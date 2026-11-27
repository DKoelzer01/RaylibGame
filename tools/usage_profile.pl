#!/usr/bin/perl
# Token usage profile of the RaylibGame Claude Code sessions (#1880, tracking #1874): where the context cost goes.
#
#   perl tools/usage_profile.pl [SINCE] [UNTIL] [transcript dirs or .jsonl files...]
#
#   SINCE / UNTIL   YYYY-MM-DD (default: the last 7 days, until now)
#   dirs            default: every ~/.claude/projects/*RaylibGame* folder (main sessions and their subagents)
#
# Prints API-equivalent cost by main session vs subagent and model (input, output, cache read, cache write), the number
# of API calls and their mean context, the context-size spread, and tool-result tokens by tool and by shell command type,
# weighted by how many later calls carried them ("carried", until a compaction or the end of the session). Tokens of
# tool results are estimated as characters / 4. Prices: API list prices per 1M tokens (see %price); a subscription is
# billed differently, but its usage limits follow the same shape.
use strict; use warnings; use JSON::PP; use File::Find; use POSIX qw(strftime);

my @dates = grep { /^\d{4}-\d\d-\d\d$/ } @ARGV;
my @paths = grep { !/^\d{4}-\d\d-\d\d$/ } @ARGV;
my $since = $dates[0] // strftime('%Y-%m-%d', gmtime(time - 7 * 86400));
my $until = $dates[1] // '9999';
if (!@paths) {
    my $home = $ENV{USERPROFILE} // $ENV{HOME}; $home =~ s{\\}{/}g;
    @paths = glob("$home/.claude/projects/*RaylibGame*");
}
my (@files, %seenf);
my ($sy, $sm, $sd) = split /-/, $since;
my $since_t = POSIX::mktime(0, 0, 0, $sd, $sm - 1, $sy - 1900);
for my $p (@paths) {
    if (-f $p) { push @files, $p unless $seenf{$p}++; next }
    find(sub { push @files, $File::Find::name if /\.jsonl$/ && (stat $_)[9] >= $since_t && !$seenf{$File::Find::name}++ }, $p) if -d $p;
}
die "usage_profile: no transcripts found\n" unless @files;

# $ per 1M tokens: input, output, cache read, cache write (5-minute writes, 1.25x input)
my %price = (opus => [4, 20, 0.20, 5], sonnet => [2, 10, 0.10, 2.5], haiku => [0.10, 0.50, 0.01, 0.125]);
my $json = JSON::PP->new;
my (%tok, %cost, %ncalls, %ctxsum, @hist, %added, %carried, %shadd, %shcar);
my $carried_total = 0; my $images = 0;

sub shell_kind {
    my $cmd = shift // return '';
    return 'quiet.sh (compact)' if $cmd =~ /quiet\.sh/;
    return 'build'              if $cmd =~ /cmake --build|\bmake\b|\bninja\b/;
    return 'game/test run'      if $cmd =~ /game(\.exe)?\b.*--|game_tests|regression\.sh|--scenario|sandbox\.sh|\btest_\w+\.sh/;
    return 'gh (CI, PRs, issues)' if $cmd =~ /\bgh\b/;
    return 'git'                if $cmd =~ /\bgit\b/;
    return 'slot'               if $cmd =~ /slot\.sh/;
    return 'shell read/search'  if $cmd =~ /\b(grep|rg|find|ls|cat|head|tail|sed -n|wc|awk)\b/;
    return 'other shell';
}

for my $file (@files) {
    my $kind = ($file =~ m{[/\\]subagents[/\\]}) ? 'subagent' : 'main';
    open my $fh, '<', $file or next;
    my (%seen, %toolname, %shell, @open);
    my $calls = 0;
    my $close = sub { for my $o (@open) { my $w = $o->[1] * ($calls - $o->[2]); $carried{$o->[0]} += $w; $shcar{$o->[3]} += $w if $o->[3]; $carried_total += $w } @open = () };
    while (my $line = <$fh>) {
        my ($day) = $line =~ /"timestamp":"(\d{4}-\d\d-\d\d)/ or next;
        my $in = $day ge $since && $day le $until;
        if ($line =~ /"subtype":"compact_boundary"/) { $close->(); next }
        if ($line =~ /^\{"parentUuid".*?"type":"assistant"|"type":"assistant".*"usage"/ && $line =~ /"usage"/) {
            my $d = eval { $json->decode($line) } or next;
            my $m = $d->{message}; ref $m eq 'HASH' or next;
            for my $c (@{ ref $m->{content} eq 'ARRAY' ? $m->{content} : [] }) {
                next unless ($c->{type} // '') eq 'tool_use';
                my $n = $c->{name} // '?';
                $n = "mcp:$1" if $n =~ /^mcp__(.+?)__/;
                $toolname{$c->{id}} = $n;
                $shell{$c->{id}} = shell_kind($c->{input}{command}) if $n eq 'Bash' || $n eq 'PowerShell';
            }
            my $u = $m->{usage} or next;
            next if $seen{$m->{id} // ''}++;
            $calls++;
            next unless $in;
            my ($model) = ($m->{model} // '') =~ /(opus|sonnet|haiku)/; $model //= 'other';
            my $k = "$kind/$model";
            my @t = map { $_ // 0 } @$u{qw(input_tokens output_tokens cache_read_input_tokens cache_creation_input_tokens)};
            $tok{$k}[$_] += $t[$_] for 0..3;
            my $p = $price{$model} // [0, 0, 0, 0];
            $cost{$k} += $t[$_] * $p->[$_] / 1e6 for 0..3;
            my $ctx = $t[0] + $t[2] + $t[3];
            $ncalls{$k}++; $ctxsum{$k} += $ctx;
            $hist[$ctx >= 1e6 ? 10 : int($ctx / 1e5)]++ if $kind eq 'main';
        } elsif ($in && $line =~ /"tool_use_id"/) {
            # tool results, sized without a JSON decode (these lines are big): only the "message" part counts (the line
            # repeats the result in "toolUseResult", sometimes with the whole file), split at each tool_result block
            my $msg = $line;
            my $cut = index($msg, '"toolUseResult":'); $msg = substr($msg, 0, $cut) if $cut > 0;
            my @starts;
            # a block starts with either key: {"type":"tool_result","tool_use_id":...} or {"tool_use_id":...,"type":...}
            while ($msg =~ /[\[,]\{"(?:tool_use_id"|type":"tool_result")/g) { push @starts, $-[0] }
            for my $i (0..$#starts) {
                my $seg = substr($msg, $starts[$i], ($i < $#starts ? $starts[$i + 1] : length $msg) - $starts[$i]);
                my ($id) = $seg =~ /"tool_use_id":"([^"]+)"/ or next;
                # images (screenshots, captures) are base64 in the transcript but cost ~1.6K tokens each to the model
                my $img = 0; my $b64 = 0;
                while ($seg =~ /"data":"([A-Za-z0-9+\/=]{1000,})"/g) { $img++; $b64 += length $1 }
                my $len = length($seg) - $b64 - 120; $len = 0 if $len < 0;
                my $tk = $len / 4 + 1600 * $img;
                $images += $img;
                my $n = $toolname{$id} // 'unknown';
                my $sk = $shell{$id} // '';
                $added{$n} += $tk; $shadd{$sk} += $tk if $sk;
                push @open, [$n, $tk, $calls, $sk];
            }
        }
    }
    $close->();
}

my $tc = 0; $tc += $_ for values %cost;
printf "Claude Code usage %s .. %s: %d transcripts, API-equivalent \$%.0f\n\n", $since, ($until eq '9999' ? 'now' : $until), scalar @files, $tc;
printf "%-16s %8s %9s %11s %14s %12s %8s %6s %8s\n", qw(who calls mean_ctx output cache_read cache_write input $ share);
for my $k (sort { $cost{$b} <=> $cost{$a} } keys %cost) {
    printf "%-16s %8d %8.0fK %11d %14d %12d %8d %6.0f %5.1f%%\n", $k, $ncalls{$k}, $ctxsum{$k} / $ncalls{$k} / 1e3,
        $tok{$k}[1], $tok{$k}[2], $tok{$k}[3], $tok{$k}[0], $cost{$k}, 100 * $cost{$k} / ($tc || 1);
}
my $mo = $cost{'main/opus'} // 0;
if ($mo) {
    my $p = $price{opus};
    printf "\nmain/opus cost split: cache read %.0f%%, cache write %.0f%%, output %.0f%%, input %.0f%%\n",
        map { 100 * $tok{'main/opus'}[$_] * $p->[$_] / 1e6 / $mo } (2, 3, 1, 0);
}
print "\nmain-session calls by context size:\n";
for my $b (0..10) { next unless $hist[$b]; printf "  %s: %d\n", ($b == 10 ? '>= 1M' : sprintf('%3d-%dK', $b * 100, $b * 100 + 100)), $hist[$b] }
my $cr = 0; $cr += $tok{$_}[2] for keys %tok;
printf "\ntool results: tokens added, and token-turns carried by later calls (%.0fM in all, %.0f%% of all cache reads; %d images at ~1.6K each)\n",
    $carried_total / 1e6, 100 * $carried_total / ($cr || 1), $images;
for my $k (sort { $carried{$b} <=> $carried{$a} } keys %carried) {
    next if $carried{$k} < $carried_total / 300;
    printf "  %-26s added %8.0fK  carried %8.1fM (%4.1f%%)\n", $k, $added{$k} / 1e3, $carried{$k} / 1e6, 100 * $carried{$k} / ($carried_total || 1);
}
print "\nshell results by command type:\n";
for my $k (sort { $shcar{$b} <=> $shcar{$a} } keys %shcar) {
    printf "  %-26s added %8.0fK  carried %8.1fM (%4.1f%%)\n", $k, $shadd{$k} / 1e3, $shcar{$k} / 1e6, 100 * $shcar{$k} / ($carried_total || 1);
}
