#!/bin/bash
# busybox (w64devkit on PATH) is not bash: re-run under Git's bash (#1477)
if [ -z "${BASH_VERSINFO:-}" ]; then for b in "C:/Program Files/Git/bin/bash.exe" "/c/Program Files/Git/bin/bash.exe" "/usr/bin/bash"; do [ -x "$b" ] && exec "$b" "$0" "$@"; done; echo "needs bash" >&2; exit 2; fi
# (#1897) Restart a topic session with a fresh context at a task boundary (user, 2026-11-26: the coordinator may stop and
# restart topic sessions itself). Long sessions cost more per call: every call re-reads the whole context (#1874).
#
#   tools/coord/cycle_session.sh status                   every background session: context tokens, state, cwd
#   tools/coord/cycle_session.sh check <id|name>          exit 0 when it is over the threshold (250K), 1 when not
#   tools/coord/cycle_session.sh restart <id|name> [--dry] [--cwd <dir>] [--handoff <file>]
#                                                         stop it, take it off the list, start it again (same name, cwd,
#                                                         permission mode and model) with --autocompact 300000
#     --cwd <dir>       start the new session there instead (and look for its handoff there) (#1928)
#     --handoff <file>  the handoff the new session reads; default: <cwd>/build/session-handoff-<slug>.md when it
#                       exists (slug = the name in lower case, dashes: "UI/UX and Adventure" -> ui-ux-and-adventure),
#                       else <cwd>/build/session-handoff.md. The prompt names its absolute path.
# Two live sessions in one cwd would read each other's build/session-handoff.md: status marks them (!), restart warns.
#
# The coordinator's steps, with its own tools in between:
#   1. `check <name>`; over the threshold: SendMessage the session "finish your task, write build/session-handoff.md
#      (what's done, PRs, what's next), then stop", with notify_when_idle: true, and wait for the idle notice.
#   2. `restart <name>`. Peers keep messaging it by name; the user opens it with `claude attach "<name>"`.
# The old conversation is kept (`claude --resume` lists it).
#
# `claude rm` deletes the session's worktree when Claude created it and it has no uncommitted changes (tested
# 2026-11-26, #1897): the restarted session would lose its folder and any unpushed commits. So the old entry is removed
# only when its cwd is not under .claude/worktrees/; otherwise it stays listed (stopped) and this script says so.
#
# Env: CYCLE_THRESHOLD (tokens, default 250000), CYCLE_AUTOCOMPACT (default 300000).
set -u
threshold=${CYCLE_THRESHOLD:-250000}
autocompact=${CYCLE_AUTOCOMPACT:-300000}
projects="$HOME/.claude/projects"

# id <TAB> sessionId <TAB> name <TAB> status <TAB> state <TAB> cwd, one line per background session
agents() {
    claude agents --json 2>/dev/null | perl -MJSON::PP -0777 -ne '
        my $a = eval { decode_json($_) } || [];
        for my $s (@$a) { next unless ($s->{kind} // "") eq "background";
            print join("\t", map { $s->{$_} // "" } qw(id sessionId name status state cwd)), "\n"; }'
}

transcript() { ls "$projects"/*/"$1".jsonl 2>/dev/null | head -1; }

# the context of the session's last call (input + cache read + cache write), and its permission mode and model
# (the last ones the transcript names): "<tokens> <mode> <model>"
profile() {
    local t; t=$(transcript "$1")
    [ -n "$t" ] || { echo "? - -"; return; }
    tail -c 4000000 "$t" | perl -ne '
        if (/"usage":\{([^}]*)\}/) { my $u = $1; my $n = 0;
            $n += $1 while $u =~ /"(?:input_tokens|cache_read_input_tokens|cache_creation_input_tokens)":(\d+)/g; $ctx = $n; }
        $mode = $1 while /"permissionMode":"([^"]+)"/g;
        $model = $1 if /"type":"assistant"/ && /"model":"(claude-[^"]+)"/;
        END { printf "%s %s %s\n", $ctx // "?", $mode // "-", $model // "-" }'
}

find_one() { # <id or name> -> the agents line of that session (by id first; a name must be unique)
    local hits; hits=$(agents | awk -F'\t' -v n="$1" '$1 == n')
    [ -n "$hits" ] || hits=$(agents | awk -F'\t' -v n="$1" '$3 == n')
    [ -n "$hits" ] || { echo "cycle: no background session with the id or name '$1' (tools/coord/cycle_session.sh status)" >&2; return 1; }
    if [ "$(printf '%s\n' "$hits" | wc -l)" -gt 1 ]; then
        # (a stopped copy kept by an earlier restart shares the name; `state` is the task's state, not the process's)
        echo "cycle: several sessions are named '$1': give the id (first column of 'status'):" >&2
        printf '%s\n' "$hits" | cut -f1,3,4,6 >&2
        return 1
    fi
    printf '%s\n' "$hits"
}

normcwd() { local p=${1//\\//}; p=${p%/}; printf '%s' "${p,,}"; }
slug() { local s=${1,,}; s=$(printf '%s' "$s" | sed -E 's/[^a-z0-9]+/-/g; s/^-+//; s/-+$//'); printf '%s' "$s"; }

# the other live (not stopped) sessions whose cwd is <cwd>, excluding <id>: "name (id)" per line
cwd_peers() { # <cwd> <id>
    local want; want=$(normcwd "$1")
    agents | while IFS=$'\t' read -r id sid name status state cwd; do
        [ "$id" = "$2" ] && continue; [ "$status" = stopped ] && continue
        [ "$(normcwd "$cwd")" = "$want" ] && echo "$name ($id)"
    done
}

case "${1:-status}" in
status)
    all=$(agents)
    printf '%-26s %-9s %9s %-8s %-8s %s\n' name id context status state cwd
    while IFS=$'\t' read -r id sid name status state cwd; do
        [ -n "$id" ] || continue
        read -r ctx _ _ <<< "$(profile "$sid")"
        mark=""; [ "$ctx" != "?" ] && [ "$ctx" -ge "$threshold" ] && mark=" *"
        if [ "$status" != stopped ]; then
            n=$(printf '%s\n' "$all" | awk -F'\t' -v c="$(normcwd "$cwd")" -v me="$id" '
                $1 != me && $4 != "stopped" { p = $6; gsub(/\\/, "/", p); sub(/\/$/, "", p); if (tolower(p) == c) k++ } END { print k + 0 }')
            [ "$n" -gt 0 ] && mark="$mark !"
        fi
        printf '%-26s %-9s %9s %-8s %-8s %s%s\n' "$name" "$id" "$ctx" "$status" "$state" "$cwd" "$mark"
    done <<< "$all"
    echo "(* over $threshold tokens: cycle it at its next task boundary)"
    echo "(! another live session shares this cwd: give each its own handoff (build/session-handoff-<slug>.md) or restart with --cwd)" ;;
check)
    line=$(find_one "${2:?usage: check <name>}") || exit 2
    IFS=$'\t' read -r id sid name status state cwd <<< "$line"
    read -r ctx _ _ <<< "$(profile "$sid")"
    echo "$name ($id): $ctx tokens of context (threshold $threshold), $status/$state"
    [ "$ctx" != "?" ] && [ "$ctx" -ge "$threshold" ] ;;
restart)
    usage="usage: restart <id|name> [--dry] [--cwd <dir>] [--handoff <file>]"
    who=${2:?$usage}; shift 2
    dry=0; newcwd=""; handoff=""
    while [ $# -gt 0 ]; do
        case "$1" in
            --dry) dry=1 ;;
            --cwd) newcwd=${2:?$usage}; shift ;;
            --handoff) handoff=${2:?$usage}; shift ;;
            *) echo "cycle: unknown option '$1' ($usage)" >&2; exit 2 ;;
        esac
        shift
    done
    line=$(find_one "$who") || exit 2
    IFS=$'\t' read -r id sid name status state cwd <<< "$line"
    read -r ctx mode model <<< "$(profile "$sid")"
    [ "$status" = busy ] && { echo "cycle: '$name' is busy: wait for its idle notice first"; exit 3; }
    oldcwd=$cwd
    if [ -n "$newcwd" ]; then
        [ -d "$newcwd" ] || { echo "cycle: --cwd '$newcwd' is not a folder"; exit 2; }
        cwd=$(cd "$newcwd" && pwd -W 2>/dev/null || pwd)
    fi
    fwd=${cwd//\\//}
    if [ -z "$handoff" ]; then
        handoff="$fwd/build/session-handoff-$(slug "$name").md"
        [ -f "$handoff" ] || handoff="$fwd/build/session-handoff.md"
    fi
    case "$handoff" in /*|[A-Za-z]:*) ;; *) handoff="$fwd/$handoff" ;; esac
    peers=$(cwd_peers "$cwd" "$id")
    if [ -n "$peers" ]; then
        echo "cycle: warning: other live sessions in $cwd: $(printf '%s' "$peers" | paste -sd, -)"
        case "$handoff" in */build/session-handoff.md)
            echo "cycle: warning: they share build/session-handoff.md: check it is $name's (or pass --handoff / --cwd)" ;; esac
    fi
    prompt="You are the $name session, restarted with a fresh context (tools/coord/cycle_session.sh, #1897). Read CLAUDE.md and your handoff $handoff, check your board_inbox, then continue the work the handoff names."
    args=(--bg -n "$name" --autocompact "$autocompact")
    [ "$mode" != - ] && args+=(--permission-mode "$mode")
    [ "$model" != - ] && args+=(--model "$model")
    # rm deletes a Claude-made worktree: judge by the OLD session's cwd, which is what rm would act on
    canrm=1; case "${oldcwd//\\//}" in */.claude/worktrees/*) canrm=0 ;; esac
    echo "cycle: $name ($id, $ctx tokens, $status/$state) in $oldcwd${newcwd:+ -> restarts in $cwd}"
    [ -f "$handoff" ] || echo "cycle: warning: no handoff at $handoff: the new session starts from CLAUDE.md alone"
    if [ $dry = 1 ]; then
        echo "  claude stop $id"
        [ $canrm = 1 ] && echo "  claude rm $id" || echo "  (no rm: the cwd is a Claude worktree; rm could delete it)"
        echo "  (cd \"$cwd\" && claude ${args[*]} \"$prompt\")"
        exit 0
    fi
    claude stop "$id" || { echo "cycle: stop failed"; exit 1; }
    if [ $canrm = 1 ]; then claude rm "$id" || echo "cycle: rm $id failed: the stopped copy stays listed as '$name'"
    else echo "cycle: kept the stopped entry $id (its cwd is a Claude worktree, which rm could delete): 'claude agents' lists '$name' twice; the new one is the running one"; fi
    (cd "$cwd" && claude "${args[@]}" "$prompt") ;;
*)
    sed -n '4,28p' "$0"; exit 2 ;;
esac
