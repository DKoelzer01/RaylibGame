#!/bin/bash
# busybox (w64devkit on PATH) is not bash: re-run under Git's bash (#1477)
if [ -z "${BASH_VERSINFO:-}" ]; then for b in "C:/Program Files/Git/bin/bash.exe" "/c/Program Files/Git/bin/bash.exe" "/usr/bin/bash"; do [ -x "$b" ] && exec "$b" "$0" "$@"; done; echo "needs bash" >&2; exit 2; fi
# The PC's machine slot (#1477): one session runs the game, tests, benches or captures at a time.
# State lives in the main checkout, shared by every worktree: .machine-lock (who holds the machine) and
# .machine-queue (who waits). A session is named by its worktree folder (override with SLOT_NAME).
#
#   tools/slot.sh request <minutes> [quiet] <what>   queue, then wait (up to 60 min) until you hold the machine
#                                                    (run it in the background; it returns when it's yours)
#   tools/slot.sh take <minutes> [quiet] <what>      take it now if it's free (or already yours)
#   tools/slot.sh extend <minutes>                   more time on your hold
#   tools/slot.sh free                               release it (and drop your queue entries)
#   tools/slot.sh status                             holder, quiet or not, until when, and the queue
#   tools/slot.sh wait-quiet                         block until no one else's quiet window is on (max 30 min)
#   tools/slot.sh holds [name]                       exit 0 if name (default: this session) holds the machine
#   tools/slot.sh quiet-other                        exit 0 (and print the holder) if someone ELSE holds a quiet slot
#   tools/slot.sh queue-rm <name> | queue-top <name> coordinator: drop or move an entry to the front
# "quiet" = nobody else compiles or runs anything (timing, A/B, screenshots that need a calm machine).
# A hold expires at its end time (+2 min grace), so a crashed session can't block the machine.

common=$(git rev-parse --path-format=absolute --git-common-dir 2>/dev/null) || { echo "slot: not in a git checkout" >&2; exit 2; }
ROOT=$(dirname "$common")
LOCK="$ROOT/.machine-lock"; QUEUE="$ROOT/.machine-queue"; MUTEX="$ROOT/.machine-lock.d"
me=${SLOT_NAME:-$(basename "$(git rev-parse --show-toplevel 2>/dev/null)")}
now() { date +%s; }
hhmm() { date -d "@$1" +%H:%M 2>/dev/null || echo "$1"; }

lockon() { local t=0; until mkdir "$MUTEX" 2>/dev/null; do sleep 0.2; t=$((t+1)); [ $t -gt 150 ] && rm -rf "$MUTEX"; done; }
lockoff() { rmdir "$MUTEX" 2>/dev/null; }

# current holder fields: holder|mode|until|what (empty when free or expired)
current() {
  [ -f "$LOCK" ] || return 1
  IFS='|' read -r h m u w < "$LOCK"
  [ -n "$h" ] && [ "$(now)" -le $(( u + 120 )) ] || return 1
  echo "$h|$m|$u|$w"
}
holder() { current | cut -d'|' -f1; }

parse() { # minutes [quiet] what...
  mins=$1; shift; mode=run
  [ "$1" = quiet ] && { mode=quiet; shift; }
  what="$*"
  [[ "$mins" =~ ^[0-9]+$ ]] || { echo "slot: minutes must be a number" >&2; exit 2; }
  # 15 min at most (the test-batch limit), except the user's playtest hold (SLOT_NAME=playtest, set by the coordinator)
  # and a full regression (SLOT_NAME=regression).
  case "$me" in playtest|regression) ;; *) [ "$mins" -le 15 ] || { echo "slot: 15 minutes at most (the test-batch limit); split longer work" >&2; exit 2; } ;; esac
}

write_hold() { echo "$me|$mode|$(( $(now) + mins*60 ))|$what" > "$LOCK"; }

cmd=$1; shift
case "$cmd" in
  status)
    if c=$(current); then IFS='|' read -r h m u w <<< "$c"; echo "held by $h ($m) until $(hhmm $u): $w"; else echo "free"; fi
    if [ -s "$QUEUE" ]; then echo "queue:"; n=1; while IFS='|' read -r qn qm qmode qw qt; do echo "  $n. $qn $qm min $qmode: $qw"; n=$((n+1)); done < "$QUEUE"; else echo "queue: empty"; fi ;;
  take)
    parse "$@"; lockon
    h=$(holder)
    if [ -z "$h" ] || [ "$h" = "$me" ]; then write_hold; lockoff; echo "slot: $me holds the machine ($mode) until $(hhmm $(( $(now)+mins*60 )))"
    else lockoff; echo "slot: busy, held by $h. Use: tools/slot.sh request $mins ${mode/run/} $what" >&2; exit 1; fi ;;
  request)
    parse "$@"
    lockon; grep -v "^$me|" "$QUEUE" 2>/dev/null > "$QUEUE.tmp"; echo "$me|$mins|$mode|$what|$(now)" >> "$QUEUE.tmp"; mv "$QUEUE.tmp" "$QUEUE"; lockoff
    end=$(( $(now) + 3600 )); last=""
    while :; do
      lockon
      h=$(holder); first=$(head -1 "$QUEUE" 2>/dev/null | cut -d'|' -f1)
      if { [ -z "$h" ] || [ "$h" = "$me" ]; } && [ "$first" = "$me" ]; then
        write_hold; tail -n +2 "$QUEUE" > "$QUEUE.tmp"; mv "$QUEUE.tmp" "$QUEUE"; lockoff
        echo "GO: $me holds the machine ($mode) until $(hhmm $(( $(now)+mins*60 ))): $what"
        echo "When done: tools/slot.sh free"; exit 0
      fi
      pos=$(grep -n "^$me|" "$QUEUE" | cut -d: -f1); lockoff
      [ -z "$pos" ] && { echo "slot: your request was removed from the queue" >&2; exit 1; }
      msg="waiting: position $pos, machine ${h:-free}"; [ "$msg" != "$last" ] && { echo "$msg"; last=$msg; }
      [ "$(now)" -gt $end ] && { lockon; grep -v "^$me|" "$QUEUE" > "$QUEUE.tmp"; mv "$QUEUE.tmp" "$QUEUE"; lockoff; echo "slot: gave up after 60 min (dropped from the queue)" >&2; exit 1; }
      sleep 10
    done ;;
  extend)
    lockon; c=$(current); IFS='|' read -r h m u w <<< "$c"
    if [ "$h" = "$me" ]; then echo "$h|$m|$(( $(now) + $1*60 ))|$w" > "$LOCK"; lockoff; echo "slot: extended to $(hhmm $(( $(now)+$1*60 )))"
    else lockoff; echo "slot: you don't hold the machine" >&2; exit 1; fi ;;
  free)
    lockon; h=$(holder)
    if [ "$h" = "$me" ] || [ "$1" = "--force" ]; then rm -f "$LOCK"; fi
    grep -v "^$me|" "$QUEUE" 2>/dev/null > "$QUEUE.tmp"; mv "$QUEUE.tmp" "$QUEUE" 2>/dev/null
    lockoff; echo "slot: free" ;;
  holds) [ "$(holder)" = "${1:-$me}" ] ;;
  quiet-other)
    c=$(current) || exit 1; IFS='|' read -r h m u w <<< "$c"
    [ "$m" = quiet ] && [ "$h" != "$me" ] && { echo "$h until $(hhmm $u)"; exit 0; }; exit 1 ;;
  wait-quiet)
    end=$(( $(now) + 1800 ))
    while q=$("$0" quiet-other); do [ "$(now)" -gt $end ] && { echo "slot: quiet window still on after 30 min ($q)" >&2; exit 1; }; sleep 15; done ;;
  queue-rm) lockon; grep -v "^$1|" "$QUEUE" > "$QUEUE.tmp"; mv "$QUEUE.tmp" "$QUEUE"; lockoff ;;
  queue-top) lockon; { grep "^$1|" "$QUEUE"; grep -v "^$1|" "$QUEUE"; } > "$QUEUE.tmp"; mv "$QUEUE.tmp" "$QUEUE"; lockoff ;;
  *) sed -n '2,20p' "$0"; exit 2 ;;
esac
