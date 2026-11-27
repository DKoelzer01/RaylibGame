#!/bin/bash
# Dedicated server (#326): --server with a config file, console commands on stdin, a log file, "stop" with a final save.
# A sandbox world (up in seconds). EXE (default build/check/bin/game.exe; no network needed: listening may be unavailable).
# Prints "=== test_server: PASS|FAIL" and exits 0 on PASS. Capped at 3 minutes.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/check/bin/game.exe}"
tmp="$(mktemp -d)"
win="$(cygpath -m "$tmp" 2>/dev/null || echo "$tmp")"
cat > "$tmp/server.cfg" <<EOF
# test_server.sh
save = $win/world/world.sav
port = ${PORT:-27971}
max_players = 4
password = secret
autosave_minutes = 0
motd = hello from test_server
seed = 42
log = $win/logs/server.log
EOF
cd "$here"
( until grep -q "server: ready" "$tmp/out.txt" 2>/dev/null || [ $SECONDS -gt 150 ]; do sleep 1; done
  sleep 2; echo help; echo status; echo list; echo "say hello there"; echo nonsense
  echo "op Alice"; echo "op Bob 0123456789abcdef"; echo "ban mallory griefing"; echo "whitelist on"; echo ops; echo save; sleep 2; echo stop ) \
    | timeout 180 "$exe" --server "$win/server.cfg" --sandbox > "$tmp/out.txt" 2>&1
rc=$?
ok=1
need() { if grep -q "$1" "$2"; then echo "  ok: $3"; else echo "  MISSING: $3"; ok=0; fi; }
need "server: config .*max players 4, password set, autosave 0 min, motd: hello from test_server" "$tmp/out.txt" "config read"
need "server: ready" "$tmp/out.txt" "ready line (Pterodactyl's startup detection)"
need "server: commands:" "$tmp/out.txt" "help"
need "server: 0 player(s) online of 4" "$tmp/out.txt" "status"
need "server: \[server\] hello there" "$tmp/out.txt" "say"
need "unknown command 'nonsense'" "$tmp/out.txt" "unknown command"
need "server: saved to .* done" "$tmp/out.txt" "save command"
need "host: saved on shutdown to .* done" "$tmp/out.txt" "final save on stop"
need "server: stopped" "$tmp/out.txt" "clean stop"
need "^\[20[0-9-]* [0-9:]*\] server: ready" "$tmp/logs/server.log" "log file, time-stamped"
need "alice hasn't played on this server" "$tmp/out.txt" "op refuses a name with no player id (#1868)"
need "server: bob is now an op" "$tmp/out.txt" "op NAME UUID from the console (#345, #1868)"
need "^op bob id:0123456789abcdef" "$tmp/world/admin.txt" "op stored with its player id (#1868)"
need "server: mallory is banned (griefing)" "$tmp/out.txt" "ban (#345)"
need "^ban mallory griefing" "$tmp/world/admin.txt" "admin.txt next to the save (#345)"
need "console: whitelist on" "$tmp/world/audit.log" "audit log (#345)"
[ -f "$tmp/world/world.sav" ] && echo "  ok: save on disk" || { echo "  MISSING: save on disk"; ok=0; }
[ $rc -eq 0 ] && echo "  ok: exit code 0" || { echo "  MISSING: exit code 0 (got $rc)"; ok=0; }
echo "=== test_server: $([ $ok -eq 1 ] && echo PASS || echo FAIL) ($tmp)"
[ $ok -eq 1 ]
