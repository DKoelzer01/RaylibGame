#!/bin/bash
# Dedicated server admin end to end (#326 part 2, #345 part 2): a --server (sandbox world) and --join clients as separate
# processes. The console ops alice; alice joins (the motd reaches her), chats, and as an op bans mallory from chat - mallory,
# online, is kicked with the reason; the console's "say" reaches alice; "list" shows who is on; mallory can't come back
# (banned); the whitelist refuses carol; the console kicks alice; the audit log has it all.
# EXE: a GNS build (default build/gns/bin/game.exe). Prints "=== test_server_admin: PASS|FAIL". Capped at 6 minutes.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27991}"
tmp="$(mktemp -d)"; win="$(cygpath -m "$tmp" 2>/dev/null || echo "$tmp")"
cat > "$tmp/server.cfg" <<EOF
save = $win/world/world.sav
port = $port
max_players = 8
autosave_minutes = 0
motd = welcome to the admin test
log = $win/world/server.log
EOF
cd "$here"
SECONDS=0
out="$tmp/server.txt"
waitfor() { until grep -q "$1" "$2" 2>/dev/null || [ $SECONDS -gt 300 ]; do sleep 1; done; }
join() { # name [chat] [log] -> background client, log in $tmp/<log or name>.txt
    RAYLIBGAME_JOIN_CHAT="${2:-}" timeout 200 "$exe" --join "127.0.0.1:$port" --name "$1" --sandbox > "$tmp/${3:-$1}.txt" 2>&1 &
}
{
    waitfor "server: ready" "$out"
    echo "op alice"
    join mallory
    waitfor "mallory joined" "$out"
    join alice "hello everyone|/ban mallory griefing|/list"
    waitfor "alice joined" "$out"
    waitfor "mallory is banned" "$out"; sleep 2
    echo "say hi from the console"; echo list; sleep 3
    join mallory "" mallory2; waitfor "refused mallory" "$out"
    echo "whitelist on"; sleep 1
    join carol; waitfor "refused carol" "$out"
    echo "kick alice enough for today"; sleep 3
    echo stop
} | timeout 360 "$exe" --server "$win/server.cfg" --sandbox > "$out" 2>&1
sleep 3
ok=1
need() { if grep -q "$1" "$2"; then echo "  ok: $3"; else echo "  MISSING: $3"; ok=0; fi; }
need "join: motd: welcome to the admin test" "$tmp/alice.txt" "the motd reaches a joining player"
need "chat: alice: hello everyone" "$out" "a player's chat on the server"
need "server: mallory is banned (griefing) (by alice)" "$out" "an op's /ban from chat"
need "server: kicked mallory (banned: griefing)" "$out" "the ban kicks mallory"
need "chat: You were kicked: banned: griefing" "$tmp/mallory.txt" "mallory hears why"
need "join: disconnected (kicked)" "$tmp/mallory.txt" "mallory's client was kicked"
need "chat: \[server\] hi from the console" "$tmp/alice.txt" "the console's say reaches alice"
need "server:   alice (session [0-9]*, op," "$out" "list: alice, an op"
need "server: refused mallory (banned: griefing)" "$out" "a banned player can't rejoin"
need "server: refused carol (not on this server's whitelist)" "$out" "the whitelist refuses carol"
need "chat: online: .*alice" "$tmp/alice.txt" "an op's /list from chat"
need "server: alice was kicked" "$out" "the console kicks alice"
need "chat: You were kicked: enough for today" "$tmp/alice.txt" "alice hears why"
need "alice: ban mallory griefing" "$tmp/world/audit.log" "audit: the op's ban"
need "refused carol" "$tmp/world/audit.log" "audit: the refusal"
need "join alice from" "$tmp/world/audit.log" "audit: joins"
echo "=== test_server_admin: $([ $ok -eq 1 ] && echo PASS || echo FAIL) after ${SECONDS} s ($tmp)"
[ $ok -eq 1 ]
