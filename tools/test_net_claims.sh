#!/bin/bash
# Land claims (#343). Two players of two teams (of one each) join a peaceful --host that gives each a chest and a claim
# banner. "owner" raises its banner and builds a chest in the claim; "other" (3 m off) hears it is inside, sees the claim
# listed (NetMsg 41), then digs by the chest, builds there, opens and breaks the chest and the banner: all refused. The
# owner allies the other team (25 s): the other breaks the chest (refunded). The owner breaks its banner (32 s): no claim.
# EXE: a GNS build (default build/gns/bin/game.exe). PORT (27969). Prints "=== test_net_claims: PASS|FAIL". Git Bash.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27969}"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
printf 'erosion 0\n' > "$tmp/host.sav.worldgen"
. "$here/tools/lib_netmem.sh"
net_wait_for_memory 12 test_net_claims || exit 0
RAYLIBGAME_HOST_DIFFICULTY=0 RAYLIBGAME_HOST_GIVE="chest 1,claim_banner 1" "$exe" --host "$tmp/host.sav" --seed 42 --bench 600 --port "$port" > "$tmp/host.log" 2>&1 &
hp=$!
stop_host() { # #795: by the Windows PID too, only while it is still our game.exe
    local w; w=$(cat "/proc/$hp/winpid" 2>/dev/null)
    kill "$hp" 2>/dev/null
    for _ in 1 2 3 4 5; do kill -0 "$hp" 2>/dev/null || break; sleep 1; done
    wait "$hp" 2>/dev/null
    if [ -n "$w" ] && tasklist //FI "PID eq $w" //FI "IMAGENAME eq game.exe" 2>/dev/null | grep -q game.exe; then taskkill //F //PID "$w" >/dev/null 2>&1; fi
    true
}
trap stop_host EXIT
until grep -q "listening on port" "$tmp/host.log" 2>/dev/null || ! kill -0 "$hp" 2>/dev/null || [ $SECONDS -gt 240 ]; do sleep 1; done
until grep -q "world ready" "$tmp/host.log" 2>/dev/null || [ $SECONDS -gt 280 ]; do sleep 1; done
join() { RAYLIBGAME_CLAIMS_SECONDS=40 RAYLIBGAME_NET_CLAIMS_ROLE=$1 timeout 300 "$exe" --join "127.0.0.1:$port" --name "$1" --bench 3 $NET_LIGHT_CLIENT --test-net-claims > "$tmp/$1.log" 2>&1; }
join other & wp=$!
sleep 2
join owner
wait $wp
stop_host
ok=1
chk() { if [ "$1" = 1 ]; then echo "  ok: $2"; else echo "  MISSING: $2"; ok=0; fi; }
o=$(grep -h "^bench: net claims (owner" "$tmp/owner.log"); t=$(grep -h "^bench: net claims (other" "$tmp/other.log"); echo "$o"; echo "$t"
chk "$(echo "$o" | grep -q "banner 1, chest placed 1; claims listed 1" && echo 1)" "the owner raised a claim banner (a claim, in its list) and built a chest in it"
chk "$(echo "$t" | grep -q "claims listed 1 (not allowed 1, then allowed 1), now 0" && echo 1)" "the other's list (NetMsg 41): the claim, not allowed, then allowed (ally), then gone"
chk "$(echo "$t" | grep -q "entered 1, tried 1" && echo 1)" "the other team's player heard it was inside the claim"
chk "$(echo "$t" | grep -q "dig refused (claim) 1" && echo 1)" "its dig there was refused (Claim)"
chk "$(echo "$t" | grep -q "build/break refused (claim) 1 (3 parts refused)" && echo 1)" "building there, breaking the chest and breaking the banner were refused (Claim)"
chk "$(echo "$t" | grep -q "open refused 1" && echo 1)" "opening the chest was refused"
chk "$(echo "$t" | grep -q "allied 1, removed as an ally 1: chests 1 -> 2, chest standing 0" && echo 1)" "allied (/claim ally), it could break the chest (refunded to it)"
chk "$(echo "$o" | grep -q "banner broken 1, standing 0" && echo "$o" | grep -q "), now 0;" && echo 1)" "the owner broke its banner: the claim went"
echo "=== test_net_claims: $([ $ok = 1 ] && echo PASS || echo FAIL) after ${SECONDS} s ($tmp)"
[ $ok = 1 ]
