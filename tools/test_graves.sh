#!/bin/bash
# Graves in multiplayer (#804). A --host with --death-drops grave places a thornback by the spawn (threats on); a --join
# "victim" (given berries) stands there, is bitten to death through the server, and everything it carried goes into a
# grave the server owns for it. The client must hear the grave (owned by "victim") and be left with nothing; after the
# respawn it walks back to the grave and the server hands it all back (the grave is gone, the items are carried again).
# EXE: a GNS build (default build/gns/bin/game.exe). PORT (27963). Prints "=== test_graves: PASS|FAIL". Run it with Git Bash.
set -u
here="$(cd "$(dirname "$0")/.." && pwd)"
exe="${EXE:-$here/build/gns/bin/game.exe}"
port="${PORT:-27963}"
tmp="$(mktemp -d)"
cd "$here"
SECONDS=0
printf 'erosion 0\n' > "$tmp/host.sav.worldgen"
RAYLIBGAME_HOST_GIVE="berries 3" RAYLIBGAME_HOST_HERD="thornback 1" RAYLIBGAME_HOST_HERD_THREATS=1 "$exe" --host "$tmp/host.sav" --seed 42 --day --bench 600 --port "$port" --death-drops grave > "$tmp/host.log" 2>&1 &
hp=$!
# #795: on Windows a signal to the msys PID can miss the native game.exe; stop it by its Windows PID as well.
stop_host() {
    local w; w=$(cat "/proc/$hp/winpid" 2>/dev/null)
    kill "$hp" 2>/dev/null
    for _ in 1 2 3 4 5; do kill -0 "$hp" 2>/dev/null || break; sleep 1; done
    wait "$hp" 2>/dev/null
    if [ -n "$w" ] && tasklist //FI "PID eq $w" //FI "IMAGENAME eq game.exe" 2>/dev/null | grep -q game.exe; then taskkill //F //PID "$w" >/dev/null 2>&1; fi
    true
}
trap stop_host EXIT
[ -d /proc ] || echo "test_graves: warning: not Git Bash (no /proc): the host may be left running; run it with Git Bash" >&2
until grep -q "test herd" "$tmp/host.log" 2>/dev/null || ! kill -0 "$hp" 2>/dev/null || [ $SECONDS -gt 240 ]; do sleep 1; done
grep -h "test herd" "$tmp/host.log" || { echo "test_graves: no thornback on the host -> FAIL"; stop_host; exit 1; }
RAYLIBGAME_VICTIM_SECONDS=45 RAYLIBGAME_NET_ANIMALS_ROLE=victim timeout 200 "$exe" --join "127.0.0.1:$port" --name victim --bench 3 --test-net-animals > "$tmp/victim.log" 2>&1
stop_host
ok=1
chk() { if [ "$1" = 1 ]; then echo "  ok: $2"; else echo "  MISSING: $2"; ok=0; fi; }
sv=$(grep "^bench: net survival (victim" "$tmp/victim.log"); echo "$sv"
gv=$(grep "^bench: net graves (victim" "$tmp/victim.log"); echo "$gv"
atDeath=$(echo "$gv" | grep -o "carried [0-9-]* at death" | awk '{print $2}')
after=$(echo "$gv" | grep -o ", [0-9-]* after it" | awk '{print $2}')
back=$(echo "$gv" | grep -o "walked back: carried [0-9-]*" | awk '{print $4}')
chk "$(echo "$sv" | grep -q "died 1" && echo 1)" "the bites killed it on the server"
chk "$([ -n "$atDeath" ] && [ "$atDeath" -gt 0 ] && echo 1)" "it carried its starting kit when it died ($atDeath of telescope, map)"
chk "$([ "$after" = 0 ] && echo 1)" "death_drops grave: the server took everything ($after of them carried while dead)"
chk "$(echo "$gv" | grep -q "own grave seen 1" && echo 1)" "the client heard of its grave (owned by victim)"
chk "$(echo "$gv" | grep -q "grave gone 1" && echo 1)" "walking back to it emptied it (the grave is gone)"
chk "$([ -n "$back" ] && [ -n "$atDeath" ] && [ "$back" -ge "$atDeath" ] && echo 1)" "and everything came back ($back of $atDeath)"
echo "=== test_graves: $([ $ok = 1 ] && echo PASS || echo FAIL) after ${SECONDS} s ($tmp)"
[ $ok = 1 ]
