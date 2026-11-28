#!/bin/bash
# Cheap static checks for CI (#1016), seconds on any Linux runner (bash, awk, git; no build):
#
#   tools/ci/static_checks.sh [BASE]      BASE = the commit to compare against (the PR's base tip, a push's previous tip);
#                                         without it only the checks that need no history run
#
#   1. NetMsg ids (src/engine/netproto.h, enum class NetMsg) are unique, fit the uint8 wire id, and parse. Against BASE: a
#      message added, removed or renumbered and kNetProtocol changes are listed in the job summary; a NetMsg change without
#      a kNetProtocol bump is a warning annotation.
#   2. The erosion bake (src/engine/erosion*.cpp/h) changed in code (not only comments) while kErosionBakeVersion did not:
#      a warning annotation (stale <save>.erosion caches would load; CLAUDE.md "bump the version on ANY bake code change").
#   3. A new `gh project` command in tools/ (forbidden by CLAUDE.md: it burns the shared GraphQL budget): an error.
#   4. renderOrigin (the drawing-only floating origin, Plan 11 #1412) used more in a file than tools/ci/render_origin_baseline.txt
#      records, or in a new file: a warning (simulation reads body-local doubles). --render-origin-baseline rewrites it.
#   5. Every scenario lands in exactly one of the two runner scenario steps (or carries a local-only tag): else an error.
#   6. kNetProtocol (src/engine/netproto.h) is reserved in tools/ci/protocol_reservations.txt (#1480): else an error, so
#      two branches can't ship the same number for different wire formats.
#   7. A kNetProtocol / kSaveVersion / kErosionBakeVersion number this branch reserves is not reserved for another issue on origin/staging,
#      origin/ladder or in an open PR into them (tools/ci/reservation_clash.sh, #2009): an error against BASE.
#   8. Every tools/ci/known_failures.txt line's note starts "[<owner>, until <date>]", at most 7 days ahead (#2122): else an
#      error; a line past its date is a warning.
#
# Exit code 1 only for errors (duplicate / unparsable ids, checks 3 and 5). Warnings never fail the job.
# Writes a Markdown summary to $GITHUB_STEP_SUMMARY when set (else to stdout) and ::warning/::error lines for annotations.
set -uo pipefail
export LC_ALL=C # sort and join must agree on the order

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
BASE="${1:-}"
if [ -n "$BASE" ] && ! git cat-file -e "$BASE^{commit}" 2>/dev/null; then
    echo "static_checks: base $BASE is not in this clone: history checks skipped"
    BASE=""
fi
SUMMARY="${GITHUB_STEP_SUMMARY:-/dev/stdout}"
errors=0
NETPROTO=src/engine/netproto.h
EROSION_H=src/engine/erosion.h

sum() { printf '%s\n' "$*" >> "$SUMMARY"; }
annot() { # level file line message
    printf '::%s file=%s,line=%s::%s\n' "$1" "$2" "$3" "$4"
}

# ---- 1. NetMsg ids ----
# Prints "name value" per enumerator of enum class NetMsg (Count left out), or "!unparsed <entry>" for a value it can't read.
netmsg_ids() { # reads the header on stdin
    tr -d '\r' | awk '
        /enum class NetMsg/ { inside = 1; next }
        inside && /^[ \t]*\};/ { inside = 0 }
        inside { line = $0; sub(/\/\/.*/, "", line); buf = buf " " line }
        END {
            n = split(buf, parts, ","); v = -1
            for (i = 1; i <= n; i++) {
                p = parts[i]; gsub(/[ \t]/, "", p)
                if (p == "") continue
                eq = index(p, "=")
                if (eq) {
                    name = substr(p, 1, eq - 1); val = substr(p, eq + 1)
                    if (val !~ /^[0-9]+$/) { print "!unparsed " p; v = -1000000; continue }
                    v = val + 0
                } else { name = p; v = v + 1 }
                if (name == "Count") continue
                print name, v
            }
        }'
}
protocol_of() { tr -d '\r' | sed -n 's/^[[:space:]]*constexpr[[:space:]]\{1,\}uint32_t[[:space:]]\{1,\}kNetProtocol[[:space:]]*=[[:space:]]*\([0-9]\{1,\}\).*/\1/p' | head -1; }
line_of() { grep -n "$2" "$1" | head -1 | cut -d: -f1; }

sum "## Static checks"
sum ""
if [ ! -f "$NETPROTO" ]; then
    annot error "$NETPROTO" 1 "static_checks: $NETPROTO not found (moved? update tools/ci/static_checks.sh)"
    errors=$((errors + 1))
else
    ids="$(netmsg_ids < "$NETPROTO")"
    enum_line="$(line_of "$NETPROTO" 'enum class NetMsg')"
    count="$(printf '%s\n' "$ids" | grep -vc '^!' || true)"
    if [ "$count" -lt 2 ]; then
        annot error "$NETPROTO" "${enum_line:-1}" "static_checks: could not read enum class NetMsg (format changed? update tools/ci/static_checks.sh)"
        errors=$((errors + 1))
    fi
    while read -r _ entry; do
        [ -n "$entry" ] || continue
        annot error "$NETPROTO" "${enum_line:-1}" "NetMsg entry '$entry': the checker reads only plain decimal ids"
        errors=$((errors + 1))
    done < <(printf '%s\n' "$ids" | grep '^!unparsed' || true)
    dups="$(printf '%s\n' "$ids" | grep -v '^!' | awk '{ seen[$2] = seen[$2] (seen[$2] ? ", " : "") $1; n[$2]++ } END { for (v in n) if (n[v] > 1) print v ": " seen[v] }')"
    if [ -n "$dups" ]; then
        while read -r d; do
            annot error "$NETPROTO" "${enum_line:-1}" "Duplicate NetMsg id $d (two messages share a wire id)"
        done <<< "$dups"
        errors=$((errors + 1))
    fi
    big="$(printf '%s\n' "$ids" | grep -v '^!' | awk '$2 > 255 || $2 < 1 { print $1 " = " $2 }')"
    if [ -n "$big" ]; then
        annot error "$NETPROTO" "${enum_line:-1}" "NetMsg ids must be 1..255 (uint8 on the wire): $(echo $big)"
        errors=$((errors + 1))
    fi
    proto="$(protocol_of < "$NETPROTO")"
    proto_line="$(line_of "$NETPROTO" 'constexpr uint32_t kNetProtocol')"
    if [ -z "$proto" ]; then
        annot error "$NETPROTO" 1 "static_checks: kNetProtocol not found in $NETPROTO"
        errors=$((errors + 1))
    fi
    sum "**Network protocol:** kNetProtocol = ${proto:-?}, $count NetMsg ids, $( [ -z "$dups" ] && echo "no duplicates" || echo "**duplicates**")."

    if [ -n "$BASE" ] && git cat-file -e "$BASE:$NETPROTO" 2>/dev/null; then
        base_ids="$(git show "$BASE:$NETPROTO" | netmsg_ids | grep -v '^!' || true)"
        base_proto="$(git show "$BASE:$NETPROTO" | protocol_of)"
        added="$(join -v2 <(printf '%s\n' "$base_ids" | sort) <(printf '%s\n' "$ids" | grep -v '^!' | sort) || true)"
        removed="$(join -v1 <(printf '%s\n' "$base_ids" | sort) <(printf '%s\n' "$ids" | grep -v '^!' | sort) || true)"
        renumbered="$(join <(printf '%s\n' "$base_ids" | sort) <(printf '%s\n' "$ids" | grep -v '^!' | sort) | awk '$2 != $3 { print $1 " " $2 " -> " $3 }' || true)"
        changed=0
        if [ -n "$added" ]; then sum ""; sum "New NetMsg ids:"; while read -r n v; do sum "- \`$n = $v\`"; done <<< "$added"; changed=1; fi
        if [ -n "$removed" ]; then sum ""; sum "Removed NetMsg ids:"; while read -r n v; do sum "- \`$n = $v\`"; done <<< "$removed"; changed=1; fi
        if [ -n "$renumbered" ]; then sum ""; sum "Renumbered NetMsg ids:"; while read -r l; do sum "- \`$l\`"; done <<< "$renumbered"; changed=1; fi
        if [ -n "$proto" ] && [ -n "$base_proto" ] && [ "$proto" != "$base_proto" ]; then
            sum ""; sum "**kNetProtocol bumped: $base_proto -> $proto**"
            if [ "$proto" -lt "$base_proto" ]; then
                annot warning "$NETPROTO" "${proto_line:-1}" "kNetProtocol went down ($base_proto -> $proto): clients and servers of the two builds would accept each other"
            fi
        elif [ "$changed" = 1 ]; then
            annot warning "$NETPROTO" "${proto_line:-1}" "NetMsg ids changed but kNetProtocol is still $proto: bump it so mismatched builds refuse each other at the handshake"
        fi
    fi
fi

# ---- 2. Erosion bake version ----
if [ -n "$BASE" ]; then
    # Code lines changed in the bake's files: added or removed lines that are not blank or // comments.
    code_changes="$(git diff -U0 "$BASE" HEAD -- 'src/engine/erosion*.cpp' 'src/engine/erosion*.h' \
        | tr -d '\r' | grep -E '^[+-]' | grep -vE '^(\+\+\+|---) ' \
        | grep -vE '^[+-][[:space:]]*(//.*)?$' | grep -v 'kErosionBakeVersion' | wc -l)"
    if [ "$code_changes" -gt 0 ] && [ -f "$EROSION_H" ]; then
        ver_of() { tr -d '\r' | sed -n 's/^[[:space:]]*constexpr[[:space:]]\{1,\}uint32_t[[:space:]]\{1,\}kErosionBakeVersion[[:space:]]*=[[:space:]]*\([0-9]\{1,\}\).*/\1/p' | head -1; }
        ver="$(ver_of < "$EROSION_H")"
        base_ver="$(git show "$BASE:$EROSION_H" 2>/dev/null | ver_of)"
        files="$(git diff --name-only "$BASE" HEAD -- 'src/engine/erosion*.cpp' 'src/engine/erosion*.h' | tr '\n' ' ')"
        if [ -n "$ver" ] && [ "$ver" = "$base_ver" ]; then
            annot warning "$EROSION_H" "$(line_of "$EROSION_H" 'constexpr uint32_t kErosionBakeVersion')" \
                "The erosion bake's code changed ($files) but kErosionBakeVersion is still $ver: if the bake's output changes, bump it or stale .erosion caches load (ignore for changes outside the bake, e.g. ErosionGullies)"
            sum ""; sum "Erosion code changed ($files) with kErosionBakeVersion unchanged ($ver): see the warning."
        elif [ -n "$ver" ]; then
            sum ""; sum "kErosionBakeVersion: ${base_ver:-?} -> $ver."
        fi
    fi
fi

# ---- 3. No new `gh project` in tools/ ----
if [ -n "$BASE" ]; then
    pat='gh[[:space:]]+'"project"'([[:space:]]|$)'
    hits="$(git diff -U0 "$BASE" HEAD -- tools/ ':(exclude)tools/ci/static_checks.sh' | tr -d '\r' | awk -v pat="$pat" '
        /^\+\+\+ / { f = substr($0, 7); next }
        /^@@/ { split($3, a, ","); ln = substr(a[1], 2) + 0; next }
        /^\+/ { t = substr($0, 2); if (t ~ pat && t !~ /^[ \t]*#/) print f ":" ln; ln++ }')"
    if [ -n "$hits" ]; then
        while IFS=: read -r f l; do
            annot error "$f" "$l" "New 'gh project' command: forbidden in tools/ (CLAUDE.md, #582: it reads the whole board over GraphQL). Use tools/board.sh or REST."
        done <<< "$hits"
        sum ""; sum "**New \`gh project\` usage in tools/** (error): $(echo $hits)"
        errors=$((errors + 1))
    fi
fi

# ---- 4. renderOrigin doesn't creep into more code (Plan 11 F4, #1412) ----
# The render origin is for drawing; simulation reads body-local doubles. tools/ci/render_origin_baseline.txt holds each
# file's count of `renderOrigin` (word); a file whose count rises, or a new file using it, is a warning (an error once F8
# has moved simulation code off it). Counts may only fall: after removing uses, refresh the baseline with
#   bash tools/ci/static_checks.sh --render-origin-baseline
RO_BASE=tools/ci/render_origin_baseline.txt
ro_counts() { grep -rcw renderOrigin src --include=*.cpp --include=*.h | awk -F: '$2 > 0 { print $2 " " $1 }' | sort -k2; }
if [ "${1:-}" = "--render-origin-baseline" ]; then ro_counts > "$RO_BASE"; echo "static_checks: wrote $RO_BASE"; exit 0; fi
if [ -f "$RO_BASE" ]; then
    rose="$(ro_counts | tr -d '\r' | awk 'NR == FNR { base[$2] = $1; next } ($2 in base) ? $1 > base[$2] : 1 { print $2 " " (($2 in base) ? base[$2] : 0) " " $1 }' <(tr -d '\r' < "$RO_BASE") -)"
    if [ -n "$rose" ]; then
        while read -r f was now; do
            annot warning "$f" 1 "renderOrigin uses rose $was -> $now (Plan 11 #1412): simulation code reads body-local doubles; drawing may, then refresh $RO_BASE"
        done <<< "$rose"
        sum ""; sum "**renderOrigin uses rose** (warning; Plan 11 #1412): $(echo "$rose" | awk '{ printf "%s %s->%s; ", $1, $2, $3 }')"
    fi
fi

# 5. (#1367) Every scenario runs on the self-hosted runner exactly once: its tags match one of the workflow's two
#    `--scenario tag:...` lists (parts 1 and 2 of the headless scenarios), not both - unless it carries a local-only tag
#    (rendering, timing and looks belong on the PC). A new tag in neither list would silently never run there: an error.
LOCAL_ONLY_TAGS="perf visual fx light local" # local: needs the PC's real-time rates (fluid tests the runner can't keep up with, #1447)
WF=.github/workflows/selfhosted.yml
part1=$(grep -o -- '--scenario tag:[a-z0-9_,:-]*' "$WF" | sed -n 1p | sed 's/--scenario //; s/tag://g; s/,/ /g')
part2=$(grep -o -- '--scenario tag:[a-z0-9_,:-]*' "$WF" | sed -n 2p | sed 's/--scenario //; s/tag://g; s/,/ /g')
part3=$(grep -o -- '--scenario tag:[a-z0-9_,:-]*' "$WF" | sed -n 3p | sed 's/--scenario //; s/tag://g; s/,/ /g')
part4=$(grep -o -- '--scenario tag:[a-z0-9_,:-]*' "$WF" | sed -n 4p | sed 's/--scenario //; s/tag://g; s/,/ /g') # (#1882) world generation metrics: the full tier and manual runs
if [ -z "$part1" ] || [ -z "$part2" ] || [ -z "$part3" ] || [ "$part4" != worldgen ]; then
    echo "::error file=$WF::static check 5: the four headless scenario tag lists (part 1, part 2, pose twins, worldgen) were not found in the workflow, in that order"
    errors=$((errors + 1))
else
    unplaced=""
    for f in tests/scenarios/*.cpp; do
        name=$(grep -o 'SCENARIO("[^"]*")' "$f" | head -1 | sed 's/SCENARIO("\(.*\)")/\1/')
        [ -n "$name" ] || continue
        tags=$(grep -o 's.tags = {[^}]*}' "$f" | head -1 | grep -o '"[a-z0-9_-]*"' | tr -d '"' | tr '\n' ' ')
        in1=0; in2=0; in3=0; in4=0; local=0
        for t in $tags; do
            for x in $part1; do [ "$t" = "$x" ] && in1=1; done
            for x in $part2; do [ "$t" = "$x" ] && in2=1; done
            for x in $part3; do [ "$t" = "$x" ] && in3=1; done
            for x in $part4; do [ "$t" = "$x" ] && in4=1; done
            for x in $LOCAL_ONLY_TAGS; do [ "$t" = "$x" ] && local=1; done
        done
        if [ $((in1 + in2 + in3 + in4)) -gt 1 ]; then
            echo "::error file=$f::scenario $name matches more than one runner scenario step (tags: $tags): it would run twice"
            errors=$((errors + 1))
        elif [ $in1 = 0 ] && [ $in2 = 0 ] && [ $in3 = 0 ] && [ $in4 = 0 ] && [ $local = 0 ]; then
            echo "::error file=$f::scenario $name matches neither runner scenario step (tags: $tags): add a tag to a list in $WF or a local-only tag ($LOCAL_ONLY_TAGS)"
            errors=$((errors + 1)); unplaced="$unplaced $name"
        fi
    done
    sum ""; sum "**Scenario steps (#1367):** part 1 \`$part1\`, part 2 \`$part2\`, pose twins (nightly) \`$part3\`, world generation (full tier, manual runs) \`$part4\`; local-only tags \`$LOCAL_ONLY_TAGS\`${unplaced:+; not placed:$unplaced}"
fi

# ---- 6. kNetProtocol is reserved (#1480) ----
# Each wire-format change takes the next number from tools/ci/protocol_reservations.txt ("<n> | <branch> | <issue> | <what>",
# reserved in the PR before the bump). A number missing there, or marked "(unused)", is an error: two branches that bumped
# independently would both claim it, and a client of one would join a server of the other with a different format.
RES=tools/ci/protocol_reservations.txt
proto=$(grep -oE 'kNetProtocol = [0-9]+' "$NETPROTO" | head -1 | grep -oE '[0-9]+$')
if [ -z "$proto" ]; then
    echo "::error file=$NETPROTO::static check 6: kNetProtocol not found"; errors=$((errors + 1))
elif [ ! -f "$RES" ]; then
    echo "::error file=$RES::static check 6: the protocol reservation list is missing"; errors=$((errors + 1))
else
    line=$(tr -d '\r' < "$RES" | grep -E "^${proto}[[:space:]]*\|" | head -1)
    if [ -z "$line" ]; then
        echo "::error file=$NETPROTO::kNetProtocol $proto is not reserved in $RES: add a line \"$proto | <branch> | #<issue> | <what>\" (#1480)"
        errors=$((errors + 1))
    elif echo "$line" | grep -q "(unused)"; then
        echo "::error file=$NETPROTO::kNetProtocol $proto is marked (unused) in $RES: reserve the next free number instead (#1480)"
        errors=$((errors + 1))
    else
        sum ""; sum "**kNetProtocol $proto** reserved: \`$(echo "$line" | sed 's/  */ /g')\`"
    fi
fi
# (#1881) kSaveVersion the same way: a "save <n> | ..." line in the same file, not "(unused)". Two branches that bumped the
# save version independently would each write a different format under one number, and the migration chain would read
# one of them wrongly.
SAVEVER=src/engine/saveversion.h
savever=$(grep -oE 'kSaveVersion = [0-9]+' "$SAVEVER" 2>/dev/null | head -1 | grep -oE '[0-9]+$')
if [ -z "$savever" ]; then
    echo "::error file=$SAVEVER::static check 6: kSaveVersion not found"; errors=$((errors + 1))
elif [ -f "$RES" ]; then
    line=$(tr -d '\r' < "$RES" | grep -E "^save[[:space:]]+${savever}[[:space:]]*\|" | head -1)
    if [ -z "$line" ]; then
        echo "::error file=$SAVEVER::kSaveVersion $savever is not reserved in $RES: add a line \"save $savever | <branch> | #<issue> | <what>\" (#1881)"
        errors=$((errors + 1))
    elif echo "$line" | grep -q "(unused)"; then
        echo "::error file=$SAVEVER::kSaveVersion $savever is marked (unused) in $RES: reserve the next free number instead (#1881)"
        errors=$((errors + 1))
    else
        sum ""; sum "**kSaveVersion $savever** reserved: \`$(echo "$line" | sed 's/  */ /g')\`"
    fi
fi
# (#2009) kErosionBakeVersion: an "erosion <n> | ..." line when the branch CHANGES the version against BASE (a branch
# that keeps its base's version needs none, so ladder's older bake still passes). Two branches that bumped to one number
# would each bake different terrain under it, and a cache baked by one loads in the other.
if [ -n "$BASE" ] && [ -f "$RES" ]; then
    ero_ver() { tr -d '\r' | grep -oE 'kErosionBakeVersion = [0-9]+' | head -1 | grep -oE '[0-9]+$'; }
    ero_new=$(ero_ver < "$EROSION_H" 2>/dev/null)
    ero_old=$(git show "$BASE:$EROSION_H" 2>/dev/null | ero_ver)
    if [ -n "$ero_new" ] && [ "$ero_new" != "$ero_old" ]; then
        line=$(tr -d '\r' < "$RES" | grep -E "^erosion[[:space:]]+${ero_new}[[:space:]]*\|" | head -1)
        if [ -z "$line" ]; then
            echo "::error file=$EROSION_H::kErosionBakeVersion ${ero_old:-?} -> $ero_new is not reserved in $RES: add a line \"erosion $ero_new | <branch> | #<issue> | <what>\" (#2009)"
            errors=$((errors + 1))
        elif echo "$line" | grep -q "(unused)"; then
            echo "::error file=$EROSION_H::kErosionBakeVersion $ero_new is marked (unused) in $RES: reserve the next free number instead (#2009)"
            errors=$((errors + 1))
        else
            sum ""; sum "**kErosionBakeVersion $ero_new** reserved: \`$(echo "$line" | sed 's/  */ /g')\`"
        fi
    fi
fi

# ---- 7. (#2009) a number this branch reserves is not reserved for another issue on the other base or in an open PR ----
# (check 6 reads only this branch's copy). RESERVATION_CLASH=0 skips it (offline).
if [ "${RESERVATION_CLASH:-1}" != 0 ]; then
    rc_out=$(bash tools/ci/reservation_clash.sh "$BASE" 2>&1)
    rc_rc=$?
    printf '%s\n' "$rc_out" | grep -E '^(::|reservation_clash:)'
    sum ""; sum "**Reservation clashes (check 7):** $(printf '%s\n' "$rc_out" | grep '^reservation_clash:' | sed 's/^reservation_clash: //')"
    [ "$rc_rc" -ne 0 ] && errors=$((errors + 1))
fi

# ---- 8. (#2122) every known-failure line names its owner and an expiry, at most 7 days ahead ----
# The note starts "[<owner>, until <YYYY-MM-DD>]" (or "until never": a baseline by design). A line past its date is a
# warning (its owner fixes it or renews it, at most 7 days ahead); a missing or malformed tag, or a date further out, is an
# error. KNOWN_TODAY=<YYYY-MM-DD> overrides today (tests).
# The null-GL shards' list (#2123) follows the same rules.
for KF in tools/ci/known_failures.txt tools/ci/known_failures_nogl.txt; do
if [ -f "$KF" ]; then
    today="${KNOWN_TODAY:-$(date -u +%F)}"
    latest=$(date -u -d "$today +7 days" +%F 2>/dev/null || echo 9999-12-31)
    kf_out=$(tr -d '\r' < "$KF" | awk -F'|' -v today="$today" -v latest="$latest" -v kf="$KF" '
        /^[[:space:]]*(#|$)/ { next }
        {
            note = $3; sub(/^[[:space:]]+/, "", note)
            if (NF < 3 || note !~ /^\[(L&G|Water|UI\/UX|Progression|Space|Tools|Coordinator), until ([0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]|never)\]/) {
                printf "::error file=%s,line=%d::known failure line without \"[<owner>, until <YYYY-MM-DD>]\" at the start of its note (owners: L&G, Water, UI/UX, Progression, Space, Tools, Coordinator; #2122)\n", kf, NR; bad++; next
            }
            d = note; sub(/^\[[^,]*, until /, "", d); sub(/\].*/, "", d)
            o = note; sub(/^\[/, "", o); sub(/,.*/, "", o)
            if (d == "never") { never++; next }
            if (d > latest) { printf "::error file=%s,line=%d::known failure expiry %s is more than 7 days ahead (latest %s, #2122)\n", kf, NR, d, latest; bad++; next }
            if (d < today) { printf "::warning file=%s,line=%d::known failure expired %s (%s): fix it, drop it, or renew it (#2122)\n", kf, NR, d, o; expired++ }
            lines++
        }
        END { printf "known_lines: %d dated, %d never, %d expired, %d bad\n", lines, never, expired, bad; exit bad > 0 }')
    kf_rc=$?
    printf '%s\n' "$kf_out" | grep -E '^::'
    sum ""; sum "**Known failures (check 8, $(basename "$KF")):** $(printf '%s\n' "$kf_out" | grep '^known_lines:' | sed 's/^known_lines: //')"
    [ "$kf_rc" -ne 0 ] && errors=$((errors + 1))
fi
done

if [ -z "$BASE" ]; then sum ""; sum "_No base commit: only the NetMsg id check ran._"; fi
if [ "$errors" -gt 0 ]; then
    echo "static_checks: $errors error(s)"
    exit 1
fi
echo "static_checks: ok"
