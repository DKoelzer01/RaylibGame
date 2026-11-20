#!/bin/bash
# CI (#1016): the FetchContent sources (raylib, googletest) cloned once into a cacheable folder, so a build passes
# -DFETCHCONTENT_SOURCE_DIR_<NAME>=DIR/<name> and downloads nothing. Repositories and tags are read from CMakeLists.txt,
# so a version bump there changes the cache key and the clone.
#
#   tools/ci/fetch_deps.sh --key        prints "raylib@5.5 googletest@v1.15.2" (the cache key's input)
#   tools/ci/fetch_deps.sh DIR          clones each one missing or at another tag into DIR/<name>
#   tools/ci/fetch_deps.sh --cmake DIR  prints the -DFETCHCONTENT_SOURCE_DIR_... arguments for DIR
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
NAMES="raylib googletest" # tracy is optional (ENABLE_TRACY) and not built in CI

# "name repo tag" for the first FetchContent_Declare of each name.
decls() {
    tr -d '\r' < "$REPO/CMakeLists.txt" | awk '
        /FetchContent_Declare\(/ { s = $0; sub(/.*FetchContent_Declare\(/, "", s); split(s, w, /[ \t)]/); cur = w[1]; next }
        cur != "" && /GIT_REPOSITORY/ { if (!(cur in repo)) repo[cur] = $2 }
        cur != "" && /GIT_TAG/ { if (!(cur in tag)) tag[cur] = $2 }
        /\)/ { cur = "" }
        END { for (n in repo) print n, repo[n], tag[n] }'
}

case "${1:-}" in
    --key)
        out=""
        for n in $NAMES; do t="$(decls | awk -v n="$n" '$1 == n { print $3 }')"; out="$out $n@$t"; done
        echo "${out# }" ;;
    --cmake)
        DIR="$2"
        for n in $NAMES; do printf -- '-DFETCHCONTENT_SOURCE_DIR_%s=%s/%s ' "$(echo "$n" | tr a-z A-Z)" "$DIR" "$n"; done
        echo ;;
    ""|-h|--help) sed -n '2,9p' "$0" ;;
    *)
        DIR="$1"; mkdir -p "$DIR"
        for n in $NAMES; do
            read -r _ url tag < <(decls | awk -v n="$n" '$1 == n')
            [ -n "${url:-}" ] && [ -n "${tag:-}" ] || { echo "fetch_deps: no FetchContent_Declare($n ...) with GIT_REPOSITORY/GIT_TAG in CMakeLists.txt"; exit 1; }
            if [ -f "$DIR/$n/.ci-tag" ] && [ "$(cat "$DIR/$n/.ci-tag")" = "$url $tag" ]; then
                echo "fetch_deps: $n $tag cached"
                continue
            fi
            rm -rf "$DIR/$n"
            echo "fetch_deps: cloning $n $tag"
            git -c advice.detachedHead=false clone -q --depth 1 --branch "$tag" "$url" "$DIR/$n"
            rm -rf "$DIR/$n/.git"
            echo "$url $tag" > "$DIR/$n/.ci-tag"
        done ;;
esac
