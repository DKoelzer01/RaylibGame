#!/bin/bash
# deweyserver host side of the server package (#1016): fetches the newest Linux server tarball that the self-hosted
# workflow (.github/workflows/selfhosted.yml) built on a branch and puts it where the Pterodactyl egg's DOWNLOAD_URL can
# fetch it. The runner container has no LAN access, so this runs on the host itself, as the `claude` user, by hand or
# from a systemd timer (see docs/ci-deweyserver.md). It never restarts or reinstalls the game server: that stays a
# deliberate step (panel or the deweyserver MCP), so nobody is kicked mid-session.
#
#   deweyserver-update-server.sh --token-file FILE --dest DIR [--branch develop] [--repo DKoelzer01/RaylibGame]
#                                [--workflow selfhosted.yml] [--artifact raylibgame-server-linux-x86_64] [--force]
#
#   --token-file FILE  a GitHub token (fine-grained, this repo, "Actions: read" only) in a file of mode 600 owned by
#                      the caller; artifacts download only with a token, even from a public repo
#   --dest DIR         the folder a web server serves (LAN / tailnet only): gets
#                        raylibgame-server-<branch>-<run>.tar.gz   each fetched build (the newest 5 are kept)
#                        raylibgame-server-<branch>.tar.gz         the newest, replaced atomically (the egg's DOWNLOAD_URL)
#                        raylibgame-server-<branch>.json           run id, commit, time, sha256
#   --branch B         the branch whose push runs count (default develop): the newest one whose package job succeeded
#   --force            fetch even when the newest run is the one already in DIR
#
# Exit codes: 0 updated or already current, 1 error, 3 no run with a successful package job and the artifact (expired after 5 days?).
# Needs: curl, jq, unzip (or bsdtar, or python3), sha256sum.
set -euo pipefail

REPO=DKoelzer01/RaylibGame
WORKFLOW=selfhosted.yml
ARTIFACT=raylibgame-server-linux-x86_64
BRANCH=develop
TOKEN_FILE=""; DEST=""; FORCE=0; KEEP=5
while [ $# -gt 0 ]; do
    case "$1" in
        --token-file) TOKEN_FILE="$2"; shift ;;
        --dest) DEST="$2"; shift ;;
        --branch) BRANCH="$2"; shift ;;
        --repo) REPO="$2"; shift ;;
        --workflow) WORKFLOW="$2"; shift ;;
        --artifact) ARTIFACT="$2"; shift ;;
        --force) FORCE=1 ;;
        -h|--help) sed -n '2,25p' "$0"; exit 0 ;;
        *) echo "update-server: unknown argument $1" >&2; exit 1 ;;
    esac
    shift
done
[ -n "$TOKEN_FILE" ] && [ -n "$DEST" ] || { echo "update-server: --token-file and --dest are required (--help)" >&2; exit 1; }
for t in curl jq sha256sum; do command -v "$t" > /dev/null || { echo "update-server: '$t' missing" >&2; exit 1; }; done
[ -f "$TOKEN_FILE" ] || { echo "update-server: no token file $TOKEN_FILE" >&2; exit 1; }
perm="$(stat -c '%a' "$TOKEN_FILE")"; owner="$(stat -c '%U' "$TOKEN_FILE")"
if [ "$perm" != 600 ] && [ "$perm" != 400 ]; then echo "update-server: $TOKEN_FILE has mode $perm: chmod 600 it" >&2; exit 1; fi
[ "$owner" = "$(id -un)" ] || { echo "update-server: $TOKEN_FILE belongs to $owner, not $(id -un)" >&2; exit 1; }
case "$BRANCH" in *[!A-Za-z0-9._/-]*|"") echo "update-server: odd branch name '$BRANCH'" >&2; exit 1 ;; esac
safe_branch="${BRANCH//\//-}"

mkdir -p "$DEST"
DEST="$(cd "$DEST" && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
# The token goes to curl through a header file, never on the command line (ps would show it).
umask 077
printf 'Authorization: Bearer %s\n' "$(tr -d '\r\n' < "$TOKEN_FILE")" > "$work/auth"
api() { curl -fsSL --retry 3 -H @"$work/auth" -H "Accept: application/vnd.github+json" -H "X-GitHub-Api-Version: 2022-11-28" "$@"; }

# The newest finished push run of the workflow on the branch whose package job succeeded (#1028: the package is its own
# job, so a failing scenario or regression test elsewhere in the run doesn't hold the server back; runs from before the
# split count when the whole run succeeded) and that still has the artifact.
runs="$(api "https://api.github.com/repos/$REPO/actions/workflows/$WORKFLOW/runs?branch=$BRANCH&event=push&status=completed&per_page=10")"
run_id=""; art_url=""; sha=""; when=""
for id in $(echo "$runs" | jq -r '.workflow_runs[].id'); do
    concl="$(echo "$runs" | jq -r --argjson i "$id" '.workflow_runs[] | select(.id == $i) | .conclusion')"
    if [ "$concl" != success ]; then
        pkg="$(api "https://api.github.com/repos/$REPO/actions/runs/$id/jobs" | jq -r '[.jobs[] | select(.name == "Linux server build, unit tests, package")][0].conclusion // empty')"
        [ "$pkg" = success ] || continue
    fi
    arts="$(api "https://api.github.com/repos/$REPO/actions/runs/$id/artifacts")"
    url="$(echo "$arts" | jq -r --arg n "$ARTIFACT" '[.artifacts[] | select(.name == $n and .expired == false)][0].archive_download_url // empty')"
    if [ -n "$url" ]; then
        run_id="$id"; art_url="$url"
        sha="$(echo "$runs" | jq -r --argjson i "$id" '.workflow_runs[] | select(.id == $i) | .head_sha')"
        when="$(echo "$runs" | jq -r --argjson i "$id" '.workflow_runs[] | select(.id == $i) | .updated_at')"
        break
    fi
done
[ -n "$run_id" ] || { echo "update-server: no $WORKFLOW push run on $BRANCH with a successful package job and an unexpired $ARTIFACT"; exit 3; }

meta="$DEST/raylibgame-server-$safe_branch.json"
if [ "$FORCE" = 0 ] && [ -f "$meta" ] && [ "$(jq -r .run_id "$meta" 2>/dev/null)" = "$run_id" ]; then
    echo "update-server: $BRANCH already at run $run_id (${sha:0:7})"
    exit 0
fi

echo "update-server: fetching run $run_id (${sha:0:7}, $when)"
api -o "$work/artifact.zip" "$art_url"
mkdir "$work/x"
if command -v unzip > /dev/null; then unzip -q "$work/artifact.zip" -d "$work/x"
elif command -v bsdtar > /dev/null; then bsdtar -xf "$work/artifact.zip" -C "$work/x"
else python3 -m zipfile -e "$work/artifact.zip" "$work/x"; fi   # deweyserver has python3 but no unzip
tgz="$work/x/$ARTIFACT.tar.gz"
[ -f "$tgz" ] || { echo "update-server: the artifact holds no $ARTIFACT.tar.gz" >&2; exit 1; }
tar -tzf "$tgz" | grep -q '^raylibgame-server/game$' || { echo "update-server: the tarball has no raylibgame-server/game" >&2; exit 1; }
sum="$(sha256sum "$tgz" | cut -d' ' -f1)"

umask 022
versioned="$DEST/raylibgame-server-$safe_branch-$run_id.tar.gz"
cp "$tgz" "$versioned.part" && mv "$versioned.part" "$versioned"
cp "$tgz" "$DEST/.raylibgame-server-$safe_branch.tar.gz.part"
mv "$DEST/.raylibgame-server-$safe_branch.tar.gz.part" "$DEST/raylibgame-server-$safe_branch.tar.gz" # atomic: a download never sees half a file
jq -n --arg repo "$REPO" --arg branch "$BRANCH" --arg run "$run_id" --arg sha "$sha" --arg when "$when" --arg sum "$sum" \
    '{repo: $repo, branch: $branch, run_id: $run, commit: $sha, built: $when, sha256: $sum,
      run_url: ("https://github.com/" + $repo + "/actions/runs/" + $run)}' > "$meta.part"
mv "$meta.part" "$meta"
# Keep the newest $KEEP versioned builds of this branch.
ls -1t "$DEST"/raylibgame-server-"$safe_branch"-*.tar.gz 2>/dev/null | tail -n +$((KEEP + 1)) | xargs -r rm -f
echo "update-server: $BRANCH -> $DEST/raylibgame-server-$safe_branch.tar.gz (run $run_id, ${sha:0:7}, sha256 ${sum:0:12}...)"
echo "update-server: the game server is NOT restarted: reinstall it from the panel (or the deweyserver MCP) when nobody is playing"
