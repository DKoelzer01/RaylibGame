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
#   --branch B         the branch whose builds count (default develop): the newest run whose package job succeeded and
#                      whose package VERSION names B (push and dispatch runs on B; on the integration branch also the
#                      schedule runs, #1966)
#   --integration B    the integration branch, whose batch runs are schedule events (default staging)
#   --force            fetch even when the newest run is the one already in DIR
#
# Exit codes: 0 updated or already current, 1 error, 3 no run with a successful package job and the artifact (expired after 5 days?).
# Needs: curl, jq, unzip (or bsdtar, or python3), sha256sum.
set -euo pipefail

REPO=DKoelzer01/RaylibGame
WORKFLOW=selfhosted.yml
ARTIFACT=raylibgame-server-linux-x86_64
BRANCH=develop
INTEGRATION=staging
TOKEN_FILE=""; DEST=""; FORCE=0; KEEP=5
while [ $# -gt 0 ]; do
    case "$1" in
        --token-file) TOKEN_FILE="$2"; shift ;;
        --dest) DEST="$2"; shift ;;
        --branch) BRANCH="$2"; shift ;;
        --integration) INTEGRATION="$2"; shift ;;
        --repo) REPO="$2"; shift ;;
        --workflow) WORKFLOW="$2"; shift ;;
        --artifact) ARTIFACT="$2"; shift ;;
        --force) FORCE=1 ;;
        -h|--help) sed -n '2,24p' "$0"; exit 0 ;;
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

# The newest finished run of the workflow that built the branch, whose package job succeeded (#1028: the package is its
# own job, so a failing scenario or regression test elsewhere in the run doesn't hold the server back; runs from before
# the split count when the whole run succeeded) and that still has the artifact. (#1966) Since the CI tiers (#1871) the
# integration branch has few push runs: its batch runs are schedule events (they fire from the default branch, so their
# head_branch is main) and by-hand dispatches. Candidates: push and dispatch runs on the branch, and on the integration
# branch also the schedule runs. A dispatch can build another branch (-f ref=...) and a schedule run's head is main, so
# the package's own VERSION ("<branch> <commit>", tools/build_linux.sh) decides: a package of another branch is skipped
# and remembered, so it isn't downloaded again every hour.
meta="$DEST/raylibgame-server-$safe_branch.json"
rejected="$DEST/.rejected-runs-$safe_branch.txt"
current="$(jq -r .run_id "$meta" 2>/dev/null || true)"
all="$(api "https://api.github.com/repos/$REPO/actions/workflows/$WORKFLOW/runs?status=completed&per_page=50")"
cands="$(echo "$all" | jq -r --arg b "$BRANCH" --arg ib "$INTEGRATION" \
    '.workflow_runs[] | select((.head_branch == $b and (.event == "push" or .event == "workflow_dispatch")) or ($b == $ib and .event == "schedule")) | .id')"
run_id=""; sha=""; when=""; tgz=""; n=0
for id in $cands; do
    if [ "$FORCE" = 0 ] && [ "$id" = "$current" ]; then
        echo "update-server: $BRANCH already at run $id ($(jq -r '.commit[0:7]' "$meta" 2>/dev/null))"
        exit 0
    fi
    grep -qx "$id" "$rejected" 2>/dev/null && continue
    concl="$(echo "$all" | jq -r --argjson i "$id" '.workflow_runs[] | select(.id == $i) | .conclusion')"
    if [ "$concl" != success ]; then
        pkg="$(api "https://api.github.com/repos/$REPO/actions/runs/$id/jobs" | jq -r '[.jobs[] | select(.name == "Linux server build, unit tests, package")][0].conclusion // empty')"
        [ "$pkg" = success ] || continue
    fi
    arts="$(api "https://api.github.com/repos/$REPO/actions/runs/$id/artifacts")"
    url="$(echo "$arts" | jq -r --arg n "$ARTIFACT" '[.artifacts[] | select(.name == $n and .expired == false)][0].archive_download_url // empty')"
    [ -n "$url" ] || continue
    n=$((n + 1)); [ $n -gt 5 ] && break # at most 5 downloads (~50 MB each) per call
    echo "update-server: checking run $id"
    x="$work/x$id"; mkdir "$x"
    api -o "$x/artifact.zip" "$url"
    if command -v unzip > /dev/null; then unzip -q "$x/artifact.zip" -d "$x"
    elif command -v bsdtar > /dev/null; then bsdtar -xf "$x/artifact.zip" -C "$x"
    else python3 -m zipfile -e "$x/artifact.zip" "$x"; fi   # deweyserver has python3 but no unzip
    t="$x/$ARTIFACT.tar.gz"
    [ -f "$t" ] || { echo "update-server: run $id: the artifact holds no $ARTIFACT.tar.gz" >&2; echo "$id" >> "$rejected"; continue; }
    # (#1966) one member listed, no pipe: "tar -tzf | grep -q" under pipefail failed on every good package (grep quits at
    # the first match, tar dies of SIGPIPE with 141), so no build was ever installed.
    tar -tzf "$t" raylibgame-server/game > /dev/null 2>&1 || { echo "update-server: run $id: the tarball has no raylibgame-server/game" >&2; echo "$id" >> "$rejected"; continue; }
    ver="$(tar -xzOf "$t" raylibgame-server/VERSION 2>/dev/null | tr -d '\r\n' || true)"
    if [ "${ver%% *}" != "$BRANCH" ]; then
        echo "update-server: run $id built '${ver:-no VERSION}', not $BRANCH: skipped"; echo "$id" >> "$rejected"; continue
    fi
    run_id="$id"; tgz="$t"
    sha="${ver#* }"
    when="$(echo "$all" | jq -r --argjson i "$id" '.workflow_runs[] | select(.id == $i) | .updated_at')"
    break
done
[ -n "$run_id" ] || { echo "update-server: no $WORKFLOW run that built $BRANCH with a successful package job and an unexpired $ARTIFACT"; exit 3; }
echo "update-server: installing run $run_id ($ver, $when)"
sum="$(sha256sum "$tgz" | cut -d' ' -f1)"

umask 022
versioned="$DEST/raylibgame-server-$safe_branch-$run_id.tar.gz"
# (#1966) chmod: cp keeps the umask-077 temp file's mode 600, which the web server can't read
cp "$tgz" "$versioned.part" && chmod 644 "$versioned.part" && mv "$versioned.part" "$versioned"
cp "$tgz" "$DEST/.raylibgame-server-$safe_branch.tar.gz.part" && chmod 644 "$DEST/.raylibgame-server-$safe_branch.tar.gz.part"
mv "$DEST/.raylibgame-server-$safe_branch.tar.gz.part" "$DEST/raylibgame-server-$safe_branch.tar.gz" # atomic: a download never sees half a file
jq -n --arg repo "$REPO" --arg branch "$BRANCH" --arg run "$run_id" --arg sha "$sha" --arg when "$when" --arg sum "$sum" \
    '{repo: $repo, branch: $branch, run_id: $run, commit: $sha, built: $when, sha256: $sum,
      run_url: ("https://github.com/" + $repo + "/actions/runs/" + $run)}' > "$meta.part"
mv "$meta.part" "$meta"
# Keep the newest $KEEP versioned builds of this branch.
ls -1t "$DEST"/raylibgame-server-"$safe_branch"-*.tar.gz 2>/dev/null | tail -n +$((KEEP + 1)) | xargs -r rm -f
echo "update-server: $BRANCH -> $DEST/raylibgame-server-$safe_branch.tar.gz (run $run_id, ${sha:0:7}, sha256 ${sum:0:12}...)"
echo "update-server: the game server is NOT restarted: reinstall it from the panel (or the deweyserver MCP) when nobody is playing"
