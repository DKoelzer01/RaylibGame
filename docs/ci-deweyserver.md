# Server builds on deweyserver (#1016)

The self-hosted workflow (`.github/workflows/selfhosted.yml`) builds the Linux server tarball on every push to
`develop` / `m84-integration` and uploads it as the workflow artifact `raylibgame-server-linux-x86_64` (kept 5 days).
The runner's container has no LAN access, so it cannot touch the game server. `tools/ci/deweyserver-update-server.sh`
runs on the host instead and makes the newest build available to the Pterodactyl egg
(`tools/pterodactyl/egg-raylibgame-linux.json`, variable `DOWNLOAD_URL`). **Nothing here is deployed yet**: each step
below needs the user's OK.

## One-time setup (as `claude`)

1. A GitHub token that can only read this repo's Actions: a fine-grained token, repository `DKoelzer01/RaylibGame`,
   permission **Actions: Read-only**, nothing else. Artifacts download only with a token, even from a public repo.
   ```
   install -m 600 /dev/null ~/.config/raylibgame/gh-actions-token   # then paste the token into it
   ```
2. A folder a web server serves on the LAN / tailnet only (never the internet), e.g. `~/srv/raylibgame`, served by the
   existing reverse proxy or a small static server bound to the LAN address:
   ```
   python3 -m http.server 8787 --bind 192.168.4.29 --directory ~/srv/raylibgame   # e.g. as a systemd user service
   ```
3. A first run (needs `curl`, `jq`, `unzip`):
   ```
   tools/ci/deweyserver-update-server.sh --token-file ~/.config/raylibgame/gh-actions-token \
       --dest ~/srv/raylibgame --branch develop
   ```
   It writes `raylibgame-server-develop.tar.gz` (the newest, replaced atomically), the last 5 builds as
   `raylibgame-server-develop-<run>.tar.gz`, and `raylibgame-server-develop.json` (run id, commit, sha256). A re-run does
   nothing until a newer run whose package job succeeded exists (exit 0; #1028: failing tests elsewhere in the run
   don't hold a build back); exit 3 means no such run with an unexpired artifact.
4. The egg's `DOWNLOAD_URL` = `http://192.168.4.29:8787/raylibgame-server-develop.tar.gz` (a pinned
   `...-develop-<run>.tar.gz` holds a server on one build).

## Keeping it current

The script never restarts or reinstalls the game server: an update means a reinstall, which kicks players. Run it by
hand, or from a user timer that only fetches (`systemctl --user enable --now raylibgame-fetch.timer`):

```
# ~/.config/systemd/user/raylibgame-fetch.service
[Service]
Type=oneshot
ExecStart=%h/RaylibGame/tools/ci/deweyserver-update-server.sh --token-file %h/.config/raylibgame/gh-actions-token --dest %h/srv/raylibgame --branch develop
SuccessExitStatus=3

# ~/.config/systemd/user/raylibgame-fetch.timer
[Timer]
OnCalendar=hourly
Persistent=true
[Install]
WantedBy=timers.target
```

Then reinstall the server from the panel (or the deweyserver MCP's game-server tools) when nobody is playing: the
reinstall keeps `server/` (the world) and `server.cfg`. Run the script from a checkout that is not the runner's work
folder (the runner cleans its own).
