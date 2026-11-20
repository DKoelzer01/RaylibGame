# Bug report relay (#169)

The game's pause menu has **Report a bug**: title, description, tags (they become issue labels when the repo has them;
others are listed in the body), severity, and attachments (the screenshot taken as the game paused, `log.txt`, the last F10
capture, optionally the save). The game posts it with `curl.exe` to this relay, which:

1. stores the attachments in the Claude media store, `bug-reports/<date>/<HHMMSS>-<title>/` (Samba:
   `\\DEWEYSERVER\files\claude-media\bug-reports\...`; in a Claude session: `media_list` / `media_get`);
2. opens a GitHub issue **as the user**, labelled `bug report` (+ `bug` and the chosen tags that exist), linking the files.
   Triage is done later by the user and Claude.

User decisions (2026-11-15): tailnet-only for now (public later); posts as the user with their token; the attachments
above; labels from the in-game options.

## Deploy (on deweyserver, by the user: the token is yours)

1. Create a **fine-grained token**: GitHub, Settings, Developer settings, Fine-grained tokens. Repository access: only
   `DKoelzer01/RaylibGame`. Permissions: **Issues: read and write** (nothing else).
2. Put it on the server, readable only by you:
   ```sh
   install -d -m 700 ~/mcp-stack/secrets
   printf '%s' 'github_pat_...' > ~/mcp-stack/secrets/github_token && chmod 600 ~/mcp-stack/secrets/github_token
   ```
3. Add the service to `~/mcp-stack/docker-compose.yml` (next to the MCP server; same media mount):
   ```yaml
     bugrelay:
       image: python:3.12-slim
       restart: unless-stopped
       command: ["python", "/app/relay.py"]
       volumes:
         - ./bugrelay/relay.py:/app/relay.py:ro
         - ./secrets/github_token:/run/secrets/github_token:ro
         - /hoststorage/files/claude-media/bug-reports:/srv/claude-media/bug-reports
       ports:
         - "8765:8765"        # LAN + tailnet; drop this and use `tailscale serve` for tailnet-only
   ```
   and copy `tools/bugrelay/relay.py` to `~/mcp-stack/bugrelay/relay.py`. Then `docker compose up -d bugrelay`.
4. Check it: `curl http://deweyserver:8765/health`, which should show `"token": true`.
5. The game's default relay is `http://deweyserver:8765/report`. Another address: `RAYLIBGAME_BUG_RELAY=...`.

## Test without posting

`BUGRELAY_DRY_RUN=1 BUGRELAY_MEDIA=/tmp/bugrelay python3 relay.py` stores the report and prints the issue it would open.
The game's `--test-bug-report` sends a scripted report (with `RAYLIBGAME_BUG_RELAY` pointing at the test relay).
