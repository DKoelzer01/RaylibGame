# RaylibGame on Pterodactyl (#326)

`egg-raylibgame.json` is a Pterodactyl egg for the dedicated server (`game.exe --server server.cfg`).
**Not deployed:** importing it on deweyserver's panel and creating a server needs the user's OK.

- **Runtime:** the build is Windows-only, so the egg runs it under Wine (`ghcr.io/parkervcp/yolks:wine_latest`) inside
  `xvfb-run`. The server still opens a hidden window for the OpenGL 4.3 context the world loaders need; in a container
  without a GPU that is Mesa's llvmpipe (software). **Untested under Wine/llvmpipe** - the first run on the server should
  check that the world generates (`host: listening on port ...` then `server: ready`) and how long the erosion bake takes
  on that CPU. A native Linux build would avoid Wine entirely (later work).
- **Install:** fetches `DOWNLOAD_URL` (a zip of `game.exe` + `resources/`) if given, otherwise upload the build through the
  panel's file manager. Writes `server.cfg` from the panel's variables.
- **Config:** the panel rewrites `server.cfg`'s `port`, `max_players`, `password`, `motd`, `autosave_minutes` and `seed`
  lines at every start (the `file` parser). The save is `server/world.sav`, the log `server/server.log`.
- **Console:** the panel console is the server's stdin: `help`, `status`, `list`, `save`, `say <text>`, `kick <name>`,
  `stop`. The panel's Stop button sends `stop` (the world is saved, every player hears "server closed").
- **Startup detection:** the panel marks the server running when it prints `server: ready`.

Local test without the panel:

```
game.exe --server server.cfg            # writes a default server.cfg the first time
game.exe --server server.cfg --sandbox  # a flat test world (no erosion bake: up in seconds)
```
