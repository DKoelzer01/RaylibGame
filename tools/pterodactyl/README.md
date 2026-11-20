# RaylibGame on Pterodactyl (#326, #792)

`egg-raylibgame-linux.json` is a Pterodactyl egg for the native Linux dedicated server (`./game --server server.cfg`).
**Not deployed:** importing it on deweyserver's panel and creating a server needs the user's OK.

It replaces the earlier Wine egg (`egg-raylibgame.json`, removed: Windows build under Wine + xvfb + llvmpipe, never
tested). Windows self-hosters run `game.exe --server server.cfg` directly; no panel needed.

## What the egg does

- **Image:** `ghcr.io/parkervcp/yolks:debian` (Debian 12). The server needs nothing else: libstdc++ is linked in, and
  the runtime folder brings libglvnd's `libEGL.so.1` (`lib/`) for the GPU path.
- **Install** (in `debian:bookworm`): either fetches `DOWNLOAD_URL` (the `raylibgame-server-linux-x86_64.tar.gz` that
  `tools/build_linux.sh` makes) or, when it is empty, clones `GIT_REPO` at `GIT_REF` and builds it (about 10 minutes:
  GameNetworkingSockets + the game). Writes `server.cfg` from the panel's variables the first time. Reinstalling updates
  the build and keeps `server/` (the world) and `server.cfg`.
- **Config:** the panel rewrites `server.cfg`'s `port`, `max_players`, `password`, `motd`, `autosave_minutes` and `seed`
  lines at every start (the `file` parser). The save is `server/world.sav`, the log `server/server.log`, the admin lists
  `server/admin.txt` and `server/audit.log`.
- **Console:** the panel console is the server's stdin: `help`, `status`, `list`, `save`, `say <text>`,
  `kick <name> [reason]`, `op`, `ban`, `whitelist`..., `stop`. The panel's Stop button sends `stop` (the world is saved,
  every player hears "server closed"); a kill (SIGTERM) takes the same path.
- **Startup detection:** `server: ready`.
- **Port:** the allocation's port is the game port (UDP).

## GPU (the NVIDIA runtime)

The server has no window and needs no X11 or Wayland. At start it looks for an OpenGL 4.3 context through EGL:

1. **NVIDIA GPU in the container** (`HEADLESS: EGL 1.5 device 0 (...): NVIDIA GeForce RTX 2080 SUPER/PCIe/SSE2, OpenGL
   4.6...` in the console): the GPU paths are on (GPU fluid, GPU terrain density).
2. **No GPU** (`HEADLESS: null GL (no GPU context: ...)`): every GL call is a stub; the fluid and terrain run on their CPU
   paths. Works in any container, costs CPU.

For 1 the container must be started by the **nvidia runtime** with the GPU visible. The egg sets
`NVIDIA_VISIBLE_DEVICES=all` and `NVIDIA_DRIVER_CAPABILITIES=graphics,compute,utility` as server variables, but Wings
has no per-server runtime setting, so Docker itself must use it. On deweyserver (nvidia-container-toolkit is installed,
runtime `nvidia` exists) that means, as root, **with the user's OK**:

```
# /etc/docker/daemon.json: add (keep the existing "runtimes" entry)
"default-runtime": "nvidia"
# then
sudo systemctl restart docker   # restarts every container on the host (Wolf, media, MCP...); Wings reconnects
```

The nvidia runtime only acts on containers that set `NVIDIA_VISIBLE_DEVICES`, so the other stacks are unaffected
(Wolf already asks for the GPU itself). Without this change the egg still works, on the CPU paths.

The GPU is shared with Wolf: the server uses a few hundred MB of VRAM and a little GPU time per tick.

Switches (set as extra variables or in the startup line): `RAYLIBGAME_NO_GL=1` forces the CPU paths,
`RAYLIBGAME_EGL_DEVICE=N` picks a GPU, `RAYLIBGAME_EGL_SOFTWARE=1` accepts Mesa's llvmpipe instead of the null GL.

## Local test without the panel

```
tools/build_linux.sh                       # on a Linux host with Docker -> build/linux/raylibgame-server-linux-x86_64.tar.gz
tar xzf build/linux/raylibgame-server-linux-x86_64.tar.gz && cd raylibgame-server
./game --server server.cfg                 # writes a default server.cfg the first time
./game --server server.cfg --sandbox       # a flat test world (no erosion bake: up in seconds)
docker run --rm -it --gpus all -e NVIDIA_DRIVER_CAPABILITIES=graphics,compute,utility -p 27015:27015/udp \
    -v "$PWD":/srv -w /srv ghcr.io/parkervcp/yolks:debian ./game --server server.cfg   # as in the panel, with the GPU
```
