# Playing together (#341)

Multiplayer is by IP address: one player (or a dedicated server) hosts a world, the others join it.

## Host from the game

Main menu, **Multiplayer**, **Host a world**:

- **World:** one of your saves, or a new world (a new world is built first: a few minutes the first time).
- **Port:** UDP, default 27015. Players outside your network need it forwarded on your router to this PC.
- **Max players** and an optional **password**.
- **Your name** and **suit colour**: how the others see you.

**Host** starts a dedicated server in the background (this same `game.exe --server`, no window) and joins it. When you
leave (Save & Quit, or Exit) the server saves the world and stops. The world is a normal save in `saves\`. Next to it:
- `<world>.sav.server.cfg`: the server's settings;
- `<world>.sav.server.log`: its log, which says why if it didn't start.

If the port is already in use (another server, or another program), the Host page says so: pick another.

## Join

Main menu, **Multiplayer**, **Join a server**:
- the host's **address**, as `ip` or `ip:port`;
- your **name** and **suit colour**;
- the **password**, if the server has one.

Your last servers are listed under **Recent**. The first join fetches the world's changes (digs, buildings) from the
server; they're cached in `saves\client_<name>.sav*`, so rejoining fetches only what changed since.

In game, hold **J** to see who is playing (the pause menu shows it too).

## A dedicated server (Windows)

To keep a world running without playing on it, run `run_server.bat` next to `game.exe`, or
`game.exe --server server.cfg` from a console.

- The first start writes `server.cfg` with the defaults: port 27015, 8 players, no password, the world in
  `server\world.sav`, a log in `server\server.log`, autosave every 30 minutes. Edit it and start again.
- Console commands: `help`, `status`, `list`, `save`, `say <text>`, `kick <name> [reason]`, `stop`. `stop` saves and
  ends the server; closing the window ends it without a final save (the last autosave stands).
- The server needs the game's `resources\` folder next to `game.exe`, and opens a hidden window: Windows still needs a
  GPU driver with OpenGL 4.3 to generate the world.

A Linux build and a Pterodactyl egg are separate work (`tools/pterodactyl/`).
