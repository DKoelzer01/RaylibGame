# Playing together (#341)

Multiplayer is by IP address. There are two ways to run a world others join:

| | Host game (from the menu) | Dedicated server |
| --- | --- | --- |
| Started by | a player, from the main menu | `run_server.bat`, or `game.exe --server <cfg>` |
| Runs | while that player plays | on its own, with or without players, until stopped |
| When the host leaves | the world saves and closes for everyone | nothing changes |
| Good for | playing with friends now | a world that is always there |

## Host game

Main menu, **Multiplayer**, **Host game**. Others can join while you play; the world closes when you leave.

- **World:** one of your saves, or a new world (a new world is built first: a few minutes the first time).
- **Port:** UDP, default 27015. Players outside your network need it forwarded on your router to this PC.
- **Max players** and an optional **password**.
- **Your name** and **suit colour**: how the others see you.

**Host** opens your world to others and puts you in it. When you leave (Save & Quit, or Exit), the world is saved and
closes, and everyone else is disconnected. The world is a normal save in `saves\`, which you can also play alone. Next to it:
- `<world>.sav.server.cfg`: the hosting settings;
- `<world>.sav.server.log`: the log, which says why if the world didn't open.

(Under the hood the game runs the world in a background copy of itself, the same code as a dedicated server, and joins it.
It is tied to your game: it closes when you leave, and if the game crashes it closes too.)

If the port is already in use (another game, or another program), the Host page says so: pick another.

## Join game

Main menu, **Multiplayer**, **Join game**:
- the **address** of the host or server, as `ip` or `ip:port`;
- your **name** and **suit colour**;
- the **password**, if the world has one.

Your last addresses are listed under **Recent**. The first join fetches the world's changes (digs, buildings); they're
cached in `saves\client_<name>.sav*`, so rejoining fetches only what changed since.

In game, hold **J** to see who is playing (the pause menu shows it too).

## Dedicated server (Windows)

A dedicated server runs the world without a player and keeps it going until you stop it: players come and go, and the
world stays up. Run `run_server.bat` next to `game.exe`, or `game.exe --server server.cfg` from a console.

- The first start writes `server.cfg` with the defaults: port 27015, 8 players, no password, the world in
  `server\world.sav`, a log in `server\server.log`, autosave every 30 minutes. Edit it and start again.
- Console commands: `help`, `status`, `list`, `save`, `say <text>`, `kick <name> [reason]`, `stop`. `stop` saves and
  ends the server; closing the window ends it without a final save (the last autosave stands).
- The server needs the game's `resources\` folder next to `game.exe`, and opens a hidden window: Windows still needs a
  GPU driver with OpenGL 4.3 to generate the world.

A Linux build and a Pterodactyl egg are separate work (`tools/pterodactyl/`).
