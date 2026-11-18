@echo off
rem A dedicated server for RaylibGame (#326, #341): the same game.exe, headless, with a console.
rem The first start writes server.cfg with the defaults (port 27015, 8 players, no password, saves in server\);
rem edit it and start again. Console commands: help, status, list, save, say <text>, kick <name> [reason], stop.
rem Type "stop" to save the world and end the server (closing this window ends it without a final save:
rem the last autosave, every 30 minutes, stands). See docs\hosting.md.
cd /d "%~dp0"
game.exe --server server.cfg
pause
