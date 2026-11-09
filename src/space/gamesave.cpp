#include "gamesave.h"
#include "inventory.h"
#include "player.h"
#include "world.h"
#include "ship.h"
#include "survival.h"

#include <cstring>

void GameSave_Register() {
    saveExtrasHook = [](FILE* f) {
        fprintf(f, "player %.6f %.6f %.6f %f %f %f %d %f\n", player.position.x, player.position.y, player.position.z,
                player.forward.x, player.forward.y, player.forward.z, player.flyMode ? 1 : 0, player.jetpackFuel);
        fprintf(f, "player_up %f %f %f\n", player.up.x, player.up.y, player.up.z);
        for (const auto& kv : inventory.all()) fprintf(f, "item %s %d\n", kv.first.c_str(), kv.second);
        fprintf(f, "ship %.6f %.6f %.6f %f %f %f %f %f %d %d\n", ship.position.x, ship.position.y, ship.position.z, ship.orientation.x,
                ship.orientation.y, ship.orientation.z, ship.orientation.w, ship.energy, ship.piloted ? 1 : 0, ship.spawned ? 1 : 0);
        fprintf(f, "survival %f %f %f %d\n", survival.oxygen, survival.power, survival.health, survival.deaths);
    };
    loadExtrasHook = [](FILE* f) {
        inventory.clear();
        char tok[64] = {};
        while (fscanf(f, " %63s", tok) == 1 && strcmp(tok, "end_extras") != 0) {
            if (strcmp(tok, "player") == 0) {
                int fly = 0;
                fscanf(f, " %lf %lf %lf %f %f %f %d %f", &player.position.x, &player.position.y, &player.position.z,
                       &player.forward.x, &player.forward.y, &player.forward.z, &fly, &player.jetpackFuel);
                player.flyMode = fly != 0;
                player.velocity = { 0, 0, 0 };
                player.grounded = false;
                player.up = { 0, 0, 0 }; // set by "player_up"; older saves: rebuilt from the body at the next update
            } else if (strcmp(tok, "player_up") == 0) {
                fscanf(f, " %f %f %f", &player.up.x, &player.up.y, &player.up.z);
            } else if (strcmp(tok, "ship") == 0) {
                int piloted = 0, spawned = 0;
                fscanf(f, " %lf %lf %lf %f %f %f %f %f %d %d", &ship.position.x, &ship.position.y, &ship.position.z, &ship.orientation.x, &ship.orientation.y,
                       &ship.orientation.z, &ship.orientation.w, &ship.energy, &piloted, &spawned);
                ship.piloted = piloted != 0;
                ship.spawned = spawned != 0;
                ship.velocity = { 0, 0, 0 };
            } else if (strcmp(tok, "survival") == 0) {
                fscanf(f, " %f %f %f %d", &survival.oxygen, &survival.power, &survival.health, &survival.deaths);
                survival.dead = false;
            } else if (strcmp(tok, "item") == 0) {
                char id[64] = {};
                int n = 0;
                if (fscanf(f, " %63s %d", id, &n) == 2) inventory.add(id, n);
            }
        }
    };
}
