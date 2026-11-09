#include "gamesave.h"
#include "inventory.h"
#include "player.h"
#include "world.h"

#include <cstring>

void GameSave_Register() {
    saveExtrasHook = [](FILE* f) {
        fprintf(f, "player %.6f %.6f %.6f %f %f %f %d %f\n", player.position.x, player.position.y, player.position.z,
                player.forward.x, player.forward.y, player.forward.z, player.flyMode ? 1 : 0, player.jetpackFuel);
        for (const auto& kv : inventory.all()) fprintf(f, "item %s %d\n", kv.first.c_str(), kv.second);
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
            } else if (strcmp(tok, "item") == 0) {
                char id[64] = {};
                int n = 0;
                if (fscanf(f, " %63s %d", id, &n) == 2) inventory.add(id, n);
            }
        }
    };
}
