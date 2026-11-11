#include "gamesave.h"
#include "inventory.h"
#include "player.h"
#include "world.h"
#include "ship.h"
#include "survival.h"
#include "hotbar.h"
#include "rocket.h"
#include "terraintool.h"
#include "navigation.h"
#include "farming.h"
#include "animals.h"
#include "suit.h"
#include "toolwear.h"

#include <algorithm>
#include <cstring>

void GameSave_Register() {
    saveExtrasHook = [](FILE* f) {
        fprintf(f, "player %.6f %.6f %.6f %f %f %f %d %f\n", player.position.x, player.position.y, player.position.z,
                player.forward.x, player.forward.y, player.forward.z, player.flyMode ? 1 : 0, player.jetpackFuel);
        fprintf(f, "player_up %f %f %f\n", player.up.x, player.up.y, player.up.z);
        fprintf(f, "player_vel %f %f %f\n", player.velocity.x, player.velocity.y, player.velocity.z);
        for (int i = 0; i < Inventory::kSlots; ++i) {
            const InvSlot& s = inventory.grid()[i];
            if (!s.empty()) fprintf(f, "slot %d %s %lld\n", i, s.id.c_str(), (long long)s.centi);
        }
        QuestBook_Save(f);
        fprintf(f, "creative %d\n", terrainTool.creative ? 1 : 0);
        Hotbar_Save(f);
        Nav_Save(f);
        Farm_Save(f);
        Animals_Save(f);
        Suit_Save(f);
        ToolWear_Save(f);
        Rocket_Save(f);
        fprintf(f, "survival %f %f %f %d\n", survival.oxygen, survival.power, survival.health, survival.deaths);
        fprintf(f, "food %f\n", survival.food);
    };
    loadExtrasHook = [](FILE* f) {
        inventory.clear();
        terrainTool.creative = false;
        Hotbar_Reset();
        Nav_Reset();
        QuestBook_Reset();
        Farm_Reset();
        Animals_Reset();
        Suit_Reset();
        ToolWear_Reset();
        char tok[64] = {};
        while (fscanf(f, " %63s", tok) == 1 && strcmp(tok, "end_extras") != 0) {
            if (QuestBook_LoadLine(tok, f)) continue;
            if (Farm_LoadLine(tok, f)) continue;
            if (Animals_LoadLine(tok, f)) continue;
            if (Suit_LoadLine(tok, f)) continue;
            if (ToolWear_LoadLine(tok, f)) continue;
            if (strcmp(tok, "food") == 0) { float v = 100.0f; if (fscanf(f, " %f", &v) == 1) survival.food = std::clamp(v, 0.0f, 100.0f); continue; }
            if (strcmp(tok, "creative") == 0) { int c = 0; if (fscanf(f, " %d", &c) == 1) terrainTool.creative = c != 0; continue; }
            if (Hotbar_LoadLine(tok, f)) continue;
            if (Nav_LoadLine(tok, f)) continue;
            if (Rocket_LoadLine(tok, f)) continue;
            if (strcmp(tok, "player") == 0) {
                int fly = 0;
                fscanf(f, " %lf %lf %lf %f %f %f %d %f", &player.position.x, &player.position.y, &player.position.z,
                       &player.forward.x, &player.forward.y, &player.forward.z, &fly, &player.jetpackFuel);
                player.flyMode = fly != 0;
                survival.lastJetpackFuel = player.jetpackFuel; // a load is not a recharge
                player.velocity = { 0, 0, 0 };
                player.grounded = false;
                player.up = { 0, 0, 0 }; // set by "player_up"; older saves: rebuilt from the body at the next update
            } else if (strcmp(tok, "player_up") == 0) {
                fscanf(f, " %f %f %f", &player.up.x, &player.up.y, &player.up.z);
            } else if (strcmp(tok, "player_vel") == 0) { // older saves: stays zero
                fscanf(f, " %f %f %f", &player.velocity.x, &player.velocity.y, &player.velocity.z);
            } else if (strcmp(tok, "ship") == 0) {
                int piloted = 0, spawned = 0;
                fscanf(f, " %lf %lf %lf %f %f %f %f %f %d %d", &ship.position.x, &ship.position.y, &ship.position.z, &ship.orientation.x, &ship.orientation.y,
                       &ship.orientation.z, &ship.orientation.w, &ship.energy, &piloted, &spawned);
                ship.piloted = false; // old saves: the placeholder ship is gone (its line is read and dropped)
                ship.spawned = false; (void)piloted; (void)spawned;
                ship.velocity = { 0, 0, 0 };
            } else if (strcmp(tok, "ship_vel") == 0) { // older saves: stays zero
                fscanf(f, " %f %f %f", &ship.velocity.x, &ship.velocity.y, &ship.velocity.z);
            } else if (strcmp(tok, "survival") == 0) {
                fscanf(f, " %f %f %f %d", &survival.oxygen, &survival.power, &survival.health, &survival.deaths);
                survival.dead = false;
            } else if (strcmp(tok, "slot") == 0) {
                char id[64] = {};
                int i = 0;
                long long c = 0;
                if (fscanf(f, " %d %63s %lld", &i, id, &c) == 3 && i >= 0 && i < Inventory::kSlots && c > 0 && itemDB.find(id)) {
                    InvSlot& s = inventory.grid()[i];
                    if (s.empty()) { s.id = id; s.centi = c; } else inventory.addAmount(id, c / 100.0);
                }
            } else if (strcmp(tok, "item") == 0) {
                char id[64] = {};
                int n = 0;
                if (fscanf(f, " %63s %d", id, &n) == 2) inventory.add(id, n);
            }
        }
    };
}
