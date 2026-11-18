#include "gamesave.h"
#include "structures.h"
#include "trapping.h"
#include "combat.h"
#include "events.h"
#include "gamerng.h"
#include "headless.h"
#include "menus.h"
#include "net.h"
#include "netid.h"
#include "inventory.h"
#include "player.h"
#include "playerrecord.h"
#include "world.h"
#include "ship.h"
#include "survival.h"
#include "hotbar.h"
#include "helpmanual.h"
#include "rocket.h"
#include "terraintool.h"
#include "navigation.h"
#include "farming.h"
#include "animals.h"
#include "particles.h"
#include "ambience.h"
#include "suit.h"
#include "toolwear.h"
#include "recipebook.h"
#include "station.h"

#include <raylib.h>
#include <algorithm>
#include <cstring>
#include <ctime>
#include <filesystem>

// Save format 13 (#346): the world file's extras hold the world (ids, the RNG, farms, animals, structures, graves, events,
// rockets, stations, difficulty, teams with their quest books); the player's own state goes to their record,
// <save>.players/<uuid>.player. Format 12 saves had the player's lines in the extras: they still load (the player lines are
// read wherever they are), and the next save writes 13.
static int savedCols = 9;          // the loading save's inventory width (inv_cols; older saves: 9) (#495)
static bool playerLoaded = false;  // the load found the player's lines (a record, or a format 12 save's)
static double recordFirstSeen = 0.0, recordPlayTime = 0.0, recordLoadedAt = 0.0;
static std::function<void()> newPlayerHook;
static std::string legacyOwner;   // a server that loaded a format 12 save: the id its player is saved under as the owner

namespace {
// This game has a player of its own to save (not a dedicated server, which has none; not a joined client, whose world
// and player are the server's).
bool LocalPlayer() { return !Headless() && gNetRole != NetRole::Client; }

void WritePlayerLines(FILE* f) {
    fprintf(f, "player %.6f %.6f %.6f %f %f %f %d %f\n", player.position.x, player.position.y, player.position.z,
            player.forward.x, player.forward.y, player.forward.z, player.flyMode ? 1 : 0, player.jetpackFuel);
    fprintf(f, "player_up %f %f %f\n", player.up.x, player.up.y, player.up.z);
    fprintf(f, "player_vel %f %f %f\n", player.velocity.x, player.velocity.y, player.velocity.z);
    fprintf(f, "inv_cols %d\n", Inventory::kCols); // before the slots: older saves (no line) had 9 columns (#495)
    for (int i = 0; i < Inventory::kCells; ++i) { // with the off hand (#165)
        const InvSlot& s = inventory.grid()[i];
        if (!s.empty()) fprintf(f, "slot %d %s %lld\n", i, s.id.c_str(), (long long)s.centi);
        if (!s.empty() && s.heat > 0.0f) fprintf(f, "slot_heat %d %.3f\n", i, s.heat); // still hot (#580)
    }
    fprintf(f, "creative %d\n", terrainTool.creative ? 1 : 0);
    Hotbar_Save(f);
    HelpManual_Save(f);
    Nav_Save(f);
    Suit_Save(f);
    ToolWear_Save(f);
    fprintf(f, "survival %f %f %f %d\n", survival.oxygen, survival.power, survival.health, survival.deaths);
    fprintf(f, "food %f\n", survival.food);
    if (survival.hasBed) fprintf(f, "bed %s %u %d %d %d\n", survival.bedBody.c_str(), survival.bedPlot, survival.bedCell[0], survival.bedCell[1], survival.bedCell[2]); // #800: the player's spawn bed
}

// One of the player's lines (a record's, or a format 12 save's extras): true if it was one.
bool ReadPlayerLine(const char* tok, FILE* f) {
    if (strcmp(tok, "food") == 0) { float v = 100.0f; if (fscanf(f, " %f", &v) == 1) survival.food = std::clamp(v, 0.0f, 100.0f); return true; }
    if (strcmp(tok, "bed") == 0) { // #800: the spawn bed
        char b[64] = {}; unsigned pl = 0; int x = 0, y = 0, z = 0;
        if (fscanf(f, " %63s %u %d %d %d", b, &pl, &x, &y, &z) == 5) { survival.hasBed = true; survival.bedBody = b; survival.bedPlot = pl; survival.bedCell[0] = x; survival.bedCell[1] = y; survival.bedCell[2] = z; }
        return true;
    }
    if (strcmp(tok, "creative") == 0) { int c = 0; if (fscanf(f, " %d", &c) == 1) terrainTool.creative = c != 0; return true; }
    if (Hotbar_LoadLine(tok, f)) return true;
    if (HelpManual_LoadLine(tok, f)) return true;
    if (Nav_LoadLine(tok, f)) return true;
    if (Suit_LoadLine(tok, f)) return true;
    if (ToolWear_LoadLine(tok, f)) return true;
    if (strcmp(tok, "player") == 0) {
        int fly = 0;
        fscanf(f, " %lf %lf %lf %f %f %f %d %f", &player.position.x, &player.position.y, &player.position.z,
               &player.forward.x, &player.forward.y, &player.forward.z, &fly, &player.jetpackFuel);
        player.flyMode = fly != 0;
        survival.lastJetpackFuel = player.jetpackFuel; // a load is not a recharge
        player.velocity = { 0, 0, 0 };
        player.grounded = false;
        player.up = { 0, 0, 0 }; // set by "player_up"; older saves: rebuilt from the body at the next update
        playerLoaded = true;
        return true;
    }
    if (strcmp(tok, "player_up") == 0) { fscanf(f, " %f %f %f", &player.up.x, &player.up.y, &player.up.z); return true; }
    if (strcmp(tok, "player_vel") == 0) { fscanf(f, " %f %f %f", &player.velocity.x, &player.velocity.y, &player.velocity.z); return true; } // older saves: stays zero
    if (strcmp(tok, "survival") == 0) {
        fscanf(f, " %f %f %f %d", &survival.oxygen, &survival.power, &survival.health, &survival.deaths);
        survival.dead = false;
        return true;
    }
    if (strcmp(tok, "inv_cols") == 0) { if (fscanf(f, " %d", &savedCols) != 1) savedCols = 9; return true; }
    if (strcmp(tok, "slot") == 0) {
        char id[64] = {};
        int i = 0;
        long long c = 0;
        const bool read = fscanf(f, " %d %63s %lld", &i, id, &c) == 3;
        i = Inventory::SavedCell(i, savedCols); // a 9-wide save: its hotbar and off hand move up (#495)
        if (read && i >= 0 && i < Inventory::kCells && c > 0 && itemDB.find(id)) {
            InvSlot& s = inventory.grid()[i];
            if (s.empty()) { s.id = id; s.centi = c; } else { ++inventory.muted; inventory.addAmount(id, c / 100.0); --inventory.muted; } // no pickup feed for a load
        }
        return true;
    }
    if (strcmp(tok, "slot_heat") == 0) { // #580: after its slot line
        int i = 0; float h = 0.0f;
        if (fscanf(f, " %d %f", &i, &h) == 2 && i >= 0 && i < Inventory::kCells && !inventory.grid()[i].empty()) inventory.grid()[i].heat = std::clamp(h, 0.0f, 1.0f);
        return true;
    }
    if (strcmp(tok, "item") == 0) {
        char id[64] = {};
        int n = 0;
        if (fscanf(f, " %63s %d", id, &n) == 2) { ++inventory.muted; inventory.add(id, n); --inventory.muted; }
        return true;
    }
    return false;
}

bool ReadRecord(const std::string& path) {
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return false;
    int version = 0;
    if (fscanf(f, " RAYLIBGAME_PLAYER %d", &version) != 1 || version < 1) { fclose(f); return false; }
    savedCols = 9;
    char tok[64] = {};
    while (fscanf(f, " %63s", tok) == 1 && strcmp(tok, "end_record") != 0) {
        char word[128] = {};
        if (strcmp(tok, "uuid") == 0 || strcmp(tok, "name") == 0) { fscanf(f, " %127s", word); continue; }
        if (strcmp(tok, "team") == 0) { int t = 0; fscanf(f, " %d", &t); continue; } // (the world's team lines are the truth)
        if (strcmp(tok, "first_seen") == 0) { fscanf(f, " %lf", &recordFirstSeen); continue; }
        if (strcmp(tok, "last_seen") == 0) { double t = 0; fscanf(f, " %lf", &t); continue; }
        if (strcmp(tok, "play_time") == 0) { fscanf(f, " %lf", &recordPlayTime); continue; }
        if (!ReadPlayerLine(tok, f)) { char rest[512]; if (!fgets(rest, sizeof(rest), f)) break; } // an unknown line: skipped
    }
    fclose(f);
    return true;
}

bool WriteRecord(const std::string& uuid) {
    std::error_code ec;
    std::filesystem::create_directories(PlayerRecord_Dir(), ec);
    const double now = (double)time(nullptr);
    return Save_WriteAtomic(PlayerRecord_Path(uuid), [&](FILE* f) {
        fprintf(f, "RAYLIBGAME_PLAYER 1\n");
        fprintf(f, "uuid %s\n", uuid.c_str());
        fprintf(f, "name %s\n", settings.mpName.empty() ? "player" : settings.mpName.c_str()); // (shown only; Steam's later)
        fprintf(f, "team %d\n", Teams_Of(uuid));
        fprintf(f, "first_seen %.0f\nlast_seen %.0f\n", recordFirstSeen > 0.0 ? recordFirstSeen : now, now);
        fprintf(f, "play_time %.0f\n", recordPlayTime + std::max(0.0, GetTime() - recordLoadedAt));
        WritePlayerLines(f);
        fprintf(f, "end_record\n");
        return true;
    });
}
} // namespace

bool GameSave_PlayerLoaded() { return playerLoaded; }
void GameSave_SetNewPlayerHook(std::function<void()> hook) { newPlayerHook = std::move(hook); }

void GameSave_Register() {
    saveExtrasHook = [](FILE* f) {
        NetId_Save(f); // first: loads allocate ids for lines from older saves after it (#317)
        GameRng_Save(f); // the gameplay RNG continues where it was (#316)
        Farm_Save(f);
        Animals_Save(f);
        Structures_Save(f); // #164
        Trapping_Save(f);   // #176
        Combat_Save(f);     // graves (#165)
        Events_Save(f);     // #177
        Rocket_Save(f);
        Station_Save(f);
        fprintf(f, "difficulty %d\n", worldDifficulty);
        // Teams and their quest books (#346: quests and tiers are a team's); this game's player is in one (a team of one to
        // start with), and the world's owner is whoever saved it first.
        // A server that loaded a format 12 save writes its player as the owner (a new id), so nothing is lost: single player
        // opening the world again continues as that owner.
        const std::string me = LocalPlayer() ? PlayerId_Local() : legacyOwner;
        if (!me.empty()) {
            Teams_Join(me, settings.mpName.empty() ? "player" : settings.mpName, settings.mpSuit);
            if (Teams_Owner().empty()) Teams_SetOwner(me);
        }
        Teams_Save(f, me, [](FILE* q) { QuestBook_Save(q); });
        // The player's own record, beside the world (the world file never holds it).
        if (!me.empty() && !WriteRecord(me)) TraceLog(LOG_WARNING, "save: could not write the player record %s", PlayerRecord_Path(me).c_str());
    };
    loadExtrasHook = [](FILE* f) {
        inventory.clear();
        terrainTool.creative = false;
        Hotbar_Reset();
        Nav_Reset();
        QuestBook_Reset();
        Farm_Reset();
        Animals_Reset();
        Particles_Clear();
        Ambience_Reset();
        Station_Reset();
        Suit_Reset();
        ToolWear_Reset();
        RecipeBook_Reset();
        Combat_Reset(); // graves: the save holds them all (#165)
        Teams_Reset();
        survival.hasBed = false; // #800: the record's (or a format 12 save's) bed line sets it
        savedCols = 9;  // a save without an inv_cols line: the 9-wide grid (#495)
        worldDifficulty = 1; // a save without the line: normal
        playerLoaded = false;
        recordFirstSeen = recordPlayTime = 0.0;
        recordLoadedAt = GetTime();
        const bool local = LocalPlayer();
        const std::string me = local ? PlayerId_Local() : std::string();
        char tok[64] = {};
        while (fscanf(f, " %63s", tok) == 1 && strcmp(tok, "end_extras") != 0) {
            if (NetId_LoadLine(tok, f)) continue;
            if (GameRng_LoadLine(tok, f)) continue;
            if (Teams_LoadLine(tok, f, me)) continue;
            if (QuestBook_LoadLine(tok, f)) continue; // (format 12: the quest book was the save's own)
            if (Farm_LoadLine(tok, f)) continue;
            if (Animals_LoadLine(tok, f)) continue;
            if (Structures_LoadLine(tok, f)) continue;
            if (Trapping_LoadLine(tok, f)) continue;
            if (Combat_LoadLine(tok, f)) continue;
            if (Events_LoadLine(tok, f)) continue;
            if (Rocket_LoadLine(tok, f)) continue;
            if (Station_LoadLine(tok, f)) continue;
            if (strcmp(tok, "difficulty") == 0) { int d = 1; if (fscanf(f, " %d", &d) == 1) worldDifficulty = std::clamp(d, 0, 2); continue; }
            if (ReadPlayerLine(tok, f)) continue; // format 12: the player's lines were the save's own
            if (strcmp(tok, "ship") == 0) { // old saves: the placeholder ship is gone (its line is read and dropped)
                char rest[512]; if (!fgets(rest, sizeof(rest), f)) break;
                ship.piloted = false; ship.spawned = false; ship.velocity = { 0, 0, 0 };
                continue;
            }
            if (strcmp(tok, "ship_vel") == 0) { char rest[512]; if (!fgets(rest, sizeof(rest), f)) break; continue; }
        }
        legacyOwner.clear();
        if (!local) { if (Headless() && playerLoaded && Teams_Owner().empty()) legacyOwner = PlayerId_New(); return; } // a server: a format 12 save's player, kept

        // The player's record: ours; else (a world made on another PC, opened here) its owner's, which we continue as.
        std::string from = PlayerRecord_Path(me);
        const std::string& own = Teams_Owner();
        if (!FileExists(from.c_str()) && !own.empty() && own != me && FileExists(PlayerRecord_Path(own).c_str())) {
            from = PlayerRecord_Path(own);
            Teams_ReplaceMember(own, me);
        }
        if (!playerLoaded && FileExists(from.c_str())) ReadRecord(from); // (a format 12 save already gave its lines)
        if (!playerLoaded && newPlayerHook) newPlayerHook(); // nobody's player here for us: start at the world's start
    };
}

// --test-save-split (#346): the world and the player's record, saved and loaded back four ways.
// 1 a save: format 13, the world file without the player's lines, the team (with its quest book) in it, the record beside it;
//   loaded back: the inventory, vitals, position and quests as they were;
// 2 a format 12 save (the player's lines in the world file, no record): loads the same; its next save writes 13 and a record;
// 3 a world made on another PC (its owner's record, not ours): we continue as its owner;
// 4 no record at all: a new player at the start with the starting kit, the world's quests kept.
bool GameSave_TestSplit(Scene& world) {
    namespace fs = std::filesystem;
    const std::string save = "saves/split_test.sav";
    const std::string oldPath = WorldSavePath();
    SetWorldSavePath(save.c_str());
    std::error_code ec;
    fs::remove(save, ec); fs::remove_all(PlayerRecord_Dir(), ec);
    auto slurp = [](const std::string& p) { std::string s; if (FILE* f = fopen(p.c_str(), "rb")) { char b[4096]; size_t n; while ((n = fread(b, 1, sizeof(b), f)) > 0) s.append(b, n); fclose(f); } return s; };
    auto spit = [](const std::string& p, const std::string& s) { if (FILE* f = fopen(p.c_str(), "wb")) { fwrite(s.data(), 1, s.size(), f); fclose(f); } };
    auto lines = [](const std::string& s, const char* prefix) { int n = 0; size_t at = 0; const size_t len = strlen(prefix);
        while (at < s.size()) { const size_t e = s.find('\n', at); if (s.compare(at, len, prefix) == 0) ++n; at = e == std::string::npos ? s.size() : e + 1; } return n; };
    const std::string me = PlayerId_Local();
    // The state to keep.
    inventory.clear();
    inventory.add("iron_ingot", 7);
    survival.health = 63.0f; survival.food = 41.0f;
    QuestBook_SetTier(2); QuestBook_Count("test_split", 7.0);
    const DVec3 pos = player.position;
    auto same = [&](const char* when) {
        const bool ok = inventory.count("iron_ingot") == 7 && fabsf(survival.health - 63.0f) < 0.01f && fabsf(survival.food - 41.0f) < 0.01f &&
                        QuestBook_Tier() == 2 && fabs(QuestBook_CountOf("test_split") - 7.0) < 0.01 && (player.position - pos).length() < 0.01;
        if (!ok) printf("bench: save split: %s: iron %d, health %.1f, food %.1f, tier %d, count %.1f, moved %.3f m\n", when, inventory.count("iron_ingot"),
                        survival.health, survival.food, QuestBook_Tier(), QuestBook_CountOf("test_split"), (player.position - pos).length());
        return ok;
    };
    auto scramble = [&] { inventory.clear(); survival.health = 100.0f; survival.food = 100.0f; QuestBook_Reset(); player.position = player.position + DVec3(5.0, 0.0, 0.0); };
    // 1
    bool ok1 = SaveWorld(world);
    const std::string w = slurp(save), rec = slurp(PlayerRecord_Path(me));
    ok1 = ok1 && w.rfind("RAYLIBGAME_SAVE 13", 0) == 0 && lines(w, "slot ") == 0 && lines(w, "player ") == 0 && lines(w, "survival ") == 0 &&
          lines(w, "team ") == 1 && lines(w, "team_quests ") == 1 && lines(w, "owner ") == 1 && rec.rfind("RAYLIBGAME_PLAYER 1", 0) == 0 &&
          lines(rec, "slot ") >= 1 && lines(rec, "survival ") == 1 && lines(rec, "quest") == 0 && !fs::exists(save + ".tmp");
    scramble();
    ok1 = LoadWorldSave(world) && GameSave_PlayerLoaded() && same("saved and loaded") && ok1;
    // 2: the same save as format 12 would have had it: the record's lines in the extras, no teams, no record.
    std::string v12 = w;
    v12.replace(0, strlen("RAYLIBGAME_SAVE 13"), "RAYLIBGAME_SAVE 12");
    { // the team lines out, the quest block's lines kept in place (format 12 had them in the extras)
        std::string out; size_t at = 0;
        while (at < v12.size()) {
            const size_t e = v12.find('\n', at), end = e == std::string::npos ? v12.size() : e + 1;
            const std::string l = v12.substr(at, end - at);
            if (l.rfind("team", 0) != 0 && l.rfind("owner ", 0) != 0) out += l;
            at = end;
        }
        v12 = out;
    }
    std::string body;
    { size_t a = rec.find("first_seen"); a = rec.find('\n', rec.find("play_time")) + 1; body = rec.substr(a, rec.find("end_record") - a); }
    { const size_t x = v12.find("\nextras"); v12.insert(v12.find('\n', x + 1) + 1, body); } // after the "extras" line (CRLF or LF)
    spit(save, v12);
    fs::remove_all(PlayerRecord_Dir(), ec);
    scramble();
    bool ok2 = LoadWorldSave(world) && GameSave_PlayerLoaded() && same("format 12 loaded");
    ok2 = ok2 && SaveWorld(world) && fs::exists(PlayerRecord_Path(me)) && lines(slurp(save), "slot ") == 0 && slurp(save).rfind("RAYLIBGAME_SAVE 13", 0) == 0;
    // 3: the owner is someone else (their record, their team): we continue as them.
    const std::string them = "0f0f0f0f-1234-4abc-8def-0123456789ab";
    std::string w3 = slurp(save);
    for (size_t at; (at = w3.find(me)) != std::string::npos;) w3.replace(at, me.size(), them);
    spit(save, w3);
    fs::rename(PlayerRecord_Path(me), PlayerRecord_Path(them), ec);
    scramble();
    bool ok3 = LoadWorldSave(world) && GameSave_PlayerLoaded() && same("another PC's world") && Teams_Of(me) != 0 && Teams_Owner() == me;
    // 4: no record anywhere: a new player with the kit at the start; the world's quests stay.
    ok3 = ok3 && SaveWorld(world);
    fs::remove_all(PlayerRecord_Dir(), ec);
    scramble();
    const bool loaded4 = LoadWorldSave(world);
    const bool ok4 = loaded4 && !GameSave_PlayerLoaded() && inventory.count("telescope") == 1 && inventory.count("iron_ingot") == 0 && QuestBook_Tier() == 2;
    // 5: a dedicated server (no player of its own) loads the format 12 save and saves it: its player becomes the owner's
    // record (with the quests in the owner's team), and single player opening it again continues as that owner.
    spit(save, v12);
    fs::remove_all(PlayerRecord_Dir(), ec);
    scramble();
    gHeadless = true;
    bool ok5 = LoadWorldSave(world) && SaveWorld(world);
    gHeadless = false;
    const std::string w5 = slurp(save);
    ok5 = ok5 && !Teams_Owner().empty() && Teams_Owner() != me && fs::exists(PlayerRecord_Path(Teams_Owner())) && lines(w5, "team ") == 1 && lines(w5, "slot ") == 0;
    scramble();
    ok5 = ok5 && LoadWorldSave(world) && GameSave_PlayerLoaded() && same("a server's resave, opened alone") && Teams_Owner() == me;
    const bool all = ok1 && ok2 && ok3 && ok4 && ok5;
    printf("bench: save split (#346): saved as 13 with the record beside it and loaded back %d; a format 12 save loaded and resaved as 13 %d; "
           "another PC's world continued as its owner %d; no record: a new player with the kit, the quests kept %d; "
           "a server's resave of a format 12 save kept its player as the owner %d -> %s\n",
           ok1, ok2, ok3, ok4, ok5, all ? "PASS" : "FAIL");
    fs::remove(save, ec); fs::remove_all(PlayerRecord_Dir(), ec); fs::remove_all(save + ".regions", ec);
    SetWorldSavePath(oldPath.c_str());
    return all;
}
