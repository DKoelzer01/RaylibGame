#include "gamesave.h"
#include "blueprint.h"
#include "fishing.h"
#include "structures.h"
#include "trapping.h"
#include "combat.h"
#include "claims.h"
#include "events.h"
#include "gamerng.h"
#include "headless.h"
#include "menus.h"
#include "net.h"
#include "netid.h"
#include "inventory.h"
#include "player.h"
#include "playerrecord.h"
#include "teamquests.h"
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
#include "flares.h"
#include "particles.h"
#include "ambience.h"
#include "suit.h"
#include "toolwear.h"
#include "playersession.h"
#include "skills.h"
#include "magic.h"      // (#1249) Magic_Enabled
#include "magicsight.h"
#include "potions.h" // (#1254)
#include "recipebook.h"
#include "station.h"
#include "saveversion.h"
#include "savemigrate.h"

#include <raylib.h>
#include <algorithm>
#include <cstring>
#include <ctime>
#include <filesystem>

// Save format 13 (#346): the world file's extras hold the world (ids, the RNG, farms, animals, structures, graves, events,
// rockets, stations, difficulty, teams with their quest books); the player's own state goes to their record,
// <save>.players/<uuid>.player. Saves from before 13 are refused (#1693); a record has its own number and migration chain.
static int savedCols = 9;          // the loading save's inventory width (inv_cols; a record without the line: 9) (#495)
static bool playerLoaded = false;  // the load found the player's record
static double recordFirstSeen = 0.0, recordPlayTime = 0.0, recordLoadedAt = 0.0;
static std::function<void()> newPlayerHook;
static bool bodyFramed = false;  // the record being read gave a player_body line (its up and velocity are rebuilt)
// Body-frame positions read during a world load (#901): the records are read in the save's extras, before the load puts the
// bodies back in their saved rotation, so the conversion to world space waits for loadDoneHook. Outside a load (a record
// read for a joining player) the bodies are as they are and the conversion at once is right.
struct PendingFrame { Player* p; unsigned long long body; DVec3 local; };
static std::vector<PendingFrame> pendingFrames;
static bool inWorldLoad = false;
static bool loadedFrameSet = false; // the local player's frame from the last load, for the loading screen's hold (once)
static unsigned long long loadedFrameBody = 0;
static DVec3 loadedFrameLocal{};
bool GameSave_TakeLoadedFrame(unsigned long long& bodyNetId, DVec3& local) {
    if (!loadedFrameSet) return false;
    loadedFrameSet = false; bodyNetId = loadedFrameBody; local = loadedFrameLocal;
    return true;
}

namespace {
// This game has a player of its own to save (not a dedicated server, which has none; not a joined client, whose world
// and player are the server's).
bool LocalPlayer() { return !Headless() && gNetRole != NetRole::Client; }

// Whose lines: a player, their inventory and survival; `local` also the game's own (hotbar, help, map pins, suit, tool
// wear, creative), which a server does not hold for a remote player: those lines are kept as text and written back as
// they were (a record made in single player keeps its hotbar when the world is hosted).
struct RecordTarget { Player& p; Inventory& inv; SurvivalState& sv; bool local; std::string& kept; SkillState* skills = nullptr; }; // skills (#920): a server's remote player's, read and written as lines

void WritePlayerLines(FILE* f, RecordTarget t) {
    fprintf(f, "player %.6f %.6f %.6f %f %f %f %d %f\n", t.p.position.x, t.p.position.y, t.p.position.z,
            t.p.forward.x, t.p.forward.y, t.p.forward.z, t.p.flyMode ? 1 : 0, t.p.jetpackFuel);
    // (#346) also in the body's own frame: a server's world time runs on while the player is away (the planet turns), so
    // the absolute position would be wrong when they come back.
    if (const Planetoid* b = t.p.nearest; b && NetId_Kind(b->netId) == NetKind::Body) {
        const DVec3 l = b->toLocalD(t.p.position);
        fprintf(f, "player_body %llu %.6f %.6f %.6f\n", (unsigned long long)b->netId, l.x, l.y, l.z);
    }
    fprintf(f, "player_up %f %f %f\n", t.p.up.x, t.p.up.y, t.p.up.z);
    fprintf(f, "player_vel %f %f %f\n", t.p.velocity.x, t.p.velocity.y, t.p.velocity.z);
    fprintf(f, "inv_cols %d\n", Inventory::kCols); // before the slots: older saves (no line) had 9 columns (#495)
    for (int i = 0; i < Inventory::kCells; ++i) { // with the off hand (#165)
        const InvSlot& s = t.inv.grid()[i];
        if (!s.empty()) fprintf(f, "slot %d %s %lld\n", i, s.id.c_str(), (long long)s.centi);
        if (!s.empty() && s.heat > 0.0f) fprintf(f, "slot_heat %d %.3f\n", i, s.heat); // still hot (#580)
    }
    if (t.local) {
        fprintf(f, "creative %d\n", terrainTool.creative ? 1 : 0);
        Hotbar_Save(f);
        HelpManual_Save(f);
        Nav_Save(f);
        Suit_Save(f);
        ToolWear_Save(f);
        QuestBook_SavePins(f); // (#825) the player's pinned quests
        fputs(Skills_SaveText(LocalSession().skills).c_str(), f); // (#883) XP and perks
        if (Magic_Enabled() || LocalSession().magic.sight) fputs(MagicSight_SaveText(LocalSession().magic).c_str(), f); // (#1249)
        if (Magic_Enabled()) fputs(Potion_SaveText(Potion_Local()).c_str(), f); // (#1254) potions running, the cleansing draught's day
    } else {
        fputs(t.kept.c_str(), f);
        if (t.skills) fputs(Skills_SaveText(*t.skills).c_str(), f); // (#920) a server's remote player's XP and perks
    }
    fprintf(f, "survival %f %f %f %d\n", t.sv.oxygen, t.sv.power, t.sv.health, t.sv.deaths);
    fprintf(f, "food %f\n", t.sv.food);
    if (t.sv.hasBed) fprintf(f, "bed %s %u %d %d %d\n", t.sv.bedBody.c_str(), t.sv.bedPlot, t.sv.bedCell[0], t.sv.bedCell[1], t.sv.bedCell[2]); // #800: the player's spawn bed
}

// One of the player's lines (a record's): true if it was one.
bool ReadPlayerLine(const char* tok, FILE* f, RecordTarget t) {
    if (strcmp(tok, "food") == 0) { float v = 100.0f; if (fscanf(f, " %f", &v) == 1) t.sv.food = std::clamp(v, 0.0f, 100.0f); return true; }
    if (strcmp(tok, "bed") == 0) { // #800: the spawn bed
        char b[64] = {}; unsigned pl = 0; int x = 0, y = 0, z = 0;
        if (fscanf(f, " %63s %u %d %d %d", b, &pl, &x, &y, &z) == 5) { t.sv.hasBed = true; t.sv.bedBody = b; t.sv.bedPlot = pl; t.sv.bedCell[0] = x; t.sv.bedCell[1] = y; t.sv.bedCell[2] = z; }
        return true;
    }
    if (!t.local && t.skills && (strcmp(tok, "skill_xp") == 0 || strcmp(tok, "perk") == 0)) { // (#920) a remote player's skills, not kept as text
        char rest[256] = {};
        if (fgets(rest, sizeof(rest), f)) Skills_LoadLine(*t.skills, (std::string(tok) + rest).c_str());
        return true;
    }
    if (t.local) {
        if (strcmp(tok, "creative") == 0) { int c = 0; if (fscanf(f, " %d", &c) == 1) terrainTool.creative = c != 0; return true; }
        if (Hotbar_LoadLine(tok, f)) return true;
        if (HelpManual_LoadLine(tok, f)) return true;
        if (Nav_LoadLine(tok, f)) return true;
        if (Suit_LoadLine(tok, f)) return true;
        if (ToolWear_LoadLine(tok, f)) return true;
        if (QuestBook_LoadPinLine(tok, f)) return true; // (#825)
        if (strcmp(tok, "skill_xp") == 0 || strcmp(tok, "perk") == 0) { // (#883) the rest of the line is the skill's or perk's
            char rest[256] = {};
            if (fgets(rest, sizeof(rest), f)) Skills_LoadLine(LocalSession().skills, (std::string(tok) + rest).c_str());
            return true;
        }
        if (strncmp(tok, "magic_", 6) == 0) { // (#1249) magic sight, the pool, studies
            char rest[256] = {};
            if (fgets(rest, sizeof(rest), f)) { const std::string l = std::string(tok) + rest; if (!Potion_LoadLine(Potion_Local(), l.c_str())) MagicSight_LoadLine(LocalSession().magic, l.c_str()); } // (#1254) potion lines too
            return true;
        }
    }
    if (strcmp(tok, "player") == 0) {
        int fly = 0;
        fscanf(f, " %lf %lf %lf %f %f %f %d %f", &t.p.position.x, &t.p.position.y, &t.p.position.z,
               &t.p.forward.x, &t.p.forward.y, &t.p.forward.z, &fly, &t.p.jetpackFuel);
        t.p.flyMode = fly != 0;
        t.sv.lastJetpackFuel = t.p.jetpackFuel; // a load is not a recharge
        t.p.velocity = { 0, 0, 0 };
        t.p.grounded = false;
        t.p.up = { 0, 0, 0 }; // set by "player_up"; older saves: rebuilt from the body at the next update
        if (t.local) playerLoaded = true;
        return true;
    }
    if (strcmp(tok, "player_body") == 0) { // (#346) the body's frame wins over the absolute position
        unsigned long long id = 0; DVec3 l{};
        if (fscanf(f, " %llu %lf %lf %lf", &id, &l.x, &l.y, &l.z) == 4)
            if (const Planetoid* b = Planetoid::findNet((NetId)id)) {
                t.p.position = b->toWorldD(l); t.p.up = { 0, 0, 0 }; bodyFramed = true;
                if (inWorldLoad) pendingFrames.push_back({ &t.p, id, l }); // again once the body has its saved rotation
            }
        return true;
    }
    if (strcmp(tok, "player_up") == 0) { Vector3 u{}; fscanf(f, " %f %f %f", &u.x, &u.y, &u.z); if (!bodyFramed) t.p.up = u; return true; } // (with player_body: rebuilt from the body)
    if (strcmp(tok, "player_vel") == 0) { fscanf(f, " %f %f %f", &t.p.velocity.x, &t.p.velocity.y, &t.p.velocity.z); if (bodyFramed) t.p.velocity = { 0, 0, 0 }; return true; } // older saves: stays zero
    if (strcmp(tok, "survival") == 0) {
        fscanf(f, " %f %f %f %d", &t.sv.oxygen, &t.sv.power, &t.sv.health, &t.sv.deaths);
        t.sv.dead = false;
        return true;
    }
    if (strcmp(tok, "inv_cols") == 0) { if (fscanf(f, " %d", &savedCols) != 1) savedCols = 9; return true; }
    if (strcmp(tok, "slot") == 0) {
        char id[64] = {};
        int i = 0;
        long long c = 0;
        const bool read = fscanf(f, " %d %63s %lld", &i, id, &c) == 3;
        if (const char* now = RetiredItem(id); *now) snprintf(id, sizeof id, "%s", now); // #1090
        i = Inventory::SavedCell(i, savedCols); // a 9-wide save: its hotbar and off hand move up (#495)
        if (read && i >= 0 && i < Inventory::kCells && c > 0 && itemDB.find(id)) {
            InvSlot& s = t.inv.grid()[i];
            if (s.empty()) { s.id = id; s.centi = c; } else { ++t.inv.muted; t.inv.addAmount(id, c / 100.0); --t.inv.muted; } // no pickup feed for a load
        }
        return true;
    }
    if (strcmp(tok, "slot_heat") == 0) { // #580: after its slot line
        int i = 0; float h = 0.0f;
        if (fscanf(f, " %d %f", &i, &h) == 2 && i >= 0 && i < Inventory::kCells && !t.inv.grid()[i].empty()) t.inv.grid()[i].heat = std::clamp(h, 0.0f, 1.0f);
        return true;
    }
    if (strcmp(tok, "item") == 0) {
        char id[64] = {};
        int n = 0;
        if (fscanf(f, " %63s %d", id, &n) == 2) { if (const char* now = RetiredItem(id); *now) snprintf(id, sizeof id, "%s", now); ++t.inv.muted; t.inv.add(id, n); --t.inv.muted; }
        return true;
    }
    return false;
}

bool ReadRecordInto(const std::string& path, RecordTarget t, PlayerRecordMeta& m) {
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return false;
    int version = 0;
    if (fscanf(f, " RAYLIBGAME_PLAYER %d", &version) != 1 || version < 1) { fclose(f); return false; }
    if (version > kPlayerRecordVersion) { // (#1693) a newer build's record: not read (the next save keeps it as a .bak)
        fclose(f);
        TraceLog(LOG_ERROR, "save: %s is a player record from a newer build (version %d; this build reads up to %d)", path.c_str(), version, kPlayerRecordVersion);
        return false;
    }
    std::string migrated;
    if (version < kPlayerRecordVersion) { // (#1693) through the record's chain, in memory, into a temporary file read below
        fclose(f);
        std::vector<std::string> lines;
        MigrateContext ctx;
        std::string error;
        migrated = path + ".migrating";
        if (!Save_ReadLines(path, lines) || !PlayerRecord_Migrate(version, lines, ctx, error) || !Save_WriteLines(migrated, lines) || !(f = fopen(migrated.c_str(), "r"))) {
            TraceLog(LOG_ERROR, "save: could not upgrade the player record %s from version %d: %s", path.c_str(), version, error.c_str());
            std::error_code ec;
            std::filesystem::remove(migrated, ec);
            return false;
        }
        int upgraded = 0;
        if (fscanf(f, " RAYLIBGAME_PLAYER %d", &upgraded) != 1) { fclose(f); std::error_code ec; std::filesystem::remove(migrated, ec); return false; }
    }
    savedCols = 9;
    bodyFramed = false;
    char tok[64] = {};
    while (fscanf(f, " %63s", tok) == 1 && strcmp(tok, "end_record") != 0) {
        char word[128] = {};
        if (strcmp(tok, "uuid") == 0) { if (fscanf(f, " %127s", word) == 1) m.uuid = word; continue; }
        if (strcmp(tok, "name") == 0) { if (fscanf(f, " %127s", word) == 1) m.name = word; continue; }
        if (strcmp(tok, "team") == 0) { int tm = 0; fscanf(f, " %d", &tm); continue; } // (the world's team lines are the truth)
        if (strcmp(tok, "first_seen") == 0) { fscanf(f, " %lf", &m.firstSeen); continue; }
        if (strcmp(tok, "last_seen") == 0) { double x = 0; fscanf(f, " %lf", &x); continue; }
        if (strcmp(tok, "play_time") == 0) { fscanf(f, " %lf", &m.playTime); continue; }
        if (ReadPlayerLine(tok, f, t)) continue;
        char rest[512];
        if (!fgets(rest, sizeof(rest), f)) break;
        if (!t.local) t.kept += std::string(tok) + rest; // another game's line (hotbar, pins, ...): kept as it was
    }
    fclose(f);
    if (!migrated.empty()) { std::error_code ec; std::filesystem::remove(migrated, ec); }
    return true;
}

bool WriteRecordFrom(const std::string& path, RecordTarget t, const PlayerRecordMeta& m, double sessionSeconds) {
    std::error_code ec;
    std::filesystem::create_directories(std::filesystem::path(path).parent_path(), ec);
    const double now = (double)time(nullptr);
    return Save_WriteAtomic(path, [&](FILE* f) {
        fprintf(f, "RAYLIBGAME_PLAYER %d\n", kPlayerRecordVersion);
        fprintf(f, "uuid %s\n", m.uuid.c_str());
        fprintf(f, "name %s\n", m.name.empty() ? "player" : m.name.c_str()); // (shown only; Steam's later)
        fprintf(f, "team %d\n", Teams_Of(m.uuid));
        fprintf(f, "first_seen %.0f\nlast_seen %.0f\n", m.firstSeen > 0.0 ? m.firstSeen : now, now);
        fprintf(f, "play_time %.0f\n", m.playTime + std::max(0.0, sessionSeconds));
        WritePlayerLines(f, t);
        fprintf(f, "end_record\n");
        return true;
    });
}

std::string noKept; // the local player keeps nothing as text
RecordTarget LocalTarget() { return RecordTarget{ player, inventory, survival, true, noKept }; }

bool ReadRecord(const std::string& path) {
    PlayerRecordMeta m;
    if (!ReadRecordInto(path, LocalTarget(), m)) return false;
    recordFirstSeen = m.firstSeen; recordPlayTime = m.playTime;
    return true;
}

// (#1110) The record file this session's player was read from ("" none: a new player, or the world file's own lines).
std::string recordReadFrom;
// A record's play time, from its file (-1: no file or no line).
double RecordPlayTime(const std::string& path) {
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return -1.0;
    double t = -1.0; char tok[64] = {};
    while (fscanf(f, " %63s", tok) == 1) { if (strcmp(tok, "play_time") == 0) { if (fscanf(f, " %lf", &t) != 1) t = -1.0; break; } char rest[512]; if (!fgets(rest, sizeof rest, f)) break; }
    fclose(f);
    return t;
}
// (#1110, data loss) Never replace a player's record silently: a session that didn't read that file (a fresh player while a
// record existed, a legacy save's own lines), or would write less play time than the file holds, first keeps it beside
// itself as <id>.player.<stamp>.bak. Returns the backup's path ("" none made).
std::string KeepRecordBeforeWrite(const std::string& path, double sessionPlayTime) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (!fs::exists(path, ec)) return std::string();
    const double onDisk = RecordPlayTime(path);
    const bool ours = recordReadFrom == path;
    if (ours && onDisk <= sessionPlayTime + 1.0) return std::string(); // (the record this session is continuing)
    char stamp[32];
    const std::time_t now = std::time(nullptr);
    std::strftime(stamp, sizeof stamp, "%Y%m%d-%H%M%S", std::localtime(&now));
    const std::string bak = path + "." + stamp + ".bak";
    fs::copy_file(path, bak, fs::copy_options::overwrite_existing, ec);
    if (ec) { TraceLog(LOG_WARNING, "save: could not back up the player record %s (%s)", path.c_str(), ec.message().c_str()); return std::string(); }
    printf("save: kept the earlier player record as %s (%s; its play time %.0f s, this session's %.0f s)\n", bak.c_str(),
           ours ? "it holds more play time" : "this session didn't load it", onDisk, sessionPlayTime);
    return bak;
}

bool WriteRecord(const std::string& uuid) {
    PlayerRecordMeta m;
    m.uuid = uuid; m.name = settings.mpName.empty() ? "player" : settings.mpName;
    m.firstSeen = recordFirstSeen; m.playTime = recordPlayTime;
    const std::string path = PlayerRecord_Path(uuid);
    const double session = GetTime() - recordLoadedAt;
    KeepRecordBeforeWrite(path, recordPlayTime + session);
    const bool ok = WriteRecordFrom(path, LocalTarget(), m, session);
    if (ok) { recordReadFrom = path; recordPlayTime += session; recordLoadedAt = GetTime(); } // (this file is ours now; its time counted)
    return ok;
}
} // namespace

// A remote player's record on a server (#346 step 2): their player, inventory and survival; the rest kept as text.
bool PlayerRecord_Read(const std::string& path, Player& p, Inventory& inv, SurvivalState& sv, PlayerRecordMeta& meta, SkillState* skills) {
    meta.kept.clear();
    return ReadRecordInto(path, RecordTarget{ p, inv, sv, false, meta.kept, skills }, meta);
}
bool PlayerRecord_Write(const std::string& path, Player& p, Inventory& inv, SurvivalState& sv, PlayerRecordMeta& meta, double sessionSeconds, SkillState* skills) {
    return WriteRecordFrom(path, RecordTarget{ p, inv, sv, false, meta.kept, skills }, meta, sessionSeconds);
}
std::vector<PlayerRecordMeta> PlayerRecord_List() {
    std::vector<PlayerRecordMeta> out;
    std::error_code ec;
    for (const auto& e : std::filesystem::directory_iterator(PlayerRecord_Dir(), ec)) {
        if (e.path().extension() != ".player") continue;
        FILE* f = fopen(e.path().string().c_str(), "r");
        if (!f) continue;
        PlayerRecordMeta m;
        char tok[64] = {}, word[128] = {};
        int version = 0;
        if (fscanf(f, " RAYLIBGAME_PLAYER %d", &version) == 1)
            while ((m.uuid.empty() || m.name.empty()) && fscanf(f, " %63s", tok) == 1 && strcmp(tok, "end_record") != 0) {
                if (strcmp(tok, "uuid") == 0 && fscanf(f, " %127s", word) == 1) m.uuid = word;
                else if (strcmp(tok, "name") == 0 && fscanf(f, " %127s", word) == 1) m.name = word;
                else { char rest[512]; if (!fgets(rest, sizeof(rest), f)) break; }
            }
        fclose(f);
        if (PlayerId_Valid(m.uuid)) out.push_back(m);
    }
    return out;
}

bool GameSave_PlayerLoaded() { return playerLoaded; }
double GameSave_PlayTime() { return recordPlayTime + std::max(0.0, GetTime() - recordLoadedAt); }
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
        Blueprint_SaveGhosts(f); // blueprint parts not built yet (#1018)
        Claims_Save(f);     // land claims (#343)
        Events_Save(f);     // #177
        Rocket_Save(f);
        Station_Save(f);
        fprintf(f, "difficulty %d\n", worldDifficulty);
        fprintf(f, "xp_rate %g\n", worldXpRate); // (#883)
        // Teams and their quest books (#346: quests and tiers are a team's); this game's player is in one (a team of one to
        // start with), and the world's owner is whoever saved it first.
        // A server has no player of its own: no team join, no record here.
        const std::string me = LocalPlayer() ? PlayerId_Local() : std::string();
        if (!me.empty()) {
            Teams_Join(me, settings.mpName.empty() ? "player" : settings.mpName, settings.mpSuit);
            if (Teams_Owner().empty()) Teams_SetOwner(me);
        }
        Teams_Save(f, me, [](FILE* q) { QuestBook_Save(q); });
        // The player's own record, beside the world (the world file never holds it).
        if (!me.empty() && !WriteRecord(me)) TraceLog(LOG_WARNING, "save: could not write the player record %s", PlayerRecord_Path(me).c_str());
    };
    loadDoneHook = [] { // (#901) body-frame positions back to world space with the bodies in their saved state
        for (const PendingFrame& pf : pendingFrames)
            if (Planetoid* b = Planetoid::findNet((NetId)pf.body)) {
                pf.p->position = b->toWorldD(pf.local);
                pf.p->nearest = b; // (#1881) a save before the next update wrote no player_body line without it

                if (pf.p == &player) { loadedFrameSet = true; loadedFrameBody = pf.body; loadedFrameLocal = pf.local; } // the hold's spot
            }
        pendingFrames.clear();
        Farm_LoadDone(); // (#1085) an old save's farmland parts become tilled ground, with their crops
        inWorldLoad = false;
    };
    loadExtrasHook = [](FILE* f) {
        pendingFrames.clear(); inWorldLoad = true; // (#901) a world load: body frames converted again at its end
        inventory.clear();
        terrainTool.creative = false;
        Hotbar_Reset();
        Nav_Reset();
        QuestBook_Reset();
        Farm_Reset();
        Animals_Reset();
        Flares_Clear(); // (#2062) flares are not saved: none survive into the loaded world
        Fishing_Reset(); // (#902) a line cast in the old world held its body
        Particles_Clear();
        Ambience_Reset();
        Station_Reset();
        Suit_Reset();
        ToolWear_Reset();
        LocalSession().skills = SkillState{}; // (#883)
        LocalSession().magic = MagicPlayerState{}; Potion_Local() = PotionState{}; // (#1249)
        RecipeBook_Reset();
        Combat_Reset(); // graves: the save holds them all (#165)
        Blueprint_ResetGhosts(); // (#1018) and the ghosts
        Claims_Reset();
        Teams_Reset();
        QuestBook_ClearPins(); // (#825) the record's pins come back with it
        TeamQuests_Reset(); // (#338) the server's team books load from the teams as they are used
        survival.hasBed = false; // #800: the record's bed line sets it
        savedCols = 9;  // a save without an inv_cols line: the 9-wide grid (#495)
        worldDifficulty = 1; // a save without the line: normal
        worldXpRate = 1.0f; // (#883) likewise: 1x
        playerLoaded = false;
        bodyFramed = false;
        recordFirstSeen = recordPlayTime = 0.0;
        recordReadFrom.clear(); // (#1110) set below when a record is read
        recordLoadedAt = GetTime();
        const bool local = LocalPlayer();
        const std::string me = local ? PlayerId_Local() : std::string();
        char tok[64] = {};
        while (fscanf(f, " %63s", tok) == 1 && strcmp(tok, "end_extras") != 0) {
            if (NetId_LoadLine(tok, f)) continue;
            if (GameRng_LoadLine(tok, f)) continue;
            if (Teams_LoadLine(tok, f, me)) continue;
            if (QuestBook_LoadLine(tok, f)) continue;
            if (Farm_LoadLine(tok, f)) continue;
            if (Animals_LoadLine(tok, f)) continue;
            if (Structures_LoadLine(tok, f)) continue;
            if (Trapping_LoadLine(tok, f)) continue;
            if (Combat_LoadLine(tok, f)) continue;
            if (Blueprint_LoadGhostLine(tok, f)) continue; // (#1018)
            if (Claims_LoadLine(tok, f)) continue;
            if (Events_LoadLine(tok, f)) continue;
            if (Rocket_LoadLine(tok, f)) continue;
            if (Station_LoadLine(tok, f)) continue;
            if (strcmp(tok, "difficulty") == 0) { int d = 1; if (fscanf(f, " %d", &d) == 1) worldDifficulty = std::clamp(d, 0, 2); continue; }
            if (strcmp(tok, "xp_rate") == 0) { float x = 1.0f; if (fscanf(f, " %f", &x) == 1) worldXpRate = std::clamp(x, 0.5f, 5.0f); continue; } // (#883)
        }
        if (!local) return; // a server: it has no player of its own (the records are its clients')

        // The player's record: ours; else (a world made on another PC, opened here) its owner's, which we continue as.
        std::string from = PlayerRecord_Path(me);
        const std::string& own = Teams_Owner();
        if (!FileExists(from.c_str()) && !own.empty() && own != me && FileExists(PlayerRecord_Path(own).c_str())) {
            from = PlayerRecord_Path(own);
            Teams_ReplaceMember(own, me);
        }
        if (FileExists(from.c_str())) { if (ReadRecord(from)) recordReadFrom = from; }
        if (!playerLoaded && newPlayerHook) newPlayerHook(); // nobody's player here for us: start at the world's start
    };
}

// --test-save-split (#346): the world and the player's record, saved and loaded back four ways.
// 1 a save: format 13, the world file without the player's lines, the team (with its quest book) in it, the record beside it;
//   loaded back: the inventory, vitals, position and quests as they were;
// (#1693: the format 12 cases are gone; saves from before 13 are refused, see case 5.)
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
    // 5 (#1693): a save from a newer build, and one from before 13, are refused with a message and nothing is loaded; then the
    // first save and record put back load as they were.
    auto refused = [&](int ver, const char* words) {
        std::string t = w;
        t.replace(0, strlen("RAYLIBGAME_SAVE 13"), "RAYLIBGAME_SAVE " + std::to_string(ver));
        spit(save, t);
        scramble();
        const bool no = !LoadWorldSave(world) && LastLoadError().find(words) != std::string::npos && inventory.count("iron_ingot") == 0 && slurp(save) == t;
        if (!no) printf("bench: save split: version %d: not refused as it should be (%s)\n", ver, LastLoadError().c_str());
        return no;
    };
    bool ok5 = refused(kSaveVersion + 1, "newer build") && refused(kOldestSaveVersion - 1, "older build");
    spit(save, w);
    fs::create_directories(PlayerRecord_Dir(), ec);
    spit(PlayerRecord_Path(me), rec);
    scramble();
    ok5 = ok5 && LoadWorldSave(world) && GameSave_PlayerLoaded() && same("put back after the refusals") && LastLoadError().empty();
    // 6 (step 3): bodies by NetId: the first body's place in the list made wrong in the file, its state still comes back.
    Planetoid* body0 = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
    bool ok6 = body0 && SaveWorld(world);
    std::string w6 = slurp(save);
    const size_t p0 = w6.find("planetoid 0 netid ");
    ok6 = ok6 && p0 != std::string::npos;
    Vector3 rot0{};
    if (ok6) { w6.replace(p0, strlen("planetoid 0"), "planetoid 7"); spit(save, w6); rot0 = body0->rotation; body0->rotation = { rot0.x + 33.0f, rot0.y, rot0.z }; }
    ok6 = ok6 && LoadWorldSave(world);
    body0 = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
    ok6 = ok6 && body0 && fabsf(body0->rotation.x - rot0.x) < 0.01f && same("bodies by NetId");
    // 7 (step 3): backups: seven saves keep the newest five, each with the world file and the player records.
    fs::remove_all(BackupDir(), ec);
    bool ok7 = true;
    for (int k = 0; k < 7 && ok7; ++k) ok7 = SaveWorld(world);
    std::vector<fs::path> backups;
    for (const auto& e : fs::directory_iterator(BackupDir(), ec)) if (e.is_directory()) backups.push_back(e.path());
    std::sort(backups.begin(), backups.end());
    const fs::path saveName = fs::path(save).filename();
    ok7 = ok7 && (int)backups.size() == std::min(5, BackupsKept()) && BackupsKept() == 5 && fs::exists(backups.back() / saveName) &&
          fs::exists(backups.back() / (saveName.string() + ".players") / (me + ".player"));
    // 8 (#1110, data loss): a record exists, but this session's player didn't come from it (here: the record was away when the
    // world loaded, so a new player started): saving keeps the old record as a .bak before writing; nothing is silently lost.
    const std::string rec1110 = PlayerRecord_Path(me), aside = rec1110 + ".aside";
    auto baks = [&]() { int n = 0; for (const auto& e : fs::directory_iterator(PlayerRecord_Dir(), ec)) if (e.path().string().find(".bak") != std::string::npos) ++n; return n; };
    bool ok8 = SaveWorld(world) && fs::exists(rec1110);
    const std::string before8 = slurp(rec1110);
    fs::rename(rec1110, aside, ec);
    ok8 = ok8 && !ec && LoadWorldSave(world) && !GameSave_PlayerLoaded(); // (a new player: no record seen)
    fs::rename(aside, rec1110, ec);
    const int baks0 = baks();
    ok8 = ok8 && !ec && SaveWorld(world) && baks() == baks0 + 1;
    for (const auto& e : fs::directory_iterator(PlayerRecord_Dir(), ec))
        if (e.path().string().find(".bak") != std::string::npos) ok8 = ok8 && slurp(e.path().string()) == before8; // (the old record, whole)
    // 9 (#1110): the round trip. A record as a server writes it (PlayerRecord_Write: a menu host's player) opens in single
    // player with its items; single player's save keeps them (no .bak: it continues that record); a server reads them back.
    for (const auto& e : fs::directory_iterator(PlayerRecord_Dir(), ec)) if (e.path().string().find(".bak") != std::string::npos) fs::remove(e.path(), ec);
    bool ok9 = false;
    {
        Player hostP = player; Inventory hostInv; SurvivalState hostSv = survival; PlayerRecordMeta meta;
        hostInv.add("copper_ingot", 5); hostInv.add("telescope", 1);
        meta.uuid = me; meta.name = "host"; meta.playTime = 5000.0;
        ok9 = PlayerRecord_Write(rec1110, hostP, hostInv, hostSv, meta, 0.0, nullptr);
        ok9 = ok9 && LoadWorldSave(world) && GameSave_PlayerLoaded() && inventory.count("copper_ingot") == 5; // MP -> SP
        const int baks9 = baks();
        ok9 = ok9 && SaveWorld(world) && baks() == baks9;                                                          // SP save: no backup needed
        Player p2; Inventory inv2; SurvivalState sv2; PlayerRecordMeta m2;
        ok9 = ok9 && PlayerRecord_Read(rec1110, p2, inv2, sv2, m2, nullptr) && inv2.count("copper_ingot") == 5 && m2.playTime >= 5000.0; // SP -> MP
    }
    const bool all = ok1 && ok3 && ok4 && ok5 && ok6 && ok7 && ok8 && ok9;
    printf("bench: save split (#1110): a record this session didn't load is kept as a .bak before a save %d; a server's record opens in single player "
           "and back, play time kept %d\n", ok8, ok9);
    printf("bench: save split (#346): saved as 13 with the record beside it and loaded back %d; a newer or pre-13 save refused, the first put back %d; "
           "another PC's world continued as its owner %d; no record: a new player with the kit, the quests kept %d; "
           "a body found by its NetId %d; seven saves left %d backups (5), "
           "the newest with the world and the record %d -> %s\n",
           ok1, ok5, ok3, ok4, ok6, (int)backups.size(), ok7, all ? "PASS" : "FAIL");
    fs::remove(save, ec); fs::remove_all(PlayerRecord_Dir(), ec); fs::remove_all(save + ".regions", ec); fs::remove_all(BackupDir(), ec);
    SetWorldSavePath(oldPath.c_str());
    return all;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#869): the player record and the world file's extras (format 13); the migration chain and refusals (#1693) ----
// The world file's own header and bodies need a Scene (GL); --test-save-split covers those. These cover what the hooks
// read and write, which is where the player and the teams live.
#include "unit_helpers.h"

namespace {
std::string Slurp(const std::string& p) {
    std::string s;
    if (FILE* f = fopen(p.c_str(), "rb")) { char b[4096]; size_t n; while ((n = fread(b, 1, sizeof b, f)) > 0) s.append(b, n); fclose(f); }
    return s;
}
void Spit(const std::string& p, const std::string& s) { if (FILE* f = fopen(p.c_str(), "wb")) { fwrite(s.data(), 1, s.size(), f); fclose(f); } }
int Lines(const std::string& s, const char* prefix) {
    int n = 0; size_t at = 0; const size_t len = strlen(prefix);
    while (at < s.size()) { const size_t e = s.find('\n', at); if (s.compare(at, len, prefix) == 0) ++n; at = e == std::string::npos ? s.size() : e + 1; }
    return n;
}

struct SaveTest : ::testing::Test {
    std::string dir, oldSave;
    void SetUp() override {
        UnitTest::EnsureItemDB();
        dir = UnitTest::ScratchDir();
        oldSave = WorldSavePath();
        SetWorldSavePath((dir + "/world.sav").c_str());
        GameSave_Register();
        Teams_Reset();
    }
    void TearDown() override { SetWorldSavePath(oldSave.c_str()); Teams_Reset(); inventory.clear(); }
    // The extras as the save writes them, read back by the load hook (an "end_extras" line closes them, as in the file).
    std::string WriteExtras() {
        const std::string p = dir + "/extras.txt";
        FILE* f = fopen(p.c_str(), "w");
        saveExtrasHook(f);
        fprintf(f, "end_extras\n");
        fclose(f);
        return Slurp(p);
    }
    void ReadExtras(const std::string& text) {
        const std::string p = dir + "/extras_in.txt";
        Spit(p, text);
        FILE* f = fopen(p.c_str(), "r");
        ASSERT_NE(f, nullptr);
        loadExtrasHook(f);
        fclose(f);
    }
};
}

// A server's record of a remote player: position, inventory (with heat), vitals and the lines it doesn't hold (kept text).
TEST_F(SaveTest, PlayerRecordRoundTrip) {
    Player p; p.position = { 1234.5, 60001.25, -77.125 }; p.forward = { 0, 0, 1 }; p.jetpackFuel = 0.4f;
    Inventory inv; inv.add("iron_ingot", 7); inv.addAmount("iron_ore", 2.5, 0.75f); inv.add("stone_pickaxe", 1);
    SurvivalState sv; sv.health = 63.0f; sv.food = 41.0f; sv.oxygen = 80.0f; sv.deaths = 3;
    sv.hasBed = true; sv.bedBody = "terra"; sv.bedPlot = 4; sv.bedCell[0] = -2; sv.bedCell[1] = 5; sv.bedCell[2] = 9;
    PlayerRecordMeta meta; meta.uuid = PlayerId_FromName("unit"); meta.name = "unit"; meta.firstSeen = 1700000000.0; meta.playTime = 60.0;
    meta.kept = "hotbar 0 45\nmap_pin 1 2 3 home\n";
    const std::string path = PlayerRecord_Path(meta.uuid);
    ASSERT_TRUE(PlayerRecord_Write(path, p, inv, sv, meta, 30.0));
    EXPECT_FALSE(std::filesystem::exists(path + ".tmp"));
    const std::string text = Slurp(path);
    EXPECT_EQ(text.rfind("RAYLIBGAME_PLAYER 1", 0), 0u);

    Player p2; Inventory inv2; SurvivalState sv2; PlayerRecordMeta m2;
    ASSERT_TRUE(PlayerRecord_Read(path, p2, inv2, sv2, m2));
    EXPECT_EQ(m2.uuid, meta.uuid);
    EXPECT_EQ(m2.name, "unit");
    EXPECT_DOUBLE_EQ(m2.playTime, 90.0);
    EXPECT_LT((p2.position - p.position).length(), 1e-5);
    EXPECT_FLOAT_EQ(p2.jetpackFuel, 0.4f);
    EXPECT_EQ(inv2.count("iron_ingot"), 7);
    EXPECT_DOUBLE_EQ(inv2.amount("iron_ore"), 2.5);
    EXPECT_EQ(inv2.count("stone_pickaxe"), 1);
    for (int i = 0; i < Inventory::kCells; ++i) {
        EXPECT_EQ(inv2.grid()[i].id, inv.grid()[i].id) << "cell " << i;
        EXPECT_NEAR(inv2.grid()[i].heat, inv.grid()[i].heat, 1e-3f) << "cell " << i;
    }
    EXPECT_FLOAT_EQ(sv2.health, 63.0f); EXPECT_FLOAT_EQ(sv2.food, 41.0f); EXPECT_FLOAT_EQ(sv2.oxygen, 80.0f); EXPECT_EQ(sv2.deaths, 3);
    EXPECT_TRUE(sv2.hasBed); EXPECT_EQ(sv2.bedBody, "terra"); EXPECT_EQ(sv2.bedPlot, 4u); EXPECT_EQ(sv2.bedCell[0], -2); EXPECT_EQ(sv2.bedCell[2], 9);
    EXPECT_NE(m2.kept.find("hotbar 0 45"), std::string::npos) << "a line the server doesn't hold is kept: " << m2.kept;
    EXPECT_NE(m2.kept.find("map_pin 1 2 3 home"), std::string::npos) << m2.kept;
}

// (#920) A remote player's skills in their record: read into their SkillState, written from it, never kept as text too.
TEST_F(SaveTest, PlayerRecordSkillsRoundTrip) {
    ASSERT_TRUE(Skills_Load("resources/data/skills.txt"));
    Player p; p.forward = { 0, 0, 1 };
    Inventory inv; SurvivalState sv;
    PlayerRecordMeta meta; meta.uuid = PlayerId_FromName("skilled"); meta.name = "skilled"; meta.kept = "hotbar 0 45\n";
    SkillState sk; sk.xp["mining"] = 2478.5; sk.xp["digging"] = 12.25; sk.ranks["quick_pick"] = 2;
    const std::string path = PlayerRecord_Path(meta.uuid);
    ASSERT_TRUE(PlayerRecord_Write(path, p, inv, sv, meta, 0.0, &sk));
    Player p2; Inventory inv2; SurvivalState sv2; PlayerRecordMeta m2; SkillState sk2;
    ASSERT_TRUE(PlayerRecord_Read(path, p2, inv2, sv2, m2, &sk2));
    EXPECT_NEAR(sk2.xp["mining"], 2478.5, 1e-3);
    EXPECT_NEAR(sk2.xp["digging"], 12.25, 1e-3);
    EXPECT_EQ(sk2.ranks["quick_pick"], 2);
    EXPECT_EQ(m2.kept.find("skill_xp"), std::string::npos) << "the skills are not kept as text as well: " << m2.kept;
    EXPECT_EQ(m2.kept.find("perk"), std::string::npos) << m2.kept;
    // Written again (the server's autosave) without losing or doubling them.
    ASSERT_TRUE(PlayerRecord_Write(path, p2, inv2, sv2, m2, 0.0, &sk2));
    SkillState sk3; PlayerRecordMeta m3; Player p3; Inventory inv3; SurvivalState sv3;
    ASSERT_TRUE(PlayerRecord_Read(path, p3, inv3, sv3, m3, &sk3));
    EXPECT_NEAR(sk3.xp["mining"], 2478.5, 1e-3);
    EXPECT_EQ(sk3.ranks["quick_pick"], 2);
}

TEST_F(SaveTest, TeamLinesRoundTrip) {
    const std::string a = PlayerId_FromName("alice"), b = PlayerId_FromName("bob");
    const int t = Teams_Join(a, "alice's team", 0xD04A2F);
    Teams_SetOwner(a);
    Teams_Move(b, t);
    Teams_SetQuests(t, "quest 2\nside_done e0f\n");
    const std::string p = dir + "/teams.txt";
    { FILE* f = fopen(p.c_str(), "w"); Teams_Save(f, "", nullptr); fclose(f); }
    Teams_Reset();
    FILE* f = fopen(p.c_str(), "r");
    char tok[64];
    while (fscanf(f, " %63s", tok) == 1) ASSERT_TRUE(Teams_LoadLine(tok, f, "")) << tok;
    fclose(f);
    ASSERT_EQ(Teams_All().size(), 1u);
    const TeamRecord& r = Teams_All()[0];
    EXPECT_EQ(r.id, t);
    EXPECT_EQ(r.colour, 0xD04A2Fu);
    EXPECT_EQ(r.name, "alice's_team") << "spaces become underscores";
    EXPECT_EQ(r.members, (std::vector<std::string>{ a, b }));
    EXPECT_EQ(Teams_OwnerOf(t), a);
    EXPECT_EQ(Teams_Owner(), a);
    EXPECT_EQ(r.quests, "quest 2\nside_done e0f\n");
}

// Format 13: the world's extras hold the team and no player lines; the player's record beside it brings them back.
TEST_F(SaveTest, Format13PutsThePlayerInItsRecord) {
    inventory.clear(); inventory.add("iron_ingot", 7);
    survival.health = 63.0f; survival.food = 41.0f;
    const std::string w = WriteExtras();
    EXPECT_EQ(Lines(w, "slot "), 0);
    EXPECT_EQ(Lines(w, "player "), 0);
    EXPECT_EQ(Lines(w, "survival "), 0);
    EXPECT_EQ(Lines(w, "team "), 1);
    EXPECT_EQ(Lines(w, "owner "), 1);
    const std::string rec = Slurp(PlayerRecord_Path(PlayerId_Local()));
    EXPECT_GE(Lines(rec, "slot "), 1);
    EXPECT_EQ(Lines(rec, "survival "), 1);

    inventory.clear(); survival.health = 100.0f; survival.food = 100.0f;
    ReadExtras(w);
    EXPECT_TRUE(GameSave_PlayerLoaded());
    EXPECT_EQ(inventory.count("iron_ingot"), 7);
    EXPECT_FLOAT_EQ(survival.health, 63.0f);
    EXPECT_FLOAT_EQ(survival.food, 41.0f);
}

// (#1693) Versions: this build opens 13 only; a newer build's save and one from before 13 are refused, each with its message.
TEST(SaveVersions, OnlyTheSupportedRangeOpens) {
    std::string why;
    EXPECT_TRUE(Save_VersionSupported(kSaveVersion, why));
    EXPECT_FALSE(Save_VersionSupported(kSaveVersion + 1, why));
    EXPECT_NE(why.find("newer build"), std::string::npos) << why;
    EXPECT_NE(why.find(std::to_string(kSaveVersion + 1)), std::string::npos) << why;
    EXPECT_FALSE(Save_VersionSupported(kOldestSaveVersion - 1, why));
    EXPECT_NE(why.find("older build"), std::string::npos) << why;
}

// The chain over text lines: the current version needs no step; one with no step to it fails with a message; the lines come
// back with the version line unchanged.
TEST(SaveVersions, TheChainRunsFromAnOlderVersionToTheCurrentOne) {
    std::vector<std::string> lines = { "RAYLIBGAME_SAVE " + std::to_string(kSaveVersion), "seed 42" };
    const std::vector<std::string> before = lines;
    MigrateContext ctx;
    std::string error;
    EXPECT_TRUE(Save_Migrate(kSaveVersion, lines, ctx, error)) << error;
    EXPECT_EQ(lines, before);
    EXPECT_TRUE(ctx.dropped.empty());
    EXPECT_EQ(Save_UpgradeNote(ctx), "Save upgraded from version " + std::to_string(kSaveVersion));
    ctx.drop("a"); ctx.drop("b");
    EXPECT_EQ(Save_UpgradeNote(ctx), "Save upgraded from version " + std::to_string(kSaveVersion) + "; dropped: a, b");
    if (kSaveVersion == kOldestSaveVersion) {
        MigrateContext c2;
        EXPECT_FALSE(Save_Migrate(kSaveVersion - 1, lines, c2, error)) << "no step from a version before the oldest";
        EXPECT_FALSE(error.empty());
    }
}

// A record from a newer build is not read (its next save keeps it as a .bak, #1110); a current one is.
TEST_F(SaveTest, APlayerRecordFromANewerBuildIsNotRead) {
    const std::string path = dir + "/newer.player";
    Spit(path, "RAYLIBGAME_PLAYER " + std::to_string(kPlayerRecordVersion + 1) + "\nuuid x\nend_record\n");
    Player p; Inventory inv; SurvivalState sv; PlayerRecordMeta m;
    EXPECT_FALSE(PlayerRecord_Read(path, p, inv, sv, m));
    Spit(path, "RAYLIBGAME_PLAYER " + std::to_string(kPlayerRecordVersion) + "\nuuid x\nend_record\n");
    EXPECT_TRUE(PlayerRecord_Read(path, p, inv, sv, m));
}
#endif
