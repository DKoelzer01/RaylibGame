#include "playerrecord.h"
#include "menus.h"
#include "survival.h"
#include "world.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <random>

namespace {
    std::vector<TeamRecord> teams;
    std::string owner;
    TeamRecord* Find(int id) { for (TeamRecord& t : teams) if (t.id == id) return &t; return nullptr; }
}

std::string PlayerId_New() {
    std::random_device rd;
    uint8_t b[16];
    for (int i = 0; i < 16; i += 4) { const uint32_t r = rd(); memcpy(b + i, &r, 4); }
    b[6] = (uint8_t)((b[6] & 0x0F) | 0x40); // version 4 (random)
    b[8] = (uint8_t)((b[8] & 0x3F) | 0x80); // the RFC 4122 variant
    char s[37];
    snprintf(s, sizeof(s), "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x",
             b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]);
    return s;
}

bool PlayerId_Valid(const std::string& id) {
    if (id.size() != 36) return false;
    for (size_t i = 0; i < id.size(); ++i) {
        const bool dash = i == 8 || i == 13 || i == 18 || i == 23;
        if (dash ? id[i] != '-' : !isxdigit((unsigned char)id[i])) return false;
    }
    return true;
}

std::string PlayerId_Local() {
    // Settings_Load makes and saves it on the first start; a bench (no settings read) gets one for the run.
    if (!PlayerId_Valid(settings.playerUuid)) settings.playerUuid = PlayerId_New();
    return settings.playerUuid;
}

std::string PlayerRecord_Dir() {
    std::string s = WorldSavePath();
    return s + ".players"; // beside <save>.regions
}
std::string PlayerRecord_Path(const std::string& uuid) { return PlayerRecord_Dir() + "/" + uuid + ".player"; }

bool Save_WriteAtomic(const std::string& path, const std::function<bool(FILE*)>& write) {
    const std::string tmp = path + ".tmp";
    FILE* f = fopen(tmp.c_str(), "w");
    if (!f) return false;
    bool ok = write(f);
    ok = fflush(f) == 0 && ok;
    ok = fclose(f) == 0 && ok;
    std::error_code ec;
    if (ok) std::filesystem::rename(tmp, path, ec); // replaces the old file (MoveFileEx with replace on Windows)
    if (!ok || ec) { std::filesystem::remove(tmp, ec); return false; }
    return true;
}

void Teams_Reset() { teams.clear(); owner.clear(); }
const std::vector<TeamRecord>& Teams_All() { return teams; }
const std::string& Teams_Owner() { return owner; }
void Teams_SetOwner(const std::string& uuid) { owner = uuid; }

int Teams_Of(const std::string& uuid) {
    for (const TeamRecord& t : teams) for (const std::string& m : t.members) if (m == uuid) return t.id;
    return 0;
}

int Teams_Join(const std::string& uuid, const std::string& name, uint32_t colour) {
    if (const int t = Teams_Of(uuid)) return t;
    int id = 1;
    for (const TeamRecord& t : teams) id = std::max(id, t.id + 1);
    TeamRecord t;
    t.id = id; t.colour = colour & 0xFFFFFF; t.name = name.empty() ? "team" : name; t.members.push_back(uuid);
    teams.push_back(t);
    return id;
}

void Teams_ReplaceMember(const std::string& from, const std::string& to) {
    for (TeamRecord& t : teams) for (std::string& m : t.members) if (m == from) m = to;
    if (owner == from) owner = to;
}

bool Teams_LoadLine(const char* tok, FILE* f, const std::string& localId) {
    char word[128] = {};
    if (strcmp(tok, "owner") == 0) { if (fscanf(f, " %127s", word) == 1 && PlayerId_Valid(word)) owner = word; return true; }
    if (strcmp(tok, "team") == 0) { // team <id> <colour hex> <name>
        int id = 0; unsigned colour = 0;
        if (fscanf(f, " %d %x %127s", &id, &colour, word) == 3 && id > 0 && !Find(id)) {
            TeamRecord t; t.id = id; t.colour = colour & 0xFFFFFF; t.name = word;
            teams.push_back(t);
        }
        return true;
    }
    if (strcmp(tok, "team_member") == 0) {
        int id = 0;
        if (fscanf(f, " %d %127s", &id, word) == 2 && PlayerId_Valid(word))
            if (TeamRecord* t = Find(id)) if (std::find(t->members.begin(), t->members.end(), word) == t->members.end()) t->members.push_back(word);
        return true;
    }
    if (strcmp(tok, "team_quests") == 0) {
        // This game's team (its member, else the owner's team when the local player has none: single player opening a world
        // made on another PC continues as its owner): into the quest book. Any other: kept as text.
        int id = 0;
        fscanf(f, " %d", &id);
        const int mine = localId.empty() ? -1 : (Teams_Of(localId) ? Teams_Of(localId) : Teams_Of(owner));
        char line[512];
        if (!fgets(line, sizeof(line), f)) return true; // the rest of the team_quests line
        TeamRecord* t = Find(id);
        if (id == mine) {
            char q[64] = {};
            while (fscanf(f, " %63s", q) == 1 && strcmp(q, "team_quests_end") != 0)
                if (!QuestBook_LoadLine(q, f) && !fgets(line, sizeof(line), f)) break; // an unknown line: skipped
        } else {
            std::string text;
            while (fgets(line, sizeof(line), f)) {
                const char* s = line; while (*s == ' ' || *s == '\t') ++s;
                if (strncmp(s, "team_quests_end", 15) == 0) break;
                text += line;
            }
            if (t) t->quests = text;
        }
        return true;
    }
    return false;
}

void Teams_Save(FILE* f, const std::string& localId, const std::function<void(FILE*)>& writeLocalQuests) {
    if (!owner.empty()) fprintf(f, "owner %s\n", owner.c_str());
    const int mine = localId.empty() ? -1 : Teams_Of(localId);
    for (const TeamRecord& t : teams) {
        std::string name = t.name;
        for (char& c : name) if (isspace((unsigned char)c)) c = '_';
        fprintf(f, "team %d %06x %s\n", t.id, t.colour, name.empty() ? "team" : name.c_str());
        for (const std::string& m : t.members) fprintf(f, "team_member %d %s\n", t.id, m.c_str());
        if (t.id == mine && writeLocalQuests) { fprintf(f, "team_quests %d\n", t.id); writeLocalQuests(f); fprintf(f, "team_quests_end\n"); }
        else if (!t.quests.empty()) { fprintf(f, "team_quests %d\n%s", t.id, t.quests.c_str()); fprintf(f, "team_quests_end\n"); }
    }
}
