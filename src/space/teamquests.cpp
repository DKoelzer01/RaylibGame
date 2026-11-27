#include "teamquests.h"
#include "survival.h"
#include "questtree.h"
#include "playerrecord.h"
#include "navigation.h"
#include <unordered_set>

#include <cstdio>
#include <cstring>
#ifdef _WIN32
#include <process.h> // _getpid
#else
#include <unistd.h>
#define _getpid getpid // (#1025) the Linux server build
#endif
#include <string>
#include <raylib.h>
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

namespace {
    std::map<int, QuestBookState> books; // loaded team books (by team id)
    // (#338 slice 3) the map tiles each team has explored, by body: saved with its book (team_tiles lines), never sent to clients
    std::map<int, std::map<std::string, std::unordered_set<uint64_t>>> tiles;
    // (#1886) the map tiles each player revealed themselves, by player id: a player who changes teams carries them along
    // (the user, on #1886); saved in the world (player_tiles blocks)
    std::map<std::string, NavTiles> playerTiles;
    // (#1886) each team's shared map pins (players place them from their maps: MsgMapPin); saved by the team store
    // (team_pins blocks) and sent to the team in its "pins" section
    std::map<int, std::vector<NavPin>> pins;

    // The book selected while `fn` runs (every QuestBook_* call works on the selected one), then the local one again.
    template <class F> void With(QuestBookState& b, F fn) {
        QuestBookState* before = &QuestBook_Selected();
        QuestBook_Select(&b);
        fn();
        QuestBook_Select(before == &QuestBook_LocalState() ? nullptr : before);
    }

    // Save lines <-> a book, through a scratch file (the quest book reads and writes FILE*, as the world file does).
    // (tmpfile() makes its file in the drive's root on Windows, which often fails: a scratch file under saves/ instead)
    // (#837) per process: a host and its clients on one machine share the working directory, and one name raced
    const std::string kScratch = "saves/.teamquests." + std::to_string(_getpid()) + ".tmp";
    FILE* Scratch(const char* mode) {
        if (!DirectoryExists("saves")) MakeDirectory("saves");
        return fopen(kScratch.c_str(), mode);
    }
    void Parse(const std::string& text, std::map<std::string, std::unordered_set<uint64_t>>* teamTiles = nullptr) {
        if (text.empty()) return;
        FILE* f = Scratch("wb");
        if (!f) return;
        fwrite(text.data(), 1, text.size(), f);
        fclose(f);
        if (!(f = Scratch("rb"))) return;
        char tok[64], rest[512];
        while (fscanf(f, " %63s", tok) == 1) {
            if (!strcmp(tok, "team_tiles")) { // (#338 slice 3) team_tiles <body> <n> <hex>... (only a team's stored block has them)
                char body[64] = {}; int n = 0;
                if (fscanf(f, " %63s %d", body, &n) != 2) break;
                std::unordered_set<uint64_t> dummy;
                auto& set = teamTiles ? (*teamTiles)[body] : dummy;
                for (int k = 0; k < n; ++k) { unsigned long long key = 0; if (fscanf(f, " %llx", &key) != 1) break; set.insert(key); }
                continue;
            }
            if (!QuestBook_LoadLine(tok, f) && !fgets(rest, sizeof rest, f)) break; // an unknown line: skipped
        }
        fclose(f);
        std::remove(kScratch.c_str());
    }
    std::string Write() {
        FILE* f = Scratch("wb");
        if (!f) return std::string();
        QuestBook_Save(f);
        fclose(f);
        std::string out;
        if (!(f = Scratch("rb"))) return out;
        char buf[4096];
        size_t n;
        while ((n = fread(buf, 1, sizeof buf, f)) > 0) out.append(buf, n);
        fclose(f);
        std::remove(kScratch.c_str());
        return out;
    }
}

QuestBookState& TeamQuests_Book(int team) {
    auto it = books.find(team);
    if (it != books.end()) return it->second;
    QuestBookState& b = books[team];
    b.ownObtained = true;
    With(b, [&] {
        QuestBook_Reset();
        for (const TeamRecord& t : Teams_All()) if (t.id == team) { Parse(t.quests, &tiles[team]); break; }
    });
    return b;
}

void TeamQuests_Count(int team, const std::string& event, double n) { With(TeamQuests_Book(team), [&] { QuestBook_Count(event, n); }); }
void TeamQuests_AtLeast(int team, const std::string& event, double v) { With(TeamQuests_Book(team), [&] { QuestBook_AtLeast(event, v); }); }
void TeamQuests_NotePlaced(int team, const std::string& item) { With(TeamQuests_Book(team), [&] { QuestBook_NotePlaced(item); }); }
void TeamQuests_Obtained(int team, const std::string& item, double n) { if (n > 0.0) TeamQuests_Book(team).obtained[item] += n; }

void TeamQuests_Update(int team, const std::function<void(const std::string&, int)>& reward, const std::function<void(const std::string&)>& news) {
    QuestBookState& b = TeamQuests_Book(team);
    auto oldReward = gQuestReward; auto oldNews = gQuestNews;
    gQuestReward = reward ? reward : [](const std::string&, int) {}; // never the local inventory: a server's own player isn't on this team
    gQuestNews = news ? news : [](const std::string&) {};
    With(b, [&] { QuestBook_Update(); });
    gQuestReward = oldReward; gQuestNews = oldNews;
}

std::string TeamQuests_Text(int team) {
    std::string s;
    With(TeamQuests_Book(team), [&] { s = Write(); });
    return s;
}

void TeamQuests_StoreAll() {
    for (auto& [team, b] : books) {
        std::string s; With(b, [&] { s = Write(); });
        for (const auto& [body, set] : tiles[team]) { // (#338 slice 3)
            s += "team_tiles " + body + " " + std::to_string(set.size());
            int k = 0;
            for (uint64_t key : set) { char hex[24]; snprintf(hex, sizeof hex, " %llx", (unsigned long long)key); s += hex; if (++k % 16 == 0 && k < (int)set.size()) s += "\n"; }
            s += "\n";
        }
        Teams_SetQuests(team, s);
    }
}

void TeamQuests_Reset() { books.clear(); tiles.clear(); playerTiles.clear(); pins.clear(); }

std::string TeamQuests_TilesText(int team) {
    TeamQuests_Book(team); // (its stored tiles load with it)
    return Nav_TilesText(tiles[team]);
}

int TeamQuests_CarryPlayer(const std::string& uuid, int team) {
    auto it = playerTiles.find(uuid);
    if (it == playerTiles.end() || !team) return 0;
    QuestBookState& dst = TeamQuests_Book(team);
    size_t added = 0, all = 0;
    for (const auto& [body, set] : it->second) {
        auto& into = tiles[team][body];
        const size_t before = into.size();
        into.insert(set.begin(), set.end());
        added += into.size() - before;
    }
    for (const auto& [body, set] : tiles[team]) all += set.size();
    if (all) dst.counts["explored"] = std::max(dst.counts["explored"], (double)all);
    return (int)added;
}

void TeamQuests_SavePlayers(FILE* f) {
    for (const auto& [uuid, t] : playerTiles) {
        const std::string text = Nav_TilesText(t);
        if (!text.empty()) fprintf(f, "player_tiles %s\n%splayer_tiles_end\n", uuid.c_str(), text.c_str());
    }
}

bool TeamQuests_LoadPlayerLine(const char* tok, FILE* f) {
    if (strcmp(tok, "player_tiles") != 0) return false;
    char uuid[64] = {};
    if (fscanf(f, " %63s", uuid) != 1) return true;
    std::string text;
    char word[64];
    while (fscanf(f, " %63s", word) == 1 && strcmp(word, "player_tiles_end") != 0) { text += word; text += ' '; } // (the tile text reads as words)
    Nav_TilesFromText(text, playerTiles[uuid]);
    return true;
}

void TeamQuests_ApplyText(const std::string& text) {
    QuestBook_Select(nullptr);
    QuestBook_Reset();
    Parse(text); // into the local book ("obtained" lines into the local inventory's totals, which the HUD reads)
}

// --test-team-quests (#338): two teams' books on a "server": a team's events finish its node (e0f: eat once) and pay its
// reward through the callback; the other team's and the local book are untouched; a store and reload keeps the progress.
bool TeamQuests_SelfTest() {
    bool ok = true;
    auto check = [&](bool c, const char* what) { printf("bench: team quests: %s: %s\n", what, c ? "yes" : "NO"); ok = ok && c; };
    TeamQuests_Reset();
    const int a = Teams_Join(PlayerId_New(), "red", 0xD04A2F), b = Teams_Join(PlayerId_New(), "blue", 0x3C8CDC);
    const bool localBefore = QuestTree_IsDone("e0f");
    std::string paid; int paidN = 0; std::string heard;
    auto reward = [&](const std::string& item, int n) { paid = item; paidN += n; };
    auto news = [&](const std::string& t) { heard = t; };
    TeamQuests_Update(a, reward, news);
    check(paid.empty(), "nothing paid before the event");
    TeamQuests_Count(a, "eat", 1.0);
    TeamQuests_Update(a, reward, news);
    bool doneA = false, doneB = false;
    With(TeamQuests_Book(a), [&] { doneA = QuestTree_IsDone("e0f"); });
    With(TeamQuests_Book(b), [&] { doneB = QuestTree_IsDone("e0f"); });
    check(doneA && paid == "berries" && paidN == 4, "team a ate once: its node is done and the berries are paid through the server");
    check(!heard.empty(), "and its news is told");
    check(!doneB, "team b's book is untouched");
    check(QuestTree_IsDone("e0f") == localBefore, "the local book is untouched");
    TeamQuests_Update(a, reward, news);
    check(paidN == 4, "a node pays once");
    TeamQuests_Obtained(a, "oak_log", 3.0);
    TeamQuests_StoreAll();
    TeamQuests_Reset();
    bool again = false; double logs = 0.0;
    With(TeamQuests_Book(a), [&] { again = QuestTree_IsDone("e0f"); logs = QuestBook_Selected().obtained["oak_log"]; });
    check(again && logs > 2.99, "stored and reloaded: the node and the team's obtained totals kept");
    // (#338 slice 3) exploring and visits: two spots 600 m apart on a 60 km body add up for the team, the same spot twice
    // adds nothing, the tiles survive a store and reload; a moon visit sets only that team's flag
    const Vector3 p1 = { 0.0f, 1.0f, 0.0f }, p2 = Vector3Normalize(Vector3{ 0.01f, 1.0f, 0.0f });
    const int e1 = TeamQuests_Explore(a, "terra", p1, 60000.0f), e2 = TeamQuests_Explore(a, "terra", p2, 60000.0f), e3 = TeamQuests_Explore(a, "terra", p1, 60000.0f);
    check(e1 > 0 && e2 > e1 && e3 == 0 && TeamQuests_Explored(a) == (double)e2, "exploring adds up for the team; the same ground twice adds nothing");
    TeamQuests_Visit(a, true);
    bool moonA = false, moonB = false;
    With(TeamQuests_Book(a), [&] { moonA = QuestBook_Met("visit", "moon", 1.0); });
    With(TeamQuests_Book(b), [&] { moonB = QuestBook_Met("visit", "moon", 1.0); });
    check(moonA && !moonB, "a moon visit counts for that team only");
    TeamQuests_StoreAll();
    TeamQuests_Reset();
    check(TeamQuests_Explore(a, "terra", p1, 60000.0f) == 0 && TeamQuests_Explore(a, "terra", p2, 60000.0f) == 0, "stored and reloaded: the team's tiles kept");
    // (#861) a joiner's book merged into a team's: the higher tier, the union of done tasks / nodes / visits / tiles, the
    // larger count (not the sum); what the team had stays
    {
        const int solo = 9001, team = 9002;
        TeamQuests_Forget(solo); TeamQuests_Forget(team);
        QuestBookState& s = TeamQuests_Book(solo);
        s.tier = 2; s.tasksDone.insert({ 1, 0 }); s.treeDone.insert("e0f"); s.counts["eat"] = 3.0; s.obtained["oak_log"] = 5.0; s.visitedMoon = true;
        QuestBookState& t = TeamQuests_Book(team);
        t.tier = 1; t.tasksDone.insert({ 0, 2 }); t.treeDone.insert("e1a"); t.counts["eat"] = 5.0; t.obtained["oak_log"] = 2.0;
        TeamQuests_Explore(solo, "terra", p1, 60000.0f);
        const int teamTiles = TeamQuests_Explore(team, "terra", p2, 60000.0f);
        TeamQuests_Merge(solo, team);
        const QuestBookState& m = TeamQuests_Book(team);
        check(m.tier == 2 && m.tasksDone.count({ 1, 0 }) && m.tasksDone.count({ 0, 2 }) && m.treeDone.count("e0f") && m.treeDone.count("e1a") &&
              m.counts.at("eat") == 5.0 && m.obtained.at("oak_log") == 5.0 && m.visitedMoon, "merged: tier, tasks, nodes, visits; counts take the larger");
        check(TeamQuests_Explored(team) > teamTiles && TeamQuests_Explore(team, "terra", p1, 60000.0f) == 0, "merged: the joiner's explored tiles are the team's");
        TeamQuests_Forget(solo); TeamQuests_Forget(team);
    }
    printf("bench: team quests: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

int TeamQuests_Explore(int team, const std::string& body, Vector3 up, float radius, const std::string& uuid) {
    TeamQuests_Book(team); // (its stored tiles load with it)
    if (!uuid.empty()) Nav_RevealTiles(up, radius, playerTiles[uuid][body]); // (#1886) theirs too, to carry on a team change
    if (Nav_RevealTiles(up, radius, tiles[team][body]) == 0) return 0;
    size_t all = 0;
    for (const auto& [b, set] : tiles[team]) all += set.size();
    TeamQuests_AtLeast(team, "explored", (double)all);
    return (int)all;
}

void TeamQuests_Visit(int team, bool moon) {
    QuestBookState& b = TeamQuests_Book(team);
    if (moon) b.visitedMoon = true; else b.visitedRing = true;
}

double TeamQuests_Explored(int team) {
    double n = 0.0;
    With(TeamQuests_Book(team), [&] { n = QuestBook_CountOf("explored"); });
    return n;
}

void TeamQuests_Copy(int from, int to) {
    if (from == to) return;
    QuestBookState copy = TeamQuests_Book(from);
    TeamQuests_Book(to) = copy;
    tiles[to] = tiles[from];
    pins[to] = pins[from];
}

void TeamQuests_Merge(int from, int into) {
    if (from == into) return;
    const QuestBookState src = TeamQuests_Book(from); // (a copy: Book(into) may add to the map)
    QuestBookState& dst = TeamQuests_Book(into);
    dst.tier = std::max(dst.tier, src.tier);
    dst.placed.insert(src.placed.begin(), src.placed.end());
    dst.visitedBodies.insert(src.visitedBodies.begin(), src.visitedBodies.end());
    dst.reachedOrbit |= src.reachedOrbit; dst.visitedMoon |= src.visitedMoon; dst.visitedRing |= src.visitedRing;
    dst.tasksDone.insert(src.tasksDone.begin(), src.tasksDone.end());
    dst.treeDone.insert(src.treeDone.begin(), src.treeDone.end());
    // Counts and totals take the larger side, not the sum: both sides' work toward one task would count twice.
    for (const auto& [k, v] : src.counts) dst.counts[k] = std::max(dst.counts[k], v);
    for (const auto& [k, v] : src.obtained) dst.obtained[k] = std::max(dst.obtained[k], v);
    size_t all = 0;
    for (const auto& [body, set] : tiles[from]) tiles[into][body].insert(set.begin(), set.end());
    for (const auto& [body, set] : tiles[into]) all += set.size();
    if (all) dst.counts["explored"] = std::max(dst.counts["explored"], (double)all);
    const std::vector<NavPin> carried = pins[from]; // (#1886) the joiner's pins come along, as their map does
    for (const NavPin& p : carried) TeamQuests_SetPin(into, true, p.body, p.dir, p.name);
}

void TeamQuests_Forget(int team) { books.erase(team); tiles.erase(team); pins.erase(team); }

// (#1886) The team's shared map pins. A name already on that body gets a number added (two members' "WP1"s both stay).
bool TeamQuests_SetPin(int team, bool add, const std::string& body, Vector3 dir, const std::string& rawName) {
    const std::string name = Nav_PinName(rawName);
    if (!team || body.empty() || Nav_PinName(body) != body || name.empty()) return false; // (one-word names: the text and the save)
    std::vector<NavPin>& list = pins[team];
    auto at = [&](const std::string& n) {
        return std::find_if(list.begin(), list.end(), [&](const NavPin& p) { return p.body == body && p.name == n; });
    };
    if (!add) {
        auto it = at(name);
        if (it == list.end()) return false;
        list.erase(it);
        return true;
    }
    if (!std::isfinite(dir.x) || !std::isfinite(dir.y) || !std::isfinite(dir.z) || Vector3Length(dir) < 0.5f) return false;
    if ((int)list.size() >= kMaxTeamPins) return false;
    std::string unique = name;
    for (int k = 2; at(unique) != list.end(); ++k) unique = name.substr(0, 44) + "_" + std::to_string(k);
    list.push_back({ body, unique, Vector3Normalize(dir) });
    return true;
}

std::string TeamQuests_PinsText(int team) { auto it = pins.find(team); return it == pins.end() ? std::string() : Nav_PinsText(it->second); }
void TeamQuests_SetPinsText(int team, const std::string& text) { pins[team] = Nav_PinsFromText(text); }
int TeamQuests_PinCount(int team) { auto it = pins.find(team); return it == pins.end() ? 0 : (int)it->second.size(); }

int TeamQuests_BodiesReached(int team) { // (#343) planets and moons with explored ground (only PlanetField bodies explore)
    TeamQuests_Book(team);
    int n = 0;
    for (const auto& [body, set] : tiles[team]) if (!set.empty()) ++n;
    return n;
}
