#include "teamquests.h"
#include "survival.h"
#include "questtree.h"
#include "playerrecord.h"

#include <cstdio>
#include <cstring>
#include <raylib.h>
#include <map>

namespace {
    std::map<int, QuestBookState> books; // loaded team books (by team id)

    // The book selected while `fn` runs (every QuestBook_* call works on the selected one), then the local one again.
    template <class F> void With(QuestBookState& b, F fn) {
        QuestBookState* before = &QuestBook_Selected();
        QuestBook_Select(&b);
        fn();
        QuestBook_Select(before == &QuestBook_LocalState() ? nullptr : before);
    }

    // Save lines <-> a book, through a scratch file (the quest book reads and writes FILE*, as the world file does).
    // (tmpfile() makes its file in the drive's root on Windows, which often fails: a scratch file under saves/ instead)
    const char* kScratch = "saves/.teamquests.tmp";
    FILE* Scratch(const char* mode) {
        if (!DirectoryExists("saves")) MakeDirectory("saves");
        return fopen(kScratch, mode);
    }
    void Parse(const std::string& text) {
        if (text.empty()) return;
        FILE* f = Scratch("wb");
        if (!f) return;
        fwrite(text.data(), 1, text.size(), f);
        fclose(f);
        if (!(f = Scratch("rb"))) return;
        char tok[64], rest[512];
        while (fscanf(f, " %63s", tok) == 1)
            if (!QuestBook_LoadLine(tok, f) && !fgets(rest, sizeof rest, f)) break; // an unknown line: skipped
        fclose(f);
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
        for (const TeamRecord& t : Teams_All()) if (t.id == team) { Parse(t.quests); break; }
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
    for (auto& [team, b] : books) { std::string s; With(b, [&] { s = Write(); }); Teams_SetQuests(team, s); }
}

void TeamQuests_Reset() { books.clear(); }

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
    printf("bench: team quests: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}
