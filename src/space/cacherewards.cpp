// (#2072) A cache's team-side rewards on the host: research units into a node, alternate recipes into the team's list.
// The team store sends both to the team's members when they change (hostserver's TeamStoreTick).
#include "caches.h"
#include "research.h"

#include <algorithm>

void Caches_ApplyTeam(int teamId, const std::vector<CacheRoll>& rolls) {
    for (const CacheRoll& r : rolls) {
        if (r.kind == CacheReward::Alt) Caches_AddTeamAlt(teamId, r.id);
        if (r.kind != CacheReward::Boost || r.units <= 0) continue;
        const ResearchNode* n = Research_Find(r.id);
        ResearchState& s = Research_Team(teamId);
        if (!n || s.Done(r.id)) continue;
        int& have = s.progress[r.id];
        have += r.units;
        if (have >= n->units) { // finished, as a lab finishes it
            s.done.insert(r.id);
            s.progress.erase(r.id);
            s.queue.erase(std::remove(s.queue.begin(), s.queue.end(), r.id), s.queue.end());
        }
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>

// (#2072) A boost adds units to the node (finishing it when they reach its units); an alternate recipe joins the team's list.
TEST(Caches, ApplyTeamAddsResearchUnitsAndAlternates) {
    ASSERT_TRUE(Research_Load());
    Research_Reset(); Caches_ResetTeams();
    const ResearchNode* n = nullptr;
    for (const ResearchNode& c : Research_Nodes()) if (c.units >= 8) { n = &c; break; }
    ASSERT_NE(n, nullptr);
    CacheRoll boost; boost.kind = CacheReward::Boost; boost.id = n->id; boost.units = n->units / 2;
    CacheRoll alt; alt.kind = CacheReward::Alt; alt.id = "packed_coke";
    Caches_ApplyTeam(5, { boost, alt });
    EXPECT_EQ(Research_Team(5).progress[n->id], n->units / 2);
    EXPECT_FALSE(Research_Team(5).Done(n->id));
    EXPECT_TRUE(Caches_TeamHasAlt(5, "packed_coke"));
    boost.units = n->units;
    Caches_ApplyTeam(5, { boost });
    EXPECT_TRUE(Research_Team(5).Done(n->id)) << "enough units finish the node";
    EXPECT_EQ(Research_Team(5).progress.count(n->id), 0u);
    Research_Reset(); Caches_ResetTeams();
}
#endif
