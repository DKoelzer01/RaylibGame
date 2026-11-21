#include "magicdiscovery.h"

#include "magicsight.h"
#include "manafield.h"

#include <raymath.h>
#include <algorithm>
#include <map>
#include <set>

// ---- The rules (pure) ----

MagicFind MagicDiscovery_OnFind(MagicPlayerState& s, bool carriesPage) {
    MagicFind f;
    if (s.found) return f;
    s.found = true;
    f.first = true;
    f.givePage = !s.sight && !carriesPage;
    return f;
}

bool MagicDiscovery_NearLine(const ManaField& f, Vector3 dir, float range) {
    dir = Vector3Normalize(dir);
    for (int k = 0; k < kManaBases; ++k)
        if (f.strength(k, dir) >= 0.3f && f.lineDistance(k, dir) <= range) return true;
    return false;
}

int MagicDiscovery_CrystalOdds(float flow) {
    if (flow < 0.3f) return 0;                            // only on a strong line's core
    return std::max(60, (int)(400.0f / std::min(flow, 2.0f))); // ~1 in 400 strokes on a full line, 1 in 200 at a node (guess)
}

int MagicDiscovery_NodeKindsStudied(const MagicPlayerState& s) {
    int n = 0;
    for (const auto& [type, count] : s.studied) n += count > 0 && type.rfind("node:", 0) == 0;
    return n;
}

// ---- The game ----
#include "magic.h"
#include "tome.h"
#include "researchtab.h"
#include "envflag.h"
#include <cstdio>
#include "gamerng.h"
#include "gameui.h"
#include "inventory.h"
#include "particles.h"
#include "playersession.h"
#include "survival.h"
#include "terraintool.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

namespace {
    constexpr const char* kPage = "first_sight_page";
    constexpr const char* kCrystal = "resonant_crystal";
    float moteClock = 0.0f;

    Planetoid* StartPlanet(Scene& world) { return world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get()); }
    Color KindColour(int kind) {
        const uint32_t c = Mana_Kind(kind >= 0 ? kind : Mana_Aether()).colour;
        return { (unsigned char)(c >> 16), (unsigned char)(c >> 8), (unsigned char)c, 255 };
    }

    void Find(MagicRoute route) {
        MagicPlayerState& s = LocalSession().magic;
        const MagicFind f = MagicDiscovery_OnFind(s, inventory.count(kPage) > 0);
        if (!f.first) return;
        QuestBook_Count("arcane_find"); // opens the Arcane line
        if (f.givePage && itemDB.find(kPage)) inventory.add(kPage, 1);
        switch (route) {
            case MagicRoute::Node: GameUI_Toast("The ground hums under your feet, and a faded page lifts out of the earth: First Sight."); break;
            case MagicRoute::Crystal: GameUI_Toast("The crystal rings in your hand. Words come with the sound, and you write them down: First Sight."); break;
            case MagicRoute::Loot: GameUI_Toast("The page hums in your hand: First Sight. Hold it and right-click to read it."); break;
            case MagicRoute::Creature: GameUI_Toast("The mote settles in your palm and leaves words behind: First Sight."); break;
        }
    }
}

bool MagicDiscovery_IsPage(const std::string& id) { return id == kPage && Magic_Enabled(); }

bool MagicDiscovery_ReadPage() {
    MagicPlayerState& s = LocalSession().magic;
    if (!Magic_Enabled() || inventory.count(kPage) <= 0) return false;
    if (s.sight) { GameUI_Toast("You know this page by heart: your sight is open to you (tap Y)."); return false; }
    MagicSight_Learn(s);
    inventory.remove(kPage, 1);
    s.found = true;
    QuestBook_Count("arcane_read");
    Tome_Team(ResearchTab_Team()).pages.insert("first_sight"); // (#1251) the team's first tome page (research's resonance tuning needs it)
    GameUI_Toast("First Sight. The world is threaded with faint currents: tap Y to see them, hold Y on something to study it.");
    return true;
}

bool MagicDiscovery_LootAllowed(const Planetoid& body, Vector3 local, const std::string& item) {
    if (item != kPage) return true;
    const ManaField* f = Magic_FieldFor(body);
    return f && Vector3Length(local) > 1.0f && MagicDiscovery_NearLine(*f, local, 40.0f); // within the first sight's 50 m
}

void MagicDiscovery_Install() {
    static bool done = false;
    if (done || !Magic_Enabled()) return;
    done = true;
    // A resonant crystal in rock dug on a strong ley line (the terrain tool asks after its gem roll).
    gRockFind = [](Planetoid& body, Vector3 hitLocal) -> const char* {
        const ManaField* f = Magic_FieldFor(body);
        if (!f || Vector3Length(hitLocal) < 1.0f) return nullptr;
        const Vector3 dir = Vector3Normalize(hitLocal);
        float flow = 0.0f;
        for (int k = 0; k < kManaBases; ++k) flow = std::max(flow, f->lineFlow(k, dir));
        const int odds = MagicDiscovery_CrystalOdds(flow);
        return odds > 0 && GameRng_Int(0, odds - 1) == 0 ? kCrystal : nullptr;
    };
}

void MagicDiscovery_Update(Scene& world, float dt) {
    if (!Magic_Enabled()) return;
    MagicDiscovery_Install();
    MagicPlayerState& s = LocalSession().magic;
    // Captures (RAYLIBGAME_MAGIC_SHOW_TOME=1, #1251): once the world is idle, a reader who has found magic, with First Sight
    // and a handful of aspects, and the game window open on the Tome.
    static bool showTome = EnvFlag("RAYLIBGAME_MAGIC_SHOW_TOME");
    if (showTome && WorldIdle()) {
        showTome = false;
        s.found = true;
        Tome_Team(ResearchTab_Team()).pages.insert("first_sight");
        for (const auto& [id, n] : std::map<std::string, float>{ { "air", 7 }, { "light", 5 }, { "order", 6 }, { "crystal", 3 }, { "earth", 4 } }) s.aspects[id] += n;
        GameUI_ShowTab("Tome");
    }
    // Finds by item: a page from a chest, a crystal from the rock.
    if (!s.found && inventory.count(kPage) > 0) Find(MagicRoute::Loot);
    if (!s.found && inventory.count(kCrystal) > 0) Find(MagicRoute::Crystal);
    const ManaField* f = Magic_Field(world);
    Planetoid* body = StartPlanet(world);
    const auto* pf = body ? dynamic_cast<const PlanetField*>(body->field.get()) : nullptr;
    if (!f || !pf) return; // (not player.nearest: a bench never sets it)
    const Vector3 local = body->toLocalOf(player.position, world.renderOrigin);
    const float r = Vector3Length(local);
    if (r < 1.0f || r > pf->radius + 2000.0f) return; // on or near the start planet
    const Vector3 dir = Vector3Scale(local, 1.0f / r);
    const float ground = pf->radius + pf->elevation(dir, 2.0f);
    // Captures (RAYLIBGAME_MAGIC_AT_NODE=1): once the world is idle, the player stands on the nearest node within 3 km, and
    // the camera watches it from 9 m off and 4 m up every frame (the subject stays in view).
    static int atNode = EnvFlag("RAYLIBGAME_MAGIC_AT_NODE") ? 0 : -1;
    static ManaNode watched;
    if (atNode == 0 && WorldIdle()) {
        std::vector<ManaNode> near;
        f->nodesNear(dir, 3000.0f, near);
        float best = 1e9f;
        for (const ManaNode& n : near) if (const float d = acosf(std::clamp(Vector3DotProduct(dir, n.dir), -1.0f, 1.0f)); d < best) { best = d; watched = n; }
        if (!near.empty()) {
            const float gn = pf->radius + std::max(pf->elevation(watched.dir, 2.0f), 0.0f);
            player.placeAt(world, *body, body->toWorld(Vector3Scale(watched.dir, gn + 1.0f)));
            printf("magic: at node: %s, %.0f m from the start\n", watched.kind >= 0 ? Mana_Kind(watched.kind).name.c_str() : "?", best * pf->radius);
        }
        atNode = near.empty() ? 2 : 1;
    }
    if (atNode == 1) {
        const float gn = pf->radius + std::max(pf->elevation(watched.dir, 2.0f), 0.0f);
        Vector3 e = Vector3CrossProduct({ 0.0f, 1.0f, 0.0f }, watched.dir);
        if (Vector3Length(e) < 1e-4f) e = { 1.0f, 0.0f, 0.0f };
        e = Vector3Normalize(e);
        const Vector3 at = Vector3Scale(watched.dir, gn + 1.0f);
        world.camera.position = body->toWorld(Vector3Add(at, Vector3Add(Vector3Scale(e, 9.0f), Vector3Scale(watched.dir, 4.0f))));
        world.camera.target = body->toWorld(at);
        world.camera.up = body->dirToWorld(watched.dir);
    }
    std::vector<ManaNode> nodes;
    f->nodesNear(dir, 40.0f, nodes);
    // Faint motes rise over a node, sight or no sight: the place itself hums (the doc's "humming fault spring, glowing pool").
    moteClock += dt;
    if (moteClock >= 0.3f) {
        moteClock = 0.0f;
        for (const ManaNode& n : nodes) {
            const float gn = pf->radius + std::max(pf->elevation(n.dir, 2.0f), 0.0f);
            Vector3 e = Vector3CrossProduct({ 0.0f, 1.0f, 0.0f }, n.dir);
            if (Vector3Length(e) < 1e-4f) e = { 1.0f, 0.0f, 0.0f };
            e = Vector3Normalize(e);
            const Vector3 nn = Vector3CrossProduct(n.dir, e);
            const float a = (float)GameRng_Int(0, 628) / 100.0f, d = (float)GameRng_Int(0, 300) / 100.0f;
            const Vector3 at = Vector3Add(Vector3Scale(n.dir, gn + 0.4f), Vector3Add(Vector3Scale(e, cosf(a) * d), Vector3Scale(nn, sinf(a) * d)));
            ParticleBurst b;
            b.kind = ParticleKind::Firefly; b.color = KindColour(n.kind);
            b.count = 3; b.size = 0.18f + 0.14f * n.strength; b.life = 5.0f; b.speed = 0.35f; b.spread = 0.3f; // (capture 2026-11-21: 0.1 m motes were lost in the grass)
            Particles_Emit(*body, at, b);
        }
    }
    // Standing on a node (within 4 m, on the ground) wakes it.
    for (const ManaNode& n : nodes) {
        const float along = acosf(std::clamp(Vector3DotProduct(dir, n.dir), -1.0f, 1.0f)) * pf->radius;
        if (along > 4.0f || r - ground > 4.0f) continue;
        if (!s.found) Find(MagicRoute::Node);
        if (s.sight) QuestBook_AtLeast("arcane_node", 1.0);
        break;
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1250): the first find, line nearness, crystal odds, node kinds studied ----
#include <gtest/gtest.h>

TEST(MagicDiscovery, OnlyTheFirstFindCountsAndThePageComesOnce) {
    MagicPlayerState s;
    MagicFind f = MagicDiscovery_OnFind(s, false);
    EXPECT_TRUE(f.first); EXPECT_TRUE(f.givePage); EXPECT_TRUE(s.found);
    f = MagicDiscovery_OnFind(s, false);
    EXPECT_FALSE(f.first) << "a second find is no news";
    EXPECT_FALSE(f.givePage);
    MagicPlayerState carrying;
    EXPECT_FALSE(MagicDiscovery_OnFind(carrying, true).givePage) << "a page from a chest is already the page";
    MagicPlayerState knows; knows.sight = true;
    EXPECT_FALSE(MagicDiscovery_OnFind(knows, false).givePage) << "no page for one who already sees";
}

TEST(MagicDiscovery, CrystalsOnlyOnStrongLines) {
    EXPECT_EQ(MagicDiscovery_CrystalOdds(0.0f), 0);
    EXPECT_EQ(MagicDiscovery_CrystalOdds(0.29f), 0);
    EXPECT_GT(MagicDiscovery_CrystalOdds(0.5f), MagicDiscovery_CrystalOdds(1.0f)) << "rarer on a weaker line";
    EXPECT_GE(MagicDiscovery_CrystalOdds(2.0f), 60) << "never common";
}

TEST(MagicDiscovery, NodeKindsStudiedCountsDistinctNodeKinds) {
    MagicPlayerState s;
    s.studied["node:storm"] = 3; s.studied["node:shard"] = 1; s.studied["creature:deer"] = 2; s.studied["rock:granite"] = 1;
    EXPECT_EQ(MagicDiscovery_NodeKindsStudied(s), 2);
}

TEST(MagicDiscovery, NearLineFindsLinesAndNotEmptyGround) {
    ASSERT_TRUE(Mana_Load());
    ManaFieldParams p; p.seed = 42; p.radius = 47000.0f;
    ManaField f(p);
    std::vector<ManaSegment> segs;
    f.lineSegments(Vector3Normalize({ 0.2f, 0.9f, 0.3f }), 4000.0f, 40.0f, 0.4f, segs); // (0.4: clear of the 0.3 bar)
    ASSERT_FALSE(segs.empty());
    const Vector3 on = Vector3Normalize(Vector3Add(segs[0].a, segs[0].b));
    EXPECT_TRUE(MagicDiscovery_NearLine(f, on, 40.0f)) << "on a line of strength >= 0.3";
    int far = 0, checked = 0;
    for (int k = 0; k < 300; ++k) {
        const Vector3 d = Vector3Normalize({ sinf(k * 0.77f), cosf(k * 0.31f), sinf(k * 0.13f + 1.0f) });
        bool anyClose = false;
        for (int b = 0; b < kManaBases; ++b) anyClose = anyClose || f.lineDistance(b, d) <= 40.0f;
        if (anyClose) continue;
        ++checked;
        far += !MagicDiscovery_NearLine(f, d, 40.0f);
    }
    EXPECT_GT(checked, 50) << "most of the ground is away from every line";
    EXPECT_EQ(far, checked) << "away from every line is never near one";
}
#endif
