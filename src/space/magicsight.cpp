#include "magicsight.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cctype>
#include <cstdlib>
#include <cstring>

// ---- The rules (pure) ----

float MagicSight_Range(const MagicPlayerState& s) {
    static const float kRange[3] = { 50.0f, 150.0f, 400.0f }; // the doc's attunement steps (guess)
    return kRange[std::clamp(s.sightTier, 0, 2)];
}

float MagicSight_Reserved(const MagicPlayerState& s) { return s.sight && s.shown ? kSightReserve * s.poolMax : 0.0f; }
float MagicSight_Available(const MagicPlayerState& s) { return std::max(0.0f, s.pool - MagicSight_Reserved(s)); }

bool MagicSight_Learn(MagicPlayerState& s) {
    if (s.sight) return false;
    s.sight = true;
    return true;
}

bool MagicSight_SetShown(MagicPlayerState& s, bool on) {
    s.shown = s.sight && on;
    return s.shown == on;
}

float MagicSight_StudyValue(const MagicPlayerState& s, const std::string& type) {
    const auto it = s.studied.find(type);
    const int n = it == s.studied.end() ? 0 : it->second;
    return std::max(0.1f, powf(0.5f, (float)n));
}

float MagicSight_Study(MagicPlayerState& s, const std::string& type) {
    const float v = MagicSight_StudyValue(s, type);
    ++s.studied[type];
    return v;
}

void MagicSight_Regen(MagicPlayerState& s, float dt, float flow) {
    const float rate = 0.5f * (1.0f + 3.0f * std::clamp(flow, 0.0f, 2.0f)); // pool a second (guess): 0.5 off the lines, up to 3.5 on a node
    s.pool = std::min(s.poolMax, s.pool + rate * std::max(0.0f, dt));
}

std::string MagicSight_SaveText(const MagicPlayerState& s) {
    std::string out;
    char buf[200];
    if (s.sight || s.sightTier || s.throughRock) {
        snprintf(buf, sizeof buf, "magic_sight %d %d %d %d\n", s.sight ? 1 : 0, s.shown ? 1 : 0, s.sightTier, s.throughRock ? 1 : 0);
        out += buf;
    }
    snprintf(buf, sizeof buf, "magic_pool %.3f %.3f\n", s.pool, s.poolMax);
    out += buf;
    if (s.found) out += "magic_found 1\n"; // (#1250)
    if (s.strain > 0.0f || s.burden > 0.0f || s.wardingCharm) { // (#1255)
        snprintf(buf, sizeof buf, "magic_burden %.3f %.3f %d\n", s.strain, s.burden, s.wardingCharm ? 1 : 0);
        out += buf;
    }
    for (const auto& [id, v] : s.aspects) { // (#1251)
        if (v <= 0.0f || id.empty() || id.find_first_of(" \t\r\n") != std::string::npos) continue;
        snprintf(buf, sizeof buf, "magic_aspect %s %.3f\n", id.c_str(), v);
        out += buf;
    }
    for (const auto& [type, n] : s.studied) {
        if (n <= 0 || type.empty() || type.find_first_of(" \t\r\n") != std::string::npos) continue;
        snprintf(buf, sizeof buf, "magic_studied %s %d\n", type.c_str(), n);
        out += buf;
    }
    return out;
}

bool MagicSight_LoadLine(MagicPlayerState& s, const char* line) {
    int a = 0, b = 0, tier = 0, rock = 0, n = 0;
    float pool = 0.0f, poolMax = 0.0f;
    char type[160] = {};
    if (sscanf(line, "magic_sight %d %d %d %d", &a, &b, &tier, &rock) == 4) {
        s.sight = a != 0; s.shown = s.sight && b != 0; s.sightTier = std::clamp(tier, 0, 2); s.throughRock = rock != 0;
        return true;
    }
    if (sscanf(line, "magic_pool %f %f", &pool, &poolMax) == 2) {
        s.poolMax = std::max(1.0f, poolMax);
        s.pool = std::clamp(pool, 0.0f, s.poolMax);
        return true;
    }
    if (sscanf(line, "magic_studied %159s %d", type, &n) == 2) { if (n > 0) s.studied[type] = n; return true; }
    if (sscanf(line, "magic_found %d", &n) == 1) { s.found = n != 0; return true; } // (#1250)
    if (float st = 0.0f, bu = 0.0f; sscanf(line, "magic_burden %f %f %d", &st, &bu, &n) == 3) { // (#1255)
        s.strain = std::clamp(st, 0.0f, 100.0f); s.burden = std::clamp(bu, 0.0f, 100.0f); s.wardingCharm = n != 0;
        return true;
    }
    if (sscanf(line, "magic_aspect %159s %f", type, &pool) == 2) { if (pool > 0.0f) s.aspects[type] = pool; return true; } // (#1251)
    return false;
}

// ---- The game: the local player's key, studies, pool, look and HUD ----

#include "magic.h"
#include "animals.h"
#include "envflag.h"
#include "gameui.h"
#include "playerinput.h"
#include "playersession.h"
#include "terraintool.h"
#include "voxel.h"
#include "world.h"
#include "magicdiscovery.h" // (#1250) the Arcane line's study steps
#include "survival.h"       // QuestBook_Count
#include "tome.h"           // (#1251) aspects from studies, the team's pages
#include "burden.h"         // (#1255) strain and burden
#include "playerrecord.h"   // Teams_Of, PlayerId_Local: the local team
#include "worldtypes.h"

#include <raylib.h>
#include <raymath.h>

namespace {
    float held = 0.0f;         // s the key has been down
    bool wasDown = false;
    std::string target, targetName; // what is being studied
    float progress = 0.0f;
    bool studiedThisHold = false;

    Planetoid* StartPlanet(Scene& world) { return world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get()); }

    // What the crosshair is on, within the sight's range: an animal, a node's glow, or rock.
    bool StudyTarget(Scene& world, const MagicPlayerState& s, std::string& type, std::string& name) {
        const float range = MagicSight_Range(s);
        const int a = Animals_AimedIndex();
        if (a >= 0) {
            const auto [id, display] = Animals_SpeciesOf(a);
            if (!id.empty()) { type = "creature:" + id; name = display.empty() ? id : display; return true; }
        }
        ManaNode node; float d = 0.0f;
        if (Magic_SightNodeAhead(world, range, node, d) && node.kind >= 0) {
            type = "node:" + Mana_Kind(node.kind).id;
            name = Mana_Kind(node.kind).name + " node";
            return true;
        }
        Planetoid* body = nullptr; Vector3 local{}, normal{};
        if (TerrainTool_LastHit(body, local, normal) && body && body == StartPlanet(world)) {
            const auto* pf = dynamic_cast<const PlanetField*>(body->field.get());
            const Vector3 eye = body->toLocal(world.camera.position);
            const float r = Vector3Length(local);
            if (pf && r > 1.0f && Vector3Distance(eye, local) <= range) {
                const char* rock = pf->stoneKind(Vector3Scale(local, (r - 0.5f) / r)); // just inside the face
                if (rock && rock[0]) {
                    type = std::string("rock:") + rock;
                    name = rock;
                    if (!name.empty()) name[0] = (char)toupper((unsigned char)name[0]);
                    return true;
                }
            }
        }
        return false;
    }
}

void MagicSight_Update(Scene& world, float dt, const PlayerInput& in) {
    if (!Magic_Enabled()) return;
    MagicPlayerState& s = LocalSession().magic;
    static const bool grant = EnvFlag("RAYLIBGAME_MAGIC_SIGHT"); // until discovery (#1250) teaches First Sight
    if (const int team = Teams_Of(PlayerId_Local())) Tome_ApplyEffects(Tome_Team(team), s); // (#1251) the team's pages: First Sight, range, through rock (no team yet: none)
    if (grant && MagicSight_Learn(s)) { // =show: and shown (captures)
        const char* e = getenv("RAYLIBGAME_MAGIC_SIGHT");
        if (e && strcmp(e, "show") == 0) MagicSight_SetShown(s, true);
        if (const char* t = getenv("RAYLIBGAME_MAGIC_SIGHT_TIER")) s.sightTier = std::clamp(atoi(t), 0, 2); // captures: 0-2 (50/150/400 m) before the attunement pages
    }
    // The pool refills, faster on a ley line.
    if (const ManaField* f = Magic_Field(world)) {
        float flow = 0.0f;
        if (Planetoid* p = StartPlanet(world); p && player.nearest == p) {
            const Vector3 local = p->toLocalOf(player.position, world.renderOrigin);
            if (Vector3Length(local) > 1.0f) for (int k = 0; k < kManaBases; ++k) flow = std::max(flow, f->lineFlow(k, Vector3Normalize(local)));
        }
        MagicSight_Regen(s, dt * Burden_RegenFactor(s.burden), flow); // (#1255) Unsettled and up: 10% slower
        Strain_Drain(s, dt, false);
    }
    // The key: a tap toggles the overlay, a hold studies what the crosshair is on (only while the sight is shown).
    const bool down = in.down(Act::MagicSight);
    if (down) {
        held += dt;
        if (held >= kSightTapSeconds && s.shown && !studiedThisHold) {
            std::string type, name;
            if (StudyTarget(world, s, type, name)) {
                if (type != target) { target = type; targetName = name; progress = 0.0f; }
                progress += dt;
                if (progress >= kStudySeconds) {
                    const float v = MagicSight_Study(s, target);
                    QuestBook_Count("arcane_study"); // (#1250) the Arcane line
                    QuestBook_AtLeast("arcane_kinds", MagicDiscovery_NodeKindsStudied(s));
                    Tome_GrantStudy(s, target, v); // (#1251) its aspects
                    std::string gained;
                    if (const auto* asp = Aspects_Of(target))
                        for (const auto& [id, n] : *asp) { const AspectDef* a = Aspects_Find(id); gained += TextFormat("%s+%.1f %s", gained.empty() ? "" : ", ", v * n, a ? a->name.c_str() : id.c_str()); }
                    if (gained.empty()) gained = "nothing you can name yet";
                    GameUI_Toast(TextFormat(v >= 1.0f ? "Studied: %s (%s)" : "Studied: %s again, less to learn (%s)", targetName.c_str(), gained.c_str()));
                    studiedThisHold = true;
                    progress = 0.0f;
                }
            } else { target.clear(); progress = 0.0f; }
        }
    } else if (wasDown) {
        if (held < kSightTapSeconds) {
            if (!s.sight) GameUI_Toast("Nothing answers. (Magic sight is learned from the First Sight page.)");
            else if (MagicSight_SetShown(s, !s.shown) && s.shown) QuestBook_Count("arcane_sight"); // (#1250)
        }
        held = 0.0f; target.clear(); progress = 0.0f; studiedThisHold = false;
    }
    wasDown = down;
}

void MagicSight_Draw3D(Scene& world) {
    if (!Magic_Enabled()) return;
    const MagicPlayerState& s = LocalSession().magic;
    if (s.sight && s.shown) Magic_DrawSight3D(world, MagicSight_Range(s), s.throughRock);
}

void MagicSight_DrawHUD(Scene& world) {
    (void)world;
    if (!Magic_Enabled()) return;
    const MagicPlayerState& s = LocalSession().magic;
    if (!s.sight) return;
    const float sc = GameUI_HudScale();
    const int W = GetScreenWidth(), H = GetScreenHeight();
    // The pool, with the part the sight holds back, above the hotbar.
    const int bw = (int)(220 * sc), bh = (int)(8 * sc), bx = W / 2 - bw / 2, by = H - (int)(132 * sc);
    DrawRectangle(bx, by, bw, bh, Fade(BLACK, 0.45f));
    DrawRectangle(bx, by, (int)(bw * std::clamp(s.pool / s.poolMax, 0.0f, 1.0f)), bh, Color{ 200, 170, 255, 200 });
    const float reserved = MagicSight_Reserved(s) / s.poolMax;
    if (reserved > 0.0f) DrawRectangle(bx + (int)(bw * (1.0f - reserved)), by, (int)(bw * reserved), bh, Color{ 120, 90, 200, 220 });
    const char* label = TextFormat("Aether %.0f / %.0f%s%s", s.pool, s.poolMax, s.shown ? TextFormat("  (sight holds %.0f)", MagicSight_Reserved(s)) : "", s.burden >= 20.0f ? TextFormat("   %s", Burden_BandName(Burden_Band(s.burden))) : ""); // (#1255) the burden band beside the pool
    DrawText(label, bx, by - (int)(16 * sc), (int)(14 * sc), Fade(RAYWHITE, 0.85f));
    // Studying: a bar under the crosshair with what it is.
    if (!target.empty() && progress > 0.0f) {
        const int sw = (int)(140 * sc), sx = W / 2 - sw / 2, sy = H / 2 + (int)(36 * sc);
        DrawRectangle(sx, sy, sw, (int)(5 * sc), Fade(BLACK, 0.5f));
        DrawRectangle(sx, sy, (int)(sw * std::clamp(progress / kStudySeconds, 0.0f, 1.0f)), (int)(5 * sc), Color{ 200, 170, 255, 230 });
        DrawText(TextFormat("Studying %s", targetName.c_str()), sx, sy + (int)(8 * sc), (int)(14 * sc), Fade(RAYWHITE, 0.9f));
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1249): range, the pool's reserve, learning and showing, studies, the record's lines ----
#include <gtest/gtest.h>
#include <sstream>

TEST(MagicSight, ShownOnlyWhenKnownAndHoldsTenPercent) {
    MagicPlayerState s;
    EXPECT_FALSE(MagicSight_SetShown(s, true)) << "not known yet";
    EXPECT_FALSE(s.shown);
    EXPECT_FLOAT_EQ(MagicSight_Reserved(s), 0.0f);
    EXPECT_TRUE(MagicSight_Learn(s));
    EXPECT_FALSE(MagicSight_Learn(s)) << "permanent: learned once";
    EXPECT_FLOAT_EQ(MagicSight_Reserved(s), 0.0f) << "free while hidden (decided)";
    EXPECT_TRUE(MagicSight_SetShown(s, true));
    EXPECT_FLOAT_EQ(MagicSight_Reserved(s), 10.0f) << "10% of the pool while shown";
    EXPECT_FLOAT_EQ(MagicSight_Available(s), 90.0f);
    s.pool = 5.0f;
    EXPECT_FLOAT_EQ(MagicSight_Available(s), 0.0f) << "never below zero";
    s.poolMax = 250.0f; s.pool = 250.0f;
    EXPECT_FLOAT_EQ(MagicSight_Reserved(s), 25.0f) << "a share of the pool, not a fixed amount";
}

TEST(MagicSight, RangeFollowsAttunement) {
    MagicPlayerState s;
    EXPECT_FLOAT_EQ(MagicSight_Range(s), 50.0f);
    s.sightTier = 1; EXPECT_FLOAT_EQ(MagicSight_Range(s), 150.0f);
    s.sightTier = 2; EXPECT_FLOAT_EQ(MagicSight_Range(s), 400.0f);
    s.sightTier = 9; EXPECT_FLOAT_EQ(MagicSight_Range(s), 400.0f);
}

TEST(MagicSight, StudyingGivesFullInsightTheFirstTimeAndLessAfter) {
    MagicPlayerState s;
    EXPECT_FLOAT_EQ(MagicSight_Study(s, "node:storm"), 1.0f);
    EXPECT_FLOAT_EQ(MagicSight_Study(s, "node:storm"), 0.5f);
    EXPECT_FLOAT_EQ(MagicSight_Study(s, "node:storm"), 0.25f);
    EXPECT_FLOAT_EQ(MagicSight_Study(s, "creature:deer"), 1.0f) << "per type";
    for (int i = 0; i < 10; ++i) MagicSight_Study(s, "node:storm");
    EXPECT_FLOAT_EQ(MagicSight_StudyValue(s, "node:storm"), 0.1f) << "never nothing";
}

TEST(MagicSight, ThePoolRefillsFasterOnLeyLines) {
    MagicPlayerState off, on;
    off.pool = on.pool = 0.0f;
    MagicSight_Regen(off, 10.0f, 0.0f);
    MagicSight_Regen(on, 10.0f, 1.0f);
    EXPECT_GT(off.pool, 0.0f);
    EXPECT_GT(on.pool, 3.0f * off.pool);
    MagicSight_Regen(on, 1e6f, 2.0f);
    EXPECT_FLOAT_EQ(on.pool, on.poolMax) << "capped";
}

TEST(MagicSight, TheRecordKeepsItThroughSaveAndLoad) {
    MagicPlayerState a;
    MagicSight_Learn(a);
    MagicSight_SetShown(a, true);
    a.sightTier = 1; a.throughRock = true; a.pool = 42.5f; a.poolMax = 120.0f;
    a.studied["node:storm"] = 3; a.studied["rock:granite"] = 1;
    MagicPlayerState b;
    std::istringstream in(MagicSight_SaveText(a));
    for (std::string line; std::getline(in, line);) EXPECT_TRUE(MagicSight_LoadLine(b, line.c_str())) << line;
    EXPECT_TRUE(b.sight); EXPECT_TRUE(b.shown); EXPECT_EQ(b.sightTier, 1); EXPECT_TRUE(b.throughRock);
    EXPECT_NEAR(b.pool, 42.5f, 1e-3f); EXPECT_NEAR(b.poolMax, 120.0f, 1e-3f);
    EXPECT_EQ(b.studied, a.studied);
    EXPECT_FALSE(MagicSight_LoadLine(b, "skill_xp mining 3"));
    MagicPlayerState c;
    EXPECT_TRUE(MagicSight_LoadLine(c, "magic_sight 0 1 0 0"));
    EXPECT_FALSE(c.shown) << "never shown without the sight";
}
#endif
