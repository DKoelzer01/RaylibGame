// terrain/place-dig-back (#1097): placed ground digs back at most what it cost, so placing and digging can't make items.
// On the sandbox, for each filler (gravel, sand: not what the sandbox ground is made of, so only the placed ground counts)
// and each digging tool (bare hands; a stone shovel with Digging's perks, which must add nothing on placed ground): three
// placing strokes at the aimed ground, then six digging strokes at the same spot; the filler's net change must be <= 0.
#include "scenario.h"
#include "playerinput.h"
#include "playersession.h"
#include "skills.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
// The camera at the player's eyes, looking at the ground 2 m ahead; then the tool's aim (its hit) from it.
bool AimAhead(ScenarioCtx& c) {
    Player& p = LocalSession().player;
    const Vector3 feet = (p.position - c.world.renderOrigin).toFloat();
    const Vector3 up = Vector3Normalize(p.up);
    Vector3 fwd = Vector3Subtract(p.forward, Vector3Scale(up, Vector3DotProduct(p.forward, up)));
    fwd = Vector3Length(fwd) < 1e-3f ? Vector3Normalize(Vector3CrossProduct(up, { 0, 0, 1 })) : Vector3Normalize(fwd);
    c.world.camera.position = Vector3Add(feet, Vector3Scale(up, 1.6f));
    c.world.camera.target = Vector3Add(feet, Vector3Scale(fwd, 2.0f));
    c.world.camera.up = up;
    TerrainTool_Update(c.world, PlayerInput::None(), true);
    Planetoid* body = nullptr; Vector3 local{}, n{};
    return TerrainTool_LastHit(body, local, n);
}
}

SCENARIO("terrain/place-dig-back") {
    s.tags = { "world" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.step("place three strokes, dig the spot back: no filler gained", 40.0, [](ScenarioCtx& c) {
        c.waitingFor("the player standing on the sandbox");
        if (c.T < 1.0) return false;
        terrainTool.creative = false;
        SkillState& sk = LocalSession().skills;
        const SkillState saved = sk;
        struct Case { const char* filler; const char* tool; bool perks; };
        const Case cases[] = { { "gravel", "", false }, { "gravel", "stone_shovel", true }, { "sand", "", false }, { "sand", "stone_shovel", true } };
        for (const Case& k : cases) {
            sk.ranks.clear();
            if (k.perks) sk.ranks = { { "quick_spade", 1 }, { "wide_scoop", 2 } }; // (Digging: wider and faster soil strokes)
            inventory.add(k.filler, 20);
            const double before = inventory.amount(k.filler);
            terrainTool.heldTool = k.filler; // placing: the filler in hand
            int placed = 0;
            for (int n = 0; n < 3; ++n) if (AimAhead(c) && TerrainTool_Stroke(false)) ++placed;
            const double afterPlace = inventory.amount(k.filler);
            terrainTool.heldTool = k.tool;   // digging: the tool in hand
            for (int n = 0; n < 6; ++n) if (AimAhead(c)) TerrainTool_Stroke(true);
            const double after = inventory.amount(k.filler);
            const std::string what = std::string(k.filler) + " with " + (k.tool[0] ? k.tool : "bare hands") + (k.perks ? " and Digging perks" : "");
            c.eq(("strokes placed: " + what).c_str(), placed, 3);
            c.gt(("placing cost some " + what).c_str(), before - afterPlace, 0.1);
            c.gt(("dug back no more than it cost: " + what + " (net change <= 0)").c_str(), before + 1e-3, after);
            c.note(TextFormat("%s: %.2f -> placed %.2f -> dug back %.2f (net %+.2f)", what.c_str(), before, afterPlace, after, after - before));
        }
        sk = saved;
        return true;
    });
}
