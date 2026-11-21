// visual/lamps and visual/torch (#870 visual tier; moved from --test-lamps #139 and --test-torch #218): local lights light
// the ground. Seed 42 starts at night, so the player stands at the start looking down at the ground ahead; the ground's
// luminance (the lower middle of the screen) is measured before and after the lights, and each lit view is a shot against
// its baseline. PASS as the flags: lamps: >= 90 lamp cubes placed and >= 90 local lights uploaded, the ground brighter by
// > 6; torch: a torch in hand lights it (the hand light on, brighter by > 4), and one set down does too (a local light,
// brighter by > 4, one torch fewer).
#include "scenario.h"
#include "hotbar.h"
#include "inventory.h"
#include "locallights.h"
#include "particles.h"
#include "player.h"
#include "playersession.h"
#include "plots.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>

namespace {
// The ground's mean luminance (0..255) in a box of the screen (fractions of its size), every 2nd pixel.
double GroundLum(float x0, float x1, float y0, float y1) {
    Image im = LoadImageFromScreen();
    ImageFormat(&im, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    const Color* px = (const Color*)im.data;
    double sum = 0.0; int n = 0;
    for (int y = (int)(im.height * y0); y < (int)(im.height * y1); y += 2)
        for (int x = (int)(im.width * x0); x < (int)(im.width * x1); x += 2) { const Color c = px[y * im.width + x]; sum += 0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b; ++n; }
    UnloadImage(im);
    return n ? sum / n : 0.0;
}
// Stand at the start on foot (not flying), as the flags did.
void StandAtStart(ScenarioCtx& c) {
    player.flyMode = false;
    player.spawnOn(c.world, *c.body, c.spawnDirWorld());
}
// Look down at the ground in front (0.6 down, 0.8 ahead), each frame.
void LookDown() {
    const Vector3 ahead = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
    player.forward = Vector3Normalize(Vector3Subtract(Vector3Scale(ahead, 0.8f), Vector3Scale(player.up, 0.6f)));
}
}

SCENARIO("visual/lamps") {
    s.tags = { "visual", "lights", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 120.0;
    s.legacyFlag = "--test-lamps";
    struct St { double lumOff = 0.0; int placed = 0; };
    static St st;
    s.step("stand at the start", 30.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); if (c.T < 1.0 || !WorldIdle()) return false; StandAtStart(c); return true; });
    s.step("look down, unlit", 10.0, [](ScenarioCtx& c) { LookDown(); return c.t > 0.5; });
    s.step("measure unlit, place 99 lamp cubes", 5.0, [](ScenarioCtx& c) {
        LookDown();
        st.lumOff = GroundLum(0.2f, 0.8f, 0.5f, 0.9f);
        c.metric("ground luminance unlit", st.lumOff, "/255", false);
        Planetoid* p = c.body;
        const Vector3 me = p->toLocal((player.position - c.world.renderOrigin).toFloat());
        const Vector3 up = Vector3Normalize(me);
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        const Vector3 t2 = Vector3CrossProduct(t1, up);
        if (!p->plots) p->plots = std::make_unique<PlotSet>();
        const int lampMat = BuildMaterialForItem("lamp");
        Plot* pl = p->plots->plotFor(me, 0.0f, true);
        for (int i = -5; i < 5; ++i)
            for (int j = -5; j < 5; ++j) {
                if (abs(2 * i + 1) <= 1 && abs(2 * j + 1) <= 1) continue; // not on the player
                const Vector3 d = Vector3Normalize(Vector3Add(me, Vector3Add(Vector3Scale(t1, 2.0f * i + 1.0f), Vector3Scale(t2, 2.0f * j + 1.0f))));
                const Vector3 at = Vector3Scale(d, c.field->radius + c.field->elevation(d) + 0.3f);
                st.placed += pl->place(Part{ PartType::Cube, pl->cellAt(pl->toPlot(at)), 0, (uint8_t)std::max(0, lampMat) }) >= 0 ? 1 : 0;
            }
        return true;
    });
    s.step("lamps lit", 10.0, [](ScenarioCtx& c) { LookDown(); c.waitingFor("the lamps' lights"); return c.t > 0.6 && LocalLights_Count() >= 90; });
    s.step("measure lit", 5.0, [](ScenarioCtx& c) {
        LookDown();
        const double on = GroundLum(0.2f, 0.8f, 0.5f, 0.9f);
        c.metric("ground luminance lit", on, "/255", false);
        c.gt("lamp cubes placed", st.placed, 89.5);
        c.gt("local lights uploaded", LocalLights_Count(), 89.5);
        c.metric("lights in the busiest tile", LocalLights_MaxPerTile());
        c.gt("ground brighter by", on - st.lumOff, 6.0);
        c.shot("lamps");
        return true;
    });
}

SCENARIO("visual/torch") {
    s.tags = { "visual", "lights", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 120.0;
    s.legacyFlag = "--test-torch";
    struct St { double lumOff = 0.0; int before = 0; bool placed = false; };
    static St st;
    s.step("stand at the start with torches", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        StandAtStart(c);
        inventory.add("torch", 3); Hotbar_Bind(7, "torch"); Hotbar_Press(8, false); // (slot 8: an empty hand)
        return true;
    });
    s.step("look down, no light", 10.0, [](ScenarioCtx& c) { LookDown(); return c.t > 0.3; });
    s.step("measure, take the torch", 5.0, [](ScenarioCtx& c) {
        LookDown();
        st.lumOff = GroundLum(0.3f, 0.7f, 0.6f, 0.9f);
        c.metric("ground luminance unlit", st.lumOff, "/255", false);
        Hotbar_Press(7, false);
        return true;
    });
    s.step("torch in hand", 10.0, [](ScenarioCtx& c) { LookDown(); return c.t > 0.4; });
    s.step("measure in hand, set it down", 5.0, [](ScenarioCtx& c) {
        LookDown();
        const double on = GroundLum(0.3f, 0.7f, 0.6f, 0.9f);
        const bool lit = c.world.lights.size() > 3 && c.world.lights[3].enabled;
        c.check("the hand light is on", lit);
        c.gt("ground brighter by (in hand)", on - st.lumOff, 4.0);
        c.shot("torch_in_hand");
        st.before = inventory.count("torch");
        st.placed = Hotbar_PlaceTorch() && inventory.count("torch") == st.before - 1;
        Hotbar_Press(8, false); // the hand away
        return true;
    });
    s.step("placed torch burning", 10.0, [](ScenarioCtx& c) { LookDown(); c.waitingFor("its local light"); return c.t > 0.5 && LocalLights_Count() >= 1; });
    s.step("measure the placed torch", 5.0, [](ScenarioCtx& c) {
        LookDown();
        const double on = GroundLum(0.3f, 0.7f, 0.6f, 0.9f);
        const bool handOff = !(c.world.lights.size() > 3 && c.world.lights[3].enabled);
        c.check("one torch set down (one fewer carried)", st.placed);
        c.check("the hand light is off", handOff);
        c.gt("local lights", LocalLights_Count(), 0.5);
        c.gt("ground brighter by (placed)", on - st.lumOff, 4.0);
        c.metric("particles alive (its flame)", Particles_Count());
        c.shot("torch_placed");
        return true;
    });
}
