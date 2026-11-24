// ui/map-live (#1795): the maps show the live world. A plot of mixed parts (a stone floor, a furnace, a chest, a generator
// and a row of conveyors) is built a few metres from the player, turned 30 degrees off north (a turned plot, as #1771 will
// place them). Within the refresh window the minimap and the big map (opened over the plot at 1/4 m a pixel) show it: each
// part in its kind's colour at its place; the conveyor row along the plot's own east, not along the body's east (the
// heading); and the ground, not the plot, just past its edge. Then a close-up at 1/16 m a pixel, and a 2.5 m pit dug beside
// the plot: within the refresh window its edit tile is re-sampled (its depth) and the map draws it as dug earth (part 2).
// Then 12 m^3 poured into it: the simulation's pond is sampled and drawn as water (part 3).
// The recording shows the HUD: the minimap, then the big map.
#include "scenario.h"
#include "fluid.h"
#include "inventory.h"
#include "mapfeatures.h"
#include "navigation.h"
#include "player.h"
#include "plots.h"
#include "profiler.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>

namespace {
Plot* plot = nullptr;
Vector3 e0{}, n0{}; // the body's east and north at the plot (unturned)
double tBuilt = 0.0, tDug = 0.0;
bool opened = false;
Vector3 pitDir{}, pitCentre{};
Color pitBefore{};
double tFlood = -1.0, nextPour = 0.0;
float poured = 0.0f;

// A body-local direction at plot cell (u, v) (continuous cells), on the floor's top.
Vector3 AtCell(float u, float v) { return Vector3Normalize(plot->toBody({ u * kCellSize, kCellSize, v * kCellSize })); }
bool Same(Color a, Color b) { return a.r == b.r && a.g == b.g && a.b == b.b; }
// The kind whose colour a map pixel has (None: no plot colour, the ground).
MapKind KindOf(Color c) {
    for (int k = 1; k < (int)MapKind::Count; ++k) if (Same(c, MapFeatures_Colour((MapKind)k))) return (MapKind)k;
    return MapKind::None;
}
MapKind MapKindAt(bool big, Vector3 dir) {
    Color c;
    return Nav_MapPixel(big, dir, c) ? KindOf(c) : MapKind::Count;
}
std::string Name(MapKind k) { return k == MapKind::Count ? "off the image" : MapKindName(k); }
}

SCENARIO("ui/map-live") {
    s.tags = { "world", "ui" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day", "--record-hud", "1" }; // the maps are HUD: the recording shows them
    s.timeout = 150.0;
    s.step("settle on the ground", 40.0, [](ScenarioCtx& c) {
        player.flyMode = false;
        c.waitingFor("grounded and an idle world");
        return c.T >= 1.0 && player.grounded && WorldIdle();
    });
    s.step("a turned plot of mixed parts", 2.0, [](ScenarioCtx& c) {
        inventory.add("map", 1);
        const Vector3 feet = c.body->toLocalD(player.position).toFloat();
        const Vector3 up = Vector3Normalize(feet);
        e0 = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, up));
        n0 = Vector3CrossProduct(up, e0);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(Vector3Add(feet, Vector3Add(Vector3Scale(n0, 6.0f), Vector3Scale(e0, 3.0f))), 0.0f, true);
        if (!plot) { c.fail("no plot"); return true; }
        const float a = 30.0f * DEG2RAD; // turned before any part is placed: the footprint uses the plot's own axes
        plot->east = Vector3Normalize(Vector3Add(Vector3Scale(e0, cosf(a)), Vector3Scale(n0, sinf(a))));
        plot->north = Vector3CrossProduct(plot->up, plot->east);
        int placed = 0, want = 0;
        for (int x = 0; x < 16; ++x)
            for (int z = 0; z < 8; ++z) { ++want; placed += plot->place({ PartType::Cube, { x, 0, z }, 0, 0 }) >= 0; }
        for (const PartType t : { PartType::Furnace, PartType::Chest, PartType::Generator }) {
            ++want;
            placed += plot->place({ t, { t == PartType::Furnace ? 2 : t == PartType::Chest ? 6 : 10, 1, 2 }, 0, 0 }) >= 0;
        }
        for (int x = 0; x < 16; ++x) { ++want; placed += plot->place({ PartType::Conveyor, { x, 1, 6 }, 0, 0 }) >= 0; }
        c.eq("parts placed", placed, want);
        tBuilt = c.T;
        opened = false;
        return true;
    });
    s.step("the minimap shows the plot", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("the furnace on the minimap");
        if (MapKindAt(false, AtCell(3.0f, 3.0f)) != MapKind::Machine) return false;
        c.metric("minimap shows the plot after", c.T - tBuilt, "s");
        c.lt("minimap refresh after the build (s)", c.T - tBuilt, 10.0);
        c.eq("minimap: the conveyor row is transport", (double)MapKindAt(false, AtCell(13.5f, 6.5f)), (double)MapKind::Transport);
        return true;
    });
    s.step("the big map shows it at its heading", 20.0, [](ScenarioCtx& c) {
        if (!opened) { Nav_OpenMapAt(AtCell(8.0f, 4.0f), 0.25f); opened = true; }
        c.waitingFor("the big map at 1/4 m a pixel");
        if (!Nav_MapImageReady() || fabsf(Nav_MapMetresPerPixel(true) - 0.25f) > 1e-3f) return false;
        if (MapKindAt(true, AtCell(3.0f, 3.0f)) != MapKind::Machine) return false;
        c.metric("big map shows the plot after", c.T - tBuilt, "s");
        struct Probe { const char* what; float u, v; MapKind want; };
        const Probe probes[] = {
            { "furnace", 3.0f, 3.0f, MapKind::Machine },  { "chest", 7.0f, 3.0f, MapKind::Storage },
            { "generator", 11.0f, 3.0f, MapKind::Power }, { "conveyor", 8.5f, 6.5f, MapKind::Transport },
            { "bare floor", 14.5f, 1.5f, MapKind::Structure }, { "past the west edge", -4.0f, 4.0f, MapKind::None },
        };
        for (const Probe& p : probes) {
            const MapKind k = MapKindAt(true, AtCell(p.u, p.v));
            c.check(std::string("big map: ") + p.what + " is " + MapKindName(p.want), k == p.want, "got " + Name(k));
        }
        // The heading: the far corner (7.75 m east, 3.75 m north of the first cell's corner along the plot's own axes) is its
        // floor; the same offsets along the body's unturned axes land 0.6 m past the plot's east end and south edge.
        const Vector3 o = plot->toBody({ 0.0f, kCellSize, 0.0f });
        const MapKind turned = MapKindAt(true, Vector3Normalize(Vector3Add(o, Vector3Add(Vector3Scale(plot->east, 7.75f), Vector3Scale(plot->north, 3.75f)))));
        const MapKind cardinal = MapKindAt(true, Vector3Normalize(Vector3Add(o, Vector3Add(Vector3Scale(e0, 7.75f), Vector3Scale(n0, 3.75f)))));
        c.check("big map: the far corner along the plot's axes is its floor", turned == MapKind::Structure, "got " + Name(turned));
        c.check("big map: the same offsets along the body's axes are off the plot", cardinal == MapKind::None, "got " + Name(cardinal));
        c.grab("big-map");
        c.keepCapture("the plot on both maps (#1795)");
        Nav_OpenMapAt(AtCell(8.0f, 4.0f), 0.0625f); // the next step's close-up
        return true;
    });
    s.step("zoomed in to 1/16 m a pixel: outlines and machine icons", 15.0, [](ScenarioCtx& c) {
        c.waitingFor("the big map at 1/16 m a pixel");
        if (!Nav_MapImageReady() || fabsf(Nav_MapMetresPerPixel(true) - 0.0625f) > 1e-4f) return false;
        const MapKind k = MapKindAt(true, AtCell(3.0f, 3.0f));
        c.check("close up: the furnace is machine", k == MapKind::Machine, "got " + Name(k));
        c.grab("big-map-close");
        // (part 2) A pit 6 m west of the plot: 2.5 m deep, carved now; the map's colour there before it.
        pitDir = AtCell(-12.0f, 4.0f);
        Nav_MapPixel(true, pitDir, pitBefore);
        TerrainEdit e{};
        e.center = Vector3Scale(pitDir, c.field->radius + c.field->elevation(pitDir, 0.5f));
        e.radius = 2.5f;
        e.carve = true;
        ApplyTerrainEdit(*c.body, e);
        pitCentre = e.center;
        tDug = c.T;
        return true;
    });
    s.step("the pit shows on the map", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("the pit re-sampled and drawn");
        const auto edits = MapFeatures_Edits(*c.body);
        float el = 0.0f, ch = 0.0f;
        if (!edits || !MapFeatures_EditedAt(*edits, pitDir, el, ch) || ch > -1.0f) return false;
        Color now;
        if (!Nav_MapPixel(true, pitDir, now)) return false;
        const int diff = std::max({ abs(now.r - pitBefore.r), abs(now.g - pitBefore.g), abs(now.b - pitBefore.b) });
        if (diff < 20) return false; // (the tile landed; the map re-renders within a second)
        c.metric("pit on the map after", c.T - tDug, "s");
        c.lt("pit drawn within the refresh window (s)", c.T - tDug, 15.0);
        c.near("the pit's sampled depth (m)", ch, -2.5, 0.6);
        c.check("the pit is drawn darker (dug earth)", now.r + now.g + now.b < pitBefore.r + pitBefore.g + pitBefore.b,
                TextFormat("before %d,%d,%d now %d,%d,%d", pitBefore.r, pitBefore.g, pitBefore.b, now.r, now.g, now.b));
        // 6 m north of its centre the ground is untouched.
        const Vector3 away = AtCell(-12.0f, 16.0f);
        float el2 = 0.0f, ch2 = 0.0f;
        const bool sampled = MapFeatures_EditedAt(*edits, away, el2, ch2);
        c.check("untouched ground beside it keeps its height", !sampled || fabsf(ch2) < 0.25f, TextFormat("change %.2f m", ch2));
        c.grab("big-map-pit");
        tFlood = -1.0;
        poured = 0.0f;
        return true;
    });
    s.step("the flooded pit shows water on the map", 40.0, [](ScenarioCtx& c) {
        // (part 3) 1 m^3 every 0.1 s into the pit's bottom, 12 m^3 in all: the simulation's water, not generated water.
        if (poured < 12.0f && c.t >= nextPour) {
            nextPour = c.t + 0.1;
            poured += Fluid::Pour(*c.body, Vector3Subtract(pitCentre, Vector3Scale(pitDir, 1.5f)), 1.0f);
            if (poured >= 12.0f) tFlood = c.T;
        }
        c.waitingFor("the pit's water sampled and drawn");
        if (tFlood < 0.0) return false;
        const auto water = MapFeatures_Water(*c.body);
        float depth = 0.0f;
        int kind = -1;
        const float ground = c.field->radius + c.field->elevation(pitDir, 0.5f) - 2.4f; // about the pit's floor
        if (!water || !MapFeatures_WaterAt(*water, Vector3Scale(pitDir, ground), depth, kind) || depth < 0.3f) return false;
        Color now;
        if (!Nav_MapPixel(true, pitDir, now) || !(now.b > now.r && now.b > now.g)) return false; // drawn as water (blue)
        c.metric("pond on the map after", c.T - tFlood, "s");
        c.lt("pond drawn within the refresh window (s)", c.T - tFlood, 25.0);
        c.eq("its fluid is water", kind, 0);
        c.gt("its sampled depth over the pit floor (m)", depth, 0.3);
        c.check("the pond is drawn as water", true, TextFormat("%d,%d,%d", now.r, now.g, now.b));
        c.grab("big-map-pond");
        return true;
    });
    s.step("hold for the recording", 3.0, [](ScenarioCtx& c) {
        if (c.t < 2.0) return false;
        Nav_CloseMap();
        // The maps' main-thread work over the whole run (plot builds landing, renders started and uploaded): never a
        // millisecond in one frame (#1795). Timing gates are only reported on the runner.
        for (const char* z : { "gu: Nav_Update", "nav: plots", "nav: edits", "nav: water", "nav: map start", "nav: map take", "nav: minimap start", "nav: minimap take" }) {
            const double ms = Profiler_ZoneMaxMs(z);
            if (ms >= 0.0) c.timingLt(std::string(z) + " max (ms)", ms, 1.0);
        }
        return true;
    });
}
