#include "profiler.h"
#include "machinemodels.h"
#include "machinetiers.h"
#include "machineact.h"
#include "buildtool.h"
#include "researchgate.h"
#include "playerinput.h"
#include "waterworks.h"
#include "survival.h"
#include "farming.h"
#include "audio.h"
#include "gamerng.h"
#include "gameui.h"
#include "hotbar.h"
#include "inventory.h"
#include "player.h"
#include "terraintool.h"
#include "world.h"
#include "worldtypes.h"
#include "simclock.h"
#include "factory.h"
#include "alchemy.h" // (#1254) phials on vessels, the run on their panel
#include "structures.h" // (#1641 slice 2c) breakable props
#include "blueprint.h"
#include <deque>
#include "claims.h"
#include "playerrecord.h"
#include "menus.h"
#include "net.h"
#include "hostserver.h" // (#1684) Host_BroadcastPartEdits
#include "multiblock.h" // (#1429 PR 3)
#include <ctime>

#include "uipalette.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>

// (buildTool: the local PlayerSession's, src/space/playersession.cpp, #313)

namespace {
    // The current aim: a placement (plot may be null = a new plot at newOrigin) and/or a part to remove.
    struct Aim {
        Planetoid* body = nullptr;
        Plot* plot = nullptr;           // existing plot to place into (null: start a new one)
        Vector3 newOrigin{};            // body-local, for a new plot
        Part part{};
        bool canPlace = false;
        Plot* removePlot = nullptr;
        int removePart = -1;
    } aim;

    // A throwaway plot frame for previewing a new plot's first part.
    Plot preview;

    Vector3 PlotDir(const Plot& p, Int3 v) {
        return Vector3Add(Vector3Add(Vector3Scale(p.east, (float)v.x), Vector3Scale(p.up, (float)v.y)), Vector3Scale(p.north, (float)v.z));
    }
}

bool IsMachineItem(const std::string& id) {
    for (int t = 0; t < (int)PartType::Count; ++t) if (id == PartCostItem((PartType)t) && IsMachine((PartType)t)) return true;
    { PartType vt; int vs = 0; if (MachineTiers_ItemPart(id, vt, vs)) return true; } // (#1427 L7a) a tiered variant
    return false;
}

// Foundation (#210): the terrain inside the part's cells is subtracted, exactly those boxes (a carve, as mining, but with
// no yield), and nothing else changes: no fill under it, no brush reaching past it (the plot carries parts above ground).
namespace {
    bool foundationBatching = false;            // (#1167) between Begin and EndFoundationBatch
    std::vector<TerrainEdit> foundationBatch;
    // (#1183) A blueprint's foundations land over the next frames: the edits wait here (with their body) and each frame
    // applies them in small groups until kFoundationFrameMs is spent; the parts are placed at once. (The terrain edit was
    // ~11.7 of a 60-part placement's ~11.8 ms in one frame.) RAYLIBGAME_BLUEPRINT_ASYNC=0: all in the placing frame (A/B).
    struct PendingFoundation { uint64_t body; TerrainEdit edit; };
    std::deque<PendingFoundation> foundationsPending;
    // A frame's share for foundations (#1372, perf/blueprint-place with the world idle, 2026-11-21): a 60-part placement all in
    // one frame was a ~12 ms frame; queued at 3 ms its worst frame after was +5.4 ms over a ~2.7 ms baseline (5 frames), at
    // 1 ms +3.9 ms (14 frames, ~0.23 s). An edit costs ~0.16 ms. RAYLIBGAME_FOUNDATION_MS overrides (A/B).
    const double kFoundationFrameMs = [] { const char* e = getenv("RAYLIBGAME_FOUNDATION_MS"); return e ? atof(e) : 1.0; }();
    constexpr size_t kFoundationGroup = 4;      // edits applied together (one snapshot each); 2 measured no better (#1183)
}
void BuildTool_BeginFoundationBatch() { // RAYLIBGAME_BLUEPRINT_BATCH=0: a batch per part, as before (A/B)
    static const bool on = !(getenv("RAYLIBGAME_BLUEPRINT_BATCH") && getenv("RAYLIBGAME_BLUEPRINT_BATCH")[0] == '0');
    foundationBatching = on; foundationBatch.clear();
}
void BuildTool_EndFoundationBatch(Planetoid& body) {
    PROFILE_SCOPE("blueprint foundations");
    foundationBatching = false;
    static const bool async = !(getenv("RAYLIBGAME_BLUEPRINT_ASYNC") && getenv("RAYLIBGAME_BLUEPRINT_ASYNC")[0] == '0');
    if (async) for (const TerrainEdit& e : foundationBatch) foundationsPending.push_back({ body.id, e });
    else if (!foundationBatch.empty()) ApplyTerrainEdits(body, foundationBatch.data(), foundationBatch.size());
    foundationBatch.clear();
}
int BuildTool_PendingFoundations() { return (int)foundationsPending.size(); }
void BuildTool_ApplyPendingFoundations(double budgetMs) {
    if (foundationsPending.empty()) return;
    PROFILE_SCOPE("blueprint foundations (queued)");
    const double t0 = GetTime();
    // (#1372) the next group only if it fits the budget by the groups' cost so far (checking only before each group overran a
    // 3 ms budget to 4); at least one group a frame, so the queue always drains
    static double groupMs = 0.0; // a running average
    std::vector<TerrainEdit> group;
    int groups = 0;
    while (!foundationsPending.empty()) {
        const double spent = (GetTime() - t0) * 1000.0;
        if (budgetMs > 0.0 && groups > 0 && spent + groupMs > budgetMs) break;
        const uint64_t id = foundationsPending.front().body;
        group.clear();
        while (!foundationsPending.empty() && foundationsPending.front().body == id && group.size() < kFoundationGroup) { group.push_back(foundationsPending.front().edit); foundationsPending.pop_front(); }
        const double g0 = GetTime();
        if (Planetoid* b = Planetoid::find(id)) ApplyTerrainEdits(*b, group.data(), group.size()); // (a body gone: its edits with it)
        const double ms = (GetTime() - g0) * 1000.0;
        groupMs = groupMs > 0.0 ? groupMs * 0.8 + ms * 0.2 : ms;
        ++groups;
    }
}

void BuildTool_LayFoundation(Planetoid& body, Plot& plot, const Part& part) {
    std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(part, fp);
    // The plot's axes in body-local space.
    const Vector3 o = plot.toBody({ 0, 0, 0 });
    const Vector3 ax = Vector3Normalize(Vector3Subtract(plot.toBody({ 1, 0, 0 }), o));
    const Vector3 ay = Vector3Normalize(Vector3Subtract(plot.toBody({ 0, 1, 0 }), o));
    const Vector3 az = Vector3Normalize(Vector3Subtract(plot.toBody({ 0, 0, 1 }), o));
    const float h = 0.5f * kCellSize;
    std::vector<TerrainEdit> edits;
    for (const auto& c : fp) {
        const Vector3 centre = plot.toBody({ (c.first.x + 0.5f) * kCellSize, (c.first.y + 0.5f) * kCellSize, (c.first.z + 0.5f) * kCellSize });
        // Anything solid in the cell? (its centre and corners)
        bool solid = body.densityAt(centre) >= ISO_LEVEL;
        for (int k = 0; k < 8 && !solid; ++k) {
            const Vector3 corner = Vector3Add(centre, Vector3Add(Vector3Scale(ax, (k & 1) ? h * 0.95f : -h * 0.95f),
                                   Vector3Add(Vector3Scale(ay, (k & 2) ? h * 0.95f : -h * 0.95f), Vector3Scale(az, (k & 4) ? h * 0.95f : -h * 0.95f))));
            solid = body.densityAt(corner) >= ISO_LEVEL;
        }
        if (!solid) continue;
        TerrainEdit e{};
        e.center = centre;
        e.carve = true;
        e.box = true;
        e.half = { h, h, h };
        e.ax = ax; e.ay = ay; e.az = az;
        e.radius = h * 1.7320508f; // its bounding reach (dirty chunks, trees, fluid)
        edits.push_back(e);
    }
    if (foundationBatching) { foundationBatch.insert(foundationBatch.end(), edits.begin(), edits.end()); return; } // (#1167)
    if (!edits.empty()) ApplyTerrainEdits(body, edits.data(), edits.size());
}

BuildToolNetEdit gBuildToolNetEdit;
std::function<void(Planetoid& body, Plot& plot, int partIndex)> gFormRequest;

// (#1429 PR 3) "Large tank formed": the toast for a controller going formed or unformed. An unformed controller that was
// broken itself says nothing (it is gone).
static void ToastFormed(const Plot& plot, int partIndex, bool formed) {
    if (partIndex < 0 || partIndex >= (int)plot.parts.size()) return;
    const MbDef* d = Multiblock_DefForController(plot.parts[partIndex]);
    if (!d) return;
    GameUI_Toast(TextFormat("%s %s", d->name.c_str(), formed ? "formed" : "unformed"));
}

void BuildTool_NoteFormed(const Planetoid& body, const Plot& plot, int partIndex, bool formed) {
    if (partIndex < 0 || partIndex >= (int)plot.parts.size() || player.nearest != &body) return;
    const Int3 o = plot.parts[partIndex].origin;
    const Vector3 at = plot.toBody({ (o.x + 0.5f) * kCellSize, (o.y + 0.5f) * kCellSize, (o.z + 0.5f) * kCellSize });
    if (Vector3Distance(at, body.toLocalD(player.position).toFloat()) > 48.0f) return; // a structure out of sight: no toast
    ToastFormed(plot, partIndex, formed);
}

// After a part was placed or broken in this process's own world (solo, or a host's own player): the structures that no longer
// validate unform; a host tells its clients (a joined client's server does this, in the order of its edits).
static void RecheckLocal(Planetoid& body, Plot& plot) {
    if (plot.formed.empty()) return;
    for (int idx : Multiblock_Recheck(plot)) {
        ToastFormed(plot, idx, false);
        if (gNetRole == NetRole::Host) Host_BroadcastForm(body, plot, idx, false);
    }
}

// The wrench on a controller, in this process's own world (solo, or a host's own player).
static void FormLocally(Planetoid& body, Plot& plot, int partIndex) {
    const MbResult r = Multiblock_CheckController(plot, partIndex);
    if (!r.formed) { GameUI_Toast(TextFormat("Can't form: %s", r.problem.c_str())); Audio_Play(Sfx::Error, 0.5f); return; }
    Multiblock_SetFormed(plot, partIndex, true);
    ToastFormed(plot, partIndex, true);
    Audio_Play(Sfx::Build, 0.5f, 1.2f);
    if (gNetRole == NetRole::Host) Host_BroadcastForm(body, plot, partIndex, true);
}

bool BuildTool_ApplyPart(Planetoid& body, Plot& plot, const Part& part) {
    if (!plot.canPlace(part)) return false;
    if (IsBeam(part.type)) { // a line of single cubes (#200): quicker to build, each removed on its own
        static std::vector<Part> cubes;
        Plot::SplitBeam(part, cubes);
        for (const Part& c : cubes) plot.place(c);
    } else plot.place(part);
    BuildTool_LayFoundation(body, plot, part);
    return true;
}

// (#1684) The part at `i` again with a new material; the edit is recorded (a remove, then the place at the index it lands on).
static void RematerialPart(Plot& plot, int i, uint8_t material, std::vector<PartEdit>* edits) {
    Part q = plot.parts[i];
    plot.remove(i);
    if (edits) edits->push_back({ true, i, {} });
    q.material = material;
    const int at = plot.place(q);
    if (edits) edits->push_back({ false, at, q });
}

int Dungeon_UnlockVault(Plot& plot, int partIndex, std::vector<PartEdit>* edits) {
    if (partIndex < 0 || partIndex >= (int)plot.parts.size() || !GateLocked(plot.parts[partIndex])) return 0;
    const Part g = plot.parts[partIndex];
    int opened = 0;
    for (int i = 0; i < (int)plot.parts.size(); ++i) { // both gates of the door (stacked: same column, 2 cells apart)
        const Part& q = plot.parts[i];
        if (!GateLocked(q) || q.origin.x != g.origin.x || q.origin.z != g.origin.z || abs(q.origin.y - g.origin.y) > 2) continue;
        RematerialPart(plot, i, 0, edits);
        ++opened;
    }
    return opened;
}

int Dungeon_PullLever(Plot& plot, int partIndex, std::vector<PartEdit>* edits) {
    if (partIndex < 0 || partIndex >= (int)plot.parts.size() || plot.parts[partIndex].type != PartType::Lever || plot.parts[partIndex].material != 0) return -1;
    const Part lever = plot.parts[partIndex];
    std::vector<int> gates; // (collected first: opening replaces the parts)
    for (size_t i = 0; i < plot.parts.size(); ++i) {
        const Part& q = plot.parts[i];
        if (!GatePortcullis(q)) continue;
        const float dx = (float)(q.origin.x - lever.origin.x), dy = (float)(q.origin.y - lever.origin.y), dz = (float)(q.origin.z - lever.origin.z);
        if (sqrtf(dx * dx + dy * dy + dz * dz) * kCellSize <= 12.0f) gates.push_back((int)i);
    }
    RematerialPart(plot, partIndex, 1, edits); // pulled: the lever is placed again, material 1
    for (int i : gates) RematerialPart(plot, i, 0, edits);
    return (int)gates.size();
}

// (#900) Claims in a game this player runs (single player, or a host's own player): the player's team (a solo player is a
// team of one: made the first time it is needed), and the map's mirror of the claims, which a joined client gets from its
// server (MsgClaims) and a game run here builds from its own list when that changes.
static int LocalTeam() {
    int team = Teams_Of(PlayerId_Local());
    if (!team) team = Teams_CreateSolo(PlayerId_Local(), settings.mpName);
    return team;
}
static void RefreshLocalClaims() {
    if (gNetRole == NetRole::Client) return;
    static uint32_t seen = 0xFFFFFFFFu;
    if (Claims_Version() == seen) return;
    seen = Claims_Version();
    const int me = Teams_Of(PlayerId_Local());
    std::vector<ClaimView> list;
    for (const Claim& c : Claims_All()) {
        const TeamRecord* t = Teams_Find(c.team);
        ClaimView v;
        v.id = c.id; v.team = c.team; v.teamName = t && !t->name.empty() ? t->name : TextFormat("team %d", c.team);
        v.colour = Claims_TeamColour(c.team, t ? t->colour : 0);
        v.body = c.body; v.centre = c.centre; v.radius = c.radius;
        v.ours = me != 0 && c.team == me; v.allowed = me != 0 && Claims_Allows(c.team, me);
        list.push_back(std::move(v));
    }
    Claims_SetMirror(list);
}

bool BuildTool_Place(Planetoid& body, Plot& plot, const Part& placed, bool charge) {
    // (#1427 L7a) a machine's material is its variant's step (0 for untiered ones: the shape material the build tool last held
    // never rides along)
    Part part = placed;
    if (IsMachine(part.type)) part.material = (uint8_t)std::clamp((int)part.material, 0, MachineTiers_MaxStep(part.type));
    if (part.type == PartType::Farmland) { GameUI_Toast("Farmland is tilled ground now: till soil with a shovel or hoe (RMB)"); return false; } // (#1085)
    if (!plot.canPlace(part)) return false;
    const PartType t = part.type;
    if (charge && ResearchGate_On() && !terrainTool.creative && !ResearchGate_PartAllowed(ResearchGate_Team(), t)) { // (#1271 L6b: a player's placing; creative and scripted placements are free)
        const ResearchNode* n = ResearchGate_NodeForPart(t);
        GameUI_Toast(TextFormat("Research %s first (Research tab)", n ? n->name.c_str() : "it"));
        Audio_Play(Sfx::Error, 0.5f);
        return false;
    }
    const char* item = PartItem(part);
    // (#900) A claim banner claims its ground for the placer's team here too, as a server does for a joined player (or it
    // is refused, nothing spent: over the budget, or over another team's claim).
    uint32_t claim = 0;
    if (const int tier = ClaimBannerTier(t); tier >= 0 && gNetRole != NetRole::Client) {
        const int team = LocalTeam();
        const Vector3 cellCentre{ (part.origin.x + 0.5f) * kCellSize, (part.origin.y + 0.5f) * kCellSize, (part.origin.z + 0.5f) * kCellSize };
        const int reached = std::max(1, (int)QuestBook_LocalState().visitedBodies.size()); // (the bodies stood on: the budget grows with them)
        std::string why;
        claim = Claims_Make(team, Claims_Budget(reached), body.name, plot.toBody(cellCentre), kClaimBannerAcross[tier] * 0.5f, (int64_t)time(nullptr), why);
        if (!claim) { GameUI_Toast(why.c_str()); Audio_Play(Sfx::Error, 0.5f); return false; }
    }
    if (charge && !terrainTool.creative && !inventory.remove(item, PartCostCount(t))) {
        GameUI_Toast(TextFormat("Need %d %s for a %s", PartCostCount(t), item, PartName(t)));
        Audio_Play(Sfx::Error, 0.5f);
        if (claim) Claims_Remove(claim);
        return false;
    }
    BuildTool_ApplyPart(body, plot, part);
    RecheckLocal(body, plot); // (#1429 PR 3) a part in a structure's air cell, or in the way
    if (claim) {
        Claims_SetBanner(claim, plot.id, plot.partAt(part.origin));
        GameUI_Toast(TextFormat("Claimed the land %.0f m across around the banner (the map shows it)", kClaimBannerAcross[ClaimBannerTier(t)]));
    }
    if (IsMachine(part.type) || IsRocketPart(part.type) || IsRotationPart(part.type)) QuestBook_NotePlaced(PartCostItem(part.type));
    if (IsMachine(part.type) && part.material > 0) QuestBook_NotePlaced(item); // (#1427) a variant counts as its base and as itself
    return true;
}

// Breaking built parts: LMB held on a part for BreakHold() seconds (in build mode or out of it; a
// part in front of the ground takes the click from the terrain tool). A short hold, so a stray click
// never deletes a wall.
namespace {
    constexpr float kBreakHold0 = 0.45f;
    float BreakHold() { return kBreakHold0 / std::max(0.2f, SkillStat("break_speed")); } // (#919) Building's perks: a shorter hold
    const Plot* breakPlot = nullptr;
    int breakPart = -1;
    float breakHeld = 0.0f;
    bool breakAiming = false;
}

float BuildTool_BreakProgress() { return breakAiming && breakHeld > 0.0f ? std::min(1.0f, breakHeld / BreakHold()) : 0.0f; }
float BuildTool_BreakHold() { return BreakHold(); }
bool BuildTool_AimingPart() { return breakAiming; }

int BuildTool_AutoOrient(int partType, Int3 target) {
    const PartType type = (PartType)partType;
    // The part's front - its local +z, where machines have their door - toward target (one of the 4 turns about up). (#1578) A
    // ladder is a thin frame against its cell's -x side: its climbing face (local +x) goes toward target, not its edge.
    const int col = type == PartType::Ladder ? 0 : 2;
    int o = 0;
    for (int k = 0; k < 4; ++k, o = RotateOrientation(o, 1)) {
        const int* m = Orientation(o);
        if (m[col] == target.x && m[3 + col] == target.y && m[6 + col] == target.z) return o;
    }
    return 0; // (straight up or down: no turn about up faces it)
}
int BuildTool_AimOrient() { return aim.plot || aim.newOrigin.x != 0.0f || aim.newOrigin.y != 0.0f || aim.newOrigin.z != 0.0f ? (int)aim.part.orient : -1; }
const Plot* BuildTool_AimPlacement(Vector3* newOrigin, Int3* origin) { if (newOrigin) *newOrigin = aim.newOrigin; if (origin) *origin = aim.part.origin; return aim.plot; }

static void UpdateBreak(Scene& world, const PlayerInput& in) {
    breakAiming = false;
    Planetoid* body = player.nearest;
    if (!body || !body->plots || player.flyMode) { breakHeld = 0.0f; breakPart = -1; return; }
    Vector3 dir, ld; DVec3 lo;
    player.aimRay(world, *body, lo, ld, dir); // (F4, #1412: the eye in double with the frame path; the camera before)
    float best = buildTool.active ? buildTool.reach : terrainTool.reach;
    Plot* hitPlot = nullptr;
    Int3 hitCell{};
    for (auto& pl : body->plots->plots) {
        Int3 c, f; float d;
        if (pl->raycast(lo, ld, best, c, f, d)) { best = d; hitPlot = pl.get(); hitCell = c; }
    }
    if (!hitPlot) { breakHeld = 0.0f; breakPart = -1; return; }
    TerrainHit th;
    if (RaycastTerrain(world, Ray{ world.camera.position, dir }, best, th) && th.planetoid == body && th.distance < best * body->scale) { breakHeld = 0.0f; breakPart = -1; return; }
    const int part = hitPlot->partAt(hitCell);
    if (part < 0) { breakHeld = 0.0f; breakPart = -1; return; }
    breakAiming = true;
    if (!gPrimaryLatch.allows(in.down(Act::Primary))) { breakHeld = std::min(breakHeld, 0.0f); breakPart = -1; return; } // (#895: not after a break, until pressed again)
    if (PartUnbreakable(hitPlot->parts[part])) { breakHeld = 0.0f; breakPart = -1; if (in.hit(Act::Primary)) GameUI_Toast(GateLocked(hitPlot->parts[part]) ? "Locked: it needs a dungeon key" : hitPlot->parts[part].type == PartType::Lever ? "The lever will not come loose" : "The portcullis is too heavy to break: look for a lever"); return; } // (#1641) a vault door, a portcullis or a lever: not broken
    if (hitPlot != breakPlot || part != breakPart) { breakPlot = hitPlot; breakPart = part; breakHeld = std::max(breakHeld, 0.0f); }
    breakHeld += GameFrameTime();
    if (breakHeld < BreakHold()) return;
    const Part removed = hitPlot->parts[part];
    if (gBuildToolNetEdit) { // joined (#330): the server removes it, refunds, and tells everyone
        gBuildToolNetEdit(*body, hitPlot, Vector3{}, removed, true, part);
        Audio_Play(Sfx::Dig, 0.5f, 1.3f);
        breakPart = -1;
        breakHeld = -0.15f;
        gPrimaryLatch.consume(); // (#895) the held button is used up
        return;
    }
    const uint32_t bannerClaim = ClaimBannerTier(removed.type) >= 0 && gNetRole != NetRole::Client ? Claims_OfBanner(body->name, hitPlot->id, part) : 0; // (#900)
    if (hitPlot->remove(part)) {
        RecheckLocal(*body, *hitPlot); // (#1429 PR 3) a structure that lost a part unforms
        if (bannerClaim) { Claims_Remove(bannerClaim); GameUI_Toast("The banner is down: the land is no longer claimed"); }
        if (removed.type == PartType::Meteorite && !terrainTool.creative) { // (#177) space rock: nickel-iron, now and then platinum or iridium ore
            const int n = GameRng_Int(6, 10), pt = GameRng_Int(0, 99) < 35 ? 2 : 0, ir = GameRng_Int(0, 99) < 12 ? 1 : 0; // the authority rolls (#316)
            inventory.add("nickel_iron_ore", n);
            if (pt) inventory.add("platinum_ore", pt);
            if (ir) inventory.add("iridium_ore", ir);
            GameUI_Toast(TextFormat("The meteorite breaks apart: +%d nickel-iron%s%s", n, pt ? ", +2 platinum ore" : "", ir ? ", +1 iridium ore" : ""));
        } else if (!terrainTool.creative) inventory.add(PartItem(removed), PartCostCount(removed.type));
        if (std::vector<std::pair<std::string, int>> drops; Structures_BreakLoot(*body, hitPlot->id, part, removed, drops)) { // (#1641 slice 2c) a dungeon urn or crate: its loot
            std::string what;
            for (const auto& [item, n] : drops) { inventory.add(item.c_str(), n); what += std::string(what.empty() ? "" : ", ") + std::to_string(n) + " " + item; }
            if (!drops.empty()) GameUI_Toast(("It breaks open: " + what).c_str());
        }
        Audio_Play(Sfx::Dig, 0.5f, 1.3f);
    }
    breakPart = -1;
    breakHeld = -0.15f;
    gPrimaryLatch.consume(); // (#895, the user: holding on mined what was behind) the next part needs a new press
}

// (#919) Building XP: a part the player placed pays once it has stood 5 minutes (placing and breaking again pays nothing).
namespace {
    struct Placed { uint64_t body; uint32_t plot; Int3 origin; PartType type; double at; }; // (the body's id: never a pointer kept)
    std::vector<Placed> placedByPlayer;
    constexpr double kPartKept = 300.0; // s
    void TickPlacedParts(Scene& world) {
        if (placedByPlayer.empty() || SimTime() - placedByPlayer.front().at < kPartKept) return;
        const Placed p = placedByPlayer.front();
        placedByPlayer.erase(placedByPlayer.begin());
        for (auto& o : world.objects) {
            auto* b = dynamic_cast<Planetoid*>(o.get());
            if (!b || b->id != p.body || !b->plots) continue;
            const Plot* plot = b->plots->byId(p.plot);
            const int idx = plot ? plot->partAt(p.origin) : -1;
            if (idx >= 0 && plot->parts[idx].type == p.type && plot->parts[idx].origin == p.origin && gSkillEvent) gSkillEvent("part_kept", 1.0, Vector3{});
            break;
        }
    }
}
void BuildTool_NotePlacedByPlayer(Planetoid& body, const Plot& plot, const Part& part) {
    if (placedByPlayer.size() < 4096) placedByPlayer.push_back({ body.id, plot.id, part.origin, part.type, SimTime() });
}

void BuildTool_BenchAimAt(Plot& plot, Int3 cell) {
    // The part's first footprint cell onto `cell`: a facing rotates the footprint about the origin cell's centre, which can
    // move it a cell off the origin (Plot::Footprint).
    aim.plot = &plot;
    aim.part.origin = { 0, 0, 0 };
    std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(aim.part, fp);
    aim.part.origin = fp.empty() ? cell : Int3{ cell.x - fp[0].first.x, cell.y - fp[0].first.y, cell.z - fp[0].first.z };
    aim.canPlace = plot.canPlace(aim.part);
}

void BuildTool_Update(Scene& world, const PlayerInput& in, bool usable) {
    TickPlacedParts(world);
    BuildTool_ApplyPendingFoundations(kFoundationFrameMs); // (#1183) a placement's foundations, a frame's share at a time
    if (!usable) { buildTool.active = false; breakAiming = false; breakHeld = 0.0f; return; }
    RefreshLocalClaims(); // (#900) the map's claims in a game run here
    UpdateBreak(world, in);
    if (in.hit(Act::BuildMode)) {
        buildTool.active = !buildTool.active;
        buildTool.orientManual = false; // parts face you again (#196)
        GameUI_Toast(buildTool.active ? "Build mode: Ctrl+wheel shape, parts face you (R/V rotate by hand), RMB place, hold LMB remove, B save blueprint, N blueprints (G to leave)" : "Build mode off");
    }
    if (buildTool.part != buildTool.autoPart) { buildTool.autoPart = buildTool.part; buildTool.orientManual = false; } // another part: automatic facing
    aim = Aim{};
    if (!buildTool.active) return;
    if (in.enabled) {
        const float wheel = in.wheel;
        // The wheel picks the shape (any shape from the selected material); machines come from the hotbar.
        const bool modWheel = IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_RIGHT_CONTROL) || IsKeyDown(KEY_LEFT_ALT) || IsKeyDown(KEY_RIGHT_ALT); // #227: the plain wheel scrolls the hotbar
        if (modWheel && IsShape((PartType)buildTool.part) && wheel != 0.0f) { // Ctrl/Alt + wheel: step through kShapes
            int k = 0;
            for (int i = 0; i < kShapeCount; ++i) if ((int)kShapes[i] == buildTool.part) k = i;
            k = (k + (wheel > 0 ? 1 : kShapeCount - 1)) % kShapeCount;
            buildTool.part = (int)kShapes[k];
        }
        if (in.hit(Act::Rotate) || in.hit(Act::Flip)) {
            if (!buildTool.orientManual) GameUI_Toast("Rotating by hand (parts face you again with the next part, or G)");
            buildTool.orientManual = true; // the orientation sticks from now on (#196); it starts from the automatic one shown
            buildTool.autoPart = buildTool.part;
            buildTool.orient = RotateOrientation(buildTool.orient, in.hit(Act::Rotate) ? 1 : 0);
        }
    }
    Planetoid* body = player.nearest;
    if (!body) return;
    Vector3 dir, ld; DVec3 lo;
    player.aimRay(world, *body, lo, ld, dir); // (F4, #1412: the eye in double with the frame path; the camera before)
    aim.body = body;
    aim.part = { (PartType)buildTool.part, { 0, 0, 0 }, (uint8_t)buildTool.orient, (uint8_t)(IsMachine((PartType)buildTool.part) || IsMultiblockPart((PartType)buildTool.part) ? buildTool.variant : buildTool.material) }; // (#1427: a machine's variant; #1429: a multiblock part's list index)

    // Nearest part face along the ray, and the terrain.
    float bestPart = buildTool.reach;
    Plot* hitPlot = nullptr;
    Int3 hitCell{}, hitFace{};
    if (body->plots)
        for (auto& pl : body->plots->plots) {
            Int3 c, f; float d;
            if (pl->raycast(lo, ld, bestPart, c, f, d)) { bestPart = d; hitPlot = pl.get(); hitCell = c; hitFace = f; }
        }
    TerrainHit th;
    const bool terrain = RaycastTerrain(world, Ray{ world.camera.position, dir }, buildTool.reach, th) && th.planetoid == body;
    // A part bigger than a cell is centred on the aimed point across the face it is placed on (#196: a 2x2 machine hung
    // off one corner of the aimed cell). `normal` is the face's axis (a plot axis); the other two axes move the
    // footprint so its middle is nearest the hit point (in cells).
    auto centreOn = [&](Plot& pl, Int3 normal, const DVec3& hitBodyLocal, bool onGround = false) {
        std::vector<std::pair<Int3, CellShape>> fp;
        Plot::Footprint(aim.part, fp);
        if (fp.size() <= 1) { aim.part.origin = Plot::SingleCellOrigin(aim.part, aim.part.origin); return; } // (#1630) its one cell where aimed, any orientation
        Int3 lo{ INT32_MAX, INT32_MAX, INT32_MAX }, hi{ INT32_MIN, INT32_MIN, INT32_MIN };
        for (const auto& c : fp) {
            lo = { std::min(lo.x, c.first.x), std::min(lo.y, c.first.y), std::min(lo.z, c.first.z) };
            hi = { std::max(hi.x, c.first.x), std::max(hi.y, c.first.y), std::max(hi.z, c.first.z) };
        }
        const Vector3 h = Vector3Scale(pl.toPlot(hitBodyLocal), 1.0f / kCellSize);
        // The shift that puts the footprint's middle nearest the hit, with hysteresis (#551): aimed near the halfway
        // point between two placements, the camera's bob flipped the box back and forth a cell every few frames. The
        // last shift is kept until the hit is 0.2 cells past the halfway point.
        // Compared in absolute cells: the aimed cell (the origin before the shift) changes as the hit moves.
        static struct { const Plot* plot = nullptr; PartType type = PartType::Count; uint8_t orient = 0; Int3 normal{}, origin{}; } last;
        const bool same = last.plot == &pl && last.type == aim.part.type && last.orient == aim.part.orient && last.normal == normal;
        constexpr float kKeep = 0.5f + 0.2f;
        auto place = [&](int base, float hv, int a, int b, int prev) {
            const float want = (float)base + hv - 0.5f * (float)(a + b + 1); // the origin that centres it exactly
            return same && fabsf(want - (float)prev) <= kKeep ? prev : base + (int)lroundf(want - (float)base);
        };
        Int3& o = aim.part.origin;
        if (normal.x == 0) o.x = place(o.x, h.x, lo.x, hi.x, last.origin.x);
        if (normal.y == 0) o.y = place(o.y, h.y, lo.y, hi.y, last.origin.y);
        if (normal.z == 0) o.z = place(o.z, h.z, lo.z, hi.z, last.origin.z);
        // On the ground (#771): the height is the cell just above the ground at the aim (the caller's o.y, from the hit
        // 0.2 cells up). On a slope that crossed a cell boundary as the aim moved, and the camera's bob flipped the box
        // up and down a cell every frame there. Keep the last height until the ground is 0.2 cells past the boundary.
        if (onGround && same && o.y != last.origin.y && abs(o.y - last.origin.y) == 1) {
            const float g = h.y + 0.2f;                                       // what o.y was taken from (cells)
            const float past = o.y > last.origin.y ? g - (float)o.y : (float)last.origin.y - g;
            if (past < 0.2f) o.y = last.origin.y;
        }
        last = { &pl, aim.part.type, aim.part.orient, normal, o };
    };
    // Automatic facing (#196, until R/V is pressed): the part's front (its local +z, where machines have their door)
    // turns toward the player - on a wall, out of the wall; on the ground or a floor, back along the view.
    auto faceToward = [&](Int3 target) {
        if (buildTool.orientManual) return;
        const int o = BuildTool_AutoOrient((int)aim.part.type, target);
        buildTool.orient = o;          // shown, and where a hand rotation starts from
        aim.part.orient = (uint8_t)o;
    };
    auto towardViewer = [&](const Plot& pl) {
        const float vx = Vector3DotProduct(ld, pl.east), vz = Vector3DotProduct(ld, pl.north);
        return fabsf(vx) > fabsf(vz) ? Int3{ vx > 0.0f ? -1 : 1, 0, 0 } : Int3{ 0, 0, vz > 0.0f ? -1 : 1 };
    };
    if (hitPlot && (!terrain || bestPart <= th.distance)) {
        aim.plot = hitPlot;
        faceToward(hitFace.y == 0 ? hitFace : towardViewer(*hitPlot));
        static const bool stackFacing = !(getenv("RAYLIBGAME_LADDER_STACK_FACING") && getenv("RAYLIBGAME_LADDER_STACK_FACING")[0] == '0'); // (A/B)
        if (stackFacing && !buildTool.orientManual && aim.part.type == PartType::Ladder && hitFace.y != 0) { // (#1640) stacked on (or under) a ladder: its facing
            const int on = hitPlot->partAt(hitCell);
            if (on >= 0 && hitPlot->parts[on].type == PartType::Ladder) { aim.part.orient = hitPlot->parts[on].orient; buildTool.orient = aim.part.orient; }
        }
        aim.part.origin = { hitCell.x + hitFace.x, hitCell.y + hitFace.y, hitCell.z + hitFace.z };
        centreOn(*hitPlot, hitFace, lo + DVec3(Vector3Scale(ld, bestPart)));
        // A part bigger than a cell grows from its origin, and its rotation can reach back into the block aimed at or stand a
        // cell off it: flush against the aimed face, every cell in front of it (#1574: it slid only outward).
        aim.part.origin = Plot::FlushAgainstFace(aim.part, hitCell, hitFace);
        aim.removePlot = hitPlot;
        aim.removePart = hitPlot->partAt(hitCell);
    } else if (terrain) {
        const Vector3 ground = th.local;
        // The part's extent (cells, from its origin): a plot whose cells come within 6 m of any of it takes it, so a big
        // part never lands on a second, unaligned grid overlapping the first one's parts (#551).
        std::vector<std::pair<Int3, CellShape>> fp;
        Plot::Footprint(aim.part, fp);
        Int3 flo{ 0, 0, 0 }, fhi{ 0, 0, 0 };
        for (const auto& c : fp) {
            flo = { std::min(flo.x, c.first.x), std::min(flo.y, c.first.y), std::min(flo.z, c.first.z) };
            fhi = { std::max(fhi.x, c.first.x), std::max(fhi.y, c.first.y), std::max(fhi.z, c.first.z) };
        }
        const float extent = (float)std::max({ fhi.x - flo.x, fhi.y - flo.y, fhi.z - flo.z }) * kCellSize;
        Plot* near = body->plots ? body->plots->plotFor(ground, 6.0f + extent, false) : nullptr;
        if (near) {
            aim.plot = near;
            faceToward(towardViewer(*near));
            // The cell just above the ground at the aim.
            aim.part.origin = near->cellAt(near->toPlot(Vector3Add(ground, Vector3Scale(near->up, 0.1f))));
            centreOn(*near, Int3{ 0, 1, 0 }, DVec3(ground), true); // on the ground: the bottom face's middle at the aim
        } else {
            // A new plot (#550): its grid starts at the aim, so the part grew from its corner there. Move the plot's
            // origin back along its east and north axes (as PlotSet::plotFor will set them) by half the footprint, so
            // the middle of the part's bottom face sits on the aim point.
            static Plot frame;
            frame.up = Vector3Normalize(ground);
            frame.east = Vector3CrossProduct({ 0, 1, 0 }, frame.up);
            if (Vector3Length(frame.east) < 0.1f) frame.east = { 1, 0, 0 };
            frame.east = Vector3Normalize(frame.east);
            frame.north = Vector3CrossProduct(frame.up, frame.east);
            faceToward(towardViewer(frame));
            aim.part.origin = { 0, 0, 0 };
            fp.clear();
            Plot::Footprint(aim.part, fp); // in the facing just chosen
            Vector3 lo{ 1e9f, 1e9f, 1e9f }, hi{ -1e9f, -1e9f, -1e9f };
            for (const auto& c : fp) {
                lo = Vector3Min(lo, Vector3{ (float)c.first.x, (float)c.first.y, (float)c.first.z });
                hi = Vector3Max(hi, Vector3{ (float)c.first.x + 1.0f, (float)c.first.y + 1.0f, (float)c.first.z + 1.0f });
            }
            const float mx = 0.5f * (lo.x + hi.x) * kCellSize, mz = 0.5f * (lo.z + hi.z) * kCellSize;
            aim.newOrigin = Vector3Subtract(ground, Vector3Add(Vector3Scale(frame.east, mx), Vector3Scale(frame.north, mz)));
            aim.newOrigin = Vector3Subtract(aim.newOrigin, Vector3Scale(frame.up, lo.y * kCellSize)); // its bottom on the ground
        }
    } else {
        return;
    }
    if (aim.plot) aim.canPlace = aim.plot->canPlace(aim.part);
    else aim.canPlace = true;

    if (!in.enabled) return;
    if (Blueprint_Keys(body, aim.removePlot)) return; // B saves the plot looked at, N cycles blueprints
    if (gBuildToolNetEdit && in.hit(Act::Secondary) && aim.canPlace) { // joined (#330): the server places and charges
        if (Blueprint_Active()) { GameUI_Toast("Blueprints can't be placed on a server yet"); return; }
        gBuildToolNetEdit(*body, aim.plot, aim.newOrigin, aim.part, false, -1);
        Audio_Play(Sfx::Build, 0.5f, 1.2f);
        return;
    }
    if (const Blueprint* bp = Blueprint_Active(); bp && in.hit(Act::Secondary)) {
        Plot* target = aim.plot;
        if (!target) {
            if (!body->plots) body->plots = std::make_unique<PlotSet>();
            target = body->plots->plotFor(aim.newOrigin, 0.0f, true);
        }
        std::string missing;
        // (#1018) What can't be paid for stays as ghosts, built from the inventory when you come by; placing the same
        // blueprint on its own ghosts again (nothing new to build) cancels them.
        const int ghostsBefore = Blueprint_GhostsOn(*target, *bp, aim.part.origin);
        const int n = Blueprint_Place(*body, *target, *bp, aim.part.origin, true, missing, true);
        if (n == 0 && ghostsBefore > 0) {
            const int gone = Blueprint_RemoveGhosts(*target, *bp, aim.part.origin);
            GameUI_Toast(TextFormat("Cancelled %d ghost part%s of %s", gone, gone == 1 ? "" : "s", bp->name.c_str()));
            return;
        }
        GameUI_Toast(missing.empty() ? TextFormat("Placed %d parts of %s", n, bp->name.c_str())
                                     : TextFormat("Placed %d parts; the rest are ghosts, built when you come by with %s (place again to cancel)", n, missing.c_str()));
        if (n > 0) Audio_Play(Sfx::Build, 0.5f, 1.0f);
        return;
    }
    if (in.hit(Act::Secondary) && aim.canPlace) {
        Plot* target = aim.plot;
        if (!target) {
            if (!body->plots) body->plots = std::make_unique<PlotSet>();
            target = body->plots->plotFor(aim.newOrigin, 0.0f, true);
        }
        if (BuildTool_Place(*body, *target, aim.part, true)) { Audio_Play(Sfx::Build, 0.5f, 1.2f); BuildTool_NotePlacedByPlayer(*body, *target, aim.part); }
    } // LMB on a part: UpdateBreak (a short hold)
}

void BuildTool_Draw3D(Scene& world) {
    if (!buildTool.active || !aim.body) return;
    Planetoid& body = *aim.body;
    const Plot* frame = aim.plot;
    if (!frame) {
        preview.origin = DVec3(aim.newOrigin);
        preview.up = Vector3Normalize(aim.newOrigin);
        preview.east = Vector3CrossProduct({ 0, 1, 0 }, preview.up);
        if (Vector3Length(preview.east) < 0.1f) preview.east = { 1, 0, 0 };
        preview.east = Vector3Normalize(preview.east);
        preview.north = Vector3CrossProduct(preview.up, preview.east);
        frame = &preview;
    }
    auto box = [&](const Plot& p, Int3 c, Color col) {
        // In plot coordinates under the plot's own double-precision model, as the blocks are drawn: float corners
        // 60 km out drifted off them and shook as the planet turned (#271).
        Vector3 v[8];
        for (int i = 0; i < 8; ++i) {
            // Inflated 2 cm so the lines stand off the faces of the blocks they outline (they z-fought, #199).
            constexpr float kPad = 0.02f;
            auto edge = [](int cell, int hi) { return (cell + hi) * kCellSize + (hi ? kPad : -kPad); };
            v[i] = { edge(c.x, i & 1), edge(c.y, (i >> 1) & 1), edge(c.z, (i >> 2) & 1) };
        }
        rlDrawRenderBatchActive();
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(p.cellModel(body, world.renderOrigin)));
        const int e[12][2] = { { 0, 1 }, { 2, 3 }, { 4, 5 }, { 6, 7 }, { 0, 2 }, { 1, 3 }, { 4, 6 }, { 5, 7 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
        for (auto& ed : e) DrawLine3D(v[ed[0]], v[ed[1]], col);
        rlPopMatrix();
    };
    // The part's real shape (#197): each footprint cell's faces (slopes, corners, beams as they will be), translucent
    // green (or red: can't place) with their edges, inflated 2 cm about the cell centre so they stand off the blocks
    // around them (#199), both sides drawn.
    // cross: can't place - each face also gets its diagonals, a cue that doesn't depend on telling the colours apart (#1432).
    auto ghost = [&](const Plot& p, const std::vector<std::pair<Int3, CellShape>>& cells, Color fill, Color edge, bool cross) {
        std::vector<std::vector<Vector3>> faces;
        rlDrawRenderBatchActive();
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(p.cellModel(body, world.renderOrigin)));
        rlDisableBackfaceCulling();
        rlDisableDepthMask();
        constexpr float kGrow = 1.0f + 0.04f / kCellSize; // 2 cm each side
        for (const auto& c : cells) {
            CellShapeFaces(c.second, faces);
            const Vector3 lo = { (float)c.first.x, (float)c.first.y, (float)c.first.z };
            auto toPlot = [&](Vector3 v) { // cell units -> plot metres, grown about the cell centre
                const Vector3 g = Vector3Add(Vector3Scale(Vector3Subtract(v, { 0.5f, 0.5f, 0.5f }), kGrow), { 0.5f, 0.5f, 0.5f });
                return Vector3Scale(Vector3Add(lo, g), kCellSize);
            };
            for (const auto& poly : faces) {
                for (size_t i = 1; i + 1 < poly.size(); ++i) DrawTriangle3D(toPlot(poly[0]), toPlot(poly[i]), toPlot(poly[i + 1]), fill);
                for (size_t i = 0; i < poly.size(); ++i) DrawLine3D(toPlot(poly[i]), toPlot(poly[(i + 1) % poly.size()]), edge);
                if (cross && poly.size() >= 3) { // corner to the opposite corner (a triangle: each corner to the far edge's middle)
                    const size_t n = poly.size();
                    if (n >= 4) for (size_t i = 0; i < n / 2; ++i) DrawLine3D(toPlot(poly[i]), toPlot(poly[i + n / 2]), edge);
                    else for (size_t i = 0; i < 3; ++i) DrawLine3D(toPlot(poly[i]), toPlot(Vector3Scale(Vector3Add(poly[(i + 1) % 3], poly[(i + 2) % 3]), 0.5f)), edge);
                }
            }
        }
        rlDrawRenderBatchActive();
        rlEnableDepthMask();
        rlEnableBackfaceCulling();
        rlPopMatrix();
    };
    std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(aim.part, fp);
    if (Blueprint_Active()) Blueprint_DrawGhost(body, *frame, aim.part.origin, true, world.renderOrigin); // the whole blueprint instead of one part
    else ghost(*frame, fp, aim.canPlace ? Fade(Ui().good, 0.28f) : Fade(Ui().bad, 0.30f), aim.canPlace ? Ui().good : Ui().bad, !aim.canPlace);
    if (aim.removePlot && aim.removePart >= 0) {
        Plot::Footprint(aim.removePlot->parts[aim.removePart], fp);
        for (const auto& c : fp) box(*aim.removePlot, c.first, Fade(ORANGE, 0.7f));
    }
}

void BuildTool_DrawHUD() {
    {   // The bar under the crosshair: a part breaking (amber, filling) or the swing recovering (white).
        const int bx = GetScreenWidth() / 2 - 30, by = GetScreenHeight() / 2 + 20, bw = 60;
        const float brk = BuildTool_BreakProgress(), cd = buildTool.active ? 0.0f : TerrainTool_Cooldown();
        if (brk > 0.0f || cd > 0.0f) {
            DrawRectangle(bx - 1, by - 1, bw + 2, 6, Fade(BLACK, 0.55f));
            if (brk > 0.0f) DrawRectangle(bx, by, (int)(bw * brk), 4, Color{ 255, 180, 60, 255 });
            else DrawRectangle(bx, by, (int)(bw * (1.0f - cd)), 4, Fade(RAYWHITE, 0.85f));
        }
    }
    if (!buildTool.active) return;
    Blueprint_DrawHUD();
    const PartType t = (PartType)buildTool.part;
    const int cx = GetScreenWidth() / 2, cy = GetScreenHeight() / 2;
    DrawLine(cx - 8, cy, cx + 8, cy, RAYWHITE);
    DrawLine(cx, cy - 8, cx, cy + 8, RAYWHITE);
    const Part p{ t, { 0, 0, 0 }, 0, (uint8_t)(IsMachine(t) || IsMultiblockPart(t) ? buildTool.variant : buildTool.material) };
    const char* item = PartItem(p);
    const char* label = IsShape(t) ? TextFormat("%s %s", GetBuildMaterial(buildTool.material).name, PartName(t)) : PartName(t);
    const char* line = TextFormat("BUILD  %s  (%d %s, have %d)", label, PartCostCount(t), item, inventory.count(item));
    DrawText(line, cx - MeasureText(line, 22) / 2, GetScreenHeight() - 150, 22, SKYBLUE);
    DrawText(IsShape(t) ? "Ctrl+wheel shape   wheel hotbar   R / V rotate   RMB place   LMB remove   G leave   1-0 hotbar"
                        : "R / V rotate   RMB place   LMB remove   G leave   1-0 hotbar",
             10, GetScreenHeight() - 30, 20, Fade(RAYWHITE, 0.7f));
}

// ---------------------------------------------------------------------------------------------
// Machines: look at one to see its state; I puts matching items in, O takes the outputs.
// ---------------------------------------------------------------------------------------------

namespace {
    struct MachineAim { uint64_t body = 0; uint32_t plot = 0; int part = -1; PartType type = PartType::Count; Int3 cell{ 0, 0, 0 }; } machineAim;
    struct FarmAim { uint64_t body = 0; Vector3 local{}; bool on = false; } farmAim; // (#1085) farmland ground or a crop, aimed
    double handPumpAt = -10.0; // the last frame the hand pump was worked (GetTime), and the litres since this go began (#571)
    float handPumpLitres = 0.0f;
    double crankAt = -10.0; // (#1269) the last frame a hand crank was turned (GetTime)
}
bool BuildTool_AimingCampfire() { return machineAim.part >= 0 && machineAim.type == PartType::Campfire; }

namespace {
// (#1716) Does RMB on the aimed part do something of its own (open a screen, work a pump, pull a lever...)? Cables and the
// look-only rotation parts (shafts, cogs, gauges) do nothing on RMB, so the held item's RMB may fire at them (a block placed
// against a shaft). In build mode RMB places parts instead, and no interactable takes it.
bool AimTakesRmb() {
    if (buildTool.active || machineAim.part < 0) return false;
    if (machineAim.type == PartType::Cable) return false;
    return !(IsRotationPart(machineAim.type) && !IsMachine(machineAim.type));
}
enum class RmbStart { None, Empty, Interactable };
RmbStart rmbStart = RmbStart::None; // where the RMB press that is still held began (as of the last Machine_Update)
void MachineInteract(Scene& world, const PlayerInput& in);
}

bool Interact_RmbConsumed() {
    if (rmbStart == RmbStart::None) return AimTakesRmb(); // not held at the last update: a press now lands on what is aimed
    return rmbStart == RmbStart::Interactable;
}

PlayerInput Interact_FilterRmb(const PlayerInput& in) {
    PlayerInput out = in;
    if (!Interact_RmbConsumed()) return out;
    const uint64_t rmb = 1ull << (int)Act::Secondary;
    out.held &= ~rmb; out.pressed &= ~rmb; out.released &= ~rmb;
    return out;
}

void Machine_Update(Scene& world, const PlayerInput& in) {
    MachineInteract(world, in);
    // Latch where a held RMB began, after the aim is known: a press on empty space keeps its item working while held (a
    // shield, a telescope zoom swept over a chest), a press on an interactable stays the interactable's until released.
    if (!in.down(Act::Secondary)) rmbStart = RmbStart::None;
    else if (rmbStart == RmbStart::None) rmbStart = AimTakesRmb() ? RmbStart::Interactable : RmbStart::Empty;
}

namespace {
void MachineInteract(Scene& world, const PlayerInput& in) {
    machineAim = MachineAim{};
    farmAim = FarmAim{};
    Planetoid* body = player.nearest;
    if (!body || player.flyMode) return;
    Vector3 dir, ld; DVec3 lo;
    player.aimRay(world, *body, lo, ld, dir); // (F4, #1412: the eye in double with the frame path; the camera before)
    float best = 6.0f;
    if (body->plots) for (auto& pl : body->plots->plots) {
        Int3 c, f; float d;
        if (pl->raycast(lo, ld, best, c, f, d)) {
            const int part = pl->partAt(c);
            const PartType pt = part >= 0 ? pl->parts[part].type : PartType::Count;
            if (part >= 0 && (IsMachine(pt) || IsRotationPart(pt) || pt == PartType::HandPump || pt == PartType::Cable || pt == PartType::Bed || GateLocked(pl->parts[part]) || pt == PartType::Lever || pt == PartType::MbController)) { // vault doors and levers (#1641); cables: their network (#570); beds (#800); shafts (#1269)
                best = d;
                machineAim = { body->id, pl->id, part, pt, c };
            }
        }
        // Chimneys stick out of the cells (#285): aimable too, when nearer.
        const Vector3 po = pl->toPlot(lo);
        const Vector3 pd = { Vector3DotProduct(ld, pl->east), Vector3DotProduct(ld, pl->up), Vector3DotProduct(ld, pl->north) };
        int dp = -1; float dd = 0.0f;
        if (Vector3Length(po) < 200.0f && Machine_RaycastDetails(*pl, po, pd, best, dp, dd)) {
            best = dd;
            machineAim = { body->id, pl->id, dp, pl->parts[dp].type };
        }
    }
    // Farmland (#1085): tilled ground, or a crop on it, nearer than any machine. I plants seeds you carry there, O harvests
    // the crop there.
    TerrainHit th;
    if (RaycastTerrain(world, Ray{ world.camera.position, dir }, best * body->scale, th) && th.planetoid == body && Farm_Aims(*body, th.local)) {
        machineAim = MachineAim{};
        farmAim = { body->id, th.local, true };
        if (!in.enabled) return;
        if (in.hit(Act::Insert)) Farm_Plant(*body, th.local);
        if (in.hit(Act::Extract)) Farm_Harvest(*body, th.local);
        return;
    }
    if (machineAim.part < 0 || !in.enabled) return;
    if (machineAim.type == PartType::MbController) { // (#1429 PR 3) the wrench (RMB, outside build mode) forms the structure around a controller; it opens nothing
        if (buildTool.active || !in.hit(Act::Secondary)) return;
        Planetoid* mb = Planetoid::find(machineAim.body);
        Plot* mp = mb && mb->plots ? mb->plots->byId(machineAim.plot) : nullptr;
        if (!mp || machineAim.part >= (int)mp->parts.size() || mp->parts[machineAim.part].type != PartType::MbController) return;
        if (Hotbar_Slot(Hotbar_Selected()) != "wrench") { GameUI_Toast("Hold a wrench and right-click the controller to form the structure"); return; }
        if (gFormRequest) gFormRequest(*mb, *mp, machineAim.part); // joined: the server checks it and tells everyone
        else FormLocally(*mb, *mp, machineAim.part);
        return;
    }
    if (machineAim.type == PartType::Gate) { // (#1641) a vault door: RMB (outside build mode) or E opens it with a dungeon_key, which it takes
        if (!((!buildTool.active && in.hit(Act::Secondary)) || in.hit(Act::Use))) return;
        Planetoid* gb = Planetoid::find(machineAim.body);
        Plot* gp = gb && gb->plots ? gb->plots->byId(machineAim.plot) : nullptr;
        if (!gp || machineAim.part >= (int)gp->parts.size() || !GateLocked(gp->parts[machineAim.part])) return;
        if (inventory.count("dungeon_key") <= 0) { GameUI_Toast("Locked: it needs a dungeon key"); Audio_Play(Sfx::Error, 0.5f); return; }
        if (MachineAct_Remote()) { gMachineRequest(machineAim.body, machineAim.plot, machineAim.part, kMachineActUnlockVault, std::string(), 0, 0, 0); return; } // (#1684) joined: the host unlocks, takes the key and tells everyone
        inventory.remove("dungeon_key", 1);
        std::vector<PartEdit> edits;
        Dungeon_UnlockVault(*gp, machineAim.part, &edits);
        if (gNetRole == NetRole::Host) Host_BroadcastPartEdits(*gb, *gp, edits); // (#1684) a host's own player: the clients get it too
        GameUI_Toast("The vault door unlocks");
        Audio_Play(Sfx::Pickup, 0.4f, 0.7f);
        return;
    }
    if (machineAim.type == PartType::Lever) { // (#1641 slice 2c) RMB (outside build mode) or E pulls it: every portcullis within 12 m of it in its plot opens
        if (!((!buildTool.active && in.hit(Act::Secondary)) || in.hit(Act::Use))) return;
        Planetoid* lb = Planetoid::find(machineAim.body);
        Plot* lp = lb && lb->plots ? lb->plots->byId(machineAim.plot) : nullptr;
        if (!lp || machineAim.part >= (int)lp->parts.size() || lp->parts[machineAim.part].type != PartType::Lever) return;
        if (lp->parts[machineAim.part].material != 0) { GameUI_Toast("The lever is already down"); return; }
        if (MachineAct_Remote()) { gMachineRequest(machineAim.body, machineAim.plot, machineAim.part, kMachineActPullLever, std::string(), 0, 0, 0); return; } // (#1684) joined: the host pulls it and tells everyone
        std::vector<PartEdit> edits;
        const int gates = Dungeon_PullLever(*lp, machineAim.part, &edits);
        if (gNetRole == NetRole::Host) Host_BroadcastPartEdits(*lb, *lp, edits); // (#1684)
        GameUI_Toast(gates <= 0 ? "The lever clicks, and nothing happens" : "Somewhere, a gate grinds open");
        Audio_Play(Sfx::Pickup, 0.5f, 0.5f);
        return;
    }
    if (machineAim.type == PartType::Cable) return; // looked at only: its panel shows the network (#570)
    if (machineAim.type == PartType::Anvil && !buildTool.active && in.down(Act::Secondary)) { // (#1271 L2) hold RMB on an anvil to hammer
        static double nextStrike = 0.0;
        if (GetTime() >= nextStrike) { Factory::Strike(machineAim.body, machineAim.plot, machineAim.part); Audio_Play(Sfx::Build, 0.45f, 1.6f); nextStrike = GetTime() + 0.45; }
        return;
    }
    if (machineAim.type == PartType::HandCrank) { // (#1269) hold RMB (outside build mode) to turn it
        if (!buildTool.active && in.down(Act::Secondary)) { Factory::Crank(machineAim.body, machineAim.plot, machineAim.part); crankAt = GetTime(); }
        return;
    }
    if (IsRotationPart(machineAim.type) && !IsMachine(machineAim.type)) return; // shafts, cogs, gearboxes, gauges: looked at only
    if (machineAim.type == PartType::HandPump) {
        // Hand pump (#280): hold RMB (outside build mode) to work it; the water goes down its pipes at the next pass.
        if (!buildTool.active && in.down(Act::Secondary)) {
            Waterworks_HandPump(machineAim.body, machineAim.plot, machineAim.part, GameFrameTime());
            if (GetTime() - handPumpAt > 1.0) handPumpLitres = 0.0f; // a new go
            handPumpAt = GetTime(); handPumpLitres += kHandPumpRate * 1000.0f * GameFrameTime(); // the HUD's feedback (#571)
            static double nextSound = 0.0;
            if (GetTime() > nextSound) { Audio_Play(Sfx::Pickup, 0.25f, 0.6f); nextSound = GetTime() + 0.6; }
        }
        return;
    }
    if (machineAim.type == PartType::Campfire && in.hit(Act::Use)) {
        // (#1103) E on a campfire lights it or puts it out, as its screen's button does (E elsewhere opens the game window)
        const MachineView v = Factory::View(machineAim.body, machineAim.plot, machineAim.part);
        const bool burning = v.status == MachineStatus::Working;
        if (burning || v.fuelSlot.second > 0 || v.energy > 0) {
            MachineAct_SetFire(machineAim.body, machineAim.plot, machineAim.part, !burning);
            GameUI_Toast(burning ? "Campfire put out" : "Campfire lit");
        } else GameUI_Toast("The campfire needs fuel to light");
        return;
    }
    if (machineAim.type == PartType::Bed) {
        // A bed (#800): RMB (outside build mode) makes it your spawn; in a joined game the server sets it.
        if (!buildTool.active && in.hit(Act::Secondary)) {
            Planetoid* bb = Planetoid::find(machineAim.body);
            Plot* bp = bb && bb->plots ? bb->plots->byId(machineAim.plot) : nullptr;
            if (bb && bp) {
                if (gBedUseRequest) gBedUseRequest(*bb, *bp, machineAim.part);
                else Survival_SetBedFor(survival, *bb, *bp, machineAim.part, true);
            }
        }
        return;
    }
    // RMB (outside build mode, where RMB places parts) opens the machine's own screen: its slots
    // beside the inventory, like the player's.
    if (!buildTool.active && in.hit(Act::Secondary)) {
        // (#1254) phials in hand on an alchemy vessel: empty ones fill from it, full ones pour in
        if (Alchemy_UsePhials(machineAim.body, machineAim.plot, machineAim.part, machineAim.type, Hotbar_Slot(Hotbar_Selected()))) { Audio_Play(Sfx::Pickup, 0.4f, 1.4f); return; }
        GameUI_OpenContainer(machineAim.body, machineAim.plot, machineAim.part, (int)machineAim.type);
        Audio_Play(Sfx::Pickup, 0.3f, 1.0f);
        return;
    }
    if (in.hit(Act::Insert) && (machineAim.type == PartType::SteamLab || (machineAim.type >= PartType::LvLab && machineAim.type <= PartType::HvLab)) && ResearchGate_HandIn() > 0) return; // (#1271 L6b-2) milestone deliveries and finds, at a lab
    if (in.hit(Act::Insert)) {
        // Everything in the inventory the machine might take: fuel for generators, recipe inputs,
        // anything for a chest (the simulation hands back what does not fit).
        // A chest takes everything (#568) but tools and what is on the hotbar, as the chest screen's "Put in all".
        if (MachineAct_Remote()) gMachineRequest(machineAim.body, machineAim.plot, machineAim.part, 2, std::string(), 0, 0, 0); // (#332) the server puts in all it takes
        const auto all = MachineAct_Remote() ? decltype(inventory.all()){} : inventory.all(); // copy: we remove while walking it
        for (const auto& kv : all) {
            if (kv.second <= 0 || !Factory::Accepts(machineAim.type, kv.first)) continue;
            if (machineAim.type == PartType::Chest) {
                bool bound = false;
                for (int s = 0; s < Inventory::kHotbar; ++s) bound = bound || Hotbar_Slot(s) == kv.first;
                if (bound || KindOf(kv.first) == ItemKind::Tool) continue;
            }
            if (inventory.remove(kv.first, kv.second)) MachineAct_Insert(machineAim.body, machineAim.plot, machineAim.part, kv.first, kv.second);
        }
        Audio_Play(Sfx::Pickup, 0.4f, 0.8f);
    }
    if (in.hit(Act::Extract)) { MachineAct_Extract(machineAim.body, machineAim.plot, machineAim.part); Audio_Play(Sfx::Pickup, 0.4f, 1.2f); }
}
} // namespace

// A ", "-separated item list as lines at most maxW px wide (font size `size`), the first starting with `prefix`; past
// maxLines the last line ends with "+N more" (#567).
std::vector<std::string> WrapItemList(const char* prefix, const std::string& list, int maxW, int size, int maxLines) {
    const auto more = [&](std::string line, int rest) { // drops entries off the end until "+N more" fits
        for (;;) {
            const std::string s = line + TextFormat("  +%d more", rest);
            const size_t cut = line.rfind(", ");
            if (MeasureText(s.c_str(), size) <= maxW || cut == std::string::npos) return s;
            line.resize(cut); ++rest;
        }
    };
    std::vector<std::string> entries;
    for (size_t a = 0; a <= list.size();) {
        size_t b = list.find(", ", a);
        if (b == std::string::npos) b = list.size();
        entries.push_back(list.substr(a, b - a));
        a = b + 2;
    }
    std::vector<std::string> lines;
    std::string cur = prefix;
    bool fresh = true; // nothing on `cur` but the prefix or indent
    for (size_t i = 0; i < entries.size(); ++i) {
        const std::string next = fresh ? cur + entries[i] : cur + ", " + entries[i];
        if (fresh || MeasureText(next.c_str(), size) <= maxW) { cur = next; fresh = false; continue; }
        if ((int)lines.size() + 1 == maxLines) { lines.push_back(more(cur, (int)(entries.size() - i))); return lines; } // full
        lines.push_back(cur + ",");
        cur = "    " + entries[i];
    }
    lines.push_back(cur);
    return lines;
}

void Machine_DrawHUD() {
    if (GameUI_ContainerOpen()) return; // the open screen watches its own machine
    if (farmAim.on) { // (#1085) farmland ground or a crop
        if (Planetoid* fb = Planetoid::find(farmAim.body)) Farm_DrawAimHUD(*fb, farmAim.local);
        return;
    }
    if (machineAim.part < 0) return;
    if (machineAim.type == PartType::MbController) { // (#1429 PR 3)
        const Planetoid* cb = Planetoid::find(machineAim.body);
        const Plot* cp = cb && cb->plots ? cb->plots->byId(machineAim.plot) : nullptr;
        const MbDef* cd = cp && machineAim.part < (int)cp->parts.size() ? Multiblock_DefForController(cp->parts[machineAim.part]) : nullptr;
        const bool on = cp && std::binary_search(cp->formed.begin(), cp->formed.end(), machineAim.part);
        const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40;
        DrawRectangle(x - 8, y - 6, 400, 52, Fade(BLACK, 0.55f));
        DrawText(TextFormat("%s controller: %s", cd ? cd->name.c_str() : "multiblock", on ? "formed" : "not formed"), x, y, 20, SKYBLUE);
        DrawText(on ? "break a part of the structure to unform it" : "wrench + RMB to form the structure", x, y + 24, 16, RAYWHITE);
        return;
    }
    if (machineAim.type == PartType::Lever) { // (#1641 slice 2c)
        const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40;
        DrawRectangle(x - 8, y - 6, 300, 52, Fade(BLACK, 0.55f));
        DrawText("lever", x, y, 20, SKYBLUE);
        DrawText("RMB / E to pull it", x, y + 24, 16, RAYWHITE);
        return;
    }
    if (machineAim.type == PartType::Gate) { // (#1641) a vault door
        const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40;
        DrawRectangle(x - 8, y - 6, 340, 52, Fade(BLACK, 0.55f));
        DrawText("vault door: locked", x, y, 20, SKYBLUE);
        DrawText(inventory.count("dungeon_key") > 0 ? "RMB / E to unlock with your dungeon key" : "it needs a dungeon key", x, y + 24, 16, RAYWHITE);
        return;
    }
    if (machineAim.type == PartType::HandPump) {
        const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40;
        const bool pumping = GetTime() - handPumpAt < 0.15; // worked this frame or the last few (#571)
        DrawRectangle(x - 8, y - 6, 380, pumping ? 76 : 52, Fade(BLACK, 0.55f));
        DrawText(pumping ? "hand pump: PUMPING" : "hand pump", x, y, 20, pumping ? Color{ 120, 200, 255, 255 } : SKYBLUE);
        DrawText(TextFormat("hold RMB to pump %.0f L/s into its pipes", kHandPumpRate * 1000.0f), x, y + 24, 16, RAYWHITE);
        if (pumping) { // a bar that runs with each stroke (a stroke a second) and the litres sent this go
            const float stroke = (float)fmod(GetTime(), 1.0);
            DrawRectangle(x, y + 48, 160, 10, Fade(RAYWHITE, 0.2f));
            DrawRectangle(x, y + 48, (int)(160 * (stroke < 0.5f ? stroke * 2.0f : 2.0f - stroke * 2.0f)), 10, Color{ 90, 170, 255, 255 });
            DrawText(TextFormat("%.0f L this go", handPumpLitres), x + 172, y + 45, 16, RAYWHITE); // (none flows with no water under it)
        }
        return;
    }
    if (IsRotationPart(machineAim.type) && !IsMachine(machineAim.type)) { // (#1269) a shaft line: its speed and its network's load
        const RotationView rv = Factory::RotationAt(machineAim.body, machineAim.plot, machineAim.part);
        const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40;
        const bool crank = machineAim.type == PartType::HandCrank;
        DrawRectangle(x - 8, y - 6, 420, crank ? 100 : 76, Fade(BLACK, 0.55f));
        DrawText(PartName(machineAim.type), x, y, 20, SKYBLUE);
        if (!rv.valid) { DrawText("...", x, y + 24, 16, RAYWHITE); return; }
        const Color amber = { 255, 190, 70, 255 };
        if (rv.jam[0]) DrawText(TextFormat("jammed: %s", rv.jam), x, y + 24, 16, Color{ 255, 110, 90, 255 });
        else DrawText(TextFormat("%.1f RPM  (%s)", fabsf(rv.rpm), rv.status), x, y + 24, 16, RAYWHITE);
        if (rv.demand > 0 || rv.capacity > 0)
            DrawText(TextFormat("stress %d%% - running at %d%%  (%d / %d SU)", rv.stressPercent, rv.speedPercent, rv.demand, rv.capacity), x, y + 44, 16, rv.stressPercent > 100 ? amber : RAYWHITE);
        else DrawText("no load", x, y + 44, 16, RAYWHITE);
        if (crank) DrawText(GetTime() - crankAt < 0.2 ? "turning" : "hold RMB to turn it", x, y + 68, 16, GetTime() - crankAt < 0.2 ? Color{ 120, 200, 255, 255 } : RAYWHITE);
        return;
    }
    if (machineAim.type == PartType::Cable) { // its power network (#570): what feeds it and what flowed in the last tick
        const CableView cv = Factory::Cable(machineAim.body, machineAim.plot, machineAim.cell);
        const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40;
        DrawRectangle(x - 8, y - 6, 400, 76, Fade(BLACK, 0.55f));
        DrawText("cable", x, y, 20, SKYBLUE);
        if (!cv.valid) { DrawText("...", x, y + 24, 16, RAYWHITE); return; }
        if (cv.voltage <= 0) { DrawText("no power: nothing feeds this cable", x, y + 24, 16, Fade(RAYWHITE, 0.75f)); return; }
        DrawText(TextFormat("%d V network: %d generator%s, %d machine%s, %d batter%s", cv.voltage, cv.generators, cv.generators == 1 ? "" : "s",
                            cv.consumers, cv.consumers == 1 ? "" : "s", cv.batteries, cv.batteries == 1 ? "y" : "ies"), x, y + 24, 16, RAYWHITE);
        const float frac = cv.demand > 0.0f ? std::min(1.0f, cv.supply / cv.demand) : (cv.supply > 0.0f ? 1.0f : 0.0f);
        DrawRectangle(x, y + 48, 160, 10, Fade(RAYWHITE, 0.2f));
        DrawRectangle(x, y + 48, (int)(160 * frac), 10, cv.supply + 0.01f >= cv.demand ? Color{ 250, 210, 70, 255 } : ORANGE);
        DrawText(TextFormat("%.0f EU/t supplied, %.0f wanted", cv.supply, cv.demand), x + 172, y + 45, 16, RAYWHITE);
        return;
    }
    const MachineView v = Factory::View(machineAim.body, machineAim.plot, machineAim.part);
    const bool chest = machineAim.type == PartType::Chest;
    const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40, w = 420;
    // The item lists ("coal x3200, stone x40, ...", totals per item) wrap inside the panel; past kMaxLines the rest is
    // counted (#567: a full chest ran off the screen).
    const std::vector<std::string> in = WrapItemList(chest ? "holds: " : "in: ", v.input, w - 16, 16, 4);
    const std::vector<std::string> out = chest ? std::vector<std::string>{} : WrapItemList("out: ", v.output, w - 16, 16, 2);
    const std::string runLine = Alchemy_IsVessel(machineAim.type) ? Alchemy_PanelLine(machineAim.body, machineAim.plot, machineAim.part) : std::string(); // (#1254)
    const int lines = (int)(in.size() + out.size()) + (runLine.empty() ? 0 : 1);
    DrawRectangle(x - 8, y - 6, w, 72 + 20 * lines, Fade(BLACK, 0.55f));
    DrawText(TextFormat("%s: %s", PartName(machineAim.type), v.valid ? MachineStatusName(v.status) : "..."), x, y, 20, SKYBLUE);
    if (v.duration > 0) {
        DrawRectangle(x, y + 26, 300, 8, Fade(RAYWHITE, 0.2f));
        DrawRectangle(x, y + 26, 300 * v.progress / v.duration, 8, ORANGE);
    }
    int ty = y + 40;
    if (v.rotation) { // (#1269) the shaft it runs on
        DrawRectangle(x - 8, y + 66 + 20 * lines, w, 22, Fade(BLACK, 0.55f));
        DrawText(TextFormat("%.1f RPM (%s), stress %d%% - running at %d%%", v.rpm, v.rotStatus, v.stressPercent, v.speedPercent), x, y + 68 + 20 * lines, 16,
                 v.stressPercent > 100 ? Color{ 255, 190, 70, 255 } : RAYWHITE);
    }
    { // (#1271 L6b) recipes of this machine the research has not unlocked, and what unlocks them
        static PartType noteFor = PartType::Count; static std::string note; static double noteAt = -10.0;
        if (noteFor != machineAim.type || GetTime() - noteAt > 1.0) { note = ResearchGate_MachineLockNote(machineAim.type); noteFor = machineAim.type; noteAt = GetTime(); }
        if (!note.empty()) {
            const int ny = y + 66 + 20 * lines + (v.rotation ? 22 : 0);
            DrawRectangle(x - 8, ny - 2, w, 22, Fade(BLACK, 0.55f));
            DrawText(note.c_str(), x, ny, 16, Color{ 255, 190, 70, 255 });
        }
    }
    for (const std::string& l : in) { DrawText(l.c_str(), x, ty, 16, RAYWHITE); ty += 20; }
    for (const std::string& l : out) { DrawText(l.c_str(), x, ty, 16, RAYWHITE); ty += 20; }
    if (!runLine.empty()) { DrawText(runLine.c_str(), x, ty, 16, runLine.find("refused") != std::string::npos ? Color{ 255, 150, 110, 255 } : Color{ 205, 190, 255, 255 }); ty += 20; }
    if (v.energyMax > 0) DrawText(TextFormat("charge %d / %d EU", v.energy, v.energyMax), x, ty, 16, YELLOW);
    else DrawText(machineAim.type == PartType::Anvil ? "hold RMB to hammer   I put work in   O take it out" : chest ? "RMB open   I put items in   O take all" : "RMB open   I put items in   O take the output", x, ty, 16, Fade(RAYWHITE, 0.7f));
}

// The world is going (#902): the aim points into its bodies and plots.
void BuildTool_ForgetWorld() { aim = {}; }

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"
// #1578 (playtest: a ladder placed against a wall showed its edge): automatic facing turns a ladder's climbing face out of the
// wall (its frame against the wall), and a machine's front (local +z) as before.
TEST(BuildFacing, LadderClimbingFaceOutOfTheWallMachinesFrontAsBefore) {
    for (const Int3 out : { Int3{ 1, 0, 0 }, Int3{ -1, 0, 0 }, Int3{ 0, 0, 1 }, Int3{ 0, 0, -1 } }) {
        Part ladder{ PartType::Ladder, { 0, 0, 0 }, (uint8_t)BuildTool_AutoOrient((int)PartType::Ladder, out), 0 };
        std::vector<std::pair<Int3, CellShape>> fp;
        Plot::Footprint(ladder, fp);
        ASSERT_EQ(fp.size(), 1u);
        const Vector3 n = fp[0].second.n[0]; // the frame's open side: the cell beyond its climbing face
        EXPECT_EQ(Int3({ (int)lroundf(n.x), (int)lroundf(n.y), (int)lroundf(n.z) }), out) << "ladder, out of the wall (" << out.x << ", " << out.z << ")";
        const int* m = Orientation(BuildTool_AutoOrient((int)PartType::Chest, out));
        EXPECT_EQ(Int3({ m[2], m[5], m[8] }), out) << "a chest's front, out of the wall (" << out.x << ", " << out.z << ")";
    }
}
// #1684: a vault unlock and a lever pull are one shared edit; the part changes it records, replayed on another plot in order
// (as a client applies the host's MsgPlotEdit), leave that plot with the same parts (indices included).
static void ReplayPartEdits(Plot& plot, const std::vector<PartEdit>& edits) {
    for (const PartEdit& e : edits) {
        if (e.remove) plot.remove(e.index);
        else EXPECT_EQ(plot.place(e.part), e.index);
    }
}
static void ExpectSameParts(const Plot& a, const Plot& b) {
    ASSERT_EQ(a.parts.size(), b.parts.size());
    for (size_t i = 0; i < a.parts.size(); ++i) {
        EXPECT_EQ((int)a.parts[i].type, (int)b.parts[i].type) << "part " << i;
        EXPECT_EQ(a.parts[i].material, b.parts[i].material) << "part " << i;
        EXPECT_EQ(a.parts[i].origin, b.parts[i].origin) << "part " << i;
    }
}
TEST(DungeonEdits, VaultUnlockOpensBothHalvesAndReplaysOnAClient) {
    Plot host, client;
    for (Plot* p : { &host, &client }) {
        p->place(Part{ PartType::Gate, { 4, 0, 4 }, 0, 1 });
        p->place(Part{ PartType::Gate, { 4, 2, 4 }, 0, 1 });
        p->place(Part{ PartType::Gate, { 20, 0, 4 }, 0, 1 }); // another vault door: stays locked
    }
    EXPECT_EQ(Dungeon_UnlockVault(host, 2 + 99), 0);          // no such part
    std::vector<PartEdit> edits;
    EXPECT_EQ(Dungeon_UnlockVault(host, 1, &edits), 2);        // the aimed half opens its other half too
    EXPECT_EQ(host.parts[3].material, 0);
    EXPECT_EQ(host.parts[4].material, 0);
    EXPECT_TRUE(GateLocked(host.parts[2]));
    EXPECT_EQ(Dungeon_UnlockVault(host, 3), 0);                // already open: not a vault door any more
    ReplayPartEdits(client, edits);
    ExpectSameParts(host, client);
}
TEST(DungeonEdits, LeverOpensPortcullisesWithin12mOnly) {
    Plot host, client;
    for (Plot* p : { &host, &client }) {
        p->place(Part{ PartType::Lever, { 0, 0, 0 }, 0, 0 });
        p->place(Part{ PartType::Gate, { 10, 0, 0 }, 0, 2 });  // 5 m: opens
        p->place(Part{ PartType::Gate, { 40, 0, 0 }, 0, 2 });  // 20 m: stays shut
        p->place(Part{ PartType::Gate, { 0, 0, 10 }, 0, 1 });  // a vault door: not a portcullis
    }
    EXPECT_EQ(Dungeon_PullLever(host, 3), -1);                 // not a lever
    std::vector<PartEdit> edits;
    EXPECT_EQ(Dungeon_PullLever(host, 0, &edits), 1);
    EXPECT_EQ(host.parts[4].material, 1);                      // the lever, down
    EXPECT_EQ(host.parts[5].material, 0);                      // the near portcullis, open
    EXPECT_TRUE(GatePortcullis(host.parts[2]));                // the far one, shut
    EXPECT_TRUE(GateLocked(host.parts[3]));
    EXPECT_EQ(Dungeon_PullLever(host, 4), -1);                 // already down
    ReplayPartEdits(client, edits);
    ExpectSameParts(host, client);
}
#endif
