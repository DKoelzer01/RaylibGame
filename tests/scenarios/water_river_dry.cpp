// water/river-dry (#1438): rivers that show as dry sand. UI/UX's river-valley slides found a wide sandy channel, with a thin
// crease along its middle and no water drawn, where waterSourceAt reports river water. This samples every 10 m within
// 2 km of the start: at each river point (kind 1) it compares the river's level with the ground, both from the height map
// and from the real (density) surface found along the vertical, and with the generated water the fluid model would give the
// column (a cell is filled only where the level stands over the ground). Reported: river points, how many stand dry (level
// at or under the ground) by each measure, how deep the dry ones are under it, and the nearest dry spots (lat/lon, to look
// at). Reported only, not gated (a diagnosis).
#include "scenario.h"
#include "rivers.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
constexpr float kSpacing = 10.0f, kReach = 2000.0f;

// The real ground under a direction: the highest solid point within +-12 m of the height map (density, all edits).
float SurfaceAt(ScenarioCtx& c, Vector3 dir, float hm) {
    for (float h = hm + 12.0f; h > hm - 12.0f; h -= 0.1f)
        if (c.body->densityAt(Vector3Scale(dir, c.field->radius + h)) >= ISO_LEVEL) return h;
    return hm - 12.0f;
}

// The channels (pieces) whose half width holds a point, and the lowest and highest of their levels there (the carve's
// cross-section uses every one of them; nearest() reports the lowest level of those).
struct Holding { int n = 0; float lo = 1e30f, hi = -1e30f; };
Holding ChannelsHolding(const RiverNetwork& net, Vector3 dir, float R) {
    Holding h;
    const Vector3 p = Vector3Scale(dir, R);
    net.forEachLink([&](const RiverNetwork::Link& l) {
        const Vector3 ab = Vector3Subtract(l.b, l.a);
        const float len2 = Vector3DotProduct(ab, ab);
        const float t = len2 > 0.0f ? std::clamp(Vector3DotProduct(Vector3Subtract(p, l.a), ab) / len2, 0.0f, 1.0f) : 0.0f;
        const float d = Vector3Distance(p, Vector3Add(l.a, Vector3Scale(ab, t)));
        const float half = 0.5f * (l.wa + (l.wb - l.wa) * t); // Link widths are full widths
        if (d >= half) return;
        const float lv = l.la + (l.lb - l.la) * t;            // (ignores a waterfall's step: a diagnosis)
        ++h.n; h.lo = std::min(h.lo, lv); h.hi = std::max(h.hi, lv);
    });
    return h;
}
}

SCENARIO("water/river-dry") {
    s.tags = { "water", "terrain", "local" }; // local: a diagnosis, ~2 min of field queries (no runner value)
    s.world = ScenarioWorld::Planet;
    s.timeout = 300.0;
    s.record = false;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        return c.T > 1.0 && WorldIdle();
    });
    s.step("sample the rivers around the start", 200.0, [](ScenarioCtx& c) {
        const Vector3 home = c.spawnDirLocal();
        const Vector3 ax = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), ay = Vector3CrossProduct(ax, home);
        const float R = c.field->radius;
        int river = 0, dryMap = 0, drySurf = 0, thin = 0;
        int dryBand[3] = {}; // dry over the real ground, by where across the channel: the centre (d < 0.5 half), the middle, the edge (> 0.85)
        double sumUnder = 0.0;
        const RiverNetwork* net = c.field->rivers();
        // Each piece's in-lake flag (forEachLink walks the pieces in index order, as Hit::piece counts them): carve() skips
        // in-lake pieces but nearest() doesn't, so a point by one outside the lake's water gets river water and no channel.
        std::vector<char> pieceInLake;
        if (net) net->forEachLink([&](const RiverNetwork::Link& l) { pieceInLake.push_back(l.inLake ? 1 : 0); });
        int dryChannelInLake = 0; // dry at the centre or middle, nearest piece in-lake
        int dryChannelTwo = 0, dryChannelLow = 0; // ...held by two channels or more; ...on ground under 0.5 m over the sea
        struct Dry { float dist, under, lat, lon, width, across, half; bool inLake; };
        std::vector<Dry> dry;
        for (float y = -kReach; y <= kReach; y += kSpacing)
            for (float x = -kReach; x <= kReach; x += kSpacing) {
                const Vector3 dir = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(ax, x / R), Vector3Scale(ay, y / R))));
                const PlanetField::WaterSource ws = c.field->waterSourceAt(dir);
                if (ws.kind != 1) continue;
                ++river;
                const float hm = c.field->elevation(dir);
                const float depthMap = ws.level - hm;              // water over the height map (m)
                if (depthMap <= 0.0f) ++dryMap;
                const float surf = SurfaceAt(c, dir, hm);
                const float depth = ws.level - surf;               // water over the real ground (m)
                if (depth > 0.0f && depth < 0.25f) ++thin;          // under one 0.25 m cell: drawn as a film at most
                if (depth <= 0.0f) {
                    ++drySurf; sumUnder += -depth;
                    // The dry point's channel width: river points along x on this row (a rough measure).
                    float w = 0.0f;
                    for (float dx = kSpacing; dx < 200.0f; dx += kSpacing) {
                        const Vector3 d2 = Vector3Normalize(Vector3Add(dir, Vector3Scale(ax, dx / R)));
                        if (c.field->waterSourceAt(d2).kind != 1) break;
                        w += kSpacing;
                    }
                    // Where across the channel (the carve's bed stands over the level in its outer ~12% by design: a bank).
                    const RiverNetwork::Hit rh = net ? net->nearest(dir, 60.0f) : RiverNetwork::Hit{};
                    const float across = rh.found && rh.halfWidth > 0.0f ? rh.dist / rh.halfWidth : -1.0f;
                    ++dryBand[across < 0.5f ? 0 : across <= 0.85f ? 1 : 2];
                    const bool inLake = rh.piece >= 0 && rh.piece < (int)pieceInLake.size() && pieceInLake[rh.piece];
                    if (inLake && across <= 0.85f) ++dryChannelInLake;
                    if (across <= 0.85f && net) { // a dry channel: two channels at once, or low ground (the carve's sea fade)?
                        const Holding hold = ChannelsHolding(*net, dir, R);
                        const float base = c.field->elevationBase(dir);
                        if (hold.n >= 2) ++dryChannelTwo;
                        if (base < 0.5f) ++dryChannelLow;
                        c.note(TextFormat("dry channel point %.0f m out: %d channels hold it (levels %.2f-%.2f m, reported %.2f), ground before rivers %.2f m over the sea, after %.2f",
                                          sqrtf(x * x + y * y), hold.n, hold.n ? hold.lo : 0.0f, hold.n ? hold.hi : 0.0f, ws.level, base, hm));
                    }
                    dry.push_back({ sqrtf(x * x + y * y), -depth, asinf(dir.y) * RAD2DEG, atan2f(dir.x, dir.z) * RAD2DEG, w, across, rh.halfWidth, inLake });
                }
            }
        // The centre and middle first (a dry channel), the edges last (banks by design); nearest first within each.
        std::sort(dry.begin(), dry.end(), [](const Dry& a, const Dry& b) {
            const bool ea = a.across > 0.85f, eb = b.across > 0.85f;
            return ea != eb ? !ea : a.dist < b.dist; });
        c.note(TextFormat("dry points across the channel: %d at its centre (< 0.5 of the half-width), %d in the middle, %d at the edge (> 0.85: the carve's bank)",
                          dryBand[0], dryBand[1], dryBand[2]));
        c.metric("dry at the channel's centre", dryBand[0], "points");
        c.metric("dry in the channel's middle", dryBand[1], "points");
        c.note(TextFormat("of the %d dry at the centre or middle, %d lie by an in-lake piece (carve() skips those; nearest() doesn't)",
                          dryBand[0] + dryBand[1], dryChannelInLake));
        c.metric("dry in the channel by an in-lake piece", dryChannelInLake, "points");
        c.note(TextFormat("of those, %d lie inside two channels or more, %d on ground under 0.5 m over the sea (the carve fades out there)",
                          dryChannelTwo, dryChannelLow));
        c.check("sampled river points near the start", river > 0);
        c.note(TextFormat("%d river points within %.0f m (every %.0f m): %d stand dry over the height map, %d over the real ground (mean %.2f m under it); %d hold under one 0.25 m cell of water",
                          river, kReach, kSpacing, dryMap, drySurf, drySurf ? sumUnder / drySurf : 0.0, thin));
        for (size_t k = 0; k < dry.size() && k < 8; ++k)
            c.note(TextFormat("dry river point %.0f m from the start: the level %.2f m under the ground, %.2f of the half-width (%.1f m) from the centre line%s, lat %.4f lon %.4f, river %.0f m wide along x",
                              dry[k].dist, dry[k].under, dry[k].across, dry[k].half, dry[k].inLake ? " (an in-lake piece)" : "", dry[k].lat, dry[k].lon, dry[k].width));
        c.metric("river points", river, "", false);
        c.metric("dry over the real ground", drySurf, "points");
        c.metric("dry over the height map", dryMap, "points");
        c.metric("under one cell of water", thin, "points");
        return true;
    });
}
