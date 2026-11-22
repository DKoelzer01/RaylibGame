#include "locate.h"
#include "ambience.h"
#include "structures.h"
#include "voxel.h"
#include "world.h"

#include <raymath.h>
#include <algorithm>
#include <cctype>
#include <cmath>

namespace {
constexpr float kDeg = 57.29578f;

std::string Lower(std::string s) { for (char& c : s) c = (char)tolower((unsigned char)c); return s; }

Vector3 Ground(const PlanetField& pf, Vector3 dir) { return Vector3Scale(dir, pf.radius + pf.elevation(dir)); }

// The tangent frame at a direction: north toward the body's +Y pole, east = up x north (longitude grows eastward).
void Frame(Vector3 up, Vector3& north, Vector3& east) {
    north = Vector3Subtract({ 0, 1, 0 }, Vector3Scale(up, up.y));
    if (Vector3Length(north) < 1e-4f) north = Vector3Subtract({ 0, 0, -1 }, Vector3Scale(up, up.z)); // at a pole
    north = Vector3Normalize(north);
    east = Vector3CrossProduct(up, north);
}

// Candidates on a square grid round `home` (metres apart, out to `range`) where `ok` holds, nearest first; points closer
// than `apart` to a nearer one are the same feature.
template <class F> std::vector<Vector3> GridSearch(const PlanetField& pf, Vector3 home, float step, float range, float apart, int n, F ok) {
    Vector3 north, east; Frame(home, north, east);
    struct C { Vector3 d; float dd; };
    std::vector<C> found;
    for (float a = -range; a <= range; a += step)
        for (float b = -range; b <= range; b += step) {
            const float dd = a * a + b * b;
            if (dd > range * range) continue;
            const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, pf.radius), Vector3Add(Vector3Scale(east, a), Vector3Scale(north, b))));
            if (ok(d)) found.push_back({ d, dd });
        }
    std::sort(found.begin(), found.end(), [](const C& x, const C& y) { return x.dd < y.dd; });
    std::vector<Vector3> out;
    for (const C& c : found) {
        bool same = false;
        for (const Vector3& o : out) if (Vector3Distance(o, c.d) * pf.radius < apart) { same = true; break; }
        if (same) continue;
        out.push_back(c.d);
        if ((int)out.size() >= n) break;
    }
    return out;
}

// Volcano summits (#593, as --test-terrain's count): where the designed relief passes 900 m, climbed to its top.
std::vector<Vector3> Summits(const PlanetField& pf, Vector3 home, int n) {
    std::vector<Vector3> peaks;
    const std::vector<Vector3> slopes = GridSearch(pf, home, 1000.0f, 80000.0f, 3000.0f, 4 * n, [&](Vector3 d) { return pf.volcano(Vector3Scale(d, pf.radius)) > 900.0f; });
    for (const Vector3& d0 : slopes) {
        Vector3 b = d0; float bv = pf.volcano(Vector3Scale(b, pf.radius));
        Vector3 north, east; Frame(b, north, east);
        for (int step = 0; step < 60; ++step) {
            bool moved = false;
            for (int k = 0; k < 8; ++k) {
                const float a = k * 0.785398f;
                const Vector3 q = Vector3Normalize(Vector3Add(b, Vector3Scale(Vector3Add(Vector3Scale(east, cosf(a)), Vector3Scale(north, sinf(a))), 100.0f / pf.radius)));
                const float qv = pf.volcano(Vector3Scale(q, pf.radius));
                if (qv > bv) { bv = qv; b = q; moved = true; }
            }
            if (!moved) break;
        }
        bool dup = false;
        for (const Vector3& q : peaks) if (Vector3Distance(q, b) * pf.radius < 3000.0f) { dup = true; break; }
        if (!dup) peaks.push_back(b);
    }
    std::sort(peaks.begin(), peaks.end(), [&](const Vector3& x, const Vector3& y) { return Vector3DotProduct(x, home) > Vector3DotProduct(y, home); });
    if ((int)peaks.size() > n) peaks.resize(n);
    return peaks;
}
}

Vector3 Locate_DirFromLatLon(float latDeg, float lonDeg) {
    const float la = latDeg / kDeg, lo = lonDeg / kDeg;
    return { cosf(la) * cosf(lo), sinf(la), cosf(la) * sinf(lo) };
}
void Locate_LatLonOf(Vector3 local, float& latDeg, float& lonDeg) {
    const Vector3 d = Vector3Normalize(local);
    latDeg = asinf(std::clamp(d.y, -1.0f, 1.0f)) * kDeg;
    lonDeg = atan2f(d.z, d.x) * kDeg;
}

const char* Locate_CompassPoint(float bearing) {
    static const char* kPoints[] = { "N", "NE", "E", "SE", "S", "SW", "W", "NW" };
    const int k = (int)floorf(fmodf(bearing + 22.5f + 360.0f, 360.0f) / 45.0f);
    return kPoints[std::clamp(k, 0, 7)];
}

std::string Locate_FeatureList() {
    return "cave, hillside, sinkhole, resurgence, chamber, karst, tube, volcano, ruin, dungeon, camp, wreck, wonder, spring, oilseep, lake, river, ore, ore:<kind>, biome:<name>";
}

bool Locate_Find(const Planetoid& body, Vector3 from, const std::string& featureIn, int n, std::vector<LocateHit>& out, std::string& why) {
    out.clear();
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    if (!pf) { why = "Nothing to find here (not a planet)"; return false; }
    n = std::clamp(n, 1, 10);
    const std::string feature = Lower(featureIn);
    const Vector3 home = Vector3Normalize(from);
    auto add = [&](Vector3 at, const std::string& what) { out.push_back({ at, what }); };
    if (feature == "cave" || feature == "hillside" || feature == "sinkhole" || feature == "resurgence" || feature == "chamber" || feature == "karst" || feature == "tube") {
        if (!PlanetField::CavesGraph()) { why = "The cave graph is off (RAYLIBGAME_CAVES=0)"; return false; }
        // nearestCaveMouth's kinds (the cave finders, #615/#1052/#1181): -1 any entrance, 0 hillside mouth, 1 shaft (sinkholes,
        // karst sinks, lava-tube skylights), 2 resurgence, 3 chamber centre. karst and tube look near their country first.
        Vector3 near = home; float range = 3000.0f;
        const int kind = feature == "cave" ? -1 : feature == "hillside" ? 0 : feature == "resurgence" ? 2 : feature == "chamber" ? 3 : 1;
        if (feature == "karst") { // (#1052) karst sinks near strong karst ground (regions ~10 km apart)
            const auto g = GridSearch(*pf, home, 500.0f, 30000.0f, 0.0f, 1, [&](Vector3 d) { return pf->karstWeight(Vector3Scale(d, pf->radius)) > 0.6f && pf->elevation(d) > 10.0f; });
            if (g.empty()) { why = "No karst country within 30 km"; return false; }
            near = g[0]; range = 1500.0f;
        } else if (feature == "tube") { // (#1063) lava tube skylights on volcano flanks (volcanoes ~20 km apart)
            const auto g = GridSearch(*pf, home, 1000.0f, 60000.0f, 0.0f, 1, [&](Vector3 d) { const float w = pf->volcanoWeight(Vector3Scale(d, pf->radius)); return w > 0.4f && w < 0.8f && pf->elevation(d) > 10.0f; });
            if (g.empty()) { why = "No volcano flank within 60 km"; return false; }
            near = g[0]; range = 3000.0f;
        }
        Vector3 at{}, outward{};
        if (!pf->nearestCaveMouth(near, range, at, outward, kind)) { why = "None within " + std::to_string((int)(range / 1000.0f)) + " km of there"; return false; }
        add(at, feature == "cave" ? "cave entrance" : feature == "hillside" ? "hillside cave mouth" : feature == "sinkhole" ? "sinkhole" : feature == "resurgence" ? "resurgence (a spring from a cave)"
              : feature == "chamber" ? "cave chamber" : feature == "karst" ? "karst sink" : "lava tube skylight");
        if (n > 1) why = "(caves: the nearest only)";
    } else if (feature == "volcano") {
        for (const Vector3& d : Summits(*pf, home, n)) add(Ground(*pf, d), "volcano summit");
    } else if (feature == "ruin" || feature == "dungeon" || feature == "camp" || feature == "wreck" || feature == "wonder" || feature == "site") {
        std::vector<SiteInfo> sites;
        Structures_SitesNear(body, from, 12000.0f, sites);
        for (const SiteInfo& s : sites) {
            if ((int)out.size() >= n) break;
            if (feature == "dungeon" ? !s.dungeon : feature != "site" && s.kind != feature) continue;
            add(s.at, s.kind + ": " + s.set + (s.dungeon ? " (a dungeon under it)" : ""));
        }
    } else if (feature == "spring" || feature == "oilseep" || feature == "seep") {
        std::vector<Vector3> found;
        if (feature == "spring") Ambience_SpringsNear(body, from, 15000.0f, found); else Ambience_SeepsNear(body, from, 15000.0f, found);
        std::sort(found.begin(), found.end(), [&](const Vector3& a, const Vector3& b) { return Vector3Distance(a, from) < Vector3Distance(b, from); });
        for (const Vector3& p : found) { if ((int)out.size() >= n) break; add(p, feature == "spring" ? "water spring" : "oil seep"); }
    } else if (feature == "lake" || feature == "river") {
        const int kind = feature == "river" ? 1 : 2;
        const float step = kind == 1 ? 100.0f : 250.0f, range = kind == 1 ? 6000.0f : 15000.0f;
        for (const Vector3& d : GridSearch(*pf, home, step, range, 1500.0f, n, [&](Vector3 d) { return pf->waterSourceAt(d).kind == kind; }))
            add(Ground(*pf, d), kind == 1 ? "river" : "lake");
        if (kind == 2) { // and the bake's ponds (small lakes)
            for (const Vector3& d : pf->bakePondsByDistance(home)) {
                if (Vector3Distance(d, home) * pf->radius > range) break;
                add(Ground(*pf, d), "pond");
            }
        }
    } else if (feature == "ore" || feature.rfind("ore:", 0) == 0) {
        // Ore deposits (#207: one blob per 96 m cell in about a third of them): the cells under the ground round the
        // player, 5-100 m down; a deposit counts when its centre is under the surface. ore:<kind> = iron, copper, tin,
        // coal, bauxite, ilmenite.
        int want = -1;
        if (feature.size() > 4) {
            const std::string name = feature.substr(4);
            for (int k = 0; k < PlanetField::kOreKindCount; ++k) {
                const std::string id = PlanetField::OreKindName(k);
                if (id == name || id == name + "_ore") want = k;
            }
            if (want < 0) {
                why = "Unknown ore '" + name + "': ";
                for (int k = 0; k < PlanetField::kOreKindCount; ++k) {
                    std::string id = PlanetField::OreKindName(k);
                    if (id.size() > 4 && id.compare(id.size() - 4, 4, "_ore") == 0) id.resize(id.size() - 4);
                    why += (k ? ", " : "") + id;
                }
                return false;
            }
        }
        std::vector<Vector3> centres;
        GridSearch(*pf, home, 96.0f, 1500.0f, 0.0f, 1, [&](Vector3 d) {
            const float ground = pf->radius + pf->elevation(d);
            for (float depth = 5.0f; depth <= 100.0f; depth += 47.5f) {
                Vector3 c; float r; int kind;
                if (!pf->depositAt(Vector3Scale(d, ground - depth), c, r, kind) || (want >= 0 && kind != want)) continue;
                const float down = pf->radius + pf->elevation(Vector3Normalize(c)) - Vector3Length(c);
                if (down < 2.0f) continue; // (in the air, or barely under: no ore there)
                bool seen = false;
                for (const Vector3& o : centres) if (Vector3Distance(o, c) < 1.0f) { seen = true; break; }
                if (seen) continue;
                centres.push_back(c);
                char what[96];
                snprintf(what, sizeof what, "%s deposit, ~%.0f m across, %.0f m down", PlanetField::OreKindName(kind), 2.0f * r, down);
                add(Ground(*pf, Vector3Normalize(c)), what);
            }
            return false; // (the deposits are collected above; the grid only walks the ground)
        });
    } else if (feature.rfind("biome:", 0) == 0 || feature.rfind("biome ", 0) == 0) {
        const std::string name = feature.substr(6);
        int b = -1;
        for (int k = 0; k < (int)Biome::Count; ++k) if (Lower(BiomeName((Biome)k)) == name) b = k;
        if (b < 0) {
            why = "Unknown biome '" + name + "': ";
            for (int k = 0; k < (int)Biome::Count; ++k) why += (k ? ", " : "") + Lower(BiomeName((Biome)k));
            return false;
        }
        const Vector3 d = pf->nearestBiomeDir(home, (Biome)b);
        if (Vector3Distance(d, home) < 1e-6f && pf->biomeAt(home, pf->elevation(home)) != (Biome)b) { why = "No " + name + " found"; return false; }
        add(Ground(*pf, Vector3Normalize(d)), "biome: " + name);
    } else {
        why = "Unknown feature '" + featureIn + "'. Features: " + Locate_FeatureList();
        return false;
    }
    if (out.empty()) { if (why.empty()) why = "No " + feature + " found near here"; return false; }
    // Distance and bearing from where the player stands; nearest first (ponds and lakes come from two searches).
    Vector3 north, east; Frame(home, north, east);
    for (LocateHit& h : out) {
        h.distance = Vector3Distance(h.at, from);
        const Vector3 t = Vector3Subtract(h.at, Vector3Scale(home, Vector3DotProduct(h.at, home))); // across the ground
        h.bearing = fmodf(atan2f(Vector3DotProduct(t, east), Vector3DotProduct(t, north)) * kDeg + 360.0f, 360.0f);
    }
    std::sort(out.begin(), out.end(), [](const LocateHit& a, const LocateHit& b) { return a.distance < b.distance; });
    if ((int)out.size() > n) out.resize(n);
    return true;
}

std::string Locate_LocateCommand(const Planetoid& body, Vector3 from, const std::vector<std::string>& args, LocateHit* nearest) {
    if (args.empty()) return "/locate <feature> [n]: " + Locate_FeatureList();
    const int n = args.size() > 1 ? std::max(1, atoi(args[1].c_str())) : 1;
    std::vector<LocateHit> hits;
    std::string why;
    if (!Locate_Find(body, from, args[0], n, hits, why)) return why;
    std::string reply;
    for (size_t i = 0; i < hits.size(); ++i) {
        float la, lo; Locate_LatLonOf(hits[i].at, la, lo);
        char dist[32];
        if (hits[i].distance < 1000.0f) snprintf(dist, sizeof dist, "%.0f m", hits[i].distance); else snprintf(dist, sizeof dist, "%.1f km", hits[i].distance / 1000.0f);
        char line[256];
        snprintf(line, sizeof line, "%s%s: %s %s (%.0f deg), lat %.4f lon %.4f", i ? "\n" : "", hits[i].what.c_str(), dist, Locate_CompassPoint(hits[i].bearing), hits[i].bearing, la, lo);
        reply += line;
    }
    if (nearest) *nearest = hits[0];
    if (!why.empty()) reply += "\n" + why;
    return reply;
}

bool Locate_TpTarget(const Planetoid& body, Vector3 from, const std::vector<std::string>& args, Vector3& feet, bool& hover, std::string& reply) {
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    if (!pf) { reply = "Nowhere to go here (not a planet)"; return false; }
    Vector3 to{}; float height = 0.0f; std::string where;
    char* end = nullptr;
    const float lat = args.empty() ? 0.0f : strtof(args[0].c_str(), &end);
    if (args.size() >= 2 && end && *end == 0) { // lat lon [height]
        if (fabsf(lat) > 90.0f) { reply = "Latitude is -90 to 90"; return false; }
        const float lon = (float)atof(args[1].c_str());
        height = args.size() > 2 ? std::max(0.0f, (float)atof(args[2].c_str())) : 0.0f;
        const Vector3 d = Locate_DirFromLatLon(lat, lon);
        to = Vector3Scale(d, pf->radius + pf->elevation(d));
        char w[64]; snprintf(w, sizeof w, "lat %.4f lon %.4f", lat, lon); where = w;
    } else if (!args.empty()) { // a feature: its nearest
        std::vector<LocateHit> hits; std::string why;
        if (!Locate_Find(body, from, args[0], 1, hits, why)) { reply = why; return false; }
        to = hits[0].at; where = hits[0].what;
        height = args.size() > 1 ? std::max(0.0f, (float)atof(args[1].c_str())) : 0.0f;
    } else { reply = "/tp <lat> <lon> [height] or /tp <feature> [height] (" + Locate_FeatureList() + ")"; return false; }
    // On the ground (a cave feature: at it, underground) 0.3 m up so the feet don't start in the soil, or height m over it.
    feet = Vector3Add(to, Vector3Scale(Vector3Normalize(to), 0.3f + height));
    hover = height > 2.0f; // (creative: hover there rather than fall)
    char h[32]; snprintf(h, sizeof h, ", %.0f m up", height);
    reply = "Teleported to " + where + (height > 0.0f ? std::string(h) : std::string());
    return true;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1171): latitude and longitude round-trip; compass points ----
#include <gtest/gtest.h>

TEST(Locate, LatLonRoundTrip) {
    for (float la : { -80.0f, -12.5f, 0.0f, 33.3f, 89.0f })
        for (float lo : { -179.0f, -45.0f, 0.0f, 90.0f, 170.0f }) {
            float la2, lo2; Locate_LatLonOf(Vector3Scale(Locate_DirFromLatLon(la, lo), 60000.0f), la2, lo2);
            EXPECT_NEAR(la2, la, 1e-3f); EXPECT_NEAR(lo2, lo, 1e-3f);
        }
    const Vector3 up = Locate_DirFromLatLon(0.0f, 0.0f);
    EXPECT_NEAR(up.x, 1.0f, 1e-6f) << "lat 0 lon 0 is +X (the tests' convention, RAYLIBGAME_MOUTH_AT)";
}

TEST(Locate, CompassPoints) {
    EXPECT_STREQ(Locate_CompassPoint(0.0f), "N"); EXPECT_STREQ(Locate_CompassPoint(359.0f), "N");
    EXPECT_STREQ(Locate_CompassPoint(44.0f), "NE"); EXPECT_STREQ(Locate_CompassPoint(90.0f), "E");
    EXPECT_STREQ(Locate_CompassPoint(181.0f), "S"); EXPECT_STREQ(Locate_CompassPoint(270.0f), "W");
}
#endif
