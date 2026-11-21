#include "worldprofile.h"
#include "voxel.h"
#include "utils/logger.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <mutex>
#include <sstream>
#include <unordered_map>

bool WorldProfile::allowsCreature(const std::string& id) const {
    return allCreatures || std::find(creatures.begin(), creatures.end(), id) != creatures.end();
}

WorldProfile WorldProfile_Default(BodyKind kind, const std::string& name) {
    WorldProfile p; // the field defaults are today's planet (Terra)
    p.name = name;
    p.source = std::string("default:") + BodyKindName(kind);
    if (kind == BodyKind::Moon) { // today's moon (Luna): airless terrain, no New World choices, bake, foliage, air or sky
        p.airless = true;
        p.newWorld = p.erosion = p.foliage = p.planetPalette = false;
        p.hasAir = false; p.airTop = 0.0f;
        p.hasSky = false;
        p.allBiomes = false; p.biomes = { BiomeName(Biome::Polar) };
        p.hazards = { "vacuum" };
        // No creatures (#1329, the user 2026-11-21: the night hunter used to come out on moons too).
        p.allCreatures = false;
        p.creatures.clear();
    }
    return p;
}

namespace {
    bool ParseYesNo(const std::string& v, bool& out) {
        if (v == "yes" || v == "on" || v == "1") { out = true; return true; }
        if (v == "no" || v == "off" || v == "0") { out = false; return true; }
        return false;
    }
    bool ParsePositive(const std::string& s, float& out) {
        char* end = nullptr;
        out = strtof(s.c_str(), &end);
        return end && end != s.c_str() && *end == '\0' && out > 0.0f;
    }
    bool KnownBiome(const std::string& n) {
        for (int b = 0; b < (int)Biome::Count; ++b) if (n == BiomeName((Biome)b)) return true;
        return false;
    }
}

bool WorldProfile_ParseLine(const std::string& line, WorldProfile& p, std::string& error) {
    std::stringstream ss(line);
    std::string key;
    if (!(ss >> key)) { error = "empty line"; return false; }
    std::vector<std::string> v;
    for (std::string w; ss >> w;) v.push_back(w);
    auto one = [&](std::string& out) {
        if (v.size() != 1) { error = "'" + key + "' takes one value"; return false; }
        out = v[0];
        return true;
    };
    auto flag = [&](bool& out) {
        std::string s;
        if (!one(s)) return false;
        if (!ParseYesNo(s, out)) { error = "'" + key + "' is yes or no, not '" + s + "'"; return false; }
        return true;
    };
    if (key == "terrain") {
        std::string s;
        if (!one(s)) return false;
        if (s == "earthlike") p.airless = false;
        else if (s == "airless") p.airless = true;
        else { error = "terrain is earthlike or airless, not '" + s + "'"; return false; }
        return true;
    }
    if (key == "radius") {
        std::string s;
        if (!one(s)) return false;
        if (!ParsePositive(s, p.radius)) { error = "radius is a length in metres > 0, not '" + s + "'"; return false; }
        p.hasRadius = true;
        return true;
    }
    if (key == "gravity") {
        std::string s;
        if (!one(s)) return false;
        if (s == "auto") p.gravity = -1.0f;
        else {
            char* end = nullptr;
            p.gravity = strtof(s.c_str(), &end);
            if (!end || end == s.c_str() || *end != '\0' || p.gravity < 0.0f) { error = "gravity is auto or m/s^2 >= 0, not '" + s + "'"; return false; }
        }
        p.hasGravity = true;
        return true;
    }
    if (key == "newworld") return flag(p.newWorld);
    if (key == "erosion") return flag(p.erosion);
    if (key == "foliage") return flag(p.foliage);
    if (key == "palette") {
        std::string s;
        if (!one(s)) return false;
        if (s == "planet") p.planetPalette = true;
        else if (s == "scene") p.planetPalette = false;
        else { error = "palette is planet or scene, not '" + s + "'"; return false; }
        return true;
    }
    if (key == "air" || key == "sky") {
        const bool air = key == "air";
        if (v.size() == 1 && v[0] == "none") {
            if (air) { p.hasAir = false; p.airTop = 0.0f; } else p.hasSky = false;
            return true;
        }
        float f[3];
        if (v.size() != 3 || !ParsePositive(v[0], f[0]) || !ParsePositive(v[1], f[1]) || !ParsePositive(v[2], f[2])) {
            error = air ? "air is none or: top scale-height density (all > 0)" : "sky is none or: shell haze clouds (heights > 0)";
            return false;
        }
        if (air) { p.hasAir = true; p.airTop = f[0]; p.airScaleHeight = f[1]; p.airDensity = f[2]; }
        else { p.hasSky = true; p.skyHeight = f[0]; p.hazeHeight = f[1]; p.cloudHeight = f[2]; }
        return true;
    }
    if (key == "biomes") {
        if (v.empty()) { error = "biomes is all or a list of biome names"; return false; }
        if (v.size() == 1 && v[0] == "all") { p.allBiomes = true; p.biomes.clear(); return true; }
        for (const std::string& b : v)
            if (!KnownBiome(b)) { error = "unknown biome '" + b + "'"; return false; }
        p.allBiomes = false;
        p.biomes = v;
        return true;
    }
    if (key == "creatures") {
        if (v.empty()) { error = "creatures is all, none or a list of species ids"; return false; }
        p.creatures.clear();
        if (v.size() == 1 && v[0] == "all") { p.allCreatures = true; return true; }
        p.allCreatures = false;
        if (!(v.size() == 1 && v[0] == "none")) p.creatures = v;
        return true;
    }
    if (key == "resources" || key == "hazards") {
        std::vector<std::string>& list = key == "resources" ? p.resources : p.hazards;
        list = (v.size() == 1 && v[0] == "none") ? std::vector<std::string>{} : v;
        return true;
    }
    if (key == "about") {
        std::string rest;
        for (size_t i = 0; i < v.size(); ++i) rest += (i ? " " : "") + v[i];
        p.about = rest;
        return true;
    }
    error = "unknown key '" + key + "'";
    return false;
}

bool WorldProfile_LoadFile(const std::string& path, BodyKind kind, const std::string& name, WorldProfile& out) {
    out = WorldProfile_Default(kind, name);
    std::ifstream f(path);
    if (!f.is_open()) return true; // no file: the kind's default
    out.source = path;
    bool ok = true;
    std::string line;
    for (int lineNo = 1; std::getline(f, line); ++lineNo) {
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        std::string error;
        if (!WorldProfile_ParseLine(line, out, error)) {
            LOG_ERROR("Worlds", "%s:%d: %s\n", path.c_str(), lineNo, error.c_str());
            ok = false;
        }
    }
    return ok;
}

namespace {
    std::mutex gProfileMutex;
    std::unordered_map<std::string, std::shared_ptr<const WorldProfile>> gProfiles; // key: kind + name
}

std::shared_ptr<const WorldProfile> WorldProfile_For(const std::string& name, BodyKind kind) {
    if (kind != BodyKind::Planet && kind != BodyKind::Moon) return nullptr;
    const std::string key = std::string(BodyKindName(kind)) + ":" + name;
    std::lock_guard<std::mutex> lock(gProfileMutex);
    auto it = gProfiles.find(key);
    if (it != gProfiles.end()) return it->second;
    auto p = std::make_shared<WorldProfile>();
    WorldProfile_LoadFile("resources/data/worlds/" + name + ".txt", kind, name, *p);
    LOG_INFO("Worlds", "%s: %s (%s terrain%s%s%s)\n", name.c_str(), p->source.c_str(), p->airless ? "airless" : "earthlike",
             p->foliage ? ", foliage" : "", p->hasAir ? ", air" : "", p->hasSky ? ", sky" : "");
    gProfiles[key] = p;
    return p;
}

void WorldProfile_ClearCache() {
    std::lock_guard<std::mutex> lock(gProfileMutex);
    gProfiles.clear();
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1274): the world files parse, and Terra's and Luna's say what the game did before them ----
#include "unit_helpers.h"
#include <cstdio>

TEST(WorldProfile, DefaultsAreTodaysPlanetAndMoon) {
    const WorldProfile t = WorldProfile_Default(BodyKind::Planet, "x");
    EXPECT_FALSE(t.airless);
    EXPECT_TRUE(t.newWorld && t.erosion && t.foliage && t.planetPalette && t.hasAir && t.hasSky && t.allBiomes && t.allCreatures);
    EXPECT_EQ(t.skyHeight, 25000.0f); // Scene::Atmosphere's values
    EXPECT_EQ(t.hazeHeight, 8000.0f);
    EXPECT_EQ(t.cloudHeight, 1500.0f);
    EXPECT_EQ(t.airTop, 35000.0f);    // rocket.cpp kAirTop
    const WorldProfile m = WorldProfile_Default(BodyKind::Moon, "y");
    EXPECT_TRUE(m.airless);
    EXPECT_FALSE(m.newWorld || m.erosion || m.foliage || m.planetPalette || m.hasAir || m.hasSky || m.allBiomes);
    EXPECT_FALSE(m.allowsCreature("dusk_stalker")); // #1329: no night hunter on moons
    EXPECT_FALSE(m.allowsCreature("woolback"));
}

// The shipped files must not change what Terra and Luna are (#1274: "moved onto it unchanged").
TEST(WorldProfile, TerraAndLunaFilesMatchTheirDefaults) {
    for (BodyKind kind : { BodyKind::Planet, BodyKind::Moon }) {
        const std::string name = kind == BodyKind::Planet ? "terra" : "luna";
        WorldProfile f;
        ASSERT_TRUE(WorldProfile_LoadFile("resources/data/worlds/" + name + ".txt", kind, name, f)) << name;
        EXPECT_EQ(f.source, "resources/data/worlds/" + name + ".txt") << name << ": the file is missing";
        const WorldProfile d = WorldProfile_Default(kind, name);
        EXPECT_EQ(f.airless, d.airless) << name;
        EXPECT_EQ(f.newWorld, d.newWorld) << name;
        EXPECT_EQ(f.erosion, d.erosion) << name;
        EXPECT_EQ(f.foliage, d.foliage) << name;
        EXPECT_EQ(f.planetPalette, d.planetPalette) << name;
        EXPECT_EQ(f.hasSky, d.hasSky) << name;
        if (d.hasSky) {
            EXPECT_EQ(f.skyHeight, d.skyHeight) << name;
            EXPECT_EQ(f.hazeHeight, d.hazeHeight) << name;
            EXPECT_EQ(f.cloudHeight, d.cloudHeight) << name;
        }
        EXPECT_EQ(f.hasAir, d.hasAir) << name;
        if (d.hasAir) {
            EXPECT_EQ(f.airTop, d.airTop) << name;
            EXPECT_EQ(f.airScaleHeight, d.airScaleHeight) << name;
            EXPECT_EQ(f.airDensity, d.airDensity) << name;
        }
        EXPECT_EQ(f.allBiomes, d.allBiomes) << name;
        EXPECT_EQ(f.biomes, d.biomes) << name;
        for (const char* sp : { "woolback", "plains_runner", "clucker", "hornhide", "thornback", "dusk_stalker" })
            EXPECT_EQ(f.allowsCreature(sp), d.allowsCreature(sp)) << name << " " << sp;
    }
}

TEST(WorldProfile, ParseLineAcceptsAndRejects) {
    WorldProfile p = WorldProfile_Default(BodyKind::Planet, "z");
    std::string e;
    EXPECT_TRUE(WorldProfile_ParseLine("terrain airless", p, e)); EXPECT_TRUE(p.airless);
    EXPECT_TRUE(WorldProfile_ParseLine("sky none", p, e)); EXPECT_FALSE(p.hasSky);
    EXPECT_TRUE(WorldProfile_ParseLine("sky 9000 2000 800", p, e)); EXPECT_TRUE(p.hasSky); EXPECT_EQ(p.hazeHeight, 2000.0f);
    EXPECT_TRUE(WorldProfile_ParseLine("air 12000 3000 0.4", p, e)); EXPECT_EQ(p.airDensity, 0.4f);
    EXPECT_TRUE(WorldProfile_ParseLine("biomes Desert Badlands", p, e)); EXPECT_FALSE(p.allBiomes); EXPECT_EQ(p.biomes.size(), 2u);
    EXPECT_TRUE(WorldProfile_ParseLine("creatures none", p, e)); EXPECT_FALSE(p.allowsCreature("woolback"));
    EXPECT_TRUE(WorldProfile_ParseLine("creatures woolback", p, e)); EXPECT_TRUE(p.allowsCreature("woolback")); EXPECT_FALSE(p.allowsCreature("clucker"));
    EXPECT_TRUE(WorldProfile_ParseLine("about a small grey moon", p, e)); EXPECT_EQ(p.about, "a small grey moon");
    EXPECT_FALSE(WorldProfile_ParseLine("terrain soupy", p, e));
    EXPECT_FALSE(WorldProfile_ParseLine("biomes Moonscape", p, e));
    EXPECT_FALSE(WorldProfile_ParseLine("sky 1 2", p, e));
    EXPECT_FALSE(WorldProfile_ParseLine("foliage maybe", p, e));
    EXPECT_FALSE(WorldProfile_ParseLine("colour red", p, e));
}

// (#1298) radius and gravity moved from bodies.txt into the world files: the loaded system, and so the lines saves
// and joining players get, must be exactly what bodies.txt said before the move.
TEST(WorldProfile, BodiesFileResolvesToTheOldDefinitions) {
    WorldProfile_ClearCache();
    SystemDef sys;
    ASSERT_TRUE(LoadSystemDef("resources/data/bodies.txt", sys));
    ASSERT_GE(sys.bodies.size(), 2u);
    EXPECT_EQ(FormatBodyDef(sys.bodies[0]), "planet terra radius=60000 position=0,0,0 day=3600");
    EXPECT_EQ(FormatBodyDef(sys.bodies[1]), "moon luna radius=10000 orbit=terra,500000,30 day=86400");
    EXPECT_EQ(sys.bodies[0].gravity, -1.0f); // auto, as before
    EXPECT_EQ(sys.bodies[1].gravity, -1.0f);
}

TEST(WorldProfile, BodiesLineOverridesTheWorldFile) {
    const char* path = "worldprofile_test_bodies.tmp";
    {
        std::ofstream f(path);
        f << "planet terra radius=500 gravity=2 position=0,0,0 spin=0\n";
        f << "moon luna orbit=terra,1000,0 spin=0\n";
        f << "moon nowhere_moon orbit=terra,2000,0 spin=0\n"; // no world file and no radius: an error
    }
    WorldProfile_ClearCache();
    SystemDef sys;
    EXPECT_FALSE(LoadSystemDef(path, sys)); // the moon without a radius
    std::ofstream(path, std::ios::trunc) << "planet terra radius=500 gravity=2 position=0,0,0 spin=0\nmoon luna orbit=terra,1000,0 spin=0\n";
    ASSERT_TRUE(LoadSystemDef(path, sys));
    std::remove(path);
    EXPECT_EQ(sys.bodies[0].radius, 500.0f);
    EXPECT_EQ(sys.bodies[0].gravity, 2.0f);
    EXPECT_EQ(sys.bodies[1].radius, 10000.0f); // from worlds/luna.txt
}

TEST(WorldProfile, AsteroidsHaveNone) {
    EXPECT_EQ(WorldProfile_For("goliath", BodyKind::Asteroid), nullptr);
    EXPECT_EQ(WorldProfile_For("rock", BodyKind::Ring), nullptr);
}
#endif
