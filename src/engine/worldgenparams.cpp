#include "worldgenparams.h"
#include "bodies.h"
#include "world.h"
#include <sstream>

#include <cstdio>
#include <cstring>
#include <random>

bool WorldGenParams::isDefault() const {
    const WorldGenParams d;
    return ocean == d.ocean && mountains == d.mountains && roughness == d.roughness && biomeSize == d.biomeSize &&
           resources == d.resources && erosion == d.erosion && erosionStrength == d.erosionStrength && cliffs == d.cliffs &&
           karst == d.karst && volcanoes == d.volcanoes && waterCycle == d.waterCycle;
}

uint64_t WorldGenParams::key() const {
    uint64_t h = 1469598103934665603ull; // FNV-1a
    auto mix = [&](const void* p, size_t n) { const unsigned char* b = (const unsigned char*)p; for (size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ull; } };
    mix(&ocean, sizeof ocean); mix(&mountains, sizeof mountains); mix(&roughness, sizeof roughness);
    mix(&biomeSize, sizeof biomeSize); mix(&resources, sizeof resources);
    const int e = erosion ? 1 : 0; mix(&e, sizeof e); mix(&erosionStrength, sizeof erosionStrength);
    mix(&cliffs, sizeof cliffs); // #589
    if (!karst) { const int k = 0x592; mix(&k, sizeof k); }       // #592 (only when off: default keys unchanged)
    if (!volcanoes) { const int v = 0x593; mix(&v, sizeof v); }   // #593
    if (waterCycle != 1.0f) mix(&waterCycle, sizeof waterCycle); // #637 (only when changed: default keys unchanged)
    return h;
}

WorldGenParams WorldGenParams_FromText(const std::string& text) {
    WorldGenParams p;
    std::istringstream in(text);
    std::string key; float v;
    while (in >> key >> v) {
        if (key == "ocean") p.ocean = v;
        else if (key == "mountains") p.mountains = v;
        else if (key == "roughness") p.roughness = v;
        else if (key == "biomeSize") p.biomeSize = v;
        else if (key == "resources") p.resources = v;
        else if (key == "erosion") p.erosion = v != 0.0f;
        else if (key == "erosionStrength") p.erosionStrength = v;
        else if (key == "waterCycle") p.waterCycle = v;
        else if (key == "cliffs") p.cliffs = v;
        else if (key == "karst") p.karst = v != 0.0f;
        else if (key == "volcanoes") p.volcanoes = v != 0.0f;
    }
    return p;
}

std::string WorldGenParams_ToText(const WorldGenParams& p) {
    // %.9g: the floats come back bit for bit (the bake's cache key hashes them).
    char buf[512];
    snprintf(buf, sizeof buf, "ocean %.9g\nmountains %.9g\nroughness %.9g\nbiomeSize %.9g\nresources %.9g\nerosion %d\nerosionStrength %.9g\ncliffs %.9g\nkarst %d\nvolcanoes %d\nwaterCycle %.9g\n",
             p.ocean, p.mountains, p.roughness, p.biomeSize, p.resources, p.erosion ? 1 : 0, p.erosionStrength, p.cliffs, p.karst ? 1 : 0, p.volcanoes ? 1 : 0, p.waterCycle);
    return buf;
}

WorldGenParams WorldGenParams_ForSave(const std::string& savePath) {
    FILE* f = fopen((savePath + ".worldgen").c_str(), "r");
    if (!f) return WorldGenParams{};
    std::string text; char buf[512]; size_t n;
    while ((n = fread(buf, 1, sizeof buf, f)) > 0) text.append(buf, n);
    fclose(f);
    return WorldGenParams_FromText(text);
}

bool WorldGenParams_Save(const std::string& savePath, const WorldGenParams& p) {
    FILE* f = fopen((savePath + ".worldgen").c_str(), "w");
    if (!f) return false;
    fputs(WorldGenParams_ToText(p).c_str(), f);
    fclose(f);
    return true;
}

bool WorldGen_StartPlanet(unsigned int worldSeed, unsigned int& terrainSeed, float& radius) {
    if (!worldGen.systemLoaded) ReloadSystemDef();
    // The same draws as generateWorld (world.cpp): rings round a parent draw nothing; a loose belt draws per member;
    // every other body draws 3 positions, its terrain seed, then its spin.
    std::mt19937 gen(worldSeed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    bool seenParent = false;
    for (const BodyDef& d : worldGen.system.bodies) {
        if (d.kind == BodyKind::Ring && !d.orbitParent.empty()) continue;
        if (d.kind == BodyKind::Ring) {
            for (int k = 0; k < d.count; ++k) { (void)unit(gen); (void)unit(gen); (void)unit(gen); (void)gen(); (void)unit(gen); (void)unit(gen); }
            continue;
        }
        for (int i = 0; i < 3; ++i) (void)unit(gen);
        const unsigned int drawn = (unsigned int)gen();
        (void)unit(gen);
        if (d.kind == BodyKind::Planet) { terrainSeed = d.hasSeed ? d.seed : drawn; radius = d.radius; return true; }
        seenParent = true;
    }
    (void)seenParent;
    return false;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#869): the world generation sidecar and the Welcome's worldgen text ----
#include "unit_helpers.h"
#include <iterator>

namespace {
// Every field set away from its default, one at a time (the name is for the failure message).
struct FieldCase { const char* name; void (*set)(WorldGenParams&); };
const FieldCase kFields[] = {
    { "ocean", [](WorldGenParams& p) { p.ocean = -0.375f; } },
    { "mountains", [](WorldGenParams& p) { p.mountains = 1.7f; } },
    { "roughness", [](WorldGenParams& p) { p.roughness = 0.31f; } },
    { "biomeSize", [](WorldGenParams& p) { p.biomeSize = 1.9f; } },
    { "resources", [](WorldGenParams& p) { p.resources = 0.6f; } },
    { "erosion", [](WorldGenParams& p) { p.erosion = false; } },
    { "erosionStrength", [](WorldGenParams& p) { p.erosionStrength = 1.3f; } },
    { "waterCycle", [](WorldGenParams& p) { p.waterCycle = 1.75f; } }, // #686: once dropped from the text
    { "cliffs", [](WorldGenParams& p) { p.cliffs = 0.1f; } },
    { "karst", [](WorldGenParams& p) { p.karst = false; } },
    { "volcanoes", [](WorldGenParams& p) { p.volcanoes = false; } },
};
}

// A new field must be added to the text, the key, isDefault and kFields above: this fails until it is.
TEST(WorldGenParams, EveryFieldIsListedHere) {
    EXPECT_EQ(sizeof(WorldGenParams), 40u) << "WorldGenParams changed size: add the new field to ToText, FromText, key(), "
                                              "isDefault() and kFields in this test";
    EXPECT_EQ(std::size(kFields), 11u);
}

TEST(WorldGenParams, DefaultsRoundTrip) {
    const WorldGenParams d;
    EXPECT_TRUE(d.isDefault());
    const WorldGenParams back = WorldGenParams_FromText(WorldGenParams_ToText(d));
    EXPECT_TRUE(back.isDefault());
    EXPECT_EQ(back.key(), d.key());
    EXPECT_TRUE(WorldGenParams_FromText("").isDefault()) << "an old save without a sidecar generates as before";
}

TEST(WorldGenParams, EveryFieldRoundTripsThroughTheText) {
    const uint64_t defaultKey = WorldGenParams{}.key();
    for (const FieldCase& f : kFields) {
        WorldGenParams p;
        f.set(p);
        EXPECT_FALSE(p.isDefault()) << f.name << ": isDefault() ignores it";
        EXPECT_NE(p.key(), defaultKey) << f.name << ": the bake's cache key ignores it";
        const std::string text = WorldGenParams_ToText(p);
        const WorldGenParams back = WorldGenParams_FromText(text);
        EXPECT_EQ(back.key(), p.key()) << f.name << " lost in the text:\n" << text;
        EXPECT_FALSE(back.isDefault()) << f.name << " lost in the text:\n" << text;
    }
}

TEST(WorldGenParams, FloatsComeBackBitForBit) {
    WorldGenParams p;
    p.ocean = 0.1f; p.mountains = 1.0f / 3.0f; p.waterCycle = 0.7f + 1e-7f; p.cliffs = 1.9999999f;
    const WorldGenParams back = WorldGenParams_FromText(WorldGenParams_ToText(p));
    EXPECT_EQ(memcmp(&back.ocean, &p.ocean, sizeof(float)), 0);
    EXPECT_EQ(memcmp(&back.mountains, &p.mountains, sizeof(float)), 0);
    EXPECT_EQ(memcmp(&back.waterCycle, &p.waterCycle, sizeof(float)), 0);
    EXPECT_EQ(back.key(), p.key());
}

TEST(WorldGenParams, UnknownKeysAreSkipped) {
    const WorldGenParams p = WorldGenParams_FromText("futureKnob 3\nmountains 1.5\nanother 0\n");
    EXPECT_FLOAT_EQ(p.mountains, 1.5f);
    EXPECT_FLOAT_EQ(p.roughness, 1.0f);
}

TEST(WorldGenParams, SidecarRoundTrip) {
    const std::string save = UnitTest::ScratchDir() + "/world.sav";
    WorldGenParams p; p.ocean = 0.25f; p.waterCycle = 0.5f; p.karst = false;
    ASSERT_TRUE(WorldGenParams_Save(save, p));
    EXPECT_EQ(WorldGenParams_ForSave(save).key(), p.key());
    EXPECT_TRUE(WorldGenParams_ForSave(save + ".missing").isDefault());
}
#endif
