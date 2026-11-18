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
