#include "worldgenparams.h"
#include "bodies.h"
#include "world.h"

#include <cstdio>
#include <cstring>
#include <random>

bool WorldGenParams::isDefault() const {
    const WorldGenParams d;
    return ocean == d.ocean && mountains == d.mountains && roughness == d.roughness && biomeSize == d.biomeSize &&
           resources == d.resources && erosion == d.erosion && erosionStrength == d.erosionStrength;
}

uint64_t WorldGenParams::key() const {
    uint64_t h = 1469598103934665603ull; // FNV-1a
    auto mix = [&](const void* p, size_t n) { const unsigned char* b = (const unsigned char*)p; for (size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ull; } };
    mix(&ocean, sizeof ocean); mix(&mountains, sizeof mountains); mix(&roughness, sizeof roughness);
    mix(&biomeSize, sizeof biomeSize); mix(&resources, sizeof resources);
    const int e = erosion ? 1 : 0; mix(&e, sizeof e); mix(&erosionStrength, sizeof erosionStrength);
    return h;
}

WorldGenParams WorldGenParams_ForSave(const std::string& savePath) {
    WorldGenParams p;
    FILE* f = fopen((savePath + ".worldgen").c_str(), "r");
    if (!f) return p;
    char key[64]; float v;
    while (fscanf(f, " %63s %f", key, &v) == 2) {
        if (!strcmp(key, "ocean")) p.ocean = v;
        else if (!strcmp(key, "mountains")) p.mountains = v;
        else if (!strcmp(key, "roughness")) p.roughness = v;
        else if (!strcmp(key, "biomeSize")) p.biomeSize = v;
        else if (!strcmp(key, "resources")) p.resources = v;
        else if (!strcmp(key, "erosion")) p.erosion = v != 0.0f;
        else if (!strcmp(key, "erosionStrength")) p.erosionStrength = v;
    }
    fclose(f);
    return p;
}

bool WorldGenParams_Save(const std::string& savePath, const WorldGenParams& p) {
    FILE* f = fopen((savePath + ".worldgen").c_str(), "w");
    if (!f) return false;
    // %.9g: the floats come back bit for bit (the bake's cache key hashes them).
    fprintf(f, "ocean %.9g\nmountains %.9g\nroughness %.9g\nbiomeSize %.9g\nresources %.9g\nerosion %d\nerosionStrength %.9g\n",
            p.ocean, p.mountains, p.roughness, p.biomeSize, p.resources, p.erosion ? 1 : 0, p.erosionStrength);
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
