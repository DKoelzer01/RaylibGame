#ifndef LOCATE_H
#define LOCATE_H

// /locate and /tp (#1171): find the world's features from where the player stands, for players (creative) and testers who
// can't find a cave, a lava tube or a ruin by walking. Each feature uses the finder its test uses (caves: the cave graph's
// mouths; volcanoes: the designed relief's summits; ruins: the structure sites; springs and seeps: the source grid...).
// Everything is computed from the seed and the bake, so a joined client's world would give the same answers.

#include <raylib.h>
#include <string>
#include <vector>

class Planetoid;

struct LocateHit {
    Vector3 at{};        // body-local: on the ground over the feature (a cave mouth or chamber: the feature itself)
    std::string what;    // "lava tube skylight", "ruin: Ruined watchtower (dungeon)", "biome: taiga"...
    float distance = 0.0f; // metres from `from`, straight
    float bearing = 0.0f;  // degrees clockwise from north (the body's +Y pole), at `from`
};

// Up to n of `feature` nearest `from` (body-local), nearest first; false with why (unknown feature, none in range).
bool Locate_Find(const Planetoid& body, Vector3 from, const std::string& feature, int n, std::vector<LocateHit>& out, std::string& why);
std::string Locate_FeatureList();               // for /help and an unknown feature's answer
const char* Locate_CompassPoint(float bearing); // "N", "NE", ...

// Latitude and longitude (degrees) <-> a body-local direction, the tests' convention (RAYLIBGAME_MOUTH_AT):
// (cos lat cos lon, sin lat, cos lat sin lon).
Vector3 Locate_DirFromLatLon(float latDeg, float lonDeg);
void Locate_LatLonOf(Vector3 local, float& latDeg, float& lonDeg);

#endif
