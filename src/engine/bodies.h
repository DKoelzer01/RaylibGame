#ifndef BODIES_H
#define BODIES_H

// Star-system definitions: what bodies exist, read from resources/data/bodies.txt instead of
// being generated in code. One body (or one procedural ring of bodies) per line:
//
//   kind name key=value ...
//
// kinds    asteroid | planet | moon   one body
//          ring                       `count` procedural asteroids scattered around the origin
// keys     radius=R          body radius in metres (ring: base size of its members)
//          position=x,y,z    world position, or
//          orbit=parent,altitude,phaseDeg   placed `altitude` above the parent's surface, at
//                            `phaseDeg` around its equator (motion arrives in milestone 1 step 6)
//          spin=deg/s  |  day=seconds        rotation about the body's Y axis (day<0: retrograde)
//          gravity=auto | m/s^2              auto: asteroids 3 + 0.03 R (max 9.81),
//                                            planets/moons 9.81 * R / 100 km
//          seed=N            terrain seed (default: drawn from the system seed)
//          ring only: count=N  size=lo,hi (radius multipliers)  extent=m  flatten=f  maxspin=deg/s
//
// Bodies draw from the system seed in a fixed order (3 position values, the terrain seed, a
// size value for ring members, a spin value), whether or not a key overrides them, so a seed
// reproduces the same system as the original hard-coded generator.

#include <raylib.h>
#include <string>
#include <vector>

enum class BodyKind { Asteroid, Planet, Moon, Ring };

struct BodyDef {
    BodyKind kind = BodyKind::Asteroid;
    std::string name;
    float radius = 50.0f;
    bool hasPosition = false;
    Vector3 position = { 0, 0, 0 };
    std::string orbitParent;          // empty: use position
    float orbitAltitude = 0.0f, orbitPhaseDeg = 0.0f;
    float spinDegPerSec = 0.0f;       // set from spin= or day=
    bool spinFromDay = false;         // remember which key was used (for writing it back)
    float dayLength = 0.0f;
    float gravity = -1.0f;            // < 0: auto
    bool hasSeed = false;
    unsigned int seed = 0;
    // ring
    int count = 0;
    float sizeMin = 0.5f, sizeMax = 3.0f;
    float extent = 1500.0f, maxSpinDegPerSec = 3.0f;
    double flatten = 0.2;             // double: the y scatter reproduces the old generator exactly
};

struct SystemDef {
    std::vector<BodyDef> bodies;
};

const char* BodyKindName(BodyKind kind);
// Parses one definition line (no comment or blank handling). On failure returns false and sets error.
bool ParseBodyDef(const std::string& line, BodyDef& out, std::string& error);
std::string FormatBodyDef(const BodyDef& def);          // round-trips through ParseBodyDef
bool LoadSystemDef(const std::string& path, SystemDef& out); // logs errors with line numbers
// The system the hard-coded generator used to make, from its old settings (old saves, fallback).
SystemDef LegacySystemDef(int count, int size, float systemRadius, float maxSpin, Vector3 first);
float AutoSurfaceGravity(BodyKind kind, float radius);

#endif
