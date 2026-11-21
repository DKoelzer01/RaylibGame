#ifndef WORLDPROFILE_H
#define WORLDPROFILE_H

// What a planet or moon is (Plan 13 step 2, #1274): one world file per body, resources/data/worlds/<body name>.txt,
// read where the game used to decide by the body's kind. bodies.txt keeps what a body does in the system (radius,
// orbit, spin, gravity: saved with the world and sent to joining players); the world file holds what it is like to be
// there (terrain style, foliage, air, sky, biomes, resources, creatures). Bodies without a file get their kind's
// default, which is today's planet (Terra) or moon (Luna), so a missing file changes nothing.
//
// One `key value...` per line, '#' starts a comment, unknown keys are errors (logged with the line):
//   terrain   earthlike | airless   PlanetField's style: earthlike = sea, rivers, erosion, caves, the 3D band, biomes;
//                                   airless = cratered regolith, no water, one cold biome
//   newworld  yes | no              the New World screen's choices (#504) shape this body (the home world)
//   erosion   yes | no              the Plan 7 erosion bake runs for it (earthlike terrain only)
//   foliage   yes | no              trees, bushes and grass (and the grass blockers under parts)
//   palette   planet | scene        the terrain's texture palette: the layer palette, or the scene's
//   air       none | top scale density   the air: its top (m above sea level), scale height (m), sea-level density (kg/m^3)
//   sky       none | shell haze clouds   an atmosphere drawn over the stars: its height, haze layer and cloud layer (m)
//   biomes    all | Name...         the biomes it has (voxel.h Biome names)
//   resources id...                 item ids it is a source of (declared; the destination steps hook them up)
//   hazards   word...               what it takes to go there (vacuum, heat, cold, radiation...; declared)
//   creatures all | none | id...    animal species (animals.txt ids) that spawn on it
//   about     text...               one line for the telescope and the help pages

#include "bodies.h"

#include <memory>
#include <string>
#include <vector>

struct WorldProfile {
    std::string name;                 // the body (and file) name
    std::string source;               // the file it came from, or "default:<kind>"
    bool airless = false;             // PlanetField's `moon` style
    bool newWorld = true;
    bool erosion = true;
    bool foliage = true;
    bool planetPalette = true;
    // Air (#1274: declared here; rocket.cpp drag and frames.cpp's "in orbit" altitude keep their own 35 km until
    // Plan 11 F5 moves them onto this, see #667).
    bool hasAir = true;
    float airTop = 35000.0f, airScaleHeight = 4000.0f, airDensity = 1.2f;
    // Sky (Scene::Atmosphere): drawn when hasSky.
    bool hasSky = true;
    float skyHeight = 25000.0f, hazeHeight = 8000.0f, cloudHeight = 1500.0f;
    bool allBiomes = true;
    std::vector<std::string> biomes;  // when !allBiomes
    std::vector<std::string> resources, hazards;
    bool allCreatures = true;
    std::vector<std::string> creatures; // when !allCreatures (empty: none)
    std::string about;

    bool allowsCreature(const std::string& id) const;
};

// Today's values for a kind (planet: Terra, moon: Luna); asteroids and rings have no world profile.
WorldProfile WorldProfile_Default(BodyKind kind, const std::string& name);
// Applies one line (no comment handling) to `p`. On failure returns false and sets error.
bool WorldProfile_ParseLine(const std::string& line, WorldProfile& p, std::string& error);
// The kind's default with the file's lines applied; errors are logged with line numbers (the bad lines are skipped).
// Returns false when the file exists but has errors; `out` is still usable.
bool WorldProfile_LoadFile(const std::string& path, BodyKind kind, const std::string& name, WorldProfile& out);
// The profile for a body: resources/data/worlds/<name>.txt over its kind's default, read once and cached.
// nullptr for asteroids and rings.
std::shared_ptr<const WorldProfile> WorldProfile_For(const std::string& name, BodyKind kind);
void WorldProfile_ClearCache();       // the next WorldProfile_For reads the files again (debug UI "Reload bodies file")

#endif
