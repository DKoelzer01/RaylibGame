#ifndef MAPFEATURES_H
#define MAPFEATURES_H

// The live world on the minimap and the big map (#1795): build plots drawn top-down, Factorio style.
//
// Each plot gets a footprint: per column of cells (plot x, z) the part shown from above, coloured by its kind. A machine
// shows over the floor or roof around it. Footprints are built on the job system from a copy of the plot's parts, only when
// the plot's version changes. The map's render workers then sample them per pixel. Zoomed in, each part has an outline.
// Zoomed out, coarser levels (2, 4, 8... cells a block, the highest-ranked kind in it) keep machines and thin rows of parts
// visible as a solid footprint. Colours come from resources/data/map_colours.txt.

#include <raylib.h>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include "dvec3.h"

class Planetoid;
enum class PartType : uint8_t;

enum class MapKind : uint8_t { None = 0, Structure, Machine, Transport, Storage, Power, Farm, Other, Count };
MapKind MapKindOf(PartType t);
const char* MapKindName(MapKind k);
Color MapFeatures_Colour(MapKind k); // (main thread, after the first MapFeatures_Update)

struct MapIcon {
    std::string item;  // the machine's item (its icon)
    Vector3 dir{};     // body-local unit direction of its centre
    float size = 1.0f; // m across
};

struct PlotMap {
    uint32_t id = 0, version = 0;
    DVec3 origin;              // body-local corner of cell (0,0,0)
    Vector3 east{}, north{};   // body-local plot x and z axes
    DVec3 centre;              // body-local centre of the footprint
    float radius = 0.0f;       // m from the centre to its farthest corner
    int x0 = 0, z0 = 0;        // the first column (plot cells)
    struct Level { int w = 0, h = 0; std::vector<uint8_t> kind; }; // level L: blocks of 2^L x 2^L columns
    std::vector<Level> levels;
    std::vector<int32_t> part; // level 0: the part shown per column (-1 none), for outlines
    std::vector<MapIcon> icons;
};

// The plot footprints one map image needs, immutable once made (the render workers read it).
struct MapPlots {
    std::vector<std::shared_ptr<const PlotMap>> plots;
    Color colour[(int)MapKind::Count]{};
};

// Main thread, each frame: starts a footprint build for each plot of `body` whose version changed. True when a build landed
// or a plot went away since the last call (the maps are due a re-render).
bool MapFeatures_Update(const Planetoid& body);
// Main thread: the footprints of `body`'s plots within `radius` m of the body-local point `centre`; nullptr if none.
std::shared_ptr<const MapPlots> MapFeatures_PlotsNear(const Planetoid& body, const DVec3& centre, float radius);
// Any thread: the colour of the plots at body-local `p` in an image of `mpp` m a pixel with body-local axes east, north.
// False where no part is.
bool MapFeatures_Sample(const MapPlots& s, const DVec3& p, Vector3 east, Vector3 north, float mpp, Color& out);
// Tests: a footprint from a plot's parts, built at once on this thread.
std::shared_ptr<const PlotMap> MapFeatures_BuildNow(const class Plot& plot);
void MapFeatures_Reset();

#endif
