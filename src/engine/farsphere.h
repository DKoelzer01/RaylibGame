#ifndef FARSPHERE_H
#define FARSPHERE_H

// Distance sphere: planets and moons seen from far away (beyond kFarSphereRadii body radii, the
// telescope's zoom taken into account) are drawn as one textured sphere instead of their terrain
// meshes. The textures are baked once per body on the worker threads from the same height map and
// biome rules as the terrain (sea shaded by depth, forest canopy painted in): an equirectangular
// colour map and an object-space normal map (from the height grid, so mountains still shade).
// The LOD keeps its coarse nodes meanwhile, so turning back toward the body needs no streaming.

#include <raylib.h>
#include <atomic>
#include <memory>
#include <vector>

class PlanetField;
class Planetoid;
struct RenderContext;

constexpr float kFarSphereRadii = 3.0f;   // the coarse meshes look rougher than the bake past here

class FarSphere {
public:
    FarSphere(std::shared_ptr<const PlanetField> field, bool foliage, unsigned int foliageSeed);
    ~FarSphere();
    bool ready();                         // main thread: uploads the bake once its rows are done
    bool isReady() const { return uploaded; }
    struct Bake;                          // pixels, shared with the worker jobs
    void draw(const Planetoid& body, const RenderContext& ctx, Shader shader);

private:
    std::shared_ptr<Bake> bake;
    Texture2D albedo{}, normals{};
    Mesh mesh{};
    bool uploaded = false;
};

#endif
