#ifndef OBJECT_H
#define OBJECT_H

#include <string>
#include <vector>
#include <memory>
#include <raylib.h>
#include <iostream>

#include "frustum.h"
#include "dvec3.h"

// Draws many meshes with one material: binds its shader and texture maps once, then per mesh only
// sets mvp / matModel and draws (raylib's DrawMesh re-binds and unbinds everything on every call,
// which dominated the frame with hundreds of terrain chunks). Mirrors DrawMesh's OpenGL 3.3 path;
// meshes must have been uploaded with a VAO (UploadMesh). Unbinds in the destructor.
struct TextureArray;

class ChunkBatch {
public:
    // `layers` (optional): a texture array bound to unit kLayersUnit and given to sampler `layersLoc`.
    explicit ChunkBatch(const Material& material, const TextureArray* layers = nullptr, int layersLoc = -1);
    static constexpr int kLayersUnit = 11; // a free 2D slot among raylib's 12 material maps
    ~ChunkBatch();
    ChunkBatch(const ChunkBatch&) = delete;
    ChunkBatch& operator=(const ChunkBatch&) = delete;
    void setUniform(int loc, const void* value, int uniformType); // the shader is already bound
    void draw(const Mesh& mesh, const Matrix& model);

private:
    const Material& material;
    const TextureArray* layers;
    Matrix view, projection, viewProjection;
};

// Everything an object needs to submit draws for one pass. Owned by the Scene; objects never
// look up shaders or uniforms themselves, they hand meshes + a material to DrawMesh, which sets
// mvp / matModel / matNormal and binds the material's texture maps.
struct RenderContext {
    Material* terrainMaterial = nullptr;  // lit terrain material (lighting shader + textures + shadow map)
    Material* depthMaterial = nullptr;    // depth-only material for the shadow pass
    const Frustum* frustum = nullptr;     // frustum of the pass's camera (main or light), for culling
    int* drawCalls = nullptr;             // optional stats counter
    int chunkOriginLoc = -1;              // lighting shader "chunkTexOrigin" (chunk origin wrapped to the texture repeat)
    DVec3 renderOrigin{};                 // the scene's floating origin (chunk placement in double)
    float texScale = 0.05f;               // terrain texture repeats per metre (wraps chunkTexOrigin)
    ChunkBatch* batch = nullptr;          // set while a body draws its chunks (see ChunkBatch)
    const TextureArray* terrainLayers = nullptr; // terrain texture array (planet palette)
    int terrainLayersLoc = -1;            // its sampler uniform
    int paletteLoc = -1;                  // lighting shader "palette"; bodies may override the scene's
    int scenePalette = 0;
    int sunOccluderLoc = -1;              // lighting shader "sunOccluder" (a planet blocking the sun)
    Material* waterMaterial = nullptr;    // sea surface (water shader); drawn in the transparent pass
    int waterTexOriginLoc = -1, waterClipLoc = -1, waterOccluderLoc = -1;
    int chunkClipLoc = -1;                // lighting shader "chunkClipMax" (chunk-local overlap clipping)
    Shader* grassShader = nullptr;        // foliage: grass tufts (drawn after the opaque terrain)
    int grassModelLoc = -1, grassOccluderLoc = -1;
    Vector3 cameraPos = { 0, 0, 0 };      // render space
    Shader* farShader = nullptr;          // distance sphere (farsphere.h)
    Shader* partShader = nullptr;         // build plots (plots.h)
    int partOccluderLoc = -1;
    Shader* treeShader = nullptr;         // foliage: trees and bushes
    int treeModelLoc = -1, treeOccluderLoc = -1, treeFadeLoc = -1;
};

extern bool frustumCullingEnabled; // debug toggle

class Object {
public:
    std::string type; // Type of the object (e.g., "cube", "sphere")
    std::string name; // Name of the object
    Vector3 position; // Position of the object
    Vector3 rotation; // Rotation of the object (degrees)
    Color color; // Color of the object
    float scale; // Scale of the object

    Object* parent = nullptr;
    std::vector<std::unique_ptr<Object>> children;

    bool isActive = true; // Whether the object is active in the scene

    // Top-level objects with hasWorldPos keep their authoritative position in doubles; `position`
    // is then derived each frame as worldPos - Scene::renderOrigin (floating origin).
    DVec3 worldPos{};
    bool hasWorldPos = false;

    Object(std::string type, std::string name, Vector3 position, Vector3 rotation, Color color, float scale);
    virtual void draw(const RenderContext& ctx);
    virtual void drawDepthOnly(const RenderContext& ctx);
    virtual void drawTransparent(const RenderContext&) {} // after every opaque object (water)
    virtual ~Object(); // Make destructor virtual for safe polymorphic deletion

    // Scale, then rotate (XYZ, degrees), then translate.
    Matrix localTransform() const;

    // Disable copying and moving: derived classes own GPU resources.
    Object(const Object&) = delete;
    Object& operator=(const Object&) = delete;
};

#endif
