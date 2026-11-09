#ifndef OBJECT_H
#define OBJECT_H

#include <string>
#include <vector>
#include <memory>
#include <raylib.h>
#include <iostream>

#include "frustum.h"
#include "dvec3.h"

// Everything an object needs to submit draws for one pass. Owned by the Scene; objects never
// look up shaders or uniforms themselves, they hand meshes + a material to DrawMesh, which sets
// mvp / matModel / matNormal and binds the material's texture maps.
struct RenderContext {
    Material* terrainMaterial = nullptr;  // lit terrain material (lighting shader + textures + shadow map)
    Material* depthMaterial = nullptr;    // depth-only material for the shadow pass
    const Frustum* frustum = nullptr;     // frustum of the pass's camera (main or light), for culling
    int* drawCalls = nullptr;             // optional stats counter
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
    virtual ~Object(); // Make destructor virtual for safe polymorphic deletion

    // Scale, then rotate (XYZ, degrees), then translate.
    Matrix localTransform() const;

    // Disable copying and moving: derived classes own GPU resources.
    Object(const Object&) = delete;
    Object& operator=(const Object&) = delete;
};

#endif
