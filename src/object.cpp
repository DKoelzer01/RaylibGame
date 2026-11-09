#include "object.h"
#include "textobject.h"
#include "gameobject.h"
#include <raymath.h>
#include <cstring>

Object::Object(std::string type,std::string name, Vector3 position, Vector3 rotation, Color color, float scale)
    : type(type), name(name), position(position), rotation(rotation), color(color), scale(scale) {}

Matrix Object::localTransform() const {
    Matrix m = MatrixScale(scale, scale, scale);
    m = MatrixMultiply(m, MatrixRotateXYZ({ DEG2RAD * rotation.x, DEG2RAD * rotation.y, DEG2RAD * rotation.z }));
    return MatrixMultiply(m, MatrixTranslate(position.x, position.y, position.z));
}

void Object::draw(const RenderContext&) {
    // Base objects do not render anything; derived classes implement their own draw logic.
}

void Object::drawDepthOnly(const RenderContext&) {
    // Base objects do not render depth.
}

Object::~Object() {}

bool frustumCullingEnabled = true;

TextObject::TextObject(std::string text, int fontSize)
    : Object("text","text", {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, WHITE, 1.0f), text(text), fontSize(fontSize) {}

TextObject::TextObject(std::string text, Vector3 position, Vector3 rotation, Color color, float scale, int fontSize)
    : Object("text","text", position, rotation, color, scale), text(text), fontSize(fontSize) {}

void TextObject::draw(const RenderContext&) {
    DrawText(text.c_str(), static_cast<int>(position.x), static_cast<int>(position.y), fontSize, color);
}

GameObject::GameObject(std::string type,std::string name, Vector3 position, Vector3 rotation, Color color, float scale)
    : Object(type, name, position, rotation, color, scale) {}

GameObject::GameObject(std::string type, std::string name, Vector3 position, Vector3 rotation, Color color, float scale, Model model)
    : Object(type, name, position, rotation, color, scale), model(model) {}

GameObject::~GameObject() {
    // UnloadModel frees meshes and material map arrays but not shaders/textures (those may be shared).
    if (model.meshCount > 0) {
        UnloadModel(model);
        model = Model{};
    }
}

void GameObject::draw(const RenderContext&) {
    if (!isActive) return;
    if (type == "cube") {
        DrawCube(position, scale, scale, scale, color);
        DrawCubeWires(position, scale, scale, scale, BLACK);
    } else if (type == "sphere") {
        DrawSphere(position, scale, color);
        DrawSphereWires(position, scale, 0, 1, BLACK);
    } else if (model.meshCount > 0) {
        const Matrix m = localTransform();
        for (int i = 0; i < model.meshCount; ++i) {
            DrawMesh(model.meshes[i], model.materials[model.meshMaterial[i]], m);
        }
    }
}

void GameObject::drawDepthOnly(const RenderContext& ctx) {
    if (!isActive || model.meshCount == 0 || !ctx.depthMaterial) return;
    const Matrix m = localTransform();
    for (int i = 0; i < model.meshCount; ++i) {
        DrawMesh(model.meshes[i], *ctx.depthMaterial, m);
    }
}

// ---------------------------------------------------------------------------------------------
// Chunk
// ---------------------------------------------------------------------------------------------

Chunk::~Chunk() { releaseMesh(); }

void Chunk::releaseMesh() {
    if (mesh.vaoId > 0) {
        UnloadMesh(mesh); // frees GPU buffers and the MemAlloc'd CPU arrays
    } else {
        MemFree(mesh.vertices);
        MemFree(mesh.normals);
        MemFree(mesh.colors);
        MemFree(mesh.indices);
    }
    mesh = Mesh{};
}

template <typename T>
static T* CopyToRaylib(const std::vector<T>& v) {
    T* p = (T*)MemAlloc((unsigned int)(v.size() * sizeof(T)));
    memcpy(p, v.data(), v.size() * sizeof(T));
    return p;
}

void Chunk::setMesh(const ChunkMeshData& data) {
    releaseMesh();
    mesh.vertexCount = (int)(data.positions.size() / 3);
    mesh.triangleCount = (int)(data.indices.size() / 3);
    mesh.vertices = CopyToRaylib(data.positions);
    mesh.normals = CopyToRaylib(data.normals);
    mesh.colors = CopyToRaylib(data.colors);
    mesh.indices = CopyToRaylib(data.indices);
    UploadMesh(&mesh, false);
}

Matrix Chunk::modelMatrix(const Matrix& parentTransform) const {
    const Vector3 o = voxelOffset();
    return MatrixMultiply(MatrixTranslate(o.x, o.y, o.z), parentTransform);
}

bool Chunk::visible(const RenderContext& ctx, const Matrix& model) const {
    if (!ctx.frustum || !frustumCullingEnabled) return true;
    const float h = CHUNK_SIZE * 0.5f;
    const Vector3 centre = Vector3Transform({ h, h, h }, model);
    // Bounding sphere of the chunk cube; account for the parent's scale.
    const float s = Vector3Length({ model.m0, model.m1, model.m2 });
    return ctx.frustum->sphereVisible(centre, h * 1.7320508f * s);
}

void Chunk::drawChunk(const RenderContext& ctx, const Matrix& parentTransform) {
    if (!isActive || mesh.vaoId == 0 || !ctx.terrainMaterial) return;
    const Matrix m = modelMatrix(parentTransform);
    if (!visible(ctx, m)) return;
    DrawMesh(mesh, *ctx.terrainMaterial, m);
    if (ctx.drawCalls) (*ctx.drawCalls)++;
}

void Chunk::drawChunkDepth(const RenderContext& ctx, const Matrix& parentTransform) {
    if (!isActive || mesh.vaoId == 0 || !ctx.depthMaterial) return;
    const Matrix m = modelMatrix(parentTransform);
    if (!visible(ctx, m)) return;
    DrawMesh(mesh, *ctx.depthMaterial, m);
}
