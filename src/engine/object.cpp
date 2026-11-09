#include "object.h"
#include "textobject.h"
#include "gameobject.h"
#include <raymath.h>
#include <rlgl.h>
#include "terrainlayers.h"
#include <algorithm>
#include <cmath>
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
    ready = true;
    releaseMesh();
    if (data.positions.empty()) return; // node has no surface
    mesh.vertexCount = (int)(data.positions.size() / 3);
    mesh.triangleCount = (int)(data.indices.size() / 3);
    mesh.vertices = CopyToRaylib(data.positions);
    mesh.normals = CopyToRaylib(data.normals);
    mesh.colors = CopyToRaylib(data.colors);
    if (!data.texcoords.empty()) mesh.texcoords = CopyToRaylib(data.texcoords);    // material weights 4-5
    if (!data.texcoords2.empty()) mesh.texcoords2 = CopyToRaylib(data.texcoords2); // material weights 6-7
    mesh.indices = CopyToRaylib(data.indices);
    UploadMesh(&mesh, false);
}

bool Chunk::visible(const RenderContext& ctx, const Matrix& model) const {
    if (!ctx.frustum || !frustumCullingEnabled) return true;
    const float h = key.size() * 0.5f * (CHUNK_SIZE + 1.0f) / CHUNK_SIZE; // includes the overlap cell
    const Vector3 centre = Vector3Transform({ h, h, h }, model);
    // Bounding sphere of the chunk cube; account for the parent's scale.
    const float s = Vector3Length({ model.m0, model.m1, model.m2 });
    return ctx.frustum->sphereVisible(centre, h * 1.7320508f * s);
}

// raylib allocates MAX_MATERIAL_MAPS (config.h, default 12) maps per Material; not exported by raylib.h.
static constexpr int kMaterialMaps = 12;

ChunkBatch::ChunkBatch(const Material& m, const TextureArray* arr, int arrLoc) : material(m), layers(arr && arr->id ? arr : nullptr) {
    view = rlGetMatrixModelview();
    projection = rlGetMatrixProjection();
    viewProjection = MatrixMultiply(view, projection);
    const int* locs = material.shader.locs;
    rlEnableShader(material.shader.id);
    if (locs[SHADER_LOC_COLOR_DIFFUSE] != -1) {
        const Color c = material.maps[MATERIAL_MAP_DIFFUSE].color;
        const float v[4] = { c.r / 255.0f, c.g / 255.0f, c.b / 255.0f, c.a / 255.0f };
        rlSetUniform(locs[SHADER_LOC_COLOR_DIFFUSE], v, SHADER_UNIFORM_VEC4, 1);
    }
    if (locs[SHADER_LOC_MATRIX_VIEW] != -1) rlSetUniformMatrix(locs[SHADER_LOC_MATRIX_VIEW], view);
    if (locs[SHADER_LOC_MATRIX_PROJECTION] != -1) rlSetUniformMatrix(locs[SHADER_LOC_MATRIX_PROJECTION], projection);
    for (int i = 0; i < kMaterialMaps; i++) {
        if (material.maps[i].texture.id == 0) continue;
        rlActiveTextureSlot(i);
        if (i == MATERIAL_MAP_IRRADIANCE || i == MATERIAL_MAP_PREFILTER || i == MATERIAL_MAP_CUBEMAP) rlEnableTextureCubemap(material.maps[i].texture.id);
        else rlEnableTexture(material.maps[i].texture.id);
        rlSetUniform(locs[SHADER_LOC_MAP_DIFFUSE + i], &i, SHADER_UNIFORM_INT, 1);
    }
    if (layers) {
        BindTextureArray(kLayersUnit, *layers);
        const int unit = kLayersUnit;
        if (arrLoc >= 0) rlSetUniform(arrLoc, &unit, SHADER_UNIFORM_INT, 1);
    }
}

ChunkBatch::~ChunkBatch() {
    if (layers) UnbindTextureArray(kLayersUnit);
    for (int i = 0; i < kMaterialMaps; i++) {
        if (material.maps[i].texture.id == 0) continue;
        rlActiveTextureSlot(i);
        if (i == MATERIAL_MAP_IRRADIANCE || i == MATERIAL_MAP_PREFILTER || i == MATERIAL_MAP_CUBEMAP) rlDisableTextureCubemap();
        else rlDisableTexture();
    }
    rlDisableVertexArray();
    rlDisableVertexBuffer();
    rlDisableVertexBufferElement();
    rlDisableShader();
    rlSetMatrixModelview(view);
    rlSetMatrixProjection(projection);
}

void ChunkBatch::setUniform(int loc, const void* value, int uniformType) {
    if (loc >= 0) rlSetUniform(loc, value, uniformType, 1);
}

void ChunkBatch::draw(const Mesh& mesh, const Matrix& model) {
    const int* locs = material.shader.locs;
    if (locs[SHADER_LOC_MATRIX_MODEL] != -1) rlSetUniformMatrix(locs[SHADER_LOC_MATRIX_MODEL], model);
    if (locs[SHADER_LOC_MATRIX_MVP] != -1) rlSetUniformMatrix(locs[SHADER_LOC_MATRIX_MVP], MatrixMultiply(model, viewProjection));
    if (!rlEnableVertexArray(mesh.vaoId)) return; // chunk meshes are always uploaded with a VAO
    if (mesh.indices != nullptr) rlDrawVertexArrayElements(0, mesh.triangleCount * 3, 0);
    else rlDrawVertexArray(0, mesh.vertexCount);
}

void Chunk::drawChunk(const RenderContext& ctx, const Matrix& m) {
    if (!isActive || mesh.vaoId == 0 || !ctx.terrainMaterial) return;
    if (!visible(ctx, m)) return;
    // Per-chunk uniforms. Both stay small far from a body's centre: the texturing origin is the
    // chunk origin wrapped to the texture's repeat (computed in double, so neighbouring chunks still
    // line up), and the clip bound is relative to the chunk.
    const Vector3 o = key.origin();
    const double period = 1.0 / std::max(ctx.texScale, 1e-6f);
    auto wrap = [period](float v) { return (float)(fmod(fmod((double)v, period) + period, period)); };
    const Vector3 texOrigin = { wrap(o.x), wrap(o.y), wrap(o.z) };
    const Vector3 clip = { clipMax.x - o.x, clipMax.y - o.y, clipMax.z - o.z }; // 1e30 stays "no clip"
    if (ctx.batch) {
        ctx.batch->setUniform(ctx.chunkOriginLoc, &texOrigin, SHADER_UNIFORM_VEC3);
        ctx.batch->setUniform(ctx.chunkClipLoc, &clip, SHADER_UNIFORM_VEC3);
        ctx.batch->draw(mesh, m);
    } else {
        // DrawMesh draws immediately, so per-draw uniforms set just before it are safe.
        if (ctx.chunkOriginLoc >= 0) SetShaderValue(ctx.terrainMaterial->shader, ctx.chunkOriginLoc, &texOrigin, SHADER_UNIFORM_VEC3);
        if (ctx.chunkClipLoc >= 0) SetShaderValue(ctx.terrainMaterial->shader, ctx.chunkClipLoc, &clip, SHADER_UNIFORM_VEC3);
        DrawMesh(mesh, *ctx.terrainMaterial, m);
    }
    if (ctx.drawCalls) (*ctx.drawCalls)++;
}

void Chunk::drawChunkDepth(const RenderContext& ctx, const Matrix& m) {
    if (!isActive || mesh.vaoId == 0 || !ctx.depthMaterial) return;
    if (!visible(ctx, m)) return;
    if (ctx.batch) ctx.batch->draw(mesh, m);
    else DrawMesh(mesh, *ctx.depthMaterial, m);
}
