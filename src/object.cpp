#include "object.h"
#include "textobject.h"
#include "gameobject.h"
#include <raymath.h>

Object::Object(std::string type,std::string name, Vector3 position, Vector3 rotation, Color color, float scale)
    : type(type), name(name), position(position), rotation(rotation), color(color), scale(scale) {}

void Object::draw(Shader* lightingShader) {
    if (!isActive) return; // Skip drawing if the object is not active
    // Base objects do not render anything by default
    // Derived classes should implement their own draw logic
    std::cout << "Drawing base Object: " << name << " of type " << type << std::endl;
}

void Object::drawDepthOnly(const Matrix& lightSpaceMatrix, Shader* depthShader) {
    // Base objects do not render depth
}

Object::~Object() {
    //std::cout << "Object of type " << type << " destroyed." << std::endl;
}

TextObject::TextObject(std::string text, int fontSize) 
    : Object("text","text", {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, WHITE, 1.0f), text(text), fontSize(fontSize) {}

TextObject::TextObject(std::string text, Vector3 position, Vector3 rotation, Color color, float scale, int fontSize)
    : Object("text","text", position, rotation, color, scale), text(text), fontSize(fontSize) {}

void TextObject::draw(Shader* lightingShader) {
    DrawText(text.c_str(), static_cast<int>(position.x), static_cast<int>(position.y), fontSize, color);
}

GameObject::GameObject(std::string type,std::string name, Vector3 position, Vector3 rotation, Color color, float scale)
    : Object(type, name, position, rotation, color, scale) {}

GameObject::GameObject(std::string type, std::string name, Vector3 position, Vector3 rotation, Color color, float scale, Model model)
    : Object(type, name, position, rotation, color, scale), model(model) {}

GameObject::~GameObject() {
    // Only unload if this is the only owner and model is valid
    if (model.meshCount > 0) {
        UnloadModel(model);
        model.meshCount = 0; // Prevent double-free if destructor called again
    }
}

ChunkObject::ChunkObject(std::string type, std::string name, Vector3 position, Vector3 rotation, Color color, float scale, Chunk chunk)
    : GameObject(type, name, position, rotation, color, scale), chunk(chunk) {
    // Do NOT assign this->model = chunk.model; Each chunk must have its own unique model/mesh.
}

ChunkObject::~ChunkObject() {
    // Do NOT unload chunk.model here, as it is a copy and may be managed elsewhere.
    // If you ever manage chunk.model separately, add safe unload logic here.
}

void GameObject::draw(Shader* lightingShader) {
    if (!isActive) return; // Skip drawing if the object is not active
    if (type == "cube") {
        DrawCube(position, scale, scale, scale, color);
        DrawCubeWires(position, scale, scale, scale, BLACK);
    } else if (type == "sphere") {
        DrawSphere(position, scale, color);
        DrawSphereWires(position, scale,0,1, BLACK);
    } else  {
        Matrix matModel = MatrixIdentity();
        matModel = MatrixMultiply(matModel, MatrixScale(scale, scale, scale));
        matModel = MatrixMultiply(matModel, MatrixRotateXYZ((Vector3){ DEG2RAD*rotation.x, DEG2RAD*rotation.y, DEG2RAD*rotation.z }));
        matModel = MatrixMultiply(matModel, MatrixTranslate(position.x, position.y, position.z));

        int matModelLoc = GetShaderLocation(*lightingShader, "matModel");
        SetShaderValueMatrix(*lightingShader, matModelLoc, matModel);
        DrawModel(model, {0,0,0}, 1.0f, WHITE);
        DrawModelWires(model, {0,0,0}, 1.0f, WHITE);
    }
}

void ChunkObject::draw(Shader* lightingShader) {
    if (!isActive) return;
    Matrix matModel = MatrixIdentity();
    matModel = MatrixMultiply(matModel, MatrixScale(scale, scale, scale));
    // Only apply chunk-local position, do NOT add parent->position (planetoid origin)
    matModel = MatrixMultiply(matModel, MatrixTranslate(position.x, position.y, position.z));
    // Verbose logging for debugging chunk transform
    // logger.logf("Drawing ChunkObject '%s' at local position (%f, %f, %f) with scale %.2f\n",
    //     name.c_str(), position.x, position.y, position.z, scale);
    // logger.logf("ChunkObject '%s' model matrix:\n", name.c_str());
    // logger.logf("  Scale: (%f, %f, %f)\n", scale, scale, scale);
    // logger.logf("  Rotation: (%f, %f, %f)\n", rotation.x, rotation.y, rotation.z);
    // logger.logf("  Position: (%f, %f, %f)\n", position.x, position.y, position.z);
    // logger.logf("  Model Matrix: \n");
    // logger.logf("  [%f, %f, %f, %f]\n", matModel.m0, matModel.m1, matModel.m2, matModel.m3);
    // logger.logf("  [%f, %f, %f, %f]\n", matModel.m4, matModel.m5, matModel.m6, matModel.m7);
    // logger.logf("  [%f, %f, %f, %f]\n", matModel.m8, matModel.m9, matModel.m10, matModel.m11);
    // logger.logf("  [%f, %f, %f, %f]\n", matModel.m12, matModel.m13, matModel.m14, matModel.m15);

    // Diagnostic: print first 10 mesh vertices for this chunk
    // if (chunk.vertices.size() > 0) {
    //     logger.logf("ChunkObject '%s' first 10 mesh vertices (chunk-local):\n", name.c_str());
    //     for (size_t i = 0; i < chunk.vertices.size() && i < 10; ++i) {
    //         const Vector3& v = chunk.vertices[i];
    //         logger.logf("  v[%zu] = (%.6f, %.6f, %.6f)\n", i, v.x, v.y, v.z);
    //     }
    // }
    // // Diagnostic: print first 10 mesh indices for this chunk
    // if (chunk.indices.size() > 0) {
    //     logger.logf("ChunkObject '%s' first 10 mesh indices:\n", name.c_str());
    //     for (size_t i = 0; i < chunk.indices.size() && i < 10; ++i) {
    //         logger.logf("  idx[%zu] = %d\n", i, chunk.indices[i]);
    //     }
    // }

    // Diagnostic: print mesh vertices on +X, +Y, +Z chunk faces (borders)
    // int borderCount = 0;
    // float epsilon = 1e-4f;
    // // Assume CHUNK_SIZE is available globally or via chunk struct
    // extern const int CHUNK_SIZE;
    // for (size_t i = 0; i < chunk.vertices.size(); ++i) {
    //     const Vector3& v = chunk.vertices[i];
    //     // Check if vertex is on +X, +Y, or +Z face (within epsilon)
    //     if (fabs(v.x - CHUNK_SIZE) < epsilon || fabs(v.y - CHUNK_SIZE) < epsilon || fabs(v.z - CHUNK_SIZE) < epsilon) {
    //         logger.logf("ChunkObject '%s' border vertex idx %zu: (%.6f, %.6f, %.6f)\n", name.c_str(), i, v.x, v.y, v.z);
    //         borderCount++;
    //         if (borderCount >= 20) break; // Limit output
    //     }
    // }
    // if (borderCount == 0) {
    //     logger.logf("ChunkObject '%s' has no mesh vertices on +X/+Y/+Z border faces.\n", name.c_str());
    // }

    if (chunk.mesh.vertexCount > 0 && chunk.model.meshCount > 0 && chunk.model.materialCount > 0 && chunk.model.materials != nullptr) {
        for (int i = 0; i < chunk.model.materialCount; i++) {
            chunk.model.materials[i].shader = *lightingShader;
        }
        BeginShaderMode(chunk.model.materials[0].shader);
        DrawMesh(chunk.mesh, chunk.model.materials[0], matModel);
        EndShaderMode();
    }
}

void GameObject::drawDepthOnly(const Matrix& lightSpaceMatrix, Shader* depthShader) {
    if (!isActive) return;
    int lightSpaceLoc = GetShaderLocation(*depthShader, "lightSpaceMatrix");
    SetShaderValueMatrix(*depthShader, lightSpaceLoc, lightSpaceMatrix);
    Matrix matModel = MatrixIdentity();
    matModel = MatrixMultiply(matModel, MatrixScale(scale, scale, scale));
    matModel = MatrixMultiply(matModel, MatrixRotateXYZ((Vector3){ DEG2RAD*rotation.x, DEG2RAD*rotation.y, DEG2RAD*rotation.z }));
    matModel = MatrixMultiply(matModel, MatrixTranslate(position.x, position.y, position.z));
    int matModelLoc = GetShaderLocation(*depthShader, "matModel");
    SetShaderValueMatrix(*depthShader, matModelLoc, matModel);
    BeginShaderMode(*depthShader);
    DrawModel(model, {0,0,0}, 1.0f, WHITE); // Use identity transform, model matrix handles all
    EndShaderMode();
}

void ChunkObject::drawDepthOnly(const Matrix& lightSpaceMatrix, Shader* depthShader) {
    if (!isActive) return;
    if (chunk.model.meshCount == 0 || chunk.model.meshes == nullptr) return;
    int lightSpaceLoc = GetShaderLocation(*depthShader, "lightSpaceMatrix");
    SetShaderValueMatrix(*depthShader, lightSpaceLoc, lightSpaceMatrix);
    int modelLoc = GetShaderLocation(*depthShader, "model");
    Matrix modelMat = MatrixMultiply(MatrixTranslate(position.x, position.y, position.z), MatrixScale(scale, scale, scale));
    SetShaderValueMatrix(*depthShader, modelLoc, modelMat);
    BeginShaderMode(*depthShader);
    chunk.model.materials[0].shader = *depthShader; // Ensure correct shader for chunk model
    DrawMesh(chunk.mesh, chunk.model.materials[0], MatrixIdentity());
    EndShaderMode();
}
