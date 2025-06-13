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

Chunk::~Chunk() {
    // Unload mesh and model resources
    if (mesh.vertices) delete[] mesh.vertices;
    if (mesh.indices) delete[] mesh.indices;
    if (mesh.normals) delete[] mesh.normals;
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

void Chunk::draw(Shader* lightingShader) {
    if (!isActive) return;
    Matrix matModel = MatrixIdentity();
    matModel = MatrixMultiply(matModel, MatrixScale(scale, scale, scale));
    matModel = MatrixMultiply(matModel, MatrixTranslate(position.x, position.y, position.z));

    if (mesh.vertexCount > 0 && model.meshCount > 0 && model.materialCount > 0 && model.materials != nullptr) {
        for (int i = 0; i < model.materialCount; i++) {
            model.materials[i].shader = *lightingShader;
        }
        BeginShaderMode(model.materials[0].shader);
        DrawMesh(mesh, model.materials[0], matModel);
        EndShaderMode();
    }

    // Draw chunk boundary for debugging
    DrawCubeWires((Vector3){
        position.x + CHUNK_SIZE/2.0f,
        position.y + CHUNK_SIZE/2.0f,
        position.z + CHUNK_SIZE/2.0f
    }, CHUNK_SIZE, CHUNK_SIZE, CHUNK_SIZE, GREEN);
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

void Chunk::drawDepthOnly(const Matrix& lightSpaceMatrix, Shader* depthShader) {
    if (!isActive) return;
    if (model.meshCount == 0 || model.meshes == nullptr) return;
    int lightSpaceLoc = GetShaderLocation(*depthShader, "lightSpaceMatrix");
    SetShaderValueMatrix(*depthShader, lightSpaceLoc, lightSpaceMatrix);
    int modelLoc = GetShaderLocation(*depthShader, "model");
    Matrix modelMat = MatrixMultiply(MatrixTranslate(position.x, position.y, position.z), MatrixScale(scale, scale, scale));
    SetShaderValueMatrix(*depthShader, modelLoc, modelMat);
    BeginShaderMode(*depthShader);
    model.materials[0].shader = *depthShader; // Ensure correct shader for chunk model
    DrawMesh(mesh, model.materials[0], MatrixIdentity());
    EndShaderMode();
}
