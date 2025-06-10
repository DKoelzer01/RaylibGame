#include "object.h"
#include "textobject.h"
#include "gameobject.h"

Object::Object(std::string type,std::string name, Vector3 position, Vector3 rotation, Color color, float scale)
    : type(type), name(name), position(position), rotation(rotation), color(color), scale(scale) {}

void Object::draw() {
}

Object::~Object() {
    //std::cout << "Object of type " << type << " destroyed." << std::endl;
}

TextObject::TextObject(std::string text, int fontSize) 
    : Object("text","text", {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, WHITE, 1.0f), text(text), fontSize(fontSize) {}

TextObject::TextObject(std::string text, Vector3 position, Vector3 rotation, Color color, float scale, int fontSize)
    : Object("text","text", position, rotation, color, scale), text(text), fontSize(fontSize) {}

void TextObject::draw() {
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
    this->model = chunk.model;
}

ChunkObject::~ChunkObject() {
    // Do NOT unload chunk.model here, as it is a copy and may be managed elsewhere.
    // If you ever manage chunk.model separately, add safe unload logic here.
}

void GameObject::draw() {
    if (!isActive) return; // Skip drawing if the object is not active
    if (type == "cube") {
        DrawCube(position, scale, scale, scale, color);
        DrawCubeWires(position, scale, scale, scale, BLACK);
    } else if (type == "sphere") {
        DrawSphere(position, scale, color);
        DrawSphereWires(position, scale,0,1, BLACK);
    } else  {
        DrawModel(model, position, scale, color);
        DrawModelWires(model, position, scale, BLACK);
    }
}

void ChunkObject::draw() {
    // logger.logf("Drawing ChunkObject: %s at position (%f, %f, %f) active:%d\n", name.c_str(), position.x, position.y, position.z, isActive);
    if (!isActive) return; // Skip drawing if the object is not active
    DrawModel(chunk.model, position, scale, color);
    // DrawModelWires(chunk.model, position, scale, BLACK);
    // DrawCube({position.x+16,position.y+16,position.z+16},32,32,32, {200,200,200,50}); // Draw a cube at the chunk position for visualization
}
