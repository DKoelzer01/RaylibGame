#include "object.h"
#include "textobject.h"
#include "gameobject.h"

Object::Object(std::string type, Vector3 position, Vector3 rotation, Color color, float scale)
    : type(type), position(position), rotation(rotation), color(color), scale(scale) {}

void Object::draw() {
}

Object::~Object() {
    std::cout << "Object of type " << type << " destroyed." << std::endl;
}

TextObject::TextObject(std::string text, int fontSize) 
    : Object("text", {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, WHITE, 1.0f), text(text), fontSize(fontSize) {}

TextObject::TextObject(std::string text, Vector3 position, Vector3 rotation, Color color, float scale, int fontSize)
    : Object("text", position, rotation, color, scale), text(text), fontSize(fontSize) {}

void TextObject::draw() {
    DrawText(text.c_str(), static_cast<int>(position.x), static_cast<int>(position.y), fontSize, color);
}

GameObject::GameObject(std::string type, Vector3 position, Vector3 rotation, Color color, float scale)
    : Object(type, position, rotation, color, scale) {}

void GameObject::draw() {
    if (type == "cube") {
        DrawCube(position, scale, scale, scale, color);
    } else if (type == "sphere") {
        DrawSphere(position, scale, color);
    } else {
        std::cout << "Unknown object type: " << type << std::endl;
    }
}
