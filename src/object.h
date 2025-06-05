#ifndef OBJECT_H
#define OBJECT_H

#include <string>
#include <raylib.h>
#include <iostream>

#pragma once
class Object {
public:
    std::string type; // Type of the object (e.g., "cube", "sphere")
    Vector3 position; // Position of the object
    Vector3 rotation; // Rotation of the object
    Color color; // Color of the object
    float scale; // Scale of the object

    Object(std::string type, Vector3 position, Vector3 rotation, Color color, float scale);
    virtual void draw();
    virtual ~Object();
};

#endif