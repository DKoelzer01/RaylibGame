#ifndef OBJECT_H
#define OBJECT_H

#include <string>
#include <vector>
#include <memory>
#include <raylib.h>
#include <iostream>

#pragma once
class Object {
public:
    std::string type; // Type of the object (e.g., "cube", "sphere")
    std::string name; // Name of the object
    Vector3 position; // Position of the object
    Vector3 rotation; // Rotation of the object
    Color color; // Color of the object
    float scale; // Scale of the object

    Object* parent;
    std::vector<std::unique_ptr<Object>> children;

    bool isActive; // Whether the object is active in the scene

    Object(std::string type, std::string name, Vector3 position, Vector3 rotation, Color color, float scale);
    virtual void draw();
    virtual void drawDepthOnly(const Matrix& lightSpaceMatrix, Shader* depthShader);
    virtual ~Object(); // Make destructor virtual for safe polymorphic deletion

    // Disable copying
    Object(const Object&) = delete;
    Object& operator=(const Object&) = delete;

    // Allow moving
    Object(Object&&) = default;
    Object& operator=(Object&&) = default;
};

#endif