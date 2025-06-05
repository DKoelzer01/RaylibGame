#ifndef gameobject_h
#define gameobject_h
#include "object.h"

#pragma once

class GameObject: public Object {
    public:
        GameObject(std::string type, Vector3 position, Vector3 rotation, Color color, float scale);
        void draw() override;
};

#endif