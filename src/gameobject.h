#ifndef gameobject_h
#define gameobject_h
#include "object.h"
#include "worldtypes.h"
#include "utils/logger.h"

#pragma once

class GameObject: public Object {
    public:
        Model model;
        Texture2D texture;
        std::string modelPath; // Path to the model file
        std::string texturePath; // Path to the texture file
        GameObject(std::string type,std::string name, Vector3 position, Vector3 rotation, Color color, float scale);
        GameObject(std::string type,std::string name, Vector3 position, Vector3 rotation, Color color, float scale, Model model);
        virtual ~GameObject();

        void draw() override;
};

class ChunkObject: public GameObject {
    public:
        Chunk chunk; // Store chunk data
        ChunkObject(std::string type,std::string name, Vector3 position, Vector3 rotation, Color color, float scale, Chunk chunk);
        virtual ~ChunkObject();

        void draw() override;
};

#endif