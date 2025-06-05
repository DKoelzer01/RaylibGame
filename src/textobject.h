#ifndef TEXTOBJECT_H
#define TEXTOBJECT_H

#include "object.h"

#pragma once
class TextObject: public Object     {
    public:
        std::string text; // Text to display
        int fontSize; // Font size of the text
        TextObject(std::string text, int fontSize);
        TextObject(std::string text, Vector3 position, Vector3 rotation, Color color, float scale, int fontSize);
        virtual void draw() override;

};

#endif