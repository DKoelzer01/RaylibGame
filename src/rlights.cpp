// rlights.c - Raylib lights implementation (public domain)
// https://github.com/raysan5/raylib/blob/master/examples/shaders/rlights.h
#include "rlights.h"

int lightCount = 0;

Light CreateLight(int type, Vector3 position, Vector3 target, Color color, Shader shader)
{
    Light light = { 0 };
    light.enabled = true;
    light.type = type;
    light.position = position;
    light.target = target;
    light.color = color;

    // Assign a unique index to the light
    light.enabledLoc = GetShaderLocation(shader, TextFormat("lights[%i].enabled", lightCount));
    light.typeLoc = GetShaderLocation(shader, TextFormat("lights[%i].type", lightCount));
    light.positionLoc = GetShaderLocation(shader, TextFormat("lights[%i].position", lightCount));
    light.targetLoc = GetShaderLocation(shader, TextFormat("lights[%i].target", lightCount));
    light.colorLoc = GetShaderLocation(shader, TextFormat("lights[%i].color", lightCount));

    UpdateLightValues(shader, light);

    lightCount++;
    return light;
}

void UpdateLightValues(Shader shader, Light light)
{
    // Send to shader
    SetShaderValue(shader, light.enabledLoc, &light.enabled, SHADER_UNIFORM_INT);
    SetShaderValue(shader, light.typeLoc, &light.type, SHADER_UNIFORM_INT);
    float pos[3] = { light.position.x, light.position.y, light.position.z };
    SetShaderValue(shader, light.positionLoc, pos, SHADER_UNIFORM_VEC3);
    float tgt[3] = { light.target.x, light.target.y, light.target.z };
    SetShaderValue(shader, light.targetLoc, tgt, SHADER_UNIFORM_VEC3);
    float col[4] = { (float)light.color.r/255.0f, (float)light.color.g/255.0f, (float)light.color.b/255.0f, (float)light.color.a/255.0f };
    SetShaderValue(shader, light.colorLoc, col, SHADER_UNIFORM_VEC4);
}