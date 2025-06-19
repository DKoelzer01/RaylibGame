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

    UpdateLightValues(shader, light, lightCount);

    lightCount++;
    return light;
}
    
void UpdateLightValues(Shader shader, Light light, int index = 0)
{
    char uniformName[32];

    // Position
    sprintf(uniformName, "lights[%d].position", index);
    int posLoc = GetShaderLocation(shader, uniformName);
    SetShaderValue(shader, posLoc, &light.position, SHADER_UNIFORM_VEC3);

    // Target
    sprintf(uniformName, "lights[%d].target", index);
    int tgtLoc = GetShaderLocation(shader, uniformName);
    SetShaderValue(shader, tgtLoc, &light.target, SHADER_UNIFORM_VEC3);

    // Color
    sprintf(uniformName, "lights[%d].color", index);
    int colLoc = GetShaderLocation(shader, uniformName);
    float color[4] = { light.color.r/255.0f, light.color.g/255.0f, light.color.b/255.0f, light.color.a/255.0f };
    SetShaderValue(shader, colLoc, color, SHADER_UNIFORM_VEC4);

    // Enabled
    sprintf(uniformName, "lights[%d].enabled", index);
    int enLoc = GetShaderLocation(shader, uniformName);
    int enabled = light.enabled ? 1 : 0;
    SetShaderValue(shader, enLoc, &enabled, SHADER_UNIFORM_INT);

    // Type
    sprintf(uniformName, "lights[%d].type", index);
    int typeLoc = GetShaderLocation(shader, uniformName);
    int type = light.type;
    SetShaderValue(shader, typeLoc, &type, SHADER_UNIFORM_INT);
}