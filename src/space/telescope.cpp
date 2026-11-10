#include "telescope.h"
#include "inventory.h"
#include "player.h"
#include "world.h"
#include "worldtypes.h"
#include "audio.h"

#include <algorithm>
#include <cmath>

TelescopeState telescope;

void Telescope_SetActive(Scene& world, bool on) {
    if (on == telescope.active) return;
    telescope.active = on;
    if (on) telescope.baseFov = world.camera.fovy;
    else {
        world.camera.fovy = telescope.baseFov;
        player.lookScale = 1.0f;
        SetLodZoom(world, 1.0f, { 0, 0, 1 }, 0.0f);
    }
}

void Telescope_Update(Scene& world, bool inputEnabled, bool usable) {
    const bool has = inventory.count("telescope") > 0;
    if (!usable || !has) { Telescope_SetActive(world, false); return; }
    if (inputEnabled && IsKeyPressed(KEY_T)) {
        Telescope_SetActive(world, !telescope.active);
        Audio_Play(Sfx::Build, 0.25f, telescope.active ? 1.6f : 1.3f);
    }
    if (!telescope.active) return;
    if (inputEnabled) {
        const float wheel = GetMouseWheelMove();
        if (wheel != 0.0f) telescope.zoom *= wheel > 0.0f ? 1.25f : 1.0f / 1.25f;
    }
    telescope.zoom = std::clamp(telescope.zoom, telescope.minZoom, telescope.maxZoom);
    // Zoom by the tangent of the half-angle, so 8x really magnifies 8 times.
    const float baseHalf = telescope.baseFov * 0.5f * DEG2RAD;
    const float half = atanf(tanf(baseHalf) / telescope.zoom);
    world.camera.fovy = 2.0f * half * RAD2DEG;
    player.lookScale = 1.0f / telescope.zoom;
    // LOD cone to the screen corners.
    const float aspect = (float)GetScreenWidth() / (float)std::max(1, GetScreenHeight());
    const float corner = atanf(tanf(half) * sqrtf(1.0f + aspect * aspect));
    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    SetLodZoom(world, telescope.zoom, dir, corner);
}

void Telescope_DrawHUD(Scene& world) {
    if (!telescope.active) return;
    const int w = GetScreenWidth(), h = GetScreenHeight();
    const Vector2 c = { w * 0.5f, h * 0.5f };
    const float r = h * 0.47f;
    // Eyepiece: black outside a circle, a thin reticle inside.
    DrawRing(c, r, sqrtf((float)(w * w + h * h)), 0.0f, 360.0f, 96, BLACK);
    DrawCircleLinesV(c, r, Fade(RAYWHITE, 0.25f));
    DrawLineV({ c.x - r, c.y }, { c.x - 12, c.y }, Fade(RAYWHITE, 0.35f));
    DrawLineV({ c.x + 12, c.y }, { c.x + r, c.y }, Fade(RAYWHITE, 0.35f));
    DrawLineV({ c.x, c.y - r }, { c.x, c.y - 12 }, Fade(RAYWHITE, 0.35f));
    DrawLineV({ c.x, c.y + 12 }, { c.x, c.y + r }, Fade(RAYWHITE, 0.35f));

    // What the reticle is on: the nearest body whose disc covers the centre of the view.
    const Vector3 eye = world.camera.position;
    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, eye));
    const Planetoid* target = nullptr;
    float targetDist = 0.0f;
    for (auto& obj : world.objects) {
        const auto* p = dynamic_cast<const Planetoid*>(obj.get());
        if (!p || !p->field) continue;
        const Vector3 to = Vector3Subtract(p->position, eye);
        const float d = Vector3Length(to);
        const float rad = p->field->maxRadius() * p->scale;
        if (d <= rad) continue; // standing on it
        if (Vector3DotProduct(to, dir) / d < cosf(asinf(std::min(1.0f, rad / d)))) continue;
        if (!target || d < targetDist) { target = p; targetDist = d; }
    }
    int y = (int)(c.y + r) - 64;
    DrawText(TextFormat("x%.0f", telescope.zoom), (int)c.x - MeasureText(TextFormat("x%.0f", telescope.zoom), 24) / 2, y, 24, RAYWHITE);
    if (target) {
        const char* label = TextFormat("%s %s  %.0f km", BodyKindName(target->kind), target->name.c_str(), targetDist / 1000.0f);
        DrawText(label, (int)c.x - MeasureText(label, 18) / 2, y + 30, 18, Fade(SKYBLUE, 0.9f));
    }
    DrawText("wheel zoom   T lower", 10, h - 30, 20, Fade(RAYWHITE, 0.7f));
}
