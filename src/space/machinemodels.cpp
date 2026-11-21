#include "machinemodels.h"
#include "rotationmodels.h"
#include "factory.h"
#include "trapping.h"
#include "instancebatch.h"
#include "plots.h"
#include "scene.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
#include <unordered_map>

namespace {

// The machine's frame: x along the plot's east, y up, z along its north, origin at the centre of the
// footprint's floor; w, h, d its size (m); front = the part's facing in that frame (unit, horizontal).
struct Frame {
    float w, h, d; Vector3 front; float yawDeg;
    // Free space over the casing and under its intake corner (m, whole cells; 1e9 = open), from the plot's other parts
    // (#221): details that would reach into a neighbouring block are dropped (stacked machines, no z-fighting).
    float roomAbove = 1e9f, roomBelow = 1e9f;
};

Color Shade(Color c, float light, float m = 1.0f) {
    return { (unsigned char)std::clamp(c.r * light * m, 0.0f, 255.0f), (unsigned char)std::clamp(c.g * light * m, 0.0f, 255.0f), (unsigned char)std::clamp(c.b * light * m, 0.0f, 255.0f), c.a };
}

struct Painter {
    float light;
    void cube(Vector3 c, Vector3 s, Color col, float m = 1.0f) const { IB_Cube(c, s, Shade(col, 1.0f, m)); } // lit by the batch shader
    void glow(Vector3 c, Vector3 s, Color col) const { IB_Cube(c, s, col, true); } // unlit (lit doors, hot coils)
    void cyl(Vector3 base, float rTop, float rBottom, float h, Color col, float m = 1.0f, int = 10) const { IB_Cylinder(base, rTop, rBottom, h, Shade(col, 1.0f, m)); }
};

constexpr float kTorchStick = 0.05f; // a placed torch's stick centre from its cell's low x/z corner (plots.cpp: 0.1 m thick, #558)
const Color kSteel = { 140, 144, 150, 255 }, kDarkSteel = { 70, 72, 78, 255 }, kBrick = { 150, 80, 60, 255 }, kCopper = { 190, 110, 60, 255 };
const Color kGlowOrange = { 255, 150, 50, 255 }, kGlass = { 120, 200, 230, 255 }, kPaint = { 200, 170, 50, 255 };

void Chimney(const Painter& p, const Frame& f, float x, float z, float r, float h, Color c) {
    p.cyl({ x, f.h, z }, r, r * 1.15f, h, c);
    p.cyl({ x, f.h + h, z }, r * 1.2f, r * 1.2f, 0.06f, kDarkSteel);
}

// Each machine's chimneys (frame metres: x east, z north from the floor centre; r, h; colour). One table for the
// model and the smoke (Machine_SmokePoints, #220).
struct Stack { float x, z, r, h; Color c; };
int StacksAll(PartType type, const Frame& f, Stack* out) {
    const float x0 = -f.w * 0.5f, z0 = -f.d * 0.5f;
    switch (type) {
    case PartType::Kiln:           out[0] = { x0 + 0.2f, z0 + 0.2f, 0.1f, 0.6f, kBrick }; return 1;
    case PartType::CokeOven:       out[0] = { x0 + 0.2f, z0 + 0.2f, 0.09f, 0.7f, kBrick }; out[1] = { -x0 - 0.2f, z0 + 0.2f, 0.09f, 0.7f, kBrick }; return 2;
    case PartType::BlastFurnace:   out[0] = { 0, 0, 0.16f, 1.2f, kBrick }; return 1;
    case PartType::Boiler:         out[0] = { x0 + 0.15f, z0 + 0.15f, 0.07f, 0.8f, kDarkSteel }; return 1;
    case PartType::SteamMacerator: out[0] = { x0 + 0.15f, z0 + 0.15f, 0.06f, 0.5f, kDarkSteel }; return 1;
    default: return 0;
    }
}
// The chimneys with room over the casing (#221): one under another block is left out (no model, smoke or hit box).
int Stacks(PartType type, const Frame& f, Stack* out) {
    int n = StacksAll(type, f, out), k = 0;
    for (int i = 0; i < n; ++i) if (out[i].h + 0.06f <= f.roomAbove) out[k++] = out[i];
    return k;
}

// The door on the front face (front = +z after the yaw): dark when idle, glowing when working.
void Door(const Painter& p, const Frame& f, bool working, float wFrac = 0.45f, float hFrac = 0.45f) {
    rlPushMatrix();
    rlRotatef(f.yawDeg, 0, 1, 0);
    const float depth = (fabsf(f.front.x) > 0.5f ? f.w : f.d) * 0.5f;
    const Vector3 c = { 0, f.h * 0.35f, depth + 0.01f }, s = { f.w * wFrac, f.h * hFrac, 0.03f };
    p.cube({ c.x, c.y, c.z - 0.005f }, { s.x + 0.08f, s.y + 0.08f, 0.02f }, kDarkSteel);
    if (working) p.glow(c, s, Color{ (unsigned char)(220 + 35 * sinf((float)GameClock() * 7.0f)), 120, 40, 255 });
    else p.cube(c, s, { 30, 28, 26, 255 });
    rlPopMatrix();
}

// How far each machine's top details (not its chimneys or pump pipe) reach over the casing (m), as DrawMachine draws them.
float TopReach(PartType type) {
    switch (type) {
    case PartType::BlastFurnace: return 0.6f;
    case PartType::Furnace: return 0.04f;
    case PartType::ArcFurnace: return 0.55f;
    case PartType::AlloySmelter: return 0.2f;
    case PartType::Boiler: return 0.46f;
    case PartType::ForgeHammer: return 0.86f;
    case PartType::Macerator: case PartType::SteamMacerator: return 0.39f;
    case PartType::Generator: case PartType::MvGenerator: case PartType::HvGenerator: return 0.08f;
    case PartType::Battery: return 0.12f;
    case PartType::OreWasher: return 0.47f;
    case PartType::Centrifuge: return 0.3f;
    case PartType::HandPump: return 0.71f;
    case PartType::WaterPump: return 0.45f;
    case PartType::Electrolyzer: return 0.5f;
    case PartType::ChemicalReactor: return 0.6f;
    case PartType::Refinery: return 1.85f;
    case PartType::Miner: return 0.93f;
    case PartType::PrecisionAssembler: return 0.95f;
    case PartType::LvLab: case PartType::MvLab: case PartType::HvLab: case PartType::SteamLab: return 0.42f; // (#1193) the dish on its mast
    case PartType::ResearchTable: return 0.25f; case PartType::ManaSiphon: return 0.7f; case PartType::ManaDynamo: return 0.3f; // (#1252)
    case PartType::FusionReactor: return 0.29f;
    case PartType::MassDriver: return 1.3f;
    case PartType::Torch: return 0.25f;
    default: return 0.0f; // Kiln, CokeOven, Chest: nothing over the top but chimneys
    }
}

void DrawMachine(PartType type, const Frame& f, const Painter& p, bool working) {
    const float t = (float)GameClock();
    const float spin = working ? t : 0.0f;
    const float x0 = -f.w * 0.5f, z0 = -f.d * 0.5f, top = f.h;
    {   // chimneys (one table with the smoke, #220)
        Stack st[4];
        for (int i = 0, n = Stacks(type, f, st); i < n; ++i) Chimney(p, f, st[i].x, st[i].z, st[i].r, st[i].h, st[i].c);
    }
    if (type == PartType::HandPump || type == PartType::WaterPump) { // the intake pipe, cut at a block below or above (#221)
        const float lo = -std::min(0.6f, f.roomBelow), hi = top + std::min(0.3f, f.roomAbove);
        p.cyl({ x0 + 0.2f, lo, z0 + 0.2f }, type == PartType::HandPump ? 0.07f : 0.08f, type == PartType::HandPump ? 0.07f : 0.08f, hi - lo, kSteel);
    }
    if (f.roomAbove + 1e-3f < TopReach(type)) { // a block sits on it (#221): only the details on the casing's sides
        switch (type) {
        case PartType::BlastFurnace: Door(p, f, working, 0.35f, 0.3f); break;
        case PartType::Furnace: Door(p, f, working); break;
        case PartType::ArcFurnace: Door(p, f, working, 0.4f, 0.3f); break;
        case PartType::Generator: case PartType::MvGenerator: case PartType::HvGenerator: {
            const Color coil = type == PartType::HvGenerator ? Color{ 200, 60, 60, 255 } : type == PartType::MvGenerator ? Color{ 200, 120, 50, 255 } : kCopper;
            rlPushMatrix(); rlTranslatef(x0 - 0.01f, top * 0.5f, 0); rlRotatef(90, 0, 0, 1);
            p.cyl({ 0, -0.02f, 0 }, 0.28f, 0.28f, 0.18f, coil);
            rlPopMatrix();
            if (working) p.glow({ x0 - 0.12f, top * 0.5f, 0 }, { 0.02f, 0.3f, 0.3f }, { 120, 200, 255, 255 });
            break;
        }
        case PartType::Battery: p.cube({ 0, top * 0.5f, f.d * 0.5f + 0.01f }, { f.w * 0.7f, 0.08f, 0.02f }, kPaint); break;
        default: break; // the rest have only top details
        }
        return;
    }
    switch (type) {
    case PartType::Kiln:
        Door(p, f, working);
        break;
    case PartType::CokeOven:
        Door(p, f, working, 0.6f, 0.35f);
        break;
    case PartType::BlastFurnace:
        p.cyl({ 0, top + 0.5f, 0 }, 0.22f, 0.22f, 0.1f, kDarkSteel); // a charging ring
        Door(p, f, working, 0.35f, 0.3f);
        break;
    case PartType::Furnace:
        Door(p, f, working);
        for (int k = 0; k < 3; ++k) p.cube({ x0 + 0.25f + k * 0.25f, top + 0.02f, 0 }, { 0.12f, 0.04f, f.d * 0.6f }, kDarkSteel);
        break;
    case PartType::ArcFurnace:
        for (int k = 0; k < 3; ++k) {
            const float a = k * 2.094f;
            p.cyl({ 0.18f * cosf(a), top, 0.18f * sinf(a) }, 0.04f, 0.04f, 0.55f, { 60, 60, 66, 255 });
            if (working) p.glow({ 0.18f * cosf(a), top + 0.02f, 0.18f * sinf(a) }, { 0.1f, 0.04f, 0.1f }, { 200, 230, 255, 255 });
        }
        Door(p, f, working, 0.4f, 0.3f);
        break;
    case PartType::AlloySmelter:
        for (int k = 0; k < 2; ++k) {
            p.cyl({ (k ? 0.22f : -0.22f), top, 0 }, 0.17f, 0.12f, 0.2f, kBrick);
            if (working) p.glow({ (k ? 0.22f : -0.22f), top + 0.19f, 0 }, { 0.2f, 0.02f, 0.2f }, kGlowOrange);
        }
        break;
    case PartType::Boiler:
        rlPushMatrix(); rlTranslatef(0, top + 0.22f, 0); rlRotatef(90, 0, 0, 1);
        p.cyl({ 0, -f.w * 0.42f, 0 }, 0.2f, 0.2f, f.w * 0.84f, kCopper);    // a horizontal tank
        rlPopMatrix();
        p.cyl({ -x0 - 0.2f, top + 0.42f, 0 }, 0.06f, 0.06f, 0.04f, { 230, 230, 230, 255 }); // gauge
        break;
    case PartType::ForgeHammer: {
        p.cube({ x0 + 0.08f, top + 0.4f, 0 }, { 0.08f, 0.8f, 0.12f }, kDarkSteel);
        p.cube({ -x0 - 0.08f, top + 0.4f, 0 }, { 0.08f, 0.8f, 0.12f }, kDarkSteel);
        p.cube({ 0, top + 0.82f, 0 }, { f.w, 0.08f, 0.14f }, kDarkSteel);
        const float lift = working ? 0.25f * (0.5f + 0.5f * sinf(t * 5.0f)) : 0.25f;
        p.cube({ 0, top + 0.2f + lift, 0 }, { 0.35f, 0.18f, 0.3f }, kSteel);   // the hammer head
        p.cube({ 0, top + 0.03f, 0 }, { 0.45f, 0.06f, 0.4f }, kDarkSteel);      // anvil
        break;
    }
    case PartType::Macerator:
    case PartType::SteamMacerator:
        p.cyl({ 0, top, 0 }, 0.38f, 0.18f, 0.35f, kSteel);                       // hopper
        p.cyl({ 0, top + 0.35f, 0 }, 0.4f, 0.4f, 0.04f, kDarkSteel);
        rlPushMatrix(); rlTranslatef(0, top + 0.12f, f.d * 0.5f + 0.02f); rlRotatef(spin * 400.0f, 0, 0, 1);
        p.cube({ 0, 0, 0 }, { 0.3f, 0.05f, 0.04f }, kDarkSteel); p.cube({ 0, 0, 0 }, { 0.05f, 0.3f, 0.04f }, kDarkSteel); // grinder wheel
        rlPopMatrix();

        break;
    case PartType::Generator:
    case PartType::MvGenerator:
    case PartType::HvGenerator: {
        const Color coil = type == PartType::HvGenerator ? Color{ 200, 60, 60, 255 } : type == PartType::MvGenerator ? Color{ 200, 120, 50, 255 } : kCopper;
        rlPushMatrix(); rlTranslatef(x0 - 0.01f, top * 0.5f, 0); rlRotatef(90, 0, 0, 1);
        p.cyl({ 0, -0.02f, 0 }, 0.28f, 0.28f, 0.18f, coil);                     // the coil drum on the side
        rlPopMatrix();
        rlPushMatrix(); rlTranslatef(0, top + 0.05f, 0); rlRotatef(spin * 720.0f, 0, 1, 0);
        p.cube({ 0, 0, 0 }, { 0.6f, 0.03f, 0.1f }, kDarkSteel); p.cube({ 0, 0, 0 }, { 0.1f, 0.03f, 0.6f }, kDarkSteel); // fan
        rlPopMatrix();
        p.cyl({ 0, top, 0 }, 0.34f, 0.34f, 0.02f, kSteel);
        if (working) p.glow({ x0 - 0.12f, top * 0.5f, 0 }, { 0.02f, 0.3f, 0.3f }, { 120, 200, 255, 255 });
        break;
    }
    case PartType::Battery:
        p.cyl({ -0.2f, top, 0 }, 0.07f, 0.07f, 0.12f, { 200, 40, 40, 255 });
        p.cyl({ 0.2f, top, 0 }, 0.07f, 0.07f, 0.12f, { 30, 30, 30, 255 });
        p.cube({ 0, top * 0.5f, f.d * 0.5f + 0.01f }, { f.w * 0.7f, 0.08f, 0.02f }, kPaint);
        break;
    case PartType::OreWasher:
        rlPushMatrix(); rlTranslatef(0, top + 0.25f, 0); rlRotatef(90, 0, 0, 1); rlRotatef(spin * 180.0f, 0, 1, 0);
        p.cyl({ 0, -f.w * 0.4f, 0 }, 0.22f, 0.22f, f.w * 0.8f, kSteel);        // the washing drum, turning
        rlPopMatrix();
        p.cube({ 0, top + 0.03f, 0 }, { f.w * 0.9f, 0.06f, 0.5f }, kDarkSteel);
        break;
    case PartType::Centrifuge:
        rlPushMatrix(); rlTranslatef(0, top, 0); rlRotatef(spin * 900.0f, 0, 1, 0);
        p.cyl({ 0, 0, 0 }, 0.4f, 0.4f, 0.3f, kSteel);
        p.cube({ 0.3f, 0.15f, 0 }, { 0.14f, 0.12f, 0.14f }, kPaint);           // a marker to see it spin
        rlPopMatrix();
        break;
    case PartType::HandPump: // (#280) a wooden post with a lever and a spout, its intake pipe down into the water
        p.cube({ 0.0f, top + 0.35f, 0.0f }, { 0.16f, 0.7f, 0.16f }, { 110, 80, 55, 255 });
        p.cube({ 0.22f, top + 0.68f, 0.0f }, { 0.55f, 0.05f, 0.06f }, kSteel);
        p.cyl({ -0.1f, top + 0.45f, 0.0f }, 0.04f, 0.04f, 0.18f, kSteel);
        break;
    case PartType::WaterPump:
        p.cyl({ 0.1f, top, 0.1f }, 0.28f, 0.28f, 0.45f, { 70, 120, 190, 255 }); // tank
        break;
    case PartType::Electrolyzer:
        p.cyl({ -0.22f, top, 0 }, 0.16f, 0.16f, 0.5f, { 230, 230, 230, 255 }); // oxygen
        p.cyl({ 0.22f, top, 0 }, 0.16f, 0.16f, 0.5f, { 120, 160, 230, 255 });  // hydrogen
        if (working) p.glow({ 0, top + 0.05f, 0 }, { 0.1f, 0.1f, 0.1f }, { 150, 220, 255, 255 });
        break;
    case PartType::ChemicalReactor:
        p.cyl({ 0, top, 0 }, 0.3f, 0.36f, 0.45f, { 120, 170, 120, 255 });
        p.cyl({ 0, top + 0.45f, 0 }, 0.12f, 0.3f, 0.15f, { 120, 170, 120, 255 }, 0.9f);
        p.cyl({ x0 + 0.12f, top, 0 }, 0.05f, 0.05f, 0.6f, kSteel);
        break;
    case PartType::Refinery:
        p.cyl({ 0, top, 0 }, 0.22f, 0.26f, 1.6f, kSteel);                      // the distillation column
        for (int k = 1; k <= 4; ++k) p.cyl({ 0, top + 0.35f * k, 0 }, 0.29f, 0.29f, 0.04f, kDarkSteel);
        p.cyl({ 0, top + 1.6f, 0 }, 0.05f, 0.05f, 0.25f, kDarkSteel);
        break;
    case PartType::Miner: {
        // A derrick and a drill bit under the casing's front edge, turning while it works.
        for (int k = 0; k < 4; ++k) p.cube({ (k & 1 ? 0.35f : -0.35f), top + 0.45f, (k & 2 ? 0.35f : -0.35f) }, { 0.05f, 0.9f, 0.05f }, kPaint);
        p.cube({ 0, top + 0.9f, 0 }, { 0.8f, 0.06f, 0.8f }, kPaint);
        rlPushMatrix(); rlTranslatef(0, top + 0.1f, 0); rlRotatef(spin * 300.0f, 0, 1, 0);
        p.cyl({ 0, 0, 0 }, 0.12f, 0.02f, 0.7f, kDarkSteel);
        rlPopMatrix();
        break;
    }
    case PartType::PrecisionAssembler: {
        p.cube({ 0, top + 0.02f, 0 }, { f.w * 0.8f, 0.04f, f.d * 0.8f }, kDarkSteel);
        const float a = working ? 40.0f * sinf(t * 1.6f) : 0.0f;
        rlPushMatrix(); rlTranslatef(x0 + 0.2f, top, 0); rlRotatef(a, 0, 1, 0);
        p.cube({ 0, 0.25f, 0 }, { 0.08f, 0.5f, 0.08f }, kPaint);
        rlTranslatef(0, 0.5f, 0); rlRotatef(-60.0f + (working ? 20.0f * sinf(t * 2.3f) : 0.0f), 0, 0, 1);
        p.cube({ 0.22f, 0, 0 }, { 0.44f, 0.06f, 0.06f }, kPaint);               // the arm
        rlPopMatrix();
        rlPushMatrix(); rlRotatef(f.yawDeg, 0, 1, 0);
        p.glow({ 0, top * 0.6f, (fabsf(f.front.x) > 0.5f ? f.w : f.d) * 0.5f + 0.01f }, { f.w * 0.4f, f.h * 0.25f, 0.02f }, working ? Color{ 80, 230, 140, 255 } : Color{ 20, 50, 40, 255 }); // screen
        rlPopMatrix();
        break;
    }
    case PartType::LvLab: case PartType::MvLab: case PartType::HvLab: case PartType::SteamLab: { // (#1193; #1271 L6: the steam lab too) a research lab: a white bench, flasks, a dish on a mast
        const Color band = type == PartType::SteamLab ? Color{ 190, 150, 90, 255 } : type == PartType::LvLab ? Color{ 90, 150, 220, 255 } : type == PartType::MvLab ? Color{ 230, 150, 60, 255 } : Color{ 180, 90, 210, 255 }; // its tier's colour
        p.cube({ 0, top * 0.55f, 0 }, { f.w * 1.01f, 0.06f, f.d * 1.01f }, band);
        for (int k = -1; k <= 1; ++k) p.cyl({ k * 0.18f, top, 0.12f }, 0.05f, 0.03f, 0.16f, kGlass);  // flasks
        rlPushMatrix(); rlTranslatef(-0.25f, top, -0.2f); rlRotatef(working ? t * 45.0f : 0.0f, 0, 1, 0);
        p.cyl({ 0, 0, 0 }, 0.025f, 0.025f, 0.3f, kDarkSteel);                                  // the mast
        p.cyl({ 0, 0.3f, 0 }, 0.02f, 0.14f, 0.1f, kSteel);                                      // the dish
        rlPopMatrix();
        rlPushMatrix(); rlRotatef(f.yawDeg, 0, 1, 0);
        p.glow({ 0, top * 0.3f, (fabsf(f.front.x) > 0.5f ? f.w : f.d) * 0.5f + 0.01f }, { f.w * 0.5f, f.h * 0.18f, 0.02f }, working ? Color{ 120, 200, 255, 255 } : Color{ 20, 40, 60, 255 }); // its screen
        rlPopMatrix();
        break;
    }
    case PartType::ResearchTable: { // (#1252, #1327) a wooden table: an open tome and a crystal on it
        p.cube({ 0.1f, top + 0.02f, 0.0f }, { 0.36f, 0.04f, 0.26f }, Color{ 225, 215, 190, 255 });                 // the open tome
        p.cube({ 0.1f, top + 0.045f, 0.0f }, { 0.02f, 0.01f, 0.26f }, Color{ 120, 70, 50, 255 });                 // its spine
        p.glow({ -0.28f, top + 0.08f, 0.15f }, { 0.07f, 0.14f, 0.07f }, Color{ 170, 230, 255, 255 });             // a resonant crystal
        break;
    }
    case PartType::ManaSiphon: { // (#1252) a stone base and a violet crystal column that pulses while it draws
        const float pulse = working ? 0.6f + 0.4f * sinf(t * 2.5f) : 0.2f;
        p.cyl({ 0, top, 0 }, 0.12f, 0.06f, 0.55f, kDarkSteel);
        p.glow({ 0, top + 0.62f, 0 }, { 0.1f, 0.16f, 0.1f }, Color{ (unsigned char)(170 * pulse), (unsigned char)(130 * pulse), (unsigned char)(255 * pulse), 255 });
        break;
    }
    case PartType::ManaDynamo: { // (#1252) a coil that spins while it turns mana into power
        rlPushMatrix(); rlTranslatef(0, top + 0.12f, 0); rlRotatef(working ? t * 120.0f : 0.0f, 0, 1, 0);
        for (int k = 0; k < 4; ++k) p.cube({ 0.16f * cosf(k * 1.5708f), 0, 0.16f * sinf(k * 1.5708f) }, { 0.08f, 0.2f, 0.08f }, kCopper);
        rlPopMatrix();
        p.glow({ 0, top + 0.26f, 0 }, { 0.08f, 0.04f, 0.08f }, working ? Color{ 160, 140, 255, 255 } : Color{ 40, 30, 70, 255 });
        break;
    }
    case PartType::FusionReactor: {
        const float pulse = working ? 0.6f + 0.4f * sinf(t * 3.0f) : 0.25f;
        for (int k = 0; k < 16; ++k) {
            const float a = k * 0.3927f;
            const Vector3 c = { 0.42f * cosf(a), top + 0.18f, 0.42f * sinf(a) };
            p.cube(c, { 0.16f, 0.22f, 0.16f }, kSteel);
            p.glow({ c.x * 0.8f, c.y, c.z * 0.8f }, { 0.06f, 0.06f, 0.06f }, Color{ (unsigned char)(200 * pulse), (unsigned char)(120 * pulse), (unsigned char)(255 * pulse), 255 });
        }
        break;
    }
    case PartType::MassDriver:
        rlPushMatrix(); rlTranslatef(0, top + 0.1f, 0); rlRotatef(f.yawDeg, 0, 1, 0); rlRotatef(-35.0f, 1, 0, 0);
        p.cube({ 0, 0, 0.8f }, { 0.22f, 0.22f, 1.8f }, kSteel);                // the launch rail, raised
        for (int k = 0; k < 5; ++k) p.cube({ 0, 0, 0.1f + 0.35f * k }, { 0.3f, 0.3f, 0.06f }, working ? kGlass : kCopper);
        rlPopMatrix();
        break;
    case PartType::Snare: { // (#176) a bent sapling and a noose on the frame; with a catch, a small furry animal in it
        const float slab = top * 0.5f; // the frame is a half-cell slab
        const Color wood = { 120, 90, 55, 255 }, cord = { 200, 180, 130, 255 };
        p.cube({ x0 + 0.06f, slab + 0.18f, 0 }, { 0.04f, 0.36f, 0.04f }, wood);               // the spring pole
        p.cube({ x0 + 0.2f, slab + 0.36f, 0 }, { 0.3f, 0.03f, 0.03f }, wood);                 // bent over
        p.cube({ 0.05f, slab + 0.2f, 0 }, { 0.015f, 0.32f, 0.015f }, cord);                   // the line
        p.cube({ 0.05f, slab + 0.03f, 0 }, { 0.14f, 0.02f, 0.14f }, cord);                    // the noose
        if (working) {
            p.cube({ 0.05f, slab + 0.08f, 0.02f }, { 0.16f, 0.1f, 0.24f }, { 130, 105, 80, 255 });   // a hare-sized body
            p.cube({ 0.05f, slab + 0.14f, 0.15f }, { 0.08f, 0.08f, 0.08f }, { 140, 112, 86, 255 });  // head
            p.cube({ 0.03f, slab + 0.21f, 0.15f }, { 0.02f, 0.08f, 0.02f }, { 150, 120, 95, 255 });  // ears
            p.cube({ 0.07f, slab + 0.21f, 0.15f }, { 0.02f, 0.08f, 0.02f }, { 150, 120, 95, 255 });
        }
        break;
    }
    case PartType::Torch: { // (#218) the wrapped head on the thin stick (its low x/z corner: placed torches aren't turned)
        // #558: the flame is particles (ambience.cpp, from Machine_SmokePoints); only a small glowing core stays here, so
        // the head reads lit where no particles are drawn.
        p.cube({ x0 + kTorchStick, top - 0.03f, z0 + kTorchStick }, { 0.13f, 0.08f, 0.13f }, kDarkSteel);
        p.glow({ x0 + kTorchStick, top + 0.03f, z0 + kTorchStick }, { 0.07f, 0.05f, 0.07f }, Color{ 255, 190, 80, 255 });
        break;
    }
    case PartType::Campfire: { // #552: a firepit - a ring of stones round the ash bed, logs leaning in a teepee (flames: ambience.cpp)
        const Color stone = { 150, 145, 138, 255 }, stoneDark = { 122, 117, 110, 255 }, bark = { 92, 64, 40, 255 }, charred = { 40, 32, 28, 255 };
        const float bed = 0.125f;
        for (int k = 0; k < 10; ++k) {
            const float a = 2.0f * PI * k / 10.0f + 0.13f * (k % 3);
            const float s = 0.13f + 0.03f * ((k * 7) % 3);                                        // a little irregular
            p.cube({ 0.4f * cosf(a), bed + s * 0.35f, 0.4f * sinf(a) }, { s, s * 0.75f, s * 0.9f }, (k & 1) ? stone : stoneDark);
        }
        for (int k = 0; k < 5; ++k) { // tilted 32 degrees in, feet on a 0.24 m circle, tops meeting over the middle
            const float a = 360.0f * k / 5.0f + 15.0f;
            rlPushMatrix();
            rlRotatef(a, 0, 1, 0);
            rlTranslatef(0.24f, bed, 0.0f);
            rlRotatef(32.0f, 0, 0, 1); // leans toward -x: the centre
            p.cyl({ 0, 0, 0 }, 0.035f, 0.045f, 0.5f, k == 0 ? charred : bark);
            p.cube({ 0, 0.02f, 0 }, { 0.1f, 0.06f, 0.1f }, charred);                              // the burnt foot
            rlPopMatrix();
        }
        if (working) p.glow({ 0, bed + 0.04f, 0 }, { 0.22f, 0.06f, 0.22f }, Color{ 255, 120, 30, 255 }); // embers
        break;
    }
    case PartType::Chest:
        p.cube({ 0, top - 0.06f, 0 }, { f.w * 1.02f, 0.06f, f.d * 1.02f }, kDarkSteel);  // a banded lid
        rlPushMatrix(); rlRotatef(f.yawDeg, 0, 1, 0);
        p.cube({ 0, top * 0.7f, (fabsf(f.front.x) > 0.5f ? f.w : f.d) * 0.5f + 0.02f }, { 0.12f, 0.14f, 0.04f }, kPaint); // the latch
        rlPopMatrix();
        break;
    default: break;
    }
}

} // namespace

// Per plot, rebuilt when its version changes: each machine's floor centre (plot metres), size and facing.
namespace {
    struct MachineSpot { int part; PartType type; Vector3 floorPlot; Frame f; float buried = 0.0f; }; // buried: 1 - sky exposure (#284)
    struct PlotCache { uint32_t version = 0xFFFFFFFFu; std::vector<MachineSpot> spots; };
    std::unordered_map<uint64_t, PlotCache> plotCache; // key: body id ^ plot id

    // A machine's frame from its footprint: the floor centre in plot metres, size and facing.
    bool SpotOf(const Plot& plot, size_t i, MachineSpot& s) {
        static thread_local std::vector<std::pair<Int3, CellShape>>* fpP = nullptr; if (!fpP) fpP = new std::vector<std::pair<Int3, CellShape>>(); auto& fp = *fpP; // (#1133, #485: MinGW frees thread_local storage before running its destructor at thread exit: never destroyed)
        const Part& part = plot.parts[i];
        Plot::Footprint(part, fp);
        if (fp.empty()) return false;
        Int3 lo = fp[0].first, hi = fp[0].first;
        for (auto& c : fp) { lo = { std::min(lo.x, c.first.x), std::min(lo.y, c.first.y), std::min(lo.z, c.first.z) }; hi = { std::max(hi.x, c.first.x), std::max(hi.y, c.first.y), std::max(hi.z, c.first.z) }; }
        s.part = (int)i; s.type = part.type;
        s.floorPlot = { (lo.x + hi.x + 1) * 0.5f * kCellSize, lo.y * kCellSize, (lo.z + hi.z + 1) * 0.5f * kCellSize };
        s.f.w = (hi.x - lo.x + 1) * kCellSize; s.f.h = (hi.y - lo.y + 1) * kCellSize; s.f.d = (hi.z - lo.z + 1) * kCellSize;
        const int* m = Orientation(part.orient);
        Vector3 front = { (float)m[2], 0.0f, (float)m[8] };              // the part's +z, flattened
        if (Vector3LengthSqr(front) < 0.5f) front = { 0, 0, 1 };
        s.f.front = Vector3Normalize(front);
        s.f.yawDeg = atan2f(s.f.front.x, s.f.front.z) * RAD2DEG;
        // Room for details (#221): whole free cell layers over the footprint (4 = open), and under its intake
        // corner (the pump pipes stand at the footprint's low x/z corner).
        auto other = [&](int x, int y, int z) { const int q = plot.partAt({ x, y, z }); return q >= 0 && q != (int)i; };
        for (int k = 1; k <= 4; ++k) {
            bool hit = false;
            for (int x = lo.x; x <= hi.x && !hit; ++x)
                for (int z = lo.z; z <= hi.z && !hit; ++z) hit = other(x, hi.y + k, z);
            if (hit) { s.f.roomAbove = (k - 1) * kCellSize; break; }
        }
        for (int k = 1; k <= 2; ++k)
            if (other(lo.x, lo.y - k, lo.z)) { s.f.roomBelow = (k - 1) * kCellSize; break; }
        return true;
    }
}

int Machine_SmokePoints(const Plot& plot, size_t part, Vector3* outBody, int max) {
    MachineSpot s;
    if (max <= 0 || part >= plot.parts.size() || !SpotOf(plot, part, s)) return 0;
    Stack st[4];
    const int n = Stacks(s.type, s.f, st);
    int k = 0;
    for (int i = 0; i < n && k < max; ++i) // the chimney tops (#220: the smoke came off the wrong corner, in the air)
        outBody[k++] = plot.toBody({ s.floorPlot.x + st[i].x, s.floorPlot.y + s.f.h + st[i].h + 0.06f, s.floorPlot.z + st[i].z });
    if (k == 0 && n == 0 && StacksAll(s.type, s.f, st) > 0) return 0; // chimneys capped by a block above: no smoke (#221)
    if (s.type == PartType::Torch) { // the flame's foot, just over the head on the stick (#558)
        outBody[0] = plot.toBody({ s.floorPlot.x - s.f.w * 0.5f + kTorchStick, s.floorPlot.y + s.f.h + 0.05f, s.floorPlot.z - s.f.d * 0.5f + kTorchStick });
        return 1;
    }
    if (k == 0) { // no chimney: off the middle of the top (the campfire: just over its flames)
        const float up = s.type == PartType::Campfire ? 0.5f : s.f.h + 0.05f;
        outBody[k++] = plot.toBody({ s.floorPlot.x, s.floorPlot.y + up, s.floorPlot.z });
    }
    return k;
}

static int drawnLastFrame = 0;
int Machines_DrawnCount() { return drawnLastFrame; }

void Machines_Draw3D(Scene& world) {
    rlDrawRenderBatchActive();
    rlDisableBackfaceCulling(); // the plot frame (east, up, north) is left-handed: mirrored winding
    drawnLastFrame = 0;
    const Vector3 sun = Vector3Normalize(world.sunDirection);
    std::vector<std::pair<Int3, CellShape>> fp;
    struct Pending { float d2; Planetoid* body; Plot* plot; const MachineSpot* spot; Vector3 e, u, n; float light; };
    static std::vector<Pending> pending;
    pending.clear();
    for (const auto& obj : world.objects) {
        Planetoid* body = dynamic_cast<Planetoid*>(obj.get());
        if (!body || !body->isActive || !body->plots) continue;
        const Vector3 camPlotBody = body->toLocal(world.camera.position);
        for (auto& plot : body->plots->plots) {
            if (Vector3Distance(plot->toBody({ 0, 0, 0 }), camPlotBody) > 500.0f) continue;
            if (!gOcclusionSkipPlots.empty() && gOcclusionSkipPlots.count(plot->id)) continue; // #134 experiment
            PlotCache& pc = plotCache[(body->id * 1315423911ull) ^ plot->id];
            if (pc.version != plot->version) {
                pc.version = plot->version;
                pc.spots.clear();
                for (size_t i = 0; i < plot->parts.size(); ++i) {
                    const Part& part = plot->parts[i];
                    if (part.type == PartType::Count || (!IsMachine(part.type) && part.type != PartType::HandPump && part.type != PartType::Torch && part.type != PartType::Snare)) continue;
                    MachineSpot s;
                    if (!SpotOf(*plot, i, s)) continue;
                    if (body->field) { // sky exposure 1 m over the floor: the terrain's rule, as plot cells (#216, #284)
                        const Vector3 at = plot->toBody({ s.floorPlot.x, s.floorPlot.y + 1.0f, s.floorPlot.z });
                        const float len = Vector3Length(at);
                        float sky = 1.0f;
                        if (len > 1.0f) body->field->materialSky(at, 1.0f, 1.0f, Vector3Scale(at, 1.0f / len), sky);
                        s.buried = 1.0f - sky;
                    }
                    pc.spots.push_back(s);
                }
            }
            if (pc.spots.empty()) continue;
            const Vector3 camPlot = plot->toPlot(camPlotBody);                     // camera in plot metres: cheap culling
            const Vector3 e = body->dirToWorld(plot->east), u = body->dirToWorld(plot->up), n = body->dirToWorld(plot->north);
            Painter p;
            p.light = 0.45f + 0.55f * std::clamp(Vector3DotProduct(sun, u) * 2.0f + 0.5f, 0.0f, 1.0f);
            for (const MachineSpot& s : pc.spots) {
                const float d2 = Vector3DistanceSqr(s.floorPlot, camPlot);
                if (d2 > 32.0f * 32.0f) continue; // the casings read on their own beyond this
                Pending q; q.d2 = d2; q.body = body; q.plot = plot.get(); q.spot = &s; q.e = e; q.u = u; q.n = n; q.light = p.light;
                pending.push_back(q);
            }
        }
    }
    // The nearest kMaxDetailed only: a huge factory in view costs a bounded amount.
    constexpr size_t kMaxDetailed = 120;
    if (pending.size() > kMaxDetailed) {
        std::nth_element(pending.begin(), pending.begin() + kMaxDetailed, pending.end(), [](const Pending& a, const Pending& b) { return a.d2 < b.d2; });
        pending.resize(kMaxDetailed);
    }
    for (const Pending& q : pending) {
        const MachineSpot& s = *q.spot;
        const bool working = s.type == PartType::Snare ? Trapping_HasCatch(q.body->id, q.plot->id, s.part) : Factory::Working(q.body->id, q.plot->id, s.part); // a snare: "working" = holding a catch (#176)
        // The plot's own model, its translation in double (#284): placed through the float toWorld ~60 km from the
        // planet's centre (~4 mm steps, float trig every frame) the details shook against the casing, which is drawn
        // with the double localModel. Plot axes in body space, then the body's rotation and the anchor's offset.
        const Plot& pl = *q.plot;
        const Vector3 off = Vector3Add(Vector3Add(Vector3Scale(pl.east, s.floorPlot.x), Vector3Scale(pl.up, s.floorPlot.y)), Vector3Scale(pl.north, s.floorPlot.z));
        Matrix axes = MatrixIdentity();
        axes.m0 = pl.east.x; axes.m4 = pl.up.x; axes.m8 = pl.north.x;
        axes.m1 = pl.east.y; axes.m5 = pl.up.y; axes.m9 = pl.north.y;
        axes.m2 = pl.east.z; axes.m6 = pl.up.z; axes.m10 = pl.north.z;
        const Matrix mm = MatrixMultiply(axes, q.body->localModel(DVec3(pl.origin.x + off.x, pl.origin.y + off.y, pl.origin.z + off.z), world.renderOrigin));
        Painter pp; pp.light = q.light;
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(mm));
        IB_SetBuried(s.buried);
        DrawMachine(s.type, s.f, pp, working);
        IB_SetBuried(0.0f);
        ++drawnLastFrame;
        rlPopMatrix();
    }
    rlDrawRenderBatchActive();
    rlEnableBackfaceCulling();
    RotationModels_Draw3D(world); // (#1269) shafts, cogs, cranks and the rotation machines' turning parts
}

bool Machine_RaycastDetails(const Plot& plot, Vector3 o, Vector3 d, float maxDist, int& part, float& dist) {
    // Chimneys stand out of the casing's cells (#285): a box round each, in plot metres (slab test).
    bool hit = false;
    for (size_t i = 0; i < plot.parts.size(); ++i) {
        if (plot.parts[i].type == PartType::Count || (!IsMachine(plot.parts[i].type) && plot.parts[i].type != PartType::HandPump)) continue;
        MachineSpot s;
        Stack st[4];
        if (!SpotOf(plot, i, s)) continue;
        const int n = Stacks(s.type, s.f, st);
        for (int k = 0; k < n; ++k) {
            const float r = st[k].r * 1.25f;
            const Vector3 lo = { s.floorPlot.x + st[k].x - r, s.floorPlot.y + s.f.h, s.floorPlot.z + st[k].z - r };
            const Vector3 hi = { s.floorPlot.x + st[k].x + r, s.floorPlot.y + s.f.h + st[k].h + 0.06f, s.floorPlot.z + st[k].z + r };
            float t0 = 0.0f, t1 = maxDist;
            const float oa[3] = { o.x, o.y, o.z }, da[3] = { d.x, d.y, d.z }, la[3] = { lo.x, lo.y, lo.z }, ha[3] = { hi.x, hi.y, hi.z };
            bool miss = false;
            for (int a = 0; a < 3 && !miss; ++a) {
                if (fabsf(da[a]) < 1e-8f) { miss = oa[a] < la[a] || oa[a] > ha[a]; continue; }
                float ta = (la[a] - oa[a]) / da[a], tb = (ha[a] - oa[a]) / da[a];
                if (ta > tb) std::swap(ta, tb);
                t0 = std::max(t0, ta); t1 = std::min(t1, tb);
                miss = t0 > t1;
            }
            if (!miss && t0 < maxDist) { maxDist = t0; part = (int)i; dist = t0; hit = true; }
        }
    }
    return hit;
}

bool Machines_TestStack() {
    // --test-machine-stack (#221): a block over a kiln caps its chimney (no model, smoke or hit box), freed again on removal;
    // a pump on a block keeps its pipe out of it.
    Plot plot;
    plot.east = { 1, 0, 0 }; plot.up = { 0, 1, 0 }; plot.north = { 0, 0, 1 };
    const int kiln = plot.place({ PartType::Kiln, { 0, 0, 0 }, 0 });
    Vector3 pts[4];
    int part = -1; float dist = 0.0f;
    auto chimneyHit = [&] { MachineSpot s; Stack st[4]; SpotOf(plot, (size_t)kiln, s); StacksAll(s.type, s.f, st);
        const Vector3 o = { s.floorPlot.x + st[0].x, s.floorPlot.y + s.f.h + 0.3f, s.floorPlot.z + st[0].z - 2.0f };
        return Machine_RaycastDetails(plot, o, { 0, 0, 1 }, 4.0f, part, dist); };
    const int free0 = Machine_SmokePoints(plot, (size_t)kiln, pts, 4);
    const bool hit0 = chimneyHit();
    const int cap = plot.place({ PartType::Battery, { 0, 3, 0 }, 0 }); // one free layer, then a block: 0.5 m < the 0.66 m chimney
    const int capped = Machine_SmokePoints(plot, (size_t)kiln, pts, 4);
    const bool hit1 = chimneyHit();
    plot.remove(cap);
    const int freed = Machine_SmokePoints(plot, (size_t)kiln, pts, 4);
    MachineSpot pump;
    const int pumpIdx = plot.place({ PartType::WaterPump, { 0, 2, 0 }, 0 });
    SpotOf(plot, (size_t)pumpIdx, pump);
    MachineSpot k2; SpotOf(plot, (size_t)kiln, k2);
    Stack kst[4];
    const bool ok = kiln >= 0 && cap >= 0 && pumpIdx >= 0 && free0 == 1 && hit0 && capped == 0 && !hit1 && freed == 1
        && pump.f.roomBelow == 0.0f && k2.f.roomAbove == 0.0f && Stacks(k2.type, k2.f, kst) == 0;
    printf("bench: machine-stack: smoke free %d capped %d freed %d, chimney hit %d/%d, pump roomBelow %.2f, kiln roomAbove %.2f: %s\n",
        free0, capped, freed, (int)hit0, (int)hit1, pump.f.roomBelow, k2.f.roomAbove, ok ? "PASS" : "FAIL");
    return ok;
}
