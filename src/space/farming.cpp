#include "farming.h"
#include "audio.h"
#include "buildtool.h"
#include "gameui.h"
#include "inventory.h"
#include "plots.h"
#include "survival.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>

namespace {
    struct CropDef { std::string id, seed, yield; int yieldCount = 1; float growSeconds = 3600.0f; Color color = GREEN; };
    std::vector<CropDef> crops;
    bool cropsLoaded = false;

    void LoadCrops() {
        if (cropsLoaded) return;
        cropsLoaded = true;
        FILE* f = fopen("resources/data/crops.txt", "r");
        if (!f) return;
        char line[512];
        while (fgets(line, sizeof line, f)) {
            if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
            char id[64] = {}, seed[64] = {}, yield[64] = {};
            int n = 1, r = 0, g = 0, b = 0;
            float secs = 3600.0f;
            // id | seed item | yield item count | daylight seconds to ripen | colour r,g,b
            if (sscanf(line, " %63[^ |] | %63[^ |] | %63s %d | %f | %d , %d , %d", id, seed, yield, &n, &secs, &r, &g, &b) == 8)
                crops.push_back({ id, seed, yield, n, secs, { (unsigned char)r, (unsigned char)g, (unsigned char)b, 255 } });
        }
        fclose(f);
    }
    const CropDef* Crop(const std::string& id) { LoadCrops(); for (const CropDef& c : crops) if (c.id == id) return &c; return nullptr; }

    // Crop state per farmland part: "<body name>|<plot id>|<part>" -> crop and growth (0..1).
    struct Growing { std::string crop; float growth = 0.0f; };
    std::map<std::string, Growing> fields;

    std::string Key(const Planetoid& body, const Plot& plot, int part) { return body.name + "|" + std::to_string(plot.id) + "|" + std::to_string(part); }
    bool Parse(const std::string& key, std::string& body, uint32_t& plot, int& part) {
        const size_t a = key.find('|'), b = key.find('|', a + 1);
        if (a == std::string::npos || b == std::string::npos) return false;
        body = key.substr(0, a);
        plot = (uint32_t)std::stoul(key.substr(a + 1, b - a - 1));
        part = std::stoi(key.substr(b + 1));
        return true;
    }
    Planetoid* BodyNamed(Scene& world, const std::string& name) {
        for (auto& o : world.objects) if (auto* p = dynamic_cast<Planetoid*>(o.get()); p && p->name == name) return p;
        return nullptr;
    }
    bool IsFarmland(const Plot& plot, int part) { return part >= 0 && part < (int)plot.parts.size() && plot.parts[part].type == PartType::Farmland; }

    // The farmland's centre on top of the soil, body-local.
    Vector3 TopCentre(const Plot& plot, const Part& p) {
        return plot.toBody({ (p.origin.x + 1.0f) * kCellSize, (p.origin.y + 1.0f) * kCellSize, (p.origin.z + 1.0f) * kCellSize });
    }
}

void Farm_Reset() { fields.clear(); }

void Farm_Update(Scene& world, float dt) {
    if (dt <= 0.0f || fields.empty()) return;
    LoadCrops();
    const Vector3 sun = Vector3Normalize(world.sunDirection);
    for (auto it = fields.begin(); it != fields.end();) {
        std::string bodyName; uint32_t plotId = 0; int part = -1;
        Planetoid* body = Parse(it->first, bodyName, plotId, part) ? BodyNamed(world, bodyName) : nullptr;
        Plot* plot = body && body->plots ? body->plots->byId(plotId) : nullptr;
        if (body && plot && !IsFarmland(*plot, part)) { it = fields.erase(it); continue; } // dug up
        Growing& g = it->second;
        const CropDef* c = Crop(g.crop);
        if (!body || !plot || !c || g.growth >= 1.0f) { ++it; continue; }
        // Daylight only, twice as fast with water within 4 m (a river, a pond, the sea shore).
        const Vector3 top = TopCentre(*plot, plot->parts[part]);
        const Vector3 up = Vector3Normalize(top);
        if (Vector3DotProduct(sun, body->dirToWorld(up)) < 0.05f) { ++it; continue; }
        float rate = 1.0f;
        if (const auto* pf = dynamic_cast<const PlanetField*>(body->field.get())) {
            const float e = pf->elevation(up);
            if (e < 3.0f) rate = 2.0f; // the sea shore
            else if (pf->inlandWaterNear(up, 4.0f) && pf->waterElevation(up) != PlanetField::kNoWater) rate = 2.0f;
        }
        g.growth = std::min(1.0f, g.growth + dt * rate / c->growSeconds);
        ++it;
    }
}

bool Farm_Plant(Planetoid& body, Plot& plot, int part, const std::string& seed) {
    LoadCrops();
    if (!IsFarmland(plot, part)) return false;
    const std::string key = Key(body, plot, part);
    if (fields.count(key)) { GameUI_Toast("Something is growing here already"); return false; }
    for (const CropDef& c : crops) {
        if (!seed.empty() && c.seed != seed) continue;
        if (inventory.count(c.seed) <= 0 || !inventory.remove(c.seed, 1)) continue;
        fields[key] = { c.id, 0.0f };
        QuestBook_Count("plant");
        Audio_Play(Sfx::Pickup, 0.4f, 0.7f);
        return true;
    }
    GameUI_Toast("No seeds (harvest bushes for them)");
    Audio_Play(Sfx::Error, 0.3f);
    return false;
}

int Farm_Harvest(Planetoid& body, Plot& plot, int part) {
    auto it = fields.find(Key(body, plot, part));
    if (it == fields.end()) return 0;
    const CropDef* c = Crop(it->second.crop);
    if (!c) { fields.erase(it); return 0; }
    if (it->second.growth < 1.0f) { GameUI_Toast(TextFormat("The %s is not ripe yet (%.0f%%)", c->id.c_str(), it->second.growth * 100.0f)); return 0; }
    const int seeds = 1 + (GetRandomValue(0, 1));
    inventory.add(c->yield, c->yieldCount);
    inventory.add(c->seed, seeds);
    GameUI_Toast(TextFormat("+%d %s  +%d seeds", c->yieldCount, c->yield.c_str(), seeds));
    QuestBook_Count("harvest");
    Audio_Play(Sfx::Pickup, 0.5f);
    fields.erase(it);
    return c->yieldCount + seeds;
}

bool Farm_Crop(const Planetoid& body, const Plot& plot, int part, std::string& crop, float& growth) {
    auto it = fields.find(Key(body, plot, part));
    if (it == fields.end()) { crop.clear(); growth = 0.0f; return false; }
    crop = it->second.crop;
    growth = it->second.growth;
    return true;
}

void Farm_DrawAimHUD(const Planetoid& body, const Plot& plot, int part) {
    std::string crop; float growth = 0.0f;
    const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40;
    DrawRectangle(x - 8, y - 6, 360, 62, Fade(BLACK, 0.55f));
    if (Farm_Crop(body, plot, part, crop, growth)) {
        DrawText(TextFormat("farmland: %s %s", crop.c_str(), growth >= 1.0f ? "- ripe" : ""), x, y, 20, Color{ 140, 220, 120, 255 });
        DrawRectangle(x, y + 26, 300, 8, Fade(RAYWHITE, 0.2f));
        DrawRectangle(x, y + 26, (int)(300 * growth), 8, growth >= 1.0f ? GOLD : GREEN);
        DrawText(growth >= 1.0f ? "O harvest" : "grows in daylight; water nearby doubles it", x, y + 38, 14, Fade(RAYWHITE, 0.7f));
    } else {
        DrawText("farmland: empty", x, y, 20, Color{ 200, 170, 120, 255 });
        DrawText("I plant seeds you carry", x, y + 26, 16, Fade(RAYWHITE, 0.7f));
    }
}

bool Farm_Till(Scene& world, Planetoid& body, Vector3 ground) {
    (void)world;
    if (!body.plots) body.plots = std::make_unique<PlotSet>();
    Plot* plot = body.plots->plotFor(ground, 6.0f, false);
    Part part{ PartType::Farmland, { 0, 0, 0 }, 0, 0 };
    if (plot) {
        part.origin = plot->cellAt(plot->toPlot(Vector3Add(ground, Vector3Scale(plot->up, 0.1f))));
        part.origin.x &= ~1; part.origin.z &= ~1; // on a 1 m grid, so tilled squares line up
    } else {
        plot = body.plots->plotFor(ground, 0.0f, true);
    }
    if (!plot->canPlace(part)) return false;
    if (!BuildTool_Place(body, *plot, part, false)) return false;
    QuestBook_Count("till");
    return true;
}

void Farm_Draw3D(Scene& world) {
    if (fields.empty()) return;
    LoadCrops();
    for (const auto& kv : fields) {
        std::string bodyName; uint32_t plotId = 0; int part = -1;
        Planetoid* body = Parse(kv.first, bodyName, plotId, part) ? BodyNamed(world, bodyName) : nullptr;
        Plot* plot = body && body->plots ? body->plots->byId(plotId) : nullptr;
        if (!plot || !IsFarmland(*plot, part)) continue;
        const CropDef* c = Crop(kv.second.crop);
        if (!c) continue;
        const Part& p = plot->parts[part];
        const Vector3 top = body->toWorld(TopCentre(*plot, p));
        if (Vector3Distance(top, world.camera.position) > 120.0f) continue;
        // Plot axes in render space: the crop stands along the plot's up.
        const Vector3 e = body->dirToWorld(plot->east), u = body->dirToWorld(plot->up), n = body->dirToWorld(plot->north);
        const float g = kv.second.growth;
        const int stage = std::min(3, (int)(g * 4.0f));
        // Four stages: sprouts, leafy, tall, ripe (heads in the crop's colour on top).
        static const float kH[4] = { 0.08f, 0.25f, 0.5f, 0.7f };
        const bool ripe = g >= 1.0f;
        const float h = kH[ripe ? 3 : std::min(stage, 2)];
        const Color col = ColorLerp({ 70, 150, 50, 255 }, { 110, 150, 60, 255 }, g);
        // Right-handed (east, up, south): the plot frame (east, up, north) is mirrored, which would
        // turn the cubes inside out under backface culling.
        Matrix m = { e.x, u.x, -n.x, top.x, e.y, u.y, -n.y, top.y, e.z, u.z, -n.z, top.z, 0, 0, 0, 1 };
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(m));
        const float w = ripe ? 0.07f : 0.05f + 0.025f * std::min(stage, 2);
        for (int k = 0; k < 9; ++k) {
            const float x = (k % 3 - 1) * 0.3f + ((k * 37) % 7 - 3) * 0.012f, z = (k / 3 - 1) * 0.3f + ((k * 53) % 7 - 3) * 0.012f;
            const float hk = h * (0.85f + 0.05f * (k % 4));
            DrawCube({ x, hk * 0.5f, z }, w, hk, w, col);
            if (ripe) DrawCube({ x, hk + 0.06f, z }, w * 1.8f, 0.12f, w * 1.8f, c->color); // the head
        }
        rlPopMatrix();
    }
}

void Farm_Save(FILE* f) {
    for (const auto& kv : fields) fprintf(f, "crop %s %s %.4f\n", kv.first.c_str(), kv.second.crop.c_str(), kv.second.growth);
}

bool Farm_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "crop") != 0) return false;
    char key[256] = {}, crop[64] = {};
    float g = 0.0f;
    if (fscanf(f, " %255s %63s %f", key, crop, &g) == 3) fields[key] = { crop, std::clamp(g, 0.0f, 1.0f) };
    return true;
}
