#include "farming.h"
#include "gamerng.h"
#include "terraintool.h"
#include "fluid.h"
#include "audio.h"
#include "buildtool.h"
#include "gameui.h"
#include "inventory.h"
#include "plots.h"
#include "survival.h"
#include "weather.h"
#include "player.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include "net.h"
#include "netproto.h"

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
    struct Growing { std::string crop; float growth = 0.0f; float rate = 0.0f; double wetChecked = -1e9; bool flowing = false; }; // rate: 1, or 2 by static water (worked out once; 0 = not yet); flowing: moving water within 4 m (checked every 2 s)
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
    if (gNetRole == NetRole::Client) return; // (#335) the server grows them; its fields come as FarmFields
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
        if (g.rate <= 0.0f) { // the terrain queries are costly: once per field
            g.rate = 1.0f;
            if (const auto* pf = dynamic_cast<const PlanetField*>(body->field.get())) {
                const float e = pf->elevation(up);
                if (e < 3.0f) g.rate = 2.0f; // the sea shore
                else if (pf->inlandWaterNear(up, 4.0f)) {
                    const float w = pf->waterElevation(up);
                    // A pond the simulation drained is dry (#360): the generated level no longer holds there.
                    if (w != PlanetField::kNoWater && Fluid::GeneratedWaterValid(body->id, Vector3Scale(up, pf->radius + w - 0.2f))) g.rate = 2.0f;
                }
            }
        }
        // Rain on the field speeds it by half again (the weather where the player is: fields
        // more than 3 km away grow at their own rate).
        const bool rainy = (Weather_Here().kind == WeatherKind::Rain || Weather_Here().kind == WeatherKind::Storm) && body == player.nearest &&
                           Vector3Distance(top, body->toLocal((player.position - world.renderOrigin).toFloat())) < 3000.0f;
        // Irrigation by moving water (a channel dug to the field, water carried in a bucket).
        if (g.rate < 2.0f && SimTime() - g.wetChecked > 2.0) { g.wetChecked = SimTime(); g.flowing = Fluid::StandingWaterNear(*body, top, 4.0f); }
        const float rate = std::max(g.rate, g.flowing ? 2.0f : 0.0f) * (rainy ? 1.5f : 1.0f);
        g.growth = std::min(1.0f, g.growth + dt * rate / c->growSeconds);
        ++it;
    }
}

bool Farm_SowFor(Inventory& inv, Planetoid& body, Plot& plot, int part, const std::string& seed, std::string* why) {
    LoadCrops();
    if (!IsFarmland(plot, part)) { if (why) *why = "Not farmland"; return false; }
    const std::string key = Key(body, plot, part);
    if (fields.count(key)) { if (why) *why = "Something is growing here already"; return false; }
    for (const CropDef& c : crops) {
        if (!seed.empty() && c.seed != seed) continue;
        if (inv.count(c.seed) <= 0 || !inv.remove(c.seed, 1)) continue;
        static const bool instant = getenv("RAYLIBGAME_FARM_INSTANT") != nullptr; // (#335 tests) sown ripe
        fields[key] = { c.id, instant ? 1.0f : 0.0f };
        return true;
    }
    if (why) *why = "No seeds (harvest bushes for them)";
    return false;
}
bool Farm_Plant(Planetoid& body, Plot& plot, int part, const std::string& seed) {
    if (gFarmRequest) { // joined (#335): the server sows from the inventory it keeps; its fields come back
        if (!IsFarmland(plot, part)) return false;
        if (fields.count(Key(body, plot, part))) { GameUI_Toast("Something is growing here already"); return false; }
        gFarmRequest(body, plot, part, false, seed);
        Audio_Play(Sfx::Pickup, 0.4f, 0.7f);
        return true;
    }
    std::string why;
    if (!Farm_SowFor(inventory, body, plot, part, seed, &why)) {
        if (!why.empty()) GameUI_Toast(why.c_str());
        if (why.rfind("No seeds", 0) == 0) Audio_Play(Sfx::Error, 0.3f);
        return false;
    }
    QuestBook_Count("plant");
    if (gSkillEvent) gSkillEvent("plant", 1.0, Vector3{}); // (#919) Farming XP
    Audio_Play(Sfx::Pickup, 0.4f, 0.7f);
    return true;
}
int Farm_HarvestFor(Inventory& inv, Planetoid& body, Plot& plot, int part, std::string* why, bool local) {
    auto it = fields.find(Key(body, plot, part));
    if (it == fields.end()) return 0;
    const CropDef* c = Crop(it->second.crop);
    if (!c) { fields.erase(it); return 0; }
    if (it->second.growth < 1.0f) { if (why) *why = TextFormat("The %s is not ripe yet (%.0f%%)", c->id.c_str(), it->second.growth * 100.0f); return 0; }
    // (#919) Farming's perks: more crop and seed (the local player's; a server's remote players: #920)
    const int seeds = local ? SkillScaled(1 + GameRng_Int(0, 1), "seed_yield") : 1 + GameRng_Int(0, 1);
    const int crops = local ? SkillScaled(c->yieldCount, "crop_yield") : c->yieldCount;
    inv.add(c->yield, crops);
    inv.add(c->seed, seeds);
    fields.erase(it);
    return crops + seeds;
}
int Farm_Harvest(Planetoid& body, Plot& plot, int part) {
    if (gFarmRequest) { // joined (#335)
        auto it = fields.find(Key(body, plot, part));
        if (it == fields.end()) return 0;
        if (it->second.growth < 1.0f) { GameUI_Toast(TextFormat("The %s is not ripe yet (%.0f%%)", it->second.crop.c_str(), it->second.growth * 100.0f)); return 0; }
        gFarmRequest(body, plot, part, true, "");
        Audio_Play(Sfx::Pickup, 0.5f);
        return 1;
    }
    std::string why;
    const int got = Farm_HarvestFor(inventory, body, plot, part, &why, true);
    if (!why.empty()) GameUI_Toast(why.c_str());
    if (got <= 0) return 0;
    QuestBook_Count("harvest");
    if (gSkillEvent) gSkillEvent("harvest", 1.0, Vector3{}); // (#919) Farming XP
    Audio_Play(Sfx::Pickup, 0.5f);
    return got;
}
void Farm_FieldsNear(const Planetoid& body, Vector3 local, float range, std::vector<NetFarmField>& out) {
    out.clear();
    const std::string prefix = body.name + "|";
    for (auto it = fields.lower_bound(prefix); it != fields.end() && it->first.compare(0, prefix.size(), prefix) == 0; ++it) {
        std::string b; uint32_t plotId = 0; int part = -1;
        if (!Parse(it->first, b, plotId, part)) continue;
        const Plot* plot = body.plots ? body.plots->byId(plotId) : nullptr;
        if (!plot || !IsFarmland(*plot, part) || Vector3Distance(TopCentre(*plot, plot->parts[part]), local) > range) continue;
        out.push_back({ plotId, part, it->second.crop, it->second.growth });
    }
}
void Farm_SetRemoteFields(const Planetoid& body, const std::vector<NetFarmField>& in) {
    const std::string prefix = body.name + "|";
    for (auto it = fields.lower_bound(prefix); it != fields.end() && it->first.compare(0, prefix.size(), prefix) == 0;) it = fields.erase(it);
    for (const NetFarmField& f : in) fields[prefix + std::to_string(f.plot) + "|" + std::to_string(f.part)] = { f.crop, f.growth };
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

bool Farm_TillPart(Planetoid& body, Vector3 ground, Plot*& plot, Part& part) {
    part = Part{ PartType::Farmland, { 0, 0, 0 }, 0, 0 };
    plot = body.plots ? body.plots->plotFor(ground, 6.0f, false) : nullptr;
    if (plot) {
        part.origin = plot->cellAt(plot->toPlot(Vector3Add(ground, Vector3Scale(plot->up, 0.1f))));
        part.origin.x &= ~1; part.origin.z &= ~1; // on a 1 m grid, so tilled squares line up
        return plot->canPlace(part);
    }
    return true; // a new plot there (anchored at the ground), the part at its origin
}
bool Farm_Till(Scene& world, Planetoid& body, Vector3 ground) {
    (void)world;
    if (gFarmTillRequest) return gFarmTillRequest(body, ground); // joined (#335): the server places the farmland
    Plot* plot = nullptr; Part part;
    if (!Farm_TillPart(body, ground, plot, part)) return false;
    if (!plot) {
        if (!body.plots) body.plots = std::make_unique<PlotSet>();
        plot = body.plots->plotFor(ground, 0.0f, true);
    }
    if (!plot->canPlace(part)) return false;
    if (!BuildTool_Place(body, *plot, part, false)) return false;
    QuestBook_Count("till");
    if (gSkillEvent) gSkillEvent("till", 1.0, Vector3{}); // (#919) Farming XP
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
