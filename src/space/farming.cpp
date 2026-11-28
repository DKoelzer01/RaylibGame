#include "farming.h"
#include "uipalette.h"
#include "tools.h"
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
    struct CropDef { std::string id, seed, yield; int yieldCount = 1; float growSeconds = 3600.0f; Color color = GREEN; float spacing = 0.5f; };
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
            float secs = 3600.0f, spacing = 0.5f;
            // id | seed item | yield item count | daylight seconds to ripen | colour r,g,b [| spacing m]
            const int got = sscanf(line, " %63[^ |] | %63[^ |] | %63s %d | %f | %d , %d , %d | %f", id, seed, yield, &n, &secs, &r, &g, &b, &spacing);
            if (got >= 8)
                crops.push_back({ id, seed, yield, n, secs, { (unsigned char)r, (unsigned char)g, (unsigned char)b, 255 }, std::clamp(spacing, 0.2f, 4.0f) });
        }
        fclose(f);
    }
    const CropDef* Crop(const std::string& id) { LoadCrops(); for (const CropDef& c : crops) if (c.id == id) return &c; return nullptr; }

    // Crops stand on farmland ground (#1085: a state of dirt, painted by a till, no grid): keyed by their body and their
    // body-local point on a 5 cm lattice (host and clients agree on it), "<body name>|<x>|<y>|<z>" -> crop and growth.
    struct Growing { std::string crop; float growth = 0.0f; float rate = 0.0f; double wetChecked = -1e9; bool flowing = false; double soilChecked = -1e9; }; // rate: 1, or 2 by static water (worked out once; 0 = not yet); flowing: moving water within 4 m (checked every 2 s)
    std::map<std::string, Growing> fields;
    // Old saves (before #1085): "<body>|<plot id>|<part>" keys of farmland parts, placed on the ground by Farm_LoadDone.
    std::vector<std::pair<std::string, Growing>> legacy;
    constexpr float kGrid = 0.05f;    // m: the crop points' lattice
    constexpr float kAimReach = 0.4f; // m: a crop is the one aimed at within this of the aimed ground

    Vector3 Snap(Vector3 p) { return { roundf(p.x / kGrid) * kGrid, roundf(p.y / kGrid) * kGrid, roundf(p.z / kGrid) * kGrid }; }
    std::string Key(const Planetoid& body, Vector3 local) {
        return body.name + "|" + std::to_string(lroundf(local.x / kGrid)) + "|" + std::to_string(lroundf(local.y / kGrid)) + "|" + std::to_string(lroundf(local.z / kGrid));
    }
    std::string Prefix(const Planetoid& body) { return body.name + "|"; }
    // The parts of a key between '|'s.
    std::vector<std::string> Split(const std::string& key) {
        std::vector<std::string> out;
        size_t a = 0;
        for (size_t b; (b = key.find('|', a)) != std::string::npos; a = b + 1) out.push_back(key.substr(a, b - a));
        out.push_back(key.substr(a));
        return out;
    }
    bool Parse(const std::string& key, std::string& body, Vector3& pos) {
        const std::vector<std::string> s = Split(key);
        if (s.size() != 4) return false;
        body = s[0];
        pos = { std::stol(s[1]) * kGrid, std::stol(s[2]) * kGrid, std::stol(s[3]) * kGrid };
        return true;
    }
    Planetoid* BodyNamed(const std::string& name) {
        for (Planetoid* p : Planetoid::All()) if (p->name == name) return p;
        return nullptr;
    }
    bool OnFarmland(const Planetoid& body, Vector3 local) { return body.materialAt(local) == VoxelMaterial::Farmland; }
    // The crop on `body` nearest `local`, within `range` (fields.end() if none).
    std::map<std::string, Growing>::iterator Nearest(const Planetoid& body, Vector3 local, float range) {
        const std::string prefix = Prefix(body);
        auto best = fields.end();
        float bestD = range;
        for (auto it = fields.lower_bound(prefix); it != fields.end() && it->first.compare(0, prefix.size(), prefix) == 0; ++it) {
            std::string b; Vector3 pos;
            if (!Parse(it->first, b, pos)) continue;
            const float d = Vector3Distance(pos, local);
            if (d <= bestD) { bestD = d; best = it; }
        }
        return best;
    }
    // Why a crop of `c` can't go at `at` on `body` ("" if it can): the ground must be farmland, and no other crop within
    // the larger of the two crops' spacing.
    std::string CantSow(const Planetoid& body, Vector3 at, const CropDef& c) {
        if (!OnFarmland(body, at)) return "Till the ground first (a shovel or hoe, right mouse)";
        const std::string prefix = Prefix(body);
        for (auto it = fields.lower_bound(prefix); it != fields.end() && it->first.compare(0, prefix.size(), prefix) == 0; ++it) {
            std::string b; Vector3 pos;
            if (!Parse(it->first, b, pos)) continue;
            const CropDef* o = Crop(it->second.crop);
            const float need = std::max(c.spacing, o ? o->spacing : 0.5f);
            if (Vector3Distance(pos, at) < need) return TextFormat("Too close to the next crop (%.1f m apart)", need);
        }
        return "";
    }
}

void Farm_Reset() { fields.clear(); legacy.clear(); }

void Farm_Update(Scene& world, float dt) {
    if (dt <= 0.0f || fields.empty()) return;
    if (gNetRole == NetRole::Client) return; // (#335) the server grows them; its fields come as FarmFields
    LoadCrops();
    const Vector3 sun = Vector3Normalize(world.sunDirection);
    for (auto it = fields.begin(); it != fields.end();) {
        std::string bodyName; Vector3 pos;
        Planetoid* body = Parse(it->first, bodyName, pos) ? BodyNamed(bodyName) : nullptr;
        Growing& g = it->second;
        // Dug up, filled over, built on: the ground is no longer farmland (checked every 2 s; #1085).
        if (body && SimTime() - g.soilChecked > 2.0) {
            g.soilChecked = SimTime();
            if (!OnFarmland(*body, pos)) { it = fields.erase(it); continue; }
        }
        const CropDef* c = Crop(g.crop);
        if (!body || !c || g.growth >= 1.0f) { ++it; continue; }
        // Daylight only, twice as fast with water within 4 m (a river, a pond, the sea shore).
        const Vector3 up = Vector3Normalize(pos);
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
                           Vector3Distance(pos, body->toLocalOf(player.position, world.renderOrigin)) < 3000.0f;
        // Irrigation by moving water (a channel dug to the field, water carried in a bucket).
        if (g.rate < 2.0f && SimTime() - g.wetChecked > 2.0) { g.wetChecked = SimTime(); g.flowing = Fluid::StandingWaterNear(*body, pos, 4.0f); }
        const float rate = std::max(g.rate, g.flowing ? 2.0f : 0.0f) * (rainy ? 1.5f : 1.0f);
        g.growth = std::min(1.0f, g.growth + dt * rate / c->growSeconds);
        ++it;
    }
}

bool Farm_SowFor(Inventory& inv, Planetoid& body, Vector3 ground, const std::string& seed, std::string* why) {
    LoadCrops();
    const Vector3 at = Snap(ground);
    std::string refusal = "No seeds (harvest bushes for them)";
    for (const CropDef& c : crops) {
        if (!seed.empty() && c.seed != seed) continue;
        if (inv.count(c.seed) <= 0) continue;
        const std::string cant = CantSow(body, at, c);
        if (!cant.empty()) { refusal = cant; continue; }
        if (!inv.remove(c.seed, 1)) continue;
        static const bool instant = getenv("RAYLIBGAME_FARM_INSTANT") != nullptr; // (#335 tests) sown ripe
        Growing g; g.crop = c.id; g.growth = instant ? 1.0f : 0.0f;
        fields[Key(body, at)] = g;
        return true;
    }
    if (why) *why = refusal;
    return false;
}
bool Farm_Plant(Planetoid& body, Vector3 ground, const std::string& seed) {
    if (gFarmRequest) { // joined (#335): the server sows from the inventory it keeps; its fields come back
        if (!OnFarmland(body, ground)) { GameUI_Toast("Till the ground first (a shovel or hoe, right mouse)"); return false; }
        gFarmRequest(body, Snap(ground), false, seed);
        Audio_Play(Sfx::Pickup, 0.4f, 0.7f);
        return true;
    }
    std::string why;
    if (!Farm_SowFor(inventory, body, ground, seed, &why)) {
        if (!why.empty()) GameUI_Toast(why.c_str());
        Audio_Play(Sfx::Error, 0.3f);
        return false;
    }
    QuestBook_Count("plant");
    if (gSkillEvent) gSkillEvent("plant", 1.0, Vector3{}); // (#919) Farming XP
    Audio_Play(Sfx::Pickup, 0.4f, 0.7f);
    return true;
}
int Farm_HarvestFor(Inventory& inv, Planetoid& body, Vector3 ground, std::string* why, const std::function<float(const char*)>& stat) {
    auto it = Nearest(body, ground, kAimReach);
    if (it == fields.end()) return 0;
    const CropDef* c = Crop(it->second.crop);
    if (!c) { fields.erase(it); return 0; }
    if (it->second.growth < 1.0f) { if (why) *why = TextFormat("The %s is not ripe yet (%.0f%%)", c->id.c_str(), it->second.growth * 100.0f); return 0; }
    // (#919) Farming's perks: more crop and seed (the harvester's: the local player, or #920 a server's remote player)
    const int seeds = stat ? SkillScaled(1 + GameRng_Int(0, 1), "seed_yield", stat) : 1 + GameRng_Int(0, 1);
    const int crops = stat ? SkillScaled(c->yieldCount, "crop_yield", stat) : c->yieldCount;
    inv.add(c->yield, crops);
    inv.add(c->seed, seeds);
    fields.erase(it);
    return crops + seeds;
}
int Farm_Harvest(Planetoid& body, Vector3 ground) {
    if (gFarmRequest) { // joined (#335)
        auto it = Nearest(body, ground, kAimReach);
        if (it == fields.end()) return 0;
        if (it->second.growth < 1.0f) { GameUI_Toast(TextFormat("The %s is not ripe yet (%.0f%%)", it->second.crop.c_str(), it->second.growth * 100.0f)); return 0; }
        gFarmRequest(body, Snap(ground), true, "");
        Audio_Play(Sfx::Pickup, 0.5f);
        return 1;
    }
    std::string why;
    const int got = Farm_HarvestFor(inventory, body, ground, &why, SkillStat);
    if (!why.empty()) GameUI_Toast(why.c_str());
    if (got <= 0) return 0;
    QuestBook_Count("harvest");
    if (gSkillEvent) gSkillEvent("harvest", 1.0, Vector3{}); // (#919) Farming XP
    Audio_Play(Sfx::Pickup, 0.5f);
    return got;
}
void Farm_FieldsNear(const Planetoid& body, Vector3 local, float range, std::vector<NetFarmField>& out) {
    out.clear();
    const std::string prefix = Prefix(body);
    for (auto it = fields.lower_bound(prefix); it != fields.end() && it->first.compare(0, prefix.size(), prefix) == 0; ++it) {
        std::string b; Vector3 pos;
        if (!Parse(it->first, b, pos) || Vector3Distance(pos, local) > range) continue;
        NetFarmField f; f.x = pos.x; f.y = pos.y; f.z = pos.z; f.crop = it->second.crop; f.growth = it->second.growth;
        out.push_back(std::move(f));
    }
}
void Farm_SetRemoteFields(const Planetoid& body, const std::vector<NetFarmField>& in) {
    const std::string prefix = Prefix(body);
    for (auto it = fields.lower_bound(prefix); it != fields.end() && it->first.compare(0, prefix.size(), prefix) == 0;) it = fields.erase(it);
    for (const NetFarmField& f : in) { Growing g; g.crop = f.crop; g.growth = f.growth; fields[Key(body, { f.x, f.y, f.z })] = g; }
}
bool Farm_Crop(const Planetoid& body, Vector3 ground, std::string& crop, float& growth) {
    auto it = Nearest(body, ground, kAimReach);
    if (it == fields.end()) { crop.clear(); growth = 0.0f; return false; }
    crop = it->second.crop;
    growth = it->second.growth;
    return true;
}
bool Farm_Aims(const Planetoid& body, Vector3 ground) { return OnFarmland(body, ground) || Nearest(body, ground, kAimReach) != fields.end(); }
int Farm_CropCount() { return (int)fields.size(); }

void Farm_DrawAimHUD(const Planetoid& body, Vector3 ground) {
    std::string crop; float growth = 0.0f;
    const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40;
    DrawRectangle(x - 8, y - 6, 360, 62, Fade(BLACK, 0.55f));
    if (Farm_Crop(body, ground, crop, growth)) {
        DrawText(TextFormat("farmland: %s %s", crop.c_str(), growth >= 1.0f ? "- ripe" : ""), x, y, 20, Color{ 140, 220, 120, 255 });
        DrawRectangle(x, y + 26, 300, 8, Fade(RAYWHITE, 0.2f));
        DrawRectangle(x, y + 26, (int)(300 * growth), 8, growth >= 1.0f ? Ui().warn : Ui().good);
        DrawText(growth >= 1.0f ? "O harvest" : "grows in daylight; water nearby doubles it", x, y + 38, 14, Fade(RAYWHITE, 0.7f));
    } else {
        DrawText("farmland", x, y, 20, Color{ 200, 170, 120, 255 });
        DrawText("I plant seeds you carry", x, y + 26, 16, Fade(RAYWHITE, 0.7f));
    }
}

// Farmland is a state of dirt (#1085): a till paints the soil around the aimed ground (a free brush, no grid) Farmland in
// the terrain's material channel - saved and replicated as any terrain edit, and grass there goes.
bool Farm_Till(Scene& world, Planetoid& body, Vector3 ground, float radius) {
    (void)world;
    if (!IsTillable(body.materialAt(ground))) return false; // rock, sand, gravel, snow: nothing to till
    if (gFarmTillRequest) return gFarmTillRequest(body, ground, radius); // joined (#335): the server tills
    TerrainEdit e{};
    e.till = true; e.center = ground; e.radius = radius;
    ApplyTerrainEdit(body, e);
    QuestBook_Count("till");
    if (gSkillEvent) gSkillEvent("till", 1.0, Vector3{}); // (#919) Farming XP
    return true;
}

void Farm_Draw3D(Scene& world) {
    if (fields.empty()) return;
    LoadCrops();
    for (const auto& kv : fields) {
        std::string bodyName; Vector3 pos;
        Planetoid* body = Parse(kv.first, bodyName, pos) ? BodyNamed(bodyName) : nullptr;
        const CropDef* c = body ? Crop(kv.second.crop) : nullptr;
        if (!c) continue;
        const Vector3 base = body->toWorld(pos);
        if (Vector3Distance(base, world.camera.position) > 120.0f) continue;
        if (!world.inView(base, 0.9f)) continue; // (#2140) a plant under 0.8 m tall, a quarter metre across
        // A frame on the ground: up from the body's centre, east and south across it (right-handed, so the cubes are not
        // turned inside out under backface culling).
        const Vector3 ul = Vector3Normalize(pos);
        Vector3 el = Vector3CrossProduct({ 0, 1, 0 }, ul);
        el = Vector3Length(el) < 1e-3f ? Vector3{ 1, 0, 0 } : Vector3Normalize(el);
        const Vector3 sl = Vector3CrossProduct(el, ul);
        const Vector3 e = body->dirToWorld(el), u = body->dirToWorld(ul), s = body->dirToWorld(sl);
        const float g = kv.second.growth;
        const int stage = std::min(3, (int)(g * 4.0f));
        // Four stages: sprouts, leafy, tall, ripe (heads in the crop's colour on top).
        static const float kH[4] = { 0.08f, 0.25f, 0.5f, 0.7f };
        const bool ripe = g >= 1.0f;
        const float h = kH[ripe ? 3 : std::min(stage, 2)];
        const Color col = ColorLerp({ 70, 150, 50, 255 }, { 110, 150, 60, 255 }, g);
        Matrix m = { e.x, u.x, s.x, base.x, e.y, u.y, s.y, base.y, e.z, u.z, s.z, base.z, 0, 0, 0, 1 };
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(m));
        const float w = ripe ? 0.05f : 0.035f + 0.015f * std::min(stage, 2);
        // One plant: five stalks in a tuft a quarter metre across (crops stand at least their spacing apart).
        static const float kTuft[5][2] = { { 0.0f, 0.0f }, { 0.09f, 0.04f }, { -0.07f, 0.08f }, { -0.05f, -0.09f }, { 0.06f, -0.07f } };
        for (int k = 0; k < 5; ++k) {
            const float hk = h * (0.85f + 0.06f * (k % 3));
            DrawCube({ kTuft[k][0], hk * 0.5f, kTuft[k][1] }, w, hk, w, col);
            if (ripe) DrawCube({ kTuft[k][0], hk + 0.05f, kTuft[k][1] }, w * 1.8f, 0.1f, w * 1.8f, c->color); // the head
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
    if (fscanf(f, " %255s %63s %f", key, crop, &g) == 3) {
        Growing gr; gr.crop = crop; gr.growth = std::clamp(g, 0.0f, 1.0f);
        if (Split(key).size() == 3) legacy.push_back({ key, gr }); // a farmland part's (before #1085): Farm_LoadDone places it
        else fields[key] = gr;
    }
    return true;
}

// Old saves (#1085): every farmland part becomes tilled ground under where it stood (a patch about its 1 m square), its
// crop at the patch's centre, and the part goes. Runs once the bodies, their terrain edits and their plots are loaded.
void Farm_LoadDone() {
    std::map<std::string, std::pair<Planetoid*, Vector3>> placed; // old key -> body and the crop's ground
    for (Planetoid* body : Planetoid::All()) {
        if (!body->plots) continue;
        for (auto& pl : body->plots->plots)
            for (int i = 0; i < (int)pl->parts.size(); ++i) {
                const Part& part = pl->parts[i];
                if (part.type != PartType::Farmland) continue;
                // The part's footprint centre at its bottom (2 x 2 cells), then the real ground in that column: searched
                // from 1.5 m over it downwards (the part stood on the density surface).
                const Vector3 bottom = pl->toBody({ (part.origin.x + 1.0f) * kCellSize, part.origin.y * kCellSize, (part.origin.z + 1.0f) * kCellSize });
                const Vector3 d = Vector3Normalize(bottom);
                const float r0 = Vector3Length(bottom);
                Vector3 ground = bottom;
                for (float r = r0 + 1.5f; r >= r0 - 1.5f; r -= 0.05f)
                    if (body->densityAt(Vector3Scale(d, r)) >= ISO_LEVEL) { ground = Vector3Scale(d, r + 0.025f); break; }
                TerrainEdit e{};
                e.till = true; e.center = ground; e.radius = 0.75f;
                ApplyTerrainEdit(*body, e);
                pl->remove(i);
                placed[body->name + "|" + std::to_string(pl->id) + "|" + std::to_string(i)] = { body, ground };
            }
    }
    int moved = 0;
    for (auto& [key, g] : legacy) {
        auto it = placed.find(key);
        if (it == placed.end()) continue; // its part was gone already
        fields[Key(*it->second.first, Snap(it->second.second))] = g;
        ++moved;
    }
    if (!placed.empty()) printf("farming: %zu farmland parts became tilled ground (%d crops kept) (#1085)\n", placed.size(), moved);
    legacy.clear();
}
