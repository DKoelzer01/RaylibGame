#include "navigation.h"
#include "survival.h"
#include "animals.h"
#include "farsphere.h"
#include "gameui.h"
#include "inventory.h"
#include "jobsystem.h"
#include "player.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

// ---------------------------------------------------------------------------------------------
// Exploration: 256 m tiles on a cube-sphere grid
// ---------------------------------------------------------------------------------------------

constexpr float kTileSize = 256.0f;
constexpr float kRevealRadius = 150.0f;
const Vector3 kN[6] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
const Vector3 kU[6] = { { 0, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 1 }, { 1, 0, 0 }, { 1, 0, 0 } };
const Vector3 kV[6] = { { 0, 0, 1 }, { 0, 0, -1 }, { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 } };

uint64_t TileKey(Vector3 dir, float radius) {
    int face = 0;
    float best = -2.0f;
    for (int f = 0; f < 6; ++f) { const float d = Vector3DotProduct(dir, kN[f]); if (d > best) { best = d; face = f; } }
    const int n = std::max(1, (int)ceilf(0.5f * PI * radius / kTileSize));
    const float da = 0.5f * PI / n;
    const float a = atanf(Vector3DotProduct(dir, kU[face]) / best), b = atanf(Vector3DotProduct(dir, kV[face]) / best);
    const int i = std::clamp((int)floorf((a + 0.25f * PI) / da), 0, n - 1), j = std::clamp((int)floorf((b + 0.25f * PI) / da), 0, n - 1);
    return ((uint64_t)face << 48) | ((uint64_t)i << 24) | (uint64_t)j;
}

using TileSet = std::unordered_set<uint64_t>;
std::unordered_map<std::string, TileSet> explored;     // body name -> explored tiles (main thread)
std::unordered_map<std::string, int> exploredVersion;  // bumped when a tile is added

struct Waypoint { std::string body, name; Vector3 dir; };
std::vector<Waypoint> waypoints;
std::string startBody;
Vector3 startDir{ 0, 0, 0 };
bool hasStart = false;

// The tangent frame at a body-local unit direction: north toward the spin axis (+Y), east the
// way the surface turns.
bool Frame(Vector3 up, Vector3& east, Vector3& north) {
    const Vector3 c = Vector3CrossProduct({ 0, 1, 0 }, up);
    if (Vector3Length(c) < 0.03f) return false; // within ~2 km of a pole on a 60 km planet
    east = Vector3Normalize(c);
    north = Vector3CrossProduct(up, east);
    return true;
}

// The body under the player, when it has a map (a planet or moon).
struct Here {
    Planetoid* body = nullptr;
    std::shared_ptr<const PlanetField> field;
    Vector3 up{};      // body-local unit direction of the player
    Vector3 forward{}; // body-local view direction
};
Here Locate(Scene& world) {
    Here h;
    Planetoid* p = player.nearest;
    if (!p || p->flat || !p->field) return h;
    auto pf = std::dynamic_pointer_cast<const PlanetField>(p->field);
    if (!pf) return h;
    const Vector3 local = p->toLocal((player.position - world.renderOrigin).toFloat());
    const float d = Vector3Length(local);
    if (d < 1.0f) return h;
    h.body = p;
    h.field = pf;
    h.up = Vector3Scale(local, 1.0f / d);
    h.forward = p->dirToLocal(player.forward);
    return h;
}

// Degrees clockwise from north of a body-local direction seen from `up` (-1 at a pole).
float BearingOf(Vector3 up, Vector3 d) {
    Vector3 east, north;
    if (!Frame(up, east, north)) return -1.0f;
    const float x = Vector3DotProduct(d, east), y = Vector3DotProduct(d, north);
    if (x * x + y * y < 1e-12f) return -1.0f;
    float a = atan2f(x, y) * RAD2DEG;
    if (a < 0.0f) a += 360.0f;
    if (a >= 360.0f) a -= 360.0f; // -0 wraps to 360
    return a;
}

// ---------------------------------------------------------------------------------------------
// Map images, rendered on the workers
// ---------------------------------------------------------------------------------------------

struct MapImage {
    std::shared_ptr<const PlanetField> field;
    std::shared_ptr<const TileSet> tiles;
    Vector3 centre{}, east{}, north{};  // body-local frame of the image centre
    float mpp = 1.0f;                   // metres per pixel
    int W = 0, H = 0;
    bool circle = false;                // alpha 0 outside the inscribed circle (the minimap)
    std::vector<unsigned char> px;
    std::atomic<int> rowsDone{ 0 };
};

unsigned char Enc(float v) { return (unsigned char)std::clamp((int)lroundf(powf(std::max(v, 0.0f), 1.0f / 2.2f) * 255.0f), 0, 255); }

Vector3 PixelDir(const MapImage& m, float x, float y) {
    const float R = m.field->radius;
    const Vector3 p = Vector3Add(Vector3Scale(m.centre, R), Vector3Add(Vector3Scale(m.east, (x - m.W * 0.5f) * m.mpp), Vector3Scale(m.north, (m.H * 0.5f - y) * m.mpp)));
    return Vector3Normalize(p);
}

void RenderRows(MapImage& m, int y0, int y1) {
    const PlanetField& f = *m.field;
    const int W = m.W;
    // Heights for the rows and one above and below (hillshade needs the neighbours).
    const int rows = y1 - y0 + 2;
    std::vector<float> e((size_t)rows * (W + 2));
    auto E = [&](int x, int y) -> float& { return e[(size_t)(y - y0 + 1) * (W + 2) + (x + 1)]; };
    for (int y = y0 - 1; y <= y1; ++y)
        for (int x = -1; x <= W; ++x) E(x, y) = f.elevation(PixelDir(m, x + 0.5f, y + 0.5f), m.mpp);
    const Vector3 light = Vector3Normalize({ -0.6f, 0.75f, 0.6f }); // from the north-west, as maps are lit
    for (int y = y0; y < y1; ++y)
        for (int x = 0; x < W; ++x) {
            unsigned char* p = &m.px[((size_t)y * W + x) * 4];
            const float cx = x + 0.5f - W * 0.5f, cy = y + 0.5f - m.H * 0.5f;
            if (m.circle && cx * cx + cy * cy > (W * 0.5f) * (W * 0.5f)) { p[0] = p[1] = p[2] = p[3] = 0; continue; }
            const Vector3 dir = PixelDir(m, x + 0.5f, y + 0.5f);
            if (!m.tiles->count(TileKey(dir, f.radius))) { p[0] = 16; p[1] = 19; p[2] = 26; p[3] = 235; continue; }
            const float h = E(x, y);
            // Slope and hillshade from the height grid (x east, y south in image space).
            const float dx = (E(x + 1, y) - E(x - 1, y)) / (2.0f * m.mpp), dy = (E(x, y - 1) - E(x, y + 1)) / (2.0f * m.mpp);
            const Vector3 n = Vector3Normalize({ -dx, 1.0f, -dy });
            const float slope = 1.0f - n.y;
            bool water = false;
            Vector3 c = SurfaceAlbedo(f, dir, h, slope, f.foliage, f.foliageSeed, water);
            if (!water && f.inlandWaterNear(dir, 0.0f, m.mpp)) {
                const float w = f.waterElevation(dir, m.mpp);
                if (w != PlanetField::kNoWater && w > h) { c = { 0.04f, 0.16f, 0.30f }; water = true; }
            }
            const float shade = water ? 1.0f : std::clamp(0.55f + 0.75f * Vector3DotProduct(n, light), 0.35f, 1.25f);
            c = Vector3Scale(c, shade * 1.6f); // a little brighter than the lit terrain: a map, not a photo
            p[0] = Enc(c.x); p[1] = Enc(c.y); p[2] = Enc(c.z); p[3] = 255;
        }
    m.rowsDone.fetch_add(y1 - y0);
}

std::shared_ptr<MapImage> StartImage(const Here& h, Vector3 centre, float mpp, int W, int H, bool circle) {
    auto m = std::make_shared<MapImage>();
    m->field = h.field;
    auto it = explored.find(h.body->name);
    m->tiles = std::make_shared<const TileSet>(it != explored.end() ? it->second : TileSet{});
    m->centre = centre;
    if (!Frame(centre, m->east, m->north)) { m->east = { 1, 0, 0 }; m->north = Vector3Normalize(Vector3CrossProduct(centre, m->east)); }
    m->mpp = mpp;
    m->W = W; m->H = H;
    m->circle = circle;
    m->px.assign((size_t)W * H * 4, 0);
    const int block = 16;
    for (int y = 0; y < H; y += block) {
        const int y1 = std::min(H, y + block);
        WorldJobs().submit(circle ? 300.0f : 200.0f, [m, y, y1] { RenderRows(*m, y, y1); });
    }
    return m;
}

// A map texture with the image it came from (its frame, for placing markers).
struct MapView {
    Texture2D tex{};
    std::shared_ptr<MapImage> shown;    // the image on the texture
    std::shared_ptr<MapImage> pending;  // being rendered
    void take() {
        if (!pending || pending->rowsDone.load() < pending->H) return;
        Image img{ pending->px.data(), pending->W, pending->H, 1, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8 };
        if (tex.id) UnloadTexture(tex);
        tex = LoadTextureFromImage(img);
        SetTextureFilter(tex, TEXTURE_FILTER_BILINEAR);
        std::vector<unsigned char>().swap(pending->px);
        shown = pending;
        pending.reset();
    }
    void reset() { if (tex.id) UnloadTexture(tex); tex = Texture2D{}; shown.reset(); pending.reset(); }
};

MapView minimap, bigMap;
std::string minimapBody;
int minimapTilesVersion = -1;
bool mapOpen = false;
Vector3 mapCentre{ 0, 0, 0 };          // body-local unit direction
std::string mapBody;
int mapZoom = 1;                       // index into kZooms (4 m a pixel)
const float kZooms[] = { 2.0f, 4.0f, 8.0f, 16.0f, 32.0f, 64.0f };
constexpr int kMapPixels = 640;
constexpr int kMiniPixels = 176;
constexpr float kMiniMetres = 300.0f;
bool dragging = false, dragMoved = false; // a drag pans; a press and release in place is a click
Vector2 dragFrom{};
double mapDirtyAt = -1.0;
double lastReveal = -1.0;

bool HasMap() { return inventory.count("map") > 0; }

// Screen rectangle of the full map.
Rectangle MapRect() {
    const float s = std::min(GetScreenHeight() * 0.86f, GetScreenWidth() * 0.7f);
    return { (GetScreenWidth() - s) * 0.5f, (GetScreenHeight() - s) * 0.5f + 10.0f, s, s };
}

// A body-local direction on a map image, in its pixels (centre-relative, metres -> pixels).
Vector2 OnImage(const MapImage& m, Vector3 dir) {
    const float R = m.field->radius;
    const Vector3 rel = Vector3Subtract(Vector3Scale(dir, R), Vector3Scale(m.centre, R));
    return { Vector3DotProduct(rel, m.east) / m.mpp, -Vector3DotProduct(rel, m.north) / m.mpp };
}

void DrawMarker(Vector2 at, Color col, const char* label, float size) {
    DrawPoly(at, 4, size, 45.0f, col);
    DrawPolyLines(at, 4, size, 45.0f, BLACK);
    if (label && *label) DrawText(label, (int)(at.x + size + 3), (int)(at.y - 7), 14, RAYWHITE);
}

} // namespace

// ---------------------------------------------------------------------------------------------

void Nav_Reset() {
    explored.clear();
    exploredVersion.clear();
    waypoints.clear();
    hasStart = false;
    minimap.reset();
    bigMap.reset();
    mapOpen = false;
}

void Nav_SetStart(const std::string& body, Vector3 dir) { startBody = body; startDir = Vector3Normalize(dir); hasStart = true; }
bool Nav_MapOpen() { return mapOpen; }
void Nav_CloseMap() { mapOpen = false; dragging = false; }

void Nav_ToggleMap() {
    if (mapOpen) { Nav_CloseMap(); return; }
    if (!HasMap()) { GameUI_Toast("You need a Map (craft one from fibre and charcoal)"); return; }
    mapOpen = true;
    mapCentre = Vector3{ 0, 0, 0 }; // centred on the player when it opens
    mapDirtyAt = 0.0;
}

int Nav_ExploredTiles(const std::string& body) { auto it = explored.find(body); return it == explored.end() ? 0 : (int)it->second.size(); }
int Nav_WaypointCount() { return (int)waypoints.size(); }
bool Nav_MinimapReady() { return minimap.tex.id != 0; }
bool Nav_MapImageReady() { return bigMap.tex.id != 0 && !bigMap.pending; }

void Nav_AddWaypoint(const std::string& body, Vector3 dir, const std::string& name) {
    waypoints.push_back({ body, name, Vector3Normalize(dir) });
    QuestBook_Count("waypoint");
}

float Nav_Heading(Scene& world) {
    const Here h = Locate(world);
    if (!h.body) return -1.0f;
    return BearingOf(h.up, h.forward);
}

float Nav_Bearing(Scene& world, Vector3 dir) {
    const Here h = Locate(world);
    if (!h.body) return -1.0f;
    return BearingOf(h.up, Vector3Subtract(Vector3Normalize(dir), h.up));
}

void Nav_Update(Scene& world, bool inputEnabled) {
    const Here h = Locate(world);
    // Reveal the ground around the player a few times a second.
    if (h.body && GetTime() - lastReveal > 0.3) {
        lastReveal = GetTime();
        TileSet& set = explored[h.body->name];
        const size_t before = set.size();
        Vector3 east, north;
        if (!Frame(h.up, east, north)) { east = { 1, 0, 0 }; north = Vector3Normalize(Vector3CrossProduct(h.up, east)); }
        const float R = h.field->radius;
        for (float a = -kRevealRadius; a <= kRevealRadius; a += 50.0f)
            for (float b = -kRevealRadius; b <= kRevealRadius; b += 50.0f) {
                if (a * a + b * b > kRevealRadius * kRevealRadius) continue;
                const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(h.up, R), Vector3Add(Vector3Scale(east, a), Vector3Scale(north, b))));
                set.insert(TileKey(d, R));
            }
        if (set.size() != before) ++exploredVersion[h.body->name];
        if (set.size() != before) { size_t all = 0; for (const auto& kv : explored) all += kv.second.size(); QuestBook_AtLeast("explored", (double)all); }
    }
    // Minimap: re-rendered when the player has moved 10 m off its centre, or new ground shows.
    minimap.take();
    if (h.body && HasMap() && !minimap.pending) {
        const int ver = exploredVersion[h.body->name];
        const bool moved = !minimap.shown || minimapBody != h.body->name ||
                           Vector3Distance(minimap.shown->centre, h.up) * h.field->radius > 10.0f;
        if (moved || ver != minimapTilesVersion) {
            minimap.pending = StartImage(h, h.up, kMiniMetres / kMiniPixels, kMiniPixels, kMiniPixels, true);
            minimapBody = h.body->name;
            minimapTilesVersion = ver;
        }
    }
    // The full map: pan (drag), zoom (wheel), waypoints (click, right-click removes).
    bigMap.take();
    if (!mapOpen) return;
    if (!h.body || !HasMap()) { Nav_CloseMap(); return; }
    if (mapBody != h.body->name || Vector3LengthSqr(mapCentre) < 0.5f) { mapBody = h.body->name; mapCentre = h.up; mapDirtyAt = 0.0; }
    if (inputEnabled) {
        const Rectangle r = MapRect();
        const float mpp = kZooms[mapZoom], scale = kMapPixels / r.width; // image pixels per screen pixel
        const float wheel = GetMouseWheelMove();
        if (wheel != 0.0f) { mapZoom = std::clamp(mapZoom + (wheel < 0 ? 1 : -1), 0, (int)(sizeof(kZooms) / sizeof(kZooms[0])) - 1); mapDirtyAt = GetTime(); }
        const Vector2 mouse = GetMousePosition();
        const bool inside = CheckCollisionPointRec(mouse, r);
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && inside) { dragging = true; dragMoved = false; dragFrom = mouse; }
        if (dragging && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            const Vector2 d = Vector2Subtract(mouse, dragFrom);
            if (Vector2Length(d) > 4.0f && bigMap.shown) {
                // Move the centre opposite to the drag (the ground follows the cursor).
                const MapImage& m = *bigMap.shown;
                const float R = m.field->radius;
                const Vector3 c = Vector3Add(Vector3Scale(mapCentre, R), Vector3Add(Vector3Scale(m.east, -d.x * scale * mpp), Vector3Scale(m.north, d.y * scale * mpp)));
                mapCentre = Vector3Normalize(c);
                dragFrom = mouse;
                mapDirtyAt = GetTime();
                dragMoved = true;
            }
        }
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT) && dragging) {
            const Vector2 d = Vector2Subtract(mouse, dragFrom);
            if (!dragMoved && Vector2Length(d) <= 4.0f && inside && bigMap.shown && !bigMap.pending && mapDirtyAt < 0.0) { // only on a settled view (the image is the view)
                // A click: a waypoint there.
                const MapImage& m = *bigMap.shown;
                const Vector3 dir = PixelDir(m, (mouse.x - r.x) * scale, (mouse.y - r.y) * scale);
                int n = 1;
                for (const Waypoint& w : waypoints) if (w.body == h.body->name) ++n;
                Nav_AddWaypoint(h.body->name, dir, TextFormat("WP%d", n));
            }
            dragging = false;
        }
        if (IsMouseButtonPressed(MOUSE_BUTTON_RIGHT) && inside && bigMap.shown) {
            // Remove the waypoint nearest the cursor (within 14 px).
            const MapImage& m = *bigMap.shown;
            const float k = m.mpp / mpp; // the image as the view draws it (see Nav_DrawMap)
            const Vector2 off = OnImage(m, mapCentre);
            int best = -1;
            float bestD = 14.0f;
            for (size_t i = 0; i < waypoints.size(); ++i) {
                if (waypoints[i].body != h.body->name) continue;
                const Vector2 p = OnImage(m, waypoints[i].dir);
                const Vector2 s = { r.x + r.width * 0.5f + (p.x - off.x) * k / scale, r.y + r.height * 0.5f + (p.y - off.y) * k / scale };
                const float dd = Vector2Distance(s, mouse);
                if (dd < bestD) { bestD = dd; best = (int)i; }
            }
            if (best >= 0) waypoints.erase(waypoints.begin() + best);
        }
        if (IsKeyPressed(KEY_C)) { mapCentre = h.up; mapDirtyAt = GetTime(); } // back to the player
    }
    // Re-render after the view settles (a short pause, so dragging doesn't queue a render a frame).
    if (mapDirtyAt >= 0.0 && GetTime() - mapDirtyAt > 0.05 && !bigMap.pending) {
        bigMap.pending = StartImage(h, mapCentre, kZooms[mapZoom], kMapPixels, kMapPixels, false);
        mapDirtyAt = -1.0;
    }
}

void Nav_DrawHUD(Scene& world) {
    const Here h = Locate(world);
    if (!h.body) return;
    const int sw = GetScreenWidth();
    // Compass strip: +-90 degrees around the heading.
    const float heading = BearingOf(h.up, h.forward);
    const int cw = 520, ch = 34, cx = sw / 2, cy = 10; // clear of the debug text at the top left
    DrawRectangle(cx - cw / 2, cy, cw, ch, Fade(BLACK, 0.45f));
    DrawRectangleLines(cx - cw / 2, cy, cw, ch, Fade(RAYWHITE, 0.25f));
    if (heading < 0.0f) {
        DrawText("POLE", cx - MeasureText("POLE", 20) / 2, cy + 7, 20, SKYBLUE);
    } else {
        const float pxPerDeg = cw / 180.0f;
        auto xOf = [&](float bearing) {
            float d = bearing - heading;
            while (d > 180.0f) d -= 360.0f;
            while (d < -180.0f) d += 360.0f;
            return d;
        };
        static const char* names[8] = { "N", "NE", "E", "SE", "S", "SW", "W", "NW" };
        for (int b = 0; b < 360; b += 15) {
            const float d = xOf((float)b);
            if (fabsf(d) > 88.0f) continue;
            const int x = cx + (int)(d * pxPerDeg);
            if (b % 45 == 0) {
                const char* t = names[b / 45];
                const int fs = b % 90 == 0 ? 20 : 15;
                DrawText(t, x - MeasureText(t, fs) / 2, cy + (b % 90 == 0 ? 7 : 10), fs, b == 0 ? Color{ 255, 110, 90, 255 } : RAYWHITE);
            } else {
                DrawLine(x, cy + 4, x, cy + 12, Fade(RAYWHITE, 0.6f));
            }
        }
        // Markers: the start and this body's waypoints.
        auto marker = [&](Vector3 dir, Color col, const char* label) {
            const float b = BearingOf(h.up, Vector3Subtract(dir, h.up));
            if (b < 0.0f) return;
            const float d = std::clamp(xOf(b), -88.0f, 88.0f);
            const int x = cx + (int)(d * pxPerDeg);
            DrawTriangle({ (float)x, (float)(cy + ch - 2) }, { (float)x + 6, (float)(cy + ch - 10) }, { (float)x - 6, (float)(cy + ch - 10) }, col);
            const float metres = Vector3Distance(dir, h.up) * h.field->radius;
            DrawText(TextFormat("%s %.0fm", label, metres), x - 20, cy + ch + 2, 12, col);
        };
        if (hasStart && startBody == h.body->name) marker(startDir, GREEN, "start");
        for (const Waypoint& w : waypoints) if (w.body == h.body->name) marker(w.dir, GOLD, w.name.c_str());
        DrawLine(cx, cy, cx, cy + ch, SKYBLUE);
        DrawText(TextFormat("%03.0f", heading), cx - MeasureText("000", 14) / 2, cy + ch + 16, 14, SKYBLUE);
    }
    // Minimap (with a Map): north-up image turned so the view points up; a circle.
    if (!HasMap() || !minimap.tex.id || !minimap.shown || minimapBody != h.body->name) return;
    const MapImage& m = *minimap.shown;
    const float size = 200.0f, rad = size * 0.5f;
    const Vector2 c = { sw - rad - 20.0f, 150.0f + rad };
    const float rot = heading >= 0.0f ? -heading : 0.0f;
    DrawCircleV(c, rad + 4.0f, Fade(BLACK, 0.55f));
    // The image is centred where it was rendered; shift it by how far the player has moved since.
    const Vector2 off = OnImage(m, h.up);
    const float k = size / m.W;
    const float sr = sinf(rot * DEG2RAD), cr = cosf(rot * DEG2RAD);
    const Vector2 shift = { -(off.x * cr - off.y * sr) * k, -(off.x * sr + off.y * cr) * k };
    DrawTexturePro(minimap.tex, { 0, 0, (float)m.W, (float)m.H }, { c.x + shift.x, c.y + shift.y, size, size }, { rad, rad }, rot, WHITE);
    DrawRing(c, rad, rad + 3.0f, 0, 360, 48, Fade(RAYWHITE, 0.6f));
    // Waypoints and the start on it.
    auto onMini = [&](Vector3 dir, Color col) {
        const Vector2 p = OnImage(m, dir);
        const Vector2 rel = { (p.x - off.x) * k, (p.y - off.y) * k };
        Vector2 s = { c.x + rel.x * cr - rel.y * sr, c.y + rel.x * sr + rel.y * cr };
        const Vector2 d = Vector2Subtract(s, c);
        if (Vector2Length(d) > rad - 6.0f) s = Vector2Add(c, Vector2Scale(Vector2Normalize(d), rad - 6.0f)); // pinned to the rim
        DrawMarker(s, col, nullptr, 5.0f);
    };
    if (hasStart && startBody == h.body->name) onMini(startDir, GREEN);
    for (const Waypoint& w : waypoints) if (w.body == h.body->name) onMini(w.dir, GOLD);
    // Animals nearby: small dots (not pinned to the rim), wild ones tan.
    for (int i = 0; i < Animals_Count(); ++i) {
        if (Animals_Body(i) != h.body->name) continue;
        const Vector3 dir = Vector3Normalize(Animals_Position(i));
        const Vector2 p = OnImage(m, dir);
        const Vector2 rel = { (p.x - off.x) * k, (p.y - off.y) * k };
        const Vector2 s = { c.x + rel.x * cr - rel.y * sr, c.y + rel.x * sr + rel.y * cr };
        if (Vector2Distance(s, c) < rad - 4.0f) DrawCircleV(s, 2.5f, Animals_Hostile(i) ? RED : Animals_Wild(i) ? Color{ 230, 180, 110, 255 } : RAYWHITE);
    }
    // The player: an arrow pointing up (the view), north marked on the rim.
    DrawTriangle({ c.x, c.y - 9 }, { c.x - 6, c.y + 7 }, { c.x + 6, c.y + 7 }, SKYBLUE);
    if (heading >= 0.0f) {
        const Vector2 n = { c.x + (rad - 10.0f) * sinf(rot * DEG2RAD), c.y - (rad - 10.0f) * cosf(rot * DEG2RAD) };
        DrawText("N", (int)n.x - 5, (int)n.y - 8, 16, Color{ 255, 110, 90, 255 });
    }
}

void Nav_DrawMap(Scene& world) {
    if (!mapOpen) return;
    const Here h = Locate(world);
    DrawRectangle(0, 0, GetScreenWidth(), GetScreenHeight(), Fade(BLACK, 0.7f));
    const Rectangle r = MapRect();
    DrawRectangleRec(r, { 16, 19, 26, 255 });
    // The last rendered image is drawn where the current view puts it (shifted while dragging, scaled
    // after a zoom), so the map follows the cursor at once and the fresh render only sharpens it.
    if (bigMap.tex.id && bigMap.shown) {
        const MapImage& m = *bigMap.shown;
        const float scale = m.W / r.width, k = m.mpp / kZooms[mapZoom];
        const Vector2 off = Vector3LengthSqr(mapCentre) > 0.5f ? OnImage(m, mapCentre) : Vector2{ 0, 0 };
        const Rectangle dst = { r.x + r.width * 0.5f - (m.W * 0.5f + off.x) * k / scale, r.y + r.height * 0.5f - (m.H * 0.5f + off.y) * k / scale, m.W * k / scale, m.H * k / scale };
        BeginScissorMode((int)r.x, (int)r.y, (int)r.width, (int)r.height);
        DrawTexturePro(bigMap.tex, { 0, 0, (float)m.W, (float)m.H }, dst, { 0, 0 }, 0.0f, WHITE);
        EndScissorMode();
    }
    DrawRectangleLinesEx(r, 2.0f, Fade(RAYWHITE, 0.5f));
    if (bigMap.shown && h.body) {
        const MapImage& m = *bigMap.shown;
        const float scale = m.W / r.width, k = m.mpp / kZooms[mapZoom];
        const Vector2 off = Vector3LengthSqr(mapCentre) > 0.5f ? OnImage(m, mapCentre) : Vector2{ 0, 0 };
        auto toScreen = [&](Vector3 dir) { const Vector2 p = OnImage(m, dir); return Vector2{ r.x + r.width * 0.5f + (p.x - off.x) * k / scale, r.y + r.height * 0.5f + (p.y - off.y) * k / scale }; };
        BeginScissorMode((int)r.x, (int)r.y, (int)r.width, (int)r.height);
        if (hasStart && startBody == h.body->name) DrawMarker(toScreen(startDir), GREEN, "Start", 7.0f);
        for (const Waypoint& w : waypoints) if (w.body == h.body->name) DrawMarker(toScreen(w.dir), GOLD, w.name.c_str(), 7.0f);
        // The player and where they look.
        const Vector2 p = toScreen(h.up);
        const float hd = BearingOf(h.up, h.forward);
        const float a = (hd < 0.0f ? 0.0f : hd) * DEG2RAD;
        const Vector2 f = { sinf(a), -cosf(a) }, sd = { -f.y, f.x };
        DrawTriangle(Vector2Add(p, Vector2Scale(f, 12.0f)), Vector2Subtract(Vector2Subtract(p, Vector2Scale(f, 7.0f)), Vector2Scale(sd, 7.0f)),
                     Vector2Add(Vector2Subtract(p, Vector2Scale(f, 7.0f)), Vector2Scale(sd, 7.0f)), SKYBLUE);
        EndScissorMode();
        // Scale bar: 100 m .. 10 km.
        const float metresPerScreenPx = kZooms[mapZoom] * scale;
        float bar = 100.0f;
        while (bar / metresPerScreenPx < 80.0f) bar *= 2.0f;
        const int bw = (int)(bar / metresPerScreenPx), bx = (int)r.x + 16, by = (int)(r.y + r.height) - 24;
        DrawRectangle(bx - 6, by - 18, bw + 12, 30, Fade(BLACK, 0.5f));
        DrawLine(bx, by, bx + bw, by, RAYWHITE);
        DrawLine(bx, by - 5, bx, by + 5, RAYWHITE);
        DrawLine(bx + bw, by - 5, bx + bw, by + 5, RAYWHITE);
        DrawText(bar >= 1000.0f ? TextFormat("%.1f km", bar / 1000.0f) : TextFormat("%.0f m", bar), bx, by - 16, 14, RAYWHITE);
    }
    if (bigMap.pending) DrawText("drawing...", (int)(r.x + r.width) - 90, (int)r.y + 8, 16, Fade(RAYWHITE, 0.7f));
    DrawText(TextFormat("MAP  %s", h.body ? h.body->name.c_str() : ""), (int)r.x, (int)r.y - 30, 22, RAYWHITE);
    DrawText("M close   wheel zoom   drag pan   click: waypoint   right-click: remove   C: centre on you", (int)r.x, (int)(r.y + r.height) + 8, 16, Fade(RAYWHITE, 0.7f));
}

void Nav_Save(FILE* f) {
    for (const auto& kv : explored) {
        if (kv.second.empty()) continue;
        fprintf(f, "nav_explored %s %zu", kv.first.c_str(), kv.second.size());
        for (uint64_t k : kv.second) fprintf(f, " %llu", (unsigned long long)k);
        fprintf(f, "\n");
    }
    for (const Waypoint& w : waypoints) fprintf(f, "nav_waypoint %s %f %f %f %s\n", w.body.c_str(), w.dir.x, w.dir.y, w.dir.z, w.name.c_str());
    if (hasStart) fprintf(f, "nav_start %s %f %f %f\n", startBody.c_str(), startDir.x, startDir.y, startDir.z);
}

bool Nav_LoadLine(const char* tok, FILE* f) {
    char body[128] = {}, name[64] = {};
    if (strcmp(tok, "nav_explored") == 0) {
        size_t n = 0;
        if (fscanf(f, " %127s %zu", body, &n) != 2) return true;
        TileSet& set = explored[body];
        for (size_t i = 0; i < n; ++i) {
            unsigned long long k = 0;
            if (fscanf(f, " %llu", &k) != 1) break;
            set.insert((uint64_t)k);
        }
        ++exploredVersion[body];
        return true;
    }
    if (strcmp(tok, "nav_waypoint") == 0) {
        Vector3 d{};
        if (fscanf(f, " %127s %f %f %f %63s", body, &d.x, &d.y, &d.z, name) == 5) waypoints.push_back({ body, name, d });
        return true;
    }
    if (strcmp(tok, "nav_start") == 0) {
        Vector3 d{};
        if (fscanf(f, " %127s %f %f %f", body, &d.x, &d.y, &d.z) == 4) Nav_SetStart(body, d);
        return true;
    }
    return false;
}
