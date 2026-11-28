#include "ruincache.h"
#include "buildtool.h"
#include "cacheopeners.h"
#include "inventory.h"
#include "worldtypes.h"
#include "research.h"
#include "uipalette.h"
#include "voxel.h"

#include <raylib.h>
#include <algorithm>
#include <cmath>

namespace {
std::vector<std::string> shown;   // the result lines on screen
CacheKind shownKind = CacheKind::Small;
double shownAt = -100.0;          // GetTime() when they were shown
constexpr double kShowSeconds = 7.0;
float hold = -1.0f;               // the hold ring's progress, <0 none
double holdAt = -100.0;           // when it was last set: a hold not kept up (E let go, looked away) stops drawing

uint64_t Mix(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdull; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull; x ^= x >> 33;
    return x;
}
} // namespace

const char* CacheKindName(CacheKind k) {
    switch (k) {
    case CacheKind::Vault: return "sealed vault";
    case CacheKind::DeepVault: return "deep vault";
    default: return "sealed cache";
    }
}

std::string CacheMissing(const Planetoid& body, const Plot& plot, int part, const Inventory& inv) {
    if (part < 0 || part >= (int)plot.parts.size()) return "";
    const CacheOpener* o = CacheOpenerOf(plot.parts[part]);
    if (!o) return "Sealed: nothing opens it"; // an empty or broken cache_openers.txt (logged at load)
    switch (o->kind) {
    case OpenerKind::Power: {
        const float charge = CacheOpener_Charge(body, plot, part);
        if (charge >= 1.0f) return "";
        if (charge > 0.0f) return TextFormat("Charging %d%%: %s", (int)(charge * 100.0f), CacheOpener_NeedText(*o).c_str());
        return "Sealed: its lock " + CacheOpener_NeedText(*o);
    }
    case OpenerKind::Tool:
    case OpenerKind::Part: return inv.count(o->item) > 0 ? "" : "Sealed: it " + CacheOpener_NeedText(*o);
    }
    return "";
}

uint64_t CacheId(const Planetoid& body, const Plot& plot, const Part& p) {
    // Where the cache's origin cell is on its body, to 0.5 m: ruins are built from the seed in the same place everywhere.
    const Vector3 at = plot.toBody({ (p.origin.x + 0.5f) * kCellSize, (p.origin.y + 0.5f) * kCellSize, (p.origin.z + 0.5f) * kCellSize });
    const PlanetField* pf = dynamic_cast<const PlanetField*>(body.field.get());
    const uint64_t seed = pf ? (uint64_t)pf->foliageSeed : 0;
    uint64_t h = Mix(seed ^ 0xCAC4E5ull);
    for (float v : { at.x, at.y, at.z }) h = Mix(h ^ (uint64_t)(int64_t)lroundf(v * 2.0f));
    return h;
}

CacheTeam CacheTeamOf(int team) {
    CacheTeam t;
    const ResearchState& rs = Research_Team(team);
    t.researched.assign(rs.done.begin(), rs.done.end());
    for (const std::string& id : rs.done)
        if (const ResearchNode* n = Research_Find(id)) t.tier = std::max(t.tier, n->tier);
    return t;
}

bool Cache_Open(Planetoid& body, Plot& plot, int part, Inventory& inv, int team, std::vector<CacheRoll>& rolls,
                std::vector<std::pair<std::string, int>>& overflow, std::vector<PartEdit>* edits, std::string* why) {
    auto no = [&](const std::string& w) { if (why) *why = w; return false; };
    if (part < 0 || part >= (int)plot.parts.size() || plot.parts[part].type != PartType::RuinCache) return no("");
    const Part p = plot.parts[part];
    if (CacheOpened(p)) return no("Opened already: empty");
    const CacheKind kind = CacheKindOf(p);
    const std::string missing = CacheMissing(body, plot, part, inv);
    if (!missing.empty()) return no(missing);
    if (const CacheOpener* o = CacheOpenerOf(p)) { // a tool is kept
        if (o->kind == OpenerKind::Power) CacheOpener_Spend(body, plot, part); // the lock's charge, emptied
        else if (o->kind == OpenerKind::Part) inv.remove(o->item, 1);
    }
    const int band = 0; // surface ruins only so far: the depth bands come with the underground templates (#1468)
    const PlanetField* pf = dynamic_cast<const PlanetField*>(body.field.get());
    rolls = Caches_Roll(pf ? (uint64_t)pf->foliageSeed : 0, CacheId(body, plot, p), kind, band, CacheTeamOf(team));
    for (const CacheRoll& r : rolls) {
        if (r.kind != CacheReward::Parts) continue; // blueprints: with the team-library slice; alts and boosts: ApplyTeam
        for (const auto& [item, n] : r.items) {
            const double left = inv.addAmount(item, n);
            if (left > 0.0) overflow.push_back({ item, (int)ceil(left) });
        }
    }
    Caches_ApplyTeam(team, rolls);
    Part q = p;
    q.material = (uint8_t)(p.material | kCacheOpened);
    plot.remove(part);
    if (edits) edits->push_back({ true, part, {} });
    const int at = plot.place(q);
    if (edits) edits->push_back({ false, at, q });
    return true;
}

void RuinCache_ShowResult(CacheKind kind, const std::vector<std::string>& lines) { shown = lines; shownKind = kind; shownAt = GetTime(); }
void RuinCache_SetHold(float progress) { hold = progress; holdAt = GetTime(); }
const std::vector<std::string>& RuinCache_LastResult() { return shown; }

void RuinCache_DrawHUD() {
    const int sw = GetScreenWidth(), sh = GetScreenHeight();
    if (hold >= 0.0f && GetTime() - holdAt < 0.2) { // the hold ring round the crosshair
        const Vector2 c{ sw * 0.5f, sh * 0.5f };
        DrawRing(c, 22.0f, 27.0f, 0.0f, 360.0f, 48, Fade(BLACK, 0.45f));
        DrawRing(c, 22.0f, 27.0f, -90.0f, -90.0f + 360.0f * std::min(hold, 1.0f), 48, Ui().warn);
    }
    const double age = GetTime() - shownAt;
    if (shown.empty() || age > kShowSeconds) return;
    const float a = (float)std::min(1.0, (kShowSeconds - age) / 0.8); // fades out over its last 0.8 s
    const int fs = 20, pad = 14, lh = fs + 8;
    std::string title = std::string("Opened the ") + CacheKindName(shownKind);
    int w = MeasureText(title.c_str(), fs + 4);
    for (const std::string& l : shown) w = std::max(w, MeasureText(l.c_str(), fs));
    const int h = pad * 2 + (fs + 12) + lh * (int)shown.size();
    const int x = (sw - w) / 2 - pad, y = (int)(sh * 0.30f);
    DrawRectangle(x, y, w + pad * 2, h, Fade(BLACK, 0.70f * a));
    DrawRectangleLines(x, y, w + pad * 2, h, Fade(Ui().warn, a));
    DrawText(title.c_str(), x + pad, y + pad, fs + 4, Fade(Ui().warn, a));
    for (size_t i = 0; i < shown.size(); ++i)
        DrawText(shown[i].c_str(), x + pad, y + pad + fs + 12 + lh * (int)i, fs, Fade(RAYWHITE, a));
}
