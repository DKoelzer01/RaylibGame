#include "statusbar.h"
#include "gameui.h"
#include "itemicons.h"
#include "simclock.h"

#include <algorithm>
#include <string>
#include <vector>

namespace {
struct Entry {
    std::string id, label, icon;
    Color color;
    double start = 0.0, seconds = 0.0;
};
std::vector<Entry> gEntries; // in the order they were first added: the stack's order, top to bottom

// Drops the entries that ran out, and any started "after now" (the sim clock went back: a new world).
void Prune() {
    const double now = SimTime();
    gEntries.erase(std::remove_if(gEntries.begin(), gEntries.end(), [&](const Entry& e) { return now >= e.start + e.seconds || now < e.start; }),
                   gEntries.end());
}
Entry* Find(const char* id) {
    for (Entry& e : gEntries) if (e.id == id) return &e;
    return nullptr;
}
}

void StatusBar_Add(const char* id, float seconds, const char* label, Color color, const char* icon) {
    if (!id || seconds <= 0.0f) return;
    Prune();
    Entry* e = Find(id);
    if (!e) { gEntries.push_back(Entry{}); e = &gEntries.back(); e->id = id; }
    e->label = label ? label : id;
    e->icon = icon ? icon : "";
    e->color = color;
    e->start = SimTime();
    e->seconds = seconds;
}

void StatusBar_Remove(const char* id) {
    if (!id) return;
    gEntries.erase(std::remove_if(gEntries.begin(), gEntries.end(), [&](const Entry& e) { return e.id == id; }), gEntries.end());
}

bool StatusBar_Active(const char* id) { Prune(); return id && Find(id) != nullptr; }

float StatusBar_Remaining(const char* id) {
    Prune();
    const Entry* e = id ? Find(id) : nullptr;
    if (!e) return 0.0f;
    return std::clamp((float)(1.0 - (SimTime() - e->start) / e->seconds), 0.0f, 1.0f);
}

int StatusBar_Count() { Prune(); return (int)gEntries.size(); }

void StatusBar_Clear() { gEntries.clear(); }

void StatusBar_Draw() {
    Prune();
    if (gEntries.empty()) return;
    const float s = GameUI_HudScale();
    const int fs = std::max(10, (int)(16 * s));
    const int icon = fs + 4, barW = (int)(160 * s), barH = std::max(4, (int)(8 * s)), gap = (int)(10 * s);
    const int rowH = std::max(icon, fs) + 3 + barH, x = (int)(16 * s);
    const int total = (int)gEntries.size() * rowH + ((int)gEntries.size() - 1) * gap;
    int y = (GetScreenHeight() - total) / 2; // the stack is centred on the left edge's middle
    const double now = SimTime();
    for (const Entry& e : gEntries) {
        const float left = std::clamp((float)(1.0 - (now - e.start) / e.seconds), 0.0f, 1.0f);
        int tx = x;
        if (!e.icon.empty() && ItemIcons_Has(e.icon)) {
            ItemIcons_Draw(e.icon, Rectangle{ (float)x, (float)y, (float)icon, (float)icon });
            tx += icon + (int)(4 * s);
        }
        const int ty = y + (std::max(icon, fs) - fs) / 2;
        DrawText(e.label.c_str(), tx + 1, ty + 1, fs, Fade(BLACK, 0.7f)); // the label carries the meaning, not the colour
        DrawText(e.label.c_str(), tx, ty, fs, RAYWHITE);
        const int by = y + std::max(icon, fs) + 3;
        DrawRectangle(x - 1, by - 1, barW + 2, barH + 2, Fade(BLACK, 0.6f));
        DrawRectangle(x, by, (int)(barW * left + 0.5f), barH, e.color); // drains toward the left edge
        y += rowH + gap;
    }
}
