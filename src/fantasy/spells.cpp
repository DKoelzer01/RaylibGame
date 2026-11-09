#include "spells.h"
#include "player.h"
#include "world.h"
#include "audio.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>

SpellState spells;

namespace {
    struct Projectile { DVec3 pos; Vector3 vel; float life; float power; Color color; int light; };
    struct TimedLight { int slot; DVec3 pos; double until; bool followPlayer; float range; Color color; };
    struct Beam { Vector3 from, to; double until; };

    std::vector<Projectile> projectiles;
    std::vector<TimedLight> timedLights;
    std::vector<Beam> beams;
    double featherfallUntil = 0.0;
    SealAnalysis preview;
    bool drawing = false;

    Vector2 CanvasCenter() { return { GetScreenWidth() * 0.5f, GetScreenHeight() * 0.5f }; }
    float CanvasRadius() { return 0.38f * (float)std::min(GetScreenWidth(), GetScreenHeight()); }

    Vector3 Forward(const Scene& w) { return Vector3Normalize(Vector3Subtract(w.camera.target, w.camera.position)); }

    // Point lights live in Scene::lights[1..3] (slot 0 is the sun). Returns the slot or -1.
    int AddLight(Scene& w, DVec3 pos, Color c, float range, float seconds, bool follow) {
        int slot = -1;
        for (int i = 1; i < 4; ++i) {
            if (i >= (int)w.lights.size()) { w.lights.push_back(Light{}); slot = i; break; }
            if (!w.lights[i].enabled) { slot = i; break; }
        }
        if (slot < 0) { // all busy: replace the one that expires first
            double soonest = 1e30;
            for (auto& t : timedLights) if (t.until < soonest) { soonest = t.until; slot = t.slot; }
            timedLights.erase(std::remove_if(timedLights.begin(), timedLights.end(), [&](const TimedLight& t) { return t.slot == slot; }), timedLights.end());
        }
        Light& l = w.lights[slot];
        l.type = LIGHT_POINT;
        l.enabled = true;
        l.color = c;
        l.target = { range, 0, 0 };   // lighting.fs reads a point light's range from target.x
        l.position = (pos - w.renderOrigin).toFloat();
        timedLights.push_back({ slot, pos, GetTime() + seconds, follow, range, c });
        return slot;
    }

    bool TerrainHitAhead(Scene& w, float reach, TerrainHit& hit) {
        return RaycastTerrain(w, Ray{ w.camera.position, Forward(w) }, reach, hit);
    }

    void Crater(Planetoid& body, Vector3 renderPoint, float radius) {
        // Density gradient is 0.25/m, so a centre strength of ~0.3*r carves about r metres deep.
        ApplyTerrainEdit(body, TerrainEdit{ body.toLocal(renderPoint), radius / body.scale, -0.3f * radius });
        spells.craters++;
    }

    bool HasModifier(const SealAnalysis& a, const char* sign) {
        for (int m : a.modifiers) if (a.glyphs[m].sign == sign) return true;
        return false;
    }
}

int Spells_ActiveProjectiles() { return (int)projectiles.size(); }

SpellResult Spells_Cast(Scene& w, const std::vector<Stroke>& strokes) {
    SpellResult r;
    r.analysis = AnalyzeSeal(strokes);
    const SealAnalysis& a = r.analysis;
    spells.lastTime = GetTime();
    if (!a.problem.empty()) { r.message = a.problem; Audio_Play(Sfx::Error, 0.5f); spells.last = r; return r; }

    const std::string core = a.glyphs[a.core].sign;
    r.power = std::clamp(a.quality * 1.4f * (1.0f + 0.25f * a.dots), 0.0f, 2.0f);
    if (a.quality < 0.45f) {
        r.name = "Fizzle";
        r.message = TextFormat("The seal wavers and fizzles (ring %.0f%% perfect). Draw a rounder, closed ring.", a.ring.quality * 100.0f);
        spells.health = std::max(0.0f, spells.health - 5.0f);
        spells.fizzles++;
        Audio_Play(Sfx::Error, 0.6f, 0.7f);
        spells.last = r;
        return r;
    }
    r.manaCost = 15.0f + 20.0f * r.power;
    if (spells.mana < r.manaCost) {
        r.message = TextFormat("Not enough mana (%.0f needed)", r.manaCost);
        Audio_Play(Sfx::Error, 0.5f);
        spells.last = r;
        return r;
    }

    const Vector3 fwd = Forward(w);
    const Vector3 eye = w.camera.position;
    Planetoid* body = NearestPlanetoid(w, eye);
    const bool arrow = HasModifier(a, "arrow"), orb = HasModifier(a, "orb"), line = HasModifier(a, "line");
    const float P = r.power;
    TerrainHit hit;
    r.success = true;

    if (core == "flame") {
        if (arrow) {
            r.name = "Fireball";
            Projectile pr{ w.renderOrigin + Vector3Add(eye, Vector3Scale(fwd, 1.5f)), Vector3Scale(fwd, 38.0f), 4.0f, P, ORANGE, -1 };
            pr.light = AddLight(w, pr.pos, ORANGE, 10.0f + 6.0f * P, 4.0f, false);
            projectiles.push_back(pr);
        } else if (orb) {
            r.name = "Flare Burst";
            if (body && TerrainHitAhead(w, 40.0f, hit)) {
                Crater(*hit.planetoid, hit.point, 3.5f + 2.5f * P);
                AddLight(w, w.renderOrigin + Vector3Add(hit.point, Vector3Scale(hit.normal, 2.0f)), ORANGE, 18.0f * P, 1.2f, false);
            } else { r.success = false; r.message = "Flare Burst needs ground within 40 m"; }
        } else {
            r.name = "Ember";
            AddLight(w, w.renderOrigin + eye, { 255, 150, 60, 255 }, 9.0f + 5.0f * P, 25.0f, true);
        }
    } else if (core == "earth") {
        if (line) {
            r.name = "Stone Wall";
            const Vector3 up = player.up;
            const Vector3 side = Vector3Normalize(Vector3CrossProduct(fwd, up));
            const Vector3 flat = Vector3Normalize(Vector3CrossProduct(up, side));
            const Vector3 centre = Vector3Add(eye, Vector3Scale(flat, 6.0f));
            int placed = 0;
            for (int i = -3; i <= 3; ++i) {
                const Vector3 above = Vector3Add(Vector3Add(centre, Vector3Scale(side, i * 1.6f)), Vector3Scale(up, 4.0f));
                if (RaycastTerrain(w, Ray{ above, Vector3Negate(up) }, 14.0f, hit)) {
                    Planetoid& b = *hit.planetoid;
                    // Two stacked brushes per column: a ~3.5 m * P high wall.
                    for (int k = 0; k < 2; ++k)
                        ApplyTerrainEdit(b, TerrainEdit{ b.toLocal(Vector3Add(hit.point, Vector3Scale(up, 1.0f + k * 1.6f * P))), 1.6f, 0.9f });
                    placed++;
                }
            }
            if (!placed) { r.success = false; r.message = "No ground to raise a wall on"; }
        } else if (arrow) {
            r.name = "Wallbreaker";
            int dug = 0;
            for (int i = 1; i <= 8; ++i) {
                const Vector3 p = Vector3Add(eye, Vector3Scale(fwd, i * 1.8f * (0.6f + 0.4f * P)));
                if (Planetoid* b = NearestPlanetoid(w, p)) { ApplyTerrainEdit(*b, TerrainEdit{ b->toLocal(p), 2.2f, -1.2f * P }); dug++; }
            }
            spells.craters += dug ? 1 : 0;
        } else {
            r.name = "Mound";
            if (TerrainHitAhead(w, 30.0f, hit)) ApplyTerrainEdit(*hit.planetoid, TerrainEdit{ hit.planetoid->toLocal(hit.point), 3.0f, 0.8f * P });
            else { r.success = false; r.message = "Mound needs ground within 30 m"; }
        }
    } else if (core == "water") {
        r.name = "Mending";
        spells.health = std::min(100.0f, spells.health + 35.0f * P);
    } else if (core == "wind") {
        player.grounded = false;
        if (arrow) {
            r.name = "Gale Step";
            player.velocity = Vector3Add(player.velocity, Vector3Add(Vector3Scale(fwd, 16.0f * P), Vector3Scale(player.up, 3.0f)));
        } else if (orb) {
            r.name = "Featherfall";
            featherfallUntil = GetTime() + 8.0 * P;
            player.gravityScale = 0.2f;
        } else {
            r.name = "Updraft";
            player.velocity = Vector3Add(player.velocity, Vector3Scale(player.up, 9.0f * P));
        }
    } else if (core == "light") {
        if (orb) {
            r.name = "Lantern Orb";
            const Vector3 at = TerrainHitAhead(w, 30.0f, hit) ? Vector3Add(hit.point, Vector3Scale(hit.normal, 1.5f)) : Vector3Add(eye, Vector3Scale(fwd, 8.0f));
            AddLight(w, w.renderOrigin + at, { 255, 240, 200, 255 }, 14.0f * P, 60.0f, false);
        } else if (arrow) {
            r.name = "Sunlance";
            if (TerrainHitAhead(w, 80.0f, hit)) {
                beams.push_back({ Vector3Add(eye, Vector3Scale(fwd, 1.0f)), hit.point, GetTime() + 0.35 });
                Crater(*hit.planetoid, hit.point, 1.5f + P);
                AddLight(w, w.renderOrigin + Vector3Add(hit.point, Vector3Scale(hit.normal, 1.0f)), { 255, 250, 200, 255 }, 12.0f, 1.0f, false);
            } else { r.success = false; r.message = "Sunlance found nothing to strike"; }
        } else {
            r.name = "Glow";
            AddLight(w, w.renderOrigin + eye, { 230, 240, 255, 255 }, 10.0f + 6.0f * P, 30.0f, true);
        }
    }

    if (r.success) {
        spells.mana -= r.manaCost;
        spells.casts++;
        if (r.message.empty())
            r.message = TextFormat("%s! power %.0f%%  (ring %.0f%%, %s %.0f%%%s)", r.name.c_str(), P * 100.0f, a.ring.quality * 100.0f,
                                   core.c_str(), a.glyphs[a.core].score * 100.0f, a.dots ? TextFormat(", %d dots", a.dots) : "");
        Audio_Play(Sfx::Craft, 0.7f, 0.8f + 0.2f * P);
    } else {
        Audio_Play(Sfx::Error, 0.5f);
    }
    spells.last = r;
    return r;
}

void Spells_Update(Scene& w, float dt, bool inputEnabled) {
    spells.mana = std::min(spells.maxMana, spells.mana + spells.manaRegenPerSec * dt);
    if (featherfallUntil > 0.0 && GetTime() > featherfallUntil) { player.gravityScale = 1.0f; featherfallUntil = 0.0; }

    // Grimoire input.
    if (inputEnabled && IsKeyPressed(KEY_G)) {
        if (!spells.grimoireOpen) { spells.grimoireOpen = true; spells.strokes.clear(); preview = SealAnalysis{}; }
        else { if (!spells.strokes.empty()) Spells_Cast(w, spells.strokes); spells.grimoireOpen = false; }
    }
    if (spells.grimoireOpen) {
        if (IsKeyPressed(KEY_ENTER)) { if (!spells.strokes.empty()) Spells_Cast(w, spells.strokes); spells.grimoireOpen = false; }
        if (IsKeyPressed(KEY_ESCAPE)) spells.grimoireOpen = false;
        if (IsKeyPressed(KEY_C)) { spells.strokes.clear(); preview = SealAnalysis{}; }
        if ((IsMouseButtonPressed(MOUSE_BUTTON_RIGHT) || IsKeyPressed(KEY_BACKSPACE)) && !spells.strokes.empty()) {
            spells.strokes.pop_back();
            preview = AnalyzeSeal(spells.strokes);
        }
        const Vector2 m = GetMousePosition();
        const bool onCanvas = Vector2Distance(m, CanvasCenter()) < CanvasRadius() * 1.08f;
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && onCanvas) { spells.strokes.push_back({ m }); drawing = true; }
        if (drawing && IsMouseButtonDown(MOUSE_BUTTON_LEFT) && Vector2Distance(m, spells.strokes.back().back()) > 2.0f) spells.strokes.back().push_back(m);
        if (drawing && IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) {
            drawing = false;
            if (spells.strokes.back().size() < 2) spells.strokes.back().push_back(Vector2Add(m, { 1, 1 })); // a tap = a dot
            preview = AnalyzeSeal(spells.strokes); // live feedback as you draw
        }
    }

    // Projectiles.
    for (auto it = projectiles.begin(); it != projectiles.end();) {
        Projectile& p = *it;
        p.life -= dt;
        const Vector3 from = (p.pos - w.renderOrigin).toFloat();
        const Vector3 to = Vector3Add(from, Vector3Scale(p.vel, dt));
        p.vel = Vector3Add(p.vel, Vector3Scale(player.up, -2.0f * dt)); // a gentle arc
        bool hitSomething = false;
        TerrainHit hit;
        const float len = Vector3Distance(from, to);
        if (len > 0.0f && RaycastTerrain(w, Ray{ from, Vector3Scale(Vector3Subtract(to, from), 1.0f / len) }, len + 0.5f, hit)) {
            Crater(*hit.planetoid, hit.point, 2.5f + 2.0f * p.power);
            AddLight(w, w.renderOrigin + Vector3Add(hit.point, Vector3Scale(hit.normal, 2.0f)), ORANGE, 16.0f * p.power, 1.0f, false);
            Audio_Play(Sfx::Dig, 1.0f, 0.5f);
            hitSomething = true;
        }
        p.pos = w.renderOrigin + to;
        if (p.light >= 0 && p.light < (int)w.lights.size()) {
            for (auto& t : timedLights) if (t.slot == p.light) t.pos = p.pos;
        }
        if (hitSomething || p.life <= 0.0f) {
            if (p.light >= 0) for (auto& t : timedLights) if (t.slot == p.light) t.until = 0.0; // light dies with it
            it = projectiles.erase(it);
        } else {
            ++it;
        }
    }

    // Timed lights.
    for (auto it = timedLights.begin(); it != timedLights.end();) {
        if (it->slot >= (int)w.lights.size()) { it = timedLights.erase(it); continue; }
        Light& l = w.lights[it->slot];
        if (GetTime() > it->until) { l.enabled = false; it = timedLights.erase(it); continue; }
        if (it->followPlayer) it->pos = w.renderOrigin + Vector3Add((player.position - w.renderOrigin).toFloat(), Vector3Scale(player.up, 2.2f));
        l.position = (it->pos - w.renderOrigin).toFloat();
        ++it;
    }
    beams.erase(std::remove_if(beams.begin(), beams.end(), [](const Beam& b) { return GetTime() > b.until; }), beams.end());
}

void Spells_Draw3D(const Scene& w) {
    for (const Projectile& p : projectiles) {
        const Vector3 at = (p.pos - w.renderOrigin).toFloat();
        DrawSphere(at, 0.35f + 0.15f * p.power, p.color);
        DrawSphereWires(at, 0.6f + 0.2f * p.power, 6, 8, Fade(YELLOW, 0.5f));
    }
    for (const TimedLight& t : timedLights) {
        if (t.followPlayer) continue;
        DrawSphere((t.pos - w.renderOrigin).toFloat(), 0.25f, t.color);
    }
    for (const Beam& b : beams) {
        DrawCylinderEx(b.from, b.to, 0.12f, 0.05f, 6, { 255, 250, 210, 255 });
    }
}

void Spells_DrawOverlay(Scene&) {
    const int W = GetScreenWidth(), H = GetScreenHeight();
    // HUD bars.
    DrawRectangle(20, H - 80, 260, 16, Fade(BLACK, 0.55f));
    DrawRectangle(20, H - 80, (int)(260 * spells.health / 100.0f), 16, { 220, 70, 70, 255 });
    DrawText(TextFormat("HP %3.0f", spells.health), 26, H - 79, 14, RAYWHITE);
    DrawRectangle(20, H - 58, 260, 16, Fade(BLACK, 0.55f));
    DrawRectangle(20, H - 58, (int)(260 * spells.mana / spells.maxMana), 16, { 110, 120, 255, 255 });
    DrawText(TextFormat("MANA %3.0f", spells.mana), 26, H - 57, 14, RAYWHITE);
    if (GetTime() - spells.lastTime < 5.0) {
        const Color c = spells.last.success ? Color{ 180, 220, 255, 255 } : Color{ 255, 190, 120, 255 };
        DrawText(spells.last.message.c_str(), 20, H - 110, 20, c);
    }
    if (!spells.grimoireOpen) {
        DrawText("G grimoire (draw a seal)   LMB dig  RMB build  F fly  F1 debug", 20, H - 30, 18, Fade(RAYWHITE, 0.75f));
        return;
    }

    // Grimoire canvas.
    const Vector2 c = CanvasCenter();
    const float R = CanvasRadius();
    DrawRectangle(0, 0, W, H, Fade(BLACK, 0.35f));
    DrawCircleV(c, R * 1.12f, Fade({ 238, 226, 196, 255 }, 0.95f));
    DrawCircleLinesV(c, R, Fade({ 120, 100, 70, 255 }, 0.25f)); // faint guide: you still draw the ring yourself
    for (const Stroke& s : spells.strokes)
        for (size_t i = 1; i < s.size(); ++i) DrawLineEx(s[i - 1], s[i], 4.0f, { 40, 30, 90, 255 });

    // Live analysis.
    if (preview.ring.found) {
        DrawCircleLinesV(preview.ring.center, preview.ring.radius, Fade(DARKGREEN, 0.6f));
        DrawText(TextFormat("ring %.0f%%  (round %.0f%%, closed %.0f%%)", preview.ring.quality * 100.0f, preview.ring.roundness * 100.0f, preview.ring.closure * 100.0f),
                 (int)(c.x - R), (int)(c.y + R * 1.15f), 18, DARKGREEN);
    }
    for (size_t i = 0; i < preview.glyphs.size(); ++i) {
        const GlyphMatch& g = preview.glyphs[i];
        const bool core = (int)i == preview.core;
        DrawText(TextFormat("%s %.0f%%%s", g.sign.c_str(), g.score * 100.0f, core ? " (core)" : ""), (int)(g.center.x + g.size * 0.4f), (int)(g.center.y - 10), 16,
                 core ? MAROON : DARKBLUE);
    }
    if (preview.dots) DrawText(TextFormat("+%d power dots", preview.dots), (int)(c.x - R), (int)(c.y + R * 1.15f + 22), 16, DARKBLUE);
    if (!preview.problem.empty() && !spells.strokes.empty()) DrawText(preview.problem.c_str(), (int)(c.x - R), (int)(c.y - R * 1.2f), 18, MAROON);
    DrawText("Draw a ring, a core sign in the middle, modifiers near the ring.  G/Enter cast   RMB undo   C clear   Esc close",
             W / 2 - MeasureText("Draw a ring, a core sign in the middle, modifiers near the ring.  G/Enter cast   RMB undo   C clear   Esc close", 18) / 2, 20, 18, RAYWHITE);
    DrawText("core: triangle=flame  square=earth  wave=water  spiral=wind  asterisk=light     modifiers: chevron=arrow  small circle=orb  stroke=line  dots=power",
             W / 2 - MeasureText("core: triangle=flame  square=earth  wave=water  spiral=wind  asterisk=light     modifiers: chevron=arrow  small circle=orb  stroke=line  dots=power", 16) / 2, 46, 16, Fade(RAYWHITE, 0.85f));
}

void Spells_ShowInGrimoire(const std::vector<Stroke>& strokes) {
    spells.grimoireOpen = true;
    spells.strokes = strokes;
    preview = AnalyzeSeal(strokes);
}
