// Creature framework (#1223): see creatures.h.
#include "creatures.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>

namespace {
std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}
std::vector<std::string> SplitOn(const std::string& s, char sep) {
    std::vector<std::string> out;
    size_t start = 0, at;
    while ((at = s.find(sep, start)) != std::string::npos) { out.push_back(Trim(s.substr(start, at - start))); start = at + 1; }
    out.push_back(Trim(s.substr(start)));
    return out;
}
std::string ReadFile(const char* path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return {};
    std::stringstream ss; ss << f.rdbuf();
    return ss.str();
}
bool ToFloat(const std::string& v, float& out) { char* end = nullptr; out = strtof(v.c_str(), &end); return end && end != v.c_str() && *end == '\0'; }
} // namespace

// ---- Factions ----
int FactionTable::find(const std::string& name) const {
    for (size_t i = 0; i < names.size(); ++i) if (names[i] == name) return (int)i;
    return -1;
}

int FactionTable::add(const std::string& name) {
    const int f = find(name);
    if (f >= 0) return f;
    const size_t n = names.size();
    std::vector<Stance> grown((n + 1) * (n + 1), Stance::Neutral);
    for (size_t a = 0; a < n; ++a) for (size_t b = 0; b < n; ++b) grown[a * (n + 1) + b] = m[a * n + b];
    m.swap(grown);
    names.push_back(name);
    return (int)n;
}

Stance FactionTable::stance(int a, int b) const {
    const int n = (int)names.size();
    if (a < 0 || b < 0 || a >= n || b >= n) return Stance::Neutral;
    if (a == b) return Stance::Friendly;
    return m[a * n + b];
}

void FactionTable::set(int a, int b, Stance s) {
    const int n = (int)names.size();
    if (a >= 0 && b >= 0 && a < n && b < n) m[a * n + b] = s;
}

bool Factions_Parse(const std::string& text, FactionTable& out, std::string* error) {
    std::istringstream in(text);
    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        const auto parts = SplitOn(line, '|');
        if (parts[0].empty()) { if (error) *error = "line " + std::to_string(lineNo) + ": no faction"; return false; }
        const int self = out.add(parts[0]);
        for (size_t k = 1; k < parts.size(); ++k) {
            const size_t colon = parts[k].find(':');
            if (colon == std::string::npos) { if (error) *error = "line " + std::to_string(lineNo) + ": '" + parts[k] + "' needs 'stance: factions'"; return false; }
            const std::string what = Trim(parts[k].substr(0, colon));
            Stance st;
            if (what == "hostile") st = Stance::Hostile;
            else if (what == "fears") st = Stance::Fears;
            else if (what == "friendly") st = Stance::Friendly;
            else if (what == "neutral") st = Stance::Neutral;
            else { if (error) *error = "line " + std::to_string(lineNo) + ": unknown stance '" + what + "'"; return false; }
            for (const std::string& other : SplitOn(parts[k].substr(colon + 1), ','))
                if (!other.empty()) { const int o = out.add(other); out.set(self, o, st); }
        }
    }
    return true;
}

const FactionTable& Factions() {
    static FactionTable t;
    static bool loaded = false;
    if (!loaded) {
        loaded = true;
        std::string err;
        if (!Factions_Parse(ReadFile("resources/data/factions.txt"), t, &err)) fprintf(stderr, "factions.txt: %s\n", err.c_str());
        t.add("player");
    }
    return t;
}

int Faction_Player() { return Factions().find("player"); }

// ---- Species data ----
bool Creatures_ParseLine(const std::string& line, CreatureSpec& c, std::string* error) {
    const auto parts = SplitOn(line, '|');
    if (parts.size() < 2 || parts[0].empty() || parts[1].empty()) { if (error) *error = "needs 'id | faction | key=value ...'"; return false; }
    c = CreatureSpec{};
    c.id = parts[0];
    c.faction = parts[1];
    auto fail = [&](const std::string& why) { if (error) *error = c.id + ": " + why; return false; };
    std::string rest;
    for (size_t k = 2; k < parts.size(); ++k) rest += parts[k] + " ";
    std::istringstream ws(rest);
    std::string tok;
    while (ws >> tok) {
        const size_t eq = tok.find('=');
        if (eq == std::string::npos) return fail("'" + tok + "' is not key=value");
        const std::string k = tok.substr(0, eq), v = tok.substr(eq + 1);
        float f = 0.0f;
        auto num = [&](float& dst) { if (!ToFloat(v, f)) return false; dst = f; return true; };
        bool ok = true;
        if (k == "aggro") {
            if (v == "passive") c.aggro = Aggro::Passive; else if (v == "defensive") c.aggro = Aggro::Defensive;
            else if (v == "territorial") c.aggro = Aggro::Territorial; else if (v == "predator") c.aggro = Aggro::Predator; else ok = false;
        } else if (k == "range") ok = num(c.aggroRange);
        else if (k == "chase") ok = num(c.chaseRange);
        else if (k == "angry") ok = num(c.angryTime);
        else if (k == "retreat") ok = num(c.retreatAt);
        else if (k == "threat") ok = num(c.threat);
        else if (k == "sight") ok = num(c.sight);
        else if (k == "darksight") ok = num(c.darkSight);
        else if (k == "hearing") ok = num(c.hearing);
        else if (k == "light") { if (v == "none") c.light = LightPref::None; else if (v == "avoids") c.light = LightPref::Avoids; else if (v == "drawn") c.light = LightPref::Drawn; else ok = false; }
        else if (k == "sleeps") c.sleepsByDay = v == "day" || v == "1";
        else if (k == "fire") c.fireScares = v == "scares" || v == "1";
        else if (k == "firesafe") c.fireSafe = !(v == "0" || v == "no");
        else if (k == "shield") ok = num(c.frontShield);
        else if (k == "shieldarc") ok = num(c.shieldArcDeg);
        else if (k == "wake") { c.wakeByLight = v == "light"; ok = v == "light" || v == "none"; }
        else if (k == "wakerange") ok = num(c.wakeRange);
        else if (k == "rare") { // item:n:chance,...
            for (const std::string& p : SplitOn(v, ',')) {
                const auto f = SplitOn(p, ':');
                float n = 1.0f, ch = 0.0f;
                if (f.size() != 3 || f[0].empty() || !ToFloat(f[1], n) || !ToFloat(f[2], ch)) { ok = false; break; }
                c.rare.push_back({ f[0], (int)n, ch });
            }
        }
        else if (k == "ceiling") c.ceiling = v == "1" || v == "yes";
        else if (k == "drop") ok = num(c.dropScale);
        else if (k == "droprange") ok = num(c.dropRange);
        else if (k == "burst") c.burst = v == "1" || v == "yes";
        else if (k == "group") { float n = 1.0f; ok = num(n); c.group = std::max(1, (int)n); }
        else if (k == "thrower") c.thrower = v == "1" || v == "yes";
        else if (k == "thief") c.thief = v == "1" || v == "yes";
        else if (k == "ranged") c.ranged = v == "1" || v == "yes";
        else if (k == "dart") c.dart = v == "1" || v == "yes";
        else if (k == "dartspeed") ok = num(c.dartSpeed);
        else if (k == "immune") { for (const std::string& s : SplitOn(v, ',')) { if (s == "bleed") c.immuneBleed = true; else if (s == "poison") c.immunePoison = true; else ok = false; } }
        else if (k == "torch") c.torch = v == "1" || v == "yes";
        else if (k == "split") c.split = v == "1" || v == "yes";
        else if (k == "wet") c.wet = v == "1" || v == "yes";
        else if (k == "corrode") { float n = 0.0f; ok = num(n); c.corrode = (int)n; }
        else if (k == "throw") ok = num(c.throwDamage);
        else if (k == "throwmin") ok = num(c.throwMin);
        else if (k == "throwrange") ok = num(c.throwRange);
        else if (k == "throwwindup") ok = num(c.throwWindUp);
        else if (k == "throwcooldown") ok = num(c.throwCooldown);
        else if (k == "damage") ok = num(c.damage);
        else if (k == "reach") ok = num(c.reach);
        else if (k == "windup") ok = num(c.windUp);
        else if (k == "cooldown") ok = num(c.cooldown);
        else if (k == "knock") ok = num(c.knockback);
        else if (k == "poise") ok = num(c.poise);
        else if (k == "bite") c.bite = v == "1" || v == "yes";
        else if (k == "type") {
            ok = false;
            for (int t = 0; t < (int)DamageType::Count; ++t) if (v == DamageType_Name((DamageType)t)) { c.type = (DamageType)t; ok = true; }
        } else if (k == "armour") ok = num(c.resist.armour);
        else if (k == "resist") { // type:share,type:share
            for (const std::string& p : SplitOn(v, ',')) {
                const size_t colon = p.find(':');
                bool found = false; float share = 0.0f;
                if (colon != std::string::npos && ToFloat(p.substr(colon + 1), share))
                    for (int t = 0; t < (int)DamageType::Count; ++t) if (p.substr(0, colon) == DamageType_Name((DamageType)t)) { c.resist.r[t] = share; found = true; }
                if (!found) { ok = false; break; }
            }
        } else if (k == "habitat") {
            if (v == "surface") c.habitat = Habitat::Surface; else if (v == "night") c.habitat = Habitat::Night; else if (v == "cave") c.habitat = Habitat::Cave;
            else if (v == "mouth") c.habitat = Habitat::Mouth; else if (v == "dungeon") c.habitat = Habitat::Dungeon; else ok = false;
        } else if (k == "depth") { // min-max (m)
            const size_t dash = v.find('-');
            ok = dash != std::string::npos && ToFloat(v.substr(0, dash), c.depthMin) && ToFloat(v.substr(dash + 1), c.depthMax);
        } else if (k == "dark") ok = num(c.darkMax);
        else if (k == "cap") { float n = 0.0f; ok = num(n); c.capPerChunk = (int)n; }
        else if (k == "respawn") ok = num(c.respawn);
        else if (k == "nests") c.nests = v == "1" || v == "yes";
        else if (k == "loot") { // item:n,item:n
            for (const std::string& p : SplitOn(v, ',')) {
                const size_t colon = p.find(':');
                float n = 1.0f;
                if (p.empty() || (colon != std::string::npos && !ToFloat(p.substr(colon + 1), n))) { ok = false; break; }
                c.loot.push_back({ p.substr(0, colon), (int)n });
            }
        } else if (k.rfind("sound_", 0) == 0) {
            static const char* kNames[5] = { "idle", "alert", "telegraph", "hurt", "death" };
            ok = false;
            for (int e = 0; e < 5; ++e) if (k.substr(6) == kNames[e]) { c.sounds[e] = v; ok = true; }
        } else return fail("unknown key '" + k + "'");
        if (!ok) return fail("bad value for " + k + ": '" + v + "'");
    }
    return true;
}

bool Creatures_On() {
    static const bool on = !(getenv("RAYLIBGAME_CREATURES") && getenv("RAYLIBGAME_CREATURES")[0] == '0');
    return on;
}

static std::vector<CreatureSpec>& AllMut() {
    static std::vector<CreatureSpec> all;
    static bool loaded = false;
    if (!loaded) {
        loaded = true;
        std::istringstream in(ReadFile("resources/data/creatures.txt"));
        std::string line;
        while (std::getline(in, line)) {
            const std::string t = Trim(line);
            if (t.empty() || t[0] == '#') continue;
            CreatureSpec c; std::string err;
            if (!Creatures_ParseLine(t, c, &err)) { fprintf(stderr, "creatures.txt: %s\n", err.c_str()); continue; }
            c.factionIdx = Factions().find(c.faction);
            if (c.factionIdx < 0) fprintf(stderr, "creatures.txt: %s: faction '%s' is not in factions.txt\n", c.id.c_str(), c.faction.c_str());
            all.push_back(c);
        }
    }
    return all;
}

const std::vector<CreatureSpec>& Creature_All() { return AllMut(); }

bool Creatures_TestAdd(const std::string& line, std::string* error) {
    CreatureSpec c;
    if (!Creatures_ParseLine(line, c, error)) return false;
    c.factionIdx = Factions().find(c.faction);
    for (CreatureSpec& o : AllMut()) if (o.id == c.id) { o = c; return true; }
    AllMut().push_back(c);
    return true;
}

const CreatureSpec* Creature_Spec(const std::string& id) {
    if (!Creatures_On()) return nullptr;
    for (const CreatureSpec& c : Creature_All()) if (c.id == id) return &c;
    return nullptr;
}

const char* CreatureEvent_Name(CreatureEvent e) {
    switch (e) {
        case CreatureEvent::Idle: return "idle";
        case CreatureEvent::Alert: return "alert";
        case CreatureEvent::Telegraph: return "telegraph";
        case CreatureEvent::Hurt: return "hurt";
        case CreatureEvent::Death: return "death";
    }
    return "?";
}

static int gEventCounts[5] = {};
void Creature_Event(const CreatureSpec& c, CreatureEvent e, uint64_t body, Vector3 local) {
    ++gEventCounts[(int)e];
    if (gCreatureEvent) gCreatureEvent(c.id, e, body, local, c.sounds[(int)e]);
}
int Creature_EventCount(CreatureEvent e) { return gEventCounts[(int)e]; }

// ---- Senses ----
float Sense_SightRange(const CreatureSpec& c, float light) {
    const float l = std::clamp(light, 0.0f, 1.0f);
    return c.sight * (c.darkSight + (1.0f - c.darkSight) * l);
}
bool Sense_Sees(const CreatureSpec& c, float dist, float light, bool los) { return los && dist <= Sense_SightRange(c, light); }
bool Sense_Hears(const CreatureSpec& c, float dist, float loudness) { return loudness > 0.0f && dist <= c.hearing * loudness; }

static std::vector<NoiseEvent> gNoise;
void Noise_Emit(uint64_t body, Vector3 local, float loudness, double now, uint64_t source) {
    if (gNoise.size() > 256) gNoise.erase(gNoise.begin(), gNoise.begin() + 128);
    gNoise.push_back({ body, local, loudness, now, source });
}
const std::vector<NoiseEvent>& Noise_Recent(double now) {
    gNoise.erase(std::remove_if(gNoise.begin(), gNoise.end(), [&](const NoiseEvent& e) { return now - e.at > NoiseTune::kKeep || e.at > now + 1.0; }), gNoise.end());
    return gNoise;
}
void Noise_Clear() { gNoise.clear(); }

// ---- The brain ----
const char* Mood_Name(Mood m) {
    static const char* n[] = { "idle", "wander", "alert", "chase", "attack", "retreat", "flee", "sleep" };
    return n[std::clamp((int)m, 0, 7)];
}

BrainOut Brain_Step(BrainState& s, const CreatureSpec& c, const BrainInput& in, float dt) {
    BrainOut o;
    const FactionTable& ft = Factions();
    // (#1227) A light sleeper: asleep until a light comes near or it is hit; then awake a minute after the last of either.
    if (c.wakeByLight) {
        if (in.lightNear || in.hurt) s.awakeUntil = in.now + 60.0;
        if (in.now >= s.awakeUntil) { s.mood = Mood::Sleep; s.windUpAt = -1.0; s.target = 0; o.mood = Mood::Sleep; return o; }
    }
    // Fire: a creature that fears it backs off (3 s away from it).
    if (c.fireScares && in.nearFire) { s.mood = Mood::Flee; s.timer = std::max(s.timer, 3.0f); s.goal = in.fire; }
    // Hurt: the passive run; one badly hurt retreats; the rest turn on whoever hit it for a while.
    if (in.hurt) {
        if (c.aggro == Aggro::Passive || (c.retreatAt > 0.0f && in.health01 < c.retreatAt)) {
            s.mood = c.aggro == Aggro::Passive ? Mood::Flee : Mood::Retreat;
            s.timer = 5.0f; s.goal = in.hurtFrom; s.windUpAt = -1.0;
        } else { s.angryUntil = in.now + c.angryTime; s.angryAt = in.hurtBy; }
    }
    // Fleeing or retreating: away from what it runs from until the timer runs out (no attacks meanwhile).
    if ((s.mood == Mood::Flee || s.mood == Mood::Retreat) && s.timer > 0.0f) {
        s.timer -= dt;
        s.windUpAt = -1.0;
        o.mood = s.mood; o.move = true; o.away = true; o.goal = s.goal; o.speed = 1.0f;
        return o;
    }
    // Something it fears in sight: away from it.
    for (const Percept& p : in.others)
        if (p.seen && ft.stance(c.factionIdx, p.faction) == Stance::Fears && p.dist < Sense_SightRange(c, 1.0f)) {
            s.mood = Mood::Flee; s.timer = 3.0f; s.goal = p.pos; s.windUpAt = -1.0;
            o.mood = s.mood; o.move = true; o.away = true; o.goal = p.pos; o.speed = 1.0f;
            return o;
        }
    // A target: the nearest it is hostile to (a predator by its range; a territorial one when it comes near and is seen;
    // defensive and territorial ones also whoever hit it, while angry), within its chase range. Never on peaceful.
    const Percept* best = nullptr;
    if (!in.peaceful && !(c.sleepsByDay && in.daylight) && c.aggro != Aggro::Passive)
        for (const Percept& p : in.others) {
            if (!p.vulnerable || p.safe) continue;
            const bool angryAt = in.now < s.angryUntil && p.id == s.angryAt;
            const bool hostile = ft.stance(c.factionIdx, p.faction) == Stance::Hostile;
            // (#1851) Taken within its range; the target it has is kept out to its chase range (it gave up at the range).
            const bool inRange = p.dist < c.aggroRange || (p.id == s.target && s.target != 0);
            bool want = false;
            if (c.aggro == Aggro::Predator) want = hostile && inRange && (p.seen || p.heard);
            else if (c.aggro == Aggro::Territorial) want = (hostile && inRange && (p.seen || p.heard)) || angryAt;
            else want = angryAt; // defensive
            if (!want || p.dist >= std::max(c.aggroRange, c.chaseRange)) continue; // (never given up within the range it takes at)
            if (!best || p.dist < best->dist) best = &p;
        }
    if (!best) {
        s.target = 0;
        s.windUpAt = -1.0; // a wind-up broken off (the target gone or safe) never lands later
        if (in.heardNoise && c.aggro != Aggro::Passive && !(c.sleepsByDay && in.daylight)) { s.mood = Mood::Alert; s.timer = 6.0f; s.goal = in.noiseAt; }
        if (s.mood == Mood::Alert && s.timer > 0.0f) {
            s.timer -= dt;
            o.mood = Mood::Alert; o.move = true; o.goal = s.goal; o.speed = 0.35f;
            return o;
        }
        s.mood = c.sleepsByDay && in.daylight ? Mood::Sleep : Mood::Idle;
        o.mood = s.mood;
        return o;
    }
    // Chase, and in reach the telegraphed attack: a wind-up, then the strike (which misses if the target stepped away).
    s.target = best->id;
    o.target = best->id;
    if (best->dist < c.reach) {
        s.mood = Mood::Attack;
        if (s.windUpAt < 0.0 && in.now - s.lastStrike > c.cooldown) { s.windUpAt = in.now; o.windUpStarted = true; }
    } else {
        s.mood = Mood::Chase;
        o.move = true; o.goal = best->pos; o.speed = 1.0f;
    }
    if (s.windUpAt >= 0.0 && in.now - s.windUpAt >= c.windUp) { s.windUpAt = -1.0; s.lastStrike = in.now; o.strike = true; }
    o.mood = s.mood;
    return o;
}


// ---- Spawning (#1223) ----
bool Spawn_Fits(const CreatureSpec& c, const SpawnSite& s) {
    if (s.kind == 2) return false;                                 // rubble is solid
    if (s.depth < c.depthMin || s.depth > c.depthMax) return false;
    if (s.light > c.darkMax + 1e-4f) return false;
    if (c.wet && !s.wet) return false;                             // (#1226) by the cave water only
    switch (c.habitat) {
        case Habitat::Cave: return !s.nearMouth && (s.kind == 0 || s.kind == 3 || s.kind == 4);
        case Habitat::Mouth: return s.nearMouth;                    // scavengers: by the entrances (decided)
        default: return false;                                     // the surface and night ones spawn as herds and hunters do; dungeons later
    }
}

int Spawn_Due(const CreatureSpec& c, SpawnChunk& ch, int alive, double now) {
    if (alive >= c.capPerChunk) return 0;
    if (!ch.filled) { ch.filled = true; ch.nextAt = now + c.respawn; return c.capPerChunk - alive; }
    if (c.nests && ch.hasNest && ch.nestDestroyed) return 0;      // the nest is gone: the chunk stays clear
    if (now < ch.nextAt) return 0;
    ch.nextAt = now + c.respawn;
    return 1;
}

// ---- Throws (#1225) ----
bool Throw_Aim(Vector3 from, Vector3 to, Vector3 up, float speed, float g, Vector3& vel) {
    const Vector3 d = Vector3Subtract(to, from);
    const float h = Vector3DotProduct(d, up);
    Vector3 flat = Vector3Subtract(d, Vector3Scale(up, h));
    const float x = Vector3Length(flat);
    flat = x > 1e-4f ? Vector3Scale(flat, 1.0f / x) : Vector3{ 0, 0, 0 };
    const float v2 = speed * speed;
    const float disc = v2 * v2 - g * (g * x * x + 2.0f * h * v2);
    if (x < 1e-4f) { vel = Vector3Scale(up, h >= 0.0f ? speed : 0.0f); return true; }
    if (disc < 0.0f || g <= 0.0f) {
        const float c = 0.70710678f;
        vel = Vector3Add(Vector3Scale(flat, speed * c), Vector3Scale(up, speed * c));
        return false;
    }
    const float ang = atanf((v2 - sqrtf(disc)) / (g * x)); // the flatter arc
    vel = Vector3Add(Vector3Scale(flat, speed * cosf(ang)), Vector3Scale(up, speed * sinf(ang)));
    return true;
}
#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
namespace {
// A small world of factions for the tests (the brain reads Factions(); these tests build specs against it).
CreatureSpec Thornback() {
    CreatureSpec c; std::string err;
    EXPECT_TRUE(Creatures_ParseLine("thornback | beast | aggro=territorial range=12 chase=30 angry=20 retreat=0.35 sight=1000 darksight=1 "
                                    "damage=8 reach=1.8 windup=0.5 cooldown=1.2 knock=6", c, &err)) << err;
    c.factionIdx = Factions().find(c.faction);
    return c;
}
Percept Player(float dist, bool safe = false) {
    Percept p; p.id = 7; p.pos = { dist, 0, 0 }; p.dist = dist; p.faction = Faction_Player(); p.seen = true; p.safe = safe; return p;
}
BrainOut Step(BrainState& s, const CreatureSpec& c, double now, std::vector<Percept> others, bool hurt = false, float health = 1.0f) {
    BrainInput in; in.now = now; in.others = std::move(others); in.hurt = hurt; in.hurtBy = 7; in.health01 = health;
    return Brain_Step(s, c, in, 1.0f / 60.0f);
}
}

TEST(Creatures, ParseAndRejects) {
    CreatureSpec c; std::string err;
    ASSERT_TRUE(Creatures_ParseLine("ooze | ooze | aggro=defensive type=heat resist=heat:0.5,cold:0.2 armour=2 depth=10-80 loot=ooze_gel:2,stone:1 sound_telegraph=ooze_gurgle cap=3", c, &err)) << err;
    EXPECT_EQ(c.aggro, Aggro::Defensive);
    EXPECT_EQ(c.type, DamageType::Heat);
    EXPECT_FLOAT_EQ(c.resist.r[(int)DamageType::Heat], 0.5f);
    EXPECT_FLOAT_EQ(c.resist.armour, 2.0f);
    EXPECT_FLOAT_EQ(c.depthMin, 10.0f); EXPECT_FLOAT_EQ(c.depthMax, 80.0f);
    ASSERT_EQ(c.loot.size(), 2u); EXPECT_EQ(c.loot[0].first, "ooze_gel"); EXPECT_EQ(c.loot[0].second, 2);
    EXPECT_EQ(c.sounds[(int)CreatureEvent::Telegraph], "ooze_gurgle");
    EXPECT_EQ(c.capPerChunk, 3);
    EXPECT_FALSE(Creatures_ParseLine("x | beast | speed=3", c, &err));     // unknown key
    EXPECT_FALSE(Creatures_ParseLine("x | beast | aggro=angry", c, &err)); // bad value
    EXPECT_FALSE(Creatures_ParseLine("x", c, &err));                       // no faction
}

TEST(Creatures, FactionTable) {
    FactionTable t; std::string err;
    ASSERT_TRUE(Factions_Parse("# test\nbeast | hostile: player, livestock\nlivestock | fears: beast\nscavenger | hostile: lurker | friendly: player\n", t, &err)) << err;
    const int beast = t.find("beast"), player = t.find("player"), stock = t.find("livestock"), lurker = t.find("lurker"), scav = t.find("scavenger");
    ASSERT_GE(lurker, 0);
    EXPECT_EQ(t.stance(beast, player), Stance::Hostile);
    EXPECT_EQ(t.stance(stock, beast), Stance::Fears);
    EXPECT_EQ(t.stance(player, beast), Stance::Neutral); // not stated: neutral
    EXPECT_EQ(t.stance(scav, lurker), Stance::Hostile);
    EXPECT_EQ(t.stance(beast, beast), Stance::Friendly);
    EXPECT_FALSE(Factions_Parse("beast | loves: player\n", t, &err));
}

TEST(Creatures, SensesLightAndNoise) {
    CreatureSpec c; c.sight = 20.0f; c.darkSight = 0.25f; c.hearing = 10.0f;
    EXPECT_FLOAT_EQ(Sense_SightRange(c, 1.0f), 20.0f);
    EXPECT_FLOAT_EQ(Sense_SightRange(c, 0.0f), 5.0f);
    EXPECT_TRUE(Sense_Sees(c, 15.0f, 1.0f, true));
    EXPECT_FALSE(Sense_Sees(c, 15.0f, 0.0f, true));  // in the dark it sees 5 m
    EXPECT_FALSE(Sense_Sees(c, 2.0f, 1.0f, false));  // behind rock
    EXPECT_TRUE(Sense_Hears(c, 9.0f, 1.0f));
    EXPECT_FALSE(Sense_Hears(c, 9.0f, 0.4f));        // a walk is quieter
    Noise_Clear();
    Noise_Emit(1, { 0, 0, 0 }, 1.0f, 10.0);
    EXPECT_EQ(Noise_Recent(10.5).size(), 1u);
    EXPECT_EQ(Noise_Recent(12.0).size(), 0u);        // kept 1.5 s
}

// The acceptance gate (#1223): a thornback on the brain behaves as the old hostile code did.
TEST(Creatures, ThornbackAsBefore) {
    if (Faction_Player() < 0 || Factions().find("beast") < 0) GTEST_SKIP() << "factions.txt not found (run from the repo root)";
    const CreatureSpec c = Thornback();
    BrainState s;
    EXPECT_EQ(Step(s, c, 0.0, { Player(13.0f) }).mood, Mood::Idle);        // outside 12 m: left alone
    BrainOut o = Step(s, c, 0.1, { Player(11.0f) });
    EXPECT_EQ(o.mood, Mood::Chase); EXPECT_TRUE(o.move);                    // inside: it charges
    EXPECT_EQ(Step(s, c, 0.2, { Player(11.0f, true) }).mood, Mood::Idle); // the player by a fire is safe
    // In reach: the wind-up, the strike 0.5 s later, then the 1.2 s cooldown.
    o = Step(s, c, 1.0, { Player(1.0f) });
    EXPECT_EQ(o.mood, Mood::Attack); EXPECT_TRUE(o.windUpStarted); EXPECT_FALSE(o.strike);
    EXPECT_FALSE(Step(s, c, 1.4, { Player(1.0f) }).strike);
    EXPECT_TRUE(Step(s, c, 1.5, { Player(1.0f) }).strike);
    EXPECT_FALSE(Step(s, c, 2.0, { Player(1.0f) }).windUpStarted);       // cooling down
    EXPECT_TRUE(Step(s, c, 2.8, { Player(1.0f) }).windUpStarted);
    // Hit from 25 m: angry for 20 s, it chases out to 30 m, not past.
    BrainState h;
    EXPECT_EQ(Step(h, c, 10.0, { Player(25.0f) }, true).mood, Mood::Chase);
    EXPECT_EQ(Step(h, c, 20.0, { Player(29.0f) }).mood, Mood::Chase);
    EXPECT_EQ(Step(h, c, 21.0, { Player(31.0f) }).mood, Mood::Idle);
    EXPECT_EQ(Step(h, c, 31.0, { Player(20.0f) }).mood, Mood::Idle);        // anger over
    // Badly hurt (under 35%): it runs for 5 s.
    BrainState r;
    EXPECT_EQ(Step(r, c, 0.0, { Player(3.0f) }, true, 0.3f).mood, Mood::Retreat);
}

// (#1851) A target taken within the range is chased out to the chase range: it gave up at the range.
TEST(Creatures, ChaseKeptToChaseRange) {
    if (Faction_Player() < 0 || Factions().find("beast") < 0) GTEST_SKIP() << "factions.txt not found (run from the repo root)";
    for (const char* aggro : { "predator", "territorial" }) {
        CreatureSpec c; std::string err;
        ASSERT_TRUE(Creatures_ParseLine(std::string("stalker | beast | aggro=") + aggro + " range=12 chase=30 sight=1000 darksight=1 damage=8 reach=1.8", c, &err)) << err;
        c.factionIdx = Factions().find(c.faction);
        BrainState s;
        EXPECT_EQ(Step(s, c, 0.0, { Player(13.0f) }).mood, Mood::Idle) << aggro;  // not taken past its range
        EXPECT_EQ(Step(s, c, 0.1, { Player(11.0f) }).mood, Mood::Chase) << aggro; // taken
        EXPECT_EQ(Step(s, c, 0.2, { Player(13.0f) }).mood, Mood::Chase) << aggro; // kept past the range
        EXPECT_EQ(Step(s, c, 0.3, { Player(29.0f) }).mood, Mood::Chase) << aggro;
        EXPECT_EQ(Step(s, c, 0.4, { Player(31.0f) }).mood, Mood::Idle) << aggro;  // given up past the chase range
        EXPECT_EQ(Step(s, c, 0.5, { Player(20.0f) }).mood, Mood::Idle) << aggro;  // and not taken again from there
    }
}

TEST(Creatures, NightHunterFireAndDay) {
    if (Faction_Player() < 0) GTEST_SKIP() << "factions.txt not found";
    CreatureSpec c; std::string err;
    ASSERT_TRUE(Creatures_ParseLine("dusk_stalker | beast | aggro=predator range=70 sight=1000 darksight=1 sleeps=day fire=scares reach=1.7", c, &err)) << err;
    c.factionIdx = Factions().find(c.faction);
    BrainState s;
    BrainInput in; in.now = 0.0; in.daylight = true; in.others = { Player(40.0f) };
    EXPECT_NE(Brain_Step(s, c, in, 0.1f).mood, Mood::Chase); // daylight: no hunt
    in.daylight = false;
    EXPECT_EQ(Brain_Step(s, c, in, 0.1f).mood, Mood::Chase); // dark, within 70 m
    in.nearFire = true; in.fire = { 5, 0, 0 };
    const BrainOut o = Brain_Step(s, c, in, 0.1f);
    EXPECT_EQ(o.mood, Mood::Flee); EXPECT_TRUE(o.away);      // a fire within 14 m: it backs off
}

TEST(Creatures, InfightingAndAlert) {
    if (Factions().find("beast") < 0 || Factions().find("livestock") < 0) GTEST_SKIP() << "factions.txt not found";
    CreatureSpec c = Thornback();
    BrainState s;
    Percept prey; prey.id = 99; prey.dist = 6.0f; prey.pos = { 6, 0, 0 }; prey.faction = Factions().find("livestock"); prey.seen = true;
    const bool beastHuntsStock = Factions().stance(c.factionIdx, prey.faction) == Stance::Hostile;
    const BrainOut o = Step(s, c, 0.0, { prey });
    EXPECT_EQ(o.mood == Mood::Chase, beastHuntsStock); // as factions.txt says
    // A noise with nothing in sight: it goes to look.
    BrainState a;
    BrainInput in; in.now = 0.0; in.heardNoise = true; in.noiseAt = { 3, 0, 4 };
    const BrainOut al = Brain_Step(a, c, in, 0.1f);
    EXPECT_EQ(al.mood, Mood::Alert); EXPECT_TRUE(al.move);
    // Peaceful: never hostile.
    BrainState p;
    BrainInput pin; pin.now = 0.0; pin.peaceful = true; pin.others = { Player(2.0f) };
    EXPECT_EQ(Brain_Step(p, c, pin, 0.1f).mood, Mood::Idle);
}

TEST(Creatures, SpawnFitsAndDue) {
    CreatureSpec c; std::string err;
    ASSERT_TRUE(Creatures_ParseLine("lurker | lurker | habitat=cave depth=8-200 dark=0.1 cap=2 respawn=600 nests=1", c, &err)) << err;
    SpawnSite s; s.kind = 0; s.depth = 20.0f; s.light = 0.0f;
    EXPECT_TRUE(Spawn_Fits(c, s));
    s.depth = 4.0f; EXPECT_FALSE(Spawn_Fits(c, s));       // too shallow
    s.depth = 20.0f; s.light = 0.5f; EXPECT_FALSE(Spawn_Fits(c, s)); // too light
    s.light = 0.0f; s.nearMouth = true; EXPECT_FALSE(Spawn_Fits(c, s)); // a cave dweller keeps off the entrances
    s.nearMouth = false; s.kind = 2; EXPECT_FALSE(Spawn_Fits(c, s));   // rubble
    CreatureSpec m; ASSERT_TRUE(Creatures_ParseLine("scav | scavenger | habitat=mouth depth=0-30", m, &err)) << err;
    SpawnSite mouth; mouth.kind = 1; mouth.depth = 3.0f; mouth.nearMouth = true;
    EXPECT_TRUE(Spawn_Fits(m, mouth));
    mouth.nearMouth = false; EXPECT_FALSE(Spawn_Fits(m, mouth));
    // A chunk: filled to the cap at once, then one per 600 s while below it, none once its nest is destroyed.
    SpawnChunk ch;
    EXPECT_EQ(Spawn_Due(c, ch, 0, 0.0), 2);
    EXPECT_EQ(Spawn_Due(c, ch, 2, 10.0), 0);   // full
    EXPECT_EQ(Spawn_Due(c, ch, 0, 300.0), 0);  // cleared, waiting
    EXPECT_EQ(Spawn_Due(c, ch, 0, 600.0), 1);
    EXPECT_EQ(Spawn_Due(c, ch, 1, 700.0), 0);
    ch.hasNest = true; ch.nestDestroyed = true;
    EXPECT_EQ(Spawn_Due(c, ch, 0, 5000.0), 0);
}

// #1224: a blind hunter (sight 0) goes by sound alone: a heard player is chased, a silent one (sneaking) is not.
TEST(Creatures, BlindHunterBySound) {
    if (Faction_Player() < 0 || Factions().find("lurker") < 0) GTEST_SKIP() << "factions.txt not found";
    CreatureSpec c; std::string err;
    ASSERT_TRUE(Creatures_ParseLine("tunnel_lurker | lurker | aggro=predator range=24 sight=0 darksight=0 hearing=22 ceiling=1 burst=1 drop=1.5", c, &err)) << err;
    c.factionIdx = Factions().find(c.faction);
    EXPECT_TRUE(c.ceiling); EXPECT_TRUE(c.burst); EXPECT_FLOAT_EQ(c.dropScale, 1.5f);
    EXPECT_FALSE(Sense_Sees(c, 3.0f, 1.0f, true));              // blind, even in full light
    EXPECT_TRUE(Sense_Hears(c, 20.0f, NoiseTune::kSprint));      // a sprint at 20 m
    EXPECT_FALSE(Sense_Hears(c, 12.0f, NoiseTune::kWalk));       // a walk only within ~9 m
    Percept p; p.id = 7; p.dist = 10.0f; p.pos = { 10, 0, 0 }; p.faction = Faction_Player();
    BrainState s; BrainInput in; in.now = 0.0; in.others = { p };
    EXPECT_EQ(Brain_Step(s, c, in, 0.1f).mood, Mood::Idle);      // silent: nothing
    in.others[0].heard = true;
    EXPECT_EQ(Brain_Step(s, c, in, 0.1f).mood, Mood::Chase);     // heard: it comes
}

// #1225: a throw's arc lands where it was aimed (12 m off, 1 m up), and out of reach it says so.
TEST(Creatures, ThrowAimLands) {
    const Vector3 up{ 0, 1, 0 }, from{ 0, 1.5f, 0 }, to{ 12.0f, 2.5f, 0 };
    Vector3 v;
    ASSERT_TRUE(Throw_Aim(from, to, up, 16.0f, 9.81f, v));
    Vector3 p = from, vel = v;
    float best = 1e9f;
    for (int i = 0; i < 600; ++i) {
        vel = Vector3Add(vel, Vector3Scale(up, -9.81f / 600.0f * 2.0f));
        p = Vector3Add(p, Vector3Scale(vel, 2.0f / 600.0f));
        best = std::min(best, Vector3Distance(p, to));
    }
    EXPECT_LT(best, 0.3f);
    EXPECT_FALSE(Throw_Aim(from, { 200.0f, 0, 0 }, up, 16.0f, 9.81f, v)); // 200 m: out of reach at 16 m/s
    CreatureSpec c; std::string err;
    ASSERT_TRUE(Creatures_ParseLine("scav | scavenger | group=3 thrower=1 thief=1 torch=1 throw=6 throwmin=5 throwrange=14 throwwindup=0.6 throwcooldown=2.8", c, &err)) << err;
    EXPECT_EQ(c.group, 3); EXPECT_TRUE(c.thrower && c.thief && c.torch); EXPECT_FLOAT_EQ(c.throwRange, 14.0f);
}

// #1227: a light sleeper sleeps through a player near in the dark, wakes at a light (or a hit), and stays awake a minute.
TEST(Creatures, LightSleeperAndShieldKeys) {
    if (Faction_Player() < 0 || Factions().find("undead") < 0) GTEST_SKIP() << "factions.txt not found";
    CreatureSpec c; std::string err;
    ASSERT_TRUE(Creatures_ParseLine("lost_miner | undead | aggro=territorial range=12 sight=1000 darksight=1 wake=light wakerange=10 shield=0.85 shieldarc=60 rare=map_fragment:1:0.25", c, &err)) << err;
    c.factionIdx = Factions().find(c.faction);
    EXPECT_TRUE(c.wakeByLight); EXPECT_FLOAT_EQ(c.frontShield, 0.85f);
    ASSERT_EQ(c.rare.size(), 1u); EXPECT_EQ(c.rare[0].item, "map_fragment"); EXPECT_FLOAT_EQ(c.rare[0].chance, 0.25f);
    BrainState s; BrainInput in; in.now = 0.0; in.others = { Player(4.0f) };
    EXPECT_EQ(Brain_Step(s, c, in, 0.1f).mood, Mood::Sleep);  // a player 4 m off in the dark: asleep
    in.lightNear = true; in.now = 1.0;
    EXPECT_EQ(Brain_Step(s, c, in, 0.1f).mood, Mood::Chase);  // a light: awake and after them
    in.lightNear = false; in.now = 30.0;
    EXPECT_EQ(Brain_Step(s, c, in, 0.1f).mood, Mood::Chase);  // still awake 29 s on
    in.now = 62.0; in.others.clear();
    EXPECT_EQ(Brain_Step(s, c, in, 0.1f).mood, Mood::Sleep);  // a minute after the light: asleep again
}

// #1226: the ooze keys, and a wet-only species keeps to wet sites.
TEST(Creatures, OozeKeysAndWetSites) {
    CreatureSpec c; std::string err;
    ASSERT_TRUE(Creatures_ParseLine("ooze | ooze | split=1 corrode=15 wet=1 habitat=cave depth=10-400 dark=0.2", c, &err)) << err;
    EXPECT_TRUE(c.split); EXPECT_TRUE(c.wet); EXPECT_EQ(c.corrode, 15);
    SpawnSite s; s.kind = 0; s.depth = 30.0f; s.light = 0.0f;
    EXPECT_FALSE(Spawn_Fits(c, s)); // dry
    s.wet = true;
    EXPECT_TRUE(Spawn_Fits(c, s));
}

// #1228: the sentinel keys; a weakness parses as a negative resistance.
TEST(Creatures, SentinelKeys) {
    CreatureSpec c; std::string err;
    ASSERT_TRUE(Creatures_ParseLine("ruin_sentinel | construct | ranged=1 dart=1 dartspeed=28 immune=bleed,poison resist=blunt:-0.5,slash:0.4", c, &err)) << err;
    EXPECT_TRUE(c.ranged); EXPECT_TRUE(c.dart); EXPECT_FLOAT_EQ(c.dartSpeed, 28.0f);
    EXPECT_TRUE(c.immuneBleed); EXPECT_TRUE(c.immunePoison);
    EXPECT_FLOAT_EQ(c.resist.r[(int)DamageType::Blunt], -0.5f);
    EXPECT_FLOAT_EQ(Damage_Apply(10.0f, DamageType::Blunt, c.resist), 15.0f);
    EXPECT_FALSE(Creatures_ParseLine("x | construct | immune=fire", c, &err));
}
#endif
