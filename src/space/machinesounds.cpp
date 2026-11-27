#include "machinesounds.h"

#include "audio.h"
#include "factory.h"
#include "player.h"
#include "profiler.h"
#include "scene.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <unordered_map>

const char* MachineSound_Hum(PartType t) {
    switch (t) {
        case PartType::Campfire: case PartType::Kiln: case PartType::CokeOven: case PartType::BlastFurnace: return "hum_fire";
        case PartType::Boiler: case PartType::SteamMacerator: case PartType::ForgeHammer: case PartType::AlloySmelter: return "hum_steam";
        case PartType::Generator: case PartType::Furnace: case PartType::Macerator: case PartType::Battery:
        case PartType::OreWasher: case PartType::Centrifuge: case PartType::WaterPump: return "hum_lv";
        case PartType::MvGenerator: case PartType::Electrolyzer: case PartType::ChemicalReactor: case PartType::Miner: return "hum_mv";
        case PartType::HvGenerator: case PartType::ArcFurnace: case PartType::Refinery: case PartType::PrecisionAssembler: return "hum_hv";
        case PartType::FusionReactor: case PartType::MassDriver: return "hum_ev";
        // (#1269, #1270) mechanical power: transmission clatter, and each rotation machine and source its own voice
        case PartType::Shaft: case PartType::Cogwheel: case PartType::LargeCogwheel: case PartType::Gearbox: case PartType::HandCrank: return "mech_gears";
        case PartType::Millstone: return "mech_millstone";
        case PartType::Saw: return "mech_saw";
        case PartType::Press: case PartType::Lathe: case PartType::Wiremill: case PartType::TripHammer: return "mech_press";
        case PartType::Mixer: case PartType::RotarySieve: case PartType::RotaryPump: return "mech_millstone"; // (#1734) GUESS: a slow grinding turn
        case PartType::WaterWheel: return "mech_waterwheel";
        case PartType::Windmill: return "mech_windmill";
        default: return nullptr; // chests, hand pumps, rocket parts, connectors, gauges
    }
}

bool MachineSound_IsTransmission(PartType t) {
    return t == PartType::Shaft || t == PartType::Cogwheel || t == PartType::LargeCogwheel || t == PartType::Gearbox || t == PartType::HandCrank;
}

float MachineSound_RpmPitch(float rpm) {
    return std::clamp(fabsf(rpm) / 64.0f, 0.5f, 2.0f); // the loops are made at 64 RPM: slower turning plays them slower
}

bool MachineSound_IdleHums(PartType t) {
    const char* h = MachineSound_Hum(t);
    return h && (!strcmp(h, "hum_lv") || !strcmp(h, "hum_mv") || !strcmp(h, "hum_hv") || !strcmp(h, "hum_ev")) && t != PartType::Generator &&
           t != PartType::MvGenerator && t != PartType::HvGenerator; // a generator idle is off; a powered consumer hums
}

MachineCue MachineSound_Cue(int was, int now) {
    if (was == now || was < 0) return MachineCue::None; // (was < 0: just seen)
    if (now == 2) return MachineCue::Stall;            // stopped for a reason: no power, output full, no fuel, no water
    if (was == 1 && now == 0) return MachineCue::Done; // ran out of work: the batch is done
    return MachineCue::None;
}

float MachineSound_BedVolume(int count) {
    if (count <= 0) return 0.0f;
    return std::min(1.0f, 0.18f * log2f(1.0f + (float)count)); // 1: 0.18, 7: 0.54, 31: 0.9, 63+: 1
}

void MachineSound_Pick(const std::vector<float>& distances, const std::vector<bool>& audible, int max, std::vector<int>& out) {
    out.clear();
    std::vector<int> idx(distances.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return distances[a] < distances[b]; });
    for (int i : idx) {
        if ((int)out.size() >= max) break;
        if (audible[i]) out.push_back(i);
    }
}

namespace {
    constexpr float kRange = 64.0f;      // machines within this of the listener are considered
    constexpr float kCueRange = 40.0f;   // state cues within this
    struct Spot { int part; PartType type; Vector3 centre; }; // centre in plot metres
    struct PlotCache { uint32_t version = ~0u; std::vector<Spot> spots; };
    std::unordered_map<uint64_t, PlotCache> cache;
    struct Live { uint64_t key; uint64_t body; uint32_t plot; int part; PartType type; DVec3 pos; float dist; int state; float pitch = 1.0f; };
    std::vector<Live> live;
    std::unordered_map<uint64_t, int> lastState;
    std::unordered_map<uint64_t, double> nextClick;
    struct Slot { int emitter = -1, sound = -1; uint64_t key = 0; bool used = false; };
    Slot slots[kMachineIndividual];
    int bedEmitter = -1;
    double nextScan = 0.0;

    uint64_t KeyOf(uint64_t body, uint32_t plot, int part) { return (body * 1315423911ull) ^ ((uint64_t)plot << 20) ^ (uint64_t)part; }

    void PlayCue(const char* name, const DVec3& at, float volume) {
        const int snd = Audio_Sound(name);
        if (snd < 0) return;
        AudioAttenuation att;
        att.refDistance = 2.0f;
        att.maxDistance = kCueRange;
        att.lowpassDistance = 30.0f;
        Audio_PlaySoundAt(snd, at, AudioBus::Machines, AudioPriority::Normal, volume, 1.0f, att);
    }

    void Scan(Planetoid& p) {
        live.clear();
        if (!p.plots) return;
        const Vector3 lis = p.toLocalD(Audio_GetListener().pos).toFloat();
        std::vector<std::pair<Int3, CellShape>> cells;
        for (auto& plotPtr : p.plots->plots) {
            Plot& plot = *plotPtr;
            if (plot.cells.empty() || Vector3Distance(plot.toBody({ 0, 0, 0 }), lis) > kRange + 200.0f) continue;
            PlotCache& pc = cache[(p.id * 1315423911ull) ^ plot.id];
            if (pc.version != plot.version) {
                pc.version = plot.version;
                pc.spots.clear();
                for (size_t i = 0; i < plot.parts.size(); ++i) {
                    const Part& part = plot.parts[i];
                    if (part.type == PartType::Count || !MachineSound_Hum(part.type)) continue;
                    Plot::Footprint(part, cells);
                    Vector3 c{ 0, 0, 0 };
                    for (const auto& cell : cells) c = Vector3Add(c, { cell.first.x + 0.5f, cell.first.y + 0.5f, cell.first.z + 0.5f });
                    if (!cells.empty()) c = Vector3Scale(c, kCellSize / (float)cells.size());
                    pc.spots.push_back({ (int)i, part.type, c });
                }
            }
            if (pc.spots.empty()) continue;
            const std::vector<int> stalled = Factory::StalledParts(p.id, plot.id);
            // (#1269) Rotation parts turn at their network's speed: transmission parts of a plot sound as one clatter at
            // their middle (a 20-shaft line is one voice, not twenty); machines and sources each by their own RPM.
            static std::vector<std::pair<int, float>> rpms;
            rpms.clear();
            bool anyRotation = false;
            for (const Spot& s : pc.spots) anyRotation = anyRotation || IsRotationPart(s.type);
            if (anyRotation) Factory::RotationSpeeds(p.id, plot.id, rpms);
            auto rpmOf = [&](int part) { for (const auto& kv : rpms) if (kv.first == part) return kv.second; return 0.0f; };
            Vector3 gearSum{ 0, 0, 0 };
            int gearN = 0;
            float gearRpm = 0.0f;
            for (const Spot& s : pc.spots) {
                const Vector3 at = plot.toBody(s.centre);
                const float d = Vector3Distance(at, lis);
                if (d > kRange) continue;
                if (IsRotationPart(s.type)) {
                    const float rpm = fabsf(rpmOf(s.part));
                    if (MachineSound_IsTransmission(s.type)) { gearSum = Vector3Add(gearSum, at); ++gearN; gearRpm += rpm; continue; }
                    live.push_back({ KeyOf(p.id, plot.id, s.part), p.id, plot.id, s.part, s.type, p.toWorldD(DVec3(at)), d, rpm > 1.0f ? 1 : 0, MachineSound_RpmPitch(rpm) });
                    continue;
                }
                const bool working = Factory::Working(p.id, plot.id, s.part);
                const bool stall = !working && std::find(stalled.begin(), stalled.end(), s.part) != stalled.end();
                live.push_back({ KeyOf(p.id, plot.id, s.part), p.id, plot.id, s.part, s.type, p.toWorldD(DVec3(at)), d, working ? 1 : stall ? 2 : 0 });
            }
            if (gearN > 0) {
                const Vector3 at = Vector3Scale(gearSum, 1.0f / gearN);
                const float rpm = gearRpm / gearN;
                live.push_back({ KeyOf(p.id, plot.id, -1), p.id, plot.id, -1, PartType::Shaft, p.toWorldD(DVec3(at)), Vector3Distance(at, lis),
                                 rpm > 1.0f ? 1 : 0, MachineSound_RpmPitch(rpm) });
            }
        }
    }
}

void MachineSounds_Reset() {
    for (Slot& s : slots) { Audio_EmitterDestroy(s.emitter); s = Slot{}; }
    Audio_EmitterDestroy(bedEmitter);
    bedEmitter = -1;
    cache.clear(); live.clear(); lastState.clear(); nextClick.clear();
}

void MachineSounds_Update(Scene& world, float dt) {
    (void)world; (void)dt;
    if (!Audio_V2()) return;
    const double now = GetTime();
    if (now < nextScan) return;
    nextScan = now + 0.25;
    Planetoid* p = player.nearest;
    if (!p) return;
    PROFILE_SCOPE("audio:machines");
    Scan(*p);

    // State cues (a machine first seen has no "before": no cue).
    std::unordered_map<uint64_t, int> seen;
    for (const Live& m : live) {
        auto it = lastState.find(m.key);
        const MachineCue cue = MachineSound_Cue(it == lastState.end() ? -1 : it->second, m.state);
        if (m.dist < kCueRange && !IsRotationPart(m.type)) { // (a shaft stopping is no finished batch)
            if (cue == MachineCue::Stall) { PlayCue("machine_stall", m.pos, 0.6f); nextClick[m.key] = now + 3.0; }
            else if (cue == MachineCue::Done) PlayCue("machine_done", m.pos, 0.45f);
            else if (m.state == 2 && m.dist < 24.0f && now >= nextClick[m.key]) { PlayCue("machine_click", m.pos, 0.5f); nextClick[m.key] = now + 3.0; }
        }
        seen[m.key] = m.state;
    }
    lastState.swap(seen);

    // The nearest audible machines get their own voices; the other working ones make the factory bed.
    std::vector<float> dist(live.size());
    std::vector<bool> audible(live.size());
    for (size_t i = 0; i < live.size(); ++i) {
        dist[i] = live[i].dist;
        audible[i] = live[i].state == 1 || (live[i].state == 0 && MachineSound_IdleHums(live[i].type));
    }
    std::vector<int> picked;
    MachineSound_Pick(dist, audible, kMachineIndividual, picked);
    for (Slot& s : slots) s.used = false;
    std::vector<int> unplaced;
    for (int i : picked) { // keep a machine on the slot it already had
        bool kept = false;
        for (Slot& s : slots)
            if (!s.used && s.emitter >= 0 && s.key == live[i].key) { s.used = true; kept = true; break; }
        if (!kept) unplaced.push_back(i);
    }
    for (int i : unplaced) {
        Slot* free = nullptr;
        for (Slot& s : slots) if (!s.used) { free = &s; break; }
        if (!free) break;
        const int snd = Audio_Sound(MachineSound_Hum(live[i].type));
        if (snd < 0) continue;
        if (free->sound != snd) { // another kind of hum: a new emitter
            Audio_EmitterDestroy(free->emitter);
            AudioAttenuation att;
            att.refDistance = 2.0f;
            att.maxDistance = kRange;
            att.lowpassDistance = 40.0f;
            free->emitter = Audio_EmitterCreate(snd, AudioBus::Machines, AudioPriority::Low, att, true);
            free->sound = snd;
        }
        free->key = live[i].key;
        free->used = free->emitter >= 0;
    }
    for (Slot& s : slots) {
        if (s.emitter < 0) continue;
        if (!s.used) { Audio_EmitterSetPlaying(s.emitter, false); continue; }
        const Live* m = nullptr;
        for (const Live& l : live) if (l.key == s.key) { m = &l; break; }
        if (!m) { Audio_EmitterSetPlaying(s.emitter, false); continue; }
        Audio_EmitterSetPosition(s.emitter, m->pos);
        Audio_EmitterSetVolume(s.emitter, m->state == 1 ? 0.8f : 0.2f); // idle: a quiet hum
        Audio_EmitterSetPitch(s.emitter, m->state == 1 ? m->pitch : 0.97f); // (rotation: by its RPM)
        Audio_EmitterSetPlaying(s.emitter, true);
    }
    // The bed: every other working machine in range, at their middle.
    int count = 0;
    DVec3 sum{};
    for (size_t i = 0; i < live.size(); ++i) {
        if (live[i].state != 1 || std::find(picked.begin(), picked.end(), (int)i) != picked.end()) continue;
        ++count;
        sum += live[i].pos;
    }
    if (bedEmitter < 0 && count > 0) {
        const int snd = Audio_Sound("factory_bed");
        AudioAttenuation att;
        att.refDistance = 12.0f;
        att.maxDistance = 140.0f;
        att.lowpassDistance = 90.0f;
        if (snd >= 0) bedEmitter = Audio_EmitterCreate(snd, AudioBus::Machines, AudioPriority::Low, att, true);
    }
    if (bedEmitter >= 0) {
        if (count > 0) Audio_EmitterSetPosition(bedEmitter, sum * (1.0 / count));
        Audio_EmitterSetVolume(bedEmitter, MachineSound_BedVolume(count));
        Audio_EmitterSetPlaying(bedEmitter, count > 0);
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1243): hums by type, cues by state change, the bed's growth, the nearest-first pick ----
#include "unit_helpers.h"

TEST(MachineSounds, HumsByTypeAndTier) {
    EXPECT_STREQ(MachineSound_Hum(PartType::Campfire), "hum_fire");
    EXPECT_STREQ(MachineSound_Hum(PartType::Boiler), "hum_steam");
    EXPECT_STREQ(MachineSound_Hum(PartType::Macerator), "hum_lv");
    EXPECT_STREQ(MachineSound_Hum(PartType::Electrolyzer), "hum_mv");
    EXPECT_STREQ(MachineSound_Hum(PartType::ArcFurnace), "hum_hv");
    EXPECT_STREQ(MachineSound_Hum(PartType::FusionReactor), "hum_ev");
    EXPECT_EQ(MachineSound_Hum(PartType::Chest), nullptr);
    // (#1269) mechanical power
    EXPECT_STREQ(MachineSound_Hum(PartType::Shaft), "mech_gears");
    EXPECT_STREQ(MachineSound_Hum(PartType::Gearbox), "mech_gears");
    EXPECT_STREQ(MachineSound_Hum(PartType::WaterWheel), "mech_waterwheel");
    EXPECT_STREQ(MachineSound_Hum(PartType::Millstone), "mech_millstone");
    EXPECT_EQ(MachineSound_Hum(PartType::StressGauge), nullptr);
    EXPECT_TRUE(MachineSound_IsTransmission(PartType::Cogwheel));
    EXPECT_FALSE(MachineSound_IsTransmission(PartType::Windmill));
    EXPECT_FALSE(MachineSound_IdleHums(PartType::Shaft));
    EXPECT_FLOAT_EQ(MachineSound_RpmPitch(64.0f), 1.0f);
    EXPECT_FLOAT_EQ(MachineSound_RpmPitch(-128.0f), 2.0f);
    EXPECT_FLOAT_EQ(MachineSound_RpmPitch(8.0f), 0.5f);
    EXPECT_FLOAT_EQ(MachineSound_RpmPitch(1000.0f), 2.0f);
    EXPECT_TRUE(MachineSound_IdleHums(PartType::Macerator));
    EXPECT_FALSE(MachineSound_IdleHums(PartType::Generator));
    EXPECT_FALSE(MachineSound_IdleHums(PartType::Campfire));
}

TEST(MachineSounds, CuesOnStateChanges) {
    EXPECT_EQ(MachineSound_Cue(-1, 2), MachineCue::None);  // just seen
    EXPECT_EQ(MachineSound_Cue(1, 2), MachineCue::Stall);
    EXPECT_EQ(MachineSound_Cue(0, 2), MachineCue::Stall);
    EXPECT_EQ(MachineSound_Cue(1, 0), MachineCue::Done);
    EXPECT_EQ(MachineSound_Cue(2, 0), MachineCue::None);   // a stall cleared by taking the outputs: no chime
    EXPECT_EQ(MachineSound_Cue(0, 1), MachineCue::None);
    EXPECT_EQ(MachineSound_Cue(1, 1), MachineCue::None);
}

TEST(MachineSounds, TheBedStaysTameAndThePickIsNearestFirst) {
    EXPECT_FLOAT_EQ(MachineSound_BedVolume(0), 0.0f);
    EXPECT_LT(MachineSound_BedVolume(1), MachineSound_BedVolume(10));
    EXPECT_LT(MachineSound_BedVolume(100) - MachineSound_BedVolume(50), 0.2f);
    EXPECT_LE(MachineSound_BedVolume(100000), 1.0f);
    const std::vector<float> d = { 30, 5, 12, 2, 50, 8 };
    const std::vector<bool> a = { true, true, false, true, true, true };
    std::vector<int> out;
    MachineSound_Pick(d, a, 3, out);
    ASSERT_EQ(out.size(), 3u);
    EXPECT_EQ(out[0], 3);
    EXPECT_EQ(out[1], 1);
    EXPECT_EQ(out[2], 5); // (index 2 is nearer but silent)
}
#endif
