// --test-rocket (milestone 7 gate): build, launch, orbit and land.
// A rocket is built from parts on a plot near the start (an engine on four landing legs, four fuel
// tanks, a guidance computer, a cockpit), fuelled from the inventory and boarded; the autopilot flies
// a gravity turn to a 45 km apoapsis, circularises (orbit: periapsis above the 35 km air), burns
// retrograde to come down, and lands with a suicide burn. Run with --time-scale (the orbit takes
// ~25 minutes of world time).

#include "rocket_test.h"
#include "rocket.h"
#include "buildtool.h"
#include "inventory.h"
#include "player.h"
#include "survival.h"
#include "plots.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include "station.h"
#include "gamesave.h"

#include <raymath.h>
#include <algorithm>
#include <cstdio>

namespace {
    int phase = 0;
    double phaseStart = 0.0, testStart = 0.0;
    Vector3 padDir{};
    bool orbited = false, savedInOrbit = false, mapShot = false, warped = false, warpDone = false;
    double warpStart = 0.0;
    RocketTelemetry warpBefore;
    double maxAlt = 0.0;
    float airLoaded = 0.0f;   // minutes of air loaded into the cockpit
    void Next() { ++phase; phaseStart = GetTime(); }
    double T() { return GetTime() - phaseStart; }
}

int RocketTest_Update(Scene& world) {
    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
    if (!pf) { printf("bench: --test-rocket needs the planet\n"); return -1; }
    const RocketTelemetry t = Rocket_Telemetry(world);
    if (phase >= 2) maxAlt = std::max(maxAlt, t.altitude);
    switch (phase) {
    case 0: {
        testStart = GetTime();
        const Vector3 home = Vector3Normalize(p->toLocal((player.position - world.renderOrigin).toFloat()));
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 }));
        padDir = Vector3Normalize(Vector3Add(home, Vector3Scale(t1, 14.0f / pf->radius)));
        if (!p->plots) p->plots = std::make_unique<PlotSet>();
        Plot* plot = p->plots->plotFor(Vector3Scale(padDir, pf->radius + pf->elevation(padDir)), 0.0f, true);
        bool ok = BuildTool_Place(*p, *plot, { PartType::RocketEngine, { 0, 0, 0 }, 0 }, false);
        const Int3 legs[4] = { { -1, 0, -1 }, { 2, 0, -1 }, { -1, 0, 2 }, { 2, 0, 2 } };
        for (const Int3& l : legs) ok &= BuildTool_Place(*p, *plot, { PartType::LandingLeg, l, 0 }, false);
        for (int k = 0; k < 4; ++k) ok &= BuildTool_Place(*p, *plot, { PartType::FuelTank, { 0, 2 + 2 * k, 0 }, 0 }, false);
        ok &= BuildTool_Place(*p, *plot, { PartType::Guidance, { 0, 10, 0 }, 0 }, false);
        ok &= BuildTool_Place(*p, *plot, { PartType::Cockpit, { 0, 11, 0 }, 0 }, false);
        inventory.add("rocket_fuel", 60);
        const int loaded = Rocket_LoadFuel(*plot, 60);
        inventory.add("oxygen_canister", 1);
        airLoaded = Rocket_LoadAir(*plot, 1) * 10.0f; // 10 min for the cockpit's life support
        const RocketStats s = Rocket_Stats(*plot, p);
        printf("bench: rocket: built %d parts (%s); fuel %d items -> %.0f kg; dry %.0f kg, thrust %.0f kN, TWR %.2f, delta-v %.0f m/s (orbit needs ~%.0f): %s\n",
               plot->partCount(), ok ? "ok" : "FAILED", loaded, s.fuel, s.dryMass, s.thrust / 1000.0f, s.twr, s.deltaV, s.deltaVNeeded, s.ready ? "READY" : s.missing.front().c_str());
        if (!ok || !s.ready || !Rocket_Board(world, *p, *plot)) return -1;
        Next();
        break;
    }
    case 1:
        if (T() > 1.0) {
            TakeScreenshot("bench_rocket_pad.png");
            Rocket_SetAutopilot(RocketAutopilot::Ascent, 45000.0);
            printf("bench: rocket: launch (autopilot: gravity turn to a 45 km apoapsis)\n");
            Next();
        }
        break;
    case 2: // ascent and circularisation
        if (t.state == RocketState::Flying && t.periapsis > 35000.0 && t.apoapsis > 35000.0) {
            orbited = true;
            TakeScreenshot("bench_rocket_orbit.png");
            printf("bench: rocket: IN ORBIT after %.0f s: apoapsis %.1f km, periapsis %.1f km, %.0f m/s, fuel %.0f kg left\n",
                   GetTime() - testStart, t.apoapsis / 1000.0, t.periapsis / 1000.0, t.orbitalSpeed, t.fuel);
            Next();
        } else if (T() > 240.0) { printf("bench: rocket: no orbit (apoapsis %.1f km, periapsis %.1f km, fuel %.0f kg)\n", t.apoapsis / 1000.0, t.periapsis / 1000.0, t.fuel); return -1; }
        break;
    case 3: // coast a while, then come down
        if (T() > 2.0 && !savedInOrbit) {
            // A save and a load in flight: the orbit must survive.
            savedInOrbit = true;
            const RocketTelemetry before = Rocket_Telemetry(world);
            const bool ok = SaveWorld(world) && LoadWorldSave(world);
            const RocketTelemetry after = Rocket_Telemetry(world);
            const bool same = ok && after.state == RocketState::Flying && fabs(after.apoapsis - before.apoapsis) < 500.0 && fabs(after.periapsis - before.periapsis) < 500.0;
            printf("bench: rocket: saved and loaded in orbit: %s (apoapsis %.1f -> %.1f km, periapsis %.1f -> %.1f km)\n", same ? "orbit kept" : "LOST",
                   before.apoapsis / 1000.0, after.apoapsis / 1000.0, before.periapsis / 1000.0, after.periapsis / 1000.0);
            if (!same) return -1;
        }
        // The orbit map, then time warp: 100x for a while (coasting), the orbit must hold.
        if (T() > 3.0 && !mapShot) { mapShot = true; Rocket_SetMapView(true); }
        if (T() > 3.5 && !warped) {
            warped = true;
            TakeScreenshot("bench_rocket_map.png");
            Rocket_SetMapView(false);
            warpBefore = Rocket_Telemetry(world);
            Rocket_SetWarp(6);
            warpStart = GetTime();
        }
        if (warped && !warpDone && GetTime() - warpStart > 8.0) {
            warpDone = true;
            const float factor = Rocket_Warp();
            Rocket_SetWarp(0);
            const RocketTelemetry a = Rocket_Telemetry(world);
            const bool held = factor >= 100.0f && fabs(a.apoapsis - warpBefore.apoapsis) < 1500.0 && fabs(a.periapsis - warpBefore.periapsis) < 1500.0;
            printf("bench: rocket: time warp %.0fx for 8 s (~%.0f min of orbit): apoapsis %.1f -> %.1f km, periapsis %.1f -> %.1f km: %s\n", factor, 8.0 * factor / 60.0,
                   warpBefore.apoapsis / 1000.0, a.apoapsis / 1000.0, warpBefore.periapsis / 1000.0, a.periapsis / 1000.0, held ? "orbit held" : "DRIFTED");
            if (!held) return -1;
        }
        if (warpDone) { Rocket_SetAutopilot(RocketAutopilot::Deorbit); printf("bench: rocket: deorbit burn\n"); Next(); }
        break;
    case 4: {
        static double lastLog = 0.0;
        if (GetTime() - lastLog > 1.0 && t.groundAltitude < 8000.0) { lastLog = GetTime(); printf("bench: rocket: descent %.0f m above ground, vertical %+.1f m/s, surface %.1f m/s, throttle %.2f, fuel %.0f\n", t.groundAltitude, t.verticalSpeed, t.surfaceSpeed, t.throttle, t.fuel); }
    }
        if (t.state == RocketState::OnPad && t.landings > 0) {
            const Vector3 landed = Vector3Normalize(p->toLocal((player.position - world.renderOrigin).toFloat()));
            const double km = acos(std::clamp(Vector3DotProduct(landed, padDir), -1.0f, 1.0f)) * pf->radius / 1000.0;
            printf("bench: rocket: LANDED on %s after %.0f s: touchdown %.1f m/s, %d parts broken, %.0f km from the pad, fuel %.0f kg left (highest %.1f km, hottest %.0f C)\n",
                   t.body.c_str(), GetTime() - testStart, t.lastTouchdown, t.damage, km, t.fuel, maxAlt / 1000.0, Rocket_MaxTemperature());
            Next();
        } else if (T() > 300.0) { printf("bench: rocket: did not land (altitude %.0f m, vertical %.1f m/s, fuel %.0f kg)\n", t.groundAltitude, t.verticalSpeed, t.fuel); return -1; }
        break;
    case 5: {
        if (T() < 2.0) break;
        TakeScreenshot("bench_rocket_landed.png");
        // Life support spent some of the cockpit's air above the atmosphere (none on the ground).
        const float air = Rocket_AirMinutes();
        const bool breathed = air < airLoaded - 0.05f && air > 0.0f;
        printf("bench: rocket: life support: %.1f of %.0f min of air left%s\n", air, airLoaded, breathed ? "" : " - FAILED: the air was not used");
        const bool out = Rocket_Exit(world);
        const bool pass = orbited && t.damage == 0 && out && survival.deaths == 0 && breathed;
        printf("bench: rocket gate (build, launch, orbit, land): %s\n", pass ? "PASS" : "FAIL");
        return pass ? 1 : -1;
    }
    }
    return 0;
}

// --test-station (plan 3, M20 gate): the rocket above with a station core on top flies to orbit,
// deploys the core as a station (its own body on a Kepler orbit), docks to it, gives it fuel,
// undocks, docks again and takes some fuel back; the station and its fuel survive a save and load.
namespace {
    int sphase = 0;
    double sStart = 0.0, sPhaseStart = 0.0;
    bool sOk = true;
    float sGiven = 0.0f, sTaken = 0.0f, sFuelSaved = 0.0f;
    double sPe0 = 0.0, sPe1 = 0.0;
    void SNext() { ++sphase; sPhaseStart = GetTime(); }
    double ST() { return GetTime() - sPhaseStart; }
    void SCheck(bool c, const char* what) { if (!c) { printf("bench: station: FAILED %s\n", what); sOk = false; } }
    float StationFuel(Scene& world) {
        float f = 0.0f;
        for (auto& o : world.objects)
            if (auto* b = dynamic_cast<Planetoid*>(o.get()); b && Station_Is(*b) && b->plots)
                for (const auto& pl : b->plots->plots) {
                    bool core = false;
                    for (const Part& part : pl->parts) core = core || part.type == PartType::StationCore;
                    if (core) f += pl->fuel;
                }
        return f;
    }
}

int StationTest_Update(Scene& world) {
    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
    if (!pf) { printf("bench: --test-station needs the planet\n"); return -1; }
    const RocketTelemetry t = Rocket_Telemetry(world);
    switch (sphase) {
    case 0: {
        sStart = GetTime();
        const Vector3 home = Vector3Normalize(p->toLocal((player.position - world.renderOrigin).toFloat()));
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 }));
        const Vector3 pad = Vector3Normalize(Vector3Add(home, Vector3Scale(t1, 14.0f / pf->radius)));
        if (!p->plots) p->plots = std::make_unique<PlotSet>();
        Plot* plot = p->plots->plotFor(Vector3Scale(pad, pf->radius + pf->elevation(pad)), 0.0f, true);
        bool ok = BuildTool_Place(*p, *plot, { PartType::RocketEngine, { 0, 0, 0 }, 0 }, false);
        const Int3 legs[4] = { { -1, 0, -1 }, { 2, 0, -1 }, { -1, 0, 2 }, { 2, 0, 2 } };
        for (const Int3& l : legs) ok &= BuildTool_Place(*p, *plot, { PartType::LandingLeg, l, 0 }, false);
        for (int k = 0; k < 4; ++k) ok &= BuildTool_Place(*p, *plot, { PartType::FuelTank, { 0, 2 + 2 * k, 0 }, 0 }, false);
        ok &= BuildTool_Place(*p, *plot, { PartType::Guidance, { 0, 10, 0 }, 0 }, false);
        ok &= BuildTool_Place(*p, *plot, { PartType::Cockpit, { 0, 11, 0 }, 0 }, false);
        ok &= BuildTool_Place(*p, *plot, { PartType::StationCore, { 0, 13, 0 }, 0 }, false);
        inventory.add("rocket_fuel", 75);
        Rocket_LoadFuel(*plot, 75);
        const RocketStats s = Rocket_Stats(*plot, p);
        printf("bench: station: rocket of %d parts with a station core: fuel %.0f of %.0f kg, TWR %.2f, delta-v %.0f m/s: %s\n",
               plot->partCount(), s.fuel, s.fuelCapacity, s.twr, s.deltaV, s.ready ? "READY" : s.missing.front().c_str());
        if (!ok || !s.ready || !Rocket_Board(world, *p, *plot)) return -1;
        SNext();
        break;
    }
    case 1:
        if (ST() > 1.0) { Rocket_SetAutopilot(RocketAutopilot::Ascent, 45000.0); SNext(); }
        break;
    case 2: // to orbit
        if (t.state == RocketState::Flying && t.periapsis > 35000.0 && t.apoapsis > 35000.0) {
            printf("bench: station: in orbit after %.0f s (apoapsis %.1f km, periapsis %.1f km, fuel %.0f kg)\n", GetTime() - sStart, t.apoapsis / 1000.0, t.periapsis / 1000.0, t.fuel);
            SNext();
        } else if (ST() > 240.0) {
            printf("bench: station: no orbit (state %d, apoapsis %.1f km, periapsis %.1f km, altitude %.1f km, fuel %.0f kg, autopilot %d)\n",
                   (int)t.state, t.apoapsis / 1000.0, t.periapsis / 1000.0, t.altitude / 1000.0, t.fuel, (int)t.autopilot);
            return -1;
        }
        break;
    case 3: { // deploy
        if (ST() < 1.0) break;
        SCheck(Rocket_DeployStation(), "the core deployed");
        double ap = 0.0;
        SCheck(Station_Count() == 1 && Station_Orbit(0, ap, sPe0) && sPe0 > 30000.0, "a station on a stable orbit");
        printf("bench: station: founded (orbit %.1f x %.1f km); station tanks hold %.0f kg\n", ap / 1000.0, sPe0 / 1000.0, StationFuel(world));
        SNext();
        break;
    }
    case 4: // dock (it is right beside us, at the same speed)
        if (ST() < 1.0) break;
        SCheck(Rocket_Dock() && Rocket_DockedAtStation(), "docked");
        SNext();
        break;
    case 5: { // give fuel, undock
        if (ST() < 0.5) break;
        const float before = StationFuel(world);
        sGiven = -Rocket_TransferFuel(-400.0f);
        SCheck(sGiven > 300.0f && StationFuel(world) > before + 300.0f, "fuel given to the station");
        SCheck(Rocket_Undock() && Rocket_Telemetry(world).state == RocketState::Flying, "undocked");
        SNext();
        break;
    }
    case 6: // drift a few seconds beside it, dock again, take fuel back
        if (ST() < 3.0) break;
        SCheck(Rocket_Dock(), "docked again");
        sTaken = Rocket_TransferFuel(200.0f);
        SCheck(sTaken > 150.0f, "fuel taken from the station");
        SNext();
        break;
    case 7: { // save and load the station
        if (ST() < 0.5) break;
        sFuelSaved = StationFuel(world);
        const std::string path = std::string(GetWorkingDirectory()) + "/station_test.tmp";
        if (FILE* f = fopen(path.c_str(), "w")) { saveExtrasHook(f); fprintf(f, "end_extras\n"); fclose(f); }
        Rocket_Exit(world);
        Station_Reset();
        SCheck(Station_Count() == 0, "stations cleared");
        if (FILE* f = fopen(path.c_str(), "r")) { loadExtrasHook(f); fclose(f); }
        remove(path.c_str());
        double ap = 0.0;
        SCheck(Station_Count() == 1 && Station_Orbit(0, ap, sPe1) && fabs(sPe1 - sPe0) < 500.0, "the station and its orbit loaded");
        SCheck(fabsf(StationFuel(world) - sFuelSaved) < 1.0f, "its fuel loaded");
        printf("bench: station: gave %.0f kg, undocked and docked again, took %.0f kg back; saved and loaded (periapsis %.1f -> %.1f km, %.0f kg in its tanks) -> %s\n",
               sGiven, sTaken, sPe0 / 1000.0, sPe1 / 1000.0, StationFuel(world), sOk ? "PASS" : "FAIL");
        return sOk ? 1 : -1;
    }
    }
    return 0;
}
