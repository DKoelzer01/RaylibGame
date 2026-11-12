#ifndef ROCKET_H
#define ROCKET_H

// Rocketry (milestone 7, design doc section 7). A rocket is a build plot with a cockpit: the
// player builds it freely from parts (rocket engines, fuel tanks, a cockpit, a guidance computer,
// landing legs, any structure) whose stats come from resources/data/ship_parts.txt.
//
// Launch check (looking at the cockpit, or aboard): a cockpit, a guidance computer, at least one
// engine, fuel loaded (I, from rocket_fuel items), thrust-to-weight above 1.2 here, and delta-v
// for orbit (tanks and engine efficiency, against this body's orbital speed plus losses). The panel
// lists what is missing.
// Flight: E boards. Throttle up (Shift/Ctrl, X cuts) and the plot detaches from the ground and flies:
// summed gravity, atmospheric drag against the turning air, fuel burnt at the engines' efficiency,
// all in double precision (coasting above the air is a real Kepler orbit). The guidance computer
// turns it (W/S pitch, A/D yaw, Z/C roll; it holds still when you let go). HUD: altitude, speeds,
// apoapsis / periapsis, fuel, thrust-to-weight.
// Landing: a gentle, upright touchdown turns it back into a plot on whatever it landed on (snapped
// to that ground's grid); a hard one breaks the landing legs, a harder one the engines too (damage,
// not destruction). E leaves the cockpit once landed. Saves keep a rocket in flight.

#include <cstdio>
#include <string>
#include <vector>

class Scene;
class Planetoid;
class Plot;

struct RocketStats {
    float dryMass = 0, fuel = 0, fuelCapacity = 0, thrust = 0, isp = 0; // kg, N, s
    int engines = 0, tanks = 0, cockpits = 0, guidance = 0, legs = 0;
    float damage = 0;                     // mean damage of the rocket parts, 0..1
    int damagedParts = 0;                 // rocket parts with any damage
    float leak = 0;                       // kg of fuel a second lost from damaged tanks
    float airCapacity = 0;                // minutes of air the cockpits (10) and cabins (30) hold
    int seats = 0;                        // a cockpit 1, a cabin 4
    float twr = 0, deltaV = 0, deltaVNeeded = 0;
    bool ready = false;
    std::vector<std::string> missing;
};
// Stats of a plot as a rocket on (or near) a body.
RocketStats Rocket_Stats(const Plot& plot, const Planetoid* body);
// Repair kits (from the hotbar): the most damaged part of the nearest rocket within 12 m (the boarded
// one, or one on a pad) is half repaired. False if there is none to repair.
bool Rocket_UseRepairKit();
// Life support: oxygen canisters loaded into a rocket (I, with the fuel), 10 minutes of air each.
int Rocket_LoadAir(Plot& plot, int canisters);
// Aboard outside breathable air, the cockpit spends its air: true while there is some.
bool Rocket_BreatheAboard(float dt);
float Rocket_AirMinutes();
// Stations (plan 3, M20; station.h). In a stable orbit with a station core aboard, K deploys it as a
// station; near a station (within 40 m, under 2 m/s) K docks. Docked: T / Shift+T move fuel between the
// rocket and the station's tanks (kg > 0: into the rocket), and launching (or U) undocks.
bool Rocket_DeployStation();
bool Rocket_Dock();
bool Rocket_Undock();
float Rocket_TransferFuel(float kg);
bool Rocket_DockedAtStation();
// Tests: damage a rocket part (part index) by an amount.
void Rocket_DamagePart(Plot& plot, int part, float amount);

enum class RocketState { None, OnPad, Flying };
enum class RocketAutopilot { Off, Ascent, Circularize, Deorbit, Land };

struct RocketTelemetry {
    RocketState state = RocketState::None;
    double altitude = 0, groundAltitude = 0;   // above the reference body's radius / the ground
    double verticalSpeed = 0, surfaceSpeed = 0, orbitalSpeed = 0;
    double apoapsis = 0, periapsis = 0;        // altitudes; apoapsis < 0: escaping
    float fuel = 0, throttle = 0, twr = 0;
    std::string body;
    RocketAutopilot autopilot = RocketAutopilot::Off;
    int landings = 0;
    float lastTouchdown = 0;                   // m/s
    int damage = 0;                            // parts broken so far
};

void Rocket_Init(Scene* world);                // ship_parts.txt; the world (for loading a rocket in flight)
bool Rocket_Boarded();
// Time warp while coasting above the air (period / comma), dropping to 1x for burns and in the air.
float Rocket_Warp();
void Rocket_SetWarp(int level);    // 0 = 1x .. 6 = 100x (tests)
bool Rocket_MapView();             // M: the orbit map is up
void Rocket_SetMapView(bool on);
float Rocket_MaxTemperature();     // hottest the leading face got (tests)
bool Rocket_AimingCockpit();                   // on foot, looking at a rocket's cockpit (E boards it)
// Once a frame. On foot: the launch-check panel, E boards, I loads fuel. Aboard: flight and camera.
void Rocket_Update(Scene& world, bool inputEnabled, float dt);
void Rocket_DrawHUD(Scene& world);
RocketTelemetry Rocket_Telemetry(Scene& world);

// Tests and scripted use.
bool Rocket_Board(Scene& world, Planetoid& body, Plot& plot);
bool Rocket_Exit(Scene& world);
void Rocket_SetThrottle(float t);
void Rocket_SetAutopilot(RocketAutopilot mode, double targetAltitude = 45000.0);
int Rocket_LoadFuel(Plot& plot, int items);    // from the inventory; returns items used

void Rocket_Save(FILE* f);                     // game save extras
bool Rocket_LoadLine(const char* tok, FILE* f);  // true if it consumed a rocket line
void Rocket_Reset();                           // a new world / load

#endif
