#ifndef SCENARIO_H
#define SCENARIO_H

// Tier 2 scenarios (#870, testing framework #868): gameplay checked in a running world, one .cpp per scenario under
// tests/scenarios/, registered with SCENARIO(...). main.cpp knows nothing about individual scenarios; it calls four hooks:
//   Scenario_PrepareArgs(argc, argv)  first thing in main(): --scenario NAME adds the scenario's own flags to argv (a bench
//                                     with fixed dt, its seed, sandbox or planet, light streaming, a temporary recording)
//   Scenario_Setup()                  after the bench save path is set: the save goes to the run's scratch dir instead
//   Scenario_Frame(world)             each frame while a bench waits to start: runs the scenario; true while it runs
//                                     (holds the bench start, as the --test-* flags do)
//   Scenario_Shutdown()               after Recorder_Shutdown: results.json final, the capture kept or deleted
//
// Command line:
//   --scenario NAME | tag:TAG | all     run one scenario (in this process) or every match (each in its own process)
//   --scenario-list [tag:TAG]          list the scenarios and exit
//   --scenario-out DIR                 where results.json and captures go (default: %TEMP%/raylibgame_scenarios/<time>)
//   --scenario-keep                    keep the recording even on PASS (evidence for a PR)
//   --scenario-no-record               no recording
//   --scenario-report IN [OUT.html]    results.json (or a dir of them) -> an HTML report; appends to the history
//                                      (RAYLIBGAME_TEST_HISTORY, default //DEWEYSERVER/files/claude-media/raylib-test-history)
//
// A scenario is a list of steps. Each step runs every frame until it returns true (done) or its deadline (sim seconds since
// the step began) passes, which fails the scenario. Waits are on state, never on the wall clock alone; the whole scenario
// has a wall-clock timeout too. A heartbeat line every 5 s says which step runs and what it waits for. Checks
// (c.eq / c.near / c.lt / c.gt / c.between / c.check) record the value, the expectation and the verdict in results.json;
// the scenario passes when every check passed and every step finished.

#include <functional>
#include <string>
#include <vector>
#include <raylib.h>

class Scene;
class Planetoid;
class PlanetField;

enum class ScenarioWorld { Sandbox, Planet };

struct ScenarioCheck { std::string name, op, expected, detail; double value = 0.0; bool pass = false; };
struct ScenarioMetric { std::string name, unit; double value = 0.0; bool lowerIsBetter = true; };

struct ScenarioCtx {
    Scene& world;
    Planetoid* body = nullptr;            // world.objects[0]
    const PlanetField* field = nullptr;   // its PlanetField (planet and sandbox worlds have one)
    double t = 0.0;                       // sim seconds since the current step began
    double T = 0.0;                       // sim seconds since the scenario began
    double wall = 0.0;                    // wall seconds since the current step began
    double clock = 0.0;                   // what the step's deadline counts: t, or the lesser of t and wall (realTime, #1972),
                                          // minus the sim time held by world work (#2119: the step waits while the world
                                          // streams or builds, or a loading screen ran)
    long frame = 0;                       // frames since the scenario began
    int stepIndex = 0;

    // Checks: recorded in results.json and printed ("scenario: check NAME = V (want ...) ok|FAIL"). Return the verdict.
    bool check(const std::string& name, bool ok, const std::string& detail = "");
    bool eq(const std::string& name, double v, double want, double tol = 0.0);
    bool near(const std::string& name, double v, double want, double tol);
    bool timingLt(const std::string& name, double v, double bound); // a timing gate (#1034): reported only where RAYLIBGAME_TIMING_GATES=0
    bool lt(const std::string& name, double v, double bound);
    bool gt(const std::string& name, double v, double bound);
    bool between(const std::string& name, double v, double lo, double hi);
    void metric(const std::string& name, double v, const std::string& unit = "", bool lowerIsBetter = true); // reported only
    void note(const std::string& text);       // printed and kept in results.json
    void fail(const std::string& reason);     // ends the scenario now: FAIL
    void waitingFor(const std::string& what); // the heartbeat says so; while the world isn't idle, the deadline pauses (#2119)
    void keepCapture(const std::string& reason); // keep the recording whatever the verdict (a look change for the user)
    void cue();                               // the recording's first frame from now (--record-on-cue)
    // Visual checks (#870 visual tier): grab the screen now as NAME.png (in the run's shots/ folder) and compare it with the
    // baseline on the share (Scenario_BaselineDir()/<scenario>/NAME.png), both shrunk to 320x180: the check fails when more
    // than tolerancePct % of the pixels differ by more than 24 of 255 in any channel. No baseline yet, or
    // RAYLIBGAME_BASELINE_UPDATE=1: the shot becomes the baseline (the check passes, noted). On a failure the shot and a
    // difference image stay in shots/ and the recording is kept.
    bool shot(const std::string& name, double tolerancePct = 3.0);
    // (#1071) A picture for the menu slideshow, not a check: only when RAYLIBGAME_SLIDE_SHOTS names a folder, the screen
    // goes there once as <scenario>-<name>.png (the scenario's '/' as '_'); otherwise nothing happens. Scenarios that
    // already frame their subject well call it at that moment (run them with --no-hud --native-window).
    void grab(const std::string& name);

    // Helpers the fluid and terrain tests share.
    Vector3 spawnDirLocal() const;            // the start's direction, body-local (--spawn-dir override, else the field's)
    Vector3 spawnDirWorld() const;            // the same in the world frame (main.cpp's StartSpawnDir)
    // Camera at body-local points, the player flying there (the subject in view: CLAUDE.md).
    void camera(Vector3 eyeLocal, Vector3 targetLocal, Vector3 upLocal);
};

struct ScenarioStep {
    std::string name;
    double deadline = 60.0;                   // on ScenarioCtx::clock: sim seconds after the step began, less any held by
                                              // world work while the step waits (#2119)
    std::function<bool(ScenarioCtx&)> run;    // true: done
};

struct ScenarioDef {
    std::string name, file;
    std::vector<std::string> tags;
    ScenarioWorld world = ScenarioWorld::Sandbox;
    unsigned seed = 42;
    double timeout = 180.0;                   // wall seconds for the whole scenario (loading excluded)
    std::vector<std::string> args;            // extra game flags (e.g. --day)
    std::vector<std::pair<std::string, std::string>> env; // environment set before the game starts
    std::string legacyFlag;                   // the --test-* flag it replaces (the report pairs their verdicts)
    bool record = true;                       // the temporary recording (perf scenarios: off, it costs frame time)
    bool needsPlanetField = true;             // body 0 must have a PlanetField (false: any body, e.g. the asteroid system, #667)
    // A legacy bench drives the run (its flag in args, e.g. --shot-graves): the scenario has no steps and never holds the
    // bench; when the game ends, every bench_*.png it wrote (or just `shots`, if listed) is compared with its baseline, as
    // c.shot does, and those checks are the verdict. The way to move a screenshot bench without rewriting it (#870).
    bool legacyDrives = false;
    std::vector<std::string> shots;           // (legacyDrives) the files to compare, e.g. "bench_graves.png"; empty: all new
    double shotTolerance = 3.0;               // (legacyDrives) % of differing pixels allowed
    // (#1442) Every scenario reports camera_inside_solid_frames / player_inside_solid_frames (cliptelemetry.h); with this
    // set, a frame with the camera or the player inside solid ground or a built cube turns a PASS into a FAIL.
    bool expectNoClip = false;
    // (#1972) The steps wait on wall-clock work (a join over real sockets, the factory's 20 Hz thread): a headless sim runs
    // up to 50x real time, so a sim-time deadline expired before the work could finish, and a loaded machine runs it slower
    // than real time. With this set, the step deadlines and the timeout count ScenarioCtx::clock, the lesser of sim and wall.
    bool realTime = false;
    std::vector<ScenarioStep> steps;
    void step(const std::string& name, double deadline, std::function<bool(ScenarioCtx&)> run) { steps.push_back({ name, deadline, std::move(run) }); }
    bool hasTag(const std::string& t) const { for (const auto& x : tags) if (x == t) return true; return false; }
};

using ScenarioBuild = void (*)(ScenarioDef&);
struct ScenarioRegistrar { ScenarioRegistrar(const char* name, const char* file, ScenarioBuild build); };

#define SCENARIO_CAT2(a, b) a##b
#define SCENARIO_CAT(a, b) SCENARIO_CAT2(a, b)
// SCENARIO("area/name") { s.tags = { ... }; s.world = ...; s.step("...", deadline, [](ScenarioCtx& c) { ...; return done; }); }
#define SCENARIO(NAME)                                                                                                  \
    static void SCENARIO_CAT(ScenarioBuild_, __LINE__)(ScenarioDef & s);                                              \
    static ScenarioRegistrar SCENARIO_CAT(ScenarioReg_, __LINE__)(NAME, __FILE__, SCENARIO_CAT(ScenarioBuild_, __LINE__)); \
    static void SCENARIO_CAT(ScenarioBuild_, __LINE__)(ScenarioDef & s)

// The hooks (see the top).
void Scenario_PrepareArgs(int& argc, char**& argv);
void Scenario_Setup();
bool Scenario_Frame(Scene& world);
void Scenario_Loading(); // (#1422) each frame of a loading screen: heartbeats and the wall backstop for a run held there
void Scenario_Shutdown();
bool Scenario_Active();
// (#1216) A scenario drives gameplay input with scripted holds (gFakeHeld): the bench run takes gameplay input while it is set
// (and the window need not have focus). Scenarios that set it clear it again.
inline bool gScenarioInput = false;
// (#1105) A scenario places the camera itself (world.camera, e.g. a view of the player from outside, for its recording):
// while set, the player's first-person view does not overwrite it each frame. Scenarios that set it clear it again.
inline bool gScenarioCamera = false;
std::string Scenario_BaselineDir(); // RAYLIBGAME_BASELINES, default //DEWEYSERVER/files/claude-media/raylib-baselines
std::string Scenario_RunDir();      // (#1693) this run's folder (results.json, shots/): where saves/golden-write puts golden-v<N>/

#endif
