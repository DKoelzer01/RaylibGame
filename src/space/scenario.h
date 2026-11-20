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
    void waitingFor(const std::string& what); // the heartbeat says so
    void keepCapture(const std::string& reason); // keep the recording whatever the verdict (a look change for the user)
    void cue();                               // the recording's first frame from now (--record-on-cue)

    // Helpers the fluid and terrain tests share.
    Vector3 spawnDirLocal() const;            // the start's direction, body-local (--spawn-dir override, else the field's)
    Vector3 spawnDirWorld() const;            // the same in the world frame (main.cpp's StartSpawnDir)
    // Camera at body-local points, the player flying there (the subject in view: CLAUDE.md).
    void camera(Vector3 eyeLocal, Vector3 targetLocal, Vector3 upLocal);
};

struct ScenarioStep {
    std::string name;
    double deadline = 60.0;                   // sim seconds after the step began
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
void Scenario_Shutdown();
bool Scenario_Active();

#endif
