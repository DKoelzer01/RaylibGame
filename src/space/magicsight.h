#pragma once

// Magic sight (#1249, magic step 2; the magic doc section 2): a permanent buff cast once on yourself (the First Sight page,
// #1250). It stays through death, sleep and reloads. A key toggles its overlay (ley lines as ribbons coloured by kind, nodes
// as pulsing glows), and while shown it holds back 10% of the mana pool (decided 2026-11-21: free while hidden). Holding the
// key while looking at a node, a creature or a rock studies it without destroying it: full insight the first time per
// type, less after. The aspect points it gives land with the tome (#1251); until then studies are counted per type.
// Attunement pages (later) widen the range (50 m, 150 m, 400 m) and show lines through rock.
//
// Behind RAYLIBGAME_MAGIC; until discovery (#1250) teaches First Sight, RAYLIBGAME_MAGIC_SIGHT=1 grants it (=show: shown too).

#include <map>
#include <string>

class Scene;
struct PlayerInput;

struct MagicPlayerState {
    bool sight = false;       // knows magic sight (permanent)
    bool shown = false;       // its overlay is on
    int sightTier = 0;        // attunement: 0 = 50 m, 1 = 150 m, 2 = 400 m
    bool throughRock = false; // attunement: lines seen through rock
    bool found = false;       // (#1250) has stumbled on magic (the first find: the Arcane line opens)
    float pool = 100.0f, poolMax = 100.0f; // the Aether pool (spells come in step 6, #1253)
    std::map<std::string, int> studied;    // study target type ("node:storm", "creature:deer", "rock:granite") -> times
};

constexpr float kSightReserve = 0.10f;   // of the pool, while the overlay is shown (decided)
constexpr float kStudySeconds = 3.0f;    // looking and holding the key (the doc: a few seconds)
constexpr float kSightTapSeconds = 0.3f; // a shorter press toggles the overlay; a longer one studies

// The rules (pure).
float MagicSight_Range(const MagicPlayerState& s);     // m
float MagicSight_Reserved(const MagicPlayerState& s);  // pool held back now
float MagicSight_Available(const MagicPlayerState& s); // pool free for spells
bool MagicSight_Learn(MagicPlayerState& s);            // false if already known
bool MagicSight_SetShown(MagicPlayerState& s, bool on); // false (and stays hidden) if not known
float MagicSight_StudyValue(const MagicPlayerState& s, const std::string& type); // 1 the first time, then halving, at least 0.1
float MagicSight_Study(MagicPlayerState& s, const std::string& type);            // counts it; the value it gave
// The pool refills over time, faster on ley lines (the doc: "regenerating faster on ley lines"); flow = the strongest
// base kind's line flow here (0..2).
void MagicSight_Regen(MagicPlayerState& s, float dt, float flow);
// The player record's lines (magic_sight, magic_pool, magic_studied).
std::string MagicSight_SaveText(const MagicPlayerState& s);
bool MagicSight_LoadLine(MagicPlayerState& s, const char* line);

// The game (the local player): the key, studying, the pool; the ribbons; the HUD.
void MagicSight_Update(Scene& world, float dt, const PlayerInput& in);
void MagicSight_Draw3D(Scene& world);
void MagicSight_DrawHUD(Scene& world);
