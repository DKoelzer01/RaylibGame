#pragma once

// Aspects and the tome (#1251, magic step 4; the magic doc sections 2 and 7). Things are made of aspects
// (resources/data/aspects.txt); studying them with magic sight (or, later, spending them at the research table) gives the
// player their aspects. Tome pages (resources/data/tome.txt) cost aspect combinations and are kept per team; aspects are
// per player (decided). Pages cross-link with the research tree both ways: a research node may need a page (research.txt
// "tome"), a page may need a research node. Behind RAYLIBGAME_MAGIC like the rest of magic.

#include <map>
#include <set>
#include <string>
#include <vector>

struct MagicPlayerState;
struct ResearchState;
struct ResearchNode;

struct AspectDef { std::string id, name, latin; unsigned colour = 0xFFFFFF; };

struct TomePage {
    std::string id, name, latin, branch, text;
    std::map<std::string, int> cost;     // aspect -> amount
    std::vector<std::string> needs;      // pages
    std::vector<std::string> research;   // research nodes (the cross link)
    std::vector<std::string> unlocks;    // "machine:<id>", "recipe:<id>" (later steps)
    int sightTier = -1;                  // effect: magic sight's range tier (-1 none)
    bool throughRock = false;            // effect: lines seen through rock
    float minBurden = 0.0f;              // (#1255) shown and writable only at this burden or more (the late Void and Rot pages: 20)
};

// The tables. Load replaces a table only on success.
bool Aspects_LoadText(const std::string& text, std::string& error);
bool Tome_LoadText(const std::string& text, std::string& error);
bool Tome_Load(const char* aspectsPath = "resources/data/aspects.txt", const char* tomePath = "resources/data/tome.txt");
const std::vector<AspectDef>& Aspects_All();
const AspectDef* Aspects_Find(const std::string& id);
const std::map<std::string, int>* Aspects_Of(const std::string& type); // "item:iron_ingot", "creature:..", "rock:..", "node:.."
const std::map<std::string, std::map<std::string, int>>& Aspects_Sources(); // (#1254) every "of" line: type -> aspect -> amount
const std::vector<TomePage>& Tome_Pages();
const TomePage* Tome_Find(const std::string& id);
// Checks: known aspects in every source and cost, unique ids, pages' needs exist and have no cycle, research links exist,
// and every research node's "tome" page exists (when the research table is loaded). The problems, one a line ("" none).
std::string Tome_Validate();

// A team's tome.
struct TomeState {
    std::set<std::string> pages;
    std::string Text() const;
    bool FromText(const std::string& text);
    bool Has(const std::string& page) const { return pages.count(page) > 0; }
};
TomeState& Tome_Team(int team);
void Tome_Reset();
std::string Tome_TeamText(int team);
void Tome_SetTeamText(int team, const std::string& text);

// The rules (pure).
// Studying a thing of this type with this value (1 the first time, less after: MagicSight_StudyValue) gives its aspects.
void Tome_GrantStudy(MagicPlayerState& me, const std::string& type, float value);
// Can the reader buy the page for the team? Not owned, its pages and research done, the aspects in hand; why not in `why`.
bool Tome_CanUnlock(const TomeState& team, const MagicPlayerState& me, const ResearchState& research, const TomePage& p, std::string* why = nullptr);
bool Tome_Unlock(TomeState& team, MagicPlayerState& me, const ResearchState& research, const std::string& page); // spends the aspects
// The team's pages' effects on a member's sight (range tier, through rock): never lowers what the member has.
void Tome_ApplyEffects(const TomeState& team, MagicPlayerState& me);
// The cross link: a research node that needs a tome page may start only once the team owns it.
bool Tome_ResearchAllowed(const ResearchNode& n, const TomeState& team);
// Spending an item for its aspects (the research table's job; until the table is a placeable machine, owning its page is
// enough): one unit leaves the inventory and gives its aspects in full. False (nothing spent) without the page, without
// the item, or for an item made of nothing known.
bool Tome_SpendItem(const TomeState& team, MagicPlayerState& me, class Inventory& inv, const std::string& item);
