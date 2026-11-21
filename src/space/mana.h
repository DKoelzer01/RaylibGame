#pragma once

// Mana kinds and their algebra (#1001, the magic doc section 3.1; #1248). The kinds come from resources/data/mana.txt:
// four base kinds, compounds of two different bases (combine), inverses of a base or Life (invert) and Aether (cancel: a
// kind and its inverse annihilate into it). One step deep in v1: no compound of compounds, at most 16 kinds.
// Everything magic is behind RAYLIBGAME_MAGIC until the first playable slice (magic.h).

#include <cstdint>
#include <string>
#include <vector>

constexpr int kManaMaxKinds = 16; // decided (2026-11-21): 16 kinds max in v1

enum class ManaOp : uint8_t { Base, Combine, Invert, Cancel };

struct ManaKind {
    std::string id, name, latin;
    ManaOp op = ManaOp::Base;
    int a = -1, b = -1;        // Combine: the two bases; Invert: the kind it inverts (a); else -1
    std::string status;        // the combat status effect ("" none)
    uint32_t colour = 0xFFFFFF; // 0xRRGGBB
    bool slice = false;        // in the first playable slice
};

// The table: loaded once from resources/data/mana.txt (or a test's text). Load replaces the table only on success.
bool Mana_LoadText(const std::string& text, std::string& error);
bool Mana_Load(const char* path = "resources/data/mana.txt");
bool Mana_Loaded();
int Mana_Count();
const ManaKind& Mana_Kind(int k);
int Mana_Find(const std::string& id); // -1 if none
// Checks the algebra: unique ids, exactly four bases, each compound two different bases and each pair at most once,
// each inverse of a base or a compound and each kind inverted at most once, exactly one cancel kind, at most 16 kinds.
// The problems, one a line ("" none).
std::string Mana_Validate();

// The algebra (kind indices; -1 = no such kind).
int Mana_Combine(int a, int b);  // A + B: the compound of two different bases; A + ~A: Aether; else -1
int Mana_Inverse(int k);         // ~A (symmetric: the inverse of Frost is Fire); -1 if it has none
bool Mana_Cancels(int a, int b); // a and b are each other's inverse
int Mana_Aether();               // the cancel kind
// Conversion rates (mana.txt "rate" lines): combine, split, invert, cancel, alone, from_aether. 0 if unknown.
float Mana_Rate(const std::string& name);
