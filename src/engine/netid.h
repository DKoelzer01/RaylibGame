#ifndef NETID_H
#define NETID_H

// Stable network ids (#317, multiplayer phase 1). Replication messages refer to things by these, never by pointers,
// vector indices or the runtime Planetoid::id (a per-run counter). A NetId is a kind in the top 8 bits and a value:
//   Body      the body's place in the system definition's creation order (= its index in the save): the same in every
//             run of a seed and on every machine that generated the world from it
//   RingRock  ring name hash (24 bits) << 32 | the rock's index in its ring (promoted ring rocks)
//   Animal, Vehicle, Station   allocated from a per-world counter (NetId_Allocate), saved with the world, never reused
//   Player    0 = the local player; a server issues ids (or SteamIDs) later
// Plots are (body, plot id) and parts (plot, part index): both already stable and saved (plots.cpp); machines are parts.
// Plants keep their own stable ids (tile + cell).

#include <cstdint>
#include <cstdio>

using NetId = uint64_t; // 0: none

enum class NetKind : uint8_t { None = 0, Body = 1, RingRock = 2, Animal = 3, Vehicle = 4, Station = 5, Player = 6 };

constexpr uint64_t kNetIdValueMask = 0x00FFFFFFFFFFFFFFull;
constexpr NetId NetId_Make(NetKind kind, uint64_t value) { return ((uint64_t)kind << 56) | (value & kNetIdValueMask); }
constexpr NetKind NetId_Kind(NetId id) { return (NetKind)(id >> 56); }
constexpr uint64_t NetId_Value(NetId id) { return id & kNetIdValueMask; }

NetId NetId_Allocate(NetKind kind);  // the next free id (Animal, Vehicle, Station)
void NetId_Note(NetId id);           // a loaded id: later allocations stay above it
void NetId_Reset();                  // a new world: allocation starts again at 1
void NetId_Save(FILE* f);            // the counter ("netid_next N")
bool NetId_LoadLine(const char* tok, FILE* f);
uint64_t NetId_Next();
uint32_t NetId_NameHash(const char* name); // 24-bit FNV-1a (ring rocks)

struct PlotRef { NetId body = 0; uint32_t plot = 0; };
struct PartRef { NetId body = 0; uint32_t plot = 0; int32_t part = -1; };

#endif
