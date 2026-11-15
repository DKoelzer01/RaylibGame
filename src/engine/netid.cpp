#include "netid.h"

#include <algorithm>
#include <cstring>

namespace {
uint64_t next = 1; // allocated values share one counter across kinds (unique per world)
}

NetId NetId_Allocate(NetKind kind) { return NetId_Make(kind, next++); }
void NetId_Note(NetId id) { next = std::max(next, NetId_Value(id) + 1); }
void NetId_Reset() { next = 1; }
uint64_t NetId_Next() { return next; }
void NetId_Save(FILE* f) { fprintf(f, "netid_next %llu\n", (unsigned long long)next); }

bool NetId_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "netid_next") != 0) return false;
    unsigned long long n = 1;
    if (fscanf(f, " %llu", &n) == 1) next = std::max(next, (uint64_t)n);
    return true;
}

uint32_t NetId_NameHash(const char* name) {
    uint32_t h = 2166136261u;
    for (const char* c = name; *c; ++c) { h ^= (uint8_t)*c; h *= 16777619u; }
    return h & 0xFFFFFFu;
}
