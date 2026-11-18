#pragma once
#include <cstdlib>
#include <cstring>

// An on/off environment switch (#751): on when set to anything but empty, "0", "false" or "off". Testing only whether a
// variable is set made FLUID_CPU=0 select the CPU fluid kernel all the same.
inline bool EnvFlag(const char* name) {
    const char* e = std::getenv(name);
    if (!e || !e[0]) return false;
    return std::strcmp(e, "0") != 0 && std::strcmp(e, "false") != 0 && std::strcmp(e, "off") != 0;
}
