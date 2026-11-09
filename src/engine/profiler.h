#ifndef PROFILER_H
#define PROFILER_H

// Tracy profiler hooks. Compiled out unless the build defines TRACY_ENABLE
// (CMake: -DENABLE_TRACY=ON or the "profile" preset).
// Run the Tracy GUI (tracy-profiler.exe, same version) and it connects to the game automatically.
#ifdef TRACY_ENABLE
    #include <tracy/Tracy.hpp>
    #define PROFILE_CONCAT_(a, b) a##b
    #define PROFILE_CONCAT(a, b) PROFILE_CONCAT_(a, b)
    // Unique variable per line, so several zones can share one scope.
    #define PROFILE_SCOPE(name) ZoneNamedN(PROFILE_CONCAT(profileZone_, __LINE__), name, true)
    #define PROFILE_FUNCTION()  ZoneScoped
    #define PROFILE_FRAME()     FrameMark
#else
    #define PROFILE_SCOPE(name)
    #define PROFILE_FUNCTION()
    #define PROFILE_FRAME()
#endif

#endif
