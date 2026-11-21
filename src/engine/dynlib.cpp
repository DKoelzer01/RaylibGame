#include "dynlib.h"

// Its own file: windows.h and raylib.h cannot share a translation unit (CloseWindow, ShowCursor, Rectangle, ...).
#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

void* DynLib_Open(const char* name) {
#if defined(_WIN32)
    return (void*)LoadLibraryA(name);
#else
    return dlopen(name, RTLD_NOW | RTLD_LOCAL);
#endif
}

void* DynLib_Symbol(void* lib, const char* name) {
    if (!lib) return nullptr;
#if defined(_WIN32)
    return reinterpret_cast<void*>(GetProcAddress((HMODULE)lib, name));
#else
    return dlsym(lib, name);
#endif
}

void DynLib_Close(void* lib) {
    if (!lib) return;
#if defined(_WIN32)
    FreeLibrary((HMODULE)lib);
#else
    dlclose(lib);
#endif
}
