#ifndef DYNLIB_H
#define DYNLIB_H

// A shared library loaded at run time (Steam Audio, #1292): LoadLibrary / dlopen. nullptr when it isn't there.
void* DynLib_Open(const char* name);
void* DynLib_Symbol(void* lib, const char* name);
void DynLib_Close(void* lib);

#endif
