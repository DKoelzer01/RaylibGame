#ifndef GPUMEM_H
#define GPUMEM_H

// VRAM accounting (#302): with --vram, the GL allocation calls (buffer/texture/renderbuffer storage and their deletes)
// are wrapped through glad's function pointers and every object is counted with its size and the subsystem tag that
// was current when it was allocated (GpuMemTag scopes at the main call sites; "other" otherwise). GpuMem_Report prints
// the totals per tag and what the driver says is in use (NVX/ATI memory info, when available). Off (no wrappers, no
// cost) unless GpuMem_Install is called.

#include <cstddef>

void GpuMem_Install();          // after the GL context exists
bool GpuMem_Active();
void GpuMem_Report(const char* label); // "bench: vram: ..." lines
size_t GpuMem_Total();

struct GpuMemTag {               // the subsystem allocations in this scope are counted under
    const char* prev;
    explicit GpuMemTag(const char* tag);
    ~GpuMemTag();
};

#endif
