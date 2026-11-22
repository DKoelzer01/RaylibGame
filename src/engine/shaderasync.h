#ifndef SHADERASYNC_H
#define SHADERASYNC_H

// (#1047) Compute shaders compiled without blocking: GL_KHR_parallel_shader_compile (or the ARB one) lets the driver compile
// and link on its own threads; GL_COMPLETION_STATUS_KHR says when it is done without waiting. A cold compile of
// fluid_step.comp took ~3 s here and ~37 s on the CI runner in ONE frame (a joining client's loading screen froze and its
// connection stopped being pumped). Started early (Compute_Prewarm, gpucompute.h), the work is done by the time it is needed.
// Without the extension Begin compiles nothing: the caller compiles as before (blocking).

#include <string>

struct AsyncProgram {
    unsigned shader = 0, program = 0;
    bool started() const { return program != 0; }
};

bool ShaderAsync_Available();                                  // the extension is there (and RAYLIBGAME_ASYNC_SHADERS is not 0)
AsyncProgram ShaderAsync_BeginCompute(const std::string& src); // compile + link started; {} if not available
bool ShaderAsync_Done(const AsyncProgram& p);                  // compile and link finished (never blocks)
// The linked program (0 on a compile or link error, logged with `name`); waits for the rest if it isn't done. Frees the shader.
unsigned ShaderAsync_Finish(AsyncProgram& p, const char* name);
void ShaderAsync_Discard(AsyncProgram& p);

#endif
