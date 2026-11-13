#ifndef GPUCOMPUTE_H
#define GPUCOMPUTE_H

// GPU compute (OpenGL 4.3 compute shaders + SSBOs through rlgl). The game is built for a 4.3 context;
// The game requires a 4.3 context (#185: no window without it), so Compute_Available() is normally true; callers keep their CPU
// paths for the off switches and for a compute shader that fails to compile.
//
// Shaders live in resources/compute/*.comp and may `#include "file"` (relative to the shader), the same
// expansion HotShader does. Buffers are SSBOs bound by index (layout(std430, binding = N)).

#include <string>
#include <vector>

bool Compute_Available();

// Expands `#include "file"` lines (relative to the including file); the included paths are appended.
std::string LoadShaderSource(const std::string& path, std::vector<std::string>& includes);

struct ComputeProgram {
    unsigned int id = 0;
    std::string path;
    bool load(const std::string& path);      // compiles; false (and logged) on error
    void unload();
    void use() const;                       // rlEnableShader
    int uniform(const char* name) const;    // location
    void setInt(const char* name, int v) const;
    void setUint(const char* name, unsigned v) const;
    void setFloat(const char* name, float v) const;
    void setVec3(const char* name, float x, float y, float z) const;
    void setVec4(const char* name, float x, float y, float z, float w) const;
    void dispatch(unsigned gx, unsigned gy = 1, unsigned gz = 1) const; // then a storage barrier
};

struct ComputeBuffer {
    unsigned int id = 0;
    unsigned int bytes = 0;
    void create(unsigned int size, const void* data = nullptr, bool dynamic = true);
    void resize(unsigned int size);         // contents lost
    void upload(const void* data, unsigned int size, unsigned int offset = 0);
    void download(void* dest, unsigned int size, unsigned int offset = 0) const; // stalls: use sparingly
    void bind(unsigned int index) const;
    void destroy();
};

// glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_BUFFER_UPDATE_BARRIER_BIT) and friends.
void Compute_Barrier();

#endif
