#include "gpucompute.h"
#include "utils/logger.h"

#include <raylib.h>
#include <rlgl.h>

extern "C" {
    typedef void (*PfnMemoryBarrier)(unsigned int);
    extern PfnMemoryBarrier glad_glMemoryBarrier;
}

bool Compute_Available() {
    return rlGetVersion() == RL_OPENGL_43 && glad_glMemoryBarrier != nullptr;
}

std::string LoadShaderSource(const std::string& path, std::vector<std::string>& includes) {
    struct Rec {
        static std::string Load(const std::string& p, std::vector<std::string>& inc, int depth) {
            char* text = LoadFileText(p.c_str());
            if (!text) return std::string();
            std::string src(text), out;
            UnloadFileText(text);
            const std::string dir = p.substr(0, p.find_last_of("/\\") + 1);
            size_t pos = 0;
            while (pos < src.size()) {
                size_t eol = src.find('\n', pos);
                if (eol == std::string::npos) eol = src.size();
                const std::string line = src.substr(pos, eol - pos);
                const size_t at = line.find("#include");
                const size_t q0 = line.find('"'), q1 = line.rfind('"');
                if (at != std::string::npos && line.find_first_not_of(" \t") == at && q0 != std::string::npos && q1 > q0 && depth < 4) {
                    const std::string file = dir + line.substr(q0 + 1, q1 - q0 - 1);
                    inc.push_back(file);
                    out += Load(file, inc, depth + 1);
                } else {
                    out += line;
                }
                out += "\n";
                pos = eol + 1;
            }
            return out;
        }
    };
    return Rec::Load(path, includes, 0);
}

bool ComputeProgram::load(const std::string& p) {
    if (!Compute_Available()) return false;
    std::vector<std::string> inc;
    const std::string src = LoadShaderSource(p, inc);
    if (src.empty()) { logger.logf("[ERROR] compute shader %s: not found\n", p.c_str()); return false; }
    const unsigned int sh = rlCompileShader(src.c_str(), RL_COMPUTE_SHADER);
    if (sh == 0) { logger.logf("[ERROR] compute shader %s failed to compile\n", p.c_str()); TraceLog(LOG_WARNING, "Compute shader %s failed to compile", p.c_str()); return false; }
    const unsigned int prog = rlLoadComputeShaderProgram(sh);
    if (prog == 0) { logger.logf("[ERROR] compute shader %s failed to link\n", p.c_str()); return false; }
    unload();
    id = prog;
    path = p;
    return true;
}
void ComputeProgram::unload() { if (id) rlUnloadShaderProgram(id); id = 0; }
void ComputeProgram::use() const { rlEnableShader(id); }
int ComputeProgram::uniform(const char* name) const { return rlGetLocationUniform(id, name); }
void ComputeProgram::setInt(const char* name, int v) const { const int l = uniform(name); if (l >= 0) rlSetUniform(l, &v, RL_SHADER_UNIFORM_INT, 1); }
void ComputeProgram::setUint(const char* name, unsigned v) const { const int l = uniform(name); if (l >= 0) rlSetUniform(l, &v, RL_SHADER_UNIFORM_UINT, 1); }
void ComputeProgram::setFloat(const char* name, float v) const { const int l = uniform(name); if (l >= 0) rlSetUniform(l, &v, RL_SHADER_UNIFORM_FLOAT, 1); }
void ComputeProgram::setVec3(const char* name, float x, float y, float z) const { const float v[3] = { x, y, z }; const int l = uniform(name); if (l >= 0) rlSetUniform(l, v, RL_SHADER_UNIFORM_VEC3, 1); }
void ComputeProgram::setVec4(const char* name, float x, float y, float z, float w) const { const float v[4] = { x, y, z, w }; const int l = uniform(name); if (l >= 0) rlSetUniform(l, v, RL_SHADER_UNIFORM_VEC4, 1); }
void ComputeProgram::dispatch(unsigned gx, unsigned gy, unsigned gz) const {
    rlComputeShaderDispatch(gx, gy, gz);
    Compute_Barrier();
}

void ComputeBuffer::create(unsigned int size, const void* data, bool dynamic) {
    destroy();
    id = rlLoadShaderBuffer(size, data, dynamic ? RL_DYNAMIC_COPY : RL_STATIC_DRAW);
    bytes = size;
}
void ComputeBuffer::resize(unsigned int size) { create(size, nullptr, true); }
void ComputeBuffer::upload(const void* data, unsigned int size, unsigned int offset) { if (id) rlUpdateShaderBuffer(id, data, size, offset); }
void ComputeBuffer::download(void* dest, unsigned int size, unsigned int offset) const { if (id) rlReadShaderBuffer(id, dest, size, offset); }
void ComputeBuffer::bind(unsigned int index) const { rlBindShaderBuffer(id, index); }
void ComputeBuffer::destroy() { if (id) rlUnloadShaderBuffer(id); id = 0; bytes = 0; }

void Compute_Barrier() {
    // GL_SHADER_STORAGE_BARRIER_BIT 0x2000 | GL_BUFFER_UPDATE_BARRIER_BIT 0x200 | GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT 0x1
    // | GL_COMMAND_BARRIER_BIT 0x40 | GL_ELEMENT_ARRAY_BARRIER_BIT 0x2
    if (glad_glMemoryBarrier) glad_glMemoryBarrier(0x2000 | 0x200 | 0x1 | 0x40 | 0x2);
}
