#include "shaderasync.h"
#include "utils/logger.h"

#include "external/glad.h"
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <vector>

extern "C" void* glfwGetProcAddress(const char* name); // raylib links GLFW (the headless platform answers it too)

namespace {
constexpr GLenum kCompletionStatus = 0x91B1; // GL_COMPLETION_STATUS_KHR (= _ARB)
typedef void (GLAD_API_PTR* PfnMaxThreads)(GLuint count);

int gState = -1; // -1 not decided, 0 off, 1 on
}

bool ShaderAsync_Available() {
    if (gState >= 0) return gState == 1;
    gState = 0;
    const char* env = getenv("RAYLIBGAME_ASYNC_SHADERS");
    if (env && env[0] == '0') { LOG_INFO("Shaders", "async compute compiles off (RAYLIBGAME_ASYNC_SHADERS=0)\n"); return false; }
#if RAYLIBGAME_HEADLESS_PLATFORM
    // The headless Linux platform (EGL; the server and the CI runner's clients): NVIDIA's driver there reports the compile
    // complete but does the work when the status is first read (fluid_step.comp: 5-12 s waited after a 4-minute "background"
    // compile on the runner), so it gains nothing; compile as before unless RAYLIBGAME_ASYNC_SHADERS=1 asks for it.
    if (!(env && env[0] == '1')) { LOG_INFO("Shaders", "async compute compiles off on the headless platform\n"); return false; }
#endif
    if (!glad_glGetStringi || !glad_glCreateShader) return false;
    bool khr = false, arb = false;
    GLint n = 0;
    glGetIntegerv(GL_NUM_EXTENSIONS, &n);
    for (GLint i = 0; i < n; ++i)
        if (const char* e = (const char*)glGetStringi(GL_EXTENSIONS, (GLuint)i)) {
            khr = khr || strcmp(e, "GL_KHR_parallel_shader_compile") == 0;
            arb = arb || strcmp(e, "GL_ARB_parallel_shader_compile") == 0;
        }
    if (!khr && !arb) { LOG_INFO("Shaders", "no parallel shader compile: compute compiles block\n"); return false; }
    PfnMaxThreads maxThreads = (PfnMaxThreads)glfwGetProcAddress(khr ? "glMaxShaderCompilerThreadsKHR" : "glMaxShaderCompilerThreadsARB");
    if (maxThreads) maxThreads(0xFFFFFFFFu); // as many threads as the driver likes
    gState = 1;
    LOG_INFO("Shaders", "parallel shader compile (%s): compute compiles run in the background\n", khr ? "KHR" : "ARB");
    return true;
}

AsyncProgram ShaderAsync_BeginCompute(const std::string& src) {
    AsyncProgram p;
    if (!ShaderAsync_Available()) return p;
    p.shader = glCreateShader(GL_COMPUTE_SHADER);
    const char* text = src.c_str();
    glShaderSource(p.shader, 1, &text, nullptr);
    glCompileShader(p.shader);
    p.program = glCreateProgram();
    glAttachShader(p.program, p.shader);
    glLinkProgram(p.program); // queued behind the compile on the driver's threads; no status asked yet (that would wait)
    return p;
}

bool ShaderAsync_Done(const AsyncProgram& p) {
    if (!p.started()) return true;
    GLint done = 0;
    glGetProgramiv(p.program, kCompletionStatus, &done);
    return done != 0;
}

unsigned ShaderAsync_Finish(AsyncProgram& p, const char* name) {
    if (!p.started()) return 0;
    GLint ok = 0;
    glGetShaderiv(p.shader, GL_COMPILE_STATUS, &ok); // waits for whatever is left
    if (!ok) {
        GLint len = 0; glGetShaderiv(p.shader, GL_INFO_LOG_LENGTH, &len);
        std::vector<char> log((size_t)std::max(len, 1));
        glGetShaderInfoLog(p.shader, (GLsizei)log.size(), nullptr, log.data());
        LOG_ERROR("Shaders", "compute shader %s failed to compile:\n%s\n", name, log.data());
        ShaderAsync_Discard(p);
        return 0;
    }
    glGetProgramiv(p.program, GL_LINK_STATUS, &ok);
    if (!ok) {
        GLint len = 0; glGetProgramiv(p.program, GL_INFO_LOG_LENGTH, &len);
        std::vector<char> log((size_t)std::max(len, 1));
        glGetProgramInfoLog(p.program, (GLsizei)log.size(), nullptr, log.data());
        LOG_ERROR("Shaders", "compute shader %s failed to link:\n%s\n", name, log.data());
        ShaderAsync_Discard(p);
        return 0;
    }
    const unsigned prog = p.program;
    glDetachShader(prog, p.shader);
    glDeleteShader(p.shader);
    p = AsyncProgram{};
    return prog;
}

void ShaderAsync_Discard(AsyncProgram& p) {
    if (p.program) glDeleteProgram(p.program);
    if (p.shader) glDeleteShader(p.shader);
    p = AsyncProgram{};
}
