// A headless raylib platform backend for the Linux dedicated server (#792, under #326).
//
// raylib's rcore.c includes one platform backend (GLFW, SDL, DRM, ...) chosen by a PLATFORM_* define; with none defined
// it includes nothing and expects a custom one. This file is that custom backend, compiled as one unit with rcore.c
// (it includes rcore.c first, so it sees raylib's CORE state, as the bundled backends do). CMake builds raylib from
// this file instead of raylib's own CMake when RAYLIBGAME_HEADLESS is on (Linux; see CMakeLists.txt).
//
// No window system: there is no X11, Wayland or GLFW in the process. "The window" is an off-screen GL context:
//   1. EGL (libEGL.so.1, loaded with dlopen, so the binary needs no GL libraries to start):
//      a. every EGL device (EGL_EXT_platform_device: the NVIDIA driver in a container started with --gpus / the
//         nvidia runtime, or Mesa on a /dev/dri render node), then
//      b. Mesa's surfaceless platform, then the default display;
//      each asked for an OpenGL 4.3 core context with a small pbuffer (or no surface: EGL_KHR_surfaceless_context).
//      Software renderers (llvmpipe, softpipe) are skipped unless RAYLIBGAME_EGL_SOFTWARE=1: a server without a GPU
//      is cheaper on the null GL below than emulating the GPU paths on its CPU.
//   2. Null GL, when no hardware context can be made (no GPU in the container): every GL function is a stub that
//      does nothing and returns a plausible value (ids, "compiled", "complete", zeroed reads). The game sees
//      RaylibHeadless_GLMode() == 2 and keeps its GPU paths off (Compute_Available() is false), so the sim runs on the
//      CPU paths: CPU fluid, CPU terrain density. Nothing is drawn on a server anyway.
//   RAYLIBGAME_NO_GL=1 forces the null GL; RAYLIBGAME_EGL_DEVICE=N picks one EGL device; RAYLIBGAME_EGL_LIB names the
//   libEGL to load. Input is empty (a server has a stdin console instead, src/space/serverconsole.*).

#define _GNU_SOURCE 1  // setenv, readlink (rcore.c asks only for _XOPEN_SOURCE 500)
#include "rcore.c"

#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

//----------------------------------------------------------------------------------
// EGL, declared here (no EGL headers needed to build; the library is loaded at run time)
//----------------------------------------------------------------------------------
typedef void *EGLDisplay_, *EGLConfig_, *EGLContext_, *EGLSurface_, *EGLDeviceEXT_;
typedef int32_t EGLint_;
typedef unsigned int EGLBoolean_, EGLenum_;

enum {
    kEGL_NONE = 0x3038, kEGL_SURFACE_TYPE = 0x3033, kEGL_PBUFFER_BIT = 0x0001, kEGL_RENDERABLE_TYPE = 0x3040,
    kEGL_OPENGL_BIT = 0x0008, kEGL_RED_SIZE = 0x3024, kEGL_GREEN_SIZE = 0x3023, kEGL_BLUE_SIZE = 0x3022,
    kEGL_ALPHA_SIZE = 0x3021, kEGL_DEPTH_SIZE = 0x3025, kEGL_WIDTH = 0x3057, kEGL_HEIGHT = 0x3056,
    kEGL_OPENGL_API = 0x30A2, kEGL_CONTEXT_MAJOR_VERSION = 0x3098, kEGL_CONTEXT_MINOR_VERSION = 0x30FB,
    kEGL_CONTEXT_OPENGL_PROFILE_MASK = 0x30FD, kEGL_CONTEXT_OPENGL_CORE_PROFILE_BIT = 0x0001, kEGL_EXTENSIONS = 0x3055,
    kEGL_VENDOR = 0x3053, kEGL_VERSION = 0x3054, kEGL_PLATFORM_DEVICE_EXT = 0x313F,
    kEGL_PLATFORM_SURFACELESS_MESA = 0x31DD, kEGL_DRM_DEVICE_FILE_EXT = 0x3233,
};

static struct {
    void *lib;
    void *(*GetProcAddress)(const char *);
    EGLint_ (*GetError)(void);
    EGLDisplay_ (*GetDisplay)(void *);
    EGLBoolean_ (*Initialize)(EGLDisplay_, EGLint_ *, EGLint_ *);
    EGLBoolean_ (*Terminate)(EGLDisplay_);
    const char *(*QueryString)(EGLDisplay_, EGLint_);
    EGLBoolean_ (*ChooseConfig)(EGLDisplay_, const EGLint_ *, EGLConfig_ *, EGLint_, EGLint_ *);
    EGLBoolean_ (*BindAPI)(EGLenum_);
    EGLContext_ (*CreateContext)(EGLDisplay_, EGLConfig_, EGLContext_, const EGLint_ *);
    EGLBoolean_ (*DestroyContext)(EGLDisplay_, EGLContext_);
    EGLSurface_ (*CreatePbufferSurface)(EGLDisplay_, EGLConfig_, const EGLint_ *);
    EGLBoolean_ (*DestroySurface)(EGLDisplay_, EGLSurface_);
    EGLBoolean_ (*MakeCurrent)(EGLDisplay_, EGLSurface_, EGLSurface_, EGLContext_);
    EGLBoolean_ (*SwapBuffers)(EGLDisplay_, EGLSurface_);
    // extensions
    EGLBoolean_ (*QueryDevicesEXT)(EGLint_, EGLDeviceEXT_ *, EGLint_ *);
    const char *(*QueryDeviceStringEXT)(EGLDeviceEXT_, EGLint_);
    EGLDisplay_ (*GetPlatformDisplayEXT)(EGLenum_, void *, const EGLint_ *);
} egl;

typedef struct {
    int glMode;                 // 0 not initialised, 1 EGL context, 2 null GL
    EGLDisplay_ display;
    EGLContext_ context;
    EGLSurface_ surface;
    char description[256];      // what was picked, for the log
} PlatformData;

extern CoreData CORE;
static PlatformData platform = { 0 };

int InitPlatform(void);
void ClosePlatform(void);

// 0 before InitWindow, 1 = a real GL context from EGL, 2 = the null GL (no GPU): the game turns its GPU paths off.
int RaylibHeadless_GLMode(void) { return platform.glMode; }
const char *RaylibHeadless_Description(void) { return platform.description; }

//----------------------------------------------------------------------------------
// Null GL: stub functions for every GL entry point (see the top of the file)
//----------------------------------------------------------------------------------
// GL enums used by the stubs (from glcorearb.h)
enum {
    kGL_VENDOR = 0x1F00, kGL_RENDERER = 0x1F01, kGL_VERSION = 0x1F02, kGL_EXTENSIONS = 0x1F03,
    kGL_SHADING_LANGUAGE_VERSION = 0x8B8C, kGL_VIEWPORT = 0x0BA2, kGL_SCISSOR_BOX = 0x0C10,
    kGL_MAX_VIEWPORT_DIMS = 0x0D3A, kGL_DEPTH_RANGE = 0x0B70, kGL_COLOR_CLEAR_VALUE = 0x0C22,
    kGL_COLOR_WRITEMASK = 0x0C23, kGL_MAX_TEXTURE_SIZE = 0x0D33, kGL_MAX_SAMPLES = 0x8D57,
    kGL_MAX_COLOR_ATTACHMENTS = 0x8CDF, kGL_MAX_DRAW_BUFFERS = 0x8824, kGL_MAX_TEXTURE_IMAGE_UNITS = 0x8872,
    kGL_MAX_COMBINED_TEXTURE_IMAGE_UNITS = 0x8B4D, kGL_MAX_VERTEX_ATTRIBS = 0x8869,
    kGL_MAX_SHADER_STORAGE_BUFFER_BINDINGS = 0x90DD, kGL_MAX_SHADER_STORAGE_BLOCK_SIZE = 0x90DE,
    kGL_MAX_UNIFORM_BLOCK_SIZE = 0x8A30, kGL_MAJOR_VERSION = 0x821B, kGL_MINOR_VERSION = 0x821C,
    kGL_NUM_EXTENSIONS = 0x821D, kGL_COMPRESSED_TEXTURE_FORMATS = 0x86A3, kGL_MAX_TEXTURE_MAX_ANISOTROPY = 0x84FF,
    kGL_MAX_COMPUTE_WORK_GROUP_COUNT = 0x91BE, kGL_MAX_COMPUTE_WORK_GROUP_SIZE = 0x91BF,
    kGL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS = 0x90EB, kGL_MAX_ARRAY_TEXTURE_LAYERS = 0x88FF,
    kGL_MAX_3D_TEXTURE_SIZE = 0x8073, kGL_MAX_CUBE_MAP_TEXTURE_SIZE = 0x851C, kGL_MAX_RENDERBUFFER_SIZE = 0x84E8,
    kGL_COMPILE_STATUS = 0x8B81, kGL_LINK_STATUS = 0x8B82, kGL_VALIDATE_STATUS = 0x8B83,
    kGL_COMPUTE_WORK_GROUP_SIZE = 0x8267, kGL_QUERY_RESULT_AVAILABLE = 0x8867,
    kGL_FRAMEBUFFER_COMPLETE = 0x8CD5, kGL_ALREADY_SIGNALED = 0x911A,
    kGL_RGBA = 0x1908, kGL_RGB = 0x1907, kGL_RG = 0x8227, kGL_BGRA = 0x80E1,
    kGL_UNSIGNED_BYTE = 0x1401, kGL_FLOAT = 0x1406, kGL_UNSIGNED_SHORT = 0x1403, kGL_HALF_FLOAT = 0x140B,
};

static intptr_t NullGL_Zero(void) { return 0; } // any GL function without a useful result (the x86-64 SysV ABI lets a
                                                // caller pass arguments a callee ignores)
static unsigned int nullNextId = 1;
static unsigned int NullGL_NewId(void) { return __atomic_fetch_add(&nullNextId, 1, __ATOMIC_RELAXED); }
static void NullGL_Gen(int n, unsigned int *ids) { for (int i = 0; ids && i < n; i++) ids[i] = NullGL_NewId(); }
static void NullGL_CreateTyped(unsigned int target, int n, unsigned int *ids) { (void)target; NullGL_Gen(n, ids); }
static unsigned int NullGL_CreateObj(void) { return NullGL_NewId(); }
static unsigned int NullGL_CreateShaderProgramv(unsigned int t, int n, const char *const *s) { (void)t; (void)n; (void)s; return NullGL_NewId(); }
static unsigned int NullGL_GetError(void) { return 0; }
static const unsigned char *NullGL_GetString(unsigned int name)
{
    switch (name)
    {
        case kGL_VENDOR: return (const unsigned char *)"RaylibGame";
        case kGL_RENDERER: return (const unsigned char *)"null GL (no GPU: nothing is drawn, GPU paths off)";
        case kGL_VERSION: return (const unsigned char *)"4.3 null";
        case kGL_SHADING_LANGUAGE_VERSION: return (const unsigned char *)"4.30";
        default: return (const unsigned char *)"";
    }
}
static const unsigned char *NullGL_GetStringi(unsigned int name, unsigned int i) { (void)name; (void)i; return (const unsigned char *)""; }
static int NullGL_IntValue(unsigned int pname, int *out) // how many values pname writes; fills them
{
    const int w = CORE.Window.screen.width, h = CORE.Window.screen.height;
    switch (pname)
    {
        case kGL_VIEWPORT: case kGL_SCISSOR_BOX: out[0] = 0; out[1] = 0; out[2] = w; out[3] = h; return 4;
        case kGL_MAX_VIEWPORT_DIMS: out[0] = out[1] = 16384; return 2;
        case kGL_DEPTH_RANGE: out[0] = 0; out[1] = 1; return 2;
        case kGL_COLOR_CLEAR_VALUE: case kGL_COLOR_WRITEMASK: out[0] = out[1] = out[2] = out[3] = 0; return 4;
        case kGL_COMPRESSED_TEXTURE_FORMATS: return 0; // GL_NUM_COMPRESSED_TEXTURE_FORMATS said 0
        case kGL_MAX_TEXTURE_SIZE: case kGL_MAX_RENDERBUFFER_SIZE: case kGL_MAX_CUBE_MAP_TEXTURE_SIZE: out[0] = 16384; return 1;
        case kGL_MAX_3D_TEXTURE_SIZE: case kGL_MAX_ARRAY_TEXTURE_LAYERS: out[0] = 2048; return 1;
        case kGL_MAX_SAMPLES: case kGL_MAX_COLOR_ATTACHMENTS: case kGL_MAX_DRAW_BUFFERS: out[0] = 8; return 1;
        case kGL_MAX_TEXTURE_IMAGE_UNITS: out[0] = 32; return 1;
        case kGL_MAX_COMBINED_TEXTURE_IMAGE_UNITS: out[0] = 192; return 1;
        case kGL_MAX_VERTEX_ATTRIBS: case kGL_MAX_SHADER_STORAGE_BUFFER_BINDINGS: out[0] = 16; return 1;
        case kGL_MAX_SHADER_STORAGE_BLOCK_SIZE: out[0] = 1 << 27; return 1;
        case kGL_MAX_UNIFORM_BLOCK_SIZE: out[0] = 65536; return 1;
        case kGL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS: out[0] = 1024; return 1;
        case kGL_MAJOR_VERSION: out[0] = 4; return 1;
        case kGL_MINOR_VERSION: out[0] = 3; return 1;
        default: out[0] = 0; return 1; // GL_NUM_EXTENSIONS and the rest: 0
    }
}
static void NullGL_GetIntegerv(unsigned int pname, int *data) { int v[4]; const int n = NullGL_IntValue(pname, v); for (int i = 0; data && i < n; i++) data[i] = v[i]; }
static void NullGL_GetInteger64v(unsigned int pname, int64_t *data) { int v[4]; const int n = NullGL_IntValue(pname, v); for (int i = 0; data && i < n; i++) data[i] = v[i]; }
static void NullGL_GetFloatv(unsigned int pname, float *data)
{
    if (!data) return;
    if (pname == kGL_MAX_TEXTURE_MAX_ANISOTROPY) { data[0] = 1.0f; return; }
    int v[4]; const int n = NullGL_IntValue(pname, v); for (int i = 0; i < n; i++) data[i] = (float)v[i];
}
static void NullGL_GetDoublev(unsigned int pname, double *data) { int v[4]; const int n = NullGL_IntValue(pname, v); for (int i = 0; data && i < n; i++) data[i] = v[i]; }
static void NullGL_GetBooleanv(unsigned int pname, unsigned char *data) { int v[4]; const int n = NullGL_IntValue(pname, v); for (int i = 0; data && i < n; i++) data[i] = (unsigned char)(v[i] != 0); }
static void NullGL_GetIntegeri_v(unsigned int target, unsigned int index, int *data)
{
    (void)index; if (!data) return;
    data[0] = target == kGL_MAX_COMPUTE_WORK_GROUP_COUNT ? 65535 : target == kGL_MAX_COMPUTE_WORK_GROUP_SIZE ? 1024 : 0;
}
static void NullGL_GetShaderiv(unsigned int obj, unsigned int pname, int *params) { (void)obj; if (params) params[0] = pname == kGL_COMPILE_STATUS ? 1 : 0; }
static void NullGL_GetProgramiv(unsigned int obj, unsigned int pname, int *params)
{
    (void)obj; if (!params) return;
    if (pname == kGL_COMPUTE_WORK_GROUP_SIZE) { params[0] = params[1] = params[2] = 1; return; }
    params[0] = (pname == kGL_LINK_STATUS || pname == kGL_VALIDATE_STATUS) ? 1 : 0;
}
static void NullGL_GetInfoLog(unsigned int obj, int bufSize, int *length, char *log) { (void)obj; if (length) *length = 0; if (log && bufSize > 0) log[0] = 0; }
static int NullGL_GetLocation(unsigned int program, const char *name) { (void)program; (void)name; return 0; }
static unsigned int NullGL_CheckFramebufferStatus(unsigned int target) { (void)target; return kGL_FRAMEBUFFER_COMPLETE; }
static unsigned int NullGL_CheckNamedFramebufferStatus(unsigned int fb, unsigned int target) { (void)fb; (void)target; return kGL_FRAMEBUFFER_COMPLETE; }
static void *NullGL_FenceSync(unsigned int cond, unsigned int flags) { (void)cond; (void)flags; return (void *)1; }
static unsigned int NullGL_ClientWaitSync(void *sync, unsigned int flags, uint64_t timeout) { (void)sync; (void)flags; (void)timeout; return kGL_ALREADY_SIGNALED; }
static unsigned char NullGL_True(void) { return 1; }
static void NullGL_GetQueryObjectiv(unsigned int id, unsigned int pname, int *p) { (void)id; if (p) p[0] = pname == kGL_QUERY_RESULT_AVAILABLE ? 1 : 0; }
static void NullGL_GetQueryObjecti64v(unsigned int id, unsigned int pname, int64_t *p) { (void)id; if (p) p[0] = pname == kGL_QUERY_RESULT_AVAILABLE ? 1 : 0; }
static void NullGL_GetOneInt(unsigned int a, unsigned int b, int *p) { (void)a; (void)b; if (p) p[0] = 0; }
static void NullGL_GetOneInt4(unsigned int a, unsigned int b, unsigned int c, int *p) { (void)a; (void)b; (void)c; if (p) p[0] = 0; }
static void NullGL_GetTexLevelParameteriv(unsigned int t, int level, unsigned int pname, int *p) { (void)t; (void)level; (void)pname; if (p) p[0] = 0; }
static void NullGL_GetBufferSubData(unsigned int target, intptr_t offset, intptr_t size, void *data) { (void)target; (void)offset; if (data && size > 0) memset(data, 0, (size_t)size); }
static void NullGL_ReadPixels(int x, int y, int w, int h, unsigned int format, unsigned int type, void *data)
{
    (void)x; (void)y;
    if ((uintptr_t)data < 65536 || w <= 0 || h <= 0) return; // a pixel-pack buffer offset, or nothing to write
    const int comps = (format == kGL_RGBA || format == kGL_BGRA) ? 4 : format == kGL_RGB ? 3 : format == kGL_RG ? 2 : 1;
    const int size = (type == kGL_FLOAT) ? 4 : (type == kGL_UNSIGNED_SHORT || type == kGL_HALF_FLOAT) ? 2 : (type == kGL_UNSIGNED_BYTE) ? 1 : 4;
    memset(data, 0, (size_t)w*h*comps*size);
}
// Mapped buffers: zeroed memory, freed on unmap (persistent maps live as long as the process).
#define NULLGL_MAPS 64
static struct { unsigned int key; void *ptr; } nullMaps[NULLGL_MAPS];
static void *NullGL_MapAlloc(unsigned int key, intptr_t length)
{
    if (length <= 0) length = 1 << 20; // glMapBuffer: the size is not known; 1 MB of zeros
    void *p = calloc(1, (size_t)length);
    for (int i = 0; i < NULLGL_MAPS; i++) if (!nullMaps[i].ptr) { nullMaps[i].key = key; nullMaps[i].ptr = p; break; }
    return p;
}
static unsigned char NullGL_MapFree(unsigned int key)
{
    for (int i = 0; i < NULLGL_MAPS; i++) if (nullMaps[i].ptr && nullMaps[i].key == key) { free(nullMaps[i].ptr); nullMaps[i].ptr = NULL; break; }
    return 1;
}
static void *NullGL_MapBuffer(unsigned int target, unsigned int access) { (void)access; return NullGL_MapAlloc(target, 0); }
static void *NullGL_MapBufferRange(unsigned int target, intptr_t offset, intptr_t length, unsigned int access) { (void)offset; (void)access; return NullGL_MapAlloc(target, length); }
static unsigned char NullGL_UnmapBuffer(unsigned int target) { return NullGL_MapFree(target); }
// (the DSA variants key on the buffer name instead of the target; ids and targets do not collide in practice)
static void *NullGL_MapNamedBufferRange(unsigned int buffer, intptr_t offset, intptr_t length, unsigned int access) { (void)offset; (void)access; return NullGL_MapAlloc(0x80000000u | buffer, length); }
static unsigned char NullGL_UnmapNamedBuffer(unsigned int buffer) { return NullGL_MapFree(0x80000000u | buffer); }

static void *NullGL_GetProcAddress(const char *name)
{
    static const struct { const char *name; void *fn; } kTable[] = {
        { "glGetError", (void *)NullGL_GetError },
        { "glGetString", (void *)NullGL_GetString }, { "glGetStringi", (void *)NullGL_GetStringi },
        { "glGetIntegerv", (void *)NullGL_GetIntegerv }, { "glGetInteger64v", (void *)NullGL_GetInteger64v },
        { "glGetFloatv", (void *)NullGL_GetFloatv }, { "glGetDoublev", (void *)NullGL_GetDoublev },
        { "glGetBooleanv", (void *)NullGL_GetBooleanv }, { "glGetIntegeri_v", (void *)NullGL_GetIntegeri_v },
        { "glGetShaderiv", (void *)NullGL_GetShaderiv }, { "glGetProgramiv", (void *)NullGL_GetProgramiv },
        { "glGetShaderInfoLog", (void *)NullGL_GetInfoLog }, { "glGetProgramInfoLog", (void *)NullGL_GetInfoLog },
        { "glGetUniformLocation", (void *)NullGL_GetLocation }, { "glGetAttribLocation", (void *)NullGL_GetLocation },
        { "glGetUniformBlockIndex", (void *)NullGL_GetLocation },
        { "glCreateShader", (void *)NullGL_CreateObj }, { "glCreateProgram", (void *)NullGL_CreateObj },
        { "glCreateShaderProgramv", (void *)NullGL_CreateShaderProgramv },
        { "glGenBuffers", (void *)NullGL_Gen }, { "glGenTextures", (void *)NullGL_Gen },
        { "glGenVertexArrays", (void *)NullGL_Gen }, { "glGenFramebuffers", (void *)NullGL_Gen },
        { "glGenRenderbuffers", (void *)NullGL_Gen }, { "glGenQueries", (void *)NullGL_Gen },
        { "glGenSamplers", (void *)NullGL_Gen }, { "glGenTransformFeedbacks", (void *)NullGL_Gen },
        { "glGenProgramPipelines", (void *)NullGL_Gen },
        { "glCreateBuffers", (void *)NullGL_Gen }, { "glCreateVertexArrays", (void *)NullGL_Gen },
        { "glCreateFramebuffers", (void *)NullGL_Gen }, { "glCreateRenderbuffers", (void *)NullGL_Gen },
        { "glCreateSamplers", (void *)NullGL_Gen },
        { "glCreateTextures", (void *)NullGL_CreateTyped }, { "glCreateQueries", (void *)NullGL_CreateTyped },
        { "glCheckFramebufferStatus", (void *)NullGL_CheckFramebufferStatus },
        { "glCheckNamedFramebufferStatus", (void *)NullGL_CheckNamedFramebufferStatus },
        { "glFenceSync", (void *)NullGL_FenceSync }, { "glClientWaitSync", (void *)NullGL_ClientWaitSync },
        { "glIsSync", (void *)NullGL_True }, { "glIsShader", (void *)NullGL_True }, { "glIsProgram", (void *)NullGL_True },
        { "glIsBuffer", (void *)NullGL_True }, { "glIsTexture", (void *)NullGL_True },
        { "glGetQueryObjectiv", (void *)NullGL_GetQueryObjectiv }, { "glGetQueryObjectuiv", (void *)NullGL_GetQueryObjectiv },
        { "glGetQueryObjecti64v", (void *)NullGL_GetQueryObjecti64v }, { "glGetQueryObjectui64v", (void *)NullGL_GetQueryObjecti64v },
        { "glGetQueryiv", (void *)NullGL_GetOneInt }, { "glGetBufferParameteriv", (void *)NullGL_GetOneInt },
        { "glGetTexParameteriv", (void *)NullGL_GetOneInt }, { "glGetRenderbufferParameteriv", (void *)NullGL_GetOneInt },
        { "glGetNamedBufferParameteriv", (void *)NullGL_GetOneInt }, { "glGetFramebufferParameteriv", (void *)NullGL_GetOneInt },
        { "glGetFramebufferAttachmentParameteriv", (void *)NullGL_GetOneInt4 }, { "glGetActiveUniformBlockiv", (void *)NullGL_GetOneInt4 },
        { "glGetProgramInterfaceiv", (void *)NullGL_GetOneInt4 }, { "glGetQueryIndexediv", (void *)NullGL_GetOneInt4 },
        { "glGetProgramStageiv", (void *)NullGL_GetOneInt4 },
        { "glGetTexLevelParameteriv", (void *)NullGL_GetTexLevelParameteriv },
        { "glGetBufferSubData", (void *)NullGL_GetBufferSubData }, { "glGetNamedBufferSubData", (void *)NullGL_GetBufferSubData },
        { "glReadPixels", (void *)NullGL_ReadPixels },
        { "glMapBuffer", (void *)NullGL_MapBuffer }, { "glMapBufferRange", (void *)NullGL_MapBufferRange },
        { "glUnmapBuffer", (void *)NullGL_UnmapBuffer },
        { "glMapNamedBufferRange", (void *)NullGL_MapNamedBufferRange }, { "glUnmapNamedBuffer", (void *)NullGL_UnmapNamedBuffer },
    };
    for (size_t i = 0; i < sizeof(kTable)/sizeof(kTable[0]); i++) if (!strcmp(kTable[i].name, name)) return kTable[i].fn;
    return (void *)NullGL_Zero;
}

//----------------------------------------------------------------------------------
// EGL context
//----------------------------------------------------------------------------------
static void *HeadlessGetProcAddress(const char *name)
{
    if (platform.glMode == 2) return NullGL_GetProcAddress(name);
    return (egl.GetProcAddress != NULL)? egl.GetProcAddress(name) : NULL;
}

// The game's engine resolves a few optional GL entry points by name (glClipControl, glMultiDrawElementsIndirectCount)
// through glfwGetProcAddress, which the GLFW backend provides. There is no GLFW here: answer with the same loader.
void *glfwGetProcAddress(const char *name) { return HeadlessGetProcAddress(name); }

static int EnvIs(const char *name, const char *value) { const char *e = getenv(name); return (e != NULL) && !strcmp(e, value); }

static int LoadEGL(void)
{
    char exeDir[1024] = { 0 };
    ssize_t n = readlink("/proc/self/exe", exeDir, sizeof(exeDir) - 1);
    if (n > 0) { exeDir[n] = 0; char *slash = strrchr(exeDir, '/'); if (slash) *slash = 0; } else exeDir[0] = 0;

    // A runtime folder may bring libglvnd's dispatcher (lib/libEGL.so.1 + lib/libGLdispatch.so.0), for images without it.
    // The vendor library itself (libEGL_nvidia.so.0) comes from the driver, mounted in by the NVIDIA container runtime.
    char path[1200];
    const char *candidates[4] = { getenv("RAYLIBGAME_EGL_LIB"), NULL, "libEGL.so.1", "libEGL.so" };
    char bundled[1200] = { 0 };
    if (exeDir[0])
    {
        snprintf(bundled, sizeof(bundled), "%s/lib/libEGL.so.1", exeDir);
        if (access(bundled, R_OK) == 0)
        {
            snprintf(path, sizeof(path), "%s/lib/libGLdispatch.so.0", exeDir);
            if (access(path, R_OK) == 0) dlopen(path, RTLD_NOW | RTLD_GLOBAL); // so the bundled libEGL finds it by soname
            candidates[1] = bundled;
        }
    }
    // Prefer the system's libEGL (it matches the system's glvnd and vendor files); the bundled one is the fallback.
    const char *order[4] = { candidates[0], candidates[2], candidates[3], candidates[1] };
    for (int i = 0; i < 4 && !egl.lib; i++) if (order[i]) egl.lib = dlopen(order[i], RTLD_NOW | RTLD_LOCAL);
    if (!egl.lib) return 0;

    // glvnd finds vendor libraries through JSON files; the runtime folder ships them for images that have none.
    if (exeDir[0] && !getenv("__EGL_VENDOR_LIBRARY_FILENAMES") && !getenv("__EGL_VENDOR_LIBRARY_DIRS")
        && access("/usr/share/glvnd/egl_vendor.d/10_nvidia.json", R_OK) != 0 && access("/etc/glvnd/egl_vendor.d/10_nvidia.json", R_OK) != 0)
    {
        snprintf(path, sizeof(path), "/etc/glvnd/egl_vendor.d:/usr/share/glvnd/egl_vendor.d:%s/lib/egl_vendor.d", exeDir);
        setenv("__EGL_VENDOR_LIBRARY_DIRS", path, 0);
    }

#define LOAD(field, sym) egl.field = dlsym(egl.lib, sym); if (!egl.field) return 0
    LOAD(GetProcAddress, "eglGetProcAddress");
    LOAD(GetError, "eglGetError");
    LOAD(GetDisplay, "eglGetDisplay");
    LOAD(Initialize, "eglInitialize");
    LOAD(Terminate, "eglTerminate");
    LOAD(QueryString, "eglQueryString");
    LOAD(ChooseConfig, "eglChooseConfig");
    LOAD(BindAPI, "eglBindAPI");
    LOAD(CreateContext, "eglCreateContext");
    LOAD(DestroyContext, "eglDestroyContext");
    LOAD(CreatePbufferSurface, "eglCreatePbufferSurface");
    LOAD(DestroySurface, "eglDestroySurface");
    LOAD(MakeCurrent, "eglMakeCurrent");
    LOAD(SwapBuffers, "eglSwapBuffers");
#undef LOAD
    egl.QueryDevicesEXT = egl.GetProcAddress("eglQueryDevicesEXT");
    egl.QueryDeviceStringEXT = egl.GetProcAddress("eglQueryDeviceStringEXT");
    egl.GetPlatformDisplayEXT = egl.GetProcAddress("eglGetPlatformDisplayEXT");
    return 1;
}

static int HasWord(const char *list, const char *word)
{
    if (!list) return 0;
    const size_t n = strlen(word);
    for (const char *p = list; (p = strstr(p, word)) != NULL; p += n)
        if ((p == list || p[-1] == ' ') && (p[n] == ' ' || p[n] == 0)) return 1;
    return 0;
}

// Tries one display: an OpenGL 4.3 core context, current on a pbuffer (or on no surface). 1 = done.
static int TryDisplay(EGLDisplay_ dpy, const char *what)
{
    if (!dpy) return 0;
    EGLint_ major = 0, minor = 0;
    if (!egl.Initialize(dpy, &major, &minor)) { TRACELOG(LOG_INFO, "HEADLESS: %s: eglInitialize failed (0x%x)", what, egl.GetError()); return 0; }
    const char *exts = egl.QueryString(dpy, kEGL_EXTENSIONS);
    if (!egl.BindAPI(kEGL_OPENGL_API)) { TRACELOG(LOG_INFO, "HEADLESS: %s: no desktop OpenGL API", what); egl.Terminate(dpy); return 0; }

    const EGLint_ pbufferAttribs[] = { kEGL_SURFACE_TYPE, kEGL_PBUFFER_BIT, kEGL_RENDERABLE_TYPE, kEGL_OPENGL_BIT,
        kEGL_RED_SIZE, 8, kEGL_GREEN_SIZE, 8, kEGL_BLUE_SIZE, 8, kEGL_ALPHA_SIZE, 8, kEGL_DEPTH_SIZE, 24, kEGL_NONE };
    const EGLint_ anyAttribs[] = { kEGL_SURFACE_TYPE, 0, kEGL_RENDERABLE_TYPE, kEGL_OPENGL_BIT, kEGL_NONE };
    EGLConfig_ config = NULL; EGLint_ count = 0;
    int pbuffer = egl.ChooseConfig(dpy, pbufferAttribs, &config, 1, &count) && count > 0;
    if (!pbuffer)
    {
        count = 0;
        if (!HasWord(exts, "EGL_KHR_surfaceless_context") || !egl.ChooseConfig(dpy, anyAttribs, &config, 1, &count) || count == 0)
        { TRACELOG(LOG_INFO, "HEADLESS: %s: no OpenGL config", what); egl.Terminate(dpy); return 0; }
    }
    const EGLint_ ctxAttribs[] = { kEGL_CONTEXT_MAJOR_VERSION, 4, kEGL_CONTEXT_MINOR_VERSION, 3,
        kEGL_CONTEXT_OPENGL_PROFILE_MASK, kEGL_CONTEXT_OPENGL_CORE_PROFILE_BIT, kEGL_NONE };
    EGLContext_ ctx = egl.CreateContext(dpy, config, NULL, ctxAttribs);
    if (!ctx) { TRACELOG(LOG_INFO, "HEADLESS: %s: no OpenGL 4.3 core context (0x%x)", what, egl.GetError()); egl.Terminate(dpy); return 0; }
    EGLSurface_ surface = NULL;
    if (pbuffer)
    {
        const EGLint_ sizeAttribs[] = { kEGL_WIDTH, CORE.Window.screen.width, kEGL_HEIGHT, CORE.Window.screen.height, kEGL_NONE };
        surface = egl.CreatePbufferSurface(dpy, config, sizeAttribs);
    }
    if (!egl.MakeCurrent(dpy, surface, surface, ctx))
    {
        TRACELOG(LOG_INFO, "HEADLESS: %s: eglMakeCurrent failed (0x%x)", what, egl.GetError());
        if (surface) egl.DestroySurface(dpy, surface);
        egl.DestroyContext(dpy, ctx); egl.Terminate(dpy); return 0;
    }
    // Which renderer is it? glGetString is safe now that a context is current.
    const unsigned char *(*getString)(unsigned int) = (const unsigned char *(*)(unsigned int))egl.GetProcAddress("glGetString");
    const char *renderer = getString ? (const char *)getString(kGL_RENDERER) : NULL;
    const char *version = getString ? (const char *)getString(kGL_VERSION) : NULL;
    if (!renderer || !version)
    {
        TRACELOG(LOG_INFO, "HEADLESS: %s: the context answers no glGetString", what);
        egl.MakeCurrent(dpy, NULL, NULL, NULL);
        if (surface) egl.DestroySurface(dpy, surface);
        egl.DestroyContext(dpy, ctx); egl.Terminate(dpy); return 0;
    }
    const int software = strstr(renderer, "llvmpipe") || strstr(renderer, "softpipe") || strstr(renderer, "SwiftShader") || strstr(renderer, "Software Rasterizer");
    if (software && !EnvIs("RAYLIBGAME_EGL_SOFTWARE", "1"))
    {
        TRACELOG(LOG_INFO, "HEADLESS: %s: skipping the software renderer %s (RAYLIBGAME_EGL_SOFTWARE=1 allows it)", what, renderer);
        egl.MakeCurrent(dpy, NULL, NULL, NULL);
        if (surface) egl.DestroySurface(dpy, surface);
        egl.DestroyContext(dpy, ctx); egl.Terminate(dpy); return 0;
    }
    platform.display = dpy; platform.context = ctx; platform.surface = surface;
    snprintf(platform.description, sizeof(platform.description), "EGL %d.%d %s: %s, OpenGL %s%s", (int)major, (int)minor, what,
             renderer, version, surface ? "" : " (surfaceless)");
    return 1;
}

static int InitEGLContext(void)
{
    if (EnvIs("RAYLIBGAME_NO_GL", "1")) { TRACELOG(LOG_INFO, "HEADLESS: RAYLIBGAME_NO_GL=1"); return 0; }
    if (!LoadEGL()) { TRACELOG(LOG_INFO, "HEADLESS: no libEGL.so.1 (install libegl1 for GPU access)"); return 0; }
    const char *clientExts = egl.QueryString(NULL, kEGL_EXTENSIONS); // EGL_EXT_client_extensions: may be NULL
    const char *pick = getenv("RAYLIBGAME_EGL_DEVICE");

    // 1. EGL devices (a GPU without a display server)
    if (egl.QueryDevicesEXT && egl.GetPlatformDisplayEXT && HasWord(clientExts, "EGL_EXT_platform_device"))
    {
        EGLDeviceEXT_ devices[16]; EGLint_ count = 0;
        if (egl.QueryDevicesEXT(16, devices, &count))
        {
            TRACELOG(LOG_INFO, "HEADLESS: %d EGL device(s)", (int)count);
            for (int i = 0; i < count; i++)
            {
                if (pick && atoi(pick) != i) continue;
                const char *file = egl.QueryDeviceStringEXT ? egl.QueryDeviceStringEXT(devices[i], kEGL_DRM_DEVICE_FILE_EXT) : NULL;
                char what[160]; snprintf(what, sizeof(what), "device %d (%s)", i, file ? file : "no DRM file");
                if (TryDisplay(egl.GetPlatformDisplayEXT(kEGL_PLATFORM_DEVICE_EXT, devices[i], NULL), what)) return 1;
            }
        }
    }
    // 2. Mesa's surfaceless platform, then the default display
    if (egl.GetPlatformDisplayEXT && HasWord(clientExts, "EGL_MESA_platform_surfaceless")
        && TryDisplay(egl.GetPlatformDisplayEXT(kEGL_PLATFORM_SURFACELESS_MESA, NULL, NULL), "surfaceless")) return 1;
    if (TryDisplay(egl.GetDisplay(NULL), "default display")) return 1;
    return 0;
}

//----------------------------------------------------------------------------------
// Module Functions Definition: Window and Graphics Device (nothing to show; state only)
//----------------------------------------------------------------------------------
bool WindowShouldClose(void) { return CORE.Window.ready ? CORE.Window.shouldClose : true; }
void ToggleFullscreen(void) { }
void ToggleBorderlessWindowed(void) { }
void MaximizeWindow(void) { }
void MinimizeWindow(void) { }
void RestoreWindow(void) { }
void SetWindowState(unsigned int flags) { CORE.Window.flags |= flags; }
void ClearWindowState(unsigned int flags) { CORE.Window.flags &= ~flags; }
void SetWindowIcon(Image image) { (void)image; }
void SetWindowIcons(Image *images, int count) { (void)images; (void)count; }
void SetWindowTitle(const char *title) { CORE.Window.title = title; }
void SetWindowPosition(int x, int y) { CORE.Window.position = (Point){ x, y }; }
void SetWindowMonitor(int monitor) { (void)monitor; }
void SetWindowMinSize(int width, int height) { CORE.Window.screenMin.width = width; CORE.Window.screenMin.height = height; }
void SetWindowMaxSize(int width, int height) { CORE.Window.screenMax.width = width; CORE.Window.screenMax.height = height; }
void SetWindowSize(int width, int height) { CORE.Window.screen.width = width; CORE.Window.screen.height = height; }
void SetWindowOpacity(float opacity) { (void)opacity; }
void SetWindowFocused(void) { }
void *GetWindowHandle(void) { return NULL; }
int GetMonitorCount(void) { return 1; }
int GetCurrentMonitor(void) { return 0; }
Vector2 GetMonitorPosition(int monitor) { (void)monitor; return (Vector2){ 0, 0 }; }
int GetMonitorWidth(int monitor) { (void)monitor; return 1920; }   // a nominal monitor: callers size windows from it
int GetMonitorHeight(int monitor) { (void)monitor; return 1080; }
int GetMonitorPhysicalWidth(int monitor) { (void)monitor; return 527; }
int GetMonitorPhysicalHeight(int monitor) { (void)monitor; return 296; }
int GetMonitorRefreshRate(int monitor) { (void)monitor; return 60; }
const char *GetMonitorName(int monitor) { (void)monitor; return "headless"; }
Vector2 GetWindowPosition(void) { return (Vector2){ (float)CORE.Window.position.x, (float)CORE.Window.position.y }; }
Vector2 GetWindowScaleDPI(void) { return (Vector2){ 1.0f, 1.0f }; }
void SetClipboardText(const char *text) { (void)text; }
const char *GetClipboardText(void) { return ""; }
Image GetClipboardImage(void) { Image image = { 0 }; return image; }
void ShowCursor(void) { CORE.Input.Mouse.cursorHidden = false; }
void HideCursor(void) { CORE.Input.Mouse.cursorHidden = true; }
void EnableCursor(void) { CORE.Input.Mouse.cursorHidden = false; }
void DisableCursor(void) { CORE.Input.Mouse.cursorHidden = true; }

void SwapScreenBuffer(void)
{
    if (platform.glMode == 1 && platform.surface) egl.SwapBuffers(platform.display, platform.surface);
}

double GetTime(void)
{
    struct timespec ts = { 0 };
    clock_gettime(CLOCK_MONOTONIC, &ts);
    const unsigned long long int nanoSeconds = (unsigned long long int)ts.tv_sec*1000000000LLU + (unsigned long long int)ts.tv_nsec;
    return (double)(nanoSeconds - CORE.Time.base)*1e-9; // since InitTimer()
}

void OpenURL(const char *url) { TRACELOG(LOG_INFO, "HEADLESS: OpenURL(%s) ignored", url ? url : ""); }
int SetGamepadMappings(const char *mappings) { (void)mappings; return 0; }
void SetMousePosition(int x, int y)
{
    CORE.Input.Mouse.currentPosition = (Vector2){ (float)x, (float)y };
    CORE.Input.Mouse.previousPosition = CORE.Input.Mouse.currentPosition;
}
void SetMouseCursor(int cursor) { CORE.Input.Mouse.cursor = cursor; }
const char *GetKeyName(int key) { (void)key; return ""; }

// No input devices: keep raylib's per-frame bookkeeping so Is*Pressed() stay false.
void PollInputEvents(void)
{
#if defined(SUPPORT_GESTURES_SYSTEM)
    UpdateGestures();
#endif
    CORE.Input.Keyboard.keyPressedQueueCount = 0;
    CORE.Input.Keyboard.charPressedQueueCount = 0;
    CORE.Input.Gamepad.lastButtonPressed = 0;
    for (int i = 0; i < MAX_KEYBOARD_KEYS; i++)
    {
        CORE.Input.Keyboard.previousKeyState[i] = CORE.Input.Keyboard.currentKeyState[i];
        CORE.Input.Keyboard.keyRepeatInFrame[i] = 0;
    }
    for (int i = 0; i < MAX_MOUSE_BUTTONS; i++) CORE.Input.Mouse.previousButtonState[i] = CORE.Input.Mouse.currentButtonState[i];
    CORE.Input.Mouse.previousWheelMove = CORE.Input.Mouse.currentWheelMove;
    CORE.Input.Mouse.currentWheelMove = (Vector2){ 0.0f, 0.0f };
    CORE.Input.Mouse.previousPosition = CORE.Input.Mouse.currentPosition;
    for (int i = 0; i < MAX_TOUCH_POINTS; i++) CORE.Input.Touch.previousTouchState[i] = CORE.Input.Touch.currentTouchState[i];
}

int InitPlatform(void)
{
    if (CORE.Window.screen.width <= 0) CORE.Window.screen.width = 320;
    if (CORE.Window.screen.height <= 0) CORE.Window.screen.height = 180;
    CORE.Window.display.width = CORE.Window.screen.width;
    CORE.Window.display.height = CORE.Window.screen.height;
    CORE.Window.render.width = CORE.Window.screen.width;
    CORE.Window.render.height = CORE.Window.screen.height;
    CORE.Window.currentFbo.width = CORE.Window.render.width;
    CORE.Window.currentFbo.height = CORE.Window.render.height;
    CORE.Window.flags |= FLAG_WINDOW_HIDDEN;
    CORE.Window.flags &= ~FLAG_MSAA_4X_HINT;

    if (InitEGLContext())
    {
        platform.glMode = 1;
        rlLoadExtensions(HeadlessGetProcAddress);
        TRACELOG(LOG_WARNING, "HEADLESS: %s", platform.description); // (warning level: the game logs from WARNING up)
    }
    else
    {
        platform.glMode = 2;
        snprintf(platform.description, sizeof(platform.description), "null GL (no GPU context: GPU paths off, the sim runs on the CPU)");
        rlLoadExtensions(NullGL_GetProcAddress);
        TRACELOG(LOG_WARNING, "HEADLESS: %s", platform.description);
    }
    CORE.Window.ready = true;

    InitTimer();
    CORE.Storage.basePath = GetWorkingDirectory();
    TRACELOG(LOG_INFO, "PLATFORM: HEADLESS: Initialized successfully");
    return 0;
}

void ClosePlatform(void)
{
    if (platform.glMode == 1)
    {
        egl.MakeCurrent(platform.display, NULL, NULL, NULL);
        if (platform.surface) egl.DestroySurface(platform.display, platform.surface);
        if (platform.context) egl.DestroyContext(platform.display, platform.context);
        egl.Terminate(platform.display);
    }
    platform.glMode = 0;
}
