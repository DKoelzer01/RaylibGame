#include "gpumem.h"
#include "external/glad.h"
#include <algorithm>
#include <cstdio>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstring>
#if defined(_WIN32)
extern "C" __declspec(dllimport) unsigned short __stdcall RtlCaptureStackBackTrace(unsigned long, unsigned long, void**, unsigned long*); // kernel32
#elif defined(__GLIBC__)
#include <execinfo.h>
#endif

namespace {
    // (#1603) The return addresses above the GL wrapper, so an untagged allocation names its caller (raylib's rlLoad* and
    // the game code that called it). Only for "other" allocations, and only with --vram (the wrappers are installed then).
    void CaptureStack(void* (&out)[6]) {
#if defined(_WIN32)
        RtlCaptureStackBackTrace(2, 6, out, nullptr); // skip CaptureStack and Set
#elif defined(__GLIBC__)
        void* f[8] = {};
        const int n = backtrace(f, 8);
        for (int i = 2; i < n && i < 8; ++i) out[i - 2] = f[i];
#else
        (void)out;
#endif
    }
    bool active = false;
    const char* currentTag = "other";
    struct Obj {
        const char* tag = "other";
        void* stack[6] = {}; // (#1603) who allocated an untagged one (return addresses: the GL caller, its callers)
        char desc[40] = "";
        std::map<int, size_t> parts; // per mip level / cube face (re-specifying one replaces it)
        size_t total() const { size_t s = 0; for (auto& p : parts) s += p.second; return s; }
    };
    std::unordered_map<GLuint, Obj> buffers, textures, renderbuffers;

    PFNGLBUFFERDATAPROC oBufferData; PFNGLBUFFERSTORAGEPROC oBufferStorage; PFNGLNAMEDBUFFERSTORAGEPROC oNamedBufferStorage;
    PFNGLDELETEBUFFERSPROC oDeleteBuffers;
    PFNGLTEXIMAGE2DPROC oTexImage2D; PFNGLTEXIMAGE3DPROC oTexImage3D; PFNGLTEXSTORAGE2DPROC oTexStorage2D; PFNGLTEXSTORAGE3DPROC oTexStorage3D;
    PFNGLTEXIMAGE2DMULTISAMPLEPROC oTexImage2DMS; PFNGLTEXSTORAGE2DMULTISAMPLEPROC oTexStorage2DMS; PFNGLCOMPRESSEDTEXIMAGE2DPROC oCompressed2D;
    PFNGLGENERATEMIPMAPPROC oGenerateMipmap; PFNGLDELETETEXTURESPROC oDeleteTextures;
    PFNGLRENDERBUFFERSTORAGEPROC oRbStorage; PFNGLRENDERBUFFERSTORAGEMULTISAMPLEPROC oRbStorageMS; PFNGLDELETERENDERBUFFERSPROC oDeleteRenderbuffers;

    GLenum BufferBinding(GLenum target) {
        switch (target) {
            case GL_ARRAY_BUFFER: return GL_ARRAY_BUFFER_BINDING;
            case GL_ELEMENT_ARRAY_BUFFER: return GL_ELEMENT_ARRAY_BUFFER_BINDING;
            case GL_COPY_READ_BUFFER: return GL_COPY_READ_BUFFER_BINDING;
            case GL_COPY_WRITE_BUFFER: return GL_COPY_WRITE_BUFFER_BINDING;
            case GL_PIXEL_PACK_BUFFER: return GL_PIXEL_PACK_BUFFER_BINDING;
            case GL_PIXEL_UNPACK_BUFFER: return GL_PIXEL_UNPACK_BUFFER_BINDING;
            case GL_UNIFORM_BUFFER: return GL_UNIFORM_BUFFER_BINDING;
            case GL_SHADER_STORAGE_BUFFER: return GL_SHADER_STORAGE_BUFFER_BINDING;
            case GL_DRAW_INDIRECT_BUFFER: return GL_DRAW_INDIRECT_BUFFER_BINDING;
            case GL_DISPATCH_INDIRECT_BUFFER: return GL_DISPATCH_INDIRECT_BUFFER_BINDING;
            case GL_TRANSFORM_FEEDBACK_BUFFER: return GL_TRANSFORM_FEEDBACK_BUFFER_BINDING;
            case GL_ATOMIC_COUNTER_BUFFER: return GL_ATOMIC_COUNTER_BUFFER_BINDING;
            default: return 0;
        }
    }
    GLenum TextureBinding(GLenum target) {
        switch (target) {
            case GL_TEXTURE_2D: return GL_TEXTURE_BINDING_2D;
            case GL_TEXTURE_2D_ARRAY: return GL_TEXTURE_BINDING_2D_ARRAY;
            case GL_TEXTURE_3D: return GL_TEXTURE_BINDING_3D;
            case GL_TEXTURE_2D_MULTISAMPLE: return GL_TEXTURE_BINDING_2D_MULTISAMPLE;
            case GL_TEXTURE_CUBE_MAP: case GL_TEXTURE_CUBE_MAP_POSITIVE_X: case GL_TEXTURE_CUBE_MAP_NEGATIVE_X: case GL_TEXTURE_CUBE_MAP_POSITIVE_Y:
            case GL_TEXTURE_CUBE_MAP_NEGATIVE_Y: case GL_TEXTURE_CUBE_MAP_POSITIVE_Z: case GL_TEXTURE_CUBE_MAP_NEGATIVE_Z: return GL_TEXTURE_BINDING_CUBE_MAP;
            case GL_TEXTURE_1D: return GL_TEXTURE_BINDING_1D;
            default: return 0;
        }
    }
    GLuint Bound(GLenum binding) { if (!binding) return 0; GLint id = 0; glGetIntegerv(binding, &id); return (GLuint)id; }
    int Face(GLenum target) { return target >= GL_TEXTURE_CUBE_MAP_POSITIVE_X && target <= GL_TEXTURE_CUBE_MAP_NEGATIVE_Z ? (int)(target - GL_TEXTURE_CUBE_MAP_POSITIVE_X) : 0; }
    size_t Bpp(GLenum f) {
        switch (f) {
            case GL_R8: case GL_RED: case GL_R8I: case GL_R8UI: return 1;
            case GL_RG8: case GL_R16F: case GL_R16: case GL_R16I: case GL_R16UI: case GL_DEPTH_COMPONENT16: return 2;
            case GL_RGB8: case GL_RGB: case GL_SRGB8: return 3;
            case GL_RGB16F: return 6;
            case GL_DEPTH32F_STENCIL8: case GL_RGBA16F: case GL_RG32F: case GL_RGBA16: case GL_RG32I: case GL_RG32UI: case GL_RGBA16I: case GL_RGBA16UI: return 8;
            case GL_RGB32F: return 12;
            case GL_RGBA32F: case GL_RGBA32I: case GL_RGBA32UI: return 16;
            default: return 4; // RGBA8, R32F/I, RG16F, depth 24/32, depth-stencil, unknown
        }
    }
    void Set(std::unordered_map<GLuint, Obj>& m, GLuint id, int part, size_t bytes, const char* desc = nullptr) {
        if (!id) id = 0xFFFFFFFFu; // unknown binding: still counted
        Obj& o = m[id];
        if (o.parts.empty()) { o.tag = currentTag; if (std::strcmp(currentTag, "other") == 0) CaptureStack(o.stack); if (desc) snprintf(o.desc, sizeof(o.desc), "%s", desc); }
        o.parts[part] = bytes;
    }
    void Erase(std::unordered_map<GLuint, Obj>& m, GLsizei n, const GLuint* ids) { for (GLsizei i = 0; i < n; ++i) m.erase(ids[i]); }

    // (#1603) what was made, for the "other" list
    char gDesc[40];
    template <typename... A> const char* Desc(const char* f, A... a) { snprintf(gDesc, sizeof(gDesc), f, a...); return gDesc; }

    void GLAD_API_PTR wBufferData(GLenum t, GLsizeiptr size, const void* d, GLenum u) { oBufferData(t, size, d, u); Set(buffers, Bound(BufferBinding(t)), 0, (size_t)size, Desc("buffer 0x%x", t)); }
    void GLAD_API_PTR wBufferStorage(GLenum t, GLsizeiptr size, const void* d, GLbitfield f) { oBufferStorage(t, size, d, f); Set(buffers, Bound(BufferBinding(t)), 0, (size_t)size, Desc("buffer 0x%x", t)); }
    void GLAD_API_PTR wNamedBufferStorage(GLuint b, GLsizeiptr size, const void* d, GLbitfield f) { oNamedBufferStorage(b, size, d, f); Set(buffers, b, 0, (size_t)size, "named buffer"); }
    void GLAD_API_PTR wDeleteBuffers(GLsizei n, const GLuint* ids) { Erase(buffers, n, ids); oDeleteBuffers(n, ids); }
    void GLAD_API_PTR wTexImage2D(GLenum t, GLint level, GLint ifmt, GLsizei w, GLsizei h, GLint b, GLenum fmt, GLenum type, const void* px) {
        oTexImage2D(t, level, ifmt, w, h, b, fmt, type, px);
        if (t == GL_PROXY_TEXTURE_2D) return;
        Set(textures, Bound(TextureBinding(t)), level * 8 + Face(t), (size_t)w * h * Bpp((GLenum)ifmt), Desc("tex2d %dx%d fmt 0x%x", w, h, ifmt));
    }
    void GLAD_API_PTR wTexImage3D(GLenum t, GLint level, GLint ifmt, GLsizei w, GLsizei h, GLsizei dd, GLint b, GLenum fmt, GLenum type, const void* px) {
        oTexImage3D(t, level, ifmt, w, h, dd, b, fmt, type, px);
        Set(textures, Bound(TextureBinding(t)), level * 8, (size_t)w * h * dd * Bpp((GLenum)ifmt), Desc("tex3d %dx%dx%d fmt 0x%x", w, h, dd, ifmt));
    }
    void GLAD_API_PTR wTexStorage2D(GLenum t, GLsizei levels, GLenum ifmt, GLsizei w, GLsizei h) {
        oTexStorage2D(t, levels, ifmt, w, h);
        size_t s = 0;
        for (int l = 0; l < levels; ++l) s += (size_t)std::max(1, w >> l) * std::max(1, h >> l) * Bpp(ifmt);
        Set(textures, Bound(TextureBinding(t)), 0, s * (t == GL_TEXTURE_CUBE_MAP ? 6 : 1), Desc("storage2d %dx%d x%d fmt 0x%x", w, h, levels, ifmt));
    }
    void GLAD_API_PTR wTexStorage3D(GLenum t, GLsizei levels, GLenum ifmt, GLsizei w, GLsizei h, GLsizei dd) {
        oTexStorage3D(t, levels, ifmt, w, h, dd);
        size_t s = 0;
        for (int l = 0; l < levels; ++l) s += (size_t)std::max(1, w >> l) * std::max(1, h >> l) * (t == GL_TEXTURE_3D ? std::max(1, dd >> l) : dd) * Bpp(ifmt);
        Set(textures, Bound(TextureBinding(t)), 0, s, Desc("storage3d %dx%dx%d fmt 0x%x", w, h, dd, ifmt));
    }
    void GLAD_API_PTR wTexImage2DMS(GLenum t, GLsizei samples, GLenum ifmt, GLsizei w, GLsizei h, GLboolean fixed) {
        oTexImage2DMS(t, samples, ifmt, w, h, fixed);
        Set(textures, Bound(TextureBinding(t)), 0, (size_t)w * h * Bpp(ifmt) * std::max(1, (int)samples), Desc("tex2d ms %dx%d fmt 0x%x", w, h, ifmt));
    }
    void GLAD_API_PTR wTexStorage2DMS(GLenum t, GLsizei samples, GLenum ifmt, GLsizei w, GLsizei h, GLboolean fixed) {
        oTexStorage2DMS(t, samples, ifmt, w, h, fixed);
        Set(textures, Bound(TextureBinding(t)), 0, (size_t)w * h * Bpp(ifmt) * std::max(1, (int)samples), Desc("storage2d ms %dx%d fmt 0x%x", w, h, ifmt));
    }
    void GLAD_API_PTR wCompressed2D(GLenum t, GLint level, GLenum ifmt, GLsizei w, GLsizei h, GLint b, GLsizei size, const void* d) {
        oCompressed2D(t, level, ifmt, w, h, b, size, d);
        Set(textures, Bound(TextureBinding(t)), level * 8 + Face(t), (size_t)size, Desc("compressed %dx%d", w, h));
    }
    void GLAD_API_PTR wGenerateMipmap(GLenum t) {
        oGenerateMipmap(t);
        auto it = textures.find(Bound(TextureBinding(t)));
        if (it != textures.end()) { const size_t base = it->second.parts.count(0) ? it->second.parts[0] : 0; it->second.parts[1000] = base / 3; } // the chain: ~1/3 more
    }
    void GLAD_API_PTR wDeleteTextures(GLsizei n, const GLuint* ids) { Erase(textures, n, ids); oDeleteTextures(n, ids); }
    void GLAD_API_PTR wRbStorage(GLenum t, GLenum ifmt, GLsizei w, GLsizei h) { oRbStorage(t, ifmt, w, h); Set(renderbuffers, Bound(GL_RENDERBUFFER_BINDING), 0, (size_t)w * h * Bpp(ifmt), Desc("rb %dx%d fmt 0x%x", w, h, ifmt)); }
    void GLAD_API_PTR wRbStorageMS(GLenum t, GLsizei samples, GLenum ifmt, GLsizei w, GLsizei h) {
        oRbStorageMS(t, samples, ifmt, w, h);
        Set(renderbuffers, Bound(GL_RENDERBUFFER_BINDING), 0, (size_t)w * h * Bpp(ifmt) * std::max(1, (int)samples), Desc("rb ms %dx%d fmt 0x%x", w, h, ifmt));
    }
    void GLAD_API_PTR wDeleteRenderbuffers(GLsizei n, const GLuint* ids) { Erase(renderbuffers, n, ids); oDeleteRenderbuffers(n, ids); }

    template <typename F> void Hook(F& slot, F& orig, F wrapper) { if (slot) { orig = slot; slot = wrapper; } }
}

void GpuMem_Install() {
    if (active) return;
    active = true;
    Hook(glad_glBufferData, oBufferData, wBufferData);
    Hook(glad_glBufferStorage, oBufferStorage, wBufferStorage);
    Hook(glad_glNamedBufferStorage, oNamedBufferStorage, wNamedBufferStorage);
    Hook(glad_glDeleteBuffers, oDeleteBuffers, wDeleteBuffers);
    Hook(glad_glTexImage2D, oTexImage2D, wTexImage2D);
    Hook(glad_glTexImage3D, oTexImage3D, wTexImage3D);
    Hook(glad_glTexStorage2D, oTexStorage2D, wTexStorage2D);
    Hook(glad_glTexStorage3D, oTexStorage3D, wTexStorage3D);
    Hook(glad_glTexImage2DMultisample, oTexImage2DMS, wTexImage2DMS);
    Hook(glad_glTexStorage2DMultisample, oTexStorage2DMS, wTexStorage2DMS);
    Hook(glad_glCompressedTexImage2D, oCompressed2D, wCompressed2D);
    Hook(glad_glGenerateMipmap, oGenerateMipmap, wGenerateMipmap);
    Hook(glad_glDeleteTextures, oDeleteTextures, wDeleteTextures);
    Hook(glad_glRenderbufferStorage, oRbStorage, wRbStorage);
    Hook(glad_glRenderbufferStorageMultisample, oRbStorageMS, wRbStorageMS);
    Hook(glad_glDeleteRenderbuffers, oDeleteRenderbuffers, wDeleteRenderbuffers);
}

bool GpuMem_Active() { return active; }

GpuMemTag::GpuMemTag(const char* tag) : prev(currentTag) { currentTag = tag; }
GpuMemTag::~GpuMemTag() { currentTag = prev; }

size_t GpuMem_Total() {
    size_t s = 0;
    for (auto* m : { &buffers, &textures, &renderbuffers }) for (auto& kv : *m) s += kv.second.total();
    return s;
}

void GpuMem_Report(const char* label) {
    if (!active) return;
    std::map<std::string, size_t> all;
    size_t nBuf = buffers.size(), nTex = textures.size() + renderbuffers.size();
    for (auto& kv : buffers) all[std::string(kv.second.tag)] += kv.second.total();
    for (auto& kv : textures) all[std::string(kv.second.tag)] += kv.second.total();
    for (auto& kv : renderbuffers) all[std::string(kv.second.tag)] += kv.second.total();
    std::vector<std::pair<size_t, std::string>> rows;
    for (auto& kv : all) rows.push_back({ kv.second, kv.first });
    std::sort(rows.rbegin(), rows.rend());
    const double MB = 1024.0 * 1024.0;
    printf("bench: vram (%s): counted %.0f MB in %zu buffers and %zu textures/renderbuffers:", label, GpuMem_Total() / MB, nBuf, nTex);
    for (auto& r : rows) printf(" %s %.1f", r.second.c_str(), r.first / MB);
    printf(" (MB)\n");
    {   // (#1603) The biggest untagged objects: what they are and who made them (return addresses in the executable; the
        // anchor is GpuMem_Report's own address, so `nm game.exe` maps one: its address = it - anchor + nm's GpuMem_Report).
        std::vector<std::pair<size_t, const Obj*>> other;
        for (auto* m : { &buffers, &textures, &renderbuffers })
            for (auto& kv : *m) if (std::string(kv.second.tag) == "other") other.push_back({ kv.second.total(), &kv.second });
        std::sort(other.begin(), other.end(), [](auto& a, auto& b) { return a.first > b.first; });
        printf("bench: vram (%s): other: %zu objects, anchor %p; the largest:", label, other.size(), (const void*)&GpuMem_Report);
        for (size_t i = 0; i < other.size() && i < 16; ++i) {
            printf(" [%.1f MB %s @", other[i].first / MB, other[i].second->desc);
            for (void* p : other[i].second->stack) if (p) printf(" %p", p);
            printf("]");
        }
        printf("\n");
    }
    // What the driver says (all processes on the GPU): NVX_gpu_memory_info, else ATI_meminfo.
    GLint total = 0, avail = 0;
    while (glGetError() != GL_NO_ERROR) {}
    glGetIntegerv(0x9048, &total);  // GL_GPU_MEMORY_INFO_TOTAL_AVAILABLE_MEMORY_NVX (KB)
    glGetIntegerv(0x9049, &avail);  // GL_GPU_MEMORY_INFO_CURRENT_AVAILABLE_VIDMEM_NVX (KB)
    if (glGetError() == GL_NO_ERROR && total > 0) {
        printf("bench: vram (%s): driver: %.0f of %.0f MB in use on the GPU (every process)\n", label, (total - avail) / 1024.0, total / 1024.0);
        return;
    }
    GLint ati[4] = { 0, 0, 0, 0 };
    while (glGetError() != GL_NO_ERROR) {}
    glGetIntegerv(0x87FC, ati); // GL_TEXTURE_FREE_MEMORY_ATI (KB)
    if (glGetError() == GL_NO_ERROR && ati[0] > 0) printf("bench: vram (%s): driver: %.0f MB free for textures\n", label, ati[0] / 1024.0);
}
