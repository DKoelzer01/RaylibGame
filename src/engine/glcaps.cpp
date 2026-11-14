#include "glcaps.h"
#include "utils/logger.h"

#include "external/glad.h"
#include <cstdlib>
#include <cstring>
#include <unordered_set>

const GLCaps& GL_Caps() {
    static GLCaps c;
    static bool done = false;
    if (done) return c;
    done = true;
    glGetIntegerv(GL_MAJOR_VERSION, &c.major);
    glGetIntegerv(GL_MINOR_VERSION, &c.minor);
    if (const char* v = (const char*)glGetString(GL_VENDOR)) c.vendor = v;
    if (const char* r = (const char*)glGetString(GL_RENDERER)) c.renderer = r;
    std::unordered_set<std::string> ext;
    GLint n = 0;
    glGetIntegerv(GL_NUM_EXTENSIONS, &n);
    for (GLint i = 0; i < n; ++i)
        if (const char* e = (const char*)glGetStringi(GL_EXTENSIONS, (GLuint)i)) ext.insert(e);
    auto at = [&](int ma, int mi) { return c.major > ma || (c.major == ma && c.minor >= mi); };
    auto has = [&](const char* e) { return ext.count(e) > 0; };
    const char* lim = getenv("RAYLIBGAME_GL_LIMIT");
    c.limited = lim && strcmp(lim, "43") == 0;
    if (!c.limited) {
        c.clipControl = at(4, 5) || has("GL_ARB_clip_control");
        c.bufferStorage = (at(4, 4) || has("GL_ARB_buffer_storage")) && glad_glBufferStorage && glad_glMapBufferRange;
        c.directStateAccess = at(4, 5) || has("GL_ARB_direct_state_access");
        c.indirectParameters = at(4, 6) || has("GL_ARB_indirect_parameters");
        c.shaderDrawParameters = at(4, 6) || has("GL_ARB_shader_draw_parameters");
        c.anisotropic = at(4, 6) || has("GL_ARB_texture_filter_anisotropic") || has("GL_EXT_texture_filter_anisotropic");
        if (c.anisotropic) glGetFloatv(0x84FF, &c.maxAnisotropy); // GL_MAX_TEXTURE_MAX_ANISOTROPY
    }
    logger.logf("[GL] %d.%d, %s / %s%s: clip control %s, buffer storage %s, DSA %s, indirect count %s, draw parameters %s, anisotropy %s (x%.0f)\n",
                c.major, c.minor, c.vendor.c_str(), c.renderer.c_str(), c.limited ? " (RAYLIBGAME_GL_LIMIT=43)" : "",
                c.clipControl ? "yes" : "no", c.bufferStorage ? "yes" : "no", c.directStateAccess ? "yes" : "no",
                c.indirectParameters ? "yes" : "no", c.shaderDrawParameters ? "yes" : "no", c.anisotropic ? "yes" : "no", c.maxAnisotropy);
    return c;
}
