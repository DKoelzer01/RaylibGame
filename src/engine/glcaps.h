#ifndef GLCAPS_H
#define GLCAPS_H

// OpenGL capabilities above the 4.3 floor (#286): detected once on the main thread after the context exists, logged
// in one line, and used by every optional path (each keeps its 4.3 fallback). RAYLIBGAME_GL_LIMIT=43 turns all of
// them off, to exercise the fallbacks on a machine that has them.
//
//   clipControl          4.5 / ARB_clip_control: reverse-Z (scenetarget)
//   bufferStorage        4.4 / ARB_buffer_storage: persistent mapped buffers (GPU terrain readback ring)
//   directStateAccess    4.5 / ARB_direct_state_access
//   indirectParameters   4.6 / ARB_indirect_parameters: glMultiDrawElementsIndirectCount
//   shaderDrawParameters 4.6 / ARB_shader_draw_parameters: gl_DrawID, gl_BaseInstance in shaders
//   anisotropic          4.6 / ARB_texture_filter_anisotropic or EXT_texture_filter_anisotropic

#include <string>

struct GLCaps {
    int major = 0, minor = 0;
    std::string vendor, renderer;
    bool limited = false;              // RAYLIBGAME_GL_LIMIT=43: everything above 4.3 reported off
    bool clipControl = false;
    bool bufferStorage = false;
    bool directStateAccess = false;
    bool indirectParameters = false;
    bool shaderDrawParameters = false;
    bool anisotropic = false;
    float maxAnisotropy = 1.0f;
};

const GLCaps& GL_Caps();               // detects on the first call (main thread, context current)

#endif
