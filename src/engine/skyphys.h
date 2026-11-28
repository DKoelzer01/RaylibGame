#pragma once

// #142: the physical sky (Hillaire 2020 look-up tables) on by default (user, 2026-11-18); RAYLIBGAME_SKY_PHYS=0 turns it off.
// Transmittance (256x64) and multiple scattering (32x32) are built once per planet radius / air top; the sky-view LUT
// (192x108) follows the camera height and the sun's elevation. All RGBA16F, linear filtered, rendered by fullscreen
// fragment passes (resources/skyphys.glsl holds the shared model, and atmosphere.fs samples the sky view).

bool SkyPhys_On();               // on unless RAYLIBGAME_SKY_PHYS=0, and the LUT shaders compiled
float SkyPhys_Tint();            // how far atmosphere.fs pulls the result toward today's sky palette (RAYLIBGAME_SKY_PHYS_TINT, 0.35)
float SkyPhys_Exposure();        // sun illuminance x exposure into display (RAYLIBGAME_SKY_PHYS_EXPOSURE, 10)

// Rebuilds what changed. Call between passes (it binds its own framebuffer and restores the framebuffer and viewport).
// camDist: the camera's distance from the planet centre (m); sunCosZenith: dot(up at the camera, toward the sun).
void SkyPhys_Update(float planetRadius, float atmosphereTop, float camDist, float sunCosZenith);
// The LUTs' Earth-sized shell (km; the sky-view lookup takes these as skyR / skyRt) and the camera's height in it.
float SkyPhys_LutR();
float SkyPhys_LutRt();
float SkyPhys_LutHeight(float planetRadius, float atmosphereTop, float camDist);
unsigned SkyPhys_SkyViewTexture();       // 0 until built
unsigned SkyPhys_HazeTexture();          // (#819) the haze colour LUT for lit surfaces (64 x 128); 0 until built
unsigned SkyPhys_TransmittanceTexture(); // 0 until built
long SkyPhys_Rebuilds();                 // (#2143) sky-view + haze rebuilds since the start (tests)
