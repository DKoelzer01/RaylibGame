// --crack-check (#725), shared by main.cpp and the visual/crack-check scenario: see crackcheck.h.
#include "crackcheck.h"
#include <vector>

// --crack-check (#725): magenta (the clear colour, nothing drawn) pixels that are holes in the terrain, not open sky. With
// the skybox off the sky is magenta too: magenta connected to the screen's top edge (4-way, through the gaps in the HUD's
// text) is sky; magenta enclosed by drawn pixels is a hole. Colour test loose enough for the tonemapped composite. The sky
// spreads through magenta-tinted pixels too (anti-aliased edges, where far terrain blends into it), or a horizon pixel pair
// sealing a bit of sky off reads as a hole (the 4-5 left after the HUD went, #725).
// step: sample every step-th pixel both ways. skyOut: the magenta samples taken as sky.
long CrackHoles(const Color* px, int w, int h, int step, long* skyOut) {
    auto magenta = [](Color c) { return c.r > 200 && c.g < 60 && c.b > 200; };
    const int gw = (w + step - 1) / step, gh = (h + step - 1) / step;
    auto tinted = [](Color c) { return (int)c.r - (int)c.g > 25 && (int)c.b - (int)c.g > 25; }; // magenta, or a blend with it
    std::vector<uint8_t> m((size_t)gw * gh), t((size_t)gw * gh), seen((size_t)gw * gh, 0);
    long total = 0;
    for (int y = 0; y < gh; ++y)
        for (int x = 0; x < gw; ++x) { const Color c = px[(size_t)(y * step) * w + x * step]; m[(size_t)y * gw + x] = magenta(c); t[(size_t)y * gw + x] = magenta(c) || tinted(c); total += m[(size_t)y * gw + x]; }
    std::vector<int> stack;
    for (int x = 0; x < gw; ++x) if (t[x]) { seen[x] = 1; stack.push_back(x); }
    long sky = 0;
    while (!stack.empty()) {
        const int i = stack.back(); stack.pop_back(); sky += m[i];
        const int x = i % gw, y = i / gw;
        const int nb[4] = { x > 0 ? i - 1 : -1, x + 1 < gw ? i + 1 : -1, y > 0 ? i - gw : -1, y + 1 < gh ? i + gw : -1 };
        for (int j : nb) if (j >= 0 && t[j] && !seen[j]) { seen[j] = 1; stack.push_back(j); }
    }
    if (skyOut) *skyOut = sky;
    return total - sky;
}
