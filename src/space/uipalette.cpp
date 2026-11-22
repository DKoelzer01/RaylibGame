#include "uipalette.h"

#include <raymath.h>
#include <algorithm>
#include <cstdlib>

namespace {
// Normal: the colours the game used. RedGreen: sky blue against vermillion (they differ in blue-yellow and in lightness,
// which protan and deutan keep). BlueYellow: bluish green against reddish purple (they differ in red-green, which tritan
// keeps). Okabe & Ito, "Color Universal Design" (2008).
const UiPalette kPalettes[] = {
    { { 0, 228, 48, 255 }, { 230, 41, 55, 255 }, { 255, 203, 0, 255 }, { 102, 191, 255, 255 } },   // Normal (raylib GREEN, RED, GOLD, SKYBLUE)
    { { 86, 180, 233, 255 }, { 213, 94, 0, 255 }, { 240, 228, 66, 255 }, { 230, 230, 230, 255 } },  // RedGreen
    { { 0, 158, 115, 255 }, { 204, 121, 167, 255 }, { 240, 240, 240, 255 }, { 0, 114, 178, 255 } }, // BlueYellow
};
// RAYLIBGAME_COLOR_VISION=0/1/2 overrides the setting (benches and captures don't read settings.txt).
int EnvMode() { static const int m = [] { const char* e = getenv("RAYLIBGAME_COLOR_VISION"); return e && *e ? atoi(e) : -1; }(); return m; }
ColorVision gMode = EnvMode() >= 0 ? (ColorVision)EnvMode() : ColorVision::Normal;
}

void UiPalette_Set(ColorVision v) {
    if (EnvMode() >= 0) v = (ColorVision)EnvMode();
    gMode = (int)v >= 0 && v < ColorVision::Count ? v : ColorVision::Normal;
}
ColorVision UiPalette_Mode() { return gMode; }
const UiPalette& UiPaletteFor(ColorVision v) { return kPalettes[(int)v >= 0 && v < ColorVision::Count ? (int)v : 0]; }
const UiPalette& Ui() { return UiPaletteFor(gMode); }
Color UiPalette_Lerp(float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    const UiPalette& p = Ui();
    return t < 0.5f ? ColorLerp(p.bad, p.warn, t * 2.0f) : ColorLerp(p.warn, p.good, t * 2.0f - 1.0f);
}
const char* ColorVisionName(ColorVision v) {
    switch (v) {
    case ColorVision::RedGreen: return "Red-green";
    case ColorVision::BlueYellow: return "Blue-yellow";
    default: return "Normal";
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
#include <cmath>

namespace {
// Colour vision deficiency, simulated (Machado, Oliveira & Fernandes 2009, severity 1.0, on linear RGB), and the difference
// a viewer sees as CIE76 delta E in Lab (D65).
enum class Cvd { Protan, Deutan, Tritan };
const float kMachado[3][9] = {
    { 0.152286f, 1.052583f, -0.204868f, 0.114503f, 0.786281f, 0.099216f, -0.003882f, -0.048116f, 1.051998f },
    { 0.367322f, 0.860646f, -0.227968f, 0.280085f, 0.672501f, 0.047413f, -0.011820f, 0.042940f, 0.968881f },
    { 1.255528f, -0.076749f, -0.178779f, -0.078411f, 0.930809f, 0.147602f, 0.004733f, 0.691367f, 0.303900f },
};
float ToLinear(unsigned char c) { const float v = c / 255.0f; return v <= 0.04045f ? v / 12.92f : powf((v + 0.055f) / 1.055f, 2.4f); }
Vector3 Seen(Color c, const float* m) {
    const float r = ToLinear(c.r), g = ToLinear(c.g), b = ToLinear(c.b);
    if (!m) return { r, g, b };
    return { std::clamp(m[0] * r + m[1] * g + m[2] * b, 0.0f, 1.0f), std::clamp(m[3] * r + m[4] * g + m[5] * b, 0.0f, 1.0f),
             std::clamp(m[6] * r + m[7] * g + m[8] * b, 0.0f, 1.0f) };
}
Vector3 Lab(Vector3 rgb) {
    const float X = (0.4124f * rgb.x + 0.3576f * rgb.y + 0.1805f * rgb.z) / 0.95047f;
    const float Y = 0.2126f * rgb.x + 0.7152f * rgb.y + 0.0722f * rgb.z;
    const float Z = (0.0193f * rgb.x + 0.1192f * rgb.y + 0.9505f * rgb.z) / 1.08883f;
    auto f = [](float t) { return t > 0.008856f ? cbrtf(t) : 7.787f * t + 16.0f / 116.0f; };
    return { 116.0f * f(Y) - 16.0f, 500.0f * (f(X) - f(Y)), 200.0f * (f(Y) - f(Z)) };
}
float DeltaE(Color a, Color b, const float* m) { return Vector3Distance(Lab(Seen(a, m)), Lab(Seen(b, m))); }
constexpr float kApart = 40.0f; // delta E: clearly different at a glance (2.3 is "just noticeable")
}

TEST(UiPalette, GoodAndBadStayApartForTheVisionTheyAreFor) {
    const UiPalette& rg = UiPaletteFor(ColorVision::RedGreen);
    EXPECT_GT(DeltaE(rg.good, rg.bad, kMachado[(int)Cvd::Protan]), kApart);
    EXPECT_GT(DeltaE(rg.good, rg.bad, kMachado[(int)Cvd::Deutan]), kApart);
    EXPECT_GT(DeltaE(rg.good, rg.bad, nullptr), kApart);
    const UiPalette& by = UiPaletteFor(ColorVision::BlueYellow);
    EXPECT_GT(DeltaE(by.good, by.bad, kMachado[(int)Cvd::Tritan]), kApart);
    EXPECT_GT(DeltaE(by.good, by.bad, nullptr), kApart);
    // And the reason for the setting: the normal green and red come much closer under deutan vision than they look normally.
    const UiPalette& n = UiPaletteFor(ColorVision::Normal);
    EXPECT_LT(DeltaE(n.good, n.bad, kMachado[(int)Cvd::Deutan]), DeltaE(n.good, n.bad, nullptr) * 0.5f);
}

TEST(UiPalette, BadIsApartFromWarnToo) {
    for (ColorVision v : { ColorVision::RedGreen, ColorVision::BlueYellow }) {
        const UiPalette& p = UiPaletteFor(v);
        const float* m = kMachado[v == ColorVision::RedGreen ? (int)Cvd::Deutan : (int)Cvd::Tritan];
        EXPECT_GT(DeltaE(p.bad, p.warn, m), 25.0f) << ColorVisionName(v);
    }
}

TEST(UiPalette, TheModeSwitchesTheColours) {
    UiPalette_Set(ColorVision::RedGreen);
    EXPECT_EQ(Ui().good.b, UiPaletteFor(ColorVision::RedGreen).good.b);
    EXPECT_EQ(UiPalette_Lerp(0.0f).r, Ui().bad.r);
    EXPECT_EQ(UiPalette_Lerp(1.0f).b, Ui().good.b);
    UiPalette_Set(ColorVision::Normal);
    EXPECT_EQ(Ui().bad.r, 230);
}
#endif
