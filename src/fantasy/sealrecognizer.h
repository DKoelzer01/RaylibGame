#ifndef SEALRECOGNIZER_H
#define SEALRECOGNIZER_H

#include <raylib.h>
#include <random>
#include <string>
#include <vector>

// Recognises hand-drawn spell seals: a closed ring containing a core sign at the centre and
// modifier signs placed toward the ring, plus small dots that add power. Works on stroke point
// lists (what the in-game canvas records), so no image processing is needed.
//
//  ring   : Kasa least-squares circle fit per stroke -> roundness (RMS radial error), angular
//           coverage and closure gap. Quality is the "perfection" of the seal.
//  glyphs : strokes inside the ring grouped by proximity; tiny groups count as power dots.
//  signs  : each glyph classified with the $P point-cloud recogniser (Vatavu, Anthony &
//           Wobbrock 2012), which ignores stroke order/direction and handles multi-stroke signs.

using Stroke = std::vector<Vector2>;

struct RingFit {
    bool found = false;
    int stroke = -1;
    Vector2 center{};
    float radius = 0.0f;
    float roundness = 0.0f;  // 1 = perfect circle
    float coverage = 0.0f;   // fraction of 360 degrees swept
    float closure = 0.0f;    // 1 = ends meet
    float quality = 0.0f;    // combined 0..1
};

struct GlyphMatch {
    std::string sign;        // template name ("flame", "arrow", ...)
    float score = 0.0f;      // $P match score 0..1
    Vector2 center{};
    float size = 0.0f;       // bounding-box diagonal
    float polarR = 0.0f;     // distance from ring centre / ring radius
    float angleDeg = 0.0f;   // around the ring, 0 = right, 90 = up (screen space flipped)
    std::vector<int> strokes;
};

struct SealAnalysis {
    RingFit ring;
    std::vector<GlyphMatch> glyphs;
    int core = -1;                 // index into glyphs
    std::vector<int> modifiers;    // indices into glyphs
    int dots = 0;
    float quality = 0.0f;          // ring quality x core match
    std::string problem;           // empty when the seal is valid
};

SealAnalysis AnalyzeSeal(const std::vector<Stroke>& strokes);

// Recognise one glyph (one or more strokes) against the templates.
GlyphMatch RecognizeGlyph(const std::vector<Stroke>& strokes);

// Names of core signs and modifier signs the templates know.
const std::vector<std::string>& CoreSigns();      // flame, water, wind, light, earth
const std::vector<std::string>& ModifierSigns();  // arrow, orb, line

// Test helper: draw a seal procedurally (ring + signs) the way a player might, with
// `sloppiness` 0 (perfect) .. 1 (very rough) noise. Coordinates are in pixels around `center`.
std::vector<Stroke> SynthesizeSeal(const std::string& core, const std::vector<std::string>& modifiers, int dots,
                                   float sloppiness, Vector2 center, float radius, std::mt19937& rng);

#endif
