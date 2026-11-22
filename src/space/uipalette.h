#ifndef UIPALETTE_H
#define UIPALETTE_H

#include <raylib.h>

// The colours that carry a meaning (#1432, part 1 of #1017): good / bad / warn / info, chosen per colour vision setting
// (Settings > Gameplay > Colour vision) so the pairs that matter stay apart for players who can't tell red from green
// (protan, deutan) or blue from yellow (tritan). Normal keeps the game's colours; the others are from the Okabe-Ito set.
// Wherever colour alone told the state, the state now has a second cue as well (the can't-place build ghost is crossed).
enum class ColorVision { Normal = 0, RedGreen = 1, BlueYellow = 2, Count };

struct UiPalette {
    Color good;   // can place, ready, safe, full
    Color bad;    // can't place, danger, hostile, empty
    Color warn;   // in between: ripe, wild, low
    Color info;   // neutral highlights (a selection, the start marker)
};

void UiPalette_Set(ColorVision v);
ColorVision UiPalette_Mode();
const UiPalette& Ui();                                 // the current mode's colours
const UiPalette& UiPaletteFor(ColorVision v);          // any mode's (tests, the settings preview)
Color UiPalette_Lerp(float t);                         // bad (0) .. good (1), through warn: a meter's colour
const char* ColorVisionName(ColorVision v);

#endif
