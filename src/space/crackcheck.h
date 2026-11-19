#pragma once
// --crack-check (#725): holes in the terrain show the magenta clear colour. CrackHoles counts the magenta samples enclosed
// by drawn pixels (holes), not those connected to the top edge (open sky; skyOut, when given). step: every step-th pixel.
#include <raylib.h>
long CrackHoles(const Color* px, int w, int h, int step, long* skyOut);
