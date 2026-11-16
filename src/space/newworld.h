#ifndef NEWWORLD_H
#define NEWWORLD_H

// The New World screen (#504): a seed (random by default), the terrain parameters (worldgenparams.h) and a slowly
// turning preview of the planet made from the base noise only (no erosion bake: the menu has no world, #503). The
// main menu opens it for New Game and calls NewWorld_Draw each frame (immediate mode, inside BeginDrawing) until the
// player creates the world or cancels; NewWorld_Apply then sets the seed and writes the parameters beside the save,
// right before the world is generated.
#include "worldgenparams.h"

struct NewWorldChoice {
    enum Kind { None, Create, Cancel } kind = None;
    unsigned int seed = 0;
    WorldGenParams params;
};

void NewWorld_Open();               // a random seed, default parameters, the preview starts building
NewWorldChoice NewWorld_Draw();     // each frame while the screen is up; None until a button is pressed
void NewWorld_Apply(const NewWorldChoice& c); // worldGen.seed, and <WorldSavePath()>.worldgen
void NewWorld_Shutdown();           // waits for a preview build in flight (program exit)
// Tests (--test-newworld): the last preview build's time (s, < 0: none yet), and set the parameters as the sliders do.
double NewWorld_PreviewBuildSeconds();
void NewWorld_TestSetParams(const WorldGenParams& p);
// --seed-survey (#532): the preview's map (same sampling and colours) at w x h for a terrain seed, no GL; the
// elevations too when `heights` is given.
#include <raylib.h>
#include <vector>
void NewWorld_RenderMap(unsigned int terrainSeed, float radius, const WorldGenParams& p, int w, int h, std::vector<Color>& out, std::vector<float>* heights);

#endif
