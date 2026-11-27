// (#2032) Per-biome tints. The climate code (biomevariants.h BiomeClimateCode): 2 + the temperature step (0-62 over
// -30..40 C) + 64 x the biome variant id + moisture as the fraction; below 1.5 = no climate. BiomeTint_Bind sets the
// arrays: one grass and one leaf multiplier per variant id (1 with RAYLIBGAME_BIOME_TINT=0).
uniform vec3 biomeGrassTint[64];
uniform vec3 biomeLeafTint[64];
uniform float biomeClimateWeight; // how much of the climate gradient stays inside a biome (1 = all, the tints off)

float biomeStepOf(float code) { float n = floor(code) - 2.0; return n - 64.0 * floor(n / 64.0); }
float biomeIdOf(float code) { return code < 1.5 ? 0.0 : floor((floor(code) - 2.0) / 64.0); }
float biomeTempOf(float code) { return biomeStepOf(code) / 62.0 * 70.0 - 30.0; }
vec3 biomeGrass(float id) { return biomeGrassTint[int(clamp(id, 0.0, 63.0))]; }
vec3 biomeLeaf(float id) { return biomeLeafTint[int(clamp(id, 0.0, 63.0))]; }
// The climate colour pulled toward the temperate one: the biome's tint carries the main signal.
vec3 biomeClimate(vec3 climateColour, vec3 temperate) { return mix(temperate, climateColour, biomeClimateWeight); }
