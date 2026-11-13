#ifndef TERRAINLAYERS_H
#define TERRAINLAYERS_H

// Terrain texture layers: one OpenGL 2D texture array (one sampler, one texture unit, any number of
// layers; mipmapped and repeating, so no atlas bleeding). raylib has no texture-array API, so this
// uses the GL loader raylib links (its glad function pointers are public symbols).
//
// Layers are placeholder programmer art generated at startup (seamless noise / cell patterns);
// a file assets/textures/terrain/<name>.png replaces a layer when present (resized to the layer size).

#include <raylib.h>

enum TerrainLayer { LayerDirt = 0, LayerSand, LayerGravel, LayerGrass, LayerSnow, LayerIce, LayerStone, LayerCount };
const char* TerrainLayerName(int layer);

// Terrain layers come with a second array, the detail maps (Plan 5 D, stylized PBR): RG = tangent-space
// normal XY, B = roughness, A = ambient occlusion; generated from each layer's albedo brightness, or
// assets/textures/terrain/<name>_detail.png. Bound with the colour array to unit kTerrainDetailUnit.
constexpr int kTerrainDetailUnit = 12;
struct TextureArray {
    unsigned int id = 0, detailId = 0;
    int size = 0, layers = 0;
};

TextureArray LoadTerrainLayers(int size = 512);   // after InitWindow
void UnloadTextureArray(TextureArray& t);
void BindTextureArray(int unit, const TextureArray& t);
void UnbindTextureArray(int unit);
Image GenerateTerrainLayerImage(int layer, int size); // the placeholder art for one layer
Image GenerateTerrainLayerDetail(int layer, const Image& albedo); // its detail map (RGBA8 albedo in)

#endif
