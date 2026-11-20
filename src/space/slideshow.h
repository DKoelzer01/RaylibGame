#ifndef SLIDESHOW_H
#define SLIDESHOW_H

// A slow slideshow behind the main menu and the loading screens (#386). Slides come from resources/menu/slides/ (PNG or
// JPG, any size: drawn to cover the window) with an optional slides.txt ("file | caption" lines, # comments) for their
// order-independent captions. Each shows ~9 s with a 1.5 s crossfade, in random order without immediate repeats, carrying
// on across the menu and loading screens; a slow pan and zoom unless Reduce motion is on. Images are decoded on a worker
// thread and only two are on the GPU at a time, so it never slows loading. No slides: nothing is drawn (the caller's
// backdrop shows).

// Draws the slideshow over the whole window, then a darkening gradient so the UI on top stays readable (dim 0..1: how
// dark the gradient's middle is). Returns false when there are no slides (nothing drawn).
bool Slideshow_Draw(float dim = 0.45f);
void Slideshow_Unload();

#endif
