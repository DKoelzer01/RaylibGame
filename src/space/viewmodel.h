#ifndef VIEWMODEL_H
#define VIEWMODEL_H

// The hands (plan 2, milestone 9): a gloved hand holding the selected tool or item, low on the
// right of the view. A fixed reference that moves with the body while the world slides past is the
// main comfort aid for first-person walking.
//
// Drawn in its own pass after the world (depth cleared, its own near plane), so it never clips into
// walls. Bobbing follows the distance walked, not time: a gentle figure-eight of ~1.5 cm up and
// ~1 cm across per step, scaled by speed and faded out when standing, airborne or flying; "Reduce
// motion" and the "Hand bob" setting scale it (Full, Reduced, Off). The camera itself never bobs.
// A swing plays on each mining stroke or chop; switching what is held dips the hand.

#include <raylib.h>

class Scene;

void Viewmodel_Update(Scene& world, float dt);
void Viewmodel_Draw(Scene& world);        // after the world, before the 2D HUD
Vector3 Viewmodel_BobOffset();            // the current bob (view space, m) - tests
bool Viewmodel_Visible();                 // drawn this frame

#endif
