#ifndef INSTANCEBATCH_H
#define INSTANCEBATCH_H

// InstanceBatch (perf: batching): boxes and tapered cylinders collected over a frame and drawn in one
// instanced draw per shape, lit like the parts (sun shadows, headlamp, ambient, haze). Replaces raylib's
// immediate DrawCube / DrawCylinder for many small objects (animals, machine details): those rebuilt every
// vertex on the CPU each frame and were drawn unlit and unshadowed.
//
// The current rlgl transform (rlPushMatrix / rlTranslatef / rlRotatef / rlMultMatrixf) applies as with
// DrawCube. Call InstanceBatch_Flush(world) once at the end of the frame's world drawing (onDraw3D).

#include <raylib.h>

class Scene;

void IB_Cube(Vector3 centre, Vector3 size, Color color, bool glow = false);
// A cylinder from `base` up (local +y) `height`, radius rBottom at the base and rTop at the top.
void IB_Cylinder(Vector3 base, float rTop, float rBottom, float height, Color color, bool glow = false);
// How buried the instances pushed from now on are (0 open .. 1 ~10 m+ underground): scales their sky light, as for
// parts and terrain (#284). Reset to 0 by the caller when done.
void IB_SetBuried(float buried);
void InstanceBatch_Flush(Scene& world);
int InstanceBatch_LastCount();   // instances drawn last flush (tests)
void InstanceBatch_Unload();

#endif
