#ifndef ROTATIONMODELS_H
#define ROTATIONMODELS_H

// Rotation models (#159, #1269): shafts, cogwheels and hand cranks drawn turning at the speed the factory reports
// (Factory::RotationSpeeds), and the turning part on top of the rotation machines (the millstone's runner stone, the saw's
// blade, the press's flywheel). A slowed network visibly turns slower. Plots within 48 m of the camera.

class Scene;
void RotationModels_Draw3D(Scene& world);
int RotationModels_DrawnCount(); // tests: parts drawn last frame

#endif
