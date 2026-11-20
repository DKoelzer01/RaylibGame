#ifndef DEBUGUI_H
#define DEBUGUI_H

class Scene;

// Dear ImGui tweak panel (rlImGui backend). F1 toggles it; while open the mouse is released and
// camera input is paused so you can use the widgets.
void DebugUI_Init();
void DebugUI_Shutdown();
bool DebugUI_IsOpen();
void DebugUI_Toggle();
// Call inside rlImGuiBegin()/rlImGuiEnd() (main loop), after the scene and HUD.
void DebugUI_Draw(Scene& world);

#endif
