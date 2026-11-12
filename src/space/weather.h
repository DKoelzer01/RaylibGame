#ifndef WEATHER_H
#define WEATHER_H

// Weather (plan 3, M23). On a planet with air, the weather is a pure function of where and when:
// a cloud field over the sphere drifts with a slowly turning wind, and the air temperature (survival:
// latitude, height, day and night) decides what falls. Thick cloud rains, or snows below 1 C; very
// thick cloud in warm air is a storm with lightning; around dawn, calm humid air lies as fog. Nothing
// is saved: the same world time gives the same weather.
//
// It greys the sky, dims the sun and thickens the haze (Scene::weatherCloud / weatherFog), draws rain
// streaks or snowflakes around the camera, and flashes lightning. Rain speeds crops (farming), storms
// send animals to shelter, rain and snow chill the air (survival), and the HUD names the weather.

#include <raylib.h>

class Scene;

enum class WeatherKind { Clear, Cloudy, Rain, Storm, Snow, Fog };

struct WeatherState {
    WeatherKind kind = WeatherKind::Clear;
    float cloud = 0.0f;       // 0..1
    float precipitation = 0.0f; // 0..1 (rain or snow)
    float fog = 0.0f;         // 0..1
    Vector3 wind{};           // body-local m/s along the surface
};

// The weather at a body-local direction on the start planet at a world time (pure; any thread).
WeatherState Weather_At(Scene& world, Vector3 dir, double time);
void Weather_Update(Scene& world, float dt);   // the weather at the player: the sky, particles, lightning
void Weather_Draw3D(Scene& world);             // rain or snow around the camera
void Weather_DrawHUD();                        // the weather's name and lightning flashes
const WeatherState& Weather_Here();            // as of the last update (smoothed)
const char* Weather_Name(WeatherKind k);
// For other systems: rain falling at a body-local direction now (crops), a storm at the player (animals).
bool Weather_RainingAt(Scene& world, Vector3 dir);
bool Weather_StormHere();
// Tests: hold a weather (Clear..Fog) at the player instead of the forecast, or -1 to release it.
void Weather_Force(int kind);

#endif
