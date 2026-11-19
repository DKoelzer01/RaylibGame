// The wind field (#650, src/engine/wind.*): the gust texture (512 m tile: rolling waves, burst cells, direction swing)
// in the plane at an anchor near the player, sampled where the air now at p was when the pattern started (minus the
// drift), so gust fronts roll downwind. Same formula as Wind_At on the CPU. windOn 0: no wind (RAYLIBGAME_WIND=0).
uniform sampler2D windTex;
uniform float windOn;
uniform vec3 windOrigin;  // render space
uniform vec3 windE1;      // render-space tangent axes over the body's scale (dot gives body metres)
uniform vec3 windE2;
uniform vec2 windShift;   // metres the air has moved along e1/e2 (wrapped to the period)
uniform vec2 windMean;    // mean wind along e1/e2 (m/s): weather, daily cycle and storms in
uniform float windStorm;  // 0 calm .. 1 storm
uniform float windPeriod; // the texture's span (m)

// xyz: the wind at p (render space, m/s), w: how much of a burst is passing (0..1).
vec4 windAt(vec3 p)
{
    if (windOn < 0.5) return vec4(0.0);
    vec3 r = p - windOrigin;
    vec2 uv = (vec2(dot(r, windE1), dot(r, windE2)) - windShift) / windPeriod;
    vec3 s = texture(windTex, uv).rgb;
    float waves = s.r * 2.0 - 1.0;
    float thr = 0.80 - 0.14 * windStorm;
    float burst = smoothstep(thr, thr + 0.1, s.g);
    float factor = max(0.15, 1.0 + 0.5 * waves + 0.8 * burst); // (#888, as wind.cpp)
    float swing = 0.35 * (s.b * 2.0 - 1.0);
    float ms = length(windMean);
    if (ms < 1e-4) return vec4(0.0, 0.0, 0.0, burst);
    vec2 d = windMean / ms;
    vec2 dr = vec2(d.x * cos(swing) - d.y * sin(swing), d.x * sin(swing) + d.y * cos(swing));
    vec3 dir = normalize(windE1) * dr.x + normalize(windE2) * dr.y;
    return vec4(dir * ms * factor, burst);
}
