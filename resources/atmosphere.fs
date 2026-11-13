#version 330

// Sky of the planet the camera is near, blended over the star skybox: blue overhead, paler at the
// horizon, dusk colours around a low sun, a sun disc, fading to the stars at night and in space.
in vec3 fragPosition;          // view direction (cube vertex, from skybox.vs)

uniform vec3 sunDir;           // towards the sun
uniform vec3 upDir;            // away from the planet centre, at the camera
uniform float density;
uniform float overcast;        // weather: 0 clear .. 1 a thick grey sky (no sun disc)
// Clouds: a shell over the planet, sampled in the planet's own frame (they stay over the ground
// and drift with the wind). Cover never drops to nothing: a few clouds are always out.
uniform vec3 camLocal;         // camera, planet-local (m)
uniform vec3 toLocalX, toLocalY, toLocalZ; // the planet-local images of the world axes
uniform float cloudShell;      // the cloud layer's radius (m); 0 = no clouds
uniform vec3 cloudDrift;       // wind drift so far (m, planet-local)
uniform float cloudCover;      // 0..1 (weather), with a floor applied here

out vec4 finalColor;

float hash3(vec3 p) {
    p = fract(p * 0.3183099 + 0.1);
    p *= 17.0;
    return fract(p.x * p.y * p.z * (p.x + p.y + p.z));
}
float vnoise(vec3 x) {
    vec3 i = floor(x), f = fract(x);
    f = f * f * (3.0 - 2.0 * f);
    return mix(mix(mix(hash3(i + vec3(0, 0, 0)), hash3(i + vec3(1, 0, 0)), f.x), mix(hash3(i + vec3(0, 1, 0)), hash3(i + vec3(1, 1, 0)), f.x), f.y),
               mix(mix(hash3(i + vec3(0, 0, 1)), hash3(i + vec3(1, 0, 1)), f.x), mix(hash3(i + vec3(0, 1, 1)), hash3(i + vec3(1, 1, 1)), f.x), f.y), f.z);
}
float fbm(vec3 p) {
    float a = 0.5, s = 0.0;
    for (int i = 0; i < 5; i++) { s += a * vnoise(p); p = p * 2.03 + vec3(1.7, 9.2, 3.1); a *= 0.5; }
    return s;
}

void main()
{
    vec3 v = normalize(fragPosition);
    float s = dot(upDir, sunDir);                       // sun height above the horizon
    float day = smoothstep(-0.18, 0.12, s);
    float h = dot(v, upDir);                            // view height above the horizon
    float horizon = pow(1.0 - clamp(h, 0.0, 1.0), 4.0);
    vec3 sky = mix(vec3(0.18, 0.42, 0.85), vec3(0.62, 0.75, 0.92), horizon);
    float mu = dot(v, sunDir);
    float sunset = clamp(1.0 - abs(s) * 4.0, 0.0, 1.0); // sun near the horizon
    sky = mix(sky, vec3(1.0, 0.45, 0.15), sunset * pow(max(mu, 0.0), 3.0) * (0.4 + 0.6 * horizon));
    sky *= day;
    sky = mix(sky, vec3(0.52, 0.55, 0.60) * day, overcast * 0.85);       // cloud cover greys the sky
    sky += vec3(1.0, 0.9, 0.7) * pow(max(mu, 0.0), 300.0) * day * density * (1.0 - overcast); // glare around the sun
    float alpha = density * (0.06 + 0.94 * day);                         // night: the stars show through

    float cloudA = 0.0;
    // Clouds: where the view ray meets the cloud shell (from below; from above, its near side).
    if (cloudShell > 0.0 && density > 0.01) {
        vec3 vl = normalize(toLocalX * v.x + toLocalY * v.y + toLocalZ * v.z);
        float b = dot(camLocal, vl), c = dot(camLocal, camLocal) - cloudShell * cloudShell, disc2 = b * b - c;
        float t = disc2 < 0.0 ? -1.0 : (c < 0.0 ? -b + sqrt(disc2) : -b - sqrt(disc2));
        if (t > 0.0) {
            vec3 p = camLocal + vl * t + cloudDrift;
            float cover = clamp(max(cloudCover, 0.22), 0.0, 1.0);
            float thr = 0.66 - 0.32 * cover;                              // a few clouds at the floor, a full deck at 1
            float n = fbm(p / 2600.0) * 0.8 + vnoise(p / 450.0) * 0.2;
            float a = smoothstep(thr, thr + 0.12, n) * exp(-t / 45000.0) * smoothstep(-0.02, 0.10, dot(vl, normalize(camLocal)));
            float lit = mix(0.72, 1.0, smoothstep(thr, thr + 0.3, n));    // denser cores read brighter on top, greyer below
            vec3 cloud = vec3(lit) * mix(vec3(1.0), vec3(0.62, 0.64, 0.68), overcast);
            cloud = mix(cloud, vec3(1.0, 0.62, 0.42), sunset * 0.45);    // dusk tints them
            cloud *= max(day, 0.04);
            cloudA = a;
            sky = mix(sky, cloud, a * 0.92);
            alpha = max(alpha, a * 0.9 * density);                       // they hide the stars at night too
        }
    }

    // The sun disc (also in space, where the atmosphere is gone).
    float disc = smoothstep(0.99985, 0.99995, mu) * (1.0 - smoothstep(0.4, 0.8, overcast)) * (1.0 - cloudA);
    vec3 col = mix(sky, vec3(1.0, 0.97, 0.9), disc);
    alpha = max(alpha, disc);
    finalColor = vec4(col, alpha); // already display (sRGB-ish) colours, like the skybox
}
