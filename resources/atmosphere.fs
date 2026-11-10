#version 330

// Sky of the planet the camera is near, blended over the star skybox: blue overhead, paler at the
// horizon, dusk colours around a low sun, a sun disc, fading to the stars at night and in space.
in vec3 fragPosition;          // view direction (cube vertex, from skybox.vs)

uniform vec3 sunDir;           // towards the sun
uniform vec3 upDir;            // away from the planet centre, at the camera
uniform float density;         // 1 at sea level .. 0 above the atmosphere

out vec4 finalColor;

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
    sky += vec3(1.0, 0.9, 0.7) * pow(max(mu, 0.0), 300.0) * day * density; // glare around the sun
    float alpha = density * (0.06 + 0.94 * day);                         // night: the stars show through

    // The sun disc (also in space, where the atmosphere is gone).
    float disc = smoothstep(0.99985, 0.99995, mu);
    vec3 col = mix(sky, vec3(1.0, 0.97, 0.9), disc);
    alpha = max(alpha, disc);
    finalColor = vec4(col, alpha); // already display (sRGB-ish) colours, like the skybox
}
