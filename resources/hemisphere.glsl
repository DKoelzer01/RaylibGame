// Hemisphere ambient (#130): sky light from above, a dimmer warm ground bounce from below, by the surface normal against
// the local up (the planet's, at the camera: the same over a few km). ambientUp = 0: off, the flat ambient (A/B).
uniform vec3 ambientUp;

bool hemiOn() { return dot(ambientUp, ambientUp) > 0.5; }

vec3 hemiAmbient(vec3 amb, vec3 n) {
    if (!hemiOn()) return amb;
    float t = 0.5 + 0.5 * dot(n, ambientUp);
    vec3 sky = amb * vec3(1.20, 1.28, 1.42);
    vec3 ground = amb * vec3(0.62, 0.55, 0.46);
    return mix(ground, sky, t);
}
