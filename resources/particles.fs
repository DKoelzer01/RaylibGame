#version 330

// Particles (#2087): the day light the CPU gives (sun angle x sky, display-encoded as before, so daylight looks unchanged)
// plus the local lights every world shader gets through lampLight(): the headlamp, a held torch, the point lights and the
// light grid (torches, lamps, flares). A puff has no surface, so it is lit as if it faced the eye.
in vec3 fragPosition;
in vec2 fragTexCoord;
in vec4 fragColor;
in float fragDay;

uniform sampler2D texture0;
uniform vec3 viewPos;

#include "shadow.glsl"

out vec4 finalColor;

void main()
{
    vec4 sprite = texture(texture0, fragTexCoord);
    vec3 n = normalize(viewPos - fragPosition);
    vec3 albedo = pow(fragColor.rgb, vec3(2.2));
    vec3 color = albedo * (vec3(pow(fragDay, 2.2)) + lampLight(fragPosition, n));
    finalColor = vec4(pow(max(color, vec3(0.0)), vec3(1.0 / 2.2)) * sprite.rgb, fragColor.a * sprite.a);
}
