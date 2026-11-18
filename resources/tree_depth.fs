#version 330
// Trees and bushes in the shadow pass (trees.vs + this). A caster fades out with the same distance fade as trees.fs
// (#541: bushes dithered nearly invisible at 225-300 m still cast full shadows, dark spots on the mid-range ground with
// no bush over them). The fade is dithered over the shadow map's texels; the lookup's PCF averages it into a lighter
// shadow.
in vec3 fragPosition;
uniform vec3 viewPos;   // the main camera (render space), not the light's
uniform float fadeEnd;  // the draw distance of this variant (0: no fade)
void main()
{
    if (fadeEnd > 0.0) {
        float dist = length(viewPos - fragPosition);
        float fade = clamp((fadeEnd - dist) / (fadeEnd * 0.25), 0.0, 1.0);
        fade = fade * fade * (3.0 - 2.0 * fade);
        vec2 cell = floor(gl_FragCoord.xy);
        if (fract(sin(dot(cell, vec2(12.9898, 78.233))) * 43758.5453) > fade) discard;
    }
}
