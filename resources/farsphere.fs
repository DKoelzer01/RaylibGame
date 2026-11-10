#version 330

// Distance sphere: baked colour (sRGB, alpha 0 = sea) and object-space normals in equirectangular
// maps, looked up per fragment (Tarini's two-UV trick hides the longitude seam), lit by the sun
// with the sphere's own night side, and a soft sun glint on the sea.
in vec3 fragLocal;
in vec3 fragPosition;

uniform sampler2D texture0;  // colour
uniform sampler2D texture2;  // normals (body-local, 0.5 + 0.5 n)
uniform mat4 matModel;
uniform vec3 sunDir;
uniform vec3 viewPos;
uniform vec3 ambient;
uniform vec3 fogColor;
uniform float fogDensity;
uniform vec4 fogShell;      // atmosphere: centre + outer radius (w = 0: fog the whole distance)

out vec4 finalColor;

void main()
{
    vec3 d = normalize(fragLocal);
    float lon = atan(d.z, d.x);
    float lat = asin(clamp(d.y, -1.0, 1.0));
    float u1 = lon / 6.2831853 + 0.5;
    float u2 = fract(u1 + 0.5) - 0.5;
    float v = lat / 3.14159265 + 0.5;
    vec2 uv = vec2(fwidth(u1) <= fwidth(u2) + 1e-5 ? u1 : u2, v);

    vec4 c = texture(texture0, uv);
    vec3 albedo = pow(c.rgb, vec3(2.2));
    bool sea = c.a < 0.5;
    vec3 nl = texture(texture2, uv).xyz * 2.0 - 1.0;
    mat3 rot = mat3(matModel);
    vec3 up = normalize(rot * d);
    vec3 n = normalize(rot * nl);
    // Light: the surface normal, cut off past the sphere's own terminator.
    float day = smoothstep(-0.08, 0.06, dot(up, sunDir));
    float ndl = max(dot(n, sunDir), 0.0) * day;
    vec3 color = albedo * (ambient + vec3(1.0, 0.96, 0.88) * ndl);
    if (sea) {
        vec3 view = normalize(viewPos - fragPosition);
        vec3 h = normalize(sunDir + view);
        color += vec3(1.0, 0.95, 0.85) * pow(max(dot(up, h), 0.0), 250.0) * 0.5 * day;
    }
    // Haze over the part of the view ray inside the air (as the terrain), e.g. the moon by day.
    float dist = length(viewPos - fragPosition), fogDist = dist;
    if (fogShell.w > 0.0) {
        vec3 rd = (fragPosition - viewPos) / max(dist, 1e-3);
        vec3 oc = viewPos - fogShell.xyz;
        float b = dot(oc, rd), cc = dot(oc, oc) - fogShell.w * fogShell.w, disc = b * b - cc;
        float t0 = disc > 0.0 ? -b - sqrt(disc) : 0.0, t1 = disc > 0.0 ? -b + sqrt(disc) : 0.0;
        fogDist = max(0.0, min(t1, fogDist) - max(t0, 0.0));
    }
    color = mix(color, fogColor, clamp(1.0 - exp(-fogDensity * fogDist), 0.0, 1.0));
    finalColor = vec4(pow(color, vec3(1.0 / 2.2)), 1.0);
}
