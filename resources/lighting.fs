#version 330

// Input vertex attributes (from vertex shader)
in vec3 fragPosition;
in vec2 fragTexCoord;
in vec4 fragColor;
in vec3 fragNormal;
in vec4 fragPosLightSpace; // Position in light space for shadow mapping

// Input uniform values
uniform sampler2D shadowMap;
uniform sampler2D colorTexture;
uniform sampler2D normalTexture;
uniform sampler2D heightTexture;
uniform sampler2D roughnessTexture;
uniform sampler2D aoTexture;

uniform vec3 objectCenter; // Center of the asteroid/object

// Output fragment color
out vec4 finalColor;
out vec4 dummyColor;

// NOTE: Add your custom variables here

#define     MAX_LIGHTS              4
#define     LIGHT_DIRECTIONAL       0
#define     LIGHT_POINT             1

float texScale = 0.01; // Scale for texture mapping

struct Light {
    int enabled;
    int type;
    vec3 position;
    vec3 target;
    vec4 color;
};

// Input lighting values
uniform Light lights[MAX_LIGHTS];
uniform vec4 ambient;
uniform vec3 viewPos;


vec4 triplanarTextureGrad(sampler2D tex, vec3 pos, vec3 normal, float scale)
{
    // Decorrelate axes: rotate and offset each projection
    mat2 rotA = mat2(cos(0.0), -sin(0.0), sin(0.0), cos(0.0));
    mat2 rotB = mat2(cos(1.0), -sin(1.0), sin(1.0), cos(1.0));
    mat2 rotC = mat2(cos(2.0), -sin(2.0), sin(2.0), cos(2.0));

    vec2 xz = rotA * (pos.yz * scale + vec2(13.1, 7.7));
    vec2 yz = rotB * (pos.xz * scale + vec2(-5.3, 2.2));
    vec2 xy = rotC * (pos.xy * scale + vec2(8.8, -11.4));

    vec3 blend = abs(normal) + 0.0001;
    blend = pow(blend, vec3(4.0)); // Softer blend
    blend = max(blend, vec3(0.05));
    blend /= (blend.x + blend.y + blend.z);

    vec3 scaledPos = pos * scale;
    vec3 dx = dFdx(scaledPos);
    vec3 dy = dFdy(scaledPos);

    vec4 xProj = textureGrad(tex, xz, dx.yz, dy.yz);
    vec4 yProj = textureGrad(tex, yz, dx.xz, dy.xz);
    vec4 zProj = textureGrad(tex, xy, dx.xy, dy.xy);

    return xProj * blend.x + yProj * blend.y + zProj * blend.z;
}

vec3 triplanarNormal(sampler2D tex, vec3 pos, vec3 normal, float scale)
{
    // Decorrelate axes: rotate and offset each projection
    mat2 rotA = mat2(cos(0.0), -sin(0.0), sin(0.0), cos(0.0));
    mat2 rotB = mat2(cos(1.0), -sin(1.0), sin(1.0), cos(1.0));
    mat2 rotC = mat2(cos(2.0), -sin(2.0), sin(2.0), cos(2.0));

    vec2 xz = rotA * (pos.yz * scale + vec2(13.1, 7.7));
    vec2 yz = rotB * (pos.xz * scale + vec2(-5.3, 2.2));
    vec2 xy = rotC * (pos.xy * scale + vec2(8.8, -11.4));

    // Blend weights
    vec3 blend = abs(normal) + 0.0001;
    blend = pow(blend, vec3(1.0));
    blend = max(blend, vec3(0.2));
    blend /= (blend.x + blend.y + blend.z);

    // Sample tangent-space normals from each axis
    vec3 nX = texture(tex, yz).xyz * 2.0 - 1.0;
    vec3 nY = texture(tex, xz).xyz * 2.0 - 1.0;
    vec3 nZ = texture(tex, xy).xyz * 2.0 - 1.0;

    // Transform to world space
    nX = vec3(nX.z, nX.y, -nX.x); // X projection
    nY = vec3(nY.x, nY.z, -nY.y); // Y projection
    nZ = vec3(nZ.x, nZ.y, nZ.z);  // Z projection

    // Blend and normalize
    vec3 blended = normalize(nX * blend.x + nY * blend.y + nZ * blend.z);
    return blended;
}

void main()
{
    vec4 diffuse;
    vec3 normal;
    float height, rough, ao;

    diffuse    = triplanarTextureGrad(colorTexture, fragPosition, fragNormal, texScale);
    normal     = triplanarNormal(normalTexture, fragPosition, fragNormal, texScale);
    height     = triplanarTextureGrad(heightTexture, fragPosition, fragNormal, texScale).r;
    rough      = triplanarTextureGrad(roughnessTexture, fragPosition, fragNormal, texScale).r;
    ao         = triplanarTextureGrad(aoTexture, fragPosition, fragNormal, texScale).r;

    vec3 lightDot = vec3(0.0);
    vec3 viewD = normalize(viewPos - fragPosition);
    vec3 specular = vec3(0.0);
    for (int i = 0; i < MAX_LIGHTS; i++)
    {
        if (lights[i].enabled == 1)
        {
            vec3 light = vec3(0.0);
            if (lights[i].type == LIGHT_DIRECTIONAL)
            {
                light = -normalize(lights[i].target - lights[i].position);
            }
            if (lights[i].type == LIGHT_POINT)
            {
                light = normalize(lights[i].position - fragPosition);
            }
            float NdotL = max(dot(normal, light), 0.0);
            lightDot += lights[i].color.rgb*NdotL;
            float specCo = 0.0;
            if (NdotL > 0.0) specCo = pow(max(0.0, dot(viewD, reflect(-(light), normal))), 16.0); // 16 refers to shine
            specular += specCo;
        }
    }

    vec3 ambientColor = (ambient.rgb * 0.02) * diffuse.rgb; // Lower ambient
    vec3 diffuseColor = lightDot * diffuse.rgb;
    vec3 specularColor = specular;

    // --- Shadow mapping ---
    vec3 projCoords = fragPosLightSpace.xyz;
    // Use the main directional light direction for bias
    vec3 shadowLightDir = -normalize(lights[0].target - lights[0].position);
    float shadow = 1.0;
    float closestDepth = 0.0;
    float currentDepth = 0.0;
    vec2 shadowUV = projCoords.xy * 0.5 + 0.5; // Convert from [-1, 1] to [0, 1]
    shadowUV.y = 1.0 - shadowUV.y; // Flip Y coordinate for OpenGL texture coordinates
    float bias = max(0.05 * (1.0 - dot(normal, shadowLightDir)), 0.0001);
    if (shadowUV.x >= 0.0 && shadowUV.x <= 1.0 && shadowUV.y >= 0.0 && shadowUV.y <= 1.0) {
        closestDepth = texture(shadowMap, shadowUV).r; 
        currentDepth = projCoords.z * 0.5 + 0.5; // Convert from [-1, 1] to [0, 1]
        shadow = (currentDepth - bias < closestDepth) ? 0.0 : 1.0;
    } else {
        finalColor = vec4(0, 0, 1, 1); // Outside shadow map bounds
    }

    // Only apply shadow to direct lighting, not ambient
    diffuseColor *= shadow;
    specularColor *= shadow;
    vec4 result = vec4(ambientColor + diffuseColor + specularColor, 1.0);
    finalColor = pow(result, vec4(1.0/2.2)); // Gamma correction
    dummyColor = vec4(height, rough, ao, result); // Output additional data for debugging
}
