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


// Output fragment color
out vec4 finalColor;
out vec4 dummyColor;

// NOTE: Add your custom variables here

#define     MAX_LIGHTS              4
#define     LIGHT_DIRECTIONAL       0
#define     LIGHT_POINT             1

float texScale = 0.2;

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

float hash(vec3 p) {
    return fract(sin(dot(p, vec3(127.1, 311.7, 74.7))) * 43758.5453);
}

vec2 randomOffset(vec3 pos) {
    float rnd = hash(pos);
    return vec2(rnd, rnd) * 0.005; // Adjust 0.05 as needed
}

vec4 triplanarTextureGrad(sampler2D tex, vec3 pos, vec3 normal, float scale)
{
    vec2 xz = pos.yz * scale;
    vec2 yz = pos.xz * scale;
    vec2 xy = pos.xy * scale;

    vec3 blend = abs(normal) + 0.0001;
    blend = pow(blend, vec3(1.0)); // Softer blend
    blend = max(blend, vec3(0.2)); // Minimum blend to reduce axis dominance
    blend /= (blend.x + blend.y + blend.z);

    float rnd = fract(sin(dot(pos, vec3(127.1, 311.7, 74.7))) * 43758.5453);
    vec2 offset = vec2(rnd, rnd) * 0.01;

    vec4 xProj = textureGrad(tex, yz + offset, dFdx(yz), dFdy(yz));
    vec4 yProj = textureGrad(tex, xz + offset, dFdx(xz), dFdy(xz));
    vec4 zProj = textureGrad(tex, xy + offset, dFdx(xy), dFdy(xy));

    return xProj * blend.x + yProj * blend.y + zProj * blend.z;
}

vec3 triplanarNormal(sampler2D tex, vec3 pos, vec3 normal, float scale)
{
    // Projected UVs
    vec2 xz = pos.yz * scale;
    vec2 yz = pos.xz * scale;
    vec2 xy = pos.xy * scale;

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
    vec4 diffuse    = triplanarTextureGrad(colorTexture, fragPosition, fragNormal, texScale);
    vec3 normal     = triplanarNormal(normalTexture, fragPosition, fragNormal, texScale);
    float height    = triplanarTextureGrad(heightTexture, fragPosition, fragNormal, texScale).r;
    float rough     = triplanarTextureGrad(roughnessTexture, fragPosition, fragNormal, texScale).r;
    float ao        = triplanarTextureGrad(aoTexture, fragPosition, fragNormal, texScale).r;

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
    // Debug: visualize shadow 
    // vec3 blendedNormal = triplanarNormal(normalTexture, fragPosition, fragNormal, texScale);
    finalColor = vec4(texture(shadowMap, shadowUV).rgb, 1.0); // Debugging shadow map texture
    // finalColor = diffuse;
    // finalColor = vec4(normal, 1.0);
    // finalColor = vec4(vec3(shadow), 1.0);
    // finalColor = vec4(currentDepth, 0, (1.0-closestDepth)-currentDepth, 1.0); // Debugging depth values
}
