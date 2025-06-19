#version 330

// Input vertex attributes (from vertex shader)
in vec3 fragPosition;
in vec2 fragTexCoord;
in vec4 fragColor;
in vec3 fragNormal;
in vec4 fragPosLightSpace; // Position in light space for shadow mapping
in vec3 vertexPos; // Original vertex position for debugging
in float modelZ;

// Input uniform values
uniform sampler2D texture0;
uniform vec4 colDiffuse;
uniform sampler2D shadowMap;

// Output fragment color
out vec4 finalColor;

// NOTE: Add your custom variables here

#define     MAX_LIGHTS              4
#define     LIGHT_DIRECTIONAL       0
#define     LIGHT_POINT             1

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

void main()
{
    // Texel color fetching from texture sampler
    vec4 texelColor = texture(texture0, fragTexCoord);
    vec3 lightDot = vec3(0.0);
    vec3 normal = normalize(fragNormal);
    vec3 viewD = normalize(viewPos - fragPosition);
    vec3 specular = vec3(0.0);
    vec4 tint = colDiffuse * fragColor;
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

    vec3 ambientColor = (ambient.rgb * 0.02) * texelColor.rgb * tint.rgb; // Lower ambient
    vec3 diffuseColor = lightDot * texelColor.rgb * tint.rgb;
    vec3 specularColor = specular * tint.rgb;

    // --- Shadow mapping ---
    vec3 shadowCoords = fragPosLightSpace.xyz;
    vec3 projCoords = fragPosLightSpace.xyz;
    vec3 testProjCoords = projCoords;
    projCoords = projCoords * 0.5 + 0.5;
    // Use the main directional light direction for bias
    vec3 shadowLightDir = -normalize(lights[0].target - lights[0].position);
    float shadow = 1.0;
    float closestDepth = 0.0;
    float currentDepth = 0.0;
    if (projCoords.x >= 0.0 && projCoords.x <= 1.0 && projCoords.y >= 0.0 && projCoords.y <= 1.0) {
        // Try flipping Y for shadow map sample
        vec2 shadowUV = clamp(vec2(projCoords.x,projCoords.y),0.0,1.0);
        closestDepth = texture(shadowMap, shadowUV).r;
        currentDepth = projCoords.z;
        float bias = 0.01;
        shadow = (currentDepth - bias > closestDepth) ? 0.0 : 1.0;
        finalColor = vec4(shadowUV, 0, 1);
    } else {
        finalColor = vec4(0, 0, 1, 1); // Outside shadow map bounds
    }
    // Debug output: R = currentDepth, G = closestDepth, B = shadow factor
    // finalColor = vec4(projCoords, 1);
    // finalColor = vec4(currentDepth,0,0,1);
    finalColor = vec4(closestDepth, 0.0, 0.0, 1.0);
    // finalColor = vec4(texture(texture0, vec2(0.6, 0.6)).r, 0, 0, 1);
    
    // finalColor = vec4(normal * 0.5 + 0.5, 1.0);
    // finalColor = vec4(fragPosition * 0.5 + 0.5, 1.0); // Debug: visualize fragment position
    // finalColor = vec4(testProjCoords, 1.0); // Debug: visualize projected coordinates
    // finalColor = vec4(clamp(projCoords, 0.0, 1.0), 1.0);
    // finalColor = vec4(fragPosLightSpace.w, fragPosLightSpace.w, fragPosLightSpace.w, 1.0); // Debug: visualize light space position
    // finalColor = vec4(fragPosLightSpace.z, fragPosLightSpace.w, 0.0, 1.0);    
    // finalColor = vec4(shadowCoords, 1.0); // Debug: visualize shadow coordinates
    // finalColor = vec4(testProjCoords, 1.0); // No clamping or remapping
    // finalColor = vec4(fragPosLightSpace.x, fragPosLightSpace.y, fragPosLightSpace.z, fragPosLightSpace.w);    
    // finalColor = vec4(0.0, 0.0, modelZ, 1.0);
    // finalColor = vec4(fragPosLightSpace.x * 0.5 + 0.5, fragPosLightSpace.y * 0.5 + 0.5, 0.0, 1.0);
    // finalColor = vec4(vertexPos.z, 0, fragPosLightSpace.z * 0.1 + 0.5, 1.0);
    // if (isnan(projCoords.x) || isinf(projCoords.x)) {
    //     finalColor = vec4(1,0,0,1); // Red for NaN/Inf
    // } else {
    //     finalColor = vec4(projCoords, 1.0);
    // }
    // finalColor = vec4(fragPosLightSpace.w, fragPosLightSpace.w, fragPosLightSpace.w, 1.0);
    // finalColor = vec4(projCoords, 1.0); // Debug: visualize projected coordinates
    // finalColor = vec4(currentDepth, closestDepth, 0.0, 1.0);    
    // finalColor = vec4(fragTexCoord, 0.0,1.0);
    // finalColor = fragPosLightSpace;

    // Only apply shadow to direct lighting, not ambient
//     diffuseColor *= shadow;
//     specularColor *= shadow;
//     vec3 result = ambientColor + diffuseColor + specularColor;
//     finalColor = vec4(result, texel.rgb * tint.a);
//     finalColor = pow(finalColor, vec4(1.0/2.2)); // Gamma correction

//     // Debug: visualize shadow factor
//     finalColor = vec4(vec3(shadow), 1.0);
}
