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
uniform sampler2D debugTexture;

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
    vec3 result = ambientColor + diffuseColor + specularColor;
    finalColor = vec4(result, texelColor.rgb * tint.a);
    finalColor = pow(finalColor, vec4(1.0/2.2)); // Gamma correction
    // finalColor = vec4((fragPosition.z + 1.0) * 0.5, 0.0, 0.0, 1.0);
    // Debug: visualize shadow factor
    // finalColor = vec4(vec3(shadow), 1.0);
    // finalColor = vec4(currentDepth, 0, (1.0-closestDepth)-currentDepth, 1.0); // Debugging depth values
}
