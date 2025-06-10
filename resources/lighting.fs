#version 330

// Input vertex attributes (from vertex shader)
in vec3 fragPosition;
in vec2 fragTexCoord;
in vec4 fragColor;
in vec3 fragNormal;

// Input uniform values
uniform sampler2D texture0;
uniform vec4 colDiffuse;
uniform sampler2D shadowMap;
uniform mat4 lightSpaceMatrix;

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
    vec4 fragPosLightSpace = lightSpaceMatrix * vec4(fragPosition, 1.0);
    vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
    projCoords = projCoords * 0.5 + 0.5;
    float shadow = 1.0;
    if (projCoords.x >= 0.0 && projCoords.x <= 1.0 && projCoords.y >= 0.0 && projCoords.y <= 1.0) {
        float closestDepth = texture(shadowMap, projCoords.xy).r;
        float currentDepth = projCoords.z;
        float bias = max(0.004 * (1.0 - dot(normal, lightDot)), 0.0003);
        shadow = (currentDepth - bias > closestDepth) ? 0.5 : 1.0;
    }
    diffuseColor *= shadow;
    specularColor *= shadow;

    vec3 result = ambientColor + diffuseColor + specularColor;
    finalColor = vec4(result, texelColor.a * tint.a);
    // finalColor = pow(finalColor, vec4(1.0/2.2)); // Gamma correction
    // Debug: visualize shadow factor
    finalColor = vec4(vec3(shadow), 1.0);
}
