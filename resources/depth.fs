#version 330 core
void main() {
    float depth = gl_FragCoord.z;
    gl_FragColor = vec4(depth, depth, depth, 1.0);
    // gl_FragColor = (vec4(gl_FragCoord.z, gl_FragCoord.z, gl_FragCoord.z, 1.0) - vec4(0.5)) * 2.0; // Convert to range [-1, 1]
    // gl_FragColor = vec4(1,0,1,1);
}