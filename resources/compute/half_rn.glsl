// Floats to IEEE half precision, nearest even, as the CPU kernel's H() (fluidkernel.cpp) rounds them (#2116). GLSL leaves
// packHalf2x16's rounding to the driver: NVIDIA rounds toward zero, so the face flows the GPU stored were up to an ulp
// below the CPU's and the two fluid paths parted on a spring's first substep. halfRN rounds first; packHalf2x16 of the
// result is then exact on any driver (the value is a half).
float halfRN(float f) {
    if (abs(f) < 6.103515625e-05) return roundEven(f * 16777216.0) / 16777216.0; // half subnormals: 2^-24 steps
    if (abs(f) > 65504.0) return f > 0.0 ? uintBitsToFloat(0x7F800000u) : uintBitsToFloat(0xFF800000u);
    uint x = floatBitsToUint(f);
    x += 0x0FFFu + ((x >> 13) & 1u); // 23 -> 10 mantissa bits
    return uintBitsToFloat(x & ~0x1FFFu);
}
uint packHalfRN(vec2 v) { return packHalf2x16(vec2(halfRN(v.x), halfRN(v.y))); }
