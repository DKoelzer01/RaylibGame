#ifndef FRUSTUM_H
#define FRUSTUM_H

#include <raylib.h>
#include <raymath.h>

// View frustum as 6 planes (a,b,c,d) with normals pointing inward: dot(n, p) + d >= 0 is inside.
struct Frustum {
    Vector4 planes[6];

    // viewProj in raymath order: MatrixMultiply(view, proj).
    static Frustum fromViewProj(const Matrix& m) {
        // raymath stores matrices so Vector3Transform computes x' = m0*x + m4*y + m8*z + m12,
        // i.e. "row i" of the math matrix is (m[i], m[i+4], m[i+8], m[i+12]).
        const Vector4 r0 = { m.m0, m.m4, m.m8,  m.m12 };
        const Vector4 r1 = { m.m1, m.m5, m.m9,  m.m13 };
        const Vector4 r2 = { m.m2, m.m6, m.m10, m.m14 };
        const Vector4 r3 = { m.m3, m.m7, m.m11, m.m15 };
        Frustum f;
        f.planes[0] = add(r3, r0);  // left
        f.planes[1] = sub(r3, r0);  // right
        f.planes[2] = add(r3, r1);  // bottom
        f.planes[3] = sub(r3, r1);  // top
        f.planes[4] = add(r3, r2);  // near
        f.planes[5] = sub(r3, r2);  // far
        for (Vector4& p : f.planes) {
            const float len = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
            if (len > 0.0f) { p.x /= len; p.y /= len; p.z /= len; p.w /= len; }
        }
        return f;
    }

    bool sphereVisible(Vector3 c, float r) const {
        for (const Vector4& p : planes) {
            if (p.x * c.x + p.y * c.y + p.z * c.z + p.w < -r) return false;
        }
        return true;
    }

private:
    static Vector4 add(Vector4 a, Vector4 b) { return { a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w }; }
    static Vector4 sub(Vector4 a, Vector4 b) { return { a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w }; }
};

#endif
