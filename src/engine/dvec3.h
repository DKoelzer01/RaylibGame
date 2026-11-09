#ifndef DVEC3_H
#define DVEC3_H

#include <raylib.h>
#include <cmath>

// Double-precision world position. Floats lose sub-metre precision ~10^4-10^5 units from the
// origin; authoritative positions are kept in doubles and converted to small camera-relative
// floats for rendering (floating origin, see Scene::renderOrigin).
struct DVec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    DVec3() = default;
    constexpr DVec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    explicit DVec3(Vector3 v) : x(v.x), y(v.y), z(v.z) {}

    DVec3 operator+(const DVec3& o) const { return { x + o.x, y + o.y, z + o.z }; }
    DVec3 operator-(const DVec3& o) const { return { x - o.x, y - o.y, z - o.z }; }
    DVec3 operator*(double s) const { return { x * s, y * s, z * s }; }
    DVec3& operator+=(const DVec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    DVec3& operator-=(const DVec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    double length() const { return std::sqrt(x * x + y * y + z * z); }
    Vector3 toFloat() const { return { (float)x, (float)y, (float)z }; }
};

inline DVec3 operator+(const DVec3& a, Vector3 b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }

#endif
