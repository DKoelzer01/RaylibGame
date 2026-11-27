#ifndef DETMATH_H
#define DETMATH_H

// (#1905) Transcendentals that give the same bits on every platform. The C libraries' expf/logf/powf/asinf differ in the
// last bit between glibc and MinGW's CRT, and the erosion bake feeds them back through hundreds of iterations, so a Linux
// server and a Windows client baked different planets from one seed. These use only + - * / (IEEE-exact everywhere with
// -ffp-contract=off, which the engine builds with), sqrt (IEEE-exact) and frexp/ldexp (exact), in double, then round to
// float: accurate to well under a float ulp, and identical everywhere. Inputs outside a function's domain give what the
// bake needs, not NaN (noted per function).
#include <cmath>

namespace det {

// Natural log, x > 0 (x <= 0: -1e30).
inline double Log(double x) {
    if (!(x > 0.0)) return -1e30;
    int e = 0;
    double m = std::frexp(x, &e); // x = m * 2^e, m in [0.5, 1)
    if (m < 0.70710678118654752440) { m *= 2.0; --e; } // m in [sqrt(1/2), sqrt(2))
    // log(m) = 2 atanh(s), s = (m - 1) / (m + 1), |s| < 0.172: the series to s^23 is below 1e-19
    const double s = (m - 1.0) / (m + 1.0), z = s * s;
    double t = 1.0 / 23.0;
    for (int k = 21; k >= 1; k -= 2) t = t * z + 1.0 / k;
    return 2.0 * s * t + e * 0.69314718055994530942;
}

// log2, x > 0 (x <= 0: -1e30); exact for powers of two, so floor(Log2(x)) classes are exact at their edges.
inline double Log2(double x) {
    if (!(x > 0.0)) return -1e30;
    int e = 0;
    double m = std::frexp(x, &e);
    if (m < 0.70710678118654752440) { m *= 2.0; --e; }
    return e + Log(m) * 1.44269504088896340736;
}

// e^x (x < -745: 0; x > 709: 1e308).
inline double Exp(double x) {
    if (x < -745.0) return 0.0;
    if (x > 709.0) return 1e308;
    const double k = std::floor(x * 1.44269504088896340736 + 0.5); // nearest multiple of ln 2
    const double r = (x - k * 6.93147180369123816490e-01) - k * 1.90821492927058770002e-10; // ln 2 in two parts: r exact to ~1e-17
    // e^r, |r| <= 0.35: Taylor to r^16 (the next term is below 1e-21)
    double t = 1.0;
    for (int n = 16; n >= 1; --n) t = 1.0 + t * r / n;
    return std::ldexp(t, (int)k);
}

// x^y for x >= 0 (0^y: 0 for y > 0, 1 for y == 0).
inline double Pow(double x, double y) {
    if (x == 0.0) return y > 0.0 ? 0.0 : 1.0;
    return Exp(y * Log(x));
}

// asin, x clamped to [-1, 1].
inline double Asin(double x) {
    const double a = x < 0.0 ? -x : x, sign = x < 0.0 ? -1.0 : 1.0;
    auto series = [](double v) { // asin(v), |v| <= 0.5: sum c_n v^(2n+1); the term past n = 30 is below 1e-20
        const double z = v * v;
        double term = v, sum = v;
        for (int n = 1; n <= 30; ++n) {
            term *= z * (2.0 * n - 1.0) * (2.0 * n - 1.0) / ((2.0 * n) * (2.0 * n + 1.0));
            sum += term;
        }
        return sum;
    };
    if (a <= 0.5) return sign * series(a);
    const double v = std::sqrt((1.0 - (a > 1.0 ? 1.0 : a)) * 0.5); // asin(a) = pi/2 - 2 asin(sqrt((1 - a) / 2))
    return sign * (1.57079632679489661923 - 2.0 * series(v));
}

namespace detail {
// sin and cos of r, |r| <= pi/4: Taylor to r^21 / r^20 (the next terms are below 1e-20).
inline double SinPoly(double r) {
    const double z = r * r;
    double t = 1.0;
    for (int n = 20; n >= 2; n -= 2) t = 1.0 - t * z / (n * (n + 1.0));
    return r * t;
}
inline double CosPoly(double r) {
    const double z = r * r;
    double t = 1.0;
    for (int n = 20; n >= 2; n -= 2) t = 1.0 - t * z / ((n - 1.0) * n);
    return t;
}
// x = q * pi/2 + r, |r| <= pi/4. pi/2 in two parts (fdlibm's pio2_1, 33 bits, so k * it is exact, and its tail).
inline double Reduce(double x, long long& q) {
    const double k = std::floor(x * 0.63661977236758134308 + 0.5);
    q = (long long)k;
    return (x - k * 1.57079632673412561417e+00) - k * 6.07710050650619224932e-11;
}
} // namespace detail

// sin and cos (accuracy falls off past |x| ~1e6, determinism doesn't).
inline double Sin(double x) {
    long long q = 0;
    const double r = detail::Reduce(x, q);
    switch (q & 3) {
    case 0: return detail::SinPoly(r);
    case 1: return detail::CosPoly(r);
    case 2: return -detail::SinPoly(r);
    default: return -detail::CosPoly(r);
    }
}
inline double Cos(double x) {
    long long q = 0;
    const double r = detail::Reduce(x, q);
    switch (q & 3) {
    case 0: return detail::CosPoly(r);
    case 1: return -detail::SinPoly(r);
    case 2: return -detail::CosPoly(r);
    default: return detail::SinPoly(r);
    }
}

// tan (cos == 0 never happens for a double x: pi/2 isn't representable).
inline double Tan(double x) { return Sin(x) / Cos(x); }

// atan: |x| > 1 folds to pi/2 - atan(1/x), above tan(pi/8) to pi/4 + atan((x - 1) / (x + 1)), then the series on
// |t| <= 0.4143 (t^2 <= 0.1716: the term past n = 26 is below 1e-20).
inline double Atan(double x) {
    const double sign = x < 0.0 ? -1.0 : 1.0;
    double a = x < 0.0 ? -x : x, base = 0.0;
    bool inv = false;
    if (a > 1.0) { a = 1.0 / a; inv = true; }
    if (a > 0.41421356237309504880) { a = (a - 1.0) / (a + 1.0); base = 0.78539816339744830962; }
    const double z = a * a;
    double t = 1.0 / 53.0;
    for (int k = 51; k >= 1; k -= 2) t = 1.0 / k - t * z;
    double r = base + a * t;
    if (inv) r = 1.57079632679489661923 - r;
    return sign * r;
}

inline float Tanf(float x) { return (float)Tan(x); }
inline float Atanf(float x) { return (float)Atan(x); }
inline float Sinf(float x) { return (float)Sin(x); }
inline float Cosf(float x) { return (float)Cos(x); }
inline float Logf(float x) { return (float)Log(x); }
inline float Log2f(float x) { return (float)Log2(x); }
inline float Expf(float x) { return (float)Exp(x); }
inline float Powf(float x, float y) { return (float)Pow(x, y); }
inline float Asinf(float x) { return (float)Asin(x); }

} // namespace det

#endif
