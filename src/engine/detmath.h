#ifndef DETMATH_H
#define DETMATH_H

// (#1905) Transcendentals that give the same bits on every platform. The C libraries' expf/logf/powf/asinf differ in the
// last bit between glibc and MinGW's CRT, and the erosion bake feeds them back through hundreds of iterations, so a Linux
// server and a Windows client baked different planets from one seed. These use only + - * / (IEEE-exact everywhere with
// -ffp-contract=off, which the engine builds with), sqrt (IEEE-exact) and frexp/ldexp (exact), in double, then round to
// float: accurate to well under a float ulp, and identical everywhere. Inputs outside a function's domain give what the
// bake needs, not NaN (noted per function).
#include <cmath>
#include <cstdint>
#include <cstring>

namespace det {

namespace detail {
// x = m * 2^e with m in [sqrt(1/2), sqrt(2)), x > 0 finite. From the bits for normal x (a library frexp call cost a
// tenth of the bake's flood pass, #1905); subnormals through std::frexp (exact as well).
inline double SplitLog(double x, int& e) {
    uint64_t b; std::memcpy(&b, &x, 8);
    const int be = (int)((b >> 52) & 0x7ff);
    double m;
    if (be == 0) { m = std::frexp(x, &e); m *= 2.0; --e; } // subnormal: m in [1, 2)
    else { e = be - 1023; b = (b & 0x000fffffffffffffull) | 0x3ff0000000000000ull; std::memcpy(&m, &b, 8); } // m in [1, 2)
    if (m >= 1.41421356237309504880) { m *= 0.5; ++e; }
    return m;
}
// 2^k as a double, k in -1022..1023 (exact).
inline double Pow2(int k) { const uint64_t b = (uint64_t)(k + 1023) << 52; double v; std::memcpy(&v, &b, 8); return v; }
// log(m), m in [sqrt(1/2), sqrt(2)): 2 atanh(s), s = (m - 1) / (m + 1), |s| < 0.172. The series to s^13; the rest is
// below 2e-12 relative, far under a float ulp (6e-8). The coefficients are constants the compiler folds exactly.
inline double LogMantissa(double m) {
    const double s = (m - 1.0) / (m + 1.0), z = s * s;
    const double t = 1.0 + z * (1.0 / 3.0 + z * (1.0 / 5.0 + z * (1.0 / 7.0 + z * (1.0 / 9.0 + z * (1.0 / 11.0 + z * (1.0 / 13.0))))));
    return 2.0 * s * t;
}
} // namespace detail

// Natural log, x > 0 (x <= 0: -1e30).
inline double Log(double x) {
    if (!(x > 0.0)) return -1e30;
    int e = 0;
    const double m = detail::SplitLog(x, e);
    return detail::LogMantissa(m) + e * 0.69314718055994530942;
}

// log2, x > 0 (x <= 0: -1e30); exact for powers of two, so floor(Log2(x)) classes are exact at their edges.
inline double Log2(double x) {
    if (!(x > 0.0)) return -1e30;
    int e = 0;
    const double m = detail::SplitLog(x, e);
    return e + detail::LogMantissa(m) * 1.44269504088896340736;
}

// e^x (x < -708: 0; x > 709: 1e308).
inline double Exp(double x) {
    if (x < -708.0) return 0.0;
    if (x > 709.0) return 1e308;
    const double k = std::floor(x * 1.44269504088896340736 + 0.5); // nearest multiple of ln 2
    const double r = (x - k * 6.93147180369123816490e-01) - k * 1.90821492927058770002e-10; // ln 2 in two parts: r exact to ~1e-17
    // e^r, |r| <= 0.347: Taylor to r^11 (the rest is below 1e-14)
    const double t = 1.0 + r * (1.0 + r * (1.0 / 2.0 + r * (1.0 / 6.0 + r * (1.0 / 24.0 + r * (1.0 / 120.0 + r * (1.0 / 720.0 + r * (1.0 / 5040.0
                     + r * (1.0 / 40320.0 + r * (1.0 / 362880.0 + r * (1.0 / 3628800.0 + r * (1.0 / 39916800.0)))))))))));
    const int ki = (int)k;
    return ki >= -1022 ? t * detail::Pow2(ki) : std::ldexp(t, ki); // the product is exact (a power of two)
}

// x^y for x >= 0 (0^y: 0 for y > 0, 1 for y == 0).
inline double Pow(double x, double y) {
    if (x == 0.0) return y > 0.0 ? 0.0 : 1.0;
    return Exp(y * Log(x));
}

// asin, x clamped to [-1, 1].
inline double Asin(double x) {
    const double a = x < 0.0 ? -x : x, sign = x < 0.0 ? -1.0 : 1.0;
    // asin(v), |v| <= 0.5: sum c_n v^(2n+1), c_n = (2n)! / (4^n n!^2 (2n + 1)); the rest past n = 20 is below 1e-14.
    // The table is built at compile time (the compiler's double arithmetic is IEEE too, so every build has the same).
    struct Coef { double c[21]; };
    static constexpr Coef kC = [] {
        Coef t{};
        double b = 1.0;
        for (int n = 0; n <= 20; ++n) { if (n > 0) b *= (2.0 * n - 1.0) / (2.0 * n); t.c[n] = b / (2.0 * n + 1.0); }
        return t;
    }();
    auto series = [](double v) {
        const double z = v * v;
        double t = kC.c[20];
        for (int n = 19; n >= 0; --n) t = t * z + kC.c[n];
        return v * t;
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
inline float Root4f(float x) { return (float)std::sqrt(std::sqrt((double)x)); } // x^0.25, x >= 0: sqrt is exact IEEE, and cheaper than Pow
inline float Asinf(float x) { return (float)Asin(x); }

} // namespace det

#endif
