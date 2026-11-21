#include "starmath.h"

#include <algorithm>
#include <cmath>

namespace {
    constexpr double kTwoPi = 6.283185307179586476925286766559;
    double Dot(const DVec3& a, const DVec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
    DVec3 Cross(const DVec3& a, const DVec3& b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }
    // Angle between two directions, accurate for tiny angles too (acos of a dot product is not).
    double AngleBetween(const DVec3& a, const DVec3& b) { return atan2(Cross(a, b).length(), Dot(a, b)); }
    double AngularRadius(double radius, double dist) { return dist <= radius ? kTwoPi / 4.0 : asin(radius / dist); }
}

double Star_YearSeconds(double days, double dayLength) { return days * dayLength; }

double Star_GMForOrbit(double radius, double period) {
    const double w = kTwoPi / period;
    return w * w * radius * radius * radius;
}

double Star_PeriodForGM(double gm, double radius) { return kTwoPi * sqrt(radius * radius * radius / gm); }

DVec3 Star_OrbitPosition(const StarOrbit& o, double t) {
    const double a = o.phase + kTwoPi / o.period * t;
    return { cos(a) * o.radius, 0.0, sin(a) * o.radius };
}

DVec3 Star_OrbitVelocity(const StarOrbit& o, double t) {
    const double w = kTwoPi / o.period, a = o.phase + w * t, v = w * o.radius;
    return { -sin(a) * v, 0.0, cos(a) * v };
}

DVec3 Star_LightDirection(DVec3 bodyPos, DVec3 starPos) {
    const DVec3 d = starPos - bodyPos;
    const double l = d.length();
    return l > 0.0 ? d * (1.0 / l) : DVec3{ 1.0, 0.0, 0.0 };
}

double Star_DiscCover(double rs, double ro, double sep) {
    if (rs <= 0.0) return 0.0;
    if (ro <= 0.0 || sep >= rs + ro) return 0.0;
    if (sep <= fabs(rs - ro)) return ro >= rs ? 1.0 : (ro * ro) / (rs * rs); // one disc inside the other
    const double d = sep;
    const double a1 = ro * ro * acos(std::clamp((d * d + ro * ro - rs * rs) / (2.0 * d * ro), -1.0, 1.0));
    const double a2 = rs * rs * acos(std::clamp((d * d + rs * rs - ro * ro) / (2.0 * d * rs), -1.0, 1.0));
    const double k = std::max(0.0, (-d + ro + rs) * (d + ro - rs) * (d - ro + rs) * (d + ro + rs));
    const double area = a1 + a2 - 0.5 * sqrt(k);
    return std::clamp(area / (kTwoPi / 2.0 * rs * rs), 0.0, 1.0);
}

double Star_Visibility(DVec3 observer, DVec3 starPos, double starRadius, DVec3 occluderPos, double occluderRadius) {
    const DVec3 toStar = starPos - observer, toOcc = occluderPos - observer;
    const double ds = toStar.length(), dOcc = toOcc.length();
    if (ds <= 0.0 || dOcc <= 0.0) return 1.0;
    if (Dot(toStar, toOcc) <= 0.0) return 1.0;   // the occluder is behind the observer
    if (dOcc - occluderRadius >= ds) return 1.0;  // the occluder is beyond the star
    const double cover = Star_DiscCover(AngularRadius(starRadius, ds), AngularRadius(occluderRadius, dOcc), AngleBetween(toStar, toOcc));
    return 1.0 - cover;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1273): the star helpers ----
#include "unit_helpers.h"

TEST(StarMath, TerraYearClosesAndMatchesGM) {
    const double year = Star_YearSeconds(24.0, 3600.0); // Plan 13: ~24 game days of one hour
    EXPECT_DOUBLE_EQ(year, 86400.0);
    StarOrbit o{ 4.0e9, year, 0.3 };
    const DVec3 p0 = Star_OrbitPosition(o, 0.0), p1 = Star_OrbitPosition(o, year);
    EXPECT_LT((p1 - p0).length(), 1e-3) << "the orbit closes after one year";
    EXPECT_NEAR(Star_OrbitPosition(o, year * 0.5).length(), o.radius, 1e-3);
    const double gm = Star_GMForOrbit(o.radius, year);
    EXPECT_NEAR(Star_PeriodForGM(gm, o.radius), year, 1e-6 * year);
    // Circular speed = sqrt(GM / r), along the orbit (perpendicular to the radius)
    const DVec3 v = Star_OrbitVelocity(o, 1234.0), p = Star_OrbitPosition(o, 1234.0);
    EXPECT_NEAR(v.length(), sqrt(gm / o.radius), 1e-6 * v.length());
    EXPECT_NEAR(Dot(v, p) / (v.length() * p.length()), 0.0, 1e-12);
}

TEST(StarMath, OrbitFollowsTheBodiesConvention) {
    // PlaceBodyOnOrbit: (cos a, 0, sin a) * r with a = phase + rate * t, velocity (-sin a, 0, cos a) * rate * r.
    StarOrbit o{ 1000.0, 100.0, 0.0 };
    const DVec3 p = Star_OrbitPosition(o, 25.0); // a quarter turn
    EXPECT_NEAR(p.x, 0.0, 1e-9); EXPECT_NEAR(p.z, 1000.0, 1e-9); EXPECT_EQ(p.y, 0.0);
    const DVec3 v = Star_OrbitVelocity(o, 0.0);
    EXPECT_NEAR(v.x, 0.0, 1e-12); EXPECT_GT(v.z, 0.0);
}

TEST(StarMath, LightPointsAtTheStar) {
    const DVec3 d = Star_LightDirection({ 4.0e9, 0.0, 0.0 }, { 0.0, 0.0, 0.0 });
    EXPECT_DOUBLE_EQ(d.x, -1.0); EXPECT_DOUBLE_EQ(d.y, 0.0); EXPECT_DOUBLE_EQ(d.z, 0.0);
    const DVec3 same = Star_LightDirection({ 1, 2, 3 }, { 1, 2, 3 });
    EXPECT_DOUBLE_EQ(same.x, 1.0);
}

TEST(StarMath, DiscCover) {
    EXPECT_DOUBLE_EQ(Star_DiscCover(1.0, 1.0, 2.5), 0.0);         // apart
    EXPECT_DOUBLE_EQ(Star_DiscCover(1.0, 1.2, 0.1), 1.0);         // total
    EXPECT_NEAR(Star_DiscCover(1.0, 0.5, 0.2), 0.25, 1e-12);      // annular: (0.5 / 1)^2
    EXPECT_NEAR(Star_DiscCover(1.0, 1.0, 0.0), 1.0, 1e-12);       // equal discs, centred
    // Equal discs half a radius apart overlap 2 acos(d/2r) - (d/2r) sqrt(1 - (d/2r)^2) times r^2... over pi r^2:
    const double h = 0.25;
    EXPECT_NEAR(Star_DiscCover(1.0, 1.0, 0.5), (2.0 * acos(h) - 2.0 * h * sqrt(1.0 - h * h)) / (kTwoPi / 2.0), 1e-12);
    // Grows as the occluder slides across
    double last = 0.0;
    for (double s = 2.0; s >= 0.0; s -= 0.05) { const double c = Star_DiscCover(1.0, 0.8, s); EXPECT_GE(c, last - 1e-12); last = c; }
}

TEST(StarMath, VisibilityEclipse) {
    // Star at the origin, observer 4e9 m out on +X, a moon between them right on the line: total if it looks bigger.
    const DVec3 star{ 0, 0, 0 }, eye{ 4.0e9, 0, 0 };
    const double starR = 2.0e7;                                   // angular radius 5 mrad from 4e9 m
    EXPECT_NEAR(Star_Visibility(eye, star, starR, { 4.0e9 - 5.0e5, 0, 0 }, 1.0e4), 0.0, 1e-12);  // 20 mrad moon: total
    EXPECT_NEAR(Star_Visibility(eye, star, starR, { 4.0e9 - 5.0e5, 0, 1.0e5 }, 1.0e4), 1.0, 1e-12); // 200 mrad aside: none
    const double partial = Star_Visibility(eye, star, starR, { 4.0e9 - 5.0e5, 0, 1.0e4 }, 1.0e4);   // overlapping edges
    EXPECT_GT(partial, 0.0); EXPECT_LT(partial, 1.0);
    EXPECT_DOUBLE_EQ(Star_Visibility(eye, star, starR, { 4.0e9 + 5.0e5, 0, 0 }, 1.0e4), 1.0);     // behind the observer
    EXPECT_DOUBLE_EQ(Star_Visibility(eye, star, starR, { -1.0e9, 0, 0 }, 1.0e8), 1.0);            // beyond the star
}
#endif
