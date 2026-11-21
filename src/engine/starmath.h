#ifndef STARMATH_H
#define STARMATH_H

// Plan 13 step 1 groundwork (#1273, design #244): pure helpers for a star that is a body. The star sits at the
// system origin; a planet orbits it on a circle (the same convention as bodies.txt orbits, world.cpp
// PlaceBodyOnOrbit: angle = phase + rate * t in the parent's XZ plane, counter-clockwise seen from +Y); sunlight
// comes from the star's position, and a moon between the star and an observer eclipses it.
// Nothing calls these yet: the star as a body waits for Plan 11 F4 (render origin) and F7 (save format). All doubles:
// star distances are far past float precision (Frames_ToGlobal is double too).

#include "dvec3.h"

struct StarOrbit {
    double radius = 0.0;   // from the star's centre (m)
    double period = 0.0;   // seconds for one orbit (the year), > 0
    double phase = 0.0;    // rad at time 0
};

// The year as game time: `days` days of `dayLength` seconds (Terra: 24 days of 3600 s, Plan 13 decision 2026-11-21).
double Star_YearSeconds(double days, double dayLength);
// The star's GM that makes a circular orbit of `radius` take `period`: GM = (2 pi / T)^2 r^3.
double Star_GMForOrbit(double radius, double period);
// The period of a circular orbit of `radius` around a body of `gm`.
double Star_PeriodForGM(double gm, double radius);
// Position and velocity relative to the star at time t.
DVec3 Star_OrbitPosition(const StarOrbit& o, double t);
DVec3 Star_OrbitVelocity(const StarOrbit& o, double t);
// Unit direction from a body towards the star (Scene::sunDirection points towards the light). {1,0,0} when they coincide.
DVec3 Star_LightDirection(DVec3 bodyPos, DVec3 starPos);
// The share (0..1) of a disc of angular radius `rs` (the star) covered by a disc of angular radius `ro` (the occluder)
// whose centre is `sep` away (all in radians, small-angle discs: the plane circle-circle overlap).
double Star_DiscCover(double rs, double ro, double sep);
// How much of the star an observer sees (1 = all, 0 = total eclipse) with one sphere (a moon or planet) that may stand in
// front of it. An occluder farther than the star, or behind the observer, covers nothing.
double Star_Visibility(DVec3 observer, DVec3 starPos, double starRadius, DVec3 occluderPos, double occluderRadius);

#endif
