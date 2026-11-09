#include "bodies.h"
#include "utils/logger.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

const char* BodyKindName(BodyKind kind) {
    switch (kind) {
        case BodyKind::Asteroid: return "asteroid";
        case BodyKind::Planet: return "planet";
        case BodyKind::Moon: return "moon";
        case BodyKind::Ring: return "ring";
    }
    return "?";
}

float AutoSurfaceGravity(BodyKind kind, float radius) {
    if (kind == BodyKind::Planet || kind == BodyKind::Moon) return 9.81f * radius / 100000.0f; // Earth-like at 100 km
    return std::min(9.81f, 3.0f + radius * 0.03f); // asteroids: walkable even when small
}

namespace {
    bool ParseFloat(const std::string& s, float& out) {
        char* end = nullptr;
        out = strtof(s.c_str(), &end);
        return end && end != s.c_str() && *end == '\0';
    }

    bool ParseFloats(const std::string& s, float* out, int n) {
        std::stringstream ss(s);
        std::string part;
        int i = 0;
        while (std::getline(ss, part, ',')) {
            if (i >= n || !ParseFloat(part, out[i])) return false;
            ++i;
        }
        return i == n;
    }

    std::string NumD(double v) {
        char buf[40];
        snprintf(buf, sizeof(buf), "%.15g", v);              // shortest that usually round-trips
        if (strtod(buf, nullptr) != v) snprintf(buf, sizeof(buf), "%.17g", v);
        return buf;
    }

    std::string Num(float v) {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.7g", v);               // shortest that usually round-trips
        if (strtof(buf, nullptr) != v) snprintf(buf, sizeof(buf), "%.9g", v);
        return buf;
    }
}

bool ParseBodyDef(const std::string& line, BodyDef& out, std::string& error) {
    std::stringstream ss(line);
    std::string kind;
    BodyDef d;
    if (!(ss >> kind >> d.name)) { error = "expected: kind name key=value ..."; return false; }
    if (kind == "asteroid") d.kind = BodyKind::Asteroid;
    else if (kind == "planet") d.kind = BodyKind::Planet;
    else if (kind == "moon") d.kind = BodyKind::Moon;
    else if (kind == "ring") d.kind = BodyKind::Ring;
    else { error = "unknown kind '" + kind + "' (asteroid, planet, moon, ring)"; return false; }

    std::string tok;
    while (ss >> tok) {
        const size_t eq = tok.find('=');
        if (eq == std::string::npos || eq == 0) { error = "expected key=value, got '" + tok + "'"; return false; }
        const std::string key = tok.substr(0, eq), val = tok.substr(eq + 1);
        bool ok = true;
        if (key == "radius") ok = ParseFloat(val, d.radius) && d.radius > 0.0f;
        else if (key == "position") { ok = ParseFloats(val, &d.position.x, 3); d.hasPosition = ok; }
        else if (key == "orbit") {
            const size_t c = val.find(',');
            float rest[2];
            ok = c != std::string::npos && c > 0 && ParseFloats(val.substr(c + 1), rest, 2);
            if (ok) { d.orbitParent = val.substr(0, c); d.orbitAltitude = rest[0]; d.orbitPhaseDeg = rest[1]; }
        }
        else if (key == "spin") { ok = ParseFloat(val, d.spinDegPerSec); d.spinFromDay = false; }
        else if (key == "day") {
            ok = ParseFloat(val, d.dayLength) && d.dayLength != 0.0f;
            if (ok) { d.spinDegPerSec = 360.0f / d.dayLength; d.spinFromDay = true; }
        }
        else if (key == "gravity") {
            if (val == "auto") d.gravity = -1.0f;
            else ok = ParseFloat(val, d.gravity) && d.gravity >= 0.0f;
        }
        else if (key == "seed") {
            char* end = nullptr;
            const unsigned long s = strtoul(val.c_str(), &end, 10);
            ok = end && *end == '\0' && !val.empty();
            d.seed = (unsigned int)s;
            d.hasSeed = ok;
        }
        else if (key == "count") { float c; ok = ParseFloat(val, c) && c >= 0.0f; d.count = (int)c; }
        else if (key == "size") { float s[2]; ok = ParseFloats(val, s, 2) && s[0] > 0.0f && s[0] <= s[1]; if (ok) { d.sizeMin = s[0]; d.sizeMax = s[1]; } }
        else if (key == "extent") ok = ParseFloat(val, d.extent) && d.extent >= 0.0f;
        else if (key == "flatten") {
            char* end = nullptr;
            d.flatten = strtod(val.c_str(), &end);
            ok = end && end != val.c_str() && *end == '\0' && d.flatten >= 0.0;
        }
        else if (key == "maxspin") ok = ParseFloat(val, d.maxSpinDegPerSec) && d.maxSpinDegPerSec >= 0.0f;
        else { error = "unknown key '" + key + "'"; return false; }
        if (!ok) { error = "bad value for '" + key + "': '" + val + "'"; return false; }
    }
    if (d.kind != BodyKind::Ring && d.hasPosition && !d.orbitParent.empty()) {
        error = "give either position= or orbit=, not both";
        return false;
    }
    out = d;
    return true;
}

std::string FormatBodyDef(const BodyDef& d) {
    std::string s = std::string(BodyKindName(d.kind)) + " " + d.name + " radius=" + Num(d.radius);
    if (d.kind == BodyKind::Ring) {
        s += " count=" + std::to_string(d.count) + " size=" + Num(d.sizeMin) + "," + Num(d.sizeMax) +
             " extent=" + Num(d.extent) + " flatten=" + NumD(d.flatten) + " maxspin=" + Num(d.maxSpinDegPerSec);
    } else {
        if (!d.orbitParent.empty()) s += " orbit=" + d.orbitParent + "," + Num(d.orbitAltitude) + "," + Num(d.orbitPhaseDeg);
        else s += " position=" + Num(d.position.x) + "," + Num(d.position.y) + "," + Num(d.position.z);
        s += d.spinFromDay ? " day=" + Num(d.dayLength) : " spin=" + Num(d.spinDegPerSec);
    }
    if (d.gravity >= 0.0f) s += " gravity=" + Num(d.gravity);
    if (d.hasSeed) s += " seed=" + std::to_string(d.seed);
    return s;
}

bool LoadSystemDef(const std::string& path, SystemDef& out) {
    std::ifstream f(path);
    if (!f.is_open()) { logger.logf("[Bodies] cannot open %s\n", path.c_str()); return false; }
    SystemDef sys;
    std::string line;
    int lineNo = 0;
    bool ok = true;
    while (std::getline(f, line)) {
        ++lineNo;
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        BodyDef d;
        std::string error;
        if (!ParseBodyDef(line, d, error)) {
            logger.logf("[Bodies] %s:%d: %s\n", path.c_str(), lineNo, error.c_str());
            ok = false;
            continue;
        }
        for (const BodyDef& prev : sys.bodies)
            if (prev.name == d.name) { logger.logf("[Bodies] %s:%d: duplicate name '%s'\n", path.c_str(), lineNo, d.name.c_str()); ok = false; }
        if (!d.orbitParent.empty()) {
            const bool found = std::any_of(sys.bodies.begin(), sys.bodies.end(),
                                           [&](const BodyDef& b) { return b.name == d.orbitParent && b.kind != BodyKind::Ring; });
            if (!found) { logger.logf("[Bodies] %s:%d: orbit parent '%s' must be a body defined above\n", path.c_str(), lineNo, d.orbitParent.c_str()); ok = false; }
        }
        sys.bodies.push_back(d);
    }
    if (sys.bodies.empty()) { logger.logf("[Bodies] %s defines no bodies\n", path.c_str()); ok = false; }
    if (ok) out = sys;
    return ok;
}

SystemDef LegacySystemDef(int count, int size, float systemRadius, float maxSpin, Vector3 first) {
    SystemDef sys;
    BodyDef start;
    start.kind = BodyKind::Asteroid;
    start.name = "planetoid_0";
    start.radius = (float)size;
    start.hasPosition = true;
    start.position = first;
    start.spinDegPerSec = 0.0f; // the start body holds still (easier to test edits)
    sys.bodies.push_back(start);
    if (count > 1) {
        BodyDef ring;
        ring.kind = BodyKind::Ring;
        ring.name = "planetoid";
        ring.radius = (float)size;
        ring.count = count - 1;
        ring.extent = systemRadius;
        ring.maxSpinDegPerSec = maxSpin;
        sys.bodies.push_back(ring);
    }
    return sys;
}
