#include "joininfo.h"
#include "world.h"
#include "bodies.h"

#include <sstream>

JoinInfo JoinInfo_FromWorld(const Scene& world) {
    JoinInfo j;
    j.seed = worldGen.lastSeed;
    for (const BodyDef& b : worldGen.system.bodies) j.bodies.push_back(FormatBodyDef(b));
    j.fingerprint = WorldFingerprint(world);
    return j;
}

std::string JoinInfo_Write(const JoinInfo& j) {
    std::ostringstream s;
    s << "RAYLIBGAME_JOIN " << j.protocol << "\n";
    s << "seed " << j.seed << "\n";
    s << "fingerprint " << std::hex << j.fingerprint << std::dec << "\n";
    s << "bodies " << j.bodies.size() << "\n";
    for (const std::string& b : j.bodies) s << b << "\n";
    return s.str();
}

bool JoinInfo_Read(const std::string& text, JoinInfo& out, std::string& error) {
    std::istringstream s(text);
    std::string tag;
    JoinInfo j;
    size_t n = 0;
    if (!(s >> tag >> j.protocol) || tag != "RAYLIBGAME_JOIN") { error = "not a join message"; return false; }
    if (j.protocol != JoinInfo::kProtocol) { error = "protocol " + std::to_string(j.protocol) + ", this build speaks " + std::to_string(JoinInfo::kProtocol); return false; }
    if (!(s >> tag >> j.seed) || tag != "seed") { error = "no seed"; return false; }
    if (!(s >> tag >> std::hex >> j.fingerprint >> std::dec) || tag != "fingerprint") { error = "no fingerprint"; return false; }
    if (!(s >> tag >> n) || tag != "bodies") { error = "no bodies"; return false; }
    std::string line;
    std::getline(s, line); // the rest of the count's line
    for (size_t i = 0; i < n; ++i) {
        if (!std::getline(s, line)) { error = "short body list"; return false; }
        BodyDef d;
        std::string e;
        if (!ParseBodyDef(line, d, e)) { error = "body " + std::to_string(i) + ": " + e; return false; }
        j.bodies.push_back(line);
    }
    out = j;
    return true;
}
