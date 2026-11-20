#ifndef JOININFO_H
#define JOININFO_H

// The join handshake's world description (#316): what a client needs to regenerate the server's world itself - the seed,
// the system definition (bodies.txt lines, as saves store them) and the world fingerprint (world.h WorldFingerprint) to
// compare after generating, so a client built differently stops at join instead of seeing different terrain. Only edits
// are replicated after that. Text, one item per line, so it can travel as one reliable message.

#include <cstdint>
#include <string>
#include <vector>

class Scene;

struct JoinInfo {
    static constexpr uint32_t kProtocol = 1;
    uint32_t protocol = kProtocol;
    uint32_t seed = 0;
    std::vector<std::string> bodies; // FormatBodyDef lines
    uint64_t fingerprint = 0;
};

JoinInfo JoinInfo_FromWorld(const Scene& world); // the server's (or a single-player game's) current world
std::string JoinInfo_Write(const JoinInfo& info);
bool JoinInfo_Read(const std::string& text, JoinInfo& out, std::string& error);

#endif
