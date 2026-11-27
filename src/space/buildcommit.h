#ifndef BUILDCOMMIT_H
#define BUILDCOMMIT_H

// (#1906) The commit in a build's version string, and whether two builds are the same commit. kGameVersion is
// "<branch> <hash>" with " +local changes" after it for a tree with edits (cmake/game_version.cmake), or whatever
// RAYLIBGAME_VERSION says ("<tag> <hash>": the Pterodactyl egg, tools/build_linux.sh). git abbreviates the hash to as many
// digits as the clone needs to keep it unique (7 in a fresh shallow clone, 8 or more in a big one), so the same commit can
// come out as e5f5f13 on one machine and e5f5f132 on another.

#include <cctype>
#include <string>

// The last word of 7 or more hex digits (the hash); without one, the last word (e.g. "unknown").
inline std::string Build_CommitOf(const std::string& version) {
    std::string last;
    size_t i = 0;
    while (i < version.size()) {
        while (i < version.size() && version[i] == ' ') ++i;
        const size_t a = i;
        while (i < version.size() && version[i] != ' ') ++i;
        if (i == a) break;
        const std::string w = version.substr(a, i - a);
        bool hex = w.size() >= 7;
        for (char c : w) hex = hex && std::isxdigit((unsigned char)c);
        if (hex) return w;
        last = w;
    }
    return last;
}

// Equal, or both hashes (7+ hex digits) with the shorter a prefix of the longer (case-insensitive).
inline bool Build_SameCommit(const std::string& a, const std::string& b) {
    if (a == b) return true;
    auto isHash = [](const std::string& s) { if (s.size() < 7) return false; for (char c : s) if (!std::isxdigit((unsigned char)c)) return false; return true; };
    if (!isHash(a) || !isHash(b)) return false;
    const size_t n = a.size() < b.size() ? a.size() : b.size();
    for (size_t k = 0; k < n; ++k) if (std::tolower((unsigned char)a[k]) != std::tolower((unsigned char)b[k])) return false;
    return true;
}

#endif
