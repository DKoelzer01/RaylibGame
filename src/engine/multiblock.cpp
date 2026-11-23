#include "multiblock.h"
#include "plots.h"
#include "utils/logger.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <set>
#include <sstream>

// Nothing here remembers which structure a part belongs to: the validator only reads cells through the lookup, so two
// structures may share a wall of casings.

namespace {
std::string Trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r')) ++a;
    while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r')) --b;
    return s.substr(a, b - a);
}
std::vector<std::string> Split(const std::string& s, char c) {
    std::vector<std::string> out;
    size_t start = 0;
    for (size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == c) { out.push_back(Trim(s.substr(start, i - start))); start = i + 1; }
    return out;
}
std::string ReadFile(const char* path) {
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}
const char* const kHatchKinds[] = { "item_in", "item_out", "fluid_in", "fluid_out", "energy", "rotation" };
bool IsHatchKind(const std::string& k) {
    for (const char* h : kHatchKinds) if (k == h) return true;
    return false;
}
// "a..b"
bool ParseRange(const std::string& s, int& a, int& b) {
    const size_t p = s.find("..");
    if (p == std::string::npos) return false;
    const std::string l = Trim(s.substr(0, p)), r = Trim(s.substr(p + 2));
    if (l.empty() || r.empty()) return false;
    char* e1 = nullptr; char* e2 = nullptr;
    a = (int)strtol(l.c_str(), &e1, 10);
    b = (int)strtol(r.c_str(), &e2, 10);
    return *e1 == 0 && *e2 == 0;
}
std::vector<std::string> Words(const std::string& s) {
    std::vector<std::string> w;
    std::istringstream is(s);
    std::string t;
    while (is >> t) w.push_back(t);
    return w;
}

struct Pending { int line; std::string casing; };
struct Block {
    MbDef def;
    bool open = false, bad = false;
    int headerLine = 0, stretchLine = 0, stars = 0;
    bool hasStretch = false;
    std::vector<int> layerLines, legendLines;
    std::vector<Pending> casingRefs;
};
} // namespace

const MbLegend* MbDef::key(char c) const {
    for (const MbLegend& l : legend) if (l.key == c) return &l;
    return nullptr;
}
const MbDef* MbCatalog::find(const std::string& id) const {
    for (const MbDef& d : defs) if (d.id == id) return &d;
    return nullptr;
}

MbCatalog Multiblock_Parse(const std::string& text, std::vector<std::string>* errors) {
    MbCatalog cat;
    std::vector<std::vector<Pending>> refs; // parallel to cat.defs: casing ids to check once every casing line has been read
    Block cur;
    auto err = [&](int line, const std::string& msg) { if (errors) errors->push_back("line " + std::to_string(line) + ": " + msg); };
    auto bad = [&](int line, const std::string& msg) { cur.bad = true; err(line, msg); };

    auto finish = [&]() {
        MbDef& d = cur.def;
        const int layerCount = (int)d.layers.size();
        // equal row counts and widths (the structure is a box)
        if (layerCount > 0) {
            const size_t rows = d.layers[0].size();
            const size_t width = rows ? d.layers[0][0].size() : 0;
            for (int y = 0; y < layerCount; ++y) {
                if (d.layers[y].size() != rows) bad(cur.layerLines[y], "layers have different row counts");
                for (const std::string& r : d.layers[y])
                    if (r.size() != width) { bad(cur.layerLines[y], "rows have different widths"); break; }
            }
        }
        int controllers = 0;
        for (int y = 0; y < layerCount; ++y)
            for (int z = 0; z < (int)d.layers[y].size(); ++z)
                for (int x = 0; x < (int)d.layers[y][z].size(); ++x) {
                    const MbLegend* l = d.key(d.layers[y][z][x]);
                    if (!l) { bad(cur.layerLines[y], std::string("unknown legend key '") + d.layers[y][z][x] + "'"); continue; }
                    if (l->role == MbRole::Controller) { ++controllers; d.controller = { x, y, z }; }
                }
        if (controllers == 0) bad(cur.headerLine, "no controller cell");
        if (controllers > 1) bad(cur.headerLine, "more than one controller cell");
        if (cur.stars > 1) bad(cur.headerLine, "more than one stretch (*) layer");
        if (cur.hasStretch && cur.stars == 0) bad(cur.stretchLine, "stretch without a * layer");
        if (!cur.hasStretch && cur.stars == 1) bad(cur.headerLine, "a * layer without a stretch line");
        if (cur.hasStretch && (d.stretchMin < 1 || d.stretchMin > d.stretchMax)) bad(cur.stretchLine, "stretch needs 1 <= min <= max");
        if (cat.find(d.id)) bad(cur.headerLine, "duplicate multiblock '" + d.id + "'");
        if (!cur.bad) { cat.defs.push_back(d); refs.push_back(cur.casingRefs); }
        cur = Block();
    };

    std::istringstream in(text);
    std::string raw;
    int lineNo = 0;
    while (std::getline(in, raw)) {
        ++lineNo;
        const size_t hash = raw.find('#');
        const std::string line = Trim(hash == std::string::npos ? raw : raw.substr(0, hash));
        if (line.empty()) continue;
        const std::vector<std::string> w = Words(line);
        const std::string& kw = w[0];

        if (kw == "multiblock") {
            if (cur.open) { err(cur.headerLine, "multiblock '" + cur.def.id + "' has no end"); cur = Block(); }
            const auto c = Split(line.substr(kw.size()), '|');
            cur = Block();
            cur.open = true; cur.headerLine = lineNo;
            if (c.size() != 3 || c[0].empty() || c[1].empty()) { bad(lineNo, "multiblock needs: <id> | <name> | tier <n>"); continue; }
            const auto t = Words(c[2]);
            if (t.size() != 2 || t[0] != "tier") { bad(lineNo, "multiblock needs: <id> | <name> | tier <n>"); continue; }
            cur.def.id = c[0]; cur.def.name = c[1]; cur.def.tier = atoi(t[1].c_str());
            continue;
        }
        if (!cur.open) {
            if (kw == "casing" && w.size() == 4 && w[2] == "tier") {
                cat.casingTier[w[1]] = atoi(w[3].c_str());
                if (std::find(cat.casings.begin(), cat.casings.end(), w[1]) == cat.casings.end()) cat.casings.push_back(w[1]);
            }
            else if (kw == "hatch" && w.size() == 6 && w[2] == "kind" && w[4] == "tier") {
                if (!IsHatchKind(w[3])) err(lineNo, "unknown hatch kind '" + w[3] + "'");
                else {
                    cat.hatches[w[1]] = { w[3], atoi(w[5].c_str()) };
                    if (std::find(cat.hatchIds.begin(), cat.hatchIds.end(), w[1]) == cat.hatchIds.end()) cat.hatchIds.push_back(w[1]);
                }
            }
            else err(lineNo, "unexpected '" + kw + "' outside a multiblock");
            continue;
        }

        MbDef& d = cur.def;
        if (kw == "end") { finish(); continue; }
        if (kw == "legend") {
            if (w.size() < 3 || w[1].size() != 1) { bad(lineNo, "legend needs: <key> <spec>"); continue; }
            MbLegend l;
            l.key = w[1][0];
            if (l.key == '/' || l.key == '*' || l.key == '#') { bad(lineNo, std::string("legend key '") + l.key + "' is reserved"); continue; }
            if (d.key(l.key)) { bad(lineNo, std::string("legend key '") + l.key + "' defined twice"); continue; }
            bool ok = true;
            for (size_t i = 2; i < w.size() && ok; ++i) {
                const std::string& s = w[i];
                const size_t colon = s.find(':');
                const std::string head = s.substr(0, colon), val = colon == std::string::npos ? "" : s.substr(colon + 1);
                if (i == 2 && s == "air") l.role = MbRole::Air;
                else if (i == 2 && s == "any") l.role = MbRole::Any;
                else if (val.empty()) ok = false;
                else if (i == 2 && head == "controller") { l.role = MbRole::Controller; l.id = val; }
                else if (i == 2 && head == "part") { l.role = MbRole::Part; l.id = val; }
                else if (i == 2 && head == "casing") { l.role = MbRole::Casing; l.id = val; cur.casingRefs.push_back({ lineNo, val }); }
                else if (i == 2 && head == "hatch") {
                    l.role = MbRole::Hatch;
                    for (const std::string& k : Split(val, ',')) {
                        if (!IsHatchKind(k)) { bad(lineNo, "unknown hatch kind '" + k + "'"); ok = false; break; }
                        l.hatchKinds.push_back(k);
                    }
                }
                else if (i == 3 && l.role == MbRole::Hatch && head == "casing") { l.fallbackCasing = val; cur.casingRefs.push_back({ lineNo, val }); }
                else ok = false;
            }
            if (cur.bad && !ok) continue;
            if (!ok) { bad(lineNo, "bad legend spec"); continue; }
            d.legend.push_back(l);
            continue;
        }
        if (kw == "hatches") {
            int a = 0, b = 0;
            if (w.size() != 3 || !ParseRange(w[2], a, b)) { bad(lineNo, "hatches needs: <kind> <min>..<max>"); continue; }
            if (!IsHatchKind(w[1])) { bad(lineNo, "unknown hatch kind '" + w[1] + "'"); continue; }
            if (a > b) { bad(lineNo, "hatch min is above max"); continue; }
            d.hatchLimits[w[1]] = { a, b };
            continue;
        }
        if (kw == "stretch") {
            int a = 0, b = 0;
            if (w.size() != 2 || !ParseRange(w[1], a, b)) { bad(lineNo, "stretch needs: <min>..<max>"); continue; }
            cur.hasStretch = true; cur.stretchLine = lineNo;
            d.stretchMin = a; d.stretchMax = b;
            continue;
        }
        if (kw == "layer") {
            std::string rest = Trim(line.substr(kw.size()));
            bool star = false;
            if (rest.size() >= 2 && rest.back() == '*' && (rest[rest.size() - 2] == ' ' || rest[rest.size() - 2] == '\t')) {
                star = true;
                rest = Trim(rest.substr(0, rest.size() - 1));
            }
            std::vector<std::string> rows;
            bool empty = false;
            for (const std::string& r : Split(rest, '/')) {
                std::string row;
                for (char ch : r) if (ch != ' ' && ch != '\t') row += ch;
                if (row.empty()) empty = true;
                rows.push_back(row);
            }
            if (empty) { bad(lineNo, "layer has an empty row"); continue; }
            if (star) { ++cur.stars; d.stretchLayer = (int)d.layers.size(); }
            d.layers.push_back(rows);
            cur.layerLines.push_back(lineNo);
            continue;
        }
        bad(lineNo, "unknown line '" + kw + "' in a multiblock");
    }
    if (cur.open) err(cur.headerLine, "multiblock '" + cur.def.id + "' has no end");

    // A casing id must be declared by a casing line (anywhere in the file).
    for (size_t i = cat.defs.size(); i-- > 0;) {
        bool ok = true;
        for (const Pending& p : refs[i])
            if (!cat.casingTier.count(p.casing)) { err(p.line, "casing '" + p.casing + "' is not declared by a casing line"); ok = false; }
        if (!ok) cat.defs.erase(cat.defs.begin() + i);
    }
    // (#1429 PR 2) the controllers, in definition order (the index is a placed part's material: saved)
    for (const MbDef& d : cat.defs)
        for (const MbLegend& l : d.legend)
            if (l.role == MbRole::Controller && std::find(cat.controllers.begin(), cat.controllers.end(), l.id) == cat.controllers.end()) cat.controllers.push_back(l.id);
    if (cat.controllers.size() > 256) err(0, "more than 256 controllers (a part's material is a byte)");
    if (cat.casings.size() > 256) err(0, "more than 256 casings (a part's material is a byte)");
    if (cat.hatchIds.size() > 256) err(0, "more than 256 hatches (a part's material is a byte)");
    return cat;
}

bool Multiblock_ItemPart(const std::string& item, PartType& type, int& material) {
    const MbCatalog& c = Multiblock_All();
    const struct { const std::vector<std::string>& list; PartType t; } roles[3] = { { c.controllers, PartType::MbController }, { c.casings, PartType::MbCasing }, { c.hatchIds, PartType::MbHatch } };
    for (const auto& r : roles) {
        const auto it = std::find(r.list.begin(), r.list.end(), item);
        if (it != r.list.end() && it - r.list.begin() < 256) { type = r.t; material = (int)(it - r.list.begin()); return true; }
    }
    return false;
}

const char* Multiblock_PartItem(PartType type, int material) {
    const MbCatalog& c = Multiblock_All();
    const std::vector<std::string>* list = type == PartType::MbController ? &c.controllers : type == PartType::MbCasing ? &c.casings : type == PartType::MbHatch ? &c.hatchIds : nullptr;
    if (!list || material < 0 || material >= (int)list->size() || material >= 256) return "";
    return (*list)[material].c_str();
}

const MbCatalog& Multiblock_All() {
    static const MbCatalog all = [] {
        std::vector<std::string> errors;
        MbCatalog c = Multiblock_Parse(ReadFile("resources/data/multiblocks.txt"), &errors);
        for (const std::string& e : errors) LOG_WARN("Multiblock", "%s\n", e.c_str());
        return c;
    }();
    return all;
}

MbInt3 Multiblock_ToWorld(const MbDef& def, MbInt3 defCell, MbInt3 controllerAt, int fx, int fz, int stretchExtra) {
    // Cells above the stretch layer move up with the repeats; so does the controller if it is up there, so the offset
    // from the controller stays right.
    auto lifted = [&](int y) { return y + (def.stretchLayer >= 0 && y > def.stretchLayer ? stretchExtra : 0); };
    const int dx = defCell.x - def.controller.x, dz = defCell.z - def.controller.z;
    const int dy = lifted(defCell.y) - lifted(def.controller.y);
    const int rx = fz, rz = -fx; // right vector: facing (0,1) gives (1,0), the identity
    return { controllerAt.x + dx * rx + dz * fx, controllerAt.y + dy, controllerAt.z + dx * rz + dz * fz };
}

MbResult Multiblock_Validate(const MbDef& def, const MbCatalog& cat, const std::function<std::string(MbInt3)>& lookup,
                             MbInt3 controllerAt, int fx, int fz) {
    MbResult r;
    const MbLegend* ctl = def.key(def.layers.empty() ? 0 : def.layers[def.controller.y][def.controller.z][def.controller.x]);
    const bool controllerOk = ctl && lookup(controllerAt) == ctl->id;

    auto casingFits = [&](const std::string& have, const std::string& want) {
        if (have == want) return true;
        const auto h = cat.casingTier.find(have), w = cat.casingTier.find(want);
        return h != cat.casingTier.end() && w != cat.casingTier.end() && h->second >= w->second;
    };
    // One layer of the definition at a world offset; appends its cells, counts hatches in `counts`.
    auto evalLayer = [&](int defY, int extra, int yAdd, int layerIdx, std::vector<MbCellResult>& out, std::map<std::string, int>& counts) {
        for (int z = 0; z < (int)def.layers[defY].size(); ++z)
            for (int x = 0; x < (int)def.layers[defY][z].size(); ++x) {
                const MbLegend* l = def.key(def.layers[defY][z][x]);
                if (!l || l->role == MbRole::Any) continue;
                MbCellResult c;
                c.world = Multiblock_ToWorld(def, { x, defY, z }, controllerAt, fx, fz, extra);
                c.world.y += yAdd;
                c.layer = layerIdx;
                const std::string p = lookup(c.world);
                c.status = MbCell::Ok;
                switch (l->role) {
                case MbRole::Air:
                    c.want = "air";
                    if (!p.empty()) c.status = MbCell::Blocked;
                    break;
                case MbRole::Controller:
                case MbRole::Part:
                    c.want = l->id;
                    if (p != l->id) c.status = p.empty() ? MbCell::Missing : MbCell::Wrong;
                    break;
                case MbRole::Casing:
                    c.want = l->id;
                    if (p.empty()) c.status = MbCell::Missing;
                    else if (!casingFits(p, l->id)) c.status = MbCell::Wrong;
                    break;
                case MbRole::Hatch: {
                    std::string kinds;
                    for (const std::string& k : l->hatchKinds) kinds += (kinds.empty() ? "" : ",") + k;
                    c.want = l->fallbackCasing.empty() ? "hatch:" + kinds : l->fallbackCasing;
                    const auto h = cat.hatches.find(p);
                    if (h != cat.hatches.end() && std::find(l->hatchKinds.begin(), l->hatchKinds.end(), h->second.first) != l->hatchKinds.end())
                        ++counts[h->second.first];
                    else if (p.empty()) c.status = MbCell::Missing;
                    else if (l->fallbackCasing.empty() || !casingFits(p, l->fallbackCasing)) c.status = MbCell::Wrong;
                    break;
                }
                default: break;
                }
                out.push_back(c);
            }
    };
    auto allOk = [](const std::vector<MbCellResult>& v) {
        for (const MbCellResult& c : v) if (c.status != MbCell::Ok) return false;
        return true;
    };

    const int layerCount = (int)def.layers.size();
    const int sl = def.stretchLayer;
    std::map<std::string, int> counts;
    std::vector<MbCellResult> sink;
    std::map<std::string, int> sinkCounts;

    // The complete repeats of the stretch layer. The controller may sit above the stretch layer and move with it, so the
    // count is settled against the extra it implies (independent of it when the controller is not above).
    int stretch = 0, extra = 0;
    if (sl >= 0) {
        for (int e = 0; e < def.stretchMax; ++e) {
            int n = 0;
            while (n < def.stretchMax) {
                sink.clear(); sinkCounts.clear();
                evalLayer(sl, e, n, sl + n, sink, sinkCounts);
                if (!allOk(sink)) break;
                ++n;
            }
            stretch = n;
            extra = e;
            if (std::max(n, def.stretchMin) - 1 == e) break;
        }
    }
    // The layers above sit over the repeats found, or over the minimum when fewer are built (the ghost's top).
    const int aboveExtra = sl >= 0 ? std::max(stretch, def.stretchMin) - 1 : 0;
    if (sl >= 0 && aboveExtra != extra) { // controller above and the guess disagreed: count again against the final extra
        int n = 0;
        while (n < def.stretchMax) {
            sink.clear(); sinkCounts.clear();
            evalLayer(sl, aboveExtra, n, sl + n, sink, sinkCounts);
            if (!allOk(sink)) break;
            ++n;
        }
        stretch = n;
    }
    r.stretch = stretch;
    const int top = sl >= 0 ? std::max(stretch, def.stretchMin) - 1 : 0;

    for (int y = 0; y < layerCount; ++y) {
        if (sl < 0 || y < sl) evalLayer(y, 0, 0, y, r.cells, counts);
        else if (y == sl) {
            for (int k = 0; k < stretch; ++k) evalLayer(sl, top, k, sl + k, r.cells, counts);
            if (stretch < def.stretchMin) evalLayer(sl, top, stretch, sl + stretch, r.cells, counts); // the next layer to build
        } else evalLayer(y, top, 0, y + top, r.cells, counts);
    }
    r.hatchCount = counts;

    r.nextLayer = -1;
    int missing = 0, wrong = 0, blocked = 0;
    for (const MbCellResult& c : r.cells)
        if (c.status != MbCell::Ok && (r.nextLayer < 0 || c.layer < r.nextLayer)) r.nextLayer = c.layer;
    if (r.nextLayer >= 0)
        for (const MbCellResult& c : r.cells)
            if (c.layer == r.nextLayer) {
                if (c.status == MbCell::Missing) ++missing;
                else if (c.status == MbCell::Wrong) ++wrong;
                else if (c.status == MbCell::Blocked) ++blocked;
            }

    if (!controllerOk) r.problem = "no controller";
    else if (r.nextLayer >= 0) {
        std::string parts;
        auto add = [&](int n, const char* what) { if (n) parts += (parts.empty() ? "" : ", ") + std::to_string(n) + " " + what; };
        add(missing, missing == 1 ? "cell missing" : "cells missing");
        add(wrong, "wrong");
        add(blocked, "blocked");
        r.problem = "layer " + std::to_string(r.nextLayer) + ": " + parts;
    } else if (sl >= 0 && stretch < def.stretchMin) {
        r.problem = "needs " + std::to_string(def.stretchMin - stretch) + " more repeating layers";
    } else {
        for (const auto& lim : def.hatchLimits) {
            const auto it = counts.find(lim.first);
            const int n = it == counts.end() ? 0 : it->second;
            if (n < lim.second.first) { r.problem = "needs at least " + std::to_string(lim.second.first) + " " + lim.first + (lim.second.first == 1 ? " hatch" : " hatches"); break; }
            if (n > lim.second.second) { r.problem = "at most " + std::to_string(lim.second.second) + " " + lim.first + (lim.second.second == 1 ? " hatch" : " hatches") + " allowed"; break; }
        }
    }
    r.formed = r.problem.empty();
    return r;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"

#include <tuple>

namespace {
const char* kDecl =
    "casing bronze_casing tier 3\n"
    "casing steel_casing tier 4\n"
    "hatch fluid_in_hatch kind fluid_in tier 3\n"
    "hatch fluid_out_hatch kind fluid_out tier 3\n";
const char* kTank =
    "multiblock large_tank | Large tank | tier 3\n"
    "  legend C controller:large_tank_controller\n"
    "  legend c casing:bronze_casing\n"
    "  legend h hatch:fluid_in,fluid_out casing:bronze_casing\n"
    "  legend g part:tank_glass\n"
    "  legend . air\n"
    "  legend - any\n"
    "  hatches fluid_in 1..2\n"
    "  hatches fluid_out 1..2\n"
    "  stretch @MIN@..6\n"
    "  layer ccc / chc / cCc\n"
    "  layer ggg / g.g / ggg  *\n"
    "  layer ccc / chc / ccc\n"
    "end\n";
std::string Tank(const char* minRepeats = "1") {
    std::string t = kTank;
    t.replace(t.find("@MIN@"), 5, minRepeats);
    return std::string(kDecl) + t;
}
std::string Bad(const std::string& body) {
    return std::string(kDecl) + "multiblock t | T | tier 1\n legend C controller:ctl\n legend c casing:bronze_casing\n" + body;
}

using World = std::map<std::tuple<int, int, int>, std::string>;
std::function<std::string(MbInt3)> Look(const World& w) {
    return [&w](MbInt3 p) { auto it = w.find({ p.x, p.y, p.z }); return it == w.end() ? std::string() : it->second; };
}
void Put(World& w, MbInt3 p, const std::string& id) { if (id.empty()) w.erase({ p.x, p.y, p.z }); else w[{ p.x, p.y, p.z }] = id; }

// Build a structure with `repeats` stretch layers. hatchIds: the part for each hatch spot in order (default: a fluid_in
// hatch, then a fluid_out hatch, then the fallback casing).
void Stamp(World& w, const MbDef& def, MbInt3 at, int fx, int fz, int repeats = 1, std::vector<std::string> hatchIds = { "fluid_in_hatch", "fluid_out_hatch" }) {
    size_t hatchNo = 0;
    const int top = std::max(repeats, 1) - 1;
    for (int y = 0; y < (int)def.layers.size(); ++y) {
        const int reps = y == def.stretchLayer ? repeats : 1;
        for (int k = 0; k < reps; ++k)
            for (int z = 0; z < (int)def.layers[y].size(); ++z)
                for (int x = 0; x < (int)def.layers[y][z].size(); ++x) {
                    const MbLegend* l = def.key(def.layers[y][z][x]);
                    std::string id;
                    switch (l->role) {
                    case MbRole::Controller: case MbRole::Casing: case MbRole::Part: id = l->id; break;
                    case MbRole::Hatch: id = hatchNo < hatchIds.size() ? hatchIds[hatchNo] : l->fallbackCasing; if (k == 0) ++hatchNo; break;
                    default: break;
                    }
                    if (id.empty()) continue;
                    MbInt3 p = Multiblock_ToWorld(def, { x, y, z }, at, fx, fz, top);
                    p.y += k;
                    Put(w, p, id);
                }
    }
}
MbInt3 Cell(const MbDef& d, MbInt3 c, MbInt3 at, int fx, int fz, int extra = 0) { return Multiblock_ToWorld(d, c, at, fx, fz, extra); }
const int kFacings[4][2] = { { 0, 1 }, { 1, 0 }, { 0, -1 }, { -1, 0 } };
} // namespace

TEST(Multiblock, ParsesTheLargeTank) {
    std::vector<std::string> errors;
    const MbCatalog cat = Multiblock_Parse(Tank(), &errors);
    EXPECT_TRUE(errors.empty());
    ASSERT_EQ(cat.defs.size(), 1u);
    const MbDef* d = cat.find("large_tank");
    ASSERT_NE(d, nullptr);
    EXPECT_EQ(d->name, "Large tank");
    EXPECT_EQ(d->tier, 3);
    EXPECT_EQ(d->legend.size(), 6u);
    ASSERT_NE(d->key('h'), nullptr);
    EXPECT_EQ(d->key('h')->role, MbRole::Hatch);
    EXPECT_EQ(d->key('h')->hatchKinds.size(), 2u);
    EXPECT_EQ(d->key('h')->fallbackCasing, "bronze_casing");
    EXPECT_EQ(d->key('g')->id, "tank_glass");
    EXPECT_EQ(d->key('z'), nullptr);
    ASSERT_EQ(d->layers.size(), 3u);
    EXPECT_EQ(d->layers[0][2], "cCc");
    EXPECT_EQ(d->layers[1][1], "g.g");
    EXPECT_TRUE(d->controller == (MbInt3{ 1, 0, 2 }));
    EXPECT_EQ(d->stretchLayer, 1);
    EXPECT_EQ(d->stretchMin, 1);
    EXPECT_EQ(d->stretchMax, 6);
    EXPECT_EQ(d->hatchLimits.at("fluid_in"), (std::pair<int, int>{ 1, 2 }));
    EXPECT_EQ(cat.casingTier.at("steel_casing"), 4);
    EXPECT_EQ(cat.hatches.at("fluid_out_hatch").first, "fluid_out");
}

TEST(Multiblock, ParseErrorsDropTheBlock) {
    const char* bodies[] = {
        " layer cCx / ccc\nend\n",                              // unknown legend key
        " layer cCc / cc\nend\n",                               // rows of unequal width
        " layer cCc / ccc\n layer ccc\nend\n",                  // layers with unequal row counts
        " layer ccc\nend\n",                                    // no controller
        " layer cCC\nend\n",                                    // two controller cells
        " stretch 1..2\n layer cCc *\n layer ccc *\nend\n",     // two * layers
        " stretch 1..2\n layer cCc\nend\n",                     // stretch without a * layer
        " layer cCc *\nend\n",                                  // a * layer without stretch
        " stretch 0..2\n layer cCc *\nend\n",                   // min < 1
        " stretch 3..2\n layer cCc *\nend\n",                   // min > max
        " legend s casing:nope\n layer cCs\nend\n",             // undeclared casing
        " legend h hatch:warp\n layer cCh\nend\n",              // not a hatch kind
        " hatches fluid_in 3..1\n layer cCc\nend\n",            // hatch min > max
        " layer cCc\n",                                         // missing end
    };
    for (const char* b : bodies) {
        std::vector<std::string> errors;
        const MbCatalog cat = Multiblock_Parse(Bad(b), &errors);
        EXPECT_TRUE(cat.defs.empty()) << b;
        ASSERT_FALSE(errors.empty()) << b;
        EXPECT_EQ(errors[0].rfind("line ", 0), 0u) << errors[0];
    }
    // A bad block does not take a good one with it.
    std::vector<std::string> errors;
    const MbCatalog cat = Multiblock_Parse(Bad(" layer ccc\nend\n") + "multiblock ok | Ok | tier 1\n legend C controller:ctl\n layer C\nend\n", &errors);
    EXPECT_EQ(cat.defs.size(), 1u);
    EXPECT_NE(cat.find("ok"), nullptr);
}

TEST(Multiblock, RotationForTheFourFacings) {
    const MbCatalog cat = Multiblock_Parse(Tank(), nullptr);
    const MbDef& d = *cat.find("large_tank");
    const MbInt3 at{ 10, 5, 10 };
    // The def cell (0,0,0) is 1 west and 2 north of the controller (1,0,2); written out by hand per facing.
    EXPECT_TRUE(Multiblock_ToWorld(d, { 0, 0, 0 }, at, 0, 1) == (MbInt3{ 9, 5, 8 }));
    EXPECT_TRUE(Multiblock_ToWorld(d, { 0, 0, 0 }, at, 1, 0) == (MbInt3{ 8, 5, 11 }));
    EXPECT_TRUE(Multiblock_ToWorld(d, { 0, 0, 0 }, at, 0, -1) == (MbInt3{ 11, 5, 12 }));
    EXPECT_TRUE(Multiblock_ToWorld(d, { 0, 0, 0 }, at, -1, 0) == (MbInt3{ 12, 5, 9 }));
    // Above the stretch layer the cell moves up by the extra repeats.
    EXPECT_TRUE(Multiblock_ToWorld(d, { 0, 2, 0 }, at, 0, 1, 3) == (MbInt3{ 9, 10, 8 }));
    for (const auto& f : kFacings) {
        World w;
        Stamp(w, d, at, f[0], f[1]);
        const MbResult r = Multiblock_Validate(d, cat, Look(w), at, f[0], f[1]);
        EXPECT_TRUE(r.formed) << r.problem << " facing " << f[0] << "," << f[1];
    }
    // Facing south, a structure built for another facing is not formed.
    World w;
    Stamp(w, d, at, 1, 0);
    EXPECT_FALSE(Multiblock_Validate(d, cat, Look(w), at, 0, 1).formed);
}

TEST(Multiblock, MissingWrongBlocked) {
    const MbCatalog cat = Multiblock_Parse(Tank(), nullptr);
    const MbDef& d = *cat.find("large_tank");
    const MbInt3 at{ 10, 5, 10 };
    World w;
    Stamp(w, d, at, 0, 1);
    ASSERT_TRUE(Multiblock_Validate(d, cat, Look(w), at, 0, 1).formed);

    World miss = w;
    Put(miss, Cell(d, { 0, 0, 0 }, at, 0, 1), "");
    MbResult r = Multiblock_Validate(d, cat, Look(miss), at, 0, 1);
    EXPECT_FALSE(r.formed);
    EXPECT_EQ(r.nextLayer, 0);
    EXPECT_NE(r.problem.find("layer 0"), std::string::npos) << r.problem;
    EXPECT_NE(r.problem.find("1 cell missing"), std::string::npos) << r.problem;
    int missing = 0;
    for (const MbCellResult& c : r.cells) if (c.status == MbCell::Missing) { ++missing; EXPECT_EQ(c.want, "bronze_casing"); }
    EXPECT_EQ(missing, 1);

    World wrong = w;
    Put(wrong, Cell(d, { 0, 1, 0 }, at, 0, 1), "bronze_casing"); // casing where the glass goes
    r = Multiblock_Validate(d, cat, Look(wrong), at, 0, 1);
    EXPECT_FALSE(r.formed);
    int nWrong = 0;
    for (const MbCellResult& c : r.cells) if (c.status == MbCell::Wrong) ++nWrong;
    EXPECT_EQ(nWrong, 1);
    EXPECT_NE(r.problem.find("layer 1"), std::string::npos) << r.problem;

    World blocked = w;
    Put(blocked, Cell(d, { 1, 1, 1 }, at, 0, 1), "tank_glass"); // a part in the air cell
    r = Multiblock_Validate(d, cat, Look(blocked), at, 0, 1);
    EXPECT_FALSE(r.formed);
    bool sawBlocked = false;
    for (const MbCellResult& c : r.cells) if (c.status == MbCell::Blocked) { sawBlocked = true; EXPECT_EQ(c.want, "air"); }
    EXPECT_TRUE(sawBlocked);

    World noCtl = w;
    Put(noCtl, at, "");
    EXPECT_EQ(Multiblock_Validate(d, cat, Look(noCtl), at, 0, 1).problem, "no controller");
}

TEST(Multiblock, CasingTierAcceptsBetterNotWorse) {
    const MbCatalog cat = Multiblock_Parse(Tank() + "multiblock box | Box | tier 4\n legend C controller:large_tank_controller\n legend s casing:steel_casing\n layer sCs\nend\n", nullptr);
    const MbDef& tank = *cat.find("large_tank");
    const MbInt3 at{ 10, 5, 10 };
    World w;
    Stamp(w, tank, at, 0, 1);
    Put(w, Cell(tank, { 0, 0, 0 }, at, 0, 1), "steel_casing");
    EXPECT_TRUE(Multiblock_Validate(tank, cat, Look(w), at, 0, 1).formed) << "steel fits a bronze cell";

    const MbDef& box = *cat.find("box");
    World b;
    Stamp(b, box, at, 0, 1);
    EXPECT_TRUE(Multiblock_Validate(box, cat, Look(b), at, 0, 1).formed);
    Put(b, { 9, 5, 10 }, "bronze_casing");
    const MbResult r = Multiblock_Validate(box, cat, Look(b), at, 0, 1);
    EXPECT_FALSE(r.formed) << "bronze does not fit a steel cell";
    EXPECT_NE(r.problem.find("1 wrong"), std::string::npos) << r.problem;
}

TEST(Multiblock, Hatches) {
    const MbCatalog cat = Multiblock_Parse(Tank() +
        "multiblock hh | HH | tier 3\n legend C controller:large_tank_controller\n legend h hatch:fluid_in casing:bronze_casing\n"
        " hatches fluid_in 1..2\n layer hhCh\nend\n", nullptr);
    const MbDef& tank = *cat.find("large_tank");
    const MbInt3 at{ 10, 5, 10 };
    World w;
    Stamp(w, tank, at, 0, 1);
    MbResult r = Multiblock_Validate(tank, cat, Look(w), at, 0, 1);
    EXPECT_TRUE(r.formed) << r.problem;
    EXPECT_EQ(r.hatchCount["fluid_in"], 1);
    EXPECT_EQ(r.hatchCount["fluid_out"], 1);

    // Both hatch spots with the fallback casing: no hatch at all.
    World none;
    Stamp(none, tank, at, 0, 1, 1, {});
    r = Multiblock_Validate(tank, cat, Look(none), at, 0, 1);
    EXPECT_FALSE(r.formed);
    EXPECT_NE(r.problem.find("fluid_in"), std::string::npos) << r.problem;
    EXPECT_NE(r.problem.find("at least"), std::string::npos) << r.problem;

    // The spec's three: the "hh" definition has three hatch spots.
    const MbDef& hh = *cat.find("hh");
    World h3;
    Stamp(h3, hh, at, 0, 1, 1, { "fluid_in_hatch", "fluid_in_hatch", "fluid_in_hatch" });
    r = Multiblock_Validate(hh, cat, Look(h3), at, 0, 1);
    EXPECT_FALSE(r.formed);
    EXPECT_EQ(r.hatchCount["fluid_in"], 3);
    EXPECT_NE(r.problem.find("at most"), std::string::npos) << r.problem;
    World h2;
    Stamp(h2, hh, at, 0, 1, 1, { "fluid_in_hatch", "fluid_in_hatch" });
    EXPECT_TRUE(Multiblock_Validate(hh, cat, Look(h2), at, 0, 1).formed);
    // A hatch of the wrong kind is not accepted.
    World wk;
    Stamp(wk, hh, at, 0, 1, 1, { "fluid_out_hatch", "fluid_in_hatch" });
    EXPECT_FALSE(Multiblock_Validate(hh, cat, Look(wk), at, 0, 1).formed);
}

TEST(Multiblock, StretchLayers) {
    const MbCatalog cat = Multiblock_Parse(Tank(), nullptr);
    const MbDef& d = *cat.find("large_tank");
    const MbInt3 at{ 10, 5, 10 };
    World w;
    Stamp(w, d, at, 0, 1, 3);
    MbResult r = Multiblock_Validate(d, cat, Look(w), at, 0, 1);
    EXPECT_EQ(r.stretch, 3);
    EXPECT_TRUE(r.formed) << r.problem;
    // The top layer is checked two above where it sits with one repeat: its stretched index is 4.
    int topCells = 0;
    for (const MbCellResult& c : r.cells) if (c.layer == 4) { ++topCells; EXPECT_EQ(c.world.y, 5 + 4); }
    EXPECT_EQ(topCells, 9);
    World broken = w;
    Put(broken, Cell(d, { 0, 2, 0 }, at, 0, 1, 2), "");
    r = Multiblock_Validate(d, cat, Look(broken), at, 0, 1);
    EXPECT_FALSE(r.formed);
    EXPECT_EQ(r.nextLayer, 4);

    // Two repeats built with min 1, and a partly built third where the top layer should be: the top reads wrong; the
    // third repeat is not reported as a layer to build.
    World partial;
    Stamp(partial, d, at, 0, 1, 2);
    for (int z = 0; z < 3; ++z) for (int x = 0; x < 3; ++x) Put(partial, Cell(d, { x, 2, z }, at, 0, 1, 1), ""); // lift the top off...
    Put(partial, { 9, 8, 8 }, "tank_glass"); // ...and start the third glass layer
    r = Multiblock_Validate(d, cat, Look(partial), at, 0, 1);
    EXPECT_EQ(r.stretch, 2);
    EXPECT_FALSE(r.formed);
    for (const MbCellResult& c : r.cells) EXPECT_FALSE(c.layer >= 3 && c.want == "tank_glass") << "layer " << c.layer;
    for (const MbCellResult& c : r.cells) { if (c.world == (MbInt3{ 9, 8, 8 })) { EXPECT_EQ(c.status, MbCell::Wrong); } }

    // Minimum 2 with one repeat: not formed, and the ghost shows the second repeat's cells as Missing.
    const MbCatalog cat2 = Multiblock_Parse(Tank("2"), nullptr);
    const MbDef& d2 = *cat2.find("large_tank");
    World one;
    Stamp(one, d2, at, 0, 1, 1);
    for (int z = 0; z < 3; ++z) for (int x = 0; x < 3; ++x) Put(one, Cell(d2, { x, 2, z }, at, 0, 1, 0), ""); // top not built yet (it would sit where the 2nd repeat goes)
    r = Multiblock_Validate(d2, cat2, Look(one), at, 0, 1);
    EXPECT_EQ(r.stretch, 1);
    EXPECT_FALSE(r.formed);
    EXPECT_EQ(r.nextLayer, 2);
    int missing = 0;
    for (const MbCellResult& c : r.cells) if (c.layer == 2 && c.status == MbCell::Missing) { ++missing; EXPECT_EQ(c.want, "tank_glass"); }
    EXPECT_EQ(missing, 8);
    World two;
    Stamp(two, d2, at, 0, 1, 2);
    EXPECT_TRUE(Multiblock_Validate(d2, cat2, Look(two), at, 0, 1).formed);
}

TEST(Multiblock, TwoStructuresMayShareAWall) {
    const MbCatalog cat = Multiblock_Parse(Tank(), nullptr);
    const MbDef& d = *cat.find("large_tank");
    World w;
    Stamp(w, d, { 10, 5, 10 }, 0, 1);
    Stamp(w, d, { 12, 5, 10 }, 0, 1); // the second tank's west wall is the first one's east wall
    EXPECT_TRUE(Multiblock_Validate(d, cat, Look(w), { 10, 5, 10 }, 0, 1).formed);
    EXPECT_TRUE(Multiblock_Validate(d, cat, Look(w), { 12, 5, 10 }, 0, 1).formed);
}

TEST(Multiblock, TheGameFileParses) {
    std::vector<std::string> errors;
    std::ifstream f("resources/data/multiblocks.txt");
    ASSERT_TRUE(f.good()) << "run from the repo root";
    std::stringstream ss;
    ss << f.rdbuf();
    const MbCatalog cat = Multiblock_Parse(ss.str(), &errors);
    for (const std::string& e : errors) ADD_FAILURE() << e;
    EXPECT_GE(cat.defs.size(), 1u);
}

// (#1429 PR 2) every id of the real catalog's three lists is a part item, and back
TEST(Multiblock, PartItemsRoundTrip) {
    const MbCatalog& cat = Multiblock_All();
    ASSERT_FALSE(cat.controllers.empty());
    ASSERT_FALSE(cat.casings.empty());
    ASSERT_FALSE(cat.hatchIds.empty());
    const struct { const std::vector<std::string>& list; PartType t; } roles[3] = { { cat.controllers, PartType::MbController }, { cat.casings, PartType::MbCasing }, { cat.hatchIds, PartType::MbHatch } };
    for (const auto& r : roles)
        for (int i = 0; i < (int)r.list.size(); ++i) {
            PartType t = PartType::Count; int m = -1;
            ASSERT_TRUE(Multiblock_ItemPart(r.list[i], t, m)) << r.list[i];
            EXPECT_EQ((int)t, (int)r.t) << r.list[i];
            EXPECT_EQ(m, i) << r.list[i];
            EXPECT_STREQ(Multiblock_PartItem(r.t, i), r.list[i].c_str());
            EXPECT_STREQ(PartItem(Part{ r.t, { 0, 0, 0 }, 0, (uint8_t)i }), r.list[i].c_str()) << "a placed part gives its item back";
        }
    PartType t = PartType::Count; int m = -1;
    EXPECT_FALSE(Multiblock_ItemPart("iron_plate", t, m));
    EXPECT_FALSE(Multiblock_ItemPart("", t, m));
    EXPECT_STREQ(Multiblock_PartItem(PartType::MbCasing, (int)cat.casings.size()), "");
    EXPECT_STREQ(Multiblock_PartItem(PartType::Cube, 0), "");
}

// (#1429 PR 2) a large tank built of real Part structs validates through PartItem
TEST(Multiblock, ValidatesRealParts) {
    const MbCatalog& cat = Multiblock_All();
    const MbDef* d = cat.find("large_tank");
    ASSERT_NE(d, nullptr);
    const int glass = BuildMaterialForItem("glass");
    ASSERT_GE(glass, 0);
    const MbInt3 at{ 10, 5, 10 };
    World ids;
    Stamp(ids, *d, at, 0, 1, 2); // two stretch layers
    std::map<std::tuple<int, int, int>, Part> parts;
    for (const auto& [k, id] : ids) {
        PartType t = PartType::Cube; int m = 0;
        if (!Multiblock_ItemPart(id, t, m)) { ASSERT_EQ(id, "glass"); m = glass; }
        parts[k] = Part{ t, { std::get<0>(k), std::get<1>(k), std::get<2>(k) }, 0, (uint8_t)m };
    }
    const auto look = [&parts](MbInt3 p) { const auto it = parts.find({ p.x, p.y, p.z }); return it == parts.end() ? std::string() : std::string(PartItem(it->second)); };
    MbResult r = Multiblock_Validate(*d, cat, look, at, 0, 1);
    EXPECT_TRUE(r.formed) << r.problem;
    EXPECT_EQ(r.hatchCount["fluid_in"], 1);
    EXPECT_EQ(r.hatchCount["fluid_out"], 1);
    // a steel casing in a bronze casing's place is fine (a better casing); a hatch part where glass goes is not
    const MbInt3 casingCell = Multiblock_ToWorld(*d, { 0, 0, 0 }, at, 0, 1, 1);
    PartType ct; int cm = 0;
    ASSERT_TRUE(Multiblock_ItemPart("steel_casing", ct, cm));
    parts[{ casingCell.x, casingCell.y, casingCell.z }] = Part{ ct, { casingCell.x, casingCell.y, casingCell.z }, 0, (uint8_t)cm };
    EXPECT_TRUE(Multiblock_Validate(*d, cat, look, at, 0, 1).formed);
    const MbInt3 glassCell = Multiblock_ToWorld(*d, { 0, 1, 0 }, at, 0, 1, 1);
    ASSERT_TRUE(Multiblock_ItemPart("fluid_in_hatch", ct, cm));
    parts[{ glassCell.x, glassCell.y, glassCell.z }] = Part{ ct, { glassCell.x, glassCell.y, glassCell.z }, 0, (uint8_t)cm };
    EXPECT_FALSE(Multiblock_Validate(*d, cat, look, at, 0, 1).formed);
}
#endif
