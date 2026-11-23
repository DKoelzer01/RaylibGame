#include "conduit.h"
#include "utils/logger.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <set>
#include <sstream>

// A run is a polyline of axis-aligned segments; a slot is owned by at most one run. The placement check only reads the
// occupancy and the runs it is given, so the game and the tests pass their own maps.

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
std::vector<std::string> Words(const std::string& s) {
    std::vector<std::string> w;
    std::istringstream is(s);
    std::string t;
    while (is >> t) w.push_back(t);
    return w;
}
std::string ReadFile(const char* path) {
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}
bool ParseInt(const std::string& s, int& v) {
    if (s.empty()) return false;
    char* e = nullptr;
    const long r = strtol(s.c_str(), &e, 10);
    if (*e != 0) return false;
    v = (int)r;
    return true;
}

int Get(const CInt3& p, int a) { return a == 0 ? p.x : a == 1 ? p.y : p.z; }
void Set(CInt3& p, int a, int v) { (a == 0 ? p.x : a == 1 ? p.y : p.z) = v; }
CInt3 Add(CInt3 p, int a, int d) { Set(p, a, Get(p, a) + d); return p; }

// The axis a segment runs along, -1 when it is zero-length or diagonal.
int AxisOf(const CInt3& a, const CInt3& b) {
    int axis = -1, n = 0;
    for (int i = 0; i < 3; ++i) if (Get(a, i) != Get(b, i)) { axis = i; ++n; }
    return n == 1 ? axis : -1;
}

// The slots c + (da, db) on the two cross axes of `axis` around one centreline slot.
void AddBlock(CInt3 c, int axis, int size, std::vector<CInt3>& out) {
    const int ca = (axis + 1) % 3, cb = (axis + 2) % 3;
    const int lo = -(size / 2), hi = lo + size - 1;
    for (int da = lo; da <= hi; ++da)
        for (int db = lo; db <= hi; ++db) out.push_back(Add(Add(c, ca, da), cb, db));
}
void SegSlots(const CInt3& a, const CInt3& b, int axis, int size, std::vector<CInt3>& out) {
    const int dir = Get(b, axis) > Get(a, axis) ? 1 : -1;
    const int len = std::abs(Get(b, axis) - Get(a, axis));
    for (int t = 0; t <= len; ++t) AddBlock(Add(a, axis, dir * t), axis, size, out);
}

int FloorDiv8(int v) { return (v >= 0 ? v : v - 7) / 8; }

const char* const kPortKinds[] = { "power_in", "power_out", "fluid_in", "fluid_out", "item_in", "item_out", "rotation" };
bool PortKindFits(ConduitKind run, ConduitPortKind port) {
    if (run == ConduitKind::Cable) return port == ConduitPortKind::PowerIn || port == ConduitPortKind::PowerOut;
    return port == ConduitPortKind::FluidIn || port == ConduitPortKind::FluidOut;
}
} // namespace

CInt3 Conduit_CellOf(CInt3 slot) { return { FloorDiv8(slot.x), FloorDiv8(slot.y), FloorDiv8(slot.z) }; }

const ConduitMaterial* ConduitCatalog::find(const std::string& id) const {
    for (const ConduitMaterial& m : mats) if (m.id == id) return &m;
    return nullptr;
}

ConduitCatalog Conduit_ParseMaterials(const std::string& text, std::vector<std::string>* errors) {
    ConduitCatalog cat;
    std::istringstream in(text);
    std::string raw;
    int lineNo = 0;
    while (std::getline(in, raw)) {
        ++lineNo;
        const std::string line = Trim(raw);
        if (line.empty() || line[0] == '#') continue;
        auto err = [&](const std::string& msg) { if (errors) errors->push_back("line " + std::to_string(lineNo) + ": " + msg); };
        const std::vector<std::string> parts = Split(line, '|');
        const std::vector<std::string> head = Words(parts[0]);
        ConduitMaterial m;
        if (head.size() != 2 || (head[0] != "cable" && head[0] != "pipe")) { err("expected 'cable <id>' or 'pipe <id>'"); continue; }
        m.kind = head[0] == "cable" ? ConduitKind::Cable : ConduitKind::Pipe;
        m.id = head[1];
        bool ok = true;
        for (size_t i = 1; i < parts.size() && ok; ++i) {
            const std::vector<std::string> kv = Words(parts[i]);
            int n = 0;
            if (kv.size() != 2) { err("expected '<key> <value>' in '" + parts[i] + "'"); ok = false; break; }
            if (kv[0] == "item") { m.item = kv[1]; continue; }
            if (!ParseInt(kv[1], n) || n < 0) { err("bad number in '" + parts[i] + "'"); ok = false; break; }
            if (kv[0] == "units") m.unitsPerItem = n;
            else if (kv[0] == "volts" && m.kind == ConduitKind::Cable) m.volts = n;
            else if (kv[0] == "flow" && m.kind == ConduitKind::Pipe) m.flow = n;
            else if (kv[0] == "maxtemp" && m.kind == ConduitKind::Pipe) m.maxTempC = n;
            else { err("unknown field '" + kv[0] + "' for this kind"); ok = false; }
        }
        if (!ok) continue;
        if (m.item.empty()) { err("missing 'item'"); continue; }
        if (m.unitsPerItem < 1) { err("'units' must be at least 1"); continue; }
        if (cat.find(m.id)) { err("duplicate material '" + m.id + "'"); continue; }
        cat.mats.push_back(m);
    }
    return cat;
}

const ConduitCatalog& Conduit_Materials() {
    static const ConduitCatalog all = [] {
        std::vector<std::string> errors;
        ConduitCatalog c = Conduit_ParseMaterials(ReadFile("resources/data/conduits.txt"), &errors);
        for (const std::string& e : errors) LOG_WARN("Conduit", "%s\n", e.c_str());
        return c;
    }();
    return all;
}

bool Conduit_ValidRun(const ConduitRun& run, std::string* why) {
    auto fail = [&](const std::string& m) { if (why) *why = m; return false; };
    if (run.points.size() < 2) return fail("a run needs at least 2 points");
    if (run.size < 1 || run.size > 4) return fail("size must be 1..4");
    int prevAxis = -1, prevDir = 0;
    for (size_t i = 1; i < run.points.size(); ++i) {
        const int axis = AxisOf(run.points[i - 1], run.points[i]);
        if (axis < 0) return fail("segment " + std::to_string(i) + " is zero-length or diagonal");
        const int dir = Get(run.points[i], axis) > Get(run.points[i - 1], axis) ? 1 : -1;
        if (axis == prevAxis && dir == -prevDir) return fail("segment " + std::to_string(i) + " reverses onto the previous one");
        prevAxis = axis;
        prevDir = dir;
    }
    return true;
}

std::vector<CInt3> Conduit_Slots(const ConduitRun& run) {
    std::vector<CInt3> out;
    if (!Conduit_ValidRun(run, nullptr)) return out;
    for (size_t i = 1; i < run.points.size(); ++i)
        SegSlots(run.points[i - 1], run.points[i], AxisOf(run.points[i - 1], run.points[i]), run.size, out);
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

int64_t Conduit_Length(const ConduitRun& run) {
    int64_t len = 0;
    for (size_t i = 1; i < run.points.size(); ++i) {
        const CInt3& a = run.points[i - 1];
        const CInt3& b = run.points[i];
        len += (int64_t)std::abs(a.x - b.x) + std::abs(a.y - b.y) + std::abs(a.z - b.z);
    }
    return len;
}

int Conduit_Cost(const ConduitRun& run, const ConduitMaterial& mat) {
    const int64_t units = mat.unitsPerItem > 0 ? mat.unitsPerItem : 1;
    const int64_t need = Conduit_Length(run) * run.size * run.size;
    return (int)((need + units - 1) / units);
}

// ---- occupancy ----

namespace {
int BitIndex(CInt3 slot) {
    const CInt3 c = Conduit_CellOf(slot);
    return (slot.x - c.x * 8) + 8 * ((slot.y - c.y * 8) + 8 * (slot.z - c.z * 8));
}
} // namespace

void ConduitOccupancy::add(const ConduitRun& run) {
    for (const CInt3& s : Conduit_Slots(run)) {
        std::array<uint64_t, 8>& m = mask[Conduit_CellOf(s)]; // value-initialised to zero on first use
        const int b = BitIndex(s);
        m[b >> 6] |= (uint64_t)1 << (b & 63);
        owner[s] = run.id;
    }
}

void ConduitOccupancy::remove(const ConduitRun& run) {
    for (const CInt3& s : Conduit_Slots(run)) {
        auto o = owner.find(s);
        if (o == owner.end() || o->second != run.id) continue;
        owner.erase(o);
        auto it = mask.find(Conduit_CellOf(s));
        if (it == mask.end()) continue;
        const int b = BitIndex(s);
        it->second[b >> 6] &= ~((uint64_t)1 << (b & 63));
        bool any = false;
        for (uint64_t w : it->second) if (w) any = true;
        if (!any) mask.erase(it);
    }
}

int ConduitOccupancy::ownerAt(CInt3 slot) const {
    auto it = owner.find(slot);
    return it == owner.end() ? 0 : it->second;
}
bool ConduitOccupancy::cellUsed(CInt3 cell) const { return mask.count(cell) != 0; }
size_t ConduitOccupancy::cellCount() const { return mask.size(); }

// ---- placement ----

ConduitPlaceResult Conduit_CheckPlace(const ConduitRun& run, const ConduitOccupancy& occ, const std::map<int, ConduitRun>& runs,
                                      const std::function<bool(CInt3 cell)>& solidCell, const std::vector<ConduitPortAt>& ports) {
    ConduitPlaceResult res;
    if (!Conduit_ValidRun(run, &res.why)) return res;

    auto compatible = [&](int id) {
        auto it = runs.find(id);
        return it != runs.end() && it->second.kind == run.kind && it->second.material == run.material;
    };

    // The two ends: the point, the axis of its segment, its cap (the size x size block at the end).
    struct End { CInt3 pt; std::set<CInt3> cap; int owner = 0; bool portOk = false; };
    End ends[2];
    for (int e = 0; e < 2; ++e) {
        const size_t n = run.points.size();
        const CInt3 pt = e == 0 ? run.points[0] : run.points[n - 1];
        const CInt3 other = e == 0 ? run.points[1] : run.points[n - 2];
        std::vector<CInt3> cap;
        AddBlock(pt, AxisOf(pt, other), run.size, cap);
        ends[e].pt = pt;
        ends[e].cap.insert(cap.begin(), cap.end());
        ends[e].owner = occ.ownerAt(pt);
        if (ends[e].owner == run.id) ends[e].owner = 0;
        if (ends[e].owner != 0 && compatible(ends[e].owner)) res.joinsRuns.push_back(ends[e].owner);
        for (size_t i = 0; i < ports.size(); ++i) {
            const ConduitPortAt& p = ports[i];
            if (!(p.slot == pt) || !PortKindFits(run.kind, p.kind) || run.size > p.maxSize) continue;
            ends[e].portOk = true;
            res.portIdx.push_back((int)i);
            break;
        }
    }

    std::set<CInt3> blocked;
    for (const CInt3& s : Conduit_Slots(run)) {
        const int o = occ.ownerAt(s);
        if (o != 0 && o != run.id) {
            bool meets = false;
            for (const End& e : ends)
                if (e.owner == o && compatible(o) && e.cap.count(s)) meets = true;
            if (!meets) { blocked.insert(s); continue; }
        }
        if (solidCell && solidCell(Conduit_CellOf(s))) {
            bool exempt = false;
            for (const End& e : ends) if (e.portOk && e.cap.count(s)) exempt = true;
            if (!exempt) blocked.insert(s);
        }
    }

    std::sort(res.joinsRuns.begin(), res.joinsRuns.end());
    res.joinsRuns.erase(std::unique(res.joinsRuns.begin(), res.joinsRuns.end()), res.joinsRuns.end());
    std::sort(res.portIdx.begin(), res.portIdx.end());
    res.portIdx.erase(std::unique(res.portIdx.begin(), res.portIdx.end()), res.portIdx.end());
    res.blocked.assign(blocked.begin(), blocked.end());
    res.ok = res.blocked.empty();
    if (!res.ok) res.why = std::to_string(res.blocked.size()) + " slot(s) blocked";
    return res;
}

// ---- ports ----

bool Conduit_ParsePorts(const std::string& spec, std::vector<ConduitPortDef>& out, std::string* why) {
    auto fail = [&](const std::string& m) { if (why) *why = m; return false; };
    std::vector<ConduitPortDef> tmp;
    for (const std::string& tok : Words(spec)) {
        const std::vector<std::string> parts = Split(tok, ':');
        if (parts.size() < 2 || parts.size() > 3) return fail("bad port '" + tok + "'");
        ConduitPortDef d;
        const size_t f = std::string("NESWUD").find(parts[0].size() == 1 ? parts[0][0] : '?');
        if (f == std::string::npos) return fail("bad face in '" + tok + "'");
        d.face = (int)f;
        std::string kind = parts[1];
        const size_t at = kind.find('@');
        if (at != std::string::npos) {
            const std::vector<std::string> uv = Split(kind.substr(at + 1), ',');
            if (uv.size() != 2 || !ParseInt(uv[0], d.u) || !ParseInt(uv[1], d.v) || d.u < 0 || d.u > 7 || d.v < 0 || d.v > 7)
                return fail("bad u,v in '" + tok + "'");
            kind = kind.substr(0, at);
        }
        int k = -1;
        for (int i = 0; i < 7; ++i) if (kind == kPortKinds[i]) k = i;
        if (k < 0) return fail("unknown port kind in '" + tok + "'");
        d.kind = (ConduitPortKind)k;
        if (parts.size() == 3) {
            if (parts[2].compare(0, 3, "max") != 0 || !ParseInt(parts[2].substr(3), d.maxSize) || d.maxSize < 1 || d.maxSize > 4)
                return fail("bad max in '" + tok + "'");
        }
        tmp.push_back(d);
    }
    out = tmp;
    return true;
}

std::vector<ConduitPortDef> Conduit_DefaultPorts(bool eu, bool generator, bool fluid) {
    std::vector<ConduitPortDef> out;
    auto add = [&](int face, ConduitPortKind k) { ConduitPortDef d; d.face = face; d.kind = k; out.push_back(d); };
    if (eu) for (int f = 0; f < 6; ++f) add(f, ConduitPortKind::PowerIn);
    if (generator) for (int f = 0; f < 6; ++f) add(f, ConduitPortKind::PowerOut);
    if (fluid) {
        add(4, ConduitPortKind::FluidIn); add(0, ConduitPortKind::FluidIn); add(3, ConduitPortKind::FluidIn); // U, N, W
        add(5, ConduitPortKind::FluidOut); add(2, ConduitPortKind::FluidOut); add(1, ConduitPortKind::FluidOut); // D, S, E
    }
    return out;
}

CInt3 Conduit_PortSlot(CInt3 partCell, int partFacing, const ConduitPortDef& port) {
    // Local slot (0..7 per axis) in the part's own frame, written out per face. u runs left to right and v bottom to top
    // seen from outside the face: N looks toward +z (its right is -x), S toward -z (right +x), E toward -x (right -z),
    // W toward +x (right +z). U and D: u = x, v = z.
    const int u = port.u, v = port.v;
    int x = 0, y = 0, z = 0;
    switch (port.face) {
        case 0: x = 7 - u; y = v; z = 0; break; // N (-z)
        case 1: x = 7; y = v; z = 7 - u; break; // E (+x)
        case 2: x = u; y = v; z = 7; break;     // S (+z)
        case 3: x = 0; y = v; z = u; break;     // W (-x)
        case 4: x = u; y = 7; z = v; break;     // U (+y)
        default: x = u; y = 0; z = v; break;    // D (-y)
    }
    // Rotate 90 degrees clockwise seen from above per facing step: relative N (-z) becomes E (+x): (x, z) -> (7 - z, x).
    for (int i = 0, n = ((partFacing % 4) + 4) % 4; i < n; ++i) {
        const int nx = 7 - z, nz = x;
        x = nx;
        z = nz;
    }
    return { partCell.x * 8 + x, partCell.y * 8 + y, partCell.z * 8 + z };
}

// ---- save text ----

std::string Conduit_SaveRuns(const std::vector<ConduitRun>& runs) {
    std::ostringstream os;
    for (const ConduitRun& r : runs) {
        os << "run " << r.id << ' ' << (r.kind == ConduitKind::Cable ? "cable" : "pipe") << ' ' << r.material << ' ' << r.size
           << ' ' << r.points.size();
        for (const CInt3& p : r.points) os << ' ' << p.x << ' ' << p.y << ' ' << p.z;
        os << '\n';
    }
    return os.str();
}

bool Conduit_LoadRuns(const std::string& text, std::vector<ConduitRun>& out, std::string* why) {
    std::vector<ConduitRun> tmp;
    std::istringstream in(text);
    std::string raw;
    int lineNo = 0;
    while (std::getline(in, raw)) {
        ++lineNo;
        const std::string line = Trim(raw);
        if (line.empty() || line[0] == '#') continue;
        auto fail = [&](const std::string& m) { if (why) *why = "line " + std::to_string(lineNo) + ": " + m; return false; };
        const std::vector<std::string> w = Words(line);
        ConduitRun r;
        int n = 0;
        if (w.size() < 6 || w[0] != "run") return fail("expected 'run <id> <kind> <material> <size> <n> ...'");
        if (!ParseInt(w[1], r.id)) return fail("bad id");
        if (w[2] == "cable") r.kind = ConduitKind::Cable;
        else if (w[2] == "pipe") r.kind = ConduitKind::Pipe;
        else return fail("kind must be cable or pipe");
        r.material = w[3];
        if (!ParseInt(w[4], r.size)) return fail("bad size");
        if (!ParseInt(w[5], n) || n < 0) return fail("bad point count");
        if (w.size() != 6 + (size_t)n * 3) return fail("expected " + std::to_string(n) + " points");
        for (int i = 0; i < n; ++i) {
            CInt3 p;
            if (!ParseInt(w[6 + i * 3], p.x) || !ParseInt(w[7 + i * 3], p.y) || !ParseInt(w[8 + i * 3], p.z)) return fail("bad coordinate");
            r.points.push_back(p);
        }
        tmp.push_back(r);
    }
    out = tmp;
    return true;
}

// ---- conversion ----

std::vector<ConduitRun> Conduit_FromCells(const std::vector<CInt3>& cellsIn, ConduitKind kind, const std::string& material, int firstId) {
    std::set<CInt3> cells(cellsIn.begin(), cellsIn.end());
    std::vector<ConduitRun> out;
    int id = firstId;
    auto centre = [](CInt3 c) { return CInt3{ c.x * 8 + 4, c.y * 8 + 4, c.z * 8 + 4 }; };
    for (int axis = 0; axis < 3; ++axis) {
        for (const CInt3& c : cells) { // std::set iterates sorted: deterministic
            if (cells.count(Add(c, axis, -1)) || !cells.count(Add(c, axis, 1))) continue; // not a chain start, or a lone cell on this axis
            CInt3 last = c;
            while (cells.count(Add(last, axis, 1))) last = Add(last, axis, 1);
            ConduitRun r;
            r.id = id++;
            r.kind = kind;
            r.material = material;
            r.size = 1;
            r.points = { centre(c), centre(last) };
            out.push_back(r);
        }
    }
    return out;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"

namespace {
ConduitRun MakeRun(int id, ConduitKind kind, const std::string& mat, int size, std::vector<CInt3> pts) {
    ConduitRun r;
    r.id = id; r.kind = kind; r.material = mat; r.size = size; r.points = std::move(pts);
    return r;
}
ConduitRun Cable(int id, int size, std::vector<CInt3> pts) { return MakeRun(id, ConduitKind::Cable, "copper", size, std::move(pts)); }
ConduitRun Pipe(int id, int size, std::vector<CInt3> pts) { return MakeRun(id, ConduitKind::Pipe, "iron", size, std::move(pts)); }

// A little world: runs by id plus their occupancy.
struct Plot {
    std::map<int, ConduitRun> runs;
    ConduitOccupancy occ;
    void put(const ConduitRun& r) { runs[r.id] = r; occ.add(r); }
};
bool Has(const std::vector<CInt3>& v, CInt3 p) { return std::find(v.begin(), v.end(), p) != v.end(); }
const std::function<bool(CInt3)> kNoSolid = [](CInt3) { return false; };
} // namespace

TEST(Conduit, TheGameFileParses) {
    std::vector<std::string> errors;
    std::ifstream f("resources/data/conduits.txt");
    ASSERT_TRUE(f.good()) << "run from the repo root";
    std::stringstream ss;
    ss << f.rdbuf();
    const ConduitCatalog cat = Conduit_ParseMaterials(ss.str(), &errors);
    for (const std::string& e : errors) ADD_FAILURE() << e;
    ASSERT_NE(cat.find("copper"), nullptr);
    EXPECT_EQ(cat.find("copper")->kind, ConduitKind::Cable);
    EXPECT_EQ(cat.find("copper")->volts, 32);
    ASSERT_NE(cat.find("iron"), nullptr);
    EXPECT_EQ(cat.find("iron")->kind, ConduitKind::Pipe);
    EXPECT_EQ(cat.find("iron")->flow, 40);
    EXPECT_EQ(cat.find("iron")->maxTempC, 250);
    EXPECT_EQ(cat.find("iron")->item, "pipe");
}

TEST(Conduit, BadMaterialLinesError) {
    std::vector<std::string> errors;
    const ConduitCatalog cat = Conduit_ParseMaterials(
        "cable a | item cable | units 8 | volts 32\n"
        "cable b | item cable | units eight\n"
        "wire c | item x\n"
        "pipe d | units 8\n"
        "cable a | item cable\n"
        "cable e | item cable | flow 4\n", &errors);
    EXPECT_EQ(cat.mats.size(), 1u);
    EXPECT_EQ(errors.size(), 5u);
    EXPECT_NE(errors[0].find("line 2"), std::string::npos);
}

TEST(Conduit, ValidRun) {
    std::string why;
    EXPECT_TRUE(Conduit_ValidRun(Cable(1, 1, { {0, 0, 0}, {5, 0, 0}, {5, 0, 4} }), &why)) << why;
    EXPECT_FALSE(Conduit_ValidRun(Cable(1, 1, { {0, 0, 0}, {3, 0, 3} }), &why));      // diagonal
    EXPECT_FALSE(Conduit_ValidRun(Cable(1, 1, { {0, 0, 0}, {0, 0, 0} }), &why));      // zero length
    EXPECT_FALSE(Conduit_ValidRun(Cable(1, 1, { {0, 0, 0} }), &why));                 // single point
    EXPECT_FALSE(Conduit_ValidRun(Cable(1, 0, { {0, 0, 0}, {3, 0, 0} }), &why));      // size 0
    EXPECT_FALSE(Conduit_ValidRun(Cable(1, 5, { {0, 0, 0}, {3, 0, 0} }), &why));      // size 5
    EXPECT_FALSE(Conduit_ValidRun(Cable(1, 1, { {0, 0, 0}, {5, 0, 0}, {2, 0, 0} }), &why)); // reversal
    EXPECT_TRUE(Conduit_ValidRun(Cable(1, 1, { {0, 0, 0}, {5, 0, 0}, {9, 0, 0} }), &why));  // collinear is fine
}

TEST(Conduit, Slots) {
    EXPECT_EQ(Conduit_Slots(Cable(1, 1, { {0, 0, 0}, {2, 0, 0} })).size(), 3u);
    EXPECT_EQ(Conduit_Slots(Cable(1, 2, { {0, 0, 0}, {2, 0, 0} })).size(), 12u);
    // an L of size 1: 4 + 3 segments sharing the corner slot
    const std::vector<CInt3> l = Conduit_Slots(Cable(1, 1, { {0, 0, 0}, {3, 0, 0}, {3, 0, 3} }));
    EXPECT_EQ(l.size(), 7u);
    EXPECT_TRUE(std::is_sorted(l.begin(), l.end()));
    // size 2 covers c and c-1 on the cross axes (lo = -1)
    const std::vector<CInt3> t = Conduit_Slots(Cable(1, 2, { {0, 5, 0}, {0, 5, 1} }));
    EXPECT_TRUE(Has(t, { -1, 4, 0 }));
    EXPECT_TRUE(Has(t, { 0, 5, 1 }));
    EXPECT_FALSE(Has(t, { 1, 5, 0 }));
    // negative slots land in the right cells
    EXPECT_EQ(Conduit_CellOf({ -1, -8, -9 }), (CInt3{ -1, -1, -2 }));
    EXPECT_EQ(Conduit_CellOf({ 7, 8, 0 }), (CInt3{ 0, 1, 0 }));
    ConduitOccupancy occ;
    occ.add(Cable(1, 1, { {-2, 0, 0}, {1, 0, 0} }));
    EXPECT_TRUE(occ.cellUsed({ -1, 0, 0 }));
    EXPECT_TRUE(occ.cellUsed({ 0, 0, 0 }));
    EXPECT_EQ(occ.ownerAt({ -1, 0, 0 }), 1);
}

TEST(Conduit, LengthAndCost) {
    const ConduitMaterial m = *Conduit_ParseMaterials("cable t | item cable | units 8 | volts 1\n", nullptr).find("t");
    EXPECT_EQ(Conduit_Length(Cable(1, 1, { {0, 0, 0}, {3, 0, 0}, {3, 4, 0} })), 7);
    EXPECT_EQ(Conduit_Cost(Cable(1, 1, { {0, 0, 0}, {16, 0, 0} }), m), 2);
    EXPECT_EQ(Conduit_Cost(Cable(1, 2, { {0, 0, 0}, {16, 0, 0} }), m), 8);
    EXPECT_EQ(Conduit_Cost(Cable(1, 1, { {0, 0, 0}, {1, 0, 0} }), m), 1);
}

TEST(Conduit, Occupancy) {
    ConduitOccupancy occ;
    EXPECT_EQ(occ.cellCount(), 0u);
    const ConduitRun a = Cable(1, 1, { {0, 0, 0}, {20, 0, 0} }); // cells 0..2
    const ConduitRun b = Cable(2, 2, { {0, 0, 4}, {3, 0, 4} });
    occ.add(a);
    occ.add(b);
    EXPECT_EQ(occ.ownerAt({ 5, 0, 0 }), 1);
    EXPECT_EQ(occ.ownerAt({ 5, 0, 1 }), 0);
    EXPECT_EQ(occ.ownerAt({ 3, 0, 4 }), 2);
    EXPECT_TRUE(occ.cellUsed({ 2, 0, 0 }));
    EXPECT_FALSE(occ.cellUsed({ 3, 0, 0 }));
    EXPECT_EQ(occ.cellCount(), 4u); // a: 3 cells; b (size 2, y -1..0) adds (0,-1,0)
    occ.remove(a);
    EXPECT_EQ(occ.ownerAt({ 5, 0, 0 }), 0);
    EXPECT_EQ(occ.ownerAt({ 3, 0, 4 }), 2);
    EXPECT_EQ(occ.cellCount(), 2u);
    occ.remove(b);
    EXPECT_EQ(occ.cellCount(), 0u);
    EXPECT_FALSE(occ.cellUsed({ 0, 0, 0 }));
}

TEST(Conduit, CheckPlace) {
    Plot plot;
    plot.put(Cable(1, 1, { {0, 0, 0}, {20, 0, 0} }));
    plot.put(Pipe(3, 1, { {0, 0, 10}, {20, 0, 10} }));
    auto check = [&](const ConduitRun& r, const std::function<bool(CInt3)>& solid = kNoSolid, const std::vector<ConduitPortAt>& ports = {}) {
        return Conduit_CheckPlace(r, plot.occ, plot.runs, solid, ports);
    };

    ConduitPlaceResult r = check(Cable(2, 1, { {0, 0, 5}, {10, 0, 5} }));
    EXPECT_TRUE(r.ok) << r.why;
    EXPECT_TRUE(r.joinsRuns.empty());

    r = check(Cable(2, 1, { {2, 0, 5}, {2, 0, 3}, {-5, 0, 3} }));
    EXPECT_TRUE(r.ok) << r.why;

    // along the existing run, ends outside: the middle overlaps
    r = check(Cable(2, 1, { {-5, 0, 0}, {30, 0, 0} }));
    EXPECT_FALSE(r.ok);
    EXPECT_TRUE(Has(r.blocked, { 5, 0, 0 }));

    // a T: the new run's end lands on the side of run 1
    r = check(Cable(2, 1, { {5, 0, 4}, {5, 0, 0} }));
    EXPECT_TRUE(r.ok) << r.why;
    ASSERT_EQ(r.joinsRuns.size(), 1u);
    EXPECT_EQ(r.joinsRuns[0], 1);

    // an end-to-end joint
    r = check(Cable(2, 1, { {20, 0, 0}, {20, 0, 6} }));
    EXPECT_TRUE(r.ok) << r.why;
    EXPECT_EQ(r.joinsRuns, std::vector<int>({ 1 }));

    // a cable end onto a pipe
    r = check(Cable(2, 1, { {5, 0, 14}, {5, 0, 10} }));
    EXPECT_FALSE(r.ok);
    EXPECT_TRUE(Has(r.blocked, { 5, 0, 10 }));
    EXPECT_TRUE(r.joinsRuns.empty());

    // crossing in the middle
    r = check(Cable(2, 1, { {5, 0, -4}, {5, 0, 4} }));
    EXPECT_FALSE(r.ok);
    EXPECT_TRUE(Has(r.blocked, { 5, 0, 0 }));

    // a corner on the run is a crossing, not a meeting
    r = check(Cable(2, 1, { {5, 0, 4}, {5, 0, 0}, {8, 0, 0} }));
    EXPECT_FALSE(r.ok);

    // into a solid cell
    const std::function<bool(CInt3)> solid = [](CInt3 c) { return c == CInt3{ 0, 0, 1 }; };
    r = check(Cable(2, 1, { {0, 0, 5}, {0, 0, 10} }), solid);
    EXPECT_FALSE(r.ok);
    EXPECT_TRUE(Has(r.blocked, { 0, 0, 8 }));
    EXPECT_FALSE(Has(r.blocked, { 0, 0, 7 }));

    // ending exactly at a fitting port slot in a solid cell (cell 3,0,0 is x 24..31)
    const std::function<bool(CInt3)> solid3 = [](CInt3 c) { return c == CInt3{ 3, 0, 0 }; };
    ConduitPortAt port; port.slot = { 24, 4, 4 }; port.face = 3; port.kind = ConduitPortKind::PowerIn; port.maxSize = 4;
    std::vector<ConduitPortAt> ports = { { { 99, 99, 99 }, 0, ConduitPortKind::PowerIn, 4 }, port };
    r = check(Cable(2, 1, { {10, 4, 4}, {24, 4, 4} }), solid3, ports);
    EXPECT_TRUE(r.ok) << r.why;
    EXPECT_EQ(r.portIdx, std::vector<int>({ 1 }));
    // one slot short of the port is not a connection and not inside the solid cell either: ok but no port
    r = check(Cable(2, 1, { {10, 4, 4}, {23, 4, 4} }), solid3, ports);
    EXPECT_TRUE(r.ok);
    EXPECT_TRUE(r.portIdx.empty());
    // size 2 into a port that takes 2: ok; into one that takes 1: blocked
    ports[1].maxSize = 2;
    r = check(Cable(2, 2, { {10, 4, 4}, {24, 4, 4} }), solid3, ports);
    EXPECT_TRUE(r.ok) << r.why;
    ports[1].maxSize = 1;
    r = check(Cable(2, 2, { {10, 4, 4}, {24, 4, 4} }), solid3, ports);
    EXPECT_FALSE(r.ok);
    EXPECT_TRUE(r.portIdx.empty());
    // a cable into a fluid port
    ports[1].maxSize = 4;
    ports[1].kind = ConduitPortKind::FluidIn;
    r = check(Cable(2, 1, { {10, 4, 4}, {24, 4, 4} }), solid3, ports);
    EXPECT_FALSE(r.ok);
    EXPECT_TRUE(Has(r.blocked, { 24, 4, 4 }));
    r = check(Pipe(4, 1, { {10, 4, 4}, {24, 4, 4} }), solid3, ports);
    EXPECT_TRUE(r.ok) << r.why;
    // an invalid run reports why
    r = check(Cable(2, 1, { {0, 0, 0} }));
    EXPECT_FALSE(r.ok);
    EXPECT_FALSE(r.why.empty());
}

TEST(Conduit, Ports) {
    std::vector<ConduitPortDef> p;
    std::string why;
    ASSERT_TRUE(Conduit_ParsePorts("N:power_in E:power_in S:fluid_out@2,6:max2", p, &why)) << why;
    ASSERT_EQ(p.size(), 3u);
    EXPECT_EQ(p[0].face, 0);
    EXPECT_EQ(p[0].kind, ConduitPortKind::PowerIn);
    EXPECT_EQ(p[0].u, 4);
    EXPECT_EQ(p[0].maxSize, 4);
    EXPECT_EQ(p[1].face, 1);
    EXPECT_EQ(p[2].face, 2);
    EXPECT_EQ(p[2].kind, ConduitPortKind::FluidOut);
    EXPECT_EQ(p[2].u, 2);
    EXPECT_EQ(p[2].v, 6);
    EXPECT_EQ(p[2].maxSize, 2);
    ASSERT_TRUE(Conduit_ParsePorts("U:rotation D:item_in@0,7 W:item_out:max1", p, &why)) << why;
    EXPECT_EQ(p.size(), 3u);
    EXPECT_EQ(p[1].v, 7);
    EXPECT_EQ(p[2].maxSize, 1);

    const char* bad[] = { "X:power_in", "N:steam", "N", "N:power_in@8,0", "N:power_in@1", "N:power_in:max5", "N:power_in:max0", "N:power_in:big2", "NN:power_in" };
    for (const char* b : bad) EXPECT_FALSE(Conduit_ParsePorts(b, p, &why)) << b;
    EXPECT_EQ(p.size(), 3u); // a failed parse leaves the output alone

    EXPECT_EQ(Conduit_DefaultPorts(true, false, false).size(), 6u);
    EXPECT_EQ(Conduit_DefaultPorts(false, true, false).size(), 6u);
    EXPECT_EQ(Conduit_DefaultPorts(false, false, true).size(), 6u);
    EXPECT_EQ(Conduit_DefaultPorts(true, false, true).size(), 12u);
    EXPECT_TRUE(Conduit_DefaultPorts(false, false, false).empty());
    int fluidIn = 0;
    for (const ConduitPortDef& d : Conduit_DefaultPorts(false, false, true)) if (d.kind == ConduitPortKind::FluidIn) { ++fluidIn; EXPECT_TRUE(d.face == 4 || d.face == 0 || d.face == 3); }
    EXPECT_EQ(fluidIn, 3);
}

TEST(Conduit, PortSlot) {
    const CInt3 cell{ 1, 2, -1 }; // slot origin (8, 16, -8)
    ConduitPortDef d;
    auto at = [&](int face, int facing, int u = 4, int v = 4) { d.face = face; d.u = u; d.v = v; return Conduit_PortSlot(cell, facing, d); };
    EXPECT_EQ(at(0, 0), (CInt3{ 11, 20, -8 }));  // N: z layer 0, x = 7 - u
    EXPECT_EQ(at(1, 0), (CInt3{ 15, 20, -5 }));  // E: x layer 7, z = 7 - u
    EXPECT_EQ(at(2, 0), (CInt3{ 12, 20, -1 }));  // S: z layer 7
    EXPECT_EQ(at(3, 0), (CInt3{ 8, 20, -4 }));   // W: x layer 0
    EXPECT_EQ(at(4, 0), (CInt3{ 12, 23, -4 }));  // U: y layer 7
    EXPECT_EQ(at(5, 0), (CInt3{ 12, 16, -4 }));  // D: y layer 0
    EXPECT_EQ(at(0, 0, 2, 6), (CInt3{ 13, 22, -8 }));
    EXPECT_EQ(at(0, 1), (CInt3{ 15, 20, -5 }));  // facing 1: relative N becomes E
    EXPECT_EQ(at(2, 1), (CInt3{ 8, 20, -4 }));   // relative S becomes W
    EXPECT_EQ(at(4, 1, 2, 6), (CInt3{ 9, 23, -6 })); // U rotates in its plane: local (2,7,6) -> (1,7,2)
    EXPECT_EQ(at(1, 4), at(1, 0));
}

TEST(Conduit, SaveLoad) {
    std::vector<ConduitRun> runs = { Cable(3, 2, { {0, 0, 0}, {5, 0, 0}, {5, -7, 0} }), Pipe(9, 1, { {-1, -2, -3}, {-1, -2, 4} }) };
    const std::string text = Conduit_SaveRuns(runs);
    EXPECT_EQ(text, "run 3 cable copper 2 3 0 0 0 5 0 0 5 -7 0\nrun 9 pipe iron 1 2 -1 -2 -3 -1 -2 4\n");
    std::vector<ConduitRun> back;
    std::string why;
    ASSERT_TRUE(Conduit_LoadRuns(text, back, &why)) << why;
    ASSERT_EQ(back.size(), 2u);
    for (size_t i = 0; i < 2; ++i) {
        EXPECT_EQ(back[i].id, runs[i].id);
        EXPECT_EQ(back[i].kind, runs[i].kind);
        EXPECT_EQ(back[i].material, runs[i].material);
        EXPECT_EQ(back[i].size, runs[i].size);
        EXPECT_TRUE(back[i].points == runs[i].points);
    }
    EXPECT_FALSE(Conduit_LoadRuns("run 1 cable copper 1 2 0 0 0 1 0 0\nrun 2 cable copper 1 2 0 0 0 1 0\n", back, &why));
    EXPECT_NE(why.find("line 2"), std::string::npos);
    EXPECT_FALSE(Conduit_LoadRuns("run 1 wire copper 1 2 0 0 0 1 0 0\n", back, &why));
    EXPECT_FALSE(Conduit_LoadRuns("run x cable copper 1 0\n", back, &why));
    EXPECT_FALSE(Conduit_LoadRuns("nonsense\n", back, &why));
}

TEST(Conduit, FromCells) {
    // an L of 5 cells: x chain (0..2, 0, 0) and z chain (2, 0, 0..2)
    const std::vector<ConduitRun> l = Conduit_FromCells({ {0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {2, 0, 1}, {2, 0, 2} }, ConduitKind::Cable, "copper", 10);
    ASSERT_EQ(l.size(), 2u);
    EXPECT_EQ(l[0].id, 10);
    EXPECT_EQ(l[1].id, 11);
    EXPECT_TRUE((l[0].points == std::vector<CInt3>{ {4, 4, 4}, {20, 4, 4} }));
    EXPECT_TRUE((l[1].points == std::vector<CInt3>{ {20, 4, 4}, {20, 4, 20} }));
    EXPECT_EQ(l[0].size, 1);
    EXPECT_EQ(l[0].material, "copper");
    // they meet: the second run's end lies in the first run's slots
    ConduitOccupancy occ;
    occ.add(l[0]);
    EXPECT_EQ(occ.ownerAt(l[1].points[0]), 10);

    // a plus: both arms of x and z through the centre cell (1,0,1): two runs crossing at the shared centre slot
    const std::vector<CInt3> plus = { {1, 0, 1}, {0, 0, 1}, {2, 0, 1}, {1, 0, 0}, {1, 0, 2} };
    const std::vector<ConduitRun> p = Conduit_FromCells(plus, ConduitKind::Pipe, "iron", 1);
    ASSERT_EQ(p.size(), 2u);
    EXPECT_EQ(p[0].kind, ConduitKind::Pipe);
    EXPECT_TRUE((p[0].points == std::vector<CInt3>{ {4, 4, 12}, {20, 4, 12} }));
    EXPECT_TRUE((p[1].points == std::vector<CInt3>{ {12, 4, 4}, {12, 4, 20} }));
    const std::vector<CInt3> sa = Conduit_Slots(p[0]), sb = Conduit_Slots(p[1]);
    EXPECT_TRUE(Has(sa, { 12, 4, 12 }));
    EXPECT_TRUE(Has(sb, { 12, 4, 12 }));

    // an isolated cell gives nothing; a vertical pair gives one run along y
    EXPECT_TRUE(Conduit_FromCells({ {5, 5, 5} }, ConduitKind::Cable, "copper", 1).empty());
    const std::vector<ConduitRun> v = Conduit_FromCells({ {0, 1, 0}, {0, 0, 0}, {0, 1, 0} }, ConduitKind::Cable, "copper", 1);
    ASSERT_EQ(v.size(), 1u);
    EXPECT_TRUE((v[0].points == std::vector<CInt3>{ {4, 4, 4}, {4, 12, 4} }));
}
#endif
