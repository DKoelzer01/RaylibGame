#ifndef CONDUIT_H
#define CONDUIT_H

// (#1620, conduits infra PR C1; design: #1429 conduits addendum) Conduit runs on a sub-grid: cables and pipes that are
// polylines of 1x1 .. 4x4 slots instead of whole cells. Pure data and logic: ports data, run validation, the slots a run
// occupies, its cost, an occupancy map, a placement check, port slot mapping, save text and the conversion of the old
// cell-based cables. No game code calls it yet (networks are PR C2, the placement UI is C4).
//
// Units: 8 slots per cell per axis. A slot coordinate is in 1/8-cell units relative to the plot origin; cell = floor(slot / 8).

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

struct CInt3 {
    int x = 0, y = 0, z = 0;
    bool operator==(const CInt3& o) const { return x == o.x && y == o.y && z == o.z; }
    bool operator<(const CInt3& o) const {
        if (x != o.x) return x < o.x;
        if (y != o.y) return y < o.y;
        return z < o.z;
    }
};

// The cell holding a slot (floor division: the cell of slot -1 is -1).
CInt3 Conduit_CellOf(CInt3 slot);

// (#1735) Belt: a thin conveyor lying on a cell floor. Its cross-section is kBeltWidth slots across (half a cell) and
// kBeltHeight slots up: the belt's own slot layer plus one open slot above for the items it carries (a quarter-cell
// placement box). Cables and pipes run in the 6 slots above. A belt run's points are its centreline on that floor layer.
enum class ConduitKind : uint8_t { Cable, Pipe, Belt };
constexpr int kBeltWidth = 4;
constexpr int kBeltHeight = 2;

struct ConduitMaterial {
    std::string id;
    ConduitKind kind = ConduitKind::Cable;
    std::string item;          // the item one unit of conduit costs
    int unitsPerItem = 8;      // slots of 1x1 conduit per item (belt: slots of belt length per item)
    int volts = 0;             // cable: EU/t per amp
    int flow = 0;              // pipe: L/t at 1x1
    int maxTempC = 0;          // pipe
    int speed = 0;             // belt: 1/64-cell units an item moves per factory tick
};

struct ConduitCatalog {
    std::vector<ConduitMaterial> mats;
    const ConduitMaterial* find(const std::string& id) const;
};

// Parse; errors appended as "line N: message". Bad lines are skipped, good ones kept.
ConduitCatalog Conduit_ParseMaterials(const std::string& text, std::vector<std::string>* errors);
const ConduitCatalog& Conduit_Materials(); // resources/data/conduits.txt, loaded once (function-local static)

struct ConduitRun {
    int id = 0;
    ConduitKind kind = ConduitKind::Cable;
    std::string material;
    int size = 1;                // 1..4: the cross-section is size x size slots (belt: always kBeltWidth)
    std::vector<CInt3> points;   // centreline corners
};

// >= 2 points; size 1..4; consecutive points differ in exactly one axis; no segment reverses straight back onto the previous.
// A belt also: size kBeltWidth, horizontal segments only, every point on one floor layer (y a multiple of 8) and on the
// quarter-cell grid (x and z even).
bool Conduit_ValidRun(const ConduitRun& run, std::string* why);
// Every occupied slot, no duplicates, sorted. For a centreline slot c: c + (a, b) on the two cross axes, a, b in
// [lo, lo + size - 1] with lo = -(size / 2). A belt: the horizontal cross axis spans [c - 2, c + 1], y spans [c.y, c.y + 1].
std::vector<CInt3> Conduit_Slots(const ConduitRun& run);
int64_t Conduit_Length(const ConduitRun& run);                          // the centreline length in slots
// ceil(length * size * size / unitsPerItem); a belt: ceil(length / unitsPerItem).
int Conduit_Cost(const ConduitRun& run, const ConduitMaterial& mat);

class ConduitOccupancy {
public:
    void add(const ConduitRun& run);     // marks its slots, owner = run id
    void remove(const ConduitRun& run);  // clears the slots this run owns
    int ownerAt(CInt3 slot) const;       // run id, 0 = free
    bool cellUsed(CInt3 cell) const;     // any slot in that cell
    size_t cellCount() const;
private:
    std::map<CInt3, std::array<uint64_t, 8>> mask; // per cell 512 bits
    std::map<CInt3, int> owner;                    // slot -> run id
};

// Ports: a machine face connection. Faces 0..5 = N E S W U D = -z +x +z -x +y -y.
enum class ConduitPortKind { PowerIn, PowerOut, FluidIn, FluidOut, ItemIn, ItemOut, Rotation };

struct ConduitPortDef {
    int face = 0;
    int u = 4, v = 4;            // 0..7 on the face; 4,4 = the face centre
    ConduitPortKind kind = ConduitPortKind::PowerIn;
    int maxSize = 4;
};

// The slot a conduit end must occupy to connect, the slot of the port inside the (solid) part cell.
struct ConduitPortAt {
    CInt3 slot;
    int face = 0;
    ConduitPortKind kind = ConduitPortKind::PowerIn;
    int maxSize = 4;
};

struct ConduitPlaceResult {
    bool ok = false;
    std::string why;
    std::vector<int> joinsRuns;   // ids of the runs met at an end, sorted, unique
    std::vector<int> portIdx;     // indices into the ports argument that an end connects to
    std::vector<CInt3> blocked;   // the slots that stop the placement
};

// An end of the new run may meet (end-to-end or end-to-side) a run of the same kind and material: only the end's own
// size x size cap may overlap that run. Any other overlap, or any overlap with another kind or material, is blocked.
// Solid cells are blocked except the cap of an end that sits exactly on a fitting port's slot.
// Belts (#1735): a belt meets any belt (any material: a side belt merges into the other); its end cap reaches
// kBeltWidth / 2 slots back from the end point, so an end on another belt's centreline overlaps that belt's near half.
// A belt end connects to a side face's item port when the end point is the port's slot moved down to the floor of the
// port's cell (the cross axis within one slot: face centres on N and E are odd slots).
// Every floor-layer slot of a belt needs a solid cell under it (the ground or a part's top), except in a port's cell.
ConduitPlaceResult Conduit_CheckPlace(const ConduitRun& run, const ConduitOccupancy& occ, const std::map<int, ConduitRun>& runs,
                                      const std::function<bool(CInt3 cell)>& solidCell, const std::vector<ConduitPortAt>& ports);

// (#1735 PR 3) The belt end rule above: a belt whose end point is `end` (the next point along it: `next`) meets `port`.
bool Conduit_BeltEndAtPort(CInt3 end, CInt3 next, const ConduitPortAt& port);

// "N:power_in E:power_in S:fluid_out@2,6:max2": whitespace-separated <face letter>:<kind>[@u,v][:max<n>].
bool Conduit_ParsePorts(const std::string& spec, std::vector<ConduitPortDef>& out, std::string* why);
// Face-centre ports. EU consumers: PowerIn on all 6 faces. Generators: PowerOut on all 6. Fluid machines: FluidIn on
// U, N, W and FluidOut on D, S, E. Flags combine.
std::vector<ConduitPortDef> Conduit_DefaultPorts(bool eu, bool generator, bool fluid);
// partFacing 0..3 rotates the faces N E S W around +y (1 = 90 degrees clockwise seen from above: N becomes E).
CInt3 Conduit_PortSlot(CInt3 partCell, int partFacing, const ConduitPortDef& port);

// "run <id> <cable|pipe|belt> <material> <size> <n> x y z x y z ..." one run per line.
std::string Conduit_SaveRuns(const std::vector<ConduitRun>& runs);
bool Conduit_LoadRuns(const std::string& text, std::vector<ConduitRun>& out, std::string* why);

// Old cell-based cables/pipes -> centred 1x1 runs (a cell's centre is slot cell*8 + 4). Maximal straight chains of
// face-adjacent cells along +x, +y, +z (>= 2 cells) become one run from the first cell's centre to the last's. A cell
// with no neighbour at all gives no run. Chains along different axes share the centre slot of their common cell (a joint).
// Order: axis x, y, z, then by the chain's first cell; ids from firstId upward.
std::vector<ConduitRun> Conduit_FromCells(const std::vector<CInt3>& cells, ConduitKind kind, const std::string& material, int firstId);

#endif
