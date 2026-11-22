#include "composed.h"
#include "itemdata.h"
#include "inventory.h"
#include "utils/logger.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace {
std::unordered_map<std::string, ComposedClass> gClasses;
std::unordered_map<std::string, PartMaterial> gMaterials;
bool gLoaded = false;

std::string Trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && isspace((unsigned char)s[a])) ++a;
    while (b > a && isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}
std::vector<std::string> Split(const std::string& s, char c) {
    std::vector<std::string> out; std::string cur;
    for (char ch : s) { if (ch == c) { out.push_back(Trim(cur)); cur.clear(); } else cur += ch; }
    out.push_back(Trim(cur));
    return out;
}
void EnsureLoaded() { if (!gLoaded) Composed_Load(); }

ToolType ToolTypeOf(const std::string& t) {
    if (t == "pickaxe") return ToolType::Pickaxe;
    if (t == "shovel") return ToolType::Shovel;
    if (t == "axe") return ToolType::Axe;
    if (t == "hoe") return ToolType::Hoe;
    return ToolType::Hands;
}
// The usual handle and binding for a head (what an old tool converts to): a wooden handle; fibre lashing up to iron,
// leather for bronze, an iron binding for steel and better.
std::string DefaultBinding(const std::string& head) {
    const PartMaterial* m = Composed_Material(head);
    if (head == "bronze") return "leather";
    if (m && m->tier >= 3) return "iron";
    return "fibre";
}
}

bool Composed_Load(const char* classesPath, const char* materialsPath) {
    gLoaded = true;
    gClasses.clear(); gMaterials.clear();
    bool ok = true;
    {
        std::ifstream f(classesPath);
        if (!f) { LOG_ERROR("Composed", "%s: not found\n", classesPath); ok = false; }
        std::string line;
        while (std::getline(f, line)) {
            line = Trim(line);
            if (line.empty() || line[0] == '#') continue;
            const auto p = Split(line, '|');
            if (p.size() < 6) continue;
            gClasses[p[0]] = ComposedClass{ p[0], p[1], p[2], p[3], p[4], p[5] };
        }
    }
    {
        std::ifstream f(materialsPath);
        if (!f) { LOG_ERROR("Composed", "%s: not found\n", materialsPath); ok = false; }
        std::string line;
        while (std::getline(f, line)) {
            line = Trim(line);
            if (line.empty() || line[0] == '#') continue;
            const auto p = Split(line, '|');
            if (p.size() < 9) continue;
            PartMaterial m;
            m.id = p[0]; m.name = p[1];
            const auto rgb = Split(p[2], ',');
            if (rgb.size() == 3) { m.r = (unsigned char)atoi(rgb[0].c_str()); m.g = (unsigned char)atoi(rgb[1].c_str()); m.b = (unsigned char)atoi(rgb[2].c_str()); }
            m.tier = atoi(p[3].c_str());
            m.hardness = (float)atof(p[4].c_str()); m.toughness = (float)atof(p[5].c_str()); m.density = (float)atof(p[6].c_str());
            m.sockets = atoi(p[7].c_str()); m.factor = (float)atof(p[8].c_str());
            gMaterials[m.id] = m;
        }
    }
    return ok && !gClasses.empty() && !gMaterials.empty();
}

bool Composed_Is(const std::string& id) { EnsureLoaded(); return gClasses.count(id) > 0; }
const ComposedClass* Composed_Class(const std::string& id) { EnsureLoaded(); auto it = gClasses.find(id); return it == gClasses.end() ? nullptr : &it->second; }
const PartMaterial* Composed_Material(const std::string& id) { EnsureLoaded(); auto it = gMaterials.find(id); return it == gMaterials.end() ? nullptr : &it->second; }

std::string Composed_MakeData(const std::string& head, const std::string& handle, const std::string& binding, int qualityPct) {
    std::string d;
    d = ItemData::Set(d, "h", head);
    d = ItemData::Set(d, "g", handle.empty() ? std::string("wood") : handle);
    d = ItemData::Set(d, "b", binding.empty() ? DefaultBinding(head) : binding);
    if (qualityPct != 100) d = ItemData::SetInt(d, "q", qualityPct);
    return d;
}

bool Composed_IsGem(const std::string& id) { return id == "ruby" || id == "sapphire" || id == "emerald" || id == "amethyst" || id == "quartz"; }

std::string Composed_GemText(const std::string& id) { // (#1219, combat doc: GUESS effects)
    if (id == "ruby") return "Ruby (burns)";
    if (id == "sapphire") return "Sapphire (chills)";
    if (id == "emerald") return "Emerald (+25% durability)";
    if (id == "amethyst") return "Amethyst (shocks)";
    if (id == "quartz") return "Quartz (+50% poise damage)";
    return id;
}

std::string Composed_FromParts(const std::string& cls, const std::vector<std::string>& ingredients, int qualityPct) {
    const ComposedClass* c = Composed_Class(cls);
    if (!c) return {};
    std::string headPart = c->headPart;
    for (char& ch : headPart) if (ch == ' ') ch = '_';
    std::string head, handle, binding;
    auto materialOf = [](const std::string& id, const std::string& part) -> std::string { // "<material>_<part>" with a known material
        if (id.size() <= part.size() + 1 || id.compare(id.size() - part.size(), part.size(), part) != 0 || id[id.size() - part.size() - 1] != '_') return {};
        const std::string m = id.substr(0, id.size() - part.size() - 1);
        return Composed_Material(m) ? m : std::string();
    };
    for (const std::string& id : ingredients) {
        std::string m;
        if (head.empty() && !(m = materialOf(id, headPart)).empty()) head = m;
        else if (handle.empty() && !(m = materialOf(id, "handle")).empty()) handle = m;
        else if (binding.empty() && !(m = materialOf(id, "binding")).empty()) binding = m;
    }
    if (head.empty()) return {};
    return Composed_MakeData(head, handle, binding, qualityPct);
}

ComposedStats Composed_Stats(const std::string& id, const std::string& data) {
    ComposedStats s;
    const ComposedClass* c = Composed_Class(id);
    if (!c) return s;
    s.cls = id;
    s.head = ItemData::Get(data, "h");
    s.handle = ItemData::Get(data, "g");
    s.binding = ItemData::Get(data, "b");
    const PartMaterial* head = Composed_Material(s.head);
    if (!head) return s;
    const PartMaterial* handle = Composed_Material(s.handle.empty() ? std::string("wood") : s.handle);
    const PartMaterial* binding = Composed_Material(s.binding.empty() ? DefaultBinding(s.head) : s.binding);
    s.valid = true;
    s.quality = std::clamp(ItemData::GetInt(data, "q", 100), 90, 110) / 100.0f;
    s.wear = std::max(0, ItemData::GetInt(data, "w", 0));
    const float hf = handle ? handle->factor : 1.0f, bf = binding ? binding->factor : 1.0f;
    s.tool.type = ToolTypeOf(c->toolType);
    s.tool.tier = s.tool.type == ToolType::Hands ? 0 : std::clamp(head->tier, 1, 4);
    s.durability = (int)std::lround(25.0f * head->toughness * hf * bf * s.quality);
    s.damage = (3.0f + 3.0f * head->hardness) * s.quality;
    // Dig speed against the tier's reference hardness (the tier's own number): flint bites a little faster than stone.
    s.digSpeed = head->tier > 0 ? std::clamp(head->hardness / (float)head->tier, 0.5f, 1.5f) : 1.0f;
    s.weight = 0.4f * head->density + 0.3f * (handle ? handle->density : 0.7f) + 0.05f * (binding ? binding->density : 0.3f);
    s.sockets = std::max(head->sockets, binding ? binding->sockets : 0);
    { // (#1219) gems set in its sockets: "j=ruby,sapphire", no more than it has sockets; an emerald: x1.25 durability (GUESS)
        const std::string j = ItemData::Get(data, "j");
        size_t a = 0;
        while (a < j.size() && (int)s.gems.size() < s.sockets) {
            size_t b = j.find(',', a); if (b == std::string::npos) b = j.size();
            const std::string g = j.substr(a, b - a);
            if (Composed_IsGem(g)) s.gems.push_back(g);
            a = b + 1;
        }
        for (const std::string& g : s.gems) if (g == "emerald") s.durability = (int)std::lround(s.durability * 1.25f);
    }
    std::string noun = c->noun;
    std::string mat = head->name;
    s.name = mat + " " + noun;
    return s;
}

bool Composed_FromLegacy(const std::string& oldId, std::string& cls, std::string& data) {
    EnsureLoaded();
    const size_t u = oldId.find('_');
    if (u == std::string::npos) return false;
    const std::string mat = oldId.substr(0, u), type = oldId.substr(u + 1);
    const std::string c = "tool_" + type;
    if (!gClasses.count(c) || !gMaterials.count(mat)) return false;
    if (gMaterials[mat].tier < 1) return false; // wood_..., fibre_...: not a head
    cls = c;
    data = Composed_MakeData(mat);
    return true;
}

int Composed_ConvertInventory(Inventory& inv, std::map<std::string, int>* wearUsed) {
    int n = 0;
    auto& cells = inv.grid();
    for (int i = 0; i < Inventory::kCells; ++i) {
        InvSlot& s = cells[i];
        if (s.empty() || s.hasData()) continue;
        std::string cls, data;
        if (!Composed_FromLegacy(s.id, cls, data) || !itemDB.find(cls)) continue;
        if (wearUsed) {
            auto it = wearUsed->find(s.id);
            if (it != wearUsed->end()) { if (it->second > 0) data = ItemData::SetInt(data, "w", it->second); wearUsed->erase(it); }
        }
        const std::string old = s.id;
        const int64_t extra = s.centi - 100; // tools hold one a cell; any more (an odd save) become copies of their own
        s.id = cls; s.data = data; s.centi = 100;
        ++n;
        for (int64_t k = extra; k >= 100; k -= 100) { ++inv.muted; inv.addInstance(cls, Composed_MakeData(ItemData::Get(data, "h"))); --inv.muted; ++n; }
        (void)old;
    }
    return n;
}

std::string Composed_DisplayName(const std::string& id, const std::string& data) {
    if (!Composed_Is(id) || data.empty()) return std::string();
    const ComposedStats s = Composed_Stats(id, data);
    return s.valid ? s.name : std::string();
}

std::string Composed_PartsText(const std::string& id, const std::string& data) {
    if (!Composed_Is(id) || data.empty()) return std::string();
    const ComposedStats s = Composed_Stats(id, data);
    if (!s.valid) return std::string();
    auto nameOf = [](const std::string& m) { const PartMaterial* p = Composed_Material(m); return p ? p->name : m; };
    std::string t = "Head: " + nameOf(s.head) + "   Handle: " + nameOf(s.handle.empty() ? "wood" : s.handle) +
                    "   Binding: " + nameOf(s.binding.empty() ? std::string("fibre") : s.binding);
    char buf[160];
    if (s.tool.type != ToolType::Hands) { snprintf(buf, sizeof buf, "\nTier %d (%s)", s.tool.tier, ToolTierName(s.tool.tier)); t += buf; }
    snprintf(buf, sizeof buf, "\nDurability %d / %d", std::max(0, s.durability - s.wear), s.durability); t += buf;
    if (s.quality != 1.0f) { snprintf(buf, sizeof buf, "   Quality %+.0f%%", (s.quality - 1.0f) * 100.0f); t += buf; }
    if (s.sockets > 0) {
        snprintf(buf, sizeof buf, "   Sockets %d/%d", (int)s.gems.size(), s.sockets); t += buf; // (#1219) set / free
        if (!s.gems.empty()) {
            t += "\nGems:";
            for (const std::string& g : s.gems) t += " " + Composed_GemText(g);
        }
    }
    return t;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1449): composed items' stats from their parts ----
#include "unit_helpers.h"

TEST(Composed, DefaultPartsReproduceTodaysTools) {
    ASSERT_TRUE(Composed_Load());
    // Every old tool converts to its class with the default parts, and keeps its type, tier and durability.
    for (const char* t : { "pickaxe", "shovel", "axe", "hoe" })
        for (const char* m : { "stone", "bronze", "iron", "steel", "titanium" }) {
            const std::string old = std::string(m) + "_" + t;
            std::string cls, data;
            ASSERT_TRUE(Composed_FromLegacy(old, cls, data)) << old;
            EXPECT_EQ(cls, std::string("tool_") + t);
            const ComposedStats s = Composed_Stats(cls, data);
            ASSERT_TRUE(s.valid) << old;
            const ToolSpec was = ToolOf(old);
            EXPECT_EQ((int)s.tool.type, (int)was.type) << old;
            EXPECT_EQ(s.tool.tier, was.tier) << old;
            EXPECT_EQ(s.durability, ToolDurability(was)) << old << ": the same strokes as before";
            EXPECT_FLOAT_EQ(s.digSpeed, 1.0f) << old;
        }
}

TEST(Composed, StatsFollowTheParts) {
    ASSERT_TRUE(Composed_Load());
    const ComposedStats stone = Composed_Stats("tool_pickaxe", Composed_MakeData("stone"));
    const ComposedStats flint = Composed_Stats("tool_pickaxe", Composed_MakeData("flint"));
    EXPECT_EQ(stone.name, "Stone pickaxe");
    EXPECT_EQ(flint.tool.tier, 1);
    EXPECT_GT(flint.digSpeed, stone.digSpeed) << "flint bites faster";
    EXPECT_LT(flint.durability, stone.durability) << "but wears sooner";
    // A tougher binding and handle last longer; quality moves it a few percent either way.
    const ComposedStats plain = Composed_Stats("tool_axe", Composed_MakeData("bronze", "wood", "leather"));
    const ComposedStats bound = Composed_Stats("tool_axe", Composed_MakeData("bronze", "bronze", "bronze"));
    EXPECT_GT(bound.durability, plain.durability);
    const ComposedStats good = Composed_Stats("tool_axe", Composed_MakeData("bronze", {}, {}, 110));
    EXPECT_NEAR(good.durability, plain.durability * 1.1f, 1.0f);
    EXPECT_NEAR(good.damage, plain.damage * 1.1f, 1e-3f);
    EXPECT_EQ(Composed_Stats("tool_axe", Composed_MakeData("bronze", {}, {}, 300)).quality, 1.1f) << "clamped to 110 %";
    // Sockets: bronze one, steel two; weapons have no tool type.
    EXPECT_EQ(plain.sockets, 1);
    const ComposedStats sword = Composed_Stats("weapon_sword", Composed_MakeData("steel"));
    EXPECT_EQ(sword.sockets, 2);
    EXPECT_EQ((int)sword.tool.type, (int)ToolType::Hands);
    EXPECT_EQ(sword.name, "Steel sword");
    EXPECT_GT(sword.damage, Composed_Stats("weapon_sword", Composed_MakeData("bronze")).damage);
    // Wear is read from the record.
    EXPECT_EQ(Composed_Stats("tool_pickaxe", ItemData::SetInt(Composed_MakeData("iron"), "w", 37)).wear, 37);
}

// (#1449 part 2) Old tools in an inventory become composed ones, their wear moved into the first copy of each.
TEST(Composed, ConvertInventoryMovesOldToolsAndTheirWear) {
    UnitTest::EnsureItemDB();
    ASSERT_TRUE(Composed_Load());
    Inventory inv;
    inv.add("stone_pickaxe", 1);
    inv.add("iron_axe", 1);
    inv.add("stone", 4);
    std::map<std::string, int> wear = { { "stone_pickaxe", 37 }, { "bronze_axe", 5 } };
    EXPECT_EQ(Composed_ConvertInventory(inv, &wear), 2);
    EXPECT_EQ(inv.count("stone_pickaxe"), 0);
    EXPECT_EQ(inv.countInstances("tool_pickaxe"), 1);
    EXPECT_EQ(inv.countInstances("tool_axe"), 1);
    EXPECT_EQ(inv.count("stone"), 4) << "everything else stays";
    const int p = inv.findInstance("tool_pickaxe", ItemData::SetInt(Composed_MakeData("stone"), "w", 37));
    EXPECT_GE(p, 0) << "the old pickaxe's 37 strokes are in its record";
    EXPECT_EQ(wear.count("stone_pickaxe"), 0u) << "its per-id wear went";
    EXPECT_EQ(wear.count("bronze_axe"), 1u) << "wear for a tool not held stays (nothing to move it into)";
    // The converted copies are tools with their old tier and durability.
    EXPECT_EQ(ToolOf("tool_axe", Composed_MakeData("iron")).tier, 2);
    EXPECT_EQ(ToolDurabilityOf("tool_axe", Composed_MakeData("iron")), 400);
    EXPECT_EQ(ToolDurabilityOf("stone_axe", ""), 150) << "old ids still answer as before";
    EXPECT_EQ(Composed_ConvertInventory(inv, &wear), 0) << "converting again does nothing";
}

TEST(Composed, AssemblyFromPartItems) {
    ASSERT_TRUE(Composed_Load());
    const std::string d = Composed_FromParts("tool_pickaxe", { "bronze_pick_head", "oak_planks", "wood_handle", "leather_binding" });
    EXPECT_EQ(ItemData::Get(d, "h"), "bronze");
    EXPECT_EQ(ItemData::Get(d, "g"), "wood");
    EXPECT_EQ(ItemData::Get(d, "b"), "leather");
    EXPECT_EQ(Composed_FromParts("tool_pickaxe", { "steel_pick_head" }), Composed_MakeData("steel")) << "missing parts: the head's defaults";
    EXPECT_EQ(ItemData::Get(Composed_FromParts("weapon_sword", { "steel_blade", "bone_handle" }), "g"), "bone");
    EXPECT_EQ(ItemData::Get(Composed_FromParts("weapon_spear", { "flint_spearhead" }), "h"), "flint");
    EXPECT_EQ(Composed_FromParts("tool_pickaxe", { "bronze_blade", "wood_handle" }), "") << "a sword's blade is no pick head";
    EXPECT_EQ(Composed_FromParts("tool_pickaxe", { "unobtainium_pick_head" }), "") << "an unknown material";
    EXPECT_EQ(Composed_FromParts("stone", { "bronze_pick_head" }), "") << "not a class";
    EXPECT_EQ(ItemData::GetInt(Composed_FromParts("tool_axe", { "iron_axe_head" }, 106), "q", 100), 106) << "the anvil's quality";
}

TEST(Composed, NamesAndPartsForTheUi) {
    ASSERT_TRUE(Composed_Load());
    const std::string d = ItemData::SetInt(Composed_MakeData("bronze", {}, {}, 105), "w", 40);
    EXPECT_EQ(Composed_DisplayName("tool_pickaxe", d), "Bronze pickaxe");
    EXPECT_EQ(Composed_DisplayName("tool_pickaxe", ""), "") << "no record: the item table's name";
    EXPECT_EQ(Composed_DisplayName("stone", "h=bronze"), "") << "not a class";
    const std::string t = Composed_PartsText("tool_pickaxe", d);
    EXPECT_NE(t.find("Head: Bronze"), std::string::npos) << t;
    EXPECT_NE(t.find("Binding: Leather"), std::string::npos) << t;
    EXPECT_NE(t.find("Tier 2"), std::string::npos) << t;
    EXPECT_NE(t.find("Durability 380 / 420"), std::string::npos) << t; // 400 x 1.05 quality, 40 used
    EXPECT_NE(t.find("Quality +5%"), std::string::npos) << t;
    EXPECT_NE(t.find("Sockets 0/1"), std::string::npos) << t; // (#1219: set / sockets)
}

TEST(Composed, UnknownsAreNotValid) {
    ASSERT_TRUE(Composed_Load());
    EXPECT_FALSE(Composed_Stats("stone", "h=stone").valid) << "not a class";
    EXPECT_FALSE(Composed_Stats("tool_pickaxe", "h=cheese").valid) << "not a material";
    EXPECT_FALSE(Composed_Stats("tool_pickaxe", "").valid) << "no head";
    std::string c, d;
    EXPECT_FALSE(Composed_FromLegacy("stone", c, d));
    EXPECT_FALSE(Composed_FromLegacy("wood_pickaxe", c, d)) << "wood isn't a head material";
    EXPECT_FALSE(Composed_FromLegacy("stone_spoon", c, d));
    EXPECT_TRUE(Composed_Is("tool_pickaxe"));
    EXPECT_FALSE(Composed_Is("stone_pickaxe"));
}
#endif
