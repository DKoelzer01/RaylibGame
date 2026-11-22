#include "cavekinds.h"

#include <cstdio>
#include <cstring>
#include <string>

namespace {
std::string F(float v) { // a GLSL float literal that reads back as exactly v
    char b[48];
    snprintf(b, sizeof b, "%.9g", (double)v);
    std::string s(b);
    if (s.find_first_of(".eEn") == std::string::npos) s += ".0";
    return s;
}
template <typename Get> std::string Array(const char* type, const char* name, Get get) {
    std::string s = std::string("const ") + type + " " + name + "[" + std::to_string(kCaveKindCount) + "] = " + type + "[](";
    for (int i = 0; i < kCaveKindCount; ++i) s += (i ? ", " : "") + get(kCaveKinds[i]);
    return s + ");\n";
}
} // namespace

std::string CaveKinds_Glsl() {
    std::string s = "// (generated from src/engine/cavekinds.h, CaveKinds_Glsl)\n";
    s += "const int kCaveKindCount = " + std::to_string(kCaveKindCount) + ";\n";
    s += Array("int", "kCK_solid", [](const CaveKindDef& k) { return std::to_string(k.solid); });
    s += Array("float", "kCK_outline", [](const CaveKindDef& k) { return F(k.outline); });
    s += Array("int", "kCK_lopsided", [](const CaveKindDef& k) { return std::to_string(k.lopsided); });
    s += Array("float", "kCK_rubble", [](const CaveKindDef& k) { return F(k.rubble); });
    s += Array("float", "kCK_fixedLevel", [](const CaveKindDef& k) { return F(k.fixedLevel); });
    s += Array("float", "kCK_floorR", [](const CaveKindDef& k) { return F(k.floorR); });
    s += Array("int", "kCK_tilt", [](const CaveKindDef& k) { return std::to_string(k.tilt); });
    s += Array("int", "kCK_scale", [](const CaveKindDef& k) { return std::to_string(k.scaleToPassage); });
    s += Array("float", "kCK_intrude", [](const CaveKindDef& k) { return F(k.intrude); });
    s += Array("float", "kCK_outcropCap", [](const CaveKindDef& k) { return F(k.outcropCap); });
    // CaveKind_Index, as in cavekinds.h
    s += "int caveKindIndex(float k) { int i = int(k + 0.5); return k < 0.5 ? 0 : i >= kCaveKindCount ? kCaveKindCount - 1 : i; }\n";
    return s;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"
TEST(CaveKinds, TodaysKindsAndTheirIndex) {
    ASSERT_EQ(kCaveKindCount, 5);
    EXPECT_STREQ(kCaveKinds[(int)CaveKind::Passage].name, "passage");
    EXPECT_STREQ(kCaveKinds[(int)CaveKind::Mouth].name, "mouth");
    EXPECT_STREQ(kCaveKinds[(int)CaveKind::Solid].name, "solid");
    EXPECT_STREQ(kCaveKinds[(int)CaveKind::Chamber].name, "chamber");
    EXPECT_STREQ(kCaveKinds[(int)CaveKind::Cavern].name, "cavern");
    // the old inline tests: mouth 0.5 < k < 1.5, solid 1.5 < k < 2.5, hall k > 2.5, else passage
    EXPECT_EQ(CaveKind_Index(0.0f), 0);
    EXPECT_EQ(CaveKind_Index(1.0f), 1);
    EXPECT_EQ(CaveKind_Index(2.0f), 2);
    EXPECT_EQ(CaveKind_Index(3.0f), 3);
    EXPECT_EQ(CaveKind_Index(4.0f), 4) << "a cavern";
    EXPECT_EQ(CaveKind_Index(7.0f), 4) << "past the last: the last";
    EXPECT_EQ(CaveKind_Index(-1.0f), 0);
}
TEST(CaveKinds, TheGlslBlockReadsBackExactly) {
    const std::string g = CaveKinds_Glsl();
    EXPECT_NE(g.find("const float kCK_outline[5] = float[](0.200000003, 0.400000006, 0.200000003, 0.200000003, 0.200000003);"), std::string::npos) << g;
    EXPECT_NE(g.find("const float kCK_fixedLevel[5] = float[](-1.0, 0.449999988, -1.0, -1.0, -1.0);"), std::string::npos) << g;
    EXPECT_NE(g.find("const int kCK_solid[5] = int[](0, 0, 1, 0, 0);"), std::string::npos) << g;
}
#endif
