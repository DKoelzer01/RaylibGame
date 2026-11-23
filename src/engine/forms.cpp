#include "forms.h"
#include "utils/logger.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

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
// "<n> = <form> <n>, <form> <n>": the output count and the inputs (forms of the same material; "item:<id>" a fixed item)
bool ParseSides(const std::string& s, int& outCount, std::vector<std::pair<std::string, int>>& in) {
    const size_t eq = s.find('=');
    if (eq == std::string::npos) return false;
    outCount = std::max(1, atoi(Trim(s.substr(0, eq)).c_str()));
    for (const std::string& part : Split(s.substr(eq + 1), ',')) {
        std::istringstream ps(part);
        std::string form; int n = 0;
        if (ps >> form >> n && n > 0) in.push_back({ form, n });
    }
    return !in.empty();
}
struct Rule { std::string shape, machine; int outCount = 1; std::vector<std::pair<std::string, int>> in; int ticks = 0, eu = 0; };

std::string ReadFile(const char* path) {
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}
std::vector<std::string> gShapes;
} // namespace

std::vector<FormRecipe> Forms_Generate(const std::string& formsText, const std::string& materialsText, std::vector<std::string>* shapesOut) {
    // forms.txt: "hand <shape> | <n> = <form> <n>, ..."  or  "machine <shape> | <machine> | <n> = <form> <n>, ... | ticks | eu"
    std::vector<Rule> rules;
    std::set<std::string> shapes;
    {
        std::istringstream in(formsText);
        for (std::string line; std::getline(in, line);) {
            line = Trim(line);
            if (line.empty() || line[0] == '#') continue;
            const auto c = Split(line, '|');
            std::istringstream head(c[0]);
            std::string kind, shape;
            head >> kind >> shape;
            Rule r;
            r.shape = shape;
            bool ok = false;
            if (kind == "hand" && c.size() >= 2) ok = ParseSides(c[1], r.outCount, r.in);
            else if (kind == "machine" && c.size() >= 5) {
                r.machine = c[1];
                ok = ParseSides(c[2], r.outCount, r.in);
                r.ticks = std::max(1, atoi(c[3].c_str()));
                r.eu = std::max(0, atoi(c[4].c_str()));
            }
            if (!ok || shape.empty()) { LOG_WARN("Forms", "bad line '%s'\n", line.c_str()); continue; }
            shapes.insert(shape);
            rules.push_back(r);
        }
    }
    if (shapesOut) shapesOut->assign(shapes.begin(), shapes.end()); // (only the caller's: the game's table keeps its own, tests pass their text)
    // materials.txt: "<id> | <Name> | r,g,b | <forms>"
    std::vector<FormRecipe> out;
    std::istringstream in(materialsText);
    for (std::string line; std::getline(in, line);) {
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        const auto c = Split(line, '|');
        if (c.size() < 4) continue;
        std::set<std::string> has;
        std::istringstream fs(c[3]);
        for (std::string f; fs >> f;) has.insert(f);
        for (const Rule& r : rules) {
            if (!has.count(r.shape)) continue;
            FormRecipe fr;
            fr.machine = r.machine;
            fr.output = c[0] + "_" + r.shape;
            fr.outCount = r.outCount;
            fr.ticks = r.ticks; fr.eu = r.eu;
            bool inputsExist = true;
            for (const auto& [form, n] : r.in) {
                if (form.rfind("item:", 0) == 0) { fr.inputs.push_back({ form.substr(5), n }); continue; }
                if (!has.count(form)) { inputsExist = false; break; } // (a material without the input's form: no recipe)
                fr.inputs.push_back({ c[0] + "_" + form, n });
            }
            if (inputsExist) out.push_back(fr);
        }
    }
    return out;
}

const std::vector<FormRecipe>& Forms_All() {
    static const std::vector<FormRecipe> all = Forms_Generate(ReadFile("resources/data/forms.txt"), ReadFile("resources/data/materials.txt"), &gShapes);
    return all;
}
const std::vector<std::string>& Forms_Shapes() { Forms_All(); return gShapes; }

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"

TEST(Forms, OneLinePerShapeGivesEveryMetalItsRecipe) {
    const auto r = Forms_Generate(
        "hand plate | 1 = ingot 2\n"
        "machine plate | mechanical_press | 1 = ingot 1 | 100 | 0\n"
        "hand screw | 2 = rod 1\n",
        "iron  | Iron  | 1,2,3 | ingot plate rod screw\n"
        "zinc  | Zinc  | 1,2,3 | ingot plate\n"
        "brass | Brass | 1,2,3 | ingot screw\n"); // brass has screws but no rods: no screw recipe
    auto has = [&](const std::string& machine, const std::string& out) {
        for (const FormRecipe& f : r) if (f.machine == machine && f.output == out) return &f;
        return (const FormRecipe*)nullptr;
    };
    ASSERT_NE(has("", "iron_plate"), nullptr);
    EXPECT_EQ(has("", "iron_plate")->inputs[0], (std::pair<std::string, int>{ "iron_ingot", 2 }));
    ASSERT_NE(has("mechanical_press", "zinc_plate"), nullptr);
    EXPECT_EQ(has("mechanical_press", "zinc_plate")->ticks, 100);
    ASSERT_NE(has("", "iron_screw"), nullptr);
    EXPECT_EQ(has("", "iron_screw")->outCount, 2);
    EXPECT_EQ(has("", "brass_screw"), nullptr) << "no rod form: no screw recipe";
    EXPECT_EQ(has("", "zinc_screw"), nullptr) << "no screw form";
    EXPECT_GE(Forms_Shapes().size(), 6u) << "a test's generate leaves the game's shapes alone (it emptied them on Windows' test order)";
}
#endif
