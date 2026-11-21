#include "itemdata.h"

#include <cstdlib>
#include <vector>

namespace {
bool Bad(char c) { return c == ' ' || c == '\t' || c == ';' || c == '=' || c == '|' || c == '\n' || c == '\r'; }

struct Pair { std::string key, value; };
std::vector<Pair> Parse(const std::string& data) {
    std::vector<Pair> out;
    size_t pos = 0;
    while (pos < data.size()) {
        size_t end = data.find(';', pos);
        if (end == std::string::npos) end = data.size();
        const std::string item = data.substr(pos, end - pos);
        const size_t eq = item.find('=');
        if (eq != std::string::npos && eq > 0) out.push_back({ item.substr(0, eq), item.substr(eq + 1) });
        pos = end + 1;
    }
    return out;
}
std::string Join(const std::vector<Pair>& pairs) {
    std::string out;
    for (const Pair& p : pairs) { if (!out.empty()) out += ';'; out += p.key + '=' + p.value; }
    return out;
}
}

namespace ItemData {

std::string Sanitize(const std::string& value) {
    std::string s = value;
    for (char& c : s) if (Bad(c)) c = '_';
    return s;
}

std::string Get(const std::string& data, const std::string& key) {
    for (const Pair& p : Parse(data)) if (p.key == key) return p.value;
    return std::string();
}

bool Has(const std::string& data, const std::string& key) {
    for (const Pair& p : Parse(data)) if (p.key == key) return true;
    return false;
}

int GetInt(const std::string& data, const std::string& key, int fallback) {
    const std::string v = Get(data, key);
    if (v.empty()) return fallback;
    char* end = nullptr;
    const long n = strtol(v.c_str(), &end, 10);
    return end && *end == '\0' ? (int)n : fallback;
}

std::string Set(const std::string& data, const std::string& key, const std::string& value) {
    std::vector<Pair> pairs = Parse(data);
    const std::string k = Sanitize(key), v = Sanitize(value);
    for (size_t i = 0; i < pairs.size(); ++i)
        if (pairs[i].key == k) {
            if (v.empty()) pairs.erase(pairs.begin() + (long)i); else pairs[i].value = v;
            return Join(pairs);
        }
    if (!v.empty() && !k.empty()) pairs.push_back({ k, v });
    return Join(pairs);
}

std::string SetInt(const std::string& data, const std::string& key, int value) { return Set(data, key, std::to_string(value)); }

bool Valid(const std::string& data) {
    if (data.empty()) return true;
    size_t pos = 0;
    while (pos <= data.size()) {
        size_t end = data.find(';', pos);
        if (end == std::string::npos) end = data.size();
        const std::string item = data.substr(pos, end - pos);
        const size_t eq = item.find('=');
        if (eq == std::string::npos || eq == 0 || eq + 1 >= item.size()) return false;
        for (char c : item) if (c != '=' && Bad(c)) return false;
        if (item.find('=', eq + 1) != std::string::npos) return false;
        if (end == data.size()) break;
        pos = end + 1;
    }
    return true;
}

} // namespace ItemData

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1446): the instance data record ----
#include "unit_helpers.h"

TEST(ItemData, SetGetAndRemove) {
    std::string d;
    d = ItemData::Set(d, "h", "bronze");
    d = ItemData::Set(d, "g", "oak");
    d = ItemData::SetInt(d, "q", 104);
    EXPECT_EQ(d, "h=bronze;g=oak;q=104");
    EXPECT_EQ(ItemData::Get(d, "g"), "oak");
    EXPECT_EQ(ItemData::GetInt(d, "q"), 104);
    EXPECT_EQ(ItemData::GetInt(d, "w", 7), 7) << "absent: the fallback";
    d = ItemData::Set(d, "h", "steel");
    EXPECT_EQ(ItemData::Get(d, "h"), "steel") << "replaced in place";
    d = ItemData::Set(d, "g", "");
    EXPECT_FALSE(ItemData::Has(d, "g")) << "an empty value removes the key";
    EXPECT_EQ(d, "h=steel;q=104");
    EXPECT_TRUE(ItemData::Valid(d));
}

TEST(ItemData, ValuesAreSanitizedForSavesAndProtocol) {
    const std::string d = ItemData::Set("", "name", "my sword;x=1|y\n");
    EXPECT_EQ(ItemData::Get(d, "name"), "my_sword_x_1_y_");
    EXPECT_TRUE(ItemData::Valid(d));
    EXPECT_EQ(d.find(' '), std::string::npos);
    EXPECT_FALSE(ItemData::Valid("a b=c"));
    EXPECT_FALSE(ItemData::Valid("a=b=c"));
    EXPECT_FALSE(ItemData::Valid("=x"));
    EXPECT_FALSE(ItemData::Valid("a="));
    EXPECT_TRUE(ItemData::Valid(""));
    EXPECT_EQ(ItemData::GetInt("q=1x", "q", -1), -1) << "not a number: the fallback";
}
#endif
