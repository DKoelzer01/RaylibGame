#ifndef ITEMDATA_H
#define ITEMDATA_H

// (#1446) Per-item instance data: what one weapon (later armour) is made of and has been through - parts, anvil quality,
// sockets, wear - as a short "key=value;key=value" record on its inventory cell (InvSlot::data). Empty for plain items.
// Items with data never stack. Keys are short (h head, g grip/handle, b binding, q quality %, s sockets, w wear used).
// Values are sanitized: no spaces, ';', '=', '|' or line breaks (the record must survive the save's and protocol's
// whitespace- and '|'-separated fields).

#include <string>

namespace ItemData {
std::string Get(const std::string& data, const std::string& key);            // "" if absent
int GetInt(const std::string& data, const std::string& key, int fallback = 0);
bool Has(const std::string& data, const std::string& key);
std::string Set(const std::string& data, const std::string& key, const std::string& value); // "" value removes the key
std::string SetInt(const std::string& data, const std::string& key, int value);
std::string Sanitize(const std::string& value);                               // what Set stores for a value
bool Valid(const std::string& data);                                          // only key=value pairs of sanitized text
}

#endif
