/*
 * Implementations of the string-operation functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "StrFunc.hpp"

#include <algorithm>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

int StrFunc::split_string(const std::string& str, std::vector<std::string>& out_vec, const std::string& separators) {
  if (str.empty()) return 0;
  out_vec.clear();
  out_vec.reserve(str.size() / 4);  // heuristic, reduce reallocation

  // A char is content iff it belongs to the printable-ASCII pool (plus '\t', '\n') and is
  // not listed in separators; anything else (control chars incl. '\r', DEL, bytes >= 0x80)
  // acts as a separator, so CRLF-terminated files parse cleanly.
  bool is_content[256] = {};
  for (unsigned char c : std::string_view(
           "`1234567890-=~!@#$%^&*()_+qwertyuiop[]\\asdfghjkl;'zxcvbnm,./QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>? \t\n"))
    is_content[c] = true;
  for (unsigned char c : separators) is_content[c] = false;

  size_t start_pos = 0;
  size_t len = 0;
  for (size_t i = 0; i < str.size(); ++i) {
    auto c = static_cast<unsigned char>(str[i]);
    if (!is_content[c]) {
      if (len > 0) {
        out_vec.emplace_back(str, start_pos, len);
        len = 0;
      }
    } else {
      if (len == 0) start_pos = i;
      ++len;
    }
  }
  if (len > 0) out_vec.emplace_back(str, start_pos, len);

  return static_cast<int>(out_vec.size());
}

void StrFunc::to_upper(char* str, int len) {
  int i = 0;
  for (i = 0; i < len; i++) {
    if (str[i] >= 'a' && str[i] <= 'z') str[i] += 'A' - 'a';
  }
}

// Uppercase ASCII, avoid `std::to_upper` locale table lookup.
void StrFunc::to_upper(std::string& str) { to_upper(str.data(), static_cast<int>(str.size())); }

void StrFunc::match(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB,
                    std::vector<int>& VecC) {
  std::unordered_map<std::string_view, int> id_map;
  id_map.reserve(VecB.size());
  VecC.clear();
  VecC.reserve(VecA.size());
  for (size_t i = 0; i < VecB.size(); i++) id_map.emplace(VecB[i], static_cast<int>(i));
  for (const auto& value : VecA) {
    auto iter = id_map.find(value);
    if (iter == id_map.end()) VecC.push_back(-9);
    else VecC.push_back(iter->second);
  }
}

/// Get intersection of `VecA` and `VecB`. `VecC` store the index of intersection in `VecB`.
void StrFunc::match_only(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB,
                         std::vector<int>& VecC) {
  std::unordered_map<std::string_view, int> id_map;
  id_map.reserve(VecB.size());
  VecC.clear();
  VecC.reserve(VecA.size());
  for (size_t i = 0; i < VecB.size(); i++) id_map.emplace(VecB[i], static_cast<int>(i));
  for (const auto& value : VecA) {
    auto iter = id_map.find(value);
    if (iter != id_map.end()) VecC.push_back(iter->second);
  }
}

void StrFunc::set_complement(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB,
                             const std::vector<int>& tmp, std::vector<int>& VecC) {
  std::unordered_set<std::string_view> to_remove(VecA.begin(), VecA.end());

  VecC.clear();
  VecC.reserve(VecB.size());
  for (size_t i = 0; i < VecB.size(); i++)
    if (to_remove.find(VecB[i]) == to_remove.end()) VecC.push_back(tmp[i]);
}

// form head

bool StrFunc::has_suffix(const std::string& str, const std::string& suffix) {
  return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void StrFunc::set_intersect(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB,
                            std::vector<std::string>& VecC) {
  std::unordered_set<std::string_view> id_set(VecB.begin(), VecB.end());
  VecC.clear();
  VecC.reserve(std::min(VecA.size(), VecB.size()));
  for (const auto& value : VecA)
    if (id_set.find(value) != id_set.end()) VecC.push_back(value);
}

void StrFunc::set_intersect(const std::vector<int>& VecA, const std::vector<int>& VecB, std::vector<int>& VecC) {
  std::unordered_set<int> id_set(VecB.begin(), VecB.end());
  VecC.clear();
  VecC.reserve(std::min(VecA.size(), VecB.size()));
  for (const auto& value : VecA)
    if (id_set.find(value) != id_set.end()) VecC.push_back(value);
}

void StrFunc::set_complement(const std::vector<int>& toRm, const std::vector<int>& source, std::vector<int>& VecC) {
  std::unordered_set<int> to_remove(toRm.begin(), toRm.end());

  VecC.clear();
  VecC.reserve(source.size());
  for (const auto& value : source)
    if (to_remove.find(value) == to_remove.end()) VecC.push_back(value);
}
