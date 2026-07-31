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

  bool look = false;
  std::string str_buf;
  std::string symbol_pool =
      "`1234567890-=~!@#$%^&*()_+qwertyuiop[]\\asdfghjkl;'zxcvbnm,./QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>? \t\n";
  std::string::size_type pos;

  // Remove seperators
  for (char separator : separators) {
    pos = symbol_pool.find(separator);
    if (pos != std::string::npos) symbol_pool.erase(symbol_pool.begin() + pos);
  }

  for (char i : str) {
    if (symbol_pool.find(i) != std::string::npos) {
      if (!look) look = true;
      str_buf += i;
    } else {
      if (look) {
        look = false;
        out_vec.push_back(str_buf);
        str_buf.erase(str_buf.begin(), str_buf.end());
      }
    }
  }
  if (look) out_vec.push_back(str_buf);

  return static_cast<int>(out_vec.size());
}

/// Both `str` and `separators` are consist of ASCII characters.
int StrFunc::split_string_fast(const std::string& str, std::vector<std::string>& out_vec,
                               const std::string& separators) {
  if (str.empty()) return 0;
  out_vec.clear();

  // heuristic, reduce reallocation
  out_vec.reserve(str.size() / 4);

  bool is_sep[256] = {};
  for (unsigned char c : separators) {
    is_sep[c] = true;
  }

  // Record the start position and length of each substring
  size_t start_pos = 0;
  size_t len = 0;

  // Traverse the input string
  for (size_t i = 0; i < str.size(); ++i) {
    // If the current character is a separator, add the previous substring to the output vector
    if (is_sep[static_cast<unsigned char>(str[i])]) {
      if (len > 0) {
        out_vec.emplace_back(str, start_pos, len);
        len = 0;
      }
    } else {
      // If the current character is not a separator, update the start position and length of the substring
      if (len == 0) {
        start_pos = i;
      }
      ++len;
    }
  }

  // Add the last substring
  if (len > 0) {
    out_vec.emplace_back(str, start_pos, len);
  }

  return out_vec.size();
}

void StrFunc::to_upper(char* str, int len) {
  int i = 0;
  for (i = 0; i < len; i++) {
    if (str[i] >= 'a' && str[i] <= 'z') str[i] += 'A' - 'a';
  }
}

// Uppercase ASCII, avoid `std::to_upper` locale table lookup.
void StrFunc::to_upper(std::string& str) {
  size_t i = 0;
  for (i = 0; i < str.size(); i++) {
    if (str[i] >= 'a' && str[i] <= 'z') str[i] += 'A' - 'a';
  }
}

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
  std::unordered_set<std::string> id_set(VecB.begin(), VecB.end());
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
