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
#include <cstdio>
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
  for (size_t i = 0; i < separators.size(); i++) {
    pos = symbol_pool.find(separators[i]);
    if (pos != std::string::npos) symbol_pool.erase(symbol_pool.begin() + pos);
  }

  for (size_t i = 0; i < str.size(); i++) {
    if (symbol_pool.find(str[i]) != std::string::npos) {
      if (!look) look = true;
      str_buf += str[i];
    } else {
      if (look) {
        look = false;
        out_vec.push_back(str_buf);
        str_buf.erase(str_buf.begin(), str_buf.end());
      }
    }
  }
  if (look) out_vec.push_back(str_buf);

  return (int)out_vec.size();
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

std::string StrFunc::first_string(const std::string& str, const char separator) {
  int pos = (int)str.find(separator);
  if (pos != -1) return std::string(str.begin(), str.begin() + pos);
  return std::string("");
}

std::string StrFunc::last_string(const std::string& str, const char separator) {
  int pos = (int)str.find_last_of(separator);
  if (pos != -1) return std::string(str.begin() + pos + 1, str.end());
  return std::string("");
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

void StrFunc::to_lower(std::string& str) {
  size_t i = 0;
  for (i = 0; i < str.size(); i++) {
    if (str[i] >= 'A' && str[i] <= 'Z') str[i] -= 'A' - 'a';
  }
}

std::string StrFunc::get_sub_str(const std::string& rst, size_t pos) {
  std::vector<std::string> vs_buf;
  StrFunc::split_string(rst, vs_buf);
  return vs_buf[pos];
}

bool StrFunc::StrEqual(const std::string& StrA, const std::string& StrB, bool NoCaseSens) {
  if (!NoCaseSens) return StrA == StrB;
  std::string StrBufA = StrA, StrBufB = StrB;
  to_upper(StrBufA);
  to_upper(StrBufB);
  return StrBufA == StrBufB;
}

bool StrFunc::StrVecEqual(const std::vector<std::string>& VsBufA, const std::vector<std::string>& VsBufB, int Pos) {
  int SizeA = (int)VsBufA.size(), SizeB = (int)VsBufB.size();
  if (SizeA != SizeB) return false;
  if (Pos >= SizeA) throw("Invalid Pos! StrFunc::StrVecEqual");

  int i = 0;
  for (i = Pos; i < SizeA; i++) {
    if (VsBufA[i] != VsBufB[i]) return false;
  }

  return true;
}

bool StrFunc::str_within_quto(const std::string& str, std::string& str_buf) {
  unsigned int begin = str.find_first_of("\"");
  unsigned int end = str.find_last_of("\"");
  if (begin == std::string::npos || end == std::string::npos || begin == end) return false;

  str_buf = "";
  str_buf.insert(str_buf.begin(), str.begin() + begin + 1, str.begin() + end);
  return true;
}

std::vector<std::string>::iterator StrFunc::find(std::vector<std::string>& target_vs, const std::string& target_str) {
  std::string str_buf = target_str;
  std::vector<std::string> vs_buf = target_vs;

  size_t i = 0;
  for (i = 0; i < vs_buf.size(); i++) to_upper(vs_buf[i]);
  to_upper(str_buf);
  return target_vs.begin() + (std::find(vs_buf.begin(), vs_buf.end(), str_buf) - vs_buf.begin());
}

std::string::iterator StrFunc::find(std::string& target_str, const char target_ch) {
  char ch_buf = target_ch;
  std::string str_buf = target_str;
  to_upper(str_buf);
  if (ch_buf > 'a' && ch_buf < 'z') ch_buf += 'A' - 'a';
  return target_str.begin() + (std::find(str_buf.begin(), str_buf.end(), ch_buf) - str_buf.begin());
}

bool StrFunc::goto_str(std::istream& in_file, const std::string& str) {
  std::string str_buf;
  std::string query_str = str;
  std::vector<std::string> vs_buf;
  StrFunc::to_upper(query_str);
  while (in_file >> str_buf) {
    if (StrFunc::split_string(str_buf, vs_buf) > 0) {
      str_buf = vs_buf[0];
    } else {
      continue;
    }

    StrFunc::to_upper(str_buf);
    if (str_buf == "#") {
      std::getline(in_file, str_buf);
      continue;
    }
    if (str_buf == query_str) return true;
  }

  return false;
}

void StrFunc::rewind_if(std::istream& in_file) {
  in_file.clear(std::ios::goodbit);
  in_file.seekg(std::ios::beg);
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

void StrFunc::match_only(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB,
                         std::vector<std::uint32_t>& VecC) {
  std::unordered_map<std::string_view, std::uint32_t> id_map;
  id_map.reserve(VecB.size());
  VecC.clear();
  VecC.reserve(VecA.size());
  for (size_t i = 0; i < VecB.size(); i++) id_map.emplace(VecB[i], static_cast<std::uint32_t>(i));
  for (const auto& value : VecA) {
    auto iter = id_map.find(value);
    if (iter != id_map.end()) VecC.push_back(iter->second);
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

void StrFunc::set_complement(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB,
                             const std::vector<int>& tmp, std::vector<std::uint32_t>& VecC) {
  std::unordered_set<std::string_view> to_remove(VecA.begin(), VecA.end());

  VecC.clear();
  VecC.reserve(VecB.size());
  for (size_t i = 0; i < VecB.size(); i++)
    if (to_remove.find(VecB[i]) == to_remove.end()) VecC.push_back(tmp[i]);
}

// form head
int StrFunc::split_string_skip(const std::string& str, std::vector<std::string>& vec_str, std::string separator,
                               int num2skip) {
  if (str.empty()) return 0;
  vec_str.clear();

  bool look = false;
  std::string str_buf;
  int count = 0;

  for (size_t i = 0; i < str.size(); i++) {
    if (separator.find(str[i]) == std::string::npos) {
      if (!look) look = true;
      str_buf += str[i];
    } else {
      if (look) {
        look = false;
        if (++count > num2skip) vec_str.push_back(str_buf);
        str_buf.erase(str_buf.begin(), str_buf.end());
      }
    }
  }
  if (look) vec_str.push_back(str_buf);

  return (int)vec_str.size();
}

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

bool StrFunc::stringNumCheck(std::string a, int num) {
  std::vector<std::string> b;
  int number = split_string(a, b, ", \t\n");
  return (number == num);
}
