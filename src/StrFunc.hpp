/*
 * Interface to the string-operation functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#pragma once

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace StrFunc {

int split_string(const std::string& str, std::vector<std::string>& out_vec, const std::string& separators = " ,\t;\n");
int split_string_fast(const std::string& str, std::vector<std::string>& out_vec,
                      const std::string& separators = " ,\t;\n");
void to_upper(char* str, int len);
void to_upper(std::string& str);

// match two vectors
void match(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB, std::vector<int>& VecC);
void match_only(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB, std::vector<int>& VecC);

void set_complement(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB,
                    const std::vector<int>& tmp, std::vector<int>& VecC);
bool has_suffix(const std::string& str, const std::string& suffix);
void set_intersect(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB,
                   std::vector<std::string>& VecC);
void set_intersect(const std::vector<int>& VecA, const std::vector<int>& VecB, std::vector<int>& VecC);
void set_complement(const std::vector<int>& toRm, const std::vector<int>& source, std::vector<int>& VecC);

}  // namespace StrFunc
