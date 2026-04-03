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

bool str_within_quto(const std::string& str, std::string& str_buf);
int split_string(const std::string& str, std::vector<std::string>& out_vec, const std::string& separators = " ,\t;\n");
int split_string_fast(const std::string& str, std::vector<std::string>& out_vec,
                      const std::string& separators = " ,\t;\n");
std::string first_string(const std::string& str, const char separator);
std::string last_string(const std::string& str, const char separator);
void to_upper(char* str, int len);
void to_upper(std::string& str);
void to_lower(std::string& str);
std::string get_sub_str(const std::string& rst, size_t pos);
bool StrEqual(const std::string& StrA, const std::string& StrB, bool NoCaseSens = true);
bool StrVecEqual(const std::vector<std::string>& VsBufA, const std::vector<std::string>& VsBufB, int Pos);

// find a string in a string vector ignoring upper or lower case
std::vector<std::string>::iterator find(std::vector<std::string>& target_vs, const std::string& target_str);

// find a char in a string ignoring upper or lower case
std::string::iterator find(std::string& target_str, const char target_ch);

// go to the postion of a give string in a stream ignoring upper or lower case
bool goto_str(std::istream& in_file, const std::string& str);

// rewind a stream
void rewind_if(std::istream& in_file);

// match two vectors
void match(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB, std::vector<int>& VecC);
void match_only(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB, std::vector<int>& VecC);

void set_complement(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB,
                    const std::vector<int>& tmp, std::vector<int>& VecC);
void set_complement(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB,
                    const std::vector<int>& tmp, std::vector<std::uint32_t>& VecC);
int split_string_skip(const std::string& str, std::vector<std::string>& vec_str, std::string separator, int num2skip);
void match_only(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB,
                std::vector<std::uint32_t>& VecC);
bool has_suffix(const std::string& str, const std::string& suffix);
void set_intersect(const std::vector<std::string>& VecA, const std::vector<std::string>& VecB,
                   std::vector<std::string>& VecC);
void set_intersect(const std::vector<int>& VecA, const std::vector<int>& VecB, std::vector<int>& VecC);
void set_complement(const std::vector<int>& toRm, const std::vector<int>& source, std::vector<int>& VecC);
bool stringNumCheck(std::string a, int num);
double pchisqd1(double x);

bool rankContrast(int n, double* Z);

}  // namespace StrFunc
