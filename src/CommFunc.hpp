/*
 * Interface to the commonly-used functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#ifndef _COMMFUNC_H
#define _COMMFUNC_H

#include <cstddef>
#include <cstdint>

#if defined _WIN64 || defined _WIN32
constexpr std::size_t MAX_LINE_SIZE = 0x10000;
#else
constexpr std::size_t MAX_LINE_SIZE = 0x80000;
#endif

// clang-format off
// uint32 + floats
constexpr std::uint32_t DENSE_FILE_TYPE_1 = 0;
// uint32 + uint64_t + uint64_ts + uint32_ts + floats
constexpr std::uint32_t SPARSE_FILE_TYPE_3F = 0x40400000;
// 16*uint32s + uint64_t + uint64_ts + uint32_ts + floats (indicator+samplesize+snpnumber+probenumber+ 6*-9s +valnumber+cols+rowids+betases) [default]
constexpr std::uint32_t SPARSE_FILE_TYPE_3 = 3;
// 16*uint32s + floats (indicator+samplesize+snpnumber+probenumber+ 6*-9s + values) [default]
constexpr std::uint32_t DENSE_FILE_TYPE_3 = 5;
constexpr int RESERVEDUNITS = 16;
constexpr std::size_t FNAMESIZE = 4096;

// clang-format on

#include <cstdio>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace Eigen;

namespace CommFunc {
const double FloatErr = std::numeric_limits<double>::epsilon();
template <typename T>
inline T ABS(T const& a) {
  return (T{} < a) ? a : -a;
}
template <typename T>
extern void free2(T** to) {
  if (*to != nullptr) {
    delete[] *to;
    *to = nullptr;
  }
}
template <typename T>
extern inline std::string atosm(T const& a) {
  if (a == -9) return ("NA");
  std::stringstream ss;
  ss << a;
  return (ss.str());
}
double Abs(const double& x);
bool FloatEqual(double lhs, double rhs);
void FileExist(const std::string& filename);
int max_abs_id(VectorXd& zsxz);
int max_abs_id(std::vector<double>& zsxz);
void getRank_norep(std::vector<int>& a, std::vector<int>& b);
void getUnique(std::vector<std::uint32_t>& a);
void match(const std::vector<std::uint32_t>& VecA, const std::vector<std::uint32_t>& VecB, std::vector<int>& VecC);
static inline unsigned int fputs_checked(const char* ss, FILE* outfile) {
  fputs(ss, outfile);
  return ferror(outfile);
}
void strcpy2(char** to, const std::string& from);
float readfloat(FILE* f);
std::uint64_t readuint64(FILE* f);
std::uint32_t readuint32(FILE* f);
int readint(FILE* f);
double cor(std::vector<double>& y, std::vector<double>& x);
void update_map_kp(const std::vector<std::string>& id_list, std::map<std::string, int>& id_map, std::vector<int>& keep);
void update_map_rm(const std::vector<std::string>& id_list, std::map<std::string, int>& id_map, std::vector<int>& keep);

template <typename T>
inline std::string atos(T const& a) {
  std::stringstream ss;
  ss << a;
  return (ss.str());
}

std::string dtos(double value);
std::string itos(int value);
void update_id_map_kp(const std::vector<std::string>& id_list, std::unordered_map<std::string, int>& id_map,
                      std::vector<int>& keep);
void update_id_map_rm(const std::vector<std::string>& id_list, std::unordered_map<std::string, int>& id_map,
                      std::vector<int>& keep);
void read_indi_list(const std::string& indi_list_file, std::vector<std::string>& indi_list);
void read_msglist(const std::string& msglistfile, std::vector<std::string>& msglist, const std::string& msg);
void progress(int cur, double& disp, int ttl);
}  // namespace CommFunc

#endif
