/*
 * Implementations of the commonly-used functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "CommFunc.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

void CommFunc::update_id_map_kp(const std::vector<std::string>& id_list, std::unordered_map<std::string, int>& id_map,
                                std::vector<int>& keep) {
  std::unordered_set<std::string> keep_ids;
  keep_ids.reserve(id_list.size());
  for (const auto& id : id_list) keep_ids.emplace(id);

  for (auto iter = id_map.begin(); iter != id_map.end();) {
    if (keep_ids.find(iter->first) == keep_ids.end()) {
      iter = id_map.erase(iter);
    } else {
      ++iter;
    }
  }

  keep.clear();
  keep.reserve(id_map.size());
  for (const auto& entry : id_map) keep.push_back(entry.second);
  std::stable_sort(keep.begin(), keep.end());
}

void CommFunc::update_id_map_rm(const std::vector<std::string>& id_list, std::unordered_map<std::string, int>& id_map,
                                std::vector<int>& keep) {
  std::unordered_set<std::string> remove_ids;
  remove_ids.reserve(id_list.size());
  for (const auto& id : id_list) remove_ids.emplace(id);

  for (auto iter = id_map.begin(); iter != id_map.end();) {
    if (remove_ids.find(iter->first) != remove_ids.end()) {
      iter = id_map.erase(iter);
    } else {
      ++iter;
    }
  }

  keep.clear();
  keep.reserve(id_map.size());
  for (const auto& entry : id_map) keep.push_back(entry.second);
  std::stable_sort(keep.begin(), keep.end());
}

void CommFunc::read_indi_list(const std::string& indi_list_file, std::vector<std::string>& indi_list) {
  std::ifstream i_indi_list(indi_list_file.c_str());
  if (!i_indi_list) throw std::runtime_error("Error: can not open the file [" + indi_list_file + "] to read.");
  std::string str_buf, id_buf;
  indi_list.clear();
  while (i_indi_list) {
    i_indi_list >> str_buf;
    if (i_indi_list.eof()) break;
    id_buf = str_buf + ":";
    i_indi_list >> str_buf;
    id_buf += str_buf;
    indi_list.push_back(id_buf);
    std::getline(i_indi_list, str_buf);
  }
  i_indi_list.close();
}

void CommFunc::read_msglist(const std::string& msglistfile, std::vector<std::string>& msglist, const std::string& msg) {
  // Read msglist file
  msglist.clear();
  std::string StrBuf;
  std::ifstream i_msglist(msglistfile.c_str());
  if (!i_msglist) throw std::runtime_error("Error: can not open the file [" + msglistfile + "] to read.");
  std::cout << "Reading a list of " << msg << " from [" + msglistfile + "]." << std::endl;
  while (i_msglist >> StrBuf) {
    msglist.push_back(StrBuf);
    std::getline(i_msglist, StrBuf);
  }
  i_msglist.close();
}

std::string CommFunc::dtos(double value) {
  std::stringstream ss;
  ss << std::scientific << value;
  return (ss.str());
}

std::string CommFunc::itos(int value) {
  std::stringstream ss;
  ss << value;
  return (ss.str());
}

double CommFunc::Abs(const double& x) {
  std::complex<double> cld(x);
  double ldAbs = abs(cld);
  return (ldAbs);
}

bool CommFunc::FloatEqual(double lhs, double rhs) {
  if (Abs(lhs - rhs) < FloatErr) return true;
  return false;
}

void CommFunc::FileExist(const std::string& filename) {
  std::ifstream ifile(filename.c_str());
  if (!ifile) throw std::runtime_error("Error: can not open the file [" + filename + "] to read.");
}

int CommFunc::max_abs_id(std::vector<double>& zsxz) {
  int id = 0;
  double tmpVal, cmpVal = abs(zsxz[0]);
  for (int i = 1; i < zsxz.size(); i++) {
    tmpVal = abs(zsxz[i]);
    if (cmpVal - tmpVal < 1e-6) {
      cmpVal = tmpVal;
      id = i;
    }
  }
  return id;
}

int CommFunc::max_abs_id(VectorXd& zsxz) {
  int id = 0;
  double tmpVal, cmpVal = fabs(zsxz[0]);
  for (int i = 1; i < zsxz.size(); i++) {
    tmpVal = fabs(zsxz[i]);
    if (cmpVal - tmpVal < 1e-6) {
      cmpVal = tmpVal;
      id = i;
    }
  }
  return (id);
}

void CommFunc::getRank_norep(std::vector<int>& a, std::vector<int>& b) {
  b.resize(a.size());
  std::map<int, int> rep_chck;
  long mapsize = 0;
  for (int i = static_cast<int>(a.size()) - 1; i >= 0; i--) {
    int count = 0;
    for (int j = 0; j < a.size(); j++)
      if (a[j] < a[i]) count++;
    rep_chck.insert(std::pair<int, int>(count, i));
    while (rep_chck.size() == mapsize) {
      count++;
      rep_chck.insert(std::pair<int, int>(count, i));
    }
    mapsize = rep_chck.size();
    b[i] = count;
  }
}

void CommFunc::getUnique(std::vector<std::uint32_t>& a) {
  std::sort(a.begin(), a.end());
  auto it = std::unique(a.begin(), a.end());
  a.erase(it, a.end());
}

void CommFunc::match(const std::vector<std::uint32_t>& VecA, const std::vector<std::uint32_t>& VecB,
                     std::vector<int>& VecC) {
  std::unordered_map<std::uint32_t, int> id_map;
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

void CommFunc::strcpy2(char** to, const std::string& from) {
  char* tmp = new char[from.size() + 1];
  std::copy(from.begin(), from.end(), tmp);
  tmp[from.size()] = '\0';
  *to = tmp;
}

float CommFunc::readfloat(FILE* f) {
  float v;
  if (fread(reinterpret_cast<void*>(&v), sizeof(v), 1, f) != 1)
    throw std::runtime_error("ERROR: unexpected end of file or read error while reading float data.");
  return v;
}

std::uint64_t CommFunc::readuint64(FILE* f) {
  std::uint64_t v;
  if (fread(reinterpret_cast<void*>(&v), sizeof(v), 1, f) != 1)
    throw std::runtime_error("ERROR: unexpected end of file or read error while reading uint64 data.");
  return v;
}

std::uint32_t CommFunc::readuint32(FILE* f) {
  std::uint32_t v;
  if (fread(reinterpret_cast<void*>(&v), sizeof(v), 1, f) != 1)
    throw std::runtime_error("ERROR: unexpected end of file or read error while reading uint32 data.");
  return v;
}

int CommFunc::readint(FILE* f) {
  int v;
  if (fread(reinterpret_cast<void*>(&v), sizeof(v), 1, f) != 1)
    throw std::runtime_error("ERROR: unexpected end of file or read error while reading int data.");
  return v;
}

double CommFunc::cor(std::vector<double>& y, std::vector<double>& x) {
  long N = x.size();
  if (N != y.size() || N < 1) {
    printf("Error: The lengths of x and y do not match.\n");
    exit(EXIT_FAILURE);
  }

  int i = 0;
  double d_buf = 0.0, y_mu = 0.0, x_mu = 0.0, x_var = 0.0, y_var = 0.0, cov = 0.0;
  for (i = 0; i < N; i++) {
    x_mu += x[i];
    y_mu += y[i];
  }
  x_mu /= static_cast<double>(N);
  y_mu /= static_cast<double>(N);
  for (i = 0; i < N; i++) {
    d_buf = (x[i] - x_mu);
    x_var += d_buf * d_buf;
    d_buf = (y[i] - y_mu);
    y_var += d_buf * d_buf;
  }
  x_var /= (N - 1.0);
  y_var /= (N - 1.0);
  for (i = 0; i < N; i++) cov += (x[i] - x_mu) * (y[i] - y_mu);
  cov /= static_cast<double>(N - 1);
  double a = 0.0, b = 0.0, sse = 0.0, a_se = 0.0, b_se = 0.0, r = 0.0;
  if (x_var > 0.0) b = cov / x_var;
  a = y_mu - b * x_mu;
  for (i = 0; i < N; i++) {
    d_buf = y[i] - a - b * x[i];
    sse += d_buf * d_buf;
  }
  if (x_var > 0.0) {
    a_se = sqrt((sse / (N - 2.0)) * (1.0 / N + x_mu * x_mu / (x_var * (N - 1.0))));
    b_se = sqrt(sse / x_var / (N - 1.0) / (N - 2.0));
  }
  if (x_var > 0.0 && y_var > 0.0) {
    r = cov / sqrt(y_var * x_var);
  }

  return (r);
}

void CommFunc::update_map_kp(const std::vector<std::string>& id_list, std::map<std::string, int>& id_map,
                             std::vector<int>& keep) {
  int i = 0;
  std::map<std::string, int> id_map_buf(id_map);
  for (i = 0; i < id_list.size(); i++) id_map_buf.erase(id_list[i]);
  std::map<std::string, int>::iterator iter;
  for (iter = id_map_buf.begin(); iter != id_map_buf.end(); iter++) id_map.erase(iter->first);

  keep.clear();
  for (iter = id_map.begin(); iter != id_map.end(); iter++) keep.push_back(iter->second);
  std::stable_sort(keep.begin(), keep.end());
}

void CommFunc::update_map_rm(const std::vector<std::string>& id_list, std::map<std::string, int>& id_map,
                             std::vector<int>& keep) {
  int i = 0;
  for (i = 0; i < id_list.size(); i++) id_map.erase(id_list[i]);

  keep.clear();
  std::map<std::string, int>::iterator iter;
  for (iter = id_map.begin(); iter != id_map.end(); iter++) keep.push_back(iter->second);
  std::stable_sort(keep.begin(), keep.end());
}

void CommFunc::progress(int& cur, double& disp, int ttl) {
  double desti = 1.0 * cur / (ttl - 1);
  if (desti >= disp) {
    printf("%3.0f%%\r", 100.0 * desti);
    fflush(stdout);
    if (disp == 0) disp += 0.05;
    else if (disp == 0.05) disp += 0.2;
    else if (disp == 0.25) disp += 0.5;
    else disp += 0.25;
  }
}
