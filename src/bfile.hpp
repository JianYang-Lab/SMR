//
//  bfile.hpp
//  SMR_CPP
//
//  Created by Futao Zhang on 5/07/2018.
//  Copyright © 2018 Futao Zhang. All rights reserved.
//

#pragma once

#include <sys/stat.h>
#include <sys/types.h>

#include <map>
#include <string>
#include <vector>

#if defined _WIN64 || defined _WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#ifndef __APPLE__
#include <omp.h>
#endif

#include <Eigen/Eigen>
#include <map>
#include <string>
#include <vector>

using namespace Eigen;

namespace SMRDATA {
// implementation of SoA(Struct of Array)
struct bInfo {
  // bim file
  int _autosome_num;
  std::vector<int> _chr;
  std::vector<std::string> _snp_name;
  std::map<std::string, int> _snp_name_map;  // snp name to snp index
  // std::unordered_map<std::string,int> _snp_name_map;

  std::vector<double> _genet_dst;
  std::vector<int> _bp;
  std::vector<std::string> _allele1;
  std::vector<std::string> _allele2;
  std::vector<std::string> _ref_A;    // reference allele
  std::vector<std::string> _other_A;  // the other allele
  int _snp_num;
  std::vector<double> _rc_rate;
  std::vector<int> _include;  // snps indices, initialized in the read_bimfile()
  VectorXd _maf;

  // fam file
  std::vector<std::string> _fid;
  std::vector<std::string> _pid;
  std::map<std::string, int> _id_map;  // "<fam_id>:<iid>" to individual index
  std::vector<std::string> _fa_id;
  std::vector<std::string> _mo_id;
  std::vector<int> _sex;
  std::vector<double> _pheno;
  int _indi_num;
  std::vector<int> _keep;  // individual indices, initialized in the read_famfile()
  MatrixXd _varcmp_Py;     // BLUP solution to the total genetic effects of individuals

  // bed file
  std::vector<std::vector<bool>> _snp_1;
  std::vector<std::vector<bool>> _snp_2;

  // imputed data
  bool _dosage_flag = false;
  std::vector<std::vector<float>> _geno_dose;
  std::vector<double> _impRsq;

  // genotypes
  MatrixXf _geno;

  std::vector<double> _mu;
};

struct ldInfo {
  std::vector<int> _esi_chr;
  std::vector<std::string> _esi_rs;
  std::vector<int> _esi_gd;
  std::vector<int> _esi_bp;
  std::vector<std::string> _esi_allele1;
  std::vector<std::string> _esi_allele2;
  std::vector<int> _esi_include;
  std::map<std::string, int> _snp_name_map;
  // std::unordered_map<std::string,int> _snp_name_map;

  std::vector<float> _esi_freq;

  std::vector<std::uint64_t> _cols;
  std::vector<float> _val;

  std::uint64_t _snpNum;
  std::uint64_t _valNum;
};

void filter_snp_maf(bInfo* bdata, double maf);
void calcu_mu(bInfo* bdata, bool ssq_flag = false);
void ld_report(char* outFileName, char* bFileName, char* indilstName, char* indilst2remove, char* snplstName,
               char* snplst2exclde, int chr, char* rs, double maf, bool ldr, bool ldr2, int ldWind);
void calcu_mean_rsq(char* outFileName, char* bFileName, char* indilstName, char* indilst2remove, char* snplstName,
                    char* snplst2exclde, int chr, double maf, bool ldr, bool ldr2, int ldWind, double rsq_cutoff);
void lookup(char* outFileName, char* bldFileName, char* snplstName, char* snplst2exclde, int chr, char* snprs,
            char* snprs2exclde, char* fromsnprs, char* tosnprs, int snpWind, bool snpWindflg, int fromsnpkb,
            int tosnpkb, int ld_wind);
void read_ld_esifile(ldInfo* ldinfo, char* esiFileName);
void ld_esi_man(ldInfo* ldinfo, char* snplstName, char* snplst2exclde, int chr, char* snprs, char* fromsnprs,
                char* tosnprs, int snpWind, bool snpwindFlag, int fromsnpkb, int tosnpkb, char* snprs2exclde);
void fetch_ld_by_id(ldInfo* ldinfo, FILE* ldfprt, int sid, std::vector<float>& ld);
void fetch_ld_by_id(ldInfo* ldinfo, FILE* ldfprt, std::vector<std::uint32_t>& curId, int sid, std::vector<float>& ld);
void fetch_ld_by_snps(ldInfo* ldinfo, FILE* ldfprt, std::string rs, std::vector<float>& ld);
void ld_calc_o2m(VectorXd& ld_v, long targetid, MatrixXd& X, bool centered = false);
}  // namespace SMRDATA
