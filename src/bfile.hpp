//
//  bfile.hpp
//  SMR_CPP
//
//  Created by Futao Zhang on 5/07/2018.
//  Copyright © 2018 Futao Zhang. All rights reserved.
//

#pragma once

#include <cstdint>
#include <cstdio>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#if defined _WIN64 || defined _WIN32
  #include <direct.h>
#else
  #include <unistd.h>
#endif

#ifndef __APPLE__
  #include <omp.h>
#endif

#include <Eigen/Eigen>

using namespace Eigen;

namespace SMRDATA {
// PLINK binary genotype data (bed/bim/fam) in Struct-of-Arrays layout:
// every vector below is indexed by SNP index k or individual index i.
struct bInfo {
  // bim file
  int _autosome_num;      // number of autosomes (22 for human); 1..22 = autosomes, 23 = X, 24 = Y, 25 = XY, 26 = MT
  std::vector<int> _chr;  // chromosome code per SNP (PLINK coding above)
  std::vector<std::string> _snp_name;
  std::unordered_map<std::string, int> _snp_name_map;  // SNP rs ID -> SNP index

  std::vector<double> _genet_dst;  // genetic distance (bim col 3); passed through to output files, not used in analysis
  std::vector<int> _bp;            // base-pair position per SNP
  std::vector<std::string> _allele1;  // A1 from the .bim file: the coded (effect) allele,
                                      // unchanged after load (used for output)
  std::vector<std::string> _allele2;  // A2 from the .bim file: the other allele, unchanged after load
  std::vector<std::string> _ref_A;    // reference allele (dosages are counted for this allele);
                                      // initialized as _allele1, swapped with _other_A by allele_check
                                      // to match the summary data's coded allele
  std::vector<std::string> _other_A;  // the allele paired with _ref_A; initialized as _allele2,
                                      // swapped together with _ref_A
  int _snp_num;                       // total number of SNPs in the .bim file
  std::vector<int> _include;          // indices of SNPs retained for analysis (subset of [0, _snp_num), maintained by
                                      // extract/exclude/filter ops)
  VectorXd _maf;                      // allele frequencies per SNP (used by plot paths; analysis uses _mu/2 instead)

  // fam file
  std::vector<std::string> _fid;                 // family ID per individual (fam col 1)
  std::vector<std::string> _pid;                 // individual ID per individual (fam col 2)
  std::unordered_map<std::string, int> _id_map;  // "<fam_id>:<iid>" to individual index
  std::vector<std::string> _fa_id;               // paternal ID per individual (fam col 3)
  std::vector<std::string> _mo_id;               // maternal ID per individual (fam col 4)
  std::vector<int> _sex;       // sex code per individual (fam col 5): 1 = male, 2 = female, other = unknown; used for X
                               // dosage weighting in calcu_mu
  std::vector<double> _pheno;  // phenotype value per individual (fam col 6); parsed but not used by SMR analyses
  int _indi_num;               // total number of individuals in the .fam file
  std::vector<int> _keep;      // individual indices, initialized in the read_famfile()

  // bed file
  std::vector<std::vector<bool>> _snp_1;  // packed genotype bitplane 1: _snp_1[k][i]
  std::vector<std::vector<bool>>
      _snp_2;  // packed genotype bitplane 2: _snp_2[k][i]; (snp1,snp2): 00 = BB, 01 = het, 11 = AA, 10 = missing

  // imputed data
  std::vector<std::vector<float>> _geno_dose;  // imputed dosage matrix [snp][individual] (populated by the imputed-data
                                               // reader; unused by current analysis paths)
  std::vector<double> _impRsq;                 // imputation quality score (R^2) per SNP

  // the mean dosage of SNP k over kept, non-missing individuals
  std::vector<double> _mu;

  bool containsSNP(const std::string& snp) const { return _snp_name_map.find(snp) != _snp_name_map.end(); }
  bool is_autosome(int chr) const { return chr >= 1 && chr <= _autosome_num; }
  bool is_x(int chr) const { return chr == _autosome_num + 1; }
};

// LD panel data (.bld + .esi) in Struct-of-Arrays layout.
struct ldInfo {
  std::vector<int> _esi_chr;         // chromosome code per SNP (from .esi)
  std::vector<std::string> _esi_rs;  // SNP rs ID
  std::vector<int> _esi_gd;          // genetic distance (esi col 3); passed through to outputs, not used in analysis
  std::vector<int> _esi_bp;          // base-pair position per SNP
  std::vector<std::string> _esi_allele1;  // coded allele (A1) from .esi
  std::vector<std::string> _esi_allele2;  // the other allele (A2) from .esi
  std::vector<int>
      _esi_include;  // indices of SNPs retained in the LD panel (subset of [0, _snpNum), maintained by ld_esi_man)
  std::map<std::string, int> _snp_name_map;  // SNP rs ID -> SNP index

  std::vector<float> _esi_freq;  // frequency of A1 per SNP (used by the allele-frequency QC)

  std::vector<std::uint64_t> _cols;  // .bld cols array: cumulative LD value counts per SNP (random access into _val)
  std::vector<float> _val;           // all LD values in .bld order (r or r^2, per the header indicator)

  std::uint64_t _snpNum;  // SNP count (from the .bld header)
  std::uint64_t _valNum;  // total number of LD values (from the .bld header)
};

void filter_snp_maf(bInfo* bdata, double maf);
void calcu_mu(bInfo* bdata, bool ssq_flag = false);
void ld_report(char* outFileName, char* bFileName, char* indilstName, char* indilst2remove, char* snplstName,
               char* snplst2exclde, int chr, char* rs, double maf, bool ldr, bool ldr2, int ldWind);
void lookup(char* outFileName, char* bldFileName, char* snplstName, char* snplst2exclde, int chr, char* snprs,
            char* snprs2exclde, char* fromsnprs, char* tosnprs, int snpWind, bool snpWindflg, int fromsnpkb,
            int tosnpkb, int ld_wind);
void read_ld_esifile(ldInfo* ldinfo, char* esiFileName);
void ld_esi_man(ldInfo* ldinfo, char* snplstName, char* snplst2exclde, int chr, char* snprs, char* fromsnprs,
                char* tosnprs, int snpWind, bool snpwindFlag, int fromsnpkb, int tosnpkb, char* snprs2exclde);
void fetch_ld_by_id(ldInfo* ldinfo, FILE* ldfprt, std::vector<std::uint32_t>& curId, int sid, std::vector<float>& ld);
void ld_calc_o2m(VectorXd& ld_v, long target, const MatrixXd& X, bool centered = false);
}  // namespace SMRDATA
