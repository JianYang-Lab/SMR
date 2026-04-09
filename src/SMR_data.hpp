//
//  SMR_data.h
//  SRM_CPP
//
//  Created by Futao Zhang on 29/06/15.
//  Copyright (c) 2015 Futao Zhang. All rights reserved.
//

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "bfile.hpp"
#include "mmap_read.hpp"

#define MAX_LEN_PRBNAME 40
#define MAX_LEN_GENENAME 40
#define MAX_LEN_PATH 512

extern int thread_num;
extern bool mute;
extern int xh;  // using when testing
extern bool forcefrqck;
extern FILE* techeQTLfile;
extern char* outFileName;

namespace SMRDATA {

struct gwasData {
  long snpNum;
  std::vector<std::string> snpName;
  std::vector<int> snpBp;
  std::vector<std::string> allele_1;
  std::vector<std::string> allele_2;
  std::vector<double> freq;
  std::vector<double> byz;
  std::vector<double> seyz;
  std::vector<double> pvalue;
  std::vector<std::uint32_t> splSize;
  std::vector<int> _include;
  std::unordered_map<std::string, int> _snp_name_map;

  void reset() {
    snpNum = 0;
    _include.clear();
    _snp_name_map.clear();
    snpName.clear();
    snpBp.clear();
    allele_1.clear();
    allele_2.clear();
    freq.clear();
    byz.clear();
    seyz.clear();
    pvalue.clear();
    splSize.clear();
  }

  bool containsSNP(const std::string& snp_name) const { return _snp_name_map.find(snp_name) != _snp_name_map.end(); }
};

struct eqtlInfo {
  std::vector<int> _esi_chr;
  std::vector<std::string> _esi_rs;
  std::vector<int> _esi_gd;
  std::vector<int> _esi_bp;
  std::vector<std::string> _esi_allele1;
  std::vector<std::string> _esi_allele2;
  std::vector<int> _esi_include;  // snp indices (lineNum in *.esi file), initialized in the readesi
  // std::map<std::string,int> _snp_name_map;
  std::unordered_map<std::string, int> _snp_name_map;

  std::vector<float> _esi_freq;

  std::vector<int> _epi_chr;
  std::vector<std::string> _epi_prbID;
  std::vector<int> _epi_gd;
  std::vector<int> _epi_bp;
  std::vector<std::string> _epi_gene;
  std::vector<char> _epi_orien;
  // initialized in the readepi, probe_idx
  std::vector<int> _include;
  std::unordered_map<std::string, int> _probe_name_map;
  std::vector<double> _epi_var;
  /* if no probe sequence region input, its size should be 0.
     for the probe not for probe sequence file, the value should be
     std::set as -9, no technical eQTL would be removed from this probe.
   */
  std::vector<int> _epi_start;
  std::vector<int> _epi_end;

  // for sparse
  std::vector<std::uint64_t> _cols;
  std::vector<std::uint32_t> _rowid;  // index of snp in _esi_include
  std::vector<float> _val;
  // for dense
  std::vector<std::vector<float>> _bxz;  // first dimension is probe, second is snp
  std::vector<std::vector<float>> _sexz;

  std::uint64_t _probNum;
  std::uint64_t _snpNum;
  std::uint64_t _valNum;

  void reset_esi() {
    _esi_chr.clear();
    _esi_rs.clear();
    _esi_gd.clear();
    _esi_bp.clear();
    _esi_allele1.clear();
    _esi_allele2.clear();
    _esi_include.clear();
    _snp_name_map.clear();
    _esi_freq.clear();
  }

  void reset_epi() {
    _epi_chr.clear();
    _epi_prbID.clear();
    _epi_gd.clear();
    _epi_bp.clear();
    _epi_gene.clear();
    _epi_orien.clear();
    _include.clear();
    _probe_name_map.clear();
  }

  void reset_besd_sparse() {
    _cols.clear();
    _rowid.clear();
    _val.clear();
    _valNum = 0;
  }

  void reset_besd_dense() {
    _bxz.clear();
    _sexz.clear();
  }

  bool containsSNP(const std::string& snp_name) const { return _snp_name_map.find(snp_name) != _snp_name_map.end(); }
};

struct probeinfolst {
  char* probeId;
  char* genename;
  char* esdpath;
  char* bfilepath;
  int probechr;
  int gd;
  int bp;
  char orien;
  int start;
  int end;
};

struct snpinfolst {
  char* snprs;
  char* a1;
  char* a2;
  int snpchr;
  int gd;
  int bp;
  float beta;
  float se;
  float freq;
  float estn;
};

struct probeinfolst2 {
  std::vector<std::string> besdpath;
  char* probeId;
  char* genename;
  int probechr;
  int gd;
  int bp;
  char orien;
  std::uint64_t vnum;
  std::uint32_t* rowid;
  float* beta_se;
  snpinfolst* sinfo;
};

struct info4trans {
  char* snprs;
  char* a1;
  char* a2;
  int snpchr;
  int gd;  // to store curID
  int bp;
  float beta;
  float se;
  float freq;
  float byz;
  float seyz;
  float pyz;
};

struct SMRWK {
  int cur_chr;
  int cur_prbidx;
  std::vector<double> bxz, sexz, freq, byz, seyz;
  std::vector<double> pyz, zxz;
  std::vector<std::uint32_t> curId;
  std::vector<int> bpsnp, snpchrom;
  std::vector<std::string> rs, allele1, allele2;

  void init(int chr, int prbidx, size_t data_size) {
    cur_chr = chr;
    cur_prbidx = prbidx;
    bxz.reserve(data_size);
    sexz.reserve(data_size);
    freq.reserve(data_size);
    byz.reserve(data_size);
    seyz.reserve(data_size);
    pyz.reserve(data_size);
    zxz.reserve(data_size);
    curId.reserve(data_size);
    bpsnp.reserve(data_size);
    snpchrom.reserve(data_size);
    rs.reserve(data_size);
    allele1.reserve(data_size);
    allele2.reserve(data_size);
  }
};

struct SMRRLT {
  std::string ProbeID;
  int ProbeChr;
  std::string Gene;
  int Probe_bp;
  std::string SNP;
  int SNP_Chr;
  int SNP_bp;
  std::string A1;
  std::string A2;
  float Freq;
  float b_GWAS;
  float se_GWAS;
  double p_GWAS;
  float b_eQTL;
  float se_eQTL;
  double p_eQTL;
  float b_SMR;
  float se_SMR;
  double p_SMR;
  double p_HET;
  int nsnp;
  char Orien;
  double p_SSMR;
};

struct smr_probeinfo {
  int* ptr;
  char* probeId;
  char* genename;
  char* esdpath;
  char* bfilepath;
  int probechr;
  int gd;
  int bp;
  char orien;
};

struct smr_snpinfo {
  int* rstr;
  bool* revs;
  char* snprs;
  char* a1;
  char* a2;
  int snpchr;
  int gd;
  int bp;
  float beta;
  float se;
  float freq;
  float estn;
};

void read_bimfile(bInfo* bdata, std::string bimfile);
void read_famfile(bInfo* bdata, std::string famfile);
void read_bedfile(bInfo* bdata, std::string bedfile);
void read_gwas_data(gwasData* gdata, char* gwasFileName, bool enableGwasComments = false);
void read_esifile(eqtlInfo* eqtlinfo, std::string esifile, bool prtscr = true);
void read_esifile_by_chr(eqtlInfo* eqtlinfo, std::string esifile, int snpchr, bool prtscr = true);

void read_epifile(eqtlInfo* eqtlinfo, std::string epifile, bool prtscr = true);
void read_besdfile(eqtlInfo* eqtlinfo, std::string besdfile, bool prtscr = true);
void read_besdfile_mmap(eqtlInfo* eqtlinfo, MappedFile mapped, bool prtscr = true);

bool has_prefix(const std::string& str, const std::string& prefix);
bool has_suffix(const std::string& str, const std::string& suffix);
void get_square_idxes(std::vector<int>& sn_ids, VectorXd& zsxz, double threshold);
void est_cov_bxy(MatrixXd& covbxy, VectorXd& _zsxz, VectorXf& _bxy, VectorXd& _seyz, VectorXd& _bxz,
                 MatrixXd& _LD_heidi);
double bxy_hetero3(VectorXd& _byz, VectorXd& _bxz, VectorXd& _seyz, VectorXd& _sexz, VectorXd& _zsxz,
                   MatrixXd& _LD_heidi, long* snp_num);
float bxy_mltheter_so(VectorXd& _byz, VectorXd& _bxz, VectorXd& _seyz, VectorXd& _sexz, VectorXd& _zsxz,
                      MatrixXd& _LD_heidi, long* snp_num, double theta);
void allele_check(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata);
void allele_check(ldInfo* ldinfo, gwasData* gdata, eqtlInfo* esdata, double maf, double& freqthresh,
                  double& percenthresh);

void allele_check_opt(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata);
void update_geIndx(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata);
void allele_check(gwasData* gdata, eqtlInfo* esdata);
void update_gwas(gwasData* gdata);

bool make_XMat(bInfo* bdata, MatrixXd& X);
void makex_eigenVector(bInfo* bdata, int j, VectorXd& x, bool resize, bool minus_2p);
// inline functions
template <typename ElemType>
void makex(bInfo* bdata, int j, std::vector<ElemType>& x, bool minus_2p = false) {
  int i = 0;
  x.resize(bdata->_keep.size());
  for (i = 0; i < bdata->_keep.size(); i++) {
    if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
      if (bdata->_allele1[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]])
        x[i] =
            (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
      else
        x[i] = 2.0 - (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] +
                      bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
    } else x[i] = bdata->_mu[bdata->_include[j]];
    if (minus_2p) x[i] -= bdata->_mu[bdata->_include[j]];
  }
}

void makeptrx(bInfo* bdata, int bsnpid, int cursnpid, float* x, bool minus_2p = false);
void keep_indi(bInfo* bdata, std::string indi_list_file);
void extract_snp(bInfo* bdata, std::string snplistfile);
void extract_snp(bInfo* bdata, int chr);

void extract_prob(eqtlInfo* eqtlinfo, std::string problstName);
void extract_eqtl_snp(eqtlInfo* eqtlinfo, std::string snplstName);
void exclude_eqtl_snp(eqtlInfo* eqtlinfo, std::string snplstName);
void exclude_prob(eqtlInfo* eqtlinfo, std::string problstName);

void free_gwas_data(gwasData* gdata);

int file_read_check(std::ifstream* in_file, const char* filename);
void init_smr_wk(SMRWK* smrwk);
long fill_smr_wk(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata, SMRWK* smrwk, const char* refSNP, int lowerbp,
                 int upperbp, bool heidioffFlag);
long fill_smr_wk(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata, SMRWK* smrwk, const char* refSNP, int cis_itvl,
                 bool heidioffFlag);
long fill_smr_wk(ldInfo* ldinfo, gwasData* gdata, eqtlInfo* esdata, SMRWK* smrwk, const char* refSNP, int cis_itvl,
                 bool heidioffFlag);

double heidi_test(bInfo* bdata, SMRWK* smrwk, long maxid, double ld_top, double threshold, int m_hetero, long& nsnp);
double heidi_test_new(bInfo* bdata, SMRWK* smrwk, double ldr2_top, double threshold, int m_hetero, long& nsnp,
                      double ld_min, int opt_hetero, bool sampleoverlap, double theta);
double heidi_test_new(ldInfo* ldinfo, FILE* ldfptr, SMRWK* smrwk, double ldr2_top, double threshold, int m_hetero,
                      long& nsnp, double ld_min, int opt_hetero, bool sampleoverlap, double theta);

void update_snidx(SMRWK* smrwk, std::vector<int>& sn_ids, int max_snp_slct, std::string forwhat);
void extract_smrwk(SMRWK* smrwk, std::vector<int>& sn_ids, SMRWK* smrwk2);
void rm_cor_sbat(MatrixXd& R, double R_cutoff, int m, std::vector<int>& rm_ID1);
void update_smrwk_x(SMRWK* smrwk, std::vector<int>& sn_ids, MatrixXd& X);
void smr_heidi_func(std::vector<SMRRLT>& smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, eqtlInfo* esdata,
                    int cis_itvl, bool heidioffFlag, double heidiskipthresh, const char* refSNP, double p_hetero,
                    double ld_top, int m_hetero, double p_smr, double threshpsmrest, bool new_het_mtd, bool opt,
                    double ld_min, int opt_hetero, bool sampleoverlap, double pmecs, int minCor,
                    std::unordered_map<std::string, std::string>& prb_snp, bool targetLstFlg);
void smr(char* outFileName, char* bFileName, char* bldFileName, char* gwasFileName, char* eqtlFileName, double maf,
         char* indilstName, char* snplstName, char* problstName, bool bFlag, double p_hetero, double ld_top,
         int m_hetero, int opt_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde, double p_smr,
         char* refSNP, bool heidioffFlag, double heidiskipthresh, int cis_itvl, char* genelistName, int chr, int prbchr,
         const char* prbname, char* fromprbname, char* toprbname, int prbWind, int fromprbkb, int toprbkb,
         bool prbwindFlag, char* genename, int snpchr, char* snprs, char* fromsnprs, char* tosnprs, int snpWind,
         int fromsnpkb, int tosnpkb, bool snpwindFlag, bool cis_flag, double threshpsmrest, bool new_het_mth, bool opt,
         char* prbseqregion, double ld_min, bool sampleoverlap, double pmecs, int minCor, char* targetsnpproblstName,
         char* snpproblstName, double afthresh, double percenthresh);

double heidi_test_ref_new(ldInfo* ldinfo, FILE* ldfptr, SMRWK* smrwk, double ldr2_top, double threshold, int m_hetero,
                          long& nsnp, int refid, double ld_min, int opt_hetero, bool sampleoverlap, double theta);

void smr_trans(char* outFileName, char* bFileName, char* gwasFileName, char* eqtlFileName, double maf,
               char* indilstName, char* snplstName, char* problstName, bool bFlag, double p_hetero, double ld_top,
               int m_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde, double p_trans,
               char* refSNP, bool heidioffFlag, int cis_itvl, int trans_itvl, char* genelistName, int chr, int prbchr,
               const char* prbname, char* fromprbname, char* toprbname, int prbWind, int fromprbkb, int toprbkb,
               bool prbwindFlag, char* genename, int snpchr, char* snprs, char* fromsnprs, char* tosnprs, int snpWind,
               int fromsnpkb, int tosnpkb, bool snpwindFlag, bool cis_flag, double threshpsmrest, bool new_het_mth,
               double p_smr, double ld_min);
void smr_trans_region(char* outFileName, char* bFileName, char* gwasFileName, char* eqtlFileName, double maf,
                      char* indilstName, char* snplstName, char* problstName, bool bFlag, double p_hetero,
                      double ld_top, int m_hetero, int opt_hetero, char* indilst2remove, char* snplst2exclde,
                      char* problst2exclde, double p_trans, char* refSNP, bool heidioffFlag, int cis_itvl,
                      int trans_itvl, char* genelistName, int chr, int prbchr, const char* prbname, char* fromprbname,
                      char* toprbname, int prbWind, int fromprbkb, int toprbkb, bool prbwindFlag, char* genename,
                      int snpchr, char* snprs, char* fromsnprs, char* tosnprs, int snpWind, int fromsnpkb, int tosnpkb,
                      bool snpwindFlag, bool cis_flag, double threshpsmrest, bool new_het_mth, double p_smr, bool opt,
                      double ld_min, double afthresh, double percenthresh);

void make_full_besd(char* outFileName, char* eqtlFileName, char* snplstName, char* problstName, bool bFlag,
                    bool make_besd_flag, char* snplst2exclde, char* problst2exclde, int cis_itvl, char* genelistName,
                    int chr, int prbchr, char* prbname, char* fromprbname, char* toprbname, int prbWind, int fromprbkb,
                    int toprbkb, bool prbwindFlag, char* genename, int snpchr, char* snprs, char* fromsnprs,
                    char* tosnprs, int snpWind, int fromsnpkb, int tosnpkb, bool snpwindFlag, bool cis_flag, int addn);

void smr_e2e(char* outFileName, char* bFileName, char* eqtlFileName, char* eqtlFileName2, double maf, char* indilstName,
             char* snplstName, char* problstName, char* eproblstName, char* mproblstName, bool bFlag, double p_hetero,
             double ld_prune, int m_hetero, int opt_hetero, char* indilst2remove, char* snplst2exclde,
             char* problst2exclde, char* eproblst2exclde, char* mproblst2exclde, double p_smr, char* refSNP,
             bool heidioffFlag, double heidiskipthresh, int cis_itvl, char* traitlstName, int op_wind, char* oprobe,
             char* eprobe, char* oprobe2rm, char* eprobe2rm, double prepsmrthres, bool new_het_mth, bool opt,
             double ld_min, bool cis2all, bool sampleoverlap, double pmecs, int minCor, bool ssmrflg, int expanWind,
             double ld_top_multi, char* targetsnpproblstName, char* snpproblstName, double afthresh,
             double percenthresh);

void ssmr_heidi_func(std::vector<SMRRLT>& smrrlts, char* outFileName, bInfo* bdata, gwasData* gdata, eqtlInfo* esdata,
                     int cis_itvl, bool heidioffFlag, double heidiskipthresh, const char* refSNP, double p_hetero,
                     double ld_top, int m_hetero, double p_smr, double threshpsmrest, double ld_min, int opt_hetero,
                     int expanWind, bool sampleoverlap, double pmecs, int minCor, double ld_top_multi);
void meta(char* outFileName, char* eqtlFileName, char* eqtlFileName2);
int comp_esi(const void* a, const void* b);
int comp(const void* a, const void* b);
int comp2(const void* a, const void* b);
int comp_estn(const void* a, const void* b);
int comp_i4tran(const void* a, const void* b);

void read_smaslist(std::vector<std::string>& smasNames, std::string eqtlsmaslstName);
void remove_indi(bInfo* bdata, std::string indi_list_file);
void extract_snp(bInfo* bdata, std::string snplistfile);
void exclude_snp(bInfo* bdata, std::string snplistfile);
void allele_check(bInfo* bdata, eqtlInfo* etrait, eqtlInfo* esdata);
void update_geIndx(bInfo* bdata, eqtlInfo* etrait, eqtlInfo* esdata);
void progress_print(float progress);
void make_XMat(bInfo* bdata, std::vector<std::uint32_t>& snpids, MatrixXd& X, bool minus_2p = false);
void get_square_ldpruning_idxes(std::vector<int>& sn_ids, VectorXd& zsxz, double threshold, VectorXd& ld_v, long maxid,
                                double ld_top);
void cor_calc(MatrixXd& LD, MatrixXd& X);
void cor_calc(MatrixXd& LD, ldInfo* ldinfo, FILE* ldfprt, const std::vector<std::uint32_t>& curId, int indicator);

void extract_prob_by_gene(eqtlInfo* eqtlinfo, std::string genelistName);
void update_freq(char* eqtlFileName, char* frqfile);
void write_besd(std::string outFileName, eqtlInfo* eqtlinfo);
void extract_region_bp(bInfo* bdata, int chr, int fromkb, int tokb);
void extract_eqtl_single_snp(eqtlInfo* eqtlinfo, std::string snprs);
void extract_eqtl_snp(eqtlInfo* eqtlinfo, std::string fromsnprs, std::string tosnprs);
void extract_eqtl_snp(eqtlInfo* eqtlinfo, int chr, int fromsnpkb, int tosnpkb);
void extract_prob(eqtlInfo* eqtlinfo, std::string prbname, int prbWind);
void extract_eqtl_single_probe(eqtlInfo* eqtlinfo, std::string prbname, bool prtscr = true);
void extract_eqtl_prob(eqtlInfo* eqtlinfo, std::string fromprbname, std::string toprbname);
void extract_eqtl_prob(eqtlInfo* eqtlinfo, int chr, int fromprbkb, int toprbkb);
void extract_prob_by_single_gene(eqtlInfo* eqtlinfo, std::string genename);
void extract_epi_by_chr(eqtlInfo* eqtlinfo, int prbchr);
void extract_eqtl_by_chr(eqtlInfo* eqtlinfo, int snpchr);
void extract_eqtl_snp(eqtlInfo* eqtlinfo, std::string snporprb, int Wind, std::string msg);
void epi_man(eqtlInfo* eqtlinfo, char* problstName, char* genelistName, int chr, int prbchr, const char* prbname,
             char* fromprbname, char* toprbname, int prbWind, int fromprbkb, int toprbkb, bool prbwindFlag,
             char* genename);
void esi_man(eqtlInfo* eqtlinfo, char* snplstName, int chr, int snpchr, char* snprs, char* fromsnprs, char* tosnprs,
             int snpWind, int fromsnpkb, int tosnpkb, bool snpwindFlag, bool cis_flag, int cis_itvl,
             const char* prbname);
void exclude_eqtl_single_probe(eqtlInfo* eqtlinfo, std::string prbname);
void slct_sparse_per_prb(std::vector<int>& slct_idx, probeinfolst* prbifo, std::vector<snpinfolst>& snpinfo,
                         long cis_itvl, long trans_itvl, double transThres, double restThres, FILE* logfile,
                         bool extract_cis_only, bool techHit = false);
void read_gene_anno(char* geneAnnoName, std::vector<int>& chr, std::vector<std::string>& genename,
                    std::vector<int>& start, std::vector<int>& end);
void read_gene_anno_strand(char* geneAnnoName, std::vector<int>& chr, std::vector<std::string>& genename,
                           std::vector<int>& start, std::vector<int>& end, std::vector<std::string>& strand);
void rm_unmatched_snp(gwasData* gdata, eqtlInfo* esdata);
void rm_unmatched_snp(eqtlInfo* etrait, eqtlInfo* esdata);
void filter_snp_null(eqtlInfo* eqtlinfo);
int max_zsmr_id(SMRWK* smrwk, double p_smr);
double heidi_test_ref_new(bInfo* bdata, SMRWK* smrwk, double ldr2_top, double threshold, int m_hetero, long& nsnp,
                          int refid, double ld_min, int opt_hetero, bool sampleoverlap, double theta);
double heidi_test_ref_new(ldInfo* ldinfo, FILE* ldfptr, SMRWK* smrwk, double ldr2_top, double threshold, int m_hetero,
                          long& nsnp, int refid, double ld_min, int opt_hetero, bool sampleoverlap, double theta);
void free_snplist(std::vector<snpinfolst>& a);
void free_probelist(std::vector<probeinfolst>& a);
void free_snplist(std::vector<info4trans>& a);
void free_probelist(std::vector<probeinfolst2>& a);
void read_epistartend(eqtlInfo* eqtlinfo, char* prbseqregion);
int shown(std::string besdfile);
double freq_check(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata, double& freqthresh, double& percenthresh);
void slct_trans_per_prb(std::vector<int>& slct_idx, std::vector<int>& regionChr, std::vector<long>& snpNumPerRegion,
                        std::vector<long>& leftbound, std::vector<long>& rightbound, probeinfolst* prbifo,
                        std::vector<info4trans>& snpinfo, long cis_itvl, long trans_itvl, double transThres);
void extract_targets(eqtlInfo* eqtlinfo, std::string snpprblistfile,
                     std::unordered_map<std::string, std::string>& prb_snp);
}  // namespace SMRDATA
