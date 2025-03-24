//
//  bfile.hpp
//  SMR_CPP
//
//  Created by Futao Zhang on 5/07/2018.
//  Copyright © 2018 Futao Zhang. All rights reserved.
//

#ifndef bfile_hpp
#define bfile_hpp
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
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
#include <bitset>
#include "CommFunc.hpp"
#include "StatFunc.hpp"
#include "StrFunc.hpp"

using namespace std;
using namespace StatFunc;
using namespace StrFunc;
using namespace CommFunc;

namespace SMRDATA
{
    //implementation of SoA(Struct of Array)
    typedef struct{
        // bim file
        int _autosome_num;
        vector<int> _chr;
        vector<string> _snp_name;
        map<string, int> _snp_name_map;
        vector<double> _genet_dst;
        vector<int> _bp;
        vector<string> _allele1;
        vector<string> _allele2;
        vector<string> _ref_A; // reference allele
        vector<string> _other_A; // the other allele
        int _snp_num;
        vector<double> _rc_rate;
        vector<int> _include; // initialized in the read_bimfile()
        VectorXd _maf;
        
        // fam file
        vector<string> _fid;
        vector<string> _pid;
        map<string, int> _id_map;
        vector<string> _fa_id;
        vector<string> _mo_id;
        vector<int> _sex;
        vector<double> _pheno;
        int _indi_num;
        vector<int> _keep; // initialized in the read_famfile()
        MatrixXd _varcmp_Py; // BLUP solution to the total genetic effects of individuals
        
        // bed file
        vector< vector<bool> > _snp_1;
        vector< vector<bool> > _snp_2;
        
        // imputed data
        bool _dosage_flag;
        vector< vector<float> > _geno_dose;
        vector<double> _impRsq;
        
        // genotypes
        MatrixXf _geno;
        
        
        vector<double> _mu;
    } bInfo;
    
    typedef struct{
        vector<int> _esi_chr;
        vector<string> _esi_rs;
        vector<int> _esi_gd;
        vector<int> _esi_bp;
        vector<string> _esi_allele1;
        vector<string> _esi_allele2;
        vector<int> _esi_include;
        map<string,int> _snp_name_map;
        vector<float> _esi_freq;
        
        vector<uint64_t> _cols;
        vector<float> _val;
        
        uint64_t _snpNum;
        uint64_t _valNum;
        
    } ldInfo;
    void filter_snp_maf(bInfo* bdata,double maf);
    void calcu_mu(bInfo* bdata, bool ssq_flag=false);
    void ld_report(char* outFileName, char* bFileName,char* indilstName, char* indilst2remove,char* snplstName, char* snplst2exclde,int chr, char* rs, double maf, bool ldr, bool ldr2, int ldWind);
    void calcu_mean_rsq(char* outFileName, char* bFileName,char* indilstName, char* indilst2remove,char* snplstName, char* snplst2exclde,int chr,double maf, bool ldr, bool ldr2, int ldWind, double rsq_cutoff);
    void lookup(char* outFileName, char* bldFileName, char* snplstName, char* snplst2exclde,int chr,char* snprs, char* snprs2exclde, char* fromsnprs, char* tosnprs,int snpWind, bool snpWindflg, int fromsnpkb, int tosnpkb, int ld_wind);
    void read_ld_esifile(ldInfo* ldinfo, char* esiFileName);
    void ld_esi_man(ldInfo* ldinfo,char* snplstName, char* snplst2exclde, int chr,char* snprs, char* fromsnprs, char* tosnprs,int snpWind,bool snpwindFlag, int fromsnpkb, int tosnpkb,char* snprs2exclde);
    void fetch_ld_by_id(ldInfo* ldinfo,FILE* ldfprt, int sid, vector<float> &ld);
    void fetch_ld_by_id(ldInfo* ldinfo,FILE* ldfprt, vector<uint32_t> &curId, int sid, vector<float> &ld);
    void fetch_ld_by_snps(ldInfo* ldinfo,FILE* ldfprt, string rs, vector<float> &ld);
    void ld_calc_o2m(VectorXd &ld_v,long targetid, MatrixXd &X, bool centered=false);
}

#endif /* bfile_hpp */
