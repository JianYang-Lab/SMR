//  SMR_data_p2.h
//  SMR_CPP
//
//  Created by Futao Zhang on 10/06/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//  Modified by fanghl 20210412

#ifndef __SMR_CPP__SMR_data_p2__
#define __SMR_CPP__SMR_data_p2__

#include "SMR_data.hpp"

namespace SMRDATA
{
    void make_besd(char * outFileName, char * syllabusName, bool gctaflag, \
        bool plinkflag,bool gemmaflag, bool merlinflag, bool boltflag, \
        bool save_dense_flag, int cis_itvl, int trans_itvl, \
        float transThres, float restThres, bool samegeno, int addn);
    void make_besd_byQfile(char * qfile, char * outFileName, \
        bool save_dense_flag, int cis_itvl, int trans_itvl, \
        float transThres, float restThres, int addn);
    void make_besd_fmat(char * fmatfileName, char * outFileName, \
        bool mateqtlflag, bool fastnflag, bool fastpflag, \
        bool qtltoolsnflag, bool qtltoolspflag, int addn);
    void make_qfile(char* outFileName, char* eFileName, char* mafFileName);
    void selective_cpy(eqtlInfo* etrait, eqtlInfo* esdata);
    void write_smr_esi(string outFileName, eqtlInfo* eqtlinfo);
    void write_smr_epi(string outFileName, eqtlInfo* eqtlinfo);
}
#endif /* defined(__SMR_CPP__SMR_data_p2__) */
