//  SMR.h
//  SMR_CPP
//
//  Created by Futao Zhang on 1/09/2015.
//  Copyright (c) 2015 Futao Zhang. All rights reserved.
//  Modified by fanghl in 20210512

#ifndef SMR_CPP_SMR_h
#define SMR_CPP_SMR_h

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "SMR_data.hpp"

void option(int option_num, char* option_str[]);

static inline bool not_in_flags(std::vector<std::string>& flags, std::string str) {
  return std::find(flags.begin(), flags.end(), str) == flags.end();
}

static inline void FLAGS_VALID_CK(int option_num, char* option_str[]) {
  const char* flgs[] = {"--bfile",
                        "--gwas-summary",
                        "--beqtl-summary",
                        "--beqtl-list",
                        "--beqtl-list-index",
                        "--maf",
                        "--keep",
                        "--remove",
                        "--extract-snp",
                        "--exclude-snp",
                        "--extract-probe",
                        "--extract-outcome-probe",
                        "--extract-exposure-probe",
                        "--exclude-probe",
                        "--exclude-outcome-probe",
                        "--exclude-exposure-probe",
                        "--eqtl-summary",
                        "--ld-upper-limit",
                        "--peqtl-heidi",
                        "--heidi-min-m",
                        "--make-besd",
                        "--out",
                        "--peqtl-smr",
                        "--smr",
                        "--cis-wind",
                        "--peqtl-trans",
                        "--peqtl-cis",
                        "--peqtl-other",
                        "--efile",
                        "--query",
                        "--heidi-off",
                        "--heidi-skip-thresh",
                        "--target-snp",
                        "--extract-trait",
                        "--thread-num",
                        "--besd-flist",
                        "--trans-wind",
                        "--plot",
                        "--eqtl-flist",
                        "--smr-format",
                        "--merlin-fastassoc-format",
                        "--plink-qassoc-format",
                        "--gemma-format",
                        "--update-freq",
                        "--genes",
                        "--smr-wind",
                        "--smr-file",
                        "--recode",
                        "--probe-var",
                        "--chr",
                        "--probe-chr",
                        "--snp-chr",
                        "--snp",
                        "--from-snp",
                        "--to-snp",
                        "--probe",
                        "--from-probe",
                        "--to-probe",
                        "--snp-wind",
                        "--probe-wind",
                        "--gene",
                        "--from-snp-kb",
                        "--to-snp-kb",
                        "--from-probe-kb",
                        "--to-probe-kb",
                        "--extract-single-exposure-probe",
                        "--extract-single-outcome-probe",
                        "--exclude-single-exposure-probe",
                        "--exclude-single-outcome-probe",
                        "--gene-list",
                        "--set-list",
                        "--set-wind",
                        "--smr-multi",
                        "--qfile",
                        "--make-besd-dense",
                        "--bolt-assoc-format",
                        "--geno-uni",
                        "--psmr",
                        "--beqtl-qc",
                        "--z-thresh",
                        "--heidi-mtd",
                        "--phet",
                        "--descriptive-trans",
                        "--descriptive-cis",
                        "--meta",
                        "--trans",
                        "--extract-cis",
                        "--rm-technical",
                        "--p-technical",
                        "--ld-lower-limit",
                        "--max_num_ld",
                        "--heidi-max-m",
                        "--extract-snp-p",
                        "--exclude-snp-p",
                        "--matrix-eqtl-format",
                        "--fastqtl-nominal-format",
                        "--fastqtl-permu-format",
                        "--add-n",
                        "--show-n",
                        "--update-epi",
                        "--update-esi",
                        "--cis-to-all",
                        "--mecs",
                        "--pmecs",
                        "--mmecs",
                        "--sample-overlap",
                        "--ld-multi-snp",
                        "--extract-target-snp-probe",
                        "--extract-snp-probe",
                        "--disable-freq-ck",
                        "--diff-freq",
                        "--diff-freq-prop",
                        "--r",
                        "--r2",
                        "--ld-wind",
                        "--bld",
                        "--make-bld",
                        "--snp-rm",
                        "--qtltools-nominal-format",
                        "--qtltools-permu-format",
                        "--nmecs",
                        "--outcome-wind",
                        "--enable-gwas-comments"};

  // define varible flags, using start pointer and end pointer.
  std::vector<std::string> flags(flgs, flgs + sizeof(flgs) / sizeof(flgs[0]));
  if (option_num < 3) {
    std::cout << "flags include:" << std::endl;
    int cur_mark = 0;
    for (int i = 0; i < flags.size(); i++) {
      int tmp = i >> 2;
      if (tmp > cur_mark) {
        std::cout << std::endl;
        cur_mark = tmp;
      }
      std::cout << flags[i] << ",";
    }
    std::cout << std::endl;
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < option_num; i++) {
    if (SMRDATA::has_prefix(option_str[i], "--"))
      if (not_in_flags(flags, option_str[i])) {
        fprintf(stderr, "%s: Invalid option\n", option_str[i]);
        exit(EXIT_FAILURE);
      }
  }
}

static inline void FLAG_VALID_CK(std::string str, char* flag) {
  if (flag == NULL || SMRDATA::has_prefix(flag, "--")) {
    fprintf(stderr, "Please verify the flag %s!: \n", str.c_str());
    exit(EXIT_FAILURE);
  }
}
#endif
