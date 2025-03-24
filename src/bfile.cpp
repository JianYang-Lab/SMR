//
//  bfile.cpp
//  SMR_CPP
//
//  Created by Futao Zhang on 5/07/2018.
//  Copyright © 2018 Futao Zhang. All rights reserved.
//

#include "bfile.hpp"

namespace SMRDATA
{
    void read_famfile(bInfo* bdata, string famfile) {
        bdata->_autosome_num = 22;
        ifstream Fam(famfile.c_str());
        if (!Fam) throw ("Error: can not open the file [" + famfile + "] to read.");
        cout << "Reading PLINK FAM file from [" + famfile + "]." << endl;
        
        int i = 0;
        string str_buf;
        bdata->_fid.clear();
        bdata->_pid.clear();
        bdata->_fa_id.clear();
        bdata->_mo_id.clear();
        bdata->_sex.clear();
        bdata->_pheno.clear();
        while (Fam) {
            Fam >> str_buf;
            if (Fam.eof()) break;
            bdata->_fid.push_back(str_buf);
            Fam >> str_buf;
            bdata->_pid.push_back(str_buf);
            Fam >> str_buf;
            bdata->_fa_id.push_back(str_buf);
            Fam >> str_buf;
            bdata->_mo_id.push_back(str_buf);
            Fam >> str_buf;
            bdata->_sex.push_back(atoi(str_buf.c_str()));
            Fam >> str_buf;
            bdata->_pheno.push_back(atoi(str_buf.c_str()));
        }
        Fam.clear();
        Fam.close();
        bdata->_indi_num = bdata->_fid.size();
        cout << bdata->_indi_num << " individuals to be included from [" + famfile + "]." << endl;
        
        // Initialize _keep
        bdata->_keep.clear();
        bdata->_keep.resize(bdata->_indi_num);
        bdata->_id_map.clear();
        int size = 0;
        for (int i = 0; i < bdata->_indi_num; i++) {
            bdata->_keep[i] = i;
            bdata->_id_map.insert(pair<string, int>(bdata->_fid[i] + ":" + bdata->_pid[i], i));
            if (size == bdata->_id_map.size()) throw ("Error: Duplicate individual ID found: \"" + bdata->_fid[i] + "\t" + bdata->_pid[i] + "\".");
            size = bdata->_id_map.size();
        }
    }
    


    void update_bim(bInfo* bdata,vector<int> &rsnp) {
        int i = 0;
        
        //update bim information
        vector<int> chr_buf, bp_buf;
        vector<string> a1_buf, a2_buf, ref_A_buf, other_A_buf;
        vector<string> snp_name_buf;
        vector<double> genet_dst_buf, impRsq_buf;
        for (i = 0; i < bdata->_snp_num; i++) {
            if (!rsnp[i]) continue;
            chr_buf.push_back(bdata->_chr[i]);
            snp_name_buf.push_back(bdata->_snp_name[i]);
            genet_dst_buf.push_back(bdata->_genet_dst[i]);
            bp_buf.push_back(bdata->_bp[i]);
            a1_buf.push_back(bdata->_allele1[i]);
            a2_buf.push_back(bdata->_allele2[i]);
            ref_A_buf.push_back(bdata->_ref_A[i]);
            other_A_buf.push_back(bdata->_other_A[i]);
            if(bdata->_impRsq.size()>0) impRsq_buf.push_back(bdata->_impRsq[i]);
        }
        bdata->_chr.clear();
        bdata->_snp_name.clear();
        bdata->_genet_dst.clear();
        bdata->_bp.clear();
        bdata->_allele1.clear();
        bdata->_allele2.clear();
        bdata->_ref_A.clear();
        bdata->_other_A.clear();
        bdata->_impRsq.clear();
        bdata->_chr = chr_buf;
        bdata->_snp_name = snp_name_buf;
        bdata->_genet_dst = genet_dst_buf;
        bdata->_bp = bp_buf;
        bdata->_allele1 = a1_buf;
        bdata->_allele2 = a2_buf;
        bdata->_ref_A = ref_A_buf;
        bdata->_other_A = other_A_buf;
        bdata->_impRsq=impRsq_buf;
        bdata->_snp_num = bdata->_chr.size();
        bdata->_include.clear();
        bdata-> _include.resize(bdata->_snp_num);
        bdata->_snp_name_map.clear();
        
        for (i = 0; i < bdata->_snp_num; i++) {
            bdata->_include[i] = i;
            bdata->_snp_name_map.insert(pair<string, int>(bdata->_snp_name[i], i));
        }
    }
    

    void update_fam(bInfo* bdata,vector<int> &rindi) {
        //update fam information
        int i = 0;
        vector<string> fid_buf, pid_buf, fa_id_buf, mo_id_buf;
        vector<int> sex_buf;
        vector<double> pheno_buf;
        for (i = 0; i < bdata->_indi_num; i++) {
            if (!rindi[i]) continue;
            fid_buf.push_back(bdata->_fid[i]);
            pid_buf.push_back(bdata->_pid[i]);
            fa_id_buf.push_back(bdata->_fa_id[i]);
            mo_id_buf.push_back(bdata->_mo_id[i]);
            sex_buf.push_back(bdata->_sex[i]);
            pheno_buf.push_back(bdata->_pheno[i]);
        }
        bdata->_fid.clear();
        bdata->_pid.clear();
        bdata->_fa_id.clear();
        bdata->_mo_id.clear();
        bdata->_sex.clear();
        bdata->_pheno.clear();
        bdata->_fid = fid_buf;
        bdata->_pid = pid_buf;
        bdata->_fa_id = fa_id_buf;
        bdata->_mo_id = mo_id_buf;
        bdata->_sex = sex_buf;
        bdata->_pheno = pheno_buf;
        
        bdata->_indi_num = bdata->_fid.size();
        bdata->_keep.clear();
        bdata->_keep.resize(bdata->_indi_num);
        bdata->_id_map.clear();
        for (i = 0; i < bdata->_indi_num; i++) {
            bdata->_keep[i] = i;
            bdata->_id_map.insert(pair<string, int>(bdata->_fid[i] + ":" + bdata->_pid[i], i));
        }
    }
    

    void read_bimfile(bInfo* bdata,string bimfile) {
        // Read bim file: recombination rate is defined between SNP i and SNP i-1
        int ibuf = 0;
        string cbuf = "0";
        double dbuf = 0.0;
        string str_buf;
        ifstream Bim(bimfile.c_str());
        if (!Bim) throw ("Error: can not open the file [" + bimfile + "] to read.");
        cout << "Reading PLINK BIM file from [" + bimfile + "]." << endl;
        bdata->_chr.clear();
        bdata->_snp_name.clear();
        bdata->_genet_dst.clear();
        bdata->_bp.clear();
        bdata->_allele1.clear();
        bdata->_allele2.clear();
        while (Bim) {
            Bim >> ibuf;
            if (Bim.eof()) break;
            bdata->_chr.push_back(ibuf);
            Bim >> str_buf;
            bdata->_snp_name.push_back(str_buf);
            Bim >> dbuf;
            bdata->_genet_dst.push_back(dbuf);
            Bim >> ibuf;
            bdata->_bp.push_back(ibuf);
            Bim >> cbuf;
            StrFunc::to_upper(cbuf);
            bdata->_allele1.push_back(cbuf.c_str());
            Bim >> cbuf;
            StrFunc::to_upper(cbuf);
            bdata->_allele2.push_back(cbuf.c_str());
        }
        Bim.close();
        bdata->_snp_num = bdata->_chr.size();
        bdata->_ref_A = bdata->_allele1;
        bdata->_other_A = bdata->_allele2;
        cout << bdata->_snp_num << " SNPs to be included from [" + bimfile + "]." << endl;
        
        // Initialize _include
        bdata->_include.clear();
        bdata->_include.resize( bdata->_snp_num);
        bdata->_snp_name_map.clear();
        
        for (int i = 0; i <  bdata->_snp_num; i++) {
            
            bdata->_include[i] = i;
            if( bdata->_snp_name_map.find(bdata->_snp_name[i]) != bdata->_snp_name_map.end()){
                cout << "Warning: Duplicated SNP ID \"" + bdata->_snp_name[i] + "\" ";
                stringstream ss;
                ss << bdata->_snp_name[i] << "_" << i + 1;
                bdata->_snp_name[i] = ss.str();
                cout<<"has been changed to \"" + bdata->_snp_name[i] + "\".\n";
            }
            bdata->_snp_name_map.insert(pair<string, int>(bdata->_snp_name[i], i));
        }
    }


    // some code are adopted from PLINK with modifications
    void read_bedfile(bInfo* bdata, string bedfile)
    {
        int i = 0, j = 0, k = 0;
        
        // Flag for reading individuals and SNPs
        vector<int> rindi, rsnp;
        //get_rindi
        rindi.clear();
        rindi.resize(bdata->_indi_num);
        for (int i = 0; i < bdata->_indi_num; i++) {
            if (bdata->_id_map.find(bdata->_fid[i] + ":" + bdata->_pid[i]) != bdata->_id_map.end()) rindi[i] = 1;
            else rindi[i] = 0;
        }
        //get_rsnp
        rsnp.clear();
        rsnp.resize(bdata->_snp_num);
        for (int i = 0; i < bdata->_snp_num; i++) {
            if (bdata->_snp_name_map.find(bdata->_snp_name[i]) != bdata->_snp_name_map.end()) rsnp[i] = 1;
            else rsnp[i] = 0;
        }
        
        if (bdata->_include.size() == 0) throw ("Error: No SNP is retained for analysis.");
        if (bdata->_keep.size() == 0) throw ("Error: No individual is retained for analysis.");
        
        // Read bed file
        char ch[1];
        bitset<8> b;
        bdata->_snp_1.resize(bdata->_include.size());
        bdata->_snp_2.resize(bdata->_include.size());
        for (i = 0; i < bdata->_include.size(); i++) {
            bdata->_snp_1[i].reserve(bdata->_keep.size());
            bdata->_snp_2[i].reserve(bdata->_keep.size());
        }
        fstream BIT(bedfile.c_str(), ios::in | ios::binary);
        if (!BIT) throw ("Error: can not open the file [" + bedfile + "] to read.");
        cout << "Reading PLINK BED file from [" + bedfile + "] in SNP-major format ..." << endl;
        for (i = 0; i < 3; i++) BIT.read(ch, 1); // skip the first three bytes
        int snp_indx = 0, indi_indx = 0;
        for (j = 0, snp_indx = 0; j < bdata->_snp_num; j++) { // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
            if (!rsnp[j]) {
                for (i = 0; i < bdata->_indi_num; i += 4) BIT.read(ch, 1);
                continue;
            }
            for (i = 0, indi_indx = 0; i < bdata->_indi_num;) {
                BIT.read(ch, 1);
                if (!BIT) throw ("Error: problem with the BED file ... has the FAM/BIM file been changed?");
                b = ch[0];
                k = 0;
                while (k < 7 && i < bdata->_indi_num) { // change code: 11 for AA; 00 for BB;
                    if (!rindi[i]) k += 2;
                    else {
                        bdata->_snp_2[snp_indx][indi_indx] = (!b[k++]);
                        bdata->_snp_1[snp_indx][indi_indx] = (!b[k++]);
                        indi_indx++;
                    }
                    i++;
                }
            }
            if (snp_indx == bdata->_include.size()) break;
            snp_indx++;
        }
        BIT.clear();
        BIT.close();
        cout << "Genotype data for " << bdata->_keep.size() << " individuals and " << bdata->_include.size() << " SNPs to be included from [" + bedfile + "]." << endl;
        
        update_fam(bdata, rindi);
        update_bim(bdata, rsnp);
    }


    void keep_indi(bInfo* bdata,string indi_list_file)
    {
        vector<string> indi_list;
        read_indi_list(indi_list_file, indi_list);
        update_id_map_kp(indi_list, bdata->_id_map, bdata->_keep);
        cout<<bdata->_keep.size()<<" individuals are kept from ["+indi_list_file+"]."<<endl;
    }
    

    void remove_indi(bInfo* bdata, string indi_list_file) {
        vector<string> indi_list;
        read_indi_list(indi_list_file, indi_list);
        int prev_size = bdata->_keep.size();
        update_id_map_rm(indi_list, bdata->_id_map, bdata->_keep);
        cout << prev_size - bdata->_keep.size() << " individuals are removed from [" + indi_list_file + "] and there are " << bdata->_keep.size() << " individuals remaining." << endl;
    }


    void extract_region_bp(bInfo* bdata, int chr, int fromkb, int tokb)
    {
        int frombp=fromkb*1000;
        int tobp=tokb*1000;
        vector<string> snplist;
        for(int i = 0; i < bdata->_include.size(); i++){
            int j = bdata->_include[i];
            if(bdata->_chr[j] == chr && bdata->_bp[j]<=tobp && bdata->_bp[j]>=frombp) snplist.push_back(bdata->_snp_name[j]);
        }
        if(snplist.empty()) throw ("Error: on SNP found in this region.");
        update_id_map_kp(snplist, bdata->_snp_name_map, bdata->_include);
        cout << bdata->_include.size() << " SNPs are extracted from SNP BP: "<<fromkb<<"Kb"<<"to SNP BP: "<<tokb<<"Kb."<< endl;
    }


    void extract_snp(bInfo* bdata, int chr)
    {
        vector<string> snplist;
        for(int i = 0; i < bdata->_include.size(); i++){
            int j = bdata->_include[i];
            if(bdata->_chr[j] == chr) snplist.push_back(bdata->_snp_name[j]);
        }
        if(snplist.empty()) throw ("Error: on SNP found in this region.");
        update_id_map_kp(snplist, bdata->_snp_name_map, bdata->_include);
        cout << bdata->_include.size() << " SNPs are extracted from chromosome "<<chr<< "."<<endl;
    }


    void extract_snp(bInfo* bdata,string snplistfile)
    {
        vector<string> snplist;
        string msg="SNPs";
        read_msglist(snplistfile, snplist,msg);
        update_id_map_kp(snplist, bdata->_snp_name_map, bdata->_include);
        cout<<bdata->_include.size()<<" SNPs are extracted from ["+snplistfile+"]."<<endl;
    }


    void exclude_snp(bInfo* bdata,string snplistfile)
    {
        vector<string> snplist;
        string msg="SNPs";
        read_msglist(snplistfile, snplist,msg);
        int prev_size = bdata->_include.size();
        update_id_map_rm(snplist, bdata->_snp_name_map, bdata->_include);
        cout << prev_size - bdata->_include.size() << " SNPs are excluded from [" + snplistfile + "] and there are " << bdata->_include.size() << " SNPs remaining." << endl;
    }


    int getMaxNum(bInfo* bdata,int ldWind, vector<uint64_t> &cols)
    {
        int n=0, loopj=0, preldnum=1;
        long window=ldWind*1000;
        int maxldnum=0;
        cols.resize(bdata->_include.size()+1);
        for(int i=0;i<bdata->_include.size();i++)
        {
             int chri=bdata->_chr[bdata->_include[i]];
            int bpi=bdata->_bp[bdata->_include[i]];
            int ldnum=0;
            for(int j=loopj+1;j<bdata->_include.size();j++)
            {
                int chrj=bdata->_chr[bdata->_include[j]];
                int bpj=bdata->_bp[bdata->_include[j]];
                if(chri==chrj && abs(bpj-bpi)<=window)
                {
                    ldnum++;
                    loopj=j;
                }
                else
                {
                    break;
                }
            }
            ldnum += preldnum-1;
            preldnum=ldnum;
            cols[i+1]=cols[i]+ldnum;
            if(ldnum>maxldnum) maxldnum=ldnum;
        }
        ++maxldnum;
        while((maxldnum>>=1) != 0) n++;
        maxldnum=2<<n;
        return maxldnum;
    }


    void mu_func(bInfo* bdata, int j, vector<double> &fac) {
        int i = 0;
        bdata->_dosage_flag = 0;
        double fcount = 0.0, f_buf = 0.0;
        if (bdata->_dosage_flag) {
            for (i = 0; i < bdata->_keep.size(); i++) {
                if (bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]] < 1e5) {
                    bdata->_mu[bdata->_include[j]] += fac[i] * bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]];
                    fcount += fac[i];
                }
            }
        } else {
            for (i = 0; i < bdata->_keep.size(); i++) {
                if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
                    f_buf = (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                    if (bdata->_allele2[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) f_buf = 2.0 - f_buf;
                    bdata->_mu[bdata->_include[j]] += fac[i] * f_buf;
                    fcount += fac[i];
                }
            }
        }
        
        if (fcount > 0.0)bdata->_mu[bdata->_include[j]] /= fcount;
    }


    void calcu_mu(bInfo* bdata, bool ssq_flag)
    {
        int i = 0, j = 0;
        
        vector<double> auto_fac(bdata->_keep.size()), xfac(bdata->_keep.size()), fac(bdata->_keep.size());
        for (i = 0; i < bdata->_keep.size(); i++)
        {
            auto_fac[i] = 1.0;
            if (bdata->_sex[bdata->_keep[i]] == 1) xfac[i] = 0.5;
            else if (bdata->_sex[bdata->_keep[i]] == 2) xfac[i] = 1.0;
            fac[i] = 0.5;
        }
        
        cout << "Calculating allele frequencies ..." << endl;
        bdata->_mu.clear();
        bdata->_mu.resize(bdata->_snp_num);
        
#pragma omp parallel for
        for (j = 0; j < bdata->_include.size(); j++)
        {
            if (bdata->_chr[bdata->_include[j]]<(bdata->_autosome_num + 1)) mu_func(bdata, j, auto_fac);
            else if (bdata->_chr[bdata->_include[j]] == (bdata->_autosome_num + 1)) mu_func(bdata,j, xfac);
            else mu_func(bdata, j, fac);
        }
    }
    

    bool make_XMat_subset(bInfo* bdata, MatrixXf &X, vector<int> &snp_indx, bool divid_by_std)
    {
        if(snp_indx.empty()) return false;
        if (bdata->_mu.empty()) calcu_mu(bdata);
        
        int i = 0, j = 0, k = 0, n = (int)bdata->_keep.size(), m = (int)snp_indx.size();
        vector<double> sd_SNP(m);
        
        X.resize(n, m);
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                k = bdata->_include[snp_indx[j]];
                if (!bdata->_snp_1[k][bdata->_keep[i]] || bdata->_snp_2[k][bdata->_keep[i]]) {
                    if (bdata->_allele1[k] == bdata->_ref_A[k]) X(i,j) = bdata->_snp_1[k][bdata->_keep[i]] + bdata->_snp_2[k][bdata->_keep[i]];
                    else X(i,j) = 2.0 - (bdata->_snp_1[k][bdata->_keep[i]] + bdata->_snp_2[k][bdata->_keep[i]]);
                    X(i,j) -= bdata->_mu[k];
                }
                else X(i,j) = 0.0;
                sd_SNP[j]+=X(i,j)*X(i,j);
            }
        }

        if(divid_by_std){
            for (j = 0; j < m; j++){
                sd_SNP[j]=sd_SNP[j]/(n-1.0);
                if (fabs(sd_SNP[j]) < 1.0e-50) sd_SNP[j] = 0.0;
                else sd_SNP[j] = sqrt(1.0 / sd_SNP[j]);
            }
            for (j = 0; j < m; j++) X.col(j) = X.col(j).array() * sd_SNP[j];
        }
        
        /*
         // alternative method. a little different
        if(divid_by_std){
            vector<double> sd_SNP(m);
            for (j = 0; j < m; j++){
                k = bdata->_include[snp_indx[j]];
                sd_SNP[j] = bdata->_mu[k]*(1.0 - 0.5 * bdata->_mu[k]);
                if (fabs(sd_SNP[j]) < 1.0e-50) sd_SNP[j] = 0.0;
                else sd_SNP[j] = sqrt(1.0 / sd_SNP[j]);
            }
            for (j = 0; j < m; j++) X.col(j) = X.col(j).array() * sd_SNP[j];
        }
        */
        return true;
    }


    void makex_xVec_subset(bInfo* bdata,int j, VectorXf &x, bool resize, bool divid_by_std)
    {
        if (resize) x.resize(bdata->_keep.size());
        int k=bdata->_include[j];
        double sd_SNP=0;
        for (int i = 0; i < bdata->_keep.size(); i++)
        {
            if (!bdata->_snp_1[k][bdata->_keep[i]] || bdata->_snp_2[k][bdata->_keep[i]])
            {
                if (bdata->_allele1[k] == bdata->_ref_A[k]) x[i] = (bdata->_snp_1[k][bdata->_keep[i]] + bdata->_snp_2[k][bdata->_keep[i]]);
                else x(i) = 2.0 - (bdata->_snp_1[k][bdata->_keep[i]] + bdata->_snp_2[k][bdata->_keep[i]]);
                x(i) -= bdata->_mu[k];
            }
            else x(i) = 0.0;
            sd_SNP += x(i)*x(i);
        }
  
        if(divid_by_std)
        {
            sd_SNP /= (bdata->_keep.size()-1.0);
            if (fabs(sd_SNP) < 1.0e-50) sd_SNP = 0.0;
            else sd_SNP = sqrt(1.0 / sd_SNP);
            x*=sd_SNP;
        }
        /*
        if(divid_by_std)
        {
            double sd_SNP = bdata->_mu[k]*(1.0 - 0.5 * bdata->_mu[k]);
            if (fabs(sd_SNP) < 1.0e-50) sd_SNP = 0.0;
            else sd_SNP = sqrt(1.0 / sd_SNP);
            x*=sd_SNP;
        }
         */
    }


    void initX(bInfo* bdata,MatrixXf &X, long snpnum)
    {
        vector<int> snpids(snpnum);
        for(int i=0;i<snpnum;i++) snpids[i]=i;
        //make_XMat(bdata, snpids,X,true); //centered
        make_XMat_subset(bdata,X,snpids,true);
        
    }


    void write_smr_esi(char* outFileName, bInfo* binfo)
    {
        string epiName=string(outFileName)+".esi";
        FILE* efile=fopen(epiName.c_str(),"w");
        if (!(efile)) {
            printf("Error: Failed to open file %s.\n",epiName.c_str());
            exit(EXIT_FAILURE);
        }

        printf("Saving SNP information ...\n");
        for(int i=0;i<binfo->_include.size();i++)
        {
            string chrstr;
            if(binfo->_chr[binfo->_include[i]]==23) chrstr="X";
            else if(binfo->_chr[binfo->_include[i]]==24) chrstr="Y";
            else chrstr=atosm(binfo->_chr[binfo->_include[i]]);
            string freqstr;
            if(binfo->_mu.empty()) freqstr="NA";
            else freqstr=atos(binfo->_mu[binfo->_include[i]] /2 );
            string str=chrstr+'\t'+binfo->_snp_name[binfo->_include[i]]+'\t'+atos(0)+'\t'+atosm(binfo->_bp[binfo->_include[i]])+'\t'+binfo->_allele1[binfo->_include[i]]+'\t'+binfo->_allele2[binfo->_include[i]]+'\t'+freqstr+'\n';
            if(fputs_checked(str.c_str(),efile))
            {
                printf("ERROR: in writing file %s .\n", epiName.c_str());
            }
        }
        fclose(efile);
        printf("%ld SNPs have been saved in the file %s .\n", binfo->_include.size(), epiName.c_str());
    }


    void filter_snp_maf(bInfo* bdata,double maf)
    {
        if(bdata->_mu.empty()) calcu_mu(bdata);
        
        cout<<"Pruning SNPs with MAF > "<<maf<<" ..."<<endl;
        map<string, int> id_map_buf(bdata->_snp_name_map);
        map<string, int>::iterator iter, end=id_map_buf.end();
        int prev_size=bdata->_include.size();
        double fbuf=0.0;
        bdata->_include.clear();
        bdata->_snp_name_map.clear();
        for(iter=id_map_buf.begin(); iter!=end; iter++){
            fbuf=bdata->_mu[iter->second]*0.5;
            if(fbuf<=maf || (1.0-fbuf)<=maf) continue;
            bdata->_snp_name_map.insert(*iter);
            bdata->_include.push_back(iter->second);
        }
        if(bdata->_include.size()==0) throw("Error: No SNP is retained for analysis.");
        else{
            stable_sort(bdata->_include.begin(), bdata->_include.end());
            cout<<"After pruning SNPs with MAF > "<<maf<<", there are "<<bdata->_include.size()<<" SNPs ("<<prev_size-bdata->_include.size()<<" SNPs with MAF < "<<maf<<")."<<endl;
        }
        
    }


    void ld_calc_o2m(VectorXd &ld_v,long targetid, MatrixXd &X, bool centered)
    {
        long size=X.cols();
        long n=X.rows();
        
        VectorXd tmpX(size);
        VectorXd tmpX2(size);
        VectorXd tmpXY(size);
        
        if(centered)
        {
#pragma omp parallel for
            for(int i=0;i<size;i++){
                tmpX2[i]=X.col(i).dot(X.col(i));
                tmpXY[i]=X.col(targetid).dot(X.col(i));
            }
            float tmpY2=X.col(targetid).dot(X.col(targetid));
            ld_v=tmpXY.array()/sqrt(tmpX2.array()*tmpY2);
        }
        else
        {
#pragma omp parallel for
            for(int i=0;i<size;i++){
                tmpX[i]=X.col(i).sum();
                tmpX2[i]=X.col(i).dot(X.col(i));
                tmpXY[i]=X.col(targetid).dot(X.col(i));
            }
            float tmpY=X.col(targetid).sum();;
            float tmpY2=X.col(targetid).dot(X.col(targetid));
            ld_v=(tmpXY*n-tmpX*tmpY).array()/sqrt((tmpX2.array()*n-tmpX.array()*tmpX.array())*(tmpY2*n-tmpY*tmpY));
        }
        
    }


    void extract_snp_kb(bInfo* bdata,string rsnames, int windInKb)
    {
        map<string, int>::iterator iter;
        iter=bdata->_snp_name_map.find(rsnames);
        if(iter==bdata->_snp_name_map.end())
        {
            printf("ERROR: Can't find SNP %s.\n",  rsnames.c_str());
            exit(EXIT_FAILURE);
        }
        long idx=iter->second;
        int snpbp=bdata->_bp[idx];
        int snpchr=bdata->_chr[idx];
        int upbound=snpbp+windInKb*1000;
        int tmpint=snpbp-windInKb*1000;
        int lowbound=tmpint>0?tmpint:0;
        vector<string> snplist;
        for(int i=0;i<bdata->_include.size();i++)
        {
            tmpint=bdata->_include[i];
            if(bdata->_chr[tmpint]==snpchr && bdata->_bp[tmpint]>=lowbound && bdata->_bp[tmpint]<=upbound) snplist.push_back(bdata->_snp_name[tmpint]);
        }
        update_id_map_kp(snplist, bdata->_snp_name_map, bdata->_include);
        printf("%ld SNPs including %s are extracted.\n",bdata->_include.size(), rsnames.c_str());
    }


    void ld_report(char* outFileName, char* bFileName,char* indilstName, char* indilst2remove,char* snplstName, char* snplst2exclde,int chr,char* rs, double maf, bool ldr, bool ldr2, int ldWind)
    {
        if(!(ldr || ldr2))
        {
            printf("Please specify --r or --r2 \n");
            exit(EXIT_FAILURE);
        }
        bool bitmod=true;
        bInfo bdata;
        MatrixXf X;
        vector<uint64_t> cols;
        vector<float> lds;
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        read_bimfile(&bdata, string(bFileName)+".bim");
        if(chr) extract_snp(&bdata, chr);
        if(rs != NULL) extract_snp_kb(&bdata,rs,ldWind);
        if(snplstName != NULL) extract_snp(&bdata, snplstName);
        if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0) filter_snp_maf(&bdata, maf);
        int m=getMaxNum(&bdata,ldWind, cols);
        if(m==1)
        {
            printf("No SNP pair included in %d Kb.\n",ldWind);
            exit(EXIT_FAILURE);
        }
        if(m>bdata._include.size())
        {
            m=(int)bdata._include.size();
            bitmod=false;
        }
        write_smr_esi(outFileName, &bdata);
        string bldname=string(outFileName)+".bld";
        FILE* outfile=fopen(bldname.c_str(), "wb");
        if (!(outfile)) {
            printf("Error: Failed to open file %s.\n",bldname.c_str());
            exit(EXIT_FAILURE);
        }
        vector<int> reserved(RESERVEDUNITS);
        if(ldr) reserved[0]=0;
        else if(ldr2) reserved[0]=1;
        reserved[1]=(int)bdata._keep.size();
        reserved[2]=(int)bdata._include.size();
        reserved[3]=ldWind;
        for(int i=4;i<RESERVEDUNITS;i++) reserved[i]=-9;
        fwrite(&reserved[0],sizeof(int), RESERVEDUNITS, outfile);
        uint64_t valnum=cols[bdata._include.size()];
        fwrite(&valnum,sizeof(uint64_t), 1, outfile);
        fwrite(&cols[0],sizeof(uint64_t), bdata._include.size()+1, outfile);
        long window=ldWind*1000, vc=0;
        double cr=0;
        for(int i=0;i<bdata._include.size();i++)
        {
            progress(i, cr, (int)bdata._include.size());

            if(X.size()==0) initX(&bdata, X,m);
            int start =-9;
            if(bitmod) start = (i & (m-1));
            else start = i % m;
           
            VectorXf x=X.col(start);
            
            int chri=bdata._chr[bdata._include[i]];
            int bpi=bdata._bp[bdata._include[i]];
            
            //clock_t begin_time = clock();
            VectorXf ldv=X.col(start)/(X.rows()-1);
            ldv=X.transpose()*ldv;
            if(ldr2) ldv=ldv.array()*ldv.array();
            int st=-9, ed=-9;
            for(int j=1;j<m && i+j<bdata._include.size();j++)
            {
                
                int chrj=bdata._chr[bdata._include[i+j]];
                int bpj=bdata._bp[bdata._include[i+j]];
                int cur = -9;
                if(bitmod) cur=((start+j) & (m-1));
                else cur= (start+j) % m ;
                if(chri==chrj && abs(bpj-bpi)<=window)
                {
                    if(st<0) st=cur;
                    ed=cur;
                    if(ed==m-1)
                    {
                        fwrite(&ldv(st),sizeof(float), ed-st+1, outfile);
                        st=-9;
                        ed=-9;
                    }
                    vc++;
                }
            }
            if(ed>=0 && st>=0) fwrite(&ldv(st),sizeof(float), ed-st+1, outfile);
            
            //printf(" cost2: %f ms.\n",float( clock () - begin_time ) /  1000);
            //begin_time = clock();
            if(i+m<bdata._include.size())
            {
                makex_xVec_subset(&bdata, i+m, x, false, true);
                X.col(start)=x;
            }
            //printf(" cost3: %f ms.\n",float( clock () - begin_time ) /  1000);
        }
        if(vc!=valnum)
        {
            printf("Error: predicted number vs observed number: %ld, %llu.\n",vc,valnum);
            printf("Please repot this bug.\n");
            exit(EXIT_FAILURE);
            
        }
        fclose(outfile);
        printf("LD information is saved in the binary file %s.\n",bldname.c_str());

    }


    void read_ld_esifile(ldInfo* ldinfo, char* esiFileName)
    {
        vector<string> strlist;
        uint32_t line_idx = 0;
        int colnum=7;
        FILE* esifile=fopen(esiFileName, "r");
        if (!(esifile)) {
            printf("Error: Failed to open file %s.\n",esiFileName);
            exit(EXIT_FAILURE);
        }
        printf("Reading SNP information from %s ...\n", esiFileName);
        ldinfo->_esi_chr.clear();
        ldinfo->_esi_rs.clear();
        ldinfo->_esi_gd.clear();
        ldinfo->_esi_bp.clear();
        ldinfo->_esi_allele1.clear();
        ldinfo->_esi_allele2.clear();
        ldinfo->_esi_include.clear();
        ldinfo->_snp_name_map.clear();
        ldinfo->_esi_freq.clear();
        
        bool chrwarning=false;
        bool feqwarning=false, allele1warning=false, allele2warning=false;
        bool orienwarning=false;
        char Tbuf[MAX_LINE_SIZE];
        while(fgets(Tbuf, MAX_LINE_SIZE, esifile))
        {
            split_string(Tbuf, strlist, ", \t\n");
            if(Tbuf[0]=='\0') {
                printf("ERROR: Line %u is blank.\n", line_idx);
                exit(EXIT_FAILURE);
            }
            if(strlist.size()<colnum-1)
            {
                printf("ERROR: Line %u has less than %d items.\n", line_idx,colnum-1);
                exit(EXIT_FAILURE);
            } else if(strlist.size()==colnum-1) {
                if(!feqwarning) {
                    printf("WARNING: Maybe this is an old .esi file which doesn't contain frequency information. \n");
                    feqwarning=true;
                }
            } else if(strlist.size()>colnum) {
                //printf("WARNING: Line %u has more than %d items. The first %d columns would be used. \n", line_idx,colnum,colnum);
            }
            
            ldinfo->_snp_name_map.insert(pair<string,int>(strlist[1],line_idx));
            if(ldinfo->_snp_name_map.size()==line_idx)
            {
                printf("ERROR: Duplicate SNP : %s.\n", strlist[1].c_str());
                exit(EXIT_FAILURE);
            }
            if(strlist[0]=="X" || strlist[0]=="x") ldinfo->_esi_chr.push_back(23);
            else if(strlist[0]=="Y" || strlist[0]=="y") ldinfo->_esi_chr.push_back(24);
            else if(strlist[0]=="NA" || strlist[0]=="na"){
                ldinfo->_esi_chr.push_back(-9);
                if(!chrwarning) {
                    printf("WARNING: At least one SNP chr is missing.\n");
                    chrwarning=true;
                }
            } else if (atoi(strlist[0].c_str())==0 ) {
                //printf("WARNING: unrecongized chromosome found. This chromosome is set to 0:\n");
                //printf("%s\n",Tbuf);
                ldinfo->_esi_chr.push_back(atoi(strlist[0].c_str()));
            } else if ( atoi(strlist[0].c_str())>24 || atoi(strlist[0].c_str())<0) {
                //printf("WARNING: abmormal chromosome found:\n");
                //printf("%s\n",Tbuf);
                ldinfo->_esi_chr.push_back(atoi(strlist[0].c_str()));
            } else ldinfo->_esi_chr.push_back(atoi(strlist[0].c_str()));
            
            if(strlist[1]=="NA" || strlist[1]=="na") {
                printf("ERROR: NA SNP ID found:\n");
                printf("%s\n",Tbuf);
                exit(EXIT_FAILURE);
            }
            ldinfo->_esi_rs.push_back(strlist[1]);
            ldinfo->_esi_gd.push_back(atoi(strlist[2].c_str()));
            if(strlist[3]=="NA" || strlist[3]=="na") ldinfo->_esi_bp.push_back(-9);
            else ldinfo->_esi_bp.push_back(atoi(strlist[3].c_str()));
            if(strlist[4]=="NA" || strlist[4]=="na") {
                if(!allele1warning) {
                    //printf("WARNING: At least one reference allele is missing.\n");
                    allele1warning=true;
                }
            }
            to_upper(strlist[4]);
            ldinfo->_esi_allele1.push_back(strlist[4].c_str());
            if(strlist[5]=="NA" || strlist[5]=="na") {
                if(!allele2warning) {
                    //printf("WARNING: At least one alternative allele is missing.\n");
                    allele2warning=true;
                }
            }
            to_upper(strlist[5]);
            ldinfo->_esi_allele2.push_back(strlist[5].c_str());
            if(strlist.size()==colnum)
            {
                if(strlist[6]=="NA" || strlist[6]=="na"){
                    if(!orienwarning){
                        //printf("WARNING: frequency is \"NA\" in one or more rows.\n");
                        orienwarning=true;
                    }
                    ldinfo->_esi_freq.push_back(-9);
                } else {
                    ldinfo->_esi_freq.push_back(atof(strlist[6].c_str()));
                }
            } else {
                ldinfo->_esi_freq.push_back(-9);
            }
            ldinfo->_esi_include.push_back(line_idx);
            line_idx++;
        }
        ldinfo->_snpNum =line_idx;
        fclose(esifile);
        printf("%llu SNPs to be included from  %s .\n", ldinfo->_snpNum, esiFileName);
    }


    void extract_ld_esi_by_chr(ldInfo* ldinfo, int snpchr)
    {
        vector<int> newIcld;
        ldinfo->_snp_name_map.clear();
        for(int i=0;i<ldinfo->_esi_include.size();i++)
        {
            int tmpint=ldinfo->_esi_include[i];
            if(ldinfo->_esi_chr[tmpint]==snpchr ) {
                newIcld.push_back(tmpint);
                ldinfo->_snp_name_map.insert(pair<string,int>(ldinfo->_esi_rs[tmpint],tmpint));
            }
        }
        ldinfo->_esi_include.swap(newIcld);
        printf("%ld SNPs are extracted from chromosome %d. \n",ldinfo->_esi_include.size(),snpchr);
    }


    void extract_ld_esi_snps(ldInfo* ldinfo, string snplstName)
    {
        vector<string> snplist;
        string msg = "SNPs";
        read_msglist(snplstName, snplist, msg);
        update_map_kp(snplist, ldinfo->_snp_name_map, ldinfo->_esi_include);
        printf("%ld SNPs are extracted from %s.\n", ldinfo->_esi_include.size() ,snplstName.c_str() );
    }


    void extract_ld_esi_snps(ldInfo* ldinfo, string snp, int Wind)
    {
        int bp=-9;
        int chr=-9;
        map<string, int>::iterator iter;
                    iter=ldinfo->_snp_name_map.find(snp);
            if(iter==ldinfo->_snp_name_map.end())
            {
                printf("ERROR: Can't find SNP %s.\n",  snp.c_str());
                exit(EXIT_FAILURE);
            }
            bp=ldinfo->_esi_bp[iter->second];
            chr=ldinfo->_esi_chr[iter->second];
            if(chr < 0) {
                printf("ERROR: Missing Chromosome found of SNP %s.\n", snp.c_str());
                exit(EXIT_FAILURE);
            }
            if(bp < 0) {
                printf("ERROR: Missing BP found of SNP %s.\n", snp.c_str());
                exit(EXIT_FAILURE);
            }
        
        int upbound=bp+Wind*1000;
        int tmpint=bp-Wind*1000;
        int lowbound=tmpint>0?tmpint:0;
        vector<int> newIcld;
        ldinfo->_snp_name_map.clear();
        for(int i=0;i<ldinfo->_esi_include.size();i++)
        {
            tmpint=ldinfo->_esi_include[i];
            if(ldinfo->_esi_chr[tmpint]==chr && ldinfo->_esi_bp[tmpint]>=lowbound && ldinfo->_esi_bp[tmpint]<=upbound) {
                newIcld.push_back(tmpint);
                ldinfo->_snp_name_map.insert(pair<string,int>(ldinfo->_esi_rs[tmpint],tmpint));
            }
        }
        ldinfo->_esi_include.swap(newIcld);
        printf("%ld SNPs are extracted from the region: %d Kb around SNP %s.\n", ldinfo->_esi_include.size() ,Wind, snp.c_str());
    }


    void extract_ld_esi_single_snp(ldInfo* ldinfo, string snprs)
    {
        map<string, int>::iterator iter;
        iter=ldinfo->_snp_name_map.find(snprs);
        if(iter==ldinfo->_snp_name_map.end())
        {
            printf("ERROR: Can't find SNP %s.\n", snprs.c_str());
            exit(EXIT_FAILURE);
        }
        long idx=iter->second;
        ldinfo->_esi_include.clear();
        ldinfo->_esi_include.push_back((int)idx);
        ldinfo->_snp_name_map.clear();
        ldinfo->_snp_name_map.insert(pair<string,int>(snprs,idx));
        printf("SNP %s is extracted.\n",  snprs.c_str());
    }


    void extract_ld_esi_snps(ldInfo* ldinfo, string fromsnprs, string tosnprs)
    {
        map<string, int>::iterator iter;
        iter=ldinfo->_snp_name_map.find(fromsnprs);
        if(iter==ldinfo->_snp_name_map.end())
        {
            printf("ERROR: Can't find probe %s.\n",  fromsnprs.c_str());
            exit(EXIT_FAILURE);
        }
        int fromsnpbp=ldinfo->_esi_bp[iter->second];
        int snpchr=ldinfo->_esi_chr[iter->second];
        
        iter=ldinfo->_snp_name_map.find(tosnprs);
        if(iter==ldinfo->_snp_name_map.end())
        {
            printf("ERROR: Can't find probe %s.\n",  tosnprs.c_str());
            exit(EXIT_FAILURE);
        }
        int tosnpbp=ldinfo->_esi_bp[iter->second];
        int tosnpchr=ldinfo->_esi_chr[iter->second];
        if(tosnpchr != snpchr)
        {
            printf("ERROR: SNP %s and SNP %s are not from the same chromosome.\n", fromsnprs.c_str(), tosnprs.c_str());
            exit(EXIT_FAILURE);
        }
        if(fromsnpbp>tosnpbp)
        {
            int tmp=fromsnpbp;
            fromsnpbp=tosnpbp;
            tosnpbp=tmp;
        }
        
        vector<string> snplst;
        for(int i=0;i<ldinfo->_esi_include.size();i++)
        {
            int tmpint=ldinfo->_esi_include[i];
            if(ldinfo->_esi_chr[tmpint]==snpchr && ldinfo->_esi_bp[tmpint]>=fromsnpbp && ldinfo->_esi_bp[tmpint]<=tosnpbp) snplst.push_back(ldinfo->_esi_rs[tmpint]);
        }
        update_map_kp(snplst, ldinfo->_snp_name_map, ldinfo->_esi_include);
        printf("%ld SNPs are extracted from SNP %s to SNP %s.\n",  ldinfo->_esi_include.size(), fromsnprs.c_str(), tosnprs.c_str());
    }


    void extract_ld_esi_snps(ldInfo* ldinfo, int chr, int fromsnpkb, int tosnpkb)
    {
        int fromsnpbp=fromsnpkb*1000;
        int tosnpbp=tosnpkb*1000;
        
        if(fromsnpbp>tosnpbp)
        {
            int tmp=fromsnpbp;
            fromsnpbp=tosnpbp;
            tosnpbp=tmp;
        }
        
        vector<int> newIcld;
        ldinfo->_snp_name_map.clear();
        for(int i=0;i<ldinfo->_esi_include.size();i++)
        {
            int tmpint=ldinfo->_esi_include[i];
            if( ldinfo->_esi_chr[tmpint]==chr &&ldinfo->_esi_bp[tmpint]>=fromsnpbp && ldinfo->_esi_bp[tmpint]<=tosnpbp) {
                newIcld.push_back(tmpint);
                ldinfo->_snp_name_map.insert(pair<string,int>(ldinfo->_esi_rs[tmpint],tmpint));
            }
        }
        ldinfo->_esi_include.swap(newIcld);
        printf("%ld SNPs are extracted from SNP BP:  %dKb to SNP BP: %dKb on chromosome %d.\n",ldinfo->_esi_include.size(),fromsnpkb,tosnpkb, chr);
    }


    void exclude_ld_esi_snp(ldInfo* ldinfo, string snplstName)
    {
        vector<string> snplist;
        string msg="SNPs";
        read_msglist(snplstName, snplist,msg);
        long pre_num=ldinfo->_esi_include.size();
        update_map_rm(snplist, ldinfo->_snp_name_map, ldinfo->_esi_include);
        printf("%ld SNPs are excluded from %s and there are %ld SNPs remaining.\n",  pre_num - ldinfo->_esi_include.size(), snplstName.c_str(),ldinfo->_esi_include.size());
    }


    void exclude_ld_esi_snps(ldInfo* ldinfo, string snprs2exclde)
    {
        vector<string> snplist;
        snplist.push_back(snprs2exclde);
        update_map_rm(snplist, ldinfo->_snp_name_map, ldinfo->_esi_include);
        printf("SNP %s are excluded and there are %ld probe remaining.\n", snprs2exclde.c_str(), ldinfo->_esi_include.size());
    }


    void ld_esi_man(ldInfo* ldinfo,char* snplstName, char* snplst2exclde, int chr,char* snprs, char* fromsnprs, char* tosnprs,int snpWind,bool snpwindFlag, int fromsnpkb, int tosnpkb,char* snprs2exclde)
    {
        string logstr;
        int flags4snp=0;
        if(snplstName != NULL) flags4snp++;
        if(snprs != NULL) flags4snp++;
        if(fromsnprs!=NULL) flags4snp++;
        if(fromsnpkb>=0) flags4snp++;
        if(flags4snp>1)
        {
            printf("WARNING: Flags for SNPs in this section are mutual exclusive. The priority order (from high to low) is: --extract-snp, --snp-wind, --snp, --from(to)--snp, --from(to)-snp-kb.\n");
        }
        if(chr>0)
        {
            extract_ld_esi_by_chr(ldinfo, chr);
        }
        
        if (snplstName != NULL) extract_ld_esi_snps(ldinfo, snplstName);
        else if (snpwindFlag)
        {
            if(snprs==NULL)
            {
                printf("ERROR: please specify the SNP name by --snp when using --snp-wind.\n");
               exit(EXIT_FAILURE);
            }
            extract_ld_esi_snps(ldinfo, snprs, snpWind);
        }
        else if(snprs!=NULL)
        {
            extract_ld_esi_single_snp(ldinfo, snprs);
        }
        else if(fromsnprs!=NULL)
        {
            if(tosnprs==NULL)
            {
                printf("ERROR: please specify the SNP name by --to-snp.\n");
                exit(EXIT_FAILURE);
            }
            extract_ld_esi_snps(ldinfo, fromsnprs, tosnprs);
        }
        else if(fromsnpkb>=0)
        {
            
            if(fromsnpkb>=0 && chr==0) {
                printf("ERROR: please specify the chromosome by --chr.\n");
                exit(EXIT_FAILURE);
            }
            
            if(tosnpkb<0)
            {
                printf("ERROR: SNP BP can't be negative.\n");
                exit(EXIT_FAILURE);
            }
            extract_ld_esi_snps(ldinfo, chr, fromsnpkb, tosnpkb);
        }
        if(snplst2exclde!=NULL) exclude_ld_esi_snp(ldinfo,snplst2exclde);
        if(snprs2exclde!=NULL) exclude_ld_esi_snps(ldinfo,snprs2exclde);
        
    }


    void get_BlddHeaders(char* blddFileName, vector<int> &headers)//
    {
        headers.resize(RESERVEDUNITS);
        FILE* bld=fopen(blddFileName,"rb");
        if(bld == NULL)
        {
            printf("Error: can't open file %s.\n",blddFileName);
            exit(EXIT_FAILURE);
        }
        if(fread(&headers[0], sizeof(int),RESERVEDUNITS, bld)<1)
        {
            printf("ERROR: File %s read failed!\n", blddFileName);
            exit(EXIT_FAILURE);
        }
        fclose(bld);
    }


    void fetch_ld_by_id(ldInfo* ldinfo,FILE* ldfprt, vector<uint32_t> &curId, int sid, vector<float> &ld)
    {
        ld.resize(curId.size());
        int testid=ldinfo->_esi_include[curId[sid]];
        uint64_t valSTART=RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) + (ldinfo->_snpNum+1)*sizeof(uint64_t);
        for(int i=0;i<curId.size();i++)
        {
            int toid=ldinfo->_esi_include[curId[i]];
            if(testid==toid) ld[i]=1;
            else if(testid>toid)
            {
                long poss=ldinfo->_cols[toid];
                fseek( ldfprt, (poss+testid-toid-1)*sizeof(float)+valSTART, SEEK_SET );
                ld[i]=readfloat(ldfprt);
            }
            else
            {
                long poss=ldinfo->_cols[testid];
                fseek( ldfprt, (poss+toid-testid-1)*sizeof(float)+valSTART, SEEK_SET );
                ld[i]=readfloat(ldfprt);
            }
        }
    }


    void fetch_ld_by_id(ldInfo* ldinfo,FILE* ldfprt, int sid, vector<float> &ld)
    {  // only the SNPs with bigger pb would be extracted!
        uint64_t valSTART=RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) + ldinfo->_snpNum*sizeof(uint64_t);
        uint64_t poss=ldinfo->_cols[sid];
        uint64_t post=ldinfo->_cols[sid+1];
        long num=post-poss;
        ld.resize(num);
        fseek( ldfprt, poss*sizeof(float)+valSTART, SEEK_SET );
        if(fread(&ld[0], sizeof(float),num, ldfprt)!=num)
        {
            printf("ERROR: File read failed!\n");
            exit(EXIT_FAILURE);
        }

    }


    void fetch_ld_by_snps(ldInfo* ldinfo,FILE* ldfprt, string rs, vector<float> &ld)
    {
        map<string,int>::iterator iter;
        iter=ldinfo->_snp_name_map.find(rs);
        if(iter==ldinfo->_snp_name_map.end())
        {
            printf("ERROR: can't find SNP %s.\n",rs.c_str());
             exit(EXIT_FAILURE);
        }
        int sid=iter->second;
        uint64_t valSTART=RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) + ldinfo->_snpNum*sizeof(uint64_t);
        uint64_t poss=ldinfo->_cols[sid];
        uint64_t post=ldinfo->_cols[sid+1];
        long num=post-poss;
        ld.resize(num);
        fseek( ldfprt, poss*sizeof(float)+valSTART, SEEK_SET );
        if(fread(&ld[0], sizeof(float),num, ldfprt)!=num)
        {
            printf("ERROR: File read failed!\n");
            exit(EXIT_FAILURE);
        }
    }


    void lookup(char* outFileName, char* bldFileName, char* snplstName, char* snplst2exclde,int chr,char* snprs, char* snprs2exclde, char* fromsnprs, char* tosnprs,int snpWind, bool snpWindflg, int fromsnpkb, int tosnpkb, int ld_wind)
    {
        ldInfo ldinfo;
        if(bldFileName == NULL)
        {
            printf("Error: please input the ld information by the option --bld.\n");
            exit(EXIT_FAILURE);
        }
        char inputname[FNAMESIZE];
        memcpy(inputname,bldFileName,strlen(bldFileName)+1);
        char* suffix=inputname+strlen(bldFileName);
        memcpy(suffix,".esi",5);
        read_ld_esifile(&ldinfo, inputname);
        
        memcpy(suffix,".bld",5);
        vector<int> headers;
        headers.resize(RESERVEDUNITS);
        FILE* bld=fopen(inputname,"rb");
        if(bld == NULL)
        {
            printf("Error: can't open file %s.\n",inputname);
            exit(EXIT_FAILURE);
        }
        if(fread(&headers[0], sizeof(int),RESERVEDUNITS, bld)<1)
        {
            printf("ERROR: File %s read failed!\n", inputname);
            exit(EXIT_FAILURE);
        }
        int ldwind=headers[3];
        if(ldwind>ld_wind) ldwind=ld_wind;
        ld_esi_man(&ldinfo, snplstName, snplst2exclde,chr,  snprs,  fromsnprs,  tosnprs, snpWind, snpWindflg, fromsnpkb,  tosnpkb,snprs2exclde);
        if(ldinfo._esi_include.size()==0)
        {
            printf("Error: no SNP included.\n");
            exit(EXIT_FAILURE);
        }
        int indicator = headers[0];
        if(indicator==0) printf("\nReading ld r from binary file %s...\n", inputname);
        else printf("\nReading ld r-squared from binary file %s...\n", inputname);
        uint64_t valnum=readuint64(bld), colNum=ldinfo._snpNum+1;
        uint64_t cur_pos = ftell( bld );
        fseek( bld, 0L, SEEK_END );
        uint64_t size_file = ftell( bld );
        fseek( bld, cur_pos, SEEK_SET );
        if( size_file - (RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valnum*sizeof(float)) != 0) {
            printf("ERROR: File %s is broken!\n", inputname);
            exit(EXIT_FAILURE);
        }
        
        ldinfo._cols.resize(colNum);
        if(fread(&ldinfo._cols[0], sizeof(uint64_t),colNum, bld)<1)
        {
            printf("ERROR: File %s read failed!\n", inputname);
            exit(EXIT_FAILURE);
        }
        long maxnum=0;
        for(int i=1;i<colNum;i++) {
            long numtmp=ldinfo._cols[i]-ldinfo._cols[i-1];
            if(numtmp>maxnum) maxnum=numtmp;
        }
        FILE* outfile=fopen(outFileName, "w");
        if (!(outfile)) {
            printf("Error: Failed to open file %s.\n",outFileName);
            exit(EXIT_FAILURE);
        }
        string tmpstr="CHR_A\tBP_A\tSNP_A\tCHR_B\tBP_B\tSNP_B\tR\n";
        if(indicator) tmpstr="CHR_A\tBP_A\tSNP_A\tCHR_B\tBP_B\tSNP_B\tR2\n";
        fputs(tmpstr.c_str(),outfile);
        uint64_t valSTART=RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) + colNum*sizeof(uint64_t);
        long wcount=0;
        float* buffer = (float*) malloc (sizeof(float)*maxnum);
        if (buffer == NULL) {
            printf("ERROR: memory allocation failed to read %s.\n", inputname);
            exit(EXIT_FAILURE);
        }
        if(ldinfo._esi_include.size()==1)
        {
            vector<int> snplist;
            int sid=ldinfo._esi_include[0];
            int chri=ldinfo._esi_chr[sid];
            int bpi=ldinfo._esi_bp[sid];
            string rsi=ldinfo._esi_rs[sid];
            for(int i=sid-1;i>=0;i--)
            {
                int chrj=ldinfo._esi_chr[i];
                int bpj=ldinfo._esi_bp[i];
                if(chri==chrj && bpi-bpj<ldwind*1000) snplist.push_back(i);
            }
            for(int i=(int)snplist.size()-1;i>=0;i--)
            {
                int sidi=snplist[i];
                int chrj=ldinfo._esi_chr[sidi];
                int bpj=ldinfo._esi_bp[sidi];
                string rsj=ldinfo._esi_rs[sidi];
                uint64_t poss=ldinfo._cols[sidi];
                fseek( bld, (poss+sid-sidi-1)*sizeof(float)+valSTART, SEEK_SET );
                float ldv=readfloat(bld);
                string tmpstr = atos(chrj)+'\t'+atos(bpj)+ '\t'+rsj+'\t'+ atos(chri)+'\t'+atos(bpi)+ '\t'+rsi+'\t'+atos(ldv)+'\n';
                fputs(tmpstr.c_str(),outfile);
                wcount++;
            }
        }
        for(int i=0;i<ldinfo._esi_include.size();i++)
        {
            
            int sid=ldinfo._esi_include[i];
            int chri=ldinfo._esi_chr[sid];
            int bpi=ldinfo._esi_bp[sid];
            string rsi=ldinfo._esi_rs[sid];
            uint64_t poss=ldinfo._cols[sid];
            uint64_t post=ldinfo._cols[sid+1];
            long num=post-poss;
            fseek( bld, poss*sizeof(float)+valSTART, SEEK_SET );
            if(fread(buffer, sizeof(float),num, bld)!=num)
            {
                printf("ERROR: File %s read failed!\n", inputname);
                exit(EXIT_FAILURE);
            }
            for(int j=0;j<num;j++)
            {
                int sid2=sid+j+1;
                int chrj=ldinfo._esi_chr[sid2];
                int bpj=ldinfo._esi_bp[sid2];
                string rsj=ldinfo._esi_rs[sid2];
                float ldv=buffer[j];
                string tmpstr = atos(chri)+'\t'+atos(bpi)+ '\t'+rsi+'\t'+atos(chrj)+'\t'+atos(bpj)+ '\t'+rsj+'\t'+atos(ldv)+'\n';
                fputs(tmpstr.c_str(),outfile);
                wcount++;
            }
        }
        free(buffer);
        fclose(bld);
        fclose(outfile);
        printf("%ld pairs of ld value are saved in the file %s.\n",wcount,outFileName);
        
    }
    

    void check_autosome(bInfo* bdata) {
        for (int i = 0; i < bdata->_include.size(); i++) {
            if (bdata->_chr[bdata->_include[i]] > bdata->_autosome_num) throw ("Error: this option is for the autosomal SNPs only. Please check the option --autosome.");
        }
    }


    void get_ld_blk_pnt(bInfo* bdata, vector<int> &brk_pnt1, vector<int> &brk_pnt2, vector<int> &brk_pnt3, int wind_bp, int wind_snp)
    {
        unsigned long i = 0, j = 0, k = 0, m = bdata->_include.size();
        
        brk_pnt1.clear();
        brk_pnt1.push_back(0);
        bool chr_start = true;
        for (i = 1, j = 0; i < m; i++) {
            if (i == (m - 1)){
                if(chr_start
                   || ((bdata->_bp[bdata->_include[i]] - bdata->_bp[bdata->_include[brk_pnt1[j]]] > 0.5 * wind_bp)
                       && (i - brk_pnt1[j] > 0.5 * wind_snp))) brk_pnt1.push_back(m - 1);
                else brk_pnt1[j - 1] = brk_pnt1[j] = m - 1;
            }
            else if (bdata->_chr[bdata->_include[i]] != bdata->_chr[bdata->_include[brk_pnt1[j]]] || bdata->_bp[bdata->_include[i]] - bdata->_bp[bdata->_include[i-1]] > 1e6) {
                if(chr_start
                   || ((bdata->_bp[bdata->_include[i-1]] - bdata->_bp[bdata->_include[brk_pnt1[j]]] > 0.5 * wind_bp)
                       && (i - 1 - brk_pnt1[j] > 0.5 * wind_snp))){
                       brk_pnt1.push_back(i - 1);
                       j++;
                       brk_pnt1.push_back(i);
                       j++;
                   }
                else{
                    brk_pnt1[j - 1] = i - 1;
                    brk_pnt1[j] = i;
                }
                chr_start = true;
            }
            else if ((bdata->_bp[bdata->_include[i]] - bdata->_bp[bdata->_include[brk_pnt1[j]]] > wind_bp) && (i - brk_pnt1[j] >= wind_snp)) {
                chr_start = false;
                brk_pnt1.push_back(i - 1);
                j++;
                brk_pnt1.push_back(i);
                j++;
            }
        }
        stable_sort(brk_pnt1.begin(), brk_pnt1.end());
        brk_pnt1.erase(unique(brk_pnt1.begin(), brk_pnt1.end()), brk_pnt1.end());
        
        brk_pnt2.clear();
        brk_pnt3.clear();
        for (i = 1; i < brk_pnt1.size() && brk_pnt1.size() > 2; i++) {
            if ((bdata->_chr[bdata->_include[brk_pnt1[i - 1]]] == bdata->_chr[bdata->_include[brk_pnt1[i]]]) && (brk_pnt1[i] - brk_pnt1[i - 1] > 1)) {
                int i_buf = (brk_pnt1[i - 1] + brk_pnt1[i]) / 2;
                brk_pnt2.push_back(i_buf);
                brk_pnt2.push_back(i_buf + 1);
                brk_pnt3.push_back(brk_pnt1[i]);
                brk_pnt3.push_back(brk_pnt1[i]);
            }
        }
    }
    

    void calcu_ld_blk_split(bInfo* bdata, int size, int size_limit, MatrixXf &X_sub, VectorXf &ssx_sqrt_i_sub, double rsq_cutoff, VectorXf &rsq_size, VectorXf &mean_rsq_sub, VectorXf &max_rsq_sub, int s1, int s2, bool second)
    {
        int i = 0, j = 0, k = 0, m = 0, n = bdata->_keep.size();
        vector<int> brk_pnt_sub;
        brk_pnt_sub.push_back(0);
        for (i = size_limit; i < size - size_limit; i += size_limit) {
            brk_pnt_sub.push_back(i - 1);
            brk_pnt_sub.push_back(i);
            j = i;
        }
        j = (size - j) / 2 + j;
        brk_pnt_sub.push_back(j - 1);
        brk_pnt_sub.push_back(j);
        brk_pnt_sub.push_back(size - 1);
        
        for (i = 0; i < brk_pnt_sub.size() - 1; i++) {
            int size_sub = brk_pnt_sub[i + 1] - brk_pnt_sub[i] + 1;
            if (size_sub < 3) continue;
            
            VectorXf ssx_sqrt_i_sub_sub = ssx_sqrt_i_sub.segment(brk_pnt_sub[i], size_sub);
            MatrixXf rsq_sub_sub = X_sub.block(0,brk_pnt_sub[i],n,size_sub).transpose() * X_sub;
            VectorXf rsq_size_sub(size_sub), mean_rsq_sub_sub(size_sub), max_rsq_sub_sub = VectorXf::Constant(size_sub, -1.0);
            
            for (j = 0; j < size_sub; j++) {
                unsigned long s = j + brk_pnt_sub[i];
                rsq_size_sub[j] = 0.0;
                mean_rsq_sub_sub[j] = 0.0;
                for (k = 0; k < size; k++) {
                    if (second) {
                        if (s <= s1 && k <= s1) continue;
                        if (s >= s2 && k >= s2) continue;
                    }
                    if (k == s) continue;
                    rsq_sub_sub(j,k) *= (ssx_sqrt_i_sub_sub[j] * ssx_sqrt_i_sub[k]);
                    rsq_sub_sub(j,k) = rsq_sub_sub(j,k) * rsq_sub_sub(j,k);
                    if (rsq_sub_sub(j,k) >= rsq_cutoff) {
                        
                        mean_rsq_sub_sub[j] += rsq_sub_sub(j,k);
                        rsq_size_sub[j] += 1.0;
                    }
                    if(rsq_sub_sub(j,k) > max_rsq_sub_sub[j]) max_rsq_sub_sub[j] = rsq_sub_sub(j,k);
                }
                
                if (rsq_size_sub[j] > 0.0) mean_rsq_sub_sub[j] /= rsq_size_sub[j];
            }
            
            for (j = 0, k = brk_pnt_sub[i]; j < size_sub; j++, k++) {
                mean_rsq_sub[k] = mean_rsq_sub_sub[j];
                rsq_size[k] = rsq_size_sub[j];
                max_rsq_sub[k] = max_rsq_sub_sub[j];
            }
        }
    }
    
    
    void calcu_ld_blk(bInfo* bdata, vector<int> &brk_pnt, vector<int> &brk_pnt3, VectorXf &mean_rsq, VectorXf &snp_num, VectorXf &max_rsq, bool second, double rsq_cutoff)
    {
        int i = 0, j = 0, k = 0, s1 = 0, s2 = 0, n = bdata->_keep.size(), m = bdata->_include.size(), size = 0, size_limit = 10000;
        int ttlsize=0;
        for (i = 0; i < brk_pnt.size() - 1; i++)
        {
            if (bdata->_chr[bdata->_include[brk_pnt[i]]] != bdata->_chr[bdata->_include[brk_pnt[i + 1]]]) continue;
            size = brk_pnt[i + 1] - brk_pnt[i] + 1;
            if (size < 3) continue;
            if (second) {
                s1 = brk_pnt3[i] - brk_pnt[i];
                s2 = s1 + 1;
            }
            else {
                s1 = 0;
                s2 = size - 1;
            }
            
            VectorXf rsq_size(size), mean_rsq_sub(size), max_rsq_sub = VectorXf::Constant(size, -1.0);
            ttlsize+=size;
            // make genotype matrix
            vector<int> snp_indx(size);
            for (j = brk_pnt[i], k = 0; j <= brk_pnt[i + 1]; j++, k++) snp_indx[k] = j;
            MatrixXf X_sub;
            make_XMat_subset(bdata,X_sub, snp_indx, true);
            VectorXf ssx_sqrt_i_sub(size);
            for (j = 0; j < size; j++){
                ssx_sqrt_i_sub[j] = X_sub.col(j).squaredNorm();
                if (ssx_sqrt_i_sub[j] < 1.0e-30) ssx_sqrt_i_sub[j] = 0.0;
                else ssx_sqrt_i_sub[j] = 1.0 / sqrt(ssx_sqrt_i_sub[j]);
            }
            
            if (size > size_limit) calcu_ld_blk_split(bdata,size, size_limit, X_sub, ssx_sqrt_i_sub, rsq_cutoff, rsq_size, mean_rsq_sub, max_rsq_sub, s1, s2, second);
            else {
                MatrixXf rsq_sub = X_sub.transpose() * X_sub;
                for (j = 0; j < size; j++) {
                    rsq_size[j] = 0.0;
                    mean_rsq_sub[j] = 0.0;
                    for (k = 0; k < size; k++) {
                        if (second) {
                            if (j <= s1 && k <= s1) continue;
                            if (j >= s2 && k >= s2) continue;
                        }
                        if (k == j) continue;
                        rsq_sub(j,k) *= (ssx_sqrt_i_sub[j] * ssx_sqrt_i_sub[k]);
                        rsq_sub(j,k) = rsq_sub(j,k) * rsq_sub(j,k);
                        if (rsq_sub(j,k) >= rsq_cutoff) {
                             mean_rsq_sub[j] += rsq_sub(j,k);
                            rsq_size[j] += 1.0;
                        }
                        if (rsq_sub(j,k) > max_rsq_sub[j]) max_rsq_sub[j] = rsq_sub(j,k);
                    }
                    if (rsq_size[j] > 0.0) mean_rsq_sub[j] /= rsq_size[j];
                }
            }
            
            for (j = 0, k = brk_pnt[i]; j < size; j++, k++) {
                if (second) {
                    if (rsq_size[j] > 0.0) {
                        mean_rsq[k] = (mean_rsq[k] * snp_num[k] + mean_rsq_sub[j] * rsq_size[j]) / (snp_num[k] + rsq_size[j]);
                        snp_num[k] = (snp_num[k] + rsq_size[j]);
                        if(max_rsq[k] < max_rsq_sub[j]) max_rsq[k] = max_rsq_sub[j];
                    }
                }
                else {
                    mean_rsq[k] = mean_rsq_sub[j];
                    snp_num[k] = rsq_size[j];
                    max_rsq[k] = max_rsq_sub[j];
                }
            }
        }
        cout<<ttlsize<<endl;
    }
    

    void calcu_mean_rsq(char* outFileName, char* bFileName,char* indilstName, char* indilst2remove,char* snplstName, char* snplst2exclde,int chr,double maf, bool ldr, bool ldr2, int ldWind, double rsq_cutoff)
    {
        bInfo bdata;
        MatrixXd X;
        vector<uint64_t> cols;
        vector<float> lds;
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        read_bimfile(&bdata, string(bFileName)+".bim");
        if(chr) extract_snp(&bdata, chr);
        if(snplstName != NULL) extract_snp(&bdata, snplstName);
        if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0) filter_snp_maf(&bdata, maf);
        check_autosome(&bdata);
        
        int i = 0, m = bdata._include.size();
        
        cout << "\nCalculating LD score for SNPs (block size of " << ldWind / 1000 << "Kb with an overlap of "<<ldWind/2000<<"Kb between blocks); LD rsq threshold = " << rsq_cutoff << ") ... " << endl;
        
        vector<int> brk_pnt1, brk_pnt2, brk_pnt3;
        get_ld_blk_pnt(&bdata,brk_pnt1, brk_pnt2, brk_pnt3, ldWind,0);
        
        VectorXf mean_rsq = VectorXf::Zero(m), snp_num = VectorXf::Zero(m), max_rsq = VectorXf::Zero(m);
        clock_t begin_time = clock();
        calcu_ld_blk(&bdata, brk_pnt1, brk_pnt3, mean_rsq, snp_num, max_rsq, false, rsq_cutoff);
        if (brk_pnt2.size() > 1) calcu_ld_blk(&bdata,brk_pnt2, brk_pnt3, mean_rsq, snp_num, max_rsq, true, rsq_cutoff);
        //printf(" cost2: %f ms.\n",float( clock () - begin_time ) /  1000);
        /*
        string mrsq_file = "";
        mrsq_file = _out + ".score.ld";
        ofstream o_mrsq(mrsq_file.data());
        o_mrsq<<"SNP chr bp MAF mean_rsq snp_num max_rsq ldscore"<<endl;
        double ldscore = 0.0;
        for (i = 0; i < m; i++){
            o_mrsq << bdata->_snp_name[bdata->_include[i]] << " " << bdata->_chr[bdata->_include[i]] << " " << bdata->_bp[bdata->_include[i]] << " ";
            double MAF = 0.5 * bdata->_mu[bdata->_include[i]];
            if(MAF > 0.5) MAF = 1.0 - MAF;
            ldscore = 1.0 + mean_rsq[i] * snp_num[i];
            o_mrsq << MAF << " " << mean_rsq[i] << " " << snp_num[i] << " " << max_rsq[i] << " " << ldscore << "\n";
        }
        o_mrsq << endl;
        cout << "LD score for " << m << " SNPs have been saved in the file [" + mrsq_file + "]." << endl;
         */
    }
    
}
