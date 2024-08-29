#ifndef EIGENSTRAT_HPP
#define EIGENSTRAT_HPP

#include <iostream>
#include <string>
#include <cassert>
#include <vector>
#include <map>
#include <math.h>
#include <RcppArmadillo.h>

#include "gzstream.hpp"

#if HAVE_ENDIAN_H
#include <endian.h>
#elif HAVE_MACHINE_ENDIAN_H
#include <machine/endian.h>
#elif HAVE_SYS_ENDIAN_H
#include <sys/endian.h>
#endif

inline bool exists (const std::string& name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }   
}

class eigenstrat{

  //class to read/write bed.
  //define with file pointer or filename and never use it to write and read.

  private:

    int N, L, num_groups;
    igzstream is_geno, is_snp; 
    int count_snps;

  public:

    std::map<std::string, int> groups; 
    std::vector<std::string> unique_groups; //unique groups
    std::vector<int> membership; //group membership of each individual

    eigenstrat(std::string pref);

    int ReadSNP(std::vector<int>& sequence, std::string& chr, int& bp, std::string& REF, std::string& ALT); //gets hap info for SNP

    int GetN(){return(N);}
    int GetL(){return(L);}

};

class plink{

  //class to read/write bed.
  //define with file pointer or filename and never use it to write and read.

  private:

#if __BYTE_ORDER == __LITTLE_ENDIAN
#include "snp_lookup_little.h"
#else
#include "snp_lookup_big.h"
#endif /* End test endianess */

    int N, L, maxind, num_groups, num_blocks;
    igzstream is_geno, is_snp; 
    FILE* fp_geno;
    int count_snps;
    unsigned char* blocks;

  public:

    std::map<std::string, int> groups; 
    std::vector<std::string> unique_groups; //unique groups
    std::vector<int> membership; //group membership of each individual

    plink(std::string pref);
    ~plink(){
      free(blocks);
    }

    void SetMaxInd(int max_ind){
      maxind = max_ind;
    }
    void ChangeFAM(std::vector<std::string>& fam);
    int ReadSNP(std::vector<int>& sequence, std::string& chr, int& bp, double& rec, std::string& REF, std::string& ALT); //gets hap info for SNP

    int GetN(){return(N);}
    int GetL(){return(L);}

};


#endif //EIGENSTRAT_HPP
