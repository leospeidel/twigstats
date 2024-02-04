#ifndef DATA_HPP
#define DATA_HPP

#include "collapsed_matrix.hpp"
#include "gzstream.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <tgmath.h>
#include <string>
#include <stdlib.h>
#include <string.h>
#include <cassert>
#include <RcppArmadillo.h>

class mgzip{

  private:

    bool is_gzipped = false;
    bool is_open    = false;

  public:

   mgzip(){};

   //function that takes in a file pointer, filename and opens the file
   FILE* open(const char* filename, const char* mode);
   
   //function that takes in a file pointer and closes the file.
   void close(FILE* fp);

};

//struct recording all the data needed for building anc
struct Data{

  int N, L; //number of sequences, number of SNPs
  int Ne; //effective population size
  double mu; //mutation rate

  std::vector<int> pos;    //vector specifying location of each SNP along the genome
  std::vector<double> r;   //vector of recombination distances from one SNP to the next
  std::vector<double> rpos; //vector of cumulative recombination distances

  ///////////

  Data(){}
  //Constructor, which only assigns values to N, L, Ne, mu
  Data(int N, int L, int Ne = 3e4, double mu = 1.25e-8);
  
};


class haps{

  //class to read/write bed.
  //define with file pointer or filename and never use it to write and read.

  private:

    int N, L;
    FILE* fp;
    mgzip g; 
    char* line;

  public:

    char chr[1024];
    char rsid[1024];
    char ancestral[1024], alternative[1024];

    haps(const char* filename_haps, const char* filename_sample){ 

      //fp = fopen(filename_sample, "r");
      fp = g.open(filename_sample, "r");
      assert(fp);
      N = 0;
      while(!feof(fp)){
        if(fgetc(fp) == '\n'){
          N++;
        }
      }
      N -= 2;
      N *= 2;
      g.close(fp);

      //fp = fopen(filename_haps, "r");
      fp = g.open(filename_haps, "r");
      assert(fp);
      L = 0;
      while(!feof(fp)){
        if(fgetc(fp) == '\n'){
          L++;
        }
      }
      g.close(fp);

      //fp = fopen(filename_haps, "r");
      fp = g.open(filename_haps, "r");
      assert(fp);
      line = (char*) malloc(2*N+10);

    }
    ~haps(){
      free(line);
    }

    void ReadSNP(std::vector<char>& sequence, int& bp); //gets hap info for SNP
    void DumpSNP(std::vector<char>& sequence, int bp, FILE* fp_out); //dumps hap info for SNP
    void CloseFile(){g.close(fp);};

    int GetN(){return(N);}
    int GetL(){return(L);}

};

class map{

  private:

    FILE* fp;
    char buffer[30];

  public:
 
    std::vector<int> bp;
		std::vector<double> recrate;
    std::vector<double> gen_pos;
    map(const char* filename);

};

struct fasta{

  std::string seq;
  void Read(const std::string filename);

};

#endif //DATA_HPP
