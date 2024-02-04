#include "eigenstrat.hpp"

////////////////////////////

eigenstrat::eigenstrat(std::string pref){

	std::string tmp;

	//read ind file
	igzstream is(pref + ".ind");
	if(is.fail()){
		is.open(pref + ".ind.gz");
	}
	if(is.fail()){
    Rcpp::stop("ERROR: unable to open ind file.");
	}
	N = 0;
	std::string id, rel, g;
	num_groups = 0;
	while(is >> id >> rel >> g){
		//check if group label exists	
		if(groups.count(g) == 0){
			groups[g] = num_groups;
			unique_groups.push_back(g);
			num_groups++;
		}
		membership.push_back(groups[g]);
		N++;
	}
	is.close();

	//now open geno and snp
	is_snp.open(pref + ".snp");
	if(is_snp.fail()){
		is_snp.open(pref + ".snp.gz");
	}
	if(is_snp.fail()){
    Rcpp::stop("ERROR: unable to open snp file.");
	}
	L = 0;
	while(getline(is_snp, tmp)){
		L++;
	}
	is_snp.close();
	is_snp.open(pref + ".snp");
	if(is_snp.fail()){
		is_snp.open(pref + ".snp.gz");
	}
	if(is_snp.fail()){
    Rcpp::stop("ERROR: unable to open snp file.");
	}

  count_snps = 0;
	is_geno.open(pref + ".geno");
	if(is_geno.fail()){
		is_geno.open(pref + ".geno.gz");
	}
	if(is_geno.fail()){
    Rcpp::stop("ERROR: unable to open geno file.");
	}

}

int
eigenstrat::ReadSNP(std::vector<int>& sequence, std::string& chr, int& bp, std::string& REF, std::string& ALT){

	if(sequence.size() < (unsigned int) N) sequence.resize(N);

  count_snps++;

	std::string tmp;
  if(count_snps <= L){
    is_snp >> tmp >> chr >> tmp >> bp >> REF >> ALT;
	  getline(is_geno, tmp);

		std::string::iterator it_tmp = tmp.begin();
		for(std::vector<int>::iterator it_seq = sequence.begin(); it_seq != sequence.end();){
      (*it_seq) = (*it_tmp) - '0';
			it_tmp++;
			it_seq++;
			if(0){
      if(it_tmp == tmp.end()){
        if(it_seq != sequence.end()){
          Rcpp::stop("ERROR while reading geno file.");
				}
			}
      }
		}

	  return 1;
	}else{
    return 0;
	}

}

////////////////////////////

//https://www.cog-genomics.org/plink2/formats#bed

plink::plink(std::string pref){

  std::string tmp;

  //read ind file
  igzstream is(pref + ".fam");
  if(is.fail()){
    is.open(pref + ".fam.gz");
  }
  if(is.fail()){
    Rcpp::stop("ERROR: unable to open ind file.");
  }
  N = 0;
  std::string id, rel, g;
  num_groups = 0;
  while(is >> g >> id >> rel >> rel >> rel >> rel){
    //check if group label exists	
    if(groups.count(g) == 0){
      groups[g] = num_groups;
      unique_groups.push_back(g);
      num_groups++;
    }
    membership.push_back(groups[g]);
    N++;
  }
  is.close();
  num_blocks = std::ceil(N/4.0);

  //now open geno and snp
  is_snp.open(pref + ".bim");
  if(is_snp.fail()){
    is_snp.open(pref + ".bim.gz");
  }
  if(is_snp.fail()){
    Rcpp::stop("ERROR: unable to open snp file.");
  }
  L = 0;
  while(getline(is_snp, tmp)){
    L++;
  }
  is_snp.close();
  is_snp.open(pref + ".bim");
  if(is_snp.fail()){
    is_snp.open(pref + ".bim.gz");
  }
  if(is_snp.fail()){
    Rcpp::stop("ERROR: unable to open snp file.");
  }

  count_snps = 0;
  fp_geno = fopen( (pref + ".bed").c_str() , "rb");
  if(fp_geno == NULL){
    Rcpp::stop("ERROR: unable to open bed file.");
  }
  unsigned char header[3];
  if(fread( header, 1, 3, fp_geno) != 3){
    Rcpp::stop("Failed to parse bed file. Possibly empty.");
  }else{
    if(header[0] != 0x6c || header[1] != 0x1b ||  header[2] != 0x01){
      Rcpp::warning("Magic numbers in bed not matching with expectation.");
    }
  }
  blocks = (unsigned char*) malloc(num_blocks);

}


void
plink::ChangeFAM(std::vector<std::string>& fam){
  int N2 = 0;
  groups.clear();
  unique_groups.clear();
  membership.clear();
  num_groups = 0;
  for(int i = 0; i < fam.size(); i++){
    if(groups.count(fam[i]) == 0){
      groups[fam[i]] = num_groups;
      unique_groups.push_back(fam[i]);
      num_groups++;
    }
    membership.push_back(groups[fam[i]]);
    N2++;
  }
  assert(N == N2);
}

int
plink::ReadSNP(std::vector<int>& sequence, std::string& chr, int& bp, double& rec, std::string& REF, std::string& ALT){

  if(sequence.size() < (unsigned int) maxind) sequence.resize(maxind);

  count_snps++;

  std::string tmp;
  if(count_snps <= L){

    //chrID variantID cM bp allele1 allele2
    is_snp >> chr >> tmp >> rec >> bp >> REF >> ALT;

    //parse bed file
    fread(blocks,1,num_blocks,fp_geno);

    //0,1,2 genotypes and 9 means missing
    int j = 0;
    std::vector<int>::iterator it_seq = sequence.begin();
    int* snp_block = (int*) malloc(4*sizeof(int));
    for(int i = 0; i < num_blocks; i++){
      snp_block = snp_lookup[ (int)blocks[i] ];
      //std::cerr << (int) blocks[i] << " ";
      int k = 0;
      int j_start = j;
      for(; j < std::min(j_start+4,maxind); j++){
        *it_seq = snp_block[k];
        //*it_seq = snp_lookup[ (int)blocks[i] ][k];
        k++;
        it_seq++;
      }
    }
    //std::cerr << std::endl;

    return 1;
  }else{
    return 0;
  }

}

