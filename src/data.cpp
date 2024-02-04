#include "data.hpp"


static double lower_bound = 1e-10;

FILE*
mgzip::open(const char* filename, const char* mode){

  //if(mode != "r" || mode != "rb"){
  //  std::cerr << "Mode " << mode << " is currently not supported" << std::endl;
  //}

  if(is_open == true){
    Rcpp::stop("Object already in use. Please define a new object for each file.");
  }

  //determine if file is gzipped
  FILE* fp_check = fopen(filename, "rb");
  unsigned char buffer;
  if(fp_check == NULL){
    fprintf(stderr, "Failed to open file %s", filename);
    Rcpp::stop("Error");
  }
  fread(&buffer, 1, 1, fp_check);
  is_gzipped = false;
  if(buffer == 0x1f){
    fread(&buffer, 1, 1, fp_check);
    if(buffer == 0x8b){
      fread(&buffer, 1, 1, fp_check);
      if(buffer == 0x08){
        is_gzipped = true;
      }
    }
  }
  fclose(fp_check);

  if(is_gzipped){

    const char prefix[] = "gunzip -c ";
    char *cmd;
    cmd = (char*) malloc(sizeof(prefix) + strlen(filename) + 10);
    if (!cmd) {
      fprintf(stderr, "%s: malloc: %s\n", filename, strerror(errno));
      Rcpp::stop("Error");
    }
    sprintf(cmd,"%s \'%s\'", prefix, filename);

    FILE* fp = popen(cmd, mode);
    if(fp == NULL){
      std::cerr << "Failed to open file " << filename << std::endl;
      Rcpp::stop("Error");
    }
    return(fp);

  }else{
  
    FILE* fp = fopen(filename, mode);
    if(fp == NULL){
      fprintf(stderr, "Failed to open file %s", filename);
      Rcpp::stop( "Error" );
    }
    return(fp);
  
  }

}

void
mgzip::close(FILE* fp){
  if(is_gzipped){
    pclose(fp);
  }else{
    fclose(fp);
  }
  is_open = false;
}

////////////////////////////
Data::Data(int N, int L, int Ne, double mu): N(N), L(L), Ne(Ne), mu(mu){};

////////////////////////////
void
haps::ReadSNP(std::vector<char>& sequence, int& bp){

  //read a line, extract bp and fp_props
  //snp;pos_of_snp;rs-id;ancestral_allele/alternative_allele;downstream_allele;upstream_allele;All;
  fscanf(fp, "%s %s %d %s %s", chr, rsid, &bp, ancestral, alternative);

  assert(sequence.size() > 0);
  //read haplotypes into sequence
  std::vector<char>::iterator it_seq = sequence.begin();  

  fgets(line, 2*N+10, fp);
  char d = line[0];
  int i  = 0;
  while(d != '\0' && it_seq != sequence.end()){
    if(d == '0'){
      *it_seq = '0';
      it_seq++;
    }else if(d == '1'){
      *it_seq = '1';
      it_seq++;
    }
    i++;
    d = line[i];
  }
  if(it_seq != sequence.end()){
    std::cerr << chr << " " << rsid << " " << bp << " " << ancestral << " " << alternative << std::endl;
    std::cerr << line << " " << i << std::endl;
    std::cerr << *it_seq << std::endl;
  }
  assert(it_seq == sequence.end());

}

void
haps::DumpSNP(std::vector<char>& sequence, int bp, FILE* fp_out){

  //read a line, extract bp and fp_props
  //snp;pos_of_snp;rs-id;ancestral_allele/alternative_allele;downstream_allele;upstream_allele;All;
  fprintf(fp_out, "%s %s %d %s %s", chr, rsid, bp, ancestral, alternative);

  //read haplotypes into sequence
  std::vector<char>::iterator it_seq = sequence.begin(); 
  for(; it_seq != sequence.end(); it_seq++){
    fprintf(fp_out, " %c", *it_seq);
  }
  fprintf(fp_out, "\n");

}

////////////////////////////
map::map(const char* filename){

  mgzip g;

  fp = g.open(filename, "r");
  assert(fp);
  int lines = 0;
  while(!feof(fp)){
    if(fgetc(fp) == '\n'){
      lines++;
    }
  }
  lines--;//don't count the header
  g.close(fp);

  fp = g.open(filename, "r");
  assert(fp);
  //skip header
  fscanf(fp, "%s", buffer); 
  fscanf(fp, "%s", buffer); 
  fscanf(fp, "%s", buffer); 

  bp.resize(lines);
	recrate.resize(lines);
  gen_pos.resize(lines);

  float dummy;
  double fbp;
  for(int snp = 0; snp < lines; snp++){
    fscanf(fp, "%lf %lf %lf", &fbp, &recrate[snp], &gen_pos[snp]);
    bp[snp] = fbp;
  }

  g.close(fp);

}

////////////////////////////
void
fasta::Read(const std::string filename){

  igzstream is(filename);
  if(is.fail()){
    is.open(filename + ".gz");
  }
  if(is.fail()){
    Rcpp::stop( "Error while opening file " + filename );
  }
  std::string line;
  getline(is,line);
  while(getline(is,line)){

    for(std::string::iterator it_c = line.begin(); it_c != line.end(); it_c++){
      *it_c = std::toupper(*it_c);
    }
    seq += line;
  }
  is.close();

}


