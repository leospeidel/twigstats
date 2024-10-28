#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <numeric>
#include <iomanip>

#include "gzstream.hpp"
#include "sample.hpp"
#include "eigenstrat.hpp"
#include "mutations.hpp"

using namespace Rcpp;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//' Chromosome painting using genealogies.
//'
//' This function outputs the first coalescence with an individual from a pre-specified group identity along the genome.
//' If the first such coalescnece involves several copying candidates, a random haplotype is chosen.
//' Output is in GLOBEtrotter format.
//'
//' @param file_anc Filename of anc file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).
//' @param file_mut Filename of mut file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).
//' @param file_map File prefix of recombination map.
//' @param file_out File prefix of output files
//' @param poplabels Filename of poplabels file
//' @param blgsize (Optional) SNP block size in Morgan. Default is 0.05 (5 cM). If blgsize is 1 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
//' @param pops (Optional) Populations for which data should be extracted. Names need to match the second column of the poplabels file
//' @param chrs (Optional) Vector of chromosome IDs
//' @return void. Write three files idfile, paint, rec to disc. 
//' @examples
//' file_anc  <- system.file("sim/msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats")
//' file_mut  <- system.file("sim/msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats")
//' poplabels <- system.file("sim/msprime_ad0.8_split250_1.poplabels", package = "twigstats")
//' file_map  <- system.file("sim/genetic_map_combined_b37_chr1.txt.gz", package = "twigstats")
//'
//' #define populations to paint against:
//' pops <- c("P1","P2","P3","P4")
//'
//' Painting(file_anc, file_mut, file_map, file_out = "test", poplabels, blgsize = 1e-5)
//' @export
// [[Rcpp::export]]
void Painting( SEXP file_anc, SEXP file_mut, SEXP file_map, SEXP file_out, SEXP poplabels, Nullable<double> blgsize = R_NilValue, Nullable<StringVector> pops = R_NilValue, Nullable<CharacterVector> chrs = R_NilValue) {

	double binsize_ = 0.05;
	if(blgsize.isNotNull()){
		binsize_ = as<double>(blgsize);
	}
	bool isM = true;
	if(binsize_ >= 1) isM = false;
	if(isM) binsize_ *= 100; //convert to cM

	std::string filename_poplabels = as<std::string>(poplabels);
	std::vector<std::string> filename_anc, filename_mut, filename_rec, chrname, filename_recout, filename_paintout;
	std::string map_ = as<std::string>(file_map);
	std::string filename_out = as<std::string>(file_out);

	if(chrs.isNotNull()){
		CharacterVector chrs_(chrs);
		for(int i = 0; i < chrs_.size(); i++){
			filename_anc.push_back( as<std::string>(file_anc) + "_chr" + as<std::string>(chrs_[i]) + ".anc" );
			filename_mut.push_back( as<std::string>(file_mut) + "_chr" + as<std::string>(chrs_[i]) + ".mut" );
			filename_rec.push_back( map_ + "_chr" + as<std::string>(chrs_[i]) + ".txt" );
			chrname.push_back(as<std::string>(chrs_[i]));
			filename_recout.push_back( filename_out + "_rec_chr" + chrname[i] + ".txt" );
			filename_paintout.push_back( filename_out + "_painting_chr" + chrname[i] + ".txt" );
		}
	}else{
		filename_anc.push_back( as<std::string>(file_anc) );
		filename_mut.push_back( as<std::string>(file_mut) );
		filename_rec.push_back( map_ );
		filename_recout.push_back( filename_out + "_rec.txt" );
		filename_paintout.push_back( filename_out + "_painting.txt" );
	}

	Sample sample;
	sample.Read(filename_poplabels);

	StringVector pops_(sample.groups.size());
	pops_ = sample.groups;
	if(pops.isNotNull()){
		pops_ = as<StringVector>(pops);
	}
	std::vector<int> group_used(sample.groups.size());
	for(int i = 0; i < pops_.size(); i++){
		bool exists = false;
		for(int j = 0; j < sample.groups.size(); j++){
			if( pops_[i] == sample.groups[j] ){
				group_used[j] = 1;
				exists = true;
				break;
			}
		}
		if(!exists){
			Rcpp::stop("Group does not exist");
		}
	}


	std::cerr << "Painting genomes using first coalescence among reference groups for " << sample.groups.size() << " groups." << std::endl;
	std::cerr << std::endl;

	MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
	Muts::iterator it_mut; //iterator for mut file
	float num_bases_tree_persists = 0.0;

	std::mt19937 rng;
	std::uniform_real_distribution<double> runif(0.0f, 1.0f);

	//For each individual, decide who they are painted by out of individuals in pops
	//Then output this per individual
	//I want to output this every 10kb or so.

	int N = sample.group_of_haplotype.size();

	bool is_haploid = false;
	if(N == sample.ind.size()){
		is_haploid = true;
	}else{
		assert(2*sample.ind.size() == N);
	}

	////////// 1. Read one tree at a time /////////

	double total = 0;
	double rec = 0.0, genpos = 0.0;
	int  irec = 0;
	for(int chr = 0; chr < filename_anc.size(); chr++){

		int percentage = 0;

		std::vector<std::vector<int>> painting_profile(N);
		std::ofstream os_rec(filename_recout[chr]);
		os_rec << "start.pos recom.rate.perbp\n"; 

		//We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
		//The mut file is read once, file is closed after constructor is called.
		AncMutIterators ancmut(filename_anc[chr], filename_mut[chr]);
		Data data(ancmut.NumTips(), ancmut.NumSnps());
		map recmap(filename_rec[chr].c_str());
		rec = 0.0;
		genpos = 0.0;
		irec = 0;

		//mtr stores the marginal tree, it_mut points to the first SNP at which the marginal tree starts.
		//If marginal tree has no SNPs mapping to it, it_mut points to the first SNP of the following tree that has a SNP mapping to it
		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

		std::vector<float> coords(2*data.N-1, 0.0);
		std::vector<Leaves> desc;


		int current_ipos = 0, current_igenpos = 0;
		int current_pos = (*it_mut).pos;
		double current_genpos = 0;
		//get genpos for first tree
		if(recmap.bp[irec] >= (*it_mut).pos){
			if(irec > 0){
				genpos = ((double)(*it_mut).pos - recmap.bp[irec-1])/(recmap.bp[irec] - recmap.bp[irec-1]) * (recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) + recmap.gen_pos[irec-1];
			}else{
				genpos = ((double)(*it_mut).pos)/(recmap.bp[irec]) * (recmap.gen_pos[irec]);         
			}
		}else{
			while( recmap.bp[irec] < (*it_mut).pos ){
				irec++;
				if(irec == recmap.bp.size()){
					irec--;
					break;
				}
			}
			if(irec == recmap.bp.size()-1){
				genpos = recmap.gen_pos[irec]/recmap.bp[irec] * ((*it_mut).pos - recmap.gen_pos[irec]) + recmap.gen_pos[irec];
			}else{
				assert(recmap.bp[irec] >= (*it_mut).pos);
				assert(irec > 0);
				genpos = ((double)(*it_mut).pos - recmap.bp[irec-1])/(recmap.bp[irec] - recmap.bp[irec-1]) * (recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) + recmap.gen_pos[irec-1];
			}
		}
		current_genpos = genpos;

		//
		if(isM){
			current_igenpos = std::ceil(current_genpos/binsize_)-1;
			//current_genpos  = current_igenpos * binsize_;
		}else{
			current_ipos = std::ceil(current_pos/binsize_)-1;
			//current_pos  = current_igenpos * binsize_;
		}

		bool first = true;
		int next_pos;

		//iterate through whole file
		while(num_bases_tree_persists >= 0.0){

			if(((*it_mut).tree % (int)(ancmut.NumTrees()/100)) == 0){
				std::cerr << "[" << percentage << "%]\r";
				Rcpp::checkUserInterrupt();
				percentage++;
				std::cerr.flush();
			}

			//get bp position of next tree
			int treeID = (*it_mut).tree;
			while(treeID == (*it_mut).tree){
				next_pos = (*it_mut).pos;
				it_mut++;
				if(it_mut == ancmut.mut_end()) break;
			}
			//std::cerr << "next_pos: " << next_pos << std::endl;

			//get genpos of next tree
			int irec_orig = irec;
			if(recmap.bp[irec] >= next_pos){
				if(irec > 0){
					genpos = ((double)next_pos - recmap.bp[irec-1])/(recmap.bp[irec] - recmap.bp[irec-1]) * (recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) + recmap.gen_pos[irec-1];
				}else{
					genpos = ((double)next_pos)/(recmap.bp[irec]) * (recmap.gen_pos[irec]);         
				}
			}else{
				while( recmap.bp[irec] < next_pos ){
					irec++;
					if(irec == recmap.bp.size()){
						irec--;
						break;
					}
				}
				if(irec == recmap.bp.size()-1){
					genpos = recmap.gen_pos[irec]/recmap.bp[irec] * (next_pos - recmap.bp[irec]) + recmap.gen_pos[irec];
				}else{
					assert(recmap.bp[irec] >= next_pos);
					assert(irec > 0);
					genpos = ((double)next_pos - recmap.bp[irec-1])/(recmap.bp[irec] - recmap.bp[irec-1]) * (recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) + recmap.gen_pos[irec-1];
				}
			}

			//assume that current_pos and current_genpos are coordinates of last bin
			//output painting every x bases or x cM
			if( (!isM && next_pos - current_pos > binsize_) || (isM && genpos - current_genpos > binsize_) ){

				//get all pos/gpos pairs until next tree
				std::vector<double> pos, gpos;
				irec = irec_orig;
				if(!isM){
					assert((current_ipos+1)*binsize_ >= current_pos);
					assert((current_ipos+1)*binsize_ <= next_pos);
					for(int bp = (current_ipos+1)*binsize_; bp <= next_pos; bp += binsize_){
						if(recmap.bp[irec] >= bp){
							if(irec > 0){
								gpos.push_back(((double)bp - recmap.bp[irec-1])/(recmap.bp[irec] - recmap.bp[irec-1]) * (recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) + recmap.gen_pos[irec-1]);
							}else{
								gpos.push_back(((double)bp)/(recmap.bp[irec]) * (recmap.gen_pos[irec]));         
							}
						}else{
							while( recmap.bp[irec] < bp ){
								irec++;
								if(irec == recmap.bp.size()){
									irec--;
									break;
								}
							}
							if(irec == recmap.bp.size()-1){
								gpos.push_back(recmap.gen_pos[irec]/recmap.bp[irec] * (bp - recmap.bp[irec]) + recmap.gen_pos[irec]);
							}else{
								assert(recmap.bp[irec] >= bp);
								assert(irec > 0);
								gpos.push_back(((double)bp - recmap.bp[irec-1])/(recmap.bp[irec] - recmap.bp[irec-1]) * (recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) + recmap.gen_pos[irec-1]);
							}
						}
						pos.push_back((int)bp);
					}
				}else{
					assert((current_igenpos+1)*binsize_ >= current_genpos);
					assert((current_igenpos+1)*binsize_ <= genpos);
					for(double gbp = (current_igenpos + 1.0)*binsize_; gbp <= genpos; gbp += binsize_){
						if(recmap.gen_pos[irec] >= gbp){
							assert(recmap.gen_pos[irec-1] <= gbp);
							if(irec > 0){
								pos.push_back((int)(((double)gbp - recmap.gen_pos[irec-1])/(recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) * (recmap.bp[irec] - recmap.bp[irec-1]) + recmap.bp[irec-1]));
							}else{
								pos.push_back((int)(((double)gbp)/(recmap.gen_pos[irec]) * (recmap.bp[irec])));         
							}
						}else{
							while( recmap.gen_pos[irec] < gbp ){
								irec++;
								if(irec == recmap.bp.size()){
									irec--;
									break;
								}
							}
							if(irec == recmap.bp.size()-1){
								pos.push_back((int)(recmap.bp[irec]/recmap.gen_pos[irec] * (gbp - recmap.gen_pos[irec]) + recmap.bp[irec]));
							}else{
								assert(recmap.gen_pos[irec] >= gbp);
								assert(recmap.gen_pos[irec-1] <= gbp);
								assert(irec > 0);
								pos.push_back((int)(((double)gbp - recmap.gen_pos[irec-1])/(recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) * (recmap.bp[irec] - recmap.bp[irec-1]) + ((double)recmap.bp[irec-1]) ));
							}
						}
						gpos.push_back(gbp);
					}
				}
				assert(pos.size() > 0);
				assert(gpos.size() > 0);

				rec = (gpos[0] - current_genpos)/(pos[0] - current_pos);
				assert(rec >= 0.0);
				if(!first) os_rec << current_pos << " " << rec/1e2 << "\n";
				first = false;

				for(int j = 1; j < pos.size(); j++){
					rec = (gpos[j] - gpos[j-1])/(pos[j] - pos[j-1]);
					assert(rec >= 0.0);
					os_rec << (int)pos[j-1] << " " << rec/1e2 << "\n";
				}
				current_ipos += pos.size();
				current_pos = pos[pos.size()-1];
				current_igenpos += gpos.size();
				current_genpos = gpos[gpos.size()-1];

				mtr.tree.GetCoordinates(coords);
				mtr.tree.FindAllLeaves(desc);
				total += num_bases_tree_persists;

				//Find nearest neighbour, get their identity, and add to corresponding row
				for(int i = 0; i < data.N; i++){
					int node = i;
					int parent = (*mtr.tree.nodes[node].parent).label;
					int sib    = (*mtr.tree.nodes[parent].child_left).label;
					int group = sample.group_of_haplotype[i];
					int group2;
					if(node == sib) sib = (*mtr.tree.nodes[parent].child_right).label;
					int nonself;

					//naive implementation of nearest neighbour 
					while(mtr.tree.nodes[node].parent != NULL){
						parent = (*mtr.tree.nodes[node].parent).label;
						sib    = (*mtr.tree.nodes[parent].child_left).label;
						if(node == sib) sib = (*mtr.tree.nodes[parent].child_right).label;
						nonself = 0;
						for(int k = 0; k < desc[sib].member.size(); k++){
							if(group_used[ sample.group_of_haplotype[desc[sib].member[k]] ] != 0) nonself++;
						}
						/*
							 if(pos[0] == 509747 && i == 0){
							 std::cerr << i << " " << nonself << " " << desc[sib].member.size() << std::endl;
							 for(int k = 0; k < group_used.size(); k++){
							 std::cerr << group_used[k] << " ";
							 }
							 std::cerr << std::endl;
							 for(int k = 0; k < desc[sib].member.size(); k++){
							 std::cerr << desc[sib].member[k] << " of group " << sample.group_of_haplotype[desc[sib].member[k]] << std::endl;
							 }
							 }
							 */
						//draw a random number between 0 and nonself
						int copyind = (int)nonself*runif(rng);

						//if(pos[0] == 500152) std::cerr << copyind << std::endl;

						if(nonself > 0){
							int count = 0;
							for(int k = 0; k < desc[sib].member.size(); k++){
								group2 = sample.group_of_haplotype[desc[sib].member[k]];    
								if( group_used[group2] != 0 ){
									if(count == copyind){
										//if(pos[0] == 500152) std::cerr << "ind i = " << i << " has " << desc[sib].member[k] << " of group " << group2 << std::endl;
										for(int j = 0; j < pos.size(); j++){
											painting_profile[i].push_back(desc[sib].member[k]+1);
										}
										break;
									}
									count++;
								}
							}
							break;
						}
						node = parent;
					}

				}

			}
			num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
		}

		rec = 0.0;
		os_rec << current_pos << " " << rec/1e2 << "\n";
		current_pos = (*it_mut).pos;
		os_rec.close();

		//output to file
		std::ofstream os(filename_paintout[chr]);
		os << "EM_iter = 0 (N_e = 0 / copy_prop = 0 / mutation = 0 / mutationGLOBAL = 0), nsamples = 1, N_e_start = 400000.000000 (divided by number of donor haplotypes), region_size = 100.000000\n";

		if(is_haploid){
			for(int i = 0; i < sample.ind.size(); i++){
				os << "HAP 1 " << sample.ind[i] << "\n";
				os << "1 ";
				for(int j = 0; j < painting_profile[i].size(); j++){
					os << 2*painting_profile[i][j] << " ";
				}
				os << "\n";

				os << "HAP 2 " << sample.ind[i] << "\n";
				os << "1 ";
				for(int j = 0; j < painting_profile[i].size(); j++){
					os << 2*painting_profile[i][j] << " ";
				}
				os << "\n";
			}
			os.close();
		}else{
			for(int i = 0; i < sample.ind.size(); i++){
				os << "HAP 1 " << sample.ind[i] << "\n";
				os << "1 ";
				for(int j = 0; j < painting_profile[2*i].size(); j++){
					os << painting_profile[2*i][j] << " ";
				}
				os << "\n";

				os << "HAP 2 " << sample.ind[i] << "\n";
				os << "1 ";
				for(int j = 0; j < painting_profile[2*i+1].size(); j++){
					os << painting_profile[2*i+1][j] << " ";
				}
				os << "\n";
			}
			os.close();
		}

	}

	if(is_haploid){
		std::ofstream os_id(filename_out + "_idfile.txt");
		for(int i = 0; i < sample.ind.size(); i++){
			os_id << sample.ind[i] << " " << sample.groups[sample.group_of_haplotype[i]] << " 1\n";
		}
		os_id.close();
	}else{
		std::ofstream os_id(filename_out + "_idfile.txt");
		for(int i = 0; i < sample.ind.size(); i++){
			os_id << sample.ind[i] << " " << sample.groups[sample.group_of_haplotype[2*i]] << " 1\n";
		}
		os_id.close();
	}

	std::cerr << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Helper function to split a string by space and convert to a vector of doubles
std::vector<double> split_string_to_doubles(const std::string &str) {
	std::stringstream ss(str);
	std::vector<double> result;
	double num;
	while (ss >> num) {
		result.push_back(num-1);
	}
	return result;
}

// Function to read ID file
std::vector<std::pair<std::string, std::string>> read_ids(const std::string &filename) {
	igzstream file(filename);
	std::vector<std::pair<std::string, std::string>> ids;
	std::string line;
	while (std::getline(file, line)) {
		std::stringstream ss(line);
		std::string id;
		std::string population;
		ss >> id >> population;
		ids.emplace_back(id, population);
		ids.emplace_back(id, population);
	}
	return ids;
}

// Function to print progress bar
void print_progress_bar(int iteration, int total, int bar_length = 50) {
	int progress = static_cast<int>(static_cast<double>(iteration) / total * bar_length);
	int percent = static_cast<int>(static_cast<double>(iteration) / total * 100);
	Rcout << "\r[";
	for (int i = 0; i < progress; ++i) Rcout << "=";
	for (int i = progress; i < bar_length; ++i) Rcout << " ";
	Rcout << "] " << percent << "%" << std::flush;
}

// Function to perform bootstrap sampling with replacement
std::vector<std::map<std::string, std::map<std::string, double>>> bootstrap_sample(
		const std::vector<std::map<std::string, std::map<std::string, double>>>& df,
		int nboot) {

	std::vector<std::map<std::string, std::map<std::string, double>>> bootstrapped;
	int seed = 1;
	std::mt19937 rng(seed);
	std::uniform_int_distribution<int> distribution(0, df.size() - 1);

	if(nboot > 1){

		// Perform 'nboot' iterations
		for (int i = 0; i < nboot; ++i) {
			std::map<std::string, std::map<std::string, double>> bootstrapped_map;
			std::map<std::string, double> bootstrapped_map_first;

			std::vector<int> selected_indices(df.size());
			// Sample df.size() entries from df with replacement
			for (size_t j = 0; j < df.size(); ++j) {
				selected_indices[j] = distribution(rng);
			}
			for(const auto& index : selected_indices){
				const auto& selected_map = df[index];
				// Count number of blocks per target
				for (const auto& [first_str, inner_map] : selected_map) {
					for (const auto& [second_str, value] : inner_map) {
						bootstrapped_map_first[first_str] += value;
					}
				}
			}
			for(const auto& index : selected_indices){
				const auto& selected_map = df[index];
				// Count number of blocks per target and reference population
				for (const auto& [first_str, inner_map] : selected_map) {
					for (const auto& [second_str, value] : inner_map) {
						bootstrapped_map[first_str][second_str] += value/bootstrapped_map_first[first_str]; 
					}
				}
			}
			bootstrapped.push_back(bootstrapped_map);
		}

	}else{

		for (int i = 0; i < nboot; ++i) {
			std::map<std::string, std::map<std::string, double>> bootstrapped_map;
			std::map<std::string, double> bootstrapped_map_first;

			std::vector<int> selected_indices(df.size());
			// Sample df.size() entries from df with replacement
			for (size_t j = 0; j < df.size(); ++j) {
				selected_indices[j] = j;
			}
			for(const auto& index : selected_indices){
				const auto& selected_map = df[index];
				// Count number of blocks per target
				for (const auto& [first_str, inner_map] : selected_map) {
					for (const auto& [second_str, value] : inner_map) {
						bootstrapped_map_first[first_str] += value;
					}
				}
			}
			for(const auto& index : selected_indices){
				const auto& selected_map = df[index];
				// Count number of blocks per target and reference population
				for (const auto& [first_str, inner_map] : selected_map) {
					for (const auto& [second_str, value] : inner_map) {
						bootstrapped_map[first_str][second_str] += value/bootstrapped_map_first[first_str]; // Sum values for identical pairs
					}
				}
			}
			bootstrapped.push_back(bootstrapped_map);
		}
	}

	return bootstrapped;
}

//' Function to compute genome-wide copying proportions.
//'
//' This function takes the output of the function Painting and computes the genome-wide 'copying vectors',
//' i.e., the proportion of the genome copied from each reference population.
//'
//' @param filename_painting A vector containing filenames of painting profiles. These are the outputs of Painting.
//' @param filename_idfile The filename of the idfile, which is an output of Painting.
//' @param nboot Number of bootstrap samples to generate.
//' @param blocksize Number of blocks to combine for the bootstrap. For example, if Painting was run with a blgsize of 1e-5 Morgans, blocksize should be 5000 to achieve a blocksize of 5cM.
//' @param use_IDs If TRUE, compute profiles for each sample. If FALSE (default), compute profiles for each group as specified in the second column of the idfile.
//' @return A DataFrame with the copying proportions per bootstrap sample.
//' @examples
//' # Example usage in R:
//' path      <- paste0(system.file("test/", package = "twigstats"),"/")
//' prefix    <- "test" # prefix of files under path
//' df <- PaintingProfile(c(paste0(path, prefix, "_painting.txt.gz")),
//'                          paste0(path, prefix, "_idfile.txt.gz"),
//'                          nboot = 10, blocksize = 5000)
//' head(df)
//' @export
// [[Rcpp::export]]
DataFrame PaintingProfile(std::vector<std::string> filename_painting, std::string filename_idfile, 
		int nboot, int blocksize, bool use_IDs = false) {
	if (nboot <= 0) {
		stop("Value of nboot needs to be at least 1.");
	}
	if (blocksize <= 0) {
		stop("Value of blocksize needs to be at least 1.");
	}

	// Read the ID file
	std::vector<std::pair<std::string, std::string>> ids = read_ids(filename_idfile);
	std::vector<std::map<std::string, std::map<std::string, double>>> df;
	int chr = 1;

	int block_id = 0;
	// Process each painting file
	for (const auto &file : filename_painting) {
		igzstream infile(file);
		std::string line;
		int indiv  = 0;
		std::string pop;

		// Skip the first two lines (header)
		std::getline(infile, line);
		//std::getline(infile, line);

		int block_id_tmp;
		while (std::getline(infile, line)) {

			std::getline(infile, line);
			std::vector<double> fields = split_string_to_doubles(line);

			if(use_IDs){
				pop = ids[indiv].first;
			}else{
				pop = ids[indiv].second;
			}

			if(block_id + (fields.size() / blocksize) + 10 > df.size()){
				df.resize( df.size() + (fields.size() / blocksize) + 10 );
			}

			block_id_tmp = block_id;
			for (int i = 1; i < fields.size(); i++) { // Skip first column

				if(i % blocksize == 0){
					block_id_tmp++;
				}
				std::string lab = ids[static_cast<int>(fields[i])].second;
				df[block_id_tmp][pop][lab]++;

			}

			indiv++;
		}

		block_id = block_id_tmp + 1;
		print_progress_bar(chr, filename_painting.size());
		chr++;
	}
	Rcout << std::endl;
  df.resize(block_id);

	//std::cerr << df.size() << std::endl;
	//std::cerr << df[0]["P1"]["P1"] << " " << df[1]["P1"]["P1"] << " " << df[2]["P1"]["P1"] << " " <<  df[3]["P1"]["P1"] << std::endl;

	//for(int k = 0; k < df.size(); k++){
	//  std::cerr << std::get<0>(df[k]) << " " << std::get<1>(df[k]) << " " << std::get<2>(df[k]) << " " << std::get<3>(df[k]) << std::endl;
	//}

	// Perform bootstrapping and summarize results
	auto summary = bootstrap_sample(df, nboot);

	// Prepare the result as a DataFrame
	std::vector<std::string> pop_vec;
	std::vector<std::string> lab_vec;
	std::vector<double> prop_vec;
	std::vector<int> bootstrap_vec;

	for(int b = 0; b < nboot; b++){
		for (const auto &entry : summary[b]) {
			const std::string& pop = entry.first;
			for (const auto &lab_entry : entry.second) {
				pop_vec.push_back(pop);
				lab_vec.push_back(lab_entry.first);
				prop_vec.push_back(lab_entry.second);
				bootstrap_vec.push_back(b); // Assuming all are part of the bootstrap
			}
		}
	}

	return DataFrame::create(
			Named("POP") = pop_vec,
			Named("labs") = lab_vec,
			Named("prop") = prop_vec,
			Named("bootstrap") = bootstrap_vec
			);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO: Implement genetic distance
//TODO: Implement bootstrap

//' Function to calculate first coalescence copying vector for each population specified in poplabels file.
//'
//' This function calculates the proportion of the genome where the first coalescence is with group x.
//' It returns a version where coalescence with own group is counted, and one where the own group is excluded.
//'
//' @param file_anc Filename of anc file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).
//' @param file_mut Filename of mut file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).
//' @param poplabels Filename of poplabels file
//' @param pops (Optional) Populations for which data should be extracted. Names need to match the second column of the poplabels file
//' @param chrs (Optional) Vector of chromosome IDs
//' @return 2d array of dimension #groups x #groups.
//' @keywords internal
//' @examples
//' file_anc  <- system.file("sim/msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats")
//' file_mut  <- system.file("sim/msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats")
//' poplabels <- system.file("sim/msprime_ad0.8_split250_1.poplabels", package = "twigstats")
//'
//' #Calculate f2s between all pairs of populations
//' nn_res <- ExpPaintingProfile(file_anc, file_mut, poplabels)
//' if(require(lsei)){
//'   pnnls(t(nn_res[c("P2","P3"),c("P1","P2","P3","P4"),1]),(nn_res["PX",c("P1","P2","P3","P4"),1]), sum = 1)$x
//' }
//'
//' @export
// [[Rcpp::export]]
NumericVector ExpPaintingProfile( SEXP file_anc, SEXP file_mut, SEXP poplabels, Nullable<StringVector> pops = R_NilValue, Nullable<CharacterVector> chrs = R_NilValue) {

	std::string filename_poplabels = as<std::string>(poplabels);
	std::vector<std::string> filename_anc, filename_mut;

	if(chrs.isNotNull()){
		CharacterVector chrs_(chrs);
		for(int i = 0; i < chrs_.size(); i++){
			filename_anc.push_back( as<std::string>(file_anc) + "_chr" + as<std::string>(chrs_[i]) + ".anc" );
			filename_mut.push_back( as<std::string>(file_mut) + "_chr" + as<std::string>(chrs_[i]) + ".mut" );
		}
	}else{
		filename_anc.push_back( as<std::string>(file_anc) );
		filename_mut.push_back( as<std::string>(file_mut) );
	}

	Sample sample;
	sample.Read(filename_poplabels);

	StringVector pops_(sample.groups.size());
	pops_ = sample.groups;
	if(pops.isNotNull()){
		pops_ = as<StringVector>(pops);
	}
	std::vector<int> group_used(sample.groups.size());
	for(int i = 0; i < pops_.size(); i++){
		bool exists = false;
		for(int j = 0; j < sample.groups.size(); j++){
			if( pops_[i] == sample.groups[j] ){
				group_used[j] = 1;
				exists = true;
				break;
			}
		}
		if(!exists){
			Rcpp::stop("Group does not exist");
		}
	}


	std::cerr << "Computing nearest neighbours for " << sample.groups.size() << " groups." << std::endl;
	std::cerr << std::endl;

	MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
	Muts::iterator it_mut; //iterator for mut file
	float num_bases_tree_persists = 0.0;

	////////// 1. Read one tree at a time /////////

	double total = 0;
	arma::cube nn_res = arma::zeros<arma::cube>(sample.groups.size(), sample.groups.size(), 2);
	for(int chr = 0; chr < filename_anc.size(); chr++){

		//We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
		//The mut file is read once, file is closed after constructor is called.
		AncMutIterators ancmut(filename_anc[chr], filename_mut[chr]);
		Data data(ancmut.NumTips(), ancmut.NumSnps());

		//mtr stores the marginal tree, it_mut points to the first SNP at which the marginal tree starts.
		//If marginal tree has no SNPs mapping to it, it_mut points to the first SNP of the following tree that has a SNP mapping to it
		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

		std::vector<float> coords(2*data.N-1, 0.0);
		std::vector<Leaves> desc;
		int current_pos = 0, percentage = 0;
		//iterate through whole file
		while(num_bases_tree_persists >= 0.0){

			if(((*it_mut).tree % (int)(ancmut.NumTrees()/100)) == 0){
				std::cerr << "[" << percentage << "%]\r";
				Rcpp::checkUserInterrupt();
				percentage++;
				std::cerr.flush();
			}

			mtr.tree.GetCoordinates(coords);
			mtr.tree.FindAllLeaves(desc);
			total += num_bases_tree_persists;

			//Find nearest neighbour, get their identity, and add to corresponding row
			for(int i = 0; i < data.N; i++){
				int node = i;
				int parent = (*mtr.tree.nodes[node].parent).label;
				int sib    = (*mtr.tree.nodes[parent].child_left).label;
				int group = sample.group_of_haplotype[i];
				int group2;
				if(node == sib) sib = (*mtr.tree.nodes[parent].child_right).label;
				int nonself;

				//naive implementation of nearest neighbour
				while(mtr.tree.nodes[node].parent != NULL){
					parent = (*mtr.tree.nodes[node].parent).label;
					sib    = (*mtr.tree.nodes[parent].child_left).label;
					if(node == sib) sib = (*mtr.tree.nodes[parent].child_right).label;
					nonself = 0;
					for(int k = 0; k < desc[sib].member.size(); k++){
						if(group_used[ sample.group_of_haplotype[desc[sib].member[k]] ] != 0) nonself++;
					}
					if(nonself > 0){
						for(int k = 0; k < desc[sib].member.size(); k++){
							group2 = sample.group_of_haplotype[desc[sib].member[k]];
							if( group_used[group2] != 0 ){
								nn_res(group, group2,0) += num_bases_tree_persists/nonself;
							}
						}
						break;
					}
					node = parent;
				}

				//naive implementation of nearest neighbour
				while(mtr.tree.nodes[node].parent != NULL){
					parent = (*mtr.tree.nodes[node].parent).label;
					sib    = (*mtr.tree.nodes[parent].child_left).label;
					if(node == sib) sib = (*mtr.tree.nodes[parent].child_right).label;
					nonself = 0;
					for(int k = 0; k < desc[sib].member.size(); k++){
						group2 = sample.group_of_haplotype[desc[sib].member[k]];
						if(group != group2 && group_used[ group2 ] != 0) nonself++;
					}
					if(nonself > 0){
						for(int k = 0; k < desc[sib].member.size(); k++){
							group2 = sample.group_of_haplotype[desc[sib].member[k]];
							if( group != group2 && group_used[ group2 ] != 0 ){
								nn_res(group, group2,1) += num_bases_tree_persists/nonself;
							}
						}
						break;
					}
					node = parent;
				}

			}

			num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
		}

	}

	//divide all entries by group size
	for(int i = 0; i < sample.groups.size(); i++){
		for(int j = 0; j < pops_.size(); j++){
			nn_res(i,j,0) /= (1e6 * sample.group_sizes[i]);
			nn_res(i,j,1) /= (1e6 * sample.group_sizes[i]);
		}
	}

	NumericVector nn_res_copy(Dimension(sample.groups.size(), sample.groups.size(),2));
	std::copy(nn_res.begin(), nn_res.end(), nn_res_copy.begin());
	std::vector<std::string> names = {"incl_self", "excl_self"};

	CharacterVector dim1 = wrap(sample.groups);
	CharacterVector dim2 = wrap(names);
	nn_res_copy.attr("dimnames") = List::create(dim1, dim1, dim2);


	return nn_res_copy;

}


