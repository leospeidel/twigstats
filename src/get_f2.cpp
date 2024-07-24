#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "gzstream.hpp"
#include "sample.hpp"
#include "eigenstrat.hpp"
#include "mutations.hpp"

using namespace Rcpp;


int findInterval(const std::vector<float>& boundaries, float value) {
    auto it = std::lower_bound(boundaries.begin(), boundaries.end(), value);
    
    if (it == boundaries.begin()) {
        // Value is less than the first boundary
			  std::cerr << "Epoch out of bounds on the lower side" << std::endl;
        return -1;  // Indicates out of bounds on the lower side
    } else if (it == boundaries.end()) {
        // Value is greater than or equal to the last boundary
        return boundaries.size()-1;  // Indicates out of bounds on the upper side
    } else {
        // Value is between boundaries
        return std::distance(boundaries.begin(), it) - 1;
    }
}

void
traversef2(const Node& n, const Sample& sample, std::vector<float>& freqs, arma::mat& f2, float factor, std::vector<float>& coords, bool correct, double time_cutoff = std::numeric_limits<double>::infinity(), double min_cutoff = 0, double minfreq = 1){
	double overlap = 1.0;
	if(coords[n.label] > time_cutoff){
		overlap = 0.0;
	}else if(n.parent != NULL){
		if(coords[(*n.parent).label] < min_cutoff){
			overlap = 0.0; 
		}else if(coords[(*n.parent).label] > time_cutoff && coords[n.label] < min_cutoff){
			overlap = (time_cutoff - min_cutoff)/(coords[(*n.parent).label] - coords[n.label]);
		}else if(coords[(*n.parent).label] > time_cutoff && coords[(*n.parent).label] > coords[n.label]){
			overlap = (time_cutoff - coords[n.label])/(coords[(*n.parent).label] - coords[n.label]);
		}else if(coords[n.label] < min_cutoff && coords[(*n.parent).label] > coords[n.label]){
			overlap = (coords[(*n.parent).label] - min_cutoff)/(coords[(*n.parent).label] - coords[n.label]);
		}
	}

	if(n.child_left == NULL){

		if(freqs.size() < sample.groups.size()+1) freqs.resize(sample.groups.size()+1);
		std::fill(freqs.begin(), freqs.end(), 0.0);

		for(int i = 0; i < sample.groups.size(); i++){
			if(sample.group_of_haplotype[n.label] == i){
				freqs[i] = 1.0/sample.group_sizes[i];

				float val;
				if(1.0 > minfreq){
					if(correct){
						float a   = freqs[i] * (1.0 - freqs[i])/std::max(sample.group_sizes[i]-1,1);
						val = overlap*factor*(freqs[i]*freqs[i] - a)*n.branch_length;
					}else{
						val = overlap*factor*freqs[i]*freqs[i]*n.branch_length;
					}
					for(int j = 0; j < sample.groups.size(); j++){
						if(i != j){
							f2(i,j) += val;
							f2(j,i) += val;
						}
					}
					f2(i,sample.groups.size()) += val;
					f2(sample.groups.size(),i) += val;
				}

				break;
			}
		}

	}else{

		std::vector<float> freqs_tmp(sample.groups.size(),0.0);
		std::vector<int> zeros(sample.groups.size()), nonzeros(sample.groups.size());
		int nz = 0, nnz = 0;
		float b;

		traversef2((*n.child_left), sample, freqs_tmp, f2, factor, coords, correct, time_cutoff, min_cutoff, minfreq);
		freqs = freqs_tmp;
		std::fill(freqs_tmp.begin(), freqs_tmp.end(), 0.0);
		traversef2((*n.child_right), sample, freqs_tmp, f2, factor, coords, correct, time_cutoff, min_cutoff, minfreq);
		double total = 0;
		for(int i = 0; i < freqs.size(); i++){
			freqs[i] += freqs_tmp[i];
			if(i < sample.group_sizes.size()) total += freqs[i]*sample.group_sizes[i];
			if(freqs[i] > 0){
				nonzeros[nnz] = i;
				nnz++;
			}else{
				zeros[nz] = i;
				nz++;
			}
		}

		if(correct){
			if(total > minfreq){
				float k = overlap*factor*n.branch_length;
				if(k > 0){
					for(std::vector<int>::iterator itnz = nonzeros.begin(); itnz != std::next(nonzeros.begin(),nnz); itnz++){
						assert(*itnz < sample.group_sizes.size());
						float f = freqs[*itnz];
						float a = f * (1.0 - f)/std::max(sample.group_sizes[*itnz]-1,1);
						f = (f*f - a)*k;
						for(std::vector<int>::iterator itz = zeros.begin(); itz != std::next(zeros.begin(),nz); itz++){
							f2(*itnz,*itz) += f;
							f2(*itz,*itnz) += f;
						}
						//f2(*itnz,sample.groups.size()) += f;
						//f2(sample.groups.size(),*itnz) += f;

						f = freqs[*itnz];
						for(std::vector<int>::iterator itnz2 = std::next(itnz,1); itnz2 != std::next(nonzeros.begin(),nnz); itnz2++){
							assert(*itnz2 < sample.group_sizes.size());
							b = freqs[*itnz2] * (1.0 - freqs[*itnz2])/std::max(sample.group_sizes[*itnz2]-1,1);
							f2(*itnz,*itnz2) += ((f-freqs[*itnz2])*(f-freqs[*itnz2])-a-b)*k;
							f2(*itnz2,*itnz) += ((f-freqs[*itnz2])*(f-freqs[*itnz2])-a-b)*k;
						}


					}

				}
			}
		}else{
			if(total > minfreq){
				float k = overlap*factor*n.branch_length;
				if(k > 0){
					for(std::vector<int>::iterator itnz = nonzeros.begin(); itnz != std::next(nonzeros.begin(),nnz); itnz++){
						float f = freqs[*itnz];
						f = f*f*k;
						for(std::vector<int>::iterator itz = zeros.begin(); itz != std::next(zeros.begin(),nz); itz++){
							f2(*itnz,*itz) += f;
							f2(*itz,*itnz) += f;
						}

						f = freqs[*itnz];
						for(std::vector<int>::iterator itnz2 = std::next(itnz,1); itnz2 != std::next(nonzeros.begin(),nnz); itnz2++){
							f2(*itnz,*itnz2) += (f-freqs[*itnz2])*(f-freqs[*itnz2])*k;
							f2(*itnz2,*itnz) += (f-freqs[*itnz2])*(f-freqs[*itnz2])*k;
						}
					}
				}
			}
		}

		if(0){
			if(overlap > 0){

				for(int i = 0; i < sample.groups.size(); i++){
					for(int j = i+1; j < sample.groups.size(); j++){
						f2(i,j) += overlap * factor * (freqs[i] - freqs[j])*(freqs[i] - freqs[j])*n.branch_length;
						f2(j,i) += overlap * factor * (freqs[i] - freqs[j])*(freqs[i] - freqs[j])*n.branch_length;
					}
					f2(i,sample.groups.size()) += overlap * factor * freqs[i] * freqs[i] * n.branch_length;
					f2(sample.groups.size(),i) += overlap * factor * freqs[i] * freqs[i] * n.branch_length;
				}
			}
		}

	}

}

void
traverseTMRCA(const MarginalTree& mtr, std::vector<Leaves>& desc, arma::mat& f2, float factor, std::vector<float>& coords, float t){
	
	float coords_value;
	for(std::vector<Node>::const_iterator it_node = std::next(mtr.tree.nodes.begin(), (mtr.tree.nodes.size() + 1.0)/2.0); it_node != mtr.tree.nodes.end(); it_node++){
	
		coords_value = coords[(*it_node).label];
		if(coords_value > t){
      coords_value = t;
		}
		for(int i = 0; i < desc[(*(*it_node).child_left).label].member.size(); i++){
			for(int j = 0; j < desc[(*(*it_node).child_right).label].member.size(); j++){
				f2(desc[(*(*it_node).child_left).label].member[i],desc[(*(*it_node).child_right).label].member[j]) += factor * coords_value;
				f2(desc[(*(*it_node).child_right).label].member[j],desc[(*(*it_node).child_left).label].member[i]) += factor * coords_value;
			}
		}
	}

}

void
traverseTMRCAdist(const MarginalTree& mtr, std::vector<Leaves>& desc, std::vector<std::vector<std::vector<float>>>& tmrcas, std::vector<float>& epochs, float factor, std::vector<float>& coords){

	int ep = 0;
	float coords_value;
	for(std::vector<Node>::const_iterator it_node = std::next(mtr.tree.nodes.begin(), (mtr.tree.nodes.size() + 1.0)/2.0); it_node != mtr.tree.nodes.end(); it_node++){

		coords_value = coords[(*it_node).label];
		ep = findInterval(epochs, coords_value);

		for(int i = 0; i < desc[(*(*it_node).child_left).label].member.size(); i++){
			for(int j = 0; j < desc[(*(*it_node).child_right).label].member.size(); j++){
				tmrcas[desc[(*(*it_node).child_left).label].member[i]][desc[(*(*it_node).child_right).label].member[j]][ep] += factor;
				tmrcas[desc[(*(*it_node).child_right).label].member[j]][desc[(*(*it_node).child_left).label].member[i]][ep] += factor;
			}
		}
	}

}


//' Function to calculate f2 statistics from Relate trees for pairs of populations specified in poplabels.
//'
//' This function will calculate f2 statistics in blocks of prespecified size for all pairs of populations specified in the poplabels file.
//' Please refer to the Relate documentation for input file formats (https://myersgroup.github.io/relate/).
//' The output is in a format that is directly accepted by the admixtools R package to calculate 
//' f3, f4, f4ratio, D statistics and more (https://uqrmaie1.github.io/admixtools/).
//'
//' @param file_anc Filename of anc file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).
//' @param file_mut Filename of mut file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).
//' @param poplabels Filename of poplabels file
//' @param chrs (Optional) Vector of chromosome IDs
//' @param blgsize (Optional) SNP block size in Morgan. Default is 0.05 (50 cM). If blgsize is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
//' @param file_map (Optional) File prefix of recombination map. Not needed if blgsize is given in base-pairs, i.e. blgsize > 100
//' @param mu (Optional) Per base per generation mutation rate to scale f2 values. Default: 1.25e-8
//' @param t (Optional) Time cutoff in generations. Default: Inf
//' @param tmin (Optional) Minimum time cutof in generations. Any lineages younger than tmin will be excluded from the analysis. Default: t = 0.
//' @param minMAF (Optional) Minimum frequency cutoff. Default: 1 (i.e. excl singletons)
//' @param use_muts (Optional) Calculate traditional f2 statistics by only using mutations mapped to Relate trees. Default: false.
//' @param transitions (Optional) Set this to FALSE to exclude transition SNPs. Only meaningful with use_muts
//' @param apply_corr (Optional) Use small sample size correction. Default: true.
//' @param dump_blockpos (Optional) Filename of blockpos file.
//' @return 3d array of dimension #groups x #groups x #blocks. Analogous to output of f2_from_geno in admixtools.
//' @examples
//' file_anc  <- system.file("sim/msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats")
//' file_mut  <- system.file("sim/msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats")
//' poplabels <- system.file("sim/msprime_ad0.8_split250_1.poplabels", package = "twigstats")
//' file_map  <- system.file("sim/genetic_map_combined_b37_chr1.txt.gz", package = "twigstats")
//'
//' #Calculate f2s between all pairs of populations
//' f2_blocks <- f2_blocks_from_Relate(file_anc, file_mut, poplabels, file_map)
//' f4_ratio(f2_blocks, popX="PX", popI="P1", pop1="P2", pop2="P3", popO="P4")
//'
//' #Use a cutoff of 500 generations
//' f2_blocks <- f2_blocks_from_Relate(file_anc, file_mut, poplabels, file_map, t = 500)
//' f4_ratio(f2_blocks, popX="PX", popI="P1", pop1="P2", pop2="P3", popO="P4")
//' @export
// [[Rcpp::export]]
NumericVector f2_blocks_from_Relate( SEXP file_anc, SEXP file_mut, SEXP poplabels, SEXP file_map = R_NilValue, Nullable<CharacterVector> chrs = R_NilValue, Nullable<double> blgsize = R_NilValue, Nullable<double> mu = R_NilValue, Nullable<double> tmin = R_NilValue, Nullable<double> t = R_NilValue, Nullable<int> transitions = R_NilValue, Nullable<int> use_muts = R_NilValue, Nullable<int> minMAF = R_NilValue, SEXP dump_blockpos = R_NilValue, Nullable<int> apply_corr = R_NilValue) {

	std::string filename_poplabels = as<std::string>(poplabels);
	std::vector<std::string> filename_anc, filename_mut, filename_rec;
	std::string map_ = "";
	if(file_map != R_NilValue){
		map_ = as<std::string>(file_map);
	}

	if(chrs.isNotNull()){
		CharacterVector chrs_(chrs);
		for(int i = 0; i < chrs_.size(); i++){
			filename_anc.push_back( as<std::string>(file_anc) + "_chr" + as<std::string>(chrs_[i]) + ".anc" );
			filename_mut.push_back( as<std::string>(file_mut) + "_chr" + as<std::string>(chrs_[i]) + ".mut" );
			filename_rec.push_back( map_ + "_chr" + as<std::string>(chrs_[i]) + ".txt" );
		}
	}else{
		filename_anc.push_back( as<std::string>(file_anc) );
		filename_mut.push_back( as<std::string>(file_mut) );
		filename_rec.push_back( map_ );
	}

	double mu_ = 1.25e-8;
	if(mu.isNotNull()){
		mu_ = as<double>(mu);
	}
	double binsize_ = 0.05;
	if(blgsize.isNotNull()){
		binsize_ = as<double>(blgsize);
	}
	bool use_muts_ = false;
	if(use_muts.isNotNull()){
		use_muts_ = as<bool>(use_muts);
	}
	double t_ = std::numeric_limits<double>::infinity();
	if(t.isNotNull()){
		t_ = as<double>(t);
	}
	double t_min_ = 0;
	if(tmin.isNotNull()){
		t_min_ = as<double>(tmin);
	}
	if(t_min_ >= t_){
		Rcpp::stop("t needs to be greater than tmin");
	}
	int transitions_ = 1;
	if(transitions.isNotNull()){
		transitions_ = as<int>(transitions);
	}
	int minmafcut = 1;
	if(minMAF.isNotNull()){
		minmafcut = as<int>(minMAF);
	}

	bool correct = true;
	if(apply_corr.isNotNull()){
		correct = as<int>(apply_corr);
	}

	bool isM = true;
	if(binsize_ > 100) isM = false;
	if(isM){
		if(map_ == ""){
			Rcpp::stop("Need to specify recombination map when blgsize < 100");
		}
	}

	Sample sample;
	sample.Read(filename_poplabels);

	std::cerr << "Computing pairwise f2 for " << sample.groups.size() + 1 << " groups." << std::endl;
	std::cerr << "blgsize: " << binsize_;
	if(isM) std::cerr << "M" << std::endl;
	if(!isM) std::cerr << "bp" << std::endl;
	std::cerr << std::endl;

	MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
	Muts::iterator it_mut; //iterator for mut file
	float num_bases_tree_persists = 0.0;

	////////// 1. Read one tree at a time /////////

	//NumericVector f2_block(Dimension(sample.groups.size(), sample.groups.size(), 10000));
	arma::cube f2_block = arma::zeros<arma::cube>(sample.groups.size()+1, sample.groups.size()+1, 1000);
	std::vector<std::string> splicenames(1000);
	int blockID = 0;
	float num_infSNPs = 0.0;

	bool blockBP = false;
	std::ofstream os;
	if(dump_blockpos != R_NilValue){
		blockBP = true;
		os.open(as<std::string>(dump_blockpos));
	}

	for(int chr = 0; chr < filename_anc.size(); chr++){

		//We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
		//The mut file is read once, file is closed after constructor is called.
		AncMutIterators ancmut(filename_anc[chr], filename_mut[chr]);
		Data data(ancmut.NumTips(), ancmut.NumSnps());

		map recmap(filename_rec[chr].c_str());
		int irec = 0;

		if( blockID+100 >= size(f2_block)[2] ){
			f2_block.resize(sample.groups.size()+1, sample.groups.size()+1, size(f2_block)[3] + 1000);
			splicenames.resize(splicenames.size() + 1000);
		}

		//mtr stores the marginal tree, it_mut points to the first SNP at which the marginal tree starts.
		//If marginal tree has no SNPs mapping to it, it_mut points to the first SNP of the following tree that has a SNP mapping to it
		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

		std::vector<float> freqs(2, 0.0), coords(2*data.N-1, 0.0);
		std::vector<Leaves> desc;
		int current_pos = 0;
		double current_rec = 0;
		double rec;
		double factor = 1.0;
		int percentage = 0;
		//iterate through whole file
		std::cerr << "Block: " << blockID << "\r";
		while(num_bases_tree_persists >= 0.0){

			mtr.tree.GetCoordinates(coords);
			mtr.tree.FindAllLeaves(desc);

			if(((*it_mut).tree % (int)(ancmut.NumTrees()/100)) == 0){
				std::cerr << "Block: " << blockID << ", " << "[" << percentage << "%]\r";
				Rcpp::checkUserInterrupt();
				percentage++;
				std::cerr.flush();
			}

			if(isM){
				if(recmap.bp[irec] >= (*it_mut).pos){
					if(irec > 0){
						rec = ((double)(*it_mut).pos - recmap.bp[irec-1])/(recmap.bp[irec] - recmap.bp[irec-1]) * (recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) + recmap.gen_pos[irec-1];
					}else{
						rec = ((double)(*it_mut).pos)/(recmap.bp[irec]) * (recmap.gen_pos[irec]);         
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
						rec = recmap.gen_pos[irec]/recmap.bp[irec] * ((*it_mut).pos - recmap.bp[irec]) + recmap.gen_pos[irec];
					}else{
						assert(recmap.bp[irec] >= (*it_mut).pos);
						assert(irec > 0);
						rec = ((double)(*it_mut).pos - recmap.bp[irec-1])/(recmap.bp[irec] - recmap.bp[irec-1]) * (recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) + recmap.gen_pos[irec-1];
					}
				}
				assert(rec >= 0.0);
			}

			//for each tree, calculate f2 stats given sample assignment to pops
			if( (!isM && (*it_mut).pos - current_pos > binsize_) || (isM && rec/100.0 - current_rec > binsize_) ){

				Rcpp::checkUserInterrupt();
				if((int)std::round(num_infSNPs) > 0){
					for(arma::cube::slice_iterator it_c = f2_block.begin_slice(blockID); it_c != f2_block.end_slice(blockID); it_c++){
						*it_c /= num_infSNPs;
					}
					splicenames[blockID] = "l" + std::to_string((int)std::round(num_infSNPs));
					num_infSNPs = 0;
					if(blockBP) os << blockID << " " << current_pos << " " << (*it_mut).pos << "\n";
					blockID++;
				}
				if( blockID+100 >= size(f2_block)[2] ){
					f2_block.resize(sample.groups.size()+1, sample.groups.size()+1, size(f2_block)[2] + 1000);
					splicenames.resize(splicenames.size() + 1000);
				}
				if(isM){
					while(rec/100.0 - current_rec > binsize_){
						current_rec += binsize_; 
					}
				}else{
					while((*it_mut).pos - current_pos > binsize_){
						current_pos += binsize_; 
					}
				}

			}

			if(!use_muts_){
				factor = mu_ * num_bases_tree_persists;

				double overlap = 1.0;	
				for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
					overlap = 1.0;
					if(coords[(*it_node).label] > t_ || desc[(*it_node).label].num_leaves <= minmafcut){
						overlap = 0.0;
					}else if((*it_node).parent != NULL){
						if(coords[(*(*it_node).parent).label] < t_min_){
							overlap = 0.0; 
						}else if(coords[(*(*it_node).parent).label] > t_ && coords[(*it_node).label] < t_min_){
							overlap = (t_ - t_min_)/(coords[(*(*it_node).parent).label] - coords[(*it_node).label]);
						}else if(coords[(*(*it_node).parent).label] > t_ && coords[(*(*it_node).parent).label] > coords[(*it_node).label]){
							overlap = (t_ - coords[(*it_node).label])/(coords[(*(*it_node).parent).label] - coords[(*it_node).label]);
						}else if(coords[(*it_node).label] < t_min_ && coords[(*(*it_node).parent).label] > coords[(*it_node).label]){
							overlap = (coords[(*(*it_node).parent).label] - t_min_)/(coords[(*(*it_node).parent).label] - coords[(*it_node).label]);
						}
					}else{
						overlap = 0.0;
					}
					num_infSNPs += overlap*(*it_node).branch_length * factor;
				}

				traversef2(*std::prev(mtr.tree.nodes.end(),1), sample, freqs, f2_block.slice(blockID), factor, coords, correct, t_, t_min_, minmafcut);
			}else{

				for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
					(*it_node).branch_length = 0.0;
				}

				int treeID = (*it_mut).tree;
				while((*it_mut).tree == treeID){

					bool use = ((*it_mut).branch.size() == 1);
					if(transitions_ == 0 && use){
						if((*it_mut).mutation_type != "C/A" && (*it_mut).mutation_type != "A/C" &&
								(*it_mut).mutation_type != "G/C" && (*it_mut).mutation_type != "C/G" &&
								(*it_mut).mutation_type != "A/T" && (*it_mut).mutation_type != "T/A" &&
								(*it_mut).mutation_type != "T/G" && (*it_mut).mutation_type != "G/T"
							){
							use = false;
						}
					}

					if(use){
						mtr.tree.nodes[(*it_mut).branch[0]].branch_length += 1.0;
					}
					it_mut++;
					if(it_mut == ancmut.mut_end()) break;
				}

				double overlap = 1.0;
				for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
					overlap = 1.0;
					if(coords[(*it_node).label] > t_ || desc[(*it_node).label].num_leaves <= minmafcut){
						overlap = 0.0;
					}else if((*it_node).parent != NULL){
						if(coords[(*(*it_node).parent).label] < t_min_){
							overlap = 0.0; 
						}else if(coords[(*(*it_node).parent).label] > t_ && coords[(*it_node).label] < t_min_){
							overlap = (t_ - t_min_)/(coords[(*(*it_node).parent).label] - coords[(*it_node).label]);
						}else if(coords[(*(*it_node).parent).label] > t_ && coords[(*(*it_node).parent).label] > coords[(*it_node).label]){
							overlap = (t_ - coords[(*it_node).label])/(coords[(*(*it_node).parent).label] - coords[(*it_node).label]);
						}else if(coords[(*it_node).label] < t_min_ && coords[(*(*it_node).parent).label] > coords[(*it_node).label]){
							overlap = (coords[(*(*it_node).parent).label] - t_min_)/(coords[(*(*it_node).parent).label] - coords[(*it_node).label]);
						}
					}else{
						overlap = 0.0;
					}

					assert(overlap >= 0);
					assert(overlap <= 1.0);
					num_infSNPs += overlap * (*it_node).branch_length;
				}

				traversef2(*std::prev(mtr.tree.nodes.end(),1), sample, freqs, f2_block.slice(blockID), factor, coords, correct, t_, t_min_, minmafcut);

			}

			num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
		}

		if((int)std::round(num_infSNPs) > 0){
			for(arma::cube::slice_iterator it_c = f2_block.begin_slice(blockID); it_c != f2_block.end_slice(blockID); it_c++){
				*it_c /= num_infSNPs;
			}
			splicenames[blockID] = "l" + std::to_string((int)std::round(num_infSNPs));
			num_infSNPs = 0;
			if(blockBP) os << blockID << " " << current_pos << " " << (*it_mut).pos << "\n";
			blockID++;
		}

	}

	sample.groups.push_back("Root");

	f2_block.reshape(sample.groups.size(), sample.groups.size(), blockID);
	splicenames.resize(blockID);

	NumericVector f2_block_copy(Dimension(sample.groups.size(), sample.groups.size(), blockID));
	std::copy(f2_block.begin(), f2_block.end(), f2_block_copy.begin());
	f2_block_copy.attr("dimnames") = List::create(wrap(sample.groups), wrap(sample.groups), wrap(splicenames));

	std::cerr << "Block: " << blockID << ", " << "[100%]\r";
	std::cerr << std::endl;

	return f2_block_copy;

}

//' Function to calculate f2 statistics from plink files, ascertained using mutation ages in Relate trees.
//'
//' This function will calculate f2 statistics in blocks of prespecified size for all pairs of populations specified in the input files.
//' Input is assumed to be in PLINK binary format https://www.cog-genomics.org/plink/1.9/formats#bed and mutation ages are in Relate mut format https://myersgroup.github.io/relate/.
//' The output is in a format that is directly accepted by the admixtools R package to calculate 
//' f2, f3, f4, f4ratio, D statistics and more (https://uqrmaie1.github.io/admixtools/).
//'
//' @param pref Prefix of PLINK binary files, assuming filenames of form $\{pref\}.bed, $\{pref\}.bim, $\{pref\}.fam.
//' @param file_mut (Optional) Prefix of filenames of mut files, assuming filenames of form $\{file_mut\}_chr1.mut(.gz). Chromosome names have to be consistent to those in the PLINK files. If no file is specified, all mutations in the PLINK file are used.
//' @param blgsize (Optional) SNP block size in Morgan. Default is 0.05 (50 cM). If blgsize is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
//' @param t (Optional) Time cutoff in generations. Any mutations older that t will be excluded from the analysis. Default: t = Inf.
//' @param tmin (Optional) Minimum time cutof in generations. Any mutations younger than tmin will be excluded from the analysis. Default: t = 0.
//' @param transitions (Optional) Set this to FALSE to exclude transition SNPs
//' @param maxmiss (Optional) Discard SNPs which are missing in a fraction of populations higher than maxmiss
//' @param pops (Optional) Populations for which data should be extracted. Names need to match the first column in the fam file (or fam option below)
//' @param fam (Optional) 1d-array assigning individuals to populations. Corresponds to the first column in the fam file and is useful if you want to change population assignments.
//' @param chrs (Optional) List chromosome names to use.
//' @param include_undated (Optional) Include mutations that are not dated. Default: false.
//' @param minMAF (Optional) minimum minor allele count. Default: 1.
//' @param apply_corr (Optional) Use small sample size correction. Default: true.
//' @param debug_mode (Optional) Prints progress used for debugging.
//' @return 3d array of dimension #groups x #groups x #blocks containing f2 statistics. Analogous to output of f2_from_geno in admixtools.
//' @examples
//' path <- paste0(system.file("sim/", package = "twigstats"),"/")
//' pref <- "msprime_ad0.8_split250_1"
//' file_plink <- paste0(path,pref) #only need prefix
//' file_mut  <- paste0(path,pref) #only need prefix (here same name as plink file but can be different)
//'
//' system(paste0("gunzip ", file_plink, ".bim.gz"))
//'
//' #Compute f2 statistics between all pairs of populations. You can use pops to only calculate f2s between specified populations.
//' f2_blocks <- f2_blocks_from_RelateAges(pref = file_plink, file_mut)
//' f4_ratio(f2_blocks, popX="PX", popI="P1", pop1="P2", pop2="P3", popO="P4")
//'
//' #Use a cutoff of 500 generations
//' f2_blocks <- f2_blocks_from_RelateAges(pref = file_plink, file_mut, t = 500)
//' f4_ratio(f2_blocks, popX="PX", popI="P1", pop1="P2", pop2="P3", popO="P4")
//' @export
// [[Rcpp::export]]
NumericVector f2_blocks_from_RelateAges( SEXP pref, SEXP file_mut, Nullable<double> blgsize = R_NilValue, Nullable<int> transitions = R_NilValue, Nullable<double> maxmiss = R_NilValue, Nullable<StringVector> fam = R_NilValue, Nullable<StringVector> pops = R_NilValue, Nullable<StringVector> chrs = R_NilValue, Nullable<double> tmin = R_NilValue, Nullable<double> t = R_NilValue, Nullable<int> include_undated = R_NilValue, Nullable<int> minMAF = R_NilValue, Nullable<int> apply_corr = R_NilValue, Nullable<int> debug_mode = 0){

	double binsize_ = 0.05;
	if(blgsize.isNotNull()){
		binsize_ = as<double>(blgsize);
	}
	double t_ = std::numeric_limits<double>::infinity();
	if(t.isNotNull()){
		t_ = as<double>(t);
	}
	double t_min_ = 0;
	if(tmin.isNotNull()){
		t_min_ = as<double>(tmin);
	}
	if(t_min_ >= t_){
		Rcpp::stop("t needs to be greater than tmin");
	}
	bool include = false;
	if(include_undated.isNotNull()){
		include = as<int>(include_undated);
	}
	bool debug = false;
	if(debug_mode.isNotNull()){
		debug = as<int>(debug_mode);
	}
	int minmafcut = 1;
	if(minMAF.isNotNull()){
		minmafcut = as<int>(minMAF);
	}

	bool correct = true;
	if(apply_corr.isNotNull()){
		correct = as<int>(apply_corr);
	}

	std::vector<std::string> chrs_;
	if(chrs.isNotNull()){
		StringVector chrs_tmp_ = as<StringVector>(chrs);
		for(int i = 0; i < chrs_tmp_.size(); i++){
			chrs_.push_back(as<std::string>(chrs_tmp_[i]));
		}
	}

	bool isM = true;
	if(binsize_ > 100) isM = false;

	std::string filename_mut = as<std::string>(file_mut);

	int ploidy_ = 2;
	//if(ploidy.isNotNull()){
	//  ploidy_ = as<int>(ploidy);
	//}
	int transitions_ = 1;
	if(transitions.isNotNull()){
		transitions_ = as<int>(transitions);
	}
	double maxmiss_ = 0.0;
	if(maxmiss.isNotNull()){
		maxmiss_ = as<double>(maxmiss);
	}

	//eigenstrat es(as<std::string>(pref));
	plink es(as<std::string>(pref));
	if(fam.isNotNull()){
		StringVector fam_;
		std::vector<std::string> vfam_;
		fam_ = as<StringVector>(fam);
		for(int i = 0; i < fam_.size(); i++){
			vfam_.push_back(as<std::string>(fam_[i]));
		}
		es.ChangeFAM(vfam_);
	}

	Data data(es.GetN(), es.GetL());

	StringVector pops_(es.unique_groups.size());
	pops_ = es.unique_groups;
	if(pops.isNotNull()){
		pops_ = as<StringVector>(pops);
		//TODO: allow to combine populations
	}
	std::cerr << pops_.size() << " " << es.unique_groups.size() << std::endl;

	//need a mapping of individual to group in pops_
	std::vector<int> groups_mapping(es.unique_groups.size(), -1);
	for(int g = 0; g < pops_.size(); g++){
		bool exists = false;
		for(int l = 0; l < es.unique_groups.size(); l++){
			if( pops_[g] == es.unique_groups[l] ){
				groups_mapping[l] = g;
				exists = true;
				break;
			}
		}
		if(!exists){
			Rcpp::stop("Group " + pops_[g] + " does not exist.");
		}
	}

	int max_ind = 0;
	std::vector<int> use_ind(data.N, -1);
	for(int s = 0; s < data.N; s++){
		//std::cerr << s << " " << es.membership[s] << " " << es.unique_groups[es.membership[s]] << " " << groups_mapping[es.membership[s]] << std::endl;
		if(groups_mapping[es.membership[s]] != -1){
			max_ind = s;
			use_ind[s] = groups_mapping[es.membership[s]];
		}
	}
	max_ind++;
	es.SetMaxInd(max_ind);

	std::cerr << "Computing pairwise f2 for " << pops_.size() << " groups." << std::endl;
	std::cerr << "blgsize: " << binsize_;
	if(isM) std::cerr << "M" << std::endl;
	if(!isM) std::cerr << "bp" << std::endl;
	std::cerr << std::endl;

	Mutations mut;
	Muts::iterator it_mut; //iterator for mut file

	arma::cube f2_block = arma::zeros<arma::cube>(pops_.size(), pops_.size(), 1000);
	arma::mat num_infSNPs_ps = arma::zeros<arma::mat>(pops_.size(), pops_.size());
	std::vector<std::string> splicenames(1000);
	int blockID = 0, current_pos = 0;
	float num_infSNPs = 0.0, num_snps = 0, num_used_snps = 0;

	std::vector<int> sequence;
	std::vector<double> freqs(pops_.size(), 0), num_inds(pops_.size(), 0);
	std::string REF, ALT, chr, current_chr;
	int bp, percentage = 0, snp = 0;
	double rec, current_rec = 0;
	bool mut_exists = false, never_exists = false;
	if(filename_mut == "") never_exists = true;

	std::vector<int> pass(pops_.size(), 0);

	//std::ofstream os("debug.txt");
	//for(int g = 0; g < es.unique_groups.size(); g++){
	//  os << es.unique_groups[g] << " ";
	//}
	//os << "\n";

	bool use_chr = true;
	int num_used_chr = 0;
	//Read EIGENSTRAT file and iterate through.
	//Look at chrID, try open the corresponding mut file and WARN if it doesn't exist. 
	while(es.ReadSNP(sequence, chr, bp, rec, REF, ALT)){

		//std::cerr << bp << std::endl;

		if((snp % (int)(data.L/100)) == 0){
			std::cerr << "Block: " << blockID << ", " << "[" << percentage << "%]\r";
			Rcpp::checkUserInterrupt();
			percentage++;
			std::cerr.flush();
		}

		if(current_chr != chr || (!isM && bp - current_pos > binsize_) || (isM && rec - current_rec > binsize_) ){

			if(current_chr != chr && chrs_.size() > 0){
				use_chr = false;
				for(int i = 0; i < chrs_.size(); i++){
					if(chr == chrs_[i]){
						use_chr = true;
						break;
					}
				}
				if(use_chr == true) num_used_chr++;
				if(use_chr == false && num_used_chr >= chrs_.size()) break;
			}

			Rcpp::checkUserInterrupt();

			//std::cerr << blockID << " " << bp << " " << rec << " " << current_rec << std::endl;

			bool all_data = true;
			for(int i = 0; i < pops_.size(); i++){
				for(int j = 0; j < pops_.size(); j++){

					//if(i != j){
					//  if(num_infSNPs_ps(i,j) != num_infSNPs){
					//    std::cerr << std::endl;
					//    std::cerr << es.unique_groups[i] << " " << es.unique_groups[j] << " " << num_infSNPs_ps(i,j) << " " << num_infSNPs << std::endl;
					//  }
					//  assert(num_infSNPs_ps(i,j) == num_infSNPs);
					//}

					if(i != j && num_infSNPs_ps(i,j) == 0){
						all_data = false;
						break;
					}
				}
				if(!all_data) break;
			}
			if(num_infSNPs > 0 && all_data){
				for(int i = 0; i < pops_.size(); i++){
					for(int j = 0; j < pops_.size(); j++){
						if(i != j) f2_block(i,j,blockID) /= num_infSNPs_ps(i,j);
					}
				}
				//std::cerr << chr << " " << blockID << " " << num_infSNPs << " " << num_used_snps << std::endl << std::endl;
				splicenames[blockID] = "l" + std::to_string(1+(int)std::round(num_infSNPs));
				//std::cerr << num_snps << std::endl;
				//splicenames[blockID] = "l100";
				num_infSNPs = 0;
				num_infSNPs_ps = arma::zeros<arma::mat>(pops_.size(), pops_.size());
				num_snps = 0;
				blockID++;
			}
			if( blockID+100 >= size(f2_block)[2] ){
				f2_block.resize(pops_.size(), pops_.size(), size(f2_block)[2] + 1000);
				splicenames.resize(splicenames.size() + 1000);
			}
			if(current_chr != chr){
				current_pos = 0.0;
				current_rec = 0.0;
			}
			if(isM){
				while(rec - current_rec > binsize_){
					current_rec += binsize_; 
				}
			}else{
				while(bp - current_pos > binsize_){
					current_pos += binsize_; 
				}
			}
		}

		//std::cerr << current_chr << " " << chr << std::endl;
		if(current_chr != chr){
			current_chr = chr;
			if(!never_exists){
				//std::cerr << filename_mut + "_chr" + chr + ".mut" << std::endl; 
				if(exists(filename_mut + "_chr" + chr + ".mut") || exists(filename_mut + "_chr" + chr + ".mut.gz")){
					mut.Read(filename_mut + "_chr" + chr + ".mut");
					it_mut = mut.info.begin();
					mut_exists = true;
				}else{
					Rcpp::warning("Warning: Failed to open mut file " + filename_mut + "_chr" + chr + ".mut(.gz)");
					mut_exists = false;
				}
			}
		}

		if(use_chr){
			bool use = true;
			double overlap = 1.0;

			if(REF != "A" && REF != "C" && REF != "G" && REF != "T" && REF != "0" && REF != "1") use = false;
			if(ALT != "A" && ALT != "C" && ALT != "G" && ALT != "T" && ALT != "1" && ALT != "0") use = false;

			if(transitions_ == 0){
				if( (REF == "C" && ALT == "T") || (REF == "T" && ALT == "C") || (REF == "G" && ALT == "A") || (REF == "A" && ALT == "G") ){
					use = false;
				}
			}

			if(debug) std::cerr << chr << " " << snp << " " << data.L << " " << blockID << " " << bp << " " << " " << (*it_mut).pos << " " << rec << " " << current_chr << " " << current_pos << " " << current_rec << " " << binsize_ << " " << REF << " " << ALT << " " << use << " " << mut_exists << std::endl;

			if(!include && mut_exists && use){
				if(it_mut != mut.info.end()){
					while((*it_mut).pos < bp){
						it_mut++;
						if(it_mut == mut.info.end()) break;
					}
					//std::cerr << (*it_mut).pos << " " << (int) bp << std::endl;

					//std::cerr << use << " " << (*it_mut).pos << " " << (int) bp << " " << prop_excl << " " << snp << " " << ((double) prop_excl/snp) << std::endl;
					//exclude or calculate overlap here
					if( (*it_mut).pos == bp ){
						std::string mut_type1 = REF + "/" + ALT;
						std::string mut_type2 = ALT + "/" + REF;	
						if( (*it_mut).branch.size() == 1 && (*it_mut).flipped == 0 && (*it_mut).age_end > (*it_mut).age_begin ){
							if( (*it_mut).mutation_type == mut_type1 || (*it_mut).mutation_type == mut_type2 ){
								double begin = (*it_mut).age_begin;
								double end   = (*it_mut).age_end;

								if(t_ < begin || end < t_min_){
									use = false;
								}
								if(t_ < end){
									end = t_;
								}
								if(begin < t_min_){
									begin = t_min_;
								}

								if(0){
									if(t_ < begin){
										use = false;
									}else if(t_ < end){
										begin = t_;
									}
									if(end < t_min_){
										use = false;
									}else if( begin < t_min_ ){
										begin = t_min_;
									}
								}
								overlap = (end - begin)/((*it_mut).age_end-(*it_mut).age_begin);
							}else{
								if(!include) use = false;
							}
						}else{
							if(!include) use = false;
						}
					}else{
						if(!include) use = false;
					}
				}else{
					use = false;
				}
			}else{
				if(!include) use = false;
			}

			if(1){
				if(use){
					std::fill(pass.begin(), pass.end(), 0);
					int s = 0;
					for(std::vector<int>::iterator it_seq = sequence.begin(); it_seq != sequence.end(); it_seq++){
						if((*it_seq) != 9 && use_ind[s] != -1){
							pass[use_ind[s]]++;
						}
						s++;
					}
					float num_miss = 0.0;
					for(int g = 0; g < pops_.size(); g++){
						if(pass[g] == 0){
							num_miss += 1.0;
							if(num_miss > maxmiss_*pops_.size()){
								use = false;
								break;
							}
						}
					}
				}
			}

			if(use){

				int s = 0;
				std::fill(freqs.begin(), freqs.end(), 0);
				std::fill(num_inds.begin(), num_inds.end(), 0);

				float total = 0.0, num_called = 0.0;
				for(std::vector<int>::iterator it_seq = sequence.begin(); it_seq != sequence.end(); it_seq++){
					if((*it_seq) != 9 && use_ind[s] != -1){
						freqs[use_ind[s]]    += (*it_seq);
						num_inds[use_ind[s]] += ploidy_;
						total      += (*it_seq);
						num_called += ploidy_;
					}
					s++;
				}

				if(total < num_called-minmafcut && total > minmafcut){

					num_used_snps++;
					std::vector<double>::iterator it_num = num_inds.begin();
					for(std::vector<double>::iterator it_freq = freqs.begin(); it_freq != freqs.end(); it_freq++){
						if(*it_num > 0){
							*it_freq /= *it_num;
							assert(*it_freq >= 0.0);
							assert(*it_freq <= 1.0);
						}
						it_num++;
					}

					if(correct){
						float a,b;
						for(int i = 0; i < pops_.size(); i++){
							a = freqs[i]*(1-freqs[i])/std::max(num_inds[i]-1,1.0);
							for(int j = i+1; j < pops_.size(); j++){
								if(num_inds[i] > 0 && num_inds[j] > 0){
									num_infSNPs_ps(i,j)   += overlap;
									num_infSNPs_ps(j,i)   += overlap;
									b = freqs[j]*(1-freqs[j])/std::max(num_inds[j]-1,1.0);

									f2_block(i,j,blockID) += overlap * ((freqs[i] - freqs[j])*(freqs[i] - freqs[j]) - a - b);
									f2_block(j,i,blockID) += overlap * ((freqs[i] - freqs[j])*(freqs[i] - freqs[j]) - a - b);
								}
							}
						}
					}else{
						for(int i = 0; i < pops_.size(); i++){
							for(int j = i+1; j < pops_.size(); j++){
								if(num_inds[i] > 0 && num_inds[j] > 0){
									num_infSNPs_ps(i,j)   += overlap;
									num_infSNPs_ps(j,i)   += overlap;
									f2_block(i,j,blockID) += overlap * (freqs[i] - freqs[j])*(freqs[i] - freqs[j]);
									f2_block(j,i,blockID) += overlap * (freqs[i] - freqs[j])*(freqs[i] - freqs[j]);
								}
							}
						}
					}
					num_infSNPs += overlap;
				}

			}

		}


		num_snps++;
		snp++;

	}

	bool all_data = true;
	for(int i = 0; i < pops_.size(); i++){
		for(int j = 0; j < pops_.size(); j++){
			if(i != j && num_infSNPs_ps(i,j) == 0){
				all_data = false;
				break;
			}
		}
		if(!all_data) break;
	}
	if(num_infSNPs > 0 && all_data){
		//for(arma::cube::slice_iterator it_c = f2_block.begin_slice(blockID); it_c != f2_block.end_slice(blockID); it_c++){
		//  *it_c /= num_infSNPs;
		//}
		for(int i = 0; i < pops_.size(); i++){
			for(int j = 0; j < pops_.size(); j++){
				if(i != j) f2_block(i,j,blockID) /= num_infSNPs_ps(i,j);
			}
		}
		//std::cerr << chr << " " << blockID << " " << num_infSNPs << std::endl << std::endl;
		splicenames[blockID] = "l" + std::to_string(1+(int)std::round(num_infSNPs));
		//splicenames[blockID] = "l100";
		num_infSNPs = 0;
		num_infSNPs_ps = arma::zeros<arma::mat>(pops_.size(), pops_.size());
		blockID++;
	}

	f2_block.reshape(pops_.size(), pops_.size(), blockID);
	splicenames.resize(blockID);

	NumericVector f2_block_copy(Dimension(pops_.size(), pops_.size(), blockID));
	std::copy(f2_block.begin(), f2_block.end(), f2_block_copy.begin());
	f2_block_copy.attr("dimnames") = List::create(wrap(pops_), wrap(pops_), wrap(splicenames));

	std::cerr << "Block: " << blockID << ", " << "[100%]\r";
	std::cerr << std::endl;
	std::cerr << "Number of SNPs used: " << num_used_snps << std::endl;
	//std::cerr << "Excluded " << std::round(10000 * prop_excl/data.L)/100 << "% of SNPs in EIGENSTRAT file." << std::endl;

	return f2_block_copy;

}


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
//' @param blgsize (Optional) SNP block size in Morgan. Default is 0.05 (50 cM). If blgsize is 1 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
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
	nn_res_copy.attr("dimnames") = List::create(wrap(sample.groups), wrap(sample.groups), wrap(names));

	return nn_res_copy;

}


//' Function to calculates mean TMRCAs from Relate trees for pairs of populations specified in poplabels.
//'
//' This function will calculate TMRCAs in blocks of prespecified size for all pairs of populations specified in the poplabels file.
//' Please refer to the Relate documentation for input file formats (https://myersgroup.github.io/relate/).
//'
//' @param file_anc Filename of anc file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).
//' @param file_mut Filename of mut file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).
//' @param file_out Filename of output file.
//' @param poplabels Filename of poplabels file
//' @param t (Optional) Time cutoff in generations. Any coalescences older that t will be set to t in the analysis. Default: t = Inf.
//' @param chrs (Optional) Vector of chromosome IDs
//' @param blgsize (Optional) SNP block size in Morgan. Default is 0.05 (50 cM). If blgsize is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
//' @param file_map (Optional) File prefix of recombination map. Not needed if blgsize is given in base-pairs, i.e. blgsize > 100
//' @return 3d array of dimension #groups x #groups x #blocks. Analogous to output of f2_from_geno in admixtools.
//' @keywords internal
//' @examples
//' file_anc  <- system.file("sim/msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats")
//' file_mut  <- system.file("sim/msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats")
//' poplabels <- system.file("sim/msprime_ad0.8_split250_1.poplabels", package = "twigstats")
//' file_map  <- system.file("sim/genetic_map_combined_b37_chr1.txt.gz", package = "twigstats")
//'
//' #Calculate f2s between all pairs of populations
//' TMRCA_from_Relate(file_anc, file_mut, poplabels, file_out = "test", file_map)
//' @export
// [[Rcpp::export]]
void TMRCA_from_Relate( SEXP file_anc, SEXP file_mut, SEXP poplabels, SEXP file_out, SEXP file_map = R_NilValue, Nullable<CharacterVector> chrs = R_NilValue, Nullable<double> t = R_NilValue, Nullable<double> blgsize = R_NilValue) {

	std::string filename_poplabels = as<std::string>(poplabels);
	std::vector<std::string> filename_anc, filename_mut, filename_rec;
	std::string map_ = "";
	if(file_map != R_NilValue){
		map_ = as<std::string>(file_map);
	}

	double t_ = std::numeric_limits<double>::infinity();
	if(t.isNotNull()){
		t_ = as<double>(t);
	}

	CharacterVector chrs_;
	if(chrs.isNotNull()){
		chrs_ = chrs;
		for(int i = 0; i < chrs_.size(); i++){
			filename_anc.push_back( as<std::string>(file_anc) + "_chr" + as<std::string>(chrs_[i]) + ".anc" );
			filename_mut.push_back( as<std::string>(file_mut) + "_chr" + as<std::string>(chrs_[i]) + ".mut" );
			filename_rec.push_back( map_ + "_chr" + as<std::string>(chrs_[i]) + ".txt" );
		}
	}else{
		chrs_ = CharacterVector::create("1");
		filename_anc.push_back( as<std::string>(file_anc) );
		filename_mut.push_back( as<std::string>(file_mut) );
		filename_rec.push_back( map_ );
	}

	double binsize_ = 0.05;
	if(blgsize.isNotNull()){
		binsize_ = as<double>(blgsize);
	}

	bool isM = true;
	if(binsize_ > 100) isM = false;
	if(isM){
		if(map_ == ""){
			Rcpp::stop("Need to specify recombination map when blgsize < 100");
		}
	}

	Sample sample;
	sample.Read(filename_poplabels);

	std::cerr << "Computing pairwise f2 for " << sample.groups.size() + 1 << " groups." << std::endl;
	std::cerr << "blgsize: " << binsize_;
	if(isM) std::cerr << "M" << std::endl;
	if(!isM) std::cerr << "bp" << std::endl;
	std::cerr << std::endl;

	MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
	Muts::iterator it_mut; //iterator for mut file
	float num_bases_tree_persists = 0.0;

	////////// 1. Read one tree at a time /////////

	//as output I want a data frame with columns blockID, chr, start, end, ind1, ind2, pop1, pop2, f2

	std::vector<float> numsnps(1000), chrom(1000), blockBP(1000);
	arma::cube f2_block = arma::zeros<arma::cube>(sample.group_of_haplotype.size(), sample.group_of_haplotype.size(), 1000);

	int blockID = 0;
	float num_infSNPs = 0.0;
	for(int chr = 0; chr < filename_anc.size(); chr++){

		//We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
		//The mut file is read once, file is closed after constructor is called.
		AncMutIterators ancmut(filename_anc[chr], filename_mut[chr]);
		Data data(ancmut.NumTips(), ancmut.NumSnps());

		map recmap(filename_rec[chr].c_str());
		int irec = 0;

		if( blockID+100 >= size(f2_block)[2] ){
			f2_block.resize(sample.group_of_haplotype.size(), sample.group_of_haplotype.size(), size(f2_block)[3] + 1000);
			numsnps.resize(numsnps.size() + 1000);
			chrom.resize(chrom.size() + 1000);
			blockBP.resize(blockBP.size() + 1000);
		}

		//mtr stores the marginal tree, it_mut points to the first SNP at which the marginal tree starts.
		//If marginal tree has no SNPs mapping to it, it_mut points to the first SNP of the following tree that has a SNP mapping to it
		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

		std::vector<float> freqs(2, 0.0), coords(2*data.N-1, 0.0);
		std::vector<Leaves> desc;
		int current_pos = 0;
		double current_rec = 0;
		double rec;
		double factor = 0.0;
		int percentage = 0;
		//iterate through whole file
		std::cerr << "Block: " << blockID << "\r";
		while(num_bases_tree_persists >= 0.0){

			mtr.tree.GetCoordinates(coords);
			mtr.tree.FindAllLeaves(desc);

			if(((*it_mut).tree % (int)(ancmut.NumTrees()/100)) == 0){
				std::cerr << "Block: " << blockID << ", " << "[" << percentage << "%]\r";
				Rcpp::checkUserInterrupt();
				percentage++;
				std::cerr.flush();
			}

			//compute genetic distance at the start of this tree
			if(isM){
				if(recmap.bp[irec] >= (*it_mut).pos){
					if(irec > 0){
						rec = ((double)(*it_mut).pos - recmap.bp[irec-1])/(recmap.bp[irec] - recmap.bp[irec-1]) * (recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) + recmap.gen_pos[irec-1];
					}else{
						rec = ((double)(*it_mut).pos)/(recmap.bp[irec]) * (recmap.gen_pos[irec]);         
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
						rec = recmap.gen_pos[irec]/recmap.bp[irec] * ((*it_mut).pos - recmap.bp[irec]) + recmap.gen_pos[irec];
					}else{
						assert(recmap.bp[irec] >= (*it_mut).pos);
						assert(irec > 0);
						rec = ((double)(*it_mut).pos - recmap.bp[irec-1])/(recmap.bp[irec] - recmap.bp[irec-1]) * (recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) + recmap.gen_pos[irec-1];
					}
				}
				assert(rec >= 0.0);
			}

			//for each tree, calculate f2 stats given sample assignment to pops
			if( (!isM && (*it_mut).pos - current_pos > binsize_) || (isM && rec/100.0 - current_rec > binsize_) ){

				Rcpp::checkUserInterrupt();
				if((int)std::round(num_infSNPs) > 0){
					for(arma::cube::slice_iterator it_c = f2_block.begin_slice(blockID); it_c != f2_block.end_slice(blockID); it_c++){
						*it_c /= factor;
					}
					numsnps[blockID] = std::round(num_infSNPs);
					chrom[blockID] = chr;
					blockBP[blockID] = current_pos;
					num_infSNPs = 0;
					factor = 0.0;
					blockID++;
				}
				if( blockID+100 >= size(f2_block)[2] ){
					f2_block.resize(sample.group_of_haplotype.size(), sample.group_of_haplotype.size(), size(f2_block)[2] + 1000);
					numsnps.resize(numsnps.size() + 1000);
					chrom.resize(chrom.size() + 1000);
					blockBP.resize(blockBP.size() + 1000);
				}
				if(isM){
					while(rec/100.0 - current_rec > binsize_){
						current_rec += binsize_; 
					}
					current_pos = (*it_mut).pos;
				}else{
					while((*it_mut).pos - current_pos > binsize_){
						current_pos += binsize_; 
					}
				}

			}

			//now populate values for this tree
			factor += num_bases_tree_persists;
			int tree = (*it_mut).tree;
			while((*it_mut).tree == tree){
				num_infSNPs++;
				it_mut++;
			}
			traverseTMRCA(mtr, desc, f2_block.slice(blockID), num_bases_tree_persists, coords, t_);

			num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
		}

		if((int)std::round(num_infSNPs) > 0){
			for(arma::cube::slice_iterator it_c = f2_block.begin_slice(blockID); it_c != f2_block.end_slice(blockID); it_c++){
				*it_c /= factor;
			}
			numsnps[blockID] = std::round(num_infSNPs);
			chrom[blockID] = chr;
			blockBP[blockID] = current_pos;
			num_infSNPs = 0;
			blockID++;
		}

	}

	//Output to file
	f2_block.reshape(sample.group_of_haplotype.size(), sample.group_of_haplotype.size(), blockID);
	numsnps.resize(blockID);
	chrom.resize(blockID);
	blockBP.resize(blockID);

	std::ofstream os(as<std::string>(file_out) + ".tmrca");

	os << "blockID\tCHR\tBP\tHAPID1\tHAPID2\tIND1\tIND2\tPOP1\tPOP2\tTMRCA\n";

	//blockID chr bp ind1 ind2 pop1 pop2 TMRCA
	for(int b = 0; b < blockID; b++){
		int k = 0;
		for(arma::cube::slice_iterator it_c = f2_block.begin_slice(b); it_c != f2_block.end_slice(b); it_c++){

			int i = k/f2_block.n_rows;
			int j = k%f2_block.n_rows;

			os << b << "\t" << chrs_[chrom[b]] << "\t" << blockBP[b] << "\t" << i << "\t" << j << "\t" << sample.ind[(int) (i/2.0)] << "\t" << sample.ind[(int) (j/2.0)] << "\t" << 
				sample.groups[sample.group_of_haplotype[i]] << "\t" << sample.groups[sample.group_of_haplotype[j]] << "\t" << *it_c << "\n";

			k++;
		}
	}

	os.close();

	std::cerr << "Block: " << blockID << ", " << "[100%]\r";
	std::cerr << std::endl;

}


//' Function to calculates TMRCA distributions from Relate trees for pairs of populations specified in poplabels.
//'
//' This function will calculate TMRCAs in blocks of prespecified size for all pairs of populations specified in the poplabels file.
//' Please refer to the Relate documentation for input file formats (https://myersgroup.github.io/relate/).
//'
//' @param file_anc Filename of anc file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).
//' @param file_mut Filename of mut file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).
//' @param file_out Filename of output file.
//' @param poplabels Filename of poplabels file
//' @param epochs Vector of epoch boundaries. Should start at 0.
//' @param chrs (Optional) Vector of chromosome IDs
//' @param blgsize (Optional) SNP block size in Morgan. Default is 0.05 (50 cM). If blgsize is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
//' @param file_map (Optional) File prefix of recombination map. Not needed if blgsize is given in base-pairs, i.e. blgsize > 100
//' @return 3d array of dimension #groups x #groups x #blocks. Analogous to output of f2_from_geno in admixtools.
//' @keywords internal
//' @examples
//' file_anc  <- system.file("sim/msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats")
//' file_mut  <- system.file("sim/msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats")
//' poplabels <- system.file("sim/msprime_ad0.8_split250_1.poplabels", package = "twigstats")
//' file_map  <- system.file("sim/genetic_map_combined_b37_chr1.txt.gz", package = "twigstats")
//'
//' #Calculate f2s between all pairs of populations
//' TMRCAdist_from_Relate(file_anc, file_mut, poplabels, file_out = "test", file_map, epochs = c(0,10^seq(3,7,length.out=49)/28))
//' @export
// [[Rcpp::export]]
void TMRCAdist_from_Relate( SEXP file_anc, SEXP file_mut, SEXP poplabels, SEXP file_out, NumericVector epochs, SEXP file_map = R_NilValue, Nullable<CharacterVector> chrs = R_NilValue, Nullable<double> t = R_NilValue, Nullable<double> blgsize = R_NilValue) {

	std::string filename_poplabels = as<std::string>(poplabels);
	std::vector<std::string> filename_anc, filename_mut, filename_rec;
	std::string map_ = "";
	if(file_map != R_NilValue){
		map_ = as<std::string>(file_map);
	}

	std::vector<float> epochs_ = as<std::vector<float> >(epochs);

	double t_ = std::numeric_limits<double>::infinity();
	if(t.isNotNull()){
		t_ = as<double>(t);
	}

	CharacterVector chrs_;
	if(chrs.isNotNull()){
		chrs_ = chrs;
		for(int i = 0; i < chrs_.size(); i++){
			filename_anc.push_back( as<std::string>(file_anc) + "_chr" + as<std::string>(chrs_[i]) + ".anc" );
			filename_mut.push_back( as<std::string>(file_mut) + "_chr" + as<std::string>(chrs_[i]) + ".mut" );
			filename_rec.push_back( map_ + "_chr" + as<std::string>(chrs_[i]) + ".txt" );
		}
	}else{
		chrs_ = CharacterVector::create("1");
		filename_anc.push_back( as<std::string>(file_anc) );
		filename_mut.push_back( as<std::string>(file_mut) );
		filename_rec.push_back( map_ );
	}

	double binsize_ = 0.05;
	if(blgsize.isNotNull()){
		binsize_ = as<double>(blgsize);
	}

	bool isM = true;
	if(binsize_ > 100) isM = false;
	if(isM){
		if(map_ == ""){
			Rcpp::stop("Need to specify recombination map when blgsize < 100");
		}
	}

	Sample sample;
	sample.Read(filename_poplabels);

	std::cerr << "Computing pairwise f2 for " << sample.groups.size() + 1 << " groups." << std::endl;
	std::cerr << "blgsize: " << binsize_;
	if(isM) std::cerr << "M" << std::endl;
	if(!isM) std::cerr << "bp" << std::endl;
	std::cerr << std::endl;

	MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
	Muts::iterator it_mut; //iterator for mut file
	float num_bases_tree_persists = 0.0;

	////////// 1. Read one tree at a time /////////

	//as output I want a data frame with columns blockID, chr, start, end, ind1, ind2, pop1, pop2, f2

	//std::vector<float> numsnps(1000), chrom(1000), blockBP(1000);
	double numsnps, chrom, blockBP;
	int N = sample.group_of_haplotype.size();
	std::vector<std::vector<std::vector<float>>> tmrcas(N);
	for(int i = 0; i < tmrcas.size(); i++){
		tmrcas[i].resize(N);
		for(int j = 0; j < N; j++){
      tmrcas[i][j].resize(epochs_.size());
			std::fill(tmrcas[i][j].begin(), tmrcas[i][j].end(), 0.0);
		}
	}

	//output
	std::ofstream os(as<std::string>(file_out) + ".tmrca");

	os << "blockID\tCHR\tBP\tHAPID1\tHAPID2\tIND1\tIND2\tPOP1\tPOP2\t";
	for(int e = 0; e < epochs_.size(); e++){
		os << epochs_[e] << "\t";
	}
	os << "\n";


	int blockID = 0;
	float num_infSNPs = 0.0;
	for(int chr = 0; chr < filename_anc.size(); chr++){

		//We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
		//The mut file is read once, file is closed after constructor is called.
		AncMutIterators ancmut(filename_anc[chr], filename_mut[chr]);
		Data data(ancmut.NumTips(), ancmut.NumSnps());

		map recmap(filename_rec[chr].c_str());
		int irec = 0;

		//mtr stores the marginal tree, it_mut points to the first SNP at which the marginal tree starts.
		//If marginal tree has no SNPs mapping to it, it_mut points to the first SNP of the following tree that has a SNP mapping to it
		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

		std::vector<float> freqs(2, 0.0), coords(2*data.N-1, 0.0);
		std::vector<Leaves> desc;
		int current_pos = 0;
		double current_rec = 0;
		double rec;
		double factor = 0.0;
		int percentage = 0;
		//iterate through whole file
		std::cerr << "Block: " << blockID << "\r";
		while(num_bases_tree_persists >= 0.0){

			mtr.tree.GetCoordinates(coords);
			mtr.tree.FindAllLeaves(desc);

			if(((*it_mut).tree % (int)(ancmut.NumTrees()/100)) == 0){
				std::cerr << "Block: " << blockID << ", " << "[" << percentage << "%]\r";
				Rcpp::checkUserInterrupt();
				percentage++;
				std::cerr.flush();
			}

			//compute genetic distance at the start of this tree
			if(isM){
				if(recmap.bp[irec] >= (*it_mut).pos){
					if(irec > 0){
						rec = ((double)(*it_mut).pos - recmap.bp[irec-1])/(recmap.bp[irec] - recmap.bp[irec-1]) * (recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) + recmap.gen_pos[irec-1];
					}else{
						rec = ((double)(*it_mut).pos)/(recmap.bp[irec]) * (recmap.gen_pos[irec]);         
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
						rec = recmap.gen_pos[irec]/recmap.bp[irec] * ((*it_mut).pos - recmap.bp[irec]) + recmap.gen_pos[irec];
					}else{
						assert(recmap.bp[irec] >= (*it_mut).pos);
						assert(irec > 0);
						rec = ((double)(*it_mut).pos - recmap.bp[irec-1])/(recmap.bp[irec] - recmap.bp[irec-1]) * (recmap.gen_pos[irec] - recmap.gen_pos[irec-1]) + recmap.gen_pos[irec-1];
					}
				}
				assert(rec >= 0.0);
			}

			//for each tree, calculate f2 stats given sample assignment to pops
			if( (!isM && (*it_mut).pos - current_pos > binsize_) || (isM && rec/100.0 - current_rec > binsize_) ){

				Rcpp::checkUserInterrupt();
				if((int)std::round(num_infSNPs) > 0){
					numsnps = std::round(num_infSNPs);
					chrom = chr;
					blockBP = current_pos;
					num_infSNPs = 0;
					factor = 0.0;
					for(int i = 0; i < N; i++){
						for(int j = 0; j < N; j++){
							os << blockID << "\t" << chrs_[chrom] << "\t" << blockBP << "\t" << i << "\t" << j << "\t" << sample.ind[(int) (i/2.0)] << "\t" << sample.ind[(int) (j/2.0)] << "\t" <<	sample.groups[sample.group_of_haplotype[i]] << "\t" << sample.groups[sample.group_of_haplotype[j]] << "\t";
							for(int e = 0; e < epochs_.size(); e++){
							  os << tmrcas[i][j][e] << "\t";
							}
							os << "\n";
						}
					}

					for(int i = 0; i < tmrcas.size(); i++){
						for(int j = 0; j < N; j++){
							std::fill(tmrcas[i][j].begin(), tmrcas[i][j].end(), 0.0);
						}
					}

					blockID++;
				}

				if(isM){
					while(rec/100.0 - current_rec > binsize_){
						current_rec += binsize_; 
					}
					current_pos = (*it_mut).pos;
				}else{
					while((*it_mut).pos - current_pos > binsize_){
						current_pos += binsize_; 
					}
				}

			}

			//now populate values for this tree
			factor += num_bases_tree_persists;
			int tree = (*it_mut).tree;
			while((*it_mut).tree == tree){
				num_infSNPs++;
				it_mut++;
			}
			traverseTMRCAdist(mtr, desc, tmrcas, epochs_, num_bases_tree_persists, coords);

			num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
		}

		if((int)std::round(num_infSNPs) > 0){
			numsnps = std::round(num_infSNPs);
			chrom = chr;
			blockBP = current_pos;
			num_infSNPs = 0;

			for(int i = 0; i < N; i++){
				for(int j = 0; j < N; j++){
					os << blockID << "\t" << chrs_[chrom] << "\t" << blockBP << "\t" << i << "\t" << j << "\t" << sample.ind[(int) (i/2.0)] << "\t" << sample.ind[(int) (j/2.0)] << "\t" <<	sample.groups[sample.group_of_haplotype[i]] << "\t" << sample.groups[sample.group_of_haplotype[j]] << "\t";
					for(int e = 0; e < epochs_.size(); e++){
						os << tmrcas[i][j][e] << "\t";
					}
					os << "\n";
				}
			}

			for(int i = 0; i < tmrcas.size(); i++){
				for(int j = 0; j < N; j++){
					std::fill(tmrcas[i][j].begin(), tmrcas[i][j].end(), 0.0);
				}
			}

			blockID++;
		}

	}

	os.close();

	std::cerr << "Block: " << blockID << ", " << "[100%]\r";
	std::cerr << std::endl;

}

