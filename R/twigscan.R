library(dplyr)

#' Function implementing TwigScan
#'
#' This function computes f2s or Fst in blocks along the genome and returns a data frame with columns blockID, pos, pop1, pop2, and the f2 or Fst values.
#'
#' @rdname TwigScan
#' @param file_anc Filename of anc file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).
#' @param file_mut Filename of mut file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).
#' @param file_map File prefix of recombination map.
#' @param file_out File prefix of output files
#' @param poplabels Filename of poplabels file
#' @param blgsize (Optional) SNP block size in Morgan. Default is 0.05 (5 cM). If blgsize is 1 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
#' @param t (Optional) Time cutoff in generations. Default: Inf
#' @param use_muts (Optional) Calculate traditional f2 statistics by only using mutations mapped to Relate trees. Default: False.
#' @param Fst (Optional) If TRUE, compute Fst. Default: FALSE
#' @return Returns a data frame with Fst or f2 values.
#' @examples
#' #These lines assign file names to variables file_anc, file_mut, poplabels, file_map.
#' #see https://myersgroup.github.io/relate/getting_started.html#Output for file formats
#' file_anc  <- system.file("sim/msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats")
#' file_mut  <- system.file("sim/msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats")
#' #see https://myersgroup.github.io/relate/input_data.html for file formats
#' poplabels <- system.file("sim/msprime_ad0.8_split250_1.poplabels", package = "twigstats")
#' file_map  <- system.file("sim/genetic_map_combined_b37_chr1.txt.gz", package = "twigstats") #recombination map (three column format)
#'
#' df <- TwigScan(file_anc = file_anc,
#'			   file_mut = file_mut,
#'			   poplabels = poplabels,
#'			   file_map = file_map,
#'			   file_out = "test",
#'			   blgsize = 10000,  #optional
#'			   use_muts = FALSE, #optional
#'			   t = 1000,	       #optional
#'			   Fst = TRUE	       #optional
#'	 )
#'
#' print(head(df))
#' @export 
TwigScan <- function(file_anc, file_mut, poplabels, file_map, file_out, blgsize = 10000, t = Inf, use_muts = FALSE, Fst = FALSE){

	if(Fst){

		f2_blocks <- Fst_blocks_from_Relate(file_anc = file_anc,
																				file_mut = file_mut,
																				poplabels = poplabels,
																				blgsize = blgsize,
																				file_map = file_map,
																				dump_blockpos = paste0(file_out,'_t',t,'.pos'),
																				use_muts = use_muts,
																				t = t)

	}else{

		f2_blocks <- f2_blocks_from_Relate(file_anc = file_anc,
																			 file_mut = file_mut,
																			 poplabels = poplabels,
																			 blgsize = blgsize,
																			 file_map = file_map,
																			 dump_blockpos = paste0(file_out,'_t',t,'.pos'),
																			 use_muts = use_muts,
																			 t = t)

	}

	pos <- read.table(paste0(file_out,'_t',t,'.pos'))
	names <- dimnames(f2_blocks)[[1]]

	df <- data.frame()
	for(i in 1:length(names)){
		for(j in i:length(names)){
			if(i != j){
				df <- rbind(df, data.frame(block = pos[,1], pos = pos[,2], val = f2_blocks[i,j,], pop1 = names[i], pop2 = names[j], cutoff = t, use_muts = use_muts))
			}
		}
	}

	df %>% filter(val > -Inf, val < Inf) %>% group_by(pop1,pop2) %>% mutate(zval = (val - median(val))/mad(val)) -> df

	return(df)

}



