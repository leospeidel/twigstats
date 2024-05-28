library(dplyr)

#' Function to compute genome-wide copying proportions.
#'
#' This function takes the output of the function Painting and computes the genome-wide 'copying vectors', 
#' i.e. the proportion of the genome copied from each other reference population..
#'
#' @param filename_painting Vector containing filenames of painting profiles. Output of Painting.
#' @param filename_idfile Filename of idfile. Output of Painting.
#' @param nboot Number of bootstrap samples.
#' @param blocksize Number of blocks to combine for the bootstrap. E.g. if Painting was run with a blgsize of 1e-5 Morgans, blocksize should be 5000 to achieve a blocksize of 5cM.
#' @param use_IDs If TRUE, compute profile for each sample. If FALSE (default), compute profile for each group as specified in the second column of the idfile.
#' @return Returns a data frame with copying proportions per bootstrap sample.
#' @examples
#' #This path stores files precomputed using Painting().
#' path      <- paste0(system.file("test/", package = "twigstats"),"/")
#' prefix    <- "test" #prefix of files under path
#'
#' #compute the painting profiles with 10 bootstrap samples and a blocksize of 5cM (5000*0.001)
#' df <- PaintingProfile(c(paste0(path,prefix,"_painting.txt.gz")), paste0(path,prefix,"_idfile.txt.gz"), nboot = 10, blocksize = 5000)
#' head(df)
#' @export
PaintingProfile <- function(filename_painting, filename_idfile, nboot, blocksize, use_IDs = FALSE){

	if (nboot <= 0) {
		stop("Value of nboot needs to be at least 1.")
	}
	if (blocksize <= 0) {
		stop("Value of blocksize needs to be at least 1.")
	}
  idcol <- 2
	if(use_IDs){
    idcol <- 1
	}

	print_progress_bar <- function(iteration, total, length = 50) {
		percent <- round(iteration / total * 100)
		progress <- round(length * iteration / total)
		bar <- paste0(c(rep("=", progress), rep(" ", length - progress)), collapse = "")
		cat("\r[", bar, "] ", percent, "%", sep = "")
		flush.console()
	}

	ids <- read.table(filename_idfile)
	ids <- ids[sort(rep(1:nrow(ids),2)),]

	df <- data.frame()
	chr <- 1
	print_progress_bar(0, length(filename_painting))
	for(file in filename_painting){

		lines <- readLines(file)[-c(1:2)]
		lines <- lines[seq(1,length(lines),2)]
		df_chr <- as.data.frame(do.call(rbind,lapply(strsplit(lines, " "), as.numeric)))
		#df_chr <- read.table(pipe(paste0("tail -n +2", file, "| sed -n 'n;p'")))
		#df_chr <- read.table(pipe(paste0("sed -n '3~2p' ", file)))

		df_chr <- df_chr[,-1]
		df_chr <- df_chr[order(ids[,2]),]
		n <- ncol(df_chr)
		df_chr <- ids[t(as.matrix(df_chr)),2] 

		#POP is the ID of the population
		#labs is the group I copied from
		df_chr <- data.frame(POP = sort(rep(ids[,idcol], n)), CHR = chr, pos_id = rep(1:n,nrow(ids)), labs = df_chr)
		df_chr$block <- cut(df_chr$pos_id, breaks = seq(1,n,blocksize))
		df_chr %>% group_by(POP, labs, CHR, block) %>% summarize(count = length(POP)) -> df_chr
		df <- bind_rows(df, df_chr)
		rm(df_chr)
		print_progress_bar(chr, length(filename_painting))
		chr <- chr + 1
	}
	cat("\n")

	#now do a single bootstrap on df
	if(nboot > 1){
		df$block <- paste0(df$CHR, "-", df$block)
		df_sum <- data.frame()
		for(b in 1:nboot){
			df_sum <- bind_rows(df_sum, df %>% sample_n(size = nrow(df), replace = TRUE) %>% group_by(POP, labs) %>% summarize(prop = sum(count)) %>% group_by(POP) %>% mutate(prop = prop/sum(prop), bootstrap = b))
		}
	}else{
		df %>% group_by(POP, labs) %>% summarize(prop = sum(count)) %>% group_by(POP) %>% mutate(prop = prop/sum(prop)) -> df_sum
	}

	return( df_sum )
}

#' Function to run an NNLS on genome-wide copying proportions to estimate admixture proportions.
#'
#' This function computes admixture proprtions using a non-negative least squares on painting profiles obtained from the function PaintingProfile
#'
#' @param df_sum Output of PaintingProfile
#' @param target_pops Vector of populations to be fitted
#' @param source_pops (Optional) Vector of putative source populations. If not provided, all remaining populations are used.
#' @return Returns a data frame with admixture proportions.
#' @examples
#' library(dplyr)
#'
#' #This path stores files precomputed using Painting().
#' path      <- paste0(system.file("test/", package = "twigstats"),"/")
#' prefix    <- "test" #prefix of files under path
#'
#' #compute the painting profiles with 100 bootstrap samples and a blocksize of 5cM (5000*0.001)
#' df <- PaintingProfile(c(paste0(path,prefix,"_painting.txt.gz")), paste0(path,prefix,"_idfile.txt.gz"), nboot = 100, blocksize = 5000)
#' 
#' #compute the NNLS to get admixture proportions
#' df <- PaintingNNLS(df, target_pops = c("PX"), source_pops = c("P2","P3"))
#'
#' #Now you can summarize the bootstrap samples 
#' df %>% group_by(target, POP) %>% summarize(mean_ancestry = mean(ancestry), sd_ancestry = sd(ancestry)) -> df
#' print(df)
#' @export
PaintingNNLS <- function(df_sum, target_pops, source_pops = NULL){

	res <- data.frame()
	if("bootstrap" %in% colnames(df_sum)){
	  for(b in unique(df_sum$bootstrap)){
			df_sum %>% filter(bootstrap == b) %>% select(-bootstrap) %>% tidyr::spread(key = "labs", value = "prop") -> mat

      mat[is.na(mat)] <- 0
			if(is.null(source_pops)){
				mat %>% filter(!POP %in% target_pops) -> sources
			}else{
				mat %>% filter(POP %in% source_pops) -> sources
			}
			mat %>% filter(POP %in% target_pops)  -> targets

			pops <- as.matrix(sources[,1])
			tpops <- as.matrix(targets[,1])
			sources <- t(as.matrix(sources[,-1]))
			targets <- t(as.matrix(targets[,-1]))

			for(i in 1:length(tpops)){
				res <- rbind(res, data.frame(bootstrap = b, target = tpops[i], POP = pops, ancestry = lsei::pnnls(sources, targets[,i], sum = 1)$x))
			}
		}
	}else{
		df_sum %>% tidyr::spread(key = "labs", value = "prop") -> mat

		if(is.null(source_pops)){
			mat %>% filter(!POP %in% target_pops) -> sources
		}else{
			mat %>% filter(POP %in% source_pops) -> sources
		}
		mat %>% filter(POP %in% target_pops)  -> targets

		pops <- as.matrix(sources[,1])
		tpops <- as.matrix(targets[,1])
		sources <- t(as.matrix(sources[,-1]))
		targets <- t(as.matrix(targets[,-1]))

		for(i in 1:length(tpops)){
			res <- rbind(res, data.frame(target = tpops[i], POP = pops, ancestry = lsei::pnnls(sources, targets[,i], sum = 1)$x))
		}
	}

	return(res)

}
