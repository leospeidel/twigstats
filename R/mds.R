library(dplyr)
library(tidyr)

#' Compute f2 stats genome-wide
#'
#' This function aggregates f2_blocks to compute a point estimate of f2 genome-wide
#'
#' @param f2_blocks A 3d array of blocked f2 statistics
#' @return Returns NxN matrix storing f2 statistics between pairs of populations
#' @keywords internal
calcf2 <- function(f2_blocks){

	N <- ncol(f2_blocks)
	L <- length(f2_blocks[1,1,])

	weights <- as.numeric(gsub(names(f2_blocks[1,1,]), pattern = "l", replace = ""))

	f2_blocks_2d <- matrix(f2_blocks, nrow = N * N, ncol = L)
	weight_2d <- matrix(weights,ncol = L, nrow = N*N, byrow = T)

	weighted_values <- f2_blocks_2d * weight_2d
	# Compute the sum along the L dimension
	sum_weights <- rowSums(weighted_values)
	# Compute the sum of weights
	total_weight <- sum(weights)
	# Calculate the final result for each pair i, j
	result <- sum_weights / total_weight
	# Reshape the result back to a matrix (N x N)
	result_matrix <- matrix(result, nrow = N, ncol = N)

	return(result_matrix)

}



#' Function that computes a PCA and MDS from f2_blocks
#'
#' This function takes f2_blocks as input, computes outgroup f3 statistics, and then computes PCA and MDS from the f3 statistics.
#'
#' @rdname mds
#' @param f2_blocks A 3d array of blocked f2 statistics
#' @param poplabels Filename of poplabels file
#' @param outgroup Name of outgroup population
#' @return Returns a data frame storing the PCA and MDS results.
#' @examples
#' #These lines assign file names to variables file_anc, file_mut, poplabels, file_map.
#' #see https://myersgroup.github.io/relate/getting_started.html#Output for file formats
#' file_anc  <- system.file("sim/msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats")
#' file_mut  <- system.file("sim/msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats")
#' #see https://myersgroup.github.io/relate/input_data.html for file formats
#' poplabels <- system.file("sim/msprime_ad0.8_split250_1.ind.poplabels", package = "twigstats")
#' file_map       <- system.file("sim/genetic_map_combined_b37_chr1.txt.gz", package = "twigstats") #recombination map (three column format)
#'
#' f2_blocks <- f2_blocks_from_Relate(file_anc = file_anc, file_mut = file_mut, poplabels = poplabels, file_map = file_map, t = 1000)
#' df <- calc_mds(f2_blocks, poplabels, "P4")
#' print(head(df))
#'
#' @export 
calc_mds <- function(f2_blocks, poplabels, outgroup){

	#load(file)
	d <- dim(f2_blocks)[1:2]
	subn <- colnames(f2_blocks)
	subn <- subn[subn != "Root"]

	poplabels <- read.table(poplabels, header = T)
	poplabels <- subset(poplabels, ID %in% subn)
	subn <- unique(as.matrix(poplabels$ID))

	n <- colnames(f2_blocks)
	oname <- outgroup
	oname <- oname[oname %in% n]
	subn <- intersect(subn,n)

	f2_blocks <- f2_blocks[c(oname,subn),c(oname,subn),]

	A <- calcf2(f2_blocks)
	colnames(A) <- colnames(f2_blocks)
	rownames(A) <- colnames(f2_blocks)
	
  #compute f3 from f2
	N <- ncol(A)
	Ak_i <- matrix(A[oname, ], nrow = N, ncol = N, byrow = F)
	Ak_j <- matrix(A[oname, ], nrow = N, ncol = N, byrow = T)
	A_i_j <- A
	result <- 0.5*(Ak_i + Ak_j - A_i_j)

	result <- cbind(result, pop2 = rownames(result))
	result <- as.data.frame(result)
	result %>% pivot_longer(cols = !pop2,names_to = "pop3", values_to = "est") -> f3_matrix
	f3_matrix     <- subset(f3_matrix, pop2 != oname & pop3 != oname)
	f3_matrix$est <- as.numeric(as.matrix(f3_matrix$est))

	f3_matrix %>% pivot_wider(id_cols = pop2, names_from = pop3, values_from = est) -> df_wider
	colnames(df_wider)[1] <- "pop1"

	names <- unique(df_wider$pop1)
	cols <- names(df_wider[1,])
	df_wider <- df_wider[order(df_wider$pop1),]
	n <- as.matrix(df_wider$pop1)
	df_wider <- df_wider[, c("pop1", n[n %in% subn])]

	df_ret <- data.frame()

	pca <- df_wider %>% 
		select(where(is.numeric)) %>% # retain only numeric columns
		prcomp(scale = T) 
	poplabels$POP <- factor(poplabels$POP, levels = n)
	poplabels <- poplabels[order(poplabels$POP),]
	df_ret <- rbind(df_ret, cbind(data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], 
																					 ID = as.matrix(poplabels$POP[!is.na(poplabels$POP)])), poplabels[,-1], method = "PCA"))

	cmd_fit <- cmdscale(1-df_wider[,-1], eig = T, k = 2)
	poplabels$POP <- factor(poplabels$POP, levels = n)
	poplabels <- poplabels[order(poplabels$POP),]
	df_ret <- rbind(df_ret, cbind(data.frame(PC1 = cmd_fit$points[,1], PC2 = cmd_fit$points[,2], PC3 = NA, 
																					 ID = as.matrix(poplabels$POP)), poplabels[,-1], method = "MDS"))

	return(df_ret)

}






