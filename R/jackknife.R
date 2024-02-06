library(dplyr)

#' Function implementing the F4-ratio statistic
#'
#' This function computes the admixture proportion given five populations. A population history following ((PI,P1),P2,PO) is assumed,
#' and the target is assumed to be a mixture of proximal sources P1 and P2, i.e. PX = alpha*P2 + (1-alpha)*P1 
#'
#' @rdname f-statistics
#' @param f2_blocks A 3d array of blocked f2 statistics
#' @param popO Name of outgroup population
#' @param popI Name of ingroup population
#' @param pop1 Name of source that clusters with ingroup
#' @param pop2 Name of other source
#' @param popX Name of target group.
#' @return Returns a data frame with admixture proportion estimates and jacknifed standard errors.
#' @examples
#' #These lines assign file names to variables file_anc, file_mut, poplabels, file_map.
#' #see https://myersgroup.github.io/relate/getting_started.html#Output for file formats
#' file_anc  <- system.file("sim/msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats")
#' file_mut  <- system.file("sim/msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats")
#' #see https://myersgroup.github.io/relate/input_data.html for file formats
#' poplabels <- system.file("sim/msprime_ad0.8_split250_1.poplabels", package = "twigstats")
#' file_map  <- system.file("sim/genetic_map_combined_b37_chr1.txt", package = "twigstats") #recombination map (three column format)
#' 
#' #Calculate regular f2s between all pairs of populations
#' f2_blocks1 <- f2_blocks_from_Relate(file_anc, file_mut, poplabels, file_map)
#' f4_ratio(f2_blocks1, popX="PX", popI="P1", pop1="P2", pop2="P3", popO="P4")
#' 
#' #Use a twigstats cutoff of 500 generations
#' f2_blocks2 <- f2_blocks_from_Relate(file_anc, file_mut, poplabels, file_map, t = 500)
#' f4_ratio(f2_blocks2, popX="PX", popI="P1", pop1="P2", pop2="P3", popO="P4") 
#' @export 
f4_ratio <- function(f2_blocks, popO,popI,pop1,pop2,popX, mode = 1){

  if(mode == 1){
    #f4(popO,popI,popX,pop1)/f4(popO,popI,pop2,pop1)
    df_jack <- data.frame(blockID = 1:dim(f2_blocks)[3],
                          hj      = as.numeric(gsub(dimnames(f2_blocks)[3][[1]], pattern = "l", replace = "")),
                          numj    = f2_blocks[popO,pop1,] + f2_blocks[popI,popX,] - f2_blocks[popO,popX,] - f2_blocks[popI,pop1,],
                          denomj  = f2_blocks[popO,pop1,] + f2_blocks[popI,pop2,] - f2_blocks[popO,pop2,] - f2_blocks[popI,pop1,]
    )

    num        <- sum(df_jack$hj * df_jack$numj)/sum(df_jack$hj)
    denom      <- sum(df_jack$hj * df_jack$denomj)/sum(df_jack$hj)
    df_jack$Dj <- (sum(df_jack$hj) * num - df_jack$hj * df_jack$numj)/(sum(df_jack$hj) * denom - df_jack$hj * df_jack$denomj)
    df_jack$hj <- sum(df_jack$hj)/df_jack$hj
  }else{
    #f4(popO,popI,popX,pop2)/f4(popO,popI,pop1,pop2)
    df_jack <- data.frame(blockID = 1:dim(f2_blocks)[3],
                          hj      = as.numeric(gsub(dimnames(f2_blocks)[3][[1]], pattern = "l", replace = "")),
                          numj    = f2_blocks[popO,pop2,] + f2_blocks[popI,popX,] - f2_blocks[popO,popX,] - f2_blocks[popI,pop2,],
                          denomj  = f2_blocks[popO,pop2,] + f2_blocks[popI,pop1,] - f2_blocks[popO,pop1,] - f2_blocks[popI,pop2,]
                          )

    num        <- sum(df_jack$hj * df_jack$numj)/sum(df_jack$hj)
    denom      <- sum(df_jack$hj * df_jack$denomj)/sum(df_jack$hj)
    df_jack$Dj <- (sum(df_jack$hj) * num - df_jack$hj * df_jack$numj)/(sum(df_jack$hj) * denom - df_jack$hj * df_jack$denomj)
    df_jack$hj <- sum(df_jack$hj)/df_jack$hj
  }

	return( cbind(popO = popO, popI = popI, pop1 = pop1, pop2 = pop2, popX = popX, jackknife(df_jack)) )

}


#' Function implementing the block jackknife
#'
#' This function implements a jackknife on the input data.
#'
#' @rdname f-statistics
#' @param data frame. Three columns called blockID, hj, Dj, storing block ID, weight of block, and statistic without that block
#' @return Returns a data table with columns val and se.
#' @export
jackknife <- function(df_jack){
  D <- mean(df_jack$Dj)

  df_jack %>% filter(!is.infinite(hj)) %>%
  mutate(pseudoj = hj*D - (hj-1)*Dj) %>%
  summarize(g = length(Dj), DJ = sum(D - Dj) + sum(Dj/hj), se2 = mean( (pseudoj - DJ)^2/(hj-1) ) ) -> df_jack

  return( data.frame(val = D, se = sqrt(df_jack$se2)) )
}

