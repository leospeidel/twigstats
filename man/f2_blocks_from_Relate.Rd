% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{f2_blocks_from_Relate}
\alias{f2_blocks_from_Relate}
\title{Function to calculate f2 statistics from Relate trees for pairs of populations specified in poplabels.}
\usage{
f2_blocks_from_Relate(
  file_anc,
  file_mut,
  poplabels,
  file_map = NULL,
  chrs = NULL,
  blgsize = NULL,
  mu = NULL,
  tmin = NULL,
  t = NULL,
  transitions = NULL,
  use_muts = NULL,
  minMAF = NULL,
  dump_blockpos = NULL,
  apply_corr = NULL
)
}
\arguments{
\item{file_anc}{Filename of anc file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).}

\item{file_mut}{Filename of mut file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).}

\item{poplabels}{Filename of poplabels file}

\item{file_map}{(Optional) File prefix of recombination map. Not needed if blgsize is given in base-pairs, i.e. blgsize > 100}

\item{chrs}{(Optional) Vector of chromosome IDs}

\item{blgsize}{(Optional) SNP block size in Morgan. Default is 0.05 (5 cM). If blgsize is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance. If blgsize is negative, every tree is its own block.}

\item{mu}{(Optional) Per base per generation mutation rate to scale f2 values. Default: 1.25e-8}

\item{tmin}{(Optional) Minimum time cutof in generations. Any lineages younger than tmin will be excluded from the analysis. Default: t = 0.}

\item{t}{(Optional) Time cutoff in generations. Default: Inf}

\item{transitions}{(Optional) Set this to FALSE to exclude transition SNPs. Only meaningful with use_muts}

\item{use_muts}{(Optional) Calculate traditional f2 statistics by only using mutations mapped to Relate trees. Default: false.}

\item{minMAF}{(Optional) Minimum frequency cutoff. Default: 1 (i.e. excl singletons)}

\item{dump_blockpos}{(Optional) Filename of blockpos file.}

\item{apply_corr}{(Optional) Use small sample size correction. Default: true.}
}
\value{
3d array of dimension #groups x #groups x #blocks. Analogous to output of f2_from_geno in admixtools.
}
\description{
This function will calculate f2 statistics in blocks of prespecified size for all pairs of populations specified in the poplabels file.
Please refer to the Relate documentation for input file formats (https://myersgroup.github.io/relate/).
The output is in a format that is directly accepted by the admixtools R package to calculate
f3, f4, f4ratio, D statistics and more (https://uqrmaie1.github.io/admixtools/).
}
\examples{
file_anc  <- system.file("sim/msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats")
file_mut  <- system.file("sim/msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats")
poplabels <- system.file("sim/msprime_ad0.8_split250_1.poplabels", package = "twigstats")
file_map  <- system.file("sim/genetic_map_combined_b37_chr1.txt.gz", package = "twigstats")

#Calculate f2s between all pairs of populations
f2_blocks <- f2_blocks_from_Relate(file_anc, file_mut, poplabels, file_map)
f4_ratio(f2_blocks, popX="PX", popI="P1", pop1="P2", pop2="P3", popO="P4")

#Use a cutoff of 500 generations
f2_blocks <- f2_blocks_from_Relate(file_anc, file_mut, poplabels, file_map, t = 500)
f4_ratio(f2_blocks, popX="PX", popI="P1", pop1="P2", pop2="P3", popO="P4")
}
