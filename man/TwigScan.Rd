% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twigscan.R
\name{TwigScan}
\alias{TwigScan}
\title{Function implementing TwigScan}
\usage{
TwigScan(
  file_anc,
  file_mut,
  poplabels,
  file_map,
  file_out,
  blgsize = 10000,
  t = Inf,
  use_muts = FALSE,
  Fst = FALSE
)
}
\arguments{
\item{file_anc}{Filename of anc file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).}

\item{file_mut}{Filename of mut file. If chrs is specified, this should only be the prefix, resulting in filenames of $\{file_anc\}_chr$\{chr\}.anc(.gz).}

\item{poplabels}{Filename of poplabels file}

\item{file_map}{File prefix of recombination map.}

\item{file_out}{File prefix of output files}

\item{blgsize}{(Optional) SNP block size in Morgan. Default is 0.05 (5 cM). If blgsize is 1 or greater, if will be interpreted as base pair distance rather than centimorgan distance.}

\item{t}{(Optional) Time cutoff in generations. Default: Inf}

\item{use_muts}{(Optional) Calculate traditional f2 statistics by only using mutations mapped to Relate trees. Default: False.}

\item{Fst}{(Optional) If TRUE, compute Fst. Default: FALSE}
}
\value{
Returns a data frame with Fst or f2 values.
}
\description{
This function computes f2s or Fst in blocks along the genome and returns a data frame with columns blockID, pos, pop1, pop2, and the f2 or Fst values.
}
\examples{
#These lines assign file names to variables file_anc, file_mut, poplabels, file_map.
#see https://myersgroup.github.io/relate/getting_started.html#Output for file formats
file_anc  <- system.file("sim/msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats")
file_mut  <- system.file("sim/msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats")
#see https://myersgroup.github.io/relate/input_data.html for file formats
poplabels <- system.file("sim/msprime_ad0.8_split250_1.poplabels", package = "twigstats")
file_map  <- system.file("sim/genetic_map_combined_b37_chr1.txt.gz", package = "twigstats") #recombination map (three column format)

df <- TwigScan(file_anc = file_anc,
		   file_mut = file_mut,
		   poplabels = poplabels,
		   file_map = file_map,
		   file_out = "test",
		   blgsize = 10000,  #optional
		   use_muts = FALSE, #optional
		   t = 1000,	       #optional
		   Fst = TRUE	       #optional
 )

print(head(df))
}
