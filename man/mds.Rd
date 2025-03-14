% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mds.R
\name{calc_mds}
\alias{calc_mds}
\title{Function that computes a PCA and MDS from f2_blocks}
\usage{
calc_mds(f2_blocks, poplabels, outgroup)
}
\arguments{
\item{f2_blocks}{A 3d array of blocked f2 statistics}

\item{poplabels}{Filename of poplabels file}

\item{outgroup}{Name of outgroup population}
}
\value{
Returns a data frame storing the PCA and MDS results.
}
\description{
This function takes f2_blocks as input, computes outgroup f3 statistics, and then computes PCA and MDS from the f3 statistics.
}
\examples{
#These lines assign file names to variables file_anc, file_mut, poplabels, file_map.
#see https://myersgroup.github.io/relate/getting_started.html#Output for file formats
file_anc  <- system.file("sim/msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats")
file_mut  <- system.file("sim/msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats")
#see https://myersgroup.github.io/relate/input_data.html for file formats
poplabels <- system.file("sim/msprime_ad0.8_split250_1.ind.poplabels", package = "twigstats")
file_map       <- system.file("sim/genetic_map_combined_b37_chr1.txt.gz", package = "twigstats") #recombination map (three column format)

f2_blocks <- f2_blocks_from_Relate(file_anc = file_anc, file_mut = file_mut, poplabels = poplabels, file_map = file_map, t = 1000)
df <- calc_mds(f2_blocks, poplabels, "P4")
print(head(df))

}
