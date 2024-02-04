# twigstats

Boost f-statitics power using genealogies. Compatible with [admixtools2](https://uqrmaie1.github.io/admixtools/index.html)

## Installation

You can install this package by running the following command in R:
```R
library(devtools)
install_github("leospeidel/twigstats")
```

Alternatively, clone this directory (https://github.com/leospeidel/twigstats) and then in R type
```R
library(devtools)
install()
```

<br/>

## Basic Usage

Twigstats computes f2 statistics for pairs of populations. The output can be directly fed into functions of the admixtools2 R package (https://uqrmaie1.github.io/admixtools/index.html) to compute any other f-statistic.<br>

We will use an example stored under ```inst/sim/```, which contains simulated data for five populations as described in our paper (Figure 1).

### Computing f2-statistics from Relate genealogies
We can compute f2-statistics on Relate genealogies (https://myersgroup.github.io/relate/) inferred for our samples of interest.<br>

The central function is called ```f2_blocks_from_Relate()``` which takes Relate output files as input.<br>
Please type ```?f2_blocks_from_Relate``` for more details:<br>

```R
library(twigstats)

#These lines assign file names to variables file_anc, file_mut, poplabels, file_map.
#see https://myersgroup.github.io/relate/getting_started.html#Output for file formats
file_anc  <- system.file("sim/msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats")
file_mut  <- system.file("sim/msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats")
#see https://myersgroup.github.io/relate/input_data.html for file formats
poplabels <- system.file("sim/msprime_ad0.8_split250_1.poplabels", package = "twigstats")
file_map  <- system.file("sim/genetic_map_combined_b37_chr1.txt", package = "twigstats") #recombination map (three column format)

#Calculate regular f2s between all pairs of populations
f2_blocks1 <- f2_blocks_from_Relate(file_anc, file_mut, poplabels, file_map)
f4_ratio(f2_blocks1, popX="PX", popI="P1", pop1="P2", pop2="P3", popO="P4")

#Now calculate f2s using a twigstats cutoff of 500 generations. 
#This should give us a big power boost.
f2_blocks2 <- f2_blocks_from_Relate(file_anc, file_mut, poplabels, file_map, t = 500)
f4_ratio(f2_blocks2, popX="PX", popI="P1", pop1="P2", pop2="P3", popO="P4")
```
- Use argument <b>use_muts</b> to compute f2-statistics on the (age-ascertained) mutations, instead of branch lengths of Relate genealogies.
- Use argument <b>blgsize</b> to change jackknife block size
- Set argument <b>transitions</b> to False to exclude transitions

#### Expected output
In the above example, PX is admixed between P2 and P3. The f4-ratio statistic computes an estimate of this admixture proportion.<br>
The expected output from running the code above is
```
#without twigstats ascertainment
  popO popI pop1 pop2 popX      val        se
1   P4   P1   P2   P3   PX 1.572514 0.3871271

#with twigtstats ascertainment of 500 generations
  popO popI pop1 pop2 popX       val         se
1   P4   P1   P2   P3   PX 0.7602321 0.02410571
```
So while we are unable to get a reliable estimate without twigstats ascertainment, we estimate PX to carry approximately 80% P2 and 20% P3 ancestry with twigstats.

#### Using admixtools2
You can now directly feed this into admixtools2 functions:
```R
library(admixtools)
#no twigstats ascertainment
f4(f2_blocks1, "P4", "P3", "P2", "PX")
#with twigstats ascertainment
f4(f2_blocks2, "P4", "P3", "P2", "PX")
```
The expected outcome of running this code is
```
# A tibble: 1 × 8
  pop1  pop2  pop3  pop4        est       se      z     p
  <chr> <chr> <chr> <chr>     <dbl>    <dbl>  <dbl> <dbl>
1 P4    P3    P2    PX    -0.000670 0.000694 -0.965 0.335
# A tibble: 1 × 8
  pop1  pop2  pop3  pop4       est        se     z         p
  <chr> <chr> <chr> <chr>    <dbl>     <dbl> <dbl>     <dbl>
1 P4    P3    P2    PX    0.000572 0.0000230  24.8 4.71e-136
```
<br/>

### Computing f2-statistics on age-ascertained SNPs
Sometimes, we are prefer not to graft samples into genealogies (e.g. due to low sequencing coverage, or if we don't trust the imputation). <br>
In this case, we can still compute f2-statistics on age ascertained mutations.<br>
The central function is called ```f2_blocks_from_RelateAges()``` which takes plink (bed/bim/fam) files and mutation ages (.mut) as input.<br>
Please type ```?f2_blocks_from_RelateAges``` for more details:<br>

```R
library(twigstats)

path <- paste0(system.file("sim/", package = "twigstats"),"/")
pref <- "msprime_ad0.8_split250_1"
file_plink <- paste0(path,pref) #only need prefix
file_mut  <- paste0(path,pref) #only need prefix (here same name as plink file but can be different)

#Compute regular f2 statistics between all pairs of populations. You can use pops to only calculate f2s between specified populations.
f2_blocks1 <- f2_blocks_from_RelateAges(pref = file_plink, file_mut)
f4_ratio(f2_blocks1, popX="PX", popI="P1", pop1="P2", pop2="P3", popO="P4")

#Use a cutoff of 500 generations
f2_blocks2 <- f2_blocks_from_RelateAges(pref = file_plink, file_mut, t = 500)
f4_ratio(f2_blocks2, popX="PX", popI="P1", pop1="P2", pop2="P3", popO="P4")
```
- Use argument <b>pops</b> to analyse a subset of populations and <b>fam</b> to change population assignment of individuals
- Use argument <b>blgsize</b> to change jackknife block size
- Set argument <b>transitions</b> to False to exclude transitions
- Use file_mut = "./1240k/1240k" to calculate f2 statistics at 1240k sites
- Use file_mut = "./1000G_ages/1000GP_Phase3_mask_prene_nosingle_chr1" to calculate f2 statistics at 1240k sites

#### Expected output
In the above example, PX is admixed between P2 and P3. The f4-ratio statistic computes an estimate of this admixture proportion.<br>
The expected output from running the code above is
```
#without twigstats ascertainment
Number of SNPs used: 396386
  popO popI pop1 pop2 popX     val        se
1   P4   P1   P2   P3   PX 1.56105 0.3807564

#with twigtstats ascertainment of 500 generations
Number of SNPs used: 78888
  popO popI pop1 pop2 popX       val         se
1   P4   P1   P2   P3   PX 0.7817311 0.02798447
```
So while we are unable to get a reliable estimate without twigstats ascertainment, we estimate PX to carry approximately 80% P2 and 20% P3 ancestry with twigstats.

#### Using admixtools2
You can now directly feed this into admixtools2 functions:
```R
library(admixtools)
#no twigstats ascertainment
f4(f2_blocks1, "P4", "P3", "P2", "PX")
#with twigstats ascertainment
f4(f2_blocks2, "P4", "P3", "P2", "PX")
```
The expected outcome of running this code is
```
# A tibble: 1 × 8
  pop1  pop2  pop3  pop4        est       se      z     p
  <chr> <chr> <chr> <chr>     <dbl>    <dbl>  <dbl> <dbl>
1 P4    P3    P2    PX    -0.000698 0.000723 -0.965 0.334
# A tibble: 1 × 8
  pop1  pop2  pop3  pop4       est        se     z        p
  <chr> <chr> <chr> <chr>    <dbl>     <dbl> <dbl>    <dbl>
1 P4    P3    P2    PX    0.000600 0.0000318  18.9 2.29e-79
```
