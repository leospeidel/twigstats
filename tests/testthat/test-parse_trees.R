
test_that("Run f2 function", {

	filename_anc <- system.file("sim", "msprime_ad0.8_split250_1_chr1.anc.gz", package = "twigstats", mustWork = T)
	filename_mut <- system.file("sim", "msprime_ad0.8_split250_1_chr1.mut.gz", package = "twigstats", mustWork = T)
	filename_popl <- system.file("sim", "msprime_ad0.8_split250_1.poplabels", package = "twigstats", mustWork = T)
	file_map <- system.file("sim", "genetic_map_combined_b37_chr1.txt.gz", package = "twigstats", mustWork = T)

	expect_error(f2_blocks_from_Relate(filename_anc, filename_mut, filename_popl, file_map), NA);

})
