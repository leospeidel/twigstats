/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include <testthat.h>
#include <RcppArmadillo.h>

#include "./data.hpp"
#include "./mutations.hpp"

std::string get_datapath_fromR(std::string filename, std::string packageName){
    Rcpp::Environment base("package:base");
    Rcpp::Function sys_file = base["system.file"];
    // "inst" field is necessary at this point
    Rcpp::StringVector file_path_sv = sys_file(
        "inst", "sim", filename,
        Rcpp::_["package"] = packageName,
        Rcpp::_["mustWork"] = true
    );
    std::string file_path = Rcpp::as<std::string>(file_path_sv);
    return file_path;
}

//Test parsing anc/mut
context("Parsing Relate files") {
	
	test_that("Open files and count number of tips") {

		AncMutIterators ancmut(get_datapath_fromR("msprime_ad0.8_split250_1_chr1.anc.gz", "twigstats"), get_datapath_fromR("msprime_ad0.8_split250_1_chr1.mut.gz", "twigstats"));
		Data data(ancmut.NumTips(), ancmut.NumSnps());
		expect_true( data.N == 100 );

	}

}

