# twigstats

Boost f-statitics power using genealogies. Compatible with [admixtools2](https://uqrmaie1.github.io/admixtools/index.html).<br/>
Full documentation at [leospeidel.com/twigstats](https://leospeidel.com/twigstats).<br/>
Each function is described on this page: [leospeidel.com/twigstats/reference](https://leospeidel.com/twigstats/reference).

## Installation

Please install this package by running the following command in R:
```R
library(remotes)
install_github("leospeidel/twigstats")
```

Alternatively, clone this directory (https://github.com/leospeidel/twigstats) and then in R type
```R
library(devtools)
install()
```

Please make sure to have an up-to-date R version (>=3.6.0) and C/C++ compiler (e.g., >=GCC v8) loaded in your environment.
If you encounter issues, it can help to create a clean R library. Assuming you place this in your home directory, this is done using
```
mkdir ~/R_libs_for_twigstats/
export R_LIBS_USER="~/R_libs_for_twigstats/"
```
Now, all required R packages will be installed into this fresh R library.

## Basic Usage

Please see [leospeidel.com/twigstats/articles/basic-usage.html](https://leospeidel.com/twigstats/articles/basic-usage.html).<br/>
For a small real data example see [leospeidel.com/twigstats/articles/real-data-example.html](https://leospeidel.com/twigstats/articles/real-data-example.html).
