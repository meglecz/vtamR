
<!-- README.md is generated from README.Rmd. Please edit that file -->

# vtamR

<!-- badges: start -->
<!-- badges: end -->

The goal of vtamR is to …

## Installation

You can install the development version of vtamR from
[GitHub](https://github.com/) with:

``` r
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pak("meglecz/vtamR")
#> ✔ Updated metadata database: 5.74 MB in 15 files.
#> ℹ Updating metadata database✔ Updating metadata database ... done
#>  
#> → Will update 1 package.
#> → Will download 1 package with unknown size.
#> + vtamR 0.2.0 → 0.0.1.0 [bld][cmp][dl] (GitHub: 4c6c288)
#> ℹ Getting 1 pkg with unknown size
#> ✔ Got vtamR 0.0.1.0 (source) (45.40 MB)
#> ℹ Packaging vtamR 0.0.1.0
#> ✔ Packaged vtamR 0.0.1.0 (25.9s)
#> ℹ Building vtamR 0.0.1.0
#> ✔ Built vtamR 0.0.1.0 (3.6s)
#> ✔ Installed vtamR 0.0.1.0 (github::meglecz/vtamR@4c6c288) (71ms)
#> ✔ 1 pkg: upd 1, dld 1 (NA B) [47.4s]
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(vtamR)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
