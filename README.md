
# nmabat

<!-- badges: start -->
<!-- badges: end -->

The `nmabat` package implements numerical methods for calculating and visualizing the bias adjustment thresholds of selected data points within a network meta-analysis (NMA). Through iterative modification of the data and reassessment of the decision function, the boundary at which bias changes the decision can be identified.  

## Installation

You can install the development version of `nmabat` from [GitHub](https://github.com/miriam-vb/nmabat) with:

``` r
# install.packages("remotes")
remotes::install_github("miriam-vb/nmabat", build_vignettes=TRUE)
```

The vignettes associated with the package contain worked examples illustrating applications of the method, and they may take a few minutes to build during installation. The vignettes can then be accessed with the following call:

``` r
browseVignettes("nmabat")
```
