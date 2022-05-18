[![CRAN status](https://www.r-pkg.org/badges/version/spots)](https://CRAN.R-project.org/package=spots) [![Build status](https://ci.appveyor.com/api/projects/status/lbdrt72vyiqwxxff/branch/main?svg=true)](https://ci.appveyor.com/project/stevexniu/spots/branch/main) [![Build Status](https://app.travis-ci.com/stevexniu/spots.svg?branch=main)](https://app.travis-ci.com/stevexniu/spots) [![codecov](https://codecov.io/gh/stevexniu/spots/branch/main/graph/badge.svg?token=7KF4D3GGUB)](https://codecov.io/gh/stevexniu/spots)

# spots: Spatial Component Analysis <img src="man/figures/logo.png" align="right" width="150"/>

The **```spots```** package is designed for spatial omics (10x Visium, etc.) data analysis. 

It performs various statistical analyses and tests, including spatial component analysis (SCA), both global and local spatial statistics, such as univariate and bivariate Moran's I, Getis-Ord Gi* statistics, etc.

See <a href="https://doi.org/10.1101/2022.03.15.484516" target="_blank">Integrated protein and transcriptome high-throughput spatial profiling (2022)</a> for more details.

Installation
-----
Install from CRAN release:

``` r
install.packages("spots")
```

Install from Github:

``` r
install.packages("devtools")
devtools::install_github("stevexniu/spots")
```

Install the light-weight ```spots-feather``` without eign decomposition (```SCA``` functionality).

``` r
install.packages("devtools")
devtools::install_github("stevexniu/spots", ref = "feather")
```

Usage
-----
See tutorials:

[SCA cs PCA](https://stevexniu.github.io/spots/articles/SCA_vs_PCA.html)

[Get Started](https://stevexniu.github.io/spots/articles/get_started.html)

[Spatial Statistics](https://stevexniu.github.io/spots/articles/spatial_statistics.html)

© [X. Steve Niu](https://github.com/stevexniu) from [Landau Lab](https://www.landaulab.org), Weill Cornell Medicine and New York Genome Center
