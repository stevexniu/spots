[![Build status](https://ci.appveyor.com/api/projects/status/lbdrt72vyiqwxxff/branch/main?svg=true)](https://ci.appveyor.com/project/stevexniu/spots/branch/main) [![codecov](https://codecov.io/gh/stevexniu/spots/branch/main/graph/badge.svg?token=7KF4D3GGUB)](https://codecov.io/gh/stevexniu/spots)

# spots: Spatial Component Analysis <img src="man/figures/logo.png" align="right" width="150"/>

The **```spots```** package is designed for spatial omics (10x Visium, etc.) data analysis. 

It performs various statistical analyses and tests, including spatial component analysis (SCA), both global and local spatial statistics, such as univariate and bivariate Moran's I, Getis-Ord Gi* statistics, etc.

See <a href="https://doi.org/10.1101/2022.03.15.484516" target="_blank">Integrated protein and transcriptome high-throughput spatial profiling (2022)</a> for more details.

The ```feather``` branch is a light-weight version without ```SCA``` function and subsequent eigen decomposition steps.

Installation
-----
Install from Github:

``` r
install.packages("devtools")
devtools::install_github("stevexniu/spots", ref = "feather")
```

Usage
-----
See tutorials:

[Get Started](https://stevexniu.github.io/spots/articles/get_started.html)

[Spatial Statistics](https://stevexniu.github.io/spots/articles/spatial_statistics.html)

© [X. Steve Niu](https://github.com/stevexniu) from [Landau Lab](https://www.landaulab.org), Weill Cornell Medicine and New York Genome Center
