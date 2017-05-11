# speaq: an R package for 1D NMR spectra data processing 

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.com/Beirnaert/speaq.svg?token=RasChTxxFxz6YzsLEhBK&branch=master)](https://travis-ci.com/Beirnaert/speaq)


The **speaq** package is meant to make Nuclear Magnetic Resonance spectroscopy (NMR spectroscopy) data analysis as easy as possible by only requiring a small set of functions to perform an entire analysis. Parallel processing is extensively implemented throughout these function to limit processing time.  


### Installation

``` r
install.packages("speaq")
```

### Features

Functions include

* Plot functions for raw spectra
* Peak picking with wavelets
* Peak grouping
* Peak Filling
* Linear model based differential analysis
* Silhouette values (check for alignment quality)
* Converting raw spectra of unequal length (unequal measurement time) to matrix of equal length spectra
* SCANT: a function to scale, normalise or transform a data matrix (included besides the standard are pareto scaling, probabilistic quotient normalization, range scaling, etc.)
* etc.



### Usage

See the vignette on CRAN for a usage case of most functions. See the GitHub repo for more extensive vignettes.



### License

The **speaq** package is licensed under the Apache License v2.0 


![Screenshot](docimages/speaq-github-workflow.png)
