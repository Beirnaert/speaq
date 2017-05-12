# speaq: an R package for 1D NMR spectra data processing 

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.com/Beirnaert/speaq.svg?token=RasChTxxFxz6YzsLEhBK&branch=master)](https://travis-ci.com/Beirnaert/speaq)


The **speaq** package is meant to make Nuclear Magnetic Resonance spectroscopy (NMR spectroscopy) data analysis as easy as possible by only requiring a small set of functions to perform an entire analysis. **speaq** offers two main possibilities for NMR spectral processing:

* Raw spectra alignment & processing (speaq version <= 1.2.3)

* Peak based processing (new since version 2.0)


A visual representation of what each method encompasses, and how both methods can be used together can be found in the image below. More specific information on both methods can be found in their respective publications (see [References](#Refs)). 

![Workflow](docimages/speaq-github-workflow.png)


## Installation


### From CRAN using R
```R
install.packages("speaq")   
library("speaq")
```

### From Github using R (with devtools):
```R
library("devtools")
install_github("speaq","beirnaert")
library("speaq")
```


## Features


Functions include

* Raw spectra 
	- spectra alignment (CluPA)
	- quantitation (BW-ratio)
	- visualizations
	- etc.

* Peak based approach (new since v2.0)
	- peak picking (wavelets)
	- grouping
	- peak filling
	- scaling & imputations
	- differential analysis (also more than 2 groups)
	- visualizations
	- etc.



### <a name="Refs"></a> References 

Vu TN, Valkenborg D, Smets K, Verwaest KA, Dommisse R, Lemière F, Verschoren A, Goethals B, Laukens K. (2011) An integrated workflow for robust alignment and simplified quantitative analysis of NMR spectrometry data. BMC Bioinformatics. 2011 Oct 20;12:405.

### License

The **speaq** package is licensed under the Apache License v2.0 


