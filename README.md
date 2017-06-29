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

### older versions

speaq versions <= 1.2.3 are still available at https://github.com/nghiavtr/speaq


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


## A minimal example

```R
library(speaq)

# get the data. (spectra in matrix format)
spectra.matrix = as.matrix(Winedata$spectra)
ppm.vector = as.numeric(Winedata$ppm)
class.factor = as.factor(Winedata$wine.color)

# plot the spectra
drawSpecPPM(Y.spec = spectra.matrix, X.ppm = ppm.vector, groupFactor = class.vector, title = 'Example spectra')

# peak detection
peaks <- getWaveletPeaks(Y.spec = spectra.matrix, X.ppm = ppm.vector)  # the default mode is parallel

# peak grouping
groups <- PeakGrouper(Y.peaks = test.peaks)

# get the feature matrix
Features <- BuildFeatureMatrix(Y.data = groups)

# this featurematrix can be processed further (scaling, transformations) or analysed with the statistical tools of interest like PCA. 
```

### <a name="Refs"></a> References 

Vu TN, Valkenborg D, Smets K, Verwaest KA, Dommisse R, Lemière F, Verschoren A, Goethals B, Laukens K. An integrated workflow for robust alignment and simplified quantitative analysis of NMR spectrometry data. BMC Bioinformatics, 2011, 12:405.

**preprint** Beirnaert C, Meysman P, Vu TN, Hermans N, Apers S, Pieters L, Covaci A, Laukens K (2017) speaq 2.0: a complete workflow for high-throughput 1D NMR spectra processing and quantification. bioRxiv preprint server, 2017, doi: https://doi.org/10.1101/138503

### License

The **speaq** package is licensed under the Apache License v2.0 


