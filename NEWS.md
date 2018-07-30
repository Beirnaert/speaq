# speaq 2.3.4
### 30/07/2018

* fixed an issue with DrawSpecPPM and detectSpecPeaks. If the ppm vector or the spectra matrix contained NA's, the function would throw an uninformative error. Now the software tries to fix it, and throws an informative error if it can't. Thanks to Precious Kwadzo Pomary for the alert.
* updated SilhouetR documentation.

# speaq 2.3.3
### 03/07/2018

* fixed an issue with the detectSpecPeaks. If no sample labels are provided it would give a confusing warning. 


# speaq 2.3.2
### 19/04/2018

* fixed an issue with the detectSpecPeaks. An error would occur if no peaks are found in a spectrum.
* fixed an issue with getWaveletPeaks. In case very short spectra where analysed an error would occur because of an optimization step that required spectra to be at least 512 or 1024 measurement points long (for the FFT). This is fixed now. Thanks to Pedro Lopez Garcia for the alert.


# speaq 2.3.1
### 27/03/2018

* fixed an issue with the SilhouetR which did not return the correct index values of the groups.
* peakGrouper now allows the use of non numeric sample labels (character vectors, factors). It would produce errors otherwise. Thanks to Guillaume Marti for the alert.
* getWaveletPeaks now contains a multiplier for the duplicate detection limit. This is usefull when working with non NMR spectra that are possibly a bit more distorted than NMR spectra. It allows a wider window to indicate peaks as multiple detections.


# speaq 2.3:
### 16/03/2018

* fixed issue of condition length 1. Thanks to Tomas Kalibera for the alert.
* fixed PeakFilling issue if MassSpecWaveet returns numeric(0).


# speaq 2.2:
### 01/03/2018

* changed the relevant.feature.p function to include multiple responses
* typo fixes
* vignette updates (now with performance analysis for a simulated dataset)
* peak picker supports raw peak height from now on
* fixed a "lostpeak" bug in dohClusterCustommedSegments(). Thanks to Manolis Matzapetakis for the alert.


# speaq 2.1:
### 14/09/2017

* HMDB searchfunction added
* general bugfixes


# speaq 2.0:
### 12/05/2017

### implemented new functionality to allow peak based analysis of NMR spectra by using wavelet based peak detection. New functions include:

* Peak picking with wavelets
* Peak grouping
* Peak Filling
* Linear model based differential analysis
* Silhouette values (check for alignment quality)
* Converting raw spectra of unequal length (unequal measurement time) to matrix of equal length spectra
* SCANT: a function to scale, normalise or transform a data matrix (included besides the standard are pareto scaling, probabilistic quotient normalization, range scaling, etc.)
* Plot functions for raw spectra
* etc.


# speaq 1.2.3:
### 25/02/2017

* Allow to automatically detect the optimal value for maxShift in function dohCluster(). The default setting (maxShift=100) usually works well for NMR spectra. However, for other types of spectra such as chromatograms, this value might be too large. In this new version, when the value of maxShift is set by NULL (maxShift=NULL), CluPA will automatically learn to select the optimal value based on the median Pearson correlation coefficent between spectra. It is worth noting that this metric is significantly effected by high peaks in the spectra, so it might not be the best measure for evaluating alignment performances. However, it is fast for the purpose of detecting the suitable *maxShift* value. A plot of correlation across the maxShift values also reported if the verbose=TRUE is applied.
* Do scaling data before Fast Fourier Transformation (FFT) cross-correlation in function findShiftStepFFT() if the input spectra are very low abundant (possible in chromatograms).
* Fix small bugs of detectSpecPeaks() when errors happen in function peakDetectionCWT() of MassSpecWavelet package.


# speaq 1.2.2:
### 15/11/2016

* Fix the issue of "if (condition) return;" might happen in function dohClusterCustommedSegments(). I acknowledge Duncan Murdoch <murdoch.duncan@gmail.com> for the alert.


# speaq 1.2.1:
### 24/02/2015

* Replace R version depends to R (>= 3.1.0) in order to remove the error of using anyNA().
* Remove the period mark in the end of the package title.
* Convert the title field to title case.
