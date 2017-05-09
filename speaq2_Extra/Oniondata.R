#' Onion Intake in Mice Dataset
#'
#' 1H-NMR data of 32 rats fed different onion diets to identify onion intake biomarkers.
#'
#' @docType data
#'
#' @usage data(Oniondata)
#'
#' @format A list with the spectra, ppm vector and onion intake percentages. The peak data (after running the speaq 2.0 function 'getWaveletPeaks') is also included since it is needed in the vignette but would otherwise require too much processing power to build the vignettes (which are run daily for quality checks).
#'
#' @keywords datasets
#'
#' @references  Winning et al. (2009)  An exploratory NMR nutri-metabonomic investigation reveals dimethyl sulfone as a dietary biomarker for onion intake. Analyst 134 (2009) 2344-2351
#' (\href{http://pubs.rsc.org/is/content/articlelanding/2009/an/b918259d/unauth#!divAbstract}{DOI: 10.1039/B918259D})
#'
#' @source \href{http://www.models.life.ku.dk/onionnmr}{University of Copenhagen, Dept. of Food Science, Quality & Technology}
#'
#' @examples
#' data(Oniondata)
#' Spectra <- Oniondata$spectra 
#' ppm.wine <- Oniondata$ppm 
#' onion.percent <- Oniondata$onion.percent 
#' onion.grouplabels <- Oniondata$grouplabels
"Oniondata"
