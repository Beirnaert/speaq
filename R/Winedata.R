#' Wine dataset
#'
#' 1H-NMR data of 40 wines, different origins and colors are included.
#'
#' @docType data
#'
#' @usage data(Winedata)
#'
#' @format A list with the spectra, ppm values, color and origin as list entries.
#'
#' @keywords datasets
#'
#' @references  Larsen et al. (2006) An exploratory chemometric study of 1H-NMR spectra of table wines. J.Chemom. 20 (2006) 198-208
#' (\href{http://onlinelibrary.wiley.com/doi/10.1002/cem.991/full}{Wiley Online Library})
#'
#' @source \href{http://www.models.life.ku.dk/datasets}{University of Copenhagen, Dept. of Food Science, Quality & Technology}
#'
#' @examples
#' data(Winedata)
#' Spectra <- Winedata$spectra 
#' ppm.wine <- Winedata$ppm 
#' wine.color <- Winedata$wine.color 
#' wine.origin <- Winedata$origin 
"Winedata"
