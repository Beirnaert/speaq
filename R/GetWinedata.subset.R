#' Get subset of Winedata for code examples
#'
#' This functions extracts a small part of the Winedata to be used in code testing and code examples
#' 
#' @return list of 2: spectra and ppm values
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' 
#' subset <- GetWinedata.subset()
#' subset.spectra = subset$Spectra
#' subset.ppm = subset$PPM.vector
#' 
#' @export
#' 
GetWinedata.subset <- function(){
    Spectra <- Winedata$spectra[1:4, 1001:2000]
    PPM.vector <- Winedata$ppm[1001:2000]
    
    return(list(Spectra = Spectra, PPM = PPM.vector))
}