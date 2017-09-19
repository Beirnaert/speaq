#' Get subset of Winedata for code examples
#'
#' This functions extracts a small part of the Winedata to be used in code testing and code examples
#' 
#' @return list of 3: spectra, ppm values, color and origin.
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
    Spectra <- speaq::Winedata$spectra[1:3, 1001:2000]
    PPM.vector <- speaq::Winedata$ppm[1001:2000]
    wine.color <- as.character(speaq::Winedata$wine.color[1:3])
    wine.origin <- as.character(speaq::Winedata$origin[1:3])
    return(list(Spectra = Spectra, PPM = PPM.vector, Color = wine.color, Origin = wine.origin))
}