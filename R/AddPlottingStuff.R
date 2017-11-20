#' Add plotting variables
#'
#' This functions adds a few variables which make plotting features easier (and more informative). Since for example every peaks keeps it original ppm value, if you want to plot the groups this function adds the group ppm value. Also sample labels can be added.
#'
#' @param Y.peaks data frame obtained from either of the 'getWaveletPeaks', 'PeakAligner' or 'PeakFilling' function.
#' @param X.ppm The vector with the ppm values (numeric vector).
#' @param groupLabels The groupLabels (numeric or factor).
#' 
#' @return Returns a data frame with added plotting variables (groupPPM for aligned features and labels for plotting).
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' subset <- GetWinedata.subset()
#' subset.spectra = as.matrix(subset$Spectra)
#' subset.ppm = as.numeric(subset$PPM)
#' 
#' test.peaks <- getWaveletPeaks(Y.spec=subset.spectra, 
#'                               X.ppm=subset.ppm,
#'                               nCPU = 2) # nCPU set to 2 for the vignette build
#'                               
#' test.peaks.plot = AddPlottingStuff(test.peaks, subset.ppm, subset$Color)
#' head(test.peaks.plot)
#' 
#' @export
AddPlottingStuff <- function(Y.peaks, X.ppm = NULL, groupLabels = NULL) {
    
    if (!is.null(X.ppm)) {
        if ("matrix" %in% class(X.ppm)) {
            stop("X.ppm is a matrix, for the plotting stuff only a numeric vector is allowed.")
        } else if (!"numeric" %in% class(X.ppm)) {
            X.ppm <- as.numeric(as.character(X.ppm))
        }
        Y.peaks$groupPPM <- as.numeric(X.ppm[Y.peaks$peakIndex])
    }
    
    if (!is.null(groupLabels)) {
        
        if (! "factor" %in% class(groupLabels)) {
            groupLabels <- as.factor(groupLabels)
        }
        if (length(groupLabels) != length(unique(Y.peaks$Sample))) {
            warning("There are not an equal amount of groupLabels as there are samples in Y.peaks")
        }
        if (length(groupLabels) < length(unique(Y.peaks$Sample))) {
            stop("There are less groupLabels than there are samples in Y.peaks")
        }
        Y.peaks$label <- groupLabels[Y.peaks$Sample]
    }
    
    return(Y.peaks)
}
