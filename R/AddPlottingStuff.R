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
#' \dontrun{
#' # This function works on a data frame resulting from the 'getWaveletPeaks' function 
#' # DetectedPeaks <- getWaveletPeaks(X.ppm= PPM.vector, Y=Y.spec,  baselineThresh = 10,nCPU  = 4)
#' Aligned.peaks <- PeakAligner(Y.peaks = DetectedPeaks)
#' }
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
        
        if ("factor" %in% class(groupLabels)) {
            groupLabels <- as.character(groupLabels)
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
