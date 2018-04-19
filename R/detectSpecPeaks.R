#' Peak detection for spectra
#'
#' Divide the whole spectra into smaller segments and detect peaks by using MassSpecWavelet package. 
#' Note that, the peak lists could be found by using other methods, this function is just a choice.
#'
#'
#' @param X The spectral dataset in matrix format in which each row contains a single sample
#' @param nDivRange The size of a single small segment after division of spectra
#' @param scales The parameter of peakDetectionCWT function of MassSpecWavelet package, look it up in the original function.
#' @param baselineThresh It will remove all peaks under an intensity set by baselineThresh.
#' @param SNR.Th The parameter of peakDetectionCWT function of MassSpecWavelet package, look it up in the original function. If you set -1, the function will itself re-compute this value.
#' @param verbose A boolean value to allow print out process information.
#'
#' @return The peak lists of the spectra
#'
#'
#' @author Trung Nghia Vu
#'
#' @examples
#' res=makeSimulatedData();
#' X=res$data;
#' groupLabel=res$label;
#' peakList <- detectSpecPeaks(X,
#'                             nDivRange = c(128),
#'                             scales = seq(1, 16, 2),
#'                             baselineThresh = 50000,
#'                             SNR.Th = -1,
#'                             verbose=FALSE
#' );
#' 
#' @export
#' 
detectSpecPeaks <- function (X, nDivRange = 128, scales = seq(1, 16, 2), baselineThresh = 50000, 
SNR.Th = -1, verbose = TRUE) 
{
    nFea <- ncol(X)
    nSamp <- nrow(X)
    noiseEsp <- 0.005
    if (SNR.Th < 0){
        SNR.Th <- max(scales) * 0.05
    }
    pList <- NULL
    for (i in 1:nSamp) {
        myPeakRes <- NULL
        mySpec <- X[i, ]
        for (k in 1:length(nDivRange)) {
            divR <- nDivRange[k]
            for (j in 1:(trunc(nFea/divR) - 3)) {
                startR <- (j - 1) * divR + 1
                if (startR >= nFea){
                    startR <- nFea
                }
                endR <- (j + 3) * divR
                if (endR > nFea){
                    endR <- nFea
                }
                xRange <- mySpec[startR:endR]
                xMean <- mean(xRange)
                xMedian <- stats::median(xRange)
                peakDetectionFlag <- TRUE
                if ((xMean == xMedian) || abs(xMean - xMedian)/((xMean + xMedian) * 2) < noiseEsp) peakDetectionFlag = FALSE
                if (peakDetectionFlag) {
                    majorPeakInfo <- NULL
                    tryCatch({
                        peakInfo = MassSpecWavelet::peakDetectionCWT(mySpec[startR:endR], 
                                                                     scales = scales, SNR.Th = SNR.Th)
                        majorPeakInfo = peakInfo$majorPeakInfo
                    },  error=function(e){if (verbose) cat("\n", "Wannings for spectrum",i, "from function peakDetectionCWT() of MassSpecWavelet package: ",conditionMessage(e), "\n")})
                    
                    if(!is.null(majorPeakInfo)){
                        if (length(majorPeakInfo$peakIndex) > 0) {
                            print(paste("j =",j," peak =" , majorPeakInfo$peakIndex," savedpeak = ", (majorPeakInfo$peakIndex + startR - 1), sep = "" ))
                            myPeakRes <- c(myPeakRes, majorPeakInfo$peakIndex + startR - 1)
                        }
                    }
                }
            }
        }
        
        pList[i] <- list(myPeakRes)
        if(!is.null(myPeakRes)){
            pList[[i]] <- sort(unique(pList[[i]]))
            pList[[i]] <- pList[[i]][which(X[i, pList[[i]]] > baselineThresh)]
            pList[[i]] <- sort(pList[[i]])
        }
        if (verbose){
            cat("\n Spectrum", i, "has", length(pList[[i]]), "peaks")
        }
        
    }
    return(pList)
}
