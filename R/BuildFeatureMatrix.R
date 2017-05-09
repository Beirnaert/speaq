#' Build a Feature matrix from the with speaq 2.0 processed data
#'
#' This function converts the grouped peak data to a matrix.
#' The matrix has features (peaks groups) in the columns and the value of the peak for every sample in the rows. 
#'
#' @param Y.data The dataset after (at least) peak detection and grouping with speaq 2.0.
#' @param var The variable to be used in the Featurematrix. This can be any of 'peakIndex', 'peakPPM', 'peakValue' (default), 'peakSNR', 'peakScale', or 'Sample'.
#' @param impute What to impute when a certain peak is missing for a certain sample and feature combo. Options are 'zero' (or 'zeros'), any other statement will produce NA's.
#' @param delete.below.threshold Whether to ignore peaks for which the 'var' variable has a value below 'baselineThresh' (default = FALSE).
#' @param baselineThresh The threshold for the 'var' variable peaks have to surpass to be included in the feature matrix.
#' @param snrThres The threshold for the signal-to-noise ratio of a peak.
#' @param thresholds.pass This variable lets users deside whether a peak has to pass all the thresholds (both snrThres and baselineThresh), or just one. (If the peak does not need to surpass any thresholds set 'delete.below.threshold' to FALSE). 
#' 
#' @return a matrix, data.matrix, with samples for rows and features for columns. The values in the matrix are thoes of the 'var' variable.
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' \dontrun{
#' # 'DetectedPeaks_grouped' is the peak data after wavelet based peak detection and grouping
#' Featurematrix <- BuildFeatureMatrix(Y.data = DetectedPeaks_grouped, var = 'peakValue')
#' }
#' @export

BuildFeatureMatrix <- function(Y.data, var = "peakValue", impute = "zero", delete.below.threshold = FALSE, 
    baselineThresh = 500, snrThres = 3, thresholds.pass = "any-to-pass") {
    if (!var %in% names(Y.data)) {
        stop("the variable of interest defined by 'var' is not in the names of the data")
    } else {
        VOI <- which(names(Y.data) == var)  # Variable Of Interest
    }
    
    if(!is.data.frame(Y.data) | length(class(Y.data)) > 1){
        Y.data <- data.frame(Y.data)
    }
    
    Features <- unique(Y.data$peakIndex)
    nFeatures <- length(Features)
    
    # if groups where all values are below threshold have to be deleted, the following loop is executed
    # to remove these groups from Y.data
    to.delete <- NULL
    if (delete.below.threshold) {
        if (!thresholds.pass %in% c("any-to-pass", "any to pass", "any", "all-to-pass", "all to pass", 
            "all")) {
            thresholds.pass <- "any"
            warning("'thresholds.pass' was not set to an appropriate value, set to 'any-to-pass' meaning that if the filled peak is higher than one of the thresholds it is kept in the results")
        } else if (thresholds.pass %in% c("any-to-pass", "any to pass", "any")) {
            thresholds.pass <- "any"
        }
        if (thresholds.pass == "any") {
            for (gg in 1:nFeatures) {
                # maxVal = max(Y.data$peakValue[Y.data$peakIndex==Features[gg] & !is.na(Y.data$peakSNR)])
                maxVal <- max(Y.data[Y.data$peakIndex == Features[gg] & !is.na(Y.data$peakSNR), VOI])
                maxSNR <- max(Y.data$peakSNR[Y.data$peakIndex == Features[gg] & !is.na(Y.data$peakSNR)])
                if (maxVal < baselineThresh & maxSNR < snrThres) {
                  to.delete <- c(to.delete, gg)
                }
                
            }
        } else {
            for (gg in 1:nFeatures) {
                # maxVal = max(Y.data$peakValue[Y.data$peakIndex==Features[gg] & !is.na(Y.data$peakSNR)])
                maxVal <- max(Y.data[Y.data$peakIndex == Features[gg] & !is.na(Y.data$peakSNR), VOI])
                maxSNR <- max(Y.data$peakSNR[Y.data$peakIndex == Features[gg] & !is.na(Y.data$peakSNR)])
                if (maxVal < baselineThresh | maxSNR < snrThres) {
                  to.delete <- c(to.delete, gg)
                }
                
            }
        }
    }
    
    if (length(to.delete) > 0) {
        Features <- Features[-to.delete]
        nFeatures <- length(Features)
    }
    
    Samples <- unique(Y.data$Sample)
    nSamples <- length(Samples)
    Samples <- Samples[order(Samples)]
    nSampl.seq <- seq(1, nSamples)
    
    if (impute == "zero" | impute == "zeros") {
        data.matrix <- matrix(0, nrow = nSamples, ncol = nFeatures)
    } else {
        data.matrix <- matrix(NA, nrow = nSamples, ncol = nFeatures)
    }
    
    for (k in 1:nFeatures) {
        curr.peak <- Y.data[Y.data$peakIndex == Features[k], , drop = FALSE]
        # data.matrix[ nSampl.seq[ Samples %in% curr.peak$Sample] ,k] =
        # curr.peak$peakValue[order(curr.peak$Sample)]
        data.matrix[nSampl.seq[Samples %in% curr.peak$Sample], k] <- curr.peak[order(curr.peak$Sample), 
            VOI]
    }
    colnames(data.matrix) <- Features
    rownames(data.matrix) <- Samples
    
    if (impute == "mean") {
        for (k in 1:nFeatures) {
            data.matrix[is.na(data.matrix[, k]), k] <- mean(data.matrix[!is.na(data.matrix[, k]), k])
        }
    }
    
    return(data.matrix)
}
