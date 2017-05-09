#' Peak filling of any missed peaks
#'
#' This functions detects which samples (after grouping) are missing from every peak group and reanalyses the raw data to verify whether this peak is actually non-existent for this sample 
#'
#' @param Y.grouped Peaks groups (output from the 'PeakGrouper' function).
#' @param Y.spec The raw NMR spectra in matrix format.
#' @param max.index.shift Maximal shift in index between a filled peak and the group it belongs to.
#' @param window.width The width of the detection window for the wavelets. Because of the Fourier transform lengths of 512 ( window.width = 'small') of 1024 ( window.width = 'large') are preferable.
#' @param nCPU The amount of cpu's to be used for peak detection. If set to '-1' all available cores minus 1 will be used.
#'
#' @return Returns a data frame with grouped peaks and possibly extra peaks obtained from the raw data (these peaks have SNR = NA).
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' \dontrun{
#' # This function works on a data frame resulting from the 'PeakGrouper' function 
#' # DetectedPeaks <- getWaveletPeaks(X.ppm= PPM.vector, Y=Y.spec,  baselineThresh = 10,nCPU  = 4)
#' # Grouped.peaks = PeakGrouper = function (Y.peaks = DetectedPeaks)
#' Filled.Peaks = SpecPeak.filling(Y.grouped = Grouped.peaks, Y.spec = Y)
#' }
#' 
#' @export
#' 
#' @importFrom MassSpecWavelet tuneInPeakInfo
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom data.table rbindlist
#' @importFrom stats median na.omit
#' @importFrom foreach %dopar% foreach
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW

#'
PeakFilling <- function(Y.grouped, Y.spec, max.index.shift = 10, window.width = "small", nCPU = -1) {
    
    
    groups <- unique(Y.grouped$peakIndex)
    
    Samples <- sort(unique(Y.grouped$Sample))
    nSamples <- length(Samples)
    print(paste("There are", as.character(nSamples), "samples in the dataset."))
    
    if ("small" %in% window.width & !"large" %in% window.width | missing(window.width)) {
        FFTwindow <- 512
    } else if ("large" %in% window.width & !"small" %in% window.width) {
        FFTwindow <- 1024
    } else {
        warning("'window.width' is defined ambiguously or wrong, set to default small.")
        FFTwindow <- 512
    }
    
    if (nCPU == -1) {
        nCPU <- parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1
    }
    if(nCPU > nSamples){
        nCPU <- nSamples
    }
    
    
    groups.to.fill <- rep(NA, length(groups))
    # determine which groups (peakIndex) do not contain a peak for every sample
    for (g in 1:length(groups)) {
        if (length(unique(Y.grouped$Sample[Y.grouped$peakIndex == groups[g]])) != nSamples) 
            groups.to.fill[g] <- groups[g]
    }
    groups.to.fill <- as.numeric(stats::na.omit(groups.to.fill))
    
    if (length(groups.to.fill) == 0) {
        print("No peak group filling performed. All samples are present in every group. If this is unexpected see whether setting the 'nSamples' parameter manually.")
        return(Y.grouped)
    } else {
        
        
        
        if(length(groups.to.fill) > nSamples){
            group.division = cut(seq(1,length(groups.to.fill)),breaks=nSamples,labels=FALSE)
        }else{
            group.division = rep(1,length(groups.to.fill))
        }
        
        nRuns = length(unique(group.division))
        
        cl <- parallel::makeCluster(nCPU)
        #doParallel::registerDoParallel(cl)
        doSNOW::registerDoSNOW(cl)
        FilledList <- list()
        Parcounter <- NULL
        pb <- txtProgressBar(max=nRuns, style=3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress=progress)
        
        FilledList <- foreach::foreach(Parcounter = 1:nRuns, .options.snow=opts, .inorder = FALSE, .packages = c("stats", "MassSpecWavelet")) %dopar% 
        {
            
            groups.to.fill.subset = groups.to.fill[group.division == Parcounter]
            Peaks.filled <- list()
            
            for (gg in 1:length(groups.to.fill.subset)) {
                # for (gg in 1:length(28)){
                grpdata <- Y.grouped[Y.grouped$peakIndex == groups.to.fill.subset[gg], ]
                
                Mn.ppm <- mean(grpdata$peakPPM)
                Med.index <- stats::median(grpdata$peakIndex)
                Med.scale <- stats::median(grpdata$peakScale)
                
                samples.to.fill <- Samples[!Samples %in% grpdata$Sample]
                filled <- matrix(ncol = 6, nrow = length(samples.to.fill))
                
                for (k in 1:length(samples.to.fill)) {
                    
                    res <- tryCatch({
                        res.in.try <- MassSpecWavelet::tuneInPeakInfo(Y.spec[samples.to.fill[k], (Med.index - 
                                                                                                      floor((FFTwindow - 1)/2)):(Med.index + ceiling((FFTwindow - 1)/2))], peakIndex = FFTwindow/2, 
                                                                      peakScale = round(Med.scale)) 
                        # This rounding of Med.scale is very important as Massspecwavelet gives warnings in 
                        # the function tuneInPeakInfo -> getLocalMaximumCWT -> localMaximum
                    }, warning = function(warn) {
                        res.in.try <- rep(NA, 6)  # this will result in FWHM=-1 and resolving_power = -1 for this run
                        return(res.in.try)
                    }, error = function(err) {
                        res.in.try <- rep(NA, 6)  # this will result in FWHM=-1 and resolving_power = -1 for this run
                        return(res.in.try)
                    })  # end trycatch
                    
                    
                    res[sapply(res, is.null)] <- NA
                    if (length(res) != 6) 
                        res <- rep(NA, 6)  # !important after sapply <- NA! : sometimes the $peakCenterIndex from tuneInPeakInfo returns 'numeric(0)', throw it away
                    res <- unlist(res)
                    filled[k, ] <- res
                    
                }
                colnames(filled) <- c("peakIndex", "peakValue", "peakCenterIndex", "peakSNR", "peakScale", 
                                      "unProcessedPeak")
                filled <- data.frame(filled)
                filled$Sample <- samples.to.fill
                filled$peakPPM <- Mn.ppm
                filled$grp.index <- Med.index
                
                # delete the new peaks who's center index is shifted too far from the original center of the group
                # (larger than +/- max.index.shift)
                filled <- filled[filled$peakCenterIndex <= (FFTwindow/2 + max.index.shift) & filled$peakCenterIndex >= 
                                     (FFTwindow/2 - max.index.shift), , drop = F]
                # delete the ones where tuneInPeakInfo returned 6 NA
                filled <- filled[rowSums(is.na(filled[, 1:6])) != 6, , drop = F]
                # delete the peaks where the wavelet search returned NaN for peakValue
                filled <- filled[!is.nan(filled$peakValue), ]
                
                Peaks.filled[[gg]] <- filled
                
                # update progress bar
                #utils::setTxtProgressBar(pb, gg)
            }
            Peaks.filled <- data.table::rbindlist(Peaks.filled)
            Peaks.filled <- data.frame(Peaks.filled$grp.index, Peaks.filled$peakPPM, Peaks.filled$peakValue, 
                                       Peaks.filled$peakSNR, Peaks.filled$peakScale, Peaks.filled$Sample)
            colnames(Peaks.filled) <- colnames(Y.grouped)
            return(Peaks.filled)   
        }
        close(pb)  # close progress bar
        parallel::stopCluster(cl)
        FilledList = data.table::rbindlist(FilledList)
        
        
        Y.filled <- rbind(Y.grouped, FilledList)
        Y.filled <- data.frame(Y.filled)
        return(Y.filled)
    }
}



