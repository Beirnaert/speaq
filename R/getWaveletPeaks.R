#' Convert raw NMR spectra to peak data by using wavelets
#'
#' This function converts phase corrected NMR spectra to peak data by using wavelet based peak detection (with the MassSpecWavelet package)
#'
#' @param X.ppm The x/ppm values of the spectra (in single vector or matrix format).
#' @param Y.spec The spectra in matrix format (rows = samples, columns = measurement points ).
#' @param sample.labels The sample labels (numeric), if not supplied these will simply be the sample numbers.
#' @param window.width The width of the detection window for the wavelets. Because of the Fourier transform lengths of 512 ( window.width = 'small') of 1024 ( window.width = 'large') are preferable.
#' @param window.split A positive, even and whole number indicating in how many parts the sliding window is split up. With every iteration the window slides one part further.
#' @param scales The scales to be used in the wavelet based peak detection, see \link[MassSpecWavelet]{peakDetectionCWT}.
#' @param baselineThresh Peaks with a peakValue lower than this threshold will be removed (default = 1000).
#' @param SNR.Th The Signal-to-noise threshold, see \link[MassSpecWavelet]{peakDetectionCWT}.
#' @param nCPU The amount of cpu's to be used for peak detection. If set to '-1' all available cores minus 1 will be used.
#' @param include_nearbyPeaks If set to TRUE small peaks in the tails of larger ones will be included in the peak data, see \link[MassSpecWavelet]{peakDetectionCWT}.
#' @param raw_peakheight (default = FALSE) Whether to use the raw peak height of a peak instead of the optimal CWT coefficient (which is a measure for AUC).
#'
#' @return The peaks detected with the wavelets.
#' 
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' subset <- GetWinedata.subset()
#' subset.spectra = as.matrix(subset$Spectra)
#' subset.ppm = as.numeric(subset$PPM)
#'
#' test.peaks <- getWaveletPeaks(Y.spec=subset.spectra, 
#'                               X.ppm=subset.ppm ,
#'                               nCPU = 1) # nCPU set to 2 for the vignette build
#'
#' @export
#' 
#' @importFrom foreach %dopar% foreach
#' @importFrom data.table rbindlist
#' @importFrom MassSpecWavelet peakDetectionCWT tuneInPeakInfo
#' @importFrom stats median
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
#' 
getWaveletPeaks <- function(Y.spec, X.ppm, sample.labels = NULL, window.width = "small", window.split = 4, 
                            scales = seq(1, 16, 1), baselineThresh = 1000, SNR.Th = -1, nCPU = -1, include_nearbyPeaks = TRUE,
                            raw_peakheight = FALSE) {
    # require(data.table) require(MassSpecWavelet) require(parallel) require(foreach)
    
    # error checks and miscellaneous parameter fixing
    for (w in 1:length(window.split)) {
        if (!window.split[w] %in% c(2, 4, 16, 32, 64)) {
            warning(paste(" 'window.split value", as.character(window.split[w]), "is not a power of 2 (max 64). It is set to the default: 4", 
                          sep = " "))
        }
    }
    
    if (nCPU == -1) {
        nCPU <- parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1
    }
    
    
    if ("small" %in% window.width & !"large" %in% window.width | missing(window.width)) {
        FFTwindow <- 512
    } else if ("large" %in% window.width & !"small" %in% window.width) {
        FFTwindow <- 1024
    } else {
        warning("'window.width' is defined ambiguously or wrong, set to default small.")
        FFTwindow <- 512
    }
    
    
    if (class(Y.spec) != "matrix") {
        print("the raw spectra, Y.spec, are not in matrix format, attempting conversion")
        warning("the raw spectra, Y.spec, are not in matrix format, conversion attempted")
        if ("numeric" %in% class(Y.spec)) {
            Y.spec <- matrix(Y.spec, nrow = 1, ncol = length(Y.spec))
        } else {
            Y.spec <- as.data.frame(Y.spec)
            Y.spec <- as.matrix(Y.spec)
        }
    }
    
    nPPM <- ncol(Y.spec)  # the number of measurement points are represented by the amount of columns in the matrix
    nSamp <- nrow(Y.spec)  # the number of samples are represented by the amount of rows in the matrix
    
    if(nCPU > nSamp){
        nCPU <- nSamp
    }
    
    if("data.frame" %in% class(X.ppm)){
        print("X.ppm is a data frame, attempting conversion." )
        warning("X.ppm was in data frame format. Conversion to numeric vector or matrix attempted")
        if(1 %in% dim(X.ppm)){
            X.ppm = as.numeric(X.ppm)
        } else if(nPPM %in% dim(X.ppm) & nSamp %in% dim(X.ppm)){
            X.ppm = as.matrix(X.ppm)
            if(ncol(X.ppm) == nSamp){
                X.ppm = t(X.ppm)
                warning("PPM matrix was transposed to have the samples in the matrix rows.")
            }
        }
    }
    
    if (is.numeric(X.ppm)) {
        if (length(X.ppm) == nPPM) {
            X.ppm.matrix <- matrix(rep(X.ppm, nSamp), ncol = nPPM, nrow = nSamp, byrow = TRUE)
        } else {
            stop("the length of the ppm vector, X.ppm, does not match with the amount of columns in the data matrix, Y.ppm.")
        }
    } else if (is.matrix(X.ppm)) {
        if (ncol(X.ppm) != nPPM | nrow(X.ppm) != nSamp) {
            stop("X.ppm is a matrix but the dimensions do not match with the Y.ppm matrix")
        }
        X.ppm.matrix <- X.ppm
    } else {
        stop("X.ppm is not a numeric nor a matrix. If X.ppm is a data frame: convert it to a matrix or a numeric vector")
    }
    
    
    # the actual function
    
    if (is.null(sample.labels)) {
        sample.labels <- seq(from = 1, to = nSamp)
    } else if (nSamp != length(sample.labels)) {
        warning("Sample labels do not match amount of rows in Y matrix, default row numbers will be used as sample labels")
        sample.labels <- seq(from = 1, to = nSamp)
    }
    
    noiseEsp <- 0.005
    if (SNR.Th < 0) 
        SNR.Th <- max(scales) * 0.05
    
    
    # library(doParallel)
    print("detecting peaks")
    cl <- parallel::makeCluster(nCPU)
    #doParallel::registerDoParallel(cl)
    doSNOW::registerDoSNOW(cl)
    peakList <- list()
    Parcounter <- NULL
    pb <- txtProgressBar(max=nSamp, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    peakList <- foreach::foreach(Parcounter = 1:nSamp, .options.snow=opts, .inorder = TRUE, .packages = c("data.table", "MassSpecWavelet")) %dopar% 
    {
        # noiselevel.estimate = NULL noisecounter = 0 library(MassSpecWavelet) library(data.table)
        teller <- 0
        WavPeaks <- list()
        currentSpec <- as.numeric(Y.spec[Parcounter, ])
        ppm_vector <- X.ppm.matrix[Parcounter, ]
        for (kw in 1:length(window.split)) {
            window.increment <- FFTwindow/window.split[kw]  # determine width of window increment
            nshifts <- ceiling(nPPM/window.increment) - window.split + 1  # determine how many increments there have to be
            for (j in 1:nshifts) {
                startR <- (j - 1) * window.increment + 1
                if (startR >= nPPM) 
                    next
                endR <- (j - 1) * window.increment + FFTwindow
                if (endR > nPPM) {
                    endR <- nPPM
                    startR <- nPPM - FFTwindow + 1
                }
                subSpec <- currentSpec[startR:endR]
                subMean <- mean(subSpec)
                subMedian <- stats::median(subSpec)
                subMax <- max(subSpec)
                subLeft <- 0
                # extra for when small spectra are analysed (zero padding left and right to make the spectrum of
                # length = FFTwindow )
                if (length(subSpec) < FFTwindow) {
                    subSpec.extended <- rep(0, FFTwindow)
                    subLeft <- floor((FFTwindow - (endR - startR + 1))/2)
                    subRight <- ceiling((FFTwindow - (endR - startR + 1))/2)
                    subSpec.extended[(subLeft + 1):(FFTwindow - subRight)] <- subSpec
                    subSpec <- subSpec.extended
                }
                if ((subMean == subMedian) || abs(subMean - subMedian)/((subMean + subMedian) * 2) < 
                    noiseEsp || subMax < baselineThresh) {
                    # noisecounter = noisecounter + 1 noiselevel.estimate[noisecounter] = subMean
                    next  # there is only noise
                } else {
                    Wavelet.Peaks <- -1
                    try({
                        peakInfo <- MassSpecWavelet::peakDetectionCWT(subSpec, scales = scales, SNR.Th = SNR.Th, 
                                                                      nearbyPeak = include_nearbyPeaks)
                        majorPeakInfo <- peakInfo$majorPeakInfo
                        
                        betterPeakInfo <- MassSpecWavelet::tuneInPeakInfo(subSpec, majorPeakInfo)
                        betterPeakInfo$sample <- sample.labels[Parcounter]
                        
                        peakIndex <- as.matrix(betterPeakInfo$peakIndex + startR - 1 - subLeft, ncol = 1)  # get the true index, not the one from the subprofile
                        peakPPM <- as.matrix(ppm_vector[peakIndex], ncol = 1)
                        peakSNR <- as.matrix(betterPeakInfo$peakSNR, ncol = 1)
                        peakValue <- as.matrix(betterPeakInfo$peakValue, ncol = 1)
                        peakScale <- as.matrix(betterPeakInfo$peakScale, ncol = 1)
                        Sample <- betterPeakInfo$sample
                        
                        Wavelet.Peaks <- data.frame(peakIndex, peakPPM, peakValue, peakSNR, peakScale, 
                                                    Sample)
                        
                    }, silent = TRUE)
                    
                    oldWarningLevel <- getOption("warn")
                    options(warn = -1)
                    
                    if (head(Wavelet.Peaks, n = 1)[1] != -1) {
                        teller <- teller + 1
                        WavPeaks[[teller]] <- Wavelet.Peaks
                    }
                    options(warn = oldWarningLevel)
                }
            }
        }
        
        WavPeaks <- data.table::rbindlist(WavPeaks)
        
        WavPeaks <- WavPeaks[!duplicated(WavPeaks$peakIndex), ]  # Remove all duplicated elements
        WavPeaks <- WavPeaks[WavPeaks$peakValue > baselineThresh, ]  # remove all elements with smaller than baseline threshold
        
        return(WavPeaks)
        
        
    }
    close(pb)
    parallel::stopCluster(cl)
    WaveletPeaks <- data.table::rbindlist(peakList)
    
    if(nrow(WaveletPeaks) != 0){
        
        ###### Fixing the duplicate detections
        print("fixing duplicate detections")
        
        window.increment <- FFTwindow/min(window.split)  # determine width of window increment (take the smalest window width)
        nshifts <- ceiling(nPPM/window.increment) - window.split + 1  # determine how many increments there where
        startR <- rep(NA,nshifts)
        endR <- rep(NA,nshifts)
        for (j in 1:nshifts) {
            startR[j] <- (j - 1) * window.increment + 1
            endR[j] <- (j - 1) * window.increment + FFTwindow
            if (startR[j] >= nPPM){
                startR[j] <- NA
                endR[j] <- NA
            }
        }
        startR <- startR[complete.cases(startR)]
        endR <- endR[complete.cases(startR)]
        
        
        cl <- parallel::makeCluster(nCPU)
        #doParallel::registerDoParallel(cl)
        doSNOW::registerDoSNOW(cl)
        to.delete <- list()
        Parcounter <- NULL
        pb <- txtProgressBar(max=(nshifts-1), style=3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress=progress)
        
        to.delete <- foreach::foreach(Parcounter = 1:(nshifts-1), .options.snow=opts, .inorder = TRUE, .packages = c("cluster")) %dopar% 
        {
            check.dat <- WaveletPeaks[WaveletPeaks$peakIndex>=(startR[Parcounter]) & WaveletPeaks$peakIndex<=endR[Parcounter+1],] # the data in two neighbouring windows
            if(nrow(check.dat) != 0){
                check.dist <- cluster::daisy(matrix(check.dat$peakPPM,ncol=1), metric = "euclid", stand = FALSE) # cluster matrix
                check.distM <- as.matrix(check.dist)
                check.distM[lower.tri(check.distM,diag=T)] <- 100 # only take the top triangle (identical to down triangle) and remove the diag because this is obviously 0 and of no interest
                # take the first quartile as a maximal distance
                indices = which(check.distM <= 0.005,arr.ind = TRUE) # which indices have a distance smaller 
                
                
                distance.data <- matrix(NA,nrow = nrow(indices), ncol = 9)
                colnames(distance.data) <- c("distance", "sample1","sample2","peakIndex1", "peakIndex2", "ppm.diff", "peakValue.diff.abs", "peakValue1", "peakValue2")
                distance.data[,1] <- check.distM[indices]
                distance.data[,2] <- check.dat$Sample[indices[,1]]
                distance.data[,3] <- check.dat$Sample[indices[,2]]
                distance.data[,4] <- check.dat$peakIndex[indices[,1]]
                distance.data[,5] <- check.dat$peakIndex[indices[,2]]
                distance.data[,6] <- abs(check.dat$peakPPM[indices[,2]]-check.dat$peakPPM[indices[,1]])
                distance.data[,7] <- check.dat$peakValue[indices[,2]]-check.dat$peakValue[indices[,1]]
                distance.data[,8] <- check.dat$peakValue[indices[,1]]
                distance.data[,9] <- check.dat$peakValue[indices[,2]]
                
                distance.data[,7] <- distance.data[,7, drop = FALSE] /apply(distance.data[,8:9, drop = FALSE], MARGIN = 1, max)
                distance.data <- distance.data[,1:7, drop = FALSE]
                colnames(distance.data) <- c("distance", "sample1","sample2","peakIndex1", "peakIndex2", "ppm.diff", "peakValue.diff.prop" )
                
                distance.data <- distance.data[order(distance.data[,1]), , drop = FALSE]
                distance.data <- distance.data[distance.data[,2]==distance.data[,3], ,drop = FALSE]
                
                left.largest  <- which(sign(distance.data[,7])==-1)
                distance.data[left.largest, c(4,5)] <- distance.data[left.largest, c(5,4)]
                distance.data[left.largest, 7] <- abs(distance.data[left.largest, 7])
                to.delete <- distance.data[,c(3,4), drop = FALSE]
            } else{
                distance.data <- matrix(NA,nrow =0, ncol = 7)
                colnames(distance.data) <- c("distance", "sample1","sample2","peakIndex1", "peakIndex2", "ppm.diff", "peakValue.diff.prop" )
                to.delete <- distance.data[,c(3,4),drop = FALSE]
            }
            
            return(to.delete)
        }
        close(pb)
        parallel::stopCluster(cl)
        
        
        to.delete <- do.call(rbind, to.delete)
        to.delete <- unique(to.delete)
        
        for(sn in 1:nSamp){
            if(nrow(to.delete[,1,drop=FALSE]) !=0){
                curr.samp.to.delete <- to.delete[to.delete[,1,drop=FALSE]==sn,2]
                WaveletPeaks$peakIndex[WaveletPeaks$Sample == sn & WaveletPeaks$peakIndex %in% curr.samp.to.delete] <- NA
            }
        }
        
        
        WaveletPeaks <- WaveletPeaks[!is.na(WaveletPeaks$peakIndex), ,drop=FALSE]
        WaveletPeaks <- data.frame(WaveletPeaks)
        
        if(raw_peakheight){
            WaveletPeaks$peakValue = Y.spec[as.matrix(WaveletPeaks[,c("Sample", "peakIndex")])]
        }
        
    } else{
        print("No peaks detected")
    } 
    
    return(WaveletPeaks)
    
}

