#' Realignment with hierarchical clustering 
#'
#' This functions aligns the data frame obtained after wavelet based peak detection with the 'getWaveletPeaks' function
#'
#' @param Y.peaks data frame obtained from the 'getWaveletPeaks' function.
#' @param grouping.window.width The width of the sliding window (in measurement points). Measurments are taken for when this sliding window is taken too small, but best set this too a value that a normal peak is comfortably in a window. Note if large shifts occur in your dataset (like in the wine dataset) it is best to set this parameter larger.
#' @param verbose If set to TRUE the window selection process is documented in real time (default = FALSE).
#' @param min.samp.grp The minimal amount of samples needed te form a group, see \link[speaq2]{hclust.grouping}.
#' @param max.dupli.prop The maximal duplication proportion allowed for a group to be considered a sigle group, see \link[speaq2]{hclust.grouping}.
#' @param maxClust The maximum number of clusters (depth of the tree), see \link[speaq2]{hclust.grouping}.
#' @param Jaccard.regroup.threshold If 2 neighbouring groups have a jaccard index smaller than this 'Jaccard.regroup.threshold' (indicating that they are quite complementary as they have little peaks samples in common), then they are merged and regrouped. This situation can occur if a group is accidentally cut in half by the window approach.
#' @param linkage The linkage to be used in the hierarchical clustering. See the 'method' argument in \link[stats]{hclust}.
#'
#' @return Returns a data frame with realigned peaks (aligned on groupIndex).
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' \dontrun{
#' # This function works on a data frame resulting from the 'getWaveletPeaks' function 
#' # DetectedPeaks <- getWaveletPeaks(X.ppm= PPM.vector, Y=Y.spec,  baselineThresh = 10,nCPU  = 4)
#' Aligned.peaks = PeakAligner = function (Y.peaks = DetectedPeaks)
#' }
#' 
#' @export
#' 
#' @importFrom data.table rbindlist
#' @importFrom utils txtProgressBar
#' @importFrom stats complete.cases
#' 
PeakAligner <- function(Y.peaks, grouping.window.width = 100, verbose = FALSE, min.samp.grp = 1, max.dupli.prop = 0.25, 
    maxClust = 10, Jaccard.regroup.threshold = 0.25, linkage = "average") {
    
    
    .Deprecated("PeakGrouper")
    
    
    # Y.peaks is the result from detectspecPeaks.charlie in list or 1 large dataframe format
    
    if ("list" %in% class(Y.peaks)) {
        Y.peaks <- data.table::rbindlist(Y.peaks)
    } else if ("data.frame" %in% class(Y.peaks)) {
        Y.peaks <- Y.peaks
    } else {
        stop("NMR.peaks format is not a list or data.frame")
    }
    
    
    
    ##### Peak aligner
    
    Samples <- unique(Y.peaks$Sample)
    
    
    user.peaks.threshold <- 5
    Y.aligned <- matrix(NA, nrow = nrow(Y.peaks), ncol = ncol(Y.peaks))
    alignindex.start <- 1
    alignindex.stop <- 0
    alignindex.counter <- 0
    
    maxFea <- max(Y.peaks$peakIndex)
    minFea <- min(Y.peaks$peakIndex)
    
    nSamp <- length(unique(Y.peaks$Sample))
    
    window.length <- grouping.window.width  # 100 by default
    window.variation.global <- 0.1
    window.breaks <- rep(NA, round(maxFea/window.length) * 2)
    # what is done here: you take a window size roughly close to window.length but it is crucial that a
    # gaussian peak distribution is not cut in half therefor vary the right edge a bit, until there is a
    # plateau where no extra peaks appear.  then take the outer egde to work with, but when removing
    # these of the remaining.indexes, only remove until the inner edge. As to make sure that the next
    # iteration, the left edge is also seccured of some blank space. minFea = 9985 - 100 maxFea = 9985 +
    # 410
    remaining.index <- seq(from = minFea, to = maxFea, by = 1)
    max.iterator <- ceiling(length(remaining.index)/window.length) * 2
    iterator <- 0
    window.variation <- window.variation.global
    
    # progress bar
    print("realigning peaks")
    maxPBindex <- length(remaining.index)
    pb <- utils::txtProgressBar(min = 0, max = maxPBindex, style = 3)
    # end progress bar
    
    while (length(remaining.index) > 0 & iterator < max.iterator) {
        
        iterator <- iterator + 1
        
        # check that there is more than twice the window.length remaining in remaining.index
        if (length(remaining.index) > 2 * window.length) {
            # loop management checks
            nsubpeaks.NotOK <- FALSE
            no.ending <- FALSE
            window.overexpanded <- FALSE
            window.overshrunk <- FALSE
            
            min.edge <- round(window.length * (1 - window.variation))
            max.edge <- round(window.length * (1 + window.variation))
            
            var.seq <- seq(from = min.edge, to = max.edge, by = 1)
            seq.peaks <- NULL
            for (l in 1:length(var.seq)) {
                min.index_b <- remaining.index[1]
                max.index_b <- remaining.index[var.seq[l]]
                
                seq.peaks[l] <- nrow(Y.peaks[Y.peaks$peakIndex >= min.index_b & Y.peaks$peakIndex <= 
                  max.index_b, , drop = F])
            }
            
            
            
            
            if (max(seq.peaks) < 5 | max(seq.peaks) < (nSamp * 0.3)) {
                remaining.index <- remaining.index[-(1:min.edge)]
                if (verbose == TRUE) 
                  print("too little samples")
                nsubpeaks.NotOK <- TRUE
                # break # hop out of this loop and therefor delete this outlier peak (it is not part of a group)
            } else {
                
                peaks.edge.seq <- data.frame(var.seq, seq.peaks)
                
                seq.table <- data.frame(table(seq.peaks))
                seq.table$seq.peaks <- as.numeric(as.character(seq.table$seq.peaks))
                # endpoint.2nd.half = seq.table$seq.peaks[ seq.table$Freq == max(seq.table$Freq[
                # round(nrow(seq.table)/2) : nrow(seq.table)] ) ]
                endpoint.max <- seq.table$seq.peaks[seq.table$Freq == max(seq.table$Freq)]  # get the highest freq of number of peaks in the range. This is indicative of a plateau where nothing changes.
                if (length(endpoint.max) > 1) {
                  endpoint.max <- endpoint.max[1]
                }
                
                endpoint.peaks <- peaks.edge.seq[peaks.edge.seq$seq.peaks == endpoint.max, ]
                
                
                min.edge <- endpoint.peaks$var.seq[1]
                max.edge <- endpoint.peaks$var.seq[nrow(endpoint.peaks)]
                
                if (min.edge < (max.edge - 2)) {
                  confirmed.end <- TRUE
                } else {
                  confirmed.end <- FALSE
                }
                
                
                
                if (is.null(confirmed.end) & (window.variation < window.variation.global * 2) & (window.variation > 
                  window.variation.global/2) & (max(seq.peaks) < 10 * nSamp)) {
                  # no window.variation.global piece found with all zeros
                  window.variation <- window.variation * 1.1
                  if (verbose == TRUE) 
                    print("no ending found, window enlarged")
                  no.ending <- TRUE
                } else if (is.null(confirmed.end) & (window.variation > window.variation.global/2) & (window.variation < 
                  window.variation.global * 2) & (max(seq.peaks) > 10 * nSamp)) {
                  # no window.variation.global piece found with all zeros
                  window.variation <- window.variation/1.1
                  if (verbose == TRUE) 
                    print("no ending found, a lot of peaks: window halved")
                  no.ending <- TRUE
                } else if (is.null(confirmed.end) & (window.variation > window.variation.global * 2)) {
                  # no window.variation.global piece found with all zeros
                  remaining.index <- remaining.index[-(1:min.edge)]
                  if (verbose == TRUE) 
                    print("window increased too much")
                  window.overexpanded <- TRUE
                } else if (is.null(confirmed.end) & (window.variation < window.variation.global/2)) {
                  # no window.variation.global piece found with all zeros
                  remaining.index <- remaining.index[-(1:min.edge)]
                  if (verbose == TRUE) 
                    print("window got too small")
                  window.overshrunk <- TRUE
                } else {
                  peaks.in.seq <- peaks.edge.seq$seq.peaks[peaks.edge.seq$var.seq == max.edge]
                  
                  if (peaks.in.seq < 3 | peaks.in.seq < (nSamp * 0.2) | peaks.in.seq < user.peaks.threshold) {
                    remaining.index <- remaining.index[-(1:min.edge)]
                    if (verbose == TRUE) 
                      print("too little samples")
                    nsubpeaks.NotOK <- TRUE
                    # break # hop out of this loop and therefor delete this outlier peak (it is not part of a group)
                  }
                }
            }
        } else {
            
            nsubpeaks.NotOK <- FALSE
            no.ending <- FALSE
            window.overexpanded <- FALSE
            window.overshrunk <- FALSE
            min.edge <- length(remaining.index)
            max.edge <- min.edge
        }
        
        if (nsubpeaks.NotOK | no.ending | window.overexpanded | window.overshrunk) {
            # next
        } else {
            startR <- remaining.index[1]
            endR <- remaining.index[max.edge]
            
            window.breaks[iterator] <- endR
            
            remaining.index <- remaining.index[-(1:min.edge)]  # remove the chosen series (but not the buffer) from the remaining indexes
            window.variation <- window.variation.global  # reset window variation
            
            
            
            current.peaks <- Y.peaks[Y.peaks$peakIndex >= startR & Y.peaks$peakIndex <= endR, , drop = FALSE]
            
            
            ## A group is said to be a group when there are almost no duplicates in it
            
            realigned.peaks <- speaq2::hclust.grouping(current.peaks, min.samp.grp = min.samp.grp, max.dupli.prop = max.dupli.prop, 
                maxClust = maxClust, linkage = linkage)
            
            
            if (nrow(realigned.peaks) > 0) {
                alignindex.stop <- alignindex.counter + nrow(realigned.peaks)
                Y.aligned[alignindex.start:alignindex.stop, ] <- as.matrix(realigned.peaks)
                alignindex.counter <- alignindex.stop
                alignindex.start <- alignindex.stop + 1
                
                if (verbose == TRUE) 
                  print(paste("Realigned", nrow(realigned.peaks), "peaks"))
            }
        }
        
        utils::setTxtProgressBar(pb, maxPBindex - length(remaining.index))
    }
    close(pb)
    
    
    Y.aligned <- data.frame(Y.aligned)
    colnames(Y.aligned) <- colnames(Y.peaks)
    Y.aligned <- Y.aligned[stats::complete.cases(Y.aligned$peakIndex), ]
    
    
    ###### Verify realignment (check for faulty groupings where window splits occured)
    
    aligned.groupindexes <- unique(Y.aligned$peakIndex)
    window.breaks <- window.breaks[!is.na(window.breaks) & window.breaks >= min(aligned.groupindexes) & 
        window.breaks < max(aligned.groupindexes)]
    
    print("verifying realignment")
    maxPBindex <- length(window.breaks) + length(aligned.groupindexes)
    pb <- utils::txtProgressBar(min = 0, max = maxPBindex, style = 3)
    
    # append Y.aligned with some empty space to continuously add the reclustered peaks without rbind
    # (this would be slow)
    empty.start <- nrow(Y.aligned) + 1
    empty.data <- data.frame(matrix(NA, ncol = ncol(Y.aligned), nrow = nrow(Y.aligned)))
    colnames(empty.data) <- colnames(Y.aligned)
    Y.aligned <- rbind(Y.aligned, empty.data)
    
    
    Jac <- rep(NA, length(aligned.groupindexes))
    reclustered.groups.accumulator <- matrix(NA, nrow = nrow(Y.aligned), ncol = ncol(Y.aligned))
    strt <- 1
    stp <- 0
    for (k in 2:(length(aligned.groupindexes) - 1)) {
        
        dat.to.regroup <- Y.aligned[Y.aligned$peakIndex %in% aligned.groupindexes[c(k - 1, k, k + 1)], 
            ]
        regroup.indexes <- aligned.groupindexes[c(k - 1, k, k + 1)]
        
        if (nrow(dat.to.regroup) > 1) {
            
            Jac[k] <- sum(duplicated(dat.to.regroup$Sample))/length(unique(dat.to.regroup$Sample))
            
            if (Jac[k] <= Jaccard.regroup.threshold) {
                reclustered.groups <- speaq2::hclust.grouping(dat.to.regroup, min.samp.grp = min.samp.grp, max.dupli.prop = max.dupli.prop, 
                  maxClust = maxClust, linkage = linkage)
                
                Y.aligned[Y.aligned$peakIndex %in% regroup.indexes, ] <- NA
                Y.aligned <- Y.aligned[order(Y.aligned[, 1]), ]
                
                empty.start <- which(is.na(Y.aligned[, 1]))[1]
                stp <- empty.start + nrow(reclustered.groups) - 1
                Y.aligned[empty.start:stp, ] <- reclustered.groups
                
            }
        }
        utils::setTxtProgressBar(pb, k)
    }
    
    Y.aligned <- Y.aligned[stats::complete.cases(Y.aligned$peakIndex), ]
    
    
    # append Y.aligned with some empty space to continuously add the reclustered peaks without rbind
    # (this would be slow)
    empty.start <- nrow(Y.aligned) + 1
    Y.aligned <- rbind(Y.aligned, empty.data)
    
    aligned.groupindexes.update <- unique(Y.aligned$peakIndex[!is.na(Y.aligned$peakIndex)])
    Jaccard.index <- rep(NA, length(window.breaks))
    strt <- 1
    stp <- 0
    for (k in 1:length(window.breaks)) {
        
        # warnings are supressed since if min and max return 0 or Inf no further action is taken
        left.group <- suppressWarnings(max(aligned.groupindexes.update[aligned.groupindexes.update <= 
            window.breaks[k]]))
        right.group <- suppressWarnings(min(aligned.groupindexes.update[aligned.groupindexes.update > 
            window.breaks[k]]))
        
        if (left.group != 0 & right.group != 0 & is.finite(left.group) & is.finite(right.group)) {
            # calculate jaccard index
            Jaccard.index[k] <- sum(duplicated(Y.aligned$Sample[Y.aligned$peakIndex %in% c(left.group, 
                right.group)]))/length(unique(Y.aligned$Sample[Y.aligned$peakIndex %in% c(left.group, 
                right.group)]))
            
            if (Jaccard.index[k] <= Jaccard.regroup.threshold) {
                
                reclustered.groups <- speaq2::hclust.grouping(Y.aligned[Y.aligned$peakIndex %in% c(left.group, 
                  right.group), ], min.samp.grp = min.samp.grp, max.dupli.prop = max.dupli.prop, maxClust = maxClust, linkage = linkage)
                
                Y.aligned[Y.aligned$peakIndex %in% c(left.group, right.group), ] <- NA
                Y.aligned <- Y.aligned[order(Y.aligned[, 1]), ]
                
                empty.start <- which(is.na(Y.aligned[, 1]))[1]
                stp <- empty.start + nrow(reclustered.groups) - 1
                Y.aligned[empty.start:stp, ] <- reclustered.groups
                
                aligned.groupindexes.update <- unique(Y.aligned$peakIndex[!is.na(Y.aligned$peakIndex)])
                # Y.aligned[Y.aligned$peakIndex %in% c(left.group,right.group), ] = NA stp = empty.start +
                # nrow(reclustered.groups) - 1 Y.aligned[empty.start:stp, ] = reclustered.groups empty.start = stp +
                # 1 aligned.groupindexes.update = unique(Y.aligned$peakIndex[!is.na(Y.aligned$peakIndex)])
            }
        }
        
        utils::setTxtProgressBar(pb, k + length(aligned.groupindexes))
    }
    close(pb)
    
    Y.aligned <- Y.aligned[stats::complete.cases(Y.aligned$peakIndex), ]
    
    return(Y.aligned)
}


