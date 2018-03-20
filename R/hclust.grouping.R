#' Grouping with hierarchical clustering (used in the PeakGrouper function)
#'
#' Internal function in the PeakGrouper function for generating the hierarchical clustering tree and cutting it.
#'
#' @param current.peaks A number of neighbouring peaks to be grouped.
#' @param min.samp.grp The minimal amount of samples needed to form a group.
#' @param max.dupli.prop The maximal duplication proportion allowed for a group to be considered a single group.
#' @param maxClust The maximum number of clusters (depth of the tree).
#' @param linkage The linkage to be used in the hierarchical clustering. See the 'method' argument in \link[stats]{hclust}.
#'
#' @return Returns a data frame with grouped peaks.
#' 
#' 
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @seealso \code{\link{PeakGrouper}}
#' 
#' @export
#' 
#' @importFrom cluster daisy
#' @importFrom stats hclust cutree kmeans
#' 
hclust.grouping <- function(current.peaks, min.samp.grp = 1, max.dupli.prop = 0.25, maxClust = 10, linkage = "average") {
    
    if (nrow(current.peaks) <= 1) {
        regrouped.peaks <- current.peaks
    } else {
        
        if (nrow(current.peaks) < maxClust) {
            maxClust <- nrow(current.peaks) - 1
        }
        
        dissE <- cluster::daisy(matrix(c(current.peaks$peakPPM), ncol = 1), metric = "euclid", stand = TRUE, warnType = FALSE)  # the variables are not standardized 
        
        fit <- stats::hclust(dissE, method = linkage)
        # silhouet = numeric(maxClust) silhouet[1]= 0
        duplicates <- numeric(maxClust)
        duplicates[1] <- sum(duplicated(current.peaks$Sample))
        groupassignment <- matrix(NA, ncol = maxClust + 1, nrow = length(current.peaks$Sample))
        grpcounter <- 1
        if (sum(duplicated(current.peaks$Sample))/length(current.peaks$Sample) < max.dupli.prop) {
            groupassignment[, 1] <- grpcounter
            grpcounter <- grpcounter + 1
        }
        
        for (k in 2:maxClust) {
            
            groups <- stats::cutree(fit, k)
            # sk <- silhouette(groups, dissE) silhouet[k] = mean(sk[,3])
            dup.sub.prop <- numeric(k)
            sub.size <- numeric(k)
            for (kk in 1:k) {
                dup.sub.prop[kk] <- sum(duplicated(current.peaks$Sample[groups == kk]))/length(current.peaks$Sample[groups == 
                  kk])
                sub.size[kk] <- length(current.peaks$Sample[groups == kk])
                if (dup.sub.prop[kk] <= max.dupli.prop & sub.size[kk] >= min.samp.grp) {
                  groupassignment[groups == kk, k] <- grpcounter
                  grpcounter <- grpcounter + 1
                }
            }
            duplicates[k] <- mean(dup.sub.prop)
        }
        
        
        
        for (s in 1:length(current.peaks$Sample)) {
            
            groupassignment[s, maxClust + 1] <- suppressWarnings(min(groupassignment[s, 1:maxClust], 
                na.rm = TRUE))
            
        }
        
        # the group with 0 as its group value (if there is one) is the garbage group
        group <- rep(0, length(current.peaks$Sample))
        grps <- unique(groupassignment[, maxClust + 1])
        grps <- grps[!is.infinite(grps)]
        if (length(grps) > 0) {
            for (cl in 1:length(grps)) {
                group[groupassignment[, maxClust + 1] == grps[cl]] <- cl
            }
            ngrps <- max(group)
            
            
            
            starts <- numeric(ngrps)
            stops <- numeric(ngrps)
            for (g in 1:ngrps) {
                starts[g] <- min(current.peaks$peakPPM[group == g])
                stops[g] <- max(current.peaks$peakPPM[group == g])
            }
            
            
            
            
            # remove duplicates
            
            for (gg in 1:ngrps) {
                
                if (any(duplicated(current.peaks$Sample[group == gg]))) {
                  duplicated.samples <- unique(current.peaks$Sample[group == gg][duplicated(current.peaks$Sample[group == 
                    gg])])  # get the samples that have 1 or more duplicated peaks
                  group.meanppm <- mean(current.peaks$peakPPM[group == gg][!current.peaks$Sample[group == 
                    gg] %in% duplicated.samples])  # get the average ppm for all the samples without duplicates
                  group.meanSNR <- mean(current.peaks$peakSNR[group == gg][!current.peaks$Sample[group == 
                    gg] %in% duplicated.samples])  # get the average SNR for all the samples without duplicates
                  
                  for (ds in 1:length(duplicated.samples)) {
                    
                    dat.mat <- matrix(c(group.meanppm, current.peaks$peakPPM[group == gg & current.peaks$Sample == 
                      duplicated.samples[ds]], group.meanSNR, current.peaks$peakSNR[group == gg & current.peaks$Sample == 
                      duplicated.samples[ds]]), ncol = 2, byrow = F)
                    nduplis <- nrow(dat.mat) - 1
                    peak.to.keep <- which.min(cluster::daisy(dat.mat, metric = "gower", stand = FALSE)[1:nduplis])  # the peak with the minimal distance to the group mean values is the one to keep
                    to.remove <- seq(1, nduplis)
                    to.remove <- to.remove[-peak.to.keep]
                    group[group == gg & current.peaks$Sample == duplicated.samples[ds]][to.remove] <- rep(0, 
                      length(to.remove))
                    
                  }
                  
                }
                
            }
            
            #### see that no 2 groups overlap ### do this by first checking if there is an overlap somewhere (it
            #### will usually not be the case) If there is an overlap, do an extensive search to find the overlap
            #### check steps: sort the starting values from low to high and sort the stop values acording to the
            #### starting values check if stop j is smaller than start j+1. If this is the case than there is no
            #### overlap
            group.order <- order(starts)
            stops <- stops[group.order]
            starts <- starts[group.order]
            
            
            if (any((starts[2:ngrps] - stops[1:ngrps - 1]) < 0)) {
                overlapping.groups <- which((starts[2:ngrps] - stops[1:ngrps - 1]) < 0)
                all.overlapping.groups <- matrix(c(overlapping.groups, overlapping.groups + 1), ncol = 2)
                # note that the way this is implemented multiple groups can overlap and they will all be merged
                for (jj in nrow(all.overlapping.groups):1) {
                  group[group == all.overlapping.groups[jj, 2]] <- all.overlapping.groups[jj, 1]
                }
                ngrps <- ngrps - nrow(all.overlapping.groups)
                # check if the groups are stil sequentially ordered, if not do this
                if (any(!(seq(1, ngrps) %in% group))) {
                  curr.grps <- sort(unique(group[group > 0]))
                  for (gg in 1:length(curr.grps)) {
                    group[group == curr.grps[gg]] <- gg
                  }
                }
                
            }
            
            
            
            regrouped.peaks <- current.peaks[group > 0, , drop = FALSE]
            group <- group[group > 0, drop = FALSE]
            
            if (nrow(regrouped.peaks) > 0) {
                for (g in 1:max(group)) {
                  theoretical.centerIndex <- round(as.numeric(stats::kmeans(regrouped.peaks$peakIndex[group == 
                    g], centers = 1)[2]))
                  regrouped.peaks$peakIndex[group == g] <- theoretical.centerIndex
                }
                regrouped.peaks <- regrouped.peaks[order(regrouped.peaks$peakIndex), ]
            }
        } else {
            regrouped.peaks <- matrix(NA, ncol = 6, nrow = 0)
        }
    }
    
    return(regrouped.peaks)
}

