#' Regroup faulty grouped peaks
#'
#' If there are peaks wrongly grouped by the peakGrouper function, they will be regrouped by using the ppm values together with the peak signal to noise ratio.
#'
#' @param grouped.peaks The grouped peaks data.
#' @param list.to.regroup The peak indices of groups to regroup (the groups, indicated by their peakIndex, in 1 list item will be merged and regrouped).
#' @param min.samp.grp The minimal amount of samples needed to form a group.
#' @param max.dupli.prop The maximal duplication proportion allowed for a group to be considered a single group.
#' @param maxClust The maximum number of clusters (depth of the tree).
#'
#' @return Returns a data frame with regrouped peaks.
#' 
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @export
#' 
#' @importFrom cluster daisy
#' @importFrom stats hclust cutree kmeans
#' @importFrom data.table rbindlist
#' 
regroupR <- function(grouped.peaks, list.to.regroup, min.samp.grp = 1, max.dupli.prop = 0.1, maxClust = 10) {
    
    
    if (class(list.to.regroup) != "list") {
        message("'list.to.regroup' is not a list. Assuming it is a list with only one element.")
        list.to.regroup <- list(list.to.regroup)
    }
    if (length(list.to.regroup) == 0) {
        stop("list.to.regroup is empty")
    }
    
    all.regrouped.peaks <- list()
    list.counter <- 0
    for (iter in 1:length(list.to.regroup)) {
        
        Indexes.to.regroup <- list.to.regroup[[iter]]
        
        current.peaks <- grouped.peaks[grouped.peaks$peakIndex %in% Indexes.to.regroup, ]
        grouped.peaks[grouped.peaks$peakIndex %in% Indexes.to.regroup, ] <- NA
        
        if (nrow(current.peaks) < maxClust) {
            maxClust <- nrow(current.peaks) - 1
        }
        
        dissE <- cluster::daisy(matrix(c(current.peaks$peakPPM, current.peaks$peakSNR), ncol = 2, byrow = FALSE), 
            metric = "gower", stand = TRUE)  # the variables are not standardized 
        
        fit <- stats::hclust(dissE, "average")
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
        
        if (nrow(regrouped.peaks) > 0) {
            list.counter <- list.counter + 1
            all.regrouped.peaks[[list.counter]] <- regrouped.peaks
        }
        
    }
    
    if (length(all.regrouped.peaks) > 0) {
        all.regrouped.peaks <- data.table::rbindlist(all.regrouped.peaks)
    }
    
    grouped.peaks <- grouped.peaks[complete.cases(grouped.peaks$peakIndex), ]
    grouped.peaks <- rbind(grouped.peaks, all.regrouped.peaks)
    grouped.peaks <- grouped.peaks[order(grouped.peaks$peakIndex), ]
    
    return(grouped.peaks)
}
