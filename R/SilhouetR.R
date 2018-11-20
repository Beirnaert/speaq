#' SilhouetR
#'
#' This function calculate Silhouette values. The function is generic, 
#' as such silhouette values can be calculated between samples of different classes or 
#' it can be used to calculate silhouette values between different groups of peaks. 
#' This is the way in which it is used for the speaq package (see the example).
#'
#' @param DataMatrix a matrix with the raw data, 1 variable per column.
#' @param GroupIndices The vector with the group indices (length must be equal to the amount of rows in DataMatrix).
#' @param distance The distance metric to be used, "euclidean" or "manhattan".
#'
#' @return Returns the silhouette values. Note if a group contains only 1 no Silhouette value can be calculated (will give NA)
#' 
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' subset <- GetWinedata.subset()
#' # to reduce the example time we only select spectra 1 & 2
#' subset.spectra = as.matrix(subset$Spectra)[1:2,] 
#' subset.ppm = as.numeric(subset$PPM)
#' 
#' test.peaks <- getWaveletPeaks(Y.spec=subset.spectra, 
#'                               X.ppm=subset.ppm,
#'                               nCPU = 1) # nCPU set to 2 for the vignette build
#'                               
#' test.grouped <- PeakGrouper(Y.peaks = test.peaks)
#'
#' Silhouette.values = SilhouetR(DataMatrix = test.grouped$peakPPM, 
#'                               test.grouped$peakIndex, 
#'                               distance = "euclidean")
#'                               
#' hist(Silhouette.values$SilhouetteValues)
#' 
#' 
#' @export
#' 
#' @importFrom Rfast dista
#' @importFrom utils head tail txtProgressBar
#' 
SilhouetR <- function(DataMatrix, GroupIndices, distance = "euclidean") {
    
    if (!"matrix" %in% class(DataMatrix) & !"data.frame" %in% class(DataMatrix)) {
        DataMatrix <- matrix(DataMatrix, ncol = 1)
        message("DataMatrix is not a matrix, attempting conversion with the assumption of only 1 variable (1 column)")
    }
    
    if(any(is.na(GroupIndices))){
        message("GroupIndices contains NA values. These will be deleted together with the corresponding DataMatrix entries.")
        DataMatrix = DataMatrix[!is.na(GroupIndices), , drop = FALSE]
        GroupIndices = GroupIndices[!is.na(GroupIndices)]
    }
    
    # order the Datamatrix according to increasing GroupIndex, this is to speed up the code later on by
    # calling indices instead of comparing to group
    DataMatrix <- DataMatrix[order(GroupIndices), , drop = FALSE]
    original.order = seq(1,nrow(DataMatrix))[order(GroupIndices)]
    GroupIndices <- GroupIndices[order(GroupIndices)]
    
    groups <- unique(GroupIndices)
    Ngroups <- length(groups)
    Nindividuals <- nrow(DataMatrix)
    individual_indices <- 1:Nindividuals
    
    
    if(Ngroups > 2){
        start_stop_indices <- matrix(NA, ncol = 2, nrow = Ngroups)
        start_stop_indices[1,1] <- 1
        start_stop_indices[Ngroups,2] <- length(GroupIndices)
        start_stop_indices[1:(Ngroups-1),2] <- which(diff(GroupIndices)>0)
        start_stop_indices[2:Ngroups,1] <- which(diff(GroupIndices)>0) + 1
    } else{
        start_stop_indices <- matrix(NA, ncol = 2, nrow = Ngroups)
        for (k in 1:Ngroups) {
            start_stop_indices[k, 1] <- utils::head(which(GroupIndices == groups[k]), n = 1)
            start_stop_indices[k, 2] <- utils::tail(which(GroupIndices == groups[k]), n = 1)
        }
    }
    
    print("Computing silhouette values")
    pb <- utils::txtProgressBar(min = 0, max = Nindividuals, style = 3, width = 100)
    
    A <- NULL
    B <- NULL
    SilhouetteValues <- rep(NA, Nindividuals)
    for (i in 1:Nindividuals) {
        ABvalues <- rep(NA, Ngroups)
        if(ncol(DataMatrix) == 1){
            dissSim <- abs(DataMatrix - DataMatrix[i,])
        } else{
            dissSim <- Rfast::dista(xnew = DataMatrix[i, , drop = FALSE], x = DataMatrix, type = distance )
        }
        for (grp in 1:Ngroups) {
            ABvals <- dissSim[start_stop_indices[grp, 1]:start_stop_indices[grp, 2]]
            ABvalues[grp] <- sum(ABvals)/length(ABvals)
        }
        individual_grp <- which(groups == GroupIndices[i])
        # A = ABmatrix[i,individual_grp] B = min(ABmatrix[i,-individual_grp])
        ABvalues[individual_grp] <- mean(dissSim[GroupIndices == GroupIndices[i] & individual_indices != 
                                                     i])
        A <- ABvalues[individual_grp]
        B <- min(ABvalues[-individual_grp])
        SilhouetteValues[i] <- (B - A)/max(c(A, B))
        utils::setTxtProgressBar(pb, i)
    }
    close(pb)
    SilhouetteValues.df <- data.frame(matrix(c(original.order[individual_indices], SilhouetteValues, GroupIndices), ncol = 3, 
                                             byrow = FALSE))
    colnames(SilhouetteValues.df) <- c("index", "SilhouetteValues", "GroupIndices")
    
    SilhouetteValues.df$SilhouetteValues[is.nan(SilhouetteValues.df$SilhouetteValues)] = NA
    
    SilhouetteValues.df = SilhouetteValues.df[order(original.order),]
    
    return(SilhouetteValues.df)
    
}