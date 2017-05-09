#' SilhouetR
#'
#' Calculate the Silhouette value for 
#'
#' @param DataMatrix a matrix with the raw data, 1 variable per column.
#' @param GroupIndices The vector with the group indices (length must be equal to the amount of rwos in DataMatrix).
#' @param distance The distance metric to be used, see \link[cluster]{daisy}.
#' @param stand whether to standardize the data before calculating the dissimilarities. See \link[cluster]{daisy}.
#'
#' @return Returns the silhouette values.
#' 
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' \dontrun{
#' # This function works on grouped peak data
#' library(ggplot2)
#' Silhouette.values = SilhouetR(DataMatrix = Grouped.peaks$peakPPM, 
#' Grouped.peaks$peakIndex, distance = 'euclid', stand = TRUE)
#' ggplot(SilhouetteValues, aes(SilhouetteValues)) +geom_freqpoly(binwidth = 0.03) +theme_bw()
#' }
#' 
#' @export
#' 
#' @importFrom cluster daisy
#' @importFrom utils head tail txtProgressBar
#' 
SilhouetR <- function(DataMatrix, GroupIndices, distance = "euclid", stand = TRUE) {
    
    if (!"matrix" %in% class(DataMatrix) & !"data.frame" %in% class(DataMatrix)) {
        DataMatrix <- matrix(DataMatrix, ncol = 1)
        message("DataMatrix is not a matrix, attempting conversion with the assumption of only 1 variable (1 column)")
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
    
    start_stop_indices <- matrix(NA, ncol = 2, nrow = Ngroups)
    for (k in 1:Ngroups) {
        start_stop_indices[k, 1] <- utils::head(which(GroupIndices == groups[k]), n = 1)
        start_stop_indices[k, 2] <- utils::tail(which(GroupIndices == groups[k]), n = 1)
    }
    
    
    
    print("Computing dissimilarity matrix")
    # dissimilarity matrix
    dissSim <- as.matrix(cluster::daisy(DataMatrix, metric = distance, stand = stand))
    
    print("Computing silhouette values")
    pb <- utils::txtProgressBar(min = 0, max = Nindividuals, style = 3, width = 100)
    
    A <- NULL
    B <- NULL
    SilhouetteValues <- rep(NA, Nindividuals)
    for (i in 1:Nindividuals) {
        ABvalues <- rep(NA, Ngroups)
        for (grp in 1:Ngroups) {
            # ABvalues[grp] = mean(dissSim[i,GroupIndices == groups[grp] & individual_indices!=i])
            ABvalues[grp] <- mean(dissSim[i, start_stop_indices[grp, 1]:start_stop_indices[grp, 2]])
        }
        individual_grp <- which(groups == GroupIndices[i])
        # A = ABmatrix[i,individual_grp] B = min(ABmatrix[i,-individual_grp])
        ABvalues[individual_grp] <- mean(dissSim[i, GroupIndices == GroupIndices[i] & individual_indices != 
            i])
        A <- ABvalues[individual_grp]
        B <- min(ABvalues[-individual_grp])
        SilhouetteValues[i] <- (B - A)/max(c(A, B))
        utils::setTxtProgressBar(pb, i)
    }
    close(pb)
    SilhouetteValues.df <- data.frame(matrix(c(individual_indices, SilhouetteValues, GroupIndices), ncol = 3, 
        byrow = FALSE))
    colnames(SilhouetteValues.df) <- c("index", "SilhouetteValues", "GroupIndices")
    
    SilhouetteValues.df = SilhouetteValues.df[order(original.order),]
    
    return(SilhouetteValues.df)
}
