#' SCAle, Normalize and Transform a data matrix
#'
#' This function allows the column-wise or row-wise scaling, normalization and transformation operations on a data matrix.
#'
#' @param data.matrix the data matrix to be scaled, normalized or transformed.
#' @param type the operations to be performed, this can be multiple and are performed sequentially. Any of 'unit', 'pareto', 'log10', 'log2', 'center', 'range', 'vast', 'prob.Q' or 'max' are accepted.
#' @param feature_orientation default = "columns". This corresponds to the default feature matrix with samples as rows and features as columns. The other option is "rows": samples as columns and different features as different rows.
#' @param what (deprecated, use feature_orientation) to specify on which to perform the operations (row or column).
#'
#' @return The scaled, normalized and/or transformed matrix.
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' Samples <- 10
#' Features <- 20
#' data.matrix <- matrix(runif(n=Features*Samples, min=0,max=100), 
#' ncol = Features, nrow = Samples) 
#' 
#' changed_matrix = SCANT(data.matrix, type=c('pareto', 'center'), feature_orientation = 'columns')
#' 
#' @references van den Berg RA, Hoefsloot HCJ, Westerhuis JA, et al. Centering, scaling, and transformations: improving the biological information content of metabolomics data. BMC Genomics 2006; 7:142.
#' 
#' @export
#' 
#' @importFrom mQTL normalise
#' @importFrom stats sd
#' 
SCANT <- function(data.matrix, type = "unit", feature_orientation = "columns", what = NA) {
    
    for (N in seq_along(type)) {
        if (!(type[N] %in% c("unit", "pareto", "log10", "log2", "center", "range", "vast", "prob.Q", "max", "TotSum"))) {
            stop("No appropriate type of normalization/scaling selected")
        }
    }
    
    if(!is.na(what)){
        warning("The use of the 'what' parameter is deprecated. Please use 'feature_orientation' instead. For now feature_orientation has been set in accordance with the supplied 'what' parameter.")
        feature_orientation = what
    }
    
    if (!(feature_orientation %in% c("columns", "rows"))) {
        stop("No appropriate 'feature_orientation' is selected")
    }
    
    if (!("matrix" %in% class(data.matrix))) {
        data.matrix <- as.matrix(data.matrix)
        print("data.matrix was not a matrix, it has subsequently been converted. To make sure no problems arise, do this in advance.")
    }
    
    trans <- FALSE
    if (feature_orientation == "rows") {
        trans <- TRUE
        data.matrix <- t(data.matrix)
    }
    
    
    scaled.data <- data.matrix
    
    for (N in seq_along(type)) {
        Curr.type <- type[N]
        
        # standard ones
        if (Curr.type == "log10") {
            scaled.data[is.na(scaled.data)] <- 1  # for the log scaling
            scaled.data[scaled.data < 0.1] <- 0.1
            scaled.data <- log10(scaled.data)
        }
        if (Curr.type == "log2") {
            scaled.data[is.na(scaled.data)] <- 1  # for the log scaling
            scaled.data[scaled.data < 0.1] <- 0.1
            scaled.data <- log2(scaled.data)
        }
        if (Curr.type == "unit") {
            scaled.data <- scale(scaled.data, center = FALSE, scale = TRUE)
        }
        if (Curr.type == "center") {
            scaled.data <- scale(scaled.data, center = TRUE, scale = FALSE)
        }
        # from R van den Berg, Centering, scaling, and transformations: improving the biological information
        # content of metabolomics data
        if (Curr.type == "pareto") {
            scaled.data <- apply(scaled.data, 2, function(x) x/sqrt(stats::sd(x)))
        }
        if (Curr.type == "range") {
            scaled.data <- apply(scaled.data, 2, function(x) (x - mean(x))/(max(x) - min(x)))
        }
        if (Curr.type == "vast") {
            scaled.data <- apply(scaled.data, 2, function(x) ((x - mean(x)) * mean(x))/(stats::sd(x)^2))
        }
        ## Quotient probabilistic normalisation
        if (Curr.type == "prob.Q") {
            scaled.data <- mQTL::normalise(scaled.data, "prob")
            scaled.data <- scaled.data[[1]]
        }
        if (Curr.type == "max") {
            scaled.data <- apply(scaled.data, 2, function(x) (x/max(x)))
        }
        if (Curr.type == "TotSum"){
            scaled.data <-apply(scaled.data, 1, function(x) (x/sum(x)))
        }
        
    }
    
    if (trans) {
        scaled.data <- t(scaled.data)
    }
    
    return(scaled.data)
    
}
