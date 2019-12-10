#' Identify features (columns in the datamatrix) which are significantly associated with the outcome.
#'
#' This function produces a p-value for every column in the datamatrix, corresponding to the null hypothesis that outcome/response is independent of that feature. 
#'
#' @param datamatrix The data matrix with a column for each feature.
#' @param response A vector or matrix of outcomes/responses (e.g. class labels). the length of this vector or the amount of rows in this matrix should match the amount of rows in datamatrix.
#' @param p.adj The adjustment method for the p-values. Any of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH' (default), 'BY', 'fdr' or 'none' are accepted.
#' @param POI Only if 'response' is a matrix! The p values of interest. This is a number indicating which column of the 'response' matrix you are interested in. POI can range from 1 (default) to the number of columns in 'response'. 
#' @param responsevector (deprecated), please use the the more general 'response' variable instead.
#'
#' @return data with the features and their (adjusted) p-values, one for every column in the datamatrix  .
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#' 
#' @examples
#' nSamples <- 10
#' nFeatures <- 20
#' data.matrix <- matrix( stats::runif(n=nFeatures*nSamples, min=0,max=100), 
#' ncol = nFeatures, nrow = nSamples)
#' 
#' responseVec <- c( rep(0,nSamples/2), rep(1,nSamples/2) )
#' p_values <- relevant.features.p(datamatrix = data.matrix, response = 
#' responseVec, p.adj = 'none')
#' p_values_adjusted <- relevant.features.p( datamatrix = data.matrix, 
#' response = responseVec, p.adj = 'bonferroni')
#'
#' @export
#' 
#' @importFrom stats lm p.adjust var
#' 
relevant.features.p <- function(datamatrix, response, p.adj = "BH", POI = 1, responsevector = NULL) {
    
    ##### error checks #####
    if (!p.adj %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) 
        stop("p-value correction method is not one of the allowed possibilities")
    
    if(!is.null(responsevector)){
        warning("'responsevector' is deprecated, please use the more general 'response' variable.")
        response <- responsevector
    }
    
    if(is(response, "matrix")){
        if (nrow(response) != nrow(datamatrix)) {
            stop("The amount of rows of the response does not match with the amount of rows in the datamatrix")
        }
        if (POI > ncol(response)) {
            stop("The p-value of interest, 'POI', was larger than the number of columns in 'response'.")
        }
        if(POI == 0){
            stop("Mininmal POI = 1")
        }
        coeff.index <- POI + 1 # +1 because of the intercept
    } else{
        if (!is(response, "numeric")) {
            print("response was transformed to numeric")
            levels(response) <- 1:length(levels(response))
            response <- as.numeric(response)
        }
        
        if (length(response) != nrow(datamatrix)) {
            stop("The length of the response does not match with the amount of rows in the datamatrix")
        }
        coeff.index <- 2 # 2 because of the intercept
    }
    
    ##### the actual function #####
    nfeat <- ncol(datamatrix)
    
    p <- rep(NA, nfeat)
    for (coln in 1:nfeat) {
        if (stats::var(datamatrix[, coln]) == 0) {
            p[coln] <- 1
        } else {
            linreg <- stats::lm(datamatrix[, coln] ~ response)
            p[coln] <- summary(linreg)$coefficients[coeff.index, "Pr(>|t|)"]
        }
    }
    
    p.adjusted <- stats::p.adjust(p, method = p.adj, n = nfeat)
    if (!is.null(colnames(datamatrix))) {
        # this is the case when submitting a speaq 2 dataset where the colnames correspond to the
        # feature/peak indexes
        features.indexes <- as.numeric(colnames(datamatrix))
    } else {
        features.indexes <- seq(1, length(p.adjusted))
    }
    
    P.data <- data.frame(matrix(NA, ncol = 2, nrow = length(p.adjusted)))
    colnames(P.data) <- c("index", "p.values")
    P.data$index <- features.indexes
    P.data$p.values <- p.adjusted
    
    
    return(P.data)
    
}