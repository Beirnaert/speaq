#' Build a raw data matrix (spectra) from spectra of unequal length
#'
#' This function can be used to build a data matrix from ill aligned spectra or of spectra of unequal length.
#' the result is a matrix whereby the first column matches (approximately) with a sinle left ppm value and the 
#' last column matches (approxiamtely) with a single right ppm value. Crucial is that the sample rates of the machine are the same
#' this should be always the case otherwise comparing intensities becomes meaningless.
#' Note that, as standard in NMR spectra, the highest ppm value is on the left
#'
#' @param spectrum.list A list of the spectra (y-values). Since by definition some of these differ in length this has to be in list form with a single spectrum per list item.
#' @param ppm.list The list of corresponding ppm values (x-values) with the highest ppm-value at the begining (left) as is the convention for NMR spectra. (This our ppm.edges.matrix has to be provided). 
#' @param ppm.edges.matrix The list with the starting and ending ppm values (highest ppm-value on the left/in the beginning). This or ppm.list has to be provided.
#'
#' @return SpectraAndPPM A list with 2 elements, the DataMatrix and the ppmMatrix.
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' # this is an example for 3 meaningless spectra
#' lengths_of_spectra <- c(100,150,120)
#' measurement_distance <- 0.01
#' starting_ppm_values <- c(8.7, 9.0, 9.0)
#' spectra <- list()
#' ppm_values <- list()
#' for (k in 1:3) {
#'     spectra[[k]] <- runif(lengths_of_spectra[k], min = 0, max = 10)
#'     
#'     # note the minus sign in the 'by' statement
#'     ppm_values[[k]] <- seq(from = starting_ppm_values[k], by = -measurement_distance, 
#'                            length.out = lengths_of_spectra[k])  
#' }
#' new.Data <- BuildRawDataMatrix(spectrum.list = spectra, ppm.list = ppm_values)
#' spectraMatrix <- new.Data$DataMatrix
#' ppmMatrix <- new.Data$ppmMatrix
#' 
#' @export
BuildRawDataMatrix <- function(spectrum.list, ppm.list = NULL, ppm.edges.matrix = NULL) {
    
    # This function can be used to build a data matrix from ill aligned spectra or of spectra of unequal
    # length.  the result is a matrix whereby the first column matches (approximately) with a sinle left
    # ppm value and the last column matches (approxiamtely) with a single right ppm value.  crucial is
    # that the sample rates of the machine are thesame this should be always the case otherwise comparing
    # intensities becomes meaningless
    
    # Note that, as standard in NMR spectra, the highest ppm value is on the left
    
    # checks and setting some variables
    if ("list" %in% class(spectrum.list)) {
        nspectra <- length(spectrum.list)
    } else stop("spectrum.list is not a list")
    
    
    if (is.null(ppm.list) & is.null(ppm.edges.matrix)) 
        stop(" at least 1 of 'ppm.list' or 'ppm.edges.matrix has to be provided'")
    
    if (is.null(ppm.list) & !is.null(ppm.edges.matrix)) {
        ppm.list <- NULL
        for (l in 1:nspectra) {
            ppm.list[[l]] <- seq(from = ppm.edges.matrix[l, 1], to = ppm.edges.matrix[l, 2], length.out = length(spectrum.list[[l]]))
        }
    }
    
    if (!is.null(ppm.list) & is.null(ppm.edges.matrix)) {
        ppm.edges.matrix <- matrix(NA, ncol = 2, nrow = nspectra)
        for (l in 1:nspectra) {
            ppm.edges.matrix[l, ] <- c(ppm.list[[l]][1], ppm.list[[l]][length(spectrum.list[[l]])])
        }
    }
    
    if (!all(ppm.edges.matrix[, 1] > ppm.edges.matrix[, 2])) 
        stop("some ppm vectors do not have the highest (most positive) value on the left")
    
    # the actual function
    lowest.left.index <- which.min(ppm.edges.matrix[, 1])
    lowest.left <- ppm.edges.matrix[lowest.left.index, 1]
    
    start.left.indexes <- rep(0, nspectra)
    for (s in 1:nspectra) {
        start.left.indexes[s] <- which.min(abs(ppm.list[[s]] - lowest.left))
    }
    
    points.left <- rep(0, nspectra)
    for (s in 1:nspectra) {
        points.left[s] <- length(ppm.list[[s]]) - start.left.indexes[s] + 1
    }
    
    Npoints <- min(points.left)
    
    DataMatrix <- matrix(NA, nrow = nspectra, ncol = Npoints)
    for (s in 1:nspectra) {
        DataMatrix[s, ] <- spectrum.list[[s]][start.left.indexes[s]:(start.left.indexes[s] + Npoints - 
            1)]
    }
    
    ppmMatrix <- matrix(NA, nrow = nspectra, ncol = Npoints)
    for (s in 1:nspectra) {
        ppmMatrix[s, ] <- ppm.list[[s]][start.left.indexes[s]:(start.left.indexes[s] + Npoints - 1)]
    }
    
    last.ppm.column <- ppmMatrix[, ncol(ppmMatrix)]
    end.ppm.diff <- max(last.ppm.column) - min(last.ppm.column)
    
    end.ppm.index.diff <- ceiling(end.ppm.diff/((ppm.list[[lowest.left.index]][1] - ppm.list[[lowest.left.index]][length(ppm.list[[lowest.left.index]])])/(length(ppm.list[[lowest.left.index]]))))
    
    print(paste("The ppm difference at the end is ", as.character(end.ppm.diff), " ppm. This is ", as.character(end.ppm.index.diff), 
        " measurement point(s)", sep = ""))
    
    SpectraAndPPM <- list()
    SpectraAndPPM$DataMatrix <- DataMatrix
    SpectraAndPPM$ppmMatrix <- ppmMatrix
    
    
    return(SpectraAndPPM)
}
