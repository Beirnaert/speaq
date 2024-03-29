% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCANT.R
\name{SCANT}
\alias{SCANT}
\title{SCAle, Normalize and Transform a data matrix}
\usage{
SCANT(data.matrix, type = "unit", feature_orientation = "columns")
}
\arguments{
\item{data.matrix}{the data matrix to be scaled, normalized or transformed.}

\item{type}{the operations to be performed, this can be multiple and are performed sequentially. Any of 'unit', 'pareto', 'log10', 'log2', 'center', 'range', 'vast', or 'max' are accepted.}

\item{feature_orientation}{default = "columns". This corresponds to the default feature matrix with samples as rows and features as columns. The other option is "rows": samples as columns and different features as different rows.}
}
\value{
The scaled, normalized and/or transformed matrix.
}
\description{
This function allows the column-wise or row-wise scaling, normalization and transformation operations on a data matrix.
}
\examples{
Samples <- 10
Features <- 20
data.matrix <- matrix(runif(n=Features*Samples, min=0,max=100), 
ncol = Features, nrow = Samples) 

changed_matrix = SCANT(data.matrix, type=c('pareto', 'center'), feature_orientation = 'columns')

}
\references{
van den Berg RA, Hoefsloot HCJ, Westerhuis JA, et al. Centering, scaling, and transformations: improving the biological information content of metabolomics data. BMC Genomics 2006; 7:142.
}
\author{
Charlie Beirnaert, \email{charlie.beirnaert@uantwerpen.be}
}
