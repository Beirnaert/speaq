#' Selecting the peaks in a segment
#'
#' This function is to find out which peaks belonging to a segment which ranges from startP to endP
#' 
#' @param peakList The peak lists of the spectra.
#' @param startP The starting point of the segment.
#' @param endP The ending point of the segment.
#' 
#' @return The list of indices of the peaks in the segment.
#'
#' @author Trung Nghia Vu
#' 
#' @seealso \code{\link{dohClusterCustommedSegments}}
#' 
#' @examples
#' res=makeSimulatedData();
#' X=res$data;
#' groupLabel=res$label;
#' peakList <- detectSpecPeaks(X,
#'                             nDivRange = c(128),
#'                             scales = seq(1, 16, 2),
#'                             baselineThresh = 50000,
#'                             SNR.Th = -1,
#'                             verbose=FALSE
#'                             );
#' cat("\n ", peakList[[1]])
#' segmentpeakList= findSegPeakList(peakList[[1]],400,600);
#' cat("\n ", segmentpeakList)     
#'                        
#' @export
#' 
findSegPeakList <-function(peakList, startP, endP){
    res=0;  
    for (i in 1:length(peakList)){    
        if (peakList[i]>startP&&peakList[i]<endP){    
            res=c(res,peakList[i]-startP+1);
        }    
    }
    return(res[-1])
}
