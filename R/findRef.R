#' Reference finding
#'
#' This function is to heuristically detect a reference spectrum.
#' 
#' @param peakList The peak lists of the spectra.
#' 
#' 
#' @return list of 2: refInd (The index of the reference spectrum found by the algorithm) and orderSpec (A sorted array of the spectra by their goodness values)
#'
#' @author Trung Nghia Vu
#' 
#' @references Vu TN, Valkenborg D, Smets K, Verwaest KA, Dommisse R, Lemi\`{e}re F, Verschoren A, Goethals B, Laukens K. (2011) An integrated workflow for robust alignment and simplified quantitative analysis of NMR spectrometry data. BMC Bioinformatics. 2011 Oct 20;12:405.
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
#' );
#' cat("\n Find the spectrum reference...")
#' resFindRef<- findRef(peakList);
#' refInd <- resFindRef$refInd;
#' cat("\n Order of spectrum for reference  \n");
#' for (i in seq_along(resFindRef$orderSpec))
#'     cat(paste(i, ":",resFindRef$orderSpec[i],sep=""), " ");
#' cat("\n The reference is: ", refInd);
#' 
#' @export
#' 
findRef <-function(peakList)
{
    disS=matrix(data=NA,ncol=length(peakList),nrow=length(peakList));
    sumDis=double(length(peakList));
    for(refInd in seq_along(peakList)){
        for(tarInd in seq_along(peakList))
            if (refInd!=tarInd)
            {
                disS[refInd,tarInd]=0;
                for (i in seq_along(peakList[[tarInd]]))
                    disS[refInd,tarInd]=disS[refInd,tarInd]+
                        min(abs(peakList[[tarInd]][i]-peakList[[refInd]]));
            }
    }
    
    for(refInd in seq_along(peakList)) {
        disS[refInd,refInd]=0;
        sumDis[refInd]=sum(disS[refInd,]);
    }
    orderSumdis=order(sumDis);
    return(list(refInd=orderSumdis[1],orderSpec=orderSumdis));
}
