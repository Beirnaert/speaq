#' BW and percentile ratios plot
#'
#' This function is used to plot BW and percentile ratios
#' 
#' @param BW An array of the BW ratios.
#' @param perc An array of the percentile ratios.
#' @param X The spectral dataset in matrix format in which each row contains a single sample.
#' @param startP The starting point of the segment. If it is -1, the starting point is from beginning of the spectra.
#' @param endP The ending point of the segment. If it is -1, the ending point is the last point of the spectra.
#' @param groupLabel The default value is NULL, it means that a single spectrum has a distinct color. Otherwise, the spectra is colored by their label.
#' @param highBound Default value is -1, that means the plot covers also the highest intensity peaks in the figure. If the users want to limit the upper height of the figure, set this parameter by the limited value.
#' @param lowBound Default value is -1, that means the plot covers also the lowest intensity peaks in the figure. If the users want to limit the under height of the figure, set this parameter by the limited value.
#' @param nAxisPos The number of ticks that will be displayed in the horizontal axis.
#' @param offside The offside of values in x-axis for display.
#'
#' @return Return a plot containing both the BW and the spectra.
#'
#' @author Trung Nghia Vu
#' 
#' @seealso \code{\link{drawSpec}}
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
#' resFindRef<- findRef(peakList);
#' refInd <- resFindRef$refInd;
#' maxShift = 50;
#' Y <- dohCluster(X,
#'                 peakList = peakList,
#'                 refInd = refInd,
#'                 maxShift  = maxShift,
#'                 acceptLostPeak = TRUE, verbose=FALSE);
#' # find the BW-statistic
#' BW = BWR(Y, groupLabel);
#' N = 100;
#' alpha = 0.05;
#' # create sampled H0 and export to file
#' H0 = createNullSampling(Y, groupLabel, N = N,verbose=FALSE)
#' #compute percentile of alpha
#' perc = double(ncol(Y));
#' alpha_corr = alpha/sum(returnLocalMaxima(Y[2,])$pkMax>50000);
#' for (i in seq_along(perc)) {
#'     perc[i] = quantile(H0[,i],1-alpha_corr, type = 3);
#' }
#' drawBW(BW, perc,Y, groupLabel = groupLabel)
#' 
#' @export
#' 
#' @importFrom graphics plot axis lines
#' 
drawBW <-function(BW, perc, X, startP=-1, endP=-1, groupLabel=NULL,
                  highBound=-1, lowBound=-1, nAxisPos=4, offside=0){
    if (startP==-1) startP=1;
    if (endP==-1) endP=ncol(X);
    GraphRange<-c(startP:endP); 
    op=par(mfrow=c(2,1))
    
    plot(BW[GraphRange],ylim=c(min(c(BW[GraphRange]),perc[GraphRange]),
                               max(c(BW[GraphRange]),perc[GraphRange])),type="n", ylab="BW",
         xlab="index", xaxt="n")
    tempVal =trunc(length(GraphRange)/nAxisPos);
    xPos=c(0:nAxisPos) * tempVal; 
    axis(1,at=xPos,labels=xPos+startP+offside);
    lines(BW[GraphRange],col= "blue")
    lines(perc[GraphRange],col= "black")
    drawSpec(X,xlab="index",ylab="intensity",startP=startP,endP=endP,groupLabel=groupLabel,
             offside=offside,highBound=highBound,lowBound=lowBound)
}
