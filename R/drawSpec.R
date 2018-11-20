#' Spectral plot
#'
#' This function allows to draw a segment or the whole spectra with limited high/low bounds of intensity.
#' 
#' @param X The spectral dataset in matrix format in which each row contains a single sample.
#' @param startP The starting point of the segment. If it is -1, the starting point is from beginning of the spectra.
#' @param endP The ending point of the segment. If it is -1, the ending point is the last point of the spectra.
#' @param groupLabel The default value is NULL, it means that a single spectrum has a distinct color. Otherwise, the spectra is colored by their label.
#' @param useLog The default value is -1, that means do not use a log transformation. If users want to transform the intensities to logarit values before plotting, set it to 1.
#' @param highBound Default value is -1, that means the plot covers also the highest intensity peaks in the figure. If the users want to limit the upper height of the figure, set this parameter by the limited value.
#' @param lowBound Default value is -1, that means the plot covers also the lowest intensity peaks in the figure. If the users want to limit the under height of the figure, set this parameter by the limited value.
#' @param xlab The default value is NULL, if so, "index" is displayed at the horizontal axis.
#' @param ylab The default value is NULL, if so, "intensity" is displayed at the vertical axis.
#' @param main The default value is NULL, if so, the title shows the values of startP and endP.
#' @param nAxisPos The number of ticks that you want to display in horizontal axis.
#' @param offside The offside of values in x-axis for display.
#' 
#' 
#' @return Return a plot of the spectra.
#'
#' @author Trung Nghia Vu
#' 
#' @seealso \code{\link{drawBW}}
#'
#' @examples
#'  res=makeSimulatedData();
#'  X=res$data;
#'  groupLabel=res$label;
#'  drawSpec(X)
#' 
#' @export
#' 
#' @importFrom graphics plot axis lines legend
#' 
drawSpec <-function(X, startP=-1, endP=-1, groupLabel=NULL, useLog=-1,
                    highBound=-1, lowBound=-1, xlab=NULL, ylab=NULL, main=NULL,
                    nAxisPos=4, offside=0)
{
    groupLabel_name=groupLabel;
    X=as.data.frame(X);
    colnames(X)=c(1:ncol(X));
    X=as.matrix(X);
    if (highBound!=-1){
        for (i in seq_len(nrow(X))){
            myIndex=which(X[i,]>highBound);
            X[i,myIndex]=highBound;
        }        
    }
    
    if (lowBound!=-1){
        for (i in seq_len(nrow(X))){
            myIndex=which(X[i,]<lowBound);
            X[i,myIndex]=lowBound;
        }        
    }
    
    if (is.null(groupLabel)){
        groupLabel=c(1:nrow(X));
        groupLabel=as.factor(groupLabel);
    }else{
        levels(groupLabel) =c(1:length(levels(groupLabel)))
    }
    
    if (startP==-1) startP=1;
    if (endP==-1) endP=ncol(X);
    
    if (is.null(xlab)){xlab='index' }
    if (is.null(ylab)){ylab='intensity' }
    if (is.null(main)){main=paste(' ',startP+offside,'-',endP+offside)}
    
    GraphRange<-c(startP:endP); 
    
    yn<-X[,GraphRange];
    if (useLog!=-1) yn=log(yn);
    
    plot(yn[1,],ylim=c(min(yn),max(yn)),type="n",
         ylab=ylab,
         xlab=xlab,
         main=main,
         xaxt="n" 
    )
    tempVal =trunc(length(GraphRange)/nAxisPos);
    xPos=c(0:nAxisPos) * tempVal; 
    axis(1,at=xPos,labels=xPos+startP+offside);
    
    for(i in seq_along(levels(groupLabel))) {
        groupLabelIdx=which(groupLabel==levels(groupLabel)[i]);
        for (j in seq_along(groupLabelIdx)){
            lines(yn[groupLabelIdx[j],],col=as.integer(levels(groupLabel)[i]))        
        }
    }
    
    if (!is.null(groupLabel_name)){
        legendPos="topleft";
        legend(legendPos,levels(groupLabel_name),col=as.integer(levels(groupLabel)),
               text.col="black",pch=c(19,19),bg='gray90')
    }
}
