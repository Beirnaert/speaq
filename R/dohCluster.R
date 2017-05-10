#' CluPA function for multiple spectra.
#'
#' Use CluPA for alignment for multiple spectra.
#'
#'
#' @param X The spectral dataset in the matrix format in which each row contains a single sample
#' @param peakList The peak lists of the spectra
#' @param refInd The index of the reference spectrum.
#' @param maxShift The maximum number of the points for a shift step.
#' @param acceptLostPeak This is an option for users, TRUE is the default value. If the users believe that all the peaks in the peak list are true positive, change it to FALSE.
#' @param verbose A boolean value to allow print out process information.
#'
#' @return The aligned spectra.
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
#' );
#' resFindRef<- findRef(peakList);
#' refInd <- resFindRef$refInd;
#' maxShift = 50;
#' Y <- dohCluster(X,
#'                 peakList = peakList,
#'                 refInd = refInd,
#'                 maxShift  = maxShift,
#'                 acceptLostPeak = TRUE, verbose=FALSE);
#' 
#' @export
#' 
#' @importFrom graphics plot points
#' 
dohCluster <-function(X, peakList, refInd=0, 
                      maxShift =100, acceptLostPeak=TRUE, verbose=TRUE){
    
    .withMaxShift <-function(X, peakList, refInd=0, 
                             maxShift =100, acceptLostPeak=TRUE, verbose=TRUE){
        
        Y=X;
        peakListNew=peakList;
        
        if (verbose) startTime=proc.time();  
        refSpec=Y[refInd,];        
        for (tarInd in 1: nrow(X))
            if (tarInd!=refInd)
            {
                if (verbose) cat("\n aligning spectrum ",tarInd);
                targetSpec=Y[tarInd,];
                myPeakList=c(peakList[[refInd]],peakList[[tarInd]]);
                myPeakLabel=double(length(myPeakList));
                for (i in 1:length(peakList[[refInd]]) ) myPeakLabel[i]=1;
                startP=1;
                endP=length(targetSpec);
                res=hClustAlign(refSpec,targetSpec,myPeakList,myPeakLabel,startP,endP,
                                maxShift=maxShift,acceptLostPeak=acceptLostPeak)
                Y[tarInd,]=res$tarSpec;
                
                peakListNew[[tarInd]]=
                    res$PeakList[(length(peakList[[refInd]])+1):length(myPeakList)]
            }
        
        if (verbose){
            Z=stats::cor(t(Y))
            newCor= stats::median(Z[lower.tri(Z)])
            cat("\n Median pearson correlation of aligned spectra:",newCor)
            
            endTime=proc.time();
            cat("\n Alignment time: ",(endTime[3]-startTime[3])/60," minutes");
        }
        return(Y);
    }
    
    .withoutMaxShift <- function (X, peakList, refInd = 0, acceptLostPeak = TRUE, 
                                  verbose = TRUE) 
    {
        if (verbose) 
            startTime = proc.time()
        
        maxShift_ladder=2^(c(1:trunc(log2(ncol(X)/2))))
        bestCor=-1
        corVec=NULL;
        bestY=NULL
        bestMaxShift=0
        for (maxShift_val in maxShift_ladder){
            
            if (verbose) cat("\n maxShift=",maxShift_val)
            Y = X
            peakListNew = peakList
            refSpec = Y[refInd, ]
            for (tarInd in 1:nrow(X)) if (tarInd != refInd) {
                #if (verbose) cat("\n aligning spectrum ", tarInd)
                targetSpec = Y[tarInd, ]
                myPeakList = c(peakList[[refInd]], peakList[[tarInd]])
                myPeakLabel = c(rep(1,length(peakList[[refInd]])),rep(0,length(peakList[[tarInd]])))
                
                startP = 1
                endP = length(targetSpec)
                res = hClustAlign(refSpec, targetSpec, myPeakList, myPeakLabel, 
                                  startP, endP, maxShift = maxShift_val, acceptLostPeak = acceptLostPeak)
                Y[tarInd, ] = res$tarSpec
                peakListNew[[tarInd]] = res$PeakList[(length(peakList[[refInd]]) + 
                                                          1):length(myPeakList)]
            }
            
            Z=stats::cor(t(Y))
            newCor= stats::median(Z[lower.tri(Z)])
            corVec=c(corVec,newCor);
            if (verbose) cat("\n Median Pearson correlation coefficent:",newCor,", the best result:",bestCor)
            if (newCor>bestCor){
                bestCor=newCor
                bestY=Y
                bestMaxShift=maxShift_val
            }
            
        }
        
        if (verbose){ 
            cat("\nOptimal maxShift=",bestMaxShift,"with median Pearson correlation of aligned spectra=",bestCor);
            plot(log2(maxShift_ladder),corVec,type="b",  xlab="log2(maxShift)",ylab="Median Pearson correlation coefficent", 
                 main=paste("Optimal maxShift=",bestMaxShift," (red star) \n with median Pearson correlation coefficent of ",round(bestCor,6), sep=""))
            points(log2(bestMaxShift),bestCor, col="red",pch=8, cex=2.0)
            
        }    
        if (verbose) {
            endTime = proc.time()
            cat("\n Alignment time: ", (endTime[3] - startTime[3])/60, 
                " minutes")
        }
        return(bestY)
    }
    
    
    if (is.null(maxShift)){
        if (verbose) {
            cat("\n --------------------------------")
            cat("\n maxShift=NULL, thus CluPA will automatically detect the optimal value of maxShift.")
            cat("\n --------------------------------\n")
        }
        res=.withoutMaxShift(X=X, peakList=peakList, refInd=refInd, 
                             acceptLostPeak=acceptLostPeak, verbose=verbose)
    }
    else{
        if (verbose) {
            cat("\n --------------------------------")
            cat("\n CluPA will run with maxShift=",maxShift)
            cat("\n If you want CluPA automatically detect the optimal maxShift,")
            cat("\n let try with dohCluster(...,maxShift=NULL,..)")
            cat("\n --------------------------------\n")
        }
        res=.withMaxShift(X=X, peakList=peakList, refInd=refInd, 
                          maxShift=maxShift, acceptLostPeak=acceptLostPeak, verbose=verbose)
    }
    
    return(res)
}
