#' Use CluPA for alignment with additional information
#'
#' This function integrates some additional information from user such as references for each specific segment, segment ignorance, maximum step size.. to align spectra using CluPA.
#'
#' @details 
#' Each row of the segmentInfoMat matrix includes 5 values. 
#' For example, it could be imported from a CSV file consisting of following content:
#' #
#' begin,end,forAlign,ref,maxShift
#' 100,200,0,0,0
#' 450,680,1,0,50
#' #
#' Each column could be explained as the following:
#' * begin: the starting point of the segment.
#' * end: the end point of the segment.
#' * forAlign: the segment is aligned (1) or not (0).
#' * ref: the index of the reference spectrum. If 0, the algorithm will select the reference found by the reference finding step.
#' * maxShift: the maximum number of points of a shift to left/right.
#' It is worth to note that only segments with forAlign=1 (column 3) will be taken into account for spectral alignment.
#'
#' @param X The spectral dataset in the matrix format in which each row contains a single sample
#' @param peakList The peak lists of the spectra
#' @param refInd The index of the reference spectrum.
#' @param segmentInfoMat The matrix containing the additional information for segments from the users. This parameter must be a matrix.
#' @param minSegSize The minimum size of the segments which could be considered for alignment.
#' @param maxShift The maximum number of the points for a shift step.
#' @param acceptLostPeak This is an option for users, TRUE is the default value. If the users believe that all the peaks in the peak list are true positive, change it to FALSE.
#' @param verbose A boolean value to allow print out process information.
#'
#' @return The aligned spectral segments.
#'
#' @author Trung Nghia Vu
#' 
#' @seealso \code{\link{dohCluster}}
#'
#' @examples
#' cat("\n Please see more examples in the vignettes file.")
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
#' segmentInfoMat=matrix(data=c(100,200,0,0,0,
#'                              50,680,1,0,50),nrow=2,ncol=5,byrow=TRUE
#' )
#' colnames(segmentInfoMat)=c("begin","end","forAlign","ref","maxShift")
#' segmentInfoMat
#' maxShift = 50;
#' Yc <- dohClusterCustommedSegments(X,
#'                                   peakList = peakList,
#'                                   refInd = refInd,
#'                                   maxShift  = maxShift,
#'                                   acceptLostPeak = TRUE,
#'                                   segmentInfoMat = segmentInfoMat,
#'                                   minSegSize = 128,
#'                                   verbose=FALSE)
#' 
#' @export
#' 
dohClusterCustommedSegments <-function(X, peakList, refInd, segmentInfoMat, 
                                       minSegSize=128, maxShift=100, acceptLostPeak=TRUE, verbose=TRUE){
    if (!is.matrix(segmentInfoMat)) {
        cat("ERROR! segmentInfoMat must be in a matrix format.")
        return(NULL)
    }
    
    if (segmentInfoMat[1,1]>minSegSize) 
        mysegments=c(1,segmentInfoMat[1,1]-1,0,0,0) else mysegments=c();
        i = 0;
        if (nrow(segmentInfoMat)>1){   
            for (i in 1:(nrow(segmentInfoMat)-1)){
                mysegments=c(mysegments,c(segmentInfoMat[i,]))
                if (segmentInfoMat[i+1,1]-segmentInfoMat[i,2]>minSegSize)
                    mysegments=c(mysegments,c(segmentInfoMat[i,2]+1,segmentInfoMat[i+1,1]-1,0,0,0))
            }
        }
        mysegments=c(mysegments,c(segmentInfoMat[i+1,]))  
        
        if (ncol(X)-segmentInfoMat[i+1,2]-1>minSegSize) 
            mysegments=c(mysegments,c(segmentInfoMat[i+1,2]+1,ncol(X),0,0,0))
        
        mysegments=matrix(data=mysegments,nrow=length(mysegments)/5,
                          ncol=5,byrow=TRUE,dimnames=NULL)
        
        mysegments[which(mysegments[,4]==0),4]=refInd;
        mysegments[which(mysegments[,5]==0),5]=maxShift;
        
        if (sum(mysegments[,3]!=0)==0){
            cat("\n No segments are set for alignment! Please set 
                at least 1 values in columnn 3 in segmentInfoMat matrix be 1 ")
            return(NULL)
        }
        
        Y=X;
        
        for (i in 1:nrow(mysegments))
            if (mysegments[i,3]!=0)
            {
                if (verbose)
                    cat("\n Doing alignment a segment from ",
                        mysegments[i,1]," to ",mysegments[i,2]," ...");
                
                segmentpeakList=peakList;
                for (j in seq_along(peakList)){
                    segmentpeakList[[j]]=
                        findSegPeakList(peakList[[j]],mysegments[i,1],mysegments[i,2]);
                }    
                Y[,c(mysegments[i,1]:mysegments[i,2])]=
                    dohCluster(X[,c(mysegments[i,1]:mysegments[i,2])],peakList=segmentpeakList,
                               refInd=mysegments[i,4],maxShift =mysegments[i,5],
                               acceptLostPeak=acceptLostPeak, verbose=verbose);    
            }else{
                if (verbose)
                    cat("\n The segment ",
                        mysegments[i,1],"-",mysegments[i,2], " is not aligned");
            }  
        return(Y)
}