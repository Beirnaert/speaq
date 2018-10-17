#' Building a null hypothesis data
#'
#' Create a null sampling data (N times) and write them to a file
#'
#'
#' @param X The spectral dataset in the matrix format in which each row contains a single sample
#' @param groupLabel Group label of samples in the dataset
#' @param N The number of iteration for creating null sample distribution
#' @param verbose A boolean value to allow print out process information
#'
#' @return A matrix with N rows containing the null distribution.
#'
#'
#' @author Trung Nghia Vu
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
#' H0 = createNullSampling(Y, groupLabel, N = 100,verbose=FALSE)
#' 
#' @export
#' 
#' 
createNullSampling <-function(X, groupLabel, N=100, 
                              verbose=TRUE){
    
    groupNum=length(levels(groupLabel));
    samplePool=X;  
    groupMean=list();
    
    for (i in 1:groupNum){
        groupLabeli=which(groupLabel==levels(groupLabel)[i]);
        Xi=X[groupLabeli,]
        mi=colMeans(Xi);
        groupMean[[i]]=mi;    
    }  
    
    for (i in seq_len(nrow(samplePool))){
        samplePool[i,]=
            X[i,]-groupMean[[which(levels(groupLabel)==groupLabel[i])]];
    }
    
    L=nrow(X);    
    H0=matrix(data=0,nrow=N,ncol=ncol(X));
    
    for(i in 1 : N){
        if (verbose) cat("\n Permutation th",i);
        index=sample(L);    
        H0[i,]=BWR(samplePool[index,],groupLabel);
    }
    return(H0);
}
