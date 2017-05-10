#' BW ratio calculation
#'
#' Compute the BW ratios from data groups
#'
#'
#' @param X The spectral dataset in the matrix format in which each row contains a single sample.
#' @param groupLabel Group label of samples in the dataset.
#' 
#' @return Return BW ratio
#'
#' @author Trung Nghia Vu
#'
#' @seealso \code{\link{createNullSampling}}
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
#' 
#' @export
#' 
#'
BWR <-function(X, groupLabel)
{
    groupNum=length(levels(groupLabel));
    m=colMeans(X);
    groupMean=list();
    groupW=list();
    B=X;
    B[,]=0;
    for (i in 1:groupNum){
        groupLabeli=which(groupLabel==levels(groupLabel)[i]);
        Xi=X[groupLabeli,]
        mi=colMeans(Xi);
        groupMean[[i]]=mi;
        tempi=matrix(data=rep(mi,nrow(Xi)),nrow=nrow(Xi),
                     ncol=length(mi),byrow=TRUE);
        
        Wi=(Xi-tempi)^2; 
        groupW[[i]]=Wi;
        temp=matrix(data=rep(m,nrow(Xi)),nrow=nrow(Xi),
                    ncol=length(m),byrow=TRUE);
        B[groupLabeli,]=(tempi-temp)^2;    
    }  
    
    BW=double(ncol(B));
    for (i in 1:length(BW)){
        BW_denominator=0;
        for (j in 1:groupNum){
            BW_denominator=BW_denominator+sum(groupW[[j]][,i]);
        }
        Bw_numerator=sum(B[,i]);
        BW[i]=Bw_numerator/BW_denominator;    
    }
    return(BW); 
}