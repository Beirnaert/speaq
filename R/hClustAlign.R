#' CluPA function for two spectra.
#'
#' This function implements the idea of the CluPA algorithm to align the target spectrum against the reference spectrum.
#' 
#' @param refSpec The reference spectrum.
#' @param tarSpec The target spectrum.
#' @param peakList List of peaks of the both reference and target spectra
#' @param peakLabel The list of the labels of the peaks
#' @param startP The starting point of the segment.
#' @param endP The ending point of the segment.
#' @param distanceMethod The distance method for the hierarchial clustering algorithm.
#' @param maxShift The maximum number of points for a shift step.
#' @param acceptLostPeak This is an option for users, TRUE is the default value. If the users believe that all the peaks in the peak list are true positive, change it to FALSE.
#' 
#' @return list of 2: tarSpec (The target spectrum after alignment) and peakList (The peak list after alignment)
#' 
#' @author Trung Nghia Vu
#' 
#' @references Vu TN, Valkenborg D, Smets K, Verwaest KA, Dommisse R, Lemi\`{e}re F, Verschoren A, Goethals B, Laukens K. (2011) An integrated workflow for robust alignment and simplified quantitative analysis of NMR spectrometry data. BMC Bioinformatics. 2011 Oct 20;12:405.
#' 
#' @seealso \code{\link{dohCluster}} 
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
#' tarInd=1;
#' refSpec=X[refInd,];
#' tarSpec=X[tarInd,];
#' mergedPeakList=c(peakList[[refInd]],peakList[[tarInd]]);
#' mergedPeakLabel=double(length(mergedPeakList));
#' for (i in 1:length(peakList[[refInd]]) ) mergedPeakLabel[i]=1;
#' startP=1;
#' endP=length(tarSpec);
#' res=hClustAlign(refSpec,tarSpec,mergedPeakList,mergedPeakLabel,startP,endP,
#'                 maxShift=50,acceptLostPeak=TRUE)
#'                        
#' @export
#' 
hClustAlign <-function(refSpec, tarSpec, peakList, peakLabel, startP, endP,
                       distanceMethod="average", maxShift=0, acceptLostPeak=FALSE){
    minpeakList=min(peakList)[1];
    maxpeakList=max(peakList)[1];
    startCheckP=startP+which.min(tarSpec[startP:(minpeakList-1)])[1]-1;
    if (is.na(startCheckP)) {
        startCheckP=startP;
    }
    if(startCheckP < 1){
        startCheckP <- startP
    }
    
    endCheckP=maxpeakList+ which.min(tarSpec[(maxpeakList+1):endP])[1];
    if (is.na(endCheckP)){
        endCheckP=endP;
    } 
    if(endCheckP > length(tarSpec)){
        endCheckP <- endP
    }
    
    if ((endCheckP-startCheckP)<2) {
        return (list(tarSpec=tarSpec,peakList=peakList));
    }
    adj=findShiftStepFFT(refSpec[startCheckP:endCheckP],
                         tarSpec[startCheckP:endCheckP],maxShift=maxShift);
    if (adj$stepAdj!=0){
        if (acceptLostPeak) isImplementShift=TRUE 
        else isImplementShift=(adj$stepAdj<0&&adj$stepAdj+
                                   minpeakList >=startCheckP )||(adj$stepAdj>0&&adj$stepAdj+
                                                                     maxpeakList<=endCheckP);
        if (isImplementShift)
        {
            newTargetSpecRegion=doShift(tarSpec[startCheckP:endCheckP],adj$stepAdj);
            tarSpec[startCheckP:endCheckP]=newTargetSpecRegion;
            
            peakListTarget=which(peakLabel==0);
            peakList[peakListTarget]=peakList[peakListTarget]+adj$stepAdj;
            
            lostPeaks <- which(peakList <= 0 | peakList > length(tarSpec))
            if (length(lostPeaks) >0){
                
                peakList=peakList[-lostPeaks];
                peakLabel=peakLabel[-lostPeaks];
            }
        }
    }
    
    if (length(peakList)<3) {return (list(tarSpec=tarSpec,peakList=peakList));}
    hc=stats::hclust(stats::dist(peakList),method=distanceMethod)
    clusterLabel=stats::cutree(hc,h=hc$height[length(hc$height)-1]);
    if (length(unique(clusterLabel))<2){ 
        return (list(tarSpec=tarSpec,peakList=peakList));
    }
    
    labelID1=which(clusterLabel==1);
    subData1=peakList[labelID1];
    subLabel1=peakLabel[labelID1];
    
    labelID2=which(clusterLabel==2);
    subData2=peakList[labelID2];
    subLabel2=peakLabel[labelID2];
    maxsubData1=max(subData1)[1];
    minsubData2=min(subData2)[1];
    
    if (maxsubData1<minsubData2){
        endP1=maxsubData1+which.min(tarSpec[(maxsubData1+1) :(minsubData2-1)])[1];
        if (is.na(endP1)) endP1=maxsubData1;
        startP2=endP1+1;
        if (length(unique(subLabel1))>1){
            res=hClustAlign(refSpec,tarSpec,subData1,subLabel1,startP,endP1,
                            distanceMethod=distanceMethod,maxShift=maxShift,
                            acceptLostPeak=acceptLostPeak);
            tarSpec=res$tarSpec;
            peakList[labelID1]=res$peakList;
        }        
        if (length(unique(subLabel2))>1){
            res=hClustAlign(refSpec,tarSpec,subData2,subLabel2,startP2,endP,
                            distanceMethod=distanceMethod,maxShift=maxShift,
                            acceptLostPeak=acceptLostPeak);
            tarSpec=res$tarSpec;
            peakList[labelID2]=res$peakList;
        }
    }else{        
        maxsubData2=max(subData2)[1];
        minsubData1=min(subData1)[1];
        endP2=maxsubData2+which.min(tarSpec[(maxsubData2+1) :(minsubData1-1)])[1];
        if (is.na(endP2)) endP2=maxsubData2;
        startP1=endP2+1;
        if (length(unique(subLabel2))>1){
            res=hClustAlign(refSpec,tarSpec,subData2,subLabel2,startP,endP2,
                            distanceMethod=distanceMethod,maxShift=maxShift,
                            acceptLostPeak=acceptLostPeak);
            tarSpec=res$tarSpec;
            peakList[labelID2]=res$peakList;
        }    
        if (length(unique(subLabel1))>1){
            res=hClustAlign(refSpec,tarSpec,subData1,subLabel1,startP1,endP,
                            distanceMethod=distanceMethod,maxShift=maxShift,
                            acceptLostPeak=acceptLostPeak);
            tarSpec=res$tarSpec;
            peakList[labelID1]=res$peakList;
        }        
    }    
    return (list(tarSpec=tarSpec,peakList=peakList));
}