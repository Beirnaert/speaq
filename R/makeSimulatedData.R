#' Create a simulated NMR spectral data
#'
#' Generate an NMR spectral data for testing.
#'
#' @details 
#' We generate a NMR spectral data sets that contains two group A and group B. 
#' One at around 300 has a single tip and the other at around 600 has double 
#' tips that intentionally contains biological variation. First, a single spectrum 
#' is created based on statistic information (mean, standard deviation of intensity) 
#' achieved from real NMR spectra. Then, we randomly shift the spectrum to maximum 50 
#' data points and add some biological and technical variations to each point intensity 
#' to the spectrum to create a new spectrum. The collection of spectra from each group 
#' is the final dataset.
#'
#'
#' @return a list with 2 elements: data (The simulated NMR spectral data matrix) and label (Group label of each spectrum)
#'
#'
#' @author Trung Nghia Vu
#'
#' @examples
#' res <- makeSimulatedData();
#' X <- res$data;
#' groupLabel <- res$label;
#' 
#' @export
#' 
#' @importFrom stats rnorm median
#' 
makeSimulatedData <- function(){
        set.seed(2015);
        groupSize1=groupSize2=15;
        
        maxShift=50;
        pointNum=1000;
        pointNum=pointNum+maxShift*2;
        
        peakNum=2;
        peakCenterPos=c(300, 600);
        peakTipNum=c(1,2);
        
        x1=stats::rnorm(pointNum, mean = 9000, sd = 5000)
        
        tipMax=600000;
        slopeSpeed=0.01
        peakWide=14;
        
        peakPos=peakCenterPos[1];
        peakInd=1;
        x1[peakPos[peakInd]]=tipMax - sample(round(slopeSpeed*tipMax),1)
        for (i in 1:peakWide) {
            diffPerc=slopeSpeed*((1-slopeSpeed)+(sample(100,1)-50)/100)
            x1[peakPos[peakInd]+i]= x1[peakPos[peakInd]+i-1]*(1-slopeSpeed)* 
                x1[peakPos[peakInd]+i-1]/tipMax+ 
                sample(round(x1[peakPos[peakInd]+i-1]*diffPerc*2),1)-
                round(x1[peakPos[peakInd]+i-1]*diffPerc)
            
            diffPerc=slopeSpeed * ((1-slopeSpeed)+(sample(100,1)-50)/100)
            x1[peakPos[peakInd]-i]=x1[peakPos[peakInd]-i+1]*(1-slopeSpeed)*
                x1[peakPos[peakInd]-i+1]/tipMax+
                sample(round(x1[peakPos[peakInd]-i+1]*diffPerc*2),1)-
                round(x1[peakPos[peakInd]-i+1]*diffPerc)
        }
        
        
        
        tipMax=500000;
        slopeSpeed=0.05
        peakWide=10;
        
        peakPos=c(peakCenterPos[2]-peakWide,peakCenterPos[2]+peakWide);
        peakInd=1;
        x1[peakPos[peakInd]]=tipMax - sample(round(slopeSpeed*tipMax),1)
        for (i in 1:peakWide) {
            diffPerc=slopeSpeed*((1-slopeSpeed)+(sample(100,1)-50)/100)
            x1[peakPos[peakInd]+i]= x1[peakPos[peakInd]+i-1]*(1-slopeSpeed)*
                x1[peakPos[peakInd]+i-1]/tipMax+
                sample(round(x1[peakPos[peakInd]+i-1]*diffPerc*2),1)-
                round(x1[peakPos[peakInd]+i-1]*diffPerc)
            diffPerc=slopeSpeed * ((1-slopeSpeed)+(sample(100,1)-50)/100)
            x1[peakPos[peakInd]-i]=x1[peakPos[peakInd]-i+1]*(1-slopeSpeed)*
                x1[peakPos[peakInd]-i+1]/tipMax+
                sample(round(x1[peakPos[peakInd]-i+1]*diffPerc*2),1)-
                round(x1[peakPos[peakInd]-i+1]*diffPerc)
        }
        peakInd=2;
        x1[peakPos[peakInd]]=tipMax - sample(round(slopeSpeed*tipMax),1)
        for (i in 1:peakWide) {
            diffPerc=slopeSpeed*((1-slopeSpeed)+(sample(100,1)-50)/100)
            x1[peakPos[peakInd]+i]= x1[peakPos[peakInd]+i-1]*(1-slopeSpeed)*
                x1[peakPos[peakInd]+i-1]/tipMax+
                sample(round(x1[peakPos[peakInd]+i-1]*diffPerc*2),1)-
                round(x1[peakPos[peakInd]+i-1]*diffPerc)
            diffPerc=slopeSpeed * ((1-slopeSpeed)+(sample(100,1)-50)/100)
            x1[peakPos[peakInd]-i]=x1[peakPos[peakInd]-i+1]*(1-slopeSpeed)*
                x1[peakPos[peakInd]-i+1]/tipMax+
                sample(round(x1[peakPos[peakInd]-i+1]*diffPerc*2),1)-
                round(x1[peakPos[peakInd]-i+1]*diffPerc)
        }
        
        peakPos=peakCenterPos[2];
        x1[peakPos]= tipMax*0.25
        peakInd=1;
        for (i in 1:peakWide) {
            x1[peakPos[peakInd]+i]= x1[peakPos] + x1[peakPos[peakInd]+i]*0.75
            x1[peakPos[peakInd]-i]= x1[peakPos] + x1[peakPos[peakInd]-i+1]*0.75
        }
        
        
        
        groupLabel=as.factor(c(rep("Group 1",groupSize1),rep("Group 2",groupSize2)))
        X=NULL;
        for (i in seq_along(groupLabel)){
            groupType=groupLabel[i];
            
            x2=doShift(x1,shiftStep=sample(maxShift*2,1)-maxShift)
            
            xMed=stats::median(x1)
            xMed=round(xMed/1)
            x2=x2+ sample(xMed*2,1)-xMed
            
            peakRanges=c(550:650)
            xMax=max(x2[peakRanges]);
            if (groupType==levels(groupLabel)[1]){
                x2[peakRanges]=ifelse (x2[peakRanges]/xMax>0.5,x2[peakRanges]/xMax,0.5)*
                    x2[peakRanges]*((120+sample(50,1))/100)
            }else{
                x2[peakRanges]=ifelse (x2[peakRanges]/xMax>0.5,x2[peakRanges]/xMax,0.5)*
                    x2[peakRanges]*((25+sample(50,1))/100)
            }
            X=rbind(X,x2)
        }
        
        X=as.matrix(X[,(maxShift+1):(ncol(X)-maxShift)]);
        return(list(data=X,label=groupLabel))
    }