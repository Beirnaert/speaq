## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(tidy = TRUE)
figwidth.out <- 600

## ----Read_data_input,fig.keep='none', tidy=FALSE, message=F, warning=F----
library(speaq)
#Generate a simulated NMR data set for this experiment
res=makeSimulatedData();
X=res$data;
groupLabel=res$label;

## ----Unaligned_spectral_plots--------------------------------------------
drawSpec(X);

## ----Peak_detection------------------------------------------------------
cat("\n detect peaks....");
startTime <- proc.time();
peakList <- detectSpecPeaks(X,
    nDivRange = c(128),                
    scales = seq(1, 16, 2),
    baselineThresh = 50000,
    SNR.Th = -1,
    verbose=FALSE
);

endTime <- proc.time();
cat("Peak detection time:",(endTime[3]-startTime[3])/60," minutes");

## ----Reference_finding---------------------------------------------------

cat("\n Find the spectrum reference...")
resFindRef<- findRef(peakList);
refInd <- resFindRef$refInd;

#The ranks of spectra
for (i in 1:length(resFindRef$orderSpec))
{
    cat(paste(i, ":",resFindRef$orderSpec[i],sep=""), " ");
    if (i %% 10 == 0) cat("\n")
}
    
cat("\n The reference is: ", refInd);

## ----Spectral_alignment--------------------------------------------------
# Set maxShift
maxShift = 50;

Y <- dohCluster(X,
                peakList = peakList,
                refInd = refInd,
                maxShift  = maxShift,
                acceptLostPeak = TRUE, verbose=FALSE);


## ----Spectral_alignment_optimal_maxShift,fig.align='center'--------------
Y <- dohCluster(X,
                peakList = peakList,
                refInd = refInd,
                maxShift  = NULL,
                acceptLostPeak = TRUE, verbose=TRUE);


## ----table, echo=FALSE---------------------------------------------------
library(knitr)

nghiaTable = matrix(c(c(100, 200, 0, 0, 0),c(450, 680, 1, 0, 50)), nrow = 2, byrow = T)
colnames(nghiaTable) = c("begin" , "end" , "forAlign" , "ref" , "maxShift")

kable(nghiaTable)

## ----Spectral_segment_alignment------------------------------------------
segmentInfoMat=matrix(data=c(100,200,0,0,0,
                      450,680,1,0,50),nrow=2,ncol=5,byrow=TRUE
                      )
colnames(segmentInfoMat)=c("begin","end","forAlign","ref","maxShift")
segmentInfoMat

Yc <- dohClusterCustommedSegments(X,
                                 peakList = peakList,
                                 refInd = refInd,
                                 segmentInfoMat = segmentInfoMat,
                                 minSegSize = 128,
                                 verbose=FALSE)
                                 

## ----AlignedSpectral_plots-----------------------------------------------
drawSpec(Y);

## ----AlignedSpectral_plots_limited_height--------------------------------
drawSpec(Y,
        startP=450,
        endP=680,
        highBound = 5e+5,
        lowBound = -100);

## ----Aligned_spectral_plots_customized-----------------------------------
drawSpec(Yc);

## ----Quantitative_analysis-----------------------------------------------
N = 100;
alpha = 0.05;

# find the BW-statistic
BW = BWR(Y, groupLabel);

# create sampled H0 and export to file
H0 = createNullSampling(Y, groupLabel, N = N,verbose=FALSE)

#compute percentile of alpha
perc = double(ncol(Y));
alpha_corr = alpha/sum(returnLocalMaxima(Y[2,])$pkMax>50000);
for (i in 1 : length(perc)){    
    perc[i] = quantile(H0[,i],1-alpha_corr, type = 3);
}

## ----drawBW_1, dpi=200, fig.width=7, fig.height=4, out.width = figwidth.out----

drawBW(BW, perc,Y, groupLabel = groupLabel)


## ----drawBW_2, dpi=200, fig.width=7, fig.height=4, out.width = figwidth.out----

drawBW(BW, perc, Y ,startP=450, endP=680, groupLabel = groupLabel)


