## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(tidy = FALSE)
figwidth.out <- 600
dpi.HQ <- 150
dpi.LQ <- 120

## ----wine data, dpi=dpi.HQ, fig.width=7, fig.height=4, out.width = figwidth.out----
library(speaq)
data(Winedata)
Spectra.wine <- as.matrix(Winedata$spectra )
ppm.wine <- as.numeric(Winedata$ppm) 
wine.color <- Winedata$wine.color 
wine.origin <- Winedata$origin 
# all spectra
speaq::drawSpecPPM(Y.spec = Spectra.wine, 
                   X.ppm = ppm.wine, 
                   title = 'Wine data spectra', 
                   groupFactor = wine.color, 
                   legend.extra.x = 1, 
                   legend.extra.y = 1.1)

## ----wine excerpt, dpi=dpi.HQ, fig.width=7, fig.height=4, out.width = figwidth.out----
# small excerpt by defining the region of interest
speaq::drawSpecPPM(Y.spec = Spectra.wine, 
                   X.ppm = ppm.wine, 
                   groupFactor = as.factor(wine.color), 
                   title = 'Raw wine data excerpt', 
                   legend.extra.x = 1.1, 
                   legend.extra.y = 1.0,
                   ROI.ppm = 3.6, 
                   ROI = NULL, 
                   roiWidth.ppm = 0.15,
                   legendpos = "topright" )

## ----detect winepeaks,  results = "hide"---------------------------------
wine.peaks <- speaq::getWaveletPeaks(Y.spec=Spectra.wine, 
                                     X.ppm=ppm.wine, 
                                     baselineThresh = 10,
                                     SNR.Th = -1, 
                                     nCPU = 2, 
                                     include_nearbyPeaks = TRUE) # nCPU set to 2 for the vignette build


wine.grouped <- speaq::PeakGrouper(Y.peaks = wine.peaks,  
                                   min.samp.grp = 5, 
                                   grouping.window.width = 200)


## ----plots base, dpi=dpi.HQ, fig.width=7, fig.height=10, fig.keep = "last", out.width = figwidth.out, warnings = FALSE----
# adding labels to the dat a for plotting and the group ppm values
library(ggplot2)
ROI.ppm <- 1.330
roiWidth.ppm <- 0.025

speaq::ROIplot(Y.spec = Spectra.wine, 
               X.ppm = ppm.wine, 
               ungrouped.peaks = wine.peaks,
               grouped.peaks = wine.grouped, 
               ROI.ppm = ROI.ppm,
               roiWidth.ppm = roiWidth.ppm, 
               groupLabels = as.factor(wine.color))

## ----silhouette values,  results = "hide", dpi=dpi.LQ, fig.width=6, fig.height=3.5, out.width = 500----
SilhouetteValues <- speaq::SilhouetR(DataMatrix = wine.grouped$peakPPM, 
                                     GroupIndices = wine.grouped$peakIndex)

Silh_plot <- ggplot(SilhouetteValues, aes(SilhouetteValues)) +
             geom_freqpoly(binwidth = 0.03) +
             theme_bw()
Silh_plot



## ----average silhouette, tidy = TRUE-------------------------------------
groups <- unique(SilhouetteValues$GroupIndices)
Ngroups <- length(groups)
sil_means <- matrix(NA, ncol = 3, nrow = Ngroups)

for(k in 1:Ngroups){
    sil_means[k,1] = groups[k]
    sil_means[k,2] = mean(SilhouetteValues$SilhouetteValues[SilhouetteValues$GroupIndices==groups[k]])
    sil_means[k,3] = mean(wine.grouped$peakSNR[wine.grouped$peakIndex==groups[k]])
}

sil_means <- sil_means[order(sil_means[,2]),]
colnames(sil_means) <- c("groupIndex", "avg_silhouette_val", "avg. SNR")
head(sil_means)


## ----wrong grouping plot,  dpi=dpi.HQ, fig.width=7, fig.height=10, fig.keep = "last", out.width = figwidth.out, warnings = FALSE----

faulty.groupIndex <- sil_means[1,1]
ROI.ppm <- ppm.wine[faulty.groupIndex]
roiWidth.ppm <- 0.1

speaq::ROIplot(Y.spec = Spectra.wine, 
               X.ppm = ppm.wine, 
               ungrouped.peaks = wine.peaks,
               grouped.peaks = wine.grouped, 
               ROI.ppm = ROI.ppm,
               roiWidth.ppm = roiWidth.ppm, 
               groupLabels = as.factor(wine.color))


## ----regroup-------------------------------------------------------------
wrong.groups <- sort(sil_means[sil_means[,1]>=sil_means[1,1],1])[1:2]

wine.regrouped <- speaq::regroupR(grouped.peaks = wine.grouped,
                                  list.to.regroup = wrong.groups, 
                                  min.samp.grp = 5,
                                  max.dupli.prop = 0.1)


## ----regroup fix plot,  dpi=dpi.HQ, fig.width=7, fig.height=10, fig.keep = "last", out.width = figwidth.out, warnings = FALSE----


faulty.groupIndex <- sil_means[1,1]
ROI.ppm <- ppm.wine[faulty.groupIndex]
roiWidth.ppm <- 0.1



speaq::ROIplot(Y.spec = Spectra.wine, 
               X.ppm = ppm.wine, 
               ungrouped.peaks = wine.peaks,
               grouped.peaks = wine.regrouped, 
               ROI.ppm = ROI.ppm,
               roiWidth.ppm = roiWidth.ppm, 
               groupLabels = as.factor(wine.color))



## ----data matrix, results = "hide", message=FALSE------------------------

wine.filled <- speaq::PeakFilling(Y.grouped = wine.regrouped, 
                                  Y.spec = Spectra.wine,  
                                  max.index.shift = 200,
                                  nCPU = 2) # nCPU set to 1 for the vignette build

wine.Features <- speaq::BuildFeatureMatrix(wine.filled)


## ----scaling-------------------------------------------------------------
wine.Features.scaled <- speaq::SCANT(data.matrix = wine.Features, 
                                     type = c("pareto", "center"))  


## ----PCA, dpi=dpi.LQ, fig.width=7, fig.height=5, out.width = 500---------


common.pca <- prcomp(wine.Features.scaled) 


loadings <- common.pca$rotation
scores <- common.pca$x
varExplained <- common.pca$sdev^2

barplot(varExplained/sum(varExplained), 
        main="Scree Plot",ylab="Proportion of variance explained", 
        xlab = "Principal comonent", 
        names.arg = as.character(seq(1,length(varExplained))))

## ----PCA2, dpi=200, fig.width=7, fig.height=5,out.width = figwidth.out, tidy = FALSE----
plot.marks <- as.numeric(wine.color)
plot.marks[plot.marks == 1] <- 8 
plot.marks[plot.marks == 2] <- 15
plot.marks[plot.marks == 3] <- 1

cp1 <- 1
cp2 <- 2 
plot(scores[,cp1]/max(scores[,cp1]), scores[,cp2]/max(scores[,cp2]),
     main=paste("score plot, PC",cp1," vs. PC",cp2,sep=""),
     xlab=paste("PC",cp1,round(varExplained[cp1]/sum(varExplained),digits=2),""),
     ylab=paste("PC",cp2,round(varExplained[cp2]/sum(varExplained),digits=2),""),
     pch = plot.marks)
text(scores[,cp1]/max(scores[,cp1]),scores[,cp2]/max(scores[,cp2]), wine.color, cex=0.5, pos=4, col="red")
lines(x = c(-100,100), y = c(0,0))
lines(x = c(0,0), y = c(-100,100))
legend("topleft", 
       legend = c("red  ","rosé  ","white      "), 
       pch = c(8,15,1), 
       y.intersp = 1.9)




## ----no rose-------------------------------------------------------------
red.white.scaled <- speaq::SCANT(wine.Features[wine.color!="rose",], 
                                 type = c("pareto", "center"))  

red.white.colors <- as.factor(as.character(wine.color)[wine.color!="rose"])


## ----relevant peaks------------------------------------------------------
p.all_bonf <- speaq::relevant.features.p(datamatrix = red.white.scaled,
                                         responsevector = red.white.colors, 
                                         p.adj = "bonferroni")

significant.features <- p.all_bonf[p.all_bonf$p.values<=0.05, ]

# order from most significant
significant.features <- significant.features[order(significant.features[,2]),]
head(significant.features)

## ----r significant features plots 1 , dpi=dpi.HQ, fig.width=5, fig.height=6, fig.keep = "all", tidy = FALSE, warnings = FALSE, out.width = 500----





peak_of_interest <- 5 # change this to the peak you want

interest.groupIndex <- significant.features$index[peak_of_interest]
interest.peakIndex <- as.numeric(rownames(significant.features))[peak_of_interest]


ROI.ppm <- ppm.wine[interest.groupIndex]
roiWidth.ppm <- 0.03


ggplot(p.all_bonf, aes(x=as.numeric(rownames(p.all_bonf)), y= -log10(p.values) )) + 
       geom_point(data = p.all_bonf[-interest.peakIndex,],  
                  aes(x=as.numeric(rownames(p.all_bonf[-interest.peakIndex,])), y= -log10(p.values) ),
                  shape = 16) +
       geom_point(data = p.all_bonf[interest.peakIndex,],  aes(x=interest.peakIndex, y= -log10(p.values) ),
                  shape = 18, 
                  size = 3, 
                  colour ="#00B0F6" ) +
       xlab("feature index") + 
       ylab("- log10 p-value") + 
       ggtitle("Bonferroni corrected p-values") +
       geom_hline(aes(yintercept= -log10(0.05), color="red"),linetype = 2) + 
       guides(color=FALSE) +
       theme_bw() + 
       theme(plot.title = element_text(lineheight = 0.8, face="bold", margin = margin(12,0,13,0),hjust = 0.5, size = 15), 
             text = element_text(size=14))






## ----r significant features plots 2 , dpi=dpi.HQ, fig.width=6, fig.height=8, fig.keep = "all", tidy = FALSE, warnings = FALSE, out.width = figwidth.out----




speaq::ROIplot(Y.spec = Spectra.wine,
               X.ppm = ppm.wine, 
               ungrouped.peaks = wine.peaks,
               grouped.peaks = wine.regrouped, 
               ROI.ppm = ROI.ppm,
               roiWidth.ppm = roiWidth.ppm, 
               groupLabels = as.factor(wine.color))



## ----plot significant features, dpi=dpi.HQ, fig.width=7, fig.height=4, out.width = figwidth.out----

peak_of_interest <- 5# change this number to any of the peaks you want to see
drawSpecPPM(Y.spec = Spectra.wine[wine.color != "rose", ], 
            X.ppm = ppm.wine, 
            groupFactor = red.white.colors, 
            title = paste("significant feature, p-value =",
                          format(significant.features$p.values[peak_of_interest], 
                                 scientific = TRUE, 
                                 digits=2),
                          sep=" "), 
            legend.extra.x = 1.1, 
            legend.extra.y = 0.9, 
            ROI = significant.features$index[peak_of_interest], 
            roiWidth = 100, 
            legendpos = "topright" )


## ----speaq 1.0,  dpi=dpi.HQ, fig.width=7, fig.height=4, message=F, out.width = figwidth.out----


peakList <- speaq::detectSpecPeaks(as.matrix(Spectra.wine),   
                                   nDivRange = 128,   
                                   scales = seq(1, 16, 2),
                                   baselineThresh = 100,
                                   SNR.Th = -1,
                                   verbose=FALSE)

resFindRef <- findRef(peakList)
refInd <- resFindRef$refInd

aligned.spectra <- speaq::dohCluster(as.matrix(Spectra.wine), 
                                     peakList = peakList,
                                     refInd = refInd,
                                     maxShift  = 200,
                                     acceptLostPeak = TRUE, verbose=FALSE)
      
                
speaq::drawSpecPPM(Y.spec = aligned.spectra[wine.color != "rose", ], 
                   X.ppm = ppm.wine, 
                   groupFactor = red.white.colors, 
                   title = paste("significant feature, p-value =",
                                 format(significant.features$p.values[peak_of_interest], 
                                        scientific = TRUE, 
                                        digits=2),
                                 sep=" "), 
                   legend.extra.x = 1.1, 
                   legend.extra.y = 0.9,
                   ROI = significant.features$index[peak_of_interest], 
                   roiWidth = 100, 
                   legendpos = "topright" )               

