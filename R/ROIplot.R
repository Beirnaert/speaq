#' Plot NMR spectra, together with raw and grouped peaks
#'
#' This function plots NMR spectra, peak plots and grouped peak plots all in figure for easy comparison.
#'
#' @param Y.spec (required) The raw spectra in matrix format (1 sample per row) or numeric vector (in case of 1 spectrum)
#' @param X.ppm (required) The vector with the ppm values
#' @param ungrouped.peaks (required) The data resulting from peak detecion with wavelets
#' @param grouped.peaks (required) The data after grouping (with PeakGrouper)
#' @param ROI If provided (with an index value, not a ppm value) only this region of interest will be plotted. (supply no ROI or ROI.ppm values, for the full spectrum, or specify only 1, either ROI or ROI.ppm).
#' @param ROI.ppm If provided (a ppm value, not an index value) only this region of interest will be plotted. (supply no ROI or ROI.ppm values, for the full spectrum, or specify only 1, either ROI or ROI.ppm).
#' @param roiWidth The width of the ROI (region of interest) plot in index points/measurement points. The plot will span from ROI/ROI.ppm - roiWidth to ROI/ROI.ppm + roiWidth. (only supply roiWidth or roiWidth.ppm if needed).
#' @param roiWidth.ppm The width of the ROI (region of interest) plot in ppm. The plot will span from ROI/ROI.ppm - roiWidth.ppm to ROI/ROI.ppm + roiWidth.ppm. (only supply roiWidth or roiWidth.ppm if needed).
#' @param groupLabels The vector with group labels (as factors)
#' @param output Whether to return a plot (default), or the individual ggplot objects (output = "ggObjects")
#' 
#' 
#' @return a plot
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#' 
#' @examples 
#' subset <- GetWinedata.subset()
#' subset.spectra = as.matrix(subset$Spectra)
#' subset.ppm = as.numeric(subset$PPM)
#'
#' test.peaks <- getWaveletPeaks(Y.spec=subset.spectra, 
#'                               X.ppm=subset.ppm,
#'                               nCPU = 1) # nCPU set to 2 for the vignette build
#'
#' test.grouped <- PeakGrouper(Y.peaks = test.peaks)
#'                            
#' ROI.ppm <- 4.9
#' roiWidth.ppm <- 0.15
#'
#' plots <- ROIplot(Y.spec = subset.spectra, 
#'                  X.ppm =subset.ppm, 
#'                  ungrouped.peaks = test.peaks,
#'                  grouped.peaks = test.grouped ,
#'                  ROI.ppm = ROI.ppm,
#'                  roiWidth.ppm = roiWidth.ppm , 
#'                  output = "ggObjects"
#'                  )
#'
#' @export
#' 
#' @import ggplot2
#' @importFrom stats median
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange
#' 
ROIplot <- function(Y.spec, X.ppm, ungrouped.peaks, grouped.peaks, ROI = NULL, ROI.ppm = NULL, roiWidth = 100, roiWidth.ppm = NULL, groupLabels = NULL, output = NULL) {
 
    if (!is.null(roiWidth.ppm)) {
        step <- stats::median(abs(diff(X.ppm)))
        roiWidth <- round(roiWidth.ppm/step)
    } else if (!is.null(roiWidth)){
        step <- stats::median(abs(diff(X.ppm)))
        roiWidth.ppm <- roiWidth*step
    }
    
    if(is.null(ROI.ppm) & is.null(ROI)){
        stop("no ROI or ROI.ppm is provided.")
    }
    
    if (!is.null(ROI.ppm)) {
        ROI <- which(abs(X.ppm - ROI.ppm) == min(abs(X.ppm - ROI.ppm)))[1]
    }
    
    if((ROI - roiWidth) < 1 | (ROI + roiWidth) > length(X.ppm) ){
        # this means roiwidth is chosen too large for the available ppm range
        roiWidth = roiWidth - abs(min(c((ROI - roiWidth), (length(X.ppm)-(ROI + roiWidth))))) - 1
        roiWidth.ppm = abs( (X.ppm[(ROI - roiWidth)] - X.ppm[(ROI + roiWidth)]) / 2 )
        print("plotting width chosen too large for available range (roiWidth or roiWidth.ppm). Rescaling.")
    }
    
    
    if (is.null(ROI)) {
        if (10 * min(abs(diff(X.ppm))) < max(abs(diff(X.ppm)))) {
            message(paste("There might be a gap in the data, the minimal and maximal differences between consecutive ppm points are respectively", 
                          as.character(min(abs(diff(X.ppm)))), "and", as.character(max(abs(diff(X.ppm))))))
        }
    } else {
        if (10 * min(abs(diff(X.ppm[(ROI - roiWidth):(ROI + roiWidth)]))) < max(abs(diff(X.ppm[(ROI - 
                                                                                                roiWidth):(ROI + roiWidth)])))) {
            message(paste("There might be a gap in the data, the minimal and maximal differences between consecutive ppm points are respectively", 
                          as.character(min(abs(diff(X.ppm)))), "and", as.character(max(abs(diff(X.ppm))))))
        }
    }
    
    if (is.null(groupLabels)) {
        groupLabels <- c(1:nrow(Y.spec))
        groupLabels <- as.factor(groupLabels)
    } else {
        if (!"factor" %in% class(groupLabels)) {
            warning("groupLabels is not a factor, attempting conversion.")
            groupLabels <- tryCatch({
                groupLabels <- as.factor(groupLabels)
            }, error = function(err) {
                groupLabels <- as.factor(as.numeric(groupLabels))
            })
        }
    }
    
    if(is.null(output)){
         
    } else if(output != "ggObjects"){
        output = NULL
    }
    
    peaks.plot <- AddPlottingStuff(Y.peaks = ungrouped.peaks, 
                                           X.ppm = X.ppm, 
                                           groupLabels = groupLabels)

    grouped.plot <- AddPlottingStuff(Y.peaks = grouped.peaks, 
                                             X.ppm = X.ppm, 
                                             groupLabels = groupLabels)
    
    peakValRange = range(ungrouped.peaks$peakValue, finite = TRUE)
    Range.extra = (peakValRange[2] - peakValRange[1]) * 0.01

    pp1 <- ggplot(peaks.plot[peaks.plot$peakPPM > (ROI.ppm - roiWidth.ppm ) & 
                             peaks.plot$peakPPM < (ROI.ppm + roiWidth.ppm ) ,], 
                  aes_string(x = "peakPPM", y = "peakValue", colour = "label") ) + 
           geom_point() + 
           theme_bw() + 
           xlim(c(ROI.ppm + roiWidth.ppm, ROI.ppm - roiWidth.ppm)) +
           ylim(c(peakValRange[1] - Range.extra, peakValRange[2]+Range.extra)) +
           labs(x = "ppm", y = "peak value") + 
           ggtitle("After peak detection") +
           theme(legend.title = element_blank(),
                  plot.title = element_text(face = "bold",hjust = 0.5)) 

    pp2 <- ggplot(grouped.plot[grouped.plot$peakPPM > (ROI.ppm - roiWidth.ppm ) &
                               grouped.plot$peakPPM < (ROI.ppm + roiWidth.ppm ) ,],
                  aes_string(x = "groupPPM", y = "peakValue", colour = "label") ) +
           geom_point() + 
           theme_bw() + 
           xlim(c(ROI.ppm + roiWidth.ppm, ROI.ppm - roiWidth.ppm)) +
           ylim(c(peakValRange[1] - Range.extra, peakValRange[2]+Range.extra)) +
           labs(x = "ppm", y = "peak value") + 
           ggtitle("After grouping") +
           theme(legend.title = element_blank(),
                 plot.title = element_text(face = "bold",hjust = 0.5)) 

    ROImatrix <- Y.spec[,X.ppm > (ROI.ppm - roiWidth.ppm ) & X.ppm < (ROI.ppm + roiWidth.ppm )]
    ROIppm <- X.ppm[X.ppm > (ROI.ppm - roiWidth.ppm ) & X.ppm < (ROI.ppm + roiWidth.ppm )]
    colnames(ROImatrix ) <- ROIppm
    ROI.df <- reshape2::melt(t(ROImatrix))
    names(ROI.df) <- c("ppm", "sample", "intensity")
    ROI.df$class <-  groupLabels[ROI.df$sample]

    pp0 <- ggplot(ROI.df, aes_string(x = "ppm", y = "intensity", group = "sample", colour = "class")) +
           scale_x_reverse()+
           theme_bw() +
           geom_line(size = 0.3) +
           ggtitle("Spectra") +
           theme(legend.title = element_blank(),
                 plot.title = element_text(face = "bold",hjust = 0.5)) 
  
    if(is.null(output)){
        gridExtra::grid.arrange(pp0, pp1, pp2, ncol=1)
    } else{
        return(list(SpectraPlot = pp0, peakPlot = pp1, GroupPlot = pp2))
    }
}

