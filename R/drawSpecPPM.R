#' Plot NMR spectra from a spectra data matrix
#'
#' This function plots NMR spectra (so with the largest ppm values on the left) with a number of plotting options
#'
#' @param Y.spec (required) The raw spectra in matrix format (1 sample per row) or numeric vector (in case of 1 spectrum)
#' @param X.ppm (required) The vector with the ppm values
#' @param LeftIndex The starting index of the ppm values for plotting. default = -1 indicates the first ppm (the largest) value is the start of the plot  
#' @param RightIndex The stopping index for plotting. default = -1 indicates the last ppm value (the smallest) is the end of the plot
#' @param groupFactor The groupFactors. If provided different colors will be used for each group.
#' @param useLog If set to 'TRUE' the spectra will be log10 transformed (default = FALSE).
#' @param maxHeight The maximal height of the plot (default = -1, this indicates no maximal value).
#' @param minHeight The minimal height of the plot (default = -1, this indicates no minimal value).
#' @param nAxisPos The number of equally spaced tickmarks.
#' @param xlab The label on the x axis.
#' @param ylab The label on the y axis.
#' @param title The title of the plot.
#' @param ticks Position tick manually by providing ppm values.
#' @param ROI If provided (with an index value, not a ppm value) only this region of interest will be plotted. (supply no ROI or ROI.ppm values, for the full spectrum, or specify only 1, either ROI or ROI.ppm).
#' @param ROI.ppm If provided (a ppm value, not an index value) only this region of interest will be plotted. (supply no ROI or ROI.ppm values, for the full spectrum, or specify only 1, either ROI or ROI.ppm).
#' @param roiWidth The width of the ROI (region of interest) plot in index points/measurement points. The plot will span from ROI/ROI.ppm - roiWidth to ROI/ROI.ppm + roiWidth. (only supply roiWidth or roiWidth.ppm if needed).
#' @param roiWidth.ppm The width of the ROI (region of interest) plot in ppm. The plot will span from ROI/ROI.ppm - roiWidth.ppm to ROI/ROI.ppm + roiWidth.ppm. (only supply roiWidth or roiWidth.ppm if needed).
#' @param legend.extra.x Increase (or decrease) the horizontal space in the legend, this is useful for exporting larger figures.
#' @param legend.extra.y Increase (or decrease) the vertical space in the legend, this is useful for exporting larger figures.
#' @param legendpos The position of the legend (standard R legend positioning, default = 'topleft').
#' @param colourstyle The colours used in the plot, either standard R or ggplot colours (default).
#' @param manual.colours Provide specific colours to be used in the plot.
#' @param lwd The linewidth.
#' @param noLegend If set to TRUE no legend will be plotted (default = FALSE).
#' 
#' 
#' @return an R plot
#'
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' data(Winedata)
#' Spectra = Winedata$spectra 
#' ppm.wine = Winedata$ppm
#' wine.color = Winedata$wine.color 
#' drawSpecPPM(Y.spec=Spectra, X.ppm=ppm.wine, groupFactor = wine.color, 
#' title = 'Raw wine data spectra')
#' 
#' 
#' @export
#' 
#' @importFrom stats median
#' @importFrom grDevices hcl
#' @importFrom graphics plot par axis lines legend
#' 
#' 
drawSpecPPM <- function(Y.spec, X.ppm, LeftIndex = -1, RightIndex = -1, groupFactor = NULL, useLog = FALSE, 
    maxHeight = -1, minHeight = -1, nAxisPos = 4, xlab = NULL, ylab = NULL, title = NULL, ticks = NULL, 
    ROI = NULL, ROI.ppm = NULL, roiWidth = 100, roiWidth.ppm = NULL, legend.extra.x = 2, legend.extra.y = 2, 
    legendpos = NULL, colourstyle = "ggplot", manual.colours = NULL, lwd = 1, noLegend = FALSE) {
    
    if ("matrix" %in% class(X.ppm)) {
        stop("X.ppm is a matrix, for the plotting function only a numeric vector is allowed.")
    } else if (!"numeric" %in% class(X.ppm)) {
        X.ppm <- as.numeric(as.character(X.ppm))
    }
    if ("numeric" %in% class(Y.spec)) {
        Y.spec <- matrix(Y.spec, nrow = 1, ncol = length(Y.spec))
    } else {
        Y.spec <- as.data.frame(Y.spec)
        colnames(Y.spec) <- c(1:ncol(Y.spec))
        Y.spec <- as.matrix(Y.spec)
    }
    
    if (any(is.na(X.ppm))) {
        warning("X.ppm contains NA's, trying to remove these.")
        Y.spec <- Y.spec[,!is.na(X.ppm)]
        X.ppm <- X.ppm[!is.na(X.ppm)]
    }
    if(any(is.na(Y.spec))){
        stop("Y.spec contains NA's. Don't know how to deal with these.")
    }
    
    groupFactor_name <- groupFactor
    
    if (X.ppm[1] < X.ppm[2]) {
        message("X.ppm is not inverted. Inverting X.ppm and Y.spec")
        Y.spec <- Y.spec[,order(X.ppm, decreasing = T)]
        X.ppm <- rev(X.ppm)
    }
    
    
    if (!is.null(roiWidth.ppm)) {
        step <- stats::median(abs(diff(X.ppm)))
        roiWidth <- round(roiWidth.ppm/step)
    }
    
    if (!is.null(ROI.ppm)) {
        ROI <- which(abs(X.ppm - ROI.ppm) == min(abs(X.ppm - ROI.ppm)))[1]
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
    
    if (maxHeight != -1) {
        for (i in seq_len(nrow(Y.spec))) {
            myIndex <- which(Y.spec[i, ] > maxHeight)
            Y.spec[i, myIndex] <- maxHeight
        }
    }
    if (minHeight != -1) {
        for (i in seq_len(nrow(Y.spec))) {
            myIndex <- which(Y.spec[i, ] < minHeight)
            Y.spec[i, myIndex] <- minHeight
        }
    }
    if (is.null(groupFactor)) {
        groupFactor <- c(1:nrow(Y.spec))
        groupFactor <- as.factor(groupFactor)
    } else {
        if (!"factor" %in% class(groupFactor)) {
            warning("groupFactor is not a factor, attempting conversion.")
            groupFactor <- tryCatch({
                groupFactor <- as.factor(groupFactor)
            }, error = function(err) {
                groupFactor <- as.factor(as.numeric(groupFactor))
            })
            groupFactor_name <- groupFactor
        }
        levels(groupFactor) <- c(1:length(levels(groupFactor)))
    }
    if (LeftIndex == -1 & RightIndex == -1 & !is.null(ROI)) {
        LeftIndex <- ROI - roiWidth
        RightIndex <- ROI + roiWidth
        if (LeftIndex < 1) 
            LeftIndex <- 1
        if (RightIndex > ncol(Y.spec)) 
            RightIndex <- ncol(Y.spec)
    }
    if (LeftIndex == -1) {
        LeftIndex <- 1
    }
    if (RightIndex == -1) {
        RightIndex <- ncol(Y.spec)
    }
    if (is.null(xlab)) {
        xlab <- "ppm"
    }
    if (is.null(ylab)) {
        ylab <- "intensity"
    }
    if (is.null(title)) {
        title <- ""
        graphics::par(mar = c(4.5, 4.1, 1.5, 2.1))
    } else {
        graphics::par(mar = c(4.5, 4.1, 3.1, 2.1))
    }
    GraphRange <- c(LeftIndex:RightIndex)
    yn <- Y.spec[, GraphRange, drop = F]
    X.ppm.sub <- X.ppm[GraphRange]
    if (useLog) {
        yn <- log10(yn)
    }
    graphics::plot(yn[1, ], ylim = c(min(yn), max(yn)), type = "n", ylab = ylab, xlab = xlab, main = title, 
        xaxt = "n")
    tempVal <- trunc(length(GraphRange)/nAxisPos)
    xPos <- c(1, c(1:nAxisPos) * tempVal)
    if (is.null(X.ppm)) {
        graphics::axis(1, at = xPos, labels = xPos + LeftIndex)
    } else if (!is.null(ticks)) {
        ticks <- sort(ticks, decreasing = TRUE)
        xPos.user <- NULL
        if (!"numeric" %in% class(ticks)) 
            ticks <- as.numeric(as.character(ticks))
        for (tk in seq_along(ticks)) {
            xPos.user[tk] <- which(abs(X.ppm - ticks[tk]) == min(abs(X.ppm - ticks[tk])))[1]
            if (min(abs(X.ppm - ticks[tk])) > min(abs(diff(ticks)))) 
                stop(paste("tickmark", as.character(ticks[tk]), "is not close to any value in the provided ppm vector"))
        }
        graphics::axis(1, at = xPos.user, labels = ticks)
        
    } else {
        # this will plot in the ppm space get the possible scales for the plot
        max.logscale <- round(log10((max(X.ppm.sub) - min(X.ppm.sub)) * 2))
        possible.ppm.ticks <- c(as.numeric(matrix(c(1, 2, 5), byrow = TRUE, ncol = 3, nrow = 1) %o% matrix(10^(seq(-6, 
            (max.logscale - 1), by = 1)), byrow = TRUE, ncol = 1)), 10^max.logscale)
        # get the spacing in ppm
        axis.ppm.diff <- abs(X.ppm.sub[xPos][1] - X.ppm.sub[xPos][2])
        ppm.tickwidth <- possible.ppm.ticks[abs(possible.ppm.ticks - as.numeric(axis.ppm.diff)) == min(abs(possible.ppm.ticks - 
            as.numeric(axis.ppm.diff)))]
        Ndigits <- ceiling(-(log10(ppm.tickwidth))) + 1
        # get the spacing in indexes
        centerindex <- round(length(X.ppm.sub)/2)  # use the center index to avoid porblems
        # ppm sequence we want
        left.tick.ppm <- floor((X.ppm.sub[xPos][1])/(ppm.tickwidth/2)) * (ppm.tickwidth/2)
        right.tick.ppm <- ceiling((X.ppm.sub[xPos][length(xPos)])/(ppm.tickwidth/2)) * (ppm.tickwidth/2)
        ppm.ticks <- seq(left.tick.ppm, right.tick.ppm, -ppm.tickwidth)
        ppm.mindiff <- min(abs(diff(ppm.ticks)))
        xPos.auto <- NULL
        for (tk in seq_along(ppm.ticks)) {
            xPos.auto[tk] <- which(abs(X.ppm.sub - ppm.ticks[tk]) == min(abs(X.ppm.sub - ppm.ticks[tk])))[1]
            if (min(abs(X.ppm.sub - ppm.ticks[tk])) > ppm.mindiff/5) {
                xPos.auto[tk] <- NA
                ppm.ticks[tk] <- NA
            }
        }
        xPos.auto <- xPos.auto[!is.na(xPos.auto)]
        ppm.ticks <- ppm.ticks[!is.na(ppm.ticks)]
        graphics::axis(1, at = xPos.auto, labels = formatC(ppm.ticks, digits = Ndigits, format = "f"))
    }
    
    if (!is.null(manual.colours)) {
        plotcols <- manual.colours
    } else if (colourstyle == "ggplot") {
        nLevels <- length(levels(groupFactor))
        angles.col <- seq(0, 360, length.out = nLevels + 1) + 15
        plotcols <- grDevices::hcl(h = angles.col, l = 65, c = 120)[1:nLevels]
    } else {
        plotcols <- as.integer(levels(groupFactor))
    }
    
    for (i in 1:length(levels(groupFactor))) {
        groupFactorIdx <- which(groupFactor == levels(groupFactor)[i])
        for (j in seq_along(groupFactorIdx)) {
            graphics::lines(yn[groupFactorIdx[j], ], col = plotcols[i], lwd = lwd)
        }
    }
    
    if (!is.null(groupFactor_name) & !noLegend) {
        if (is.null(legendpos)) 
            legendpos <- "topleft"
        graphics::legend(legendpos, levels(groupFactor_name), col = plotcols, text.col = "black", lwd = lwd + 
            1, y.intersp = legend.extra.y, x.intersp = legend.extra.x)
    }
    
}
