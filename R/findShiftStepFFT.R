#' Finding the shift-step by using Fast Fourier Transform cross- correlation
#'
#' This function uses Fast Fourier Transform cross-correlation to find out the shift step between two spectra.
#' 
#' @param refSpec The reference spectrum.
#' @param tarSpec The target spectrum which needs to be aligned.
#' @param maxShift The maximum number of points for a shift step. If this value is zero, the algo- rithm will check on the whole length of the spectra.
#' @param scale Boolean value (TRUE/FALSE) for scaling data before Fast Fourier Transform cross-correlation step. If scale=NULL but mean/median of absolute data is too small (<1), the scaling will be applied. This might happen for very low abundant spectra like chromatograms. For normal NMR spectra, the scaling is usually not applied.
#' 
#' @return list of 2: corValue (The best correlation value) and stepAdj (The shift step found by the algorithm)
#' 
#' @author Trung Nghia Vu
#' 
#' @seealso \code{\link{hClustAlign}}
#' 
#' @examples
#' res=makeSimulatedData();
#' X=res$data;
#' groupLabel=res$label;
#' maxShift=50;
#' refSpec=X[1,];
#' tarSpec=X[2,];
#' adj=findShiftStepFFT(refSpec, tarSpec,maxShift=maxShift);
#'                        
#' @export
#' 
findShiftStepFFT<-function (refSpec, tarSpec, maxShift = 0, scale=NULL) 
{
    #do scaling if the spectra are low abundant
    refScale=stats::median(abs(refSpec)); if (refScale==0) refScale=mean(abs(refSpec))
    tarScale=stats::median(abs(tarSpec)); if (tarScale==0) tarScale=mean(abs(tarSpec))
    scaleFactor=10^round(-log10(min(tarScale,refScale)))
    if (is.null(scale) & scaleFactor > 1) scale=TRUE
    if (is.null(scale)) scale=FALSE
    if (scale){
        refScale=refScale*scaleFactor
        tarScale=tarScale*scaleFactor
    }
    
    
    M = length(refSpec)
    zeroAdd = 2^ceiling(log2(M)) - M
    r = c(refSpec*1e6, double(zeroAdd))
    s = c(tarSpec*1e6, double(zeroAdd))
    M = M + zeroAdd
    fftR = stats::fft(r)
    fftS = stats::fft(s)
    R = fftR * Conj(fftS)
    R = R/M
    rev = stats::fft(R, inverse = TRUE)/length(R)
    vals = Re(rev)
    maxi = -1
    maxpos = 1
    lenVals <- length(vals)
    if ((maxShift == 0) || (maxShift > M)) 
        maxShift = M
    if (anyNA(vals)) {
        lag = 0
        return(list(corValue = maxi, stepAdj = lag))
    }
    for (i in 1:maxShift) {
        if (vals[i] > maxi) {
            maxi = vals[i]
            maxpos = i
        }
        if (vals[lenVals - i + 1] > maxi) {
            maxi = vals[lenVals - i + 1]
            maxpos = lenVals - i + 1
        }
    }
    
    if (maxi < 0.1) {
        lag = 0
        return(list(corValue = maxi, stepAdj = lag))
    }
    if (maxpos > lenVals/2) 
        lag = maxpos - lenVals - 1
    else lag = maxpos - 1
    return(list(corValue = maxi, stepAdj = lag))
}
