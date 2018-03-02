#' Local maximum detection
#'
#' Find and return local maximum of a single spectrum.
#' 
#' @param spectrum A spectral sample in the vector format.
#' 
#' 
#' @return list of 2: locMax (Locations of the found local maximum peaks) and pkMax (Intensities of the found local maximum peaks)
#' 
#' @author Trung Nghia Vu
#' 
#' @examples
#' res=makeSimulatedData();
#' X=res$data;
#' groupLabel=res$label;
#' returnLocalMaxima(X[2,])
#'                        
#' @export
#' 
returnLocalMaxima <-function(spectrum){
    l=length(spectrum);
    diffs=diff(spectrum);
    diffl=diffs[1:(l-1)];
    diffr=diffs[2:l];
    locMax=unique(c(which((diffl>=0&diffr<0))+1,which((diffl>0&diffr<=0))+1));
    pkMax=spectrum[locMax];
    return(list(pkMax=pkMax,locMax=locMax))
}
