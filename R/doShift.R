#' Segment shift
#'
#' Move a spectral segment of a sample shiftStep points to right or left
#'
#'
#' @param specSeg The segment which needs to be shifted
#' @param shiftStep The shift step for moving. If it is a negative (positive) value, the segment is moved to left (right).
#'
#' @return The new segment after shifting.
#'
#' @author Trung Nghia Vu
#' 
#' @seealso \code{\link{hClustAlign}}, \code{\link{findShiftStepFFT}}
#'
#' @examples
#' res=makeSimulatedData();
#' X=res$data;
#' groupLabel=res$label;
#' maxShift=50;
#' refSpec=X[1,];
#' tarSpec=X[2,];
#' adj=findShiftStepFFT(refSpec, tarSpec,maxShift=maxShift);
#' newTarSpec=doShift(tarSpec,adj$stepAdj);
#' 
#' @export
#' 
doShift <-function(specSeg, shiftStep){
    nFea=length(specSeg);
    newSegment=double(nFea);
    for (j in 1:nFea)
        if (shiftStep+j>0&&shiftStep+j<=nFea) 
            newSegment[shiftStep+j]=specSeg[j];
        if (shiftStep>0){
            for (j in 1: shiftStep) newSegment[j]=newSegment[shiftStep+1];
        }
        else{
            for (j in (nFea+shiftStep): nFea) 
                newSegment[j]=newSegment[(nFea+shiftStep-1)];
        }
        return (newSegment);
}
