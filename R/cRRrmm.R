#' @title Create metadata objects from changeRangeR
#' @description Creates and populates a \code{rangeModelMetadata} object from the output of \code{changeRangeR}.
#' See Merow \emph{et al.} (2019) for more details on the nature of the metadata and the \code{rangeModelMetadata} package.
#' To improve reproducibility of the study, this metadata object can be used as supplemental information for a manuscript, shared with collaborators, etc.
#' @param
#' @export
#' @author pgalante@@amnh.org
#' @return raster object showing species endemism.
#' @references Merow, C., Maitner, B. S., Owens, H. L., Kass, J. M., Enquist, B. J., Jetz, W., & Guralnick, R. (2019). Species' range model metadata standards: RMMS. \emph{Global Ecology and Biogeography}, \bold{28}: 1912-1924. \doi{10.1111/geb.12993}
#' @examples
#'


buildRMM <- function(rmm=NULL){
  if(is.null(rmm)) {
    rmm <- rangeModelMetadata::rmmTemplate()
  }

}

