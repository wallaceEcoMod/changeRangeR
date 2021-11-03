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


buildRMM <- function(binaryRange = NULL, rmm=NULL, locs = NULL, AOOarea=NULL, PE = NULL, PhyloTree = NULL,
                     complementarity = NULL, complementarity.of = NULL, complementarity.mask = NULL,
                     enChange = NULL, envChange.rStack = NULL, envChange.binaryRange = NULL, envChange.threshold=NULL,
                     envChange.bound = NULL, envChange.correlation = NULL, futureOverlap.binRasters = NULL,
                     futureOverlap.futures = NULL, mcp = NULL, mcpSDM = NULL, ratioOverlap = NULL, ratioOverlap.shape = NULL,
                     ratioOverlap.field = NULL, SE = NULL, SE.ranges = NULL){
  if(is.null(rmm)) {
    rmm <- rangeModelMetadata::rmmTemplate()
  }

  # AOOarea
  rmm$postprocess$inputs$binarySDM <- r
  rmm$postprocess$inputs$occurrences <- locs
  rmm$postprocess$AOO$Area <- AOOarea$area
  rmm$postprocess$AOO$raster <- AOOarea$aooRaster
  rmm$postprocess$AOO$withOccs <- AOOarea$aooPixels
  # phylogenetic endemism
  rmm$multispecies$PE <- PE
  rmm$multispecies$inputs$Phylogenetic_tree <- PhyloTree
  # complementarity
  rmm$postprocess$inputs$complementarity$type <- complementarity.of
  rmm$postprocess$inputs$complementarity$mask <- complementarity.mask
  rmm$postprocess$complementarity$withinMask <- complementarity$Percent_of_Total
  rmm$postprocess$complementarity$outsideMask <- complementarity$Percent_unique_values
  # envChange
  rmm$postprocess$inputs$TemporalSDM$rStack <- envChange.rStack
  rmm$postprocess$inputs$TemporalSDM$binaryRange <- envChange.binaryRange
  rmm$postprocess$inputs$TemporalSDM$threshold <- envChange.threshold
  rmm$postprocess$inputs$TemporalSDM$bound <- envChange.bound
  rmm$postprocess$inputs$TemporalSDM$correlation <- envChange.correlation
  rmm$postprocess$TemporalSDM$Area <- envChange$allAreas
  rmm$postprocess$TemporalSDM$masks <- envChange$masks
  # futureOverlap
  rmm$postprocess$inputs$TemporalOverlap$binaryRasters <- futureOverlap.binRasters
  rmm$postprocess$inputs$TemporalOverlap$futures <- futureOverlap.futures
  rmm$postprocess$TemporalOverlap$areaOverTime
  # mcp
  rmm$postprocess$mcp$mcpLocs <- mcp
  rmm$postprocess$mcp$mcpSDM <- mcpSDM
  # ratioOVerlap
  rmm$postprocess$inputs$binaryRange <- binaryRange
  rmm$postprocess$inputs$overlap.shape <- ratioOverlap.shape
  rmm$postprocess$inputs$overlap.field <- ratioOverlap.field
  rmm$postprocess$overlap$maskedRange <- ratioOverlap$maskedRange
  rmm$postprocess$overlap$rangeInCategory <- ratioOverlap$ratio
  rmm$postprocess$overlap$correlation <- ratioOverlap$correlation
  #Species endemism
  rmm$multispecies$inputs$ranges <- SE.ranges
  rmm$multispecies$SE <- SE
}

