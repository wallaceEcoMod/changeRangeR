#' @title Create metadata objects from changeRangeR
#' @description Creates and populates a \code{rangeModelMetadata} object from the output of \code{changeRangeR}.
#' See Merow \emph{et al.} (2019) for more details on the nature of the metadata and the \code{rangeModelMetadata} package.
#' To improve reproducibility of the study, this metadata object can be used as supplemental information for a manuscript, shared with collaborators, etc.
#' @param binaryRange Raster layer of a binary SDM. Must be either unprojected in the WGS84 datum, or projected in a UTM projection measured in meters.
#' @param rmm A `rangemodelmetadata` object.
#' @param locs `data.frame` of occurrence records: Longitude and latitude.
#' @param AOOarea The output of the function `AOOarea`.
#' @param PE The output of the function `calc_PE`.
#' @param PhyloTree class phylo object of phylogenetic tree. The names on tree's tip labels need to match the column names on the matrix.
#' @param complementarity The out put of the function `complementarity`.
#' @param complementarity.of `raster` of categorical values (e.g., a species richness map).
#' @param complementarity.mask `raster` object representing areas of interest (e.g., protected areas).
#' @param envChange The output of the function `envChange`
#' @param envChange.rStack `rasterStack` of environmental variable to measure within binary SDM through time.
#' @param envChange.binaryRange `raster` object or shapefile of binary range (SDM, AOO, EOO) with same projection as rStack.
#' @param envChange.threshold integer (or integers if bound = "both") of where rStack layers should be thresholded.
#' @param envChange.bound character string characterizing the way the threshold should happen. "upper" removes values above
#' the threshold (e.g., maximum human footprint)."lower" removes values below the threshold (e.g., minimum forest cover).
#' "neither" does not threshold at any point. "both" thresholds at both threshold values (if provided; e.g.,
#' minimum and maximum temperature).
#' @param envChange.correlation boolean. If FALSE, environmental variable will be converted to a binary map and used as a mask.
#' If TRUE, environmental variable is only thresholded by bounds, but left continuous. Then, Pearson's correlation
#' coefficient with SDM will be computed for overlapping areas.
#' @param futureOverlap The output of the function `futureOverlap`.
#' @param futureOverlap.binRasters `list` of rasters of binary SDMs..
#' @param futureOverlap.futures `list` of `SpatialPolygons*` objects with same CRS as futureOverlap.binRasters.
#' @param mcp The output of the function `mcp`
#' @param mcpSDM The output of the function `mcpSDM`
#' @param ratioOverlap The output of the function `ratioOverlap`
#' @param ratioOverlap.shape Either 1) a `shapefile` of land cover features or 2) a continuous `raster`. Must be in same projection as r parameter. If shp is a raster, then the number of cells within each quantile are calculated.
#' @param ratioOverlap.field A character string representing the shapefile field attribute containing the features to compare (i.e., the column name).
#' @param SE The output of the function `SE`
#' @param SE.ranges a `rasterStack` of binary species presences
#' @export
#' @author pgalante@@amnh.org
#' @return Populated rmm object
#' @references Merow, C., Maitner, B. S., Owens, H. L., Kass, J. M., Enquist, B. J., Jetz, W., & Guralnick, R. (2019). Species' range model metadata standards: RMMS. \emph{Global Ecology and Biogeography}, \bold{28}: 1912-1924. \doi{10.1111/geb.12993}
#' @examples
#' buildRMM(binaryRange = NULL)
#'


buildRMM <- function(binaryRange = NULL, rmm=NULL, locs = NULL, AOOarea=NULL, PE = NULL, PhyloTree = NULL,
                     complementarity = NULL, complementarity.of = NULL, complementarity.mask = NULL,
                     envChange = NULL, envChange.rStack = NULL, envChange.binaryRange = NULL, envChange.threshold=NULL,
                     envChange.bound = NULL, envChange.correlation = NULL, futureOverlap = NULL, futureOverlap.binRasters = NULL,
                     futureOverlap.futures = NULL, mcp = NULL, mcpSDM = NULL, ratioOverlap = NULL, ratioOverlap.shape = NULL,
                     ratioOverlap.field = NULL, SE = NULL, SE.ranges = NULL){

  #require(rangeModelMetadata)

  if(is.null(rmm)) {
    rmm <- rangeModelMetadata::rmmTemplate()
  }

  # AOOarea
  rmm$postprocess$inputs$binarySDM <- binaryRange
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
  rmm$postprocess$TemporalOverlap$areaOverTime <- futureOverlap
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

