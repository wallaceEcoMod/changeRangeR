#' @title Calculate land cover overlapping with thresholded prediction in sq km
#' @description Calculate area of land cover of choice for thresholded prediction
#' @param r Raster layer of a binary SDM
#' @param proj Character string of r's projected coordinate system proj4string
#' @param locs data.frame of occurrence records: Longitude and latitude
#' @export

aooArea <- function(r, proj, locs=NULL) {
  ## If locs is NULL
  if (is.null(locs)){
  # Project everything to WGS84
    r.2km <- raster::projectRaster(from = r, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", method = "ngb")
  # Create new dummy raster with correct 2km resolution
    dummy <- r.2km
    res(dummy) <- 0.01666667
    # resample SDM
    r.2km <- raster::resample(x = r.2km, y = dummy)
    # Calculate number of cells
    fc.cells <- cellStats(!is.na(r.2km), stat = sum) * 2
    return(paste0("AOO:", fc.cells, "km^2"))
  }
  ## If interested in # of cells with occurrence points
  # Project everything to WGS84
  r.2km <- raster::projectRaster(from = r, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", method = "ngb")
  # Create new dummy raster with correct 2km resolution
  dummy <- r.2km
  res(dummy) <- 0.01666667
  # resample SDM
  r.2km <- raster::resample(x = r.2km, y = dummy)
  locCells <- extract(r.2km, locs)
  return(paste0("AOO of cells with occurrence records:", length(locCells) * 2, "km^2"))
}
