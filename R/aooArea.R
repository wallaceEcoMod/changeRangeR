#' @title Calculate AOO
#' @description Calculate area of occupancy measured in 2km resolution using a binary SDM
#' @param r Raster layer of a binary SDM. Must be either unprojected in the WGS84 datum, or projected in a UTM projection measured in meters.
#' @param locs (optional) data.frame of occurrence records: Longitude and latitude. If provided, AOO of cells containing occurrence points
#' is returned. If NULL, AOO of SDM is returned.
#' @export


aooArea <- function (r, locs = NULL){
  if (is.null(locs)) {
    if (raster::isLonLat(r)) {
      dummy <- r
      raster::res(dummy) <- 0.01666667
      r.2km <- raster::resample(x = r, y = dummy, method = "ngb")
      fc.cells <- raster::cellStats(!is.na(r.2km), stat = sum) * 4
      rasArea <- paste0("AOO:", fc.cells, " km^2")
    }
    else {
      r.dummy <- r
      agg <- 2000/raster::res(r)[1]
      r.2km <- raster::aggregate(x = r, fact = agg, fun = "max")
      fc.cells <- raster::cellStats(r.2km, na.rm = T, stat = sum) *
        4
      rasArea <- paste0("AOO: ", fc.cells, " km^2")
    }
  }
  else {
    r.2km <- raster::projectRaster(from = r, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                                   method = "ngb")
    dummy <- r.2km
    raster::res(dummy) <- 0.01666667
    r.2km <- raster::resample(x = r.2km, y = dummy, method = "ngb")
    locCells <- raster::extract(r.2km, locs)
    rasArea <- paste0("AOO of cells with occurrence records:",
                      length(locCells) * 4, "km^2")
  }
  return(list(area = rasArea, aooRaster = r.2km))
}
