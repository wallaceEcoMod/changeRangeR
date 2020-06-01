#' @title Calculate AOO
#' @description Calculate area of occupancy measured in 2km resolution using a binary SDM
#' @param r Raster layer of a binary SDM
#' @param proj Character string of r's projected coordinate system proj4string
#' @param locs (optional) data.frame of occurrence records: Longitude and latitude. If provided, AOO of cells containing occurrence points
#' is returned. If NULL, AOO of SDM is returned.
#' @export

#### FIX SO IT CAN TAKE ANY RESOLUTION or projection. Use native projection if available.


aooArea <- function(r, proj, locs=NULL) {
  ## If locs is NULL
  if (is.null(locs)){
    if (crs(r) == "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"){
      # Project everything to WGS84
      #  r.2km <- raster::projectRaster(from = r, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", method = "ngb")
      # Create new dummy raster with correct 2km resolution
      dummy <- r
      res(dummy) <- 0.01666667
      # resample SDM
      r.2km <- raster::resample(x = r, y = dummy, method = 'ngb')
      # Calculate number of cells
      fc.cells <- cellStats(!is.na(r.2km), stat = sum) * 4
      return(paste0("AOO:", fc.cells, "km^2"))
    } else {
      r.dummy <- r
      agg <- 2000 /res(r)[1]
      r.resam <- raster::aggregate(x = r, fact = agg, fun = 'max')
      #fc.cells<- cellStats(!is.na(r.resam), stat=sum)
      fc.cells <- cellStats(r.resam, na.rm=T,stat=sum) * 4
      return(paste0("AOO: ", fc.cells))
    } else {
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
  }
}
