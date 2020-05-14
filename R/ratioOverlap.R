#' @title Calculate the proportion of range area that is contained by landcover categories
#' @description Calculate the proportion of the species' range (e.g., a thresholded SDM) that is contained by landcover categories
#' taken from a shapefile. Example shapefile categories include protected areas, threatened areas.
#' @param r Binary Raster
#' @param shp a shapefile of land cover features
#' @param category a list of the names of shapefile features to compare.
#' @param proj proj4string of crs of landcover layer
#' @export
#'
#'

ratioOverlap <- function(r, shp, field, category){
  setClass("ratioOverlap", slots = list(maskedRange = "RasterLayer", ratio = "character"))
  require(sf)
  require(raster)
  require(dplyr)

  shp <- st_as_sf(shp)
  # if (crs == NULL){
  #   fc <- filter(shp, grepl(category, field))
  #   out<- mask(r, fc)
  # } else {
  r <- projectRaster(r, crs = crs(shp))
  fc <- lapply(category, function(x) filter(shp, shp[[field]]==x))
  fc <- do.call("rbind", fc)
  #fc <- filter(shp, shp[[field]]==category)
  maskedRange <- mask(r, fc)
  #  }
  #  return(out)
  ratio <- ncell(maskedRange[!is.na(maskedRange)]) / ncell(r[!is.na(r)]) * 100
  ratio <- paste0("Percentage of range within shape is ", ratio, "%")

  out <- new("ratioOverlap", maskedRange = maskedRange, ratio = ratio)
}



