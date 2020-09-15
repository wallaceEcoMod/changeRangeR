#' @title Calculate the proportion of range area that is contained by landcover categories
#' @description Calculate the proportion of the species' range (e.g., a thresholded SDM) that is contained by landcover categories
#' taken from a shapefile. Example shapefile categories include protected areas, threatened areas. ratioOverlap returns an s4 object
#' containing the masked raster layer and the percent of the total range that lies within the shapefile polygons specified.
#' @param r Binary Raster. Must be in same projection as shp parameter
#' @param shp a shapefile of land cover features.
#' @param field The shapefile field attribute containing the features to compare (i.e., the column name).
#' @param category a list of the names of shapefile features to compare. If all features are to be used, input "All".
#' @param proj character string proj4string of crs of landcover layer.
#' @export
#'
#'

### FIX RESOLUTION AND PROJECTION STUFF . Beth projected shapefile to the raster
## if one is projected in utm, project other to utm


ratioOverlap <- function(r, shp, field, category){
  #setClass("ratioOverlap", slots = list(maskedRange = "RasterLayer", ratio = "character"))
  require(sf)
  require(rgdal)
  require(raster)
  require(dplyr)

  if (category == "All"){
    shp <- sf::st_as_sf(shp)
    r <- raster::projectRaster(r, crs = crs(shp)@projargs, method = 'ngb')
    maskedRange <- raster::mask(r, shp)
  }
  else {
    shp <- sf::st_as_sf(shp)
    # if (crs == NULL){
    #   fc <- filter(shp, grepl(category, field))
    #   out<- mask(r, fc)
    # } else {
    r <- raster::projectRaster(r, crs = crs(shp))
    fc <- lapply(category, function(x) dplyr::filter(shp, shp[[field]]==x))
    fc <- do.call("rbind", fc)
    #fc <- filter(shp, shp[[field]]==category)
    maskedRange <- raster::mask(r, fc)
    #  }
    #  return(out)
  }
  ratio <- raster::ncell(maskedRange[!is.na(maskedRange)]) / raster::ncell(r[!is.na(r)]) * 100
  ratio <- paste0("Percentage of range within shape is ", ratio, "%")

  #out <- new("ratioOverlap", maskedRange = maskedRange, ratio = ratio)
  return(list(maskedRange = maskedRange, ratio = ratio))

}

