#' @title Calculate the proportion of a range area that is either 1: contained by landcover categories, or 2: correlated with
#' a continuous environmental layer.
#' @description Calculate the proportion of the species' range (e.g., a thresholded SDM) that is contained by landcover categories
#' taken from a shapefile. Example shapefile categories include protected areas, threatened areas. ratioOverlap returns an s4 object
#' containing the masked raster layer and the percent of the total range that lies within the shapefile polygons specified.
#' @param r Binary Raster. Must be in same projection as shp parameter
#' @param shp (optional) a shapefile of land cover features.
#' @param rasMask (optional) a raster layer to calculate the relationship with the object r.
#' @param field The shapefile field attribute containing the features to compare (i.e., the column name).
#' @param category a list of the names of shapefile features to compare. If all features are to be used, input "All".
#' @param proj character string proj4string of crs of landcover layer.
#' @export
#'
#'



ratioOverlap <- function(r, shp = NULL, rasMask = NULL, field, category){
  #setClass("ratioOverlap", slots = list(maskedRange = "RasterLayer", ratio = "character"))
  require(sf)
  require(rgdal)
  require(raster)
  require(dplyr)

  if(!is.null(shp)){
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
  }

  if(!is.null(rasMask)){
    rasMask.resam <- raster::resample(rasMask, r, method = "bilinear")
    rRasmask <- raster::stack(rasMask.resam, r)
    layerCorrs <- raster::layerStats(rRasmask, stat = "pearson")
    correlation <- layerCorrs$`pearson correlation coefficient`[[2]]
  }

  #out <- new("ratioOverlap", maskedRange = maskedRange, ratio = ratio)
  return(list(maskedRange = maskedRange, ratio = ratio, correlation = correlation))

}

