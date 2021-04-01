#' @title Calculate the ratio of future overlap of SDMs with shapefile categories
#' @description Calculate future overlap of SDMs with shapefile categories
#' @param r list of Rasters of binary SDMs
#' @param r.names list of character values of the names representing each raster in r
#' @param futures List of SpatialPolygons* objects with same CRS as r
#' @param futures.names list of character values of the names representing each SpatialPolygons* object in futures.
#' @param field The shapefile field attribute containing the features to compare (i.e., the column name).
#' @param category a list of the names of shapefile features to compare. If all features are to be used, input "All".
#' @export

futureOverlap <- function(r, futures, field, category, r.names, futures.names){
  # setClass("ratioOverlap", slots = list(maskedRange = "RasterLayer", ratio = "character"))
  require(sf)
  require(rgdal)
  require(raster)
  require(dplyr)
  require(sp)

  #rat.io <- lapply(futures, function(shp){
  if (category == "All"){
    shp <- lapply(futures, sf::st_as_sf)
    #r <- lapply(r, function(x) raster::projectRaster(x, crs = crs(shp[[1]]), method = "ngb"))
    # shp <- lapply(shp, function(x) sf::st_transform(x, crs = crs(r[[1]])))
    maskedRange <- mapply(raster::mask, r, shp)
    # maskedRange <- lapply(r, function(x) lapply(futures, function(y) raster::mask(x ,y)))
  } else {
    shp <- lapply(futures, sf::st_as_sf)
    # if (crs == NULL){
    #   fc <- filter(shp, grepl(category, field))
    #   out<- mask(r, fc)
    # } else {
    #r <- lapply(raster::projectRaster(r, crs = crs(shp)))
    #shp <- lapply(shp, function(x) sf::st_transform(x, crs = crs(r[[1]])))
    fc <- lapply(category, function(x) dplyr::filter(shp, shp[[field]]==x))
    fc <- do.call("rbind", fc)
    #fc <- filter(shp, shp[[field]]==category)
    maskedRange <- raster::mask(r, fc)
    #  }
    #  return(out)
  }
  #})

  # rat.io.s <- mapply(maskedRange, function(x){
  #   raster::ncell(maskedRange[!is.na(maskedRange)]) / raster::ncell(r[!is.na(r)]) * 100
  # })

  rat.io.s <- mapply(function(x, y){raster::ncell(x[!is.na(x)]) / raster::ncell(y[!is.na(y)]) * 100}, x=maskedRange, y=r)

  #ratio <- ncell(maskedRange[!is.na(maskedRange)]) / ncell(r[!is.na(r)]) * 100
  #ratio <- paste0("Percentage of range within shape is ", ratio, "%")
  ratioValues <- mapply(function(x,y,z){paste0("Overlap between ", x, " and ", y, " is ", z)}, x=r.names, y=futures.names, z=rat.io.s)
  ratioValues <- data.frame(strsplit(ratioValues, "is "))
  return(t(ratioValues))
  # out <- new("ratioOverlap", maskedRange = maskedRange, ratio = ratio)
}
