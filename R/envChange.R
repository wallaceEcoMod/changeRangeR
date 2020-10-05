#' @title Calculate change in suitable SDM area through time
#' @description Calculate SDM area after masking for environmental variables through time
#' @param SDM raster object of binary SDM with same projection as rStack
#' @param rStack rasterStack of environmental variable to measure within binary SDM through time
#' @param threshold integer of where rStack layers should be thresholded
#' @export


#SDM <- raster::raster("inst/extdata/DemoData/SDM/olinguito/Climatically_suitable_projected1.tif")
#rStack <- raster::stack(list.files(path = "inst/extdata/DemoData/MODIS", pattern = "\\.tif$", full.names = T))
#rStack <- raster::projectRaster(rStack, SDM, method = 'bilinear')
#threshold <- 50.086735

envChange <- function(rStack, SDM, threshold){
  require(raster)
  require(rgdal)

  rStack[rStack < threshold] <- NA
  rStack[rStack > threshold] <- 1
  masks <- lapply(raster::unstack(rStack), function(x) raster::mask(x, SDM))
  maskStack <- raster::stack(masks)

  if (!isLonLat(maskStack)){
    areas <- lapply(masks, function(x) res(x)[[1]] *ncell(x[!is.na(x)]))
  } else {
    area <- lapply(masks, raster::area)
    areas.1 <- lapply(area, function(x) x[!is.na(x)])
    areas <- lapply(areas.1, function(x) length(x) * median(x))
  }
  return(list(Area = areas, masks = maskStack))
}
