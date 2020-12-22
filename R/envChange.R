#' @title Calculate change in suitable SDM area through time
#' @description Calculate SDM area after masking for environmental variables through time
#' @param SDM raster object of binary SDM with same projection as rStack
#' @param rStack rasterStack of environmental variable to measure within binary SDM through time
#' @param threshold integer (or integers if bound = "both") of where rStack layers should be thresholded
#' @param bound character string characterizing the way the threshold should happen. "upper" removes values above
#' the threshold (e.g., maximum human footprint)."lower" removes values below the threshold (e.g., minimum forest cover).
#' "neither" does not threshold at any point. "both" thresholds at both threshold values (if provided; e.g.,
#' minimum and maximum temperature).
#' @param correlation boolean. If FALSE, environmental variable will be converted to a binary map and used as a mask.
#' If TRUE, environmental variable is only thresholded by bounds, but left continuous. Then#' Pearson's correlation
#' coefficient with SDM will be computed for overlapping areas.
#' @export


#SDM <- raster::raster("inst/extdata/DemoData/SDM/olinguito/Climatically_suitable_projected1.tif")
#rStack <- raster::stack(list.files(path = "inst/extdata/DemoData/MODIS", pattern = "\\.tif$", full.names = T))
#rStack <- raster::projectRaster(rStack, SDM, method = 'bilinear')
#threshold <- 50.086735

### NEED TO ADD IN CONTINUOUS OPTION OVER TIME PLOTTING

envChange <- function(rStack, SDM, threshold, bound, correlation = F){
  require(raster)
  require(rgdal)

  if(bound == "lower"){
#    if(correlation = FALSE){
      rStack[rStack < threshold] <- NA
      rStack[rStack > threshold] <- 1
#    } else {
#      rStack[rStack < threshold] <- NA
#    }
  }

  if(bound == "upper"){
#    if(correlation = FALSE){
      rStack[rStack < threshold] <- 1
      rStack[rStack > threshold] <- NA
#    } else {
#      rStack[rStack > threshold] <- NA
#    }
  }

  if(bound == "neither"){
    rStack = rStack
  }

  if(bound == "both"){
#    if(correlation = FALSE){
      rStack[rStack < min(thredhold)] <- NA
      rStack[rStack > max(threshold)] <- NA
      rStack[!is.na(rStack)] <- 1
#    } else {
#      rStack[rStack < min(thredhold)] <- NA
#      rStack[rStack > max(threshold)] <- NA
#    }
  }

  masks <- lapply(raster::unstack(rStack), function(x) raster::mask(x, SDM))
  maskStack <- raster::stack(masks)

  if (!isLonLat(maskStack)){
    areas <- lapply(masks, function(x) res(x)[[1]] *ncell(x[!is.na(x)]))
  } else {
    area <- lapply(masks, raster::area, na.rm = T)
    areas.1 <- lapply(area, function(x) x[!is.na(x)])
    areas <- lapply(areas.1, function(x) length(x) * median(x))
    #areas <- lapply(masks, raster::area)
  }
  allAreas <- cbind(unlist(areas))
  colnames(allAreas) <- "Area"
  rownames(allAreas) <- names(rStack)


  return(list(Area = allAreas, masks = maskStack))
}
