#' @title Calculate change in suitable SDM area through time
#' @description Calculate SDM area after masking for environmental variables through time
#' @param binaryRange raster object or shapefile of binary range (SDM, AOO, EOO) with same projection as rStack
#' @param rStack rasterStack of environmental variable to measure within binary SDM through time
#' @param threshold integer (or integers if bound = "both") of where rStack layers should be thresholded
#' @param bound character string characterizing the way the threshold should happen. "upper" removes values above
#' the threshold (e.g., maximum human footprint)."lower" removes values below the threshold (e.g., minimum forest cover).
#' "neither" does not threshold at any point. "both" thresholds at both threshold values (if provided; e.g.,
#' minimum and maximum temperature).
#' @param correlation boolean. If FALSE, environmental variable will be converted to a binary map and used as a mask.
#' If TRUE, environmental variable is only thresholded by bounds, but left continuous. Then, Pearson's correlation
#' coefficient with SDM will be computed for overlapping areas.
#' @return A list two objects. Area is a data.frame showing the total area of each masked raster. masks is a rasterStack
#' of each masked raster.
#' @examples
#' # create rStack
#' r1 <- raster::raster(nrows=108, ncols=108, xmn=-50, xmx=50)
#' raster::values(r1)<- runif(n = (108*108))
#' r2 <-  raster::raster(nrows=108, ncols=108, xmn=-50, xmx=50)
#' raster::values(r2)<- runif(n = (108*108))
#' rStack <-  raster::stack(r1,r2)
#' # create binaryRange
#' binaryRange <- raster::crop(r1, raster::extent(c(-50, 50, 0, 90)))
#' binaryRange <- raster::extend(binaryRange, r1)
#' binaryRange[!is.na(binaryRange)] <- 1
#' # set threshold
#' threshold <- 0.5
#' # set bound
#' bound <- "upper"
#' # Run function
#' envChange(rStack = rStack, binaryRange = binaryRange, threshold = threshold, bound = bound)
#' @author pgalante@@amnh.org
#' @export


#SDM <- raster::raster("inst/extdata/DemoData/SDM/olinguito/olinguitoSDM_coarse.tif")
#binaryRange <- raster::raster("inst/extdata/DemoData/SDM/olinguito/Climatically_suitable_projected_coarse.tif")
#binaryRange <- raster::projectRaster(binaryRange, SDM, method = "bilinear")
#rStack <- raster::stack(list.files(path = "inst/extdata/DemoData/MODIS", pattern = "\\.tif$", full.names = T))
#rStack <- raster::projectRaster(rStack, SDM, method = 'bilinear')
#threshold <- 50.086735
#test <- envChange(rStack, binaryRange, threshold, bound = "upper")
#test2 <- envChange(rStack, binaryRange, threshold, bound = "lower")

envChange <- function(rStack, binaryRange, threshold, bound, correlation = F){
  #require(raster)
  #require(rgdal)

  # if binaryRange is a shapefile, convert to raster then run like normal
  if(!("RasterLayer" %in% class(binaryRange))){
    binaryRange <- raster::rasterize(binaryRange, rStack)
  }

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
    rStack[rStack > threshold] <- NA
    rStack[rStack < threshold] <- 1
    #    } else {
    #      rStack[rStack > threshold] <- NA
    #    }
  }

  if(bound == "neither"){
    rStack = rStack
  }

  if(bound == "both"){
    #    if(correlation = FALSE){
    rStack[rStack < min(threshold)] <- NA
    rStack[rStack > max(threshold)] <- NA
    rStack[!is.na(rStack)] <- 1
    #    } else {
    #      rStack[rStack < min(threshold)] <- NA
    #      rStack[rStack > max(threshold)] <- NA
    #    }
  }

  masks <- lapply(raster::unstack(rStack), function(x) raster::mask(x, binaryRange))
  maskStack <- raster::stack(masks)

  if (!raster::isLonLat(maskStack)){
    areas <- lapply(masks, function(x) raster::res(x)[[1]] * raster::ncell(x[!is.na(x)]))
  } else {
    area <- lapply(masks, raster::area, na.rm = T)
    areas.1 <- lapply(area, function(x) x[!is.na(x)])
    areas <- lapply(areas.1, function(x) length(x) * stats::median(x))
    #areas <- lapply(masks, raster::area)
  }
  allAreas <- cbind(unlist(areas))
  colnames(allAreas) <- "Area"
  rownames(allAreas) <- names(rStack)


  return(list(Area = allAreas, masks = maskStack))
}
