#' @title Calculate species endemism
#' @description Calculate species endemism as the number of species in a place divided by the total number of places in which those species are found.
#' @param rStack a rasterStack of binary species presences
#' @export
#' @author pgalante@@amnh.org
#' @return raster object showing species endemism.
#' @examples
#' # create binary raster
#' r1 <- raster::raster(nrows=108, ncols=108, xmn=-50, xmx=50)
#' raster::values(r1)<- runif(n = (108*108))
#' r1[r1 < 0.5] <- NA
#' r1[r1 > 0.5] <- 1
#' r2 <- raster::raster(nrows=108, ncols=108, xmn=-50, xmx=50)
#' raster::values(r2)<- runif(n = (108*108))
#' r2[r2 < 0.5] <- NA
#' r2[r2 > 0.5] <- 1
#' r3 <- raster::raster(nrows=108, ncols=108, xmn=-50, xmx=50)
#' raster::values(r3)<- runif(n = (108*108))
#' r3[r3 < 0.5] <- NA
#' r3[r3 > 0.5] <- 1
#' rStack <- stack(r1, r2, r3)
#' # calculate SE
#' SE(rStack)
#' @export


SE <- function (rStack){
  require(raster)
  rStack[rStack == 0] <- NA
  p.df <- as.data.frame(raster::rasterToPoints(rStack))
  sspTotal <- colSums(p.df, na.rm = T)
  ssp.df <- t(p.df[-c(1, 2)]) * (sspTotal[-c(1, 2)])
  ssp.PixSum <- rowSums(t(ssp.df), na.rm = T)
  spSum <- rowSums(p.df[-c(1, 2)], na.rm = T)
  SEvals <- spSum/ssp.PixSum
  stackTotal <- sum(rStack, na.rm = T)
  stackTotal[stackTotal == 0] <- NA
  stackTotal[!is.na(stackTotal)] <- SEvals
  return(stackTotal)
}
