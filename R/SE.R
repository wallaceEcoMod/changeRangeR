#' @title Calculate species endemism
#' @description Calculate species endemism as the number of species in a place divided by the total number of places in which those species are found.
#' @param rStack a rasterStack of binary species presences
#' @export
#' @author pgalante@@amnh.org
#'


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
