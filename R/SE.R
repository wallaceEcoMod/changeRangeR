#' @title Calculate species endemism
#' @description Calculate species endemism as the number of species in a place divided by the total number of places in which those species are found.
#' @param rStack a rasterStack of binary species presences
#' @export
#' @author pgalante@@amnh.org
#'


SE <- function(rStack){
  # Convert rasterstack to points dataframe
  p.df <- as.data.frame(rasterToPoints(rStack))
  # For each species, sum total pixels
  sspTotal <- colSums(p.df, na.rm = T)
  # Multiply presence in each pixel by total for that species
  ssp.df <- t(p.df[-c(1,2)]) * (sspTotal[-c(1,2)])
  # For every pixel, calculate the sum of areas of species found there
  ssp.PixSum <- rowSums(t(ssp.df), na.rm = T)
  # for every pixel, calculate sum of species found
  spSum <- rowSums(p.df[-c(1,2)], na.rm=T)
  # Calculate SE by dividing
  SEvals <- spSum/ssp.PixSum
  # put back into raster
  stackTotal <- sum(rStack, na.rm = F)
  stackTotal[!is.na(stackTotal)] <- SEvals
  return(stackTotal)
}
