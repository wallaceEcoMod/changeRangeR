#' @name complementarity
#' @title Compare continuous raster values within and outside of a mask
#' @param ras1 Raster object of categorical values
#' @param ras1mask Raster object of ras1 that has been previously masked
#' @export
#'

complementarity <- function(ras1, ras1mask){
  ## calculate percentage of values protected
  prop <- cellStats(ras1mask, stat = "sum", na.rm = T) / cellStats(ras1, stat = "sum", na.rm = T) * 100
  ## Get percentage of values that fall in polygons
  SRdf <- as.data.frame(table(as.matrix(ras1)))
  df <- as.data.frame(table(as.matrix(ras1mask)))
  # merge datafames
  dfMerged <- merge(SRdf, df, by = "Var1", all = T)
  dfMerged[is.na(dfMerged)] <- 0
  dfMerged$percent <- dfMerged$Freq.y / dfMerged$Freq.x * 100

  out <- list(Percent_of_Total = prop, Percent_per_Category = dfMerged)
  return(out)
}
