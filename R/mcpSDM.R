#' @name mcpSDM
#' @title SDM-based Minimum Convex Polygon
#' @description Generates a minimum convex polygon (MCP; i.e., convex hull) that is delineated from a thresholded SDM. 
#' For each increment of 0.01 between a user-specified threshold and the maximum SDM prediction value, the prediction is
#' thresholded to this value to make a binary raster. This raster is then converted to points, which are used to delineate
#' a trial MCP. Each trial MCP is spatially intersected with the original MCP (based on the occurrence coordinates) and
#' the original occurrence points. The Jaccard similarity index is calculated to determine geographic similarity between
#' the trial and observed MCP. The trial MCP is also spatially intersected with the original occurrence points to determine
#' how many were omitted. The "best" MCP is the one that has the highest JSI and also omits the least original occurrence points. 
#' @export
#' 

mcpSDM <- function(p, xy, ch.orig, thr) {
  vals.p <- getValues(p)
  x <- seq(thr, max(vals.p, na.rm=TRUE), 0.01)
  jsi.vec <- numeric(length(x))
  ov.pts.vec <- numeric(length(x))
  ch.vec <- list()
  
  for(i in 1:length(x)) {
    th <- x[i]
    p.i <- p >= th
    p.i[p.i == 0] <- NA
    p.i.xy <- rasterToPoints(p.i)
    if(nrow(p.i.xy) > 1) {
      ch.i <- mcp(p.i.xy[,1:2], crs = wgs)
      ov.xy <- rgeos::gIntersection(ch.i, SpatialPoints(xy, proj4string = CRS(wgs)))
      if(!is.null(ov.xy)) ov.pts.vec[i] <- nrow(ov.xy@coords)
      ch.vec[[i]] <- ch.i
      ov <- rgeos::gIntersection(ch.i, ch.orig)
      if(!is.null(ov)) {
        A <- raster::area(ch.i)
        B <- raster::area(ch.orig)
        C <- raster::area(ov)
        jsi.vec[i] <- C / (A + B - C)  
      }  
    }
  }

  jsi.vec.allPts <- jsi.vec[which(ov.pts.vec == max(ov.pts.vec))]
  i.bestfit <- which(jsi.vec.allPts == max(jsi.vec.allPts))
  ch.bestfit <- ch.vec[[i.bestfit]]
  
  return(list(jsi = jsi.vec, thr = x, ov.pts = ov.pts.vec, best.fit = ch.bestfit, best.fit.ind = i.bestfit))
}
