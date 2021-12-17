#' @title Minimum Convex Hull Polygon
#' @param xy Matrix or Data frame of occurrence coordinates
#' @param crs Character of coordinate reference system for minimum convex hull
#' @description Generates a minimum convex polygon (MCP; i.e., convex hull) that
#' is delineated from occurrence locality coordinates.
#' This function is a wrapper for `chull()` that makes a SpatialPolygons object.
#' @return a SpatialPolygons object of the minimum convex hull around occurrences.
#' @examples
#' # generate occurrences
#' ras1 <- raster::raster(nrows=108, ncols=108, xmn=-50, xmx=50)
#' raster::values(ras1)<- runif(n = (108*108))
#' occs <- dismo::randomPoints(ras1, 4)
#' # create mcp
#' mcp(occs)
#' @export

# make a minimum convex polygon as SpatialPolygons object
mcp <- function(xy, crs = NULL) {
  #require(sp)
  xy <- as.data.frame(sp::coordinates(xy))
  coords.t <- grDevices::chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  poly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(as.matrix(xy.bord))), 1)))
  if(!is.null(crs)) {
    poly@proj4string <- sp::CRS(crs)
  } else {
    message("WARNING: this minimum convex polygon has no coordinate reference system.")
  }
  return(poly)
}
