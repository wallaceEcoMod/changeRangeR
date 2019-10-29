#' @export

# make a minimum convex polygon as SpatialPolygons object
mcp <- function(xy, crs = NULL) {
  xy <- as.data.frame(sp::coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
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