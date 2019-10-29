#' @title Calculate land cover overlapping with thresholded prediction in sq km
#' @description Calculate area of land cover of choice for thresholded prediction
#' @param r Binary Raster object
#' @param proj Projected coordinate system as CRS string
#' @export

aooArea <- function(r, proj) {
	r.2km <- projectRaster(lc, crs = proj, res = 2000, method = "ngb")
	fc.cells <- cellStats(!is.na(r.2km), stat = sum)
	# multiply by area of cell (2 km * 2 km)
	return(4 * fc.cells)
}
