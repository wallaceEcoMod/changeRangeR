### These are the to functions that Mary wanted us to see from the BioModelos github...
# I translated what the functions do, but the objects that go into the functions are not defined anywhere.

### Historic forest loss: Calculates the percentage of occurrence of the species in the forest for5 different years
# 6. `Historico perdida de bosque`: Calcula el porcentaje de ocurrencia de la especie en bosque para 5 años diferentes (1990, 2000, 2005, 2010, 2012).
ForestLoss <- function(r, area.raster, f90, f00, f05, f10, f12, rast=TRUE) {
  r <- r *area.raster
  com90 <- r * f90 # Los valores f90, f00, f05, f10, f12, corresponden al area en km2 de bosque en el pixel. El rango es de 0 a 1, siendo 1km2 el area total del pixel.
  com00 <- r * f00
  com05 <- r * f05
  com10 <- r * f10
  com12 <- r * f12

  forestArea <- data.frame(y1990 = cellStats(com90, sum),
                           y2000 = cellStats(com00, sum),
                           y2005 = cellStats(com05, sum),
                           y2010 = cellStats(com10, sum),
                           y2012 = cellStats(com12, sum))
  return(forestArea)
}

### Climate scenarios: Calculates the percentage of occurrence of the species in the forest in 3 climate scenarios: historic, conservationist (best case?), and industrialized (worst case I think?)
# 7. Escenarios 2030: Calcula el porcentaje de ocurrencia de la especie en bosque en trese escenarios para el 2030; histórico, conservacionista e industrializado.
Escenarios2030 <- function (r, area.raster, con2030_res,  des2030_res, his2030_res, sppFolder, rast) {
  inRaster <- r *area.raster
  escenarios2030 <- data.frame (con_2030 = NA, des_20302 = NA, his_2030 = NA)

  #con2030
  inRaster_2030c <- con2030_res * inRaster
  escenarios2030[1,1] <- cellStats (inRaster_2030c, stat = "sum")

  #des2030
  inRaster_2030d <- des2030_res * inRaster
  escenarios2030[1,2] <- cellStats (inRaster_2030d, stat = "sum")

  #his2030
  inRaster_2030h <- his2030_res * inRaster
  escenarios2030[1,3] <-cellStats (inRaster_2030h, stat = "sum")

  return (escenarios2030)
}
