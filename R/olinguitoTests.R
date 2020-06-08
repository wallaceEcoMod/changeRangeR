# library(rgdal)
# library(raster)
# ### EOO occurrences
# locs <- read.csv("C:/Users/pgalante/OneDrive - AMNH/NASA-Wallace/changeRangeR/singleSpecies/Olinguito_data/Olinguito_RS_data/Occurrence_data/thinned_records/10KM_thin_2017.csv")
# head(locs)
# eoo <- mcp(locs[,1:2])
# area <- area(eoo)
# ### EOO SDM
# p <-raster('C:/Users/pgalante/Projects/NASA/maskRangeR/dataDriven/olinguito/olinguito_sdm.tif')
# thr <- 0.3380209
# p[p<thr] <- NA
# p.pts <- rasterToPoints(p)
# eooSDM <- mcp(p.pts[,1:2])
# aeoosdm <- area(eooSDM)
#
# ### AOO prethreshold
# # AOO if raster projected
# r.dummy <- r
# agg <- 2000 /res(r.dummy)[1]
# r.resam <- raster::aggregate(r, agg)
# plot(r)
#
#
# locs <- read.csv("C:/Users/pgalante/OneDrive - AMNH/NASA-Wallace/changeRangeR/singleSpecies/Olinguito_data/Olinguito_RS_data/Occurrence_data/thinned_records/10KM_thin_2017.csv")
# locs <-locs[,1:2]
# p <-raster('C:/Users/pgalante/Projects/NASA/maskRangeR/dataDriven/olinguito/olinguito_sdm.tif')
# # Sum 2x2km gridcells with occurrence points
# p[!is.na(p)] <- 1
# AOOlocs<-aooArea(r = p,proj = crs(p), locs = locs) # 36 km^2
#
# # Sum 2x2km gridcells from premasked thresholded SDM
# # convert to binary
# p[!is.na(p)] <- 1
# AOO<-aooArea(r = p,proj = crs(p)) # 87942 km^2
#
# # Sum 2x2km gridcells from masked thresholded SDM
# ## Function is same as above - Need to know which input raster to use
#
#
# ### Change in env values over time (forest cover) - and into future
# lc80 <- raster("C:/Users/pgalante/layers/cropsuitability/cropsuitability_rainfed_and_irrigated/1961-1990/overall_cropsuit_i_1961-1990/overall_cropsuit_i_1961-1990.bil")
# lc90 <- raster("C:/Users/pgalante/layers/cropsuitability/cropsuitability_rainfed_and_irrigated/1981-2010/overall_cropsuit_i_1981-2010/overall_cropsuit_i_1981-2010.bil")
# lc20 <- raster('C:/Users/pgalante/layers/cropsuitability/cropsuitability_rainfed_and_irrigated/2011-2040/overall_cropsuit_i_2011-2040/overall_cropsuit_i_2011-2040.bil')
# lc280 <- raster("C:/Users/pgalante/layers/cropsuitability/cropsuitability_rainfed_and_irrigated/2071-2100/overall_cropsuit_i_2071-2100/overall_cropsuit_i_2071-2100.bil")
# p <-raster('C:/Users/pgalante/Projects/NASA/maskRangeR/dataDriven/olinguito/olinguito_sdm.tif')
# lc80 <- crop(lc80, extent(p))
# lc90 <- crop(lc90, extent(p))
# lc20 <- crop(lc20, extent(p))
# lc280 <- crop(lc280, extent(p))
#
#
#
# ### Optimized model threshold
# p <-raster('C:/Users/pgalante/Projects/NASA/maskRangeR/dataDriven/olinguito/olinguito_sdm.tif')
# xy <- read.csv("C:/Users/pgalante/OneDrive - AMNH/NASA-Wallace/changeRangeR/singleSpecies/Olinguito_data/Olinguito_RS_data/Occurrence_data/thinned_records/10KM_thin_2017.csv")
# ch.orig <- mcp(locs[,1:2])
# thr <- 0.3380209
# SDMeoo <- mcpSDM(p = p, xy = xy[,1:2], ch.orig = ch.orig, thr = thr)
#
#
#
# ### ratio overlap test
# r <- raster('inst/extdata/DemoData/SDM/Forest_suitable_projected1.tif')
# # shp <- readOGR("C:/Users/pgalante/OneDrive - AMNH/NASA-Wallace/changeRangeR/singleSpecies/Olinguito_data/Olinguito_RS_data/WDPA_protected_areas/WDPA_May2020_COL-shapefile", "WDPA_May2020_COL-shapefile-polygons")
# shp <- readOGR("C:/Users/pgalante/OneDrive - AMNH/NASA-Wallace/changeRangeR/singleSpecies/Olinguito_data/Olinguito_RS_data/WDPA_protected_areas/WDPA_May2020_COL-shapefile", "WDPA_May2020_COL-shapefile-polygons")
# field <- "DESIG_ENG"
# category <- "All"
# test <- ratioOverlap(r = r, shp = shp, field = field, category = category) # 8.52774000069742%
