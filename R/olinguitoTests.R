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
#
# ### Optimized model threshold
# p <-raster('C:/Users/pgalante/Projects/NASA/maskRangeR/dataDriven/olinguito/olinguito_sdm.tif')
# xy <- read.csv("C:/Users/pgalante/OneDrive - AMNH/NASA-Wallace/changeRangeR/singleSpecies/Olinguito_data/Olinguito_RS_data/Occurrence_data/thinned_records/10KM_thin_2017.csv")
# ch.orig <- mcp(locs[,1:2])
# thr <- 0.3380209
# SDMeoo <- mcpSDM(p = p, xy = xy[,1:2], ch.orig = ch.orig, thr = thr)
#
# ### ratio overlap test
# r <- raster('C:/Users/pgalante/Projects/NASA/maskRangeR/dataDriven/olinguito/olinguito_sdm.tif')
# # shp <- readOGR("C:/Users/pgalante/OneDrive - AMNH/NASA-Wallace/changeRangeR/singleSpecies/Olinguito_data/Olinguito_RS_data/WDPA_protected_areas/WDPA_May2020_COL-shapefile", "WDPA_May2020_COL-shapefile-polygons")
# shp <- readOGR("C:/Users/pgalante/OneDrive - AMNH/NASA-Wallace/changeRangeR/singleSpecies/Olinguito_data/Olinguito_RS_data/WDPA_protected_areas/WDPA_May2020_Ecuador-shapefile", "WDPA_May2020_Ecuador-polygons")
# field <- "DESIG_ENG"
# category <- "All"
# test <- ratioOverlap(r = r, shp = shp, field = field, category = category) # 8.52774000069742%
