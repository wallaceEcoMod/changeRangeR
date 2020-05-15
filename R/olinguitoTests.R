# ### EOO occurrences
# locs <- read.csv("C:/Users/pgalante/OneDrive - AMNH/NASA-Wallace/changeRangeR/singleSpecies/Olinguito_data/Olinguito_RS_data/Occurrence_data/thinned_records/10KM_thin_2017.csv")
# head(locs)
# eoo <- mcp(locs[,1:2])
# area <- area(eoo)
# ### EOO SDM
# p <-raster('C:/Users/pgalante/Projects/NASA/maskRangeR/dataDriven/olinguito/olinguito_sdm.tif')
# xy <- read.csv("C:/Users/pgalante/OneDrive - AMNH/NASA-Wallace/changeRangeR/singleSpecies/Olinguito_data/Olinguito_RS_data/Occurrence_data/thinned_records/10KM_thin_2017.csv")
# ch.orig <- mcp(locs[,1:2])
# thr <- 0.3380209
# SDMeoo <- mcpSDM(p = p, xy = xy[,1:2], ch.orig = ch.orig, thr = thr)#### throws ERROR
# ### AOO prethreshold
#
# ### Change in env values over time (forest cover) - and into future
#
#
#
# ### ratio overlap test
# r <- raster('C:/Users/pgalante/Projects/NASA/maskRangeR/dataDriven/olinguito/olinguito_sdm.tif')
# # shp <- readOGR("C:/Users/pgalante/OneDrive - AMNH/NASA-Wallace/changeRangeR/singleSpecies/Olinguito_data/Olinguito_RS_data/WDPA_protected_areas/WDPA_May2020_COL-shapefile", "WDPA_May2020_COL-shapefile-polygons")
# shp <- readOGR("C:/Users/pgalante/OneDrive - AMNH/NASA-Wallace/changeRangeR/singleSpecies/Olinguito_data/Olinguito_RS_data/WDPA_protected_areas/WDPA_May2020_Ecuador-shapefile", "WDPA_May2020_Ecuador-polygons")
# field <- "DESIG_ENG"
# category <- "All"
# test <- ratioOverlap(r = r, shp = shp, field = field, category = category)
