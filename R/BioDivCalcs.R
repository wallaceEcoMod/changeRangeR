# library(raster)
# library(ape)
# ## Using Nivel 2 binary maps (34 species) calculate: PE, PD using trees 1, 300, 500, 1000
# #### Among trees, get a measure of variation
# # load rasters and edit names
# primRas <- stack(list.files(path = "C:/Users/pgalante/Projects/NASA/changeRangeR/Nivel2", pattern = "\\.tif$", full.names = T))
# names(primRas) <- gsub("*_veg", "", names(primRas))
# # Drop models that are resolved into one group; Cebus_albifrons
# primRas <- dropLayer(primRas, c("Cebus_versicolor", "Cebus_leucocephalus", "Cebus_malitiosus", "Cebus_albifrons"))
# # load trees
# primTree <- read.nexus("C:/Users/pgalante/Projects/NASA/changeRangeR/calcComparisons/output.nex")
# # Rename albifrons_all created from sum of "Cebus_versicolor", "Cebus_leucocephalus", "Cebus_malitiosus", "Cebus_albifrons"
# names(primRas[[15]]) <- "Cebus_albifrons"
# ##################  Change rasternames to match tree tip names
# #        (1) Cheracebus_lugens to Callicebus_lugens, (2) Cheracebus_medemi to Callicebus_medemi, (3) Lagothrix_lagothricha to Lagothrix_lagotricha,
# #        (4) Leontocebus_fuscus to Leontocebus_fuscicollis, (5) Plecturocebus_caquetensis to Callicebus_caquetensis, (6) Plecturocebus_discolor to Callicebus_discolor,
# #        (7) Plecturocebus_ornatus to Callicebus_ornatus, (8) Saimiri_cassiquiarensis to Saimiri_sciureus
# sum(names(primRas) %in% primTree[[500]]$tip.label)
# oldnames <- names(primRas)
# namesToChange <- c("Cheracebus_lugens", "Cheracebus_medemi", "Lagothrix_lagothricha", "Leontocebus_fuscus", "Plecturocebus_caquetensis", "Plecturocebus_discolor", "Plecturocebus_ornatus", "Saimiri_cassiquiarensis")
# newNames <- c("Callicebus_lugens", "Callicebus_medemi", "Lagothrix_lagotricha", "Leontocebus_fuscicollis", "Callicebus_caquetensis", "Callicebus_discolor", "Callicebus_ornatus", "Saimiri_sciureus")
# oldnames[oldnames %in% namesToChange] <- newNames
# names(primRas) <- oldnames
# names(primRas[[15]]) <- "Cebus_albifrons"
#
#
# ###  Crop rasters to minimum extent
# allRas <- sum(primRas, na.rm=T)
# #allRas[allRas == 0] <- NA
# #allRas <- trim(allRas)
# e <- extent(c(-79, -66, -4.25, 11.8))
# allRas <- crop(allRas, e)
# plot(allRas)
# primRas <- crop(primRas, allRas)
# rm(allRas)
#
# ## Species richness
# SR <- sum(primRas, na.rm=T)
# writeRaster(SR, "/home/pgalante/Projects/NASA/changeRangeR/calcComparisons/SR.tif")
#
# ## Species endemism
# source("C:\\Users\\pgalante\\Git\\changeRangeR\\R/calc_PE.R")
# library(phylobase)
# library(dplyr)
# #primRas2 <- stack(aggregate(primRas, fact = 10, fun = max))
# Allxy <- rasterToPoints(primRas)
# sites <- as.data.frame(Allxy[,1:ncol(Allxy)])
# sites[is.na(sites)] <- 0
#
# pEprims <- calc_PE(tree = primTree[[1]], sites_x_tips = sites, presence = "presence")
#
# cd <- read.csv('C:/Users/pgalante/Projects/NASA/changeRangeR/matrix_primate_composition.csv')
# cdColnames <- colnames(cd)
# cdColnames <- recode(cdColnames, Cheracebus_lugens = "Callicebus_lugens", Cheracebus_medemi = "Callicebus_medemi", Lagothrix_lagothricha = "Lagothrix_lagotricha",
#                      Leontocebus_fuscus = "Leontocebus_fuscicollis", Plecturocebus_discolor = "Callicebus_discolor", Plecturocebus_ornatus = "Callicebus_ornatus",
#                      Saimiri_cassiquiarensis = "Saimiri_sciureus", Plecturocebus_caquetensis = "Callicebus_caquetensis")
# colnames(cd) <- cdColnames
#
# pEprims <- calc_PE(tree = primTree[[1]], sites_x_tips = cd[,4:ncol(cd)], presence = "presence")
#
# # SR <- sum(primRas, na.rm=T)
# # SR[SR == 0] <- NA
# # for (i in 1:length(SR[SR > 0])){
# #   SR[SR > 0][i]<-pEprims$PE[i]
# # }
#
# sum(pEprims$PE)
#
# rxyz <- cbind(Allxy[,1:2], pEprims$PE)
# dfr <- rasterFromXYZ(rxyz)
# plot(dfr)
# dfr[dfr == 0] <- NA
#
# writeRaster(dfr, "C:/Users/pgalante/Garbage/PETree1TestR.tif")
#
#
# #############################################################################################
# ######  Calculate SE
# library(SSDM)
# library(raster)
# library(ape)
# library(dplyr)
# ## Using Nivel 2 binary maps (34 species) calculate: PE, PD using trees 1, 300, 500, 1000
# #### Among trees, get a measure of variation
# # load rasters and edit names
# primRas <- stack(list.files(path = "C:/Users/pgalante/Projects/NASA/changeRangeR/Nivel2", pattern = "\\.tif$", full.names = T))
# names(primRas) <- gsub("*_veg", "", names(primRas))
# # Drop models that are resolved into one group; Cebus_albifrons
# primRas <- dropLayer(primRas, c("Cebus_versicolor", "Cebus_leucocephalus", "Cebus_malitiosus", "Cebus_albifrons"))
# # load trees
# primTree <- read.nexus("C:/Users/pgalante/Projects/NASA/changeRangeR/PD/output.nex")
# # Rename albifrons_all created from sum of "Cebus_versicolor", "Cebus_leucocephalus", "Cebus_malitiosus", "Cebus_albifrons"
# names(primRas[[15]]) <- "Cebus_albifrons"
# primRas$Cebus_albifrons[primRas$Cebus_albifrons == 0] <- NA
# ##################  Change rasternames to match tree tip names
# #        (1) Cheracebus_lugens to Callicebus_lugens, (2) Cheracebus_medemi to Callicebus_medemi, (3) Lagothrix_lagothricha to Lagothrix_lagotricha,
# #        (4) Leontocebus_fuscus to Leontocebus_fuscicollis, (5) Plecturocebus_caquetensis to Callicebus_caquetensis, (6) Plecturocebus_discolor to Callicebus_discolor,
# #        (7) Plecturocebus_ornatus to Callicebus_ornatus, (8) Saimiri_cassiquiarensis to Saimiri_sciureus
# sum(names(primRas) %in% primTree[[500]]$tip.label)
# oldnames <- names(primRas)
# namesToChange <- c("Cheracebus_lugens", "Cheracebus_medemi", "Lagothrix_lagothricha", "Leontocebus_fuscus", "Plecturocebus_caquetensis", "Plecturocebus_discolor", "Plecturocebus_ornatus", "Saimiri_cassiquiarensis")
# newNames <- c("Callicebus_lugens", "Callicebus_medemi", "Lagothrix_lagotricha", "Leontocebus_fuscicollis", "Callicebus_caquetensis", "Callicebus_discolor", "Callicebus_ornatus", "Saimiri_sciureus")
# oldnames[oldnames %in% namesToChange] <- newNames
# names(primRas) <- oldnames
# primRas[primRas == 0] <- NA
#
# plot(primRas$Alouatta_seniculus)
#
# ### For each cell, species richness / total # cells of each species present
# ## Calculate SR
# SR <- sum(primRas, na.rm=T)
#
#  # Calculate # pixels per species
# areaList <- data.frame(matrix(ncol = 2, nrow = nlayers(primRas)))
# colnames(areaList) <- c("species", "area")
# for(i in 1:nlayers(primRas)){
#   areaPrimate <- ncell(primRas[[i]][!is.na(primRas[[i]])])
#   areaList$species[[i]] <- names(primRas[[i]])
#   areaList$area[[i]] <- areaPrimate
# }
# # Get species per pixel
# rasDF <- rasterToPoints(primRas)
# SEval <- list()
# for (i in 1:nrow(rasDF)){
#   primsPixel <- as.matrix(Filter(function(x)!all(is.na(x)), rasDF[30,3:nlayers(primRas)]))
#   areaList2 <- filter(areaList, areaList$species %in% rownames(primsPixel))
#   SEval[i] <- values(SR)[i] / sum(areaList2$area)
# }
#
