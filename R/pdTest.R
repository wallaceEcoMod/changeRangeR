#
# library(picante)
# library(raster)
#
# # samp: community data matrix
# # tree: Phylo tree object
# # include.root = T
#
# ### Read in rasters and adjust names to match tree including some collapsed nodes
# primRas <- stack(list.files(path = "/home/pgalante/Projects/NASA/changeRangeR/Nivel2", pattern = "\\.tif$", full.names = T))
# names(primRas) <- gsub("*_veg", "", names(primRas))
# # Drop models that are resolved into one group; Cebus_albifrons
# primRas <- dropLayer(primRas, c("Cebus_versicolor", "Cebus_leucocephalus", "Cebus_malitiosus", "Cebus_albifrons"))
# oldnames <- names(primRas)
# namesToChange <- c("Cheracebus_lugens", "Cheracebus_medemi", "Lagothrix_lagothricha", "Leontocebus_fuscus", "Plecturocebus_caquetensis", "Plecturocebus_discolor", "Plecturocebus_ornatus", "Saimiri_cassiquiarensis")
# newNames <- c("Callicebus_lugens", "Callicebus_medemi", "Lagothrix_lagotricha", "Leontocebus_fuscicollis", "Callicebus_caquetensis", "Callicebus_discolor", "Callicebus_ornatus", "Saimiri_sciureus")
# oldnames[oldnames %in% namesToChange] <- newNames
# names(primRas) <- oldnames
# ### Read in tree
# primTree <- read.nexus("/home/pgalante/Projects/NASA/changeRangeR/calcComparisons/output.nex")
# ## Get mask raster of colombia
# col <- raster::getData(name = "alt", country = "COL")
# col <- crop(col, primRas)
# primRas <- crop(primRas, col)
# ### Create community data
# plot(primRas[[1]])
# # set all cells =NA to 0, then mask by Colombia
# for (i in 1:nlayers(primRas)){
#   primRas[[i]][is.na(primRas[[i]])] <- 0
# }
# primRas <- mask(primRas, col)
# # convert raster to dataframe
# commDat <- rasterToPoints(primRas)
# rownames(commDat) <- 1:nrow(commDat)
#
# ### pd
# phydiv <- pd(commDat[,3:ncol(commDat)], primTree[[500]], include.root=F)
# colnames(commDat[,3:ncol(commDat)])
# primTree[[1]]$tip.label
# sum(phydiv$SR)
#
# phydiv$x <- commDat[,1]
# phydiv$y <- commDat[,2]
# phydiv <- phydiv[,c(3,4,1,2)]
# phyRas <- rasterFromXYZ(phydiv)
# plot(phyRas[[2]])
#
# writeRaster(phyRas[[1]], '/home/pgalante/Garbage/PD/pdtree500R.tif')
#
# r1 <- raster('/home/pgalante/Garbage/PD/pdtree1R.tif')
# r300 <- raster("/home/pgalante/Garbage/PD/pdtree300R.tif")
# r500 <- raster("/home/pgalante/Garbage/PD/pdtree500R.tif")
#
# library(dismo)
# nicheOverlap(r300, r500, stat = "D")
