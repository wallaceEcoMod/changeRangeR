source("R/SE.R")

# Load binary maps/SDMs
primRas <- stack(list.files(path = paste0(system.file(package="changeRangeR"), "/extdata/DemoData/SDM/colPrimates_binary"), pattern = "\\.tif$", full.names = T))
# Because of species' name resolutions in the phylogeny we use here, we need to rename some of the species in this rasterStack
names(primRas) <- gsub("*_veg", "", names(primRas))
# Drop models that are resolved into one group; Cebus_albifrons
primRas <- dropLayer(primRas, c("Cebus_versicolor", "Cebus_leucocephalus", "Cebus_malitiosus", "Cebus_albifrons"))
oldnames <- names(primRas)
namesToChange <- c("Cebus_all_albifrons", "Cheracebus_lugens", "Cheracebus_medemi", "Lagothrix_lagothricha", "Leontocebus_fuscus", "Plecturocebus_caquetensis", "Plecturocebus_discolor", "Plecturocebus_ornatus", "Saimiri_cassiquiarensis")
newNames <- c("Cebus_albifrons", "Callicebus_lugens", "Callicebus_medemi", "Lagothrix_lagotricha", "Leontocebus_fuscicollis", "Callicebus_caquetensis", "Callicebus_discolor", "Callicebus_ornatus", "Saimiri_sciureus")
oldnames[oldnames %in% namesToChange] <- newNames
names(primRas) <- oldnames

# crop for speed
e <- extent(c(-75, -70, 0, 2.5))
primRas <- crop(primRas, e)

# run function
out <- SE(primRas)

## TESTS
test_that("output type checks", {
  expect_is(out, "raster")
})

