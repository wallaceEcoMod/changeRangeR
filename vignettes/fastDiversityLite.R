## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache=T)

## ----message=FALSE, warning=FALSE----------------------------------------
library(changeRangeR)
library(raster)
library(rasterVis)
library(rgdal)
library(Matrix.utils)
library(tidyverse)
if(Sys.info()["sysname"]== "Windows") library (parallelsugar)
mc.cores=6

## ------------------------------------------------------------------------
### determine where you want outputs
summaryBaseDir='/Volumes/cm2/changeRangerDemos/trees190_test3'
#summaryBaseDir='/Volumes/cm2/changeRangerDemos/trees190/selec'
if(!file.exists(summaryBaseDir)) dir.create(summaryBaseDir)

# Indicate scenarios. these names will be used throughout to structure folders 
allScen=c('present','8580')
#load environnent with the reprojected projection (typically the one you used for modeling). You only need the raster grid that , not the actual layer values
envGrid=raster::stack(system.file("extdata/treeDemo/envGrid.tif",package='changeRangeR'))

# folder of binary range rasters
myDir=system.file("extdata/treeDemo/BinaryMaps",package='changeRangeR')

# shapefiles for plotting. This one comes preinstallted
world.shp=readOGR(system.file(
  "extdata/treeDemo/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS_SIMPL-0.3.shp",
  package='changeRangeR'),'TM_WORLD_BORDERS_SIMPL-0.3',verbose=F)
world.shp2=spTransform(world.shp,projection(envGrid))

## ------------------------------------------------------------------------
sumDirs=setupSummaryDirectories(summaryBaseDir, optionalSubDirs=c('funcDiv','phyloDiv'))
str(sumDirs)

## ----message=F-----------------------------------------------------------
allSpeciesMaps=tibble(rasterFiles=list.files(paste0(myDir,'/present'),
                                             recursive=T, full.names=T)) %>% 
  mutate(sp.names= rasterFiles %>% basename %>% file_path_sans_ext) %>%
  separate(sp.names,into=c('g','sp','t','s'),sep='_') %>% select(-t,-s) %>%
  unite(sp.names,g,sp)
# species index table. columns: species name, integer index
sp.ind=speciesIndexTable(allSpeciesMaps,sumDirs)
# cell Index Table. columns: long, lat, cellid
cell.ind=cellIndexTable(envGrid,nCellChunks=10,sumDirs)

for (scn in allScen){
  allSpeciesMaps=tibble(rasterFiles=list.files(paste0(myDir,'/',scn),
                                               recursive=T, full.names=T)) %>% 
    mutate(sp.names= rasterFiles %>% basename %>% file_path_sans_ext) %>%
    separate(sp.names,into=c('g','sp','t','s'),sep='_') %>% select(-t,-s) %>%
    unite(sp.names,g,sp)
  cellBySpeciesMatrices(sumDirs$cbsDir,
	                      allSpeciesMaps=allSpeciesMaps,
	                      scenario=scn,
	                      envGrid=envGrid,
	                      sp.ind=sp.ind,
                        cell.ind=cell.ind,
	                      nCellChunks=10,
	                      mc.cores=mc.cores,
	                      overwrite=T)
}

## ------------------------------------------------------------------------
cell.ind=readRDS(paste0(sumDirs$myBaseDir,'/cellIndexTable.rds'))
head(cell.ind)

## ------------------------------------------------------------------------
sp.ind=readRDS(paste0(sumDirs$myBaseDir,'/speciesIndexTable.rds'))
head(sp.ind)

## ------------------------------------------------------------------------
cell.ind=readRDS(paste0(sumDirs$myBaseDir,'/cellIndexTable.rds'))
chunks.r=sparseToRaster(cell.ind,envGrid,'chunkID')
fdMapPlot(chunks.r,shp=world.shp2)

## ------------------------------------------------------------------------
rich=lapply(allScen,function(scn){
  r=richnessFromCBS(cbsDir=sumDirs$cbsDir,
                    scenario=scn,env=envGrid,
                    mc.cores=mc.cores, outDir=sumDirs$richDir)
  fdMapPlot(r,paste0(sumDirs$figDir,'/Richness_',scn,'.pdf'),shp=world.shp2)
  r
})
fdMapPlot(rich[[1]],shp=world.shp2)

## ------------------------------------------------------------------------
sp.ind=readRDS(paste0(sumDirs$sumBaseDir,'/speciesIndexTable.rds'))
ra=lapply(allScen,function(scn){
	    rangeArea(cbsDir=sumDirs$cbsDir,scenario=scn,
	              sp.ind=sp.ind,outDir=sumDirs$rangeSizeDir,mc.cores=mc.cores)
})
str(ra)
hist(ra[[1]]$rangeArea)

## ------------------------------------------------------------------------
sumDirs$rarityDir=file.path(sumDirs$rangeSizeDir,'Rarity')
if(!file.exists(sumDirs$rarityDir)) dir.create(sumDirs$rarityDir)

rar=lapply(allScen,function(scn){
	# generate the value of 1/range size for each species in an attribute table
	raritySpAttr=sumDirs$rangeSizeDir %>% 
	  list.files(full.names=T,pattern=scn) %>% 
	  readRDS %>% 
	  mutate(rarity=1/rangeArea) %>% 
	  select(-rangeArea)
	r=speciesAttributeByCell(cbsDir=sumDirs$cbsDir,scenario=scn,
	                         attrTable=raritySpAttr, method='mean', 
	                         env=envGrid, outDir=sumDirs$rarityDir)	
	fdMapPlot(log(r),plotFile=paste0(sumDirs$figDir,'/Rarity_',scn,'.pdf'), 
	          shp=world.shp2,legend.args=list(text='log(rarity)',line=2,side=4))
	r
})
rar.st=stack(rar)
fdMapPlot(log(rar.st[[1]]),shp=world.shp2,legend.args=list(text='log(rarity)'))

