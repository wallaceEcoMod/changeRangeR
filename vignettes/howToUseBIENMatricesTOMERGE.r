library(parallel)
library(raster)
library(rgdal)
library(Matrix)
library(textTinyR)
mc.cores=1 # choose whatever your computer can handle. i use 1 less than available.
# set this for whereever you download the data
myBaseDir='/Users/ctg/Documents/SDMs/BIEN41/NWPlants_BinaryOnly/BIEN41_Summaries'

# there are 3 main tables that connect species to the locations they occur in
# 1. cell index table
cell.ind=read.csv(paste0(myBaseDir,'/cellIndexTable.csv'))
# 2. species index table
sp.ind=read.csv(paste0(myBaseDir,'/speciesIndexTable.csv'))
# 3. cell by species matrices
	# these are stored in chunks of cells because most computers wont handle the whole cell by species matrix at once. So we'll read in the chunks below and subset them before keeping them in memory 
	# this first one lists the scenarios avaliable. we'll focus on the present
scenarios=list.files(paste0(myBaseDir,'/CellSpeciesListsBeforeBug'),full.names=T)
scenarios=scenarios[-grep('csv',scenarios)]
# file names for the cell by species matrice for just the present 
scenario.files=list.files(scenarios[1],full.names=T)

# an additional table can be used to subset species (or cells, analogously)
# here's an example of how to subset a group of species. i want to split the present day ranges by algorithm (PPM, RangeBagging, Points) because i only want to compare the PPM models to future models (since those are the only species we projected into the future.)
myAttrTable=read.csv(paste0(myBaseDir,'/AttributeTables/algorithmType.csv'))
keep=subset(myAttrTable,PPM==1)$index

# the matrix with all species and cells might fit in memory on your computer. But if you're doing a summary, such as summing over columns to get richness, it's better to do that in the loop and avoid putting it all in memory. First I'll show an example of getting richness, and then we'll load the whole matrix in and we can see if it breaks your machine
# roughly 1 min with 1 core 
richByChunk=mclapply(seq_along(scenario.files),function(x){
	print(x)
	tmp=readRDS(scenario.files[x])
	# column names correspond to species indices
	colsToKeep=which(colnames(tmp) %in% keep)
	tmp1=tmp[,colsToKeep]
	data.frame(cellID=as.numeric(rownames(tmp)),rich=textTinyR:: sparse_Sums(tmp, rowSums = T)) # rowSums=F gives column sums if you want to make a statistic for a species
},mc.cores=mc.cores)
rich=do.call('rbind',richByChunk)
# ok, so this structure is kind of a pain in the ass, but you also just made a richness map of 50k+ species in a few seconds....

# plot the richness map
# this env.template is a raster where the grid cell IDs correspond the cell indices stored in cell.ind, so we can put any values in a raster and plot it
env=stack(paste0(myBaseDir,'/Convenience/Env/AllEnv.tif'))
ln=read.csv(paste0(myBaseDir,'/Convenience/Env/layerNames.csv'),header=F)
names(env)=ln[,1]
env.template=env[[1]]; values(env.template)=NA
rich.r=env.template
values(rich.r)[rich$cellID]=rich$rich
# cory colors
cm.cols1=function(x,bias=1) { colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x) }
# it nice to add a shp file especially if you zoom in.
world.shp=readOGR(paste0(myBaseDir,'/Convenience/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS_SIMPL-0.3.shp'))
world.shp=spTransform(world.shp,projection(env))

# make the plot
plot.f='~/Desktop/tmpRich.png'
png(plot.f,h=1480*2, w=1095*2)
	plot(rich.r,col=cm.cols1(100),axis.args=list(cex.axis=3),smallplot= c(.9,.94,.1,.9),axes=F,xlab="",ylab="",xaxt='n',yaxt='n',box=F)
	plot(world.shp,add=T)
dev.off()
system(paste0('open ',plot.f))

# DANGER!
# here's where you're on your own. you can load in the whole cell by species matrix at once and hope that you don't run out of ram. this can take a few minutes, and I have 32gb of ram, so only for the brave or foolish. Really you should figure out a way to subset it for your analysis as above to be practical about working with it.
out=mclapply(seq_along(scenario.files),function(x){
	print(x)
	tmp=readRDS(scenario.files[x])
	# column names correspond to species indices
	colsToKeep=which(colnames(tmp) %in% keep)
	tmp[,colsToKeep]
},mc.cores=mc.cores)

# This makes the species by cell map for all species and cells
# you can do your summary statistics on this object
cbs=do.call('rbind',out)
rm(out) # cuz it's big

# see what a single species map looks like
mySp=env.template
values(mySp)[as.numeric(rownames(cbs))]=cbs[,100]
mySp=raster::trim(mySp)
plot(mySp)


