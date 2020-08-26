#============================================================
#============================================================
#' Initialize a set of directories to store summary outputs
#' @description Make a standardized set of directories
#' @param summaryBaseDir a file path
#' @param optionalSubDirs a vector of names (strings) of optional directories the user would like to include.
#' @examples
#' setupSummaryDirectories('~/tmpChangeRanger')
#' @return a list of file paths
#' @author Cory Merow <cory.merow@@gmail.com>
# @note 
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation
#' @export
setupSummaryDirectories=function(summaryBaseDir,
																 optionalSubDirs=c('funcDiv','phlyoDiv')){
	# this is a string that can be sub'd in so you can share code among computers
	myBaseDir=summaryBaseDir
	if(!file.exists(summaryBaseDir)) dir.create(summaryBaseDir)
	figDir=paste0(summaryBaseDir,'/Figures')
	if(!file.exists(figDir)) dir.create(figDir)
	attrDir=paste0(summaryBaseDir,'/AttributeTables')
	if(!file.exists(attrDir)) dir.create(attrDir)
	cbsDir=paste0(summaryBaseDir,'/CellSpeciesLists')
	if(!file.exists(cbsDir)) dir.create(cbsDir)
	rangeSizeDir=paste0(summaryBaseDir,'/RangeSize')
	if(!file.exists(rangeSizeDir)) dir.create(rangeSizeDir)
	richDir=paste0(summaryBaseDir,'/Richness/')
	if(!file.exists(richDir)) dir.create(richDir)
	envDir=paste0(summaryBaseDir,'/Env/')
	if(!file.exists(envDir)) dir.create(envDir)
	miscDir=paste0(summaryBaseDir,'/Misc/')
	if(!file.exists(miscDir)) dir.create(miscDir)
	optional=lapply(optionalSubDirs,function(x){
		d=paste0(summaryBaseDir,'/',x)
		if(!file.exists(d)) dir.create(d)
		d
	})
	names(optional)=paste0(optionalSubDirs,'Dir')
	#uniqueComDir=paste0(summaryBaseDir,'/uniqueCom/')
	#if(!file.exists(uniqueComDir)) dir.create(uniqueComDir)
	
	out=list(myBaseDir=myBaseDir,sumBaseDir=summaryBaseDir, figDir=figDir,attrDir=attrDir,cbsDir=cbsDir, rangeSizeDir=rangeSizeDir,  richDir=richDir,envDir=envDir,miscDir=miscDir,#uniqueComDir=uniqueComDir,
		optional)
	out=.flattenlist(out)
	return(out)
}

#https://stackoverflow.com/questions/16300344/how-to-flatten-a-list-of-lists
.flattenlist <- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  }
}



#============================================================
#============================================================
#============================================================
#' @title Convert directory structure to work on different paths or computers
#' @description change the `dirname` for a previously saved set of directories
#' @param sumDirs a list of directories as created by `setupSummaryDirectories`
#' @param newBaseDir a file path
#' @examples
#' sumDirs=setupSummaryDirectories('~/tmpChangeRanger')
#' convertSummaryBaseDir(sumDirs,'~/tmpChangeRanger2')
#' @return a list of file paths
#' @author Cory Merow <cory.merow@@gmail.com>
#' @export

convertSummaryBaseDir=function(sumDirs,newBaseDir){
	out=lapply(sumDirs,function(x){sub(sumDirs$myBaseDir,newBaseDir,x)})
}



#=============================================================
#=============================================================
# inputs (in preparation for a function producing all summaries)
	# intermediate directories may be needed, with 1 map for each speceis (e.g., the best map)
	# you have to get this set up for each workflow run on your own before using this package. I may have some utilities functions to help for common use cases
	# file tree should look like this
# -  folderOfMaps
#		- 	Scenario 1 (E.g, present)
#				-  Algorithm 1 (e.g., PPM)
# 			-  Algorithm 2	(e.g., Range bagging)
# 	- 	Scenario 2 (e.g., Future 2050)
#				-  Algorithm 1 (e.g., PPM)
# 			-  Algorithm 2	(e.g., Range bagging)
# 	- 	Scenario 3 (e.g., Future 2070)
#				-  Algorithm 1 (e.g., PPM)
# 			-  Algorithm 2	(e.g., Range bagging)

# inputsFromSDMWorkflow=list(
# 	scenarios=c('Present'),
# 	# can add other directories, e.g., from other algorithms. names can be anything you like but scenarios should match those above. 
# 	binaryMapDirs=list( 
# 		Present=list(sdm= list.files('/Volumes/cm2/Aus_Fire/BinaryMaps', recursive=T,full.names=T,pattern='__noBias')
# )))


# this will be a custom way to unpack the species name
	# z=inputsFromSDMWorkflow$binaryMapDirs[[1]]; z1=z[[1]]; x=z1[[1]]
# inputsFromSDMWorkflow$sp.names=mclapply(
# 	inputsFromSDMWorkflow$binaryMapDirs,function(z){
# 		lapply(z,function(y){getSpNamesFromDirs(y)})
# 		}, mc.cores=length(inputsFromSDMWorkflow$binaryMapDirs$Present))
		

# simplify to apply this to 1 single scenario at a time. that removes dependenceies on inputFromSDMWorklfow. just need the file  names and the species names. and by default they can be equal.

# allSpeciesMaps=list(rasterFiles=list.files('/Volumes/cm2/Aus_Fire/BinaryMaps', recursive=T,full.names=T,pattern='__noBias'))
# allSpeciesMaps$sp.names=getSpNamesFromDirs(allSpeciesMaps$rasterFiles)
	
	# could choose to concatentate other directories, but they better have the same way to get the names!


#---------------------------------------------------------------------
#' @title Build sparse matrices from rasters
#' @description Convert a folder of species range maps to a sparse matrix
#' @param sumDirs a list of directories as created by `setupSummaryDirectories`
#' @param scenario character; any user-defined terms for scenarios are acceptable. For example, ‘Present’, ’2050’, ’2070’
#' @param env a raster object defining the grid associated with cell indices
#' @param mc.cores number of cores to use
#' @param allSpeciesMaps
#' @param nCellChunks integer
#' @param reprojectToEnv
#' @param myTempPath
#' @param overwrite logical; overwrite existing files?
# @examples
#
#' @return NULL
#' @author Cory Merow <cory.merow@@gmail.com>, Pep Serra-Diaz
# @note 
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the docum
#' @export
setupSparseSpeciesMatrixStuff=function(sumDirs,
																			 scenario='Present',
																			 env,
																			 allSpeciesMaps,
																			 nCellChunks=10,
																			 mc.cores=1,
																			 reprojectToEnv=FALSE,
																			 myTempPath=rasterOptions()$tmpdir,
																			 overwrite=F){

	#  for testing
	#  scenario='Present'; nCellChunks=10; myTempPath=rasterOptions()$tmpdir; overwrite=F
	t1=proc.time()
	#=============================================================
	# species index table
		# columns: species name, integer index
	sp.ind=.speciesIndexTable(allSpeciesMaps,sumDirs)
	message('wrote species index table')
	#================================================================
	# cell Index Table
		# columns: long, lat, cellid
	cell.ind=.cellIndexTable(env,nCellChunks,sumDirs)
	message('wrote cell index table')
	#================================================================
	# cell by species list	
		# generate cell by species matrices, to make everything downstream faster
	# # write out cell ids for each species as an intermediate product
	# make a looooong table with columns of species index and cell index
		# make multiple in chunks and concatenate them later

	message('writing cell by species matrices; this can be slow...')
	out=.cellBySpeciesMatrices(outDir=sumDirs$cbsDir, allSpeciesMaps=allSpeciesMaps,scenario=scenario,env=env, sp.ind=sp.ind,cell.ind=cell.ind, nCellChunks=nCellChunks, removeTempFiles=FALSE, mc.cores=mc.cores, verbose=F, myTempPath=myTempPath,overwrite=overwrite)
	
	t2=proc.time()-t1
	message(paste0(round(t2[3],2),' s'))
	out
}

