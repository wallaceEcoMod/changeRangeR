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

#' @export
setupSparseSpeciesMatrixStuff=function(sumDirs,
																			 scenario='Present',
																			 env,
																			 allSpeciesMaps,
																			 nCellChunks=10,
																			 mc.cores=1,
																			 myTempPath=rasterOptions()$tmpdir){

	#  for testing
	#  scenario='Present'; nCellChunks=10; myTempPath=rasterOptions()$tmpdir
	
	# make summaryBaseDir inside? no, because you may use that without running the cellBySpeciesMatrics
	
	#=============================================================
	# species index table
		# columns: species name, integer index
	sp.ind=speciesIndexTable(allSpeciesMaps,sumDirs)
	message('wrote species index table')
	#================================================================
	# cell Index Table
		# columns: long, lat, cellid
	cell.ind=cellIndexTable(env,nCellChunks,sumDirs)
	message('wrote cell index table')

	#================================================================
	# cell by species list	
		# generate cell by species matrices, to make everything downstream faster
	# # write a single file for each cell NOPE
	# 	# write out cell ids for each species as an intermediate product
	# 	# list, with each element containing the indices of species occurring there
	# 	# all cells on the grid have a list element, even those that are NA, to make the indexing easeier
	# try a looooong table with columns of species index and cell index
		# make multiple in chunks and concatenate them later


	#mapDir=paste0(summaryBaseDir,'/SpeciesConcensusMaps')
	#scenarioDirs=c(list.files(mapDir,full.names=T),'/Users/ctg/Documents/SDMs/BIEN41/NWPlants/BinaryMaps/Present/BinaryMaps')

	#for(i in 1:length(inputsFromSDMWorkflow$scenarios)){
		#scenario=inputsFromSDMWorkflow$scenarios[i]
		#print(scenario)
		#if(length(inputsFromSDMWorkflow$binaryMapDirs[[scenario]])>0){
		message('writing cell by species matrices; this can be slow...')

			cellBySpeciesMatrices(outDir=sumDirs$cbsDir, allSpeciesMaps=allSpeciesMaps,scenario=scenario, sp.ind=sp.ind,cell.ind=cell.ind, nCellChunks= nCellChunks, removeTempFiles=FALSE, mc.cores=mc.cores, myTempPath=myTempPath,overwrite=F)
		#}
	#}

}