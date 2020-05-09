#===================================================================
#===================================================================
#===================================================================
#' @description Conver rasters to long format data frames. Optionally reprojectt rasters first.
#' @note The reason to reproject internally is that it may be inconvenient to store tons of maps in a different projection, so this just handles the reprojection internally
#' @export
speciesByCellLongFormat=function(scenario,
																 f,
																 env,
																 nChunks,
																 sp.names,
																 mc.cores,
																 outDir,
																 sp.ind,
																 myTempPath,
																 reprojectToEnv=F,
																 verbose=T){

	#  for testing
	#  f=allSpeciesMaps$rasterFiles; nChunks=nCellChunks; scenario='Present'; sp.names=allSpeciesMaps$sp.names; reprojectToEnv=F; outDir=sumDirs$cbsDir

	t1=proc.time()

	# set up directories
	message(paste0(length(f),' species'))
	outDir2=paste0(outDir,'/',scenario)
	if(!file.exists(outDir2)) dir.create(outDir2)
	problem.f=paste0(outDir,'/',scenario,'_problems.csv')
	if(!file.exists(problem.f)) 	file.create(problem.f)
	# !! somehow not printing the right name for errors

	# split tasks
	if(nChunks>1) {
		tasks=chopTasks(f,nChunks)
		tasks.sp.name=chopTasks(sp.names,nChunks)
	} else{
		tasks=list(f)
		tasks.sp.name=list(sp.names)
	}

	message(paste0(nChunks ,' chunks each with ',length(tasks[[1]]),' species'))

	lapply(seq_along(tasks),function(x){ # loop over nChunks
	  if(verbose) message('chunk ',x,' \n')
		defaultRasterTmpDir=raster::rasterOptions(overwrite=T)$tmpdir #reset later
		myRasterTmpDir=paste0(myTempPath,'/tmp_',x)
		if(!file.exists(myRasterTmpDir)) dir.create(myRasterTmpDir)
		raster::rasterOptions(tmpdir=myRasterTmpDir)

		# internal function  -------------------------------------------------
		# reproject if needed, and make the matrix for a single species
		longFun=function(ii,tasks.sp.name,sp.ind){
			#if(verbose) cat(ii,' ')
			sp.name=tasks.sp.name[[x]][ii]
			this.sp.ind=sp.ind$index[which(sp.ind$species %in% sp.name)]
			if(length(this.sp.ind)>1){
				warning(paste0('There are multiple rasters for this species: ',sp.name,". Because we use names to match files, multiple files with the same basename are not supported. As a (possibly bad) default I'm selecting just the first instance of this species"))
				this.sp.ind=this.sp.ind[1]
			}

			# error handling
			if(length(this.sp.ind)==0) {
				print(paste0(sp.name,' ',x,' ',ii,' does not have a species index '))
				sink(problem.f,append=T)
				cat(sp.name, 'speciesName', sep=','); cat('\n')
				sink()
				return(NULL)
			}

			map=try(raster(tasks[[x]][ii]))
			if(class(map)=='try-error'){
				print(paste0(sp.name,' ',x,' ',ii,' cannot be read by raster '))
				sink(problem.f,append=T)
				cat(sp.name, 'speciesName', sep=','); cat('\n')
				sink()
				return(NULL)
			}

			if(suppressWarnings(is.na(maxValue(map)))){
				print(paste0(sp.name,' ',x,' ',ii,' map is all NA'))
				sink(problem.f,append=T)
				cat(sp.name, 'speciesName', sep=','); cat('\n')
				sink()
				return(NULL)
			}

			if(maxValue(map)==0){
				print(paste0(sp.name,' ',x,' ',ii,' map is all zeroes '))
				sink(problem.f,append=T)
				cat(sp.name, 'speciesName', sep=','); cat('\n')
				sink()
				return(NULL)
			}

			if(reprojectToEnv){
				map2=raster::projectRaster(map,env)
				# For some reason, ngb doesn't work for small rasters and just sets the values to zero. so using bilinear to get nonzero values, but than choosing the top N cells to set to 1 so range size doesn't change
				originalRangeSize=sum(values(map)==1,na.rm=T)
				fuck=sort(values(map2),decreasing=T)
				keep=min(originalRangeSize,length(fuck))
				map2=map2>=fuck[keep]
			} else {# this extend doesn't guarantee the same extent
				map2=raster::extend(map,extent(env))
			}

			cellID=which(values(map2)>0)

			if(length(cellID)==0) {# if no cells
				print(sp.name);
				sink(problem.f,append=T)
				cat(sp.name, 'cellsNA',sep=','); cat('\n')
				sink()
				return()
			}

			rm(map2)
			if(length(c(cellID,rep(this.sp.ind,length(cellID))))==3) stop()
			cc=matrix(c(cellID,rep(this.sp.ind,length(cellID))),ncol=2)
			cc
		} # end internal function

	 # parallelize ------------------------------------------------------
	 if(Sys.info()["sysname"]== "Windows"){
		# foreach is just too finicky for my patience...
		# 			cl = parallel::makeCluster(mc.cores)
		#     	doSNOW::registerDoSNOW(cl)
		# 			pb = txtProgressBar(max = nrow(sp.ind), style = 3)
		# 			progress = function(n) setTxtProgressBar(pb, n)
		#     	out1= foreach::foreach(ii = 1:nrow(sp.ind),
		#     		#.export= c('sp.ind','tasks.sp.name'),
		# 				.options.snow=list(progress = progress),.combine='rbind') %dopar% {
		# 				# !!PEP!!
		# 					# for some reason tasks.sp.name is not being passed. its in the environmetn and it works interactively. just not when i run the function speciesListLongFormat
		# 					print(ii)
		# 					return(longFun(ii,tasks.sp.name=tasks.sp.name,sp.ind=sp.ind))
		# 			}
		# 			parallel::stopCluster(cl)
		#     	foreach::registerDoSEQ()
			out=parallelsugar::mclapply(seq_along(tasks[[x]]),function(ii){
				longFun(ii,tasks.sp.name,sp.ind)
			},mc.cores=mc.cores)
		} else {	# mac
			t3=proc.time()
			out=parallel::mclapply(seq_along(tasks[[x]]),function(ii){
			#for(ii in seq_along(tasks[[x]])){ # for testing
				#print(ii)
				longFun(ii,tasks.sp.name,sp.ind)
			},mc.cores=mc.cores)
			(t4=proc.time()-t3)
			# tried this too, but doesn't matter on a mac
			#  t1=proc.time()
			#  subTasks=chopTasks(seq_along(tasks[[x]]),mc.cores)
			#  out=mclapply(seq_along(subTasks),function(ii){
			#  #for(ii in seq_along(tasks[[x]])){ # for testing
			# 	 #print(ii)
			# 	 a=lapply(subTasks[[ii]],function(yy){longFun(yy,tasks.sp.name, sp.ind)})
			# 	 do.call('rbind',a)
			#  },mc.cores=mc.cores)
			#  (t2=proc.time()-t1)
		}
		out1=do.call('rbind',out)


		# outputs  ------------------------------------------------------
		colnames(out1)=c('cellID','spID')
		saveRDS(out1,file=paste0(outDir2,'/temp_long_chunk_',x,'.rds'))
		unlink(raster::rasterOptions(overwrite=T)$tmpdir,recursive=T) # clean up
		raster::rasterOptions(tmpdir=defaultRasterTmpDir)
		gc()
	})

	t2=(proc.time()['elapsed']-t1['elapsed'])/60
	message(paste('Ran in ',round(t2,2),' minutes'))
}


#===================================================================
#===================================================================
#===================================================================
#' @export
cellBySpeciesMatrices=function(outDir,
														   allSpeciesMaps,
														   scenario,
														   sp.ind,
														   cell.ind,
														   env,
														   nCellChunks,
														   removeTempFiles=FALSE,
														   mc.cores=mc.cores,
														   overwrite=FALSE,
														   myTempPath=rasterOptions()$tmpdir,
														   reprojectToEnv=F,
														   verbose=T){

	# Strategy
		# internally maintain an easy and quick version when chunks aren't needed
		# separately maintain a more elaborate version with chunks
		# to make parallelization optimized for windows and mac, write a single function and apply it with foreach or mclapply

	#  for testing
	#  scenario='Present'; outDir=sumDirs$cbsDir; myTempPath='~/Desktop'; overwrite=F; removeTempFiles=F; verbose=T; nCellChunks=6

	# intermediate step of long format ---------------------------------
	scenDir=paste0(outDir,'/',scenario)
	if(!file.exists(scenDir)) dir.create(scenDir)

	chunk.f=list.files(scenDir,full.names=T,pattern='temp_long')
	if(length(chunk.f)==0 | overwrite ){
		# testing for ordering of names
			# 		f.sp=unlist(lapply(seq_along(f),function(x) strsplit(basename(f[x]),'__')[[1]][2]))
			# 		f.ind=which(sp.ind$species %in% f.sp)
		#save(scenario,allSpeciesMaps, cell.ind, env, mc.cores,outDir,sp.ind,myTempPath, reprojectToEnv, verbose=verbose,file='~/Desktop/forPep.rdata')
		speciesByCellLongFormat(scenario,f=allSpeciesMaps$rasterFiles, nChunks=max(cell.ind$chunkID), env=env, sp.names=allSpeciesMaps$sp.names,mc.cores=mc.cores,outDir=outDir,sp.ind=sp.ind,myTempPath=myTempPath, reprojectToEnv=reprojectToEnv, verbose=verbose)
	}
	message('done making intermediate long format')

	# 		# for testing that range sizes are correctly recorded
	# 	if(scenario=='Present'){
	# 		chunk.f=list.files(paste0(outDir,scenario), full.names=T,pattern='temp')
	# 		range.f=list.files( '/Users/ctg/Documents/SDMs/BIEN41/NWPlants_BinaryOnly/BIEN41_outputs/PPM/BinaryMaps',full.names=T,recursive=T,pattern='TP05')
	# 		sp.ind=read.csv(paste0(summaryBaseDir, '/speciesIndexTable.csv'),stringsAsFactors=F)
	# 		aa=readRDS(chunk.f[80])
	# 		keep=sample(unique(aa[,1]),50)
	# 		all.sp=sp.ind[sapply(keep,function(yy) which(sp.ind$index==yy)),]
	# 		out=mclapply(1:nrow(all.sp),function(x){
	# 			sp=all.sp$species[x]
	# 			r=raster(range.f[grep(sp,basename(range.f))])
	# 			true=cellStats(r,sum,na.rm=T)
	# 			rsInLong=aa[aa[,1]==all.sp$index[x],]
	# 			keep1=grep(sp, rangeSize$species)
	# 			print(x)
	# 			data.frame(sp=sp,true=true,inLong=nrow(rsInLong), cbsRangesize=rangeSize[keep1,3])
	# 		},mc.cores=5)
	# 		(out1=do.call(rbind,out))
	# 	}
	#

	# make sparse matrices -------------------------------------------

	cbsFun=function(chunkF,overwrite){
		if(verbose) message(chunkF)

		#check if output file exists
		outFile=paste0(outDir,'/',scenario,'/chunk_',chunkF,'.rds')
		# may add this back in; wasn't working
		#if (file.exists(outFile) & !overwrite) {message (paste('chunk',chunkF,'already exists. Check it out')); return(NULL)}

		#Select cells that appear in this specific chunk, chunkF
		cellsIDinChunk = cell.ind$cellID[which(cell.ind$chunkID==chunkF)]

		#select species occurences in cells of this chunk, we remove files for memory purposes
		#sco= readRDS(paste0(outDir2,'/temp_long_chunk_1.rds'))
		outDir2=paste0(outDir,'/',scenario)
		scoDT = readRDS(paste0(outDir2,'/temp_long_chunk_1.rds')) %>% data.table::as.data.table()
		#rm(sco); gc(verbose = F)
		spOccCellChunk = scoDT[scoDT$cellID %in% cellsIDinChunk,]
		if(nrow(spOccCellChunk)==0) saveRDS(NULL,outFile)

		rm(scoDT); gc(verbose = F)
		#these are the indices of the rows in the cell id just for this chunk (cellsIDinChunk) (not the rows of values(env), because NAs were removed)
		cellsIDinChunk.row.index=match(spOccCellChunk$cellID,cellsIDinChunk)
		# this would be the row in values(env)
		#cellChunkID = cell.ind$cellID[match (spOccCellChunk$cellID,cell.ind$cellID)]
		#build a matrix of species occurrence with the cellID OF THIS SPECIFIC CHUNK
		#matrixChunkLocations =  matrix (c(cellChunkID, spOccCellChunk$spID),ncol = 2,byrow = F)
		matrixChunkLocations =  matrix (c(cellsIDinChunk.row.index, spOccCellChunk$spID),ncol = 2,byrow = F)

		#Alternatively you could use merge, but it is very time consuming compared to match
		#we can apply match above bc cell.ind works as a dictiornary (one entry per cellID)
		# matrixChunkLocations =  merge (spOccCellChunk,
		#                                cell.ind[,c('cellID','cellChunkID')],
		#                                by='cellID',all.x=T)
		# m = as.matrix (matrixChunkLocations[,c('chunkcellID','spID')])


		#create sparse matrix
		sM <- Matrix::Matrix(data = 0,
												 nrow = length(cellsIDinChunk),
												 ncol = nrow(sp.ind),
												 dimnames = list(as.character(cellsIDinChunk), as.character(sp.ind[,2])),
												 sparse = T)

		#add 1 where locations of occurrence in matrix
		sM[matrixChunkLocations] <- 1
		rm(matrixChunkLocations); gc(verbose = F)

		#write output
		saveRDS(sM,outFile)
	}

	#if(nCellChunks==1){
	message ('Building sparse matrices')

	# parallelize ---------------------------------------------
	#using mclapply
	#if (Sys.info()["sysname"]== "Windows") {mclapply <- parallelsugar::mclapply}
	if(Sys.info()["sysname"]== "Windows"){
		t1=proc.time()
		cl = parallel::makeCluster(mc.cores)
		doSNOW::registerDoSNOW(cl)
		pb = txtProgressBar(max = nCellChunks, style = 3)
		progress = function(n) setTxtProgressBar(pb, n)
		opts = list(progress = progress)
		lp = (.packages())
		foreach::foreach(chunkF = 1:nCellChunks,.options.snow=opts,.packages = lp) %dopar% { cbsFun(chunkF,overwrite=overwrite)}
		snow::stopCluster(cl)
		(t2=proc.time()-t1)
	} else {
		t1=proc.time()
		parallel::mclapply(1:nCellChunks,cbsFun,overwrite=overwrite)
		t2=(proc.time()['elapsed']-t1['elapsed'])/60
		message(paste('Ran in ',round(t2,2),' minutes'))
	}

	#}	# end nCellChunk==1

	#=============================================================
	#=============================================================
	# if(nCellChunks>1){
#
# 				parallel::mclapply(1:nCellChunks,cbsFun)
#
#
#
# 	# Cory's old way
# 		# for each chunk of cells, read in each chunk of species and grab just the ones associated with the cells. this is a little inefficient because you have to duplicate reading in the species nChunks times, but its easier to code and never requires reading all data in at once.
# 		chunk.f=list.files(paste0(outDir,'/',scenario),full.names=T, pattern='temp_long')
# 		mclapply(1:nCellChunks,function(x){
# 			message(paste0('chunk ',x,' of ',nCellChunks))
# 			outFile=paste0(outDir,'/',scenario,'/chunk_',x,'.rds')
# 			if(!overwrite) { if(file.exists(outFile)) return()}
# 			cellsThisChunk=cell.ind$cellID[cell.ind$chunkID==x]
# 			# make an empty matrix with cells x row
# 				# the rows will be a subset of all cells, but the columns will include all species
# 			cellBySp=matrix(0,nrow=length(cellsThisChunk),ncol=nrow(sp.ind))
# 			colnames(cellBySp)=sp.ind[,2]
# 			rownames(cellBySp)=cellsThisChunk
# 			for(j in seq_along(chunk.f)){
# 				if(verbose) cat(j,' ')
# 				tmp.dat=readRDS(chunk.f[j])
# 				keep=which(c(tmp.dat[,2]) %in% cellsThisChunk)
# 				tmp.dat1=matrix(tmp.dat[keep,],ncol=2)
# 				# becuase there's no multiple match function.
# 				fuck=data.frame(tmp.dat1, cellIndex=do.call(c,lapply(tmp.dat1[,2],function(x) which(x==cellsThisChunk))))
# 				spIDs=unique(fuck[,1])
# 				for(ii in spIDs){
# 					print(ii)
# 					tmp10=fuck[fuck[,1]==ii,]
# 					cellBySp[tmp10$cellIndex,as.numeric(ii)]=1
# 				}
# 			}
# 			cbs=Matrix(cellBySp, sparse = TRUE)
# 			saveRDS(cbs,outFile)
# 			message(paste0('chunk ',x,' done'))
# 		}, mc.cores=mc.cores)
# 		if(removeTempFiles) file.remove(chunk.f)
# 	} # end multichunk
#

}
