#===================================================================
#===================================================================
#===================================================================
#' @export

# species index table
	# columns: species name, integer index
speciesIndexTable=function(allSpeciesMaps,sumDirs){
	# old way with lists
	#sp.ind=data.frame(species=unique(unlist(inputsFromSDMWorkflow$sp.names, recursive=T)))
	sp.ind=data.frame(species=allSpeciesMaps$sp.names)
	sp.ind$index=1:nrow(sp.ind)
	write.csv(sp.ind,file=paste0(sumDirs$sumBaseDir,'/speciesIndexTable.csv'), row.names=F)
	sp.ind
}

#===================================================================
#===================================================================
#===================================================================
#' @export
# cell Index Table
	# columns: long, lat, cellid
cellIndexTable=function(env,nCellChunks,sumDirs){
	co=coordinates(env)
	keep=apply(values(env),1,function(x) any(!is.na(x)))
	co=co[keep,]
	co=apply(co,2,as.integer)
	# chunks can be based on regions (e.g. if cell values depend on one another) or just random
	chunks=cut(1:nrow(co), nCellChunks, labels = FALSE)
	cell.ind=data.frame(co,cellID=as.integer(cellFromXY(env,co)), chunkID=as.integer(chunks))
	write.csv(cell.ind,file=paste0(sumDirs$sumBaseDir, '/cellIndexTable.csv'),row.names=F)
	cell.ind
}

#===================================================================
#===================================================================
#===================================================================
#' @export
speciesByCellLongFormat=function(scenario,
																 f,
																 nChunks,
																 sp.names,
																 mc.cores,
																 outDir,
																 sp.ind,
																 myTempPath,
																 verbose=T){

	#  for testing
	#  nChunks=max(cell.ind$chunkID)

	print(paste0(length(f),' species'))
	#outDir2=paste0(outDir,gsub('BinaryMaps','Present', basename(scenarioDirs[i])))
	outDir2=paste0(outDir,'/',scenario)
	if(!file.exists(outDir2)) dir.create(outDir2)
	print(outDir2)
	t1=proc.time()
	problem.f=paste0(outDir,'/',scenario,'_problems.csv')
	if(!file.exists(problem.f)) 	file.create(problem.f)
	# !! somehow not printing the right name for errors

	tasks=chopTasks(f,nChunks)
	tasks.sp.name=chopTasks(sp.names,nChunks)
	message(paste0(nChunks ,' chunks each with ',length(tasks[[1]]),' species'))

	lapply(seq_along(tasks),function(x){
		message('chunk ',x,' ====================================\n')
		defaultRasterTmpDir=suppressMessages(rasterOptions()$tmpdir)
		myRasterTmpDir=paste0(myTempPath,'/tmp_',x)
		if(!file.exists(myRasterTmpDir)) dir.create(myRasterTmpDir)
		suppressMessages(rasterOptions(tmpdir=myRasterTmpDir))
		out=mclapply(seq_along(tasks[[x]]),function(ii){
			if(verbose) cat(ii,' ')
			#if(length(grep('concensus',basename(tasks[[x]][ii])))>0){
			#sp.name=strsplit(basename(tasks[[x]][ii]),'__')[[1]][1]
			#} else {


			# UNCOMMENT FOR WORKFLOW!!!!
			# sp.name=strsplit(basename(tasks[[x]][ii]),'__')[[1]][2] #}

			#sp.name=file_path_sans_ext(basename(tasks[[x]][ii]))
			sp.name=tasks.sp.name[[x]][ii]

			this.sp.ind=sp.ind$index[which(sp.ind$species %in% sp.name)]

			if(length(this.sp.ind)==0) {
				print(sp.name)
				sink(problem.f,append=T)
				cat(sp.name, 'speciesName', sep=',')
				cat('\n')
				sink()
				return(NULL)
			}
			map=try(raster(tasks[[x]][ii]))
			if(class(map)=='try-error'){
				print(sp.name)
				sink(problem.f,append=T)
				cat(sp.name, 'speciesName', sep=',')
				cat('\n')
				sink()
				return(NULL)
			}
			map2=raster::extend(map,extent(env))
			cellID=which(values(map2)>0)
			# cellStats(map,sum,na.rm=T)
# 			cellStats(map2,sum,na.rm=T)
# 			length(cellID)
# 			temp=env[[1]]
# 			values(temp)=NA
# 			values(temp)[cc[,2]]= 1
			if(length(cellID)==0) {
				print(sp.name);
				sink(problem.f,append=T)
				cat(sp.name, 'cellsNA',sep=',')
				cat('\n')
				sink()
				return()
			}

			rm(map2)
			cc=matrix(c(rep(this.sp.ind,length(cellID)),cellID),ncol=2)
			cc
		},mc.cores=mc.cores)

		out1=do.call('rbind',out)
		colnames(out1)=c('spID','cellID')
		saveRDS(out1,file=paste0(outDir2,'/temp_chunk_',x,'.rds'))
				# clean up tmp files
		unlink(rasterOptions()$tmpdir,recursive=T)
		suppressMessages(rasterOptions(tmpdir=defaultRasterTmpDir))
		gc()
	})

	t2=(proc.time()-t1)/60
	message(t2)
	#Sys.sleep(180)
}

# aa=out1#readRDS(paste0(outDir2,'/temp_chunk_',x,'.rds'))
# bb=subset(aa,aa[,1]==this.sp.ind)
# bb[,2] %in% cc
#===================================================================
#===================================================================
#===================================================================
#' @export
cellBySpeciesMatrices=function(outDir,
														   allSpeciesMaps,
														   scenario,
														   sp.ind,
														   cell.ind,
														   nCellChunks,
														   removeTempFiles=FALSE,
														   mc.cores=mc.cores,
														   overwrite=FALSE,
														   myTempPath=rasterOptions()$tmpdir,
														   verbose=T){

	#  for testing
	#  scenario='Present'; outDir=sumDirs$cbsDir; myTempPath='~/Desktop'; overwrite=F; removeTempFiles=F; verbose=T

	# intermediate step of long format
	scenDir=paste0(outDir,'/',scenario)
	if(!file.exists(scenDir)) dir.create(scenDir)

	chunk.f=list.files(scenDir,full.names=T,pattern='temp')
	if(length(chunk.f)==0 | overwrite ){

		#f=do.call(c,inputsFromSDMWorkflow$binaryMapDirs[[scenario]])
		#sp.names=do.call(c,inputsFromSDMWorkflow$sp.names[[scenario]])
		#names(f)=NULL
		# testing for ordering of names
			# 		f.sp=unlist(lapply(seq_along(f),function(x) strsplit(basename(f[x]),'__')[[1]][2]))
			# 		f.ind=which(sp.ind$species %in% f.sp)
		speciesByCellLongFormat(scenario,f=allSpeciesMaps$rasterFiles, nChunks=max(cell.ind$chunkID), sp.names=allSpeciesMaps$sp.names,mc.cores=mc.cores,outDir=outDir,sp.ind=sp.ind,myTempPath=myTempPath, verbose=verbose)
		gc()
	}

# 		# for testing
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
	message('done making intermediate long format')

	# concatenate into cell matrices
	# for each chunk of cells, read in each chunk of species and grab just the ones associated with the cells. this is a little inefficient because you have to duplicate reading in the species nChunks times, but its easier to code and never requires reading all data in at once.
	chunk.f=list.files(paste0(outDir,'/',scenario),full.names=T,pattern='temp')
	mclapply(1:nCellChunks,function(x){
		message(paste0('chunk ',x,' of ',nCellChunks))
		outFile=paste0(outDir,'/',scenario,'/chunk_',x,'.rds')
		if(!overwrite) { if(file.exists(outFile)) return()}
		cellsThisChunk=cell.ind$cellID[cell.ind$chunkID==x]
		# make an empty matrix with cells x row
			# the rows will be a subset of all cells, but the columns will include all species
		cellBySp=matrix(0,nrow=length(cellsThisChunk),ncol=nrow(sp.ind))
		colnames(cellBySp)=sp.ind[,2]
		rownames(cellBySp)=cellsThisChunk
		for(j in seq_along(chunk.f)){
			if(verbose) cat(j,' ')
			tmp.dat=readRDS(chunk.f[j])
			keep=which(c(tmp.dat[,2]) %in% cellsThisChunk)
			tmp.dat1=matrix(tmp.dat[keep,],ncol=2)
			# becuase there's no multiple match function.
			fuck=data.frame(tmp.dat1, cellIndex=do.call(c,lapply(tmp.dat1[,2],function(x) which(x==cellsThisChunk))))
			spIDs=unique(fuck[,1])
			for(ii in spIDs){
				tmp10=fuck[fuck[,1]==ii,]
				cellBySp[tmp10$cellIndex,ii]=1
			}
		}
		cbs=Matrix(cellBySp, sparse = TRUE)
		saveRDS(cbs,outFile)
		message(paste0('chunk ',x,' done'))
	}, mc.cores=mc.cores)
	if(removeTempFiles) file.remove(chunk.f)

}


#===================================================================
#===================================================================
#===================================================================

#' @export
richnessFromCBS=function(cbsDir,scenario,env,mc.cores){
	cbs.f=list.files(paste0(cbsDir,'/',scenario),full.names=T)
	toss=grep('temp_',cbs.f)
	cbs.f=cbs.f[-toss]
	richByCell=parallel::mclapply(seq_along(cbs.f), function(x){
		message(x)
		cbs=readRDS(cbs.f[x])
		data.frame(cellID=as.numeric(rownames(cbs)),
		           rich=textTinyR:: sparse_Sums(cbs, rowSums = T))
	},mc.cores=mc.cores)
	rich.vec=do.call('rbind',richByCell)
	rich.r=raster(env[[1]])
	values(rich.r)[rich.vec$cellID]= rich.vec$rich
	rich.r
}

#===================================================================
#===================================================================
#===================================================================
#' @notes optionally uses an antribute table to subset by each column not called `species` or `index` and creates a richness file based on the name of the attribute table column.
#' @export
richnessFromCBSAttr=function(cbsDir,
														 scenario,
														 env,
														 mc.cores,
														 attrTable=NULL,
														 outFileBase){
	cbs.f=list.files(paste0(cbsDir,scenario),full.names=T)
	toss=grep('temp_',cbs.f)
	cbs.f=cbs.f[-toss]
	attrNames=names(attrTable)
	attrNames=attrNames[-mapply(function(x){grep(x,attrNames)}, c('species','index'))]
	out=lapply(seq_along(attrNames),function(y){
		message(attrNames[y])
		keep=attrTable$index[attrTable[attrNames[y]]==1]
		richByCell=mclapply(seq_along(cbs.f), function(x){
			message(x)
			cbs.tmp=readRDS(cbs.f[x])
			cbs=cbs.tmp[,keep]
			#fuck=data.frame(spID=as.numeric(colnames(cbs)),rich=textTinyR:: sparse_Sums(cbs, rowSums = F))
			#print(fuck[fuck[,2]>0,])
			data.frame(cellID=as.numeric(rownames(cbs)),rich=textTinyR:: sparse_Sums(cbs, rowSums = T))
		},mc.cores=mc.cores)
		rich.vec=do.call('rbind',richByCell)
		rich.r=env[[1]]
		values(rich.r)=NA
		values(rich.r)[rich.vec$cellID]= rich.vec$rich
		writeRaster(rich.r,file=paste0(outFileBase,'/', attrNames[y],'.tif'), overwrite=T)
		rich.r
	})
	stack(out)
}

#ongspeciesInPolygon=function(cbs,myPolygon)

