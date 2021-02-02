
#===================================================================
#===================================================================
#===================================================================
#' Create a table linking cell IDs with their longitude and latitude which can optoionally be used as a template for a	`cellAttributeTable`.
#' @notes Columns are long, lat, cellid
#' @export
cellIndexTable=function(env,nCellChunks,sumDirs=NULL,toInteger=T){
	co=coordinates(env)
	#keep=apply(values(env),1,function(x) any(!is.na(x))) # not sure why i did this
	keep=complete.cases(values(env))
	co=co[keep,]
	if(toInteger) co=apply(co,2,as.integer)
	# chunks can be based on regions (e.g. if cell values depend on one another) or just random
	if(nCellChunks>1) { chunks=cut(1:nrow(co), nCellChunks, labels = FALSE)
	} else { chunks=rep(1,nrow(co)) }
	cell.ind=data.frame(co,cellID=as.integer(cellFromXY(env,co)), chunkID=as.integer(chunks))
	cell.ind=sapply(cell.ind,as.integer) %>% data.frame
	if(!is.null(sumDirs))	saveRDS(cell.ind,file=paste0(sumDirs$sumBaseDir, '/cellIndexTable.rds'))
	cell.ind
}

#===================================================================
#===================================================================
#===================================================================
# used internally to read in only the necessary chunks for a given analysis
# This seems to assume that `someRaster` is already exactly matched to the env grid used for cell.ind
#' @export
.chunkFinder=function(someRaster,cell.ind,chunkFiles){
	notNAs=which(!is.na(values(someRaster)))
	keep=which(cell.ind$cellID %in% notNAs )
	chunks=unique(cell.ind[keep,'chunkID'])
	keep=mapply(function(x){grep(paste0('chunk_',x), basename(chunkFiles))},chunks)
	chunkFiles[keep]
}

