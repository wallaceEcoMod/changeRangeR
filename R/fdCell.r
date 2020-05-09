
#===================================================================
#===================================================================
#===================================================================
#' @export
# cell Index Table
	# columns: long, lat, cellid
cellIndexTable=function(env,nCellChunks,sumDirs,toInteger=T){
	co=coordinates(env)
	#keep=apply(values(env),1,function(x) any(!is.na(x))) # not sure why i did this 
	keep=complete.cases(values(env))
	co=co[keep,]
	if(toInteger) co=apply(co,2,as.integer)
	# chunks can be based on regions (e.g. if cell values depend on one another) or just random 
	if(nCellChunks>1) { chunks=cut(1:nrow(co), nCellChunks, labels = FALSE) 
	} else { chunks=rep(1,nrow(co)) }
	cell.ind=data.frame(co,cellID=as.integer(cellFromXY(env,co)), chunkID=as.integer(chunks))
	write.csv(cell.ind,file=paste0(sumDirs$sumBaseDir, '/cellIndexTable.csv'),row.names=F)
	cell.ind
}

#===================================================================
#===================================================================
#===================================================================
#' @export
chunkFinder=function(someRaster,cell.ind,chunkFiles){
	notNAs=which(!is.na(values(someRaster)))
	keep=which(cell.ind$cellID %in% notNAs )
	chunks=unique(cell.ind[keep,'chunkID'])
	keep=mapply(function(x){grep(paste0('chunk_',x), basename(chunkFiles))},chunks)
	chunkFiles[keep]
}

#===================================================================
#===================================================================
#===================================================================
