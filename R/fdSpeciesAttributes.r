
#' @title title
#' @description
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
# @examples
#
#' @return
#' @author Cory Merow <cory.merow@@gmail.com>
#' @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the docum
#' @export

#' @notes Currently just handles the mean, but will do other moments some day. This is only for species attributes that don't depend on the other species (e.g., good for traits, not nearest phylogenic neighbor)


speciesAttributeByCell=function(cbsDir,
													      attrTable,
													      scenario,
													      method='mean',
													      env,
																outDir=NULL,
																richnessRaster=NULL,
																verbose=F){

	t1=proc.time()
	message(paste0('starting ',scenario))
	if(Sys.info()["sysname"]== "Windows") mclapply=parallelsugar::mclapply
	if(Sys.info()["sysname"]!= "Windows") mclapply=parallel::mclapply
	cbs.f=changeRangeR:::.getCBS(cbsDir,scenario)

	attrNames=names(attrTable)
	# if species and index were included, toss them
	toss=unlist(mapply(function(x){grep(x,attrNames)}, c('species','index')))
	if(length(toss) > 0 ) attrNames=attrNames[-toss]

	out=lapply(seq_along(attrNames),function(y){
		if(verbose) message(attrNames[y])
		#keep=attrTable$index[attrTable[attrNames[y]]==1]
		if(method=='mean'){
			attrByCell=mclapply(seq_along(cbs.f), function(x){
				if(verbose) message(x)
				cbs=readRDS(cbs.f[x])
				b=cbs %*% as.matrix(attrTable[attrNames[y]],ncol=1)
				rich=textTinyR::sparse_Sums(cbs, rowSums = T)
				data.frame(cellID=as.numeric(rownames(cbs)),thisAttr=as.numeric(b/rich))
			},mc.cores=mc.cores)
		} else{
			stop('sorry, only mean values supported at this point')
				# Var(X) = Σ ( Xi - X_mean )2 / N
		}

		attr.vec=do.call('rbind',attrByCell)
		attr.r=raster(env[[1]])
		values(attr.r)[attr.vec$cellID]= attr.vec$thisAttr
		if(!is.null(outDir)) writeRaster(attr.r,file=paste0(outDir,'/',attrNames[y],'_', scenario,'.tif'), overwrite=T)
		attr.r
	})
	out1=stack(out)
	names(out1)=attrNames

	t2=proc.time()-t1
	message(paste0(round(t2[3],2),' s'))
  out1
}
