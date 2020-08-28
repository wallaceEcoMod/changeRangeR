


#============================================================
#============================================================
#============================================================

#============================================================
#============================================================
#============================================================

#' Color palletes
#' Useful contrasts for plotting

#' @export
cm.cols1=function(x,bias=1) { colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x)
}

#' @export
cm.cols.dif3=function(x,bias=1){ colorRampPalette(c('navyblue','steelblue4','steelblue1', 'white', 'red1','red3','darkred'))(x)
}

#' @export
cm.cols.dif2=function(x,bias=1){ colorRampPalette(c('steelblue4','steelblue1','steelblue1', 'white', 'red1','red1','red4'))(x)
}


#============================================================
#============================================================
#============================================================
<<<<<<< HEAD
#' @title
=======
#' @title title
>>>>>>> f844c8a487e5ed5b4b3c878a5ea65be608980558
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
getSpNamesFromDirs=function(x){basename(dirname(x))}

#============================================================
#============================================================
#============================================================
#' @title Split up jobs for parallelization
#'
#' @description See examples
#'
#' @details
#' See Examples.
#'
#' @param x a vector or dataframe with 1 column
#' @param ntasks number of cores to parallelize over
# @keywords
#' @export
#'
#' @examples
#' x=1:100
#' chopTasks(x,7)
#' @return Returns a list with \code{ntasks} elements, useful for sending to a \code{foreach} loop
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.
#' @export

.chopTasks=function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))



#============================================================
#============================================================
#============================================================
#' @title Aggregate cells and connect cell IDs at different scales
#'
#' @description See examples
#'
#' @details
#' See Examples.
#'
#' @param x a vector or dataframe with 1 column
#' @param ntasks number of cores to parallelize over
# @keywords
#' @export
#'
#' @examples
#' x=1:100
#' chopTasks(x,7)
#' @return Returns a list with \code{ntasks} elements, useful for sending to a \code{foreach} loop
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.
#' @export

aggregateCells=function(cell.ind,facts,mc.cores=mc.cores,outDir){
	#  for testing
	#  outDir=sumDirs$envDir

	coordinates(cell.ind)=c('x','y')
	if (Sys.info()["sysname"]== "Windows") {mclapply <- parallelsugar::mclapply}
	if (Sys.info()["sysname"]!= "Windows") {mclapply <- parallel::mclapply}
	coarser=mclapply(facts,function(x){
		coarseEnv=raster::aggregate(env[[1]],fact=x)
		if(!is.null(outDir)) writeRaster(coarseEnv,file=paste0(outDir,'/env_Agg',x,'.tif'),datatype = "INT2S",overwrite=T)
		cellFromXY(coarseEnv,cell.ind)
	},mc.cores=mc.cores)
	# merge with cell.ind
	c1=do.call(cbind,coarser) %>% data.frame
	names(c1)=paste0('Agg',facts)
	cellPerChunk=cell.ind@data %>% dplyr::filter(chunkID==1) %>% nrow
	tmp=lapply(1:ncol(c1),function(x){
		coarseCellID=unique(c1[,x])
		nchunks=ceiling(length(coarseCellID)/cellPerChunk)
		ch=data.frame(coarseCellID, rep(1:nchunks,each=cellPerChunk)[1:length(coarseCellID)])
		names(ch)=c(names(c1)[x],paste0('chunkID_',names(c1)[x]))
		c1 %>% dplyr::left_join(ch) %>% dplyr::select_if(grepl('chunk',names(.)))
	})
	tmp1=do.call('cbind',tmp)
	out=data.frame(cell.ind,c1,tmp1) %>% dplyr::select(-optional)
	out=sapply(out,as.integer)
}

#============================================================
#============================================================
#============================================================
#' @note gets them in the right order
#' @export
.getCBS=function(cbsDir,scenario){
	cbs.f=list.files(paste0(cbsDir,'/',scenario),full.names=T,pattern='chunk')
<<<<<<< HEAD
	fuck=cbs.f %>% basename %>% file_path_sans_ext %>% data.frame %>%
	  separate('.',c('c','ind')) %>% dplyr::select(ind)
=======
	fuck=cbs.f %>% basename %>% file_path_sans_ext %>% data.frame %>% separate('.',c('c','ind')) %>% select(ind)
>>>>>>> f844c8a487e5ed5b4b3c878a5ea65be608980558
	ord=order(as.numeric(fuck$ind))
	cbs.f[ord]
}


#============================================================
#============================================================
#============================================================
#' @note gets them in the right order
#' @export
.getUCom=function(cbsDir,scenario){
  uc.f=list.files(paste0(cbsDir,'/',scenario),full.names=T,pattern='unique')
  fuck=uc.f %>% basename %>% file_path_sans_ext %>% data.frame %>%
    separate('.',c('c','ind')) %>% dplyr::select(ind)
  ord=order(as.numeric(fuck$ind))
  uc.f[ord]
}

#============================================================
#============================================================
#============================================================
#' @title Put sparse matrix values into raster
#'
#' @description See examples
#'
#' @details
#' See Examples.
#'
#' @param
#' @param
# @keywords
#' @export
#'
# @examples
#'
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.
#' @export

sparseToRaster=function(cell.ind,envGrid,colName){
	out=raster(envGrid[[1]])
	out[cell.ind$cellID]=cell.ind[,colName]
	out
}
