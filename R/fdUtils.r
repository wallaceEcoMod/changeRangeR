


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
#' @title title
#'
# @description
# @param
# @param
# # @examples#@return
# @author Cory Merow <cory.merow@@gmail.com>
# @note
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
	fuck=cbs.f %>% basename %>% file_path_sans_ext %>% data.frame %>%
	  separate('.',c('c','ind')) %>% dplyr::select(ind)
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

cellIDToRaster=function(cell.ind,envGrid,colName){
	out=raster(envGrid[[1]])
	out[cell.ind$cellID]=cell.ind[,colName]
	out
}

#============================================================
#============================================================
#============================================================
#' @title Find corresponding cell IDs between two matrices
#'
#' @description See examples
#'
#' @details
#' See Examples.
#'
#' @param r1 a raster
#' @param r2 a raster
# @keywords
#' @export
#'
# @examples
#'
#' @author Cory Merow <cory.merow@@gmail.com>
# @note this is more reliable they're in the same projection. Otherwise, we just check which cell centers of r2 are an a cell of r1
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.
#' @export

# this only works if they're in the same projection!!
rasterDictionary <- function(r1, r2){ 
  # This finds the cell where the center of each cell is located. Is that guaranteed to be the cell where the majority of a cell is?  
  coords <- as.data.frame(raster::xyFromCell(r1, 1:ncell(r1))) 
  coords$cellID_r1 <- 1:ncell(r1) 
  # avoid huge slow operations when cells are NA
  coords=coords[complete.cases(values(r1)),]
  if(projection(r1)==projection(r2)){
    coords$cellID_r2 <- raster::cellFromXY(r2, as.matrix(coords[,1:2])) %>% as.integer
  } else {
  	coords.projr2=coords %>% select(x,y) %>% as.matrix %>% SpatialPoints(proj4string=CRS(projection(r1))) %>% spTransform(projection(r2))
  	coords$cellID_r2 <- raster::cellFromXY(r2, coords.projr2) %>% as.integer
  }
  return(coords[complete.cases(coords),])
} 

#============================================================
#============================================================
#============================================================
#' @title Can you put a sparse matrix in Memory
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
# @note based on code stolen from raster::canProcessInMemory
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.
#' @export

# > canProcessInMemory
# function (x, n = 4, verbose = FALSE) 
# {
#     if (.toDisk()) {
#         return(FALSE)
#     }
#     nc <- ncell(x)
#     n <- n * nlayers(x)
#     memneed <- nc * n * 8
#     maxmem <- .maxmemory()
#     memavail <- .availableRAM(maxmem)
#     if (verbose) {
#         gb <- 1073741824
#         cat("memory stats in GB")
#         cat(paste("\nmem available:", round(memavail/gb, 2)))
#         cat(paste0("\n        ", round(100 * .memfrac()), "%  : ", 
#             round(.memfrac() * memavail/gb, 2)))
#         cat(paste("\nmem needed   :", round(memneed/gb, 2)))
#         cat(paste("\nmax allowed  :", round(maxmem/gb, 2), " (if available)\n"))
#     }
#     if (nc > (2^31 - 1)) 
#         return(FALSE)
#     memavail <- .memfrac() * memavail
#     memavail <- min(memavail, maxmem)
#     if (memneed > memavail) {
#         options(rasterChunk = min(.chunksize(), memavail * 0.25))
#         return(FALSE)
#     }
#     else {
#         return(TRUE)
#     }
# }
