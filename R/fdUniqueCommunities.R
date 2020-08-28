#' @title
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


#' @param cbsDir directory where the CBS.rds are stored
#' @param mc.cores allowed cores for parallelization, default = 1
#' @param
#' @note it will write the uniqueCommunities for each chunk in the cbsDir


makeUniqueCommunities <- function (cbsDir,scenario,mc.cores,overwrite=F){
  #list cbs rds files
  message(scenario)
  t1=proc.time()
	cbs.f=changeRangeR:::.getCBS(cbsDir,scenario)
  if (length(cbs.f)==0) stop ('No cellBySpecies.RDS files found, cannot compute unique communities')

  #check corresponding cbs unique communities
  uc.f=list.files(paste0(cbsDir,'/',scenario),full.names=T,pattern = 'uniqueComm')

  #check whether they already exist
  if(length(uc.f) == length(cbs.f) & !overwrite) {
    message('Unique communities already computed. Choose overwrite=T to remake.')
    return (uc.f)
  }
  #check corresponding cbs and unique communities match if not compute them
  #if(length(uc.f) == 0 | (length(uc.f) != length(cbs.f))){
    # make it Windows compatible :/
	if (Sys.info()["sysname"]== "Windows") {mclapply <- parallelsugar::mclapply}
	if (Sys.info()["sysname"]!= "Windows") {mclapply <- parallel::mclapply}

	mclapply(seq_along(cbs.f), function(x){
		message(x)
		cbs=readRDS(cbs.f[x])
		uc.cbs = .duplicated.dgCMatrix(cbs,MARGIN = 1,include.all.zero.vectors = T)
		uc.cbs$cellind = as.integer(rownames(cbs))
		names (uc.cbs)=c('duplicated','comChunkID','cellind')
		uc.filename = basename(cbs.f[x])
		uc.filename = paste0(cbsDir,'/',scenario,'/uniqueComm',strsplit (uc.filename,split = 'chunk')[[1]][2])
		saveRDS (uc.cbs,uc.filename)
	},mc.cores =mc.cores)
  #}

  uc.f=list.files(paste0(cbsDir,'/',scenario),full.names=T,pattern = 'uniqueComm_')
  t2=proc.time()-t1
  message(paste0(round(t2[3],2),' s'))
  uc.f
}


.duplicated.dgCMatrix <- function (dgCMat, MARGIN, include.all.zero.vectors = TRUE) {
  MARGIN <- as.integer(MARGIN)
  J <- rep(1:ncol(dgCMat), diff(dgCMat@p))
  I <- dgCMat@i + 1
  x <- dgCMat@x
  if (MARGIN == 1L) {
    ## check duplicated rows
    names(x) <- J
    if (include.all.zero.vectors) {
      RowLst <- split(x, factor(I, levels = 1:nrow(dgCMat)))
    } else {
      RowLst <- split(x, I)  ## will do `factor(I)` internally in `split`
    }
    dupsLogical <- duplicated.default(RowLst)
    uniques <- unique (RowLst)
    dupsIdx <- match (RowLst,uniques)
    list (dupsLogical,dupsIdx)


  } else if (MARGIN == 2L) {
    ## check duplicated columns
    names(x) <- I
    if (include.all.zero.vectors) {
      ColLst <- split(x, factor(J, levels = 1:ncol(dgCMat)))
    } else {
      ColLst <- split(x, J)  ## will do `factor(J)` internally in `split`
    }
    dupsLogical <- duplicated.default(RowLst)
    uniques <- unique (RowLst)
    dupsIdx <- match (RowLst,uniques)
    list (dupsLogical,dupsIdx)

  } else {
    stop("invalid MARGIN; return NULL")
    return(NULL)
  }

}
