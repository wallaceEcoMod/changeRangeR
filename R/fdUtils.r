
#============================================================
#============================================================
#============================================================
# general
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

#' @export
setupSummaryDirectories=function(summaryBaseDir){
	if(!file.exists(summaryBaseDir)) dir.create(summaryBaseDir)
	figDir=paste0(summaryBaseDir,'/Figures')
	if(!file.exists(figDir)) dir.create(figDir)
	attrDir=paste0(summaryBaseDir,'/AttributeTables')
	if(!file.exists(attrDir)) dir.create(attrDir)
	cbsDir=paste0(summaryBaseDir,'/CellSpeciesLists')
	if(!file.exists(cbsDir)) dir.create(cbsDir)
	rangeSizeDir=paste0(summaryBaseDir,'/RangeSize')
	if(!file.exists(rangeSizeDir)) dir.create(rangeSizeDir)
	richDir=paste0(summaryBaseDir,'/SpeciesRichness/')
	if(!file.exists(richDir)) dir.create(richDir)

	return(list(sumBaseDir=summaryBaseDir, figDir=figDir,attrDir=attrDir,cbsDir=cbsDir, rangeSizeDir=rangeSizeDir,  richDir=richDir))
}
#============================================================
#============================================================
#============================================================

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

chopTasks=function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))


#============================================================
#============================================================
#============================================================

duplicated.dgCMatrix <- function (dgCMat, MARGIN, include.all.zero.vectors = TRUE) {
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
