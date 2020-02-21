
#' @export
cm.cols.dif2=function(x,bias=1){ colorRampPalette(c('steelblue4','steelblue1','steelblue1', 'white', 'red1','red1','red4'))(x)
}

# general
cm.cols.dif3=function(x,bias=1){ colorRampPalette(c('navyblue','steelblue4','steelblue1', 'white', 'red1','red3','darkred'))(x)
}

#' @export
setupSummaryDirectories=function(summaryBaseDir){
	if(!file.exists(summaryBaseDir)) dir.create(summaryBaseDir)
	figDir=paste0(summaryBaseDir,'/Figures')
	if(!file.exists(figDir)) dir.create(figDir)
	attrDir=paste0(summaryBaseDir,'/AttributeTables')
	if(!file.exists(attrDir)) dir.create(attrDir)
	cbsDir=paste0(summaryBaseDir,'/CellSpeciesLists')
	if(!file.exists(cbsDir)) dir.create(cbsDir)
	rangeSizeDir=paste0(summaryBaseDir,'/RangeSizeDir')
	if(!file.exists(rangeSizeDir)) dir.create(rangeSizeDir)
	richDir=paste0(summaryBaseDir,'/SpeciesRichness/')
	if(!file.exists(richDir)) dir.create(richDir)

	return(list(sumBaseDir=summaryBaseDir, figDir=figDir,attrDir=attrDir,cbsDir=cbsDir, rangeSizeDir=rangeSizeDir,  richDir=richDir))
}

#' @export
getSpNamesFromDirs=function(x){basename(dirname(x))}


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
