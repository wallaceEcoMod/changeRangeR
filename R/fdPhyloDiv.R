
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

#### SMALLER FUNCTION applied to
phyloDiv  <- function (uniqueCom,tree,fullMatch=T,verbose=T,...){
  # for testing
  #  uniqueCom; tree; fullMatch=T; verbose=T
  #require(ape)
  #require(PhyloMeasures)
  #check that tips of the tree and species name match
  percMatch = (sum (colnames(uniqueCom)%in%tree$tip.label) / length(colnames(uniqueCom))) *100

  if (fullMatch & percMatch!=100){stop('Some species not in the phylogenetic tree')}
  if (!fullMatch & percMatch < 100){
    if (verbose)    warning( paste('Some species not in the phylogenetic tree. Check output. Percentage match =',percMatch))
    selecSp = which (colnames(uniqueCom)%in%tree$tip.label)
    spNotIn = colnames(uniqueCom) [selecSp * (-1)]
    uniqueCom = uniqueCom[,selecSp]
  }

  #just in case the inpute has been passed as a sparse matrix
  if (!is.matrix (uniqueCom)) uniqueCom = as.matrix (uniqueCom)


  suppressMessages(base::suppressWarnings({outPhyloMetric = PhyloMeasures::pd.query(tree = tree,matrix = uniqueCom)}))


  outPhyloMetric

}
