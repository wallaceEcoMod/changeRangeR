#### SMALLER FUNCTION applied to
phyloDiv  <- function (uniqueCom,tree,fullMatch=T,verbose=T,...){

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

  PhyloMeasures::pd.query(tree = tree,matrix = uniqueCom)
}
