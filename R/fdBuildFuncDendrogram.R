
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

buildFuncDivTree=function(#GENERAL INPUTS
													sumDirs,
													speciesAttributeTable=NULL,
													colFuncTraits,#COLUMNS FUNCTIONAL TRAITS
													verbose=T,
													#FUNCTION DETAILS
													transformTraits=T,
													rulePCchoice='cumulativeProp', #how to choose PC? other options varianceProp
													rulePCthreshold= 0.8, #numeric, whatever number you like to apply the rule.
													writeOut=F) {

  #  forTesting
  # colFuncTraits=2:9; speciesAttributeTable= readRDS("./speciesAttributes.rds")
  #require(ape)
  #require(fastcluster)
  #require(tidyverse)


  #get the traits [the first line is just cuz I dont know if
  #we want to read the attr table as part of the upper function
  #also, one consolidated speciesAttrTable (?) or serveral? I think the first
  # CM: I think there could be multiple attribute tables so we should automatically read one in. it would be clean to just have 1 master table for all species attributes, but that could be too prescriptive for some people
  #what format for the table? I suggest one conslidated table and then user saves it as an RDSfile.
  # CM: yes agree; and i've now fixed to save as an rds

  #if(is.null(speciesAttributeTable)) speciesAttributeTable = paste0(sumDirs$attrDir, "/speciesAttributes.rds")
  if(is.character(speciesAttributeTable)) speciesAttributeTable = readRDS(file = speciesAttributeTable)

  #check is df-like object. We need that cuz the 99% we will have numeric and other things
  # CM: yes, it'll always be a df because it should have species names and indices.
  stopifnot(is.data.frame(speciesAttributeTable))

  #spTraits <- read.csv ("./TC_species_8Traits_mean_final.csv",header = T, sep = ';')
  species  = speciesAttributeTable[,1]
  spTraits = speciesAttributeTable[,colFuncTraits]

  #here is where I suggest B.Maitner could start splitting the PCAs or PCO depending on the data
  #a good expample in function FD::dbFD, the initial checks (e.g. half of the code)
  #are meant for guessing the type of data you have and apply the correct distance measure
  #most of the time this will be gowdis not euclidean.
  if (all (apply(spTraits,2,function (x)is.numeric(x))))
    (message('all traits are numeric.Show must go on!'))
  if (!all (apply(spTraits,2,function (x)is.numeric(x))))
    (stop('Sorry, right now this finicky function only takes numeric tratis with no NAs'))


  #log-transformed the traits and run PCA
  if (transformTraits) traits.scaled = spTraits %>% dplyr::mutate_if(is.numeric,log1p) %>% dplyr::mutate_if(is.numeric, scale)
  if (!transformTraits) traits.scaled= spTraits

  #subset to get the trait table
  traits.pc = princomp(traits.scaled)
  if (verbose) {message('Principal Components Results:'); print(summary(traits.pc))}

  #choice OF PC [underdevelopped, but fine for now]
  #also thought of as an indep function, but it will be so much dependent whteht PCA or PCO...not sure
  vars = apply(traits.pc$scores, 2, var)
  props = vars / sum(vars)
  cumprop = cumsum(props)
  if (rulePCchoice == 'cumulativeProp') {pcSelection = cumprop < rulePCthreshold}
  if (rulePCchoice == 'varianceProp') {pcSelection = props > rulePCthreshold}

  traits.pc.scores = as.data.frame(traits.pc$scores[, pcSelection])
  #traits.pc.spscores <- cbind(species, traits.pc.scores )

  rownames(traits.pc.scores) = species
  #traits.pc.spscores <- traits.pc.spscores[, -1]

  #attention BM this also could change depening on the data, but for all numeric this is right
  #however, this is only for a small amoutn of hte cases
  traits.pc.dist.mat = dist(traits.pc.scores, method = "euclidean")

  #quick solution, get FD using the phylo function
  traits.dendro = hclust(traits.pc.dist.mat, method = "average") #= UPGMA
  traits.phy <- ape::as.phylo(traits.dendro)
  traits.phy$root.edge <- 0

  if (writeOut==T) {
    #i dont know where you want to put it
    if (is.null(sumDirs$Misc)) dir.create (paste0(sumDirs$sumBaseDir,'/Misc'))
    ape::write.tree(traits.phy,file = paste0(sumDirs$Misc,'/functionaTraitDendrogram.tree'))
  }

  return (traits.phy)


  #### TO DEVELOP W BMAITNER
  # #Pep
  # #OKK so all above it is to compute a dendrogram of functional distance
  # # then pass it to an efficient packages of Phylogenetic Diversity
  # # the package FD does this an more, the thing is that some of this calculations take a long time
  # # and it should be tweaked to allow for BigData, PCComp or the convex hull part
  # # they recommend seting it to FALSE when computer crashes.
  #
  # THIS CODE IS A FIRST STAB AT IT...IN WHICH I TRIED TO RECHUNK EVERYTHING
  # library (FD)
  # source ("E:/Pep/TREECHANGE/SDM_ForDivFate/TCCC_WORLD_v6_10km__Summaries_sdmTP05rbg0165thr/workingSCripts/F_dbTRY.R")
  # str (traits.pc.dist.mat)
  # #FD is very phiniky so if the dist matrix is bigger than your species list, then you are fck
  # #When does this happen ?
  # #this happens if you used a matrix of all plants, but you do not have SDMs for all plants for which you have traits
  # #yet you want to depict your results in the plant functional space not your SDM subset
  # #I opted for recomputing the dist based on a subset
  #
  #
  # #eliminate all allZero Communities and species without NAs
  # commWithSp = which (rowSums(uniqueComm)!=0)
  # spOccurring = which (colSums(uniqueComm) != 0)
  # uniqueComm.subset= uniqueComm[commWithSp,spOccurring]
  #
  # #subset to common species
  # commonSp = dplyr::intersect(colnames(uniqueComm.subset),rownames (traits.pc.scores))
  # uniqueComm.subset = uniqueComm.subset[,commonSp]
  # traits.pc.spscores.subset = subset (traits.pc.scores,rownames(traits.pc.scores) %in% commonSp)
  # traits.pc.dist.mat.subset = dist(traits.pc.spscores.subset, method = "euclidean")
  #
  # #transform to ucmatrix
  # uc.matrix = as.matrix (uniqueComm.subset)
  #
  # div = dbFD.browse(x = traits.pc.dist.mat.subset,a = uc.matrix,calc.FRic = F)
  #

}
