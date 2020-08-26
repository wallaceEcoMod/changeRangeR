
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

#### build species based maps
spRangeChangeFromCBS=function(sumDirs,
															scn1,
															scn2,
															env,
															mc.cores){
  #for testing
  # scn1='present'
  # scn2='8580'
  # mc.cores=15

  cbsDir=sumDirs$cbsDir

  cbsList = lapply (list(scn1,scn2), function (x){
    cbs.f=list.files(paste0(cbsDir,'/',x),pattern = "chunk_",full.names=T)
    toss=grep('temp_|spCellOcc',cbs.f)
    if(length(toss>0)) cbs.f=cbs.f[-toss]
    cbs.f
  })

  #check they have the same number of cbs
  L <- lapply (cbsList, function (x) length(x))
  if (!all(sapply(L, identical, L[[1]]))) stop ('scenarios should have the sam nubmer of chunks')

  #get number of chunks
  nchunks = length(cbsList[[1]])
  cl = parallel::makeCluster(min (mc.cores,nchunks))
  doSNOW::registerDoSNOW(cl)
  iterations = length(nchunks)
  pb = txtProgressBar(max = iterations, style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)

  require(doParallel)
  srcBySP=foreach::foreach(x= 1:nchunks,
                             .options.snow=opts,
                             .packages = c((.packages()))
                              )%dopar%{

                               message(x)

                               cbsScenarios = lapply (cbsList, function (l,i=x){
                                 readRDS(l[i])

                               })

                               cbsSpDyn = cbsScenarios[[2]]-cbsScenarios[[1]]

                               Lost <-  Matrix::colSums(cbsSpDyn == -1,sparseResult = T)
                               outL = Matrix(data = Lost,nrow=length(Lost),ncol=1)

                               Gained <-  Matrix::colSums (cbsSpDyn == 1,sparseResult = T)
                               outG = Matrix(data = Gained,nrow=length(Gained),ncol=1)

                               Maintained  <-  Matrix::colSums((cbsScenarios[[2]]+cbsScenarios[[1]]) == 2, sparseResult = T)
                               outM = Matrix(data = Maintained,nrow=length(Maintained),ncol=1)

                               outTotal = cbind(outL,outG,outM)
                               colnames (outTotal)=c('Lost','Gained','Maintained')
                               rownames (outTotal)= colnames (cbsSpDyn)

                               outTotal

                             }
  registerDoSEQ()
  gc()

  srcBySP = Reduce(`+`,srcBySP)

  #load the species index table
  spIndTable = data.table::fread(paste0(sumDirs$sumBaseDir,'/speciesIndexTable.csv'))

  srcBySP= as.matrix(srcBySP)
  srcBySP=as_tibble(srcBySP)
  srcBySP$spID= spIndTable$species
  srcBySP


}


