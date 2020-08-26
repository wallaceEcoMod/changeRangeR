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



#' @param FUN function defining the metric to be calculated.
#' @param outputFUNnames desired output name of the metric
#' @param ... other arguments to be passed to FUN
#' @note most functions `FUN` will need to have the colnames to be the name of the species

metricFromCBS=function(cbsDir,
											 scenario,
											 env,
											 spIndTable,
											 cellAttributeTable=NULL, speciesAttributeTable=NULL,# I think we need to add this when wanting to subset or give weights
											 mc.cores=1,
											 #specify the function that takes dense regular matrix of locationBySpecies for computation
											 #bear in mind that the first argument of this function should be the matrix
											 FUN=phyloDiv,
											 outputFUNnames ,
											 ...){

  t1=proc.time()
  message(scenario)

	cbs.f=changeRangeR:::.getCBS(cbsDir,scenario)
	if (Sys.info()["sysname"]== "Windows") {mclapply <- parallelsugar::mclapply}
  if (Sys.info()["sysname"]!= "Windows") {mclapply <- parallel::mclapply}

  #Pep
  #check if unique communities files exists otherwise compute them.
  #CM: i think we have users explicitly make these, as in the demo
  #uc.f =makeUniqueCommunities(cbsDir =paste0(cbsDir,'/',scenario),mc.cores = mc.cores )

  #create list of pairs of uc and cbs
  inputCBSlist = lapply(seq_along(cbs.f), function (x) c(cbs.f[x],uc.f[x]))

  #function to apply over a pairs of cbs.f and uc.f a given function
  #function that converts inputs to unique communities dense matrices and pass them to the target function.
  #if sparse -> dense throws a memory issue, it computes over blocks that converts inputs to unique communities dense matrices and passs them to the target function.

  #the aim of this function is to extract unique communities, pass them into a dense matrix and then compute the metricFunction FUN
  applyFUNtoDenseMatrixComm = function(x,spind=spIndTable$species, ...){ #uniqueCom = uc.cbs.matrix,fullMatch = F
    message(basename(x[1]))
    cbs=readRDS(x[1])
    ucList=readRDS(x[2])
    colnames(cbs)=spind
    #cellBySp Matrix for unique Communites
    uc.cbs = cbs[!ucList$duplicated,]
    #somewhere here be able to subset by cellAttr or spAttr

    #convert to dense matrix ()if you can
# 		#=== End CM Edit ====================================
#     # CM: lets guess the matrix size and skip this if its likely too big
#     mat.size=round(nrow(uc.cbs)) * ncol(uc.cbs) * 8e-9
#     if(mat.size < memSizeGb){
#     	uc.cbs.matrix = try ({as.matrix (uc.cbs)},silent=T)
# 			gc()
#     	#if memory holds, go directly to function
#     	#if (!class(uc.cbs.matrix)=='try-error')
#     	outMetric = try(FUN(uc.cbs.matrix,...))
#     	#outMetric = try(FUN(uc.cbs.matrix,fullMatch=F,tree=phyloTree)) # for testing
#     } else {uc.cbs.matrix=outMetric=NULL}
# 		#=== End CM Edit ====================================
#
# 		# CM: also getting memory errors when running FUN even after the matrix works
#     #fuck loop across if too much memory (this should be lapply not a for)
#     if (class(uc.cbs.matrix)=='try-error' | class(outMetric)=='try-error' | mat.size >= memSizeGb){

      # THIS IS A CODE CHUNK TO BE ACTIVATED if uc.cbs.matrix is too big. Nut sure it does the jub
      #bc w/ lapply it keeps the chunks in memory while computing, I think
      # CM: this update to mclapply seems to work fine with memory
      levelsRows = seq (from=1, to = nrow (uc.cbs), by = 1000)
      if(length(levelsRows)==1){
      	rowGroups=rep(1,nrow(uc.cbs))
      } else {
	      rowGroups = cut(1:nrow(uc.cbs), length(levelsRows), labels = FALSE)
	    }
      #outMetric = vector('list',length(levelsRows))
      #for (g in seq_along(levelsRows)) {
      outMetric=mclapply(seq_along(levelsRows),function(g){
        message(paste0('Computing group : ',g,'/',max(rowGroups)))
        uniqueComSubset = uc.cbs[which(rowGroups == g),]
        uc.cbs.matrix = as.matrix(uniqueComSubset)
        #out = FUN(uc.cbs.matrix,...) # comment for testing
        out = FUN(uc.cbs.matrix,fullMatch=F,tree=phyloTree) # uncomment for testing
        gc()
        out
      },mc.cores=mc.cores)
      #if (is.vector(outMetric[[1]])) outMetric = do.call (c, outMetric)
      #if (!is.vector(outMetric[[1]])) outMetric = do.call (cbind, outMetric)
			#CM:  previous cases won't happen
			outMetric2 = do.call('c', outMetric)
			gc()
    #}

    #reconstitute cellID from community ID
    #  CM: i think this will be the only case now, right?
    uc.id = tibble::tibble(comChunkID =ucList$comChunkID[!ucList$duplicated], outMetric2)
#     if (is.vector (outMetric2)){ uc.id = tibble::tibble(comChunkID =ucList$comChunkID[!ucList$duplicated], outMetric2)}
#     if (!is.vector (outMetric)) {
#       uc.id = tibble::as_tibble(outMetric)
#       uc.id$comChunkID =ucList$comChunkID [!ucList$duplicated]
#     }

    uc.cellid = tibble::tibble (comChunkID= ucList$comChunkID,cellind = ucList$cellind)

    outmetric.cellid = dplyr::left_join(uc.cellid,uc.id,by = "comChunkID")
    outmetric.cellid
  }

  kk =lapply(seq_along(inputCBSlist), function (x) {
  	#out=applyFUNtoDenseMatrixComm(inputCBSlist[[x]],...)
  	out=applyFUNtoDenseMatrixComm(inputCBSlist[[x]],fullMatch=F,tree=phyloTree) # for testing
  	gc()
  	out
  })
  pdByCell=do.call('rbind',kk)
  #this is not elegant
  nvar = ncol (pdByCell)-2
  listRas = lapply(3:ncol(pdByCell), function(outmetricID){
    pd.r = raster(env[[1]])
    pd.r[as.numeric(pdByCell$cellind)]= pdByCell%>% dplyr::pull(outmetricID)
    pd.r
  })
  ras = do.call(stack,listRas)
  names (ras) = outputFUNnames

  t2=proc.time()-t1
	message(paste0(round(t2[3],2),' s'))
  return(ras)


}
