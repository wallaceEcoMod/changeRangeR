#' @title cbs scenario filtering
#' @description modify an output scenario based on binary species rasters to get a new scenario that supports a condition on species and cell charactristics
#' @param cbsDir directory of cell by species matrix
#' @param scenario name of the target scenario to modify
#' @param suffixScenario suffix to be added to scenario for the output (the modified cbs)
#' @param spAttrTable species attribute table
#' @param spAttr character. The name of the field with the species attribute of interest
#' @param cellAttrTable site (cell) attribute table
#' @param cellAttr character. The name of the field with the species attribute of interest
#' @param condition character. Logical condition to be used to compare spAttr to cellAttr
#' @param NAfill numeric or NULL. Should NA in cellAttr or spAttr be filled with a number
#' @param mc.cores numeric. Number of cores to be used.
#' @param g numeric. Number of species chunks to process the data. Higher numbers imply less RAM utilization. Default to 10.

# @examples
#
#' @return saves a cbs file
#' @author Cory Merow <cory.merow@@gmail.com>, Pep Serra
#' @note Mask cbs by a condition involving species characteristics and cell characteristcis.
#'Example when you have a species that can't live near cities: input a species with a city tolorance and a table with cell urbanized
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the docum
#' @export

fdCBSfilterBySpByCell=function(cbsDir,
															 scenario,
															 suffixScenario,
															 spAttrTable,
															 spAttr,
															 cellAttrTable ,
															 cellAttr ,
															 condition,
															 NAfill = NULL,
															 mc.cores=1,
															 useForeach=F,
															 g = 10,verbose=T ){

  #  for testing
  # cbsDir = sumDirs$cbsDir,
  # scenario = '4580_200myr',
  # suffixScenario= 'LUIfilter2',
  # spAttrTable = spLUI,
  # spAttr='intQ90',
  # cellAttrTable = cellLUI,
  # cellAttr ='intensity',
  # condition = '>=',
  # NAfill = 0,
  # mc.cores=15,
  # g = 10

  t1=proc.time()
  message(paste0('starting ',scenario))
  cbs.f=changeRangeR:::.getCBS(cbsDir,scenario)
  if(Sys.info()["sysname"]== "Windows") mclapply=parallelsugar::mclapply
  if(Sys.info()["sysname"]!= "Windows") mclapply=parallel::mclapply

  #tibblize
  cellAttrTable <- tibble::as_tibble(cellAttrTable)
  spAttrTable <- tibble::as_tibble(spAttrTable)

  #build output folder
  outDir = paste0(cbsDir,'/',scenario,'_',suffixScenario)
  dir.create (outDir)

  ### with foreach
  #get number of chunks
  nchunks = length(cbs.f)

	if(!useForeach){
		cbs.new=mclapply(1:nchunks,function (x){
			if(verbose) message(x)
			#read in cbs
			cbs=readRDS(cbs.f[x])
			#get cell ID and the attribute values of cell interest
			chunkCellID = as.numeric(rownames(cbs))
			chunkCellAttr = cellAttrTable [which (cellAttrTable$cellind %in%chunkCellID),]
			if (is.numeric(NAfill)) {chunkCellAttr[,cellAttr][is.na(chunkCellAttr[,cellAttr])] <- NAfill }
			cellValue = chunkCellAttr[,cellAttr] %>% pull
			#get species attributes
			if (is.numeric(NAfill)) {spAttrTable [,spAttr][is.na(spAttrTable [,spAttr])] <- NAfill }
			spValue   = spAttrTable [,spAttr] %>%pull

			#prepare cbs out
			cbs.mod = Matrix (data = c(0),nrow = nrow(cbs),ncol=ncol(cbs))

			#create species chunks (default is 10)
			#g = 10
			byArg = round (ncol (cbs.mod)/g,digits = 0)
			a =seq(from = 1,to = ncol (cbs.mod), by= byArg)
			b = a+ (byArg)-1
			b[length(b)] <- ncol(cbs.mod)

			#loop by species chunks
			for (colChunk in 1:g) {
				tictoc::tic()
				spInfoMatrixChunk = Matrix(data = spValue[a[colChunk]:b[colChunk]],
																	 ncol = length(a[colChunk]:b[colChunk]),
																	 nrow = nrow(cbs.mod),
																	 byrow = T)

				cellInfoMatrixChunk = Matrix(data = cellValue,
																		 ncol = length(a[colChunk]:b[colChunk]),
																		 nrow = nrow(cbs.mod),
																		 byrow = F)

				cbs.Chunk = cbs[,a[colChunk]:b[colChunk]]
				binaryMatrixCond = eval(expr = parse(text = paste('spInfoMatrixChunk' ,condition ,'cellInfoMatrixChunk')))

				out = cbs.Chunk * binaryMatrixCond
				cbs.mod [,a[colChunk]:b[colChunk]] <- out
				tictoc::toc()

			}

			#write output
			dimnames (cbs.mod)<- dimnames(cbs)
			saveRDS(object = cbs.mod, file = paste0(outDir,'/',basename(cbs.f[x])))

		})
	} 
	#---------------------------------------------------------------
	# if your life sucks and you use a PC
	if(useForeach){
	
	 #prepare parallel
  cl = parallel::makeCluster(mc.cores)
  doSNOW::registerDoSNOW(cl)
  iterations = length(nchunks)
  pb = txtProgressBar(max = iterations, style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  library(doParallel)

  cbs.new =foreach::foreach(x = 1:nchunks,
		.options.snow=opts,
		.packages = c((.packages()))) %dopar% {
			if(verbose) message(x)
			#read in cbs
			cbs=readRDS(cbs.f[x])
			#get cell ID and the attribute values of cell interest
			chunkCellID = as.numeric(rownames(cbs))
			chunkCellAttr = cellAttrTable [which (cellAttrTable$cellind %in%chunkCellID),]
			if (is.numeric(NAfill)) {chunkCellAttr[,cellAttr][is.na(chunkCellAttr[,cellAttr])] <- NAfill }
			cellValue = chunkCellAttr[,cellAttr] %>% pull
			#get species attributes
			if (is.numeric(NAfill)) {spAttrTable [,spAttr][is.na(spAttrTable [,spAttr])] <- NAfill }
			spValue   = spAttrTable [,spAttr] %>%pull

			#prepare cbs out
			cbs.mod = Matrix (data = c(0),nrow = nrow(cbs),ncol=ncol(cbs))

			#create species chunks (default is 10)
			#g = 10
			byArg = round (ncol (cbs.mod)/g,digits = 0)
			a =seq(from = 1,to = ncol (cbs.mod), by= byArg)
			b = a+ (byArg)-1
			b[length(b)] <- ncol(cbs.mod)

			#loop by species chunks
			for (colChunk in 1:g) {
				tictoc::tic()
				spInfoMatrixChunk = Matrix(data = spValue[a[colChunk]:b[colChunk]],
																	 ncol = length(a[colChunk]:b[colChunk]),
																	 nrow = nrow(cbs.mod),
																	 byrow = T)

				cellInfoMatrixChunk = Matrix(data = cellValue,
																		 ncol = length(a[colChunk]:b[colChunk]),
																		 nrow = nrow(cbs.mod),
																		 byrow = F)

				cbs.Chunk = cbs[,a[colChunk]:b[colChunk]]
				binaryMatrixCond = eval(expr = parse(text = paste('spInfoMatrixChunk' ,condition ,'cellInfoMatrixChunk')))

				out = cbs.Chunk * binaryMatrixCond
				cbs.mod [,a[colChunk]:b[colChunk]] <- out
				tictoc::toc()

			}

			#write output
			dimnames (cbs.mod)<- dimnames(cbs)
			saveRDS(object = cbs.mod, file = paste0(outDir,'/',basename(cbs.f[x])))

		}

  	registerDoSEQ()
  	stopCluster(cl)
	}

  message ('done!')

}


#=========================================================================
# CM: commented out because its now an option in the main function. This can be deleted by PEP whenever he's sure its not needed.
# # 
# # #' @title cbs scenario filtering
# # #' @description modify an output scenario based on binary species rasters to get a new scenario that supports a condition on species and cell charactristics
# # #' @param cbsDir directory of cell by species matrix
# # #' @param scenario name of the target scenario to modify
# # #' @param suffixScenario suffix to be added to scenario for the output (the modified cbs)
# # #' @param spAttrTable species attribute table
# # #' @param spAttr character. The name of the field with the species attribute of interest
# # #' @param cellAttrTable site (cell) attribute table
# # #' @param cellAttr character. The name of the field with the species attribute of interest
# # #' @param condition character. Logical condition to be used to compare spAttr to cellAttr
# # #' @param NAfill numeric or NULL. Should NA in cellAttr or spAttr be filled with a number
# # #' @param mc.cores numeric. Number of cores to be used.
# # #' @param g numeric. Number of species chunks to process the data. Higher numbers imply less RAM utilization. Default to 10.
# # 
# # # @examples
# # #
# # #' @return saves a cbs file
# # #' @author Cory Merow <cory.merow@@gmail.com>, Pep Serra
# # #' @note Mask cbs by a condition involving species characteristics and cell characteristcis.
# # #'Example when you have a species that can't live near cities: input a species with a city tolorance and a table with cell urbanized
# # # @seealso
# # # @references
# # # @aliases - a list of additional topic names that will be mapped to
# # # this documentation when the user looks them up from the command
# # # line.
# # # @family - a family name. All functions that have the same family tag will be linked in the docum
# # #' @export
# # 
# # .fdCBSfilterBySpByCell_foreach=function(cbsDir,
# # 																				 scenario,
# # 																				 suffixScenario,
# # 																				 spAttrTable,
# # 																				 spAttr,
# # 																				 cellAttrTable ,
# # 																				 cellAttr ,
# # 																				 condition,
# # 																				 NAfill = NULL,
# # 																				 mc.cores=1,
# # 																				 g = 10,verbose=T ){
# # 
# #   #for testing
# #   # cbsDir = sumDirs$cbsDir,
# #   # scenario = '4580_200myr',
# #   # suffixScenario= 'LUIfilter2',
# #   # spAttrTable = spLUI,
# #   # spAttr='intQ90',
# #   # cellAttrTable = cellLUI,
# #   # cellAttr ='intensity',
# #   # condition = '>=',
# #   # NAfill = 0,
# #   # mc.cores=15,
# #   # g = 10
# # 
# #   t1=proc.time()
# #   message(paste0('starting ',scenario))
# #   cbs.f=changeRangeR:::.getCBS(cbsDir,scenario)
# #   if(Sys.info()["sysname"]== "Windows") mclapply=parallelsugar::mclapply
# #   if(Sys.info()["sysname"]!= "Windows") mclapply=parallel::mclapply
# # 
# #   #tibblize
# #   cellAttrTable <- tibble::as_tibble(cellAttrTable)
# #   spAttrTable <- tibble::as_tibble(spAttrTable)
# # 
# #   #build output folder
# #   outDir = paste0(cbsDir,'/',scenario,'_',suffixScenario)
# #   dir.create (outDir)
# # 
# #   ### with foreach
# #   #get number of chunks
# #   nchunks = length(cbs.f)
# # 
# #   #prepare parallel
# #   cl = parallel::makeCluster(mc.cores)
# #   doSNOW::registerDoSNOW(cl)
# #   iterations = length(nchunks)
# #   pb = txtProgressBar(max = iterations, style = 3)
# #   progress = function(n) setTxtProgressBar(pb, n)
# #   opts = list(progress = progress)
# #   library(doParallel)
# # 
# #   cbs.new =foreach::foreach(x = 1:nchunks,
# #                             .options.snow=opts,
# #                             .packages = c((.packages()))) %dopar% {
# #                               if(verbose) message(x)
# #                               #read in cbs
# #                               cbs=readRDS(cbs.f[x])
# #                               #get cell ID and the attribute values of cell interest
# #                               chunkCellID = as.numeric(rownames(cbs))
# #                               chunkCellAttr = cellAttrTable [which (cellAttrTable$cellind %in%chunkCellID),]
# #                               if (is.numeric(NAfill)) {chunkCellAttr[,cellAttr][is.na(chunkCellAttr[,cellAttr])] <- NAfill }
# #                               cellValue = chunkCellAttr[,cellAttr] %>% pull
# #                               #get species attributes
# #                               if (is.numeric(NAfill)) {spAttrTable [,spAttr][is.na(spAttrTable [,spAttr])] <- NAfill }
# #                               spValue   = spAttrTable [,spAttr] %>%pull
# # 
# #                               #prepare cbs out
# #                               cbs.mod = Matrix (data = c(0),nrow = nrow(cbs),ncol=ncol(cbs))
# # 
# #                               #create species chunks (default is 10)
# #                               #g = 10
# #                               byArg = round (ncol (cbs.mod)/g,digits = 0)
# #                               a =seq(from = 1,to = ncol (cbs.mod), by= byArg)
# #                               b = a+ (byArg)-1
# #                               b[length(b)] <- ncol(cbs.mod)
# # 
# #                               #loop by species chunks
# #                               for (colChunk in 1:g) {
# #                                 tictoc::tic()
# #                                 spInfoMatrixChunk = Matrix(data = spValue[a[colChunk]:b[colChunk]],
# #                                                            ncol = length(a[colChunk]:b[colChunk]),
# #                                                            nrow = nrow(cbs.mod),
# #                                                            byrow = T)
# # 
# #                                 cellInfoMatrixChunk = Matrix(data = cellValue,
# #                                                              ncol = length(a[colChunk]:b[colChunk]),
# #                                                              nrow = nrow(cbs.mod),
# #                                                              byrow = F)
# # 
# #                                 cbs.Chunk = cbs[,a[colChunk]:b[colChunk]]
# #                                 binaryMatrixCond = eval(expr = parse(text = paste('spInfoMatrixChunk' ,condition ,'cellInfoMatrixChunk')))
# # 
# #                                 out = cbs.Chunk * binaryMatrixCond
# #                                 cbs.mod [,a[colChunk]:b[colChunk]] <- out
# #                                 tictoc::toc()
# # 
# #                               }
# # 
# #                               #write output
# #                               dimnames (cbs.mod)<- dimnames(cbs)
# #                               saveRDS(object = cbs.mod, file = paste0(outDir,'/',basename(cbs.f[x])))
# # 
# #                             }
# # 
# # 
# #   registerDoSEQ()
# #   stopCluster(cl)
# # 
# #   message ('done!')
# # 
# # 
# # 
# # }
