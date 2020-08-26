
#============================================================
#============================================================
#============================================================
#' Plot a raster to a pdf


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



#' @param r a raster
#' @param plotFile optiona; path to for writing pdf
#' @param shp optional shapefile to plot
#' @param legend.args see `plot.raster`
#' @param 									 axis.args=list(cex.axis=1.1),
#' @param  ... arguments to be passed to plot
#' @export
fdMapPlot=function(r,
									 plotFile=NULL,
								   shp=NULL,
									 legend.args=list(text='# species',line=1,side=3,cex=1),
									 axis.args=list(cex.axis=1.1),
									 open=F,
									 ...){

	if(!is.null(plotFile)) pdf(plotFile,h=4*nlayers(r),w=8)
	par(mfrow=c(nlayers(r),1),mar=c(0,0,2,0),oma=c(0,0,0,0))

	zlims=c(min(minValue(r)),max(maxValue(r)))
	for(i in 1:nlayers(r)){
		plot(r[[i]],zlim=zlims,axes=F,xlab="",ylab="",xaxt='n',yaxt='n', box=FALSE,col= c('white',cm.cols1(100)),legend.args=legend.args, axis.args=axis.args,smallplot= c(.89,.91,.1,.9),...)
		mtext(names(r)[i],3,line=0)
		if(!is.null(shp)) plot(shp,add=T,lwd=1)
	}
	if(!is.null(plotFile)){
		dev.off()
		if(open) system(paste0('open ',plotFile))
	}
}

#============================================================
#============================================================
#============================================================
#' @title Map a species based on cbs
#' @description Quickly map species stored in sparse matrices
#' @param
#' @param species a vector of species names or integers refering to species indices
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

mapSpecies=function(cbsDir,species,scenario,sp.ind,cell.ind,envGrid,shp,...){
	#  for testing
	#  cbsDir=sumDirs$cbsDir; species=1:3; shp=world.shp
	if(length(species)>12) warning("probably shouldn't try to plot more than 12 species at once")
	if(!is.numeric(species)) species = which(sp.ind$species %in% species)
	sp.names=sp.ind$species[species]

	cbs.f=changeRangeR:::.getCBS(cbsDir,scenario)

	if(Sys.info()["sysname"]== "Windows") mclapply=parallelsugar::mclapply
	if(Sys.info()["sysname"]!= "Windows") mclapply=parallel::mclapply

	spMaps=mclapply(seq_along(cbs.f), function(x){
		readRDS(cbs.f[x])[,species]
	},mc.cores=mc.cores)
	sp.mat=do.call('rbind',spMaps)
	# this probably isn't fast
	sp.r=mclapply(1:ncol(sp.mat),function(x){
		r=raster(envGrid)
		values(r)[cell.ind$cellID]= sp.mat[,x]
		r
	},mc.cores=mc.cores) %>% stack
	names(sp.r)=sp.names

	plot(sp.r)
	sp.r
}



#============================================================
#============================================================
#============================================================
#' FUNCTIONS TO COMPARE SCENARIOS

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

#' @param
#' @export

#compare richness rasters
compareRichness = function (sumDirs,scn1, scn2,plotFig=T){
  #for testint
  # scn1 = "present"
  #scn2 = '8580'
  f1 = list.files (sumDirs$rich,pattern=scn1,full.names = T)
  f2 = list.files (sumDirs$rich,pattern=scn2,full.names = T)

  if (length(f1)>1){stop('more than one raster for scn1')}
  if (length(f1)== 0){stop('no raster for scn1')}
  if (length(f2)>1){stop('more than one raster for scn2')}
  if (length(f2)== 0){stop('no raster for scn2')}

  r1 = raster(f1) ; r2 = raster (f2)

  #species richness change
  rdiff = r2-r1

  #Positive vs negative
  rDyn = rdiff
  rDyn[values (rDyn) >0 ] = 1
  rDyn[values (rDyn) <0 ] = -1
  if (plotFig==T){
    rasterVis::levelplot(rdiff)
  }

  raster::stack (list (deltaRichness=rdiff,binRichnessChange=rDyn))
}


#============================================================
#============================================================
#============================================================
#' FUNCTIONS TO COMPARE SCENARIOS


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

#calculate turnover
turnoverFromCBS=function(cbsDir,
												 scn1,
												 scn2,
												 env,
												 mc.cores,
												 betaDivChange=F,
												 outputTable=F){
  #  for testing
  # cbsDir=sumDirs$cbsDir
  # scn1='present'
  # scn2='8580'
  # mc.cores=15

  cbsList = lapply (list(scn1,scn2), function (x){
    cbs.f=list.files(paste0(cbsDir,'/',x),full.names=T)
    toss=grep('temp_|spCellOcc|uniqueComm_',cbs.f)
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
  srcByCell=foreach::foreach(x= 1:nchunks,
                             .options.snow=opts,
                             .packages = c((.packages())),
                             #.export =c('cbsList'),
                             .combine = rbind)%dopar%{

    message(x)

    cbsScenarios = lapply (cbsList, function (l,i=x){
      readRDS(l[i])

    })

    cbsSpDyn = cbsScenarios[[2]]-cbsScenarios[[1]]

    Lost <-  Matrix::rowSums(cbsSpDyn == -1)
    Gained <-  Matrix::rowSums (cbsSpDyn == 1)
    Maintained  <-  Matrix::rowSums((cbsScenarios[[2]]+cbsScenarios[[1]]) == 2)
    Turnover = round (((Lost+Gained)/(Lost+Maintained+Gained))*100,digits=0)

    tibble::tibble (cellID=as.numeric(rownames(cbsSpDyn)),
               lost=Lost,
               gained=Gained,
               maintained=Maintained,
               turnover=Turnover)
                              }
  stopCluster(cl)
  registerDoSEQ()


  #build temporal beta diversity
  if (betaDivChange==T){

    betaTempByCell=lapply (1:nchunks, function (x) {

      message(x)

      cbsScenarios = lapply (cbsList, function (l,i=x){
        readRDS(l[i])

      })
      cbsT0 = cbsScenarios[[1]]
      cbsT1 = cbsScenarios[[2]]
      #get the partitioning
      levelsRows = seq (from=1, to = nrow (cbsT0), by = 1000)
      if(length(levelsRows)==1) {rowGroups = rep(1,cbsT0)}
      if(length(levelsRows)> 1) {rowGroups = cut(1:nrow(cbsT0), length(levelsRows), labels = FALSE)}

      #apply beta temporal diversity in chunks of 1000
      # make it Windows compatible :/
      if (Sys.info()["sysname"]== "Windows") {mclapply <- parallelsugar::mclapply}
      if (Sys.info()["sysname"]!= "Windows") {mclapply <- parallel::mclapply}

      betaTemp=mclapply(seq_along(levelsRows),function(g,rg=rowGroups){
        message (paste0('Computing group : ',g,'/',max(rg)))


        cbsSubset = lapply (cbsScenarios, function (x,...){
          y = x[which(rg == g),]
          y = as.matrix(y)
          y
        })

        col2Reduce = which (colSums (cbsSubset[[1]]+cbsSubset[[2]])==0)
        if (length(col2Reduce)==ncol(cbsSubset[[1]])) {
          out=data.frame(beta.sim=rep(NA,nrow(cbsSubset[[1]]),row.names = rownames(cbsSubset)))
          out$beta.sne=NA ; out$beta.sor=NA
          return (out)
        }
        cbsSubset = lapply (cbsSubset, function (x,...){
          x[,col2Reduce*(-1)]
        })

        out = betapart::beta.temp(cbsSubset[[1]],cbsSubset[[2]])
        out
      },mc.cores=mc.cores)

      betaTemp= do.call (rbind,betaTemp)
      dplyr::as_tibble(betaTemp)


    })
    betaTempByCell = do.call (rbind, betaTempByCell)
    #this should be ok but I added for assurance
    betaTempByCell =as_tibble(betaTempByCell)

  }

  #Build raster and assign values
  if (betaDivChange==T){srcByCell = dplyr::bind_cols(srcByCell,betaTempByCell)}
  rDyn=raster(env[[1]])
  rDynStack = lapply (2:ncol(srcByCell), function (j){
    r=rDyn
    r[]=NA
    r[srcByCell$cellID] <-srcByCell[[j]]
    r
  })
  rDynStack = raster::stack(rDynStack)
  names (rDynStack) = names (srcByCell)[2:ncol(srcByCell)]

  #return objects
  #ouptu tibble?
  if(outputTable==T) {list (srcByCell,rDynStack)} else {rDynStack}

}



