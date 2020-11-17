fdBetaTemp=function(scen1path,scen2path,outDir=NULL,indexFam,verbose=T){
	#  for testing
	#  scen1path= '/Volumes/cm/BIEN_42_Diversity/SparseMatrices/CellBySpecies/present_OccEco_TP05'; scen2path='/Volumes/cm/BIEN_42_Diversity/SparseMatrices/CellBySpecies/bc8590_OccEco_TP05'; outDir=sumDirs$betaTempDir ; indexFam='sorensen'
	
	t1=proc.time()
	message(paste0('starting ',basename(scen1path),' and ',basename(scen2path)))
	cbs.f1=changeRangeR:::.getCBS(dirname(scen1path),basename(scen1path))
	cbs.f2=changeRangeR:::.getCBS(dirname(scen2path),basename(scen2path))
	if(Sys.info()["sysname"]== "Windows") mclapply=parallelsugar::mclapply
	if(Sys.info()["sysname"]!= "Windows") mclapply=parallel::mclapply

	out.tmp=mclapply(seq_along(cbs.f1),function(i){
	  if(verbose) message(i)
		cbs1=readRDS(cbs.f1[i])
		cbs2=readRDS(cbs.f2[i])
		if(!nrow(cbs1)==nrow(cbs2)) stop(paste0('chunk',i,'has different numbers of rows'))
		# remove columns of species not in either because these objects are too large to make data frames from
		cbs1.keep=textTinyR::sparse_Sums(cbs1, rowSums = F)
		cbs2.keep=textTinyR::sparse_Sums(cbs2, rowSums = F)
		keepCols=which(cbs1.keep>0 & cbs2.keep>0)
		cbs1.keep2=textTinyR::sparse_Sums(cbs1, rowSums = T)
		cbs2.keep2=textTinyR::sparse_Sums(cbs2, rowSums = T)
		keepRows=which(cbs1.keep2>0 & cbs2.keep2>0)
		x=as.data.frame(as.matrix(cbs1[keepRows,keepCols]))
		y=as.data.frame(as.matrix(cbs2[keepRows,keepCols]))
		
		t1=proc.time()
		# stolen from betapart::beta.temp to speed up
		ai <- apply(x & y, 1, sum)
    bi <- apply(x & !y, 1, sum)
    ci <- apply(!x & y, 1, sum)
    switch(indexFam, sorensen = {
        beta.sor <- (bi + ci)/(2 * ai + bi + ci)
        beta.sim <- pmin(bi, ci)/(ai + pmin(bi, ci))
        beta.sne <- beta.sor - beta.sim
        result <- data.frame(beta.sim, beta.sne, beta.sor)
    }, jaccard = {
        beta.jac <- (bi + ci)/(ai + bi + ci)
        beta.jtu <- 2 * pmin(bi, ci)/(ai + (2 * pmin(bi, ci)))
        beta.jne <- beta.jac - beta.jtu
        result <- data.frame(beta.jtu, beta.jne, beta.jac)
    })
		#out=betapart::beta.temp(x,y,indexFam)	
		(t2=proc.time()-t1)
		message(paste0('chunk ',i,' ',round(t2[3],1),' s'))
		out.df=data.frame(matrix(NA,nrow=nrow(cbs1),ncol=3))
		names(out.df)=names(out)
		out.df[keepRows,]=result
		out.df
	},mc.cores=mc.cores)
	# should join this with cell.ind to keep cell IDs stored
	out2=do.call(rbind,out.tmp)
	if(!is.null(outDir)){
		saveRDS(out2,file=paste0(outDir,'/',basename(scen1path),'_', basename(scen2path),'.rds'))
	}
	out2
}
