#' @title cbs scenario masking
#' @description modify an output scenario based on binary species rasters
#' @param sumDirs
#' @param allSpeciesMasks
#' @param scenario
#' @param maskName
#' @param nCellChunks
#' @param mc.cores
#' @param
#' @param
# @examples
#
#' @return saves a cbs file
#' @author Cory Merow <cory.merow@@gmail.com>
#' @note This is the typical use for dispersal constrained scenario
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the docum
#' @export
# FUNCTION PARAMETERS


maskCBSscenario <- function (sumDirs,
                             sp.ind,
                             allSpeciesMasks,
                             scenario,
                             maskName,
                             nCellChunks=NULL,
                             mc.cores){

 # for testing
 # spMaskDir=system.file("extdata/treeDemo/Migrationmask/120km/",package='changeRangeR')
 # allSpeciesMasks=tibble(rasterFiles=list.files(spMaskDir, recursive=T, full.names=T)) %>%
 #   mutate(sp.names= rasterFiles %>% basename %>% file_path_sans_ext)
 # sumDirs= sumDirs; allSpeciesMasks = allSpeciesMasks; scenario = '8580';
 # maskName = 'dispers120';nCellChunks = NULL;mc.cores=6
  t1=proc.time()

  #initial congruence tests
  if (length(allSpeciesMasks$sp.names) != nrow(sp.ind)) {warning ('not the same number of species and masks')}
  if (!all(allSpeciesMasks$sp.names %in% as.character (sp.ind$species))) {
    stop('some species masks are not part of your analysis')
    print ( allSpeciesMasks$sp.names [which (!allSpeciesMasks$sp.names %in% sp.ind$species)] )
  }
  exampleRas = raster (allSpeciesMasks$rasterFiles[1])
  if (!compareRaster(exampleRas,envGrid)) stop ('the envGrid and the masks do not coincide')
  rm (exampleRas)
  outDir = sumDirs$cbsDir

  #check if need to build masks with raster[]<- 1 for species w/o masks ]
  ###Pep: not sure if we need/want this.
  buildFakeMasks = ifelse (!all(sp.ind$species %in% spMaskRaster),F ,T )
  if (buildFakeMasks){print('Some species modelled do not have a mask. Automatic mask generation not implemented yet, Talk to Cory'); stop()}

  #get info on nCellChunks from the scenario
  if (is.null (nCellChunks)){
    chunkFiles = textclean::drop_element(x = list.files(cbsScnDir,pattern = '.rds'),pattern = 'temp')
    nCellChunks = length (chunkFiles)
  }

  # check if the masks are done
  ###Pep: not sure if we should have a default sparseDir$CellBySpeciesMasks in setup setupSummaryDirectories, your call.
  sumDirs$cbsMaskDir = paste0(sumDirs$sparseDir,'/CellBySpeciesMasks')
  dir.create (sumDirs$cbsMaskDir,showWarnings = F)
  myMaskDir =  paste0 (sumDirs$cbsMaskDir,'/',maskName)

  if (dir.exists(myMaskDir) ){
    spCellMaskCBS = textclean::drop_element(pattern = 'temp_',list.files (path = myMaskDir , full.names = T))
    doCBSmasks = ifelse (length(spCellMaskCBS) == nCellChunks,F,T)
  }
  if (!dir.exists(myMaskDir))  doCBSmasks = T
  if (!doCBSmasks) (message ('CellBySpeciesMasks sparse matrices are already computed'))
  if (doCBSmasks){
    # build sparse matrices for masks
    message ('Building CellBySpeciesMasks. This could take a while.')
    cellBySpeciesMatrices(outDir = sumDirs$cbsMaskDir,
                          allSpeciesMaps=allSpeciesMasks,
                          scenario=maskName,
                          envGrid=envGrid,
                          sp.ind=sp.ind,
                          cell.ind=cell.ind,
                          nCellChunks=nCellChunks,
                          mc.cores=mc.cores,
                          overwrite=T)
  }

  # prepare output files and folders
  cbsScnDir =  paste0 (sumDirs$cbsDir,'/',scenario)
  spCellCBS = textclean::drop_element(pattern = 'temp_',list.files (path = cbsScnDir,full.names = T))
  spCellMaskCBS = textclean::drop_element(pattern = 'temp_',list.files (path = myMaskDir , full.names = T))
  outputScenario = paste0(sumDirs$cbsDir,'/', scenario,'_',maskName)
  dir.create(outputScenario,showWarnings = F)

  # sparse matrix multiplication and saving
  lapply(1:nCellChunks, function (i){
    message (paste('multiplying chunk', i))
    sdmCBS = readRDS(spCellCBS[i])
    maskCBS = readRDS(spCellMaskCBS[i])
    out = sdmCBS * maskCBS
    saveRDS(out, paste0(outputScenario,'/chunk_',i,'.rds'))
  })

  t2=proc.time()-t1
  message( paste0(round(t2[3],2),' s') )
}

