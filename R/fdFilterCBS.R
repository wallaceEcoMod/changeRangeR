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

filterCBS <- function (sumDir,sumDirOut, spSelec=NULL, cellIdSelec=NULL) {

  #determine input and output directory
  #sumDir = "E:/Pep/TREECHANGE/SDM_ForDivFate/TCCC_WORLD_v6_10km__Summaries_sdmTP05rbg0165thr"
  cbsDir = paste0 (sumDir,'/CellSpeciesLists/')
  #sumDirOut = "E:/Pep/TREECHANGE/SDM_ForDivFate/treeSelec/"
  if (!dir.exists(sumDirOut)) dir.create(sumDirOut,recursive = T)

  #list cbs rds files
  cbs.f=list.files(cbsDir,full.names=T,pattern = 'chunk_',recursive = T)
  if (length(cbs.f)==0) stop ('No cellBySpecies.RDS files found, cannot compute unique communities')

  #determine input and output directory


  if (!is.null(spSelec)){
    spId = readr::read_csv (paste0 (sumDir,'/speciesIndexTable.csv'))
    spId = subset (spId, species %in% spSelec)
    spIdNewIndex = spId
    spIdNewIndex$index = 1:nrow (spIdNewIndex)
    readr::write_csv(x =spIdNewIndex,path = paste0 (sumDirOut,'/speciesIndexTable.csv') )
  }

  if (!is.null(cellIdSelec)){
     print ("subseting by cellInd not implemented yet ")
  }

  #select and copy cbs
  lapply (cbs.f, function (cbsFile,colSelec =!is.null(spSelec), rowSelec= !is.null(cellIdSelec)){

    cbs = readRDS(cbsFile)

    if (colSelec) {
      cbsSelec = cbs[,spId$index]
      cbsFileOut = strsplit(cbsFile,split=sumDir)[[1]][2]
      cbsFileOut = paste0(sumDirOut,cbsFileOut)
      dir.create (dirname(cbsFileOut),recursive = T)
      saveRDS(cbsSelec,cbsFileOut)
    }
    if (rowSelec) {
      print ("subseting by cellInd not implemented yet ")
    }



  })


}


#### example

# filterCBS(sumDir ="E:/Pep/TREECHANGE/SDM_ForDivFate/TCCC_WORLD_v6_10km__Summaries_sdmTP05rbg0165thr",
#           sumDirOut = "E:/Pep/TREECHANGE/SDM_ForDivFate/treeSelec/",
#           spSelec =  oneSelectionSp )
#

