oldNames=list.files('/Users/ctg/Dropbox/Projects/Wallace/changeRangeR/inst/extdata/DemoData/BinaryMaps',full.names=T,recursive=T)
newNames=paste0(dirname(oldNames),'.tif')
file.copy(oldNames,newNames)
unlink(dirname(oldNames),recursive=T)