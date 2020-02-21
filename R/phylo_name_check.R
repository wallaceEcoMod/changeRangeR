# ### make this into a function that checks names of species against tip labels
#
# ### Notes and code on primate phylogenies for PD/PE calculations on changeRangeR package
# setwd("~/Dropbox/Wallace_project_Andrea/PD_primates/Primates_all/")
# library(ape)
# trees<-read.nexus("Primates_col_all_1k.nex")
# setwd("~/Dropbox/Wallace_project_Andrea/PD_primates")
# names<-read.table("primate_names_for_phylo_feb2020.txt",head=F)
# names$V1<-as.character(names$V1)
# library(geiger)
# ###To get list of species not in tree
# geiger::name.check(trees[[1]], names, data.names=names$V1)$data_not_tree
# geiger::name.check(trees[[1]], names, data.names=names$v1$tree_not_data
# ## Print a warning of the names that do not match. Give option to drop all non-analogs, or rename outside of function
