mm=matrix(c(1,3,5.6,7,9,11),nrow=2)
class(mm[1,1])
object.size(mm[1,1])
object.size(mm[1,2])
object.size(mm[2,2])

56*33000^2*1e-9

forceMatrixToInteger <- function(m){
    apply (m, c (1, 2), function (x) {
         (as.integer(x))
    })
}

b=forceMatrixToInteger(mm)
object.size(b[1,1])
# 
# 1.7e6 * 3e5 *56 *1e-9 *.001
# 
# object.size(mm)
# textTinyR::matrix_sparsity(as(mm,"dgCMatrix"))
# 
# rs=read.csv('/Users/ctg/Dropbox/Projects/BIEN/Modeling41/Summaries/RangeSizes/RangeSize_SpAT_Present.csv')
# mean(rs$rangeSize)
# 1-(1.1e4 / 1.6e6)
