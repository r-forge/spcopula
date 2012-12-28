########################################
##                                    ##
## an empirical copula representation ##
##                                    ##
########################################

# constructor
empiricalCopula <- function (sample=NULL, copula) {
  if(is.null(sample) && missing(copula))
    stop("At least one parameter of copula or sample must be provided.")
  if(is.null(sample))
    return(genEmpCop(copula))
  if(missing(copula))
    return(new("empiricalCopula", dimension = as.integer(ncol(sample)), 
               parameters = as.numeric(NA), param.names = "unknown", 
               param.lowbnd = as.numeric(NA), param.upbnd = as.numeric(NA), 
               fullname = "Unkown empirical copula based on a sample.",
               sample=sample))
  new("empiricalCopula", dimension = copula@dimension, 
      parameters = copula@parameters, param.names = copula@param.names, 
      param.lowbnd = copula@param.lowbnd, param.upbnd = copula@param.upbnd, 
      fullname = paste("Empirical copula derived from",copula@fullname),
      sample=sample)
}

# simplified constructor
genEmpCop <- function(copula, sample.size=1e5) {
  cat("Note: the copula will be empirically represented by a sample of size:", sample.size, "\n")
  empiricalCopula(rCopula(sample.size,copula), copula)
}


## density, not yet needed and hence not implemented ##

## jcdf ##
# pempCop <- function(u, copula) {
#   t_data <- t(copula@sample)
#   
#   u <- matrix(u,ncol=copula@dimension)
#     
#   # --/-- make this a C-function?
#   res <- NULL
#   for(i in 1:nrow(u)) {
#     bool <- t_data <= u[i,]
#     for (i in 2:nrow(t_data)) bool[1,] <- bool[1,] * bool[i,]
#     res <- c(res,sum(bool[1,]))
#   }
#   # --//--
#   
#   return(res/ncol(t_data))
# }

# from package copula
pempCop.C <- function(u, copula) {
  .C("RmultCn", as.double(copula@sample), as.integer(nrow(copula@sample)),
     copula@dimension, as.double(u), as.integer(nrow(u)), as.double(u[,1]),
     PACKAGE="copula")[[6]]
}

##
# us <- matrix(runif(100),ncol=2)
# copula <- genEmpCop(normalCopula(.3),1e6)
# 
# hist(pCopula(us,normalCopula(.3)) - pempCop.C(us,copula),main="C")
# hist(pCopula(us,normalCopula(.3)) - pempCop(us,copula),main="R")
# 
# system.time(pempCop.C(us,copula))
# system.time(pempCop(us,copula))

setMethod("pCopula", signature("numeric", "empiricalCopula"),
          function(u, copula, ...) {
            pempCop.C(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix", "empiricalCopula"), pempCop.C)


tauempCop <- function(copula){
  CDVine:::fasttau(copula@sample[,1],copula@sample[,2])
}

setMethod("tau",signature("asCopula"),tauempCop)


rhoempCop <- function(copula){
  cor(copula@sample,method="spearman")
}

setMethod("rho",signature("asCopula"),rhoempCop)