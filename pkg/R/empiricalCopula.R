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
  empiricalCopula(rCopula(sample.size, copula), copula)
}


## density, not yet needed and hence not implemented ##

## jcdf ##
# from package copula
pempCop.C <- function(u, copula) {
    return(Cn(copula@sample,u))
}

setMethod("pCopula", signature("numeric", "empiricalCopula"),
          function(u, copula, ...) {
            pempCop.C(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix", "empiricalCopula"), pempCop.C)


tauempCop <- function(copula){
  TauMatrix(copula@sample)[1,2]
}

setMethod("tau",signature("empiricalCopula"), tauempCop)


rhoempCop <- function(copula){
  cor(copula@sample,method="spearman")
}

setMethod("rho",signature("empiricalCopula"), rhoempCop)

setMethod("lambda", signature("empiricalCopula"), 
          function(copula, ...) stop("No evaluation possible, try to plot 'empBivJointDepFun' for a visual assessment."))

# Vine Copula - empirical evaluation
## jcdf ##
pvineCopula <- function(u, copula) {
  empCop <- genEmpCop(copula, 1e5)
  
  return(pCopula(u, empCop))
}

setMethod("pCopula", signature("numeric","vineCopula"), 
          function(u, copula) {
            pvineCopula(matrix(u, ncol=copula@dimension), copula)
          })
setMethod("pCopula", signature("data.frame","vineCopula"), 
          function(u, copula) pvineCopula(as.matrix(u), copula))
setMethod("pCopula", signature("matrix","vineCopula"), pvineCopula)