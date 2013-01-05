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
# from package copula
pempCop.C <- function(u, copula) {
#   Cn(u, copula@sample, do.pobs=F, method="C") # preferred use instead of direct C-code from copula <=0.999-5

  # annoying hack, to be removed after release of copula 0.999-6 with the line above
  if("RmultCn" %in% names(getDLLRegisteredRoutines(getLoadedDLLs()[["copula"]][["path"]])[[1]])) {
    if(getDLLRegisteredRoutines(getLoadedDLLs()[["copula"]][["path"]])[[1]][["RmultCn"]]$numParameters == 6)
       .C("RmultCn", as.double(copula@sample), as.integer(nrow(copula@sample)),
         copula@dimension, as.double(u), as.integer(nrow(u)), as.double(u[,1]),
         PACKAGE="copula")[[6]]
    else
      .C("RmultCn", as.double(copula@sample), as.integer(nrow(copula@sample)),
         copula@dimension, as.double(u), as.integer(nrow(u)), as.double(u[,1]),
         as.double(0), PACKAGE="copula")[[6]]
  } else # copula > 0.999-5
    .C("Cn_C", as.double(copula@sample), as.integer(nrow(copula@sample)),
       copula@dimension, as.double(u), as.integer(nrow(u)), as.double(u[,1]),
       as.double(0), PACKAGE="copula")[[6]]
}

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