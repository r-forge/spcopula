####################
##                ##
## the Joe copula ##
##                ##
####################
# Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall. 

validJoeCopula = function(object) {
  if (object@dimension != 2)
    return("Only Joe copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param >= upper | param <= lower ))
    return("Parameter value out of bound.")
  else return (TRUE)
}

setClass("JoeCopula",
  representation = representation("copula", family="numeric"),
  validity = validJoeCopula,
  contains = list("copula")
)

# constructor
JoeCopula <- function (param) {
  new("JoeCopula", dimension = as.integer(2), parameters = param, param.names = c("theta"),
      param.lowbnd = 1, param.upbnd = Inf, family=6, 
      fullname = "Joe copula family. Number 6 in CDVine.")
}

## density ##
setMethod("dCopula", signature("numeric","JoeCopula"), 
          function(u, copula, log) {
            linkCDVine.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","JoeCopula"), function(u, copula, log) linkCDVine.PDF(u, copula, log))

## jcdf ##
setMethod("pCopula", signature("numeric","JoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","JoeCopula"), linkCDVine.CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","JoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","JoeCopula"), linkCDVine.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","JoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","JoeCopula"), linkCDVine.ddv)

## random number generater
setMethod("rCopula", signature("numeric","JoeCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("JoeCopula"), 
          function(copula, tau) {
            if(tau <= 0) warning("The Joe copula can only represent positive dependence!")
            linkCDVine.calibKendallsTau(copula, max(1e-6,abs(tau)))
          })

## kendall distribution/measure, taken from CDVine:::obs.stat
kendall.Joe <- function(copula, t){
  par = copula@parameters[1]
  
  kt <- rep(NA,length(t))
  kt <- t - (log(1 - (1 - t)^par) * (1 - (1 - t))^par)/(par * (1 - t)^(par - 1))
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("JoeCopula"), kendall.Joe)

setMethod("getKendallDistr", signature("JoeCopula"), 
          function(copula) return(function(t) kendall.Joe(copula, t)))

#########################
## Joe survival copula ##
#########################

setClass("surJoeCopula",
  representation = representation("copula", family="numeric"),
  validity = validJoeCopula,
  contains = list("copula")
)

# constructor
surJoeCopula <- function (param) {
  new("surJoeCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"),
      param.lowbnd = 1, param.upbnd = Inf, family=16, 
      fullname = "Survival Joe copula family. Number 16 in CDVine.")
}

## density ##
setMethod("dCopula", signature("numeric","surJoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.PDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dCopula", signature("matrix","surJoeCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","surJoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","surJoeCopula"), linkCDVine.surCDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","surJoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","surJoeCopula"), linkCDVine.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","surJoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","surJoeCopula"), linkCDVine.ddv)

## random number generater
setMethod("rCopula", signature("numeric","surJoeCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("surJoeCopula"), 
          function(copula, tau) {
            if(tau <= 0) warning("The survival Joe copula can only represent positive dependence!")
            linkCDVine.calibKendallsTau(copula, max(1e-6,abs(tau)))
          })

###################
## Joe copula 90 ##
###################

validRotJoeCopula = function(object) {
  if (object@dimension != 2)
    return("Only Joe copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param >= upper | param <= lower))
    return("Parameter value out of bound.")
  else return (TRUE)
}

setClass("r90JoeCopula",
  representation = representation("copula", family="numeric"),
  validity = validRotJoeCopula,
  contains = list("copula")
)

# constructor
r90JoeCopula <- function (param) {
  new("r90JoeCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"),
      param.lowbnd = -Inf, param.upbnd = -1, family=26, 
      fullname = "90 deg rotated Joe copula family. Number 26 in CDVine.")
}

## density ##
setMethod("dCopula", signature("numeric","r90JoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.PDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dCopula", signature("matrix","r90JoeCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r90JoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r90JoeCopula"), linkCDVine.r90CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90JoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r90JoeCopula"), linkCDVine.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r90JoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r90JoeCopula"), linkCDVine.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90JoeCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("r90JoeCopula"),
          function(copula, tau) {
            if(tau >= 0) warning("The rotated Joe copula can only represent negative dependence!")
            linkCDVine.calibKendallsTau(copula, min(-1e-6,-abs(tau)))
          })

####################
## Joe copula 270 ##
####################

setClass("r270JoeCopula",
  representation = representation("copula", family="numeric"),
  validity = validRotJoeCopula,
  contains = list("copula")
)

# constructor
r270JoeCopula <- function (param) {
  new("r270JoeCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), 
      param.lowbnd = -Inf, param.upbnd = -1, family=36, 
      fullname = "270 deg rotated Joe copula family. Number 36 in CDVine.")
}

## density ##
setMethod("dCopula", signature("numeric","r270JoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.PDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dCopula", signature("matrix","r270JoeCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r270JoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r270JoeCopula"), linkCDVine.r270CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270JoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r270JoeCopula"), linkCDVine.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","r270JoeCopula"), 
          function(u, copula, ...) {
            linkCDVine.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r270JoeCopula"), linkCDVine.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270JoeCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("r270JoeCopula"), 
          function(copula, tau) {
            if(tau >= 0) warning("The rotated Joe copula can only represent negative dependence!")
            linkCDVine.calibKendallsTau(copula, min(-1e-6,-abs(tau)))
          })