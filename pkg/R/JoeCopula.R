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
  if (any(is.na(param) | param >= upper | param < lower ))
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
  new("JoeCopula", dimension = 2, parameters = param, param.names = c("theta", "delta"),
      param.lowbnd = 1, param.upbnd = Inf, family=6, 
      message = "Joe copula family. Number 6 in CDVine.")
}

## density ##
setMethod("dcopula", signature("JoeCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("JoeCopula"), linkCDVine.CDF)

## partial derivatives ##
# ddu
setMethod("dducopula", signature("JoeCopula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("JoeCopula"), linkCDVine.ddv)

## random number generater ??
setMethod("rcopula", signature("JoeCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("JoeCopula"), linkCDVine.calibKendallsTau)

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
  new("surJoeCopula", dimension = 2, parameters = param, param.names = c("theta", "delta"),
      param.lowbnd = 1, param.upbnd = Inf, family=16, 
      message = "Survival Joe copula family. Number 16 in CDVine.")
}

## density ##
setMethod("dcopula", signature("surJoeCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("surJoeCopula"), linkCDVine.surCDF)
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("surJoeCopula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("surJoeCopula"), linkCDVine.ddv)

## random number generater ??
setMethod("rcopula", signature("surJoeCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("surJoeCopula"), linkCDVine.calibKendallsTau)

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
    return("Parameter value out of bound")
  else return (TRUE)
}

setClass("r90JoeCopula",
  representation = representation("copula", family="numeric"),
  validity = validRotJoeCopula,
  contains = list("copula")
)

# constructor
r90JoeCopula <- function (param) {
  new("r90JoeCopula", dimension = 2, parameters = param, param.names = c("theta", "delta"),
      param.lowbnd = -Inf, param.upbnd = -1, family=26, 
      message = "90 deg rotated Joe copula family. Number 26 in CDVine.")
}

## density ##
setMethod("dcopula", signature("r90JoeCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r90JoeCopula"), linkCDVine.r90CDF)

## partial derivatives ##
# ddu
setMethod("dducopula", signature("r90JoeCopula"), linkCDVine.ddu)

## ddv
setMethod("ddvcopula", signature("r90JoeCopula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r90JoeCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("r90JoeCopula"), linkCDVine.calibKendallsTau)

#####################
## Joe copula 270� ##
#####################

setClass("r270JoeCopula",
  representation = representation("copula", family="numeric"),
  validity = validRotJoeCopula,
  contains = list("copula")
)

# constructor
r270JoeCopula <- function (param) {
  new("r270JoeCopula", dimension = 2, parameters = param, param.names = c("theta", "delta"), 
      param.lowbnd = -Inf, param.upbnd = -1, family=36, 
      message = "270 deg rotated Joe copula family. Number 36 in CDVine.")
}

## density ##
setMethod("dcopula", signature("r270JoeCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r270JoeCopula"), linkCDVine.r270CDF)
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("r270JoeCopula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("r270JoeCopula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r270JoeCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("r270JoeCopula"), linkCDVine.calibKendallsTau)