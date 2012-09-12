#####################
##                 ##
## the BB8 copulas ##
##                 ##
#####################
# Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall. 

validBB8Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB8 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param)) | param[1] >= upper[1] | param[2] > upper[2] | param[1] < lower[1] | param[2] <= lower[2])
    return("Parameter value out of bound.")
  else return (TRUE)
}

setClass("BB8Copula",
  representation = representation("copula",family="numeric"),
  validity = validBB8Copula,
  contains = list("copula")
)

# constructor
BB8Copula <- function (param) {
  val <- new("BB8Copula", dimension = as.integer(2), parameters = param, 
  param.names = c("theta", "delta"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, 1), family=10, fullname = "BB8 copula family. Number 10 in CDVine.")
  val
}

## density ##
setMethod("dCopula", signature("numeric","BB8Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","BB8Copula"), linkCDVine.CDF)

## partial derivatives ##
## ddu
setMethod("dduCopula", signature("numeric","BB8Copula"),linkCDVine.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","BB8Copula"),linkCDVine.ddv)

## random number generator
setMethod("rCopula", signature("numeric","BB8Copula"),linkCDVine.r)

## kendall distribution/measure, taken from CDVine:::obs.stat
kendall.BB8 <- function(copula, t){
  theta = copula@parameters[1]
  delta = copula@parameters[2]
  
  kt <- rep(NA,length(t))
  kt <- t + log(((1 - t * delta)^theta - 1)/((1 - delta)^theta - 1)) * (1 - t * delta - (1 - t * delta)^(-theta) + (1 - t * delta)^(-theta) * t * delta)/ (theta * delta)
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("BB8Copula"), kendall.BB8)

setMethod("getKendallDistr", signature("BB8Copula"), 
          function(copula) return(function(t) kendall.BB8(copula, t)))

#########################
## BB8 survival copula ##
#########################

setClass("surBB8Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB8Copula,
  contains = list("copula")
)

# constructor
surBB8Copula <- function (param) {
  val <- new("surBB8Copula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, 1), family=20, fullname = "Survival BB8 copula family. Number 20 in CDVine.")
  val
}

## density ##
setMethod("dCopula", signature("numeric","surBB8Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","surBB8Copula"), linkCDVine.surCDF)

## partial derivatives ##
## ddu
setMethod("dduCopula", signature("numeric","surBB8Copula"), linkCDVine.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","surBB8Copula"), linkCDVine.ddv)

## random number generator
setMethod("rCopula", signature("numeric","surBB8Copula"), linkCDVine.r)

#######################
## BB8 copula 90 deg ##
#######################

validRotBB8Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB8 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param[1] > upper[1] | param[2] >= upper[2] | param[1] <= lower[1] | param[2] < lower[2]))
    return("Parameter value out of bound")
  else return (TRUE)
}

setClass("r90BB8Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB8Copula,
  contains = list("copula")
)

# constructor
r90BB8Copula <- function (param) {
  val <- new("r90BB8Copula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -1), param.upbnd = c(-1, 0), family=30, fullname = "90 deg rotated BB8 copula family. Number 30 in CDVine.")
  val
}

## density ##
setMethod("dCopula", signature("numeric","r90BB8Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r90BB8Copula"), linkCDVine.r90CDF)
  
## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90BB8Copula"), linkCDVine.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r90BB8Copula"), linkCDVine.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90BB8Copula"), linkCDVine.r)

#####################
## BB8 copula 270� ##
#####################

setClass("r270BB8Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB8Copula,
  contains = list("copula")
)

# constructor
r270BB8Copula <- function (param) {
  val <- new("r270BB8Copula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -1), param.upbnd = c(-1, 0), family=40, fullname = "270 deg rotated BB8 copula family. Number 40 in CDVine.")
  val
}

## density ##
setMethod("dCopula", signature("numeric","r270BB8Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r270BB8Copula"), linkCDVine.r270CDF)
  
## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270BB8Copula"), linkCDVine.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","r270BB8Copula"), linkCDVine.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270BB8Copula"), linkCDVine.r)