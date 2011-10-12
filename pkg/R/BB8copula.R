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
  val <- new("BB8Copula", dimension = 2, parameters = param, 
  param.names = c("theta", "delta"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, 1), family=10, message = "BB8 copula family. Number 10 in CDVine.")
  val
}

## density ##
setMethod("dcopula", signature("BB8Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("BB8Copula"), linkCDVine.CDF)

# persp(BB8Copula(c(5.329995,0.1201476)),dcopula)

## partial derivatives ##
## ddu
setMethod("dducopula", signature("BB8Copula"),linkCDVine.ddu)

## ddv
setMethod("ddvcopula", signature("BB8Copula"),linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("BB8Copula"),linkCDVine.r)
# rcopula(BB8Copula(c(5.329995,0.1201476)),500)

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
  val <- new("surBB8Copula", dimension = 2, parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, 1), family=20, message = "Survival BB8 copula family. Number 20 in CDVine.")
  val
}

## density ##
setMethod("dcopula", signature("surBB8Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("surBB8Copula"), linkCDVine.surCDF)
# persp(surBB8Copula(c(5.329995, 0.9201476)),pcopula)  

## partial derivatives ##
## ddu
setMethod("dducopula", signature("surBB8Copula"), linkCDVine.ddu)

## ddv
setMethod("ddvcopula", signature("surBB8Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("surBB8Copula"), linkCDVine.r)
# rcopula(surBB8Copula(c(5.329995,0.9201476)),500)

####################
## BB8 copula 90° ##
####################

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
  val <- new("r90BB8Copula", dimension = 2, parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -1), param.upbnd = c(-1, 0), family=30, message = "90° rotated BB8 copula family. Number 30 in CDVine.")
  val
}

## density ##
setMethod("dcopula", signature("r90BB8Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r90BB8Copula"), linkCDVine.r90CDF)
# persp(r90BB8Copula(c(-5.329995,-0.1201476)), dcopula)
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("r90BB8Copula"), linkCDVine.ddu)

## ddv
setMethod("ddvcopula", signature("r90BB8Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r90BB8Copula"), linkCDVine.r)
# rcopula(r90BB8Copula(c(-5.329995,-0.1201476)),500)

#####################
## BB8 copula 270° ##
#####################

setClass("r270BB8Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB8Copula,
  contains = list("copula")
)

# constructor
r270BB8Copula <- function (param) {
  val <- new("r270BB8Copula", dimension = 2, parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -1), param.upbnd = c(-1, 0), family=40, message = "270° rotated BB8 copula family. Number 40 in CDVine.")
  val
}

## density ##
setMethod("dcopula", signature("r270BB8Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r270BB8Copula"), linkCDVine.r270CDF)

# persp(r270BB8Copula(c(-5.329995,-0.1201476)), pcopula)
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("r270BB8Copula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("r270BB8Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r270BB8Copula"), linkCDVine.r)