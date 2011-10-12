#####################
##                 ##
## the BB7 copulas ##
##                 ##
#####################
# Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall. 

validBB7Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB7 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param >= upper | param[1] < lower[1] | param[2] <= lower[2]))
    return("Parameter value out of bound.")
  else return (TRUE)
}

setClass("BB7Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB7Copula,
  contains = list("copula")
)

# constructor
BB7Copula <- function (param) {
    val <- new("BB7Copula", dimension = 2, parameters = param, 
        param.names = c("theta", "delta"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, Inf), family=9, message = "BB7 copula family. Number 9 in CDVine.")
    val
}

## density ##
setMethod("dcopula", signature("BB7Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("BB7Copula"), linkCDVine.CDF)

## partial derivatives ##
# ddu
setMethod("dducopula", signature("BB7Copula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("BB7Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("BB7Copula"), linkCDVine.r)


#########################
## BB7 survival copula ##
#########################

setClass("surBB7Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB7Copula,
  contains = list("copula")
)

# constructor
surBB7Copula <- function (param) {
  val <- new("surBB7Copula", dimension = 2, parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(1, 0), param.upbnd = c(Inf, Inf), family= 19, message = "Survival BB7 copula family. Number 19 in CDVine.")
  return(val)
}

## density ##
setMethod("dcopula", signature("surBB7Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("surBB7Copula"), linkCDVine.surCDF)
# persp(surBB7Copula(c(5.329995,2.1201476)),pcopula,zlim=c(0,1))
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("surBB7Copula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("surBB7Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("surBB7Copula"), linkCDVine.r)

####################
## BB7 copula 90� ##
####################

validRotBB7Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB7 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param[1] > upper[1] | param[2] >= upper[2] | param <= lower))
    return("Parameter value out of bound")
  else return (TRUE)
}

setClass("r90BB7Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB7Copula,
  contains = list("copula")
)

# constructor
r90BB7Copula <- function (param) {
  val <- new("r90BB7Copula", dimension = 2, parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -Inf), param.upbnd = c(-1, 0), family=29, message = "90° rotated BB7 copula family. Number 29 in CDVine.")
  val
}

## density ##
setMethod("dcopula", signature("r90BB7Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r90BB7Copula"), linkCDVine.r90CDF)
# persp(r90BB7Copula(c(-1.329995,-1.1201476)), dcopula)
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("r90BB7Copula"), linkCDVine.ddu)

## ddv
setMethod("ddvcopula", signature("r90BB7Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r90BB7Copula"), linkCDVine.r)
# rcopula(r90BB7Copula(c(-5.329995,-1.1201476)),500)

#####################
## BB7 copula 270° ##
#####################

setClass("r270BB7Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB7Copula,
  contains = list("copula")
)

# constructor
r270BB7Copula <- function (param) {
  val <- new("r270BB7Copula", dimension = 2, parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -Inf), param.upbnd = c(-1, -1), family=39, message = "270° rotated BB7 copula family. Number 39 in CDVine.")
  val
}

## density ##
setMethod("dcopula", signature("r270BB7Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r270BB7Copula"), linkCDVine.r270CDF)
# persp(r270BB7Copula(c(-5.329995,-1.1201476)), dcopula, zlim=c(0,20))
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("r270BB7Copula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("r270BB7Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r270BB7Copula"), linkCDVine.r)