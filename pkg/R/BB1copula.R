#####################
##                 ##
## the BB1 copulas ##
##                 ##
#####################
# Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall. 

validBB1Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB1 copulas of dimension 2 are supported.")
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

setClass("BB1Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB1Copula,
  contains = list("copula")
)

# constructor
BB1Copula <- function (param) {
    val <- new("BB1Copula", dimension = as.integer(2), parameters = param, 
        param.names = c("theta", "delta"), param.lowbnd = c(0, 1), param.upbnd = c(Inf, Inf), family=7, message = "BB1 copula family. Number 7 in CDVine.")
    val
}

## density ##
setMethod("dcopula", signature("BB1Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("BB1Copula"), linkCDVine.CDF)

## partial derivatives ##
# ddu
setMethod("dducopula", signature("BB1Copula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("BB1Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("BB1Copula"), linkCDVine.r)

## kendall distribution/measure
kendall.BB1 <- function(copula, t){
  theta = copula@parameters[1]
  delta = copula@parameters[2]
  
  kt <- rep(NA,length(t))
  kt <- t + 1/(theta * delta) * (t^(-theta) - 1)/(t^(-1 - theta))
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("BB1Copula"), kendall.BB1)

setMethod("getKendallDistr", signature("BB1Copula"), function(copula) return(function(t) kendall.BB1(copula, t)) )

#########################
## BB1 survival copula ##
#########################

setClass("surBB1Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB1Copula,
  contains = list("copula")
)

# constructor
surBB1Copula <- function (param) {
  val <- new("surBB1Copula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(0, 1), param.upbnd = c(Inf, Inf), family=17, message = "Survival BB1 copula family. Number 17 in CDVine.")
  val
}

## density ##
setMethod("dcopula", signature("surBB1Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("surBB1Copula"), linkCDVine.surCDF)
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("surBB1Copula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("surBB1Copula"), linkCDVine.ddv)

## random number generater ??
setMethod("rcopula", signature("surBB1Copula"), linkCDVine.r)

#######################
## BB1 copula 90 deg ##
#######################

validRotBB1Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB1 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param[1] >= upper[1] | param[2] > upper[2] | param <= lower))
    return("Parameter value out of bound")
  else return (TRUE)
}

setClass("r90BB1Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB1Copula,
  contains = list("copula")
)

# constructor
r90BB1Copula <- function (param) {
  val <- new("r90BB1Copula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -Inf), param.upbnd = c(0, -1), family=27, message = "90 deg rotated BB1 copula family. Number 27 in CDVine.")
  val
}

## density ##
setMethod("dcopula", signature("r90BB1Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r90BB1Copula"), linkCDVine.r90CDF)

## partial derivatives ##
# ddu
setMethod("dducopula", signature("r90BB1Copula"), linkCDVine.ddu)

## ddv
setMethod("ddvcopula", signature("r90BB1Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r90BB1Copula"), linkCDVine.r)

########################
## BB1 copula 270 deg ##
########################

setClass("r270BB1Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB1Copula,
  contains = list("copula")
)

# constructor
r270BB1Copula <- function (param) {
  val <- new("r270BB1Copula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -Inf), param.upbnd = c(0, -1), family=37, message = "270 deg rotated BB1 copula family. Number 37 in CDVine.")
  val
}

## density ##
setMethod("dcopula", signature("r270BB1Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r270BB1Copula"), linkCDVine.r270CDF)
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("r270BB1Copula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("r270BB1Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r270BB1Copula"), linkCDVine.r)