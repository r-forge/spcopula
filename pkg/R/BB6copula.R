#####################
##                 ##
## the BB6 copulas ##
##                 ##
#####################
# Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall. 

validBB6Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB6 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param >= upper | param < lower))
    return("Parameter value out of bound.")
  else return (TRUE)
}

setClass("BB6Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB6Copula,
  contains = list("copula")
)

# constructor
BB6Copula <- function (param) {
    val <- new("BB6Copula", dimension = 2, parameters = param, 
        param.names = c("theta", "delta"), param.lowbnd = c(1, 1), param.upbnd = c(Inf, Inf), family=8, message = "BB6 copula family. Number 8 in CDVine.")
    val
}

## density ##
setMethod("dcopula", signature("BB6Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("BB6Copula"), linkCDVine.CDF)

## partial derivatives ##
# ddu
setMethod("dducopula", signature("BB6Copula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("BB6Copula"), linkCDVine.ddv)

## random number generater ??
setMethod("rcopula", signature("BB6Copula"), linkCDVine.r)

## kendall distribution/measure, taken from CDVine:::obs.stat
kendall.BB6 <- function(copula, t){
  theta = copula@parameters[1]
  delta = copula@parameters[2]
  
  kt <- rep(NA,length(t))
  kt <- t + log(-(1 - t)^theta + 1) * (1 - t - (1 - t)^(-theta) + (1 - t)^(-theta) * t)/(delta * theta)
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("BB6Copula"), kendall.BB6)

setMethod("getKendallDistr", signature("BB6Copula"), 
          function(copula) return(function(t) kendall.BB6(copula, t)))

#########################
## BB6 survival copula ##
#########################

setClass("surBB6Copula",
  representation = representation("copula", family="numeric"),
  validity = validBB6Copula,
  contains = list("copula")
)

# constructor
surBB6Copula <- function (param) {
  val <- new("surBB6Copula", dimension = 2, parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(1, 1), param.upbnd = c(Inf, Inf), family=18, message = "Survival BB6 copula family. Number 18 in CDVine.")
  val
}

## density ##
setMethod("dcopula", signature("surBB6Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("surBB6Copula"), linkCDVine.surCDF)
# persp(surBB6Copula(c(5.329995,2.1201476)),dcopula,zlim=c(0,20))
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("surBB6Copula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("surBB6Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("surBB6Copula"), linkCDVine.r)

#######################
## BB6 copula 90 deg ##
#######################

validRotBB6Copula = function(object) {
  if (object@dimension != 2)
    return("Only BB6 copulas of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param > upper | param <= lower))
    return("Parameter value out of bound")
  else return (TRUE)
}

setClass("r90BB6Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB6Copula,
  contains = list("copula")
)

# constructor
r90BB6Copula <- function (param) {
  val <- new("r90BB6Copula", dimension = 2, parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -Inf), param.upbnd = c(-1, -1), family=28, message = "90 deg rotated BB6 copula family. Number 28 in CDVine.")
  val
}

## density ##
setMethod("dcopula", signature("r90BB6Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r90BB6Copula"), linkCDVine.r90CDF)
# persp(r90BB6Copula(c(-1.329995,-1.1201476)), pcopula)
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("r90BB6Copula"), linkCDVine.ddu)

## ddv
setMethod("ddvcopula", signature("r90BB6Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r90BB6Copula"), linkCDVine.r)
# rcopula(r90BB6Copula(c(-5.329995,-1.1201476)),500)

#####################
## BB6 copula 270� ##
#####################

setClass("r270BB6Copula",
  representation = representation("copula", family="numeric"),
  validity = validRotBB6Copula,
  contains = list("copula")
)

# constructor
r270BB6Copula <- function (param) {
  val <- new("r270BB6Copula", dimension = 2, parameters = param, param.names = c("theta", "delta"), param.lowbnd = c(-Inf, -Inf), param.upbnd = c(-1, -1), family=38, message = "270 deg rotated BB6 copula family. Number 38 in CDVine.")
  val
}

## density ##
setMethod("dcopula", signature("r270BB6Copula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r270BB6Copula"), linkCDVine.r270CDF)
# persp(r270BB6Copula(c(-5.329995,-1.1201476)), dcopula)
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("r270BB6Copula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("r270BB6Copula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r270BB6Copula"), linkCDVine.r)