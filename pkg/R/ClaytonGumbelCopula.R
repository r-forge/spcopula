#####################################
##                                 ##
## additions to the Clayton copula ##
##                                 ##
#####################################

#############################
## Clayton survival copula ##
#############################

validClaytonCopula = function(object) {
  if (object@dimension != 2)
    return("Only Clayton copulas of dimension 2 are supported.")
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

setClass("surClaytonCopula",
  representation = representation("copula", family="numeric"),
  validity = validClaytonCopula,
  contains = list("copula")
)

# constructor
surClaytonCopula <- function (param) {
  new("surClaytonCopula", dimension = as.integer(2), parameters = param, param.names = c("theta"),
      param.lowbnd = 0, param.upbnd = Inf, family=13, 
      message = "Survival Clayton copula family. Number 13 in CDVine.")
}

## density ##
setMethod("dcopula", signature("surClaytonCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("surClaytonCopula"), linkCDVine.surCDF)
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("surClaytonCopula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("surClaytonCopula"), linkCDVine.ddv)

## random number generater ??
setMethod("rcopula", signature("surClaytonCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("surClaytonCopula"), 
          function(copula, tau) {
            if(tau <= 0) warning("The survival Clayton copula can only represent positive dependence!")
            linkCDVine.calibKendallsTau(copula, max(1e-6,abs(tau)))
          })

#######################
## Clayton copula 90 ##
#######################

validRotClaytonCopula = function(object) {
  if (object@dimension != 2)
    return("Only Clayton copulas of dimension 2 are supported.")
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

setClass("r90ClaytonCopula",
  representation = representation("copula", family="numeric"),
  validity = validRotClaytonCopula,
  contains = list("copula")
)

# constructor
r90ClaytonCopula <- function (param) {
  new("r90ClaytonCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"),
      param.lowbnd = -Inf, param.upbnd = 0, family=23, 
      message = "90 deg rotated Clayton copula family. Number 23 in CDVine.")
}

## density ##
setMethod("dcopula", signature("r90ClaytonCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r90ClaytonCopula"), linkCDVine.r90CDF)

## partial derivatives ##
# ddu
setMethod("dducopula", signature("r90ClaytonCopula"), linkCDVine.ddu)

## ddv
setMethod("ddvcopula", signature("r90ClaytonCopula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r90ClaytonCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("r90ClaytonCopula"),
          function(copula, tau) {
            if(tau >= 0) warning("The rotated Clayton copula can only represent negative dependence!")
            linkCDVine.calibKendallsTau(copula, min(-1e-6,-abs(tau)))
          })

########################
## Clayton copula 270 ##
########################

setClass("r270ClaytonCopula",
  representation = representation("copula", family="numeric"),
  validity = validRotClaytonCopula,
  contains = list("copula")
)

# constructor
r270ClaytonCopula <- function (param) {
  new("r270ClaytonCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), 
      param.lowbnd = -Inf, param.upbnd = 0, family=33, 
      message = "270 deg rotated Clayton copula family. Number 33 in CDVine.")
}

## density ##
setMethod("dcopula", signature("r270ClaytonCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r270ClaytonCopula"), linkCDVine.r270CDF)
  
## partial derivatives ##
# ddu
setMethod("dducopula", signature("r270ClaytonCopula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("r270ClaytonCopula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r270ClaytonCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("r270ClaytonCopula"), 
          function(copula, tau) {
            if(tau >= 0) warning("The rotated Clayton copula can only represent negative dependence!")
            linkCDVine.calibKendallsTau(copula, min(-1e-6,-abs(tau)))
          })

####################################
##                                ##
## additions to the Gumbel copula ##
##                                ##
####################################

validGumbelCopula = function(object) {
  if (object@dimension != 2)
    return("Only Gumbel copulas of dimension 2 are supported.")
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

############################
## Gumbel survival copula ##
############################

setClass("surGumbelCopula",
         representation = representation("copula", family="numeric"),
         validity = validGumbelCopula,
         contains = list("copula")
         )

# constructor
surGumbelCopula <- function (param) {
  new("surGumbelCopula", dimension = as.integer(2), parameters = param, param.names = c("theta"),
      param.lowbnd = 1, param.upbnd = Inf, family=14, 
      message = "Survival Gumbel copula family. Number 14 in CDVine.")
}

## density ##
setMethod("dcopula", signature("surGumbelCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("surGumbelCopula"), linkCDVine.surCDF)

## partial derivatives ##
# ddu
setMethod("dducopula", signature("surGumbelCopula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("surGumbelCopula"), linkCDVine.ddv)

## random number generater ??
setMethod("rcopula", signature("surGumbelCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("surGumbelCopula"), 
          function(copula, tau) {
            if(tau < 0) warning("The survival Gumbel copula can only represent non-negative dependence!")
            linkCDVine.calibKendallsTau(copula, max(0,abs(tau)))
          })

#######################
## Gumbel copula 90 ##
#######################

validRotGumbelCopula = function(object) {
  if (object@dimension != 2)
    return("Only Gumbel copulas of dimension 2 are supported.")
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

setClass("r90GumbelCopula",
         representation = representation("copula", family="numeric"),
         validity = validRotGumbelCopula,
         contains = list("copula")
         )

# constructor
r90GumbelCopula <- function (param) {
  new("r90GumbelCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"),
      param.lowbnd = -Inf, param.upbnd = -1, family=24, 
      message = "90 deg rotated Gumbel copula family. Number 24 in CDVine.")
}

## density ##
setMethod("dcopula", signature("r90GumbelCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r90GumbelCopula"), linkCDVine.r90CDF)

## partial derivatives ##
# ddu
setMethod("dducopula", signature("r90GumbelCopula"), linkCDVine.ddu)

## ddv
setMethod("ddvcopula", signature("r90GumbelCopula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r90GumbelCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("r90GumbelCopula"),
          function(copula, tau) {
            if(tau > 0) warning("The rotated Gumbel copula can only represent non-positive dependence!")
            linkCDVine.calibKendallsTau(copula, min(0,-abs(tau)))
          })

########################
## Gumbel copula 270 ##
########################

setClass("r270GumbelCopula",
         representation = representation("copula", family="numeric"),
         validity = validRotGumbelCopula,
         contains = list("copula")
         )

# constructor
r270GumbelCopula <- function (param) {
  new("r270GumbelCopula", dimension = as.integer(2), parameters = param, param.names = c("theta", "delta"), 
      param.lowbnd = -Inf, param.upbnd = -1, family=34, 
      message = "270 deg rotated Gumbel copula family. Number 34 in CDVine.")
}

## density ##
setMethod("dcopula", signature("r270GumbelCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pcopula", signature("r270GumbelCopula"), linkCDVine.r270CDF)

## partial derivatives ##
# ddu
setMethod("dducopula", signature("r270GumbelCopula"), linkCDVine.ddu)

# ddv
setMethod("ddvcopula", signature("r270GumbelCopula"), linkCDVine.ddv)

## random number generator
setMethod("rcopula", signature("r270GumbelCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("calibKendallsTau", signature("r270GumbelCopula"), 
          function(copula, tau) {
            if(tau >= 0) warning("The rotated Gumbel copula can only represent negative dependence!")
            linkCDVine.calibKendallsTau(copula, min(-1e-6,-abs(tau)))
          })