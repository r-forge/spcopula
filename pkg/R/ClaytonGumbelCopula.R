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
      fullname = "Survival Clayton copula family. Number 13 in CDVine.")
}

## density ##
setMethod("dCopula", signature("numeric","surClaytonCopula"), 
          function(u, copula, log) {
            linkCDVine.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","surClaytonCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","surClaytonCopula"), 
          function(u, copula, log) {
            linkCDVine.surCDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("pCopula", signature("matrix","surClaytonCopula"), linkCDVine.surCDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","surClaytonCopula"), 
          function(u, copula, log) {
            linkCDVine.ddu(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dduCopula", signature("matrix","surClaytonCopula"), linkCDVine.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","surClaytonCopula"), 
          function(u, copula, log) {
            linkCDVine.ddv(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("ddvCopula", signature("matrix","surClaytonCopula"), linkCDVine.ddv)

## random number generater ??
setMethod("rCopula", signature("numeric","surClaytonCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("surClaytonCopula"), 
          function(copula, tau) {
            if(tau <= 0) warning("The survival Clayton copula can only represent positive dependence!")
            linkCDVine.iTau(copula, max(1e-6,abs(tau)))
          })

setMethod("tau",signature("surClaytonCopula"),linkCDVine.tau)
setMethod("tailIndex",signature("surClaytonCopula"),linkCDVine.tailIndex)

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
      fullname = "90 deg rotated Clayton copula family. Number 23 in CDVine.")
}

## density ##
setMethod("dCopula", signature("numeric","r90ClaytonCopula"), 
          function(u, copula, log) {
            linkCDVine.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r90ClaytonCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r90ClaytonCopula"), 
          function(u, copula) {
            linkCDVine.r90CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r90ClaytonCopula"), linkCDVine.r90CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90ClaytonCopula"), 
          function(u, copula) {
            linkCDVine.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r90ClaytonCopula"), linkCDVine.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r90ClaytonCopula"), 
          function(u, copula) {
            linkCDVine.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r90ClaytonCopula"), linkCDVine.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90ClaytonCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("r90ClaytonCopula"),
          function(copula, tau) {
            if(tau >= 0) warning("The rotated Clayton copula can only represent negative dependence!")
            linkCDVine.iTau(copula, min(-1e-6,-abs(tau)))
          })

setMethod("tau",signature("r90ClaytonCopula"),linkCDVine.tau)

setMethod("tailIndex",signature("r90ClaytonCopula"),linkCDVine.tailIndex)

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
      fullname = "270 deg rotated Clayton copula family. Number 33 in CDVine.")
}

## density ##
setMethod("dCopula", signature("numeric","r270ClaytonCopula"), 
          function(u, copula, log) {
            linkCDVine.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r270ClaytonCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r270ClaytonCopula"), 
          function(u, copula) {
            linkCDVine.r270CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r270ClaytonCopula"), linkCDVine.r270CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270ClaytonCopula"), 
          function(u, copula) {
            linkCDVine.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r270ClaytonCopula"), linkCDVine.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r270ClaytonCopula"), 
          function(u, copula) {
            linkCDVine.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r270ClaytonCopula"), linkCDVine.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270ClaytonCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("r270ClaytonCopula"), 
          function(copula, tau) {
            if(tau >= 0) warning("The rotated Clayton copula can only represent negative dependence!")
            linkCDVine.iTau(copula, min(-1e-6,-abs(tau)))
          })

setMethod("tau",signature("r270ClaytonCopula"),linkCDVine.tau)

setMethod("tailIndex",signature("r270ClaytonCopula"),linkCDVine.tailIndex)

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
      fullname = "Survival Gumbel copula family. Number 14 in CDVine.")
}

## density ##
setMethod("dCopula", signature("numeric","surGumbelCopula"), 
          function(u, copula, log) {
            linkCDVine.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","surGumbelCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","surGumbelCopula"), 
          function(u, copula) {
            linkCDVine.surCDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","surGumbelCopula"), linkCDVine.surCDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","surGumbelCopula"), 
          function(u, copula) {
            linkCDVine.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","surGumbelCopula"), linkCDVine.ddu)

# ddv
setMethod("ddvCopula", signature("numeric","surGumbelCopula"), 
          function(u, copula) {
            linkCDVine.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","surGumbelCopula"), linkCDVine.ddv)

## random number generater ??
setMethod("rCopula", signature("numeric","surGumbelCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("surGumbelCopula"), 
          function(copula, tau) {
            if(tau < 0) warning("The survival Gumbel copula can only represent non-negative dependence!")
            linkCDVine.iTau(copula, max(0,abs(tau)))
          })

setMethod("tau",signature("surGumbelCopula"),linkCDVine.tau)

setMethod("tailIndex",signature("surGumbelCopula"),linkCDVine.tailIndex)

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
      fullname = "90 deg rotated Gumbel copula family. Number 24 in CDVine.")
}

## density ##
setMethod("dCopula", signature("numeric","r90GumbelCopula"), 
          function(u, copula, log) {
            linkCDVine.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r90GumbelCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r90GumbelCopula"), 
          function(u, copula) {
            linkCDVine.r90CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r90GumbelCopula"), linkCDVine.r90CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r90GumbelCopula"), 
          function(u, copula) {
            linkCDVine.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r90GumbelCopula"), linkCDVine.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r90GumbelCopula"), 
          function(u, copula) {
            linkCDVine.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r90GumbelCopula"), linkCDVine.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r90GumbelCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("r90GumbelCopula"),
          function(copula, tau) {
            if(tau > 0) warning("The rotated Gumbel copula can only represent non-positive dependence!")
            linkCDVine.iTau(copula, min(0,-abs(tau)))
          })

setMethod("tau",signature("r90GumbelCopula"),linkCDVine.tau)

setMethod("tailIndex",signature("r90GumbelCopula"),linkCDVine.tailIndex)

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
      fullname = "270 deg rotated Gumbel copula family. Number 34 in CDVine.")
}

## density ##
setMethod("dCopula", signature("numeric","r270GumbelCopula"), 
          function(u, copula, log) {
            linkCDVine.PDF(matrix(u,ncol=copula@dimension),copula, log)
          })
setMethod("dCopula", signature("matrix","r270GumbelCopula"), linkCDVine.PDF)

## jcdf ##
setMethod("pCopula", signature("numeric","r270GumbelCopula"), 
          function(u, copula) {
            linkCDVine.r270CDF(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix","r270GumbelCopula"), linkCDVine.r270CDF)

## partial derivatives ##
# ddu
setMethod("dduCopula", signature("numeric","r270GumbelCopula"), 
          function(u, copula) {
            linkCDVine.ddu(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","r270GumbelCopula"), linkCDVine.ddu)

## ddv
setMethod("ddvCopula", signature("numeric","r270GumbelCopula"), 
          function(u, copula) {
            linkCDVine.ddv(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","r270GumbelCopula"), linkCDVine.ddv)

## random number generator
setMethod("rCopula", signature("numeric","r270GumbelCopula"), linkCDVine.r)

## Kendalls tau to parameter conversion
setMethod("iTau", signature("r270GumbelCopula"), 
          function(copula, tau) {
            if(tau >= 0) warning("The rotated Gumbel copula can only represent negative dependence!")
            linkCDVine.iTau(copula, min(-1e-6,-abs(tau)))
          })

setMethod("tau",signature("r270GumbelCopula"),linkCDVine.tau)

setMethod("tailIndex",signature("r270GumbelCopula"),linkCDVine.tailIndex)
