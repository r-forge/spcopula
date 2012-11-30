## an asymmetric copula with cubic and quadratic sections

validAsCopula = function(object) {
  if (object@dimension != 2)
    return("Only copulas with cubic quadratic sections of dimension 2 are supported.")
  param <- object@parameters
  upper <- object@param.upbnd
  lower <- object@param.lowbnd
  if (length(param) != length(upper))
    return("Parameter and upper bound have non-equal length")
  if (length(param) != length(lower))
    return("Parameter and lower bound have non-equal length")
  if (any(is.na(param) | param > upper | param < lower))
    return("Parameter value out of bound")
  else return (TRUE)
}

# the lower bound of the parameter a dependening on the parameter b
limA <- function (b) {
  (b-3-sqrt(9+6*b-3*b^2))/2
}

setClass("asCopula",
  representation = representation("copula"),
  validity = validAsCopula,
  contains = list("copula")
)

####
## a symmetric copula with cubic and quadratic sections

validCqsCopula <- validAsCopula

setClass("cqsCopula",
  representation = representation("copula"),
  validity = validCqsCopula,
  contains = list("copula")
)

## 
## the spatial copula
##
## realized as a distance dependent convex combination of biv copulas

# dimension = "numeric"     set to 2
# parameters = "numeric"    set of parameters
# param.names = "character" appropriate names
# param.lowbnd = "numeric"  appropriate lower bounds
# param.upbnd = "numeric"   appropriate upper bounds
# message = "character"     messgae printed with "show"
# components="list"         list of copulas 
# distances="numeric"       the linking distances
# unit="character"          measurement unit of distance
# depFun="function"         an optional dependence function; depFun(NULL)
#                             has to return either "spearman" or "kendall" 
#                             dependening on the moa used. Make sure depFun
#                             assings valid parameters to the copulas involved

validSpCopula <- function(object) {
  if (length(object@components) != length(object@distances)) return("Length of components + 1 does not equal length of distances. \n Note: The last distance must give the range and it is automatically associated with the indepenence copula.")
  check.upper <- NULL
  check.lower <- NULL
  
  if(!is.null(object@calibMoa(normalCopula(0),0))) {
    for (i in 1:(length(object@components)-1)) {
      check.upper <- c(check.upper, is.na(object@calibMoa(object@components[[i]], object@distances[i+1])))
      check.lower <- c(check.lower, is.na(object@calibMoa(object@components[[i]], c(0,object@distances)[i])))
    }
    if(sum(check.upper>0)) return(paste("Reconsider the upper boundary conditions of the following copula(s): \n",
                                        paste(sapply(object@components[check.upper], function(x) x@message), 
                                              "at", object@distances[check.upper],collapse="\n")))
    if(sum(check.lower>0)) return(paste("Reconsider the lower boundary conditions of the following copula(s): \n",
                                        paste(sapply(object@components[check.lower], function(x) x@message), 
                                              "at", object@distances[check.lower],collapse="\n")))
  }
  return(TRUE)
}

setClass("spCopula", representation = representation("copula", 
                                                     components="list",
                                                     distances="numeric", 
                                                     calibMoa="function", 
                                                     unit="character"),
         validity = validSpCopula, contains = list("copula"))

############################
## Spatio-Temporal Copula ##
############################

validStCopula <- function(object) {
  if(length(object@t.lags) != length(object@spCopList)) return("The length of the temporal distance vector must equal the number of spatial copulas.")
  return(TRUE) # validity of any spCopula in spCopList is tested by the constructor, I believe
}

setClass("stCopula", representation = representation("copula", 
                                                     spCopList="list", 
                                                     t.lags="numeric",
                                                     t.res="character"),
         validity = validStCopula, contains = list("copula"))



########################################
## spatial classes providing the data ##
########################################

## neighbourhood:

sizeLim <- 25 #  a constant
# setSizeLim <- function(x) {
#   env <- parent.env(environment())
#   unlockBinding("neighbourLim",env)
#   assign("neighbourLim", x,envir=env)
#   lockBinding("neighbourLim",env)
# }

# a class combining two matrices holding the data and the corresponding 
# distances as well a slot for the coordinates refernce system and an attribute
# if the data is already transformed to uniform on [0,1] distributed variables
# data:		a list of data.frames holding the data per neighbour. each neighbour needs to have the same number of variables in the same order
# sp: an optional slot providing the coordinates of locations
# index: a matrix linking the data entries with the coordinates of the locations
validNeighbourhood <- function(object) {
  sizeN <- ncol(object@distances)+1
  nVars <- length(object@varNames)
  if (sizeN > sizeLim) return("The limting size of the neighbourhood is exceeded. Increase the constant sizeLim if needed.")
  if (nrow(object@data) != nrow(object@distances)) return("Data and distances have unequal number of rows.")
  if (ncol(object@data) %% sizeN != 0) return("Data and distances have non matching number of columns.")
  if (nrow(object@data) != nrow(object@coords) ) return("Data and sp@coordinates have unequal number of rows.")
  if (nrow(object@data) != nrow(object@index)) return("Data and index have unequal number of rows.")
  if (sizeN != ncol(object@index)) return("Data and index have unequal number of columns.")
  if (ncol(object@data) != sizeN * nVars) return(paste("Number of columns in data does not equal the product of the neighbourhood's size (",sizeN,") with number of variables (",nVars,").",sep=""))
  else return(TRUE)
}

setClass("neighbourhood",
  representation = representation(data = "data.frame", distances="matrix", "SpatialPoints", index="matrix", varNames="character"),
  validity = validNeighbourhood,
  contains = list("SpatialPoints"))

