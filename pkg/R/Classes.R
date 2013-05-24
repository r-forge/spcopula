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
  stopifnot(abs(b) <= 1)
  0.5*(-sqrt(-3*b^2+6*b+9)+b-3)
}

# the lower and upper bound of the parameter b dependening on the parameter a
limB <- function (a) {
  stopifnot(a <=1 & a >= -3)
  if(a>-2)
    return(c(-1,1))
  pmax(pmin(0.5*(c(-1,1)*(sqrt(3)*sqrt(-a^2-2*a+3))+a+3),1),-1)
}

setClass("asCopula",
  representation = representation("copula"),
  validity = validAsCopula,
  contains = list("copula")
)

####
## a symmetric copula with cubic and quadratic sections

validCqsCopula <- function(object) {
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
  if (object@fixed != ""){
    if(!("a" %in% object@fixed | "b" %in% object@fixed))
      return("The slot fixed may only refer to \"a\" or \"b\".")
    if ("a" %in% object@fixed & "b" %in% object@fixed)
      return("Only one of the parameters may be kept fixed.")
  }
  else return (TRUE)
}

setClass("cqsCopula",
  representation = representation("copula",fixed="character"),
  validity = validCqsCopula,
  contains = list("copula")
)

####
## an empirical copula representation

validEmpCopula <- function(object) {
  if(ncol(object@sample) != object@dimension)
    return("Dimension of the copula and the sample do not match.")
  else
    return(TRUE)
}

setClass("empiricalCopula",
         representation = representation("copula", sample="matrix"),
         validity = validEmpCopula,
         contains = list("copula")
)

####
## the leaf copula

validLeafCopula <- function(object) {
  if (object@dimension != 2)
    return("The leaf copula only supports two dimensions.")
  
  if (any(is.na(object@parameters)))
    return("Parameter value is \"NA\".")
  else return (TRUE)
}

setClass("leafCopula",
         representation = representation("copula"),
         validity = validLeafCopula,
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
  if (length(object@components) != length(object@distances)) 
    return("Length of components does not equal length of distances. \n Note: The last distance must give the range and it is automatically associated with the indepenence copula.")
  check.upper <- NULL
  check.lower <- NULL
  
  nComp <- length(object@components)
  if(!is.null(object@calibMoa(normalCopula(0),0))) {
    nonIndep <- sapply(object@components[-nComp], function(x) class(x) != "indepCopula")
    for (i in (1:(nComp-1))[nonIndep]) {
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

####################
##  vine copulas  ##
####################

validVineCopula = function(object) {
  dim <- object@dimension
  if( dim <= 2)
    return("Number of dimension too small (>2).")
  if(length(object@copulas)!=(dim*(dim-1)/2))
    return("Number of provided copulas does not match given dimension.")
  if(!any(unlist(lapply(object@copulas,function(x) is(x,"copula")))))
    return("Not all provided copulas in your list are indeed copulas.")
  return (TRUE)
}

setOldClass("RVineMatrix")

setClass("vineCopula",
         representation = representation(copulas="list", dimension="integer", 
                                         RVM="RVineMatrix"),
         prototype = prototype(RVM=structure(list(),class="RVineMatrix")),
         validity = validVineCopula,
         contains = list("copula")
)

#########################
## Spatial Vine Copula ##
#########################

validMixedSpVineCopula <- function(object) {
  return(all(sapply(object@spCop,validSpCopula) & validObject(object@topCop)))
}

setClass("mixedSpVineCopula", representation("copula", spCop="list", topCop="copula"),
         validity = validMixedSpVineCopula, contains=list("copula"))

validPureSpVineCopula <- function(object) {
  return(all(sapply(object@spCop,validSpCopula)))
}

setClass("pureSpVineCopula", representation("copula", spCop="list"),
         validity = validPureSpVineCopula, contains=list("copula"))

setClassUnion("spVineCopula",c("mixedSpVineCopula","pureSpVineCopula"))

#################################
## Spatio-temporal Vine Copula ##
#################################

validStVineCopula <- function(object) {
  return(validStCopula(object@stCop) & validObject(object@topCop))
}

setClass("stVineCopula", representation("copula", stCop="stCopula", topCop="copula"),
         validity = validStVineCopula, contains=list("copula"))

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
  nVars <- length(object@var)
  if (object@prediction & is.null(object@dataLocs))
    return("The locations of the data have to provided for the estimation procedure.")
  if (nrow(object@data) != nrow(object@distances)) 
    return("Data and distances have unequal number of rows.")
  if (ncol(object@data) %% (sizeN-object@prediction) != 0) 
    return("Data and distances have non matching number of columns.")
  if (nrow(object@data) != nrow(object@index)) 
    return("Data and index have unequal number of rows.")
  if (ncol(object@distances) != ncol(object@index)) 
    return("Data and index have unequal number of columns.")
  if (ncol(object@data) != (sizeN-object@prediction) * nVars) 
    return(paste("Number of columns in data does not equal the product of the neighbourhood's size (",sizeN,") with number of variables (",nVars,").",sep=""))
  else 
    return(TRUE)
}

setClassUnion("optionalDataLocs",c("NULL","Spatial"))

setClass("neighbourhood",
         representation = representation(data = "data.frame", 
                                         distances="matrix", 
                                         index="matrix",
                                         locations="Spatial",
                                         dataLocs="optionalDataLocs",
                                         var="character", 
                                         prediction="logical"),
         validity = validNeighbourhood, contains = list("Spatial"))

## ST neighbourhood

validStNeighbourhood <- function(object) {
  sizeN <- nrow(object@data)
  if (object@prediction & is.null(object@dataLocs))
    return("The spatio-temporal locations of the data have to be provided for the estimation procedure.")
  dimDists <- dim(object@distances)
  if (nrow(object@data) != dimDists[1]) 
    return("Data and distances have unequal number of rows.")
  dimInd <- dim(object@index)
  if (nrow(object@data) != dimInd[1]) 
    return("Data and index have unequal number of rows.")
  if (dimDists[2] != dimInd[2]) 
    return("Data and index have unequal number of columns.")
  else 
    return(TRUE)
}

setClassUnion("optionalST",c("NULL","ST"))

setClass("stNeighbourhood",
         representation = representation(data = "data.frame", 
                                         distances="array", 
                                         index="array",
                                         locations="ST",
                                         dataLocs="optionalST",
                                         var="character", 
                                         prediction="logical"),
         validity = validStNeighbourhood, contains = list("ST"))