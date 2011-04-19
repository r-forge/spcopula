#################################################################################
##
##   R package spcopula by Benedikt Gräler Copyright (C) 2011
##
##   This file is part of the R package spcopula.
##
##   The R package spcopula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package spcopula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################
## some additional bivariate copulas extending the set of copulas in the package copula

####
## an asymmetric copula with cubic and qudratic sections

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

####
## partial derivatives

setGeneric("dducopula", function(copula, pair) standardGeneric("dducopula"))
setGeneric("ddvcopula", function(copula, pair) standardGeneric("ddvcopula"))

## inverse partial derivatives 
setGeneric("invdducopula", function(copula, u, y) standardGeneric("invdducopula"))
setGeneric("invddvcopula", function(copula, v, y) standardGeneric("invddvcopula"))

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
  if (length(object@components) != length(object@distances)) return("Length of components + 1 does not equal length of distances. \n Note: The last distance must give the range and its automatically associated with the indepenence copula.")
  if (is.na(match(object@depFun(NULL),c("kendall","spearman","id","none")))) return("depFun(NULL) must return 'spearman', 'kendall' or 'id'.")
  else return(TRUE)
}

setClass("spCopula",
  representation = representation("copula", components="list", distances="numeric", unit="character", depFun="function"),
  validity = validSpCopula,
  contains = list("copula")
)

########################################
## spatial classes providing the data ##
########################################

## neighbourhood:

sizeLim <- 25 #  a constant
setSizeLim <- function(x) {
  env <- parent.env(environment())
  unlockBinding("neighbourLim",env)
  assign("neighbourLim", x,envir=env)
  lockBinding("neighbourLim",env)
}

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

