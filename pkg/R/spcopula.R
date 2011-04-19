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

# constructor
# dimension = "numeric"     set to 2
# parameters = "numeric"    set of parameters
# param.names = "character" appropriate names
# param.lowbnd = "numeric"  appropriate lower bounds
# param.upbnd = "numeric"   appropriate upper bounds
# message = "character"     messgae printed with "show"
# components="list"         list of copulas (will be automatically supplemented 
#			      by the independent copula)
# distances="numeric"       the linking distances + the range (will be assigned
#			      to the independent copula)
# unit="character"          measurement unit of distance
# depFun="function"         a optinal spatial dependence function providing 
#                             Kendalls tau or Spearman's rho to calib* or exact 
#                             parameters


spCopula <- function(components, distances, unit="m", depFun=NULL) {
  if (is.null(depFun)){ depFun <- function(x) return("none")}
  else cat("The parameters of the components will be recalculated according to 
    the provided function.")

  components <- append(components,indepCopula())
  param <- NULL
  param.names <- NULL
  param.low <- NULL
  param.up<- NULL
  for (cmpCop in components) {
    param <- c(param, cmpCop@parameters)
    param.names <- c(param.names,cmpCop@param.names)
    param.low   <- c(param.low,cmpCop@param.lowbnd)
    param.up  <- c(param.up, cmpCop@param.upbnd)
  }
 
  new("spCopula", dimension=2, parameters=param, param.names=param.names,
      param.lowbnd=param.low, param.upbnd=param.up,
      message="Spatial Copula: distance dependent convex combination of 
	bivariate copulas",
      components=components, distances=distances, unit=unit, depFun=depFun)
}

## show method
showCopula <- function(object) {
  cat(object@message, "\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Copulas:\n")
  for (i in 1:length(object@components)) {
    cmpCop <- object@components[[i]]
    cat("  ", cmpCop@message, "at", object@distances[i], 
      paste("[",object@unit,"]",sep=""), "\n")
    if (length(cmpCop@parameters) > 0) {
      for (i in (1:length(cmpCop@parameters))) 
        cat("    ", cmpCop@param.names[i], " = ", cmpCop@parameters[i], "\n")
    }
  }
  if(object@depFun(NULL)!="none") cat("A spatial dependence function is used.")
}

setMethod("show", signature("spCopula"), showCopula)

## spatial copula jcdf ##

# u 
#   three column matrix providing the transformed pairs and their respective 
#   separation distances

pSpCopula <- function (copula, u) {
if (!is.matrix(u)) pair <- matrix(u,ncol=3)
if (ncol(u)!=3) stop("point pairs need top be provided with their separating 
		  distance")

# ascending sorted pairs allow for easy evaluation
pairs <- u[order(u[,3]),1:2,drop=FALSE] 
depFun <- copula@depFun
method <- depFun(NULL)

h <- sort(u[,3])
n.dists <- length(distances)
res <- numeric(0)
if (method != "none"){
  calibMoa <- switch(method, kendall=calibKendallsTau,
                 spearmann=calibSpearmansRho, 
                 id=function(copula, tau) return(tau))
}

sel <- which(h < distances[1])
if(sum(sel)>0) {
  tmpH <- h[sel]
  tmpCop <- copula@components[[1]]
  tmpPairs <- pairs[sel,,drop=FALSE]
  if (method != "none"){
    for (j in 1:length(tmpH)) {
      tmpCop@parameters <- calibMoa(tmpCop,depFun(tmpH[j]))
      res <- c(res, pcopula(tmpCop,tmpPairs[j,]))
    }
  } else res <- pcopula(tmpCop,tmpPairs)
}

if (n.dists > 2) {
for ( i in 2:n.dists ) {
  low <- distances[i-1]
  high<- distances[i]
  sel <- which(h >= low & h < high)
  if (sum(sel)>0) {
    tmpH <- h[sel]
    tmpPairs <- pairs[sel,,drop=FALSE]
    lowerCop <- copula@components[[i-1]]
    upperCop <- copula@components[[i]]
    if (method != "none") {
      lowerVals <- numeric(0)
      upperVals <- numeric(0)
      for (j in 1:length(tmpH)) {
        lowerCop@parameters <- calibMoa(lowerCop,depFun(tmpH[j]))
        upperCop@parameters <- calibMoa(upperCop,depFun(tmpH[j]))
        lowerVals <- c(lowerVals, pcopula(lowerCop,tmpPairs[j,]))
        upperVals <- c(upperVals, pcopula(upperCop,tmpPairs[j,]))
      }
      res <- c( res,
        (high-tmpH)/(high-low) * lowerVals
        +(tmpH-low)/(high-low) * upperVals)
    } else {
      res <- c( res,
        (high-tmpH)/(high-low) * pcopula(lowerCop,tmpPairs)
        +(tmpH-low)/(high-low) * pcopula(upperCop,tmpPairs))
    }
  }
}}

sel <- which(h >= distances[n.dists])
if(sum(sel)>0) {
  res <- c(res,pcopula(copula@components[[n.dists]],
	       pairs[which(h >= distances[n.dists]),]))
}

# reordering the values
return(res[order(u[,3])])
}

setMethod("pcopula", signature("spCopula"), pSpCopula)

## spatial Copula density ##
# the non dynamic version: without depFun
dSpCopula <- function (copula, u) {
if (!is.matrix(u)) pair <- matrix(u,ncol=3)
if (ncol(u)!=3) stop("point pairs need top be provided with their separating 
		  distance")

# ascending sorted pairs allow for easy evaluation
pairs <- u[order(u[,3]),1:2,drop=FALSE] 
depFun <- copula@depFun
method <- depFun(NULL)

if (method != "none"){
  calibMoa <- switch(method, kendall=calibKendallsTau,
                 spearmann=calibSpearmansRho, 
                 id=function(copula, tau) return(tau))
}

h <- sort(u[,3])
n.dists <- length(distances)
res <- numeric(0)

sel <- which(h < distances[1])
if(sum(sel)>0) {
  tmpH <- h[sel]
  tmpCop <- copula@components[[1]]
  tmpPairs <- pairs[sel,,drop=FALSE]
  if (method != "none"){
    for (j in 1:length(tmpH)) {
      tmpCop@parameters <- calibMoa(tmpCop,depFun(tmpH[j]))
      res <- c(res, dcopula(tmpCop,tmpPairs[j,]))
    }
  } else res <- dcopula(tmpCop,tmpPairs)
}

if (n.dists > 2) {
for ( i in 2:n.dists ) {
  low <- distances[i-1]
  high<- distances[i]
  sel <- which(h >= low & h < high)
  if (sum(sel)>0) {
    tmpH <- h[sel]
    tmpPairs <- pairs[sel,,drop=FALSE]
    lowerCop <- copula@components[[i-1]]
    upperCop <- copula@components[[i]]
    if (method != "none") {
      lowerVals <- numeric(0)
      upperVals <- numeric(0)
      for (j in 1:length(tmpH)) {
        lowerCop@parameters <- calibMoa(lowerCop,depFun(tmpH[j]))
        upperCop@parameters <- calibMoa(upperCop,depFun(tmpH[j]))
        lowerVals <- c(lowerVals, dcopula(lowerCop,tmpPairs[j,]))
        upperVals <- c(upperVals, dcopula(upperCop,tmpPairs[j,]))
      }
      res <- c( res,
        (high-tmpH)/(high-low) * lowerVals
        +(tmpH-low)/(high-low) * upperVals)
    } else {
      res <- c( res,
        (high-tmpH)/(high-low) * dcopula(lowerCop,tmpPairs)
        +(tmpH-low)/(high-low) * dcopula(upperCop,tmpPairs))
    }
  }
}}

sel <- which(h >= distances[n.dists])
if(sum(sel)>0) {
  res <- c(res,dcopula(copula@components[[n.dists]],
	       pairs[which(h >= distances[n.dists]),]))
}

# reordering the values
return(res[order(u[,3])])
}

setMethod("dcopula", signature("spCopula"), dSpCopula)

## partial derivatives ##

## dduSpCopula
###############

dduSpCopula <- function (copula, pair) {
if (!is.matrix(pair)) pair <- matrix(pair,ncol=3)
if (ncol(pair)!=3) stop("point pairs need top be provided with their separating 
		  distance")

# ascending sorted pairs allow for easy evaluation
pairs <- pair[order(pair[,3]),1:2,drop=FALSE] 
depFun <- copula@depFun
method <- depFun(NULL)

if (method != "none"){
  calibMoa <- switch(method, kendall=calibKendallsTau,
                 spearmann=calibSpearmansRho, 
                 id=function(copula, tau) return(tau))
}

h <- sort(pair[,3])
n.dists <- length(distances)
res <- numeric(0)

sel <- which(h < distances[1])
if(sum(sel)>0) {
  tmpH <- h[sel]
  tmpCop <- copula@components[[1]]
  tmpPairs <- pairs[sel,,drop=FALSE]
  if (method != "none"){
    for (j in 1:length(tmpH)) {
      tmpCop@parameters <- calibMoa(tmpCop,depFun(tmpH[j]))
      res <- c(res, dducopula(tmpCop,tmpPairs[j,]))
    }
  } else res <- dducopula(tmpCop,tmpPairs)
}

if (n.dists > 2) {
for ( i in 2:n.dists ) {
  low <- distances[i-1]
  high<- distances[i]
  sel <- which(h >= low & h < high)
  if (sum(sel)>0) {
    tmpH <- h[sel]
    tmpPairs <- pairs[sel,,drop=FALSE]
    lowerCop <- copula@components[[i-1]]
    upperCop <- copula@components[[i]]
    if (method != "none") {
      lowerVals <- numeric(0)
      upperVals <- numeric(0)
      for (j in 1:length(tmpH)) {
        lowerCop@parameters <- calibMoa(lowerCop,depFun(tmpH[j]))
        upperCop@parameters <- calibMoa(upperCop,depFun(tmpH[j]))
        lowerVals <- c(lowerVals, dducopula(lowerCop,tmpPairs[j,]))
        upperVals <- c(upperVals, dducopula(upperCop,tmpPairs[j,]))
      }
      res <- c( res,
        (high-tmpH)/(high-low) * lowerVals
        +(tmpH-low)/(high-low) * upperVals)
    } else {
      res <- c( res,
        (high-tmpH)/(high-low) * dducopula(lowerCop,tmpPairs)
        +(tmpH-low)/(high-low) * dducopula(upperCop,tmpPairs))
    }
  }
}}

sel <- which(h >= distances[n.dists])
if(sum(sel)>0) {
  res <- c(res,dducopula(copula@components[[n.dists]],
	       pairs[which(h >= distances[n.dists]),]))
}

# reordering the values
return(res[order(pair[,3])])
}

setMethod("dducopula", signature("spCopula"), dduSpCopula)

## ddvSpCopula
###############

ddvSpCopula <- function (copula, pair) {
if (!is.matrix(pair)) pair <- matrix(pair,ncol=3)
if (ncol(pair)!=3) stop("point pairs need top be provided with their separating 
		  distance")

# ascending sorted pairs allow for easy evaluation
pairs <- pair[order(pair[,3]),1:2,drop=FALSE] 
depFun <- copula@depFun
method <- depFun(NULL)

if (method != "none"){
  calibMoa <- switch(method, kendall=calibKendallsTau,
                 spearmann=calibSpearmansRho, 
                 id=function(copula, tau) return(tau))
}

h <- sort(pair[,3])
n.dists <- length(distances)
res <- numeric(0)

sel <- which(h < distances[1])
if(sum(sel)>0) {
  tmpH <- h[sel]
  tmpCop <- copula@components[[1]]
  tmpPairs <- pairs[sel,,drop=FALSE]
  if (method != "none"){
    for (j in 1:length(tmpH)) {
      tmpCop@parameters <- calibMoa(tmpCop,depFun(tmpH[j]))
      res <- c(res, ddvcopula(tmpCop,tmpPairs[j,]))
    }
  } else res <- ddvcopula(tmpCop,tmpPairs)
}

if (n.dists > 2) {
for ( i in 2:n.dists ) {
  low <- distances[i-1]
  high<- distances[i]
  sel <- which(h >= low & h < high)
  if (sum(sel)>0) {
    tmpH <- h[sel]
    tmpPairs <- pairs[sel,,drop=FALSE]
    lowerCop <- copula@components[[i-1]]
    upperCop <- copula@components[[i]]
    if (method != "none") {
      lowerVals <- numeric(0)
      upperVals <- numeric(0)
      for (j in 1:length(tmpH)) {
        lowerCop@parameters <- calibMoa(lowerCop,depFun(tmpH[j]))
        upperCop@parameters <- calibMoa(upperCop,depFun(tmpH[j]))
        lowerVals <- c(lowerVals, ddvcopula(lowerCop,tmpPairs[j,]))
        upperVals <- c(upperVals, ddvcopula(upperCop,tmpPairs[j,]))
      }
      res <- c( res,
        (high-tmpH)/(high-low) * lowerVals
        +(tmpH-low)/(high-low) * upperVals)
    } else {
      res <- c( res,
        (high-tmpH)/(high-low) * ddvcopula(lowerCop,tmpPairs)
        +(tmpH-low)/(high-low) * ddvcopula(upperCop,tmpPairs))
    }
  }
}}

sel <- which(h >= distances[n.dists])
if(sum(sel)>0) {
  res <- c(res,ddvcopula(copula@components[[n.dists]],
	       pairs[which(h >= distances[n.dists]),]))
}

# reordering the values
return(res[order(pair[,3])])
}

setMethod("ddvcopula", signature("spCopula"), ddvSpCopula)

##
## 

## testing ##

# spCop <- spCopula(
# 	   list(normalCopula(runif(1,min=-1,max=1)),
# 	     normalCopula(runif(1,min=-1,max=1)),
# 	     normalCopula(runif(1,min=-1,max=1)),
# 	     normalCopula(runif(1,min=-1,max=1)) ),
# 	   distances=c(200,300,500,750,1000),depFun= function(x) {if(is.null(x))return("id");return(.4)}) # NULL)#  function(x) {if(is.null(x))return("kendall");return(.4)})
# 
# pcopula(spCop,matrix(c(.03, 1, 1000, 1, .4, 600, 0, 0, 400, 1, 1, 200),ncol=3,byrow=T))
# dcopula(spCop,matrix(c(.03, .04, 100, .03, .04, 400, .03, .4, 250, .03, .4, 600),ncol=3,byrow=T))
# dducopula(spCop,matrix(c(.03, .04, 100, .03, .04, 400, .03, .4, 250, .03, .4, 600),ncol=3,byrow=T))
# ddvcopula(spCop,matrix(c(.03, .04, 100, .03, .04, 400, .03, .4, 250, .03, .4, 600),ncol=3,byrow=T))
# 
# foo <- normalCopula(.5)
# foo@parameters <- .3
# 
# system.time(
# sum(dcopula(spCop,cbind(lrgUnitSq,rep(120,10000))))/10000
# )