#################################################################################
##
##   R package spcopula by Benedikt Gräler Copyright (C) 2011
##
##   This file is part of the R package spcopula.
##
##   The R package spcopula is free software: you can redistribute it and/or 
##   modify it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package spcopula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package spcopula. If not, see <http://www.gnu.org/licenses/>.
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
# depFun="function"         a optional spatial dependence function providing 
#                             Kendalls tau or Spearman's rho to calib* or exact 
#                             parameters

stCopula <- function(components, distances, t.lags, stDepFun, unit="m", t.res="day") {
  spCopList <- list()
  getSpCop <- function(comp,dist,time) spCopula(comp, dist, 
                                                spDepFun=function(h) stDepFun(h,time), unit)
  for(i in 1:length(t.lags)){
    spCopList <- append(spCopList, getSpCop(components[[i]], distances[[i]], i))
  }
  
  param       <- unlist(lapply(spCopList, function(x) x@parameters))
  param.names <- unlist(lapply(spCopList, function(x) x@param.names))
  param.low   <- unlist(lapply(spCopList, function(x) x@param.lowbnd))
  param.up    <- unlist(lapply(spCopList, function(x) x@param.upbnd))
  
  new("stCopula", dimension=2, parameters=param, param.names=param.names,
      param.lowbnd=param.low, param.upbnd=param.up,
      message="Spatio-Temporal Copula: distance and time dependent convex combination of bivariate copulas",
      spCopList=spCopList, t.lags=t.lags, t.res=t.res)
}

## show method
showStCopula <- function(object) {
  cat(object@message, "\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Copulas:\n")
  for (i in 1:length(object@spCopList)) {
    cmpCop <- object@spCopList[[i]]
    cat("  ", cmpCop@message, "at", object@t.lags[i], 
      paste("[",object@t.res,"]",sep=""), "\n")
    show(cmpCop)
    
  }
}

setMethod("show", signature("stCopula"), showStCopula)

## spatial copula jcdf ##

# u 
#   list containing two column matrix providing the transformed pairs,  their respective 
#   separation distances and time steps
pStCopula <- function (copula, u) {
  if (!is.list(u) || !length(u)==3) stop("Point pairs need to be provided with their separating spatial and temproal distances as a list.")
  
  if(!is.matrix(u[[1]])) u[[1]] <- matrix(u[[1]],ncol=2)
  n <- nrow(u[[1]])
  h <- u[[2]]
  t.dist <- u[[3]]

  if(any(is.na(match(t.dist,copula@t.lags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  if(length(h)>1 & length(h)!=n) 
    stop("The spatial distance vector must either be of same length as rows in the data pairs or a single value.")
  if(length(t.dist)>1 & length(t.dist)!=n) 
    stop("The temporal distances vector must either be of same length as rows in the data pairs or a single value.")

  if (length(t.dist)==1) {
    res <- spcopula:::pSpCopula(copula@spCopList[[match(t.dist,copula@t.lags)]], 
                                list(u[[1]], u[[2]]))
  } else {
    if(length(h)==1) h <- rep(h,n)
    res <- NULL
    for(i in 1:n) {
      cat("fooFoo",match(t.dist[i],copula@t.lags),"\n")
      res <- rbind(res, spDepFunCopSnglDist(pcopula,
                                            copula@spCopList[[match(t.dist[i],copula@t.lags)]], 
                                            u[[1]][i,], h[i]))
    }
  }
  
  return(res)
}

setMethod("pcopula", signature("stCopula"), pStCopula)


## spatial Copula density ##

# u 
#   three column matrix providing the transformed pairs and their respective 
#   separation distances
dSpCopula <- function (copula, u) {
  if (!is.list(u) || !length(u)==2) stop("Point pairs need to be provided with their separating distance as a list.")

  h <- u[[2]]
  if(length(h)>1) {
    if(length(h)!=nrow(u[[1]])) stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
    ordering <- order(h)
    
    # ascending sorted pairs allow for easy evaluation
    pairs <- u[[1]][ordering,,drop=FALSE] 
    h <- h[ordering]
  } else {
    pairs <- u[[1]]
  }

  if(is.null(copula@calibMoa)) res <- spConCop(dcopula, copula, pairs, h)
  else {
    if(length(h)>1) {
      res <- spDepFunCop(dcopula, copula, pairs, h)
      
      # reordering the values
      res <- res[order(ordering)]
    } else {
      res <- spDepFunCopSnglDist(dcopula, copula, pairs, h)
    }
  }

return(res)
}

setMethod("dcopula", signature("spCopula"), dSpCopula)


## partial derivatives ##

## dduSpCopula
###############

dduSpCopula <- function (copula, pair) {
  if (!is.list(pair) || !length(pair)==2) stop("Point pairs need to be provided with their separating distance as a list.")
  
  h <- pair[[2]]
  if(length(h)>1) {
    if(length(h)!=nrow(pair[[1]])) stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
    ordering <- order(h)
    
    # ascending sorted pairs allow for easy evaluation
    pairs <- pair[[1]][ordering,,drop=FALSE] 
    h <- h[ordering]
  } else {
    pairs <- pair[[1]]
  }
  
  if(is.null(copula@calibMoa)) res <- spConCop(dducopula, copula, pairs, h)
  else {
    if(length(h)>1) {
      res <- spDepFunCop(dducopula, copula, pairs, h)
      
      # reordering the values
      res <- res[order(ordering)]
    } else {
      res <- spDepFunCopSnglDist(dducopula, copula, pairs, h)
    }
  }
  
  return(res)
}

setMethod("dducopula", signature("spCopula"), dduSpCopula)

## ddvSpCopula
###############

ddvSpCopula <- function (copula, pair) {
  if (!is.list(pair) || !length(pair)==2) stop("Point pairs need to be provided with their separating distance as a list.")
  
  h <- pair[[2]]
  if(length(h)>1) {
    if(length(h)!=nrow(pair[[1]])) stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
    ordering <- order(h) 
    
    # ascending sorted pairs allow for easy evaluation
    pairs <- pair[[1]][ordering,,drop=FALSE] 
    h <- h[ordering]
  } else {
    pairs <- pair[[1]]
  }
  
  if(is.null(copula@calibMoa)) res <- spConCop(ddvcopula, copula, pairs, h)
  else {
    if(length(h)>1) {
      res <- spDepFunCop(ddvcopula, copula, pairs, h)
      
      # reordering the values
      res <- res[order(ordering)]
    } else {
      res <- spDepFunCopSnglDist(ddvcopula, copula, pairs, h)
    }
  }
  
  return(res)
}

setMethod("ddvcopula", signature("spCopula"), ddvSpCopula)

#############
##         ##
## FITTING ##
##         ##
#############

# two models: 
# 1) Kendall's tau driven:
#    fit curve through emp. Kendall's tau values, identify validity ranges for
#    copula families deriving parameters from the fit, fade from one family to 
#    another at borders
# 2) convex-linear combination of copulas: 
#    fit one per lag, fade from one to another

# towards the first model:

# INPUT: the stBinning
# steps
# a) fit a curve
# b) estimate bivariate copulas per lag (limited to those with some 1-1-relation 
#    to Kendall's tau')
# INTERMEDIATE RESULT
# c) select best fits based on ... e.g. log-likelihood, visual inspection
# d) compose bivariate copulas to one spatial copula
# OUTPUT: a spatial copula parametrixued by distance through Kendall's tau

# towards a)
# bins   -> typically output from calcBins
# type   -> the type of curve (by now only polynominals are supported)
# degree -> the degree of the polynominal
# cutoff -> maximal distance that should be considered for fitting
# bounds -> the bounds of the correlation function (typically c(0,1))
# method -> the measure of association, either "kendall" or "spearman"
fitCorFun <- function(bins, type="poly", degree=3, cutoff=NA, bounds=c(0,1), method="kendall") {
  bins <- as.data.frame(bins[1:2])
  if(!is.na(cutoff)) bins <- bins[which(bins[[1]] <= cutoff),]
  
  fitCor <- lm(lagCor ~ poly(meanDists, degree), data = bins)
  
  print(fitCor)
  cat("Sum of squared residuals:",sum(fitCor$residuals^2),"\n")
  
  function(x) {
    if (is.null(x)) return(method)
    return(pmin(bounds[2], pmax(bounds[1], eval(predict(fitCor, data.frame(meanDists=x))))))
  }
}


# towards b)
# bins     -> typically output from calcBins
# calcTau  -> a function on distance providing Kendall's tau estimates, 
# families -> a vector of dummy copula objects of each family to be considered
#             DEFAULT: c(normal, t_df=4, clayton, frank, gumbel
loglikByCopulasLags <- function(bins, calcTau, families=c(normalCopula(0), tCopula(0,dispstr="un"),
                                                          claytonCopula(0), frankCopula(1), gumbelCopula(1))) {
  loglik <- NULL
  for (cop in families) {
    print(cop)
    tmploglik <- NULL
    for(i in 1:length(bins$meanDists)) {
      cop@parameters[1] <- calibKendallsTau(cop,tau=calcTau(bins$meanDists[i]))
      tmploglik <- c(tmploglik, sum(log(dcopula(cop,bins$lagData[[i]]))))
    }
    loglik <- cbind(loglik, tmploglik)
  }

  colnames(loglik) <- sapply(families, function(x) class(x)[1])

  return(loglik)
}

# towards d)
composeSpCop <- function(bestFit, families, bins, calcCor) {
  nfits <- length(bestFit)
  gaps <- which(diff(bestFit)!=0)

  if(missing(calcCor)) noCor <- nfits
  else noCor <- min(which(calcCor(bins$meanDists)<=0), nfits)
  
  breaks <- sort(c(gaps, gaps+1, noCor))
  breaks <- breaks[breaks<noCor]
  
  cops <- as.list(families[bestFit[breaks]])
  
  breaks <- unique(c(breaks, min(nfits,rev(breaks)[1]+1)))
  distances <- bins$meanDists[breaks]
  
  if(length(breaks)==length(cops)) {
    warning("Lags do not cover the entire range with positive correlation. ", 
             "An artifical one fading towards the indepCopula has been added.")
    distances <- c(distances, rev(distances)[1]+diff(bins$meanDists[nfits+c(-1,0)]))
  }

  if(missing(calcCor)) return(spCopula(components=cops, distances=distances, unit="m"))
  else return(spCopula(components=cops, distances=distances, unit="m", spDepFun=calcCor))
}

# in once

# bins   -> typically output from calcBins
# cutoff -> maximal distance that should be considered for fitting
# families -> a vector of dummy copula objects of each family to be considered
#             DEFAULT: c(normal, t_df=4, clayton, frank, gumbel
# ...
# type   -> the type of curve (by now only polynominals are supported)
# degree -> the degree of the polynominal
# bounds -> the bounds of the correlation function (typically c(0,1))
# method -> the measure of association, either "kendall" or "spearman"
fitSpCopula <- function(bins, cutoff=NA, families=c(normalCopula(0), tCopula(0,dispstr="un"),
                                                    claytonCopula(0), frankCopula(1), gumbelCopula(1)), ...) {
  calcTau <- fitCorFun(bins, cutoff=cutoff, ...)
  loglik <- loglikByCopulasLags(bins, calcTau=calcTau, families=families)
  
  bestFit <- apply(apply(loglik, 1, rank),2,function(x) which(x==length(families)))
  
  return(composeSpCop(bestFit, families, bins, calcTau))
}
