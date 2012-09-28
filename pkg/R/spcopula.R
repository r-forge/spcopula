#################################################################################
##
##  R package spCopula by Benedikt Gräler Copyright (C) 2011
##
##  This file is part of the R package spCopula.
##
##  The R package spCopula is free software: you can redistribute it and/or 
##  modify it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  The R package spCopula is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with the R package spCopula. If not, see <http://www.gnu.org/licenses/>
##
#################################################################################

# constructor
# dimension = "integer"     set to 2
# parameters = "numeric"    set of parameters
# param.names = "character" appropriate names
# param.lowbnd = "numeric"  appropriate lower bounds
# param.upbnd = "numeric"   appropriate upper bounds
# fullname = "character"     messgae printed with "show"
# components="list"         list of copulas (will be automatically supplemented 
#			      by the independent copula)
# distances="numeric"       the linking distances + the range (will be assigned
#			      to the independent copula)
# unit="character"          measurement unit of distance
# depFun="function"         a optional spatial dependence function providing 
#                             Kendalls tau or Spearman's rho to calib* or exact 
#                             parameters


spCopula <- function(components, distances, spDepFun, unit="m") {
  if (missing(spDepFun)) { 
    calibMoa <- function(copula, h) return(NULL)
  } else {
    if (is.na(match(spDepFun(NULL),c("kendall","spearman","id","none")))) stop("spDepFun(NULL) must return 'spearman', 'kendall' or 'id'.")
    cat("The parameters of the components will be recalculated according to the provided spDepFun. \n")
    calibMoa <- switch(spDepFun(NULL), 
                       kendall=function(copula, h) iTau(copula, spDepFun(h)),
                       spearman=function(copula, h) iRho(copula, spDepFun(h)),
                       id=function(copula, h) return(h))
    
    for (i in 1:length(components)) {
      components[[i]]@parameters[1] <- calibMoa(components[[i]], distances[i]) # take care of non single parameter copulas
    }
  }

  components <- append(components,indepCopula())
  
  param       <- unlist(lapply(components, function(x) x@parameters))
  param.names <- unlist(lapply(components, function(x) x@param.names))
  param.low   <- unlist(lapply(components, function(x) x@param.lowbnd))
  param.up    <- unlist(lapply(components, function(x) x@param.upbnd))
     
  new("spCopula", dimension=as.integer(2), parameters=param, param.names=param.names,
      param.lowbnd=param.low, param.upbnd=param.up,
      fullname="Spatial Copula: distance dependent convex combination of bivariate copulas",
      components=components, distances=distances, calibMoa=calibMoa, unit=unit)
}

## show method
showCopula <- function(object) {
  cat(object@fullname, "\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Copulas:\n")
  for (i in 1:length(object@components)) {
    cmpCop <- object@components[[i]]
    cat("  ", cmpCop@fullname, "at", object@distances[i], 
      paste("[",object@unit,"]",sep=""), "\n")
    if (length(cmpCop@parameters) > 0) {
      for (i in (1:length(cmpCop@parameters))) 
        cat("    ", cmpCop@param.names[i], " = ", cmpCop@parameters[i], "\n")
    }
  }
  if(!is.null(object@calibMoa(normalCopula(0),0))) cat("A spatial dependence function is used. \n")
}

setMethod("show", signature("spCopula"), showCopula)

## spatial copula jcdf ##


# for spatial copulas with a spatial dependece function
spDepFunCop <- function(fun, copula, pairs, h) {
  dists <- copula@distances
  n.dists <- length(dists)
  calibPar <- copula@calibMoa
  
  res <- numeric(0)
  sel <- which(h < dists[1])
  if(sum(sel)>0) {
    tmpH <- h[sel]
    tmpCop <- copula@components[[1]]
    tmpPairs <- pairs[sel,,drop=FALSE]
    for (j in 1:length(tmpH)) {
      tmpCop@parameters[1] <- calibPar(tmpCop, tmpH[j])
      res <- c(res, fun(tmpCop, tmpPairs[j,]))
    }
  }
  
  if (n.dists >= 2) {
    for ( i in 2:n.dists ) {
      low  <- dists[i-1]
      high <- dists[i]
      sel <- which(h >= low & h < high)
      if (sum(sel)>0) {
        tmpH <- h[sel]
        tmpPairs <- pairs[sel,,drop=FALSE]
        lowerCop <- copula@components[[i-1]]
        upperCop <- copula@components[[i]]
        if (class(lowerCop) != class(upperCop)) {
          lowerVals <- numeric(0)
          upperVals <- numeric(0)
          for (j in 1:length(tmpH)) {
            lowerCop@parameters[1] <- calibPar(lowerCop,  tmpH[j])
            upperCop@parameters[1] <- calibPar(upperCop, tmpH[j])
            lowerVals <- c(lowerVals, fun(lowerCop, tmpPairs[j,]))
            upperVals <- c(upperVals, fun(upperCop, tmpPairs[j,]))
          }
          res <- c(res,(high-tmpH)/(high-low)*lowerVals+(tmpH-low)/(high-low)*upperVals)
        } else {
          newVals <- numeric(0)
          for (j in 1:length(tmpH)) {
            lowerCop@parameters <- calibPar(lowerCop, tmpH[j])
            newVals <- c(newVals, fun(lowerCop, tmpPairs[j,]))
          }
          res <- c(res, newVals)
        }
      }
    }
  }
  
  sel <- which(h >= dists[n.dists])
  if(sum(sel)>0) {
    res <- c(res, fun(copula@components[[n.dists]], 
                      pairs[which(h >= dists[n.dists]),]))
  }

  return(res)
}

# for spatial copulas with a spatial dependece function and a single distance but many pairs
spDepFunCopSnglDist <- function(fun, copula, pairs, h) {
  dists <- copula@distances
  n.dists <- length(dists)
  calibPar <- copula@calibMoa

  if(h < dists[1]) {
    tmpCop <- copula@components[[1]]
    tmpCop@parameters[1] <- calibPar(tmpCop, h)
    res <- fun(pairs, tmpCop)
  }
  
  if (n.dists >= 2) {
    for ( i in 2:n.dists ) {
      low  <- dists[i-1]
      high <- dists[i]
      if (h >= low & h < high) {
        lowerCop <- copula@components[[i-1]]
        upperCop <- copula@components[[i]]
        if (class(lowerCop) != class(upperCop)) {
          lowerCop@parameters[1] <- calibPar(lowerCop, h)
          upperCop@parameters[1] <- calibPar(upperCop, h)

          lowerVals <- fun(pairs, lowerCop)
          upperVals <- fun(pairs, upperCop)

          res <- (high-h)/(high-low)*lowerVals + (h-low)/(high-low)*upperVals
        } else {
          lowerCop@parameters <- calibPar(lowerCop, h)
          res <- fun(pairs, lowerCop)
        }
      }
    }
  }
  
  if(h >= dists[n.dists]) {
    res <- fun(pairs, copula@components[[n.dists]])
  }
  
  return(res)
}


# for static convex combinations of copulas
spConCop <- function(fun, copula, pairs, h) {
  dists <- copula@distances
  n.dists <- length(dists)
  
  res <- numeric(nrow(pairs))
  sel <- which(h < dists[1])
  if(sum(sel)>0) {
    res[sel] <- fun(pairs[sel,,drop=FALSE],copula@components[[1]])
  }
  
  if (n.dists >= 2) {
    for ( i in 2:n.dists ) {
      low  <- dists[i-1]
      high <- dists[i]
      sel <- which(h >= low & h < high)
      if (sum(sel)>0) {
        tmpH <- h[sel]
        tmpPairs <- pairs[sel,,drop=FALSE]

        lowerVals <- fun(tmpPairs[,], copula@components[[i-1]])
        upperVals <- fun(tmpPairs[,], copula@components[[i]])

        res[sel] <- (high-tmpH)/(high-low)*lowerVals+(tmpH-low)/(high-low)*upperVals
      }
    }
  }
  
  sel <- which(h >= dists[n.dists])
  if(sum(sel)>0) {
    res[sel] <- fun(pairs[which(h >= dists[n.dists]),],
                    copula@components[[n.dists]])
  }

  return(res)
}


# u 
#   three column matrix providing the transformed pairs and their respective 
#   separation distances
pSpCopula <- function (u, copula) {
  if (!is.list(u) || !length(u)>=2) stop("Point pairs need to be provided with their separating distance as a list.")
  
  pairs <- u[[1]]
  if(!is.matrix(pairs)) pairs <- matrix(pairs,ncol=2)
  n <- nrow(pairs)
  
  if(length(u)==3) {
    block <- u[[3]]
    if (n%%block != 0) stop("The block size is not a multiple of the data length:",n)
  } else block <- 1
  
  h <- u[[2]]
  if(length(h)>1 && length(h)!=nrow(u[[1]])) {
    stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
  }
  
  
  if(is.null(copula@calibMoa(normalCopula(0),0))) {
    res <- spConCop(pCopula, copula, pairs, rep(h,length.out=nrow(pairs)))
  } else {
    if(length(h)>1) {
      if (block == 1){
        ordering <- order(h)
        
        # ascending sorted pairs allow for easy evaluation
        pairs <- pairs[ordering,,drop=FALSE] 
        h <- h[ordering]
        
        res <- spDepFunCop(pCopula, copula, pairs, h)
        
        # reordering the values
        res <- res[order(ordering)]
      } else {
        res <- NULL
        for(i in 1:(n%/%block)) {
          res <- c(res, spDepFunCopSnglDist(pCopula, copula, 
                                            pairs[((i-1)*block+1):(i*block),], 
                                            h[i*block]))
        }
      }
    } else {
      res <- spDepFunCopSnglDist(pCopula, copula, pairs, h)
    }
  }
  
  return(res)
}

setMethod("pCopula", signature("list","spCopula"), pSpCopula)

## spatial Copula density ##

# u 
#   three column matrix providing the transformed pairs and their respective 
#   separation distances
dSpCopula <- function (u, copula) {
  if (!is.list(u) || !length(u)>=2) stop("Point pairs need to be provided with their separating distance as a list.")
  
  pairs <- u[[1]]
  n <- nrow(pairs)
  
  if(length(u)==3) {
    block <- u[[3]]
    if (n%%block != 0) stop("The block size is not a multiple of the data length:",n)
  } else block <- 1
  
  h <- u[[2]]
  if(length(h)>1 && length(h)!=nrow(u[[1]])) {
    stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
  }
  
  if(is.null(copula@calibMoa(normalCopula(0),0))){
    res <- spConCop(dCopula, copula, pairs, rep(h, length.out=nrow(pairs)))
  }
  else {
    if(length(h)>1) {
      if (block == 1){
        ordering <- order(h)
        
        # ascending sorted pairs allow for easy evaluation
        pairs <- pairs[ordering,,drop=FALSE] 
        h <- h[ordering]
        
        res <- spDepFunCop(dCopula, copula, pairs, h)
        
        # reordering the values
        res <- res[order(ordering)]
      } else {
        res <- NULL
        for(i in 1:(n%/%block)) {
          res <- c(res, spDepFunCopSnglDist(dCopula, copula, 
                                            pairs[((i-1)*block+1):(i*block),], 
                                            h[i*block]))
        }
      }
    } else {
      res <- spDepFunCopSnglDist(dCopula, copula, pairs, h)
    }
  }
  
  return(res)
}

setMethod("dCopula", signature("list","spCopula"), dSpCopula)


## partial derivatives ##

## dduSpCopula
###############

dduSpCopula <- function (u, copula) {
  if (!is.list(u) || !length(u)>=2) stop("Point pairs need to be provided with their separating distance as a list.")
  
  pairs <- u[[1]]
  n <- nrow(pairs)
  
  if(length(u)==3) {
    block <- u[[3]]
    if (n%%block != 0) stop("The block size is not a multiple of the data length:",n)
  } else block <- 1
  
  h <- u[[2]]
  if(length(h)>1 && length(h)!=nrow(u[[1]])) {
    stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
  }

  if(is.null(copula@calibMoa(normalCopula(0),0))) res <- spConCop(dduCopula, copula, pairs, 
                                                 rep(h,length.out=nrow(pairs)))
  else {
    if(length(h)>1) {
      if (block == 1){
        ordering <- order(h)
        
        # ascending sorted pairs allow for easy evaluation
        pairs <- pairs[ordering,,drop=FALSE] 
        h <- h[ordering]
        
        res <- spDepFunCop(dduCopula, copula, pairs, h)
        
        # reordering the values
        res <- res[order(ordering)]
      } else {
        res <- NULL
        for(i in 1:(n%/%block)) {
          res <- c(res, spDepFunCopSnglDist(dduCopula, copula, 
                                            pairs[((i-1)*block+1):(i*block),],
                                            h[i*block]))
        }
      }
    } else {
      res <- spDepFunCopSnglDist(dduCopula, copula, pairs, h)
    }
  }
  
  return(res)
}

setMethod("dduCopula", signature("numeric","spCopula"), dduSpCopula)

## ddvSpCopula
###############

ddvSpCopula <- function (u, copula) {
  if (!is.list(u) || !length(u)>=2) stop("Point pairs need to be provided with their separating distance as a list.")
  
  pairs <- u[[1]]
  n <- nrow(pairs)
  
  if(length(u)==3) {
    block <- u[[3]]
    if (n%%block != 0) stop("The block size is not a multiple of the data length:",n)
  } else block <- 1
  
  h <- u[[2]]
  if(length(h)>1 && length(h)!=nrow(u[[1]])) {
    stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
  }
  
  
  if(is.null(copula@calibMoa(normalCopula(0),0))) res <- spConCop(ddvCopula, copula, pairs, 
                                                 rep(h,length.out=nrow(pairs)))
  else {
    if(length(h)>1) {
      if (block == 1){
        ordering <- order(h)
        
        # ascending sorted pairs allow for easy evaluation
        pairs <- pairs[ordering,,drop=FALSE] 
        h <- h[ordering]
        
        res <- spDepFunCop(ddvCopula, copula, pairs, h)
        
        # reordering the values
        res <- res[order(ordering)]
      } else {
        res <- NULL
        for(i in 1:(n%/%block)) {
          res <- c(res, spDepFunCopSnglDist(ddvCopula, copula, 
                                            pairs[((i-1)*block+1):(i*block),],
                                            h[i*block]))
        }
      }
    } else {
      res <- spDepFunCopSnglDist(ddvCopula, copula, pairs, h)
    }
  }
  
  return(res)
}

setMethod("ddvCopula", signature("numeric","spCopula"), ddvSpCopula)


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
fitCorFun <- function(bins, type="poly", degree=3, cutoff=NA, bounds=c(0,1), 
                      method="kendall") {
  bins <- as.data.frame(bins[1:2])
  if(!is.na(cutoff)) bins <- bins[which(bins[[1]] <= cutoff),]
  
  fitCor <- lm(lagCor ~ poly(meanDists, degree), data = bins)
  
  print(fitCor)
  cat("Sum of squared residuals:",sum(fitCor$residuals^2),"\n")
  
  function(x) {
    if (is.null(x)) return(method)
    return(pmin(bounds[2], pmax(bounds[1], 
                                eval(predict(fitCor, 
                                             data.frame(meanDists=x))))))
  }
}


# towards b)
# bins     -> typically output from calcBins
# calcTau  -> a function on distance providing Kendall's tau estimates, 
# families -> a vector of dummy copula objects of each family to be considered
#             DEFAULT: c(normal, t_df=4, clayton, frank, gumbel
loglikByCopulasLags <- function(bins, calcTau, 
                                families=c(normalCopula(0), 
                                           tCopula(0,dispstr="un"),
                                           claytonCopula(0), frankCopula(1), 
                                           gumbelCopula(1))) {
  loglik <- NULL
  for (cop in families) {
    print(cop)
    tmploglik <- NULL
    for(i in 1:length(bins$meanDists)) {
      cop@parameters[1] <- iTau(cop, tau=calcTau(bins$meanDists[i]))
      tmploglik <- c(tmploglik, sum(log(dCopula(bins$lagData[[i]],cop))))
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
    distances <- c(distances, 
                   rev(distances)[1]+diff(bins$meanDists[nfits+c(-1,0)]))
  }

  if(missing(calcCor)) return(spCopula(components=cops, distances=distances, 
                                       unit="m"))
  else return(spCopula(components=cops, distances=distances, unit="m", 
                       spDepFun=calcCor))
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
fitSpCopula <- function(bins, cutoff=NA, 
                        families=c(normalCopula(0), 
                                   tCopula(0,dispstr="un"), claytonCopula(0), 
                                   frankCopula(1), gumbelCopula(1)), ...) {
  calcTau <- fitCorFun(bins, cutoff=cutoff, ...)
  loglik <- loglikByCopulasLags(bins, calcTau=calcTau, families=families)
  
  bestFit <- apply(apply(loglik, 1, rank),2, 
                   function(x) which(x==length(families)))
  
  return(composeSpCop(bestFit, families, bins, calcTau))
}
