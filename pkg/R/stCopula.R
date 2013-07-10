######################################
## Spatio-Temporal Bivariate Copula ##
######################################

## constructor ##
#################

stCopula <- function(components, distances, t.lags, stDepFun, unit="m", t.res="day") {
  spCopList <- list()
  
  if(!missing(stDepFun)) {
    getSpCop <- function(comp,dist,time) spCopula(comp, dist,
                                                  spDepFun=function(h) stDepFun(h,time), unit)
    for(i in 1:length(t.lags)){
      spCopList <- append(spCopList, getSpCop(components[[i]], distances[[i]], i))
    }
  } else {
    for(i in 1:length(t.lags)){
      spCopList <- append(spCopList, spCopula(components[[i]], distances[[i]], unit=unit))
    }
  }
  
  param       <- unlist(lapply(spCopList, function(x) x@parameters))
  param.names <- unlist(lapply(spCopList, function(x) x@param.names))
  param.low   <- unlist(lapply(spCopList, function(x) x@param.lowbnd))
  param.up    <- unlist(lapply(spCopList, function(x) x@param.upbnd))
  
  new("stCopula", dimension=as.integer(2), parameters=param, param.names=param.names,
      param.lowbnd=param.low, param.upbnd=param.up,
      fullname="Spatio-Temporal Copula: distance and time dependent convex combination of bivariate copulas",
      spCopList=spCopList, t.lags=t.lags, t.res=t.res)
}

## show method ##
#################

showStCopula <- function(object) {
  cat(object@fullname, "\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Copulas:\n")
  for (i in 1:length(object@spCopList)) {
    cmpCop <- object@spCopList[[i]]
    cat("  ", cmpCop@fullname, "at", object@t.lags[i], 
      paste("[",object@t.res,"]",sep=""), "\n")
    show(cmpCop)
  }
}

setMethod("show", signature("stCopula"), showStCopula)

## spatial copula cdf ##
########################

pStCopula <- function (u, copula, h) {
  stopifnot(ncol(h)==2)
  stopifnot(nrow(h)==1 || nrow(h)==nrow(u))
  
  n <- nrow(u)
  tDist <- unique(h[,2])
  
  if(any(is.na(match(tDist,copula@t.lags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  
  if (length(tDist)==1) {
    res <- pSpCopula(u, copula@spCopList[[match(tDist, copula@t.lags)]], h[,1])
  } else {
    res <- numeric(n)
    for(t in tDist) {
      tmpInd <- h[,2]==t
      tmpCop <- copula@spCopList[[match(t, copula@t.lags)]]
      res[tmpInd] <- pSpCopula(u[tmpInd,,drop=F], tmpCop, h[tmpInd,1])
    }
  }
  res
}

setMethod(pCopula, signature("numeric","stCopula"), 
          function(u, copula, log, ...) pStCopula(matrix(u,ncol=2), copula, ...))
setMethod(pCopula, signature("matrix","stCopula"), pStCopula)

## spatial Copula density ##
############################

dStCopula <- function (u, copula, log, h) {
  stopifnot(ncol(h)==2)
  stopifnot(nrow(h)==1 || nrow(h)==nrow(u))
  
  n <- nrow(u)
  tDist <- unique(h[,2])
  
  if(any(is.na(match(tDist,copula@t.lags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  
  if (length(tDist)==1) {
    res <- dSpCopula(u, copula@spCopList[[match(tDist, copula@t.lags)]], log, h[,1])
  } else {
    res <- numeric(n)
    for(t in tDist) {
      tmpInd <- h[,2]==t
      tmpCop <- copula@spCopList[[match(t, copula@t.lags)]]
      res[tmpInd] <- dSpCopula(u[tmpInd,,drop=F], tmpCop, log, h[tmpInd,1])
    }
  }
  res
}

setMethod(dCopula, signature("numeric","stCopula"), 
          function(u, copula, log, ...) dStCopula(matrix(u,ncol=2), copula, log=log, ...))
setMethod(dCopula, signature("matrix","stCopula"), dStCopula)


## partial derivatives ##

## dduSpCopula ##
#################

dduStCopula <- function (u, copula, h) {
  stopifnot(ncol(h)==2)
  stopifnot(nrow(h)==1 || nrow(h)==nrow(u))
  
  n <- nrow(u)
  tDist <- unique(h[,2])
  
  if(any(is.na(match(tDist,copula@t.lags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  
  if (length(tDist)==1) {
    res <- dduSpCopula(u, copula@spCopList[[match(tDist, copula@t.lags)]], h[,1])
  } else {
    res <- numeric(n)
    for(t in tDist) {
      tmpInd <- h[,2]==t
      tmpCop <- copula@spCopList[[match(t, copula@t.lags)]]
      res[tmpInd] <- dduSpCopula(u[tmpInd,,drop=F], tmpCop, h[tmpInd,1])
    }
  }
  res
}

setMethod("dduCopula", signature("numeric","stCopula"), 
          function(u, copula, ...) dduStCopula(matrix(u,ncol=2), copula, ...))
setMethod("dduCopula", signature("matrix","stCopula"), dduStCopula)


## ddvSpCopula ##
#################

ddvStCopula <- function (u, copula, h) {
  stopifnot(ncol(h)==2)
  stopifnot(nrow(h)==1 || nrow(h)==nrow(u))
  
  n <- nrow(u)
  tDist <- unique(h[,2])
  
  if(any(is.na(match(tDist,copula@t.lags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  
  if (length(tDist)==1) {
    res <- ddvSpCopula(u, copula@spCopList[[match(tDist,copula@t.lags)]], h[,1])
  } else {
    res <- numeric(n)
    for(t in tDist) {
      tmpInd <- h[,2]==t
      tmpCop <- copula@spCopList[[match(t, copula@t.lags)]]
      res[tmpInd] <- ddvSpCopula(u[tmpInd,,drop=F], tmpCop, h[tmpInd,1])
    }
  }
  res
}

setMethod("ddvCopula", signature("numeric","stCopula"), 
          function(u, copula, ...) ddvStCopula(matrix(u,ncol=2), copula, ...))
setMethod("ddvCopula", signature("matrix","stCopula"), ddvStCopula)