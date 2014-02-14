######################################
## Spatio-Temporal Bivariate Copula ##
######################################

## constructor ##
#################

stCopula <- function(components, t.lags, distances=NA, stDepFun, unit="m", t.res="day") {
  if(all(sapply(components, function(x) class(x)=="spCopula"))) {
    if(length(unique(sapply(components, function(x) x@unit))) >1 )
      stop("All spatial copulas need to have the same distance unit.")
    stopifnot(length(t.lags) == length(components))
    spCopList <- components
  } else {
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

# dropping a sptio-temporal tree
dropStTree <- function (stNeigh, dataLocs, stCop) {
  stopifnot(class(stNeigh) == "stNeighbourhood")
  
  u0 <- as.matrix(stNeigh@data)
  h0 <- stNeigh@distances
  u1 <- matrix(NA, nrow(u0), ncol(u0)-1-length(stNeigh@coVar))
  h1 <- array(dim = c(nrow(u0), ncol(h0)-1, 2))
  
  pb <- txtProgressBar(0,dim(h0)[2],style=3)
  for (i in 1:dim(h0)[2]) {
    u1[,i] <- dduCopula(u0[, c(1, i + 1)], stCop, h = h0[, i, ])
    if (i < ncol(h0)) {
      h1[,i,1] <- apply(stNeigh@index[, c(1, i + 1), 1], 1, 
                        function(x) spDists(dataLocs@sp[x, ])[1, 2])
      h1[,i,2] <- apply(stNeigh@index[, c(1, i + 1), 2], 1, 
                        function(x) diff(x))
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
#   # add covariate to the conditioned neighbourhood?
#   if (length(stNeigh@coVar) > 0)
#     u1[,ncol(u0)-(1:length(stNeigh@coVar))] <- u0[,ncol(u0) + 1 - (1:length(stNeigh@coVar))]
  
  varSplit <- strsplit(stNeigh@var, "|", fixed = TRUE)[[1]]
  cond <- suppressWarnings(as.numeric(varSplit[length(varSplit)]))
  
  if (is.na(cond)) {
#     coVar <- paste(stNeigh@coVar, "|0", sep = "")
    cond <- paste(stNeigh@var, "|0", sep = "")
  }
  else {
#     coVar <- stNeigh@coVar
    cond <- paste(stNeigh@var, cond + 1, sep = "")
  }
  
  return(stNeighbourhood(data = u1, distances = h1, index = stNeigh@index[, -1, ],
                         var = cond, prediction = stNeigh@prediction))
}