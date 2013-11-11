#########################################
## methods for the spatial vine copula ##
#########################################

## constructor ##
#################

stVineCopula <- function(stCop, topCop) {
  stopifnot(is(stCop,"stCopula"))
  
  new("stVineCopula", dimension = as.integer(topCop@dimension+1),
      parameters=numeric(), param.names = character(), param.lowbnd = numeric(), 
      param.upbnd = numeric(), 
      fullname = paste("Spatio-temporal vine copula family with 1 spatio-temporal tree."),
      stCop=stCop, topCop=topCop)
}

## show ##
##########

showStVineCopula <- function(object) {
  dim <- object@dimension
  cat(object@fullname, "\n")
  cat("Dimension: ", dim, "\n")
}

setMethod("show", signature("stVineCopula"), showStVineCopula)

## density ##
#############

dstVine <- function(u, stCop, topCop, log, h) {
  l0 <- rep(0,nrow(u)) # level 0 spatio-temporal density
  dimDists <- dim(h)
  stopifnot(length(dimDists)==3)
  
  u0 <- u # previous level's conditional data
  u1 <- NULL # current level of conditional data
  for(i in 1:dimDists[2]) { # i <- 1
    l0 <- l0 + dCopula(u0[,c(1,i+1)], stCop, h=matrix(h[,i,],ncol=2), log=T)
    u1 <- cbind(u1, dduCopula(u0[,c(1,i+1)], stCop, h=matrix(h[,i,],ncol=2)))
  }
  u0 <- u1
  
  if(!is.null(topCop))
    l1 <- dCopula(u0, topCop, log=T)
  else 
    l1 <- 0
  
  if(log)
    return(l0+l1)
  else(exp(l0+l1))
}

setMethod("dCopula",signature=signature("matrix","stVineCopula"),
          function(u, copula, log, ...) {
            if("topCop" %in% slotNames(copula))
              dstVine(u, copula@stCop, copula@topCop, log=log, ...)
            else
              dstVine(u, copula@stCop, NULL, log=log, ...)
          })

setMethod("dCopula",signature=signature("numeric","stVineCopula"),
          function(u, copula, log, ...) {
            if("topCop" %in% slotNames(copula))
              dstVine(matrix(u,ncol=copula@dimension), copula@stCop, copula@topCop, log=log, ...)
            else
              dstVine(matrix(u,ncol=copula@dimension), copula@stCop, NULL, log=log, ...)
          })

setMethod("dCopula",signature=signature("data.frame","stVineCopula"),
          function(u, copula, log, ...) {
            if("topCop" %in% slotNames(copula))
              dstVine(as.matrix(u), copula@stCop, copula@topCop, log=log, ...)
            else
              dstVine(as.matrix(u), copula@stCop, NULL, log=log, ...)
          })

# fiiting the spatial vine for a given list of spatial copulas
fitStVine <- function(copula, data, method, estimate.variance=F) {
  stopifnot(class(data)=="stNeighbourhood")
  stopifnot(copula@dimension == ncol(data@data))
  
  u0 <- as.matrix(data@data) # previous level's (conditional) data
  h0 <- data@distances # previous level's distances
  l0 <- rep(0,nrow(u0)) # spatial density
  u1 <- NULL # current level of conditional data
  cat("[Margin ")
  for(i in 1:dim(h0)[2]) { # i <- 1
    l0 <- l0 + dCopula(u0[,c(1,i+1)], copula@stCop, h=h0[,i,], log=T)
    cat(i,", ", sep="")
    u1 <- cbind(u1, dduCopula(u0[,c(1,i+1)], copula@stCop, h=h0[,i,]))
  }
  u0 <- u1
  cat("]\n")
  
  cat("[Estimating a",ncol(u0),"dimensional copula at the top.]\n")
  vineCopFit <- fitCopula(copula@topCop, u0, method, estimate.variance) 
    
  stVineCop <- stVineCopula(copula@stCop, vineCopFit@copula)
  loglik <- vineCopFit@loglik
  
  return(new("fitCopula", estimate = stVineCop@parameters, var.est = matrix(NA), 
             method = method, 
             loglik = sum(l0)+loglik,
             fitting.stats=list(convergence = as.integer(NA)),
             nsample = nrow(data@data), copula=stVineCop))
}

setMethod("fitCopula",signature=signature("stVineCopula"),fitStVine)

# conditional spatial vine
condStVine <- function (condVar, dists, stVine, n = 1000) {
  stopifnot(is.array(dists))
  
  # add some points in the tails
  rat <- 50:1%x%c(1e-6,1e-5,1e-4,1e-3)
  xVals <- unique(sort(c(rat, 1 - rat, 1:(n - 1)/n)))
  nx <- length(xVals)
  
  repCondVar <- matrix(condVar, ncol = length(condVar), nrow = nx, byrow = T)
  density <- dCopula(cbind(xVals, repCondVar), stVine, h = dists)
  
  # the 1-e6 corners linearily to [0,1], but ensure non-negative
  density <- c(max(0,2*density[1]-density[2]),
               density, max(0,2*density[nx]-density[nx-1]))
  linAppr <- approxfun(c(0, xVals, 1), density)
  
  # sum up the denstiy to rescale
  int <- sum(diff(c(0,xVals,1))*(0.5*diff(density)+density[-(nx+2)]))
  condVineFun <- function(u) linAppr(u)/int
  attr(condVineFun,"xVals") <- c(0,xVals,1)
  return(condVineFun)
}

## interpolation ##
###################

stCopPredict.expectation <- function(predNeigh, stVine, margin, ..., stop.on.error=F) {
  stopifnot(is.function(margin$q))
  dists <- predNeigh@distances
  
  predMean <- NULL
  for(i in 1:nrow(predNeigh@data)) { # i <-1
    cat("[Predicting location ",i,".]\n", sep="")
    condSecVine <- condStVine(as.numeric(predNeigh@data[i,]), dists[i,], stVine)
    
    condExp <-  function(x) {
      margin$q(x)*condSecVine(x)
    }
    
    ePred <- integrate(condExp,0,1,subdivisions=10000L,stop.on.error=stop.on.error, ...)
    if(ePred$abs.error > 0.01)
            warning("Numerical integration in predExpectation performed at a level of absolute error of only ",
                    ePred$abs.error, " for location ",i,".")
    predMean <- c(predMean, ePred$value)
  }
  if ("data" %in% slotNames(predNeigh@locations)) {
    res <- predNeigh@locations
    res@data[["expect"]] <- predMean
    return(res)
  } else {
    predMean <- data.frame(predMean)
    colnames(predMean) <- "expect"
    return(addAttrToGeom(predNeigh@locations, predMean, match.ID=FALSE))
  }
}

stCopPredict.quantile <- function(predNeigh, stVine, margin, p=0.5) {
  stopifnot(is.function(margin$q))
  dists <- predNeigh@distances
  
  predQuantile <- NULL
  for(i in 1:nrow(predNeigh@data)) { # i <-1
    cat("[Predicting location ",i,".]\n", sep="")
    condSecVine <- condStVine(as.numeric(predNeigh@data[i,]), dists[i,,,drop=F], stVine)
    
    xVals <- attr(condSecVine,"xVals")
    density <- condSecVine(xVals)
    nx <- length(xVals)
    int <- cumsum(c(0,diff(xVals)*(0.5*diff(density)+density[-nx])))
    lower <- max(which(int <= p))
    m <- (density[lower+1]-density[lower])/(xVals[lower+1]-xVals[lower])
    b <- density[lower]
    xRes <- -b/m+sign(m)*sqrt(b^2/m^2+2*(p-int[lower])/m)
    
    predQuantile <- c(predQuantile, margin$q(xVals[lower]+xRes))
  }
  
  if ("data" %in% slotNames(predNeigh@locations)) {
    res <- predNeigh@locations
    res@data[[paste("quantile.",p,sep="")]] <- predQuantile
    return(res)
  } else {
    predQuantile <- data.frame(predQuantile)
    colnames(predQuantile) <- paste("quantile.",p,sep="")
    return(addAttrToGeom(predNeigh@locations, predQuantile, match.ID=FALSE))
  }
}

stCopPredict <- function(predNeigh, stVine, margin, method="quantile", p=0.5, ...) {
  switch(method,
         quantile=stCopPredict.quantile(predNeigh, stVine, margin, p),
         expectation=stCopPredict.expectation(predNeigh, stVine, margin, ...))
}

dropStTree <- function(neigh, stCop) {
    stopifnot(class(neigh)=="stNeighbourhood")
    
    u0 <- as.matrix(neigh@data) # previous level's (conditional) data
    h0 <- neigh@distances # previous level's distances
    u1 <- NULL # current level of conditional data
    h1s <- NULL # upcoming distances
    h1t <- NULL # upcoming distances
    cat("[Margin ")
    for(i in 1:dim(h0)[2]) { # i <- 1
      cat(i,", ")
      u1 <- cbind(u1, dduCopula(u0[,c(1,i+1)], stCop, h=h0[,i,]))
      if (i < ncol(neigh@distances)) {
        h1s <- cbind(h1s, apply(neigh@index[, c(1, i + 1),1], 1, 
                                function(x) spDists(neigh@locations@sp[x, ])[1, 2]))
        h1t <- cbind(h1t, apply(neigh@index[, c(1, i + 1),2], 1, 
                                function(x) diff(x)))
      }
    }
    h1 <- array(dim=c(dim(h1s),2))
    h1[,,1] <- h1s
    h1[,,2] <- h1t
    
    varSplit <- strsplit(neigh@var, "|", fixed = TRUE)[[1]]
    cond <- suppressWarnings(as.numeric(varSplit[length(varSplit)]))
    if (is.na(cond)) 
      cond <- paste(neigh@var, "|0", sep = "")
    else cond <- paste(neigh@var, cond + 1, sep = "")
    return(stNeighbourhood(data=u1, distances=h1, STxDF=neigh@locations, 
                           ST=neigh@dataLocs, index=neigh@index[, -1,], 
                           prediction=neigh@prediction, var=cond))
}