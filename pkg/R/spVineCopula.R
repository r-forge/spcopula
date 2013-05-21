#########################################
## methods for the spatial vine copula ##
#########################################

# constructor
spVineCopula <- function(spCop, topCop=NULL) {
  if(!is.list(spCop))
    spCop <- list(spCop)
  
  if(is.null(topCop))
    new("pureSpVineCopula", dimension = as.integer(length(spCop)+1),
        parameters=numeric(), param.names = character(), param.lowbnd = numeric(), 
        param.upbnd = numeric(), 
        fullname = paste("Spatial vine copula family with only spatial tree(s)."),
        spCop=spCop)
  else
    new("mixedSpVineCopula", dimension = as.integer(topCop@dimension+length(spCop)),
      parameters=numeric(), param.names = character(), param.lowbnd = numeric(), 
      param.upbnd = numeric(), 
      fullname = paste("Spatial vine copula family with",length(spCop),
                       "spatial tree(s)."),
      spCop=spCop, topCop=topCop)
}

# show
showSpVineCopula <- function(object) {
  dim <- object@dimension
  cat(object@fullname, "\n")
  cat("Dimension: ", dim, "\n")
}

setMethod("show", signature("spVineCopula"), showSpVineCopula)

# density
dspVine <- function(u, spCop, topCop, log, h) {
  l0 <- rep(0,nrow(u)) # level 0 spatial density
  stopifnot(is.list(h))
  stopifnot(length(spCop)==length(h))
  
  u0 <- u # previous level's conditional data
  for(spTree in 1:length(spCop)) {
#     cat("[Evaluating spatial tree",spTree,"]\n")
    tmpH <- h[[spTree]]
    if(!is.matrix(tmpH)) 
      tmpH <- matrix(tmpH, ncol=length(tmpH))

    u1 <- NULL # current level of conditional data
    for(i in 1:ncol(tmpH)) { # i <- 1
      l0 <- l0 + dCopula(u0[,c(1,i+1)], spCop[[spTree]], h=tmpH[,i], log=T)
      u1 <- cbind(u1, dduCopula(u0[,c(1,i+1)], spCop[[spTree]], h=tmpH[,i]))
    }
    u0 <- u1
  }
  
  if(!is.null(topCop))
    l1 <- dCopula(u0, topCop, log=T)
  else 
    l1 <- 0
  
  if(log)
    return(l0+l1)
  else(exp(l0+l1))
}

setMethod("dCopula",signature=signature("matrix","spVineCopula"),
          function(u, copula, log, ...) {
            if("topCop" %in% slotNames(copula))
              dspVine(u, copula@spCop, copula@topCop, log=log, ...)
            else
              dspVine(u, copula@spCop, NULL, log=log, ...)
          })

setMethod("dCopula",signature=signature("numeric","spVineCopula"),
          function(u, copula, log, ...) {
            if("topCop" %in% slotNames(copula))
              dspVine(matrix(u,ncol=copula@dimension), copula@spCop, copula@topCop, log=log, ...)
            else
              dspVine(matrix(u,ncol=copula@dimension), copula@spCop, NULL, log=log, ...)
          })

setMethod("dCopula",signature=signature("data.frame","spVineCopula"),
          function(u, copula, log, ...) {
            if("topCop" %in% slotNames(copula))
              dspVine(as.matrix(u), copula@spCop, copula@topCop, log=log, ...)
            else
              dspVine(as.matrix(u), copula@spCop, NULL, log=log, ...)
          })

# fiiting the spatial vine for a given list of spatial copulas
fitSpVine <- function(copula, data, method, estimate.variance=F) {
  stopifnot(class(data)=="neighbourhood")
  stopifnot(copula@dimension == ncol(data@data))
  
  u0 <- as.matrix(data@data) # previous level's (contitional) data
  h0 <- data@distances # previous level's distances
  l0 <- rep(0,nrow(u0)) # spatial density
  for(spTree in 1:length(copula@spCop)) {
    cat("[Dropping ", spTree, ". spatial tree.]\n",sep="")
    u1 <- NULL # current level of conditional data
    h1 <- NULL # current level's distances
    for(i in 1:ncol(h0)) { # i <- 1
      l0 <- l0 + dCopula(u0[,c(1,i+1)], copula@spCop[[spTree]], h=h0[,i], log=T)
      u1 <- cbind(u1, dduCopula(u0[,c(1,i+1)], copula@spCop[[spTree]], h=h0[,i]))
      if (i < ncol(h0)) {
        h1 <- cbind(h1,apply(data@index[,c(spTree,spTree+i)],1, 
                             function(x) spDists(data@locations[x,])[1,2]))
      }
    }
    u0 <- u1
    h0 <- h1
  }
  
  if (ncol(u0)==1) {
    cat("[No copula to be estimated at the top.]\n")
    
    spVineCop <- spVineCopula(copula@spCop)
    loglik <- 0
  }
  if (ncol(u0)==2) {
    cat("[Estimating a single bivariate copula at the top.]\n")
    bivCop <- BiCopSelect(u0[,1],u0[,2])
    topCop <- copulaFromFamilyIndex(bivCop$family, bivCop$par, bivCop$par2)
    
    spVineCop <- spVineCopula(copula@spCop, topCop)
    loglik <- dCopula(topCop,u0,log=TRUE)
  } else {
    cat("[Estimating a",ncol(u0),"dimensional copula at the top.]\n")
    vineCopFit <- fitCopula(copula@topCop, u0, method, estimate.variance) 
    
    spVineCop <- spVineCopula(copula@spCop, vineCopFit@copula)
    loglik <- vineCopFit@loglik
  }
  
  return(new("fitCopula", estimate = spVineCop@parameters, var.est = matrix(NA), 
             method = method, 
             loglik = sum(l0)+loglik,
             fitting.stats=list(convergence = as.integer(NA)),
             nsample = nrow(data@data), copula=spVineCop))
}

setMethod("fitCopula",signature=signature("spVineCopula"),fitSpVine)

# deriving all spatial tree distances
calcSpTreeDists <- function(neigh, n.trees) {
  if(!neigh@prediction)
    data <- neigh@locations
  else
    data <- neigh@dataLocs
  
  condDists <- list(n.trees)
  condDists[[1]] <- neigh@distances
  if(n.trees==1)
    return(condDists)
  for (spTree in 1:(n.trees-1)) {
    h1 <- NULL
    for(i in 1:(ncol(neigh@distances)-spTree)) {
      h1 <- cbind(h1,apply(neigh@index[,c(spTree,i+spTree),drop=F],1, 
                           function(x) spDists(data[x,])[1,2]))
      dimnames(h1) <- NULL
    }
    condDists[[spTree+1]] <- h1
  }
  return(condDists)
}

# conditional spatial vine
condSpVine <- function (condVar, dists, spVine, n = 1000) {
  stopifnot(is.list(dists))
  stopifnot(length(spVine@spCop)==length(dists))
  
  # add some points in the tails
  rat <- 50:1%x%c(1e-6,1e-5,1e-4,1e-3)
  xVals <- unique(sort(c(rat, 1 - rat, 1:(n - 1)/n)))
  nx <- length(xVals)
  
  repCondVar <- matrix(condVar, ncol = length(condVar), nrow = nx, byrow = T)
  density <- dCopula(cbind(xVals, repCondVar), spVine, h = dists)
  
  # the 1-e6 corners linearily to [0,1], but keep non-negative
  density <- c(max(0,2*density[1]-density[2]),
               density, max(0,2*density[nx]-density[nx-1]))
  linAppr <- approxfun(c(0, xVals, 1), density)
  
  # sum up the denstiy to rescale
  int <- sum(diff(c(0,xVals,1))*(0.5*diff(density)+density[-(nx+2)]))
  condVineFun <- function(u) linAppr(u)/int
  attr(condVineFun,"xVals") <- c(0,xVals,1)
  return(condVineFun)
}

# interpolation

spCopPredict.expectation <- function(predNeigh, spVine, margin, ..., stop.on.error=F) {
#   stopifnot(!is.null(range))
#   stopifnot(is.function(margin$d))
#   stopifnot(is.function(margin$p))
  stopifnot(is.function(margin$q))
  
  dists <- calcSpTreeDists(predNeigh,length(spVine@spCop))
  
  predMean <- NULL
  for(i in 1:nrow(predNeigh@data)) { # i <-1
    cat("[Predicting location ",i,".]\n", sep="")
    condSecVine <- condSpVine(as.numeric(predNeigh@data[i,]), 
                              lapply(dists,function(x) x[i,]), spVine)
    
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

spCopPredict.quantile <- function(predNeigh, spVine, margin, p=0.5) {
  stopifnot(is.function(margin$q))
  dists <- calcSpTreeDists(predNeigh,length(spVine@spCop))
  
  predQuantile <- NULL
  for(i in 1:nrow(predNeigh@data)) { # i <-1
    cat("[Predicting location ",i,".]\n", sep="")
    condSecVine <- condSpVine(as.numeric(predNeigh@data[i,]), 
                              lapply(dists,function(x) x[i,]), spVine)
    
    xVals <- attr(condSecVine,"xVals")
    density <- condSecVine(xVals)
    nx <- length(xVals)
    int <- cumsum(c(0,diff(xVals)*(0.5*diff(density)+density[-nx])))
    lower <- max(which(int <= p))
    m <- (density[lower+1]-density[lower])/(xVals[lower+1]-xVals[lower])
    b <- density[lower]
    xRes <- -b/m+sign(m)*sqrt(b^2/m^2+2*(p-int[lower])/m)
    
#     pPred <- optimise(function(x) abs(integrate(condSecVine, 0, x, 
#                                                 subdivisions=10000L, 
#                                                 abs.tol=1e-6)$value-p), c(0,1))
#     if(pPred$objective > 1e-4)
#       warning("Numerical evaluation in predQuantile achieved an objective of only ",
#               pPred$objective, " where 0 has been sought for location ",i,".")
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

spCopPredict <- function(predNeigh, spVine, margin, method="quantile", p=0.5, ...) {
  switch(method,
         quantile=spCopPredict.quantile(predNeigh, spVine, margin, p),
         expectation=spCopPredict.expectation(predNeigh, spVine, margin, ...))
}

# draw from a spatial vine
# Algorithm 1 from Aas et al. (2006): Pair-copula constructions of multiple dependence

r.spVineCop <- function(n, spVine, h) {
  spVineDim <- spVine@dimension
  
  sims <- NULL
  for(runs in 1:n) {
    init <- runif(spVineDim)
    res <- init[1]
    v <- matrix(NA,spVineDim,spVineDim)
    v[1,1] <- init[1]
    for (i in 2:spVineDim) { # i <- 2
      v[i,1] <- init[i]
      for (k in (i-1):1) { # k <- i-1
        v[i,1] <- uniroot(function(u) {
                            v[i,1] - ddvCopula(cbind(u,v[k,k]), spVine@spCop[[k]],
                                               h=h[[k]][i-k])
                          }, c(0,1))$root
      }  
      res <- c(res,v[i,1])
      if(i==spVineDim)
        break()
      for(j in 1:(i-1)) {
        v[i,j+1] <- ddvCopula(cbind(v[i,j],v[j,j]),spVine@spCop[[k]], h=h[[j]][i-j])
      }
    }
    sims <- rbind(sims,res)
  }
  
  rownames(sims) <- NULL
  sims
}

setMethod("rCopula", signature("numeric","spVineCopula"), 
          function(n, copula, ...) r.spVineCop(n, copula, ...))