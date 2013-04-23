#########################################
## methods for the spatial vine copula ##
#########################################

# constructor
spVineCopula <- function(spCop, vineCop=vineCopula()) {
  new("spVineCopula", dimension = as.integer(vineCop@dimension+1), parameters=numeric(),
      param.names = character(), param.lowbnd = numeric(), 
      param.upbnd = numeric(), fullname = "Spatial vine copula family.",
      spCop=spCop, vineCop=vineCop)
}

# show
showSpVineCopula <- function(object) {
  dim <- object@dimension
  cat(object@fullname, "\n")
  cat("Dimension: ", dim, "\n")
}

setMethod("show", signature("spVineCopula"), showSpVineCopula)

# density
dspVine <- function(u, spCop, vine, log, h) {
  l0 <- rep(0,nrow(u)) # level 0 (spatial) density
  u0 <- NULL # level 0 conditional data
  
  if(!is.matrix(h)) 
    h <- matrix(h, ncol=length(h))
  
  for(i in 1:(ncol(u)-1)) { # i <- 1
    l0 <- l0 + dCopula(u[,c(1,i+1)], spCop, h=h[,i], log=T)
    u0 <- cbind(u0, dduCopula(u[,c(1,i+1)], spCop, h=h[,i]))
  }

  l1 <- dCopula(u0, vine, log=T)
  if(log)
    return(l0+l1)
  else(exp(l0+l1))
}

setMethod("dCopula",signature=signature("matrix","spVineCopula"),
          function(u, copula, log, ...) {
            dspVine(u, copula@spCop, copula@vineCop, log=log, ...)
          })

setMethod("dCopula",signature=signature("numeric","spVineCopula"),
          function(u, copula, log, ...) {
            dspVine(matrix(u,ncol=copula@dimension), copula@spCop, copula@vineCop, log=log, ...)
          })

setMethod("dCopula",signature=signature("data.frame","spVineCopula"),
          function(u, copula, log, ...) {
            dspVine(as.matrix(u), copula@spCop, copula@vineCop, log=log, ...)
          })

# fiiting the spatial vine for a given spatial copula

fitSpVine <- function(copula, data, method) {
  stopifnot(class(data)=="neighbourhood")
  stopifnot(copula@dimension == ncol(data@data))
  
  secLevel <- NULL
  for (i in 1:(copula@dimension-1)) { # i <- 1
    secLevel <- cbind(secLevel, 
                      dduCopula(u=as.matrix(data@data[,c(1,i+1)]), 
                                copula=copula@spCop, h=data@distances[,i]))
  }
  
  cat(summary(as.data.frame(secLevel)))
    
  vineCopFit <- fitCopula(copula@vineCop, secLevel, method) 
  
  spVineCop <- spVineCopula(copula@spCop, vineCopFit@copula)
    
  return(new("fitCopula", estimate = spVineCop@parameters, var.est = matrix(NA), 
             method = method, 
             loglik = sum(dCopula(data@data, spVineCop, h=data@distances, log=T)),
             fitting.stats=list(convergence = as.integer(NA)),
             nsample = nrow(data@data), copula=spVineCop))
}

setMethod("fitCopula",signature=signature("spVineCopula"),fitSpVine)

# conditional spatial vine
condSpVine <- function (condVar, dists, spVine, n = 1000) {
  # add some points in the tails
  rat <- 29:1%x%c(1e-6,1e-5,1e-4,1e-3)
  xVals <- unique(sort(c(rat, 1 - rat, 1:(n - 1)/(n))))
  nx <- length(xVals)
  
  repCondVar <- matrix(condVar, ncol = length(condVar), nrow = nx, byrow = T)
  density <- dCopula(cbind(xVals, repCondVar), spVine, h = dists)
  
  # the 1-e7 corners linearily to [0,1], but keep non-negative
  density <- c(max(0,2*density[1]-density[2]),
               density, max(0,2*density[nx]-density[nx-1]))
  linAppr <- approxfun(c(0, xVals, 1), density)
  
  # sum up the denstiy to rescale
  int <- sum(diff(c(0,xVals,1))*(0.5*diff(density)+density[-(nx+2)]))
  return(function(u) linAppr(u)/int)
}

# interpolation

spCopPredict.expectation <- function(predNeigh, spVine, margin, range) {
  stopifnot(!is.null(range))
  stopifnot(is.function(margin$d))
  stopifnot(is.function(margin$p))
  
  predMean <- NULL
  for(i in 1:nrow(predNeigh@data)) { # i <-1
    condSecVine <- condSpVine(as.numeric(predNeigh@data[i,]), predNeigh@distances[i,], spVine)
    
    condExp <-  function(x) {
      condSecVine(margin$p(x))*margin$d(x)*x
    }
    
    predMean <- c(predMean, integrate(condExp,range[1],range[2],subdivisions=1e6)$value)
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
  
  predQuantile <- NULL
  for(i in 1:nrow(predNeigh@data)) { # i <-1
    condSecVine <- condSpVine(as.numeric(predNeigh@data[i,]), predNeigh@distances[i,], spVine)  
    pPred <- optimise(function(x) abs(integrate(condSecVine, 0, x, 
                                                subdivisions=10000L, 
                                                abs.tol=1e-6)$value-p),
                      c(0,1))
    if(pPred$objective > 1e-6)
      warning("Numerical evaluation in predQuantile achieved an obkective of only ",
              pPred$objective, " where 0 has been sought.")
    predQuantile <- c(predQuantile, margin$q(pPred$minimum))
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

spCopPredict <- function(predNeigh, spVine, margin, method="quantile", p=0.5, range=NULL) {
  switch(method,
         quantile=spCopPredict.quantile(predNeigh, spVine, margin, p),
         expectation=spCopPredict.expectation(predNeigh, spVine, margin, range))
}