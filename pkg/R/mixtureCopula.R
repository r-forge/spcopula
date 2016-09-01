##############################
##                          ##
## a general mixture copula ##
##                          ##
##############################

# class
setClass("mixtureCopula", contains = "copula", slots = list(memberCops= "list"))

# constructor
mixtureCopula <- function (param = c(0.2, 0.2, 0.5), memberCops = c(normalCopula(), claytonCopula())) {
  stopifnot(length(memberCops) == 2)
  stopifnot(memberCops[[1]]@dimension == memberCops[[2]]@dimension)
  
  cop1.nPar <- length(memberCops[[1]]@parameters)
  cop2.nPar <- length(memberCops[[2]]@parameters)
  
  if (missing(param))
    param <- 0.5
  if (length(param) == 1)
    param <- c(memberCops[[1]]@parameters, memberCops[[2]]@parameters, 0.5)
  else {
    stopifnot(length(param) == cop1.nPar + cop2.nPar + 1)
  
    memberCops[[1]]@parameters <- param[1:cop1.nPar]
    memberCops[[2]]@parameters <- param[(1:cop2.nPar)+cop1.nPar]
  }
  
  new("mixtureCopula", dimension = memberCops[[1]]@dimension, parameters = param, memberCops = memberCops,
      param.names = c(memberCops[[1]]@param.names, memberCops[[2]]@param.names, "mixLambda"),
      param.lowbnd = c(memberCops[[1]]@param.lowbnd, memberCops[[2]]@param.lowbnd, 0),
      param.upbnd = c(memberCops[[1]]@param.upbnd, memberCops[[2]]@param.upbnd, 1), 
      fullname = paste("mixture of a", memberCops[[1]]@fullname, "and a", memberCops[[2]]@fullname))
}

## density ##
setMethod("dCopula", signature(copula = "mixtureCopula"), 
          function(u, copula, log, ...) {
            mixLambda <- tail(copula@parameters, 1)
            res <- (1-mixLambda) * dCopula(u, copula@memberCops[[1]], ...) + mixLambda * dCopula(u, copula@memberCops[[2]], ...)
            if (log)
              return(log(res))
            else 
              return(res)
          })

## jcdf ##
setMethod("pCopula", signature( copula = "mixtureCopula"),
          function(u, copula, ...) {
            mixLambda <- tail(copula@parameters, 1)
            (1-mixLambda) * pCopula(u, copula@memberCops[[1]]) + mixLambda * pCopula(u, copula@memberCops[[2]])
          })

## partial derivatives ##
## ddu

setMethod("dduCopula", signature(copula = "mixtureCopula"),
          function(u, copula, ...) {
            mixLambda <- tail(copula@parameters, 1)
            (1-mixLambda) * dduCopula(u, copula@memberCops[[1]]) + mixLambda * dduCopula(u, copula@memberCops[[2]])
          })

# ddv
setMethod("ddvCopula", signature(copula = "mixtureCopula"),
          function(u, copula, ...) {
            mixLambda <- tail(copula@parameters, 1)
            (1-mixLambda) * ddvCopula(u, copula@memberCops[[1]]) + mixLambda * ddvCopula(u, copula@memberCops[[2]])
          })

## inverse partial derivative 
# invddu
invdduMixCop <- function (u, copula, y) {
  stopifnot(length(u) == length(y)) 
 
  opti <- function(ind) {
    optFun <- function(v) {
      (dduCopula(cbind(u[ind], v), copula) - y[ind])^2
    }
    optimise(optFun, c(0,1))$minimum
  }

  sapply(1:length(y), opti)  
}

setMethod("invdduCopula", 
          signature("numeric", "mixtureCopula", "numeric"), 
          invdduMixCop)

# invddv
invddvMixCop <- function (v, copula, y) {
  stopifnot(length(v) == length(y)) 
  
  opti <- function(ind) {
    optFun <- function(u) {
      (dduCopula(cbind(u, v[ind]), copula) - y[ind])^2
    }
    optimise(optFun, c(0,1))$minimum
  }
  
  sapply(1:length(y), opti)  
}

setMethod("invddvCopula", 
          signature("numeric", "mixtureCopula", "numeric"),
          invddvMixCop)

## random number generator

rMixCop <- function(n, copula, ...) {
  u <- runif(n)
  y <- runif(n)
  
  cbind(u, invdduCopula(u, copula, y))
}

setMethod("rCopula", signature(copula = "mixtureCopula"), rMixCop)

## fitment
fitMixCop <- function(copula, data, start, method="mpl",
                      lower = NULL, upper = NULL, 
                      optim.method = "BFGS", optim.control = list(maxit = 1000), 
                      estimate.variance = FALSE, ...){
  if (missing(start))
    start <- copula@parameters
  stopifnot(method %in% c("ml", "mpl"))
  
  if(is.null(lower))
    lower <- copula@param.lowbnd
  if(is.null(upper))
    upper <- copula@param.lowbnd
    
  copula:::fitCopula.ml(copula, data, start = start, method = method,  
                              lower = lower, upper = upper, 
                              optim.method = optim.method, 
                              optim.control = optim.control, 
                              estimate.variance = estimate.variance , ...)
}

setMethod(fitCopula, 
          signature = c(copula = "mixtureCopula"), 
          fitMixCop)

mixCop <- mixtureCopula(c(0.2,0.5,0.3))
fitCopula(mixCop, rCopula(300, mixCop))

fitMixCop(mixCop, rCopula(300, mixCop))

# 
# fitCopulaASC2 <- function (copula, data, method = "ml", start=c(0,0),
#                            lower=c(-3,-1), upper=c(1,1), 
#                            optim.method="L-BFGS-B", optim.control=list(),
#                            estimate.variance = FALSE) {
#   fit <- switch(method, 
#                 ml=fitASC2.ml(copula, data, start, lower, upper, optim.control, optim.method),
#                 itau=fitASC2.itau(copula, data, estimate.variance),
#                 irho=fitASC2.irho(copula, data, estimate.variance),
#                 stop("Implemented methods for copulas in the spCopula package are: ml, itau, and irho."))
#   return(fit)
# }
# 
# setMethod("fitCopula", signature("asCopula"), fitCopulaASC2)

# setMethod("tau",signature("asCopula"),tauASC2)
# setMethod("rho", signature("asCopula"), rhoASC2)
# setMethod("lambda", signature("asCopula"), 
#           function(copula, ...) c(lower = 0, upper = 0))