######################################################
##                                                  ##
## a symmetric copula with cubic quadratic sections ##
##                                                  ##
######################################################

cqsCopula <- function (param) {
    val <- new("cqsCopula", dimension = as.integer(2), parameters = param, 
      param.names = c("a", "b"), param.lowbnd = c(limA(param[2]),-1),
      param.upbnd = c(1, 1), fullname = "copula family with cubic quadratic sections")
    val
}

## density ##

dCQSec <- function (u, copula) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  
  u1 <- u[, 1]
  u2 <- u[, 2]
  
  return(pmax(1-b*(1-2*u2)*(1-2*u1)+(b-a)*(1-u2)*(1-3*u2)*(1-u1)*(1-3*u1),0))
}

setMethod("dCopula", signature("numeric", "cqsCopula"),
          function(u, copula, ...) {
            dCQSec(matrix(u,ncol=copula@dimension), copula)
          })
setMethod("dCopula", signature("matrix", "cqsCopula"), dCQSec)

## jcdf ##

pCQSec <- function (u, copula) {
    a <- copula@parameters[1]
    b <- copula@parameters[2]
    if (!is.matrix(u)) 
        u <- matrix(u, ncol = 2)
    u1 <- u[, 1]
    u2 <- u[, 2]
return(u1*u2*(1- b*(1-u1)*(1-u2) + (b-a)*(1-u2)^2*(1-u1)^2))
}
setMethod("pCopula", signature("numeric", "cqsCopula"),
          function(u, copula, ...) {
            pCQSec(matrix(u,ncol=copula@dimension), copula)
          })

setMethod("pCopula", signature("matrix","cqsCopula"), pCQSec)

## partial derivatives ##

# solves a*x^3 + b*x^2 + c*x + d = 0
solveCubicEq <- function(a,b,c,d){
eps <- .Machine$double.eps

# using the reduced equation z^3 + 3 * p * z + q = 0 with:
  p <- 3*a*c-b^2
  q <- 2*b^3-9*a*b*c+27*a^2*d
  D <- q^2+4*p^3

  z <- matrix(NA,nrow=length(D),ncol=3)

  ind <- abs(D) <= eps
  if(any(ind)){
    z[ind,1] <- 0.5*(-4*q[ind])^(1/3)
    z[ind,2] <- -z[ind,1]
  }

  ind <- D > eps
  if(any(ind)){
    cubeRad <- -4*q[ind]+4*sqrt(D[ind])
    r1 <- sign(cubeRad)*abs(cubeRad)^(1/3)
    cubeRad <- -4*q[ind]-4*sqrt(D[ind])
    r2 <- sign(cubeRad)*abs(cubeRad)^(1/3)
    z[ind,1] <- 0.5*(r1+r2)
  }

  ind <- D < eps
  if(any(ind)){
    phi <- acos(-q[ind]/(2*sqrt(-p[ind]^3)))
    triple <- NULL
    triple <- sqrt(-p[ind])*2*cos(phi/3)
    triple <- cbind(triple,sqrt(-p[ind])*2*cos((phi+2*pi)/3))
    triple <- cbind(triple,sqrt(-p[ind])*2*cos((phi+4*pi)/3))
    z[ind,] <- triple
  }

return((z-b)/(3*a))
}

## partial derivative ddu pCQSec

dduCQSec <- function (u, copula) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]

  u1 <- u[, 1]
  u2 <- u[, 2]

  return(u2-b*(u2-u2^2-2*u1*u2+2*u1*u2^2)+(b-a)*(u2-4*u1*u2+3*u1^2*u2-2*u2^2+8*u1*u2^2-6*u1^2*u2^2+u2^3-4*u1*u2^3+3*u1^2*u2^3))
}

setMethod("dduCopula", signature("numeric","cqsCopula"),
          function(u, copula, ...) {
            dduCQSec(matrix(u,ncol=copula@dimension), copula)
          }) 
setMethod("dduCopula", signature("matrix","cqsCopula"), dduCQSec)

## inverse partial derivative ddu

## inverse partial derivative ddu
# seems to be accurate (1.4e-05 is the max out of 1000 random CQSec-copulas for 1000 random pairs (u,v) each.)
invdduCQSec <- function (u, copula, y) {
    if (length(u)!=length(y)) 
        stop("Length of u and y differ!")

    a <- copula@parameters[1]
    b <- copula@parameters[2]

# solving the cubic equation: u^3 * c3 + u^2 * c2 + u * c1 + c0 = 0
    usq <- u^2
    c3 <- (b-a)*(1-4*u+3*usq)
    c2 <- (b-a)*(-2+8*u-6*u^2)-b*(-1+2*u)
    c1 <- (b-a)*(1-4*u+3*u^2)-b*(1-2*u)+1
    c0 <- -y

v <- solveCubicEq(c3,c2,c1,c0)

filter <- function(vec){
  vec <- vec[!is.na(vec)]
  return(vec[vec >= 0 & vec <= 1])
}

return(apply(v,1,filter))
}

setMethod("invdduCopula", signature("numeric","cqsCopula","numeric"), invdduCQSec)

## partial derivative ddv

ddvCQSec <- function (u, copula) {
    a <- copula@parameters[1]
    b <- copula@parameters[2]
    if (!is.matrix(u)) u <- matrix(u, ncol = 2)

    u1 <- u[, 1]
    u2 <- u[, 2]

return(u1-b*(u1-2*u1*u2-u1^2+2*u1^2*u2)+(b-a)*(u1-2*u1^2+u1^3-4*u1*u2+8*u1^2*u2-4*u1^3*u2+3*u1*u2^2-6*u1^2*u2^2+3*u1^3*u2^2))
}

setMethod("ddvCopula", signature("numeric","cqsCopula"),
          function(u, copula, ...) {
            ddvCQSec(matrix(u,ncol=copula@dimension), copula)
          })
setMethod("ddvCopula", signature("matrix","cqsCopula"), ddvCQSec)

## inverse partial derivative ddv
# seems to be accurate (1e-05 is the max out of 5000 random CQSec-copulas for 1000 random pairs (u,v) each. Very most are below 10*.Machine$double.eps)
invddvCQSec <- function (v, copula, y) {
    if (length(v)!=length(y)) 
        stop("Length of v and y differ!")

    a <- copula@parameters[1]
    b <- copula@parameters[2]

# solving the cubic equation: u^3 * c3 + u^2 * c2 + u * c1 + c0 = 0
    vsq <- v^2
    c3 <- (b-a)*(1-4*v+3*vsq)
    c2 <- (b-a)*(-2+8*v-6*vsq)-b*(-1+2*v)
    c1 <- (b-a)*(1-4*v+3*vsq)-b*(1-2*v)+1
    c0 <- -y

u <- solveCubicEq(c3,c2,c1,c0)

filter <- function(vec){
  vec <- vec[!is.na(vec)]
  return(vec[vec >= 0 & vec <= 1])
}

return(apply(u,1,filter))
}

setMethod("invddvCopula", signature("numeric","cqsCopula","numeric"), invddvCQSec)

## random number generator

rCQSec <- function (n, copula) {
  u <- runif(n, min = 0, max = 1)
  y <- runif(n, min = 0, max = 1)
    
  res <- cbind(u, invdduCQSec(u, copula, y))
  colnames(res) <- c("u","v")
    
  return(res)
}

setMethod("rCopula", signature("numeric","cqsCopula"), rCQSec)
## fitment

fitCopula.cqs <- function (copula, data, method = "ml", start=c(0,0), 
                           lower=c(-3,-1), upper=c(1,1), 
                           optim.method="L-BFGS-B", optim.control=list(),
                           estimate.variance = FALSE) {
  fit <- switch(method,
                ml=fitCQSec.ml(copula, data, start, lower, upper, optim.control, optim.method),
                itau=fitCQSec.itau(copula, data, estimate.variance),
                irho=fitCQSec.irho(copula, data, estimate.variance),
                stop("Implemented methods for copulas in the spCopula package are: ml, itau, and irho."))
  return(fit)
}

setMethod("fitCopula", signature("cqsCopula"), fitCopula.cqs)

## Fits the copula with cubic and quadratic sections acoording to a measure of association.
## It performs a maximum likelihood evaluation over all possible pairs of 
## parameters a and b generating a copula with the given mesaure of 
## association.
#
# moa
#  measure of association, according to method
# data
#  the bivariate data set as a 2-column matrix within the unitsquare
# method
#  one of kendall or spearman according to the calculation of moa

fitCQSec.itau <- function(copula, data, estimate.variance) {
tau <- cor(data,method="kendall")[1,2]
esti <- fitCQSec.moa(tau, data, method="itau")
copula <- cqsCopula(esti)
return(new("fitCopula",
  estimate = esti, 
  var.est = matrix(NA), 
  method = "Inversion of Kendall's tau and MLE",
  loglik = sum(log(dCopula(data, copula))),
  fitting.stats=list(convergence = as.integer(NA)),
  nsample = nrow(data),
  copula=copula
))
}

fitCQSec.irho <- function(copula, data, estimate.variance){
rho <- cor(data,method="spearman")[1,2]
esti <- fitCQSec.moa(rho, data, method="irho")
copula <- cqsCopula(esti)
return(new("fitCopula",
  estimate = esti, 
  var.est = matrix(NA), 
  method = "Inversion of Spearman's rho and MLE",
  loglik = sum(log(dCopula(data, copula))),
  fitting.stats=list(convergence = as.integer(NA)),
  nsample = nrow(data),
  copula=copula
))
}

fitCQSec.moa <- function(moa, data, method="itau", tol=.Machine$double.eps^.5) {
smpl <- as.matrix(data)

iTau <- function(p) {
  iTauCQSec(p,moa)
}

iRho <- function(p) {
  iRhoCQSec(p,moa)
}

iFun <- switch(method, itau=iTau, irho=iRho)

sec <- function (parameters) {
res <- NULL
for(param in parameters) {
  res <- rbind(res, -sum(log( dCQSec(smpl, cqsCopula(c(iFun(param),param))) )))
}
return(res)
}

b <- optimize(sec,c(-1,1), tol=tol)$minimum

param <- c(iFun(b),b)

return(param)
}

# maximum log-likelihood estimation of a and b using optim

fitCQSec.ml <- function(copula, data, start, lower, upper, optim.control, optim.method) { 
  if(length(start)!=2) stop("Start values need to have same length as parameters.")
  
  optFun <- function(param=c(0,0)) {
    if(any(param > 1) | param[2] < -1 | param[1] < limA(param[2])) return(1)
    return(-sum(log( dCQSec(data, cqsCopula(param)) )))
  }
  
  optimized <- optim(par=start, fn=optFun, method = optim.method, 
                     lower=lower, upper=upper, control = optim.control)
  
  return(new("fitCopula", estimate = optimized$par, var.est = matrix(NA),
             method = "Numerical MLE over the full range.",
             loglik = -optimized$value, fitting.stats= optimized,
             nsample = nrow(data), copula=cqsCopula(optimized$par)))
}

####

iTauCQSec <- function(b,tau=0) {
return(min(max(limA(b),(b^2 + 75*b + 450*tau)/(b - 25)),1))
}

####

tauCQSec <- function(copula){
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  return( (a*b - 25*a - b^2 - 75*b)/450 )
}

setMethod("tau",signature("cqsCopula"),tauCQSec)

####
# find parameter "a" for parameter "b" under a given measure of association "rho" 
# it may return a value exceeding the limit of "a" which may result in an invalid copula.

iRhoCQSec <- function(b, rho=0) {
  return(min(max(limA(b),-3*b - 12*rho),1))
}

####

rhoCQSec <- function(copula){
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  return( -(a+3*b)/12 )
}

setMethod("rho",signature("cqsCopula"),rhoCQSec)