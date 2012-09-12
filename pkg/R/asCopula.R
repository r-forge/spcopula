#################################################################################
##
##   R package spCopula by Benedikt Gräler Copyright (C) 2011
##
##   This file is part of the R package spCopula.
##
##   The R package spCopula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package spCopula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package spCopula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################

##########################
##                      ##
## an asymmetric copula ##
##                      ##
##########################
# (see Example 3.16 in: Nelsen, Roger B. (2006): An Introduction to Copulas, second edition, Springer)

# constructor
asCopula <- function (param) {
  val <- new("asCopula", dimension = as.integer(2), parameters = param, 
             param.names = c("a", "b"), param.lowbnd = c(limA(param[2]), -1),
             param.upbnd = c(1, 1), fullname = "asymmetric copula family with cubic and quadratic sections")
  return(val)
}

## density ##

dASC2 <- function (u, copula) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  
  u1 <- u[, 1]
  u2 <- u[, 2]
  
  return(pmax(a * u2 * (((12 - 9 * u1) * u1 - 3) * u2 + u1 * (6 * u1 - 8) + 2) + b * (u2 * ((u1 * (9 * u1 - 12) + 3) * u2 + (12 - 6 * u1) * u1 - 4) - 2 * u1 + 1) + 1,0))
}

setMethod("dCopula", signature("numeric","asCopula"), 
          function(u, copula, ...) {
            dASC2(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dCopula", signature("matrix","asCopula"), dASC2)

## jcdf ##
pASC2 <- function (u, copula) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]

  u1 <- u[, 1]
  u2 <- u[, 2]
  return( u1 * u2 + u1 * u2 * (1 - u1) * (1 - u2) * ((a - b) * u2 * (1 - u1) + b) )
}

setMethod("pCopula", signature("numeric", "asCopula"),
          function(u, copula, ...) {
            pASC2(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("pCopula", signature("matrix", "asCopula"), pASC2)

## partial derivatives ##
## ddu

dduASC2 <- function (u, copula) {
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  u1 <- u[, 1]
  u2 <- u[, 2]

  return(u2*(1 + b*(-1 + 2*u1)*(-1 + u2) - (a - b)*(1 - 4*u1 + 3*u1^2)*(-1 + u2)*u2))
}

setMethod("dduCopula", signature("numeric", "asCopula"),
          function(u, copula, ...) {
            dduASC2(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix", "asCopula"), dduASC2)

## ddv
ddvASC2 <- function (u, copula){
  a <- copula@parameters[1]
  b <- copula@parameters[2]

  u1 <- u[, 1]
  u2 <- u[, 2]

  return( u1 + b*(-1 + u1)*u1*(-1 + 2*u2) - (a - b)*(-1 + u1)^2*u1*u2*(-2 + 3*u2))
}

setMethod("ddvCopula", signature("numeric", "asCopula"),
          function(u, copula, ...) {
            ddvASC2(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix", "asCopula"),ddvASC2)

## random number generater
# incorporating the inverse of the partial derivative that is solved numerically using optimize

## inverse partial derivative 

invdduASC2 <- function (u, copula, y) {
    if (length(u)!=length(y)) 
        stop("Length of u and y differ!")

    a <- copula@parameters[1]
    b <- copula@parameters[2]

# solving the cubic equation: u^3 * c3 + u^2 * c2 + u * c1 + c0 = 0
    usq <- u^2
    c3 <- (a-b)*(-3*usq+4*u-1)
    c2 <- (a-b)*(1-4*u+3*usq)+b*(- 1 + 2*u)
    c1 <- 1+b*(1-2*u)
    c0 <- -y

v <- solveCubicEq(c3,c2,c1,c0) # from cqsCopula.R

filter <- function(vec){
  vec <- vec[!is.na(vec)]
  return(vec[vec >= 0 & vec <= 1])
}

return(apply(v,1,filter))
}

setMethod("invdduCopula", signature("numeric","asCopula","numeric"),invdduASC2)

## inverse partial derivative ddv
invddvASC2 <- function (v, copula, y) {
    if (length(v)!=length(y)) 
        stop("Length of v and y differ!")

    a <- copula@parameters[1]
    b <- copula@parameters[2]

# solving the cubic equation: u^3 * c3 + u^2 * c2 + u * c1 + c0 = 0
    vsq <- v^2
    c3 <- (a-b)*(2*v-3*vsq)
    c2 <- (a-b)*(-4*v+6*vsq)+b*(-1+2*v)
    c1 <- 1+(a-b)*(2*v - 3*vsq)+b*(1-2*v)
    c0 <- -y

u <- solveCubicEq(c3,c2,c1,c0) # from cqsCopula.R

filter <- function(vec){
  vec <- vec[!is.na(vec)]
  return(vec[vec >= 0 & vec <= 1])
}

return(apply(u,1,filter))
}

setMethod("invddvCopula", signature("numeric","asCopula","numeric"),invddvASC2)

## random number generator
rASC2 <- function (n, copula) {
    u <- runif(n, min = 0, max = 1)
    y <- runif(n, min = 0, max = 1)
    return(cbind(u, invdduASC2(u, copula, y) ))
}

setMethod("rCopula", signature("numeric", "asCopula"), rASC2)

## fitment

fitCopulaASC2 <- function (copula, data, method = "ml", start=c(0,0),
                           lower=c(-3,-1), upper=c(1,1), 
                           optim.method="L-BFGS-B", optim.control=list(),
                           estimate.variance = FALSE) {
  fit <- switch(method, 
                ml=fitASC2.ml(copula, data, start, lower, upper, optim.control, optim.method),
                itau=fitASC2.itau(copula, data, estimate.variance),
                irho=fitASC2.irho(copula, data, estimate.variance),
                stop("Implemented methods for copulas in the spCopula package are: ml, itau, and irho."))
  return(fit)
}

setMethod("fitCopula", signature("asCopula"), fitCopulaASC2)

## Fits the type 2 asymmetric copula acoording to a measure of association.
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

fitASC2.itau <- function(copula, data, estimate.variance) {
tau <- cor(data,method="kendall")[1,2]
esti <- fitASC2.moa(tau, data, method="itau")
copula <- asCopula(esti)
return(new("fitCopula",
  estimate = esti, 
  var.est = matrix(NA), 
  method = "Inversion of Kendall's tau and MLE",
  loglik = sum(log(dCopula(data, copula))),
  convergence = as.integer(NA),
  nsample = nrow(data),
  copula=copula
))
}

fitASC2.irho <- function(copula, data, estimate.variance){
rho <- cor(data,method="spearman")[1,2]
esti <- fitASC2.moa(rho, data, method="itau")
copula <- asCopula(esti)
return(new("fitCopula",
  estimate = esti, 
  var.est = matrix(NA), 
  method = "Inversion of Spearman's rho and MLE",
  loglik = sum(log(dCopula(data, copula))),
  convergence = as.integer(NA),
  nsample = nrow(data),
  copula=copula
))
}

fitASC2.moa <- function(moa, data, method="itau") {
  smpl <- as.matrix(data)

  iTau <- function(p) {
    iTauASC2(p,moa)
  }

  iRho <- function(p) {
    iRhoASC2(p,moa)
  }

  iFun <- switch(method, itau=iTau, irho=iRho)

  sec <- function (parameters) {
    res <- NULL
    for(param in parameters) {
      res <- rbind(res, -sum(log( dASC2(asCopula(c(iFun(param),param)),u=smpl) )))
    }
    return(res)
  }

  b <- optimize(sec,c(-1,1))$minimum

  param <- c(iFun(b),b)

  return(param)
}

# maximum log-likelihood estimation of a and b using optim

fitASC2.ml <- function(copula, data, start, lower, upper, optim.control, optim.method) { 
  if(length(start)!=2) stop("Start values need to have same length as parameters:")
  
  optFun <- function(param=c(0,0)) {
    if(any(param > 1) | param[2] < -1 | param[1] < limA(param[2])) return(1)
    return(-sum(log( dASC2(asCopula(param),u=data))))
  }
  
  optimized <- optim(par=start, fn=optFun, method = optim.method, 
                     lower=lower, upper=upper, control = optim.control)

  return(new("fitCopula", estimate = optimized$par, var.est = matrix(NA), 
             method = "Numerical MLE over the full range.", 
             loglik = -optimized$value, convergence = optimized$convergence,
             nsample = nrow(data), copula=asCopula(optimized$par)))
}

####

iTauASC2 <- function(b,tau=0) {
return(min(max(limA(b),(450*tau-75*b+b^2)/(25-b)),1))
}

####

tauASC2 <- function(copula){
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  return((75*b-b^2+a*(25-b))/450)
}

setMethod("kendallsTau",signature("asCopula"),tauASC2)

####
# find parameter "a" for parameter "b" under a given measure of association "rho" 
# it may return a value exceeding the limit of "a" which may result in an invalid copula.

iRhoASC2 <- function(b,rho=0) {
  return(min(max(limA(b),12*rho-3*b),1))
}

####

rhoASC2 <- function(copula){
  a <- copula@parameters[1]
  b <- copula@parameters[2]
  
  return((a+3*b)/12)
}

setMethod("spearmansRho",signature("asCopula"),tauASC2)