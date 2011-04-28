#################################################################################
##
##   R package spcopula by Benedikt Gräler Copyright (C) 2011
##
##   This file is part of the R package spcopula.
##
##   The R package spcopula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
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

##########################
##                      ##
## an asymmetric copula ##
##                      ##
##########################
# (see Example 3.16 in: Nelsen, Roger B. (2006): An Introduction to Copulas, second edition, Springer)

# constructor
asCopula <-
function (param)
{
    val <- new("asCopula", dimension = 2, parameters = param, 
        param.names = c("a", "b"), param.lowbnd = c(limA(param[2]), 
            -1), param.upbnd = c(1, 1), message = "asymmetric copula family with cubic and quadratic sections")
    val
}

## density ##

dASC2 <-
function (copula, u) 
{
    a <- copula@parameters[1]
    b <- copula@parameters[2]
    if (!is.matrix(u)) 
        u <- matrix(u, ncol = 2)
    u1 <- u[, 1]
    u2 <- u[, 2]
    return(pmax(a * u2 * (((12 - 9 * u1) * u1 - 3) * u2 + u1 * 
        (6 * u1 - 8) + 2) + b * (u2 * ((u1 * (9 * u1 - 12) + 
        3) * u2 + (12 - 6 * u1) * u1 - 4) - 2 * u1 + 1) + 1, 
        0))
}

setMethod("dcopula", signature("asCopula"), dASC2)

## jcdf ##

pASC2 <-
function (copula, u) 
{
    a <- copula@parameters[1]
    b <- copula@parameters[2]
    if (!is.matrix(u)) 
        u <- matrix(u, ncol = 2)
    u1 <- u[, 1]
    u2 <- u[, 2]
    return(u1 * u2 + u1 * u2 * (1 - u1) * (1 - u2) * ((a - b) * 
        u2 * (1 - u1) + b))
}

setMethod("pcopula", signature("asCopula"), pASC2)

## partial derivatives ##
## ddu

dduASC2 <-
function (copula, pair) 
{
    a <- copula@parameters[1]
    b <- copula@parameters[2]
    if (!is.matrix(pair)) pair <- matrix(pair, ncol = 2)

    u1 <- pair[, 1]
    u2 <- pair[, 2]

    return(u2 * (u2 * (u2 * (u1 * (u1 * (3 * b - 3 * a) + 4 * 
        a - 4 * b) - a + b) + u1 * (u1 * (3 * a - 3 * b) - 4 * 
        a + 6 * b) + a - 2 * b) - 2 * b * u1 + b + 1))
}

setMethod("dducopula", signature("asCopula"),dduASC2)

## ddv

ddvASC2 <-
function (copula, pair) 
{
    a <- copula@parameters[1]
    b <- copula@parameters[2]
    if (!is.matrix(pair)) pair <- matrix(pair, ncol = 2)

    u1 <- pair[, 1]
    u2 <- pair[, 2]

    return(u1 * (-a * (u1 - 1)^2 * u2 * (3 * u2 - 2) + b * (u1 - 1) * (3 * (u1 - 1) * u2^2 - 2 * (u1 - 2) * u2 - 1) + 1))
}

setMethod("ddvcopula", signature("asCopula"),ddvASC2)

## random number generater
# incorporating the inverse of the partial derivative that is solved numerically using optimize

## inverse partial derivative 

invdduASC2 <-
function (copula, u, y) 
{
    if (length(u) != length(y)) 
        stop("Length of u and y differ!")
    res <- NULL
    for (i in 1:length(u)) {
        res <- rbind(res, optimize(function(x) (dduASC2(copula, 
            cbind(rep(u[i], length(x)), x)) - y[i])^2, 
            interval = c(0, 1))$minimum)
    }
    return(res)
}

setMethod("invdducopula", signature("asCopula"),invdduASC2)

## inverse partial derivative ddv
invddvASC2 <-
function (copula, v, y) 
{
    if (!is.matrix(u)) 
        u <- matrix(u, ncol = 2)
    res <- NULL
    for (i in 1:nrow(u)) {
        res <- rbind(res, optimize(function(x) (dduASC2(copula, 
            cbind(x, rep(v[i], length(x)))) - y[i])^2, 
            interval = c(0, 1))$minimum)
    }
    return(res)
}

setMethod("invddvcopula", signature("asCopula"),invddvASC2)


## random number generator
rASC2 <-
function (copula, n) 
{
    u <- matrix(runif(2 * n, min = 0, max = 1), ncol = 2)
    return(cbind(u[, 1], invdduASC2(copula, u)))
}

setMethod("rcopula", signature("asCopula"), rASC2)

## fitment


fitCopulaASC2 <- function (copula, data, method = "itau", start=c(0,0), lower=c(-3.-1), upper=c(1.1), optim.control=list(), optim.method="L-BFGS-B", 
    estimate.variance = FALSE) 
{
    if (method == "ml") 
        fit <- fitASC2.ml(copula, data, start, lower, upper, optim.control, optim.method)
    else if (method == "itau") 
        fit <- fitASC2.itau(copula, data, estimate.variance)
    else if (method == "irho") 
        fit <- fitASC2.irho(copula, data, estimate.variance)
    else stop("Implemented methods for copulas in the spCopula package are: ml, itau, and irho.")
    fit
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
  loglik = sum(log(dcopula(copula,data))),
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
  loglik = sum(log(dcopula(copula,data))),
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

if(length(start)!=2) stop(paste("Start values need to have same length as parameters:",copula@dimension))

optFun <- function(param=c(0,0)) {
      b <- max(min(param[2],1),-1)
      a <- max(min(param[1],1),limA(b))
      return(-sum(log( dASC2(asCopula(c(a,b)),u=data) )))
    }

optimized <- optim(par=start, fn=optFun, method = optim.method, lower=lower, upper=upper,
       control = list())

return(new("fitCopula",
  estimate = optimized$par, 
  var.est = matrix(NA), 
  method = "Numerical MLE over the full range.",
  loglik = -optimized$value,
  convergence = as.integer(NA),
  nsample = nrow(data),
  copula=asCopula(optimized$par)
))
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