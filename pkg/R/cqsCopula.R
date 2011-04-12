#################################################################################
##
##   R package spCopula by Benedikt Gräler Copyright (C) 2009
##
##   This file is part of the R package spCopula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################

######################################################
##                                                  ##
## a symmetric copula with cubic quadratic sections ##
##                                                  ##
######################################################

cqsCopula <-
function (param)
{
    val <- new("cqsCopula", dimension = 2, parameters = param, 
      param.names = c("a", "b"), param.lowbnd = c(limA(param[2]),-1),
      param.upbnd = c(1, 1), message = "copula family with cubic quadratic sections")
    val
}

## density ##

dCQSec <- function (copula, u) 
{
    a <- copula@parameters[1]
    b <- copula@parameters[2]
    if (!is.matrix(u)) 
        u <- matrix(u, ncol = 2)
    u1 <- u[, 1]
    u2 <- u[, 2]
return(pmax(1-b*(1-2*u2)*(1-2*u1)+(b-a)*(1-u2)*(1-3*u2)*(1-u1)*(1-3*u1),0))
}

setMethod("dcopula", signature("cqsCopula"), dCQSec)

## jcdf ##

pCQSec <- function (copula, u) {
    a <- copula@parameters[1]
    b <- copula@parameters[2]
    if (!is.matrix(u)) 
        u <- matrix(u, ncol = 2)
    u1 <- u[, 1]
    u2 <- u[, 2]
return(u1*u2*(1- b*(1-u1)*(1-u2) + (b-a)*(1-u2)^2*(1-u1)^2))
}

setMethod("pcopula", signature("cqsCopula"), pCQSec)

## partial derivatives ##
## ddu

dduCQSec <- function (copula, pair) 
{
    a <- copula@parameters[1]
    b <- copula@parameters[2]
    if (!is.matrix(pair)) pair <- matrix(pair, ncol = 2)

    u1 <- pair[, 1]
    u2 <- pair[, 2]

return( u2*(-a*(u1-1)*(3*u1-1)*(u2-1)^2+b*(u2-1)*(u1*(3*u1*(u2-1)-4*u2+2)+u2)+1) )
}

setMethod("dducopula", signature("cqsCopula"),dduCQSec)

## ddv

ddvCQSec <- function (copula, pair) 
{
    a <- copula@parameters[1]
    b <- copula@parameters[2]
    if (!is.matrix(pair)) pair <- matrix(pair, ncol = 2)

    u1 <- pair[, 1]
    u2 <- pair[, 2]

return( u1*(-a*(u1-1)^2*(u2-1)*(3*u2-1)+b*(u1-1)*(3*(u1-1)*u2^2-4*u1*u2+u1+2*u2)+1) )
}

setMethod("ddvcopula", signature("cqsCopula"),ddvCQSec)

## random number generater
# incorporating the inverse of the partial derivative that is solved numerically using optimize

## inverse partial derivative 

invdduCQSec <-
function (copula, u) 
{
    a <- copula@parameters[1]
    b <- copula@parameters[2]
    if (!is.matrix(u)) 
        u <- matrix(u, ncol = 2)
    u1 <- u[, 1]
    u2 <- u[, 2]
    res <- NULL
    for (i in 1:nrow(u)) {
        res <- rbind(res, optimize(function(x) (dduCQSec(copula, 
            cbind(rep(u[i, 1], length(x)), x)) - u[i, 2])^2, 
            interval = c(0, 1))$minimum)
    }
    return(res)
}

## random generator

rCQSec <-
function (copula, n) 
{
    u <- matrix(runif(2 * n, min = 0, max = 1), ncol = 2)
    return(cbind(u[, 1], invdduCQSec(copula, u)))
}

setMethod("rcopula", signature("cqsCopula"), rCQSec)

## fitment

fitCopula.cqs <- function (copula, data, method = "itau", start=c(0,0), lower=c(-3.-1), upper=c(1.1), optim.control=list(), optim.method="L-BFGS-B", 
    estimate.variance = FALSE) 
{
    if (method == "ml") 
        fit <- fitCQSec.ml(copula, data, start, lower, upper, optim.control, optim.method)
    else if (method == "itau") 
        fit <- fitCQSec.itau(copula, data, estimate.variance)
    else if (method == "irho") 
        fit <- fitCQSec.irho(copula, data, estimate.variance)
    else stop("Implemented methods for copulas in the spCopula package are: ml, itau, and irho.")
    fit
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
  loglik = sum(log(dcopula(copula,data))),
  convergence = as.integer(NA),
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
  loglik = sum(log(dcopula(copula,data))),
  convergence = as.integer(NA),
  nsample = nrow(data),
  copula=copula
))
}

fitCQSec.moa <- function(moa, data, method="itau") {
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
  res <- rbind(res, -sum(log( dCQSec(cqsCopula(c(iFun(param),param)),u=smpl) )))
}
return(res)
}

b <- optimize(sec,c(-1,1))$minimum

param <- c(iFun(b),b)

return(param)
}

####
# maximum likelihood estimation using optim

fitCQSec.ml <- function(copula, data, start, lower, upper, optim.control, optim.method) { 
if(length(start)!=2) stop(paste("Start values need to have same length as parameters:",copula@dimension))

optFun <- function(param=c(0,0)) {
      b <- max(min(param[2],1),-1)
      a <- max(min(param[1],1),limA(b))
      return(-sum(log( dCQSec(cqsCopula(c(a,b)),u=data) )))
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
  copula=cqsCopula(optimized$par)
))
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

setMethod("kendallsTau",signature("cqsCopula"),tauCQSec)

####
# find parameter "a" for parameter "b" under a given measure of association "rho" 
# it may return a value exceeding the limit of "a" which may result in an invalid copula.

iRhoCQSec <- function(b,rho=0) {
return(min(max(limA(b),-3*b - 12*rho),1))
}

####

rhoCQSec <- function(copula){
    a <- copula@parameters[1]
    b <- copula@parameters[2]
return( -(a+3*b)/12 )
}

setMethod("spearmansRho",signature("cqsCopula"),rhoCQSec)