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
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################

## partial derivatives
#######################
# of some copulas from the copula package
# new defined copulas store their partial derivative separately

## Normal Copula
#################

dduNorm <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- qnorm(pair[,1])
  v <- qnorm(pair[,2])

  return(pnorm(v,mean=rho*u,sd=1-rho^2))
}

setMethod("dducopula", signature("normalCopula"),dduNorm)

ddvNorm <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- qnorm(pair[,1])
  v <- qnorm(pair[,2])

  return(pnorm(u,mean=rho*v,sd=1-rho^2))
}

setMethod("ddvcopula", signature("normalCopula"),ddvNorm)

## independent copula
######################

dduIndep <- function(copula, pair){
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)
  return(pair[,2])
}

setMethod("dducopula", signature("indepCopula"),dduIndep)

ddvIndep <- function(copula, pair){
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)
  return(pair[,1])
}

setMethod("ddvcopula", signature("indepCopula"),ddvIndep)

## Clayton Copula
##################

dduClayton <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- pair[,1]
  v <- pair[,2]

  pmax(u^(-rho)+v^(-rho)-1,0)^((-1-rho)/rho)*u^(-rho-1)
}

setMethod("dducopula", signature("claytonCopula"),dduClayton)

ddvClayton <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- pair[,1]
  v <- pair[,2]

  pmax(v^(-rho)+u^(-rho)-1,0)^((-1-rho)/rho)*v^(-rho-1)
}

setMethod("ddvcopula", signature("claytonCopula"),ddvClayton)

## Gumbel Copula 
#################

dduGumbel <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- pair[,1]
  v <- pair[,2]

  pcopula(gumbelCopula(rho),pair) * ((-log(u))^rho+(-log(v))^rho)^(1/rho-1) * (-log(u))^(rho-1)/u
}

setMethod("dducopula", signature("gumbelCopula"),dduGumbel)

ddvGumbel <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- pair[,1]
  v <- pair[,2]

  pcopula(gumbelCopula(rho),pair) * ((-log(v))^rho+(-log(u))^rho)^(1/rho-1) * (-log(v))^(rho-1)/v
}

setMethod("ddvcopula", signature("gumbelCopula"),ddvGumbel)


## Frank Copula 
################

dduFrank <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- pair[,1]
  v <- pair[,2]

  exp(-rho*u)*(exp(-rho*v)-1) / ( (exp(-rho)-1) + (exp(-rho*u)-1)*(exp(-rho*v)-1) )
}

setMethod("dducopula", signature("frankCopula"),dduFrank)

ddvFrank <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- pair[,1]
  v <- pair[,2]

  exp(-rho*v)*(exp(-rho*u)-1) / ( (exp(-rho)-1) + (exp(-rho*v)-1)*(exp(-rho*u)-1) )
}

setMethod("ddvcopula", signature("frankCopula"),ddvFrank)