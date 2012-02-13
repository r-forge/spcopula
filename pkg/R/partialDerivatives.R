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

# partial derivatives and their inverse of some copulas from the copula package
# new defined copulas store their partial derivative separately

setGeneric("dducopula", function(copula, pair) standardGeneric("dducopula"))
setGeneric("ddvcopula", function(copula, pair) standardGeneric("ddvcopula"))

## inverse partial derivatives 
# numerical standard function
invdducopula <- function(copula, u, y) {
    if (length(u) != length(y)) 
        stop("Length of u and y differ!")
    res <- NULL
    for (i in 1:length(u)) {
        res <- rbind(res, optimize(function(x) (dducopula(copula, 
            cbind(rep(u[i], length(x)), x)) - y[i])^2, 
            interval = c(0, 1))$minimum)
    }
    return(res)
}

setGeneric("invdducopula")

invddvcopula <- function(copula, v, y) {
  standardGeneric("invddvcopula")
    if (length(v) != length(y)) 
        stop("Length of v and y differ!")
    res <- NULL
    for (i in 1:length(v)) {
        res <- rbind(res, optimize(function(x) (ddvcopula(copula, 
            cbind(x, rep(v[i], length(x)))) - y[i])^2, 
            interval = c(0, 1))$minimum)
    }
    return(res)
}

setGeneric("invddvcopula")

###################
## Normal Copula ##
###################

## partial derivative d/du
##########################

dduNorm <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- qnorm(pair[,1]) # u ~ N(0,1)
  v <- qnorm(pair[,2]) # v ~ N(0,1)

  return(pnorm(v,mean=rho*u,sd=1-rho^2))
}

setMethod("dducopula", signature("normalCopula"),dduNorm)

## inverse of the partial derivative d/du
#########################################

invdduNorm <- function(copula, u, y){
cat("This might be wrong!")
  rho <- copula@parameters
  return(pnorm(qnorm(y,mean=rho*u,sd=1-rho^2))) # in doubt
}

setMethod("invdducopula", signature("normalCopula"), invdduNorm)


## partial derivative d/dv
##########################

ddvNorm <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- qnorm(pair[,1])
  v <- qnorm(pair[,2])

  return(pnorm(u,mean=rho*v,sd=1-rho^2))
}

setMethod("ddvcopula", signature("normalCopula"),ddvNorm)

## inverse of the partial derivative d/dv
#########################################

invddvNorm <- function(copula, v, y){
  rho <- copula@parameters
  return(qnorm(y,mean=rho*v,sd=1-rho^2))
}

setMethod("invddvcopula", signature("normalCopula"), invddvNorm)

########################
## independent copula ##
########################

## partial derivative d/du
##########################

dduIndep <- function(copula, pair){
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)
  return(pair[,2])
}

setMethod("dducopula", signature("indepCopula"),dduIndep)

## inverse of the partial derivative d/du
#########################################

invdduIndep <- function(copula, u, y){
  return(y)
}

setMethod("invdducopula", signature("indepCopula"), invdduIndep)

## partial derivative d/dv
##########################

ddvIndep <- function(copula, pair){
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)
  return(pair[,1])
}

setMethod("ddvcopula", signature("indepCopula"),ddvIndep)

## inverse of the partial derivative d/dv
#########################################

invddvIndep <- function(copula, v, y){
  return(y)
}

setMethod("invddvcopula", signature("indepCopula"), invddvIndep)


####################
## Clayton Copula ##
####################

## partial derivative d/du
##########################

dduClayton <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- pair[,1]
  v <- pair[,2]

  pmax(u^(-rho)+v^(-rho)-1,0)^((-1-rho)/rho)*u^(-rho-1) # (u^(-rho)+v^(-rho)-1)^((-1-rho)/rho) * u^(-rho-1)
}

setMethod("dducopula", signature("claytonCopula"), dduClayton)

## inverse of the partial derivative d/du
#########################################

invdduClayton <- function(copula, u, y){
    rho <- copula@parameters[1]
    if (length(u)!=length(y)) 
        stop("Length of u and y differ!")
    return(((y^(rho/(-1-rho))-1)*u^(-rho)+1)^(-1/rho)) # by DL
}

setMethod("invdducopula", signature("claytonCopula"), invdduClayton)

## partial derivative d/dv
##########################

ddvClayton <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- pair[,1]
  v <- pair[,2]

  pmax(v^(-rho)+u^(-rho)-1,0)^((-1-rho)/rho)*v^(-rho-1)
}

setMethod("ddvcopula", signature("claytonCopula"),ddvClayton)

## inverse of the partial derivative d/dv
#########################################

invddvClayton <- function(copula, v, y){
    rho <- copula@parameters[1]
    if (length(v)!=length(y)) 
        stop("Length of v and y differ!")
    return(((y^(rho/(-1-rho))-1)*v^(-rho)+1)^(-1/rho))
}

setMethod("invddvcopula", signature("claytonCopula"), invddvClayton)


###################
## Gumbel Copula ##
###################

## partial derivative d/du
##########################

dduGumbel <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- pair[,1]
  v <- pair[,2]

  pcopula(gumbelCopula(rho),pair) * ((-log(u))^rho+(-log(v))^rho)^(1/rho-1) * (-log(u))^(rho-1)/u
}

setMethod("dducopula", signature("gumbelCopula"),dduGumbel)

## partial derivative d/dv
##########################

ddvGumbel <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- pair[,1]
  v <- pair[,2]

  pcopula(gumbelCopula(rho),pair) * ((-log(v))^rho+(-log(u))^rho)^(1/rho-1) * (-log(v))^(rho-1)/v
}

setMethod("ddvcopula", signature("gumbelCopula"),ddvGumbel)


##################
## Frank Copula ##
##################

## partial derivative d/du
##########################

dduFrank <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- pair[,1]
  v <- pair[,2]

  exp(-rho*u)*(exp(-rho*v)-1) / ( (exp(-rho)-1) + (exp(-rho*u)-1)*(exp(-rho*v)-1) )
}

setMethod("dducopula", signature("frankCopula"),dduFrank)

## inverse of the partial derivative d/du
#########################################

invdduFrank <- function(copula, u, y){
    rho <- copula@parameters[1]
    if (length(u)!=length(y)) 
        stop("Length of u and y differ!")
    return( (-1/rho) * log( y*( exp(-rho)-1)/(exp(-rho*u)-y*(exp(-rho*u)-1)) +1) ) # by DL
}

setMethod("invdducopula", signature("frankCopula"), invdduFrank)

## partial derivative d/dv
##########################

ddvFrank <- function(copula, pair){
  rho <- copula@parameters
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)

  u <- pair[,1]
  v <- pair[,2]

  exp(-rho*v)*(exp(-rho*u)-1) / ( (exp(-rho)-1) + (exp(-rho*v)-1)*(exp(-rho*u)-1) )
}

setMethod("ddvcopula", signature("frankCopula"),ddvFrank)

## inverse of the partial derivative d/dv
#########################################

invddvFrank <- function(copula, v, y){
    rho <- copula@parameters[1]
    if (length(v)!=length(y)) 
        stop("Length of v and y differ!")
    return( (-1/rho) * log( y*( exp(-rho)-1)/(exp(-rho*v)-y*(exp(-rho*v)-1)) +1) )
}

setMethod("invddvcopula", signature("frankCopula"), invddvFrank)

##################
## student Copula ##
##################

## partial derivative d/du
##########################

dduStudent <- function(copula, pair){
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)
  
  df <- copula@df
  v <- qt(pair,df=df)
  
  rho <- copula@parameters[1]
  
  return(pt(sqrt((df+1)/(df+v[,1]^2)) / sqrt(1 - rho^2) * (v[,2] - rho * v[,1]), df=df+1))
}

setMethod("dducopula", signature("tCopula"), dduStudent)

## partial derivative d/dv
##########################

ddvStudent <- function(copula, pair){
  if (!is.matrix(pair)) pair <- matrix(pair,ncol=2)
  
  df <- copula@df
  v <- qt(pair,df=df)
  
  rho <- copula@parameters[1]
  
  return(pt(sqrt((df+1)/(df+v[,2]^2)) / sqrt(1 - rho^2) * (v[,1] - rho * v[,2]), df=df+1))
}

setMethod("ddvcopula", signature("tCopula"), ddvStudent)

## kendall distribution

# empirical default
getKendallDistr <- function(copula, sample=NULL) {
  standardGeneric("getKendallDistr")
  if(is.null(sample)) sample <- rcopula(copula,1e6)
  empCop <- genEmpCop(sample)
  ken <- empCop(sample) # takes really long, any suggestions? Comparring a 1e6x3/1e6x2 matrix by 1e6 pairs/triplets values
  
  empKenFun <- function(tlevel) {
    res <- NULL
    for(t in tlevel) {
      res <- c(res,sum(ken<=t))
    }
    return(res/nrow(sample))
  }
  return(empKenFun)
}

setGeneric("getKendallDistr")

## 

kendallDistribution <- function(copula, t) {
  stop("There is no analytical expression implemented for this copula family. See 'getKendallDstr' for a numerical solution instead.")
}

setGeneric("kendallDistribution")

## Clayton
## kendall distribution/measure, taken from CDVine:::obs.stat
kendall.Clayton <- function(copula, t){
  par = copula@parameters
    
  kt <- rep(NA,length(t))
  kt <- t + t * (1 - t^par)/par
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("claytonCopula"), kendall.Clayton)

setMethod("getKendallDistr", signature("claytonCopula"), 
          function(copula) return(function(t) kendall.Clayton(copula, t)))

## Gumbel
## kendall distribution/measure, taken from CDVine:::obs.stat
kendall.Gumbel <- function(copula, t){
  par = copula@parameters
    
  kt <- rep(NA,length(t))
  kt <- t - t * log(t)/(par)
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("gumbelCopula"), kendall.Gumbel)

setMethod("getKendallDistr", signature("gumbelCopula"), 
          function(copula) return(function(t) kendall.Gumbel(copula, t)))

## Frank
## kendall distribution/measure, taken from CDVine:::obs.stat
kendall.Frank <- function(copula, t){
  par = copula@parameters
    
  kt <- rep(NA,length(t))
  kt <- t + log((1 - exp(-par))/(1 - exp(-par * t))) * (1 - exp(-par * t))/(par * exp(-par * t))
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("frankCopula"), kendall.Frank)

setMethod("getKendallDistr", signature("frankCopula"), 
          function(copula) return(function(t) kendall.Frank(copula, t)))

