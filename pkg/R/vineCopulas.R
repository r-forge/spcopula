####################
##                ##
##  vine copulas  ##
##                ##
####################

validVineCopula = function(object) {
  dim <- object@dimension
  if( dim <= 2)
    return("Number of dimension too small (>2).")
  if(length(object@copulas)!=(dim*(dim-1)/2))
    return("Number of provided copulas does not match given dimension.")
  if(!any(unlist(lapply(object@copulas,function(x) is(x,"copula")))))
    return("Not all provided copulas in your list are indeed copulas.")
  if(!(object@type == "c-vine" | object@type == "d-vine"))
    return("Only c-vines and d-vines are implemented.")
  else return (TRUE)
}

setClass("vineCopula",
  representation = representation(copulas="list", dimension="numeric", type="character", pdf="numeric"),
  validity = validVineCopula,
  contains = list("copula")
)

# constructor
vineCopula <- function (copulas, dim, type) {
    val <- new("vineCopula", copulas=copulas, dimension = dim, parameters = numeric(), 
        param.names = character(), param.lowbnd = numeric(), param.upbnd = numeric(), type=type, pdf=numeric(), message = paste(type, "copula family."))
    val
}

showVineCopula <- function(object) {
  cat(object@message, "\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Copulas:\n")
  for (i in (1:length(object@copulas))) cat("  ", class(object@copulas[[i]]), "with parameter(s)", object@copulas[[i]]@parameters, "\n")
}

setMethod("show", signature("vineCopula"), showVineCopula)

## density ##
dvineCopula <- function(copula, u) { 
  dim <- copula@dimension
  tmp <- u
  u <- NULL
  u[[1]] <- matrix(tmp,ncol=dim)

  den <- rep(1,nrow(u[[1]]))

# d-vine method:
  # Only the d-vine structure is implemented by now!
  newU <- NULL
  for (i in 1:(dim-1)) {
    tmpCop <- copula@copulas[[i]]
    tmpU <- u[[1]][,i:(i+1)]
    den <- den*dcopula(tmpCop,tmpU)
    if (i == 1) {
      newU <- cbind(newU, ddvcopula(tmpCop,tmpU))
    } else {
      newU <- cbind(newU, dducopula(tmpCop,tmpU))
    }
    if (1<i & i<(dim-1)) { 
      newU <- cbind(newU, ddvcopula(tmpCop,tmpU))
    }
  }
  u[[2]] <- newU

  used <- dim-1

  for (l in 2:(dim-1)) {
    newU <- NULL
    for (i in 1:(dim-l)) {
      tmpCop <- copula@copulas[[used+i]]
      tmpU <- u[[l]][,(i*2-1):(i*2)]
      den <- den*dcopula(tmpCop, tmpU)
      if (i == 1) {
	newU <- cbind(newU,ddvcopula(tmpCop,tmpU))
      } else {
	newU <- cbind(newU,dducopula(tmpCop,tmpU))
      }
      if (1<i & i<(dim-1)) { 
	newU <- cbind(newU,ddvcopula(tmpCop,tmpU))
      }
    }
    u[[l+1]] <- newU
    used <- used + dim - l + 1
  }

return(den)
} 

setMethod("dcopula", signature("vineCopula"), dvineCopula)


## jcdf ##
# calcGrid <- function(dim,p=0) {
#   steps <- round(1e7^(1/dim),0)
#   origSeq <- seq(p+(1-p)/steps,1,length.out=steps)-(1-p)/2/steps
#   grid <- matrix(origSeq,ncol=1)
#   for (i in 1:(dim-1)) {
#     grid <- grid[rep(1:nrow(grid),steps),]
#     grid <- cbind(grid, rep(origSeq,each=steps^i))
#   }
# return(grid)
# }
# 
# setGeneric("setPDF<-", function(copula, value) standardGeneric("setPDF<-"))
# 
# repPDF <- function(copula, value) {
#   cat("Note: The calculation is a simple numerical approximation based on 1e6 regular (hyper)cubes.\n")
#   copula@pdf <- dcopula(vine3d,calcGrid(vine3d@dimension))
#   gc()
#   return(copula)
# }

# setReplaceMethod("setPDF", "vineCopula", repPDF)
# 
# setGeneric("forcePDF<-", function(copula, value) standardGeneric("forcePDF<-"))
# 
# forcePDF <- function(copula, value) {
#   copula@pdf <- value
#   return(copula)
# }
# 
# setReplaceMethod("forcePDF", "vineCopula", forcePDF)

genEmpCop <- function(data) {
  empCop <- function(u) {
    u <- matrix(u,ncol=ncol(data))

# make this a C-function
    res <- NULL
    for(i in 1:nrow(u)) {
      bool <- t(t(data) <= u[i,])
      for (i in 2:ncol(data)) bool[,1] <- bool[,1] * bool[,i]
      res <- c(res,sum(bool[,1]))
    }
    return(res/nrow(data))
  }
  return(empCop)
}

# genEmpCop <- function(data) {
#   empCop <- function(u) {
#     u <- matrix(u,ncol=ncol(data))
#     res <- NULL
#     fun <- function(x){
#       bool <- t(t(data) <= x)
#       for (i in 2:ncol(data)) bool[,1] <- bool[,1] * bool[,i]
#       return(sum(bool[,1]))
#     }
#     return(apply(u,1,fun)/nrow(data))
#   }
#   return(empCop)
# }

pvineCopula <- function(copula, u) {
  cat("Note: the copula is empirically evaluated from 1M samples.")
  empCop <- genEmpCop(rcopula(copula,1e6))

  return(empCop(u))
}

setMethod("pcopula", signature("vineCopula"), pvineCopula)

#   dim <- copula@dimension
# #  grid <- calcGrid(dim)
#   pvineCubed <- copula@pdf
#   if (length(pvineCubed) == 0) {
#     stop("The Copula density surface has not yet been approximated. Call setPDF(yourCopula) <- 1 once.\n")
#   }
# 
#   if(is.null(ncol(u)) || ncol(u) != dim) u <- matrix(u,ncol=dim)
#   res <- NULL
#   sumTrues <- function(u) {
#     bools <- (grid[,1] <= u[1] & grid[,2] <= u[2] & grid[,3] <= u[3])
#     return(sum(pvineCubed[bools],rm.na=T))
#   }
#   res <- apply(u,1,sumTrues)*diff(grid[1:2,1])^dim # adjust rescaling
#   return(res)

## random numbers

linkCDVineSim <- function(copula, n) {
  if (copula@type == "c-vine") { 
    numType <- 1
  } else {
    numType <- 2
  }
  par1 <- unlist(lapply(copula@copulas,function(x) x@parameters[1]))
  par2 <- unlist(lapply(copula@copulas,function(x) x@parameters[2]))
  par2[is.na(par2)] <- 0
  return(CDVine::CDVineSim(n,unlist(lapply(copula@copulas,function(x) x@family)),par1,par2,numType))
}

setMethod("rcopula", signature("vineCopula"), linkCDVineSim)

setGeneric("qcopula_u",function(copula,p,u,...) {standardGeneric("qcopula_u")})

qCopula_u <- function(copula,p,u,sample=NULL) {
  dim <- copula@dimension
  if(length(p) != length(u)) stop("Length of p and u differ!")
  if(is.null(sample)) sample <- rcopula(copula,1e6)
  empCop <- genEmpCop(sample)
  params <- NULL
  
  for(i in 1:length(p)) {
    if (u[i] < p[i]) {
      params <- rbind(params,rep(NA,dim-1))
    } else {
    if (dim == 2) {
      params <- rbind(params,optimize(function(v) (empCop(cbind(u[i],v))-p[i])^2,c(p,1))$minimum)
    } else {
      opt <- optim(par=rep(p[i],dim-1), function(vw) (empCop(c(u[i],vw))-p[i])^2, lower=rep(p[i],dim-1), upper=rep(1,dim-1), method="L-BFGS-B")
      params <- rbind(params, opt$par)
    }
    }
  }

  return(cbind(u,params))
}

setMethod("qcopula_u",signature("copula"),qCopula_u)

## kendall function (empirical) -> spcopula

genEmpKenFun <- function(copula, sample=NULL) {
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

## inverse kendall function
# Kendall return period:
# K_x(t)= \mu / (1-K_C(t))
# 
# solve K_x(t)=1/1000
# 
# KRP = \mu / (1-K_C(t))
# <=> (1-K_C(t)) = \mu / KRP
# <=> K_C(t) = 1 - \mu /KRP
# <=> t = K_C^{-1}(1 - \mu /KRP)

genInvKenFun <- function(kenFun, ...) {
  invKenFun <- function(k){
    res <- NULL
    for(i in 1:length(k)) {
      res <- c(res, optimize(function(x) (kenFun(x)-k[i])^2,c(0,1))$minimum)
    }
    return(res)
  }
  return(invKenFun)
}

# qcopula_u(vine3d,.9,.9)
# 
# curve(empCop( cbind(c(x),.99,.99)),from=0,to=1)
# 
# empCop( cbind(runif(10), .7,.4))
# 
# empVine <- rcopula(vine3d,1e6)
# empCop <- empCopula(empVine)
# 
# 
# 
