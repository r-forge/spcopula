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
  representation = representation(copulas="list", dimension="integer", type="character", pdf="numeric"),
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

## num type

getNumType <- function(copula) {
  if (copula@type == "c-vine") return(1)
  else return(2)  
}

## density ##

## d-vine structure

dDvine <- function(copula, u){
  dim <- copula@dimension
  tmp <- u
  u <- NULL
  u[[1]] <- matrix(tmp,ncol=dim)
  
  den <- rep(1,nrow(u[[1]]))
  
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
#       cat(used+i,"\n")
      tmpCop <- copula@copulas[[used+i]]
      tmpU <- u[[l]][,(i*2-1):(i*2)]
      den <- den*dcopula(tmpCop, tmpU)
      if (l < dim-1) {
        if (i == 1) {
          newU <- cbind(newU,ddvcopula(tmpCop,tmpU))
        } else {
          newU <- cbind(newU,dducopula(tmpCop,tmpU))
        }
        if (1<i & i<(dim-1)) { 
          newU <- cbind(newU,ddvcopula(tmpCop,tmpU))
        }
      } 
    }
    u[[l+1]] <- newU
    used <- used + dim - l
  }
  
  return(den)
}

## c-vine structure

dCvine <- function(copula, u) {
#   cat("c-vine \n")
  dim <- copula@dimension
  tmp <- u
  u <- NULL
  u[[1]] <- matrix(tmp,ncol=dim)
  
  den <- rep(1,nrow(u[[1]]))

  used <- 0 # already used copulas
  
  for (l in 1:(dim-1)) {
    newU <- NULL
    for (i in 1:(dim-l)) {
#       cat(used+i,"\n")
      tmpCop <- copula@copulas[[used+i]]
      tmpU <- u[[l]][,c(1,(i+1))]
      den <- den*dcopula(tmpCop, tmpU)
      if(l < (dim-1)) newU <- cbind(newU,dducopula(tmpCop,tmpU))
    }
    if(l < (dim-1)) {
      u[[l+1]] <- newU
      used <- used + dim - l
    }
  }
  
  return(den)
}

##

dvineCopula <- function(copula, u) { 
  den <- switch(getNumType(copula),dCvine ,dDvine)
  return(den(copula, u))
} 

setMethod("dcopula", signature("vineCopula"), dvineCopula)


## jcdf ##
genEmpCop <- function(data) {
  t_data <- t(data)

  empCop <- function(u) {
    u <- matrix(u,ncol=nrow(t_data))

# --/-- make this a C-function?
    res <- NULL
    for(i in 1:nrow(u)) {
      bool <- t_data <= u[i,]
      for (i in 2:nrow(t_data)) bool[1,] <- bool[1,] * bool[i,]
      res <- c(res,sum(bool[1,]))
    }
# --//--

    return(res/ncol(t_data))
  }
  return(empCop)
}

pvineCopula <- function(copula, u) {
  cat("Note: the copula is empirically evaluated from 100.000 samples.")
  empCop <- genEmpCop(rcopula(copula,1e5))

  return(empCop(u))
}

setMethod("pcopula", signature("vineCopula"), pvineCopula)


## random numbers
linkCDVineSim <- function(copula, n) {
  numType <- getNumType

  getFamily <- function(copula) {
    if("family" %in% slotNames(copula)) numFam <- copula@family
    else {
      numFam <- switch(class(copula)[1], normalCopula=1, tCopula=2, claytonCopula=3, gumbelCopula=4, frankCopula=5)
    }
  }

  par1 <- unlist(lapply(copula@copulas,function(x) x@parameters[1]))
  par2 <- unlist(lapply(copula@copulas,function(x) x@parameters[2]))
  par2[is.na(par2)] <- 0
  numFam <- unlist(lapply(copula@copulas,getFamily))
  tcops <- which(numFam==2) #? length(which(5==3))
  if(length(tcops)>0) par2[tcops] <- unlist(lapply(copula@copulas[tcops], function(x) x@df))
  return(CDVine::CDVineSim(n,numFam,par1,par2,numType))
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

