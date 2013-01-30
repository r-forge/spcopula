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
  representation = representation(copulas="list", dimension="integer", 
                                  type="character"),
  validity = validVineCopula,
  contains = list("copula")
)

# constructor
vineCopula <- function (copulas, dim, type) {
  new("vineCopula", copulas=copulas, dimension = as.integer(dim), parameters = numeric(),
      param.names = character(), param.lowbnd = numeric(), 
      param.upbnd = numeric(), type=type, 
      fullname = paste(type, "copula family."))
}

showVineCopula <- function(object) {
  cat(object@fullname, "\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Copulas:\n")
  for (i in (1:length(object@copulas))) {
    cat("  ", class(object@copulas[[i]]), "with parameter(s)", 
        object@copulas[[i]]@parameters, "\n")
  }
}

setMethod("show", signature("vineCopula"), showVineCopula)

## num type

getNumType <- function(copula) {
  if (copula@type == "c-vine") return(1)
  else return(2)  
}

## density ##

## d-vine structure

# copula <- vineFit
# u <- empVine
#   empCopVine

# dDvine(vineFit, empVine,log=T)

dDvine <- function(copula, u, log=FALSE){
  dim <- copula@dimension
  tmp <- u
  u <- NULL
  u[[1]] <- matrix(tmp,ncol=dim)
  
  den <- rep(1,nrow(u[[1]]))
  
  newU <- NULL
  for (i in 1:(dim-1)) {
    tmpCop <- copula@copulas[[i]]
    tmpU <- u[[1]][,i:(i+1)]
    if(log)
      den <- den + dCopula(tmpU, tmpCop,log=T)
    else
      den <- den*dCopula(tmpU,tmpCop,log=F)
    if (i == 1) {
      newU <- cbind(newU, ddvCopula(tmpU, tmpCop))
    } else {
      newU <- cbind(newU, dduCopula(tmpU, tmpCop))
    }
    if (1<i & i<(dim-1)) { 
      newU <- cbind(newU, ddvCopula(tmpU, tmpCop))
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
      if(log)
        den <- den + dCopula(tmpU, tmpCop,log=T)
      else
        den <- den*dCopula(tmpU, tmpCop, log=F)
      if (l < dim-1) {
        if (i == 1) {
          newU <- cbind(newU,ddvCopula(tmpU, tmpCop))
        } else {
          newU <- cbind(newU,dduCopula(tmpU, tmpCop))
        }
        if (1<i & i<(dim-1)) { 
          newU <- cbind(newU,ddvCopula(tmpU, tmpCop))
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
      den <- den*dCopula(tmpU, tmpCop)
      if(l < (dim-1)) newU <- cbind(newU,dduCopula(tmpU, tmpCop))
    }
    if(l < (dim-1)) {
      u[[l+1]] <- newU
      used <- used + dim - l
    }
  }
  
  return(den)
}

##

dvineCopula <- function(u, copula, log=F) { 
  den <- switch(getNumType(copula),dCvine ,dDvine)
  return(den(copula, u, log))
} 

setMethod("dCopula", signature("numeric","vineCopula"), dvineCopula)
setMethod("dCopula", signature("matrix","vineCopula"), dvineCopula)

## jcdf ##
pvineCopula <- function(u, copula) {
  empCop <- genEmpCop(copula,1e5)

  return(pCopula(u, empCop))
}

setMethod("pCopula", signature("numeric","vineCopula"), pvineCopula)
setMethod("pCopula", signature("matrix","vineCopula"), pvineCopula)


## random numbers
linkVineCopSim <- function(n, copula) {
  numType <- getNumType(copula)

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
  if(length(tcops)>0) 
    par2[tcops] <- unlist(lapply(copula@copulas[tcops], function(x) x@df))
  
  return(RVineSim(n, C2RVine(1:copula@dimension, numFam, par1, par2)))
}

setMethod("rCopula", signature("numeric","vineCopula"), linkVineCopSim)