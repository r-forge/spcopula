## kendall function (empirical) -> spcopula
genEmpKenFun <- function(copula, sample=NULL) {
  if(is.null(sample)) sample <- rCopula(1e6,copula)
  # as empirical copula:
  # copula <- genEmpCop(copula, sample)
  ken <- pCopula(sample, copula)
  
  empKenFun <- function(tlevel) {
    sapply(tlevel,function(t) sum(ken<=t))/nrow(sample)
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

genInvKenFun <- function(kenFun, tol=.Machine$double.eps^.5) {
  invKenFun <- function(k){
    res <- NULL
    for(i in 1:length(k)) {
      res <- c(res, optimize(function(x) (kenFun(x)-k[i])^2, c(0,1), tol=tol)$minimum)
    }
    return(res)
  }
  return(invKenFun)
}

## return periods
kendallRP <- function(kendallFun=NULL, cl=c(.99,.999), mu=1, copula=NULL) {
  if(is.null(kendallFun) & is.null(copula)) stop("Either the kendall distribution function or the copula must be provided. Note that the calculation of the kendall distribution function from the copula is pretty time consuming. Saving them separately might be advantougous.")
  if(is.null(kendallFun)) kendallFun <- genEmpKenFun(copula)
  if(length(mu)>1 & length(cl) > 1) stop("Either the critial level (cl) or mu may be of length larger than 1!")
  return(mu/(1-kendallFun(cl)))
}   

criticalLevel <- function(kendallFun=NULL, KRP=c(100,1000), mu=1, copula=NULL) {
  if(is.null(kendallFun) & is.null(copula)) stop("Either the kendall distribution function or the copula must be provided. Note that the calculation of the kendall distribution function from the copula is pretty time consuming. Saving them separately might be advantougous.")
  if(is.null(kendallFun)) kendallFun <- genEmpKenFun(copula)
  if(length(mu)>1 & length(KRP) > 1) stop("Either the kendall return period or mu may be of length larger than 1!")
  invKenFun <- genInvKenFun(kendallFun)
  return(invKenFun(1-mu/KRP))
}

## next: calculating critical layer, sampling from the layer, selecting "typical" points
# calculate critical layer (ONLY 2D by now)
criticalPair <- function(copula, cl, u, ind, tol=sqrt(.Machine$double.eps)) {
  
  optimFun <- function(x, u, ind) {
    pair <- cbind(x,x)
    pair[,ind] <- u
    return(abs(pCopula(pair,copula)-cl))
  }
  
  sapply(u, function(uRow) {
              upper <- cl+(1-uRow)
              if (upper<cl | uRow < cl) 
                return(NA)
              if (upper == cl) 
                return(cl)
              optimize(function(x) optimFun(x, uRow, ind),
                       interval=c(cl,upper), tol=tol)$minimum
            })
}


# calculate critical layer (ONLY 3D by now)
criticalTriple <- function(copula, cl, u, ind, tol=sqrt(.Machine$double.eps)) {
  if(!is.matrix(u)) u <- matrix(u,ncol=2)
    
  optimFun <- function(x, u, ind) {
    x <- matrix(rep(x,3),ncol=3)
    x[,ind[1]] <- u[1]
    x[,ind[2]] <- u[2]
    return(abs(pCopula(x,copula)-cl))
  }
  
  apply(u, 1, 
        function(uRow) {
          upper <- min(1,2+cl-sum(uRow)) # hyperplane in the hypercube
          if (upper < cl | any(uRow < cl)) 
            return(NA)
          if (upper == cl) 
            return(cl)
          optimize(function(x) optimFun(x, uRow, ind), 
                   interval=c(cl,upper),tol=tol)$minimum
        })
}


setGeneric("qCopula_u",function(copula,p,u,...) {standardGeneric("qCopula_u")})

qCopula_u.def <- function(copula,p,u, tol=.Machine$double.eps^.5) { # sample=NULL
  dim <- copula@dimension
  if(length(p) != length(u)) stop("Length of p and u differ!")
  
  params <- NULL
  for(i in 1:length(p)) { # i <- 1
    if (u[i] < p[i]) {
      params <- rbind(params,rep(NA,dim-1))
    } else {
      if (dim == 2) {
        params <- rbind(params, 
                        optimize(function(v) abs(pCopula(cbind(rep(u[i],length(v)),v),copula)-p[i]),
                                 c(p,1), tol=tol)$minimum)
      } else {
        opt <- optim(par=rep(p[i],dim-1), 
                     function(vw) abs(pCopula(c(u[i],vw), copula)-p[i]), 
                     lower=rep(p[i],dim-1), upper=rep(1,dim-1), method="L-BFGS-B")
        params <- rbind(params, opt$par)
      }
    }
  }
  
  return(cbind(u,params))
}

setMethod("qCopula_u",signature("copula"),qCopula_u.def)