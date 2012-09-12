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

# calculate critical layer (ONLY 3D by now)
criticalTriple <- function(empCop, cl, u, ind, eps=10e-6) {
  if(!is.matrix(u)) u <- matrix(u,ncol=2)
    
  optimFun <- function(x, u, ind) {
    x <- matrix(rep(x,3),ncol=3)
    x[,ind[1]] <- u[1]
    x[,ind[2]] <- u[2]
    return((empCop(x)-cl)^2*1e6)
  }
  
  res <- apply(u, 1, 
               function(uRow) {
                  upper <- min(1,2+cl-sum(uRow)) # hyperplane in the hypercube
                  if (upper < cl | any(uRow < cl)) return(NA)
                  if (upper == cl) return(cl)
                  opt <- optimize(function(x) optimFun(x, uRow, ind), interval=c(cl,upper))
                  if(opt$objective < eps*1e6) return(opt$minimum)
                  else return(NA)
               })
  return(res)
}


setGeneric("qCopula_u",function(copula,p,u,...) {standardGeneric("qcopula_u")})

qCopula_u.def <- function(copula,p,u,sample=NULL) {
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

setMethod("qCopula_u",signature("copula"),qCopula_u.def)