# adopted from http://www.r-bloggers.com/copulas-and-tail-dependence-part-1/, 04.11.2013

lowerEmpBivTailDepFun <- function(u) {
  stopifnot(ncol(u) == 2)
  empFun <- function(x) sum((u[,1]<=x)&(u[,2]<=x))/sum(u[,1]<=x)
  function(x) sapply(x,empFun)
}

upperEmpBivTailDepFun <- function(u) {
  stopifnot(ncol(u) == 2)
  empFun <- function(x) sum((u[,1]>=x)&(u[,2]>=x))/sum(u[,1]>=x)
  function(x) sapply(x,empFun)
}

empBivTailDepFun <- function(u) {
  stopifnot(ncol(u) == 2)
  
  function(z) {
    res <- z
    res[z>0.5] <- upperEmpBivTailDepFun(u)(z[z>0.5])
    res[z<=0.5] <- lowerEmpBivTailDepFun(u)(z[z<=0.5])
    return(res)
  }
}

##

lowerBivTailDepFun <- function(copula) {
  stopifnot(copula@dimension == 2)
  function(z) pCopula(cbind(z,z),copula)/z
}

upperBivTailDepFun <- function(copula) {
  stopifnot(copula@dimension == 2)
  function(z) (1-2*z+pCopula(cbind(z,z),copula))/(1-z)
}

bivTailDepFun <- function(copula) {
  function(z) {
    res <- z
    res[z>0.5] <- upperBivTailDepFun(copula)(z[z>0.5])
    res[z<=0.5] <- lowerBivTailDepFun(copula)(z[z<=0.5])
    return(res)
  }
}