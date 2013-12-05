# adopted from http://www.r-bloggers.com/copulas-and-tail-dependence-part-1/, 04.11.2013

lowerEmpTailDepFun <- function(u) {
  empFun <- function(x) sum((u[,1]<=x)&(u[,2]<=x))/sum(u[,1]<=x)
  function(x) sapply(x,empFun)
}

upperEmpTailDepFun <- function(u) {
  empFun <- function(x) sum((u[,1]>=x)&(u[,2]>=x))/sum(u[,1]>=x)
  function(x) sapply(x,empFun)
}

empTailDepFun <- function(u) {
  function(z) {
    res <- z
    res[z>0.5] <- upperEmpTailDepFun(u)(z[z>0.5])
    res[z<=0.5] <- lowerEmpTailDepFun(u)(z[z<=0.5])
    return(res)
  }
}

##

lowerTailDepFun <- function(copula) {
  function(z) pCopula(cbind(z,z),copula)/z
}

upperTailDepFun <- function(copula) {
  function(z) (1-2*z+pCopula(cbind(z,z),copula))/(1-z)
}

tailDepFun <- function(copula) {
  function(z) {
    res <- z
    res[z>0.5] <- upperTailDepFun(copula)(z[z>0.5])
    res[z<=0.5] <- lowerTailDepFun(copula)(z[z<=0.5])
    return(res)
  }
}