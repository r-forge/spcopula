## ##################################################
## generic wrapper functions to the CDVine package ##
#####################################################


# density from BiCopPDF
linkCDVine.PDF <- function (copula, u) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  n <- nrow(u)
  u1 <- u[,1]
  u2 <- u[,2]
  fam <- copula@family

# workaround, check updates in CDVine using
#   family <- 27
# simData <- BiCopSim(500, family, -1.5, -1)
# fit <- BiCopSelect(simData[,1],simData[,2], familyset=family+c(-20,-10,0,10))
# sum(log(BiCopPDF(simData[,1],simData[,2],family=fit$family, fit$par, fit$par2))) # 190
# sum(log(BiCopPDF(simData[,1],simData[,2],family=fit$family+10, fit$par, fit$par2))) # 111
  if (fam %in% c(23,24,26:30,33,34,36:40)) {
    u1 <- u[,2]
    u2 <- u[,1]
  }

  coplik = .C("LL_mod_seperate", as.integer(fam), as.integer(n), as.double(u1), as.double(u2), as.double(param[1]), as.double(param[2]), as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
  return(exp(coplik))
}

# cdf from BiCopCDF

# for "standard" copulas: family %in% c(3:10)
linkCDVine.CDF <- function (copula, u) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  n <- nrow(u)
  fam <- copula@family
  
  res <- .C("archCDF", as.double(u[,1]), as.double(u[,2]), as.integer(n), as.double(param), as.integer(fam), as.double(rep(0, n)), PACKAGE = "CDVine")[[6]]
  return(res)
}

# for survival copulas: family %in% c(13, 14, 16:20)
linkCDVine.surCDF <- function (copula, u) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  u1 <- u[,1]
  u2 <- u[,2]
  n <- nrow(u)
  fam <- copula@family

  res <-  u1 + u2 - 1 + .C("archCDF", as.double(1 - u1), as.double(1 - u2), as.integer(n), as.double(param), as.integer(fam - 10), as.double(rep(0, n)), PACKAGE = "CDVine")[[6]]
  return(res)
}

# for 90° rotated copulas: family %in% c(23, 24, 26:30)
linkCDVine.r90CDF <- function (copula, u) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  u1 <- u[,1]
  u2 <- u[,2]
  n <- nrow(u)
  fam <- copula@family
  
  res <-  u2 - .C("archCDF", as.double(1 - u1), as.double(u2), as.integer(n), as.double(-param), as.integer(fam - 20), as.double(rep(0, n)), PACKAGE = "CDVine")[[6]]
  return(res)
}

# for 270° rotated copulas: family %in% c(33, 34, 36:40)
linkCDVine.r270CDF <- function (copula, u) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  u1 <- u[,1]
  u2 <- u[,2]
  n <- nrow(u)
  fam <- copula@family
  
  res <- u1 - .C("archCDF", as.double(u1), as.double(1 - u2), as.integer(n), as.double(-param), as.integer(fam - 30), as.double(rep(0, n)), PACKAGE = "CDVine")[[6]]
  return(res)
}

## derivtives/h-function  from BiCopHfunc
# ddu
linkCDVine.ddu <- function (copula, pair) {
  param <- copula@parameters
  u <- matrix(pair, ncol = 2)
  n <- nrow(pair)
  fam <- copula@family
  
  res <- .C("Hfunc1", as.integer(fam), as.integer(n), as.double(u[,2]), as.double(u[,1]), as.double(param[1]), as.double(param[2]), as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
  return(res)
}

# ddv
linkCDVine.ddv <- function (copula, pair) {
  param <- copula@parameters
  u <- matrix(pair, ncol = 2)
  n <- nrow(pair)
  fam <- copula@family
  
  res <- .C("Hfunc2", as.integer(fam), as.integer(n), as.double(u[,1]), as.double(u[,2]), as.double(param[1]), as.double(param[2]), as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
  return(res)
}


## random numbers from CDVineSim
linkCDVine.r <- function(copula, n){
  param <- copula@parameters
  fam <- copula@family
  if(is.na(param[2])) param <- c(param,0)
  
  tmp <- .C("pcc", as.integer(n), as.integer(2), as.integer(fam), as.integer(1), as.double(param[1]), as.double(param[2]), as.double(rep(0, n * 2)), PACKAGE = "CDVine")[[7]]
  return(matrix(tmp, ncol = 2))
}


