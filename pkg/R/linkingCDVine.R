#####################################################
## generic wrapper functions to the CDVine package ##
#####################################################


# density from BiCopPDF
linkCDVine.PDF <- function (u, copula, log=FALSE) {
  param <- copula@parameters
  if(length(param)==1) param <- c(param,0)
  n <- nrow(u)
  fam <- copula@family

  coplik = .C("LL_mod_seperate", as.integer(fam), as.integer(n), as.double(u[,1]), 
              as.double(u[,2]), as.double(param[1]), as.double(param[2]), 
              as.double(rep(0, n)), PACKAGE = "CDVine")[[7]]
  if(log) return(coplik)
  else return(exp(coplik))
}

# cdf from BiCopCDF

# for "standard" copulas: family %in% c(3:10)
linkCDVine.CDF <- function (u, copula) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  n <- nrow(u)
  fam <- copula@family
  
  res <- .C("archCDF", as.double(u[,1]), as.double(u[,2]), as.integer(n), as.double(param),
            as.integer(fam), as.double(rep(0, n)), PACKAGE = "CDVine")[[6]]
  return(res)
}

# for survival copulas: family %in% c(13, 14, 16:20)
linkCDVine.surCDF <- function (u, copula) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  u1 <- u[,1]
  u2 <- u[,2]
  n <- nrow(u)
  fam <- copula@family

  res <-  u1 + u2 - 1 + .C("archCDF", as.double(1 - u1), as.double(1 - u2), as.integer(n),
                           as.double(param), as.integer(fam - 10), as.double(rep(0, n)),
                           PACKAGE = "CDVine")[[6]]
  return(res)
}

# for 90 deg rotated copulas: family %in% c(23, 24, 26:30)
linkCDVine.r90CDF <- function (u, copula) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  u1 <- u[,1]
  u2 <- u[,2]
  n <- nrow(u)
  fam <- copula@family
  
  res <-  u2 - .C("archCDF", as.double(1 - u1), as.double(u2), as.integer(n), 
                  as.double(-param), as.integer(fam - 20), as.double(rep(0, n)), 
                  PACKAGE = "CDVine")[[6]]
  return(res)
}

# for 270 deg rotated copulas: family %in% c(33, 34, 36:40)
linkCDVine.r270CDF <- function (u, copula) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  u1 <- u[,1]
  u2 <- u[,2]
  n <- nrow(u)
  fam <- copula@family
  
  res <- u1 - .C("archCDF", as.double(u1), as.double(1 - u2), as.integer(n), 
                 as.double(-param), as.integer(fam - 30), as.double(rep(0, n)), 
                 PACKAGE = "CDVine")[[6]]
  return(res)
}

## derivtives/h-function  from BiCopHfunc
# ddu
linkCDVine.ddu <- function (u, copula) {
  param <- copula@parameters
  u <- matrix(u, ncol = 2)
  n <- nrow(u)
  fam <- copula@family
  
  res <- .C("Hfunc1", as.integer(fam), as.integer(n), as.double(u[,2]), as.double(u[,1]), 
            as.double(param[1]), as.double(param[2]), as.double(rep(0, n)), 
            PACKAGE = "CDVine")[[7]]
  return(res)
}

# ddv
linkCDVine.ddv <- function (u, copula) {
  param <- copula@parameters
  u <- matrix(u, ncol = 2)
  n <- nrow(u)
  fam <- copula@family
  
  res <- .C("Hfunc2", as.integer(fam), as.integer(n), as.double(u[,1]), as.double(u[,2]), 
            as.double(param[1]), as.double(param[2]), as.double(rep(0, n)), 
            PACKAGE = "CDVine")[[7]]
  return(res)
}


## random numbers from CDVineSim
linkCDVine.r <- function (n, copula){
  param <- copula@parameters
  fam <- copula@family
  if(is.na(param[2])) param <- c(param,0)
  
  tmp <- .C("pcc", as.integer(n), as.integer(2), as.integer(fam), as.integer(1), 
            as.double(param[1]), as.double(param[2]), as.double(rep(0, n * 2)), 
            PACKAGE = "CDVine")[[7]]
  return(matrix(tmp, ncol = 2))
}

## transform a fit from CDVine to a list of copula objects
castCDvine <- function(cdvEst) {
  copulas <- NULL
  for(i in 1: length(cdvEst$family)) {
    par1 <- cdvEst$par[i]
    par2 <- cdvEst$par2[i]
    cop <- switch(paste("fam",cdvEst$family[i],sep=""), fam0=indepCopula(dim=2), fam1=normalCopula(par1), fam2=tCopula(par1,df=par2), 
                  fam3=claytonCopula(par1),         fam4=gumbelCopula(par1),          fam5=frankCopula(par1),           fam6=joeBiCopula(par1), 
                  fam7=BB1Copula(c(par1,par2)),     fam8=BB6Copula(c(par1,par2)),     fam9=BB7Copula(c(par1,par2)),    fam10=BB8Copula(c(par1,par2)),
                  fam13=surClaytonCopula(par1),     fam14=surGumbelCopula(par1),      fam16=surJoeBiCopula(par1),
                  fam17=surBB1Copula(c(par1,par2)), fam18=surBB6Copula(c(par1,par2)), fam19=surBB7Copula(c(par1,par2)), fam20=surBB8Copula(c(par1,par2)),
                  fam23=r90ClaytonCopula(par1),     fam24=r90GumbelCopula(par1),      fam26=r90JoeBiCopula(par1),
                  fam27=r90BB1Copula(c(par1,par2)), fam28=r90BB6Copula(c(par1,par2)), fam29=r90BB7Copula(c(par1,par2)), fam30=r90BB8Copula(c(par1,par2)),
                  fam33=r270ClaytonCopula(par1),    fam34=r270GumbelCopula(par1),     fam36=r270JoeBiCopula(par1),
                  fam37=r270BB1Copula(c(par1,par2)),fam38=r270BB6Copula(c(par1,par2)),fam39=r270BB7Copula(c(par1,par2)),fam40=r270BB8Copula(c(par1,par2)))
    
    copulas <- append(copulas, cop)
  }
  if(length(copulas) ==1) copulas <- copulas[[1]]
  return(copulas)
}

## Kendall's tau
linkCDVine.tau <- function(copula) {
  par <- copula@parameters
  BiCopPar2Tau(copula@family, par[1], par[2])
}

## get parameter from Kendall's tau (only for one parameter families)
linkCDVine.iTau <- function(copula, tau) {
  BiCopTau2Par(copula@family, tau)
}

## tailIndex
linkCDVine.tailIndex <- function(copula) {
  par <- copula@parameters
  unlist(BiCopPar2TailDep(copula@family,par[1],par[2]))
}