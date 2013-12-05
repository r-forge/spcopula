copulaFromFamilyIndex <- function(family, par, par2=0) {
  constr <- switch(paste("fam",family,sep=""),
                   fam0 = function(par) indepCopula(), 
                   fam1 = function(par) normalCopula(par[1]),
                   fam2 = function(par) tCopula(par[1],df=par[2]),
                   fam3 = function(par) claytonCopula(par[1]),
                   fam4 = function(par) gumbelCopula(par[1]),
                   fam5 = function(par) frankCopula(par[1]), 
                   fam6 = function(par) joeBiCopula(par[1]),
                   fam7 = BB1Copula,
                   fam8 = BB6Copula, 
                   fam9 = BB7Copula, 
                   fam10 = BB8Copula, 
                   fam13 = function(par) surClaytonCopula(par[1]), 
                   fam14 = function(par) surGumbelCopula(par[1]),
                   fam16 = function(par) surJoeBiCopula(par[1]),
                   fam17 = surBB1Copula, 
                   fam18 = surBB6Copula, 
                   fam19 = surBB7Copula, 
                   fam20 = surBB8Copula, 
                   fam23 = function(par) r90ClaytonCopula(par[1]),
                   fam24 = function(par) r90GumbelCopula(par[1]),
                   fam26 = function(par) r90JoeBiCopula(par[1]),
                   fam27 = r90BB1Copula,
                   fam28 = r90BB6Copula,
                   fam29 = r90BB7Copula, 
                   fam30 = r90BB8Copula, 
                   fam33 = function(par) r270ClaytonCopula(par[1]),
                   fam34 = function(par) r270GumbelCopula(par[1]),
                   fam36 = function(par) r270JoeBiCopula(par[1]),
                   fam37 = r270BB1Copula, 
                   fam38 = r270BB6Copula, 
                   fam39 = r270BB7Copula, 
                   fam40 = r270BB8Copula, 
                   fam104 = tawnT1Copula,
                   fam114 = surTawnT1Copula,
                   fam124 = r90TawnT1Copula,
                   fam134 = r270TawnT1Copula,
                   fam204 = tawnT2Copula,
                   fam214 = surTawnT2Copula,
                   fam224 = r90TawnT2Copula,
                   fam234 = r270TawnT2Copula)
  constr(c(par,par2))
}

#########################################################
## generic wrapper functions to the VineCopula package ##
#########################################################

# density from BiCopPDF
linkVineCop.PDF <- function (u, copula, log=FALSE) {
  param <- copula@parameters

  if(length(param)==1) 
    param <- c(param,0)
  n <- nrow(u)
  fam <- copula@family

  coplik = RLL_mod_separate(fam, n, u, param)[[7]]
#   coplik = .C("LL_mod_seperate", as.integer(fam), as.integer(n), as.double(u[,1]), 
#               as.double(u[,2]), as.double(param[1]), as.double(param[2]), 
#               as.double(rep(0, n)), PACKAGE = "VineCopula")[[7]]
  if(log) 
    return(coplik)
  else 
    return(exp(coplik))
}

# cdf from BiCopCDF

# for "standard" copulas: family %in% c(3:10)
linkVineCop.CDF <- function (u, copula) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  n <- nrow(u)
  fam <- copula@family
  
  res <- RarchCDF(fam, n, u, param)[[6]]
#   res <- .C("archCDF", as.double(u[,1]), as.double(u[,2]), as.integer(n), as.double(param),
#             as.integer(fam), as.double(rep(0, n)), PACKAGE = "VineCopula")[[6]]
  return(res)
}

# for survival copulas: family %in% c(13, 14, 16:20)
linkVineCop.surCDF <- function (u, copula) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  u1 <- u[,1]
  u2 <- u[,2]
  n <- nrow(u)
  fam <- copula@family

  res <- u1 + u2 - 1 + RarchCDF(fam-10, n, cbind(1-u1,1-u2), param)[[6]]
#   res <-  u1 + u2 - 1 + .C("archCDF", as.double(1 - u1), as.double(1 - u2), as.integer(n),
#                            as.double(param), as.integer(fam - 10), as.double(rep(0, n)),
#                            PACKAGE = "VineCopula")[[6]]
  return(res)
}

# for 90 deg rotated copulas: family %in% c(23, 24, 26:30)
linkVineCop.r90CDF <- function (u, copula) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  u1 <- u[,1]
  u2 <- u[,2]
  n <- nrow(u)
  fam <- copula@family
  
  res <- u2 - RarchCDF(fam - 20, n, cbind(1-u1,u2), -param)[[6]]
#   u2 - .C("archCDF", as.double(1 - u1), as.double(u2), as.integer(n), 
#                   as.double(-param), as.integer(fam - 20), as.double(rep(0, n)), 
#                   PACKAGE = "VineCopula")[[6]]
  return(res)
}

# for 270 deg rotated copulas: family %in% c(33, 34, 36:40)
linkVineCop.r270CDF <- function (u, copula) {
  param <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol = 2)
  u1 <- u[,1]
  u2 <- u[,2]
  n <- nrow(u)
  fam <- copula@family
  
  res <- u1 - RarchCDF(fam-30, n, cbind(u1,1-u2), -param)[[6]]
#     u1 - .C("archCDF", as.double(u1), as.double(1 - u2), as.integer(n), 
#                  as.double(-param), as.integer(fam - 30), as.double(rep(0, n)), 
#                  PACKAGE = "VineCopula")[[6]]
  return(res)
}

## derivtives/h-function  from BiCopHfunc
# ddu
linkVineCop.ddu <- function (u, copula) {
  param <- copula@parameters
  u <- matrix(u, ncol = 2)
  n <- nrow(u)
  fam <- copula@family
  
  res <- RHfunc1(fam, n, u, param)[[7]]
#     .C("Hfunc1", as.integer(fam), as.integer(n), as.double(u[,2]), as.double(u[,1]), 
#             as.double(param[1]), as.double(param[2]), as.double(rep(0, n)), 
#             PACKAGE = "VineCopula")[[7]]
  return(res)
}

# ddv
linkVineCop.ddv <- function (u, copula) {
  param <- copula@parameters
  u <- matrix(u, ncol = 2)
  n <- nrow(u)
  fam <- copula@family
  
  res <- RHfunc2(fam, n, u, param)[[7]]
#     .C("Hfunc2", as.integer(fam), as.integer(n), as.double(u[,1]), as.double(u[,2]), 
#             as.double(param[1]), as.double(param[2]), as.double(rep(0, n)), 
#             PACKAGE = "VineCopula")[[7]]
  return(res)
}


## random numbers from VineCopulaSim
linkVineCop.r <- function (n, copula){
  param <- copula@parameters
  fam <- copula@family
  if(is.na(param[2])) param <- c(param,0)
  
  tmp <- Rpcc(fam, n, param)[[7]]
#     .C("pcc", as.integer(n), as.integer(2), as.integer(fam), as.integer(1), 
#             as.double(param[1]), as.double(param[2]), as.double(rep(0, n * 2)), 
#             PACKAGE = "VineCopula")[[7]]
  return(matrix(tmp, ncol = 2))
}


# ## transform a fit from VineCopula to a list of copula objects
# castVineCopula <- function(cdvEst) {
#   copulas <- NULL
#   for(i in 1: length(cdvEst$family)) {
#     par1 <- cdvEst$par[i]
#     par2 <- cdvEst$par2[i]
#     cop <- switch(paste("fam",cdvEst$family[i],sep=""), fam0=indepCopula(dim=2), fam1=normalCopula(par1), fam2=tCopula(par1,df=par2), 
#                   fam3=claytonCopula(par1),         fam4=gumbelCopula(par1),          fam5=frankCopula(par1),           fam6=joeBiCopula(par1), 
#                   fam7=BB1Copula(c(par1,par2)),     fam8=BB6Copula(c(par1,par2)),     fam9=BB7Copula(c(par1,par2)),    fam10=BB8Copula(c(par1,par2)),
#                   fam13=surClaytonCopula(par1),     fam14=surGumbelCopula(par1),      fam16=surJoeBiCopula(par1),
#                   fam17=surBB1Copula(c(par1,par2)), fam18=surBB6Copula(c(par1,par2)), fam19=surBB7Copula(c(par1,par2)), fam20=surBB8Copula(c(par1,par2)),
#                   fam23=r90ClaytonCopula(par1),     fam24=r90GumbelCopula(par1),      fam26=r90JoeBiCopula(par1),
#                   fam27=r90BB1Copula(c(par1,par2)), fam28=r90BB6Copula(c(par1,par2)), fam29=r90BB7Copula(c(par1,par2)), fam30=r90BB8Copula(c(par1,par2)),
#                   fam33=r270ClaytonCopula(par1),    fam34=r270GumbelCopula(par1),     fam36=r270JoeBiCopula(par1),
#                   fam37=r270BB1Copula(c(par1,par2)),fam38=r270BB6Copula(c(par1,par2)),fam39=r270BB7Copula(c(par1,par2)),fam40=r270BB8Copula(c(par1,par2)))
#     
#     copulas <- append(copulas, cop)
#   }
#   if(length(copulas) ==1) copulas <- copulas[[1]]
#   return(copulas)
# }

## Kendall's tau
linkVineCop.tau <- function(copula) {
  par <- copula@parameters
  BiCopPar2Tau(copula@family, par[1], par[2])
}

## get parameter from Kendall's tau (only for one parameter families)
linkVineCop.iTau <- function(copula, tau) {
  BiCopTau2Par(copula@family, tau)
}

## tailIndex
linkVineCop.tailIndex <- function(copula) {
  par <- copula@parameters
  unlist(BiCopPar2TailDep(copula@family,par[1],par[2]))
}