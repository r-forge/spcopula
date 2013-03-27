## utilities


# ranks are automatically removed and NAs are by default randomly distributed
rankTransform <- function(u,v=NULL, ties.method="average") {
  if(!(is.matrix(u) | is.data.frame(u))) {
    if (is.null(v))
      stop("u must either be a matrix with at least 2 columns or u and v must be given.")
    else
      u <- cbind(u,v)
  }

  bool <- apply(u,1,function(row) !any(is.na(row)))
  res <- apply(u[bool,],2,rank,ties.method=ties.method)/(sum(bool)+1)
  if(is.data.frame(u))
    return(as.data.frame(res))
  return(res)
}

##
dependencePlot <- function(var=NULL, smpl, bandwidth=0.075, 
                           main="Stength of dependece", 
                           transformation=function (x) x, ...) {
  if(is.null(var)) {
    if (ncol(smpl)>2) {
      smpl <- smpl[,1:2]
    }
  } else {
    smpl <- smpl[,var]
  }
  
  smoothScatter(smpl,bandwidth=bandwidth, asp=1, xlim=c(0,1), ylim=c(0,1), 
                nrpoints=0, main=main,
                transformation=transformation, ...)
}

##
unitScatter <- function(var=NULL, smpl, ...) {
  
  if(is.null(var)) {
    if (ncol(smpl)>2) {
      smpl <- smpl[,1:2]
    }
  } else {
    smpl <- smpl[,var]
  }

  for(variable in var){
    if( min(smpl[,variable])<0 | max(smpl[,variable])>1) {
      smpl[,variable] <- rank(smpl[,variable])/(length(smpl[,variable])+1)
      warning("The variable ",variable," seems to exceed [0,1] and has been transformed using the rank order transformation.")
    }
  }

  plot(smpl, asp=1, xlim=c(0,1), ylim=c(0,1), ...)
}

univScatter <- function(formula=NULL, smpl) {
  .Deprecated("unitScatter")
  unitScatter(formula, smpl)
}

copulaFromFamilyIndex <- function(family, par, par2=0) {
  constr <- switch(family+1, function(par) indepCopula(), 
                   function(par) normalCopula(par[1]), function(par) tCopula(par[1],df=par[2]),
                   function(par) claytonCopula(par[1]), function(par) gumbelCopula(par[1]),
                   function(par) frankCopula(par[1]), function(par) joeBiCopula(par[1]),
                   BB1Copula, BB6Copula, BB7Copula, BB8Copula, 
                   NULL, NULL,
                   function(par) surClaytonCopula(par[1]), function(par) surGumbelCopula(par[1]),
                   NULL, function(par) surJoeBiCopula(par[1]),
                   surBB1Copula, surBB6Copula, surBB7Copula, surBB8Copula, 
                   NULL, NULL,
                   function(par) r90ClaytonCopula(par[1]), function(par) r90GumbelCopula(par[1]),
                   NULL, function(par) r90JoeBiCopula(par[1]),
                   r90BB1Copula, r90BB6Copula, r90BB7Copula, r90BB8Copula, 
                   NULL, NULL,
                   function(par) r270ClaytonCopula(par[1]), function(par) r270GumbelCopula(par[1]),
                   NULL, function(par) r270JoeBiCopula(par[1]),
                   r270BB1Copula, r270BB6Copula, r270BB7Copula, r270BB8Copula)
  constr(c(par,par2))
}