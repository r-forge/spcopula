# constructor
# dimension = "integer"     set to 2
# parameters = "numeric"    set of parameters
# param.names = "character" appropriate names
# param.lowbnd = "numeric"  appropriate lower bounds
# param.upbnd = "numeric"   appropriate upper bounds
# fullname = "character"     messgae printed with "show"
# components="list"         list of copulas (will be automatically supplemented 
#			      by the independent copula)
# distances="numeric"       the linking distances + the range (will be assigned
#			      to the independent copula)
# unit="character"          measurement unit of distance
# depFun="function"         a optional spatial dependence function providing 
#                             Kendalls tau or Spearman's rho to calib* or exact 
#                             parameters

stCopula <- function(components, distances, t.lags, stDepFun, unit="m", t.res="day") {
  spCopList <- list()
  
  if(!missing(stDepFun)) {
    getSpCop <- function(comp,dist,time) spCopula(comp, dist,
                                                  spDepFun=function(h) stDepFun(h,time), unit)
    for(i in 1:length(t.lags)){
      spCopList <- append(spCopList, getSpCop(components[[i]], distances[[i]], i))
    }
  } else {
    for(i in 1:length(t.lags)){
      spCopList <- append(spCopList, spCopula(components[[i]], distances[[i]], unit=unit))
    }
  }
  
  param       <- unlist(lapply(spCopList, function(x) x@parameters))
  param.names <- unlist(lapply(spCopList, function(x) x@param.names))
  param.low   <- unlist(lapply(spCopList, function(x) x@param.lowbnd))
  param.up    <- unlist(lapply(spCopList, function(x) x@param.upbnd))
  
  new("stCopula", dimension=as.integer(2), parameters=param, param.names=param.names,
      param.lowbnd=param.low, param.upbnd=param.up,
      fullname="Spatio-Temporal Copula: distance and time dependent convex combination of bivariate copulas",
      spCopList=spCopList, t.lags=t.lags, t.res=t.res)
}

## show method
showStCopula <- function(object) {
  cat(object@fullname, "\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Copulas:\n")
  for (i in 1:length(object@spCopList)) {
    cmpCop <- object@spCopList[[i]]
    cat("  ", cmpCop@fullname, "at", object@t.lags[i], 
      paste("[",object@t.res,"]",sep=""), "\n")
    show(cmpCop)
  }
}

setMethod("show", signature("stCopula"), showStCopula)

## spatial copula jcdf ##

# TODO: add again block support to the spatio-temporal copula
# u 
#   list containing two column matrix providing the transformed pairs,  their respective 
#   separation distances and time steps
pStCopula <- function (u, copula) {
  if (!is.list(u) || !length(u)>=3) stop("Point pairs need to be provided with their separating spatial and temproal distances as a list.")
  
  if(!is.matrix(u[[1]])) u[[1]] <- matrix(u[[1]],ncol=2)
  n <- nrow(u[[1]])
  h <- u[[2]]
  t.dist <- u[[3]]
  
  if(length(u)==4) {
    block <- u[[4]]
    if (n%%block != 0) stop("The block size is not a multiple of the data length:",n)
  } else block <- 1

  if(any(is.na(match(t.dist,copula@t.lags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  if(length(h)>1 & length(h)!=n) 
    stop("The spatial distance vector must either be of same length as rows in the data pairs or a single value.")
  if(length(t.dist)>1 & length(t.dist)!=n) 
    stop("The temporal distances vector must either be of same length as rows in the data pairs or a single value.")

  if (length(t.dist)==1) {
    res <- pSpCopula(copula@spCopList[[match(t.dist,copula@t.lags)]], 
                     list(u[[1]], h))
  } else {
    if(length(h)==1) h <- rep(h,n)
    res <- NULL
    for(i in 1:(n%/%block)) {
      res <- rbind(res, pSpCopula(copula@spCopList[[match(t.dist[i],copula@t.lags)]],
                                  list(u[[1]][((i-1)*block+1):(i*block),], h[i*block])))
    }
  }
  
  return(res)
}

setMethod("pCopula", signature("numeric","stCopula"), pStCopula)


## spatial Copula density ##

# u 
#   three column matrix providing the transformed pairs and their respective 
#   separation distances
dStCopula <- function (u, copula) {
  if (!is.list(u) || !length(u)>=3) stop("Point pairs need to be provided with their separating spatial and temproal distances as a list.")
  
  if(!is.matrix(u[[1]])) u[[1]] <- matrix(u[[1]],ncol=2)
  n <- nrow(u[[1]])
  h <- u[[2]]
  t.dist <- u[[3]]
  
  if(length(u)==4) {
    block <- u[[4]]
    if (n%%block != 0) stop("The block size is not a multiple of the data length:",n)
  } else block <- 1
  
  if(any(is.na(match(t.dist,copula@t.lags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  if(length(h)>1 & length(h)!=n) 
    stop("The spatial distance vector must either be of same length as rows in the data pairs or a single value.")
  if(length(t.dist)>1 & length(t.dist)!=n) 
    stop("The temporal distances vector must either be of same length as rows in the data pairs or a single value.")
  
  if (length(t.dist)==1) {
    res <- dSpCopula(copula@spCopList[[match(t.dist,copula@t.lags)]], 
                     list(u[[1]], h))
  } else {
    if(length(h)==1) h <- rep(h,n)
    res <- NULL
    for(i in 1:(n%/%block)) {
      res <- rbind(res, dSpCopula(copula@spCopList[[match(t.dist[i],copula@t.lags)]],
                                  list(u[[1]][((i-1)*block+1):(i*block),], h[i*block])))
    }
  }
  
  return(res)
}

setMethod("dCopula", signature("list","stCopula"), dStCopula)


## partial derivatives ##

## dduSpCopula
###############

dduStCopula <- function (u, copula) {
  if (!is.list(u) || !length(u)>=3) stop("Point pairs need to be provided with their separating spatial and temproal distances as a list.")
  
  if(!is.matrix(u[[1]])) u[[1]] <- matrix(u[[1]],ncol=2)
  n <- nrow(u[[1]])
  h <- u[[2]]
  t.dist <- u[[3]]
  
  if(length(u)==4) {
    t.block <- u[[4]]
    if (n%%t.block != 0) stop("The block size is not a multiple of the data length:",n)
  } else t.block <- 1
  
  if(any(is.na(match(t.dist,copula@t.lags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  if(length(h)>1 & length(h)!=n) 
    stop("The spatial distance vector must either be of same length as rows in the data pairs or a single value.")
  if(length(t.dist)>1 & length(t.dist)!=n) 
    stop("The temporal distances vector must either be of same length as rows in the data pairs or a single value.")
  
  if (length(t.dist)==1) {
    res <- dduSpCopula(copula@spCopList[[match(t.dist,copula@t.lags)]],
                       list(u[[1]], h, block=t.block))
  } else {
    if(length(h)==1) h <- rep(h,n)
    res <- NULL
    for(i in 1:(n%/%t.block)) {
      cop <- copula@spCopList[[match(t.dist[i*t.block],copula@t.lags)]]
      tmpPair <- u[[1]][((i-1)*t.block+1):(i*t.block),]
      res <- rbind(res, dduSpCopula(cop, list(tmpPair,h[i*t.block])))
    }
  }
  
  return(res)
}

setMethod("dduCopula", signature("list","stCopula"), dduStCopula)


## ddvSpCopula
###############
ddvStCopula <- function (u, copula) {
  if (!is.list(u) || !length(u)>=3) stop("Point pairs need to be provided with their separating spatial and temproal distances as a list.")
  
  if(!is.matrix(u[[1]])) u[[1]] <- matrix(u[[1]],ncol=2)
  n <- nrow(u[[1]])
  h <- u[[2]]
  t.dist <- u[[3]]
  
  if(length(u)==4) {
    t.block <- u[[4]]
    if (n%%t.block != 0) stop("The block size is not a multiple of the data length:",n)
  } else t.block <- 1
  
  if(any(is.na(match(t.dist,copula@t.lags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  if(length(h)>1 & length(h)!=n) 
    stop("The spatial distance vector must either be of same length as rows in the data pairs or a single value.")
  if(length(t.dist)>1 & length(t.dist)!=n) 
    stop("The temporal distances vector must either be of same length as rows in the data pairs or a single value.")
  
  if (length(t.dist)==1) {
    res <- ddvSpCopula(copula@spCopList[[match(t.dist,copula@t.lags)]],
                       list(u[[1]], h, block=t.block))
  } else {
    if(length(h)==1) h <- rep(h,n)
    res <- NULL
    for(i in 1:(n%/%t.block)) {
      cop <- copula@spCopList[[match(t.dist[i*t.block],copula@t.lags)]]
      tmpPair <- u[[1]][((i-1)*t.block+1):(i*t.block),]
      res <- rbind(res, ddvSpCopula(cop, list(tmpPair,h[i*t.block])))
    }
  }
  
  return(res)
}

setMethod("ddvCopula", signature("list","stCopula"), ddvStCopula)

# #############
# ##         ##
# ## FITTING ##
# ##         ##
# #############
# 
# # two models: 
# # 1) Kendall's tau driven:
# #    fit curve through emp. Kendall's tau values, identify validity ranges for
# #    copula families deriving parameters from the fit, fade from one family to 
# #    another at borders
# # 2) convex-linear combination of copulas: 
# #    fit one per lag, fade from one to another
# 
# # towards the first model:
# 
# # INPUT: the stBinning
# # steps
# # a) fit a curve
# # b) estimate bivariate copulas per lag (limited to those with some 1-1-relation 
# #    to Kendall's tau')
# # INTERMEDIATE RESULT
# # c) select best fits based on ... e.g. log-likelihood, visual inspection
# # d) compose bivariate copulas to one spatial copula
# # OUTPUT: a spatial copula parametrixued by distance through Kendall's tau
# 
# 
# # towards b)
# # bins     -> typically output from calcBins
# # calcTau  -> a function on distance providing Kendall's tau estimates, 
# # families -> a vector of dummy copula objects of each family to be considered
# #             DEFAULT: c(normal, t_df=4, clayton, frank, gumbel
# loglikByCopulasLags <- function(bins, calcTau, families=c(normalCopula(0), tCopula(0,dispstr="un"),
#                                                           claytonCopula(0), frankCopula(1), gumbelCopula(1))) {
#   loglik <- NULL
#   for (cop in families) {
#     print(cop)
#     tmploglik <- NULL
#     for(i in 1:length(bins$meanDists)) {
#       cop@parameters[1] <- iTau(cop,tau=calcTau(bins$meanDists[i]))
#       tmploglik <- c(tmploglik, sum(log(dCopula(bins$lagData[[i]],cop))))
#     }
#     loglik <- cbind(loglik, tmploglik)
#   }
# 
#   colnames(loglik) <- sapply(families, function(x) class(x)[1])
# 
#   return(loglik)
# }
# 
# # towards d)
# composeSpCop <- function(bestFit, families, bins, calcCor) {
#   nfits <- length(bestFit)
#   gaps <- which(diff(bestFit)!=0)
# 
#   if(missing(calcCor)) noCor <- nfits
#   else noCor <- min(which(calcCor(bins$meanDists)<=0), nfits)
#   
#   breaks <- sort(c(gaps, gaps+1, noCor))
#   breaks <- breaks[breaks<noCor]
#   
#   cops <- as.list(families[bestFit[breaks]])
#   
#   breaks <- unique(c(breaks, min(nfits,rev(breaks)[1]+1)))
#   distances <- bins$meanDists[breaks]
#   
#   if(length(breaks)==length(cops)) {
#     warning("Lags do not cover the entire range with positive correlation. ", 
#              "An artifical one fading towards the indepCopula has been added.")
#     distances <- c(distances, rev(distances)[1]+diff(bins$meanDists[nfits+c(-1,0)]))
#   }
# 
#   if(missing(calcCor)) return(spCopula(components=cops, distances=distances, unit="m"))
#   else return(spCopula(components=cops, distances=distances, unit="m", spDepFun=calcCor))
# }
# 
# # in once
# 
# # bins   -> typically output from calcBins
# # cutoff -> maximal distance that should be considered for fitting
# # families -> a vector of dummy copula objects of each family to be considered
# #             DEFAULT: c(normal, t_df=4, clayton, frank, gumbel
# # ...
# # type   -> the type of curve (by now only polynominals are supported)
# # degree -> the degree of the polynominal
# # bounds -> the bounds of the correlation function (typically c(0,1))
# # method -> the measure of association, either "kendall" or "spearman"
# fitSpCopula <- function(bins, cutoff=NA, families=c(normalCopula(0), tCopula(0,dispstr="un"),
#                                                     claytonCopula(0), frankCopula(1), gumbelCopula(1)), ...) {
#   calcTau <- fitCorFun(bins, cutoff=cutoff, ...)
#   loglik <- loglikByCopulasLags(bins, calcTau=calcTau, families=families)
#   
#   bestFit <- apply(apply(loglik, 1, rank),2,function(x) which(x==length(families)))
#   
#   return(composeSpCop(bestFit, families, bins, calcTau))
# }
