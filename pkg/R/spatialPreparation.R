##################################################################
##                                                              ##
## dedicated functions based on sp preparing the use of copulas ##
##                                                              ##
##################################################################

## neighbourhood
# constructor
# data = "matrix"	a matrix or array providing the data
# distances="matrix"	a matrix providing the distances
# sp="SpatialPoints"	SpatialPoints object providing the coordinates
# index="matrix"	linking the obs. in data to the coordinates

neighbourhood <- function(data, distances, sp, index, prediction, var){
  sizeN <- ncol(distances)+1
  data <- as.data.frame(data)
  colnames(data) <- paste(paste("N", rep((0+prediction):(sizeN-1), each=length(var)), sep=""),
                          rep(var,(sizeN-prediction)),sep=".")
  new("neighbourhood", data=data, distances=distances, locations=sp, 
      bbox=sp@bbox, proj4string=sp@proj4string, index=index, var=var, 
      prediction=prediction)
}

## show
showNeighbourhood <- function(object){
  cat("A set of neighbourhoods consisting of", ncol(object@distances)+1, "locations each \n")
  cat("with",nrow(object@data),"rows of observations for:\n")
  cat(object@var,"\n")
}

setMethod(show,signature("neighbourhood"),showNeighbourhood)

## names (from sp)

setMethod(names, signature("neighbourhood"), namesNeighbourhood <- function(x) x@var)

## spplot ##
spplotNeighbourhood <- function(obj, zcol=names(obj), ..., column=0) {
  stopifnot(all(column<ncol(obj@data)))
  pattern <- paste(paste("N", column, ".", sep=""), zcol, sep="")
  spdf <- SpatialPointsDataFrame(coords=obj@locations, 
                                 data=obj@data[,pattern,drop=FALSE], 
                                 proj4string=obj@proj4string, bbox=obj@bbox)
  spplot(spdf, ...)
}

setMethod(spplot, signature("neighbourhood"), spplotNeighbourhood)

## calculate neighbourhood from SpatialPointsDataFrame

# returns an neighbourhood object
# spData	    spatialPointsDataFrame
# locations   the prediction locations, for fitting, locations=spData
# var 		    one or multiple variable names, all is the default
# size		    the size of the neighbourhood, note that for prediction the size 
#             is one less than for the copula estimation (default of 5)
# min.dist    the minimum distance between neighbours, must be positive

getNeighbours <- function(spData, locations, var=names(spData)[1], size=5, 
                          prediction=FALSE, min.dist=0.01) {
  
  stopifnot((!prediction && missing(locations)) || (prediction && !missing(locations)))
  stopifnot(min.dist>0 || prediction)
  
  if(missing(locations) && !prediction)
    locations=spData
  
  stopifnot(is(locations,"Spatial"))
  
  nLocs <- length(locations)
  
  if(any(is.na(match(var,names(spData)))))
    stop("At least one of the variables is unkown or is not part of the data.")

  size <- min(size,length(spData)+prediction)

  allDists <- NULL
  allLocs <- NULL
  allData <- NULL
  
  for(i in 1:length(locations)) { # i <- 1
    tempDists <- spDistsN1(spData,locations[i,])
    tempDists[tempDists < min.dist] <- Inf
    spLocs <- order(tempDists)[1:(size-1)]
    allLocs <- rbind(allLocs, spLocs)
    allDists <- rbind(allDists, tempDists[spLocs])
    
    if(!prediction)
      spLocs <- c(i,spLocs)
    allData <- rbind(allData, as.vector(spData[spLocs, var, drop=F]@data[[1]]))
  }
  
  return(neighbourhood(allData, allDists, locations, 
                       allLocs, prediction, var))
}

# data(meuse)
# coordinates(meuse) <- ~x+y
# meuseNeigh <- getNeighbours(meuse,var="zinc",size=5)
# str(meuseNeigh)
# 
# meuseNeigh <- addAttrToGeom(meuseNeigh@locations,data.frame(rnd=runif(155)),match.ID=F)
# str(meuseNeigh)

#############
## BINNING ##
#############

# calculates lag indicies for a Spatial object and stores the respective separating distances
# 
# boundaries  -> are the right-side limits in the dimenssion as provided by spDists
# data --------> a spatial object that can be handled by spDists()        
calcSpLagInd <- function(data, boundaries) {
  lags <- vector("list",length(boundaries))
  
  dists <- spDists(data)
  nlocs <- length(data)
  
  for (i in 1:(nlocs-1)) {
    for (j in (i+1):nlocs) {
      d <- dists[i,j]
      for ( k in 1:length(boundaries)) {
        if (d < boundaries[k]) {
          lags[[k]] <- rbind(lags[[k]],c(i,j,d))
          break()
        }
      }
    }
  }
  return(lags)
}

# the generic calcBins, calculates bins for spatial and spatio-temporal data

setGeneric("calcBins", function(data, var, nbins=15, boundaries=NA, cutoff=NA, cor.method="kendall", plot=T, ...) standardGeneric("calcBins") )

# calculating the spatial bins
# 
# data denotes the spatial data object
# var denotes the only variable name used
# cor.method is passed on to cor() (default="kendall")
# if plot=TRUE (default), the correlation measures are plotted agaisnt the mean lag separation distance
# 
calcSpBins <- function(data, var=names(data), nbins=15, boundaries=NA, cutoff=NA, cor.method="kendall", plot=TRUE) {

  if(is.na(boundaries)) {
    diagonal <- spDists(coordinates(t(data@bbox)))[1,2]
    boundaries <- ((1:nbins) * min(cutoff,diagonal/3,na.rm=T) / nbins)
  }
  
  lags <- calcSpLagInd(data, boundaries)
    
  mDists <- sapply(lags,function(x) mean(x[,3]))
  lagData <- lapply(lags, function(x) as.matrix((cbind(data[x[,1],var]@data, data[x[,2],var]@data))))
  
  if(cor.method == "fasttau")
    lagCor <- sapply(lagData, function(x) VineCopula:::fasttau(x[,1], x[,2]))
  else 
    lagCor <- sapply(lagData, function(x) cor(x,method=cor.method)[1,2])
  
  if(plot) { 
    plot(mDists, lagCor, xlab="distance",ylab=paste("correlation [",cor.method,"]",sep=""), 
         ylim=1.05*c(-abs(min(lagCor)),max(lagCor)), xlim=c(0,max(mDists)))
    abline(h=c(-min(lagCor),0,min(lagCor)),col="grey")
  }
  
  res <- list(meanDists = mDists, lagCor=lagCor, lagData=lagData, lags=lags)
  attr(res,"cor.method") <- cor.method
  return(res)
}

setMethod(calcBins, signature("Spatial"), calcSpBins)

# instances: number  -> number of randomly choosen temporal intances
#            NA      -> all observations
#            other   -> temporal indexing as in spacetime/xts, the parameter t.lags is set to 0 in this case.
# t.lags:    numeric -> temporal shifts between obs
calcStBins <- function(data, var, nbins=15, boundaries=NA, cutoff=NA, instances=10, t.lags=c(0), cor.method="kendall", plot=TRUE) {

  if(is.na(cutoff)) cutoff <- spDists(coordinates(t(data@sp@bbox)))[1,2]/3
  if(is.na(boundaries)) boundaries <- ((1:nbins) * cutoff / nbins)
  if(is.na(instances)) instances=length(data@time)
  
  spIndices <- calcSpLagInd(data@sp, boundaries)
    
  mDists <- sapply(spIndices,function(x) mean(x[,3]))
  
  lengthTime <- length(data@time)
  if (!is.numeric(instances) | !length(instances)==1) {
    tempIndices <- cbind(instances, instances)
  } 
  else {
    tempIndices <- NULL
    for (t.lag in rev(t.lags)) {
      smplInd <- sample(x=max(1,1-t.lag):min(lengthTime,lengthTime-t.lag), size=min(instances,lengthTime-max(abs(t.lags))))
      tempIndices <- cbind(smplInd+t.lag, tempIndices)
      tempIndices <- cbind(tempIndices[,1]-t.lag, tempIndices)
    }
  }
    
  retrieveData <- function(spIndex, tempIndices) {
    binnedData <- NULL
    for (i in 1:(ncol(tempIndices)/2)) {
      binnedData <- cbind(binnedData, 
                          as.matrix((cbind(data[spIndex[,1], tempIndices[,2*i-1], var]@data, 
                                           data[spIndex[,2], tempIndices[,2*i], var]@data))))
    }
    return(binnedData)
  }
  
  lagData <- lapply(spIndices, retrieveData, tempIndices=tempIndices)
  
  calcStats <- function(binnedData) {
    cors <- NULL
    for(i in 1:(ncol(binnedData)/2)) {
      cors <- c(cors, cor(binnedData[,2*i-1], binnedData[,2*i], method=cor.method, use="pairwise.complete.obs"))
    }
    return(cors)
  }
  
  calcTau <- function(binnedData) {
    cors <- NULL
    for(i in 1:(ncol(binnedData)/2)) {
      cors <- c(cors, VineCopula:::fasttau(binnedData[,2*i-1], binnedData[,2*i]))
    }
    return(cors)
  }
  
  calcCor <- switch(cor.method, fasttau=calcTau, calcStats)
  
  lagCor <- sapply(lagData, calcCor)
  
  if(plot) { 
    plot(mDists, as.matrix(lagCor)[,1], xlab="distance",ylab=paste("correlation [",cor.method,"]",sep=""), 
         ylim=1.05*c(-abs(min(lagCor)),max(lagCor)), xlim=c(0,max(mDists)))
    abline(h=c(-min(lagCor),0,min(lagCor)),col="grey")
  }
  
  res <- list(meanDists = mDists, lagCor=lagCor, lagData=lagData, lags=list(sp=spIndices, time=tempIndices))
  attr(res,"cor.method") <- cor.method
  return(res)
}

setMethod(calcBins, signature("STFDF"), calcStBins)
