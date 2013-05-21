########################################################
##                                                    ##
## functions based on sp preparing the use of copulas ##
##                                                    ##
########################################################

## neighbourhood constructor
############################

neighbourhood <- function(data, distances, sp, dataLocs=NULL, index, 
                          prediction, var) {
  sizeN <- ncol(distances)+1
  data <- as.data.frame(data)
  colnames(data) <- paste(paste("N", rep((0+as.numeric(prediction)):(sizeN-1), 
                                         each=length(var)), sep=""),
                          rep(var,(sizeN-prediction)),sep=".")
  if (anyDuplicated(rownames(data))>0)
    rownames <- 1:length(rownames)
  new("neighbourhood", data=data, distances=distances, locations=sp, 
      dataLocs=dataLocs, bbox=sp@bbox, proj4string=sp@proj4string, index=index,
      var=var, prediction=prediction)
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

selectFromNeighbourhood <- function(x, i) {
  newSp <- x@locations[i,]
  new("neighbourhood", data=x@data[i,,drop=F], 
      distances=x@distances[i,,drop=F], 
      locations=newSp, dataLocs=x@dataLocs, bbox=newSp@bbox, 
      proj4string=x@proj4string, index=x@index[i,,drop=F], var=x@var, 
      prediction=x@prediction)
}

setMethod("[[", "neighbourhood", selectFromNeighbourhood) 

## calculate neighbourhood from SpatialPointsDataFrame

# returns an neighbourhood object
##################################

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
  if (prediction)
    dataLocs <- spData
  else 
    dataLocs <- NULL
  return(neighbourhood(allData, allDists, locations, dataLocs,
                       allLocs, prediction, var))
}

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
setGeneric("calcBins", function(data, var, nbins=15, boundaries=NA, cutoff=NA,
                                cor.method="kendall", plot=TRUE, ...) {
                         standardGeneric("calcBins") 
                         })

## calculating the spatial bins
################################

calcSpBins <- function(data, var=names(data), nbins=15, boundaries=NA, 
                       cutoff=NA, cor.method="kendall", plot=TRUE) {

  if(is.na(cutoff)) {
    cutoff <- spDists(coordinates(t(data@bbox)))[1,2]/3
  }
  if(any(is.na(boundaries))) {
    boundaries <- ((1:nbins) * cutoff/nbins)
  }
  
  lags <- calcSpLagInd(data, boundaries)
    
  mDists <- sapply(lags, function(x) mean(x[,3]))
  np <- sapply(lags, function(x) length(x[,3]))
  lagData <- lapply(lags, function(x) as.matrix((cbind(data[x[,1],var]@data, data[x[,2],var]@data))))
  
  if(cor.method == "fasttau")
    lagCor <- sapply(lagData, function(x) VineCopula:::fasttau(x[,1], x[,2]))
  if(cor.method %in% c("kendall","spearman","perarson"))
    lagCor <- sapply(lagData, function(x) cor(x,method=cor.method)[1,2])
  if(cor.method == "normVariogram")  
    lagCor <- sapply(lagData, function(x) 1-cor(x,method="pearson")[1,2])
  if(cor.method == "variogram")  
    lagCor <- sapply(lagData, function(x) 0.5*mean((x[,1]-x[,2])^2,na.rm=T))
    
  if(plot) { 
    plot(mDists, lagCor, xlab="distance",ylab=paste("correlation [",cor.method,"]",sep=""), 
         ylim=1.05*c(-abs(min(lagCor)),max(lagCor)), xlim=c(0,max(mDists)))
    abline(h=c(-min(lagCor),0,min(lagCor)),col="grey")
  }
  
  res <- list(np=np, meanDists = mDists, lagCor=lagCor, lagData=lagData, lags=lags)
  attr(res,"cor.method") <- cor.method
  return(res)
}

setMethod(calcBins, signature("Spatial"), calcSpBins)

# calc bins from a (conditional) neighbourhood

calcNeighBins <- function(data, var=names(data), nbins=9, boundaries=NA, 
                          cutoff=NA, cor.method="kendall", plot=TRUE) {
  dists <- data@distances
  
  corFun <- switch(cor.method,
                   fasttau=function(x) VineCopula:::fasttau(x[,1],x[,2]),
                   function(x) cor(x,method=cor.method)[1,2])
  
  if (any(is.na(boundaries))) 
    boundaries <- quantile(as.vector(dists), probs=c(1:nbins/nbins))
  if(!is.na(cutoff)) {
    boundaries <- boundaries[boundaries < cutoff]
    boundaries <- unique(c(0,boundaries,cutoff))
  } else {
    boundaries <- unique(c(0,boundaries))
  }
  
  np <- numeric(0)
  moa <- numeric(0)
  lagData <- NULL
  meanDists <- numeric(0)

  data <- as.matrix(data@data)
  
  for ( i in 1:nbins) {
    bools <- (dists <= boundaries[i+1] & dists > boundaries[i])
    
    pairs <- NULL
    for(col in 1:(dim(bools)[2])) {
      pairs <- rbind(pairs, data[bools[,col],c(1,1+col)])
    }
    
    lagData <- append(lagData, list(pairs))
    moa <- c(moa, corFun(pairs))
    meanDists <- c(meanDists, mean(dists[bools]))
    np <- c(np, sum(bools))
  }
  
  if(plot) { 
    plot(meanDists, moa, xlab="distance", ylab=paste("correlation [",cor.method,"]",sep=""), 
         ylim=1.05*c(-abs(min(moa)),max(moa)), xlim=c(0,max(meanDists)))
    abline(h=c(-min(moa),0,min(moa)),col="grey")
  }
  
  res <- list(np=np, meanDists = meanDists, lagCor=moa, lagData=lagData)
  attr(res,"cor.method") <- switch(cor.method, fasttau="kendall", cor.method)
  return(res)
}
  
setMethod(calcBins, signature="neighbourhood", calcNeighBins)

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
