#################################################################################
##
##  R package spcopula by Benedikt Graeler Copyright (C) 2011
##
##  This file is part of the R package spcopula.
##
##  The R package spcopula is free software: you can redistribute it and/or 
##  modify it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  The R package spcopula is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with the R package spcopula. If not, see <http://www.gnu.org/licenses/>
##
#################################################################################

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

neighbourhood <- function(data, distances, sp, index){
  varNames <- names(data[[1]])
  sizeN <- ncol(distances)+1
  data <- as.data.frame(data)
  colnames(data) <- paste(paste("N",rep(0:(sizeN-1),each=length(varNames)),sep=""),rep(varNames,sizeN),sep=".")
  new("neighbourhood", data=data, distances=distances, coords=sp@coords, bbox=sp@bbox, proj4string=sp@proj4string, index=index, varNames=varNames)
}

## show
showNeighbourhood <- function(object){
  cat("A set of neighbourhoods consisting of", ncol(object@distances)+1, "locations each \n")
  cat("with",nrow(object@data),"rows of observations for:\n")
  cat(object@varNames,"\n")
}

setMethod(show,signature("neighbourhood"),showNeighbourhood)

## names (from sp)

setMethod(names, signature("neighbourhood"), namesNeighbourhood <- function(x) x@varNames)

## spplot ##

spplotNeighbourhood <- function(obj, zcol=names(obj), ..., column=0) {
  pattern <- paste(paste("N",column,".",sep=""),zcol,sep="")
  spdf <- SpatialPointsDataFrame(coords=obj@coords, data=obj@data[,pattern,drop=FALSE], proj4string=obj@proj4string, bbox=obj@bbox)
  spplot(spdf, ...)
}

setMethod(spplot, signature("neighbourhood"), spplotNeighbourhood)

## calculate neighbourhood from SpatialPointsDataFrame

# returns an neighbourhood object
# spData	spatialPointsDataFrame
# variable 		one or multiple variable names, all is the default
# size		the size of the neighbourhood, default of 5
# dep		denoting a subset of dependent locations (default NULL: all locations will be used)
# indep		denoting a subset of independent locations (default NULL: all locations will be used)
#		no location will be paired with itself
getNeighbours <- function(spData,var=names(spData),size=4,dep=NULL,indep=NULL,min.dist=10){
nLocs <- length(spData)
distMat <- spDists(spData)
if(min.dist>0) distMat[distMat<min.dist] <- Inf
if ( any(is.na( match(var,names(spData)) )) ) 
  stop("At least one of the variables is unkown is not part of the data.")
if(is.null(dep) & !is.null(indep))   dep <- 1:nLocs[-indep]
if(!is.null(dep) & is.null(indep)) indep <- 1:nLocs[-dep]
if(!is.null(dep) & !is.null(indep)) {
  cat("Reduced distance matrix is used: (",dep,") x (",indep,")",sep="")
} else {
  dep <- 1:nLocs
  indep <- 1:nLocs
}

size <- min(size,length(indep)-1)
if (size > sizeLim) {
  stop(paste("Evaluation of copulas might take a long time for more than",
         sizeLim," neighbours. Increase sizeLim if you want to evaluate neighbourhoods with",
	 size,"locations."))
}

lData <- vector("list",size)
index <- NULL
dists <- NULL

for (i in dep) {
  nbrs <- matrix(Inf,ncol=2,nrow=size-1)
  ind <- logical(nLocs)
  ind[indep] <- TRUE
  ind[i] <- FALSE
  for (j in (1:nLocs)[ind]) {
    tmpDist <- distMat[i,j]
    if (any(tmpDist < nbrs[,1])) {
      nbrs[size-1,] <- c(tmpDist,j)
      nbrs <- nbrs[order(nbrs[,1]),]
    }
  }
  lData[[1]] <- rbind(lData[[1]],spData@data[i,var,drop=FALSE])
  for (nbr in 2:size) {
    lData[[nbr]] <- rbind(lData[[nbr]],spData@data[nbrs[nbr-1,2],var,drop=FALSE])
  }
  index <- rbind(index, c(i,nbrs[,2]))
  dists <- rbind(dists, c(nbrs[,1]))
}

return(neighbourhood(lData, dists, SpatialPoints(spData), index))
}

## testing ## 
# data(meuse)
# coordinates(meuse) <- ~x+y
# str(SpatialPoints(meuse))
# 
# neighbourSet <- getNeighbours(meuse,size=5)
# 
# str(neighbourSet@index)
# 
# library(lattice)
# spplot(neighbourSet,"zinc",col.regions=bpy.colors())

# 
# as.data.frame(array(1:18,dim=c(2,3,3)))
# 
# ?spplot
# 
# array(matrix(runif(3*155),ncol=3),dim=c(155,3))
# 
# names(neighb)
# 
# ## transformation of the sample by local neighborhoods ##
# 

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

# the generic calcBins, calculates bins for spatiaql and spatio-temporal data

setGeneric("calcBins", function(data, nbins=15, boundaries=NA, cutoff=500000, ...) standardGeneric("calcBins") )

# calculating the spatial bins
# 
# data denotes the spatial data object
# var denotes the only variable name used
# cor.method is passed on to cor() (default="kendall")
# if plot=TRUE (default), the correlation measures are plotted agaisnt the mean lag separation distance
# 
calcSpBins <- function(data, var, nbins=15, boundaries=NA, cutoff=NA, cor.method="kendall", plot=TRUE) {

  if(is.na(boundaries)) {
    diagonal <- spDists(coordinates(t(data@bbox)))[1,2]
    boundaries <- ((1:nbins) * min(cutoff,diagonal/3,na.rm=T) / nbins)
  }
  
  lags <- calcSpLagInd(data, boundaries)
    
  mDists <- sapply(lags,function(x) mean(x[,3]))
  lagData <- lapply(lags, function(x) as.matrix((cbind(data[x[,1],var]@data, data[x[,2],var]@data))))
  
  lagCor <- sapply(lagData,function(x) cor(x,method=cor.method)[1,2])
  
  if(plot) { 
    plot(mDists, lagCor, xlab="distance",ylab=paste("correlation [",cor.method,"]",sep=""), 
         ylim=1.05*c(-abs(min(lagCor)),max(lagCor)), xlim=c(0,max(mDists)))
    abline(h=c(-min(lagCor),0,min(lagCor)),col="grey")
  }
  
  return(list(meanDists = mDists, lagCor=lagCor, lagData=lagData, lags=lags))
}

setMethod(calcBins, signature("Spatial"), calcSpBins)

# instances: number  -> number of randomly choosen temporal intances
#            NA      -> all observations
#            other   -> temporal indexing as in spacetime/xts, the parameter t.lags is set to 0 in this case.
# t.lags:    numeric -> temporal shifts between obs
calcStBins <- function(data, variable="PM10", nbins=15, boundaries=NA, cutoff=NA, instances=10, t.lags=c(0), cor.method="kendall", plot=TRUE) {

  if(is.na(boundaries)) {
    diagonal <- spDists(coordinates(t(data@sp@bbox)))[1,2]
    boundaries <- ((1:nbins) * min(cutoff,diagonal/3,na.rm=T) / nbins)
  }

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
                          as.matrix((cbind(data[spIndex[,1], tempIndices[,2*i-1], variable]@data, 
                                           data[spIndex[,2], tempIndices[,2*i], variable]@data))))
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
      cors <- c(cors, CDVine:::fasttau(binnedData[,2*i-1], binnedData[,2*i]))
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
  
  return(list(meanDists = mDists, lagCor=lagCor, lagData=lagData, lags=list(sp=spIndices, time=tempIndices)))
}

setMethod(calcBins, signature("STFDF"), calcStBins)
