#################################################################################
##
##   R package spCopula by Benedikt Gräler Copyright (C) 2009
##
##   This file is part of the R package spCopula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
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
# uniform="logical"	is the data distributed on [0,1]

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
# var 		one or multiple variable names, all is the default
# size		the size of the neighbourhood, default of 5
# dep		denoting a subset of dependent locations (default NULL: all locations will be used)
# indep		denoting a subset of independent locations (default NULL: all locations will be used)
#		no location will be paired with itself
getNeighbours <- function(spData,var=names(spData),size=4,dep=NULL,indep=NULL){
nLocs <- length(spData)
distMat <- spDists(spData)
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
data(meuse)
coordinates(meuse) <- ~x+y
str(SpatialPoints(meuse))

neighbourSet <- getNeighbours(meuse,size=5)

str(neighbourSet@index)

library(lattice)
spplot(neighbourSet,"zinc",col.regions=bpy.colors())

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
