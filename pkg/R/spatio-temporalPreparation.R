###############################################################
##                                                           ##
## functions based on spacetime preparing the use of copulas ##
##                                                           ##
###############################################################

## spatio-temporal neighbourhood constructor
############################################

stNeighbourhood <- function(data, distances, STxDF, ST=NULL,index, 
                          prediction, var) {
  data <- as.data.frame(data)
  sizeN <- nrow(data)
  dimDists <- dim(distances)
  
  stopifnot(dimDists[1]==sizeN)
  stopifnot(dimDists[2]==ncol(data)-(!prediction))
  stopifnot(dimDists[3]==2)
  colnames(data) <- paste(paste("N", (0+prediction):dimDists[2], sep=""),var,sep=".")
  if (anyDuplicated(rownames(data))>0)
    rownames <- 1:length(rownames)
  new("stNeighbourhood", data=data, distances=distances, locations=STxDF, 
      dataLocs=ST, index=index, prediction=prediction, var=var,
      sp=as(STxDF@sp, "Spatial"), time=STxDF@time[1], 
      endTime=STxDF@endTime[length(STxDF@endTime)])
}

## show
showStNeighbourhood <- function(object){
  cat("A set of spatio-temporal neighbourhoods consisting of", dim(object@distances)[2]+1, "locations each \n")
  cat("with",nrow(object@data),"rows of observations for:\n")
  cat(object@var,"\n")
}

setMethod(show,signature("stNeighbourhood"),showStNeighbourhood)


## calculate neighbourhood from SpatialPointsDataFrame

# returns an neighbourhood object
##################################
getStNeighbours <- function(stData, ST, var=names(stData@data)[1], spSize=4, 
                            t.lags=-(0:2), timeSteps=NA, prediction=FALSE, min.dist=0.01) {
  stopifnot((!prediction && missing(ST)) || (prediction && !missing(ST)))
  stopifnot(min.dist>0 || prediction)
  
  timeSpan <- min(t.lags)
  if(missing(ST) && !prediction)
    ST=stData
  if(is.na(timeSteps))
    timeSteps <- length(stData@time)+timeSpan
  
  stopifnot(is(ST,"ST"))
  
  nLocs <- length(ST@sp)*timeSteps
  
  if(any(is.na(match(var,names(stData@data)))))
    stop("At least one of the variables is unkown or is not part of the data.")

  if(prediction)
    nghbrs <- getNeighbours(stData[,1], ST, var, spSize, prediction, min.dist)
  else
    nghbrs <- getNeighbours(stData[,1], var=var, size=spSize, min.dist=min.dist)
  
  stNeighData <- NULL
  stDists <- array(NA,c(nLocs,(spSize-1)*length(t.lags),2))
  stInd <- array(NA,c(nLocs,(spSize-1)*length(t.lags),2))
  for(i in 1:nrow(nghbrs@index)){ # i <- 1
    tmpInst <- sample((1-timeSpan):length(stData@time), timeSteps) # draw random time steps for each neighbourhood
    tmpData <- matrix(stData[c(i, nghbrs@index[i,]),  tmpInst,  var]@data[[1]],
                      ncol=spSize, byrow=T) # retrieve the top level data
    tmpInd <- matrix(rep(tmpInst,spSize-1),ncol=spSize-1)
    for(t in t.lags[-1]) {
      tmpData <- cbind(tmpData, matrix(stData[nghbrs@index[i,], 
                                              tmpInst+t,var]@data[[1]],
                                       ncol=spSize-1, byrow=T))
      tmpInd <- cbind(tmpInd, matrix(rep(tmpInst+t,spSize-1),ncol=spSize-1))
    }
    stNeighData <- rbind(stNeighData, tmpData) # bind data row-wise
    stDists[(i-1)*timeSteps+1:timeSteps,,1] <- matrix(rep(nghbrs@distances[i,],
                                                          timeSteps*(spSize-1)),
                                                      byrow=T, ncol=length(t.lags)*(spSize-1))  # store sp distances
    stDists[(i-1)*timeSteps+1:timeSteps,,2] <- matrix(rep(rep(t.lags,each=spSize-1),
                                                          timeSteps),
                                                      byrow=T, ncol=length(t.lags)*(spSize-1))  # store tmp distances
    stInd[(i-1)*timeSteps+1:timeSteps,,1] <- matrix(rep(nghbrs@index[i,],
                                                    timeSteps*(spSize-1)),
                                                    byrow=T, ncol=length(t.lags)*(spSize-1))
    stInd[(i-1)*timeSteps+1:timeSteps,,2] <- tmpInd
  }

  if (prediction)
    dataLocs <- stData
  else 
    dataLocs <- NULL
  return(stNeighbourhood(as.data.frame(stNeighData), stDists, stData, ST, 
                         stInd, prediction, var))
}