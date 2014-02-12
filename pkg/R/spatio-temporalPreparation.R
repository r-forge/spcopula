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
  new("stNeighbourhood", data=data, distances=distances, locations=ST, 
      dataLocs=STxDF, index=index, prediction=prediction, var=var,
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


## calculate neighbourhood from ST

# returns an neighbourhood object
##################################

getStNeighbours <- function(stData, ST, var=names(stData@data)[1], spSize=4, 
                            t.lags=-(0:2), timeSteps=NA, prediction=FALSE, min.dist=0.01) {
  stopifnot((!prediction && missing(ST)) || (prediction && !missing(ST)))
  stopifnot(min.dist>0 || prediction)
  
  timeSpan <- min(t.lags)
  if(missing(ST) && !prediction)
    ST=stData
  
  stopifnot(is(ST,"ST"))
  
  if(any(is.na(match(var,names(stData@data)))))
    stop("At least one of the variables is unkown or is not part of the data.")
  
  if(!prediction) {
    if(is.na(timeSteps)) {
      timeSteps <- length(stData@time)+timeSpan
      reSample <- function() (1-timeSpan):length(stData@time)
    } else {
      reSample <- function() sort(sample((1-timeSpan):length(stData@time), timeSteps))
    }
    nLocs <- length(ST@sp)*timeSteps
    nghbrs <- getNeighbours(stData[,1], var=var, size=spSize, min.dist=min.dist)
  } else {
    nLocs <- length(ST)
    nghbrs <- getNeighbours(stData[,1], ST@sp, var, spSize, prediction, min.dist)
    timeNghbrs <- sapply(index(ST@time), function(x) which(x == index(stData@time)))
    reSample <- function() timeNghbrs
    timeSteps <- length(stData@time)+timeSpan
  }
  
  stNeighData <- matrix(NA, nLocs, (spSize-1)*length(t.lags)+1)
  stDists <- array(NA,c(nLocs,(spSize-1)*length(t.lags),2))
  stInd <- array(NA,c(nLocs,(spSize-1)*length(t.lags),2))
  
  nTimeInst <- length(reSample())
  
  for(i in 1:nrow(nghbrs@index)){ # i <- 1
    timeInst <- reSample() # draw random time steps for each neighbourhood
    stNeighData[(i-1)*timeSteps+(1:timeSteps),
                1:spSize] <- matrix(stData[nghbrs@index[i,], timeInst,
                                           var, drop=F]@data[[1]],
                                    ncol=spSize, byrow=T) # retrieve the top level data
    tmpInd <- matrix(rep(timeInst, spSize-1), ncol=spSize-1)
    for(j in 2:length(t.lags)) {
      t <- t.lags[j]
      stNeighData[(i-1)*timeSteps+(1:timeSteps),
                  (j-1)*(spSize-1)+2:(spSize)] <- matrix(stData[nghbrs@index[i,][-1],
                                                                timeInst+t,
                                                                var, drop=F]@data[[1]],
                                                         ncol=spSize-1, byrow=T)
      tmpInd <- cbind(tmpInd, matrix(rep(timeInst+t,spSize-1),ncol=spSize-1))
    }

    stDists[(i-1)*timeSteps+1:timeSteps,,1] <- matrix(rep(nghbrs@distances[i,],
                                                          timeSteps*length(t.lags)),
                                                      byrow=T, ncol=length(t.lags)*(spSize-1))   # store sp distances
    stDists[(i-1)*timeSteps+1:timeSteps,,2] <- matrix(rep(rep(t.lags,each=spSize-1),
                                                          timeSteps),
                                                      byrow=T, ncol=length(t.lags)*(spSize-1))  # store tmp distances
    stInd[(i-1)*timeSteps+1:timeSteps,,1] <- matrix(rep(nghbrs@index[i,][-1],
                                                    timeSteps*length(t.lags)),
                                                    byrow=T, ncol=length(t.lags)*(spSize-1))
    stInd[(i-1)*timeSteps+1:timeSteps,,2] <- tmpInd
  }

  if (prediction) {
    dataLocs <- stData
    stNeighData <- stNeighData[,-1]
  } else {
    dataLocs <- NULL
  }
  return(stNeighbourhood(as.data.frame(stNeighData), stDists, stData, ST, 
                         stInd, prediction, var))
}


## reduction of a larger neigbopurhood based on correlation strengths
reduceNeighbours <- function(stNeigh, stDepFun, n) {
  stopifnot(n>0)
  
  # transform distances into correlations to detect the strongest correlated ones
  dimStNeigh <- dim(stNeigh@distances)
  corMat <- matrix(NA, dimStNeigh[1], dimStNeigh[2])
  
  for (i in 1:dimStNeigh[2]) {
    boolNA <- is.na(stNeigh@data[[1]]) | is.na(stNeigh@data[[1+i]])
    stNeigh@distances[boolNA,i,] <- c(NA,NA)
    tLag <- -1*stNeigh@distances[!boolNA,i,2][1]+1
    corMat[!boolNA,i] <- stDepFun(stNeigh@distances[!boolNA,i,1], tLag)
  }
  
  highCorMat <- t(apply(corMat, 1, function(x) order(x, na.last=TRUE, decreasing=TRUE)[1:n]))
  
  stNeighDataRed <- matrix(NA, nrow=nrow(highCorMat), ncol=n+1)
  stNeighDistRed <- array(NA, dim=c(nrow(highCorMat), n, 2))
  stNeighIndeRed <- array(NA, dim=c(nrow(highCorMat), n, 2))
  for (i in 1:nrow(highCorMat)) {
    stNeighDataRed[i,] <- as.numeric(stNeigh@data[i,c(1,highCorMat[i,]+1)])
    stNeighDistRed[i,,] <- stNeigh@distances[i,highCorMat[i,],]
    stNeighIndeRed[i,,] <- stNeigh@index[i,highCorMat[i,],]
  }
  
  stNeighDataRed <- stNeighDataRed[!is.na(stNeigh@data[[1]]),]
  stNeighDistRed <- stNeighDistRed[!is.na(stNeigh@data[[1]]),,]
  stNeighIndeRed <- stNeighIndeRed[!is.na(stNeigh@data[[1]]),,]
  
  return(stNeighbourhood(stNeighDataRed,stNeighDistRed, stNeigh@dataLocs, 
                         ST=stNeigh@dataLocs, stNeighIndeRed, prediction=F, 
                         var=stNeigh@var))
}

## to be redone
calcStNeighBins <- function(data, var="uniPM10", nbins=9, t.lags=-(0:2),
                            boundaries=NA, cutoff=NA, cor.method="fasttau") {
#   dists <- data@distances[,,1]
#   
#   corFun <- switch(cor.method,
#                    fasttau=function(x) VineCopula:::fasttau(x[,1],x[,2]),
#                    function(x) cor(x,method=cor.method)[1,2])
#   
#   if (any(is.na(boundaries))) 
#     boundaries <- quantile(as.vector(dists), probs=c(1:nbins/nbins))
#   if(!is.na(cutoff)) {
#     boundaries <- boundaries[boundaries < cutoff]
#     boundaries <- unique(c(0,boundaries,cutoff))
#   } else {
#     boundaries <- unique(c(0,boundaries))
#   }
#   
#   lagData <- NULL
#   for(t.lag in t.lags) { # t.lag <- 0
#     tBool <- data@distances[,,2]==t.lag
#     tmpLagData <- NULL
#     for(i in 1:nbins) { # i <- 1
#       sBool <- (dists <= boundaries[i + 1] & dists > boundaries[i])
#       bool <- tBool & sBool
#       pairs <- NULL
#       for (col in 1:(dim(tBool)[2])) { # col <- 1
#         if(!any(bool[, col]))
#           next
#         sInd <- data@index[bool[, col], c(1, 1 + col),1]
#         tInd <- data@index[bool[, col], c(1, 1 + col),2]
#         p1 <- apply(cbind(sInd[,1], tInd[,1]),1,
#                     function(x) data@locations[x[1], x[2],var])
#         p2 <- apply(cbind(sInd[,2], tInd[,2]),1,
#                     function(x) data@locations[x[1], x[2],var])
#         pairs <- rbind(pairs, cbind(p1,p2))
#       }
#       tmpLagData <- append(tmpLagData,list(pairs))
#     }
#     lagData <- append(lagData,list(tmpLagData))
#     
#   }
#   
#   lagData <- lapply(spIndices, retrieveData, tempIndices = tempIndices)
#   calcStats <- function(binnedData) {
#     cors <- NULL
#     for (i in 1:(ncol(binnedData)/2)) {
#       cors <- c(cors, cor(binnedData[, 2 * i - 1], binnedData[, 2 * i], method = cor.method, use = "pairwise.complete.obs"))
#     }
#     return(cors)
#   }
#   calcTau <- function(binnedData) {
#     cors <- NULL
#     for (i in 1:(ncol(binnedData)/2)) {
#       cors <- c(cors, VineCopula:::fasttau(binnedData[, 2 * i - 1], binnedData[, 2 * i]))
#     }
#     return(cors)
#   }
#   calcCor <- switch(cor.method, fasttau = calcTau, calcStats)
#   lagCor <- sapply(lagData, calcCor)
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   np <- numeric(0)
#   moa <- numeric(0)
#   lagData <- NULL
#   meanDists <- numeric(0)
#   
#   data <- as.matrix(data@data)
#   
#   for ( i in 1:nbins) {
#     bools <- (dists <= boundaries[i+1] & dists > boundaries[i])
#     
#     pairs <- NULL
#     for(col in 1:(dim(bools)[2])) {
#       pairs <- rbind(pairs, data[bools[,col],c(1,1+col)])
#     }
#     
#     lagData <- append(lagData, list(pairs))
#     moa <- c(moa, corFun(pairs))
#     meanDists <- c(meanDists, mean(dists[bools]))
#     np <- c(np, sum(bools))
#   }
#   
#   if(plot) { 
#     plot(meanDists, moa, xlab="distance", ylab=paste("correlation [",cor.method,"]",sep=""), 
#          ylim=1.05*c(-abs(min(moa)),max(moa)), xlim=c(0,max(meanDists)))
#     abline(h=c(-min(moa),0,min(moa)),col="grey")
#   }
#   
#   res <- list(np=np, meanDists = meanDists, lagCor=moa, lagData=lagData)
#   attr(res,"cor.method") <- switch(cor.method, fasttau="kendall", cor.method)
#   return(res)
}

setMethod(calcBins, signature="stNeighbourhood", calcStNeighBins)