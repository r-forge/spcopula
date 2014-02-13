###############################################################
##                                                           ##
## functions based on spacetime preparing the use of copulas ##
##                                                           ##
###############################################################

## spatio-temporal neighbourhood constructor
############################################

stNeighbourhood <- function(data, distances, index, var, coVar=character(), prediction=FALSE) {
  data <- as.data.frame(data)
  sizeN <- nrow(data)
  
  dimDists <- dim(distances)
  dimInd <- dim(index)
  
  stopifnot(length(dimDists) == 3)
  stopifnot(length(dimInd) == 3)
  
  stopifnot(dimDists[1] == sizeN)
  stopifnot(dimInd[1] == dimDists[1])
  
  stopifnot(((dimDists[2] + !prediction) + length(coVar)) == ncol(data))
  stopifnot(dimInd[2] == dimDists[2]+1)
  
  stopifnot(dimDists[3] == 2)
  stopifnot(dimInd[3] == dimDists[3])
  
  colnames(data) <- paste(paste("N", (0+prediction):dimDists[2], sep=""), var, sep=".")
  
  if(length(coVar)>0)
    colnames(data)[dimDists[2] + 1:length(coVar)] <- paste("N0", coVar)
  
  if (anyDuplicated(rownames(data))>0)
    rownames <- 1:length(rownames)
  
  new("stNeighbourhood", data=data, distances=distances, index=index,
      var=var, coVar=coVar, prediction=prediction)
}

## show
showStNeighbourhood <- function(object) {
  cat("A set of spatio-temporal neighbourhoods consisting of", dim(object@distances)[2]+1, "locations each \n")
  cat("with",nrow(object@data),"rows of observations for:\n")
  cat(object@var,"\n")
  if(length(object@coVar > 0))
    cat("with covariate", object@coVar, "\n")
}

setMethod(show, signature("stNeighbourhood"), showStNeighbourhood)

# select
selectFromStNeighbourhood <- function(x, i) {
  new("stNeighbourhood", data=x@data[i,,drop=F], 
      distances=x@distances[i,,,drop=F], index=x@index[i,,,drop=F], 
      var=x@var, coVar=x@coVar, prediction=x@prediction)
}

setMethod("[", signature("stNeighbourhood","numeric"), selectFromStNeighbourhood) 

## calculate neighbourhood from ST

# returns an neighbourhood object
##################################

getStNeighbours <- function(stData, ST, spSize=4, t.lags=-(0:2),
                            var=names(stData@data)[1], coVar=character(),
                            timeSteps=NA, prediction=FALSE, min.dist=0.01) {
  stopifnot((!prediction && missing(ST)) || (prediction && !missing(ST)))
  stopifnot(min.dist>0 || prediction)
  
  timeSpan <- min(t.lags)
  if(missing(ST) && !prediction)
    ST=geometry(stData)
  
  stopifnot(is(ST,"ST"))
  
  if(any(is.na(match(var,names(stData@data)))))
    stop("The variables is not part of stData.")
  if(length(coVar)>0)
    if(any(is.na(match(coVar,names(stData@data)))))
      stop("The covariate is not part of stData.")
  
  if(!prediction) {
    if(is.na(timeSteps)) {
      timeSteps <- length(stData@time)+timeSpan
      reSample <- function() (1-timeSpan):length(stData@time)
    } else {
      reSample <- function() sort(sample((1-timeSpan):length(stData@time), timeSteps))
    }
    nLocs <- length(ST@sp)*timeSteps
    nghbrs <- getNeighbours(dataLocs=geometry(stData@sp), var=character(), size=spSize,
                            min.dist=min.dist)
  } else {
    nLocs <- length(ST)
    nghbrs <- getNeighbours(dataLocs=geometry(stData@sp), predLocs=geometry(ST@sp),
                            size=spSize, var=character(), prediction=prediction, 
                            min.dist=min.dist)
    timeNghbrs <- sapply(index(ST@time), function(x) which(x == index(stData@time)))
    reSample <- function() timeNghbrs
    timeSteps <- length(stData@time)+timeSpan
  }
  
  nStNeighs <- (spSize-1)*length(t.lags)
  
  stNeighData <- matrix(NA, nLocs, nStNeighs + 1 + length(coVar))
  stDists <- array(NA,c(nLocs, nStNeighs, 2))
  stInd <- array(NA,c(nLocs, nStNeighs + 1, 2))
  
  nTimeInst <- length(reSample())
  
  for (i in 1:nrow(nghbrs@index)) {
    timeInst <- reSample() # draw random time steps for each neighbourhood
    spInd <- (i-1)*timeSteps+(1:timeSteps)
    
    stNeighData[spInd, 1:spSize] <- matrix(stData[nghbrs@index[i,], timeInst,
                                                  var, drop=F]@data[[1]],
                                           ncol=spSize, byrow=T)
    # add covariate(s) to the last column(s)
    if (length(coVar) > 0) {
      coVarCols <- nStNeighs + 1 + (1:length(coVar))
      stNeighData[spInd, coVarCols] <- matrix(stData[nghbrs@index[i,1], timeInst,
                                                     coVar, drop=F]@data[[1]],
                                              ncol=length(coVar), byrow=T)
    }
    
    tmpInd <- matrix(rep(timeInst, spSize), ncol=spSize)
    
    for (j in 2:length(t.lags)) {
      t <- t.lags[j]
      stNeighData[spInd, (j-1)*(spSize-1)+2:(spSize)] <- matrix(stData[nghbrs@index[i,][-1],
                                                                       timeInst+t, var, drop=F]@data[[1]],
                                                                ncol=spSize-1, byrow=T)
      tmpInd <- cbind(tmpInd, matrix(rep(timeInst+t,spSize-1), ncol=spSize-1))
    }
    
    # store spatial distances
    stDists[spInd,,1] <- matrix(rep(nghbrs@distances[i,], timeSteps*length(t.lags)),
                                byrow=T, ncol=nStNeighs)
    
    # store temporal distances
    stDists[spInd,,2] <- matrix(rep(rep(t.lags,each=spSize-1), timeSteps),
                                byrow=T, ncol=nStNeighs)  
    
    # store space indices
    stInd[spInd,,1] <- matrix(rep(c(nghbrs@index[i, ], rep(nghbrs@index[i, -1], length(t.lags)-1)),
                                  timeSteps), ncol = nStNeighs + 1, byrow = T)
    
    # store time indices
    stInd[spInd,,2] <- tmpInd
  }

  if (prediction) {
    stNeighData <- stNeighData[,-1]
  } else {
    dataLocs <- NULL
  }
  return(stNeighbourhood(as.data.frame(stNeighData), stDists, stInd, var, coVar, prediction))
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
  nrCM <- nrow(highCorMat)
  
  stNeighDataRed <- matrix(NA, nrow=nrCM, ncol=n+1+length(stNeigh@coVar))
  stNeighDistRed <- array(NA, dim=c(nrCM, n, 2))
  stNeighIndeRed <- array(NA, dim=c(nrCM, n+1, 2))
  if (length(stNeigh@coVar) > 0) {
    for (i in 1:nrCM) {
      selCol <- c(1,highCorMat[i,]+1, ncol(stNeigh@data)-((length(stNeigh@coVar)-1):0))
      stNeighDataRed[i,] <- as.numeric(stNeigh@data[i,selCol])
      stNeighDistRed[i,,] <- stNeigh@distances[i,highCorMat[i,],]
      stNeighIndeRed[i,,] <- stNeigh@index[i,c(1,highCorMat[i,]+1),]
    }
  } else {
    for (i in 1:nrCM) {
      stNeighDataRed[i,] <- as.numeric(stNeigh@data[i,c(1,highCorMat[i,]+1)])
      stNeighDistRed[i,,] <- stNeigh@distances[i,highCorMat[i,],]
      stNeighIndeRed[i,,] <- stNeigh@index[i,c(1,highCorMat[i,]+1),]
    }
  }
  
  boolNA <- !is.na(stNeigh@data[[1]])
  stNeighDataRed <- stNeighDataRed[boolNA,]
  stNeighDistRed <- stNeighDistRed[boolNA,,]
  stNeighIndeRed <- stNeighIndeRed[boolNA,,]
  
  return(stNeighbourhood(stNeighDataRed, stNeighDistRed, stNeighIndeRed,
                         var=stNeigh@var, coVar=stNeigh@coVar,
                         prediction=stNeigh@prediction))
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


# instances: number  -> number of randomly choosen temporal intances
#            NA      -> all observations
#            other   -> temporal indexing as in spacetime/xts, the parameter t.lags is set to 0 in this case.
# t.lags:    numeric -> temporal shifts between obs
calcStBins <- function(data, var, nbins=15, boundaries=NA, cutoff=NA, 
                       instances=NA, t.lags=-(0:2), ...,
                       cor.method="fasttau", plot=FALSE) {
  if(is.na(cutoff)) 
    cutoff <- spDists(coordinates(t(data@sp@bbox)))[1,2]/3
  if(is.na(boundaries)) 
    boundaries <- ((1:nbins) * cutoff / nbins)
  if(is.na(instances)) 
    instances=length(data@time)
  
  spIndices <- calcSpLagInd(data@sp, boundaries)
  
  mDists <- sapply(spIndices,function(x) mean(x[,3]))
  
  lengthTime <- length(data@time)
  if (!is.numeric(instances) | !length(instances)==1) {
    tempIndices <- cbind(instances, instances)
  } 
  else {
    tempIndices <- NULL
    for (t.lag in rev(t.lags)) {
      #       smplInd <- max(1,1-min(t.lags)):min(lengthTime,lengthTime-min(t.lags))
      smplInd <- sample(x=max(1,1-min(t.lags)):min(lengthTime,lengthTime-min(t.lags)),
                        size=min(instances,lengthTime-max(abs(t.lags))))
      tempIndices <- cbind(smplInd+t.lag, tempIndices)
      tempIndices <- cbind(smplInd, tempIndices)
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
      tmpData <- binnedData[,2*i+c(-1,0)]
      tmpData <- tmpData[!apply(tmpData, 1, function(x) any(is.na(x))),]
      cors <- c(cors, TauMatrix(tmpData)[1,2])
    }
    return(cors)
  }
  
  calcCor <- switch(cor.method, fasttau=calcTau, calcStats)
  
  lagCor <- sapply(lagData, calcCor)
  
  if(plot) { 
    plot(mDists, as.matrix(lagCor)[1,], xlab="distance",ylab=paste("correlation [",cor.method,"]",sep=""), 
         ylim=1.05*c(-abs(min(lagCor)),max(lagCor)), xlim=c(0,max(mDists)))
    abline(h=c(-min(lagCor),0,min(lagCor)),col="grey")
  }
  
  res <- list(meanDists = mDists, lagCor=lagCor, lagData=lagData, lags=list(sp=spIndices, time=tempIndices))
  attr(res,"cor.method") <- cor.method
  return(res)
}

setMethod(calcBins, signature(data="STFDF"), calcStBins)