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
                                                          timeSteps*length(t.lags)),
                                                      byrow=T, ncol=length(t.lags)*(spSize-1))   # store sp distances
    stDists[(i-1)*timeSteps+1:timeSteps,,2] <- matrix(rep(rep(t.lags,each=spSize-1),
                                                          timeSteps),
                                                      byrow=T, ncol=length(t.lags)*(spSize-1))  # store tmp distances
    stInd[(i-1)*timeSteps+1:timeSteps,,1] <- matrix(rep(nghbrs@index[i,],
                                                    timeSteps*length(t.lags)),
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