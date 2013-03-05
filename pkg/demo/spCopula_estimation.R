## librarys ##
library(spcopula)
library(evd)

## meuse - spatial poionts data.frame ##
data(meuse)
coordinates(meuse) = ~x+y

spplot(meuse,"zinc", col.regions=bpy.colors(5))

## margins ##
hist(meuse[["zinc"]],freq=F,n=30,ylim=c(0,0.0035), 
     main="Histogram of zinc", xlab="zinc concentration")
gevEsti <- fgev(meuse[["zinc"]])$estimate
meanLog <- mean(log(meuse[["zinc"]]))
sdLog <- sd(log(meuse[["zinc"]]))
curve(dgev(x,gevEsti[1], gevEsti[2], gevEsti[3]),add=T,col="red")
curve(dlnorm(x,meanLog,sdLog),add=T,col="green")

ks.test(meuse[["zinc"]],pgev,gevEsti[1], gevEsti[2], gevEsti[3]) # p: 0.07
ks.test(meuse[["zinc"]],plnorm,meanLog,sdLog) # p: 0.03

pMar <- function(q) plnorm(q, meanLog, sdLog)
qMar <- function(p) qlnorm(p, meanLog, sdLog)
dMar <- function(x) dlnorm(x, meanLog, sdLog)

# pMar <- function(q) pgev(q, gevEsti[1], gevEsti[2], gevEsti[3])
# qMar <- function(p) qgev(p, gevEsti[1], gevEsti[2], gevEsti[3])
# dMar <- function(x) dgev(x, gevEsti[1], gevEsti[2], gevEsti[3])

## lag classes ##
bins <- calcBins(meuse,var="zinc",nbins=10,cutoff=800)

# transform data to the unit interval
bins$lagData <- lapply(bins$lagData, rankTransform)

## calculate parameters for Kendall's tau function ##
# either linear
calcKTauLin <- fitCorFun(bins, degree=1, cutoff=600)
curve(calcKTauLin,0, 1000, col="red",add=TRUE)

# or polynomial (used here)
calcKTauPol <- fitCorFun(bins, degree=3)
curve(calcKTauPol,0, 1000, col="purple",add=TRUE)

## find best fitting copula per lag class
loglikTau <- loglikByCopulasLags(bins, calcKTauPol,
                                 families=c(normalCopula(0), tCopula(0),
                                            claytonCopula(0), frankCopula(1), 
                                            gumbelCopula(1), joeBiCopula(1.5),
                                            indepCopula()))
bestFitTau <- apply(apply(loglikTau, 1, rank, na.last=T), 2, 
                    function(x) which(x==7))
bestFitTau

## set-up a spatial Copula ##
spCop <- spCopula(components=list(normalCopula(0), tCopula(0),
                                  frankCopula(1), normalCopula(0), 
                                  claytonCopula(0), claytonCopula(0), 
                                  claytonCopula(0), claytonCopula(0),
                                  claytonCopula(0), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")

## compare spatial copula loglik by lag:
spLoglik <- NULL
for(i in 1:length(bins$lags)) { # i <- 8
  spLoglik <- c(spLoglik,
                sum(dCopula(u=bins$lagData[[i]], spCop,log=T,
                            h=bins$lags[[i]][,3])))
}

plot(spLoglik, ylab="log-likelihood", xlim=c(1,11)) 
points(loglikTau[cbind(1:10,bestFitTau)], col="green", pch=16)
points(loglikTau[,1], col="red", pch=5)
legend(6, 50,c("Spatial Copula", "best copula per lag", "Gaussian Copula",
               "number of pairs"), 
       pch=c(1,16,5,50), col=c("black", "green", "red"))
text(x=(1:10+0.5),y=spLoglik,lapply(bins$lagData,length))

##
# spatial vine
vineDim <- 5L
meuseNeigh <- getNeighbours(meuse,"zinc",vineDim)
meuseNeigh@data <- rankTransform(meuseNeigh@data)

meuseSpVine <- fitCopula(spVineCopula(spCop, vineCopula(as.integer(vineDim-1))),
                         meuseNeigh)

meuseSpVine@vineCop

##
# leave-one-out x-validation

condVine <- function(condVar, dists, n=100) {
  rat <- 0.2/(1:(n/2))-(0.1/((n+1)/2))
  xVals <- unique(sort(c(rat,1-rat,1:(n-1)/(n))))
  xLength <- length(xVals)
  repCondVar <- matrix(condVar, ncol=length(condVar), nrow=xLength, byrow=T)
  density <- dCopula(cbind(xVals, repCondVar), meuseSpVine, h=dists)
  
  linAppr <- approxfun(c(0,xVals,1), density[c(1,1:xLength,xLength)] ,yleft=0, yright=0)
  int <- integrate(linAppr,lower=0, upper=1)$value
  
  return(function(u) linAppr(u)/int)
}

time <- proc.time()  # ~30 s
predMedian <- NULL
predMean <- NULL
for(loc in 1:nrow(meuseNeigh@data)) { # loc <- 429  predNeigh$data[loc,1]
  cat("Location:",loc,"\n")
  condSecVine <- condVine(condVar=as.numeric(meuseNeigh@data[loc,-1]), 
                          dists=meuseNeigh@distances[loc,,drop=F])
  
  predMedian <- c(predMedian, qMar(optimise(function(x) abs(integrate(condSecVine,0,x)$value-0.5),c(0,1))$minimum))
  
  condExp <-  function(x) {
    condSecVine(pMar(x))*dMar(x)*x
  }
  
  predMean <- c(predMean, integrate(condExp,0,3000,subdivisions=1e6)$value)
}
proc.time()-time

mean(abs(predMean-meuse$zinc))
mean(predMean-meuse$zinc)
sqrt(mean((predMean-meuse$zinc)^2))

mean(abs(predMedian-meuse$zinc))
mean(predMedian-meuse$zinc)
sqrt(mean((predMedian-meuse$zinc)^2))

plot(predMean,meuse$zinc)
abline(0,1)

plot(predMedian,meuse$zinc)
abline(0,1)

## kriging results:
# same neighbourhood size 5L:
# MAE:  158.61
# BIAS:  -4.24
# RMSE: 239.85
#
# global kriging:
# MAE:  148.85
# BIAS:  -3.05
# RMSE: 226.15