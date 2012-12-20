## librarys ##
library(evd)

## dataset - spatial poionts data.frame ##
data(meuse)
coordinates(meuse) = ~x+y
dataSet <- meuse

spplot(meuse,"zinc", col.regions=bpy.colors(5))

## margins ##
hist(dataSet[["zinc"]],freq=F,n=30,ylim=c(0,0.0035))
gevEsti <- fgev(dataSet[["zinc"]])$estimate
loc <- gevEsti[1]
scale <- gevEsti[2]
shape  <- gevEsti[3]
meanLog <- mean(log(meuse[["zinc"]]))
sdLog <- sd(log(meuse[["zinc"]]))
curve(dgev(x,loc,scale,shape),add=T,col="red")
curve(dlnorm(x,meanLog,sdLog),add=T,col="green")

ks.test(dataSet[["zinc"]],pgev,loc,scale,shape) # p: 0.07
ks.test(dataSet[["zinc"]],plnorm,meanLog,sdLog) # p: 0.03

## lag classes ##
bins <- calcBins(dataSet,var="zinc",nbins=10,cutoff=800)

# transform data to the unit interval
bins$lagData <- lapply(bins$lagData, function(x) cbind(rank(x[,1])/(nrow(x)+1),rank(x[,2])/(nrow(x)+1)))

## calculate parameters for Kendall's tau function ##
# either linear
calcKTauLin <- fitCorFun(bins, degree=1, cutoff=600)
curve(calcKTauLin,0, 1000, col="red",add=TRUE)

# or polynomial (used here)
calcKTauPol <- fitCorFun(bins, degree=3)
curve(calcKTauPol,0, 1000, col="purple",add=TRUE)

## find best fitting copula per lag class
loglikTau <- loglikByCopulasLags(bins, calcKTauPol,
                                 families=c(normalCopula(0), tCopula(0,dispstr = "un"),
                                            claytonCopula(0), frankCopula(1), 
                                            gumbelCopula(1), joeBiCopula(1.5),
                                            indepCopula()))
bestFitTau <- apply(apply(loglikTau, 1, rank),2,function(x) which(x==7))

## set-up a spatial Copula ##
spCop <- spCopula(components=list(normalCopula(0), tCopula(0, dispstr = "un"),
                                  frankCopula(1), normalCopula(0), claytonCopula(0),
                                  claytonCopula(0), claytonCopula(0), claytonCopula(0),
                                  claytonCopula(0), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")