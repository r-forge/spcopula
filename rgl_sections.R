library(copula)
library(VineCopula)
library(lattice)

grid <-  cbind(rep(1:99/100,99), rep(1:99/100,each=99))

plotData <- as.data.frame(grid)
plotData$sec1 <- dCopula(grid, joeBiCopula(param = 1.05))
plotData$sec2 <- dCopula(grid, normalCopula(.4))
plotData$sec3 <- dCopula(grid, claytonCopula(0.8))

colBreaks <- quantile(c(plotData$sec1,
                        plotData$sec2,
                        plotData$sec3), probs = 0:100/100)

p1 <- levelplot(sec1~V1+V2, plotData, at=colBreaks,
          col.regions = terrain.colors(102))
p2 <- levelplot(sec2~V1+V2, plotData, at=colBreaks,
          col.regions = terrain.colors(102))
p3 <- levelplot(sec3~V1+V2, plotData, at=colBreaks,
          col.regions = terrain.colors(102))

print(p1, position=c(0.2,0.65,0.8,1.05), more=T)
print(p2, position=c(0.2,0.3,0.8,0.7), more=T)
print(p3, position=c(0.2,-0.05,0.8,0.35))

##
library(rgl)
plotData$sl1 <- qnorm(0.2)
plotData$sl2 <- qnorm(0.5)
plotData$sl3 <- qnorm(0.8)

plot3d(plotData$sl1, plotData$V1, plotData$V2, 
       col=terrain.colors(102)[findInterval(plotData$sec1, colBreaks)])
plot3d(plotData$sl2, plotData$V1, plotData$V2, 
       col=terrain.colors(102)[findInterval(plotData$sec2, colBreaks)], add=T)
plot3d(plotData$sl3, plotData$V1, plotData$V2, 
       col=terrain.colors(102)[findInterval(plotData$sec3, colBreaks)], add=T)