library(spcopula)
data(simulatedTriples)

rtPair <- 1-as.matrix(rankTransform(triples[,c(1,3)]))

plot(rtPair,asp=1)

tdfEmp <- empTailDepFun(rtPair)
plot(tdfEmp,ylim=c(0,1), ylab="tail index", xlab="u")
abline(v=0.5, col="grey")

gaussCop <- fitCopula(normalCopula(0), rtPair)@copula
tdfGauss <- tailDepFun(gaussCop)
curve(tdfGauss, add=T,col="green",n=500)

gumbelCop <- fitCopula(gumbelCopula(2),rtPair)@copula
tdfGumbel <- tailDepFun(gumbelCop)
curve(tdfGumbel,add=T, col="blue",n=500)

BB6Cop <- fitCopula(BB6Copula(), rtPair)@copula
tdfBB6 <- tailDepFun(BB6Cop)
curve(tdfBB6, add=T,col="red",n=500)

legend("bottomright",
       c("empirical", "Gaussian", "Gumbel", "BB6"),
       col=c("black", "green", "blue", "red"), lty=1)