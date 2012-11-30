## get the data
data(simulatedTriples)

## rank order transformation
peakVol <- rankTransform(triples[,1],triples[,3])
colnames(peakVol) <- c("Qp","Vp")
plot(peakVol, asp=1)

# Kendall's tau correlation
cor(triples,method="kendall")

# estiamte the BB7 copula by means of maximum likelihood
copQV <- fitCopula(BB7Copula(param=c(2,14)), peakVol, method="ml",start=c(2,14), estimate.variance=F)@copula
copQV

# we use a design return period of 100 years
# the MAR-case: given the 0.99 quantile of one marginal, what is the 
# corresponding quantile of the second variable with the given joint return 
# period of 1/(1-0.99)
v_MAR <- c(0.99,invdduCopula(0.99,copQV,0.99))
v_MAR

## the anlytical kendall distribution
kendallFunQV <- getKendallDistr(copQV)

# the Kendall distribution value (fraction of pairs having a smaller copula value than "t")
kendallFunQV(t=0.99)

curve(kendallFunQV, from=0, to=1, asp=1)

# the critical level of the KEN2-RP for 100 years
t_KEN2 <- criticalLevel(kendallFunQV,100,mu=1)
t_KEN2

# the corresponmding KEN2-RP for the OR-RP
kendallRP(kendallFunQV,cl=0.99,mu=1,copula=copQV)

# illustrating the critical lines (empirically)
contour(copQV,pCopula,levels=c(0.99,t_KEN2),
        xlim=c(0.98,1), ylim=c(0.98,1), n=1000, asp=1, col="blue")