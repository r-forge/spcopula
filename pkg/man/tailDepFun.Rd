\name{tailDepFun}
\alias{tailDepFun}
\alias{lowerTailDepFun}
\alias{upperTailDepFun}

\alias{empTailDepFun}
\alias{lowerEmpTailDepFun}
\alias{upperEmpTailDepFun}

\title{
Tail dependence functions
}
\description{
Functions returning a (empirical) tail dependence function for a copula (sample). The tail dependence functions can be upper, lower or joint, where values below 0.5 are calculated from the lower tail dependence function and values larger 0.5 for the upper tail dependence function.
}
\usage{
tailDepFun(copula)
lowerTailDepFun(copula)
upperTailDepFun(copula)

empTailDepFun(u)
lowerEmpTailDepFun(u)
upperEmpTailDepFun(u)
}

\arguments{
  \item{copula}{
an object of class \code{\linkS4class{copula}}
}
  \item{u}{
a bivariate sample on (0,1)
}
}

\value{
A function taking arguments from the unit interval (0,1) and returning the corresponding tail index.
}
\references{
Inspired by: \url{http://freakonometrics.blog.free.fr/index.php?post/2012/09/11/Copulas-and-statistical-inference}
}
\author{
Benedikt Graeler
}

\examples{
data(simulatedTriples)
X <- rankTransform(triples[,c(1,3)])
  
tdfEmp <- empTailDepFun(X)
plot(tdfEmp,ylim=c(0,1),
     ylab="tail dependence index")  
abline(v=0.5, col="grey")

smplTau <- cor(X,method="kendall")[1,2]

# Gauss
tdfGauss <- tailDepFun(normalCopula(sin(smplTau*pi/2)))
curve(tdfGauss,add=TRUE,col="blue")

# survival Gumbel
tdfGumbel <- tailDepFun(surGumbelCopula(1/(1-smplTau)))
curve(tdfGumbel,add=TRUE,col="darkgreen")

# survival BB6 copula
tdfBB6 <- tailDepFun(surBB6Copula(c(4.65,2.28)))
curve(tdfBB6,add=TRUE,col="red")

legend("bottomleft",c("empircal","Gauss","surv. Gumbel","surv. BB6"),
       col=c("black","blue","darkgreen","red"),lty=1)
} 