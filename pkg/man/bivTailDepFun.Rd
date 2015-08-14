\name{bivTailDepFun}
\alias{bivTailDepFun}
\alias{lowerBivTailDepFun}
\alias{upperBivTailDepFun}

\alias{empBivTailDepFun}
\alias{lowerEmpBivTailDepFun}
\alias{upperEmpBivTailDepFun}

\title{
Tail dependence functions
}
\description{
Functions returning a (empirical) tail dependence function for a bivariate copula (sample). The tail dependence functions can be upper, lower or joint, where values below 0.5 are calculated from the lower tail dependence function and values larger 0.5 for the upper tail dependence function. The definition follows the one by Nelsen (2006, Theorem 5.4.2.). For tail dependence of multivariate copulas, see the papers by Joe et al. (2010) and Nikoloulopoulos et al. (2012).
}
\usage{
bivTailDepFun(copula)
lowerBivTailDepFun(copula)
upperBivTailDepFun(copula)

empBivTailDepFun(u)
lowerEmpBivTailDepFun(u)
upperEmpBivTailDepFun(u)
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
Definition: 

Nelsen, Roger B. An introduction to copulas. Vol. 139. Springer Science & Business Media, 2013.

Further reading on multivariate tail dependence:

Joe, Harry, Haijun Li, and Aristidis K. Nikoloulopoulos. "Tail dependence functions and vine copulas." Journal of Multivariate Analysis 101.1 (2010): 252-270.

Nikoloulopoulos, Aristidis K., Harry Joe, and Haijun Li. "Vine copulas with asymmetric tail dependence and applications to financial return data." Computational Statistics & Data Analysis 56.11 (2012): 3659-3673.

Plots inspired by:

\url{http://freakonometrics.blog.free.fr/index.php?post/2012/09/11/Copulas-and-statistical-inference}
}
\author{
Benedikt Graeler
}

\examples{
library("VineCopula")
data("simulatedTriples")
X <- rankTransform(triples[,c(1,3)])
  
tdfEmp <- empBivTailDepFun(X)
plot(tdfEmp,ylim=c(0,1),
     ylab="tail dependence index")  
abline(v=0.5, col="grey")

smplTau <- cor(X,method="kendall")[1,2]

# Gauss
tdfGauss <- bivTailDepFun(normalCopula(sin(smplTau*pi/2)))
curve(tdfGauss,add=TRUE,col="blue")

# survival Gumbel
tdfGumbel <- bivTailDepFun(surGumbelCopula(1/(1-smplTau)))
curve(tdfGumbel,add=TRUE,col="darkgreen")

# survival BB6 copula
tdfBB6 <- bivTailDepFun(surBB6Copula(c(4.65,2.28)))
curve(tdfBB6,add=TRUE,col="red")

legend("bottomleft",c("empircal","Gauss","surv. Gumbel","surv. BB6"),
       col=c("black","blue","darkgreen","red"),lty=1)
} 