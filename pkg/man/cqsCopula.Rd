\name{cqsCopula}
\alias{cqsCopula}
\title{
Constructor of a symmetric copula with cubic quadratic sections. 
}
\description{
A copula from the copula family of a symmetric copula with cubic and quadratic sections (\code{\linkS4class{cqsCopula}}) for a given pair of parameters \code{param} (Nelsen 2006).
}
\usage{
cqsCopula(param=c(0,0), fixed=character(0))
}

\arguments{
  \item{param}{
The copula is defined by the two-place parameter \code{param}. The default represents the product copula. 
}
\item{fixed}{a single character referring to "a" or "b" fixing one of the parameters. This allows to estimate the remaining parameter through inversion of Kendall's tau \code{\link{iTau}} or Spearman's rho \code{\link{iRho}}.}
}

\value{
An \code{\linkS4class{cqsCopula}} object that behaves (almost) like any other \code{\linkS4class{copula}} from the \code{\link{copula-package}}.
}
\references{
Nelsen, Roger B. (2006): An Introduction to Copulas, second edition, Springer
}
\author{
Benedikt Graeler
}
\seealso{
See also \code{\link{asCopula}}.
}
\examples{
persp(cqsCopula(c(-2,1)),dCopula)
}
\keyword{copula}
\keyword{cubic quadratic sections}
