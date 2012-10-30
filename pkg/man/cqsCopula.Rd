\name{cqsCopula}
\alias{cqsCopula}
\title{
Constructor of a symmetric copula with cubic quadratic sections. 
}
\description{
A copula from the copula family of a symmetric copula with cubic and quadratic sections (\code{\linkS4class{cqsCopula}}) for a given pair of parameters \code{param} (Nelsen 2006).
}
\usage{
cqsCopula(param)
}

\arguments{
  \item{param}{
The copula is defined by the two-place parameter \code{param}.
}
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
