\name{asCopula}
\alias{asCopula}

\title{Constructor of an asymmetric copula with cubic and quadratic sections (Nelsen 2006).}
\description{A copula from the copula family of an asymmetric copula with cubic and quadratic sections for a given pair of parameters \code{param} (Nelsen 2006).}
\usage{
asCopula(param)
}
\arguments{
  \item{param}{
The copula is defined by the two-place parameter \code{param}.
}
}

\value{
An \code{\linkS4class{asCopula}} object that behaves (almost) like any other \code{\linkS4class{copula}} from the \code{copula} package.
}
\references{
Example 3.16 in: Nelsen, Roger B. (2006): An Introduction to Copulas, second edition, Springer
}
\author{
Benedikt Graeler
}

\seealso{ 
See also \code{\link{cqsCopula}}.
}
\examples{
persp(asCopula(c(-2,1)),dCopula)
}

\keyword{asymmetric}
\keyword{copula}
\keyword{cubic quadratic sections}
