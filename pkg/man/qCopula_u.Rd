\name{qCopula_u}
\alias{qCopula_u}
\alias{qCopula_u,copula-method}

\title{
The inverse of a bivariate copula given u
}
\description{
The inverse of a bivariate copula is calculated for a given u.
}
\usage{
qCopula_u(copula, p, u, ...)
}
\arguments{
  \item{copula}{
the copula of interest
}
  \item{p}{
the probability level p
}
  \item{u}{
the conditioning variable u
}
  \item{\dots}{
Passed on to \code{\link{optim}} in the background.
}
}
\details{
The evaluation is done numerically using either \code{\link{optim}} or \code{\link{optimise}}.
}
\value{
A matrix having the same number of rows as \"u\" providing u and the other arguments.
}
\author{
Benedikt Graeler
}

\examples{
uv <- qCopula_u(asCopula(c(-1,1)), p=rep(0.9,10), u=runif(10,0.9,1))
pCopula(uv,asCopula(c(-1,1)))-0.9
}

\keyword{ distribution }
\keyword{ multivariate }
