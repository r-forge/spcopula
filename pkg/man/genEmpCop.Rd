\name{genEmpCop}
\alias{genEmpCop}
\title{
Generate an empirical copula
}
\description{
Generates an empirical copula from a sample and returns the corresponding function.
}
\usage{
genEmpCop(copula, sample.size=1e5)
}
\arguments{
  \item{copula}{
The theoretical copula to be represnted numerically.
}
\item{sample.size}{
The length of the sample to be used. The default is 1e5.
}
}
\value{
An empirical copula class \code{\linkS4class{empiricalCopula}}.
}
\author{
Benedikt Graeler
}

\examples{
empCop <- genEmpCop(frankCopula(0.7), 500)

# the empirical value
pCopula(c(0.3, 0.5), empCop)

# the theoretical value
pCopula(c(0.3, 0.5), frankCopula(0.7))
}
\keyword{ multivariate }