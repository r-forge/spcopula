\name{genEmpCop}
\alias{genEmpCop}
\title{
Generate an empirical copula
}
\description{
Generates an empirical copula from a sample and returns the corresponding function.
}
\usage{
genEmpCop(data)
}
\arguments{
  \item{data}{
The sample to be used as input for the empirical copula.
}
}
\value{
The empirical copula as a function (the multivariate cdf).
}
\author{
Benedikt Graeler
}

\examples{
empCop <- genEmpCop(rCopula(500, frankCopula(0.7)))

# the empirical value
empCop(c(0.3, 0.5))

# the theoretical value
pCopula(c(0.3, 0.5), frankCopula(0.7))
}
\keyword{ multivariate }