\name{loglikByCopulasLags}
\alias{loglikByCopulasLags}
\title{
Log-Lieklihood per lag calss
}
\description{
The functio calculates the log-lieklihood for a set of provided copula families per lag class.
}
\usage{
loglikByCopulasLags(bins, calcTau, families = c(normalCopula(0), tCopula(0, dispstr = "un"), claytonCopula(0), frankCopula(1), gumbelCopula(1)))
}

\arguments{
  \item{bins}{
Typically, the output of the \code{\link{calcBins}} function providing columns \code{meanDists} to tune the copula's parameter and \code{lagData} to evaluate the log-likelihood.
}
  \item{calcTau}{
The tuning function for the copulas' parameters. Typically, the output of the \code{\link{fitCorFun}} function. 
}
  \item{families}{
A vector of representatives of the bivariate copula families that should be investigated.
}
}
\value{
A matrix of spatial lags (rows) and copula family names (columns) holding the calculated log-likelihood values.
}
\author{
Benedikt Graeler
}

\seealso{
See also \code{\link{calcBins}}, \code{\link{fitCorFun}} or \code{\link{spCopula}}.
}

\examples{
data(spCopDemo) # load the workspace from demo(spcopula_estimation)

loglikTau <- loglikByCopulasLags(bins, calcKTauPol)
loglikTau
}

\keyword{spcopula}
\keyword{log-likelihood}