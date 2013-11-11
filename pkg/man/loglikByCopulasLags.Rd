\name{loglikByCopulasLags}
\alias{loglikByCopulasLags}
\title{
Log-likelihood per lag class
}
\description{
This function calculates the log-likelihood for a set of provided copula families per lag class. The copulas' parameters are either fitted by a provided distance dependent function, or through the function \code{\link{fitCopula}} per lag and copula family. 
}
\usage{
loglikByCopulasLags(bins, families = c(normalCopula(0), tCopula(0, dispstr = "un"),
                    claytonCopula(0), frankCopula(1), gumbelCopula(1)), calcCor)
}

\arguments{
  \item{bins}{
Typically, the output of the \code{\link{calcBins}} function providing columns \code{meanDists} to tune the copula's parameter and \code{lagData} to evaluate the log-likelihood.
}
  \item{families}{
A vector of representatives of the bivariate copula families that should be investigated.
}
  \item{calcCor}{
The tuning function for the copulas' parameters. Typically, the output of the \code{\link{fitCorFun}} function. Its either based on Kendall's tau, Spearman's rho or the identity. This function needs to return one of "kendall", "spearman" or "id" when passed "NULL". This parameter might be "missing" for static convex spatial copulas.
}
}
\value{
A list containing a matrix (\code{loglik}) of spatial lags (rows) and copula family names (columns) holding the calculated log-likelihood value and a list o the corresponding copula fits.
}
\author{
Benedikt Graeler
}

\seealso{
See also \code{\link{calcBins}}, \code{\link{fitCorFun}} or \code{\link{spCopula}}.
}

\examples{
data(spCopDemo) # load the workspace from demo(spcopula_estimation)

calcKTauPol <- fitCorFun(bins, degree=3)

loglikTau <- loglikByCopulasLags(bins, calcCor=calcKTauPol)
loglikTau$loglik
}

\keyword{spcopula}
\keyword{log-likelihood}