\name{stCopula}
\alias{stCopula}
\title{
Spatio-Temporal Copula
}
\description{
Constructor of a bivariate spatio-temporal copula \code{\linkS4class{stCopula}}.
}
\usage{
stCopula(components, distances, t.lags, stDepFun, unit="m", t.res="day")
}

\arguments{
  \item{components}{
A list of bivariate spatial copulas (\code{\linkS4class{spCopula}}) to be used at each temporal lag.
}
  \item{distances}{
A vector of the mean distances of the spatial lag classes.
}
\item{t.lags}{The temporal lags used in the spatio-temporal copula.}
  \item{stDepFun}{
A list of spatial dependence functions; one per temporal lag.
}
  \item{unit}{
The spatial unit, default: m (meters)
}
\item{t.res}{The temporal resolution, default: day}
}
\value{
An instance of the spatio-temporal Copula class \code{\linkS4class{stCopula}}.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{spCopula}}
}
\examples{
data(spCopDemo)

calcKTauPol <- fitCorFun(bins, degree=3)

spCop <- spCopula(components=list(normalCopula(0), tCopula(0, dispstr = "un"),
                                  frankCopula(1), normalCopula(0), claytonCopula(0),
                                  claytonCopula(0), claytonCopula(0), claytonCopula(0),
                                  claytonCopula(0), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")
}
\keyword{spcopula}
\keyword{copula}
