\name{spCopula}
\alias{spCopula}
\title{
Spatial Copula
}
\description{
Constructor of a bivariate spatial copula \code{\linkS4class{spCopula}}.
}
\usage{
spCopula(components, distances, spDepFun, unit = "m")
}

\arguments{
  \item{components}{
A list of bivariate copulas to be used at each lag.Will be supplemented by the independeCopula.
}
  \item{distances}{
A vecotr of the mean distances of the lag classes PLUS the range which will define the use of the independece copula.
}
  \item{spDepFun}{
The spatial dependece function. Typically, the result of \code{\link{fitCorFun}}. If no spatial dependence function is provided, the spatial copula will be a pure linear convex combination of the provided components (i.e. no further fit of their parameters will take place).
}
  \item{unit}{
The spatial unit, defualt: m (meters)
}
}
\value{
An instance of the spatial Copula class.
}
\references{
Graeler, B. and E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. Procedia Environmental Sciences (Vol. 7, pp. 206 - 211), Elsevier. 
}
\author{
Benedikt Graeler
}

\seealso{
See also \code{\link{calcBins}} or \code{\link{loglikByCopulasLags}} and \code{demo(spcopula_estimation)}.
}
\examples{
data(spCopDemo)

spCop <- spCopula(components=list(normalCopula(0), tCopula(0, dispstr = "un"),
                                  frankCopula(1), normalCopula(0), claytonCopula(0),
                                  claytonCopula(0), claytonCopula(0), claytonCopula(0),
                                  claytonCopula(0), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")
}
\keyword{spcopula}
\keyword{copula}
