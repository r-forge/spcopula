\name{stCopPredict}
\alias{stCopPredict}

\title{
spatio-temporal prediction based on a spatio-temporal vine copula
}
\description{
A spatio-temporal vine copula is used to predict values at unobserved spatio-temporal locations conditioned on observations of a local spatio-temporal neighbourhood.
}
\usage{
stCopPredict(predNeigh, stVine, margin, method = "quantile", p = 0.5, ...)
}
\arguments{
  \item{predNeigh}{the \code{\linkS4class{neighbourhood}} used for prediction, its slot \code{prediction} must be \code{TRUE}.}
  \item{stVine}{the spatio-temporal vine copula describing the spatio-temporal dependence}
  \item{margin}{the marginal distribution as a list with entries named "d" for the density function (PDF), "q" for the quantile function and "p" for cumulative distribution function (CDF).}
  \item{method}{one of \code{"quantile"} or \code{"expectation"} denoting the type of predictor.}
  \item{p}{only used for the quantile predictor indicating the desired fraction the quantile should correspond to.}
  \item{...}{Further arguments passed to \code{\link{integrate}}.}
}
\details{
Predictions are done based on \code{\link{condStVine}} through numerical integration/optimisation.
}
\value{
A \code{\link[spacetime:ST-class]{ST}} object of the same type as provided in the slot \code{locations} of the argument \code{predNeigh}.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{condStVine}} and \code{\link{spCopPredict}} for the spatial version.
}
\examples{
## the spatial version
data(spCopDemo)

calcKTauPol <- fitCorFun(bins, degree=3)

spCop <- spCopula(components=list(normalCopula(0), tCopula(0, dispstr = "un"),
                                  frankCopula(1), normalCopula(0), claytonCopula(0),
                                  claytonCopula(0), claytonCopula(0), claytonCopula(0),
                                  claytonCopula(0), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")

spVineCop <- spVineCopula(spCop, vineCopula(4L))

data(meuse.grid)
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- TRUE

data(meuse)
coordinates(meuse) <- ~x+y

predMeuseNeigh <- getNeighbours(meuse[1:4,], meuse.grid[c(9:12,15:19,24:28,34:38),],"zinc",5L,TRUE,-1)
predMeuseNeigh@data <- rankTransform(predMeuseNeigh@data)

qMar <- function(x) {
  qlnorm(x,mean(log(meuse$zinc)),sd(log(meuse$zinc)))
}

predMedian <- spCopPredict(predMeuseNeigh, spVineCop, list(q=qMar), "quantile", p=0.5)

spplot(predMedian,"quantile.0.5", 
       sp.layout=list("sp.points", meuse, pch = 19, col = "red"),
       col.regions=bpy.colors())
}

\keyword{ distribution }
\keyword{ prediction }
\keyword{ spatial }