\name{condSpVine}
\alias{condSpVine}

\title{
Conditions a spatial vine copula for conditional prediction
}
\description{
A spatial vine copula is conditioned under the observations of all but one neighbour generating a conditional uivariate distribution used ofr prediction.
}
\usage{
condSpVine(condVar, dists, spVine, n = 100)
}
\arguments{
  \item{condVar}{
the conditional variables
}
  \item{dists}{
the distances between the neighbours to adjust the spatial vine copula parameters.
}
  \item{spVine}{
the spatial vine copula
}
  \item{n}{
a proxy to the number of intervals to be used in the numerical integration. The intervals will not be split uniform with a shorter intervals at the corners of the copula.
}
}
\value{
A function describing the conditional density.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{spVineCopula}}
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

spVineCop <- spVineCopula(spCop, vineCopula(4L))

dists <- runif(4,0,800)
condVar <- runif(4)
condDensity <- condSpVine(condVar,dists,spVineCop)

curve(condDensity)
mtext(paste("Dists:",paste(round(dists,0),collapse=", ")),line=0)
mtext(paste("Cond.:",paste(round(condVar,2),collapse=", ")),line=1)
}

\keyword{ distribution }
