\name{getStNeighbours}
\alias{getStNeighbours}

\title{
Creating Local Spatio-Temporal Neighbourhoods
}
\description{
This function calculates local spatio-temporal neighbourhoods to be used for fitting of spatio-temporal vine copulas and for prediction using spatio-temporal vine copulas.
}
\usage{
getStNeighbours(stData, ST, var = names(stData@data)[1], spSize = 4, t.lags=-(0:2),
                timeSteps=NA, prediction=FALSE, min.dist = 0.01)
}
\arguments{
  \item{stData}{
some spatio-temporal data frame holding the data used for estimation/prediction
}
  \item{ST}{
A spatio-temporal object defining the prediction locations, might be missing if the spatio-temporal neighbourhood is used for fitting.
}
  \item{var}{
the variable name of interest, by default the first variable is used.
}
  \item{spSize}{
The spatial size of the neighbourhood including the location of interest (for fitting as well for prediction).
}
\item{t.lags}{The temporal lags to be used in the spatio-temporal neighbourhood.}
\item{timeSteps}{The number of tiem instances that should randomly be selected from \code{stData}. The defualt, \code{NA}, selects all locations.}
  \item{prediction}{whether the neighbourhood should be used for prediction (TRUE) or spatial/Spatio-temporal vine copula fitting.}
  \item{min.dist}{
the minimal distance for a location to be included. Must be larger than 0 for fitting purposes and might be 0 for prediction.
}
}
\value{
An object of \code{\linkS4class{stNeighbourhood}}.
}
\author{
Benedikt Graeler
}

\seealso{
See \code{\link{stNeighbourhood}} for the native constructor of a \code{\linkS4class{stNeighbourhood}} class. The pure spatial version can be found at \code{\link{getNeighbours}}.
}
\examples{
## the spatial version:
spdf <- data.frame(x=c(112,154,212,289,345),y=c(124,198,85,168,346),measure=rlnorm(5))
coordinates(spdf) <- ~x+y

getNeighbours(spdf,size=4)
}
\keyword{ spatial }