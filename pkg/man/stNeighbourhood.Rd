\name{stNeighbourhood}
\alias{stNeighbourhood}
\title{
Constructor of the \code{\linkS4class{stNeighbourhood}} class.
}
\description{An object that stores the nearest k spatio-temporal neighbours for a set of locations.}

\usage{stNeighbourhood(data, distances, STxDF, ST=NULL, index, prediction, var)}

\arguments{
  \item{data}{A data.frame that contains the variables of interest for each local spatio-temporal neighbourhood.}
  \item{distances}{A matrix with all the separating distances per local spatio-temporal neighbourhood.}
  \item{STxDF}{A \code{\link[sp:Spatial-class]{Spatial}} object describing the central locations the neighbourhood; must be of same length as rows in \code{data} and \code{distances}.}
  \item{ST}{The locations used for the prediction. May be left empty (NULL) for estimation of the spatial (vine) copula.}
  \item{index}{An index pointing from each neighbourhood to the corresponding spatial IDs.}
  \item{prediction}{whether this neighbourhood is intended to be used for prediction or fitting purpose.}
  \item{var}{The variable name(s).}
}

\value{A \code{\linkS4class{stNeighbourhood}} object.}

\author{Benedikt Graeler}

\seealso{
\code{\linkS4class{stNeighbourhood}}, \code{\link{getStNeighbours}}
}
\examples{
library(sp)
library(spacetime)

sp <- SpatialPoints(matrix(c(181000,181100,333500,333600),2))
time <- Sys.time()+60*60*24*c(0,1,2)
data <- data.frame(var1=runif(6))

stData <- STFDF(sp, time, data)
stQuerry <- STF(SpatialPoints(matrix(c(181000,181200,333600,333600),2)),
                time[2:3])

getStNeighbours(stData=stData, ST=stQuerry, prediction=TRUE, spSize=3,
                t.lags=-(0:1))
}

\keyword{spatio-temporal}