\name{neighbourhood}
\alias{neighbourhood}
\title{
Constructor of the \code{\linkS4class{neighbourhood}} class.
}
\description{An object that stores the nearest k-neighbours for a set of locations.}

\usage{neighbourhood(data, distances, index, dataLocs, predLocs=NULL, prediction, var)}

\arguments{
  \item{data}{A data.frame that contains the variables of interest for each local neighbourhood.}
  \item{distances}{A matrix with all the separating distances per local neighbourhood.}
  \item{index}{An index pointing from each neighbourhood to the corresponding spatial IDs.}
  \item{dataLocs}{A \code{\link[sp:Spatial-class]{Spatial}} object describing the central locations the neighbourhood; must be of same length as rows in \code{data} and \code{distances}.}
  \item{predLocs}{The locations used for the prediction. May be left empty (NULL) for estimation of the spatial (vine) copula.}
  \item{prediction}{whether this neighbourhood is intended to be used for prediction or fitting purpose.}
  \item{var}{The variable name(s).}
}

\value{A \code{\linkS4class{neighbourhood}} object.}

\author{Benedikt Graeler}

\seealso{
\code{\linkS4class{neighbourhood}}, \code{\link{getNeighbours}}
}
\examples{
spdf <- data.frame(x=c(112,154,212,289),y=c(124,198,85,168),measure=rlnorm(4))
coordinates(spdf) <- ~x+y

neigh <- getNeighbours(spdf,size=4)
neigh

# rebuilding neigh
neighbourhood(neigh@data, neigh@distances, neigh@index, spdf, NULL, neigh@prediction, neigh@var)
}

\keyword{spatial}