\name{neighbourhood}
\alias{neighbourhood}
\title{
Constructor of the \code{\linkS4class{neighbourhood}} class.
}
\description{
An object that stores the nearest k-neighbours for each location.
}
\usage{
neighbourhood(data, distances, sp, index)
}
\arguments{
  \item{data}{
A data.frame that contains the variables of interest for each local neighbourhood.
}
  \item{distances}{
A matrix with all the separating distances per local neighbourhood.
}
  \item{sp}{
The coordiantes of the original spatial data.frame.
}
  \item{index}{
An index pointing from each neighbourhood to the corresponding spatial IDs.
}
}
\value{
A \code{\linkS4class{neighbourhood}} object.
}

\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{neighbourhood}}, \code{\link{getNeighbours}}
}
\examples{
spdf <- data.frame(x=c(112,154,212,289),y=c(124,198,85,168),measure=rlnorm(4))
coordinates(spdf) <- ~x+y

neigh <- getNeighbours(spdf,size=4)

colnames(neigh@data) <- rep("measure",3)

# rebuilding neigh
neighbourhood(data=neigh@data, neigh@distances, SpatialPoints(neigh@coords), neigh@index)
neigh
}

\keyword{ spatial }