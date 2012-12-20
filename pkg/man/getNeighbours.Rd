\name{getNeighbours}
\alias{getNeighbours}

\title{
Creating Local Neighbourhoods
}
\description{
This function calculates a local neighbourhood to be used in the spatial and spatio-temporal vine copulas.
}
\usage{
getNeighbours(spData, var = names(spData), size = 4, dep = NULL, indep = NULL, min.dist = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spData}{
some spatial data frame
}
  \item{var}{
the variable names of interest
}
  \item{size}{
The size of the neighbourhood including the location of interest.
}
  \item{dep}{
Index of dependent locations.
}
  \item{indep}{
Index of independent locations.
}
  \item{min.dist}{
the minimal distance for a location to be included.
}
}
\value{
An object of \code{\linkS4class{neighbourhood}}.
}
\author{
Benedikt Graeler
}

\seealso{
See \code{\link{neighbourhood}} for the native constructor of a \code{\linkS4class{neighbourhood}} class.
}
\examples{
spdf <- data.frame(x=c(112,154,212,289,345),y=c(124,198,85,168,346),measure=rlnorm(5))
coordinates(spdf) <- ~x+y

neigh <- getNeighbours(spdf,size=4)

neigh
}
\keyword{ spatial }