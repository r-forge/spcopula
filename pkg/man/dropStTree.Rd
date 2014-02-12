\name{dropStTree}
\alias{dropStTree}

\title{
Drops a spatio-temporal tree
}
\description{
A spatio-temporal tree is dropped in order to fit the following copula.
}
\usage{dropStTree(neigh, stCop)}

\arguments{
  \item{neigh}{
the current spatio-temporal \code{\linkS4class{stNeighbourhood}}
}
  \item{stCop}{
the current spatio-temporal copula performing the conditioning
}
}

\value{
A conditioned spatio-temporal \code{\linkS4class{stNeighbourhood}} of dimension 1 less.
}

\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{stNeighbourhood}}
}

\keyword{ spatio-temporal }