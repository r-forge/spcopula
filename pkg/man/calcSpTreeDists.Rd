\name{calcSpTreeDists}
\alias{calcSpTreeDists}

\title{
Calculating the conditional spatial distances
}
\description{
All spatail distances of all involved neighbours are calculated according to multiple spatial copulas forming a C-Vine.
}
\usage{calcSpTreeDists(neigh, n.trees)}

\arguments{
  \item{neigh}{
the initial \code{\linkS4class{neighbourhood}}
}
  \item{n.trees}{
the number of spatial trees
}
}

\value{
A \code{list} of the spatial distance of the pairs involved in the conditinal data pairs.
}

\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{neighbourhood}}
}

\keyword{ spatial }
\keyword{ multivariate }
\keyword{distribution}
