\name{dropSpTree}
\alias{dropSpTree}

\title{
Drops a spatial tree
}
\description{
A spatial tree is dropped in order to fit the following spatial copula.
}
\usage{dropSpTree(neigh, spCop)}

\arguments{
  \item{neigh}{
the current \code{\linkS4class{neighbourhood}}
}
  \item{spCop}{
the current spatial copula performing the conditioning
}
}

\value{
A conditioned \code{\linkS4class{neighbourhood}} of dimesnion 1 less.
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
