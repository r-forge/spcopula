\name{rankTransform}
\alias{rankTransform}
\title{
rank order transformation of margins
}
\description{
performs the rank order transfroamtion
}
\usage{
rankTransform(u, v = NULL, ties.method = "average")
}
\arguments{
  \item{u}{
a matrix or data.frame with at least two columns holding the data or a vector holding the first column of the data
}
  \item{v}{
a vector holding the second column of the data
}
  \item{ties.method}{
How should ties be treated in \code{\link{rank}}?
}
}
\value{
A matrix or data.frame (as provided) with the transformed data. Rows containing any \code{NA} will be dropped.
}

\author{
Benedikt Graeler
}
\examples{
data(simulatedTriples)
str(rankTransform(triples))
}
\keyword{ function}