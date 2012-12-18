\name{unitScatter}
\alias{unitScatter}
\alias{univScatter}
\title{
A scatterplot on the unit-sqaure
}
\description{
Plots a sample on the unit-square. If needed, it will transform the sample using the rank order transformation as in \code{\link{rankTransform}}.
}
\usage{
unitScatter(formula = NULL, smpl, ...)
}
\arguments{
  \item{formula}{
A standard plotting formula to select the corresponding columns.
}
  \item{smpl}{
The data set to be used.
}
\item{...}{
Passed on to \code{\link{xyplot}}.
}
}
\value{
A lattice plot.
}
\author{
Benedikt Graeler
}
\examples{
data(loss)
rt_loss <- rankTransform(loss[,1:2])
unitScatter(loss~alae,rt_loss)
}
\keyword{ hplot }

