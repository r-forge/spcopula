\name{dependencePlot}
\alias{dependencePlot}
\title{
Kernel smoothed scatterplot
}
\description{
Plots a kernel smoothed scatterplot of the provided rank-transformed sample. The work is done by the function \code{\link{panel.smoothScatter}}.
}
\usage{
dependencePlot(formula = NULL, smpl, cuts = 15, bandwidth = 0.075, transformation = function(x) x, ..., range.x = list(c(0, 1), c(0, 1)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
the formula as in common lattice plotting
}
  \item{smpl}{
can be a two-column matrix holding the data, no formula is needed
}
  \item{cuts}{
the cuts to change colors
}
  \item{bandwidth}{
the bandwidth passed to the smoothing kernel
}
  \item{transformation}{
a transformation passed to the kernel
}
  \item{\dots}{
passed on to the function \code{\link{panel.smoothScatter}}
}
  \item{range.x}{
passed on to the function \code{\link{panel.smoothScatter}}
}
}
\details{
see \code{\link{panel.smoothScatter}}
}
\value{
A plot.
}
\author{
Benedikt Graeler
}
\note{
This is simple wrapper to \code{\link{panel.smoothScatter}}.
}

\seealso{
\code{\link{panel.smoothScatter}}
}
\examples{
dependencePlot(smpl=rCopula(500,asCopula(c(-1,1))))
}
\keyword{plot}
