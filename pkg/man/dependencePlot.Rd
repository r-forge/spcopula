\name{dependencePlot}
\alias{dependencePlot}
\title{
Kernel smoothed scatterplot
}
\description{
Plots a kernel smoothed scatterplot of the provided rank-transformed sample. The work is done by the function \code{\link{panel.smoothScatter}}.
}
\usage{
dependencePlot(var = NULL, smpl, bandwidth = 0.075,  main="Stength of dependece", transformation = function(x) x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var}{
Column IDs or variable names to be used. If not provided, the first two columns will be used.
}
  \item{smpl}{
a matrix (two-columns at least) holding the data
}
  \item{bandwidth}{
the bandwidth passed to the smoothing kernel
}
\item{main}{
the tilte of the plot
}
  \item{transformation}{
a transformation passed to the kernel
}
  \item{\dots}{
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
This is simple wrapper to \code{\link{smoothScatter}}.
}

\seealso{
\code{\link{smoothScatter}}
}
\examples{
dependencePlot(smpl=rCopula(500,asCopula(c(-1,1))))
}
\keyword{plot}
