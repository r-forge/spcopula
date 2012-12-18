\name{kendallRP}
\alias{kendallRP}
\title{
calculating the Kendall Return Period
}
\description{
The Kendall Return Period corresponding to a given critical level and a given 
kendall distribution or its underlying copula.
}
\usage{
kendallRP(kendallFun = NULL, cl = c(0.99, 0.999), mu = 1, copula = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kendallFun}{
The kendall distribution function.
}
  \item{cl}{
The desired critical level(s): the default are 0.99 and 0.999.
}
  \item{mu}{
The reoccurence of an event per time unit; the default is 1.
}
  \item{copula}{
The copula describing the dependencies between the investigated variables. Only necessary if \code{kendallFun} is not provided.
}
}
\value{
A vector of the same length as \code{cl} with the corresponding endall Return Perios.}
\references{
Salvadori, G., De Michele, C., and Durante, F.: On the return period and design in a multivariate framework, Hydrol. Earth Syst. Sci., 15, 3293-3305, doi:10.5194/hess-15-3293-2011, 2011.
}
\author{
Benedikt Graeler
}
\seealso{
\code{\link{criticalLevel}} for the reverse direction: what is the critical level for a given Kendall Return Period?
}
\examples{
kendallRP(getKendallDistr(frankCopula(7)), cl=c(0.9,0.99,0.999))
}
\keyword{ survival }
\keyword{ multivariate }