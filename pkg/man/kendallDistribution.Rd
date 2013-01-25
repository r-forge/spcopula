\name{kendallDistribution}
\alias{kendallDistribution}
\alias{kendallDistribution,claytonCopula-method}
\alias{kendallDistribution,frankCopula-method}
\alias{kendallDistribution,gumbelCopula-method}
\title{
The Kendall distribution
}
\description{
The Kendall distribution for a given copula.
}
\usage{
kendallDistribution(copula, t)
}
\arguments{
  \item{copula}{
A copula of which the Kendall Distribution should be calculated.
}
  \item{t}{
The quantile of the copula for which the Kendall distri bution should be calculated.
}
}
\value{
A vector of Kendall distribution values the same length as \code{t}.
}

\author{
Benedikt Graeler
}

\seealso{
\code{\link{getKendallDistr}}}
\examples{
kendallDistribution(frankCopula(.8),.8)
}