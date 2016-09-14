\name{mixtureCopula}
\alias{mixtureCopula}

\title{
Constructor of a mixture copula}
\description{
An instance of a convex mixture of two copulas is generated.
}
\usage{
mixtureCopula(param = c(0.2, 0.2, 0.5), memberCops = c(normalCopula(0), claytonCopula(1)))
}

\arguments{
  \item{param}{A numeric vector: the set of parameters for the first and second copula appended by the mixture coefficient.}
  \item{memberCops}{A list of two copulas.}
}
\value{
An object of class \linkS4class{mixtureCopula}}
\author{
Benedikt Graeler}

\examples{

mixCop <- mixtureCopula()

persp(mixCop, dCopula)
}
