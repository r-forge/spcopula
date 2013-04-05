\name{copulaFromFamilyIndex}
\alias{copulaFromFamilyIndex}
\title{
Construct a copula object from a \code{\link{VineCopula-package}} family index
}
\description{
A \code{\linkS4class{copula}} object is constructed from the family index used in the package \code{\link{VineCopula-package}} for the provided parameters.
}
\usage{
copulaFromFamilyIndex (family, par, par2 = 0) 
}
\arguments{
  \item{family}{
The number identifying the desired copula family.
}
  \item{par}{
the first parameter as used in \code{\link{VineCopula-package}}
}
  \item{par2}{
the second parameter as used in \code{\link{VineCopula-package}}. The default is \code{par2=0} for single parameter families. 
}
}
\value{
A \code{\linkS4class{copula}} object of the desired family with the provided parameters. 
}

\author{
Benedikt Graeler
}
\examples{
# the survival Joe Copula
cop <- copulaFromFamilyIndex(16, 3)
class(cop)
}
\keyword{ function}