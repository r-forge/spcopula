\name{vineCopula}
\alias{vineCopula}
\title{
Constructor of the class \code{\linkS4class{vineCopula}}.
}
\description{
Cosntruct an instance of the \code{\linkS4class{vineCopula}} class.
}
\usage{
vineCopula(copulas, dim, type)
}
\arguments{
  \item{copulas}{
A list of copulas that build up the vine.
}
  \item{dim}{
The dimension of the vine copula.
}
  \item{type}{
The type of the vine, either \"c-vine\" or \"d-vine\".
}
}
\value{
An instance of the \code{\linkS4class{vineCopula}} class.
}
\references{
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009). Pair-copula constructions of multiple dependence Insurance: Mathematics and Economics 44 (2), 182-198. 
}
\author{
Benedikt Graeler
}
\examples{
vine <- vineCopula(list(frankCopula(.7), gumbelCopula(3), gumbelCopula(1)),
                   3, "c-vine")
\dontrun{cloud(V1~V2+V3, as.data.frame(rCopula(500,vine)))}
}
\keyword{ mulitvariate }
\keyword{ distribution }
