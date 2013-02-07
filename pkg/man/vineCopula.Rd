\name{vineCopula}
\alias{vineCopula}
\title{
Constructor of the class \code{\linkS4class{vineCopula}}.
}
\description{
Construct an instance of the \code{\linkS4class{vineCopula}} class.
}
\usage{
vineCopula(RVM)
}
\arguments{
  \item{RVM}{
An object of class \code{RVineMatrix} generated from \code{\link{RVineMatrix}} in the package \code{\link{VineCopula-package}}.
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
# an C-vine of independent copulas
RVM <- RVineMatrix(matrix(c(1,0,0,2,2,0,3,3,3),3,byrow=TRUE))

vine <- vineCopula(RVM)
\dontrun{cloud(V1~V2+V3, as.data.frame(rCopula(500,vine)))}
}
\keyword{ mulitvariate }
\keyword{ distribution }
