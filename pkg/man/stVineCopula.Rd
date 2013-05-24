\name{stVineCopula}
\alias{stVineCopula}
\title{
Constructor of the class \code{\linkS4class{stVineCopula}}.
}
\description{
Construct an instance of the \code{\linkS4class{stVineCopula}}.
}
\usage{
stVineCopula(stCop, topCop)
}
\arguments{
  \item{stCop}{Object of class \code{\linkS4class{stCopula}}. A bivariate spatio-temporal copula. }
  \item{topCop}{Object of class \code{\linkS4class{copula}}. A bivariate, multivariate or R-Vine copula joining the spatio-temporal components (neighbours). In case all trees are spatio-temporal copulas, this argument might be \code{NULL}.}
}
\value{
An instance of the \code{\linkS4class{stVineCopula}} class.
}
\author{
Benedikt Graeler
}
\examples{
## the spatial version
# a spatial C-vine copula (with independent dummy copulas in the upper vine)
data(spCopDemo)

calcKTauPol <- fitCorFun(bins, degree=3)

spCop <- spCopula(components=list(normalCopula(0), tCopula(0, dispstr = "un"),
                                  frankCopula(1), normalCopula(0), claytonCopula(0),
                                  claytonCopula(0), claytonCopula(0), claytonCopula(0),
                                  claytonCopula(0), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")

RVM <- RVineMatrix(matrix(c(1,0,0,2,2,0,3,3,3),3,byrow=TRUE))
spVineCopula(spCop, vineCopula(RVM))
}
\keyword{ mulitvariate }
\keyword{ distribution }