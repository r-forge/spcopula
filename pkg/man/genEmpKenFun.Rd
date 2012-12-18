\name{genEmpKenFun}
\alias{genEmpKenFun}
\title{
Generates an emprical Kendall distribution function
}
\description{
Generates an emprical Kendall distribution function in cases where no closed form exists. When a sample is provided (\code{smpl}) this exact sample is used to derive the empirical Kendall distribution function. Otherweise, a sample will be generated from the porovided copula.
}
\usage{
genEmpKenFun(copula, sample = NULL)
}
\arguments{
  \item{copula}{
The underlying copula from which 1e6 samples are drawn by default.
}
  \item{sample}{
Optional, a specific sample that should be used to generate the emopircal Kendall distribution function.
}
}
\value{
A function returning the corresponding empirical  Kendalll distribution value.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{getKendallDistr}} for the closed form of the Kendall distribution function.
}
\examples{
empKenFun <- genEmpKenFun(frankCopula(0.7), rCopula(500, frankCopula(0.7)))
empKenFun(.9)
}
\keyword{ multivariate }
\keyword{ distribution }
