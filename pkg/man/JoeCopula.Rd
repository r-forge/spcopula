\name{JoeCopula}
\alias{JoeCopula}
\alias{surJoeCopula}
\alias{r90JoeCopula}
\alias{r270JoeCopula}

\title{
Constructor of the Joe family and its derivates
}
\description{
Constructs an object of the \code{\linkS4class{JoeCopula}} (survival \code{sur}, 90 degree rotated \code{r90} and 270 degree rotated \code{r270}) family for a given parameter.
}
\usage{
JoeCopula(param)
surJoeCopula(param)
r90JoeCopula(param)
r270JoeCopula(param)
}

\arguments{
  \item{param}{
The parameter \code{param} defines the copula through \code{theta} and \code{delta}.
}
}
\value{
One of the respective Joe copula classes (\code{\linkS4class{JoeCopula}}, \code{\linkS4class{surJoeCopula}}, \code{\linkS4class{r90JoeCopula}}, \code{\linkS4class{r270JoeCopula}}).
}
\references{
Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
}
\author{
Benedikt Graeler
}
\seealso{
See also \code{\link{BB6Copula}}, \code{\link{BB7Copula}}, \code{\link{BB8Copula}} and \code{\link{JoeCopula}} for further wrapper functions to the \code{\link{CDVine-package}}. 
}
\examples{
persp(JoeCopula(1.5),dCopula, zlim=c(0,10))
persp(surJoeCopula(1.5),dCopula, zlim=c(0,10))
persp(r90JoeCopula(-1.5),dCopula, zlim=c(0,10))
persp(r270JoeCopula(-1.5),dCopula, zlim=c(0,10))
}