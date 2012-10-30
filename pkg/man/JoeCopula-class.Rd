\name{JoeCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{JoeCopula-class}
\alias{dduCopula,numeric,JoeCopula-method}
\alias{ddvCopula,numeric,JoeCopula-method}
\alias{getKendallDistr,JoeCopula-method}
\alias{kendallDistribution,JoeCopula-method}

\alias{surJoeCopula-class}
\alias{dduCopula,numeric,surJoeCopula-method}
\alias{ddvCopula,numeric,surJoeCopula-method}


\alias{r90JoeCopula-class}
\alias{dduCopula,numeric,r90JoeCopula-method}
\alias{ddvCopula,numeric,r90JoeCopula-method}

\alias{r270JoeCopula-class}
\alias{dduCopula,numeric,r270JoeCopula-method}
\alias{ddvCopula,numeric,r270JoeCopula-method}

\title{Classes \code{"JoeCopula"}, \code{"surJoeCopula"}, \code{"r90JoeCopula"} and \code{"r270JoeCopula"}}
\description{
Wrapper classes representing the Joe, survial Joe, 90 degree and 270 degree rotated Joe copula families (Joe 1997) from \code{\link{CDVine-package}}.
}
\section{Objects from the Classes}{
Objects can be created by calls of the form \code{new("JoeCopula", ...)}, \code{new("surJoeCopula", ...)}, \code{new("r90JoeCopula", ...)} and \code{new("r270JoeCopula", ...)} or by the functions \code{\link{JoeCopula}}, \code{\link{surJoeCopula}}, \code{\link{r90JoeCopula}} and \code{\link{r270JoeCopula}}.
}
\section{Slots}{
  \describe{
    \item{\code{family}:}{Object of class \code{"numeric"} defining the family number in \code{\link{CDVine-package}} }
    \item{\code{dimension}:}{Object of class \code{"integer"} defining the dimension of the copula }
    \item{\code{parameters}:}{Object of class \code{"numeric"} the single parameter}
    \item{\code{param.names}:}{Object of class \code{"character"}, parameter name.}
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"}, lower bound of the copula parameter}
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"}, upper bound of the copula parameter}
    \item{\code{fullname}:}{Object of class \code{"character"}, family name of the copula.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
  \describe{
    \item{dduCopula}{\code{signature(u = "matrix", copula = "JoeCopula")}: ... }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "JoeCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "matrix", copula = "JoeCopula")}: ... }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "JoeCopula")}: ... }
    \item{getKendallDistr}{\code{signature(copula = "JoeCopula")}: ... }
    \item{kendallDistribution}{\code{signature(copula = "JoeCopula")}: ... }
   }
}
\references{
Joe, H., (1997). Multivariate Models and Dependence Concepts. Monogra. Stat. Appl. Probab. 73, London: Chapman and Hall.
}
\author{
Benedikt Graeler
}

\seealso{
See also \code{\linkS4class{BB6Copula}}, \code{\linkS4class{BB7Copula}}, \code{\linkS4class{BB8Copula}} and \code{\linkS4class{JoeCopula}} for further wrapper classes to the \code{\link{CDVine-package}}.
}
\examples{
showClass("JoeCopula")
}
\keyword{classes}