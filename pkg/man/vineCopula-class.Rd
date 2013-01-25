\name{vineCopula-class}
\Rdversion{1.1}
\docType{class}
\alias{vineCopula-class}

\title{Class \code{"vineCopula"}}
\description{
A class representing vine copulas in a object oriented implementations. Many functions go back to the package \code{\link{CDVine-package}}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("vineCopula", ...)} or through the function \code{\link{vineCopula}}.
}
\section{Slots}{
  \describe{
    \item{\code{copulas}:}{Object of class \code{"list"} holding all copulas. }
    \item{\code{dimension}:}{Object of class \code{"integer"}; the vines dimension. }
    \item{\code{type}:}{Object of class \code{"character"} defining the type of the vine (either \"c-vine\" or \"d-vine\").}
    \item{\code{parameters}:}{Object of class \code{"numeric"}: empty }
    \item{\code{param.names}:}{Object of class \code{"character"}: empty }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"}: empty }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"}: empty }
    \item{\code{fullname}:}{Object of class \code{"character"} providing a descriptive name of the vine copula.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
No additional methods yet, but uses e.g. \code{\link{dCopula}}, \code{\link{pCopula}}, \code{\link{rCopula}} as any other copula.
}
\references{
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009). Pair-copula constructions of multiple dependence Insurance: Mathematics and Economics 44 (2), 182-198. 
}
\author{
Benedikt Graeler
}

\examples{
showClass("vineCopula")
}
\keyword{classes}
