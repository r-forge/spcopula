\name{neighbourhood-class}
\Rdversion{1.1}
\docType{class}
\alias{neighbourhood-class}
\alias{names,neighbourhood-method}
\alias{[[,neighbourhood,ANY,ANY-method}
\alias{[[,neighbourhood,numeric,missing-method}

\title{Class \code{neighbourhood}}
\description{A class representing a local spatial neighbourhood.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("neighbourhood", ...)}, from the constructor \code{\link{neighbourhood}} or through the easier function \code{\link{getNeighbours}}.
}
\section{Slots}{
  \describe{
    \item{data}{Object of class \code{data.frame} holding the data for each local neighbourhood.}
    \item{distances}{Object of class \code{matrix} providing the separating distances per neighbourhood.}
    \item{index}{Object of class \code{matrix} pointing to the IDs of the underlying \code{\link[sp:Spatial-class]{Spatial}} object.}
    \item{var}{Object of class \code{character} giving the variable's name.}
    \item{coVar}{Object of class \code{character} giving the covariate's name.}
    \item{prediction}{of class \code{logical}, whether this neighbourhood is suited for fitting or prediction purposes.}
  }
}
\section{Methods}{
  \describe{
    \item{names}{\code{signature(x = "neighbourhood")}: provides the variable names of the neighbourhood. }
    \item{show}{\code{signature(object = "neighbourhood")}: a brief description of the characteristics of the neighbourhood.}
    \item{[[}{subsets the selection of neighbourhoods returning a subset of these ("column wise")}
	 }
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{getNeighbours}}
}
\examples{
library(sp)
spdf <- data.frame(x=c(112,154,212,289),y=c(124,198,85,168),measure=rlnorm(4))
coordinates(spdf) <- ~x+y

getNeighbours(spdf,size=4)
}
\keyword{classes}