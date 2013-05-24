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
    \item{dataLocs}{The locations holding the data for the prediction. May be left empty (NULL) for estimation of the spatial (vine) copula.}
    \item{index}{Object of class \code{matrix} pointing to the IDs of the underlying \code{\link[sp:Spatial-class]{Spatial}} object.}
    \item{locations}{Object of class \code{\link[sp:Spatial-class]{Spatial}} denoting the central locations of the neighbourhoods.}
    \item{var}{Object of class \code{character} giving the variable names of interest.}
    \item{prediction}{of class \code{logical}, whether this neighbourhood is suited for fitting or prediction purposes.}
    \item{bbox}{Object of class \code{matrix} the bounding box of the central neighbourhood locations.}
    \item{proj4string}{Object of class \code{CRS} the coordinate reference system of the underlying \code{\link[sp:Spatial-class]{Spatial}} object.}
  }
}
\section{Extends}{
Class \code{\link[sp:Spatial-class]{Spatial}}, directly.
}
\section{Methods}{
  \describe{
    \item{names}{\code{signature(x = "neighbourhood")}: provides the variable names of the neighbourhood. }
    \item{show}{\code{signature(object = "neighbourhood")}: a brief description of the characteristics of the neighbourhood.}
    \item{spplot}{\code{signature(obj = "neighbourhood")}: plots the values of the one or more selected columns (0="central location", 1="closest neighbours", ...) at the coordinates of the central location.}
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
spdf <- data.frame(x=c(112,154,212,289),y=c(124,198,85,168),measure=rlnorm(4))
coordinates(spdf) <- ~x+y

neigh <- getNeighbours(spdf,size=4)
spplot(neigh, column=0:3)
}
\keyword{classes}