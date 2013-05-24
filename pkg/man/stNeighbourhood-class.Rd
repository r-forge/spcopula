\name{stNeighbourhood-class}
\Rdversion{1.1}
\docType{class}
\alias{stNeighbourhood-class}

\title{Class \code{"stNeighbourhood"}}
\description{
An object representing a set of spatio-temporal neighbourhoods including data, spatio-temporal indices and spatio-temporal distances.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("stNeighbourhood", ...)} or through the simplified call to \code{\link{getStNeighbours}}.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"data.frame"} holding the data of spatio-temporal neighbourhoods. }
    \item{\code{distances}:}{Object of class \code{"array"} with the following three dimensions: number of neighbourhoods, dimension of each neighbourhood, 2 (spatial and temporal distance).}
    \item{\code{index}:}{Object of class \code{"array"} with the following three dimensions: number of neighbourhoods, dimension of each neighbourhood, 2 (spatial and temporal index). }
    \item{\code{locations}:}{Object of class \code{"ST"} the conditioning locations.}
    \item{\code{dataLocs}:}{Object of class \code{"optionalST"} the prediction locations, may be of class \code{"ST"} or \code{"NULL"}. }
    \item{\code{var}:}{Object of class \code{"character"} the name of the variable that is estimated.}
    \item{\code{prediction}:}{Object of class \code{"logical"} whether prediction or fitting takes place. }
    \item{\code{sp}:}{Object of class \code{"Spatial"} the spatial extend represented as a bounding box with CRS.}
    \item{\code{time}:}{Object of class \code{"xts"} the earliest time point of all spatio-temporal neighbourhoods. }
    \item{\code{endTime}:}{Object of class \code{"POSIXct"} the latest endTime of all spatio-temporal  neighbourhoods. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{ST}"}, directly.
}
\section{Methods}{
No methods defined with class "stNeighbourhood" in the signature.
}
\author{
Benedikt Graeler
}

\seealso{\code{\linkS4class{neighbourhood}}, \code{\link{getStNeighbours}}}

\examples{
showClass("stNeighbourhood")
}
\keyword{classes}
