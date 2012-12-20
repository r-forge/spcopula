\name{neighbourhood-class}
\Rdversion{1.1}
\docType{class}
\alias{neighbourhood-class}
\alias{names,neighbourhood-method}

\title{Class \code{"neighbourhood"}}
\description{
A class representing a local spatial neighbourhood.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("neighbourhood", ...)} and from the constructor \code{\link{neighbourhood}} or through the easier function \code{\link{getNeighbours}}.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"data.frame"} a data.frame holding the data for each local neighbourhood. }
    \item{\code{distances}:}{Object of class \code{"matrix"} a matrix providing the separating distances per neighbourhood. }
    \item{\code{index}:}{Object of class \code{"matrix"} a matrix pointing to the IDs of the underlying spatial data.frame }
    \item{\code{varNames}:}{Object of class \code{"character"} the variable names of interest. }
    \item{\code{coords}:}{Object of class \code{"matrix"} the coordinates of the underlying spatial data frame. }
    \item{\code{bbox}:}{Object of class \code{"matrix"} the bounding box of the underlying spatial data frame. }
    \item{\code{proj4string}:}{Object of class \code{"CRS"} the coordinate reference system of the underlying spatial data frame. }
  }
}
\section{Extends}{
Class \code{"\link[sp:SpatialPoints-class]{SpatialPoints}"}, directly.
Class \code{"\link[sp:Spatial-class]{Spatial}"}, by class "SpatialPoints", distance 2.
}
\section{Methods}{
  \describe{
    \item{names}{\code{signature(x = "neighbourhood")}: ... }
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

colnames(neigh@data) <- rep("measure",3)

# rebuilding neigh
neighbourhood(data=neigh@data, neigh@distances, SpatialPoints(neigh@coords), neigh@index)
neigh
}
\keyword{classes}