\name{reduceNeighbours}
\alias{reduceNeighbours}

\title{
Selecting the strongest correlated neighbours
}
\description{
A function selecting the strongest correlated neighbours from a larger set of neighbours
}
\usage{
reduceNeighbours(stNeigh, stDepFun, n)
}

\arguments{
  \item{stNeigh}{
the proxy neighbourhood to be investigated
}
  \item{stDepFun}{
a spatio-temporal dependence function that return correlation estimates based on a spatial and temporal distance
}
  \item{n}{
the number of neighbours to be selected.
}
}
\value{A spatio-temporal neighbourhood \code{\linkS4class{stNeighbourhood}} with fewer neighbours.}

\author{
Benedikt Graeler
}


\seealso{
\code{\link{getStNeighbours}}
}

\examples{
library(sp)
library(spacetime)

sp <- SpatialPoints(matrix(c(181000,181100,333500,333600),2))
time <- Sys.time()+60*60*24*c(0,1,2,3,4)
data <- data.frame(var1=runif(10))

stData <- STFDF(sp, time, data)

stNeigh <- getStNeighbours(stData, spSize=2, t.lags=-(0:2))

reduceNeighbours(stNeigh, function(h,delta) return(1/h/delta), 2)
}