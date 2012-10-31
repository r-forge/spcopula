\name{calcBins}
\alias{calcBins}

\alias{calcBins-methods}
\alias{calcBins,Spatial-method}
\alias{calcBins,STFDF-method}

\title{
A function calculating the spatial/spatio-temporal bins
}
\description{
The (spatio-temporal) space is subdivided into pairs of observations that belkopng to certain spatial/spatio-temporal distanxce classes. For each distance class, the mean separating distance of all pairs involved is calcualted alongside a correaltion measure. The spatial/spatio-temporal correlogram is plotted by default.
}
\usage{
calcBins(data, var, nbins = 15, boundaries = NA, cutoff = NA, cor.method="kendall", plot=T, ...)
}

\arguments{
  \item{data}{
Either some spatial \code{data.frame} extending \code{\link[sp:Spatial-class]{Spatial}} or a spatio-temporal \code{data.frame} extending \code{\linkS4class{ST}}.
}
  \item{var}{
A vector containing the variable names of interset to select a subset of the variables present in the data slot.
}
  \item{nbins}{
the number of spatial bins.
}
  \item{boundaries}{
The left hand-side limits of the lag classes to be investigated, if set to \code{NA}, the \code{cutoff} will be split in equal parts according to \code{nbins}.
}
  \item{cutoff}{
The largest spatial distance to be investiagted, if set to \code{NA}, one third of the boundingbox's diagonal will be used.
}
  \item{cor.method}{
  defining the correlation method that should be used. Possible choices are \code{kendall} (default), \code{pearson}, \code{spearman} and \code{fasttau} that re-uses a very fast implementation of Kendall's tau from \code{\link{CDVine-package}}.
  }
  \item{plot}{Whether the correlogram should be plotted.}
  \item{\dots}{Additioanl arguments for the spatio-temporal case:
  \describe{
    \item{instances}{To reduce the data size or circumvent unwanted autocorrelation effects, one might provide a number of randomly selected time instances from the spatio-temporal \code{data.frame}. If this parameter is set to \code{NA}, th complete time series will be used, if different from a signle number, \code{instances} will be passed on as to index time.}
    \item{t.lags}{a vector indicating the time lags to be investigated}
}
}
}
\value{A list containg the following entries
\item{meanDists}{A vector holding the mean separating distance for each lag.}
\item{lagCor}{A vector holding the calculated correlation measure for each lag.}
\item{lagData}{the data pairs for each lag}
\item{lags}{In the spatial case, a vector holding the indicies of loactions used for each lag. In the spatio-temporal case, a list with a vector of spatail and temporal indicies used in the lag classes each.}

}
\author{
Benedikt Graeler
}
\seealso{
\code{\link{CDVine-package}}
}
\examples{
data(meuse)
coordinates(meuse) = ~x+y
dataSet <- meuse

## lag classes ##
bins <- calcBins(dataSet,var="zinc",nbins=10,cutoff=800)
}

\keyword{spatial preparation}
\keyword{spatio-temporal preparation}
