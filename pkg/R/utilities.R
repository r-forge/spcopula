## utilities

# calib Kendall's tau and Spearman's Rho for the indepCopula
setMethod("calibKendallsTau", signature("indepCopula"), function(copula, tau) return(numeric(0)) )
setMethod("calibSpearmansRho", signature("indepCopula"), function(copula, rho) return(numeric(0)) )

# ranks are automatically removed and NAs are by default randomly distributed
rankTransform <- function(u,v=NULL, ties.method="average") {
  if(is.matrix(u) && ncol(u)==2) {
    v <- u[,2]
    u <- u[,1]
  } else {
    if(is.null(v)) stop("u must either be a matrix with 2 columns or u and v must be given.")
  }
  bool <- !(is.na(u)*is.na(v))
  cbind(rank(u[bool], na.last=NA, ties.method=ties.method),rank(v[bool], na.last=NA, ties.method=ties.method))/(sum(bool)+1)
}

# strength of dependence scatterplot
myPanel.smoothScatter <- function(x, y = NULL, nbin = 64, cuts = 255, bandwidth, colramp, 
  transformation = function(x) x, pch = ".", cex = 1, col = "black", range.x, ..., raster = FALSE, subscripts) {
  if (missing(colramp)) 
    colramp <- function(x) rev(heat.colors(x))
  x <- as.numeric(x)
  y <- as.numeric(y)
  if(min(x)<0 | max(x)>1 | min(y) < 0 | max(y) > 1) {
    x <- rankTransform(x)
    y <- rankTransform(y)
    warning("At least one of the margins seems to exceed [0,1] and has been transformed using rankTransform.")
  }
  xy <- xy.coords(x, y)
  x <- cbind(xy$x, xy$y)[!(is.na(xy$x) | is.na(xy$y)), , drop = FALSE]
  if (nrow(x) < 1) return()
  map <- lattice:::.smoothScatterCalcDensity(x, nbin, bandwidth, range.x)
  xm <- map$x1
  ym <- map$x2
  dens <- map$fhat
  dens <- array(transformation(dens), dim = dim(dens))
  PFUN <- if (raster) 
    panel.levelplot.raster
  else panel.levelplot
  PFUN(x = rep(xm, length(ym)), y = rep(ym, each = length(xm)), 
       z = as.numeric(dens), subscripts = TRUE, at = seq(from = 0, 
       to = 1.01 * max(dens), length = cuts + 2), col.regions = colramp(cuts + 1), ...)
}

dependencePlot <- function(formula=NULL, smpl, cuts=15, bandwidth=.075, transformation=function (x) x, ..., range.x=list(c(0,1),c(0,1))) {
  smpl <- as.data.frame(smpl)
  if(is.null(formula)) {
    if (ncol(smpl)>2) {
      warning("smpl contains more than 2 columns and no formula is given. The first two columns are plotted.")
      smpl <- as.data.frame(smpl[,1:2])
    }
    colnames(smpl) <- c("u","v")
    formula <- u~v
  }

  xyplot(formula, smpl, panel=myPanel.smoothScatter,
  aspect="iso", xlim=c(0,1), ylim=c(0,1), cuts=cuts, bandwidth=bandwidth, transformation=transformation, range.x=range.x, ...)
}

##
unitScatter <- function(formula=NULL, smpl, cuts=15, bandwidth=.075, transformation=function (x) x, ...) {
  smpl <- as.data.frame(smpl)
  if(is.null(formula)) {
    if (ncol(smpl)>2) {
      warning("smpl contains more than 2 columns and no formula is given. The first two columns are plotted.")
      smpl <- as.data.frame(smpl[,1:2])
    }
    colnames(smpl) <- c("u","v")
    formula <- v~u
  }

  for(variable in all.vars(formula)){
    if( min(smpl[,variable])<0 | max(smpl[,variable])>1) {
      smpl[,variable] <- rankTransform(smpl[,variable])
      warning("The variable ",variable," seems to exceed [0,1] and has been transformed using rankTransform.")
    }
  }

  xyplot(formula, smpl, aspect="iso", xlim=c(0,1), ylim=c(0,1), ...)
}

univScatter <- function(formula=NULL, smpl, cuts=15, bandwidth=.075, transformation=function (x) x, ...) {
  warning("Use unitScatter instead!")
  unitScatter(formula, smpl, cuts, bandwidth, transformation, ...)
}