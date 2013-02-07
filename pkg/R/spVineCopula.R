#########################################
## methods for the spatial vine copula ##
#########################################

# constructor
spVineCopula <- function(spCop, vineCop) {
  new("spVineCopula", dimension = as.integer(vineCop@dimension+1), parameters=numeric(),
      param.names = character(), param.lowbnd = numeric(), 
      param.upbnd = numeric(), fullname = "Spatial vine copula family.",
      spCop=spCop, vineCop=vineCop)
}

# show
showSpVineCopula <- function(object) {
  dim <- object@dimension
  cat(object@fullname, "\n")
  cat("Dimension: ", dim, "\n")
}

setMethod("show", signature("spVineCopula"), showSpVineCopula)

# density
dspVine <- function(u, spCop, vine, log, h) {
  l0 <- rep(0,nrow(u)) # level 0 (spatial) density
  u0 <- NULL # level 0 conditional data
  
  if(!is.matrix(h)) h <- matrix(h, ncol=length(h))
  
  for(i in 1:(ncol(u)-1)) { # i <- 1
    l0 <- l0+dCopula(as.matrix(u[,c(1,i+1)]), spCop, h=h[,i], log=T)
    u0 <- cbind(u0, dduCopula(as.matrix(u[,c(1,i+1)]), spCop, h=h[,i]))
  }
  
  l1 <- dCopula(u0, vine, log=T)
  if(log)
    return(l0+l1)
  else(exp(l0+l1))
}

setMethod("dCopula",signature=signature("matrix","spVineCopula"),
          function(u, copula, ...) {
            dspVine(u, copula@spCop, copula@vineCop, ...)
          })

setMethod("dCopula",signature=signature("numeric","spVineCopula"),
          function(u, copula, ...) {
            dspVine(matrix(u,ncol=copula@dimension), copula@spCop, copula@vineCop, ...)
          })