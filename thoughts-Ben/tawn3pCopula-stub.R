#################
## tawn copula ##
#################

library(copula)

Atawn3p <- function(t, param = c(0.9302082, 1, 8.355008)) {
  alpha <- param[1]
  beta  <- param[2]
  theta <- param[3]
  (1-beta)*(t) + (1-alpha)*(1-t) + ((alpha*(1-t))^theta+(beta*t)^theta)^(1/theta)
#   1-beta + (beta-alpha)*t + ((alpha*t)^theta + (beta*(1-t))^theta)^(1/theta)
}

curve(Atawn3p)

tawn3pCopula <- function (param = c(0.5, 0.5, 2)) {
  # A(t) = (1-beta)*(1-t) + (1-alpha)*t + ((alpha*(1-t))^theta+(beta*t)^theta)^(1/theta)
  # C(u1,u2) = exp(log(u1*u2) * A(log(u2)/log(u1*u2)))

  cdf <- expression(exp(log(u1*u2)*((1-beta)*(log(u2)/log(u1*u2)) +
                                    (1-alpha)*(1-log(u2)/log(u1*u2)) +
                                    ((alpha*(1-log(u2)/log(u1*u2)))^theta+(beta*log(u2)/log(u1*u2))^theta)^(1/theta))))
  dCdU1 <- D(cdf, "u1")
  pdf <- D(dCdU1, "u2")
  
  new("tawn3pCopula", dimension = 2L, exprdist = c(cdf = cdf, pdf = pdf),
      parameters = param, param.names = c("alpha", "beta", "theta"), 
      param.lowbnd = c(0,0,1), param.upbnd = c(1,1,Inf), 
      fullname = "Tawn copula family with three parameters; Extreme value copula")
}

setClass("tawn3pCopula", representation(exprdist = "expression"),
         contains = "evCopula")


dtawn3pCopula <- function(u, copula, log=FALSE, ...) {
  dim <- copula@dimension
  for (i in 1:dim) {
    assign(paste("u", i, sep=""), u[,i])
  }
  alpha <- copula@parameters[1]
  beta <-  copula@parameters[2]
  theta <-copula@parameters[3]
  
  val <- c(eval(copula@exprdist$pdf))
  ## FIXME: improve log-case
  if(log) log(val) else val
}

setMethod("dCopula", signature("matrix", "tawn3pCopula"), dtawn3pCopula)
setMethod("dCopula", signature("numeric", "tawn3pCopula"),dtawn3pCopula)

ptawn3pCopula <- function(u, copula, ...) {
  dim <- copula@dimension
  for (i in 1:dim) {
    assign(paste("u", i, sep=""), u[,i])
  }
  alpha <- copula@parameters[1]
  beta <-  copula@parameters[2]
  theta <-copula@parameters[3]
  
  val <- c(eval(copula@exprdist$cdf))
}

setMethod("pCopula", signature("matrix", "tawn3pCopula"),  ptawn3pCopula)
setMethod("pCopula", signature("numeric", "tawn3pCopula"), ptawn3pCopula)


persp(tawn3pCopula(c(0.25, 0.75, 2)), dCopula)
persp(tawn3pCopula(c(0.5, 1, 20)), pCopula)

plot(1-rtTriples[,c(1,3)])

dCopula(c(0.15,0.95), tawn3pCopula(c(0.5, 0.5, 2)))


tawnFit <- fitCopula(tawn3pCopula(), 1-as.matrix(rtTriples[,c(1,3)]), hideWarnings=F,estimate.variance=F,
                     start=c(0.9, 1, 8), method="mpl", lower=c(0,0,1), upper=c(1, 1, 10),
                     optim.method="L-BFGS-B",)
tawnFit@loglik # 742
tawnFit@copula

fitCopula(gumbelCopula(5),1-as.matrix(rtTriples[,c(1,3)]))@loglik # 723

dLeaf <- dCopula(as.matrix(rtTriples[,c(1,3)]), spcopula:::leafCopula()) 
sum(log(dLeaf[dLeaf>0]))

persp(tawnFit@copula, dCopula)
contour(tawnFit@copula, dCopula, levels=c(0,0.5,1,2,4,8,100), asp=1)

sum(dCopula(as.matrix(rtTriples[,c(1,3)]), cop13, log=T))
sum(dCopula(1-as.matrix(rtTriples[,c(1,3)]), tawnFit@copula, log=T))

plot(1-as.matrix(rtTriples[,c(1,3)]),asp=1,cex=0.5)
curve(x^(tawnFit@copula@parameters[1]),add=T, col="red")
abline(0,1,col="grey")

copula:::fitCopula.ml


tawn3pCopula()

persp(tawn3pCopula(),dCopula)


###
# h(t)

hist(log(1-rtTriples[,3])/log((1-rtTriples[,1])*(1-rtTriples[,3])),n=20,
     xlim=c(0,1), freq=F, add=T, col="blue")

tSmpl <- log(1-rtTriples[,3])/log((1-rtTriples[,1])*(1-rtTriples[,3]))
  ((1-rtTriples[,1])-(1-rtTriples[,3]))*(0.5)+0.5

dlogNorm <- function(x) dlnorm(x, mean(log(tSmpl)), sd(log(tSmpl)))
aGevPar <- fgev(tSmpl)$estimate
dGev <- function(x) dgev(x, 0.46938, 0.05057, 0.01720)
dGamma <- function(x) dgamma(x, 60.40556, 120.92807)
  
hist(tSmpl, freq=F, xlim=c(0,1), n=20, add=F)
curve(dlogNorm, add=T)
curve(dGev, add=T, col="red")
curve(dGamma, add=T, col="purple")

sum(log(dGev(tSmpl))) # 690
sum(log(dlogNorm(tSmpl))) # 666
sum(log(dGamma(tSmpl))) # 658


Afit <- function(t) {
  res <- t
  res[res == 0] <- 1
  intFun <- function(z) (pgev(z, 0.46938, 0.05057, 0.01720)-z)/(z-z^2)
  for(i in which(res != 1)) {
    res[i] <- exp(integrate(intFun,0,t[i])$value)
  }
  return(res)
}

curve(Afit,ylim=c(0.5,1))
abline(0, 1, col="grey")
abline(1, -1, col="grey")
curve(Atawn3p,col="red",add=T)


optFun <- function(param) {
  -sum(log(dgamma(tSmpl,param[1],param[2])))
}

optim(c(1,0.5),optFun)


library(evd)

?dnorm

##
AFun <- function(u, v, cdf=function(z) pgev(z, 0.46938, 0.05057, 0.01720)) {
  t <- log(v)/log(u*v)
  
  # rough limit case corrections
  t[v == 0 | u == 1 | v == 1 | u == 0] <- 1 
  
  res <- t
  intFun <- function(z) (cdf(z)-z)/(z-z^2)
  for(i in which(res < 1 & res > 0)) {
    res[i] <- exp(integrate(intFun,0,t[i],stop.on.error=F,subdivisions=1000L)$value)
  }
  return(res) #pmax(pmax(res,t),1-t))
}

plot(u,AFun(u),typ="l",ylim=c(0,1))
points(u,log(0.3)/log(0.3*u),typ="l",col="red") # slope less than this one
points(u,1-log(0.3)/log(0.3*u),typ="l",col="red") # slope larger than this one

AFunInT <- function(t, cdf=function(z) pgev(z, 0.5, 
                                            0.05,
                                            0.017)) {
  res <- t
  res[res == 0] <- 1
  intFun <- function(z) (cdf(z)-z)/(z-z^2)

  for(i in which(res < 1 & res > 0)) {
    res[i] <- exp(integrate(intFun,0,t[i])$value)
  }
  return(res)
}

curve(AFunInT)
abline(0,1)
abline(1,-1)

fooFun <- function(z) pgev(z, 0.5, 0.05, 0.017)
curve(fooFun,add=T)
abline(v=0.5)

aGevPar
newFit@estimate

curve(AFunInT)
abline(0,1)
abline(1,-1)


dduAFun <- function(u, v, cdf=function(z) pgev(z, 0.46938, 0.05057, 0.01720)) {
  logV  <- log(v)
  logUV <- log(u*v)
  
  t <- logV/logUV
  # limit cases: 
  # v == 0 -> t=1
  # v == 1 -> t=0
  # u == 1 -> t=1
  # u == 0 -> t=0
  t[v == 0] <- 1 
  t[v == 1] <- 0
  t[u == 1] <- 1 
  t[u == 0] <- 0 
  
  res <- t
  res[t == 0] <- -10
  res[t == 1] <- (-logV[t == 1])/u[t == 1]/logUV[t == 1]^2
  
  boolDo <- t > 0 & t < 1
  tDo <- t[boolDo]
  uDo <- u[boolDo]
  vDo <- v[boolDo]

  res[boolDo] <- AFun(uDo, vDo, cdf)*(cdf(tDo)-tDo)/(tDo-tDo^2)*(-logV[boolDo])/(u[boolDo]*logUV[boolDo]^2)
  
  return(res)
}

plot(u, dduAFun(u),typ="l")
points(u,-log(0.3)/log(0.3*u)^2/u, typ="l",col="red") # slope less than this one
points(u, log(0.3)/log(0.3*u)^2/u, typ="l",col="red") # slope larger than this one

curve(AFun,0,1,add=F,ylim=c(-2,1))
curve(dduAFun, add=T, col="green")
abline(h=0,col="grey")
abline(v=uniroot(dduAFun,c(0.2,0.4))$root, col="grey")

ddvAFun <- function(u, v=rep(0.3,length(u)), cdf=function(z) pgev(z, 0.46938, 0.05057, 0.01720)) {
  logV  <- log(v)
  logUV <- log(u*v)
  
  t <- logV/logUV
  # limit cases: 
  # v == 0 -> t=1
  # v == 1 -> t=0
  # u == 1 -> t=1
  # u == 0 -> t=0
  t[v == 0] <- 1 
  t[v == 1] <- 0
  t[u == 1] <- 1 
  t[u == 0] <- 0 
  
  res <- t
  res[t == 0] <- -(logUV[t == 0]-logV[t == 0])/(logUV[t == 0]^2*v[t == 0])
  res[t == 1] <- (logUV[t == 1]-logV[t == 1])/(logUV[t == 1]^2*v[t == 1])
  
  boolDo <- t > 0 & t < 1
  tDo <- t[boolDo]
  uDo <- u[boolDo]
  vDo <- v[boolDo]
  
  res[boolDo] <- AFun(uDo, vDo, cdf)*(cdf(tDo)-tDo)/(tDo-tDo^2)*(logUV[boolDo]-logV[boolDo])/(vDo*logUV[boolDo]^2)
  
  return(res)
}


plot(v, ddvAFun(v),typ="l")
points(v,-(log(0.3*v)-log(v))/(log(0.3*v)^2*v), typ="l",col="red") # slope less than this one
points(v, (log(0.3*v)-log(v))/(log(0.3*v)^2*v), typ="l",col="red") # slope larger than this one

curve(AFun,0,1,add=F,ylim=c(-2,1))
curve(ddvAFun, add=T, col="green")
abline(h=0,col="grey")
abline(v=uniroot(dduAFun,c(0.2,0.4))$root, col="grey")

dduvAFun <- function(u, v=rep(0.3,length(u)), 
                     cdf=function(z) pgev(z, 0.46938, 0.05057, 0.01720),
                     pdf=function(z) dgev(z, 0.46938, 0.05057, 0.01720)) {
  logV  <- log(v)
  logUV <- log(u*v)
  
  .g <- function(z) {
    (cdf(z)-z)/(z-z^2)
  }
  
  .ddzg <- function(z) {
    ((2*z-1)*cdf(z)-z*((z-1)*pdf(z)+z))/((z-1)^2*z^2)
  }
  
  t <- logV/logUV
  # limit cases
  t[v == 0] <- 1 
  t[v == 1] <- 0
  t[u == 1] <- 1 
  t[u == 0] <- 0 
  
  p1 <- dduAFun(u, v, cdf)*.g(t)*u*(logUV-logV)*logUV^2
  p2 <- AFun(u, v, cdf)*(logV*(logV-logUV) * .ddzg(t) + (2*logV-logUV)*logUV * .g(t) )
  
  (p1+p2)/(u*v*logUV^4)
}

curve(ddvAFun,0,1,ylim=c(-4,2),asp=1, add=T)
curve(dduvAFun,0.01,0.99,add=F)
abline(h=0, col="grey")

curve(AFun, add=T, col="red")
curve(dduAFun)
# curve(ddvAFun,add=T, col="red")
# curve(dduvAFun, add=F)

CDFevCopula <- function(u, v=rep(0.3,length(u)), 
                        cdf=function(z) pgev(z, 0.46938, 0.05057, 0.01720)) {
  res <- u*v
  boolDo <- res > 0 & res < 1 

  res[boolDo] <- exp(log((u*v)[boolDo])*AFun(u[boolDo], v[boolDo], cdf))
  return(res)
}

dduCDFevCop <- function(u, v=rep(0.3,length(u)),
                        cdf=function(z) pgev(z, 0.46938, 0.05057, 0.01720)) {
  res <- v
  
  # limit cases
  res[u == 0 & (res > 0 | res < 1)] <- 1
  res[u == 1 & (res > 0 | res < 1)] <- 0
  boolDo <- res > 0 & res < 1
  uDo <- u[boolDo]
  vDo <- v[boolDo]
  res[boolDo] <- CDFevCopula(uDo, vDo, cdf)/uDo
 
  res[boolDo] <- res[boolDo] * (uDo * log(uDo*vDo) * dduAFun(uDo, vDo, cdf)
                                + AFun(uDo, vDo, cdf))
  return(res)    
}

u <- 1:99/100

CDFevCopula(u)
dduCDFevCop(u)

curve(CDFevCopula,add=T)
curve(dduCDFevCop,add=F, col="red")
abline(0,1,col="grey")
abline(h=0.3, col="grey")

ddvCDFevCop <- function(u, v=rep(0.8,length(u)),
                        cdf=function(z) pgev(z, 0.46938, 0.05057, 0.01720)) {
  res <- v
  
  # limit cases
  res[v == 1 & (res > 0 | res < 1)] <- 0
  res[v == 0 & (res > 0 | res < 1)] <- 1
  boolDo <- res > 0 & res < 1
  uDo <- u[boolDo]
  vDo <- v[boolDo]
  res[boolDo] <- CDFevCopula(uDo, vDo, cdf)/vDo
  
  res[boolDo] <- res[boolDo] * (vDo * log(uDo*vDo) * ddvAFun(uDo, vDo, cdf)
                                + AFun(uDo, vDo, cdf))
  return(res)    
}

CDFevCopula(u)
ddvCDFevCop(u)
curve(ddvCDFevCop,add=F, col="red")
curve(CDFevCopula,add=T)
abline(0,1,col="grey")
abline(h=0.3, col="grey")

PDFevCopula <- function(u, v=rep(0.7,length(u)),
                        cdf=function(z) pgev(z, 0.46938, 0.05057, 0.01720),
                        pdf=function(z) dgev(z, 0.46938, 0.05057, 0.01720)) {
  res <- u
  
  # limit cases
  res[v == 1 | v == 0 | u == 1 | u == 0] <- 0
  boolDo <- res > 0 & res < 1
  uDo <- u[boolDo]
  vDo <- v[boolDo]
  logUVDo <- log(uDo*vDo)
  AFunDo <- AFun(uDo, vDo, cdf)
  dduAFunDo <- dduAFun(uDo, vDo, cdf)
  ddvAFunDo <- ddvAFun(uDo, vDo, cdf)
  
  res[boolDo] <- CDFevCopula(uDo, vDo, cdf)/(vDo*uDo)
  
  p1 <- vDo*ddvAFunDo + uDo*dduAFunDo
  p2 <- (vDo*logUVDo * ddvAFunDo + AFunDo)*(uDo*logUVDo*dduAFunDo+AFunDo)
  p3 <- uDo*vDo*logUVDo*dduvAFun(uDo, vDo, cdf)
  
  res[boolDo] <- res[boolDo] * (p1 + p2 + p3)
  return(res)    
}

PDFevCopula(u)

curve(PDFevCopula,0.01,0.99, add=F)
curve(ddvCDFevCop, add=T)

dgev()

cdfAFunCopula <- function (param = c(0.46938, 0.05057, 0.01720)) {
  new("cdfAFunCopula", dimension = 2L,  
      cdf=function(z) pgev(z, param[1], param[2], param[3]),
      pdf=function(z) pgev(z, param[1], param[2], param[3]),
      parameters = param, param.names = c("location", "scale", "shape"), 
      param.lowbnd = c(-Inf,0,-Inf), param.upbnd = c(Inf,Inf,Inf), 
      fullname = "Extreme value copula family defining A(t) through an univariate evd")
}

setClass("cdfAFunCopula", representation(cdf = "function", pdf = "function"),
         contains = "copula")

setMethod("dCopula", signature("matrix",  "cdfAFunCopula"),  
          function(u, copula, log=FALSE, ...) {
            copula <- cdfAFunCopula(copula@parameters)
            res <- PDFevCopula(u[,1], u[,2], copula@cdf, copula@pdf)
            res <- pmax(res,1e-27)
            if(log)
              return(log(res))
            return(res)
          })

setMethod("dCopula", signature("numeric", "cdfAFunCopula"),
          function(u, copula) PDFevCopula(u[1], u[2], copula@cdf, copula@pdf))

persp(cdfAFunCopula(), dCopula, zlim=c(0,10))

setMethod("pCopula", signature("matrix",  "cdfAFunCopula"),
          function(u, copula) CDFevCopula(u[,1], u[,2], copula@cdf))
setMethod("pCopula", signature("numeric", "cdfAFunCopula"),
          function(u, copula) CDFevCopula(u[1], u[2], copula@cdf))

persp(cdfAFunCopula(), pCopula)

###

library(rgl)

u=rep((1:200)/201,200)
v=rep((1:200)/201,each=200)

pGevCop <- data.frame(u, v, p=pCopula(cbind(u,v), cdfAFunCopula()))
dGevCop <- data.frame(u, v, d=pmin(dCopula(cbind(u,v), cdfAFunCopula()),10))

surface3d((1:200)/201, (1:200)/201, as.matrix(dGevCop$d/10,nrow=200),
          col=terrain.colors(11)[round(as.matrix(dGevCop$d,nrow=200),0)+1])

surface3d((1:200)/201, (1:200)/201, as.matrix(pLeafCop$p,nrow=200),
          col=terrain.colors(51)[round(as.matrix(pLeafCop$p,nrow=200)*50,0)+1])

sum(log(dCopula(as.matrix(1-rtTriples[,c(1,3)]), tawnFit@copula)))  # Tawn 3p:  742
sum(log(dCopula(as.matrix(1-rtTriples[,c(1,3)]), cdfAFunCopula(aGevPar)))) # cdf AFun: 766
sum(log(dCopula(as.matrix(rtTriples[,c(1,3)]), cop13)))             # BB7:      766
fitCopula(gumbelCopula(5),1-as.matrix(rtTriples[,c(1,3)]))@loglik   # Gumbel:   723


newFit <- fitCopula(cdfAFunCopula(aGevPar), 1-as.matrix(rtTriples[,c(1,3)]),method="mpl",
                    lower=aGevPar*0.8, upper=1.2*aGevPar, optim.method="L-BFGS-B",
                    start=aGevPar, estimate.variance=F, hideWarnings=F)

newFit@loglik
newFit@estimate-aGevPar*0.9
newFit@estimate-aGevPar*1.1
newFit@estimate

log(1e-10)
