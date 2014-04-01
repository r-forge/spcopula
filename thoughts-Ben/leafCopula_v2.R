

ellipse <- function(x, a, b) {
  -a*(1-(x/b)^2)^(1/2)+a
}

ddxEllipse <- function(x, a, b) {
  a*x/b^(2)/sqrt(1-(x/b)^2)
}

## double checking
curve(ellipse,0,0.4,asp=1)
curve(ddxEllipse,add=T,col="red")
abline(h=0)
abline(ellipse(0.2)-ddxEllipse(0.2)*0.2,ddxEllipse(0.2), col="green")
abline(ellipse(0.3)-ddxEllipse(0.3)*0.3,ddxEllipse(0.3), col="green")
abline(v=c(0.2,0.3), col="grey")
## okay

# t to x scale:

# rotate t into nice ellipse coordinates, works only with k=2:
# i.e. 45° counter clockwise
rotate_t <- function(t, a=0.2, b=0.4) {
  c <- sqrt(2)
  (sqrt(2*c*a^3*b^2*t+a^2*b^4-2*a^2*b^2*t^2)-a*b^2+c*b^2*t)/(a^2+b^2)
}

ddtRotate_t <- function(t, a=0.2, b=0.4) {
  c <- sqrt(2)
  cRoot <- sqrt(2*t*(a*c-t)+b^2)
  (b*(a^2*c+b*c*cRoot-2*a*t))/((a^2+b^2)*cRoot)
}

## double check
curve(rotate_t,0,0.5)
curve(ddtRotate_t,0,0.5,add=T)
abline(rotate_t(0.2)-ddtRotate_t(0.2)*0.2, ddtRotate_t(0.2), col="green")
abline(rotate_t(0.3)-ddtRotate_t(0.3)*0.3, ddtRotate_t(0.3), col="green")
abline(v=c(0.2,0.3), col="grey")
## okay

# scale_t <- function(t, t1=0.5-sqrt(0.5)*b, t3=0.5+sqrt(0.5)*a) {
#   (t-t1)/(t3-t1)
# }

# rotation of ellipse from "nice" coordinates to original A(t) in t:
# i.e.: 45° clockwise
ellipseInT <- function(t, a=0.2, b=0.4) {
  x <- rotate_t(t, a, b)
  -sqrt(0.5)*(x-ellipse(x, a, b))
}

curve(ellipseInT,0,0.4,asp=1)
abline(0,-1,col="grey")

ddtEllipseInT <- function(t, a=0.2, b=0.4) {
  ddtx <- ddtRotate_t(t, a, b)
  
  -sqrt(0.5)*(ddtx-ddxEllipse(rotate_t(t, a, b), a, b)*ddtx)
}

## double check
curve(ellipseInT,0,0.5)
curve(ddtEllipseInT,0,0.5,add=T)
abline(ellipseInT(0.2)-ddtEllipseInT(0.2)*0.2, ddtEllipseInT(0.2), col="green")
abline(ellipseInT(0.3)-ddtEllipseInT(0.3)*0.3, ddtEllipseInT(0.3), col="green")
abline(v=c(0.2,0.3), col="grey")
## okay


##

# Aellipse <- function(t, aRate=0.7, bRate=0.3) {
Aellipse <- function(copula, w) {
  t <- w
  aRate=0.7
  bRate=0.3
  a <- aRate*sqrt(0.5)
  b <- bRate*sqrt(0.5)
  t1=0.5-sqrt(0.5)*b
  t3=0.5+sqrt(0.5)*a

  res <- t
  boolOut <- (t <= t1 | t >= t3)
  res[boolOut] <- pmax(1-t[boolOut],t[boolOut])
  res[!boolOut] <- ellipseInT(t[!boolOut]-t1, a, b)+1-t1
  res
}

curve(1-x,0,0.5,col="grey",xlim=c(0,1),ylim=c(0.5,1),asp=1)
curve((x),0.5,1,add=T,col="grey")
segments(0,1,1,1,col="grey")
curve(Aellipse,asp=1,add=T)

##

# ddtAellipse <- function(t, aRate=0.7, bRate=0.3) {

ddtAellipse <- function(copula, w) {
  t <- w
  aRate=0.7
  bRate=0.3
  a <- aRate*sqrt(0.5)
  b <- bRate*sqrt(0.5)
  t1=0.5-sqrt(0.5)*b
  t3=0.5+sqrt(0.5)*a
  
  res <- t
  res[t <= t1] <- -1
  res[t >= t3] <- 1
  
  boolInner <- (t > t1 & t < t3)
  res[boolInner] <- ddtEllipseInT(t[boolInner]-t1, a, b)
  res
}


## double check
curve(Aellipse,0,1,ylim=c(0.5,1))
curve(ddtAellipse,0,,add=T)
abline(Aellipse(0.2)-ddtAellipse(0.2)*0.2, ddtAellipse(0.2), col="green")
abline(Aellipse(0.5)-ddtAellipse(0.5)*0.5, ddtAellipse(0.5), col="green")
abline(Aellipse(0.7)-ddtAellipse(0.7)*0.7, ddtAellipse(0.7), col="green")
abline(v=c(0.2,0.5,0.7), col="grey")
## okay

plot(ecdf(tSmpl),xlim=c(0,1))

plickandsCDF <- function(t, t1=0.35, t3=0.75) {
  stopifnot(t1 < 0.5)
  stopifnot(t3 > 0.5)
  aRate <- (t3-0.5)/0.5
  bRate <- (0.5-t1)/0.5
  t+t*(1-t)*ddtAellipse(t,aRate, bRate)/Aellipse(t,aRate, bRate)
}

curve(plickandsCDF, add=T, col="green")



setMethod("dAdu",signature("tawn3pCopula"),ddtAellipse)
setMethod("A",signature("tawn3pCopula"),Aellipse)

plot(rCopula(500,tawn3pCopula()))















#####
cdf.leaf <- function(uvab) {
  wBor <- weakBorderPoly(uvab[1])
  if( uvab[2] >= wBor & uvab[2] <= strongBorderPoly(uvab[1]))
    return(wBor+ellipseOrig(uvab[2]-wBor,uvab[3],uvab[4]))
  return(min(uvab[1:2]))
}

ddv.cdf.leaf <- function(uvab) {
  wBor <- weakBorderPoly(uvab[1])
  sBor <- strongBorderPoly(uvab[1])
  if (uvab[2] <= wBor)
    return(1)
  if (uvab[2] >= sBor)
    return(0)
  return(ddvellipseOrig(uvab[2]-wBor, uvab[3], uvab[4]))
}

invddvLeafCop <- function(vy) {
  
  optFun <- function(u) {
    ret <- apply(cbind(u,rep(vy[1],length(u)),solveA(u),solveB(u)),1,ddv.cdf.leaf)
    sqrt(mean((ret-vy[2])^2))
  }
  invWBor <- invWeakBor(vy[1])
  invSBor <- invStrongBor(vy[1])
  if(invSBor == invWBor)
    return(invSBor)
  optimise(optFun,c(min(invSBor,invWBor),max(invSBor,invWBor)))$minimum
}

ddu.cdf.leaf <- function(uvab) {
  wBor <- weakBorderPoly(uvab[1])
  sBor <- strongBorderPoly(uvab[1])
  if (uvab[2] <= wBor)
    return(0)
  if (uvab[2] >= sBor)
    return(1)
  return((ddxweakBorderPoly(uvab[1])*(1-ddvellipseOrig(uvab[2]-wBor, uvab[3], uvab[4]))
          + dduellipseOrig(uvab[2]-wBor, uvab[3], uvab[4], ddusolveA(uvab[1]), ddusolveB(uvab[1]))))
}

pdf.leaf <- function(uvab) {
  wBor <- weakBorderPoly(uvab[1])
  sBor <- strongBorderPoly(uvab[1])
  if (uvab[2] <= wBor)
    return(0)
  if (uvab[2] >= sBor)
    return(0)
  return(ddvuEllipseOrig(uvab[2]-wBor, uvab[3], uvab[4], ddusolveA(uvab[1]), ddusolveB(uvab[1]))
         - ddxweakBorderPoly(uvab[1])*ddvvEllipseOrig(uvab[2]-wBor, uvab[3], uvab[4]))
}

## random number generator
r.leaf <- function(n) {
  v <- runif(n, min = 0, max = 1)
  y <- runif(n, min = 0, max = 1)
    
  res <- cbind(apply(cbind(v,y),1,invddvLeafCop),v)
  colnames(res) <- c("u","v")
    
  return(res)
}

pLeafCopula <- function(u) {
aVec <- solveA(u[,1])
bVec <- solveB(u[,1])

apply(cbind(u,aVec,bVec),1,cdf.leaf)
}

ddvLeafCopula <- function(u) {
  aVec <- solveA(u[,1])
  bVec <- solveB(u[,1])
  
  apply(cbind(u,aVec,bVec),1,ddv.cdf.leaf)
}

dduLeafCopula <- function(u) {
  aVec <- solveA(u[,1])
  bVec <- solveB(u[,1])
  
  apply(cbind(u,aVec,bVec),1,ddu.cdf.leaf)
}

dLeafCopula <- function(u) {
  aVec <- solveA(u[,1])
  bVec <- solveB(u[,1])
  
  apply(cbind(u,aVec,bVec),1,pdf.leaf)
}

####

pLeafCop <- data.frame(u=rep((1:100)/101,each=100), v=rep((1:100)/101,100),
                       p=pCopula(cbind(rep((1:100)/101,each=100),rep((1:100)/101,100)),
                                 spcopula:::leafCopula()))

u=rep((1:200)/201,200)
v=c(rep((1:200)/201,each=200))#,strongBorderPoly((1:100)/101))




dLeafCop <- data.frame(x=u, y=v, d=pmin(dCopula(cbind(u,v), spcopula:::leafCopula()),50))
dLeafCop <- 
  dLeafCop[dLeafCop$d<=0,] <- NA

library(rgl)
surface3d((1:100)/101, (1:100)/101, as.matrix(pLeafCop$p,nrow=100,byrow=T),col=terrain.colors(100))
surface3d((1:200)/201, (1:200)/201, as.matrix(dLeafCop$d/50,nrow=200),
          col=terrain.colors(51)[round(as.matrix(dLeafCop$d,nrow=200),0)+1])

summary(dLeafCop$d)

plot3d(dLeafCop,zlim=c(5,50))

rgl.surface(matrix(u,101),t(matrix(v,100)),
            y=matrix(pmin(dCopula(cbind(u,v), spcopula:::leafCopula()),50)/100,101))

persp(spcopula:::leafCopula(),dCopula,zlim=c(0,20))

z <- 2 * volcano        # Exaggerate the relief

x <- 10 * (1:nrow(z))   # 10 meter spacing (S to N)
y <- 10 * (1:ncol(z))   # 10 meter spacing (E to W)

zlim <- range(z)
zlen <- zlim[2] - zlim[1] + 1

colorlut <- terrain.colors(zlen) # height color lookup table

col <- colorlut[ z-zlim[1]+1 ] # assign colors to heights for each point

open3d()
surface3d(x, y, z, color=col, back="lines")


plot_rgl_model_a(dLeafCop)

foo <- dLeafCop[which(dLeafCop$d<0),][1,]

filled.contour(matrix(dLeafCop$d,100),asp=1,key=F)
points(rtTriples[,c(3,1)])

wireframe(d~x+y,dLeafCop,zlim=c(0,50), drape=T)

dLeafCopula(foo[1:2])

pdf.leaf(uvab cbind(foo[1:2],solveA(foo[1]),solveB(foo[1])))

dduellipseOrig(foo[2]-weakBorderPoly(foo[1]),solveA(foo[1]),solveB(foo[1]),ddusolveA(foo[1]),ddusolveB(foo[1]))

ddvuEllipseOrig(uvab[2]-weakBorderPoly(uvab[1]), uvab[3], uvab[4], 
                ddusolveA(uvab[1]), ddusolveB(uvab[1]))

aVec <- solveA(foo[,1])
bVec <- solveB(foo[,1])

apply(cbind(foo[1:2],aVec,bVec),1,pdf.leaf)

plot(ddvLeafCopula(cbind(rep(.9,100), (1:100)/101)),typ="l")
plot(ddvLeafCopula(cbind((1:100)/101,rep(.9,100))),typ="l")

plot(dduLeafCopula(cbind((1:100)/101,rep(.9,100))),typ="l")
plot(dduLeafCopula(cbind(rep(.2,100),(1:100)/101)),typ="l")

plot(dduLeafCopula(cbind(rep(.85,100),(1:100)/101)),typ="l")
plot(pLeafCopula(cbind(rep(.85,100),(1:100)/101)),typ="l")
plot(dduLeafCopula(cbind((1:100)/101,rep(93/101,100))),typ="l")

plot(pLeafCopula(cbind((750:1000)/1001,rep(96/101,251))),typ="l", ylim=c(0.5,1))
lines(pLeafCopula(cbind((750:1000)/1001,rep(95/101,251))),typ="l",col=2)
lines(pLeafCopula(cbind((750:1000)/1001,rep(94/101,251))),typ="l",col=3)
lines(pLeafCopula(cbind((750:1000)/1001,rep(93/101,251))),typ="l",col=4)

plot(dduLeafCopula(cbind((750:1000)/1001,rep(96/101,251))),typ="l")
lines(dduLeafCopula(cbind((750:1000)/1001,rep(95/101,251))),typ="l",col=2)
lines(dduLeafCopula(cbind((750:1000)/1001,rep(94/101,251))),typ="l",col=3)
lines(dduLeafCopula(cbind((750:1000)/1001,rep(93/101,251))),typ="l",col=4)
abline(v=750/1001*101)

## double check
plot((750:1000)/1001,pLeafCopula(cbind((750:1000)/1001,rep(0.95,251))),typ="l")
abline(pLeafCopula(matrix(c(0.9,0.95),ncol=2))-dduLeafCopula(cbind(0.9,0.95))*0.9,dduLeafCopula(cbind(0.9,0.95)))
abline(pLeafCopula(matrix(c(0.84,0.95),ncol=2))-dduLeafCopula(cbind(0.84,0.95))*0.84,dduLeafCopula(cbind(0.84,0.95)))
abline(v=c(0.84,0.9), col="grey")
## okay

plot(r.leaf(500),asp=1)
points(rtTriples[,c(3,1)],col="green",cex=1)

plot_rgl_model_a(

plot(rtTriples[,c(3,1)],asp=1)
poin