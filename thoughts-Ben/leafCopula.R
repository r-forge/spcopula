# assume to be present:
weakBorderPoly <- myPol # par[1]*x^3+par[2]*x^2+x
ddxweakBorderPoly <- function(x) {
  par <- optRet$par
  3*par[1]*x^2+2*par[2]*x+1
}

invWeakBor <- function(v) {
  optFun <- function(u) {
    sqrt(mean((weakBorderPoly(u)-v)^2))
  }
  
  optimise(optFun,c(0,1))$minimum
}

strongBorderPoly <- myBorderPoly
# function(x) {
#   a <- -1
#   b <- -1/2*(1+3*a)
#   a*x^3+b*x^2+x
# }
ddxstrongBorderPoly <- function(x) {
  a <- -1
  b <- -1/2*(1+3*a)
  3*a*x^2+2*b*x+1
}

invStrongBor <- function(v) {
  optFun <- function(u) {
    sqrt(mean((strongBorderPoly(u)-v)^2))
  }
  
  optimise(optFun,c(0,1))$minimum
}


# non visible functions
# precalculate ellipse parameters
solveQ <- function(u) {
  sqrt(0.5*(strongBorderPoly(u)-u)^2)
}

ddxsolveQ <- function(u) {
  sBor <- strongBorderPoly(u)
  1/(2*sqrt(0.5*(sBor-u)^2))*(sBor-u)*(ddxstrongBorderPoly(u)-1)
}

## double check
curve(solveQ, ylim=c(-0.6,0.3)) 
curve(ddxsolveQ, add=T) 
abline(h=0)
abline(solveQ(0.5)-ddxsolveQ(0.5)*0.5,ddxsolveQ(0.5))
abline(solveQ(0.9)-ddxsolveQ(0.9)*0.9,ddxsolveQ(0.9))
abline(v=c(0.5,0.9),col="grey")
##

solveXb <- function(u) {
  sqrt(2*(u-weakBorderPoly(u))^2)+solveQ(u)
}

ddxsolveXb <- function(u) {
  wBor <- weakBorderPoly(u)
  -(sqrt(2)*(u-wBor)*(ddxweakBorderPoly(u)-1))/sqrt((u-wBor)^2)+ddxsolveQ(u)
}

## double check
curve(solveXb,ylim=c(-1,1.2))
curve(ddxsolveXb, add=T)
abline(h=0)
abline(solveXb(0.5)-ddxsolveXb(0.5)*0.5,ddxsolveXb(0.5))
abline(solveXb(0.9)-ddxsolveXb(0.9)*0.9,ddxsolveXb(0.9))
abline(v=c(0.5,0.9), col="grey")
## okay

###
# ellipse(xb)=b-q:
# b = -q/(sqrt((a^2-x^2)/a^2)-1)
# 
# ellipse'(xb)=1:
# b = (a^2 sqrt((a^2-x^2)/a^2))/x
#
# solve for a:
# -q/(sqrt((a^2-x^2)/a^2)-1) = (a^2 sqrt((a^2-x^2)/a^2))/x
###

solveA <- function(u) { #(xb,q) {
  xb <- solveXb(u)
  q <- solveQ(u)
  sqrt((-xb^3+2*q*xb^2-q^2*xb)/(-2*xb+2*q+xb))
}

ddusolveA <- function(u) { # (xb,q,dduXb,dduQ) {
  xb <- solveXb(u)
  q <- solveQ(u)
  dduXb <- ddxsolveXb(u)
  dduQ <- ddxsolveQ(u)
  1/(2*solveA(u)) * (-2*(q-xb)*(q*xb*(dduQ-3*dduXb)+q^2*dduXb+xb^2*dduXb))/(-2*q + xb)^2
}

## double check
curve(solveA,ylim=c(-2,2))
curve(ddusolveA,add=T)
abline(h=0)
abline(solveA(0.5)-ddusolveA(0.5)*0.5,ddusolveA(0.5))
abline(solveA(0.9)-ddusolveA(0.9)*0.9,ddusolveA(0.9))
abline(v=c(0.5,0.9), col="grey")
## okay

solveB <- function(u) {
  a <- solveA(u)
  xb <- solveXb(u)
  a^2*sqrt(1-(xb/a)^2)/xb
} 

ddusolveB <- function(u) {
  a <- solveA(u)
  xb <- solveXb(u)
  dduA <- ddusolveA(u)
  dduXb <- ddxsolveXb(u)
  (2*a^2*xb*dduA - xb^3*dduA - a^3*dduXb)/(a*xb^2*sqrt(1 - xb^2/a^2)) 
}

## double check
curve(solveB,ylim=c(-0.6,0.3))
curve(ddusolveB,add=T)
abline(h=0)
abline(solveB(0.5)-ddusolveB(0.5)*0.5,ddusolveB(0.5))
abline(solveB(0.9)-ddusolveB(0.9)*0.9,ddusolveB(0.9))
abline(v=c(0.5,0.9), col="grey")
## okay

ellipse <- function(v,a,b) {
  b*sqrt(1-(v/a)^2)-b
}

## double check
ellipseInV <- function(v) ellipse(v,solveA(0.9),solveB(0.9))

curve(ellipseInV,asp=2)
abline(v=solveXb(0.9)+c(0,-solveQ(0.9)))
abline(solveXb(0.9)-solveQ(0.9),-1)
abline(h=0)
## okay

dduellipse <- function(v, a=solveA(0.9), b=solveB(0.9),
                       dduA=ddusolveA(0.9), dduB=ddusolveB(0.9)) {
  dduB*sqrt(1-(v/a)^2)+b/(2*sqrt(1-(v/a)^2))*2*v^2/a^3*dduA-dduB
}

## double checking
dduellipseInU <- function(u) dduellipse(0.5,solveA(u),solveB(u),ddusolveA(u),ddusolveB(u))
ellipseInU <- function(u) ellipse(0.5,solveA(u),solveB(u))

curve(ellipseInU,0.5,1,ylim=c(-0.2,.5))
curve(dduellipseInU,0.5,1,add=T)
abline(h=0)
abline(ellipseInU(0.6)-dduellipseInU(0.6)*0.6,dduellipseInU(0.6))
abline(ellipseInU(0.9)-dduellipseInU(0.9)*0.9,dduellipseInU(0.9))
abline(v=c(0.6,0.9), col="grey")
## okay

ddvellipse <- function(v, a=solveA(0.9),b=solveB(0.9)) {
  -b*v/sqrt(1-(v/a)^2)/a^2
}

## double checking
curve(ellipseInV,0,0.6,ylim=c(-.2,0))
curve(ddvellipse,0,0.6,asp=2,add=T)
abline(h=0)
abline(ellipseInV(0.2)-ddvellipseInV(0.2)*0.2,ddvellipseInV(0.2))
abline(ellipseInV(0.4)-ddvellipseInV(0.4)*0.4,ddvellipseInV(0.4))
abline(v=c(0.2,0.4), col="grey")
## okay

ddvuEllipse <- function(v, a=solveA(0.9),b=solveB(0.9),
                        dduA=ddusolveA(0.9), dduB=ddusolveB(0.9)) {
  (2*a^2*dduA*b*v - a^3*dduB*v - dduA*b*v^3 + a*dduB*v^3)/(a^3*(a^2 - v^2)*sqrt(1 - v^2/a^2))
}

## double checking
curve(dduellipse,ylim=c(-1,.5))
curve(ddvuEllipse,add=T) #,0,0.6,ylim=c(-.2,0))
abline(h=0)
abline(dduellipse(0.2)-ddvuEllipse(0.2)*0.2,ddvuEllipse(0.2))
abline(dduellipse(0.4)-ddvuEllipse(0.4)*0.4,ddvuEllipse(0.4))
abline(v=c(0.2,0.4), col="grey")
## okay


ddvvEllipse <- function(v, a=solveA(0.9),b=solveB(0.9)) {
  -b/((a^2 - v^2)*sqrt(1 - v^2/a^2))
}

## double checking
curve(ddvellipse,ylim=c(-1,.5))
curve(ddvvEllipse,add=T) #,0,0.6,ylim=c(-.2,0))
abline(h=0)
abline(ddvellipse(0.2)-ddvvEllipse(0.2)*0.2,ddvvEllipse(0.2))
abline(ddvellipse(0.4)-ddvvEllipse(0.4)*0.4,ddvvEllipse(0.4))
abline(v=c(0.2,0.4), col="grey")
## okay

vNorm <- function(vo, a=0.5428281, b=0.06500491) {
  c <- sqrt(0.5)
  (-a^2*b + 2*a^2*c*vo + 2*a*b*sqrt(0.25*a^2 + b*c*vo - 0.5*vo^2))/(a^2 + b^2)
}

dduvNorm <- function(vo, a=solveA(0.9), b=solveB(0.9), 
                     dduA=ddusolveA(0.9), dduB=ddusolveB(0.9)) {
  cRoot <- sqrt(-0.5*vo^2 + 0.25*a^2 + sqrt(0.5)*vo*b) 
  
  low <- a^2+b^2
  high <- -a^2*b+sqrt(2)*a^2*vo+2*a*b*cRoot
  
  dduLow <- 2*a*dduA + 2*b*dduB
  dduHigh <- (-2*a*dduA*b - a^2*dduB + 2*sqrt(2)*vo*a*dduA
              + 2*(dduA*b + a*dduB)*cRoot
              + a*b/cRoot*(0.5*a*dduA+sqrt(0.5)*vo*dduB))

  (low*dduHigh - high*dduLow)/low^2
}

## double check
vNormInU <- function(u) vNorm(.5,solveA(u),solveB(u))
dduvNormInU <- function(u) dduvNorm(.5,solveA(u),solveB(u),ddusolveA(u),ddusolveB(u))

curve(vNormInU,0.6,0.85,ylim=c(0.5,.7))
curve(dduvNormInU,add=T)
abline(vNormInU(0.65)-dduvNormInU(0.65)*0.65,dduvNormInU(0.65))
abline(vNormInU(0.8)-dduvNormInU(0.8)*0.8,dduvNormInU(0.8))
abline(v=c(0.65,0.8), col="grey")
## okay

ddvvNorm <- function(vo, a=solveA(0.9),b=solveB(0.9)) {
  c <- sqrt(0.5)
  (2*a^2*(c + (b^2*(0.5*b*c - 0.5*vo))/sqrt(a^2*b^2*(0.25*a^2 + (b*c - 0.5*vo)*vo))))/(a^2 + b^2)
}

## double check
curve(vNorm,0,0.5,ylim=c(0,1))
curve(ddvvNorm,add=T)
abline(vNorm(0.4)-ddvvNorm(0.4)*0.4,ddvvNorm(0.4))
abline(vNorm(0.1)-ddvvNorm(0.1)*0.1,ddvvNorm(0.1))
abline(v=c(0.1,0.4), col="grey")
## okay

ddvuvNorm <-  function(vo, a=solveA(0.9),b=solveB(0.9),
                       dduA=ddusolveA(0.9), dduB=ddusolveB(0.9)) {

  cRoot <- sqrt(-0.5*vo^2 + 0.25*a^2 + sqrt(0.5)*vo*b) 
  
  low <- a^2+b^2
  
  dduLow <- 2*a*dduA + 2*b*dduB
  ddvHigh <- sqrt(2)*a^2+(a*b*(-vo + sqrt(0.5)*b))/cRoot
  
  ddvuHigh <- (2*sqrt(2)*a*dduA 
               + ((dduA*b + a*dduB)*(-vo+sqrt(0.5)*b)+sqrt(1/2)*a*b*dduB)/cRoot
               - (a*b*(0.5*a*dduA+sqrt(0.5)*dduB*vo)*(-vo+sqrt(0.5)*b))/(2*cRoot^3))
    
  # ddv (low*dduHigh - high*dduLow):
  (low * ddvuHigh - ddvHigh*dduLow)/low^2
}

## double check
curve(dduvNorm,0,0.5,ylim=c(-2,1))
curve(ddvuvNorm,add=T)
abline(h=0)
abline(dduvNorm(0.2)-ddvuvNorm(0.2)*0.2,ddvuvNorm(0.2))
abline(dduvNorm(0.4)-ddvuvNorm(0.4)*0.4,ddvuvNorm(0.4))
abline(v=c(0.2,0.4), col="grey")
## okay

ddvvvNorm <- function(vo, a=solveA(0.9),b=solveB(0.9)) {
  c <- sqrt(0.5)
  (-2*a^4*b^4)/((a^2*b^2*(a^2 + 2*(2*b*c - vo)*vo))^(3/2))
}

## double check
curve(ddvvNorm,0,0.5)#,ylim=c(0,1))
curve(ddvvvNorm,add=T)
abline(ddvvNorm(0.4)-ddvvvNorm(0.4)*0.4,ddvvvNorm(0.4))
abline(ddvvNorm(0.1)-ddvvvNorm(0.1)*0.1,ddvvvNorm(0.1))
abline(v=c(0.1,0.4), col="grey")
## okay

yOrig <- function(vn, a=solveA(0.9), b=solveB(0.9)) {
  sqrt(0.5)*(vn+ellipse(vn, a, b))
}

dduyOrig <- function(vn, a=solveA(0.9), b=solveB(0.9),
                     dduA=ddusolveA(0.9), dduB=ddusolveB(0.9)) {
  sqrt(0.5)*dduellipse(vn, a, b, dduA, dduB)
}

## double check
yOrigInU <- function(u) yOrig(.5,solveA(u),solveB(u))
dduyOrigInU <- function(u) dduyOrig(.5,solveA(u),solveB(u),ddusolveA(u),ddusolveB(u))

curve(yOrigInU,0.5,1,ylim=c(0.2,.4))
curve(dduyOrigInU,add=T)
abline(yOrigInU(0.6)-dduyOrigInU(0.6)*0.6,dduyOrigInU(0.6))
abline(yOrigInU(0.9)-dduyOrigInU(0.9)*0.9,dduyOrigInU(0.9))
abline(v=c(0.6,0.9), col="grey")
## okay

ddvyOrig <- function(vn, a=solveA(0.9), b=solveB(0.9)) {
  sqrt(0.5)*(1+ddvellipse(vn, a, b))
}

## double check
curve(yOrig,0,0.5,ylim=c(0,1))
curve(ddvyOrig,add=T)
abline(yOrig(0.475)-ddvyOrig(0.475)*0.475,ddvyOrig(0.475))
abline(yOrig(0.1)-ddvyOrig(0.1)*0.1,ddvyOrig(0.1))
abline(v=c(0.1,0.475), col="grey")
## okay

ddvuyOrig <- function(vn, a=solveA(0.9), b=solveB(0.9),
                      dduA=ddusolveA(0.9), dduB=ddusolveB(0.9)) {
  sqrt(0.5)*ddvuEllipse(vn, a, b, dduA, dduB)
}

## double check
curve(dduyOrig,0,0.5,ylim=c(-0.2,0.4))
curve(ddvuyOrig,add=T)
abline(dduyOrig(0.475)-ddvuyOrig(0.475)*0.475,ddvuyOrig(0.475))
abline(dduyOrig(0.1)-ddvuyOrig(0.1)*0.1,ddvuyOrig(0.1))
abline(v=c(0.1,0.475), col="grey")
## okay

ddvvyOrig <- function(vn, a=solveA(0.9), b=solveB(0.9)) {
  sqrt(0.5)*(ddvvEllipse(vn, a, b))
}

## double check
curve(ddvyOrig,0,0.5)#,ylim=c(-0.2,0.4))
# curve(ddvvyOrig,add=F)
abline(ddvyOrig(0.475)-ddvvyOrig(0.475)*0.475,ddvvyOrig(0.475))
abline(ddvyOrig(0.1)-ddvvyOrig(0.1)*0.1,ddvvyOrig(0.1))
abline(v=c(0.1,0.475), col="grey")
## okay

ellipseOrig <- function(vo, a=solveA(0.9), b=solveB(0.9)) {
  yOrig(vNorm(vo, a, b), a, b)
}

dduellipseOrig <- function(vo, a=solveA(0.9), b=solveB(0.9),
                           dduA=ddusolveA(0.9), dduB=ddusolveB(0.9)) {
  ddvyOrig(vNorm(vo, a, b), a, b)*dduvNorm(vo, a, b, dduA, dduB)+dduyOrig(vNorm(vo, a, b), a, b, dduA, dduB)
}

## double check
ellipseOrigInU <- function(u) ellipseOrig(.5,solveA(u),solveB(u))
dduellipseOrigInU <- function(u) dduellipseOrig(.5,solveA(u),solveB(u),ddusolveA(u),ddusolveB(u))

curve(ellipseOrigInU,0.6,0.85,ylim=c(0.2,0.4))
curve(dduellipseOrigInU,add=T)
abline(ellipseOrigInU(0.8)-dduellipseOrigInU(0.8)*0.8,dduellipseOrigInU(0.8))
abline(ellipseOrigInU(0.65)-dduellipseOrigInU(0.65)*0.65,dduellipseOrigInU(0.65))
abline(v=c(0.8,0.65), col="grey")
## okay

ddvellipseOrig <- function(vo, a=solveA(0.9), b=solveB(0.9)) {
  ddvyOrig(vNorm(vo, a, b), a, b)*ddvvNorm(vo, a, b)
}

## double check
curve(ellipseOrig,0,0.6,ylim=c(0,1))
curve(ddvellipseOrig,add=T,col="green")
abline(ellipseOrig(0.4)-ddvellipseOrig(0.4)*0.4,ddvellipseOrig(0.4))
abline(ellipseOrig(0.1)-ddvellipseOrig(0.1)*0.1,ddvellipseOrig(0.1))
abline(v=c(0.1,0.4), col="grey")
abline(v=strongBorderPoly(foo[1])-weakBorderPoly(foo[1]),col="blue")
abline(h=0)
## okay

ddvuEllipseOrig <- function(vo, a=solveA(0.9), b=solveB(0.9),
                            dduA=ddusolveA(0.9), dduB=ddusolveB(0.9)) {
  cvNorm <- vNorm(vo, a, b)
  (ddvuvNorm(vo, a, b, dduA, dduB)*ddvyOrig(cvNorm, a, b)
   + ddvvNorm(vo, a, b)*(dduvNorm(vo, a, b, dduA, dduB)*ddvvyOrig(cvNorm, a, b)
                         + ddvuyOrig(cvNorm, a, b, dduA, dduB)))
}

## double check
curve(dduellipseOrig,0,0.5,ylim=c(-1,1.5))
curve(ddvuEllipseOrig,add=T)
abline(dduellipseOrig(0.4)-ddvuEllipseOrig(0.4)*0.4,ddvuEllipseOrig(0.4))
abline(dduellipseOrig(0.1)-ddvuEllipseOrig(0.1)*0.1,ddvuEllipseOrig(0.1))
abline(v=c(0.1,0.4), col="grey")
## okay

ddvvEllipseOrig <- function(vo, a=solveA(0.9), b=solveB(0.9)) {
  cvNorm <- vNorm(vo, a, b)
  (ddvvNorm(vo, a, b)^2*ddvvyOrig(cvNorm, a, b)
   +ddvvvNorm(vo, a, b)*ddvyOrig(cvNorm, a, b))
}

## double check
curve(ddvellipseOrig,0,0.5,ylim=c(-1,1.5))
curve(ddvvEllipseOrig,add=T)
abline(ddvellipseOrig(0.4)-ddvvEllipseOrig(0.4)*0.4,ddvvEllipseOrig(0.4))
abline(ddvellipseOrig(0.1)-ddvvEllipseOrig(0.1)*0.1,ddvvEllipseOrig(0.1))
abline(v=c(0.1,0.4), col="grey")
## okay

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