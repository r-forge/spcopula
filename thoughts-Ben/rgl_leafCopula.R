## rgl leafCopula
library(spcopula)
library(rgl)

u=rep((1:200)/201,200)
v=rep((1:200)/201,each=200)

pLeafCop <- data.frame(u, v, p=pCopula(cbind(u,v),spcopula:::leafCopula()))
dLeafCop <- data.frame(u, v, d=pmin(dCopula(cbind(u,v), spcopula:::leafCopula()),50))

surface3d((1:200)/201, (1:200)/201, as.matrix(dLeafCop$d/50,nrow=200),
          col=terrain.colors(51)[round(as.matrix(dLeafCop$d,nrow=200),0)+1])

surface3d((1:200)/201, (1:200)/201, as.matrix(pLeafCop$p,nrow=200),
          col=terrain.colors(51)[round(as.matrix(pLeafCop$p,nrow=200)*50,0)+1])