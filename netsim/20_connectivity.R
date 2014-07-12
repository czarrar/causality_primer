#' For now only subject 1
#' Goal here is to compare the accuracy of estimating a connection using NetSim
#' I will compare the regular correlation with a host of partial correlations
#' 
#' Let's read in the data.
#+ read
net   <- as.matrix(read.table("sub01_net.txt")) # reference
ts    <- as.matrix(read.table("sub01_ts.txt"))

#+ compute the correlation
upper <- function(x) x[upper.tri(x)]
ztrue <- function(x, net) upper(x)[upper(net)!=0]
zfalse <- function(x, net) upper(x)[upper(net)==0]
cmat  <- cor(ts)
mean(ztrue(cmat, net))
mean(zfalse(cmat, net))
table(upper(cmat)>0.2, upper(net)!=0)

#' compute the inverse covariance
#+ glasso
library(glasso)
s <- cov(ts)
a <- glasso(s, rho=0.01)
cmat2 <- -cov2cor(a$wi)
#round(cmat2, 2); round(net, 2)
table(me=upper(cmat2)>0.2, ref=upper(net)!=0)

#' below is an example from the package help
#' i'm a little confused as to what passint the a$w and a$wi exactly do
#' is it different then doing that all at once?
#+ glasso-example
set.seed(100)

x<-matrix(rnorm(50*20),ncol=20)
s<- var(x)
a<-glasso(s, rho=.01)
aa<-glasso(s,rho=.02, w.init=a$w, wi.init=a$wi)
a2<-glasso(s,rho=0.2)

#' other. this does a cross-validation...
#+ parcor
library(parcor)
a <- adalasso.net(ts,k=5)
table(me=upper(a$pcor.lasso)>0.15, ref=upper(net)!=0)
table(me=upper(a$pcor.adalasso)>0.15, ref=upper(net)!=0)

a <- pls.net(ts,ncomp=10,k=5)
table(me=upper(a$pcor)>0.16, ref=upper(net)!=0)

a <- ridge.net(ts,k=5)
table(me=upper(a$pcor)>0.16, ref=upper(net)!=0)

#' not sure
#+ quic
library(QUIC)
a <- QUIC(s, rho=0.01)
cmat2 <- -cov2cor(a$X)
table(me=upper(cmat2)>0.2, ref=upper(net)!=0)

#' corpcor
#+ corpcor
library(corpcor)
a <- -cov2cor(invcov.shrink(ts))
table(me=upper(a)>0.15, ref=upper(net)!=0)

a <- -cov2cor(invcor.shrink(ts))
table(me=upper(a)>0.15, ref=upper(net)!=0)

#' might be useful for larger matrices
#+ clime
library(clime)
a <- clime(ts)
cmat2 <- -cov2cor(a$Omegalist[[70]])
table(me=upper(cmat2)>0.15, ref=upper(net)!=0)

#' also good for larger matrices?
#+ scio
library(scio)
a <- scio.cv(s, 0.2) # doesn't work well
cmat2 <- cov2cor(a$w)
table(me=upper(cmat2)>0.16, ref=upper(net)!=0)

#' gives some kind of p-value
#+ lassoscore
library(lassoscore)

a <- glassoscore(ts, 0.2)
table(me=upper(-cov2cor(a$wi))>0, ref=upper(net)!=0)
table(me=upper((a$p.model < 0.0001)), ref=upper(net)!=0)

a <- mbscore(ts, 0.2)
table(me=upper((a$p.model < 0.001)), ref=upper(net)!=0)

#' large matrices?
#+ huge
library(huge)
out.mb = huge(ts)
out.ct = huge(ts, method = "ct")
out.glasso = huge(ts, method = "glasso")

out.select = huge.select(out.mb)
out.select$refit

out.select = huge.select(out.ct)
out.select$refit

out.select = huge.select(out.glasso)
out.select$refit


#' threshold
#+ threshold
# consider using fdrtool
library(fdrtool)
fdrtool(upper(cmat), statistic="correlation")
#qcor0...
# also consider wavestrap
# can use this to determine the threshold
library(wmtsa)
z <- wavBootstrap(ts[,1], n.realization=500)
z2 <- wavBootstrap(ts[,2], n.realization=500)


# tawny
library(tawny)
a1 <- cov_shrink(ts)
-solve(a)[1:5,1:5]

library(BurStFin)
a2 <- var.shrink.eqcor(ts)
a[1:5,1:5]
-solve(a)[1:5,1:5]

a1[1:2,1:2]
a2[1:2,1:2]
all.equal(as.numeric(a1),as.numeric(a2))
