source("FrontMatter.R")

################
## Chapter 11 ##
################
plotDIRch11 <- paste (plotDIR, "chapter11", "figures", sep="/")
packages(bootstrap)
packages(rv)

## Simulation of the ``nonparametric'' change point model of Qian et al (2003)
## With a normal response and equal variance

chngp <- function(infile)
{ ## infile is a data frame with two columns
  ##  Y and X
    temp <- na.omit(infile)
    yy <- temp$Y
    xx <- temp$X
    mx <- sort(unique(xx))
    m <- length(mx)
    vi <- numeric()
    vi [m] <- sum((yy - mean(yy))^2)
    for(i in 1:(m-1))
            vi[i] <- sum((yy[xx <= mx[i]] - mean(yy[xx <=
                mx[i]]))^2) + sum((yy[xx > mx[i]] - mean(
                yy[xx > mx[i]]))^2)
    thr <- mean(mx[vi == min(vi)])
    return(thr)
}

set.seed(123)
reject <- 0
n.sims <- 50000
for (i in 1:n.sims){
    temp <- data.frame(X=runif(30), Y=rnorm(30))
    split <- chngp(temp)
    if (split==min(temp$X) | split==max(temp$X)) p.value=0.5
    else p.value <- t.test(Y~I(X<split), data=temp, var.equal=T)$p.value
    reject <- reject + (p.value<0.05)/n.sims
}
print(reject)


tikz(paste(plotDIRch11, "chngpms.tex", sep="/"),
     height=3, width=6, standAlone=F)
par(mfrow=c(1,4), mar=c(1.25,0.125,0.5,0.125), mgp=c(1.25,.25,0), tck=0.01)
plot(c(0,1), c(0,1), type="n", xlab="", ylab="", axes=F)
box()
axis(1, at=0.5, label="$\\phi$")
segments(x0=c(0, 0.5), x1=c(0.5,1), y0=c(0.3, 0.75), y1=c(0.3, 0.75))
text(0.05, 0.95, "SF")

plot(c(0,1), c(0,1), type="n", xlab="", ylab="", axes=F)
box()
axis(1, at=0.5, label="$\\phi$")
segments(x0=c(0, 0.5), x1=c(0.5,1), y0=c(0.5, 0.7), y1=c(0.7, 0.25))
text(0.05, 0.95, "HS")

plot(c(0,1), c(0, 1), type="n", xlab="", ylab="", axes=F)
box()
axis(1, at=0.5, label="$\\phi$")
segments(x0=c(0, 0.5), x1=c(0.5,1), y0=c(0.25, 0.7), y1=c(0.5, 0.3))
text(0.05, 0.95, "dBS")

plot(c(0,1), c(0, 1), type="n", xlab="", ylab="", axes=F)
box()
axis(1, at=0.5, label="$\\phi$")
curve(1/(1+exp(10-20*x)), add=T)
text(0.05, 0.95, "SM")

dev.off()


### Function for calculating IVs along a gradient
###  Corrected a coding error in B&K for handling ties in
###  gradient values

T1 <- function(grad, spAbun, minsplit=5, Log=F){
## calculating IVs using species abundance along a gradient
    if (Log) spAbun <- log10(spAbun+1)
    ## in some cases, B&K used log10(abundance + 1)
    oo <- order(grad)
    spAbun <- spAbun[oo] ## sorting
    grad <- grad[oo]
    splits <- unique(grad)
    n <- length(splits)
    nsp <- n-2*minsplit+1 ## number of potential splits
    b <- as.numeric(spAbun>0)
    IV <- numeric()
    for (i in 1:nsp){
        x1 <- sum(spAbun[grad<=splits[minsplit+i-1]])
        x2 <- sum(spAbun[grad>=splits[minsplit+i]])
        x1 <- x1/(minsplit+i-1)
        x2 <- x2/(n-minsplit-i+1)
        tabun <- x1+x2
        if(tabun==0)tabun <- 1
        Ra1 <- x1/(tabun)
        Ra2 <- 1-Ra1
        Rf1 <- sum(b[grad<=splits[minsplit+i-1]])
        Rf2 <- sum(b[grad>=splits[minsplit+i]])
        Rf1 <- Rf1/length(b[grad<=splits[minsplit+i-1]])
        Rf2 <- Rf2/length(b[grad>=splits[minsplit+i]])
        IV[i] <- 100*max(c(Ra1*Rf1, Ra2*Rf2))
    }
    new.grad <- splits[minsplit:(n-minsplit)]
    split <- seq(1,length(new.grad))[IV==max(IV)]
    if (min(split)==1 | max(split)==length(new.grad)) no.chngp <- 1
    else no.chngp <- 0
    if (max(split)==length(new.grad))
        split<- new.grad[max(split)]
    else
        split <- mean(new.grad[c(split, split+1)])
    return (list(IVs=data.frame(splits=splits[minsplit:(n-minsplit)],
                     IndVal=IV),
                 split=c(split, no.chngp)))
}

### Function for permutation test for a specific split
### See notes above
perm <- function(numperm = 250, grad, spAbun, minsplit=5, Log=FALSE){
  oo <- order(grad)
  grad<- grad[oo]
  spAbun <- spAbun[oo]
  b <- as.numeric(spAbun >0)
  if (Log) spAbun <- log10(spAbun+1)
  split <- T1(grad, spAbun, minsplit, Log=Log)
  split.obs<-max(split$IVs$IndVal)
  cp <- split$split[1]
  n1 <- sum(grad<=cp)
  n2 <- sum(grad> cp)
  split.perm <- numeric()
  grad.count <- table(grad)
  grad.order <- as.numeric(ordered(grad))
  for (i in 1:numperm){
    perms <- sample(rep(1:2, c(n1, n2)), size=max(grad.order))
    while (sum(perms==1)==0 | sum(perms==2)==0)
      perms <- sample(rep(1:2, c(n1, n2)), size=max(grad.order))
    grad.od <- rep(perms, grad.count)
    i1=which(grad.od==1)
    i2=which(grad.od==2)
    x1 <- sum(spAbun[i1])
    x2 <- sum(spAbun[i2])
    x1 <- x1/length(unique(grad[i1]))
    x2 <- x2/length(unique(grad[i2]))
    tatl <- x1+x2
    if (tatl == 0) tatl<-1
    Ra1 <- x1/(tatl)
    Ra2 <- 1-Ra1
    Rf1 <- sum(b[i1])
    Rf2 <- sum(b[i2])
    Rf1 <- Rf1/length(i1)
    Rf2 <- Rf2/length(i2)
    split.perm[i] <- 100*max(c(Ra1*Rf1, Ra2*Rf2))
  }
  return(c(p.value=mean(split.perm >= split.obs),
              z=(split.obs-mean(split.perm))/sd(split.perm),
              perm.mu = mean(split.perm), perm.sd = sd(split.perm)))
}

## simulation of a permutation with a given number of permutations.

perm.sim <- function(numperm = 250, grad, spAbun, minsplit=5, Log=F){
  oo <- order(grad)
  grad.order <- as.numeric(ordered(grad))
  grad.unique <- sort(unique(grad))
  grad<- grad[oo]
  spAbun <- spAbun[oo]
  b <- as.numeric(spAbun >0)
  if (Log) spAbun <- log10(spAbun+1)
  splitperm <- numeric()
  k <- 0
  permsim.out <- matrix(0, ncol=6, nrow=max(grad.order)-2*minsplit+1)
  for (i in minsplit:(max(grad.order)-minsplit)){
    splitobs.cp<-0.5*(grad.unique[i]+grad.unique[i+1])
    perms <- c(rep(1, i), rep(2, length(grad)-i))
    n1 <- i
    n2 <- length(grad)-i
    split.perm <- numeric()
    grad.count <- table(grad)
    grad.od <- rep(perms, grad.count)
    i1 <- which(grad.od==1)
    i2 <- which(grad.od==2)
    x1 <- sum(spAbun[i1])
    x2 <- sum(spAbun[i2])
    x1 <- x1/length(unique(grad[i1]))
    x2 <- x2/length(unique(grad[i2]))
    tatl <- x1+x2
    if (tatl == 0) tatl<-1
    Ra1 <- x1/(tatl)
    Ra2 <- 1-Ra1
    Rf1 <- sum(b[i1])
    Rf2 <- sum(b[i2])
    Rf1 <- Rf1/length(i1)
    Rf2 <- Rf2/length(i2)
    splitobs.IndVal <- 100*max(c(Ra1*Rf1, Ra2*Rf2))
    for (j in 1:numperm){
      perms <- sample(rep(1:2, c(n1, n2)), size=max(grad.order))
      while (sum(perms==1)==0 | sum(perms==2)==0)
        perms <- sample(rep(1:2, c(n1, n2)), size=max(grad.order))
      grad.od <- rep(perms, grad.count)
      i1=which(grad.od==1)
      i2=which(grad.od==2)
      x1 <- sum(spAbun[i1])
      x2 <- sum(spAbun[i2])
      x1 <- x1/length(unique(grad[i1]))
      x2 <- x2/length(unique(grad[i2]))
      tatl <- x1+x2
      if (tatl == 0) tatl<-1
      Ra1 <- x1/(tatl)
      Ra2 <- 1-Ra1
      Rf1 <- sum(b[i1])
      Rf2 <- sum(b[i2])
      Rf1 <- Rf1/length(i1)
      Rf2 <- Rf2/length(i2)
      splitperm[j] <- 100*max(c(Ra1*Rf1, Ra2*Rf2))
    }
    k<- k+1
    permsim.out[k,] <- c(cp = splitobs.cp,
                         IndVal=splitobs.IndVal,
                         p.value=mean(splitperm >= splitobs.IndVal),
                         z=(splitobs.IndVal-mean(splitperm))/sd(splitperm),
                         perm.mu = mean(splitperm),
                         perm.sd = sd(splitperm))
  }
  colnames(permsim.out) <- c("cp","IndVal","pval","z","mu","sd")
  return(permsim.out)
}

## function to draw random response variable values

drawY <- function(x, phi1, beta0, beta1, sigma, model=M1, Log=T, Plot=F, alpha0=0, alpha1=0, ...){
  mu <- model(x, b0=beta0, b1=beta1, pt=Plot, a0=alpha0, a1=alpha1, ... )
  n <- length(x)
  if (Log) return(list(x=x, y=ifelse(mu>0, exp(rnorm(sum(mu>0), mu, sigma)), 0),
              z=exp(mu+0.5*sigma^2)))
  else return(list(x=x, y=ifelse(mu>0, rnorm(sum(mu>0), mu, sigma), 0),
              z=mu))
}


## the step function
M1 <-  function(x=runif(100), phi1=0.5, b0, b1, pt=T, a0=0, a1=0,...){
  if (pt) plot(x, b0+ (x>phi1)*b1)
  return(b0+(x>phi1) * b1)
}

## linear model
M2 <- function(x=runif(100), phi1=0,  b0, b1,pt=T, a0=0, a1=0, ...){
  y <- b0 + b1 * x
  if (pt) plot(x, ifelse(y>0, y, 0))
  return (ifelse(y>0, y, 0))
}

## piecewise linear model
M3 <- function(x=runif(100), phi1=0.5,  b0, b1, pt=T, a0=0, a1=0, ...){
  y <- ifelse(x <=phi1, b0, b0+b1*(x-phi1))
  if (pt) plot(x, ifelse (y>0,y,0))
  return(ifelse (y>0,y,0))
}

## M4 SM model
M4 <- function (x=runif(100), b0=-10, b1=20, a=20, b=-10, pt=T,a0=0, a1=0, ...){
    y <- a+b*invlogit(b0+b1*x)
    if(pt) plot(x, y)
  return(ifelse(y>0, y, 0))
}

## dBS model
M5 <- function(x=runif(100), phi1=0.5, b0, b1, pt=T, a0=0, a1=0, ...){
    y <- ifelse(x<=phi1, b0+b1*(x-phi1), a0+a1*(x-phi1))
    if(pt) plot(x, y)
    return(ifelse(y>0, y, 0))
}
## Simulation without error
## SF model
tikz(file=paste(plotDIRch11, "titansim1.tex", sep="/"), height=5,
     width=5, standAlone=T)
par(mfrow=c(2,2), mar=c(1.5, 1, 0.5,1), mgp=c(0.5,0.125,0), las=1, tck=0.01)
## step function (both sides > 0)
temp1 <- drawY (x=seq(0,1,,101), phi1=0.5, beta0=log(45),
                     beta1=-log(10), sigma=0.)
stepT1 <-  T1(temp1$x, ifelse(temp1$y<0,0,temp1$y))
plot(temp1, xlab="gradient", ylab=" ",
     ylim=range(c(temp1$y, stepT1$IVs$IndVal)),
     axes=F)
box()
abline(v=0.5, col="gray")
lines(stepT1$IVs$splits, stepT1$IVs$IndVal)

## HS
temp1 <- drawY (x=seq(0,1,,101), phi1=0.5, beta0=log(40),
                beta1=-log(40), sigma=0., Log=F, model=M3)
hockeyT1 <-  T1(temp1$x, ifelse(temp1$y<0,0,temp1$y))
plot(temp1$x, temp1$y*10, xlab="gradient", ylab=" ",
     ylim=range(c(temp1$y*10, hockeyT1$IVs$IndVal/1.25)),
     axes=F)
box()
abline(v=0.5, col="gray")
lines(hockeyT1$IVs$splits, hockeyT1$IVs$IndVal/1.25)

## dBS model
temp1 <- drawY (x=seq(0,1, ,101), beta0=25, beta1=-2, alpha0=22, alpha1=-12, sigma=0,
                Log=F, model=M5)
dsbT1 <-  T1(temp1$x, ifelse(temp1$y<0,0,temp1$y))
plot(temp1$x, temp1$y, xlab="gradient", ylab=" ",
     ylim=range(c(temp1$y, dsbT1$IVs$IndVal/1.75)),
          axes=F)
box()
abline(v=0.5, col="gray")
lines(dsbT1$IVs$splits, dsbT1$IVs$IndVal/1.75)

## SM
temp1 <- drawY (x=seq(0,1,,101), beta0=-10,
                beta1=20, sigma=0., Log=F, model=M4)
sigmoT1 <-  T1(temp1$x, ifelse(temp1$y<0,0,temp1$y))
plot(temp1$x, temp1$y, xlab="gradient", ylab=" ",
     ylim=range(c(temp1$y, sigmoT1$IVs$IndVal/2.25)),
          axes=F)
box()
abline(v=0.5, col="gray")
lines(sigmoT1$IVs$splits, sigmoT1$IVs$IndVal/2.25)

dev.off()

## Type I error probability 
typeIsim <- function(n.sims=2500, ns=21){
    x <- seq(0,1,,ns)
    fr.rej <- 0
    fr.z <- 0
    for (i in 1:n.sims){
        if (i%%100==0) print(paste(i, "of", n.sims, fr.rej, fr.z))
        y <- rpois(length(x), 20)
        fr.rej <- fr.rej +
            (perm(numperm = 1000, grad=x, spAbun=y)[1]<0.05)/n.sims
        fr.z <- fr.z +
            (abs(perm(numperm = 1000, grad=x, spAbun=y)[2])> 1.96)/n.sims
    }
    return(c(fr.rej, fr.z))
}

## 1. # sampling points = 21
ns <- 15
prI.15 <- typeIsim(n.sims=5000, ns=ns)

ns <- 25
prI.25 <- typeIsim(n.sims=5000, ns=ns)
## 2. # sampling points = 51
ns <- 51
prI.51 <- typeIsim(n.sims=5000, ns=ns)
## 3. # sampling points =101
ns <- 101
prI.101 <- typeIsim(n.sims=5000, ns=ns)
## 3. # sampling points =101
ns <- 201
prI.201 <- typeIsim(n.sims=5000, ns=ns)
save(prI.15, prI.25, prI.51, prI.101, prI.201, file="prI.RData")
load("prI.RData")


## the z-score
## 1 IndVal under the null

xg <- T1(seq(0,1,,101), rpois(101, 20))$IVs[,1]
indvals <- NULL
for (i in 1:1000)
   indvals <- rbind(indvals, T1(seq(0,1,,101), rpois(101, 20))$IVs[,2])

indvals.rv <- rvsims(indvals)

tikz(paste(plotDIRch11, "indvalsims.tex", sep="/"), height=3, width=4,
     standAlone=F)
par(mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(indvals.rv, axes=F, ylim=c(49.5,57), xlab="gradient", ylab="$IndVal$")
axis(1, at=1:length(indvals.rv), labels=xg)
axis(2)
box()
dev.off()

    ## z-score

sims.fun <- function(grd = seq(0,1,,101), txnAb=20, n.sims=1000){
    zz <- mm <- sig <- indvals <- NULL
    permsim <- ifelse(n.sims>=1000, 1000, n.sims)
    print("starting ...")
    for (i in 1:n.sims){
        if(i%%10 == 0) print(paste("iteration ", i, " of ", n.sims))
        temp <- as.data.frame(perm.sim(numperm = permsim,
                                       grad = grd,
                                       spAbun = rpois(length(grd), txnAb)))
        zz <- rbind(zz, temp$z)
        mm <- rbind(mm, temp$mu)
        sig <- rbind(sig, temp$sd)
        indvals <- rbind(indvals, temp$IndVal)
    }
    indvals.rv <- rvsims(indvals)
    zz.rv <- rvsims(zz)
    mm.rv <- rvsims(mm)
    sig.rv <- rvsims(sig)
    return(list(z=zz.rv, mu=mm.rv, sig=sig.rv, indV=indvals.rv))
}

gd <- seq(0,1,,101)

## the null model
##sims.NULL <- sims.fun(n.sims=5000)
##save(sims.NULL, file="simNull.rda")
load("simNull.rda")

indV <- sims.NULL$indV
z <- sims.NULL$z
mu <- sims.NULL$mu
sig <- sims.NULL$sig

tikz(paste(plotDIRch11, "indVz0.tex", sep="/"),
     width=4, height=2, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(rvgamma(1, rep(20, 101)), axes=F, , xlab="gradient", ylab="Abundance",
     ylim=c(10, 40))
axis(1, at=1:length(gd), labels=gd)
axis(2)
box()
text(5, 35, "(a)")

plot(indV, axes=F, xlab="gradient", ylab="$IV$", ylim=c(50,60))
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
text(5, 58, "(b)")

plot(z, axes=F, xlab="gradient", ylab="$z$", ylim=c(-1.5, 4))
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
text(5, 3, "(c)")

dev.off()

tikz(paste(plotDIRch11, "musigmaNull.tex", sep="/"), height=1.75, width=3.75,
     standAlone=F)
par(mfrow=c(1,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(mu, axes=F, xlab="gradient", ylab="$\\mu$")
axis(1, at=1:length(indV), labels=xg)
axis(2, at=51:52)
box()
text(6, 52.5, "(a)")

plot(sig, axes=F, xlab="gradient", ylab="$\\sigma$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
text(6, 2, "(b)")
dev.off()

tikz(paste(plotDIRch11, "indvalsims0.tex", sep="/"), height=3, width=4,
     standAlone=F)
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(indV, axes=F, ylim=c(50, 56), xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(mu, axes=F, xlab="gradient", ylab="$\\mu$")
axis(1, at=1:length(indV), labels=xg)
axis(2, at=51:53)
box()

plot(sig, axes=F, xlab="gradient", ylab="$\\sigma$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

dev.off()

tikz(paste(plotDIRch11, "indvalsims0.tex", sep="/"), height=3, width=4,
     standAlone=F)
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(indV, axes=F, ylim=c(50, 56), xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(mu, axes=F, xlab="gradient", ylab="$\\mu$")
axis(1, at=1:length(indV), labels=xg)
axis(2, at=51:53)
box()

plot(sig, axes=F, xlab="gradient", ylab="$\\sigma$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

dev.off()

tikz(paste(plotDIRch11, "indvalsims0.tex", sep="/"), height=3, width=4,
     standAlone=F)
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(indV, axes=F, ylim=c(50, 56), xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(mu, axes=F, xlab="gradient", ylab="$\\mu$")
axis(1, at=1:length(indV), labels=xg)
axis(2, at=51:53)
box()

plot(sig, axes=F, xlab="gradient", ylab="$\\sigma$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

dev.off()

tikz(paste(plotDIRch11, "indvalsims0.tex", sep="/"), height=3, width=4,
     standAlone=F)
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(indV, axes=F, ylim=c(50, 56), xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(mu, axes=F, xlab="gradient", ylab="$\\mu$")
axis(1, at=1:length(indV), labels=xg)
axis(2, at=51:53)
box()

plot(sig, axes=F, xlab="gradient", ylab="$\\sigma$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

dev.off()

## Power

powersim <- function(x=seq(0,1,,101), spAbun=20+(gd>0.5)*10, n.sims=2500){
    if (length(x) != length(spAbun)) stop("Check Data!")
    n.perms <- ifelse(n.sims > 1000, 1000, n.sims)
    fr.rej <- 0
    fr.z <- 0
    for (i in 1:n.sims){
        if (i%%100==0) print(paste(i, "of", n.sims, fr.rej, fr.z))
        y <- rpois(length(x), spAbun)
        fr.rej <- fr.rej +
            (perm(numperm = n.perms, grad=x, spAbun=y)[1]<0.05)/n.sims
        fr.z <- fr.z +
            (abs(perm(numperm = n.perms, grad=x, spAbun=y)[2])> 1.96)/n.sims
    }
    return(c(fr.rej, fr.z))
}

gd <- seq(0,1,,101)

## SF model
power.SF <- powersim(spAbun= 20+(gd>0.5)*10, n.sims=5000)
save(power.SF, file="powerSF.rda")
load("powerSF.rda")

## HS model
power.HS <- powersim(spAbun= 20+(gd>0.5)*(20*(gd-0.5)), n.sims=5000)
save(power.HS, file="powerHS.rda")
load("powerHS.rda")

## SM model
power.SM <- powersim(spAbun= 20+10*invlogit(-5+10*gd), n.sims=5000)
save(power.SM, file="powerSM.rda")
load("powerSM.rda")

## LM model
power.LM <- powersim(spAbun= 20+10*gd, n.sims=5000)
save(power.LM, file="powerLM.rda")
load("powerSM.rda")

## alternative 1: linear
##sims.ALT1 <- sims.fun(txnAb= 20+10*gd, n.sims=5000)
##save(sims.ALT1, file="simAlt1.rda")
load("simAlt1.rda")

indV <- sims.ALT1$indV
z <- sims.ALT1$z
mu <- sims.ALT1$mu
sig <- sims.ALT1$sig

tikz(paste(plotDIRch11, "indVz1.tex", sep="/"),
     width=4, height=2, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(rvgamma(1, 20+10*gd), axes=F,  xlab="gradient", ylab="Abundance",
     ylim=c(10,50))
axis(1, at=1:length(gd), labels=gd)
axis(2)
box()
text(5, 45, "(a)")

plot(indV, axes=F, xlab="gradient", ylab="$IV$", ylim=c(50,65))
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
text(5, 63, "(b)")

plot(z, axes=F, xlab="gradient", ylab="$z$", ylim=c(-1,10))
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
text(5, 8.5, "(c)")

dev.off()


tikz(paste(plotDIRch11, "indvalsims1.tex", sep="/"), height=3, width=4,
     standAlone=T)
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(indV, ylim=c(50,61), axes=F, xlab="gradient", ylab="$IV$")
axis(1, at=1:length(gd), labels=gd)
axis(2)
box()

plot(mu, axes=F, xlab="gradient", ylab="$\\mu$")
axis(1, at=1:length(indV), labels=xg)
axis(2, at=51:53)
box()

plot(sig, axes=F, xlab="gradient", ylab="$\\sigma$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

dev.off()

## Alternative 1, reduced sample size
xg1 <- T1(seq(0,1,,21), rpois(21, 20))$IVs[,1]
gd1 <- seq(0,1,,21)
#sims.ALT11 <- sims.fun(grd=gd1, txnAb= 20+10*gd1, n.sims=5000)
#save(sims.ALT11, file="simAlt11.rda")
load("simAlt11.rda")

indV <- sims.ALT11$indV
z <- sims.ALT11$z
mu <- sims.ALT11$mu
sig <- sims.ALT11$sig

tikz(paste(plotDIRch11, "indVz11.tex", sep="/"),
     width=4, height=2, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(rvgamma(1, 20+10*gd1), axes=F, xlab="gradient", ylab="Abundance",
     ylim=c(10,50))
axis(1, at=1:length(gd1), labels=gd1)
axis(2)
box()
text(2, 45, "(a)")

plot(indV, axes=F, xlab="gradient", ylab="$IV$", ylim=c(50,65))
axis(1, at=1:length(indV), labels=xg1)
axis(2)
box()
text(5, 63, "(b)")

plot(z, axes=F, xlab="gradient", ylab="$z$", ylim=c(-1,10))
axis(1, at=1:length(indV), labels=xg1)
axis(2)
box()
text(5, 8.5, "(c)")

dev.off()


tikz(paste(plotDIRch11, "indvalsims11.tex", sep="/"), height=3, width=4,
     standAlone=F)
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(indV, ylim=c(50,61), axes=F, xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(mu, axes=F, xlab="gradient", ylab="$\\mu$")
axis(1, at=1:length(indV), labels=xg)
axis(2, at=51:53)
box()

plot(sig, axes=F, xlab="gradient", ylab="$\\sigma$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

dev.off()


## alternative 2: HS model
##sims.ALT2 <- sims.fun(txnAb= ifelse(gd<=0.5, 20, 20+20*(gd-0.5)),
##                     n.sims=5000)
##save(sims.ALT2, file="simAlt2.rda")
load("simAlt2.rda")

indV <- sims.ALT2$indV
z <- sims.ALT2$z
mu <- sims.ALT2$mu
sig <- sims.ALT2$sig

tikz(paste(plotDIRch11, "indVz2.tex", sep="/"),
     width=4, height=2, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(rvgamma(1, ifelse(gd<=0.5, 20, 20+20*(gd-0.5))), axes=F,
     xlab="gradient", ylab="Abundance")
axis(1, at=1:length(gd), labels=gd)
axis(2)
box()
text(5, 52.5, "(a)")

plot(indV, axes=F, xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
text(5, 62, "(b)")

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
abline(v=46.5, col="gray")
text(5, 9, "(c)")

dev.off()


tikz(paste(plotDIRch11, "indvalsims2.tex", sep="/"), height=3, width=4,
     standAlone=F)
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(indV, ylim=c(50,61), axes=F, xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(mu, axes=F, xlab="gradient", ylab="$\\mu$")
axis(1, at=1:length(indV), labels=xg)
axis(2, at=51:53)
box()

plot(sig, axes=F, xlab="gradient", ylab="$\\sigma$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
abline(v=46.5, col="gray")

dev.off()

## alternative 3: SF model
##sims.ALT3 <- sims.fun(txnAb=ifelse(gd<=0.5, 20, 20+10), n.sims=5000)
##save(sims.ALT3, file="simAlt3.rda")
load("simAlt3.rda")

indV <- sims.ALT3$indV
z <- sims.ALT3$z
mu <- sims.ALT3$mu
sig <- sims.ALT3$sig

tikz(paste(plotDIRch11, "indVz3.tex", sep="/"),
     width=4, height=2, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(rvgamma(1, ifelse(gd<=0.5, 20, 20+10)), axes=F,
     xlab="gradient", ylab="Abundance")
axis(1, at=1:length(gd), labels=gd)
axis(2)
box()
text(5, 52, "(a)")

plot(indV, axes=F, xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
text(5, 63, "(b)")

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
abline(v=46.5, col="gray")
text(5, 11, "(c)")

dev.off()

tikz(paste(plotDIRch11, "indvalsims3.tex", sep="/"), height=3, width=4,
     standAlone=F)
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(indV, axes=F, xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(mu, axes=F, xlab="gradient", ylab="$\\mu$")
axis(1, at=1:length(indV), labels=xg)
axis(2, at=51:53)
box()

plot(sig, axes=F, xlab="gradient", ylab="$\\sigma$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
abline(v=46.5, col="gray")

dev.off()


## alternative 4: SM model
##sims.ALT4 <- sims.fun(txnAb=20+10*invlogit(-5+10*gd), n.sims=5000)
##save(sims.ALT4, file="simAlt4.rda")
load("simAlt4.rda")

indV <- sims.ALT4$indV
z <- sims.ALT4$z
mu <- sims.ALT4$mu
sig <- sims.ALT4$sig

tikz(paste(plotDIRch11, "indVz4.tex", sep="/"),
     width=4, height=2, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(rvgamma(1, 20+10*invlogit(-5+10*gd)), axes=F,
     xlab="gradient", ylab="Abundance")
axis(1, at=1:length(gd), labels=gd)
axis(2)
box()
text(5, 55, "(a)")

plot(indV, axes=F, xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
text(5, 62, "(b)")

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
abline(v=46.5, col="gray")
text(5, 10, "(c)")
dev.off()

tikz(paste(plotDIRch11, "indvalsims4.tex", sep="/"), height=3, width=4,
     standAlone=F)
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(indV, axes=F, xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(mu, axes=F, xlab="gradient", ylab="$\\mu$")
axis(1, at=1:length(indV), labels=xg)
axis(2, at=51:53)
box()

plot(sig, axes=F, xlab="gradient", ylab="$\\sigma$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
abline(v=46.5, col="gray")
dev.off()

## alternative 5: SM model
##sims.ALT5 <- sims.fun(txnAb=20+10*invlogit(-10+20*gd), n.sims=5000)
##save(sims.ALT5, file="simAlt5.rda")
load("simAlt5.rda")

indV <- sims.ALT5$indV
z <- sims.ALT5$z
mu <- sims.ALT5$mu
sig <- sims.ALT5$sig

tikz(paste(plotDIRch11, "indVz5.tex", sep="/"),
     width=4, height=2, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(rvgamma(1, 20+10*invlogit(-10+20*gd)), axes=F,
     xlab="gradient", ylab="Abundance")
axis(1, at=1:length(gd), labels=gd)
axis(2)
box()
text(5, 52, "(a)")

plot(indV, axes=F, xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
text(5, 63, "(b)")

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
abline(v=46.5, col="gray")
text(5, 11, "(c)")
dev.off()

tikz(paste(plotDIRch11, "indvalsims5.tex", sep="/"), height=3, width=4,
     standAlone=T)
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(indV, axes=F, xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(mu, axes=F, xlab="gradient", ylab="$\\mu$")
axis(1, at=1:length(indV), labels=xg)
axis(2, at=51:53)
box()

plot(sig, axes=F, xlab="gradient", ylab="$\\sigma$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
abline(v=46.5, col="gray")
dev.off()

## alternative 6: SM model
##sims.ALT6 <- sims.fun(txnAb=20+10*invlogit(-20+40*gd), n.sims=5000)
##save(sims.ALT6, file="simAlt6.rda")
load("simAlt6.rda")

indV <- sims.ALT6$indV
z <- sims.ALT6$z
mu <- sims.ALT6$mu
sig <- sims.ALT6$sig

tikz(paste(plotDIRch11, "indVz6.tex", sep="/"),
     width=4, height=2, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(rvgamma(1, 20+10*invlogit(-20+40*gd)), axes=F,
     xlab="gradient", ylab="Abundance")
axis(1, at=1:length(gd), labels=gd)
axis(2)
box()
text(5, 57, '(a)')

plot(indV, axes=F, xlab="gradient", ylab="$IV$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
text(5, 65, '(b)')

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
abline(v=46.5, col="gray")
text(5, 12, '(c)')
dev.off()

tikz(paste(plotDIRch11, "indvalsims6.tex", sep="/"), height=3, width=4,
     standAlone=F)
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(indV, axes=F, xlab="gradient", ylab="$IndVal$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(mu, axes=F, xlab="gradient", ylab="$\\mu$")
axis(1, at=1:length(indV), labels=xg)
axis(2, at=51:53)
box()

plot(sig, axes=F, xlab="gradient", ylab="$\\sigma$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()

plot(z, axes=F, xlab="gradient", ylab="$z$")
axis(1, at=1:length(indV), labels=xg)
axis(2)
box()
abline(v=46.5, col="gray")
dev.off()

###################################################




size <- length(temp1$x)
perm(2500, temp1$x, ifelse(temp1$y<0,0,temp1$y))



perm(2500, temp1$x, ifelse(temp1$y<0,0,temp1$y))

size <- length(temp1$x)
#results1 <- bootstrap(x=1:size, nboot=2000,
#                      theta=function(x, infile) {
#                        temp <- T1(infile[x,1], infile[x,2])$split
#                      },
#                      infile=data.frame(x=temp1$x, y=temp1$y))
#CI1 <- as.vector(round(quantile(results1$thetastar[1,],
#                               prob=c(0.025,0.975), na.rm=T), 2))#

plot(temp2 <- drawY (x=seq(0,1,,100), phi1=0.5, phi2=0., beta0=log(25),
                    beta1=-log(5), sigma=0.5))
perm(2500, temp2$x, ifelse(temp2$y<0,0,temp2$y))

stepT2 <-  T1(temp2$x, ifelse(temp2$y<0,0,temp2$y))
size <- length(temp2$x)
results2 <- bootstrap(x=1:size, nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])$split
                      },
                      infile=data.frame(x=temp2$x, y=temp2$y))
CI2 <- as.vector(round(quantile(results2$thetastar[1,],
                               prob=c(0.025,0.975), na.rm=T), 2))

plot(temp3 <- drawY (x=seq(0,1,,100), phi1=0.5, phi2=0., beta0=log(25),
                    beta1=-log(5), sigma=1))
perm(2500, temp3$x, ifelse(temp3$y<0,0,temp3$y))

stepT3 <-  T1(temp3$x, ifelse(temp3$y<0,0,temp3$y))
size <- length(temp3$x)
results3 <- bootstrap(x=1:size, nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])$split
                      },
                      infile=data.frame(x=temp3$x, y=temp3$y))
CI3 <- as.vector(round(quantile(results3$thetastar[1,],
                               prob=c(0.025,0.975), na.rm=T), 2))

#pdf(file=paste(plotDIR, "stepSims.pdf", sep="/"),height=6.5, width=6.5)
tikz(file=paste(plotDIR, "stepSimsSA.tex", sep="/"),height=6.5, width=6.5, standAlone=T)
par(mfrow=c(3,3), mar=c(3,3,0.5,0.25), mgp=c(1.5,0.25,0), tck=0.02,las=1)
plot(temp1$x,temp1$y, xlab="Gradient", ylab="Abundance")
text(x=0.8, y=floor(max(temp1$y))/1.2, "$\\sigma=0.0$")
#text(x=0.8, y=floor(max(temp1$y))/1.2, expression(sigma == 0.0))
lines(x=c(0,0.5),y=rep(25, 2))
lines(x=c(0.5,1),y=rep(exp(log(25)-log(5)), 2))
plot(stepT1$IVs[,1], stepT1$IVs[,2], xlab="Gradient", ylab="IndVal", xlim=c(0,1))
hist(results1$thetastar[1,], xlim=c(0,1), prob=T, main="", xlab="Gradient")
#text(x=0.2,y=8, paste(CI1[1], CI1[2]))
abline(v=0.5, col="red", lwd=2.5)

plot(temp2$x,temp2$y, xlab="Gradient", ylab="Abundance")
text(x=0.8, y=floor(max(temp2$y))/1.2, "$\\sigma=0.5$")
#text(x=0.8, y=floor(max(temp2$y))/1.2, expression(sigma == 0.5))
lines(x=c(0,0.5),y=rep(25, 2))
lines(x=c(0.5,1),y=rep(exp(log(25)-log(5)), 2))
plot(stepT2$IVs[,1], stepT2$IVs[,2], xlab="Gradient", ylab="IndVal", xlim=c(0,1))
hist(results2$thetastar[1,], xlim=c(0,1), prob=T, main="", xlab="Gradient")
# text(x=0.2,y=4, paste(CI2[1], CI2[2]))
abline(v=0.5, col="red", lwd=2.5)

plot(temp3$x,temp3$y, xlab="Gradient", ylab="Abundance")
text(x=0.8, y=floor(max(temp3$y))/1.2, "$\\sigma=1.0$")
#text(x=0.8, y=floor(max(temp3$y))/1.2, expression(sigma == 1.0))
lines(x=c(0,0.5),y=rep(25, 2))
lines(x=c(0.5,1),y=rep(exp(log(25)-log(5)), 2))
plot(stepT3$IVs[,1], stepT3$IVs[,2], xlab="Gradient", ylab="IndVal", xlim=c(0,1))
hist(results3$thetastar[1,], xlim=c(0,1), prob=T, main="", xlab="Gradient")
#text(x=0.2,y=4, paste(CI3[1], CI3[2]))
abline(v=0.5, col="red", lwd=2.5)
dev.off()


## find the mode of the change point distribution:
binned.split<-hist(results1$thetastar, breaks=seq(0,1,,50), plot=F)
split.mode <- binned.split$mids[binned.split$counts==max(binned.split$counts)]


## linear function
plot(temp1 <- drawY (x=seq(0,1,,100), phi1=0, phi2=0., beta0=log(5),
                    beta1=3, sigma=0., model=M2, Log=F))

perm(2500, temp1$x, ifelse(temp1$y<0,0,temp1$y))
linearT1 <-  T1(temp1$x, ifelse(temp1$y<0,0,temp1$y))
size <- length(temp1$x)
results1 <- bootstrap(x=1:size, nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])$split
                      },
                      infile=data.frame(x=temp1$x, y=temp1$y))

plot(temp2 <- drawY (x=seq(0,1,,100), phi1=0, phi2=0., beta0=log(5),
                    beta1=3, sigma=0.5, model=M2, Log=F))

perm(2500, temp2$x, ifelse(temp2$y<0,0,temp2$y))
linearT2 <-  T1(temp2$x, ifelse(temp2$y<0,0,temp2$y))
size <- length(temp2$x)
results2 <- bootstrap(x=1:size, nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])$split
                      },
                      infile=data.frame(x=temp2$x, y=temp2$y))

plot(temp3 <- drawY (x=seq(0,1,,100), phi1=0, phi2=0., beta0=log(5),
                    beta1=3, sigma=1, model=M2, Log=F))

perm(2500, temp3$x, ifelse(temp3$y<0,0,temp3$y))
linearT3 <-  T1(temp3$x, ifelse(temp3$y<0,0,temp3$y))
size <- length(temp3$x)
results3 <- bootstrap(x=1:size, nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])$split
                      },
                      infile=data.frame(x=temp3$x, y=temp3$y))


#pdf(file=paste(plotDIR, "linearSims.pdf", sep="/"),height=6.5, width=6.5)
tikz(file=paste(plotDIR, "linearSimsSA.tex", sep="/"),height=6.5, width=6.5, standAlone=T)
par(mfrow=c(3,3), mar=c(3,3,0.5,0.25), mgp=c(1.5,0.25,0), tck=0.02,las=1)
plot(temp1$x,temp1$y, xlab="Gradient", ylab="Abundance", ylim=c(0,6.5))
text(x=0.2, y=6, "$\\sigma=0.0$")
#text(x=0.2, y=6, expression(sigma==0.0))
abline(log(5), 3)
plot(linearT1$IVs[,1], linearT1$IVs[,2], xlab="Gradient",
     ylab="IndVal", xlim=c(0,1))
hist(results1$thetastar[1,], xlim=c(0,1), prob=T, main="", xlab="Gradient")

plot(temp2$x,temp2$y, xlab="Gradient", ylab="Abundance", ylim=c(0,6.5))
abline(log(5), 3)
text(x=0.2, y=6, "$\\sigma=0.5$")
#text(x=0.2, y=6, expression(sigma==0.5))
plot(linearT2$IVs[,1], linearT2$IVs[,2], xlab="Gradient",
     ylab="IndVal", xlim=c(0,1))
hist(results2$thetastar[1,], xlim=c(0,1), prob=T, main="", xlab="Gradient")

plot(temp3$x,temp3$y, xlab="Gradient", ylab="Abundance", ylim=c(0,6.5))
abline(log(5), 3)
text(x=0.2, y=6, "$\\sigma=1.0$")
#text(x=0.2, y=6, expression(sigma==1.0))
plot(linearT3$IVs[,1], linearT3$IVs[,2], xlab="Gradient",
     ylab="IndVal", xlim=c(0,1))
hist(results3$thetastar[1,], xlim=c(0,1), prob=T, main="", xlab="Gradient")
dev.off()



##
## piecewise linear function (one side slope = 0, sigma=0)
plot(temp1 <- drawY (x=seq(0,1,,100), phi1=0.5, phi2=0., beta0=log(40),
                    beta1=-log(40), sigma=0., Log=F, model=M3))
hockeyT1 <-  T1(temp1$x, ifelse(temp1$y<0,0,temp1$y))
size <- length(temp1$x)
perm(2500, temp1$x, ifelse(temp1$y<0,0,temp1$y))
 results1 <- bootstrap(x=1:size, nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])$split
                      },
                      infile=data.frame(x=temp1$x, y=temp1$y))

plot(temp2 <- drawY (x=seq(0,1,,100), phi1=0.5, phi2=0., beta0=log(40),
                    beta1=log(40), beta2=-log(5), sigma=0., Log=F, model=M3.1))
hockeyT2 <-  T1(temp2$x, ifelse(temp2$y<0,0,temp2$y))
size <- length(temp2$x)
perm(2500, temp2$x, ifelse(temp2$y<0,0,temp2$y))
results2 <- bootstrap(x=1:size, nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])$split
                      },
                      infile=data.frame(x=temp2$x, y=temp2$y))

plot(temp3 <- drawY (x=seq(0,1,,100), phi1=0.5, phi2=0., beta0=log(40),
                     beta1=-log(5), beta2=-log(40), d0=-1, sigma=0.,
                     Log=F, model=M3.2))

hockeyT3 <-  T1(temp3$x, ifelse(temp3$y<0,0,temp3$y))
size <- length(temp3$x)
perm(2500, temp3$x, ifelse(temp3$y<0,0,temp3$y))
results3 <- bootstrap(x=1:size, nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])$split
                      },
                      infile=data.frame(x=temp3$x, y=temp3$y))


#pdf(file=paste(plotDIR, "hockeySims.pdf", sep="/"),height=6.5, width=6.5)
tikz(file=paste(plotDIR, "hockeySimsSA.tex", sep="/"),height=6.5, width=6.5, standAlone=T)
par(mfrow=c(3,3), mar=c(3,3,0.5,0.25), mgp=c(1.5,0.25,0), tck=0.02,las=1)
plot(temp1$x,temp1$y, xlab="Gradient", ylab="Abundance", ylim=c(1,5))
plot(hockeyT1$IVs[,1], hockeyT1$IVs[,2], xlab="Gradient",
     ylab="IndVal", xlim=c(0,1))
hist(results1$thetastar[1,], xlim=c(0,1), prob=T, main="", xlab="Gradient")

plot(temp2$x,temp2$y, xlab="Gradient", ylab="Abundance", ylim=c(1,5))
plot(hockeyT2$IVs[,1], hockeyT2$IVs[,2], xlab="Gradient",
     ylab="IndVal", xlim=c(0,1))
hist(results2$thetastar[1,], xlim=c(0,1), prob=T, main="", xlab="Gradient")

plot(temp3$x,temp3$y, xlab="Gradient", ylab="Abundance", ylim=c(1,5))
plot(hockeyT3$IVs[,1], hockeyT3$IVs[,2], xlab="Gradient",
     ylab="IndVal", xlim=c(0,1))
hist(results3$thetastar[1,], xlim=c(0,1), prob=T, main="", xlab="Gradient")
dev.off()



## reversing x: x = 1-x

stepT <-  T1(1-temp$x, ifelse(temp$y<0,0,temp$y))
size <- length(temp$x)
 results1 <- bootstrap(x=1:size, nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])$split
                      },
                      infile=data.frame(x=1-temp$x, y=temp$y))
par(mfrow=c(1,3))
 plot(1-temp$x,temp$y)
 plot(stepT$IVs[,1], stepT$IVs[,2])
hist(results1$thetastar[1,], xlim=c(0,1),
      xlab=round(quantile(results1$thetastar[1,], prob=c(0.05,0.95)), 2))
mean(results1$thetastar[2,])


## dose-response model (one side slope = 0)
plot(temp <- drawY (x=runif(100), phi1=0.35, phi2=0.65, beta0=log(30),
                    beta1=log(2), sigma=0., model=M4))

stepT <-  T1(temp$x, ifelse(temp$y<0,0,temp$y))
size <- length(temp$x)
 results1 <- bootstrap(x=1:size, nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])$split
                      },
                      infile=data.frame(x=temp$x, y=temp$y))
par(mfrow=c(1,3))
 plot(temp$x,temp$y)
 plot(stepT$IVs[,1], stepT$IVs[,2])
hist(results1$thetastar[1,], xlim=c(0,1),
      xlab=round(quantile(results1$thetastar[1,], prob=c(0.05,0.95)), 2))

## flip x: x <- 1-x
stepT <-  T1(1-temp$x, ifelse(temp$y<0,0,temp$y))
size <- length(temp$x)
 results1 <- bootstrap(x=1:size, nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])$split
                      },
                      infile=data.frame(x=1-temp$x, y=temp$y))
par(mfrow=c(1,3))
 plot(1-temp$x,temp$y)
oo <- order(temp$x)
lines(1-temp$x[oo], temp$z[oo])
plot(stepT$IVs[,1], stepT$IVs[,2], col="gray")
hist(results1$thetastar[1,], xlim=c(0,1),
      xlab=round(quantile(results1$thetastar[1,], prob=c(0.05,0.95)), 2))

## binary

plot(temp <- drawBY (x=runif(100), phi1=0, phi2=0, beta0=logit(0.25),
                    beta1=-7))

stepT <-  T1(temp$x, ifelse(temp$y<0,0,temp$y))
stepT [stepT[,3]==max(stepT[,3]),1]
size <- length(temp$x)
 results1 <- bootstrap(x=1:size, nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])
                        temp [temp[,3]==max(temp[,3]),1]
                      },
                      infile=data.frame(x=temp$x, y=temp$y))
par(mfrow=c(1,3))
 plot(temp$x,temp$y)
 plot(stepT[,1], stepT[,3])
hist(unlist(results1$thetastar), xlim=c(0,1),
      xlab=round(quantile(unlist(results1$thetastar), prob=c(0.025,0.975), na.rm=T), 2))


everg.env <- read.table("glades.env.txt", header=T)
everg.txa <- read.table("glades.taxa.txt", header=T)

perm(2500, everg.env[,1], everg.txa[,6])
stepT <-  T1(everg.env[,1], everg.txa[,7])

results1 <- bootstrap(x=1:dim(everg.env)[1], nboot=2000,
                      theta=function(x, infile) {
                        temp <- T1(infile[x,1], infile[x,2])$split
                      },
                      infile=data.frame(x=everg.env[,1], y=everg.txa[,7]))
par(mfrow=c(1,3))
 plot(everg.env[,1], everg.txa[,7])
 plot(stepT$IVs[,1], stepT$IVs[,3])
hist(unlist(results1$thetastar),
      xlab=round(quantile(unlist(results1$thetastar), prob=c(0.025,0.975), na.rm=T), 2))

## simulation of the permutation "z"'s

# 2 linear models, a step function model and a broken stick model
#
plot(temp11 <- drawY (x=seq(0,1,,100), phi1=0, phi2=0., beta0=3,
                    beta1=0, sigma=0.75, model=M2, Log=F))
linear1 <- perm.sim(250, temp11$x, ifelse(temp11$y<0,0,temp11$y))

plot(temp2 <- drawY (x=seq(0,1,,100), phi1=0, phi2=0., beta0=log(5),
                    beta1=5, sigma=0.5, model=M2, Log=F))
linear2 <- perm.sim(250, temp2$x, ifelse(temp2$y<0,0,temp2$y))

plot(temp3 <- drawY (x=seq(0,1,,100), phi1=0.5, phi2=0., beta0=log(40),
                     beta1=-log(5), beta2=-log(40), d0=-1, sigma=0.25,
                     Log=F, model=M3.2))
hockey.perm <- perm.sim(250, temp3$x, ifelse(temp3$y<0,0,temp3$y))

plot(temp4 <- drawY (x=seq(0,1,,100), phi1=0.5, phi2=0., beta0=log(25),
                    beta1=-log(5), sigma=0.5))
step.perm <- perm.sim(250, temp4$x, ifelse(temp4$y<0,0,temp4$y))

#pdf(file=paste(plotDIR, "permSims.pdf", sep="/"),height=6.5, width=7)
tikz(file=paste(plotDIR, "permSimsSA.tex", sep="/"),height=6.5, width=7, standAlone=T)
par(mfrow=c(4,5), mar=c(3,3,0.2,0.2), mgp=c(1.5,0.25,0.), tck=0.02,las=1)
plot(temp11$x,temp11$y, xlab="Gradient", ylab="Abundance", ylim=c(0,7))
abline(h=3, col="blue")
plot(linear1[,"cp"], linear1[,"IndVal"], xlab="Gradient",
     ylab="$IndVal$", xlim=c(0,1), ylim=c(40,60))
plot(linear1[,"cp"], linear1[,"z"], xlab="Gradient",
     ylab="$z$", xlim=c(0,1), ylim=c(-3,3))
plot(linear1[,"cp"], linear1[,"mu"], xlab="Gradient",
     ylab="$\\mu$", xlim=c(0,1), ylim=c(45,65))
#     ylab=expression(mu), xlim=c(0,1), ylim=c(45,65))
plot(linear1[,"cp"], linear1[,"sd"], xlab="Gradient",
     ylab="$\\sigma$", xlim=c(0,1), ylim=c(0,2))
#     ylab=expression(sigma), xlim=c(0,1))

plot(temp2$x,temp2$y, xlab="Gradient", ylab="Abundance")#, ylim=c(0,7))
abline(log(5), 5, col="blue")
plot(linear2[,"cp"], linear2[,"IndVal"], xlab="Gradient",
     ylab="$IndVal$", xlim=c(0,1))#, ylim=c(40,60))
#     ylab="IndVal", xlim=c(0,1))#, ylim=c(40,60))
plot(linear2[,"cp"], linear2[,"z"], xlab="Gradient",
     ylab="$z$", xlim=c(0,1))#, ylim=c(-3,3))
#     ylab="z", xlim=c(0,1))#, ylim=c(-3,3))
plot(linear2[,"cp"], linear2[,"mu"], xlab="Gradient",
     ylab="$\\mu$", xlim=c(0,1), ylim=c(45,65))
#     ylab=expression(mu), xlim=c(0,1), ylim=c(45,65))
plot(linear2[,"cp"], linear2[,"sd"], xlab="Gradient",
     ylab="$\\sigma$", xlim=c(0,1), ylim=c(0.5,3))
#     ylab=expression(sigma), xlim=c(0,1))

plot(temp3$x,temp3$y, xlab="Gradient", ylab="Abundance")#, ylim=c(0,7))
segments(x0=c(0,0.5), x1=c(0.5, 1),
         y0=c(log(40)-log(5)*(0-0.5),log(40)-1),
         y1=c(log(40),log(40)-1-log(40)*(1-0.5)), col="blue")
plot(hockey.perm[,"cp"], hockey.perm[,"IndVal"], xlab="Gradient",
     ylab="$IndVal$", xlim=c(0,1))#, ylim=c(40,60))
#     ylab="IndVal", xlim=c(0,1))#, ylim=c(40,60))
plot(hockey.perm[,"cp"], hockey.perm[,"z"], xlab="Gradient",
     ylab="$z$", xlim=c(0,1))#, ylim=c(-3,3))
#     ylab="z", xlim=c(0,1))#, ylim=c(-3,3))
plot(hockey.perm[,"cp"], hockey.perm[,"mu"], xlab="Gradient",
     ylab="$\\mu$", xlim=c(0,1), ylim=c(45,65))
#     ylab=expression(mu), xlim=c(0,1), ylim=c(45,65))
plot(hockey.perm[,"cp"], hockey.perm[,"sd"], xlab="Gradient",
     ylab="$\\sigma$", xlim=c(0,1), ylim=c(1,4))
#     ylab=expression(sigma), xlim=c(0,1))

plot(temp4$x,temp4$y, xlab="Gradient", ylab="Abundance")
segments(x0=c(0, 0.5), x1=c(0.5,1),
         y0=c(25, 5), y1=c(25, 5), col="blue")
plot(step.perm[,"cp"], step.perm[,"IndVal"], xlab="Gradient",
     ylab="$IndVal$", xlim=c(0,1))#, ylim=c(40,60))
#     ylab="IndVal", xlim=c(0,1))#, ylim=c(40,60))
plot(step.perm[,"cp"], step.perm[,"z"], xlab="Gradient",
     ylab="$z$", xlim=c(0,1))#, ylim=c(-3,3))
#     ylab="z", xlim=c(0,1))#, ylim=c(-3,3))
plot(step.perm[,"cp"], step.perm[,"mu"], xlab="Gradient",
     ylab="$\\mu$", xlim=c(0,1), ylim=c(45,65))
#     ylab=expression(mu), xlim=c(0,1), ylim=c(45,65))
plot(step.perm[,"cp"], step.perm[,"sd"], xlab="Gradient",
     ylab="$\\sigma$", xlim=c(0,1), ylim=c(2,7))
#     ylab=expression(sigma), xlim=c(0,1))
dev.off()


## skewed x -- synchonicity

x <- rlnorm(100)
x <- x/max(x)

tikz(file=paste(plotDIR, "histXSA.tex", sep="/"),height=3, width=3, standAlone=T)
par(mar=c(3, 3, 1, 1), mgp=c(1.5,0.25,0), tck=0.01,las=1)
hist(x, xlab="gradient", main="")
dev.off()

hist(x, log="x")

plot(temp1 <- drawY (x=x, phi1=0, phi2=0., beta0=10,
                    beta1=-log(5), sigma=0.0, model=M2, Log=F))
linear1 <- perm.sim(250, temp1$x, ifelse(temp1$y<0,0,temp1$y))

plot(temp2 <- drawY (x=x, phi1=0, phi2=0., beta0=10,
                    beta1=-log(5), sigma=0.06, model=M2, Log=F))
linear2 <- perm.sim(250, temp2$x, ifelse(temp2$y<0,0,temp2$y))

plot(temp3 <- drawY (x=x, phi1=0, phi2=0., beta0=10,
                    beta1=-log(5), sigma=0.5, model=M2, Log=F))
linear3 <- perm.sim(250, temp3$x, ifelse(temp3$y<0,0,temp3$y))

plot(temp4 <- drawY (x=x, phi1=0, phi2=0., beta0=10,
                    beta1=-log(5), sigma=0.75, model=M2, Log=F))
linear4 <- perm.sim(250, temp4$x, ifelse(temp4$y<0,0,temp4$y))


tikz(file=paste(plotDIR, "permSims2SA.tex", sep="/"),height=6.5, width=7, standAlone=FALSE)
par(mfrow=c(4,5), mar=c(0,3,3,0.2), mgp=c(1.55,0.125,0.), tck=0.02,las=1)
plot(temp1$x,temp1$y, xlab=" ", ylab="Abundance", ylim=c(7,12), xlim=c(0,1), axes=F)
axis(2)
axis(3)
box()
text(x=0.8, y=11, "$cv=0$")
#text(x=0.5, y=11.9, "Abundance")
plot(linear1[,"cp"], linear1[,"IndVal"], xlab="",
     ylab="$IndVal$", xlim=c(0,1), ylim=c(50,55), axes=F)
axis(2)
axis(3)
box()
#text(x=0.5, y=54, "$IndVal$")
plot(linear1[,"cp"], linear1[,"z"], xlab=" ",
     ylab="$z$", xlim=c(0,1), axes=F)
axis(2)
axis(3)
box()
#text(x=0.5, y=range(linear1[,"z"])[1]+0.05*diff(range(linear1[,"z"])), "$z$")
plot(linear1[,"cp"], linear1[,"mu"], xlab=" ",
     ylab="$\\mu$", xlim=c(0,1), ylim=c(50,51), axes=F)
axis(2)
axis(3)
box()
#text(x=0.5, y=50.8, "$\\mu$")
plot(linear1[,"cp"], linear1[,"sd"], xlab=" ",
     ylab="$\\sigma$", xlim=c(0,1), axes=F)
axis(2)
axis(3)
box()
#text(x=0.25, y=range(linear1[,"sd"])[2]-0.025*diff(range(linear1[,"sd"])), "$\\sigma$")

par(mar=c(1.5,3,1.5,0.2))
plot(temp2$x,temp2$y, xlab=" ", ylab="Abundance", ylim=c(7,12), xlim=c(0,1))
text(x=0.8, y=11, "$cv=25\\%$")
#text(x=0.5, y=11.9, "Abundance")
plot(linear2[,"cp"], linear2[,"IndVal"], xlab=" ",
     ylab="$IndVal$", xlim=c(0,1), ylim=c(50,55))
#text(x=0.5, y=54, "$IndVal$")
plot(linear2[,"cp"], linear2[,"z"], xlab=" ",
     ylab="$z$", xlim=c(0,1))
#text(x=0.5, y=range(linear2[,"z"])[1]+0.05*diff(range(linear2[,"z"])), "$z$")
plot(linear2[,"cp"], linear2[,"mu"], xlab=" ",
     ylab="$\\mu$", xlim=c(0,1), ylim=c(50,51))
#text(x=0.5, y=50.8, "$\\mu$")
plot(linear2[,"cp"], linear2[,"sd"], xlab=" ",
     ylab="$\\sigma$", xlim=c(0,1))
#text(x=0.25, y=range(linear2[,"sd"])[2]-0.025*diff(range(linear2[,"sd"])), "$\\sigma$")

plot(temp3$x,temp3$y, xlab=" ", ylab="Abundance", ylim=c(7,12), xlim=c(0,1), axes=F)
axis(2)
box()
text(x=0.8, y=11, "$cv=80\\%$")
#text(x=0.5, y=11.9, "Abundance")
plot(linear3[,"cp"], linear3[,"IndVal"], xlab=" ",
     ylab="$IndVal$", xlim=c(0,1), ylim=c(50,55), axes=F)
axis(2)
box()
#text(x=0.5, y=54, "$IndVal$")
plot(linear3[,"cp"], linear3[,"z"], xlab=" ",
     ylab="$z$", xlim=c(0,1), axes=F)
axis(2)
box()
#text(x=0.5, y=range(linear3[,"z"])[1]+0.05*diff(range(linear3[,"z"])), "$z$")
plot(linear3[,"cp"], linear3[,"mu"], xlab=" ",
     ylab="$\\mu$", xlim=c(0,1), ylim=c(50,51), axes=F)
axis(2)
box()
#text(x=0.5, y=50.8, "$\\mu$")
plot(linear3[,"cp"], linear3[,"sd"], xlab=" ",
     ylab="$\\sigma$", xlim=c(0,1), axes=F)
axis(2)
box()
#text(x=0.25, y=range(linear3[,"sd"])[2]-0.025*diff(range(linear3[,"sd"])), "$\\sigma$")

par(mar=c(3,3,0,0.2))
plot(temp4$x,temp4$y, xlab="gradient", ylab="Abundance", ylim=c(7,12), xlim=c(0,1))
text(x=0.8, y=11, "$cv=100\\%$")
#text(x=0.5, y=11.9, "Abundance")
plot(linear4[,"cp"], linear4[,"IndVal"], xlab="gradient",
     ylab="$IndVal$", xlim=c(0,1), ylim=c(50,55))
#text(x=0.5, y=54, "$IndVal$")
plot(linear4[,"cp"], linear4[,"z"], xlab="gradient",
     ylab="$z$", xlim=c(0,1))
#text(x=0.5, y=range(linear4[,"z"])[1]+0.05*diff(range(linear4[,"z"])), "$z$")
plot(linear4[,"cp"], linear4[,"mu"], xlab="gradient",
     ylab="$\\mu$", xlim=c(0,1), ylim=c(50,51))
#text(x=0.5, y=50.8, "$\\mu$")
plot(linear4[,"cp"], linear4[,"sd"], xlab="gradient",
     ylab="$\\sigma$", xlim=c(0,1))
#text(x=0.25, y=range(linear4[,"sd"])[2]-0.025*diff(range(linear4[,"sd"])), "$\\sigma$")

dev.off()

