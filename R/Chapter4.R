source("FrontMatter.R")

################
## Chapter 4  ##
################
plotDIRch4 <- paste(plotDIR, "chapter4","figures", sep="/")

    y <- rnorm(20, 0.2, 1)
    n <- length (y)
    y.bar <- mean(y)
    se <- sd(y)/sqrt(n)
    int.50 <- y.bar + qt(c(0.25, 0.75), df=n-1)*se
    int.95 <- y.bar + qt(c(.025, .975), df=n-1)*se
    print( c(y.bar, int.95))

 n.sims <- 1000
    n.size <- 30
    inside <- 0
    for (i in 1:n.sims){  ## looping through n.sims iterations
        y <- rnorm(n.size, mean=2.05, sd=0.34)
            ## random samples from N(2.05, 0.34)
        se <- sd(y)/sqrt(n.size)
        int.95 <- mean(y) + qt(c(.025, .975), n.size-1)*se
        inside <- inside + sum(int.95[1]<2.05 & int.95[2]>2.05)
        }
    inside/n.sims  ## fraction of times true mean inside int.95

 n.sims <- 1000
    n.size <- 30
    inside <- 0
    for (i in 1:n.sims){  ## looping through n.sims iterations
        y <- runif(n.size, 1.05, 3.05)
            ## random samples from unif(1.05, 3.05)
        se <- sd(y)/sqrt(n.size)
        int.95 <- mean(y) + qt(c(.025, .975), n.size-1)*se
        inside <- inside + sum(int.95[1]<2.05 & int.95[2]>2.05)
        }
    inside/n.sims  ## fraction of times true mean inside int.95

## central limit theorem simulation
two.prob()
central.sim(mux=1, vx=1, n=c(5, 20, 100))
#postscript(file=paste(plotDIR, "cltSims.eps", sep="/"),
#           width=4.75, height=3, horizontal=F)
central.sim(mux=1, vx=1, n=c(5, 20, 100))
#dev.off()

## Everglades example

table(TP.reference$Year, TP.reference$Month)

histogram(~log(TP) | SITE, data=TP.reference)

#postscript(file=paste(plotDIR,"EvergQQs.eps", sep="/"),
#           width=4, height=3.5, horizontal=F)
tikz(file=paste(plotDIRch4, "EvergQQs.tex", sep="/"),
     width=4, height=3.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.set(list(layout.widths=list(left.padding=0, right.padding=0,
                       ylab.axis.padding=0, axis.right=0, key.ylab.padding=0)))
print(
qqmath(~log(TP) | factor(Year),
  panel=function(x, ...){
    panel.grid()
    panel.qqmathline(x, ...)
    panel.qqmath(x, ...)
    }, data=TP.reference,
    subset=SITE!="E5"&SITE!="F5",
    par.strip.text=list(cex=0.5),
    xlab="Unit Normal Quantile",
    ylab="Log TP Concentration (ppb)")
)
dev.off()

qqmath(~log(TP) | SITE,
  panel=function(x, ...){
    panel.qqmathline(x, ...)
    panel.qqmath(x, ...)
    }, data=TP.reference,
    subset=SITE!="E5"&SITE!="F5")


qqmath(~log(TP) | Month,
  panel=function(x, ...){
    panel.qqmathline(x, ...)
    panel.qqmath(x, ...)
    }, data=TP.reference)

## precipitation
EvergPrecip <- read.table(paste(dataDIR, "EvergPrecip.txt", sep="/"),
                          header=F)
EvergPrecip2 <- read.table(paste(dataDIR, "EvergladesP.txt", sep="/"),
                           header=T)

EvergPrecip2$ANN2 <- apply(EvergPrecip2[,2:13], 1, sum)
EvergPrecip2$Msd <- apply(EvergPrecip2[,2:13], 1, sd)

Y.lim <- c(0.99, 1.1)*range(c(EvergPrecip2$ANN2+2*EvergPrecip2$Msd),
                            c(EvergPrecip2$ANN2-2*EvergPrecip2$Msd), na.rm=T)
#postscript(file=paste(plotDIR, "EvergPrecip.eps", sep="/"),
#           width=3.75, height=2.5, horizontal=F)
tikz(file=paste(plotDIRch4, "EvergPrecip.tex", sep="/"),
           width=3.75, height=2.5, standAlone=F)
par(mar=c(2.5,2.5,0.5,0.25), mgp=c(1.5, 0.5, 0))
plot(EvergPrecip2$ANN2/12 ~ EvergPrecip2$YEAR, type="n",
     xlab="Year", ylab="Annual Precipitation (in)", ylim= Y.lim, xlim=c(1990,2003))
segments(x0= EvergPrecip2$YEAR, x1=EvergPrecip2$YEAR,
         y0=EvergPrecip2$ANN2-2*EvergPrecip2$Msd,
         y1=EvergPrecip2$ANN2+2*EvergPrecip2$Msd)
segments(x0= EvergPrecip2$YEAR, x1=EvergPrecip2$YEAR,
         y0=EvergPrecip2$ANN2-EvergPrecip2$Msd,
         y1=EvergPrecip2$ANN2+EvergPrecip2$Msd, lwd=3)
points(EvergPrecip2$YEAR, EvergPrecip2$ANN2)
segments(x0=c(1993.5, 1999.5), x1=c(1993.5, 1999.5),
         y0=rep(Y.lim[1], 2), y1=rep(Y.lim[2], 2), lwd=2, col="grey")
segments(x0=rep(1993.5, 2), x1=rep(1999.5, 2),
         y0=c(Y.lim[1], Y.lim[2]), y1=c(Y.lim[1], Y.lim[2]), lwd=2, col="grey")
abline(h=mean(EvergPrecip2$ANN2, na.rm=T), col=gray(0.5) )
dev.off()

by(TP.reference$Month, TP.reference$Year, table)

table(TP.reference$Year)

TP.impacted <- wca2tp[wca2tp$Type=="I",]
subI <- (TP.impacted$SITE=="E4"|TP.impacted$SITE=="F4")&TP.impacted$Year==1994
subS <- TP.reference$Year==1994&TP.reference$SITE!="E5"&TP.reference$SITE!="F5"
two.sample <- data.frame(TP = c(TP.impacted$RESULT[subI],
                           TP.reference$RESULT[subS]),
                         Type=c(rep("I", sum(subI)), rep("R", sum(subS))))
qq(Type~TP, data=two.sample)
qqmath(~log(TP)|SITE, data=TP.impacted)

x<-log(TP.impacted$TP[subI])
y<-log(TP.reference$TP[subS])
t.test(x, y, alternative="t", var.equal=T)

y.bar <- mean(log(TP.reference$TP[subS]))
stdv<-sd(log(TP.reference$TP[subS]))
n <- sum(subS)
se <- sd(log(TP.reference$TP[subS]))/sqrt(n)
int.50 <- y.bar + qt(c(0.25, 0.75), n-1)*se
int.95 <- y.bar + qt(c(.025, .975), n-1)*se

n.sims <- 1000
n.size <- 30
inside <- 0
for (i in 1:n.sims){  ## looping through 5000 iterations
  y <- rnorm(n.size, mean=2.05, sd=0.34)
  ## random samples from N(2.05, 0.34)
  se <- sd(y)/sqrt(n.size)
  int.95 <- mean(y) + qt(c(.025, .975), n.size-1)*se
  inside <- inside + sum(int.95[1]<2.05 & int.95[2]>2.05)
}
inside/n.sims  ## fraction of times true mean inside int.95


n.sims <- 10000
n.size <- 30
inside <- 0
for (i in 1:n.sims){  ## looping through 5000 iterations
  y <- runif(n.size, 1, 3)
  ## random samples from uniform (1,3) with a mean of 2
  se <- sd(y)/sqrt(n.size)
  int.95 <- mean(y) + qt(c(.025, .975), n.size-1)*se
  inside <- inside + sum(int.95[1]<2 & int.95[2]>2)
}
inside/n.sims  ## fraction of times true mean inside int.95

##
## Simulation of central limit theorem
##


## uncertainty of sigma
n.sims <- 1000
subS <- TP.reference$Year==1994&TP.reference$SITE!="E5"&TP.reference$SITE!="F5"
y.bar <- mean(log(TP.reference$TP[subS]))
n <- sum(subS)
se <- sd(log(TP.reference$TP[subS]))
X <- rchisq (n.sims, df=n-1)
sigma.chi2 <- se * sqrt((n-1)/X)
##postscript(file="sigmaChi2.eps", width=3, height=2, horizontal=F)
tikz(file=paste(plotDIRch4, "sigmaChi2.tex", sep="/"),
     width=3, height=2, standAlone=F)
par(mar=c(2.5,.75,2.5,0.25), mgp=c(1.5, 0.5, 0), las=1, tck=0.01)
hist(sigma.chi2, axes=F, xlab="$\\hat{\\sigma}$", main="")
axis(1, cex=0.5)
abline(v=se, lwd=3, col="grey")
dev.off()

sample.mean <- rnorm(n.sims, y.bar, sigma.chi2/sqrt(n))
q.75 <- qnorm(0.75, sample.mean, sigma.chi2)

##postscript(file=paste(plotDIR, "q75.eps", sep="/"),
##           width=3, height=2, horizontal=F)
tikz(file=paste(plotDIRch4, "q75.tex", sep="/"),
           width=3, height=2, standAlone=F)
par(mar=c(2.5,.75,2.5,0.25), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
hist(exp(q.75), axes=F, xlab="0.75 Quantile Distribution", main="")
axis(1)
dev.off()

quantile(exp(q.75), prob=c(0.025,0.975))
## bootstraping
y <- log(TP.reference$TP[subS])
x <- c(94, 38, 23, 197, 99, 16, 141)
boot.theta <- numeric()
quantile(exp(y), prob=0.75)
B <- 10000
for (i in 1:B){
  boot.sample <- sample(y, size=length(y), T)
  boot.theta[i] <- quantile(boot.sample, prob=0.75)
}
sd(boot.theta)

packages(bootstrap)

y <- log(TP.reference$TP[subS])
results <- bootstrap(y, 2000, quantile, prob=0.75)
CI.t <- c(mean(results $ thetastar) - qt(0.975, 29),
          mean(results $ thetastar) + qt(0.975, 29))
CI.t
##[1] 0.1997 4.2902

CI.percent <- quantile(results$thetastar, prob=c(0.025, 0.975))


bca.results <- bcanon(y,2000,theta=quantile, prob=0.75, alpha=c(0.025, 0.975))
bca.results$confpoints

### t and normal ###
##postscript(file=paste(plotDIR, "normalT.eps", sep="/"),
##           height=2, width=2.5, horizontal=F)
tikz(file=paste(plotDIRch4, "normalT.tex", sep="/"),
           height=2, width=2.5, standAlone=F)
par(mar=c(2.5,.75,2.5,0.25), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
plot(seq(-3,3,,100), dnorm(seq(-3,3,,100)), type="l",
     col="gray", axes=F, xlab="",ylab="")
curve(dt(x, 3),add=T)
axis(1)
box()
dev.off()

### rejection region ###
##postscript(file=paste(plotDIR, "rejectR1.eps", sep="/"),
##           height=3.75, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch4, "rejectR1.tex", sep="/"),
           height=3.75, width=4.5, standAlone=F)
par(mfrow=c(2,1),mar=c(2,0,0,0), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
## the null
plot(seq(-3,3,,100), dnorm(seq(-3,3,,100)), type="l",
     xlab="", ylab="", axes=F, xlim=c(-3,6))
polygon(c(1.75, seq(1.75, 3,,50), 3), c(0, dnorm(seq(1.75, 3,,50)), 0),
        col=gray(.8))
polygon(c(2.25, seq(2.25, 3,,50), 3), c(0, dnorm(seq(2.25, 3,,50)), 0),
        density=20, angle=45)
text(1.9, dnorm(1.9)/2, "$\\alpha$")
text(2.25, dnorm(0)/1.8, "$t_{obs}$", adj=c(0,1))
segments(x0=2.25, x1=2.25, y0=0, y1=dnorm(0)/2.1)
text(2.75, dnorm(2.5)*4.2, "$p$-value", adj=c(0,0))
text(6, dnorm(0)*(2/3), "$H_0$", adj=1)
arrows(x0=2.5, y0=dnorm(2.5)/2, x1=2.75, y1=dnorm(2.5)*4, length=0.1,
       angle=25, code=1)
abline(v=1.75)
axis(1, at=seq(-2,6,2)[-3], label=rep("", 4))
axis(1, at=c(1.75), label="$t_{cutoff}$")
axis(1, at=c(0), label="$\\mu_0$")
## the alternative
plot(seq(0,6,,100), dnorm(seq(0,6,,100), mean=3), type="l",
     xlab="", ylab="", axes=F, xlim=c(-3,6))
#polygon(c(-1, seq(-1, 1.3,,50), 1.3), c(0, dnorm(seq(-1, 1.3,,50)), 0),
#col=gray(.65))
axis(1, at=c(3), label="$\\mu_a$")
axis(1, at=seq(-2,6,2)[-3], label=rep("", 4))
abline(v=1.75)
#polygon(c(1.75, seq(1.75, 6,,50), 6), c(0, dnorm(seq(1.75, 6,,50),3 ), 0),
#density=20, angle=45)
polygon(c(0, seq(0, 1.75,,50), 1.75), c(0, dnorm(seq(0, 1.75,,50),3 ), 0),
        density=15, angle=-45)
text(3, dnorm(0)/2, "$1-\\beta$")
text(3, dnorm(0)/4, "(power)")
text(1, dnorm(1.75, 3)/2, "$\\beta$")
text(-2, dnorm(0)*(2/3), "$H_a$", adj=0)
dev.off()

### two-sided p-value and alpha
##postscript(file=paste(plotDIR, "rejectR2.eps", sep="/"),
##           height=2., width=4.5, horizontal=F)
tikz(file=paste(plotDIRch4, "rejectR2.tex", sep="/"),
           height=2., width=4.5, standAlone=F)
par(mar=c(2,0,0,0), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
## the null
plot(seq(-3,3,,100), dnorm(seq(-3,3,,100)), type="l",
     xlab="", ylab="", axes=F)
polygon(c(1.75, seq(1.75, 3,,50), 3), c(0, dnorm(seq(1.75, 3,,50)), 0),
        col=gray(.8))
polygon(c(-3, seq(-3, -1.75,,50), -1.75), c(0, dnorm(seq(-3, -1.75,,50)), 0),
        col=gray(.8))
polygon(c(2.25, seq(2.25, 3,,50), 3), c(0, dnorm(seq(2.25, 3,,50)), 0),
        density=20, angle=45)
polygon(c(-3, seq(-3,-2.25,,50), -2.25), c(0, dnorm(seq(-3, -2.25,,50)), 0),
        density=20, angle=-45)

text(3, dnorm(2.5)*4.2, "$p$-value/2", adj=c(1,0), cex=0.75)
text(-3, dnorm(2.5)*4.2, "$p$-value/2", adj=c(0,0), cex=0.75)

arrows(x0=2.5, y0=dnorm(2.5)/2, x1=2.75, y1=dnorm(2.5)*4, length=0.1,
       angle=25, code=1)
arrows(x0=-2.5, y0=dnorm(2.5)/2, x1=-2.75, y1=dnorm(2.5)*4, length=0.1,
       angle=25, code=1)

axis(1, at=c(1.75), label="$t_{cutoff}$")
axis(1, at=c(2.25), label="$t_{obs}$")
axis(1, at=c(-1.75), label="$-t_{cutoff}$")
axis(1, at=c(-2.25), label="$-t_{obs}$")
axis(1, at=c(0), label="$\\mu_0$")
text( 1.8, dnorm(1.9)/2.2, "$\\alpha/2$", adj=c(0, 0.5), cex=0.75)
text(-1.8, dnorm(1.9)/2.2, "$\\alpha/2$", adj=c(1, 0.5), cex=0.75)
axis(1, at=c(-3,3), label=c("", ""))
dev.off()


##postscript(file="power2.ps", height=4, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch4, "power2.tex", sep="/"),
           height=4, width=4.5, standAlone=F)
par(mfrow=c(2,2), mar=c(1.5, .25, 0.5, 0.), mgp=c(1.5, 0.5, 0))
z <- qnorm(1-0.05/2)
plot(seq(7, 17,,100), dnorm(seq(7,17,,100), 10, 2/sqrt(10)), type="l",
     xlab="", ylab="", ylim=c(0, 1.1), axes=F)
lines(seq(7, 17,,100), dnorm(seq(7,17,,100), 12, 2/sqrt(10)), type="l",
      col=gray(0.5))
abline(v=10+z*2/sqrt(10), lty=5, col=gray(0.5))
text(10, dnorm(10,10, 2/sqrt(10)), "$H_0$", adj=c(0.5,0))
text(12, dnorm(12,12, 2/sqrt(10)), "$H_a$", adj=c(0.5,0))
text(13, 0.6, "$n=10$", adj=c(0,0.5))
text(13, 0.5, "$\\sigma=2$", adj=c(0,0.5))
text(13, 0.4, "$\\alpha=0.05$", adj=c(0,0.5))
text(13, 0.75, paste("Power = ",
                     round(1-pnorm(10+z*2/sqrt(10), 12, 2/sqrt(10)), digit=3),
                     sep=""), adj=c(0,0.5))
axis(1)

plot(seq(7, 17,,100), dnorm(seq(7,17,,100), 10, 3/sqrt(10)), type="l",
     xlab="", ylab="", ylim=c(0, 1.1), axes=F)
lines(seq(7, 17,,100), dnorm(seq(7,17,,100), 12, 3/sqrt(10)), type="l",
      col=gray(0.5))
abline(v=10+z*3/sqrt(10), lty=5, col=gray(0.5))
#abline(v=5:10, col=1:6)
text(10, dnorm(10, 10, 3/sqrt(10)), "$H_0$", adj=c(0.5,0))
text(12, dnorm(12, 12, 3/sqrt(10)), "$H_a$", adj=c(0.5,0))
text(12.75, 0.6, "$n=10$", adj=c(0,0.5))
text(12.75, 0.5, "$\\sigma=3$", adj=c(0,0.5))
text(12.75, 0.4, "$\\alpha=0.05$", adj=c(0,0.5))
text(12.75, 0.75,
     paste("Power = ", round(1-pnorm(10+z*3/sqrt(10), 12, 4/sqrt(10)), digit=3),
           sep=""), adj=c(0,0.5))
axis(1)

plot(seq(7, 17,,100), dnorm(seq(7,17,,100), 10, 2/sqrt(20)), type="l",
     xlab="", ylab="", ylim=c(0, 1.1), axes=F)
lines(seq(7, 17,,100), dnorm(seq(7,17,,100), 12, 2/sqrt(20)), type="l",
      col=gray(0.5))
abline(v=10+z*2/sqrt(20), lty=5, col=gray(0.5))
text(10, dnorm(10, 10, 2/sqrt(20)), expression(H[0]), adj=c(0.5,0))
text(12, dnorm(12, 12, 2/sqrt(20)), expression(H[a]), adj=c(0.5,0))
text(12.75, 0.6, "$n=20$", adj=c(0,0.5))
text(12.75, 0.5, "$\\sigma=2$", adj=c(0,0.5))
text(12.75, 0.4, "$\\alpha=0.05$", adj=c(0,0.5))
text(12.75, 0.75, paste("Power = ",
                        round(1-pnorm(10+z*2/sqrt(20), 12, 2/sqrt(20)), digit=3),
                        sep=""), adj=c(0,0.5))
axis(1)

z <- qnorm(1-0.1/2)
plot(seq(7, 17,,100), dnorm(seq(7,17,,100), 10, 3/sqrt(10)), type="l",
     xlab="", ylab="", ylim=c(0, 1), axes=F)
lines(seq(7, 17,,100), dnorm(seq(7,17,,100), 12, 3/sqrt(10)), type="l",
      col=gray(0.5))
abline(v=10+z*3/sqrt(10), lty=5, col=gray(0.5))
text(10, dnorm(10, 10, 3/sqrt(10)), expression(H[0]), adj=c(0.5, 0))
text(12, dnorm(12, 12, 3/sqrt(10)), expression(H[a]), adj=c(0.5, 0))
text(12.75, 0.6, "$n=10$", adj=c(0,0.5))
text(12.75, 0.5, "$\\sigma=3$", adj=c(0,0.5))
text(12.75, 0.4, "$\\alpha=0.1$", adj=c(0,0.5))
text(12.75, 0.75, paste("Power = ",
                        round(1-pnorm(10+z*3/sqrt(10), 12, 4/sqrt(10)), digit=3),
                        sep=""), adj=c(0,0.5))
axis(1)
dev.off()

## comment on nonparametric
packages(exactRankTests)
hypo.sim <- function(n.sims, rdistF, theta, ...){
    reject.t1<-0; reject.t2<-0; reject.w1<-0; reject.w2<-0
    for (i in 1:n.sims){
        u <- rdistF(20, ...)
        y <- u
        for (j in 2:20)
            y[j] <- u[j] - theta*u[j-1]
        samp1 <- data.frame(x=y, g=sample(1:2, 20, TRUE))
### randomized sample
        samp2 <- data.frame(x=y, g=rep(c(1,2), each=10))
### correlated sample
        reject.t1 <- reject.t1 +
            (t.test(x~g, data=samp1, var.equal=T)$p.value<0.05)
        reject.t2 <- reject.t2 +
            (t.test(x~g, data=samp2, var.equal=T)$p.value<0.05)
        reject.w1 <- reject.w1 +
            (wilcox.exact(x~g, data=samp1)$p.value<0.05)
        reject.w2 <- reject.w2 +
            (wilcox.exact(x~g, data=samp2)$p.value<0.05)
    }
    return(rbind(c(reject.t2,reject.t1),
                 c(reject.w2,reject.w1))/n.sims)
}

hypo.sim(n.sims=1000,rdistF=rnorm,theta=-0.4,mean=2,sd=4)
hypo.sim(n.sims=1000,rdistF=rpois,theta=-0.4,lambda=3)
hypo.sim(n.sims=1000,rdistF=runif,theta=-0.4,max=3,min=-3)

hypo.sim(n.sims=1000,rdistF=rnorm,theta=0.4,mean=2,sd=4)
hypo.sim(n.sims=1000,rdistF=rpois,theta=0.4,lambda=3)
hypo.sim(n.sims=1000,rdistF=runif,theta=0.4,max=3,min=-3)


### Examples ###
## Everglades example

#postscript(file=paste(plotDIR,"EvergQQs.eps", sep="/"),
#           width=4, height=3.5, horizontal=F)
tikz(file=paste(plotDIRch4, "EvergQQs.tex", sep="/"),
     width=4, height=3.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.set(list(layout.widths=list(left.padding=0, right.padding=0,
                       ylab.axis.padding=0, axis.right=0, key.ylab.padding=0)))
print(
qqmath(~log(RESULT) | factor(Year),
  panel=function(x, ...){
    panel.grid()
    panel.qqmathline(x, ...)
    panel.qqmath(x, ...)
    }, data=TP.reference,
    subset=SITE!="E5"&SITE!="F5",
    par.strip.text=list(cex=0.5),
    xlab="Unit Normal Quantile",
    ylab="Log TP Concentration (ppb)")
)
dev.off()

qqmath(~log(RESULT) | SITE,
  panel=function(x, ...){
    panel.qqmathline(x, ...)
    panel.qqmath(x, ...)
    }, data=TP.reference,
    subset=SITE!="E5"&SITE!="F5")


qqmath(~log(RESULT) | Month,
  panel=function(x, ...){
    panel.qqmathline(x, ...)
    panel.qqmath(x, ...)
    }, data=TP.reference)

## precipitation
EvergPrecip <- read.table(paste(dataDIR, "EvergPrecip.txt", sep="/"),
                          header=F)
EvergPrecip2 <- read.table(paste(dataDIR, "EvergladesP.txt", sep="/"),
                           header=T)

EvergPrecip2$ANN2 <- apply(EvergPrecip2[,2:13], 1, sum)
EvergPrecip2$Msd <- apply(EvergPrecip2[,2:13], 1, sd)

Y.lim <- c(0.99, 1.1)*2.54*range(c(EvergPrecip2$ANN2+2*EvergPrecip2$Msd),
                            c(EvergPrecip2$ANN2-2*EvergPrecip2$Msd), na.rm=T)
#postscript(file=paste(plotDIR, "EvergPrecip.eps", sep="/"),
#           width=3.75, height=2.5, horizontal=F)
tikz(file=paste(plotDIRch4, "EvergPrecip.tex", sep="/"),
           width=3.75, height=2.5, standAlone=F)
par(mar=c(2.5,2.5,0.5,0.25), mgp=c(1.5, 0.5, 0))
plot(2.54*EvergPrecip2$ANN2/12 ~ EvergPrecip2$YEAR, type="n",
     xlab="Year", ylab="Annual Precipitation (cm)", ylim= Y.lim, xlim=c(1990,2003))
segments(x0= EvergPrecip2$YEAR, x1=EvergPrecip2$YEAR,
         y0=2.54*EvergPrecip2$ANN2-2*2.54*EvergPrecip2$Msd,
         y1=2.54*EvergPrecip2$ANN2+2*2.54*EvergPrecip2$Msd)
segments(x0= EvergPrecip2$YEAR, x1=EvergPrecip2$YEAR,
         y0=2.54*EvergPrecip2$ANN2-2.54*EvergPrecip2$Msd,
         y1=2.54*EvergPrecip2$ANN2+2.54*EvergPrecip2$Msd, lwd=3)
points(EvergPrecip2$YEAR, 2.54*EvergPrecip2$ANN2)
segments(x0=c(1993.5, 1999.5), x1=c(1993.5, 1999.5),
         y0=rep(Y.lim[1], 2), y1=rep(Y.lim[2], 2), lwd=2, col="grey")
segments(x0=rep(1993.5, 2), x1=rep(1999.5, 2),
         y0=c(Y.lim[1], Y.lim[2]), y1=c(Y.lim[1], Y.lim[2]), lwd=2, col="grey")
abline(h=mean(2.54*EvergPrecip2$ANN2, na.rm=T), col=gray(0.5) )
dev.off()


#### 303(d) listing power ####
sample.size <- 10:50


reject <- qbinom(1-0.05, size=sample.size, prob=0.1) + 1
decision.table <- data.frame(n=sample.size, reject=reject,
                             typeI = 1-pbinom(reject-1, sample.size, 0.1),
                             power=1-pbinom(reject-1,
                               size=sample.size, prob=0.25))
##postscript(file="cwa303d.eps", width=4.5, height=2, horizontal=F)
tikz(file=paste(plotDIRch4, "cwa303d.tex", sep="/"),
           width=4.5, height=2, standAlone=F)
par(mfrow=c(1,2),mar=c(3, 3, 0.25, 0.25), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
plot(typeI ~ n,data=decision.table, type="l", ylim=c(0,0.1),
     xlab="Sample Size", ylab="$\\alpha$")
plot(power~n, data=decision.table, type="l",
     xlab="Sample Size", ylab="$1-\\beta$")
dev.off()

Everg.aov <- aov(log(TP) ~ factor(Year), data=TP.reference)
##postscript(file=paste(plotDIR, "aovResids1.eps", sep="/"),
##           width=2.5, height=2.5, horizontal=F)
tikz(file=paste(plotDIRch4, "aovResids1.tex", sep="/"),
     width=2.5, height=2.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
qqmath(~resid(Everg.aov),
       panel = function(x,...) {
           panel.grid()
           panel.qqmath(x,...)
           panel.qqmathline(x,...)
       }, ylab="Residuals", xlab="Unit Normal Quantile"
       )
dev.off()

##postscript(file=paste(plotDIR, "aovResids2.eps", sep="/"),
##           width=3, height=2.5, horizontal=F)
tikz(file=paste(plotDIRch4, "aovResids2.tex", sep="/"),
           width=3, height=2.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
xyplot(sqrt(abs(resid(Everg.aov)))~fitted(Everg.aov),
       panel=function(x,y,...){
           panel.grid()
           panel.xyplot(x, y,...)
           panel.loess(x, y, span=1, col="grey",...)
       }, ylab="Sqrt. Abs. Residualt", xlab="Fitted")
dev.off()

##postscript(file=paste(plotDIR, "aovResids3.eps", sep="/"),
##           width=3, height=2.5, horizontal=F)
tikz(file=paste(plotDIRch4, "aovResids3.tex", sep="/"),
           width=3, height=2.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
xyplot(resid(Everg.aov)~fitted(Everg.aov),
       panel=function(x,y,...){
         panel.grid()
         panel.xyplot(x, y,...)
         panel.abline(0, 0)
                    #        panel.loess(x, y, span=1, col="grey",...)
       }, ylab="Residualt", xlab="Fitted")
dev.off()

Everg.aov2 <- aov((TP)^-0.75 ~ factor(Year), data=TP.reference)
##postscript(file=paste(plotDIR, "aovResids4.eps", sep="/"),
##width=2.5, height=2.5, horizontal=F)
tikz(file=paste(plotDIRch4, "aovResids4.tex", sep="/"),
           width=2.5, height=2.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
qqmath(~resid(Everg.aov2),
       panel = function(x,...) {
         panel.grid()
         panel.qqmath(x,...)
         panel.qqmathline(x,...)
       }, ylab="Residuals", xlab="Unit Normal Quantile"
       )
dev.off()

## anova as a linear model (dummy variables)
anova.data <- data.frame(y=TP.reference$TP,
                         x2=ifelse(TP.reference$Year==1995, 1, 0),
                         x3=ifelse(TP.reference$Year==1996, 1, 0),
                         x4=ifelse(TP.reference$Year==1997, 1, 0),
                         x5=ifelse(TP.reference$Year==1998, 1, 0),
                         x6=ifelse(TP.reference$Year==1999, 1, 0))
anova.lm <- lm(log(y) ~ x2+x3+x4+x5+x6, data=anova.data)
summary(anova.lm)

anova.lm <- lm(log(TP) ~ factor(Year), data=TP.reference)
summary(anova.lm)

anova.lm <- lm(log(TP) ~ factor(Year)-1, data=TP.reference)
summary(anova.lm)

anova.lm <- lm((TP)^0.25 ~ factor(Year), data=TP.reference)
summary(anova.lm)

powerT(TP.reference$TP)

## family versus testwise alpha
anova.p <- t.p <- numeric()
    for (i in 1:1000){
        data.sim <- data.frame(y=rnorm(120),
                               g=rep(1:6, each=20))
        sample.mean <- tapply(data.sim$y, data.sim$g, mean)
        data.sim$g <- ordered(data.sim$g,
                        levels=names(sort(sample.mean)))
        data.sim$g <- as.numeric(data.sim$g)
        anova.p[i] <- summary(aov(y~factor(g),
                              data=data.sim))[[1]][1,5] < 0.05
        t.p[i] <- t.test(y~g, data=data.sim,
                         subset=g==1|g==6)$p.value < 0.05
    }
    print(c(mean(anova.p), mean(t.p)))

anova.p <- t.p <- numeric()
    for (i in 1:1000){
        data.sim <- data.frame(y=rnorm(120),
                               g=rep(1:6, each=20))
        sample.mean <- tapply(data.sim$y, data.sim$g, mean)
        data.sim$g <- ordered(data.sim$g,
                        levels=names(sort(sample.mean)))
        data.sim$g <- as.numeric(data.sim$g)
        anova.p[i] <- summary(aov(y~factor(g),
                              data=data.sim))[[1]][1,5] < 0.05
        t.p[i] <- t.test(y~g, data=data.sim,
                         subset=g==1|g==6)$p.value < 0.0033
    }
    print(c(mean(anova.p), mean(t.p)))


### ANOVA and multiple comparisons ###
# Ellison et al 1996
mangrove.sponge <- read.table(paste(dataDIR, "completespongedata.txt",
                                    sep="/"), header=T)

##postscript(file=paste(plotDIR, "mangroveData1.eps", sep="/"),
##           width=4, height=3, horizontal=F)
tikz(file=paste(plotDIRch4, "mangroveData1.tex", sep="/"),
           width=4, height=3, standAlone=F)
par(mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(RootGrowthRate ~ Treatment, data=mangrove.sponge,
     xlab="Treatment", ylab="Root Growth Rate (mm/day)", cex.axis=0.75)
dev.off()

##trellis.device(postscript, file=paste(plotDIR, "mangroveData2.eps", sep="/"),
##               width=4, height=4, horizontal=F)
tikz(file=paste(plotDIRch4, "mangroveData2.tex", sep="/"),
     width=4, height=4, standAlone=F)
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
trellis.par.set(theme = canonical.theme("postscript", col=FALSE),
                trellis.par.temp)
qqmath(~RootGrowthRate|Treatment, data=mangrove.sponge,
       panel=function(x, ...){
         panel.qqmath(x, ...)
         panel.qqmathline(x, ...)
         panel.grid()
       }, xlab="Unit Norma Quantile", ylab="Root Growth Rate (mm/day)")
dev.off()

mangrove.lm <- lm(RootGrowthRate ~ Treatment, data=mangrove.sponge)
summary.aov(mangrove.lm)
mangrove.lm2 <- lm(RootGrowthRate ~ Treatment*Location, data=mangrove.sponge)
summary.aov(mangrove.lm2)

lm.plots(mangrove.lm)
mangrove.aov <- aov(RootGrowthRate ~ Treatment, data=mangrove.sponge)
mangrove.multcomp <- TukeyHSD(mangrove.aov)

par(mar=c(3, 10, 3, 1), las=1)
plot(TukeyHSD(mangrove.aov))

packages(multcomp)
q2<-glht(mangrove.aov, linfct=mcp(Treatment="Tukey"))
summary(q2)
##postscript(file="mangroveHSD.eps", width=4, height=3, horizontal=F)
tikz(file=paste(plotDIRch4, "mangroveHSD.tex", sep="/"),
     width=4, height=3, standAlone=F)
par(mar=c(3, 10, 3, 1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(q2, main="95\\% family-wise confidence level")
dev.off()

contr <- rbind("F - C" = c(-1, 1, 0, 0),
                 "H - C" = c(-1, 0, 1, 0),
                 "T - C" = c(-1, 0, 0, 1),
                 "S - F" = c(0, -1, 1/2, 1/2),
                 "S - C" = c(-1, 0, 1/2, 1/2))
q3 <-  glht(mangrove.aov, linfct = mcp(Treatment = contr))

plot(q3)

summary(q3, test=adjusted(type=c("none")))

contr2 <- rbind(  "S - F" = c(0, -1, 1/2, 1/2),
                 "S - C" = c(-1, 0, 1/2, 1/2))
q4 <-  glht(mangrove.aov, linfct = mcp(Treatment = contr2))

plot(q4)

summary(q4, test=adjusted(type=c("none")))


mangrove.aov2 <- aov(RootGrowthRate ~ Treatment*Location, data=mangrove.sponge)
q4 <-  glht(mangrove.aov2, linfct = mcp(Treatment = contr))
summary(q4, p.adjust.methods="none")

q5 <- glht(mangrove.aov2, linfct=mcp(Treatment="Tukey"))
summary(q5, p.adjust.methods="none")

p.adjust(0.0033, method="bon", 15)

require(agricolae)
LSD.test,

HSD.test,

pairwise.t.test(mangrove.sponge$RootGrowthRate,
                mangrove.sponge$Treatment, "none")


## removing the two potential outliers
mangrove2 <- mangrove.sponge[-c(14, 39),]
mangrove.aov2 <- aov(RootGrowthRate ~ Treatment*Location, data=mangrove2)
q4 <-  glht(mangrove.aov2, linfct = mcp(Treatment = contr))
summary(q4, test=adjusted(type=c("none")))
pairwise.t.test(mangrove2$RootGrowthRate, mangrove2$Treatment, "none")

contr2 <- rbind(
                "S - F" = c(0, -1, 1/2, 1/2),
                "Trt - C" = c(-1, 1/3, 1/3, 1/3))
q6 <-  glht(mangrove.aov, linfct = mcp(Treatment = contr2))
summary(q6, test=adjusted(type="none"))

t.test(RootGrowthRate ~ Treatment, data=mangrove.sponge,
       subset= Treatment == "Foam"|Treatment=="Control", var.equal=T)

##

#> table(mangrove.sponge$Treatment)
#
#  Control      Foam Haliclona   Tedania
#       21        20        17        14

pooled.sd <- 0.462
contr.sd <- 0.462*sqrt(1/21+1/20)
delta <- 0.59150 - 0.23714

t.st <- delta/(contr.sd)

p.value <- 2*(1-pt(t.st, 20+21-2))
