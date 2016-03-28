source("FrontMatter.R")

## Chapter 9
plotDIRch9 <- paste(plotDIR, "chapter9", "figures", sep="/")
## simulation
packages(rv)
## residual

##postscript(file=paste(plotDIR, "residsim.eps", sep="/"), height=2, width=3, horizontal=F)
tikz(file=paste(plotDIRch9, "residsim.tex", sep="/"),
     height=2, width=3, standAlone=F)
par(mgp=c(1.2,0.5,0), tck=-0.02, mar=c(3,3,1,1), las=1, tck=0.01)
curve(dnorm(x, 0, 0.25), from=-2.5, to=2.5, lty=2, xlab="", ylab="")
curve(dnorm(x,0,1), from=-2.5, to=2.5, add=T)
box()
abline(v=0.5, col="gray", lwd=2)
dev.off()

## Bayesian p-value
#laketrout <- laketrout[laketrout$pcb>exp(-2)&laketrout$length>0,]
#lake.lm1 <- lm(log(pcb) ~ I(year-1974), data=laketrout)
#display(lake.lm1, 4)

laketrout <- read.csv(paste(dataDIR,"laketrout2.csv", sep="/"),
                      header=T)
## dividing by 2.54 converts 60 cm to inches
## single predictor
laketrout$length <- laketrout$length*2.54
laketrout$size<- "small"
laketrout$size[laketrout$length>60] <- "large"

laketrout$large<- 0
laketrout$large[laketrout$length>60] <- 1

laketrout$len.c <- laketrout$length-mean(laketrout$length, na.rm=T)
laketrout <- laketrout[laketrout$pcb>exp(-2)&laketrout$length>0,]

lake.lm1 <- lm(log(pcb) ~ I(year-1974), data=laketrout)
display(lake.lm1, 4)

tikz(file=paste(plotDIRch9, "lengthYear.tex",sep="/"),
           height=3, width=4.5, standAlone=F)
par(mar=c(3,3,1,0.5), las=1, tck=0.01, mgp=c(1.25,0.125,0))
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
xyplot(I(length)~year, data=laketrout, subset=pcb>exp(-2)&length>0,
       ylab="Length (cm)", xlab="Year")
dev.off()

## illustrating Bayesian p-value
tikz(file=paste(plotDIRch9, "illusims.tex", sep="/"),
     height=2.75, width=3.75, standAlone=F)
par(mar=c(3, 3, 2, 0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
curve(1+2*x, xlim=c(0,10), xlab="$x$", ylab="$y$")
abline(v=c(2,5,8), col="gray")
y1 <- seq(-3,3,,100)+1+2*2
x1 <- dnorm(y1, 1+2*2, 1)
lines(2+x1/max(x1), y1)
points(2, qnorm(0.95, 1+2*2, 1))
y11 <- seq(qnorm(0.95, 1+2*2, 1),3+1+2*2, ,25)
x11 <- 2+dnorm(y11, 1+2*2, 1)/max(x1)
polygon(x=c(2, x11, 2), y=c(y11[1], y11, y11[25]), col="gray") 

y2 <- seq(-3,3,,100)+1+2*5
x2 <- dnorm(y2, 1+2*5, 1)
lines(5+x2/max(x2), y2)
points(5, qnorm(1-0.43, 1+2*5, 1))
y22 <- seq(qnorm(1-0.43, 1+2*5, 1),3+1+2*5, ,25)
x22 <- 5+dnorm(y22, 1+2*5, 1)/max(x2)
polygon(x=c(5, x22, 5), y=c(y22[1], y22, y22[25]), col="gray") 

y3 <- seq(-3,3,,100)+1+2*8
x3 <- dnorm(y3, 1+2*8, 1)
lines(8+x3/max(x3), y3)
points(8, qnorm(1-0.91, 1+2*8, 1))
y33 <- seq(qnorm(1-0.91, 1+2*8, 1),3+1+2*8, ,25)
x33 <- 8+dnorm(y33, 1+2*8, 1)/max(x3)
polygon(x=c(8, x33, 8), y=c(y33[1], y33, y33[25]), col="gray") 

axis(3, at=c(2,5,8), label=c(0.05,0.43, 0.91))
dev.off()

packages(rv)
setnsims(5000)
lake1.sim <- posterior(lake.lm1)

pred.lake1 <- rvnorm(1, lake1.sim$beta[1]+lake1.sim$beta[2]*((1974:2000) - 1974),
                     lake1.sim$sigma)

B.pvalue <- list()
for (i in 1974:2000){
 B.pvalue[[i-1973]] <- Pr(pred.lake1[i-1973] >= log(laketrout$pcb[laketrout$year==i]))
}

n.sims<-5000
sim.results <- sim(lake.lm1, n.sims=1000)
predict.PCB07 <- exp(sim.results@coef[,1] +
                     sim.results@coef[,2]*(2007-1974) +
                     0.5* sim.results@sigma)

predict.PCB00 <- exp(sim.results@coef[,1] +
                     sim.results@coef[,2]*(2000-1974) +
                     0.5 * sim.results@sigma)
percentages <- 1-predict.PCB07/predict.PCB00
tikz(file=paste(plotDIRch9, "pcbreduction.tex", sep="/"),
     width=3, height=2.5, standAlone=F)
par(mar=c(4,3,1,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
hist(percentages*100, xlab="Predicted \\% PCB reduction", ylab="", prob=T, main="")
dev.off()


##postscript(paste(plotDIR, "lakesimBPhist.eps", sep="/"),
##    height=2.25, width=3, horizontal=F)
tikz(paste(plotDIRch9, "lakesimBPhist.tex", sep="/"),
    height=2.25, width=3, standAlone=F)
par(mgp=c(1.25,0.125,0), mar=c(3,3,1,1), las=1, tck=0.01)
hist(unlist(B.pvalue), xlab="Tail Area", ylab="", main="")
dev.off()

##postscript(paste(plotDIR, "lakesimBPbox.eps", sep="/"), height=2.5, width=4, horizontal=F)
tikz(paste(plotDIRch9, "lakesimBPbox.tex", sep="/"),
     height=2.5, width=4, standAlone=F)
par(mgp=c(1.25,0.125,0), tck= 0.01, mar=c(3,3,1,1), las=1)
boxplot(B.pvalue, names=1974:2000, xlab="Year", ylab="Tail Area")
dev.off()

pred.lake2 <- rvnorm(1, lake1.sim$beta[1]+lake1.sim$beta[2]*(na.omit(laketrout$year) - 1974),
                     lake1.sim$sigma)

##postscript(paste(plotDIR, "lakesimStats.eps", sep="/"), width=4.25, height=4., horizontal=F)
tikz(paste(plotDIRch9, "lakesimStats.tex", sep="/"),
     width=4.25, height=4., standAlone=F)
par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.25,0.125,0), tck= 0.01, las=1)
rvhist(quantile(pred.lake2, prob=0.95), xlab="95th percentile", ylab="Frequency", breaks=15,
    main=paste("Tail area: ", round(Pr(quantile(pred.lake2, prob=0.95)>quantile(log(laketrout$pcb), prob=0.95, na.rm=T)), 2)))
abline(v=quantile(log(laketrout$pcb), prob=0.95, na.rm=T), lwd=2)

rvhist(quantile(pred.lake2, prob=0.05), xlab="5th percentile", ylab="Frequancy", breaks=15,
    main=paste("Tail area: ", round(Pr(quantile(pred.lake2, prob=0.05)>quantile(log(laketrout$pcb), prob=0.05, na.rm=T)), 2)))
abline(v=quantile(log(laketrout$pcb), prob=0.05, na.rm=T), lwd=2)

rvhist(mean(pred.lake2), xlab="Mean", ylab="Frequency", breaks=15,
    main=paste("Tail area: ", round(Pr(mean(pred.lake2)>mean(log(laketrout$pcb), na.rm=T)), 2)))
abline(v=mean(log(laketrout$pcb), na.rm=T), lwd=2)

rvhist(median(pred.lake2), xlab="Median", ylab="Frequency", breaks=15,
    main=paste("Tail area: ", round(Pr(median(pred.lake2)>median(log(laketrout$pcb), na.rm=T)), 2)))
abline(v=median(log(laketrout$pcb), na.rm=T), lwd=2)
dev.off()

########################


sparrow <- read.table(paste(dataDIR, "sparrow.txt", sep="/"), header=T)

##postscript(file=paste(plotDIR, "cssparrow1.eps", sep="/"), height=2.5, width=3, horizontal=F)
tikz(file=paste(plotDIRch9, "cssparrow1.tex", sep="/"),
     height=2.5, width=3, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(as.vector(by(sparrow$Bird.Count, sparrow$year, mean)), axes=F, xlab="Year", ylab="Average Bird Counts")
axis(1, at=2:14,labels=sort(unique(sparrow$year))[-1])
axis(2)
box()
dev.off()

spar.glm1 <- glm(Bird.Count ~ factor(year), data=sparrow, family=poisson)
display(spar.glm1)

y <- sparrow$Bird.Count
n <- dim(sparrow)[1]
y.rep <- rpois(n, fitted(spar.glm1))
print(mean(y == 0)) #(0.6929774)
print(mean(y.rep==0)) #(0.5089959)


y.mean <- mean(sparrow$Bird.Count==0)
y.rep.mean <- numeric()
n <- dim(sparrow)[1]
n.sims<-5000
for (i in 1:n.sims){
    y.rep <- rpois(n, predict(spar.glm1,type="response"))
    y.rep.mean [i] <- mean(y.rep==0)
}
##postscript(file=paste(plotDIR, "cssparrow2.eps", sep="/"), height=2.25, width=3, horizontal=F)
tikz(file=paste(plotDIRch9, "cssparrow2.tex", sep="/"),
     height=2.25, width=3, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
hist(y.rep.mean*100, xlab="\\% zeroes",
    xlim=100*range(c(y.rep.mean, y.mean)),
    main="", ylab="")
abline(v=y.mean*100, lwd=2, col="gray")
dev.off()


print(mean(sparrow$Bird.Count==0))
print(mean(y.rep==0))



spar.glm2 <- glm(Bird.Count ~ factor(year), data=sparrow, family=quasipoisson)
display(spar.glm2)

for (i in 1:n.sims){
    y.rep <- rpois(n, exp(predict(spar.glm2)))
    y.rep.mean [i] <- mean(y.rep==0)
}
##postscript(file=paste(plotDIR, "sparrow3.eps", sep="/"), height=4, width=5.5, horizontal=F)
tikz(file=paste(plotDIRch9, "sparrow3.tex", sep="/"),
     height=4, width=5.5, standAlone=F)
hist(y.rep.mean*100, xlab="\\% zeroes",
    xlim=100*range(c(y.rep.mean, y.mean)),
    main="", ylab="")
abline(v=y.mean*100)
dev.off()

sparrow.sims <- sim (spar.glm1, n.sims)
y.rep.reg <- numeric()
    for (i in 1:n.sims){
    yr <- sample(1:14, size=1000, replace=T)
    mu <- sparrow.sims@coef[i, yr]
    y.rep.reg[i] <-
        mean(rpois(1000, exp(mu))==0)
}

##postscript(file=paste(plotDIR, "sparrow3.eps", sep="/"), height=4, width=5.5, horizontal=F)
tikz(file=paste(plotDIRch9, "sparrow3.tex", sep="/"),
     height=4, width=5.5, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
hist(y.rep.reg*100, xlab="\\% zeroes",
    xlim=100*range(c(y.rep.reg, y.mean)),
    main="", ylab="")
abline(v=y.mean*100, col=2)
dev.off()

## Figure 9.9
tikz(file=paste(plotDIRch9, "letterFig.tex", sep="/"),
     height=5.5, width=5, standAlone=F)
par(mfrow=c(3,2), mar=c(1.5, 1.5, 0.5, 0.5), oma=c(2, 2, 0.5, 0.5),
    mgp=c(1.25,0.125,0), tck=0.01)
for (size in c(10, 20, 30, 50, 100, 500)){
    x.temp <- runif(size, 5, 45)
    temp <- data.frame(x=x.temp, y=ifelse(x.temp<25, rnorm(sum(x.temp<25), -1, 0.5), rnorm(sum(x.temp>=25), 0.5, 0.5)))
    plot(temp, ylim=c(-2,2), xlim=c(4,45), xlab="", ylab="")
    text(10, 1.25, paste("n = ", size, sep=""))
    }
mtext("x", side=1, line=1, outer=T)    
mtext("y", side=2, line=1, outer=T)    
dev.off()

## ELISA example

mc <- c(0.167, 0.444, 01.110, 2.220, 5.550)
rOD <- c(0.784, 0.588, 0.373, 0.270,0.202)
stdcrv <- lm(log(mc) ~ rOD)
stdcrv.sims <- posterior(stdcrv, n.sims=5000)

ttt <- sim(stdcrv, 1000)

mc.sims <- exp(rvnorm(1, stdcrv.sims$beta[1]+stdcrv.sims$beta[2]* 0.261, stdcrv.sims$sigma))
mc.sims

## Toledo water crisis data
## standard solution MC concentration
stdConc8.1<- rep(c(0,0.167,0.444,1.11,2.22,5.55), each=2)

## measured OD
Abs8.1.0<-c(1.082,1.052,0.834,0.840,0.625,0.630,
            0.379,0.416,0.28,0.296,0.214,0.218)
Abs8.1.1<-c(1.265,1.153,0.94,0.856,0.591,0.643,
            0.454,0.442,0.454,0.447,0.291,0.29)
Abs8.1.2<-c(1.051,1.143,0.679,0.936,0.657,0.662,
            0.464,0.429,0.32,0.35,0.241,0.263)
Abs8.2.0<-c(1.139,1.05,0.877,0.914,0.627,0.705,
            0.498,0.495,0.289,0.321,0.214,0.231)
Abs8.2.1<-c(1.153,1.149,0.947,0.896,0.627,0.656,
            0.465,0.435,0.33,0.328,0.218,0.226)
Abs8.2.2<-c(1.124,1.109,0.879,0.838,0.61,0.611,
            0.421,0.428,0.297,0.308,0.19,0.203)

toledo <- data.frame(stdConc=rep(stdConc8.1, 6),
                     Abs=c(Abs8.1.0,Abs8.1.1,Abs8.1.2,
                         Abs8.2.0,Abs8.2.1,Abs8.2.2),
                     Test=rep(1:6, each=12))
toledo1 <- toledo[toledo$stdConc>0,]

TM1 <- nls(Abs ~ (A-D)/(1+(stdConc/C)^B)+D,
           control=list(maxiter=200), data=toledo[toledo$Test==1,],
           start=list(A=0.2,B=-1.,C=0.5,D=1.))
TM2 <- nls(Abs ~ (A-D)/(1+(stdConc/C)^B)+D, 
           control=list(maxiter=200), data=toledo[toledo$Test==2,],
           start=list(A=0.2,B=-1.,C=0.45,D=1.))
TM3 <- nls(Abs ~ (A-D)/(1+(stdConc/C)^B)+D, 
           control=list(maxiter=200), data=toledo[toledo$Test==3,],
           start=list(A=0.16,B=-1.12,C=0.45,D=1.06))
TM4 <- nls(Abs ~ (A-D)/(1+(stdConc/C)^B)+D, 
           control=list(maxiter=200), data=toledo[toledo$Test==4,],
           start=list(A=0.16,B=-1.12,C=0.45,D=1.06))
TM5 <- nls(Abs ~ (A-D)/(1+(stdConc/C)^B)+D, 
           control=list(maxiter=200), data=toledo[toledo$Test==5,],
           start=list(A=0.16,B=-1.12,C=0.45,D=1.06))
TM6 <- nls(Abs ~ (A-D)/(1+(stdConc/C)^B)+D, 
           control=list(maxiter=200), data=toledo[toledo$Test==6,],
           start=list(A=0.16,B=-1.12,C=0.45,D=1.06))

sim.nls <- function (object, n.sims=100){
    ## sim.nls:  get posterior simulations of sigma and beta from an nls object
    ##           modified from the function sim in package arm.  The following
    ##           are nearly verbatim from the initial sim function:
    ## Arguments:
    ##
    ##     object:  the output of a call to "nls"
    ##              with n data points and k predictors
    ##     n.sims:  number of independent simulation draws to create
    ##
    ## Output is a list (sigma.sim, beta.sim):
    ##
    ##     sigma.sim:  vector of n.sims random draws of sigma
    ##       (for glm's, this just returns a vector of 1's or else of the
    ##       square root of the overdispersion parameter if that is in the model)
    ##     beta.sim:  matrix (dimensions n.sims x k) of n.sims random draws of beta
    ##
    
    object.class <- class(object)[[1]]
    if (object.class!="nls") stop("not a nls object")
    
    summ <- summary (object)
    coef <- summ$coef[,1:2,drop=FALSE]
    dimnames(coef)[[2]] <- c("coef.est","coef.sd")
    sigma.hat <- summ$sigma
    beta.hat <- coef[,1]
    V.beta <- summ$cov.unscaled
    n <- summ$df[1] + summ$df[2]
    k <- summ$df[1]
    sigma <- rep (NA, n.sims)
    beta <- array (NA, c(n.sims,k))
    dimnames(beta) <- list (NULL, names(beta.hat))
    for (s in 1:n.sims){
        sigma[s] <- sigma.hat*sqrt((n-k)/rchisq(1,n-k))
        beta[s,] <- mvrnorm (1, beta.hat, V.beta*sigma[s]^2)
    }
    return (list (beta=beta, sigma=sigma))
}

test1.sim <- sim.nls(TM1, 4000)
test1.beta <- rvsims(test1.sim$beta)
test1.sigma <- rvsims(test1.sim$sigma)
test1.coef <- coef(TM1)
test1.pred <-
    summary(rvnorm(mean=(test1.beta[1]-test1.beta[4])/
                        (1+(seq(0,5.55,,500)/test1.beta[3])^test1.beta[2])+
                        test1.beta[4], sd=test1.sigma))
    
    test3.sim <- sim.nls(TM3, 4000)
    test3.beta <- rvsims(test3.sim$beta)
    test3.sigma <- rvsims(test3.sim$sigma)
    test3.coef <- coef(TM3)
    test3.pred <-
        summary(rvnorm(mean=(test3.beta[1]-test3.beta[4])/
                         (1+(seq(0,5.55,,500)/test3.beta[3])^test3.beta[2])+
                         test3.beta[4], sd=test3.sigma))
    
    tikz(file=paste(plotDIRch9, "tolUncas.tex", sep="/"), 
         height=6.5, width=4.75, standAlone=F)
    par(mfrow=c(3,1), mgp=c(1.25,0.25,0), las=1, tck=0.01)
    par(mar=c(0,3,3,2.5))
    plot( Abs ~ stdConc, type="n", data=toledo, xlab=" ",
         ylab=" ", axes=F)
    polygon(x=c(seq(0,5.55,,500), rev(seq(0,5.55,,500))),
            y=c(test1.pred[,4],rev(test1.pred[,8])), 
            col=grey(0.75), border=NA) 
    polygon(x=c(seq(0,5.55,,500), rev(seq(0,5.55,,500))),
            y=c(test1.pred[,5],rev(test1.pred[,7])), 
            col=grey(0.5), border=NA) 
    curve((test1.coef[1]-test1.coef[4])/(1+(x/test1.coef[3])^test1.coef[2])+
           test1.coef[4], add=T)
    points(toledo$stdConc[toledo$Test==1],
           toledo$Abs[toledo$Test==1], cex=0.75, pch=16)
axis(1, at=0:5, label=rep("", 6))
axis(2, at=seq(0.2,1.2,0.2), label=rep("", 6))
    axis(3)
    axis(4)
    box()
    text(5, 1.125, "(a)")
    abline(h=0.35, col=gray(0.4), lty=5)
    abline(v=1.5, col=gray(0.4), lty=5)
    par(mar=c(1.5,3,1.5,2.5))
    plot( Abs ~ stdConc, type="n", data=toledo, xlab=" ",
         ylab="OD", axes=F)
    polygon(x=c(seq(0,5.55,,500), rev(seq(0,5.55,,500))),
            y=c(test3.pred[,4],rev(test3.pred[,8])), 
            col=grey(0.75), border=NA) 
    polygon(x=c(seq(0,5.55,,500), rev(seq(0,5.55,,500))),
            y=c(test3.pred[,5],rev(test3.pred[,7])),
            col=grey(0.5), border=NA) 
    curve((test3.coef[1]-test3.coef[4])/(1+(x/test3.coef[3])^test3.coef[2])+
              test3.coef[4], add=T)
    points(toledo$stdConc[toledo$Test==3],
           toledo$Abs[toledo$Test==3], cex=0.75, pch=16)
axis(3, at=0:5, label=rep("", 6))
axis(1, at=0:5, label=rep("", 6))
    axis(2)
axis(4, at=seq(0.2,1.2,0.2), label=rep("", 6))
    box()
    text(5, 1.125, "(b)")
    
    par(mar=c(3,3,0,2.5))
    plot( Abs ~ stdConc, type="p", data=toledo, 
          xlab="MC Concentration ($\\mu$g/L)",
          ylab=" ", axes=F)
    aa <- coef(TM1)
    curve((aa[1]-aa[4])/(1+(x/aa[3])^aa[2])+aa[4], add=T, lwd=2)
    aa <- coef(TM2)
    curve((aa[1]-aa[4])/(1+(x/aa[3])^aa[2])+aa[4], add=T, col=grey(0.5))
    aa <- coef(TM3)
    curve((aa[1]-aa[4])/(1+(x/aa[3])^aa[2])+aa[4], add=T, col=grey(0.5))
    aa <- coef(TM4)
    curve((aa[1]-aa[4])/(1+(x/aa[3])^aa[2])+aa[4], add=T, col=grey(0.5))
    aa <- coef(TM5)
    curve((aa[1]-aa[4])/(1+(x/aa[3])^aa[2])+aa[4], add=T, col=grey(0.5))
    aa <- coef(TM6)
    curve((aa[1]-aa[4])/(1+(x/aa[3])^aa[2])+aa[4], add=T, col=grey(0.5))
    text(5, 1.125, "(c)")
axis(1)
axis(2, at=seq(0.2,1.2,0.2), label=rep("", 6))
axis(3, at=0:5, label=rep("", 6))
    axis(4)
    box()
    
    dev.off()

## the changepoint program for normal response:
chngp.nonpar <- function(infile)
{
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
    chngp <- mean(mx[vi == min(vi)])
    return(chngp)
}

my.bootCIs<-
function (x, nboot, theta, ..., alpha = c(0.05, 0.95)) 
{ # calculating confidence interval using both percentiles and the BCa method
  # modified from Efron's function 'bacnon'
    n <- length(x)
    thetahat <- theta(x, ...)
    bootsam <- matrix(sample(x, size = n * nboot, replace = TRUE), 
        nrow = nboot)
    thetastar <- apply(bootsam, 1, theta, ...)
    confpoints.percent <- quantile(thetastar, alpha)
    if (sum(thetastar < thetahat)==0) return(rep(as.vector(confpoints.percent), 2))
    else {
    z0 <- qnorm(sum(thetastar < thetahat)/nboot)
        u <- rep(0, n)
        for (i in 1:n) {
        u[i] <- theta(x[-i], ...)
        }
        uu <- mean(u) - u
        acc <- sum(uu * uu * uu)/(6 * (sum(uu * uu))^1.5)
        zalpha <- qnorm(alpha)
        tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))
        ooo <- trunc(tt * nboot)
        if (ooo[1]>0) confpoints <- sort(thetastar)[ooo]
        else confpoints <- rep(NA, length(alpha))
        return(as.vector(c(confpoints, confpoints.percent)))
    }
}
#################
## try vector/matrix

size.v<-c(10, 20, 30, 50, 100,500)
n.sims <- 5000
n.boot <- 1500
percent1<-percent2<-numeric()

date()
for (i in 1:length(size.v)){
print(i)
   size <- size.v[i]
   x.data <- matrix(0, nrow=n.sims, ncol=size, byrow=F)
   x.data <- t(apply(x.data, 1, FUN=function(x)runif(length(x), 5, 45)))
   y.data <- t(apply(x.data, 1, FUN=function(x)ifelse(x<25, rnorm(sum(x<25), -1), rnorm(sum(x>=25), 0.5))))
   ## note: need only loop through sample size.  Use apply to generate confidence intervals

   all.data <- cbind(x.data, y.data)
   results <- apply(all.data, 1, FUN=function(x, sz){
       temp <- data.frame(X=x[1:(length(x)/2)], 
                    Y=x[(length(x)/2+1):length(x)])
       my.bootCIs(1:sz, nboot=n.boot, theta=function(x, infile){
           chngp.nonpar(infile[x, ])}, 
               infile=temp, alpha=c(0.025, 0.975))
           }, sz=size)
   percent1[i] <- mean(apply(results[1:2,], 2, 
       FUN=function(x)x[1]<25 & x[2]>25), na.rm=T)
   percent2[i] <- mean(apply(results[3:4,], 2, 
       FUN=function(x)x[1]<25 & x[2]>25))
    
}
date()
