source("FrontMatter.R")

## Chapter 6
## nonlinear regression

plotDIRch6 <- paste(plotDIR, "chapter6", "figures", sep="/")

# the PCB in fish example
pcb.exp <- nls(pcb ~ pcb0 * exp(-k * (year-1974)), data=laketrout,
                  start=list(pcb0=10, k=0.08))

summary(pcb.exp)
nlm.coef <- coef(pcb.exp)

tikz(file=paste(plotDIRch6, "pcbnls1.tex", sep="/"),
     width=4.5, height=3, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(pcb~year, data=laketrout, col="gray", xlab="Year", ylab="PCB (mg/kg)", cex=0.5)
curve(nlm.coef[1]*exp(-nlm.coef[2]*(x-1974)), lwd=2, add=T)
dev.off()


obj2 <- qqmath(~resid(pcb.exp),
            panel = function(x,...) {
                panel.grid()
                panel.qqmath(x,...)
                panel.qqmathline(x,...)
            }, ylab="Residuals", xlab="Standard Normal Quantile",
#            scales=list(y=list(alternating=2)),
)

tikz(file=paste(plotDIRch6, "pcbnlsDiag2.tex", sep="/"),
         height=3.5, width=4, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
obj2
dev.off()

obj3 <- xyplot(resid(pcb.exp)~fitted(pcb.exp), panel=function(x,y,...){
        panel.grid()
        panel.xyplot(x, y,...)
        panel.abline(0, 0)
        panel.loess(x, y, span=1, col="gray",...)
        }, ylab="Residuals", xlab="Fitted")
## checking for patterns in residuals (independence)
tikz(file=paste(plotDIRch6, "pcbnlsDiag3.tex", sep="/"),
     height=3.5, width=4, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
obj3
dev.off()

obj4 <- xyplot(sqrt(abs(resid(pcb.exp)))~fitted(pcb.exp), panel=function(x,y,...){
        panel.grid()
        panel.xyplot(x, y,...)
        panel.loess(x, y, span=1, col="gray",...)
        }, ylab="Sqrt. Abs. Residuals", xlab="Fitted"
#        scales=list(y=list(alternating=2)),
)
## checking whether the residuals have a constant variance
tikz(file=paste(plotDIRch6, "pcbnlsDiag4.tex", sep="/"),
     height=3.5, width=4, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
obj4
dev.off()
## checking whether residuals are normally distributed

tikz(file=paste(plotDIRch6, "pcbnlsDiag5.tex", sep="/"),
     height=3.5, width=4, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.5,0.125,0), las=1, tck=0.01)
hist(resid(pcb.exp), xlab="Residuals", main="")
dev.off()

### Stow et al 2004-- 4 alternative models

pcb.exp <- nls(log(pcb) ~ log(pcb0) -k * (year-1974), data=laketrout,
                  start=list(pcb0=10, k=0.08))

summary(pcb.exp)
nlm.coef <- coef(pcb.exp)

## alternative 2
pcb.exp2 <- nls(log(pcb) ~ log(pcb0*exp(-k*(year-1974))+pcba),
                data=laketrout, start=list(pcb0=10, k=0.08, pcba=1))
nlm.coef2 <- coef(pcb.exp2)
summary(pcb.exp2)

## alternative 3
pcb.exp3 <- nls(log(pcb) ~ log(pcb01*exp(-k1*(year-1974))+pcb02*exp(-k2*(year-1974))),
                data=laketrout, start=list(pcb01=10, pcb02=2, k1=0.24, k2=0.00002),
                algorithm="port", lower = rep(0, 4))
nlm.coef3 <- coef(pcb.exp3)
summary(pcb.exp3)

## alternative 4
mixedorder <- function(x, b0, k, theta){
    LP1 <- LP2 <- 0
    if(theta==1){
        LP1 <- log(b0) - k*x
    } else {
        LP2 <- log(b0^(1-theta) - k*x*(1-theta))/(1-theta)
    }
    return( LP1 + LP2)
    }
pcb.exp4 <- nls(log(pcb) ~ mixedorder(x=year-1974, pcb0, k, phi), data=laketrout, start=list(pcb0=10, k=0.0024, phi=3.5))
nlm.coef4 <- coef(pcb.exp4)
summary(pcb.exp4)

###postscript(file="pcbnls2.eps", width=4.25, height=3, horizontal=F)
tikz(file=paste(plotDIRch6, "pcbnls2.tex", sep="/"),
     height=3, width=4.5, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(pcb~year, data=laketrout, col=grey(0.5), xlab="Year", ylab="PCB (mg/kg)", cex=0.5)
curve(nlm.coef [1]*exp(-nlm.coef[2]*(x-1974)), lwd=2, add=T, lty=1)
curve(nlm.coef2[1]*exp(-nlm.coef2[2]*(x-1974))+nlm.coef2[3], lwd=2, add=T, lty=2)
curve(nlm.coef3[1]*exp(-nlm.coef3[3]*(x-1974))+nlm.coef3[2]*exp(-nlm.coef3[4]*(x-1974)), lwd=2, add=T, lty=3)
curve((nlm.coef4[1]^(1-nlm.coef4[3]) - nlm.coef4[2]*(x-1974)*(1-nlm.coef4[3]))^(1/(1-nlm.coef4[3])), lwd=2, add=T, lty=4)
legend (x=1995, y=40, legend=1:4, lty=1:4, lwd=2, cex=0.5, bty="n")
dev.off()

exp1.pred <- exp(predict(pcb.exp,  new=data.frame(year=c(2000, 2007))))
exp2.pred <- exp(predict(pcb.exp2, new=data.frame(year=c(2000, 2007))))
exp3.pred <- exp(predict(pcb.exp3, new=data.frame(year=c(2000, 2007))))
exp4.pred <- exp(predict(pcb.exp4, new=data.frame(year=c(2000, 2007))))

logdiff1 <- diff(exp1.pred)

exp1.sim <- sim.nls (pcb.exp, 1000)
betas <- exp1.sim$beta
pred.00<-betas[,1]*exp(-betas[,2]*(2000-1974))
pred.07<-betas[,1]*exp(-betas[,2]*(2007-1974))
percent1 <- 1-pred.07/pred.00
precent1 <- percent1[!is.na(percent1)]

exp2.sim <- sim.nls (pcb.exp2, 1000)
betas <- exp2.sim$beta
pred.00<-betas[,1]*exp(-betas[,2]*(2000-1974))+betas[,3]
pred.07<-betas[,1]*exp(-betas[,2]*(2007-1974))+betas[,3]
percent2 <- 1-pred.07/pred.00
precent2 <- percent2[!is.na(percent2)]

exp3.sim <- sim.nls (pcb.exp3, 1000)
betas <- exp3.sim$beta
pred.00<-betas[,1]*exp(-ifelse(betas[,3]<0, 0, betas[,3])*(2000-1974))+betas[,2]*exp(-ifelse(betas[,4]<0, 0, betas[,4])*(2000-1974))
pred.07<-betas[,1]*exp(-ifelse(betas[,3]<0, 0, betas[,3])*(2007-1974))+betas[,2]*exp(-ifelse(betas[,4]<0, 0, betas[,4])*(2007-1974))
percent3 <- 1-pred.07/pred.00
precent3 <- percent3[!is.na(percent3)]

mixedorder2 <- function(x, b0, k, theta){
    ifelse(theta==1, exp(log(b0) - k*x), exp(log(b0^(1-theta) - k*x*(1-theta))/(1-theta)))
    }

exp4.sim <- sim.nls (pcb.exp4, 1000)
betas <- exp4.sim$beta
pred.00<-mixedorder2 (2000, betas[,1], betas[,2], betas[,3])
pred.07<-mixedorder2 (2007, betas[,1], betas[,2], betas[,3])
percent4 <- 1-pred.07/pred.00
precent4 <- percent4[!is.nan(percent4)]

tikz(file=paste(plotDIRch6, "pcb00to07.tex", sep="/"),
     width=4, height=2.5, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(y=c(1-0.25, 4+0.25), x=c(0,0.5)*100, type="n",
     xlab="\\% change, 2000 to 2007", ylab="Models", axes=F)
                 segments(y0=1:4, y1=1:4,
                          x0=100*c(quantile(percent1, prob=0.025, na.rm=T), quantile(percent2, prob=0.025, na.rm=T),
                              quantile(percent3, prob=0.025, na.rm=T), quantile(percent4, prob=0.025, na.rm=T)),
                         x1=100*c(quantile(percent1, prob=0.975, na.rm=T), quantile(percent2, prob=0.975, na.rm=T),
                              quantile(percent3, prob=0.975, na.rm=T), quantile(percent4, prob=0.975, na.rm=T)))
segments(y0=1:4, y1=1:4, x0=100*c(quantile(percent1, prob=0.25, na.rm=T), quantile(percent2, prob=0.25, na.rm=T),
                              quantile(percent3, prob=0.25, na.rm=T), quantile(percent4, prob=0.25, na.rm=T)),
                         x1=100*c(quantile(percent1, prob=0.75, na.rm=T), quantile(percent2, prob=0.75, na.rm=T),
                              quantile(percent3, prob=0.75, na.rm=T), quantile(percent4, prob=0.75, na.rm=T)),
                         lwd=2.5)
points(y=1:4, x=100*c(quantile(percent1, prob=0.5, na.rm=T), quantile(percent2, prob=0.5, na.rm=T),
                              quantile(percent3, prob=0.5, na.rm=T), quantile(percent4, prob=0.5, na.rm=T)))
abline(v=25)
axis(1)
axis(2, at=1:4)
box()
dev.off()


lake.lm7 <- lm(log(pcb) ~ len.c*I(year-1974), data=laketrout)
display(lake.lm7, 4)

#laketrout <- laketrout[!is.na(laketrout$length),]
#lake.nlm1 <- nls(lnpcb ~ beta0 + beta1*(year-1974) + beta2*len.c + delta*pmax(0, len.c-lth),
#start=list(beta0=1.6, beta1= - 0.08, beta2=0.07, delta=0.01, lth=0), data=laketrout, algorithm="port",
#na.action=na.omit)


## Figure 6.8
## piecewise linear model
  beta1 <- 0.75
  beta2 <- 2.5
  eps <- 5
  alpha1 <- 2
        x1 <- -eps
        x2 <- +eps
        b <- (x2*beta1-x1*beta2)/(x2-x1)
        cc <- (beta2-b)/(2*x2)
        a <- alpha1+beta1*x1-b*x1-cc*x1^2
        alpha2 <- - beta2*x2 +(a + b*x2 + cc*x2^2)

## postscript(file="hockeystick.eps", height=3.5, width=4, horizontal=F)
tikz(file=paste(plotDIRch6, "hockeystick.tex", sep="/"),
     height=3.25, width=4.125, standAlone=F)
plot(y=c(2.5*(35-20)+2, 0.75*(5-20)+2), x=c(35, 5), type="n", xlab="x", ylab="y")
segments(x0=c(5, 25), x1=c(15, 35), y0=c(0.75*(5-20)+2, 2.5*(25-20)+2), y1=c(0.75*(15-20)+2, 2.5*(35-20)+2))
lines(seq(15,25,,100), a + b*(seq(15,25,,100)-20) + cc*(seq(15,25,,100)-20)^2)
points(x=c(15, 25), y=c(0.75*(15-20)+2, 2.5*(25-20)+2))
segments(x0=c(15, 25), x1=c(20, 20), y0=c(0.75*(15-20)+2, 2.5*(25-20)+2), y1=c(2, 2), lty=4)
dev.off()



#hockey <- function(x,alpha1,beta1,delta,brk,eps=diff(range(x))/100) {
#
#       ## alpha1 is the intercept of the left line segment
#       ## beta1 is the slope of the left line segment
#       ## beta2 is the slope of the right line segment
#       ## brk is location of the break point
#       ## 2*eps is the length of the connecting quadratic piece
#
#       ## reference: Bacon & Watts "Estimating the Transition Between
#       ## Two Intersecting Straight Lines", Biometrika, 1971
#
#        beta2 <- beta1 + delta
#        x1 <- brk-eps
#        x2 <- brk+eps
#        b <- (x2*beta1-x1*beta2)/(x2-x1)
#        cc <- (beta2-b)/(2*x2)
#        a <- alpha1+beta1*x1-b*x1-cc*x1^2
#        alpha2 <- - beta2*x2 +(a + b*x2 + cc*x2^2)
#
#        lebrk <- (x <= brk-eps)
#        gebrk <- (x >= brk+eps)
#        eqbrk <- (x > brk-eps & x < brk+eps)
#
#        result <- rep(0,length(x))
#        result[lebrk] <- alpha1 + beta1*x[lebrk]
#        result[eqbrk] <- a + b*x[eqbrk] + cc*x[eqbrk]^2
#        result[gebrk] <- alpha2 + beta2*x[gebrk]
#        result
#}


lake.nlm1 <- nls(log(pcb) ~  hockey(len.c, beta0, beta1, delta, theta),
                 start=list(beta0=1.6, beta1=0.07, delta=0.03, theta=0),
                 data=laketrout,
                 na.action=na.omit)
summary(lake.nlm1)

lake.nlm1 <- nls(log(pcb) ~  hockey(length, beta0, beta1, delta, theta),
start=list(beta0=.6, beta1=0.07, delta=0.03, theta=60), data=laketrout,
na.action=na.omit)
summary(lake.nlm1)
lake1.coef<-coef(lake.nlm1)

lake1.sim <- sim.nls (lake.nlm1, 1000)
betas <- lake1.sim$beta

logPCB.mean  <- betas[,1] +
    (betas[,2] + betas[,3]*(60>betas[,4]))*(60-betas[,4])
pred.PCB<-exp(rnorm(1000, logPCB.mean, lake1.sim$sigma))
pred.log.95CI  <- quantile(log(pred.PCB), prob=c(0.025, 0.975))
theta <- quantile(lake1.sim$beta[,4], prob=c(0.025, 0.975))

##postscript(paste(base, "lake1sim.eps", sep="/"), height=3.5, width=4.75, horizontal=F)
tikz(file=paste(plotDIRch6, "lake1sim.tex", sep="/"),
     height=3.5, width=4.75, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,.125,0), las=1, tck=0.01)
plot(log(pcb)~length, data=laketrout, xlab="Length (cm)",ylab="log PCB", type="n")
for (i in 1:200)
curve(hockey(x, betas[i,1], betas[i,2], betas[i,3], betas[i,4]), col=grey(0.6) , add=T)
curve(hockey(x, lake1.coef[1], lake1.coef[2], lake1.coef[3], lake1.coef[4]),lwd=3, add=T)
points(laketrout$length, log(laketrout$pcb), col=grey(0.4), cex=0.5)
segments(x0=60, x1=60, y0=pred.log.95CI[1], y1=pred.log.95CI[2])
segments(y0=-1.5, y1=-1.5, x0=theta[1], x1=theta[2], col=grey(0.5))
dev.off()


hist(pred.PCB)


lake.nlm2 <- nls(log(pcb) ~  beta1*(year-1974) + hockey(length, beta0, beta2, delta, theta),
start=list(beta0=.6, beta1= - 0.08, beta2=0.07, delta=0.03, theta=60), data=laketrout,
na.action=na.omit)
summary(lake.nlm2)

lake2.coef<- coef(lake.nlm2)

##postscript(paste(base, "lake2sim.eps", sep="/"), height=3.5, width=4.75, horizontal=F)
#pdf(paste(base, "lake2sim.pdf", sep="/"), height=3.5, width=4.75)
tikz(file=paste(plotDIRch6, "lake2sim.tex", sep="/"),
     height=3.5, width=4.75, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,.125,0), las=1, tck=0.01)
plot(log(pcb)~length, data=laketrout, xlab="Length (cm)",ylab="log PCB", type="n")
curve(lake2.coef[2]*(0)+hockey(x, lake2.coef[1], lake2.coef[3], lake2.coef[4], lake2.coef[5]), col=1 , add=T)
points(laketrout$length[laketrout$year<1984], log(laketrout$pcb[laketrout$year<1984]), col=1, cex=0.5)
curve(lake2.coef[2]*(10)+hockey(x, lake2.coef[1], lake2.coef[3], lake2.coef[4], lake2.coef[5]), col=2 , add=T)
points(laketrout$length[laketrout$year>=1984&laketrout$year<1994],
       log(laketrout$pcb[laketrout$year>=1984&laketrout$year<1994]), col=2, cex=0.5)
curve(lake2.coef[2]*(20)+hockey(x, lake2.coef[1], lake2.coef[3], lake2.coef[4], lake2.coef[5]), col=3 , add=T)
points(laketrout$length[laketrout$year>=1994],
       log(laketrout$pcb[laketrout$year>=1994]), col=3, cex=0.5)
curve(lake2.coef[2]*(30)+hockey(x, lake2.coef[1], lake2.coef[3], lake2.coef[4], lake2.coef[5]), col=4 , add=T)
legend(x=30, y=3.5, col=1:4, legend=seq(1974,2004, 10), lty=1, cex=0.5)
dev.off()

##postscript(paste(base, "lake2bw.eps", sep="/"), height=3.5, width=4.75, horizontal=F)
#pdf(paste(base, "lake2sim.pdf", sep="/"), height=3.5, width=4.75)

tikz(file=paste(plotDIRch6, "lake2bw.tex", sep="/"),
     height=3.5, width=4.75, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,.125,0), las=1, tck=0.01)
plot(log(pcb)~length, data=laketrout, xlab="Length (cm)",ylab="log PCB", type="n")
curve(lake2.coef[2]*(0)+hockey(x, lake2.coef[1], lake2.coef[3], lake2.coef[4], lake2.coef[5]),
      col=1, lty=1, add=T, lwd=2)
points(laketrout$length[laketrout$year<1984], log(laketrout$pcb[laketrout$year<1984]),
      col=1, pch="1", cex=0.5)
curve(lake2.coef[2]*(10)+hockey(x, lake2.coef[1], lake2.coef[3], lake2.coef[4], lake2.coef[5]),
      col=grey(0.5), lty=2, add=T, lwd=2)
points(laketrout$length[laketrout$year>=1984&laketrout$year<1994],
       log(laketrout$pcb[laketrout$year>=1984&laketrout$year<1994]),
       col=grey(0.5), pch="2", cex=0.5)
curve(lake2.coef[2]*(20)+hockey(x, lake2.coef[1], lake2.coef[3], lake2.coef[4], lake2.coef[5]),
       col=1, lty=3, add=T, lwd=2)
points(laketrout$length[laketrout$year>=1994],
       log(laketrout$pcb[laketrout$year>=1994]),
       col=1, pch="3", cex=0.5)
curve(lake2.coef[2]*(30)+hockey(x, lake2.coef[1], lake2.coef[3], lake2.coef[4], lake2.coef[5]),
       col=1, lty=4, add=T, lwd=2)
legend(x=30, y=3.5, col=c(1, grey(0.5), 1, 1), legend=seq(1974,2004,
                                                   10), lty=1:4,
       pch=c(1:3, ""), cex=0.75, bty="n")
dev.off()

## US lilac data
USLilac <- read.csv(paste(dataDIR, "NAmlilac.csv", sep="/"))
#> USLilac[1:10,]
#    STID Year Ptype FirstLeaf FirstBloom
#1  20309 1957     2       999         67
#2  20309 1958     2       999         90

USLilac$type <- "Syringa chinensis clone"
USLilac$type[USLilac$Ptype==2] <-"Syringa vulgaris"


USLilac$FirstLeaf[USLilac$FirstLeaf==999] <- NA
USLilac$FirstBloom[USLilac$FirstBloom==999] <- NA

## Figure 6.11


## keep stations with at least 30 years of data

keep <- !is.na(match(USLilac$STID,
                     as.numeric(names(table(USLilac$STID))[table(USLilac$STID)>=30])))
uslilacs <- USLilac[keep,]
stlist <- sort(unique(uslilacs$STID))
for (i in 1:length(stlist)){
    stid <- stlist[i]
    obj <- xyplot(FirstBloom~Year, data=USLilac, panel=function(x,y,...){
        panel.xyplot(x,y,col="gray",...)
        panel.loess(x,y,span=0.75,lwd=2,...)},
        xlab="Year", ylab="First Bloom", subset=STID==stid,
        sub=paste("Station", stid))
   # postscript(file=paste("st",stid,".eps", sep=""), width=4.5, height=3, horizontal=F)
    print(obj)
   # dev.off()
}

stlist2 <- c(354147, 456974, 456624, 426357)
obj2<-list()
for (i in 1:4){
    stid <- stlist2[i]
    obj2[[i]] <- xyplot(FirstBloom~Year, data=USLilac, panel=function(x,y,...){
        panel.xyplot(x,y,col="gray",...)
        panel.loess(x,y,span=0.75,lwd=2,...)},
        xlab=paste("Station", stid), ylab="", subset=STID==stid)
#   postscript(file=paste("st",stid,".eps", sep=""), width=4.5, height=3, horizontal=F)
#   print(obj)
#   dev.off()
}
##postscript(paste(plotDIR, "uslilacs1.eps", sep="/"),
##           height=5, width=5.5, horizontal=F)
tikz(paste(plotDIRch6, "uslilacs1.tex", sep="/"),
     height=5, width=5.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75

print(obj2[[1]], position=c(0, 0, 0.525, 0.525), more=T)
print(obj2[[2]], position=c(0.475, 0, 1, 0.525), more=T)
print(obj2[[3]], position=c(0, 0.475, 0.525, 1), more=T)
print(obj2[[4]], position=c(0.475, 0.475, 1, 1), more=F)
dev.off()

## Figure 6.12
postscript(file=paste(plotDIRch6,"uslilacsData.eps", sep="/"),
           width=4.5, height=3, horizontal=F)
##tikz(file=paste(plotDIRch6, "uslilacsData.tex", sep="/"),
##     height=3, width=4.75, standAlone=F)
xyplot(FirstBloom~Year, data=USLilac, panel=function(x,y,...){
    panel.xyplot(x,y, col="gray",...)
    panel.loess(x,y,span=0.5,lwd=3,...)},
    xlab="Year", ylab="First Bloom")
dev.off()

## other plots not used in the book
#postscript(file="USlilac2.eps", width=5.5, height=3, horizontal=F)
xyplot(FirstBloom~Year|type, data=USLilac, panel=function(x,y,...){
    panel.xyplot(x,y,...)
    panel.loess(x,y,col="red",span=0.5,lwd=3,...)},
    xlab="Year", ylab="First Bloom")
#dev.off()

#postscript(file="USlilac3.eps", width=4, height=3, horizontal=F)
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0))
plot(1956:2003, table(USLilac$Year), xlab="Year", ylab="Sample Size", type="h", lwd=3)
#dev.off()

#postscript(file="JPlilacN.eps", width=4, height=3, horizontal=F)
#par(mar=c(3,3,1,1), mgp=c(1.5,.5,0))
#plot(1996:2006, table(Lilac$Year), xlab="Year", ylab="Sample Size", type="h", lwd=3)
#dev.off()


temp <- USLilac[USLilac$STID==stlist2[1],]
lilacs.lm1 <- nls( FirstBloom ~ hockey(Year, beta0, beta1, delta, theta),
        start=list(beta0=100, beta1=  0,
                delta=-0.1, theta=1980),
        data=temp, na.action=na.omit)
temp <- USLilac[USLilac$STID==stlist2[2],]
lilacs.lm2 <- nls( FirstBloom ~ hockey(Year, beta0, beta1, delta, theta),
        start=list(beta0=150, beta1=  0,
                delta=-0.1, theta=1980),
        data=temp, na.action=na.omit, control=list(maxiter=100))
temp <- USLilac[USLilac$STID==stlist2[3],]
lilacs.lm3 <- nls( FirstBloom ~ hockey(Year, beta0, beta1, delta, theta),
        start=list(beta0=120, beta1=0.1,
                delta=-1, theta=1980),
        data=temp, na.action=na.omit, control=list(maxiter=100))
temp <- USLilac[USLilac$STID==stlist2[4],]
lilacs.lm4 <- nls( FirstBloom ~ hockey(Year, beta0, beta1, delta, theta),
        start=list(beta0=120, beta1=  0,
                delta=-0.5, theta=1980),
        data=temp, na.action=na.omit, control=list(maxiter=100))

#uslilacs.ST <- read.csv(paste(dataDIR, "NAmlilacSTID.csv", sep="/"),
#                        header=T)

#uslilacs.ST[uslilacs.ST$ID%in% stlist2,]

#354147	IMNAHA	        	OR	45.34	-116.5	564.02
#456974	REPUBLIC		WA	48.39	-118.44	792.68
#456624	PORT ANGELES		WA	48.07	-123.26	60.98
#426357	OAK CITY		UT	39.23	-112.2	1547.26

### Starting values

## examples of fpl

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

tikz(file=paste(plotDIRch6, "toledo8.1.tex", sep="/"),
     height=3, width=4.25, standAlone=F)
par(mar=c(3, 3, 0.5,0.25), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
plot(Abs8.1.1 ~ stdConc8.1, xlab="MC Concentration",
     ylab="Optical Density")
dev.off()

rng <- range(Abs8.1.0)
drng <- diff(rng)
prop <- (Abs8.1.0 - rng[1] + 0.05*drng)/(1.1*drng) 
coef(lm(I(log(prop/(1-prop))) ~ log(stdConc8.1+0.001)))

coef( nls(Abs8.1.0 ~ (al1-al4)/(1+(stdConc8.1/al3)^al2)+al4, start=list(al1=1.2, al2=0.64, al3=0.16, al4=.15)))
summary( nls(Abs8.1.0 ~ cbind(1, 1/(1+(stdConc8.1/al3)^al2)), start=list(al2=0.64, al3=0.16), algorithm="plinear"))

nlm.plots <- function(nls.obj){
    obj1<-xyplot(fitted(nls.obj)~(fitted(nls.obj)+resid(nls.obj)),
                 panel = function(x, y,...) {
                     panel.xyplot(x, y,...)
                     panel.abline(0,1, lty=1,...)
                     panel.loess(x,y, span=1.0,lty=2,...)
                     panel.grid()
                 },ylab="Fitted",xlab="Observed")
    ## checking whether the predicted is in greement with the observed

    obj2<-qqmath(~resid(nls.obj),
                 panel = function(x,...) {
                     panel.grid()
                     panel.qqmath(x,...)
                     panel.qqmathline(x,...)
                 }, ylab="Residuals", xlab="Standard Normal Quantile"
                 )
    ## checking whether residuals are normally distributed

    obj3<-xyplot(resid(nls.obj)~fitted(nls.obj), panel=function(x,y,...){
        panel.grid()
        panel.xyplot(x, y,...)
        panel.abline(0, 0, col="gray")
        panel.loess(x, y, span=1, lwd=2,...)
    }, ylab="Residuals", xlab="Fitted")
    ## checking for patterns in residuals (independence)

    obj4<-xyplot(sqrt(abs(resid(nls.obj)))~fitted(nls.obj), panel=function(x,y,...){
        panel.grid()
        panel.xyplot(x, y,...)
        panel.loess(x, y, span=1, lwd=2,...)
    }, ylab="Sqrt. Abs. Residuals", xlab="Fitted")
    ## checking whether the residuals have a constant variance
    print(obj1, position = c(0.0, 0.0, 0.5, 0.5), more = T)
    print(obj2, position = c(0.5, 0.0, 1.0, 0.5), more = T)
    print(obj3, position = c(0.0, 0.5, 0.5, 1), more = T)
    print(obj4, position = c(0.5, 0.5, 1.0, 1), more = F)
    invisible()
}

## self-start function of a 4-parameter logistic function:
## y = al4+(al1-al4)/(1+(x/al3)^al2)
## the function SSfpl is for a different form of fpl:
## A+(B-A)/(1+exp((xmid-input)/scal))
## these two functions are the same (input = log(x))

## the mean function
fplModel <- function(input, al1, al2, al3, al4){
    .x <- input+0.0001
    .expr1 <- (.x/al3)^al2
    .expr2 <- al1-al4
    .expr3 <- 1 + .expr1
    .expr4 <- .x/al3
    .value <- al4 + .expr2/.expr3
    .grad <- array(0, c(length(.value), 4L),
                   list(NULL, c("al1","al2","al3","al4")))
    .grad[,"al1"] <- 1/.expr3
    .grad[,"al2"] <- -.expr2*.expr1*log(.expr4)/.expr3^2
    .grad[,"al3"] <- .expr1*.expr2*(al2/al3)/.expr3^2
    .grad[,"al4"] <- .expr1/(1+.expr1)
    attr(.value, "gradient") <- .grad
    .value
}

## initial values
fplModelInit <- function(mCall, LHS, data){
    xy <- sortedXyData(mCall[["input"]], LHS, data)
   if (nrow(xy) < 5) {
        stop("too few distinct input values to fit a four-parameter logistic")
    }
    rng <- range(xy$y)
    drng <- diff(rng)
    xy$prop <- (xy$y-rng[1]+0.05*drng)/(1.1*drng)
    xy$logx <- log(xy$x+0.0001)
    ir <- as.vector(coef(lm(I(log(prop/(1-prop))) ~ logx, data=xy)))
    pars <- as.vector(coef(nls(y~cbind(1, 1/(1+(x/exp(lal3))^al2)),
                               data=xy,
                               start=list(al2=-ir[2], lal3=-ir[1]/ir[2]),
                               algorithm="plinear")))
    value <- c(pars[4]+pars[3], pars[1], exp(pars[2]), pars[3])
    names(value) <- mCall[c("al1","al2","al3","al4")]
    value
}

SSfpl2 <- selfStart(fplModel, fplModelInit, c("al1","al2","al3","al4"))




toledo <- data.frame(stdConc=rep(stdConc8.1, 6),
                     Abs=c(Abs8.1.0,Abs8.1.1,Abs8.1.2,
                         Abs8.2.0,Abs8.2.1,Abs8.2.2),
                     Test=rep(1:6, each=12))

tm1 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==1,])

tm2 <- nls(Abs ~ SSfpl(log(stdConc), A, B, xmid, scal),
           data=toledo[toledo$Test==1&toledo$stdConc>0,])


tikz(file=paste(plotDIRch6, "toledoTM1resid.tex", sep="/"),
     height=6, width=6, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
nlm.plots(tm1)
dev.off()

tikz(file=paste(plotDIRch6, "toledoTM2resid.tex", sep="/"),
     height=6, width=6, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
nlm.plots(tm2)
dev.off()

tm12 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==2&toledo$stdConc>0,])

tm22 <- nls(Abs ~ SSfpl(log(stdConc), A, B, xmid, scal),
           data=toledo[toledo$Test==2&toledo$stdConc>0,])

tm13 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==3&toledo$stdConc>0,])

tm23 <- nls(Abs ~ SSfpl(log(stdConc), A, B, xmid, scal),
           data=toledo[toledo$Test==3&toledo$stdConc>0,])

tm14 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==4&toledo$stdConc>0,])

tm24 <- nls(Abs ~ SSfpl(log(stdConc), A, B, xmid, scal),
           data=toledo[toledo$Test==4&toledo$stdConc>0,])

tm15 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==5&toledo$stdConc>0,])

tm25 <- nls(Abs ~ SSfpl(log(stdConc), A, B, xmid, scal),
           data=toledo[toledo$Test==5&toledo$stdConc>0,])

tm16 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==6&toledo$stdConc>0,])

tm26 <- nls(Abs ~ SSfpl(log(stdConc), A, B, xmid, scal),
           data=toledo[toledo$Test==6&toledo$stdConc>0,])


## old model fit
TM1 <- nls(Abs ~ (A-D)/(1+(stdConc/C)^B)+D,
           control=list(maxiter=200), data=toledo[toledo$Test==1,],
           start=list(A=0.1,B=-1.2,C=0.16,D=1.2))
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



#### Smoothing ####

##postscript(file="smoother1.eps", width=3.5, height=3, horizontal=F)
tikz(paste(plotDIRch6, "smoother1.tex", sep="/"),
     width=3.5, height=3, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
X<-sort(c(35, seq(30, 90, 5)))
Y<-numeric()
for (i in 1:length(X)){
    Y[i]=mean(log(laketrout$pcb[laketrout$length>X[i]-5 & laketrout$length<X[i]+5]), na.rm=T)
}
plot(log(pcb)~length,data=laketrout, type="n", xlab="Fish Length (cm)", ylab="Log PCB")
points(log(pcb)~length,data=laketrout, col="gray", cex=0.5)
polygon(x=c(40, 40, 50, 50), y=c(-2, 4, 4, -2), col="gray", border=NA)
points(log(pcb)~length,data=laketrout, subset=length>40 & length < 50, cex=0.5)
points(x=45, y=mean(log(laketrout$pcb[laketrout$length>40 & laketrout$length < 50]), na.rm=T), pch=16, cex=1.5)
lines(X, Y)
for (i in 1:length(X)){
    Y[i]=mean(log(laketrout$pcb[laketrout$length>X[i]-10 & laketrout$length<X[i]+10]), na.rm=T)
}
lines(X, Y, lty=5)
box()
dev.off()

## loess example

pcb.data <- laketrout[!is.na(laketrout$length)&laketrout$pcb>0,]
oo <- order(pcb.data$length)
pcb.data <- pcb.data[oo,]
## at length=60 the 255th number in length
## n = 646, the upper bound is 255+646*0.25 = 417 or 68.5cm
##
x.range <- c(60-9, 60+9)
pcb.subset <- pcb.data[170:291,]
local.lm <- lm(log(pcb) ~ length, data=pcb.subset)
pcb.loess <- loess.smooth (y=log(pcb.data$pcb), x=pcb.data$length, data=pcb.data, degree=1, span=0.5)

##postscript(file="smoother2.eps", width=3.5, height=3, horizontal=F)
tikz(paste(plotDIRch6, "smoother2.tex", sep="/"),
     width=3.5, height=3, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(log(pcb.data$pcb)~pcb.data$length, xlab="Fish Length (cm)", ylab="Log PCB", type="n")
points(log(pcb.data$pcb)~pcb.data$length, col=grey(0.5), cex=0.5)
polygon(x=x.range[c(1,1,2,2)], y=c(-2, 4, 4, -2), col="gray", border=NA)
lines(pcb.loess$y~pcb.loess$x, lwd=2)
curve(coef(local.lm)[1] + coef(local.lm)[2]*x, add=T, xlim=x.range)
abline(v=60)
points(log(pcb)~length,data=laketrout, subset=length>x.range[1] & length < x.range[2],
       cex=0.25, col=grey(0.3))
points(x=60, y=predict(local.lm, new=data.frame(length=60)), pch=16, cex=1.25)
dev.off()

#######################
### additive models ###
#######################
##postscript(file="lmgraph1.eps", height=2.25, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch6, "lmgraph1.tex", sep="/"),
     height=2.25, width=4.5, standAlone=F)
par(mfrow=c(1,2), mgp=c(1.15,0.125,0), las=1, tck=0.01, mar=c(3,3,1,1))
plot(c(0,2), c(-3, 1), type="l", xlab="$x_1$", ylab="$\\beta_0 + \\beta_1 x_1$")
abline(h=seq(-3,1,0.1), v=seq(0,2,0.1), col="gray")
plot(c(0,4), c(0, -3), type="l", xlab="$x_2$", ylab="$\\beta_2 x_2$")
abline(v=seq(0,4,0.1), h=seq(-3,0,0.1), col="gray")
dev.off()

##postscript(file="lmgraph2.eps", height=2.25, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch6, "lmgraph2.tex", sep="/"),
     height=2.25, width=4.5, standAlone=F)
par(mfrow=c(1,2), mgp=c(1.25,0.125,0), las=1, tck=0.01, mar=c(3,3,1,1))
plot(c(0,2), c(-3, 1), type="l", xlab="$x_1$", ylab="$\\beta_0+\\beta_1x_1$")
abline(h=seq(-3,1,0.1), v=seq(0,2,0.1), col="gray")
plot(log(seq(0.1,100,,100)), log(seq(0.1,100,,100)), type="l",
     xlab="$x_2$", ylab="$\\beta_2\\log(x_2)$", axes=F)
axis(2)
axis(1, at=log(c(0.1,1,10,50,100)), labels=c("0.1","1","10","50","100"))
abline(v=log(c(seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10))), h=seq(-3,6,0.25), col="gray")
box()
dev.off()

## postscript(file="lmgraph3.eps", height=2.25, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch6, "lmgraph3.tex", sep="/"),
     height=2.25, width=4.5, standAlone=F)
par(mfrow=c(1,2), mgp=c(1.25,0.125,0), las=1, tck=0.01, mar=c(3,3,1,1))
plot(c(0,2), c(-3, 1), type="l", xlab="$x_1$", ylab="$\\beta_0+\\beta_1x_1$")
abline(h=seq(-3,1,0.1), v=seq(0,2,0.1), col="gray")
plot((seq(0.1,100,,100)), log(seq(0.1,100,,100)), type="l",
     xlab="$x_2$", ylab="$\\beta_2\\log(x_2)$", axes=F)
axis(2)
axis(1)
abline(v= seq(10,100,10), h=seq(-3,6,0.25), col="gray")
box()
dev.off()

## using PCB in fish
packages(gam)
pcbGamlo <- gam(log(pcb) ~ lo(length, span=0.75, degree=1)+lo(year, span=0.75), data=laketrout)

plot(pcbGamlo, se=T, rug=T, resid=T, pch=1, cex=0.25, ask=T)
plot(pcbGamlo, se=T, rug=T, resid=F, scale=0, pch=1, cex=0.25, ylab="s(Year)", xlab="Year")
detach(package:gam)

packages(mgcv)
pcbGam1 <- gam(log(pcb) ~ s(length)+s(year, fx=T, k=4), data=laketrout)
summary(pcbGam1)

 ##postscript("pcbGam1.eps", height=2.25, width=4.5, horizontal=F)
tikz(paste(plotDIRch6, "pcbGam1.tex", sep="/"),
     height=2.25, width=4.5, standAlone=F)
par(mfrow=c(1,2), mar=c(3,3,1,0.5), mgp=c(1.5,0.5,0))
 plot(pcbGam1, select=1, se=T, rug=T, resid=F, scale=0, pch=1, cex=0.25,ylab="s(Length)", xlab="Length")
 plot(pcbGam1, select=2, se=T, rug=T, resid=F, scale=0, pch=1, cex=0.25,ylab="s(Year)", xlab="Year")
 dev.off()
detach(package:mgcv)

packages(gam)
pcbGam2 <- gam(log(pcb)~s(length) + s(year), data=laketrout)
pcbGam3 <- gam(log(pcb)~s(length) + s(year, 8), data=laketrout)
pcbGam4 <- gam(log(pcb)~s(length) + s(year, 2), data=laketrout)
detach(package:gam)

tikz(paste(plotDIRch6, "pcbGam2.tex", sep="/"),
     height=2, width=5, standAlone=F)
 par(mfrow=c(1,3), mar=c(3,3,0.5,0.25), mgp=c(1.25,0.125,0))
 plot.gam(pcbGam4, se=T, rug=T, terms="s(year, 2)", ylab="s(year, df=2)", xlab="Year")
 plot.gam(pcbGam2, se=T, rug=T, terms="s(year)", ylab="s(year, df=4)", xlab="Year")
 plot.gam(pcbGam3, se=T, rug=T, terms="s(year, 8)", ylab="s(year, df=8)", xlab="Year")
dev.off()

nadb <- read.csv(paste(dataDIR, "nadb.csv", sep="/"), header=T)
nadb$logPLI <- log10(nadb$PLI)
nadb$logTPIn<- log10(nadb$TPIn)
nadb$logTPOut<- log10(nadb$TPOut+1)
nadb$logHLR <- log10(nadb$HLR+1)

tikz(paste(plotDIRch6, "nadbGAM1.tex", sep="/"),
     height=5, width=5, standAlone=F)
par(mar=c(2,2,1,1))
pairs(nadb[,1:5])
dev.off()


tikz(paste(plotDIRch6, "nadbGAM2.tex", sep="/"),
     height=3.5, width=5.5, standAlone=F)
par(mgp=c(1.25,0.125,0), mar=c(3,3,1,1), las=1, tck=0.01)
plot(TPOut ~ PLI, data=nadb, log="x", axes=F,
     xlab="Input TP Loading", ylab="Output TP Concentration")
axis(1, at = c(0.01, 0.1,1, 10, 100, 1000),
     labels=c("0.01","0.1","1","10","100","1000"))
axis(2)
box()
dev.off()

##postscript(paste(plotDIR, "nadbGAM3.eps", sep="/"), height=3.5, width=5.5, horizontal=F)
tikz(paste(plotDIRch6, "nadbGAM3.tex", sep="/"),
     height=3.5, width=5.5, standAlone=F)
par(mgp=c(1.5,0.5,0), mar=c(3,3,1,1))
plot(TPOut ~ PLI, data=nadb, log="xy", axes=F,
     xlab="Input TP Loading", ylab="Output TP Concentration")
axis(1, at = c(0.01, 0.1,1, 10, 100, 1000),
     labels=c("0.01","0.1","1","10","100","1000"))
axis(2)
box()
dev.off()

coplot(logTPOut ~ logTPIn|logHLR, data=nadb, given.v=co.intervals(nadb$logHLR, 4, 0.25), rows=1, panel=panel.smooth)
coplot(logTPOut ~ logHLR|logTPIn, data=nadb, given.v=co.intervals(nadb$logTPIn, 4, 0.25), rows=1, panel=panel.smooth)

require(mgcv)

nadbGam1 <- gam(logTPOut ~ s(logPLI)+s(logTPIn)+s(logHLR), data=nadb)
tikz(paste(plotDIRch6, "nadbGam4.tex", sep="/"),
     height=2, width=5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,0.5,0.25), mgp=c(1.15,0.125,0), las=1, tck=0.01)
plot(nadbGam1, select=1, se=T, rug=T, resid=T, scale=0, pch=16, cex=0.25)
plot(nadbGam1, select=2, se=T, rug=T, resid=T, scale=0, pch=16, cex=0.25)
plot(nadbGam1, select=3, se=T, rug=T, resid=T, scale=0, pch=16, cex=0.25)
dev.off()


nadbGam1.5 <- gam(logTPOut ~ s(logPLI, fx=T, k=4)+s(logTPIn, fx=T, k=4)+s(logHLR, fx=T, k=4), data=nadb)
summary(nadbGam1.5)

tikz(paste(plotDIRch6, "nadbGam4_5.tex", sep="/"),
     height=2, width=5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,0.5,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(nadbGam1.5, select=1, se=T, rug=T, resid=T, scale=0, pch=16, cex=0.25)
plot(nadbGam1.5, select=2, se=T, rug=T, resid=T, scale=0, pch=16, cex=0.25)
plot(nadbGam1.5, select=3, se=T, rug=T, resid=T, scale=0, pch=16, cex=0.25)
dev.off()


#postscript(paste(plotDIR, "nadbGam4.eps", sep="/"), height=3, width=5, horizontal=F)
##tikz(paste(plotDIRch6, "nadbGam4.tex", sep="/"),
#     height=3, width=5, standAlone=F)
par(mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(1.5,0.5,0), pty="s")
plot(nadbGam1.5, select=1, se=T, rug=T, resid=T, scale=0, pch=1)
plot(nadbGam1.5, select=2, se=T, rug=T, resid=T, scale=0, pch=1)
#dev.off()

nadbGam2 <- gam(log(TPOut) ~ s(log(PLI))+s(log(TPIn)), data=nadb, subset=TPOut>0)

par(mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(1.5,0.5,0))
plot(nadbGam2, select=1, se=T, rug=T, resid=T)
plot(nadbGam2, select=2, se=T, rug=T, resid=T)

nadbGam3 <- gam(logTPOut ~ s(var1=logTPIn, var2=logHLR), data=nadb)
summary(nadbGam3)

tikz(paste(plotDIRch6, "nadbGam5.tex", sep="/"),
     height=3.5, width=3.5, standAlone=F)
par(mar=c(4,4,3,1), pty="s")
plot(nadbGam3, select=1, resid=T, mgp=c(1.25,0.125,0), pch=1, cex=0.5, se=F)
dev.off()

tikz(paste(plotDIRch6, "nadbGam6.tex", sep="/"),
     height=3.5, width=3.5, standAlone=F)
par(mar=c(0.5,0.5,0.5,0.5), pty="s", mgp=c(1.25, 0.125,0), las=1, tck=0.01)
plot(nadbGam3, select=1, rug=T, resid=T, pers=T)
dev.off()

nadbGam4 <- gam(logTPOut ~ s(logPLI), data=nadb)
summary(nadbGam4)


tikz(paste(plotDIRch6, "nadbGam7.tex", sep="/"),
     height=2.5, width=3, standAlone=F)
par(mar=c(4,4,1,1), mgp=c(1.5,0.125,0), las=1, tck=0.01)
plot(nadbGam4, select=1, resid=T, mgp=c(1.5,0.5,0), pch=1, cex=0.5, se=T)
dev.off()

## STL and median polishing

## median polish
median.polish.ts <- function(data.ts, ylab="", plt=T){
# median polishing for missing value imputation
medpolish(matrix(data.ts, ncol=12, byrow=T), eps=0.001, na.rm=T)->temp.2w
print(names(temp.2w))
year.temp <- rep(seq(start(data.ts)[1], end(data.ts)[1]), each=12)
month.temp <- rep(1:12, length(seq(start(data.ts)[1], end(data.ts)[1])))
# plotting median polishing results
if (plt){
  par(mfrow=c(2,1))
  plot(seq(start(data.ts)[1], end(data.ts)[1]),
       temp.2w$overall+temp.2w$row, type="l",
       xlab="Year", ylab=ylab, main="De-seasonalized Trend")
  plot(seq(1,12), temp.2w$overall+temp.2w$col, type="l",
       xlab="Month", ylab=ylab, main="Seasonal Changes")
}
data.ts[is.na(data.ts)]<-temp.2w$overall +
  temp.2w$row[year.temp[is.na(data.ts)]-start(data.ts)[1]+1]+
    temp.2w$col[month.temp[is.na(data.ts)]]
invisible(data.ts)
}
#### end of median polishing for missing value imputation


co2 <- read.table(paste(dataDIR, "co2.txt", sep="/"),
                  header=T)

## the CO2 data sets needs additional manipulation
##
## the original data file has two extra columns that are not necessary
co2 <- co2[,2:13]  ## keeping only the monthly CO2 data
co2[co2<0] <- NA   ## replacing missing values (coded as -99) with NA

co2 <- ts(as.vector(t(as.matrix(co2))), start=c(1959, 1), end=c(2003,12), freq=12)
## converting the data into a time series
co2<-median.polish.ts(co2)
## imputing missing values.
## running median polish generates two figures -- long term trend and seasonal pattern


## with a time series object, you can plot it using plot:
#postscript("co2ts.eps", height=3.5, width=5, horizontal=F)

## Figure 6.27
tikz(file=paste(plotDIRch6, "co2HAts.tex", sep="/"),
           height=3.5, width=4.5, standAlone=F)
par(mar=c(3,3,0.25,0.25), mgp=c(1.5,.125,0), las=1, tck=0.01)  ## set up graphical parameters
                                            ## (margins and locations of ticks and labels)
plot(co2, ylab="CO$_2$ (ppm)", xlab="Year")
dev.off()


stl.rfs <- 
function(data = carbon.dioxide, ss.w = 25, ss.d = 1, fc.w = 120, fc.d = 1, ylab = "Carbon Dioxide (ppm)", aspect = "xy", ...)
{
 strip.background <- trellis.par.get("strip.background")
 strip.background$col <- 0.
 trellis.par.set("strip.background", strip.background)
 strip.shingle <- trellis.par.get("strip.shingle")
 strip.shingle$col <- 0.
 trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
 trellis.par.set("strip.shingle", strip.shingle)
    the.fit <- stl(data, s.window = ss.w, s.degree = ss.d, t.window = fc.w, t.degree = fc.d, ...)
    sfit <- the.fit$time.series[,1]
    tfit <- the.fit$time.series[,2]
    fit.time <- time(data)
    car.subseries <- factor(cycle(data), label = month.abb)

    obj1 <- xyplot(sfit ~ fit.time | car.subseries, layout = c(12, 1), panel = function(x, y)
    {
        panel.xyplot(x, y, type = "l")
        panel.abline(h = mean(y))
    }
    , aspect = aspect, xlab = "Year", ylab = ylab)

    obj2 <- xyplot(tfit ~ fit.time, panel = function(x, y)
    panel.xyplot(x, y, type = "l"), xlab = "", aspect = "xy", ylab = "")

    n <- length(data)
    the.fit.trend <- the.fit$time.series[,2] - mean(the.fit$time.series[,2])
    fit.components <- c(the.fit.trend, the.fit$time.series[,1], the.fit$time.series[,3])
    fit.time <- rep(time(data), 3)
    fit.names <- ordered(rep(c("Trend", "Seasonality", "Residuals"), c(n, n, n)), c("Trend", 
        "Seasonality", "Residuals"))

    obj3 <- xyplot(fit.components ~ fit.time | fit.names, panel = function(x, y)
    {
        panel.grid(h = 5)
        panel.xyplot(x, y, type = "l")
    }
    , aspect=0.75, layout = c(3, 1), ylim = c(-1, 1) * max(abs(fit.components)), xlab = "", ylab = ylab)
    print(obj1, position = c(0, 0, 1, 0.4), more = T)
    print(obj2, position = c(0, 0.33, 1, 0.67), more = T)
    print(obj3, position = c(0, 0.6, 1, 1), more = F)
    invisible(data.frame(trend=the.fit$time.series[,2], season=the.fit$time.series[,1], residual=the.fit$time.series[,3]))
}


stl.rfs2 <- 
function(data = carbon.dioxide, ss.w = 25, ss.d = 1, fc.w = 120, fc.d = 1, ylab = "Carbon Dioxide (ppm)", aspect = "xy", ...)
{
 strip.background <- trellis.par.get("strip.background")
 strip.background$col <- 0.
 trellis.par.set("strip.background", strip.background)
 strip.shingle <- trellis.par.get("strip.shingle")
 strip.shingle$col <- 0.
 trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
 trellis.par.set("strip.shingle", strip.shingle)
    the.fit <- stl(data, s.window = ss.w, s.degree = ss.d, t.window = fc.w, t.degree = fc.d, ...)
    sfit <- the.fit$time.series[,1]
    tfit <- the.fit$time.series[,2]
    fit.time <- time(data)
    car.subseries <- factor(cycle(data), label = month.abb)

    obj1 <- xyplot(sfit ~ fit.time | car.subseries, layout = c(12, 1), panel = function(x, y)
    {
        panel.xyplot(x, y, type = "l")
        panel.abline(h = mean(y))
    }
    , aspect = aspect, xlab = "Year", ylab = ylab)

    obj2 <- xyplot(tfit ~ fit.time, panel = function(x, y)
    panel.xyplot(x, y, type = "l"), xlab = "", aspect = "xy", ylab = "")

    n <- length(data)
    the.fit.trend <- the.fit$time.series[,2] - mean(the.fit$time.series[,2])
    fit.components <- c(the.fit.trend, the.fit$time.series[,1], the.fit$time.series[,3])
    fit.time <- rep(time(data), 3)
    fit.names <- ordered(rep(c("Trend", "Seasonality", "Residuals"), c(n, n, n)), c("Trend", 
        "Seasonality", "Residuals"))

    obj3 <- xyplot(fit.components ~ fit.time | fit.names, panel = function(x, y)
    {
        panel.grid(h = 5)
        panel.xyplot(x, y, type = "l")
    }
    , aspect=0.75, layout = c(3, 1), ylim = c(-1, 1) * max(abs(fit.components)), xlab = "", ylab = ylab)
    print(obj1, position = c(0, 0, 1, 0.6), more = T)
    print(obj3, position = c(0, 0.4, 1, 1), more = F)
    invisible(data.frame(trend=the.fit$time.series[,2], season=the.fit$time.series[,1], residual=the.fit$time.series[,3]))
}

## STORET data at 
storetDIR <- paste(dataDIR, "storet", sep="/")

upper.neuse <- read.table(paste(storetDIR,
                                "Data_upp_20080509_220930_RegResults.txt", sep="/"),
                          sep="~", header=T)
middl.neuse <- read.table(paste(storetDIR,
                                "Data_mid_20080509_221058_RegResults.txt", sep="/"),
                          sep="~", header=T)
lower.neuse <- read.table(paste(storetDIR,
                                "Data_low_20080509_221156_RegResults.txt", sep="/"),
                          sep="~", header=T)
J417 <- upper.neuse[upper.neuse$Station.ID=="J4170000       " , ]

J417.alt <- read.table(paste(storetDIR,
                             "Data_417_20080513_002701_RegResults.txt", sep="/"),
                       sep="~", header=T)
### fecal coliform has a long record

J417.FecalColiform <- J417[J417$Character=="Fecal Coliform",]
J417.FecalColiform$Date  <- as.Date(J417.FecalColiform$Activity.Start,  format="%Y-%m-%d %H:%M:%S")
names(J417.FecalColiform)[24] <- "Value"
plot(Value ~ Date, data=J417.FecalColiform)

J417.TP <- J417[J417$Character=="Phosphorus as P",]
J417.TP$Date <- as.Date(J417.TP$Activity.Start,  format="%Y-%m-%d %H:%M:%S")
 names(J417.TP)[24] <- "Value"
plot(Value ~ Date, data=J417.TP)

J417.TKN <- J417[J417$Character=="Nitrogen, Kjeldahl",]
J417.TKN$Date <- as.Date(J417.TKN$Activity.Start, format="%Y-%m-%d %H:%M:%S")
names(J417.TKN)[24] <- "Value"
plot(Value ~ Date, data=J417.TKN)

## calculating monthly means and median polishing

TKN <- TP <- FecalColiform <- rep(NA, 12*(2007-1970))
k <- 0
for (i in 1971:2007){ ## year
  for (j in 1:12){ ## month
     k <- k+1
     temp <- as.numeric(format(J417.TKN$Date, "%m"))==j & as.numeric(format(J417.TKN$Date, "%Y"))==i
     if (sum(temp)>0) TKN[k] <- mean(J417.TKN$Value[temp], na.rm=T)

     temp <- as.numeric(format(J417.TP$Date, "%m"))==j & as.numeric(format(J417.TP$Date, "%Y"))==i
     if (sum(temp)>0) TP[k] <- mean(J417.TP$Value[temp])

     temp <- as.numeric(format(J417.FecalColiform$Date, "%m"))==j &
                      as.numeric(format(J417.FecalColiform$Date, "%Y"))==i
     if (sum(temp)>0) FecalColiform[k] <- mean(J417.FecalColiform$Value[temp])
     }
  }
  

TKN.ts <- ts(TKN, start=c(1971,1), end=c(2007,12), freq=12)
TP.ts <- ts(TP, start=c(1971,1), end=c(2007,12), freq=12)
FecalColiform.ts <- ts(FecalColiform, start=c(1971,1), end=c(2007,12), freq=12)
FCtable <- data.frame(matrix(round(FecalColiform.ts,1), ncol=12, byrow=T))
names(FCtable) <- month.abb
row.names(FCtable) <- 1971:2007

## Data table on page 210
write.table(FCtable, file="FCts.txt", sep="&")

medpolish(matrix(TP.ts, ncol=12, byrow=T), eps=0.001, na.rm=T)

TKN.2w<-median.polish.ts(TKN.ts, plt=F)
TKN.STL <- stl.rfs2(data=na.omit(TKN.2w),
                    aspect=2.25, ylab="TKN",
                    ss.w = 21, ss.d = 1, fc.w = 101, fc.d = 1)

TP.2w<-median.polish.ts(TP.ts, plt=F)
TP.STL <- stl.rfs2(data=na.omit(TP.2w),
                   aspect=2.25, ylab="TKN",
                   ss.w = 25, ss.d = 1, fc.w = 121, fc.d = 1)

FecalColiform.2w<-median.polish.ts(FecalColiform.ts, plt=F)
FecalColiform.2w[is.na(FecalColiform.2w)] <-
  median(FecalColiform.2w[!is.na(FecalColiform.2w)])
FecalColiform.STL <- stl.rfs2(data=(FecalColiform.2w),
                              aspect=2.25, ylab="Fecal Coliform",
                              ss.w = 21, ss.d = 1, fc.w = 101, fc.d = 1,
                              na.action=na.omit, robust=T)

plot(stl(FecalColiform.2w, s.w = "per", robust = TRUE, na.action=na.omit))
## in logarithmic scale
logTKN <- logTP <- logFecalColiform <- rep(NA, 12*(2007-1970))
k <- 0
for (i in 1971:2007){ ## year
  for (j in 1:12){ ## month
     k <- k+1
     temp <- as.numeric(format(J417.TKN$Date, "%m"))==j & as.numeric(format(J417.TKN$Date, "%Y"))==i
     if (sum(temp)>0) logTKN[k] <- mean(log(J417.TKN$Value[temp]))

     temp <- as.numeric(format(J417.TP$Date, "%m"))==j & as.numeric(format(J417.TP$Date, "%Y"))==i
     if (sum(temp)>0) logTP[k] <- mean(log(J417.TP$Value[temp]))

     temp <- as.numeric(format(J417.FecalColiform$Date, "%m"))==j &
                       as.numeric(format(J417.FecalColiform$Date, "%Y"))==i
     if (sum(temp)>0) logFecalColiform[k] <- mean(log(J417.FecalColiform$Value[temp]+1))
     }
  }
  

logTKN.ts <- ts(logTKN, start=c(1971,1), end=c(2007,12), freq=12)
logTP.ts <- ts(logTP, start=c(1971,1), end=c(2007,12), freq=12)
logFecalColiform.ts <- ts(logFecalColiform, start=c(1971,1), end=c(2007,12), freq=12)

logTKN.2w<-median.polish.ts(logTKN.ts, plt=F)
logTKN.STL <- stl.rfs2(data=na.omit(logTKN.2w), aspect=2.25,
                       ylab="Log TKN", ss.w = 21, ss.d = 1, fc.w = 111, fc.d = 1,
                       robust=T)

logTP.2w<-median.polish.ts(logTP.ts, plt=F)


logFecalColiform.2w<-median.polish.ts(logFecalColiform.ts, plt=F)
logFecalColiform.2w[is.na(logFecalColiform.2w)] <-
  mean(logFecalColiform.2w[!is.na(logFecalColiform.2w)])

## Figure 6.28
##postscript(file=paste(plotDIR, "fecalTS.eps", sep="/"),
##           height=2.5, width=3.75, horizontal=F)
tikz(file=paste(plotDIRch6, "fecalTS.tex", sep="/"),
           height=2.5, width=3.75, standAlone=F)
par(mar=c(3,3,.25,.25), mgp=c(1.5,.5,0), tck=-0.025, las=1)
plot.ts(logFecalColiform.2w, xlab="Year", ylab="Log F. Coliform (\\#/100ml)",
        cex=0.75, lwd=0.5)
dev.off()

## Figure 6.29
##postscript(file=paste(plotDIR, "fecalSTLlog.eps", sep="/"),
##           height=4.5, width=5, horizontal=F)
tikz(file=paste(plotDIRch6, "fecalSTLlog.tex", sep="/"),
           height=4.5, width=5, standAlone=F)
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.75 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
trellis.par.temp$axis.text$cex <- 0.7
trellis.par.temp$axis.components$bottom$tck <- 0.5
trellis.par.temp$axis.components$left$tck <- 0.5
trellis.par.temp$axis.components$top$tck <- 0.5
trellis.par.temp$axis.components$right$tck <- 0.5
trellis.par.set(trellis.par.temp)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
logFecalColiform.STL <- stl.rfs2(data=na.omit(logFecalColiform.2w),
                                 aspect=2.25, ylab="Log F. Coliform",
                                 ss.w = 21, ss.d = 1, fc.w = 225, fc.d = 0)
dev.off()

## Figure 6.30
##postscript(file=paste(plotDIR, "TPSTLlog.eps", sep="/"),
##           height=4.5, width=5, horizontal=F)
tikz(file=paste(plotDIRch6, "TPSTLlog.tex", sep="/"),
           height=4.5, width=5, standAlone=F)
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.75 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
trellis.par.temp$axis.text$cex <- 0.7
trellis.par.temp$axis.components$bottom$tck <- 0.5
trellis.par.temp$axis.components$left$tck <- 0.5
trellis.par.temp$axis.components$top$tck <- 0.5
trellis.par.temp$axis.components$right$tck <- 0.5
trellis.par.set(trellis.par.temp)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
logTP.STL <- stl.rfs2(data=na.omit(logTP.2w), aspect=2.25,
                      ylab="Log TP", ss.w = 17, ss.d = 1, fc.w = 101, fc.d = 1)
dev.off()

### compare to ES&T paper
logTKN.short <- rep(NA, 12*(1998-1970))
k <- 0
for (i in 1971:1998){ ## year
  for (j in 1:12){ ## month
     k <- k+1
     temp <- as.numeric(format(J417.TKN$Date, "%m"))==j & as.numeric(format(J417.TKN$Date, "%Y"))==i
     if (sum(temp)>0) logTKN.short[k] <- mean(log(J417.TKN$Value[temp]))
     }
  }
  
logTKN.tsS <- ts(logTKN.short, start=c(1971,1), end=c(1998,12), freq=12)

logTKN.2wS<-median.polish.ts(logTKN.tsS, plt=F)
logTKN.STLshort <- stl.rfs2(data=na.omit(logTKN.2wS), aspect=2.25, ylab="Log TKN", ss.w = 21, ss.d = 1, fc.w = 111, fc.d = 1, robust=T)

   
## Figure 6.31

##postscript(file=paste(plotDIR, "tknTRND.eps", sep="/"), height=5, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch6, "tknTRND.tex", sep="/"),
     height=5, width=4.5, standAlone=F)
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.75 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
trellis.par.temp$axis.text$cex <- 0.7
trellis.par.temp$axis.components$bottom$tck <- 0.5
trellis.par.temp$axis.components$left$tck <- 0.5
trellis.par.temp$axis.components$top$tck <- 0.5
trellis.par.temp$axis.components$right$tck <- 0.5
trellis.par.set(trellis.par.temp)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))

the.fit.S <- stl(na.omit(logTKN.2wS), s.window = 21, s.degree = 1,
                 t.window = 111, t.degree = 1, robust=T)
tfit <- the.fit.S$time.series[,2]
fit.time <- time(na.omit(logTKN.2wS))
obj2 <- xyplot(tfit ~ fit.time, xlim=c(1971, 2002.5), panel = function(x, y)
               panel.xyplot(x, y, type = "l"), xlab = "", aspect = 0.5, ylab = "")

the.fit.L <- stl(na.omit(logTKN.2w), s.window = 21, s.degree = 1,
                 t.window = 111, t.degree = 1, robust=T)
tfit <- the.fit.L$time.series[,2]
fit.time <- time(na.omit(logTKN.2w))
obj1 <- xyplot(tfit ~ fit.time, xlim=c(1971, 2002.5), panel = function(x, y)
               panel.xyplot(x, y, type = "l"), xlab = "", aspect = 0.5, ylab = "")

obj3 <- xyplot(na.omit(logTKN.2w) ~ fit.time, xlim=c(1971, 2002.5), panel = function(x, y)
               panel.xyplot(x, y, type = "l"), xlab = "", aspect = 0.5, ylab = "")


print(obj1, position = c(0, 0, 1, 0.4), more = T)
print(obj2, position = c(0, 0.3, 1, 0.7), more = T)
print(obj3, position = c(0, 0.6, 1, 1), more = F)
dev.off()

NChurricanes <- read.csv(paste(dataDIR, "NChurricanes.csv", sep="/"), header=T)
NChurricanes$Year<-ordered(NChurricanes$Year, levels=1857:2005)

