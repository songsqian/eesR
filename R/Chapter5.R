source("FrontMatter.R")

### linear models ####
##################################
## Chapter 5 linear regression  ##
##################################
plotDIRch5 <- paste(plotDIR, "chapter5", "figures", sep="/")

## PCB in fish

##### Lake trout ####
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
laketrout$inches<-laketrout$length/2.54

laketrout <- laketrout[laketrout$pcb>exp(-2)&laketrout$length>0,]

## 60 cm is a threshold,
##

pcb.loess <- loess.smooth (y=log(laketrout$pcb), x=laketrout$length, degree=1, span=0.75)

tikz(file=paste(plotDIRch5, "pcbfishsize.tex", sep="/"),
     height=2.5, width=4.5, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.75,0.125,0), tck=0.01,las=1)
plot(pcb ~ length, data=laketrout, log="y", ylab="PCB Concentrations (mg/kg)",
     xlab="Fish Length (cm)")
lines(exp(pcb.loess$y)~pcb.loess$x, lwd=2, col="gray")
dev.off()

tikz(file=paste(plotDIRch5, "pcbfishyear.tex", sep="/"),
     height=2.5, width=4.5, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.75,0.125,0), tck=0.01,las=1)
plot(pcb ~ year, data=laketrout, log="y", ylab="PCB Concentrations (mg/kg)",
     xlab="Year")
dev.off()


temp <- data.frame(PCB = c(laketrout$pcb, log(laketrout$pcb)),
                   Large = rep(laketrout$large, 2),
                   Log=rep(c("PCB", "log PCB"), each =dim(laketrout)[1]))
tikz(file=paste(plotDIRch5, "pcbQQ.tex", sep="/"),
     height=2.5, width=5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
qq(Large ~ PCB|Log, data=temp, xlab="Small", ylab="Large",
   scales=list(x=list(relation="free", tck=-0.3),
       y=list(relation="free", tck=-0.3)), aspect=1)
dev.off()

## t-test
pcb.ttest <- t.test(log(pcb) ~ large, data=laketrout)
print(pcb.ttest)
diff(pcb.ttest$estimate)
pcb.ttest$p.value

laketrout$yr <- laketrout$year - 1974
## simple linear regression
lake.lm1 <- lm(log(pcb) ~ I(year-1974), data=laketrout)
display(lake.lm1, 4)
summary(lake.lm1)

lm1.coef<-coef(lake.lm1)
##postscript(paste(plotDIR, "PCBlm1Pred.eps", sep="/"),
##           width=4, height=2.5, horizontal=F)
tikz(paste(plotDIRch5, "PCBlm1Pred.tex", sep="/"),
           width=4.5, height=2.5, standAlone=F)
par(mar=c(3,3,0.25,0.25), mgp=c(1.25,.125,0), las=1, tck=0.01)
plot(pcb ~ year, data=laketrout, log="y",
     ylab="PCB (mg/kg)", xlab="Year", cex=0.5, col=grey(0.5))
curve(exp( lm1.coef[1] + lm1.coef[2]*(x-1974) +0.87^2/2), add=T, lwd=2)
curve(qlnorm(0.025, lm1.coef[1] + lm1.coef[2]*(x-1974), 0.87 ), add=T, lty=2)
curve(qlnorm(0.975, lm1.coef[1] + lm1.coef[2]*(x-1974), 0.87 ), add=T, lty=2)
dev.off()

lake.lm2 <- lm(log(pcb) ~ length+I(year-1974), data=laketrout)
display(lake.lm2, 3)

lake.lm3 <- lm(log(pcb) ~ len.c+I(year-1974), data=laketrout)
display(lake.lm3, 3)
lm3.coef<-coef(lake.lm3)
##postscript(file=paste(plotDIR, "PCBlm2Pred.eps", sep="/"),
##           width=4, height=2.5, horizontal=F)
tikz(file=paste(plotDIRch5, "PCBlm2Pred.tex", sep="/"),
           width=4.5, height=2.5, standAlone=F)
par(mar=c(3,3,0.25,0.25), mgp=c(1.25,.125,0), las=1, tck=0.01)
plot(pcb ~ year, data=laketrout,log="y",
     ylab="PCB (mg/kg)", xlab="Year", cex=0.5, col=grey(0.5))
curve(exp( lm3.coef[1] + lm3.coef[2]*0 + lm3.coef[3]*(x-1974) +0.543^2/2),
      add=T, lwd=2)
curve(exp( lm3.coef[1] + lm3.coef[2]*(-2.6)+lm3.coef[3]*(x-1974) +0.543^2/2),
      add=T, lty=2, lwd=2)
curve(exp( lm3.coef[1] + lm3.coef[2]*(3.4) + lm3.coef[3]*(x-1974) +0.543^2/2),
      add=T, lty="13", lwd=2)
legend(x=1990, y=40, legend=c("large fish","average fish","small fish"),
       lty=c(3, 1, 2), lwd=rep(2, 3), cex=0.75, bty="n")
dev.off()

mn.length <- mean(laketrout$length, na.rm=T)

lake.lm4 <- lm(log(pcb) ~ len.c*I(year-1974), data=laketrout)
display(lake.lm4, 4)


#lake.lm9 <- lm(log(pcb) ~ I(year-1974) + len.c * factor(size), data=laketrout)
#display(lake.lm9, 4)


#lake.lm9 <- lm(log(pcb) ~ I(year-1974) + len.c * factor(size)-1-len.c,
 #              data=laketrout)
#display(lake.lm9, 4)

#lake.lm10 <- lm(log(pcb) ~ I(year-1974)*factor(size) + len.c * factor(size),
 #               data=laketrout)
#display(lake.lm10, 4)

#lake.lm11 <- lm(log(pcb) ~ len.c + factor(size), data=laketrout)
#display(lake.lm11, 4)

#lake.lm12 <- lm(log(pcb) ~ len.c * factor(size), data=laketrout)
#display(lake.lm12, 4)

## plots
## data
##library(lattice)

postscript(file=paste(plotDIR, "PCByear.eps", sep="/"),
           height=3, width=4.5, horizontal=F, family="Times")
#pdf(file=paste(plotDIR, "PCByear.pdf", sep="/"), height=3, width=4.5)
xyplot(log(pcb)~year, data=laketrout,
    panel=function(x,y,...){
        panel.xyplot(x, y, ...)
        panel.lmline(x,y,...)}, subset=pcb>exp(-2)&length>0,
        xlab="Year",ylab="Log PCB")
dev.off()

postscript(file=paste(plotDIR, "PCBlength.eps", sep="/"),
           height=3, width=4.5, horizontal=F)
#pdf(file=paste(plotDIR, "PCBlength.pdf", sep="/"), height=3, width=4.5)
xyplot(log(pcb)~length, data=laketrout,
    panel=function(x,y,...){
        panel.xyplot(x, y, col=grey(0.5), ...)
        panel.lmline(x,y,...)
        panel.loess(x, y, span=3/4, lty=2, ...)}
        , subset=pcb>exp(-2)&length>0, xlab="Length (in)", ylab="Log PCB")
dev.off()

postscript(file=paste(plotDIR, "lengthYear.eps",sep="/"),
           height=3, width=4.5, horizontal=F)
xyplot(I(length)~year, data=laketrout, subset=pcb>exp(-2)&length>0,
       ylab="Length (cm)", xlab="Year")
dev.off()
## residuals

##postscript(file=paste(plotDIR, "PCBresidQQN.eps", sep="/"),
##          height=3, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch5, "PCBresidQQN.tex", sep="/"),
     height=3, width=4.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
qqmath(~resid(lake.lm4),
            panel = function(x,...) {
              panel.grid()
              panel.qqmath(x,...)
              panel.qqmathline(x,...)
            }, ylab="Residuals", xlab="Standard Normal Quantile"
       )
dev.off()
## checking whether residuals are normally distributed

##postscript(file=paste(plotDIR, "PCBresid2.eps", sep="/"),
##          height=3, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch5, "PCBresid2.tex", sep="/"),
           height=3, width=4.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
xyplot(resid(lake.lm4)~fitted(lake.lm4),
       panel=function(x,y,...){
         panel.grid()
         panel.xyplot(x, y, col=grey(0.5),cex=0.5, ...)
         panel.abline(0, 0, lty=2)
         panel.loess(x, y, span=1, ...)
       }, ylab="Residuals", xlab="Fitted")
dev.off()
## checking for patterns in residuals (independence)


##postscript(file=paste(pltoDIR,"PCBresidSL.eps", sep="/"),
##           height=3, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch5,"PCBresidSL.tex", sep="/"),
           height=3, width=4.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
xyplot(sqrt(abs(resid(lake.lm4)))~fitted(lake.lm4),
       panel=function(x,y,...){
         panel.grid()
         panel.xyplot(x, y, col=grey(0.5), cex=0.5, ...)
         panel.loess(x, y, span=1,...)
       }, ylab="Sqrt. Abs. Residualt", xlab="Fitted")
## checking whether the residuals have a constant variance
dev.off()

##postscript(file="PCBcook.eps", height=3, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch5,"PCBcook.tex", sep="/"),
           height=3, width=4.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
xyplot(cooks.distance(lake.lm4) ~ fitted(lake.lm4),
    panel=function(x,y,...){
      panel.xyplot(x,y,...)
      panel.grid()},
    ylab="Cook's Distance", xlab="Fitted")
## checking for influential data points
dev.off()

##postscript(file="PCBrfs.eps", height=3, width=4.75, horizontal=F)
tikz(file=paste(plotDIRch5, "PCBrfs.tex", sep="/"),
     height=3, width=4.75, standAlone=F)
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.75 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
trellis.par.set(trellis.par.temp)
rfs(lake.lm4, aspect=1)
dev.off()

lake.lm5 <- lm(log(pcb) ~ I(year-1974)*len.c + I(len.c^2), 
               data=laketrout)
display(lake.lm5, 4)

lake.lm6 <- lm(log(pcb) ~ I(year-1974)*factor(size) + 
                 len.c * factor(size),
               data=laketrout)
display(lake.lm6, 4)

lake.lm7 <- lm(log(pcb) ~ I(year-1974) + len.c * factor(size),
              data=laketrout)
display(lake.lm7, 4)

lake.lm8 <- lm(log(pcb) ~ I(year-1974) + 
                 len.c * factor(size)-1-len.c,
              data=laketrout)
display(lake.lm8, 4)

tikz(file=paste(plotDIRch5, "PCBresid3.tex", sep="/"),
           height=3, width=4.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
xyplot(resid(lake.lm7)~fitted(lake.lm7),
       panel=function(x,y,...){
         panel.grid()
         panel.xyplot(x, y, col=grey(0.5),cex=0.5, ...)
         panel.abline(0, 0, lty=2)
         panel.loess(x, y, span=1, ...)
       }, ylab="Residuals", xlab="Fitted")
dev.off()
## checking for patterns in residuals (independence)

lm.plots(lake.lm4)

#### 
lake.lm0 <- lm(pcb ~ I(year-1974) + len.c, data=laketrout,
               y=TRUE, qr=TRUE)
        ##postscript(file="pcbBoxCox.eps", height=3, width=3.5, horizontal=F)

## boxcox function extracted with:
## methods("boxcox")
## getAnywhere("boxcox.default")
boxcox.tex <-
    function (object, lambda = seq(-2, 2, 1/10), plotit = TRUE,
              interp = (plotit && (m < 100)), eps = 1/50,
              xlab = expression(lambda), ylab = "log-Likelihood", 
              ...) 
{ ## this is necessary only when using tikz 
    if (is.null(y <- object$y) || is.null(xqr <- object$qr)) 
        stop(gettextf("%s does not have both 'qr' and 'y' components", 
            sQuote(deparse(substitute(object)))), domain = NA)
    if (any(y <= 0)) 
        stop("response variable must be positive")
    n <- length(y)
    y <- y/exp(mean(log(y)))
    logy <- log(y)
    xl <- loglik <- as.vector(lambda)
    m <- length(xl)
    for (i in 1L:m) {
        if (abs(la <- xl[i]) > eps) 
            yt <- (y^la - 1)/la
        else yt <- logy * (1 + (la * logy)/2 * (1 + (la * logy)/3 * 
            (1 + (la * logy)/4)))
        loglik[i] <- -n/2 * log(sum(qr.resid(xqr, yt)^2))
    }
    if (interp) {
        sp <- spline(xl, loglik, n = 100)
        xl <- sp$x
        loglik <- sp$y
        m <- length(xl)
    }
    if (plotit) {
        mx <- (1L:m)[loglik == max(loglik)][1L]
        Lmax <- loglik[mx]
        lim <- Lmax - qchisq(19/20, 1)/2
        dev.hold()
        on.exit(dev.flush())
        plot(xl, loglik, xlab = xlab, ylab = ylab, type = "l", 
            ylim = range(loglik, lim))
        plims <- par("usr")
        abline(h = lim, lty = 3)
        y0 <- plims[3L]
        scal <- (1/10 * (plims[4L] - y0))/par("pin")[2L]
        scx <- (1/10 * (plims[2L] - plims[1L]))/par("pin")[1L]
        text(xl[1L] + scx, lim + scal, " 95\\%", xpd = TRUE)
        la <- xl[mx]
        if (mx > 1 && mx < m) 
            segments(la, y0, la, Lmax, lty = 3)
        ind <- range((1L:m)[loglik > lim])
        if (loglik[1L] < lim) {
            i <- ind[1L]
            x <- xl[i - 1] + ((lim - loglik[i - 1]) * (xl[i] - 
                xl[i - 1]))/(loglik[i] - loglik[i - 1])
            segments(x, y0, x, lim, lty = 3)
        }
        if (loglik[m] < lim) {
            i <- ind[2L] + 1
            x <- xl[i - 1] + ((lim - loglik[i - 1]) * (xl[i] - 
                xl[i - 1]))/(loglik[i] - loglik[i - 1])
            segments(x, y0, x, lim, lty = 3)
        }
    }
    invisible(list(x = xl, y = loglik))
}

tikz(file=paste(plotDIRch5, "pcbBoxCox.tex", sep="/"),
     height=3, width=3.5, standAlone=F)
par(mgp=c(1.25,0.125,0), mar=c(3,3,0.25,0.25), las=1, tck=0.02, cex=0.75)
### boxcox(lake.lm0) ## when using pdf
boxcox.tex(lake.lm0)
dev.off()

#packages(alr3)
#powerTransform(pcb ~ year + length, data=laketrout)

## the meaning of interaction
##postscript(file=paste(plotDIR, "interaction.eps", sep="/"), width=4.5, height=4, horizontal=F)
tikz(file=paste(plotDIRch5, "interaction.tex", sep="/"),
     width=4.5, height=4, standAlone=F)
par(mar=c(4,4,0.25,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
pred.lm1<-predict(lake.lm4,newdata=data.frame(year=1974:2000,
                               len.c=rep(-2.6,2000-1974+1)),type="response")
pred.lm2<-predict(lake.lm4,newdata=data.frame(year=1974:2000,
                               len.c=rep(0,2000-1974+1)),type="response")
pred.lm3<-predict(lake.lm4,newdata=data.frame(year=1974:2000,
                               len.c=rep(9,2000-1974+1)),type="response")
pred.lm4<-predict(lake.lm1,newdata=data.frame(year=1974:2000),type="response")

plot(c(1974,2000),range((laketrout$pcb), na.rm=T),type="n", xlab="Year", ylab="PCB")
points(jitter(laketrout$year), laketrout$pcb, col="gray")
lines(1974:2000,exp(pred.lm1),col="green", lwd=2)
lines(1974:2000,exp(pred.lm2),col="red", lwd=2)
lines(1974:2000,exp(pred.lm3),col="blue", lwd=2)
lines(1974:2000,exp(pred.lm4),col="black", lwd=2)
dev.off()


### The Finnish Lakes example ###

summer.All <- read.table(paste(dataDIR, "summerAll.csv", sep="/"), sep=",",header=T)

#names(summer.All)
# [1] "totp"  "chla"  "type"  "lake"  "year"  "totn"  "month" "depth" "surfa"
#[10] "color"

summer.All <- summer.All[log(summer.All$chla) > -20 ,]
#> names(summer.All)
# [1] "totp"  "chla"  "type"  "lake"  "year"  "totn"  "month" "depth" "surfa"
#[10] "color"

summer.All$y <- log(summer.All$chla)
summer.All$lxp<- scale(log(summer.All$totp), scale=F)
summer.All$lxn<- scale(log(summer.All$totn), scale=F)
summer.All$type.lake <- paste(summer.All$type, summer.All$lake)
summer.All$npr <- scale(log(summer.All$totn/summer.All$totp), scale=F)

### collinearity example ####
lake1 <- summer.All[summer.All$lake==19600 & summer.All$type==2,]
lake1$lxn <- scale(log(lake1$totn), scale=F)
lake1$lxp <- scale(log(lake1$totp), scale=F)
lake1$npr <- scale(lake1$npr, scale=F)
lake2 <- summer.All[summer.All$lake==1070,]
lake2$lxn <- scale(log(lake2$totn), scale=F)
lake2$lxp <- scale(log(lake2$totp), scale=F)
lake2$npr <- scale(lake2$npr, scale=F)

## lake 2 pairs

pairs(log(chla)~log(totn)+log(totp), data=lake2)
## Fgure 5.11
##postscript(file=paste(plotDIR, "FinnData.eps", sep="/"),
##           height=4, width=4, horizontal=F)
tikz(file=paste(plotDIRch5, "FinnData.tex", sep="/"),
           height=4.5, width=4.5, standAlone=F)
par( mgp=c(1.25,.125,0), las=1, tck=0.01)
pairs(y~lxn+lxp+npr, data=lake2, cex=0.5,
      labels=c("log Chla","log TN","log TP","log N:P"))
dev.off()

postscript(file=paste(plotDIR, "FinnCorr.eps", sep="/"),
           width=3, height=3, horizontal=F)
par(mar=c(3.5,3.5,.5,.25), mgp=c(1.5,0.125,0), tck=0.01)
plot(totn~totp, data=lake2, log="xy", xlab="TP", ylab="TN", cex=0.5)
#title(main="Lake 1070", cex=0.75)
dev.off()

### coplots ###
##postscript(file=paste(plotDIR, "Finncoplot1.eps", sep="/"),
##           width=4.5, height=4, horizontal=F)
tikz(file=paste(plotDIRch5, "Finncoplot1.tex", sep="/"),
     width=4.5, height=4, standAlone=F)
par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), tck=0.02)
given.tn <- co.intervals(lake1$lxn, number=4, overlap=.1)
 coplot(y ~ lxp | lxn, data = lake1, given.v=given.tn, rows=1,
        xlab=c("log TP", "Given: Centered log TN"), ylab="log Chla",
        panel=panel.smooth, col="gray", col.smooth=1)
dev.off()

##postscript(file=paste(plotDIR, "Finncoplot2.eps", sep="/"),
##           width=4.5, height=4, horizontal=F)
tikz(file=paste(plotDIRch5, "Finncoplot2.tex", sep="/"),
     width=4.5, height=4, standAlone=F)
par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), tck=0.02)
given.tp <- co.intervals(lake1$lxp, number=4, overlap=.1)
coplot(y ~ lxn | lxp, data = lake1, given.v=given.tp, rows=1,
       xlab=c("log TN", "Given: Centered log TP"), ylab="log Chla",
       panel=panel.smooth, col="gray", col.smooth=1)
dev.off()

## Figure 5.12
##postscript(file=paste(plotDIR, "Finncoplot3.eps", sep="/"),
##           width=4.5, height=4, horizontal=F)
tikz(file=paste(plotDIRch5, "Finncoplot3.tex", sep="/"),
     width=4.5, height=4, standAlone=F)
par(mar=c(3.5,3.5,3.5,.25), mgp=c(1.5,0.5,0), tck=0.02)
given.tn <- co.intervals(lake2$lxn, number=4, overlap=.1)
coplot(y ~ lxp | lxn, data = lake2, given.v=given.tn, rows=1,
       xlab=c("log TP", "Given: Centered log TN"), ylab="log Chla",
       panel=panel.smooth, col="gray", col.smooth=1)
dev.off()

## Figure 5.13
##postscript(file=paste(plotDIR, "Finncoplot4.eps", sep="/"),
##           width=4.5, height=4, horizontal=F)
tikz(file=paste(plotDIRch5, "Finncoplot4.tex", sep="/"),
     width=4.5, height=4, standAlone=F)
par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), tck=0.02)
given.tp <- co.intervals(lake2$lxp, number=4, overlap=.1)
coplot(y ~ lxn | lxp, data = lake2, given.v=given.tp, rows=1,
       xlab=c("log TN", "Given: Centered log TP"), ylab="log Chla",
       panel=panel.smooth, col="gray", col.smooth=1)
dev.off()

two.lakes <- rbind(lake1, lake2)
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.75 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
trellis.par.set(trellis.par.temp)
obj1 <- xyplot(y~lxp|factor(lake), data=two.lakes, panel=function(x, y, ...){
    panel.xyplot(x,y, col="gray",...)
    panel.lmline(x, y,...)
    panel.grid()
}, aspect=1, ylab="log Chla", xlab="log TP")
obj2 <- xyplot(y~lxn|factor(lake), data=two.lakes, panel=function(x, y, ...){
    panel.xyplot(x,y, col="gray", cex=0.5, ...)
    panel.lmline(x, y, cex=0.5, ...)
    panel.grid()
}, aspect=1, ylab="log Chla", xlab="log TN")

##postscript(file=paste(plotDIR, "finnlake11.eps", sep="/"),
##           height=5, width=5, horizontal = F)
tikz(file=paste(plotDIRch5, "finnlake11.tex", sep="/"),
     height=5, width=5, standAlone = F)
print(obj1, position=c(0,0,1,0.55), more=T)
print(obj2, position=c(0,0.45,1,1), more=F)
dev.off()

lake2PN <- data.frame(tptn = c(lake2$lxp, lake2$lxn),
chla=rep(lake2$y, 2),
nutrients=rep(c("TP","TN"), each=dim(lake2)[1]))

##postscript(file=paste(plotDIR, "finnlake1.eps", sep="/"),
##           height=2.5, width=4.5, horizontal=F)

tikz(file=paste(plotDIRch5, "finnlake1.tex", sep="/"),
     height=2.5, width=4.5, standAlone=F)
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.75 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
trellis.par.set(trellis.par.temp)
xyplot(chla ~ tptn|nutrients, data=lake2PN, panel=function(x,y,...){
    panel.xyplot(x, y, col="gray", ...)
    panel.lmline(x, y, ...)
    panel.grid()
}, aspect=1,  ylab="log Chla", xlab="log nutrient concentration",
       scales=list(x=list(relation="free")))
dev.off()

Finn.lm1 <- lm(y ~ lxp + lxn, data=lake1)
display(Finn.lm1)

Finn.lm2 <- lm(y ~ lxp + lxn, data=lake2)
display(Finn.lm2)

Finn.lm3 <- lm(y ~ lxp * lxn, data=lake1)
display(Finn.lm3)

Finn.lm4 <- lm(y ~ lxp * lxn, data=lake2)
display(Finn.lm4)

coef.lm3 <- coef(Finn.lm4)
quant.n <- quantile(lake2$lxn, prob=c(0.025, 0.25, 0.5,0.75, 0.975))
quant.p <- quantile(lake2$lxp, prob=c(0.025, 0.25, 0.5,0.75, 0.975))

## Figure 5.14
##postscript(file=paste(plotDIR, "finnInt.eps", sep="/"),
##           height= 2.75, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch5, "finnInt.tex", sep="/"),
           height= 2.75, width=4.5, standAlone=F)
par(mfrow=c(1,2), mgp=c(1.25, 0.125,0), mar=c(3,3,1, 0.25), las=1, tck=0.02)
plot(y~lxp, data=lake2, axes=F, xlab="TP", ylab="Chla", cex=0.5)
axis(1, at=log(c(1, 2, 5, 10,20,50,100))-attr(lake2$lxp, which="scaled:center"),
     labels=c(1, 2, 5, 10,20,50,100))
axis(2, at=log(c(1,2,5,10,20,50)),
     labels=c(1,2,5,10,20,50))
box()
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[1] +
      coef.lm3[4]*x*quant.n[1], add=T, lty=1)
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[2] +
      coef.lm3[4]*x*quant.n[2], add=T, lty=2)
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[3] +
      coef.lm3[4]*x*quant.n[3], add=T, lty=3)
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[4] +
      coef.lm3[4]*x*quant.n[4], add=T, lty=4)
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[5] +
      coef.lm3[4]*x*quant.n[5], add=T, lty=5)
legend(x=log(30)-attr(lake2$lxp, which="scaled:center"),
       y=log(5),
       legend=c("2.5\\%","25\\%","50\\%","75\\%","97.5\\%"),
       cex=0.5, lty=1:5, bty="n")
plot(y~lxn, data=lake2, axes=F, xlab="TN", ylab="Chla", cex=0.5)
axis(1, at=log(c(500,1000))-attr(lake2$lxn, which="scaled:center"),
     labels=c(500,1000))
axis(2, at=log(c(1,2,5,10,20,50)),
     labels=c(1,2,5,10,20,50))
box()
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[1] +
      coef.lm3[4]*x*quant.p[1], add=T, lty=1)
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[2] +
      coef.lm3[4]*x*quant.p[2], add=T, lty=2)
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[3] +
      coef.lm3[4]*x*quant.p[3], add=T, lty=3)
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[4] +
      coef.lm3[4]*x*quant.p[4], add=T, lty=4)
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[5] +
      coef.lm3[4]*x*quant.p[5], add=T, lty=5)
dev.off()

lake2 $ npr <- log(lake2$totn/lake2$totp)
Finn.lm5 <- lm(y ~ lxp + lxn + npr, data=lake1)
display(Finn.lm5)


## positive interaction,

lake3 <- summer.All[summer.All$lake==14700,]
lake3$lxn <- scale(log(lake3$totn), scale=F)
lake3$lxp <- scale(log(lake3$totp), scale=F)
lake4 <- summer.All[summer.All$lake==21800,]
lake4$lxn <- scale(log(lake4$totn), scale=F)
lake4$lxp <- scale(log(lake4$totp), scale=F)

Finn.lm5 <- lm(y ~ lxp * lxn, data=lake3)
display(Finn.lm5)
##### conditional plot

## Figure 5.15
##postscript(file=paste(plotDIR, "Finncoplotlake31.eps", sep="/"),
##           width=4.5, height=4, horizontal=F)
tikz(file=paste(plotDIRch5, "Finncoplotlake31.tex", sep="/"),
           width=4.5, height=4, standAlone=F)

par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), las=1, tck=0.02)
given.tn <- co.intervals(lake3$lxn, number=4, overlap=.1)
 coplot(y ~ lxp | lxn, data = lake3, given.v=given.tn, rows=1,
        xlab=c("log TP", "Given: Centered log TN"), ylab="log Chla",
        panel=panel.smooth, col="gray", col.smooth=1)
dev.off()

## Figure 5.16
##postscript(file=paste(plotDIR, "Finncoplotlake32.eps", sep="/"),
##           width=4.5, height=4, horizontal=F)
tikz(file=paste(plotDIRch5, "Finncoplotlake32.tex", sep="/"),
           width=4.5, height=4, standAlone=F)
par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), tck=0.02)
given.tp <- co.intervals(lake3$lxp, number=4, overlap=.1)
 coplot(y ~ lxn | lxp, data = lake3, given.v=given.tp, rows=1,
        xlab=c("log TN", "Given: Centered log TP"), ylab="log Chla",
        panel=panel.smooth, col="gray", col.smooth=1)
dev.off()
##### interaction plot

coef.lm3 <- coef(Finn.lm5)
quant.n <- quantile(lake3$lxn, prob=c(0.025, 0.25, 0.5,0.75, 0.975))
quant.p <- quantile(lake3$lxp, prob=c(0.025, 0.25, 0.5,0.75, 0.975))

## Figure 5.17

##postscript(file=paste(plotDIR, "finnIntlake3.eps", sep="/"),
##           height= 2.75, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch5, "finnIntlake3.tex", sep="/"),
     height= 2.75, width=4.5, standAlone=F)
par(mfrow=c(1,2), mgp=c(1.25, 0.125,0), mar=c(3,3,1, 0.25), las=1, tck=0.02)
plot(y~lxp, data=lake3, axes=F, xlab="TP", ylab="Chla", cex=0.5)
axis(1, at=log(c(1, 2, 5, 10,20,50,100,200))-attr(lake3$lxp, which="scaled:center"),
     labels=c(1, 2, 5, 10,20,50,100,200))
axis(2, at=log(c(1,2,5,10,20,50,100)),
     labels=c(1,2,5,10,20,50,100))
box()
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[1] +
      coef.lm3[4]*x*quant.n[1], add=T, lty=1)
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[2] +
      coef.lm3[4]*x*quant.n[2], add=T, lty=2)
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[3] +
      coef.lm3[4]*x*quant.n[3], add=T, lty=3)
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[4] +
      coef.lm3[4]*x*quant.n[4], add=T, lty=4)
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[5] +
      coef.lm3[4]*x*quant.n[5], add=T, lty=5)
legend(x=log(20)-attr(lake3$lxp, which="scaled:center"),
       y=log(4),
       legend=c("2.5\\%","25\\%","50\\%","75\\%","97.5\\%"),
       cex=0.5, lty=1:5, bty="n")
plot(y~lxn, data=lake3, axes=F, xlab="TN", ylab="Chla", cex=0.5)
axis(1, at=log(c(200,500,1000,1500))-attr(lake3$lxn, which="scaled:center"),
     labels=c(200,500,1000,1500))
axis(2, at=log(c(1,2,5,10,20,50,100)),
     labels=c(1,2,5,10,20,50,100))
box()
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[1] +
      coef.lm3[4]*x*quant.p[1], add=T, lty=1)
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[2] +
      coef.lm3[4]*x*quant.p[2], add=T, lty=2)
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[3] +
      coef.lm3[4]*x*quant.p[3], add=T, lty=3)
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[4] +
      coef.lm3[4]*x*quant.p[4], add=T, lty=4)
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[5] +
      coef.lm3[4]*x*quant.p[5], add=T, lty=5)
dev.off()



## negative interaction
Finn.lm6 <- lm(y ~ lxp * lxn, data=lake4)
display(Finn.lm6)

##### conditional plot
##postscript(file=paste(plotDIR, "Finncoplotlake41.eps", sep="/"),
##           width=4.5, height=4, horizontal=F)
tikz(file=paste(plotDIRch5, "Finncoplotlake41.tex", sep="/"),
           width=4.5, height=4, standAlone=F)
par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), tck=0.02)
given.tn <- co.intervals(lake4$lxn, number=4, overlap=.1)
 coplot(y ~ lxp | lxn, data = lake4, given.v=given.tn, rows=1,
        xlab=c("log TP", "Given: Centered log TN"), ylab="log Chla",
        panel=panel.smooth, col="gray", col.smooth=1)
dev.off()

##postscript(file=paste(plotDIR, "Finncoplotlake42.eps", sep="/"),
##           width=4.5, height=4, horizontal=F)
tikz(file=paste(plotDIRch5, "Finncoplotlake42.tex", sep="/"),
           width=4.5, height=4, standAlone=F)
par(mar=c(3,3,3,.25), mgp=c(1.5,0.5,0), tck=0.02)
given.tp <- co.intervals(lake4$lxp, number=4, overlap=.1)
 coplot(y ~ lxn | lxp, data = lake4, given.v=given.tn, rows=1,
        xlab=c("log TN", "Given: Centered log TP"), ylab="log Chla",
        panel=panel.smooth, col="gray", col.smooth=1)
dev.off()
##### interaction plot
coef.lm3 <- coef(Finn.lm6)
quant.n <- quantile(lake4$lxn, prob=c(0.025, 0.25, 0.5,0.75, 0.975))
quant.p <- quantile(lake4$lxp, prob=c(0.025, 0.25, 0.5,0.75, 0.975))

## Figure 5.18
##postscript(file=paste(plotDIR, "finnIntlake4.eps", sep="/"),
##           height= 2.75, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch5, "finnIntlake4.tex", sep="/"),
           height= 2.75, width=4.5, standAlone=F)
par(mfrow=c(1,2), mgp=c(1.25, 0.125,0), mar=c(3,3,1, 0.25), las=1, tck=0.02)
plot(y~lxp, data=lake4, axes=F, xlab="TP", ylab="Chla", cex=0.5)
axis(1, at=log(c(1, 2, 5, 10,20,50,100,200))-attr(lake4$lxp, which="scaled:center"),
     labels=c(1, 2, 5, 10,20,50,100,200))
axis(2, at=log(c(1,2,5,10,20,50,100)),
     labels=c(1,2,5,10,20,50,100))
box()
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[1] +
      coef.lm3[4]*x*quant.n[1], add=T, lty=1)
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[2] +
      coef.lm3[4]*x*quant.n[2], add=T, lty=2)
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[3] +
      coef.lm3[4]*x*quant.n[3], add=T, lty=3)
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[4] +
      coef.lm3[4]*x*quant.n[4], add=T, lty=4)
curve(coef.lm3[1]+coef.lm3[2]*x+coef.lm3[3]*quant.n[5] +
      coef.lm3[4]*x*quant.n[5], add=T, lty=5)
legend(x=log(100)-attr(lake4$lxp, which="scaled:center"),
       y=log(10),
       legend=c("2.5\\%","25\\%","50\\%","75\\%","97.5\\%"),
       cex=0.5, lty=1:5, bty="n")
plot(y~lxn, data=lake4, axes=F, xlab="TN", ylab="Chla", cex=0.5)
axis(1, at=log(c(200,500,1000,1500))-attr(lake4$lxn, which="scaled:center"),
     labels=c(200,500,1000,1500))
axis(2, at=log(c(1,2,5,10,20,50,100)),
     labels=c(1,2,5,10,20,50,100))
box()
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[1] +
      coef.lm3[4]*x*quant.p[1], add=T, lty=1)
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[2] +
      coef.lm3[4]*x*quant.p[2], add=T, lty=2)
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[3] +
      coef.lm3[4]*x*quant.p[3], add=T, lty=3)
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[4] +
      coef.lm3[4]*x*quant.p[4], add=T, lty=4)
curve(coef.lm3[1]+coef.lm3[3]*x+coef.lm3[2]*quant.p[5] +
      coef.lm3[4]*x*quant.p[5], add=T, lty=5)
dev.off()


### Uncertainty ###

lm1.pred <- predict(lake.lm1, new=data.frame(year=2007), se.fit=T)

### ELISA example

mc <- c(0.167, 0.444, 01.110, 2.220, 5.550)
rOD <- c(0.784, 0.588, 0.373, 0.270,0.202)
stdcrv <- lm(log(mc) ~ rOD)
betas <- coef(stdcrv)
aug01 <- predict(stdcrv, newdata=data.frame(rOD=0.261), se.fit=T, interval="prediction")      
tikz(paste(plotDIRch5, "elisalm.tex", sep="/"),
     width=4, height=3, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(mc ~ rOD, xlab="$rOD$",
     ylab="Microcystin Concentration ($\\mu$g/L)",
     ylim=range(mc, exp(aug01$fit)))
curve(exp(betas[1] + betas[2]*x), add=T, lty=2)
segments(x0=0.261, x1=0.261,
         y0=exp(aug01$fit[2]), y1=exp(aug01$fit[3]))
segments(x0=0.261, x1=0.261,
         y0=exp(aug01$fit[1]-qt(0.975, 3)*aug01$se.fit),
         y1=exp(aug01$fit[1]+qt(0.975,3)*aug01$se.fit), lwd=3)
points(0.261, exp(aug01$fit[1]), pch=16)
abline(h=1, col="gray")
dev.off()

## simulations
n.sims<-5000
sim.results <- sim(lake.lm1, n.sims=1000)
predict.PCB07 <- exp(sim.results@coef[,1] +
                     sim.results@coef[,2]*(2007-1974) +
                     0.5* sim.results@sigma)

predict.PCB00 <- exp(sim.results@coef[,1] +
                     sim.results@coef[,2]*(2000-1974) +
                     0.5 * sim.results@sigma)
percentages <- 1-predict.PCB07/predict.PCB00
postscript(file="pcbreduction.eps", width=3, height=2.5, horizontal=F)
par(mar=c(4,3,1,0.5), mgp=c(1.5,0.5,0))
hist(percentages*100, xlab="Predicted % PCB reduction", ylab="", prob=T, main="")
dev.off()


sim.2 <- sim(lake.lm7, 1000)
predict2.PCB07 <- exp(sim.2$beta[,1] +
                      sim.2$beta[,3]*(2007-1974) +
                      0.5* sim.2$sigma)

predict2.PCB00 <- exp(sim.2$beta[,1] +
                      sim.2$beta[,3]*(2000-1974) +
                      0.5 * sim.2$sigma)
percentages2 <- 1-predict2.PCB07/predict2.PCB00
postscript(file="pcbreduction2.eps", width=3, height=2.5, horizontal=F)
par(mar=c(4,3,1,0.5), mgp=c(1.5,0.5,0))
hist(percentages2*100, xlab="Predicted % PCB reduction", ylab="", prob=T, main="")
dev.off()

#### Two Way ANOVA ####
attach(mangrove.sponge)
mangrove.sponge$Control <- as.numeric(Treatment=="Control")
mangrove.sponge$Foam <- as.numeric(Treatment=="Foam")
mangrove.sponge$PurpleS <- as.numeric(Treatment=="Haliclona")
mangrove.sponge$RedS <- as.numeric(Treatment=="Tedania")
detach()
mangrove.lmDM <- lm(RootGrowthRate ~ Foam+PurpleS+RedS, data=mangrove.sponge)
display(mangrove.lmDM, 4)

attach(mangrove.sponge)
mangrove.sponge$bbs <- as.numeric(Location=="bbs")
mangrove.sponge$etb <- as.numeric(Location=="etb")
mangrove.sponge$lcn <- as.numeric(Location=="lcn")
mangrove.sponge$lcs <- as.numeric(Location=="lcs")
detach()
mangrove.lmDM2 <- lm(RootGrowthRate ~ Foam+PurpleS+RedS +
                                      etb+lcn+lcs , data=mangrove.sponge)
display(mangrove.lmDM2, 4)

mangrove.lm2w <- lm(RootGrowthRate ~ Treatment+Location, data=mangrove.sponge)


#### cross-validation ####


B <- 500
pred.resid <- numeric()
pred.bias <- numeric()
for (i in 1:B){
    smpl <- unique(sample(dim(laketrout)[1], dim(laketrout)[1], replace=T))
    asmpl <- seq(1,dim(laketrout)[1])[-smpl]
    lm.temp <- lm(log(pcb) ~ I(year-1974) + len.c * factor(size), data=laketrout[smpl,])
    pred.temp <- predict(lm.temp, new=laketrout[asmpl,], type="response")
    pred.resid [i] <- sd(log(laketrout[asmpl,"pcb"]) - pred.temp)
    pred.bias[i] <- mean(log(laketrout[asmpl,"pcb"]) - pred.temp)
}

postscript("c:/users/song/teaching/env210/fall2007/notes/PCBxvalid.eps", width=4.5, height=2.5, horizontal=F)
par(mfrow=c(1,2), mar=c(3,1.5,0.25,0.25), mgp=c(1.5,.5,0))
hist(pred.bias, main="", xlab="Prediction Bias")
abline(v=0, col="red", lwd=3)
hist(pred.resid, main="", xlab="Prediction Standard Deviation")
abline(v=sd(lake.lm9$residuals), col="red", lwd=3)
dev.off()
