source("FrontMatter.R")

##################
### Chapter 3  ###
##################
plotDIRch3 <- paste (plotDIR, "chapter3", "figures", sep="/")

##### Everglades reference sites
wca2tp <- read.csv(paste(dataDIR, "WCA2TP.csv", sep="/"), header=T)
wca2tp$TP <- 1000*wca2tp$RESULT
TP.reference <- wca2tp[wca2tp$Type=="R",]
TP.reference$SITE <- as.vector( TP.reference$SITE)
TP.reference$Month <- ordered(months(as.Date(TP.reference$Date, "%m/%d/%y"), T),
  levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

## figure \ref{fig:stdnorm} (3.1)
#postscript(file=paste(plotDIR, "normal1.eps", sep="/"),
#           height=1.5, width=2.5, horizontal=F)
tikz(file=paste(plotDIRch3, "normal1.tex", sep="/"),
     height=1.5, width=2.5, standAlone=F)
par(mar=c(2.5,2.5,.5,1), mgp=c(1.25, 0.125, 0), tck=0.01, las=1)
plot(seq(-3,3,,100), dnorm(seq(-3,3,,100)), type="l",
     xlab="$Y$", ylab="", axes=F, cex=0.75)
#polygon(c(-1, seq(-1, 1.3,,50), 1.3), c(0, dnorm(seq(-1, 1.3,,50)), 0),
#col=gray(.65))
polygon(c(-3, seq(-3, qnorm(0.25),,50), qnorm(0.25)),
        c(0, dnorm(seq(-3, qnorm(0.25),,50)), 0), col=gray(.45))
polygon(c(2, seq(2, 3,,50), 3), c(0, dnorm(seq(2, 3,,50)), 0),
        col=gray(.8))
text(-2, 0.2, "0.25", cex=0.75)
axis(1, at=c(-3,-2, qnorm(0.25),0, 1,2,3), cex=0.75,
     label=c("-3", "-2","$y$", "0", "1", "2", "3"))
axis(1, at=-1, label="")
axis(2, cex=0.75)
box()
dev.off()

dnorm(0, 0, 1)
pnorm(2, 0, 1)
1-pnorm(2,0,1)
qnorm(0.25, 0, 1)
qnorm(0.025, 0, 1)
qnorm(0.95, 0, 1)

## Figure 3.2
#postscript(file=paste(plotDIR,"referCHANCE.eps", sep="/"),
#           height=2, width=3.5, horizontal=F)
tikz(file=paste(plotDIRch3, "referCHANCE.tex", sep="/"),
     height=2, width=3.5, standAlone=F)
par(mar=c(2.5,2.5,.5,1), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
hist(log(TP.reference$TP), axes=F, prob=T, dens=-1,col=grey(0.65), border=grey(0.85),
     xlab="TP (ppb)", ylab="Density", main="")
lines(seq(0, log(100),,100),
      dnorm(seq(0, log(100),,100),
            mean(log(TP.reference$TP)),
            sd(log(TP.reference$TP))))
axis(1, at=log(c(1, 5, 10, 20, 50)), label=c("1", "5","10","20","50"))
axis(2)
box()
dev.off()

c(mean(log(TP.reference$TP)), sd(log(TP.reference$TP)))

quantile((TP.reference$TP),
         prob=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9,0.92, 0.95))
qlnorm(c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9,0.92, 0.95),
       mean(log(TP.reference$TP)), sd(log(TP.reference$TP)))
1-plnorm(15, mean(log(TP.reference$TP)), sd(log(TP.reference$TP)))
mean((TP.reference$TP)>15)

y <- rnorm(100)
n <- length(y)
yq <- ((1:n) - 0.5)/n
zq <- qnorm(yq, mean=0, sd=1)
plot(zq, sort(y), xlab="Standard Normal Quantile", ylab="Data")
abline(mean(y), sd(y))

my.qqplot <- function (y=rnorm(100), Ylab="Data"){
  n <- length(y)
  yq <- ((1:n) - 0.5)/n
  zq <- qnorm(yq, mean=0, sd=1)
  plot(zq, sort(y), xlab="Standard Normal Quantile", ylab=Ylab)
  abline(mean(y), sd(y))
  invisible()
}

my.qqplot()

tikz(file=paste(plotDIRch3, "evergnormQQ.tex", sep="/"),
     height=3, width=3, standAlone=F)
par(mgp=c(1.25, 0.125, 0), mar=c(3,3,3,1), las=1, tck=0.01)
qqnorm(log(TP.reference$TP))
qqline(log(TP.reference$TP))
dev.off()

qqmath(~log(TP), data=TP.reference, aspect=1,
       panel=function(x, ...){
           panel.qqmathline(x, ...)
           panel.qqmath(x, ...)
           panel.grid()
       }
      )

qqmath(~log(TP)|SITE, data=TP.reference, #aspect="xy",
       panel=function(x, ...){
           panel.qqmathline(x, ...)
           panel.qqmath(x, ...)
       })

## the Lake Erie data 
TPcomplete<- read.table(paste(dataDIR, "TPcomplete.csv", sep="/"),
                        header=T, sep=",")
TPcomplete <- TPcomplete[TPcomplete$TP > 0.01,]
tikz(file=paste(plotDIRch3, "ErieTP_Dist.tex", sep="/"),
     height=2.5, width=4, standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(log(TP)~log(DISTANCE), data=TPcomplete,
     xlab="Log distance to Maumee",
     ylab="Log TP") 
dev.off()


## Figure \ref{fig:evergBoxSL} 3.3
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
subset <- TP.reference$SITE=="E5" | TP.reference$SITE=="F5"
Evg.box <- bwplot((TP) ~ SITE, data=TP.reference[subset,],
                  xlab="Reference Sites", ylab="TP (ppb)")
tpref.median <- oneway(TP ~ SITE, data=TP.reference[subset,],
                       location=mean, spread=1)
par(mar=c(2.5,2.5,.5,1), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
Evg.Sl  <-
  xyplot(sqrt(abs(residuals(tpref.median)))~
           jitter(fitted.values(tpref.median), factor=.125),
#                                             aspect=2,
         panel=function(x,y){
           panel.xyplot(x,y, col=gray(0.4), cex=0.5)
           srmads <- sqrt(tapply(abs(residuals(tpref.median)),
                                 TP.reference[subset,"SITE"], 
                                 median))
           oo <- order(tpref.median$location)
           panel.lines(tpref.median$location[oo],srmads[oo])
         },
#        sub = list("Everglades Data",cex=.8),
         xlab="Mean TP",
         ylab="Sqrt Abs Residuals")

#postscript(file=paste(plotDIR, "evergBoxSL.eps", sep="/"),
#           height=2, width=5, horizontal=F)
tikz(file=paste(plotDIRch3, "evergBoxSL.tex", sep="/"),
     height=3, width=5.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
print(Evg.box, position = c(0, 0, 0.5, 1.0), more = T)
print(Evg.Sl , position = c(0.5, 0, 1, 1.0), more = F)
dev.off()


#postscript(file=paste(plotDIR,"evergSLlog.eps", sep="/"),
#           height=3, width=3, horizontal=F)
tikz(file=paste(plotDIRch3, "evergSLlog.tex", sep="/"),
         height=3, width=3, standAlone=F)
tpref.median <- oneway(log(TP) ~ SITE, data=TP.reference,
                       location=mean, spread=1)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
par(mar=c(2.5,2.5,.5,1), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
xyplot(sqrt(abs(residuals(tpref.median)))~jitter(fitted.values(tpref.median),
                                                 factor=.5),
       aspect=0.75,
       panel=function(x,y){
         panel.xyplot(x,y)
         srmads <- sqrt(tapply(abs(residuals(tpref.median)),
                               TP.reference$SITE, median))
         oo <- order(tpref.median$location)
         panel.lines(tpref.median$location[oo],srmads[oo])
       },
#      sub = list("Everglades Data",cex=.8),
       xlab="Jittered Mean log TP",
       ylab="Square Root Absolute Residual")
dev.off()

### Figure \ref{fig:evergHists} 3.4

#postscript(file=paste(plotDIR,"evergHist.eps", sep="/"),
#           height=2, width=5, horizontal=F)
tikz(file=paste(plotDIRch3, "evergHist.tex", sep="/"),
     height=2.5, width=5.5, standAlone=F)
par(mfrow=c(1,2), mar=c(2.5,2.5,.5,1), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
hist(log(TP.reference$TP), dens=-1, xlab="TP (ppb)",
     ylab="", main="", nclass=20)
hist(log(TP.reference$TP), dens=-1, xlab="TP (ppb)",
     ylab="", main="", nclass=10)
dev.off()

### Figure \ref{fig:qplot1} 3.5
TP.example<- c(0.21, 0.35, 0.50, 0.64, 0.79,0.90, 1.00, 1.01, 1.12, 5.66)
#postscript(file=paste(plotDIR,"evergQplot.eps", sep="/"),
#           height=2, width=3, horizontal=F)
tikz(file=paste(plotDIRch3,"evergQplot.tex", sep="/"),
           height=2, width=3, standAlone=F)
par(mar=c(2.5,2.5,.5,1), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
plot((1:length(TP.example) -0.5)/length(TP.example), TP.example, xlab="$f$", ylab="TP (ppb)")
dev.off()

## Figure \ref{fig:clevelandBox} 3.6
data <-
  c(0.9, 1.6, 2.26305, 2.55052, 2.61059, 2.69284, 2.78511, 2.80955,
    2.94647, 2.96043, 3.05728, 3.15748, 3.18033, 3.20021,
    3.20156, 3.24435, 3.33231, 3.34176, 3.3762, 3.39578, 3.4925,
    3.55195, 3.56207, 3.65149, 3.72746, 3.73338, 3.73869,
    3.80469, 3.85224, 3.91386, 3.93034, 4.02351, 4.03947,
    4.05481, 4.10111, 4.26249, 4.28782, 4.37586, 4.48811,
    4.6001, 4.65677, 4.66167, 4.73211, 4.80803, 4.9812, 5.17246,
    5.3156, 5.35086, 5.36848, 5.48167, 5.68, 5.98848, 6.2, 7.1,
    7.4)

## explaining the box plot
uq <- quantile(data,.75)
lq <- quantile(data,.25)
r <- 1.5*(uq-lq)
h <- c(lq-r,1.6,lq,uq,6.2,uq+r)
writing <- c("lower quartile - 1.5 r",
             "lower adjacent value",
             "lower quartile",
             "upper quartile",
             "upper adjacent value",
             "upper quartile + 1.5 r")

n <- length(data)

#postscript(file=paste(plotDIR, "explainBox.eps", sep="/"),
#           width=4.75, height=2.25, horizontal=F)
tikz(file=paste(plotDIRch3, "explainBox.tex", sep="/"),
           width=4.75, height=2.25, standAlone=F)
par(mfrow=c(1,2), mar=c(2.5,2.5,.25,0.125), mgp=c(1.25, 0.125, 0),
    las=1, tck=0.01)
            #layout(matrix(c(1,2), nrow=1), width=c(1,1.5))
boxplot(data, rep(NA, length(data)), rep(NA, length(data)), ylab = "Data")
usr <- par("usr")
x <- usr[1] + (usr[2] - usr[1]) * 1/3
at <- c(0.9, 1.6, 3.2, 3.8, 4.65, 6.2, 7.2)
arrows(rep(x * 1.15, 7), at, rep(x, 7), at, length=0.075)
mtext("Main Title",1,1,cex=.8)
text(rep(x * 1.2, 7), at, adj = 0, cex=0.7,
     labels = c("outside value", "lower adjacent value",
       "lower quartile", "median", "upper quartile",
       "upper adjacent value", "outside values"))

plot(((1:n)-0.5)/n, data, xlab="f-value", ylab="Data")
abline(h = h, lwd =2, col = gray(0.85))
text(rep(0,3), h[4:6], writing[4:6], adj=0, cex=0.75)
text(rep(1,3), h[1:3], writing[1:3], adj=1, cex=0.75)
points(((1:n)-0.5)/n, data)

dev.off()

## Figure \ref{fig:addVmlt} 3.7
### Q-Q plot for additive and multiplicative shifts
x.data <- rnorm(500)
y.data <- rnorm(500, 2)
#postscript(file=paste(plotDIR,"additive.eps", sep="/"),
#           height=4, width=3.5, horizontal=F)
tikz(file=paste(plotDIRch3,"additive.tex", sep="/"),
     height=4, width=3.5, standAlone=F)
layout(matrix(c(1,3,2,3), nrow=2), heights=c(1,2), respect=T)
par(mar=c(2.5,2.5,.5,0.5), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
hist(x.data, xlim=range(c(x.data,y.data)), xlab="x", ylab="", main="", prob=T)
curve(dnorm(x), add=T, from=-3, to=3)
hist(y.data, xlim=range(c(x.data,y.data)), xlab="y", ylab="", main="", prob=T)
curve(dnorm(x, mean=2), add=T, from=-1, to=5)
qqplot(x.data, y.data, xlim=range(c(x.data,y.data)),
       ylim=range(c(x.data,y.data)), xlab="x", ylab="y")
abline(0,1, col=gray(0.5))
abline(2,1, col=gray(0.5))
dev.off()

x.data2 <- exp(rnorm(100))
y.data2 <- exp(rnorm(100, 1))

#postscript(file=paste(plotDIR, "multiplicative.eps", sep="/"),
#           height=4, width=3.5, horizontal=F)
tikz(file=paste(plotDIRch3, "multiplicative.tex", sep="/"),
     height=4, width=3.5, standAlone=F)
layout(matrix(c(1,3,2,3), nrow=2), heights=c(1,2), respect=T)
par(mar=c(2.5,2.5,.5,0.5), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
hist(x.data2, xlim=range(c(x.data2,y.data2)),
     xlab="x", ylab="", main="", prob=T, nclass=10)
curve(dlnorm(x), add=T, from=exp(-3), to=exp(3))
hist(y.data2, xlim=range(c(x.data2,y.data2)),
     xlab="y", ylab="", main="", prob=T, nclass=10)
curve(dlnorm(x, mean=1), add=T, from=exp(-2), to=exp(4))
qqplot(log(x.data2), log(y.data2),  xlab="x", ylab="y")
abline(0,1, col=gray(0.5))
dev.off()

### Figure 3.7
#postscript(file=paste(plotDIR, "addVmult.eps", sep="/"),
#           height=3.75, width=4.75, horizontal=F)
tikz(file=paste(plotDIRch3, "addVmult.tex", sep="/"),
     height=3.75, width=4.75, standAlone=F)
layout(rbind(c(1,2,0,4,5), c(3,3,0,6,6)), widths=c(1,1,lcm(0.5), 1,1),
       heights=c(1,2), respect=T)
par(mar=c(2.5,2.5,.5,0.5), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)

hist(x.data, xlim=range(c(x.data,y.data)), xlab="x", ylab="", main="", prob=T)
curve(dnorm(x), add=T, from=-3, to=3)
hist(y.data, xlim=range(c(x.data,y.data)), xlab="y", ylab="", main="", prob=T)
curve(dnorm(x, mean=2), add=T, from=-1, to=5)
qqplot(x.data, y.data, xlim=range(c(x.data,y.data)),
       ylim=range(c(x.data,y.data)),xlab="x", ylab="y")
abline(0,1, col=gray(0.5))
abline(2,1, col=gray(0.5))

hist(x.data2, xlim=range(c(x.data2,y.data2)),
     xlab="x", ylab="", main="", prob=T, nclass=10)
curve(dlnorm(x), add=T, from=exp(-3), to=exp(3))
hist(y.data2, xlim=range(c(x.data2,y.data2)),
     xlab="y", ylab="", main="", prob=T, nclass=10)
curve(dlnorm(x, mean=1), add=T, from=exp(-2), to=exp(4))
qqplot(x.data2, y.data2, xlim=range(c(x.data2,y.data2)),
       ylim=range(c(x.data2,y.data2)), xlab="x", ylab="y")
abline(0,1, col=gray(0.5))
dev.off()

## Figure \ref{fig:cars.test} 3.8
## Car.test.frame in package rpart
packages(rpart)
names(car.test.frame)
#postscript(file=paste(plotDIR, "carsTest.eps", sep="/"),
#           width=4.75, height=2, horizontal=F)
tikz(file=paste(plotDIRch3, "carsTest.tex", sep="/"),
           width=4.75, height=2, standAlone=F)
par(mfrow=c(1,2), mar=c(2.5,2.5,.5,0.5),
    mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
plot(Mileage ~ Weight, data=car.test.frame,
     xlab="Weight (lb)", ylab="Mileage (mpg)")
abline(lsfit(car.test.frame$Weight, car.test.frame$Mileage))

plot(1/Mileage ~ Weight, data=car.test.frame,
     xlab="Weight (lb)", ylab="Mileage (mpg)")
abline(lsfit(car.test.frame$Weight, car.test.frame$Mileage),
       col=gray(0.5), lty=2)
cars.lo <- loess(I(1/Mileage) ~ Weight, data = car.test.frame)
oo <- order(cars.lo$x)
lines(cars.lo$x[oo], cars.lo$fitted[oo])
dev.off()

## Figure \ref{fig:airquality} 3.9
## Car.test.frame in package rpart

panel.hist = function(x, col.hist="grey", ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1],y, col=col.hist)
  }

#postscript(file=paste(plotDIR, "airquality.eps", sep="/"),
#           height=4.5, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch3, "airquality.tex", sep="/"),
           height=4.5, width=4.5, standAlone=F)
pairs(Ozone~Solar.R+Wind+Temp, data=airquality,
      panel=function(x, y, ...){
        points(x, y,  ...)
        panel.smooth(x, y, col.smooth=1, ...)
      }, col=gray(0.5), cex=0.5, diag.panel=panel.hist)
dev.off()

### Figure \ref{fig:irisPairs} 3.10
### the iris data
Snames <- dimnames(iris3)[[3]]
iris.df <- rbind(iris3[,,1],iris3[,,2],iris3[,,3])
iris.df <- as.data.frame(iris.df)
iris.df$Species <- factor(rep(Snames,rep(50,3)))
#pairs(iris.df[1:4],main = "Anderson?s Iris Data",
#pch = c("<", "+", ">")[unclass(iris$Species)],
#col = c(gray(0.25),gray(0.5),gray(0.25))[unclass(iris$Species)])

#postscript(file=paste(plotDIR, "irisPairs.eps", sep="/"),
#           height=4.5, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch3, "irisPairs.tex", sep="/"),
           height=4.5, width=4.5, standAlone=F)
pairs(iris.df[1:4],main = "The Iris Data",
      pch = c(2, 4, 6)[unclass(iris$Species)],
      col = c("#BFBFBF", "#808080", "#404040")[unclass(iris$Species)])
dev.off()

### Figure \ref{fig:piecewise} 3.11
#postscript(file=paste(plotDIR,"nadbplot.eps", sep="/"),
#           width=4.75, height=1.75, horizontal=F)
tikz(file=paste(plotDIRch3,"nadbplot.tex", sep="/"),
           width=4.75, height=1.75, standAlone=F)
par(mfrow=c(1,2), mar=c(2.5,2.5,.5,0.5), mgp=c(1.5, 0.5, 0))
plot(TPOut ~ PLI, data=nadb,
     xlab="P Loading", ylab="P Concentration", cex=0.5)
plot(TPOut ~ PLI, data=nadb,
     xlab="P Loading", ylab="P Concentration", log="x", axes=F, cex=0.5)
axis(1, at=c(0.01, 1, 100), label=c("0.01", "1", "100"))
axis(2)
box()
dev.off()

### Figure \ref{fig:nadbpowers} 3.12
## power transformation
powerT <- function(y, lambda1=c(-1,-1/2,-1/4), lambda2 = c(1/4, 1/2, 1), layout1=2){
  nt <- length(lambda1)+length(lambda2)+1
  transformed <- cbind(outer(y,lambda1,"^"),log(y),(outer(y,lambda2,"^")))
  y.power <- data.frame(transformed=c(transformed),
                        lambda = factor(rep(round(c(lambda1,0,lambda2), 2),
                          rep(length(y),nt))))
  ans <- qqmath(~transformed | lambda,
                data=y.power,
                prepanel = prepanel.qqmathline,
                panel = function(x, ...) {
                  panel.grid(h = 0)
                  panel.qqmath(x, col=gray(0.5))
                  panel.qqmathline(x, distribution = qnorm)
                }, aspect=1, scale = list(y = "free"),
                layout=c(layout1,ceiling(nt/layout1)),
                xlab = "Unit Normal Quantile",
                ylab = "y")
  ans
}

#postscript(file=paste(plotDIR, "nadbpowers.eps", sep="/"),
#           width= 4., height=5, horizontal=F)
tikz(file=paste(plotDIRch3, "nadbpowers.tex", sep="/"),
           width= 4., height=6, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
powerT(y=nadb$PLI)
dev.off()

### conditional plot
### Figure \ref{fig:pm2.5-1} 3.13
#trellis.device(postscript, file=paste(plotDIR, "pm25plot1.eps", sep="/"),
#               height=2, width=3, horizontal=F)
tikz(file=paste(plotDIRch3, "pm25plot1.tex", sep="/"),
               height=2, width=3, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
xyplot(log.value~AvgTemp, panel=function(x,y,...){
  panel.grid()
  panel.xyplot(x,y, col=grey(0.65), cex=0.5, ...)
  panel.loess(x,y,span=1, degree=1,col=1,...)
}, scales=list(x=list(cex=0.75), y=list(cex=0.75)),
       par.settings=trellis.par.temp,
       data=pmdata, xlab="Average Daily Temperature (F)", ylab="Log PM2.5")
dev.off()

## Figure \ref{fig:pm2.5-2}, 3.14
#postscript(file=paste(plotDIR,"pm25coplot2.eps", sep="/"),
#           height=1.25, width=4.75, horizontal=F)
tikz(file=paste(plotDIRch3,"pm25coplot2.tex", sep="/"),
           height=2.25, width=4.75, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
xyplot(log.value~AvgTemp|Month, panel=function(x,y,...){
  panel.grid()
  panel.xyplot(x,y, col=grey(0.65), cex=0.5, ...)
  panel.loess(x,y,span=1, degree=1,col=1,...)
  }, layout=c(12, 1),
       scales=list(x=list(relation="free", cex=0.4,
           alternating=c(1,2))),
       ## x-axis relation and font size
  par.settings=trellis.par.temp,
  data=pmdata, xlab="Average Daily Temperature (F)", ylab="Log PM2.5")
dev.off()

## the air quality data

### Figure \ref{fig:airQcop1}
## traditional coplot
coplot(sqrt(Ozone) ~ Solar.R|Wind*Temp , data=airquality,
       given.values=list(co.intervals(airquality$Wind, 3, 0.25),
         co.intervals(airquality$Temp, 3, 0.25)),
       panel=function(x, y, ...)
       panel.smooth(x, y, span=1, ...),
       ylab="Log Ozone",
       xlab=c("Solar Radiation", "Wind","Temperature"))

WindSpeed <- equal.count(airquality$Wind, 3, 0.25)
Temperature <- equal.count(airquality$Temp, 3, 0.25)
SolarR <- equal.count(airquality$Solar.R, 3, 0.25)


#postscript(file=paste(plotDIR, "airQcoplot1.eps", sep="/"),
#           width=3.75, height=4, horizontal=F)
tikz(file=paste(plotDIRch3, "airQcoplot1.tex", sep="/"),
           width=3.75, height=4, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.set(list(layout.widths=list(left.padding=0, right.padding=0,
                       ylab.axis.padding=0, axis.right=0, key.ylab.padding=0)))
print(
xyplot(sqrt(Ozone) ~ SolarR|WindSpeed*Temperature,
       data=airquality,
       panel=function(x,y,...){
#            panel.loess(x, y, span=1, degree=1, ...)
            panel.grid()
            panel.lmline(x, y, col="grey",...)
            panel.xyplot(x, y, col=1, cex=0.5, ...)
       },
       ylab=list(label="Sqrt Ozone", cex=0.6),
       xlab=list(label="Solar Radiation", cex=0.6),
       scales=list(x=list(alternating=c(1, 2, 1))),
#       between=list(y=1),
       par.strip.text=list(cex=0.4), aspect=1,
       par.settings=list(axis.text=list(cex=0.4)))
)
dev.off()

#postscript(file=paste(plotDIR, "airQcoplot2.eps", sep="/"),
#           width=3.75, height=4, horizontal=F)
tikz(file=paste(plotDIRch3, "airQcoplot2.tex", sep="/"),
           width=3.75, height=4, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.set(list(layout.widths=list(left.padding=0, right.padding=0,
                       ylab.axis.padding=0, axis.right=0, key.ylab.padding=0)))
print(
      xyplot(sqrt(Ozone) ~ Wind|SolarR*Temperature,
             data=airquality,
             panel=function(x,y,...){
#            panel.loess(x, y, span=1, degree=1, ...)
               panel.grid()
               panel.lmline(x, y, col="grey",...)
               panel.xyplot(x, y, col=1, cex=0.5, ...)
             },
             ylab=list(label="Sqrt Ozone", cex=0.6),
             xlab=list(label="Wind Speed", cex=0.6),
             scales=list(x=list(alternating=c(1, 2, 1))),
                                        #       between=list(y=1),
             par.strip.text=list(cex=0.4), aspect=1,
             par.settings=list(axis.text=list(cex=0.4)))
      )
dev.off()

#postscript(file=paste(plotDIR, "airQcoplot3.eps", sep="/"),
#           width=3.75, height=4, horizontal=F)
tikz(file=paste(plotDIRch3, "airQcoplot3.tex", sep="/"),
           width=3.75, height=4, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.set(list(layout.widths=list(left.padding=0, right.padding=0,
                       ylab.axis.padding=0, axis.right=0, key.ylab.padding=0)))
print(
xyplot(sqrt(Ozone) ~ Temp|WindSpeed*SolarR,
       data=airquality,
       panel=function(x,y,...){
#            panel.loess(x, y, span=1, degree=1, ...)
            panel.grid()
            panel.lmline(x, y, col="grey",...)
            panel.xyplot(x, y, col=1, cex=0.5, ...)
       },
       ylab=list(label="Sqrt Ozone", cex=0.6),
       xlab=list(label="Temperature", cex=0.6),
       scales=list(x=list(alternating=c(1, 2, 1))),
#       between=list(y=1),
       par.strip.text=list(cex=0.4), aspect=1,
       par.settings=list(axis.text=list(cex=0.4)))
)
dev.off()

