source("FrontMatter.R")

## chapter 10
plotDIRch10 <- paste(plotDIR, "chapter10", "figures", sep="/")


## seaweed grazers
seaweed <- read.csv(paste(dataDIR, "seaweed.csv", sep="/"), header=T)
seaweed$TREAT <- ordered(seaweed$TREAT, levels=c("CONTROL","f","fF","L","Lf","LfF"))
#seaweed
#   COVER   BLOCK   TREAT

## EDA --
tikz(file=paste(plotDIRch10, "seaweedEDA.tex", sep="/"),
     height=4, width=5.5, standAlone=F)
par(mfrow=c(1,2), mar=c(3,3,0.5,0.25), mgp=c(1.25,0.125,0), tck=0.015, las=1)
plot.design(seaweed)
boxplot(COVER~TREAT, data=seaweed)
dev.off()

slplot <- function(y, x, Xlab=NULL, Ylab=NULL){
    oo <- !is.na(y) & !is.na(x)
    y<-y[oo]
    x <- x[oo]
    mads.ref <- factor(x)
    oneway.output <- oneway(y ~ factor(x), location=mean, spread=1)
    obj <- xyplot(sqrt(abs(residuals(oneway.output)))~jitter(fitted.values(oneway.output), factor=.5),
        aspect=0.75,
        panel=function(x,y){
            panel.xyplot(x,y)
            srmads <- sqrt(tapply(abs(residuals(oneway.output)),
                mads.ref, median))
            oo <- order(oneway.output$location)
            panel.lines(oneway.output$location[oo],srmads[oo])
        },
        xlab=Xlab,
        ylab=Ylab)
    return(obj)
}

obj1 <- slplot(seaweed$COVER, seaweed$TREAT, "Treatment", "\\% Recover")
obj2 <- slplot(log(seaweed$COVER), seaweed$TREAT, "Treatment", "Log Recover")
obj3 <- slplot(logit(seaweed$COVER), seaweed$TREAT, "Treatment", "Logit Recover")
obj4 <- slplot(1/seaweed$COVER, seaweed$TREAT, "Treatment", "Inverse Recover")

tikz(file=paste(plotDIRch10, "seaweedSLs.tex", sep="/"),
     height=4, width=5.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
print(obj1, position=c(0,0,0.5,0.575), more=T)
print(obj2, position=c(0.5,0,1,0.575), more=T)
print(obj3, position=c(0,0.425,0.5,1), more=T)
print(obj4, position=c(0.5,0.425,1,1), more=F)
dev.off()

seaweed$y <- logit(seaweed$COVER/100)
names(seaweed)[2:3] <- c("Block","Treatment")
seaweed.lm <- lm(y ~ factor(Treatment)-1, data=seaweed)
summary(seaweed.lm)
seaweed.aov <- aov(y ~ Treatment, data=seaweed)
summary(seaweed.aov)

tikz(paste(plotDIRch10, "seaweedResd.tex", sep="/"),
     height=4, width=5, standAlone=F)
par(mfrow=c(2,2), tck=0.01, mgp=c(1.25,.125,0), mar=c(3,3,4,0.25))
 plot(seaweed.aov)
dev.off()

require(arm)
seaweed.lmer <- lmer(y ~ 1+(1|Treatment), data=seaweed)
summary(seaweed.lmer)
ranef(seaweed.lmer)
se.ranef(seaweed.lmer)

line.plots <- function(est, se, Ylabel, Xlab, yaxis=2){
    n <- length(est)
    if(n != length(se))stop("lengths not match")
    plot(1:n, 1:n, xlim=range(c(est+2*se, est-2*se)), ylim=c(0.75, n+0.25), type="n", axes=F, xlab=Xlab, ylab="")
    axis(1)
    axis(yaxis, at=1:n, labels=Ylabel, las=1)
    segments(y0=1:n, y1=1:n, x0=est-2*se, x1=est+2*se)
    segments(y0=1:n, y1=1:n, x0=est-1*se, x1=est+1*se, lwd=2.5)
    points(est, 1:n)
    abline(v=0, col="gray")
    invisible()
}

line.plots.compare <- function(est1, se1, est2, se2, Ylabel, Xlab, yaxis=2, V=NULL){
    n <- length(est1)
    if(n != length(se1) | n !=length(se2) | n != length(est2) )stop("lengths not match")
    plot(1:n, 1:n, xlim=range(c(est1+2*se1, est1-2*se1, est2+2*se2, est2-2*se2)), ylim=c(0.75, n+0.25), type="n", axes=F, xlab=Xlab, ylab="")
    axis(1)
    axis(yaxis, at=1:n, labels=Ylabel, las=1)
    segments(y0=(1:n)-0.125, y1=(1:n)-0.125, x0=est1-2*se1, x1=est1+2*se1)
    segments(y0=(1:n)-0.125, y1=(1:n)-0.125, x0=est1-1*se1, x1=est1+1*se1, lwd=2.5)
    points(est1, (1:n)-0.125, pch=16, cex=0.5)

    segments(y0=(1:n)+0.125, y1=(1:n)+0.125, x0=est2-2*se2, x1=est2+2*se2, col="gray")
    segments(y0=(1:n)+0.125, y1=(1:n)+0.125, x0=est2-1*se2, x1=est2+1*se2, lwd=2.5, col="gray")
    points(est2, (1:n)+0.125, cex=0.5, col="gray")
    if(is.null(V))
        abline(v=0, col="gray")
    else abline(v=V, col="gray")
    invisible()
}

lmer1.int<-as.data.frame(cbind(ranef(seaweed.lmer)[[1]][,1]+fixef(seaweed.lmer)[1],
                               se.ranef(seaweed.lmer)[[1]][,1]))

ls2.int<-as.data.frame(summary(seaweed.lm)$coef[, 1:2])

tikz(file=paste(plotDIRch10, "seaweed1.tex", sep="/"),
     width=3, height=2.5, standAlone=F)
par(mar=c(4,5,0.5, 0.75), mgp=c(1.25,0.25,0), tck=0.01)
line.plots.compare(ls2.int[,1], ls2.int[,2],
                   lmer1.int[,1], lmer1.int[,2],
                   levels(seaweed$Treatment), "Treatment effects",
                   V=fixef(seaweed.lmer)[[1]][1])
dev.off()

seaweed.lmer2 <- lmer(y ~ 1+(1|Treatment)+(1|Block), data=seaweed)
summary(seaweed.lmer2)

ranef(seaweed.lmer2)

seaweed.aov2 <- lm(y~Treatment+Block, data=seaweed)
summary(seaweed.aov2)

## mcmcsamp no longer available
##sims.M2 <- mcmcsamp(seaweed.lmer2, n=10000, saveb=T)
##hist(apply(sims.M2@ranef[1:8, ], 2, sd))

##block.mcmc <- sims.M2@ranef[1:8, 5001:10000]
##treat.mcmc <- sims.M2@ranef[9:14,5001:10000]

##sigma.block <- apply(block.mcmc, 2, sd)
##sigma.treat <- apply(treat.mcmc, 2, sd)
##sigma <- sims.M2@sigma[5001:10000]

##s.sum <- rbind(
##               quantile(sigma, prob=c(0.025,0.25,0.5,0.75,0.975)),
##               quantile(sigma.treat, prob=c(0.025,0.25,0.5,0.75,0.975)),
##               quantile(sigma.block, prob=c(0.025,0.25,0.5,0.75,0.975)))
##postscript(file=paste(plotDIR, "seaweedAOV2.eps", sep="/"), width=3.5, height=2.75, horizontal=F)
##par(mar=c(3.5, 7, 0.25,0.25), mgp=c(1.5,0.5,0), tck=-0.02)
##plot(c(0,1), c(0.75,3.25), xlim=range(s.sum), type="n",
##     xlab="standard deviation", ylab="", axes=F)
##abline(h=0, col="gray")
##segments(x0=s.sum[,1], x1=s.sum[,5], y0=1:3, y1=1:3)
##segments(x0=s.sum[,2], x1=s.sum[,4], y0=1:3, y1=1:3, lwd=3)
##axis(1)
##axis(2, at=1:3, labels=c("Residuals","Treatment","Block"), las=1)
##points(x=s.sum[,3], y=1:3, pch=16, cex=1.25)
##dev.off()

seaweed.aov3 <- lm(y~Treatment*Block, data=seaweed)

seaweed.lmer3 <- lmer(y~1 + (1|Treatment)+(1|Block)+(1|Treatment:Block), data=seaweed)

##sims.M3 <- mcmcsamp(seaweed.lmer3, n=10000, saveb=T)
## 1:48 -- interaction
## 49:56 -- block
## 57:62 -- treatment

## ranef
##ranef.int <- sims.M3@ranef[1:48,5001:10000]
##ranef.blk <- sims.M3@ranef[48:56,5001:10000]
##ranef.trt <- sims.M3@ranef[57:62,5001:10000]

##sigma.int <- apply(ranef.int, 2, sd)
##sigma.blk <- apply(ranef.blk, 2, sd)
##sigma.trt <- apply(ranef.trt, 2, sd)
##sigma <- sims.M3@sigma[5001:10000]

##s.sum <- rbind(
##               quantile(sigma, prob=c(0.025,0.25,0.5,0.75,0.975)),
##               quantile(sigma.int, prob=c(0.025,0.25,0.5,0.75,0.975)),
##               quantile(sigma.trt, prob=c(0.025,0.25,0.5,0.75,0.975)),
##               quantile(sigma.blk, prob=c(0.025,0.25,0.5,0.75,0.975)))
##postscript(file=paste(plotDIR, "seaweedAOV3.eps", sep="/"), width=3.5, height=2.75, horizontal=F)
##par(mar=c(3.5, 7, 0.25,0.25), mgp=c(1.5,0.5,0), tck=-0.02)
##plot(c(0,1), c(0.75,4.25), xlim=range(s.sum), type="n",
##     xlab="standard deviation", ylab="", axes=F)
##abline(h=0, col="gray")
##segments(x0=s.sum[,1], x1=s.sum[,5], y0=1:4, y1=1:4)
##segments(x0=s.sum[,2], x1=s.sum[,4], y0=1:4, y1=1:4, lwd=3)
##axis(1)
##axis(2, at=1:4, labels=c("Residuals","Interaction","Treatment","Block"), las=1)
##points(x=s.sum[,3], y=1:4, pch=16, cex=1.25)
##dev.off()

##sigma.int.plot <- list()
##sigma.int.plot$summary <- t(apply(ranef.int, 1, FUN=function(x){
##    return(c(mean(x), sd(x), quantile(x, prob=c(0.025,0.25,0.5,0.75,0.975))))
##}))

##postscript(file=paste(plotDIR, "seaweedInt.eps", sep="/"), width=5, height=3.5, horizontal=F)
##par(mgp=c(1.5,0.5,0), tck=-0.02)
##summary.plot.Interaction(sigma.int.plot, rows=1:48,
##                         treatment = levels(seaweed$Treatment),
##                         block=1:8)
##dev.off()
##summary.plot.Interaction<-
##function(bugs.out=bugs.out.S, rows, ylab = NULL, xlab=" ",
##         treatment, block, reverse=F, ymar=7, ...){
##
## Graphical presentation of interaction effect from bugs output
## 'bugs.out' is generated by function 'bugs' from linrary(R2WinBUGS)
##
##    Plot.data <- bugs.out$summary[rows,]
##    plotting.region <- range(Plot.data[,c(3,7)])
#    treatment  <- levels(seaweed$TREAT)
##    nt <- length(treatment)
#    block <- paste("Block", 1:8)
##    nb <- length(block)
##    if(reverse){
##    par(mfrow=c(1, nt), oma=c(0.5,ymar,1,1), mar=c(5, 0, 0, 0))
##    for (i in 1:nt){
##    plot(seq(plotting.region[1], plotting.region[2],,5),
##         seq(1-0.1, nb+.1, ,5), type="n",
##         xlab=treatment[i], ylab=" ", axes=F, ...)
##    axis(1)
##    if (i==1)
##        axis(2, at=1:nb, labels=block, las=1, outer=T, ...)
##    segments(x0=Plot.data[((i-1)*nb+1):(i*nb),3],
##             x1=Plot.data[((i-1)*nb+1):(i*nb),7], y0=1:nb, y1=1:nb)
##    segments(x0=Plot.data[((i-1)*nb+1):(i*nb),4],
##             x1=Plot.data[((i-1)*nb+1):(i*nb),6], y0=1:nb, y1=1:nb, lwd=3)
##    abline(v=0, col="gray")
##    points(Plot.data[((i-1)*nb+1):(i*nb),1], 1:nb, cex=0.75)
##    }
##    } else {
##    par(mfrow=c(1, nb), oma=c(5,7,4,4), mar=c(3.5, 0, 0, 0))
##    for (i in 1:nb){
##    plot(seq(plotting.region[1], plotting.region[2],,5),
##         seq(1-0.1, nt+.1, ,5), type="n",
##        xlab=block[i], ylab=" ", axes=F)
##    axis(1)
##    if (i==1)
##    axis(2, at=1:nt, labels=treatment, las=1, outer=T)
##    segments(x0=Plot.data[seq(i,nt*nb,nb),3],
##             x1=Plot.data[seq(i,nt*nb,nb),7], y0=1:nt, y1=1:nt)
##    segments(x0=Plot.data[seq(i,nt*nb,nb),4],
##             x1=Plot.data[seq(i,nt*nb,nb),6], y0=1:nt, y1=1:nt, lwd=3)
##    abline(v=0, col="gray")
##    points(Plot.data[seq(i,nt*nb,nb),1], 1:nt, cex=0.75)
##    }
##    }
##    mtext("Interaction Effects", side=1, outer=T, line=-.5)
##    invisible()
##}

seaweed.lmer4 <- lmer(COVER ~ 1 + (1|Treatment)+(1|Block)+(1|Treatment:Block), data=seaweed)

## block comments: mark the block and c-x r t
##sims.M3 <- mcmcsamp(seaweed.lmer4, n=10000, saveb=T)


## 1:48 -- interaction
## 49:56 -- block
## 57:62 -- treatment

## ranef
##ranef.int <- sims.M3@ranef[1:48,5001:10000]
##ranef.blk <- sims.M3@ranef[48:56,5001:10000]
##ranef.trt <- sims.M3@ranef[57:62,5001:10000]

##sigma.int <- apply(ranef.int, 2, sd)
##sigma.blk <- apply(ranef.blk, 2, sd)
##sigma.trt <- apply(ranef.trt, 2, sd)
##sigma <- sims.M3@sigma[5001:10000]
##
##rbind(
##      quantile(sigma, prob=c(0.025,0.25,0.5,0.75,0.975)),
##      quantile(sigma.int, prob=c(0.025,0.25,0.5,0.75,0.975)),
##      quantile(sigma.trt, prob=c(0.025,0.25,0.5,0.75,0.975)),
##      quantile(sigma.blk, prob=c(0.025,0.25,0.5,0.75,0.975)))
##t(file=paste(plotDIR, "seaweedAOV4.eps", sep="/"), width=3.5, height=2.75, horizontal=F)
##(3.5, 7, 0.25,0.25), mgp=c(1.5,0.5,0), tck=-0.02)
##1), c(0.75,4.25), xlim=range(s.sum), type="n",
##="standard deviation", ylab="", axes=F)
##0, col="gray")
##x0=s.sum[,1], x1=s.sum[,5], y0=1:4, y1=1:4)
##x0=s.sum[,2], x1=s.sum[,4], y0=1:4, y1=1:4, lwd=3)
##
##t=1:4, labels=c("Residuals","Interaction","Treatment","Block"), las=1)
##s.sum[,3], y=1:4, pch=16, cex=1.25)
##
##
##.plot <- list()
##.plot$summary <- t(apply(ranef.int, 1, FUN=function(x){
##n(c(mean(x), sd(x), quantile(x, prob=c(0.025,0.25,0.5,0.75,0.975))))
##
##
##t(file=paste(plotDIR, "seaweedInt2.eps", sep="/"), width=5, height=3.5, horizontal=F)
##(1.5,0.5,0), tck=-0.02)
##lot.Interaction(sigma.int.plot, rows=1:48,
##                treatment = levels(seaweed$Treatment),
##                block=1:8)
##

#
################ N2O emission data #####################
N2O.data <- read.csv(paste(dataDIR, "N2Oemission.csv", sep="/"), header=T)
N2O.data$y <- log(N2O.data$emission+1)
N2O.data$x <- log(N2O.data$n.input+1)
N2O.data$gf <- paste(N2O.data$group, N2O.data$fertilizer)
N2O.data$fertilizer <- as.character(N2O.data$fertilizer)
N2O.data$x2 <- N2O.data$sample.period/30  ## time in month

N2O.background <- N2O.data[N2O.data$n.input==0 & !is.na(N2O.data$carbon),]
N2O.background$group <- as.numeric(ordered(N2O.background$group))

bckg.lm <- lm(log(emission/x2) ~ factor(group)-1, data=N2O.background)
bckg.aov <- aov(log(emission/x2) ~ factor(group), data=N2O.background)

bckg.lmer <- lmer(log(emission/x2) ~ 1 + (1|group), data=N2O.background)

bckg.size <- as.vector(table(N2O.background$group))
size <- bckg.size + runif(length(bckg.size), -0.1,0.1)

bckg.lmCoef <- summary(bckg.lm)$coef

## Figure 10.2
##postscript(file=paste(plotDIR, "bckgCompare.eps", sep="/"), width=5, height=2.5, horizontal=F)
tikz(file=paste(plotDIRch10, "bckgCompare.tex", sep="/"),
     width=5, height=2.5, standAlone=F)
lower <- bckg.lmCoef[,1] - bckg.lmCoef[,2]
upper <- bckg.lmCoef[,1] + bckg.lmCoef[,2]
Ylim <- range(lower, upper)

par(mfrow=c(1,2), mgp=c(1.25,0.125,0), mar=c(3,3,0.5,0.25), las=1, tck=0.01)
plot(size, rnorm(length(size)), ylim=Ylim, type="n", xlab="Sample size",
  ylab="Log Mean (No Pooling)", log="x")
abline(h=mean(bckg.lmCoef[,1]))
segments(x0=size, x1=size, y0=lower, y1=upper)
points(size, bckg.lmCoef[,1], pch=16, cex=0.5)
#title(main="No Pooling", cex=0.75)

lmer.mean <- fixef(bckg.lmer)[1] + ranef(bckg.lmer)[[1]][,1]
lmer.se <- sqrt(se.fixef(bckg.lmer)[1]^2 + se.ranef(bckg.lmer)[[1]][,1]^2)
lower <- lmer.mean-lmer.se
upper <- lmer.mean+lmer.se

plot(size, rnorm(length(size)), ylim=Ylim, type="n", xlab="Sample size",
  ylab="Log Mean (Partial P)", log="x")
abline(h=fixef(bckg.lmer)[[1]][1])
segments(x0=size, x1=size, y0=lower, y1=upper)
points(size, lmer.mean, pch=16,cex=0.5)
dev.off()


###### group level predictor ######
carbon.group <- tapply(N2O.background$carbon/100, N2O.background$group, mean, na.rm=T)
carbon.full <- carbon.group[N2O.background$group]
carb <- logit(carbon.group)

## Figure 10.3
##postscript (paste(plotDIR, "n2ologitcarb.eps", sep="/"),
##            width=4, height=2.75, horizontal=F)
##
tikz (paste(plotDIRch10, "n2ologitcarb.tex", sep="/"),
            width=4, height=2.75, standAlone=F)
par (mar=c(3,3,3,0.25), mgp=c(1.25,.125,0), mfrow=c(1,2), las=1, tck=0.01)
hist(carbon.group*100, xlab="Soil C (\\%)", main="")
hist(carb, xlab="Logit soil C", main="")
dev.off()

bckg.lmer2 <- lmer(log(emission/x2) ~ 1 + logit(carbon.full) + (1|group), data=N2O.background)

b0.hat.M2 <- coef(bckg.lmer2)$group[,1] + coef(bckg.lmer2)$group[,2] * logit(carbon.group)
b0.se.M2 <- se.ranef(bckg.lmer2)$group[,1]
lower <- b0.hat.M2 -b0.se.M2
upper <- b0.hat.M2 + b0.se.M2

## Figure 10.4
##postscript (paste(plotDIR, "bckgGroup.eps", sep="/"), width=3.5, height=3, horizontal=F)
tikz(paste(plotDIRch10, "bckgGroup.tex", sep="/"),
     width=3.5, height=3, standAlone=F)
par (mar=c(3,3,.25,0.25), mgp=c(1.25,.125,0), tck=0.01, las=1)
plot (logit(carbon.group), b0.hat.M2, xlab="group-level soil carbon", ylab="log mean emission", pch=20,
      ylim=range(lower, upper), cex=0.5)
curve (fixef(bckg.lmer2)[1] + fixef(bckg.lmer2)[2]*x, lwd=1, add=TRUE)
segments(x0=logit(carbon.group), x1=logit(carbon.group), y0=lower, y1=upper, lwd=0.5)
dev.off()


#### EUSE ####

rtol2 <- read.csv(file=paste(dataDIR, "rtolforMS.csv", sep="/"), header=T)

packages(car)

## Figure 10.10
##postscript(file=paste(plotDIR, "richtol1.eps", sep="/"), width=4.5, height=4.5, horizontal=F)
tikz(file=paste(plotDIRch10, "richtol1.tex", sep="/"),
     width=4.5, height=4.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
xyplot(richtol~nuii|city, data=rtol2, panel=function(x,y,...){
    panel.xyplot(x, y, ...)
    panel.lmline(x,y,...)}, ylab="TOLr", xlab="NUII")
dev.off()

## environmental data
euse.env <- read.csv(paste(dataDIR, "EUSE_NAT_ENV.csv", sep="/"), header=T)
names(euse.env)[13]<-"MAX.ELEV"
names(euse.env)[12]<-"MIN.ELEV"

AvePrec <- tapply(euse.env$AnnMeanP, euse.env$CITY, mean)
AveTemp <- tapply(euse.env$AnnMeanT, euse.env$CITY, mean)
AveElev <- tapply(euse.env$MEANELEV, euse.env$CITY, mean)
AveMaxT <- tapply(euse.env$AnnMaxT, euse.env$CITY, mean)
AveMaxP <- tapply(euse.env$AnnMaxP, euse.env$CITY, mean)
AveMinP <- tapply(euse.env$AnnMinP, euse.env$CITY, mean)
AvePdif <- tapply(euse.env$AnnMaxP-euse.env$AnnMinP, euse.env$CITY, mean)

site <- as.numeric(ordered(rtol2$city))
prec.full <- as.vector(AvePrec[site])
temp.full <- as.vector(AveTemp[site])

prep.full <- AvePrec[site]
temp.full <- AveTemp[site]
elev.full <- AveElev[site]
maxp.full <- AveMaxP[site]
minp.full <- AveMinP[site]
maxt.full <- AveMaxT[site]
pdif.full <- AvePdif[site]

citygrp<-sort(unique(rtol2$city))
## background agriculture land use
city_ag<-read.csv(paste(dataDIR, "City_AG_Grassland.csv", sep="/"), header=T, na.strings = ".")
city_ag<-cbind(city_ag[,1:2],city_ag[,3:5]/100)
city_ag[order(city_ag[ ,1]), ]
ag<-city_ag[order(city_ag[ ,1]), ][,5]
ag.full <- as.vector(ag[site])
ag.cat <- ag.full>0.5

## complete pooling
euse.lm1 <- lm(richtol ~ nuii, data=rtol2)
display(euse.lm1, 4)

## no pooling
euse.lm2 <- lm(richtol ~ nuii*factor(city)-1-nuii, data=rtol2)
display(euse.lm2, 4)

ls2.int<-as.data.frame(rbind(summary(euse.lm2)$coef[1:9, 1:2], summary(euse.lm1)$coef[1, 1:2]))
ls2.slp<-as.data.frame(rbind(summary(euse.lm2)$coef[10:18, 1:2], summary(euse.lm1)$coef[2, 1:2]))

## Figure 10.11
## postscript(file=paste(plotDIR, "EUSEcoefLM.eps", sep="/"), width=4.5, height=4, horizontal=F)
tikz(file=paste(plotDIRch10, "EUSEcoefLM.tex", sep="/"),
     width=4.5, height=4, standAlone=F)
par(mfrow=c(1,2), mar=c(4,4,0.5, 0.75), mgp=c(1.25,0.125,0), tck=0.01)
line.plots(ls2.int[,1], ls2.int[,2], c(levels(rtol2$city), "All"), "Intercepts", 2)
par(mar=c(4,0.75,0.5,4))
line.plots(ls2.slp[,1], ls2.slp[,2], c(levels(rtol2$city), "All"), "Slopes", 4)
dev.off()

## partial pooling
euse.lmer1 <- lmer(richtol ~ nuii+(1+nuii|city), data=rtol2)
summary(euse.lmer1)

ranef(euse.lmer1)
fixef(euse.lmer1)

lmer1.int<-as.data.frame(rbind(cbind(ranef(euse.lmer1)[[1]][,1]+fixef(euse.lmer1)[1], se.ranef(euse.lmer1)[[1]][,1]),
                               c(fixef(euse.lmer1)[1], se.fixef(euse.lmer1)[1])))
lmer1.slp<-as.data.frame(rbind(cbind(ranef(euse.lmer1)[[1]][,2]+fixef(euse.lmer1)[2], se.ranef(euse.lmer1)[[1]][,2]),
                               c(fixef(euse.lmer1)[2], se.fixef(euse.lmer1)[2])))

##postscript(file=paste(plotDIR, "EUSEcoefLMER.eps", sep="/"), width=4.5, height=4, horizontal=F)
tikz(file=paste(plotDIRch10, "EUSEcoefLMER.tex", sep="/"),
     width=4.5, height=4, standAlone=F)
par(mfrow=c(1,2), mar=c(4,4,0.5, 0.75), mgp=c(1.25,0.125,0), tck=0.01)
line.plots(lmer1.int[,1], lmer1.int[,2], c(levels(rtol2$city), "All"), "Intercepts", 2)
par(mar=c(4,0.75,0.5,4))
line.plots(lmer1.slp[,1], lmer1.slp[,2], c(levels(rtol2$city), "All"), "Slopes", 4)
dev.off()

## Figure 10.12
tikz(file=paste(plotDIRch10, "EUSEcompare.tex", sep="/"),
     width=4.5, height=3.5, standAlone=F)
par(mfrow=c(1,2), mar=c(4,4,0.5, 0.75), mgp=c(1.25,0.125,0), tck=0.01)
line.plots.compare(ls2.int[,1], ls2.int[,2], lmer1.int[,1], lmer1.int[,2], c(levels(rtol2$city), "All"), "Intercepts", 2)
par(mar=c(4,0.75,0.5,4))
line.plots.compare(ls2.slp[,1], ls2.slp[,2], lmer1.slp[,1], lmer1.slp[,2], c(levels(rtol2$city), "All"), "Slopes", 4)
dev.off()

euse.lmer2 <- lmer(richtol ~ nuii + temp.full + nuii:temp.full + (1+nuii|city), data=rtol2)
M2.coef <- coef (euse.lmer2)
a.hat.M2 <- M2.coef[[1]][,1] + M2.coef[[1]][,3]*AveTemp
b.hat.M2 <- M2.coef[[1]][,2] + M2.coef[[1]][,4]*AveTemp
a.se.M2 <- se.ranef(euse.lmer2)[[1]][,1]
b.se.M2 <- se.ranef(euse.lmer2)[[1]][,2]

# plot estimated intercepts and slopes

## Figure 10.13
tikz(file=paste(plotDIRch10, "richtol2Temp.tex", sep="/"),
     width=4.75, height=3, standAlone=F)
par (mfrow=c(1,2), mar=c(3,3,3,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
lower <- a.hat.M2 - a.se.M2
upper <- a.hat.M2 + a.se.M2
plot (AveTemp, a.hat.M2, ylim=range(lower,upper), cex.lab=0.75, cex.axis=0.75,
      xlab="average temperature", ylab="regression intercept", pch=20)
curve (fixef(euse.lmer2)["(Intercept)"] + fixef(euse.lmer2)["temp.full"]*x, lwd=1, col="black", add=TRUE)
segments (AveTemp, lower, AveTemp, upper, lwd=.5, col="gray10")
text(AveTemp, lower, levels(rtol2$city), adj=c(.5,1),cex=0.5)

lower <- b.hat.M2 - b.se.M2
upper <- b.hat.M2 + b.se.M2

plot (AveTemp, b.hat.M2, ylim=range(lower,upper), cex.lab=0.75, cex.axis=0.75,
      xlab="average temperature", ylab="regression slope", pch=20)
curve (fixef(euse.lmer2)["nuii"] + fixef(euse.lmer2)["nuii:temp.full"]*x, lwd=1, col="black", add=TRUE)
segments (AveTemp, lower, AveTemp, upper, lwd=.5, col="gray10")
text(AveTemp, lower, levels(rtol2$city), adj=c(.5,1),cex=0.5)

dev.off()

## adding background agriculture land use
euse.lmer2 <- lmer(richtol ~ nuii+ag.full+nuii:ag.full+(1+nuii|site), data=rtol2)
summary(euse.lmer2)  ## display no plonger works for lmer

M2.coef <- coef (euse.lmer2)
a.hat.M2 <- M2.coef[[1]][,1] + M2.coef[[1]][,3]*ag
b.hat.M2 <- M2.coef[[1]][,2] + M2.coef[[1]][,4]*ag
a.se.M2 <- se.ranef(euse.lmer2)[[1]][,1]
b.se.M2 <- se.ranef(euse.lmer2)[[1]][,2]

# plot estimated intercepts and slopes

## Figure 10.14
##postscript(file=paste(plotDIR, "richtol2.eps", sep="/"), width=4.75, height=3, horizontal=F)
tikz(file=paste(plotDIRch10, "richtol2.tex", sep="/"),
     width=4.75, height=3, standAlone=F)
par (mfrow=c(1,2), mar=c(3,3,3,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
lower <- a.hat.M2 - a.se.M2
upper <- a.hat.M2 + a.se.M2
plot (ag, a.hat.M2, ylim=range(lower,upper), cex.lab=0.75, cex.axis=0.75,
      xlab="Background ag", ylab="regression intercept", pch=20)
curve (fixef(euse.lmer2)["(Intercept)"] + fixef(euse.lmer2)["ag.full"]*x, lwd=1, col="black", add=TRUE)
segments (ag, lower, ag, upper, lwd=.5, col="gray10")
text(ag, lower, levels(rtol2$city), adj=c(.5,1),cex=0.5)

lower <- b.hat.M2 - b.se.M2
upper <- b.hat.M2 + b.se.M2

plot (ag, b.hat.M2, ylim=range(lower,upper), cex.lab=0.75, cex.axis=0.75,
      xlab="Background ag", ylab="regression slope", pch=20)
curve (fixef(euse.lmer2)["nuii"] + fixef(euse.lmer2)["nuii:ag.full"]*x, lwd=1, col="black", add=TRUE)
segments (ag, lower, ag, upper, lwd=.5, col="gray10")
text(ag, lower, levels(rtol2$city), adj=c(.5,1),cex=0.5)

dev.off()

euse.lmer3 <- lmer(richtol ~ nuii+temp.full+nuii:temp.full+(1+nuii|site)+(1+nuii|ag.cat), data=rtol2)
display(euse.lmer3, 4)

M3.fixef <- fixef (euse.lmer3)
M3.ranef <- ranef (euse.lmer3)
M3.sefixef <- se.fixef (euse.lmer3)
M3.seranef <- se.ranef (euse.lmer3)

a.hat.M3 <- M3.fixef[1] + M3.ranef[[1]][,1] + M3.ranef[[2]][c(1,1,1,2,2,2,1,1,1),1] + M3.fixef[3]*AveTemp
b.hat.M3 <- M3.fixef[2] + M3.ranef[[1]][,2] + M3.ranef[[2]][c(1,1,1,2,2,2,1,1,1),2] + M3.fixef[4]*AveTemp

a.se.M3 <- M3.seranef[[1]][,1]
b.se.M3 <- M3.seranef[[1]][,2]

# plot estimated intercepts and slopes

## Figure 10.15
tikz(file=paste(plotDIRch10, "richtol3.tex", sep="/"),
     width=4.75, height=3, standAlone=F)
par (mfrow=c(1,2), mar=c(3,3,3,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
lower <- a.hat.M3 - a.se.M3
upper <- a.hat.M3 + a.se.M3
plot (AveTemp, a.hat.M3, ylim=range(lower,upper), cex.lab=0.75, cex.axis=0.75,
      xlab="Ave Temp", ylab="regression intercept", pch=20)
curve (fixef(euse.lmer3)[1] + fixef(euse.lmer3)[3]*x + M3.ranef[[2]][1,1], lwd=1, col="black", add=TRUE)
curve (fixef(euse.lmer3)[1] + fixef(euse.lmer3)[3]*x + M3.ranef[[2]][2,1], lwd=1, col="black", add=TRUE)
segments (AveTemp, lower, AveTemp, upper, lwd=.5, col="gray10")
text(AveTemp, lower, levels(rtol2$city), adj=c(.5,1),cex=0.5)

lower <- b.hat.M3 - b.se.M3
upper <- b.hat.M3 + b.se.M3

plot (AveTemp, b.hat.M3, ylim=range(lower,upper), cex.lab=0.75, cex.axis=0.75,
      xlab="Ave Temp", ylab="regression slope", pch=20)
curve (fixef(euse.lmer3)[2] + fixef(euse.lmer3)[4]*x + M3.ranef[[2]][1,2], lwd=1, col="black", add=TRUE)
curve (fixef(euse.lmer3)[2] + fixef(euse.lmer3)[4]*x + M3.ranef[[2]][2,2], lwd=1, col="black", add=TRUE)
segments (AveTemp, lower, AveTemp, upper, lwd=.5, col="gray10")
text(AveTemp, lower, levels(rtol2$city), adj=c(.5,1),cex=0.5)

dev.off()

## alternative models with both temperature and background ag as group level predictors

## Alt 1: using ag.cat as a dummy variable -- additive group level predictors
euse.lmer3.25 <-
    lmer(richtol ~ nuii+temp.full+nuii:temp.full+
         as.numeric(ag.cat)+ as.numeric(ag.cat):nuii +
         (1+nuii|site), data=rtol2)

## Alt 2: using ag.cat as a dummy variable -- consider interaction among group level predictors
euse.lmer3.26 <-
    lmer(richtol ~ nuii+temp.full+nuii:temp.full+
         as.numeric(ag.cat)+ as.numeric(ag.cat):nuii +
         as.numeric(ag.cat):temp.full +
         as.numeric(ag.cat):temp.full:nuii +
         (1+nuii|site), data=rtol2)

## Alt 3: same as Alt 2 but using non-nested groups
euse.lmer3.5 <- lmer(richtol ~ nuii+temp.full+nuii:temp.full+
                     (1+nuii|site)+
                     (1+nuii+temp.full+nuii:temp.full|ag.cat), data=rtol2)
summary(euse.lmer3.5)

M3.fixef <- fixef (euse.lmer3.5)
M3.ranef <- ranef (euse.lmer3.5)
M3.sefixef <- se.fixef (euse.lmer3.5)
M3.seranef <- se.ranef (euse.lmer3.5)

a.hat.M3 <- M3.fixef[1] + M3.ranef[[1]][,1] + M3.ranef[[2]][c(1,1,1,2,2,2,1,1,1),1] + (M3.fixef[3]+ M3.ranef[[2]][c(1,1,1,2,2,2,1,1,1), 3])*AveTemp
b.hat.M3 <- M3.fixef[2] + M3.ranef[[1]][,2] + M3.ranef[[2]][c(1,1,1,2,2,2,1,1,1),2] + (M3.fixef[4]+ M3.ranef[[2]][c(1,1,1,2,2,2,1,1,1), 4])*AveTemp

a.se.M3 <- M3.seranef[[1]][,1]
b.se.M3 <- M3.seranef[[1]][,2]

## Figure 10.16
tikz(file=paste(plotDIRch10, "richtol35.tex", sep="/"),
     width=4.75, height=3, standAlone=F)
par (mfrow=c(1,2), mar=c(3,3,3,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
lower <- a.hat.M3 - a.se.M3
upper <- a.hat.M3 + a.se.M3
plot (AveTemp, a.hat.M3, ylim=range(lower,upper), cex.lab=0.75, cex.axis=0.75,
      xlab="Ave Temp", ylab="regression intercept", pch=20)
curve (M3.fixef[1] + (M3.fixef[3]+M3.ranef[[2]][1,3])*x + M3.ranef[[2]][1,1], lwd=1, col="black", add=TRUE)
curve (M3.fixef[1] + (M3.fixef[3]+M3.ranef[[2]][2,3])*x + M3.ranef[[2]][2,1], lwd=1, col="black", add=TRUE)
segments (AveTemp, lower, AveTemp, upper, lwd=.5, col="gray10")
text(AveTemp, lower, levels(rtol2$city), adj=c(.5,1),cex=0.5)

lower <- b.hat.M3 - b.se.M3
upper <- b.hat.M3 + b.se.M3

plot (AveTemp, b.hat.M3, ylim=range(lower,upper), cex.lab=0.75, cex.axis=0.75,
      xlab="Ave Temp", ylab="regression slope", pch=20)
curve (M3.fixef[2] + (M3.fixef[4]+M3.ranef[[2]][1,4])*x + M3.ranef[[2]][1,2], lwd=1, col="black", add=TRUE)
curve (M3.fixef[2] + (M3.fixef[4]+M3.ranef[[2]][2,4])*x + M3.ranef[[2]][2,2], lwd=1, col="black", add=TRUE)
segments (AveTemp, lower, AveTemp, upper, lwd=.5, col="gray10")
text(AveTemp, lower, levels(rtol2$city), adj=c(.5,1),cex=0.5)

dev.off()



### The Finnish Lakes example in Chapter 10 Multilevel Regression###
#dataDIR <- "c:/users/Oldfiles/Finn/updating/Data"

summer.All <- read.table(paste(dataDIR, "summerAll.csv", sep="/"),
                         sep=",",header=T)

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
## lmer
Finn.M1 <- lmer(y ~ lxp + lxn + (1|type), data=summer.All)
ranef(Finn.M1)

Finn.M2 <- lmer(y ~ lxp + lxn + (1+lxp + lxn|type), data=summer.All)

Finn.M3 <- lmer(y ~ lxp+lxn+lxp:lxn+(1+lxp+lxn+lxp:lxn|type), data=summer.All)
summary(Finn.M3)

ranef(Finn.M3)$type
se.ranef(Finn.M3)$type

fixef(Finn.M3)
se.fixef(Finn.M3)

est <- t(fixef(Finn.M3) + t(as.matrix(ranef(Finn.M3)$type)))
se <- sqrt(t(se.fixef(Finn.M3)^2+t(as.matrix(se.ranef(Finn.M3)$type))^2))

line.plots <- function(est, se, Ylabel, Xlab, yaxis=NULL, hline=0, oo=NULL){
    n <- length(est)
    if (!is.null(oo)) {
        est<-est[oo]
        se <-se[oo]
    }
    if(n != length(se))stop("lengths not match")
    plot(1:n, 1:n, xlim=range(c(est+2*se, est-2*se)), ylim=c(0.75, n+0.25),
         type="n", axes=F, xlab=Xlab, ylab="")
    axis(1)
    if (!is.null(yaxis))
      axis(yaxis, at=1:n, labels=Ylabel, las=1)
    segments(y0=1:n, y1=1:n, x0=est-2*se, x1=est+2*se)
    segments(y0=1:n, y1=1:n, x0=est-1*se, x1=est+1*se, lwd=2.5)
    points(est, 1:n)
    abline(v=hline, col="gray")
    invisible()
}

line.plots.compare <- function(est1, se1, est2, se2, Ylabel, Xlab, yaxis=2){
    n <- length(est1)
    if(n != length(se1) | n !=length(se2) | n != length(est2) )stop("lengths not match")
    plot(1:n, 1:n, xlim=range(c(est1+se1, est1-se1, est2+se2, est2-se2)), ylim=c(0.75, n+0.25), type="n", axes=F, xlab=Xlab, ylab="")
    axis(1)
    axis(yaxis, at=1:n, labels=Ylabel, las=1)
#    segments(y0=(1:n)-0.125, y1=(1:n)-0.125, x0=est1-2*se1, x1=est1+2*se1)
    segments(y0=(1:n)-0.125, y1=(1:n)-0.125, x0=est1-1*se1, x1=est1+1*se1, lwd=2)
    points(est1, (1:n)-0.125, pch=16, cex=0.5)

#    segments(y0=(1:n)+0.125, y1=(1:n)+0.125, x0=est2-2*se2, x1=est2+2*se2, lty=2)
    segments(y0=(1:n)+0.125, y1=(1:n)+0.125, x0=est2-1*se2, x1=est2+1*se2, lwd=2, col="gray")
    points(est2, (1:n)+0.125, cex=0.5, col="gray")

    abline(v=0, col="gray")
    invisible()
}

## Figure 10.17
tikz(file=paste(plotDIRch10, "finnmultcoef.tex", sep="/"),
           width=4.5, height=3, standAlone=F)
par(mfrow=c(1,4), mgp=c(1.25,0.125,0), tck=0.01, las=1)
par(mar=c(3, 3, 0.5,0))
line.plots(est[,1], se[,1], Ylabel=1:9, Xlab="$\\beta_0$",
           yaxis=2, hline=fixef(Finn.M3)[1])
box(col=grey(0.3))
par(mar=c(3, 1.5,0.5,1.5))
line.plots(est[,2], se[,2], Ylabel=1:9, Xlab="$\\beta_1$",
           yaxis=4, hline=fixef(Finn.M3)[2])
box(col=grey(0.3))
line.plots(est[,3], se[,3], Ylabel=1:9, Xlab="$\\beta_2$",
           yaxis=2, hline=fixef(Finn.M3)[3])
box(col=grey(0.3))
par(mar=c(3,0, 0.5, 3))
line.plots(est[,4], se[,4], Ylabel=1:9, Xlab="$\\beta_3$",
           yaxis=4)
box(col=grey(0.3))
dev.off()

## Figure 10.18
lake3 <- summer.All[summer.All$type==1,]
tikz(file=paste(plotDIRch10, "Finncoplottype11.tex", sep="/"),
           width=4.5, height=4, standAlone=F)
par(mar=c(3.,3.,3.,.25), mgp=c(1.25,0.125,0), tck=0.01)
given.tn <- co.intervals(lake3$lxn, number=4, overlap=.1)
 coplot(y ~ lxp | lxn, data = lake3, given.v=given.tn, rows=1,
        xlab="log TP", ylab="log Chla")
dev.off()

## Figure 10.19
tikz(file=paste(plotDIRch10, "Finncoplottype12.tex", sep="/"),
           width=4.5, height=4, standAlone=F)
par(mar=c(3.,3.,3.,.25), mgp=c(1.25,0.125,0), tck=0.01)
given.tp <- co.intervals(lake3$lxp, number=4, overlap=.1)
 coplot(y ~ lxn | lxp, data = lake3, given.v=given.tp, rows=1,
        xlab="log TN", ylab="log Chla")
dev.off()

lake4 <- summer.All[summer.All$type==2&summer.All$lxn>-2,]
##tikz(file=paste(plotDIRch10, "Finncoplottype21.tex", sep="/"),
##           width=4.5, height=4, standAlone=F)
postscript(paste(plotDIRch10,  "Finncoplottype21.eps", sep="/"),
          width=4.5, height=4, horizontal=F)
par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), tck=0.01)
given.tn <- co.intervals(lake4$lxn, number=4, overlap=.1)
 coplot(y ~ lxp | lxn, data = lake4, given.v=given.tn, rows=1,
        xlab="log TP", ylab="log Chla")
dev.off()

## Figure 10.24
##tikz(file=paste(plotDIRch10, "Finncoplottype22.tex", sep="/"),
##           width=4.5, height=4, standAlone=F)
postscript(paste(plotDIRch10,  "Finncoplottype22.eps", sep="/"),
          width=4.5, height=4, horizontal=F)
par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), tck=0.01)
given.tp <- co.intervals(lake4$lxp, number=4, overlap=.1)
 coplot(y ~ lxn | lxp, data = lake4, given.v=given.tp, rows=1,
        xlab="log TN", ylab="log Chla")
dev.off()


## Figure 10.25
lake4 <- summer.All[summer.All$type==3,]
tikz(file=paste(plotDIRch10, "Finncoplottype31.tex", sep="/"),
           width=4.5, height=4, standAlone=F)
par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), tck=0.01)
given.tn <- co.intervals(lake4$lxn, number=4, overlap=.1)
 coplot(y ~ lxp | lxn, data = lake4, given.v=given.tn, rows=1,
        xlab="log TP", ylab="log Chla")
dev.off()

tikz(file=paste(plotDIRch10, "Finncoplottype32.tex", sep="/"),
           width=4.5, height=4, standAlone=F)
par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), tck=0.01)
given.tp <- co.intervals(lake4$lxp, number=4, overlap=.1)
 coplot(y ~ lxn | lxp, data = lake4, given.v=given.tp, rows=1,
        xlab="log TN", ylab="log Chla")
dev.off()

## Figure 10.20
lake4 <- summer.All[summer.All$type==6,]
tikz(file=paste(plotDIRch10, "Finncoplottype61.tex", sep="/"),
           width=4.5, height=4, standAlone=F)
par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), tck=0.01)
given.tn <- co.intervals(lake4$lxn, number=4, overlap=.1)
 coplot(y ~ lxp | lxn, data = lake4, given.v=given.tn, rows=1,
        xlab="log TP", ylab="log Chla")
dev.off()

## Figure 10.21
tikz(file=paste(plotDIRch10, "Finncoplottype62.tex", sep="/"),
           width=4.5, height=4, standAlone=F)
par(mar=c(3,3,3,.25), mgp=c(1.5,0.125,0), tck=0.01)
given.tp <- co.intervals(lake4$lxp, number=4, overlap=.1)
 coplot(y ~ lxn | lxp, data = lake4, given.v=given.tp, rows=1,
        xlab="log TN", ylab="log Chla")
dev.off()

plot(lxn~lxp, data=lake4)

## Figure 10.22
lake4 <- summer.All[summer.All$type==7,]
tikz(file=paste(plotDIRch10, "Finncoplottype71.tex", sep="/"),
           width=4.5, height=4, standAlone=F)
par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), tck=0.01)
given.tn <- co.intervals(lake4$lxn, number=4, overlap=.1)
 coplot(y ~ lxp | lxn, data = lake4, given.v=given.tn, rows=1,
        xlab="log TP", ylab="log Chla")
dev.off()

## Figure 10.23
tikz(file=paste(plotDIRch10, "Finncoplottype72.tex", sep="/"),
           width=4.5, height=4, standAlone=F)
par(mar=c(3,3,3,.25), mgp=c(1.25,0.125,0), tck=0.01)
given.tp <- co.intervals(lake4$lxp, number=4, overlap=.1)
 coplot(y ~ lxn | lxp, data = lake4, given.v=given.tp, rows=1,
        xlab="log TN", ylab="log Chla")
dev.off()

## nonlinear lmer (nlmer)
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

## See Chapter 6 script for SSfpl2
tm1 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==1,])

tm2 <- nls(Abs ~ SSfpl(log(stdConc), A, B, xmid, scal),
           data=toledo1[toledo1$Test==1,])


tm1.nlmer <- nlmer(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4)~
                       (al1+al2+al3+al4|Test), data=toledo,
                   start=c(al1=1, al2=1, al3=0.5, al4=0.2))
print(tm1.nlmer)

tm2.nlmer <- nlmer(Abs ~ SSfpl(log(stdConc), A, B, xmid, scal) ~
                               (A+B+xmid+scal|Test), data=toledo1,
                               start=c(A=1,B=0.2, xmid=-0.6, scal=0.75))
print(tm2.nlmer)

tikz(file=paste(plotDIRch10, "elisaTOranef1.tex", sep="/"),
     height=3.5, width=5.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
temp <- ranef(tm1.nlmer, condVar=T)
names(temp$Test) <- c("$\\alpha_1$", "$\\alpha_2$", "$\\alpha_3$",
                      "$\\alpha_4$")
dotplot(temp, layout=c(4,1), main=F)
dev.off()

tikz(file=paste(plotDIRch10, "elisaTOranef2.tex", sep="/"),
     height=3.5, width=5.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
temp <- ranef(tm2.nlmer, condVar=T)
names(temp$Test) <- c("$A$", "$B$", "$x_{mid}$",
                      "$scal$")
dotplot(temp, layout=c(4,1), main=F)
dev.off()


TM1 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==1,])
TM2 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==2,])
TM3 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==3,])
TM4 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==4,])
TM5 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==5,])
TM6 <- nls(Abs ~ SSfpl2(stdConc, al1, al2, al3, al4),
           data=toledo[toledo$Test==6,])

oldTM2 <- nls(Abs ~ (A-D)/(1+(stdConc/C)^B)+D,
           control=list(maxiter=200), data=toledo[toledo$Test==2,],
           start=list(A=0.2,B=1.,C=0.45,D=1.))
oldTM3 <- nls(Abs ~ (A-D)/(1+(stdConc/C)^B)+D,
           control=list(maxiter=200), data=toledo[toledo$Test==3,],
           start=list(A=0.16,B=1.12,C=0.45,D=1.06))
oldTM4 <- nls(Abs ~ (A-D)/(1+(stdConc/C)^B)+D,
           control=list(maxiter=200), data=toledo[toledo$Test==4,],
           start=list(A=0.16,B=1.12,C=0.45,D=1.06))
oldTM5 <- nls(Abs ~ (A-D)/(1+(stdConc/C)^B)+D,
           control=list(maxiter=200), data=toledo[toledo$Test==5,],
           start=list(A=0.16,B=1.12,C=0.45,D=1.06))
oldTM6 <- nls(Abs ~ (A-D)/(1+(stdConc/C)^B)+D,
           control=list(maxiter=200), data=toledo[toledo$Test==6,],
           start=list(A=0.16,B=1.12,C=0.45,D=1.06))

## difference in logit scale is more or less the same.


## Galax
## 2010 data
survey.data <- read.csv(paste(dataDIR, "galax_2010.csv", sep="/"),
                        header=T)

## check for unique species:

levels(survey.data$Site)
 table(survey.data$Site)

## model 1: Large leaf counts -- PLOT effect
lmer1 <- glmer(CountL~1+(1|Site), family="poisson",
               offset=log(TotalPoints),  ## using offset to model density
               data=survey.data)
display(lmer1)

tikz(file=paste(plotDIRch10, "galaxRanef2010.tex", sep="/"),
     height=3.25, width=4, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.75 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=1 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.55
trellis.par.set(trellis.par.temp)
dotplot(ranef(lmer1, condVar=T))
dev.off()


new.data <- data.frame(x = fixef(lmer1)[1] + ranef(lmer1)$Site[,1],
                       y=row.names(ranef(lmer1)$Site),
                       sd = sqrt(se.ranef(lmer1)$Site[,1]^2+
                                     se.fixef(lmer1)[1]^2))
new.data$y <- ordered(new.data$y, levels=new.data$y[order(new.data$x)])

tikz(file=paste(plotDIRch10, "galaxleafLD.tex", sep="/"), height=3.5,
               width=5, standAlone=F)
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.75 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=1 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
trellis.par.set(trellis.par.temp)
dotplot(y ~ exp(x), data = new.data,
        aspect = 0.8,
        xlim=c(0,1.2*range(exp(new.data$x-new.data$sd),
          exp(new.data$x+new.data$sd))[2]),
        panel = function (x, y) {
          panel.xyplot(x, y, pch = 16, col = "black")
          panel.segments(exp(new.data$x-new.data$sd), as.numeric(y),
                         exp(new.data$x+new.data$sd), as.numeric(y),
                         lty = 1, col = "black")},
        xlab="Density")
dev.off()

## Model 2 -- large leaf ratio
lmer3 <- glmer(cbind(CountL, CountS)~1+(1|Site),
               family="binomial",data=survey.data)
display(lmer3)

tikz(file=paste(plotDIRch10, "galaxlogitranef.tex", sep="/"),
     height=3.25, width=4, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.75 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=1 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.55
trellis.par.set(trellis.par.temp)
dotplot(ranef(lmer3, condVar=T))
dev.off()


lmer3.sim <- sim(lmer3, 2000)
## predicted ratio:
pred3 <- rvsims(lmer3.sim@fixef[,1] + lmer3.sim@ranef$Site[,,1])
pred3.sum <- summary(invlogit(pred3))

new.data <- data.frame(x = pred3.sum[,2],
                       y = labels(pred3),
                       q2.5=pred3.sum[,5],
                       q25= pred3.sum[,6],
                       q50= pred3.sum[,7],
                       q75 =pred3.sum[,8],
                       q975=pred3.sum[,9])
new.data$y <- ordered(new.data$y, levels=new.data$y[order(new.data$x)])
tikz(file=paste(plotDIRch10, "leafFR.tex", sep="/"),
     height=3.5, width=5, standAlone=F)
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.75 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=1 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
trellis.par.set(trellis.par.temp)
dotplot(y~x, data=new.data,xlim=c(0,1),
       # xlim=c(0.8,1.2)*range(new.data$q2.5, new.data$q975),
        aspect=0.8,
        panel=function(x,y){
##          panel.xyplot(x, y, pch=16, col="black")
          panel.points(x=new.data$q50, y=as.numeric(y), pch=16)
          panel.segments(new.data$q2.5, as.numeric(y),
                         new.data$q975, as.numeric(y), lty=1)
          panel.segments(new.data$q25, as.numeric(y),
                         new.data$q75, as.numeric(y), lty=1, lwd=2)
        }, xlab="large leaf fraction")
dev.off()

## Crypto in drinking water
icr <- read.table(paste(dataDIR,"CryptoData2.txt",sep="/"), header=T,
        na.string = ".", sep = "\t")

xy <- data.frame(x=stdConc8.1, y=Abs8.1.0)

rm(xy)
## example of Poisson multilevel model

levels(icr$MSrcCat) -> Lmsrc
levels(icr$M.WTP.Type)->Lwtptype

dcts.data <- icr[icr$MSrcCat==Lmsrc[3] | icr$MSrcCat==Lmsrc[4],]
dcts.data <- dcts.data[dcts.data$M.WTP.Type=="Y"|dcts.data$M.WTP.Type=="N",]
dcts.data$M.WTP.Type <- ordered(as.vector( dcts.data$M.WTP.Type ))
dcts.data$MSrcCat <- ordered(as.vector( dcts.data$MSrcCat ))
dcts.data$ICR.PWSID <- ordered(dcts.data$ICR.PWSID)

icr.glm <- glm(n.cT ~ factor(ICR.PWSID)-1, data=dcts.data,
               family="poisson", offset=log(volume*0.44))

display(icr.glm)

icr.lmer1 <- glmer(n.cT ~ 1+(1|ICR.PWSID),
                   data=dcts.data, family="poisson",
                   offset=log(volume*0.44))

icr.size <- as.vector(table(dcts.data$ICR.PWSID))
size <- icr.size + runif(length(icr.size), -0.1,0.1)

icr.glmCoef <- summary(icr.glm)$coef

## Figure 10.29
##postscript(file=paste(plotDIR, "cryptolmer1.eps", sep="/"),
##           width=3.5, height=2.5, horizontal=F)
tikz(file=paste(plotDIRch10, "cryptolmer1.tex", sep="/"),
     width=3.5, height=2.5, standAlone=F)
lmer.mean <- fixef(icr.lmer1)[1] + ranef(icr.lmer1)[[1]][,1]
lmer.se <- sqrt(se.fixef(icr.lmer1)[1]^2 + se.ranef(icr.lmer1)[[1]][,1]^2)
lower <- lmer.mean-lmer.se
upper <- lmer.mean+lmer.se

par(mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), tck=0.01, las=1)
plot(size, rnorm(length(size)),type="n", xlab="Sample size",
  ylab="Log Mean", log="x", ylim=range(lower,upper))
abline(h=fixef(icr.lmer1)[[1]][1])
segments(x0=size, x1=size, y0=lower, y1=upper)
points(size, lmer.mean, pch=16,cex=0.5)
dev.off()

##icr.lmer2 <- glmer(n.cT ~ 1+(1|ICR.PWSID),
##                   data=dcts.data, family="quasipoisson",
##                   offset=log(volume*0.44))

cs <- coef(icr.lmer1)[[1]][,1]
n.cs <- length(cs)

## Figure 10.30
##postscript(file=paste(plotDIR, "cryptocdf.eps", sep="/"),
##           width=3.5, height=2.5, horizontal=F)
tikz(file=paste(plotDIRch10, "cryptocdf.tex", sep="/"),
     width=3.5, height=2.5, standAlone=F)
par(mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(sort(cs), ((1:n.cs)-0.5)/n.cs, axes=F,
     xlab="Crypto Concentration (oocyst/L)", ylab="CDF",
     type="n")
points(sort(cs), ((1:n.cs)-0.5)/n.cs, pch=1,cex=0.5, col="gray")
curve(pnorm(x, -5.384, 2.08), add=T)
axis(1, at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),
     labels=c("0.0001","0.001","0.01","0.1","1","10","100"))
axis(2)
box()
dev.off()

## simulation fraction of 0s (at system level)
n.sys <- 884
n.sims <- 10000
zeros <- sys.mean <- sys.25<-sys.75<-sys.95 <-
    sys.99 <- numeric()
sys.means <- matrix(0, n.sims, n.sys)
for (i in 1:n.sims){
    zeros[i] <- 0
    for (j in 1:n.sys){
        mu <- rnorm(1, -5.384, 0.103)
        sigma <- 2.08*sqrt((13103-884)/rchisq(1, 13103-884))
        y <- rpois(icr.size[j],
                   0.44*10*exp(rnorm(icr.size[j], mu, sigma)))
        sys.means[i,j] <- mean(y)/10
        zeros[i] <- zeros[i] + (sum(y!=0)==0)/n.sys
    }
}
dcts.sys.mean <- tapply(dcts.data$n.cT/dcts.data$volume,
                        dcts.data$ICR.PWSID,
                        mean)

## postscript(file=paste(plotDIR, "cryptoSim.eps", sep="/"),
##           height=2.5, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch10, "cryptoSim.tex", sep="/"),
     height=2.5, width=4.5, standAlone=F)
par(mfrow=c(1,2), mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0), tck=0.01)
hist(zeros*100, xlab="\\% of all zero systems", main="",
     xlim=c(40,70), cex=0.75)
abline(v=mean(dcts.sys.mean==0)*100, col="gray", lwd=2)
## simulation of the 99th percentile of the system mean

hist(as.vector(apply(sys.means, 1, quantile, prob=0.99)),
     main="", xlab="99th percentile", cex=0.75)
abline(v=quantile(dcts.sys.mean, 0.99), col="gray", lwd=2)
dev.off()

