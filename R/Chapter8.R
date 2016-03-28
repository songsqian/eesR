source("FrontMatter.R")
packages(rv)

## Chapter 8
plotDIRch8 <- paste(plotDIR, "chapter8", "figures", sep="/")
## GLM

### Crypto data: data combined from Finch et at 93 and Korich et al 2000
#
crypto.data <- read.table(paste(dataDIR, "cryptoDATA.csv", sep="/"), header=T) ## but not comma separated
crypto.data$Y <- round(crypto.data$Y)


crypto.lm1 <- lm(I(logit(Y/N))~log10(Dose), data=crypto.data, subset=Y/N!=0 & Y/N!=1)
display(crypto.lm1)

crypto.glm1 <- glm(cbind(Y, N-Y)~log10(Dose), data=crypto.data, family=binomial(link="logit"))
display(crypto.glm1, 3)
## overdispersion
z <- (crypto.data$Y-crypto.data$N*fitted(crypto.glm1))/sqrt(crypto.data$N*fitted(crypto.glm1)*(1-fitted(crypto.glm1)))
z.chisq <- sum(z^2)
orverD <- z.chisq/summary(crypto.glm1)$df[2]
p.value <- 1-pchisq(z.chisq, df=summary(crypto.glm1)$df[2])



### creating new data to convert N Y to 0 and 1
temp <- NULL
for (i in 1:dim(crypto.data)[1]){
    temp <- c(temp, rep(c(0, 1), c(crypto.data$N[i]-crypto.data$Y[i], crypto.data$Y[i])))
}
new.crypto <- data.frame (Dose = rep(crypto.data$Dose, crypto.data$N), Y = temp)

## Figure 8.1
## postscript(file=paste(plotDIR, "cryptoGLM1.eps", sep="/"), height=4.5, width=6, horizontal=F)
tikz(file=paste(plotDIRch8, "cryptoGLM1.tex", sep="/"),
     height=4.5, width=6, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25, 0.125,0), las=1, tck=0.01)
plot(jitter.binary(Y)~jitter(Dose), data=new.crypto, log="x", xlab="Oocyst Dose", col="gray10", cex=0.5, ylab="Prob. of Infection")
curve(invlogit(coef(crypto.glm1)[1] + coef(crypto.glm1)[2]*log10(x)), add=T)
points(Y/N~Dose, data=crypto.data, col="gray", pch=16)
dev.off()

## Figure 8.2

##postscript(file=paste(plotDIR, "logitFIG.eps", sep="/"),
##           height=3, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch8, "logitFIG.tex", sep="/"),
           height=3, width=4.5, standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.25,.125,0), las=1, tck=0.01)
curve(invlogit(x), xlim=c(-4,4), ylab="$logit^{-1}(x)$")
dev.off()


## Figure 8.3
##postscript(file=paste(plotDIR, "cryptodata.eps", sep="/"), width=6, height=3.5, horizontal=F)
tikz(file=paste(plotDIRch8, "cryptodata.tex", sep="/"),
     width=6, height=3.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
 xyplot(logit(Y/N)~log10(Dose)|Source, data=crypto.data,
    panel=function(x,y,...){
        panel.xyplot(x,y,...)
        panel.lmline(x,y,col="gray",...)
        panel.grid()
        }, aspect=1, xlab="Log Oocyst Dose", ylab="Logit(f)",
    layout=c(4,1))
dev.off()


##postscript(file=paste(plotDIR, "cryptoGLM2.eps", sep="/"), height=4.5, width=6, horizontal=F)
tikz(file=paste(plotDIRch8, "cryptoGLM2.tex", sep="/"),
     height=4.5, width=6, standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.25,.125,0), las=1, tck=0.01)
plot(jitter.binary(Y)~jitter(Dose), data=new.crypto, log="x",
     xlab="Oocyst Dose", col="gray10", cex=0.5,
     ylab="Prob. of Infection")
curve(invlogit(coef(crypto.glm1)[1] + coef(crypto.glm1)[2]*log10(x)), add=T, col="blue")
curve(invlogit(coef(crypto.lm1)[1] + coef(crypto.lm1)[2]*log10(x)), add=T, col="red")
points(Y/N~Dose, data=crypto.data, col="gray", pch=16)
dev.off()

##
crypto.glm2 <- glm(cbind(Y, N-Y)~log10(Dose)+factor(Source)-1,
                   data=crypto.data, family=binomial(link="logit"))
display(crypto.glm2)
crypto.glm3 <- glm(cbind(Y, N-Y)~log10(Dose)*factor(Source)-1-log10(Dose),
                   data=crypto.data, family=binomial(link="logit"))
display(crypto.glm3)
crypto.glm4 <- glm(cbind(Y, N-Y)~log10(Dose)*factor(Source), data=crypto.data,
                   family=binomial(link="logit"))
display(crypto.glm4)

crypto.lm2 <- lm(I(logit(Y/N))~log10(Dose)+factor(Source)-1,
                 data=crypto.data, subset=Y/N!=0 & Y/N!=1)
display(crypto.lm2)
crypto.lm3 <- lm(I(logit(Y/N))~log10(Dose)*factor(Source)-1-log10(Dose),
                 data=crypto.data, subset=Y/N!=0 & Y/N!=1)
display(crypto.lm3)
crypto.lm4 <- lm(I(logit(Y/N))~log10(Dose)*factor(Source), data=crypto.data,
                 subset=Y!=0 & Y!=N)
display(crypto.lm4)


pred.4 <- predict(crypto.glm4, newdata=data.frame(Dose = rep(crypto.data$Dose, crypto.data$N),
                                   Source=rep(crypto.data$Source, crypto.data$N)),
                  type="response")

## Figure 8.4
##postscript (paste(plotDIR, "cryptoRes.eps", sep="/"), height=3.5, width=4.5, horizontal=F)
tikz (paste(plotDIRch8, "cryptoRes.tex", sep="/"),
      height=3.5, width=4.5, standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.25,.125,0), las=1, tck=0.01)
plot(c(0,1), c(-1,1), xlab="Estimated  Pr (Infection)", ylab="Observed - estimated",
     type="n", main="Residual plot", mgp=c(2,.5,0))
abline (0,0, col="gray", lwd=.5)
points (pred.4, new.crypto$Y-pred.4, pch=20, cex=.2)
dev.off ()

## Figure 8.5
##postscript (paste(plotDIR, "cryptoResB.eps", sep="/"), height=3.5, width=6.5, horizontal=F)
tikz (paste(plotDIRch8, "cryptoResB.tex", sep="/"),
      height=3.5, width=6.5, standAlone=F)
par(mfrow=c(1,2), mar=c(3,3,1.5,0.5), mgp=c(1.25,.125,0), las=1, tck=0.01)
plot(c(0,1), c(-1,1), xlab="Estimated  Pr (Infection)",
     ylab="Observed - estimated", type="n", main="Residual plot")
abline (0,0, col="gray", lwd=.5)
points (pred.4, new.crypto$Y-pred.4, pch=20, cex=.2)

br.8 <- binned.resids (pred.4, new.crypto$Y-pred.4, nclass=30)$binned
plot(range(br.8[,1]), range(br.8[,2],br.8[,6],-br.8[,6]),
     xlab="Estimated  Pr (Infection)", ylab="Average residual", type="n",
     main="Binned residual plot")
abline (0,0, col="gray", lwd=.5)
lines (br.8[,1], br.8[,6], col="gray", lwd=.5)
lines (br.8[,1], -br.8[,6], col="gray", lwd=.5)
points (br.8[,1], br.8[,2], pch=20, cex=.5)
dev.off()

## Another Example of binary data

## Another Example of binary data

# data
seedbank <- read.csv(paste(dataDIR, "seedbank.csv", sep="/"),
                     header=T, na.string="-")

attach(seedbank)
 plot(jitter.binary(Predation)~as.numeric(species), data=seedbank)

seedbank.glm1 <- glm(Predation ~ log(seed.weight), data=seedbank,
                     family=binomial(link="logit"))
M1.coef<-coef(seedbank.glm1)

## Figure 8.6
##postscript(file="seeddata.eps", height=4, width=5.5, horizontal=F)
tikz(file=paste(plotDIRch8, "seeddata.tex", sep="/"),
     height=4, width=5.5, standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
plot(jitter.binary(Predation)~log(seed.weight), data=seedbank,
     col="gray", ylab="Predation", xlab="Seed Weight (g/1000 seeds)",
     cex=0.5, axes=F)
axis(2, las=1)
axis(1, at=log(c(1, 5, 10, 50, 100, 500, 1000, 5000)),
     label=as.character(c(1, 5, 10, 50, 100, 500, 1000, 5000)))
box()
points(log(sort(unique(seedbank$seed.weight))),
       as.vector(by(seedbank$Predation, factor(seedbank$seed.weight),
                    mean, na.rm=T)), cex=1.5)
curve(invlogit(M1.coef[1]+M1.coef[2]*x), add=T)
dev.off()

seedbank$species <- ordered(seedbank$species, levels=rev(1:8))
seedbank$logW.c <- scale(log(seedbank$seed.weight),scale= F)
names(seedbank)

#[1] "species"     "seed.weight" "time"        "topo"        "ground"
#[6] "Predation"
#>
#seedbank.glm <- glm(Predation ~ factor(species)+factor(time)+factor(topo)+factor(ground), data=seedbank, family=binomial(link="logit"))
#summary(seedbank.glm)

plot(seed.weight~species, data=seedbank)

seedbank.glm1 <- glm(Predation ~ log(seed.weight),
                     data=seedbank, family=binomial(link="logit"))
display(seedbank.glm1)

seedbank.glm2 <- glm(Predation ~ log(seed.weight)+factor(species)-1,
                     data=seedbank, family=binomial(link="logit"))
display(seedbank.glm2)

seedbank.glm2B1 <- bayesglm(Predation ~ log(seed.weight)+factor(species)-1,
                            data=seedbank, family=binomial(link="logit"))
display(seedbank.glm2B1)

seedbank.glm2 <- glm(Predation ~ species-1, data=seedbank, family=binomial(link="logit"))
display(seedbank.glm2)

seedbank.glm2B <- bayesglm(Predation ~ species-1, data=seedbank, family=binomial(link="logit"))
display(seedbank.glm2B)

seedbank.glm2.1 <- glm(Predation ~ species, data=seedbank, family=binomial(link="logit"))
display(seedbank.glm2.1)

seedbank.glm3 <- glm(Predation ~ factor(time)+factor(species), data=seedbank, family=binomial(link="logit"))
display(seedbank.glm3)

seedbank.glm4 <- glm(Predation ~ factor(time)+log(seed.weight), data=seedbank, family=binomial(link="logit"))
display(seedbank.glm4)

seedbank.glm4.2 <- glm(Predation ~ factor(time)+log(seed.weight)-1,
                       data=seedbank, family=binomial(link="logit"))
display(seedbank.glm4.2)
invlogit(coef(seedbank.glm4.2)[1:6] + coef(seedbank.glm4.2)[7] * log(28.5))

x.range <- range(c(summary(seedbank.glm4.2)$coef[1:6,1]+2*summary(seedbank.glm4.2)$coef[1:6,2],
                   summary(seedbank.glm4.2)$coef[1:6,1]-2*summary(seedbank.glm4.2)$coef[1:6,2] ))
plot(x.range, c(1,6), type="n", xlab="Intercept", ylab="Time")
points(summary(seedbank.glm4.2)$coef[1:6,1], 1:6, col="red")
segments(x0=summary(seedbank.glm4.2)$coef[1:6,1]+2*summary(seedbank.glm4.2)$coef[1:6,2], y0=1:6,
                   x1=summary(seedbank.glm4.2)$coef[1:6,1]-2*summary(seedbank.glm4.2)$coef[1:6,2],y1=1:6, col="blue" )

## Figure 8.7
##postscript(file=paste(plotDIR, "seedbankPlot.eps", sep="/"),
##           width=4, height=2, horizontal=F)  ## for book
tikz(file=paste(plotDIRch8, "seedbankPlot.tex", sep="/"),
           width=4, height=2, standAlone=F)  ## for book
par(mar=c(3,3,0.25,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
x.range <- range(c(summary(seedbank.glm4.2)$coef[1:6,1]+
                   2*summary(seedbank.glm4.2)$coef[1:6,2],
                   summary(seedbank.glm4.2)$coef[1:6,1]-
                   2*summary(seedbank.glm4.2)$coef[1:6,2]))
plot(x.range, c(0.75,6.25), type="n",
     xlab="Intercept", ylab="Time", mgp=c(1.5,.5,0), cex=0.75)
points(summary(seedbank.glm4.2)$coef[1:6,1], 1:6, cex=1.5)
segments(x0=summary(seedbank.glm4.2)$coef[1:6,1]+2*summary(seedbank.glm4.2)$coef[1:6,2],
         y0=1:6,
         x1=summary(seedbank.glm4.2)$coef[1:6,1]-2*summary(seedbank.glm4.2)$coef[1:6,2],
         y1=1:6, col="gray10" )
segments(x0=summary(seedbank.glm4.2)$coef[1:6,1]+summary(seedbank.glm4.2)$coef[1:6,2],
         y0=1:6,
         x1=summary(seedbank.glm4.2)$coef[1:6,1]-summary(seedbank.glm4.2)$coef[1:6,2],
         y1=1:6, lwd=3 )
dev.off()



### simulation
betas <- sim(seedbank.glm4.2,5000)@coef

see.weight <- sort(unique(seedbank$seed.weight))

## Figure 8.9
##postscript(file=paste(plotDIR, "seedbankSim.eps", sep="/"),
##           width=5.5, height=7, horizontal=F)  ## for book
tikz(file=paste(plotDIRch8, "seedbankSim.tex", sep="/"),
           width=5.5, height=7, standAlone=F)  ## for book
par(mfrow=c(4,2), mar=c(3,3,0.5,1), las=1, tck=0.01, mgp=c(1.25,0.125, 0))
for (i in 1:8){
  coef.sim<-t(apply(invlogit(betas[, 1:6] + betas[,7] * log(see.weight[i])), 2,
                    FUN=function(x){
                      return(c(mean(x), quantile(x, prob=c(0.025,0.25,0.75, 0.975))))}))
  x.range <- range(coef.sim)
  plot(x.range, c(0.75,6.25), xlim=c(0,0.8),type="n",
       xlab="Probability", ylab="Time", mgp=c(1.5,.5,0), cex=0.75)
  points(coef.sim[1:6,1], 1:6, cex=1)
  segments(x0=coef.sim[1:6,2], y0=1:6, x1=coef.sim[1:6,5], y1=1:6, col="gray10" )
  segments(x0=coef.sim[1:6,3], y0=1:6, x1=coef.sim[1:6,4], y1=1:6, lwd=3)
  if(i <8)    text(0.6, 6, paste("Seed Weight = ",see.weight[i]))
  else text(0.2, 6, paste("Seed Weight = ",see.weight[i]))
}
dev.off()


seedbank.glm4.1 <- glm(Predation ~ factor(time)*log(seed.weight),
                       data=seedbank, family=binomial(link="logit"))
display(seedbank.glm4.1)

##postscript(paste(plotDIR, "seedbankIDprob.eps", sep="/"), height=4, width=6, horizontal=F)
tikz(paste(plotDIRch8, "seedbankIDprob.tex", sep="/"),
     height=4, width=6, standAlone=F)
xyplot(jitter.binary(Predation)~log10(seed.weight)|factor(time),
       data=seedbank, xlab="log seed weight", ylab="Predation",aspect=1)
dev.off()

seedbank.glm5 <- glm(Predation ~ factor(time)+factor(topo)+log(seed.weight),
                     data=seedbank, family=binomial(link="logit"))
display(seedbank.glm5)

seedbank.glm6 <- glm(Predation ~ factor(time)+factor(topo)*log(seed.weight),
                     data=seedbank, family=binomial(link="logit"))
display(seedbank.glm6)

seedbank.glm7 <- glm(Predation ~ factor(time)+factor(topo)*log(seed.weight)+factor(ground),
                     data=seedbank, family=binomial(link="logit"))
display(seedbank.glm7)

## overdispersion ##
z <- (na.omit(seedbank$Predation)-pred.7)/sqrt(pred.7*(1-pred.7))
sum.z <- sum(z^2)
p.value <- 1-pchisq(sum.z, length(pred.7)-length(coef(seedbank.glm7)))
est.overD <- sum.z/(length(pred.7)-length(coef(seedbank.glm7)))

#### Plotting the fitted models

seedbank.glm1 <- glm(Predation ~ log(seed.weight), data=seedbank, family=binomial(link="logit"))
display(seedbank.glm1)

##postscript(file="seedbankM1.eps", width=4, height=3, horizontal=F)
tikz(file=paste(plotDIRch8, "seedbankM1.tex", sep="/"),
     width=4, height=3, standAlone=F)
betas <- coef(seedbank.glm1)
plot(jitter.binary(Predation) ~ log(seed.weight), type="n", data=seedbank, xlab="log seed weight", ylab="prob. of predation")
points(jitter.binary(Predation) ~ log(seed.weight), col="gray")
curve(invlogit(betas[1]+betas[2]*x), add=T)
dev.off()

seedbank.glm4 <- glm(Predation ~ factor(time)+log(seed.weight), data=seedbank, family=binomial(link="logit"))
display(seedbank.glm4)

#(Intercept) factor(time)2 factor(time)3 factor(time)4 factor(time)5
#     -4.81590       2.37404       2.70389       3.11282       2.98472
#factor(time)6        logW.c
#      3.09617       0.45543

## Figure 8.8
##postscript(file=paste(plotDIR,"seedbankM4.eps",sep="/"),
##           width=4, height=3, horizontal=F)
tikz(file=paste(plotDIRch8,"seedbankM4.tex",sep="/"),
           width=4, height=3, standAlone=F)
betas <- coef(seedbank.glm4)
par(mgp=c(1.25,.125,0), mar=c(3,3,1,0.25), las=1, tck=0.01)
plot(jitter.binary(Predation) ~ log(seed.weight), type="n", data=seedbank,
     xlab="log seed weight", ylab="prob. of predation")
points(jitter.binary(Predation) ~ log(seed.weight), col="gray")
curve(invlogit(betas[1]+betas[7]*x), add=T, col=gray(.1))
curve(invlogit(betas[1]+betas[2]+betas[7]*x), add=T, lty=2, col=gray(.2))
curve(invlogit(betas[1]+betas[3]+betas[7]*x), add=T, lty=3, col=gray(.3))
curve(invlogit(betas[1]+betas[4]+betas[7]*x), add=T, lty=4, col=gray(.4))
curve(invlogit(betas[1]+betas[5]+betas[7]*x), add=T, lty=5, col=gray(.5))
curve(invlogit(betas[1]+betas[6]+betas[7]*x), add=T, lty=6, col=gray(.6))
legend(x=0, y=0.9, legend=month.name[seq(1,11,2)], lty=1:6, col=gray((1:6)/10),
       cex=0.5, bty="n")
dev.off()

seedbank.glm5 <- glm(Predation ~ factor(time)+factor(topo)+log(seed.weight),
                     data=seedbank, family=binomial(link="logit"))
display(seedbank.glm5)

#  (Intercept) factor(time)2 factor(time)3 factor(time)4 factor(time)5
#     -3.93925       2.68160       3.05653       3.52552       3.37838
#factor(time)6 factor(topo)2 factor(topo)3 factor(topo)4        logW.c
#      3.58596      -2.21355      -1.78835      -2.27402       0.53575

topog <- c("Hilltop","Shady Slope","Sunny Slope","Valley")

## Figure 8.10
##postscript(file="seedbankM5.eps", width=5, height=5, horizontal=F)
tikz(file=paste(plotDIRch8, "seedbankM5.tex", sep="/"),
     width=5, height=5, standAlone=F)
betas <- coef(seedbank.glm5)
par(mfrow=c(2,2),mgp=c(1.25, 0.125,0), las=1, tck=0.01, mar=c(3,3,3,1))
plot(jitter.binary(Predation) ~ log(seed.weight), type="n", data=seedbank, xlab="log seed weight", ylab="prob. of predation")
points(jitter.binary(Predation) ~ log(seed.weight), col="gray", subset=topo==1)
curve(invlogit(betas[1]+betas[10]*x), add=T, col=gray(.1))
curve(invlogit(betas[1]+betas[2]+betas[10]*x), add=T, lty=2, col=gray(.2))
curve(invlogit(betas[1]+betas[3]+betas[10]*x), add=T, lty=3, col=gray(.3))
curve(invlogit(betas[1]+betas[4]+betas[10]*x), add=T, lty=4, col=gray(.4))
curve(invlogit(betas[1]+betas[5]+betas[10]*x), add=T, lty=5, col=gray(.5))
curve(invlogit(betas[1]+betas[6]+betas[10]*x), add=T, lty=6, col=gray(.6))
legend(x=0, y=0.9, legend=month.name[seq(1,11,2)], lty=1:6, col=gray((1:6)/10), cex=0.5, bty="n")
title(main=topog[1], cex=0.75)
for (i in c(3, 2, 4)){
    plot(jitter.binary(Predation) ~ log(seed.weight), type="n", data=seedbank,
        xlab="centered log seed weight", ylab="prob. of predation")
    points(jitter.binary(Predation) ~ log(seed.weight), col="gray", subset=topo==i)
    curve(invlogit(betas[1]+betas[10]*x), add=T, col=gray(.1))
    curve(invlogit(betas[1]+betas[2]+betas[i+5]+betas[10]*x), add=T, lty=2, col=gray(.2))
    curve(invlogit(betas[1]+betas[3]+betas[i+5]+betas[10]*x), add=T, lty=3, col=gray(.3))
    curve(invlogit(betas[1]+betas[4]+betas[i+5]+betas[10]*x), add=T, lty=4, col=gray(.4))
    curve(invlogit(betas[1]+betas[5]+betas[i+5]+betas[10]*x), add=T, lty=5, col=gray(.5))
    curve(invlogit(betas[1]+betas[6]+betas[i+5]+betas[10]*x), add=T, lty=6, col=gray(.6))
    title(main=topog[i], cex=0.75)
}
dev.off()

seedbank.glm6 <- glm(Predation ~ factor(time)+factor(topo)*log(seed.weight),
                     data=seedbank, family=binomial(link="logit"))
display(seedbank.glm6)

#> betas
#         (Intercept)        factor(time)2        factor(time)3
#            -4.52327              3.20512              3.59974
#       factor(time)4        factor(time)5        factor(time)6
#             4.07728              3.92898              4.13750
#       factor(topo)2        factor(topo)3        factor(topo)4
#            -1.99616             -1.61685             -2.06389
#              logW.c factor(topo)2:logW.c factor(topo)3:logW.c
#             0.76652             -0.30976             -0.29951
#factor(topo)4:logW.c
#            -0.30458

topog <- c("Hilltop","Shady Slope","Sunny Slope","Valley")

## Figure 8.11
##postscript(file="seedbankM6.eps", width=5, height=5, horizontal=F)
tikz(file=paste(plotDIRch8, "seedbankM6.tex", sep="/"),
     width=5, height=5, standAlone=F)
betas <- coef(seedbank.glm6)
par(mfrow=c(2,2),mgp=c(1.25, 0.125,0), mar=c(3,3,3,1), las=1, tck=0.01)
plot(jitter.binary(Predation) ~ log(seed.weight), type="n", data=seedbank, xlab="centered log seed weight", ylab="prob. of predation")
points(jitter.binary(Predation) ~ log(seed.weight), col="gray", subset=topo==1)
curve(invlogit(betas[1]+betas[10]*x), add=T, col=gray(.1))
curve(invlogit(betas[1]+betas[2]+betas[10]*x), add=T, lty=2, col=gray(.2))
curve(invlogit(betas[1]+betas[3]+betas[10]*x), add=T, lty=3, col=gray(.3))
curve(invlogit(betas[1]+betas[4]+betas[10]*x), add=T, lty=4, col=gray(.4))
curve(invlogit(betas[1]+betas[5]+betas[10]*x), add=T, lty=5, col=gray(.5))
curve(invlogit(betas[1]+betas[6]+betas[10]*x), add=T, lty=6, col=gray(.6))
legend(x=0, y=0.9, legend=month.name[seq(1, 11, 2)], lty=1:6, col=gray((1:6)/10), cex=0.5, bty="n")
title(main=topog[1], cex=0.75)
for (i in c(3, 2, 4)){
    plot(jitter.binary(Predation) ~ log(seed.weight), type="n", data=seedbank,
        xlab="centered log seed weight", ylab="prob. of predation")
    points(jitter.binary(Predation) ~ log(seed.weight), col="gray", subset=topo==i)
    curve(invlogit(betas[1]+betas[10]*x), add=T, col=gray(.1))
    curve(invlogit(betas[1]+betas[2]+betas[i+5]+(betas[10]+betas[i+9])*x), add=T, lty=2, col=gray(.2))
    curve(invlogit(betas[1]+betas[3]+betas[i+5]+(betas[10]+betas[i+9])*x), add=T, lty=3, col=gray(.3))
    curve(invlogit(betas[1]+betas[4]+betas[i+5]+(betas[10]+betas[i+9])*x), add=T, lty=4, col=gray(.4))
    curve(invlogit(betas[1]+betas[5]+betas[i+5]+(betas[10]+betas[i+9])*x), add=T, lty=5, col=gray(.5))
    curve(invlogit(betas[1]+betas[6]+betas[i+5]+(betas[10]+betas[i+9])*x), add=T, lty=6, col=gray(.6))
    title(main=topog[i], cex=0.75)
}
dev.off()


pred.7 <- predict(seedbank.glm7, type="response")

## Figure 8.12
##postscript(paste(plotDIR, "seedResB.eps", sep="/"), height=4, width=5.5, horizontal=F)
tikz(paste(plotDIRch8, "seedResB.tex", sep="/"),
     height=4, width=5.5, standAlone=F)
par(mar=c(3,3,3,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
br.8 <- binned.resids (pred.7, na.omit(seedbank$Predation)-pred.7, nclass=20)$binned
plot(range(br.8[,1]), range(br.8[,2],br.8[,6],-br.8[,6]), xlab="Estimated  Pr (Predation)", ylab="Average residual", type="n", main="Binned residual plot", mgp=c(2,.5,0))
abline (0,0, col="gray", lwd=.5)
lines (br.8[,1], br.8[,6], col="gray", lwd=.5)
lines (br.8[,1], -br.8[,6], col="gray", lwd=.5)
points (br.8[,1], br.8[,2], pch=20, cex=.5)
dev.off()

error.rate <- mean((pred.7>0.5 & na.omit(seedbank$Predation)==0) |
                   (pred.7<0.5 & na.omit(seedbank$Predation)==1))


## arsenic in drinking water

arsenic <- read.csv(paste (dataDIR, "arsenic.csv", sep="/"),header=TRUE)
arsenic$Gender<- "Female"
arsenic$Gender[arsenic$gender==1]<- "Male"

arsenic$Type<- "Bladder"
arsenic$Type[arsenic$type==1]<- "Lung"

#### Poisson regression
## Arsenic in drinking water

As.m1 <- glm(events ~ conc, data=arsenic, family="poisson")
display(As.m1, 4)

As.m2 <- glm(events ~ log(conc+1), data=arsenic, family="poisson")
display(As.m2, 4)

As.m3 <- glm(events ~ conc + gender+type, data=arsenic, family="poisson")
display(As.m3, 4)

##trellis.device(postscript, file="Ascancer1.eps", height=4, width=4, horizontal=F)
tikz( file=paste(plotDIRch8, "Ascancer1.tex", sep="/"),
     height=4, width=4, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.5 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.5
xyplot(events~conc|Type*Gender, xlab="As Concentration (ppb)", ylab="Cancer Deaths", data=arsenic, par.settings=trellis.par.temp)
dev.off()

##trellis.device(postscript, file="Ascancer2.eps", height=4, width=4, horizontal=F)
tikz( file=paste(plotDIRch8, "Ascancer2.tex", sep="/"),
     height=4, width=4, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.5 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.5
xyplot(log(events+1)~log(conc+1)|Type*Gender, xlab="Log As Concentration (ppb)",
       ylab="Log Cancer Deaths", data=arsenic, par.settings=trellis.par.temp)
dev.off()

##trellis.device(postscript, file="Ascancer3.eps", height=4, width=4, horizontal=F)
tikz( file=paste(plotDIRch8, "Ascancer3.tex", sep="/"),
     height=4, width=4, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.5 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.5
xyplot(log(at.risk)~log(conc+1)|Type*Gender, xlab="Log As Concentration (ppb)",
       ylab="Log At Risk", data=arsenic, par.settings=trellis.par.temp)
dev.off()

##trellis.device(postscript, file="Ascancer4.eps", height=4, width=4, horizontal=F)
tikz( file=paste(plotDIRch8, "Ascancer4.tex", sep="/"),
     height=4, width=4, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.5 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.5
xyplot((events/at.risk)~log(conc+1)|Type*Gender, xlab="Log As Concentration (ppb)", ylab="Cancer Deaths/At Risk", data=arsenic, par.settings=trellis.par.temp)
dev.off()

As.m4 <- glm(events ~ log(conc+1) + gender + type, data=arsenic, offset=log(at.risk/100000),
             family="poisson")
display(As.m4, 4)

pred.bush <- predict(As.m4, newdata=data.frame(conc=c(rep(10,4), rep(50,4)),
                                gender=rep(c(0,0,1,1), 2), type=rep(c(0,1), 4),
                                at.risk=rep(100000, 8)),
    type="response")
data.frame(effect=pred.bush, conc=c(rep(10,4), rep(50,4)),
           gender=rep(c(0,0,1,1), 2), type=rep(c(0,1), 4), at.risk=rep(1, 8))

 matrix(pred.bush[1:4], 2, 2)
 matrix(pred.bush[5:8], 2, 2)

### overdispersion

As.yhat <- predict(As.m4, type="response")
As.z <- (arsenic$events - As.yhat)/sqrt(As.yhat)
overD <- sum(As.z^2)/summary(As.m4)$df[2]
p.value <- 1-pchisq(sum(As.z^2), summary(As.m4)$df[2])

##postscript(file="Asm4OverD.eps", width=5, height=3, horizontal=F)
tikz(file=paste(plotDIRch8, "Asm4OverD.tex", sep="/"),
     width=5, height=3, standAlone=F)
par(mfrow=c(1,2), mar=c(3,3,3,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(As.yhat, arsenic$events - As.yhat, xlab="Predicted Values", ylab="Residuals", main="Raw Residuals")
plot(As.yhat, As.z, xlab="Predicted Values", ylab="Residuals", main="Standardized Residuals")
abline(h=c(-2,2), lty=2)
dev.off()

As.m5 <- glm(events ~ log(conc+1) + gender + type, data=arsenic,
             offset=log(at.risk), family="quasipoisson")
display(As.m5, 4)

require(MASS)
As.m5nb <- glm.nb(events ~ log(conc+1) + gender + type+offset(log(at.risk)),
             data=arsenic)
summary(As.m5nb)

As.m5.5 <- glm(events ~ conc + gender + type, data=arsenic, offset=log(at.risk), family="quasipoisson")
display(As.m5.5, 4)

As.m6 <- glm(events ~ log(conc+1) * gender  * type, data=arsenic, offset=log(at.risk), family="poisson")
display(As.m6, 4)

As.m6nb<- glm.nb(events~log(conc+1)*gender*type+offset(log(at.risk)),
                 data=arsenic)

As.m7 <- update(As.m6, .~. -log(conc+1):gender:type)
display(As.m7)

As.m7nb <- update(As.m6nb, .~.-log(conc+1):gender:type)
summary(As.m7nb)

As.m8 <- update(As.m7, .~., family="quasipoisson")
display(As.m8, 4)

As.m9 <- update(As.m8, .~.-gender:type)
display(As.m9, 4)

m9.coef <- coef(As.m9)
#        (Intercept)        log(conc + 1)               gender
#          -10.526475             0.469123             0.585474
#                type log(conc + 1):gender   log(conc + 1):type
#            1.478854            -0.101225            -0.187271

##postscript(file="AsModel9.eps", width=4, height=3, horizontal=F)
tikz(file=paste(plotDIRch8, "AsModel9.tex", sep="/"),
           width=4, height=3, standAlone=F)
par(mgp=c(1.25,0.125,0), las=1, tck=0.01, mar=c(3,3,1,0.5))
with(arsenic, {
plot(log(conc+1), 100000*events/at.risk, axes=F,ylim=c(0,90),type="n",
     xlab="As concentration (ppb)", ylab="Cancer Deaths per 100,000")
axis(1, at=log(c(0, 10, 50, 100, 500, 1000)+1), label=c("0","10","50","100","500","1000"))
axis(2)
box()
curve(100000*exp(m9.coef[1]+m9.coef[2]*x), add=T, lty=1)       ##female, baldder
curve(100000*exp(m9.coef[1]+m9.coef[2]*x+m9.coef[3]+m9.coef[5]*x), add=T, lty=2) ## male bladder
curve(100000*exp(m9.coef[1]+m9.coef[2]*x+m9.coef[4]+m9.coef[6]*x), add=T, lty=3) ## female lung
curve(100000*exp(m9.coef[1]+m9.coef[2]*x+m9.coef[3]+m9.coef[4]+m9.coef[5]*x+m9.coef[6]*x), add=T, lty=4) ## male lung
abline(v=log(c(10,50)+1), col="gray", lwd=2)
legend(x=0, y=80, legend=c("female/bladder","male/bladder","female/lung","male/lung"), lty=1:4, cex=0.5)
})
dev.off()

As.pred <- predict(As.m9, type="response")

arsenic$age.c1 <- arsenic$age - mean(arsenic$age)
As.m10<-update(As.m9, .~.+age.c1*gender+age.c1:type)
display(As.m10, 4)

### not used
arsenic$age.c2<-arsenic$age/mean(arsenic$age)
As.m11<-update(As.m9, .~.+log(age.c2)*gender+log(age.c2):type)
display(As.m11, 4)

xyplot(log((events+1)/at.risk)~log(age)|type*gender, data=arsenic)

range(arsenic$age.c1)

pred.data <- data.frame(conc=rep(seq(0, 1000, 10), 4),
                        type=rep(rep(c(0,1), each=101), 2),
                        gender=rep(c(0,1), each=202))
pred.age1 <- predict(As.m10, newdata=data.frame(pred.data, age.c1= rep(-15, 404),
                                                at.risk=rep(100000, 404)), type="response")
pred.age2 <- predict(As.m10, newdata=data.frame(pred.data, age.c1= rep(0  , 404),
                                                at.risk=rep(100000, 404)), type="response")
pred.age3 <- predict(As.m10, newdata=data.frame(pred.data, age.c1= rep( 15, 404),
                                                at.risk=rep(100000, 404)), type="response")
plot.data1 <- data.frame(events=c(pred.age1, pred.age2, pred.age3),
                         rbind(pred.data, pred.data, pred.data),
                         age=rep(c(-15+52.5, 52.5, 52.5+15), each=404))
plot.data1$Type <- "Lung Cancer"
plot.data1$Type[plot.data1$type==0] <- "Bladder Cancer"
plot.data1$Gender <- "Male"
plot.data1$Gender[plot.data1$gender==0] <- "Female"

##trellis.device(postscript, file="AsModel10.eps", width=4.5, height=5, horizontal=F)
tikz(file=paste(plotDIRch8, "AsModel10.tex", sep="/"),
     height=5, width=4.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.set(list(fontsize=list(text=8),
                 par.xlab.text=list(cex=1.25),
                     add.text=list(cex=1.25),
                     superpose.symbol=list(cex=1)))
key <- simpleKey(unique(as.character(plot.data1$age)), lines=T, points=F, space = "top", columns=3)
key$text$cex <- 1.25
xyplot(events~log(conc+1)|Type*Gender, data=plot.data1, type="l",
    group=plot.data1$age, key=key, xlab="As concentration (ppb)", ylab="Cancer deaths per 100,000",
    panel=function(x,y,...){
        panel.xyplot(x,y,lwd=1.5,...)
#        panel.abline(v=log(c(10,50)+1), col="gray")
        panel.grid()
    },
    scales=list(x=list(at=log(c(0, 10, 50, 100, 500, 1000)+1), labels=as.character(c(0, 10, 50, 100, 500, 1000))))
    )
dev.off()
As.yhat <- predict(As.m10, type="response")
As.z <- (arsenic$events - As.yhat)/sqrt(As.yhat)
overD <- sum(As.z^2)/summary(As.m10)$df[2]
p.value <- 1-pchisq(sum(As.z^2), summary(As.m10)$df[2])

##postscript(file="Asm10OverD.eps", width=5, height=3, horizontal=F)
tikz(file=paste(plotDIRch8, "Asm10OverD.tex", sep="/"),
     width=5, height=3, standAlone=F)
par(mfrow=c(1,2), mar=c(3,3,3,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(As.yhat, arsenic$events - As.yhat, xlab="Predicted Values", ylab="Residuals", main="Raw Residuals")
plot(As.yhat, As.z, xlab="Predicted Values", ylab="Residuals", main="Standardized Residuals")
abline(h=c(-2,2), lty=2)
dev.off()

#### Multinomial regression
packages(nnet)

euseSPdata <- read.csv(paste(dataDIR, "usgsmultinomial.csv", sep="/"),
                       header=TRUE)
names(euseSPdata)[33] <- "Unknown" ## was named "Other"
## fit a multinomial logit model where PID are column numbers of the
## response variable, and the predictor is % developed land (NLCD2):

PID <- c(33, 30:32)
multinom.BOS1 <- multinom(as.matrix(euseSPdata[,PID])~NLCD2,
                          data=euseSPdata, subset=City=="BOS")

 ## fitting with logit transformed x:
multinom.BOS2 <- multinom(as.matrix(euseSPdata[,PID]) ~ logit(NLCD2),
                          data=euseSPdata, subset=City=="BOS")

 ## print summary
 summary(multinom.BOS1, corr=FALSE)
 summary(multinom.BOS2, corr=FALSE)

anova(multinom.BOS1, multinom.BOS2)

## prediction

pp <- predict(multinom.BOS1, type="probs",
              newdata=data.frame(NLCD2=0:100), se.fit=T)
p0 <- pp[,2]
p1 <- pp[,3]
p2 <- pp[,4]
p3 <- pp[,1]

## plots without error bands
tikz(file=paste(plotDIRch8, "tolerance1.tex", sep="/"),
           height=4, width=4, standAlone=F)
par(mfrow=c(2,2), mar=c(0.,0.,0.,0.),oma=c(3,3,3,3), mgp=c(1.25,0.25,0),
    tck=-0.015, las=1)
plot(0:100, p0, xlab="", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,1], col="gray", cex=0.5)
lines(0:100, p0)
text(50, 0.9, "Intolerant", cex=0.75)
axis(3, outer=T)
axis(2,labels=FALSE)
box()
plot(0:100, p1, xlab="", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,2], col="gray", cex=0.5)
lines(0:100, p1)
text(50, 0.9, "Moderate Tol", cex=0.75)
box()
axis(4, outer=T)
axis(3, labels=FALSE)
plot(0:100, p2, xlab="NUII", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,3], col="gray", cex=0.5)
lines(0:100, p2)
text(50, 0.9, "Tolerant", cex=0.75)
box()
axis(2, outer=T)
axis(1, labels=FALSE)
plot(0:100, p3, xlab="NUII", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,4], col="gray", cex=0.5)
lines(0:100, p3)
text(50, 0.9, "Unknown", cex=0.75)
box()
axis(1, outer=T)
axis(4, labels=FALSE)
mtext(side=1, "\\% developed land", outer=T, line=1.5)
mtext(side=2, "relative abundance", outer=T, line=1.75, las=0)
dev.off()


### simulation

## function for generating random samples of model coefficients
## assuming estimated model coefficients have a multivariate
## normal distribution

sim.multinom <- function(M, n.sims=NULL){
    ## M: a multinomial model object of class "multinom"
    ## n.sims: number of Monte Carlo somulations
    packages(rv)
    ## a package for random variate simulation and calculation
    if (is.null(n.sims)) n.sims <- getnsims()
    else setnsims(n.sims)
    ## setting simulation numbers to be either user supplied
    ## or rv package default (2500)
    object.class <- class(M)
    if(object.class[1]!="multinom") stop ("Not a multinom object")
    
    summ <- summary(M)
    beta.hat <- as.vector(t(coef(M)))
    V.beta <- vcov(M)
    k <- length(beta.hat)
    beta <- array(NA, c(n.sims, k))
    lbs <- labels(coef(M))
    dmnm <- character()
    for (i in 1:length(lbs[[1]]))
        dmnm <- c(dmnm, paste(lbs[[1]][i], lbs[[2]], sep=":"))
    dimnames(beta) <- list(NULL, dmnm)
    beta <- mvrnorm(n.sims, beta.hat, V.beta)
    return(beta)
}

sim.BOS <- rvsims(sim.multinom(multinom.BOS1, 2500))
## generating random samples of model coefficients and store them
## as an rv object

sim.BOS <- rvmatrix(sim.BOS, nrow=3,ncol=2, byrow=T)
X <- cbind(1,seq(0,100,1))
Xb1 <- X[,1]*sim.BOS[1,1]+X[,2]*sim.BOS[1,2]
Xb2 <- X[,1]*sim.BOS[2,1]+X[,2]*sim.BOS[2,2]
Xb3 <- X[,1]*sim.BOS[3,1]+X[,2]*sim.BOS[3,2]

denomsum <- 1+exp(Xb1) + exp(Xb2) + exp(Xb3)
p3 <- 1/denomsum
p0 <- exp(Xb1)/denomsum
p1 <- exp(Xb2)/denomsum
p2 <- exp(Xb3)/denomsum

rows <- euseSPdata$City=="BOS" 
## producing a postscript file of figure 1
tikz(file=paste(plotDIRch8, "tolerance1rv.tex", sep="/"),
     height=4, width=4, standAlone=F)
par(mfrow=c(2,2), mar=c(0.,0.,0.,0.),oma=c(3,3,3,3),
    mgp=c(1.25,0.25,0), tck=0.015, las=1)
plot(p0, xlab="", ylab="", ylim=c(0,1),
     col=1, cex=0.25, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,1], col="gray", cex=0.5)
text(50, 0.9, "Intolerant", cex=0.75)
axis(3, outer=T)
axis(2,labels=FALSE)
box()
plot(p1, xlab="", ylab="", ylim=c(0,1), col=1,
     cex=0.25, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,2], col="gray", cex=0.5)
text(50, 0.9, "Moderate Tol", cex=0.75)
box()
axis(4, outer=T)
axis(3, labels=FALSE)
plot(p2, xlab="", ylab="", ylim=c(0,1), col=1,
     cex=0.25, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,3], col="gray", cex=0.5)
text(50, 0.9, "Tolerant", cex=0.75)
box()
axis(2, outer=T)
axis(1, labels=FALSE)
plot(p3, xlab="", ylab="", ylim=c(0,1), col=1,
     cex=0.25, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,4], col="gray", cex=0.5)
text(50, 0.9, "Unknown", cex=0.75)
box()
axis(1, outer=T)
axis(4, labels=FALSE)
mtext(side=1, "\\% developed land", outer=T, line=1.5)
mtext(side=2, "probability of occurrence", outer=T, line=1.75, las=0)
dev.off()

## residuals
cities <- levels(euseSPdata$City)
dataP <- t(apply(euseSPdata[,30:33], 1, function(x) return(x/sum(x))))
####### assessing model's fit
txs <- names(euseSPdata)[30:33]
Rows <- euseSPdata$City=="BOS"
M <- multinom(as.matrix(euseSPdata[,30:33])~NLCD2,data=euseSPdata, subset=City=="BOS")
resid.mat <- dataP[Rows,] - fitted(M)
Resid <- unlist(resid.mat)
Fitted <- unlist(fitted(M))
Group <- rep(txs, each=dim(resid.mat)[1])
X <- rep(euseSPdata$NLCD2[Rows], dim(resid.mat)[2])

plot(euseSPdata$NLCD2[Rows], resid.mat[,1], pch=1, cex=0.5, ylim=range(resid.mat))
abline(0,0)
points(euseSPdata$NLCD2[Rows], resid.mat[,2], pch=2, cex=0.5, col=2)
points(euseSPdata$NLCD2[Rows], resid.mat[,3], pch=3, cex=0.5, col=3)
points(euseSPdata$NLCD2[Rows], resid.mat[,4], pch=4, cex=0.5, col=4)

plot(fitted(M)[,1], resid.mat[,1], pch=1, cex=0.5,
     ylim=range(resid.mat), xlim=range(fitted(M)))
abline(0,0)
points(fitted(M)[,2], resid.mat[,2], pch=2, cex=0.5, col=2)
points(fitted(M)[,3], resid.mat[,3], pch=3, cex=0.5, col=3)
points(fitted(M)[,4], resid.mat[,4], pch=4, cex=0.5, col=4)

stdResid <- as.vector(Resid/sqrt(Fitted*(1-Fitted)))
FittedV <- as.vector(Fitted)
GroupV <- rep(c("Intol","Modtol","Tol","Unknown"), each=dim(Fitted)[1])
tikz(file=paste(plotDIRch8, "residVfitted.tex", sep="/"),
               height=3.5, width=4.25, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=F))
trellis.par.set(list(fontsize=list(text=8),
                 par.xlab.text=list(cex=1.25),
                     add.text=list(cex=1.25),
                     superpose.symbol=list(cex=1)))
key <- simpleKey(c("Intol","Modtol","Tol","Unknown"), lines=F, points=T,
                 space = "top", columns=4)
key$text$cex <- 1.2
xyplot(stdResid~FittedV, group=GroupV,
       ylab="standardized residuals",xlab="Fitted Relative Abundances", key=key)
dev.off()


### Poisson -- multinomial connection
BOS.data <- euseSPdata[euseSPdata$City=="BOS", ]
pp <- predict(multinom.BOS1, type="probs",
              newdata=data.frame(NLCD2=0:100), se.fit=T)
p0 <- pp[,2]
p1 <- pp[,3]
p2 <- pp[,4]
p3 <- pp[,1]

## Independent Poisson models
BOS.p1 <- BOS.p2 <- list()
for (i in 1:4){
  BOS.p1[[i]] <- glm(BOS.data[,30-1+i] ~ NLCD2, data=BOS.data, family="poisson")
  BOS.p2[[i]] <- glm(BOS.data[,30-1+i] ~ NLCD2+I(log(NLCD2)), data=BOS.data, family="poisson")
}

## reproducing Qian et al (2012) figure 1
x <- 1:100
predP2.mat <- predP1.mat <- matrix(0, nrow=4, ncol=length(x))
for (i in 1:4){
  predP1.mat[i,] <- predict(BOS.p1[[i]], new=data.frame(NLCD2=x))
  predP2.mat[i,] <- predict(BOS.p2[[i]], new=data.frame(NLCD2=x))
}
pi.p1 <- pi.p2 <- matrix(0, nrow=4, ncol=length(x))
for (i in 1:length(x)){
  pi.p1[,i] <- exp(predP1.mat[,i])/(sum(exp(predP1.mat[,i])))
  pi.p2[,i] <- exp(predP2.mat[,i])/(sum(exp(predP2.mat[,i])))
}

## figure 1
rows <- euseSPdata$City=="BOS"
tikz(file=paste(plotDIRch8, "tolgroups1.tex", sep="/"),
     height=4, width=4, standAlone=F)
#postscript(file=paste(plotDIR, "PoisMultnFig1.eps", sep="/"),
#           height=4, width=4, horizontal=F)
par(mfrow=c(2,2), mar=c(0.,0.,0.,0.),oma=c(3,3,3,3), mgp=c(1.25,0.25,0),
    tck=-0.015, las=1)
plot(0:100, p0, xlab="", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,1], col="gray", cex=0.5)
lines(0:100, p0)
lines(1:100, pi.p1[1,], lty=2)
text(50, 0.9, "Intolerant", cex=0.75)
axis(3, outer=T)
axis(2,labels=FALSE)
box()
plot(0:100, p1, xlab="", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,2], col="gray", cex=0.5)
lines(0:100, p1)
lines(1:100, pi.p1[2,], lty=2)
text(50, 0.9, "Moderate Tol", cex=0.75)
box()
axis(4, outer=T)
axis(3, labels=FALSE)
plot(0:100, p2, xlab="NUII", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,3], col="gray", cex=0.5)
lines(0:100, p2)
lines(1:100, pi.p1[3,], lty=2)
text(50, 0.9, "Tolerant", cex=0.75)
box()
axis(2, outer=T)
axis(1, labels=FALSE)
plot(0:100, p3, xlab="NUII", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,4], col="gray", cex=0.5)
lines(0:100, p3)
lines(1:100, pi.p1[4,], lty=2)
text(50, 0.9, "Unknown", cex=0.75)
box()
axis(1, outer=T)
axis(4, labels=FALSE)
mtext(side=1, "\\% developed land", outer=T, line=1.5)
mtext(side=2, "relative abundance", outer=T, line=1.75, las=0)
dev.off()

lgnd <- c("Intolerant", "Moderate Tol", "Tolerant", "Unknown")

## total abundance plot (figure 2)
tikz(file=paste(plotDIR, "tolpoisson.tex", sep="/"),
     height=4, width=4, standAlone=F)
#postscript(file=paste(plotDIR, "PoisMultnFig2.eps", sep="/"),
#           height=4, width=4, horizontal=F)
par(mfrow=c(2,2), mar=c(0.,0.,0.,0.),oma=c(3,3,3,3), mgp=c(1.25,0.25,0),
    tck=-0.015, las=1)
for (i in 1:4){
  plot(BOS.data$NLCD2, BOS.data[,30-1+i], xlab=" ", ylab="", ylim=range(BOS.data[,30:33]),
       xlim=c(0,100), axes=F)
  if (i ==1) {
    axis(2, labels=FALSE)
    axis(3)
  }
  if (i == 2){
    axis(3, labels=FALSE)
    axis(4)
  }
  if (i == 3){
    axis(1, labels=FALSE)
    axis(2)
  }
  if (i ==4) {
    axis(4, labels=FALSE)
    axis(1)
  }
  lines(x, exp(predP1.mat[i,]))
  ##  lines(x, exp(predP2.mat[i,]), col="blue")
  text(50, 27, lgnd[i], cex=0.75)
  box()
}
mtext(side=1, "\\% developed land", outer=T, line=1.5)
mtext(side=2, "total abundance", outer=T, line=1.75, las=0)
dev.off()


### mayfly taxa
rthE <- read.csv(paste(dataDIR, "rthE.csv",sep="/"), header=T)
NN <- (1:dim(rthE)[2])[names(rthE)=="MANUII"] -1

non0row <- apply(rthE[,1:NN], 1, sum) != 0
BOS.data <- rthE[rthE$City=="BOS" & non0row,]
non0col<- apply(BOS.data[,1:NN], 2, FUN=function(x) return(sum(x>0)>2))
## removing columns with 2 or fewer non-zero counts

m0 <- multinom(as.matrix(BOS.data[,1:NN][,non0col])~NLCD2,data=BOS.data,maxit=1000)
      ## the model
m1 <- multinom(as.matrix(BOS.data[,1:NN][,non0col])~logit(NLCD2),
               data=BOS.data,maxit=1000)
anova(m0,m1)

## independent Poisson models
non0BOS <- cbind(BOS.data[,1:NN][,non0col], BOS.data$NLCD2)
names(non0BOS) <- c(paste("e", 1:sum(non0col), sep=""), "NLCD2")

BOSe1 <- BOSe2 <- list()
BOSe1[[1]] <- glm(e1~logit(NLCD2), family="poisson", data=non0BOS)
BOSe1[[2]] <- glm(e2~logit(NLCD2), family="poisson", data=non0BOS)
BOSe1[[3]] <- glm(e3~logit(NLCD2), family="poisson", data=non0BOS)
BOSe1[[4]] <- glm(e4~logit(NLCD2), family="poisson", data=non0BOS)
BOSe1[[5]] <- glm(e5~logit(NLCD2), family="poisson", data=non0BOS)
BOSe1[[6]] <- glm(e6~logit(NLCD2), family="poisson", data=non0BOS)
BOSe1[[7]] <- glm(e7~logit(NLCD2), family="poisson", data=non0BOS)
BOSe1[[8]] <- glm(e8~logit(NLCD2), family="poisson", data=non0BOS)
BOSe1[[9]] <- glm(e9~logit(NLCD2), family="poisson", data=non0BOS)
BOSe1[[10]] <- glm(e10~logit(NLCD2), family="poisson", data=non0BOS)

BOSe2[[1]] <- glm(e1~NLCD2+I(log(NLCD2)), family="poisson", data=non0BOS)
BOSe2[[2]] <- glm(e2~NLCD2+I(log(NLCD2)), family="poisson", data=non0BOS)
BOSe2[[3]] <- glm(e3~NLCD2+I(log(NLCD2)), family="poisson", data=non0BOS)
BOSe2[[4]] <- glm(e4~NLCD2+I(log(NLCD2)), family="poisson", data=non0BOS)
BOSe2[[5]] <- glm(e5~NLCD2+I(log(NLCD2)), family="poisson", data=non0BOS)
BOSe2[[6]] <- glm(e6~NLCD2+I(log(NLCD2)), family="poisson", data=non0BOS)
BOSe2[[7]] <- glm(e7~NLCD2+I(log(NLCD2)), family="poisson", data=non0BOS)
BOSe2[[8]] <- glm(e8~NLCD2+I(log(NLCD2)), family="poisson", data=non0BOS)
BOSe2[[9]] <- glm(e9~NLCD2+I(log(NLCD2)), family="poisson", data=non0BOS)
BOSe2[[10]] <- glm(e10~NLCD2+I(log(NLCD2)), family="poisson", data=non0BOS)

##
x <- seq(0, 60, 1)
pred2.mat <- pred1.mat <- matrix(0, nrow=10, ncol=length(x))
for (i in 1:10){
  pred1.mat[i,] <- predict(BOSe1[[i]], new=data.frame(NLCD2=x))
  pred2.mat[i,] <- predict(BOSe2[[i]], new=data.frame(NLCD2=x))
}
pi.e1 <- pi.e2 <- matrix(0, nrow=10, ncol=length(x))
for (i in 1:length(x)){
  pi.e1[,i] <- exp(pred1.mat[,i])/(sum(exp(pred1.mat[,i])))
  pi.e2[,i] <- exp(pred2.mat[,i])/(sum(exp(pred2.mat[,i])))
}
dataP <- t(apply(non0BOS[,-11], 1, function(x) return(x/sum(x))))

tikz(paste(plotDIRch8, "PoisMultnP.tex", sep="/"),
     height=2.75, width=5.75, standAlone=F)
par(mfrow=c(2, 5), mar=c(0.05, 0.05, 0.05, 0.05), oma=c(3,3.25,1,1), mgp=c(1.5,0.25,0),
    las=1, tck= -0.015)
for (i in 1:10){
  plot(non0BOS$NLCD2, non0BOS[,i], xlab=" ", ylab="", ylim=range(non0BOS[,-11]),
       xlim=c(0,60), axes=F)
  if (i ==1 | i==6) axis(2)
  if (i !=1 & i!=6) axis(2, label=F)
  if (i < 6 ) axis(1, label=F)
  if (i >= 6 ) axis(1)
  lines(x, exp(pred1.mat[i,]), lty=1, lwd=2)
  lines(x, exp(pred2.mat[i,]), lty=2, lwd=2)
  text(5, 900, paste("E", i, sep=""), cex=0.75)
#  box()
}
mtext("\\% developed land", side=1, outer=T, line=1.5)
mtext("Abundance", side=2, outer=T, line=2, las=0)
dev.off()

EPTeDataplot <- function(data, NN, city, x=cbind(1, 0:100)){
  non0r<-apply(data[,1:NN], 1, sum)!=0  ## all 0 rows
  dataB <- data[data$City==city & non0r,]
  non0c<-apply(dataB[,1:NN], 2, FUN=function(x)return(sum(x>0)>2)) ## taxa appear in 2 or less sites
  sp <- names(dataB)[1:NN][non0c]        ## remaining taxa names
  m <- multinom(as.matrix(dataB[,1:NN][,non0c])~NLCD2,data=dataB,
                subset=City==city, maxit=1000)  ## the model
  beta <- as.matrix(coef(m))
  n.taxa <- sum(non0c)
  Xb <- matrix(0, nrow=dim(x)[1], ncol=n.taxa-1)
  for (i in 1:(n.taxa-1)) Xb[,i] <- x %*% beta[i,]
  denomsum <- apply(exp(Xb), 1, sum)
  pp <-matrix(0, nrow=dim(x)[1], ncol=n.taxa)
  pp[,1] <- 1/(1+denomsum)
  for (i in 2:n.taxa) pp[,i] <- exp(Xb[,i-1])/(1+denomsum)
  dataP <- t(apply(dataB[,1:NN][,non0c], 1, function(x) return(x/sum(x))))
  np <- sum(non0c)
  par(mfrow=rep(floor(sqrt(np))+1, 2), mar=c(0,0,0,0),oma=c(3,3,3,3), mgp=c(1.25,0.25,0),
      tck=-0.015, las=1)
  for (i in 1:np){
    plot(x[,2], pp[,1], xlab="", ylab="", ylim=c(0,1),
       type="n", col=1, pch=16, axes=F)
  points(dataB$NLCD2, dataP[,i], col="red", cex=0.5)
    lines(x[,2], pp[,i])
  text(sum(range(x[,2]))/2, 0.9, names(dataB)[non0c][i], cex=0.75)
  axis(1, outer=T)
    axis(2, outer=T)
    box()
  }
 mtext(city, outer=T, cex=0.75)
 invisible(pp)
}
multi.pred <- EPTeDataplot(data=rthE, NN=NN, city="BOS", x=cbind(1, x))

tikz(file=paste(plotDIRch8, "PoisMultnPM.tex", sep="/"),
     height=2.75, width=5.75, standAlone=F)
par(mfrow=c(2, 5), mar=c(0.05, 0.05, 0.05, 0.05), oma=c(3,3.25,1,1), mgp=c(1.5,0.25,0),
    las=1, tck= -0.015)
for (i in 1:10){
  plot(non0BOS$NLCD2, dataP[,i], xlab=" ", ylab="",
       ylim=c(0,1), xlim=c(0,60), axes=F)
  if (i ==1 | i==6) axis(2)
  if (i !=1 & i!=6) axis(2, label=F)
  if (i < 6 ) axis(1, label=F)
  if (i >= 6 ) axis(1)
  lines(x, pi.e1[i,], lty=1, lwd=2)
  lines(x, pi.e2[i,], lty=2, lwd=2)
  lines(x, multi.pred[,i], lty=3, lwd=2)
  text(5, 0.900, paste("E", i, sep=""), cex=0.75)
}
mtext("\\% developed land", side=1, outer=T, line=1.5)
mtext("Relative Abundance", side=2, outer=T, line=1.75, las=0)
dev.off()


#### GAM ####
### cart&gam.r
whale.data <- read.csv(paste(dataDIR, "antarctic.csv", sep="/"),
                       header=T)
packages(maps)


my.box<-function(xlim, ylim, ...){
    segments(x0=xlim, y0=rep(ylim[1],2), x1=xlim, y1=rep(ylim[2], 2), ...)
    segments(y0=ylim, x0=rep(xlim[1],2), y1=ylim, x1=rep(xlim[2], 2), ...)
}

## Figure 8.21
##postscript(file=paste(plotDIR, "maps.eps", sep="/"), horizontal=F, height=5, width=4.5)
tikz(file=paste(plotDIRch8, "maps.tex", sep="/"),
     standAlone=F, height=5, width=4.5)
par(mar=c(0,0,0.5,0))
map("world", xlim = c(-90,-60), ylim = c(-74,-63), col="grey80", fill=T)
#map("nz", fill=TRUE, col="grey80")
# mini world map as guide
points(whale.data$Lon[(whale.data$MN+whale.data$BA)==0], whale.data$Lat[(whale.data$MN+whale.data$BA)==0], pch=3, col="gray", cex=0.5)
points(whale.data$Lon[(whale.data$MN+whale.data$BA)!=0], whale.data$Lat[(whale.data$MN+whale.data$BA)!=0], pch=1, col=1, cex=0.5)
maplocs <- map(projection="sp_mercator", wrap=TRUE, lwd=0.1,
               col="grey", xlim=c(-180, 0),
               interior=FALSE, orientation=c(90, 180, 0), add=TRUE,
               plot=FALSE)
xrange <- range(maplocs$x, na.rm=TRUE)
yrange <- range(maplocs$y, na.rm=TRUE)
aspect <- abs(diff(yrange))/abs(diff(xrange))
# customised to 6.5 by 4.5 figure size
par(fig=c(0.01, 0.99 - 0.5, 0.99 - 0.5*aspect*4.5/6.5, 0.99),
    mar=rep(0, 4), new=TRUE)
plot.new()
plot.window(xlim=xrange,
            ylim=yrange)
map(projection="sp_mercator", wrap=TRUE, lwd=0.25,
    interior=FALSE, orientation=c(90, 180, 0), add=TRUE)
my.box(xlim= c(1.5,2.25), ylim = c(-2.,-1.35))
#symbols(-75, -67.5, squares=c(30, 15), inches=0.1, add=TRUE)

dev.off()

## a second version of the map
##postscript(file=paste(plotDIR, "maps2.eps", sep="/"), horizontal=F, height=5, width=4)
tikz(file=paste(plotDIRch8, "maps2.tex", sep="/"),
     standAlone=F, height=5, width=4)
par(mfrow=c(2,1), mar=c(1,1,1,1))
map("world", xlim=c(-180,0), col="gray")
my.box(xlim= c(-90,-60), ylim = c(-75,-60))

map("world", xlim = c(-90,-60), ylim = c(-75,-60))
points(whale.data$Lon[(whale.data$MN+whale.data$BA)!=0], whale.data$Lat[(whale.data$MN+whale.data$BA)!=0], pch=16, cex=0.5)
points(whale.data$Lon[(whale.data$MN+whale.data$BA)==0], whale.data$Lat[(whale.data$MN+whale.data$BA)==0], pch=1,  col="gray", cex=0.5)
dev.off()

## Figure 8.22
##postscript(file=paste(plotDIR, "whaledata.eps", sep="/"), height=5, width=4.5, horizontal=F)
tikz(file=paste(plotDIRch8, "whaledata.tex", sep="/"),
     height=5, width=4.5, standAlone=F)
par(mfrow=c(3,2), mar=c(3,3,0.25,0.25), mgp=c(1.25,0.125,0), tck=0.01)
attach(whale.data)
plot(bathy, TW, xlab="Bathymetry (m)", ylab="No. whales", pch=1, cex=0.5)
plot(chla, TW, xlab="Chla (mg/L)", ylab="No. whales", pch=1, cex=0.5)
plot(S.bathy, TW, xlab="Bathymetry slope (\\%)", ylab="No. whales", pch=1, cex=0.5)
plot(A.v100, TW, xlab="Acoustic (dB)", ylab="No. whales", pch=1, cex=0.5)
#plot(W.mass, TW, xlab=expression(paste(" Water Temperature ( ", degree, "C)")), ylab="No. Whales", pch=1, cex=0.5)
plot(D.ice, TW, xlab="Dist. to ice edge (km)", ylab="No. whales", pch=1, cex=0.5)
plot(D.coast, TW, xlab="Dist. to shore (km)", ylab="No. whales", pch=1, cex=0.5)

detach()
dev.off()

whale.data$D.coast <- whale.data$D.coast/1000
whale.data$D.ice <- whale.data$D.ice/1000

packages(rpart)
set.seed(1230)
TW.rpart<-rpart(TW~bathy+chla+D.coast+D.ice+D.inswb+D.slp+S.bathy+W.mass+A.v100+A.v300.2,
            data=whale.data, method="poisson", cp=0.00)
printcp(TW.rpart)

## Figure 8.23
##postscript(file=paste(plotDIR, "whaleplotcp.eps", sep="/"), height=2.75, width=4, horizontal=F)
tikz(file=paste(plotDIRch8, "whaleplotcp.tex", sep="/"),
     height=2.75, width=4, standAlone=F)
par(mar=c(3,3,0.25,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plotcp(TW.rpart)
dev.off()

## Figure 8.24
TW.rpart.prune <- prune.rpart(TW.rpart, cp=0.037)
##postscript(file=paste(plotDIR, "whalecart.eps", sep="/"), height=2.75, width=4, horizontal=F)
tikz(file=paste(plotDIRch8, "whalecart.tex", sep="/"),
     height=2.75, width=4, standAlone=F)
par(mar=c(3,3,0.25,0.25), mgp=c(1.25,0.125,0), tck=0.01, las=1)
plot(TW.rpart.prune, margin=0.2, branch=0.75)
text(TW.rpart.prune, pretty=T, use.n=T, cex=0.75)
dev.off()


set.seed(1230)

TW.rpart2 <- rpart(I(TW>0)~bathy+chla+D.coast+D.ice+D.inswb+D.slp+S.bathy+W.mass+A.v100+A.v300.2,
            data=whale.data, method="class", parms=list(prior=c(0.5,0.5)), cp=0.00)
printcp(TW.rpart2)

##postscript(file=paste(plotDIR, "whaleplotcp2.eps", sep="/"), height=3, width=4, horizontal=F)
tikz(file=paste(plotDIRch8, "whaleplotcp2.tex", sep="/"),
     height=3, width=4, standAlone=F)
par(mar=c(3,3,0.25,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plotcp(TW.rpart2)
dev.off()

TW.rpart2.prune <- prune.rpart(TW.rpart2, cp=0.025)
##postscript(file=paste(plotDIR, "whalecart2.eps", sep="/"), height=3, width=4, horizontal=F)
tikz(file=paste(plotDIRch8, "whalecart2.tex", sep="/"),
     height=3, width=4, standAlone=F)
par(mar=c(3,3,0.25,0.25), mgp=c(1.3,0.4,0), tck=-0.02)
plot(TW.rpart2.prune, margin=0.2, branch=0.75)
text(TW.rpart2.prune, pretty=F, use.n=T, cex=0.75, digits=0)
dev.off()

## Figure 8.25
##postscript(file=paste(plotDIR, "whaleplot2.eps", sep="/"), height=2.5, width=5, horizontal=F)
tikz(file=paste(plotDIRch8, "whaleplot2.tex", sep="/"),
     height=2.5, width=5, standAlone=F)
layout(matrix(c(1,2), nrow=1), c(2, 3))
par(mar=c(3,3,3,1), mgp=c(1.25,0.125,0), tck=0.01, las=1)
plotcp(TW.rpart2)
par(mar=c(1,1,1,1))
plot(TW.rpart2.prune, margin=0.2, branch=0.75)
text(TW.rpart2.prune, pretty=T, use.n=T, cex=0.5, digits=0)
dev.off()

#### gam ####
## using gam for possible variable transformations

## variables picked from CART: A.v50, A.v100, chla, bathy, S.bathy, D.ice, D.inswb
packages(mgcv)
std.chla <- scale(whale.data$chla) #-mean(antarctic.data$chla, na.rm=T))/sd(antarctic.data$chla, na.rm=T)
std.av100 <- scale(whale.data$A.v100) #-mean(antarctic.data$A.v100, na.rm=T))/sd(antarctic.data$A.v100, na.rm=T)
std.av300 <- scale(whale.data$A.v300.2) #-mean(antarctic.data$A.v300.2, na.rm=T))/sd(antarctic.data$A.v300.2, na.rm=T)
std.bathy <- scale(whale.data$bathy) ##-mean(antarctic.data$bathy, na.rm=T))/sd(antarctic.data$bathy, na.rm=T)
std.sbathy <- scale(log(whale.data$S.bathy)) #-mean(log(antarctic.data$S.bathy), na.rm=T))/sd(log(antarctic.data$S.bathy), na.rm=T)
std.Dice <- scale(whale.data$D.ice) #-mean(antarctic.data$D.ice, na.rm=T))/sd(antarctic.data$D.ice, na.rm=T)
std.Dinswb <- scale(whale.data$D.inswb) #-mean(antarctic.data$D.inswb, na.rm=T))/sd(antarctic.data$D.inswb, na.rm=T)
std.Dcoast <- scale(whale.data$D.coast) #-mean(antarctic.data$D.inswb, na.rm=T))/sd(antarctic.data$D.inswb, na.rm=T)

subdata <- data.frame(TW=whale.data$TW, Lat=whale.data$Lat, Lon=whale.data$Lon, std.av100, std.chla, std.bathy,  std.Dice, std.Dinswb, std.Dcoast)
subdata<- na.omit(subdata)

whale.gam1 <- gam(TW~s(std.av100,bs="ts")+s(std.chla,bs="ts")+s(std.bathy,bs="ts")+ s(std.Dice,bs="ts")+s(std.Dcoast,bs="ts"),
    data=subdata, family="poisson")

## Figure 8.26
##postscript(file=paste(plotDIR, "whalegam1.eps", sep="/"), height=5.5, width=5, horizontal=F)
tikz(file=paste(plotDIRch8, "whalegam1.tex", sep="/"),
     height=5.5, width=5, standAlone=F)
par(mfrow=c(3,2), mar=c(4,4,0.5,0.5), mgp=c(1.25,.125,0), las=1, tck=0.01)
plot(whale.gam1, scale=0, pages=0, select=1, xlab="Backscatter 25-100m", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam1, scale=0, pages=0, select=2, xlab="Chlorophyll a", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam1, scale=0, pages=0, select=3, xlab="Bathymetry", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam1, scale=0, pages=0, select=4, xlab="Dist. to ice edge", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam1, scale=0, pages=0, select=5, xlab="Dist. to shore", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
dev.off()

## Standardize predictors did not make much difference
whale.gam11 <- gam(TW~s(A.v100,bs="ts")+s(chla,bs="ts")+s(bathy,bs="ts")+ s(D.ice,bs="ts")+s(D.coast,bs="ts"),
    data=whale.data, family="poisson")
#postscript(file=paste(plotDIR, "whalegam1.eps", sep="/"), height=5.5, width=5, horizontal=F)
par(mfrow=c(3,2), mar=c(4,4,0.5,0.5), mgp=c(1.5,.5,0), tck=-0.02)
plot(whale.gam11, scale=0, pages=0, select=1, xlab="Backscatter 25-100m (dB)", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam11, scale=0, pages=0, select=2, xlab="Chlorophyll a", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam11, scale=0, pages=0, select=3, xlab="Bathymetry (m)", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam11, scale=0, pages=0, select=4, xlab="Dist. to ice edge (km)", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam11, scale=0, pages=0, select=5, xlab="Dist. to shore (km)", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
#dev.off()


As.yhat <- predict(whale.gam1, type="response")
As.z <- (subdata$TW - As.yhat)/sqrt(As.yhat)
overD <- sum(As.z^2)/summary(whale.gam1)$residual.df
p.value <- 1-pchisq(sum(As.z^2), summary(whale.gam1)$residual.df)


## Figure 8.27
##postscript(file=paste(plotDIR, "whaleGAM1resid.eps", sep="/"), width=5, height=3, horizontal=F)
tikz(file=paste(plotDIRch8, "whaleGAM1resid.tex", sep="/"),
     width=5, height=3, standAlone=F)
par(mfrow=c(1,2), mar=c(3,3,3,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(As.yhat, subdata$TW - As.yhat, xlab="Predicted Values", ylab="Residuals")
title(main="Raw Residuals", cex.main=0.75)
abline(h=c(-2,2), lty=2)
plot(As.yhat, As.z, xlab="Predicted Values", ylab="Residuals")
title(main="Standardized Residuals", cex.main=0.75)
abline(h=c(-2,2), lty=2)
dev.off()

whale.gam3 <- gam(I(TW>0)~s(A.v100,bs="ts")+s(chla,bs="ts")+s(bathy,bs="ts")+ s(D.ice,bs="ts")+s(D.coast,bs="ts"),
    data=whale.data, family="binomial")

## Figure 8.28
##postscript(file=paste(plotDIR, "whalegam3.eps", sep="/"), height=5.5, width=5, horizontal=F)
tikz(file=paste(plotDIRch8, "whalegam3.tex", sep="/"),
     height=5.5, width=5, standAlone=F)
par(mfrow=c(3,2), mar=c(4,4,0.5,0.5), mgp=c(1.25,.125,0), las=1, tck=0.01)
plot(whale.gam3, scale=0, pages=0, select=1, xlab="Backscatter 25-100m (dB)", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam3, scale=0, pages=0, select=2, xlab="Chlorophyll a", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam3, scale=0, pages=0, select=3, xlab="Bathymetry (m)", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam3, scale=0, pages=0, select=4, xlab="Dist. to ice edge (km)", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam3, scale=0, pages=0, select=5, xlab="Dist. to shore (km)", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
dev.off()



### Quasi-Poisson would not change the shape of the fitted GAM

whale.gam2 <- gam(TW~s(std.av100,bs="ts")+s(std.chla,bs="ts")+s(std.bathy,bs="ts")+ s(std.Dice,bs="ts")+s(std.Dcoast,bs="ts"),
    data=subdata, family="quasipoisson")

##postscript(file=paste(plotDIR, "whalegam2.eps", sep="/"), height=5.5, width=5, horizontal=F)
tikz(file=paste(plotDIRch8, "whalegam2.tex", sep="/"),
     height=5.5, width=5, standAlone=F)
par(mfrow=c(3,2), mar=c(4,4,0.5,0.5), mgp=c(1.25,.125,0), las=1, tck=0.01)
plot(whale.gam2, scale=0, pages=0, select=1, xlab="Backscatter 25-100m", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam2, scale=0, pages=0, select=2, xlab="Chlorophyll a", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam2, scale=0, pages=0, select=3, xlab="Bathymetry", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam2, scale=0, pages=0, select=4, xlab="Dist. to ice edge", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam2, scale=0, pages=0, select=5, xlab="Dist. to shore", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
dev.off()


whale.gam1 <- gam(TW~s(std.av100,bs="ts")+s(std.chla,bs="ts")+s(std.sbathy,bs="ts")+s(std.av300,bs="ts")+ s(std.Dice,bs="ts")+s(std.dinswb,bs="ts"),
    data=subdata, family="poisson")

plot(subdata$TW, fitted(whale.gam1), ylab="Fitted",xlab="Observed")
abline(0,1)
    ## residuals
par(mfrow=c(2,1), mar=c(4,4,0.5, 1), mgp=c(1.5,0.5,0))
plot(subdata$Lat, (subdata$TW-fitted(whale.gam1))/sqrt(fitted(whale.gam1)), ylab="Residuals", xlab="Latitude")
plot(subdata$Lon, (subdata$TW-fitted(whale.gam1))/sqrt(fitted(whale.gam1)), ylab="Residuals", xlab="Longitude")

par(mfrow=c(2,1), mgp=c(1.5,0.5,0), mar=c(4,4,0.5, 1))
plot(fitted(whale.gam1), (subdata$TW-fitted(whale.gam1)), ylab="Residuals", xlab="Fitted")
plot(fitted(whale.gam1), (subdata$TW-fitted(whale.gam1))/sqrt(fitted(whale.gam1)), ylab="Residuals", xlab="Fitted")

plot(fitted(whale.gam1), (subdata$TW-fitted(whale.gam1))/sqrt(fitted(whale.gam1)), ylab="Residuals", xlab="Fitted")


whale.gam2 <- gam(TW~s(std.av100,bs="ts")+s(std.chla,bs="ts")+s(std.sbathy,bs="ts")+s(std.av300,bs="ts")+ s(std.Dice,bs="ts")+s(std.dinswb,bs="ts"),
    data=subdata, family="quasipoisson")
summary(whale.gam2)
plot(subdata$TW, fitted(whale.gam2), ylab="Fitted",xlab="Observed")
abline(0,1)

par(mfrow=c(2,1), mgp=c(1.5,0.5,0))
plot(fitted(whale.gam2), (subdata$TW-fitted(whale.gam2)), ylab="Residuals", xlab="Fitted")
plot(fitted(whale.gam2), (subdata$TW-fitted(whale.gam2))/sqrt(fitted(whale.gam2)), ylab="Residuals", xlab="Fitted")


par(mfrow=c(3,2), mar=c(4,4,0.5,0.5))
plot(whale.gam1, scale=0, pages=0, select=1, xlab="Backscatter 25-100m", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam1, scale=0, pages=0, select=2, xlab="Chlorophyll a", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam1, scale=0, pages=0, select=3, xlab="Bathymetry slope", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam1, scale=0, pages=0, select=4, xlab="Backscatter 100-300m", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam1, scale=0, pages=0, select=5, xlab="Dist. to ice edge", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam1, scale=0, pages=0, select=6, xlab="Dist. to inner shelf water", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)

par(mfrow=c(3,2), mar=c(4,4,0.5,0.5))
plot(whale.gam2, scale=0, pages=0, select=1, xlab="Backscatter 25-100m", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam2, scale=0, pages=0, select=2, xlab="Chlorophyll a", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam2, scale=0, pages=0, select=3, xlab="Bathymetry slope", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam2, scale=0, pages=0, select=4, xlab="Backscatter 100-300m", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam2, scale=0, pages=0, select=5, xlab="Dist. to ice edge", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)
plot(whale.gam2, scale=0, pages=0, select=6, xlab="Dist. to inner shelf water", ylab="f(x)", residuals=T, shade=T, lwd=2, pch=1, cex=0.5)

plot(TW~std.av100, data=subdata)
plot(TW~std.av300, data=subdata)
