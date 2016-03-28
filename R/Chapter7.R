source("FrontMatter.R")
#### CART -- using USGS pesticides data from Anderson et al 1996

## Chapter 7

plotDIRch7 <- paste(plotDIR, "chapter7", "figures", sep="/")

#### Willamette River Data

willamette.data <- read.csv(paste(dataDIR, "willamette.csv", sep="/"), header=T)

willamette.data$Diuron <- "Below MDL"
willamette.data$Diuron[willamette.data$P49300>=7.08] <- "High"
willamette.data$Diuron[willamette.data$P49300<7.08 & willamette.data$P49300>=0.83] <- "Medium"
willamette.data$Diuron[willamette.data$P49300<0.83 & willamette.data$P49300>0.02] <- "Low"
willamette.data$Diuron <- ordered(willamette.data$Diuron, levels=c("Below MDL","Low","Medium","High"))
willamette.data$Diuron[is.na(willamette.data$P49300)] <- NA
willamette.data$Month <- ordered(willamette.data$Month, levels=c("Apr","May","Jul","Oct","Nov"))
Snames <- dimnames(iris3)[[3]]
iris.df <- rbind(iris3[,,1],iris3[,,2],iris3[,,3])
iris.df <- as.data.frame(iris.df)
iris.df$Species <- factor(rep(Snames,rep(50,3)))
packages(rpart)
names(iris.df) <- c("Sepal.L", "Sepal.W", "Petal.L", "Petal.W", "Species")
iris.rpart <- rpart(Species ~ Sepal.L+Sepal.W+Petal.L+Petal.W,
                    data=iris.df, method="class")

## Figure 7.1
##postscript(file=paste(plotDIR, "irisCART.eps", sep="/"), height=3, width=3.4,
##           horizontal=F)
tikz(file=paste(plotDIRch7, "irisCART.tex", sep="/"),
     height=3, width=3.4, standAlone=F)
par(mgp=c(1.5,0.5,0), mar=c(1,1,1,1))
plot(iris.rpart, margin=0.1)
text(iris.rpart)
dev.off()

## Figure 7.2
##postscript(file=paste(plotDIR, "iris2D.eps", sep="/"),
##           height=3, width=3, horizontal=F)
tikz(file=paste(plotDIRch7, "iris2D.tex", sep="/"),
           height=3, width=3, standAlone=F)
trellis.par.set(bw.whitebg())
xyplot(Petal.W~Petal.L,iris.df,groups=Species,
    pch=1:3,col=1,
    panel=function(x,y,groups,...){
     panel.superpose(x,y,groups,...)
     panel.abline(v=2.45,lty=2)
     panel.segments(2.45,1.75,max(x)*2,1.75,lty=2)
#     panel.segments(4.95,min(y)*-2,4.95,1.75,lty=2)
#     panel.segments(2.45,1.65,4.95,1.65,lty=2)
#     panel.segments(4.95,1.55,max(x)*2,1.55,lty=2)
   },xlab="Petal Length",ylab="Petal Width",
  key=list(points=T, columns=3,pch=1:3,col=1,cex=0.5,
    text=list(c("Setosa","Versicolor","Virginica"))))
dev.off()

## Diuron -- 49300: R.62 & P49300
##postscript(file="diuronQplot.eps", width=4, height=4.5, horizontal=F)
tikz(file=paste(plotDIRch7, "diuronQplot.tex", sep="/"),
     width=4, height=4.5, standAlone=F)
par(mar=c(3, 3, 1, 1), mgp=c(1.5,0.125,0), las=1, tck=0.01)
plot(((1:94)-0.5)/94, sort(willamette.data$P49300), log="y", xlab="Quantile",
    ylab=expression(paste("Diuron ","(",mu, "g/L)", sep="")))
abline(h=c(0.83, 7.08), lty=2)
text(x=c(0.2,0.2), y=c(1.15, 9), c("0.83","7.08"))
dev.off()

##postscript(file="diurondata.eps", width=5, height=3.5, horizontal=F)
tikz(file=paste(plotDIRch7, "diurondata.tex", sep="/"),
     width=5, height=3.5, standAlone=F)
trellis.par.set(theme = canonical.theme("postscript", col=FALSE))
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
dotplot(Month ~ jitter(log(P49300)) | Class, data=willamette.data, xlab="Log Diuron")
dev.off()

packages(rpart)

## regression model
diuron.rpart <- rpart(log(P49300) ~ NH4+NO2+TKN+N2.3+TOTP+SRP+BOD+ECOL+FECAL+wtemp+bpres+flow+cond+pH+SSC+SSF+
    Longitude+Latitude+Size+LU.Ag+LU.For+LU.Resid+LU.Other+NumCrops+Month, data=willamette.data)
 plot(diuron.rpart, margin=0.1)
 text(diuron.rpart, pretty=T, use.n=T)


## a note on the differences

set.seed(12345)
diuron.rpart <- rpart(log(P49300) ~ NH4+NO2+TKN+N2.3+TOTP+SRP+BOD+ECOL+FECAL+
#    wtemp+bpres+flow+cond+pH+SSC+SSF+
    Longitude+Latitude+Size+LU.Ag+LU.For+LU.Resid+LU.Other+NumCrops+Month,
    data=willamette.data, control=rpart.control(minsplit=4,cp=0.005))
##postscript(file="diuronCART1.eps", height=7, width=5, horizontal=F)
tikz(file=paste(plotDIRch7, "diuronCART1.tex", sep="/"),
     height=7, width=5, standAlone=F)
 plot(diuron.rpart, margin=0.1)
 text(diuron.rpart, cex=0.5)
dev.off()

 printcp (diuron.rpart)
# summary(diuron.rpart)
##postscript(file="diuronCARTcp.eps", width=5, height=3.75, horizontal=F)
tikz(file=paste(plotDIRch7, "diuronCARTcp.tex", sep="/"),
     width=5, height=3.75, standAlone=F)
par(mar=c(3, 3, 1, 1), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
plotcp(diuron.rpart)
dev.off()

diuron.rpart.prune <- prune(diuron.rpart, cp=0.05)
##postscript(file=paste(base, "diuronCARTprune.eps", sep="/"), height=8, width=6, horizontal=F)
tikz(file=paste(plotDIRch7, "diuronCARTprune.tex", sep="/"),
     height=8, width=6, standAlone=F)
nf <- layout(matrix(c(1,2), nrow=2, ncol=1), 1, c(2,1))
# layout.show(nf)
par(mar=c(0,4,1,2))
plot(diuron.rpart.prune, compress=F, branch=0.4, margin=0.1)
text(diuron.rpart.prune, pretty=T, cex=0.55, use.n=T)
title(main="log diuron Concentration")
par(mar=c(0.5,4,0.5,2))
boxplot(split(predict(diuron.rpart.prune)+resid(diuron.rpart.prune),
    round(predict(diuron.rpart.prune), digits=4)),
        ylab="Diuron Concentrations",
        xlab=" ", axes=F, ylim=log(c(0.01, 50)))
axis(2, at=log(c(0.01, 0.1, 1, 10, 50)), labels=c("0.01","0.1","1","10","50"), las=1)
box()
dev.off()

diuron.rpart.prune2 <- prune(diuron.rpart, cp=0.08)
##postscript(file=paste(base, "diuronCARTprune2.eps", sep="/"), height=8, width=6, horizontal=F)
tikz(file=paste(plotDIRch7, "diuronCARTprune2.tex", sep="/"),
     height=8, width=6, standAlone=F)
nf <- layout(matrix(c(1,2), nrow=2, ncol=1), 1, c(2,1))
# layout.show(nf)
par(mar=c(0,4,1,2))
plot(diuron.rpart.prune2, compress=F, branch=0.4, margin=0.1)
text(diuron.rpart.prune2, pretty=T, cex=0.55, use.n=T)
title(main="log diuron Concentration")
par(mar=c(0.5,4,0.5,2))
boxplot(split(predict(diuron.rpart.prune2)+resid(diuron.rpart.prune2),
    round(predict(diuron.rpart.prune2), digits=4)),
        ylab="Diuron Concentrations",
        xlab=" ", axes=F, ylim=log(c(0.01, 50)))
axis(2, at=log(c(0.01, 0.1, 1, 10, 50)), labels=c("0.01","0.1","1","10","50"), las=1)
box()
dev.off()



## classification model


set.seed(12345)
diuron.rpart2 <- rpart(Diuron ~ NH4+NO2+TKN+N2.3+TOTP+SRP+BOD+ECOL+FECAL+
#    wtemp+bpres+flow+cond+pH+SSC+SSF+
    Longitude+Latitude+Size+LU.Ag+LU.For+LU.Resid+LU.Other+NumCrops+Month,
    data=willamette.data, method="class", parms=list(prior=rep(1/4, 4), split="information"),
    control=rpart.control(minsplit=4,cp=0.005))

##postscript(file="diuronCART2.eps", height=7, width=5, horizontal=F)
tikz(file=paste(plotDIRch7, "diuronCART2.tex", sep="/"),
     height=7, width=5, standAlone=F)
plot(diuron.rpart2, margin=0.1)
text(diuron.rpart2, cex=0.55)
dev.off()
#title(main="log diuron Concentration")

tikz(file=paste(plotDIRch7, "diuronCART2cp.tex", sep="/"),
     width=5, height=3.75, standAlone=F)
par(mgp=c(1.5, 0.5, 0))
plotcp(diuron.rpart2)
dev.off()
### the above is the version used in Qian and Anderson (1999)

## dsicuss the use of prior, the use of different split method, and other issues
## design a series of possible models

# 1. using default ## infor and data frequencies
set.seed(123456)
diuron.rpart3 <- rpart(Diuron ~ NH4+NO2+TKN+N2.3+TOTP+SRP+BOD+ECOL+FECAL+
#    wtemp+bpres+flow+cond+pH+SSC+SSF+
    Longitude+Latitude+Size+LU.Ag+LU.For+LU.Resid+LU.Other+NumCrops+Month,
    data=willamette.data, method="class",
    control=rpart.control(minsplit=4,cp=0.000001))#, parms=list(prior=rep(1/4, 4), split="gini"))
plotcp(diuron.rpart3)
diuron.rpart3.prune <- prune(diuron.rpart3, cp=0.06)
 plot(diuron.rpart3.prune, margin=0.1)
 text(diuron.rpart3.prune, pretty=T, use.n=T)
# 2. default -- changing split method to gini:  gini and data frequencies
set.seed(123456)
diuron.rpart4 <- rpart(Diuron ~ NH4+NO2+TKN+N2.3+TOTP+SRP+BOD+ECOL+FECAL+
#    wtemp+bpres+flow+cond+pH+SSC+SSF+
    Longitude+Latitude+Size+LU.Ag+LU.For+LU.Resid+LU.Other+NumCrops+Month,
    data=willamette.data, method="class",
    control=rpart.control(minsplit=4,cp=0.000001), parms=list(split="gini"))
plotcp(diuron.rpart4)
printcp(diuron.rpart4)
diuron.rpart4.prune <- prune(diuron.rpart4, cp=0.06)
 plot(diuron.rpart4.prune, margin=0.1)
 text(diuron.rpart4.prune, pretty=T, use.n=T)

# use equal prior and gini:
set.seed(123456)
diuron.rpart5 <- rpart(Diuron ~ NH4+NO2+TKN+N2.3+TOTP+SRP+BOD+ECOL+FECAL+
#    wtemp+bpres+flow+cond+pH+SSC+SSF+
    Longitude+Latitude+Size+LU.Ag+LU.For+LU.Resid+LU.Other+NumCrops+Month,
    data=willamette.data, method="class",
    control=rpart.control(minsplit=4,cp=0.000001), parms=list(prior=rep(1/4, 4), split="gini"))
plotcp(diuron.rpart5)
printcp(diuron.rpart5)
diuron.rpart5.prune <- prune(diuron.rpart5, cp=0.05)
 plot(diuron.rpart5.prune, margin=0.1)
 text(diuron.rpart5.prune, pretty=T, use.n=T)

## used first
diuron.rpart2.prune <- prune(diuron.rpart2, cp=0.06)
##postscript(file="diuronCART2prune.eps", width=6, height=6, horizontal=F)
tikz(file=paste(plotDIRch7, "diuronCART2prune.tex", sep="/"),
     width=6, height=6, standAlone=F)
par(mar=c(0,4,1,2))
plot(diuron.rpart2.prune, compress=F, branch=0.4, margin=0.1)
text(diuron.rpart2.prune, pretty=T, cex=0.55, use.n=T)
title(main="Diuron Concentration Level")
dev.off()

##postscript(file="diuron4models.eps", height=5, width=5, horizontal=F)
tikz(file=paste(plotDIRch7, "diuron4models.tex", sep="/"),
     height=5, width=5, standAlone=F)
par(mar=c(0,4,1,2), mfrow=c(2,2))
plot(diuron.rpart2.prune, compress=F, branch=0.4, margin=0.1)
text(diuron.rpart2.prune, pretty=T, cex=0.55, use.n=T)
title(main="Equal Prior \\& Infor", cex=0.75)

plot(diuron.rpart3.prune, compress=F, branch=0.4, margin=0.1)
text(diuron.rpart3.prune, pretty=T, cex=0.55, use.n=T)
title(main="Data Prior \\& Infor", cex=0.75)

plot(diuron.rpart5.prune, compress=F, branch=0.4, margin=0.1)
text(diuron.rpart5.prune, pretty=T, cex=0.55, use.n=T)
title(main="Equal Prior \\& Gini", cex=0.75)

plot(diuron.rpart4.prune, compress=F, branch=0.4, margin=0.1)
text(diuron.rpart4.prune, pretty=T, cex=0.55, use.n=T)
title(main="Data Prior \\& Gini", cex=0.75)

dev.off()

### plotting options
##postscript(file=paste(base, "diuronPlots1.eps", sep="/"), height=4, width=4, horizontal=F)
tikz(file=paste(plotDIRch7, "diuronPlots1.tex", sep="/"),
     height=4, width=4, standAlone=F)
## layout.show(nf)
par(mar=c(1,2,2,1))#, mfrow=c(2,2))
plot(diuron.rpart5.prune,margin=0.1,
    main="a. default")
text(diuron.rpart5.prune, cex=0.75)
dev.off()

##postscript(file=paste(base, "diuronPlots2.eps", sep="/"), height=4, width=4, horizontal=F)
tikz(file=paste(plotDIRch7, "diuronPlots2.tex", sep="/"),
           height=4, width=4, standAlone=F)
plot(diuron.rpart5.prune,uniform=T,branch=0.25,
margin=0.1, main="b. uniform with branching")
text(diuron.rpart5.prune,pretty=1,use.n=T, cex=0.75)
dev.off()

##postscript(file=paste(base, "diuronPlots3.eps", sep="/"), height=6, width=5, horizontal=F)
tikz(file=paste(plotDIRch7, "diuronPlots3.tex", sep="/"),
     height=6, width=5, standAlone=F)
plot(diuron.rpart5.prune,uniform=T,branch=0.,
    margin=0.1,main="c. fancy")
text(diuron.rpart5.prune,pretty=1,
    all=T,use.n=T,fancy=T, cex=0.7, fwidth=0.4, fheight=0.8)
dev.off()

