source("FrontMatter.R")

################
## Chapter 11 ##
################
plotDIRch11 <- paste (plotDIR, "chapter11", "figures", sep="/")
packages(bootstrap)
packages(rv)
packages(maptools)
packages(maps)
packages(mapproj)
packages(rpart)

## EDA

## Reading data and data clean up

  ## non-detect (nd, bd, below detect) are replaced with 0
  ## secchi >x is replaced with x (x = 7 or 5)

eriedata <- read.csv(paste(dataDIR, "ErieData.csv", sep="/"), header=T)
eriedata <- eriedata[, substring(names(eriedata), 1, 1)!="X"]

  ### using maps:
my.box<-function(xlim, ylim, ...){
    segments(x0=xlim, y0=rep(ylim[1],2), x1=xlim, y1=rep(ylim[2], 2), ...)
    segments(y0=ylim, x0=rep(xlim[1],2), y1=ylim, x1=rep(xlim[2], 2), ...)
}

plot(Latitude~Longitude, data=eriedata)
## Five data entry errors:
##  1. Latitude for WE8 (10/21/13) was 875.7333 --
##     replaced mean latitude of other WE8 with Latitude < 48
##  2. Latitude for WE8 (8/13/12) was 49.8369 --
##     replaced mean latitude of other WE8 with Latitude < 48
##  3. Longitude for WE4 (10/15/14) was 88.1940 -- should be 83.1940?
##  4. Latitude for WE7 (7/6/10) was 40.7649 -- should be 41.6749?
##  5. Latitude for WE 2 (5/15/12) was 41.0127 -- 41.7627?

eriedata$Longitude[eriedata$Longitude>87] <-
    eriedata$Longitude[eriedata$Longitude>87]-5
eriedata$Latitude[eriedata$Latitude>48] <-
    mean(eriedata$Latitude[eriedata$Station=="WE8" &
                               eriedata$Latitude < 48])
eriedata$Latitude[eriedata$Latitude<41.1 &
                      eriedata$Station == "WE2"] <- 41.7622
eriedata$Latitude[eriedata$Latitude<41.1 &
                      eriedata$Station == "WE7"] <- 41.6749
eriedata$Rdate<-as.Date(eriedata$Date, format="%m/%d/%y")
## format(eriedata$Rdate, format="%a %A %b %j")
## cool: http://www.r-bloggers.com/date-formats-in-r/
eriedata$Year <- as.POSIXlt(eriedata$Rdate)[[6]]+1900
eriedata$Month <- as.POSIXlt(eriedata$Rdate)[[5]]+1

eriedata$Longitude <- -eriedata$Longitude
erieLOC <- eriedata[,c("Latitude","Longitude")]
coordinates(erieLOC) <- c("Longitude","Latitude")


tikz(file=paste(plotDIRch11, "sampleLOC.tex", sep="/"),
     height=7, width=7,standAlone=F)
par(mar=rep(0, 4))
map("usa", fill=TRUE, col="grey80", xlim=c(-83.5,-82.5),
    ylim=c(41.4, 42.1))
plot(erieLOC, pch=2, col="blue", add=T)

maplocs <- map(projection="sp_mercator", wrap=TRUE, lwd=0.1,
               col="grey", xlim=c(-180, 0),
               interior=FALSE, orientation=c(90, 180, 0), add=TRUE,
               plot=FALSE)
xrange <- range(maplocs$x, na.rm=TRUE)
yrange <- range(maplocs$y, na.rm=TRUE)
aspect <- abs(diff(yrange))/abs(diff(xrange))
                                        # customised to 6.5 by 4.5 figure size
par(fig=c(0.5, 0.99, 0.99 - 0.5*aspect*4.5/6.5, 0.99),
    mar=rep(0, 4), new=TRUE)
plot.new()
plot.window(xlim=c(1,2.00),
            ylim=c(0.45,1))
map(projection="sp_mercator", wrap=TRUE, lwd=0.25, fill=F,
    col=gray(0.25), interior=TRUE, orientation=c(90, 180, 0),
    add=TRUE)
my.box(xlim=c(1.7-0.015,1.725-0.015), ylim=c(0.79, 0.81))
dev.off()
                                        #+END_SRC

## Exploratory plots

tikz(paste(plotDIRch11, "pairs1.tex", sep="/"), 
     height=7, width=7, standAlone=F)
pairs(pMC~dMC+TP+SRP+I(NH4+NO3)+TDP+chla, data=eriedata)
dev.off()
##  tikz(paste(plotDIRch11, "pairs2.tex", sep="/"), 
##       height=7, width=7, standAlone=T)
pdf(paste(plotDIRch11, "pairs2.tex", sep="/"), 
    height=7, width=7)
pairs(pMC~TP+SRP+I(NH4/1000+NO3)+chla+Temp+Secchi+POC+NH4+NO3,
      data=eriedata,log="xy")
dev.off()

tikz(file=paste(plotDIRch11, "MCPOCvN.tex", sep="/"),
     height=5, width=5, standAlone=F)
par(mfrow=c(2,2), mar=c(3,3,1,0.5), mgp=c(1.5,0.125,0), las=1, tck=0.01)
plot(pMC ~ NH4, log="xy", data=eriedata, xlab="NH$_4$",
     ylab="MC")
plot(pMC ~ NO3, log="xy", data=eriedata, xlab="NO$_3$",
     ylab="MC")
plot(POC ~ NH4, log="xy", data=eriedata, xlab="NH$_4$",
     ylab="POC")
plot(POC ~ NO3, log="xy", data=eriedata, xlab="NO$_3$",
     ylab="POC")
dev.off()

  tikz(file=paste(plotDIRch11, "TempMonth.tex", sep="/"),
       height=5, width=6, standAlone=F)
  xyplot(Temp ~ jitter(Month), data=eriedata, 
         xlab="Month", ylab="Temperature ($^\\circ$C)")
  dev.off()

  tikz(file=paste(plotDIRch11, "pMCvMonth.tex", sep="/"),
       height=5, width=6, standAlone=F)
  xyplot(log(pMC) ~ jitter(Month), data=eriedata,
         xlab="Month", ylab="log(MC)")
  dev.off()

  tikz(file=paste(plotDIRch11, "POCvMonth.tex", sep="/"),
       height=5, width=6, standAlone=F)
  xyplot(log(POC) ~ jitter(Month), data=eriedata,
         xlab="Month", ylab="log(POC)")
  dev.off()
  ## explore within month temperature effects

## Using CART:
  pMCrpart <- rpart(pMC~dMC+SRP+NH4+NO3+PON+I(NH4+NO3+PON)+
                        TDP+chla+Temp+Secchi+TP+POC+DOC+Turbidity,
                    data=eriedata)
  tikz(file=paste(plotDIRch11, "pMCcart.tex", sep="/"),
       height=5, width=6, standAlone=F)
  plot(pMCrpart, margin=0.1)
  text(pMCrpart, use.n=T)
  dev.off()

  tikz(file=paste(plotDIRch11, "pMCvPOC.tex", sep="/"),
       height=5, width=6, standAlone=F)
  plot(pMC~ POC, log="xy", data=eriedata)
  abline(v=c(2,5), col="gray")
  dev.off()

  tikz(file=paste(plotDIRch11, "SDvPOC.tex", sep="/"),
       height=5, width=6, standAlone=F)
  plot(POC~Secchi, log="xy",data=eriedata)
  dev.off()

## Conditional plots
NH4.int <- equal.count(log(eriedata$NH4), 4, 0.25)
NO3.int <- equal.count(log(eriedata$NO3), 4, 0.25)
POC.int <- equal.count(log(eriedata$POC), 9, 0.)
## first look

tikz(file=paste(plotDIRch11, "pMCNH4cPOC.tex", sep="/"),
     height=7, width=7, standAlone=F)
xyplot(log(pMC) ~ log(NH4)|POC.int, data=eriedata)
dev.off()

xyplot(log(pMC) ~ log(NO3)|POC.int, data=eriedata)
xyplot(log(pMC) ~ log(NH4/1000+NO3)|POC.int, data=eriedata)

xyplot(log(pMC) ~ log(NH4)|POC.int*factor(Month), data=eriedata)
xyplot(log(pMC) ~ log(NH4)|factor(Month), data=eriedata)

tikz(file=paste(plotDIRch11, "pMCvNH4cYear.tex", sep="/"),
       height=6.5, width=6.5, standAlone=F)
xyplot(log(pMC) ~ log(NH4)|factor(Year), data=eriedata,
       subset=POC>3)
dev.off()

tikz(file=paste(plotDIRch11, "pMCvNH4cMonth.tex", sep="/"),
     height=6.5, width=6.5, standAlone=F)
xyplot(log(pMC) ~ log(NH4)|factor(Month), data=eriedata,
       subset=POC>3)
dev.off()

summer <- eriedata[eriedata$Month==7 | eriedata$Month==8 |
                       eriedata$Month==9,]
NH4.int <- equal.count(log(summer$NH4), 4, 0.25)
  NO3.int <- equal.count(log(summer$NO3), 4, 0.25)
POC.int <- equal.count(log(summer$POC), 9, 0.)
## first look
xyplot(log(pMC) ~ log(NH4)|POC.int, data=summer)
xyplot(log(pMC) ~ log(NO3)|POC.int, data=summer)
pdf(file=paste(plotDIR, "MCvNH4cond.pdf", sep="/"),
    height=6, width=8.5)
xyplot(log(pMC) ~ log(NH4)|POC.int*factor(Month), data=summer)
dev.off()
  xyplot(log(pMC) ~ log(NH4)|factor(Month), data=summer)
## pay attention to September and August

pdf(file=paste(plotDIRch11, "POCvCHLA.pdf", sep="/"), height=5, width=6)
  plot(POC~chla, log="xy",data=eriedata)
  dev.off()

  pdf(file=paste(plotDIRch11, "POCvTP.pdf", sep="/"), height=5, width=6)
  plot(POC~TP, log="xy",data=eriedata)
  dev.off()

  pdf(file=paste(plotDIRch11, "POCvSRP.pdf", sep="/"), height=5, width=6)
  plot(POC~SRP, log="xy",data=eriedata)
  dev.off()

  pdf(file=paste(plotDIRch11, "POCvTPTN.pdf", sep="/"), height=5, width=6)
  pairs(POC~TP+NH4+NO3+PON+I(NH4/1000+NO3)+chla, log="xy",data=eriedata)
  dev.off()

  POCrpart1 <- rpart(POC~SRP+NH4+NO3+I(NH4+NO3)+TDP+chla+
                         Temp+Secchi+TP+DOC+Turbidity,
                     data=eriedata)
  pdf(file=paste(plotDIRch11, "POCcart1.pdf", sep="/"), height=5, width=6)
  plot(POCrpart1, margin=0.1)
  text(POCrpart1, use.n=T)
  dev.off()
  POCrpart2 <- rpart(POC~SRP+NH4+NO3+I(NH4+NO3)+TDP+chla+Temp+TP+DOC,
                     data=eriedata)
  POCrpart3 <- rpart(log(POC)~SRP+NH4+NO3+I(NH4+NO3)+TDP+chla+Temp+TP+DOC,
                     data=eriedata)
  tikz(file=paste(plotDIRch11, "POCcart2.tex", sep="/"),
       height=5, width=6, standAlone=F)
  plot(POCrpart2, margin=0.1)
  text(POCrpart2, use.n=T)
  dev.off()
  summary(POCrpart2)

POCrpart3.p <- prune(POCrpart3, cp=0.031)
tikz(file=paste(plotDIRch11, "POCcart3.tex", sep="/"),
       height=5, width=6, standAlone=F)
  plot(POCrpart3.p, margin=0.1)
  text(POCrpart3.p, use.n=T)
  dev.off()
  summary(POCrpart3)

  tikz(file=paste(plotDIRch11, "TPvMonth.tex", sep="/"),
       height=5, width=6, standAlone=F)
  xyplot(log(TP) ~ jitter(Month), data=eriedata,
         xlab="Month", ylab="log(TP)")
  dev.off()

## Spatial Variation (Focusing on main sites)

  sort(temp <- table(eriedata$Station))
  ## keep stations with at least 30 data points:
  eriedata$Station <- as.character(eriedata$Station)
  keep <- rep(temp >=30, temp)
  eriedata2 <- eriedata[order(eriedata$Station),][keep,]


  erieLOC2 <- eriedata2[,c("Latitude","Longitude")]
  coordinates(erieLOC2) <- c("Longitude","Latitude")

  tikz(file=paste(plotDIRch11, "sampleLOC2.tex", sep="/"), 
       height=7, width=7, standAlone=T)
  par(mar=rep(0, 4))
  map("usa", fill=TRUE, col="grey80", xlim=c(-83.5,-82.5), 
       ylim=c(41.4, 42.1))
  plot(erieLOC2, pch=2, col="blue", add=T)

  maplocs <- map(projection="sp_mercator", wrap=TRUE, lwd=0.1,
                 col="grey", xlim=c(-180, 0),
                 interior=FALSE, orientation=c(90, 180, 0), add=TRUE,
                 plot=FALSE)
  xrange <- range(maplocs$x, na.rm=TRUE)
  yrange <- range(maplocs$y, na.rm=TRUE)
  aspect <- abs(diff(yrange))/abs(diff(xrange))
  # customised to 6.5 by 4.5 figure size
  par(fig=c(0.5, 0.99, 0.99 - 0.5*aspect*4.5/6.5, 0.99),
      mar=rep(0, 4), new=TRUE)
  plot.new()
  plot.window(xlim=c(1,2.00),
              ylim=c(0.45,1))
  map(projection="sp_mercator", wrap=TRUE, lwd=0.25, fill=F, col=gray(0.25),
      interior=TRUE, orientation=c(90, 180, 0), add=TRUE)
  my.box(xlim=c(1.7-0.015,1.725-0.015), ylim=c(0.79, 0.81))
  dev.off()

  pdf(paste(plotDIRch11, "MCvPOCbySite.pdf", sep="/"), height=6, width=6)
  xyplot(log(pMC) ~ log(POC)|Station, data=eriedata2)
  dev.off()
  pdf(paste(plotDIRch11, "POCvTPbySite.pdf", sep="/"), height=6, width=6)
  xyplot(log(POC) ~ log(TP)|Station, data=eriedata2)
  dev.off()

  pdf(paste(plotDIRch11, "MCvPOCbyMonth.pdf", sep="/"), height=5, width=6)
  xyplot(log(pMC) ~ log(POC)|factor(Month), data=eriedata2)
dev.off()

## make a point here:
  pdf(paste(plotDIRch11, "POCvTPbyMonth.pdf", sep="/"), height=6, width=6)
  xyplot(log(POC) ~ log(TP)|Month, data=eriedata2)
  dev.off()

  pdf(paste(plotDIRch11, "MCvPOCbyYear.pdf", sep="/"), height=5, width=6)
  xyplot(log(pMC) ~ log(POC)|factor(Year), data=eriedata2)
  dev.off()
  pdf(paste(plotDIRch11, "POCvTPbyYear.pdf", sep="/"), height=5, width=6)
  xyplot(log(POC) ~ log(TP)|factor(Year), data=eriedata2)
  dev.off()

  pdf(paste(plotDIRch11, "MCvPOCbyYearSep.pdf", sep="/"), height=5, width=6)
  xyplot(log(pMC) ~ log(POC)|factor(Year), data=eriedata2, 
         subset=Month==9, main="September")
  dev.off()
  pdf(paste(plotDIRch11, "MCvPOCbyYearAug.pdf", sep="/"), height=5, width=6)
  xyplot(log(pMC) ~ log(POC)|factor(Year), data=eriedata2, 
         subset=Month==8, main="August")
  dev.off()

  pdf(paste(plotDIRch11, "POCvTPbyYearSep.pdf", sep="/"), height=5, width=6)
  xyplot(log(POC) ~ log(TP)|factor(Year), data=eriedata2, 
         subset=Month==9, main="September")
  dev.off()
  pdf(paste(plotDIR, "POCvTPbyYearAug.pdf", sep="/"), height=5, width=6)
  xyplot(log(POC) ~ log(TP)|factor(Year), data=eriedata2, 
         subset=Month==8, main="August")
  dev.off()
