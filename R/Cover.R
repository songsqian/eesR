## Cover ##

source("FrontMatter.R")

##### PM2.5 data #####
pmdata<-read.table(paste (dataDIR, "PM-RAW-DATA.txt", sep="/"),header=TRUE)

pmdata$rain<-pmdata$Precip > 0
pmdata$sin.date<-sin(2*pi*pmdata$date/365)
pmdata$log.wind<-log(pmdata$AvgWind)
pmdata$z.log.wind<-as.vector(scale(pmdata$log.wind))
pmdata$z.temp<-as.vector(scale(pmdata$AvgTemp))
pmdata$log.value<-log(pmdata$value)
pmdata$group1<-rep(c(1,2),1096/2)
pmdata$group2<-c(rep(1:3,365),1)
pmdata$Dates <-  as.Date("2003-01-01") + pmdata$date-1
pmdata$Weekday <- weekdays(pmdata$Dates, abbreviate=T)
pmdata$Month <- ordered(months(pmdata$Dates, T), levels=month.abb)

tikz(file=paste(plotDIR,"cover.tex", sep="/"),
           height=3.75, width=4.75, standAlone=T)
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75

obj1 <- xyplot(log.value~AvgTemp, panel=function(x,y,...){
  panel.grid()
  panel.xyplot(x,y, col=grey(0.65), cex=0.25, ...)
  panel.loess(x,y,span=1, degree=1,col=1,...)
},    scales=list(x=list(cex=0.75, tck=0.2), y=list(cex=0.75, tck=0.2)),
               par.settings=trellis.par.temp,
       data=pmdata, xlab="", ylab="Log PM2.5")

obj2 <- xyplot(log.value~AvgTemp|Month, panel=function(x,y,...){
  panel.grid()
  panel.xyplot(x,y, col=grey(0.65), cex=0.25, ...)
  panel.loess(x,y,span=1, degree=1,col=1,...)
  }, layout=c(12, 1),
       scales=list(y=list(tck=0.2),
         x=list(relation="free", cex=0.3,tck=0.2,
           alternating=c(1,2))),
       ## x-axis relation and font size
  par.settings=trellis.par.temp,
  data=pmdata, xlab="Average Daily Temperature (F)", ylab="Log PM2.5")

print(obj1, position=c(1/4, 0.35, 3/4, 1), more=T)
print(obj2, position=c(0, 0, 1,0.4), more=F)

dev.off()

