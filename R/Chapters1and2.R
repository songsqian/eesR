##################
### Chapter 1  ###
##################
## no script
source("FrontMatter.R")

plotDIRch1 <- paste (plotDIR, "chapter1", "figures", sep="/")

##################
### Chapter 2  ###
##################
##
plotDIRch2 <- paste (plotDIR, "chapter2", "figures", sep="/")
4+8*9
a <- 4+8*9
hi <- "hello, world"
3>4
3<5
Logic <- 3<5
mode(hi)
(3<4) + (3>4)
TP <- c(8.91,  4.76, 10.30,  2.32, 12.47,  4.49,  3.11,  9.61,  6.35,
        5.84,  3.30, 12.38,  8.99,  7.79,  7.58,  6.70,  8.13,  5.47,
        5.27,  3.52)
violation <- TP >10
mean(violation)
sum(TP)
length(TP)

sum(TP)/length(TP)

my.mean <- function(x){
  total <- sum(x)
  n <- length(x)
  total/n
}
my.mean(TP)
mean(TP)
example(mean)
my.mean <- function(x)
  return(mean(x, na.rm=T))

everg.precip <- read.table(paste(dataDIR, "EvergladesP.txt",
                                 sep="/"), header=T)

EvergData <- data.frame(Precip = unlist(everg.precip[,-c(1,14)]),
                        Year = rep(everg.precip[,1], 12),
                        Month = rep(1:12, each=86))
head(EvergData)

nn<- dim(everg.precip)
EvergData <- data.frame(Precip = unlist(everg.precip[,-c(1,nn[2])]),
                        Year = rep(everg.precip[,1], nn[2]-2),
                        Month = rep(names(everg.precip[-c(1, nn[2])]), each=nn[1]))

samp1 <- rnorm(10, 2, 0.75)
samp1

set.seed(123)
Viol <- numeric() ## creating an empty numeric vector
 for (i in 1:100000){
  Viol[i] <- sum(rnorm(10, 2, 0.75) > 3) > 1
}
mean(Viol)

set.seed(123)
Viol2 <- numeric() ## creating an empty numeric vector
 for (i in 1:1000){
    Viol2[i] <- sum(rnorm(10, 2, 1) > 3) > 1
}
1-mean(Viol2)

viol.sim <- function(n=10, nsims=1000, mu=2, sigma=0.75, cr=3){
    temp <- numeric()
    for (i in 1:nsims)
        temp[i] <- sum(rnorm(n, mu, sigma) > cr) > 0.1*n
    return(mean(temp))
}

pr.sim1 <- viol.sim(n=24)
pr.sim2 <- viol.sim(n=100)
pr.sim3 <- viol.sim(sigma=1)

## Data Cleaning
## add NOAA data from CILER project

## Data aggregation and Reshaping

TPdata <- matrix(c(20.1, 21.5, 30, 15.2, 31, 12, 20, 25, 19, 11, 14, 21),
                 nrow=4, ncol=3, byrow=T)

site.mean <- apply(TPdata, 1, mean)
day.mean <- apply(TPdata, 2, mean)

TPdataframe <- data.frame(TP=as.vector(TPdata),
                          Day=rep(paste("Day", 1:3, sep=" "), each=4),
                          Site=rep(paste("Site", 1:4, sep=" "), 3))

site.mean <- tapply(TPdataframe$TP, TPdataframe$Site, mean)
day.mean <- tapply(TPdataframe$TP, TPdataframe$Day, mean)

## SPARROW example
GISdata <- data.frame(DWNSTID = c(1,1,2,2,2),
                      Load=c(3,NA, 10, NA, NA),
                      X1=c(10, 14, 20, 40, 10),
                      X2=c(3,5,1,2,3),
                      Z1=c(.2,.7,.4,.3,.2))

oo <- order(GISdata$DWNSTID) ## sort by monitoring site
GISdata <- GISdata[oo,]
Y <- GISdata$Load[!is.na(GISdata$Load)] ## Load data
id <- as.numeric(ordered(GISdata$DWNSTID))
idtbl <- table(id)
ns <- max(idtbl)
nr <- max(id)

my.unstack <- function(X, ID, nc, x.names){
  temp <- tapply(X, ID, FUN = function(x, ns=nc){
                    tt <- as.vector(x)
                    if (length(tt < ns))
                       tt <- c(tt, rep(0, ns-length(tt)))
                    return(tt)})
  temp <- as.data.frame(matrix(unlist(temp), nrow=nr, ncol=ns,
           byrow=T))
  names(temp) <- x.names
  return(temp)
}
X1.names <- paste("X1", 1:ns, sep="_")
X2.names <- paste("X2", 1:ns, sep="_")
Z1.names <- paste("Z1", 1:ns, sep="_")
X1 <- my.unstack(GISdata$X1, GISdata$DWNSTID, nc=ns, X1.names)
X2 <- my.unstack(GISdata$X2, GISdata$DWNSTID, nc=ns, X2.names)
Z1 <- my.unstack(GISdata$Z1, GISdata$DWNSTID, nc=ns, Z1.names)
GISdata_reshaped <- cbind (Y, X1,X2, Z1)


## Date objects
first.date <- as.Date("5/27/2000", format="%m/%d/%Y")
second.date <- as.Date("December 31 2003", format="%B %d %Y")
second.date - first.date

first.d <- strptime("5/27/2000 22:15:00", format="%m/%d/%Y %H:%M:%S")
second.d <- strptime("December 31, 2003, 4:25:00", format="%B %d, %Y, %H:%M:%S")
second.d - first.d

mytime <- data.frame(x = rnorm(100),
                     date=as.Date(round(runif(100)*5000),
                         origin="1970-01-01"))

mytime$Month <- format(mytime$date, "%b")
mytime$weekday <- format(mytime$date, "%a")

mytime$D <- as.POSIXlt(mytime$date)
mytime$D[[7]]+1

mytime$Julian <- as.POSIXlt(mytime$D)[[8]]+1
mytime$month <- as.POSIXlt(mytime$D)[[5]]+1
mytime$Year <- as.POSIXlt(mytime$D)[[6]]+1900

