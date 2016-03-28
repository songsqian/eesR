### Environmental and Ecological Statistics with R ###
### 2nd Ed
### Copyright Song S. Qian
### R scripts
### December 15, 2015
###

## Loading and installing (if not already installed)
##  packages

packages<-function(x, repos="http://cran.r-project.org", ...){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x, repos=repos, ...)
    require(x,character.only=TRUE)
  }
}

os <- "Mac"
if (os=="win"){
    base <- "C:/Users/Song/Google Drive/EESwithR/2ndEd"
} else {
    base <- "~/Google Drive/EESwithR/2ndEd"
}
RHome <- paste(base, "R", sep="/")
dataDIR <- paste(base, "R","data", sep="/")
plotDIR <- paste(base, "chapters", sep="/")

setwd(RHome)

packages(arm)
packages(lattice)
packages(tikzDevice)
source("eesrfuns.r")


