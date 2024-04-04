library(tidyverse)
tmp <- read_csv("//Users/liujiaqi/Downloads/Math 297/tmp-20160413.csv")
tmp = as.data.frame(tmp)
names(tmp) = c("xx","yy")
dat = tmp
attach(dat)
n = 973

regressores <- function(x, y, nbins){
  xy <- data.frame(x=x, y=y)
  xy <- xy[order(xy$x),]
  z <- cut(xy$x,breaks=seq(min(xy$x),max(xy$x),length=round(nbins)+1),
           labels=1:round(nbins),include.lowest=TRUE)
  xyz <- data.frame(xy,z=z)
  MEANS <- c(by(xyz$y,xyz$z,FUN=mean))
  new = split(xyz, xyz$z)
  res = 0
  for (i in 1:round(nbins)){
    res = res + sum((new[[i]][["y"]]-MEANS[i])^2)
  }
  return(res)
}

k1 = round(n/log(n))
sigmahat = regressores(xx, yy, k1)/(n-k1)
regressrisk = rep(0, times = k1)
for (i in 1:k1){
  regressrisk[i] = regressores(xx, yy, i)/n +2*sigmahat*i/n
}
optimal = which.min(regressrisk)