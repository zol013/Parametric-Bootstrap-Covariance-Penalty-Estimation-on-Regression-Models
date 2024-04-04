library(tidyverse)
require(lspline)
require(splines)
require(MASS)
tmp <- read_csv("//Users/liujiaqi/Library/Mobile Documents/com~apple~CloudDocs/Math 297/tmp-20160413.csv")
tmp = as.data.frame(tmp)
names(tmp) = c("xx","yy")
dat = tmp
attach(dat)
n = length(xx)

#Regressogram
xy <- data.frame(x=xx, y=yy)
xy <- xy[order(xy$x),]
z <- cut(xy$x,breaks=seq(min(xy$x),max(xy$x),length=42+1),
           labels=1:42,include.lowest=TRUE)
xyz <- data.frame(xy,z=z)
MEANS <- c(by(xyz$y,xyz$z,FUN=mean))
newz = MEANS[z]
par(mar = c(3, 3, 3, 3)) 
plot(dat, main = 'Regressogram', xlab = 'x', ylab = 'y')
lines(xy$x, newz,col = "red",lwd = 2)
legend("topleft", legend=c("Data", "Regressogram fit"), col=c("black", "red"),lty =c(NA, 1),pch = c(1, NA), cex=0.8)

  

# Linear spline
k=5
knots = seq(min(xx) ,max(xx), length=k+2)
knots = knots[2:k+1]
spmodel = lm(yy ~ lspline(xx, knots), data=dat)
xlims<-range(xx)
#Generating Test Data
x.grid<-seq(from=xlims[1], to = xlims[2], length = 1000)
newdat = data.frame(x.grid)
names(newdat) = "xx"

plot(dat, main = 'Linear spline', xlab = 'x', ylab = 'y')
points(x.grid,predict(spmodel,newdat),col="red",lwd=2,type="l")
#Plotting the Regression Line to the scatterplot   
#adding cutpoints
abline(v=knots,lty=2,col="red")
legend("topleft", legend=c("Data", "Linear spline"), col=c("black", "red"),lty =c(NA, 1),pch = c(1, NA), cex=0.8)


# Kernel regression
sortdat = dat[order(xx),]
knmodel = ksmooth(sortdat$xx, sortdat$yy, 'normal', bandwidth = 2, x.points=sortdat$xx)
plot(dat, main = 'Kernel regression', xlab = 'x', ylab = 'y')
lines(knmodel, lwd = 2, col = 2)
legend("topleft", legend=c("Data", "Kernel regression"), col=c("black", "red"),lty =c(NA, 1),pch = c(1, NA), cex=0.8)


#Lowess
lowessmd = loess(yy~xx, dat, span = 0.1)
plot(dat, main = 'Lowess', xlab = 'x', ylab = 'y')
j = order(dat$xx)
lines(dat$x[j], lowessmd$fitted[j],col='red', lwd=2)
legend("topleft", legend=c("Data", "Lowess"), col=c("black", "red"),lty =c(NA, 1),pch = c(1, NA), cex=0.8)

