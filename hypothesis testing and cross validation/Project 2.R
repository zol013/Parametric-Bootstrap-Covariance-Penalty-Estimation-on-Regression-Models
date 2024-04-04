# Import data
library(tidyverse)
require(lspline)
require(splines)
require(MASS)
tmp <- read_csv("//Users/liujiaqi/Library/Mobile Documents/com~apple~CloudDocs/Math 297/Project 1/tmp-20160413.csv")
tmp = as.data.frame(tmp)
names(tmp) = c("xx","yy")
dat = tmp
#attach(dat)
xx = dat$xx
yy = dat$yy
n = length(xx)

#\Section{Parametric bootstrap}
# Just like what we did in Project I, we first consider a linear spline model with a very large number of knots $k_1=n/\log n$. This model will give us the estimated residual variance $\widehat{\sigma}^2$ and serve as the big model for the bootstrap.


# Let $\{\widehat{y}_i^{big}, i=1,...,\Sexpr{n}\}$ be the estimated value under this big elspline model with $k_1$ number of bins. We obtain $\widehat\sigma^2 = \Sexpr{sigmahat}$ and $\{\widehat{y}^{big}_i, i=1,...,\Sexpr{n}\}$. 
k_big = round(n/log(n))
es_model = lm(yy ~ elspline(xx, k_big), data=dat) #bs means B-spline basis
sigmahat = sqrt(sum(resid(es_model)^2)/(n-k_big))
yhat_big = fitted(es_model)

# With the help of $\widehat{\sigma}^2$ and $\widehat{y}^{big}$, we are going to generate $B=100$ psuedo data set. To be more precise, for every pair $(x_i,y_i)$, we generate a large number $B$ of simulated observations according to the approximated normal distribution $N(\widehat{y}^{big}_i,\widehat\sigma^2)$. On each pseudo data set, we are going to estimate risks for each models.
# We use the $\Sexpr{n}\times \Sexpr{B}$ matrix to record the pseudo data set denoted as $(x_i, y_i^{*p})$ for $i=1,...n$ and $p=1,...,B$.
B = 100
pseudo = matrix(0, nrow = n, ncol = B)
for (i in 1:B){
  pseudo[,i] = yhat_big + rnorm(n, 0, sigmahat)
}

# For each pseudo data set $(x_i, y_i^{*p})_{i=1}^{n}$, we will apply Efron's procedure to estimate risks. This requires us to generate $B=100$ simulated observations $(x_i, y_i^{*p,b})$ for $i=1,...n$ and $b=1,...,B$ according to the parametric bootstrap $N(\widehat{y}_i^{*p, big}, \widehat{sigma}^2)$. Here $(\widehat{y}_i^{*p,big})_{i=1}^{n}$ is the fitted value for the pseudo data under the big model and $\widehat{\sigma}$ is the estimated variance as before.
# We first calculate $(\widehat{y}_i^{*p,big})_{i=1}^{n}$
pseudo_hat = matrix(0, nrow = n, ncol = B)
for (i in 1:B){
  pseudo_dat = data.frame(x=xx, y=pseudo[,i])
  cur_esmodel = lm(y ~ elspline(x, k_big), data=pseudo_dat)
  pseudo_hat[,i] = fitted(cur_esmodel)
}

# Now we are able to estimate risks for each model. Using the results in Project I, we consider the regressogram with 42 number of bins, linear spline model with 5 number of knots, kernel model with bandwidth 0.25, and lowess model with bandwidth 0.25. Risks are estimated using Efron's procedure.
# Estimated risk for regressogram with 42 number of bins
regressores <- function(x, y, nbins){
  xy <- data.frame(x=x, y=y)
  xy <- xy[order(xy$x),]
  z <- cut(xy$x,breaks=seq(min(xy$x),max(xy$x),length=round(nbins)+1),
           labels=1:round(nbins),include.lowest=TRUE)
  xyz <- data.frame(xy,z=z)
  MEANS <- c(by(xyz$y,xyz$z,FUN=mean))
  new = split(xyz, xyz$z)
  ans = rep(0, times = n+1)
  ans[1:n] = MEANS[z]
  for (i in 1:round(nbins)){
    ans[n+1] = ans[n+1] + sum((new[[i]][["y"]]-MEANS[i])^2)
  }
  return(ans)
}

k1 = 42
Errhat_reg = rep(0, times = B)
covhat_reg = rep(0, times = B)
boot = matrix(0, nrow = n, ncol = B)
meanboot = rep(0, times = n)
for (i in 1:B){
  err_reg = regressores(xx, pseudo[,i], k1)[n+1]
  for (j in 1:B){
    boot[,j] = pseudo_hat[,i] + rnorm(n, 0, sigmahat)
  }
  for (j in 1:n){
    meanboot[j]= mean(boot[j,])
  }
  for (j in 1:B){
    muhatstar = regressores(xx, boot[,j], k1)[1:n]
    covhat_reg[j] = sum(muhatstar*(boot[,j]-meanboot))/(B-1)
  }
  Errhat_reg[i] = err_reg + 2*sum(covhat_reg)
}

# Estimated risk for linear spline model with 5 number of knots.
k2 = 5
knots = seq(min(xx) ,max(xx), length=k2+2)
knots = knots[2:k2+1]
Errhat_sp = rep(0, times = B)
covhat_sp = rep(0, times = B)
for (i in 1:B){
  pseudo_dat = data.frame(x=xx, y=pseudo[,i])
  spmodel = lm(y ~ lspline(xx, knots), data=pseudo_dat)
  err_sp = sum(resid(spmodel)^2)
  for (j in 1:B){
    boot[,j] = pseudo_hat[,i] + rnorm(n, 0, sigmahat)
  }
  for (j in 1:n){
    meanboot[j]= mean(boot[j,])
  }
  for (j in 1:B){
    newdat = data.frame(x=xx, y=boot[,j])
    spmodel_new = lm(y ~ lspline(x,knots), data=newdat)
    muhatstar = fitted(spmodel_new)
    covhat_sp[j] = sum(muhatstar*(boot[,j]-meanboot))/(B-1)
  }
  Errhat_sp[i] = err_sp + 2*sum(covhat_sp)
}

# Estimated risk for kernel model with bandwidth 0.25
k3 = 0.25
Errhat_ker = rep(0, times = B)
covhat_ker = rep(0, times = B)
for (i in 1:B){
  pseudo_dat = data.frame(x=xx, y=pseudo[,i])
  sort_pseudo_dat = pseudo_dat[order(pseudo_dat$x),]
  knmodel = ksmooth(sort_pseudo_dat$x, sort_pseudo_dat$y, 'normal', bandwidth = k3, x.points=sort_pseudo_dat$x)
  yhat_ker = knmodel$y
  err_ker = sum((yhat_ker-sort_pseudo_dat$y)^2)
  for (j in 1:B){
    boot[,j] = pseudo_hat[,i] + rnorm(n, 0, sigmahat)
  }
  for (j in 1:n){
    meanboot[j]= mean(boot[j,])
  }
  for (j in 1:B){
    newdat = data.frame(x=xx, y=boot[,j], z=meanboot)
    sortnewdat = newdat[order(newdat$x),]
    knmodelnew = ksmooth(sortnewdat$x, sortnewdat$y, 'normal', bandwidth = k3, x.points=sortnewdat$x)
    muhatstar = knmodelnew$y
    covhat_ker[j] = sum(muhatstar*(sortnewdat$y-sortnewdat$z))/(B-1)
  }
  Errhat_ker[i] = err_ker + 2*sum(covhat_ker)
}


# Estimated risk for lowess model with bandwidth 0.25
k4 = 0.25
Errhat_low = rep(0, times = B)
covhat_low = rep(0, times = B)
for (i in 1:B){
  pseudo_dat = data.frame(x=xx, y=pseudo[,i])
  lowessmd = loess(y~x, pseudo_dat, span = k4)
  err_low = sum(lowessmd$residuals^2)
  for (j in 1:B){
    boot[,j] = pseudo_hat[,i] + rnorm(n, 0, sigmahat)
  }
  for (j in 1:n){
    meanboot[j]= mean(boot[j,])
  }
  for (j in 1:B){
    newdat = data.frame(x=xx, y=boot[,j])
    lowessmd_new = loess(y~x, newdat, span = k4)
    muhatstar = fitted(lowessmd_new)
    covhat_low[j] = sum(muhatstar*(boot[,j]-meanboot))/(B-1)
}
  Errhat_low[i] = err_low + 2*sum(covhat_low)
}

# With four sequences of risk estimates under each model, we are going to use the hypothesis testing to diagonose the best model.
# The idea is as follows. Suppose we have two models A and B, and their corresponding sequences of estimates risks $(Err^A_{p})_{p=1}^{B}$ and $(Err^B_{p})_{p=1}^{B}$.
# Then the risk differences between model A and B will be the sequence $Diff^{A,B}_{p=1}^{B}=(Err^A_{p}-Err^B_{p})_{p=1}^{B}$. The null hypothesis is model A is no better than B and model B is no better than B, i.e. $H_0: Err^A=Err^B$.
# For every $\alpha\in (0,1)$, we can find values $a$ and $b$ such that $(1-\alpha)B$ of the estimated risk differences stay in the interval $[a,b]$. Then $[a,b]$ will be our $(1-alpha)*100\%$ confidence interval.
# According to the duality between the confidence interval and the hypothesis testing, the largest $\alpha$ for which $[a,b]$ contains 0 will be the p value for our hypothesis testing. If the p value is very small, we reject the null hypothesis and have good reason to believe that either model A is better than model B or model B is better than model A. 
# If the confidence interval falls in the positive real line, we conclude that model B is better than model A while if the confidence interval falls in the negative real line, we conclude that model A is better than model B.
# We first write a function which has inputs riskA, riskB and threshold, and outputs the better model by calculating p values. Specifically, under certain threshold, output 0 means that the null hypothesis was not rejected, 1 means model A is better and 2 means model B is better.
hypotest_p <- function(riskA, riskB, threshold){
  diff = riskA - riskB
  diff = sort(diff)
  if (diff[1]>0){
    return(2)
  }
  if (diff[B]<0){
    return(1)
  }
  i = 1
  j = B
  while (diff[i+1] < 0){
    i = i + 1
  }
  while (diff[j-1] > 0){
    j = j - 1
  }
  alpha1 = 2*i/B
  alpha2 = 1 - 2*j/B
  alpha = min(alpha1,alpha2)
  if (alpha > threshold){
    return(0)
  }
  if (B - j > i){
    return(2)
  }
  if (B - j < i){
    return(1)
  }
}

p_reg_sp = hypotest_p(Errhat_reg, Errhat_sp, 0.1)
p_reg_ker = hypotest_p(Errhat_reg, Errhat_ker, 0.1)
p_reg_low = hypotest_p(Errhat_reg, Errhat_low, 0.1)
p_sp_ker = hypotest_p(Errhat_sp, Errhat_ker, 0.1)
p_sp_low = hypotest_p(Errhat_sp, Errhat_low, 0.1)
p_ker_low = hypotest_p(Errhat_ker, Errhat_low, 0.1)
# The output is summarized in the following table
#\begin{center}
#\begin{tabular}{|c | c | c | c | c | c |} 
#\hline
#p_reg_sp & p_reg_ker & p_reg_low & p_sp_ker & p_sp_low & p_ker_low \\ 
#\hline
#\Sexpr{p_reg_sp} & \Sexpr{p_reg_ker} & \Sexpr{p_reg_low} & \Sexpr{p_sp_ker} & \Sexpr{p_sp_low} & \Sexpr{p_ker_low}\\
#\hline
#\end{tabular}
#\end{center}

# Therefore, according to the hypothesis testing, at the level of $5\%$, Regressogram is better than all the other three models, linear spline model and lowess are better than kernel regression, and linear spline model is no better than lowess model and lowess model is no better than linear spline model. Therefore, we have good reasons to believe that the regressogram with 42 number of bins is the best regression model.
# We draw the graph of the regressogram with 42 number of bins below.
#\begin{figure}
#\begin{center}
#\caption{Regressogram with 42 number of bins is the best model}
#\includegraphics[scale=0.50]{figure1}
#\end{center}
#\end{figure}

#\Section{Cross Validation}
# In this section, we apply 7-fold cross-validation on each models -- regressograms, linear splines, kernel regression and lowess, to estimate the optimal smoothing parameters.
library(caret)
# Create folds
set.seed(5)
flds = createFolds(dat$yy, k=7, list = TRUE, returnTrain = FALSE)

# Regressogram
predict_err_reg = rep(0, times = 100)
for (j in 1:100){
  predict_err = rep(0, times = 7)
  for (i in 1:7){
    validat = dat[flds[[i]], ]
    traindat = dat[-flds[[i]],]
    xy = data.frame(x = traindat$xx, y = traindat$yy)
    xy = xy[order(xy$x),]
    z = cut(xy$x,breaks=seq(min(dat$xx),max(dat$xx),length=j+1),labels=1:j,include.lowest=TRUE)
    xyz = data.frame(xy,z=z)
    MEANS = c(by(xyz$y,xyz$z,FUN=mean))
    vali_z = cut(validat$xx,breaks=seq(min(dat$xx),max(dat$xx),length=j+1),labels=1:j,include.lowest=TRUE)
    predict_err[i] = sum((MEANS[vali_z]-validat$yy)^2)/7
  }
  predict_err_reg[j] = sum(predict_err)
}

optimal_nbins = which.min(predict_err_reg)
par(mar = c(3, 3, 3, 3)) 
plot(1:100, predict_err_reg,pch=16,xlab = "Number of bins", ylab = "Predicted error")
#The optimal number of bins is $\Sexpr{optimal_nbins}$, which is the same as the one obtained in Project I. Furthermore, the graph of the predicted errors for regressograms with different number of bins has the same trend as the graph of the estimated risks of the regressograms with different number of bins plotted in Project I.
# The fitted curve of the regressogram with the optimal number of bins is as follows.
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
legend("topleft", legend=c("Data", "Regressogram fit"), col=c("black", "red"),lty =c(NA, 1),pch = c(1, NA))



# Linear spline model
predict_err_sp = rep(0, times = 20)
for (j in 1:20){
  predict_err = rep(0, times = 7)
  for (i in 1:7){
    validat = dat[flds[[i]], ]
    traindat = dat[-flds[[i]],]
    knots = seq(min(xx), max(xx), length = j+2)
    knots = knots[2:j+1]
    spmodel = lm(yy ~ lspline(xx, knots), data = traindat)
    vali_pred = predict(spmodel, newdata = validat) 
    predict_err[i] = sum((vali_pred - validat$yy)^2)/7
  }
  predict_err_sp[j] = sum(predict_err)
}

optimal_knot = which.min(predict_err_sp)

par(mar = c(3, 3, 3, 3)) 
plot(1:20, predict_err_sp,pch=16,
     xlab = "Number of knots", ylab = "Predicted error")
#The optimal number of knots is $\Sexpr{optimal_knot}$, which is different from the one obtained using Efron's procedure. However, if we plot the predicted errors for linear spline models with different number of knots, we can see that the predicted error at 5 and 10 are similar. This might indicate that if we choose a different random sampling in the cross validation fold, there is a chance that the optimal number of knots may be 5.
# Furthermore, if we compare this graph with the plot of the estimated errors of the linear spline models with different number of knots plotted in Project I, they have the same trend.
# The fitted curve of the linear spline model with the optimal number of knots is as follows.

knots = seq(min(xx) ,max(xx), length=optimal_knot+2)
knots = knots[2:optimal_knot+1]
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
predict_err_ker = rep(0, times = 20)
h = seq(from = 0.05, to = 1, length.out = 20)
for (j in 1:20){
  predict_err = rep(0, times = 7)
  for (i in 1:7){
    validat = dat[flds[[i]], ]
    sort_validat = validat[order(validat$xx),]
    traindat = dat[-flds[[i]],]
    sort_traindat = traindat[order(traindat$xx),]
    knmodel = ksmooth(sort_traindat$xx, sort_traindat$yy, 'normal', bandwidth = h[j], x.points=sort_validat$xx)
    predict_err[i] = sum((knmodel$y - sort_validat$yy)^2)/7
  }
  predict_err_ker[j] = sum(predict_err)
}

optimal_ker = h[which.min(predict_err_ker)]
par(mar = c(3, 3, 3, 3)) 
plot(h, predict_err_ker,pch=16,xlab = "Bandwidth", ylab = "Predicted error")
# The optimal bandwidth is $\Sexpr{optimal_ker}$. Notice that there isn't a big difference between the prediction errors for 0.1 and 0.15.
sortdat = dat[order(xx),]
knmodel = ksmooth(sortdat$xx, sortdat$yy, 'normal', bandwidth = 0.2, x.points=sortdat$xx)
plot(dat, main = 'Kernel regression', xlab = 'x', ylab = 'y')
lines(knmodel, lwd = 2, col = 2)
legend("topleft", legend=c("Data", "Kernel regression"), col=c("black", "red"),lty =c(NA, 1),pch = c(1, NA), cex=0.8)
# If we look at the graph for the kernel regression with 0.1 bandwidth, it is a littlbe bit over fit. The optimal bandwidth using Efron's procedure is 0.25, which is better compared with the one obtained using cross validation.

# Lowess model
predict_err_low = rep(0, times = 30)
h = seq(from = 0.05, to = 0.95, length.out = 30)
for (j in 1:30){
  predict_err = rep(0, times = 7)
  for (i in 1:7){
    validat = dat[flds[[i]], ]
    traindat = dat[-flds[[i]],]
    lowessmd = loess(yy~xx, traindat, span = h[j], control=loess.control(surface="direct"))
    vali_pred = predict(lowessmd, newdata = validat)
    predict_err[i] = sum((vali_pred - validat$yy)^2)/7
  }
  predict_err_low[j] = sum(predict_err)
}
# In the lowess regression funciton ``loess", we add ``control=loess.control(surface="direct")" to avoid NA returns in prediction.
optimal_span = h[which.min(predict_err_low)]
# The optimal span is $\Sexpr{optimal_span}$, which is similar to the optimal span 0.25 we obtained using Efron's procedure. The plot for predicted errors of lowess model with different spans is as follows.
par(mar = c(3, 3, 3, 3)) 
plot(h, predict_err_low,pch=16, xlab = "Span", ylab = "Predicted error")
# The fitted curve for the optimal loess model is as follows.
lowessmd = loess(yy~xx, dat, span = optimal_span)
plot(dat, main = 'Lowess', xlab = 'x', ylab = 'y')
j = order(dat$xx)
lines(dat$x[j], lowessmd$fitted[j],col='red', lwd=2)
legend("topleft", legend=c("Data", "Lowess"), col=c("black", "red"),lty =c(NA, 1),pch = c(1, NA), cex=0.8)



#\Section{Summary}

#\Section{Gaussian process regression model}
# In this section, we are going to use ``GauPro" package to fit a Gaussian process regression model. We plot the fitted curve estimate in red and the 95% credible intervals in blue. 
library(GauPro)
x = dat$xx
y = dat$yy
gp = GauPro(x, y, parallel = FALSE)
plot(dat, main = 'Gaussian process regression', xlab = 'x', ylab = 'y')
curve(gp$predict(x), add=T, col='red',lwd=2)
curve(gp$predict(x)+1.96*gp$predict(x, se=T)$se, add=T, col='blue', lwd=2,lty='dashed')
curve(gp$predict(x)-1.96*gp$predict(x, se=T)$se, add=T, col='blue', lwd=2,lty='dashed')
legend("topleft", legend=c("Fitted curve", "Boundaries of credible interval"), col=c("red", "blue"),lty =c(1, 2), cex=0.8)
# We see that the Gaussian process regression wroks well and most of data points are within the $95\%$ credible interval.


