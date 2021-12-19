## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=TRUE,results='hold'-------------------------------------------------
fit <- lm(weight ~ height, data=women)
summary(fit)$coef
summary(fit)$r.squared

## -----------------------------------------------------------------------------
par(mfrow=c(1,1)) 
plot(fit)

## -----------------------------------------------------------------------------
curve(dnorm(x,0,1),-5,5,lty=1,add=F,col="blue")
curve(dnorm(x,0,2),-5,5,lty=1,add=T,col="red")
curve(dnorm(x,0,3),-5,5,lty=1,add=T,col="green")
curve(dnorm(x,0,4),-5,5,lty=1,add=T,col="orange")
#legend("topleft", inset=0.01, title="norm curve", c("N(0,1)","N(0,2)","N(0,3)","N(0,4)"),lty=c(1, 2), pch=c(15, 17), col=c("blue", "red","green","orange"))

## ---- results='asis', message=TRUE, echo=TRUE---------------------------------
library(kableExtra)
library(plyr)
table1<-kable(head(mtcars))
print(table1)

## ---- results='asis', message=TRUE, echo=TRUE---------------------------------
x_html <- knitr:: kable(head(mtcars), "html")#generating a table of the head of mtcars.
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = FALSE,font_size=15)#set the size of words in the table

## ---- results='asis', message=TRUE, echo=TRUE---------------------------------
library(RColorBrewer)
library(kableExtra)
n <- 9
mycolors <- brewer.pal(n, "Set1")#use the function of brewer.pal to return a vector with ten different colors,each color is represented by a number with Hexadecimal.
head(mtcars)%>%
kable()%>% 
kable_styling("hover", full_width = T,font_size=15) %>%# set the size of words in the table.
column_spec(3:5, bold = T,color = "#E41A1C", background ="#FFFF33") %>%#chose the  three to five columns to add the colors of words and background.the parameter color represents the color of words. 
row_spec(5:6, bold = T, color = "#A65628", background ="#999999")#chose the  five to six columns to add the colors of words and background.

## ----echo=TRUE,results='hide'-------------------------------------------------
inverse_function<-function(theta,y)
{
  return(sqrt(-2*theta^2*log(y)))
}
y<-runif(1000)
inverse_function(2,y)#theta is 2

## ----echo=TRUE,results='hold'-------------------------------------------------
inverse_function<-function(theta,y)
{
  return(sqrt(-2*theta^2*log(y)))
}
y<-runif(1000)
rd<-inverse_function(2,y)#the random numbers of Rayleigh distribution,theta equals to 2
hist(rd,prob=TRUE,col="yellow")#density histogram of sample with theta equals to 2
x<-seq(0,1,.01)
curve((x/4)*exp(-x^2/8),col="black",add=TRUE)#density curve of density function of Rayleigh distribution with theta equals to 2

## ----echo=TRUE,results='hold'-------------------------------------------------
#set.seed(1234)
n<-1000
p1<-0.25
norm_1<-rnorm(n,0,1)
norm_2<-rnorm(n,3,1)
weight<-sample(c(0,1),n,replace=TRUE,prob=c(p1,1-p1))#generate the weight of norm_1 and norm_2.
mixture<-weight*norm_1+(1-weight)*norm_2#the mixture distribution
hist(mixture,prob=TRUE,col="yellow")
lines(density(mixture),col="black",lty=1)

## ----echo=TRUE,results='hold'-------------------------------------------------
#set.seed(1234)
n<-1000
norm_1<-rnorm(n,0,1)
norm_2<-rnorm(n,3,1)
mixture_function<-function(p1)
{
weight<-sample(c(0,1),n,replace=TRUE,prob=c(p1,1-p1))#generate the weight of norm_1 and norm_2.
mixture<-weight*norm_1+(1-weight)*norm_2#the mixture distribution
hist(mixture,prob=TRUE,col="yellow")
lines(density(mixture),col="black",lty=1)
}
par(mfrow = c(2, 3))
i<-1
p1<-seq(0,1,0.1)
for(i in 1:11)
{
mixture_function(p1[i])
}

## ----echo=TRUE,results='hold'-------------------------------------------------
compound_Poisson_process<-function(lamada,t,shape,rate)#t represents the time of compound Poisson process
{
  Poisson_process<-function(lamada,t)
  {
    y<-rpois(1,t*lamada)
    return(y)
  }
  u<-rgamma(Poisson_process(lamada,t),shape,rate)
  return(sum(u))#sum(u)represents sum of random variables
}
results<-replicate(30000,compound_Poisson_process(2,10,4,5))
x1<-mean(results)
y1<-var(results)
print(c(x1,y1))

## ----echo=TRUE,results='hold'-------------------------------------------------
beta_rd<-function(N,shape,rate)#use the random numbers of gamma distributions to generate the random numbers of beta distribution.
{
  set.seed(1234)
  u<-rgamma(N,shape,rate)
  v<-rgamma(N,shape,rate)
  rd<-u/(u+v)
}
empirical_beta<-function(x,N)
{
  u<-beta_rd(N,3,3)
  v<-as.integer(u<x)
  return(sum(v)/N)
}
X<-seq(0,1,0.1)
u<-numeric(length(X))
for(i in 1:length(X))
{
u[i]<-empirical_beta(X[i],100000)
}
F<-pbeta(X,3,3)
Fn<-u
library(kableExtra)
x_html <- knitr::kable((rbind(X,Fn,F)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)#set the size of words in the table
#data.frame('x'=X,'Fn'=u,'F'=pbeta(X,3,3))

## ----echo=TRUE----------------------------------------------------------------
anti <- function(N,sigma,antithetic = TRUE) 
  {
  #set.seed(1234)
  v <- runif(N)#using antithetic method
  if (!antithetic)
  {
  u <- runif(N)
  Y1<- sqrt(-2*sigma^2*log(1-v))#inverse function method to generate random numbers of Rayleigh(σ)distribution.
  Y2<- sqrt(-2*sigma^2*log(1-u))
  Y<-(Y1+Y2)/2
  }
  else
  {
  u <- 1 - v
  Y1<- sqrt(-2*sigma^2*log(1-v))#inverse function method to generate random numbers of Rayleigh(σ)distribution.
  Y2<- sqrt(-2*sigma^2*log(1-u))
  Y<-(Y1+Y2)/2
  }
  return(Y)
}
sigma<-seq(1,10,1)
sd<-numeric(10)
sd_anti<-numeric(10)
for(i in 1:10)
{
  sd_anti[i]<-sd(anti(10000,sigma[i],TRUE))
}
for(j in 1:10)
{
  sd[j]<-sd(anti(10000,sigma[j],FALSE))
}
reduction<-1-sd_anti/sd
rbind(sigma,sd,sd_anti,reduction)

## -----------------------------------------------------------------------------
curve(x^2*exp(-x^2/2)/sqrt(2*pi),0,3,lty=2,add=F,col="blue",main="importance function")
curve(exp(-x),0,3,lty=2,add=T,col="green")
curve(exp(-x^2/2)/sqrt(2*pi),0,3,lty=2,add=T,col="red")
curve(3*x^2,0,3,lty=2,add=T,col="yellow")
#curve(x^2*exp(-x^2/2)/sqrt(2*pi),0,1,lty=2,add=T,col="blue")
#legend( "topright",inset=.01, title="importance function", c("f","f1","f2","f3"),lty=c(1, 2), pch=c(15, 17), col=c("blue", "green","red","yellow"))

## -----------------------------------------------------------------------------
set.seed(12345)
m<-10000
theta.hat<-numeric(4)
se<-numeric(4)
f<-function(x)
{
  y<-(x^2*exp(-x^2/2)/sqrt(2*pi))*(x>0)*(x<1)
  return(y)
}

x<-runif(m)# using f
f_0<-f(x)#f_0 represents f(x)/1
theta.hat[1]<-mean(f_0)
se[1]<-sd(f_0)

x<-rexp(m,1)#using f1
f_1<-f(x)/exp(-x)#f_1 represents f(x)/f1
theta.hat[2]<-mean(f_1)
se[2]<-sd(f_1)

x<-rnorm(m,0,1)#using f2
f_2<-f(x)/(exp(-x^2/2)/sqrt(2*pi))#f_2 represents f(x)/f2
theta.hat[3]<-mean(f_2)
se[3]<-sd(f_2)

x<-runif(m)#using f3
y<-x^(1/3)# inverse transform method
f_3<-f(y)/(3*(y^2))
theta.hat[4]<-mean(f_3)
se[4]<-sd(f_3)
calculus<-0.5-theta.hat
library(kableExtra)
x_html <- knitr::kable((rbind(theta.hat,calculus,se)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)#set the size of words in the table

## ----results='hold'-----------------------------------------------------------
f<-function(x)
{
  y<-(x^2*exp(-x^2/2)/sqrt(2*pi))*(x>0)*(x<1)
  return(y)
}
y<-integrate(f,0,1)
print(y)

## -----------------------------------------------------------------------------
n<-20#the number of samples of R-square distribution
alpha<-0.05
set.seed(1234)
func<-function(n,alpha)
{
  sp<-rchisq(n,df=2)
  es1<-mean(sp)-sd(sp)*qt(alpha/2,df=n-1)/sqrt(n)
  es2<-mean(sp)-sd(sp)*qt(1-alpha/2,df=n-1)/sqrt(n)
  return(c(es1,es2))
}
level<-replicate(100000,func(n,alpha))
upper<-level[1,]
lower<-level[2,]
interval<-c(mean(lower),mean(upper))
real_alpha<-sum(lower<2&upper>2)/100000
interval
real_alpha

## -----------------------------------------------------------------------------
# chi-square distribution
set.seed(1234)
N<-100
alpha<-0.05
mu_0<-1 # the expectation of χ2(1).
f1<-function(N)
{
sp<-rchisq(N,df=1)
qt<-(mean(sp)-mu_0)*sqrt(N)/sd(sp)
return(qt)
}
p1<-replicate(10000,f1(N))#the p_value of 100000 simulations
x1<-qt(alpha/2,df=N-1)
x2<-qt(1-alpha/2,df=N-1)
p_chisq<-sum(p1<x1|p1>x2)/10000
se.hat1<-sqrt((p_chisq)*(1-(p_chisq))/10000)

# uniform(0,2) distribution
set.seed(1234)
N<-100
alpha<-0.05
mu_0<-1 # the expectation of χ2(1).
f2<-function(N)
{
sp<-runif(N,0,2)
qt<-(mean(sp)-mu_0)*sqrt(N-1)/sd(sp)
return(qt)
}
p2<-replicate(100000,f2(N))#the p_value of 10000 simulations
x1<-qt(alpha/2,df=N-1)
x2<-qt(1-alpha/2,df=N-1)
p_unif<-sum(p2<x1|p2>x2)/100000
se.hat2<-sqrt((p_unif)*(1-(p_unif))/100000)

#exp(1) distribution
set.seed(1234)
N<-100
alpha<-0.05
mu_0<-1 # the expectation of exp(1).
f3<-function(N)
{
sp<-rexp(N)
qt<-(mean(sp)-mu_0)*sqrt(N-1)/sd(sp)
return(qt)
}
p3<-replicate(100000,f3(N))#the p_value of 100000 simulations
x1<-qt(alpha/2,df=N-1)
x2<-qt(1-alpha/2,df=N-1)
p_exp<-sum(p3<x1|p3>x2)/100000
se.hat3<-sqrt((p_exp)*(1-(p_exp))/100000)

#make the table
x_html <- knitr::kable((cbind(p_chisq,p_unif,p_exp)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)
y_html <- knitr::kable((cbind(se.hat1,se.hat2,se.hat3)),"html")
kableExtra::kable_styling(y_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

## -----------------------------------------------------------------------------
library(MASS)
d<-2#dimension of random vectors
n<-10000#numbers of experiment 
N<- c(10, 20, 30, 50,100,500)#numbers of samples
cv<-qchisq(0.95,(d*(d+1)*(d+2)/6))#0.95 quantile of chi-square distribution
p_value <- numeric(length(N))
mn<-function(n,d)#generate random vectors of standard multivariate normal distribution,r represents number of vectors.
{
  mu<-integer(d)
  sigma<-matrix(data=0,d,d)
  sigma<-diag(d)
  y<-mvrnorm(n,mu,sigma)
  return(y)
}
statistic<-function(X)#compute the value of estimator of skewness.
{
  sigma=cov(X,X)
  n=nrow(X)
  mat=X%*%(solve(sigma))%*%(t(X))
  return(sum(mat^3)/(6*n))
}
for(i in 1:length(N))#compute the rate of rejection of test.
{
  sktests<-numeric(n)
  for(j in 1:n)
  {
    x<-scale((mn(N[i],d)),center=T,scale=T)
    sktests[j]<-as.integer(statistic(x)>cv)
    j<-j+1
  }
  p_value[i]<-mean(sktests)
  i<-i+1
}
#make the table
y_html <- knitr::kable((cbind(N,p_value)),"html")
kableExtra::kable_styling(y_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

## -----------------------------------------------------------------------------
library(MASS)
set.seed(1234)
d<-2#dimension of random vectors
n<-5000#numbers of experiment 
k<-250#numbers of samples
epsilon <- c(seq(0,0.2,0.01),seq(0.2,1,0.05))
N <- length(epsilon)
powe<- numeric(N)#power of test
cv<-qchisq(0.95,(d*(d+1)*(d+2)/6))#0.95 quantile of chi-square distribution
statistic<-function(X)#compute the value of estimator of skewness.
{
  sigma=cov(X,X)
  n=nrow(X)
  mat=X%*%(solve(sigma))%*%(t(X))
  return(sum(mat^3)/(6*n))
}
for (j in 1:N)
{ 
weight<-epsilon[j]
sktests<-numeric(n)
for (i in 1:n)
{
sigma<-sample(c(1,10),replace=TRUE,size=k*d,prob=c(1-weight,weight))
x<-array(rnorm(k*d,0,sigma),dim=c(k,d))
y<-scale(x,center=T,scale=T)
sktests[i] <- as.integer((statistic(y))>=cv)
}
powe[j] <- mean(sktests)
}
print(powe)
plot(epsilon, powe, type = "b",
xlab = bquote(epsilon), ylim = c(0,1),main="power and epsilon")
abline(h = 0.05, lty = 3)

## -----------------------------------------------------------------------------
library(bootstrap)
set.seed(1234)
B<-2000 # number of replicates
results<-numeric(B)
ev<-eigen(cov(scor))
theta<-max(ev$values)/sum(ev$values)
statis<-function(X)
{
index<-sample(1:nrow(X),size=nrow(X),replace=TRUE)
sp<-numeric(nrow(X))
for(i in 1:nrow(X))
{
sp<-rbind(sp,X[index[i],])
i<-i+1
}
sp<-sp[-1,]
ev<-eigen(cov(sp))
stat<-max(ev$values)/sum(ev$values)
return(stat)
}
results<-replicate(B,statis(scor))
#mean(results)
b_boot<-mean(results)-theta
sd_boot<-sd(results)
round(c(b_boot,sd_boot),5)

## -----------------------------------------------------------------------------
#Jackknife
ev<-eigen(cov(scor))
theta<-max(ev$values)/sum(ev$values)
theta_star<-numeric(nrow(scor))
for(i in 1:nrow(scor))
{
  X<-scor[-i,]
  ev<-eigen(cov(X))
  theta_star[i]<-max(ev$values)/sum(ev$values)
  i<i+1
}
b_jack<-(nrow(scor)-1)*(mean(theta_star)-theta)
sd_jack<-(nrow(scor)-1)*sqrt(var(theta_star)/nrow(scor))
round(c(b_jack,sd_jack),8)

## -----------------------------------------------------------------------------
library(bootstrap)
set.seed(1234)
B<-2000 # number of replicates
alpha<-0.05#quantiles
results<-numeric(B) 
ev<-eigen(cov(scor))
theta<-max(ev$values)/sum(ev$values)
statis<-function(X)
{
  index<-sample(1:nrow(X),size=nrow(X),replace=TRUE)
  sp<-numeric(nrow(X))
  for(i in 1:nrow(X))
  {
    sp<-rbind(sp,X[index[i],])
    i<-i+1
  }
  sp<-sp[-1,]
  ev<-eigen(cov(sp))
  stat<-max(ev$values)/sum(ev$values)
  return(stat)
}
results<-replicate(B,statis(scor))
results<-sort(results)
lower<-results[B*alpha/2]
upper<-results[B*(1-alpha/2)]
interval_perc<-c(lower,upper)
interval_perc

## -----------------------------------------------------------------------------
library(bootstrap)
library(boot)
set.seed(1234)
B<-2000 # number of replicates
alpha<-0.05#quantiles
results<-numeric(B)
ev<-eigen(cov(scor))
theta<-max(ev$values)/sum(ev$values)
statis<-function(X)
{
  index<-sample(1:nrow(X),size=nrow(X),replace=TRUE)
  sp<-numeric(nrow(X))
  for(i in 1:nrow(X))
  {
    sp<-rbind(sp,X[index[i],])
    i<-i+1
  }
  sp<-sp[-1,]
  ev<-eigen(cov(sp))
  stat<-max(ev$values)/sum(ev$values)
  return(stat)
}
boot_Bca<-function(X)
{
results<-replicate(B,statis(X))
k<-as.integer(results<theta)
z0<-qnorm((sum(k)/B),0,1)
L<-mean(results)-results
a<-sum(L^3)/(6 * sum(L^2)^1.5)
alpha_1<-pnorm(z0+(z0+qnorm(alpha/2))/(1-a*(z0+qnorm(alpha/2))))
alpha_2<-pnorm(z0+(z0+qnorm(1-alpha/2))/(1-a*(z0+qnorm(1-alpha/2))))
results<-sort(results)
lower<-results[B*alpha_1]
upper<-results[B*alpha_2]
interval_bca<-c(lower,upper)
interval_bca
}
boot_Bca(scor)

## -----------------------------------------------------------------------------
set.seed(1234)
stat <- function(data, index)
 {
 x <- data[index,]
 ev <- eigen(cov(x))$values
 theta <- ev[1] / sum(ev)
 return(theta)
 }
 bootstrap_result <- boot(data=scor,statistic = stat, R = B)
 interval_boot<-boot.ci(bootstrap_result, type=c("perc", "bca"))
 y_html <- knitr::kable((cbind('perc'=interval_perc,'bca'=boot_Bca(scor),'boot_perc'=interval_boot$percent[4:5],'boot_bca'=interval_boot$bca[4:5])),"html")
kableExtra::kable_styling(y_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

## -----------------------------------------------------------------------------
#N(0,1)
#set.seed(1)
options(warn=-1)
library(boot)
library(moments)
N<-1000#numbers of experiment
n<-100#numbers of random numbers
B<-2000
normal<-numeric(N)
basic<-numeric(N)
percent<-numeric(N)
func<-function(x,indices)
{
  mu2<-var(x[indices,])
  mu3<-mean((x[indices,]-mean(x[indices,]))^3)
  sk<-mu3/mu2^(3/2)
  return(sk)
}
for(j in 1:N)
{
x<-as.matrix(rnorm(100))
sk<-boot(x,statistic=func, R=B)
ci<-boot.ci(sk,type="all")
normal[j]<-as.integer((ci$normal[2]<0&&ci$normal[3]>0))
basic[j]<-as.integer((ci$basic[4]<0&&ci$basic[5]>0))
percent[j]<-as.integer((ci$percent[4]<0&&ci$percent[5]>0))
j<-j+1
}
normal_rate1<-sum(normal)/N
basic_rate1<-sum(basic)/N
percent_rate1<-sum(percent)/N
y_html <- knitr::kable((cbind('perc'=percent_rate1,'basic'=basic_rate1,'normal'=normal_rate1)),"html")
kableExtra::kable_styling(y_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

## -----------------------------------------------------------------------------
#chi-square
#set.seed(1)
options(warn=-1)
library(boot)
library(moments)
N<-1000#numbers of experiment
n<-100#numbers of random numbers
B<-2000
normal<-numeric(N)
basic<-numeric(N)
percent<-numeric(N)
func<-function(x,indices)
{
  mu2<-var(x[indices,])
  mu3<-mean((x[indices,]-mean(x[indices,]))^3)
  sk<-mu3/mu2^(3/2)
  return(sk)
}
for(j in 1:N)
{
x<-as.matrix(rchisq(100,5))
sk<-boot(x,statistic=func, R=B)
ci<-boot.ci(sk)
normal[j]<-as.integer((ci$normal[2]<2*sqrt(2/5)&&ci$normal[3]>2*sqrt(2/5)))
basic[j]<-as.integer((ci$basic[4]<2*sqrt(2/5)&&ci$basic[5]>2*sqrt(2/5)))
percent[j]<-as.integer((ci$percent[4]<2*sqrt(2/5)&&ci$percent[5]>2*sqrt(2/5)))
j<-j+1
}
normal_rate2<-sum(normal)/N
basic_rate2<-sum(basic)/N
percent_rate2<-sum(percent)/N
y_html <- knitr::kable((cbind('perc'=percent_rate2,'basic'=basic_rate2,'normal'=normal_rate2)),"html")
kableExtra::kable_styling(y_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

## -----------------------------------------------------------------------------
set.seed(0)
data<-iris
R<-999
x <- as.numeric(data[1:50,2])
y <- as.numeric(data[1:50,3])
z<-c(x,y)
N<-length(z)
per<-length(x)
results<-numeric(R)
K<-1:per
t<-(cor.test(x,y))$estimate
for(i in 1:R)
{
  per<-sample(K,replace=FALSE)
  sp_x<-z[per]
  sp_y<-z[-per]
  results[i]<-cor(sp_x,sp_y,method="spearman")
  i<-i+1
}
pvalue<-mean(abs(c(t,results))>abs(t))
round(c(pvalue,(cor.test(x,y))$p.value),3)

## -----------------------------------------------------------------------------
library(Ball)
library(RANN)
library(energy)
library(boot)
set.seed(1234)
##NN test
k<-2
n1<-50
n2<-50
alpha<-0.05
n<-500#number of experiments
NN<-numeric(n)
energy<-numeric(n)
Ball<-numeric(n)
Tn<-function(z, index,length,k)
{
  n1<-length[1] 
  n2<-length[2] 
  n<-n1+n2
  if(is.vector(z))
  z<-data.frame(z,0)
  z<-z[index, ]
  NN_result<-nn2(data=z, k=k+1)
  part_1<- NN_result$nn.idx[1:n1,-1]
  part_2<- NN_result$nn.idx[(n1+1):n,-1]
  result1<- sum(part_1<=n1+0.5)
  result2<- sum(part_2>=n1+0.5)
  return((result1+result2)/(k * n))
}
#N<-c(nrow(x),nrow(y))
NN_test<-function(data,stat,k)
{
  boot_results<-boot(data,Tn,R=999,sim="permutation",length=N,k=2)
  p_value<-mean(c(boot_results$t0,boot_results$t)>boot_results$t0)
  return(p_value)
}
for ( i in 1:n)#equal expectation and unequal variance
{
x <- matrix(rnorm(n1*k,0,1),ncol=k)
y <- matrix(rnorm(n2*k,0,1.5),ncol=k)
z <- rbind(x,y)
N<-c(nrow(x),nrow(y))
NN[i]<-NN_test(z,Tn,k)
energy[i]<-eqdist.etest(z, sizes=N, R=999)$p.value
Ball[i]<-bd.test(x=x, y=y, num.permutations=999)$p.value
i<-i+1
}
NN_pvalue_1<-mean(NN<alpha)
energy_pvalue_1<-mean(energy<alpha)
Ball_pvalue_1<-mean(Ball<alpha)
x_html<-knitr::kable((cbind('NN1'=NN_pvalue_1,'energy_1'=energy_pvalue_1,'Ball_1'=Ball_pvalue_1)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)


for( j in 1:n)#equal expectation and unequal variance
{
x <- matrix(rnorm(n1*k,0,1),ncol=k)
y <- matrix(rnorm(n2*k,1.5,1.5),ncol=k)
z <- rbind(x,y)
N<-c(nrow(x),nrow(y))
NN[i]<-NN_test(z,Tn,k)
energy[i]<-eqdist.etest(z, sizes=N, R=999)$p.value
Ball[i]<-bd.test(x=x, y=y, num.permutations=999)$p.value
i<-i+1
}
NN_pvalue_2<-mean(NN<alpha)
energy_pvalue_2<-mean(energy<alpha)
Ball_pvalue_2<-mean(Ball<alpha)
y_html<-knitr::kable((cbind('NN2'=NN_pvalue_2,'energy_2'=energy_pvalue_2,'Ball_2'=Ball_pvalue_2)),"html")
kableExtra::kable_styling(y_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

p<-0.5
for( m in 1:n)#t(1) and bimodel distribution
{
  x <- rt(n1,1)
  weight<-sample(c(0,1),n2,replace=TRUE,prob=c(p,1-p))
  y <- weight*rnorm(n2,0,1)+(1-weight)*rnorm(n2,0,2)
  z <- c(x,y)
  N<-c(length(x),length(y))
  NN[m]<-NN_test(z,Tn,k)
  energy[m]<-eqdist.etest(z, sizes=N, R=999)$p.value
  Ball[m]<-bd.test(x=x, y=y, num.permutations=999)$p.value
  m<-m+1
}
NN_pvalue_3<-mean(NN<alpha)
energy_pvalue_3<-mean(energy<alpha)
Ball_pvalue_3<-mean(Ball<alpha)
z_html<-knitr::kable((cbind('NN3'=NN_pvalue_3,'energy_3'=energy_pvalue_3,'Ball_3'=Ball_pvalue_3)),"html")
kableExtra::kable_styling(z_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

for( l in 1:n)#Unbalanced samples
{
x <- matrix(rnorm(10*k,0,1),ncol=k)
y <- matrix(rnorm(100*k,1.5,1.5),ncol=k)
z <- rbind(x,y)
N<-c(nrow(x),nrow(y))
NN[i]<-NN_test(z,Tn,k)
energy[i]<-eqdist.etest(z, sizes=N, R=999)$p.value
Ball[i]<-bd.test(x=x, y=y, num.permutations=999)$p.value
i<-i+1
}
NN_pvalue_4<-mean(NN<alpha)
energy_pvalue_4<-mean(energy<alpha)
Ball_pvalue_4<-mean(Ball<alpha)
z_html<-knitr::kable((cbind('NN4'=NN_pvalue_4,'energy_4'=energy_pvalue_4,'Ball_4'=Ball_pvalue_4)),"html")
kableExtra::kable_styling(z_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

## -----------------------------------------------------------------------------
#set.seed(0)
k<-4
N<-20000
burn<<-2000
x0<-c(-10,-5,5,10)
cauthy_chain<-function(sigma,N,mu)#sigma is the parameter of proposal function,N is numbers of chains.
{
set.seed(0)
sd_cauthy<-function(x)
{
 y<-1/(pi*(1+x^2))
 return(y)
}
xt<-numeric(N)
yt<-numeric(N)
xt[1]<-rnorm(1,mu,sigma)
for(i in 2:N)
{
  yt[i]<-rnorm(1,xt[i-1],sigma)
  u<-runif(1)
  r<-pmin(1,((sd_cauthy(yt[i])*dnorm(xt[i-1],yt[i],sigma))/(sd_cauthy(xt[i-1])*dnorm(yt[i],xt[i-1],sigma))))
  if(u<r)
  {
    xt[i]<-yt[i]
  }
  else
  {
    xt[i]<-xt[i-1]
  }
  i<-i+1
}
return(xt)
}
xt<-cauthy_chain(1,N,0)[-(1:1000)]
p<-seq(0.1,0.9,0.1)
deciles<-quantile(xt,probs=p)
decile_real<-qt(seq(0.1,0.9,0.1),1)
x_html<-knitr::kable(round(cbind(decile_real,deciles),3),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

#the function of Gelman.Rubin
Gelman.Rubin<-function(psi)
{
psi<-as.matrix(psi)
n<-ncol(psi)
k<-nrow(psi)
psi.means<-rowMeans(psi)# row means
B<-n*var(psi.means)#between variance est.
psi.w<-apply(psi,1,"var")#within variances
W<-mean(psi.w)#within est.
v.hat<-W*(n-1)/n+(B/n)#upper variance est.
r.hat<-v.hat/W#G-R statistic
return(r.hat)
}
#compute the Gelman-Rubin statistics and make a plot of R
k<-4
N<-20000
burn<<-2000
x0<-c(-10,-5,5,10)
sp<- matrix(0, nrow=k, ncol=N)
for (i in 1:k)
{
sp[i, ] <- cauthy_chain(1,N,x0[i])
i<-i+1
}
psi <- t(apply(sp, 1, cumsum))
for (j in 1:nrow(psi))
{
psi[j,] <- psi[j,] / (1:ncol(psi))
j<-j+1
}
R<-Gelman.Rubin(psi)
rhat_cauthy<-numeric(N)
for(k in(burn+1):N)
{
rhat_cauthy[k]<-Gelman.Rubin(psi[,1:k])
k<-k+1
}
plot(rhat_cauthy[(burn+1):N],type="l",xlab="N",ylab="R")
abline(h=1.2,lty=2)

## -----------------------------------------------------------------------------
set.seed(0)
N<-15000#length of the chain
discard<-1500# length of burn
a<-2
b<-3
n<-100
x0<-c(30,40,50,60)
y0<-c(0.3,0.4,0.5,0.6)
sp<-matrix(0,N,2)
sp[1,]<-c(x0[1],y0[1])
bivarate_chain<-function(a,b,n,x,y)
{
for( i in 2:N)
{
  #x<-sp[i-1,2]
  sp[i,1]<-rbinom(1,n,sp[i-1,2])
  sp[i,2]<-rbeta(1,(sp[i,1]+a),(n-x0+b))
  i<-i+1
}
return(sp)
}
chain<-bivarate_chain(a,b,n,x0,y0)[((discard+1):N),]
x<-colMeans(chain)
x_html<-knitr::kable((cbind('mean_binom'=x[1],'mean_beta'=x[2])),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

#the function of Gelman.Rubin
Gelman.Rubin<-function(psi)
{
psi<-as.matrix(psi)
n<-ncol(psi)
k<-nrow(psi)
psi.means<-rowMeans(psi)# row means
B<-n*var(psi.means)#between variance est.
psi.w<-apply(psi,1,"var")#within variances
W<-mean(psi.w)#within est.
v.hat<-W*(n-1)/n+(B/n)#upper variance est.
r.hat<-v.hat/W#G-R statistic
return(r.hat)
}
bivarate_chain1<-bivarate_chain(a,b,n,x0[1],y0[1])
bivarate_chain2<-bivarate_chain(a,b,n,x0[2],y0[2])
bivarate_chain3<-bivarate_chain(a,b,n,x0[3],y0[3])
bivarate_chain4<-bivarate_chain(a,b,n,x0[4],y0[4])
xt1=rbind(bivarate_chain1[,1],bivarate_chain2[,1] ,bivarate_chain3[,1] ,bivarate_chain4[,1])
xt2=rbind(bivarate_chain1[,2],bivarate_chain2[,2] ,bivarate_chain3[,2],bivarate_chain4[,2])
psi1 <- t(apply(xt1, 1, cumsum))
psi2 <- t(apply(xt2, 1, cumsum))
for (i in 1:nrow(psi1))
{
psi1[i,] <- psi1[i,] / (1:ncol(psi1))
psi2[i,] <- psi2[i,] / (1:ncol(psi2))
i<-i+1
}
print(c(Gelman.Rubin(psi1),Gelman.Rubin(psi2)))
rhat_binomal<-numeric(N)
rhat_beta<-numeric(N)
for(k in(discard+1):N)
{
rhat_binomal[k]<-Gelman.Rubin(psi1[,1:k])
rhat_beta[k]<-Gelman.Rubin(psi2[,1:k])
k<-k+1
}
par(mfrow=c(1,2)) 
plot(rhat_binomal[(discard+1):N],type="l",xlab="N",ylab="R")
abline(h=1.2,lty=2)
plot(rhat_beta[(discard+1):N],type="l",xlab="N",ylab="R")
abline(h=1.2,lty=2)


## -----------------------------------------------------------------------------
a<-c(1,2)
error<-10^(-100)
k<-1
fn<-function(a,d,m)
{
  N<-100
  series<-numeric(N)
  k_norm<-function(a,k)#calculate the Kth norm of a vector
  {
    s<-a^(2)
    knorm<-(sum(s))^(1/(2))
    knorm<-knorm^k
    return(knorm)
  }
  f<-function(a,d,k)#series
  {
    f1<-((-1)^k)/prod(1:k)/2^k
    f2<-k_norm(a,2*k+2)/(2*k+1)/(2*k+2)
    f3<-exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1))
    return(f1*f2*f3)
  }
  series[1]<-(k_norm(a,2))*exp(lgamma((d+1)/2)+lgamma(3/2)-lgamma(d/2+1))/2
  while(abs(series[k])>error)
  {
    series[k+1]<-f(a,length(a),k)
    k<-k+1
  }
series<-series[-(k:N)]
sum<-cumsum(series)
print(sum)
print(series)
}
fn(a,2,10)

## -----------------------------------------------------------------------------
a<-c(1,2)
eplision<-10^(-100)
k<-1
fn<-function(a,d,m)
{
  N<-100
  series<-numeric(N)
  k_norm<-function(a,k)#calculate the Kth norm of a vector
  {
    s<-a^(2)
    knorm<-sum(s)^(1/2)
    knorm<-knorm^k
    return(knorm)
  }
  f<-function(a,d,k)#series
  {
    f1<-((-1)^k)/prod(1:k)/2^k
    f2<-k_norm(a,2*k+2)/(2*k+1)/(2*k+2)
    f3<-exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1))
    return(f1*f2*f3)
  }
  series[1]<-(k_norm(a,2))*exp(lgamma((d+1)/2)+lgamma(3/2)-lgamma(d/2+1))/2
  while(abs(series[k])>eplision)
  {
    series[k+1]<-f(a,length(a),k)
    k<-k+1
  }
series<-series[-(k:N)]
return(series[m])
}
m<-seq(1,15,1)
series<-numeric(length(m))
for( i in 1:length(m))
{
  series[i]<-fn(a,length(a),m[i])
  i<-i+1
}
sum<-cumsum(series)
x_html<-knitr::kable(rbind(sum,series),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

## -----------------------------------------------------------------------------
k<-c(4,10,15,20,25)
#t-distribution
ck<-function(a,k)
{
  y<-sqrt(a^2*k/(k+1-a^2))
  return(y)
}
equation_k1<-function(a)
{
 eq<-pt(ck(a,k[1]-1),k[1]-1)-pt(ck(a,k[1]),k[1])
 return(eq)
}
equation_k2<-function(a)
{
 eq<-pt(ck(a,k[2]-1),k[2]-1)-pt(ck(a,k[2]),k[2])
 return(eq)
}
equation_k3<-function(a)
{
 eq<-pt(ck(a,k[3]-1),k[3]-1)-pt(ck(a,k[3]),k[3])
 return(eq)
}
equation_k4<-function(a)
{
 eq<-pt(ck(a,k[4]-1),k[4]-1)-pt(ck(a,k[4]),k[4])
 return(eq)
}
equation_k5<-function(a)
{
 eq<-pt(ck(a,k[5]-1),k[5]-1)-pt(ck(a,k[5]),k[5])
 return(eq)
}
solution_k1=uniroot(equation_k1,c(0.01,sqrt(k[1])),tol=1e-10)$root
solution_k2=uniroot(equation_k2,c(0.01,sqrt(k[2])-0.1),tol=1e-10)$root
solution_k3=uniroot(equation_k3,c(0.01,sqrt(k[3])-0.1),tol=1e-10)$root
solution_k4=uniroot(equation_k4,c(0.01,sqrt(k[4]-0.1)),tol=1e-10)$root
solution_k5=uniroot(equation_k5,c(0.01,sqrt(k[5]-1)),tol=1e-10)$root
x_html<-knitr::kable(rbind(sum,series),"html")
x_html<-knitr::kable((cbind('k=4'=solution_k1,'k=10'=solution_k2,'k=15'=solution_k3,'k=20'=solution_k4,'k=25'=solution_k5)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

## -----------------------------------------------------------------------------
#MLE
sp<-c(0.54,0.48,0.33,0.43,1.00,1.00,0.91,1.00,0.21,0.85)
n1<-sum(as.integer(sp<1))
n2<-sum(as.integer(sp==1))
mlogL<-function(theta=1)
{#minus log-likelihood of exp.density,rate1/theta
y<-(-(sum(sp)/theta+n1*log(theta)))
return(-y)
}
library(stats4)
fit<-mle(mlogL)
summary(fit)
#EM
ep1<-10^(-6)
k<-1
lambda<-numeric(50)
lambda[1]<-1
lambda[2]<-(sum(sp)+n2*lambda[1])/length(sp)
while(abs(lambda[k+1]-lambda[k])>10^(-6))
{
  lambda[k+2]=(sum(sp)+n2*lambda[k+1])/length(sp)
  k<-k+1
}
print(c(lambda[k+1],k))

## -----------------------------------------------------------------------------
set.seed(1234)
trims<-c(0,0.1,0.2,0.5)
x<-rcauchy(100)
lapply(trims,function(trim)mean(x,trim=trim))
lapply(trims, mean,x=x)

## -----------------------------------------------------------------------------
#exercise 3
#lapply
formulas<-list(mpg ~disp,mpg ~I(1/disp),mpg ~disp +wt,mpg ~I(1/disp) +wt)
lapply(formulas,lm,mtcars)
# for loop
results<-vector(mode="list",length=length(formulas))
for ( i in 1:length(formulas))
{
  results[[i]]<-lm(formulas[[i]],data=mtcars)
  i<-i+1
}
results

#exercise 4
bootstraps<-lapply(1:10, function(i){
  rows<-sample(1:nrow(mtcars),rep =TRUE)
  mtcars[rows, ]})
#lapply
lapply(bootstraps,lm,mtcars)
# for loop
res<-vector(mode="list",length=length(bootstraps))
for ( i in 1:length(bootstraps))
{
  res[[i]]<-lm(bootstraps[[i]],data=mtcars)
  i<-i+1
}
rsq<-function(mod)summary(mod)$r.squared
reg<-cbind(c(results,res))
rsquare<-lapply(reg,rsq)
x_html<-knitr::kable((rbind('mpg~disp[1]'=rsquare[[1]],'mpg ~I(1/disp)'=rsquare[[2]],'mpg ~disp +wt'=rsquare[[3]],'mpg ~I(1/disp) +wt'=rsquare[[4]],'mpg~disp[2]'=rsquare[[5]],'mpg~disp[3]'=rsquare[[6]],'mpg~disp[4]'=rsquare[[7]],'mpg~disp[5]'=rsquare[[8]],'mpg~disp[6]'=rsquare[[9]],'mpg~disp[7]'=rsquare[[10]],'mpg~disp[8]'=rsquare[[11]],'mpg~disp[9]'=rsquare[[12]],'mpg~disp[10]'=rsquare[[13]],'mpg~disp[11]'=rsquare[[14]])),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

## -----------------------------------------------------------------------------
#(a)
set.seed(1234)
df<-data.frame(replicate(6,sample(c(1:10),10,rep=T)))
sd<-function(x)#compute the standard deviation
{
  y<-sqrt(mean((x-mean(x))^2))
  return(y)
}
round(vapply(df,sd,FUN.VALUE=0),3)
#(b)
name=c("zhang", "Li", "Wang", "Zhao", "Ding")
sex=c("F", "F", "M", "M", "M")
age=c(16, 17, 18, 16, 19)
height=c(167.5, 156.3, 177.3, 167.5, 170.0)
weight=c(55.0, 60.0, 63.0, 53.0, 69.5)
mydataframe<-data.frame(name,sex,age,height,weight)
logit<-vapply(mydataframe,is.numeric,FUN.VALUE=c(FALSE))
index<-seq(along=logit)[logit==TRUE]
mydataframe<-mydataframe[,index]
round(vapply(mydataframe,sd,FUN.VALUE=0),3)

## -----------------------------------------------------------------------------
#pure r code
set.seed(0)
pure_R<-function(N,discard)#N represents length of the chain,discard represents the length of the burn.
{
a<-2
b<-3
n<-100
x0<-30
y0<-0.3
sp<-matrix(0,N,2)
sp[1,]<-c(x0,y0)
for( i in 2:N)
{
  #x<-sp[i-1,2]
  sp[i,1]<-rbinom(1,n,sp[i-1,2])
  sp[i,2]<-rbeta(1,(sp[i,1]+a),(n-x0+b))
  i<-i+1
}
return(sp)
}
N<-15000
discard<-1500
chain_R<-pure_R(N,discard)[((discard+1):N),]

#Rcpp function
library(Rcpp)
library(microbenchmark)
func<-'NumericMatrix chain(int N)
{
int a=2;
int b=3;
double n=100;
double x0=30;
double y0=0.3;
NumericMatrix sp(N,2);
sp(0,0)=x0;
sp(0,1)=y0;
for(int j=1;j<N;j++)
{
  sp(j,0)=rbinom(1,n,sp(j-1,1))[0];
  sp(j,1)=rbeta(1,(sp(j,0)+a),(n-x0+b))[0];
}
return(sp);
}'
cppFunction( func )
chain_C<-chain(N)
chain_C<-chain_C[((discard+1):N),]
qqplot(chain_C[,1],chain_R[,1])
qqplot(chain_C[,2],chain_R[,2])
ts_1<-microbenchmark(meanR=mean(pure_R(15000)),meancpp=mean(chain(15000)))
ts_1<-summary(ts_1)[,c(1,3,5,6)]
x_html<-knitr::kable(ts_1,"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=10)

