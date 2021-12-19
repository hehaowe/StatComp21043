## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
op_bd<-function(X,kernel)
{
  n<-length(X)
  Kf<-function(x,h,kernel)#the kernel function
  {
    if(kernel=="Gaussian")
    {
      y<-dnorm(x/h)/h
    }
    
    if(kernel=="Triangle")
    {
      y<-(1-abs(x/h))*(x<h)*(x>-h)/h
    }
    
    if(kernel=="Epanechnikow")
    {
      y<-3/4*3/4*(1-(x/h)^2)*(x<h)*(x>-h)/h
    }
    if(kernel=="Quartic")
    {
      y<-15/16*(1-(x/h)^2)^2/h
    }
    if(kernel=="Triweight")
    {
      y<-35/32*(1-(x/h)^2)^3/h
    }
    return(y)
  }
  f1<-function(h)#the first_term in the cross-validation score function.
  {
    delta<-0.001
    temp<-0
    for(i in 1:length(X))
    {
      for(j in 1:length(X))
      {
        sum0<-0
        if(kernel=="Gaussian")
        {
          for (k in 1:(10/delta))
          {
            sum0<-sum0+2*delta*Kf(h*(delta*k),h,kernel)*Kf(h*(delta*k-(X[i]-X[j])/h),h,kernel)
          }
        }
        else
        {
          for (k in 1:(1/delta))
          {
            sum0<-sum0+2*delta*Kf(h*(delta*k),h,kernel)*Kf(h*(delta*k-(X[i]-X[j])/h),h,kernel)
          }
        }
        temp<-temp+sum0
      }
    }
    y<-sum(temp)/((length(X))^2*h)
    return(y)
  }
  f2<- function(h)#the second term in the cross-validation score function.
  {
    temp=0
    n<-length(X)
    for (i in 1:n)
    {
      for (j in 1:n)
      {
        if(i!=j)
        {
          temp=temp+Kf((X[i]-X[j]),h,kernel)
        }
      }
    }
    result=temp*2/(n*(n-2))
    return(result)
  }
  CV<-function(h)
  {
    y<-f1(h)-f2(h)
    return(y)
  }
  fh_hat<- function(h,x)
  {
    temp=0
    for (i in 1:n) {
      temp=temp + Kf((X[i]-x),h,kernel)
    }
    y<-temp/n
    return(y)
  }
  curve(expr=CV(x),from = 2, to = 10, yaxs="i", xaxs="i", type='l')
  h_opt<-optimize(CV,c(0,10),maximum=F,tol=0.001)$minimum
  h_opt
  dev.new()
  #curve(expr=fh_hat(h_opt,x), from = 5,to = 40, yaxs="i", xaxs="i", type='l' )
}
op_bd(mtcars$mpg,kernel="Triangle")


## -----------------------------------------------------------------------------
X=mtcars$wt
Y=mtcars$mpg
#local polynomial regression
lpe<-function(X,Y,p,kernel)
{
  if(kernel=="Triangle")
  {
    Kf=function(x,h){
      return((1-abs(x/h))*(x<h)*(x>-h)/h)
    }
    print("Use Triangle kernel function")
  }
  if(kernel=="Epanechnikow")
  {
    Kf=function(x,h){
      return(3/4*(1-(x/h)^2)*(x<h)*(x>-h)/h)
    }
    print("Use Epanechnikow kernel function")
  }
  if(kernel=="Gaussian")
  {
    Kf=function(x,h){
      return(dnorm(x/h)/h)
  }
    print("Use Gaussian kernel function")
  }
  if(kernel=="Quartic")
  {
    Kf=function(x,h){
      return(15/16*(1-(x/h)^2)^2/h)
  }
    print("Use Quartic kernel function")
  }
  if(kernel=="Triweight")
  {
    Kf=function(x,h){
      return(35/32*(1-(x/h)^2)^3/h)
  }
    print("Use Triweight kernel function")
  }
beta_hat<-function(x,h,X,Y,p)#calculate the value of beta hat.
{
  inx<-matrix(data=0,nrow=length(X),ncol=p)
  for(i in 1:length(Y)){
    for(j in 1:p)
    {
      inx[i,j] = (X[i]-x)^(j-1)
    }
  }
  Wx<-matrix(data=0,nrow=length(X),ncol=length(Y))
  for(i in 1:length(X))
  {
    for(j in 1:length(Y))
    {
      if(i==j)
      {
        Wx[i,j]=Kf((X[i]-x),h)
      }
      else
      {
        Wx[i,j] = 0
      }
    }
  }
  #print(k)
  result1<-solve(t(inx)%*%Wx%*%inx)
  result2<-t(inx)%*%Wx%*%Y
  beta=result1%*%result2
  return(beta)
}  
#fx(30,10,X,Y)
#plotting the figure of f_hat(x)
f<-function(h,X,Y){
  xrange<-range(X)#find the min and max point of the interval.
  d = ceiling(xrange[2])-floor(xrange[1])#find a interval containing these two values.
  x1<-c()#create a null vector.
  y1<-c()#同上
  #dividing this interval into 100 parts.
  for(i in 1:100){
    a1 = floor(xrange[1])+d/100*i 
    x1 = c(x1,a1)
    a2 = beta_hat(a1,h,X,Y,p)
    y1 = c(y1,a2[1])
  }
  plot(x1,y1,col="orange") #make a plot of x1,y1
  lines(x1,y1,type = "l",col= "red")
  title(xlab="x", ylab="fhatx")  
}

#create a new matrix to prepare for cross_validation.
m=length(Y)
n=length(X)
x0 = matrix(1:m*(m-1),m,m-1)
y0 = matrix(1:m*(m-1),m,m-1)
for(i in 1:m){
  t=1
  for(j in 1:m)
  {
    if(j!=i){
      x0[i,t]=X[j]
      y0[i,t]=Y[j]
      t=t+1
    }
  }
}
#use the method cross validation to get the optimal h
CV<- function(h)
{
  sum=0
  for(i in 1:m){
    l<-c()
    l<-beta_hat(X[i],h,x0[i,],y0[i,],p)#get beta_hat（-i)X
    temp=(Y[i]-l[1])^2
    sum=temp+sum
  }
  return(sum)
}
#use the function optimize()to get h_opt
h1<-optimize(CV,c(1,5),lower=0.1,upper=1,maximum=FALSE,tol=0.1)$minimum  
h2<-optimize(CV,c(1,5),lower=1,upper=2,maximum=FALSE,tol=0.1)$minimum 

#The following are plots
dev.new()
par(mfrow=c(2:1))
plot(X,Y)
lines(lowess(X,Y),col="blue",lwd=2,lty=2)
#f(h1,X,Y)
#title(main=paste("h_pot=", as.character(h1), sep = "") )
f(h2,X,Y)
title(main=paste("h_pot=", as.character(h2), sep = ""))
dev.new()
c<-seq(0.1,2,by=0.1)
opar <- par(no.readonly=TRUE) 
par(pin=c(0.6,0.6)) 
par(mfrow=c(4:5)) 
for(i in 1:20)
{
  f(c[i],X,Y)
  title(main=paste("when h=", as.character(c[i]), sep = ""))
}
par(opar)
}
lpe(X,Y,4,kernel="Quartic")

