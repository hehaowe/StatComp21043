#' @title Kf
#' @description helper function(choose the kernel function) 
#' @param x the argument
#' @param h the bandwidth
#' @param kernel the kernel function
#' @return The function value
#' @examples
#' \dontrun{
#' Kf(x,h,"Triangle")
#' } 
#' @export
Kf<-function(x,h,kernel)
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

#' @title CV
#' @description helper function(calculate the optimal bandwidth) 
#' @param h the bandwidth
#' @examples
#' \dontrun{
#' CV(h)
#' } 
#' @export
CV<- function(h)
{
  X<-rnorm(0,1,32)
  Y<-rnorm(0,2,32)
  m=length(X)
  x0 = matrix(1:m*(m-1),m,m-1)
  y0 = matrix(1:m*(m-1),m,m-1)
  sum=0
  for(i in 1:m){
    l<-c()
    l<-beta_hat(X[i],h,x0[i,],y0[i,],3)#get beta_hat（-i)X
    temp=(Y[i]-l[1])^2
    sum=temp+sum
  }
  return(sum)
}


#' @title beta_hat
#' @description helper function(calculate the beta_hat of local polynomial estimation.)
#' @param x the argument
#' @param h the bandwidth
#' @param X the X of multivariate sample point (X,Y)
#' @param Y the X of multivariate sample point (X,Y)
#' @param p the order of polynomial 
#' @examples
#' \dontrun{
#' beta_hat(x,h,X,Y,3)
#' } 
#' @export
beta_hat<-function(x,h,X,Y,p)#calculate the value of beta hat.
{
  inx<-matrix(data=0,nrow=length(X),ncol=p)
  for(i in 1:length(Y))
  {
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
        Wx[i,j]=Kf((X[i]-x),h,kernel)
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


#' @title f
#' @description helper function(plotting the figure of f_hat(x))
#' @param h the bandwidth
#' @param X the X of multivariate sample point (X,Y)
#' @param Y the X of multivariate sample point (X,Y)
#' @examples
#' \dontrun{
#' f(h,X,Y)
#' } 
#' @export
f<-function(h,X,Y){
  xrange<-range(X)#find the min and max point of the interval.
  d = ceiling(xrange[2])-floor(xrange[1])#find a interval containing these two values.
  x1<-c()#create a null vector.
  y1<-c()#同上
  #dividing this interval into 100 parts.
  for(i in 1:100){
    a1 = floor(xrange[1])+d/100*i 
    x1 = c(x1,a1)
    a2 = beta_hat(a1,h,X,Y,3)
    y1 = c(y1,a2[1])
  }
  plot(x1,y1,col="orange") #make a plot of x1,y1
  lines(x1,y1,type = "l",col= "red")
  title(xlab="x", ylab="beta_hatx")  
}


#' @title local polynomial estimation
#' @description main function(Use R to do local polynomial estimation and make a plot of the beta_hat.)
#' @param X the X of multivariate sample point (X,Y)
#' @param Y the X of multivariate sample point (X,Y)
#' @param p the order of polynomial 
#' @param kernel the type of kernel function(6 types of kernel functions in this function)
#' @return three plots of the polynomial estimation(the Plot of X and estimated function,the plot of beta_hat with optimal bandwidth,plots of beta_hat with the bandwidths of the interval (0,1))
#' @examples
#' \dontrun{
#' lpe(X,Y,P,"Triangle")
#' } 
#' @importFrom grDevices dev.new
#' @importFrom graphics lines par title
#' @importFrom stats dnorm rnorm  kernel  lowess optimize
#' @export
lpe<-function(X,Y,p,kernel)
{
  Kf<-function(x,h,kernel)
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
  beta_hat<-function(x,h,X,Y,p)#calculate the value of beta hat.
  {
    inx<-matrix(data=0,nrow=length(X),ncol=p)
    for(i in 1:length(Y))
    {
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
          Wx[i,j]=Kf((X[i]-x),h,kernel)
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
    title(xlab="x", ylab="beta_hatx")  
  }
  
  #create a new matrix to prepare for cross_validation.
  m=length(X)
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
  return(par(opar))
}

