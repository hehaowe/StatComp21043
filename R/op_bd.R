
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

#' @title f1
#' @description helper function(the first_term in the cross-validation score function.) 
#' @param h the bandwidth
#' @return The function value
#' @examples
#' \dontrun{
#' f1(h)
#' } 
#' @importFrom stats rnorm
#' @export
f1<-function(h)#the first_term in the cross-validation score function.
{
  X<-rnorm(32)
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
          sum0<-sum0+2*delta*Kf(delta*k,h,kernel)*Kf((delta*k-(X[i]-X[j])/h),h,kernel)
        }
      }
      else
      {
        for (k in 1:(1/delta))
        {
          sum0<-sum0+2*delta*Kf(delta*k,h,kernel)*Kf((delta*k-(X[i]-X[j])/h),h,kernel)
        }
      }
      temp<-temp+sum0
    }
  }
  y<-sum(temp)/(length(X)*h)
  return(y)
}

#' @title f2
#' @description helper function(the second term in the cross-validation score function.) 
#' @param h the bandwidth
#' @return The function value
#' @examples
#' \dontrun{
#' f2(h)
#' } 
#' @importFrom stats rnorm
#' @export

f2<- function(h)#the second term in the cross-validation score function.
{
  X<-rnorm(32)
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

#' @title fh_hat
#' @description $fh_hat(x)$ 
#' @param h the bandwidth
#' @param x the argument
#' @return The function value
#' @examples
#' \dontrun{
#' fh_hat(h,x)
#' } 
#' @importFrom stats rnorm
#' @export

fh_hat<- function(x,h)
{
  X<-rnorm(32)
  n<-length(X)
  temp=0
  for (i in 1:n) {
    temp=temp + kernel((X[i]-x)/h)
  }
  y<-temp/(n*h)
  return(y)
}


#' @title CV
#' @description helper function(f1-f2) 
#' @param h the bandwidth
#' @return The function value
#' @examples
#' \dontrun{
#' CV(h)
#' } 
#' @export
CV<-function(h)
{
  y<-f1(h)-f2(h)
  return(y)
}



#' @title optimal bandwidth
#' @description main function(get the optimal bandwidth and make a plot of f_hat(x) with optimal bandwidth.)
#' @param X several samples
#' @param kernel the type of kernel function(6 types of kernel functions in this function)
#' @return two plots (the Plot of CV_score function,the plot of f_hat(x) with optimal bandwidth)
#' @examples
#' \dontrun{
#' X<-mtcars$wt
#' op_bd(X,"Triangle")
#' } 
#' @importFrom graphics curve
#' @export

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
  curve(expr=CV(x),from=2,to = 10, yaxs="i", xaxs="i", type='l')
  h_opt<-optimize(CV,c(0,10),maximum=F,tol=0.001)$minimum
  h_opt
  dev.new()
  curve(expr=fh_hat(h_opt,x), from = 5,to = 40, yaxs="i", xaxs="i", type='l' )
}
