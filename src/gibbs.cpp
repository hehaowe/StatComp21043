#include <Rcpp.h>
using namespace Rcpp;

//' @name gibbs
//' @title gibbs sample by Rcpp
//' @description Generate the gibbs sample by Rcpp,using the method of gibbs sampling.
//' @param N the sample size
//' @return a Numeric Matrix,sp stands for random sample of multivariate distribution.
//' @examples
//' \dontrun{
//' res <- gibbs_chain(N)
//' }
//' @export
//[[Rcpp::export]]
NumericMatrix gibbs_chain(int N)
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
}