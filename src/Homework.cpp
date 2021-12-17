#include <Rcpp.h>
using namespace Rcpp;

//' @title An Rcpp function for Exercise 9.8 (page 278, Statistical Computing with R).
//' @description An Rcpp function for Exercise 9.8 (page 278, Statistical Computing with R).
//' @param a A parameter the bivariate density.
//' @param b A parameter the bivariate density.
//' @param n A parameter the bivariate density.
//' @param M Length of the chain.
//' @return A chain.
//' @export
// [[Rcpp::export]]
NumericMatrix CppFun(double a, double b, int n, int M) {
  NumericMatrix A(M,2); 
  A(0, 0) = A(0, 1) = 0;
  double x; double y;
  for(int j = 1; j < M; j++)
  {
    y = A(j-1, 1);
    A(j, 0) = rbinom(1, n, y)[0];
    x = A(j, 0);
    A(j, 1) = rbeta(1, x+a, n-x+b)[0];
  }
  return A;
}