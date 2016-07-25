#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// mat test(vec y){
// 	arma::mat w = zeros(3, 1);
// 	w(0) = log(mean(y));
// 	return w;
// }


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List l0araC(arma::mat x, arma::vec y, String family, double lam, int maxit, double eps) {
  // initialization
  int n = x.n_rows;
  int m = x.n_cols;
  arma::mat w = zeros(m, 1);
  arma::mat old_w(m, 1);
  arma::mat Xt(n, m);
  arma::mat xw(n, 1);
  arma::mat s1(n, 1);
  arma::mat A(n, n);
  arma::mat a = zeros(m, 1);
  arma::mat z(n, 1);

  int iter = 1;
  if(family=="gaussian" || family=="inv.gaussian" || family=="gamma"){
    w = inv(trans(x)*x+lam*eye(m,m))*(trans(x)*y);
  }
  if(family=="logit" || family=="inv.gaussian"){
    w(0) = log(mean(y));
  }

  Xt = x;
  while(iter < maxit) {
    old_w = w;
    xw = x*w;
    if(family=="gaussian"){
      s1 = xw;
      A = eye(n, n);
      Xt = repmat(trans(w % w), n, 1) % x;
    }
    if(family=="poisson"){
      s1 = exp(xw);
      A = diagmat(vec(s1));
    }
    if(family=="logit"){
      s1 = 1/(1+exp(-xw));
      a = s1 % (1-s1);
      A = diagmat(vec(a));
    }
    if(family=="gamma"){
      s1 = 1/xw;
      a = -s1 % s1;
      A = diagmat(vec(a));
    }
    if(family=="inv.gaussian"){
      s1 = 1/sqrt(xw);
      a = -pow(s1, 3)/2;
      A = diagmat(vec(a));
    }
    z = A*xw + y-s1;
    if(n > m) {
      w = solve(trans(Xt)*A*x+lam*eye(m,m), trans(Xt)*z);
    } else {
      w = trans(Xt)*solve(A*x*trans(Xt)+lam*eye(n,n), z);
    }
    if(family!="gaussian"){
      Xt = repmat(trans(w % w), n, 1) % x;
    }
    iter += 1;

    // check for convergence
    if(iter >= maxit){
      warning("Did not converge. Increase maxit.");
    }
    if(norm(w-old_w, 2) < eps) {
      break;
    }
  }
  for(int i=0; i<m; i++){
    if(std::abs(w(i, 0)) < 1e-3){
      w(i,0) = 0;
    }
  }
  return List::create(Named("beta") = w, Named("iter") = iter);
}





