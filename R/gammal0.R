#' fit a gamma regression model with l0 penalty
#' @description l0 iterative gamma regression
#' @usage gammal0(X, y, lam, maxit, eps)
#' @param X Input matrix, of dimension nobs x nvars; each row is an observation vector
#' @param y Response variable(positive quantitative variable)
#' @param lam A user supplied \code{lambda} value
#' @param maxit Maximum number of passes over the data for \code{lambda}
#' @param eps Convergence threshold
#' @return
#'  \item{beta}{a vector of coefficients}
#'  \item{iter}{number of iterations}
#' @author Wenchuan Guo <wguo007@ucr.edu>
gammal0 <- function(X, y, lam, maxit, eps){
  n <- nrow(X)
  m <- ncol(X)
  w <- solve(t(X)%*%X+lam*diag(m))%*%(t(X)%*%y)
  i <- 1
  Xt <- X
  while(i < maxit){
    old_w <- w
    xw <- X%*%w
    s1 <- xw
    a <- xw
    s1[xw > 10000] <- .0001
    s1[xw < -10000] <- -.0001
    a[xw > 10000] <- -.0001
    a[xw < -10000] <- .0001
    id <- xw < 10000 & xw > -10000
    s1[id] <- 1/xw[id]
    a[id] <- -s1[id]^2
    A <- diag(as.vector(a))
    z <- A%*%xw + (y-s1)
    if(n > m) {
      w <- solve(t(Xt)%*%A%*%X+lam*diag(m))%*%(t(Xt)%*%z)
    }else{
      w <- t(Xt)%*%solve(A%*%X%*%t(Xt) + lam*diag(n))%*%z
    }
    w2 <- as.vector(w^2)
    Xt <- t(t(X)*w2)
    i <- i+1
    if(i >= maxit) warning("Did not converge. Increase maxit.")
    if(norm(w-old_w,"F") < eps) break
  }
  w[abs(w)<10^(-3)] <- 0
  return(list(beta=w, iter=i, family="gamma"))
}
