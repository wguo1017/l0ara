#' fit a linear regression model with l0 penalty
#' @description l0 iterative regression
#' @usage  gaussl0(X, y, lam, maxit, eps)
#' @param X Input matrix, of dimension nobs x nvars; each row is an observation vector
#' @param y Response variable(quantitative variable)
#' @param lam A user supplied \code{lambda} value
#' @param maxit Maximum number of passes over the data for \code{lambda}
#' @param eps Convergence threshold
#' @return
#'  \item{beta}{a vector of coefficients}
#'  \item{iter}{number of iterations}
#' @author Wenchuan Guo <wguo007@ucr.edu>
gaussl0 <- function(X, y, lam, maxit, eps) {
	n <- nrow(X)
	m <- ncol(X)
	w <- solve(t(X)%*%X+lam*diag(m))%*%(t(X)%*%y)
	i <- 0
	while(i < maxit){
		i <- i+1
		u <- w
		up <- abs(u)^2
		Xu <- t(t(X)*as.vector(up))
		XXfull <- t(Xu)%*%X
		Xyfull <- t(Xu)%*%y
		w <- solve(XXfull+lam*diag(m))%*%Xyfull
		if(norm(w-u,"F")<eps) break
		if(i >= maxit) warning("Did not converge. Increase maxit.")
	}
	w[abs(w)<10^(-3)] <- 0
	return(list(beta=w, iter=i, family="gaussian"))
}




