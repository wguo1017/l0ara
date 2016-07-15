#' fit a generalized linear model with l0 penalty
#' @description An efficient algorithm for feature selection with L0 penalty.
#' @usage  l0ara(x, y, family, lam, maxit, eps)
#' @param x Input matrix, of dimension nobs x nvars; each row is an observation vector
#' @param y Response variable. Quantitative for \code{family="gaussian"}; positive quantitative for \code{family="gamma"} or \code{family="inv.gaussian"} ; a factor with two levels for \code{family="logit"}; non-negative counts for \code{family="poisson"}
#' @param family Response type(see above)
#' @param lam A user supplied \code{lambda} value
#' @param maxit Maximum number of passes over the data for \code{lambda}
#' @param eps Convergence threshold
#' @return
#'  \item{beta}{a vector of coefficients}
#'  \item{df}{number of nonzero coefficients}
#'  \item{iter}{number of iterations}
#'  \item{lambda}{The lambda used}
#' @author Wenchuan Guo <wguo007@ucr.edu>
#' @seealso \code{predict}, \code{coef} methods.
#' @examples
#' # Linear regression
#' # Generate design matrix and response variable
#' n <- 100
#' p <- 40
#' x <- matrix(rnorm(n*p), n, p)
#' beta <- c(1,0,2,3,rep(0,p-4))
#' noise <- rnorm(n)
#' y <- x%*%beta+noise
#' # fit sparse linear regression
#' res.gaussian <- l0ara(x, y, family="gaussian", 1)
#'
#' # predict for new observations
#' print(res.gaussian)
#' predict(res.gaussian, newx=matrix(rnorm(3,p),3,p))
#' coef(res.gaussian)
#'
#' # Logistic regression
#' # Generate design matrix and response variable
#' n <- 100
#' p <- 40
#' x <- matrix(rnorm(n*p), n, p)
#' beta <- c(1,0,2,3,rep(0,p-4))
#' prob <- exp(x%*%beta)/(1+exp(x%*%beta))
#' y <- rbinom(n, rep(1,n), prob)
#' # fit sparse logistic regression
#' res.logit <- l0ara(x, y, family="logit", 0.7)
#'
#' # predict for new observations
#' print(res.logit)
#' predict(res.logit, newx=matrix(rnorm(3,p),3,p))
#' coef(res.logit)
#'
#' # Poisson regression
#' # Generate design matrix and response variable
#' n <- 100
#' p <- 40
#' x <- matrix(rnorm(n*p), n, p)
#' beta <- c(1,0,0.5,0.3,rep(0,p-4))
#' mu <- exp(x%*%beta)
#' y <- rpois(n, mu)
#' # fit sparse logistic regression
#' res.pois <- l0ara(x, y, family="pois", 0.7)
#'
#' # predict for new observations
#' print(res.pois)
#' predict(res.pois, newx=matrix(rnorm(3,p),3,p))
#' coef(res.pois)

#' @export
l0ara <- function(x, y, family = c("gaussian", "logit", "gamma", "poisson", "inv.gaussian"), lam, maxit = 10^3, eps = 1e-04){
  # error checking
  family <- match.arg(family)
  if(missing(lam)){
    stop("'lam' was not found(must be a postiive number). ")
  }
  if(length(lam)>1){
    stop("Require lenght 1 for 'lam'.")
  }
  if(lam < 0) {
    warning("lambda < 0; set to 0")
    lam <- 0
  }
  np = dim(x)
  if (is.null(np) | (np[2] <= 1)){
    stop("x should be a matrix with 2 or more columns")
  }
  if (length(y)!=np[1]){
    stop("length of y not equal to the number of rows of x")
  }

#   # standardize
#   if (standardize) {
#     xx <- scale(x)
#     center <- colMeans(x)
#     scale <- apply(x,2,sd)
#     nd <- which(scale > 1e-6)
#     if (length(nd) != np[2]) xx <- xx[ ,nd, drop=FALSE]
#   } else {
#     xx <- x
#   }
#   if (family=="gaussian") {
#     yy <- y - mean(y)
#   } else {
#     yy <- y
#   }

  # fitting model
  if (family == "gaussian") {
    out <- gaussl0(x, y, lam, maxit, eps)
  }
  if (family == "gamma") {
    out <- gammal0(x, y, lam, maxit, eps)
  }
  if (family == "inv.gaussian") {
    out <- invl0(x, y, lam, maxit, eps)
  }
  if (family == "logit") {
    out <- logitl0(x, y, lam, maxit, eps)
  }
  if (family == "poisson") {
    out <- poisl0(x, y, lam, maxit, eps)
  }
  # output
  res <- list(beta = out$beta, df = sum(out$beta!=0), lambda = lam, iter = out$iter, family = out$family)
  class(res) <- "l0ara"
  return(res)
}
