#' fit a generalized linear model with l0 penalty
#' @description An adaptive ridge algorithm for feature selection with L0 penalty.
#' @usage  l0ara(x, y, family, lam, maxit, eps)
#' @param x Input matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y Response variable. Quantitative for \code{family="gaussian"}; positive quantitative for \code{family="gamma"} or \code{family="inv.gaussian"} ; a factor with two levels for \code{family="logit"}; non-negative counts for \code{family="poisson"}.
#' @param family Response type(see above).
#' @param lam A user supplied \code{lambda} value. Try use \code{cv.l0ara} first to select optimal tunning and then refit with \code{lam.min}.
#' @param maxit Maximum number of passes over the data for \code{lambda}. Default value is \code{1e3}.
#' @param eps Convergence threshold. Default value is \code{1e-4}.
#' @details The sequence of models indexed by the parameter lambda is fit using adptive ridge algorithm. The objective function for generalized linear models (including \code{family} above) is defined to be \deqn{-(log likelihood)+(\lambda/2)*|\beta|_0} This adaptive ridge algorithm is developed to approximate L0 penalized generalized linear models with sequential optimization and is efficient for high-dimensional data.
#' @return An object with S3 class "l0ara" containing:
#'  \item{beta}{A vector of coefficients}
#'  \item{df}{Number of nonzero coefficients}
#'  \item{iter}{Number of iterations}
#'  \item{lambda}{The lambda used}
#'  \item{x}{Design matrix}
#'  \item{y}{Response variable}
#' @author
#' Wenchuan Guo <wguo007@ucr.edu>, Zhenqiu Liu <Zhenqiu.Liu@cshs.org>
#' @seealso \code{cv.l0ara}, \code{predict}, \code{coef} methods.
#' @useDynLib l0ara
#' @importFrom Rcpp evalCpp
#' @import stats
#' @import graphics
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
#' # fit sparse Poisson regression
#' res.pois <- l0ara(x, y, family="poisson", 0.7)
#'
#' # predict for new observations
#' print(res.pois)
#' predict(res.pois, newx=matrix(rnorm(3,p),3,p))
#' coef(res.pois)

#' @export
l0ara <- function(x, y, family = c("gaussian", "logit", "gamma", "poisson", "inv.gaussian"), lam, maxit = 10^3, eps = 1e-04){
  # error checking
  if (class(x) != "matrix") {
    tmp <- try(x <- model.matrix(~0 + ., data = x), silent = TRUE)
    if (class(tmp)[1] == "try-error")
      stop("x must be a matrix or able to be coerced to a matrix")
  }
  if (class(y) != "numeric") {
    tmp <- try(y <- as.numeric(y), silent = TRUE)
    if (class(tmp)[1] == "try-error")
      stop("y must numeric or able to be coerced to numeric")
  }
  y <- drop(y)

  family <- match.arg(family)
  if(missing(lam)){
    stop("'lam' was not found(must be a postiive number)")
  }
  if(length(lam)>1){
    stop("Require lenght 1 for 'lam'")
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

#   # fitting model
#   if (family == "gaussian") {
#     out <- gaussl0(x, y, lam, maxit, eps)
#   }
#   if (family == "gamma") {
#     out <- gammal0(x, y, lam, maxit, eps)
#   }
#   if (family == "inv.gaussian") {
#     out <- invl0(x, y, lam, maxit, eps)
#   }
#   if (family == "logit") {
#     out <- logitl0(x, y, lam, maxit, eps)
#   }
#   if (family == "poisson") {
#     out <- poisl0(x, y, lam, maxit, eps)
#   }
  out <- l0araC(x, y, family, lam, maxit, eps)
  # output
  res <- list(beta = drop(out$beta), df = sum(out$beta!=0), lambda = lam, iter = out$iter, family = family, x = x, y = y)
  class(res) <- "l0ara"
  return(res)
}
