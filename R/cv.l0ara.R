#' Cross validation for l0ara
#' @description Does k-fold cross-validation for l0ara, produces a plot, and returns the optimal \code{lambda}
#' @usage  cv.l0ara(x, y, family, lam, measure, nfolds, maxit, eps, seed)
#' @param x Input matrix as in \code{l0ara}.
#' @param y Response variable as in \code{l0ara}.
#' @param family Response type as in \code{l0ara}.
#' @param lam A user supplied \code{lambda} sequence in descending or asecending order.
#' @param measure Loss function used for corss validation. \code{measurer="mse"} or \code{"mae"} for all models. \code{"measure"="class"} only for logsitic regression.
#' @param nfolds Number of folds. Default value is 10. Smallest value is 3.
#' @param maxit Maximum number of passes over the data for \code{lambda}. Default value is \code{1e3}.
#' @param eps Convergence threshold. Default value is \code{1e-4}.
#' @param seed Seed of random number generator.
#' @details This function calls \code{l0ara} \code{nfolds} times, each time leaving out \code{1/nfolds} of the data. The cross-validation error is based on etiher mean square error (\code{mse}) or mean absolute error (\code{mae}).
#' @return An object with S3 class "cv.l0ara" containing:
#'  \item{cv.error}{The mean cross validated error for given lambda sequence}
#'  \item{lam.min}{The lambda gives min cv.error}
#'  \item{lambda}{The lambda used}
#'  \item{measure}{Type of measure}
#'  \item{family}{Model used}
#'  \item{x}{Design matrix}
#'  \item{y}{Response variable}
#' @author
#' Wenchuan Guo <wguo007@ucr.edu>, Zhenqiu Liu <Zhenqiu.Liu@cshs.org>
#' @seealso \code{l0ara}.
#' @examples
#' #' # Linear regression
#' # Generate design matrix and response variable
#' n <- 100
#' p <- 40
#' x <- matrix(rnorm(n*p), n, p)
#' beta <- c(1,0,2,3,rep(0,p-4))
#' noise <- rnorm(n)
#' y <- x%*%beta+noise
#' lam <- c(0.1, 0.3, 0.5)
#' fit <- cv.l0ara(x, y, family="gaussian", lam, measure = "mse")
#' @export

cv.l0ara <- function(x, y, family = c("gaussian", "logit"), lam, measure = c("mse", "mae","class"), nfolds=10,  maxit = 10^3, eps = 1e-04, seed){
  measure <- match.arg(measure)
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
  if (!is.null(lam) && length(lam) < 2) {
    stop("Need a sequence values for lambda")
  }
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3")
  }
  if (!missing(seed)) {
    set.seed(seed)
  }

  #initialization
  n <- nrow(x)
  p <- ncol(x)
  res <- as.list(1:nfolds)
  error <- matrix(NA, length(lam), nfolds)
  id <- sample(rep(1:nfolds, length = n))
  for (i in 1:nfolds){
    which <- id == i
    yy <- y[!which]
    xx <- x[!which, ,drop=FALSE]
    res[[i]] <- lapply(lam, function(ll) l0ara(xx,yy,lam=ll,family=family,maxit=maxit,eps=eps))
    pred <- lapply(res[[i]], predict,x[which,],type="response")
    if(measure == "mse") {
      error[,i] <- do.call(cbind, lapply(pred, function(pp) mean((pp-y[which])^2)))
    }
    if(measure == "mae") {
      error[,i] <- do.call(lapply(pred, function(pp) mean(abs(pp-y[which]))))
    }
  }
  # cv.std <- apply(error, 1, sd)/sqrt(n)
  cv.error <- rowMeans(error)
  lam.min <- lam[which.min(cv.error)]
  out <- list(cv.error=cv.error, lam.min=lam.min, measure=measure, lam=lam, family=family, x=x, y=y)
  class(out) <- "cv.l0ara"
  return(out)
}
