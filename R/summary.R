#' make predictions from a "l0ara" object.
#' @description This function makes predictions from the model.
#' @param object Fitted "l0ara" object.
#' @param newx Matrix of new values for x at which predictions are to be made. Must be a matrix.
#' @param type Type of prediction required. "link" gives the linear predictors(for "gaussian" models it gives the fitted values). "response" gives the fitted probabilities for "logit" and fitted mean for "poisson". "coefficients" gives the coefficients which is same as "coef" function. "class" (applies only to "logit") produces the class label corresponding to the maximum probability.
#' @param ... Not used argument.
#' @details This function makes it easier to use the results to make a prediction or to see the fitted model.
#' @return The object returned depends the functions.
#' @author Wenchuan Guo <wguo007@ucr.edu>
#' @seealso \code{coef} method and \code{l0ara} function.
#' @export
predict.l0ara <- function(object, newx, type=c("link", "response", "coefficients"), ...){
  type <- match.arg(type)
  beta <- coef.l0ara(object)
  if (type=="coefficients") return(beta)
  eta <- newx%*%beta
  if (type=="link" || object$family=="gaussian") return(drop(eta))
  response <- switch(object$family, logit = exp(eta)/(1+exp(eta)), poisson = exp(eta), gamma = 1/eta, inv.gaussian = 1/sqrt(eta))
  if (type=="response") return(drop(response))
  if (type=="class") {
    if (object$family=="logit") {
      return(drop(1*(eta>0)))
    } else {
      stop("type='class' can only be used with family='binomial'")
    }
  }
}

#' print coefficients from a "l0ara" object.
#' @description This function print the coefficients from the model.
#' @param object Fitted "l0ara" object.
#' @param ... Not used argument.
#' @details This function makes it easier to use the results to make a prediction or to see the fitted model.
#' @return The object returns the coefficients.
#' @author Wenchuan Guo <wguo007@ucr.edu>
#' @seealso \code{predict} method and \code{l0ara} function.
#' @export
coef.l0ara <- function(object, ...){
  coefs <- as.vector(object$beta)
  p <- length(coefs)
  names(coefs)[1] <- "Intercept"
  names(coefs)[2:p] <- paste0("X",1:(length(coefs)-1))
  return(coefs)
}

#' summarizing the fits from a "l0ara" object.
#' @description This function print the general information of the fit.
#' @param x Fitted "l0ara" object.
#' @param ... Not used argument.
#' @details This function makes it easier to see the fitted model.
#' @author Wenchuan Guo <wguo007@ucr.edu>
#' @seealso \code{predict}, \code{coef} methods and \code{l0ara} function.
#' @export
print.l0ara <- function(x, ...){
  cat("Lambda used = ", x$lam, "\n")
  cat("Model =", x$family, "\n")
  cat("Iterations =", x$iter, "\n")
  cat("Degree of freedom =", x$df,"\n")
}
