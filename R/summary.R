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
predict.l0ara <- function(object, newx, type=c("link", "response", "coefficients", "class"), ...){
  type <- match.arg(type)
  beta <- coef.l0ara(object)
  if(missing(newx)){
    newx <- object$x
  }
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

#' print coefficients from a "cv.l0ara" object.
#' @description This function print the coefficients from the model with the optimal \code{lambda}.
#' @param object Fitted "cv.l0ara" object.
#' @param ... Not used argument.
#' @details This function fit the model with the optimal \code{lambda} first and then print the coefficients. This function makes it easier to use the results to make a prediction or to see the fitted model.
#' @return The object returns the coefficients.
#' @author Wenchuan Guo <wguo007@ucr.edu>
#' @seealso \code{predict} method and \code{l0ara} function.
#' @export
coef.cv.l0ara <- function(object, ...){
  fit <- l0ara(object$x, object$y, object$family, object$lam.min)
  coefs <- as.vector(fit$beta)
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
  family <- switch(x$family, gaussian = "Linear regression", logit = "Logistic regression", poisson = "Poisson regression", inv.gaussian = "Inverse gaussian regression", gamma = "Gamma regression" )
  cat("Lambda used : ", x$lam, "\n")
  cat("Model : ", family, "\n")
  cat("Iterations : ", x$iter, "\n")
  cat("Degree of freedom : ", x$df,"\n")
}

#' summarizing the fits from a "cv.l0ara" object.
#' @description This function print the general information of the cross validated fit.
#' @param x Fitted "cv.l0ara" object.
#' @param ... Not used argument.
#' @details This function makes it easier to see the cross validation results.
#' @author Wenchuan Guo <wguo007@ucr.edu>
#' @seealso \code{predict}, \code{coef} methods and \code{l0ara} function.
#' @export
print.cv.l0ara <- function(x, ...){
  measure <- switch(x$measure, mse = "Mean square error (MSE)", mae = "Mean absolute error (MAE)")
  family <- switch(x$family, gaussian = "Linear regression", logit = "Logistic regression", poisson = "Poisson regression", inv.gaussian = "Inverse gaussian regression", gamma = "Gamma regression" )
  cat("Number of Lambda used : ", length(x$lam), "\n")
  cat("Optimal Lambda : ", x$lam.min, "\n")
  cat("Model : ", family, "\n")
  cat("Measure : ", measure, "\n")
  cat("Minimumn error : ",min(x$cv.error), "\n")
}

#' plot for an "l0ara" object
#' @description Two plots are availiable: a plot of fitted value against linear predictor; \code{roc}(\code{auc}) curve for \code{family="logit"}.
#' @param x Fitted "l0ara" object.
#' @param auc logical; if \code{TRUE}, produces \code{auc} curve for \code{family=logit}.
#' @param split logical; if if \code{TRUE}, produces seperate plots.
#' @param ... Not used argument.
#' @author Wenchuan Guo <wguo007@ucr.edu>
#' @seealso \code{predict}, \code{coef} methods and \code{l0ara} function.
#' @export
plot.l0ara <- function(x, auc = FALSE, split = FALSE, ...){
  nplots <- ifelse(auc, 2, 1)
  if(!split){
    par(mfrow=c(1,nplots))
  }
  resp <- predict(x, type="response")
  lp <-predict(x, type="link")
  plot(lp, resp, xlab="Linear predictor", ylab="Fitted value", pch=20, main="Linear predictor v.s. Fitted", ...)
  points(lp,x$y,col=3, pch=20)
  legend("bottomright", legend = c("Fitted","Truth"), col=c(1,3), pch=rep(20,2))

  if(x$family=="logit" & auc){
    pred <- predict(x, x$x, type="class")
    fp <- mean(pred[x$y==0]==1)
    tp <- mean(pred[x$y==1]==1)
    area <- fp*tp/2+(tp+1)*(1-fp)/2
    plot(NA, xlab="False positive rate", ylab="True positive rate", xlim=c(0,1), ylim=c(0,1), main=paste("ROC curve for lambda =", x$lam,";", "AUC =", area), ...)
    lines(c(0,fp,1), c(0,tp,1), col=4)
    lines(c(0,1), c(0,1), lty=2)
  }
}
