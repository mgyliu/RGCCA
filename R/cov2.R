# cov2() is similar to cov() but has an additional argument.
# The denominator \eqn{n} (bias = TRUE) can be used (instead of \eqn{n-1})
# to give a biased estimator of the (co)variance.
# Note that if values are missing the "pairwise.complete.obs"
# option is used resulting in (co)-variance matrix that are not necessarily
# positive definite.
# @param x A numeric vector, matrix or data.frame.
# @param y A numeric vector, matrix or data.frame.
# @param bias A logical value. If bias = TRUE, \eqn{n} is used to give a biased estimator of the (co)variance.
# If bias = FALSE, \eqn{n-1} is used.
# @return \item{C}{Estimation of the variance (resp. covariance) of x (resp. x and y).}
# @title Variance and Covariance (Matrices)
# @importFrom stats cov

cov2 = function (x, y = NULL, bias = TRUE, alternateCov = NULL)
{
  n = NROW(x)
  options = c("mcd", "cellwise", "mve")

  stopifnot(alternateCov %in% options)

  # Annoying workaround because R sucks at handling NULL
  if (sum(alternateCov == "cellwise") > 0) {
    locScale.x <- cellWise::estLocScale(x)
    # the wrapped data is stored in $Xw. covariances get computed on this.
    Xw.x <- cellWise::wrap(x, locScale.x$loc, locScale.x$scale)$Xw    
    if (is.null(y)) {
      Xw.x.cov <- cov(Xw.x)
      return(ifelse(bias, Xw.x.cov, (n-1)/n * Xw.x.cov))
    } else {
      locScale.y <- cellWise::estLocScale(y)
      Xw.y <- cellWise::wrap(y, locScale.y$loc, locScale.y$scale)$Xw
      Xw.xy.cov <- cov(Xw.x, Xw.y)
      return(ifelse(bias, Xw.xy.cov, (n-1)/n * Xw.xy.cov))
    }
  }

  if (sum(alternateCov == "mcd") > 0) {    
    if (is.null(y)) {
      cov <- robustbase::covMcd(x)$cov  
    } else {
      cov <- robustbase::covMcd(cbind(x,y))$cov  
    }
    return(ifelse(bias, cov, (n-1)/n * cov))
  }

  if (sum(alternateCov == "mve") > 0) {
    if (is.null(y)) {
      cov <- rrcov::CovMve(x)$cov
    } else {
      cov <- rrcov::CovMve(cbind(x,y))$cov
    }
    # browser()
    return(ifelse(bias, cov, (n-1)/n * cov))
  }

  # Original

  if (is.null(y)) {
    x = as.matrix(x)
    if(bias){
      C = ((n - 1)/n) * stats::cov(x, use = "pairwise.complete.obs")
    }
    else{
      C = stats::cov(x, use = "pairwise.complete.obs")
    }
  }

  else{
    if(bias){
      C = ((n - 1)/n) * stats::cov(x, y, use = "pairwise.complete.obs")
    }
    else{
      C = stats::cov(x, y, use = "pairwise.complete.obs")
    }
  }
  # browser()
  return(C)
}
