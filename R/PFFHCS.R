#' Simulating a progressively first-failure type-II hybrid censored sample from a univariate distribution
#'
#' @param n Total number of groups placed on test at time zero
#' @param m Number of required failures
#' @param k Number of items in a group
#' @param Time Prefixed time point
#' @param R Prefixed censoring scheme
#' @param non_zero Logical. {TRUE} if non-zero number of failures are required (default).
#' @param prob Probability of removal
#' @param QF Quantile function of the univariate distribution
#' @param ... Further arguments to be passed to \code{CDF} and \code{QF}
#' @param set_seed Random number generation seed
#'
#' @return later
#' @export
#'
#' @examples
rpffhc2 = function(n, m = NULL, k = 1, Time, R = NULL, non_zero = T, prob = NULL, QF, ..., set_seed = NULL)
{
  if(!is.null(set_seed)) set.seed(set_seed)
  stopifnot("QF must be a quantile function" = is.function(QF))
  if(non_zero)
  {
    repeat
    {
      S = rpffc2(n = n, m = m, k = k, R = R, prob = prob, QF, ..., set_seed = set_seed)
      if(S[1,1] < Time)
        break
    }
  }else
  {
    S = rpffc2(n = n, m = m, k = k, R = R, prob = prob, QF, ..., set_seed = set_seed)
  }
  X = S[1,]
  R = S[2,]

  stopifnot("No observations are recorded before completion time" = (X[1] < Time))

  if(X[length(X)] < Time)
  {
    X = X
    R = R
  }else{
    X1 = c(X[-c(which(X > Time))], Time)
    R1 = c(R[1 : (length(X1) - 1)], (n - (length(X1) - 1) - sum(R[1:(length(X1)-1)])))
    X = X1
    R = R1}
  return(rbind(X = X, R = R))
}
