#' Simulating a progressively first-failure type-II censored sample from a univariate distribution
#'
#' @param n Total number of groups placed on test at time zero
#' @param m Number of required failures
#' @param k Number of items in a group
#' @param R Prefixed censoring scheme
#' @param prob Probability of removal
#' @param QF Quantile function of the univariate distribution
#' @param ... Further arguments to be passed to \code{QF}
#' @param set_seed Random number generation seed
#'
#' @return A dataframe with progressively type-II first-failure censored sample from a univariate distribution and the censoring scheme
#' @export
#'
#' @examples
rpffc2 = function(n, m = NULL, k = 1, R = NULL, prob = NULL, QF, ..., set_seed = NULL)
{
  stopifnot("QF must be a quantile function" = is.function(QF))
  U = rpc2unif(n = n, m = m, R = R, prob = prob, set_seed = set_seed)
  X = QF((1-(1-U[1,])^(1/k)), ...)
  return(rbind(X = X, R = U[2,]))
}
