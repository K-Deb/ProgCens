#' Simulating a generalized adaptive progressively first-failure type-II censored sample from a univariate distribution
#'
#' @param n Total number of groups placed on test at time zero
#' @param m Number of required failures
#' @param k Number of items in a group
#' @param Time1 First prefixed time point
#' @param Time2 Second prefixed time point
#' @param R Prefixed censoring scheme
#' @param non_zero If \code{TRUE}, only obtain samples with at least one failure
#' @param prob Probability of removal
#' @param CDF Cumulative Distribution function of the univariate distribution
#' @param QF Quantile function of the univariate distribution
#' @param ... Further arguments to be passed to \code{CDF} and \code{QF}
#' @param set_seed Random number generation seed
#'
#' @return LAter
#' @export
#'
#' @examples
rgapffc2 = function(n, m = NULL, k = 1, Time1 = Inf, Time2 = Inf, R = NULL, non_zero = T, prob = NULL,
                    CDF, QF, ..., set_seed = NULL)
{
  if (!is.null(set_seed))
    set.seed(set_seed)
  stopifnot(`CDF must be a quantile function` = is.function(CDF))
  stopifnot(`QF must be a quantile function` = is.function(QF))

  if(non_zero)
  {
    repeat
    {
      S = rapffc2(n = n, m = m, k = k, Time = Time1, R = R, prob = prob, CDF = CDF, QF = QF, ..., set_seed = set_seed)
      if(S[4, 1] < Time2)
        break
    }
  }else {
    S = rapffc2(n = n, m = m, k = k, Time = Time1, R = R, prob = prob, CDF = CDF, QF = QF, ..., set_seed = set_seed)
  }
  J = S[3, ]
  X = S[4, ]
  R = S[5, ]

  stopifnot(`No observations are recorded before the final completion time` = (X[1] < Time2))
  if (X[length(X)] < Time2) {
    J = J
    X = X
    R = R
  }
  else {
    X1 = c(X[- c(which(X > Time2))], Time2)
    R1 = c(R[1:(length(X1) - 1)], (n - (length(X1) - 1) - sum(R[1:(length(X1) - 1)])))
    J = J[1 : length(X1)]
    X = X1
    R = R1
  }
  return(rbind(J = J, X = X, R = R))

}
