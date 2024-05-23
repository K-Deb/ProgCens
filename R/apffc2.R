#' Simulating an adaptive progressively first-failure type-II censored sample from a univariate distribution
#'
#' @param n Total number of groups placed on test at time zero
#' @param m Number of required failures
#' @param k Number of items in a group
#' @param Time Prefixed time point
#' @param R Prefixed censoring scheme
#' @param prob Probability of removal
#' @param CDF Cumulative Distribution function of the univariate distribution
#' @param QF Quantile function of the univariate distribution
#' @param ... Further arguments to be passed to \code{CDF} and \code{QF}
#' @param set_seed Random number generation seed
#'
#' @return Laetr
#' @export
#'
#' @examples
rapffc2 = function(n, m = NULL, k = 1, Time = Inf, R = NULL, prob = NULL, CDF, QF, ..., set_seed = NULL)
{
  if(!is.null(set_seed)) set.seed(set_seed)
  stopifnot("CDF must be a quantile function" = is.function(CDF))
  stopifnot("QF must be a quantile function" = is.function(QF))
  S = rpffc2(n = n, m = m, k = k, R = R, prob = prob, QF = QF, ...)
  X1 = S[1,]
  R1 = S[2,]
  J = ifelse(!any(X1<Time), 0, max(which(X1<Time|X1 == Time)))
  if(J == 0)
  {
    R2 = c(rep(0, (m - 1)), n - m)
    U = rpc2unif(n = n - 1, m = m - 1, R = R2[(J+2) : m], prob = NULL)
    ORD = 1 - (1 - U[1,])^(1/k)
    X2 = X1
    X2[(J+2):m] = QF((CDF(X2[(J+1)], ...) +
                        ORD * (1-CDF(X2[(J+1)], ...))), ...)
  }else
  {
    if(J + 1 == m)
    {
      X2 = X1
      R2 = c(R1[1:J], n-m-sum(R1[1:J]))
    }else
    {
      if(J == m)
      {
        X2 = X1
        R2 = R1
      }else
      {
        if(J + 1 < m)
        {
          R2 = c(R1[1:J], rep(0,(m-J-1)), n-m-sum(R1[1:J]))
          U = rpc2unif(n = (n-sum(R2[1:J])-J-1), m = (m-J-1), R = R2[(J+2) : m], prob = NULL)
          ORD = 1 - (1 - U[1,])^(1/k)
          X2 = X1
          X2[(J+2):m] = QF((CDF(X2[(J+1)], ...) +
                               ORD * (1-CDF(X2[(J+1)], ...))), ...)
        }
      }
    }
  }
  return(rbind(X1 = X1, R1 = R1, J = rep(J, m), X2 = X2, R2 = R2))
}
