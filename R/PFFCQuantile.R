#' Quantiles of Progressive First-failure Censored Order Statistics
#'
#' @param p Vector of quantiles
#' @param n Number of experimental groups
#' @param m Number of failures
#' @param k group size
#' @param R Progressive censoring scheme
#' @param order order of the progressive first-failure censored order statistics
#' @param QF Quantile function of the distribution
#' @param ... Additional arguments to be passed from the quantile of the distribution
#'
#' @return later
#' @importFrom stats uniroot na.omit
#' @export
#'
#' @examples
qpffc2 <- function(p = c(0.5, 0.75), n, m, k = 1, R, order, QF, ...)
{
  gam <- NA
  Cr <- NA
  for(i in 1 : order)
  {
    gam[i] <- m - i + 1 + sum(R[i : m])
  }
  for(i in 1 : order)
  {
    Cr[i] <- prod(gam[1 : i])
  }
  air <- array(dim = c(order, order))
  for(i in 1 : order)
  {
    for (j in 1 : order) {

      if(i != j)
      {
        air[i, j] <- 1/(gam[j] - gam[i])
      }
    }
  }
  A <- NA
  for(i in 1 : order)
  {
    A[i] = prod(na.omit(air[i,]))
  }
  progU_CDF <- function(u)
  {
    1 - (Cr[order] * sum((A/gam) * ((1 - u)^(k * gam))))
  }
  Obj  <- function(x, p)
  {
    return(progU_CDF(x) - p)
  }
  quan <- NA
  for(i in 1 : length(p))
  {
    quan[i] <- uniroot(Obj, interval = c(0, 1), p = p[i])$root
  }
  return(QF(quan, ...))
}
