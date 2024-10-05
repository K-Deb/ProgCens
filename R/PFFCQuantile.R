#' Quantile Function of Progressive First-failure Censored Order Statistics
#'
#' @param p probabilities
#' @param n number of experimental groups
#' @param m desired number of failures
#' @param k group sizes
#' @param R progressive censoring scheme
#' @param order order of the progressive first-failure censored order statistic
#' @param QF quantile function of the distribution
#' @param ... additional arguments to be passed to \code{QF}
#' @param precBits same as \code{Rmpfr}
#'
#' @return Later
#' @importFrom stats uniroot
#' @importFrom Rmpfr mpfr outer colSums apply
#' @export
#'
#' @examples
qpffc2 <- function(p = seq(0, 1, 0.25), n, m, k = 1, R, order, QF, ..., precBits = 1024)
{
  gam <- mpfr(rev(cumsum(R[m : 1]))[1 : order] + (m : (m - order + 1)), precBits = precBits)
  Cr <- cumprod(gam)
  air <- t(1/outer(gam, gam, "-"))
  diag(air) <- 1
  A <- apply(air, 1, prod)
  diag(air) <- NA
  progU_CDF <- function(u)
  {
    as.numeric(1 - (Cr[order] * colSums((A/gam) * outer(k * gam, 1 - u, \(x, y) y ^ x))))
  }
  Obj  <- function(x, p)
  {
    return(progU_CDF(x) - p)
  }
  quan <- NA
  for(i in 1 : length(p))
  {
    quan[i] <- uniroot(Obj, interval = c(0, 1), p = p[i], extendInt = "yes", maxiter = 1e5)$root
  }
  return(QF(quan, ...))
}
