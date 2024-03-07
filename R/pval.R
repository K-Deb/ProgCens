#' Monte Carlo simulation of p-value
#'
#' @param data A 2 by m matrix, containing failure times in first row and the progressive first-failure censoring scheme in second row
#' @param CDF Cumulative distribution function
#' @param Time Prefixed time point
#' @param k Number of items in a group
#' @param B number of Monte Carlo iteration
#' @param ... further arguments to be passed to CDF. Typically this is a vector of parameter estimates
#'
#' @return later
#' @export
#'
#' @examples
sim.pval = function(data, CDF, Time = Inf, k = 1, B = 10000, ...)
{
  R = data[2,]
  m = length(R)
  n = sum(R) + m
  stat_KS = AdapGoF(data = data, CDF = CDF, k = k, norm_approx = T, ...)$KS
  stat_CVM = AdapGoF(data = data, CDF = CDF, k = k, norm_approx = T, ...)$CVM
  stat_AD = AdapGoF(data = data, CDF = CDF, k = k, norm_approx = T, ...)$AD
  stat = c(stat_KS, stat_CVM, stat_AD)

  tmat = t(replicate(B, AdapGoF(data = rapffc2(n = n, m = m, k = k, Time = Time,
                                               R = R, CDF = pnorm, QF = qnorm, mean = 0,
                                               sd = 1)[4 : 5, ],
                                CDF = pnorm, k = k)))

  p_val = colMeans(tmat > matrix(rep(stat, B), ncol = 3, byrow = T))

  d = rbind(stat, p_val)
  colnames(d) = c("KS", "CVM", "AD")
  rownames(d)= c("Statistic", "p-value")
  return(d)

}
