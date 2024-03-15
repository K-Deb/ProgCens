#' Simulating quantiles of EDF based statistics
#'
#' @param alp level of significance
#' @param data A 2 by m matrix, containing failure times in first row and the progressive first-failure censoring scheme in second row
#' @param k Number of items in a group
#' @param B number of Monte Carlo iteration
#' @importFrom stats quantile
#' @return later
#' @export
#'
#' @examples
sim.qEDF = function(alp = 0.05, data, k = 1, B = 1e4)
{
  stopifnot("data must be a 2 by m matrix" = (is.matrix(data) & (nrow(data) == 2)))
  R = data[2,]
  m = length(R)
  n = m + sum(R)
  outmat = t(replicate(B, unlist(AdapGoF(rapffc2(n = n, m = m, k = k, Time = Inf, R = R,
                                                 CDF = pnorm, QF = qnorm, mean = 0, sd = 1)
                                         [4:5,], pnorm, k = k))))
  out = data.frame(KS = quantile(outmat[,1], probs = 1 - alp),
                   CVM = quantile(outmat[,2], probs = 1 - alp),
                   AD = quantile(outmat[,3], probs = 1 - alp))
  rownames(out) = (paste("alpha = ", alp))
  return(out)
}
