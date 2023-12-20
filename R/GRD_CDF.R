

#' The Generalized Rayleigh Distribution
#'
#' @param q Vector of quantiles
#' @param shape shape parameter
#' @param scale scale parameter
#' @param lower.tail logical
#' @param log.p logical
#'
#' @return
#' @export
#'
#' @examples
pgr = function(q, shape, scale = 1, lower.tail=TRUE, log.p=FALSE)
{
  cdf = (1-exp(-scale^2*q^2))^shape
  if(isTRUE(lower.tail)==T & isTRUE(log.p)==F)
  {
    return(cdf)
  }else if (isTRUE(lower.tail)==F & isTRUE(log.p)==F)
  {
    return(1-cdf)
  }else if (isTRUE(lower.tail)==F & isTRUE(log.p)==T)
  {
    return(log(1-cdf))
  }else { return(log(cdf))}
}

