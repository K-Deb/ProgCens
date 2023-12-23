

#' The Generalized Rayleigh Distribution
#'
#' @param x,q Vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param shape shape parameter.
#' @param scale scale parameter.
#' @param rate  an alternative way to specify the scale.
#' @param lower.tail logical; if \code{TRUE}
#' @param log logical
#' @param log.p logical
#'
#' @return
#' @export
#'
#' @examples


dGR = function(x, shape, scale = 1, rate = 1/scale, log = FALSE)
{
  scale = 1/rate
  dens = 2*shape*scale^2*x*exp(-scale^2*x^2)*(1-exp(-scale^2*x^2))^(shape-1)
  if(log)
  {
    return(log(dens))
  }else { return(dens)}
}
#' @rdname dGR
pGR = function(q, shape, scale = 1, rate = 1/scale, lower.tail=TRUE, log.p=FALSE)
{
  scale = 1/rate
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
#' @rdname dGR
qGR = function(p, shape, scale = 1, rate = 1/scale, lower.tail=TRUE)
{
  scale = 1/rate
  if(!lower.tail) {p = 1 - p}
  qf = sqrt((1/scale ^2) * log((1 - p ^ (1/shape)) ^ (-1)))
  return(qf)
}
#' @rdname dGR
rGR = function(n, shape, scale = 1, rate = 1/scale)
{
  scale = 1/rate
  out = qGR(runif(n), shape = shape, scale = scale, rate = rate)
  return(out)
}
