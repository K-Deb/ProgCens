

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
#' @name GR
NULL
#> NULL
#' @return
#' @export
#'
#' @examples

#' @rdname GR
#' @export
dGR = function(x, shape, rate = 1, scale = 1/rate, log = FALSE)
{
  if (!missing(rate) && !missing(scale)) {
    if (abs(rate * scale - 1) < 1e-15)
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  rate = 1/scale
  dens = 2 * shape * rate ^ 2 * x * exp(- (rate * x) ^ 2) * (- expm1(- (rate * x) ^ 2)) ^ (shape - 1)
  if(log)
  {
    return(log(dens))
  }else { return(dens)}
}
#' @rdname GR
#' @export
pGR = function(q, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
{
  if (!missing(rate) && !missing(scale)) {
    if (abs(rate * scale - 1) < 1e-15)
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  rate = 1/scale
  cdf = (- expm1(- (rate * q) ^ 2)) ^ shape
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
#' @rdname GR
#' @export
qGR = function(p, shape, rate = 1, scale = 1/rate, lower.tail=TRUE)
{
  if (!missing(rate) && !missing(scale)) {
    if (abs(rate * scale - 1) < 1e-15)
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  rate = 1/scale
  if(!lower.tail) {p = 1 - p}
  qf = sqrt((1/rate ^2) * log((1 - p ^ (1/shape)) ^ (-1)))
  return(qf)
}
#' @rdname GR
#' @export
rGR = function(n, shape, rate = 1, scale = 1/rate)
{
  if (!missing(rate) && !missing(scale)) {
    if (abs(rate * scale - 1) < 1e-15)
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  rate = 1/scale
  out = qGR(runif(n), shape = shape, rate = rate)
  return(out)
}
