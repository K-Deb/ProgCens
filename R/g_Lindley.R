#' The Generalized Lindley Distribution
#'
#' @param x,q Vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param shape shape parameter.
#' @param rate  rate parameter
#' @param lower.tail logical; if \code{TRUE}
#' @param log,log.p logical
#' @name gLindley
NULL
#> NULL
#' @return later
#' @importFrom stats uniroot
#' @importFrom lamW lambertWm1
#' @export
#'
#' @examples

#' @rdname gLindley
#' @export
dgLindley = function (x, shape = 1, rate, log = FALSE)
{
  denslog = log(shape) + (2 * log(rate)) - log1p(rate) + log1p(x) +
    ((shape - 1) * log1p(- ((1 + rate + (rate * x)) / (1 + rate)) * exp(- rate * x))) - (rate * x)
  denslog[x < 0] = -Inf
  denslog[(shape < 0) | (rate < 0)] = NaN
  denslog[(shape == 0) | (x == 0) | (rate == 0)] = -Inf
  if (log) {
    return((denslog))
  }
  else {
    return(exp((denslog)))
  }
}
#' @rdname gLindley
#' @export
pgLindley = function (q, shape = 1, rate, lower.tail = TRUE, log.p = FALSE)
{
  if (lower.tail) {
    if (log.p) {
      distfun = shape * log1p(- ((1 + rate + (rate * q)) / (1 + rate)) * exp(- rate * q))
      distfun[q <= 0] = -Inf
    }
    else {
      distfun = (1 - (((1 + rate + (rate * q)) / (1 + rate)) * exp(- rate * q))) ^ shape
      distfun[q <= 0] = 0
    }
  }
  else {
    if (log.p) {
      distfun = log1p(- (1 - (((1 + rate + (rate * q)) / (1 + rate)) * exp(- rate * q))) ^ shape)
      distfun[q <= 0] = 0
    }
    else {
      distfun = 1 - (1 - (((1 + rate + (rate * q)) / (1 + rate)) * exp(- rate * q))) ^ shape
      distfun[q <= 0] = 1
    }
  }
  distfun[(shape < 0) | (rate < 0)] = NaN
  return(distfun)
}
#' @rdname gLindley
#' @export
qgLindley = function (p, shape = 1, rate, lower.tail = TRUE, log.p = FALSE)
{
  if (lower.tail) {
    if (log.p) {
      qfun = - 1 - (1/rate) - ((1/rate) * lambertWm1(- (1 + rate) * exp(- 1 - rate) * (1 - (exp(p)) ^ (1/shape))))
      qfun[p > 0] = NaN
    }
    else {
      qfun = - 1 - (1/rate) - ((1/rate) * lambertWm1(- (1 + rate) * exp(- 1 - rate) * (1 - (p) ^ (1/shape))))
      qfun[p < 0] = NaN
      qfun[p > 1] = NaN
    }
  }
  else {
    if (log.p) {
      qfun = - 1 - (1/rate) - ((1/rate) * lambertWm1(- (1 + rate) * exp(- 1 - rate) * (1 - (1 - exp(p)) ^ (1/shape))))
      qfun[p > 0] = NaN
    }
    else {
      qfun = - 1 - (1/rate) - ((1/rate) * lambertWm1(- (1 + rate) * exp(- 1 - rate) * (1 - (1 - p) ^ (1/shape))))
      qfun[p < 0] = NaN
      qfun[p > 1] = NaN
    }
  }
  qfun[(shape < 0) | (rate < 0)] = NaN
  return(qfun)
}
#' @rdname gLindley
#' @export
rgLindley = function(n, shape = 1, rate)
{
  out = qgLindley(runif(n), shape = shape, rate = rate)
  return(out)
}
