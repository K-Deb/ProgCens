#' The Generalized Rayleigh Distribution
#'
#' @param x,q Vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param shape shape parameter.
#' @param scale scale parameter.
#' @param rate  an alternative way to specify the scale.
#' @param lower.tail logical; if \code{TRUE}
#' @param log,log.p logical
#' @name GR
NULL
#> NULL
#' @return
#' @export
#'
#' @examples

#' @rdname GR
#' @export
dgenrayl = function(x, shape, rate = 1, scale = 1/rate, log = FALSE)
{
  if (!missing(rate) && !missing(scale)) {
    if (abs(rate * scale - 1) < 1e-15)
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  rate = 1/scale
  denslog = log(2) + log(shape) + (2 * log(rate)) + log(x) - (rate * x)^2 + ((shape - 1) * log(- expm1(- (rate * x)^2)))
  denslog[x < 0] = -Inf
  denslog[(shape < 0) | (rate < 0)] = NaN
  denslog[shape == 0] = -Inf
  denslog[(shape == 0) & ((rate == 0) | (x == 0))] = Inf
  denslog[(x == 0) & (shape < 1)] = Inf
  denslog[(x == 0) & (shape > 1)] = -Inf
  if(log)
  {
    return((denslog))
  }else{
    return(exp((denslog)))
  }

}
#' @rdname GR
#' @export
pgenrayl = function(q, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
{
  if (!missing(rate) && !missing(scale)) {
    if (abs(rate * scale - 1) < 1e-15)
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  rate = 1/scale
  if(lower.tail)
  {
    if(log.p)
    {
      distfun = shape * log(- expm1(- (rate * q)^2))
      distfun[q <= 0] = -Inf
    }else
    {
      distfun = (- expm1(- (rate * q)^2)) ^ shape
      distfun[q <= 0] = 0
    }
  }else{
    if(log.p)
    {
      distfun = log(1 - (- expm1(- (rate * q)^2)) ^ shape)
      distfun[q <= 0] = 0
    }else{
      distfun = 1 - (- expm1(- (rate * q)^2)) ^ shape
      distfun[q <= 0] = 1
    }
  }
  distfun[(shape < 0) | (rate < 0)] = NaN
  return(distfun)
}
#' @rdname GR
#' @export
qgenrayl = function(p, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
{
  if (!missing(rate) && !missing(scale)) {
    if (abs(rate * scale - 1) < 1e-15)
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  rate = 1/scale
  if(lower.tail)
  {
    if(log.p)
    {
      qfun = (1/rate) * sqrt(- log1p(- (exp(p) ^ (1/shape))))
      qfun[p > 0] = NaN
    }else
    {
      qfun = (1/rate) * sqrt(- log1p(- ((p) ^ (1/shape))))
      qfun[p < 0] = NaN
      qfun[p > 1] = NaN
    }
  }else{
    if(log.p)
    {
      qfun = (1/rate) * sqrt(- log1p(- ((- expm1(p)) ^ (1/shape))))
      qfun[p > 0] = NaN
    }else{
      qfun = (1/rate) * sqrt(- log1p(- (1 - p) ^ (1/shape)))
      qfun[p < 0] = NaN
      qfun[p > 1] = NaN
    }
  }
  qfun[(shape < 0) | (rate < 0)] = NaN
  return(qfun)
}
#' @rdname GR
#' @export
rgenrayl = function(n, shape, rate = 1, scale = 1/rate)
{
  if (!missing(rate) && !missing(scale)) {
    if (abs(rate * scale - 1) < 1e-15)
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  rate = 1/scale
  out = qgenrayl(runif(n), shape = shape, rate = rate)
  return(out)
}
