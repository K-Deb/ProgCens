#' The Chen Distribution
#'
#' @param x,q Vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param shape shape parameter.
#' @param rate  rate parameter
#' @param lower.tail logical; if \code{TRUE}
#' @param log,log.p logical
#' @name Chen
NULL
#> NULL
#' @return
#' @export
#'
#' @examples

#' @rdname Chen
#' @export
dchen = function(x, shape, rate = 1, log = FALSE)
{
  denslog = log(rate) + log(shape) + ((shape - 1) * log(x)) + ((rate * (- expm1(x ^ shape))) + (x ^ shape))
  denslog[x < 0] = -Inf
  denslog[(shape < 0) | (rate < 0)] = NaN
  denslog[shape == 0] = -Inf
  denslog[(shape == 0) & (x == 0)] = Inf
  denslog[(x == 0) & (shape < 1)] = Inf
  denslog[(x == 0) & (shape > 1)] = -Inf
  if(log)
  {
    return((denslog))
  }else{
    return(exp((denslog)))
  }
}
#' @rdname Chen
#' @export
pchen = function(q, shape, rate = 1, lower.tail = TRUE, log.p = FALSE)
{
  if(lower.tail)
  {
    if(log.p)
    {
      distfun = log1p(- exp(rate * (- expm1(q ^ shape))))
      distfun[q <= 0] = -Inf
    }else
    {
      distfun = 1 + (- exp(rate * (- expm1(q ^ shape))))
      distfun[q <= 0] = 0
    }
  }else{
    if(log.p)
    {
      distfun = rate * (- expm1(q ^ shape))
      distfun[q <= 0] = 0
    }else{
      distfun = exp(rate * (- expm1(q ^ shape)))
      distfun[q <= 0] = 1
    }
  }
  distfun[(shape < 0) | (rate < 0)] = NaN
  return(distfun)
}
#' @rdname Chen
#' @export
qchen = function(p, shape, rate = 1, lower.tail = TRUE, log.p = FALSE)
{
  if(lower.tail)
  {
    if(log.p)
    {
      qfun = (log1p(- ((1/rate) * log1p(- exp(p))))) ^ (1/shape)
      qfun[p > 0] = NaN
    }else
    {
      qfun = (log1p(- ((1/rate) * log1p(- p)))) ^ (1/shape)
      qfun[p < 0] = NaN
      qfun[p > 1] = NaN
    }
  }else{
    if(log.p)
    {
      qfun = (log1p(- (p/rate))) ^ (1/shape)
      qfun[p > 0] = NaN
    }else{
      qfun = (log1p(- ((1/rate) * log(p)))) ^ (1/shape)
      qfun[p < 0] = NaN
      qfun[p > 1] = NaN
    }
  }
  qfun[(shape < 0) | (rate < 0)] = NaN
  return(qfun)
}
#' @rdname Chen
#' @export
rchen = function(n, shape, rate = 1)
{
  out = qchen(runif(n), shape = shape, rate = rate)
  return(out)
}
