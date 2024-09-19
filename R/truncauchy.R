#' The Truncated Cauchy Distribution
#'
#' @param x,q Vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param location location parameter.
#' @param scale  scale parameter
#' @param upper Upper truncation point
#' @param lower Lower truncation point
#' @param lower.tail logical; if \code{TRUE}
#' @param log,log.p logical
#' @name trunc
NULL
#> NULL
#' @return later
#' @importFrom stats dcauchy pcauchy qcauchy
#' @export
#'
#' @examples

#' @rdname trunc
#' @export
dtruncauchy <- function(x, location = 0, scale = 1, lower = 0, upper = Inf, log = F)
{
  denslog = dcauchy(x, location = location, scale = scale, log = T) - log(pcauchy(upper, location = location, scale = scale) - pcauchy(lower, location = location, scale = scale))
  denslog[x < 0] <- -Inf
  if(log)
  {
    return((denslog))
  }else{
    return(exp((denslog)))
  }
}
#' @rdname trunc
#' @export
ptruncauchy <- function(q, location = 0, scale = 1, lower = 0, upper = Inf, lower.tail = T, log.p = F)
{
  if(lower.tail)
  {
    if(log.p)
    {
      distfun = log(pcauchy(q, location = location, scale = scale) - pcauchy(lower, location = location, scale = scale)) - log(pcauchy(upper, location = location, scale = scale) - pcauchy(lower, location = location, scale = scale))
      distfun[q <= 0] = -Inf
    }else
    {
      distfun = (pcauchy(q, location = location, scale = scale) - pcauchy(lower, location = location, scale = scale))/(pcauchy(upper, location = location, scale = scale) - pcauchy(lower, location = location, scale = scale))
      distfun[q <= 0] = 0
    }
  }else{
    if(log.p)
    {
      distfun = log(pcauchy(upper, location = location, scale = scale) - pcauchy(q, location = location, scale = scale)) - log(pcauchy(upper, location = location, scale = scale) - pcauchy(lower, location = location, scale = scale))
      distfun[q <= 0] = 0
    }else{
      distfun = (pcauchy(upper, location = location, scale = scale) - pcauchy(q, location = location, scale = scale))/(pcauchy(upper, location = location, scale = scale) - pcauchy(lower, location = location, scale = scale))
      distfun[q <= 0] = 1
    }
  }
  return(distfun)
}
#' @rdname trunc
#' @export
qtruncauchy <- function(p, location = 0, scale = 1, lower = 0, upper = Inf, lower.tail = T, log.p = F)
{
  if (lower.tail) {
    if (log.p) {
      qfun = qcauchy(pcauchy(lower, location = location, scale = scale) + (exp(p) * (pcauchy(upper, location = location, scale = scale) - pcauchy(lower, location = location, scale = scale))), location = location, scale = scale)
      qfun[p > 0] = NaN
    }
    else {
      qfun = qcauchy(pcauchy(lower, location = location, scale = scale) + (p * (pcauchy(upper, location = location, scale = scale) - pcauchy(lower, location = location, scale = scale))), location = location, scale = scale)
      qfun[p < 0] = NaN
      qfun[p > 1] = NaN
    }
  }
  else {
    if (log.p) {
      qfun = qcauchy(pcauchy(upper, location = location, scale = scale) - (exp(p) * (pcauchy(upper, location = location, scale = scale) - pcauchy(lower, location = location, scale = scale))), location = location, scale = scale)
      qfun[p > 0] = NaN
    }
    else {
      qfun = qcauchy(pcauchy(upper, location = location, scale = scale) - (p * (pcauchy(upper, location = location, scale = scale) - pcauchy(lower, location = location, scale = scale))), location = location, scale = scale)
      qfun[p < 0] = NaN
      qfun[p > 1] = NaN
    }
  }
  return(qfun)
}
#' @rdname trunc
#' @export
rtruncauchy <- function(n, location = 0, scale = 1)
{
  out = qtruncauchy(runif(n), location = location, scale = scale)
  return(out)
}
