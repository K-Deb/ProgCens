#' The Log Cauchy Distribution
#'
#' @param x,q Vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param location location parameter.
#' @param scale  scale parameter
#' @param lower.tail logical; if \code{TRUE}
#' @param log,log.p logical
#' @name lc
NULL
#> NULL
#' @return later
#' @importFrom stats dcauchy pcauchy qcauchy
#' @export
#'
#' @examples

#' @rdname lc
#' @export
dlcauchy = function(x, location = 0, scale = 1, log = FALSE)
{
  denslog = dcauchy(log(x), location = location, scale = scale,
                    log = TRUE) - log(x)
  denslog[x < 0] = -Inf
  denslog[x == 0] = Inf
  if(log)
  {
    return((denslog))
  }else{
    return(exp((denslog)))
  }
}
#' @rdname lc
#' @export
plcauchy = function(q, location = 0, scale = 1, lower.tail = T, log.p = FALSE)
{
  if(lower.tail)
  {
    if(log.p)
    {
      distfun = pcauchy(log(q), location = location, scale = scale, log.p = T)
      distfun[q <= 0] = -Inf
    }else
    {
      distfun = pcauchy(log(q), location = location, scale = scale)
      distfun[q <= 0] = 0
    }
  }else{
    if(log.p)
    {
      distfun = pcauchy(log(q), location = location, scale = scale, lower.tail = F, log.p = T)
      distfun[q <= 0] = 0
    }else{
      distfun = pcauchy(log(q), location = location, scale = scale, lower.tail = F, log.p = F)
      distfun[q <= 0] = 1
    }
  }
  return(distfun)
}
#' @rdname lc
#' @export
qlcauchy = function(p, location = 0, scale = 1, lower.tail = T, log.p = FALSE)
{
  if(lower.tail)
  {
    if(log.p)
    {
      qfun = exp(qcauchy(p, location = location, scale = scale, log.p = T))
      qfun[p > 0] = NaN
    }else
    {
      qfun = exp(qcauchy(p, location = location, scale = scale, log.p = F))
      qfun[p < 0] = NaN
      qfun[p > 1] = NaN
    }
  }else{
    if(log.p)
    {
      qfun = exp(qcauchy(p, location = location, scale = scale, lower.tail = F, log.p = T))
      qfun[p > 0] = NaN
    }else{
      qfun = exp(qcauchy(p, location = location, scale = scale, lower.tail = F, log.p = F))
      qfun[p < 0] = NaN
      qfun[p > 1] = NaN
    }
  }
  return(qfun)
}
#' @rdname lc
#' @export
rlcauchy = function(n, location = 0, scale = 1)
{
  out = qlcauchy(runif(n), location = location, scale = scale)
  return(out)
}
