

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
#' @importFrom VGAM dgenray pgenray qgenray rgenray
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
  return(dgenray(x, scale = 1/rate, shape = shape, log = log))
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
  return(pgenray(q, scale = 1/rate, shape = shape, lower.tail = lower.tail, log.p = log.p))
}
#' @rdname GR
#' @export
qGR = function(p, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
{
  if (!missing(rate) && !missing(scale)) {
    if (abs(rate * scale - 1) < 1e-15)
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  rate = 1/scale
  return(qgenray(p, scale = 1/rate, shape = shape, lower.tail = lower.tail, log.p = log.p))
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
