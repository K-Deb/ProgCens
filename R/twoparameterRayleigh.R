#' The Two-parameter Rayleigh Distribution
#'
#' @param x,q Vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param location Location Parameter
#' @param scale Scale parameter
#' @param lower.tail logical; if \code{TRUE}
#' @param log,log.p logical
#' @name tworayl
NULL
#> NULL
#' @return
#' @export
#'
#' @examples

#' @rdname tworayl
#' @export
d2rayleigh <- function(x, location = 0, scale, log = F)
{
  denslog <- log(x - location) - (2 * log(scale)) - ((x - location) ^ 2/(2 * scale ^ 2))
  denslog[x < 0] <- -Inf
  denslog[(x < location) | (scale < 0)] <- NaN
  denslog[scale == 0] <- -Inf
  denslog[(scale == 0) & (x == 0)] <- Inf
  if (log) {
    return((denslog))
  }
  else {
    return(exp((denslog)))
  }
}
#' @rdname tworayl
#' @export
p2rayleigh <- function(q, location = 0, scale, lower.tail = T, log.p = F)
{
  if (lower.tail) {
    if (log.p) {
      distfun <- log1p(- exp(- (q - location) ^ 2/(2 * scale ^ 2)))
      distfun[q <= 0] = -Inf
    }
    else {
      distfun = - expm1(- ((q - location) ^ 2/(2 * scale ^ 2)))
      distfun[q <= 0] = 0
    }
  }
  else {
    if (log.p) {
      distfun = - ((q - location) ^ 2/(2 * scale ^ 2))
      distfun[q <= 0] = 0
    }
    else {
      distfun = exp(- ((q - location) ^ 2/(2 * scale ^ 2)))
      distfun[q <= 0] = 1
    }
  }
  distfun[(q < location) | (scale < 0)] = NaN
  return(distfun)
}
#' @rdname tworayl
#' @export
q2rayleigh <- function(p, location = 0, scale, lower.tail = T, log.p = F)
{
  if (lower.tail) {
    if (log.p) {
      qfun = location + sqrt(- (2 * (scale ^ 2)) * log1p(- exp(p)))
      qfun[p > 0] = NaN
    }
    else {
      qfun = location + sqrt(- (2 * (scale ^ 2)) * log1p(- p))
      qfun[p < 0] = NaN
      qfun[p > 1] = NaN
    }
  }
  else {
    if (log.p) {
      qfun = location + sqrt(- (2 * (scale ^ 2)) * log(exp(p)))
      qfun[p > 0] = NaN
    }
    else {
      qfun = location + sqrt(- (2 * (scale ^ 2)) * log(p))
      qfun[p < 0] = NaN
      qfun[p > 1] = NaN
    }
  }
  qfun[(scale < 0)] = NaN
  return(qfun)
}
#' @rdname tworayl
#' @export
r2rayleigh = function(n, location = 0, scale)
{
  out = q2rayleigh(runif(n), location = location, scale = scale)
  return(out)
}
