#' Non-parametric estimate of the acceleration factor for $BC_a$ confidence intervals
#'
#' @param bootdata A numeric vector or matrix of bootstrap/ jackknife samples
#'
#' @return later
#' @export
#'
#' @examples
a_hat = function(bootdata)
{
  acc = function(x)
  {
    (1/6) * (sum((x - mean(x))^3))/
      ((sum((x - mean(x))^2))^(3/2))
  }
  if(is.matrix(bootdata))
  {
    OUT = apply(bootdata, MARGIN = 2, FUN = acc)
  }else
  {
    OUT = acc(bootdata)
  }

  return(OUT)
}
