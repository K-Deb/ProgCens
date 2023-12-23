rprog = function(n, m, R, p = NULL, QF, ...)
{
  stopifnot("QF must be a quantile function" = is.function(QF))
  if(sum(R)!=n-m) stop("The provided censoring scheme does not conform to a progressive censoring design")

}
