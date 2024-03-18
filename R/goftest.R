#' Approximate Goodness-of-fit test for progressive first-failure censored data
#'
#' @param data A 2 by m matrix, containing failure times in first row and the progressive first-failure censoring scheme in second row
#' @param CDF Cumulative distribution function
#' @param k number of units per group
#' @param norm_approx logical. Whether automatically calculate approximated goodness of fit test statistics
#' @param ... further arguments to be passed to CDF. Typically this is a vector of parameter estimates
#'
#' @return later
#' @importFrom stats qnorm dnorm pnorm optim sd
#' @export
#'
#' @examples
AdapGoF = function(data, CDF, k = 1, norm_approx = T, ...)
{
  if (! is.matrix(data) & (nrow(data) != 2))
    stop("Data must be matrix with 2 rows")
  stopifnot("CDF must be a distribution function" = is.function(CDF))
  S = data[1,]
  m = length(S)
  R = data[2,]
  n = sum(R) + m
  if (norm_approx)
  {
    if (as.character(substitute(CDF)) != "pnorm")
    {
      Y = qnorm(CDF(S, ...))
    }else
    {
      Y = S
    }

    negLLY = function(theta, x)
    {
      dd = (k * (R + 1)) - 1
      LL = sum(dnorm(x, theta[1], theta[2], log = T), na.rm = T)+
        sum(dd * pnorm(x, theta[1], theta[2], lower.tail = F, log.p = T), na.rm = T)
      return(-LL)
    }
    mleY = optim(c(mean(Y[!is.infinite(Y)], na.rm = T), sd(Y[!is.infinite(Y)], na.rm = T)), negLLY, x = Y)$par
    U = pnorm((Y - mleY[1])/mleY[2])
  }else
  {
    U = CDF(S, ...)
  }

  stat_A <- stat_C <- Alpha <- rep(NA, m - 1)
  for (i in 1 : (m - 1)) {
    Prod1 <- 1
    for (j in (m - i + 1) : m) Prod1 <- Prod1 * beta((j + sum(R[(m - j + 1) : (m)]) + (1/k)), 1)/beta((j + sum(R[(m - j + 1) : (m)])), 1)
    Alpha[i] <- 1 - Prod1
    stat_A[i] <- Alpha[i]^2 * log(U[i + 1] * (1 - U[i])/(U[i] * (1 - U[i + 1]))) + 2 * Alpha[i] * log((1 - U[i + 1])/(1 - U[i]))
    stat_C[i] <- Alpha[i] * (U[i + 1] - U[i]) * (Alpha[i] - U[i + 1] - U[i])
  }
  Anderson <- n * sum(stat_A) - n * log(1 - U[m]) - n * U[m]
  Cramer <- n * sum(stat_C) + n/3 * U[m]^3
  Prod2 <- 1
  for (j in 1:m) Prod2 <- Prod2 * beta((j + sum(R[(m - j + 1):(m)]) + (1/k)), 1)/beta((j + sum(R[(m - j + 1):(m)])), 1)

  KS <- max(max(c(Alpha, 1-Prod2)-U),max(U-c(0,Alpha)))
  return(list(KS = KS, CVM = Cramer, AD = Anderson))
}
