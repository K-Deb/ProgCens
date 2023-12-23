#' Simulating a progressively type-II censored sample from standard uniform distribution
#'
#' @param n Total number of items placed on test at time zero
#' @param m Number of required failure
#' @param R Prefixed censoring scheme
#' @param prob Probability of removal
#' @param set_seed Random number generation seed
#'
#' @return A dataframe with progressively type-II censored uniform sample and the censoring scheme
#' @importFrom stats rbinom
#' @importFrom stats runif
#' @export
#'
#' @examples
rpc2unif = function(n, m = NULL, R = NULL, prob = NULL, set_seed = NULL)
{
  if(!is.null(set_seed)) set.seed(set_seed)
  if(is.null(m) & !is.null(R)) m = length(R)
  if(!is.null(R) & sum(R) != n-m) stop(paste("Total number of removed units must be", n-m))
  if(!is.null(prob) & is.null(R))
  {
    R[1] = rbinom(1, size = n-m, prob = prob)
    for(i in 2:(m-1)) R[i] = rbinom(1, size=n-m-sum((R[1:(i-1)])), prob = prob)
    if((n-m-sum(R[1:(m-1)]))>0)
    {
      R[m] = n-m-sum(R[1:(m-1)])
    }else{
      R[m] = 0
    }
  }
  V = vector(mode = "numeric",length = m)
  U = vector(mode = "numeric",length = m)
  W = runif(m,0,1)
  for(i in 1:m) V[i] = W[i]^(1/(i+sum((R[(m-i+1):m]))))

  for (i in 1:m) U[i] = 1-prod((V[(m-i+1):m]))
  return(data.frame(U = U, R = round(R)))
}

