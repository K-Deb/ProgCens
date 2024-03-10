#' Generating an adapted progressive first-failure dataset from a given dataset
#'
#' @param m number of desired failure
#' @param R Prefixed censoring scheme
#' @param Time Prefixed adaptation time
#' @param data a matrix, vector or a data frame. If the input is a vector, a matrix is created by taking the smallest factor of the length of the vector as columns.
#'
#' @return A list
#' @export
#'
#' @examples
Adap_gen = function(m, R, Time, data)
{
  if(is.data.frame(data)) data = as.matrix(data)
  is.prime <- function(n){
    ifelse(sum(n %% (1:n) == 0) > 2, FALSE, TRUE)
  }
  if(is.vector(data) & is.prime(length(data))) stop("The size of the given data is a prime number")
  if(is.vector(data) & ! is.prime(length(data)))
  {
    repeat
    {
      rem = 2
      if((length(data)%%rem) == 0)
      {
        break
      }else
      {
        rem = rem + 1
      }
    }
    data = matrix(data, ncol = rem)
  }else
  {
    data = data
  }
  Orig_data = data
  R1 = R
  group = nrow(data)
  newsamp = NA

  for(i in 1 : (m - 1))
  {
    newsamp[i] = min(data)
    data = data[- (which(data == newsamp[i], arr.ind = T)[1,1]),]
    if((newsamp[i] < Time) | (newsamp[i] == Time))
    {
      if(R[i] == 0)
      {
        data = data
      }else{
        data = data[- (sample(1 : nrow(data), R[i], replace = F)),]
      }
    }else
    {
      if(i != (m - 1))
      {
        R[i : (m - 1)] = rep(0, (m - i))
        for( j in (i + 1) : (m - 1))
        {
          newsamp[j] = min(data)
          data = data[- (which(data == newsamp[j], arr.ind = T)[1,1]),]
        }
      }else
      {
        R[i] = 0
        newsamp[i] = min(data)
        data = data[- (which(data == newsamp[m - 1], arr.ind = T)[1,1]),]
      }
      break
    }
  }
  newsamp[m] = min(data)
  if(newsamp[1] > Time)
  {
    R[m] = group - m
  }else
  {
    if(newsamp[m] < Time)
    {
      R[m] = group - m - sum(R[1 : (length(newsamp[!newsamp>Time]) - 1)])
    }else
    {
      R[m] = group - m - sum(R[1 : (length(newsamp[!newsamp>Time]))])
    }

  }
  R2 = R
  cen_samp = newsamp
  out = list("Orig" = Orig_data, "Prim" = R1, "Cen" = cen_samp, "Adap" = R2)
  return(out)
}
