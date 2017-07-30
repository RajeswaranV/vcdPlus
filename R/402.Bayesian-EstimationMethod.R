#' Bayes estimation with different hyper priors
#' @param x - Vector of positive integers
#' @param d - Number of divisions
#' @return A list of dataframes
#'  \item{BE.diff.prior.vol.df}{ Bayesian estimation Volume}
#'  \item{BE.diff.prior.LL.df}{  Dataframe of Bayesian prior Lower Limits}
#'  \item{BE.diff.prior.UL.df}{  Dataframe of Bayesian prior Upper Limits}
#'  \item{BE.diff.prior.mean.df}{  Dataframe of Bayesian prior mean}
#' @family Confidence Interval for Multinomial Proportion
#' @examples
#' x = c(56,72,73,59,62,87,68,99,98)
#' d=3
#' ci.BMDU(x,d)
#' @seealso
#'    \code{\link[rBeta2009]{rdirichlet} with additional n and shape paramaters. This function calls rBeta2009 function.}
#' @export
ci.BMDU <- function(x,d)
{
  if (missing(x)) stop("Categorical data  is missing")
  if (missing(d)) stop("'d' - the size of division is missing")
  if ((class(x) != "integer") &
      (class(x) != "numeric")) stop("x has to be numeric")
  if (d<1 || d>length(x))
    stop("Size of the division (d) should be less than the size of the input matrix")

  k=length(x)
  for(m in 1:k)
  {
    if(x[m]<0)
    {stop("Arguments must be non-negative integers")
    }
  }

    m=0
    l=0
    u=0
    diff=0
    s=sum(x)
    s1=floor(k/d)
    d1=runif(s1,0,1)###First half of the vector
    d2=runif(k-s1,1,2)###Second half of the vector
    a=c(d1,d2)
    p=x+a###Prior for Dirichlet
    dr=rdirichlet(10000, p)###Posterior
    for(j in 1:k)
    {
      l[j]=round(quantile(dr[,j],0.025),4)###Lower Limit
      u[j]=round(quantile(dr[,j],0.975),4)###Upper Limit
      m[j]=round(mean(dr[,j]),4)###Point Estimate
      diff[j]=u[j]-l[j]
    }
    p=prod(diff)
    BE.diff.prior.UL.df = data.frame(BE.diff.prior.UpLim=u)
    BE.diff.prior.LL.df = data.frame(BE.diff.prior.LowLim=l)
    BE.diff.prior.mean.df = data.frame(BE.diff.prior.Mean=m)
    BE.diff.prior.vol.df = data.frame(BE.diff.prior.Volume=p)

    result <- list(BE.diff.prior.vol.df,
                   BE.diff.prior.LL.df,
                   BE.diff.prior.UL.df,
                   BE.diff.prior.mean.df
                   )
    return(result)

}
