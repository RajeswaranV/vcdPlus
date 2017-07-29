#'  This partition calcluation
#' @param a - description
#' @param b - description
#' @param ie - description
#' @param je - description
#' @param d - description
#' @return A dataframe with
#'  \item{x1}{Sequence}
#' @family Internal methods
#' @keywords internal
gmat<-function(a,b,ie,je,d){
  x1=matrix(0,2,2)
  for(i in 2:ie)
  {
    for(j in 2:je)
    {
      x1[1,1]=sum(d[which(a<i),which(b<j)])
      x1[1,2]=sum(d[which(a<i),j])
      x1[2,1]=sum(d[i,which(b<j)])
      x1[2,2]=d[i,j]
    }
  }
  return(x1)

}
#'  This G squared metric
#' @param catd - description
#' @return A list with
#'  \item{x1}{Sequence}
#' @family Internal methods
#' @keywords internal
GSquared <- function(catd)
{
  # computes the LR (G square) test of independence
  chsq <- chisq.test(catd)
  est <- chsq$expected
  df <- chsq$parameter
  LRterm <- 2*catd*log(catd/est)
  LRterm[catd==0] <- 0
  G2 <- sum(LRterm)
  p <- 1-pchisq(G2, df)
  LR_RES=sign(catd-est)*sqrt(abs(LRterm))
  return(list(GSq=G2, df=df, p.value=p,GsqRes=LR_RES))
}
