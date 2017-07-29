#'  This function  implements 9 methods for testing 2x2 tables and outputs the p-value
#'  In hypergemetric notation out of N = m+n items r size is chosen; x is observed success
#' @param x - Numeric value representing no of success in attribute 1.
#' @param m - Numeric value representing Total no of attribute 1
#' @param n - Numeric value representing Total no of attribute 2
#' @param r - Numeric value representing Total number of success
#' @details  Nine methods - IRWIN, Central, Blaker.exact, MP_IR, MP_BK, Exact.less,
#' Exact.greater, Pearson.Chisq, Pearson.Chisq.Yates.cc are implimented from
#' exact2x2 package (7 methods) and 2 methods from base stats
#' @return A dataframe with
#'  \item{IRWIN}{Fisher exact test with 'minlike'}
#'  \item{Central}{Fisher exact test with 'central'}
#'  \item{Blaker.exact}{Blaker exact test }
#'  \item{MP_IR}{Fisher exact test with 'minlike', hypergeometric}
#'  \item{MP_BK}{Blaker exact test, hypergeometric}
#'  \item{Exact.less}{Fisher exact test if one sided is less}
#'  \item{Exact.greater}{Fisher exact test if one sided is greater}
#'  \item{Pearson.Chisq}{Pearson Chisquared test}
#'  \item{Pearson.Chisq.Yates.cc}{Pearson Chisquared test with Yates continuity correction}
#' @family Test methods
#' @seealso
#'    \code{\link[exact2x2]{fisher.exact}{blaker.exact}} with additional parameters. This function calls exact2x2 functions.
#'    \code{\link[stats]{chisq.test}} with additional parameters. This function calls base stats function.
#' @examples
#'## Example taken from -p28 Ex3.1 of Upton, G. J. (2016).
#'## Values are x,m,n,r: related to a hypergeometric distribution (?dhyper, for further help)
#' x=1; m=3; n=10; r=4
#' exact2x2tests.hypergeom(x,m,n,r)
#' @references
#' [1] Upton, G. J. (2016).
#' Categorical Data Analysis by Example.
#' John Wiley & Sons.
#' @export
exact2x2tests.hypergeom<-function(x,m,n,r)
{
  if ((class(x) != "numeric") & (class(x) != "integer") ||
      (class(m) != "numeric")  & (class(m) != "integer")   ||
      (class(n) != "numeric")  & (class(n) != "integer")  ||
      (class(r) != "numeric")  & (class(r) != "integer"))   stop("Inputs x,m,n and r have to be numeric")
  if (m < 0 || n<0) stop("m and n have to be  greater or equal to zero")
  if (x > m || x<0) stop("x has to be  greater or equal to zero and less than or equal to m")
  if (r >= n+x || r<0) stop("r has to be  greater or equal to zero and less than or equal to n+x")

  dat=matrix(c(x,m-x,r-x,n-r+x),2,2,byrow=TRUE)
  IRWIN=fisher.exact(dat,tsmethod="minlike")$p.value
  CDT=fisher.exact(dat,tsmethod="central")$p.value
  BKR=blaker.exact(dat)$p.value
  MP_IR=fisher.exact(dat,tsmethod="minlike")$p.value-0.5*dhyper(x, m, n, r, log = FALSE)
  MP_BK=blaker.exact(dat)$p.value-0.5*dhyper(x, m, n, r, log = FALSE)
  EXA_L=fisher.exact(dat,alternative = "less")$p.value
  EXA_G=fisher.exact(dat,alternative = "greater")$p.value
  PEAR=chisq.test(dat)$p.value
  PEAR_YCC=chisq.test(dat,correct = TRUE)$p.value
  ans.tests.df=data.frame(IRWIN=round(IRWIN,4),
                    Central=round(CDT,4),
                    Blaker.exact=round(BKR,4),
                    MP_IR= round(MP_IR,4),
                    MP_BK= round(MP_BK,4),
                    Exact.less=round(EXA_L,4),
                    Exact.greater=round(EXA_G,4),
                    Pearson.Chisq=round(PEAR,4),
                    Pearson.Chisq.Yates.cc=round(PEAR_YCC,4))
  return(ans.tests.df)
}

####################################################################################################
#'  This function  implements 9 methods for testing 2x2 tables and outputs the p-value
#'  This is the regular function taking the 2x2 table as input
#' @param ai - Numeric value representing True Positives
#' @param bi - Numeric value representing False Positives
#' @param ci - Numeric value representing False Negatives
#' @param di - Numeric value representing True Negatives
#' @details  Nine methods - IRWIN, Central, Blaker.exact, MP_IR, MP_BK, Exact.less,
#' Exact.greater, Pearson.Chisq, Pearson.Chisq.Yates.cc are implimented from
#' exact2x2 package (7 methods) and 2 methods from base stats
#' @return A dataframe with
#'  \item{IRWIN}{Fisher exact test with 'minlike'}
#'  \item{Central}{Fisher exact test with 'central'}
#'  \item{Blaker.exact}{Blaker exact test }
#'  \item{MP_IR}{Fisher exact test with 'minlike', hypergeometric}
#'  \item{MP_BK}{Blaker exact test, hypergeometric}
#'  \item{Exact.less}{Fisher exact test if one sided is less}
#'  \item{Exact.greater}{Fisher exact test if one sided is greater}
#'  \item{Pearson.Chisq}{Pearson Chisquared test}
#'  \item{Pearson.Chisq.Yates.cc}{Pearson Chisquared test with Yates continuity correction}
#' @family Test methods
#' @seealso
#'    \code{\link[exact2x2]{fisher.exact}{blaker.exact}} with additional parameters. This function calls exact2x2 functions.
#'    \code{\link[stats]{chisq.test}} with additional parameters. This function calls base stats function.
#' @examples
#'## Example taken from -p28 Ex3.1 of Upton, G. J. (2016).
#' ai=1; bi=3; ci=10; di=4
#' exact2x2tests.regular(ai,bi,ci,di)
#' @references
#' [1] Upton, G. J. (2016).
#' Categorical Data Analysis by Example.
#' John Wiley & Sons.
#' @export
exact2x2tests.regular<-function(ai,bi,ci,di)
{
  if ((class(ai) != "numeric") & (class(ai) != "integer") ||
      (class(bi) != "numeric")  & (class(bi) != "integer")   ||
      (class(ci) != "numeric")  & (class(ci) != "integer")  ||
      (class(di) != "numeric")  & (class(di) != "integer"))   stop("Inputs ai,bi,ci and di have to be numeric")

  x=c(ai,bi,ci,di)
  dat=matrix(x,2,2)
  IRWIN=fisher.exact(dat,tsmethod="minlike")$p.value
  CDT=fisher.exact(dat,tsmethod="central")$p.value
  BKR=blaker.exact(dat)$p.value
  MP_IR=fisher.exact(dat,tsmethod="minlike")$p.value-0.5*dhyper(dat[1,1], dat[1,1]+dat[1,2], dat[2,1]+dat[2,2], dat[1,1]+dat[2,1], log = FALSE)
  MP_BK=blaker.exact(dat)$p.value-0.5*dhyper(dat[1,1], dat[1,1]+dat[1,2], dat[2,1]+dat[2,2], dat[1,1]+dat[2,1], log = FALSE)
  #if one sided is less
  EXA_L=fisher.exact(dat,alternative = "less")$p.value
  #if one sided is less
  EXA_G=fisher.exact(dat,alternative = "greater")$p.value
  PEAR=chisq.test(dat)$p.value
  PEAR_YCC=chisq.test(dat,correct = TRUE)$p.value
  ans.tests.df=data.frame(IRWIN=round(IRWIN,4),
                          Central=round(CDT,4),
                          Blaker.exact=round(BKR,4),
                          MP_IR= round(MP_IR,4),
                          MP_BK= round(MP_BK,4),
                          Exact.less=round(EXA_L,4),
                          Exact.greater=round(EXA_G,4),
                          Pearson.Chisq=round(PEAR,4),
                          Pearson.Chisq.Yates.cc=round(PEAR_YCC,4))
  return(ans.tests.df)
}
