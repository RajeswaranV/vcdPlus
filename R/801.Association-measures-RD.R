#'  This function  calls pairwiseCI package and 3 methods from that are implemented with consistant output view.
#' @param ai - Numeric value representing True Positives
#' @param bi - Numeric value representing False Positives
#' @param ci - Numeric value representing False Negatives
#' @param di - Numeric value representing True Negatives
#' @param alp  - Numeric value representing required significance level
#' @details  Three methods from pairwiseCI package is implemented with consistant output view.
#' Methods implemented are Newcombes Hybrid score, Wald continuity corrected and Agresti Caffo.
#' @return A dataframe with
#'  \item{Method}{Three methods are implemented - Newcombes Hybrid score, Wald continuity corrected and Agresti Caffo}
#'  \item{Estimate}{The estimated value}
#'  \item{Lower.Limit}{The lower limit }
#'  \item{Upper.Limit}{The upper limit }
#' @family Association Measures methods
#' @seealso
#'    \code{\link[pairwiseCI]{Prop.diff}}  provides this functionality with additional parameters as well
#' @examples
#' ai=7; bi=11; ci=1; di=17; alp=0.05 # Example from [reference 1]
#' suppressWarnings(Association.Measures.RD(ai,bi,ci,di,alp))
#' @references
#' [1] Fagerland, M. W., & Newcombe, R. G. (2013).
#' Confidence intervals for odds ratio and relative risk based on the inverse hyperbolic sine transformation.
#' Statistics in medicine, 32(16), 2823-2836.
#' [2] Upton, G. J. (2016).
#' Categorical Data Analysis by Example.
#' John Wiley & Sons.
#' @export
Association.Measures.RD<-function(ai,bi,ci,di,alp)
{
  if ((class(ai) != "numeric") & (class(ai) != "integer") ||
      (class(bi) != "numeric")  & (class(bi) != "integer")   ||
      (class(ci) != "numeric")  & (class(ci) != "integer")  ||
      (class(di) != "numeric")  & (class(di) != "integer"))   stop("Inputs ai,bi,ci and di have to be numeric")
  if ((class(alp) != "numeric") & (class(alp) != "integer") || alp>1 || alp<0) stop("Alpha (alp) has to be between 0 and 1")

  RD1=Prop.diff(c(ai,bi), c(ci,di), conf.level=1-alp, alternative="two.sided", CImethod="NHS")
  RES_RD1=data.frame(Method="NewcombesHybridScore",Estimate=RD1$estimate, Lower.Limit=RD1$conf.int[1],Upper.Limit=RD1$conf.int[2])
  RD2=Prop.diff(c(ai,bi), c(ci,di), conf.level=1-alp, alternative="two.sided", CImethod="CC")
  RES_RD2=data.frame(Method="WaldCC" ,Estimate=RD2$estimate, Lower.Limit=RD2$conf.int[1],Upper.Limit=RD2$conf.int[2])
  RD3=Prop.diff(c(ai,bi), c(ci,di), conf.level=1-alp, alternative="two.sided", CImethod="AC")
  RES_RD3=data.frame(Method="AgrestiCaffo",Estimate=RD3$estimate, Lower.Limit=RD3$conf.int[1],Upper.Limit=RD3$conf.int[2])
  RES_RD=rbind(RES_RD1,RES_RD2,RES_RD3)

  return(RES_RD)

} # End of the function

