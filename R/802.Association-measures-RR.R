#'  This function  calculates Risk Ratio and has 6 methods implemented leveraging pairwiseCI and metafor packags;
#'  Further inverse sinh base and adjusted are made available.
#' @param ai - Numeric value representing True Positives
#' @param bi - Numeric value representing False Positives
#' @param ci - Numeric value representing False Negatives
#' @param di - Numeric value representing True Negatives
#' @param alp  - Numeric value representing required significance level
#' @param c1 - Numeric value representing adjustment factor for inverse Sinh adjusted
#' @param c2 - Numeric value representing adjustment factor for inverse Sinh adjusted
#' @param c3 - Numeric value representing adjustment factor for inverse Sinh adjusted
#' @param c4 - Numeric value representing adjustment factor for inverse Sinh adjusted
#' @details  Two methods from pairwiseCI package and two methods from metafor package is implemented
#' along with inverese Sinh and inverse Sinh adjusted with consistant output view.
#' Methods implemented are Woolf, Gart, MI_NU, MOVER, Inverse Sinh and Inverse Sinh Adjusted.
#' @return A dataframe with
#'  \item{Method}{Six methods are implemented - Woolf, Gart, MI_NU, MOVER, ISH and ISH.Adjusted}
#'  \item{Estimate}{The estimated value}
#'  \item{Lower.Limit}{The lower limit }
#'  \item{Upper.Limit}{The upper limit }
#' @family Association Measures methods
#' @seealso
#'    \code{\link[pairwiseCI]{Prop.ratio}}  provides Score, MNScore, MOVER and GNC. This function call pairwiseCI functions.
#'    \code{\link[metafor]{escalc}}  provides comprenhisive functionality from which RR alone is used in this function
#' @examples
#'## Example taken from Fagerland and Newcombe [reference 1]
#' ai=7; bi=11; ci=1; di=17; alp=0.05
#' c1=0; c2=0.2; c3=0; c4=0.8
#' suppressWarnings(Association.Measures.RR(ai,bi,ci,di,alp, c1, c2, c3, c4))
#' @references
#' [1] Fagerland, M. W., & Newcombe, R. G. (2013).
#' Confidence intervals for odds ratio and relative risk based on the inverse hyperbolic sine transformation.
#' Statistics in medicine, 32(16), 2823-2836.
#' [2] Upton, G. J. (2016).
#' Categorical Data Analysis by Example.
#' John Wiley & Sons.
#' @export
Association.Measures.RR<-function(ai,bi,ci,di,alp, c1, c2, c3, c4)
{
  if ((class(ai) != "numeric") & (class(ai) != "integer") ||
      (class(bi) != "numeric")  & (class(bi) != "integer")   ||
      (class(ci) != "numeric")  & (class(ci) != "integer")  ||
      (class(di) != "numeric")  & (class(di) != "integer"))   stop("Inputs ai,bi,ci and di have to be numeric")
  if ((class(alp) != "numeric") & (class(alp) != "integer") || alp>1 || alp<0) stop("Alpha (alp) has to be between 0 and 1")
  if ((class(c1) != "numeric") & (class(c1) != "integer") ||
      (class(c2) != "numeric")  & (class(c2) != "integer")   ||
      (class(c3) != "numeric")  & (class(c3) != "integer")  ||
      (class(c4) != "numeric")  & (class(c4) != "integer"))   stop("Inputs c1,c2,c3 and c4 have to be numeric")

  rt1=ai+bi
  rt2=ci+di
  RR=(ai/rt1)/(ci/rt2)
  #####################metafor functions ##################################
  #Kart-Woolf
  RR1=escalc(measure = "RR",ai=ai,bi=bi,ci=ci,di=di)
  RES_RR1=summary(RR1,level=1-alp,transf=exp)
  R1=data.frame(Method="Woolf", Estimate=RES_RR1$yi,Lower.Limit=round(RES_RR1$ci.lb,4), Upper.Limit=round(RES_RR1$ci.ub,4))

  x=c(ai,bi,ci,di)
  #Adjusted log-Gart
  x2=x+0.5
  RR2=escalc(measure = "RR",ai=x2[1],bi=x2[2],ci=x2[3],di=x2[4])
  RES_RR2=summary(RR2,level=1-alp,transf=exp)
  R2=data.frame(Method="Gart",Estimate=RES_RR2$yi,Lower.Limit=round(RES_RR2$ci.lb,4), Upper.Limit=round(RES_RR2$ci.ub,4))

  #####################pairwiseCI functions ##################################
  RR3=Prop.ratio(c(ai,bi), c(ci,di), conf.level=1-alp, CImethod= "MNScore")
  R3=data.frame(Method="MI_NU",Estimate=RR3$estimate[[1]],Lower.Limit=round(RR3$conf.int[1],4), Upper.Limit=round(RR3$conf.int[2],4))
  RR4=Prop.ratio(c(ai,bi), c(ci,di), conf.level=1-alp, CImethod= "MOVER")
  R4=data.frame(Method="MOVER",Estimate=RR4$estimate[[1]],Lower.Limit=round(RR4$conf.int[1],4), Upper.Limit=round(RR4$conf.int[2],4))

  ##################### Own functions ##################################
  #Inverse Sinh
  z=qnorm(alp/2,lower.tail = FALSE)
  term=2*asinh((z/2)*sqrt((1/ai)+(1/ci)-(1/rt1)-(1/rt2)))
  LLRR_SiH= log(RR)-term
  ULRR_SiH=log(RR)+term
  R5=data.frame(Method="ISH",Estimate=RR,Lower.Limit=round(exp(LLRR_SiH),4), Upper.Limit=round(exp(ULRR_SiH),4))
  #Inverse Sinh - Adjusted
  RR_A=((ai+c1)/(rt1+c1+c2))/((ci+c1)/(rt2+c1+c2))
  z=qnorm(alp/2,lower.tail = FALSE)
  term_A=2*asinh((z/2)*sqrt((1/(ai+c3))+(1/(ci+c3))-(1/(rt1+c3+c4))-(1/(rt2+c3+c4))))
  LLRR_A_SiH= log(RR_A)-term_A
  ULRR_A_SiH=log(RR_A)+term_A
  R6=data.frame(Method="ISH.Adjusted",Estimate=RR_A,Lower.Limit=round(exp(LLRR_A_SiH),4), Upper.Limit=round(exp(ULRR_A_SiH),4))
  res_RR=rbind(R1,R2,R3,R4,R5,R6)
  rownames(res_RR)<-NULL

  return(res_RR)

} # End of the function

