#'  This function  implements 10 methods for OR confidence interaval - Woolf, Gart, Agresti,
#'  MOVER, MI_NU, CornEX, CornMP, BaPiEx, BaPiMP, Inverse Sinh and Inverse Sinh Adjusted
#' @param ai - Numeric value representing True Positives
#' @param bi - Numeric value representing False Positives
#' @param ci - Numeric value representing False Negatives
#' @param di - Numeric value representing True Negatives
#' @param alp  - Numeric value representing required significance level
#' @param e1 - Numeric value representing adjustment factor for inverse Sinh adjusted
#' @param e2 - Numeric value representing adjustment factor for inverse Sinh adjusted
#' @details  Ten methods - Woolf, Gart, Agresti, MOVER, MI_NU, CornEX, CornMP,
#' BaPiEx, BaPiMP, Inverse Sinh, Inverse Sinh Adjusted are implimented from
#' ORCI package (9 methods) and PropCIs package (1 method-Miettinen-Nurminen).
#' @return A dataframe with
#'  \item{Method}{Ten methods are implemented - Woolf, Gart, Agresti, MOVER, MI_NU, CornEX,
#'  CornMP, BaPiEx, BaPiMP, Inverse Sinh, Inverse Sinh Adjusted}
#'  \item{Estimate}{The estimated value}
#'  \item{Lower.Limit}{The lower limit }
#'  \item{Upper.Limit}{The upper limit }
#' @family Association Measures methods
#' @seealso
#'    \code{\link[ORCI]{Woolf.CI}{Gart.CI}{Agrestiind.CI}{MOVER.CI}{Cornfieldexact.CI}{Cornfieldmidp.CI}{BPexact.CI}{BPmidp.CI}{Invsinh.CI}} with additional parameters. This function calls ORCI functions.
#'    \code{\link[PropCIs]{orscoreci}} with additional parameters. This function calls orscoreci function.
#' @examples
#' ai=7; bi=11; ci=1; di=17; alp=0.05; e1=0.5; e2=0.5
#' Association.Measures.OR(ai, bi, ci, di, alp, e1, e2)
#' @references
#' [1] Fagerland, M. W., & Newcombe, R. G. (2013).
#' Confidence intervals for odds ratio and relative risk based on the inverse hyperbolic sine transformation.
#' Statistics in medicine, 32(16), 2823-2836.
#' [2] Upton, G. J. (2016).
#' Categorical Data Analysis by Example.
#' John Wiley & Sons.
#' @export
Association.Measures.OR<-function(ai,bi,ci,di,alp, e1, e2)
{
  if ((class(ai) != "numeric") & (class(ai) != "integer") ||
      (class(bi) != "numeric")  & (class(bi) != "integer")   ||
      (class(ci) != "numeric")  & (class(ci) != "integer")  ||
      (class(di) != "numeric")  & (class(di) != "integer"))   stop("Inputs ai,bi,ci and di have to be numeric")
  if ((class(alp) != "numeric") & (class(alp) != "integer") || alp>1 || alp<0) stop("Alpha (alp) has to be between 0 and 1")
  if ((class(e1) != "numeric") & (class(e1) != "integer") || e1>1 || e1<0) stop("Inverse Sinh limit (e1) has to be between 0 and 1")
  if ((class(e2) != "numeric") & (class(e2) != "integer") || e2>1 || e2<0) stop("Inverse Sinh limit (e2) has to be between 0 and 1")

  rt1=ai+bi
  rt2=ci+di
  ct1=ai+ci
  ct2=bi+di
  x=c(ai,bi,ci,di)
  N=sum(x)
  ODR=(ai/bi)/(ci/di)
  x1=x+0.5
  ODR1=(x1[1]/x1[2])/(x1[3]/x1[4])
  x2=c(ai+(rt1*ct1/N^2),bi+(rt1*ct2/N^2),ci+(rt2*ct1/N^2),di+(rt2*ct2/N^2))
  ODR2=(x2[1]/x2[2])/(x2[3]/x2[4])
  R11=Woolf.CI(ai, rt1, ci, rt2, conf = 1-alp)
  R1=data.frame(Method="Woolf", Estimate=ODR,Lower.Limit=round(R11[1],4), Upper.Limit=round(R11[2],4))

  R21=Gart.CI(ai, rt1, ci, rt2,conf = 1-alp)
  R2=data.frame(Method="Gart", Estimate=(ODR1),Lower.Limit=round(R21[1],4), Upper.Limit=round(R21[2],4))

  R31=Agrestiind.CI(x2[1], rt1, x2[3], rt2, conf = 1-alp)
  R3=data.frame(Method="Agresti", Estimate=(ODR2),Lower.Limit=round(R31[1],4), Upper.Limit=round(R31[2],4))

  R41=MOVER.CI(ai, rt1, ci, rt2, conf = 1-alp)
  R4=data.frame(Method="MOVER", Estimate=(ODR),Lower.Limit=round(R41[1],4), Upper.Limit=round(R41[2],4))

  R51=orscoreci(ai, rt1, ci, rt2, conf.level = 1-alp)
  R5=data.frame(Method="MI_NU", Estimate=(ODR),Lower.Limit=round(R51$conf.int[1],4), Upper.Limit=round(R51$conf.int[2],4))

  R61=Cornfieldexact.CI(ai, rt1, ci, rt2, conf = 1-alp)
  R6=data.frame(Method="Cornfield.Exact", Estimate=(ODR),Lower.Limit=round(R61[1],4), Upper.Limit=round(R61[2],4))

  R71=Cornfieldmidp.CI(ai, rt1, ci, rt2, conf = 1-alp)
  R7=data.frame(Method="Cornfield.midp", Estimate=(ODR),Lower.Limit=round(R71[1],4), Upper.Limit=round(R71[2],4))

  R81=BPexact.CI(ai, rt1, ci, rt2, conf = 1-alp)
  R8=data.frame(Method="BaPiEx", Estimate=(ODR),Lower.Limit=round(R81[1],4), Upper.Limit=round(R81[2],4))

  R91=BPmidp.CI(ai, rt1, ci, rt2, conf = 1-alp)
  R9=data.frame(Method="BaPiMP", Estimate=(ODR),Lower.Limit=round(R91[1],4), Upper.Limit=round(R91[2],4))

  R101=Invsinh.CI(ai, rt1, ci, rt2,e1,e2,conf = 1-alp)#Adjusted
  R10=data.frame(Method="ISH.Adjusted", Estimate=(ODR),Lower.Limit=round(R101[1],4), Upper.Limit=round(R101[2],4))

  #Inverse Sinh
  z=qnorm(alp/2,lower.tail = FALSE)
  term=2*asinh((z/2)*sqrt((1/ai)+(1/bi)+(1/ci)+(1/di)))
  LL_SiH= log(ODR)-term
  UL_SiH=log(ODR)+term
  R11=data.frame(Method="ISH", Estimate=(ODR),Lower.Limit=round(exp(LL_SiH),4), Upper.Limit=round(exp(UL_SiH),4))

  res_OR=rbind(R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11)

  return(res_OR)

} # End of the function

