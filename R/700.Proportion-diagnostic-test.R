#'  This function  provides proportion estimates based on the
#'  given 2 x 2 table. This is an array of estimates other than association measures such as odds
#'  ratio, relative risk and risk difference. Main aim is to have efficiency of diagonstic test for
#'   a disease. Data is available in a 2 x 2 table which classifies the cases of
#'   correct / incorrect diagonsis of a test.
#' @param ai - Numeric value representing True Positives
#' @param bi - Numeric value representing False Positives
#' @param ci - Numeric value representing False Negatives
#' @param di - Numeric value representing True Negatives
#' @param p  - Numeric value representing prevelance
#' @details  The function takes True Positives, False Positives , False Negatives , True Negatives  and
#' prevalance as input and then calculates the efficiency of diagonstic test.
#' @return A dataframe with
#'  \item{Sensitivity}{The proportion of positives that are correctly identified as such }
#'  \item{Specificity}{The proportion of negatives that are correctly identified as such }
#'  \item{JIndex}{Sensitivity + Specificty -1}
#'  \item{Posivitve.Predictive.Value}{Proportions of positive results that are true positive}
#'  \item{Negetive.Predictive.Value}{Proportions of negative results that are  true negative}
#'  \item{Positive.likelihood.ratio}{Likelihood ratio for positive results}
#'  \item{Negative.likelihood.ratio}{Likelihood ratio for negetive results}
#'  \item{PPV.adjusted.for.Prevalance}{Posivitve predictive value adjusted for prevelance}
#'  \item{NPV.adjusted.for.Prevalance}{Negetive predictive value adjusted for prevelance}
#' @family Test methods
#' @examples
#'## Example taken from Breast cancer data [reference 1]
#' ai=31; bi=12; ci=3; di=32; p=0.07
#' Proportion.diagnostic.test(ai,bi,ci,di,p)
#'########Interval estimation
#' @references
#' [1] Mercaldo, N. D., Lau, K. F., & Zhou, X. H. (2007).
#' Confidence intervals for predictive values with an emphasis to case-control studies.
#' Statistics in medicine, 26(10), 2170-2183.
#' [2] Upton, G. J. (2016).
#' Categorical Data Analysis by Example.
#' John Wiley & Sons.
#' @export
Proportion.diagnostic.test<-function(ai,bi,ci,di,p)
{
  if ((class(ai) != "numeric") & (class(ai) != "integer") ||
      (class(bi) != "numeric")  & (class(bi) != "integer")   ||
      (class(ci) != "numeric")  & (class(ci) != "integer")  ||
      (class(di) != "numeric")  & (class(di) != "integer"))   stop("Inputs ai,bi,ci and di have to be numeric")
  if ((class(p) != "numeric") & (class(p) != "integer") || p>1 || p<0) stop("prevelance (p) has to be between 0 and 1")

  m = ai + bi
  n = ci + di
  r = ai + ci
  s =  bi + di
  SENST=ai/r
  SPECI=di/s
  JINDX=SENST+SPECI-1
  PPV=ai/m
  NPV=di/n
  PLR=SENST/(1-SPECI)
  NLR=(1-SENST)/SPECI
  PPV_p = SENST*p/(SENST*p+(1-SPECI)*(1-p))
  NPV_p = SPECI*(1-p)/((1-SENST)*p+SPECI*(1-p))
  #ans=round(rbind(SENST,SPECI,JINDX,PPV,NPV,PLR,NLR,PPV_p,NPV_p),4)
  ans.df=data.frame(Sensitivity=round(SENST,4),
                    Specificity=round(SPECI,4),
                    JIndex=round(JINDX,4),
                    Posivitve.Predictive.Value= round(PPV,4),
                    Negetive.Predictive.Value=  round(NPV,4),
                    Positive.likelihood.ratio=  round(PLR,4),
                    Negative.likelihood.ratio=  round(NLR,4),
                    PPV.adjusted.for.Prevalance=round(PPV_p,4),
                    NPV.adjusted.for.Prevalance=round(NPV_p,4))

  return(ans.df)

} # End of the function

