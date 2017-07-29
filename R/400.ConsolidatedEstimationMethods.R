#' The simultaneous confidence interval for multinomial proportions based on the method proposed in Quesenberry and Hurst (1964)
#' @param inpmat - The input matrix
#' @param alpha - Alpha value (significance level required)
#' @return A list of dataframes
#'  \item{QH.Volume}{ QH Volume}
#'  \item{QH.UpLim}{  Dataframe of QH Upper Limits}
#'  \item{QH.LowLim}{  Dataframe of QH Lower Limits}
#'  \item{QH.Length}{  Dataframe of QH Lengths}
#' @family Confidence Interval for Multinomial Proportion
#' @examples
#' x = c(56,72,73,59,62,87,68,99,98)
#' inpmat = cbind(x[1:3],x[4:6],x[7:9])
#' alpha=0.05
#' ci.QH(inpmat,alpha)
#' @references
#' [1] Quesensberry, C.P. and Hurst, D.C. (1964).
#' Large Sample Simultaneous Confidence Intervals for Multinational Proportions.
#' Technometrics, 6: 191-195.
#' @export
ci.QH<-function(inpmat, alpha)
{
  if (missing(inpmat)) stop("Categorical data table is missing")
  if (missing(alpha)) stop("'alpha' is missing")
  if ((class(inpmat) != "matrix") || sum(inpmat < 0) > 0) stop("Check and correct Categorical data table")
  if (alpha>1 || alpha<0 || length(alpha) >1) stop("'alpha' has to be between 0 and 1")

  k = length(inpmat)
  s = sum(inpmat)
  chi = qchisq(1-alpha, df=k-1)
  pi = inpmat/s

  QH.UL = round((chi + 2*inpmat + sqrt(chi*chi + 4*inpmat*chi*(1 - pi)))/(2*(chi+s)),3)
  QH.LL = round((chi + 2*inpmat - sqrt(chi*chi + 4*inpmat*chi*(1 - pi)))/(2*(chi+s)),3)


  QH.WI = QH.Length=QH.UL - QH.LL 	#Length of the interval
  QH.VL = prod(QH.WI)


  QH.UL.df = data.frame(QH.UpLim=QH.UL)
  QH.LL.df = data.frame(QH.LowLim=QH.LL)
  QH.WI.df = data.frame(QH.Length=QH.WI)
  QH.VL.df = data.frame(QH.Volume=QH.VL)

  newList <- list(QH.VL.df,QH.UL.df,QH.LL.df,QH.WI.df)
  return(newList)
}
##############################################################################
#' The simultaneous confidence interval for multinomial proportions based on the method proposed in Goodman (1965)
#' @param inpmat - The input matrix
#' @param alpha - Alpha value (significance level required)
#' @return A list of dataframes
#'  \item{GM.Volume}{ GM Volume}
#'  \item{GM.UpLim}{  Dataframe of GM Upper Limits}
#'  \item{GM.LowLim}{  Dataframe of GM Lower Limits}
#'  \item{GM.Length}{  Dataframe of GM Lengths}
#' @family Confidence Interval for Multinomial Proportion
#' @examples
#' x = c(56,72,73,59,62,87,68,99,98)
#' inpmat = cbind(x[1:3],x[4:6],x[7:9])
#' alpha=0.05
#' ci.GM(inpmat,alpha)
#' @references
#' [1] Goodman, L.A. (1965).
#' On Simultaneous Confidence Intervals for Multinomial Proportions.
#' Technometrics 7: 247-254.
#' @export
ci.GM<-function(inpmat, alpha)
{
  if (missing(inpmat)) stop("Categorical data table is missing")
  if (missing(alpha)) stop("'alpha' is missing")
  if ((class(inpmat) != "matrix") || sum(inpmat < 0) > 0) stop("Check and correct Categorical data table")
  if (alpha>1 || alpha<0 || length(alpha) >1) stop("'alpha' has to be between 0 and 1")

  k = length(inpmat)
  s = sum(inpmat)
  chi = qchisq(1-(alpha/k), df=1)
  pi = inpmat/s
  GM.UL = round((chi + 2*inpmat + sqrt(chi*chi + 4*inpmat*chi*(1 - pi)))/(2*(chi+s)),3)
  GM.LL = round((chi + 2*inpmat - sqrt(chi*chi + 4*inpmat*chi*(1 - pi)))/(2*(chi+s)),3)
  GM.WI = GM.UL - GM.LL 	#Length of the interval
  GM.VL = prod(GM.WI)

  GM.UL.df = data.frame(GM.UpLim=GM.UL)
  GM.LL.df = data.frame(GM.LowLim=GM.LL)
  GM.WI.df = data.frame(GM.Length=GM.WI)
  GM.VL.df = data.frame(GM.Volume=GM.VL)

  newList <- list(GM.VL.df,GM.UL.df,GM.LL.df,GM.WI.df)
  return(newList)
}
##############################################################################
#' The simple Wald type interval for multinomial proportions which is symmetrical about the sample
#' proportions. In this method no continuity corrections are made to avoid zero width intervals when
#' the sample proportions are at extreme.
#' @param inpmat - The input matrix
#' @param alpha - Alpha value (significance level required)
#' @return A list of dataframes
#'  \item{Wald.Volume}{ Wald Volume}
#'  \item{Wald.UpLim}{  Dataframe of Wald Upper Limits}
#'  \item{Wald.LowLim}{  Dataframe of Wald Lower Limits}
#'  \item{Wald.Length}{  Dataframe of Wald Lengths}
#' @family Confidence Interval for Multinomial Proportion
#' @examples
#' x = c(56,72,73,59,62,87,68,99,98)
#' inpmat = cbind(x[1:3],x[4:6],x[7:9])
#' alpha=0.05
#' ci.Wald(inpmat,alpha)
#' @references
#' [1] Wald,
#' A Tests of statistical hypotheses concerning several parameters when the number of observations is large,
#' Trans. Am. Math. Soc. 54 (1943) 426-482.
#' @export
ci.Wald<-function(inpmat, alpha)
{
  if (missing(inpmat)) stop("Categorical data table is missing")
  if (missing(alpha)) stop("'alpha' is missing")
  if ((class(inpmat) != "matrix") || sum(inpmat < 0) > 0) stop("Check and correct Categorical data table")
  if (alpha>1 || alpha<0 || length(alpha) >1) stop("'alpha' has to be between 0 and 1")

  k = length(inpmat)
  s = sum(inpmat)
  chi = qchisq(1-alpha, df=1)
  pi = inpmat/s
  Wald.LL = pi - (sqrt(chi*(pi)*(1-pi)/s))
  Wald.UL = pi + (sqrt(chi*(pi)*(1-pi)/s))
  Wald.WI = Wald.UL-Wald.LL
  Wald.VL = prod(Wald.WI)

  Wald.UL.df = data.frame(Wald.UpLim=Wald.UL)
  Wald.LL.df = data.frame(Wald.LowLim=Wald.LL)
  Wald.WI.df = data.frame(Wald.Length=Wald.WI)
  Wald.VL.df = data.frame(Wald.Volume=Wald.VL)

  newList <- list(Wald.VL.df,Wald.UL.df,Wald.LL.df,Wald.WI.df)
  return(newList)
}
###############################################################################
#' The simple Wald type interval with continuity corrections for multinomial proportions which
#' is symmetrical about the sample proportions.
#' @param inpmat - The input matrix
#' @param alpha - Alpha value (significance level required)
#' @return A list of dataframes
#'  \item{WaldCC.Volume}{ Wald Volume}
#'  \item{WaldCC.UpLim}{  Dataframe of Wald Upper Limits}
#'  \item{WaldCC.LowLim}{  Dataframe of Wald Lower Limits}
#'  \item{WaldCC.Length}{  Dataframe of Wald Lengths}
#' @family Confidence Interval for Multinomial Proportion
#' @examples
#' x = c(56,72,73,59,62,87,68,99,98)
#' inpmat = cbind(x[1:3],x[4:6],x[7:9])
#' alpha=0.05
#' ci.WaldCC(inpmat,alpha)
#' @export
ci.WaldCC<-function(inpmat, alpha)
{
  if (missing(inpmat)) stop("Categorical data table is missing")
  if (missing(alpha)) stop("'alpha' is missing")
  if ((class(inpmat) != "matrix") || sum(inpmat < 0) > 0) stop("Check and correct Categorical data table")
  if (alpha>1 || alpha<0 || length(alpha) >1) stop("'alpha' has to be between 0 and 1")

  k = length(inpmat)
  s = sum(inpmat)
  chi = qchisq(1-alpha, df=1)
  pi = inpmat/s
  WALDCC.LL = pi - (sqrt(chi*(pi)*(1-pi)/s))-(1/(2*s))
  WALDCC.UL = pi + (sqrt(chi*(pi)*(1-pi)/s))+(1/(2*s))
  WALDCC.WI=WALDCC.UL-WALDCC.LL
  WALDCCTable <- round(cbind(WALDCC.LL, WALDCC.UL,WALDCC.WI),3)
  WALDCC.VL = prod(WALDCC.WI)

  WALDCC.UL.df = data.frame(WaldCC.UpLim=WALDCC.UL)
  WALDCC.LL.df = data.frame(WaldCC.LowLim=WALDCC.LL)
  WALDCC.WI.df = data.frame(WaldCC.Length=WALDCC.WI)
  WALDCC.VL.df = data.frame(WaldCC.Volume=WALDCC.VL)

  newList <- list(WALDCC.VL.df,WALDCC.UL.df,WALDCC.LL.df,WALDCC.WI.df)
  return(newList)
}
###############################################################
#' The simultaneous confidence interval for multinomial proportions based on the
#' method proposed in Fitzpatrick and Scott (1987)
#' @param inpmat - The input matrix
#' @param alpha - Alpha value (significance level required)
#' @return A list of dataframes
#'  \item{WaldCC.Volume}{ Wald Volume}
#'  \item{WaldCC.UpLim}{  Dataframe of Wald Upper Limits}
#'  \item{WaldCC.LowLim}{  Dataframe of Wald Lower Limits}
#'  \item{WaldCC.Length}{  Dataframe of Wald Lengths}
#' @family Confidence Interval for Multinomial Proportion
#' @examples
#' x = c(56,72,73,59,62,87,68,99,98)
#' inpmat = cbind(x[1:3],x[4:6],x[7:9])
#' alpha=0.05
#' ci.FS(inpmat,alpha)
#' @references
#' [1] Fitzpatrick, S. and Scott, A. (1987).
#' Quick simultaneous confidence interval for multinomial proportions.
#' Journal of American Statistical Association 82(399): 875-878.
#' @export
ci.FS<-function(inpmat, alpha)
{
  if (missing(inpmat)) stop("Categorical data table is missing")
  if (missing(alpha)) stop("'alpha' is missing")
  if ((class(inpmat) != "matrix") || sum(inpmat < 0) > 0) stop("Check and correct Categorical data table")
  if (alpha>1 || alpha<0 || length(alpha) >1) stop("'alpha' has to be between 0 and 1")

  k = length(inpmat)
  s = sum(inpmat)
  zval = abs(qnorm(1-(alpha/2)))
  pi = inpmat/s
  FS.LL = pi - (zval/(2*sqrt(s)))
  FS.UL = pi + (zval/(2*sqrt(s)))
  FS.WI = FS.UL-FS.LL
  FS.VL = prod(FS.WI)

  FS.UL.df = data.frame(FS.UpLim=FS.UL)
  FS.LL.df = data.frame(FS.LowLim=FS.LL)
  FS.WI.df = data.frame(FS.Length=FS.WI)
  FS.VL.df = data.frame(FS.Volume=FS.VL)

  newList <- list(FS.VL.df,FS.UL.df,FS.LL.df,FS.WI.df)
  return(newList)
}
###################################################################
#' The simultaneous confidence interval for multinomial proportions based on the method
#' proposed in Wilson (1927)
#' @param inpmat - The input matrix
#' @param alpha - Alpha value (significance level required)
#' @return A list of dataframes
#'  \item{WaldCC.Volume}{ Wald Volume}
#'  \item{WaldCC.UpLim}{  Dataframe of Wald Upper Limits}
#'  \item{WaldCC.LowLim}{  Dataframe of Wald Lower Limits}
#'  \item{WaldCC.Length}{  Dataframe of Wald Lengths}
#' @family Confidence Interval for Multinomial Proportion
#' @examples
#' x = c(56,72,73,59,62,87,68,99,98)
#' inpmat = cbind(x[1:3],x[4:6],x[7:9])
#' alpha=0.05
#' ci.WS(inpmat,alpha)
#' @references
#' [1] E.B. Wilson,
#' Probable inference, the law of succession and statistical inference,
#' J.Am. Stat. Assoc.22 (1927) 209-212.
#' @export
ci.WS<-function(inpmat, alpha)
{
  if (missing(inpmat)) stop("Categorical data table is missing")
  if (missing(alpha)) stop("'alpha' is missing")
  if ((class(inpmat) != "matrix") || sum(inpmat < 0) > 0) stop("Check and correct Categorical data table")
  if (alpha>1 || alpha<0 || length(alpha) >1) stop("'alpha' has to be between 0 and 1")

  k = length(inpmat)
  s = sum(inpmat)
  chi = qchisq(1-alpha, df=1)
  pi = inpmat/s
  WS.UL = (chi + 2*inpmat + sqrt(chi*chi + 4*inpmat*chi*(1 - pi)))/(2*(chi+s))
  WS.LL = (chi + 2*inpmat - sqrt(chi*chi + 4*inpmat*chi*(1 - pi)))/(2*(chi+s))
  WS.WI = WS.UL-WS.LL
  WS.VL = prod(WS.WI)
  WS.VL = prod(WS.WI)

  WS.UL.df = data.frame(WS.UpLim=WS.UL)
  WS.LL.df = data.frame(WS.LowLim=WS.LL)
  WS.WI.df = data.frame(WS.Length=WS.WI)
  WS.VL.df = data.frame(WS.Volume=WS.VL)

  newList <- list(WS.VL.df,WS.UL.df,WS.LL.df,WS.WI.df)
  return(newList)
}
