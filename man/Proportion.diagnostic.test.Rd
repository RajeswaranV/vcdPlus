% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/700.Proportion-diagnostic-test.R
\name{Proportion.diagnostic.test}
\alias{Proportion.diagnostic.test}
\title{This function  provides proportion estimates based on the
 given 2 x 2 table. This is an array of estimates other than association measures such as odds
 ratio, relative risk and risk difference. Main aim is to have efficiency of diagonstic test for
  a disease. Data is available in a 2 x 2 table which classifies the cases of
  correct / incorrect diagonsis of a test.}
\usage{
Proportion.diagnostic.test(ai, bi, ci, di, p)
}
\arguments{
\item{ai}{- Numeric value representing True Positives}

\item{bi}{- Numeric value representing False Positives}

\item{ci}{- Numeric value representing False Negatives}

\item{di}{- Numeric value representing True Negatives}

\item{p}{- Numeric value representing prevelance}
}
\value{
A dataframe with
 \item{Sensitivity}{The proportion of positives that are correctly identified as such }
 \item{Specificity}{The proportion of negatives that are correctly identified as such }
 \item{JIndex}{Sensitivity + Specificty -1}
 \item{Posivitve.Predictive.Value}{Proportions of positive results that are true positive}
 \item{Negetive.Predictive.Value}{Proportions of negative results that are  true negative}
 \item{Positive.likelihood.ratio}{Likelihood ratio for positive results}
 \item{Negative.likelihood.ratio}{Likelihood ratio for negetive results}
 \item{PPV.adjusted.for.Prevalance}{Posivitve predictive value adjusted for prevelance}
 \item{NPV.adjusted.for.Prevalance}{Negetive predictive value adjusted for prevelance}
}
\description{
This function  provides proportion estimates based on the
 given 2 x 2 table. This is an array of estimates other than association measures such as odds
 ratio, relative risk and risk difference. Main aim is to have efficiency of diagonstic test for
  a disease. Data is available in a 2 x 2 table which classifies the cases of
  correct / incorrect diagonsis of a test.
}
\details{
The function takes True Positives, False Positives , False Negatives , True Negatives  and
prevalance as input and then calculates the efficiency of diagonstic test.
}
\examples{
## Example taken from Breast cancer data [reference 1]
ai=31; bi=12; ci=3; di=32; p=0.07
Proportion.diagnostic.test(ai,bi,ci,di,p)
########Interval estimation
}
\references{
[1] Mercaldo, N. D., Lau, K. F., & Zhou, X. H. (2007).
Confidence intervals for predictive values with an emphasis to case-control studies.
Statistics in medicine, 26(10), 2170-2183.
[2] Upton, G. J. (2016).
Categorical Data Analysis by Example.
John Wiley & Sons.
}
\seealso{
Other Test methods: \code{\link{Partition.table}},
  \code{\link{exact2x2tests.hypergeom}},
  \code{\link{exact2x2tests.regular}}
}
