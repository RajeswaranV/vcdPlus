% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/400.ConsolidatedEstimationMethods.R
\name{ci.Wald}
\alias{ci.Wald}
\title{The simple Wald type interval for multinomial proportions which is symmetrical about the sample
proportions. In this method no continuity corrections are made to avoid zero width intervals when
the sample proportions are at extreme.}
\usage{
ci.Wald(inpmat, alpha)
}
\arguments{
\item{inpmat}{- The input matrix}

\item{alpha}{- Alpha value (significance level required)}
}
\value{
A list of dataframes
 \item{Wald.Volume}{ Wald Volume}
 \item{Wald.UpLim}{  Dataframe of Wald Upper Limits}
 \item{Wald.LowLim}{  Dataframe of Wald Lower Limits}
 \item{Wald.Length}{  Dataframe of Wald Lengths}
}
\description{
The simple Wald type interval for multinomial proportions which is symmetrical about the sample
proportions. In this method no continuity corrections are made to avoid zero width intervals when
the sample proportions are at extreme.
}
\examples{
x = c(56,72,73,59,62,87,68,99,98)
inpmat = cbind(x[1:3],x[4:6],x[7:9])
alpha=0.05
ci.Wald(inpmat,alpha)
}
\references{
[1] Wald,
A Tests of statistical hypotheses concerning several parameters when the number of observations is large,
Trans. Am. Math. Soc. 54 (1943) 426-482.
}
\seealso{
Other Confidence Interval for Multinomial Proportion: \code{\link{ci.BMDU}},
  \code{\link{ci.FS}}, \code{\link{ci.GM}},
  \code{\link{ci.QH}}, \code{\link{ci.SG}},
  \code{\link{ci.WS}}, \code{\link{ci.WaldCC}}
}
\author{
Subbiah and Balakrishna S Kesavan
}
