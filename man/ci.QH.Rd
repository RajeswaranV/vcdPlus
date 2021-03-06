% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/400.ConsolidatedEstimationMethods.R
\name{ci.QH}
\alias{ci.QH}
\title{The simultaneous confidence interval for multinomial proportions based on the method proposed in Quesenberry and Hurst (1964)}
\usage{
ci.QH(inpmat, alpha)
}
\arguments{
\item{inpmat}{- The input matrix}

\item{alpha}{- Alpha value (significance level required)}
}
\value{
A list of dataframes
 \item{QH.Volume}{ QH Volume}
 \item{QH.UpLim}{  Dataframe of QH Upper Limits}
 \item{QH.LowLim}{  Dataframe of QH Lower Limits}
 \item{QH.Length}{  Dataframe of QH Lengths}
}
\description{
The simultaneous confidence interval for multinomial proportions based on the method proposed in Quesenberry and Hurst (1964)
}
\examples{
x = c(56,72,73,59,62,87,68,99,98)
inpmat = cbind(x[1:3],x[4:6],x[7:9])
alpha=0.05
ci.QH(inpmat,alpha)
}
\references{
[1] Quesensberry, C.P. and Hurst, D.C. (1964).
Large Sample Simultaneous Confidence Intervals for Multinational Proportions.
Technometrics, 6: 191-195.
}
\seealso{
Other Confidence Interval for Multinomial Proportion: \code{\link{ci.BMDU}},
  \code{\link{ci.FS}}, \code{\link{ci.GM}},
  \code{\link{ci.SG}}, \code{\link{ci.WS}},
  \code{\link{ci.WaldCC}}, \code{\link{ci.Wald}}
}
\author{
Subbiah and Balakrishna S Kesavan
}
