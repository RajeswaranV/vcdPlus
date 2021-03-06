% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/500.SummaryMeasuresCategoricalData20170518.R
\name{SumMedianSD}
\alias{SumMedianSD}
\title{Calculates the median, SD, upper limit, lower limit and cell values within these limits}
\usage{
SumMedianSD(inpmat)
}
\arguments{
\item{inpmat}{- The input matrix}
}
\value{
A dataframe that provides the median, SD, upper limit, lower limit and cell
values within these limits
 \item{Feature}{Feature description -
 "Median of values in the input" showing the median of the input cells
 "SD of values in the input" showing the standard deviation (SD) of the input cells
 "Lower Limit of values in the input" showing the value of (median - 3*SD)
 "Upper Limit of values in the input" showing the maximum of (median + 3*SD)
 "Percentage of cells in the input whose values lie between the upper and lower limits" showing the percentage of the cells which fall within this range"}
 \item{Value}{ Calculated value of features}
}
\description{
Calculates the median, SD, upper limit, lower limit and cell values within these limits
}
\details{
The function calculates the Median, SD, Lower limit (defined as median - 3SD), Upper limit (defined as median + 3*SD) and the percentage of the cells which fall within this range
}
\examples{
x=c(7,61,55,129,489,570,475,431,293,154,38,12)
inpmat = cbind(x[1:4],x[5:8],x[9:12])
SumMedianSD(inpmat)
}
\seealso{
Other The summary of measures for categorical data: \code{\link{SumCount}},
  \code{\link{SumMinMaxRange}}, \code{\link{SumZeroData}},
  \code{\link{SummaryData}}
}
\author{
Subbiah and Balakrishna  S Kesavan
}
