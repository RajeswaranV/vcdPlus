% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/500.SummaryMeasuresCategoricalData20170518.R
\name{SumCount}
\alias{SumCount}
\title{The summary of measures for categorical data}
\usage{
SumCount(inpmat, min = NULL, max = NULL)
}
\arguments{
\item{inpmat}{- The input matrix}

\item{min}{- The minimum value for the range}

\item{max}{- The maximum value for the range}
}
\value{
A dataframe that provides a count of cells in a given matrix with a value
of zero or value within the range specified as min and max
 \item{Feature}{Feature description -
 "# of cells with zero value" showing the count of cells with all zeros
 "# of cells with value in specified range" showing the count of cells which have value within the specified range
 "Percentage of cells with zero value" showing the percentage of cells with all zeros
 "Percentage of cells with value in specified range" showing the percentage of cells which have value within the specified range}
 \item{Value}{ Calculated value of features}
}
\description{
The summary of measures for categorical data
}
\details{
The function calculates the number and percentage of cells which are all zeros within the given range
}
\examples{
x=c(7,61,55,129,489,570,475,431,293,154,38,12)
inpmat = cbind(x[1:4],x[5:8],x[9:12])
min = 0
max = 150
SumCount(inpmat, min, max)
}
\seealso{
Other The summary of measures for categorical data: \code{\link{SumMedianSD}},
  \code{\link{SumMinMaxRange}}, \code{\link{SumZeroData}},
  \code{\link{SummaryData}}
}
\author{
Subbiah and Balakrishna  S Kesavan
}
