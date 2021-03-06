% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/110.Sub-table-analysis.R
\name{generate.heatmap.matrix}
\alias{generate.heatmap.matrix}
\title{Given an input matrix, we can split it into smaller sub-matrix (min 2x2) and then find
 the Chi-squared test for each sub-matrix. The smaller matrix can "support"  or "oppose"
 (have a different conclusion at 95% confidence interaval) compared with the overall Chi-squared
 test value of the full input matrix. We count the number of times each cell supports or
 opposes the overall Chi-squared test. We also generate the possible list of sub-matrix.}
\usage{
generate.heatmap.matrix(mat, details = FALSE)
}
\arguments{
\item{mat}{- matrix for which the sub-matrix is to be generated}

\item{details}{- If this is set to TRUE, the the return value includes the full list of
sub-matrix}
}
\value{
A list  with
 \item{Hot.df}{  Dataframe with the difference between the supporting matrix and the opposing matrix}
 \item{Suport.df }{ Dataframe of the cell counts for support of table level Chi-sqaured }
 \item{Oppose.df }{ Dataframe of the cell counts for opposing of table level Chi-sqaured }
 \item{sub-matrix }{ list of sub-matrix - this is returned only if the details flag is set to TRUE }
}
\description{
Given an input matrix, we can split it into smaller sub-matrix (min 2x2) and then find
 the Chi-squared test for each sub-matrix. The smaller matrix can "support"  or "oppose"
 (have a different conclusion at 95% confidence interaval) compared with the overall Chi-squared
 test value of the full input matrix. We count the number of times each cell supports or
 opposes the overall Chi-squared test. We also generate the possible list of sub-matrix.
}
\details{
This can be used as an outlier detection method as well as observing the individual
cells within an IxJ table
}
\examples{
Drills=c(2,	10,	4,	2 ) # Example data from [reference 1]
Pots= c(3,	8,	4,	6)
Grinding.Stones=c( 13, 5, 3, 9)
Point.Fragments=c(20, 36, 19, 20)
mat=rbind(Drills,Pots,Grinding.Stones,Point.Fragments)
generate.heatmap.matrix(mat,details=FALSE)
}
\references{
[1] Mosteller F, Parunak A (2006)
Identifying extreme cells in a sizable contingency table: Probabilistic and exploratory approaches.
In: Hoaglin DC, Mosteller F, Tukey JW (eds) Exploring Data Tables, Trends, and Shapes,
John Wiley & Sons, pp 189-224
}
\seealso{
Other IxJ Inference methods: \code{\link{Local.Odds.Ratio}},
  \code{\link{Nominal.Odds.Ratio}},
  \code{\link{Plot.heatmap}}, \code{\link{Reversal.point}},
  \code{\link{Subtable.Odds.Ratio}}
}
