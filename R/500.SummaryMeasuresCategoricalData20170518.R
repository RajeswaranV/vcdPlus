#' Calculates the count of rows and columns that contain all zeros, in a given matrix
#' @param inpmat - The input matrix
#' @param start.col - The minimum of the range of interest. Starting column number
#' @param end.col - The minimum of the range of interest. Ending column number
#' @return A dataframe that provides a count of rows and columns that contain all
#' zeros, in a given matrix
#'  \item{Feature}{Feature description -
#'  '# of all-zero rows' showing the count of rows with all zeros
#'  '# of all-zero columns' showing the count of columns with all zeros}
#'  \item{Value}{ Calculated value of counts}
#' @details  The function calculates the number of rows and columns which are all zeros, in the given range
#' @family The summary of measures for categorical data
#' @examples
#' x=c(7,61,55,129,489,570,475,431,293,154,38,12)
#' inpmat = cbind(x[1:4],x[5:8],x[9:12])
#' SumZeroData(inpmat)
#'###Example with some zeros
#' inpmat = matrix(c(0, 0, 0, 1, 0, 7), nrow=3, ncol=2)
#' SumZeroData(inpmat)
#' @export
SumZeroData <- function(inpmat, start.col = NULL, end.col = NULL){

  if (missing(inpmat)) stop("Categorical data table is missing")
  if ((class(inpmat) != "matrix") || sum(inpmat < 0) > 0) stop("Check and correct Categorical data table")
  if(is.null(start.col)) { start.col=1}
  if(is.null(end.col)) { end.col=ncol(inpmat)}
  if ( start.col < 0 ) stop("start.col value has to be  >= 0")
  if ( end.col < 0 ) stop("end.col value has to be  >= 0")
  if (end.col <= start.col) stop("end.col has to be greater than start.col")

  if ( start.col > ncol(inpmat) ||
       end.col > ncol(inpmat) )  stop("'start.col' and 'end.col' has to be between 1 and max number of columns")

  RowZero <- cbind.data.frame("# of all-zero rows =",sum(rowSums(inpmat[,start.col:end.col])==0))
  ColZero <- cbind.data.frame("# of all-zero columns =",sum(colSums(inpmat[,start.col:end.col])==0))

  colnames(RowZero) = c("Feature", "Value")
  colnames(ColZero) = c("Feature", "Value")

  SumZeroData <- rbind(RowZero, ColZero)
  return(SumZeroData)
}

###################################################################
#' The summary of measures for categorical data
#' @param inpmat - The input matrix
#' @param min - The minimum value for the range
#' @param max - The maximum value for the range
#' @return A dataframe that provides a count of cells in a given matrix with a value
#' of zero or value within the range specified as min and max
#'  \item{Feature}{Feature description -
#'  "# of cells with zero value" showing the count of cells with all zeros
#'  "# of cells with value in specified range" showing the count of cells which have value within the specified range
#'  "Percentage of cells with zero value" showing the percentage of cells with all zeros
#'  "Percentage of cells with value in specified range" showing the percentage of cells which have value within the specified range}
#'  \item{Value}{ Calculated value of features}
#' @details  The function calculates the number and percentage of cells which are all zeros within the given range
#' @family The summary of measures for categorical data
#' @examples
#' x=c(7,61,55,129,489,570,475,431,293,154,38,12)
#' inpmat = cbind(x[1:4],x[5:8],x[9:12])
#' min = 0
#' max = 150
#' SumCount(inpmat, min, max)
#' @export
SumCount <- function(inpmat, min = NULL, max = NULL){

  if (missing(inpmat)) stop("Categorical data table is missing")
  if ((class(inpmat) != "matrix") || sum(inpmat < 0) > 0) stop("Check and correct Categorical data table")
  if (is.null(min)) { min=0}
  if (is.null(max)) { max = Inf}
  if ((class(min) != "numeric") & (class(min) != "integer")) stop("Minimum Value has to be a number")
  if ((class(max) != "numeric") & (class(max) != "integer")) stop("Maximum Value has to be a number")
  if (min >= max) stop("Minimum value has to be less than maximum value")

  zeroValueCellVal <- sum(apply(inpmat, 2, function(c){sum(c==0)}))
  zeroValueCell <- cbind.data.frame("# of cells with zero value",zeroValueCellVal)
  specificRangeValueCellVal <- sum(apply(inpmat, 2, function(c){sum(c>=min & c<= max)}))
  specificRangeValueCell <- cbind.data.frame("# of cells with value in specified range",specificRangeValueCellVal)
  percZeroValueCell <- cbind.data.frame("Percentage of cells with zero value",round(zeroValueCellVal*100/(nrow(inpmat)*ncol(inpmat)),2))
  percSpecificRangeValueCell <- cbind.data.frame("Percentage of cells with value in specified range",round(specificRangeValueCellVal*100/(nrow(inpmat)*ncol(inpmat)),2))

  colnames(zeroValueCell) = c("Feature", "Value")
  colnames(specificRangeValueCell) = c("Feature", "Value")
  colnames(percZeroValueCell) = c("Feature", "Value")
  colnames(percSpecificRangeValueCell) = c("Feature", "Value")

  SumCount = rbind(zeroValueCell, specificRangeValueCell, percZeroValueCell, percSpecificRangeValueCell)

  return(SumCount)
}

###################################################################
#' Calculates the min, max and range of values in a given matrix
#' @param inpmat - The input matrix
#' @return A dataframe that provides the min, max and range of values in a given matrix
#'  \item{Feature}{ Feature description -
#'  " Max value in table" showing the maximum of the input cells
#'  "Min value in table" showing the minimum of the input cells
#'  "Range of values in table" showing the difference between the maximum and minimum"}
#'  \item{Value}{ Calculated value of features}
#' @details  The function calculates the min, max and range of values in a given matrix
#' @family The summary of measures for categorical data
#' @examples
#' x=c(7,61,55,129,489,570,475,431,293,154,38,12)
#' inpmat = cbind(x[1:4],x[5:8],x[9:12])
#' SumMinMaxRange(inpmat)
#' @export
SumMinMaxRange <- function(inpmat){

  if (missing(inpmat)) stop("Categorical data table is missing")
  if ((class(inpmat) != "matrix") || sum(inpmat < 0) > 0) stop("Check and correct Categorical data table")

  maxValueAcrossRowsVal <- max(apply(inpmat, 1, function(c){max(c)}))
  minValueAcrossRowsVal <- min(apply(inpmat, 1, function(c){min(c)}))
  maxValueAcrossRows <- cbind.data.frame("Max value in table  =",maxValueAcrossRowsVal)
  minValueAcrossRows <- cbind.data.frame("Min value in table  =",minValueAcrossRowsVal)
  rangeValue <- cbind.data.frame("Range of values in table  =",(maxValueAcrossRowsVal-minValueAcrossRowsVal))

  colnames(maxValueAcrossRows) = c("Feature", "Value")
  colnames(minValueAcrossRows) = c("Feature", "Value")
  colnames(rangeValue) = c("Feature", "Value")

  SumMinMaxRange <- rbind(maxValueAcrossRows, minValueAcrossRows, rangeValue)

  return(SumMinMaxRange)
}

###################################################################
#' Calculates the median, SD, upper limit, lower limit and cell values within these limits
#' @param inpmat - The input matrix
#' @return A dataframe that provides the median, SD, upper limit, lower limit and cell
#' values within these limits
#'  \item{Feature}{Feature description -
#'  "Median of values in the input" showing the median of the input cells
#'  "SD of values in the input" showing the standard deviation (SD) of the input cells
#'  "Lower Limit of values in the input" showing the value of (median - 3*SD)
#'  "Upper Limit of values in the input" showing the maximum of (median + 3*SD)
#'  "Percentage of cells in the input whose values lie between the upper and lower limits" showing the percentage of the cells which fall within this range"}
#'  \item{Value}{ Calculated value of features}
#' @details  The function calculates the Median, SD, Lower limit (defined as median - 3SD), Upper limit (defined as median + 3*SD) and the percentage of the cells which fall within this range
#' @family The summary of measures for categorical data
#' @examples
#' x=c(7,61,55,129,489,570,475,431,293,154,38,12)
#' inpmat = cbind(x[1:4],x[5:8],x[9:12])
#' SumMedianSD(inpmat)
#' @export
SumMedianSD <- function(inpmat){

  if (missing(inpmat)) stop("Categorical data table is missing")
  if ((class(inpmat) != "matrix") || sum(inpmat < 0) > 0) stop("Check and correct Categorical data table")

  medianValue <- median(as.vector(inpmat))
  sdValue <- round((sd(as.vector(inpmat))),2)
  lowerLimitValue <- round((medianValue-3*sdValue),2)
  upperLimitValue <- round((medianValue+3*sdValue),2)
  withinLimitCellsCount <- 100*sum(apply(inpmat, 2, function(c){sum(c>=lowerLimitValue & c<= upperLimitValue)}))/(nrow(inpmat)*ncol(inpmat))
  medianValue.df <- cbind.data.frame("Median of values in the input =", medianValue)
  sdValue.df <- cbind.data.frame("SD of values in the input  =", sdValue)
  lowerLimit <- cbind.data.frame("Lower Limit of values in the input  =", lowerLimitValue)
  upperLimit <- cbind.data.frame("Upper Limit of values in the input  =", upperLimitValue)
  withinLimitCells <- cbind.data.frame("Percentage of cells in the input whose values lie between the upper and lower limits =", withinLimitCellsCount)

  colnames(medianValue.df) = c("Feature", "Value")
  colnames(sdValue.df) = c("Feature", "Value")
  colnames(lowerLimit) = c("Feature", "Value")
  colnames(upperLimit) = c("Feature", "Value")
  colnames(withinLimitCells) = c("Feature", "Value")

  SumMedianSD <- rbind(medianValue.df, sdValue.df, lowerLimit, upperLimit, withinLimitCells)

  return(SumMedianSD)
}

###################################################################
#' The summary of measures for categorical data
#' @param inpmat - The input matrix
#' @param min - The minimum of the range of interest
#' @param max - The minimum of the range of interest
#' @return A list of dataframes
#'  \item{SumZeroData}{ Count of rows and columns that contain all zeros, in a given matrix}
#'  \item{SumCount}{ Count of cells in a given matrix with value of zero or value within the range specified as min and max}
#'  \item{SumMinMaxRange}{  Calculates the min, max and range of values in a given matrix}
#'  \item{SumMedianSD}{  Calculates median, SD, upper limit, lower limit and cell values within these limits}
#' @details  The function calculates the summary using
#' @family The summary of measures for categorical data
#' @examples
#' x=c(7,61,55,129,489,570,475,431,293,154,38,12)
#' inpmat = cbind(x[1:4],x[5:8],x[9:12])
#' min = 0
#' max = 150
#' SummaryData(inpmat, min, max)
#' SummaryData(inpmat)
#' @export
SummaryData <- function(inpmat, min=NULL, max=NULL){

  if (missing(inpmat)) stop("Categorical data table is missing")
  if ((class(inpmat) != "matrix") || sum(inpmat < 0) > 0) stop("Check and correct Categorical data table")
  if (is.null(min)) { min=0}
  if (is.null(max)) { max = Inf}
  if ((class(min) != "numeric") & (class(min) != "integer")) stop("Minimum Value has to be a number")
  if ((class(max) != "numeric") & (class(max) != "integer")) stop("Maximum Value has to be a number")
  if (min >= max) stop("Minimum value has to be less than maximum value")

  result <- list(SumZeroData=SumZeroData(inpmat),
                 SumCount=SumCount(inpmat, min, max),
                 SumMinMaxRange=SumMinMaxRange(inpmat),
                 SumMedianSD=SumMedianSD(inpmat))
  return(result)
}


