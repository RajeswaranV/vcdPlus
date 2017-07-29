#'  Given an input matrix, we can split it into smaller sub-matrix (min 2x2) and then find
#'  the Chi-squared test for each sub-matrix. The smaller matrix can "support"  or "oppose"
#'  (have a different conclusion at 95% confidence interaval) compared with the overall Chi-squared
#'  test value of the full input matrix. We count the number of times each cell supports or
#'  opposes the overall Chi-squared test. We also generate the possible list of sub-matrix.
#' @param mat - matrix for which the sub-matrix is to be generated
#' @details  This can be used as an outlier detection method as well as observing the individual
#' cells within an IxJ table
#' @return A dataframe with
#'  \item{Dimention.row}{Number of rows}
#'  \item{Dimention.col}{Number of columns}
#'  \item{Chi-sq-test.value}{Chi squared test value}
#'  \item{Significance at .95}{Test for significance at .95 confidence level}
#' @family I x J table summary functions
#' @examples
#'## Most influential school of Psychiatric thought and ascribed origin of schizophrenia- Agresti 1992
#' Eclectic=c(90,	12,	78) # Example from [reference 4]
#' Medical=c(13,	1,	6)
#' Psychoanalytic=c(19,	13,	50)
#' mat=rbind(Eclectic,Medical,Psychoanalytic)
#' Reversal.Association.Finder(mat)
#' @references
#' [1] J.Berkson
#' Some difficulties of interpretation encountered in the application of the chi-square test
#' Journal of the American Statistical Association. 33,  1938, 526-536.
#' [2] H.W,Norton
#' Calculation of chi-square for complex contingency tables
#' Journal of the American Statistical Association. 40, 1945, 251-258.
#' [3] C.R. Blyth
#' On Simpsons paradox and the sure thing principle.
#' Journal of the American Statistical Association. 67, 1972, 364-366.
#' [4] A.Agresti
#' Categorical Data Analysis
#' (New York: Wiley & Sons 1990) pp 51-54.
#' @export
Reversal.Association.Finder <-function(mat)
{
  if (missing(mat)) stop("mat is missing")
  if ((class(mat) != "matrix"))  stop("mat has to be a matrix with minimum 2x2")
  if ((dim(mat)[1] < 2) || (dim(mat)[2] <2 )) stop("Matrix has to be minimum of 2x2")

  out.df = NULL

  Total.columns=ncol(mat)
  Total.rows=nrow(mat)
  Overall.Chi.Sq=(chisq.test(mat)$p.value < 0.05)
  Big.counter.mat=mat-mat
  Positive.mat=Negetive.mat=Big.counter.mat

  # Take cols first and then find the row elements for each col combination
  # Column.combos=Total.columns-2
  for(Column.iterator in 2:Total.columns){
    #  print(Column.iterator)
    Current.Column.Combo=combn(Total.columns,Column.iterator)

    # Now loop over the combinations of columns
    # This loop below will generate the data for the various column combinations
    # Using this we go over each row combination
    for(Loop.col in 1:ncol(Current.Column.Combo)){

      Input.col.mat = mat[,c(Current.Column.Combo[,Loop.col])]

      # Now generating the row combinations for the given columns
      for(Row.iterator in 2:Total.rows){

        Current.Row.Combo=combn(Total.rows,Row.iterator)
        for(Row.lowest.selection in 1:ncol(Current.Row.Combo)){
          Useful.mat = Input.col.mat[c(Current.Row.Combo[,Row.lowest.selection]),]

          out.full = as.data.frame(t(c(dim(Useful.mat),
                                       chisq.test(Useful.mat)$p.value,
                                       (chisq.test(Useful.mat)$p.value < 0.05))))
          out.df = rbind(out.df, out.full)

        } # End of for loop for lowest Row selection

      } # End of for loop in Row.iterator for the current col combination

    } # End of for loop in Loop.col for the current col combination

  } # End of for-loop Column.iterator

  names(out.df)=c("Dimension.row","Dimention.col", "Chi-sq-test.value", "Significance.at.95")
  return(out.df)
} # End of function
