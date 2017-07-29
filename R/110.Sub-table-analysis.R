#'  Given an input matrix, we can split it into smaller sub-matrix (min 2x2) and then find
#'  the Chi-squared test for each sub-matrix. The smaller matrix can "support"  or "oppose"
#'  (have a different conclusion at 95% confidence interaval) compared with the overall Chi-squared
#'  test value of the full input matrix. We count the number of times each cell supports or
#'  opposes the overall Chi-squared test. We also generate the possible list of sub-matrix.
#' @param mat - matrix for which the sub-matrix is to be generated
#' @param details - If this is set to TRUE, the the return value includes the full list of
#' sub-matrix
#' @details  This can be used as an outlier detection method as well as observing the individual
#' cells within an IxJ table
#' @return A list  with
#'  \item{Hot.df}{  Dataframe with the difference between the supporting matrix and the opposing matrix}
#'  \item{Suport.df }{ Dataframe of the cell counts for support of table level Chi-sqaured }
#'  \item{Oppose.df }{ Dataframe of the cell counts for opposing of table level Chi-sqaured }
#'  \item{sub-matrix }{ list of sub-matrix - this is returned only if the details flag is set to TRUE }
#' @family IxJ Inference methods
#' @examples
#' Drills=c(2,	10,	4,	2 ) # Example data from [reference 1]
#' Pots= c(3,	8,	4,	6)
#' Grinding.Stones=c( 13, 5, 3, 9)
#' Point.Fragments=c(20, 36, 19, 20)
#' mat=rbind(Drills,Pots,Grinding.Stones,Point.Fragments)
#' generate.heatmap.matrix(mat,details=FALSE)
#' @references
#' [1] Mosteller F, Parunak A (2006)
#' Identifying extreme cells in a sizable contingency table: Probabilistic and exploratory approaches.
#' In: Hoaglin DC, Mosteller F, Tukey JW (eds) Exploring Data Tables, Trends, and Shapes,
#' John Wiley & Sons, pp 189-224
#' @export
generate.heatmap.matrix<-function(mat,details=FALSE)
{
  if (missing(mat)) stop("'mat' is missing")
  if ((class(mat) != "matrix"))  stop("'mat' has to be a matrix with minimum 2x2")
  if ((dim(mat)[1] < 2) || (dim(mat)[2] <2 )) stop("Matrix has to be minimum of 2x2")
  if ((class(details) != "logical")) stop("details has to be a logical vector")

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
    #print(Input.col.mat)

    # Now generating the row combinations for the given columns
    for(Row.iterator in 2:Total.rows){
      #  print(Row.iterator)
      Current.Row.Combo=combn(Total.rows,Row.iterator)
        for(Row.lowest.selection in 1:ncol(Current.Row.Combo)){
        Useful.mat = Input.col.mat[c(Current.Row.Combo[,Row.lowest.selection]),]
      #  print(Useful.mat)
      #  Big.counter.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])]=
       #   Big.counter.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])]+1
        if(chisq.test(Useful.mat)$p.value < 0.05){
          Positive.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])] = Positive.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])]+1
          } else { Negetive.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])] = Negetive.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])]+1
            } # End of chi-squared test and matrix update

        out.full = as.data.frame(t(c(dim(Useful.mat),chisq.test(Useful.mat)$p.value,(chisq.test(Useful.mat)$p.value < 0.05))))
        out.df = rbind(out.df, out.full)
        #  write.table(out.full,"D:/Research/CDA/Testing/Full.out.csv",append=TRUE,row.names=FALSE, col.names = FALSE, sep=",")
        } # End of for loop for lowest Row selection

    } # End of for loop in Row.iterator for the current col combination

    } # End of for loop in Loop.col for the current col combination

} #End of for-loop Column.iterator


# write.table(out.df,"D:/Research/CDA/Testing/Full.out.csv",append=TRUE,row.names=FALSE, col.names = FALSE, sep=",")


# We need to correct the chi-sq value by 1 to ensure we dont use the full matrix count
if(chisq.test(mat)$p.value < 0.05){
  Positive.mat = Positive.mat-1
  Support.mat = Positive.mat
  Oppose.mat = Negetive.mat
} else { Negetive.mat= Negetive.mat-1
Support.mat = Negetive.mat
Oppose.mat = Positive.mat
} # End of chi-squared test correction to ensure that we dont take full matrix value

Support.df=data.frame(Support=Support.mat)
Oppose.df=data.frame(Oppose=Oppose.mat)
Hot.mat= Support.mat-Oppose.mat
Hot.df = data.frame(Heatmap=Hot.mat)

names(out.df)=c("# of rows", "# of columns", "Chisq", "Significant at 95%")

if(details){st.list =list(Hot.df,Support.df,Oppose.df, out.df)
}else {st.list=list(Hot.df,Support.df,Oppose.df)}

return(st.list)
} # End of function
