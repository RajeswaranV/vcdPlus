#'  Given an input matrix, we can split it into smaller sub-matrix (min 2x2) and then find
#'  the Chi-squared test for each sub-matrix. The smaller matrix can "support"  or "oppose"
#'  (have a different conclusion at 95% confidence interaval) compared with the overall Chi-squared
#'  test value of the full input matrix. We count the number of times each cell supports or
#'  opposes the overall Chi-squared test. We also generate the possible list of sub-matrix.
#' @param mat - matrix for which the sub-matrix is to be generated
#' @details  This can be used as an outlier detection method as well as observing the individual
#' cells within an IxJ table
#' @return A of dataframes with
#'  \item{Hot.df}{  Dataframe with the difference between the supporting matrix and the opposing matrix}
#'  \item{Suport.df }{ Dataframe of the cell counts for support of table level Chi-sqaured }
#'  \item{Oppose.df }{ Dataframe of the cell counts for opposing of table level Chi-sqaured }
#'  \item{sub-matrix }{ list of sub-matrix - this is returned only if the details flag is set to TRUE }
#' @family IxJ Inference methods
#' @examples
#'##  Example data from [reference 1]
#' Drills=c(2,	10,	4,	2 )
#' Pots= c(3,	8,	4,	6)
#' Grinding.Stones=c( 13, 5, 3, 9)
#' Point.Fragments=c(20, 36, 19, 20)
#' mat=rbind(Drills,Pots,Grinding.Stones,Point.Fragments)
#' Local.Odds.Ratio(mat)
#' @references
#' [1] Mosteller F, Parunak A (2006)
#' Identifying extreme cells in a sizable contingency table: Probabilistic and exploratory approaches.
#' In: Hoaglin DC, Mosteller F, Tukey JW (eds) Exploring Data Tables, Trends, and Shapes,
#' John Wiley & Sons, pp 189-224
#' @export
Local.Odds.Ratio<-function(mat)
{
  if (missing(mat)) stop("'mat' is missing")
  if ((class(mat) != "matrix"))  stop("'mat' has to be a matrix with minimum 2x2")
  if ((dim(mat)[1] < 2) || (dim(mat)[2] <2 )) stop("Matrix has to be minimum of 2x2")

  out.df = NULL

  Total.columns=ncol(mat)
  Total.rows=nrow(mat)
  Big.counter.mat=mat-mat
  Index.mat=mat-mat

  # Take cols first and then find the row elements for each col combination
  # Column.combos=Total.columns-2
  Current.Column.Combo=combn(Total.columns,2)

  # Now loop over the combinations of columns
  # This loop below will generate the data for the various column combinations
  # Using this we go over each row combination
  for(Loop.col in 1:ncol(Current.Column.Combo)){

    Input.col.mat = mat[,c(Current.Column.Combo[,Loop.col])]
    #print(Input.col.mat)

    # Now generating the row combinations for the given columns
    Current.Row.Combo=combn(Total.rows,2)
    for(Row.lowest.selection in 1:ncol(Current.Row.Combo)){
      Useful.mat = Input.col.mat[c(Current.Row.Combo[,Row.lowest.selection]),]
      #  print(Useful.mat)
      #  Big.counter.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])]=
      #   Big.counter.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])]+1

      Index.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])]=
        Index.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])]+1

      if(dim(Useful.mat)[1]==2 & dim(Useful.mat)[2]==2){
        # Only 2x2 matrix have been filtered. Within these 2x2 matrix we find the OR.
        # Just doing a count of the cells taking only the 2x2 matrix.
        Index.4.cells=which(Index.mat !=0, arr.ind = T)

        # Resetting the Index mat to zero's - we have the cell information in Index.4.cells
        Index.mat = Index.mat - Index.mat

        # Only 2x2 matrix have been filtered. Within these 2x2 matrix we find the OR applying
        # 0.5 cc if there are zero's

        a=Useful.mat[1,1]
        b=Useful.mat[1,2]
        c=Useful.mat[2,1]
        d=Useful.mat[2,2]
        if(a==0 || b==0 || c==0 || d==0){
          a.all=a+0.5
          b.all=b+0.5
          c.all=c+0.5
          d.all=d+0.5
        } else {
          a.all=a
          b.all=b
          c.all=c
          d.all=d
        }# End of adding cc

        ora.all  = (a.all*d.all ) / (b.all*c.all)

        out.full =    data.frame(    Odds.Ratio=ora.all,
                                     a=a.all,
                                     b=b.all,
                                     c=c.all,
                                     d=d.all,
                                     Index1RC=paste(Index.4.cells[1,1],Index.4.cells[1,2],sep=","),
                                     Index2RC=paste(Index.4.cells[3,1],Index.4.cells[3,2],sep=","),
                                     Index3RC=paste(Index.4.cells[2,1],Index.4.cells[2,2],sep=","),
                                     Index4RC=paste(Index.4.cells[4,1],Index.4.cells[4,2],sep=","))
        out.df = rbind(out.df, out.full)
      } # End of if condition checking for 2x2 matrix

    } # End of for loop for lowest Row selection

  } # End of for loop in Loop.col for the current col combination


  return(out.df)
} # End of function
