#' Given an input matrix, we can split it into smaller 2x2 sub-matrix and then find
#'  the Nominal Odds Ratio  for each sub-matrix (subtable) keeping a reference row and column
#'  which can be defined by the user. By default the last row and column is taken as the reference.
#' @param mat - matrix for which the sub-matrix is to be generated.
#' @param Reference.Row.ID - The row to be used as reference. By default, the last row is taken as the reference.
#' @param Reference.Col.ID - The column to be used as reference. By default, the last column is taken as the reference.
#' @details This can be used as input to models.
#'   If any of the cells a,b,c or d is zero, a continuity correction of 0.5 is added
#'   to all the cells and the odds ratio is calculated. For nominal classification
#'   variables this set of basic 2x2 tables is defined in terms of a reference category,
#'   usually the cell (I,J).
#'   Then the 2x2 tables formed have in their upper diagonal cell the (i, j) cell of the initial table,
#'   for i= 1 to I -1, j=1 to J-1 and in the lower diagonal cell always the reference cell (I,J). The
#'   non-diagonal cells are the cells of the initial table that share one classification
#'   variable index with each diagonal cell, i.e., they are the cells (i,J) and (I, j).
#' @return A of dataframes with
#'  \item{Odds.Ratio}{  Odds Ratio of the 2x2 subtable}
#'  \item{a }{ The value of a taken from the input table}
#'  \item{b }{ The value of b taken from the input table }
#'  \item{c }{ The value of c taken from the input table }
#'  \item{d }{ The value of d taken from the input table }
#'  \item{Index1RC }{ Index reference from the input matrix for a (Row,Column) }
#'  \item{Index2RC }{ Index reference from the input matrix for b (Row,Column) }
#'  \item{Index3RC }{ Index reference from the input matrix for c (Row,Column) }
#'  \item{Index4RC }{ Index reference from the input matrix for d (Row,Column) }
#' @family IxJ Inference methods
#' @examples
#'##  Here we use the default last column as the reference.
#' mat=matrix(c(1:9), nrow=3, ncol=3)
#' Nominal.Odds.Ratio(mat,Reference.Row.ID=3, Reference.Col.ID=3 ) # Same as Nominal.Odds.Ratio(mat)
#'##  Now we change the reference row and column to 2
#' Nominal.Odds.Ratio(mat, 2, 2)
#' @references
#' [1] Kateri, Maria.
#' Contingency Table Analysis.
#' Springer New York, 2014.
#' @export
Nominal.Odds.Ratio<-function(mat,Reference.Row.ID=NULL, Reference.Col.ID=NULL )
{
  if (missing(mat)) stop("'mat' is missing")
  if ((class(mat) != "matrix"))  stop("'mat' has to be a matrix with minimum 2x2")
  if ((dim(mat)[1] < 2) || (dim(mat)[2] <2 )) stop("Matrix has to be minimum of 2x2")
  if (is.null(Reference.Row.ID)) { Reference.Row.ID = nrow(mat)}
  if (is.null(Reference.Col.ID)) { Reference.Col.ID = ncol(mat)}
  if (((class(Reference.Row.ID) != "numeric" & (class(Reference.Row.ID) != "integer")) ||
       (class(Reference.Col.ID) != "numeric" & (class(Reference.Col.ID) != "integer"))  ||
       Reference.Row.ID%%1!=0 ||
       Reference.Col.ID%%1!=0 ))
    stop("'Reference.Row.ID' and 'Reference.Col.ID' have to be an integer")
  if ( Reference.Row.ID > nrow(mat) ||
       Reference.Row.ID < 1 )  stop("'Reference.Col.ID' has to be between 1 and max number of rows")
  if ((Reference.Col.ID) > ncol(mat) ||
      (Reference.Col.ID) < 1)  stop("'Reference.Col.ID' has to be between 1 and max number of columns")

  if(is.null(Reference.Row.ID)){Reference.Row.ID=nrow(mat)}
  if(is.null(Reference.Col.ID)){Reference.Col.ID=ncol(mat)}

  out.df = NULL

  Total.columns=ncol(mat)
  Total.Rows=nrow(mat)
  Index.mat=mat-mat

  # Take cols first and then find the row elements for each col combination

  Current.Column.Combo.row1=seq(from=1, to=Total.columns)
  Current.Column.Combo.row2=rep(Reference.Col.ID,Total.columns)
  Current.Column.Combo=rbind(Current.Column.Combo.row1,Current.Column.Combo.row2)


  # Now loop over the combinations of columns
  # This loop below will generate the data for the various column combinations
  # Using this we go over each row combination

  for(Loop.col in 1:ncol(Current.Column.Combo)){
    if(Loop.col!=Reference.Col.ID){
      Input.col.mat = mat[,c(Current.Column.Combo[,Loop.col])]
      #print(Input.col.mat)

      # Now generating the row combinations for the given columns
      Current.Row.Combo.row1=seq(from=1,to=Total.Rows)
      Current.Row.Combo.row2=rep(Reference.Row.ID,Total.Rows)
      Current.Row.Combo=rbind(Current.Row.Combo.row1,Current.Row.Combo.row2)

      for(Row.lowest.selection in 1:ncol(Current.Row.Combo)){
        if(Row.lowest.selection!=Reference.Row.ID){
          Useful.mat = Input.col.mat[c(Current.Row.Combo[,Row.lowest.selection]),]
          #  print(Useful.mat)
          #  Big.counter.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])]=
          #   Big.counter.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])]+1

          Index.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])]=
            Index.mat[c(Current.Row.Combo[,Row.lowest.selection]),c(Current.Column.Combo[,Loop.col])]+1
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

          RR.ID.Flag=0

          if(Reference.Row.ID!=Total.Rows & Row.lowest.selection>Reference.Row.ID &
             Reference.Col.ID!=Total.columns & Loop.col>Reference.Col.ID) {
            Index1RC=paste(Index.4.cells[4,1],Index.4.cells[4,2],sep=",")
            Index2RC=paste(Index.4.cells[2,1],Index.4.cells[2,2],sep=",")
            Index3RC=paste(Index.4.cells[3,1],Index.4.cells[3,2],sep=",")
            Index4RC=paste(Index.4.cells[1,1],Index.4.cells[1,2],sep=",")
            RR.ID.Flag=1}

          if(Reference.Row.ID!=Total.Rows & Row.lowest.selection>Reference.Row.ID & RR.ID.Flag<1)
          {
            Index1RC=paste(Index.4.cells[2,1],Index.4.cells[2,2],sep=",")
            Index2RC=paste(Index.4.cells[4,1],Index.4.cells[4,2],sep=",")
            Index3RC=paste(Index.4.cells[1,1],Index.4.cells[1,2],sep=",")
            Index4RC=paste(Index.4.cells[3,1],Index.4.cells[3,2],sep=",")
            RR.ID.Flag=1
          } else {
            if(RR.ID.Flag<1){
              Index1RC=paste(Index.4.cells[1,1],Index.4.cells[1,2],sep=",")
              Index2RC=paste(Index.4.cells[3,1],Index.4.cells[3,2],sep=",")
              Index3RC=paste(Index.4.cells[2,1],Index.4.cells[2,2],sep=",")
              Index4RC=paste(Index.4.cells[4,1],Index.4.cells[4,2],sep=",")
            }
          }

          if(Reference.Col.ID!=Total.columns & Loop.col>Reference.Col.ID & RR.ID.Flag<1)
          {
            Index1RC=paste(Index.4.cells[3,1],Index.4.cells[3,2],sep=",")
            Index2RC=paste(Index.4.cells[1,1],Index.4.cells[1,2],sep=",")
            Index3RC=paste(Index.4.cells[4,1],Index.4.cells[4,2],sep=",")
            Index4RC=paste(Index.4.cells[2,1],Index.4.cells[2,2],sep=",")
          } else {
            if(RR.ID.Flag<1){
              Index1RC=paste(Index.4.cells[1,1],Index.4.cells[1,2],sep=",")
              Index2RC=paste(Index.4.cells[3,1],Index.4.cells[3,2],sep=",")
              Index3RC=paste(Index.4.cells[2,1],Index.4.cells[2,2],sep=",")
              Index4RC=paste(Index.4.cells[4,1],Index.4.cells[4,2],sep=",")
            }
          }



          out.full =    data.frame(    Odds.Ratio=ora.all,
                                       a=a.all,
                                       b=b.all,
                                       c=c.all,
                                       d=d.all,
                                       Index1RC,
                                       Index2RC,
                                       Index3RC,
                                       Index4RC)

          out.df = rbind(out.df, out.full)

        } # End of if checking for Row.lowest.selection!=Reference.Row.ID
      } # End of for loop for lowest Row selection

    } # End of if checking for Ref.Col!=Loop.col
  } # End of for loop in Loop.col for the current col combination

  return(out.df)
} # End of function
