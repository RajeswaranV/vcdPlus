#'  A partitioning may show that an association reflects primarily differences between certain
#'  categories or groupings of categories. In case of IxJ table, independent chi-squared components
#'  result from comparing columns 1 and 2 and then combining them and comparing them to column 3,
#'  and so on. This generates (I-1)(J-1) unique partioned tables for a given IxJ input.
#'  The Chi-squared and G-sqaured test for each partition is calculated.
#'  Comparision of overall Chi-sqaured and G-squared value is done with the sumation of the
#'  Chi-squared and G-sqaured values of the partioned tables
#' @param mat - matrix for which the sub-matrix is to be generated
#' @param details - If this is set to TRUE, the the return value includes the full list of
#' partioned tables with its Chi-sqaured and G-squared values
#' @details  A partitioning may show that an association reflects primarily differences between certain
#'  categories or groupings of categories. This is based on 2 x 2 tables so that df will be 1 for each table.
#'  In total we generate (I-1)(J-1) tables so that total
#'  df = (I-1)(J-1) which is df of original table
#' @return A list of dataframes with
#'  \item{Chisq.compare}{  Dataframe with the difference between the supporting matrix and the opposing matrix}
#'  \item{Gsq.compare }{ Dataframe of the cell counts for support of table level Chi-sqaured }
#'  \item{Partitioning.df }{ Dataframe of partioned tables with its Chi-squared and G-squared values
#'  - this is returned only if the details flag is set to TRUE }
#' @family Test methods
#' @examples
#'## Example  - Agresti CDA 2002 p82 Section 3.3.3 - Example 3.3.4 is illustrated here.
#' mat=matrix(c(90,12,78,13,1,6,19,13,50),3,3,byrow = TRUE)
#' Partition.table(mat,details=FALSE)
#' @references
#' [1] Alan Agresti
#' Categorical Data Analysis, 2nd Edition
#' John Wiley & Sons
#' @export
Partition.table<-function(mat,details=FALSE){
if (missing(mat)) stop("'mat' is missing")
if ((class(mat) != "matrix"))  stop("'mat' has to be a matrix with minimum 2x2")
if ((dim(mat)[1] < 2) || (dim(mat)[2] <2 )) stop("Matrix has to be minimum of 2x2")
if ((class(details) != "logical")) stop("details has to be a logical vector")

Partitioning.df = NULL

Total.columns=ncol(mat)
Total.rows=nrow(mat)
a=1:Total.rows
b=1:Total.columns


for(Row.iterator in 2:Total.rows){

  for(Column.iterator in 2:Total.columns){

     otemp=gmat(a,b,Row.iterator,Column.iterator,mat)


        out.full = data.frame(a=otemp[1,1],b=otemp[1,2],c=otemp[2,1],d=otemp[2,2],
                              ChiSq=chisq.test(otemp)$statistic,
                              G.Squared=GSquared(otemp)$GSq)
        Partitioning.df = rbind(Partitioning.df, out.full)
        #  write.table(out.full,"D:/Research/CDA/Testing/Full.out.csv",append=TRUE,row.names=FALSE, col.names = FALSE, sep=",")
      } #End of for-loop Column.iterator

    } # End of for loop in Row.iterator

names(Partitioning.df)=c("a","b","c","d","Chisq","Gsq")
Chisq.compare = data.frame(Original.table.Chisq=chisq.test(mat)$statistic, Partition.table.Chisq.sum=sum(Partitioning.df$Chisq))
Gsq.compare = data.frame(Original.table.Gsq=GSquared(mat)$GSq, Partition.table.Gsq.sum=sum(Partitioning.df$Gsq))

rownames(Partitioning.df) <- NULL
if(details)
    {Partitioning.list =list(Chisq.compare,Gsq.compare, Partitioning.df)
    }else {Partitioning.list=list(Chisq.compare,Gsq.compare)}

return(Partitioning.list)
} # End of function
