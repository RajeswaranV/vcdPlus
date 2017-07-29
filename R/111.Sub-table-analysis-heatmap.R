#'  Given an input matrix, we can split it into smaller sub-matrix (min 2x2) and then find
#'  the Chi-squared test for each sub-matrix. The smaller matrix can "support"  or "oppose"
#'  (have a different conclusion at 95% confidence interaval) compared with the overall Chi-squared
#'  test value of the full input matrix. We count the number of times each cell supports or
#'  opposes the overall Chi-squared test. We also generate the possible list of sub-matrix.
#' @param mat - matrix for which the sub-matrix is to be generated
#' @details  This can be used as an outlier detection method as well as observing the individual
#' cells within an IxJ table
#' @return A heatmap plot
#' @family IxJ Inference methods
#' @examples
#' Drills=c(2,	10,	4,	2 ) # Example data from [reference 1]
#' Pots= c(3,	8,	4,	6)
#' Grinding.Stones=c( 13, 5, 3, 9)
#' Point.Fragments=c(20, 36, 19, 20)
#' mat=rbind(Drills,Pots,Grinding.Stones,Point.Fragments)
#' Plot.heatmap(mat)
#' @references
#' [1] Mosteller F, Parunak A (2006)
#' Identifying extreme cells in a sizable contingency table: Probabilistic and exploratory approaches.
#' In: Hoaglin DC, Mosteller F, Tukey JW (eds) Exploring Data Tables, Trends, and Shapes,
#' John Wiley & Sons, pp 189-224
#' [2] Kuhnt, S., Rapallo, F. & Rehage,
#' Outlier detection in contingency tables based on minimal patterns
#' A. Stat Comput (2014) 24: 481. doi:10.1007/s11222-013-9382-8
#' @export
Plot.heatmap<-function(mat)
{
  if (missing(mat)) stop("'mat' is missing")
  if ((class(mat) != "matrix"))  stop("'mat' has to be a matrix with minimum 2x2")
  if ((dim(mat)[1] < 2) || (dim(mat)[2] <2 )) stop("Matrix has to be minimum of 2x2")

  Var1=Var2=value=NULL

  juju= generate.heatmap.matrix(mat,details = FALSE)
  ju = as.matrix(juju[[1]])
  dat2 <- ju %>%
  tbl_df() %>%
  rownames_to_column('Var1') %>%
  gather(Var2, value, -Var1) %>%
  mutate(
    Var1 = factor(Var1, levels=1:ncol(ju)),
    Var2 = factor(gsub("Heatmap.", "", Var2), levels=1:nrow(ju))
  )

## plot data
ggplot(dat2, aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white", high = "red")

} # End of function
