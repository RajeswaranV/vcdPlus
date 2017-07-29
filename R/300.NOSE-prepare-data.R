#'  Given an 2x2 table with ai, bi, ci and di, if there are any zero cells, a continuity correction
#'  (CC) is added upto 0.5. This function checks if this addition of continuity correction causes
#'  the reversing the Odds Ratio, Relative Risk and Risk Difference.
#'  We provide 2 ways to add this CC. Adding CC to all the cells or only to the zero cells.
#'  The function will identify if there is a reversal and also give the point of reversal.
#' @param ai - Numeric value
#' @param bi - Numeric value
#' @param ci - Numeric value
#' @param di - Numeric value
#' @details  The function only analyzes 2x2 tables where there is at least 1 zero cell.
#' By default it calculates the reversal point (if it exists) for various values of continuity correction
#' (CC) starting from 10^(-8) (suggested by reference [5]) all the way to 0.5. If there is a reversal due to CC
#' the function returns the value above which there will be a reversal for OR, RR and RD.
#' The function adds the CC using 2 methods.
#' In Method 1 - it adds CC to all the cells (called All) which is a widespread practice
#' In Method 2 - it adds CC to only the zero cells (called Only.Zero)
#' If the rerutn values with the dataframe are zero it means that there was no reversal for CC upto 0.5.
#' @return A dataframe with
#'  \item{Reversal.Flag}{  Text description -
#'  "None" if there is no reversal,
#'  "Both" if there is reversal when adding CC to all cells and adding CC to zero cells only,
#'  "Only.All" if there is reversal when adding CC to all cells only and no reversal when adding to zero cells,
#'  "Only.Zero" if there is reversal when adding CC to zero cells only and no reversal when adding CC to all cells }
#'  \item{OR.all.reversal.point}{ Value above which there will be reversal of Odds Ratio when adding CC to all the cells. If this is zero then there is no reversal}
#'  \item{RR.all.reversal.point }{ Value above which there will be reversal of Relative Risk when adding CC to all the cells. If this is zero then there is no reversal }
#'  \item{RD.all.reversal.point}{  Value above which there will be reversal of Risk Difference when adding CC to all the cells. If this is zero then there is no reversal}
#'  \item{OR.zero.reversal.point }{ Value above which there will be reversal of Odds Ratio when adding CC to only the zero cells. If this is zero then there is no reversal}
#'  \item{RR.zero.reversal.point}{ Value above which there will be reversal of Relative Risk when adding CC to only the zero cells. If this is zero then there is no reversal}
#'  \item{RD.zero.reversal.point }{ Value above which there will be reversal of Risk Difference when adding CC to only the zero cells. If this is zero then there is no reversal}
#' @family IxJ Inference methods
#' @seealso \code{\link[metafor]{rma}}  provides OR, RR and RD values for various inputs, with some control
#' of the cc via the "to=" option. However, we suggest using metafor by passing the continuity corrected values
#' to replicate the results of this function, since adding only to zero cell is not an option.
#' @examples
#'## Example 1 - From reference [2] - reversal only when CC is added to the zero cells alone
#' ai=80; bi=1; ci=39; di=0
#' Reversal.point(ai,bi,ci,di)
#'##  Example 2 - From reference [3] -reversal when CC is added
#' ai=1; bi=33; ci=0; di=7
#' Reversal.point(ai,bi,ci,di)
#'##  Example 3 - From reference [4] - no reversal when CC is added
#' ai=0; bi=282; ci=1; di=249
#' Reversal.point(ai,bi,ci,di)
#' @references
#' [1] J Sweeting, Michael, Alexander J Sutton, and Paul C Lambert.
#' "What to add to nothing? Use and avoidance of continuity corrections in meta-analysis of sparse data."
#' Statistics in medicine 23.9 (2004): 1351-1375.
#' [2] Song H, Zhu J, Lu D. Molecular-targeted first-line therapy for advanced gastric cancer.
#' Cochrane Database of Systematic Reviews 2016, Issue 7. Art. No.: CD011461.
#' [3] Lensen SF, Manders M, Nastri CO, Gibreel A, Martins WP, Templer GE, Farquhar C.
#' Endometrial injury for pregnancy following sexual intercourse or intrauterine insemination.
#' Cochrane Database of Systematic Reviews 2016, Issue 6. Art. No.: CD011424.
#' [4] Nidorf SM, Eikelboom JW, Budgeon CA, Thompson PL.
#' Low-dose colchicine for secondary prevention of cardiovascular disease.
#' Journal of the American College of Cardiology 2013;61(4):404-10.
#' @export
Reversal.point<-function(ai,bi,ci,di)
{
  if ((class(ai) != "numeric") || class(bi) != "numeric" ||
    class(ci) != "numeric" || class(di) != "numeric" )  stop("Inputs ai,bi,ci and di have to be numeric")
  if (ai*bi*ci*di>0)  stop("At least one of the inputs ai,bi,ci or di have to be zero")

mrt = data.frame(a=ai,b=bi,c=ci,d=di)

tpm8=1/100000000
tpm4=1/10000

cc.vals = c(tpm8,tpm4,0.0011,0.011,0.0151,0.11,.151,.21,.251,.31,.351,.41,.451,.5)

Tpm8.OR=Calculate.cc.OR(mrt,tpm8)
Point5.OR=Calculate.cc.OR(mrt,0.5)
if(!(is.between(1,Tpm8.OR$ora,Point5.OR$ora)) &
   (!(is.between(1,Tpm8.OR$orz,Point5.OR$orz)))){
  Reversal.point.return=data.frame(Reversal.Flag="None",
                            OR.all.reversal.point=0,
                            RR.all.reversal.point=0,
                            RD.all.reversal.point=0,
                            OR.zero.reversal.point=0,
                            RR.zero.reversal.point=0,
                            RD.zero.reversal.point=0)
  return(Reversal.point.return)
}

  vbn= In.Each.Row(mrt[1,],cc.vals)
  svbn=cbind(mrt[1,],vbn)
  #new.De=ibv.check(svbn)
  Gobject=Find.RS.metric(svbn)

  if(Gobject$RS.all.or ==0 & Gobject$RS.zero.zor==0){
    Reversal.point.return=data.frame(Reversal.Flag="None",
                              OR.all.reversal.point=Gobject$RS.all.or,
                              RR.all.reversal.point=Gobject$RS.all.rr,
                              RD.all.reversal.point=Gobject$RS.all.rd,
                              OR.zero.reversal.point=Gobject$RS.zero.zor,
                              RR.zero.reversal.point=Gobject$RS.zero.zrr,
                              RD.zero.reversal.point=Gobject$RS.zero.zrd)}

  if(Gobject$RS.all.or >0 & Gobject$RS.zero.zor>0){
  Reversal.point.return=data.frame(Reversal.Flag="Both",
                            OR.all.reversal.point=Gobject$RS.all.or,
                            RR.all.reversal.point=Gobject$RS.all.rr,
                            RD.all.reversal.point=Gobject$RS.all.rd,
                            OR.zero.reversal.point=Gobject$RS.zero.zor,
                            RR.zero.reversal.point=Gobject$RS.zero.zrr,
                            RD.zero.reversal.point=Gobject$RS.zero.zrd)}

    if(Gobject$RS.all.or >0 & Gobject$RS.zero.zor==0){
    Reversal.point.return=data.frame(Reversal.Flag="Only.All",
                              OR.all.reversal.point=Gobject$RS.all.or,
                              RR.all.reversal.point=Gobject$RS.all.rr,
                              RD.all.reversal.point=Gobject$RS.all.rd,
                              OR.zero.reversal.point=Gobject$RS.zero.zor,
                              RR.zero.reversal.point=Gobject$RS.zero.zrr,
                              RD.zero.reversal.point=Gobject$RS.zero.zrd)}

  if(Gobject$RS.all.or ==0 & Gobject$RS.zero.zor>0){
    Reversal.point.return=data.frame(Reversal.Flag="Only.Zero",
                              OR.all.reversal.point=Gobject$RS.all.or,
                              RR.all.reversal.point=Gobject$RS.all.rr,
                              RD.all.reversal.point=Gobject$RS.all.rd,
                              OR.zero.reversal.point=Gobject$RS.zero.zor,
                              RR.zero.reversal.point=Gobject$RS.zero.zrr,
                              RD.zero.reversal.point=Gobject$RS.zero.zrd)}


  return(Reversal.point.return)

} # End of the function

