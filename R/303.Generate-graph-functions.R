#'  Given an 2x2 table with ai, bi, ci and di, if there are any zero cells, a continuity correction
#'  (CC) is added upto 0.5. This function checks if this addition of continuity correction causes
#'  the reversing the Odds Ratio, Relative Risk and Risk Difference.
#'  We provide 2 ways to add this CC. Adding CC to all the cells or only to the zero cells.
#'  The function will identify if there is a reversal and also give the point of reversal.
#'  The continuity correction can be also graphed if the option is selected.
#' @param ai - Numeric value
#' @param bi - Numeric value
#' @param ci - Numeric value
#' @param di - Numeric value
#' @param With.CI - Logical TRUE/FALSE. Default is FALSE
#' @details  The function only analyzes 2x2 tables where there is at least 1 zero cell.
#' By default it calculates the reversal point (if it exists) for various values of continuity correction
#' (CC) starting from 10^(-8) (suggested by reference [5]) all the way to 0.5. If there is a reversal due to CC
#' the function returns the graph showing the reversal for Odds Ratio.
#' The function adds the CC using 2 methods.
#' In Method 1 - it adds CC to all the cells (called All) which is a widespread practice
#' In Method 2 - it adds CC to only the zero cells (called Only.Zero)
#' Option is provided to graph the confidence interaval (CI).
#' @return A graph showing the reversal if there is reversal. Or text if there is no reversal
#' @family 2x2 Summary methods
#' @examples
#'## Example 1 - From reference [2] - reversal only when CC is added to the zero cells alone
#' ai=80; bi=1; ci=39; di=0
#' Plot.reversals.OR(ai,bi,ci,di, With.CI=FALSE)
#'##  Example 2 - From reference [3] -reversal when CC is added
#' ai=1; bi=33; ci=0; di=7
#' Plot.reversals.OR(ai,bi,ci,di, With.CI=FALSE)
#'## Example 3 - Above example with CI
#' ai=1; bi=33; ci=0; di=7
#' Plot.reversals.OR(ai,bi,ci,di, With.CI=TRUE)
#' @references
#' [1] J Sweeting, Michael, Alexander J Sutton, and Paul C Lambert.
#' What to add to nothing? Use and avoidance of continuity corrections in meta-analysis of sparse data.
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
Plot.reversals.OR<-function(ai,bi,ci,di,With.CI=FALSE){
#   if ((class(ai) != "numeric") || class(bi) != "numeric" ||
#        class(ci) != "numeric"  || class(di) != "numeric" )  stop("Inputs ai, bi, ci and di have to be numeric")
# if (ai*bi*ci*di>0)  stop("At least one of the inputs ai,bi,ci or di have to be zero")
# if ((class(With.CI)) != "logical") stop("With.CI has to be either TRUE or FALSE")

mrt = data.frame(a=ai,b=bi,c=ci,d=di)
tpm8=1/100000000
tpm4=1/10000

cc.vals = c(tpm8,tpm4,0.0011,0.011,0.0151,0.11,.151,.21,.251,.31,.351,.41,.451,.5)


vbn= In.Each.Row(mrt[1,],cc.vals)
svbn=cbind(mrt[1,],vbn)
ttrange.df=Find.RS.metric(svbn)
if(ttrange.df$RS.all.or ==0 & ttrange.df$RS.zero.zor==0){
return("There is no reversal to Plot")
  }

RS.all.or=ttrange.df$RS.all.or[1]
RS.zero.or=ttrange.df$RS.zero.zor[1]


# Generating sequence and then calling rma for OR/RR/RD with CI
if(RS.all.or>0){RS.all.or.seq=Find.sequence.rma(RS.all.or); orm.all=Find.rma.metric(ai,bi,ci,di,RS.all.or.seq, 1)}
if(RS.zero.or>0){RS.zero.or.seq=Find.sequence.rma(RS.zero.or); orm.zero=Find.rma.metric(ai,bi,ci,di,RS.zero.or.seq, 4)}

cc=Metric=Low=High=NULL
###############################################################################################
######## Taking care of conditions where there is no reversal #################################
if(RS.all.or==0){
  if(With.CI==TRUE){
    #### With CI plot
    pl2 = data.frame( orm.zero)
    names(pl2)=c("cc","Metric", "Low", "High")
    pl2.noci.df=data.frame(orm.zero)
    names(pl2.noci.df)=c("cc","Metric", "Low", "High")
    pl2.noci.df$Low=pl2.noci.df$Metric
    pl2.noci.df$High=pl2.noci.df$Metric
    pl2.noci.df$cigroup="Without CI"
    pl2$cigroup="With CI"
    pl2.noci.df$allz="CC added only to 0 cells"
    pl2$allz="CC added only to 0 cells"
    or.combined.df=rbind(pl2,pl2.noci.df)

    p<-ggplot(or.combined.df,aes(cc,Metric))+geom_line()+
      geom_point(data=or.combined.df,aes(x=cc,y=Metric))+
      geom_line(data=or.combined.df)+
      geom_ribbon(aes(ymin=Low, ymax=High), alpha=0.2) +
      labs(title = "Odds Ratio",subtitle = "No reversal when CC added to all cells") +
      ylab("Odds Ratio") +
      facet_wrap(allz~cigroup, scales = "free" )
    # End of if condition with CI
  print(p)
  } else {
    pl2 = data.frame( orm.zero)
names(pl2)=c("cc","Metric", "Low", "High")

l2=as.character(pl2$cc[6])
cpoint.pl2=pl2$cc[6]

pl2g.orz=pl2[6,]
p<-ggplot(pl2,aes(cc,Metric))+
  geom_line(data=pl2)+
  geom_point(data=pl2,aes(x=cc,y=Metric))+
  geom_point(data=pl2g.orz, colour="blue") +  # this adds a blue point
  geom_text(data=pl2g.orz, label=l2, hjust=0,vjust=0)+ # this adds a label for the blue point
  labs(title = "Odds Ratio",subtitle = "No reversal when CC added to all cells") +
  ylab("Odds Ratio") +
  geom_hline(yintercept=1)
print(p)
  } #End of plot without CI

} # End of if condition when there is only 1 reversal

###############################################################################################
###########Creates a plot for both OR.all and OR.zero without CI ##############################
if(RS.all.or>0 & RS.zero.or>0){
  if(With.CI==TRUE){

    pl = data.frame( orm.all)
    names(pl)=c("cc","Metric", "Low", "High")
    pl2 = data.frame( orm.zero)
    names(pl2)=c("cc","Metric", "Low", "High")

    pl2g.orz=pl2[6,]
    pl1g.or=pl[6,]
    #### With CI plot -pl
    pl.noci.df=data.frame(orm.all)
    names(pl.noci.df)=c("cc","Metric", "Low", "High")
    pl.noci.df$Low=pl.noci.df$Metric
    pl.noci.df$High=pl.noci.df$Metric
    pl.noci.df$cigroup="Without CI"
    pl$cigroup="With CI"
    pl.noci.df$allz="CC added to all cells"
    pl$allz="CC added to all cells"

    #### With CI plot - pl2
    pl2.noci.df=data.frame(orm.zero)
    names(pl2.noci.df)=c("cc","Metric", "Low", "High")
    pl2.noci.df$Low=pl2.noci.df$Metric
    pl2.noci.df$High=pl2.noci.df$Metric
    pl2.noci.df$cigroup="Without CI"
    pl2$cigroup="With CI"
    pl2.noci.df$allz="CC added only to 0 cells"
    pl2$allz="CC added only to 0 cells"

    or.combined.df=rbind(pl,pl.noci.df,pl2,pl2.noci.df)

    p<-ggplot(or.combined.df,aes(cc,Metric))+geom_line()+
      geom_point(data=or.combined.df,aes(x=cc,y=Metric))+
      geom_line(data=or.combined.df)+
      geom_ribbon(aes(ymin=Low, ymax=High), alpha=0.2) +
      labs(title = "Odds Ratio") +
      ylab("Odds Ratio") +
      facet_grid(cigroup~allz, scales = "free" )
    # End of if condition with CI
    print(p)
    } else {
    pl = data.frame( orm.all)
    names(pl)=c("cc","Metric", "Low", "High")
    pl2 = data.frame( orm.zero)
    names(pl2)=c("cc","Metric", "Low", "High")

    l1=as.character(pl$cc[6])
    l2=as.character(pl2$cc[6])
    cpoint.pl=pl$cc[6]
    cpoint.pl2=pl2$cc[6]

    pl2g.orz=pl2[6,]
    pl1g.or=pl[6,]
    p<-ggplot(pl,aes(cc,Metric))+geom_line(aes(color="CC added to all cells"))+
      geom_point(data=pl,aes(x=cc,y=Metric))+
      geom_line(data=pl2,aes(color="CC added to only the 0 cells"))+
      geom_point(data=pl2,aes(x=cc,y=Metric))+
      geom_point(data=pl2g.orz, colour="blue") +  # this adds a blue point
      geom_text(data=pl2g.orz, label=l2, hjust=0,vjust=0)+ # this adds a label for the blue point
      geom_point(data=pl1g.or, colour="red") +  # this adds a red point
      geom_text(data=pl1g.or, label=l1, hjust=0,vjust=0)+ # this adds a label for the red point
      labs(title = "Odds Ratio") +
      ylab("Odds Ratio") +
      geom_hline(yintercept=1)

    print(p)
    } #End of plot without CI


  } # End of if condition when there is both reversal

} # End of OR function

#############################################################################################################
#'  Given an 2x2 table with ai, bi, ci and di, if there are any zero cells, a continuity correction
#'  (CC) is added upto 0.5. This function checks if this addition of continuity correction causes
#'  the reversing the Odds Ratio, Relative Risk and Risk Difference.
#'  We provide 2 ways to add this CC. Adding CC to all the cells or only to the zero cells.
#'  The function will identify if there is a reversal and also give the point of reversal.
#'  The continuity correction can be also graphed if the option is selected.
#' @param ai - Numeric value
#' @param bi - Numeric value
#' @param ci - Numeric value
#' @param di - Numeric value
#' @param With.CI - Logical TRUE/FALSE. Default is FALSE
#' @details  The function only analyzes 2x2 tables where there is at least 1 zero cell.
#' By default it calculates the reversal point (if it exists) for various values of continuity correction
#' (CC) starting from 10^(-8) (suggested by reference [5]) all the way to 0.5. If there is a reversal due to CC
#' the function returns the graph showing the reversal for Relative Risk.
#' The function adds the CC using 2 methods.
#' In Method 1 - it adds CC to all the cells (called All) which is a widespread practice
#' In Method 2 - it adds CC to only the zero cells (called Only.Zero)
#' Option is provided to graph the confidence interaval (CI).
#' @return A graph showing the reversal if there is reversal. Or text if there is no reversal
#' @family 2x2 Summary methods
#' @examples
#'## Example 1 - From reference [2] - reversal only when CC is added to the zero cells alone
#' ai=80; bi=1; ci=39; di=0
#' Plot.reversals.RR(ai,bi,ci,di, With.CI=FALSE)
#'##  Example 2 - From reference [3] -reversal when CC is added
#' ai=1; bi=33; ci=0; di=7
#' Plot.reversals.RR(ai,bi,ci,di, With.CI=FALSE)
#'## Example 3 - Above example with CI
#' ai=1; bi=33; ci=0; di=7
#' Plot.reversals.RR(ai,bi,ci,di, With.CI=TRUE)
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
Plot.reversals.RR<-function(ai,bi,ci,di,With.CI=FALSE){
  if ((class(ai) != "numeric") || class(bi) != "numeric" ||
      class(ci) != "numeric"  || class(di) != "numeric" )  stop("Inputs ai,bi,ci and di have to be numeric")
  if (ai*bi*ci*di>0)  stop("At least one of the inputs ai,bi,ci or di have to be zero")
  if ((class(With.CI)) != "logical") stop("With.CI has to be either TRUE or FALSE")

  mrt = data.frame(a=ai,b=bi,c=ci,d=di)
  tpm8=1/100000000
  tpm4=1/10000
  cc=Metric=Low=High=NULL

  cc.vals = c(tpm8,tpm4,0.0011,0.011,0.0151,0.11,.151,.21,.251,.31,.351,.41,.451,.5)


  vbn= In.Each.Row(mrt[1,],cc.vals)
  svbn=cbind(mrt[1,],vbn)
  #new.De=ibv.check(svbn)
  ttrange.df=Find.RS.metric(svbn)
  if(ttrange.df$RS.all.rr ==0 & ttrange.df$RS.zero.zrr==0){
    return("There is no reversal to Plot")
  }

  RS.all.rr=ttrange.df$RS.all.rr[1]
  RS.zero.rr=ttrange.df$RS.zero.zrr[1]

  # Generating sequence and then calling rma for OR/RR/RD with CI
  if(RS.all.rr>0){RS.all.rr.seq=Find.sequence.rma(RS.all.rr); rrm.all=Find.rma.metric(ai,bi,ci,di,RS.all.rr.seq, 2)}
  if(RS.zero.rr>0){RS.zero.rr.seq=Find.sequence.rma(RS.zero.rr); rrm.zero=Find.rma.metric(ai,bi,ci,di,RS.zero.rr.seq, 5)}

###############################################################################################
######## Taking care of conditions where there is no reversal #################################
  if(RS.all.rr==0){
    if(With.CI==TRUE){
      #### With CI plot
      rr2=data.frame(rrm.zero)
      names(rr2)=c("cc","Metric", "Low", "High")
      rr2.noci.df=data.frame(rrm.zero)
      names(rr2.noci.df)=c("cc","Metric", "Low", "High")
      rr2.noci.df$Low=rr2.noci.df$Metric
      rr2.noci.df$High=rr2.noci.df$Metric
      rr2.noci.df$cigroup="Without CI"
      rr2$cigroup="With CI"
      rr2.noci.df$allz="CC added only to 0 cells"
      rr2$allz="CC added only to 0 cells"
      rr.combined.df=rbind(rr2,rr2.noci.df)

      p<-ggplot(rr.combined.df,aes(cc,Metric))+geom_line()+
        geom_point(data=rr.combined.df,aes(x=cc,y=Metric))+
        geom_line(data=rr.combined.df)+
        geom_ribbon(aes(ymin=Low, ymax=High), alpha=0.2) +
        labs(title = "Relative Risk", subtitle = "No reversal when CC added to all cells") +
        ylab("Relative Risk") +
        facet_wrap(allz~cigroup, scales = "free" )
      # End of if condition with CI
      print(p)
      } else {

    rr2 = data.frame( rrm.zero)
    names(rr2)=c("cc","Metric", "Low", "High")

    l2.rr=as.character(rr2$cc[6])
    point.rr2=rr2$cc[6]
    rr2g.rrz=rr2[6,]
    p<-ggplot(rr2,aes(cc,Metric))+
      geom_line(data=rr2)+
      geom_point(data=rr2,aes(x=cc,y=Metric))+
      geom_point(data=rr2g.rrz, colour="blue") +  # this adds a blue point
      geom_text(data=rr2g.rrz, label=l2.rr, hjust=0,vjust=0)+ # this adds a label for the blue point
      labs(title = "Relative Risk", subtitle = "No reversal when CC added to all cells") +
      ylab("Relative Risk") +
      geom_hline(yintercept=1)
    print(p)

    } #End of plot without CI

} # End of if condition when there is only 1 reversal


  ###############################################################################################
  ###########Creates a plot for both RR.all and RR.zero without CI ##############################
if(RS.all.rr>0 & RS.zero.rr>0){
  if(With.CI==TRUE){

    #### With CI plot - rr
    rr = data.frame( rrm.all)
    names(rr)=c("cc","Metric", "Low", "High")
    rr.noci.df=data.frame(rrm.all)
    names(rr.noci.df)=c("cc","Metric", "Low", "High")
    rr.noci.df$Low=rr.noci.df$Metric
    rr.noci.df$High=rr.noci.df$Metric
    rr.noci.df$cigroup="Without CI"
    rr$cigroup="With CI"
    rr.noci.df$allz="CC added to all cells"
    rr$allz="CC added to all cells"

    #### With CI plot - rr2
    rr2=data.frame(rrm.zero)
    names(rr2)=c("cc","Metric", "Low", "High")
    rr2.noci.df=data.frame(rrm.zero)
    names(rr2.noci.df)=c("cc","Metric", "Low", "High")
    rr2.noci.df$Low=rr2.noci.df$Metric
    rr2.noci.df$High=rr2.noci.df$Metric
    rr2.noci.df$cigroup="Without CI"
    rr2$cigroup="With CI"
    rr2.noci.df$allz="CC added only to 0 cells"
    rr2$allz="CC added only to 0 cells"

    rr.combined.df=rbind(rr,rr.noci.df,rr2,rr2.noci.df)

    p<-ggplot(rr.combined.df,aes(cc,Metric))+geom_line()+
      geom_point(data=rr.combined.df,aes(x=cc,y=Metric))+
      geom_line(data=rr.combined.df)+
      geom_ribbon(aes(ymin=Low, ymax=High), alpha=0.2) +
      labs(title = "Relative Risk") +
      ylab("Relative Risk") +
      facet_grid(cigroup~allz, scales = "free" )
    # End of if condition with CI
    print(p)
    } else {

  rr = data.frame( rrm.all)
  names(rr)=c("cc","Metric", "Low", "High")
  rr2 = data.frame( rrm.zero)
  names(rr2)=c("cc","Metric", "Low", "High")

  l1.rr=as.character(rr$cc[6])
  l2.rr=as.character(rr2$cc[6])
  cpoint.rr=rr2$cc[6]
  cpoint.rr2=rr2$cc[6]
  rr2g.rrz=rr2[6,]
  rr1g.rr=rr[6,]
  p<-ggplot(rr,aes(cc,Metric))+geom_line(aes(color="CC added to all cells"))+
    geom_point(data=rr,aes(x=cc,y=Metric))+
    geom_line(data=rr2,aes(color="CC added to only the 0 cells"))+
    geom_point(data=rr2,aes(x=cc,y=Metric))+
    geom_point(data=rr2g.rrz, colour="blue") +  # this adds a blue point
    geom_text(data=rr2g.rrz, label=l2.rr, hjust=0,vjust=0)+ # this adds a label for the blue point
    geom_point(data=rr1g.rr, colour="red") +  # this adds a red point
    geom_text(data=rr1g.rr, label=l1.rr, hjust=0,vjust=0)+ # this adds a label for the red point
    labs(title = "Relative Risk") +
    ylab("Relative Risk") +
    geom_hline(yintercept=1)

  print(p)
    } #End of plot without CI

  } # End of if condition when there is both reversal

} # End of RR function


####################################################################################################
#'  Given an 2x2 table with ai, bi, ci and di, if there are any zero cells, a continuity correction
#'  (CC) is added upto 0.5. This function checks if this addition of continuity correction causes
#'  the reversing the Odds Ratio, Relative Risk and Risk Difference.
#'  We provide 2 ways to add this CC. Adding CC to all the cells or only to the zero cells.
#'  The function will identify if there is a reversal and also give the point of reversal.
#'  The continuity correction can be also graphed if the option is selected.
#' @param ai - Numeric value
#' @param bi - Numeric value
#' @param ci - Numeric value
#' @param di - Numeric value
#' @param With.CI - Logical TRUE/FALSE. Default is FALSE
#' @details  The function only analyzes 2x2 tables where there is at least 1 zero cell.
#' By default it calculates the reversal point (if it exists) for various values of continuity correction
#' (CC) starting from 10^(-8) (suggested by reference [5]) all the way to 0.5. If there is a reversal due to CC
#' the function returns the graph showing the reversal for Risk Difference.
#' The function adds the CC using 2 methods.
#' In Method 1 - it adds CC to all the cells (called All) which is a widespread practice
#' In Method 2 - it adds CC to only the zero cells (called Only.Zero)
#' Option is provided to graph the confidence interaval (CI).
#' @return A graph showing the reversal if there is reversal. Or text if there is no reversal
#' @family 2x2 Summary methods
#' @examples
#'## Example 1 - From reference [2] - reversal only when CC is added to the zero cells alone
#' ai=80; bi=1; ci=39; di=0
#' Plot.reversals.RD(ai,bi,ci,di, With.CI=FALSE)
#'##  Example 2 - From reference [3] -reversal when CC is added
#' ai=1; bi=33; ci=0; di=7
#' Plot.reversals.RD(ai,bi,ci,di, With.CI=FALSE)
#'## Example 3 - Above example with CI
#' ai=1; bi=33; ci=0; di=7
#' Plot.reversals.RD(ai,bi,ci,di, With.CI=TRUE)
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
Plot.reversals.RD<-function(ai,bi,ci,di,With.CI=FALSE){
  if ((class(ai) != "numeric") || class(bi) != "numeric" ||
      class(ci) != "numeric"  || class(di) != "numeric" )  stop("Inputs ai,bi,ci and di have to be numeric")
  if (ai*bi*ci*di>0)  stop("At least one of the inputs ai,bi,ci or di have to be zero")
  if ((class(With.CI)) != "logical") stop("With.CI has to be either TRUE or FALSE")

  mrt = data.frame(a=ai,b=bi,c=ci,d=di)
  tpm8=1/100000000
  tpm4=1/10000
  cc=Metric=Low=High=NULL

  cc.vals = c(tpm8,tpm4,0.0011,0.011,0.0151,0.11,.151,.21,.251,.31,.351,.41,.451,.5)


  vbn= In.Each.Row(mrt[1,],cc.vals)
  svbn=cbind(mrt[1,],vbn)
  ttrange.df=Find.RS.metric(svbn)
  if(ttrange.df$RS.all.rd ==0 & ttrange.df$RS.zero.zrd==0){
    return("There is no reversal to Plot")
  }

  RS.all.rd=ttrange.df$RS.all.rd[1]
  RS.zero.rd=ttrange.df$RS.zero.zrd[1]

  # Generating sequence and then calling rma for OR/RR/RD with CI
  if(RS.all.rd>0){RS.all.rd.seq=Find.sequence.rma(RS.all.rd); rdm.all=Find.rma.metric(ai,bi,ci,di,RS.all.rd.seq, 3)}
  if(RS.zero.rd>0){RS.zero.rd.seq=Find.sequence.rma(RS.zero.rd); rdm.zero=Find.rma.metric(ai,bi,ci,di,RS.zero.rd.seq, 6)}

if(RS.all.rd==0){
  if(With.CI==TRUE){
    #### With CI plot
    rd2 = data.frame( rdm.zero)
    names(rd2)=c("cc","Metric", "Low", "High")
    rd2.noci.df=data.frame(rdm.zero)
    names(rd2.noci.df)=c("cc","Metric", "Low", "High")
    rd2.noci.df$Low=rd2.noci.df$Metric
    rd2.noci.df$High=rd2.noci.df$Metric
    rd2.noci.df$cigroup="Without CI"
    rd2$cigroup="With CI"
    rd2.noci.df$allz="CC added only to 0 cells"
    rd2$allz="CC added only to 0 cells"
    rd.combined.df=rbind(rd2,rd2.noci.df)

    p<-ggplot(rd.combined.df,aes(cc,Metric))+geom_line()+
      geom_point(data=rd.combined.df,aes(x=cc,y=Metric))+
      geom_line(data=rd.combined.df)+
      geom_ribbon(aes(ymin=Low, ymax=High), alpha=0.2) +
      labs(title = "Risk Difference", subtitle = "No reversal when CC added to all cells") +
      ylab("Risk Difference") +
      facet_wrap(allz~cigroup, scales = "free" )
  # End of if condition with CI
    print(p)
    } else {
  rd2 = data.frame( rdm.zero)
  names(rd2)=c("cc","Metric", "Low", "High")

  l2.rd=as.character(rd2$cc[6])
  cpoint.rd2=rd2$cc[6]

  rd2g.rdz=rd2[6,]

  p<-ggplot(rd2,aes(cc,Metric))+
    geom_line(data=rd2)+
    geom_point(data=rd2,aes(x=cc,y=Metric))+
    geom_point(data=rd2g.rdz, colour="blue") +  # this adds a blue point
    geom_text(data=rd2g.rdz, label=l2.rd,  hjust=0,vjust=0)+ # this adds a label for the blue point
    labs(title = "Risk Difference", subtitle = "No reversal when CC added to all cells") +
    ylab("Risk Difference") +
    geom_hline(yintercept=0)

  print(p)
  } #End of plot without CI
} # End of if condition when there is only 1 reversal



###############################################################################################
###########Creates a plot for both RD.all and RD.zero without CI ##############################
if(RS.all.rd>0 & RS.zero.rd>0){
  if(With.CI==TRUE){
    #### With CI plot -rd
    rd = data.frame( rdm.all)
    names(rd)=c("cc","Metric", "Low", "High")
    rd.noci.df=data.frame(rdm.all)
    names(rd.noci.df)=c("cc","Metric", "Low", "High")
    rd.noci.df$Low=rd.noci.df$Metric
    rd.noci.df$High=rd.noci.df$Metric
    rd.noci.df$cigroup="Without CI"
    rd$cigroup="With CI"
    rd.noci.df$allz="CC added to all cells"
    rd$allz="CC added to all cells"

    #### With CI plot -rd2
    rd2 = data.frame( rdm.zero)
    names(rd2)=c("cc","Metric", "Low", "High")
    rd2.noci.df=data.frame(rdm.zero)
    names(rd2.noci.df)=c("cc","Metric", "Low", "High")
    rd2.noci.df$Low=rd2.noci.df$Metric
    rd2.noci.df$High=rd2.noci.df$Metric
    rd2.noci.df$cigroup="Without CI"
    rd2$cigroup="With CI"
    rd2.noci.df$allz="CC added only to 0 cells"
    rd2$allz="CC added only to 0 cells"

    rd.combined.df=rbind(rd,rd.noci.df,rd2,rd2.noci.df)

    p<-ggplot(rd.combined.df,aes(cc,Metric))+geom_line()+
      geom_point(data=rd.combined.df,aes(x=cc,y=Metric))+
      geom_line(data=rd.combined.df)+
      geom_ribbon(aes(ymin=Low, ymax=High), alpha=0.2) +
      labs(title = "Risk Difference") +
      ylab("Risk Difference") +
      facet_grid(cigroup~allz, scales = "free" )
    # End of if condition with CI
    print(p)
    } else {

  rd = data.frame( rdm.all)
  names(rd)=c("cc","Metric", "Low", "High")
  rd2 = data.frame( rdm.zero)
  names(rd2)=c("cc","Metric", "Low", "High")
  # rd2$group<-"NI"
  # rd2$group[6]<-"Important"

  l1.rd=as.character(rd$cc[6])
  l2.rd=as.character(rd2$cc[6])
  cpoint.rd=rd$cc[6]
  cpoint.rd2=rd2$cc[6]

  rd2g.rdz=rd2[6,]
  rd1g.rd=rd[6,]

  p<-ggplot(rd,aes(cc,Metric))+geom_line(aes(color="CC added to all cells"))+
    geom_point(data=rd,aes(x=cc,y=Metric))+
    geom_line(data=rd2,aes(color="CC added to only the 0 cells"))+
    geom_point(data=rd2,aes(x=cc,y=Metric))+
    geom_point(data=rd2g.rdz, colour="blue") +  # this adds a blue point
    geom_text(data=rd2g.rdz, label=l2.rd,  hjust=0,vjust=0)+ # this adds a label for the blue point
    geom_point(data=rd1g.rd, colour="red") +  # this adds a red point
    geom_text(data=rd1g.rd, label=l1.rd,  hjust=0,vjust=0)+ # this adds a label for the red point
    labs(title = "Risk Difference") +
    ylab("Risk Difference") +
    geom_hline(yintercept=0)

  print(p)
    } #End of plot without CI

  } # End of if condition when there is both reversal

} # End of RD function

