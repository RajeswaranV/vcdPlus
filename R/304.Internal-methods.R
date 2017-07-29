####################################################################################################
#'  This function  finds the sequence for graph
#' @param RS.seq - Sequence of values for cc
#' @return A dataframe with
#'  \item{Ret.seq}{Sequence}
#' @family Internal methods
#' @keywords internal
Find.sequence.rma<- function(RS.val) {
  small=RS.val-0.1
  big=RS.val+0.1
  New.seq=seq(small, big, length.out = 10)
  Ret.seq=sort(c(RS.val,New.seq),decreasing = FALSE)
  return(Ret.seq)
}

####################################################################################################
#'  This function finds the rma OR/RR/RD for the given  range using metafor
#' @param a - description
#' @param b - description
#' @param c - description
#' @param d - description
#' @param RS.seq - description
#' @param metric - description
#' @return A dataframe with
#'  \item{output}{Sequence}
#' @family Internal methods
#' @keywords internal
Find.rma.metric<- function(a,b,c,d,RS.seq, metric) {
  iterations = length(RS.seq)
  variables = 4

  output <- matrix(ncol=variables, nrow=iterations)

  for(gseq in 1:iterations){
    cc=RS.seq[gseq]
    a.all = a + cc
    b.all = b + cc
    c.all = c + cc
    d.all = d + cc

    if(a==0) {a.zero = a + cc} else {a.zero=a}
    if(b==0) {b.zero = b + cc} else {b.zero=b}
    if(c==0) {c.zero = c + cc} else {c.zero=c}
    if(d==0) {d.zero = d + cc} else {d.zero=d}

    if(metric==1){  output[gseq,] <- c(cc,
                                       exp(metafor::rma(ai=a.all,bi=b.all,ci=c.all,di=d.all,measure="OR")$b),
                                       exp(metafor::rma(ai=a.all,bi=b.all,ci=c.all,di=d.all,measure="OR")$ci.lb),
                                       exp(metafor::rma(ai=a.all,bi=b.all,ci=c.all,di=d.all,measure="OR")$ci.ub))
    }
    if(metric==2){
      output[gseq,] <- c(cc,
                         exp(metafor::rma(ai=a.all,bi=b.all,ci=c.all,di=d.all,measure="RR")$b),
                         exp(metafor::rma(ai=a.all,bi=b.all,ci=c.all,di=d.all,measure="RR")$ci.lb),
                         exp(metafor::rma(ai=a.all,bi=b.all,ci=c.all,di=d.all,measure="RR")$ci.ub))
    }
    if(metric==3){
      output[gseq,] <- c(cc,
                         metafor::rma(ai=a.all,bi=b.all,ci=c.all,di=d.all,measure="RD")$b,
                         metafor::rma(ai=a.all,bi=b.all,ci=c.all,di=d.all,measure="RD")$ci.lb,
                         metafor::rma(ai=a.all,bi=b.all,ci=c.all,di=d.all,measure="RD")$ci.ub)
    } # End of if condition testing for metric type

    if(metric==4){  output[gseq,] <- c(cc,
                                       exp(metafor::rma(ai=a.zero,bi=b.zero,ci=c.zero,di=d.zero,measure="OR")$b),
                                       exp(metafor::rma(ai=a.zero,bi=b.zero,ci=c.zero,di=d.zero,measure="OR")$ci.lb),
                                       exp(metafor::rma(ai=a.zero,bi=b.zero,ci=c.zero,di=d.zero,measure="OR")$ci.ub))
    }
    if(metric==5){
      output[gseq,] <- c(cc,
                         exp(metafor::rma(ai=a.zero,bi=b.zero,ci=c.zero,di=d.zero,measure="RR")$b),
                         exp(metafor::rma(ai=a.zero,bi=b.zero,ci=c.zero,di=d.zero,measure="RR")$ci.lb),
                         exp(metafor::rma(ai=a.zero,bi=b.zero,ci=c.zero,di=d.zero,measure="RR")$ci.ub))
    }
    if(metric==6){
      output[gseq,] <- c(cc,
                         metafor::rma(ai=a.zero,bi=b.zero,ci=c.zero,di=d.zero,measure="RD")$b,
                         metafor::rma(ai=a.zero,bi=b.zero,ci=c.zero,di=d.zero,measure="RD")$ci.lb,
                         metafor::rma(ai=a.zero,bi=b.zero,ci=c.zero,di=d.zero,measure="RD")$ci.ub)
    } # End of if condition testing for metric type

  } # End of inner for loop iterating over gseq (number of cc to run through)

  return(output)

} # End of function

####################################################################################################
#'  This Function to calculate OR/RR/RD and the LHS/RHS of change
#' @param odf.row - description
#' @param cc.vals - description
#' @return A dataframe with
#'  \item{Full.Slevel.df}{Sequence}
#' @family Internal methods
#' @keywords internal
In.Each.Row <-function(odf.row,cc.vals)
{

  iterations = length(cc.vals)

  output.OR <- matrix(ncol=2, nrow=iterations)
  output.RR <- matrix(ncol=2, nrow=iterations)
  output.RD <- matrix(ncol=2, nrow=iterations)

  for(j in 1:iterations){

    cc.value = cc.vals[j]

    ora.results=Calculate.cc.OR(odf.row, cc.value )
    output.OR[j,] <- as.matrix(ora.results)

    rra.results=Calculate.cc.RR(odf.row, cc.value )
    output.RR[j,] <- as.matrix(rra.results)

    rda.results=Calculate.cc.RD(odf.row, cc.value )
    output.RD[j,] <- as.matrix(rda.results)

  }

  odf.or=data.frame(t(output.OR[,1]),t(output.OR[,2]))
  odf.rr=data.frame(t(output.RR[,1]),t(output.RR[,2]))
  odf.rd=data.frame(t(output.RD[,1]),t(output.RD[,2]))

  names(odf.or)<-c("or.cc.all.tpm8","or.cc.all.tpm4",
                   "or.cc.all.001","or.cc.all.01","or.cc.all.015","or.cc.all.1","or.cc.all.15",
                   "or.cc.all.2","or.cc.all.25","or.cc.all.3","or.cc.all.35",
                   "or.cc.all.4","or.cc.all.45","or.cc.all.5",
                   "or.cc.zero.tpm8","or.cc.zero.tpm4",
                   "or.cc.zero.001","or.cc.zero.01","or.cc.zero.015","or.cc.zero.1","or.cc.zero.15",
                   "or.cc.zero.2","or.cc.zero.25","or.cc.zero.3","or.cc.zero.35",
                   "or.cc.zero.4","or.cc.zero.45","or.cc.zero.5")
  Smallest.val.all.or= min(output.OR[,1])
  Largest.val.all.or= max(output.OR[,1])
  Smallest.val.zero.or= min(output.OR[,2])
  Largest.val.zero.or= max(output.OR[,2])

  names(odf.rr)<-c("rr.cc.all.tpm8","rr.cc.all.tpm4",
                   "rr.cc.all.001","rr.cc.all.01","rr.cc.all.015","rr.cc.all.1","rr.cc.all.15",
                   "rr.cc.all.2","rr.cc.all.25","rr.cc.all.3","rr.cc.all.35",
                   "rr.cc.all.4","rr.cc.all.45","rr.cc.all.5",
                   "rr.cc.zero.tpm8","rr.cc.zero.tpm4",
                   "rr.cc.zero.001","rr.cc.zero.01","rr.cc.zero.015","rr.cc.zero.1","rr.cc.zero.15",
                   "rr.cc.zero.2","rr.cc.zero.25","rr.cc.zero.3","rr.cc.zero.35",
                   "rr.cc.zero.4","rr.cc.zero.45","rr.cc.zero.5")
  Smallest.val.all.rr= min(output.RR[,1])
  Largest.val.all.rr= max(output.RR[,1])
  Smallest.val.zero.rr= min(output.RR[,2])
  Largest.val.zero.rr= max(output.RR[,2])

  names(odf.rd)<-c("rd.cc.all.tpm8","rd.cc.all.tpm4",
                   "rd.cc.all.001","rd.cc.all.01","rd.cc.all.015","rd.cc.all.1","rd.cc.all.15",
                   "rd.cc.all.2","rd.cc.all.25","rd.cc.all.3","rd.cc.all.35",
                   "rd.cc.all.4","rd.cc.all.45","rd.cc.all.5",
                   "rd.cc.zero.tpm8","rd.cc.zero.tpm4",
                   "rd.cc.zero.001","rd.cc.zero.01","rd.cc.zero.015","rd.cc.zero.1","rd.cc.zero.15",
                   "rd.cc.zero.2","rd.cc.zero.25","rd.cc.zero.3","rd.cc.zero.35",
                   "rd.cc.zero.4","rd.cc.zero.45","rd.cc.zero.5")
  Smallest.val.all.rd= min(output.RD[,1])
  Largest.val.all.rd= max(output.RD[,1])
  Smallest.val.zero.rd= min(output.RD[,2])
  Largest.val.zero.rd= max(output.RD[,2])


  Full.Slevel.df= data.frame(odf.or,Max.or.all= Largest.val.all.or ,Min.or.all= Smallest.val.all.or,
                             Max.or.zero= Largest.val.zero.or ,Min.or.zero= Smallest.val.zero.or,
                             odf.rr,Max.rr.all= Largest.val.all.rr ,Min.rr.all= Smallest.val.all.rr,
                             Max.rr.zero= Largest.val.zero.rr ,Min.rr.zero= Smallest.val.zero.rr,
                             odf.rd,Max.rd.all= Largest.val.all.rd ,Min.rd.all= Smallest.val.all.rd,
                             Max.rd.zero= Largest.val.zero.rd ,Min.rd.zero= Smallest.val.zero.rd)
  return(Full.Slevel.df)
}

############################################################################################
#'  This Function to calculate OR and the point of change
#' @param cc.df - description
#' @param cc.value - description
#' @return A dataframe with
#'  \item{r.df}{Sequence}
#' @family Internal methods
#' @keywords internal
Calculate.cc.OR <-function(cc.df,cc.value)
{
  a=cc.df$a; b=cc.df$b;
  c=cc.df$c; d=cc.df$d;

  cc=cc.value

  a.all = a + cc
  b.all = b + cc
  c.all = c + cc
  d.all = d + cc

  ################ Calculat only for cells with zero

  if(a==0) {a.zero = a + cc} else {a.zero=a}
  if(b==0) {b.zero = b + cc} else {b.zero=b}
  if(c==0) {c.zero = c + cc} else {c.zero=c}
  if(d==0) {d.zero = d + cc} else {d.zero=d}

  ora.all  = (a.all*d.all)/(b.all*c.all)
  ora.zero  = (a.zero*d.zero)  / (b.zero*c.zero)

  r.df=data.frame(ora=ora.all,orz=ora.zero)
  return(r.df)
}


############################################################################################
#'  This Function to calculate RR and the point of change
#' @param cc.df - description
#' @param cc.value - description
#' @return A dataframe with
#'  \item{r.df}{Sequence}
#' @family Internal methods
#' @keywords internal
Calculate.cc.RR <-function(cc.df,cc.value)
{
  a=cc.df$a; b=cc.df$b;
  c=cc.df$c; d=cc.df$d;

  cc=cc.value
  nt=a+b
  nc=c+d

  a.all = a + cc
  b.all = b + cc
  c.all = c + cc
  d.all = d + cc

  nc.all = nc + (2 * cc)
  nt.all = nt + (2 * cc)

  ################ Calculat only for cells with zero

  if(a==0) {a.zero = a + cc} else {a.zero=a}
  if(b==0) {b.zero = b + cc} else {b.zero=b}
  if(c==0) {c.zero = c + cc} else {c.zero=c}
  if(d==0) {d.zero = d + cc} else {d.zero=d}

  nc.zero = c.zero+d.zero
  nt.zero = a.zero+b.zero

  rr.all = (a.all / nt.all) / (c.all / nc.all)
  rr.zero = (a.zero / nt.zero) / (c.zero / nc.zero)

  r.rr.df=data.frame(rra=rr.all,rrz=rr.zero)
  return(r.rr.df)
}


############################################################################################
#'  This Function to calculate RD and the point of change
#' @param cc.df - description
#' @param cc.value - description
#' @return A dataframe with
#'  \item{r.df}{Sequence}
#' @family Internal methods
#' @keywords internal
Calculate.cc.RD <-function(cc.df,cc.value)
{
  a=cc.df$a; b=cc.df$b;
  c=cc.df$c; d=cc.df$d;

  cc=cc.value
  nt=a+b
  nc=c+d

  a.all = a + cc
  b.all = b + cc
  c.all = c + cc
  d.all = d + cc

  nc.all = nc + (2 * cc)
  nt.all = nt + (2 * cc)

  ################ Calculat only for cells with zero

  if(a==0) {a.zero = a + cc} else {a.zero=a}
  if(b==0) {b.zero = b + cc} else {b.zero=b}
  if(c==0) {c.zero = c + cc} else {c.zero=c}
  if(d==0) {d.zero = d + cc} else {d.zero=d}

  nc.zero = c.zero+d.zero
  nt.zero = a.zero+b.zero

  rd.all = (a.all / nt.all) - (c.all / nc.all)
  rd.zero = (a.zero / nt.zero) - (c.zero / nc.zero)

  r.rd.df=data.frame(rda=rd.all,rdz=rd.zero)
  return(r.rd.df)
}

############################################################################################
#'  This Is between utility function
#' @param x - description
#' @param a - description
#' @param b - description
#' @return A dataframe with
#'  \item{r.df}{Sequence}
#' @family Internal methods
#' @keywords internal
is.between <- function(x, a, b) {
  small= min(a,b)
  big  = max(a,b)
  x > small & x < big
}

