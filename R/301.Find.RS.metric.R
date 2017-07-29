#'  Find the point of reversal of OR/RR/RD
#' @param Figure.df - description
#' @return A dataframe with
#'  \item{x1}{Sequence}
#' @family Internal methods
#' @keywords internal
Find.RS.metric<-function(Figure.df)
{

LHS.Flag.all.or=0
RHS.Flag.all.or=0

LHS.Flag.all.rr=0
RHS.Flag.all.rr=0

LHS.Flag.all.rd=0
RHS.Flag.all.rd=0

LHS.Flag.zero.or=0
RHS.Flag.zero.or=0

LHS.Flag.zero.rr=0
RHS.Flag.zero.rr=0

LHS.Flag.zero.rd=0
RHS.Flag.zero.rd=0

Figure.df$LHS.Flag.all.or=NULL
Figure.df$RHS.Flag.all.or=NULL

Figure.df$LHS.Flag.all.rr=NULL
Figure.df$RHS.Flag.all.rr=NULL

Figure.df$LHS.Flag.all.rd=NULL
Figure.df$RHS.Flag.all.rd=NULL

Figure.df$LHS.Flag.zero.or=NULL
Figure.df$RHS.Flag.zero.or=NULL

Figure.df$LHS.Flag.zero.rr=NULL
Figure.df$RHS.Flag.zero.rr=NULL

Figure.df$LHS.Flag.zero.rd=NULL
Figure.df$RHS.Flag.zero.rd=NULL

ibv.Rows=nrow(Figure.df)

for(i in 1:ibv.Rows){

  if(is.between(1,Figure.df$or.cc.all.tpm8[i],Figure.df$or.cc.all.tpm4[i])){LHS.Flag.all.or=1/100000000;RHS.Flag.all.or=1/10000}
  if(is.between(1,Figure.df$or.cc.all.tpm4[i],Figure.df$or.cc.all.001[i])){LHS.Flag.all.or=1/10000;RHS.Flag.all.or=0.0011}
  if(is.between(1,Figure.df$or.cc.all.001[i],Figure.df$or.cc.all.01[i])){LHS.Flag.all.or=0.0011;RHS.Flag.all.or=0.011}
  if(is.between(1,Figure.df$or.cc.all.01[i],Figure.df$or.cc.all.015[i])){LHS.Flag.all.or=0.011;RHS.Flag.all.or=0.0151}
  if(is.between(1,Figure.df$or.cc.all.015[i],Figure.df$or.cc.all.1[i])){LHS.Flag.all.or=0.0151;RHS.Flag.all.or=0.11}
  if(is.between(1,Figure.df$or.cc.all.1[i],Figure.df$or.cc.all.15[i])){LHS.Flag.all.or=0.11;RHS.Flag.all.or=0.151}
  if(is.between(1,Figure.df$or.cc.all.15[i],Figure.df$or.cc.all.2[i])){LHS.Flag.all.or=0.151;RHS.Flag.all.or=0.21}
  if(is.between(1,Figure.df$or.cc.all.2[i],Figure.df$or.cc.all.25[i])){LHS.Flag.all.or=0.21;RHS.Flag.all.or=0.251}
  if(is.between(1,Figure.df$or.cc.all.25[i],Figure.df$or.cc.all.3[i])){LHS.Flag.all.or=0.251;RHS.Flag.all.or=0.31}
  if(is.between(1,Figure.df$or.cc.all.3[i],Figure.df$or.cc.all.35[i])){LHS.Flag.all.or=0.31;RHS.Flag.all.or=0.351}
  if(is.between(1,Figure.df$or.cc.all.35[i],Figure.df$or.cc.all.4[i])){LHS.Flag.all.or=0.351;RHS.Flag.all.or=0.41}
  if(is.between(1,Figure.df$or.cc.all.4[i],Figure.df$or.cc.all.45[i])){LHS.Flag.all.or=0.41;RHS.Flag.all.or=0.451}
  if(is.between(1,Figure.df$or.cc.all.45[i],Figure.df$or.cc.all.5[i])){LHS.Flag.all.or=0.451;RHS.Flag.all.or=0.5}

  if(is.between(1,Figure.df$or.cc.zero.tpm8[i],Figure.df$or.cc.zero.tpm4[i])){LHS.Flag.zero.or=1/100000000;RHS.Flag.zero.or=1/10000}
  if(is.between(1,Figure.df$or.cc.zero.tpm4[i],Figure.df$or.cc.zero.001[i])){LHS.Flag.zero.or=1/10000;RHS.Flag.zero.or=0.0011}
  if(is.between(1,Figure.df$or.cc.zero.001[i],Figure.df$or.cc.zero.01[i])){LHS.Flag.zero.or=0.0011;RHS.Flag.zero.or=0.011}
  if(is.between(1,Figure.df$or.cc.zero.01[i],Figure.df$or.cc.zero.015[i])){LHS.Flag.zero.or=0.011;RHS.Flag.zero.or=0.0151}
  if(is.between(1,Figure.df$or.cc.zero.015[i],Figure.df$or.cc.zero.1[i])){LHS.Flag.zero.or=0.0151;RHS.Flag.zero.or=0.11}
  if(is.between(1,Figure.df$or.cc.zero.1[i],Figure.df$or.cc.zero.15[i])){LHS.Flag.zero.or=0.11;RHS.Flag.zero.or=0.151}
  if(is.between(1,Figure.df$or.cc.zero.15[i],Figure.df$or.cc.zero.2[i])){LHS.Flag.zero.or=0.151;RHS.Flag.zero.or=0.21}
  if(is.between(1,Figure.df$or.cc.zero.2[i],Figure.df$or.cc.zero.25[i])){LHS.Flag.zero.or=0.21;RHS.Flag.zero.or=0.251}
  if(is.between(1,Figure.df$or.cc.zero.25[i],Figure.df$or.cc.zero.3[i])){LHS.Flag.zero.or=0.251;RHS.Flag.zero.or=0.31}
  if(is.between(1,Figure.df$or.cc.zero.3[i],Figure.df$or.cc.zero.35[i])){LHS.Flag.zero.or=0.31;RHS.Flag.zero.or=0.351}
  if(is.between(1,Figure.df$or.cc.zero.35[i],Figure.df$or.cc.zero.4[i])){LHS.Flag.zero.or=0.351;RHS.Flag.zero.or=0.41}
  if(is.between(1,Figure.df$or.cc.zero.4[i],Figure.df$or.cc.zero.45[i])){LHS.Flag.zero.or=0.41;RHS.Flag.zero.or=0.451}
  if(is.between(1,Figure.df$or.cc.zero.45[i],Figure.df$or.cc.zero.5[i])){LHS.Flag.zero.or=0.451;RHS.Flag.zero.or=0.5}


  Figure.df$LHS.Flag.all.or[i]=LHS.Flag.all.or
  Figure.df$RHS.Flag.all.or[i]=RHS.Flag.all.or
  Figure.df$LHS.Flag.zero.or[i]=LHS.Flag.zero.or
  Figure.df$RHS.Flag.zero.or[i]=RHS.Flag.zero.or

  if(is.between(1,Figure.df$rr.cc.all.tpm8[i],Figure.df$rr.cc.all.tpm4[i])){LHS.Flag.all.rr=1/100000000;RHS.Flag.all.rr=1/10000}
  if(is.between(1,Figure.df$rr.cc.all.tpm4[i],Figure.df$rr.cc.all.001[i])){LHS.Flag.all.rr=1/10000;RHS.Flag.all.rr=0.0011}
  if(is.between(1,Figure.df$rr.cc.all.001[i],Figure.df$rr.cc.all.01[i])){LHS.Flag.all.rr=0.0011;RHS.Flag.all.rr=0.011}
  if(is.between(1,Figure.df$rr.cc.all.01[i],Figure.df$rr.cc.all.015[i])){LHS.Flag.all.rr=0.011;RHS.Flag.all.rr=0.0151}
  if(is.between(1,Figure.df$rr.cc.all.015[i],Figure.df$rr.cc.all.1[i])){LHS.Flag.all.rr=0.0151;RHS.Flag.all.rr=0.11}
  if(is.between(1,Figure.df$rr.cc.all.1[i],Figure.df$rr.cc.all.15[i])){LHS.Flag.all.rr=0.11;RHS.Flag.all.rr=0.151}
  if(is.between(1,Figure.df$rr.cc.all.15[i],Figure.df$rr.cc.all.2[i])){LHS.Flag.all.rr=0.151;RHS.Flag.all.rr=0.21}
  if(is.between(1,Figure.df$rr.cc.all.2[i],Figure.df$rr.cc.all.25[i])){LHS.Flag.all.rr=0.21;RHS.Flag.all.rr=0.251}
  if(is.between(1,Figure.df$rr.cc.all.25[i],Figure.df$rr.cc.all.3[i])){LHS.Flag.all.rr=0.251;RHS.Flag.all.rr=0.31}
  if(is.between(1,Figure.df$rr.cc.all.3[i],Figure.df$rr.cc.all.35[i])){LHS.Flag.all.rr=0.31;RHS.Flag.all.rr=0.351}
  if(is.between(1,Figure.df$rr.cc.all.35[i],Figure.df$rr.cc.all.4[i])){LHS.Flag.all.rr=0.351;RHS.Flag.all.rr=0.41}
  if(is.between(1,Figure.df$rr.cc.all.4[i],Figure.df$rr.cc.all.45[i])){LHS.Flag.all.rr=0.41;RHS.Flag.all.rr=0.451}
  if(is.between(1,Figure.df$rr.cc.all.45[i],Figure.df$rr.cc.all.5[i])){LHS.Flag.all.rr=0.451;RHS.Flag.all.rr=0.5}

  if(is.between(1,Figure.df$rr.cc.zero.tpm8[i],Figure.df$rr.cc.zero.tpm4[i])){LHS.Flag.zero.rr=1/100000000;RHS.Flag.zero.rr=1/10000}
  if(is.between(1,Figure.df$rr.cc.zero.tpm4[i],Figure.df$rr.cc.zero.001[i])){LHS.Flag.zero.rr=1/10000;RHS.Flag.zero.rr=0.0011}
  if(is.between(1,Figure.df$rr.cc.zero.001[i],Figure.df$rr.cc.zero.01[i])){LHS.Flag.zero.rr=0.0011;RHS.Flag.zero.rr=0.011}
  if(is.between(1,Figure.df$rr.cc.zero.01[i],Figure.df$rr.cc.zero.015[i])){LHS.Flag.zero.rr=0.011;RHS.Flag.zero.rr=0.0151}
  if(is.between(1,Figure.df$rr.cc.zero.015[i],Figure.df$rr.cc.zero.1[i])){LHS.Flag.zero.rr=0.0151;RHS.Flag.zero.rr=0.11}
  if(is.between(1,Figure.df$rr.cc.zero.1[i],Figure.df$rr.cc.zero.15[i])){LHS.Flag.zero.rr=0.11;RHS.Flag.zero.rr=0.151}
  if(is.between(1,Figure.df$rr.cc.zero.15[i],Figure.df$rr.cc.zero.2[i])){LHS.Flag.zero.rr=0.151;RHS.Flag.zero.rr=0.21}
  if(is.between(1,Figure.df$rr.cc.zero.2[i],Figure.df$rr.cc.zero.25[i])){LHS.Flag.zero.rr=0.21;RHS.Flag.zero.rr=0.251}
  if(is.between(1,Figure.df$rr.cc.zero.25[i],Figure.df$rr.cc.zero.3[i])){LHS.Flag.zero.rr=0.251;RHS.Flag.zero.rr=0.31}
  if(is.between(1,Figure.df$rr.cc.zero.3[i],Figure.df$rr.cc.zero.35[i])){LHS.Flag.zero.rr=0.31;RHS.Flag.zero.rr=0.351}
  if(is.between(1,Figure.df$rr.cc.zero.35[i],Figure.df$rr.cc.zero.4[i])){LHS.Flag.zero.rr=0.351;RHS.Flag.zero.rr=0.41}
  if(is.between(1,Figure.df$rr.cc.zero.4[i],Figure.df$rr.cc.zero.45[i])){LHS.Flag.zero.rr=0.41;RHS.Flag.zero.rr=0.451}
  if(is.between(1,Figure.df$rr.cc.zero.45[i],Figure.df$rr.cc.zero.5[i])){LHS.Flag.zero.rr=0.451;RHS.Flag.zero.rr=0.5}


  Figure.df$LHS.Flag.all.rr[i]=LHS.Flag.all.rr
  Figure.df$RHS.Flag.all.rr[i]=RHS.Flag.all.rr
  Figure.df$LHS.Flag.zero.rr[i]=LHS.Flag.zero.rr
  Figure.df$RHS.Flag.zero.rr[i]=RHS.Flag.zero.rr

  if(is.between(0,Figure.df$rd.cc.all.tpm8[i],Figure.df$rd.cc.all.tpm4[i])){LHS.Flag.all.rd=1/100000000;RHS.Flag.all.rd=1/10000}
  if(is.between(0,Figure.df$rd.cc.all.tpm4[i],Figure.df$rd.cc.all.001[i])){LHS.Flag.all.rd=1/10000;RHS.Flag.all.rd=0.0011}
  if(is.between(0,Figure.df$rd.cc.all.001[i],Figure.df$rd.cc.all.01[i])){LHS.Flag.all.rd=0.0011;RHS.Flag.all.rd=0.011}
  if(is.between(0,Figure.df$rd.cc.all.01[i],Figure.df$rd.cc.all.015[i])){LHS.Flag.all.rd=0.011;RHS.Flag.all.rd=0.0151}
  if(is.between(0,Figure.df$rd.cc.all.015[i],Figure.df$rd.cc.all.1[i])){LHS.Flag.all.rd=0.0151;RHS.Flag.all.rd=0.11}
  if(is.between(0,Figure.df$rd.cc.all.1[i],Figure.df$rd.cc.all.15[i])){LHS.Flag.all.rd=0.11;RHS.Flag.all.rd=0.151}
  if(is.between(0,Figure.df$rd.cc.all.15[i],Figure.df$rd.cc.all.2[i])){LHS.Flag.all.rd=0.151;RHS.Flag.all.rd=0.21}
  if(is.between(0,Figure.df$rd.cc.all.2[i],Figure.df$rd.cc.all.25[i])){LHS.Flag.all.rd=0.21;RHS.Flag.all.rd=0.251}
  if(is.between(0,Figure.df$rd.cc.all.25[i],Figure.df$rd.cc.all.3[i])){LHS.Flag.all.rd=0.251;RHS.Flag.all.rd=0.31}
  if(is.between(0,Figure.df$rd.cc.all.3[i],Figure.df$rd.cc.all.35[i])){LHS.Flag.all.rd=0.31;RHS.Flag.all.rd=0.351}
  if(is.between(0,Figure.df$rd.cc.all.35[i],Figure.df$rd.cc.all.4[i])){LHS.Flag.all.rd=0.351;RHS.Flag.all.rd=0.41}
  if(is.between(0,Figure.df$rd.cc.all.4[i],Figure.df$rd.cc.all.45[i])){LHS.Flag.all.rd=0.41;RHS.Flag.all.rd=0.451}
  if(is.between(0,Figure.df$rd.cc.all.45[i],Figure.df$rd.cc.all.5[i])){LHS.Flag.all.rd=0.451;RHS.Flag.all.rd=0.5}

  if(is.between(0,Figure.df$rd.cc.zero.tpm8[i],Figure.df$rd.cc.zero.tpm4[i])){LHS.Flag.zero.rd=1/100000000;RHS.Flag.zero.rd=1/10000}
  if(is.between(0,Figure.df$rd.cc.zero.tpm4[i],Figure.df$rd.cc.zero.001[i])){LHS.Flag.zero.rd=1/10000;RHS.Flag.zero.rd=0.0011}
  if(is.between(0,Figure.df$rd.cc.zero.001[i],Figure.df$rd.cc.zero.01[i])){LHS.Flag.zero.rd=0.0011;RHS.Flag.zero.rd=0.011}
  if(is.between(0,Figure.df$rd.cc.zero.01[i],Figure.df$rd.cc.zero.015[i])){LHS.Flag.zero.rd=0.011;RHS.Flag.zero.rd=0.0151}
  if(is.between(0,Figure.df$rd.cc.zero.015[i],Figure.df$rd.cc.zero.1[i])){LHS.Flag.zero.rd=0.0151;RHS.Flag.zero.rd=0.11}
  if(is.between(0,Figure.df$rd.cc.zero.1[i],Figure.df$rd.cc.zero.15[i])){LHS.Flag.zero.rd=0.11;RHS.Flag.zero.rd=0.151}
  if(is.between(0,Figure.df$rd.cc.zero.15[i],Figure.df$rd.cc.zero.2[i])){LHS.Flag.zero.rd=0.151;RHS.Flag.zero.rd=0.21}
  if(is.between(0,Figure.df$rd.cc.zero.2[i],Figure.df$rd.cc.zero.25[i])){LHS.Flag.zero.rd=0.21;RHS.Flag.zero.rd=0.251}
  if(is.between(0,Figure.df$rd.cc.zero.25[i],Figure.df$rd.cc.zero.3[i])){LHS.Flag.zero.rd=0.251;RHS.Flag.zero.rd=0.31}
  if(is.between(0,Figure.df$rd.cc.zero.3[i],Figure.df$rd.cc.zero.35[i])){LHS.Flag.zero.rd=0.31;RHS.Flag.zero.rd=0.351}
  if(is.between(0,Figure.df$rd.cc.zero.35[i],Figure.df$rd.cc.zero.4[i])){LHS.Flag.zero.rd=0.351;RHS.Flag.zero.rd=0.41}
  if(is.between(0,Figure.df$rd.cc.zero.4[i],Figure.df$rd.cc.zero.45[i])){LHS.Flag.zero.rd=0.41;RHS.Flag.zero.rd=0.451}
  if(is.between(0,Figure.df$rd.cc.zero.45[i],Figure.df$rd.cc.zero.5[i])){LHS.Flag.zero.rd=0.451;RHS.Flag.zero.rd=0.5}


  Figure.df$LHS.Flag.all.rd[i]=LHS.Flag.all.rd
  Figure.df$RHS.Flag.all.rd[i]=RHS.Flag.all.rd
  Figure.df$LHS.Flag.zero.rd[i]=LHS.Flag.zero.rd
  Figure.df$RHS.Flag.zero.rd[i]=RHS.Flag.zero.rd

} # End of for loop extracting ibv range

ibv.ident.Rows=nrow(Figure.df)

for(kk in 1:ibv.ident.Rows){
  Bet.all.or=Figure.df$RHS.Flag.all.or[kk]-Figure.df$LHS.Flag.all.or[kk]
  if(Bet.all.or>0){
    BAll.seq.or= Find.sequence(Figure.df$RHS.Flag.all.or[kk],Figure.df$LHS.Flag.all.or[kk])
    BAll.values.or= pass2.In.Each.Row.new.cc.or(Figure.df[kk,],BAll.seq.or)
  }
  if(Bet.all.or==0){
    BAll.values.or= data.frame( RS.all.or=0)
  }

  Bet.zero.or=Figure.df$RHS.Flag.zero.or[kk]-Figure.df$LHS.Flag.zero.or[kk]
  if(Bet.zero.or>0){
    BAll.seq.or= Find.sequence(Figure.df$RHS.Flag.zero.or[kk],Figure.df$LHS.Flag.zero.or[kk])
    BAll.values.zor= pass2.In.Each.Row.new.cc.or(Figure.df[kk,],BAll.seq.or)
  }
  if(Bet.zero.or==0){
    BAll.values.zor= data.frame( RS.zero.or=0)
  }

  ##################### RR #################################################################################
  Bet.all.rr=Figure.df$RHS.Flag.all.rr[kk]-Figure.df$LHS.Flag.all.rr[kk]
  if(Bet.all.rr>0){
    BAll.seq.rr= Find.sequence(Figure.df$RHS.Flag.all.rr[kk],Figure.df$LHS.Flag.all.rr[kk])
    BAll.values.rr= pass2.In.Each.Row.new.cc.rr(Figure.df[kk,],BAll.seq.rr)
  }
  if(Bet.all.rr==0){
    BAll.values.rr= data.frame( RS.all.rr=0)
  }

  Bet.zero.rr=Figure.df$RHS.Flag.zero.rr[kk]-Figure.df$LHS.Flag.zero.rr[kk]
  if(Bet.zero.rr>0){
    BAll.seq.rr= Find.sequence(Figure.df$RHS.Flag.zero.rr[kk],Figure.df$LHS.Flag.zero.rr[kk])
    BAll.values.zrr= pass2.In.Each.Row.new.cc.rr(Figure.df[kk,],BAll.seq.rr)
  }
  if(Bet.zero.rr==0){
    BAll.values.zrr= data.frame( RS.zero.rr=0)
  }

  #################### RD #################################################################################
  Bet.all.rd=Figure.df$RHS.Flag.all.rd[kk]-Figure.df$LHS.Flag.all.rd[kk]
  if(Bet.all.rd>0){
    BAll.seq.rd= Find.sequence(Figure.df$RHS.Flag.all.rd[kk],Figure.df$LHS.Flag.all.rd[kk])
    BAll.values.rd= pass2.In.Each.Row.new.cc.rd(Figure.df[kk,],BAll.seq.rd)
  }
  if(Bet.all.rd==0){
    BAll.values.rd= data.frame(RS.all.rd=0)
  }

  Bet.zero.rd=Figure.df$RHS.Flag.zero.rd[kk]-Figure.df$LHS.Flag.zero.rd[kk]
  if(Bet.zero.rd>0){
    BAll.seq.rd= Find.sequence(Figure.df$RHS.Flag.zero.rd[kk],Figure.df$LHS.Flag.zero.rd[kk])
    BAll.values.zrd= pass2.In.Each.Row.new.cc.rd(Figure.df[kk,],BAll.seq.rd)
  }
  if(Bet.zero.rd==0){
    BAll.values.zrd= data.frame(RS.zero.rd=0)
  }


  tvar=cbind(
    BAll.values.or$RS.all.or,
    BAll.values.rr$RS.all.rr,
    BAll.values.rd$RS.all.rd,
    BAll.values.zor$RS.zero.or,
    BAll.values.zrr$RS.zero.rr,
    BAll.values.zrd$RS.zero.rd)

  tvar.df=as.data.frame(tvar)
  names(tvar.df)=   c( "RS.all.or",
                    "RS.all.rr",
                    "RS.all.rd",
                    "RS.zero.zor"      ,
                    "RS.zero.zrr",
                    "RS.zero.zrd" )

  rcc=cbind(Figure.df[kk,],RS.all.or=tvar.df$RS.all.or, RS.all.rr=tvar.df$RS.all.rr, RS.all.rd=tvar.df$RS.all.rd,
            RS.zero.zor=tvar.df$RS.zero.zor,RS.zero.zrr=tvar.df$RS.zero.zrr,RS.zero.zrd=tvar.df$RS.zero.zrd )

}
return(rcc)
} # End of function Find.RS.Metric

######################################################################################################
#'  Find the sequence for the given  range
#' @param RHS/value - description
#' @param LHS.value - description
#' @return A dataframe with
#'  \item{x1}{Sequence}
#' @family Internal methods
#' @keywords internal

Find.sequence <- function(RHS.value,LHS.value){
  small=min(LHS.value, RHS.value)
  #if(small>0.01){small=small+0.001}
  big=max(LHS.value, RHS.value)
  Ret.seq=seq(small, big, length.out = 10)
  return(Ret.seq)
}


######################################################################################################
#'  Find the OR,RR and RD for each row with new CC range
#' @param odf.row - description
#' @param cc.vals - description
#' @return A dataframe with
#'  \item{x1}{Sequence}
#' @family Internal methods
#' @keywords internal
pass2.In.Each.Row.new.cc.or <-function(odf.row,cc.vals)
{
  iterations = length(cc.vals)

  previous.ora=data.frame(ora=0,orz=0)
  previous.cc = 0
  RS.all.or=0
  RS.zero.or=0

  for(i in 1:iterations){

    cc.value = cc.vals[i]

    ora.results=Calculate.cc.OR(odf.row, cc.value )
    if(i>1){
      if(previous.ora$ora==1){RS.all.or=previous.cc}
      if(is.between(1,previous.ora$ora,ora.results$ora)){RS.all.or=previous.cc}
      if(previous.ora$orz==1){RS.zero.or=previous.cc}
      if(is.between(1,previous.ora$orz,ora.results$orz)){RS.zero.or=previous.cc}
    }
    previous.ora=ora.results
    previous.cc = cc.value

  }

  RS.metric.or.df = data.frame(
    RS.all.or=RS.all.or,
    RS.zero.or=RS.zero.or)
  return(RS.metric.or.df)
}


######################################################################################################
#'  Find the OR,RR and RD for each row with new CC range
#' @param odf.row - description
#' @param cc.vals - description
#' @return A dataframe with
#'  \item{x1}{Sequence}
#' @family Internal methods
#' @keywords internal
pass2.In.Each.Row.new.cc.rr <-function(odf.row,cc.vals)
{

  iterations = length(cc.vals)

  previous.rra=data.frame(rra=0,rrz=0)
  previous.cc=0
  RS.all.rr=0
  RS.zero.rr=0

  for(i in 1:iterations){

    cc.value = cc.vals[i]

    rra.results=Calculate.cc.RR(odf.row, cc.value )
    if(i>1){
      if(previous.rra$rra==1){RS.all.rr=previous.cc}
      if(previous.rra$rrz==1){RS.zero.rr=previous.cc}
      if(is.between(1,previous.rra$rra,rra.results$rra)){RS.all.rr=previous.cc}
      if(is.between(1,previous.rra$rrz,rra.results$rrz)){RS.zero.rr=previous.cc}
    }
    previous.rra=rra.results
    previous.cc=cc.value

  }

  RS.metric.rr.df = data.frame(
    RS.all.rr=RS.all.rr,
    RS.zero.rr=RS.zero.rr)

  return(RS.metric.rr.df)
}
############################################################################################


######################################################################################################
#'  Find the OR,RR and RD for each row with new CC range
#' @param odf.row - description
#' @param cc.vals - description
#' @return A dataframe with
#'  \item{x1}{Sequence}
#' @family Internal methods
#' @keywords internal
pass2.In.Each.Row.new.cc.rd <-function(odf.row,cc.vals)
{

  iterations = length(cc.vals)

  previous.rda=data.frame(rda=0,rdz=0)
  previous.cc=0
  RS.all.rd=0
  RS.zero.rd=0

  for(i in 1:iterations){

    cc.value = cc.vals[i]
    rda.results=Calculate.cc.RD(odf.row, cc.value )

    if(i>1){
      if(previous.rda$rda==0){RS.all.rd=previous.cc}
      if(previous.rda$rdz==0){RS.zero.rd=previous.cc}
      if(is.between(0,previous.rda$rda,rda.results$rda)){RS.all.rd=previous.cc}
      if(is.between(0,previous.rda$rdz,rda.results$rdz)){RS.zero.rd=previous.cc}
    }
    previous.rda=rda.results
    previous.cc=cc.value
  }

  RS.metric.rd.df = data.frame(
    RS.all.rd=RS.all.rd,
    RS.zero.rd=RS.zero.rd)
  return(RS.metric.rd.df)
}

############################################################################################



