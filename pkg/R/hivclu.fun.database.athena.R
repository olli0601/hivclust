#' @export
HIVC.db.locktime	<-  as.Date("30/03/2013", format="%d/%m/%Y")

######################################################################################
hivc.db.resetNegTbyPoslRNA_T1<- function(df.all,verbose=1)
{
	if(verbose)	cat(paste("\nreset NegT to NA for the following seq"))
	if(verbose) print(subset( df.all,!is.na(PoslRNA_T1) & PoslRNA_T1<NegT ))
	tmp		<- which(df.all[, NegT_Acc=="Yes" & !is.na(PoslRNA_T1) & PoslRNA_T1<NegT])
	if(verbose) cat(paste("\nnumber of seq with NegT_Acc=='Yes' & !is.na(PoslRNA_T1) & NegT > time of VL test, n=",length(tmp)))	
	set(df.all, tmp, "NegT", df.all[tmp, PoslRNA_T1-6*30])
	set(df.all, tmp, "NegT_Acc", 'No')
	tmp		<- which(df.all[, NegT_Acc=="No" & !is.na(PoslRNA_T1) & PoslRNA_T1<NegT])
	if(verbose) cat(paste("\nnumber of seq with NegT_Acc=='No' & !is.na(PoslRNA_T1) & NegT > time of VL test, n=",length(tmp)))			
	set(df.all, tmp, "NegT", df.all[tmp, PoslRNA_T1-6*30])
	set(df.all, tmp, "NegT_Acc", 'No')
	df.all
}
######################################################################################
hivc.db.resetNegTbyAnyPosT<- function(df.all,verbose=1)
{	
	tmp		<- df.all[, which(NegT_Acc=="Yes" & NegT>AnyPos_T1)]
	cat(paste("\nnumber of seq with NegT_Acc=='Yes' & NegT> pos test or pos seq, n=",length(tmp)))
	print(as.character(unique(df.all[tmp,Patient])))
	set(df.all, tmp, "NegT", NA_real_)
	set(df.all, tmp, "NegT_Acc", NA_character_)
	tmp		<- df.all[, which(NegT_Acc=="No" & NegT>AnyPos_T1)]
	cat(paste("\nnumber of seq with NegT_Acc=='No' & NegT> pos test or pos seq, n=",length(tmp)))
	print(as.character(unique(df.all[tmp,Patient])))
	set(df.all, tmp, "NegT", NA)
	set(df.all, tmp, "NegT_Acc", NA_character_)
	df.all
}		
######################################################################################
hivc.db.resetCD4byNegT<- function(df.cross, with.NegT_Acc.No=0, verbose=1)
{
	tmp			<- df.cross[, which(NegT_Acc=="Yes" & 	NegT>PosCD4)]
	if(verbose)		cat(paste("\nnumber of entries with NegT_Acc=='Yes' & 	NegT>PosCD4, n=", length(tmp),"SETTING CD4 to NA"))
	set(df.cross, tmp, "PosCD4", NA)
	set(df.cross, tmp, "CD4", NA)
	df.cross	<- subset(df.cross,!is.na(CD4))
	if(!with.NegT_Acc.No)
		return(df.cross)
	tmp			<- df.cross[,which(NegT_Acc=="No" & 	NegT>PosCD4)]
	if(verbose)		cat(paste("\nnumber of entries with NegT_Acc=='No' & 	NegT>PosCD4, n=", length(tmp),"SETTING CD4 to NA"))
	set(df.cross, tmp, "PosCD4", NA)
	set(df.cross, tmp, "CD4", NA)
	df.cross	<- subset(df.cross,!is.na(CD4))
	df.cross
}	
######################################################################################
hivc.db.getTrIMo<- function(df.cross, verbose=1)
{
	ans		<- df.cross[, 	{ 
				x		<- data.table( PosSeqT, StartTime, StopTime, TrI )
				x		<- subset(x, !is.na(StartTime))		#discard first time period if unknown StartTime					
				if(!nrow(x) || is.na(x[1,PosSeqT]))
				{
					TrImo_bTS		<- TrImo_aTS	<- NA_real_
				}
				else if(x[1,PosSeqT]<=x[1,StartTime])
				{
					TrImo_bTS		<- 0
					TrImo_aTS		<- sum(as.numeric(subset(x, TrI=="Yes")[,difftime(StopTime,StartTime,units="days")/30]))
				}
				else
				{
					z				<- which( as.numeric( x[,difftime(PosSeqT,StartTime,units="days")] )>0 )	#all rows in which PosSeqT>=StartTime						
					z				<- z[length(z)]
					z2				<- seq.int(1,z)
					xb				<- x[z2,]
					set(xb,z,"StopTime",x[z,PosSeqT])
					z2				<- seq.int(z,nrow(x))
					xa				<- x[z2,]
					set(xa,1L,"StartTime",x[z,PosSeqT])		
					z				<- as.numeric(xb[,difftime(StopTime,StartTime,units="days")/30])
					z2				<- as.numeric(xa[,difftime(StopTime,StartTime,units="days")/30])
					TrImo_bTS		<- sum(z[ which(xb[, TrI=="Yes"]) ])
					TrImo_aTS		<- sum(z2[ which(xa[, TrI=="Yes"]) ])
				}					
				list( 	TrImo_bTS	= TrImo_bTS, 
						TrImo_aTS	= TrImo_aTS			)																					
			}, by=idx]
	set(ans,NULL,"TrImo_bTS",round(ans[,TrImo_bTS],d=1))
	set(ans,NULL,"TrImo_aTS",round(ans[,TrImo_aTS],d=1))
	if(verbose)	cat(paste("\nnumber of entries, n=",nrow(ans)))
	ans		
}
######################################################################################
hivc.db.getlRNA.T1andTS<- function(df.cross, lRNA.bTS.quantile= 0.75, lRNA.aTS.quantile= 0.25, lRNAi.min= log10(1e4), db.locktime= HIVC.db.locktime, verbose=1)
{						
	if(verbose)	cat(paste("\nlRNA.aTS.quantile is",lRNA.aTS.quantile))
	if(verbose)	cat(paste("\nlRNA.bTS.quantile is",lRNA.bTS.quantile))
	if(verbose)	cat(paste("\nlRNAi.min is",lRNAi.min))
	if(verbose)	cat(paste("\nnumber of entries in cross product is, n=",nrow(df.cross)))
	df.cross[, lRNAi:=lRNA>=lRNAi.min]			
	ans		<-	df.cross[, 	{
				z	<- which.min(PosRNA)
				if(!is.na(PosSeqT[1]))
				{
					z2			<- which.min(abs(difftime(PosSeqT, PosRNA, units="weeks")))
					PoslRNA_TS	<- PosRNA[z2]
					lRNA_TS		<- round( mean( lRNA[seq.int(z2-1,z2+1)], na.rm=1 ), d=1 )				
					lRNA_bTS	<- round( quantile(lRNA[seq_len(z2+1)], probs = lRNA.bTS.quantile, na.rm=T, names = F), d=1 )
					lRNA_aTS	<- round( quantile(lRNA[seq.int(z2-1,length(lRNA))], probs = lRNA.aTS.quantile, na.rm=T, names = F), d=1 )							
					lRNAit		<- as.numeric( difftime(c(PosRNA[-1],db.locktime), PosRNA,units="days") )
					#print(data.table(Patient, FASTASampleCode, AnyPos_T1, PosSeqT, PosRNA, lRNA, lRNAi, lRNA.hb4tr_LT, lRNA.early)); print(lRNAit); print(z2); print(lRNAit[ seq.int(z2,length(lRNA)) ]); print(lRNAit[ seq_len(z2) ][  lRNAi[ seq_len(z2) ]  ])
					lRNAi_bTS	<- sum( lRNAit[ seq_len(z2) ][  lRNAi[ seq_len(z2) ]  ] ) / sum( lRNAit[ seq_len(z2) ] ) 							
					lRNAi_aTS	<- sum( lRNAit[seq.int(z2,length(lRNAit))][  lRNAi[ seq.int(z2,length(lRNAi)) ]	] ) /	sum( lRNAit[seq.int(z2,length(lRNAit))] )
				}
				else
				{
					PoslRNA_TS	<- lRNAi_bTS<- lRNAi_aTS<- as.Date(NA)
					lRNA_TS 	<- lRNA_bTS <- lRNA_aTS <- NA_real_
				}
				list(	PoslRNA_T1		= PosRNA[z], 
						lRNA_T1			= lRNA[z], 
						PoslRNA_TS		= PoslRNA_TS, 
						lRNA_TS			= lRNA_TS, 
						lRNA_bTS		= lRNA_bTS, 
						lRNA_aTS		= lRNA_aTS,
						lRNAi_bTS		= lRNAi_bTS,
						lRNAi_aTS		= lRNAi_aTS,
						lRNA.early		= lRNA.early[1]				) 							 	
			},by=idx]
	if(verbose)	cat(paste("\nnumber of seq with PosCD4_T1 CD4_T1  PosCD4_TS CD4_TS is n=",nrow(ans)))
	df.cross[,lRNAi:=NULL]
	ans	
}
######################################################################################
hivc.db.getcoverage<- function(df)
{
	ans	<- paste( c("df[,list(", paste( 	sapply(colnames(df), function(x) paste(x,"= length(which(!is.na(",x,")))",sep='')), collapse=",", sep='' ),")]") , collapse='',sep='')
	ans	<- eval(parse(text=ans))	/ nrow(df)
	ans
}
######################################################################################
hivc.db.getCD4.T1andTS<- function(df.cross, verbose=1, CD4.HIVNeg.min= 500, CD4.bTS.quantile= 0.75, CD4.aTS.quantile= 0.25)
{
	if(verbose)	cat(paste("\nCD4.HIVNeg.min is",CD4.HIVNeg.min))
	if(verbose)	cat(paste("\nCD4.aTS.quantile is",CD4.aTS.quantile))
	if(verbose)	cat(paste("\nCD4.bTS.quantile is",CD4.bTS.quantile))
	if(verbose)	cat(paste("\nnumber of entries in cross product is, n=",nrow(df.cross)))
	tmp	<- which( df.cross[, AnyPos_T1>PosCD4 & CD4>CD4.HIVNeg.min] )
	if(verbose)	cat(paste("\nnumber of patients with AnyPos_T1>PosCD4 & CD4>CD4.HIVNeg.min, n=",length(unique(df.cross[tmp,Patient])),"SETTING PosCD4 to NA since Patients could be healthy for corresponding CD4"))
	set(df.cross, tmp, "CD4", NA)		
	df.cross		<- subset(df.cross, !is.na(CD4))
	if(verbose)	cat(paste("\nnumber of entries in cross product is, n=",nrow(df.cross)))
	
	ans		<-	df.cross[, 	{
				z	<- which.min(PosCD4)
				if(!is.na(PosSeqT[1]))
				{
					z2			<- which.min(abs(difftime(PosSeqT, PosCD4, units="weeks")))
					PosCD4_TS	<- PosCD4[z2]
					CD4_TS		<- round( mean( CD4[seq.int(z2-1,z2+1)], na.rm=1 ), d=0 )				
					CD4_bTS		<- round( quantile(CD4[seq_len(z2+1)], probs = CD4.bTS.quantile, na.rm=T, names = F), d=0 )
					CD4_aTS		<- round( quantile(CD4[seq.int(z2-1,length(CD4))], probs = CD4.aTS.quantile, na.rm=T, names = F), d=0 )
				}
				else
				{
					PosCD4_TS	<- as.Date(NA)
					CD4_TS <- CD4_bTS <- CD4_aTS <- NA_real_
				}
				list(PosCD4_T1=PosCD4[z], CD4_T1=CD4[z], PosCD4_TS=PosCD4_TS, CD4_TS=CD4_TS, CD4_bTS=CD4_bTS, CD4_aTS=CD4_aTS	 ) 	
			},by=idx]
	
	if(verbose)	cat(paste("\nnumber of seq with PosCD4_T1 CD4_T1  PosCD4_TS CD4_TS is n=",nrow(ans)))
	ans
}
######################################################################################
hivc.db.reset.inaccuratePosT<- function(df, nacc.dy.dy= 30, nacc.mody.mo= 11, nacc.mody.dy= 31, verbose=1)
{	
	nacc.dy		<- which( df[, PosT_Acc=="No" & !is.na(PosT) & as.POSIXlt(PosT)$mday==15] )
	nacc.mody	<- which( df[, PosT_Acc=="No" & !is.na(PosT) & as.POSIXlt(PosT)$mon==6 & as.POSIXlt(PosT)$mday==1] )
	if(verbose) cat(paste("\nnumber of uncertain PosT, day only, n=",length(nacc.dy)))
	if(verbose) cat(paste("\nnumber of uncertain PosT, month & day, n=",length(nacc.mody)))
	# reset nacc.dy
	tmp							<- as.POSIXlt(df[nacc.dy,PosT] )
	tmp$mday					<- nacc.dy.dy
	set(df, nacc.dy, "PosT", as.Date(tmp))
	# reset nacc.mody
	tmp							<- as.POSIXlt(df[nacc.mody,PosT] )
	tmp$mday					<- nacc.mody.dy
	tmp$mon						<- nacc.mody.mo
	set(df, nacc.mody, "PosT", as.Date(tmp))
	df
}
######################################################################################
hivc.db.reset.inaccurateNegT<- function(df, nacc.dy.dy= 1, nacc.mody.mo= 0, nacc.mody.dy= 1, verbose=1)
{	
	serocon.nacc.dy		<- which( df[, NegT_Acc=="No" & !is.na(NegT) & as.POSIXlt(NegT)$mday==15] )
	serocon.nacc.mody	<- which( df[, NegT_Acc=="No" & !is.na(NegT) & as.POSIXlt(NegT)$mon==6 & as.POSIXlt(NegT)$mday==1] )
	if(verbose) cat(paste("\nnumber of uncertain NegT, day only, n=",length(serocon.nacc.dy)))
	if(verbose) cat(paste("\nnumber of uncertain NegT, month & day, n=",length(serocon.nacc.mody)))
	# reset serocon.nacc.dy
	tmp							<- as.POSIXlt(df[serocon.nacc.dy,NegT] )
	tmp$mday					<- nacc.dy.dy
	set(df, serocon.nacc.dy, "NegT", as.Date(tmp))
	#set(df, serocon.nacc.dy, "NegT_Acc", "Yes")
	# reset serocon.nacc.mody
	tmp							<- as.POSIXlt(df[serocon.nacc.mody,NegT] )
	tmp$mday					<- nacc.mody.dy
	tmp$mon						<- nacc.mody.mo
	set(df, serocon.nacc.mody, "NegT", as.Date(tmp))
	#set(df, serocon.nacc.mody, "NegT_Acc", "Yes")
	df
}
######################################################################################
hivc.db.reset.PosT1byCD4T1<- function(df.all, verbose=1)
{	
	tmp	<- df.all[, which(PosCD4_T1 < AnyPos_T1 & as.POSIXlt(AnyPos_T1)$mon==11 & as.POSIXlt(AnyPos_T1)$mday==31)] 
	if(verbose) cat(paste("\nnumber of PosCD4_T1 < AnyPos_T1 with AnyPos_T1==XX-12-31, reset, n=",length(tmp)))
	set(df.all, tmp, 'AnyPos_T1', df.all[tmp, PosCD4_T1])
	
	tmp	<- df.all[, which(PosCD4_T1 < AnyPos_T1 & as.POSIXlt(AnyPos_T1)$mon==10 & as.POSIXlt(AnyPos_T1)$mday==11)] 
	if(verbose) cat(paste("\nnumber of PosCD4_T1 < AnyPos_T1 with AnyPos_T1==XX-11-11, reset, n=",length(tmp)))
	set(df.all, tmp, 'AnyPos_T1', df.all[tmp, PosCD4_T1])
	
	tmp	<- df.all[, which(PosCD4_T1 < AnyPos_T1 & as.POSIXlt(AnyPos_T1)$mon==6 & as.POSIXlt(AnyPos_T1)$mday==15)] 
	if(verbose) cat(paste("\nnumber of PosCD4_T1 < AnyPos_T1 with AnyPos_T1==XX-07-15, reset, n=",length(tmp)))
	set(df.all, tmp, 'AnyPos_T1', df.all[tmp, PosCD4_T1])
		
	#tmp	<- subset(df.all, PosCD4_T1 < AnyPos_T1, c(FASTASampleCode, Patient, AnyPos_T1, PosSeqT, PosT, PosT_Acc, PoslRNA_T1, lRNA_T1,  PosCD4_T1, CD4_T1))
	#tmp[, diff:= difftime(PosCD4_T1, AnyPos_T1, units='days')]
	
	#rest unclear - reset
	tmp	<- df.all[, which(PosCD4_T1 < AnyPos_T1)]
	if(verbose) cat(paste("\nnumber of PosCD4_T1 < AnyPos_T1 with other patterns, reset, n=",length(tmp)))
	set(df.all, tmp, 'AnyPos_T1', df.all[tmp, PosCD4_T1])
	
	df.all
}
######################################################################################
hivc.db.getplot.newdiagnosesbyCD4<- function(df, plot.file=NULL, plot.file.p=NULL, plot.ylab=NULL, verbose=1)
{
	set(df, NULL, "AnyPos_T1", hivc.db.Date2numeric(df[,AnyPos_T1]))
	df[,AnyPos_yr:=	floor(df[,AnyPos_T1])]
	df[,CD4_T1bin:= my.aggregate(df[,CD4_T1], c(-Inf, 200, 350, 500, Inf))]
	#
	t.newdiag			<- table( df[,AnyPos_yr,CD4_T1bin] )
	rownames(t.newdiag)	<- c("<200","200-349","350-499",">=500")	
	if(!is.null(plot.file))
	{
		cols				<- brewer.pal(nrow(t.newdiag),"Set1")
		if(verbose)	cat(paste("\nplot file prop to ",plot.file))
		pdf(file=plot.file, width=5, height=4)
		par(mar=c(4,6,0.5,0.5))
		my.barplot.table(t.newdiag, 0.4, "year", plot.ylab, cols, x.baseline=0, xax=as.numeric(colnames(t.newdiag)), legend.loc="topleft")
		dev.off()
	}
	#
	p.newdiag			<- t.newdiag / matrix(rep(apply(t.newdiag,2,sum), each=nrow(t.newdiag)), nrow(t.newdiag), ncol(t.newdiag))
	if(!is.null(plot.file.p))
	{
		if(verbose)	cat(paste("\nplot file prop to ",plot.file.p))
		pdf(file=plot.file.p, width=5, height=4)
		par(mar=c(4,6,0.5,0.5))
		cols				<- brewer.pal(nrow(p.newdiag),"Set1")
		my.barplot.table(p.newdiag, 0.4, "year", paste("%",plot.ylab), cols, x.baseline=0, xax=as.numeric(colnames(p.newdiag)), legend.loc=NULL)
		dev.off()
	}
	list( t.newdiag=t.newdiag, p.newdiag=p.newdiag)
}
######################################################################################
hivc.db.getplot.livingbyCD4<- function(df, df.immu, plot.file, plot.file.p, plot.ylab, db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)
{
	#
	set(df, NULL, "AnyPos_T1", hivc.db.Date2numeric(df[,AnyPos_T1]))
	set(df, NULL, "DateDied", hivc.db.Date2numeric(df[,DateDied]))
	set(df, NULL, "DateLastContact", hivc.db.Date2numeric(df[,DateLastContact]))
	set(df, NULL, "AnyT_T1", hivc.db.Date2numeric(df[,AnyT_T1]))
	df[,DateEnd:= DateDied]
	#compute DateEnd, either death or lost contact
	tmp				<- which( df[, (DateDied-DateLastContact)>db.diff.lastcontact2died] )
	if(verbose)	cat(paste("\n(DateDied-DateLastContact)>db.diff.lastcontact2died for n=",length(tmp)))
	set(df, tmp, "DateEnd", df[tmp,DateLastContact])
	tmp				<- which(df[, is.na(DateDied) & (db.endtime-DateLastContact)>db.diff.lastcontact2now])
	if(verbose)	cat(paste("\n(db.endtime-DateLastContact)>db.diff.lastcontact2now for n=",length(tmp)))
	set(df, tmp, "DateEnd", df[tmp,DateLastContact])
	tmp				<- which(is.na(df[,DateLastContact]))
	if(verbose)	cat(paste("\nnumber patients with missing DateLastContact n=",length(tmp)))
	tmp				<- which(is.na(df[,DateEnd]))
	if(verbose)	cat(paste("\nnumber patients still alive n=",length(tmp)))
	set(df, tmp,"DateEnd",db.endtime)
	df[, AnyPos_yr:= floor(AnyPos_T1)]
	df[, DateEnd_yr:= floor(DateEnd)]
	#
	set(df, which(is.na(df[,AnyT_T1])), "AnyT_T1", max(df[,DateEnd_yr])+1)
	#
	df				<- subset(df, select=c(Patient, AnyPos_T1, CD4_T1, isAcute, AnyT_T1, DateEnd, AnyPos_yr, DateEnd_yr))	
	#compute patients alive in year Db_yr 
	tmp				<- seq.int(min(df[, AnyPos_yr]),max(df[, DateEnd_yr]))
	tmp				<- lapply(tmp, function(x)
			{				
				tmp	<- subset(df, AnyPos_yr<=x & DateEnd_yr>=x)
				cbind(tmp,data.table(Db_yr=rep(x,nrow(tmp))))
			})
	df				<- rbindlist(tmp)
	setkey(df, Patient)	
	#compute proportions of patients in recent, undiagnosed by CD4 count, treated per year	
	df	<- df[,	{
				x	<- data.table(Patient, AnyPos_T1, isAcute,  AnyT_T1, Db_yr, DateEnd_yr)
				#x	<- subset(df, Patient=="M10833", select=c(Patient, AnyPos_T1, isAcute,  AnyT_T1, Db_yr, DateEnd_yr))
				z	<- merge(x, df.immu, by="Patient")
				z	<- subset(z, floor(PosCD4)==Db_yr)
				tmp	<- setdiff( seq.int(x[,floor(AnyPos_T1)[1]], x[,DateEnd_yr[1]]), unique(z[,Db_yr]) )	#years for which no CD4 count available
				if(length(tmp))
					x	<- rbind(z,data.table(Patient=x[,Patient[1]], AnyPos_T1= x[,AnyPos_T1[1]], isAcute=x[,isAcute[1]], AnyT_T1=x[,AnyT_T1[1]], Db_yr=tmp, DateEnd_yr=x[,DateEnd_yr[1]],   PosCD4= tmp+.5, CD4=NA)) 
				else
					x	<- z
				
				tmp<- x[, 	{									
							tmp		<- min(Db_yr[1]+1,AnyT_T1[1])-max(Db_yr[1],AnyPos_T1[1])
							Acute_p	<- ifelse(!is.na(isAcute[1]) && isAcute[1]=="Yes" && floor(AnyPos_T1[1])==Db_yr, tmp, 0)
							AnyT_p	<- min(1,max(0,Db_yr[1]+1-AnyT_T1[1]))
							CD4_p	<- ifelse((is.na(isAcute[1]) || isAcute[1]!="Yes") && floor(AnyT_T1[1])>=Db_yr, tmp, 0)
							#print(Acute_p); print(AnyT_p); print(CD4_p)
							tmp2	<- which(PosCD4<=AnyT_T1)
							CD4_med	<- ifelse(length(tmp2),median(CD4[tmp2]), NA_real_)
							#print(CD4_yr)
							list(Patient=Patient[1], Acute_p=Acute_p, CD4_med=CD4_med, AnyT_p=AnyT_p, NotYetT_p=CD4_p )
						},by=Db_yr]
				tmp
			},by=Patient]
	df[,CD4_medbin:= my.aggregate(df[,CD4_med], c(-Inf, 200, 350, 500, Inf))]	
	set(df, which(is.na(df[,CD4_medbin])), "CD4_medbin", "unknown")
	#take sum for each stratification of interest
	t.living	<- df[,	{				
				CD4_b200_n<- which(CD4_medbin=="-Inf,200")
				CD4_200_n<- which(CD4_medbin=="200,350")
				CD4_350_n<- which(CD4_medbin=="350,500")
				CD4_500_n<- which(CD4_medbin=="500,Inf")
				CD4_NA_n<- which(CD4_medbin=="unknown")
				list(Acute_n=sum(Acute_p), T_n=sum(AnyT_p), CD4_b200_n=sum(NotYetT_p[CD4_b200_n]), CD4_200_n=sum(NotYetT_p[CD4_200_n]), CD4_350_n=sum(NotYetT_p[CD4_350_n]), CD4_500_n=sum(NotYetT_p[CD4_500_n]), CD4_NA_n=sum(NotYetT_p[CD4_NA_n]))	
			},by=Db_yr]
	setkey(t.living, Db_yr)
	tmp					<- t.living[,Db_yr]
	t.living			<- t( as.matrix(subset(t.living, select=c(Acute_n, CD4_b200_n, CD4_200_n, CD4_350_n, CD4_500_n, CD4_NA_n, T_n))) )
	colnames(t.living)	<- tmp
	rownames(t.living)	<- c("Acute","<200","200-349","350-499",">=500","unknown","Treated")
	
	if(!is.null(plot.file))
	{
		cols				<- brewer.pal(nrow(t.living),"Set1")
		if(verbose)	cat(paste("\nplot file prop to ",plot.file))
		pdf(file=plot.file, width=5, height=4)
		par(mar=c(4,6,0.5,0.5))
		my.barplot.table(t.living, 0.4, "year", plot.ylab, cols, x.baseline=0, xax=as.numeric(colnames(t.living)), legend.loc="topleft")
		dev.off()
	}
	p.living			<- t.living / matrix(rep(apply(t.living,2,sum), each=nrow(t.living)), nrow(t.living), ncol(t.living))
	if(!is.null(plot.file.p))
	{
		if(verbose)	cat(paste("\nplot file prop to ",plot.file.p))
		pdf(file=plot.file.p, width=5, height=4)
		par(mar=c(4,6,0.5,0.5))
		cols				<- brewer.pal(nrow(p.living),"Set1")
		my.barplot.table(p.living, 0.4, "year", paste("%",plot.ylab), cols, x.baseline=0, xax=as.numeric(colnames(p.living)), legend.loc=NULL)
		dev.off()
	}
	list(t.living=t.living,p.living=p.living)
}
######################################################################################
hivc.db.getplot.livingbyexposure<- function(df, plot.file, plot.file.p, plot.ylab, db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)
{
	#simplify risk group 
	set(df,which(df[,Trm=="SXCH"]),"Trm","OTH")
	set(df,which(df[,Trm=="PREG"]),"Trm","OTH")
	set(df,which(df[,Trm=="NEEACC"]),"Trm","OTH")
	set(df,which(df[,Trm=="BLOOD"]),"Trm","OTH")
	set(df,which(df[,Trm=="HETfa"]),"Trm","HET")	
	set(df,which(is.na(df[,Trm])),"Trm","unknown")	
	set(df, NULL, "Trm", factor(df[,Trm]))
	#
	set(df, NULL, "AnyPos_T1", hivc.db.Date2numeric(df[,AnyPos_T1]))
	set(df, NULL, "DateDied", hivc.db.Date2numeric(df[,DateDied]))
	set(df, NULL, "DateLastContact", hivc.db.Date2numeric(df[,DateLastContact]))
	df[,DateEnd:= DateDied]
	#compute DateEnd, either death or lost contact
	tmp				<- which( df[, (DateDied-DateLastContact)>db.diff.lastcontact2died] )
	if(verbose)	cat(paste("\n(DateDied-DateLastContact)>db.diff.lastcontact2died for n=",length(tmp)))
	set(df, tmp, "DateEnd", df[tmp,DateLastContact])
	tmp				<- which(df[, is.na(DateDied) & (db.endtime-DateLastContact)>db.diff.lastcontact2now])
	if(verbose)	cat(paste("\n(db.endtime-DateLastContact)>db.diff.lastcontact2now for n=",length(tmp)))
	set(df, tmp, "DateEnd", df[tmp,DateLastContact])
	tmp				<- which(is.na(df[,DateLastContact]))
	if(verbose)	cat(paste("\nnumber patients with missing DateLastContact n=",length(tmp)))
	tmp				<- which(is.na(df[,DateEnd]))
	if(verbose)	cat(paste("\nnumber patients still alive n=",length(tmp)))
	set(df, tmp,"DateEnd",db.endtime)
	df				<- subset(df,select=c(Patient, AnyPos_T1, Trm, DateEnd))	
	df[, AnyPos_yr:= floor(AnyPos_T1)]
	df[, DateEnd_yr:= floor(DateEnd)]
	#compute patients alive in year Db_yr 
	tmp				<- seq.int(min(df[, AnyPos_yr]),max(df[, DateEnd_yr]))
	tmp				<- lapply(tmp, function(x)
			{				
				tmp	<- subset(df, AnyPos_yr<=x & DateEnd_yr>=x)
				cbind(tmp,data.table(Db_yr=rep(x,nrow(tmp))))
			})
	df				<- rbindlist(tmp)
	#
	t.living			<- table(df[,Db_yr,Trm])
	rownames(t.living)[rownames(t.living)=="BI"]<- "Bi"
	rownames(t.living)[rownames(t.living)=="HET"]<- "Het"
	rownames(t.living)[rownames(t.living)=="IDU"]<- "DU"
	rownames(t.living)[rownames(t.living)=="OTH"]<- "other"	
	t.living	<- t.living[c("MSM","Bi",setdiff(rownames(t.living),c("MSM","Bi","unknown")),"unknown"),]
	#	
	if(!is.null(plot.file))
	{
		cols				<- brewer.pal(nrow(t.living),"Set2")
		if(verbose)	cat(paste("\nplot file prop to ",plot.file))
		pdf(file=plot.file, width=5, height=4)
		par(mar=c(4,6,0.5,0.5))
		my.barplot.table(t.living, 0.4, "year", plot.ylab, cols, x.baseline=0, xax=as.numeric(colnames(t.living)), legend.loc="topleft")
		dev.off()
	}		
	p.living			<- t.living / matrix(rep(apply(t.living,2,sum), each=nrow(t.living)), nrow(t.living), ncol(t.living))
	if(!is.null(plot.file.p))
	{
		if(verbose)	cat(paste("\nplot file prop to ",plot.file.p))
		pdf(file=plot.file.p, width=5, height=4)
		par(mar=c(4,6,0.5,0.5))
		cols				<- brewer.pal(nrow(p.living),"Set2")
		my.barplot.table(p.living, 0.4, "year", paste("%",plot.ylab), cols, x.baseline=0, xax=as.numeric(colnames(p.living)), legend.loc=NULL)
		dev.off()
	}
	list(t.living=t.living,p.living=p.living)
}
######################################################################################
hivc.db.reset.Date<- function(df, col, NA.time, date.format="%Y-%m-%d")
{
	for(x in col)
	{
		cat(paste("\nprocess column", x))
		for(y in NA.time)
		{			
			nok.idx				<- which( df[,x]==y )
			cat(paste("\nentries with format ",y,", n=", length(nok.idx)))						
			if(length(nok.idx))	
				df[nok.idx,x]	<- NA			
		}
		df[,x]	<- as.Date(df[,x], format=date.format)	
	}		
	df
}
######################################################################################
hivc.db.reset.ARTStartFromAccurate<- function(df, col="AnyT_T1")
{
	nacc				<- which(df[, AnyT_T1_Acc=="NAccD"])		
	tmp					<- as.POSIXlt(unclass(df[nacc, col, with=FALSE])[[col]])
	tmp$mon				<- tmp$mon+1
	tmp$mday			<- 1
	set(df, nacc, col, as.Date(tmp))
	nacc				<- which(df[, AnyT_T1_Acc%in%c("NAccMD","NAccYMD")])
	tmp					<- as.POSIXlt(unclass(df[nacc, col, with=FALSE])[[col]])
	tmp$mday			<- 31
	tmp$mon				<- 11
	set(df, nacc, col, as.Date(tmp))
	df
}
######################################################################################
hivc.db.Date2numeric<- function( x )
{
	if(!class(x)%in%c('Date','character'))	return( x )
	x	<- as.POSIXlt(x)
	tmp	<- x$year + 1900
	x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
	x	
}
