######################################################################################
eval.seq.first.sequences.numbers<- function()
{
	infile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	outfile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	load(infile)
	#
	#	sequence numbers per individual	
	#
	dfs		<- subset(df.all, !is.na(FASTASampleCode) & !is.na(PosSeqT))	
	dfs		<- dfs[, list(FASTASampleCode=FASTASampleCode[which.min(PosSeqT)], PosSeqT=min(PosSeqT), AnyPos_T1=AnyPos_T1[1]), by=c('Patient','Trm','Region_first','SUBTYPE_C')]
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unkown')
	set(dfs, dfs[, which(is.na(Region_first))], 'Region_first', 'Unkown')		
	set(dfs, dfs[, which(SUBTYPE_C%in%c('A1','A2'))], 'SUBTYPE_C', 'A')
	set(dfs, dfs[, which(SUBTYPE_C%in%c('F1','F2'))], 'SUBTYPE_C', 'F')
	set(dfs, dfs[, which(SUBTYPE_C%in%c('D-check'))], 'SUBTYPE_C', 'D')
	set(dfs, dfs[, which(SUBTYPE_C%in%c('H','J','K'))], 'SUBTYPE_C', 'Other')
	set(dfs, NULL, 'Time2Seq', dfs[,PosSeqT-AnyPos_T1])
	set(dfs, NULL, 'AnyPos_T1_Y', dfs[,floor(AnyPos_T1)])
	#
	#	plot by risk group
	#
	ggplot( dfs, aes(x= floor(PosSeqT), fill=Trm) ) + geom_bar() +
			labs(x='sampling time of first HIV sequence per patient', y='number / year', fill='risk group') +			
			scale_x_continuous(breaks=seq(1980, 2020, 2)) +
			facet_grid(Region_first~.) +
			scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(legend.position = "bottom") 
	ggsave(file=paste(outfile,'seqno_by_region.pdf'), w=6, h=10)
	#
	#	plot by subtype and risk group
	#
	ggplot( dfs, aes(x= floor(PosSeqT), fill=SUBTYPE_C) ) + geom_bar() +
			labs(x='sampling time of first HIV sequence per patient', y='number / year', fill='risk group') +			
			scale_x_continuous(breaks=seq(1980, 2020, 2)) +
			facet_grid(Trm~., scales='free_y') +
			scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(legend.position = "bottom") 
	ggsave(file=paste(outfile,'seqno_by_riskgroup.pdf'), w=6, h=6)	
}
######################################################################################
eval.seq.first.sequences.timetosequencing<- function()
{
	infile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	outfile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	load(infile)
	#
	#	sequence numbers per individual	
	#
	dfs		<- subset(df.all, !is.na(FASTASampleCode) & !is.na(PosSeqT))	
	dfs		<- dfs[, list(FASTASampleCode=FASTASampleCode[which.min(PosSeqT)], PosSeqT=min(PosSeqT), AnyPos_T1=AnyPos_T1[1]), by=c('Patient','Trm','Region_first','SUBTYPE_C')]
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unkown')
	set(dfs, dfs[, which(is.na(Region_first))], 'Region_first', 'Unkown')		
	set(dfs, dfs[, which(SUBTYPE_C%in%c('A1','A2'))], 'SUBTYPE_C', 'A')
	set(dfs, dfs[, which(SUBTYPE_C%in%c('F1','F2'))], 'SUBTYPE_C', 'F')
	set(dfs, dfs[, which(SUBTYPE_C%in%c('D-check'))], 'SUBTYPE_C', 'D')
	set(dfs, dfs[, which(SUBTYPE_C%in%c('H','J','K'))], 'SUBTYPE_C', 'Other')
	set(dfs, NULL, 'Time2Seq', dfs[,PosSeqT-AnyPos_T1])
	set(dfs, NULL, 'AnyPos_T1_Y', dfs[,floor(AnyPos_T1)])
	#
	#	time to first sequence taken by region
	#
	tmp		<- dfs[, list(mean=mean(Time2Seq), median=median(Time2Seq)), by=c('AnyPos_T1_Y','Region_first')]
	tmp		<- melt(tmp, id.vars=c('AnyPos_T1_Y','Region_first'))
	tmp		<- subset(tmp, AnyPos_T1_Y>1995)
	ggplot(tmp, aes(x=AnyPos_T1_Y, y=value, linetype=variable, pch=variable)) +
			geom_point() + geom_line() +
			scale_x_continuous(breaks=seq(1980, 2020, 2)) +
			scale_y_continuous(breaks=seq(0,20,1)) +
			facet_grid(Region_first~.) +
			labs(x='sampling time of first HIV sequence per patient', y='time to first sequence taken after diagnosis\n(years)', pch='', linetype='') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'_seqtime_from_diagnosis.pdf'), w=6, h=8)		
}
######################################################################################
#
######################################################################################
eval.migrants.WTprop.160618<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	indir		<- "~/Dropbox (Infectious Disease)/2015_ATHENA_May_Update"
	infile		<- file.path(indir,"ATHENA_1502_All_PatientKeyCovariates.R")
	outdir		<- "/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2016/2016_GuidingTransmissionElimination"
	load(infile)	
	df			<- copy(df.all)
	#	convert dates to numeric
	set(df, NULL, "DateBorn", df[,hivc.db.Date2numeric( DateBorn )])
	set(df, NULL, "DateLastContact", df[,hivc.db.Date2numeric( DateLastContact )])	
	set(df, NULL, "DateDied", df[,hivc.db.Date2numeric( DateDied )])
	set(df, NULL, "NegT", df[,hivc.db.Date2numeric( NegT )])
	set(df, NULL, "NegT_Crude", df[,hivc.db.Date2numeric( NegT_Crude )])
	set(df, NULL, "AnyPos_T1", df[,hivc.db.Date2numeric( AnyPos_T1 )])
	set(df, NULL, "GGD_RecordTime", df[,hivc.db.Date2numeric( GGD_RecordTime )])
	set(df, NULL, "DateInCare", df[,hivc.db.Date2numeric( DateInCare )])
	set(df, NULL, "FirstMed", df[,hivc.db.Date2numeric( FirstMed )])
	set(df, NULL, "PoslRNA_T1", df[,hivc.db.Date2numeric( PoslRNA_T1 )])
	set(df, NULL, "PoslRNAg500_T1", df[,hivc.db.Date2numeric( PoslRNAg500_T1 )])	
	set(df, NULL, "PosCD4_T1", df[,hivc.db.Date2numeric( PosCD4_T1 )])
	set(df, NULL, "AnyT_T1", df[,hivc.db.Date2numeric( AnyT_T1 )])
	set(df, NULL, "AnyT_T1_Crude", df[,hivc.db.Date2numeric( AnyT_T1_Crude )])
	set(df, NULL, "DateAIDS", df[,hivc.db.Date2numeric( DateAIDS )])
	set(df, NULL, "PosSeqT", df[,hivc.db.Date2numeric( PosSeqT )])
	#	redefine RegionOrigin
	set(df, df[, which(RegionOrigin%in%c("Central_EU"))], "RegionOrigin", "Eastern_EU_stans")
	set(df, df[, which(RegionOrigin%in%c("Austr_NZ","Sout_SouthEast_Asia","Oceania_Pacific","North_Africa_Middle_East"))], "RegionOrigin", "Other")	
	set(df, NULL, "RegionOrigin", df[, factor(RegionOrigin)])
	#	redefine Tansmission
	set(df, df[, which(Trm%in%c("MSM",'BI'))], "Trm", "MSM")
	set(df, df[, which(Trm%in%c("HET",'HETfa'))], "Trm", "HET")
	set(df, df[, which(Trm%in%c("IDU","BLOOD","NEEACC","PREG","BREAST","SXCH"))], "Trm", "OTH")
	set(df, df[,which(is.na(Trm))], "Trm", "Unknown")
	set(df, NULL, "Trm", df[, factor(Trm)])
	#	age at diagnosis
	set(df, NULL, "Age_AnyPosT1", df[, AnyPos_T1-DateBorn])	
	#	reduce to patients with one sequence
	df[, DUMMY:=seq_len(nrow(df))]
	tmp			<- df[, list(DUMMY=ifelse(all(is.na(PosSeqT)), DUMMY, DUMMY[which.min(PosSeqT)] )), by='Patient']
	df			<- merge(df, tmp, by=c('Patient','DUMMY'))
	df[,DUMMY:=NULL]
	#	define Amsterdam
	df[, AMST:= NA_character_]
	set(df, df[, which( Region_RegT1=='Amst')], 'AMST', 'Y')
	set(df, df[, which( Region_RegT1!='Amst')], 'AMST', 'N')
	set(df, df[, which( is.na(Region_RegT1) & Region_now=='Amst')], 'AMST', 'Y')
	set(df, df[, which( is.na(Region_RegT1) & Region_now!='Amst')], 'AMST', 'N')
	set(df, df[, which( !is.na(CountryInfection) & CountryInfection!='NL')], 'AMST', 'N')	
	stopifnot( !nrow(subset(df, is.na(AMST) & AnyPos_T1>2010)) )	
	#	define time period
	set(df, NULL, 'TP', df[, cut(AnyPos_T1, breaks=c(-Inf, 2010, 2011, 2012, 2013, 2014, 2015, 2016), labels=c('<2010','2010','2011','2012','2013','2014','2015'))])
	#	define young
	set(df, NULL, 'YOUNG', df[, cut(Age_AnyPosT1, breaks=c(0,28,100), labels=c('16-27','28-80'))])
	#	define migrant
	set(df, NULL, 'MIGRANT', df[, RegionOrigin])
	#
	df			<- subset(df, !is.na(AnyPos_T1))	
	
	subset(df, AMST=='Y')[, table(Trm, MIGRANT)]
	
	#
	#	proportion of migrants in Amsterdam diagnosed early
	#
	pea		<- subset(df, !is.na(isAcute) & AMST=='Y' & MIGRANT!='NL' & TP!='<2010' )[, list(PROP_EARLY=mean(isAcute=='Yes')), by="Trm"]
	pead	<- subset(df, !is.na(isAcute) & AMST=='Y' & MIGRANT=='NL' & TP!='<2010' )[, list(PROP_EARLY=mean(isAcute=='Yes')), by="Trm"]
	peo		<- subset(df, !is.na(isAcute) & AMST=='N' & MIGRANT!='NL' & TP!='<2010' )[, list(PROP_EARLY=mean(isAcute=='Yes')), by="Trm"]
	peod	<- subset(df, !is.na(isAcute) & AMST=='N' & MIGRANT=='NL' & TP!='<2010' )[, list(PROP_EARLY=mean(isAcute=='Yes')), by="Trm"]
	#
	#	proportion of migrants in Amsterdam diagnosed late, of those with a first CD4 count in the first year after diagnosis and not on ART before 
	#	with first CD4<350
	#
	pla		<- subset(df, AMST=='Y' & MIGRANT!='NL' & TP!='<2010' & !is.na(PosCD4_T1) & AnyPos_T1+1>PosCD4_T1 & PosCD4_T1<=AnyT_T1)[, list(PROP_LATE=mean(CD4_T1<=350)), by='Trm']
	#
	#	numbers of sequences and isAcute
	tmp			<- subset(df, TP!='<2010' & TP!='2015')
	tmp[, ADJ:=1]
	tmp2		<- subset(df, TP=='2014')
	tmp2[, ADJ:=0.85]
	tmp2[, TP:='2015']
	tmp			<- rbind(tmp, tmp2)
	tmp2[, ADJ:=0.87^2]
	tmp2[, TP:='2016']
	tmp			<- rbind(tmp, tmp2)
	tmp2[, ADJ:=0.87^3]
	tmp2[, TP:='2017']
	tmp			<- rbind(tmp, tmp2)
	tmp2[, ADJ:=0.87^4]
	tmp2[, TP:='2018']
	tmp			<- rbind(tmp, tmp2)
	tmp			<- tmp[, {
				z										<- list()	
				z[['A_Migrant_Any_Any_MSM']]			<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL')])
				z[['O_Migrant_Any_Any_MSM']]			<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL')])
				z[['A_Dutch_Any_Any_MSM']]				<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
				z[['O_Dutch_Any_Any_MSM']]				<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
				
				z[['A_Migrant_Seq_Any_MSM']]			<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & !is.na(PosSeqT))])
				z[['O_Migrant_Seq_Any_MSM']]			<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & !is.na(PosSeqT))])
				z[['A_Dutch_Seq_Any_MSM']]				<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL' & !is.na(PosSeqT))])
				z[['O_Dutch_Seq_Any_MSM']]				<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL' & !is.na(PosSeqT))])
				
				z[['A_Migrant_Any_Rec_MSM']]			<- z[['A_Migrant_Any_Any_MSM']] * subset(pea, Trm=='MSM')[, PROP_EARLY]	#sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & isAcute=='Yes')])
				z[['O_Migrant_Any_Rec_MSM']]			<- z[['O_Migrant_Any_Any_MSM']]	* subset(peo, Trm=='MSM')[, PROP_EARLY]#sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & isAcute=='Yes')])
				z[['A_Dutch_Any_Rec_MSM']]				<- z[['A_Dutch_Any_Any_MSM']] * subset(pead, Trm=='MSM')[, PROP_EARLY]	#sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL' & isAcute=='Yes')])
				z[['O_Dutch_Any_Rec_MSM']]				<- z[['O_Dutch_Any_Any_MSM']] * subset(peod, Trm=='MSM')[, PROP_EARLY]	#sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL' & isAcute=='Yes')])
				
				
				z[['A_Migrant_Seq_Rec_MSM']]			<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				z[['O_Migrant_Seq_Rec_MSM']]			<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				z[['A_Dutch_Seq_Rec_MSM']]				<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				z[['O_Dutch_Seq_Rec_MSM']]				<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				
				z[['A_Migrant_Any_Any_HET']]			<- sum(ADJ[which(AMST=='Y' & Trm=='HET' & MIGRANT!='NL')])
				z[['O_Migrant_Any_Any_HET']]			<- sum(ADJ[which(AMST=='N' & Trm=='HET'& MIGRANT!='NL')])
				z[['A_Dutch_Any_Any_HET']]				<- sum(ADJ[which(AMST=='Y' & Trm=='HET'& MIGRANT=='NL')])
				z[['O_Dutch_Any_Any_HET']]				<- sum(ADJ[which(AMST=='N' & Trm=='HET'& MIGRANT=='NL')])
				
				z[['A_Migrant_Seq_Any_HET']]			<- sum(ADJ[which(AMST=='Y' & Trm=='HET'& MIGRANT!='NL' & !is.na(PosSeqT))])
				z[['O_Migrant_Seq_Any_HET']]			<- sum(ADJ[which(AMST=='N' & Trm=='HET'& MIGRANT!='NL' & !is.na(PosSeqT))])
				z[['A_Dutch_Seq_Any_HET']]				<- sum(ADJ[which(AMST=='Y' & Trm=='HET'& MIGRANT=='NL' & !is.na(PosSeqT))])
				z[['O_Dutch_Seq_Any_HET']]				<- sum(ADJ[which(AMST=='N' & Trm=='HET'& MIGRANT=='NL' & !is.na(PosSeqT))])
				
				z[['A_Migrant_Any_Rec_HET']]			<- z[['A_Migrant_Any_Any_HET']]	* subset(pea, Trm=='HET')[, PROP_EARLY]#sum(ADJ[which(AMST=='Y' & MIGRANT!='NL' & isAcute=='Yes')])
				z[['O_Migrant_Any_Rec_HET']]			<- z[['O_Migrant_Any_Any_HET']]	* subset(peo, Trm=='HET')[, PROP_EARLY]#sum(ADJ[which(AMST=='N' & MIGRANT!='NL' & isAcute=='Yes')])
				z[['A_Dutch_Any_Rec_HET']]				<- z[['A_Dutch_Any_Any_HET']]	* subset(pead, Trm=='HET')[, PROP_EARLY]#sum(ADJ[which(AMST=='Y' & MIGRANT=='NL' & isAcute=='Yes')])
				z[['O_Dutch_Any_Rec_HET']]				<- z[['O_Dutch_Any_Any_HET']]	* subset(peod, Trm=='HET')[, PROP_EARLY]#sum(ADJ[which(AMST=='N' & MIGRANT=='NL' & isAcute=='Yes')])
				
				
				z[['A_Migrant_Seq_Rec_HET']]			<- sum(ADJ[which(AMST=='Y' & Trm=='HET'& MIGRANT!='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				z[['O_Migrant_Seq_Rec_HET']]			<- sum(ADJ[which(AMST=='N' & Trm=='HET'& MIGRANT!='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				z[['A_Dutch_Seq_Rec_HET']]				<- sum(ADJ[which(AMST=='Y' & Trm=='HET'& MIGRANT=='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				z[['O_Dutch_Seq_Rec_HET']]				<- sum(ADJ[which(AMST=='N' & Trm=='HET'& MIGRANT=='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				
				z
			}, by='TP']
	cnt			<- melt(tmp, id.vars='TP')
	cnt[, LOC:= factor(substr(variable, 1, 1), levels=c('A','O'), labels=c('Amst','Other'))]
	cnt[, ORIGIN:= cnt[, sapply(strsplit(as.character(variable), '_'), '[[', 2)]]
	cnt[, SEQ:= cnt[, sapply(strsplit(as.character(variable), '_'), '[[', 3)]]	
	cnt[, REC:= cnt[, sapply(strsplit(as.character(variable), '_'), '[[', 4)]]
	cnt[, TRM:= cnt[, sapply(strsplit(as.character(variable), '_'), '[[', 5)]]
	#set(cnt, NULL, 'value', cnt[, round(value)])
	set(cnt, NULL, 'variable', NULL)
	#	aggregate
	cnt			<- rbind(	subset(cnt, TP%in%c('2010','2011','2012','2013','2014'))[, list(TP='2010-2014', value=round(sum(value))), by=c('LOC','ORIGIN','SEQ','REC','TRM')],
			subset(cnt, TP%in%c('2015','2016','2017','2018'))[, list(TP='2015-2018', value=round(sum(value))), by=c('LOC','ORIGIN','SEQ','REC','TRM')],
			subset(cnt, !TP%in%c('<2010'))[, list(TP='2010-2018', value=round(sum(value))), by=c('LOC','ORIGIN','SEQ','REC','TRM')]	
	)
	set(cnt, NULL, 'TP', cnt[, factor(TP, levels=c('2010-2014','2015-2018','2010-2018'))])	
	
	dcast.data.table(subset(cnt, REC=='Any' & SEQ=='Any' & LOC=='Amst'), TP+ORIGIN~TRM, value.var='value')
	dcast.data.table(subset(cnt, REC=='Rec' & SEQ=='Any' & LOC=='Amst'), TP+ORIGIN~TRM, value.var='value')
	
	#
	# 	expected number of recipients and actual transmission pairs captured, 2010-2018:
	#	baseline:
	if(1)
	{
		m.A			<- 0.99 * 0.9*(0.7*.99+.3*.45)
		m.O			<- 0.45 * 0.9*(0.3*.99+.7*.45)		
	}
	#	rather than RITA, let s get 50% coverage outside Amsterdam
	if(0)
	{
		m.A			<- 0.95 * 0.9*(0.7*.95+.3*.5)
		m.O			<- 0.5 * 0.9*(0.3*.95+.7*.5)		
	}
	#	for young, let s get 75% coverage outside Amsterdam
	if(0)
	{
		m.A			<- 0.95 * 0.9*(0.7*.95+.3*.75)
		m.O			<- 0.75 * 0.9*(0.3*.95+.7*.75)		
	}
	tmp			<- subset(df, TP!='<2010' & TP!='2015')
	tmp[, ADJ:=1]
	tmp2		<- subset(df, TP=='2014')
	tmp2[, ADJ:=0.85]
	tmp2[, TP:='2015']
	tmp			<- rbind(tmp, tmp2)
	tmp2[, ADJ:=0.87^2]
	tmp2[, TP:='2016']
	tmp			<- rbind(tmp, tmp2)
	tmp2[, ADJ:=0.87^3]
	tmp2[, TP:='2017']
	tmp			<- rbind(tmp, tmp2)
	tmp2[, ADJ:=0.87^4]
	tmp2[, TP:='2018']
	tmp			<- rbind(tmp, tmp2)
	tmp			<- tmp[, {
				z			<- list()												
				z[['A_NL_MSM_ALL_1YE']]				<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(pead, Trm=='MSM')[, PROP_EARLY]
				z[['A_NL_HETM_ALL_1YE']]			<- sum(ADJ[which(AMST=='Y' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]
				z[['A_NL_HETF_ALL_1YE']]			<- sum(ADJ[which(AMST=='Y' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]				
				z[['A_Migrant_ALL_MSM_1YE']]		<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(pea, Trm=='MSM')[, PROP_EARLY]
				z[['A_Migrant_HETM_ALL_1YE']]		<- sum(ADJ[which(AMST=='Y' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
				z[['A_Migrant_HETF_ALL_1YE']]		<- sum(ADJ[which(AMST=='Y' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
				
				z[['A_NL_MSM_ALL_1Y']]				<- sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
				z[['A_NL_HETM_ALL_1Y']]				<- sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
				z[['A_NL_HETF_ALL_1Y']]				<- sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
				z[['A_Migrant_MSM_ALL_1Y']]			<- sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
				z[['A_Migrant_HETM_ALL_1Y']]		<- sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				z[['A_Migrant_HETF_ALL_1Y']]		<- sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				
				z[['A_NL_MSM_28_1YE']]				<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(pead, Trm=='MSM')[, PROP_EARLY]
				z[['A_NL_HETM_28_1YE']]				<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]
				z[['A_NL_HETF_28_1YE']]				<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]				
				z[['A_Migrant_28_MSM_1YE']]			<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(pea, Trm=='MSM')[, PROP_EARLY]
				z[['A_Migrant_HETM_28_1YE']]		<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
				z[['A_Migrant_HETF_28_1YE']]		<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
				
				z[['A_NL_MSM_28_1Y']]				<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
				z[['A_NL_HETM_28_1Y']]				<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
				z[['A_NL_HETF_28_1Y']]				<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
				z[['A_Migrant_MSM_28_1Y']]			<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
				z[['A_Migrant_HETM_28_1Y']]			<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				z[['A_Migrant_HETF_28_1Y']]			<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				
				if(1)	#since recently infected, in-country infection more likely
				{
					z[['A_NL_MSM_ALL_PE']]			<- m.A*sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(pead, Trm=='MSM')[, PROP_EARLY]
					z[['A_NL_HETM_ALL_PE']]			<- m.A*sum(ADJ[which(AMST=='Y' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]
					z[['A_NL_HETF_ALL_PE']]			<- m.A*sum(ADJ[which(AMST=='Y' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]				
					z[['A_Migrant_MSM_ALL_PE']]		<- m.A*sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(pea, Trm=='MSM')[, PROP_EARLY]
					z[['A_Migrant_HETM_ALL_PE']]	<- m.A*sum(ADJ[which(AMST=='Y' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					z[['A_Migrant_HETF_ALL_PE']]	<- m.A*sum(ADJ[which(AMST=='Y' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					
					z[['A_NL_MSM_ALL_P']]			<- m.A*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) 
					z[['A_NL_HETM_ALL_P']]			<- m.A*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
					z[['A_NL_HETF_ALL_P']]			<- m.A*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
					z[['A_Migrant_MSM_ALL_P']]		<- m.A*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETM_ALL_P']]		<- m.A*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETF_ALL_P']]		<- m.A*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					
					z[['A_NL_MSM_28_PE']]			<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(pead, Trm=='MSM')[, PROP_EARLY]
					z[['A_NL_HETM_28_PE']]			<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]
					z[['A_NL_HETF_28_PE']]			<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]				
					z[['A_Migrant_MSM_28_PE']]		<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(pea, Trm=='MSM')[, PROP_EARLY]
					z[['A_Migrant_HETM_28_PE']]		<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					z[['A_Migrant_HETF_28_PE']]		<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					
					z[['A_NL_MSM_28_P']]			<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) 
					z[['A_NL_HETM_28_P']]			<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
					z[['A_NL_HETF_28_P']]			<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
					z[['A_Migrant_MSM_28_P']]		<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETM_28_P']]		<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETF_28_P']]		<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
				}
				if(0)	#prob in-country infection from aMASE
				{
					z[['A_NL_MSM_ALL_PE']]			<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(pead, Trm=='MSM')[, PROP_EARLY]
					z[['A_NL_HETM_ALL_PE']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]
					z[['A_NL_HETF_ALL_PE']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]				
					z[['A_Migrant_MSM_ALL_PE']]		<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(pea, Trm=='MSM')[, PROP_EARLY]
					z[['A_Migrant_HETM_ALL_PE']]	<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					z[['A_Migrant_HETF_ALL_PE']]	<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					
					z[['A_NL_MSM_ALL_P']]			<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) 
					z[['A_NL_HETM_ALL_P']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
					z[['A_NL_HETF_ALL_P']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
					z[['A_Migrant_MSM_ALL_P']]		<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETM_ALL_P']]		<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETF_ALL_P']]		<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					
					z[['A_NL_MSM_28_PE']]			<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(pead, Trm=='MSM')[, PROP_EARLY]
					z[['A_NL_HETM_28_PE']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]
					z[['A_NL_HETF_28_PE']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]				
					z[['A_Migrant_MSM_28_PE']]		<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(pea, Trm=='MSM')[, PROP_EARLY]
					z[['A_Migrant_HETM_28_PE']]		<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					z[['A_Migrant_HETF_28_PE']]		<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					
					z[['A_NL_MSM_28_P']]			<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) 
					z[['A_NL_HETM_28_P']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
					z[['A_NL_HETF_28_P']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
					z[['A_Migrant_MSM_28_P']]		<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETM_28_P']]		<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETF_28_P']]		<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
				}				
				#z[['O']]						<- sum(ADJ[which(AMST=='N')])				
				z[['O_NL_MSM_ALL_1YE']]				<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(peod, Trm=='MSM')[, PROP_EARLY]
				z[['O_NL_HETM_ALL_1YE']]			<- sum(ADJ[which(AMST=='N' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]
				z[['O_NL_HETF_ALL_1YE']]			<- sum(ADJ[which(AMST=='N' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]								
				z[['O_Migrant_MSM_ALL_1YE']]		<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(peo, Trm=='MSM')[, PROP_EARLY]
				z[['O_Migrant_HETM_ALL_1YE']]		<- sum(ADJ[which(AMST=='N' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
				z[['O_Migrant_HETF_ALL_1YE']]		<- sum(ADJ[which(AMST=='N' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
				
				z[['O_NL_MSM_ALL_1Y']]				<- sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
				z[['O_NL_HETM_ALL_1Y']]				<- sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
				z[['O_NL_HETF_ALL_1Y']]				<- sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
				z[['O_Migrant_MSM_ALL_1Y']]			<- sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
				z[['O_Migrant_HETM_ALL_1Y']]		<- sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				z[['O_Migrant_HETF_ALL_1Y']]		<- sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				
				z[['O_NL_MSM_28_1YE']]				<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(peod, Trm=='MSM')[, PROP_EARLY]
				z[['O_NL_HETM_28_1YE']]				<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]
				z[['O_NL_HETF_28_1YE']]				<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]								
				z[['O_Migrant_MSM_28_1YE']]			<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(peo, Trm=='MSM')[, PROP_EARLY]
				z[['O_Migrant_HETM_28_1YE']]		<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
				z[['O_Migrant_HETF_28_1YE']]		<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
				
				z[['O_NL_MSM_28_1Y']]				<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
				z[['O_NL_HETM_28_1Y']]				<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
				z[['O_NL_HETF_28_1Y']]				<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
				z[['O_Migrant_MSM_28_1Y']]			<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
				z[['O_Migrant_HETM_28_1Y']]			<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				z[['O_Migrant_HETF_28_1Y']]			<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				
				if(1)	#since recently infected, in-country infection more likely
				{
					z[['O_NL_MSM_ALL_PE']]			<- m.O*sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(peod, Trm=='MSM')[, PROP_EARLY]
					z[['O_NL_HETM_ALL_PE']]			<- m.O*sum(ADJ[which(AMST=='N' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]
					z[['O_NL_HETF_ALL_PE']]			<- m.O*sum(ADJ[which(AMST=='N' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]								
					z[['O_Migrant_MSM_ALL_PE']]		<- m.O*sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(peo, Trm=='MSM')[, PROP_EARLY]
					z[['O_Migrant_HETM_ALL_PE']]	<- m.O*sum(ADJ[which(AMST=='N' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					z[['O_Migrant_HETF_ALL_PE']]	<- m.O*sum(ADJ[which(AMST=='N' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					
					z[['O_NL_MSM_ALL_P']]			<- m.O*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
					z[['O_NL_HETM_ALL_P']]			<- m.O*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])
					z[['O_NL_HETF_ALL_P']]			<- m.O*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])								
					z[['O_Migrant_MSM_ALL_P']]		<- m.O*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')])
					z[['O_Migrant_HETM_ALL_P']]		<- m.O*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					z[['O_Migrant_HETF_ALL_P']]		<- m.O*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					
					z[['O_NL_MSM_28_PE']]			<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(peod, Trm=='MSM')[, PROP_EARLY]
					z[['O_NL_HETM_28_PE']]			<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]
					z[['O_NL_HETF_28_PE']]			<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]								
					z[['O_Migrant_MSM_28_PE']]		<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(peo, Trm=='MSM')[, PROP_EARLY]
					z[['O_Migrant_HETM_28_PE']]		<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					z[['O_Migrant_HETF_28_PE']]		<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					
					z[['O_NL_MSM_28_P']]			<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
					z[['O_NL_HETM_28_P']]			<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])
					z[['O_NL_HETF_28_P']]			<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])								
					z[['O_Migrant_MSM_28_P']]		<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')])
					z[['O_Migrant_HETM_28_P']]		<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					z[['O_Migrant_HETF_28_P']]		<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
				}				
				if(0)	#prob in-country infection from aMASE
				{
					z[['O_NL_MSM_ALL_PE']]			<- 0.5*0.54*sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(peod, Trm=='MSM')[, PROP_EARLY]
					z[['O_NL_HETM_ALL_PE']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]
					z[['O_NL_HETF_ALL_PE']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]								
					z[['O_Migrant_MSM_ALL_PE']]		<- 0.5*0.54*sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(peo, Trm=='MSM')[, PROP_EARLY]
					z[['O_Migrant_HETM_ALL_PE']]	<- 0.5*0.36*sum(ADJ[which(AMST=='N' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					z[['O_Migrant_HETF_ALL_PE']]	<- 0.5*0.36*sum(ADJ[which(AMST=='N' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					
					z[['O_NL_MSM_ALL_P']]			<- 0.5*0.54*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
					z[['O_NL_HETM_ALL_P']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])
					z[['O_NL_HETF_ALL_P']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])								
					z[['O_Migrant_MSM_ALL_P']]		<- 0.5*0.54*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')])
					z[['O_Migrant_HETM_ALL_P']]		<- 0.5*0.36*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					z[['O_Migrant_HETF_ALL_P']]		<- 0.5*0.36*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					
					z[['O_NL_MSM_28_PE']]			<- 0.5*0.54*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(peod, Trm=='MSM')[, PROP_EARLY]
					z[['O_NL_HETM_28_PE']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]
					z[['O_NL_HETF_28_PE']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]								
					z[['O_Migrant_MSM_28_PE']]		<- 0.5*0.54*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(peo, Trm=='MSM')[, PROP_EARLY]
					z[['O_Migrant_HETM_28_PE']]		<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					z[['O_Migrant_HETF_28_PE']]		<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					
					z[['O_NL_MSM_28_P']]			<- 0.5*0.54*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
					z[['O_NL_HETM_28_P']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])
					z[['O_NL_HETF_28_P']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])								
					z[['O_Migrant_MSM_28_P']]		<- 0.5*0.54*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')])
					z[['O_Migrant_HETM_28_P']]		<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					z[['O_Migrant_HETF_28_P']]		<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])					
				}				
				z
			}, by='TP']
	cnt			<- melt(tmp, id.vars='TP')
	cnt[, LOC:= factor(substr(variable, 1, 1), levels=c('A','O'), labels=c('Amst','Other'))]
	cnt[, WHAT:= regmatches(variable,regexpr('[^_]*$', variable))]
	cnt[, add_RITA:= factor(grepl('1YE|PE', WHAT), levels=c(TRUE,FALSE), labels=c('Y','N'))]
	set(cnt, NULL, 'WHAT', cnt[, factor(grepl('1Y',WHAT), levels=c(TRUE,FALSE),labels=c('recent', 'pair'))])	
	set(cnt, NULL, 'YOUNG', cnt[, factor(grepl('28',variable), levels=c(TRUE,FALSE),labels=c('<28', 'all'))])	
	cnt[, TRM:= cnt[, sapply(strsplit(as.character(variable), '_'), '[[', 3)]]
	cnt[, MIGRANT:= cnt[, factor(sapply(strsplit(as.character(variable), '_'), '[[', 2), levels=c('Migrant','NL'), labels=c('Y','N'))]]	
	
	subset(cnt, add_RITA=='Y' & YOUNG=='all' & LOC=='Amst' & WHAT=='pair')[, list(SC='RITA_1018', N=round(sum(value))), by=c('WHAT','YOUNG','MIGRANT','TRM','LOC')]
	
	
	cnt			<- rbind( 	subset(cnt,SC_RITA_1016=='Y')[, list(SC='RITA_1016', N=round(sum(value))), by=c('WHAT','YOUNG','MIGRANT','TRM','LOC')],
			subset(cnt,SC_RITA_1315=='Y')[, list(SC='RITA_1315', N=round(sum(value))), by=c('WHAT','YOUNG','MIGRANT','TRM','LOC')]	)
	setkey(cnt, SC, WHAT, MIGRANT, TRM, LOC)				
	dcast.data.table( subset(cnt, WHAT=='pair' & YOUNG=='all' & SC=='RITA_1016'), MIGRANT+LOC~TRM, value.var='N')
	dcast.data.table( subset(cnt, WHAT=='pair' & YOUNG=='all' & SC=='RITA_1315'), MIGRANT+LOC~TRM, value.var='N')
	dcast.data.table( subset(cnt, WHAT=='pair' & YOUNG=='<28' & SC=='RITA_1016'), MIGRANT+LOC~TRM, value.var='N')
	dcast.data.table( subset(cnt, WHAT=='pair' & YOUNG=='<28' & SC=='RITA_1315'), MIGRANT+LOC~TRM, value.var='N')
}
######################################################################################
#
######################################################################################
eval.migrants.WTprop.160228<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	indir		<- "~/Dropbox (Infectious Disease)/2015_ATHENA_May_Update"
	infile		<- file.path(indir,"ATHENA_1502_All_PatientKeyCovariates.R")
	outdir		<- "/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2016/2016_GuidingTransmissionElimination"
	load(infile)	
	df			<- copy(df.all)
	#	convert dates to numeric
	set(df, NULL, "DateBorn", df[,hivc.db.Date2numeric( DateBorn )])
	set(df, NULL, "DateLastContact", df[,hivc.db.Date2numeric( DateLastContact )])	
	set(df, NULL, "DateDied", df[,hivc.db.Date2numeric( DateDied )])
	set(df, NULL, "NegT", df[,hivc.db.Date2numeric( NegT )])
	set(df, NULL, "NegT_Crude", df[,hivc.db.Date2numeric( NegT_Crude )])
	set(df, NULL, "AnyPos_T1", df[,hivc.db.Date2numeric( AnyPos_T1 )])
	set(df, NULL, "GGD_RecordTime", df[,hivc.db.Date2numeric( GGD_RecordTime )])
	set(df, NULL, "DateInCare", df[,hivc.db.Date2numeric( DateInCare )])
	set(df, NULL, "FirstMed", df[,hivc.db.Date2numeric( FirstMed )])
	set(df, NULL, "PoslRNA_T1", df[,hivc.db.Date2numeric( PoslRNA_T1 )])
	set(df, NULL, "PoslRNAg500_T1", df[,hivc.db.Date2numeric( PoslRNAg500_T1 )])	
	set(df, NULL, "PosCD4_T1", df[,hivc.db.Date2numeric( PosCD4_T1 )])
	set(df, NULL, "AnyT_T1", df[,hivc.db.Date2numeric( AnyT_T1 )])
	set(df, NULL, "AnyT_T1_Crude", df[,hivc.db.Date2numeric( AnyT_T1_Crude )])
	set(df, NULL, "DateAIDS", df[,hivc.db.Date2numeric( DateAIDS )])
	set(df, NULL, "PosSeqT", df[,hivc.db.Date2numeric( PosSeqT )])
	#	redefine RegionOrigin
	set(df, df[, which(RegionOrigin%in%c("Central_EU"))], "RegionOrigin", "Eastern_EU_stans")
	set(df, df[, which(RegionOrigin%in%c("Austr_NZ","Sout_SouthEast_Asia","Oceania_Pacific","North_Africa_Middle_East"))], "RegionOrigin", "Other")	
	set(df, NULL, "RegionOrigin", df[, factor(RegionOrigin)])
	#	redefine Tansmission
	set(df, df[, which(Trm%in%c("MSM",'BI'))], "Trm", "MSM")
	set(df, df[, which(Trm%in%c("HET",'HETfa'))], "Trm", "HET")
	set(df, df[, which(Trm%in%c("IDU","BLOOD","NEEACC","PREG","BREAST","SXCH"))], "Trm", "OTH")
	set(df, df[,which(is.na(Trm))], "Trm", "Unknown")
	set(df, NULL, "Trm", df[, factor(Trm)])
	#	age at diagnosis
	set(df, NULL, "Age_AnyPosT1", df[, AnyPos_T1-DateBorn])	
	#	reduce to patients with one sequence
	df[, DUMMY:=seq_len(nrow(df))]
	tmp			<- df[, list(DUMMY=ifelse(all(is.na(PosSeqT)), DUMMY, DUMMY[which.min(PosSeqT)] )), by='Patient']
	df			<- merge(df, tmp, by=c('Patient','DUMMY'))
	df[,DUMMY:=NULL]
	#	define Amsterdam
	df[, AMST:= NA_character_]
	set(df, df[, which( Region_RegT1=='Amst')], 'AMST', 'Y')
	set(df, df[, which( Region_RegT1!='Amst')], 'AMST', 'N')
	set(df, df[, which( is.na(Region_RegT1) & Region_now=='Amst')], 'AMST', 'Y')
	set(df, df[, which( is.na(Region_RegT1) & Region_now!='Amst')], 'AMST', 'N')
	set(df, df[, which( !is.na(CountryInfection) & CountryInfection!='NL')], 'AMST', 'N')	
	stopifnot( !nrow(subset(df, is.na(AMST) & AnyPos_T1>2010)) )	
	#	define time period
	set(df, NULL, 'TP', df[, cut(AnyPos_T1, breaks=c(-Inf, 2010, 2011, 2012, 2013, 2014, 2015, 2016), labels=c('<2010','2010','2011','2012','2013','2014','2015'))])
	#	define young
	set(df, NULL, 'YOUNG', df[, cut(Age_AnyPosT1, breaks=c(0,28,100), labels=c('16-27','28-80'))])
	#	define migrant
	set(df, NULL, 'MIGRANT', df[, RegionOrigin])
	#
	df			<- subset(df, !is.na(AnyPos_T1))	
	#
	#	proportion of early/late/unknown early by migrants and Dutch origin in the NL
	#
	tmp		<- subset(df, !is.na(MIGRANT) & Trm%in%c('HET','MSM') & TP!='<2010' )
	set(tmp, NULL, 'MIGRANT', tmp[,factor(as.character(MIGRANT)=='NL', levels=c(TRUE,FALSE), labels=c('Dutch origin','Migrant'))])
	set(tmp, tmp[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')
	set(tmp, NULL, 'Trm', tmp[, factor(Trm, levels=c('HET','HETM','MSM'), labels=c('Heterosexual women','Heterosexual men', 'MSM'))])
	tmp[, STAGE:=NA_character_]
	set(tmp, tmp[, which(isAcute=='Yes')], 'STAGE', 'Acute')
	set(tmp, tmp[, which(!is.na(PosCD4_T1) & AnyPos_T1+1>PosCD4_T1 & PosCD4_T1-1/12<=AnyT_T1 & CD4_T1<350)], 'STAGE', 'Late')
	set(tmp, tmp[, which(is.na(STAGE) & is.na(isAcute))], 'STAGE', 'Incomplete')
	set(tmp, tmp[, which(is.na(STAGE) & is.na(PosCD4_T1))], 'STAGE', 'Incomplete')
	set(tmp, tmp[, which(is.na(STAGE) & !is.na(PosCD4_T1) & (AnyPos_T1+1<=PosCD4_T1 | PosCD4_T1-1/12>AnyT_T1))], 'STAGE', 'Incomplete')
	set(tmp, tmp[, which(is.na(STAGE))], 'STAGE', 'NoInd')
	set(tmp, NULL, 'STAGE', tmp[, factor(STAGE, levels=c('Acute','NoInd','Late','Incomplete'), labels=c('Confirmed recent infection','No indication for recent infection\nor late presentation','Late presenter','Incomplete data'))])	
	col		<- c("#66C2A5", "#FC8D62", "#8DA0CB", "grey80")  
	ggplot(tmp, aes(x=MIGRANT, fill=STAGE)) + geom_bar(position='fill') +
			facet_grid(~Trm) + theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=col) +
			scale_y_continuous(labels=percent) +
			labs(x='', y='New HIV diagnoses in the Netherlands\n2010-2014', fill='Infection stage\nat diagnosis') +
			guides(fill=guide_legend(nrow=2,byrow=TRUE))
	ggsave(file=file.path(outdir,'ATHENA1502_NewDiag1014NL_PropsEarlyLate1.pdf'), w=8, h=5)
	#
	#	proportion of early/late/unknown early by migrants and Dutch origin in Amsterdam
	#
	tmp		<- subset(df, !is.na(MIGRANT) & AMST=='Y' & Trm%in%c('HET','MSM') & TP!='<2010' )
	set(tmp, NULL, 'MIGRANT', tmp[,factor(as.character(MIGRANT)=='NL', levels=c(TRUE,FALSE), labels=c('Dutch origin','Migrant'))])
	set(tmp, tmp[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')
	set(tmp, NULL, 'Trm', tmp[, factor(Trm, levels=c('HET','HETM','MSM'), labels=c('Heterosexual women','Heterosexual men', 'MSM'))])
	tmp[, STAGE:=NA_character_]
	set(tmp, tmp[, which(isAcute=='Yes')], 'STAGE', 'Acute')
	set(tmp, tmp[, which(!is.na(PosCD4_T1) & AnyPos_T1+1>PosCD4_T1 & PosCD4_T1-1/12<=AnyT_T1 & CD4_T1<350)], 'STAGE', 'Late')
	set(tmp, tmp[, which(is.na(STAGE) & is.na(isAcute))], 'STAGE', 'Incomplete')
	set(tmp, tmp[, which(is.na(STAGE) & is.na(PosCD4_T1))], 'STAGE', 'Incomplete')
	set(tmp, tmp[, which(is.na(STAGE) & !is.na(PosCD4_T1) & (AnyPos_T1+1<=PosCD4_T1 | PosCD4_T1-1/12>AnyT_T1))], 'STAGE', 'Incomplete')
	set(tmp, tmp[, which(is.na(STAGE))], 'STAGE', 'NoInd')
	set(tmp, NULL, 'STAGE', tmp[, factor(STAGE, levels=c('Acute','NoInd','Late','Incomplete'), labels=c('Confirmed recent infection','No indication for recent infection\nor late presentation','Late presenter','Incomplete data'))])	
	col		<- c("#66C2A5", "#FC8D62", "#8DA0CB", "grey80")  
	ggplot(tmp, aes(x=MIGRANT, fill=STAGE)) + geom_bar(position='fill') +
			facet_grid(~Trm) + theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=col) +
			scale_y_continuous(labels=percent) +
			labs(x='', y='New HIV diagnoses in Amsterdam\n2010-2014', fill='Infection stage\nat diagnosis') +
			guides(fill=guide_legend(nrow=2,byrow=TRUE))
	ggsave(file=file.path(outdir,'ATHENA1502_NewDiag1014AMST_PropsEarlyLate1.pdf'), w=8, h=5)
	#
	#	proportion of early/late/unknown early by migrant origin 2010-2014 in the Netherlands
	#
	tmp			<- subset(df, !is.na(MIGRANT) & Trm%in%c('MSM','HET') & AnyPos_T1>2010 & AnyPos_T1<2015)	
	set(tmp, tmp[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')
	set(tmp, NULL, 'Trm', tmp[, factor(Trm, levels=c('HET','HETM','MSM'), labels=c('Heterosexual women','Heterosexual men', 'MSM'))])
	set(tmp, NULL, 'MIGRANT', tmp[, factor(as.character(MIGRANT), 	levels=rev(c("NL","Western_EU", "Latin_South_America", "Eastern_EU_stans", "Caribbean", "Sub_Saharan_Africa", "Other")),
							labels=rev(c("Dutch origin","Western Europe", "Latin & South America", "Eastern Europe", "Caribbean", "Sub-Saharan Africa", "Other")))])
	tmp[, STAGE:=NA_character_]
	set(tmp, tmp[, which(isAcute=='Yes')], 'STAGE', 'Acute')
	set(tmp, tmp[, which(!is.na(PosCD4_T1) & AnyPos_T1+1>PosCD4_T1 & PosCD4_T1-1/12<=AnyT_T1 & CD4_T1<350)], 'STAGE', 'Late')
	set(tmp, tmp[, which(is.na(STAGE) & is.na(isAcute))], 'STAGE', 'Incomplete')
	set(tmp, tmp[, which(is.na(STAGE) & is.na(PosCD4_T1))], 'STAGE', 'Incomplete')
	set(tmp, tmp[, which(is.na(STAGE) & !is.na(PosCD4_T1) & (AnyPos_T1+1<=PosCD4_T1 | PosCD4_T1-1/12>AnyT_T1))], 'STAGE', 'Incomplete')
	set(tmp, tmp[, which(is.na(STAGE))], 'STAGE', 'NoInd')
	set(tmp, NULL, 'STAGE', tmp[, factor(STAGE, levels=c('Acute','NoInd','Late','Incomplete'), labels=c('Confirmed recent infection','No indication for recent infection\nor late presentation','Late presenter','Incomplete data'))])	
	col		<- c("#66C2A5", "#FC8D62", "#8DA0CB", "grey80")  
	ggplot(tmp, aes(x=MIGRANT, fill=STAGE)) + geom_bar(position='fill') +
			theme_bw() + theme(legend.position='bottom') +
			coord_flip() + facet_grid(~Trm) + 
			scale_fill_manual(values=col) +
			scale_y_continuous(labels=percent) +
			labs(x='', y='\nNew HIV diagnoses in the Netherlands\n2010-2014', fill='Infection stage\nat diagnosis') +
			guides(fill=guide_legend(nrow=2,byrow=TRUE))
	ggsave(file=file.path(outdir,'ATHENA1502_NewDiag1014NL_PropEarlyLateByOrigin1.pdf'), w=9.5, h=4)
	#
	#	proportion of early/late/unknown early by migrant origin 2010-2014 in Amsterdam
	#
	tmp			<- subset(df, !is.na(MIGRANT) & AMST=='Y' & Trm%in%c('MSM','HET') & AnyPos_T1>2010 & AnyPos_T1<2015)	
	set(tmp, tmp[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')
	set(tmp, NULL, 'Trm', tmp[, factor(Trm, levels=c('HET','HETM','MSM'), labels=c('Heterosexual women','Heterosexual men', 'MSM'))])
	set(tmp, NULL, 'MIGRANT', tmp[, factor(as.character(MIGRANT), 	levels=rev(c("NL","Western_EU", "Latin_South_America", "Eastern_EU_stans", "Caribbean", "Sub_Saharan_Africa", "Other")),
							labels=rev(c("Dutch origin","Western Europe", "Latin & South America", "Eastern Europe", "Caribbean", "Sub-Saharan Africa", "Other")))])
	tmp[, STAGE:=NA_character_]
	set(tmp, tmp[, which(isAcute=='Yes')], 'STAGE', 'Acute')
	set(tmp, tmp[, which(!is.na(PosCD4_T1) & AnyPos_T1+1>PosCD4_T1 & PosCD4_T1-1/12<=AnyT_T1 & CD4_T1<350)], 'STAGE', 'Late')
	set(tmp, tmp[, which(is.na(STAGE) & is.na(isAcute))], 'STAGE', 'Incomplete')
	set(tmp, tmp[, which(is.na(STAGE) & is.na(PosCD4_T1))], 'STAGE', 'Incomplete')
	set(tmp, tmp[, which(is.na(STAGE) & !is.na(PosCD4_T1) & (AnyPos_T1+1<=PosCD4_T1 | PosCD4_T1-1/12>AnyT_T1))], 'STAGE', 'Incomplete')
	set(tmp, tmp[, which(is.na(STAGE))], 'STAGE', 'NoInd')
	set(tmp, NULL, 'STAGE', tmp[, factor(STAGE, levels=c('Acute','NoInd','Late','Incomplete'), labels=c('Confirmed recent infection','No indication for recent infection\nor late presentation','Late presenter','Incomplete data'))])	
	col		<- c("#66C2A5", "#FC8D62", "#8DA0CB", "grey80")  
	ggplot(tmp, aes(x=MIGRANT, fill=STAGE)) + geom_bar(position='fill') +
			theme_bw() + theme(legend.position='bottom') +
			coord_flip() + facet_grid(~Trm) + 
			scale_fill_manual(values=col) +
			scale_y_continuous(labels=percent) +
			labs(x='', y='\nNew HIV diagnoses in Amsterdam\n2010-2014', fill='Infection stage\nat diagnosis') +
			guides(fill=guide_legend(nrow=2,byrow=TRUE))
	ggsave(file=file.path(outdir,'ATHENA1502_NewDiag1014AMST_PropEarlyLateByOrigin1.pdf'), w=9.5, h=4)	
	#
	#	proportion of migrants in Amsterdam diagnosed early
	#
	pea		<- subset(df, !is.na(isAcute) & AMST=='Y' & MIGRANT!='NL' & TP!='<2010' )[, list(PROP_EARLY=mean(isAcute=='Yes')), by="Trm"]
	pead	<- subset(df, !is.na(isAcute) & AMST=='Y' & MIGRANT=='NL' & TP!='<2010' )[, list(PROP_EARLY=mean(isAcute=='Yes')), by="Trm"]
	peo		<- subset(df, !is.na(isAcute) & AMST=='N' & MIGRANT!='NL' & TP!='<2010' )[, list(PROP_EARLY=mean(isAcute=='Yes')), by="Trm"]
	peod	<- subset(df, !is.na(isAcute) & AMST=='N' & MIGRANT=='NL' & TP!='<2010' )[, list(PROP_EARLY=mean(isAcute=='Yes')), by="Trm"]
	#
	#	proportion of migrants in Amsterdam diagnosed late, of those with a first CD4 count in the first year after diagnosis and not on ART before 
	#	with first CD4<350
	#
	pla		<- subset(df, AMST=='Y' & MIGRANT!='NL' & TP!='<2010' & !is.na(PosCD4_T1) & AnyPos_T1+1>PosCD4_T1 & PosCD4_T1<=AnyT_T1)[, list(PROP_LATE=mean(CD4_T1<=350)), by='Trm']
	#
	#	number new diagnoses by region origin and sex,orientation
	#
	tmp			<- subset(df, !is.na(MIGRANT) & AMST=='Y' & Trm%in%c('MSM','HET') & AnyPos_T1>2009 & AnyPos_T1<2015)	
	set(tmp, tmp[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')
	#tmp[, YR:= tmp[, cut(AnyPos_T1, breaks=seq(2001, 2015, 2), labels=c('2001-2002','2003-2004','2005-2006','2007-2008','2009-2010','2011-2012','2013-2014'))]]
	set(tmp, NULL, 'Trm', tmp[, factor(Trm, levels=c('HET','HETM','MSM'), labels=c('Heterosexual women','Heterosexual men', 'MSM'))])	
	set(tmp, NULL, 'MIGRANT', tmp[, factor(as.character(MIGRANT), 	levels=c("NL","Western_EU", "Latin_South_America", "Eastern_EU_stans", "Caribbean", "Sub_Saharan_Africa", "Other"),
							labels=c("Dutch","Western Europe", "Latin & South America", "Eastern Europe", "Caribbean", "Sub-Saharan Africa", "Other"))])
	ggplot(tmp, aes(x=Trm, fill=MIGRANT)) + geom_bar(position='stack') +
			theme_bw() + theme(legend.position='bottom') +
			coord_flip() + 
			scale_y_continuous() +		
			scale_fill_manual(values=c("Dutch"='grey50',"Western Europe"="#D53E4F", "Latin & South America"="#FC8D59", "Eastern Europe"="#FEE08B", "Caribbean"="#E6F598", "Sub-Saharan Africa"="#99D594", "Other"="#3288BD")) +				
			labs(x='', fill='Region of origin', y='\nNew HIV diagnoses among migrants in Amsterdam') 
	ggsave(file=file.path(outdir,'ATHENA1502_NewDiag0114AMST_NOrigin.pdf'), w=8, h=3)
	#
	#	cascade of care
	#
	dc	<- data.table(	TYPE=	c('Undiagnosed','Diagnosed','In Active Follow Up','ART started','Virus suppressed'),
			NC=		c(6100,5669,5266,4992,4693))
	set(dc, NULL, 'TYPE', dc[,factor(TYPE, levels=rev(c('Undiagnosed','Diagnosed','In Active Follow Up','ART started','Virus suppressed')), labels=rev(c('Undiagnosed','Diagnosed','In Active Follow Up','ART started','Virus suppressed')))])			
	dc[, N:= c(NC[-length(NC)]-NC[-1], NC[length(NC)]) ]
	dc[, P:= round(100*N/sum(N),d=1)]
	dc[, POS:=N/2]
	ggplot(dc, aes(x=TYPE, label=paste(P,'%',sep=''))) + geom_bar(aes(y=N), stat='identity', fill='black') + 
			coord_flip() + labs(x='', y='') + theme_bw() +
			geom_text(aes(y=POS), colour='white', size=2) 
	ggsave(file=file.path(outdir,'ATHENA1502_Cascade.pdf'), w=6, h=4)
	#
	#	care cascade results
	#
	dc	<- data.table(TYPE=c('U','D','T','S','L'), C=c(0.71,.224,.057,.016,.01), L=c(0.66,.207,.052,.011,.007), U=c(0.725,.262,.078,.025,.016))
	ggplot(dc, aes(x=TYPE, y=C, ymin=L, ymax=U, fill=TYPE)) + 
			geom_bar(stat='identity') + geom_errorbar(width=.3) + coord_flip() + labs(x='', y='') + 
			scale_y_continuous(labels=percent) + #scale_x_discrete(expand=c(0,0)) +
			#scale_fill_colour() +
			theme_bw()
	ggsave(file=file.path(outdir,'ATHENA1502_Sources.pdf'), w=4, h=3)
	#
	#	number new diagnoses in Amsterdam by Trm and MIGRANT
	#
	#	table
	subset(df, AMST=='Y' & MIGRANT!='NL' & TP!='<2010')[, table(RegionOrigin, interaction(Trm, Sex))]
	#	table by year
	tmp		<- subset(df, AMST=='Y' & Trm=='MSM' & AnyPos_T1>2000)
	tmp[, YR:= tmp[, cut(AnyPos_T1, breaks=seq(2000, 2016, 2))]]
	tmp		<- tmp[, list(DIAG_N=length(Patient)), by=c('YR','RegionOrigin')]	
	dcast.data.table(tmp, YR~RegionOrigin, value.var='DIAG_N')
	#	plot by year
	tmp			<- subset(df, !is.na(MIGRANT) &  AMST=='Y' & Trm%in%c('MSM','HET') & AnyPos_T1>2001 & AnyPos_T1<2015)	
	set(tmp, NULL, 'MIGRANT', tmp[,factor(as.character(MIGRANT)=='NL', levels=c(TRUE,FALSE), labels=c('Dutch origin','Migrant'))])
	set(tmp, tmp[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')
	set(tmp, NULL, 'Trm', tmp[, factor(Trm, levels=c('HET','HETM','MSM'), labels=c('Heterosexual women','Heterosexual men', 'MSM'))])	
	tmp[, YR:= tmp[, cut(AnyPos_T1, breaks=seq(2001, 2015, 2), labels=c('2001-2002','2003-2004','2005-2006','2007-2008','2009-2010','2011-2012','2013-2014'))]]
	tmp			<- tmp[, list(DIAG_N=length(Patient)), by=c('YR','MIGRANT','Trm')]
	col			<- c("#CC4C02","#225EA8","#238443")
	names(col)	<- tmp[, levels(Trm)]
	ggplot( tmp, aes(x=YR, y=DIAG_N, group=interaction(MIGRANT,Trm), fill=Trm, alpha=MIGRANT)) + 
			geom_bar(stat='identity', position='dodge') +  
			scale_fill_manual(values=col) +
			scale_alpha_manual(values=c('Dutch origin'=0.5, 'Migrant'=1)) +			
			labs(x='', y='new HIV diagnoses\nin Amsterdam', fill='Colour:', alpha='Translucency:') + 
			theme_bw() + theme(legend.position='bottom')
	ggsave(file=file.path(outdir,'ATHENA1502_NewDiag0114AMST_Total1.pdf'), w=8, h=5)
	#
	#	number new diagnoses in NL by Trm and MIGRANT
	#
	#	table by sex
	subset(df, MIGRANT!='NL' & TP!='<2010')[, table(MIGRANT, interaction(Trm, Sex))]
	#	table by year
	tmp		<- subset(df, MIGRANT!='NL' & Trm=='MSM' & AnyPos_T1>2000)
	tmp[, YR:= tmp[, cut(AnyPos_T1, breaks=seq(2001, 2017, 2))]]
	tmp		<- tmp[, list(DIAG_N=length(Patient)), by=c('YR','MIGRANT')]	
	dcast.data.table(tmp, YR~MIGRANT, value.var='DIAG_N')
	#	plot by year
	tmp			<- subset(df, !is.na(MIGRANT) &  Trm%in%c('MSM','HET') & AnyPos_T1>2001 & AnyPos_T1<2015)	
	set(tmp, NULL, 'MIGRANT', tmp[,factor(as.character(MIGRANT)=='NL', levels=c(TRUE,FALSE), labels=c('Dutch origin','Migrant'))])
	set(tmp, tmp[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')
	set(tmp, NULL, 'Trm', tmp[, factor(Trm, levels=c('HET','HETM','MSM'), labels=c('Heterosexual women','Heterosexual men', 'MSM'))])	
	tmp[, YR:= tmp[, cut(AnyPos_T1, breaks=seq(2001, 2015, 2), labels=c('2001-2002','2003-2004','2005-2006','2007-2008','2009-2010','2011-2012','2013-2014'))]]
	tmp			<- tmp[, list(DIAG_N=length(Patient)), by=c('YR','MIGRANT','Trm')]
	col			<- c("#CC4C02","#225EA8","#238443")
	names(col)	<- tmp[, levels(Trm)]
	ggplot( tmp, aes(x=YR, y=DIAG_N, group=interaction(MIGRANT,Trm), fill=Trm, alpha=MIGRANT)) + 
			geom_bar(stat='identity', position='dodge') +  
			scale_fill_manual(values=col) +
			scale_alpha_manual(values=c('Dutch origin'=0.5, 'Migrant'=1)) +			
			labs(x='', y='new HIV diagnoses\nin the Netherlands', fill='Colour:', alpha='Translucency:') + 
			theme_bw() + theme(legend.position='bottom')
	ggsave(file=file.path(outdir,'ATHENA1502_NewDiag0114NL_Total1.pdf'), w=8, h=5)
	#
	#	migrant origin 2010-2014 in the Netherlands
	#
	tmp			<- subset(df, !is.na(MIGRANT) & MIGRANT!='NL' & Trm%in%c('MSM','HET') & AnyPos_T1>2001 & AnyPos_T1<2015)	
	set(tmp, tmp[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')
	tmp[, YR:= tmp[, cut(AnyPos_T1, breaks=seq(2001, 2015, 2), labels=c('2001-2002','2003-2004','2005-2006','2007-2008','2009-2010','2011-2012','2013-2014'))]]
	set(tmp, NULL, 'Trm', tmp[, factor(Trm, levels=c('HET','HETM','MSM'), labels=c('Heterosexual women','Heterosexual men', 'MSM'))])	
	set(tmp, NULL, 'MIGRANT', tmp[, factor(as.character(MIGRANT), 	levels=c("Western_EU", "Latin_South_America", "Eastern_EU_stans", "Caribbean", "Sub_Saharan_Africa", "Other"),
							labels=c("Western Europe", "Latin & South America", "Eastern Europe", "Caribbean", "Sub-Saharan Africa", "Other"))])
	ggplot(tmp, aes(x=YR, fill=MIGRANT)) + geom_bar(position='fill') +
			theme_bw() + theme(legend.position='bottom') +
			coord_flip() + facet_grid(~Trm) +
			scale_y_continuous(labels=percent) +
			scale_alpha_manual(values=c('Heterosexual women'=0.4, 'Heterosexual men'=0.7, 'MSM'=1)) +
			scale_fill_brewer(palette='Spectral') +	
			labs(x='', fill='Region of origin', y='\nNew HIV diagnoses among migrants in the Netherlands') 
	ggsave(file=file.path(outdir,'ATHENA1502_NewDiag0114NL_PropOrigin1.pdf'), w=8, h=3.5)
	#
	#	migrant origin 2010-2014 in Amsterdam
	#
	tmp			<- subset(df, !is.na(MIGRANT) & AMST=='Y' & MIGRANT!='NL' & Trm%in%c('MSM','HET') & AnyPos_T1>2001 & AnyPos_T1<2015)	
	set(tmp, tmp[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')
	tmp[, YR:= tmp[, cut(AnyPos_T1, breaks=seq(2001, 2015, 2), labels=c('2001-2002','2003-2004','2005-2006','2007-2008','2009-2010','2011-2012','2013-2014'))]]
	set(tmp, NULL, 'Trm', tmp[, factor(Trm, levels=c('HET','HETM','MSM'), labels=c('Heterosexual women','Heterosexual men', 'MSM'))])	
	set(tmp, NULL, 'MIGRANT', tmp[, factor(as.character(MIGRANT), 	levels=c("Western_EU", "Latin_South_America", "Eastern_EU_stans", "Caribbean", "Sub_Saharan_Africa", "Other"),
							labels=c("Western Europe", "Latin & South America", "Eastern Europe", "Caribbean", "Sub-Saharan Africa", "Other"))])
	ggplot(tmp, aes(x=YR, fill=MIGRANT)) + geom_bar(position='fill') +
			theme_bw() + theme(legend.position='bottom') +
			coord_flip() + facet_grid(~Trm) +
			scale_y_continuous(labels=percent) +
			scale_alpha_manual(values=c('Heterosexual women'=0.4, 'Heterosexual men'=0.7, 'MSM'=1)) +
			scale_fill_brewer(palette='Spectral') +	
			labs(x='', fill='Region of origin', y='\nNew HIV diagnoses among migrants in Amsterdam') 
	ggsave(file=file.path(outdir,'ATHENA1502_NewDiag0114AMST_PropOrigin1.pdf'), w=8, h=3.5)
	#
	#	migrant age 2001-2014 in NL
	#
	require(viridis)
	tmp			<- subset(df, !is.na(MIGRANT) & Trm%in%c('MSM','HET') & AnyPos_T1>2010 & AnyPos_T1<2015)	
	set(tmp, tmp[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')	
	set(tmp, NULL, 'Trm', tmp[, factor(Trm, levels=c('HET','HETM','MSM'), labels=c('Heterosexual women','Heterosexual men', 'MSM'))])		
	set(tmp, NULL, 'AGE', tmp[,cut(Age_AnyPosT1, breaks=c(0,28,38,48,58,100), labels=c('<28 yrs','28-37 yrs','38-47 yrs','48-57 yrs','>57 yrs'))])	
	set(tmp, NULL, 'MIGRANT', tmp[, factor(as.character(MIGRANT), 	levels=rev(c("NL","Western_EU", "Latin_South_America", "Eastern_EU_stans", "Caribbean", "Sub_Saharan_Africa", "Other")),
							labels=rev(c("Dutch origin","Western Europe", "Latin & South America", "Eastern Europe", "Caribbean", "Sub-Saharan Africa", "Other")))])
	ggplot(tmp, aes(x=MIGRANT, fill=AGE)) + geom_bar(position='fill') +
			theme_bw() + theme(legend.position='bottom') +
			coord_flip() + facet_grid(~Trm) +
			scale_y_continuous(labels=percent) +
			scale_alpha_manual(values=c('Heterosexual women'=0.4, 'Heterosexual men'=0.7, 'MSM'=1)) +
			#scale_fill_brewer(palette='Spectral') +	
			scale_fill_viridis(option="viridis", discrete=TRUE) +
			labs(x='', fill='Age at diagnosis among\nnew diagnoses in the Netherlands,   \n2010-2014', y='') 
	ggsave(file=file.path(outdir,'ATHENA1502_NewDiag1014NL_PropAgeByOrigin.pdf'), w=10, h=3.5)
	#
	#	migrant age 2001-2014 in Amsterdam
	#	
	tmp			<- subset(df, AMST=='Y' & !is.na(MIGRANT) & Trm%in%c('MSM','HET') & AnyPos_T1>2010 & AnyPos_T1<2015)	
	set(tmp, tmp[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')	
	set(tmp, NULL, 'Trm', tmp[, factor(Trm, levels=c('HET','HETM','MSM'), labels=c('Heterosexual women','Heterosexual men', 'MSM'))])		
	set(tmp, NULL, 'AGE', tmp[,cut(Age_AnyPosT1, breaks=c(0,28,38,48,58,100), labels=c('<28 yrs','28-37 yrs','38-47 yrs','48-57 yrs','>57 yrs'))])	
	set(tmp, NULL, 'MIGRANT', tmp[, factor(as.character(MIGRANT), 	levels=rev(c("NL","Western_EU", "Latin_South_America", "Eastern_EU_stans", "Caribbean", "Sub_Saharan_Africa", "Other")),
							labels=rev(c("Dutch origin","Western Europe", "Latin & South America", "Eastern Europe", "Caribbean", "Sub-Saharan Africa", "Other")))])
	ggplot(tmp, aes(x=MIGRANT, fill=AGE)) + geom_bar(position='fill') +
			theme_bw() + theme(legend.position='bottom') +
			coord_flip() + facet_grid(~Trm) +
			scale_y_continuous(labels=percent) +
			scale_alpha_manual(values=c('Heterosexual women'=0.4, 'Heterosexual men'=0.7, 'MSM'=1)) +
			#scale_fill_brewer(palette='Spectral') +	
			scale_fill_viridis(option="viridis", discrete=TRUE) +
			labs(x='', fill='Age at diagnosis among\nnew diagnoses in Amsterdam,   \n2010-2014', y='') 
	ggsave(file=file.path(outdir,'ATHENA1502_NewDiag1014AMST_PropAgeByOrigin.pdf'), w=10, h=3.5)
	#
	#	subtype B in the Netherlands
	#
	tmp			<- subset(df, !is.na(Subtype) & !is.na(MIGRANT) & Trm%in%c('MSM','HET') & AnyPos_T1>2010 & AnyPos_T1<2015)	
	set(tmp, tmp[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')
	set(tmp, NULL, 'Subtype', tmp[, factor(Subtype=='B', levels=c(TRUE, FALSE), labels=c('B','non-B'))])	
	set(tmp, NULL, 'Trm', tmp[, factor(Trm, levels=c('HET','HETM','MSM'), labels=c('Heterosexual women','Heterosexual men', 'MSM'))])		
	set(tmp, NULL, 'AGE', tmp[,cut(Age_AnyPosT1, breaks=c(0,28,38,48,58,100), labels=c('<28 yrs','28-37 yrs','38-47 yrs','48-57 yrs','>57 yrs'))])	
	set(tmp, NULL, 'MIGRANT', tmp[, factor(as.character(MIGRANT), 	levels=rev(c("NL","Western_EU", "Latin_South_America", "Eastern_EU_stans", "Caribbean", "Sub_Saharan_Africa", "Other")),
							labels=rev(c("Dutch origin","Western Europe", "Latin & South America", "Eastern Europe", "Caribbean", "Sub-Saharan Africa", "Other")))])
	ggplot(tmp, aes(x=MIGRANT, fill=Subtype)) + geom_bar(position='fill') +
			theme_bw() + theme(legend.position='bottom') +
			coord_flip() + facet_grid(~Trm) +
			scale_y_continuous(labels=percent) +			
			scale_fill_manual(values=c('B'="#66C2A5", 'non-B'="#F46D43")) +	
			#scale_fill_viridis(option="magma", discrete=TRUE) +
			labs(x='', fill='Subtype distribution among\nnew diagnoses in the Netherlands,\n2010-2014', y='') 
	ggsave(file=file.path(outdir,'ATHENA1502_NewDiag1014NL_PropSubtypeByOrigin.pdf'), w=10, h=3.5)
	
	# 	for each time period: 
	#	new diagnoses (all; MSM; MSM/young-old; HET/young-old; Migrant; Migrant/young-old)	
	tmp			<- subset(df, TP!='<2010')
	tmp[, TP:='2010-15']
	tmp			<- rbind(df, tmp)	
	tmp[, {
				z			<- list()
				z[['Reg_NA']]							<- length(which(AMST=='Maybe'))
				z[['Migrant_1YNAM']]					<- length(which(MIGRANT!='NL' & (is.na(isAcute) | isAcute=='Maybe')))
				z[['A']]								<- length(which(AMST=='Y'))
				z[['A_MSM']]							<- length(which(AMST=='Y' & Trm%in%c('MSM','BI')))
				z[['A_MSM_16-27']]						<- length(which(AMST=='Y' & Trm%in%c('MSM','BI') & YOUNG=='16-27'))
				z[['A_MSM_28-80']]						<- length(which(AMST=='Y' & Trm%in%c('MSM','BI') & YOUNG=='28-80'))
				z[['A_MSM_Migrant_16-27']]				<- length(which(AMST=='Y' & MIGRANT!='NL'& Trm%in%c('MSM','BI') & YOUNG=='16-27'))
				z[['A_MSM_Migrant_28-80']]				<- length(which(AMST=='Y' & MIGRANT!='NL'& Trm%in%c('MSM','BI') & YOUNG=='28-80'))
				z[['A_MSM_NL_16-27']]					<- length(which(AMST=='Y' & MIGRANT=='NL'& Trm%in%c('MSM','BI') & YOUNG=='16-27'))
				z[['A_MSM_NL_28-80']]					<- length(which(AMST=='Y' & MIGRANT=='NL'& Trm%in%c('MSM','BI') & YOUNG=='28-80'))				
				z[['A_HET']]							<- length(which(AMST=='Y' & Trm%in%c('HET','HETfa')))
				z[['A_Migrant']]						<- length(which(AMST=='Y' & MIGRANT!='NL'))
				z[['A_Migrant_M']]						<- length(which(AMST=='Y' & MIGRANT!='NL' & Sex=='M'))
				z[['A_Migrant_F']]						<- length(which(AMST=='Y' & MIGRANT!='NL' & Sex=='F'))
				z[['A_Migrant_1YNAM']]					<- length(which(AMST=='Y' & MIGRANT!='NL' & (is.na(isAcute) | isAcute=='Maybe')))
				z[['A_Migrant_MSM']]					<- length(which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL'))				
				z[['A_Migrant_MSM_1Y']]					<- length(which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & isAcute=='Yes'))
				z[['A_Migrant_MSM_1YE']]				<- round( length(which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL')) * subset(pea, Trm=='MSM')[, PROP_EARLY], d=0 )
				z[['A_Migrant_HET']]					<- length(which(AMST=='Y' & Trm%in%c('HET','HETfa') & MIGRANT!='NL'))
				z[['A_Migrant_16-27']]					<- length(which(AMST=='Y' & MIGRANT!='NL' & YOUNG=='16-27'))
				z[['A_Migrant_28-80']]					<- length(which(AMST=='Y' & MIGRANT!='NL' & YOUNG=='28-80'))
				z[['A_NL_16-27']]						<- length(which(AMST=='Y' & MIGRANT=='NL' & YOUNG=='16-27'))
				z[['A_NL_28-80']]						<- length(which(AMST=='Y' & MIGRANT=='NL' & YOUNG=='28-80'))
				z
			}, by='TP']
	#
	#	numbers of sequences and isAcute
	tmp			<- subset(df, TP!='<2010' & TP!='2015')
	tmp[, ADJ:=1]
	tmp2		<- subset(df, TP=='2014')
	tmp2[, ADJ:=0.8]
	tmp2[, TP:='2015']
	tmp			<- rbind(tmp, tmp2)
	tmp2[, ADJ:=0.8]
	tmp2[, TP:='2016']
	tmp			<- rbind(tmp, tmp2)
	tmp			<- tmp[, {
				z										<- list()	
				z[['A_Migrant_Any_Any_MSM']]			<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL')])
				z[['O_Migrant_Any_Any_MSM']]			<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL')])
				z[['A_Dutch_Any_Any_MSM']]				<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
				z[['O_Dutch_Any_Any_MSM']]				<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
				
				z[['A_Migrant_Seq_Any_MSM']]			<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & !is.na(PosSeqT))])
				z[['O_Migrant_Seq_Any_MSM']]			<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & !is.na(PosSeqT))])
				z[['A_Dutch_Seq_Any_MSM']]				<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL' & !is.na(PosSeqT))])
				z[['O_Dutch_Seq_Any_MSM']]				<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL' & !is.na(PosSeqT))])
				
				z[['A_Migrant_Any_Rec_MSM']]			<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & isAcute=='Yes')])
				z[['O_Migrant_Any_Rec_MSM']]			<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & isAcute=='Yes')])
				z[['A_Dutch_Any_Rec_MSM']]				<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL' & isAcute=='Yes')])
				z[['O_Dutch_Any_Rec_MSM']]				<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL' & isAcute=='Yes')])
				
				z[['A_Migrant_Seq_Rec_MSM']]			<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				z[['O_Migrant_Seq_Rec_MSM']]			<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				z[['A_Dutch_Seq_Rec_MSM']]				<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				z[['O_Dutch_Seq_Rec_MSM']]				<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				
				z[['A_Migrant_Any_Any_Any']]			<- sum(ADJ[which(AMST=='Y' & MIGRANT!='NL')])
				z[['O_Migrant_Any_Any_Any']]			<- sum(ADJ[which(AMST=='N' & MIGRANT!='NL')])
				z[['A_Dutch_Any_Any_Any']]				<- sum(ADJ[which(AMST=='Y' & MIGRANT=='NL')])
				z[['O_Dutch_Any_Any_Any']]				<- sum(ADJ[which(AMST=='N' & MIGRANT=='NL')])
				
				z[['A_Migrant_Seq_Any_Any']]			<- sum(ADJ[which(AMST=='Y' & MIGRANT!='NL' & !is.na(PosSeqT))])
				z[['O_Migrant_Seq_Any_Any']]			<- sum(ADJ[which(AMST=='N' & MIGRANT!='NL' & !is.na(PosSeqT))])
				z[['A_Dutch_Seq_Any_Any']]				<- sum(ADJ[which(AMST=='Y' & MIGRANT=='NL' & !is.na(PosSeqT))])
				z[['O_Dutch_Seq_Any_Any']]				<- sum(ADJ[which(AMST=='N' & MIGRANT=='NL' & !is.na(PosSeqT))])
				
				z[['A_Migrant_Any_Rec_Any']]			<- sum(ADJ[which(AMST=='Y' & MIGRANT!='NL' & isAcute=='Yes')])
				z[['O_Migrant_Any_Rec_Any']]			<- sum(ADJ[which(AMST=='N' & MIGRANT!='NL' & isAcute=='Yes')])
				z[['A_Dutch_Any_Rec_Any']]				<- sum(ADJ[which(AMST=='Y' & MIGRANT=='NL' & isAcute=='Yes')])
				z[['O_Dutch_Any_Rec_Any']]				<- sum(ADJ[which(AMST=='N' & MIGRANT=='NL' & isAcute=='Yes')])
				
				
				z[['A_Migrant_Seq_Rec_Any']]			<- sum(ADJ[which(AMST=='Y' & MIGRANT!='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				z[['O_Migrant_Seq_Rec_Any']]			<- sum(ADJ[which(AMST=='N' & MIGRANT!='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				z[['A_Dutch_Seq_Rec_Any']]				<- sum(ADJ[which(AMST=='Y' & MIGRANT=='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				z[['O_Dutch_Seq_Rec_Any']]				<- sum(ADJ[which(AMST=='N' & MIGRANT=='NL' & !is.na(PosSeqT) & isAcute=='Yes')])
				
				z
			}, by='TP']
	cnt			<- melt(tmp, id.vars='TP')
	cnt[, LOC:= factor(substr(variable, 1, 1), levels=c('A','O'), labels=c('Amst','Other'))]
	cnt[, ORIGIN:= cnt[, sapply(strsplit(as.character(variable), '_'), '[[', 2)]]
	cnt[, SEQ:= cnt[, sapply(strsplit(as.character(variable), '_'), '[[', 3)]]	
	cnt[, REC:= cnt[, sapply(strsplit(as.character(variable), '_'), '[[', 4)]]
	cnt[, TRM:= cnt[, sapply(strsplit(as.character(variable), '_'), '[[', 5)]]
	set(cnt, NULL, 'value', cnt[, round(value)])
	set(cnt, NULL, 'variable', NULL)
	cnt			<- rbind(cnt, cnt[, list(TP='2010-2016', value=sum(value)), by=c('LOC','ORIGIN','SEQ','REC','TRM')], use.names=TRUE)	
	
	#	tables for Amsterdam
	tmp			<- merge( 	dcast.data.table(subset(cnt, REC=='Any' & TRM=='Any' & LOC=='Amst'), TP~SEQ+ORIGIN, value.var='value'),
			dcast.data.table(subset(cnt, SEQ=='Any' & TRM=='Any' & LOC=='Amst'), TP~REC+ORIGIN, value.var='value'),
			by=c('TP','Any_Dutch','Any_Migrant'))
	tmp[, Seq_Dutch_P:= round(Seq_Dutch/Any_Dutch, d=2)]
	tmp[, Seq_Migrant_P:= round(Seq_Migrant/Any_Migrant, d=2)]
	tmp[, Rec_Dutch_P:= round(Rec_Dutch/Any_Dutch, d=2)]
	tmp[, Rec_Migrant_P:= round(Rec_Migrant/Any_Migrant, d=2)]
	tmp[, LOC:='Amsterdam']
	ans			<- subset(tmp, select=c(LOC, TP, Any_Dutch, Any_Migrant, Seq_Dutch, Seq_Dutch_P, Seq_Migrant, Seq_Migrant_P, Rec_Dutch, Rec_Dutch_P, Rec_Migrant, Rec_Migrant_P))
	#	tables for Netherlands
	tmp			<- merge( 	dcast.data.table(subset(cnt, REC=='Any' & TRM=='Any' & LOC=='Other'), TP~SEQ+ORIGIN, value.var='value'),
			dcast.data.table(subset(cnt, SEQ=='Any' & TRM=='Any' & LOC=='Other'), TP~REC+ORIGIN, value.var='value'),
			by=c('TP','Any_Dutch','Any_Migrant'))
	tmp[, Seq_Dutch_P:= round(Seq_Dutch/Any_Dutch, d=2)]
	tmp[, Seq_Migrant_P:= round(Seq_Migrant/Any_Migrant, d=2)]
	tmp[, Rec_Dutch_P:= round(Rec_Dutch/Any_Dutch, d=2)]
	tmp[, Rec_Migrant_P:= round(Rec_Migrant/Any_Migrant, d=2)]
	tmp[, LOC:='Other']	
	ans			<- rbind( ans, subset(tmp, select=c(LOC, TP, Any_Dutch, Any_Migrant, Seq_Dutch, Seq_Dutch_P, Seq_Migrant, Seq_Migrant_P, Rec_Dutch, Rec_Dutch_P, Rec_Migrant, Rec_Migrant_P)) )
	
	write.csv(ans, row.names=FALSE, file=file.path(outdir,'ATHENA1502_Diagnoses.csv'))
	
	#
	# 	expected number of recipients and actual transmission pairs captured, 2010-2016:
	#	baseline:
	if(0)
	{
		m.A			<- 0.95 * 0.9*(0.7*.95+.3*.35)
		m.O			<- 0.35 * 0.9*(0.3*.95+.7*.35)		
	}
	#	rather than RITA, let s get 50% coverage outside Amsterdam
	if(0)
	{
		m.A			<- 0.95 * 0.9*(0.7*.95+.3*.5)
		m.O			<- 0.5 * 0.9*(0.3*.95+.7*.5)		
	}
	#	for young, let s get 75% coverage outside Amsterdam
	if(1)
	{
		m.A			<- 0.95 * 0.9*(0.7*.95+.3*.75)
		m.O			<- 0.75 * 0.9*(0.3*.95+.7*.75)		
	}
	tmp			<- subset(df, TP!='<2010' & TP!='2015')
	tmp[, ADJ:=1]
	tmp2		<- subset(df, TP=='2014')
	tmp2[, ADJ:=0.8]
	tmp2[, TP:='2015']
	tmp			<- rbind(tmp, tmp2)
	tmp2[, ADJ:=0.8]
	tmp2[, TP:='2016']
	tmp			<- rbind(tmp, tmp2)
	tmp			<- tmp[, {
				z			<- list()												
				z[['A_NL_MSM_ALL_1YE']]				<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(pead, Trm=='MSM')[, PROP_EARLY]
				z[['A_NL_HETM_ALL_1YE']]			<- sum(ADJ[which(AMST=='Y' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]
				z[['A_NL_HETF_ALL_1YE']]			<- sum(ADJ[which(AMST=='Y' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]				
				z[['A_Migrant_ALL_MSM_1YE']]		<- sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(pea, Trm=='MSM')[, PROP_EARLY]
				z[['A_Migrant_HETM_ALL_1YE']]		<- sum(ADJ[which(AMST=='Y' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
				z[['A_Migrant_HETF_ALL_1YE']]		<- sum(ADJ[which(AMST=='Y' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
				
				z[['A_NL_MSM_ALL_1Y']]				<- sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
				z[['A_NL_HETM_ALL_1Y']]				<- sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
				z[['A_NL_HETF_ALL_1Y']]				<- sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
				z[['A_Migrant_MSM_ALL_1Y']]			<- sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
				z[['A_Migrant_HETM_ALL_1Y']]		<- sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				z[['A_Migrant_HETF_ALL_1Y']]		<- sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				
				z[['A_NL_MSM_28_1YE']]				<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(pead, Trm=='MSM')[, PROP_EARLY]
				z[['A_NL_HETM_28_1YE']]				<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]
				z[['A_NL_HETF_28_1YE']]				<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]				
				z[['A_Migrant_28_MSM_1YE']]			<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(pea, Trm=='MSM')[, PROP_EARLY]
				z[['A_Migrant_HETM_28_1YE']]		<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
				z[['A_Migrant_HETF_28_1YE']]		<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
				
				z[['A_NL_MSM_28_1Y']]				<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
				z[['A_NL_HETM_28_1Y']]				<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
				z[['A_NL_HETF_28_1Y']]				<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
				z[['A_Migrant_MSM_28_1Y']]			<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
				z[['A_Migrant_HETM_28_1Y']]			<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				z[['A_Migrant_HETF_28_1Y']]			<- sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				
				if(1)	#since recently infected, in-country infection more likely
				{
					z[['A_NL_MSM_ALL_PE']]			<- m.A*sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(pead, Trm=='MSM')[, PROP_EARLY]
					z[['A_NL_HETM_ALL_PE']]			<- m.A*sum(ADJ[which(AMST=='Y' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]
					z[['A_NL_HETF_ALL_PE']]			<- m.A*sum(ADJ[which(AMST=='Y' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]				
					z[['A_Migrant_MSM_ALL_PE']]		<- m.A*sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(pea, Trm=='MSM')[, PROP_EARLY]
					z[['A_Migrant_HETM_ALL_PE']]	<- m.A*sum(ADJ[which(AMST=='Y' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					z[['A_Migrant_HETF_ALL_PE']]	<- m.A*sum(ADJ[which(AMST=='Y' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					
					z[['A_NL_MSM_ALL_P']]			<- m.A*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) 
					z[['A_NL_HETM_ALL_P']]			<- m.A*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
					z[['A_NL_HETF_ALL_P']]			<- m.A*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
					z[['A_Migrant_MSM_ALL_P']]		<- m.A*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETM_ALL_P']]		<- m.A*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETF_ALL_P']]		<- m.A*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					
					z[['A_NL_MSM_28_PE']]			<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(pead, Trm=='MSM')[, PROP_EARLY]
					z[['A_NL_HETM_28_PE']]			<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]
					z[['A_NL_HETF_28_PE']]			<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]				
					z[['A_Migrant_MSM_28_PE']]		<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(pea, Trm=='MSM')[, PROP_EARLY]
					z[['A_Migrant_HETM_28_PE']]		<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					z[['A_Migrant_HETF_28_PE']]		<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					
					z[['A_NL_MSM_28_P']]			<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) 
					z[['A_NL_HETM_28_P']]			<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
					z[['A_NL_HETF_28_P']]			<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
					z[['A_Migrant_MSM_28_P']]		<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETM_28_P']]		<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETF_28_P']]		<- m.A*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
				}
				if(0)	#prob in-country infection from aMASE
				{
					z[['A_NL_MSM_ALL_PE']]			<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(pead, Trm=='MSM')[, PROP_EARLY]
					z[['A_NL_HETM_ALL_PE']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]
					z[['A_NL_HETF_ALL_PE']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]				
					z[['A_Migrant_MSM_ALL_PE']]		<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(pea, Trm=='MSM')[, PROP_EARLY]
					z[['A_Migrant_HETM_ALL_PE']]	<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					z[['A_Migrant_HETF_ALL_PE']]	<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					
					z[['A_NL_MSM_ALL_P']]			<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) 
					z[['A_NL_HETM_ALL_P']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
					z[['A_NL_HETF_ALL_P']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
					z[['A_Migrant_MSM_ALL_P']]		<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETM_ALL_P']]		<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETF_ALL_P']]		<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					
					z[['A_NL_MSM_28_PE']]			<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(pead, Trm=='MSM')[, PROP_EARLY]
					z[['A_NL_HETM_28_PE']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]
					z[['A_NL_HETF_28_PE']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(pead, Trm=='HET')[, PROP_EARLY]				
					z[['A_Migrant_MSM_28_PE']]		<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(pea, Trm=='MSM')[, PROP_EARLY]
					z[['A_Migrant_HETM_28_PE']]		<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					z[['A_Migrant_HETF_28_PE']]		<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(pea, Trm=='HET')[, PROP_EARLY]
					
					z[['A_NL_MSM_28_P']]			<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) 
					z[['A_NL_HETM_28_P']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
					z[['A_NL_HETF_28_P']]			<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
					z[['A_Migrant_MSM_28_P']]		<- 0.95*0.54*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETM_28_P']]		<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
					z[['A_Migrant_HETF_28_P']]		<- 0.95*0.36*sum(ADJ[which(AMST=='Y' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
				}				
				#z[['O']]						<- sum(ADJ[which(AMST=='N')])				
				z[['O_NL_MSM_ALL_1YE']]				<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(peod, Trm=='MSM')[, PROP_EARLY]
				z[['O_NL_HETM_ALL_1YE']]			<- sum(ADJ[which(AMST=='N' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]
				z[['O_NL_HETF_ALL_1YE']]			<- sum(ADJ[which(AMST=='N' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]								
				z[['O_Migrant_MSM_ALL_1YE']]		<- sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(peo, Trm=='MSM')[, PROP_EARLY]
				z[['O_Migrant_HETM_ALL_1YE']]		<- sum(ADJ[which(AMST=='N' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
				z[['O_Migrant_HETF_ALL_1YE']]		<- sum(ADJ[which(AMST=='N' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
				
				z[['O_NL_MSM_ALL_1Y']]				<- sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
				z[['O_NL_HETM_ALL_1Y']]				<- sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
				z[['O_NL_HETF_ALL_1Y']]				<- sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
				z[['O_Migrant_MSM_ALL_1Y']]			<- sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
				z[['O_Migrant_HETM_ALL_1Y']]		<- sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				z[['O_Migrant_HETF_ALL_1Y']]		<- sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				
				z[['O_NL_MSM_28_1YE']]				<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(peod, Trm=='MSM')[, PROP_EARLY]
				z[['O_NL_HETM_28_1YE']]				<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]
				z[['O_NL_HETF_28_1YE']]				<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]								
				z[['O_Migrant_MSM_28_1YE']]			<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(peo, Trm=='MSM')[, PROP_EARLY]
				z[['O_Migrant_HETM_28_1YE']]		<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
				z[['O_Migrant_HETF_28_1YE']]		<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
				
				z[['O_NL_MSM_28_1Y']]				<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
				z[['O_NL_HETM_28_1Y']]				<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 
				z[['O_NL_HETF_28_1Y']]				<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) 				
				z[['O_Migrant_MSM_28_1Y']]			<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) 
				z[['O_Migrant_HETM_28_1Y']]			<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				z[['O_Migrant_HETF_28_1Y']]			<- sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) 
				
				if(1)	#since recently infected, in-country infection more likely
				{
					z[['O_NL_MSM_ALL_PE']]			<- m.O*sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(peod, Trm=='MSM')[, PROP_EARLY]
					z[['O_NL_HETM_ALL_PE']]			<- m.O*sum(ADJ[which(AMST=='N' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]
					z[['O_NL_HETF_ALL_PE']]			<- m.O*sum(ADJ[which(AMST=='N' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]								
					z[['O_Migrant_MSM_ALL_PE']]		<- m.O*sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(peo, Trm=='MSM')[, PROP_EARLY]
					z[['O_Migrant_HETM_ALL_PE']]	<- m.O*sum(ADJ[which(AMST=='N' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					z[['O_Migrant_HETF_ALL_PE']]	<- m.O*sum(ADJ[which(AMST=='N' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					
					z[['O_NL_MSM_ALL_P']]			<- m.O*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
					z[['O_NL_HETM_ALL_P']]			<- m.O*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])
					z[['O_NL_HETF_ALL_P']]			<- m.O*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])								
					z[['O_Migrant_MSM_ALL_P']]		<- m.O*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')])
					z[['O_Migrant_HETM_ALL_P']]		<- m.O*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					z[['O_Migrant_HETF_ALL_P']]		<- m.O*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					
					z[['O_NL_MSM_28_PE']]			<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(peod, Trm=='MSM')[, PROP_EARLY]
					z[['O_NL_HETM_28_PE']]			<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]
					z[['O_NL_HETF_28_PE']]			<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]								
					z[['O_Migrant_MSM_28_PE']]		<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(peo, Trm=='MSM')[, PROP_EARLY]
					z[['O_Migrant_HETM_28_PE']]		<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					z[['O_Migrant_HETF_28_PE']]		<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					
					z[['O_NL_MSM_28_P']]			<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
					z[['O_NL_HETM_28_P']]			<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])
					z[['O_NL_HETF_28_P']]			<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])								
					z[['O_Migrant_MSM_28_P']]		<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')])
					z[['O_Migrant_HETM_28_P']]		<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					z[['O_Migrant_HETF_28_P']]		<- m.O*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
				}				
				if(0)	#prob in-country infection from aMASE
				{
					z[['O_NL_MSM_ALL_PE']]			<- 0.5*0.54*sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(peod, Trm=='MSM')[, PROP_EARLY]
					z[['O_NL_HETM_ALL_PE']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]
					z[['O_NL_HETF_ALL_PE']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]								
					z[['O_Migrant_MSM_ALL_PE']]		<- 0.5*0.54*sum(ADJ[which(AMST=='N' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(peo, Trm=='MSM')[, PROP_EARLY]
					z[['O_Migrant_HETM_ALL_PE']]	<- 0.5*0.36*sum(ADJ[which(AMST=='N' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					z[['O_Migrant_HETF_ALL_PE']]	<- 0.5*0.36*sum(ADJ[which(AMST=='N' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					
					z[['O_NL_MSM_ALL_P']]			<- 0.5*0.54*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
					z[['O_NL_HETM_ALL_P']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])
					z[['O_NL_HETF_ALL_P']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])								
					z[['O_Migrant_MSM_ALL_P']]		<- 0.5*0.54*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')])
					z[['O_Migrant_HETM_ALL_P']]		<- 0.5*0.36*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					z[['O_Migrant_HETF_ALL_P']]		<- 0.5*0.36*sum(ADJ[which(AMST=='N' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					
					z[['O_NL_MSM_28_PE']]			<- 0.5*0.54*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT=='NL')]) * subset(peod, Trm=='MSM')[, PROP_EARLY]
					z[['O_NL_HETM_28_PE']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]
					z[['O_NL_HETF_28_PE']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')]) * subset(peod, Trm=='HET')[, PROP_EARLY]								
					z[['O_Migrant_MSM_28_PE']]		<- 0.5*0.54*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Trm%in%c('MSM','BI') & MIGRANT!='NL')]) * subset(peo, Trm=='MSM')[, PROP_EARLY]
					z[['O_Migrant_HETM_28_PE']]		<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					z[['O_Migrant_HETF_28_PE']]		<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')]) * subset(peo, Trm=='HET')[, PROP_EARLY]
					
					z[['O_NL_MSM_28_P']]			<- 0.5*0.54*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT=='NL')])
					z[['O_NL_HETM_28_P']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])
					z[['O_NL_HETF_28_P']]			<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT=='NL')])								
					z[['O_Migrant_MSM_28_P']]		<- 0.5*0.54*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Trm%in%c('MSM','BI') & MIGRANT!='NL')])
					z[['O_Migrant_HETM_28_P']]		<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='M' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])
					z[['O_Migrant_HETF_28_P']]		<- 0.5*0.36*sum(ADJ[which(AMST=='N' & YOUNG=='16-27' & isAcute=='Yes' & Sex=='F' & Trm%in%c('HET','HETfa') & MIGRANT!='NL')])					
				}				
				z
			}, by='TP']
	cnt			<- melt(tmp, id.vars='TP')
	cnt[, LOC:= factor(substr(variable, 1, 1), levels=c('A','O'), labels=c('Amst','Other'))]
	cnt[, WHAT:= regmatches(variable,regexpr('[^_]*$', variable))]
	cnt[, add_RITA:= factor(grepl('1YE|PE', WHAT), levels=c(TRUE,FALSE), labels=c('Y','N'))]
	set(cnt, NULL, 'WHAT', cnt[, factor(grepl('1Y',WHAT), levels=c(TRUE,FALSE),labels=c('recent', 'pair'))])	
	set(cnt, NULL, 'YOUNG', cnt[, factor(grepl('28',variable), levels=c(TRUE,FALSE),labels=c('<28', 'all'))])	
	cnt[, TRM:= cnt[, sapply(strsplit(as.character(variable), '_'), '[[', 3)]]
	cnt[, MIGRANT:= cnt[, factor(sapply(strsplit(as.character(variable), '_'), '[[', 2), levels=c('Migrant','NL'), labels=c('Y','N'))]]	
	cnt[, SC_RITA_1016:= add_RITA]
	cnt[, SC_RITA_1315:= factor( (TP%in%c('2010','2011','2012','2016') & add_RITA=='N') |  (TP%in%c('2013','2014','2015') & add_RITA=='Y'), levels=c(TRUE,FALSE), labels=c('Y','N'))]
	cnt			<- rbind( 	subset(cnt,SC_RITA_1016=='Y')[, list(SC='RITA_1016', N=round(sum(value))), by=c('WHAT','YOUNG','MIGRANT','TRM','LOC')],
			subset(cnt,SC_RITA_1315=='Y')[, list(SC='RITA_1315', N=round(sum(value))), by=c('WHAT','YOUNG','MIGRANT','TRM','LOC')]	)
	setkey(cnt, SC, WHAT, MIGRANT, TRM, LOC)				
	dcast.data.table( subset(cnt, WHAT=='pair' & YOUNG=='all' & SC=='RITA_1016'), MIGRANT+LOC~TRM, value.var='N')
	dcast.data.table( subset(cnt, WHAT=='pair' & YOUNG=='all' & SC=='RITA_1315'), MIGRANT+LOC~TRM, value.var='N')
	dcast.data.table( subset(cnt, WHAT=='pair' & YOUNG=='<28' & SC=='RITA_1016'), MIGRANT+LOC~TRM, value.var='N')
	dcast.data.table( subset(cnt, WHAT=='pair' & YOUNG=='<28' & SC=='RITA_1315'), MIGRANT+LOC~TRM, value.var='N')
}
######################################################################################
eval.diag.newdiagnoses.by.migration.region<- function()
{
	require(viridis)
	infile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	outfile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	load(infile)
	#
	#	age group
	#	
	dfs		<- df.all[, {
				z	<- NA_character_
				zz	<- NA_real_
				if(any(!is.na(PosSeqT)))
				{
					z	<- FASTASampleCode[which.min(PosSeqT)]
					zz	<- min(PosSeqT, na.rm=TRUE)
				}
				list(FASTASampleCode=z, PosSeqT=zz)
			}, by=c('Patient')]
	tmp		<- copy(df.all)
	set(tmp, NULL, c('PosSeqT','FASTASampleCode','SUBTYPE_C','SEQ_L'), NULL)
	dfs		<- merge(dfs, unique(tmp), by='Patient')
	
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG','SXCH'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unkown')
	set(dfs, dfs[, which(is.na(Region_first))], 'Region_first', 'Unknown')
	set(dfs, dfs[, which(Region_first%in%c('North','South','East'))], 'Region_first', 'Other')
	set(dfs, dfs[, which(RegionOrigin%in%c("Central_EU","Eastern_EU_stans"))], "RegionOrigin", "Cen_East_Stans_EU")
	set(dfs, dfs[, which(RegionOrigin%in%c("Austr_NZ","Sout_SouthEast_Asia","Oceania_Pacific","North_Africa_Middle_East"))], "RegionOrigin", "Other")
	set(dfs, dfs[, which(is.na(RegionOrigin))], "RegionOrigin", "Unknown")
	set(dfs, NULL, "RegionOrigin", dfs[, as.character(factor(RegionOrigin))])
	set(dfs, NULL, 'Age_AnyPosT1', dfs[, AnyPos_T1-DateBorn])
	set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, breaks=c(-.01,15,seq(28,48,5),100), labels=c('0-15','16-27','28-32','33-37','38-42','43-47','48-80'))])
	set(dfs, NULL, 'YOUNG', dfs[, cut(Age_AnyPosT1, breaks=c(-.01,15,28,100), labels=c('0-15', '16-27','28-80'))])
	set(dfs, NULL, 'OLD', dfs[, cut(Age_AnyPosT1, breaks=c(-.01,48,100), labels=c('0-47', '48-80'))])
	#
	#
	#
	ggplot(subset(dfs, Region_first!='Unknown'), aes(x=floor(AnyPos_T1), fill=RegionOrigin)) +
			geom_bar() +
			scale_x_continuous(breaks=seq(1980, 2020, 5)) +
			#scale_y_continuous(breaks=seq(0,1,.1), labels=percent) +
			scale_fill_brewer(palette='Spectral') +
			labs(x='\nyear of diagnosis', y='new HIV+ diagnoses\n') +
			facet_grid(Region_first~Trm) +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiag_by_region_trm_migrant.pdf'), w=15, h=15)
}
######################################################################################
eval.diag.newdiagnoses.by.age<- function()
{
	require(viridis)
	infile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	outfile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	load(infile)
	
	
	#	first sequences
	dfs		<- df.all[, {
				z	<- NA_character_
				zz	<- NA_real_
				if(any(!is.na(PosSeqT)))
				{
					z	<- FASTASampleCode[which.min(PosSeqT)]
					zz	<- min(PosSeqT, na.rm=TRUE)
				}
				list(FASTASampleCode=z, PosSeqT=zz)
			}, by=c('Patient')]
	tmp		<- copy(df.all)
	set(tmp, NULL, c('PosSeqT','FASTASampleCode','SUBTYPE_C','SEQ_L'), NULL)
	dfs		<- merge(dfs, unique(tmp), by='Patient')
	#	redefine factors
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG','SXCH'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unkown')	
	#set(dfs, dfs[, which(Region_first%in%c('North','South','East'))], 'Region_first', 'Other')
	set(dfs, dfs[, which(RegionOrigin%in%c("Central_EU","Eastern_EU_stans"))], "RegionOrigin", "Cen_East_Stans_EU")
	set(dfs, dfs[, which(RegionOrigin%in%c("Austr_NZ","Sout_SouthEast_Asia","Oceania_Pacific","North_Africa_Middle_East"))], "RegionOrigin", "Other")
	set(dfs, dfs[, which(is.na(RegionOrigin))], "RegionOrigin", "Unknown")
	set(dfs, NULL, "RegionOrigin", dfs[, as.character(factor(RegionOrigin))])
	set(dfs, NULL, 'Age_AnyPosT1', dfs[, AnyPos_T1-DateBorn])
	set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,seq(15,65,10),100), labels=c('0-14','15-24','25-34','35-44','45-54','55-64','65-99'))])
	set(dfs, NULL, 'AGE', dfs[, factor(as.character(AGE), levels=rev(levels(dfs$AGE)))])
	set(dfs, NULL, 'YOUNG', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,15,30,100), labels=c('0-14', '15-29','30-99'))])
	set(dfs, NULL, 'OLD', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,15,50,70,100), labels=c('0-14', '15-49', '50-69', '70-99'))])	
	#
	#	age
	#
	ggplot(subset(dfs, Region_first!='Unknown'), aes(x=floor(AnyPos_T1), fill=AGE)) +
			geom_bar() +
			scale_x_continuous(breaks=seq(1980, 2020, 5)) +
			scale_fill_viridis(option="magma", discrete=TRUE, end=0.85) +
			#scale_y_continuous(breaks=seq(0,1,.1), labels=percent) +			
			labs(x='\nyear of diagnosis', y='new HIV+ diagnoses\n', fill='age group') +
			facet_grid(Region_first~Trm, scales='free_y') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiag_by_age_count.pdf'), w=15, h=15)
	ggplot(subset(dfs, Region_first!='Unknown' & !Trm%in%c('IDU','OTH','Unkown')), aes(x=floor(AnyPos_T1), fill=AGE)) +
			geom_bar(position='fill') +
			scale_x_continuous(breaks=seq(1980, 2020, 5)) +
			scale_fill_viridis(option="magma", discrete=TRUE, end=0.85) +
			scale_y_continuous(breaks=seq(0,1,.1), labels=percent) +			
			labs(x='\nyear of diagnosis', y='new HIV+ diagnoses\n', fill='age group') +
			facet_grid(Trm~Region_first, scales='free_y') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiag_by_age_percent.pdf'), w=25, h=8)
	#
	#	young
	#
	ggplot(subset(dfs, Region_first!='Unknown'), aes(x=floor(AnyPos_T1), fill=YOUNG)) +
			geom_bar() +
			scale_x_continuous(breaks=seq(1980, 2020, 5)) +
			scale_fill_viridis(option="magma", discrete=TRUE, end=0.85) +
			#scale_y_continuous(breaks=seq(0,1,.1), labels=percent) +			
			labs(x='\nyear of diagnosis', y='new HIV+ diagnoses\n', fill='age group') +
			facet_grid(Region_first~Trm, scales='free_y') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiag_by_young_count.pdf'), w=15, h=15)
	ggplot(subset(dfs, Region_first!='Unknown' & !Trm%in%c('IDU','OTH','Unkown')), aes(x=floor(AnyPos_T1), fill=YOUNG)) +
			geom_bar(position='fill') +
			scale_x_continuous(breaks=seq(1980, 2020, 5)) +
			scale_fill_viridis(option="magma", discrete=TRUE, end=0.85) +
			scale_y_continuous(breaks=seq(0,1,.1), labels=percent) +			
			labs(x='\nyear of diagnosis', y='new HIV+ diagnoses\n', fill='age group') +
			facet_grid(Trm~Region_first, scales='free_y') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiag_by_young_percent.pdf'), w=25, h=8)
	#
	#	old
	#
	ggplot(subset(dfs, Region_first!='Unknown'), aes(x=floor(AnyPos_T1), fill=OLD)) +
			geom_bar() +
			scale_x_continuous(breaks=seq(1980, 2020, 5)) +
			scale_fill_viridis(option="magma", discrete=TRUE, end=0.85) +
			#scale_y_continuous(breaks=seq(0,1,.1), labels=percent) +			
			labs(x='\nyear of diagnosis', y='new HIV+ diagnoses\n', fill='age group') +
			facet_grid(Region_first~Trm, scales='free_y') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiag_by_old_count.pdf'), w=15, h=15)
	ggplot(subset(dfs, Region_first!='Unknown' & !Trm%in%c('IDU','OTH','Unkown')), aes(x=floor(AnyPos_T1), fill=OLD)) +
			geom_bar(position='fill') +
			scale_x_continuous(breaks=seq(1980, 2020, 5)) +
			scale_fill_viridis(option="magma", discrete=TRUE, end=0.85) +
			scale_y_continuous(breaks=seq(0,1,.1), labels=percent) +			
			labs(x='\nyear of diagnosis', y='new HIV+ diagnoses\n', fill='age group') +
			facet_grid(Trm~Region_first, scales='free_y') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiag_by_old_percent.pdf'), w=25, h=8)
}
######################################################################################
eval.diag.inrecentinfection<- function()
{
	require(viridis)
	infile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	outfile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	load(infile)
	
	dfs		<- df.all[, {
				z	<- NA_character_
				zz	<- NA_real_
				if(any(!is.na(PosSeqT)))
				{
					z	<- FASTASampleCode[which.min(PosSeqT)]
					zz	<- min(PosSeqT, na.rm=TRUE)
				}
				list(FASTASampleCode=z, PosSeqT=zz)
			}, by=c('Patient')]
	tmp		<- copy(df.all)
	set(tmp, NULL, c('PosSeqT','FASTASampleCode','SUBTYPE_C','SEQ_L'), NULL)
	dfs		<- merge(dfs, unique(tmp), by='Patient')	
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG','SXCH'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unkown')
	set(dfs, dfs[, which(Trm%in%c('Het') & Sex=='M')], 'Trm', 'male heterosexual')
	set(dfs, dfs[, which(Trm%in%c('Het') & Sex=='F')], 'Trm', 'female heterosexual')
	set(dfs, dfs[, which(is.na(Region_now))], 'Region_now', 'Unkown')
	#set(dfs, dfs[, which(Region_now=='Rott')], 'Region_now', 'Other')	
	set(dfs, NULL, 'ACUTE_L', NA_character_)
	set(dfs, dfs[, which(is.na(Acute_Spec))],'ACUTE_L', 'missing')
	set(dfs, dfs[, which(isAcute=='No')],'ACUTE_L', 'not indicated')
	set(dfs, dfs[, which(isAcute=='Yes')],'ACUTE_L', 'yes')
	set(dfs, dfs[, which(isAcute=='Maybe')],'ACUTE_L', 'not indicated')
	set(dfs, NULL, 'AnyPosT1_Y', dfs[, floor(AnyPos_T1)])
		
	ggplot(subset(dfs, AnyPos_T1>2000 & (Trm=='MSM' | grepl('heterosexual',Trm))), aes(x=AnyPos_T1, fill=ACUTE_L)) + 			
			geom_bar(binwidth=1, position='fill', alpha=0.8) +
			scale_y_continuous(breaks=seq(0,1,0.2),labels=percent, expand=c(0,0)) +
			scale_x_continuous(breaks=seq(1980,2020,2), minor_breaks=NULL, expand=c(0,0)) +
			scale_fill_viridis(option="viridis", discrete=TRUE) +			
			labs(x='', y='Infection status at diagnosis\n( % )',fill='in recent infection at diagnosis\n(< 1yr)') + 
			facet_grid(Trm~.) +					
			theme_bw() +
			theme(legend.position='bottom') +	
			guides(fill=guide_legend(ncol=2))
	ggsave(file=paste(outfile,'newdiag_by_recentinfection.pdf'), w=10, h=8)	
	
	ggplot(subset(dfs, AnyPos_T1>2000 & (Trm=='MSM' | grepl('heterosexual',Trm))), aes(x=AnyPos_T1, fill=ACUTE_L)) + 			
			geom_bar(binwidth=1, position='fill', alpha=0.8) +
			scale_y_continuous(breaks=seq(0,1,0.2),labels=percent, expand=c(0,0)) +
			scale_x_continuous(breaks=seq(1980,2020,2), minor_breaks=NULL, expand=c(0,0)) +
			scale_fill_viridis(option="viridis", discrete=TRUE) +			
			labs(x='', y='Infection status at diagnosis\n( % )',fill='in recent infection at diagnosis\n(< 1yr)') + 
			facet_grid(Trm~Region_now) +					
			theme_bw() +
			theme(legend.position='bottom') +	
			guides(fill=guide_legend(ncol=2))
	ggsave(file=paste(outfile,'newdiag_by_recentinfection_region.pdf'), w=12, h=12)
}
######################################################################################
newdiag.characteristics.table<- function(dc)
{
	#	% last neg test in year before diagnosis
	dci	<- dc[, {
				z<- table(!is.na(NegT) & (AnyPos_T1-NegT)<=1)
				list(STAT='Last_Neg_Annual_Test', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])
			}, by='GROUP']
	#	% in recent infection at diagnosis
	tmp	<- dc[, {
				z<- table(isAcute[!is.na(isAcute)]=='Yes')
				list(STAT='In_Recent_When_Not_Missing', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])
			}, by='GROUP']
	dci	<- rbind(dci, tmp)
	#	CD4>500 at diagnosis when immediately available (at most 6mo after diagnosis) and before ART (grace 1 month)
	tmp	<- dc[, {
				z<- table(CD4_T1[(PosCD4_T1-AnyPos_T1)<.5 & PosCD4_T1<=(AnyT_T1+(1/12))]>500)
				list(STAT='CD4g500', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])
			}, by='GROUP']
	dci	<- rbind(dci, tmp)
	#	region origin
	tmp	<- dc[, {
				z<- table(RegionOrigin)
				list(STAT='RegionOrigin', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])
			}, by='GROUP']
	dci	<- rbind(dci, tmp)
	#	ever sequenced
	tmp	<- dc[, {
				z<- table(!is.na(PosSeqT))
				list(STAT='With_Sequence', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])
			}, by='GROUP']
	dci	<- rbind(dci, tmp)
	#	self reported infection abroad
	tmp	<- dc[, {
				z<- table(CountryInfection[!is.na(CountryInfection)]!='NL')
				list(STAT='Slf_Infection_Abroad', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])
			}, by='GROUP']
	dci	<- rbind(dci, tmp)
	#	add proportions
	set(dci, dci[, which(is.na(N))], 'N', 0)
	tmp	<- dci[, {
				z	<- binconf(N, sum(N))
				list(LABEL=LABEL, PROP=round(100*z[,1], d=1), PROP_L=round(100*z[,2], d=1), PROP_U=round(100*z[,3], d=1))	
			}, by=c('GROUP','STAT')]
	dci	<- merge(dci, tmp, by=c('GROUP','STAT','LABEL'))
	dci[, PROP_CI:= paste('(',PROP_L,'-',PROP_U,')',sep='')]
	#	two columns
	set(dci, NULL, 'GROUP', dci[, paste('Age',GROUP,'PC',sep='_')])
	tmp	<- dcast.data.table(dci, STAT+LABEL~GROUP, value.var='PROP')
	set(dci, NULL, 'GROUP', dci[, gsub('_PC','_CI',GROUP)])
	tmp	<- merge(tmp, dcast.data.table(dci, STAT+LABEL~GROUP, value.var='PROP_CI'), by=c('STAT','LABEL'))
	set(dci, NULL, 'GROUP', dci[, gsub('_CI','_N',GROUP)])
	dci	<- dcast.data.table(dci, STAT+LABEL~GROUP, value.var='N')
	dci	<- merge(dci, tmp, by=c('STAT','LABEL'))
	dci
}
######################################################################################
eval.diag.characteristics.by.age<- function()
{
	require(viridis)
	require(scales)
	require(Hmisc)
	
	infile.hiv	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	infile.pop	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/CBS_1612_Population_GGD_Age.rda'
	outfile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	load(infile.hiv)
	#
	#	reduce to patients with / without first sequence
	#	
	dfs		<- df.all[, {
				z	<- NA_character_
				zz	<- NA_real_
				if(any(!is.na(PosSeqT)))
				{
					z	<- FASTASampleCode[which.min(PosSeqT)]
					zz	<- min(PosSeqT, na.rm=TRUE)
				}
				list(FASTASampleCode=z, PosSeqT=zz)
			}, by=c('Patient')]
	tmp		<- copy(df.all)
	set(tmp, NULL, c('PosSeqT','FASTASampleCode','SUBTYPE_C','SEQ_L'), NULL)
	dfs		<- merge(dfs, unique(tmp), by='Patient')	
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG','SXCH'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(Sex=='M' & Trm=='Het')], 'Trm', 'HetM')
	set(dfs, dfs[, which(Sex=='F' & Trm=='Het')], 'Trm', 'HetF')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unknown')
	set(dfs, dfs[, which(RegionOrigin%in%c("Central_EU","Eastern_EU_stans"))], "RegionOrigin", "Cen_East_Stans_EU")
	set(dfs, dfs[, which(RegionOrigin%in%c("Austr_NZ","Sout_SouthEast_Asia","Oceania_Pacific","North_Africa_Middle_East"))], "RegionOrigin", "Other")
	set(dfs, dfs[, which(is.na(RegionOrigin))], "RegionOrigin", "Unknown")
	set(dfs, NULL, "RegionOrigin", dfs[, as.character(factor(RegionOrigin))])
	set(dfs, NULL, 'Age_AnyPosT1', dfs[, AnyPos_T1-DateBorn])
	set(dfs, NULL, 'YEAR', dfs[, floor(AnyPos_T1)])		
	set(dfs, dfs[, which(Region_first=='Amst')], 'Region_first', 'North')
	set(dfs, NULL, 'GGD_first', dfs[, gsub('Hulpverlening_','',as.character(GGD_first))])
	set(dfs, dfs[, which(is.na(GGD_first))], "GGD_first", "Unknown")
	set(dfs, NULL, 'GGD_now', dfs[, gsub('Hulpverlening_','',as.character(GGD_now))])
	set(dfs, dfs[, which(is.na(GGD_now))], "GGD_now", "Unknown")		
	#set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, breaks=c(-.01,15,seq(28,48,5),100), labels=c('0-15','16-27','28-32','33-37','38-42','43-47','48-80'))])
	set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,seq(15,65,10),100), labels=c('0-14','15-24','25-34','35-44','45-54','55-64','65-99'))])
	set(dfs, NULL, 'YOUNG', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,15,30,100), labels=c('0-14', '15-29','30-99'))])
	set(dfs, NULL, 'OLD', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,15,50,70,100), labels=c('0-14', '15-49', '50-69', '70-99'))])	
	#
	#	focus on young vs others among MSM West and South (across all regions)
	#
	dc	<- subset(dfs, Trm=='MSM' & YOUNG!='0-14' & Region_first%in%c('Den_Haag','Rotterdam_Rijnmond','South','West','Amsterdam') & AnyPos_T1>2010)
	setnames(dc, 'YOUNG','GROUP')
	tmp	<- newdiag.characteristics.table(dc)
	write.csv(tmp, row.names=FALSE, file=paste(outfile,'table_newdiagMSMWest_by_young.csv'))
	dc	<- subset(dfs, Trm=='MSM' & YOUNG!='0-14' & AnyPos_T1>2010)
	setnames(dc, 'YOUNG','GROUP')
	tmp	<- newdiag.characteristics.table(dc)
	write.csv(tmp, row.names=FALSE, file=paste(outfile,'table_newdiagMSM_by_young.csv'))
	dc	<- subset(dfs, Trm=='MSM' & YOUNG!='0-14' & GGD_first%in%c('Den_Haag','Rotterdam_Rijnmond','Amsterdam','Utrecht') & AnyPos_T1>2010)
	setnames(dc, 'YOUNG','GROUP')
	tmp	<- newdiag.characteristics.table(dc)
	write.csv(tmp, row.names=FALSE, file=paste(outfile,'table_newdiagMSM4Cities_by_young.csv'))
	#
	#	focus on old vs others among MSM 
	#
	dc	<- subset(dfs, Trm=='MSM' & !OLD%in%c('0-14','70-99') & AnyPos_T1>2010)
	setnames(dc, 'OLD','GROUP')
	tmp	<- newdiag.characteristics.table(dc)
	write.csv(tmp, row.names=FALSE, file=paste(outfile,'table_newdiagMSM_by_old.csv'))
	dc	<- subset(dfs, Trm=='MSM' & !OLD%in%c('0-14','70-99') & AnyPos_T1>2010 & GGD_first%in%c('Den_Haag','Rotterdam_Rijnmond','Amsterdam','Utrecht') )
	setnames(dc, 'OLD','GROUP')
	tmp	<- newdiag.characteristics.table(dc)
	write.csv(tmp, row.names=FALSE, file=paste(outfile,'table_newdiagMSM4Cities_by_old.csv'))
	
}
######################################################################################
eval.spatial.crudetimetrends.allForeign<- function()
{
	require(maptools)
	require(maps)
	require(spdep)
	require(spacetime)
	library(geoR)
	require(xts)
	require(INLA)
	require(data.table)
	library(RColorBrewer)
	infile			<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/ATHENA_1610_INLA_preprocessed.rda'	
	outfile.base	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	#
	#	load: 
	#		dd 		(preprocessed HIV diagnoses + population counts for 2010-2015 in Netherlands)
	#		ggd.shp (GGD shape file)
	#	
	load(infile)
	#
	#	prepare shape file  
	#
	ggd.shp@data$GGD_INLA_IDX	<- seq_len(nrow(ggd.shp@data))		
	#
	#	plot diagnoses since 2008 for every GGD among all foreign-born individuals
	#	
	da		<- subset(dm, STAT=='ORIGIN_TRM' & grepl('Foreign',GROUP) & YEAR>=2004)
	da		<- da[, list(Foreign_DIAG=sum(NEW_DIAG)), by=c('GGD','REGION','GGD_ID','GGD_INLA_IDX','YEAR')]
	tmp		<- subset(da, YEAR>=2004 & YEAR<=2007)[, list(YEAR=2007, Foreign_DIAG=mean(Foreign_DIAG)), by=c('GGD','REGION','GGD_ID','GGD_INLA_IDX')]
	da		<- rbind(tmp, subset(da, YEAR>2007))	
	tmp		<- unique(subset(da, select=c(GGD, REGION, GGD_ID, GGD_INLA_IDX)))
	tmp2	<- as.data.table(expand.grid(GGD=da[, unique(GGD)], YEAR=da[, unique(YEAR)]))
	tmp		<- merge(tmp, tmp2, by=c('GGD'))	
	da		<- merge(tmp, da, by=c('GGD','REGION','GGD_ID','GGD_INLA_IDX','YEAR'), all.x=1)
	set(da, da[, which(is.na(Foreign_DIAG))], 'Foreign_DIAG', 0)			
	setkey(da, GGD_INLA_IDX, YEAR)
	da[, YEAR_Dt:=as.Date(paste0(YEAR,'-12-31'),format="%Y-%m-%d")]	
	cols	<- colorRampPalette(brewer.pal(9, "RdYlBu"))
	#	plot by time
	tmp		<- dcast.data.table(da, GGD_INLA_IDX~YEAR_Dt, value.var='Foreign_DIAG')	
	tmp		<- STFDF(	as(ggd.shp,"SpatialPolygons"),
			xts(1:ncol(tmp[,-1]), as.Date(colnames(tmp[,-1]),format="%Y-%m-%d")), # xts(num_time_points,date)
			data.frame(Foreign_MSM=unlist(tmp[,-1]),row.names=1:((dim(tmp[,-1])[1])*(dim(tmp[,-1])[2])))) # data.frame(vector values district order)
	pdf(file=paste(outfile.base, 'map_diagnumber_GGD_ForeignAll.pdf'), w=9, h=6)	
	stplot(tmp, 			da[, unique(YEAR_Dt)], 
			at= c(-1,5,10,25,50,100,200),
			col.regions=rev(cols(8)),		
			main='annual HIV diagnoses among foreign-born individuals',
			animate=0 )
	dev.off()
	#
	#	average p.a. 2013-2015
	#	
	da		<- subset(dm, STAT=='ORIGIN_TRM' & grepl('Foreign',GROUP) & YEAR>=2013)
	da		<- da[, list(Foreign_DIAG=sum(NEW_DIAG)), by=c('GGD','REGION','GGD_ID','GGD_INLA_IDX','YEAR')]
	tmp		<- da[, list(Foreign_DIAG=mean(Foreign_DIAG)), by=c('GGD','REGION','GGD_ID','GGD_INLA_IDX')]
	#1:    Amsterdam          Amsterdam   3406           15    84.000000
	#3:    Den_Haag           Den_Haag   3906           13    24.333333
	#20:   Rotterdam_Rijnmond Rotterdam_Rijnmond   4607           10    44.000000

	#	plot 
	tmp		<- merge(ggd.shp@data, subset(tmp, select=c(GGD_INLA_IDX, Foreign_DIAG)), by='GGD_INLA_IDX')
	ggd.out	<- copy(ggd.shp)
	attr(ggd.out, "data")	<- tmp
	pdf(file=paste(outfile.base, 'map_diagnumber_GGD_ForeignAll_20132015.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='Foreign_DIAG', asp=1, col.regions=rev(cols(8)), at=c(-1,5,10,20,40,80,100))
	dev.off()	
	pdf(file=paste(outfile.base, 'map_INLA_logzeta_contrastage_youngMSM.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='LOG_ZETA_CONTRAST', asp=1)
	dev.off()	
	pdf(file=paste(outfile.base, 'map_INLA_logzetap0_contrastage_youngMSM.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='LOG_ZETA_CONTRAST_P0', asp=1)
	dev.off()	
}
######################################################################################
eval.spatial.crudetimetrends.allMSM<- function()
{
	require(maptools)
	require(maps)
	require(spdep)
	require(spacetime)
	library(geoR)
	require(xts)
	require(INLA)
	require(data.table)
	library(RColorBrewer)
	infile			<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/ATHENA_1610_INLA_preprocessed.rda'	
	outfile.base	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	#
	#	load: 
	#		dd 		(preprocessed HIV diagnoses + population counts for 2010-2015 in Netherlands)
	#		ggd.shp (GGD shape file)
	#	
	load(infile)
	#
	#	prepare shape file  
	#
	ggd.shp@data$GGD_INLA_IDX	<- seq_len(nrow(ggd.shp@data))		
	#
	#	plot change in overall diagnosis rates since 2008 for every GGD
	#	
	da		<- subset(dd, STAT=='AGE' & GROUP!='0-14' & YEAR>=2004)[, list(M=sum(M), F=sum(F), POP=sum(POP), DIAG_MSM=sum(DIAG_MSM)), by=c('GGD','GGD_ID','GGD_INLA_IDX','REGION','YEAR')]
	set(da, NULL, 'DIAG_MSM', da[, as.double(DIAG_MSM)])
	da[, R_MSM:= DIAG_MSM/M]
	tmp		<- subset(da, YEAR>=2004 & YEAR<=2007)[, list(BASELINE_DIAG_MSM=mean(DIAG_MSM), BASELINE_R_MSM=mean(DIAG_MSM)/mean(M)), by='GGD']
	da		<- merge(da, tmp, by='GGD')
	da		<- subset(da, YEAR>=2007)
	tmp		<- da[, which(YEAR==2007)]
	set(da, tmp, 'DIAG_MSM', da[tmp,BASELINE_DIAG_MSM])
	set(da, tmp, 'R_MSM', da[tmp,BASELINE_R_MSM])
	da[, RCH_DIAG_MSM:= DIAG_MSM/BASELINE_DIAG_MSM]
	da[, RCH_R_MSM:= R_MSM/BASELINE_R_MSM]
	setkey(da, GGD_INLA_IDX, YEAR)
	da[, YEAR_Dt:=as.Date(paste0(YEAR,'-12-31'),format="%Y-%m-%d")]	
	cols	<- colorRampPalette(brewer.pal(9, "RdYlBu"))
	#	plot
	tmp		<- dcast.data.table(da, GGD_INLA_IDX~YEAR_Dt, value.var='RCH_DIAG_MSM')	
	tmp		<- STFDF(	as(ggd.shp,"SpatialPolygons"),
			xts(1:ncol(tmp[,-1]), as.Date(colnames(tmp[,-1]),format="%Y-%m-%d")), # xts(num_time_points,date)
			data.frame(RCH_DIAG_MSM=unlist(tmp[,-1]),row.names=1:((dim(tmp[,-1])[1])*(dim(tmp[,-1])[2])))) # data.frame(vector values district order)
	pdf(file=paste(outfile.base, 'map_chdiagnumber_GGD_MSM.pdf'), w=9, h=6)	
	stplot(tmp, 			da[, unique(YEAR_Dt)], 
			at= c(-0.001,1/seq(2,1,-.2), seq(1,2,.2),15),
			col.regions=rev(cols(24)),		
			main='change in new HIV diagnoses among MSM\nrelative to 2004-2007',
			animate=0 )
	dev.off()
	tmp		<- dcast.data.table(da, GGD_INLA_IDX~YEAR_Dt, value.var='RCH_R_MSM')	
	tmp		<- STFDF(	as(ggd.shp,"SpatialPolygons"),
			xts(1:ncol(tmp[,-1]), as.Date(colnames(tmp[,-1]),format="%Y-%m-%d")), # xts(num_time_points,date)
			data.frame(RCH_DIAG_MSM=unlist(tmp[,-1]),row.names=1:((dim(tmp[,-1])[1])*(dim(tmp[,-1])[2])))) # data.frame(vector values district order)
	pdf(file=paste(outfile.base, 'map_chdiagrates_GGD_MSM.pdf'), w=9, h=6)	
	stplot(tmp, 			da[, unique(YEAR_Dt)], 
			at= c(-0.001,1/seq(2,1,-.2), seq(1,2,.2),15),
			col.regions=rev(cols(24)),		
			main='change in new HIV diagnosis rates among MSM\nrelative to 2004-2007',
			animate=0 )
	dev.off()
}
######################################################################################
eval.spatial.crudetimetrends.youngMSM<- function()
{
	require(maptools)
	require(maps)
	require(spdep)
	require(spacetime)
	library(geoR)
	require(xts)
	require(INLA)
	require(data.table)
	library(RColorBrewer)
	infile			<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/ATHENA_1610_INLA_preprocessed.rda'	
	outfile.base	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	#
	#	load: 
	#		dd 		(preprocessed HIV diagnoses + population counts for 2010-2015 in Netherlands)
	#		ggd.shp (GGD shape file)
	#	
	load(infile)
	#
	#	prepare shape file  
	#
	ggd.shp@data$GGD_INLA_IDX	<- seq_len(nrow(ggd.shp@data))		
	#
	#	plot change in overall diagnosis rates since 2008 for every GGD
	#	
	dy		<- subset(dd, STAT=='YOUNG' & GROUP=='15-29' & YEAR>=2004)[, list(M=sum(M), F=sum(F), POP=sum(POP), DIAG_MSM=sum(DIAG_MSM)), by=c('GGD','GGD_ID','GGD_INLA_IDX','REGION','YEAR')]
	set(dy, NULL, 'DIAG_MSM', dy[, as.double(DIAG_MSM)])
	dy[, R_MSM:= DIAG_MSM/M]
	tmp		<- subset(dy, YEAR>=2004 & YEAR<=2007)[, list(BASELINE_DIAG_MSM=mean(DIAG_MSM), BASELINE_R_MSM=mean(DIAG_MSM)/mean(M)), by='GGD']
	dy		<- merge(dy, tmp, by='GGD')
	dy		<- subset(dy, YEAR>=2007)
	tmp		<- dy[, which(YEAR==2007)]
	set(dy, tmp, 'DIAG_MSM', dy[tmp,BASELINE_DIAG_MSM])
	set(dy, tmp, 'R_MSM', dy[tmp,BASELINE_R_MSM])
	dy[, RCH_DIAG_MSM:= DIAG_MSM/BASELINE_DIAG_MSM]
	dy[, RCH_R_MSM:= R_MSM/BASELINE_R_MSM]
	setkey(dy, GGD_INLA_IDX, YEAR)
	dy[, YEAR_Dt:=as.Date(paste0(YEAR,'-12-31'),format="%Y-%m-%d")]	
	cols	<- colorRampPalette(brewer.pal(9, "RdYlBu"))
	#	plot
	tmp		<- dcast.data.table(dy, GGD_INLA_IDX~YEAR_Dt, value.var='RCH_DIAG_MSM')	
	tmp		<- STFDF(	as(ggd.shp,"SpatialPolygons"),
			xts(1:ncol(tmp[,-1]), as.Date(colnames(tmp[,-1]),format="%Y-%m-%d")), # xts(num_time_points,dyte)
			data.frame(RCH_DIAG_MSM=unlist(tmp[,-1]),row.names=1:((dim(tmp[,-1])[1])*(dim(tmp[,-1])[2])))) # data.frame(vector values district order)
	pdf(file=paste(outfile.base, 'map_chdiagnumber_GGD_youngMSM.pdf'), w=9, h=6)	
	stplot(tmp, 			dy[, unique(YEAR_Dt)], 
			at= c(-0.001,1/seq(2,1,-.2), seq(1,2,.2),15),
			col.regions=rev(cols(24)),		
			main='change in new HIV diagnoses among MSM 15-29 yrs\nrelative to 2004-2007',
			animate=0 )
	dev.off()
	tmp		<- dcast.data.table(dy, GGD_INLA_IDX~YEAR_Dt, value.var='RCH_R_MSM')	
	tmp		<- STFDF(	as(ggd.shp,"SpatialPolygons"),
			xts(1:ncol(tmp[,-1]), as.Date(colnames(tmp[,-1]),format="%Y-%m-%d")), # xts(num_time_points,date)
			data.frame(RCH_DIAG_MSM=unlist(tmp[,-1]),row.names=1:((dim(tmp[,-1])[1])*(dim(tmp[,-1])[2])))) # data.frame(vector values district order)
	pdf(file=paste(outfile.base, 'map_chdiagrates_GGD_youngMSM.pdf'), w=9, h=6)	
	stplot(tmp, 			dy[, unique(YEAR_Dt)], 
			at= c(-0.001,1/seq(2,1,-.2), seq(1,2,.2),15),
			col.regions=rev(cols(24)),		
			main='change in new HIV diagnosis rates among MSM 15-29 yrs\nrelative to 2004-2007',
			animate=0 )
	dev.off()
}
######################################################################################
eval.spatial.crudetimetrends.midMSM<- function()
{
	require(maptools)
	require(maps)
	require(spdep)
	require(spacetime)
	library(geoR)
	require(xts)
	require(INLA)
	require(data.table)
	library(RColorBrewer)
	infile			<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/ATHENA_1610_INLA_preprocessed.rda'	
	outfile.base	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	#
	#	load: 
	#		dd 		(preprocessed HIV diagnoses + population counts for 2010-2015 in Netherlands)
	#		ggd.shp (GGD shape file)
	#	
	load(infile)
	#
	#	prepare shape file  
	#
	ggd.shp@data$GGD_INLA_IDX	<- seq_len(nrow(ggd.shp@data))		
	#
	#	plot change in overall diagnosis rates since 2008 for every GGD
	#	
	dmid	<- subset(dd, STAT=='MID' & GROUP=='30-49' & YEAR>=2004)[, list(M=sum(M), F=sum(F), POP=sum(POP), DIAG_MSM=sum(DIAG_MSM)), by=c('GGD','GGD_ID','GGD_INLA_IDX','REGION','YEAR')]
	set(dmid, NULL, 'DIAG_MSM', dmid[, as.double(DIAG_MSM)])
	dmid[, R_MSM:= DIAG_MSM/M]
	tmp		<- subset(dmid, YEAR>=2004 & YEAR<=2007)[, list(BASELINE_DIAG_MSM=mean(DIAG_MSM), BASELINE_R_MSM=mean(DIAG_MSM)/mean(M)), by='GGD']
	dmid	<- merge(dmid, tmp, by='GGD')
	dmid	<- subset(dmid, YEAR>=2007)
	tmp		<- dmid[, which(YEAR==2007)]
	set(dmid, tmp, 'DIAG_MSM', dmid[tmp,BASELINE_DIAG_MSM])
	set(dmid, tmp, 'R_MSM', dmid[tmp,BASELINE_R_MSM])
	dmid[, RCH_DIAG_MSM:= DIAG_MSM/BASELINE_DIAG_MSM]
	dmid[, RCH_R_MSM:= R_MSM/BASELINE_R_MSM]
	setkey(dmid, GGD_INLA_IDX, YEAR)
	dmid[, YEAR_Dt:=as.Date(paste0(YEAR,'-12-31'),format="%Y-%m-%d")]	
	cols	<- colorRampPalette(brewer.pal(9, "RdYlBu"))
	#	plot
	tmp		<- dcast.data.table(dmid, GGD_INLA_IDX~YEAR_Dt, value.var='RCH_DIAG_MSM')	
	tmp		<- STFDF(	as(ggd.shp,"SpatialPolygons"),
			xts(1:ncol(tmp[,-1]), as.Date(colnames(tmp[,-1]),format="%Y-%m-%d")), # xts(num_time_points,dyte)
			data.frame(RCH_DIAG_MSM=unlist(tmp[,-1]),row.names=1:((dim(tmp[,-1])[1])*(dim(tmp[,-1])[2])))) # data.frame(vector values district order)
	pdf(file=paste(outfile.base, 'map_chdiagnumber_GGD_midMSM.pdf'), w=9, h=6)	
	stplot(tmp, 			dmid[, unique(YEAR_Dt)], 
			at= c(-0.001,1/seq(2,1,-.2), seq(1,2,.2),20),
			col.regions=rev(cols(24)),		
			main='change in new HIV diagnoses among MSM 30-49 yrs\nrelative to 2004-2007',
			animate=0 )
	dev.off()
	tmp		<- dcast.data.table(dmid, GGD_INLA_IDX~YEAR_Dt, value.var='RCH_R_MSM')	
	tmp		<- STFDF(	as(ggd.shp,"SpatialPolygons"),
			xts(1:ncol(tmp[,-1]), as.Date(colnames(tmp[,-1]),format="%Y-%m-%d")), # xts(num_time_points,date)
			data.frame(RCH_DIAG_MSM=unlist(tmp[,-1]),row.names=1:((dim(tmp[,-1])[1])*(dim(tmp[,-1])[2])))) # data.frame(vector values district order)
	pdf(file=paste(outfile.base, 'map_chdiagrates_GGD_midMSM.pdf'), w=9, h=6)	
	stplot(tmp, 			dmid[, unique(YEAR_Dt)], 
			at= c(-0.001,1/seq(2,1,-.2), seq(1,2,.2),20),
			col.regions=rev(cols(24)),		
			main='change in new HIV diagnosis rates among MSM 30-49 yrs\nrelative to 2004-2007',
			animate=0 )
	dev.off()
}
######################################################################################
eval.spatial.crudetimetrends.oldMSM<- function()
{
	require(maptools)
	require(maps)
	require(spdep)
	require(spacetime)
	library(geoR)
	require(xts)
	require(INLA)
	require(data.table)
	library(RColorBrewer)
	infile			<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/ATHENA_1610_INLA_preprocessed.rda'	
	outfile.base	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	#
	#	load: 
	#		dd 		(preprocessed HIV diagnoses + population counts for 2010-2015 in Netherlands)
	#		ggd.shp (GGD shape file)
	#	
	load(infile)
	#
	#	prepare shape file  
	#
	ggd.shp@data$GGD_INLA_IDX	<- seq_len(nrow(ggd.shp@data))		
	#
	#	plot change in overall diagnosis rates since 2008 for every GGD
	#	
	do		<- subset(dd, STAT=='OLD' & GROUP=='50-69' & YEAR>=2004)[, list(M=sum(M), F=sum(F), POP=sum(POP), DIAG_MSM=sum(DIAG_MSM)), by=c('GGD','GGD_ID','GGD_INLA_IDX','REGION','YEAR')]
	set(do, NULL, 'DIAG_MSM', do[, as.double(DIAG_MSM)])
	do[, R_MSM:= DIAG_MSM/M]
	tmp		<- subset(do, YEAR>=2004 & YEAR<=2007)[, list(BASELINE_DIAG_MSM=mean(DIAG_MSM), BASELINE_R_MSM=mean(DIAG_MSM)/mean(M)), by='GGD']
	do		<- merge(do, tmp, by='GGD')
	do		<- subset(do, YEAR>=2007)
	tmp		<- do[, which(YEAR==2007)]
	set(do, tmp, 'DIAG_MSM', do[tmp,BASELINE_DIAG_MSM])
	set(do, tmp, 'R_MSM', do[tmp,BASELINE_R_MSM])
	do[, RCH_DIAG_MSM:= DIAG_MSM/BASELINE_DIAG_MSM]
	do[, RCH_R_MSM:= R_MSM/BASELINE_R_MSM]
	setkey(do, GGD_INLA_IDX, YEAR)
	do[, YEAR_Dt:=as.Date(paste0(YEAR,'-12-31'),format="%Y-%m-%d")]	
	cols	<- colorRampPalette(brewer.pal(9, "RdYlBu"))
	#	plot
	tmp		<- dcast.data.table(do, GGD_INLA_IDX~YEAR_Dt, value.var='RCH_DIAG_MSM')	
	tmp		<- STFDF(	as(ggd.shp,"SpatialPolygons"),
			xts(1:ncol(tmp[,-1]), as.Date(colnames(tmp[,-1]),format="%Y-%m-%d")), # xts(num_time_points,dyte)
			data.frame(RCH_DIAG_MSM=unlist(tmp[,-1]),row.names=1:((dim(tmp[,-1])[1])*(dim(tmp[,-1])[2])))) # data.frame(vector values district order)
	pdf(file=paste(outfile.base, 'map_chdiagnumber_GGD_oldMSM.pdf'), w=9, h=6)	
	stplot(tmp, 			do[, unique(YEAR_Dt)], 
			at= c(-0.001,1/seq(2,1,-.2), seq(1,2,.2),20),
			col.regions=rev(cols(24)),		
			main='change in new HIV diagnoses among MSM 50-69 yrs\nrelative to 2004-2007',
			animate=0 )
	dev.off()
	tmp		<- dcast.data.table(do, GGD_INLA_IDX~YEAR_Dt, value.var='RCH_R_MSM')	
	tmp		<- STFDF(	as(ggd.shp,"SpatialPolygons"),
			xts(1:ncol(tmp[,-1]), as.Date(colnames(tmp[,-1]),format="%Y-%m-%d")), # xts(num_time_points,date)
			data.frame(RCH_DIAG_MSM=unlist(tmp[,-1]),row.names=1:((dim(tmp[,-1])[1])*(dim(tmp[,-1])[2])))) # data.frame(vector values district order)
	pdf(file=paste(outfile.base, 'map_chdiagrates_GGD_oldMSM.pdf'), w=9, h=6)	
	stplot(tmp, 			do[, unique(YEAR_Dt)], 
			at= c(-0.001,1/seq(2,1,-.2), seq(1,2,.2),20),
			col.regions=rev(cols(24)),		
			main='change in new HIV diagnosis rates among MSM 50-69 yrs\nrelative to 2004-2007',
			animate=0 )
	dev.off()
}
######################################################################################
eval.spatial.INLA.allMSM.2010to2015<- function()
{
	require(maptools)
	require(maps)
	require(spdep)
	require(spacetime)
	library(geoR)
	require(xts)
	require(INLA)
	require(data.table)
	library(RColorBrewer)
	infile			<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/ATHENA_1610_INLA_preprocessed.rda'
	ggd.inla.file	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/GGD_2012_INLA.graph"
	outfile.base	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	#
	#	load: 
	#		dd 		(preprocessed HIV diagnoses + population counts for 2010-2015 in Netherlands)
	#		ggd.shp (GGD shape file)
	#	
	load(infile)
	#
	#	prepare shape file for INLA 
	#
	ggd.shp@data$GGD_INLA_IDX	<- seq_len(nrow(ggd.shp@data))	
	nb2INLA(ggd.inla.file, poly2nb(ggd.shp))
	ggd.inla		<- inla.read.graph(filename=ggd.inla.file)
	image(inla.graph2matrix(ggd.inla), xlab='', ylab='')
	
	#	plot GGD names
	pdf(file=paste(outfile.base, 'map_GGD_names.pdf'), w=6, h=6)
	par(mar=c(0,0,0,0))
	plot(ggd.shp)
	text(getSpPPolygonsLabptSlots(ggd.shp), labels=ggd.shp$naam, cex=0.6)
	dev.off()
	#	plot population 2012
	tmp		<- colorRampPalette(brewer.pal(9, "Blues")[-1])
	pdf(file=paste(outfile.base, 'map_GGD_pop2012.pdf'), w=6, h=6)
	par(mar=c(0,0,0,0))
	spplot(ggd.shp,"inwonertal",cuts=10, col.regions=rev(tmp(11)))
	plot(ggd.shp)	
	dev.off()
	
	#	GGD with 1 neighbour
	plot(ggd.shp)
	tmp	<- which(1==sapply( poly2nb(ggd.shp), length ))
	for(i in tmp)
		plot(ggd.shp[i,], col='blue', add=T)	
	#	TODO: add index to poly2nb(ggd.shp) if you want to force other neighbourhoods
	
	#
	#
	#	AIM1	excess diagnosis rates among all MSM by GGD 2010-2015
	#	AIM2	changes in diagnosis rates among all MSM by GGD 2000-2004, 2005-2009, 2010-2015
	#	AIM3	contrast in spatially structured rndm effect in diagnosis rates among MSM 15-29 from age-average
	#
	#
	#	AIM1
	#	
	da		<- subset(dd, STAT=='AGE' & GROUP!='0-14')[, list(M=sum(M), F=sum(F), POP=sum(POP), DIAG_MSM=sum(DIAG_MSM)), by=c('GGD','GGD_ID','GGD_INLA_IDX','REGION','YEAR')]
	da		<- subset(da, YEAR>=2010)
	#	calculate E
	da[, E_NL:= sum(DIAG_MSM)/length(unique(YEAR))]	
	tmp		<- da[, list(M=mean(M)), by=c('GGD')]
	tmp[, Mp:= M/sum(M)]		
	da		<- merge(da, subset(tmp, select=c(GGD, Mp)), by='GGD')
	set(da, NULL, 'E_GGD', da[, E_NL*Mp])	
	set(da, NULL, 'Mp', NULL)
	#	add crude diag rate in period 2010-2015
	tmp2	<- tmp[, sum(M)]
	tmp		<- da[, list(R_MSM=mean(DIAG_MSM)/mean(M)), by=c('GGD')]
	da		<- merge(da, tmp, by='GGD')
	da[, RR_MSM:= R_MSM/(E_NL/tmp2)]
	#
	da.f	<- DIAG_MSM ~ 1 + f(GGD_INLA_IDX, model='bym', 
									graph=ggd.inla.file,
									scale.model=TRUE)
	da.m	<- inla(da.f, family='poisson', data=da, E=da$E_GGD, control.compute(dic=TRUE), control.predictor=list(compute=TRUE))
	#
	#	posterior mean for sp str rnd effect
	#
	ggd.n	<- da[, length(unique(GGD))]
	tmp		<- da.m$marginals.random[[1]][1:ggd.n]
	zeta	<- lapply(tmp, function(x) inla.emarginal(exp,x))
	tmp		<- da.m$marginals.random[[1]][(ggd.n+1):(2*ggd.n)]
	us		<- lapply(tmp, function(x) inla.emarginal(exp,x))
	da.p	<- merge(unique(da, by='GGD'),data.table(GGD_INLA_IDX=1:ggd.n, ZETA=as.numeric(unlist(zeta))), by='GGD_INLA_IDX')
	da.p	<- merge(da.p,data.table(GGD_INLA_IDX=1:ggd.n, U=as.numeric(unlist(us))), by='GGD_INLA_IDX')
	tmp		<- merge(ggd.shp@data, subset(da.p, select=c(GGD_INLA_IDX, R_MSM, RR_MSM, ZETA)), by='GGD_INLA_IDX')
	ggd.out	<- copy(ggd.shp)
	attr(ggd.out, "data")	<- tmp
	pdf(file=paste(outfile.base, 'map_INLA1.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='ZETA', asp=1)
	dev.off()	
	pdf(file=paste0(outfile.base, 'map_rcrude_diagnoses_MSM_2010_2015.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='RR_MSM', asp=1)
	dev.off()	
	#
	#	variance explained through spatial str rnd effect
	#
	vars	<- data.table(GGD_INLA_IDX=1:ggd.n, DUMMY=(ggd.n+1):(2*ggd.n))[, {
							list(VAR_U= var(inla.rmarginal(1e5, da.m$marginals.random[['GGD_INLA_IDX']][[DUMMY]])))
						}, by=c('GGD_INLA_IDX','DUMMY')]
	vars[, VAR_V:= var(inla.rmarginal(1e5, inla.tmarginal(function(x) 1/x, da.m$marginals.hyper[['Precision for GGD_INLA_IDX (iid component)']])))]	
	vars[, mean(VAR_U/(VAR_U+VAR_V))]
	#	0.300615
}
######################################################################################
eval.spatial.INLA.youngMSM.2010to2015<- function()
{
	require(maptools)
	require(maps)
	require(spdep)
	require(spacetime)
	library(geoR)
	require(xts)
	require(INLA)
	require(data.table)
	library(RColorBrewer)
	infile			<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/ATHENA_1610_INLA_preprocessed.rda'
	ggd.inla.file	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/GGD_2012_INLA.graph"
	outfile.base	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	#
	#	load: 
	#		dd 		(preprocessed HIV diagnoses + population counts for 2010-2015 in Netherlands)
	#		ggd.shp (GGD shape file)
	#	
	load(infile)
	#
	#	prepare shape file for INLA 
	#
	ggd.shp@data$GGD_INLA_IDX	<- seq_len(nrow(ggd.shp@data))	
	nb2INLA(ggd.inla.file, poly2nb(ggd.shp))
	ggd.inla		<- inla.read.graph(filename=ggd.inla.file)	
	#
	#	AIM3	contrast in spatially structured rndm effect in diagnosis rates among MSM 15-29 from age-average
	#
	da		<- subset(dd, STAT=='YOUNG' & GROUP!='0-14' & YEAR>=2010)	
	da[, R_MSM:=NULL]
	#	calculate E: percent pop in age group in GGD in 2010-2015 	
	tmp		<- da[, list(M=mean(M)), by=c('GGD','GROUP')]
	tmp		<- tmp[, list(GGD=GGD, Mp=M/sum(M)), by='GROUP']			
	da		<- merge(da, subset(tmp, select=c(GGD, GROUP, Mp)), by=c('GGD','GROUP'))
	tmp		<- da[, list(E=sum(DIAG_MSM), M_GROUP=sum(M)), by=c('YEAR','GROUP')]
	tmp		<- tmp[, list(E=mean(E), M_GROUP=mean(M_GROUP)), by='GROUP']
	da		<- merge(da, tmp, by=c('GROUP'))	
	#	add crude diag rate in period 2010-2015	
	tmp		<- da[, list(R_MSM=mean(DIAG_MSM)/mean(M)), by=c('GGD','GROUP')]
	da		<- merge(da, tmp, by=c('GGD','GROUP'))
	#	add crude relative diag rate (relate to pop average by group)
	da[, RR_MSM:= R_MSM/(E/M_GROUP)]	
	set(da, NULL, 'E', da[, E*Mp])	
	set(da, NULL, 'Mp', NULL)	
	set(da, NULL, 'GROUP', da[, factor(as.character(GROUP))])
	#
	da[, GGD_INLA_IDX_YOUNG:= GGD_INLA_IDX]
	da[, GGD_INLA_IDX_OTHER:= GGD_INLA_IDX]
	#da.f	<- DIAG_MSM ~ GROUP + 	f(GGD_INLA_IDX*I(GROUP=='15-29'), model='bym', graph=ggd.inla.file,scale.model=TRUE) +
	#								f(GGD_INLA_IDX2*I(GROUP=='30-69'), model='bym', graph=ggd.inla.file,scale.model=TRUE)	
	set(da, da[, which(GROUP=='15-29')],'GGD_INLA_IDX_OTHER', NA_integer_)
	set(da, da[, which(GROUP=='30-69')],'GGD_INLA_IDX_YOUNG', NA_integer_)
	da.f	<- DIAG_MSM ~ GROUP + 	f(GGD_INLA_IDX_YOUNG, model='bym', graph=ggd.inla.file,scale.model=TRUE) +
									f(GGD_INLA_IDX_OTHER, model='bym', graph=ggd.inla.file,scale.model=TRUE)
	
	da.m	<- inla(da.f, family='poisson', data=da, E=da$E, control.compute(dic=TRUE), control.predictor=list(compute=TRUE))
	#	contrast in posterior means for sp str rnd effect
	ggd.n	<- da[, length(unique(GGD))]		
	lzeta.c	<- da.m$summary.random[["GGD_INLA_IDX_YOUNG"]][1:ggd.n,'mean'] - da.m$summary.random[["GGD_INLA_IDX_OTHER"]][1:ggd.n,'mean']	
	zeta.o	<- as.numeric(unlist(lapply(da.m$marginals.random[["GGD_INLA_IDX_OTHER"]][1:ggd.n], function(x) inla.emarginal(exp,x))))
	us.o	<- as.numeric(unlist(lapply(da.m$marginals.random[["GGD_INLA_IDX_OTHER"]][(ggd.n+1):(2*ggd.n)], function(x) inla.emarginal(exp,x))))
	zeta.y	<- as.numeric(unlist(lapply(da.m$marginals.random[["GGD_INLA_IDX_YOUNG"]][1:ggd.n], function(x) inla.emarginal(exp,x))))
	us.y	<- as.numeric(unlist(lapply(da.m$marginals.random[["GGD_INLA_IDX_YOUNG"]][(ggd.n+1):(2*ggd.n)], function(x) inla.emarginal(exp,x))))
	#	quantiles of the contrast spatial effect YOUNG-OTHER 
	#	on log scale, the marginals are Laplace approximated normal distributions
	#	simulate from YOUNG and OTHER marginal posterior
	#	take the difference
	tmp			<- cbind( da.m$summary.random[["GGD_INLA_IDX_YOUNG"]][1:ggd.n,c('mean','sd')], da.m$summary.random[["GGD_INLA_IDX_OTHER"]][1:ggd.n,c('mean','sd')] )
	tmp			<- cbind( tmp[,1]-tmp[,3], sqrt(tmp[,2]*tmp[,2]+tmp[,4]*tmp[,4]) )
	lzeta.c.p0	<- pnorm(0, mean=tmp[,1], sd=tmp[,2], lower.tail=FALSE)
	#	prepare output data.frame	
	da.p		<- melt(da, id.vars=c('GGD','GGD_ID', 'GGD_INLA_IDX', 'REGION', 'STAT', 'GROUP', 'YEAR'))
	da.p		<- subset(da.p, YEAR==2010 & variable%in%c('M','DIAG_MSM','E','R_MSM','RR_MSM'))	
	da.p		<- dcast.data.table(da.p, GGD+GGD_ID+GGD_INLA_IDX+REGION~variable+GROUP, value.var='value')
	tmp			<- data.table(GGD_INLA_IDX=1:ggd.n, ZETA_15_29=zeta.y, ZETA_30_69=zeta.o, U_15_29=us.y, U_30_69=us.o, LOG_ZETA_CONTRAST=lzeta.c, LOG_ZETA_CONTRAST_P0=lzeta.c.p0)	
	da.p		<- merge(da.p, tmp, by='GGD_INLA_IDX')
	da.p[, ZETA_CONTRAST:=ZETA_15_29-ZETA_30_69]
	cols		<- colorRampPalette(brewer.pal(9, "RdYlBu"))
	#	plot 
	tmp		<- merge(ggd.shp@data, subset(da.p, select=c(GGD_INLA_IDX, LOG_ZETA_CONTRAST, LOG_ZETA_CONTRAST_P0, ZETA_CONTRAST)), by='GGD_INLA_IDX')
	ggd.out	<- copy(ggd.shp)
	attr(ggd.out, "data")	<- tmp
	pdf(file=paste(outfile.base, 'map_INLA_zeta_contrastage_youngMSM.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='ZETA_CONTRAST', asp=1)
	dev.off()	
	pdf(file=paste(outfile.base, 'map_INLA_logzeta_contrastage_youngMSM.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='LOG_ZETA_CONTRAST', asp=1)
	dev.off()	
	pdf(file=paste(outfile.base, 'map_INLA_logzetap0_contrastage_youngMSM.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='LOG_ZETA_CONTRAST_P0', asp=1)
	dev.off()	
	#
	#	variance explained through spatial str rnd effect
	#	
	vars	<- as.data.table(expand.grid(GGD_INLA_IDX=1:ggd.n, EFFECT= c('GGD_INLA_IDX_YOUNG','GGD_INLA_IDX_OTHER'), stringsAsFactors=FALSE)) 
	vars	<- merge(vars, data.table(GGD_INLA_IDX=1:ggd.n, DUMMY=(ggd.n+1):(2*ggd.n)), by='GGD_INLA_IDX')	
	vars	<- vars[, list(VAR_U= var(inla.rmarginal(1e5, da.m$marginals.random[[EFFECT]][[DUMMY]]))), by=c('EFFECT','GGD_INLA_IDX','DUMMY')]
	vars[, DUMMY:=NULL]
	tmp		<- unique(vars, by='EFFECT')
	tmp[, DUMMY:= paste0("Precision for ",EFFECT," (iid component)")]
	tmp		<- tmp[, list(VAR_V= var(inla.rmarginal(1e5, inla.tmarginal(function(x) 1/x, da.m$marginals.hyper[[DUMMY]])))), by='EFFECT']
	vars	<- merge(vars, tmp, by='EFFECT')
	vars[, list(VAR_U_PC= mean(VAR_U/(VAR_U+VAR_V))), by='EFFECT']
	#	           EFFECT  VAR_U_PC
	#1: GGD_INLA_IDX_OTHER 0.1716102
	#2: GGD_INLA_IDX_YOUNG 0.1641316

	#
	#	can I calculate the variance explained on the contrast?
	#		1- 	I have mean sd for log.u.contrast. 
	#			I can use this to obtain quantiles (q=0.01 .. 0.99) and pdfs at these quantiles for (x,y).
	#			on this, I can run 'inla.rmarginal' as before
	#		2- 	I can calculate log.v from log.zeta-log.u for both age groups
	#			from this, I can calculate log.v.contrast
	#			from this, I can follow the procedure under (1)
}
######################################################################################
eval.spatial.INLA.oldMSM.2010to2015<- function()
{
	require(maptools)
	require(maps)
	require(spdep)
	require(spacetime)
	library(geoR)
	require(xts)
	require(INLA)
	require(data.table)
	library(RColorBrewer)
	infile			<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/ATHENA_1610_INLA_preprocessed.rda'
	ggd.inla.file	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/GGD_2012_INLA.graph"
	outfile.base	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	#
	#	load: 
	#		dd 		(preprocessed HIV diagnoses + population counts for 2010-2015 in Netherlands)
	#		ggd.shp (GGD shape file)
	#	
	load(infile)
	#
	#	prepare shape file for INLA 
	#
	ggd.shp@data$GGD_INLA_IDX	<- seq_len(nrow(ggd.shp@data))	
	nb2INLA(ggd.inla.file, poly2nb(ggd.shp))
	ggd.inla		<- inla.read.graph(filename=ggd.inla.file)	
	#
	#	AIM3	contrast in spatially structured rndm effect in diagnosis rates among MSM 15-29 from age-average
	#
	da		<- subset(dd, STAT=='OLD' & GROUP!='0-14' & YEAR>=2010)	
	da[, R_MSM:=NULL]
	#	calculate E: percent pop in age group in GGD in 2010-2015 	
	tmp		<- da[, list(M=mean(M)), by=c('GGD','GROUP')]
	tmp		<- tmp[, list(GGD=GGD, Mp=M/sum(M)), by='GROUP']			
	da		<- merge(da, subset(tmp, select=c(GGD, GROUP, Mp)), by=c('GGD','GROUP'))
	tmp		<- da[, list(E=sum(DIAG_MSM), M_GROUP=sum(M)), by=c('YEAR','GROUP')]
	tmp		<- tmp[, list(E=mean(E), M_GROUP=mean(M_GROUP)), by='GROUP']
	da		<- merge(da, tmp, by=c('GROUP'))	
	#	add crude diag rate in period 2010-2015	
	tmp		<- da[, list(R_MSM=mean(DIAG_MSM)/mean(M)), by=c('GGD','GROUP')]
	da		<- merge(da, tmp, by=c('GGD','GROUP'))
	#	add crude relative diag rate (relate to pop average by group)
	da[, RR_MSM:= R_MSM/(E/M_GROUP)]	
	set(da, NULL, 'E', da[, E*Mp])	
	set(da, NULL, 'Mp', NULL)	
	set(da, NULL, 'GROUP', da[, factor(as.character(GROUP))])
	#
	da[, GGD_INLA_IDX_OLD:= GGD_INLA_IDX]
	da[, GGD_INLA_IDX_OTHER:= GGD_INLA_IDX]
	set(da, da[, which(GROUP=='50-69')],'GGD_INLA_IDX_OTHER', NA_integer_)
	set(da, da[, which(GROUP=='15-49')],'GGD_INLA_IDX_OLD', NA_integer_)
	da.f	<- DIAG_MSM ~ GROUP + 	f(GGD_INLA_IDX_OLD, model='bym', graph=ggd.inla.file,scale.model=TRUE) +
									f(GGD_INLA_IDX_OTHER, model='bym', graph=ggd.inla.file,scale.model=TRUE)
	
	da.m	<- inla(da.f, family='poisson', data=da, E=da$E, control.compute(dic=TRUE), control.predictor=list(compute=TRUE))
	#	contrast in posterior means for sp str rnd effect
	ggd.n	<- da[, length(unique(GGD))]		
	lzeta.c	<- da.m$summary.random[["GGD_INLA_IDX_OLD"]][1:ggd.n,'mean'] - da.m$summary.random[["GGD_INLA_IDX_OTHER"]][1:ggd.n,'mean']	
	zeta.oth<- as.numeric(unlist(lapply(da.m$marginals.random[["GGD_INLA_IDX_OTHER"]][1:ggd.n], function(x) inla.emarginal(exp,x))))
	us.oth	<- as.numeric(unlist(lapply(da.m$marginals.random[["GGD_INLA_IDX_OTHER"]][(ggd.n+1):(2*ggd.n)], function(x) inla.emarginal(exp,x))))
	zeta.old<- as.numeric(unlist(lapply(da.m$marginals.random[["GGD_INLA_IDX_OLD"]][1:ggd.n], function(x) inla.emarginal(exp,x))))
	us.old	<- as.numeric(unlist(lapply(da.m$marginals.random[["GGD_INLA_IDX_OLD"]][(ggd.n+1):(2*ggd.n)], function(x) inla.emarginal(exp,x))))
	#	quantiles of the contrast spatial effect YOUNG-OTHER 
	#	on log scale, the marginals are Laplace approximated normal distributions
	#	simulate from YOUNG and OTHER marginal posterior
	#	take the difference
	tmp			<- cbind( da.m$summary.random[["GGD_INLA_IDX_OLD"]][1:ggd.n,c('mean','sd')], da.m$summary.random[["GGD_INLA_IDX_OTHER"]][1:ggd.n,c('mean','sd')] )
	tmp			<- cbind( tmp[,1]-tmp[,3], sqrt(tmp[,2]*tmp[,2]+tmp[,4]*tmp[,4]) )
	lzeta.c.p0	<- pnorm(0, mean=tmp[,1], sd=tmp[,2], lower.tail=FALSE)
	#	prepare output data.frame	
	da.p		<- melt(da, id.vars=c('GGD','GGD_ID', 'GGD_INLA_IDX', 'REGION', 'STAT', 'GROUP', 'YEAR'))
	da.p		<- subset(da.p, YEAR==2010 & variable%in%c('M','DIAG_MSM','E','R_MSM','RR_MSM'))	
	da.p		<- dcast.data.table(da.p, GGD+GGD_ID+GGD_INLA_IDX+REGION~variable+GROUP, value.var='value')
	tmp			<- data.table(GGD_INLA_IDX=1:ggd.n, ZETA_50_69=zeta.old, ZETA_15_49=zeta.oth, U_50_69=us.old, U_15_49=us.oth, LOG_ZETA_CONTRAST=lzeta.c, LOG_ZETA_CONTRAST_P0=lzeta.c.p0)	
	da.p		<- merge(da.p, tmp, by='GGD_INLA_IDX')
	da.p[, ZETA_CONTRAST:=ZETA_50_69-ZETA_15_49]
	cols		<- colorRampPalette(brewer.pal(9, "RdYlBu"))
	#	plot 
	tmp		<- merge(ggd.shp@data, subset(da.p, select=c(GGD_INLA_IDX, LOG_ZETA_CONTRAST, LOG_ZETA_CONTRAST_P0, ZETA_CONTRAST)), by='GGD_INLA_IDX')
	ggd.out	<- copy(ggd.shp)
	attr(ggd.out, "data")	<- tmp
	pdf(file=paste(outfile.base, 'map_INLA_zeta_contrastage_oldMSM.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='ZETA_CONTRAST', asp=1)
	dev.off()	
	pdf(file=paste(outfile.base, 'map_INLA_logzeta_contrastage_oldMSM.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='LOG_ZETA_CONTRAST', asp=1)
	dev.off()	
	pdf(file=paste(outfile.base, 'map_INLA_logzetap0_contrastage_oldMSM.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='LOG_ZETA_CONTRAST_P0', asp=1)
	dev.off()		
}
######################################################################################
eval.spatial.INLA.midMSM.2010to2015<- function()
{
	require(maptools)
	require(maps)
	require(spdep)
	require(spacetime)
	library(geoR)
	require(xts)
	require(INLA)
	require(data.table)
	library(RColorBrewer)
	infile			<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/ATHENA_1610_INLA_preprocessed.rda'
	ggd.inla.file	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/GGD_2012_INLA.graph"
	outfile.base	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	#
	#	load: 
	#		dd 		(preprocessed HIV diagnoses + population counts for 2010-2015 in Netherlands)
	#		ggd.shp (GGD shape file)
	#	
	load(infile)
	#
	#	prepare shape file for INLA 
	#
	ggd.shp@data$GGD_INLA_IDX	<- seq_len(nrow(ggd.shp@data))	
	nb2INLA(ggd.inla.file, poly2nb(ggd.shp))
	ggd.inla		<- inla.read.graph(filename=ggd.inla.file)	
	#
	#	AIM3	contrast in spatially structured rndm effect in diagnosis rates among MSM 15-29 from age-average
	#
	da		<- subset(dd, STAT=='MID' & YEAR>=2010)	
	da[, R_MSM:=NULL]
	#	calculate E: percent pop in age group in GGD in 2010-2015 	
	tmp		<- da[, list(M=mean(M)), by=c('GGD','GROUP')]
	tmp		<- tmp[, list(GGD=GGD, Mp=M/sum(M)), by='GROUP']			
	da		<- merge(da, subset(tmp, select=c(GGD, GROUP, Mp)), by=c('GGD','GROUP'))
	tmp		<- da[, list(E=sum(DIAG_MSM), M_GROUP=sum(M)), by=c('YEAR','GROUP')]
	tmp		<- tmp[, list(E=mean(E), M_GROUP=mean(M_GROUP)), by='GROUP']
	da		<- merge(da, tmp, by=c('GROUP'))	
	#	add crude diag rate in period 2010-2015	
	tmp		<- da[, list(R_MSM=mean(DIAG_MSM)/mean(M)), by=c('GGD','GROUP')]
	da		<- merge(da, tmp, by=c('GGD','GROUP'))
	#	add crude relative diag rate (relate to pop average by group)
	da[, RR_MSM:= R_MSM/(E/M_GROUP)]	
	set(da, NULL, 'E', da[, E*Mp])	
	set(da, NULL, 'Mp', NULL)	
	set(da, NULL, 'GROUP', da[, factor(as.character(GROUP))])
	#
	da[, GGD_INLA_IDX_MID:= GGD_INLA_IDX]
	da[, GGD_INLA_IDX_OTHER:= GGD_INLA_IDX]
	set(da, da[, which(GROUP=='30-49')],'GGD_INLA_IDX_OTHER', NA_integer_)
	set(da, da[, which(GROUP=='15-29,50-69')],'GGD_INLA_IDX_MID', NA_integer_)
	da.f	<- DIAG_MSM ~ GROUP + 	f(GGD_INLA_IDX_MID, model='bym', graph=ggd.inla.file,scale.model=TRUE) +
									f(GGD_INLA_IDX_OTHER, model='bym', graph=ggd.inla.file,scale.model=TRUE)	
	da.m	<- inla(da.f, family='poisson', data=da, E=da$E, control.compute(dic=TRUE), control.predictor=list(compute=TRUE))
	#	contrast in posterior means for sp str rnd effect
	ggd.n	<- da[, length(unique(GGD))]		
	lzeta.c	<- da.m$summary.random[["GGD_INLA_IDX_MID"]][1:ggd.n,'mean'] - da.m$summary.random[["GGD_INLA_IDX_OTHER"]][1:ggd.n,'mean']	
	zeta.oth<- as.numeric(unlist(lapply(da.m$marginals.random[["GGD_INLA_IDX_OTHER"]][1:ggd.n], function(x) inla.emarginal(exp,x))))
	us.oth	<- as.numeric(unlist(lapply(da.m$marginals.random[["GGD_INLA_IDX_OTHER"]][(ggd.n+1):(2*ggd.n)], function(x) inla.emarginal(exp,x))))
	zeta.mid<- as.numeric(unlist(lapply(da.m$marginals.random[["GGD_INLA_IDX_MID"]][1:ggd.n], function(x) inla.emarginal(exp,x))))
	us.mid	<- as.numeric(unlist(lapply(da.m$marginals.random[["GGD_INLA_IDX_MID"]][(ggd.n+1):(2*ggd.n)], function(x) inla.emarginal(exp,x))))
	#	quantiles of the contrast spatial effect YOUNG-OTHER 
	#	on log scale, the marginals are Laplace approximated normal distributions
	#	simulate from YOUNG and OTHER marginal posterior
	#	take the difference
	tmp			<- cbind( da.m$summary.random[["GGD_INLA_IDX_MID"]][1:ggd.n,c('mean','sd')], da.m$summary.random[["GGD_INLA_IDX_OTHER"]][1:ggd.n,c('mean','sd')] )
	tmp			<- cbind( tmp[,1]-tmp[,3], sqrt(tmp[,2]*tmp[,2]+tmp[,4]*tmp[,4]) )
	lzeta.c.p0	<- pnorm(0, mean=tmp[,1], sd=tmp[,2], lower.tail=FALSE)
	#	prepare output data.frame	
	da.p		<- melt(da, id.vars=c('GGD','GGD_ID', 'GGD_INLA_IDX', 'REGION', 'STAT', 'GROUP', 'YEAR'))
	da.p		<- subset(da.p, YEAR==2010 & variable%in%c('M','DIAG_MSM','E','R_MSM','RR_MSM'))	
	da.p		<- dcast.data.table(da.p, GGD+GGD_ID+GGD_INLA_IDX+REGION~variable+GROUP, value.var='value')
	tmp			<- data.table(GGD_INLA_IDX=1:ggd.n, ZETA_30_49=zeta.mid, ZETA_15_29_50_69=zeta.oth, U_30_49=us.mid, U_15_29_50_69=us.oth, LOG_ZETA_CONTRAST=lzeta.c, LOG_ZETA_CONTRAST_P0=lzeta.c.p0)	
	da.p		<- merge(da.p, tmp, by='GGD_INLA_IDX')
	da.p[, ZETA_CONTRAST:=ZETA_30_49-ZETA_15_29_50_69]
	cols		<- colorRampPalette(brewer.pal(9, "RdYlBu"))
	#	plot 
	tmp		<- merge(ggd.shp@data, subset(da.p, select=c(GGD_INLA_IDX, LOG_ZETA_CONTRAST, LOG_ZETA_CONTRAST_P0, ZETA_CONTRAST)), by='GGD_INLA_IDX')
	ggd.out	<- copy(ggd.shp)
	attr(ggd.out, "data")	<- tmp
	pdf(file=paste(outfile.base, 'map_INLA_zeta_contrastage_midMSM.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='ZETA_CONTRAST', asp=1)
	dev.off()	
	pdf(file=paste(outfile.base, 'map_INLA_logzeta_contrastage_midMSM.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='LOG_ZETA_CONTRAST', asp=1)
	dev.off()	
	pdf(file=paste(outfile.base, 'map_INLA_logzetap0_contrastage_midMSM.pdf'), w=6, h=6)
	spplot(obj=ggd.out, zcol='LOG_ZETA_CONTRAST_P0', asp=1)
	dev.off()		
}
######################################################################################
eval.spatial.INLA.allMSM.timetrends<- function()
{
	require(maptools)
	require(maps)
	require(spdep)
	require(spacetime)
	library(geoR)
	require(xts)
	require(INLA)
	require(data.table)
	library(RColorBrewer)
	infile			<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/ATHENA_1610_INLA_preprocessed.rda'
	ggd.inla.file	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/GGD_2012_INLA.graph"
	outfile.base	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	#
	#	load: 
	#		dd 		(preprocessed HIV diagnoses + population counts for 2010-2015 in Netherlands)
	#		ggd.shp (GGD shape file)
	#	
	load(infile)
	#
	#	prepare shape file for INLA 
	#
	ggd.shp@data$GGD_INLA_IDX	<- seq_len(nrow(ggd.shp@data))	
	nb2INLA(ggd.inla.file, poly2nb(ggd.shp))
	ggd.inla		<- inla.read.graph(filename=ggd.inla.file)
	image(inla.graph2matrix(ggd.inla), xlab='', ylab='')		 
	#
	da		<- subset(dd, STAT=='AGE' & GROUP!='0-14' & YEAR>=2004)[, list(M=sum(M), F=sum(F), POP=sum(POP), DIAG_MSM=sum(DIAG_MSM)), by=c('GGD','GGD_ID','GGD_INLA_IDX','REGION','YEAR')]
	set(da, NULL, 'DIAG_MSM', da[, as.double(DIAG_MSM)])
	da[, R_MSM:= DIAG_MSM/M]
	tmp		<- subset(da, YEAR>=2004 & YEAR<=2007)[, list(BASELINE_DIAG_MSM=mean(DIAG_MSM), BASELINE_R_MSM=mean(DIAG_MSM)/mean(M)), by='GGD']
	da		<- merge(da, tmp, by='GGD')
	da		<- subset(da, YEAR>=2007)
	tmp		<- da[, which(YEAR==2007)]
	set(da, tmp, 'DIAG_MSM', da[tmp,BASELINE_DIAG_MSM])
	set(da, tmp, 'R_MSM', da[tmp,BASELINE_R_MSM])
	da[, RCH_DIAG_MSM:= DIAG_MSM/BASELINE_DIAG_MSM]
	da[, RCH_R_MSM:= R_MSM/BASELINE_R_MSM]
	setkey(da, GGD_INLA_IDX, YEAR)
	#
	#	AIM2	changes in diagnosis rates among all MSM by GGD 2000-2004, 2005-2009, 2010-2015
	#
		
	#	calculate E: diagnoses for 2004-2007 across entire Netherlands times proportion per year
	da[, E:=subset(da, YEAR==2007)[, sum(DIAG_MSM)]]
	tmp		<- da[, list(GGD=GGD, Mp=M/sum(M)), by=c('YEAR')]
	da		<- merge(da, tmp, by=c('GGD','YEAR'))
	set(da, NULL, 'E_GGD', da[, E*Mp])	
	set(da, NULL, c('Mp','E'), NULL)
	#
	da.f	<- DIAG_MSM ~ 1 + f(GGD_INLA_IDX, model='bym',																		
			graph=ggd.inla.file,
			group=as.numeric(as.factor(YEAR)),
			control.group=list(model='rw2'),
			scale.model=TRUE)
	da.f	<- DIAG_MSM ~ 1 + f(GGD_INLA_IDX, model='iid',  																		
			graph=ggd.inla.file,
			group=as.numeric(as.factor(YEAR)),
			control.group=list(model='rw1'))
	
	da.f	<- DIAG_MSM ~ 1 + f(GGD_INLA_IDX, model='iid',  																		
									graph=ggd.inla.file)
	da.m	<- inla(da.f, family='poisson', data=da, E=da$E_GGD, verbose=TRUE,																		
									control.compute(dic=TRUE), control.predictor=list(compute=TRUE))
	#	posterior mean for sp str rnd effect
	#	we want to plot log zeta
	ggd.n	<- da[, length(unique(GGD))]	
	lzeta	<- as.data.table(da.m$summary.random[["GGD_INLA_IDX"]])	
	tmp		<- as.data.table(expand.grid(GGD_INLA_IDX=seq_len(ggd.n),YEAR=da[, sort(unique(YEAR))]))	
	lzeta	<- cbind(tmp, subset(lzeta, ID<=ggd.n))
	da.p	<- merge(da, lzeta, by=c('GGD_INLA_IDX','YEAR'))
	setnames(da.p, c('0.025quant','0.5quant','0.975quant'), c('ql','median','qu'))
	ggplot(da.p, aes(x=YEAR)) +
			geom_errorbar(aes(ymin=ql, ymax=qu)) + geom_point(aes(y=mean)) +
			facet_grid(GGD~.) +
			labs(x='', y='log zeta', colour='public health region')
			theme_bw() + theme()
	ggsave(paste0(outfile.base, 'map_logzeta_allMSM_time_2007_2015.pdf'), w=6, h=25)
	
	
}
######################################################################################
eval.spatial.prepare<- function()
{
	require(maptools)
	require(maps)
	require(spdep)
	require(INLA)
	infile.ggd		<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/original_latest/GGD_2012.shp"
	infile.hiv		<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	infile.pop		<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/CBS_1612_Population_GGD_Age.rda'			
	outfile			<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/diagrates_INLA/ATHENA_1610_INLA_preprocessed.rda'
	#
	#	prepare INLA graph file
	#
	ggd.shp			<- readShapeSpatial(infile.ggd)	
	#	ggd.shp@data
	#
	#	load and prepare HIV cases
	#
	load(infile.hiv)
	#	reduce to patients with / without first sequence
	dfs		<- df.all[, {
				z	<- NA_character_
				zz	<- NA_real_
				if(any(!is.na(PosSeqT)))
				{
					z	<- FASTASampleCode[which.min(PosSeqT)]
					zz	<- min(PosSeqT, na.rm=TRUE)
				}
				list(FASTASampleCode=z, PosSeqT=zz)
			}, by=c('Patient')]
	tmp		<- copy(df.all)
	set(tmp, NULL, c('PosSeqT','FASTASampleCode','SUBTYPE_C','SEQ_L'), NULL)
	dfs		<- merge(dfs, unique(tmp), by='Patient')
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')	
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG','SXCH'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(Sex=='M' & Trm=='Het')], 'Trm', 'HetM')
	set(dfs, dfs[, which(Sex=='F' & Trm=='Het')], 'Trm', 'HetF')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unknown')
	set(dfs, dfs[, which(RegionOrigin%in%c("Caribbean","Latin_South_America"))], "RegionOrigin", "Carib_Southern_America")
	set(dfs, dfs[, which(RegionOrigin%in%c("Central_EU","Eastern_EU_stans"))], "RegionOrigin", "Cen_East_Stans_EU")
	set(dfs, dfs[, which(RegionOrigin%in%c("Austr_NZ","Sout_SouthEast_Asia","Oceania_Pacific","North_Africa_Middle_East"))], "RegionOrigin", "Other")
	set(dfs, dfs[, which(is.na(RegionOrigin))], "RegionOrigin", "Unknown")
	set(dfs, NULL, "RegionOrigin", dfs[, as.character(factor(RegionOrigin))])
	set(dfs, NULL, 'ORIGIN_TRM', NA_character_)
	set(dfs, dfs[, which(RegionOrigin=='NL' & Trm=='MSM')], 'ORIGIN_TRM', 'NL_MSM')
	set(dfs, dfs[, which(RegionOrigin=='NL' & Trm!='MSM')], 'ORIGIN_TRM', 'NL_Oth')
	set(dfs, dfs[, which(RegionOrigin!='NL' & Trm=='MSM')], 'ORIGIN_TRM', 'Foreign_MSM')
	set(dfs, dfs[, which(RegionOrigin!='NL' & Trm!='MSM')], 'ORIGIN_TRM', 'Foreign_Oth')
	set(dfs, dfs[, which(RegionOrigin=='Unknown' & Trm=='MSM')], 'ORIGIN_TRM', 'UnknownOri_MSM')
	set(dfs, dfs[, which(RegionOrigin=='Unknown' & Trm!='MSM')], 'ORIGIN_TRM', 'UnknownOri_Oth')	
	set(dfs, NULL, 'Age_AnyPosT1', dfs[, AnyPos_T1-DateBorn])
	set(dfs, NULL, 'YEAR', dfs[, floor(AnyPos_T1)])		
	set(dfs, NULL, 'GGD_first', dfs[, gsub('-','_',gsub('Hulpverlening_','',as.character(GGD_first)))])
	set(dfs, dfs[, which(is.na(GGD_first))], "GGD_first", "Unknown")
	set(dfs, dfs[, which(is.na(Region_first))], 'Region_first', 'Unknown')
	set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,seq(15,65,10),100), labels=c('0-14','15-24','25-34','35-44','45-54','55-64','65-99'))])
	set(dfs, NULL, 'YOUNG', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,15,30, 70, 100), labels=c('0-14', '15-29','30-69', '70-99'))])
	set(dfs, NULL, 'OLD', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,15,50,70,100), labels=c('0-14', '15-49', '50-69', '70-99'))])
	set(dfs, NULL, 'MID', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,15,30,50,70,100), labels=c('0-14', '15-29,50-69', '30-49', '50-69', '70-99'))])
	set(dfs, dfs[, which(MID=='50-69')], 'MID', '15-29,50-69')
	set(dfs, NULL, 'MID', dfs[, factor(as.character(MID))])
	#
	#	focus on 2000-2015 
	#
	dfs		<- subset(dfs, Region_first!='Unknown' & AnyPos_T1>2000 & AnyPos_T1<2016)	
	#stopifnot( !nrow(subset(dfs, is.na(GGD_INLA_IDX))))
	#
	#	get population by AGE and GGD 
	#
	load(infile.pop)	
	set(dpop, NULL, 'GGD', dpop[, gsub('GooienVechtstreek','GooiVechtstreek',as.character(GGD))])
	tmp		<- data.table(GGD=ggd.shp@data$naam, GGD_INLA_IDX= seq_len(nrow(ggd.shp@data)))
	set(tmp, NULL, 'GGD', tmp[, gsub('[^[:alnum:]]','',gsub('^ +','',gsub('â','a',gsub('GGD','',gsub('GG en GD| en','',gsub('Regio|Hulpverlening','',GGD))))))])
	#	merge(unique(subset(dpop, select=GGD)), tmp, all=1, by='GGD')	
	dpop	<- merge(dpop, tmp, all.x=1, by='GGD')
	#
	dpop	<- dpop[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD_INLA_IDX','GGD','AGE','YEAR')]	
	tmp		<- dpop[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD_INLA_IDX','GGD','AGE','YEAR')]
	set(tmp, tmp[, which(AGE%in%c('0-4','5-9','10-14'))], 'AGE', '0-14')
	set(tmp, tmp[, which(AGE%in%c("15-19","20-24"))], 'AGE', '15-24')
	set(tmp, tmp[, which(AGE%in%c("25-29","30-34"))], 'AGE', '25-34')
	set(tmp, tmp[, which(AGE%in%c("35-39","40-44"))], 'AGE', '35-44')
	set(tmp, tmp[, which(AGE%in%c("45-49","50-54"))], 'AGE', '45-54')
	set(tmp, tmp[, which(AGE%in%c("55-59","60-64"))], 'AGE', '55-64')
	set(tmp, tmp[, which(AGE%in%c('65-69',"70-74","75-79","80-84","85-89","90-94","95-99"))], 'AGE', '65-99')
	tmp		<- tmp[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD_INLA_IDX','GGD','AGE','YEAR')]
	tmp[, STAT:='AGE']
	setnames(tmp, 'AGE','GROUP')	
	#	add population by YOUNG and GGD
	tmp2	<- subset(dpop, AGE%in%c("30-34","35-39","40-44","45-49","50-54","55-59","60-64",'65-69'))
	set(tmp2, NULL, 'YOUNG', '30-69')
	tmp3	<- subset(dpop, AGE%in%c("15-19","20-24","25-29"))
	set(tmp3, NULL, 'YOUNG', '15-29')	
	tmp2	<- rbind(tmp2, tmp3)
	tmp2	<- tmp2[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD_INLA_IDX','GGD','YOUNG','YEAR')]
	tmp2[, STAT:='YOUNG']
	setnames(tmp2, 'YOUNG','GROUP')
	tmp		<- rbind(tmp, tmp2)
	#	add population by OLD and GGD
	tmp2	<- subset(dpop, AGE%in%c("15-19","20-24","25-29","30-34","35-39","40-44","45-49"))
	set(tmp2, NULL, 'OLD', '15-49')
	tmp3	<- subset(dpop, AGE%in%c("50-54","55-59","60-69"))
	set(tmp3, NULL, 'OLD', '50-69')	
	tmp2	<- rbind(tmp2, tmp3)
	tmp2	<- tmp2[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD_INLA_IDX','GGD','OLD','YEAR')]
	tmp2[, STAT:='OLD']
	setnames(tmp2, 'OLD','GROUP')
	tmp		<- rbind(tmp, tmp2)
	#	add population by MID and GGD
	tmp2	<- subset(dpop, AGE%in%c("15-19","20-24","25-29","50-54","55-59","60-69"))
	set(tmp2, NULL, 'MID', '"15-29,50-69"')
	tmp3	<- subset(dpop, AGE%in%c("30-34","35-39","40-44","45-49"))
	set(tmp3, NULL, 'MID', '30-49')	
	tmp2	<- rbind(tmp2, tmp3)
	tmp2	<- tmp2[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD_INLA_IDX','GGD','MID','YEAR')]
	tmp2[, STAT:='MID']
	setnames(tmp2, 'MID','GROUP')
	tmp		<- rbind(tmp, tmp2)	
	#	adjust GGD names, add REGION
	dpop	<- copy(tmp)
	dpop[, POP:=M+F]		
	tmp		<- data.table(GGD_= setdiff(dfs[,sort(unique(GGD_first))],'Unknown'), GGD=dpop[,sort(unique(GGD))])
	dpop	<- merge(dpop, tmp, by='GGD')
	set(dpop, NULL, 'GGD', NULL)
	setnames(dpop, 'GGD_', 'GGD')	
	tmp		<- unique(subset(dfs, select=c(Region_first, GGD_first)))
	setnames(tmp, c('Region_first','GGD_first'), c('REGION','GGD'))
	dpop	<- merge(dpop, tmp, by='GGD')
	#
	#	from list of individual-level diagnoses above get total new HIV diagnoses
	#	by age group
	#
	tmp		<- dfs[, list(NEW_DIAG= length(Patient)), by=c('GGD_first','YEAR','AGE','Trm')]
	tmp		<- rbind(tmp, tmp[, list(Trm='ALL', NEW_DIAG=sum(NEW_DIAG)), by=c('GGD_first','YEAR','AGE')])
	tmp[, STAT:='AGE']
	setnames(tmp, 'AGE','GROUP')
	#	... by young
	tmp2	<- dfs[, list(NEW_DIAG= length(Patient)), by=c('GGD_first','YEAR','YOUNG','Trm')]
	tmp2	<- rbind(tmp2, tmp2[, list(Trm='ALL', NEW_DIAG=sum(NEW_DIAG)), by=c('GGD_first','YEAR','YOUNG')])
	tmp2[, STAT:='YOUNG']
	setnames(tmp2, 'YOUNG','GROUP')
	tmp		<- rbind(tmp, tmp2)
	#	... by old
	tmp2	<- dfs[, list(NEW_DIAG= length(Patient)), by=c('GGD_first','YEAR','OLD','Trm')]
	tmp2	<- rbind(tmp2, tmp2[, list(Trm='ALL', NEW_DIAG=sum(NEW_DIAG)), by=c('GGD_first','YEAR','OLD')])
	tmp2[, STAT:='OLD']
	setnames(tmp2, 'OLD','GROUP')
	tmp		<- rbind(tmp, tmp2)
	#	... by mid
	tmp2	<- dfs[, list(NEW_DIAG= length(Patient)), by=c('GGD_first','YEAR','MID','Trm')]
	tmp2	<- rbind(tmp2, tmp2[, list(Trm='ALL', NEW_DIAG=sum(NEW_DIAG)), by=c('GGD_first','YEAR','MID')])
	tmp2[, STAT:='MID']
	setnames(tmp2, 'MID','GROUP')
	dd		<- rbind(tmp, tmp2)
	#
	dd		<- dcast.data.table(dd, GGD_first+YEAR+STAT+GROUP~Trm, value.var='NEW_DIAG')
	#dd		<- subset(dd, YEAR>1999 & YEAR<2016 & GGD_first!='Unknown')
	setnames(dd, c('GGD_first','MSM','IDU','OTH','HetM','HetF','Unknown','ALL'), c('GGD','DIAG_MSM','DIAG_IDU','DIAG_OTH','DIAG_HetM','DIAG_HetF','DIAG_Unknown','DIAG_ALL'))
	dd		<- merge(dpop, dd, by=c('GGD','YEAR','STAT','GROUP'),all.x=1)
	for(x in colnames(dd)[ grepl('^DIAG_',colnames(dd)) ])
		set(dd, which(is.na(dd[[x]])), x, 0L)
	dd[, R_ALL:= DIAG_ALL/POP]
	dd[, R_MSM:= DIAG_MSM/M]
	#
	#	from list of individual-level diagnoses above get total new HIV diagnoses
	#	by region origin
	#
	tmp		<- dfs[, list(NEW_DIAG= length(Patient)), by=c('GGD_first','YEAR','RegionOrigin')]
	tmp[, STAT:='RegionOrigin']
	setnames(tmp, 'RegionOrigin','GROUP')
	#	... by ORIGIN_TRM
	tmp2	<- dfs[, list(NEW_DIAG= length(Patient)), by=c('GGD_first','YEAR','ORIGIN_TRM')]	
	tmp2[, STAT:='ORIGIN_TRM']
	setnames(tmp2, 'ORIGIN_TRM','GROUP')
	dm		<- rbind(tmp, tmp2)	
	setnames(dm, c('GGD_first'), c('GGD'))	
	dm		<- merge(unique(subset(dpop, select=c(GGD, REGION, GGD_ID, GGD_INLA_IDX))), dm, by=c('GGD'), all.x=1)
	
	#	save
	save(dd, dpop, ggd.shp, dm, file=outfile)
}
######################################################################################
eval.diag.rates.by.age<- function()
{
	require(viridis)
	require(scales)
	
	infile.hiv	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	infile.pop	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/CBS_1612_Population_GGD_Age.rda'
	outfile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	load(infile.hiv)
	#
	#	reduce to patients with / without first sequence
	#	
	dfs		<- df.all[, {
				z	<- NA_character_
				zz	<- NA_real_
				if(any(!is.na(PosSeqT)))
				{
					z	<- FASTASampleCode[which.min(PosSeqT)]
					zz	<- min(PosSeqT, na.rm=TRUE)
				}
				list(FASTASampleCode=z, PosSeqT=zz)
			}, by=c('Patient')]
	tmp		<- copy(df.all)
	set(tmp, NULL, c('PosSeqT','FASTASampleCode','SUBTYPE_C','SEQ_L'), NULL)
	dfs		<- merge(dfs, unique(tmp), by='Patient')	
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG','SXCH'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(Sex=='M' & Trm=='Het')], 'Trm', 'HetM')
	set(dfs, dfs[, which(Sex=='F' & Trm=='Het')], 'Trm', 'HetF')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unknown')
	set(dfs, dfs[, which(RegionOrigin%in%c("Central_EU","Eastern_EU_stans"))], "RegionOrigin", "Cen_East_Stans_EU")
	set(dfs, dfs[, which(RegionOrigin%in%c("Austr_NZ","Sout_SouthEast_Asia","Oceania_Pacific","North_Africa_Middle_East"))], "RegionOrigin", "Other")
	set(dfs, dfs[, which(is.na(RegionOrigin))], "RegionOrigin", "Unknown")
	set(dfs, NULL, "RegionOrigin", dfs[, as.character(factor(RegionOrigin))])
	set(dfs, NULL, 'Age_AnyPosT1', dfs[, AnyPos_T1-DateBorn])
	set(dfs, NULL, 'YEAR', dfs[, floor(AnyPos_T1)])		
	set(dfs, dfs[, which(Region_first=='Amst')], 'Region_first', 'North')
	set(dfs, NULL, 'GGD_first', dfs[, gsub('Hulpverlening_','',as.character(GGD_first))])
	set(dfs, dfs[, which(is.na(GGD_first))], "GGD_first", "Unknown")
	set(dfs, NULL, 'GGD_now', dfs[, gsub('Hulpverlening_','',as.character(GGD_now))])
	set(dfs, dfs[, which(is.na(GGD_now))], "GGD_now", "Unknown")
	set(dfs, dfs[, which(is.na(Region_first))], 'Region_first', 'Unknown')
	set(dfs, dfs[, which(Region_first%in%c('Rott','Den_Haag'))], 'Region_first', 'West')
	set(dfs, dfs[, which(GGD_first%in%c("Groningen", "Fryslan", "Drenthe"))], 'Region_first', 'North')
	set(dfs, dfs[, which(GGD_first%in%c("Flevoland","Twente","Gelre_IJssel","Midden_Nederland","Gelderland_Midden","Nijmegen","Rivierenland","IJsselland"))], 'Region_first', 'East')
	set(dfs, dfs[, which(GGD_first%in%c("Amsterdam", 'Den_Haag', "Zuid_Holland_West", "Zeeland", "Gooi_Vechtstreek", "Hollands_Noorden", "Rotterdam_Rijnmond", "Hollands_Midden", "Zaanstreek_Waterland", "Kennemerland", "Utrecht", "Zuid_Holland_Zuid"))], 'Region_first', 'West')
	set(dfs, dfs[, which(GGD_first%in%c("Hart_voor_Brabant","Zuid_Limburg","Brabant_Zuidoost","Limburg-Noord","West_Brabant"))], 'Region_first', 'South')
	#set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, breaks=c(-.01,15,seq(28,48,5),100), labels=c('0-15','16-27','28-32','33-37','38-42','43-47','48-80'))])
	set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,seq(15,65,10),100), labels=c('0-14','15-24','25-34','35-44','45-54','55-64','65-99'))])
	set(dfs, NULL, 'YOUNG', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,15,30,100), labels=c('0-14', '15-29','30-99'))])
	set(dfs, NULL, 'OLD', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,15,50,70,100), labels=c('0-14', '15-49', '50-69', '70-99'))])	
	#
	#	get population by age year ggd 
	#
	load(infile.pop)	
	dpop	<- dpop[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD','AGE','YEAR')]	
	tmp		<- dpop[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD','AGE','YEAR')]
	set(tmp, tmp[, which(AGE%in%c('0-4','5-9','10-14'))], 'AGE', '0-14')
	set(tmp, tmp[, which(AGE%in%c("15-19","20-24"))], 'AGE', '15-24')
	set(tmp, tmp[, which(AGE%in%c("25-29","30-34"))], 'AGE', '25-34')
	set(tmp, tmp[, which(AGE%in%c("35-39","40-44"))], 'AGE', '35-44')
	set(tmp, tmp[, which(AGE%in%c("45-49","50-54"))], 'AGE', '45-54')
	set(tmp, tmp[, which(AGE%in%c("55-59","60-64"))], 'AGE', '55-64')
	set(tmp, tmp[, which(AGE%in%c('65-69',"70-74","75-79","80-84","85-89","90-94","95-99"))], 'AGE', '65-99')
	tmp		<- tmp[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD','AGE','YEAR')]
	tmp[, STAT:='AGE']
	setnames(tmp, 'AGE','GROUP')	
	#
	tmp2	<- subset(dpop, AGE%in%c("15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64",'65-69'))
	set(tmp2, NULL, 'YOUNG', '15-69')
	tmp3	<- subset(dpop, AGE%in%c("15-19","20-24","25-29"))
	set(tmp3, NULL, 'YOUNG', '15-29')	
	tmp2	<- rbind(tmp2, tmp3)
	tmp2	<- tmp2[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD','YOUNG','YEAR')]
	tmp2[, STAT:='YOUNG']
	setnames(tmp2, 'YOUNG','GROUP')
	tmp		<- rbind(tmp, tmp2)
	#
	tmp2	<- subset(dpop, AGE%in%c("15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64",'65-69'))
	set(tmp2, NULL, 'OLD', '15-69')
	tmp3	<- subset(dpop, AGE%in%c("50-54","55-59","60-69"))
	set(tmp3, NULL, 'OLD', '50-69')	
	tmp2	<- rbind(tmp2, tmp3)
	tmp2	<- tmp2[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD','OLD','YEAR')]
	tmp2[, STAT:='OLD']
	setnames(tmp2, 'OLD','GROUP')
	tmp		<- rbind(tmp, tmp2)
	#
	dpop	<- copy(tmp)
	dpop[, POP:=M+F]		
	tmp		<- data.table(GGD_= setdiff(dfs[,sort(unique(GGD_first))],'Unknown'), GGD=dpop[,sort(unique(GGD))])
	dpop	<- merge(dpop, tmp, by='GGD')
	set(dpop, NULL, 'GGD', NULL)
	setnames(dpop, 'GGD_', 'GGD')	
	tmp		<- unique(subset(dfs, select=c(Region_first, GGD_first)))
	setnames(tmp, c('Region_first','GGD_first'), c('REGION','GGD'))
	dpop	<- merge(dpop, tmp, by='GGD')
	#
	#	HIV by age year ggd
	#
	tmp		<- dfs[, list(NEW_DIAG= length(Patient)), by=c('GGD_first','YEAR','AGE','Trm')]
	tmp		<- rbind(tmp, tmp[, list(Trm='ALL', NEW_DIAG=sum(NEW_DIAG)), by=c('GGD_first','YEAR','AGE')])
	tmp[, STAT:='AGE']
	setnames(tmp, 'AGE','GROUP')
	#	young
	tmp2	<- subset(dfs, Age_AnyPosT1>=15 & Age_AnyPosT1<70)
	set(tmp2, NULL, 'YOUNG', '15-69')
	tmp3	<- subset(dfs, YOUNG=='15-29')
	set(tmp3, NULL, 'YOUNG', '15-29')
	tmp2	<- rbind(tmp2, tmp3)	
	tmp2	<- tmp2[, list(NEW_DIAG= length(Patient)), by=c('GGD_first','YEAR','YOUNG','Trm')]
	tmp2	<- rbind(tmp2, tmp2[, list(Trm='ALL', NEW_DIAG=sum(NEW_DIAG)), by=c('GGD_first','YEAR','YOUNG')])
	tmp2[, STAT:='YOUNG']
	setnames(tmp2, 'YOUNG','GROUP')
	tmp		<- rbind(tmp, tmp2)
	#	old
	tmp2	<- subset(dfs, Age_AnyPosT1>=15 & Age_AnyPosT1<70)
	set(tmp2, NULL, 'OLD', '15-69')
	tmp3	<- subset(dfs, OLD=='50-69')
	set(tmp3, NULL, 'OLD', '50-69')
	tmp2	<- rbind(tmp2, tmp3)	
	tmp2	<- tmp2[, list(NEW_DIAG= length(Patient)), by=c('GGD_first','YEAR','OLD','Trm')]
	tmp2	<- rbind(tmp2, tmp2[, list(Trm='ALL', NEW_DIAG=sum(NEW_DIAG)), by=c('GGD_first','YEAR','OLD')])
	tmp2[, STAT:='OLD']
	setnames(tmp2, 'OLD','GROUP')
	tmp		<- rbind(tmp, tmp2)
	#
	tmp		<- dcast.data.table(tmp, GGD_first+YEAR+STAT+GROUP~Trm, value.var='NEW_DIAG')
	tmp		<- subset(tmp, YEAR>1999 & YEAR<2016 & GGD_first!='Unknown')
	setnames(tmp, c('GGD_first','MSM','IDU','OTH','HetM','HetF','Unknown','ALL'), c('GGD','DIAG_MSM','DIAG_IDU','DIAG_OTH','DIAG_HetM','DIAG_HetF','DIAG_Unknown','DIAG_ALL'))
	dpop	<- merge(dpop, tmp, by=c('GGD','YEAR','STAT','GROUP'),all.x=1)
	for(x in colnames(dpop)[ grepl('^DIAG_',colnames(dpop)) ])
		set(dpop, which(is.na(dpop[[x]])), x, 0L)
	dpop[, R_ALL:= DIAG_ALL/POP]
	dpop[, R_MSM:= DIAG_MSM/M]
	#
	#	plot rates by GGD and AGE
	#
	dp		<- subset(dpop, STAT=='AGE')
	set(dp, NULL, 'GROUP', dp[, factor(as.character(GROUP), levels=rev(levels(dpop$GROUP)))])
	tmp		<- subset(dp, GROUP!='0-14')[, list(POP=sum(POP), M=sum(M), DIAG_ALL=sum(DIAG_ALL), DIAG_MSM=sum(DIAG_MSM), DIAG_HET=sum(DIAG_HetM+DIAG_HetF)),by=c('REGION','GGD','YEAR')]
	tmp[, R_ALL:= DIAG_ALL/POP]
	tmp[, R_MSM:= DIAG_MSM/M]	
	tmp[, R_HET:= DIAG_HET/POP]
	tmp2	<- subset(tmp, YEAR==2001)	
	ggplot(tmp, aes(x=YEAR, y=R_ALL)) +
			geom_line(aes(colour=GGD)) + geom_text(data=tmp2, aes(label=GGD), size=2, hjust=0) +
			facet_grid(~REGION) +		
			labs(x='', y='HIV diagnosis rate\n(new diagnoses / population aged 15 or higher)', colour='Public health region') +
			theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(outfile,'diagrate.pdf'), w=10, h=10)
	ggplot(tmp, aes(x=YEAR, y=R_MSM)) +
			geom_line(aes(colour=GGD)) + geom_text(data=tmp2, aes(label=GGD), size=2, hjust=0) +
			facet_grid(~REGION) +		
			labs(x='', y='HIV diagnosis rate\n(new MSM diagnoses / male population aged 15 or higher)', colour='Public health region') +
			theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(outfile,'diagrate_MSM.pdf'), w=10, h=10)
	ggplot(tmp, aes(x=YEAR, y=R_HET)) +
			geom_line(aes(colour=GGD)) + geom_text(data=tmp2, aes(label=GGD), size=2, hjust=0) +
			facet_grid(~REGION) +		
			labs(x='', y='HIV diagnosis rate\n(new HET diagnoses / population aged 15 or higher)', colour='Public health region') +
			theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(outfile,'diagrate_HET.pdf'), w=10, h=10)
	#
	ggplot(subset(dp, GROUP!='0-14'), aes(x=YEAR, y=R_ALL, colour=GROUP)) + geom_line() +
			scale_color_viridis(option="magma", discrete=TRUE, end=0.85) +
			scale_x_continuous(breaks=seq(2000,2020,5), minor_breaks=seq(2000,2020,1)) +
			facet_wrap(~REGION+GGD, ncol=5) +
			labs(x='', y='HIV diagnosis rate\n(new HET diagnoses / population aged 15 or higher)', colour='age group') +
			theme_bw() + theme(legend.position='bottom') 
	ggsave(file=paste(outfile,'diagrate_by_age.pdf'), w=15, h=15)
	ggplot(subset(dp, GROUP!='0-14'), aes(x=YEAR, y=R_MSM, colour=GROUP)) + geom_line() +
			scale_color_viridis(option="magma", discrete=TRUE, end=0.85) +
			scale_x_continuous(breaks=seq(2000,2020,5), minor_breaks=seq(2000,2020,1)) +
			facet_wrap(~REGION+GGD, ncol=5) +
			labs(x='', y='HIV diagnosis rate\n(new diagnoses MSM / male population aged 15 or higher)', colour='age group') +
			theme_bw() + theme(legend.position='bottom') 
	ggsave(file=paste(outfile,'diagrate_by_age_MSM.pdf'), w=15, h=15)
	#
	#	single out prop YOUNG 
	#
	dp		<- subset(dpop, STAT=='YOUNG')
	set(dp, NULL, 'GROUP', dp[, factor(as.character(GROUP), levels=rev(levels(dpop$GROUP)))])
	ggplot(subset(dp, GROUP!='0-14'), aes(x=YEAR, y=R_ALL, colour=GROUP)) + geom_line() +
			scale_color_viridis(option="magma", discrete=TRUE, end=0.85) +
			scale_x_continuous(breaks=seq(2000,2020,5), minor_breaks=seq(2000,2020,1)) +
			facet_wrap(~REGION+GGD, ncol=5) +
			labs(x='', y='HIV diagnosis rate\n(new diagnoses / population aged 15-69)', colour='age group') +
			theme_bw() + theme(legend.position='bottom') 
	ggsave(file=paste(outfile,'diagrate_by_young.pdf'), w=15, h=15)
	ggplot(subset(dp, GROUP!='0-14'), aes(x=YEAR, y=R_MSM, colour=GROUP)) + geom_line() +
			scale_color_viridis(option="magma", discrete=TRUE, end=0.85) +
			scale_x_continuous(breaks=seq(2000,2020,5), minor_breaks=seq(2000,2020,1)) +
			facet_wrap(~REGION+GGD, ncol=5) +
			labs(x='', y='HIV diagnosis rate\n(new diagnoses MSM / male population aged 15-69)', colour='age group') +
			theme_bw() + theme(legend.position='bottom') 
	ggsave(file=paste(outfile,'diagrate_by_young_MSM.pdf'), w=15, h=15)
	#
	#	single out prop OLD 
	#
	dp		<- subset(dpop, STAT=='OLD')
	set(dp, NULL, 'GROUP', dp[, factor(as.character(GROUP), levels=rev(levels(dpop$GROUP)))])
	ggplot(subset(dp, GROUP!='0-14'), aes(x=YEAR, y=R_ALL, colour=GROUP)) + geom_line() +
			scale_color_viridis(option="magma", discrete=TRUE, end=0.85) +
			scale_x_continuous(breaks=seq(2000,2020,5), minor_breaks=seq(2000,2020,1)) +
			facet_wrap(~REGION+GGD, ncol=5) +
			labs(x='', y='HIV diagnosis rate\n(new diagnoses / population aged 15-69)', colour='age group') +
			theme_bw() + theme(legend.position='bottom') 
	ggsave(file=paste(outfile,'diagrate_by_old.pdf'), w=15, h=15)
	ggplot(subset(dp, GROUP!='0-14'), aes(x=YEAR, y=R_MSM, colour=GROUP)) + geom_line() +
			scale_color_viridis(option="magma", discrete=TRUE, end=0.85) +
			scale_x_continuous(breaks=seq(2000,2020,5), minor_breaks=seq(2000,2020,1)) +
			facet_wrap(~REGION+GGD, ncol=5) +
			labs(x='', y='HIV diagnosis rate\n(new diagnoses MSM / male population aged 15-69)', colour='age group') +
			theme_bw() + theme(legend.position='bottom') 
	ggsave(file=paste(outfile,'diagrate_by_old_MSM.pdf'), w=15, h=15)
	
}
######################################################################################
eval.pop.by.age.migration<- function()
{
	require(viridis)
	require(scales)
	
	infile.hiv	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	infile.pop	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/CBS_1612_Population_GGD_Age.rda'
	outfile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	load(infile.hiv)
	#
	#	age group
	#	
	dfs		<- df.all[, {
				z	<- NA_character_
				zz	<- NA_real_
				if(any(!is.na(PosSeqT)))
				{
					z	<- FASTASampleCode[which.min(PosSeqT)]
					zz	<- min(PosSeqT, na.rm=TRUE)
				}
				list(FASTASampleCode=z, PosSeqT=zz)
			}, by=c('Patient')]
	tmp		<- copy(df.all)
	set(tmp, NULL, c('PosSeqT','FASTASampleCode','SUBTYPE_C','SEQ_L'), NULL)
	dfs		<- merge(dfs, unique(tmp), by='Patient')	
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG','SXCH'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(Sex=='M' & Trm=='Het')], 'Trm', 'HetM')
	set(dfs, dfs[, which(Sex=='F' & Trm=='Het')], 'Trm', 'HetF')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unknown')
	set(dfs, dfs[, which(RegionOrigin%in%c("Central_EU","Eastern_EU_stans"))], "RegionOrigin", "Cen_East_Stans_EU")
	set(dfs, dfs[, which(RegionOrigin%in%c("Austr_NZ","Sout_SouthEast_Asia","Oceania_Pacific","North_Africa_Middle_East"))], "RegionOrigin", "Other")
	set(dfs, dfs[, which(is.na(RegionOrigin))], "RegionOrigin", "Unknown")
	set(dfs, NULL, "RegionOrigin", dfs[, as.character(factor(RegionOrigin))])
	set(dfs, NULL, 'Age_AnyPosT1', dfs[, AnyPos_T1-DateBorn])
	set(dfs, NULL, 'YEAR', dfs[, floor(AnyPos_T1)])		
	set(dfs, dfs[, which(Region_first=='Amst')], 'Region_first', 'North')
	set(dfs, NULL, 'GGD_first', dfs[, gsub('Hulpverlening_','',as.character(GGD_first))])
	set(dfs, dfs[, which(is.na(GGD_first))], "GGD_first", "Unknown")
	set(dfs, NULL, 'GGD_now', dfs[, gsub('Hulpverlening_','',as.character(GGD_now))])
	set(dfs, dfs[, which(is.na(GGD_now))], "GGD_now", "Unknown")
	set(dfs, dfs[, which(is.na(Region_first))], 'Region_first', 'Unknown')
	set(dfs, dfs[, which(Region_first%in%c('Rott','Den_Haag'))], 'Region_first', 'West')
	set(dfs, dfs[, which(GGD_first%in%c("Groningen", "Fryslan", "Drenthe"))], 'Region_first', 'North')
	set(dfs, dfs[, which(GGD_first%in%c("Flevoland","Twente","Gelre_IJssel","Midden_Nederland","Gelderland_Midden","Nijmegen","Rivierenland","IJsselland"))], 'Region_first', 'East')
	set(dfs, dfs[, which(GGD_first%in%c("Amsterdam", 'Den_Haag', "Zuid_Holland_West", "Zeeland", "Gooi_Vechtstreek", "Hollands_Noorden", "Rotterdam_Rijnmond", "Hollands_Midden", "Zaanstreek_Waterland", "Kennemerland", "Utrecht", "Zuid_Holland_Zuid"))], 'Region_first', 'West')
	set(dfs, dfs[, which(GGD_first%in%c("Hart_voor_Brabant","Zuid_Limburg","Brabant_Zuidoost","Limburg-Noord","West_Brabant"))], 'Region_first', 'South')
	#set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, breaks=c(-.01,15,seq(28,48,5),100), labels=c('0-15','16-27','28-32','33-37','38-42','43-47','48-80'))])
	set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,seq(15,65,10),100), labels=c('0-14','15-24','25-34','35-44','45-54','55-64','65-99'))])
	#
	#	get population by age year ggd 
	#
	load(infile.pop)	
	dpop	<- dpop[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD','AGE','YEAR')]
	set(dpop, dpop[, which(AGE%in%c('0-4','5-9','10-14'))], 'AGE', '0-14')
	set(dpop, dpop[, which(AGE%in%c("15-19","20-24"))], 'AGE', '15-24')
	set(dpop, dpop[, which(AGE%in%c("25-29","30-34"))], 'AGE', '25-34')
	set(dpop, dpop[, which(AGE%in%c("35-39","40-44"))], 'AGE', '35-44')
	set(dpop, dpop[, which(AGE%in%c("45-49","50-54"))], 'AGE', '45-54')
	set(dpop, dpop[, which(AGE%in%c("55-59","60-64"))], 'AGE', '55-64')
	set(dpop, dpop[, which(AGE%in%c('65-69',"70-74","75-79","80-84","85-89","90-94","95-99"))], 'AGE', '65-99')
	dpop	<- dpop[, list(M=sum(M), F=sum(F)), by=c('GGD_ID','GGD','AGE','YEAR')]
	dpop[, POP:=M+F]		
	tmp		<- data.table(GGD_= setdiff(dfs[,sort(unique(GGD_first))],'Unknown'), GGD=dpop[,sort(unique(GGD))])
	dpop	<- merge(dpop, tmp, by='GGD')
	set(dpop, NULL, 'GGD', NULL)
	setnames(dpop, 'GGD_', 'GGD')	
	tmp		<- unique(subset(dfs, select=c(Region_first, GGD_first)))
	setnames(tmp, c('Region_first','GGD_first'), c('REGION','GGD'))
	dpop	<- merge(dpop, tmp, by='GGD')
	#
	#	HIV by age year ggd
	#
	tmp		<- dfs[, list(NEW_DIAG= length(Patient)), by=c('GGD_first','YEAR','AGE','Trm')]
	tmp		<- rbind(tmp, tmp[, list(Trm='ALL', NEW_DIAG=sum(NEW_DIAG)), by=c('GGD_first','YEAR','AGE')])
	tmp		<- dcast.data.table(tmp, GGD_first+YEAR+AGE~Trm, value.var='NEW_DIAG')
	tmp		<- subset(tmp, YEAR>1999 & YEAR<2016 & GGD_first!='Unknown')
	setnames(tmp, c('GGD_first','MSM','IDU','OTH','HetM','HetF','Unknown','ALL'), c('GGD','DIAG_MSM','DIAG_IDU','DIAG_OTH','DIAG_HetM','DIAG_HetF','DIAG_Unknown','DIAG_ALL'))
	dpop	<- merge(dpop, tmp, by=c('GGD','YEAR','AGE'),all.x=1)
	for(x in colnames(dpop)[ grepl('^DIAG_',colnames(dpop)) ])
		set(dpop, which(is.na(dpop[[x]])), x, 0L)
	
	dpop[, DIAG_ALL_R:= DIAG_ALL/POP]	
	set(dpop, NULL, 'AGE', dpop[, factor(as.character(AGE), levels=rev(levels(dpop$AGE)))])
	
	#
	#	plot pop counts by GGD
	#
	tmp		<- dpop[, list(POP=sum(POP)),by=c('REGION','GGD','YEAR')]
	tmp2	<- subset(tmp, YEAR==2001)
	ggplot(tmp, aes(x=YEAR, y=POP)) +
			geom_line(aes(colour=GGD)) + geom_text(data=tmp2, aes(label=GGD), size=3, hjust=-.5) +
			facet_grid(~REGION) +		
			labs(x='', y='population\n(at Dec-31)', colour='Public health region') +
			theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(outfile,'pop_count.pdf'), w=10, h=10)
	
	ggplot(dpop, aes(x=YEAR, y=POP, fill=AGE)) +
			geom_bar(stat='identity') +
			scale_fill_viridis(option="magma", discrete=TRUE, end=0.85) +
			facet_wrap(~REGION+GGD, ncol=5) +
			labs(x='', y='population\n(at Dec-31)', fill='age group') +
			theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(outfile,'pop_by_age_count.pdf'), w=15, h=15)
	ggplot(dpop, aes(x=YEAR, y=POP, fill=AGE)) +
			geom_bar(stat='identity', position='fill') +
			scale_y_continuous(labels=percent, expand=c(0,0)) + 
			scale_fill_viridis(option="magma", discrete=TRUE, end=0.85) +
			facet_wrap(~REGION+GGD, ncol=5) +
			labs(x='', y='population\n(at Dec-31)', fill='age group') +
			theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(outfile,'pop_by_age_percent.pdf'), w=15, h=15)	
}
######################################################################################
eval.diag.newdiagnoses.by.migration<- function()
{
	require(viridis)
	require(data.table)
	require(ggplot2)
	infile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	outfile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	load(infile)
	#
	#	get first sequence
	#	
	dfs		<- df.all[, {
				z	<- NA_character_
				zz	<- NA_real_
				if(any(!is.na(PosSeqT)))
				{
					z	<- FASTASampleCode[which.min(PosSeqT)]
					zz	<- min(PosSeqT, na.rm=TRUE)
				}
				list(FASTASampleCode=z, PosSeqT=zz)
			}, by=c('Patient')]
	tmp		<- copy(df.all)
	set(tmp, NULL, c('PosSeqT','FASTASampleCode','SUBTYPE_C','SEQ_L'), NULL)
	dfs		<- merge(dfs, unique(tmp), by='Patient')
	
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	#set(dfs, dfs[, which(Sex=='M' & Trm=='Het')], 'Trm', 'HetM')
	#set(dfs, dfs[, which(Sex=='F' & Trm=='Het')], 'Trm', 'HetF')	
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG','SXCH'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unknown')
	set(dfs, dfs[, which(RegionOrigin%in%c("Austr_NZ","Western_EU"))], "RegionOrigin", "Western")
	set(dfs, dfs[, which(RegionOrigin%in%c("Caribbean","Latin_South_America"))], "RegionOrigin", "Carib_Southern_America")
	set(dfs, dfs[, which(RegionOrigin%in%c("Central_EU","Eastern_EU_stans"))], "RegionOrigin", "Cen_East_Stans_EU")
	set(dfs, dfs[, which(RegionOrigin%in%c("Sout_SouthEast_Asia"))], "RegionOrigin", "Other")	
	set(dfs, dfs[, which(RegionOrigin%in%c("Oceania_Pacific","North_Africa_Middle_East"))], "RegionOrigin", "Other")
	set(dfs, dfs[, which(is.na(RegionOrigin))], "RegionOrigin", "Unknown")
	set(dfs, NULL, "RegionOrigin", dfs[, as.character(factor(RegionOrigin))])
	set(dfs, NULL, 'ORIGIN_TRM', NA_character_)
	set(dfs, dfs[, which(RegionOrigin=='NL' & Trm=='MSM')], 'ORIGIN_TRM', 'NL_MSM')
	set(dfs, dfs[, which(RegionOrigin=='NL' & Trm!='MSM')], 'ORIGIN_TRM', 'NL_Oth')
	set(dfs, dfs[, which(RegionOrigin!='NL' & Trm=='MSM')], 'ORIGIN_TRM', 'Foreign_MSM')
	set(dfs, dfs[, which(RegionOrigin!='NL' & Trm!='MSM')], 'ORIGIN_TRM', 'Foreign_Oth')
	set(dfs, dfs[, which(RegionOrigin=='Unknown' & Trm=='MSM')], 'ORIGIN_TRM', 'UnknownOri_MSM')
	set(dfs, dfs[, which(RegionOrigin=='Unknown' & Trm!='MSM')], 'ORIGIN_TRM', 'UnknownOri_Oth')		
	set(dfs, NULL, 'Age_AnyPosT1', dfs[, AnyPos_T1-DateBorn])
	set(dfs, NULL, 'YEAR', dfs[, floor(AnyPos_T1)])		
	set(dfs, NULL, 'GGD_first', dfs[, gsub('-','_',gsub('Hulpverlening_','',as.character(GGD_first)))])
	set(dfs, dfs[, which(is.na(GGD_first))], "GGD_first", "Unknown")
	set(dfs, dfs[, which(is.na(Region_first))], 'Region_first', 'Unknown')
	
	#
	#	% diagnoses among migrants 2010-2015 in Amsterdam
	#
	dfa	<- subset(dfs, !Trm%in%c('Unknown','OTH','IDU') & RegionOrigin!='Unknown' & AnyPos_T1>=2010 & AnyPos_T1<2016 & GGD_first=='Amsterdam')
	dfa[, mean(RegionOrigin!='NL')]
	#
	#	% diagnoses among migrants 2010-2015 in Amsterdam, The Hague, Rotterdam
	#
	dfa	<- subset(dfs, !Trm%in%c('Unknown','OTH','IDU') & !RegionOrigin%in%c('Unknown','NL') & AnyPos_T1>=2010 & AnyPos_T1<2016)
	dfa[, mean(GGD_first%in%c('Amsterdam','Rotterdam_Rijnmond','Den_Haag'))]
	
	#
	#	cumulative numbers by sex orientation and origin in Amsterdam since 2010
	#
	dfa	<- subset(dfs, !Trm%in%c('Unknown','OTH','IDU') & RegionOrigin!='Unknown' & AnyPos_T1>=2010 & AnyPos_T1<2016 & GGD_first=='Amsterdam')
	dfa	<- dfa[, list(ND=length(Patient)), by=c('YEAR','Trm','RegionOrigin')]
	tmp	<- as.data.table(expand.grid(YEAR=dfa[, sort(unique(YEAR))], Trm=dfa[, sort(unique(Trm))], RegionOrigin=dfa[, sort(unique(RegionOrigin))], stringsAsFactors=FALSE))  
	dfa	<- merge(tmp, dfa, by=c('YEAR','Trm','RegionOrigin'), all.x=1)
	set(dfa, dfa[, which(is.na(ND))], 'ND', 0)
	ggplot(dfa, aes(x=YEAR, y=ND, colour=RegionOrigin)) + geom_line() + facet_grid(~Trm)
	#	extrapolate numbers for 2016 - 2020
	#	simply take 2014-2015 numbers for all but Dutch MSM
	tmp	<- subset(dfa, YEAR>2013)[,list(ND=0.7*mean(ND)),by=c('Trm','RegionOrigin')]
	set(tmp, tmp[, which(Trm=='MSM' & RegionOrigin=='NL')],'ND',40)
	tmp	<- merge(as.data.table(expand.grid(YEAR=2016:2020, Trm=dfa[, sort(unique(Trm))], RegionOrigin=dfa[, sort(unique(RegionOrigin))], stringsAsFactors=FALSE)), tmp, by=c('Trm','RegionOrigin'))
	dfa	<- rbind(dfa, tmp)
	ggplot(dfa, aes(x=YEAR, y=ND, colour=RegionOrigin)) + geom_line() + facet_grid(~Trm)
	#	add sequencing probs to calculate trm events that are in data
	probs<- data.table(	Trm=c('MSM','Het'), 
						IN_NL=c(0.75,0.5), 
						IN_A=0.7,
						SA=0.95,
						SO=0.5)
	dfa	<- merge(dfa, probs, by=c('Trm'))
	set(dfa, dfa[, which(RegionOrigin=='NL')], 'IN_NL', 0.9)
	dfa[, TRM_INDATA:= ND*IN_NL*IN_A*SA + ND*IN_NL*(1-IN_A)*SO]
	
	
	setkey(dfa, Trm, RegionOrigin, YEAR)
	tmp	<- dfa[, list(YEAR=YEAR, NDC=round(cumsum(ND)), TRM_INDATA_C=round(cumsum(TRM_INDATA))), by=c('Trm','RegionOrigin')]
	ggplot(tmp, aes(x=YEAR, y=TRM_INDATA_C, colour=RegionOrigin)) + geom_line() + facet_grid(~Trm) 
	tmp	<- dcast.data.table(melt(tmp, measure.vars=c('NDC','TRM_INDATA_C')), Trm+YEAR~RegionOrigin+variable, value.var='value')
	
	write.csv(tmp, row.names=FALSE, file='~/Dropbox (Infectious Disease)/OR_Work/2017/2017_AIDSFonds/trm_events_numbers.csv')
	
}
######################################################################################
eval.diag.newdiagnoses.by.age.migration<- function()
{
	require(viridis)
	infile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	outfile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	load(infile)
	#
	#	age group
	#	
	dfs		<- df.all[, {
				z	<- NA_character_
				zz	<- NA_real_
				if(any(!is.na(PosSeqT)))
				{
					z	<- FASTASampleCode[which.min(PosSeqT)]
					zz	<- min(PosSeqT, na.rm=TRUE)
				}
				list(FASTASampleCode=z, PosSeqT=zz)
			}, by=c('Patient')]
	tmp		<- copy(df.all)
	set(tmp, NULL, c('PosSeqT','FASTASampleCode','SUBTYPE_C','SEQ_L'), NULL)
	dfs		<- merge(dfs, unique(tmp), by='Patient')
	
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG','SXCH'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unkown')
	set(dfs, dfs[, which(is.na(Region_first))], 'Region_first', 'Unknown')
	set(dfs, dfs[, which(Region_first%in%c('Rott','Den_Haag'))], 'Region_first', 'West')
	set(dfs, dfs[, which(RegionOrigin%in%c("Central_EU","Eastern_EU_stans"))], "RegionOrigin", "Cen_East_Stans_EU")
	set(dfs, dfs[, which(RegionOrigin%in%c("Austr_NZ","Sout_SouthEast_Asia","Oceania_Pacific","North_Africa_Middle_East"))], "RegionOrigin", "Other")
	set(dfs, dfs[, which(is.na(RegionOrigin))], "RegionOrigin", "Unknown")
	set(dfs, NULL, "RegionOrigin", dfs[, as.character(factor(RegionOrigin))])
	set(dfs, NULL, 'Age_AnyPosT1', dfs[, AnyPos_T1-DateBorn])
	set(dfs, NULL, 'AnyPosT1_Y', dfs[, floor(AnyPos_T1)])	
	#set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, breaks=c(-.01,15,seq(28,48,5),100), labels=c('0-15','16-27','28-32','33-37','38-42','43-47','48-80'))])
	set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,seq(15,65,10),100), labels=c('0-14','15-24','25-34','35-44','45-54','55-64','65-80'))])	
	set(dfs, NULL, 'YOUNG', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,16, 28,100), labels=c('0-15', '16-27','28-80'))])
	set(dfs, NULL, 'OLD', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,50,100), labels=c('0-49', '50-80'))])
	
	tmp	<- subset(dfs, AnyPos_T1>2005 & AnyPos_T1<2016)
	set(tmp, NULL, 'TIME', tmp[, cut(AnyPos_T1, right=FALSE, breaks=c(0,2010,2020), labels=c('2005-2009','2010-2015'))])
	set(tmp, tmp[, which(RegionOrigin%in%c("Sub_Saharan_Africa","Unknown"))], "RegionOrigin", "Other")
	set(tmp, NULL, 'RegionOrigin', tmp[, factor(RegionOrigin, levels=c('NL','Caribbean','Cen_East_Stans_EU','Latin_South_America','Other','Western_EU'))])
	ggplot(tmp, aes(x=Age_AnyPosT1, fill=RegionOrigin)) + geom_histogram() +
			facet_grid(~Trm) +
			scale_x_continuous(breaks=seq(0, 90, 10)) +
			scale_y_continuous(breaks=seq(0,1e3,1e2)) +
			scale_fill_brewer(palette='PRGn') +			
			labs(x='\nage at diagnosis\namong newly diagnosed 2005-2015', y='new HIV+ diagnoses\n') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiag_by_ageatdiag_migrant_histo.pdf'), w=15, h=8)
	#
	#	how many young individuals are HIV+ by region origin?
	#
	ggplot(dfs, aes(x=floor(AnyPos_T1), fill=YOUNG)) +
			geom_bar() +
			scale_x_continuous(breaks=seq(1980, 2020, 5)) +
			#scale_y_continuous(breaks=seq(0,1,.1), labels=percent) +
			scale_fill_brewer(palette='Set1') +
			labs(x='\nyear of diagnosis', y='new HIV+ diagnoses\n') +
			facet_grid(RegionOrigin~.) +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiag_by_young_migrant.pdf'), w=5, h=15)
	#
	#	how many old individuals are HIV+ by region origin?
	#	
	ggplot(dfs, aes(x=floor(AnyPos_T1), fill=OLD)) +
			geom_bar() +
			scale_x_continuous(breaks=seq(1980, 2020, 5)) +
			#scale_y_continuous(breaks=seq(0,1,.1), labels=percent) +
			scale_fill_brewer(palette='Set1') +
			labs(x='\nyear of diagnosis', y='new HIV+ diagnoses\n') +
			facet_grid(RegionOrigin~.) +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiag_by_old_migrant.pdf'), w=5, h=15)
	#
	#	plot 10% 25% 50% 75% 90% quantiles over time by region origin for MSM
	#
	tmp	<- subset(dfs, AnyPos_T1>2000 & AnyPos_T1<2016 & Trm=='MSM')
	set(tmp, tmp[, which(RegionOrigin%in%c("Sub_Saharan_Africa","Unknown"))], "RegionOrigin", "Other")
	tmp	<- tmp[, list(qu20= quantile(Age_AnyPosT1, p=0.2), qu80= quantile(Age_AnyPosT1, p=0.8), median=median(Age_AnyPosT1)), by=c('RegionOrigin','AnyPosT1_Y')]
	ggplot(tmp, aes(x=AnyPosT1_Y)) +			
			geom_ribbon(aes(ymin=qu20, ymax=qu80,  fill=RegionOrigin), alpha=0.6, colour='transparent') +
			geom_line(aes(y=median, group=RegionOrigin)) + geom_point(aes(y=median, colour=RegionOrigin)) +
			scale_x_continuous(breaks=seq(1980, 2020, 5)) +
			facet_grid(~RegionOrigin) +
			labs(x='\nyear of diagnosis MSM', y='Age at diagnosis\n', linetype='quantile') +			
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiag_by_ageatdiag_migrant.pdf'), w=15, h=7)
	#
	#	cumulative new diagnoses by age group
	#
	tmp	<- subset(dfs, AnyPos_T1>2008)[, list(AnyPos_T1=sort(AnyPos_T1), CM=seq_along(AnyPos_T1)), by=c('Region_first','AGE')]
	ggplot(tmp, aes(x=AnyPos_T1, y=CM, colour=Region_first)) +
			geom_line() + 
			scale_x_continuous(breaks=seq(1980, 2020, 2)) +
			facet_grid(AGE~.) +
			labs(x='\nyear of diagnosis MSM', y='cumulative diagnoses since 2008\n', fill='age at diagnosis') +			
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'cumuldiag_by_ageatdiag.pdf'), w=8, h=12)	
	#subset(dfs, AnyPos_T1>2008 & GGD_first=='Rotterdam_Rijnmond' & AGE=='35-44')	
	tmp	<- subset(dfs, AnyPos_T1>2008 & Region_first=='West')[, list(AnyPos_T1=sort(AnyPos_T1), CM=seq_along(AnyPos_T1)), by=c('GGD_first','AGE')]
	ggplot(tmp, aes(x=AnyPos_T1, y=CM, colour=GGD_first)) +
			geom_line() + 
			scale_x_continuous(breaks=seq(1980, 2020, 2)) +
			facet_grid(AGE~.) +
			labs(x='\nyear of diagnosis MSM', y='cumulative diagnoses since 2008\n', fill='age at diagnosis') +			
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'cumuldiag_by_ageatdiag_ggd_West.pdf'), w=8, h=12)
	#	
	#	disproportional increases in 16-27 and 48-80 in West: any particular group?
	#
	ggplot(subset(dfs, AnyPos_T1>2008 & Region_first!='Unknown' ), aes(x=AGE, fill=RegionOrigin)) +
			geom_bar(position='fill') +			
			scale_y_continuous(breaks=seq(0,1,.2), labels=percent) +			
			scale_fill_brewer(palette='Spectral') +
			labs(x='\nage at diagnosis', y='new HIV+ diagnoses since 2008\n') +		
			facet_grid(Region_first~.) +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'diagsince2008_by_ageatdiag_region_percentage.pdf'), w=8, h=12)	
	ggplot(subset(dfs, AnyPos_T1>2008 & !is.na(Region_first)), aes(x=AGE, fill=RegionOrigin)) +
			geom_bar() +			
			scale_fill_brewer(palette='Spectral') +
			labs(x='\nage at diagnosis', y='new HIV+ diagnoses since 2008\n') +		
			facet_grid(Region_first~.) +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'diagsince2008_by_ageatdiag_region_number.pdf'), w=8, h=12)
	
	ggplot(subset(dfs, AnyPos_T1>2008 & !is.na(GGD_AnyPosT1) & Region_first=='West'), aes(x=AGE, fill=RegionOrigin)) +
			geom_bar() +			
			scale_fill_brewer(palette='Spectral') +
			labs(x='\nage at diagnosis', y='new HIV+ diagnoses since 2008\n') +		
			facet_grid(GGD_AnyPosT1~Region_first) +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'diagsince2008_by_ageatdiag_region_number_West.pdf'), w=8, h=12)
}
######################################################################################
eval.diag.newdiagnoses.by.age.migration.city<- function()
{
	require(viridis)
	infile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	outfile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	load(infile)
	#
	#	age group
	#	
	dfs		<- df.all[, {
				z	<- NA_character_
				zz	<- NA_real_
				if(any(!is.na(PosSeqT)))
				{
					z	<- FASTASampleCode[which.min(PosSeqT)]
					zz	<- min(PosSeqT, na.rm=TRUE)
				}
				list(FASTASampleCode=z, PosSeqT=zz)
			}, by=c('Patient')]
	tmp		<- copy(df.all)
	set(tmp, NULL, c('PosSeqT','FASTASampleCode','SUBTYPE_C','SEQ_L'), NULL)
	dfs		<- merge(dfs, unique(tmp), by='Patient')
	
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG','SXCH'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unkown')
	set(dfs, dfs[, which(is.na(Region_first))], 'Region_first', 'Unknown')	
	set(dfs, dfs[, which(RegionOrigin%in%c("Central_EU","Eastern_EU_stans"))], "RegionOrigin", "Cen_East_Stans_EU")
	set(dfs, dfs[, which(RegionOrigin%in%c("Austr_NZ","Sout_SouthEast_Asia","Oceania_Pacific","North_Africa_Middle_East"))], "RegionOrigin", "Other")
	set(dfs, dfs[, which(is.na(RegionOrigin))], "RegionOrigin", "Unknown")
	set(dfs, NULL, "RegionOrigin", dfs[, as.character(factor(RegionOrigin))])
	set(dfs, NULL, 'Age_AnyPosT1', dfs[, AnyPos_T1-DateBorn])
	set(dfs, NULL, 'AnyPosT1_Y', dfs[, floor(AnyPos_T1)])	
	#set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, breaks=c(-.01,15,seq(28,48,5),100), labels=c('0-15','16-27','28-32','33-37','38-42','43-47','48-80'))])
	set(dfs, NULL, 'Age_AnyPosT1', dfs[, AnyPos_T1-DateBorn])
	set(dfs, NULL, 'AGE', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,seq(15,65,10),100), labels=c('0-14','15-24','25-34','35-44','45-54','55-64','65-99'))])
	#set(dfs, NULL, 'AGE', dfs[, factor(as.character(AGE), levels=rev(levels(dfs$AGE)))])
	set(dfs, NULL, 'YOUNG', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,15,30,100), labels=c('0-14', '15-29','30-99'))])
	set(dfs, NULL, 'OLD', dfs[, cut(Age_AnyPosT1, right=FALSE, breaks=c(-.01,15,50,70,100), labels=c('0-14', '15-49', '50-69', '70-99'))])	
	#	
	#	disproportional increases in 16-27 and 48-80 in West: any particular group?
	#
	ggplot(subset(dfs, AnyPos_T1>2010 & AGE!='0-14' & Region_first!='Unknown' & Trm%in%c('MSM','Het')), aes(x=AGE, fill=RegionOrigin)) +
			geom_bar(position='fill') +			
			scale_y_continuous(breaks=seq(0,1,.2), labels=percent) +			
			scale_fill_brewer(palette='Spectral') +
			labs(x='\nage at diagnosis', y='new HIV+ diagnoses since 2010\n') +		
			facet_grid(Region_first~Trm) +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiagsince2010_by_age_regionorigin_percent.pdf'), w=12, h=12)	
	ggplot(subset(dfs, AnyPos_T1>2010 & AGE!='0-14' & Region_first!='Unknown' & Trm%in%c('MSM','Het')), aes(x=AGE, fill=RegionOrigin)) +
			geom_bar() +			
			scale_fill_brewer(palette='Spectral') +
			labs(x='\nage at diagnosis', y='new HIV+ diagnoses since 2008\n') +		
			facet_grid(Region_first~Trm) +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiagsince2010_by_age_regionorigin_number.pdf'), w=12, h=12)
	
	ggplot(subset(dfs, AnyPos_T1>2010 & !is.na(GGD_first) & Trm%in%c('MSM','Het') & Region_first%in%c('West','Amsterdam','Den_Haag','Rotterdam_Rijnmond','Utrecht')), aes(x=AGE, fill=RegionOrigin)) +
			geom_bar() +			
			scale_fill_brewer(palette='Spectral') +
			labs(x='\nage at diagnosis', y='new HIV+ diagnoses since 2008\n') +		
			facet_grid(GGD_first~Trm, scales='free_y') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'newdiagsince2010_by_age_regionorigin_number_West.pdf'), w=12, h=12)
}
######################################################################################
eval.seq.sequence.coverage.among.diagnosed<- function()
{
	infile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_PatientKeyCovariates_Numeric.rda'
	outfile	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/info/ATHENA_1610_'
	load(infile)
	#
	#	sequence coverage	
	#
	dfs		<- copy(df.all)	
	dfs		<- dfs[, {
				z	<- NA_character_
				zz	<- NA_real_
				if(any(!is.na(PosSeqT)))
				{
					z	<- FASTASampleCode[which.min(PosSeqT)]
					zz	<- min(PosSeqT, na.rm=TRUE)
				}
				list(FASTASampleCode=z, PosSeqT=zz, AnyPos_T1=AnyPos_T1[1], DateDied=DateDied[1])
			}, by=c('Patient','Trm','Region_first','SUBTYPE_C')]
	set(dfs, dfs[, which(Trm%in%c('MSM','BI'))], 'Trm', 'MSM')
	set(dfs, dfs[, which(Trm%in%c('HET','HETfa'))], 'Trm', 'Het')
	set(dfs, dfs[, which(Trm%in%c('BLOOD','NEEACC','PREG'))], 'Trm', 'OTH')
	set(dfs, dfs[, which(is.na(Trm))], 'Trm', 'Unkown')
	set(dfs, dfs[, which(is.na(Region_first))], 'Region_first', 'Unkown')		
	
	
	dfd	<- as.data.table(expand.grid(t=seq(1996.5, 2016, 0.125), Region=dfs[, unique(as.character(Region_first))], stringsAsFactors=FALSE))
	dfd	<- dfd[, list( 	N_DIAG=nrow(subset(dfs, AnyPos_T1<t & Region_first==Region & (is.na(DateDied) | DateDied>t))),
						N_SEQ=nrow(subset(dfs, PosSeqT<t & Region_first==Region & (is.na(DateDied) | DateDied>t)))		)
			, by=c('t','Region')]
	#
	#	sequence coverage by region
	#
	dfd[, NOT_SEQ:= N_DIAG-N_SEQ]
	dfd[, SEQ_C:= N_SEQ/N_DIAG]
	ggplot(dfd, aes(x=t, y=SEQ_C, colour=Region)) + geom_line() +
			scale_x_continuous(breaks=seq(1980, 2020, 2)) +
			scale_y_continuous(breaks=seq(0,1,.1), labels=percent) +			
			labs(x='\ndate of diagnosis', y='proportion of diagnosed and alive with sequence\n(years)\n', pch='Region', linetype='Region') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'seqcoverage_by_region.pdf'), w=5, h=5)	
	#
	#	sequence coverage overall
	#
	tmp	<- dfd[, list(SEQ_C= sum(N_SEQ)/sum(N_DIAG)), by='t']
	ggplot(tmp, aes(x=t, y=SEQ_C)) + geom_line() +
			scale_x_continuous(breaks=seq(1980, 2020, 2)) +
			scale_y_continuous(breaks=seq(0,1,.1), labels=percent) +			
			labs(x='\ndate of diagnosis', y='proportion of diagnosed and alive with sequence\n(years)\n') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'seqcoverage.pdf'), w=5, h=5)	
	#
	#	sequence coverage among those newly diagnosed in year t by region
	#
	dfd	<- as.data.table(expand.grid(t=seq(1996, 2016, 1), Region=dfs[, unique(as.character(Region_first))], stringsAsFactors=FALSE))
	dfd	<- dfd[, list( 	N_DIAG=nrow(subset(dfs, floor(AnyPos_T1)==t & Region_first==Region)),
						N_SEQ=nrow(subset(dfs, floor(AnyPos_T1)==t & !is.na(PosSeqT) & Region_first==Region))		)
			, by=c('t','Region')]
	dfd[, SEQ_C:= N_SEQ/N_DIAG]
	ggplot(subset(dfd, Region!='Other'), aes(x=factor(t), y=SEQ_C, fill=Region)) + 
			geom_bar(stat='identity') +
			scale_x_discrete(breaks=seq(1980, 2020, 2)) +
			scale_y_continuous(breaks=seq(0,1,.2), labels=percent) +
			scale_fill_brewer(palette='Set1') +
			facet_grid(Region~.) +
			labs(x='\nyear', y='proportion of individuals with a sequence\namong diagnosed HIV+ that were alive in given year\n', pch='Region', linetype='Region') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'seqcoverage_by_yeardiagnosed.pdf'), w=5, h=12)
	#
	#	seq coverage among individuals diagnosed after 2005 
	#
	dfd	<- as.data.table(expand.grid(t=seq(2005, 2016, 0.125), Region=dfs[, unique(as.character(Region_first))], stringsAsFactors=FALSE))
	dfd	<- dfd[, list( 	N_DIAG=nrow(subset(dfs, AnyPos_T1>=2005 & AnyPos_T1<t & Region_first==Region & (is.na(DateDied) | DateDied>t))),
					N_SEQ=nrow(subset(dfs, AnyPos_T1>=2005 & PosSeqT<t & Region_first==Region & (is.na(DateDied) | DateDied>t)))		)
			, by=c('t','Region')]
	dfd[, SEQ_C:= N_SEQ/N_DIAG]
	ggplot(dfd, aes(x=t, y=SEQ_C, colour=Region)) + geom_line() +
			scale_x_continuous(breaks=seq(1980, 2020, 2)) +
			scale_y_continuous(breaks=seq(0,1,.1), labels=percent) +			
			labs(x='\nyear', y='proportion of individuals with a sequence\namong HIV+ that were diagnosed after 2005 and were alive in given year\n', pch='Region', linetype='Region') +
			theme_bw() + theme(legend.position = "bottom")
	ggsave(file=paste(outfile,'seqcoverage2005_by_region.pdf'), w=5, h=5)	
	
}
