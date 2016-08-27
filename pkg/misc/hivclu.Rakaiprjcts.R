######################################################################################
project.Rakai.checkMissingRakai.150307<- function()
{
	png.f	<- '~/Dropbox (Infectious Disease)/PANGEA_data/2016-01-20_PANGEA_stats_by_sample.csv'	
	#infile	<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/CheckPangeaId.csv'
	#infile	<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/CheckPangeaId2.csv'
	infile	<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/SequenceDataNOTAvailable_Documented.csv'
	infile	<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/SequenceDataAvailable_NotDocumented.csv'
	
	#
	#
	#
	dpng	<- as.data.table(read.csv(png.f))
	setnames(dpng, 'ProjectID', 'PANGEA_ID')
	df		<- as.data.table(read.csv(infile))
	setnames(df, 'x', 'PANGEA_ID')
	df		<- merge(dpng, df, by='PANGEA_ID')	
	write.csv(df, row.names=FALSE, file=gsub('.csv','_info.csv',infile))
	
	png.f	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/PANGEA_orig/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	load(png.f)	#sq, sqi, si
	set(si, NULL, 'PANGEA_ID', si[, regmatches(PANGEA_ID, regexpr('PG[0-9]+-[^-]+',PANGEA_ID))])
	merge(si, df, by='PANGEA_ID')		
}
#
#	plots of individual histories
#	
RakaiCirc.circ.timelines.plots<- function(rt, wdir)
{
#	plot total circumcision by year among individuals in follow up
	setkey(rt, RID, ROUND)
	ggplot(unique(rt), aes(x=factor(ROUND_ST), fill=CIRC)) + geom_bar() + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1e4, 5e2)) +
			facet_grid(SEX~.) +
			labs(x='\nsurvey rounds\n(start date)', y='HIV infected RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Circumcision\nstatus')
	ggsave(file.path(wdir, 'RCCSinf_by_circ.pdf'), w=10, h=8)
	#	plot proportion circumcision by year among individuals in follow up
	ggplot(subset(unique(rt), SEX=='male'), aes(x=factor(ROUND_ST), fill=CIRC)) + geom_bar(position='fill') + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1, 0.1), labels = scales::percent) +
			labs(x='\nsurvey rounds\n(start date)', y='HIV infected male RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Circumcision\nstatus')
	ggsave(file.path(wdir, 'RCCSinf_by_circ_proportions_men.pdf'), w=10, h=5)	
	#	plot sequence coverage among HIV infected participants
	tmp		<- unique(rt)
	set(tmp, NULL, 'SEQ_TYPE', tmp[, factor(SEQ_TYPE, levels=c('gag_or_partial_gag_PANGEA','gag_or_partial_gag_historic_only','Sanger_completed_with_IVA','Sanger_completed_without_IVA','Sanger_not_started_by_Jul2016','other_gene_historic_only','Sanger_failed','no_sequence'))])
	ggplot(tmp, aes(x=factor(ROUND_ST), fill=SEQ_TYPE)) + geom_bar() + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1e4, 5e2)) +
			scale_fill_brewer(palette='Set2') +
			labs(x='\nsurvey rounds\n(start date)', y='HIV infected RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Sequence\nstatus')
	ggsave(file.path(wdir, 'RCCSinf_by_sequence_coverage.pdf'), w=10, h=5)
	#	plot proportion sequence coverage among HIV infected participants
	ggplot(tmp, aes(x=factor(ROUND_ST), fill=SEQ_TYPE)) + geom_bar(position='fill') + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1, 0.1), labels = scales::percent) +
			scale_fill_brewer(palette='Set3') +
			labs(x='\nsurvey rounds\n(start date)', y='HIV infected RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Sequence\nstatus')
	ggsave(file.path(wdir, 'RCCSinf_by_sequence_coverage_proportions.pdf'), w=10, h=5)	
	#	plot circumcised by religion
	tmp		<- subset(unique(rt), EVERCIRC=='Y')
	ggplot(tmp, aes(x=factor(ROUND_ST), fill=RELIGION)) + geom_bar(position='fill') + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1, 0.1), labels = scales::percent) +
			scale_fill_brewer(palette='Set2') +
			labs(x='\nsurvey rounds\n(start date)', y='circumcised male \nHIV infected RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Faith')
	ggsave(file.path(wdir, 'RCCScircmen_by_religion_proportions.pdf'), w=10, h=5)		
}
#
#	add first sequences
#	
RakaiCirc.circ.timelines.addfirstseq.160816<- function(rt, rs)
{
	tmp		<- rs[, {
						z	<- which.min(DATE)
						ps	<- NA_character_						
						if(any(SEQTYPE=='Sanger not started by Jul2016'))
							ps	<- 'Sanger not started by Jul2016'												
						if(any(SEQTYPE=='Sanger failed'))
							ps	<- 'Sanger failed'												
						if(any(SEQTYPE=='Sanger completed without IVA'))
							ps	<- 'Sanger completed without IVA'												
						if(any(SEQTYPE=='Sanger completed with IVA'))
							ps	<- 'Sanger completed with IVA'						
						if(any(SEQTYPE=='PANGEA'))
						{
							ps	<- 'PANGEA'
							z	<- which.min(DATE)
						}													
						list(	VISIT=VISIT[z], 
								DATE=DATE[z], 
								SEQ_HIST= any(SEQTYPE=='HISTORIC'),
								SEQ_HIST_GAG= ifelse(any(SEQTYPE=='HISTORIC'), any(GENE_REGION[SEQTYPE=='HISTORIC']=='Region1'), FALSE),
								SEQ_PANGEA_NOW= any(SEQTYPE=='PANGEA'),
								SEQ_PANGEA_NOW_GAG= ifelse(any(SEQTYPE=='PANGEA'), any(GENE_REGION[SEQTYPE=='PANGEA']=='Region1'), FALSE),
								SEQ_N=length(DATE),
								PANGEA_SEQTYPE=ps,								 
								SEQ_NOW_GAG= any(GENE_REGION=='Region1'))
					}, by=RID]
	#	complete to all available RIDs
	tmp		<- merge(unique(subset(rt, select=RID)), tmp, all.x=1, by='RID')			
	#	define SEQ_TYPE of interest
	tmp[, SEQ_TYPE:= 'no_sequence']
	tmp2	<- tmp[, which(grepl('Sanger',PANGEA_SEQTYPE))]
	set(tmp, tmp2, 'SEQ_TYPE', gsub(' ','_',tmp[tmp2, PANGEA_SEQTYPE]))	
	set(tmp, tmp[, which(!SEQ_HIST_GAG & SEQ_HIST & !SEQ_PANGEA_NOW)], 'SEQ_TYPE', "other_gene_historic_only")
	set(tmp, tmp[, which(SEQ_HIST_GAG & SEQ_HIST & !SEQ_PANGEA_NOW)], 'SEQ_TYPE', "gag_or_partial_gag_historic_only")
	set(tmp, tmp[, which(!SEQ_PANGEA_NOW_GAG & SEQ_PANGEA_NOW)], 'SEQ_TYPE', "other_gene_PANGEA")
	set(tmp, tmp[, which(SEQ_PANGEA_NOW_GAG & SEQ_PANGEA_NOW)], 'SEQ_TYPE', "gag_or_partial_gag_PANGEA")	
	#	
	tmp[, SEQ_NOW:= factor(SEQ_TYPE%in%c('gag_or_partial_gag_historic_only','gag_or_partial_gag_PANGEA','other_gene_historic_only'), levels=c(TRUE,FALSE), labels=c('Y','N'))]
	set(tmp, NULL, 'SEQ_NOW_GAG', tmp[, factor(SEQ_NOW=='Y' & !is.na(SEQ_NOW_GAG) & SEQ_NOW_GAG, levels=c(TRUE,FALSE), labels=c('Y','N'))])
	set(tmp, NULL, c('VISIT','SEQ_HIST','SEQ_PANGEA_NOW','SEQ_N', 'PANGEA_SEQTYPE','SEQ_HIST_GAG','SEQ_PANGEA_NOW_GAG'), NULL)
	setnames(tmp, 'DATE', 'FIRST_SEQ_DATE')
	rt		<- merge(rt, tmp, by='RID')
	rt
}
#
#	make individual timeline over visits
#
RakaiCirc.circ.timelines.init.160816<- function(rh, rd)
{	
	rt		<- rbind( subset(rh, select=c(RID, VISIT)), subset(rd, select=c(RID, VISIT)) )
	setkey(rt, RID, VISIT)
	rt		<- unique(rt)		
	cat('\nall visit data', nrow(rt), '\tall data in history',nrow(rh))
	rt[, ROUND:=VISIT]
	set(rt, rt[, which(ROUND==15.1)], 'ROUND', 15)
	#	expand to all possible rounds after entry into cohort
	tmp		<- as.data.table(expand.grid(RID=rt[, unique(RID)], ROUND=rt[, seq.int(1,max(ROUND))]))
	rt		<- merge(rt, tmp, by=c('RID','ROUND'),all=1)	
	#	add basic demographic data
	tmp		<- subset(rd, select=c(RID, SEX, RELIGION, CAUSE_OF_DEATH, BIRTHDATE, LASTNEGDATE, FIRSTPOSDATE, EST_DATEDIED))
	setkey(tmp, RID)
	tmp		<- unique(tmp)
	set(tmp, tmp[, which(is.na(RELIGION))], 'RELIGION', 'missing')
	set(tmp, NULL, 'SEX', tmp[, factor(SEX, levels=c('M','F'), labels=c('male','female'))])
	rt		<- merge(tmp, rt, by='RID')
	#	add start and end of visit rounds	
	tmp		<- rd[, list(ROUND_ST= round(min(DATE),d=1), ROUND_E= round(max(DATE),d=1) ), by='VISIT']
	setnames(tmp, 'VISIT','ROUND')
	#	TODO: rounds are overlapping? Here, I cut the end times of each period so that the periods do not overlap (hack).
	setkey(tmp, ROUND)
	set(tmp, NULL, 'ROUND_E2', c(tmp[, ROUND_ST-.1][-1], Inf) )
	tmp2	<- tmp[, which(ROUND_E2<ROUND_E)]
	set(tmp, tmp2, 'ROUND_E', tmp[tmp2,ROUND_E2])
	set(tmp, NULL, 'ROUND_E2', NULL)
	#	TODO: confirm no data from round 5?
	tmp		<- rbind(tmp, data.table(ROUND=5, ROUND_ST=1998.4, ROUND_E=1999.2))	
	rt		<- merge(rt, tmp, by='ROUND')
	setkey(rt, RID, ROUND)
	#	delete rounds before diagnosis or known visit
	rt		<- subset(rt, !is.na(VISIT) | (is.na(VISIT) & FIRSTPOSDATE<ROUND_E))
	stopifnot( nrow(subset(rt[, all(is.na(VISIT)), by=RID], V1))==0 )
	#	delete rounds after recorded death event
	rt		<- subset(rt, is.na(EST_DATEDIED) | (!is.na(EST_DATEDIED) & ROUND_ST<EST_DATEDIED))
	stopifnot( nrow(subset(rt[, all(is.na(VISIT)), by=RID], V1))==0 )
	#	define missing rounds as 
	#		'no_follow_up'			temporarily not in follow up, visits before and after the current round
	#		'before_first_visit'	already estimated to be HIV+, first visit after the current round
	#		'no_obs_est_died'		after last visit, estimated to have died in the current round
	#		'not_seen_again'		after last visit and no record of death
	tmp		<- rt[, list(ROUND_FIRST_ST=min(ROUND_ST[!is.na(VISIT)]), ROUND_LAST_E=max(ROUND_E[!is.na(VISIT)])), by='RID']	
	rt		<- merge(rt, tmp, by='RID')
	set(rt, NULL, 'VISIT', rt[, as.character(VISIT)])
	set(rt, rt[, which(is.na(VISIT) & ROUND_ST>ROUND_FIRST_ST & ROUND_E<ROUND_LAST_E)], 'VISIT', 'no_follow_up')	
	set(rt, rt[, which(is.na(VISIT) & ROUND_ST<ROUND_FIRST_ST)], 'VISIT', 'before_first_visit')
	set(rt, rt[, which(is.na(VISIT) & ROUND_E==ROUND_LAST_E)], 'VISIT', 'no_obs_est_died')
	set(rt, rt[, which(is.na(VISIT) & ROUND_E>ROUND_LAST_E)], 'VISIT', 'not_seen_again')
	stopifnot( nrow(subset(rt, is.na(VISIT)))==0 )
	#	define first circumcision 
	tmp		<- subset(rh, CIRC=='Y')[, list(FIRSTCIRC_V= min(VISIT) ), by='RID']
	tmp2	<- subset(rh, CIRC=='N')[, list(LASTNOTCIRC_V= max(VISIT) ), by='RID']
	tmp		<- merge(tmp, tmp2, by='RID', all=1)	
	#	TODO individuals with conflicting circumcision data. For now ignore conflicts (hack).
	tmp2	<- subset(tmp, FIRSTCIRC_V<=LASTNOTCIRC_V)
	write.csv(tmp2, row.names=FALSE, file=file.path(wdir, 'check_conflictingCircumcisionStatus.csv'))
	cat('\nWarning: Found conflicting circumcision data for n=',nrow(tmp2))
	set(tmp, tmp[, which(FIRSTCIRC_V<=LASTNOTCIRC_V)], 'LASTNOTCIRC_V', NA_real_)
	#	define ever circumcised and all missing
	tmp2	<- rh[, list(EVERCIRC=as.integer(any(CIRC=='Y', na.rm=1)), NODATA=all(is.na(CIRC))), by='RID']
	tmp		<- merge(tmp, tmp2, by='RID', all=1)	
	set(tmp, tmp[, which(NODATA)], 'EVERCIRC', 2L)
	#	tmp[, table(EVERCIRC, NODATA, useNA='if')]
	set(tmp, NULL, 'EVERCIRC', tmp[, factor(EVERCIRC, levels=c(0,1,2), labels=c('N','Y','missing_all'))])
	set(tmp, NULL, 'NODATA', NULL)
	rt		<- merge(rt, tmp, by='RID')
	#	define circumcision timeline
	rt[, CIRC:=NA_character_]	
	set(rt, rt[, which(FIRSTCIRC_V<=ROUND)], 'CIRC', 'Y')
	set(rt, rt[, which(ROUND<=LASTNOTCIRC_V)], 'CIRC', 'N')
	set(rt, rt[, which(EVERCIRC=='missing_all')], 'CIRC', 'missing_all')
	set(rt, rt[, which(EVERCIRC=='N')], 'CIRC', 'N')
	set(rt, rt[, which(is.na(CIRC) & FIRSTCIRC_V>ROUND)], 'CIRC', 'missing_before_circ')
	set(rt, rt[, which(is.na(CIRC) & LASTNOTCIRC_V<ROUND)], 'CIRC', 'missing_after_notcirc')
	stopifnot( nrow(subset(rt, is.na(CIRC)))==0 )
	set(rt, NULL, c('FIRSTCIRC_V','LASTNOTCIRC_V'), NULL)
	rt
}

RakaiCirc.recipient.get.info<- function(wdir=NA)
{
	infile				<- "~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda"
	indir.historicseq	<- "~/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree"
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	#
	#	prepare RCCS data
	#	
	load(infile)
	#	a bit of clean up 
	rd		<- as.data.table(rccsData)
	setnames(rd, colnames(rd), gsub('\\.','_',toupper(colnames(rd))))	
	for(x in colnames(rd))
		if(class(rd[[x]])=='Date')
			set(rd, NULL, x, hivc.db.Date2numeric(rd[[x]]))
	rh		<- as.data.table(rccsHistory)
	setnames(rh, colnames(rh), gsub('\\.','_',toupper(colnames(rh))))
	#rd[, table(VISIT)]
	#	make shorter
	setnames(rd, 'RCCS_STUDYID', 'RID')
	setnames(rd, 'PANGEA_ID', 'PID')
	setnames(rd, 'STUDYID', 'SID')
	setnames(rh, 'RCCS_STUDYID', 'RID')	
	
	rrec	<- subset(rd, SEX=='F') 
	stopifnot( nrow(subset(rrec, is.na(LASTNEGVIS) & !is.na(LASTNEGDATE)))==0 )
	stopifnot( nrow(subset(rrec, !is.na(LASTNEGVIS) & is.na(LASTNEGDATE)))==0 )
	rrec	<- subset(rrec, !is.na(LASTNEGVIS), select=c(RID, SEX, VISIT, DATE, REGION, COMM_NUM, HH_NUM, LASTNEGVIS, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE))
	stopifnot( nrow(subset(rrec, FIRSTPOSDATE<LASTNEGDATE))==0 )
	rrec[, SC_WINDOW:= FIRSTPOSDATE-LASTNEGDATE]		
	#	add first sequences. load sequence info. expect "rs".
	#	(if not exist, run RakaiCirc.seq.get.info() )
	load(file.path(wdir,'RCCS_SeqInfo_160816.rda'))	
	rs		<- subset(rs, !is.na(VISIT))	
	rrec	<- RakaiCirc.circ.timelines.addfirstseq.160816(rrec, rs)	
	#	rm B009659 (check) has missing first pos date, was sequenced, but Sanger_completed_without_IVA
	rrec	<- subset(rrec, RID!='B009659')
	#	rm B106184 (check) sequence date before last neg date 
	rrec	<- subset(rrec, RID!='B106184')
	stopifnot( nrow(subset(rrec, is.na(FIRSTPOSDATE)))==0 )
	stopifnot( nrow(subset(rrec, FIRSTPOSDATE<=LASTNEGDATE))==0 )
	#	make unique by RID
	set(rrec, NULL, c('VISIT','DATE'), NULL)
	rrec	<- unique(rrec)
	if(!is.na(wdir))
	{
		#	plot recipient sequence coverage
		tmp		<- subset(rrec, SC_WINDOW<3.4)
		tmp[, SC_WINDOWc:= cut(SC_WINDOW, breaks=seq(0,3.4,.2))]		
		set(tmp, NULL, 'SEQ_TYPE', tmp[, factor(SEQ_TYPE, levels=c('gag_or_partial_gag_PANGEA','gag_or_partial_gag_historic_only','Sanger_completed_with_IVA','Sanger_completed_without_IVA','Sanger_not_started_by_Jul2016','other_gene_historic_only','Sanger_failed','no_sequence'))])                           
		tmp		<- tmp[, list(N=length(RID)), by=c('SC_WINDOWc','SEQ_TYPE')]
		tmp		<- merge( tmp, data.table(expand.grid(SC_WINDOWc= tmp[, unique(SC_WINDOWc)], SEQ_TYPE= tmp[, unique(SEQ_TYPE)])), by=c('SC_WINDOWc','SEQ_TYPE'), all=1)
		set(tmp, tmp[, which(is.na(N))], 'N', 0)	
		setkey(tmp, SEQ_TYPE, SC_WINDOWc)
		tmp		<- tmp[, list(CN= cumsum(N), SC_WINDOWc=SC_WINDOWc), by='SEQ_TYPE']	
		ggplot(tmp, aes(x=SC_WINDOWc, y=CN, fill=SEQ_TYPE)) + geom_bar(position='stack', stat='identity') + 			
				scale_y_continuous(breaks=seq(0,5e3,2e2)) +
				scale_fill_brewer(palette='Set3') + 
				theme_bw() + theme(legend.position='bottom') + guides(fill=guide_legend(ncol=2)) +
				labs(x='\nseroconversion window\n(years)', y='female RCCS participants\nwith last negative test\n(cumulated counts)\n')
		ggsave(file.path(wdir, 'RCCSfemalerecipient_by_sequencecoverage.pdf'), w=14, h=8)
	}	
	rrec
}

RakaiCirc.seq.get.info<- function()
{
	#	TODO:
	#GENE_REGION
	#SEQTYPE                         Region1 Region2 Region3 Region4
	#HISTORIC                         2381       0       0    2607
	#PANGEA                           2826    2826    2826       0
	
	infile				<- "~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda"	
	infile.sangerstats	<- "~/Dropbox (Infectious Disease)/PANGEA_data/2016-07-07_PANGEA_stats_by_sample.csv"
	infile.relabel		<- "~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional phylogenetic analyses/Region 1 gag analysis/SummaryofGagSequenceData.rda"
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	#
	#	read all processed RCCS sequences 
	#	
	infile.region1		<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree/Prelim_RakaiPangeaSqnAndMetaData_IqTree_Region1.rda'
	infile.region2		<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree/Prelim_RakaiPangeaSqnAndMetaData_IqTree_Region2.rda'
	infile.region3		<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree/Prelim_RakaiPangeaSqnAndMetaData_IqTree_Region3.rda'
	infile.region4		<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree/Prelim_RakaiPangeaSqnAndMetaData_IqTree_Region4.rda'
	
	load(infile.region1)
	tmp		<- as.character(gag.sqn)
	tmp		<- data.table(	relabel= rownames(tmp),
							BASE_N= apply(tmp, 1, function(x) sum(as.numeric(x%in%c('a','t','g','c')))) )
	tmp		<- merge(tmp, as.data.table(subset(summaryData, select=c('seqid', 'relabel'))), by='relabel', all.x=1)
	#	get Pangea data from region 1
	z2		<- subset(as.data.table(pangeaMetaData), select=c(RCCS_studyid, visit, date, seqid, CometSubtype, Pangea.id))	
	z2		<- merge(z2, tmp, by='seqid', all.x=1)
	stopifnot( nrow(subset(z2, is.na(Pangea.id)))==0 )
	z2		<- subset(z2, BASE_N>0)
	z2[, SEQTYPE:='PANGEA']
	z2[, GENE_REGION:= 'Region1']
	rs		<- copy(z2)
	#	add historical data from region 1
	z2		<- subset(as.data.table(rakhistMetaData), select=c(RCCS_studyid, visit, date, seqid, CometSubtype))	
	z2		<- merge(z2, tmp, by='seqid', all.x=1)
	stopifnot( nrow(subset(z2, is.na(RCCS_studyid)))==0 )
	z2		<- subset(z2, BASE_N>0)
	z2[, SEQTYPE:='HISTORIC']
	z2[, GENE_REGION:= 'Region1']
	rs		<- rbind(rs, z2, use.names=TRUE, fill=TRUE)
	#	
	load(infile.region2)
	tmp		<- as.character(pol.sqn)
	tmp		<- data.table(	relabel= rownames(tmp),
			BASE_N= apply(tmp, 1, function(x) sum(as.numeric(x%in%c('a','t','g','c')))) )
	tmp		<- merge(tmp, as.data.table(subset(summaryData, select=c('seqid', 'relabel'))), by='relabel', all.x=1)
	#	get Pangea data from region 2
	z2		<- subset(as.data.table(pangeaMetaData), select=c(RCCS_studyid, visit, date, seqid, CometSubtype, Pangea.id))	
	z2		<- merge(z2, tmp, by='seqid', all.x=1)
	stopifnot( nrow(subset(z2, is.na(Pangea.id)))==0 )
	z2		<- subset(z2, BASE_N>0)
	z2[, SEQTYPE:='PANGEA']
	z2[, GENE_REGION:= 'Region2']
	rs		<- rbind(rs, z2, use.names=TRUE, fill=TRUE)
	#	no historical data from region 2
	#	
	load(infile.region3)
	#	TODO: there are 729 duplicates in env.sqn ie RCCS*A108646*vis15*reg13*com774*hh75*F*17*prev*2012.37
	tmp		<- as.character(env.sqn)
	tmp		<- data.table(	relabel= rownames(tmp),
							BASE_N= apply(tmp, 1, function(x) sum(as.numeric(x%in%c('a','t','g','c')))) )
	tmp		<- subset(tmp, BASE_N>0)	
	setkey(tmp, relabel)	
	tmp		<- merge(unique(tmp), as.data.table(subset(summaryData, select=c('seqid', 'relabel'))), by='relabel', all.x=1)
	#	get Pangea data from region 3
	z2		<- subset(as.data.table(pangeaMetaData), select=c(RCCS_studyid, visit, date, seqid, CometSubtype, Pangea.id))	
	z2		<- merge(z2, tmp, by='seqid', all.x=1)
	stopifnot( nrow(subset(z2, is.na(Pangea.id)))==0 )
	z2		<- subset(z2, BASE_N>0)
	z2[, SEQTYPE:='PANGEA']
	z2[, GENE_REGION:= 'Region3']
	rs		<- rbind(rs, z2, use.names=TRUE, fill=TRUE)
	#	no historical data from region 2
	#
	load(infile.region4)
	tmp		<- as.character(gp.sqn)
	tmp		<- data.table(	relabel= rownames(tmp),
							BASE_N= apply(tmp, 1, function(x) sum(as.numeric(x%in%c('a','t','g','c')))) )
	tmp		<- merge(tmp, as.data.table(subset(summaryData, select=c('seqid', 'relabel'))), by='relabel', all.x=1)
	#
	#	TODO: no PANGEA data from region 4
	#
	#	add historical data from region 4
	z2		<- subset(as.data.table(rakhistMetaData), select=c(RCCS_studyid, visit, date, seqid, CometSubtype))	
	z2		<- merge(z2, tmp, by='seqid', all.x=1)
	stopifnot( nrow(subset(z2, is.na(RCCS_studyid)))==0 )
	z2		<- subset(z2, BASE_N>0)
	z2[, SEQTYPE:='HISTORIC']
	z2[, GENE_REGION:= 'Region4']
	rs		<- rbind(rs, z2, use.names=TRUE, fill=TRUE)
	#
	#	reading done 
	#	clean up names etc
	#
	setnames(rs, colnames(rs), gsub('\\.','_',toupper(colnames(rs))))
	setnames(rs, c('RCCS_STUDYID','PANGEA_ID'), c('RID','PID'))
	set(rs, NULL, 'DATE', hivc.db.Date2numeric(rs[['DATE']]))
	set(rs, NULL, 'RELABEL', NULL)
	rs		<- subset(rs, !is.na(RID))		
	#
	#	add unprocessed but shipped PANGEA sequences
	#	
	load(infile)
	#	a bit of clean up 
	rd		<- as.data.table(rccsData)
	setnames(rd, colnames(rd), gsub('\\.','_',toupper(colnames(rd))))	
	for(x in colnames(rd))
		if(class(rd[[x]])=='Date')
			set(rd, NULL, x, hivc.db.Date2numeric(rd[[x]]))
	#	make shorter
	setnames(rd, 'RCCS_STUDYID', 'RID')
	setnames(rd, 'PANGEA_ID', 'PID')	
	#	
	tmp		<- subset(rd, !is.na(PID), select=c(RID, VISIT, DATE, PID))
	tmp2	<- setdiff(tmp[, PID], subset(rs, !is.na(PID))[, unique(PID)])
	cat('\nWarning: Found PANGEA sequences not in rd (error)? n=', length(setdiff(subset(rs, !is.na(PID))[, unique(PID)], tmp[, PID] )))	
	cat('\nWarning: PANGEA sequences not yet processed n=', length(tmp2))	
	setkey(tmp, RID)	
	#	TODO out of curiosity why are these 349 duplicated?
	cat('\nFound participants with more than one PANGEA sequence. n=',tmp[duplicated(tmp),][, length(unique(RID))],'\nids=',tmp[duplicated(tmp),][, paste(unique(RID), collapse=' ')])
	#	TODO there are 2674 unprocessed PANGEA sequences
	tmp		<- merge(tmp, data.table(PID=tmp2), by='PID')	
	tmp[, SEQTYPE:='PANGEA_not_yet_processed']	
	rs		<- rbind(rs, tmp, fill=TRUE, use.names=TRUE)
	#
	#	add 160120 statistics from Sanger
	#
	tmp		<- as.data.table(read.csv(infile.sangerstats, stringsAsFactors=FALSE))	
	setnames(tmp, colnames(tmp), gsub('\\.','_',toupper(colnames(tmp))))
	setnames(tmp, c('PROJECTID','STATUS'), c('PID','SANGER_STATUS'))
	tmp		<- subset(tmp, select=c(PID, SANGER_STATUS, ASSEMBLED, HIVCONTIG))
	tmp		<- tmp[-nrow(tmp),]
	set(tmp, tmp[, which(HIVCONTIG=='HIVcontig')],  'SANGER_STATUS',  'Sanger completed with IVA')
	set(tmp, tmp[, which(SANGER_STATUS=='Assume no HIV')], 'SANGER_STATUS', 'Sanger completed without IVA')
	set(tmp, tmp[, which(SANGER_STATUS=='Assume assembly failed')], 'SANGER_STATUS', 'Sanger failed')
	set(tmp, tmp[, which(SANGER_STATUS=='Assume sequencing failed')], 'SANGER_STATUS', 'Sanger failed')
	set(tmp, tmp[, which(SANGER_STATUS=='No sample')], 'SANGER_STATUS', 'Sanger failed')
	set(tmp, tmp[, which(SANGER_STATUS=='ID not in system yet')], 'SANGER_STATUS', 'Sanger not started by Jul2016')
	rs		<- merge(rs, subset(tmp, select=c(PID, SANGER_STATUS)), by='PID', all.x=1)
	tmp		<- rs[, which(SEQTYPE=='PANGEA_not_yet_processed')]
	set(rs, tmp, 'SEQTYPE', rs[tmp, SANGER_STATUS])
	set(rs, rs[, which(is.na(SEQTYPE))], 'SEQTYPE', 'Sanger not started by Jul2016')
	set(rs, NULL, 'SANGER_STATUS', NULL)
	#
	#	add Phylo ID for gag tree
	#
	load(infile.relabel)
	tmp		<- subset(as.data.table(summaryData), select=c(seqid,relabel))
	setnames(tmp, c('seqid','relabel'),c('SEQID','SEQIDb'))
	rs		<- merge(rs, tmp, by='SEQID',all.x=1)
	#	> rs[, table(SEQTYPE, GENE_REGION)]
	#                               GENE_REGION
	#	SEQTYPE                         Region1 Region2 Region3 Region4 <NA>
  	#	HISTORIC                         2381       0       0    2608    0
  	#	PANGEA                           2905    2401    1302       0    0
  	#	Sanger completed with IVA           0       0       0       0  598
  	#	Sanger completed without IVA        0       0       0       0  286
  	#	Sanger failed                       0       0       0       0  292
  	#	Sanger not started by Jul2016       0       0       0       0 1501
	#
	#	save to file
	#
	save(rs, file=file.path(wdir,'RCCS_SeqInfo_160816.rda'))
}

RakaiCirc.various<- function()
{
	require(ape)
	require(data.table)
	
	wdir	<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	wdir	<- '/work/or105/Gates_2014/Rakai'
	load(file.path(wdir,'RCCS_PhInfo_160825.rda'))
	
	#	add index for speed
	tmp			<- data.table(SEQIDc=rownames(sq), IDX=seq_len(nrow(sq)))
	gdgag		<- merge(tmp, gdgag, by='SEQIDc')
	setnames(tmp, c('SEQIDc','IDX'), c('SEQIDc2','IDX2'))
	gdgag		<- merge(tmp, gdgag, by='SEQIDc2')
	#	
	tmp			<- (as.character(sq)!='-')
	setkey(gdgag, IDX, IDX2)
	tmp2		<- gdgag[, list(OVERLAP= sum(apply(tmp[c(IDX,IDX2),],2,all)) ), by=c('IDX','IDX2')]
	gdgag		<- merge(gdgag, tmp2, by=c('SEQIDc','SEQIDc2'))
	
	save(phf, php, phfi, phpi, sq, gdgag, file=file.path(wdir,'RCCS_PhInfo_2_160825.rda'))
}

RakaiCirc.seq.get.phylogenies<- function()
{
	require(ape)
	require(data.table)
	wdir	<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	#	start with Kate s GP24 FastTree
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional phylogenetic analyses/Region 1 gag analysis/gagAllSmall.nwk'
	phf		<- read.tree(infile)
	phfi	<- data.table(SEQIDb= phf$tip.label)
	#	get patristic distances
	tmp								<- cophenetic.phylo(phf)
	tmp								<- as.matrix(tmp)
	tmp[upper.tri(tmp, diag=TRUE)]	<- NA_real_
	tmp								<- as.data.table(melt(tmp))
	setnames(tmp, c('Var1','Var2','value'), c('SEQIDb','SEQIDb2','PD'))
	phfi							<- subset(tmp, !is.na(PD))
	#	also use the PANGEA ExaML tree
	infile	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_fasttree.rda'
	load(infile)	#loads "ph", "dist.brl", "ph.gdtr", "ph.mrca", "clustering"
	php		<- ph
	tmp		<- copy(ph.gdtr)
	tmp[upper.tri(tmp, diag=TRUE)]	<- NA_real_
	tmp								<- as.data.table(melt(tmp))
	setnames(tmp, c('Var1','Var2','value'), c('PID','PID2','PD'))
	phpi							<- subset(tmp, !is.na(PD))
	#	get raw genetic distances on latest Region1 alignment
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1UG_codonaligned_p_PANGEA151113_p_HXB2.rda'
	load(infile)
	
	sq			<- as.character(sq)
	sq[sq=='?']	<- '-'
	sq			<- as.DNAbin(sq)
	tmp			<- dist.dna(sq, model='raw', pairwise.deletion=TRUE)
	tmp								<- as.matrix(tmp)
	tmp[upper.tri(tmp, diag=TRUE)]	<- NA_real_
	tmp								<- as.data.table(melt(tmp))
	setnames(tmp, c('Var1','Var2','value'), c('SEQIDc','SEQIDc2','GDRW'))
	gdgag		<- subset(tmp, !is.na(GDRW))
	set(gdgag, NULL, 'SEQIDc', gdgag[, as.character(SEQIDc)])
	set(gdgag, NULL, 'SEQIDc2', gdgag[, as.character(SEQIDc2)])
	#	add index for speed
	tmp			<- data.table(SEQIDc=rownames(sq), IDX=seq_len(nrow(sq)))
	gdgag		<- merge(tmp, gdgag, by='SEQIDc')
	setnames(tmp, c('SEQIDc','IDX'), c('SEQIDc2','IDX2'))
	gdgag		<- merge(tmp, gdgag, by='SEQIDc2')
	#	add overlap
	tmp			<- (as.character(sq)!='-')
	setkey(gdgag, IDX, IDX2)
	tmp2		<- gdgag[, list(OVERLAP= sum(apply(tmp[c(IDX,IDX2),],2,all)) ), by=c('IDX','IDX2')]
	gdgag		<- merge(gdgag, tmp2, by=c('SEQIDc','SEQIDc2'))
	
	
	save(phf, php, phfi, phpi, sq, gdgag, file=file.path(wdir,'RCCS_PhInfo_160825.rda'))	
}

RakaiCirc.circ.dev160815<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"	
	infile				<- "~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda"
	load(infile)
	#	a bit of clean up 
	rd		<- as.data.table(rccsData)
	setnames(rd, colnames(rd), gsub('\\.','_',toupper(colnames(rd))))	
	for(x in colnames(rd))
		if(class(rd[[x]])=='Date')
			set(rd, NULL, x, hivc.db.Date2numeric(rd[[x]]))
	rh		<- as.data.table(rccsHistory)
	setnames(rh, colnames(rh), gsub('\\.','_',toupper(colnames(rh))))
	#rd[, table(VISIT)]
	#	make shorter
	setnames(rd, 'RCCS_STUDYID', 'RID')
	setnames(rd, 'PANGEA_ID', 'PID')
	setnames(rd, 'STUDYID', 'SID')
	setnames(rh, 'RCCS_STUDYID', 'RID')	
	#	data checks
	setkey(rh, VISIT, RID)
	stopifnot(nrow(rh)==nrow(unique(rh)))
	#	define circumcision	
	set(rh, rh[, which(!CIRC%in%c(1,2))], 'CIRC', NA_integer_)
	set(rh, NULL, 'CIRC', rh[, factor(CIRC, levels=c(1,2), labels=c('Y','N'))])
	#	TODO Warning: found 1 female circumcised A006734
	tmp		<- rh[, which(CIRC=='Y' & SEX=='F')]
	cat('\nWarning: found female circumcised',rh[tmp, paste(RID, collapse=' ')])
	set(rh, tmp, 'CIRC', NA_integer_)
	#
	#	make individual timeline over visits
	#
	rt	<- RakaiCirc.circ.timelines.init.160816(rh, rd)
	#
	#	add first sequences
	#	load sequence info. expect "rs".
	#	(if not exist, run RakaiCirc.seq.get.info() )
	#
	load(file.path(wdir,'RCCS_SeqInfo_160816.rda'))	
	#	TODO: has not visit: E21593M2
	rs		<- subset(rs, !is.na(VISIT))
	rt		<- RakaiCirc.circ.timelines.addfirstseq.160816(rt, rs)
	setkey(rt, RID, ROUND)
	#RakaiCirc.circ.timelines.plots(rt, wdir)
	
	rrec	<- RakaiCirc.recipient.get.info(wdir=NA)
	
	#	check how many rec in tree, how many prob transmitters, and in which sequ sample
	tmp		<- RakaiCirc.seq.get.phylogenies()
	phf		<- tmp$phf
	php		<- tmp$php
	phfi	<- tmp$phfi
	phpi	<- tmp$phpi
	
	pd.max	<- 0.02
	#	may have multiple seqs per female recipient, keep all
	rrecs	<- merge(rrec, unique(subset(rs, select=c(RID, SEQID, SEQIDb))), by='RID', all.x=1)
	#	phfi is only lower triangular, complete 
	tmp		<- subset(phfi, PD<=pd.max)	
	tmp2	<- copy(tmp)	
	setnames(tmp2, c('SEQIDb','SEQIDb2'), c('SEQIDb2','SEQIDb'))	
	tmp		<- rbind(tmp, tmp2, use.names=TRUE)
	set(tmp, NULL, 'SEQIDb', tmp[, as.character(SEQIDb)])
	set(tmp, NULL, 'SEQIDb2', tmp[, as.character(SEQIDb2)])
	
	rrecs	<- merge(rrecs, tmp, all=1, by='SEQIDb')
	rrecs	<- subset(rrecs, !is.na(RID))
	#	add basic info on transmitters
	setnames(rrecs, 'SEQIDb2', 'SEQIDb_TR') 
	tmp		<- unique(subset(rs, !is.na(SEQIDb), select=c(RID, SEQIDb, SEQTYPE)))
	tmp2	<- subset(rd, select=c(RID, SEX))
	setkey(tmp2, RID)
	tmp		<- merge(tmp, unique(tmp2), by='RID')
	setnames(tmp, c('RID', 'SEQIDb', 'SEQTYPE', 'SEX'), c('RID_TR', 'SEQIDb_TR', 'SEQTYPE_TR', 'SEX_TR'))	
	rp		<- merge(rrecs, tmp, by='SEQIDb_TR', all.x=1)
	
	#	plot tree with pairs highlighted
	tmp		<- data.table(SEQIDb=phf$tip.label, IDX=seq_len(Ntip(phf)), TIP_COL='black')
	tmp2	<- data.table(SEQIDb=subset(rp, !is.na(SEQIDb_TR))[, unique(na.omit(c(SEQIDb, SEQIDb_TR))) ], PAIR='Y')	
	tmp		<- merge(tmp, tmp2, by='SEQIDb', all.x=1)
	set(tmp, tmp[, which(PAIR=='Y')], 'TIP_COL', 'red')
	setkey(tmp, IDX)
	file	<- file.path(wdir, paste('150825_FastTree_gagAllSmall_pairs_patristic',pd.max*100,'pc.pdf',sep=''))
	invisible(hivc.clu.plot(phf, tip.color=tmp[,TIP_COL], file=file, pdf.scaley=100, show.tip.label=TRUE, pdf.width=30))
	
	#	for every female recipient, count prob transmitters by SEQ_TYPE
	tmp		<- subset(rp, SC_WINDOW<=2)
	rpi		<- tmp[, list(	FIRSTPOSDATE=			FIRSTPOSDATE[1],
							SEQ_R_N= 				ifelse(any(!is.na(SEQID)), length(unique(SEQID)), 0L),
							SEQ_R_PANGEA= 			any(grepl('PANGEA', SEQ_TYPE)),
							PHP24_R_N= 				ifelse(any(!is.na(SEQIDb)), length(unique(SEQIDb)), 0L),
							PHP24_TR_N= 			ifelse(any(!is.na(SEQIDb_TR)), length(unique(SEQIDb_TR)), 0L),
							PHP24_MTR_N= 			ifelse(any(!is.na(SEQIDb_TR[SEX_TR=='M'])), length(unique(SEQIDb_TR[SEX_TR=='M'])), 0L),
							PHP24_MTR_MINPD=		ifelse(any(!is.na(PD[SEX_TR=='M'])), min(PD[SEX_TR=='M'], na.rm=TRUE), NA_real_),
							PHP24_MTR_PANGEA= 		length(which(grepl('PANGEA', SEQTYPE_TR[SEX_TR=='M']))),
							PHP24_MTR_PANGEA_MINPD=	ifelse(any(!is.na(PD[SEX_TR=='M' & SEQTYPE_TR=='PANGEA'])), min(PD[SEX_TR=='M' & SEQTYPE_TR=='PANGEA'], na.rm=TRUE), NA_real_)
							), by='RID']
	#	
	set(rpi, NULL, 'FIRSTPOSDATEc', rpi[, cut(FIRSTPOSDATE, breaks=c(1994,1999, 2011, 2016), labels=c('<1999','<2011','2011-2015'))])
	rpi[, REC_S:= as.integer(SEQ_R_N>0)]
	set(rpi, rpi[, which(PHP24_R_N>0)], 'REC_S', 2L)
	set(rpi, rpi[, which(SEQ_R_PANGEA & REC_S==1L)], 'REC_S', 3L)
	set(rpi, rpi[, which(SEQ_R_PANGEA & REC_S==2L)], 'REC_S', 4L)	
	set(rpi, NULL, 'REC_S', rpi[, factor(REC_S, levels=c(0L,1L,2L,3L,4L), labels=c('recipient has no sequence','recipient sequenced HISTORIC only','recipient in p24 tree HISTORIC only','recipient sequenced PANGEA','recipient in p24 tree PANGEA'))])
	#
	rpi[, TR_S:= as.integer(PHP24_R_N>0)]
	set(rpi, rpi[, which(PHP24_TR_N>0)], 'TR_S', 1L)	
	set(rpi, rpi[, which(PHP24_MTR_N>0)], 'TR_S', 2L)
	set(rpi, rpi[, which(PHP24_MTR_PANGEA>0)], 'TR_S', 3L)
	set(rpi, rpi[, which(!is.na(PHP24_MTR_MINPD) & !is.na(PHP24_MTR_PANGEA_MINPD) & PHP24_MTR_MINPD==PHP24_MTR_PANGEA_MINPD)], 'TR_S', 4L)
	set(rpi, NULL, 'TR_S', rpi[, factor(TR_S, levels=c(0L,1L,2L,3L,4L), labels=c('recipient not in p24 tree','recipient in p24 tree has only female pr tr','recipient in p24 tree has male pr tr in HISTORIC only','recipient in p24 tree with male pr tr in PANGEA that is not closest','recipient in p24 tree with male pr tr in PANGEA that is closest'))])
	
	 
	ggplot(rpi, aes(x=REC_S, fill=as.factor(FIRSTPOSDATEc))) + geom_bar() + coord_flip() +
			theme_bw() + theme(legend.position='bottom') +
			labs(y='count',x='',fill='date first positive', title='RCCS female seroconverter\n(seroconversion window <= 2yrs)\n')
	ggsave(file=file.path(wdir,paste('150825_FastTree_Recipients_SCw2_patristic',pd.max*100,'.pdf',sep='')), h=7, w=7)
	
	ggplot(rpi, aes(x=TR_S, fill=as.factor(FIRSTPOSDATEc))) + geom_bar() + coord_flip() +
			theme_bw() + theme(legend.position='bottom') + facet_grid(.~REC_S) +
			scale_y_continuous(minor_breaks=seq(0,1e3,10)) +
			labs(y='count',x='',fill='date first positive', title=paste('pairs with RCCS female seroconverter\n(seroconversion window <= 2yrs)\n(patristic distance <',100*pd.max,'%)\n'))
	ggsave(file=file.path(wdir,paste('150825_FastTree_Pairs_SCw2_patristic',pd.max*100,'.pdf',sep='')), h=7, w=14)
	
	
	rpi[, table(REC_S, TR_S)]	
	subset(rpi, REC_S=='recipient in p24 tree PANGEA' & TR_S=='recipient in p24 tree with male pr tr in PANGEA that is closest')
	#	34		at 4%
	#	24	 	at 2%
	subset(rpi, REC_S%in%c('recipient in p24 tree HISTORIC only','recipient in p24 tree PANGEA') &
				TR_S%in%c('recipient in p24 tree with male pr tr in PANGEA that is not closest','recipient in p24 tree with male pr tr in PANGEA that is closest','recipient in p24 tree has male pr tr in HISTORIC only'))
	#	191		at 4%
	#	121		at 2%


	#	check whereabouts of those not yet processed
	tmp		<- unique(subset(rt, SEQ_TYPE=='Sanger not started by Jul2016', select=c(RID, ROUND, SEX)))
	tmp		<- merge(tmp, subset(rd, select=c(RID, PID)), by='RID')
	write.table(tmp, file.path(wdir, 'RCCSparticipants_sequencingunclear.csv'))
	infile.sangerstats	<- "~/Dropbox (Infectious Disease)/PANGEA_data/2016-07-07_PANGEA_stats_by_sample.csv"
	tmp2	<- as.data.table(read.csv(infile.sangerstats, stringsAsFactors=FALSE))	
	setnames(tmp2, colnames(tmp2), gsub('\\.','_',toupper(colnames(tmp2))))
	setnames(tmp2, c('PROJECTID','STATUS'), c('PID','SANGER_STATUS'))
	merge(tmp, tmp2, by='PID')	# they are not there!
}
######################################################################################
#
######################################################################################
project.Rakai.aliRegion1.add.HXB2.RCCSmissing<- function()
{
	wdir	<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	#	load Susanna s codon alignment
	susa.f	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah3/Region1_codon aligned.fasta'
	susa.s	<- read.dna(susa.f, format='fasta')
	susa.d	<- data.table(	ID=gsub('\\*.*','',rownames(susa.s)),
							DATA= factor(grepl('^PG', rownames(susa.s)), levels=c(TRUE,FALSE), label=c('PNG','LNL')))
	infile.relabel	<- "~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional phylogenetic analyses/Region 1 gag analysis/SummaryofGagSequenceData.rda"				
	load(infile.relabel)
	s		<- as.data.table(summaryData)
	s[, TAXA:= gsub('a\\||b\\|','',seqid)]		#taxa names have been tinkered with. does not merge.
	s[, ID:=gsub('\\*.*','',TAXA)]	
	stopifnot( nrow(subset(s, grepl('\\|',TAXA)))==0 )
	stopifnot( nrow(subset(s, is.na(ID)))==0 )
	#
	#	see what s missing
	#
	ssu		<- merge(s, susa.d, by='ID', all=1)
	cat('\nseqs in codon alignment that are not in summaryData', nrow(subset(ssu, is.na(seqid))))
	#	0, yay
	cat('\nseqs in summaryData that are not in codon alignment', nrow(subset(ssu, is.na(DATA))))
	#	3094 
	#	subset(ssu, is.na(DATA))[, table(dataSource, useNA='if')]
	#	select missing RCCS historical and HXB2 
	tmp		<- subset(ssu, is.na(DATA) & dataSource=='RCCS historical' )
	tmp		<- rbind(tmp,subset(ssu, is.na(DATA) & dataSource=='reference' & grepl('HXB2',TAXA)))
	tmp[, relabel2:= gsub('^[^\\*]+\\*','',relabel)]
	# 
	#	get missing sequences
	#
	infile.gag	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree/gag.sqn.fasta'
	seq.m		<- read.dna(infile.gag, format='fasta')
	seq.m		<- seq.m[ tmp[, relabel2], ]
	file.seq.m	<- file.path('/Users/Oliver/duke/tmp', 'tmp_missingseq_150825.fasta')
	write.dna(seq.m, file=file.seq.m, format='fasta')
	#
	#	add missing sequences to codon alignment, potentially clipping missing seqs
	#
	outfile		<- file.path('/Users/Oliver/duke/tmp', '150825_Region1_codonaligned_p_missing_p_HXB2.fasta')
	cmd			<- paste('mafft --anysymbol',' --keeplength --mapout',' --add "',file.seq.m,'" --auto "',file.path('/Users/Oliver/duke/tmp',gsub(' ','_',basename(susa.f))),'" > ', outfile, '', sep='')
	system(cmd)
	#	inspecting mapout: HXB2 no internal sites removed, only end clipped 
	#
	#	move to '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments'
	#
	file.rename(outfile, file.path('~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments',basename(outfile)))
}
######################################################################################
#
######################################################################################
project.Rakai.aliRegion1.merge.with.PANGEA.160825<- function()
{	
	#	load Susanna s codon alignment and remove any PANGEA seqs
	susa.f	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151129.fasta'
	susa.f	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_missing_p_HXB2.fasta'
	susa.s	<- read.dna(susa.f, format='fasta')
	susa.d	<- data.table(TAXA=rownames(susa.s), DATA= factor(grepl('^PG', rownames(susa.s)), levels=c(TRUE,FALSE), label=c('PNG','LNL')))	
	in.s	<- susa.s[subset(susa.d, DATA!='PNG')[, TAXA],]	
	#	load PANGEA alignment	
	png.f	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/PANGEA_orig/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	load(png.f)	#sq, sqi, si
	sq		<- sq[grepl('^PG[0-9]+|HXB2',rownames(sq)),]
	#	required: HXB2 in alignment
	sqn		<- seq.align.based.on.common.reference(in.s, sq, return.common.sites=TRUE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_HXB2.fasta'
	write.dna(sqn, file=outfile, format='fasta', colsep='', nbcol=-1)
	
	#
	#	subset to Rakai and relabel
	#
	
	#	load summaryData
	infile.relabel	<- "~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional phylogenetic analyses/Region 1 gag analysis/SummaryofGagSequenceData.rda"				
	load(infile.relabel)
	s		<- as.data.table(summaryData)
	s[, TAXA:= gsub('a\\||b\\|','',seqid)]		#taxa names have been tinkered with. does not merge.
	s[, ID:=gsub('\\*.*','',TAXA)]	
	stopifnot( nrow(subset(s, grepl('\\|',TAXA)))==0 )
	stopifnot( nrow(subset(s, is.na(ID)))==0 )	
	#	get info on new alignment
	sqni		<- data.table(SEQID=rownames(sqn), IDX=seq_len(nrow(sqn)))		
	sqni[, SRC:= factor(grepl('^PG',SEQID),levels=c(TRUE,FALSE),labels=c('PNG','LANL'))]
	set(sqni, NULL, "TAXA_ID", sqni[, gsub('RCCS\\*|Ref\\*', '', SEQID)])
	
	sqni[, ID:= gsub('\\*.*','',TAXA_ID)]	
	sqni		<- merge(sqni, s, by='ID', all.x=1)
	#
	#	collect info on taxa missing in summaryData
	#
	stopifnot( nrow(subset(sqni, is.na(IDX)))==0 )
	not.in.summaryData	<- subset(sqni, is.na(relabel) & SRC!='PNG' & !grepl('HXB2', SEQID))
	#	deselect PANGEA ZA and BW
	sqni		<- subset(sqni, !(is.na(relabel) & SRC=='PNG' & !grepl('UG', SEQID)))
	#	
	set(sqni, sqni[, which(is.na(relabel) & SRC=='PNG' & grepl('UG', SEQID))], 'dataSource', 'UG-MRC')
	set(sqni, sqni[, which(grepl('HXB2', SEQID))], 'dataSource', 'reference')
	set(sqni, sqni[, which(is.na(relabel) & SRC!='PNG' & grepl('RCCS', SEQID))], 'dataSource', 'RCCS historical')
	set(sqni, sqni[, which(is.na(relabel) & SRC!='PNG' & !grepl('RCCS', SEQID))], 'dataSource', 'reference')
	#
	setkey(sqni, IDX)
	set(sqni, NULL, c('IDX','TAXA_ID','ID','SRC','TAXA'), NULL)
	set(sqni, NULL, 'index', 1:nrow(sqni))
	sq			<- sqn[sqni[, SEQID],]
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1UG_codonaligned_p_PANGEA151113_p_HXB2.fasta'
	write.dna(sqn, file=outfile, format='fasta', colsep='', nbcol=-1)
	summaryData	<- copy(sqni)
	save(sq, summaryData, file=gsub('\\.fasta','.rda',outfile))
	
	
	#
	#	old code perhaps useful at some point
	#
	
	tmp			<- sqni[, which(SRC=='LANL' & IDX<300)]
	set(sqni, tmp, 'SRC', 'COMPENDIUM')	
	tmp			<- as.character(sqn)
	tmp			<- data.table(	TAXA=rownames(tmp), 
								FIRST= apply( tmp, 1, function(x) which(!x%in%c('-','?'))[1] ),
								LAST= ncol(tmp)-apply( tmp, 1, function(x) which(!rev(x)%in%c('-','?'))[1] ) + 1L,						
								COV=apply(tmp, 1, function(x) length(which(!x%in%c('-','?')))	))
	sqni		<- merge(sqni,tmp,by='TAXA')
	tmp			<- subset(sqni, SRC=='LANL')[, max(LAST)]
	tmp			<- as.character(sqn[,seq_len(tmp)])
	tmp			<- data.table(	TAXA=rownames(tmp), 
								COV_REGION=apply(tmp, 1, function(x) length(which(!x%in%c('-','?')))	))
	sqni		<- merge(sqni,tmp,by='TAXA')
	sqni[, SUBT:=NA_character_]
	tmp			<- sqni[, which(SRC=='COMPENDIUM')]
	set(sqni, tmp, 'SUBT', sqni[tmp, sapply(strsplit(TAXA,'.',fixed=1),'[[',1)])
	tmp			<- sqni[, which(SRC=='LANL')]
	set(sqni, tmp, 'SUBT', sqni[tmp, sapply(strsplit(TAXA,'*',fixed=1),'[[',3)])
	tmp			<- sqni[, which(SUBT%in%c('-','U'))]
	set(sqni, tmp, 'SUBT', NA_character_)
	sqni[, ACCN:=NA_character_]
	tmp			<- sqni[, which(SRC=='COMPENDIUM')]
	set(sqni, tmp, 'ACCN', sqni[tmp, sapply(strsplit(TAXA,'.',fixed=1),function(x) rev(x)[1])])
	tmp			<- sqni[, which(SRC=='LANL')]
	set(sqni, tmp, 'ACCN', sqni[tmp, sapply(strsplit(TAXA,'*',fixed=1),'[[',1)])
	write.dna( sqn, format='fasta', file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201.fasta', colsep='', nbcol=-1)
	save(sqni, sqn, file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201.rda')
	
	sqn				<- read.dna(file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201_rm99gps.fasta',format='fasta')
	rownames(sqn)	<- gsub(' (stripped)','',rownames(sqn),fixed=1)
	#	get info
	sqni		<- data.table(TAXA=rownames(sqn), IDX=seq_len(nrow(sqn)))		
	sqni[, SRC:= factor(grepl('^PG',TAXA),levels=c(TRUE,FALSE),labels=c('PNG','LANL'))]
	tmp			<- sqni[, which(SRC=='LANL' & IDX<300)]
	set(sqni, tmp, 'SRC', 'COMPENDIUM')	
	tmp			<- as.character(sqn)
	tmp			<- data.table(	TAXA=rownames(tmp), 
			FIRST= apply( tmp, 1, function(x) which(!x%in%c('-','?'))[1] ),
			LAST= ncol(tmp)-apply( tmp, 1, function(x) which(!rev(x)%in%c('-','?'))[1] ) + 1L,						
			COV=apply(tmp, 1, function(x) length(which(!x%in%c('-','?')))	))
	sqni		<- merge(sqni,tmp,by='TAXA')
	tmp			<- subset(sqni, SRC=='LANL')[, max(LAST)]
	tmp			<- as.character(sqn[,seq_len(tmp)])
	tmp			<- data.table(	TAXA=rownames(tmp), 
			COV_REGION=apply(tmp, 1, function(x) length(which(!x%in%c('-','?')))	))
	sqni		<- merge(sqni,tmp,by='TAXA')
	sqni[, SUBT:=NA_character_]
	tmp			<- sqni[, which(SRC=='COMPENDIUM')]
	set(sqni, tmp, 'SUBT', sqni[tmp, sapply(strsplit(TAXA,'.',fixed=1),'[[',1)])
	tmp			<- sqni[, which(SRC=='LANL')]
	set(sqni, tmp, 'SUBT', sqni[tmp, sapply(strsplit(TAXA,'*',fixed=1),'[[',3)])
	tmp			<- sqni[, which(SUBT%in%c('-','U'))]
	set(sqni, tmp, 'SUBT', NA_character_)
	sqni[, ACCN:=NA_character_]
	tmp			<- sqni[, which(SRC=='COMPENDIUM')]
	set(sqni, tmp, 'ACCN', sqni[tmp, sapply(strsplit(TAXA,'.',fixed=1),function(x) rev(x)[1])])
	tmp			<- sqni[, which(SRC=='LANL')]
	set(sqni, tmp, 'ACCN', sqni[tmp, sapply(strsplit(TAXA,'*',fixed=1),'[[',1)])
	save(sqni, sqn, file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201_rm99gps.rda')	
	write.dna( sqn, format='fasta', file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201_rm99gps.fasta', colsep='', nbcol=-1)
}
######################################################################################
project.Rakai.checkforARVs.150910<- function()
{
	require(data.table)
	#	Kate, I generally use data.table rather than data.frame since it s faster and more versatile. 
	#	warning: it s a bit of a learning curve though to use it
	
	require(big.phylo)
	#	library(devtools)
	#	install_github("olli0601/big.phylo")
	
	f.arv	<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_150909.rda'
	f.rccsid<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/150625_PangeaBoxes.csv'
	f.sid	<- '~/Dropbox (Infectious Disease)/PANGEA_data/SangerUpdates/2015-07-20_PANGEA_3.csv'
	f.seq	<- '~/Dropbox (Infectious Disease)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908.fasta'
	#
	#	get IDs into OK format
	#
	load(f.arv)	#	loads arvdat
	arvdat	<- as.data.table(arvdat)
	d.rccsid<- as.data.table(read.csv(f.rccsid, stringsAsFactors=0))
	d.rccsid<- subset(d.rccsid, select=c(Pangea.ID, Rakai.Study.ID))
	setnames(d.rccsid, c('Rakai.Study.ID','Pangea.ID'), c('RCCS_ID','PNG_ID'))
	set(d.rccsid, NULL, 'RCCS_ID', d.rccsid[, substr(RCCS_ID, 1, nchar(RCCS_ID)-3)])
	set(d.rccsid, NULL, 'PNG_ID', d.rccsid[, gsub('-Sy*','',PNG_ID)])
	setnames(arvdat, 'RCCS_studyid', 'RCCS_ID')
	d.sid	<- as.data.table(read.csv(f.sid, stringsAsFactors=0, header=0))
	setnames(d.sid, c('V1','V2','V3'), c('PNG_ID_FULL', 'SNG_ID', 'x'))
	d.sid[, PNG_ID:= gsub('-S[0-9]+','',PNG_ID_FULL)]
	set(d.sid, NULL, 'SNG_ID', d.sid[, gsub('#','_',SNG_ID)])
	d.sid[, x:=NULL]
	#
	#	merge
	#	
	subset(arvdat, any.pangea==1)	
	#	126		
	arvdat	<- merge( arvdat, d.rccsid, by='RCCS_ID' )
	#	85
	arvdat	<- merge( arvdat, d.sid, by='PNG_ID')
	#	76
	
	#
	#	check for DRMs
	#
	seq				<- read.dna(file=f.seq, format='fasta')	
	seq				<- seq[c("B.FR.83.HXB2_LAI_IIIB_BRU.K03455",arvdat[, PNG_ID_FULL]),]
	rownames(seq)[1]	<- 'HXB2'
	outfile			<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_PotentialDRMs_OR_150910.R'
	tmp				<- big.phylo:::seq.rm.drugresistance(seq, outfile=outfile)
	nodr.info		<- tmp$nodr.info
	nodr.seq		<- tmp$nodr.seq
}
######################################################################################
project.Rakai.checkforARVs.150911<- function()
{
	require(data.table)
	#	Kate, I generally use data.table rather than data.frame since it s faster and more versatile. 
	#	warning: it s a bit of a learning curve though to use it
	
	require(big.phylo)
	#	library(devtools)
	#	install_github("olli0601/big.phylo")
	
	f.arv			<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_150911.rda'
	f.rccsid		<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/150625_PangeaBoxes.csv'
	f.sid			<- '~/Dropbox (Infectious Disease)/PANGEA_data/SangerUpdates/2015-07-20_PANGEA_3.csv'
	f.seq			<- '~/Dropbox (Infectious Disease)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908.fasta'
	#
	load(f.arv)	#	loads rakaidat
	rakaidat		<- as.data.table(rakaidat)
	setnames(rakaidat, 'ProjectID', 'PNG_ID')
	#	482
	seq				<- read.dna(file=f.seq, format='fasta')	
	rownames(seq)	<- gsub('-S[0-9]+','',rownames(seq))
	d.seq			<- merge( data.table(PNG_ID=rownames(seq)), rakaidat, by='PNG_ID' )
	#	459 with PANGEA sequence
	
	
	#
	#	check for DRMs
	#
	seq				<- seq[c("B.FR.83.HXB2_LAI_IIIB_BRU.K03455",d.seq[, PNG_ID]),]		
	rownames(seq)[1]	<- 'HXB2'
	outfile			<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_PotentialDRMs_OR_150911.R'
	tmp				<- big.phylo:::seq.rm.drugresistance(seq, outfile=outfile)
	nodr.info		<- tmp$nodr.info
	nodr.seq		<- tmp$nodr.seq
}