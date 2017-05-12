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

RakaiCirc.recipient.female.get.info<- function(wdir=NA)
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
	rrec	<- subset(rrec, !is.na(LASTNEGVIS), select=c(RID, SEX, VISIT, DATE, REGION, COMM_NUM, HH_NUM, BIRTHDATE, LASTNEGVIS, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE))
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

RakaiCirc.epi.get.info.170120<- function()
{
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
	#	define sexual relationships
	set(rh, rh[, which(is.na(RLTN1))],'RLTN1',99L)
	set(rh, rh[, which(is.na(RLTN2))],'RLTN2',99L)
	set(rh, rh[, which(is.na(RLTN3))],'RLTN3',99L)
	set(rh, rh[, which(is.na(RLTN4))],'RLTN4',99L)
	setnames(rh, c('RLTN1','RLTN2','RLTN3','RLTN4'), c('RLTN1_','RLTN2_','RLTN3_','RLTN4_'))
	warning('undocumented RLTN codes 15 16 88 -> set to Unknown')
	tmp		<- as.data.table(matrix(data=c(	'1','Current wife (at the time)',
							'2','Current consensual partner (at the time)',
							'3','Former wife/consensual partner',
							'4','Girlfriend',
							'5','Occasional or casual friend',
							'6','Visitor (incl. wedding/funeral)',
							'7','Stranger',
							'8','Workmate',
							'9','Boss/work supervisor',
							'10','Employee',
							'11','Fellow student',
							'12','Sugar mummy',
							'13','Relative other than spouse',
							'14','Other non relative',
							'15','Unknown',
							'16','Unknown',
							'88','Unknown',
							'0','Unknown',
							'97','Unknown',
							'98','Unknown',
							'99','Unknown'
					), ncol=2, byrow=TRUE))
	setnames(tmp, colnames(tmp), c('RLTN1_','RLTN1'))
	set(tmp, NULL,'RLTN1_',tmp[, as.integer(RLTN1_)])
	rh		<- merge(rh,tmp,by='RLTN1_')
	setnames(tmp, c('RLTN1_','RLTN1'), c('RLTN2_','RLTN2'))	
	rh		<- merge(rh,tmp,by='RLTN2_')
	setnames(tmp, c('RLTN2_','RLTN2'), c('RLTN3_','RLTN3'))	
	rh		<- merge(rh,tmp,by='RLTN3_')
	setnames(tmp, c('RLTN3_','RLTN3'), c('RLTN4_','RLTN4'))	
	rh		<- merge(rh,tmp,by='RLTN4_')
	set(rh, NULL, c('RLTN1_','RLTN2_','RLTN3_','RLTN4_'), NULL)	
	tmp		<- rh[, which(RLTN1=='Unknown' & RLTN2!='Unknown')]
	set(rh, tmp, 'RLTN1', rh[tmp, RLTN2])
	set(rh, tmp, 'RLTN2', 'Unknown')
	tmp		<- rh[, which(RLTN2=='Unknown' & RLTN3!='Unknown')]
	set(rh, tmp, 'RLTN2', rh[tmp, RLTN3])
	set(rh, tmp, 'RLTN3', 'Unknown')
	tmp		<- rh[, which(RLTN3=='Unknown' & RLTN4!='Unknown')]
	set(rh, tmp, 'RLTN3', rh[tmp, RLTN4])
	set(rh, tmp, 'RLTN4', 'Unknown')
	#	define number named sexual relations last year
	tmp		<- melt(rh, id.vars=c('RID','VISIT'), measure.vars=c('RLTN1','RLTN2','RLTN3','RLTN4'))
	tmp		<- tmp[, list(RLTN_NAMED= length(which(value!='Unknown'))), by=c('RID','VISIT')]
	rh		<- merge(rh, tmp, by=c('RID','VISIT'))
	#	define occuption
	set(rh, rh[, which(is.na(OCCUP1))],'OCCUP1',99L)
	set(rh, rh[, which(is.na(OCCUP2))],'OCCUP2',99L)
	setnames(rh, c('OCCUP1','OCCUP2'), c('OCCUP1_','OCCUP2_'))
	tmp		<- as.data.table(matrix(data=c(	'1','Agriculture for home use/barter',
							'2','Agriculture for selling',
							'3','Housework in your own home',
							'4','Housekeeper',
							'5','Home brewing',
							'6','Government/clerical/teaching',
							'7','Fishing',
							'8','Student',
							'9','Military/police',
							'10','Shopkeeper',
							'11','Trading/vending',
							'12','Bar worker or owner',
							'13','Trucker',
							'14','Unemployed',
							'15','Other',
							'88','No additional occupation',
							'16','Medical worker',
							'17','Casual laborer',
							'18','Waitress/Waiter/restaurant owner',
							'19','Hair dresser/Salon owner',
							'20','Construction (brick maker, builder, porter, painter, roofing)',
							'21','Mechanic (automobiles, bicycles, electronics)',
							'22','Boda Boda',
							'23','Client/Sex worker',
							'0','Unknown',
							'98','Unknown',
							'99','Unknown'), ncol=2, byrow=TRUE))
	setnames(tmp, colnames(tmp), c('OCCUP1_','OCCUP1'))
	set(tmp, NULL,'OCCUP1_',tmp[, as.integer(OCCUP1_)])
	rh		<- merge(rh,tmp,by='OCCUP1_')
	setnames(tmp, c('OCCUP1_','OCCUP1'), c('OCCUP2_','OCCUP2') )
	rh		<- merge(rh,tmp,by='OCCUP2_')
	set(rh, NULL, c('OCCUP1_','OCCUP2_'), NULL)
	tmp		<- rh[, which(OCCUP1=='Unknown' & OCCUP2!='Unknown')]
	set(rh, tmp, 'OCCUP1', rh[tmp, OCCUP2])
	set(rh, tmp, 'OCCUP2', 'Unknown')
	#	define SEXYEAR
	set(rh, rh[, which(is.na(SEXYEAR))],'SEXYEAR', 99)
	set(rh, NULL, 'SEXYEAR', rh[, gsub('_[0-9]$','',as.character(factor(SEXYEAR, levels=c(0,1,2,8,99), labels=c('Unknown_1','Y','N','Unknown_2','Unknown_3'))))])
	stopifnot( !nrow(subset(rh, is.na(SEXYEAR))) )	
	#	define SEXP1YR 
	setnames(rh, c('SEXP1YR'), c('SEXP1YR_'))	
	rh[, SEXP1YR:= as.character(SEXP1YR_)]
	set(rh, rh[, which(SEXP1YR_==92)],'SEXP1YR','<3')
	set(rh, rh[, which(SEXP1YR_==93)],'SEXP1YR','3+')
	set(rh, rh[, which(SEXP1YR_%in%c(97,98,99) | is.na(SEXP1YR_))],'SEXP1YR','Unknown')
	set(rh, NULL, 'SEXP1YR_', NULL)
	#	revisit SEXP1YR based on relationships
	set(rh, rh[, which(SEXP1YR%in%c('Unknown') & RLTN_NAMED==1)], 'SEXP1YR', '1')
	set(rh, rh[, which(SEXP1YR%in%c('Unknown') & RLTN_NAMED==2)], 'SEXP1YR', '2')
	set(rh, rh[, which(SEXP1YR%in%c('Unknown') & RLTN_NAMED==3)], 'SEXP1YR', '3')
	set(rh, rh[, which(SEXP1YR%in%c('Unknown') & RLTN_NAMED==4)], 'SEXP1YR', '3+')
	tmp		<- rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED>0)]
	if( nrow(rh[tmp,]) )
		warning("found SEXP1YR%in%c('0') & RLTN_NAMED>0  --> set to had sex last year, n=", length(tmp))
	set(rh, rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED==1)], 'SEXP1YR', '1')
	set(rh, rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED==2)], 'SEXP1YR', '2')
	set(rh, rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED==3)], 'SEXP1YR', '3')
	set(rh, rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED==4)], 'SEXP1YR', '3+')
	#	revisit SEXYEAR based on SEXP1YR
	set(rh, rh[, which(SEXYEAR%in%c('Unknown') & !SEXP1YR%in%c('0','Unknown'))], 'SEXYEAR', 'Y')
	tmp		<- rh[, which(SEXYEAR%in%c('N') & !SEXP1YR%in%c('0','Unknown'))]
	if( nrow(rh[tmp,]) )
		warning("found SEXYEAR%in%c('N') & !SEXP1YR%in%c('0','Unknown') --> set to had sex last year, n=", length(tmp))	
	set(rh, tmp, 'SEXYEAR', 'Y')	
	#	define SEXP1OUT 
	setnames(rh, c('SEXP1OUT'), c('SEXP1OUT_'))	
	rh[, SEXP1OUT:= as.character(SEXP1OUT_)]
	set(rh, rh[, which(SEXP1OUT_==92)],'SEXP1OUT','<3')
	set(rh, rh[, which(SEXP1OUT_==93)],'SEXP1OUT','3+')
	set(rh, rh[, which(SEXP1OUT_%in%c(97,98,99) | is.na(SEXP1OUT_))],'SEXP1OUT','Unknown')
	set(rh, NULL, 'SEXP1OUT_', NULL)
	#	revisit SEXYEAR based on SEXP1OUT
	set(rh, rh[, which(SEXYEAR%in%c('Unknown') & !SEXP1OUT%in%c('0','Unknown'))], 'SEXYEAR', 'Y')
	tmp		<- rh[, which(SEXYEAR%in%c('N') & !SEXP1OUT%in%c('0','Unknown'))]
	if( nrow(rh[tmp,]) )
		warning("found SEXYEAR%in%c('N') & !SEXP1OUT%in%c('0','Unknown') --> set to had sex last year, n=", length(tmp))
	set(rh, tmp, 'SEXYEAR', 'Y')
	#	revisit SEXP1YR based on SEXP1OUT
	set(rh, rh[, which(SEXP1OUT=='3+' & SEXP1YR%in%c('0','1','2','<3','Unknown'))], 'SEXP1YR', '3+')
	set(rh, rh[, which(SEXP1OUT=='<3' & SEXP1YR%in%c('0','Unknown'))], 'SEXP1YR', '<3')
	rh[, DUMMY:= seq_len(nrow(rh))]
	warning("set(rh, rh[, which(RID=='G013746' & VISIT==14)], 'SEXP1OUT', 1) --> think this is typo")
	set(rh, rh[, which(RID=='G013746' & VISIT==14)], 'SEXP1OUT', '1')	
	tmp		<- rh[, which(	!SEXP1OUT%in%c('3+','<3','Unknown') & !SEXP1YR%in%c('3+','<3','Unknown'))]  
	tmp		<- subset(rh[tmp, ], as.numeric(SEXP1OUT)>as.numeric(SEXP1YR))[, DUMMY]
	warning("as.numeric(SEXP1OUT)>as.numeric(SEXP1YR), set to SEXP1OUT, n=", length(tmp))
	set(rh, tmp, 'SEXP1YR', rh[tmp, SEXP1OUT])
	#	revisit SEXP1YR based on SEXYEAR
	set(rh, rh[, which(SEXYEAR=='Unknown' & SEXP1YR=='0')], 'SEXP1YR', 'Unknown')
	set(rh, rh[, which(SEXYEAR=='Y' & SEXP1YR=='0')], 'SEXP1YR', 'Unknown')
	#	define SEXPEVER
	setnames(rh, c('SEXPEVER'), c('SEXPEVER_'))	
	rh[, SEXPEVER:= as.character(SEXPEVER_)]
	set(rh, rh[, which(SEXPEVER_==92)],'SEXPEVER','<3')
	set(rh, rh[, which(SEXPEVER_==93)],'SEXPEVER','3+')
	set(rh, rh[, which(SEXPEVER_%in%c(97,98,99) | is.na(SEXPEVER_))],'SEXPEVER','Unknown')
	set(rh, NULL, 'SEXPEVER_', NULL)
	stopifnot( !nrow(subset(rh, is.na(SEXPEVER))) )
	#	revisit SEXPEVER based on SEXP1YR
	set(rh, rh[, which(SEXP1YR=='3+' & SEXPEVER%in%c('0','1','2','<3','Unknown'))], 'SEXPEVER', '3+')
	set(rh, rh[, which(SEXP1YR=='<3' & SEXPEVER%in%c('0','Unknown'))], 'SEXPEVER', '<3')	
	tmp		<- rh[, which(	!SEXP1YR%in%c('3+','<3','Unknown') & !SEXPEVER%in%c('3+','<3','Unknown'))]  
	tmp		<- subset(rh[tmp, ], as.numeric(SEXP1YR)>as.numeric(SEXPEVER))[, DUMMY]
	warning("as.numeric(SEXP1YR)>as.numeric(SEXPEVER), set to SEXP1YR, n=", length(tmp))
	set(rh, tmp, 'SEXPEVER', rh[tmp, SEXP1YR])
	#	define EVERSEX 
	set(rh, rh[, which(is.na(EVERSEX))],'EVERSEX', 99)
	set(rh, NULL, 'EVERSEX', rh[, gsub('_[0-9]$','',as.character(factor(EVERSEX, levels=c(0,1,2,3,8,99), labels=c('Unknown_1','Y','N','Unknown_2','Unknown_3','Unknown_4'))))])
	stopifnot( !nrow(subset(rh, is.na(EVERSEX))) )	
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & SEXYEAR=='Y')], 'EVERSEX', 'Y')
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & CURRMARR==1)], 'EVERSEX', 'Y')
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & EVERMARR==1)], 'EVERSEX', 'Y')
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & RLTN_NAMED>0)], 'EVERSEX', 'Y')
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & !is.na(AGE1STSEX)==1)], 'EVERSEX', 'Y')
	#	revisit EVERSEX based on SEXPEVER
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & !SEXPEVER%in%c('0','Unknown'))], 'EVERSEX', 'Y')
	tmp		<- rh[, which(EVERSEX%in%c('N') & !SEXPEVER%in%c('0','Unknown'))]
	if( nrow(rh[tmp,]) )
		warning("found EVERSEX%in%c('N') & !SEXPEVER%in%c('0','Unknown') --> set to had sex ever, n=", length(tmp))
	set(rh, tmp, 'EVERSEX', 'Y')
	#	revisit SEXPEVER based on EVERSEX
	set(rh, rh[, which(EVERSEX=='Unknown' & SEXPEVER=='0')], 'SEXP1YR', 'Unknown')
	set(rh, rh[, which(EVERSEX=='Y' & SEXPEVER=='0')], 'SEXP1YR', 'Unknown')	
	#	check SEXPEVER	
	tmp		<- subset(rh, SEXPEVER=='0' & HIVPREV==1)
	if( nrow(tmp) )
		warning("found SEXPEVER=='0' & HIVPREV==1 --> report only, n=", nrow(tmp))	
	#	define sexual contacts
	rh[, SEXC:= as.character(factor(CURRMARR,levels=c(1,2,8),labels=c('currently married','CURRMARR_No','CURRMARR_NA')))]	
	set(rh, rh[, which(grepl('currently married',SEXC) & !POLYMAR%in%c(0,1,98,99))], 'SEXC', 'currently married >1 wives')
	set(rh, rh[, which(grepl('CURRMARR',SEXC) & EVERMARR==1)], 'SEXC', 'previously married')
	set(rh, rh[, which(grepl('CURRMARR',SEXC) & SEXYEAR=='Y')], 'SEXC', 'sex contact last year')
	set(rh, rh[, which(grepl('CURRMARR',SEXC) & EVERSEX=='Y')], 'SEXC', 'sex contact before last year')
	set(rh, rh[, which(grepl('CURRMARR',SEXC) & EVERSEX=='N')], 'SEXC', 'no sex contact ever')
	set(rh, rh[, which(grepl('CURRMARR',SEXC))], 'SEXC', 'Unknown')	
	tmp		<- rh[, which(grepl(	'currently married', SEXC) & 
							(	!grepl('Current wife|Current consensual partner|Unknown', RLTN1) | 
								!grepl('Current wife|Current consensual partner|Unknown', RLTN2) | 
								!grepl('Current wife|Current consensual partner|Unknown', RLTN3) | 
								!grepl('Current wife|Current consensual partner|Unknown', RLTN4))	)]
	set(rh, tmp, 'SEXC', 'currently married and other partner')			
	#	check circumcision
	tmp		<- rh[, which(CIRC=='Y' & SEX=='F')]
	cat('\nWarning: found female circumcised --> set to NA' ,rh[tmp, paste(RID, collapse=' ')])	
	set(rh, tmp, 'CIRC', NA_integer_)
	#	set to NULL
	set(rh, NULL, c('VDEX','AGEYRS','ORALHC','INJHC','EVERMARR','MARORDER','CURRMARR','POLYMAR','AGE1STSEX','WHOFSTSEX','CIRCUM','CIRCY1','CIRCY2','DUMMY','RLTN1','RLTN2','RLTN3','RLTN4','RLTN_NAMED'), NULL)
	list(rd=rd, rh=rh)
}


RakaiCirc.epi.get.info.170208<- function()
{
	hivc.db.Date2numeric<- function( x )
	{
		if(!class(x)%in%c('Date','character'))	return( x )
		x	<- as.POSIXlt(x)
		tmp	<- x$year + 1900
		x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
		x	
	}
	#
	infile				<- "~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/alldat_r15tor17.rda"
	load(infile)
	ra		<- as.data.table(alldat)
	ra		<- subset(ra, select=c(RCCS_studyid,REGION,COMM_NUM,HH_NUM,SEX,AGEYRS,visdate,visit,lastNegDate,hiv,firstPosDate,eversex, evermarr, currmarr, polymar, sexpever, sexp1yr, sexp1out, sexgift, sexyear, religion, educate, educyrs, edcat, occ, occ2))
	#	a bit of clean up 	
	setnames(ra, colnames(ra), gsub('\\.','_',toupper(colnames(ra))))	
	set(ra, NULL, 'VISDATE', ra[, as.Date(VISDATE, format='%d/%m/%Y')])	
	for(x in colnames(ra))
		if(class(ra[[x]])=='Date')
			set(ra, NULL, x, hivc.db.Date2numeric(ra[[x]]))
	#	
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
	rn		<- as.data.table(neuroData)
	setnames(rn, colnames(rn), gsub('\\.','_',toupper(colnames(rn))))
	for(x in colnames(rn))
		if(class(rn[[x]])=='Date')
			set(rn, NULL, x, hivc.db.Date2numeric(rn[[x]]))
	#rd[, table(VISIT)]
	#	make shorter
	setnames(ra, 'RCCS_STUDYID', 'RID')
	setnames(rd, 'RCCS_STUDYID', 'RID')
	setnames(rd, 'PANGEA_ID', 'PID')
	setnames(rd, 'STUDYID', 'SID')
	setnames(rh, 'RCCS_STUDYID', 'RID')	
	setnames(rn, 'STUDYID', 'RID')
	setnames(rn, 'PANGEA_ID', 'PID')
	#	get neuro into format of rd
	setnames(rn, c('SAMPLEDATE','GENDER','CD4','VIRALLOAD'), c('DATE','SEX','RECENTCD4','RECENTVL'))
	set(rn, NULL, c('RECENTVLDATE','RECENTCD4DATE'), rn[, DATE])
	set(rn, NULL, c('TIMESINCEVL','TIMESINCECD4'), 0)
	#	data checks
	setkey(rh, VISIT, RID)
	stopifnot(nrow(rh)==nrow(unique(rh)))
	#	define circumcision	
	set(rh, rh[, which(!CIRCUM%in%c(1,2))], 'CIRCUM', NA_integer_)
	set(rh, NULL, 'CIRCUM', rh[, factor(CIRCUM, levels=c(1,2), labels=c('Y','N'))])
	#	define sexual relationships
	set(rh, rh[, which(is.na(RLTN1))],'RLTN1',99L)
	set(rh, rh[, which(is.na(RLTN2))],'RLTN2',99L)
	set(rh, rh[, which(is.na(RLTN3))],'RLTN3',99L)
	set(rh, rh[, which(is.na(RLTN4))],'RLTN4',99L)
	setnames(rh, c('RLTN1','RLTN2','RLTN3','RLTN4'), c('RLTN1_','RLTN2_','RLTN3_','RLTN4_'))
	warning('undocumented RLTN codes 15 16 88 -> set to Unknown')
	tmp		<- as.data.table(matrix(data=c(	'1','Current wife (at the time)',
											'2','Current consensual partner (at the time)',
											'3','Former wife/consensual partner',
											'4','Girlfriend',
											'5','Occasional or casual friend',
											'6','Visitor (incl. wedding/funeral)',
											'7','Stranger',
											'8','Workmate',
											'9','Boss/work supervisor',
											'10','Employee',
											'11','Fellow student',
											'12','Sugar mummy',
											'13','Relative other than spouse',
											'14','Other non relative',
											'15','Unknown',
											'16','Unknown',
											'88','Unknown',
											'0','Unknown',
											'97','Unknown',
											'98','Unknown',
											'99','Unknown'
									), ncol=2, byrow=TRUE))
	setnames(tmp, colnames(tmp), c('RLTN1_','RLTN1'))
	set(tmp, NULL,'RLTN1_',tmp[, as.integer(RLTN1_)])
	rh		<- merge(rh,tmp,by='RLTN1_')
	setnames(tmp, c('RLTN1_','RLTN1'), c('RLTN2_','RLTN2'))	
	rh		<- merge(rh,tmp,by='RLTN2_')
	setnames(tmp, c('RLTN2_','RLTN2'), c('RLTN3_','RLTN3'))	
	rh		<- merge(rh,tmp,by='RLTN3_')
	setnames(tmp, c('RLTN3_','RLTN3'), c('RLTN4_','RLTN4'))	
	rh		<- merge(rh,tmp,by='RLTN4_')
	set(rh, NULL, c('RLTN1_','RLTN2_','RLTN3_','RLTN4_'), NULL)	
	tmp		<- rh[, which(RLTN1=='Unknown' & RLTN2!='Unknown')]
	set(rh, tmp, 'RLTN1', rh[tmp, RLTN2])
	set(rh, tmp, 'RLTN2', 'Unknown')
	tmp		<- rh[, which(RLTN2=='Unknown' & RLTN3!='Unknown')]
	set(rh, tmp, 'RLTN2', rh[tmp, RLTN3])
	set(rh, tmp, 'RLTN3', 'Unknown')
	tmp		<- rh[, which(RLTN3=='Unknown' & RLTN4!='Unknown')]
	set(rh, tmp, 'RLTN3', rh[tmp, RLTN4])
	set(rh, tmp, 'RLTN4', 'Unknown')
	#	define number named sexual relations last year
	tmp		<- melt(rh, id.vars=c('RID','VISIT'), measure.vars=c('RLTN1','RLTN2','RLTN3','RLTN4'))
	tmp		<- tmp[, list(RLTN_NAMED= length(which(value!='Unknown'))), by=c('RID','VISIT')]
	rh		<- merge(rh, tmp, by=c('RID','VISIT'))
	#	define occupation
	#	set back OCCUP codes according to patient sheet
	setnames(rh, c('OCCUP1','OCCUP2','OCCUP3','OCCUP4'), c('OCCUP11','OCCUP12','OCCUP13','OCCUP14'))
	setnames(rh, c('OCC','OCC2'), c('OCCUP1','OCCUP2'))	
	set(rh, rh[, which(is.na(OCCUP1))],'OCCUP1',99L)
	set(rh, rh[, which(is.na(OCCUP2))],'OCCUP2',99L)
	set(rh, NULL, 'OCCUP1', rh[, as.integer(OCCUP1)])	
	setnames(rh, c('OCCUP1','OCCUP2'), c('OCCUP1_','OCCUP2_'))
	#	handle ra
	setnames(ra, c('OCC','OCC2'), c('OCCUP1','OCCUP2'))	
	set(ra, ra[, which(is.na(OCCUP1))],'OCCUP1',99L)
	set(ra, ra[, which(is.na(OCCUP2))],'OCCUP2',99L)
	set(ra, NULL, 'OCCUP1', ra[, as.integer(OCCUP1)])	
	setnames(ra, c('OCCUP1','OCCUP2'), c('OCCUP1_','OCCUP2_'))	
	tmp		<- as.data.table(matrix(data=c(	'1','Agriculture for home use/barter',
											'2','Agriculture for selling',
											'3','Housework in your own home',
											'4','Housekeeper',
											'5','Home brewing',
											'6','Government/clerical/teaching',
											'7','Fishing',
											'8','Student',
											'9','Military/police',
											'10','Shopkeeper',
											'11','Trading/vending',
											'12','Bar worker or owner',
											'13','Trucker',
											'14','Unemployed',
											'15','Other',
											'88','No additional occupation',
											'16','Medical worker',
											'17','Casual laborer',
											'18','Waitress/Waiter/restaurant owner',
											'19','Hair dresser/Salon owner',
											'20','Construction (brick maker, builder, porter, painter, roofing)',
											'21','Mechanic (automobiles, bicycles, electronics)',
											'22','Boda Boda',
											'23','Client/Sex worker',
											'0','Unknown',
											'98','Unknown',
											'99','Unknown'), ncol=2, byrow=TRUE))
	setnames(tmp, colnames(tmp), c('OCCUP1_','OCCUP1'))
	set(tmp, NULL,'OCCUP1_',tmp[, as.integer(OCCUP1_)])
	rh		<- merge(rh,tmp,by='OCCUP1_')
	setnames(tmp, c('OCCUP1_','OCCUP1'), c('OCCUP2_','OCCUP2') )
	rh		<- merge(rh,tmp,by='OCCUP2_')	
	set(rh, NULL, c('OCCUP1_','OCCUP2_'), NULL)
	setnames(tmp, c('OCCUP2_','OCCUP2'), c('OCCUP1_','OCCUP1') )
	#	handle ra
	ra		<- merge(ra,tmp,by='OCCUP1_')
	setnames(tmp, c('OCCUP1_','OCCUP1'), c('OCCUP2_','OCCUP2') )
	ra		<- merge(ra,tmp,by='OCCUP2_')	
	set(ra, NULL, c('OCCUP1_','OCCUP2_'), NULL)
	tmp		<- rh[, which(OCCUP1=='Unknown' & OCCUP2!='Unknown')]
	set(rh, tmp, 'OCCUP1', rh[tmp, OCCUP2])
	set(rh, tmp, 'OCCUP2', 'Unknown')
	tmp		<- ra[, which(OCCUP1=='Unknown' & OCCUP2!='Unknown')]
	set(ra, tmp, 'OCCUP1', ra[tmp, OCCUP2])
	set(ra, tmp, 'OCCUP2', 'Unknown')	
	#	condense OCCUP1 and OCCUP2
	for(x in c('OCCUP1','OCCUP2'))
	{
		set(rh, which(rh[[x]]%in%c('Boda Boda','Trucker')), x, 'Boda/Trucking')
		set(rh, which(rh[[x]]%in%c('Government/clerical/teaching','Military/police','Medical worker')), x, 'Government/clerical/related')
		set(rh, which(rh[[x]]%in%c('Trading/vending','Shopkeeper','Hair dresser/Salon owner')), x, 'Trading/Shopkeeper/Hair')
		set(rh, which(rh[[x]]%in%c('Agriculture for home use/barter','Agriculture for selling','Housekeeper','Housework in your own home','Home brewing')), x, 'Agro/House')
		set(rh, which(rh[[x]]%in%c('Waitress/Waiter/restaurant owner','Bar worker or owner')), x, 'Bar/waitress')
		set(rh, which(rh[[x]]%in%c('Casual laborer','Unemployed')), x, 'Casual laborer/unemployed')
		set(rh, which(rh[[x]]%in%c('Construction (brick maker, builder, porter, painter, roofing)','Mechanic (automobiles, bicycles, electronics)')), x, 'Construction/Mechanic')
		
		set(ra, which(ra[[x]]%in%c('Boda Boda','Trucker')), x, 'Boda/Trucking')
		set(ra, which(ra[[x]]%in%c('Government/clerical/teaching','Military/police','Medical worker')), x, 'Government/clerical/related')
		set(ra, which(ra[[x]]%in%c('Trading/vending','Shopkeeper','Hair dresser/Salon owner')), x, 'Trading/Shopkeeper/Hair')
		set(ra, which(ra[[x]]%in%c('Agriculture for home use/barter','Agriculture for selling','Housekeeper','Housework in your own home','Home brewing')), x, 'Agro/House')
		set(ra, which(ra[[x]]%in%c('Waitress/Waiter/restaurant owner','Bar worker or owner')), x, 'Bar/waitress')
		set(ra, which(ra[[x]]%in%c('Casual laborer','Unemployed')), x, 'Casual laborer/unemployed')
		set(ra, which(ra[[x]]%in%c('Construction (brick maker, builder, porter, painter, roofing)','Mechanic (automobiles, bicycles, electronics)')), x, 'Construction/Mechanic')
	}
	set(rh, rh[, which(OCCUP1=='No additional occupation' & OCCUP2=='Unknown')], 'OCCUP1', 'Unknown')
	#	refine OCCUP1 OCCUP2 when OCAT==Student
	set(rh, rh[, which(OCAT%in%c('Student'))], 'OCCUP1', 'Student')
	set(rh, rh[, which(OCAT%in%c('Student'))], 'OCCUP2', 'Student')
	#	define own OCCUP at time of RCCS visit
	rh[, OCCUP_OLLI:= OCCUP1]
	for(x in c('Student','Boda/Trucking','Bar/waitress','Client/Sex worker'))
		set(rh, rh[, which(OCCUP1==x | OCCUP2==x )],'OCCUP_OLLI', x)
	#	handle ra
	set(ra, ra[, which(OCCUP1=='No additional occupation' & OCCUP2=='Unknown')], 'OCCUP1', 'Unknown')
	ra[, OCCUP_OLLI:= OCCUP1]
	for(x in c('Student','Boda/Trucking','Bar/waitress','Client/Sex worker'))
		set(ra, ra[, which(OCCUP1==x | OCCUP2==x )],'OCCUP_OLLI', x)
	#	OK this is getting complicated: Kate is using all OCC codes to override OCCUP1 
	#	as deemed sensible
	#	just OCAT for simplicity
	#subset(rh, OCAT=='Bar/waitress' & OCCUP1!='Bar/waitress')
	#subset(rh, RID=='G030852')
	#	define SEXWORK
	set(rh, NULL, 'SEXWORK', rh[, as.character(factor(SEXWORK,levels=c(0,1),labels=c('N','Y')))])
	#	define SEXBAR (either work in bar or have sexual contact with barworker)
	set(rh, NULL, 'SEXBAR', rh[, as.character(factor(SEXBAR,levels=c(0,1),labels=c('N','Y')))])	
	#	extend MARSTAT for rh
	tmp		<- rh[, which((is.na(MARSTAT)|MARSTAT=='Never Married') & (RLTN1=='Current wife (at the time)'|RLTN2=='Current wife (at the time)'|RLTN3=='Current wife (at the time)'|RLTN4=='Current wife (at the time)'))]
	set(rh, tmp, 'EVERMARR', 1L)
	set(rh, tmp, 'CURRMARR', 1L)
	set(rh, tmp, 'MARSTAT', 'Monogamous')
	tmp		<- rh[, which(is.na(MARSTAT) & (RLTN1=='Former wife/consensual partner'|RLTN2=='Former wife/consensual partner'|RLTN3=='Former wife/consensual partner'|RLTN4=='Former wife/consensual partner'))]
	set(rh, tmp, 'EVERMARR', 1L)
	set(rh, tmp, 'PREVMAR', 1L)
	set(rh, tmp, 'MARSTAT', 'Previously Married')
	tmp		<- rh[, which(is.na(MARSTAT) & (!RLTN1%in%c('Unknown','Current consensual partner (at the time)')))]
	set(rh, tmp, 'EVERMARR', 0L)
	set(rh, tmp, 'PREVMAR', 0L)
	set(rh, tmp, 'CURRMARR', 0L)
	set(rh, tmp, 'MARSTAT', 'Never Married')
	set(rh, rh[, which(is.na(MARSTAT))], 'MARSTAT', 'Unknown')
	#	define MARSTAT for ra
	ra[, MARSTAT:='Unknown']
	set(ra, ra[, which(EVERMARR==2)], 'MARSTAT', 'Never Married')
	set(ra, ra[, which(EVERMARR==1 & CURRMARR==2)], 'MARSTAT', 'Previously Married')
	set(ra, ra[, which(EVERMARR==1 & CURRMARR==1)], 'MARSTAT', 'Monogamous')
	set(ra, ra[, which(POLYMAR>1 & !POLYMAR%in%c(97,98))], 'MARSTAT', 'Polygamous')
	#	define SEXYEAR
	set(rh, rh[, which(is.na(SEXYEAR))],'SEXYEAR', 99)
	set(rh, NULL, 'SEXYEAR', rh[, gsub('_[0-9]$','',as.character(factor(SEXYEAR, levels=c(0,1,2,8,99), labels=c('Unknown_1','Y','N','Unknown_2','Unknown_3'))))])
	stopifnot( !nrow(subset(rh, is.na(SEXYEAR))) )
	set(ra, ra[, which(is.na(SEXYEAR))],'SEXYEAR', 99)
	set(ra, NULL, 'SEXYEAR', ra[, gsub('_[0-9]$','',as.character(factor(SEXYEAR, levels=c(0,1,2,8,9,99), labels=c('Unknown_1','Y','N','Unknown_2','Unknown_3','Unknown_4'))))])
	stopifnot( !nrow(subset(ra, is.na(SEXYEAR))) )		
	#	define SEXP1YR 
	setnames(rh, c('SEXP1YR'), c('SEXP1YR_'))	
	rh[, SEXP1YR:= as.character(SEXP1YR_)]
	set(rh, rh[, which(SEXP1YR_==92)],'SEXP1YR','<3')
	set(rh, rh[, which(SEXP1YR_==93)],'SEXP1YR','3+')
	set(rh, rh[, which(SEXP1YR_%in%c(97,98,99) | is.na(SEXP1YR_))],'SEXP1YR','Unknown')
	set(rh, NULL, 'SEXP1YR_', NULL)
	setnames(ra, c('SEXP1YR'), c('SEXP1YR_'))	
	ra[, SEXP1YR:= as.character(SEXP1YR_)]
	set(ra, ra[, which(SEXP1YR_==92)],'SEXP1YR','<3')
	set(ra, ra[, which(SEXP1YR_==93)],'SEXP1YR','3+')
	set(ra, ra[, which(SEXP1YR_%in%c(97,98,99) | is.na(SEXP1YR_))],'SEXP1YR','Unknown')
	set(ra, NULL, 'SEXP1YR_', NULL)	
	#	revisit SEXP1YR based on relationships
	set(rh, rh[, which(SEXP1YR%in%c('Unknown') & RLTN_NAMED==1)], 'SEXP1YR', '1')
	set(rh, rh[, which(SEXP1YR%in%c('Unknown') & RLTN_NAMED==2)], 'SEXP1YR', '2')
	set(rh, rh[, which(SEXP1YR%in%c('Unknown') & RLTN_NAMED==3)], 'SEXP1YR', '3')
	set(rh, rh[, which(SEXP1YR%in%c('Unknown') & RLTN_NAMED==4)], 'SEXP1YR', '3+')
	tmp		<- rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED>0)]
	if( nrow(rh[tmp,]) )
		warning("found SEXP1YR%in%c('0') & RLTN_NAMED>0  --> set to had sex last year, n=", length(tmp))
	set(rh, rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED==1)], 'SEXP1YR', '1')
	set(rh, rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED==2)], 'SEXP1YR', '2')
	set(rh, rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED==3)], 'SEXP1YR', '3')
	set(rh, rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED==4)], 'SEXP1YR', '3+')
	#	revisit SEXYEAR based on SEXP1YR
	set(rh, rh[, which(SEXYEAR%in%c('Unknown') & !SEXP1YR%in%c('0','Unknown'))], 'SEXYEAR', 'Y')
	tmp		<- rh[, which(SEXYEAR%in%c('N') & !SEXP1YR%in%c('0','Unknown'))]
	if( nrow(rh[tmp,]) )
		warning("found SEXYEAR%in%c('N') & !SEXP1YR%in%c('0','Unknown') --> set to had sex last year, n=", length(tmp))	
	set(rh, tmp, 'SEXYEAR', 'Y')		
	set(ra, ra[, which(SEXYEAR%in%c('Unknown') & !SEXP1YR%in%c('0','Unknown'))], 'SEXYEAR', 'Y')
	tmp		<- ra[, which(SEXYEAR%in%c('N') & !SEXP1YR%in%c('0','Unknown'))]
	if( nrow(ra[tmp,]) )
		warning("found SEXYEAR%in%c('N') & !SEXP1YR%in%c('0','Unknown') --> set to had sex last year, n=", length(tmp))	
	set(ra, tmp, 'SEXYEAR', 'Y')		
	#	define SEXP1OUT 
	setnames(rh, c('SEXP1OUT'), c('SEXP1OUT_'))	
	rh[, SEXP1OUT:= as.character(SEXP1OUT_)]
	set(rh, rh[, which(SEXP1OUT_==92)],'SEXP1OUT','<3')
	set(rh, rh[, which(SEXP1OUT_==93)],'SEXP1OUT','3+')
	set(rh, rh[, which(SEXP1OUT_%in%c(97,98,99) | is.na(SEXP1OUT_))],'SEXP1OUT','Unknown')
	set(rh, NULL, 'SEXP1OUT_', NULL)	
	setnames(ra, c('SEXP1OUT'), c('SEXP1OUT_'))	
	ra[, SEXP1OUT:= as.character(SEXP1OUT_)]
	set(ra, ra[, which(SEXP1OUT_==92)],'SEXP1OUT','<3')
	set(ra, ra[, which(SEXP1OUT_==93)],'SEXP1OUT','3+')
	set(ra, ra[, which(SEXP1OUT_%in%c(97,98,99) | is.na(SEXP1OUT_))],'SEXP1OUT','Unknown')
	set(ra, NULL, 'SEXP1OUT_', NULL)	
	#	revisit SEXYEAR based on SEXP1OUT
	set(rh, rh[, which(SEXYEAR%in%c('Unknown') & !SEXP1OUT%in%c('0','Unknown'))], 'SEXYEAR', 'Y')
	tmp		<- rh[, which(SEXYEAR%in%c('N') & !SEXP1OUT%in%c('0','Unknown'))]
	if( nrow(rh[tmp,]) )
		warning("found SEXYEAR%in%c('N') & !SEXP1OUT%in%c('0','Unknown') --> set to had sex last year, n=", length(tmp))
	set(rh, tmp, 'SEXYEAR', 'Y')
	set(ra, ra[, which(SEXYEAR%in%c('Unknown') & !SEXP1OUT%in%c('0','Unknown'))], 'SEXYEAR', 'Y')
	tmp		<- ra[, which(SEXYEAR%in%c('N') & !SEXP1OUT%in%c('0','Unknown'))]
	if( nrow(ra[tmp,]) )
		warning("found SEXYEAR%in%c('N') & !SEXP1OUT%in%c('0','Unknown') --> set to had sex last year, n=", length(tmp))
	set(ra, tmp, 'SEXYEAR', 'Y')	
	#	revisit SEXP1YR based on SEXP1OUT
	set(rh, rh[, which(SEXP1OUT=='3+' & SEXP1YR%in%c('0','1','2','<3','Unknown'))], 'SEXP1YR', '3+')
	set(rh, rh[, which(SEXP1OUT=='<3' & SEXP1YR%in%c('0','Unknown'))], 'SEXP1YR', '<3')
	rh[, DUMMY:= seq_len(nrow(rh))]
	warning("set(rh, rh[, which(RID=='G013746' & VISIT==14)], 'SEXP1OUT', 1) --> think this is typo")
	set(rh, rh[, which(RID=='G013746' & VISIT==14)], 'SEXP1OUT', '1')	
	tmp		<- rh[, which(	!SEXP1OUT%in%c('3+','<3','Unknown') & !SEXP1YR%in%c('3+','<3','Unknown'))]  
	tmp		<- subset(rh[tmp, ], as.numeric(SEXP1OUT)>as.numeric(SEXP1YR))[, DUMMY]
	warning("as.numeric(SEXP1OUT)>as.numeric(SEXP1YR), set to SEXP1OUT, n=", length(tmp))
	set(rh, tmp, 'SEXP1YR', rh[tmp, SEXP1OUT])	
	set(ra, ra[, which(SEXP1OUT=='3+' & SEXP1YR%in%c('0','1','2','<3','Unknown'))], 'SEXP1YR', '3+')
	set(ra, ra[, which(SEXP1OUT=='<3' & SEXP1YR%in%c('0','Unknown'))], 'SEXP1YR', '<3')
	ra[, DUMMY:= seq_len(nrow(ra))]
	tmp		<- ra[, which(	!SEXP1OUT%in%c('3+','<3','Unknown') & !SEXP1YR%in%c('3+','<3','Unknown'))]  
	tmp		<- subset(ra[tmp, ], as.numeric(SEXP1OUT)>as.numeric(SEXP1YR))[, DUMMY]
	warning("as.numeric(SEXP1OUT)>as.numeric(SEXP1YR), set to SEXP1OUT, n=", length(tmp))
	set(ra, tmp, 'SEXP1YR', ra[tmp, SEXP1OUT])	
	#	revisit SEXP1YR based on SEXYEAR
	set(rh, rh[, which(SEXYEAR=='Unknown' & SEXP1YR=='0')], 'SEXP1YR', 'Unknown')
	set(rh, rh[, which(SEXYEAR=='Y' & SEXP1YR=='0')], 'SEXP1YR', 'Unknown')
	set(ra, ra[, which(SEXYEAR=='Unknown' & SEXP1YR=='0')], 'SEXP1YR', 'Unknown')
	set(ra, ra[, which(SEXYEAR=='Y' & SEXP1YR=='0')], 'SEXP1YR', 'Unknown')	
	#	define SEXPEVER
	setnames(rh, c('SEXPEVER'), c('SEXPEVER_'))	
	rh[, SEXPEVER:= as.character(SEXPEVER_)]
	set(rh, rh[, which(SEXPEVER_==92)],'SEXPEVER','<3')
	set(rh, rh[, which(SEXPEVER_==93)],'SEXPEVER','3+')
	set(rh, rh[, which(SEXPEVER_%in%c(97,98,99) | is.na(SEXPEVER_))],'SEXPEVER','Unknown')
	set(rh, NULL, 'SEXPEVER_', NULL)
	stopifnot( !nrow(subset(rh, is.na(SEXPEVER))) )	
	setnames(ra, c('SEXPEVER'), c('SEXPEVER_'))	
	ra[, SEXPEVER:= as.character(SEXPEVER_)]
	set(ra, ra[, which(SEXPEVER_==92)],'SEXPEVER','<3')
	set(ra, ra[, which(SEXPEVER_==93)],'SEXPEVER','3+')
	set(ra, ra[, which(SEXPEVER_%in%c(97,98,99) | is.na(SEXPEVER_))],'SEXPEVER','Unknown')
	set(ra, NULL, 'SEXPEVER_', NULL)
	stopifnot( !nrow(subset(ra, is.na(SEXPEVER))) )	
	#	revisit SEXPEVER based on SEXP1YR
	set(rh, rh[, which(SEXP1YR=='3+' & SEXPEVER%in%c('0','1','2','<3','Unknown'))], 'SEXPEVER', '3+')
	set(rh, rh[, which(SEXP1YR=='<3' & SEXPEVER%in%c('0','Unknown'))], 'SEXPEVER', '<3')	
	tmp		<- rh[, which(	!SEXP1YR%in%c('3+','<3','Unknown') & !SEXPEVER%in%c('3+','<3','Unknown'))]  
	tmp		<- subset(rh[tmp, ], as.numeric(SEXP1YR)>as.numeric(SEXPEVER))[, DUMMY]
	warning("as.numeric(SEXP1YR)>as.numeric(SEXPEVER), set to SEXP1YR, n=", length(tmp))
	set(rh, tmp, 'SEXPEVER', rh[tmp, SEXP1YR])	
	set(ra, ra[, which(SEXP1YR=='3+' & SEXPEVER%in%c('0','1','2','<3','Unknown'))], 'SEXPEVER', '3+')
	set(ra, ra[, which(SEXP1YR=='<3' & SEXPEVER%in%c('0','Unknown'))], 'SEXPEVER', '<3')	
	tmp		<- ra[, which(	!SEXP1YR%in%c('3+','<3','Unknown') & !SEXPEVER%in%c('3+','<3','Unknown'))]  
	tmp		<- subset(ra[tmp, ], as.numeric(SEXP1YR)>as.numeric(SEXPEVER))[, DUMMY]
	warning("as.numeric(SEXP1YR)>as.numeric(SEXPEVER), set to SEXP1YR, n=", length(tmp))
	set(ra, tmp, 'SEXPEVER', ra[tmp, SEXP1YR])	
	#	define EVERSEX 
	set(rh, rh[, which(is.na(EVERSEX))],'EVERSEX', 99)
	set(rh, NULL, 'EVERSEX', rh[, gsub('_[0-9]$','',as.character(factor(EVERSEX, levels=c(0,1,2,3,8,99), labels=c('Unknown_1','Y','N','Unknown_2','Unknown_3','Unknown_4'))))])
	stopifnot( !nrow(subset(rh, is.na(EVERSEX))) )	
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & SEXYEAR=='Y')], 'EVERSEX', 'Y')
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & CURRMARR==1)], 'EVERSEX', 'Y')
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & EVERMARR==1)], 'EVERSEX', 'Y')
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & RLTN_NAMED>0)], 'EVERSEX', 'Y')	
	set(ra, ra[, which(is.na(EVERSEX))],'EVERSEX', 99)
	set(ra, NULL, 'EVERSEX', ra[, gsub('_[0-9]$','',as.character(factor(EVERSEX, levels=c(0,1,2,3,8,99), labels=c('Unknown_1','Y','N','Unknown_2','Unknown_3','Unknown_4'))))])
	stopifnot( !nrow(subset(ra, is.na(EVERSEX))) )	
	set(ra, ra[, which(EVERSEX%in%c('Unknown') & SEXYEAR=='Y')], 'EVERSEX', 'Y')
	set(ra, ra[, which(EVERSEX%in%c('Unknown') & CURRMARR==1)], 'EVERSEX', 'Y')
	set(ra, ra[, which(EVERSEX%in%c('Unknown') & EVERMARR==1)], 'EVERSEX', 'Y')		
	#	revisit EVERSEX based on SEXPEVER
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & !SEXPEVER%in%c('0','Unknown'))], 'EVERSEX', 'Y')
	tmp		<- rh[, which(EVERSEX%in%c('N') & !SEXPEVER%in%c('0','Unknown'))]
	if( nrow(rh[tmp,]) )
		warning("found EVERSEX%in%c('N') & !SEXPEVER%in%c('0','Unknown') --> set to had sex ever, n=", length(tmp))
	set(rh, tmp, 'EVERSEX', 'Y')
	set(ra, ra[, which(EVERSEX%in%c('Unknown') & !SEXPEVER%in%c('0','Unknown'))], 'EVERSEX', 'Y')
	tmp		<- ra[, which(EVERSEX%in%c('N') & !SEXPEVER%in%c('0','Unknown'))]
	if( nrow(ra[tmp,]) )
		warning("found EVERSEX%in%c('N') & !SEXPEVER%in%c('0','Unknown') --> set to had sex ever, n=", length(tmp))
	set(ra, tmp, 'EVERSEX', 'Y')	
	#	revisit SEXPEVER based on EVERSEX
	set(rh, rh[, which(EVERSEX=='Unknown' & SEXPEVER=='0')], 'SEXP1YR', 'Unknown')
	set(rh, rh[, which(EVERSEX=='Y' & SEXPEVER=='0')], 'SEXP1YR', 'Unknown')
	set(ra, ra[, which(EVERSEX=='Unknown' & SEXPEVER=='0')], 'SEXP1YR', 'Unknown')
	set(ra, ra[, which(EVERSEX=='Y' & SEXPEVER=='0')], 'SEXP1YR', 'Unknown')		
	#	check SEXPEVER	
	tmp		<- subset(rh, SEXPEVER=='0' & SEXACTIVE==1)
	if( nrow(tmp) )
		warning("found SEXPEVER=='0' & SEXACTIVE==1 --> report only, n=", nrow(tmp))	
	#	add extra-marital partner to MARSTAT
	tmp	<- rh[, which(MULTIPART>0 & MARSTAT!='Previously Married')]
	set(rh, tmp, 'MARSTAT', rh[tmp, paste0(MARSTAT,' + casual partner')])
	tmp	<- rh[, which(MARSTAT%in%c('Monogamous') & !SEXP1YR%in%c('1','Unknown'))]
	set(rh, tmp, 'MARSTAT', rh[tmp, paste0(MARSTAT,' + casual partner')])	
	tmp	<- rh[, which(MARSTAT%in%c('Polygamous') & !SEXP1YR%in%c('1','Unknown') & POLYMAR<as.numeric(gsub('Unknown','0',gsub('+','',SEXP1YR,fixed=1))) ) ]
	set(rh, tmp, 'MARSTAT', rh[tmp, paste0(MARSTAT,' + casual partner')])
	tmp	<- ra[, which(MARSTAT%in%c('Monogamous') & !SEXP1YR%in%c('1','Unknown'))]
	set(ra, tmp, 'MARSTAT', ra[tmp, paste0(MARSTAT,' + casual partner')])	
	tmp	<- ra[, which(MARSTAT%in%c('Polygamous') & !SEXP1YR%in%c('1','Unknown') & POLYMAR<as.numeric(gsub('Unknown','0',gsub('+','',SEXP1YR,fixed=1))) ) ]
	set(ra, tmp, 'MARSTAT', ra[tmp, paste0(MARSTAT,' + casual partner')])
	#	define ever alcohol use during sex
	set(rh, NULL, 'ALC', rh[, as.character(factor(ALC, levels=c(0,1), labels=c('N','Y')))])
	#	redefine SEXP1YR
	set(rh, rh[, which(!SEXP1YR%in%c('0','1','2','Unknown'))], 'SEXP1YR','3+')
	set(ra, ra[, which(!SEXP1YR%in%c('0','1','2','Unknown'))], 'SEXP1YR','3+')
	#	redefine SEXP1OUT
	set(rh, rh[, which(!SEXP1OUT%in%c('0','1','2','Unknown'))], 'SEXP1OUT','3+')
	set(ra, ra[, which(!SEXP1OUT%in%c('0','1','2','Unknown'))], 'SEXP1OUT','3+')
	#	check circumcision
	tmp		<- rh[, which(CIRCUM=='Y' & SEX=='F')]
	cat('\nWarning: found female circumcised --> set to NA' ,rh[tmp, paste(RID, collapse=' ')])	
	set(rh, tmp, 'CIRCUM', NA_integer_)
	#	define COMM_TYPE
	set(rh, NULL, 'COMM_NUM', rh[, as.character(COMM_NUM)])
	set(ra, NULL, 'COMM_NUM', ra[, as.character(COMM_NUM)])
	tmp		<- data.table(	COMM_NUM=	c("1","2","3","4","5","6","7","8","9","14","15","16","18","19","22","23","24","25","29","32","33","34","35","36","38","40","44","45","46","51","52","53","54","55", "56","57","58","59","60","61","62","65","67","74","77","81","84","89","94","95","103","106","107","108","109","120","177", "183", "256", "370","391","401","451", "468","602", "754", "755", "760", "770","771","772","773","774","776"),
							COMM_TYPE=	c("T","A","A","T","A","A","A","A","A", "A", "A", "T", "A", "A", "T", "A", "T", "A", "A", "A", "T", "A", "A", "A", "F", "A", "A", "A", "A", "T", "A", "A", "A", "A",  "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",  "A","A",   "A",  "T",  "A",  "A",  "A",  "A",   "A",   "A",   "A",  "A", "A",  "A",    "A",  "A",  "A",    "A",   "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])		
	stopifnot(!length(setdiff( rh[, COMM_NUM], tmp[, COMM_NUM] )))
	stopifnot(!length(setdiff( ra[, COMM_NUM], tmp[, COMM_NUM] )))
	rh		<- merge(rh, tmp, by='COMM_NUM')
	ra		<- merge(ra, tmp, by='COMM_NUM')	
	#	merge community numbers for same community
	set(rh, NULL, 'COMM_NUM', rh[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',COMM_NUM))))])
	set(ra, NULL, 'COMM_NUM', ra[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',COMM_NUM))))])
	set(rd, NULL, 'COMM_NUM', rd[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',COMM_NUM))))])	
	#	set to NULL
	set(rh, NULL, c('VDEX','EVERMARR','CURRMARR','RELIGION','POLYMAR','DUMMY','RLTN1','RLTN2','RLTN3','RLTN4','RLTN_NAMED'), NULL)
	set(rh, NULL, c('BVDEX','EVERSEX','SEXGIFT','SEXYEAR','EDUCATE','EDUCYRS','EDCAT','ARVMED','CNDEVER1','RNYRCON1','CNDEVER2','RNYRCON2','CNDEVER3','RNYRCON3','CNDEVER4','RNYRCON4','RLTNCON1'),NULL)
	set(rh, NULL, c('RLTNCON2','RLTNCON3','RLTNCON4','ALC1B','ALC2B','ALC3B','ALC4B','ALC1F','ALC2F','ALC3F','ALC4F','OCCUP11','OCCUP12','OCCUP13','OCCUP14','OCCUP21','OCCUP22','OCCUP23','OCCUP24','SEXHIGH','SEXOUT'),NULL)
	set(rh, NULL, c('SEXCAT','PREVMAR','AGECAT','AGECAT2','HIVPREV2','UNDER25','AGE15TO19','AGE20TO24','AGE25TO29','AGE30TO34','AGE35TO39','AGE40TO44','AGE45TO49','OCCLAG1','SUM_ALC'),NULL)
	set(rh, NULL, c('SEXACTIVE','MULTIPART','CAS','SUMCON','CONCON','NEVERSEX','OCCUP1','OCCUP2','SEXPEVER'),NULL)
	set(rn, NULL, c('SAMPLEREASON'), NULL)
	rd[, COHORT:= 'RCCS']
	rh[, COHORT:= 'RCCS']
	ra[, COHORT:= 'RCCS']
	#
	list(rd=rd, rh=rh, ra=ra, rn=rn)
}

RakaiCirc.seq.get.info<- function()
{
	infile				<- "~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda"	
	infile.sangerstats	<- "~/Dropbox (Infectious Disease)/PANGEA_data/2016-07-07_PANGEA_stats_by_sample.csv"
	infile.relabel		<- "~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional phylogenetic analyses/Region 1 gag analysis/SummaryofGagSequenceData.rda"
	infile.assembly		<- "~/Dropbox (Infectious Disease)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
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
	#	add Sanger IDs
	#
	tmp		<- as.data.table(read.csv(infile.assembly, stringsAsFactors=FALSE))
	setnames(tmp, c('Sanger.ID','PANGEA.ID'), c('SID','PID'))
	set(tmp, NULL, 'PID', tmp[, gsub('-S[0-9]+','',PID)])
	set(tmp, NULL, 'PID', tmp[, gsub('^\\s','',PID)])
	rs		<- merge(rs, subset(tmp, select=c(PID, SID)), by='PID', all.x=1)
	#	NOTE: there are a few PID with multiple SIDs !
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

RakaiCirc.seq.get.info.PANGEA.170505<- function()
{
	indir	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301'
	#	start with latest Sanger IDs
	dc		<- as.data.table(read.csv('~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/WTSI_PANGEA_InfoFind_2017-02-14.csv', header=TRUE, stringsAsFactors=FALSE))
	setnames(dc, c('Lane','Public'), c('SID','PIDF'))
	dc		<- subset(dc, select=c(SID,PIDF))
	set(dc, NULL, 'SID', dc[, gsub('#','_',SID)])
	set(dc, NULL, 'PID', dc[, gsub('-S[0-9]+$','',PIDF)])
	#
	#	list of files that Chris has
	#
	infiles	<- data.table(F=list.files(indir, pattern='^CW_PANGEA',full.names=TRUE))
	tmp		<- do.call('rbind',lapply(infiles[, F], function(file){
						tmp<- as.data.table(read.csv(file, header=FALSE, col.names='PID', stringsAsFactors=FALSE))
						tmp[, F:=file]
						tmp
					}))
	set(tmp, NULL, 'PROC_STATUS', tmp[, gsub('CW_PANGEA_Rakai_','',gsub('\\.txt','',basename(F)))])	
	set(tmp, NULL, 'F', NULL)
	#	all files for which we have some PID
	dc		<- merge(dc, tmp, by='PID',all=1)	
	stopifnot(!nrow(subset(dc, is.na(SID) & PROC_STATUS!='ThoseWithoutFastqs')))
	dc		<- subset(dc, is.na(SID) | SID!='15351_1_1')		#	remove 15351_1_1  with no PANGEA info in WTSI file
	dc		<- subset(dc, is.na(SID) | SID!='15430_1_75')	#	remove 15430_1_75  with no PANGEA info in WTSI file
	stopifnot(!nrow(subset(dc, is.na(PIDF) & !is.na(SID))))
	dc[, PART:= as.numeric(gsub('^[0-9]+_([0-9])_[0-9]+','\\1',SID))]
	dc[, DUMMY:= gsub('^([0-9]+)_[0-9]_([0-9]+)','\\1_x_\\2',SID)]	
	dc		<- merge(dc, dc[, list(N_PART=length(PART), S_PART=sum(PART)), by='DUMMY'], by='DUMMY')	
	#	check if we always have _1_ and _2_
	stopifnot(!nrow(subset(dc, N_PART==2 & S_PART!=3)))
	#	assume Tanya merges to _3_
	tmp		<- dc[, which(N_PART==2)]
	set(dc, tmp, 'SID', dc[tmp, gsub('^([0-9]+)_[0-9]_([0-9]+)','\\1_3_\\2',SID)])
	set(dc, NULL, c('PART','N_PART','S_PART','DUMMY'),NULL)
	dc		<- unique(dc)
	tmp		<- subset(dc, !is.na(SID))[, list(N_SID=length(SID)), by='PIDF']
	dc		<- merge(dc, tmp, by='PIDF',all.x=1)
	#	define controls
	set(dc, dc[, which(grepl('neg',PID))], 'PROC_STATUS', 'NegControl')
	#	define not in Chris census
	set(dc, dc[, which(is.na(PROC_STATUS))], 'PROC_STATUS', 'NotTrackedByChris')	
	#	add extra category to PROC_STATUS: not processed by Kate
	tmp		<- as.data.table(read.table("~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/KG_PANGEA_Processed_597.txt", header=TRUE,stringsAsFactors=FALSE))
	setnames(tmp, c('SampleID','LaneID'), c('PID','SID'))
	tmp[, KATE_PROC:='Y']
	dc		<- merge(dc, tmp, by=c('PID','SID'), all.x=1)	
	set(dc, dc[, which(PROC_STATUS=='ThoseWithFastqs_WithKateShiverOutput' & is.na(KATE_PROC))], 'PROC_STATUS','ThoseWithFastqs_KateNotProcessed') 
	set(dc, NULL, 'KATE_PROC', NULL)	
	#	add RIDs
	load('~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda')
	tmp		<- subset(as.data.table(rccsData), select=c(RCCS_studyid,Pangea.id, batch, date, SEX))
	setnames(tmp, c('RCCS_studyid','Pangea.id','batch','date'), c('RID','PID','RCCS_SHIP_BATCH','SAMPLE_DATE'))
	tmp		<- subset(tmp, !is.na(PID))
	tmp		<- unique(tmp, by=c('RID','PID'))
	tmp2	<- subset(as.data.table(neuroData), select=c(studyid, Pangea.id, sampleDate, gender))
	setnames(tmp2, c('studyid','Pangea.id','sampleDate','gender'), c('RID','PID','SAMPLE_DATE','SEX'))
	tmp2[, RCCS_SHIP_BATCH:='neuro']
	tmp		<- rbind(tmp, tmp2, use.names=TRUE)
	set(tmp, NULL, 'RID', tmp[, as.character(RID)])
	set(tmp, NULL, 'PID', tmp[, as.character(PID)])
	set(tmp, NULL, 'SEX', tmp[, as.character(SEX)])
	dc		<- merge(dc, tmp, by='PID',all.x=1)
	#	flag test plate
	set(dc, dc[, which(grepl('PG14-UG9000[0-9][0-9]',PID))], 'RCCS_SHIP_BATCH', 'test')
	#	see if on HPC	
	tmp		<- as.data.table(read.csv('~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/HPC_census_bams.txt', header=FALSE, col.names='SID', stringsAsFactors=FALSE))
	tmp[, HPC_BAM:='Y']
	set(tmp, NULL, 'SID', tmp[, gsub('\\.bam','',SID)])
	tmp		<- subset(tmp, grepl('^[0-9]+_[0-9]_[0-9]+',SID))	#	reduce to bams from SANGER
	stopifnot( !length(setdiff( tmp[, SID], dc[, SID] )) )
	dc		<- merge(dc, tmp, by='SID',all.x=TRUE)
	set(dc, dc[, which(is.na(HPC_BAM))], 'HPC_BAM', 'N')	
	tmp		<- as.data.table(read.csv('~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/HPC_census_refs.txt', header=FALSE, col.names='SID', stringsAsFactors=FALSE))
	tmp[, HPC_REF:='Y']
	set(tmp, NULL, 'SID', tmp[, gsub('_ref.fasta','',SID)])
	tmp		<- subset(tmp, grepl('^[0-9]+_[0-9]_[0-9]+',SID))	#	reduce to refs from SANGER	
	stopifnot( !length(setdiff( tmp[, SID], dc[, SID] )) )	
	dc		<- merge(dc, tmp, by='SID',all.x=TRUE)
	set(dc, dc[, which(is.na(HPC_REF))], 'HPC_REF', 'N')
	tmp		<- as.data.table(read.csv('~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/HPC_census_fastq.txt', header=FALSE, col.names='SID', stringsAsFactors=FALSE))
	tmp[, HPC_FASTQ:='Y']
	set(tmp, NULL, 'SID', tmp[, gsub('_[0-9].fastq.gz','',SID)])	
	tmp		<- subset(tmp, grepl('^[0-9]+_[0-9]_[0-9]+',SID))	#	reduce to fastqz's from SANGER	
	dc		<- merge(dc, unique(tmp), by='SID',all.x=TRUE)
	set(dc, dc[, which(is.na(HPC_FASTQ))], 'HPC_FASTQ', 'N')	
	#	add latest PANGEA stats
	tmp		<- as.data.table(read.csv('~/Dropbox (Infectious Disease)/PANGEA_data/2016-07-07_PANGEA_stats_by_sample.csv', stringsAsFactors=FALSE))
	setnames(tmp, 	c('Status','Submitted','DayDiff','ProjectID','Cohort','Submitted.1','Sequenced','Assembled','HIVcontig'), 
			c('WTSI_STATUS','WTSI_SUBMITTED_DATE','DayDiff','PID','Cohort','WTSI_SUBMITTED','WTSI_SEQUENCED','WTSI_ASSEMBLED','WTSI_HIVCONTIG'))
	tmp		<- subset(tmp, PID!='')
	set(tmp, NULL, 'WTSI_SUBMITTED', tmp[, as.character(factor(WTSI_SUBMITTED=='',levels=c(TRUE,FALSE),labels=c('N','Y')))])
	set(tmp, NULL, 'WTSI_SEQUENCED', tmp[, as.character(factor(WTSI_SEQUENCED=='',levels=c(TRUE,FALSE),labels=c('N','Y')))])
	set(tmp, NULL, 'WTSI_ASSEMBLED', tmp[, as.character(factor(WTSI_ASSEMBLED=='',levels=c(TRUE,FALSE),labels=c('N','Y')))])
	set(tmp, NULL, 'WTSI_HIVCONTIG', tmp[, as.character(factor(WTSI_HIVCONTIG=='',levels=c(TRUE,FALSE),labels=c('N','Y')))])
	set(tmp, NULL, c('DayDiff','Cohort'), NULL)
	dc		<- merge(dc, tmp, by='PID', all.x=1)
	#	check what Dan assembled HISEQ
	tmp		<- as.data.table(read.csv('~/Dropbox (Infectious Disease)/PANGEA_data/PANGEA_UCL_Feb2017_collated_stats_all_genomes_UCL_release_Feb2017_allhitodate.csv', stringsAsFactors=FALSE))
	setnames(tmp, c('PANGEA_ID','PGID_full','WTSI_ID','Length', 'Cohort'), c('PID','PIDF','SID','UCL_LEN', 'COHORT'))
	tmp		<- subset(tmp, select=c('PID','PIDF','SID','UCL_LEN', 'COHORT'))
	#	convert SID to _3_ to match our convention
	set(tmp, NULL, 'SID', tmp[, gsub('^([0-9]+)_[0-9]_([0-9]+)','\\1_3_\\2',SID)])
	#	check what Dan assembled MISEQ
	tmp2	<- as.data.table(read.csv('~/Dropbox (Infectious Disease)/PANGEA_data/PANGEA_UCL_Feb2017_collated_stats_all_genomes_UCL_release_Feb2017_allmitodate.csv', stringsAsFactors=FALSE))
	setnames(tmp2, c('PANGEA_ID','PGID_full','WTSI_lane_ID','Length', 'Cohort'), c('PID','PIDF','SID','UCL_LEN', 'COHORT'))
	tmp2	<- subset(tmp2, select=c('PID','PIDF','SID','UCL_LEN', 'COHORT'))	
	tmp		<- rbind(tmp, tmp2)
	stopifnot( !length(setdiff(tmp[, SID], dc[, SID])) )	
	#	n=6 plus test plate
	if(0)
	{
		tmp		<- subset(tmp, grepl('Rakai',COHORT))
		tmp2	<- setdiff( tmp[, unique(sort(PID))], dc[, unique(sort(PID))] )
		write.csv(data.table(ID=tmp2), file=file.path(indir,'IDs_that_Dan_flags_from_Rakai_but_not_in_Chris_census.csv'), row.names=FALSE)			
	}
	set(tmp, NULL, 'COHORT', NULL)		
	dc		<- merge(dc, tmp, by=c('PID','PIDF','SID'),all.x=1)	
	#	resolve 60 from Chris that he confirmed he has not processed -- at least one SID from said individual has been processed
	set(dc, dc[, which(PROC_STATUS=='ThoseWithFastqs_WithChrisShiverOutput' & HPC_BAM=='N')], 'PROC_STATUS','ThoseWithFastqs_ChrisNotProcessed')
	#
	#	subset to RAKAI
	#
	dc		<- subset(dc, !is.na(RID))
	#	define RIDs from whom all SIDs are complete and on HPC
	tmp		<- dc[, list(	HPC_ALL_SID_FOR_RID= all(HPC_BAM=='Y') & all(HPC_REF=='Y')), by='RID']	
	set(tmp, NULL, 'HPC_ALL_SID_FOR_RID', tmp[,as.character(factor(HPC_ALL_SID_FOR_RID, levels=c(TRUE,FALSE), labels=c('Y','N')))])
	dc		<- merge(dc, tmp, by='RID',all.x=1)
	#	define Sampling Time
	set(dc, NULL, 'WTSI_SUBMITTED_DATE', dc[, as.Date(WTSI_SUBMITTED_DATE)])	
	stopifnot(!nrow(subset(dc, HPC_BAM!=HPC_REF)))
	#	reduce to currently on HPC 
	rs		<- subset(dc, !is.na(SID), select=c(RID, PID, PIDF, SID, RCCS_SHIP_BATCH, SAMPLE_DATE))
	set(rs, NULL, 'SAMPLE_DATE', rs[, hivc.db.Date2numeric(SAMPLE_DATE)])
	#	add VISIT for those individuals in pop surveillance
	infile	<- "~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda"
	load(infile)
	rd		<- as.data.table(rccsData)
	setnames(rd, colnames(rd), gsub('\\.','_',toupper(colnames(rd))))	
	for(x in colnames(rd))
		if(class(rd[[x]])=='Date')
			set(rd, NULL, x, hivc.db.Date2numeric(rd[[x]]))
	setnames(rd, c('RCCS_STUDYID','PANGEA_ID','DATE'), c('RID','PID','SAMPLE_DATE'))
	rd		<- subset(rd, !is.na(PID), select=c(RID,PID,VISIT,SAMPLE_DATE))
	#	round sample dates since we need to merge by those	
	set(rs, NULL, 'SAMPLE_DATE', rs[, round(SAMPLE_DATE,d=3)])
	set(rd, NULL, 'SAMPLE_DATE', rd[, round(SAMPLE_DATE,d=3)])
	rs		<- merge(rs, rd, by=c('RID','PID','SAMPLE_DATE'), all=1)
	rs		<- subset(rs, !is.na(SID))
	#
	#	save to file
	#
	save(rs, file=file.path(wdir,'RCCS_SeqInfo_170505.rda'))
}

RakaiCirc.various<- function()
{
	require(ape)
	require(data.table)
	if(0)
	{
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
	if(1)
	{
		require(ape)
		require(data.table)
		
		wdir	<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
		wdir	<- '/work/or105/Gates_2014/Rakai'
		infile	<- file.path(wdir, "PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda")
		load(infile)
		#	prepare genetic distance matrix
		sq			<- as.character(sq)
		sq[sq=='?']	<- '-'
		sq			<- as.DNAbin(sq)
		sq.gd		<- dist.dna(sq, pairwise.deletion=TRUE)
		sq.gd		<- as.data.table(melt(as.matrix(sq.gd), varnames=c('TAXA1','TAXA2')))
		setnames(sq.gd, 'value', 'PD')
		sq.gd		<- subset(sq.gd, TAXA1!=TAXA2)
		set(sq.gd, NULL, 'TAXA1', sq.gd[, as.character(TAXA1)])
		set(sq.gd, NULL, 'TAXA2', sq.gd[, as.character(TAXA2)])		
		#	add overlap	
		setkey(sq.gd, TAXA1, TAXA2)
		tmp			<- as.character(sq)
		tmp			<- !( tmp=='-' | tmp=='?' | tmp=='n' )
		tmp[]		<- as.integer(tmp) 
		tmp2		<- sq.gd[, list(OVERLAP= sum(bitwAnd(tmp[TAXA1,], tmp[TAXA2,])) ), by=c('TAXA1','TAXA2')]	
		sq.gd		<- merge(sq.gd, tmp2, by=c('TAXA1','TAXA2'))
		
		save(sq.gd, file=gsub('\\.rda','_gd.rda',infile))		
	}	
}

RakaiCirc.seq.get.phylogenies<- function()
{
	require(ape)
	require(data.table)
	wdir	<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	#
	#	start with Kate s GP24 FastTree
	#
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
	#
	#	also use the PANGEA ExaML tree
	#
	infile	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_fasttree.rda'
	load(infile)	#loads "ph", "dist.brl", "ph.gdtr", "ph.mrca", "clustering"
	php		<- ph
	tmp		<- copy(ph.gdtr)
	tmp[upper.tri(tmp, diag=TRUE)]	<- NA_real_
	tmp								<- as.data.table(melt(tmp))
	setnames(tmp, c('Var1','Var2','value'), c('PID','PID2','PD'))
	phpi							<- subset(tmp, !is.na(PD))
	#
	#	get raw genetic distances on latest Region1 alignment
	#
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
	#	this takes about 2 days.. ..ran on cluster.
	tmp			<- (as.character(sq)!='-')
	setkey(gdgag, IDX, IDX2)
	tmp2		<- gdgag[, list(OVERLAP= sum(apply(tmp[c(IDX,IDX2),],2,all)) ), by=c('IDX','IDX2')]
	gdgag		<- merge(gdgag, tmp2, by=c('SEQIDc','SEQIDc2'))
	save(phf, php, phfi, phpi, sq, gdgag, file=file.path(wdir,'RCCS_PhInfo_160825.rda'))	
	#
	#	plot nucleotide frequencies (unfinished)
	#
	if(0)
	{
		sqf			<- data.table(SITE=seq_len(ncol(sq)))
		sqf			<- sqf[, {
					tmp	<- base.freq(sq[,SITE], freq=TRUE, all=TRUE)
					list(NT=paste('NT',names(tmp),sep=''), N=tmp)
				}, by='SITE']
		sqf			<- subset(sqf, N>0)
		tmp			<- subset(sqf,!NT%in%c('NTn','NT-','NT?'))[, list(TOTAL_SEQ=sum(N)), by='SITE']
		sqf			<- merge(sqf, tmp, all.x=1, by='SITE')
		sqf[, NTc:= 'NTother']
		tmp			<- sqf[, which(NT%in%c('NTa','NTt','NTg','NTc','NT-','NT?'))]
		set(sqf, tmp, 'NTc', sqf[tmp, NT])
		set(sqf, sqf[, which(NT=='NTn')], 'NTc', 'NT?')
		ggplot(subset(sqf, SITE<1000 & NTc!='NT?' & NTc!='NT-'), aes(x=SITE, fill=NTc)) + 
				geom_bar(position='fill') +
				scale_fill_brewer(palette='Set1')
	}
}


RakaiCirc.circ.dev160907<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	
	#	get epi info
	tmp		<- RakaiCirc.epi.get.info()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#	get sequence info and add to recipients
	load(file.path(wdir,'RCCS_SeqInfo_160816.rda'))		
	rs		<- subset(rs, !is.na(VISIT))	
	#	load female recipients
	rrec	<- RakaiCirc.recipient.female.get.info()
	#	add sequence info to recipients
	rrec	<- merge(rrec, unique(subset(rs, select=c(RID, PID, SID))), by='RID', all.x=1)
	#	select female recipients with at least one PANGEA sequence
	rrec	<- subset(rrec, !is.na(SID))
	
	
	#	select likely pairs, use same selection as before just for consistency
	select.discsib	<- 0.65	
	#
	#	collect runs
	#
	infiles	<- data.table(	F_TRM= c(	'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w200/RCCS_160902_w200_trmStats.rda',
										'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w220/RCCS_160902_w220_trmStats.rda',
										'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w250/RCCS_160902_w250_trmStats.rda',
										'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w270/RCCS_160902_w270_trmStats.rda',
										'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w280/RCCS_160902_w280_trmStats.rda'
							))
	infiles[, F_PH:= gsub('trmStats.rda','trees.rda', F_TRM)]
	infiles[, DIR:= dirname(F_TRM)]
	infiles[, RUN:= gsub('_trmStats.rda','',basename(F_TRM))]				
	#
	#	for each run: get list of pairs
	#	
	rp		<- infiles[, {
				#F_TRM	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w270/RCCS_160902_w270_trmStats.rda'
				#F_PH	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w270/RCCS_160902_w270_trees.rda'
				load(F_TRM)	#loads df
				dlkl	<- copy(df)
				load(F_PH)	#loads phs and dfr
				#	select likely pairs
				dlkl[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]
				tmp		<- dcast.data.table(dlkl, PAIR_ID~TYPE, value.var='WIN_OF_TYPE_P')
				set(tmp, tmp[, which(is.na(disconnected))], 'disconnected', 0)
				set(tmp, tmp[, which(is.na(sib))], 'sib', 0)
				tmp		<- subset(tmp, disconnected+sib < select.discsib, PAIR_ID)
				dlkl	<- merge(dlkl, tmp, by='PAIR_ID')
				#	get likely pairs that involve female recipients
				tmp		<- copy(dlkl)
				setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
				set(tmp, tmp[, which(TYPE=='anc_12')], 'TYPE', 'anc')
				set(tmp, tmp[, which(TYPE=='anc_21')], 'TYPE', 'anc_12')
				set(tmp, tmp[, which(TYPE=='anc')], 'TYPE', 'anc_21')
				tmp		<- rbind(tmp, dlkl)
				setnames(tmp, 'ID1', 'SID')
				rp		<- merge(tmp, rrec, by='SID')
				setnames(rp, 'ID2', 'TR_SID')
				cat('\nFound female recipients among likely pairs, n=', rp[, length(unique(RID))])
				#	get info on transmitters: RID and PID
				tmp		<- subset(rs, !is.na(SID), select=c(RID, PID, SID))
				setkey(tmp, SID)
				tmp		<- unique(tmp)
				setnames(tmp, colnames(tmp), paste('TR_',colnames(tmp),sep=''))	
				#
				#	TODO check whereabouts of 12559_1_9 15065_1_32 15034_1_79 15081_1_71 15430_1_73
				#
				tmp2	<- setdiff( rp[, unique(TR_SID)], tmp[, TR_SID] )
				if(length(tmp2))
					cat('\nWarning: Sanger ID in rp that is not in dictionary', tmp2)
				rp		<- merge(rp, tmp, by='TR_SID')
				#	get info on transmitters: demographic stuff	
				tmp		<- subset(rd, !is.na(PID), select=c(RID, SEX, REGION, COMM_NUM, HH_NUM, BIRTHDATE, LASTNEGVIS, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, RELIGION))
				setkey(tmp, RID)
				tmp		<- unique(tmp)
				setnames(tmp, colnames(tmp), paste('TR_',colnames(tmp),sep=''))
				rp		<- merge(rp,tmp,by='TR_RID')
				#	get info on sequences: date of sampling
				tmp		<- unique(subset(rs, !is.na(SID), select=c(SID, DATE)))
				setnames(tmp, 'DATE', 'SEQ_DATE')
				rp		<- merge(rp, tmp, by='SID')
				setnames(tmp, c('SID','SEQ_DATE'), c('TR_SID','TR_SEQ_DATE'))
				rp		<- merge(rp, tmp, by='TR_SID')				
				
				rp
			}, by=c('RUN','DIR','F_TRM','F_PH')]
	#
	#	for each run: plot pairs & trees	
	#
	setkey(rp, RUN, PAIR_ID, SID, TR_SID)
	for( run in rp[, unique(RUN)] )
	{		
		#run		<- 'RCCS_160902_w270'
		dir		<- subset(rp, RUN==run)[1,DIR]
		cat('\ndir is',dir,'\trun is',run)
		df		<- subset(rp, RUN==run)
		rpok	<- subset(df, TR_SEX=='M')	
		rpff	<- subset(df, TR_SEX=='F')
		rpoku	<- unique(rpok)
		rpffu	<- unique(rpff)
		
		#	for every F-F pair, keep only one instance
		#	and select different individuals
		rpffu		<- merge(rpffu, rpffu[, list(SID=SID[1],TR_SID=TR_SID[1]), by='PAIR_ID'], by=c('PAIR_ID','SID','TR_SID'))	
		rpffu		<- subset(rpffu, RID!=TR_RID)		
		#
		#	plot evidence
		#		
		tmp		<- copy(rpoku)
		tmp		<- tmp[order(-PAIR_ID),]
		tmp[, LABEL:= factor(PAIR_ID, levels=PAIR_ID, labels=paste('Pair',PAIR_ID, 	'\nPerson1 ', RID, ' ', SID,' -sex:',SEX,' -loc:',REGION,',',COMM_NUM,',',HH_NUM,' -birth:',BIRTHDATE,' -neg:',LASTNEGDATE,' -pos:',FIRSTPOSDATE,' -seq:',SEQ_DATE,
								'\n<->', 
								'\nPerson2 ', TR_RID, ' ', TR_SID,' -sex:',TR_SEX,' -loc:',TR_REGION,',',TR_COMM_NUM,',',TR_HH_NUM,' -birth:',TR_BIRTHDATE,' -neg:',TR_LASTNEGDATE,' -pos:',TR_FIRSTPOSDATE,' -seq:',TR_SEQ_DATE,																				
								'\n',sep=''))]
		tmp		<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
		set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('anc_12','anc_21','sib','int','disconnected'), labels=c('1 ancestral to 2','2 ancestral to 1','1, 2 are siblings','1, 2 are intermingled','1, 2 are disconnected'))])
		ggplot(tmp, aes(x=LABEL, y=WIN_OF_TYPE, fill=TYPE)) +
				geom_bar(stat='identity', position='stack') +
				coord_flip() +
				labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
				scale_fill_manual(values=c('1 ancestral to 2'="#9E0142",'2 ancestral to 1'="#F46D43",'1, 2 are siblings'="#ABDDA4",'1, 2 are intermingled'="#3288BD",'1, 2 are disconnected'='grey50')) +
				theme_bw() + theme(legend.position='top') +
				guides(fill=guide_legend(ncol=2))
		ggsave(file=file.path(dir, paste(run,'_RCCS_lkltransmissionpairs_withfemaleseroconverter_discsiblt',100*select.discsib,'.pdf',sep='')), w=15, h=max(3,0.23*nrow(tmp)), limitsize = FALSE)	
		#	look at F-F pairs
		tmp			<- copy(rpffu)
		tmp			<- tmp[order(-PAIR_ID),]
		tmp[, LABEL:= factor(PAIR_ID, levels=PAIR_ID, labels=paste('Pair',PAIR_ID, 	'\nPerson1 ', RID, ' ', SID,' -sex:',SEX,' -loc:',REGION,',',COMM_NUM,',',HH_NUM,' -birth:',BIRTHDATE,' -neg:',LASTNEGDATE,' -pos:',FIRSTPOSDATE,' -seq:',SEQ_DATE,
								'\n<->', 
								'\nPerson2 ', TR_RID, ' ', TR_SID,' -sex:',TR_SEX,' -loc:',TR_REGION,',',TR_COMM_NUM,',',TR_HH_NUM,' -birth:',TR_BIRTHDATE,' -neg:',TR_LASTNEGDATE,' -pos:',TR_FIRSTPOSDATE,' -seq:',TR_SEQ_DATE,																				
								'\n',sep=''))]	
		tmp		<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
		set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('anc_12','anc_21','sib','int','disconnected'), labels=c('1 ancestral to 2','2 ancestral to 1','1, 2 are siblings','1, 2 are intermingled','1, 2 are disconnected'))])
		ggplot(tmp, aes(x=LABEL, y=WIN_OF_TYPE, fill=TYPE)) +
				geom_bar(stat='identity', position='stack') +
				coord_flip() +
				labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
				scale_fill_manual(values=c('1 ancestral to 2'="#9E0142",'2 ancestral to 1'="#F46D43",'1, 2 are siblings'="#ABDDA4",'1, 2 are intermingled'="#3288BD",'1, 2 are disconnected'='grey50')) +
				theme_bw() + theme(legend.position='top') +
				guides(fill=guide_legend(ncol=2))
		ggsave(file=file.path(dir, paste(run,'_RCCS_lkltransmissionpairs_femalefemale_discsiblt',100*select.discsib,'.pdf',sep='')), w=15, h=max(3,0.2*nrow(tmp)), limitsize = FALSE)				
	}
	#
	#	re-examine phylogenies
	#
	setkey(rp, RUN, PAIR_ID, SID, TR_SID)
	for( run in rp[, unique(RUN)] )
	{
		dir		<- subset(rp, RUN==run)[1,DIR]
		cat('\ndir is',dir,'\trun is',run)
		df		<- subset(rp, RUN==run)
		rpok	<- subset(df, TR_SEX=='M')	
		rpff	<- subset(df, TR_SEX=='F')
		rpoku	<- unique(rpok)
		rpffu	<- unique(rpff)
		rpffu	<- merge(rpffu, rpffu[, list(SID=SID[1],TR_SID=TR_SID[1]), by='PAIR_ID'], by=c('PAIR_ID','SID','TR_SID'))	
		rpffu	<- subset(rpffu, RID!=TR_RID)		
		
		load( df[1, F_PH] )	#loads phs and dfr
		setkey(rpok, PAIR_ID, SID, TR_SID)
		rpoku	<- unique(rpok)
		invisible(sapply(seq_len(nrow(rpoku)), function(ii)
						{								
							pair.id		<- rpoku[ii, PAIR_ID]
							pty.run		<- rpoku[ii, PTY_RUN]
							dfs			<- subset(dfr, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
							dfs[, TITLE:= dfs[, paste('pair', pair.id, 'ids', rpoku[ii, TR_SID], rpoku[ii, SID], '\nrun', pty.run, '\nwindow', W_FROM,'-', W_TO)]]			
							plot.file	<- file.path(dir, paste(run,'_RCCS_lkltransmissionpairs_withfemaleseroconverter_discsiblt65_pair',pair.id,'.pdf',sep=''))			
							invisible(phsc.plot.selected.pairs(phs, dfs, rpoku[ii, TR_SID], rpoku[ii, SID], plot.file=plot.file, pdf.h=50, pdf.rw=10))
						}))		
		tmp		<- copy(rpffu)
		invisible(sapply(seq_len(nrow(tmp)), function(ii)
						{
							pair.id		<- tmp[ii, PAIR_ID]
							pty.run		<- tmp[ii, PTY_RUN]
							dfs			<- subset(dfr, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
							dfs[, TITLE:= dfs[, paste('pair', pair.id, 'ids', tmp[ii, TR_SID], tmp[ii, SID], '\nrun', pty.run, '\nwindow', W_FROM,'-', W_TO)]]			
							plot.file	<- file.path(dir, paste(run,'_RCCS_lkltransmissionpairs_femalefemale_discsiblt65_pair',pair.id,'.pdf',sep=''))			
							invisible(phsc.plot.selected.pairs(phs, dfs, tmp[ii, TR_SID], tmp[ii, SID], plot.file=plot.file, pdf.h=50, pdf.rw=10))
						}))
	}

	
	#	the 200 + 270 phylogenies look identical!?!?!
	subset(rp, grepl('w270',RUN) & PAIR_ID==3)
	subset(dfr, grepl('w270',RUN), select=c(PTY_RUN, W_FROM, W_TO, IDX))
	
	subset(rp, grepl('w200',RUN) & PAIR_ID==5)
	
}	

RakaiCirc.circ.dev160901<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	
	#	get epi info
	tmp		<- RakaiCirc.epi.get.info()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#	get sequence info and add to recipients
	load(file.path(wdir,'RCCS_SeqInfo_160816.rda'))		
	rs		<- subset(rs, !is.na(VISIT))	
	
	
	#	load likely transmissions summary from phyloscanner
	infile.phsc.trms	<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/pty_Rakai_160825/RCCS_run160825_lkltrms_summary.rda"
	load(infile.phsc.trms)
	dlkl				<- copy(df)
	#	load trees from phyloscanner
	infile.phsc.trees	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/pty_Rakai_160825/RCCS_run160825_all_trees.rda'	
	load(infile.phsc.trees)
	phs					<- tmp$phs
	stat.phs			<- tmp$dfr
	
	
	
	#	select likely pairs
	select.discsib	<- 0.65
	dlkl[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]
	tmp		<- dcast.data.table(dlkl, PAIR_ID~TYPE, value.var='WIN_OF_TYPE_P')
	set(tmp, tmp[, which(is.na(disconnected))], 'disconnected', 0)
	set(tmp, tmp[, which(is.na(sib))], 'sib', 0)
	tmp		<- subset(tmp, disconnected+sib < select.discsib, PAIR_ID)
	dlkl	<- merge(dlkl, tmp, by='PAIR_ID')
	#	load female recipients
	rrec	<- RakaiCirc.recipient.female.get.info()
	#	add sequence info to recipients
	rrec	<- merge(rrec, unique(subset(rs, select=c(RID, PID, SID))), by='RID', all.x=1)
	#	select female recipients with at least one PANGEA sequence
	rrec	<- subset(rrec, !is.na(SID))
	#	get likely pairs that involve female recipients
	tmp		<- copy(dlkl)
	setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
	set(tmp, tmp[, which(TYPE=='anc_12')], 'TYPE', 'anc')
	set(tmp, tmp[, which(TYPE=='anc_21')], 'TYPE', 'anc_12')
	set(tmp, tmp[, which(TYPE=='anc')], 'TYPE', 'anc_21')
	tmp		<- rbind(tmp, dlkl)
	setnames(tmp, 'ID1', 'SID')
	rp		<- merge(tmp, rrec, by='SID')
	setnames(rp, 'ID2', 'TR_SID')
	cat('\nFound female recipients among likely pairs, n=', rp[, length(unique(RID))])
	#	get info on transmitters: RID and PID
	tmp		<- subset(rs, !is.na(SID), select=c(RID, PID, SID))
	setkey(tmp, SID)
	tmp		<- unique(tmp)
	setnames(tmp, colnames(tmp), paste('TR_',colnames(tmp),sep=''))	
	#
	#	TODO check whereabouts of 12559_1_9 15065_1_32 15034_1_79 15081_1_71 15430_1_73
	#
	tmp2	<- setdiff( rp[, unique(TR_SID)], tmp[, TR_SID] )
	if(length(tmp2))
		cat('\nWarning: Sanger ID in rp that is not in dictionary', tmp2)
	rp		<- merge(rp, tmp, by='TR_SID')
	#	get info on transmitters: demographic stuff	
	tmp		<- subset(rd, !is.na(PID), select=c(RID, SEX, REGION, COMM_NUM, HH_NUM, BIRTHDATE, LASTNEGVIS, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, RELIGION))
	setkey(tmp, RID)
	tmp		<- unique(tmp)
	setnames(tmp, colnames(tmp), paste('TR_',colnames(tmp),sep=''))
	rp		<- merge(rp,tmp,by='TR_RID')
	#	get info on sequences: date of sampling
	tmp		<- unique(subset(rs, !is.na(SID), select=c(SID, DATE)))
	setnames(tmp, 'DATE', 'SEQ_DATE')
	rp		<- merge(rp, tmp, by='SID')
	setnames(tmp, c('SID','SEQ_DATE'), c('TR_SID','TR_SEQ_DATE'))
	rp		<- merge(rp, tmp, by='TR_SID')
	setkey(rp, SID, TR_SID)
	rpok	<- subset(rp, TR_SEX=='M')	
	rpff	<- subset(rp, TR_SEX=='F')
	rpoku	<- unique(rpok)
	rpffu	<- unique(rpff)
	#	for every F-F pair, keep only one instance
	#	and select different individuals
	rpffu		<- merge(rpffu, rpffu[, list(SID=SID[1],TR_SID=TR_SID[1]), by='PAIR_ID'], by=c('PAIR_ID','SID','TR_SID'))	
	rpffu		<- subset(rpffu, RID!=TR_RID)		
	#
	#	plot evidence
	#		
	tmp		<- copy(rpoku)
	tmp		<- tmp[order(-PAIR_ID),]
	tmp[, LABEL:= factor(PAIR_ID, levels=PAIR_ID, labels=paste('Pair',PAIR_ID, 	'\nPerson1 ', RID, ' ', SID,' -sex:',SEX,' -loc:',REGION,',',COMM_NUM,',',HH_NUM,' -birth:',BIRTHDATE,' -neg:',LASTNEGDATE,' -pos:',FIRSTPOSDATE,' -seq:',SEQ_DATE,
							'\n<->', 
							'\nPerson2 ', TR_RID, ' ', TR_SID,' -sex:',TR_SEX,' -loc:',TR_REGION,',',TR_COMM_NUM,',',TR_HH_NUM,' -birth:',TR_BIRTHDATE,' -neg:',TR_LASTNEGDATE,' -pos:',TR_FIRSTPOSDATE,' -seq:',TR_SEQ_DATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('anc_12','anc_21','sib','int','disconnected'), labels=c('1 ancestral to 2','2 ancestral to 1','1, 2 are siblings','1, 2 are intermingled','1, 2 are disconnected'))])
	ggplot(tmp, aes(x=LABEL, y=WIN_OF_TYPE, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('1 ancestral to 2'="#9E0142",'2 ancestral to 1'="#F46D43",'1, 2 are siblings'="#ABDDA4",'1, 2 are intermingled'="#3288BD",'1, 2 are disconnected'='grey50')) +
			theme_bw() + theme(legend.position='top') +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(wdir, paste('150901_RCCS_lkltransmissionpairs_withfemaleseroconverter_discsiblt',100*select.discsib,'.pdf',sep='')), w=15, h=0.23*nrow(tmp), limitsize = FALSE)	
	#	look at F-F pairs
	tmp			<- copy(rpffu)
	tmp			<- tmp[order(-PAIR_ID),]
	tmp[, LABEL:= factor(PAIR_ID, levels=PAIR_ID, labels=paste('Pair',PAIR_ID, 	'\nPerson1 ', RID, ' ', SID,' -sex:',SEX,' -loc:',REGION,',',COMM_NUM,',',HH_NUM,' -birth:',BIRTHDATE,' -neg:',LASTNEGDATE,' -pos:',FIRSTPOSDATE,' -seq:',SEQ_DATE,
							'\n<->', 
							'\nPerson2 ', TR_RID, ' ', TR_SID,' -sex:',TR_SEX,' -loc:',TR_REGION,',',TR_COMM_NUM,',',TR_HH_NUM,' -birth:',TR_BIRTHDATE,' -neg:',TR_LASTNEGDATE,' -pos:',TR_FIRSTPOSDATE,' -seq:',TR_SEQ_DATE,																				
							'\n',sep=''))]	
	tmp		<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('anc_12','anc_21','sib','int','disconnected'), labels=c('1 ancestral to 2','2 ancestral to 1','1, 2 are siblings','1, 2 are intermingled','1, 2 are disconnected'))])
	ggplot(tmp, aes(x=LABEL, y=WIN_OF_TYPE, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('1 ancestral to 2'="#9E0142",'2 ancestral to 1'="#F46D43",'1, 2 are siblings'="#ABDDA4",'1, 2 are intermingled'="#3288BD",'1, 2 are disconnected'='grey50')) +
			theme_bw() + theme(legend.position='top') +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(wdir, paste('150901_RCCS_lkltransmissionpairs_femalefemale_discsiblt',100*select.discsib,'.pdf',sep='')), w=15, h=0.2*nrow(tmp), limitsize = FALSE)	
	#
	#	re-examine phylogenies
	#	
	setkey(rpok, PAIR_ID, SID, TR_SID)
	rpoku	<- unique(rpok)
	invisible(sapply(seq_len(nrow(rpoku)), function(ii)
					{
						pair.id		<- rpoku[ii, PAIR_ID]
						pty.run		<- rpoku[ii, PTY_RUN]
						dfs			<- subset(dfr, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, TITLE:= dfs[, paste('pair', pair.id, 'ids', rpoku[ii, TR_SID], rpoku[ii, SID], '\nrun', pty.run, '\nwindow', W_FROM,'-', W_TO)]]			
						plot.file	<- file.path(wdir, paste('150901_RCCS_lkltransmissionpairs_withfemaleseroconverter_discsiblt65_pair',pair.id,'.pdf',sep=''))			
						invisible(phsc.plot.selected.pairs(phs, dfs, rpoku[ii, TR_SID], rpoku[ii, SID], plot.file=plot.file, pdf.h=50, pdf.rw=10))
					}))
	
	tmp		<- copy(rpffu)
	invisible(sapply(seq_len(nrow(tmp)), function(ii)
					{
						pair.id		<- tmp[ii, PAIR_ID]
						pty.run		<- tmp[ii, PTY_RUN]
						dfs			<- subset(dfr, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, TITLE:= dfs[, paste('pair', pair.id, 'ids', tmp[ii, TR_SID], tmp[ii, SID], '\nrun', pty.run, '\nwindow', W_FROM,'-', W_TO)]]			
						plot.file	<- file.path(wdir, paste('150901_RCCS_lkltransmissionpairs_femalefemale_discsiblt65_pair',pair.id,'.pdf',sep=''))			
						invisible(phsc.plot.selected.pairs(phs, dfs, tmp[ii, TR_SID], tmp[ii, SID], plot.file=plot.file, pdf.h=50, pdf.rw=10))
					}))
}	
	
RakaiCirc.circ.dev160815<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"	
	tmp		<- RakaiCirc.epi.get.info()
	rh		<- tmp$rh
	rd		<- tmp$rd
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
	
	rrec	<- RakaiCirc.recipient.female.get.info(wdir=NA)		
	#
	#	load info from p24 tree, PANGEA tree, genetic distances
	#
	if(0)
	{
		tmp		<- RakaiCirc.seq.get.phylogenies()
		phf		<- tmp$phf
		php		<- tmp$php
		phfi	<- tmp$phfi
		phpi	<- tmp$phpi		
	}
	load(file=file.path(wdir,'RCCS_PhInfo_160825.rda'))
	#
	#	load info on phylotype runs
	#
	load("~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda")	
	ptyi	<- subset(pty.runs, select=c(TAXA, FILE_ID, PTY_RUN))
	set(ptyi,NULL,'TAXA', ptyi[,gsub('_','-',gsub('_S[0-9]+','',TAXA))])
	setnames(ptyi, c('TAXA','FILE_ID'), c('PID','SID'))
	#
	#	check how many rec in tree, how many prob transmitters, and in which sequ sample
	#
	#pd.max	<- 0.02
	pd.max	<- 0.04
	#	may have multiple seqs per female recipient, keep all
	rrecs	<- merge(rrec, unique(subset(rs, select=c(RID, PID, SEQID, SEQIDb))), by='RID', all.x=1)
	rrecs	<- merge(rrecs, ptyi, by='PID', all.x=1)
	#	add taxa in p24 tree that are close
	#	phfi is only lower triangular, complete 
	tmp		<- subset(phfi, PD<=pd.max)	
	tmp2	<- copy(tmp)	
	setnames(tmp2, c('SEQIDb','SEQIDb2'), c('SEQIDb2','SEQIDb'))	
	tmp		<- rbind(tmp, tmp2, use.names=TRUE)
	set(tmp, NULL, 'SEQIDb', tmp[, as.character(SEQIDb)])
	set(tmp, NULL, 'SEQIDb2', tmp[, as.character(SEQIDb2)])
	#	now add
	rrecs	<- merge(rrecs, tmp, all=1, by='SEQIDb')
	rrecs	<- subset(rrecs, !is.na(RID))
	#	add basic epi info and IDs for transmitters
	setnames(rrecs, 'SEQIDb2', 'TR_SEQIDb') 
	tmp		<- unique(subset(rs, !is.na(SEQIDb), select=c(RID, PID, SEQIDb, SEQTYPE)))
	tmp2	<- subset(rd, select=c(RID, SEX))
	setkey(tmp2, RID)
	tmp		<- merge(tmp, unique(tmp2), by='RID')
	#	add basic info about phylotype runs for transmitters	
	tmp		<- merge(tmp, ptyi, by='PID', all.x=1)
	setnames(tmp, c('RID', 'PID', 'SEQIDb', 'SID', 'PTY_RUN','SEQTYPE', 'SEX'), c('TR_RID', 'TR_PID', 'TR_SEQIDb', 'TR_SID', 'TR_PTY_RUN', 'TR_SEQTYPE', 'TR_SEX'))	
	rp		<- merge(rrecs, tmp, by='TR_SEQIDb', all.x=1)	
	tmp		<- subset(rp, PTY_RUN!=TR_PTY_RUN & TR_SEX=='M')	
	if(nrow(tmp))
	{
		cat('\nWarning: potential M->F transmission pair not in phylotype groupings\n', paste( tmp[, paste(TR_PID, '->', PID, sep='')], collapse='\n'))
		save(tmp, file=file.path(wdir, paste('150830_FastTree_gagAllSmall_pairs_patristic',pd.max*100,'_missedpotentialpairs.rda',sep='')))
	}
	#	plot tree with pairs highlighted
	if(0)
	{
		tmp		<- data.table(SEQIDb=phf$tip.label, IDX=seq_len(Ntip(phf)), TIP_COL='black')
		tmp2	<- data.table(SEQIDb=subset(rp, !is.na(TR_SEQIDb))[, unique(na.omit(c(SEQIDb, TR_SEQIDb))) ], PAIR='Y')	
		tmp		<- merge(tmp, tmp2, by='SEQIDb', all.x=1)
		set(tmp, tmp[, which(PAIR=='Y')], 'TIP_COL', 'red')
		setkey(tmp, IDX)
		file	<- file.path(wdir, paste('150825_FastTree_gagAllSmall_pairs_patristic',pd.max*100,'pc.pdf',sep=''))
		invisible(hivc.clu.plot(phf, tip.color=tmp[,TIP_COL], file=file, pdf.scaley=100, show.tip.label=TRUE, pdf.width=30))		
	}
	#
	#	get sample counts
	#
	#	for every female recipient, count prob transmitters by SEQ_TYPE
	tmp		<- subset(rp, SC_WINDOW<=2)
	rpi		<- tmp[, list(	FIRSTPOSDATE=			FIRSTPOSDATE[1],
							SEQ_R_N= 				ifelse(any(!is.na(SEQID)), length(unique(SEQID)), 0L),
							SEQ_R_PANGEA= 			any(grepl('PANGEA', SEQ_TYPE)),
							PHP24_R_N= 				ifelse(any(!is.na(SEQIDb)), length(unique(SEQIDb)), 0L),
							PHP24_TR_N= 			ifelse(any(!is.na(TR_SEQIDb)), length(unique(TR_SEQIDb)), 0L),
							PHP24_MTR_N= 			ifelse(any(!is.na(TR_SEQIDb[TR_SEX=='M'])), length(unique(TR_SEQIDb[TR_SEX=='M'])), 0L),
							PHP24_MTR_MINPD=		ifelse(any(!is.na(PD[TR_SEX=='M'])), min(PD[TR_SEX=='M'], na.rm=TRUE), NA_real_),
							PHP24_MTR_PANGEA= 		length(which(grepl('PANGEA', TR_SEQTYPE[TR_SEX=='M']))),
							PHP24_MTR_PANGEA_MINPD=	ifelse(any(!is.na(PD[TR_SEX=='M' & TR_SEQTYPE=='PANGEA'])), min(PD[TR_SEX=='M' & TR_SEQTYPE=='PANGEA'], na.rm=TRUE), NA_real_)
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
	#	save
	save(rp, rpi, file=file.path(wdir, paste('150825_FastTree_gagAllSmall_pairs_patristic',pd.max*100,'pc.rda',sep='')))	

	 #	plot
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

	#	extract phylotype candidates
	pty.rp	<- subset(rp, !is.na(PID) & !is.na(TR_PID) & !is.na(TR_SID) & !is.na(SID) & TR_SEX=='M')
	#	add geographic info on transmitters
	#	TODO have these been moving in the RCCS?
	tmp2	<- subset(rd, !is.na(PID) & SEX=='M', select=c(RID, REGION, COMM_NUM, HH_NUM, BIRTHDATE))
	setkey(tmp2, RID)
	tmp2	<- unique(tmp2)
	setnames(tmp2, c('RID','REGION','COMM_NUM','HH_NUM','BIRTHDATE'), c('TR_RID','TR_REGION','TR_COMM_NUM','TR_HH_NUM','TR_BIRTHDATE'))
	pty.rp	<- merge(pty.rp, tmp2, by='TR_RID')	
	setkey(pty.rp, PTY_RUN, RID)
	write.csv(pty.rp, row.names=FALSE, file=file.path(wdir, paste('150830_FastTree_Pairs_SCw2_patristic',pd.max*100,'_PANGEA_M2F_hypothetical_pairs.csv',sep='')))
	save(pty.rp, file=file.path(wdir, paste('150830_FastTree_Pairs_SCw2_patristic',pd.max*100,'_PANGEA_M2F_hypothetical_pairs.rda',sep='')))
	#
	#	TODO add circumcision status on transmitters (in rt --> now move to timelines)
	#	TODO add ART trajectories
	#
	
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
######################################################################################
hivc.db.Date2numeric<- function( x )
{
	if(!class(x)%in%c('Date','character'))	return( x )
	x	<- as.POSIXlt(x)
	tmp	<- x$year + 1900
	x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
	x	
}
######################################################################################
######################################################################################
project.Rakai.aliRegion1.597<- function()
{	
	require(big.phylo)
	require(plyr)
	
	#	load 597 new files
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/PANGEA_UG full alignment/PANGEA_HIV_n597_Imperial_v160916_RakaiPlates.fasta'
	sqn		<- read.dna(infile, format='fa')
	sqni	<- data.table(TAXA=rownames(sqn))
	sqni	<- subset(sqni, grepl('HXB2|consensus',TAXA))
	sqni[, SID:= gsub('_consensus','',TAXA)]
	set(sqni, sqni[, which(grepl('HXB2',TAXA))],'SID',NA_character_)
	#merge(sqni, unique(subset(rs, !is.na(SID), select=c(RID,SEQTYPE,SID,PID))), by='SID',all.x=1)
	sqn		<- sqn[ sqni[,TAXA], ]
	
	#	load PANGEA alignment
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_HXB2.fasta'
	sqp		<- read.dna(infile,format='fa')
	
	#	required: HXB2 in alignment
	ans		<- seq.align.based.on.common.reference(sqn, sqp, return.common.sites=TRUE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2.fasta'
	write.dna(ans, file=outfile, format='fasta', colsep='', nbcol=-1)
	
	#	get genetic distance matrix
	sq			<- ans
	sq			<- as.character(sq)
	sq[sq=='?']	<- '-'
	sq			<- as.DNAbin(sq)
	sq.gd		<- dist.dna(sq, pairwise.deletion=TRUE)
	sq.gd		<- as.data.table(melt(as.matrix(sq.gd), varnames=c('TAXA1','TAXA2')))
	setnames(sq.gd, 'value', 'PD')
	sq.gd		<- subset(sq.gd, TAXA1!=TAXA2)
	set(sq.gd, NULL, 'TAXA1', sq.gd[, as.character(TAXA1)])
	set(sq.gd, NULL, 'TAXA2', sq.gd[, as.character(TAXA2)])		
	#	add overlap	
	setkey(sq.gd, TAXA1, TAXA2)
	tmp			<- as.character(sq)
	tmp			<- !( tmp=='-' | tmp=='?' | tmp=='n' )
	tmp[]		<- as.integer(tmp) 
	tmp2		<- sq.gd[, list(OVERLAP= sum(bitwAnd(tmp[TAXA1,], tmp[TAXA2,])) ), by=c('TAXA1','TAXA2')]	
	sq.gd		<- merge(sq.gd, tmp2, by=c('TAXA1','TAXA2'))
	save(sq, sq.gd, file=gsub('\\.fasta','\\.rda',outfile))
	
	#
	#	resolve SANGER IDs
	#
	load('~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2.rda')
	sqi		<- data.table(TAXA=rownames(sq), ID=seq_len(nrow(sq)))
	
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_data/2016-09-18_PAN_SANGER_IDs.txt'
	pni		<- as.data.table(read.table(infile, header=TRUE, sep='\t', comment.char='%'))
	setnames(pni, colnames(pni), toupper(colnames(pni)))
	set(pni, NULL, 'SANGERID', pni[, gsub('#','_',SANGERID)])
	pni[, TAXA:= paste(SANGERID,'_consensus',sep='')]
	sqi		<- merge(sqi, pni, by='TAXA',all.x=1)
	stopifnot( nrow(subset(sqi, grepl('consensus',TAXA) & is.na(SANGERID)))==0 )
	#	set new taxa names
	sqi[, TAXA_NEW:=TAXA]
	tmp		<- sqi[, which(!is.na(PANGEA_ID))]
	set(sqi, tmp, 'TAXA_NEW', sqi[tmp, PANGEA_ID])
	setkey(sqi, ID)
	rownames(sq)	<- sqi[, TAXA_NEW]
	setnames(sqi, c('TAXA','TAXA_NEW'), c('TAXA1','TAXA1_NEW'))
	sq.gd			<- merge(sq.gd, subset(sqi, select=c(TAXA1,TAXA1_NEW)), by='TAXA1')
	setnames(sqi, c('TAXA1','TAXA1_NEW'), c('TAXA2','TAXA2_NEW'))
	sq.gd			<- merge(sq.gd, subset(sqi, select=c(TAXA2,TAXA2_NEW)), by='TAXA2')
	set(sq.gd, NULL, c('TAXA1','TAXA2'), NULL)
	setnames(sq.gd, c('TAXA1_NEW','TAXA2_NEW'), c('TAXA1','TAXA2'))
	#	write to file
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2.fasta'
	write.dna(sq, file=outfile, format='fasta', colsep='', nbcol=-1)	
	#
	#	read COMETv0.5
	#
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2_COMETv0.5.txt'
	sqi		<- as.data.table(read.table(infile, header=TRUE, sep='\t',stringsAsFactors=FALSE))
	setnames(sqi, c('name','subtype','length'), c('TAXA','COMET_ST','COMET_N'))
	sqi[, COMET_V:='0.5']
	#
	#	read COMETv2.1
	#
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2_COMETv2.1.txt'
	tmp		<- as.data.table(read.table(infile, header=TRUE, sep='\t',stringsAsFactors=FALSE))
	setnames(tmp, c('name','subtype','bootstrap.support'), c('TAXA','COMET_ST','COMET_BOOTSTRAP'))
	tmp[, COMET_V:='2.1']
	sqi		<- rbind(sqi, tmp, use.names=TRUE, fill=TRUE)
	#	set to factors to save a bit of space
	set(sq.gd, NULL, 'TAXA1', sq.gd[,factor(TAXA1)])
	set(sq.gd, NULL, 'TAXA2', sq.gd[,factor(TAXA2)])
	save(sq, sqi, sq.gd, file=gsub('\\.fasta','\\.rda',outfile))
}
######################################################################################
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