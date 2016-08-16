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

RakaiCirc.prepare.seqinfo<- function()
{
	infile				<- "~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda"
	indir.historicseq	<- "~/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree"
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	#
	#	prepare RCCS data (need towards end)
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
	#	read all processed RCCS sequences 
	#	
	rs		<- data.table(FILE=list.files(indir.historicseq, full.names=TRUE, pattern='Prelim_RakaiPangeaSqnAndMetaData_IqTree_Region[0-9].rda'))
	rs[, GENE_REGION:= regmatches(FILE, regexpr('Region[0-9]', FILE))]		
	tmp		<- lapply(subset(rs, GENE_REGION%in%c('Region1','Region2','Region3'))[, FILE], function(z){
				load(z)
				z2	<- subset(as.data.table(pangeaMetaData), select=c(RCCS_studyid, visit, date, seqid, CometSubtype, Pangea.id))
				z2[, SEQTYPE:='PANGEA']
				z2[, GENE_REGION:= regmatches(z, regexpr('Region[0-9]', z))]
				z2
			})
	tmp2	<- lapply(subset(rs, GENE_REGION%in%c('Region1','Region4'))[, FILE], function(z){
				load(z)
				z2	<- subset(as.data.table(rakhistMetaData), select=c(RCCS_studyid, visit, date, seqid, CometSubtype))
				z2[, SEQTYPE:='HISTORIC']
				z2[, GENE_REGION:= regmatches(z, regexpr('Region[0-9]', z))]
				z2
			})
	rs		<- do.call('rbind',tmp) 
	tmp2	<- do.call('rbind',tmp2)
	rs		<- rbind(rs, tmp2, use.names=TRUE, fill=TRUE)	
	setnames(rs, colnames(rs), gsub('\\.','_',toupper(colnames(rs))))
	setnames(rs, c('RCCS_STUDYID','PANGEA_ID'), c('RID','PID'))
	set(rs, NULL, 'DATE', hivc.db.Date2numeric(rs[['DATE']]))
	rs		<- subset(rs, !is.na(RID))		
	#	add unprocessed but shipped PANGEA sequences
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
	#	save to file
	#
	save(rs, file=file.path(wdir,'RCCS_SeqInfo_160816.rda'))
}

RakaiCirc.dev160815<- function()
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
	#
	#	add first sequences
	#	
	#	load sequence info
	#	(if not exist, run RakaiCirc.prepare.seqinfo)
	load(file.path(wdir,'RCCS_SeqInfo_160816.rda'))	
	#	TODO: has not visit: E21593M2
	rs		<- subset(rs, !is.na(VISIT))	
	tmp		<- rs[, {
				z<- which.min(DATE)
				list(VISIT=VISIT[z], DATE=DATE[z], SEQ_HIST= any(SEQTYPE=='HISTORIC'), SEQ_PANGEA_NOW= any(SEQTYPE=='PANGEA'), SEQ_PANGEA_FUTURE= any(SEQTYPE=='PANGEA_not_yet_processed'), SEQ_NOW_GAG= any(GENE_REGION=='Region1'))
			}, by=RID]
	#	complete to all available RIDs
	tmp		<- merge(unique(subset(rt, select=RID)), tmp, all.x=1)
	tmp[, SEQ_TYPE:= 'no_seq']
	set(tmp, tmp[, which(SEQ_HIST & !SEQ_PANGEA_NOW & !SEQ_PANGEA_FUTURE)], 'SEQ_TYPE', "hist_only")
	set(tmp, tmp[, which(SEQ_PANGEA_NOW | SEQ_PANGEA_FUTURE)], 'SEQ_TYPE', "pangea")
	set(tmp, tmp[, which(SEQ_HIST & (SEQ_PANGEA_NOW | SEQ_PANGEA_FUTURE))], 'SEQ_TYPE', "hist_pangea")
	tmp[, SEQ_NOW:= factor(!is.na(SEQ_HIST) & (SEQ_HIST | SEQ_PANGEA_NOW), levels=c(TRUE,FALSE), labels=c('Y','N'))]
	set(tmp, NULL, 'SEQ_NOW_GAG', tmp[, factor(!is.na(SEQ_NOW_GAG) & SEQ_NOW_GAG, levels=c(TRUE,FALSE), labels=c('Y','N'))])
	set(tmp, NULL, c('VISIT','SEQ_HIST','SEQ_PANGEA_NOW','SEQ_PANGEA_FUTURE'), NULL)
	setnames(tmp, 'DATE', 'FIRST_SEQ_DATE')
	rt		<- merge(rt, tmp, by='RID')
	
	
	#	plot total circumcision by year among individuals in follow up
	ggplot(rt, aes(x=factor(ROUND_ST), fill=CIRC)) + geom_bar() + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1e4, 5e2)) +
			facet_grid(SEX~.) +
			labs(x='\nsurvey rounds\n(start date)', y='HIV infected RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Circumcision\nstatus')
	ggsave(file.path(wdir, 'RCCSinf_by_circ.pdf'), w=10, h=8)
	#	plot proportion circumcision by year among individuals in follow up
	ggplot(subset(rt, SEX=='male'), aes(x=factor(ROUND_ST), fill=CIRC)) + geom_bar(position='fill') + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1, 0.1), labels = scales::percent) +
			labs(x='\nsurvey rounds\n(start date)', y='HIV infected male RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Circumcision\nstatus')
	ggsave(file.path(wdir, 'RCCSinf_by_circ_proportions_men.pdf'), w=10, h=5)	
	#	plot sequence coverage among HIV infected participants
	tmp		<- copy(rt)
	tmp[, DUMMY:= as.integer(SEQ_TYPE!='no_seq')]
	set(tmp, tmp[, which(SEQ_NOW_GAG=='Y')], 'DUMMY', 2L)
	set(tmp, NULL, 'DUMMY', tmp[, as.character(factor(DUMMY, levels=c(0L,1L, 2L), labels=c('No sequence','sequence available','sequence forthcoming')))])
	ggplot(tmp, aes(x=factor(ROUND_ST), fill=DUMMY)) + geom_bar() + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1e4, 5e2)) +
			labs(x='\nsurvey rounds\n(start date)', y='HIV infected RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Sequence\nstatus')
	ggsave(file.path(wdir, 'RCCSinf_by_sequence_coverage.pdf'), w=10, h=5)
	#	plot circumcised by religion
	tmp		<- subset(rt, EVERCIRC=='Y')
	ggplot(tmp, aes(x=factor(ROUND_ST), fill=RELIGION)) + geom_bar(position='fill') + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1, 0.1), labels = scales::percent) +
			scale_fill_brewer(palette='Set2') +
			labs(x='\nsurvey rounds\n(start date)', y='circumcised male \nHIV infected RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Faith')
	ggsave(file.path(wdir, 'RCCScircmen_by_religion_proportions.pdf'), w=10, h=5)	
	
}


######################################################################################
project.Rakai.processRegion1PANGEAalignment.151129<- function()
{
	susa.f	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151129.fasta'
	susa.s	<- read.dna(susa.f, format='fasta')
	susa.d	<- data.table(TAXA=rownames(susa.s), DATA= factor(grepl('^PG', rownames(susa.s)), levels=c(TRUE,FALSE), label=c('PNG','LNL')))
	
	png.f	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/PANGEA_orig/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	load(png.f)	#sq, sqi, si
	
	tmp		<- subset(sqi, is.na(SITE))[, PANGEA_ID]
	susa.d[, tmp%in%TAXA]	#Susanna did not keep the references?
	tmp		<- sapply(strsplit(tmp, '.', fixed=1), function(x) tail(x,1))
	sapply(tmp, function(x) susa.d[, any(grepl(x, TAXA))]	)	#no...
	
	#
	# find PANGEA sequence with greatest coverage in the PANGEA and Su alignment overlap
	#
	sqiu	<- subset(sqi, SITE=='UG')	
	tmp		<- subset(sqiu, COV==max(COV))[1, TAXA]
	susa.s[tmp,]	# empty!	
	#	3628 UG seqs in PANGEA alignment
	subset(susa.d, DATA=='PNG')
	#	2923 UG PANGEA in SU alignment	
	tmp		<- merge(susa.d, sqiu, by='TAXA')
	#	2903 of 3628 in SU alignment
	setdiff( subset(susa.d, DATA=='PNG')[, TAXA], sqiu[, TAXA])	# Susanna gave a different name to the duplicates TODO: Chris to change taxa names with ( and space
	setdiff( sqiu[, TAXA], subset(susa.d, DATA=='PNG')[, TAXA])	# quite a few more missing!
	#	find start of SU alignment in PANGEA alignment
	tmp2	<- susa.s[tmp[, TAXA],]
	tmp2	<- as.character(tmp2[,1:50])
	tmp2	<- data.table(TAXA=rownames(tmp2), COV50=apply(tmp2, 1, function(x) length(which(!x%in%c('-','?')))	))	
	mxcv	<- tmp2[which.max(COV50),TAXA]						# mxcov for key
	tmp2	<- gsub('?','',gsub('-','',paste(as.vector(as.character(susa.s[mxcv,1:50])), collapse='')),fixed=1)		#build SU alignment key
	tmp2	<- paste(strsplit(tmp2, '')[[1]],'-*',sep='',collapse='')
	startat	<- regexpr(tmp2, gsub('?','-',paste(as.vector(as.character(sq[mxcv,])), collapse=''),fixed=1))	
	#	TODO: double check that start of region1 is at PANGEA 1
	#	TODO: startat ignores - at start, so pos will not be right if startat != 1
	#	code below assumes that startat is 1
	
	#
	# calculate offset
	#	
	lanl.s	<- susa.s[subset(susa.d, DATA!='PNG')[, TAXA],]
	tmp2	<- as.character(lanl.s)
	tmp2	<- data.table(	TAXA= rownames(lanl.s), 
			FIRST= apply( tmp2, 1, function(x) which(x!='-')[1] ),
			LAST= ncol(lanl.s)-apply( tmp2, 1, function(x) which(rev(x)!='-')[1] )+1L		)
	lanl.d	<- merge(susa.d, tmp2, by='TAXA')	
	lanl.l	<- lanl.d[, max(LAST)]			#this is the final pos with partial seqs in the new alignment
	sqm		<- sq[tmp[, TAXA],seq.int(startat,lanl.l)]
	tmp2	<- as.character(sqm)
	tmp2	<- data.table(TAXA=rownames(tmp2), COV=apply(tmp2, 1, function(x) length(which(!x%in%c('-','?')))	))
	mxcv	<- tmp2[which.max(COV),TAXA]						# mxcov for offset calculations
	# calculate offset	
	y		<- susa.s[mxcv,]		#y must be the longer sequence with extra gaps
	x		<- sq[mxcv,]	
	x		<- as.character(x)	
	y		<- as.character(y)	
	stopat	<- 9e3
	#substring(paste(as.vector(x), collapse=''),280,400)
	#k<- 285
	#rbind( z[1, ((k-5):(k+20))+offset.in.x[(k-5):(k+20)]], z[2, ((k-5):(k+20))+offset.in.y[(k-5):(k+20)]] )
	#z[1, ((k-5):(k+20))]
	if(0)	#DEV
	{
		x<- matrix(strsplit('tcaa---------aaattt','')[[1]], nrow=1)
		y<- matrix(strsplit('-tcaaaa---a-ttt----','')[[1]], nrow=1)		
	}	
	stopifnot( gsub('?','',gsub('-','',paste(as.vector(x), collapse='')),fixed=1)==gsub('?','',gsub('-','',paste(as.vector(y), collapse='')),fixed=1) )
	z							<- rbind.fill.matrix(x,y)
	z[1,which(is.na(z[1,]))]	<- '-'
	z[2,which(is.na(z[2,]))]	<- '-'
	zz							<- unname(z)	
	offset.to.x	<- rep(0, ncol(z))
	offset.to.y	<- rep(0, ncol(z))			
	k			<- seq_len(ncol(z))+offset.to.y
	k			<- k[ k>0 & k<=ncol(z)]
	k			<- which( z[1, seq_along(k)]!=z[2, k] )[1]
	while(!is.na(k) & k<stopat)
	{
		#print(k)		
		done	<- 0
		if( !zz[1,k]%in%c('-','?') )
		{
			kn		<- k-offset.to.x[min(k,length(offset.to.x))]
			offset.to.x[ seq.int(kn,length(offset.to.x)) ] <- offset.to.x[ seq.int(kn,length(offset.to.x)) ]+1
			if(k==1)
				zz	<- rbind( c('-',zz[1,]), c(zz[2,],'-') )
			if(k>1)
				zz	<- rbind( c(zz[1, 1:(k-1)],'-',zz[1,k:ncol(zz)]), c(zz[2,],'-')	)
			done	<- 1
		}
		if( !done & zz[1,k]%in%c('-','?') )
		{
			kn		<- k-offset.to.y[min(k,length(offset.to.x))]
			offset.to.y[ seq.int(kn,length(offset.to.x)) ] <- offset.to.y[ seq.int(kn,length(offset.to.y)) ]+1
			if(k==1)
				zz	<- rbind( c(zz[1,],'-'), c('-',zz[2,]) )
			if(k>1)
				zz	<- rbind( c(zz[1,],'-'), c(zz[2, 1:(k-1)],'-',zz[2,k:ncol(zz)])	)							
		}
		k	<- unname(which(zz[1,]!=zz[2,])[1])		
	}
	if(0)
	{
		#	DEV
		# to x add offset.to.x
		k			<- offset.to.x+seq_along(offset.to.x)
		xx			<- matrix('-',ncol=max(k),nrow=nrow(x))
		xx[,k]		<- x[, seq_along(k)]
		# to y add offset.to.y
		k			<- offset.to.y+seq_along(offset.to.y)
		yy			<- matrix('-',ncol=max(k),nrow=nrow(y))
		yy[,k]		<- y[, seq_along(k)]
		rbind(xx, yy)
	}
	# add gaps to SU alignment
	tmp			<- as.character(lanl.s[, seq_len(lanl.l)])
	lanl.s2		<- rbind(tmp,as.character(susa.s[mxcv, seq_len(lanl.l)]))	#ONLY TO CHECK		
	k			<- offset.to.y[seq_len(lanl.l)]+seq_len(lanl.l)
	tmp			<- matrix('-',ncol=max(k),nrow=nrow(lanl.s2), dimnames=dimnames(lanl.s2))
	tmp[,k]		<- lanl.s2[, seq_along(k)]
	lanl.s2		<- as.DNAbin(tmp)
	# add gaps to PNGEA alignment 
	tmp			<- offset.to.x+seq_along(offset.to.x)
	k			<- tmp[tmp<=ncol(lanl.s2)]
	tmp			<- as.character(sq[,seq_along(k)])			#take PANGEA alignment and stretch
	sqn			<- matrix('-',nrow=nrow(tmp), ncol=max(k), dimnames=dimnames(tmp))
	sqn[,k]		<- tmp[,seq_along(k)]							
	sqn			<- as.DNAbin(sqn)
	# now add non-stretched part of PANGEA alignment
	sqn			<- cbind(sqn,sq[,seq.int(length(k)+1,ncol(sq))])
	# now bring LANL to same size as new alignment
	tmp			<- as.DNAbin( matrix('-',nrow=nrow(lanl.s2), ncol=ncol(sqn)-ncol(lanl.s2), dimnames=dimnames(lanl.s2)) )
	lanl.s2		<- cbind(lanl.s2, tmp)
	# now merge!
	sqn			<- rbind(sqn,lanl.s2)
	if(0)	
	{	
		#CHECK
		write.dna( sqn[which(mxcv==rownames(sqn)),], format='fasta', file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/check.fasta' )
		# regmatches(paste(as.vector(x), collapse=''),regexpr('ttttttagggagactt.*',paste(as.vector(x), collapse='')))		
	}
	#	rm Check
	sqn			<- sqn[-which(mxcv==rownames(sqn))[2],]
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