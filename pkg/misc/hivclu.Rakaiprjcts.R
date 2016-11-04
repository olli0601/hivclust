RakaiCouples.setup.phyloscan.runs<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ape)
	
	wdir	<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples"
	#
	#	load table of sanger ids and pangea ids
	load("~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda")
	#	we need from load: dm,  sqi
	#
	#	load couples data
	infile	<- '~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/Pangea_Couples.csv'
	rc		<- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	setnames(rc, c('female.PANGEA.ID','male.PANGEA.ID'), c('female.TAXA','male.TAXA'))
	setnames(rc, colnames(rc), gsub('\\.','_',toupper(colnames(rc))))
	setnames(rc, c('MALE_RCCS_STUDYID','FEMALE_RCCS_STUDYID'), c('MALE_RID','FEMALE_RID'))
	set(rc, NULL, 'MALE_DATE', rc[, hivc.db.Date2numeric(as.Date(MALE_DATE))])
	set(rc, NULL, 'FEMALE_DATE', rc[, hivc.db.Date2numeric(as.Date(FEMALE_DATE))])
	set(rc, NULL, 'MALE_LASTNEGDATE', rc[, hivc.db.Date2numeric(as.Date(MALE_LASTNEGDATE))])
	set(rc, NULL, 'FEMALE_LASTNEGDATE', rc[, hivc.db.Date2numeric(as.Date(FEMALE_LASTNEGDATE))])	
	set(rc, NULL, 'MALE_FIRSTPOSDATE', rc[, hivc.db.Date2numeric(as.Date(MALE_FIRSTPOSDATE))])
	set(rc, NULL, 'FEMALE_FIRSTPOSDATE', rc[, hivc.db.Date2numeric(as.Date(FEMALE_FIRSTPOSDATE))])
	#
	#	see if couples in same phylotype runs that were already set up
	infile	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda'
	load(infile)
	#	prepare patristic distance matrix
	ph.gdtr	<- as.data.table(melt(ph.gdtr, varnames=c('TAXA1','TAXA2')))
	setnames(ph.gdtr, 'value', 'PD')
	ph.gdtr	<- subset(ph.gdtr, TAXA1!=TAXA2)
	set(ph.gdtr, NULL, 'TAXA1', ph.gdtr[, gsub('_','-',as.character(TAXA1))])
	set(ph.gdtr, NULL, 'TAXA2', ph.gdtr[, gsub('_','-',as.character(TAXA2))])
	tmp		<- subset(pty.runs, select=c(PTY_RUN, TAXA))
	setnames(tmp, c('TAXA','PTY_RUN'), c('TAXA2','PTY_RUN'))
	ph.gdtr	<- merge(ph.gdtr, tmp, all.x=1, by='TAXA2')
	#	load genetic distance matrix with overlap
	infile		<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_gd.rda'
	load(infile)	#loads sq.gd
	#	reduce genetic distances to Rakai seqs
	tmp			<- subset(dm, COHORT=='RCCS', TAXA) 
	setnames(tmp, 'TAXA', 'TAXA1')
	sq.gd		<- merge(sq.gd, tmp, by='TAXA1')
	setnames(tmp, 'TAXA1', 'TAXA2')
	sq.gd		<- merge(sq.gd, tmp, by='TAXA2')
	#	add PTY_RUN
	tmp			<- subset(pty.runs, select=c(PTY_RUN, TAXA))
	setnames(tmp, c('TAXA','PTY_RUN'), c('TAXA2','PTY_RUN'))
	set(tmp, NULL, 'TAXA2', tmp[, gsub('_','-',as.character(TAXA2))])
	sq.gd		<- merge(sq.gd, tmp, all.x=1, by='TAXA2')
	#
	#
	#
	tmp		<- copy(pty.runs)
	setnames(tmp, c('FILE_ID','PTY_RUN','IDX'), c('FEMALE_SANGER_ID','FEMALE_PTY_RUN','FEMALE_PHIDX'))
	tmp		<- subset(tmp, select=c('FEMALE_SANGER_ID','FEMALE_PTY_RUN','FEMALE_PHIDX'))
	setkey(tmp, FEMALE_SANGER_ID)
	rc		<- merge(rc, unique(tmp), by='FEMALE_SANGER_ID', all.x=1)
	setnames(tmp, colnames(tmp), gsub('^FE','',colnames(tmp)))
	rc		<- merge(rc, unique(tmp), by='MALE_SANGER_ID', all.x=1)
	#	get all possible unique pairings of SANGER_IDs from couples
	rp		<- subset(rc, !is.na(MALE_SANGER_ID), c(COUPID, MALE_RID, MALE_SANGER_ID, MALE_TAXA, MALE_PTY_RUN ))
	tmp		<- subset(rc, !is.na(FEMALE_SANGER_ID), c(COUPID, FEMALE_RID, FEMALE_SANGER_ID, FEMALE_TAXA, FEMALE_PTY_RUN ))	
	rp		<- merge(rp, tmp, by='COUPID')
	setkey(rp, FEMALE_SANGER_ID, MALE_SANGER_ID)
	rp		<- unique(rp)
	#	exclude couple F108560:F108560
	rp		<- subset(rp, FEMALE_SANGER_ID!=MALE_SANGER_ID)
	#
	cat('\npairings that were shipped', nrow(rp))
	#
	#	drop pairings including SANGER_IDs for which I don t have bam files 
	rp		<- merge(rp, data.table(MALE_TAXA=sqi[, TAXA]), by='MALE_TAXA')	
	rp		<- merge(rp, data.table(FEMALE_TAXA=sqi[, TAXA]), by='FEMALE_TAXA')
	#
	cat('\npairings for which I have bam files right now', nrow(rp))
	#
	#	assing run id to partners with no run id
	tmp		<- rp[, which(!is.na(MALE_PTY_RUN) & is.na(FEMALE_PTY_RUN))]
	set(rp, tmp, 'FEMALE_PTY_RUN', rp[tmp, MALE_PTY_RUN])
	tmp		<- rp[, which(is.na(MALE_PTY_RUN) & !is.na(FEMALE_PTY_RUN))]
	set(rp, tmp, 'MALE_PTY_RUN', rp[tmp, FEMALE_PTY_RUN])
	#	fixup inconsistent phylotype run assignments
	tmp		<- subset(rp, FEMALE_PTY_RUN!=MALE_PTY_RUN)
	set(tmp, NULL, 'FEMALE_PTY_RUN', tmp[, MALE_PTY_RUN] )
	tmp2	<- rp[, which(FEMALE_PTY_RUN!=MALE_PTY_RUN)]
	set(rp, tmp2, 'MALE_PTY_RUN', rp[tmp2, FEMALE_PTY_RUN])
	rp		<- rbind(rp, tmp)
	#	assign run id to pairs with no run id; discard pairs with no reasonable overlap to anything else
	tmp		<- subset(rp, is.na(MALE_PTY_RUN))
	tmp2	<- subset(sq.gd, is.finite(PD))
	tmp		<- tmp[, {
				z	<- subset(tmp2, (TAXA1==MALE_TAXA|TAXA1==FEMALE_TAXA) & OVERLAP>400 & !is.na(PTY_RUN))
				if(nrow(z)>0)
					ans	<- z[which.min(PD), PTY_RUN]				
				if(nrow(z)==0 & nrow(subset(tmp2, (TAXA1==MALE_TAXA|TAXA1==FEMALE_TAXA) & OVERLAP>400)))
					ans	<- 0L
				if(nrow(z)==0)
					ans	<- -1L
				list(PTY_RUN= ans)
			}, by=c('COUPID','MALE_TAXA','FEMALE_TAXA')]
	rp		<- merge(rp, subset(tmp, PTY_RUN>0), by=c('COUPID','MALE_TAXA','FEMALE_TAXA'), all.x=1)
	tmp		<- rp[, which(is.na(MALE_PTY_RUN) & !is.na(PTY_RUN))]
	set(rp, tmp, c('FEMALE_PTY_RUN','MALE_PTY_RUN'), rp[tmp,PTY_RUN])
	set(rp, NULL, 'PTY_RUN', NULL)
	cat('\nnot enough sequence overlap to find good taxa (ie couple cons sequence too short)',subset(rp, is.na(MALE_PTY_RUN))[, paste(COUPID, collapse=', ')])
	rp		<- subset(rp, !is.na(MALE_PTY_RUN))
	setkey(rp, COUPID, FEMALE_TAXA, MALE_TAXA, MALE_PTY_RUN)
	rp		<- unique(rp)	
	#	make sure that multiple sequences from same male are never in same phylotype run
	tmp		<- rp[, list(N=length(MALE_SANGER_ID)), by=c('MALE_RID','MALE_PTY_RUN')]
	tmp		<- subset(tmp, N>1)
	tmp		<- merge(tmp, rp, by=c('MALE_RID','MALE_PTY_RUN'))
	tmp		<- tmp[, {
				#MALE_TAXA<- 'PG14-UG501971-S03975'; MALE_PTY_RUN<- c(64)
				cat('\n',MALE_RID)
				z	<- subset(tmp2, TAXA1%in%MALE_TAXA & OVERLAP>400 & !is.na(PTY_RUN) & !PTY_RUN%in%unique(MALE_PTY_RUN))
				if(nrow(z)==0)
					z	<- subset(tmp2, TAXA1%in%MALE_TAXA & OVERLAP>200 & !is.na(PTY_RUN) & !PTY_RUN%in%unique(MALE_PTY_RUN))
				if(nrow(z)==0)
					z	<- subset(tmp2, TAXA1%in%MALE_TAXA & OVERLAP>150 & !is.na(PTY_RUN) & !PTY_RUN%in%unique(MALE_PTY_RUN))				
				stopifnot(nrow(z)>0)
				z	<- z[, list(PD=min(PD)), by='PTY_RUN']
				setkey(z, PD)
				ans						<- MALE_PTY_RUN
				ans[duplicated(ans)]	<- z[seq_len( length(MALE_PTY_RUN)-length(unique(ans)) ), PTY_RUN]
				list(MALE_TAXA=MALE_TAXA, FEMALE_TAXA=FEMALE_TAXA, MALE_PTY_RUN=MALE_PTY_RUN, PTY_RUN= ans)
			}, by=c('MALE_RID')]
	rp		<- merge(rp, tmp, by=c('MALE_RID','MALE_TAXA','FEMALE_TAXA','MALE_PTY_RUN'), all.x=1)
	tmp		<- rp[, which(!is.na(PTY_RUN))]
	set(rp, tmp, c('FEMALE_PTY_RUN','MALE_PTY_RUN'), rp[tmp,PTY_RUN])
	set(rp, NULL, 'PTY_RUN', NULL)
	#	make sure that multiple sequences from same female are never in same phylotype run
	tmp		<- rp[, list(N=length(FEMALE_SANGER_ID)), by=c('FEMALE_RID','FEMALE_PTY_RUN')]
	tmp		<- subset(tmp, N>1)
	tmp		<- merge(tmp, rp, by=c('FEMALE_RID','FEMALE_PTY_RUN'))
	tmp		<- tmp[, {
				#FEMALE_TAXA<- 'PG14-UG501488-S03492'; FEMALE_PTY_RUN<- c(21)
				cat('\n',FEMALE_RID)
				z	<- subset(tmp2, TAXA1%in%FEMALE_TAXA & OVERLAP>400 & !is.na(PTY_RUN) & !PTY_RUN%in%unique(FEMALE_PTY_RUN))
				if(nrow(z)==0)
					z	<- subset(tmp2, TAXA1%in%FEMALE_TAXA & OVERLAP>200 & !is.na(PTY_RUN) & !PTY_RUN%in%unique(FEMALE_PTY_RUN))
				if(nrow(z)==0)
					z	<- subset(tmp2, TAXA1%in%FEMALE_TAXA & OVERLAP>149 & !is.na(PTY_RUN) & !PTY_RUN%in%unique(FEMALE_PTY_RUN))				
				if(nrow(z)==0)
					z	<- subset(tmp2, TAXA1%in%FEMALE_TAXA & !is.na(PTY_RUN) & !PTY_RUN%in%unique(FEMALE_PTY_RUN))
				#	exclude also run ids of the male partner
				ans		<- unique(rp[['MALE_PTY_RUN']][ which(rp[['MALE_RID']]%in%MALE_RID) ])
				z		<- subset(z, !PTY_RUN%in%ans)
				#	find min distance for each run
				z	<- z[, list(PD=min(PD)), by='PTY_RUN']
				setkey(z, PD)
				ans						<- FEMALE_PTY_RUN
				ans[duplicated(ans)]	<- z[seq_len( length(FEMALE_PTY_RUN)-length(unique(ans)) ), PTY_RUN]
				list(MALE_TAXA=MALE_TAXA, FEMALE_TAXA=FEMALE_TAXA, FEMALE_PTY_RUN=FEMALE_PTY_RUN, PTY_RUN= ans)
			}, by=c('FEMALE_RID')]
	rp		<- merge(rp, tmp, by=c('FEMALE_RID','MALE_TAXA','FEMALE_TAXA','FEMALE_PTY_RUN'), all.x=1)
	tmp		<- rp[, which(!is.na(PTY_RUN))]
	set(rp, tmp, c('FEMALE_PTY_RUN','MALE_PTY_RUN'), rp[tmp,PTY_RUN])
	set(rp, NULL, 'PTY_RUN', NULL)
	#	check
	stopifnot( nrow(subset(rp[, list(N=length(FEMALE_SANGER_ID)), by=c('FEMALE_RID','FEMALE_PTY_RUN')], N>1))==0 )
	stopifnot( nrow(subset(rp[, list(N=length(MALE_SANGER_ID)), by=c('MALE_RID','MALE_PTY_RUN')], N>1 ))==0 )
	#
	#	done with pair assignments
	#	
	#	condense pairs to list
	tmp		<- subset(rp, select=c(COUPID, MALE_TAXA, MALE_SANGER_ID, MALE_PTY_RUN))
	tmp[, PARTNER:='Male']
	setnames(tmp, c('MALE_TAXA','MALE_SANGER_ID','MALE_PTY_RUN'),c('TAXA','FILE_ID','PTY_RUN'))
	rp		<- subset(rp, select=c(COUPID, FEMALE_TAXA, FEMALE_SANGER_ID, FEMALE_PTY_RUN))
	rp[, PARTNER:='Female']
	setnames(rp, c('FEMALE_TAXA','FEMALE_SANGER_ID','FEMALE_PTY_RUN'),c('TAXA','FILE_ID','PTY_RUN'))
	rp		<- rbind(rp, tmp)
	setkey(rp, COUPID, FILE_ID, PTY_RUN)	
	# 	complete each run with individuals from outside couple so we have 22 individuals in each run
	tmp2	<- subset(sq.gd, is.finite(PD))
	ans		<- rp[, {
				cat('\n',PTY_RUN)
				#	select study participants not in couple 
				tmp	<- unlist(strsplit(COUPID,':',fixed=1))
				tmp	<- subset(dm, COHORT=='RCCS' & !STUDY_ID%in%tmp, TAXA)
				setnames(tmp, 'TAXA', 'TAXA2')
				tmp	<- merge(tmp2, tmp, by='TAXA2')				
				#	get those that are close to individuals in the current run
				tmp	<- subset(tmp, TAXA1%in%TAXA)
				setkey(tmp, PD)
				z	<- subset(tmp, OVERLAP>400)
				if(nrow(z)==0)
					z	<- subset(tmp, OVERLAP>200)
				if(nrow(z)==0)
					z	<- subset(tmp, OVERLAP>149)
				stopifnot(nrow(z)>0)
				z	<- z[, list(PD=min(PD)), by='TAXA2']
				setkey(z, PD)				
				tmp	<- integer(0)
				if(length(TAXA)<22L)
					tmp	<- seq_len(22L-length(TAXA))
				list(	TAXA= c(TAXA, z[tmp, TAXA2]), 
						PARTNER= c(PARTNER, rep('Other',length(tmp))), 
						COUPID=c(COUPID,rep('Other',length(tmp))) )				
			}, by='PTY_RUN']
	ans			<- merge(ans, subset(dm, select=c(TAXA, SANGER_ID)), by='TAXA')
	setnames(ans, 'SANGER_ID', 'FILE_ID')	
	#	to each run, add one individual with super coverage: 15172_1_43, PG14-UG500280-S02284
	ans			<- rbind(ans, as.data.table(expand.grid(TAXA='PG14-UG500280-S02284', FILE_ID='15172_1_43', PARTNER='Other', COUPID='Other', PTY_RUN=ans[, sort(unique(PTY_RUN))])))
	setkey(ans, PTY_RUN, TAXA)	
	pty.runs	<- copy(ans)
	save(pty.runs, file=file.path('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples','Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda'))
	#	add second dummy individual
	ans			<- rbind(ans, as.data.table(expand.grid(TAXA='PG14-UG501310-S03314', FILE_ID='15777_1_5', PARTNER='Other', COUPID='Other', PTY_RUN=ans[, sort(unique(PTY_RUN))])))
	setkey(ans, PTY_RUN, TAXA)
	pty.runs	<- copy(ans)
	save(pty.runs, file=file.path('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples','Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns_2dummy.rda'))		
}
######################################################################################
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

RakaiCirc.epi.get.info<- function()
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
	#	TODO Warning: found 1 female circumcised A006734
	tmp		<- rh[, which(CIRC=='Y' & SEX=='F')]
	cat('\nWarning: found female circumcised',rh[tmp, paste(RID, collapse=' ')])
	set(rh, tmp, 'CIRC', NA_integer_)
	list(rd=rd, rh=rh)
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

RakaiCouples.process.couples.161007<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	setkey(rp, COUPID)
	unique(rp)[, table(COUP_SC)]
	#	total of couples with assigned SANGER_ID in RCCS
	#   F->M    M->F seroinc seropos 
	# 	10      19      52     266 	
	merge(unique(rp), unique(subset(pty.runs, COUPID!='Other', COUPID)), by='COUPID')[, table(COUP_SC)]
	#	total of couples with SANGER_ID for which I have data
	#	F->M    M->F seroinc seropos 
	#	7       16      45     235 
	#
	#	collect runs
	#
	infiles	<- data.table(	FILE= c(	'~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161007/RCCS_161007_w270_phscout.rda' ))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_phscout.rda','',basename(FILE))]				
	#
	#	for each run: get list of pairs
	#	
	rps		<- infiles[, {
				#FILE	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161007/RCCS_161007_w270_phscout.rda'
				load(FILE)	#loads phs dtrms dtrees
				#	select likely pairs -- these are ordered
				dtrms[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]				
				tmp		<- dcast.data.table(dtrms, PAIR_ID~TYPE, value.var='WIN_OF_TYPE_P')
				set(tmp, tmp[, which(is.na(disconnected))], 'disconnected', 0)
				set(tmp, tmp[, which(is.na(int))], 'int', 0)
				set(tmp, tmp[, which(is.na(unint))], 'unint', 0)
				set(tmp, tmp[, which(is.na(cher))], 'cher', 0)
				set(tmp, tmp[, which(is.na(trans_12))], 'trans_12', 0)
				set(tmp, tmp[, which(is.na(trans_21))], 'trans_21', 0)
				tmp		<- melt.data.table(tmp, id.vars='PAIR_ID', value.name='WIN_OF_TYPE_P', variable.name='TYPE')				
				stopifnot( tmp[, list(CHECK=sum(WIN_OF_TYPE_P)),by='PAIR_ID'][, all(abs(CHECK-1)<=2*.Machine$double.eps)] )				
				dtrms	<- merge(unique(subset(dtrms, select=c(PAIR_ID, ID1, ID2, PTY_RUN, WIN_TOTAL, SCORE))), tmp, by='PAIR_ID')
				#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
				tmp		<- copy(dtrms)
				setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
				set(tmp, tmp[, which(TYPE=='trans_12')], 'TYPE', 'trans')
				set(tmp, tmp[, which(TYPE=='trans_21')], 'TYPE', 'trans_12')
				set(tmp, tmp[, which(TYPE=='trans')], 'TYPE', 'trans_21')
				dtrms	<- rbind(tmp, dtrms)
				#	first individual is always male	
				setnames(dtrms, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
				set(dtrms, dtrms[, which(TYPE=='trans_12')], 'TYPE', 'trans_mf')
				set(dtrms, dtrms[, which(TYPE=='trans_21')], 'TYPE', 'trans_fm')
				set(dtrms, NULL, 'TYPE', dtrms[, as.character(TYPE)])
				dtrms	<- merge(dtrms, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))				
				dtrms
			}, by=c('RUN','DIR','FILE')]
	save(rps, file= '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161007/RCCS_161007_w270_assignments.rda')
	#
	cat('\nNumber of couples',rps[, length(unique(COUPID))])
	setkey(rps, RUN, COUPID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	cat('\nNumber of sequence pairs from couples',nrow(unique(rps)))
	setkey(rps, RUN, COUPID, PAIR_ID)
	cat('\nNumber of pairings (including repeated sequence pairs)',nrow(unique(rps)))	
	#
	#	for each run: plot pairs	
	#	
	#for( run in rps[, unique(RUN)] )
	#{			
	run		<- 'RCCS_161007_w270'
	dir		<- subset(rps, RUN==run)[1,DIR]
	cat('\ndir is',dir,'\trun is',run)
	df		<- subset(rps, RUN==run)
	setkey(df, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)					
	#
	#	plot evidence
	#		
	tmp		<- unique(df)
	tmp[, PLOT_ID:= as.numeric(gsub('-','\\.',PAIR_ID))]
	tmp		<- tmp[order(-PLOT_ID),]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair', PAIR_ID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('trans_mf','trans_fm','unint','int','cher','disconnected'), labels=c('M transmit to F','F transmit to M','M, F are unint','M, F are intermingled','M, F are a cherry','M, F are disconnected'))])
	tmp		<- subset(tmp, select=c(PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, WIN_OF_TYPE_P, WIN_TOTAL, TYPE, LABEL))
	tmp[, WIN_OF_TYPE_N:=WIN_OF_TYPE_P*WIN_TOTAL]
	#tmp		<- melt(tmp, measure.vars=c('WIN_OF_TYPE_N', 'WIN_OF_TYPE_P'))
	#set(tmp, tmp[, which(variable=='WIN_OF_TYPE_N')],'variable','number of read windows')
	#set(tmp, tmp[, which(variable=='WIN_OF_TYPE_P')],'variable','proportion of read windows')	
	tmp2	<- subset(tmp, COUP_SC=='F->M')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC, ncol=2) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_F2M.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	tmp2	<- subset(tmp, COUP_SC=='M->F')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_grid(~COUP_SC) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_M2F.pdf',sep='')), w=25, h=max(3,0.23*nrow(tmp2)), limitsize = FALSE)
	tmp2	<- subset(tmp, COUP_SC=='seroinc')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_grid(~COUP_SC) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroInc.pdf',sep='')), w=25, h=max(3,0.23*nrow(tmp2)), limitsize = FALSE)
	tmp2	<- subset(tmp, COUP_SC=='seropos')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_grid(~COUP_SC) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroPos.pdf',sep='')), w=25, h=max(3,0.23*nrow(tmp2)), limitsize = FALSE)
	#
	#	tabulate correct classifications with phyloscanner for serodiscordant couples
	#
	#	correct: trm between M and F.
	rpa		<- subset(df, COUP_SC=='F->M' | COUP_SC=='M->F')[, list(CLASS='ancestral in either direction\nor intermingled', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='trans_mf'|TYPE=='trans_fm'|TYPE=='int'|TYPE=='cher'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(df, COUP_SC=='F->M')[, list(CLASS='ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='trans_fm'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp)
	tmp		<- subset(df, COUP_SC=='M->F')[, list(CLASS='ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='trans_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp)
	#	rank each couple
	tmp		<- rpa[order(CLASS, -CLASS_PROP),][, list(PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)), by='CLASS']
	rpa		<- merge(rpa, tmp, by=c('CLASS','PAIR_ID'))	
	setkey(rpa, CLASS, CLASS_RANK)
	#	plot by rank
	ggplot(rpa, aes(x=CLASS_PROP, y=CLASS_RANK, colour=CLASS)) + 
			geom_point() + geom_step() +
			geom_text(aes(label=PAIR_ID), size=2, nudge_x=-.05, nudge_y=.5) +
			scale_y_continuous(breaks=seq(0,50,5)) +
			scale_x_reverse(labels = scales::percent, breaks=seq(0,1,0.1)) +
			scale_colour_brewer(palette='Set1') +
			facet_grid(~CLASS) +
			labs(	x= '\nminimum proportion\n(proportion of ancestral windows out of all windows\nthat have reads from both individuals is at least x%)', 
					y='sequence pairs\n(#)\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-correctancestral.pdf',sep='')), w=12, h=7)
	#	write to file
	tmp			<- subset(rpa, CLASS=='ancestral in either direction\nor intermingled', select=c(COUPID, PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, CLASS, CLASS_PROP, CLASS_RANK))
	write.csv(tmp, row.names=FALSE, file=file.path(dir, paste(run,'-phsc-serodiscpairs-assignments.csv',sep='')) )
	#	numbers
	tmp			<- subset(rpa, CLASS=='ancestral in either direction\nor intermingled')
	cat('\nNumber of couples',tmp[, length(unique(COUPID))])
	setkey(tmp, COUPID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	cat('\nNumber of sequence pairs from couples',nrow(unique(tmp)))
	setkey(tmp, RUN, COUPID, PAIR_ID)
	cat('\nNumber of pairings (including repeated sequence pairs)',nrow(unique(tmp)))	
	#
	#	plot on proportion of assignments in epidemiologically possible direction
	#
	rpb		<- subset(df, COUP_SC=='F->M')[, list(CLASS='prop ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='trans_fm'])/sum(WIN_OF_TYPE_P[TYPE=='trans_fm'|TYPE=='trans_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]	
	tmp		<- subset(df, COUP_SC=='M->F')[, list(CLASS='prop ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='trans_mf'])/sum(WIN_OF_TYPE_P[TYPE=='trans_fm'|TYPE=='trans_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpb		<- rbind(rpb, tmp)
	tmp		<- rpb[order(CLASS, -CLASS_PROP),][, list(PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)), by='CLASS']
	rpb		<- merge(rpb, tmp, by=c('CLASS','PAIR_ID'))	
	setkey(rpb, CLASS, CLASS_RANK)
	ggplot(rpb, aes(x=CLASS_RANK, y=cumsum(CLASS_PROP))) + geom_line() + geom_point() +
			coord_cartesian(xlim=c(0, max(rpb[,CLASS_RANK])), ylim=c(0,max(rpb[,CLASS_RANK]))) +
			geom_abline(intercept=0, slope=0.5, colour='blue') +
			geom_abline(intercept=0, slope=1, colour='blue') +
			geom_text(aes(label=PAIR_ID), size=2, nudge_x=-.2, nudge_y=.8) +
			scale_y_continuous(expand=c(0,0)) +
			scale_x_continuous(expand=c(0,0)) +
			theme_bw() +
			labs(	x='\nsequence pairs\nof serodiscordant couples whose uninfected partners turns positive\n(cumulated)',
					y='# ancestral assignments in direction that is epidemiologically possible\nout of all ancestral assignments\n(cumulated)')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-only_possible_direction_assigned.pdf',sep='')), w=7, h=7)
	#
	#	plot on number of ancestral windows 
	#
	df[, WIN_OF_TYPE_N:= WIN_OF_TYPE_P*WIN_TOTAL]
	rpa		<- subset(df, COUP_SC=='F->M')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='trans_fm']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='trans_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(df, COUP_SC=='M->F')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='trans_mf']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='trans_fm'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= IN_DIR+AGAINST_DIR]
	rpa		<- melt(rpa, measure.vars=c('IN_DIR', 'AGAINST_DIR'))
	set(rpa, rpa[, which(variable=='IN_DIR')], 'variable', 'trm assignment in the only epidemiologically possible direction')
	set(rpa, rpa[, which(variable=='AGAINST_DIR')], 'variable', 'trm assignment against the only epidemiologically possible direction')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(rpa)[order(-WIN_TRM),][, list(COUPID=COUPID, PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(PAIR_ID, ' (', COUPID, ')', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PAIR_ID','COUPID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			#scale_x_discrete(labels=rpa[,PAIR_ID]) +
			labs(	x= '\nsequence pairs of couples', 
					y='number of transmission windows with direction\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trm_windows.pdf',sep='')), w=10, h=7)
	#
	#	plot on all windows 
	#	
	rpa		<- subset(df, COUP_SC=='F->M')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int'|TYPE=='cher']), NO_ANC=sum(WIN_OF_TYPE_N[!(TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int'|TYPE=='cher')])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(df, COUP_SC=='M->F')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_mf'|TYPE=='trans_mf'|TYPE=='int'|TYPE=='cher']), NO_ANC=sum(WIN_OF_TYPE_N[!(TYPE=='trans_mf'|TYPE=='trans_mf'|TYPE=='int'|TYPE=='cher')])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= ANY_ANC]
	rpa		<- melt(rpa, measure.vars=c('ANY_ANC', 'NO_ANC'))
	set(rpa, rpa[, which(variable=='ANY_ANC')], 'variable', 'transmission assignments or intermingled or cherry')
	set(rpa, rpa[, which(variable=='NO_ANC')], 'variable', 'unint or disconnected')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(rpa)[order(-WIN_TRM),][, list(COUPID=COUPID, PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(PAIR_ID, ' (', COUPID, ')', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PAIR_ID','COUPID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set2') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			#scale_x_discrete(labels=rpa[,PAIR_ID]) +
			labs(	x= '\nsequence pairs of couples', 
					y='number of windows\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_other_windows.pdf',sep='')), w=10, h=7)
	
	
	
	
	#
	#	inconsistent pairings of same sequences
	#	inconsistent pairings of same couples
	
	#}
	#
	#	re-examine phylogenies for all sero-discordant couples
	#	
	tmp		<- subset(df, COUP_SC=='F->M' | COUP_SC=='M->F')
	setkey(tmp, RUN, PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	dir		<- tmp[1,DIR]
	cat('\ndir is',dir,'\trun is',run)			
	setkey(tmp, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)				
	rpoku	<- unique(tmp)
	load( tmp[1, FILE] )	#loads phs dtrms dtrees
	invisible(sapply(seq_len(nrow(rpoku)), function(ii)
					{								
						pair.id		<- rpoku[ii, PAIR_ID]
						pty.run		<- rpoku[ii, PTY_RUN]
						dfs			<- subset(dtrees, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, TITLE:= dfs[, paste('pair', pair.id, '\n', rpoku[ii, COUP_SC], '\nid M: ', rpoku[ii, MALE_RID], ' (', rpoku[ii, MALE_SANGER_ID], ')\nid F: ', rpoku[ii, FEMALE_RID], ' (', rpoku[ii, FEMALE_SANGER_ID], ')\nrun ', pty.run, '\nwindow ', W_FROM,'-', W_TO,sep='')]]			
						plot.file	<- file.path(dir, paste(run,'-phsc-serodiscpairs-',rpoku[ii, COUP_SC],'-M-', rpoku[ii, MALE_RID],'-F-',rpoku[ii, FEMALE_RID],'-', pair.id,'.pdf',sep=''))			
						invisible(phsc.plot.selected.pairs(phs, dfs, rpoku[ii, MALE_SANGER_ID], rpoku[ii, FEMALE_SANGER_ID], plot.file=plot.file, pdf.h=150, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40))
					}))				
	
}

RakaiCouples.process.couples.161027<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	setkey(rp, COUPID)
	unique(rp)[, table(COUP_SC)]
	#	total of couples with assigned SANGER_ID in RCCS
	#   F->M    M->F seroinc seropos 
    # 	10      19      52     266 	
	tmp		<- merge(unique(rp), unique(subset(pty.runs, COUPID!='Other', COUPID)), by='COUPID')	
	tmp[, table(COUP_SC)]
	#	total of couples with SANGER_ID for which I have data
	#	F->M    M->F seroinc seropos 
	#	7       19   50      251 			total: 327
	#
	#	collect runs
	#
	infiles	<- data.table(	FILE= c(	'~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/RCCS_161027_w270_d20_phscout.rda',
										'~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/RCCS_161027_w270_d50_phscout.rda',
										'~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/RCCS_161027_w270_d200_phscout.rda'))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_phscout.rda','',basename(FILE))]				
	#
	#	for each run: get list of pairs
	#	
	rps		<- infiles[, {
				#FILE	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/RCCS_161027_w270_d20_phscout.rda'
				cat('\n',FILE)
				load(FILE)	#loads phs dtrms dtrees
				#	select likely pairs -- these are ordered
				dtrms[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]			
				tmp		<- dtrms[, list(CH=WIN_TOTAL-sum(WIN_OF_TYPE)), by=c('PTY_RUN','ID1','ID2')]
				stopifnot( tmp[, all(CH==0)] )
				tmp		<- dcast.data.table(dtrms, PAIR_ID~TYPE, value.var='WIN_OF_TYPE_P')
				set(tmp, tmp[, which(is.na(disconnected))], 'disconnected', 0)
				set(tmp, tmp[, which(is.na(int))], 'int', 0)
				set(tmp, tmp[, which(is.na(unint))], 'unint', 0)
				set(tmp, tmp[, which(is.na(cher))], 'cher', 0)
				set(tmp, tmp[, which(is.na(trans_12))], 'trans_12', 0)
				set(tmp, tmp[, which(is.na(trans_21))], 'trans_21', 0)
				tmp		<- melt.data.table(tmp, id.vars='PAIR_ID', value.name='WIN_OF_TYPE_P', variable.name='TYPE')				
				#
				dtrms	<- merge(unique(subset(dtrms, select=c(PAIR_ID, ID1, ID2, PTY_RUN, WIN_TOTAL, SCORE))), tmp, by='PAIR_ID')
				#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
				tmp		<- copy(dtrms)
				setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
				set(tmp, tmp[, which(TYPE=='trans_12')], 'TYPE', 'trans')
				set(tmp, tmp[, which(TYPE=='trans_21')], 'TYPE', 'trans_12')
				set(tmp, tmp[, which(TYPE=='trans')], 'TYPE', 'trans_21')
				dtrms	<- rbind(tmp, dtrms)
				#	first individual is always male	
				setnames(dtrms, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
				set(dtrms, dtrms[, which(TYPE=='trans_12')], 'TYPE', 'trans_mf')
				set(dtrms, dtrms[, which(TYPE=='trans_21')], 'TYPE', 'trans_fm')
				set(dtrms, NULL, 'TYPE', dtrms[, as.character(TYPE)])
				dtrms	<- merge(dtrms, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))				
				dtrms
			}, by=c('RUN','DIR','FILE')]
	save(rps, file= '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/RCCS_161027_w270_assignments.rda')
	rpsn	<- copy(rps)
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161007/RCCS_161007_w270_assignments.rda')
	rps		<- rbind(rpsn, rps)
	#
	cat('\nNumber of couples',paste(rps[, length(unique(COUPID)), by='RUN'], collapse=''))
	setkey(rps, RUN, COUPID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	cat('\nNumber of sequence pairs from couples',nrow(unique(rps)))
	setkey(rps, RUN, COUPID, PAIR_ID)
	cat('\nNumber of pairings (including repeated sequence pairs)',nrow(unique(rps)))
	#
	#	check inconsistent
	#
	if(0)
	{
		rpsn	<- copy(rps)
		rps.new<- unique(subset(rpsn, RUN=='RCCS_161027_w270_d200', select=c(COUPID,MALE_TAXA,FEMALE_TAXA)))
		rps.new[, TYPE:='NEW']	
		rps.old<- unique(subset(rps, RUN=='RCCS_161007_w270', select=c(COUPID,MALE_TAXA,FEMALE_TAXA)))
		rps.old[, TYPE:='OLD']
		tmp		<- merge(rps.new, rps.old, by=c('COUPID','MALE_TAXA','FEMALE_TAXA'),all=1)
		subset(tmp, is.na(TYPE.x))
		subset(pty.runs, TAXA=='PG14-UG503118-S05122' | TAXA=='PG14-UG503081-S05085')
		subset(pty.runs, TAXA=='PG14-UG500525-S02529' | TAXA=='PG14-UG500526-S02530')
		
		subset(rps, MALE_TAXA=='PG14-UG503118-S05122' | FEMALE_TAXA=='PG14-UG503081-S05085')
	}
	set(rps, NULL, 'RUN', rps[, factor(RUN, levels=c("RCCS_161007_w270","RCCS_161027_w270_d20","RCCS_161027_w270_d50","RCCS_161027_w270_d200"))])
	#
	#	for each run: plot pairs	
	#	
	run		<- 'RCCS_161027_w270_dxxx'
	dir		<- rps$DIR[1]
	rpp		<- subset(rps, RUN==rps$RUN[1])
	setkey(rpp, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)
	tmp		<- unique(rpp)	
	tmp[, PLOT_ID:= as.numeric(gsub('-','\\.',PAIR_ID))]
	tmp		<- tmp[order(-PLOT_ID),]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair', PAIR_ID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL)), rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('trans_mf','trans_fm','unint','int','cher','disconnected'), labels=c('M transmit to F','F transmit to M','M, F are unint','M, F are intermingled','M, F are a cherry','M, F are disconnected'))])
	tmp		<- subset(tmp, select=c(RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, WIN_OF_TYPE_P, WIN_TOTAL, TYPE, LABEL))
	tmp[, WIN_OF_TYPE_N:=WIN_OF_TYPE_P*WIN_TOTAL]
	setkey(tmp, COUP_SC, LABEL, RUN, TYPE)
	#	F->M
	tmp2	<- subset(tmp, COUP_SC=='F->M')
	ggplot(tmp2, aes(x=RUN, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC+LABEL, ncol=1) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_F2M.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	#	M->F
	tmp2	<- subset(tmp, COUP_SC=='M->F')
	ggplot(tmp2, aes(x=RUN, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC+LABEL, ncol=1) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_M2F.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	#	seroinc
	tmp2	<- subset(tmp, COUP_SC=='seroinc')
	ggplot(tmp2, aes(x=RUN, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC+LABEL, ncol=1) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroInc.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	#	seropos
	tmp2	<- subset(tmp, COUP_SC=='seropos')
	ggplot(tmp2, aes(x=RUN, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC+LABEL, ncol=1) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroPos.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	#
	#	plot number of directional trm assignments in only possible direction 
	#
	rps[, WIN_OF_TYPE_N:= WIN_OF_TYPE_P*WIN_TOTAL]
	rpa		<- subset(rps, COUP_SC=='F->M')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='trans_fm']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='trans_mf'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(rps, COUP_SC=='M->F')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='trans_mf']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='trans_fm'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= IN_DIR+AGAINST_DIR]
	rpa		<- melt(rpa, measure.vars=c('IN_DIR', 'AGAINST_DIR'))
	set(rpa, rpa[, which(variable=='IN_DIR')], 'variable', 'trm assignment in the only epidemiologically possible direction')
	set(rpa, rpa[, which(variable=='AGAINST_DIR')], 'variable', 'trm assignment against the only epidemiologically possible direction')	
	setkey(rpa, PAIR_ID)		
	tmp		<- unique(subset(rpa, RUN==rpa$RUN[1]))[order(-WIN_TRM),][, list(COUPID=COUPID, MALE_SANGER_ID=MALE_SANGER_ID, FEMALE_SANGER_ID=FEMALE_SANGER_ID, PTY_RUN=PTY_RUN, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PTY_RUN','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~RUN) +
			coord_flip() +
			labs(	x= '\nsequence pairs of couples', 
					y='number of transmission windows with direction\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trm_windows.pdf',sep='')), w=20, h=10)
	#
	#
	#
	rpa		<- subset(rps, COUP_SC=='F->M')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), CHER=WIN_OF_TYPE_N[TYPE=='cher'], NO_ANC=sum(WIN_OF_TYPE_N[!(TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int'|TYPE=='cher')])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(rps, COUP_SC=='M->F')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_mf'|TYPE=='trans_mf'|TYPE=='int']), CHER=WIN_OF_TYPE_N[TYPE=='cher'], NO_ANC=sum(WIN_OF_TYPE_N[!(TYPE=='trans_mf'|TYPE=='trans_mf'|TYPE=='int'|TYPE=='cher')])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= ANY_ANC]
	rpa		<- melt(rpa, measure.vars=c('ANY_ANC', 'NO_ANC','CHER'))
	set(rpa, rpa[, which(variable=='ANY_ANC')], 'variable', 'transmission assignments or intermingled')
	set(rpa, rpa[, which(variable=='NO_ANC')], 'variable', 'unint or disconnected')
	set(rpa, rpa[, which(variable=='CHER')], 'variable', 'cherry')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(subset(rpa, RUN==rpa$RUN[1]))[order(-WIN_TRM),][, list(COUPID=COUPID, MALE_SANGER_ID=MALE_SANGER_ID, FEMALE_SANGER_ID=FEMALE_SANGER_ID, PTY_RUN=PTY_RUN, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PTY_RUN','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set2') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~RUN) +
			coord_flip() +
			labs(	x= '\nsequence pairs of couples', 
					y='number of windows\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_other_windows.pdf',sep='')), w=20, h=10)
	
	#
	#	inspect rogue case by window
	#
	dr		<- as.data.table(read.csv('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/couple_374_rogues_OR.csv'))
	setnames(dr, c('window.start','J110207','G105508'), c('W_FROM','15778_1_82','15758_1_76'))
	tmp		<- melt(subset(dr, select=c('W_FROM','15778_1_82','15758_1_76')), id.vars='W_FROM',value.name='BRL',variable.name='ID')  
	dr		<- melt(subset(dr, select=c('W_FROM','rogue_15778_1_82','rogue_15758_1_76')), id.vars='W_FROM',value.name='ROGUE',variable.name='ID')
	set(dr, NULL, 'ID', dr[, gsub('rogue_','',ID)])
	dr		<- merge(tmp, dr, by=c('W_FROM','ID'))
	
	infiles<- c(	'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d20_rerun/ptyr66_trmStatsPerWindow.rda',
					'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d50_rerun/ptyr66_trmStatsPerWindow.rda',
					'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d200_rerun/ptyr66_trmStatsPerWindow.rda'	)
	id.m	<- '15778_1_82'
	id.f	<- '15758_1_76'	
	df		<- phsc.get.assignments.by.window.for.couple(id1=paste(id.m,'.bam',sep=''), id2=paste(id.f,'.bam',sep=''), infiles)
	set(df, NULL, 'ID1', df[, gsub('\\.bam','',ID1)])
	set(df, NULL, 'ID2', df[, gsub('\\.bam','',ID2)])
	setnames(df, colnames(df), gsub('_rerun','',gsub('couples_','',gsub('Rakai_ptoutput_','',colnames(df)))))
	
	setnames(dr, c('ID','BRL','ROGUE'), c('ID1','BRL1','ROGUE1'))
	df		<- merge(df,dr,by=c('ID1','W_FROM'))
	setnames(dr, c('ID1','BRL1','ROGUE1'), c('ID2','BRL2','ROGUE2'))
	df		<- merge(df,dr,by=c('ID2','W_FROM'))
	df		<- subset(df, !is.na(BRL1) & !is.na(BRL2) & !is.na(ROGUE1) & !is.na(ROGUE2))
	#
	#	branches among rogues and non rogues
	#
	tmp		<- subset(df, select=c(BRL1,ROGUE1))
	setnames(tmp, c('BRL1','ROGUE1'), c('BRL','ROGUE'))
	tmp2	<- subset(df, select=c(BRL2,ROGUE2))
	setnames(tmp2, c('BRL2','ROGUE2'), c('BRL','ROGUE'))	
	tmp		<- rbind(tmp, tmp2)
	ggplot(tmp, aes(x=BRL, fill=factor(ROGUE))) + geom_histogram() + facet_grid(~ROGUE)
	
	subset(df, (BRL1>0.05 & !ROGUE1) | (BRL2>0.05 & !ROGUE2))
	subset(df, (BRL1<0.07 & ROGUE1) | (BRL2<0.07 & ROGUE2))
	#	calculate patristic distance..
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/RCCS_161027_w270_d200_phscout.rda')
	ph		<- phs[[ subset(dtrees, PTY_RUN==66 & W_FROM==1225)[, IDX] ]]
	tmp		<- as.matrix(cophenetic.phylo(ph))
	max( tmp['15758_1_76.bam_read_7_count_1', grepl('15758_1_76', colnames(tmp))] )
	
	#	collect branch lengths from individuals that I classified manually
	tmp		<- merge(subset(dtrees, PTY_RUN==66, c(W_FROM, IDX)), subset(df, select=c(ID1, W_FROM, ROGUE1)),by='W_FROM')
	setnames(tmp, c('ID1','ROGUE1'), c('ID','ROGUE'))
	dr		<- tmp[, {
				#	IDX<- 26801; ID1<- '15778_1_82'
				ph	<- phs[[ IDX ]]
				z	<- ph$tip.label[!grepl(ID,ph$tip.label)]
				ph	<- drop.tip(ph,z)
				list(BRL=ph$edge.length) 
			}, by=c('W_FROM','ID','ROGUE')]	
	tmp		<- merge(subset(dtrees, PTY_RUN==66, c(W_FROM, IDX)), subset(df, select=c(ID2, W_FROM, ROGUE2)),by='W_FROM')
	setnames(tmp, c('ID2','ROGUE2'), c('ID','ROGUE'))
	tmp		<- tmp[, {
				#	IDX<- 26801; ID1<- '15778_1_82'
				ph	<- phs[[ IDX ]]
				z	<- ph$tip.label[!grepl(ID,ph$tip.label)]
				ph	<- drop.tip(ph,z)
				list(BRL=ph$edge.length) 
			}, by=c('W_FROM','ID','ROGUE')]
	dr		<- rbind(dr, tmp)	
	#	calculate prob that the max BRL is > x under the "null" that all distances are from the same distribution
	#	using Weibull as "null" model because (1) for x positive and (2) it is easy to calculate 1-CDF(max X_i))
	require(fitdistrplus)
	dp		<- dr[, {
				x			<- max(BRL)
				p			<- NA_real_
				z			<- BRL[ BRL>1e-5 ]
				if(length(z)<=1)
					cat('\n', W_FROM, ID)
				if(length(z)>5)
				{	
					cat('\n', W_FROM, ID)
					w.mle		<- fitdist(z, 'weibull')	
					#	denscomp(w.mle, addlegend=FALSE, xlab='weibull') 	#use this to compare fit to data
					w.k			<- w.mle$estimate['shape']
					w.l			<- w.mle$estimate['scale']				
					p			<- 1 - ( 1 - exp( -( max(z)/w.l )^w.k ) )^length(z)					
				}
				list(P=p, BRL.mx=max(BRL))
			}, by=c('W_FROM','ID','ROGUE')]
	
	require(gamlss)
	dp		<- dr[, {
				x			<- max(BRL)
				p			<- NA_real_
				z			<- BRL[ BRL>1e-5 ]
				if(length(z)>3)	# don't think makes sense to fit Weibull to 3 data points
				{	
					# cat('\n', W_FROM, ID) # for debugging
					w		<- gamlss(data=data.table(BRL=z), formula=BRL~1, family=WEI, trace=FALSE)
					w.l		<- exp(coef(w, what='mu'))
					w.k		<- exp(coef(w, what='sigma'))
					p			<- 1 - ( 1 - exp( -( max(z)/w.l )^w.k ) )^length(z)					
				}
				list(P=p, BRL.mx=max(BRL))
			}, by=c('W_FROM','ID','ROGUE')]	
	dp[, BRL.mx.c:= cut(BRL.mx, breaks=c(0,0.04,Inf), labels=c('<4%','>=4%'))]		
	subset(dp, !is.na(P))[, table(ROGUE, BRL.mx.c)]
	subset(dp, !is.na(P))[, table(ROGUE, P<.05)]
	subset(dp, !is.na(P))[, table(ROGUE, P<.0001)]	
	ggplot(subset(dp, !is.na(P)), aes(x=factor(ROGUE), y=P)) + 
			geom_boxplot() + 
			facet_grid(.~BRL.mx.c) +
			labs(x='manual classification rogues',y='prob that max BRL is > observed max\nunder Weibull model')
	
	#
	#	inspect couple K061956:J061939 ( M:15103_1_74 F:15861_1_22 run:91 )
	#
	infiles<- c(	'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d20_rerun/ptyr91_trmStatsPerWindow.rda',
					'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d50_rerun/ptyr91_trmStatsPerWindow.rda',
					'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d200_rerun/ptyr91_trmStatsPerWindow.rda'	)
	id.m	<- '15103_1_74'
	id.f	<- '15861_1_22'	
	df		<- phsc.get.assignments.by.window.for.couple(id1=paste(id.m,'.bam',sep=''), id2=paste(id.f,'.bam',sep=''), infiles)
	set(df, NULL, 'ID1', df[, gsub('\\.bam','',ID1)])
	set(df, NULL, 'ID2', df[, gsub('\\.bam','',ID2)])
	setnames(df, colnames(df), gsub('_rerun','',gsub('couples_','',gsub('Rakai_ptoutput_','',colnames(df)))))

	infiles<- c(	'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d20_rerun/ptyr62_trmStatsPerWindow.rda',
					'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d50_rerun/ptyr62_trmStatsPerWindow.rda',
					'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d200_rerun/ptyr62_trmStatsPerWindow.rda'	)
	id.m	<- '15861_1_26'
	id.f	<- '15861_1_22'	
	df		<- phsc.get.assignments.by.window.for.couple(id1=paste(id.m,'.bam',sep=''), id2=paste(id.f,'.bam',sep=''), infiles)
	set(df, NULL, 'ID1', df[, gsub('\\.bam','',ID1)])
	set(df, NULL, 'ID2', df[, gsub('\\.bam','',ID2)])
	setnames(df, colnames(df), gsub('_rerun','',gsub('couples_','',gsub('Rakai_ptoutput_','',colnames(df)))))
}


RakaiCouples.process.couples.160930<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	setkey(rp, COUPID)
	unique(rp)[, table(COUP_SC)]
	#	total of couples with assigned SANGER_ID in RCCS
	#   F->M    M->F seroinc seropos 
	# 	10      19      52     266 	
	merge(unique(rp), unique(subset(pty.runs, COUPID!='Other', COUPID)), by='COUPID')[, table(COUP_SC)]
	#	total of couples with SANGER_ID for which I have data
	#	F->M    M->F seroinc seropos 
	#	7       16      45     235 
	#
	#	collect runs
	#
	infiles	<- data.table(	FILE= c(	'~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/160930/RCCS_160930_w270_phscout.rda' ))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_phscout.rda','',basename(FILE))]				
	#
	#	for each run: get list of pairs
	#	
	rpso		<- infiles[, {
				#F_TRM	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/RCCS_160919_w270_trmStats.rda'; F_PH	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/RCCS_160919_w270_trees.rda'
				load(FILE)	#loads phs dtrms dtrees
				#	select likely pairs -- these are ordered
				dtrms[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]				
				tmp		<- dcast.data.table(dtrms, PAIR_ID~TYPE, value.var='WIN_OF_TYPE_P')
				set(tmp, tmp[, which(is.na(disconnected))], 'disconnected', 0)
				set(tmp, tmp[, which(is.na(sib))], 'sib', 0)
				set(tmp, tmp[, which(is.na(int))], 'int', 0)
				set(tmp, tmp[, which(is.na(anc_12))], 'anc_12', 0)
				set(tmp, tmp[, which(is.na(anc_21))], 'anc_21', 0)
				tmp		<- melt.data.table(tmp, id.vars='PAIR_ID', value.name='WIN_OF_TYPE_P', variable.name='TYPE')
				dtrms	<- merge(unique(subset(dtrms, select=c(PAIR_ID, ID1, ID2, PTY_RUN, WIN_TOTAL, SCORE))), tmp, by='PAIR_ID')
				#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
				tmp		<- copy(dtrms)
				setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
				set(tmp, tmp[, which(TYPE=='anc_12')], 'TYPE', 'anc')
				set(tmp, tmp[, which(TYPE=='anc_21')], 'TYPE', 'anc_12')
				set(tmp, tmp[, which(TYPE=='anc')], 'TYPE', 'anc_21')
				dtrms	<- rbind(tmp, dtrms)
				#	first individual is always male	
				setnames(dtrms, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
				set(dtrms, dtrms[, which(TYPE=='anc_12')], 'TYPE', 'anc_mf')
				set(dtrms, dtrms[, which(TYPE=='anc_21')], 'TYPE', 'anc_fm')
				set(dtrms, NULL, 'TYPE', dtrms[, as.character(TYPE)])
				dtrms	<- merge(dtrms, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))				
				dtrms
			}, by=c('RUN','DIR','FILE')]
	#
	cat('\nNumber of couples',rps[, length(unique(COUPID))])
	setkey(rps, RUN, COUPID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	cat('\nNumber of sequence pairs from couples',nrow(unique(rps)))
	setkey(rps, RUN, COUPID, PAIR_ID)
	cat('\nNumber of pairings (including repeated sequence pairs)',nrow(unique(rps)))	
	#
	#	for each run: plot pairs	
	#	
	#for( run in rps[, unique(RUN)] )
	#{		
	run		<- 'RCCS_160919_w270'
	run		<- 'RCCS_160930_w270'
	dir		<- subset(rps, RUN==run)[1,DIR]
	cat('\ndir is',dir,'\trun is',run)
	df		<- subset(rps, RUN==run)
	setkey(df, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)					
	#
	#	plot evidence
	#		
	tmp		<- unique(df)
	tmp[, PLOT_ID:= as.numeric(gsub('-','\\.',PAIR_ID))]
	tmp		<- tmp[order(-PLOT_ID),]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair', PAIR_ID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('anc_mf','anc_fm','sib','int','disconnected'), labels=c('M ancestral to F','F ancestral to M','M, F are siblings','M, F are intermingled','M, F are disconnected'))])
	tmp		<- subset(tmp, select=c(PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, WIN_OF_TYPE_P, WIN_TOTAL, TYPE, LABEL))
	tmp[, WIN_OF_TYPE_N:=WIN_OF_TYPE_P*WIN_TOTAL]
	#tmp		<- melt(tmp, measure.vars=c('WIN_OF_TYPE_N', 'WIN_OF_TYPE_P'))
	#set(tmp, tmp[, which(variable=='WIN_OF_TYPE_N')],'variable','number of read windows')
	#set(tmp, tmp[, which(variable=='WIN_OF_TYPE_P')],'variable','proportion of read windows')	
	tmp2	<- subset(tmp, COUP_SC=='F->M')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M ancestral to F'="#9E0142",'F ancestral to M'="#F46D43",'M, F are siblings'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC, ncol=2) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_F2M.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	tmp2	<- subset(tmp, COUP_SC=='M->F')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M ancestral to F'="#9E0142",'F ancestral to M'="#F46D43",'M, F are siblings'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_grid(~COUP_SC) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_M2F.pdf',sep='')), w=25, h=max(3,0.23*nrow(tmp2)), limitsize = FALSE)
	tmp2	<- subset(tmp, COUP_SC=='seroinc')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M ancestral to F'="#9E0142",'F ancestral to M'="#F46D43",'M, F are siblings'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_grid(~COUP_SC) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroInc.pdf',sep='')), w=25, h=max(3,0.23*nrow(tmp2)), limitsize = FALSE)
	tmp2	<- subset(tmp, COUP_SC=='seropos')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M ancestral to F'="#9E0142",'F ancestral to M'="#F46D43",'M, F are siblings'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_grid(~COUP_SC) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroPos.pdf',sep='')), w=25, h=max(3,0.23*nrow(tmp2)), limitsize = FALSE)
	#
	#	tabulate correct classifications with phyloscanner for serodiscordant couples
	#
	#	correct: trm between M and F.
	rpa		<- subset(df, COUP_SC=='F->M' | COUP_SC=='M->F')[, list(CLASS='ancestral in either direction\nor intermingled', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='anc_mf'|TYPE=='anc_fm'|TYPE=='int'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(df, COUP_SC=='F->M')[, list(CLASS='ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='anc_fm'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp)
	tmp		<- subset(df, COUP_SC=='M->F')[, list(CLASS='ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='anc_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp)
	#	rank each couple
	tmp		<- rpa[order(CLASS, -CLASS_PROP),][, list(PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)), by='CLASS']
	rpa		<- merge(rpa, tmp, by=c('CLASS','PAIR_ID'))	
	setkey(rpa, CLASS, CLASS_RANK)
	#	plot by rank
	ggplot(rpa, aes(x=CLASS_PROP, y=CLASS_RANK, colour=CLASS)) + 
			geom_point() + geom_step() +
			geom_text(aes(label=PAIR_ID), size=2, nudge_x=-.05, nudge_y=.5) +
			scale_y_continuous(breaks=seq(0,50,5)) +
			scale_x_reverse(labels = scales::percent, breaks=seq(0,1,0.1)) +
			scale_colour_brewer(palette='Set1') +
			facet_grid(~CLASS) +
			labs(	x= '\nminimum proportion\n(proportion of ancestral windows out of all windows\nthat have reads from both individuals is at least x%)', 
					y='sequence pairs\n(#)\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-correctancestral.pdf',sep='')), w=12, h=7)
	#	write to file
	tmp			<- subset(rpa, CLASS=='ancestral in either direction\nor intermingled', select=c(COUPID, PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, CLASS, CLASS_PROP, CLASS_RANK))
	write.csv(tmp, row.names=FALSE, file=file.path(dir, paste(run,'-phsc-serodiscpairs-assignments.csv',sep='')) )
	#	numbers
	tmp			<- subset(rpa, CLASS=='ancestral in either direction\nor intermingled')
	cat('\nNumber of couples',tmp[, length(unique(COUPID))])
	setkey(tmp, COUPID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	cat('\nNumber of sequence pairs from couples',nrow(unique(tmp)))
	setkey(tmp, RUN, COUPID, PAIR_ID)
	cat('\nNumber of pairings (including repeated sequence pairs)',nrow(unique(tmp)))	
	#
	#	plot on proportion of assignments in epidemiologically possible direction
	#
	rpb		<- subset(df, COUP_SC=='F->M')[, list(CLASS='prop ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='anc_fm'])/sum(WIN_OF_TYPE_P[TYPE=='anc_fm'|TYPE=='anc_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]	
	tmp		<- subset(df, COUP_SC=='M->F')[, list(CLASS='prop ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='anc_mf'])/sum(WIN_OF_TYPE_P[TYPE=='anc_fm'|TYPE=='anc_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpb		<- rbind(rpb, tmp)
	tmp		<- rpb[order(CLASS, -CLASS_PROP),][, list(PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)), by='CLASS']
	rpb		<- merge(rpb, tmp, by=c('CLASS','PAIR_ID'))	
	setkey(rpb, CLASS, CLASS_RANK)
	ggplot(rpb, aes(x=CLASS_RANK, y=cumsum(CLASS_PROP))) + geom_line() + geom_point() +
			coord_cartesian(xlim=c(0, max(rpb[,CLASS_RANK])), ylim=c(0,max(rpb[,CLASS_RANK]))) +
			geom_abline(intercept=0, slope=0.5, colour='blue') +
			geom_abline(intercept=0, slope=1, colour='blue') +
			geom_text(aes(label=PAIR_ID), size=2, nudge_x=-.2, nudge_y=.8) +
			scale_y_continuous(expand=c(0,0)) +
			scale_x_continuous(expand=c(0,0)) +
			theme_bw() +
			labs(	x='\nsequence pairs\nof serodiscordant couples whose uninfected partners turns positive\n(cumulated)',
					y='# ancestral assignments in direction that is epidemiologically possible\nout of all ancestral assignments\n(cumulated)')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-only_possible_direction_assigned.pdf',sep='')), w=7, h=7)
	#
	#	plot on number of ancestral windows 
	#
	df[, WIN_OF_TYPE_N:= WIN_OF_TYPE_P*WIN_TOTAL]
	rpa		<- subset(df, COUP_SC=='F->M')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='anc_fm']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='anc_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(df, COUP_SC=='M->F')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='anc_mf']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='anc_fm'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= IN_DIR+AGAINST_DIR]
	rpa		<- melt(rpa, measure.vars=c('IN_DIR', 'AGAINST_DIR'))
	set(rpa, rpa[, which(variable=='IN_DIR')], 'variable', 'trm assignment in the only epidemiologically possible direction')
	set(rpa, rpa[, which(variable=='AGAINST_DIR')], 'variable', 'trm assignment against the only epidemiologically possible direction')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(rpa)[order(-WIN_TRM),][, list(COUPID=COUPID, PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(PAIR_ID, ' (', COUPID, ')', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PAIR_ID','COUPID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			#scale_x_discrete(labels=rpa[,PAIR_ID]) +
			labs(	x= '\nsequence pairs of couples', 
					y='number of ancestral windows\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trm_windows.pdf',sep='')), w=10, h=7)
	#
	#	plot on all windows 
	#	
	rpa		<- subset(df, COUP_SC=='F->M')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='anc_fm'|TYPE=='anc_mf'|TYPE=='int']), NO_ANC=sum(WIN_OF_TYPE_N[!(TYPE=='anc_fm'|TYPE=='anc_mf'|TYPE=='int')])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(df, COUP_SC=='M->F')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='anc_mf'|TYPE=='anc_mf'|TYPE=='int']), NO_ANC=sum(WIN_OF_TYPE_N[!(TYPE=='anc_fm'|TYPE=='anc_mf'|TYPE=='int')])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= ANY_ANC]
	rpa		<- melt(rpa, measure.vars=c('ANY_ANC', 'NO_ANC'))
	set(rpa, rpa[, which(variable=='ANY_ANC')], 'variable', 'ancestral assignments or intermingled')
	set(rpa, rpa[, which(variable=='NO_ANC')], 'variable', 'sibling or disconnected')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(rpa)[order(-WIN_TRM),][, list(COUPID=COUPID, PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(PAIR_ID, ' (', COUPID, ')', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PAIR_ID','COUPID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set2') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			#scale_x_discrete(labels=rpa[,PAIR_ID]) +
			labs(	x= '\nsequence pairs of couples', 
					y='number of ancestral windows\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_other_windows.pdf',sep='')), w=10, h=7)
	
	
	
	
	#
	#	inconsistent pairings of same sequences
	#	inconsistent pairings of same couples
	
	#}
	#
	#	re-examine phylogenies for all sero-discordant couples
	#	
	tmp		<- subset(df, COUP_SC=='F->M' | COUP_SC=='M->F')
	setkey(tmp, RUN, PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	dir		<- tmp[1,DIR]
	cat('\ndir is',dir,'\trun is',run)			
	setkey(tmp, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)				
	rpoku	<- unique(tmp)
	load( tmp[1, FILE] )	#loads phs dtrms dtrees
	invisible(sapply(seq_len(nrow(rpoku)), function(ii)
					{								
						pair.id		<- rpoku[ii, PAIR_ID]
						pty.run		<- rpoku[ii, PTY_RUN]
						dfs			<- subset(dtrees, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, TITLE:= dfs[, paste('pair', pair.id, '\n', rpoku[ii, COUP_SC], '\nid M: ', rpoku[ii, MALE_RID], ' (', rpoku[ii, MALE_SANGER_ID], ')\nid F: ', rpoku[ii, FEMALE_RID], ' (', rpoku[ii, FEMALE_SANGER_ID], ')\nrun ', pty.run, '\nwindow ', W_FROM,'-', W_TO,sep='')]]			
						plot.file	<- file.path(dir, paste(run,'-phsc-serodiscpairs-',rpoku[ii, COUP_SC],'-M-', rpoku[ii, MALE_RID],'-F-',rpoku[ii, FEMALE_RID],'-', pair.id,'.pdf',sep=''))			
						invisible(phsc.plot.selected.pairs(phs, dfs, rpoku[ii, MALE_SANGER_ID], rpoku[ii, FEMALE_SANGER_ID], plot.file=plot.file, pdf.h=150, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40))
					}))				
	
}

RakaiCouples.save.couples.to.rda<- function()
{
	require(data.table)
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples"
	#	get epi info
	tmp		<- RakaiCirc.epi.get.info()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#	get sequence info and add to recipients
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_160816.rda')		
	rs		<- subset(rs, !is.na(VISIT))	
	
	#	load combination of all couples sequences
	infile	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda'
	load
	infile	<- '~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/Pangea_Couples.csv'
	rc		<- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	setnames(rc, c('female.PANGEA.ID','male.PANGEA.ID'), c('female.TAXA','male.TAXA'))
	setnames(rc, colnames(rc), gsub('\\.','_',toupper(colnames(rc))))
	setnames(rc, c('MALE_RCCS_STUDYID','FEMALE_RCCS_STUDYID'), c('MALE_RID','FEMALE_RID'))
	set(rc, NULL, 'MALE_DATE', rc[, hivc.db.Date2numeric(as.Date(MALE_DATE))])
	set(rc, NULL, 'FEMALE_DATE', rc[, hivc.db.Date2numeric(as.Date(FEMALE_DATE))])
	set(rc, NULL, 'MALE_LASTNEGDATE', rc[, hivc.db.Date2numeric(as.Date(MALE_LASTNEGDATE))])
	set(rc, NULL, 'FEMALE_LASTNEGDATE', rc[, hivc.db.Date2numeric(as.Date(FEMALE_LASTNEGDATE))])	
	set(rc, NULL, 'MALE_FIRSTPOSDATE', rc[, hivc.db.Date2numeric(as.Date(MALE_FIRSTPOSDATE))])
	set(rc, NULL, 'FEMALE_FIRSTPOSDATE', rc[, hivc.db.Date2numeric(as.Date(FEMALE_FIRSTPOSDATE))])
	#
	#	get all possible unique pairings of SANGER_IDs from couples
	#
	rp		<- subset(rc, !is.na(MALE_SANGER_ID), c(COUPID, MALE_RID, MALE_SANGER_ID, MALE_TAXA ))
	tmp		<- subset(rc, !is.na(FEMALE_SANGER_ID), c(COUPID, FEMALE_RID, FEMALE_SANGER_ID, FEMALE_TAXA ))	
	rp		<- merge(rp, tmp, by='COUPID')
	setkey(rp, FEMALE_SANGER_ID, MALE_SANGER_ID)
	rp		<- unique(rp)
	#
	#	add epi info 
	#	for males
	tmp		<- subset(rd, !is.na(PID), select=c(RID, SEX, REGION, COMM_NUM, HH_NUM, BIRTHDATE, LASTNEGVIS, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, RELIGION))
	setkey(tmp, RID)
	tmp		<- unique(tmp)
	setnames(tmp, colnames(tmp), paste('MALE_',colnames(tmp),sep=''))
	rp		<- merge(rp,tmp,by='MALE_RID')
	#	for females
	tmp		<- subset(rd, !is.na(PID), select=c(RID, SEX, REGION, COMM_NUM, HH_NUM, BIRTHDATE, LASTNEGVIS, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, RELIGION))
	setkey(tmp, RID)
	tmp		<- unique(tmp)
	setnames(tmp, colnames(tmp), paste('FEMALE_',colnames(tmp),sep=''))
	rp		<- merge(rp,tmp,by='FEMALE_RID')	
	#
	#	get info on sequences: date of sampling
	#	for males
	tmp		<- unique(subset(rs, !is.na(SID), select=c(SID, DATE)))
	setnames(tmp, c('SID','DATE'), c('MALE_SANGER_ID','MALE_SEQDATE'))
	rp		<- merge(rp, tmp, by='MALE_SANGER_ID')
	tmp		<- unique(subset(rs, !is.na(SID), select=c(SID, DATE)))
	setnames(tmp, c('SID','DATE'), c('FEMALE_SANGER_ID','FEMALE_SEQDATE'))	
	rp		<- merge(rp, tmp, by='FEMALE_SANGER_ID')				
	#
	#	define direction based on seroconversion dates
	#
	rp[, COUP_SC:='seropos']
	set(rp, rp[, which(MALE_LASTNEGDATE>=FEMALE_FIRSTPOSDATE)],'COUP_SC','F->M')
	set(rp, rp[, which(FEMALE_LASTNEGDATE>=MALE_FIRSTPOSDATE)],'COUP_SC','M->F')
	set(rp, rp[, which(FEMALE_FIRSTPOSDATE==MALE_FIRSTPOSDATE)],'COUP_SC','seroinc')
	
	save(rp, file="~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")	
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