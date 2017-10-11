######################################################################################

project.Bezemer.VLIntros<- function()
{
	vli.bez.LSD()
	#vli.bez.FastTrees()
}

vli.asr.get.tips.in.NL.subtrees<- function(ph, da, dm, id.regex, select.descendants.regex)
{
	#	determine max parsimony state transitions into NL 
	setnames(da, colnames(da), gsub('\\.','_',toupper(colnames(da))))
	da	<- subset(da, HOSTS=='NL')
	ans	<- da[, list(TAXA=ph$tip.label[ unlist(Descendants(ph, ROOT_NOS, type='tips')) ] ), by='UNIQUE_SPLITS']
	ans[, ID:=gsub(id.regex,'\\1',TAXA)]
	ans	<- subset(ans, grepl(select.descendants.regex,TAXA))
	ans	<- merge(dm, ans, by='ID')	# we will need this below too
	set(ans, NULL, 'TAXA', NULL)
	setnames(ans, 'UNIQUE_SPLITS', 'SUBTREE_ID')
	ans	
}

vli.estimate.proportion.introductions.overall <- function(ds, sampling.prob=0.43, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), seed=42, bs.len=1e3)
{
	set.seed(seed)	
	dsa		<- subset(ds, select=c(SUBTREE_ID, SUBTREE_SIZE))
	tmp		<- dsa[, list(BS=seq_len(bs.len), SUBTREE_UNOBS=rnbinom(bs.len, SUBTREE_SIZE, prob=sampling.prob)), by='SUBTREE_ID']
	dsa		<- merge(dsa, tmp, by='SUBTREE_ID')
	ans		<- dsa[, list(DUMMY=1, STAT='within_NL', N=sum(SUBTREE_SIZE)+sum(SUBTREE_UNOBS)-1), by='BS']
	tmp		<- ds[, list(DUMMY=1, N=length(SUBTREE_ID)), by='ANCESTRAL_STATE']
	tmp		<- melt(tmp, measure.vars='ANCESTRAL_STATE', value.name='STAT')
	tmp[, variable:=NULL]
	tmp		<- merge(expand.grid(BS=seq_len(bs.len), DUMMY=1),tmp,by='DUMMY')
	ans		<- rbind(ans, tmp)
	ans[, DUMMY:=NULL]
	ans		<- merge(ans, ans[, list(TOTAL= sum(N)), by='BS'], by='BS')
	ans[, PROP:= N/TOTAL]
	ans		<- ans[, list(P=ans.quantiles, Q=quantile(PROP, p=ans.quantiles)), by='STAT']
	ans	
}

vli.estimate.subtree.distribution.by.sexual.orientation <- function(dss, sampling.prob=0.43, resample.intro.n=3, resample.missing.n=3, sxo.select= c('MSM','hsx'), size.breaks=c(0,1,2,5,10,500), size.labels=c('1','2','3-5','6-10','>10'), intro.year=2000, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), seed=42)	
{
	#sampling.prob=0.43; resample.intro.n=3; resample.missing.n=3; sxo.select= c('MSM','hsx'); intro.year=2000; ans.quantiles=c(0.025,0.25,0.5,0.75,0.975); seed=42
	set.seed(seed)
		
	ans	<- subset(dss, select=c(ST, ASR_TYPE, REP, SUBTREE_ID, ANCESTRAL_STATE, MIN_SAMPLING_DATE, SUBTREE_SIZE))
	ans[, INTRO_SMPL:=1L]		
	ans	<- subset(ans, floor(MIN_SAMPLING_DATE)>=intro.year)
	tmp	<- ans[, list(MISSING_SMPL=seq_len(resample.missing.n), SUBTREE_UNOBS=rnbinom(resample.missing.n, SUBTREE_SIZE, prob=sampling.prob)), by=c('SUBTREE_ID','INTRO_SMPL')]
	ans	<- merge(ans, tmp, by=c('SUBTREE_ID','INTRO_SMPL'))
		
	tmp	<- copy(dss)
	tmp	<- melt(tmp, id.vars=c('SUBTREE_ID','SUBTREE_SIZE','REP'), measure.vars=colnames(tmp)[grepl('SUBTREE_SXO_',colnames(tmp))], variable.name='SXO')
	set(tmp, NULL, 'SXO', tmp[, gsub('SUBTREE_SXO_','',SXO)])
	tmp[, SXO_PROP:= value/SUBTREE_SIZE]
	tmp	<- subset(tmp, select=c(SUBTREE_ID,REP,SXO,SXO_PROP))
	tmp	<- subset(tmp, SXO%in%sxo.select & SXO_PROP>0)
	
	ans	<- merge(ans, tmp, by=c('SUBTREE_ID','REP'),allow.cartesian=TRUE)	
	tmp	<- ans[, list(SXO=SXO, NL_TRMS_SXO_RND=SUBTREE_SIZE*SXO_PROP+as.vector(rmultinom(1, SUBTREE_UNOBS, SXO_PROP))), by=c('SUBTREE_ID','INTRO_SMPL','MISSING_SMPL')]
	ans	<- merge(ans, tmp, by=c('SUBTREE_ID','INTRO_SMPL','MISSING_SMPL','SXO'))
	
	tmp	<- ans[, list(N=length(SUBTREE_ID)), by=c('SXO','MISSING_SMPL','INTRO_SMPL','REP','NL_TRMS_SXO_RND')]
	ans	<- subset(ans, MISSING_SMPL==1)[, list(N=length(SUBTREE_ID)), by=c('SXO','MISSING_SMPL','INTRO_SMPL','REP','SUBTREE_SIZE')]
	ans	<- rbind( melt(ans, measure.vars='SUBTREE_SIZE'), melt(tmp, measure.vars='NL_TRMS_SXO_RND') )
	
	ans[, VALUE_C:= ans[, cut(value, breaks=size.breaks, labels=size.labels)]]
	ans	<- ans[, list(N=sum(N)), by=c('SXO','MISSING_SMPL','INTRO_SMPL','REP','variable','VALUE_C')]
	tmp	<- ans[, list(VALUE_C=VALUE_C, P=N/sum(N)), by=c('SXO','MISSING_SMPL','INTRO_SMPL','REP','variable')]
	ans	<- merge(ans, tmp, by=c('SXO','MISSING_SMPL','INTRO_SMPL','REP','variable','VALUE_C'))
		
	tmp	<- subset(ans, variable=='SUBTREE_SIZE')
	total.replicates	<- tmp[, length(unique(MISSING_SMPL))*length(unique(INTRO_SMPL))*length(unique(REP))]
	tmp	<- tmp[, list(	P=ans.quantiles, 
						QN=quantile(c(N, rep(0, total.replicates-length(N))), p=ans.quantiles),
						QF=quantile(c(P, rep(0, total.replicates-length(P))), p=ans.quantiles)
						), by=c('SXO','variable','VALUE_C')]
	ans	<- subset(ans, variable=='NL_TRMS_SXO_RND')
	total.replicates	<- ans[, length(unique(MISSING_SMPL))*length(unique(INTRO_SMPL))*length(unique(REP))]
	ans	<- ans[, list(	P=ans.quantiles, 
						QN=quantile(c(N, rep(0, total.replicates-length(N))), p=ans.quantiles),
						QF=quantile(c(P, rep(0, total.replicates-length(P))), p=ans.quantiles)
						), by=c('SXO','variable','VALUE_C')]
	ans	<- rbind(ans, tmp)
	
	ans
}

vli.estimate.subtree.distribution.by.subtype <- function(dss, resample.intro.n=1, sxo.select= c('MSM','hsx'), intro.year=2000, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), size.breaks=c(0,1,2,5,10,500), size.labels=c('1','2','3-5','6-10','>10'), seed=42)	
{
	#	resample.intro.n=1; sxo.select= c('MSM','hsx'); intro.year=2000; ans.quantiles=c(0.025,0.25,0.5,0.75,0.975); seed=42
	set.seed(seed)
	
	ans	<- subset(dss, select=c(ST, ASR_TYPE, REP, SUBTREE_ID, ANCESTRAL_STATE, MIN_SAMPLING_DATE, SUBTREE_SIZE))
	ans[, INTRO_SMPL:=1L]		
	ans	<- subset(ans, floor(MIN_SAMPLING_DATE)>=intro.year)
	
	tmp	<- copy(dss)
	tmp	<- melt(tmp, id.vars=c('SUBTREE_ID','SUBTREE_SIZE','REP'), measure.vars=colnames(tmp)[grepl('SUBTREE_SXO_',colnames(tmp))], variable.name='SXO')
	set(tmp, NULL, 'SXO', tmp[, gsub('SUBTREE_SXO_','',SXO)])
	tmp[, SXO_PROP:= value/SUBTREE_SIZE]
	tmp	<- subset(tmp, select=c(SUBTREE_ID,REP,SXO,SXO_PROP))
	tmp	<- subset(tmp, SXO%in%sxo.select & SXO_PROP>0)	
	ans	<- merge(ans, tmp, by=c('SUBTREE_ID','REP'),allow.cartesian=TRUE)	
	ans	<- ans[, list(N=length(SUBTREE_ID)), by=c('ST','SXO','REP','SUBTREE_SIZE')]
	
	
	ans[, SUBTREE_SIZE_C:= ans[, cut(SUBTREE_SIZE, breaks=size.breaks, labels=size.labels)]]
	ans	<- ans[, list(N=sum(N)), by=c('ST','SXO','REP','SUBTREE_SIZE_C')]
	tmp	<- ans[, list(SUBTREE_SIZE_C=SUBTREE_SIZE_C, P=N/sum(N)), by=c('ST','SXO','REP')]
	ans	<- merge(ans, tmp, by=c('ST','SXO','REP','SUBTREE_SIZE_C'))
	
	total.replicates	<- ans[, length(unique(REP))]
	ans	<- ans[, list(	P=ans.quantiles, 
						QN=quantile(c(N, rep(0, total.replicates-length(N))), p=ans.quantiles),
						QF=quantile(c(P, rep(0, total.replicates-length(P))), p=ans.quantiles)
			), by=c('ST','SXO','SUBTREE_SIZE_C')]
	ans	
}

vli.estimate.proportion.introductions.by.sexual.orientation <- function(df, sampling.prob=0.43, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), seed=42, bs.len=1e2)
{
	set.seed(seed)
	df	<- subset(ds, ASR_TYPE=='inf')
	
	dsa		<- subset(df, select=c(SUBTREE_ID, SUBTREE_SIZE, ANCESTRAL_STATE, REP))
	tmp		<- dsa[, list(BS=seq_len(bs.len), SUBTREE_UNOBS=rnbinom(bs.len, SUBTREE_SIZE, prob=sampling.prob)), by='SUBTREE_ID']
	dsa		<- merge(dsa, tmp, by='SUBTREE_ID')	
	tmp		<- melt(df, id.vars=c('SUBTREE_ID','SUBTREE_SIZE','REP'), measure.vars=colnames(df)[grepl('SUBTREE_SXO_',colnames(df))], variable.name='SXO')
	set(tmp, NULL, 'SXO', tmp[, gsub('SUBTREE_SXO_','',SXO)])
	tmp[, SXO_PROP:= value/SUBTREE_SIZE]
	tmp	<- subset(tmp, select=c(SUBTREE_ID,REP,SXO,SXO_PROP))
	dsa	<- merge(dsa, tmp, by=c('SUBTREE_ID','REP'),allow.cartesian=TRUE)	
	dsa[, NL_TRMS_SXO_WEIGHTED:= SXO_PROP*(SUBTREE_SIZE+SUBTREE_UNOBS-1)]
	setnames(dsa, 'SXO_PROP', 'ANCESTRAL_TRMS_SXO_WEIGHTED')
	ans	<- dsa[, list(	ANCESTRAL_TRMS_SXO_WEIGHTED=sum(ANCESTRAL_TRMS_SXO_WEIGHTED), 
						NL_TRMS_SXO_WEIGHTED=sum(NL_TRMS_SXO_WEIGHTED)), 
						by=c('SXO','BS','REP')]
	ans[, PROP_INTRODUCTIONS:= ANCESTRAL_TRMS_SXO_WEIGHTED/(ANCESTRAL_TRMS_SXO_WEIGHTED+NL_TRMS_SXO_WEIGHTED)]
	ans		<- ans[, list(P=ans.quantiles, Q=quantile(PROP_INTRODUCTIONS, p=ans.quantiles)), by='SXO']
	ans	
}

vli.asr.summary<- function(ph, da, dm, id.regex, select.descendants.regex, is.dated.tree)
{
	#	determine max parsimony state transitions into NL 
	setnames(da, colnames(da), gsub('\\.','_',toupper(colnames(da))))
	da	<- subset(da, HOSTS=='NL')
	
	#	determine time of nodes at either end of branch with state transition	
	tmp	<- da[, list(EDGE_ID= which(ph$edge[,2]==ROOT_NOS) ), by='UNIQUE_SPLITS']
	da	<- merge(da, tmp, by='UNIQUE_SPLITS')
	da[, BRL:= ph$edge.length[ EDGE_ID] ]	#this should be same as LENGTHS
	tmp	<- node.depth.edgelength(ph)
	da[, ROOT_HEIGHT:= max(tmp) - tmp[ ROOT_NOS ]]	# this is the height between the inner node and the most recent tip
	if(is.dated.tree)
	{
		tmp2	<- gsub(id.regex,'\\1',ph$tip.label[which.max(tmp)])	# this is the ID of the most recent tip
		tmp2	<- subset(dm, ID==tmp2)[, SAMPLING_DATE]					# this is the corresponding sampling date
		da[, ROOT_DATE:= tmp2-ROOT_HEIGHT]
		da[, PARENT_ROOT_DATE:= ROOT_DATE-LENGTHS]
	}
	if(!is.dated.tree)
	{
		set(da, NULL, c('ROOT_DATE','PARENT_ROOT_DATE'), NA_real_)
	}
	
	#	determine number of tips per subtree
	da[, SUBTREE_SIZE:=sapply(Descendants(ph, da[, ROOT_NOS], type='tips'),function(x) length(which(grepl(select.descendants.regex, ph$tip.label[x]))))]
	
	#	determine composition of subtrees by sexual orientation 
	tmp	<- da[, list(TAXA=ph$tip.label[ unlist(Descendants(ph, ROOT_NOS, type='tips')) ] ), by='UNIQUE_SPLITS']
	tmp[, ID:=gsub(id.regex,'\\1',TAXA)]
	tmp	<- subset(tmp, grepl(select.descendants.regex,TAXA))	
	tmp	<- merge(dm, tmp, by='ID')	# we will need this below too
	setnames(tmp, 'SXO', 'L')
	tmp2<- as.data.table(expand.grid(UNIQUE_SPLITS=unique(tmp$UNIQUE_SPLITS), L=unique(subset(dm, SXO!='LANL', SXO))[, SXO] ))	
	tmp2<- merge(tmp2, tmp[, list(N=length(ID)), by=c('UNIQUE_SPLITS','L')], by=c('UNIQUE_SPLITS','L'), all.x=TRUE)
	set(tmp2, tmp2[, which(is.na(N))],'N',0L)	
	set(tmp2, NULL, 'L', tmp2[,paste0('SUBTREE_SXO_',L)])
	tmp2<- dcast.data.table(tmp2, UNIQUE_SPLITS~L, value.var='N')
	da	<- merge(da, tmp2, by='UNIQUE_SPLITS')
	setnames(tmp, 'L', 'SXO')
	
	#	determine composition of subtrees by birth location
	setnames(tmp, 'BORN_LOC', 'L')
	tmp2<- as.data.table(expand.grid(UNIQUE_SPLITS=unique(tmp$UNIQUE_SPLITS), L=unique(subset(dm, select=BORN_LOC))[, BORN_LOC] ))	
	tmp2<- merge(tmp2, tmp[, list(N=length(ID)), by=c('UNIQUE_SPLITS','L')], by=c('UNIQUE_SPLITS','L'), all.x=TRUE)
	set(tmp2, tmp2[, which(is.na(N))],'N',0L)	
	set(tmp2, NULL, 'L', tmp2[,paste0('BORN_LOC_',tolower(L))])
	tmp2<- dcast.data.table(tmp2, UNIQUE_SPLITS~L, value.var='N')
	da	<- merge(da, tmp2, by='UNIQUE_SPLITS')
	setnames(tmp, 'L', 'BORN_LOC')
	
	#	clean up 
	set(da, NULL, c('PARENT_SPLITS','HOSTS','BRL','LENGTHS','ROOT_HEIGHT','ROOT_NOS','EDGE_ID'), NULL)
	setnames(da, c('UNIQUE_SPLITS','PARENT_HOSTS','ROOT_DATE','PARENT_ROOT_DATE'), c('SUBTREE_ID','ANCESTRAL_STATE','INTRO_LATEST','INTRO_EARLIEST')	)
	da
}

vli.bez.explore.ancestral.state.reconstruction.170828<- function()
{
	require(data.table)
	require(phangorn)
	require(gpplot2)
	
	
	infile.meta	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo_OR.csv'
	infile.asr	<- '~/Dropbox (SPH Imperial College)/scout_out_Matthew/workspace_F_withA1_lsd_infLoc.rda'
	outfile.base<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/170825_'
	load(infile.asr)
	dm	<- as.data.table(read.csv(infile.meta))
	set(dm, NULL, 'SXO', dm[, gsub('HTF|HTM','hsx',gsub('BL|CH|CHS|DU','other',as.character(GENDER_TRMGROUP)))])
		
	#	read variables from workspace
	ph	<- all.tree.info$only.tree$tree
	da	<- as.data.table(all.tree.info$only.tree$classification.results$collapsed)
	
	#	summarise introductions into NL
	id.regex	<- '^([^_]*).*'
	is.dated.tree	<- TRUE
	ds	<- vli.ancestral.state.reconstruction.summary(ph, da, dm, id.regex, is.dated.tree)	
	ds[, INTRO_MID:= (INTRO_EARLIEST+INTRO_LATEST)/2]
	setkey(ds, INTRO_MID)
	ds[, INTRO_ID:= seq_len(nrow(ds))]
	
	#	get meta data on tips that descend from viral introductions
	dti	<- vli.asr.get.tips.in.NL.subtrees(ph, da, dm, id.regex)
	dti	<- merge(unique(subset(ds, select=c(INTRO_ID, SUBTREE_ID))), dti, by='SUBTREE_ID')
	setkey(dti, INTRO_ID, SAMPLING_DATE)
	dti[, SAMPLING_YR:= floor(SAMPLING_DATE)]
	
			
	#	plot state transitions 
	#	these are defined by self reported location of infection, not length of branches
	#	the phylogenetics only gives the min/max date of the state transition
	ggplot(ds, aes(x=INTRO_MID, y=INTRO_ID, colour=ANCESTRAL_STATE)) +
			geom_point() +
			geom_errorbarh(aes(xmin=INTRO_EARLIEST, xmax=INTRO_LATEST, y=INTRO_ID)) +
			theme_bw() +
			labs(x='time of state transition\n(assumed to occur on phylogenetic branch)', y='branch with state transition into NL\n(ID)', colour='ancestral state\n(world region)')
	ggsave(file=paste0(outfile.base,'F_statetransitions.pdf'), w=6, h=4)
	
	#	plot composition of subtrees by sexual orientation
	dc	<- melt(ds, id.vars=c('INTRO_ID','INTRO_MID'), measure.vars=colnames(ds)[grepl('SUBTREE_SXO_',colnames(ds))])
	set(dc, NULL, 'variable', dc[, gsub('SUBTREE_SXO_','',variable)])
	ggplot(dc, aes(x=INTRO_ID, y=value, fill=variable)) +
			geom_bar(stat='identity') +
			theme_bw() +
			scale_y_continuous(expand=c(0,0)) +
			labs(x='subtree ID', y='SHM patients\ndescendant from introduction',fill='sexual orientation')
	ggsave(file=paste0(outfile.base,'F_subtree_composition_by_sxo.pdf'), w=10, h=6)
	
	#
	#	plot outbreak size from midpoint until now
	#
	dtc	<- dti[, list(N=length(ID)), by=c('INTRO_ID','SAMPLING_YR')]	
	tmp	<- as.data.table(expand.grid(INTRO_ID=unique(dtc$INTRO_ID), SAMPLING_YR=seq(min(dtc$SAMPLING_YR),max(dtc$SAMPLING_YR))))
	dtc	<- merge(tmp, dtc, by=c('INTRO_ID','SAMPLING_YR'), all.x=TRUE)
	set(dtc, dtc[, which(is.na(N))],'N',0L)
	setkey(dtc, INTRO_ID, SAMPLING_YR)
	tmp	<- dtc[, list(SAMPLING_YR=SAMPLING_YR, CUM_N=cumsum(N)), by='INTRO_ID']
	dtc	<- merge(dtc, tmp, c('INTRO_ID','SAMPLING_YR'))
	dtc	<- merge(dtc, subset(ds, select=c(INTRO_ID, ANCESTRAL_STATE)),by='INTRO_ID')	
	ggplot(dtc, aes(x=SAMPLING_YR, y=CUM_N, fill=factor(INTRO_ID))) + 
			geom_bar(stat='identity') +
			theme_bw() +		
			scale_y_continuous(expand=c(0,0)) +
			labs(x='year of diagnosis', y='SHM patients\ndescendant from introductions', fill='subtree ID')
	ggsave(file=paste0(outfile.base,'F_cumulative_cases_since_intro_by_subtree.pdf'), w=6, h=10)		
	ggplot(dtc, aes(x=SAMPLING_YR, y=CUM_N, fill=ANCESTRAL_STATE)) + 
			geom_bar(stat='identity') +
			theme_bw() +		
			scale_y_continuous(expand=c(0,0)) +
			labs(x='year of diagnosis', y='SHM patients\ndescendant from introductions', fill='ancestral state of introduction\n(world region)')
	ggsave(file=paste0(outfile.base,'F_cumulative_cases_since_intro_by_region.pdf'), w=6, h=10)
	
	
	#
	#	estimate proportion of viral introductions among all transmissions
	#
	sampling.prob	<- 0.43
	dp	<- vli.estimate.proportion.introductions.by.sexual.orientation(ds, sampling.prob=sampling.prob)
	set(dp, NULL, 'P', dp[, paste0('qu_',P)])
	dp	<- dcast.data.table(dp, SXO~P, value.var='Q')
	ggplot(dp, aes(x=SXO, fill=SXO, ymin=qu_0.025, lower=qu_0.25, middle=qu_0.5, upper=qu_0.75, ymax=qu_0.975)) +
			geom_boxplot(stat='identity') +
			theme_bw() +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.2), expand=c(0,0), limits=c(0,1)) +
			labs(x='sexual orientation', y='proportion of viral introductions\namong transmissions', fill='')
	ggsave(file=paste0(outfile.base,'F_estimated_proportion_of_viral_introductions.pdf'), w=8, h=5)
	
}

vli.bez.explore.ancestral.state.reconstruction.170901<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	outfile.base	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/170831_'
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/170831_asr_summaries.rda'
	load(infile)	
	set(ds, NULL, 'SUBTREE_ID', ds[,paste0(ST,'-',REP,'-',SUBTREE_ID)])
	
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/170831_asr_tipdata.rda'
	load(infile)	
	set(dti, NULL, 'SUBTREE_ID', dti[,paste0(ST,'-',REP,'-',SUBTREE_ID)])
	
	asr.type<- 'inf'	
	tmp	<- unique(subset(ds, ASR_TYPE==asr.type & REP==0, select=c(SUBTREE_ID, ANCESTRAL_STATE, INTRO_EARLIEST, INTRO_MID, INTRO_LATEST, SUBTREE_SIZE)))
	dp	<- merge(tmp, subset(dti, ASR_TYPE==asr.type & REP==0), by='SUBTREE_ID')
	dp	<- merge(dp, dp[, list(MIN_SAMPLING_DATE=min(SAMPLING_DATE), MAX_SAMPLING_DATE=max(SAMPLING_DATE)), by='SUBTREE_ID'], by='SUBTREE_ID')
	dp	<- subset(dp, MIN_SAMPLING_DATE>1995)
	#dp	<- subset(dp, ANCESTRAL_STATE=='EUW')
	#	order subtree_ids by intro_mid
	do	<- unique(subset(dp, select=c(SUBTREE_ID, ANCESTRAL_STATE, MIN_SAMPLING_DATE, MAX_SAMPLING_DATE)))
	setkey(do, ANCESTRAL_STATE, MIN_SAMPLING_DATE)
	do	<- do[, list(SUBTREE_ID=SUBTREE_ID, MIN_SAMPLING_DATE=MIN_SAMPLING_DATE, MAX_SAMPLING_DATE=MAX_SAMPLING_DATE, SUBTREE_ORDER=seq_len(length(SUBTREE_ID))), by='ANCESTRAL_STATE']
	dp	<- merge(dp, subset(do, select=c(SUBTREE_ID, SUBTREE_ORDER)), by='SUBTREE_ID')	
	tmp	<- unique(subset(ds, ASR_TYPE==asr.type & REP==0, select=c(SUBTREE_ID, SUBTREE_SIZE, SUBTREE_SXO_MSM, SUBTREE_SXO_U, SUBTREE_SXO_hsx, SUBTREE_SXO_other)))	
	tmp	<- melt(tmp, id.vars=c('SUBTREE_ID','SUBTREE_SIZE'), measure.vars=colnames(tmp)[grepl('SUBTREE_SXO_',colnames(tmp))], variable.name='SXO')	
	tmp[, SXO_PROP:= value/SUBTREE_SIZE]
	tmp	<- dcast.data.table(tmp, SUBTREE_ID+SUBTREE_SIZE~SXO, value.var='SXO_PROP')
	do	<- merge(tmp, do, by='SUBTREE_ID')	
	tmp	<- subset(melt(subset(do,SUBTREE_SIZE==1), measure.vars=colnames(do)[grepl('SUBTREE_SXO_',colnames(do))]), value>0)
	if(0)
	{
		ggplot(dp) +			
				geom_point(aes(x=SAMPLING_DATE, y=SUBTREE_ORDER), colour='grey50') +
				geom_scatterpie(data=subset(do,SUBTREE_SIZE>1), aes(x=MIN_SAMPLING_DATE, y=SUBTREE_ORDER, group=SUBTREE_ORDER, r=log(SUBTREE_SIZE)/5+0.2), cols=colnames(do)[grepl('SUBTREE_SXO_',colnames(do))]) +
				#geom_point(data=tmp, aes(x=MIN_SAMPLING_DATE, y=SUBTREE_ORDER, colour=variable)) +			
				theme_bw() +
				scale_x_continuous(breaks=seq(1990,2020,2)) +
				labs(x='', y='new in-country subtrees') +
				facet_grid(ANCESTRAL_STATE~., scales='free_y', space='free_y')
		ggsave(file=paste0(outfile.base,asr.type,'_ALLST_newtrees_over_time.pdf'),w=8,h=30)
		
	}
	if(1)
	{
		ggplot(dp) +			
				geom_segment(data=subset(do, SUBTREE_SIZE>1), aes(x=MIN_SAMPLING_DATE, xend=MAX_SAMPLING_DATE, y=SUBTREE_ORDER, yend=SUBTREE_ORDER), colour='black', size=1) +
				geom_point(aes(x=SAMPLING_DATE, y=SUBTREE_ORDER, colour=SXO)) +
				theme_bw() +
				theme(panel.spacing.y=unit(0.1, "lines"), strip.text.y=element_text(angle=0)) +
				scale_x_continuous(breaks=seq(1990,2020,2)) +			
				labs(x='', y='new in-country subtrees', fill='sexual\norientation') +
				facet_grid(ANCESTRAL_STATE~., scales='free_y', space='free_y')
		ggsave(file=paste0(outfile.base,asr.type,'_ALLST_newtrees_over_time.pdf'),w=8,h=16)		
	}
	if(1)
	{
		ggplot(dp) +			
				geom_segment(data=subset(do, SUBTREE_SIZE>1), aes(x=MIN_SAMPLING_DATE, xend=MAX_SAMPLING_DATE, y=SUBTREE_ORDER, yend=SUBTREE_ORDER), colour='black', size=1) +
				geom_point(aes(x=SAMPLING_DATE, y=SUBTREE_ORDER, colour=ST)) +
				theme_bw() +
				theme(panel.spacing.y=unit(0.1, "lines"), strip.text.y=element_text(angle=0)) +
				scale_x_continuous(breaks=seq(1990,2020,2)) +			
				labs(x='', y='new in-country subtrees', colour='HIV-1 subtype') +
				facet_grid(ANCESTRAL_STATE~., scales='free_y', space='free_y')
		ggsave(file=paste0(outfile.base,asr.type,'_ALLST_newtrees_over_time_colour_subtype.pdf'),w=8,h=16)		
	}
}


vli.bez.explore.ancestral.state.reconstruction.170831<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	outfile.base	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/170831_'
	infile	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/170831_asr_summaries.rda'
	load(infile)	
	set(ds, NULL, 'SUBTREE_ID', ds[,paste0(ST,'-',REP,'-',SUBTREE_ID)])
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/170831_asr_tipdata.rda'
	load(infile)	
	set(dti, NULL, 'SUBTREE_ID', dti[,paste0(ST,'-',REP,'-',SUBTREE_ID)])
	tmp			<- dti[, list(MIN_SAMPLING_DATE=min(SAMPLING_DATE)), by=c('SUBTREE_ID','ASR_TYPE')]
	ds			<- merge(ds, tmp, by=c('SUBTREE_ID','ASR_TYPE'))
	tmp			<- dti[, list( BORN_LOC_FIRST=BORN_LOC[sample(which.min(SAMPLING_DATE)[1],1)] ), by=c('SUBTREE_ID','ASR_TYPE')]
	ds			<- merge(ds, tmp, by=c('SUBTREE_ID','ASR_TYPE'))
	
	#	plot number of state transitions by //INTRO_MID// MIN_SAMPLING_DATE
	asr.type	<- 'inf'
	asr.type	<- 'sam'
	asr.type	<- 'bir'
	dss	<- subset(ds, ASR_TYPE==asr.type)	
	dss[, INTRO_MID_C:= cut(INTRO_MID, breaks=c(1930,1950,1970,1990,seq(1992,2014,2)))]
	dss[, MIN_SAMPLING_DATE_C:= cut(MIN_SAMPLING_DATE, breaks=c(1930,1950,1970,1990,seq(1992,2016,2)))]	
	tmp	<- subset(dss, REP==0)
	ggplot(tmp, aes(x=MIN_SAMPLING_DATE_C, fill=ANCESTRAL_STATE)) + 
			geom_bar() + 
			theme_bw() + 
			labs(x='earliest date of diagnosis\n(of first subject with a sequence that is part of the subtree)',y='new in-country subtrees', fill='geographic origin')
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_histogram_counts.pdf'), w=12, h=5)
	ggplot(tmp, aes(x=MIN_SAMPLING_DATE_C, fill=ANCESTRAL_STATE)) + 
			geom_bar(position='fill') + 
			scale_y_continuous(labels=scales:::percent, expand=c(0,0)) +
			theme_bw() + 
			labs(x='earliest date of diagnosis\n(of first subject with a sequence that is part of the subtree)',y='new in-country subtrees', fill='geographic origin')
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_histogram_pc.pdf'), w=12, h=5)	
	ggplot(tmp, aes(x=MIN_SAMPLING_DATE_C, fill=ANCESTRAL_STATE)) + 
			geom_bar() + 
			theme_bw() + 
			facet_wrap(~ST, ncol=2) +
			labs(x='earliest date of diagnosis\n(of first subject with a sequence that is part of the subtree)',y='new in-country subtrees', fill='geographic origin')
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_histogram_counts_by_st.pdf'), w=24, h=10)
	#
	#	with uncertainty on proportions and counts
	#	by looking at all 100 tree replicates
	#	
	dc	<- subset(dss, MIN_SAMPLING_DATE>=1995)
	dc	<- dc[, list(N=length(SUBTREE_ID)), by=c('REP','MIN_SAMPLING_DATE_C','ANCESTRAL_STATE')]
	dc	<- merge(dc, dc[, list(ANCESTRAL_STATE=ANCESTRAL_STATE, P=N/sum(N)), by=c('MIN_SAMPLING_DATE_C','REP')], by=c('MIN_SAMPLING_DATE_C','REP','ANCESTRAL_STATE'))
	dc	<- merge(as.data.table(expand.grid(MIN_SAMPLING_DATE_C=unique(dc$MIN_SAMPLING_DATE_C), REP=unique(dc$REP), ANCESTRAL_STATE=unique(dc$ANCESTRAL_STATE))), dc, by=c('MIN_SAMPLING_DATE_C','REP','ANCESTRAL_STATE'), all.x=TRUE)
	set(dc, dc[, which(is.na(N))], c('N','P'), 0)	
	#	numbers by ancestral state
	dc	<- subset(dss, MIN_SAMPLING_DATE>=1995)
	dc	<- dc[, list(N=length(SUBTREE_ID)), by=c('REP','ANCESTRAL_STATE')]
	dc	<- merge(dc, dc[, list(ANCESTRAL_STATE=ANCESTRAL_STATE, P=N/sum(N)), by=c('REP')], by=c('REP','ANCESTRAL_STATE'))
	dc	<- merge(as.data.table(expand.grid(REP=unique(dc$REP), ANCESTRAL_STATE=unique(dc$ANCESTRAL_STATE))), dc, by=c('REP','ANCESTRAL_STATE'), all.x=TRUE)
	set(dc, dc[, which(is.na(N))], c('N','P'), 0)	
	tmp	<- dc[, list(QU=quantile(N, prob=c(0.025,0.25,0.5,0.75,0.975)), P=paste0('p',c(0.025,0.25,0.5,0.75,0.975))), by=c('ANCESTRAL_STATE')]
	tmp	<- dcast.data.table(tmp, ANCESTRAL_STATE~P, value.var='QU')	
	ggplot(tmp, aes(x=ANCESTRAL_STATE, y=p0.5, fill=ANCESTRAL_STATE, ymin=p0.025, ymax=p0.975)) +
			geom_bar(stat='identity', position=position_dodge(width=0.8)) +
			geom_errorbar(position=position_dodge(width=0.8), width=0.5, colour="black")	+
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous() +
			labs(x='geographic origin',y='new in-country subtrees', fill='geographic origin')			
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_counts_origin_overall.pdf'), w=6, h=5)
	write.csv(tmp, file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_counts_origin_overall.csv'))
	#	proportions by ancestral state
	tmp	<- dc[, list(QU=quantile(P, prob=c(0.025,0.25,0.5,0.75,0.975)), P=paste0('p',c(0.025,0.25,0.5,0.75,0.975))), by=c('ANCESTRAL_STATE')]
	tmp	<- dcast.data.table(tmp, ANCESTRAL_STATE~P, value.var='QU')	
	ggplot(tmp, aes(x=ANCESTRAL_STATE, y=p0.5, fill=ANCESTRAL_STATE, ymin=p0.025, ymax=p0.975)) +
			geom_bar(stat='identity', position=position_dodge(width=0.8)) +
			geom_errorbar(position=position_dodge(width=0.8), width=0.5, colour="black")	+
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1), breaks=seq(0,1,0.2)) +
			labs(x='earliest date of diagnosis\n(of first subject with a sequence that is part of the subtree)',y='new in-country subtrees', fill='geographic origin')			
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_counts_origin_overall.pdf'), w=6, h=5)
	write.csv(tmp, file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_counts_origin_overall.csv'))
	#	numbers overall
	dc	<- subset(dss, MIN_SAMPLING_DATE>=1995)
	dc	<- dc[, list(N=length(SUBTREE_ID)), by=c('REP')]		
	set(dc, dc[, which(is.na(N))], c('N'), 0)	
	tmp	<- dc[, list(QU=quantile(N, prob=c(0.025,0.25,0.5,0.75,0.975)), P=paste0('p',c(0.025,0.25,0.5,0.75,0.975)))]
	write.csv(tmp, file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_counts_overall.csv'))
	#	numbers over time
	dc	<- subset(dss, MIN_SAMPLING_DATE>=1995)
	dc	<- dc[, list(N=length(SUBTREE_ID)), by=c('REP','MIN_SAMPLING_DATE_C')]
	dc	<- merge(as.data.table(expand.grid(REP=unique(dc$REP), MIN_SAMPLING_DATE_C=unique(dc$MIN_SAMPLING_DATE_C))), dc, by=c('REP','MIN_SAMPLING_DATE_C'), all.x=TRUE)
	set(dc, dc[, which(is.na(N))], c('N'), 0)	
	tmp	<- dc[, list(QU=quantile(N, prob=c(0.025,0.25,0.5,0.75,0.975)), P=paste0('p',c(0.025,0.25,0.5,0.75,0.975))), by=c('MIN_SAMPLING_DATE_C')]
	tmp	<- dcast.data.table(tmp, MIN_SAMPLING_DATE_C~P, value.var='QU')	
	ggplot(tmp, aes(x=MIN_SAMPLING_DATE_C, y=p0.5, ymin=p0.025, ymax=p0.975)) +
			geom_bar(stat='identity', colour='grey70') +
			geom_errorbar(colour="black",width=0.5)	+
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous() +
			labs(x='earliest date of diagnosis\n(of first subject with a sequence that is part of the subtree)',y='new in-country subtrees', fill='geographic origin')			
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_counts_time.pdf'), w=6, h=5)
	write.csv(tmp, file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_counts_time.csv'))
	#	proportions by ancestral state and time
	dc	<- subset(dss, MIN_SAMPLING_DATE>=1995)
	dc	<- dc[, list(N=length(SUBTREE_ID)), by=c('REP','MIN_SAMPLING_DATE_C','ANCESTRAL_STATE')]
	dc	<- merge(dc, dc[, list(ANCESTRAL_STATE=ANCESTRAL_STATE, P=N/sum(N)), by=c('REP','MIN_SAMPLING_DATE_C')], by=c('REP','MIN_SAMPLING_DATE_C','ANCESTRAL_STATE'))
	dc	<- merge(as.data.table(expand.grid(REP=unique(dc$REP), ANCESTRAL_STATE=unique(dc$ANCESTRAL_STATE), MIN_SAMPLING_DATE_C=unique(dc$MIN_SAMPLING_DATE_C))), dc, by=c('REP','ANCESTRAL_STATE','MIN_SAMPLING_DATE_C'), all.x=TRUE)
	set(dc, dc[, which(is.na(N))], c('N','P'), 0)	
	tmp	<- dc[, list(QU=quantile(P, prob=c(0.025,0.25,0.5,0.75,0.975)), P=paste0('p',c(0.025,0.25,0.5,0.75,0.975))), by=c('MIN_SAMPLING_DATE_C','ANCESTRAL_STATE')]
	tmp	<- dcast.data.table(tmp, ANCESTRAL_STATE+MIN_SAMPLING_DATE_C~P, value.var='QU')	
	tmp	<- subset(tmp, !MIN_SAMPLING_DATE_C%in%c('(1990,1992]','(1992,1994]','(1996,1998]','(1998,2000]'))
	ggplot(tmp, aes(x=MIN_SAMPLING_DATE_C, y=p0.5, fill=ANCESTRAL_STATE, ymin=p0.025, ymax=p0.975)) +
			geom_bar(stat='identity', position=position_dodge(width=0.8)) +
			geom_errorbar(position=position_dodge(width=0.8), width=0.5, colour="black")	+
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1), breaks=seq(0,1,0.2)) +
			labs(x='earliest date of diagnosis\n(of first subject with a sequence that is part of the subtree)',y='new in-country subtrees', fill='geographic origin')			
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_proportion_origin_over_time.pdf'), w=10, h=5)
	write.csv(tmp, file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_proportion_origin_over_time.csv'))
	
	#
	#	proportion of first cases from the Netherlands
	#
	dc	<- subset(dss, MIN_SAMPLING_DATE>=1995 & MIN_SAMPLING_DATE<2015)
	dc[, MIN_SAMPLING_DATE_C:= cut(MIN_SAMPLING_DATE, breaks=seq(1995,2015,5))]
	dc[, BORN_LOC_FIRST_NL:= as.character(factor(BORN_LOC_FIRST=='NL', levels=c(TRUE,FALSE),labels=c('Netherlands',"foreign-born")))]
	dc	<- dc[, list(N=length(SUBTREE_ID)), by=c('REP','MIN_SAMPLING_DATE_C','BORN_LOC_FIRST_NL')]	
	dc	<- merge(dc, dc[, list(BORN_LOC_FIRST_NL=BORN_LOC_FIRST_NL, P=N/sum(N)), by=c('REP','MIN_SAMPLING_DATE_C')], by=c('REP','MIN_SAMPLING_DATE_C','BORN_LOC_FIRST_NL'))
	dc	<- merge(as.data.table(expand.grid(REP=unique(dc$REP), BORN_LOC_FIRST_NL=unique(dc$BORN_LOC_FIRST_NL), MIN_SAMPLING_DATE_C=unique(dc$MIN_SAMPLING_DATE_C))), dc, by=c('REP','BORN_LOC_FIRST_NL','MIN_SAMPLING_DATE_C'), all.x=TRUE)
	set(dc, dc[, which(is.na(N))], c('N','P'), 0)	
	tmp	<- dc[, list(QU=quantile(N, prob=c(0.025,0.25,0.5,0.75,0.975)), P=paste0('p',c(0.025,0.25,0.5,0.75,0.975))), by=c('MIN_SAMPLING_DATE_C','BORN_LOC_FIRST_NL')]
	tmp	<- dcast.data.table(tmp, BORN_LOC_FIRST_NL+MIN_SAMPLING_DATE_C~P, value.var='QU')		
	ggplot(tmp, aes(x=MIN_SAMPLING_DATE_C, y=p0.5, fill=BORN_LOC_FIRST_NL, ymin=p0.025, ymax=p0.975)) +
			geom_bar(stat='identity', position=position_dodge(width=0.9)) +
			geom_errorbar(position=position_dodge(width=0.9), width=0.5, colour="black")	+
			theme_bw() + theme(legend.position='bottom') +			
			labs(	x='date of diagnosis\nof first subject in subtree',
					y='new in-country subtrees', 
					fill='location of birth\nof first subject in subtree')			
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_birthfirstsubject_counts_over_time.pdf'), w=6, h=5)
	write.csv(tmp, file=paste0(outfile.base, asr.type, '_ALLST_birthfirstsubject_counts_over_time.csv'))
	tmp	<- dc[, list(QU=quantile(P, prob=c(0.025,0.25,0.5,0.75,0.975)), P=paste0('p',c(0.025,0.25,0.5,0.75,0.975))), by=c('MIN_SAMPLING_DATE_C','BORN_LOC_FIRST_NL')]
	tmp	<- dcast.data.table(tmp, BORN_LOC_FIRST_NL+MIN_SAMPLING_DATE_C~P, value.var='QU')		
	ggplot(tmp, aes(x=MIN_SAMPLING_DATE_C, y=p0.5, fill=BORN_LOC_FIRST_NL, ymin=p0.025, ymax=p0.975)) +
			geom_bar(stat='identity', position=position_dodge(width=0.9)) +
			geom_errorbar(position=position_dodge(width=0.9), width=0.5, colour="black")	+
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1), breaks=seq(0,1,0.2)) +
			labs(	x='date of diagnosis\nof first subject in subtree',
					y='new in-country subtrees', 
					fill='location of birth\nof first subject in subtree')			
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_birthfirstsubject_proportion_over_time.pdf'), w=6, h=5)
	write.csv(tmp, file=paste0(outfile.base, asr.type, '_ALLST_birthfirstsubject_proportion_over_time.csv'))
	
	#
	#	subtree sizes, of subtrees introduced since 1995
	#
	asr.type	<- 'inf'; asr.label	<- 'subtrees by location of infection'
	asr.type	<- 'sam'; asr.label	<- 'subtrees by sampling location'	
	asr.type	<- 'bir'; asr.label	<- 'subtrees by location of birth'
	size.breaks	<- c(0,1,2,5,10,500); size.labels<- c('1','2','3-5','6-10','>10')	
	dss			<- subset(ds, ASR_TYPE==asr.type)
	df			<- vli.estimate.subtree.distribution.by.sexual.orientation(dss, sampling.prob=0.43, resample.intro.n=1, resample.missing.n=100, sxo.select= c('MSM','hsx'), intro.year=2000, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), size.breaks=size.breaks, size.labels=size.labels, seed=42)	
	set(df, NULL, 'P', df[, paste0('p',P)])
	tmp			<- dcast.data.table(df, SXO+variable+VALUE_C~P, value.var='QF')
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('SUBTREE_SIZE','NL_TRMS_SXO_RND'), labels=c(paste0('subtrees\nas observed in phylogeny'),paste0('subtrees\nadjusted for missing sequences')))])	
	ggplot(tmp, aes(x=VALUE_C, y=p0.5, ymin=p0.025, ymax=p0.975, fill=SXO)) + 
			geom_bar(position=position_dodge(), stat='identity') +	
			geom_errorbar(position=position_dodge(width=0.9), width=0.5, colour="black")	+
			facet_grid(~variable, scales='free') +
			scale_y_continuous(labels=scales:::percent) + 
			theme_bw() +
			labs(x='\nsubtree size', y=paste0(asr.label,'\n'), fill='sexual\norientation')
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_subtreesize_distribution.pdf'), w=9, h=4)
	write.csv(tmp, file=paste0(outfile.base, asr.type, '_ALLST_subtreesize_distribution.csv'))
	size.breaks	<- c(0,2,500); size.labels<- c('1-2','>2')
	dss			<- subset(ds, ASR_TYPE==asr.type)
	df			<- vli.estimate.subtree.distribution.by.sexual.orientation(dss, sampling.prob=0.43, resample.intro.n=1, resample.missing.n=100, sxo.select= c('MSM','hsx'), intro.year=1995, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), size.breaks=size.breaks, size.labels=size.labels, seed=42)	
	set(df, NULL, 'P', df[, paste0('p',P)])
	tmp			<- dcast.data.table(df, SXO+variable+VALUE_C~P, value.var='QF')
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('SUBTREE_SIZE','NL_TRMS_SXO_RND'), labels=c(paste0('subtrees\nas observed in phylogeny'),paste0('subtrees\nadjusted for missing sequences')))])	
	ggplot(tmp, aes(x=VALUE_C, y=p0.5, ymin=p0.025, ymax=p0.975, fill=SXO)) + 
			geom_bar(position=position_dodge(), stat='identity') +	
			geom_errorbar(position=position_dodge(width=0.9), width=0.5, colour="black")	+
			facet_grid(~variable, scales='free') +
			scale_y_continuous(labels=scales:::percent) + 
			theme_bw() +
			labs(x='\nsubtree size', y=paste0(asr.label,'\n'), fill='sexual\norientation')
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_subtreesize_distribution_breaks12.pdf'), w=9, h=4)
	write.csv(tmp, file=paste0(outfile.base, asr.type, '_ALLST_subtreesize_distribution_breaks12.csv'))
	
	#
	#	subtree sizes by subtype
	#	
	size.breaks	<- c(0,2,500); size.labels<- c('1-2','>2')
	dss			<- subset(ds, ASR_TYPE==asr.type)
	df			<- vli.estimate.subtree.distribution.by.subtype(dss, resample.intro.n=1, size.breaks=size.breaks, size.labels=size.labels, sxo.select= c('MSM','hsx'), intro.year=1995, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), seed=42)
	set(df, NULL, 'P', df[, paste0('p',P)])
	tmp			<- dcast.data.table(df, ST+SXO+SUBTREE_SIZE_C~P, value.var='QF')		
	ggplot(tmp, aes(x=SUBTREE_SIZE_C, y=p0.5, ymin=p0.025, ymax=p0.975, fill=ST)) + 
			geom_bar(position=position_dodge(), stat='identity') +	
			geom_errorbar(position=position_dodge(width=0.9), width=0.5, colour="black")	+
			facet_grid(~SXO, scales='free') +
			scale_y_continuous(labels=scales:::percent) + 
			theme_bw() +
			labs(x='\nsubtree size\nas observed in phylogeny', y=paste0(asr.label,'\n'), fill='sexual\norientation')
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_subtreesize_distribution_by_subtype_breaks12.pdf'), w=12, h=4)
	write.csv(tmp, file=paste0(outfile.base, asr.type, '_ALLST_subtreesize_distribution_by_subtype_breaks12.csv'))
	
	#
	#	plot composition of subtrees by sexual orientation
	#
	dc	<- subset(ds, REP==0 & ASR_TYPE==asr.type)
	dc	<- melt(dc, id.vars=c('SUBTREE_ID','INTRO_MID'), measure.vars=colnames(ds)[grepl('SUBTREE_SXO_',colnames(ds))])
	set(dc, NULL, 'variable', dc[, gsub('SUBTREE_SXO_','',variable)])
	ggplot(dc, aes(x=SUBTREE_ID, y=value, fill=variable)) +
			geom_bar(stat='identity') +
			theme_bw() + theme(legend.position='bottom', axis.text.x=element_blank()) +
			scale_y_continuous(expand=c(0,0)) +
			labs(x='subtree ID', y='diagnosed cases\ndescendant from new in-country HIV acquisition',fill='sexual orientation')
	ggsave(file=paste0(outfile.base, asr.type,'_ALLST_subtree_composition_by_sxo.pdf'), w=12, h=6)
	
}

vli.bez.summarize.tip.data<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	outfile.save	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/170831_asr_tipdata.rda'
	
	infile.meta	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo_OR.csv'
	dm	<- as.data.table(read.csv(infile.meta))
	set(dm, NULL, 'SXO', dm[, gsub('HTF|HTM','hsx',gsub('BL|CH|CHS|DU','other',as.character(GENDER_TRMGROUP)))])
	
	indir.asr	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results'
	infiles		<- data.table(F=list.files(indir.asr, pattern='^workspace', full.names=TRUE))
	infiles[, ASR_TYPE:= gsub('.*_annotated_([a-z]+).*','\\1',basename(F))]
	infiles[, REP:= as.integer(gsub('.*_([0-9]+)_annotated.*','\\1',basename(F)))]
	infiles[, ST:= gsub('^workspace_([^_]+)_.*','\\1',basename(F))]
	
	select.descendants.regex<- '.*_infNL_.*'
	id.regex	<- '^([^_]*).*'
	is.dated.tree	<- TRUE
	dti			<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results/workspace_A1_withC_lsd_000_annotated_bir.rda'
				cat('\nprocessing file ',F)
				load(F)
				#	read variables from workspace
				ph	<- all.tree.info$only.tree$tree
				da	<- as.data.table(all.tree.info$only.tree$classification.results$collapsed)	
				#	summarise introductions into NL
				dti	<- vli.asr.get.tips.in.NL.subtrees(ph, da, dm, id.regex, select.descendants.regex)
				dti
			}, by=c('ST','ASR_TYPE','REP')]	
	save(dti, file=outfile.save)
}

vli.bez.summarize.state.transitions<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	outfile.save	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/170831_asr_summaries.rda'
	
	infile.meta	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo_OR.csv'
	dm	<- as.data.table(read.csv(infile.meta))
	set(dm, NULL, 'SXO', dm[, gsub('HTF|HTM','hsx',gsub('BL|CH|CHS|DU','other',as.character(GENDER_TRMGROUP)))])
	
	indir.asr	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results'
	infiles		<- data.table(F=list.files(indir.asr, pattern='^workspace', full.names=TRUE))
	infiles[, ASR_TYPE:= gsub('.*_annotated_([a-z]+).*','\\1',basename(F))]
	infiles[, REP:= as.integer(gsub('.*_([0-9]+)_annotated.*','\\1',basename(F)))]
	infiles[, ST:= gsub('^workspace_([^_]+)_.*','\\1',basename(F))]
	
	select.descendants.regex<- '.*_infNL_.*'
	id.regex	<- '^([^_]*).*'
	is.dated.tree	<- TRUE
	ds			<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results/workspace_C_withD_lsd_000_annotated_inf.rda'
				cat('\nprocessing file ',F)
				load(F)
				#	read variables from workspace
				ph	<- all.tree.info$only.tree$tree
				da	<- as.data.table(all.tree.info$only.tree$classification.results$collapsed)	
				#	summarise introductions into NL
				ds	<- vli.asr.summary(ph, da, dm, id.regex, select.descendants.regex, is.dated.tree)	
				ds[, INTRO_MID:= (INTRO_EARLIEST+INTRO_LATEST)/2]
				setkey(ds, INTRO_MID)
				ds
			}, by=c('ST','ASR_TYPE','REP')]	
	save(ds, file=outfile.save)
}

######################################################################################

vli.bez.DataFile<- function()
{
	require(data.table)
	require(ggplot2)
	infile	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo.csv'
	indir	<- dirname(infile)
	basename<- file.path(indir,gsub('\\.csv','',basename(infile)))
	
	df		<- as.data.table(read.csv(infile))
	df1		<- subset(df, select=c(patient, sex, routesex, NL, subtype, date_LSD, SampleRegion_V1, COMBIREGIONDIST_V1, infectionregion_V1, country, country_infection, countryborn))
	set(df1, df1[, which(infectionregion_V1=='ZAz')], 'infectionregion_V1', 'ZAZ')
	set(df1, df1[, which(countryborn=='LA')], 'countryborn', 'LANL')
	set(df1, NULL, 'COMBIREGIONDIST_V1', df1[, gsub('_SHM|_LANL','',COMBIREGIONDIST_V1)])
	setnames(df1, 	c('patient','sex','routesex','NL','subtype','date_LSD','SampleRegion_V1','COMBIREGIONDIST_V1','infectionregion_V1','country','country_infection','countryborn'), 
			c('ID','GENDER','GENDER_TRMGROUP','DATA_FROM_SHM','SUBTYPE','SAMPLING_DATE','SAMPLING_LOC','BORN_LOC','INFECTED_LOC','SAMPLING_COUNTRY','INFECTED_COUNTRY','BORN_COUNTRY'))
	#	write to file
	#
	write.csv(df1, row.names=FALSE, file=gsub('\\.csv','_OR.csv',infile))	
	#
	#	plot born_loc vs infected_loc	
	tmp	<- df1[, list(N=length(ID)), by=c('BORN_LOC','INFECTED_LOC')]
	ggplot(tmp, aes(x=BORN_LOC, y=INFECTED_LOC, size=N)) + geom_point()
	ggsave(file=paste0(basename,'_plot_BORNLOC_vs_INFECTEDLOC.pdf'), h=6, w=6)
	ggplot(subset(tmp, BORN_LOC=='NL'), aes(x=INFECTED_LOC, y=N)) + geom_bar(stat='identity')
	ggsave(file=paste0(basename,'_plot_BORNNL_INFECTEDLOC_numbers.pdf'), h=5, w=8)
	ggplot(subset(tmp, BORN_LOC=='NL'), aes(x=INFECTED_LOC, y= N/sum(N))) + geom_bar(stat='identity') + scale_y_continuous(labels = scales:::percent, breaks=seq(0,1,0.1)) 
	ggsave(file=paste0(basename,'_plot_BORNNL_INFECTEDLOC_percent.pdf'), h=5, w=8)
	#
	#	plot sampling locations over time
	tmp	<- copy(df1)
	tmp[, SAMPLING_DATE2:= cut(SAMPLING_DATE, seq(1981,2016,5))]	
	ggplot(tmp, aes(x=SAMPLING_DATE2, fill=SAMPLING_LOC)) + geom_bar() + scale_fill_brewer(palette='Set3')
	ggsave(file=paste0(basename,'_plot_SAMPLINGDATE_by_SAMPLINGLOC_numbers.pdf'), h=5, w=8)
	ggplot(tmp, aes(x=SAMPLING_DATE2, fill=SAMPLING_LOC)) + geom_bar(position='fill') + scale_fill_brewer(palette='Set3') + scale_y_continuous(labels = scales:::percent, breaks=seq(0,1,0.1))
	ggsave(file=paste0(basename,'_plot_SAMPLINGDATE_by_SAMPLINGLOC_percent.pdf'), h=5, w=8)
}

######################################################################################

vli.bez.LSD<- function()
{	
	require(data.table)
	require(big.phylo)
	#
	#	write the overall dates file for NON-B subtypes
	#
	if(0)
	{
		infile.dates	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo_OR.csv'
		df				<- as.data.table(read.csv(infile.dates))
		df				<- subset(df, select=c(ID, SUBTYPE, SAMPLING_DATE))
		setnames(df, c('ID','SAMPLING_DATE'), c('TAXA','DATE'))
		tmp				<- copy(df)
		set(tmp, NULL, 'TAXA', tmp[,gsub('[0-9]$','',paste0(TAXA,'_subtype',SUBTYPE))])
		df				<- rbind(df, tmp)	
		infile.dates	<- gsub('\\.csv','_lsddates.csv',infile.dates)
		write.csv(df, row.names=FALSE, infile.dates)
		#
		#	write subtype specific dates files
		#
		indir.ft		<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/trees_ft'
		infiles			<- data.table(F=list.files(indir.ft, pattern='*newick$', full.names=TRUE))
		infiles[, DATES_FILE:= gsub('_ft_bs100\\.newick','_lsddates.csv',F)]
		invisible(infiles[, {
							cmd		<- cmd.lsd.dates(infile.dates, F, DATES_FILE, run.lsd=FALSE)
							system(cmd)
							NULL
						}, by='F'])
	}
	#
	#	write the overall dates file for subtype B
	#
	if(0)
	{
		infile.dates	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo_plusB.csv'
		df				<- as.data.table(read.csv(infile.dates))
		df				<- subset(df, select=c(ID, subtype, SAMPLING_DATE))
		setnames(df, c('ID','subtype','SAMPLING_DATE'), c('TAXA','SUBTYPE','DATE'))
		tmp				<- copy(df)
		set(tmp, NULL, 'TAXA', tmp[,gsub('[0-9]$','',paste0(TAXA,'_subtype',SUBTYPE))])
		df				<- rbind(df, tmp)	
		infile.dates	<- gsub('\\.csv','_lsddates.csv',infile.dates)
		write.csv(df, row.names=FALSE, infile.dates)
		#
		#	write subtype specific dates files
		#
		indir.ft		<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/trees_ft'
		infiles			<- data.table(F=list.files(indir.ft, pattern='polB.*newick$', full.names=TRUE))
		infiles[, DATES_FILE:= gsub('_ft_bs100\\.newick','_lsddates.csv',F)]
		invisible(infiles[, {
							cmd		<- cmd.lsd.dates(infile.dates, F, DATES_FILE, run.lsd=FALSE)
							system(cmd)
							NULL
						}, by='F'])
	}	
	
	#
	#	run LSD 
	#
	#indir.ft		<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/trees_ft_rerooted'
	#indir.dates		<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/trees_lsd'
	#outdir			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/trees_lsd'
	
	indir.ft		<- '/work/or105/ATHENA_2016/vlintros/trees_ft_rerooted'
	indir.dates		<- '/work/or105/ATHENA_2016/vlintros/trees_lsd'
	outdir			<- '/work/or105/ATHENA_2016/vlintros/trees_lsd'
	
	ali.len			<- 1000				
	lsd.args		<- '-v 1 -c -b 10'			# no rooting, no re-estimation of rates				
	#	get files	
	infiles	<- data.table(F=list.files(indir.ft, pattern='*newick$', full.names=TRUE, recursive=TRUE))
	infiles	<- subset(infiles, !grepl('RAxML',F))
	infiles[, SUBTYPE:= gsub('pol','',gsub('_.*','',basename(F)))]
	tmp		<- data.table(FD=list.files(indir.dates, pattern='*_lsddates.csv$', full.names=TRUE, recursive=TRUE))
	tmp[, SUBTYPE:= gsub('pol','',gsub('_.*','',basename(FD)))]
	infiles	<- merge(infiles, tmp, by='SUBTYPE')
	infiles[, FL:= file.path(outdir,gsub('_rr\\.|_rr_','_lsd_',gsub('\\.newick','',basename(F))))]	
	#	build LSD commands
	dlsd	<- infiles[, {
				cmd		<- cmd.lsd(F, FD, ali.len, outfile=FL, pr.args=lsd.args)
				list(CMD=cmd)
			}, by='F']
	tmp		<- dlsd[, paste(CMD, collapse='\n')]
	#	qsub
	cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=30, hpc.q="pqeelab", hpc.mem="5800mb",  hpc.nproc=1, hpc.load='module load R/3.3.3')							
	cmd		<- paste(cmd,tmp,sep='\n')
	#cat(cmd)					
	outfile	<- paste("scRAr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
	cmd.hpccaller(indir.ft, outfile, cmd)	
}

######################################################################################

vli.bez.addHXB2<- function()
{
	require(ape)
	require(data.table)
	outfile.hxb2	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/SequenceAlignments/hxb2.fasta'
	load('~/git/hivclust/pkg/data/refseq_hiv1_hxb2.rda')
	tmp				<- as.DNAbin(as.matrix(t(as.character(hxb2[2000:4000, HXB2.K03455]))))
	rownames(tmp)	<- 'HXB2'
	write.dna(tmp, file=outfile.hxb2, format='fa', colsep='', nbcol=-1)
	
	file.seq.m		<- '/Users/Oliver/Dropbox\\ \\(SPH\\ Imperial\\ College\\)/2017_NL_Introductions/SequenceAlignments/polB_SHM_LANL_reduced149_withD.fasta'
	outfile			<- '/Users/Oliver/Dropbox\\ \\(SPH\\ Imperial\\ College\\)/2017_NL_Introductions/SequenceAlignments/polB_SHM_LANL_reduced149_withD_withHXB2.fasta'
	outfile.hxb2	<- '/Users/Oliver/Dropbox\\ \\(SPH\\ Imperial\\ College\\)/2017_NL_Introductions/SequenceAlignments/hxb2.fasta'
	cmd				<- paste('mafft --anysymbol --add ',outfile.hxb2,' --auto ',file.seq.m,' > ', outfile, '', sep='')
	system(cmd)	
}

vli.bez.rmDRM<- function()
{
	require(ape)
	require(data.table)
	require(big.phylo)
	file				<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/SequenceAlignments/polB_SHM_LANL_reduced149_withD_withHXB2_e.fasta'
	outfile				<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/SequenceAlignments/polB_SHM_LANL_reduced149_withD_withHXB2_e_noDRM.fasta'
	s					<- read.dna(file, format='fa')
	tmp					<- which(grepl("HXB2",rownames(s)))
	rownames(s)[tmp]	<- 'HXB2'
	tmp					<- big.phylo:::seq.rm.drugresistance(s)
	nodr.info			<- tmp$nodr.info
	seq					<- tmp$nodr.seq
	write.dna(seq, file= outfile, format='fasta', colsep='', nbcol=-1)
	save(seq, nodr.info, file= gsub('fasta','rda',outfile))
}

vli.bez.FastTrees<- function()
{	
	require(big.phylo)
	#	run FastTree on HPC
	if(0)
	{
		#indir			<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle'
		indir			<- '/work/or105/ATHENA_2016/vlintros/trees_ft'
		infiles			<- list.files(indir, pattern='fasta$')
		for(infile in infiles)
		{
			infile.fasta	<- file.path(indir,infile)
			bs.dir			<- gsub('.fasta','_bootstrap_trees',infile.fasta)
			bs.n			<- 100	
			dir.create(bs.dir)	
			outfile.ft		<- gsub('\\.fasta',paste0('_ft_bs',bs.n,'.newick'),infile.fasta)
			tmp				<- cmd.fasttree.many.bootstraps(infile.fasta, bs.dir, bs.n, outfile.ft, pr.args='-nt -gtr -gamma', opt.bootstrap.by='nucleotide')
			
			#	run on HPC
			cmd				<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=998, hpc.q="pqeelab", hpc.mem="5800mb",  hpc.nproc=1, hpc.load='module load R/3.3.3')
			cmd				<- paste(cmd,tmp,sep='\n')
			cat(cmd)					
			outfile.cmd		<- paste("bez",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			cmd.hpccaller(indir, outfile.cmd, cmd)	
		}
	}
	#	re-root at random taxon with name 'subtree'
	if(1)
	{
		require(ape)
		require(adephylo)
		require(phytools)		
		#indir			<- '/Users/Oliver/Dropbox (SPH Imperial College)/2017_NL_Introductions/trees_ft'
		indir			<- '/work/or105/ATHENA_2016/vlintros/trees_ft'
		infiles			<- data.table(F=list.files(indir, pattern='newick$', recursive=TRUE, full.names=TRUE))
		#infiles			<- subset(infiles, grepl('cpx06', F))
		invisible(infiles[, {
							#F	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2017_NL_Introductions/trees_ft/06cpx_withD_bootstrap_trees/06cpx_withD_ft.000.newick'
							ph	<- read.newick(F)		
							tmp				<- which(grepl('subtype',ph$tip.label))[1]		
							ph				<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
							ph				<- ladderize(ph)	
							write.tree(ph, file=gsub('_ft\\.','_rr.',F))
							pdf(file=paste0(gsub('_ft\\.','_rr.',F),'.pdf'), w=20, h=10+Ntip(ph)/10)
							plot(ph, show.node.label=TRUE, cex=0.3)
							dev.off()
						}, by='F'])
		
	}
}
