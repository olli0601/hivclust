cr.various<- function()
{
	indir				<- '/work/or105/ATHENA_2016/master_examples'
	par.base.pattern	<- 'm3.RR5.n1250_seed123'	
	if(1)
	{
		par.s				<- 0.5
		cr.master.ex3.runcoalreg.using.TYPE.BFGS2(indir, par.base.pattern, par.s)	
	}	
	if(0)
	{
		par.s				<- 0.5
		cr.master.ex3.runcoalreg.using.TYPE.ETFI.vanilla.BFGS2(indir, par.base.pattern, par.s)	
	}
	if(0)
	{
		par.s				<- 0.5
		par.maxNodeDepth	<- 3
		cr.master.ex3.runcoalreg.using.TYPE.ETFI.vanilla.BFGS3(indir, par.base.pattern, par.s, par.maxNodeDepth)	
	}
	if(0)
	{
		par.s				<- 0.5
		par.tsimb			<- 0.5
		cr.master.ex3.runcoalreg.using.TYPE.ETFI.mbias.BFGS2(indir, par.base.pattern, par.s, par.tsimb)		
	}
	if(0)
	{
		par.s				<- 0.5
		par.tsimb			<- 1
		par.tsimn			<- 1
		cr.master.ex3.runcoalreg.using.TYPE.ETFI.lnnoise.BFGS2(indir, par.base.pattern, par.s, par.tsimb, par.tsimn)	
	}	
	if(0)
	{
		par.s				<- 0.5
		par.tsimb			<- 0.5
		par.tsimn			<- 1
		cr.master.ex3.runcoalreg.using.TYPE.ETFI.lnnoise.BFGS2(indir, par.base.pattern, par.s, par.tsimb, par.tsimn)	
	}			
}

cr.png.generate.data<- function()
{
	require(data.table)
	require(ape)
	require(phangorn)
	infiles	<- data.table(FR= c(	'~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43/150129_HPTN071_scHN_INTERNAL/150129_HPTN071_scHN_SIMULATED_INTERNAL.R',
									'~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov22.4-seed43/150129_HPTN071_scHN_INTERNAL/150129_HPTN071_scHN_SIMULATED_INTERNAL.R',
									'~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov41.3-seed43/150129_HPTN071_scHN_INTERNAL/150129_HPTN071_scHN_SIMULATED_INTERNAL.R'),
						  FT= c(	'~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43/150129_HPTN071_scHN_SIMULATED_TREE/150129_HPTN071_scHN_DATEDTREE.newick',
								  	'~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov22.4-seed43/150129_HPTN071_scHN_SIMULATED_TREE/150129_HPTN071_scHN_DATEDTREE.newick',
								  	'~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov41.3-seed43/150129_HPTN071_scHN_SIMULATED_TREE/150129_HPTN071_scHN_DATEDTREE.newick'))
	infiles[, SC:= basename(dirname(dirname(FR)))]
	outfile	<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-'
	t.end	<- 2020
	t.adult	<- 14
	t.early	<- 1
	t.acute	<- 3/12	
	#
	#	prepare tree files
	#
	infiles[, {
				#FR		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43/150129_HPTN071_scHN_INTERNAL/150129_HPTN071_scHN_SIMULATED_INTERNAL.R'
				#FT		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43/150129_HPTN071_scHN_SIMULATED_TREE/150129_HPTN071_scHN_DATEDTREE.newick'
				cat('\n',FR)
				load(FR)
				#	collect all meta data that we want
				tmp		<- subset(df.inds, !is.na(TIME_SEQ), c(IDPOP, TIME_TR, GENDER, DOB, RISK, DIAG_T, DIAG_CD4, ART1_T, TIME_SEQ))
				tmp[, DIAG_IN_RECENT:= as.character(factor(DIAG_T-TIME_TR<t.early, levels=c(TRUE,FALSE),labels=c('Y','N')))]
				tmp[, DIAG_IN_ACUTE:= as.character(factor(DIAG_T-TIME_TR<t.acute, levels=c(TRUE,FALSE),labels=c('Y','N')))]
				tmp[, ETSI:= TIME_SEQ-TIME_TR]
				tmp2	<- subset(df.trms, select=c(IDREC, IDTR, IDTR_TIME_INFECTED, SAMPLED_TR))
				setnames(tmp, 'IDPOP', 'IDREC')
				tmp		<- merge(tmp, tmp2, by='IDREC')
				tmp[, TRM_FROM_RECENT:= as.character(factor(TIME_TR-IDTR_TIME_INFECTED<t.early, levels=c(TRUE,FALSE),labels=c('Y','N')))]
				tmp[, TRM_FROM_ACUTE:= as.character(factor(TIME_TR-IDTR_TIME_INFECTED<t.acute, levels=c(TRUE,FALSE),labels=c('Y','N')))]
				tmp[, SAMPLED_TR:= as.character(factor(!is.na(SAMPLED_TR), levels=c(TRUE,FALSE),labels=c('Y','N')))]	
				tmp[, SEQ_COV_2020:= subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & (is.na(DOD) | DOD>t.end) & HIV=='Y' & TIME_TR<t.end)[, round(mean(!is.na(TIME_SEQ)),d=2)]]
				#	read tree and change tip labels
				cat('\n',FT)
				ph		<- read.tree(FT)
				phi		<- data.table(TAXA= ph$tip.label, DUMMY=seq_along(ph$tip.label))
				phi[, IDREC:= as.integer(gsub('^IDPOP_([0-9]+).*','\\1',TAXA))]
				phi		<- merge(phi, tmp, by='IDREC')
				set(phi, NULL, 'TAXA', phi[, paste(IDREC,GENDER,DOB,RISK,round(TIME_TR,d=3),round(DIAG_T,d=3),DIAG_CD4,DIAG_IN_RECENT,DIAG_IN_ACUTE,round(TIME_SEQ,d=3),round(ETSI,d=3),SAMPLED_TR,TRM_FROM_RECENT,TRM_FROM_ACUTE,SEQ_COV_2020,sep='|')])
				setkey(phi, DUMMY)
				ph$tip.label	<- phi[, TAXA]
				#	save new tree
				write.tree(ph, paste0(dirname(dirname(FT)),'.newick'))
			}, by='FR']	
	#
	#	calculate empirical risk ratios p(T=1 | X=x) / p(T=0 | X=x)  
	#	focus on:
	#	- risk group, with recent infection at diagnosis, transmission from recent, transmission from acute
	#	
	dfr	<- infiles[, {
			#FR		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43/150129_HPTN071_scHN_INTERNAL/150129_HPTN071_scHN_SIMULATED_INTERNAL.R'			
			cat('\n',FR)
			load(FR)
			#
			#	denominator population: all adults by 2020 that were infected before 2020
			dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR<t.end)
			#
			#	transmission rate given risk category: 
			#	what I can calculate is the #transmissions per year when in stage X
			#	so I need the duration of stage X and the #transmissions while in stage X
			#
			dfd[, DUR_RISK:= DOD]
			set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
			set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-(DOB+t.adult)])
			tmp		<- subset(df.trms, TIME_TR<t.end)[, list(TR_N=length(IDREC)), by='IDTR']
			setnames(tmp, 'IDTR', 'IDPOP')
			dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
			set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
			dfr		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='RISK'], measure.vars='RISK', variable.name='GROUP', value.name='FACTOR') 
			#	for early infection < 1 year vs late infection > 1 year
			#	(of course this will never be known, but include as ideal case)
			set(dfd, NULL, 'TR_N', NULL)
			tmp		<- as.data.table(expand.grid(TR_STAGE=c('early','late'),IDPOP=unique(dfd$IDPOP)))
			dfd		<- merge(dfd,tmp,by='IDPOP')			
			tmp		<- dfd[, which(TR_STAGE=='early')]
			set(dfd, tmp, 'DUR_RISK', dfd[tmp, pmin(DUR_RISK,rep(1,length(tmp)))])
			tmp		<- dfd[, which(TR_STAGE=='late')]
			set(dfd, tmp, 'DUR_RISK', dfd[tmp, pmax(DUR_RISK-1,rep(0,length(tmp)))])
			df.trms[, TR_STAGE:= as.character(factor(TIME_TR-IDTR_TIME_INFECTED<=t.early,levels=c(TRUE,FALSE),labels=c('early','late')))]
			tmp		<- subset(df.trms, TIME_TR<t.end)[, list(TR_N=length(IDREC)), by=c('IDTR','TR_STAGE')]
			df.trms[, TR_STAGE:=NULL]
			setnames(tmp, 'IDTR', 'IDPOP')
			dfd		<- merge(dfd, tmp, by=c('IDPOP','TR_STAGE'), all.x=1)
			set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
			#	marginal trm risk
			tmp		<- melt(subset(dfd, DUR_RISK>0)[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='TR_STAGE'], measure.vars='TR_STAGE', variable.name='GROUP', value.name='FACTOR')
			dfr		<- rbind(dfr, tmp)
			#	conditional trm risk for RISK group= M
			tmp		<- melt(subset(dfd, DUR_RISK>0 & RISK=='M')[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='TR_STAGE'], measure.vars='TR_STAGE', variable.name='GROUP', value.name='FACTOR')
			tmp[, GROUP:='TR_STAGE_RISKM']
			dfr		<- rbind(dfr, tmp)
			#	for diagnosed w CD4<350 vs diagnosed w CD4>350
			dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR<t.end & !is.na(DIAG_T))
			dfd[, DIAG_CD4_STAGE:= cut(DIAG_CD4, breaks=c(-1,350,1e4),labels=c('l350','g350'))]			
			dfd[, DUR_RISK:= DOD]
			set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
			set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-(DOB+t.adult)])
			tmp		<- subset(df.trms, TIME_TR<t.end)[, list(TR_N=length(IDREC)), by='IDTR']
			setnames(tmp, 'IDTR', 'IDPOP')
			dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
			set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
			tmp		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='DIAG_CD4_STAGE'], measure.vars='DIAG_CD4_STAGE', variable.name='GROUP', value.name='FACTOR') 
			dfr		<- rbind(dfr, tmp)
			tmp		<- melt(subset(dfd, RISK=='M')[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='DIAG_CD4_STAGE'], measure.vars='DIAG_CD4_STAGE', variable.name='GROUP', value.name='FACTOR')
			tmp[, GROUP:='DIAG_CD4_STAGE_RISKM']
			dfr		<- rbind(dfr, tmp)
			#	for diagnosed w recent infection vs late infection
			dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR<t.end & !is.na(DIAG_T))
			setnames(dfd, 'RECENT_TR','DIAG_IN_RECENT')
			dfd[, DUR_RISK:= DOD]
			set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
			set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-(DOB+t.adult)])
			tmp		<- subset(df.trms, TIME_TR<t.end)[, list(TR_N=length(IDREC)), by='IDTR']
			setnames(tmp, 'IDTR', 'IDPOP')
			dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
			set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
			tmp		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='DIAG_IN_RECENT'], measure.vars='DIAG_IN_RECENT', variable.name='GROUP', value.name='FACTOR') 
			dfr		<- rbind(dfr, tmp)
			tmp		<- melt(subset(dfd, RISK=='M')[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='DIAG_IN_RECENT'], measure.vars='DIAG_IN_RECENT', variable.name='GROUP', value.name='FACTOR')
			tmp[, GROUP:='DIAG_IN_RECENT_RISKM']
			dfr		<- rbind(dfr, tmp)
			dfr
		}, by='SC']
	save(dfr, file=paste0(outfile,'empirical_transmission_rates.rda'))	
}

cr.master.ex3.generate.data<- function()
{
	require(data.table)
	require(ape)
	require(hivclust)
	require(phangorn)
	#
	#	generate MASTER simulations based on Erik's XML file
	#
	#infile.xml	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150.xml'
	infile.xml	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250.xml'
	pr.beast2	<- '/Applications/BEAST 2.4.4/lib/beast.jar'
	pr.wdir		<- '/Users/Oliver/duke/tmp'
	pr.seed		<- 123
	cmd			<- hivc.cmd.beast2.runxml(	indir=dirname(infile.xml), 
											infile=basename(infile.xml), 
											outfile=paste0(gsub('\\.xml','',basename(infile.xml)),'_seed',pr.seed), 
											insignat=NULL, 
											prog.beast=pr.beast2, prog.wdir=pr.wdir, prog.seed=pr.seed, prog.opt.Xms="64m", prog.opt.Xmx="400m", hpc.ncpu=1)
	cat(cmd)
	#
	#	read out tip state info 
	#	
	#file	<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123.nex'
	file	<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123.nex'	 
	tmp		<- hivc.beast2out.read.nexus.and.stats(file, method.node.stat='any.node')
	phs		<- tmp$tree
	phi		<- copy(tmp$node.stat)
	set(phi, NULL, 'VALUE', phi[,gsub('"','',VALUE)])
	#
	#	write newick files with state (I0/I1) and exact infection time up to 4 digits
	#
	invisible(lapply(seq_along(phs), function(i)
			{
				#i		<- 1
				ph		<- phs[[i]]
				#	in MASTER simulation: 	all coalescent events are infection events
				#							so we only need the tip branch length
				stopifnot( !any(subset(phi, TREE_ID==names(phs)[i] & NODE_ID<=Ntip(ph) & STAT=='reaction')[, unique(VALUE)]%in%c('i00','i10')) )
				stopifnot( all(subset(phi, TREE_ID==names(phs)[i] & NODE_ID>Ntip(ph) & STAT=='reaction')[, unique(VALUE)]%in%c('i00','i10')) )				
				#	get state of individual at sampling time
				tmp		<- subset(phi, TREE_ID==names(phs)[i] & NODE_ID<=Ntip(ph) & STAT=='type', select=c(NODE_ID, STAT, VALUE))
				tmp		<- dcast.data.table(tmp, NODE_ID~STAT, value.var='VALUE')	
				#	get infection time of individual
				tmp2	<- data.table(	NODE_ID=seq_len(Ntip(ph)),
										TSI=sapply(seq_len(Ntip(ph)), function(x) ph$edge.length[which(ph$edge[,2]==x)]))
				tmp		<- merge(tmp, tmp2, by='NODE_ID')
				setnames(tmp, colnames(tmp), toupper(colnames(tmp)))
				set(tmp, NULL, 'TSI', tmp[, round(TSI, d=4)])
				setkey(tmp, NODE_ID)
				#	set taxon labels				
				ph$tip.label	<- tmp[, paste0(NODE_ID,'_',TYPE,'_',TSI)]
				write.tree(ph, file.path(dirname(file), gsub('\\.nex',paste0('_rep',i,'.nwk'),basename(file))))
			}))
}

cr.master.ex3.runcoalreg.using.TYPE.BFGS2<- function(indir, par.base.pattern, par.s)
{
	require(coalreg)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples'
		par.base.pattern	<- 'm3.RR5.n150_seed123'
		par.s				<- 0.5
	}
		
	#
	#	run coalreg	run 
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#	with 50% sampling
	#	
	set.seed(42)	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))
	infiles[, {
				#F		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep1.nwk'
				ph		<- read.tree( F )					
				ph 		<- drop.tip(ph, sample(ph$tip.label, replace=FALSE, size=length(ph$tip.label)*par.s))
				tmp 	<- setNames(dist.nodes( ph )[(Ntip(ph)+1), seq_len(Ntip(ph))], ph$tip.label)
				dph		<- DatedTree(ph, sampleTimes=tmp)
				#	create data.table with infection type
				phi		<- data.table(TAXA=ph$tip.label, TYPE=sapply(strsplit(ph$tip.label,'_',),'[[',2))
				set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
				rownames(phi)	<- phi[, TAXA]
				set(phi, NULL, 'TAXA', NULL)
				tmp		<- data.matrix(phi)
				fit 	<- trf.lasso(dph, 	tmp, trf_names = c( 'TYPE'), aoi_names = c( 'TYPE' ), lasso_threshold=5, method = 'BFGS',
											lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )	
				tmp		<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPE_BFGSargs2_s',par.s*100,'.rda'),basename(F)))
				save( fit, fci, pci, file=tmp)
			}, by='F']					
}

cr.master.ex3.runcoalreg.using.TYPE.ETFI.compare<- function()
{
	require(coalreg)
	require(viridis)
	indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples'
	infiles	<- data.table(F=list.files(indir, pattern='rda$',full.names=TRUE))
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep99_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
				load(F)
				list(	THETA= names(fit$bestfit$par),
						MLE= unname(fit$bestfit$par),
						MLE_CONVERGED= fit$bestfit$convergence==0,
						PROF_MEDIAN=unname(apply(pci$sample, 2, median)))				
			}, by='F']
	res[, REP:= as.integer(gsub('.*_rep([0-9]+)_.*','\\1',F))]
	res[, N:= as.integer(gsub('.*\\.n([0-9]+)_.*','\\1',F))]
	res[, TAXA_SAMPLED:= 1]
	tmp	<- res[, which(grepl('.*_s([0-9]+).*',F))]
	set(res, tmp, 'TAXA_SAMPLED', res[tmp, as.numeric(gsub('.*_s([0-9]+).*','\\1',F))/100])
	res[, TR_COEFF:= NA_character_]
	set(res, res[, which(grepl('using_TYPEtrf|using_TYPE',F))],'TR_COEFF', 'I0 vs I1')
	res[, TSI_COEFF:= NA_character_]
	set(res, res[, which(grepl('using_TYPE_|using_TYPE\\.',F))],'TSI_COEFF', 'I0 vs I1')
	set(res, res[, which(grepl('ETFIaoi',F))],'TSI_COEFF', 'exact time since infection')
	set(res, res[, which(grepl('ETFIaoi',F) & grepl('_mb100',F) & grepl('_en100|_ln100',F))],'TSI_COEFF', 'noisy time since infection\nno bias')
	set(res, res[, which(grepl('ETFIaoi',F) & grepl('_mb50',F) & !grepl('_en|_ln',F))],'TSI_COEFF', 'biased time since infection\nno noise')
	set(res, res[, which(grepl('ETFIaoi',F) & grepl('_mb50',F) & grepl('_en100|_ln100',F))],'TSI_COEFF', 'biased time since infection\nnoisy time since infection')
	res[, MAXNODE:=Inf]
	set(res, res[, which(grepl('BFGSargs3',F))], 'MAXNODE', 2)
	tmp	<- res[, which(grepl('maxNodeDepth[0-9]+',F))]
	set(res, tmp, 'MAXNODE', res[tmp,as.numeric(gsub('.*maxNodeDepth([0-9]+)_.*','\\1',F))])	
	res[, BFGSargs:= 'no constraints']
	set(res, res[, which(grepl('BFGS',F))], 'BFGSargs', 'r < -1')
	set(res, res[, which(grepl('BFGSargs2|BFGSargs3',F))], 'BFGSargs', 'r in -4,2')	
	res[, NOISE_MODEL:='none']
	set(res, res[, which(grepl('_en100',F))],'NOISE_MODEL', 'exp')
	set(res, res[, which(grepl('_ln100',F))],'NOISE_MODEL', 'lognormal')	
	set(res, NULL, 'THETA', res[, gsub('ETSI_NOISE','ETSI', THETA)])
	set(res, NULL, 'TSI_COEFF', res[, factor(TSI_COEFF, levels=c('I0 vs I1','exact time since infection','noisy time since infection\nno bias','biased time since infection\nno noise'))])
	#
	#
	resp	<- melt(res, measure.vars=c('MLE','PROF_MEDIAN'), variable.name='STAT', value.name='V')
	#
	#	 compare MLEs with constraint r< -1
	#
	subset(resp, BFGSargs=='r < -1' & STAT=='MLE' & TAXA_SAMPLED==.5 & NOISE_MODEL%in%c('exp','none') & THETA=='b_TYPE')[, table(N, TSI_COEFF)]
	ggplot(subset(resp, BFGSargs=='r < -1' & STAT=='MLE' & TAXA_SAMPLED==.5 & NOISE_MODEL%in%c('exp','none')), aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +	
			coord_cartesian(ylim=c(-5,5)) +
			scale_y_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			facet_grid(N~TSI_COEFF) 
	ggsave(file= file.path(indir,'compare_MLEs_BFGSargs1.pdf'), w=10,h=7)
	#
	#	compare MLEs with constraint r in -4,2
	#
	tmp	<- subset(resp, BFGSargs=='r in -4,2' & STAT=='MLE' & TAXA_SAMPLED==.5 & NOISE_MODEL%in%c('lognormal','none') & MAXNODE==Inf)
	subset(tmp,THETA=='b_TYPE')[, table(N, TSI_COEFF)]
	ggplot(tmp, aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +	
			coord_cartesian(ylim=c(-5,5)) +
			scale_y_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			facet_grid(N~TSI_COEFF) 
	ggsave(file= file.path(indir,'compare_MLEs_BFGSargs2.pdf'), w=10,h=7)	
	#
	#	compare maxNodeDepth results
	#
	tmp	<- subset(resp, BFGSargs=='r in -4,2' & grepl('BFGSargs3',F) & STAT=='MLE' & TAXA_SAMPLED==.5 & NOISE_MODEL%in%c('lognormal','none'))
	subset(tmp,THETA=='b_TYPE')[, table(N, MAXNODE)]
	tmp[, LABEL:=paste0('exact time since infection\nmaxNodeDepth\n',MAXNODE)]
	ggplot(tmp, aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +	
			coord_cartesian(ylim=c(-5,5)) +
			scale_y_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			facet_grid(N~LABEL) 
	ggsave(file= file.path(indir,'compare_MLEs_BFGSargs3_maxNodeDepth.pdf'), w=10,h=7)	
	

	#
	#	compare bias in MLE for trm risk
	#
	ggplot(subset(res, BFGSargs & STAT=='MLE' & THETA=='b_TYPE'), aes(x=TSI_COEFF, y=V-log(5))) +
			geom_hline(yintercept=0, colour='grey50', size=2) +
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='\nbias in \nlog risk ratio I1 vs baseline I0') +	
			coord_cartesian(ylim=c(-5,5)) +
			scale_y_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			facet_grid(~N) + 
			coord_flip()
	ggsave(file= file.path(indir,'compare_bias_trmrisk.pdf'), w=10,h=7)
	#
	#	compare BFGSargs3 to BFGSargs2
	#
subset(resp, BFGSargs=='r in -4,2' & STAT=='MLE' & TAXA_SAMPLED==.5 & NOISE_MODEL%in%c('lognormal','none') & THETA=='b_TYPE' & MAXNODE==Inf)[, table(N, TSI_COEFF)]

	tmp		<- subset(resp, BFGSargs%in%c('r < -1','r in -4,2','r in -4,2, maxDepth 2') & STAT=='MLE' & TSI_COEFF=='exact time since infection' & TAXA_SAMPLED==.5 & NOISE_MODEL%in%c('none'))
	subset(tmp, THETA=='b_TYPE')[, table(N, BFGSargs)]
	ggplot(tmp, aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +	
			coord_cartesian(ylim=c(-5,5)) +
			scale_y_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			facet_grid(N~BFGSargs) 
	
}


cr.master.ex3.runcoalreg.using.TYPE.ETFI.vanilla.BFGS1<- function(indir, par.base.pattern, par.s)
{
	require(coalreg)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples'			
		par.base.pattern	<- 'm3.RR5.n150_seed123'
		par.s				<- 0.5
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -3, lnrLimits = c(-Inf, -1), scale=F, lasso_threshold=5, method = 'BFGS'
	#	with 50% sampling	
	set.seed(42)	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))
	infiles[, {
				#F		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep1.nwk'
				ph		<- read.tree( F )					
				ph 		<- drop.tip(ph, sample(ph$tip.label, replace=FALSE, size=length(ph$tip.label)*par.s))
				tmp 	<- setNames(dist.nodes( ph )[(Ntip(ph)+1), seq_len(Ntip(ph))], ph$tip.label)
				dph		<- DatedTree(ph, sampleTimes=tmp)
				#	create data.table with infection type
				phi		<- data.table(	TAXA=ph$tip.label, 
						TYPE=sapply(strsplit(ph$tip.label,'_',),'[[',2),
						ETSI=as.numeric(sapply(strsplit(ph$tip.label,'_',),'[[',3)))
				set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
				#	zero mean ETSI so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI', phi[, ETSI-mean(ETSI)])				
				rownames(phi)	<- phi[, TAXA]
				set(phi, NULL, 'TAXA', NULL)
				tmp		<- data.matrix(phi)
				fit 	<- trf.lasso(dph, 	tmp, trf_names = c( 'TYPE'), aoi_names = c( 'ETSI' ), lasso_threshold=5, method = 'BFGS',
						lnr0 = -3, lnrLimits = c(-Inf, -1), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )	
				tmp		<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s',par.s*100,'.rda'),basename(F)))
				save( fit, fci, pci, file=tmp)
			}, by='F']		
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'),full.names=TRUE))
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep99_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
				load(F)
				list(	THETA= names(fit$bestfit$par),
						MLE= unname(fit$bestfit$par),
						PROF_MEDIAN=unname(apply(pci$sample, 2, median)))				
			}, by='F']
	res		<- melt(res, measure.vars=c('MLE','PROF_MEDIAN'), variable.name='STAT', value.name='V')
	ggplot(res, aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +			
			coord_cartesian(ylim=c(-7.5,7.5)) +
			scale_y_continuous(breaks=seq(-10,10,1)) +
			facet_grid(.~STAT) 
	ggsave(file= gsub('\\.rda','_violin.pdf',gsub('_rep','',infiles[1,F])), w=5,h=4)
	
}

cr.master.ex3.runcoalreg.using.TYPE.ETFI.vanilla.BFGS2<- function(indir, par.base.pattern, par.s)
{
	require(coalreg)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples'			
		par.base.pattern	<- 'm3.RR5.n150_seed123'	
		par.s				<- 0.5		
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#	with 50% sampling	 
	set.seed(42)	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))
	infiles[, {
				#F		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep1.nwk'
				ph		<- read.tree( F )					
				ph 		<- drop.tip(ph, sample(ph$tip.label, replace=FALSE, size=length(ph$tip.label)*par.s))
				tmp 	<- setNames(dist.nodes( ph )[(Ntip(ph)+1), seq_len(Ntip(ph))], ph$tip.label)
				dph		<- DatedTree(ph, sampleTimes=tmp)
				#	create data.table with infection type
				phi		<- data.table(	TAXA=ph$tip.label, 
										TYPE=sapply(strsplit(ph$tip.label,'_',),'[[',2),
										ETSI=as.numeric(sapply(strsplit(ph$tip.label,'_',),'[[',3)))
				set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
				#	zero mean ETSI so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI', phi[, ETSI-mean(ETSI)])				
				rownames(phi)	<- phi[, TAXA]
				set(phi, NULL, 'TAXA', NULL)
				tmp		<- data.matrix(phi)
				fit 	<- trf.lasso(dph, 	tmp, trf_names = c( 'TYPE'), aoi_names = c( 'ETSI' ), lasso_threshold=5, method = 'BFGS', lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )	
				tmp		<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs2_s',par.s*100,'.rda'),basename(F)))
				save( fit, fci, pci, file=tmp)
			}, by='F']		
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs2_s50.rda'),full.names=TRUE))
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep99_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
				load(F)
				list(	THETA= names(fit$bestfit$par),
						MLE= unname(fit$bestfit$par),
						PROF_MEDIAN=unname(apply(pci$sample, 2, median)))				
			}, by='F']
	res		<- melt(res, measure.vars=c('MLE','PROF_MEDIAN'), variable.name='STAT', value.name='V')
	ggplot(res, aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +			
			coord_cartesian(ylim=c(-7.5,7.5)) +
			scale_y_continuous(breaks=seq(-10,10,1)) +
			facet_grid(.~STAT) 
	ggsave(file= gsub('\\.rda','_violin.pdf',gsub('_rep','',infiles[1,F])), w=5,h=4)	
}

cr.png.runcoalreg.using.TRSTAGE.ETSI.vanilla.BFGS2<- function(indir, par.base.pattern, par.s)
{
	require(coalreg)
	require(data.table)
	require(ape)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_results'
		par.base.pattern	<- 'PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43'		
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#	with 50% sampling	 
	set.seed(42)	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'.*.newick'),full.names=TRUE))
	infiles[, {
				#F		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43.newick'
				ph		<- read.tree( F )
				ph		<- multi2di(ladderize(ph),random=FALSE)
				#	create data.table with infection type				
				#	paste(	IDREC,GENDER,DOB,RISK,round(TIME_TR,d=3),round(DIAG_T,d=3),DIAG_CD4,
				#			DIAG_IN_RECENT,DIAG_IN_ACUTE,round(TIME_SEQ,d=3),round(ETSI,d=3),SAMPLED_TR,TRM_FROM_RECENT,TRM_FROM_ACUTE,SEQ_COV_2020,sep='|')				
				phi		<- data.table(	TAXA=ph$tip.label,
										IDPOP=as.integer(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
										SEX=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',2),
										DOB=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',3)),
										RISK=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',4),
										TIME_TR=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',5)),
										DIAG_T=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',6)),
										DIAG_CD4=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',7)),										
										DIAG_IN_RECENT=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',8),
										DIAG_IN_ACUTE=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',9),
										TIME_SEQ=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',10)),
										ETSI=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',11)),
										SAMPLED_TR=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',12),
										TRM_FROM_RECENT=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',13),
										TRM_FROM_ACUTE=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',14),
										SEQ_COV_2020=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',15))	
										)				
				set(phi, NULL, 'TRM_FROM_RECENT', phi[,factor(TRM_FROM_RECENT,levels=c('N','Y'))])
				#	zero mean ETSI so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI', phi[, ETSI-mean(ETSI)])
				#	prepare coalreg input
				ph$tip.label	<- phi[, IDPOP]
				tmp				<- data.matrix(subset(phi, select=c(TRM_FROM_RECENT,ETSI)))
				rownames(tmp)	<- phi[, IDPOP]
				stopifnot(!any(is.na(tmp)))
				#	The rigorous way to choose lasso_threshold would be to use cross-validation which I do not have coded up yet; 
				#	in lieu of that I set it to approximately (number free parameters) * (maximum expected effect size). 
				#	You don't need to worry much about it for the BD sims, but for pangea it will help prevent over fitting.
				#		
				#	BFGS is fast, but i find is less robust than Nelder-Mead and gets stuck in local optima. 
				#	Feel free to experiment with that. 
				#		
				#	I think your suggested lnr limits are good. 
				#	Do you see that lnr is correlated with other parameter estimates ? 
				#	If outlier lnr's correspond to outlier parameter estimates, that is good cause to put limits on it. 
				#
				#	I set scale=F so that parameter estimates could be interpreted as log-odds for binary covariates.. 
				#	You will want this to be TRUE for PANGEA 			
				fit 	<- trf.lasso(ph, tmp, trf_names=c('TRM_FROM_RECENT'), aoi_names=c( 'ETSI' ), 
											lasso_threshold=5, method = 'BFGS', lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )	
				tmp		<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs2_s',par.s*100,'.rda'),basename(F)))
				save( fit, fci, pci, file=tmp)
			}, by='F']		
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs2_s50.rda'),full.names=TRUE))
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep99_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
				load(F)
				list(	THETA= names(fit$bestfit$par),
						MLE= unname(fit$bestfit$par),
						PROF_MEDIAN=unname(apply(pci$sample, 2, median)))				
			}, by='F']
	res		<- melt(res, measure.vars=c('MLE','PROF_MEDIAN'), variable.name='STAT', value.name='V')
	ggplot(res, aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +			
			coord_cartesian(ylim=c(-7.5,7.5)) +
			scale_y_continuous(breaks=seq(-10,10,1)) +
			facet_grid(.~STAT) 
	ggsave(file= gsub('\\.rda','_violin.pdf',gsub('_rep','',infiles[1,F])), w=5,h=4)	
}

cr.master.ex3.runcoalreg.using.TYPE.ETFI.vanilla.BFGS3<- function(indir, par.base.pattern, par.s, par.maxNodeDepth)
{
	require(coalreg)
	require(viridis)
	require(data.table)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples'			
		par.base.pattern	<- 'm3.RR5.n150_seed123'	
		par.maxNodeDepth	<- 2
		par.s				<- 0.5		
	}
	#
	#	run coalreg	run using exact time to infection
	#	with maxNodeDepth=2, extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#	with 50% sampling	 
	set.seed(42)	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))
	infiles[, {
				#F		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep1.nwk'
				ph		<- read.tree( F )					
				ph 		<- drop.tip(ph, sample(ph$tip.label, replace=FALSE, size=length(ph$tip.label)*par.s))
				tmp 	<- setNames(dist.nodes( ph )[(Ntip(ph)+1), seq_len(Ntip(ph))], ph$tip.label)
				dph		<- DatedTree(ph, sampleTimes=tmp)
				#	create data.table with infection type
				phi		<- data.table(	TAXA=ph$tip.label, 
						TYPE=sapply(strsplit(ph$tip.label,'_',),'[[',2),
						ETSI=as.numeric(sapply(strsplit(ph$tip.label,'_',),'[[',3)))
				set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
				#	zero mean ETSI so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI', phi[, ETSI-mean(ETSI)])				
				rownames(phi)	<- phi[, TAXA]
				set(phi, NULL, 'TAXA', NULL)
				tmp		<- data.matrix(phi)
				fit 	<- trf.lasso(dph, 	tmp, trf_names = c( 'TYPE'), aoi_names = c( 'ETSI' ), maxNodeDepth=par.maxNodeDepth, lasso_threshold=5, method = 'BFGS', lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )	
				tmp		<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs3_maxNodeDepth',par.maxNodeDepth,'_s',par.s*100,'.rda'),basename(F)))
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.master.ex3.runcoalreg.using.TYPE.ETFI.mbias.BFGS1<- function(indir, par.base.pattern, par.s, par.tsimb)
{
	require(coalreg)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples'
		par.base.pattern	<- 'm3.RR5.n150_seed123'
		par.s				<- 0.5
		par.tsimb			<- 0.5		
	}	
	#
	#	run coalreg	run using time to infection under multiplicative bias
	#	with extra args lnr0 = -3, lnrLimits = c(-Inf, -1), scale=F, lasso_threshold=5, method = 'BFGS'
	#	with 50% sampling
	set.seed(42)	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))
	infiles[, {
				#F		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11.nwk'
				ph		<- read.tree( F )					
				ph 		<- drop.tip(ph, sample(ph$tip.label, replace=FALSE, size=length(ph$tip.label)*par.s))
				tmp 	<- setNames(dist.nodes( ph )[(Ntip(ph)+1), seq_len(Ntip(ph))], ph$tip.label)
				dph		<- DatedTree(ph, sampleTimes=tmp)
				#	create data.table with infection type
				phi		<- data.table(	TAXA=ph$tip.label, 
						TYPE=sapply(strsplit(ph$tip.label,'_',),'[[',2),
						ETSI=as.numeric(sapply(strsplit(ph$tip.label,'_',),'[[',3)))
				set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
				#	multiplicative bias
				set(phi, NULL, 'ETSI', phi[, ETSI*par.tsimb])
				#	zero mean ETSI so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI', phi[, ETSI-mean(ETSI)])				
				rownames(phi)	<- phi[, TAXA]
				set(phi, NULL, 'TAXA', NULL)
				tmp		<- data.matrix(phi)
				fit 	<- trf.lasso(dph, 	tmp, trf_names = c( 'TYPE'), aoi_names = c( 'ETSI' ), lasso_threshold=5, method = 'BFGS',
						lnr0 = -3, lnrLimits = c(-Inf, -1), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )	
				tmp		<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s',par.s*100,'_mb',par.tsimb*100,'.rda'),basename(F)))
				save( ph, fit, fci, pci, file=tmp)
			}, by='F']	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s',par.s*100,'_mb',par.tsimb*100,'.rda'),full.names=TRUE))
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
				load(F)
				list(	THETA= names(fit$bestfit$par),
						MLE= unname(fit$bestfit$par),
						PROF_MEDIAN=unname(apply(pci$sample, 2, median)))				
			}, by='F']	
	res		<- melt(res, measure.vars=c('MLE','PROF_MEDIAN'), variable.name='STAT', value.name='V')
	ggplot(res, aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +			
			coord_cartesian(ylim=c(-7.5,7.5)) +
			scale_y_continuous(breaks=seq(-10,10,1)) +
			facet_grid(.~STAT) 
	ggsave(file= gsub('\\.rda','_violin.pdf',gsub('_rep','',infiles[1,F])), w=5,h=4)	
}

cr.master.ex3.runcoalreg.using.TYPE.ETFI.mbias.BFGS2<- function(indir, par.base.pattern, par.s, par.tsimb)
{
	require(coalreg)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples'
		par.base.pattern	<- 'm3.RR5.n150_seed123'
		par.s				<- 0.5
		par.tsimb			<- 0.5		
	}
	#
	#	run coalreg	run using time to infection under multiplicative bias
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#	with 50% sampling
	set.seed(42)	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))
	infiles[, {
				#F		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11.nwk'
				ph		<- read.tree( F )					
				ph 		<- drop.tip(ph, sample(ph$tip.label, replace=FALSE, size=length(ph$tip.label)*par.s))
				tmp 	<- setNames(dist.nodes( ph )[(Ntip(ph)+1), seq_len(Ntip(ph))], ph$tip.label)
				dph		<- DatedTree(ph, sampleTimes=tmp)
				#	create data.table with infection type
				phi		<- data.table(	TAXA=ph$tip.label, 
						TYPE=sapply(strsplit(ph$tip.label,'_',),'[[',2),
						ETSI=as.numeric(sapply(strsplit(ph$tip.label,'_',),'[[',3)))
				set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
				#	multiplicative bias
				set(phi, NULL, 'ETSI', phi[, ETSI*par.tsimb])
				#	zero mean ETSI so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI', phi[, ETSI-mean(ETSI)])				
				rownames(phi)	<- phi[, TAXA]
				set(phi, NULL, 'TAXA', NULL)
				tmp		<- data.matrix(phi)
				fit 	<- trf.lasso(dph, 	tmp, trf_names = c( 'TYPE'), aoi_names = c( 'ETSI' ), lasso_threshold=5, method = 'BFGS', lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )	
				tmp		<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs2_s',par.s*100,'_mb',par.tsimb*100,'.rda'),basename(F)))
				save( ph, fit, fci, pci, file=tmp)
			}, by='F']	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs2_s',par.s*100,'_mb',par.tsimb*100,'.rda'),full.names=TRUE))
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
				load(F)
				list(	THETA= names(fit$bestfit$par),
						MLE= unname(fit$bestfit$par),
						PROF_MEDIAN=unname(apply(pci$sample, 2, median)))				
			}, by='F']	
	res		<- melt(res, measure.vars=c('MLE','PROF_MEDIAN'), variable.name='STAT', value.name='V')
	ggplot(res, aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +			
			coord_cartesian(ylim=c(-7.5,7.5)) +
			scale_y_continuous(breaks=seq(-10,10,1)) +
			facet_grid(.~STAT) 
	ggsave(file= gsub('\\.rda','_violin.pdf',gsub('_rep','',infiles[1,F])), w=5,h=4)	
}

cr.master.ex3.runcoalreg.using.TYPE.ETFI.lnnoise.BFGS2<- function(indir, par.base.pattern, par.s, par.tsimb, par.tsimn)
{
	require(coalreg)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples'
		par.base.pattern	<- 'm3.RR5.n150_seed123'
		par.s				<- 0.5
		par.tsimb			<- 1
		par.tsimn			<- 1		
	}
	#
	#	run coalreg	run using time to infection under noise of 1 std dev
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#	with 50% sampling
	set.seed(42)	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))
	infiles[, {
				#F		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11.nwk'
				ph		<- read.tree( F )					
				ph 		<- drop.tip(ph, sample(ph$tip.label, replace=FALSE, size=length(ph$tip.label)*par.s))
				tmp 	<- setNames(dist.nodes( ph )[(Ntip(ph)+1), seq_len(Ntip(ph))], ph$tip.label)
				dph		<- DatedTree(ph, sampleTimes=tmp)
				#	create data.table with infection type
				phi		<- data.table(	TAXA=ph$tip.label, 
						TYPE=sapply(strsplit(ph$tip.label,'_',),'[[',2),
						ETSI=as.numeric(sapply(strsplit(ph$tip.label,'_',),'[[',3)))
				set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
				#	add noise: 	lognormal model with mean= ETSI and std dev= ETSI*par.tsimn
				tmp		<- phi[, list(ETSI_NOISE=rlnorm(1, meanlog= log(ETSI)-0.5*log(par.tsimn*par.tsimn+1), sdlog= sqrt(log(par.tsimn*par.tsimn+1)))), by='TAXA']				
				phi		<- merge(phi, tmp, by='TAXA')
				#	add multiplicative bias			
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE*par.tsimb])
				#	zero mean ETSI_NOISE so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE-mean(ETSI_NOISE)])				
				rownames(phi)	<- phi[, TAXA]
				set(phi, NULL, 'TAXA', NULL)
				tmp		<- data.matrix(phi)
				fit 	<- trf.lasso(dph, tmp, trf_names = c( 'TYPE'), aoi_names = c( 'ETSI_NOISE' ), lasso_threshold=5, method = 'BFGS',lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )	
				tmp		<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs2_s',par.s*100,'_mb',par.tsimb*100,'_ln',par.tsimn*100,'.rda'),basename(F)))
				save( ph, fit, fci, pci, file=tmp)
			}, by='F']	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs2_s',par.s*100,'_mb',par.tsimb*100,'_ln',par.tsimn*100,'.rda'),full.names=TRUE))
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
				load(F)
				list(	THETA= names(fit$bestfit$par),
						MLE= unname(fit$bestfit$par),
						PROF_MEDIAN=unname(apply(pci$sample, 2, median)))				
			}, by='F']	
	res		<- melt(res, measure.vars=c('MLE','PROF_MEDIAN'), variable.name='STAT', value.name='V')
	ggplot(res, aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +			
			coord_cartesian(ylim=c(-7.5,7.5)) +
			scale_y_continuous(breaks=seq(-10,10,1)) +
			facet_grid(.~STAT) 
	ggsave(file= gsub('\\.rda','_violin.pdf',gsub('_rep','',infiles[1,F])), w=5,h=4)	
}

cr.master.ex3.runcoalreg.using.TYPE.ETFI.expnoise.BFGS1<- function(indir, par.base.pattern, par.s, par.tsimb, par.tsimn)
{
	require(coalreg)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples'
		par.base.pattern	<- 'm3.RR5.n150_seed123'
		par.s				<- 0.5
		par.tsimb			<- 1
		par.tsimn			<- 1		
	}
	#
	#	run coalreg	run using time to infection under noise of 1 std dev
	#	with extra args lnr0 = -3, lnrLimits = c(-Inf, -1), scale=F, lasso_threshold=5, method = 'BFGS'
	#	with 50% sampling
	set.seed(42)	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))
	infiles[, {
				#F		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11.nwk'
				ph		<- read.tree( F )					
				ph 		<- drop.tip(ph, sample(ph$tip.label, replace=FALSE, size=length(ph$tip.label)*par.s))
				tmp 	<- setNames(dist.nodes( ph )[(Ntip(ph)+1), seq_len(Ntip(ph))], ph$tip.label)
				dph		<- DatedTree(ph, sampleTimes=tmp)
				#	create data.table with infection type
				phi		<- data.table(	TAXA=ph$tip.label, 
						TYPE=sapply(strsplit(ph$tip.label,'_',),'[[',2),
						ETSI=as.numeric(sapply(strsplit(ph$tip.label,'_',),'[[',3)))
				set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
				#	add noise: exponential model has mean==std dev, so:
				tmp		<- phi[, list(ETSI_NOISE=rexp(1,rate=1/(ETSI*par.tsimn))-ETSI), by='TAXA']	
				phi		<- merge(phi, tmp, by='TAXA')
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI+ETSI_NOISE])
				tmp		<- phi[, which(ETSI_NOISE<0)]
				set(phi, tmp, 'ETSI_NOISE', phi[tmp, ETSI])				
				#	add multiplicative bias			
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE*par.tsimb])
				#	zero mean ETSI_NOISE so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE-mean(ETSI_NOISE)])				
				rownames(phi)	<- phi[, TAXA]
				set(phi, NULL, 'TAXA', NULL)
				tmp		<- data.matrix(phi)
				fit 	<- trf.lasso(dph, 	tmp, trf_names = c( 'TYPE'), aoi_names = c( 'ETSI_NOISE' ), lasso_threshold=5, method = 'BFGS',lnr0 = -3, lnrLimits = c(-Inf, -1), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )	
				tmp		<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s',par.s*100,'_mb',par.tsimb*100,'_en',par.tsimn*100,'.rda'),basename(F)))
				save( ph, fit, fci, pci, file=tmp)
			}, by='F']	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s',par.s*100,'_mb',par.tsimb*100,'_en',par.tsimn*100,'.rda'),full.names=TRUE))
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
				load(F)
				list(	THETA= names(fit$bestfit$par),
						MLE= unname(fit$bestfit$par),
						PROF_MEDIAN=unname(apply(pci$sample, 2, median)))				
			}, by='F']	
	res		<- melt(res, measure.vars=c('MLE','PROF_MEDIAN'), variable.name='STAT', value.name='V')
	ggplot(res, aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +			
			coord_cartesian(ylim=c(-7.5,7.5)) +
			scale_y_continuous(breaks=seq(-10,10,1)) +
			facet_grid(.~STAT) 
	ggsave(file= gsub('\\.rda','_violin.pdf',gsub('_rep','',infiles[1,F])), w=5,h=4)	
}

cr.master.ex3.test<- function()
{
	require(coalreg)
	pkg.dir	<- '/Users/Oliver/git/coalreg'	
	X 		<- readRDS( file.path(pkg.dir,paste0('tests/', 'X-scenJ.rds')) ) 
	trees 	<- readRDS( file.path(pkg.dir,paste0('tests/', 'trees-scenJ.rds')) ) 
	
	fit 	<- trf.lasso( 	trees, 
							X,
							trf_names = c( 'AIDS', 'DIAG_CD4' , 'RISKH', 'RISKM', 'ISMALE'),
							aoi_names = c( 'AIDS', 'DIAG_CD4' ),
							lasso_threshold= 9, 
							method = 'Nelder-Mead'
							) 	
	fit_fisherci 	<-  fisher.ci( fit ) 			# ad hoc monte carlo
	fitci 			<- trf.confint( fit ) 			# ad hoc monte carlo	
	ci_aids 		<- prof.ci.vn2( fit, 'b_AIDS' )	# lik profile ci
}
