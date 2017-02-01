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
	indir	<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples'
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
				tmp		<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPE_BFGSargs_s',par.s*100,'.rda'),basename(F)))
				save( fit, fci, pci, file=tmp)
			}, by='F']	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs2_s50.rda'),full.names=TRUE))
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep99_coalreg_using_TYPE.rda'
				load(F)
				list(	THETA= names(fit$bestfit$par),
						MLE= unname(fit$bestfit$par),
						PROF_MEDIAN=unname(apply(pci$sample, 2, median)))				
			}, by='F']
	res		<- melt(res, measure.vars=c('MLE','PROF_MEDIAN'), variable.name='STAT', value.name='V')
	ggplot(res, aes(x=THETA, y=V)) + 
			geom_violin() + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +			
			coord_cartesian(ylim=c(-7.5,7.5)) +
			scale_y_continuous(breaks=seq(-10,10,1)) +
			facet_grid(.~STAT) 
	ggsave(file= gsub('\\.rda','_violin.pdf',gsub('_rep','',infiles[1,F])), w=5,h=4)			
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
	set(res, tmp, 'TAXA_SAMPLED', res[tmp, as.numeric(gsub('.*_s([0-9]+).*','\\1',F))])
	res[, TR_COEFF:= NA_character_]
	set(res, res[, which(grepl('using_TYPEtrf|using_TYPE',F))],'TR_COEFF', 'I0 vs I1')
	res[, TSI_COEFF:= NA_character_]
	set(res, res[, which(grepl('using_TYPE_|using_TYPE\\.',F))],'TSI_COEFF', 'I0 vs I1')
	set(res, res[, which(grepl('ETFIaoi',F))],'TSI_COEFF', 'exact time since infection')
	set(res, res[, which(grepl('ETFIaoi',F) & grepl('_mb100',F) & grepl('_en100',F))],'TSI_COEFF', 'noisy time since infection\nno bias')
	set(res, res[, which(grepl('ETFIaoi',F) & grepl('_mb50',F) & !grepl('_en',F))],'TSI_COEFF', 'biased time since infection\nno noise')
	res[, BFGSargs:= grepl('BFGS',F)]
	set(res, NULL, 'THETA', res[, gsub('ETSI_NOISE','ETSI', THETA)])
	set(res, NULL, 'TSI_COEFF', res[, factor(TSI_COEFF, levels=c('I0 vs I1','exact time since infection','noisy time since infection\nno bias','biased time since infection\nno noise'))])
	#
	#
	res		<- melt(res, measure.vars=c('MLE','PROF_MEDIAN'), variable.name='STAT', value.name='V')
	#
	#	 compare MLEs
	#
	ggplot(subset(res, BFGSargs & STAT=='MLE'), aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +	
			coord_cartesian(ylim=c(-5,5)) +
			scale_y_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			facet_grid(N~TSI_COEFF) 
	ggsave(file= file.path(indir,'compare_MLEs.pdf'), w=10,h=7)
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
