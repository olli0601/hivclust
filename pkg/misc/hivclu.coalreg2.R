cr.master.ex3.dev.MCMC2<- function()
{
	require(coalreg)	
	require(data.table)
	require(ape)
	require(coda)
	
	par.s	<- 0.5	
	F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep1.nwk'
	#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep2.nwk'
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_50100.rda'),basename(F)))
	
	set.seed(42)
	phylo	<- read.tree( F )					
	phylo 	<- multi2di(drop.tip(phylo, sample(phylo$tip.label, replace=FALSE, size=length(phylo$tip.label)*par.s)), random=FALSE)
	#	create data.table with infection type
	phi		<- data.table(	TAXA=phylo$tip.label, 
							TYPE=sapply(strsplit(phylo$tip.label,'_',),'[[',2),
							ETSI=as.numeric(sapply(strsplit(phylo$tip.label,'_',),'[[',3)))
	set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
	#set(phi, NULL, 'ETSI', phi[,ETSI-mean(ETSI)])
	#	coreg
	X				<- data.matrix(phi[,2:3,with=FALSE])
	rownames(X)	<- phi[, TAXA]
	colnames(X)	<- colnames(phi)[2:3]
	# 	make sure baseline is coded as x=0, and that acute stage is coded x=1 
	#	hence exp(beta*x)=exp(beta) for individuals with x=1 and 1 for individuals with x=0
	X[, 'TYPE']	<- X[, 'TYPE']-1 	 
	fit				<- coreg.adaptiveMHwithinGibbs.or170710(	X, 
																phylo,
																transmission=~TYPE,
																infection=~ETSI,
																adapt.batch=function(accn){  ifelse(accn<1e3, 20, 50) },
																adapt.schedule=function(b){ 1+1*(1/b)^(1/3) }, 
																mhsteps=3e3,													
																maxHeight=10,
																maxNodeDepth=Inf,
																mincladesize=100,
																coef_logprior_sd=10,			
																hetInflation_logprior=NA,
																debug.nolkl=FALSE
																)	
	save(fit, file=outfile)	
	
	fit.mcmc	<- cbind(fit$trace, fit$trace_tr, fit$trace_inf)				
	if(is.na(par.hetInflation_logprior))
		fit.mcmc<- fit.mcmc[,-which(colnames(fit.mcmc)=='logHetInflation')]
	fit.mcmc	<- mcmc(fit.mcmc)
	pdf(file.path(dirname(outfile), gsub('\\.rda',paste0('_coeff_trace.pdf'),basename(outfile))), w=10, h=7)
	plot(fit.mcmc)
	dev.off()
	fit.mcmc.a	<- mcmc(fit$trace_qsd)
	pdf(file.path(dirname(outfile), gsub('\\.rda',paste0('_sdadapt.pdf'),basename(outfile))), w=5, h=10)
	plot(fit.mcmc.a)
	dev.off()
}

cr.master.ex3.dev.MCMC3<- function()
{
	require(coalreg)	
	require(data.table)
	require(ape)
	require(coda)
	
	par.s	<- 1	
	F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep30.nwk'
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_50100.rda'),basename(F)))
	
	set.seed(42)
	phylo	<- read.tree( F )
	if(par.s<1)
		phylo 	<- multi2di(drop.tip(phylo, sample(phylo$tip.label, replace=FALSE, size=length(phylo$tip.label)*par.s)), random=FALSE)
	#	create data.table with infection type
	phi		<- data.table(	TAXA=phylo$tip.label, 
			TYPE=sapply(strsplit(phylo$tip.label,'_',),'[[',2),
			ETSI=as.numeric(sapply(strsplit(phylo$tip.label,'_',),'[[',3)))
	set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
	#set(phi, NULL, 'ETSI', phi[,ETSI-mean(ETSI)])
	#	coreg
	X				<- data.matrix(phi[,2,with=FALSE])
	rownames(X)	<- phi[, TAXA]
	colnames(X)	<- colnames(phi)[2]
	# 	make sure baseline is coded as x=0, and that acute stage is coded x=1 
	#	hence exp(beta*x)=exp(beta) for individuals with x=1 and 1 for individuals with x=0
	X[, 'TYPE']	<- X[, 'TYPE']-1 	 
	
	
	n 	<- length( phylo$tip.label )
	sts <- setNames( node.depth.edgelength( phylo )[1:n] , phylo$tip.label )
	bdt <- DatedTree( phylo, sts )
	maxHeight <- 25
	mincladesize <- 100
	bdts <- .slice.and.stitch.DatedTrees( bdt, maxHeight, mincladesize )
	
	
	fit				<- coreg.adaptiveMHwithinGibbs.or170710(	X, 
			phylo,
			transmission=~TYPE,
			infection=~TYPE,
			adapt.batch=function(accn){  ifelse(accn<1e3, 20, 50) },
			adapt.schedule=function(b){ 1+1*(1/b)^(1/3) }, 
			mhsteps=3e3,													
			maxHeight=25,
			maxNodeDepth=Inf,
			mincladesize=100,
			coef_logprior_sd=10,			
			hetInflation_logprior=NA,
			debug.nolkl=FALSE
	)	
	save(fit, file=outfile)	
	
	fit.mcmc	<- cbind(fit$trace, fit$trace_tr, fit$trace_inf)				
	if(is.na(par.hetInflation_logprior))
		fit.mcmc<- fit.mcmc[,-which(colnames(fit.mcmc)=='logHetInflation')]
	fit.mcmc	<- mcmc(fit.mcmc)
	pdf(file.path(dirname(outfile), gsub('\\.rda',paste0('_coeff_trace.pdf'),basename(outfile))), w=10, h=7)
	plot(fit.mcmc)
	dev.off()
	fit.mcmc.a	<- mcmc(fit$trace_qsd)
	pdf(file.path(dirname(outfile), gsub('\\.rda',paste0('_sdadapt.pdf'),basename(outfile))), w=5, h=10)
	plot(fit.mcmc.a)
	dev.off()
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
	#infile.xml	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150.xml'
	infile.xml	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250.xml'
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
	#file	<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123.nex'
	file	<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123.nex'	 
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

cr.master.ex3.check.likelihoods.aoiModel<- function()
{
	require(coalregaoiModel1)
	require(data.table)
	require(ape)
	require(viridis)
	
	indir		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
	infile		<- file.path(indir, 'm3.RR5.n1250_seed123_rep5.nwk')	
	tr			<- read.tree( infile )
	tr			<- multi2di(ladderize(tr),random=FALSE)	
	tmp 		<- setNames(dist.nodes( tr )[(Ntip(tr)+1), seq_len(Ntip(tr))], tr$tip.label)
	bdt 		<- DatedTree( tr, tmp )	
	X 			<- data.frame( TYPE=grepl('I0', tr$tip.label) )
	rownames(X) <- tr$tip.label		
	X			<- data.matrix(X)	 
	fit 		<- trf.lasso(	bdt, 
								X, 
								trf_names='TYPE', 
								aoi_names='TYPE', 
								maxNodeDepth=Inf,
								maxHeight=10,
								lasso_threshold=5, 
								method='BFGS', 
								lnr0 = -2, 
								lnrLimits = c(-4, 2), 
								scale=FALSE)	
	outfile	<- file.path(dirname(infile), gsub('\\.nwk',paste0('_coalreg_trf_lasso.rda'),basename(infile)))
	save(fit, file=outfile)
}

cr.master.ex3.check.likelihoods.mincladesize100.maxHeight10<- function()
{
	require(coalreg)	
	require(data.table)
	require(ape)
	require(coda)
	
	indir		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'	
	coefs0		<- as.data.table(expand.grid( 	LOGSCALE=c(0.11, -1.04),
												LOGHET=0,
												COEF_INF=c(1.033, 1.29),
												DUMMY=1))	
	coefs0		<- coefs0[c(1,4),]
	coefs0[, FIXED_PARAMS:= c('mean MCMC','trf.lasso best.fit')]
	coef_trs	<- as.data.table(expand.grid( COEF_TR=c(1,log(5),2,2.5,3,3.5,4,5,6,7,8,10), 
					                          SAMPLE=c(1,0.8,0.5),
							  				  DUMMY=1))
	coef_trs	<- merge(coef_trs, coefs0, by='DUMMY',allow.cartesian=TRUE)
	tmp			<- data.table(DUMMY=1, REP=1:10) 
	coef_trs	<- merge(coef_trs, tmp, by='DUMMY',allow.cartesian=TRUE)	
	setkey(coef_trs, COEF_INF, SAMPLE, COEF_TR)
	infile		<- file.path(indir, 'm3.RR5.n1250_seed123_rep5.nwk')
	
	coef_trs	<- coef_trs[, {
				coef_tr		<- COEF_TR
				coef_inf	<- COEF_INF
				logscale	<- LOGSCALE
				loghet		<- LOGHET
				par.s		<- SAMPLE
				set.seed(REP)
				cat('\nprocess',coef_tr, ', ',par.s)	
				tryCatch({
							tr 			<- read.tree( infile )
							if(par.s<1)
								tr 		<- multi2di(drop.tip(tr, sample(tr$tip.label, replace=FALSE, size=length(tr$tip.label)*par.s)), random=FALSE)
							n			<- length( tr$tip.label )
							X 			<- data.frame( TYPE=grepl('I0', tr$tip.label) )
							rownames(X) <- tr$tip.label		
							sts 		<- setNames( node.depth.edgelength( tr )[1:n] , tr$tip.label )
							bdt 		<- DatedTree( tr, sts )
							X 			<- as.data.frame(X)
							maxHeight	<- 10
							maxNodeDepth<- Inf
							mincladesize<- 100
							transmission<- ~TYPE
							infection	<- ~TYPE
							bdts 		<- .slice.and.stitch.DatedTrees( bdt, maxHeight, mincladesize )		
							X_tr 		<- as.data.frame( model.matrix( formula( transmission ) , data = X) )
							X_inf 		<-  as.data.frame( model.matrix( formula( infection), data = X ) ) 
							X_tr_mat 	<- data.matrix( X_tr )
							X_inf_mat 	<- data.matrix( X_inf )
							names_coef_tr	<- colnames( X_tr )[-1]
							names_coef_inf	<- colnames( X_inf )[-1]	
							.h <- function( coef_tr, .bdt = bdt ){
								if (length( names_coef_tr ) == 0) {
									rv <-  ( rep(1, n ) )
								} else if (length( names_coef_tr ) == 1)  { 
									rv <- exp(  X_tr[.bdt$tip.label, names_coef_tr] * coef_tr) ;
								} else{
									rv <- exp( as.vector( X_tr_mat[.bdt$tip.label, names_coef_tr] %*% coef_tr ) )
								}
								setNames( rv / mean ( rv ) , .bdt$tip.label )
							}
							.s <- function( coef_inf, .bdt = bdt ){
								if (length( names_coef_inf ) == 0) {
									rv <- ( rep(1, n ) )
								} else if (length( names_coef_inf ) == 1) {
									rv <-   exp( X_inf[.bdt$tip.label, names_coef_inf] * coef_inf ) ; 
								} else {
									rv <- exp( as.vector( X_inf_mat[.bdt$tip.label, names_coef_inf] %*% coef_inf) )
								}
								setNames( rv / mean ( rv ) , .bdt$tip.label)
							}
							h0			<- .h(coef_tr )
							s0			<- .s(coef_inf )
							loglkl		<- sum(sapply( bdts, 
											function(.bdt ){
												.loglik(.bdt, 
														h0[.bdt$tip.label], 
														s0[.bdt$tip.label], 
														exp(logscale), 
														maxHeight = maxHeight, 
														maxNodeDepth = maxNodeDepth, 
														hetInflation = exp(loghet))
											}))
							
						}, error=function(e){ print(e$message); loglkl<<- NA_real_})								
				list(LKL=loglkl)
			}, by=c('COEF_TR','SAMPLE','COEF_INF','LOGSCALE','LOGHET','FIXED_PARAMS','REP')]
	coef_trs<- subset(coef_trs, !is.na(LKL))
	set(coef_trs, NULL, 'SAMPLE', coef_trs[,paste0('sampling prob\n',SAMPLE)])
	set(coef_trs, NULL, 'REP', coef_trs[,paste0('seed\n',REP)])
	outfile	<- file.path(dirname(infile), gsub('\\.nwk',paste0('_coalreg_likelihoods_for_mean_other_params.rda'),basename(infile)))
	save(infile, coef_trs, file=outfile)	
	
	ggplot(coef_trs, aes(x=COEF_TR, y=LKL, colour=FIXED_PARAMS, group=interaction(FIXED_PARAMS,REP))) + 
			geom_line() + geom_point(size=0.8) + 
			geom_vline(xintercept=log(5)) +
			theme_bw() + facet_wrap(SAMPLE~REP, scales='free', ncol=9)
	outfile	<- file.path(dirname(infile), gsub('\\.nwk',paste0('_coalreg_likelihoods_for_mean_other_params_1.pdf'),basename(infile)))
	ggsave(file=outfile, w=15, h=15)
}


cr.master.ex3.check.likelihoods.mincladesize200.maxHeight30<- function()
{
	require(coalreg)	
	require(data.table)
	require(ape)
	require(coda)
	
	indir		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'	
	coefs0		<- as.data.table(expand.grid( 	LOGSCALE=c(0.11, -1.04),
					LOGHET=0,
					COEF_INF=c(1.033, 1.29),
					DUMMY=1))	
	coefs0		<- coefs0[c(1,4),]
	coefs0[, FIXED_PARAMS:= c('mean MCMC','trf.lasso best.fit')]
	coef_trs	<- as.data.table(expand.grid( COEF_TR=c(1,log(5),2,2.5,3,3.5,4,5,6,7,8,10), 
					SAMPLE=c(1,0.8,0.5),
					DUMMY=1))
	coef_trs	<- merge(coef_trs, coefs0, by='DUMMY',allow.cartesian=TRUE)
	tmp			<- data.table(DUMMY=1, REP=1:10) 
	coef_trs	<- merge(coef_trs, tmp, by='DUMMY',allow.cartesian=TRUE)	
	setkey(coef_trs, COEF_INF, SAMPLE, COEF_TR)
	infile		<- file.path(indir, 'm3.RR5.n1250_seed123_rep5.nwk')
	
	coef_trs	<- coef_trs[, {
				#COEF_TR<-1; SAMPLE<- 0.5; LOGSCALE<- 0.11; LOGHET<- 0; COEF_INF<- 1.033; FIXED_PARAMS<- 'mean MCMC'; REP<-1
				coef_tr		<- COEF_TR
				coef_inf	<- COEF_INF
				logscale	<- LOGSCALE
				loghet		<- LOGHET
				par.s		<- SAMPLE

				set.seed(REP)
				cat('\nprocess',coef_tr, ', ',par.s)	
				tryCatch({
							tr 			<- read.tree( infile )
							if(par.s<1)
								tr 		<- multi2di(drop.tip(tr, sample(tr$tip.label, replace=FALSE, size=length(tr$tip.label)*par.s)), random=FALSE)
							n			<- length( tr$tip.label )
							X 			<- data.frame( TYPE=grepl('I0', tr$tip.label) )
							rownames(X) <- tr$tip.label		
							sts 		<- setNames( node.depth.edgelength( tr )[1:n] , tr$tip.label )
							bdt 		<- DatedTree( tr, sts )
							X 			<- as.data.frame(X)
							maxHeight	<- 30
							maxNodeDepth<- Inf
							mincladesize<- 200
							transmission<- ~TYPE
							infection	<- ~TYPE
							bdts 		<- .slice.and.stitch.DatedTrees( bdt, maxHeight, mincladesize )		
							X_tr 		<- as.data.frame( model.matrix( formula( transmission ) , data = X) )
							X_inf 		<-  as.data.frame( model.matrix( formula( infection), data = X ) ) 
							X_tr_mat 	<- data.matrix( X_tr )
							X_inf_mat 	<- data.matrix( X_inf )
							names_coef_tr	<- colnames( X_tr )[-1]
							names_coef_inf	<- colnames( X_inf )[-1]	
							.h <- function( coef_tr, .bdt = bdt ){
								if (length( names_coef_tr ) == 0) {
									rv <-  ( rep(1, n ) )
								} else if (length( names_coef_tr ) == 1)  { 
									rv <- exp(  X_tr[.bdt$tip.label, names_coef_tr] * coef_tr) ;
								} else{
									rv <- exp( as.vector( X_tr_mat[.bdt$tip.label, names_coef_tr] %*% coef_tr ) )
								}
								setNames( rv / mean ( rv ) , .bdt$tip.label )
							}
							.s <- function( coef_inf, .bdt = bdt ){
								if (length( names_coef_inf ) == 0) {
									rv <- ( rep(1, n ) )
								} else if (length( names_coef_inf ) == 1) {
									rv <-   exp( X_inf[.bdt$tip.label, names_coef_inf] * coef_inf ) ; 
								} else {
									rv <- exp( as.vector( X_inf_mat[.bdt$tip.label, names_coef_inf] %*% coef_inf) )
								}
								setNames( rv / mean ( rv ) , .bdt$tip.label)
							}
							h0			<- .h(coef_tr )
							s0			<- .s(coef_inf )
							loglkl		<- sum(sapply( bdts, 
											function(.bdt ){
												.loglik(.bdt, 
														h0[.bdt$tip.label], 
														s0[.bdt$tip.label], 
														exp(logscale), 
														maxHeight = maxHeight, 
														maxNodeDepth = maxNodeDepth, 
														hetInflation = exp(loghet))
											}))
							
						}, error=function(e){ print(e$message); loglkl<<- NA_real_})								
				list(LKL=loglkl)
			}, by=c('COEF_TR','SAMPLE','COEF_INF','LOGSCALE','LOGHET','FIXED_PARAMS','REP')]
	coef_trs<- subset(coef_trs, !is.na(LKL))
	set(coef_trs, NULL, 'SAMPLE', coef_trs[,paste0('sampling prob\n',SAMPLE)])
	set(coef_trs, NULL, 'REP', coef_trs[,paste0('seed\n',REP)])
	outfile	<- file.path(dirname(infile), gsub('\\.nwk',paste0('_coalreg_likelihoods_for_mean_other_params_mincladesize200_maxHeight30.rda'),basename(infile)))
	save(infile, coef_trs, file=outfile)	
	
	ggplot(coef_trs, aes(x=COEF_TR, y=LKL, colour=FIXED_PARAMS, group=interaction(FIXED_PARAMS,REP))) + 
			geom_line() + geom_point(size=0.8) + 
			geom_vline(xintercept=log(5)) +
			theme_bw() + facet_wrap(SAMPLE~REP, scales='free', ncol=9)
	outfile	<- file.path(dirname(infile), gsub('\\.nwk',paste0('_coalreg_likelihoods_for_mean_other_params_mincladesize200_maxHeight30_1.pdf'),basename(infile)))
	ggsave(file=outfile, w=15, h=15)
}