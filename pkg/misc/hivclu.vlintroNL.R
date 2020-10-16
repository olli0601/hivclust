######################################################################################

project.Bezemer.VLIntros<- function()
{
	vli.bez.LSD()
	#vli.bez.FastTrees()
}

vli.estimate.proportion.introductions.overall <- function(df, sampling.prob=0.43, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), seed=42, bs.len=1e3)
{
	set.seed(seed)	
	dsa		<- subset(df, select=c(SUBTREE_ID, SUBTREE_SIZE))
	tmp		<- dsa[, list(BS=seq_len(bs.len), SUBTREE_UNOBS=rnbinom(bs.len, SUBTREE_SIZE, prob=sampling.prob)), by='SUBTREE_ID']
	dsa		<- merge(dsa, tmp, by='SUBTREE_ID',allow.cartesian=TRUE)
	ans		<- dsa[, list(DUMMY=1, STAT='within_NL', N=sum(SUBTREE_SIZE)+sum(SUBTREE_UNOBS)-1), by='BS']
	tmp		<- df[, list(DUMMY=1, N=length(SUBTREE_ID)), by='ANCESTRAL_STATE']
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

vli.estimate.subtree.distribution.by.sexual.orientation <- function(df, sampling.prob=0.43, resample.intro.n=3, resample.missing.n=3, sxo.select= c('MSM','hsx'), size.breaks=c(0,1,2,5,10,500), size.labels=c('1','2','3-5','6-10','>10'), intro.year=2000, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), seed=42)	
{
	#	df<- copy(dpo); sampling.prob=0.43; resample.intro.n=3; resample.missing.n=3; sxo.select= c('MSM','hsx'); intro.year=2000; ans.quantiles=c(0.025,0.25,0.5,0.75,0.975); seed=42
	#	df<- copy(dpb); sampling.prob=0.43; resample.intro.n=1; resample.missing.n=100; sxo.select= c('MSM','hsx'); intro.year=1995; ans.quantiles=c(0.025,0.25,0.5,0.75,0.975); seed=42
	#	size.breaks=size.breaks, size.labels=size.labels, 
	set.seed(seed)
	df[, SUBTREE_ID2:= seq_len(nrow(df))]	
	setnames(df, c('SUBTREE_ID','SUBTREE_ID2'), c('SUBTREE_ID2','SUBTREE_ID'))
	
	ans	<- subset(df, select=c(ST, ASR_TYPE, REP, SUBTREE_ID, ANCESTRAL_STATE, MIN_SAMPLING_DATE, SUBTREE_SIZE))
	ans[, INTRO_SMPL:=1L]		
	ans	<- subset(ans, floor(MIN_SAMPLING_DATE)>=intro.year)
	tmp	<- ans[, list(MISSING_SMPL=seq_len(resample.missing.n), SUBTREE_UNOBS=rnbinom(resample.missing.n, SUBTREE_SIZE, prob=sampling.prob)), by=c('SUBTREE_ID','INTRO_SMPL')]
	ans	<- merge(ans, tmp, by=c('SUBTREE_ID','INTRO_SMPL'), allow.cartesian=TRUE)
		
	tmp	<- copy(df)
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
	
	ans[, VALUE_C:= ans[, cut(value, breaks=size.breaks, labels=size.labels, right=TRUE)]]
	ans	<- ans[, list(N=sum(N)), by=c('SXO','MISSING_SMPL','INTRO_SMPL','REP','variable','VALUE_C')]
	tmp	<- ans[, list(VALUE_C=VALUE_C, P=N/sum(N)), by=c('SXO','MISSING_SMPL','INTRO_SMPL','REP','variable')]
	ans	<- merge(ans, tmp, by=c('SXO','MISSING_SMPL','INTRO_SMPL','REP','variable','VALUE_C'))
		
	tmp	<- subset(ans, variable=='SUBTREE_SIZE')
	total.replicates	<- tmp[, length(unique(MISSING_SMPL))*length(unique(INTRO_SMPL))*length(unique(REP))]
	tmp2<- tmp[, list(	N=length(N) ), by=c('SXO','variable','VALUE_C')]
	stopifnot(all(tmp2$N==total.replicates))
	#	this is weird, not sure why I would need total.replicates-length(N)??
	tmp	<- tmp[, list(	P=ans.quantiles, 
						QN=quantile(c(N, rep(0, total.replicates-length(N))), p=ans.quantiles),
						QF=quantile(c(P, rep(0, total.replicates-length(P))), p=ans.quantiles)
						), by=c('SXO','variable','VALUE_C')]
	ans	<- subset(ans, variable=='NL_TRMS_SXO_RND')
	total.replicates	<- ans[, length(unique(MISSING_SMPL))*length(unique(INTRO_SMPL))*length(unique(REP))]
	tmp2<- ans[, list(	N=length(N) ), by=c('SXO','variable','VALUE_C')]
	stopifnot(all(tmp2$N==total.replicates))
	ans	<- ans[, list(	P=ans.quantiles, 
						QN=quantile(c(N, rep(0, total.replicates-length(N))), p=ans.quantiles),
						QF=quantile(c(P, rep(0, total.replicates-length(P))), p=ans.quantiles)
						), by=c('SXO','variable','VALUE_C')]
	ans	<- rbind(ans, tmp)	
	ans
}

vli.estimate.branchingprocess.incountryacquisition.200606 <- function()
{
	require(ggplot2)
	require(data.table)
	
	indir <- '~/Box/OR_Work/2017/2017_NL_Introductions/analysis'
	fname_origins <- '190801_subgraph_origins_N.csv'
	fnames <- c('190801_hsx_nonB_stananalysis.rda','190801_hsx_B_stananalysis.rda','190801_msm_nonB_stananalysis.rda','190801_msm_B_stananalysis.rda')
	fnames <- data.table(F= fnames)
	fnames[, ST:= gsub('^[0-9]+_([a-z]+)_([a-zA-Z]+)_.*$','\\2',F)]
	fnames[, RISK:= gsub('^[0-9]+_([a-z]+)_([a-zA-Z]+)_.*$','\\1',F)]
	fnames[, TRM_GROUP_SUB:= paste0(gsub('non','N',ST),'_NL',toupper(RISK))]
	
	do <- as.data.table( read.csv( file.path(indir, fname_origins), stringsAsFactor=FALSE) )
	do <- subset(do, PARENT_STATE!='Unknown')
	do[, PARENT_STATE2:= as.character(factor(grepl('^NL',PARENT_STATE), levels=c(TRUE,FALSE), labels=c('NL','EXT')))]
	do <- do[, list(N=sum(N)), by=c('TRM_GROUP_SUB','PARENT_STATE2')]
	do <- merge(do, do[, list(TOTAL=sum(N)), by='TRM_GROUP_SUB'], by='TRM_GROUP_SUB')
	
	ans<- list()	
	for(i in seq_len(nrow(fnames)))
	{
		fname <- fnames[i,F]
		load( file.path(indir,fname) )		
		tmp <- subset(do, TRM_GROUP_SUB==fnames$TRM_GROUP_SUB[i])
		tmp <- rbeta( length(sources_external), tmp$N[tmp$PARENT_STATE2=='EXT']+0.5, tmp$N[tmp$PARENT_STATE2=='NL']+0.5)		
		ans[[ i ]] <- data.table(	TRM_GROUP_SUB= fnames$TRM_GROUP_SUB[i],
									P_INTROS= sources_external, 
									P_NOINTROS= 1-sources_external, 
									P_EXT_ORIGIN= tmp,
									P_EXT_INTROS= sources_external * tmp,
									P_INCOUNTRY_ACQU= 1-sources_external * tmp) 
	}
	ans <- do.call('rbind',ans)
	ans <- melt(ans, id.vars='TRM_GROUP_SUB')
	ans <- ans[, list(V=quantile(value, p=c(0.025,0.5,0.975)), STAT=c('CL','M','CU')), by=c('TRM_GROUP_SUB','variable')]
	ans <- dcast.data.table(ans, TRM_GROUP_SUB+variable~STAT, value.var='V')
	ans[, L:= paste0( round(M*100, d=1), '% [',  round(CL*100, d=1),'-', round(CU*100, d=1),'%]')]
	ans[, XLAB:= factor(TRM_GROUP_SUB, levels=c('NB_NLHSX','NB_NLMSM','B_NLHSX','B_NLMSM'), c('heterosexuals\nnon-B subtypes','MSM\nnon-B subtypes','heterosexuals\nsubtype B','MSM\nsubtype B'))]
	ggplot(subset(ans, variable=='P_INCOUNTRY_ACQU'), aes(x=XLAB)) +
			geom_errorbar(aes(ymin=CL, ymax=CU), width=0.5) +
			geom_point(aes(y=M)) +			
			scale_y_continuous(labels=scales::percent, expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05)) +
			coord_cartesian(ylim=c(0,1)) +
			theme_bw() +
			labs(y='in-country HIV acquisitions\n(posterior median and 95% credibility intervals)\n', x='')
	ggsave(file=file.path(indir,'190801_incountryacquisitions.pdf'), w=6, h=6)
}

vli.estimate.branchingprocess.transmissionchains.200606 <- function()
{
	require(rstan)
	require(bayesplot)
	require(hexbin)
	
	args <- list()
	args$indir <- '~/Box/OR_Work/2017/2017_NL_Introductions/analysis'	
	args$outdir <- '~/Box/OR_Work/2017/2017_NL_Introductions/analysis'	
	args$infile.subgraphs <- file.path(args$indir, '190801_subgraph_data.csv')
	#	Dutch HSX non-B
	if(0)
	{
		ds <- as.data.table(read.csv(args$infile.subgraphs))
		ds <- subset(ds, REP==0 & TRM_GROUP=='NLHSX' & ST!='B')
		ds <- ds[, list(N=length(NAME)), by='SIZE']
		tmp <- vector('numeric', max(ds$SIZE))
		tmp[ds$SIZE] <- ds$N
		args$infile.subgraph.sizes <- '190801_subgraphsizes_hsx_nonB.csv'
		args$outfile.base <- '190801_hsx_nonB_'
		write.csv(tmp, file=file.path(args$indir,args$infile.subgraph.sizes))		
		args$size.inf.pop <- round( 6983*1402/(1402+1064), d=0 )
		args$size.inf.pop.sampled <- 1402					
	}
	#	Dutch HSX B
	if(0)
	{
		ds <- as.data.table(read.csv(args$infile.subgraphs))
		ds <- subset(ds, REP==0 & TRM_GROUP=='NLHSX' & ST=='B')
		ds <- ds[, list(N=length(NAME)), by='SIZE']
		ds <- ds[1:11,]
		tmp <- vector('numeric', max(ds$SIZE))
		tmp[ds$SIZE] <- ds$N
		args$infile.subgraph.sizes <- '190801_subgraphsizes_hsx_B.csv'
		args$outfile.base <- '190801_hsx_B_'
		write.csv(tmp, file=file.path(args$indir,args$infile.subgraph.sizes))		
		args$size.inf.pop <- round( 6983*1064/(1402+1064), d=0 )
		args$size.inf.pop.sampled <- 1064					
	}
	#	Dutch MSM B
	if(0)
	{
		ds <- as.data.table(read.csv(args$infile.subgraphs))
		ds <- subset(ds, REP==0 & TRM_GROUP=='NLMSM' & ST=='B')
		ds <- ds[, list(N=length(NAME)), by='SIZE']
		tmp <- vector('numeric', max(ds$SIZE))
		tmp[ds$SIZE] <- ds$N
		args$infile.subgraph.sizes <- '190801_subgraphsizes_msm_B.csv'
		args$outfile.base <- '190801_msm_B_'
		write.csv(tmp, file=file.path(args$indir,args$infile.subgraph.sizes))		
		args$size.inf.pop <- round( 13340*5113/(328+5113), d=0 )
		args$size.inf.pop.sampled <- 5113					
	}
	#	Dutch MSM non-B
	if(1)
	{
		ds <- as.data.table(read.csv(args$infile.subgraphs))
		ds <- subset(ds, REP==0 & TRM_GROUP=='NLMSM' & ST!='B')
		ds <- ds[, list(N=length(NAME)), by='SIZE']
		tmp <- vector('numeric', max(ds$SIZE))
		tmp[ds$SIZE] <- ds$N
		args$infile.subgraph.sizes <- '190801_subgraphsizes_msm_nonB.csv'
		args$outfile.base <- '190801_msm_nonB_'
		write.csv(tmp, file=file.path(args$indir,args$infile.subgraph.sizes))		
		args$size.inf.pop <- round( 13340*328/(328+5113), d=0 )
		args$size.inf.pop.sampled <- 328			
	}	
	args$upper.bound.multiplier <- 10
	
	stan.code <- "
			data{
			int<lower=1> N_cs_obs;							// max obs chain size
			int<lower=1> N_cs_actual;						// max size of actual chain
			row_vector<lower=0>[N_cs_obs] cs_obs;			// index i holds number of observed chains of size i
			real<lower=0> sampling_n;						// sampling total
			real<lower=0,upper=sampling_n> sampling_k;		// sampling success	
			}
			
			transformed data{
			//	get Binomial sampling successes 0,...,N_cs_obs.
			vector[N_cs_obs+1] bin_successes;	
			bin_successes[1] = 0;
			for(i in 1:N_cs_obs)
			{
			bin_successes[i+1] = bin_successes[i]+1; 
			}
			}
			
			parameters{
			real<lower=1e-10, upper=1> r0;	// R0
			real<lower=0> vmr_minus_one;	// variance to mean ratio of NegBin offspring distribution minus one, 1+R0/kappa-1= R0/kappa
			real<lower=0, upper=1> rho;		// sampling probability					
			}
			
			transformed parameters{
			real<lower=0> kappa;
			vector[N_cs_actual] cs_actual_lpmf;
			vector[N_cs_obs+1] cs_obs_lpmf;	
			
			// use local scoping to declare variables that don t need to be tracked
			{	
			real log_rho;
			real log_one_minus_rho;	
			real log_r0_div_kappa;
			real log_vmr;	
			matrix[N_cs_actual, 5] tmp;
			vector[5] ones;
			matrix[N_cs_obs+1, N_cs_actual+1] bin_lpmf;	
			int tmp_int;	// must declare index integer inside block
			
			// define transformed parameters
			kappa = r0/vmr_minus_one;
			log_one_minus_rho = log( 1 - rho );
			log_rho = log( rho );	
			log_r0_div_kappa = log( vmr_minus_one );
			log_vmr = log( 1+vmr_minus_one );
			
			// calculate lpmf of actual chain sizes given R0 and kappa
			ones = rep_vector(1., 5);
			for(i in 1:N_cs_actual)
			{
			tmp[i, 1]= kappa*i + i -1;
			tmp[i, 2]= kappa*i;
			tmp[i, 3]= i+1;
			tmp[i, 4]= i-1;
			tmp[i, 5]= kappa*i + i -1;
			}
			tmp[,1] = lgamma( tmp[,1] );
			tmp[,2] = -lgamma( tmp[,2] );
			tmp[,3] = -lgamma( tmp[,3] );
			tmp[,4] *= log_r0_div_kappa;
			tmp[,5] *= -log_vmr;
			cs_actual_lpmf = tmp * ones;		
			
			// calculate lpmf of sampling probabilities given rho		
			for(j in 1:(N_cs_actual+1))
			{		
			tmp_int= min(j, N_cs_obs+1);
			bin_lpmf[1:tmp_int,j] = bin_successes[1:tmp_int] * log_rho;
			bin_lpmf[1:tmp_int,j] += ( (rep_vector(j-1, tmp_int) - bin_successes[1:tmp_int]) * log_one_minus_rho );
			for(i in 1:tmp_int)
			{
			bin_lpmf[i,j] += lchoose( 1.*(j-1), bin_successes[i] );
			} 
			if(tmp_int<(N_cs_obs+1))
			{
			bin_lpmf[ (tmp_int+1):(N_cs_obs+1), j ] = rep_vector( negative_infinity(), N_cs_obs+1-tmp_int);
			}
			}
			
			//	calculate lpmf of observed chain sizes given R0 and kappa
			cs_obs_lpmf[1] = log_sum_exp( cs_actual_lpmf[1:N_cs_actual] + (bin_lpmf[1, 2:(N_cs_actual+1)])' );
			for(i in 1:N_cs_obs)
			{
			cs_obs_lpmf[i+1] = log_sum_exp( cs_actual_lpmf[i:N_cs_actual] + (bin_lpmf[i+1, (i+1):(N_cs_actual+1)])' );
			}
			
			//	renormalise conditional on 1 - prob nothing sampled
			cs_obs_lpmf[ 2:(N_cs_obs+1) ] -= log( 1-exp(cs_obs_lpmf[1]) );
			}		
			}
			
			model{
			// priors
			target+= beta_lpdf(r0 | 2, 2);
			target+= exponential_lpdf(vmr_minus_one | 1);
			target+= beta_lpdf(rho | sampling_k+0.5, sampling_n-sampling_k+0.5); 
			// likelihood
			target+= cs_obs * cs_obs_lpmf[ 2:(N_cs_obs+1) ];
			}
			"	
	
	stan.model <- stan_model(model_name= 'Blumberg2013', model_code = gsub('\t',' ',stan.code))
	
	infile.subgraph.sizes <- file.path(args$indir, args$infile.subgraph.sizes)
	outdir <- args$outdir
	
	#	debugging
	if(0)
	{
		stan.data <- list()		
		stan.data$cs_obs <- c(4,2,0,1)
		stan.data$N_cs_obs <- length(stan.data$cs_obs)
		stan.data$N_cs_actual <- 10 
		stan.data$sampling_n <- 100
		stan.data$sampling_k <- 80
		fit <- rstan::sampling(stan.model, data=stan.data, iter=10, warmup=5, 
				chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999),
				init= list(list(r0=0.5, vmr_minus_one=0.05, rho=0.5)))
		
	}
	
	stan.data <- list()
	tmp <- read.csv( infile.subgraph.sizes )
	colnames(tmp) <- c('SIZE','N')
	stan.data$cs_obs <- tmp$N
	stan.data$N_cs_obs <- length(stan.data$cs_obs)
	stan.data$N_cs_actual <- min(500,ceiling(stan.data$N_cs_obs / (args$size.inf.pop.sampled/args$size.inf.pop) * args$upper.bound.multiplier))	#	set upper bound for infinite sum approximation 
	stan.data$sampling_n <- args$size.inf.pop				# replace with actual number infected HSX Amsterdam
	stan.data$sampling_k <- args$size.inf.pop.sampled 		# replace with actual number sequenced infected HSX Amsterdam	
	fit <- rstan::sampling(stan.model, data=stan.data, iter=5e3, warmup=5e2, 
			chains=3, control = list(max_treedepth= 15, adapt_delta= 0.999),
			init= list(list(r0=0.5, vmr_minus_one=0.05, rho=0.5), list(r0=0.25, vmr_minus_one=0.05, rho=0.5), list(r0=0.75, vmr_minus_one=0.05, rho=0.5)))
	
	
	save(fit, file=file.path(outdir, paste0( args$outfile.base,'stanfit.rda' )))
	#	runs in 2 minutes
	
	
	#	examine neff and rhat
	fit.target.pars <- c('r0','vmr_minus_one','rho','kappa')
	rstan::monitor(rstan::extract(fit, pars=fit.target.pars, permuted = FALSE, inc_warmup = TRUE))
	#	neff too low, should ideally be 2e3
	
	#	traces	
	color_scheme_set("mix-blue-red")
	p <- rstan::traceplot(fit, pars=fit.target.pars, inc_warmup=TRUE, ncol = 1)
	pdf(file=file.path(outdir, paste0( args$outfile.base,'mcmc_traces.pdf' )), w=10, h=8)
	print(p)
	dev.off()
	
	#	keep only what we need to avoid memory meltdown
	fit.po <- rstan::extract(fit, permuted=TRUE, inc_warmup=FALSE)			
	
	#	pair plots	
	p <- mcmc_pairs(rstan::extract(fit, pars=fit.target.pars, permuted=FALSE, inc_warmup=FALSE), diag_fun = "dens", off_diag_fun = "hex")
	pdf(file=file.path(outdir, paste0( args$outfile.base,'mcmc_pairs.pdf' )), w=10, h=10)
	print(p)
	dev.off()
	
	# generate posterior predictive samples of actual chain sizes until we reach sampling_n individuals
	cs_actual_cdf <- exp(fit.po$cs_actual_lpmf)
	cs_actual_cdf <- t(apply( cs_actual_cdf, 1, cumsum))	
	samples.n <- 5e3
	cs_actual_postpred <- lapply( seq_len(nrow(cs_actual_cdf)), function(i){
				#IDX <<- i
				tmp <- runif(samples.n)
				#TMP <<- tmp
				cs_actual_postpred <- sapply(tmp, function(x){  head( which( x < c(cs_actual_cdf[i,],1) ), 1 )  })
				tmp2 <- which( args$size.inf.pop <= cumsum(cs_actual_postpred) )[1] 
				cs_actual_postpred[1:tmp2] 								
			})
		
	# posterior estimate of transmissions originating from outside Amsterdam
	sources_external <- sapply(cs_actual_postpred, function(x) length(x)/sum(x) )
	quantile(sources_external, p=c(0.5, 0.025, 0.975))
	save(sources_external, cs_actual_postpred, cs_actual_cdf, file=file.path(outdir, paste0( args$outfile.base,'stananalysis.rda' )))
	
	#	for non-B HSX:
	#	50%       2.5%      97.5% 
	#	0.5821159 0.5378998 0.6269218 

	#	for B HSX:
	#	      50%      2.5%     97.5% 
	#	0.5267176 0.4640007 0.5887040  

	#	for non-B MSM
	#	50%       2.5%      97.5%
	#	0.4054726 0.2311924 0.5581683 

	#	for B MSM
	#	50%       2.5%      97.5%
	#	0.2156988 0.1775749 0.2538293
}

vli.estimate.patients.in.small.subtrees.by.sexual.orientation <- function(df, sampling.prob=0.43, resample.intro.n=3, resample.missing.n=3, sxo.select= c('MSM','hsx'), size.breaks=c(0,1,2,5,10,500), size.labels=c('1','2','3-5','6-10','>10'), intro.year=2000, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), seed=42)	
{
	#	df<- copy(dpo); sampling.prob=0.43; resample.intro.n=3; resample.missing.n=3; sxo.select= c('MSM','hsx'); intro.year=2000; ans.quantiles=c(0.025,0.25,0.5,0.75,0.975); seed=42
	#	df<- copy(dpo); sampling.prob=0.43; resample.intro.n=1; resample.missing.n=3; sxo.select= c('MSM','hsx'); intro.year=2000; ans.quantiles=c(0.025,0.25,0.5,0.75,0.975); seed=42
	#	size.breaks=size.breaks, size.labels=size.labels, 
	set.seed(seed)
	df[, SUBTREE_ID2:= seq_len(nrow(df))]	
	setnames(df, c('SUBTREE_ID','SUBTREE_ID2'), c('SUBTREE_ID2','SUBTREE_ID'))
	
	ans	<- subset(df, select=c(ST, ASR_TYPE, REP, SUBTREE_ID, ANCESTRAL_STATE, MIN_SAMPLING_DATE, SUBTREE_SIZE))
	ans[, INTRO_SMPL:=1L]		
	ans	<- subset(ans, floor(MIN_SAMPLING_DATE)>=intro.year)
	tmp	<- ans[, list(MISSING_SMPL=seq_len(resample.missing.n), SUBTREE_UNOBS=rnbinom(resample.missing.n, SUBTREE_SIZE, prob=sampling.prob)), by=c('SUBTREE_ID','INTRO_SMPL')]
	ans	<- merge(ans, tmp, by=c('SUBTREE_ID','INTRO_SMPL'), allow.cartesian=TRUE)
	
	tmp	<- copy(df)
	tmp	<- melt(tmp, id.vars=c('SUBTREE_ID','SUBTREE_SIZE','REP'), measure.vars=colnames(tmp)[grepl('SUBTREE_SXO_',colnames(tmp))], variable.name='SXO', value.name='SXO_N')
	set(tmp, NULL, 'SXO', tmp[, gsub('SUBTREE_SXO_','',SXO)])
	tmp[, SXO_PROP:= SXO_N/SUBTREE_SIZE]	
	tmp	<- subset(tmp, select=c(SUBTREE_ID,REP,SXO,SXO_N,SXO_PROP))
	tmp	<- subset(tmp, SXO%in%sxo.select & SXO_N>0)
	
	ans	<- merge(ans, tmp, by=c('SUBTREE_ID','REP'),allow.cartesian=TRUE)	
	tmp	<- ans[, list(	SXO=SXO, 
						NL_TRMS_SXO=SUBTREE_SIZE*SXO_PROP,
						NL_TRMS_SXO_RND=SUBTREE_SIZE*SXO_PROP+as.vector(rmultinom(1, SUBTREE_UNOBS, SXO_PROP))), 
				by=c('SUBTREE_ID','REP','INTRO_SMPL','MISSING_SMPL')]
	ans	<- merge(ans, tmp, by=c('SUBTREE_ID','REP','INTRO_SMPL','MISSING_SMPL','SXO'))
	tmp	<- ans[, list(NL_TRMS_RND=sum(NL_TRMS_SXO_RND)), by=c('SUBTREE_ID','INTRO_SMPL','MISSING_SMPL')]
	ans	<- merge(ans, tmp, by=c('SUBTREE_ID','INTRO_SMPL','MISSING_SMPL'))	
	tmp	<- melt(subset(ans, MISSING_SMPL==1), id.vars=c('SXO','MISSING_SMPL','INTRO_SMPL','REP','SUBTREE_ID','SUBTREE_SIZE','NL_TRMS_RND'), measure.vars='NL_TRMS_SXO')
	ans	<- melt(ans, id.vars=c('SXO','MISSING_SMPL','INTRO_SMPL','REP','SUBTREE_ID','SUBTREE_SIZE','NL_TRMS_RND'), measure.vars='NL_TRMS_SXO_RND')
	ans	<- rbind(ans, tmp)
	
	ans[, NETWORK_SIZE_C:= NA_character_]
	tmp	<- ans[, which(variable=='NL_TRMS_SXO_RND')]
	set(ans, tmp, 'NETWORK_SIZE_C', ans[tmp, as.character(cut(NL_TRMS_RND, breaks=size.breaks, labels=size.labels, right=TRUE))])
	tmp	<- ans[, which(variable=='NL_TRMS_SXO')]
	set(ans, tmp, 'NETWORK_SIZE_C', ans[tmp, as.character(cut(SUBTREE_SIZE, breaks=size.breaks, labels=size.labels, right=TRUE))])
	
	ans	<- ans[,list(IND_IN_SUBTREE=sum(value)), by=c('SXO','MISSING_SMPL','INTRO_SMPL','REP','NETWORK_SIZE_C','variable')]	
	ans	<- ans[, list(NETWORK_SIZE_C=NETWORK_SIZE_C, N=IND_IN_SUBTREE, P=IND_IN_SUBTREE/sum(IND_IN_SUBTREE)), by=c('SXO','MISSING_SMPL','INTRO_SMPL','REP','variable')]
	
	ans	<- ans[, list(	P=ans.quantiles, 				
						QN=round(quantile(N, p=ans.quantiles)),
						QF=quantile(P, p=ans.quantiles)
						), by=c('SXO','variable','NETWORK_SIZE_C')]	
	setkey(ans, variable, SXO, NETWORK_SIZE_C, P)
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
	#	df<- copy(dpo); sampling.prob=0.43; resample.intro.n=1; resample.missing.n=100; sxo.select= c('MSM','hsx'); intro.year=2000; ans.quantiles=c(0.025,0.25,0.5,0.75,0.975); seed=42; bs.len=1e2
	#	size.breaks=size.breaks, size.labels=size.labels, 
	set.seed(seed)	
	#df		<- subset(ds, ASR_TYPE=='inf')	
	dsa		<- subset(df, select=c(SUBTREE_ID, SUBTREE_SIZE, ANCESTRAL_STATE, REP))
	tmp		<- dsa[, list(BS=seq_len(bs.len), SUBTREE_UNOBS=rnbinom(bs.len, SUBTREE_SIZE, prob=sampling.prob)), by='SUBTREE_ID']
	dsa		<- merge(dsa, tmp, by='SUBTREE_ID',allow.cartesian=TRUE)	
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

vli.asr.get.tips.in.NL.subtrees<- function(ph, da, dm, id.regex, select.descendants.regex="^.*_infNL_.*$", deadend.descendants.regex=NA)
{
	#	determine max parsimony state transitions into NL 
	setnames(da, colnames(da), gsub('\\.','_',toupper(colnames(da))))
	#	find all descendants of interest that are in NL subtrees
	#	these are in subtrees that are in NL state
	#	plus potentially other taxa of individuals that were sampled in NL, but not infected in NL, and don t show up in the ancestral state reconstrution
	if(is.na(deadend.descendants.regex))
	{
		subtrees.NL						<- subset(da, HOSTS=='NL')
	}
	if(!is.na(deadend.descendants.regex))
	{
		tmp								<- subset(da, HOSTS=='NL')[, ROOT_NOS]
		descendants.in.NL.subtrees		<- sapply(Descendants(ph, tmp, type='tips'),function(x) ph$tip.label[x][grepl(select.descendants.regex, ph$tip.label[x])])
		descendants.in.NL.subtrees		<- unlist(descendants.in.NL.subtrees)
		descendants.in.NL.subtrees		<- match(descendants.in.NL.subtrees, ph$tip.label)
		#	check if other subtrees contain any other descendants of interest
		#	if yes, create new splits from these
		tmp								<- subset(da, HOSTS!='NL')
		tmp								<- tmp[, 	{
					#ROOT_NOS<- 6441L
					z<- unlist(Descendants(ph, ROOT_NOS, type='tips'))
					z<- setdiff(z[grepl(deadend.descendants.regex, ph$tip.label[z])],descendants.in.NL.subtrees)
					if(!length(z))
						z<- NA_integer_
					list(DESC_NOTIN_NLSUBTREES=z)
				}, by='ROOT_NOS']									
		subtrees.NL						<- subset(tmp, !is.na(DESC_NOTIN_NLSUBTREES))
		#subtrees.NL[, TAXA:= ph$tip.label[DESC_NOTIN_NLSUBTREES]]		
		#merge(subtrees.NL, z, by='TAXA', all.x=TRUE)
		#
		#	get depth from root_no to tip
		subtrees.NL[, DEPTH:=NA_integer_]
		if(nrow(subtrees.NL))
			subtrees.NL					<- subtrees.NL[, {
							z<- c(DESC_NOTIN_NLSUBTREES, Ancestors(ph, DESC_NOTIN_NLSUBTREES, type="all"))
							z<- sapply(ROOT_NOS, function(x) which(z==x))
							list(ROOT_NOS=ROOT_NOS, DEPTH=z-1L)
						}, by='DESC_NOTIN_NLSUBTREES']
		#	select parent subtree with min depth, which can be self		
		subtrees.NL						<- subtrees.NL[, list(ROOT_NOS=ROOT_NOS[which.min(DEPTH)]), by='DESC_NOTIN_NLSUBTREES']
		#
		subtrees.NL[, SPLITS:= paste0('NL-DEADEND-SPLIT',seq_len(nrow(subtrees.NL)))]
		tmp								<- subset(merge(da, unique(subset(subtrees.NL, select=ROOT_NOS)), by='ROOT_NOS'), select=c(ROOT_NOS, UNIQUE_SPLITS, HOSTS))
		tmp[, PARENT_HOSTS:=HOSTS]
		subtrees.NL						<- merge(subtrees.NL, tmp, by='ROOT_NOS')
		setnames(subtrees.NL, c('ROOT_NOS','DESC_NOTIN_NLSUBTREES','SPLITS','UNIQUE_SPLITS'), c('PARENT_ROOT_NOS','ROOT_NOS','UNIQUE_SPLITS','PARENT_SPLITS'))
		subtrees.NL[, LENGTHS:=NA_real_]
		if(nrow(subtrees.NL))
			set(subtrees.NL, NULL, 'LENGTHS', ph$edge.length[sapply(subtrees.NL$ROOT_NOS, function(x) which(ph$edge[, 2]==x))])
		subtrees.NL[, PARENT_ROOT_NOS:=NULL]
		subtrees.NL						<- rbind(subset(da, HOSTS=='NL'),subtrees.NL)    		
	}
	#	add depth of subtree MRCAs
	tmp			<- subtrees.NL[, list(DEPTH=1L+length(Ancestors(ph, ROOT_NOS, type='all'))) , by='ROOT_NOS']
	subtrees.NL	<- merge(subtrees.NL, tmp, by='ROOT_NOS')
	#	get taxa
	ans	<- subtrees.NL[, list(IDX=unlist(Descendants(ph, ROOT_NOS, type='tips'))), by=c('UNIQUE_SPLITS','ROOT_NOS','DEPTH')]
	ans[, TAXA:= ph$tip.label[IDX]]
	ans[, ID:=gsub(id.regex,'\\1',TAXA)]
	ans	<- subset(ans, grepl(select.descendants.regex,TAXA) | (grepl('DEADEND',UNIQUE_SPLITS) & grepl(deadend.descendants.regex,TAXA)))	
	ans	<- merge(dm, ans, by='ID')	# we will need this below too
	#	make sure every taxon appears exactly once, keep the subtree with highest depth
	#	the fundamental problem is that subtree are not specified fully by MRCA
	#	so we could list for each subtree all inner nodes that constitute it, but that s not really needed for summarising the tips
	tmp	<- ans[, list(ROOT_NOS=ROOT_NOS[which.max(DEPTH)]), by='IDX']
	ans	<- merge(ans, tmp, by=c('ROOT_NOS','IDX'))
	set(ans, NULL, c('TAXA','DEPTH','ROOT_NOS','IDX'), NULL)
	setnames(ans, 'UNIQUE_SPLITS', 'SUBTREE_ID')
	ans	
}

vli.asr.summary<- function(ph, da, dm, id.regex, select.descendants.regex="^.*_infNL_.*$", deadend.descendants.regex=NA, is.dated.tree=TRUE)
{
	#	determine max parsimony state transitions into NL 
	setnames(da, colnames(da), gsub('\\.','_',toupper(colnames(da))))
	#	find all descendants of interest that are in NL subtrees
	#	these are in subtrees that are in NL state
	#	plus potentially other taxa of individuals that were sampled in NL, but not infected in NL, and don t show up in the ancestral state reconstrution
	if(is.na(deadend.descendants.regex))
	{
		subtrees.NL						<- subset(da, HOSTS=='NL')
	}
	if(!is.na(deadend.descendants.regex))
	{
		tmp								<- subset(da, HOSTS=='NL')[, ROOT_NOS]
		descendants.in.NL.subtrees		<- sapply(Descendants(ph, tmp, type='tips'),function(x) ph$tip.label[x][grepl(select.descendants.regex, ph$tip.label[x])])
		descendants.in.NL.subtrees		<- unlist(descendants.in.NL.subtrees)
		descendants.in.NL.subtrees		<- match(descendants.in.NL.subtrees, ph$tip.label)
		#	check if other subtrees contain any other descendants of interest
		#	if yes, create new splits from these
		tmp								<- subset(da, HOSTS!='NL')
		tmp								<- tmp[, 	{
					z<- unlist(Descendants(ph, ROOT_NOS, type='tips'))
					z<- setdiff(z[grepl(deadend.descendants.regex, ph$tip.label[z])],descendants.in.NL.subtrees)
					if(!length(z))
						z<- NA_integer_					
					list(DESC_NOTIN_NLSUBTREES=z)
				}, by='ROOT_NOS']									
		subtrees.NL						<- subset(tmp, !is.na(DESC_NOTIN_NLSUBTREES))
		#	get depth from root_no to tip
		subtrees.NL[, DEPTH:=NA_integer_]
		if(nrow(subtrees.NL))
			subtrees.NL					<- subtrees.NL[, {
						z<- c(DESC_NOTIN_NLSUBTREES, Ancestors(ph, DESC_NOTIN_NLSUBTREES, type="all"))
						z<- sapply(ROOT_NOS, function(x) which(z==x))
						list(ROOT_NOS=ROOT_NOS, DEPTH=z-1L)
					}, by='DESC_NOTIN_NLSUBTREES']
		#	select parent subtree with min depth, which can be self		
		subtrees.NL						<- subtrees.NL[, list(ROOT_NOS=ROOT_NOS[which.min(DEPTH)]), by='DESC_NOTIN_NLSUBTREES']
		#
		subtrees.NL[, SPLITS:= paste0('NL-DEADEND-SPLIT',seq_len(nrow(subtrees.NL)))]
		tmp								<- subset(merge(da, unique(subset(subtrees.NL, select=ROOT_NOS)), by='ROOT_NOS'), select=c(ROOT_NOS, UNIQUE_SPLITS, HOSTS))		
		tmp[, PARENT_HOSTS:=HOSTS]
		subtrees.NL						<- merge(subtrees.NL, tmp, by='ROOT_NOS')
		setnames(subtrees.NL, c('ROOT_NOS','DESC_NOTIN_NLSUBTREES','SPLITS','UNIQUE_SPLITS'), c('PARENT_ROOT_NOS','ROOT_NOS','UNIQUE_SPLITS','PARENT_SPLITS'))
		subtrees.NL[, LENGTHS:=NA_real_]
		if(nrow(subtrees.NL))
			set(subtrees.NL, NULL, 'LENGTHS', ph$edge.length[sapply(subtrees.NL$ROOT_NOS, function(x) which(ph$edge[, 2]==x))])
		subtrees.NL[, PARENT_ROOT_NOS:=NULL]
		subtrees.NL						<- rbind(subset(da, HOSTS=='NL'),subtrees.NL)    		
	}
	
	
	#	determine time of nodes at either end of branch with state transition	
	tmp			<- subtrees.NL[, list(EDGE_ID= which(ph$edge[,2]==ROOT_NOS) ), by='UNIQUE_SPLITS']
	subtrees.NL	<- merge(subtrees.NL, tmp, by='UNIQUE_SPLITS')
	subtrees.NL[, BRL:= ph$edge.length[ EDGE_ID] ]	#this should be same as LENGTHS
	tmp			<- node.depth.edgelength(ph)
	subtrees.NL[, ROOT_HEIGHT:= max(tmp) - tmp[ ROOT_NOS ]]	# this is the height between the inner node and the most recent tip
	if(is.dated.tree)
	{
		tmp2	<- gsub(id.regex,'\\1',ph$tip.label[which.max(tmp)])	# this is the ID of the most recent tip
		tmp2	<- subset(dm, ID==tmp2)[, SAMPLING_DATE]					# this is the corresponding sampling date
		subtrees.NL[, ROOT_DATE:= tmp2-ROOT_HEIGHT]
		subtrees.NL[, PARENT_ROOT_DATE:= ROOT_DATE-LENGTHS]
	}
	if(!is.dated.tree)
	{
		set(subtrees.NL, NULL, c('ROOT_DATE','PARENT_ROOT_DATE'), NA_real_)
	}
	
	#	determine TAXA that are in this particular tree, and not any other tree
	#	add depth of subtree MRCAs
	tmp			<- subtrees.NL[, list(DEPTH=1L+length(Ancestors(ph, ROOT_NOS, type='all'))) , by='ROOT_NOS']
	subtrees.NL	<- merge(subtrees.NL, tmp, by='ROOT_NOS')
	#	get taxa
	ans	<- subtrees.NL[, list(IDX=unlist(Descendants(ph, ROOT_NOS, type='tips'))), by=c('UNIQUE_SPLITS','ROOT_NOS','DEPTH')]
	ans[, TAXA:= ph$tip.label[IDX]]
	ans[, ID:=gsub(id.regex,'\\1',TAXA)]
	ans	<- subset(ans, grepl(select.descendants.regex,TAXA) | (grepl('DEADEND',UNIQUE_SPLITS) & grepl(deadend.descendants.regex,TAXA)))	
	tmp	<- ans[, list(ROOT_NOS=ROOT_NOS[which.max(DEPTH)]), by='IDX']
	ans	<- merge(ans, tmp, by=c('ROOT_NOS','IDX'))
	
	
	#	determine number of tips per subtree
	ans			<- merge(ans, ans[, list(SUBTREE_SIZE=length(IDX)), by='ROOT_NOS'], by='ROOT_NOS')
	#	determine composition of subtrees by sexual orientation
	tmp			<- merge(dm, subset(ans, select=c(ID,UNIQUE_SPLITS)), by='ID')	# we will need this below too
	setnames(tmp, 'SXO', 'L')
	tmp2		<- as.data.table(expand.grid(UNIQUE_SPLITS=unique(tmp$UNIQUE_SPLITS), L=unique(subset(dm, SXO!='LANL', SXO))[, SXO] ))	
	tmp2		<- merge(tmp2, tmp[, list(N=length(ID)), by=c('UNIQUE_SPLITS','L')], by=c('UNIQUE_SPLITS','L'), all.x=TRUE)
	set(tmp2, tmp2[, which(is.na(N))],'N',0L)	
	set(tmp2, NULL, 'L', tmp2[,paste0('SUBTREE_SXO_',L)])
	tmp2		<- dcast.data.table(tmp2, UNIQUE_SPLITS~L, value.var='N')
	ans	<- merge(ans, tmp2, by='UNIQUE_SPLITS')
	setnames(tmp, 'L', 'SXO')
	
	
	#	determine composition of subtrees by birth location
	setnames(tmp, 'BORN_LOC', 'L')
	tmp2		<- as.data.table(expand.grid(UNIQUE_SPLITS=unique(tmp$UNIQUE_SPLITS), L=unique(subset(dm, select=BORN_LOC))[, BORN_LOC] ))	
	tmp2		<- merge(tmp2, tmp[, list(N=length(ID)), by=c('UNIQUE_SPLITS','L')], by=c('UNIQUE_SPLITS','L'), all.x=TRUE)
	set(tmp2, tmp2[, which(is.na(N))],'N',0L)	
	set(tmp2, NULL, 'L', tmp2[,paste0('BORN_LOC_',tolower(L))])
	tmp2		<- dcast.data.table(tmp2, UNIQUE_SPLITS~L, value.var='N')
	ans			<- merge(ans, tmp2, by='UNIQUE_SPLITS')
	setnames(tmp, 'L', 'BORN_LOC')
	
	#	remove individuals and keep summary
	set(ans, NULL, c('ROOT_NOS','IDX','DEPTH','TAXA','ID'), NULL)
	ans			<- unique(ans)
	subtrees.NL	<- merge(subtrees.NL, ans, by='UNIQUE_SPLITS')
	
	#	clean up 
	set(subtrees.NL, NULL, c('PARENT_SPLITS','HOSTS','BRL','LENGTHS','ROOT_HEIGHT','ROOT_NOS','EDGE_ID','DEPTH'), NULL)
	setnames(subtrees.NL, c('UNIQUE_SPLITS','PARENT_HOSTS','ROOT_DATE','PARENT_ROOT_DATE'), c('SUBTREE_ID','ANCESTRAL_STATE','INTRO_LATEST','INTRO_EARLIEST')	)
	subtrees.NL[, DEAD_END:= as.integer(grepl('DEADEND',SUBTREE_ID))]
	subtrees.NL
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

vli.bez.explore.ancestral.state.reconstruction.170901.plot.distinct.networks<- function()
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

vli.bez.explore.ancestral.state.reconstruction.171113.origin.of.subtrees<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	
	asr.type		<- 'inf'; asr.label	<- 'subtrees by location of infection'
	#asr.type	<- 'sam'; asr.label	<- 'subtrees by sampling location'	
	#asr.type	<- 'bir'; asr.label	<- 'subtrees by location of birth'
	
	outfile.base	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_'
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_summaries.rda'
	#
	#	get subtrees subst/site, inf ancestral state reconstruction, non-bootstrap 
	load(infile)	
	dss				<- subset(ds, TREE_TYPE!='lsd' & ASR_TYPE==asr.type)
	set(dss, NULL, 'SUBTREE_ID', dss[,paste0(ST,'-',REP,'-',SUBTREE_ID)])	
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_tips.rda'
	load(infile)	
	dti				<- subset(dti, TREE_TYPE!='lsd' & ASR_TYPE==asr.type)
	set(dti, NULL, 'ST', dti[,gsub('pol','',ST)])
	set(dti, NULL, 'SUBTREE_ID', dti[,paste0(ST,'-',REP,'-',SUBTREE_ID)])
	dp				<- unique(subset(dss, select=c(SUBTREE_ID, ASR_TYPE, ST, REP, ANCESTRAL_STATE, SUBTREE_SIZE, SUBTREE_SXO_MSM, SUBTREE_SXO_U, SUBTREE_SXO_hsx, SUBTREE_SXO_other)))	
	tmp				<- dti[, list(MIN_SAMPLING_DATE=min(SAMPLING_DATE), MAX_SAMPLING_DATE=max(SAMPLING_DATE)), by=c('SUBTREE_ID')]
	dp				<- merge(dp, tmp, by='SUBTREE_ID')
	dp				<- subset(dp, MIN_SAMPLING_DATE>1995)
	dp[, MIN_SAMPLING_DATE_C:= cut(MIN_SAMPLING_DATE, breaks=c(1930,1950,1970,1990,seq(1992,2016,2)))]
	dp[, STB:= as.character(factor(ST=='B', levels=c(TRUE,FALSE), labels=c('B','non-B')))]
	
	dc	<- dp[, list(N_TREE=length(SUBTREE_ID), N_PATIENT=sum(SUBTREE_SIZE)), by=c('STB','REP','ANCESTRAL_STATE')]
	dc	<- merge(as.data.table(expand.grid(STB=unique(dc$STB), REP=unique(dc$REP), ANCESTRAL_STATE=unique(dc$ANCESTRAL_STATE))), dc, by=c('STB','REP','ANCESTRAL_STATE'), all.x=TRUE)
	set(dc, dc[, which(is.na(N_TREE))], 'N_TREE', 0L)
	set(dc, dc[, which(is.na(N_PATIENT))], 'N_PATIENT', 0L)	
	dc	<- merge(dc, dc[, list(ANCESTRAL_STATE=ANCESTRAL_STATE, P_TREE=N_TREE/sum(N_TREE), P_PATIENT=N_PATIENT/sum(N_PATIENT)), by=c('STB','REP')], by=c('STB','REP','ANCESTRAL_STATE'))
	dc	<- melt(dc, id.vars=c('STB','REP','ANCESTRAL_STATE'))
	dc	<- dc[, list(	P=paste0('p',c(0.025,0.25,0.5,0.75,0.975)),
						Q=quantile(value, prob=c(0.025,0.25,0.5,0.75,0.975))
						), by=c('STB','variable','ANCESTRAL_STATE')]
	dc	<- dcast.data.table(dc, variable+STB+ANCESTRAL_STATE~P, value.var='Q')	
	dc[, LABEL:= NA_character_]
	tmp	<- dc[, which(substr(variable,1,1)=='P')]
	set(dc, tmp, 'LABEL', dc[tmp, paste0(round(100*p0.5, d=1),'% [',round(100*p0.025, d=1),'%-',round(100*p0.975, d=1),'%]')])
	tmp	<- dc[, which(substr(variable,1,1)=='N')]
	set(dc, tmp, 'LABEL', dc[tmp, paste0(round(p0.5, d=1),' [',round(p0.025, d=1),'-',round(p0.975, d=1),']')])
	set(dc, NULL, 'ANCESTRAL_STATE', dc[, factor(ANCESTRAL_STATE,	levels=c("SSA","EUC","EUO","EUW","LAT","OAP","ZAZ","NAM","USCAU","unsampled_region"),
																	labels=c("Sub-Saharan Africa","East EU","Former USSR","West EU","Latin America","China, Japan, Korea","South and Southeast\nAsia","North Africa, Middle East","US, Canada","unresolved"))
																	])
	save(dc, file=paste0(outfile.base, asr.type, '_ALLST_origin.rda'))
	write.csv(dc, row.names=FALSE, file=paste0(outfile.base, asr.type, '_ALLST_origin.csv'))
	ggplot(subset(dc, variable=='P_TREE'), aes(x=ANCESTRAL_STATE, y=p0.5, fill=ANCESTRAL_STATE, ymin=p0.025, ymax=p0.975)) +
			geom_bar(stat='identity', position=position_dodge(width=0.8)) +
			geom_errorbar(position=position_dodge(width=0.8), width=0.5, colour="black")	+
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,0.66)) +
			facet_grid(~STB) + coord_flip() +
			labs(x='geographic origin\n',y='\nphylogenetically identified transmission chains', fill='geographic origin')
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_origin.pdf'), w=10, h=6)
	
	tmp	<- subset(dc, variable=='N_TREE')[, list(p0.025=sum(p0.025), p0.975=sum(p0.975)), by='ANCESTRAL_STATE']
	ggplot(subset(dc, variable=='N_TREE'), aes(x=ANCESTRAL_STATE)) +
			geom_bar(aes(y=p0.5, fill=STB), stat='identity') +
			geom_errorbar(data=tmp, aes(ymin=p0.025, ymax=p0.975), width=0.5, colour="black")	+
			coord_flip() +			
			theme_bw() +
			labs(x='geographic origin\n',y='\nphylogenetically identified transmission chains', fill='subtype')
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_origin2.pdf'), w=10, h=8)
}

vli.bez.explore.ancestral.state.reconstruction.171113.debug<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	asr.type		<- 'inf'; asr.label	<- 'subtrees by location of infection'
	#asr.type	<- 'sam'; asr.label	<- 'subtrees by sampling location'	
	#asr.type	<- 'bir'; asr.label	<- 'subtrees by location of birth'
	
	outfile.base	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_'
	
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_summaries.rda'
	load(infile)	
	subset(ds, ST=='AG' & SUBTREE_ID=='NL-SPLIT387' & REP==0)
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_tips.rda'
	load(infile)	
	subset(dti, ST=='AG' & SUBTREE_ID=='NL-SPLIT387' & REP==0)
	subset(do, ST=='AG' & SUBTREE_ID=='NL-SPLIT387' & REP==0)
	
	infile.meta	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo_plusB_OR.csv'
	dm	<- as.data.table(read.csv(infile.meta))
	set(dm, NULL, 'SXO', dm[, gsub('HTF|HTM','hsx',gsub('BL|CH|CHS|DU','other',as.character(GENDER_TRMGROUP)))])	
	indir.asr	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results'
	infiles		<- data.table(F=list.files(indir.asr, pattern='^workspace', full.names=TRUE, recursive=TRUE))
	infiles[, ASR_TYPE:= gsub('.*_annotated_([a-z]+).*','\\1',basename(F))]
	infiles[, REP:= as.integer(gsub('^.*_([0-9]+)_annotated.*$|^.*_rr([0-9]+)_annotated.*$','\\1',basename(F)))]
	tmp			<- infiles[, which(is.na(REP))]
	set(infiles, tmp, 'REP', infiles[tmp,as.integer(gsub('^.*_rr([0-9]+)_annotated.*$','\\1',basename(F)))])	
	infiles[, ST:= gsub('pol','',gsub('^workspace_([^_]+)_.*','\\1',basename(F)))]
	infiles[, TREE_TYPE:=gsub('^.*(lsd|rr).*$','\\1',basename(F))]
	set(infiles, infiles[, which(TREE_TYPE!='lsd')],'TREE_TYPE','substsite')
	infiles[, SELECT_DESC:= as.character(factor(ASR_TYPE, levels=c('bir','sam','inf'), labels=c('^.*_bornNL_.*$','^.*_sampNL_.*$','^.*_infNL_.*$')))]
	infiles[, DEADEND_DESC:= NA_character_]
	infiles[, IS_DATED_TREE:= TREE_TYPE=='lsd']
	set(infiles, infiles[, which(ASR_TYPE=='inf')], 'DEADEND_DESC',  "^.*_sampNL_.*$|^.*_infNL_.*$")		
	infiles		<- subset(infiles, !IS_DATED_TREE)
	id.regex	<- '^([^_]*).*'
	subset(infiles, REP==0 & ST=='AG')
	
	do<- copy(dti)
	
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results/AG/sam/workspace_AG_withD_rr000_annotated_sam.rda')
	#	F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results/AG/sam/workspace_AG_withD_rr000_annotated_sam.rda'
	select.descendants.regex	<- '^.*_sampNL_.*$'
	is.dated.tree				<- 0	
	deadend.descendants.regex	<- NA_character_						
	cat('\nprocessing file ',F,'\n',select.descendants.regex,' ',deadend.descendants.regex,' ',is.dated.tree)
	load(F)								
	#	read variables from workspace
	ph	<- all.tree.info$only.tree$tree
	da	<- as.data.table(all.tree.info$only.tree$classification.results$collapsed)	
	#	summarise introductions into NL
	dti	<- vli.asr.get.tips.in.NL.subtrees(ph, da, dm, id.regex, select.descendants.regex=select.descendants.regex, deadend.descendants.regex=deadend.descendants.regex)				
	subset(dti, SUBTREE_ID=='NL-SPLIT387')
	
	
	ph$edge[ph$edge[,2]==which(grepl('M35476',ph$tip.label)),]
	#8885 4409
	ph$edge[ph$edge[,2]==8885,]
	#8106 8107
	attr(ph,"SUBGRAPH_MRCA")[8885]
	attr(ph,"SUBGRAPH_MRCA")[8884]
	attr(ph,"BRANCH_COLOURS")[8885]
	attr(ph,"BRANCH_COLOURS")[8884]
	
	
	
	ph$edge[ph$edge[,2]==which(grepl('M10041',ph$tip.label)),]
	#8877 4394
	ph$edge[ph$edge[,2]==8877,]
	#8876 8877
	subset(da, ROOT_NOS==8877)
	ph$tip.label[ unlist(Descendants(ph, 8877, type='tips')) ]	#correct: only 2 taxa
	
	any(Descendants(ph, 4481, type='all')==8877)	#8877 is a descendant
	attr(ph,"BRANCH_COLOURS")[8876]
	#	so within the subtree called 4481, there are subtrees of other colours
	
}

vli.bez.explore.ancestral.state.reconstruction.171113.size.networks<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	asr.type		<- 'inf'; asr.label	<- 'subtrees by location of infection'
	#asr.type	<- 'sam'; asr.label	<- 'subtrees by sampling location'	
	#asr.type	<- 'bir'; asr.label	<- 'subtrees by location of birth'
	
	outfile.base	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_'
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_summaries.rda'
	#
	#	get subtrees subst/site, inf ancestral state reconstruction, non-bootstrap 
	load(infile)	
	dss				<- subset(ds, TREE_TYPE!='lsd' & ASR_TYPE==asr.type)
	set(dss, NULL, 'SUBTREE_ID', dss[,paste0(ST,'-',REP,'-',SUBTREE_ID)])	
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_tips.rda'
	load(infile)	
	dti				<- subset(dti, TREE_TYPE!='lsd' & ASR_TYPE==asr.type)
	set(dti, NULL, 'ST', dti[,gsub('pol','',ST)])
	set(dti, NULL, 'SUBTREE_ID', dti[,paste0(ST,'-',REP,'-',SUBTREE_ID)])		
	#subset(dti, REP==0 & ST!='B')[, table(grepl('DEADEND',SUBTREE_ID))]	
	dp				<- unique(subset(dss, select=c(SUBTREE_ID, ASR_TYPE, ST, REP, ANCESTRAL_STATE, SUBTREE_SIZE, SUBTREE_SXO_MSM, SUBTREE_SXO_U, SUBTREE_SXO_hsx, SUBTREE_SXO_other)))	
	tmp				<- dti[, list(MIN_SAMPLING_DATE=min(SAMPLING_DATE), MAX_SAMPLING_DATE=max(SAMPLING_DATE)), by=c('SUBTREE_ID')]
	dp				<- merge(dp, tmp, by='SUBTREE_ID')
	dp				<- subset(dp, MIN_SAMPLING_DATE>1995)
	dp[, MIN_SAMPLING_DATE_C:= cut(MIN_SAMPLING_DATE, breaks=c(1930,1950,1970,1990,seq(1992,2016,2)))]	
	dpb			<- subset(dp, ST=='B')
	dpo			<- subset(dp, ST!='B')
	
	#
	#	subtrees of sizes <3, of subtrees introduced since 1995
	#
	size.breaks	<- c(0,2,50000); size.labels<- c('1-2','>2')
	df			<- vli.estimate.subtree.distribution.by.sexual.orientation(dpo, sampling.prob=0.43, resample.intro.n=1, resample.missing.n=100, sxo.select= c('MSM','hsx'), intro.year=1995, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), size.breaks=size.breaks, size.labels=size.labels, seed=42)	
	set(df, NULL, 'P', df[, paste0('p',P)])
	df[, STB:='non-B']
	ans			<- copy(df)	
	df			<- vli.estimate.subtree.distribution.by.sexual.orientation(dpb, sampling.prob=0.43, resample.intro.n=1, resample.missing.n=100, sxo.select= c('MSM','hsx'), intro.year=1995, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), size.breaks=size.breaks, size.labels=size.labels, seed=42)	
	set(df, NULL, 'P', df[, paste0('p',P)])
	df[, STB:='B']
	ans			<- rbind(ans, df)
	setnames(ans, 'VALUE_C', 'NETWORK_SIZE_C')
	setkey(ans, variable, STB, SXO, NETWORK_SIZE_C, P)
	save(ans, file=paste0(outfile.base, asr.type, '_ALLST_subtreesize_distribution_breaks12.rda'))	
	tmp			<- dcast.data.table(ans, STB+variable+SXO+NETWORK_SIZE_C~P, value.var='QF')
	tmp[, LABEL:= paste0(round(100*p0.5, d=1),'% [',round(100*p0.025, d=1),'%-',round(100*p0.975, d=1),'%]')]
	tmp			<- subset(tmp, select=c(STB, variable, SXO, NETWORK_SIZE_C, LABEL))
	tmp			<- dcast.data.table(tmp, variable+STB+SXO~NETWORK_SIZE_C, value.var='LABEL')	
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('NL_TRMS_SXO','NL_TRMS_SXO_RND'), labels=c(paste0('subtrees\nas observed in phylogeny'),paste0('subtrees\nadjusted for missing sequences')))])
	write.csv(tmp, row.names=FALSE, file=paste0(outfile.base, asr.type, '_ALLST_subtreesize_distribution_breaks12.csv'))
	
	size.breaks	<- c(0,3,50000); size.labels<- c('1-3','>3')
	df			<- vli.estimate.subtree.distribution.by.sexual.orientation(dpo, sampling.prob=0.43, resample.intro.n=1, resample.missing.n=100, sxo.select= c('MSM','hsx'), intro.year=1995, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), size.breaks=size.breaks, size.labels=size.labels, seed=42)	
	set(df, NULL, 'P', df[, paste0('p',P)])
	df[, STB:='non-B']
	ans			<- copy(df)	
	df			<- vli.estimate.subtree.distribution.by.sexual.orientation(dpb, sampling.prob=0.43, resample.intro.n=1, resample.missing.n=100, sxo.select= c('MSM','hsx'), intro.year=1995, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), size.breaks=size.breaks, size.labels=size.labels, seed=42)	
	set(df, NULL, 'P', df[, paste0('p',P)])
	df[, STB:='B']
	ans			<- rbind(ans, df)
	setnames(ans, 'VALUE_C', 'NETWORK_SIZE_C')
	setkey(ans, variable, STB, SXO, NETWORK_SIZE_C, P)
	save(ans, file=paste0(outfile.base, asr.type, '_ALLST_subtreesize_distribution_breaks12.rda'))	
	tmp			<- dcast.data.table(ans, STB+variable+SXO+NETWORK_SIZE_C~P, value.var='QF')
	tmp[, LABEL:= paste0(round(100*p0.5, d=1),'% [',round(100*p0.025, d=1),'%-',round(100*p0.975, d=1),'%]')]
	tmp			<- subset(tmp, select=c(STB, variable, SXO, NETWORK_SIZE_C, LABEL))
	tmp			<- dcast.data.table(tmp, variable+STB+SXO~NETWORK_SIZE_C, value.var='LABEL')	
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('NL_TRMS_SXO','NL_TRMS_SXO_RND'), labels=c(paste0('subtrees\nas observed in phylogeny'),paste0('subtrees\nadjusted for missing sequences')))])
	write.csv(tmp, row.names=FALSE, file=paste0(outfile.base, asr.type, '_ALLST_subtreesize_distribution_breaks13.csv'))
}

vli.bez.explore.ancestral.state.reconstruction.171113.patients.in.small.networks<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	asr.type		<- 'inf'; asr.label	<- 'subtrees by location of infection'
	#asr.type	<- 'sam'; asr.label	<- 'subtrees by sampling location'	
	#asr.type	<- 'bir'; asr.label	<- 'subtrees by location of birth'
	
	outfile.base	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_'
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_summaries.rda'
	#
	#	get subtrees subst/site, inf ancestral state reconstruction, non-bootstrap 
	load(infile)	
	dss				<- subset(ds, TREE_TYPE!='lsd' & ASR_TYPE==asr.type)
	set(dss, NULL, 'SUBTREE_ID', dss[,paste0(ST,'-',REP,'-',SUBTREE_ID)])	
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_tips.rda'
	load(infile)	
	dti				<- subset(dti, TREE_TYPE!='lsd' & ASR_TYPE==asr.type)
	set(dti, NULL, 'ST', dti[,gsub('pol','',ST)])
	set(dti, NULL, 'SUBTREE_ID', dti[,paste0(ST,'-',REP,'-',SUBTREE_ID)])	
	dp				<- unique(subset(dss, select=c(SUBTREE_ID, ASR_TYPE, ST, REP, ANCESTRAL_STATE, SUBTREE_SIZE, SUBTREE_SXO_MSM, SUBTREE_SXO_U, SUBTREE_SXO_hsx, SUBTREE_SXO_other)))	
	tmp				<- dti[, list(MIN_SAMPLING_DATE=min(SAMPLING_DATE), MAX_SAMPLING_DATE=max(SAMPLING_DATE)), by=c('SUBTREE_ID')]
	dp				<- merge(dp, tmp, by='SUBTREE_ID')
	dp				<- subset(dp, MIN_SAMPLING_DATE>1995)
	dp[, MIN_SAMPLING_DATE_C:= cut(MIN_SAMPLING_DATE, breaks=c(1930,1950,1970,1990,seq(1992,2016,2)))]	
	dpb			<- subset(dp, ST=='B')
	dpo			<- subset(dp, ST!='B')
	
	#
	#	subtrees of sizes <3, of subtrees introduced since 1995
	#
	size.breaks	<- c(0,2,5e4); size.labels<- c('1-2','>2')	
	df			<- vli.estimate.patients.in.small.subtrees.by.sexual.orientation(dpo, sampling.prob=0.43, resample.intro.n=1, resample.missing.n=100, sxo.select= c('MSM','hsx'), intro.year=1995, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), size.breaks=size.breaks, size.labels=size.labels, seed=42)	
	set(df, NULL, 'P', df[, paste0('p',P)])
	df[, STB:='non-B']
	ans			<- copy(df)
	df			<- vli.estimate.patients.in.small.subtrees.by.sexual.orientation(dpb, sampling.prob=0.43, resample.intro.n=1, resample.missing.n=100, sxo.select= c('MSM','hsx'), intro.year=1995, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), size.breaks=size.breaks, size.labels=size.labels, seed=42)	
	set(df, NULL, 'P', df[, paste0('p',P)])
	df[, STB:='B']
	ans			<- rbind(ans, df)
	setkey(ans, variable, STB, SXO, NETWORK_SIZE_C, P)
	save(ans, file=paste0(outfile.base, asr.type, '_ALLST_patientsinsmallnetworks_breaks12.rda'))	
	tmp			<- dcast.data.table(ans, STB+variable+SXO+NETWORK_SIZE_C~P, value.var='QF')
	set(tmp, NULL, 'NETWORK_SIZE_C', tmp[, paste0('subtree_size_',NETWORK_SIZE_C)])
	tmp[, LABEL:= paste0(round(100*p0.5, d=1),'% [',round(100*p0.025, d=1),'%-',round(100*p0.975, d=1),'%]')]
	tmp			<- subset(tmp, select=c(STB, variable, SXO, NETWORK_SIZE_C, LABEL))
	tmp			<- dcast.data.table(tmp, variable+STB+SXO~NETWORK_SIZE_C, value.var='LABEL')	
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('NL_TRMS_SXO','NL_TRMS_SXO_RND'), labels=c(paste0('patients\nas observed in phylogeny'),paste0('patients\nadjusted for missing sequences')))])
	write.csv(tmp, row.names=FALSE, file=paste0(outfile.base, asr.type, '_ALLST_patientsinsmallnetworks_breaks12.csv'))
	tmp			<- dcast.data.table(ans, STB+variable+SXO+NETWORK_SIZE_C~P, value.var='QN')
	set(tmp, NULL, 'NETWORK_SIZE_C', tmp[, paste0('subtree_size_',NETWORK_SIZE_C)])
	tmp[, LABEL:= paste0(round(p0.5, d=0),' [',round(p0.025, d=0),'-',round(p0.975, d=0),']')]
	tmp			<- subset(tmp, select=c(STB, variable, SXO, NETWORK_SIZE_C, LABEL))
	tmp			<- dcast.data.table(tmp, variable+STB+SXO~NETWORK_SIZE_C, value.var='LABEL')	
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('NL_TRMS_SXO','NL_TRMS_SXO_RND'), labels=c(paste0('patients\nas observed in phylogeny'),paste0('patients\nadjusted for missing sequences')))])
	write.csv(tmp, row.names=FALSE, file=paste0(outfile.base, asr.type, '_ALLST_patientsinsmallnetworks_breaks12_N.csv'))
	
	#
	#	subtrees of sizes <=3, of subtrees introduced since 1995
	#
	size.breaks	<- c(0,3,5e4); size.labels<- c('1-3','>3')	
	df			<- vli.estimate.patients.in.small.subtrees.by.sexual.orientation(dpo, sampling.prob=0.43, resample.intro.n=1, resample.missing.n=100, sxo.select= c('MSM','hsx'), intro.year=1995, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), size.breaks=size.breaks, size.labels=size.labels, seed=42)	
	set(df, NULL, 'P', df[, paste0('p',P)])
	df[, STB:='non-B']
	ans			<- copy(df)
	df			<- vli.estimate.patients.in.small.subtrees.by.sexual.orientation(dpb, sampling.prob=0.43, resample.intro.n=1, resample.missing.n=100, sxo.select= c('MSM','hsx'), intro.year=1995, ans.quantiles=c(0.025,0.25,0.5,0.75,0.975), size.breaks=size.breaks, size.labels=size.labels, seed=42)	
	set(df, NULL, 'P', df[, paste0('p',P)])
	df[, STB:='B']
	ans			<- rbind(ans, df)
	setkey(ans, variable, STB, SXO, NETWORK_SIZE_C, P)
	save(ans, file=paste0(outfile.base, asr.type, '_ALLST_patientsinsmallnetworks_breaks13.rda'))	
	tmp			<- dcast.data.table(ans, STB+variable+SXO+NETWORK_SIZE_C~P, value.var='QF')
	set(tmp, NULL, 'NETWORK_SIZE_C', tmp[, paste0('subtree_size_',NETWORK_SIZE_C)])
	tmp[, LABEL:= paste0(round(100*p0.5, d=1),'% [',round(100*p0.025, d=1),'%-',round(100*p0.975, d=1),'%]')]
	tmp			<- subset(tmp, select=c(STB, variable, SXO, NETWORK_SIZE_C, LABEL))
	tmp			<- dcast.data.table(tmp, variable+STB+SXO~NETWORK_SIZE_C, value.var='LABEL')	
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('NL_TRMS_SXO','NL_TRMS_SXO_RND'), labels=c(paste0('patients\nas observed in phylogeny'),paste0('patients\nadjusted for missing sequences')))])
	write.csv(tmp, row.names=FALSE, file=paste0(outfile.base, asr.type, '_ALLST_patientsinsmallnetworks_breaks13.csv'))
	tmp			<- dcast.data.table(ans, STB+variable+SXO+NETWORK_SIZE_C~P, value.var='QN')
	set(tmp, NULL, 'NETWORK_SIZE_C', tmp[, paste0('subtree_size_',NETWORK_SIZE_C)])
	tmp[, LABEL:= paste0(round(p0.5, d=0),' [',round(p0.025, d=0),'-',round(p0.975, d=0),']')]
	tmp			<- subset(tmp, select=c(STB, variable, SXO, NETWORK_SIZE_C, LABEL))
	tmp			<- dcast.data.table(tmp, variable+STB+SXO~NETWORK_SIZE_C, value.var='LABEL')	
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('NL_TRMS_SXO','NL_TRMS_SXO_RND'), labels=c(paste0('patients\nas observed in phylogeny'),paste0('patients\nadjusted for missing sequences')))])
	write.csv(tmp, row.names=FALSE, file=paste0(outfile.base, asr.type, '_ALLST_patientsinsmallnetworks_breaks13_N.csv'))
}

vli.bez.explore.ancestral.state.reconstruction.171113.check.trees<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	asr.type		<- 'inf'
	#asr.type		<- 'sam'
	outfile.base	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_'
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_summaries.rda'
	#
	#	get subtrees subst/site, inf ancestral state reconstruction, non-bootstrap 
	load(infile)	
	dss				<- subset(ds, TREE_TYPE!='lsd' & ASR_TYPE==asr.type)
	set(dss, NULL, 'SUBTREE_ID', dss[,paste0(ST,'-',REP,'-',SUBTREE_ID)])	
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_tips.rda'
	load(infile)	
	dti				<- subset(dti, TREE_TYPE!='lsd' & ASR_TYPE==asr.type)
	set(dti, NULL, 'ST', dti[,gsub('pol','',ST)])
	set(dti, NULL, 'SUBTREE_ID', dti[,paste0(ST,'-',REP,'-',SUBTREE_ID)])
	tmp				<- unique(subset(dss, select=c(SUBTREE_ID, ANCESTRAL_STATE, SUBTREE_SIZE)))
	dp				<- merge(tmp, dti, by='SUBTREE_ID')
	dp				<- merge(dp, dp[, list(MIN_SAMPLING_DATE=min(SAMPLING_DATE), MAX_SAMPLING_DATE=max(SAMPLING_DATE)), by='SUBTREE_ID'], by='SUBTREE_ID')
	dp				<- subset(dp, MIN_SAMPLING_DATE>1995)
	dp[, MIN_SAMPLING_DATE_C:= cut(MIN_SAMPLING_DATE, breaks=c(1930,1950,1970,1990,seq(1992,2016,2)))]	
	dp[, STB:= as.character(factor(ST=='B', levels=c(TRUE,FALSE),labels=c('B','non-B')))]		
	#
	
	subset(dp, ST=='AG')[, length(unique(REP))]
	subset(dp, ST=='C')[, length(unique(REP))]
	subset(dp, ST=='cpx06')[, length(unique(REP))]
	subset(dp, ST=='D')[, length(unique(REP))]
	subset(dp, ST=='B')[, length(unique(REP))]
	subset(dp, ST=='F')[, length(unique(REP))]
	subset(dp, ST=='G')[, length(unique(REP))]
	
	unique(dp, by=c('ST','REP'))[, table(ST)]
	subset(dp, ST=='A1')[, unique(REP)]	
	
	dc	<- dp[, list(N=length(unique(SUBTREE_ID))), by=c('ST','REP')]			
	ggplot(dc, aes(x=N)) + geom_histogram() + facet_wrap(~ST, scales='free')
	
	
	dc	<- dp[, list(N=length(unique(SUBTREE_ID))), by=c('STB','REP')]			
	ggplot(subset(dc, STB=='non-B'), aes(x=N)) + geom_histogram() 
	
	#subset(dp, REP==0 & ST!='B')[, table(grepl('DEADEND',SUBTREE_ID))]
	tmp	<- dp[, list(NP=length(ID)), by=c('ST','REP')]
	subset(tmp, ST=='A1')[, table(NP)]	#	ok
	subset(tmp, ST=='AG')[, table(NP)]	#	many odd
	subset(tmp, ST=='C')[, table(NP)]	# ok
	subset(tmp, ST=='cpx06')[, table(NP)]	#ok
	subset(tmp, ST=='D')[, table(NP)]	#	1 odd, REP is 53
	subset(tmp, ST=='B')[, table(NP)]	#	each run different!!!
	subset(tmp, ST=='F')[, table(NP)]	#	many odd
	subset(tmp, ST=='G')[, table(NP)]	#	ok
	
	
	ggplot( tmp, aes(x=NP)) + facet_wrap(~ST, scales='free') + geom_histogram(binwidth=10)
}

vli.bez.explore.ancestral.state.reconstruction.171113.number.networks<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	asr.type		<- 'inf'
	#asr.type		<- 'sam'
	outfile.base	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_'
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_summaries.rda'
	#
	#	get subtrees subst/site, inf ancestral state reconstruction, non-bootstrap 
	load(infile)	
	dss				<- subset(ds, TREE_TYPE!='lsd' & ASR_TYPE==asr.type)
	set(dss, NULL, 'SUBTREE_ID', dss[,paste0(ST,'-',REP,'-',SUBTREE_ID)])	
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_tips.rda'
	load(infile)	
	dti				<- subset(dti, TREE_TYPE!='lsd' & ASR_TYPE==asr.type)
	set(dti, NULL, 'ST', dti[,gsub('pol','',ST)])
	set(dti, NULL, 'SUBTREE_ID', dti[,paste0(ST,'-',REP,'-',SUBTREE_ID)])
	tmp				<- unique(subset(dss, select=c(SUBTREE_ID, ANCESTRAL_STATE, SUBTREE_SIZE)))
	dp				<- merge(tmp, dti, by='SUBTREE_ID')
	dp				<- merge(dp, dp[, list(MIN_SAMPLING_DATE=min(SAMPLING_DATE), MAX_SAMPLING_DATE=max(SAMPLING_DATE)), by='SUBTREE_ID'], by='SUBTREE_ID')
	dp				<- subset(dp, MIN_SAMPLING_DATE>1995)
	dp[, MIN_SAMPLING_DATE_C:= cut(MIN_SAMPLING_DATE, breaks=c(1930,1950,1970,1990,seq(1992,2016,2)))]
	dp[, STB:= as.character(factor(ST=='B', levels=c(TRUE,FALSE),labels=c('B','non-B')))]
	
	
	#	numbers overall for B and non-B		
	dc	<- dp[, list(N=length(unique(SUBTREE_ID))), by=c('STB','REP')]		
	tmp	<- dc[, list(QU=round(quantile(N, prob=c(0.025,0.25,0.5,0.75,0.975))), P=paste0('p',c(0.025,0.25,0.5,0.75,0.975))), by='STB']
	write.csv(tmp, row.names=FALSE, file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_counts_overall.csv'))
	
	#	numbers over time for B and non-B
	dc	<- dp[, list(N=length(unique(SUBTREE_ID))), by=c('STB','REP','MIN_SAMPLING_DATE_C')]
	dc	<- merge(as.data.table(expand.grid(STB=unique(dc$STB), REP=unique(dc$REP), MIN_SAMPLING_DATE_C=unique(dc$MIN_SAMPLING_DATE_C))), dc, by=c('STB','REP','MIN_SAMPLING_DATE_C'), all.x=TRUE)
	set(dc, dc[, which(is.na(N))], c('N'), 0)	
	tmp	<- dc[, list(QU=quantile(N, prob=c(0.025,0.25,0.5,0.75,0.975)), P=paste0('p',c(0.025,0.25,0.5,0.75,0.975))), by=c('STB','MIN_SAMPLING_DATE_C')]
	tmp	<- dcast.data.table(tmp, STB+MIN_SAMPLING_DATE_C~P, value.var='QU')	
	ggplot(tmp, aes(x=MIN_SAMPLING_DATE_C, y=p0.5, ymin=p0.025, ymax=p0.975)) +
			geom_bar(stat='identity', colour='grey70') +
			geom_errorbar(colour="black",width=0.5)	+
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous() +
			facet_grid(STB~., scales='free') +
			labs(x='earliest date of diagnosis\n(of first subject with a sequence that is part of the subtree)',y='new in-country subtrees', fill='geographic origin')			
	ggsave(file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_counts_time.pdf'), w=6, h=5)
	write.csv(tmp, file=paste0(outfile.base, asr.type, '_ALLST_statetransitions_counts_time.csv'))
}

vli.bez.explore.ancestral.state.reconstruction.171113.plot.distinct.networks<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	asr.type		<- 'inf'
	#asr.type		<- 'sam'
	outfile.base	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_'
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_summaries.rda'
	#
	#	get subtrees subst/site, inf ancestral state reconstruction, non-bootstrap 
	load(infile)	
	ds				<- subset(ds, TREE_TYPE!='lsd' & ASR_TYPE==asr.type & REP==0)
	set(ds, NULL, 'SUBTREE_ID', ds[,paste0(ST,'-',REP,'-',SUBTREE_ID)])	
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_tips.rda'
	load(infile)	
	dti				<- subset(dti, TREE_TYPE!='lsd' & ASR_TYPE==asr.type & REP==0)
	set(dti, NULL, 'ST', dti[,gsub('pol','',ST)])
	set(dti, NULL, 'SUBTREE_ID', dti[,paste0(ST,'-',REP,'-',SUBTREE_ID)])
	
	#
	#	
	tmp	<- unique(subset(ds, select=c(SUBTREE_ID, ANCESTRAL_STATE, INTRO_EARLIEST, INTRO_MID, INTRO_LATEST, SUBTREE_SIZE)))
	dp	<- merge(tmp, dti, by='SUBTREE_ID')
	dp	<- merge(dp, dp[, list(MIN_SAMPLING_DATE=min(SAMPLING_DATE), MAX_SAMPLING_DATE=max(SAMPLING_DATE)), by='SUBTREE_ID'], by='SUBTREE_ID')
	dp	<- subset(dp, MIN_SAMPLING_DATE>1995)
	set(dp, NULL, 'SXO', dp[, factor(SXO, levels=c('hsx','MSM','other','U'), labels=c('Heterosexual','MSM','Other','Unknown'))])
	set(dp, NULL, 'ANCESTRAL_STATE', dp[, factor(ANCESTRAL_STATE,	levels=c("SSA","EUC","EUO","EUW","LAT","OAP","ZAZ","NAM","USCAU","unsampled_region"),
																	labels=c("Sub-\nSaharan\nAfrica","East\nEU","Former\nUSSR","West\nEU","Latin\nAmerica","China\nJapan\nKorea","Sout and\nSoutheast\nAsia","North\nAfrica,\nMiddle East","US, Canada","unresolved"))
											])
	dpb	<- subset(dp, ST=='B')
	dpo	<- subset(dp, ST!='B')
	#		
	#	order non-B and B subtrees
	do	<- unique(subset(dpb, select=c(SUBTREE_ID, SUBTREE_SIZE, ST, ANCESTRAL_STATE, MIN_SAMPLING_DATE, MAX_SAMPLING_DATE)))
	setkey(do, ANCESTRAL_STATE, MIN_SAMPLING_DATE)
	dob	<- do[, list(SUBTREE_ID=SUBTREE_ID, SUBTREE_SIZE=SUBTREE_SIZE, MIN_SAMPLING_DATE=MIN_SAMPLING_DATE, MAX_SAMPLING_DATE=MAX_SAMPLING_DATE, ST=ST, SUBTREE_ORDER=seq_len(length(SUBTREE_ID))), by=c('ANCESTRAL_STATE')]
	dpb	<- merge(dpb, subset(dob, select=c(SUBTREE_ID, SUBTREE_ORDER)), by=c('SUBTREE_ID'))
	do	<- unique(subset(dpo, select=c(SUBTREE_ID, SUBTREE_SIZE, ST, ANCESTRAL_STATE, MIN_SAMPLING_DATE, MAX_SAMPLING_DATE)))
	setkey(do, ANCESTRAL_STATE, MIN_SAMPLING_DATE)
	doo	<- do[, list(SUBTREE_ID=SUBTREE_ID, SUBTREE_SIZE=SUBTREE_SIZE, MIN_SAMPLING_DATE=MIN_SAMPLING_DATE, MAX_SAMPLING_DATE=MAX_SAMPLING_DATE, ST=ST, SUBTREE_ORDER=seq_len(length(SUBTREE_ID))), by=c('ANCESTRAL_STATE')]
	dpo	<- merge(dpo, subset(doo, select=c(SUBTREE_ID, SUBTREE_ORDER)), by=c('SUBTREE_ID'))	
	
	if(0)
	{
		ggplot(dp) +			
				geom_point(aes(x=SAMPLING_DATE, y=SUBTREE_ORDER), colour='grey50') +
				geom_scatterpie(data=subset(do,SUBTREE_SIZE>1), aes(x=MIN_SAMPLING_DATE, y=SUBTREE_ORDER, group=SUBTREE_ORDER, r=log(SUBTREE_SIZE)/5+0.2), cols=colnames(do)[grepl('SUBTREE_SXO_',colnames(do))]) +
				#geom_point(data=tmp, aes(x=MIN_SAMPLING_DATE, y=SUBTREE_ORDER, colour=variable)) +			
				theme_bw() +
				scale_x_continuous(breaks=seq(1990,2020,2)) +
				labs(x='\ndate of diagnosis\n(patients with sequence)', y='phylogenetically observed transmission networks\n') +
				facet_grid(ANCESTRAL_STATE~., scales='free_y', space='free_y')
		ggsave(file=paste0(outfile.base,asr.type,'_ALLST_newtrees_over_time.pdf'),w=8,h=30)
		
	}
	if(1)	#non-B, colour = SXO
	{
		ggplot(dpo) +			
				geom_segment(data=subset(doo, SUBTREE_SIZE>1), aes(x=MIN_SAMPLING_DATE, xend=MAX_SAMPLING_DATE, y=SUBTREE_ORDER, yend=SUBTREE_ORDER), colour='black', size=1) +
				geom_point(aes(x=SAMPLING_DATE, y=SUBTREE_ORDER, colour=SXO)) +
				theme_bw() +
				theme(panel.spacing.y=unit(0.1, "lines"), strip.text.y=element_text(angle=0)) +
				scale_x_continuous(breaks=seq(1990,2020,2)) +	
				scale_y_continuous(expand=c(0,0)) +
				labs(title='Phylogenetically identified transmission networks among non-B sequences\n', x='\ndate of diagnosis\n(patients with sequence)', y='', colour='Sexual orientation of HIV-infected individuals\nin the Netherlands') +
				theme(legend.position='bottom') +
				facet_grid(ANCESTRAL_STATE~., scales='free_y', space='free_y')
		ggsave(file=paste0(outfile.base,asr.type,'_ALLST_newtrees_over_time_nonB.pdf'),w=8,h=20)
		ggplot(dpb) +			
				geom_segment(data=subset(dob, SUBTREE_SIZE>1), aes(x=MIN_SAMPLING_DATE, xend=MAX_SAMPLING_DATE, y=SUBTREE_ORDER, yend=SUBTREE_ORDER), colour='black', size=1) +
				geom_point(aes(x=SAMPLING_DATE, y=SUBTREE_ORDER, colour=SXO)) +
				theme_bw() +
				theme(panel.spacing.y=unit(0.1, "lines"), strip.text.y=element_text(angle=0)) +
				scale_x_continuous(breaks=seq(1990,2020,2)) +
				scale_y_continuous(expand=c(0,0)) +
				labs(title='Phylogenetically identified transmission networks among B sequences\n', x='\ndate of diagnosis\n(patients with sequence)', y='', colour='Sexual orientation of HIV-infected individuals\nin the Netherlands') +
				theme(legend.position='bottom') +
				facet_grid(ANCESTRAL_STATE~., scales='free_y', space='free_y')
		ggsave(file=paste0(outfile.base,asr.type,'_ALLST_newtrees_over_time_B.pdf'),w=8,h=20)
	}
	if(0)
	{
		ggplot(dp) +			
				geom_segment(data=subset(do, SUBTREE_SIZE>1), aes(x=MIN_SAMPLING_DATE, xend=MAX_SAMPLING_DATE, y=SUBTREE_ORDER, yend=SUBTREE_ORDER), colour='black', size=1) +
				geom_point(aes(x=SAMPLING_DATE, y=SUBTREE_ORDER, colour=ST)) +
				theme_bw() +
				theme(panel.spacing.y=unit(0.1, "lines"), strip.text.y=element_text(angle=0)) +
				scale_x_continuous(breaks=seq(1990,2020,2)) +
				scale_y_continuous(expand=c(0,0)) +
				labs(x='\ndate of diagnosis\n(patients with sequence)', y='phylogenetically observed transmission networks\n', colour='HIV-1 subtype') +				
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

vli.bez.figure.for.Daniela <- function()
{
	infile <- '~/Box/OR_Work/2019/2019_Daniela_paper/Fig2_amase_OR_1.csv'
	di <- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	
	#	data type
	tmp <- di[, which(grepl('ancestral state', data.set))]
	set(di, tmp, 'data_type', 'phylogenetic estimate, empirical')
	tmp <- di[, which(grepl('incomplete', data.set))]
	set(di, tmp, 'data_type', 'phylogenetic estimate, adjusted for partial sequence sampling')
	tmp <- di[, which(grepl('aMASE', data.set))]
	set(di, tmp, 'data_type', 'clinic survey, aMASE')
	tmp <- di[, which(grepl('Self-reported', data.set))]
	set(di, tmp, 'data_type', 'patient data, self-reported likely location of infection')
	
	#	sample population
	tmp <- di[, which(grepl('non-B', subtype))]
	set(di, tmp, 'sample_pop', 'patients in ATHENA cohort,\nsequenced,\nnon-B subtype')
	tmp <- di[, which(grepl('^B', subtype))]
	set(di, tmp, 'sample_pop', 'patients in ATHENA cohort,\nsequenced,\nsubtype B')
	tmp <- di[, which(grepl('Sub-Saharan', subtype))]
	set(di, tmp, 'sample_pop', 'aMASE\nsurvey participants, born in Sub-Saharan Africa')
	tmp <- di[, which(grepl('foreign born', subtype))]
	set(di, tmp, 'sample_pop', 'aMASE\nsurvey participants, foreign-born')
	
	#	facet_label
	di[, facet_label:=paste0(risk.group,'\n',paste0('foreign-born',gsub(', born in Sub-Saharan Africa|, foreign-born','',sample_pop)))]
	di[, facet_label:=gsub(', born in Sub-Saharan Africa|, foreign-born','',sample_pop)]
	tmp <- unique(subset(di, select=facet_label))
	tmp[, facet_label_2:= c(2,3,1)]	
	setkey(tmp, facet_label_2)
	tmp[, facet_label_2:= factor(facet_label_2, levels=facet_label_2, labels=facet_label)]
	di <- merge(di, tmp, by='facet_label')
	
	#	x label
	tmp <- di[, which(grepl('ancestral state', data.set))]
	set(di, tmp, 'x_label', 'phylogenetic estimate,\nempirical')
	set(di, tmp, 'x_label_idx', 1L)
	tmp <- di[, which(grepl('incomplete', data.set))]
	set(di, tmp, 'x_label', 'phylogenetic estimate,\nadjusted for partial sequence sampling')
	set(di, tmp, 'x_label_idx', 2L)
	tmp <- di[, which(grepl('aMASE', data.set) & grepl('foreign-born', sample_pop))]
	set(di, tmp, 'x_label', 'all participants')
	set(di, tmp, 'x_label_idx', 4L)
	tmp <- di[, which(grepl('aMASE', data.set) & grepl('Sub-Saharan', sample_pop))]
	set(di, tmp, 'x_label', 'born in Sub-Saharan Africa')
	set(di, tmp, 'x_label_idx', 5L)
	tmp <- di[, which(grepl('Self-reported', data.set))]
	set(di, tmp, 'x_label', 'self-reported patient data')
	set(di, tmp, 'x_label_idx', 3L)
	tmp <- unique(subset(di, select=c(x_label,x_label_idx)))
	setkey(tmp, x_label_idx)
	set(tmp, NULL, 'x_label_2', tmp[, factor(x_label_idx, levels=x_label_idx, labels=x_label)])
	di <- merge(di, tmp, by=c('x_label_idx','x_label'))
	
	set(di, NULL, 'mid.point', di[, mid.point/100])
	set(di, NULL, 'lower', di[, lower/100])
	set(di, NULL, 'upper', di[, upper/100])
	
	ggplot(di) +
			geom_bar(aes(x=x_label_2, y= mid.point, fill= x_label), stat='identity') +
			geom_errorbar(aes(x=x_label_2, ymin=lower, ymax=upper), width=0.5) +
			scale_y_continuous(label=scales::percent, limits=c(0,1), expand=c(0,0)) +
			scale_fill_brewer(palette='Set2') +
			theme_bw() +		
			theme(legend.position = "none",
					panel.spacing = unit(1, "lines"),
					axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
			facet_grid(risk.group~facet_label_2, scales='free_x', space = "free") +
			labs(fill='sample population', y='HIV acquisitions originating in the Netherlands', x='')
	ggsave(file='~/Box/OR_Work/2019/2019_Daniela_paper/Fig2_amase_OR_1.pdf', w=6, h=7)
	
}

vli.bez.summarize.tip.data<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	outfile.save	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_tips.rda'	
	infile.meta	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo_plusB_OR.csv'
	dm	<- as.data.table(read.csv(infile.meta))
	set(dm, NULL, 'SXO', dm[, gsub('HTF|HTM','hsx',gsub('BL|CH|CHS|DU','other',as.character(GENDER_TRMGROUP)))])	
	indir.asr	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results'
	infiles		<- data.table(F=list.files(indir.asr, pattern='^workspace', full.names=TRUE, recursive=TRUE))
	infiles[, ASR_TYPE:= gsub('.*_annotated_([a-z]+).*','\\1',basename(F))]
	infiles[, REP:= as.integer(gsub('^.*_([0-9]+)_annotated.*$|^.*_rr([0-9]+)_annotated.*$','\\1',basename(F)))]
	tmp			<- infiles[, which(is.na(REP))]
	set(infiles, tmp, 'REP', infiles[tmp,as.integer(gsub('^.*_rr([0-9]+)_annotated.*$','\\1',basename(F)))])	
	infiles[, ST:= gsub('pol','',gsub('^workspace_([^_]+)_.*','\\1',basename(F)))]
	infiles[, TREE_TYPE:=gsub('^.*(lsd|rr).*$','\\1',basename(F))]
	set(infiles, infiles[, which(TREE_TYPE!='lsd')],'TREE_TYPE','substsite')
	infiles[, SELECT_DESC:= as.character(factor(ASR_TYPE, levels=c('bir','sam','inf'), labels=c('^.*_bornNL_.*$','^.*_sampNL_.*$','^.*_infNL_.*$')))]
	infiles[, DEADEND_DESC:= NA_character_]
	infiles[, IS_DATED_TREE:= TREE_TYPE=='lsd']
	set(infiles, infiles[, which(ASR_TYPE=='inf')], 'DEADEND_DESC',  "^.*_sampNL_.*$|^.*_infNL_.*$")		
	infiles		<- subset(infiles, !IS_DATED_TREE)
	id.regex	<- '^([^_]*).*'
			
	dti			<- infiles[, {
				#	IS_DATED_TREE<- 1; DEADEND_DESC<- '^.*_sampNL_.*$|^.*_infNL_.*$'; SELECT_DESC<- '^.*_infNL_.*$'; F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results/B/sam/workspace_polB_SHM_LANL_reduced149_withD_withHXB2_e_noDRM_rr_095_annotated_sam.rda'
				#	IS_DATED_TREE<- 0; DEADEND_DESC<- '^.*_sampNL_.*$|^.*_infNL_.*$'; SELECT_DESC<- '^.*_infNL_.*$'; F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results/A1/inf/workspace_A1_withC_rr000_annotated_inf.rda'
				#	IS_DATED_TREE<- 1; DEADEND_DESC<- NA_character_; SELECT_DESC<- '^.*_saNL_.*$'; F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results/B/sam/workspace_polB_SHM_LANL_reduced149_withD_withHXB2_e_noDRM_rr_095_annotated_sam.rda'
				select.descendants.regex	<- SELECT_DESC
				is.dated.tree				<- IS_DATED_TREE	
				deadend.descendants.regex	<- DEADEND_DESC						
				cat('\nprocessing file ',F,'\n',select.descendants.regex,' ',deadend.descendants.regex,' ',is.dated.tree)
				load(F)								
				#	read variables from workspace
				ph	<- all.tree.info$only.tree$tree
				da	<- as.data.table(all.tree.info$only.tree$classification.results$collapsed)	
				#	summarise introductions into NL
				dti	<- vli.asr.get.tips.in.NL.subtrees(ph, da, dm, id.regex, select.descendants.regex=select.descendants.regex, deadend.descendants.regex=deadend.descendants.regex)				
				dti
			}, by=c('TREE_TYPE','ST','ASR_TYPE','REP')]	
	save(dti, file=outfile.save)
	write.csv(dti, row.names=FALSE, file=gsub('rda$','csv',outfile.save))
}

vli.bez.evaluate.infection.loc.missing<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	infile			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/Table1_infoyesno.csv'
	df				<- as.data.table(read.csv(infile))
	df				<- merge(df, df[, list(EQ=EQ, P=COUNT/sum(COUNT)), by=c('CAT', 'haveinfo')], by=c('CAT', 'haveinfo', 'EQ'))
	dcast.data.table(df, CAT~haveinfo+EQ, value.var='P')
	
}

vli.bez.central.subgraphs.v190723<- function()
{
	require(data.table)
	require(Gmedian)
	require(pracma)
	require(OjaNP)
	
	#	write to file
	indir <- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis'
	infile.subgraphs <- file.path(indir,'190801_subgraph_data.csv')
	infile.subgraphtaxa <- file.path(indir,'190801_individuals_in_subgraph_data.csv')
	infile.subgraphs <- file.path(indir,'190801_subgraph_data_national.csv')
	infile.subgraphtaxa <- file.path(indir,'190801_individuals_in_subgraph_data_national.csv')
	
	dsubgraphs <- as.data.table(read.csv(infile.subgraphs, stringsAsFactors=FALSE))
	dsubgraphtaxa <- as.data.table(read.csv(infile.subgraphtaxa, stringsAsFactors=FALSE))
	
	#	count number of subgraphs by size
	dsizes <- dsubgraphs[, list(SIZE_N= length(NAME)), by=c('ST','REP','TRM_GROUP','SIZE')]
	tmp <- dsubgraphs[, list(MAX_SIZE=max(SIZE)), by=c('ST','TRM_GROUP')]
	dsizes <- merge(dsizes, tmp, by=c('ST','TRM_GROUP'))
	
	#	for subtype and trm group, put into matrix
	dsizes <- split(dsizes, by=c('ST','TRM_GROUP'))
	dsmed <- lapply( seq_along(dsizes), function(i)
			{
				#i<- 1
				print(i)
				ds <- dsizes[[i]]
				ans <- data.table(NAME= names(dsizes)[i], REP=0)
				if(ds$MAX_SIZE[1]>1)
				{
					ds <- ds[, {
								ss <- rep(0, MAX_SIZE[1])
								ss[SIZE] <- SIZE_N
								list(SIZE=seq_along(ss), SIZE_N=ss)
							}, by='REP']
					ds <- as.matrix(dcast.data.table(ds, REP~SIZE, value.var='SIZE_N'))[,-1]
					dsmed <- Gmedian(ds)
					dist.to.med <- apply( abs(ds - matrix(dsmed, nrow=nrow(ds), ncol=ncol(ds), byrow=TRUE)), 1, sum)
					median.weiszfeld <- which.min(dist.to.med)
					ans <- data.table(NAME= names(dsizes)[i], REP=median.weiszfeld-1L)	
				}
				ans
			})
	dsmed <- do.call('rbind',dsmed)
	dsmed[, ST:= gsub('^([A-Za-z0-9]+)\\.([A-Za-z0-9]+)$','\\1', NAME)]
	dsmed[, TRM_GROUP:= gsub('^([A-Za-z0-9]+)\\.([A-Za-z0-9]+)$','\\2', NAME)]
	set(dsmed, NULL, 'NAME', NULL)
	
	#	write results of central subgraphs to file
	dsubgraphs.central <- merge(dsubgraphs, dsmed, by=c('ST','TRM_GROUP','REP'))
	dsubgraphtaxa.central <- merge(dsubgraphtaxa, dsmed, by=c('ST','TRM_GROUP','REP'))
	
	#	write to file
	outfile.base <- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/190801_'
	write.csv(dsubgraphs.central, row.names=FALSE, file=paste0(outfile.base,'centralsubgraph_data_national.csv'))
	write.csv(dsubgraphtaxa.central, row.names=FALSE, file=paste0(outfile.base,'individuals_in_centralsubgraph_data_national.csv'))
}

vli.bez.summarize.state.transitions.v190723<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	require(reshape)
	require(phyloscannerR)
	
	#	working directory with phyloscanner output		
	indir.phsc	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results_190723_bytrmgroup'
	indir.phsc	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results_190723_national'
	#	meta data
	infile.meta	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo_NLbyTrmGroup.csv'
	#	file with individuals with identical sequence
	infile.identical <- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/addtocluster_190801.txt'
	
	#
	#	extract subgraphs
	#	(only needs to be done once, skip thereafter)
	#
	infiles		<- data.table(F=list.files(indir.phsc, pattern='workspace.rda$', full.names=TRUE, recursive=TRUE))
	#infiles[, SELECT_DESC_HSX:= 'NLHSX']
	#infiles[, SELECT_DESC_MSM:= 'NLMSM']
	#infiles[, SELECT_DESC_IDU:= 'NLIDU']
	infiles[, SELECT_DESC_ALL:= 'NL']
	infiles 	<- melt(infiles, measure.vars=c('SELECT_DESC_ALL'), value.name='SELECT', variable.name='SELECT_TYPE')
	#infiles 	<- melt(infiles, measure.vars=c('SELECT_DESC_HSX','SELECT_DESC_MSM','SELECT_DESC_IDU'), value.name='SELECT', variable.name='SELECT_TYPE')
	#infiles 	<- subset(infiles, grepl('^G', basename(F)))
	for(i in seq_len(nrow(infiles)))
	{
		#i<- 1
		cat('process', i,'\n')
		infile <- infiles[i, F]
		host <- infiles[i,SELECT]
		load(infile)	
		ph <- phyloscanner.trees[[1]][['tree']]		
		mrcas <- which( attr(ph, 'SUBGRAPH_MRCA') )
		mrcas <- mrcas[ attr(ph, 'INDIVIDUAL')[mrcas]==host ]	
		# convert tree to class simmap
		ph <- phyloscanner.to.simmap(ph)
		# extract subgraphs
		subgraphs <- lapply(mrcas, function(mrca) extract.subgraph(ph, mrca))
		# save
		outfile <- gsub('__workspace',paste0('__subgraphs_',host),infile)
		save(subgraphs, file=outfile)
	}
	
	#
	#	get origin for all subgraphs 
	#
	infiles		<- data.table(F=list.files(indir.phsc, pattern='^.*_subgraphs_.*$', full.names=TRUE, recursive=TRUE))
	infiles[, ST:= gsub('pol','',gsub('^([^_]+)_.*','\\1',basename(F)))]	
	infiles[, REP:= as.integer(gsub('^.*([0-9]{3})_([a-z]+)__.*$','\\1',basename(F)))]
	infiles[, TRM_GROUP:= gsub('^.*_subgraphs_([A-Z]+)\\.rda$','\\1',basename(F))]	
	dsubgraphs <- infiles[, {
				#i<- 1
				#infile <- infiles[i, F]			
				infile <- F
				cat('Process',infile,'\n')
				load(infile)				
				subgraph.names <- sapply( subgraphs, function(subgraph)  subgraph$subgraph.name )
				subgraph.parent.state <- sapply( subgraphs, function(subgraph)  subgraph$subgraph.parent.state )
				subgraph.parent.state[is.na(subgraph.parent.state)] <- 'Unknown'
				list(	NAME=subgraph.names, 						 
						PARENT_STATE=subgraph.parent.state)				
			}, by=c('ST','REP','TRM_GROUP')]
	#set(dsubgraphs, which(dsubgraphs[, is.na(PARENT_STATE)]), 'PARENT_STATE', 'Unknown')
	
	#
	#	get taxa for all subgraphs 
	#
	dsubgraphtaxa <- infiles[, {
			#i<- 1
			#infile <- infiles[i, F]			
			infile <- F
			cat('Process',infile,'\n')
			load(infile)			
			subgraph.names <- unlist(sapply( subgraphs, function(subgraph)  rep(subgraph$subgraph.name, length(subgraph$tip.label)) ))
			subgraph.taxa <- unlist(sapply( subgraphs, function(subgraph)  subgraph$tip.label))
			subgraph.taxa <- gsub('^([A-Za-z0-9]+)_.*$','\\1', subgraph.taxa)
			list(	NAME=subgraph.names, 
					ID= subgraph.taxa 
					)				
		}, by=c('ST','REP','TRM_GROUP')]	
	#	add individuals with identical sequence that were left out from phylo analysis
	tmp <- as.data.table(read.table(infile.identical, header=TRUE, stringsAsFactors=FALSE))
	tmp <- rbind( tmp, data.table(In_cluster=tmp$In_cluster, Add_to_cluster=tmp$In_cluster) )
	tmp <- unique(tmp)
	setnames(tmp, c('In_cluster','Add_to_cluster'), c('ID','ID2'))
	dsubgraphtaxa <- merge(dsubgraphtaxa, tmp, by='ID',all.x=TRUE,allow.cartesian=TRUE)
	tmp <- which(is.na(dsubgraphtaxa[, ID2]))
	set(dsubgraphtaxa, tmp, 'ID2', dsubgraphtaxa[tmp,ID])
	set(dsubgraphtaxa, NULL, 'ID', NULL)
	setnames(dsubgraphtaxa,'ID2','ID') 	
	#	add meta data
	dm	<- as.data.table(read.csv(infile.meta, stringsAsFactors=FALSE))
	set(dm, dm[, which(grepl('^NL',SAMPLING_LOC))], 'SAMPLING_LOC', 'NL')
	set(dm, dm[, which(grepl('^NL',INFECTED_LOC))], 'INFECTED_LOC', 'NL')
	set(dm, dm[, which(grepl('^NL',BORN_LOC))], 'BORN_LOC', 'NL')
	dsubgraphtaxa <- merge(dsubgraphtaxa, dm, by='ID')
	
	#	add first/last sampling date and size for each subgraph to 'dsubgraphs'
	tmp <- dsubgraphtaxa[,  list(MIN_SAMPLING_DATE=min(SAMPLING_DATE), MAX_SAMPLING_DATE=max(SAMPLING_DATE), SIZE=length(ID)), by=c('ST','REP','TRM_GROUP','NAME')]
	dsubgraphs <- merge(dsubgraphs,tmp,by=c('ST','REP','TRM_GROUP','NAME'))
	
	#	write to file
	outfile.base <- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/190801_'
	write.csv(dsubgraphs, row.names=FALSE, file=paste0(outfile.base,'subgraph_data_national.csv'))
	write.csv(dsubgraphtaxa, row.names=FALSE, file=paste0(outfile.base,'individuals_in_subgraph_data_national.csv'))
}

vli.bez.summarize.state.transitions.v171019<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	
	outfile.save	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/analysis/171109_asr_summaries.rda'	
	infile.meta	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo_plusB_OR.csv'
	dm	<- as.data.table(read.csv(infile.meta))
	set(dm, NULL, 'SXO', dm[, gsub('HTF|HTM','hsx',gsub('BL|CH|CHS|DU','other',as.character(GENDER_TRMGROUP)))])
	
	indir.asr	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results'
	infiles		<- data.table(F=list.files(indir.asr, pattern='^workspace', full.names=TRUE, recursive=TRUE))
	infiles[, ASR_TYPE:= gsub('.*_annotated_([a-z]+).*','\\1',basename(F))]
	infiles[, REP:= as.integer(gsub('^.*_([0-9]+)_annotated.*$|^.*_rr([0-9]+)_annotated.*$','\\1',basename(F)))]
	tmp			<- infiles[, which(is.na(REP))]
	set(infiles, tmp, 'REP', infiles[tmp,as.integer(gsub('^.*_rr([0-9]+)_annotated.*$','\\1',basename(F)))])
	infiles[, ST:= gsub('pol','',gsub('^workspace_([^_]+)_.*','\\1',basename(F)))]
	infiles[, TREE_TYPE:=gsub('.*_(lsd|rr)_.*','\\1',basename(F))]
	set(infiles, infiles[, which(TREE_TYPE!='lsd')],'TREE_TYPE','substsite')
	infiles[, SELECT_DESC:= as.character(factor(ASR_TYPE, levels=c('bir','sam','inf'), labels=c('^.*_bornNL_.*$','^.*_sampNL_.*$','^.*_infNL_.*$')))]
	infiles[, DEADEND_DESC:= NA_character_]
	infiles[, IS_DATED_TREE:= TREE_TYPE=='lsd']
	set(infiles, infiles[, which(ASR_TYPE=='inf')], 'DEADEND_DESC',  "^.*_sampNL_.*$|^.*_infNL_.*$")		
	infiles		<- subset(infiles, !IS_DATED_TREE)
	id.regex	<- '^([^_]*).*'
	
	ds			<- infiles[, {
				#	IS_DATED_TREE<- 0; DEADEND_DESC<- "^.*_sampNL_.*$|^.*_infNL_.*$"; SELECT_DESC<- '^.*_infNL_.*$'; F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2017_NL_Introductions/phyloscanner_results/B/inf/workspace_polB_SHM_LANL_reduced149_withD_withHXB2_e_noDRM_rr_015_annotated_inf.rda'				
				select.descendants.regex	<- SELECT_DESC
				is.dated.tree				<- IS_DATED_TREE	
				deadend.descendants.regex	<- DEADEND_DESC						
				cat('\nprocessing file ',F,'\n',select.descendants.regex,' ',deadend.descendants.regex,' ',is.dated.tree)
				load(F)				
				#	read variables from workspace
				ph	<- all.tree.info$only.tree$tree
				da	<- as.data.table(all.tree.info$only.tree$classification.results$collapsed)	
				#	summarise introductions into NL
				ds	<- vli.asr.summary(ph, da, dm, id.regex, select.descendants.regex=select.descendants.regex, deadend.descendants.regex=deadend.descendants.regex, is.dated.tree=is.dated.tree)	
				ds[, INTRO_MID:= (INTRO_EARLIEST+INTRO_LATEST)/2]
				setkey(ds, INTRO_MID)
				ds
			}, by=c('TREE_TYPE','ST','ASR_TYPE','REP')]	
	save(ds, file=outfile.save)
	write.csv(ds, row.names=FALSE, file=gsub('rda$','csv',outfile.save))
	
}

######################################################################################

vli.bez.DataFile.170901<- function()
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

vli.bez.DataFile.171019<- function()
{
	require(data.table)
	require(ggplot2)
	infile	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo_plusB.csv'
	indir	<- dirname(infile)
	basename<- file.path(indir,gsub('\\.csv','',basename(infile)))
	
	df		<- as.data.table(read.csv(infile))
	df1		<- subset(df, select=c(ID, GENDER, GENDER_TRMGROUP, DATA_FROM_SHM, subtype, SAMPLING_DATE, SAMPLING_LOC_V1, BORN_LOC_V1, INFECTED_LOC_V1, SAMPLING_COUNTRY, INFECTED_COUNTRY, BORN_COUNTRY))
	setnames(df1, 	c('subtype','SAMPLING_LOC_V1','BORN_LOC_V1','INFECTED_LOC_V1'), 
					c('SUBTYPE','SAMPLING_LOC','BORN_LOC','INFECTED_LOC'))
	
	stopifnot( 	df1[, !any(is.na(SAMPLING_LOC))], 
				df1[, all(unique(SAMPLING_LOC)%in%c('EUC','EUO','EUW','LAT','NAM','NL','OAP','SSA','USCAU','ZAZ','XX'))] )
	stopifnot( 	df1[, !any(is.na(BORN_LOC))], 
				df1[, all(unique(BORN_LOC)%in%c('EUC','EUO','EUW','LAT','NAM','NL','OAP','SSA','USCAU','ZAZ','XX'))] )
	stopifnot( 	df1[, !any(is.na(INFECTED_LOC))], 
				df1[, all(unique(INFECTED_LOC)%in%c('EUC','EUO','EUW','LAT','NAM','NL','OAP','SSA','USCAU','ZAZ','XX'))] )
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
