amsterdam.200306.alignment.make.bootstraps <- function()
{   
  home <- '/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223'
  indir <- file.path(home,'alignments')
  outdir <- file.path(home,'alignments_bs')
  bsn <- 10
  seed <- 42L
  #
  df <- data.table(F=list.files(indir, pattern='fasta$',full.names=TRUE))
  df <- subset(df, grepl('noDRM',F))
  for(i in seq_len(nrow(df)))
  {
    infile.fasta <- df[i,F]
    cat('\nprocess ',infile.fasta)
    sq                  <- read.dna(infile.fasta, format='fa')  
    set.seed(seed)
    for(b in 0:bsn)
    {
      outfile <- gsub('\\.fasta',paste0('_',sprintf("%03d",b),'.fasta'),infile.fasta)
      outfile <- file.path(outdir,basename(outfile))
      bs <- 1:ncol(sq)
      if(b>0)
      {
        bs<- sample(1:ncol(sq), ncol(sq), replace=TRUE) 
      }           
      write.dna(sq[,bs], file=outfile, format='fa', colsep='', nbcol=-1)
    }       
  }   
}

amsterdam.200306.fastree<- function()
{   
  #require(big.phylo)
  require(data.table)
  
  #home <- '/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223'
  home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'

  indir <- file.path(home,'alignments_bs')
  outdir <- file.path(home,'fasttree')
  
  #   get UNIX commands for all trees to be processed
  df <- data.table(FIN=list.files(indir, pattern='fasta$',full.names=TRUE))
  df <- subset(df, grepl('noDRM',FIN))
  df[, ST:= gsub('.*_subtype_([A-Z0-9a-z]+)_.*','\\1',basename(FIN))]
  df[, FOUT:= file.path(outdir,gsub('\\.fasta','_ft.newick',basename(FIN)))]
  cmds <- vector('list', nrow(df))
  for(i in seq_len(nrow(df)))
  {
    infile.fasta <- df[i,FIN]
    outfile <- df[i,FOUT]           
    tmp <- cmd.fasttree(infile.fasta, outfile=outfile, pr.args='-nt -gtr -gamma', check.binary=TRUE)
    cmds[[i]] <- tmp
  }
  df[, CMD:= unlist(cmds)]    
  
  #   submit jobs like this one:
  cat(df[1,CMD])
  
  #   run on HPC as array job
  #df <- subset(df, ST=='B')
  df[, CASE_ID:= 1:nrow(df)]
  
  #   make PBS header
  #hpc.load    <- "module load R/3.4.0"
  hpc.load    <- "module load anaconda3/personal"
  r.activate  <- "source activate Renv"
  hpc.select  <- 1                        # number of nodes
  hpc.nproc   <- 1                        # number of processors on node
  hpc.walltime<- 72                      # walltime
  hpc.q       <- NA                       # PBS queue
  hpc.mem     <- "6gb"                    # RAM
  hpc.array   <- length(unique(df$CASE_ID))   # number of runs for job array
  pbshead     <- "#!/bin/sh"
  tmp         <- paste("#PBS -l walltime=", hpc.walltime, ":00:00,pcput=", hpc.walltime, ":00:00", sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  tmp         <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  pbshead     <- paste(pbshead, "#PBS -j oe", sep = "\n") 
  if(!is.na(hpc.array))
    pbshead <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')        
  if(!is.na(hpc.q)) 
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead     <- paste(pbshead, hpc.load, r.activate, sep = "\n")         
  #   make array job
  cmd     <- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
  cmd     <- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]         
  cmd     <- paste(pbshead,cmd,sep='\n')  
  #   submit job
  outfile     <- gsub(':','',paste("trs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
  outfile     <- file.path(outdir, outfile)
  cat(cmd, file=outfile)
  cmd         <- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))
}

amsterdam.200312.sequence.labels<- function()
{
  require(data.table)
  require(ape)
  require(tidyverse)
  
  home <- '/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223'
  home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
  infile.indinfo <- file.path(home,'data_191223_Amsterdam','SHM_1902_ROADMAP_191223_tblBAS.csv')
  infile.nlinfo <- file.path(home,'data_200316_Netherlands','SHM_1902_ROADMAP_200316_tblBAS.csv')
  infile.subtypes <- file.path(home,'misc','ROADMAP_200319_All_Taxa_withsubtype.rda')
  infile.seqinfo <- file.path(home,'data_191223_Amsterdam','SHM_1902_ROADMAP_191223_tblLAB_seq.rda')   
  infile.nlseqinfo <- file.path(home,'data_200316_Netherlands','SHM_1902_ROADMAP_200316_tblLAB_seq.rda')   
  infile.georeg <- file.path(home,'misc','UN_geographic_locations.csv')
  infile.countrycodes <- file.path(home,'misc','ISO_country_codes.csv')
  infiles.lanl <- file.path(home,'alignments','data_200305_LANL_alignment.fasta')
  outfile.base <- file.path(home,'misc','200319_')
  
  #
  # collect country codes - from UN (geographic regions) and ISO.org (for country codes to match)
  #
  headers = read.csv(infile.georeg, header = F, nrows = 1, as.is = T)
  dgeo <- read.csv(infile.georeg,header=F,skip=26)
  colnames(dgeo)= headers
  names(dgeo)[5] <- 'Alpha.3.code'
  dgeo <- dgeo[,c('Location','Alpha.3.code','GeoRegName')]
  iso <- read.csv(infile.countrycodes,header=T)
  iso <- merge(iso,dgeo,by='Alpha.3.code',all=T)
  iso <- iso[!is.na(iso[,'Alpha.2.code']),]  
  iso  <- iso %>% mutate(Alpha.2.code:=as.character(Alpha.2.code),
                         WRLD:= as.character(GeoRegName),
                         CNTRY:= as.character(Location)) %>%
    mutate(WRLD:= case_when(WRLD=='Latin America and the Caribbean'~'LaAmCarr',
                            WRLD=='Northern America'~'NorthAm',
                            WRLD!='Latin America and the Caribbean'&WRLD!='Northern America'~WRLD),
          CNTRY:= case_when(CNTRY=='United States of America'~'USA',
                            CNTRY=='United Kingdom'~'UK',
                            CNTRY!='United States of America'&CNTRY!='United Kingdom'~CNTRY))  
  
  #
  # read Amsterdam individual data and merge with Dutch data
  #   
  dind <- read.csv(infile.indinfo,header=T) 
  dnl <- read.csv(infile.nlinfo,header=T) 
  dind$city <- 'Amsterdam'
  dnl$city <- 'Non-Amsterdam'
  dind <- merge(dind,dnl,all=T)

  # Recode baseline variables
  dind <- dind %>% 
    select(PATIENT,GENDER,MODE,MODE_OTH,ORIGIN,MIG_D,MIG_D_pq,MIG_D_aq,INF_NL,INF_COUNTRY_1,city) %>%
    mutate( gender2:= case_when(GENDER==1~'Male',
                                GENDER==2~'Female'),
            transm:= case_when(MODE==1~'MSM',
                              MODE==2~'IDU',
                              MODE==4~'Other',
                              MODE==5~'Other',
                              MODE==6~'HSX',
                              MODE==8~'Other',
                              MODE==90~'Other',
                              MODE==99~'Other'),
            country:= case_when(ORIGIN==""~'Unknown',
                                ORIGIN=="NL"~'Dutch',
                                ORIGIN!='NL'&ORIGIN!=""~'Other'),
            infcountry:= case_when(INF_COUNTRY_1=='NL'~'Dutch',
                                   INF_COUNTRY_1=='-1'~'Unknown',
                                   INF_COUNTRY_1==""~'Unknown',
                                   INF_COUNTRY_1!='NL'&INF_COUNTRY_1!='-1'&INF_COUNTRY_1!=""~'Non-NL'
                               )) %>%
    mutate( Alpha.2.code:= as.character(ORIGIN))
  
  # Obtain world region for country of origin and recode a few manually not in geog regions data
  dind <- dind %>% left_join(iso,by='Alpha.2.code') %>%
          mutate(CNTRY:= case_when(Alpha.2.code=='AN'~'Netherlands Antilles',
                                   Alpha.2.code=='YU'~'Yugoslavia',
                                   Alpha.2.code=='CS'~'Serbia and Montenegro',
                                   country=='Unknown'~'Unknown',
                                   country==NA~'Unknown',
                                   TRUE ~ CNTRY),
                 WRLD:= case_when(Alpha.2.code=='AN'~'LaAmCarr',
                                  Alpha.2.code=='YU'~'Europe',
                                  Alpha.2.code=='CS'~'Europe',
                                  country=='Unknown'~'Unknown',
                                  country==NA~'Unknown',
                                  TRUE ~ CNTRY))
  # Read in subtype file
  st <- load(infile.subtypes)
  st.l <- ds %>% select(TAXA_L,SUBTYPE_L) %>% distinct()
  st.a <- ds %>% select(FASTASampleCode,SUBTYPE) %>% distinct()
  colnames(st.a)[1] <- 'SEQ_LABEL'
  
  #
  # read sequence labels for Amsterdam and NL and add in subtype and geographical data
  #   
  load(infile.seqinfo)
  dseq <- ds
  load(infile.nlseqinfo)
  dseq <- rbind(dseq,ds)
  dseq <- dseq %>% inner_join(st.a, by='SEQ_LABEL')  
  
  dseq <- dind %>% inner_join(dseq, by='PATIENT')
  dseq <- dseq %>% select(PATIENT, SEQ_LABEL, SEQ_ID, SEQ_D, SEQ_L, SUBTYPE, ORIGIN, gender2, transm, city, CNTRY, WRLD) %>%
          rename(SEQ_DATE=SEQ_D, SEQ_LENGTH=SEQ_L, GENDER=gender2, TRANSM=transm, CITY=city, BIRTHCOUNTRY=CNTRY)

  #
  # collect LANL labels
  #
  dl <- read.dna(infiles.lanl,format='fa')
  dl <- tibble(TAXA=unlist(rownames(dl))) %>% 
    filter(!grepl('-.-.-.',TAXA))  %>% 
    mutate( SUBTYPE:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',1),
            Alpha.2.code:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',2),
            YEAR:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',3),
            GENBANK:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',5)) %>% 
    distinct()
  dl$TAXA[grepl('K03455',dl$TAXA)] <- 'HXB2'
  
  outfile <- paste0(outfile.base,'sequence_labels.rda')
  save(dseq, dind, dl, iso, file= outfile)
  
  #
  #   make sequence labels for AmsMSM
  #
  tmp <- dseq %>% rename(TAXA:= SEQ_LABEL) %>%
    mutate(GRP:= case_when(CITY=='Amsterdam'&TRANSM=='MSM'~'AmsMSM',
                           CITY=='Amsterdam'&TRANSM!='MSM'~'AmsnonMSM',
                           CITY!='Amsterdam'~'NL')) %>%
    mutate(TAXA_NEW:= paste0(GRP,'___',SEQ_ID,'_',PATIENT,'_',CITY,'_',ORIGIN,'_',GENDER,'_',TRANSM,'_',SEQ_DATE,'_',SEQ_LENGTH)) %>%
    select(SUBTYPE,WRLD,TAXA,TAXA_NEW)
  tmp <- dl %>% left_join(iso, by='Alpha.2.code') %>%
    mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA)) %>%
    select(SUBTYPE,WRLD,TAXA,TAXA_NEW)  %>%
    rbind(tmp)  %>%
    arrange(WRLD)
  outfile <- paste0(outfile.base,'sequence_labels_AmsMSM.csv')
  write.csv(tmp, file=outfile)
  #
  #   make sequence labels for AmsHSX
  #
  tmp <- dseq %>% rename(TAXA:= SEQ_LABEL) %>%
    mutate(GRP:= case_when(CITY=='Amsterdam'&TRANSM=='HSX'~'AmsHSX',
                            CITY=='Amsterdam'&TRANSM!='HSX'~'AmsnonHSX',
                            CITY!='Amsterdam'~'NL')) %>%
    mutate(TAXA_NEW:= paste0(GRP,'___',SEQ_ID,'_',PATIENT,'_',CITY,'_',ORIGIN,'_',GENDER,'_',TRANSM,'_',SEQ_DATE,'_',SEQ_LENGTH)) %>%
    select(SUBTYPE,WRLD,TAXA,TAXA_NEW)
  tmp <- dl %>% left_join(iso, by='Alpha.2.code') %>%
    mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA)) %>%
    select(SUBTYPE,WRLD,TAXA,TAXA_NEW) %>%
    rbind(tmp) %>%
    arrange(WRLD)
  outfile <- paste0(outfile.base,'sequence_labels_AmsHSX.csv')
  write.csv(tmp, file=outfile)
  #
  #   make sequence labels for Amsterdam overall
  #
  tmp <- dseq %>% rename(TAXA:= SEQ_LABEL) %>%
    mutate(GRP:= case_when(CITY=='Amsterdam'~'Ams',                            
                            CITY!='Amsterdam'~'NL')) %>%
    mutate(TAXA_NEW:= paste0(GRP,'___',SEQ_ID,'_',PATIENT,'_',CITY,'_',ORIGIN,'_',GENDER,'_',TRANSM,'_',SEQ_DATE,'_',SEQ_LENGTH)) %>%
    select(SUBTYPE,WRLD,TAXA,TAXA_NEW)
  tmp <- dl %>% left_join(iso, by='Alpha.2.code') %>%
    mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA)) %>%
    select(SUBTYPE,WRLD,TAXA,TAXA_NEW) %>%
    rbind(tmp) %>%
    arrange(WRLD)
  outfile <- paste0(outfile.base,'sequence_labels_Ams.csv')
  write.csv(tmp, file=outfile)
}


roadmap.200203.estimate.transmissions.from.external.branchingprocess <- function()
{
	require(rstan)
	require(bayesplot)
	require(hexbin)
	
	args <- list()
	args$size.inf.pop <- 5e3
	args$size.inf.pop.sampled <- 3e3
	args$upper.bound.multiplier <- 10
	args$indir <- '~/Box/OR_Work/AIDSFonds/analysis_200407/subgraphs'
	args$infile.subgraph.sizes <- '200512_subgraphsizes_hsx.csv'
	args$outdir <- '~/Box/OR_Work/AIDSFonds/analysis_200407/subgraphs'
	
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
	real log_rho;
	real log_one_minus_rho;	
	real log_r0_div_kappa;
	real log_vmr;	
	matrix[N_cs_actual, 5] tmp;
	vector[5] ones;
	vector[N_cs_actual] cs_actual_lpmf;
	vector[N_cs_obs+1] cs_obs_lpmf;	
	matrix[N_cs_obs+1, N_cs_actual+1] bin_lpmf;	

	// define transformed parameters
	kappa = r0/vmr_minus_one;
	log_one_minus_rho = log( 1 - rho );
	log_rho = log( rho );	
	log_r0_div_kappa = log( vmr_minus_one );
	log_vmr = log( 1+vmr_minus_one );

	// calculate lpmf of actual chain sizes given R0 and kappa
	// TODO I am not clear if Prob( size of actual chain=0 ) is well defined. 
	//		Code below is assuming this is zero.
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
	{
		int tmp_int;	// must declare index integer inside block
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

model{
	// priors
	target+= beta_lpdf(r0 | 2, 2);
	target+= exponential_lpdf(vmr_minus_one | 1);
	target+= beta_lpdf(rho | sampling_k+0.5, sampling_n-sampling_k+0.5); 
	// likelihood
	target+= cs_obs * cs_obs_lpmf[ 2:(N_cs_obs+1) ];
}

generated quantities{
}
"	
	
	stan.code.samplefromprior <- "
data{
	int<lower=1> N_cs_obs;							// max obs chain size
	int<lower=1> N_cs_actual;						// max size of actual chain
	row_vector<lower=0>[N_cs_obs] cs_obs;			// index i holds number of observed chains of size i
	real<lower=0> sampling_n;						// sampling total
	real<lower=0,upper=sampling_n> sampling_k;		// sampling success	
}

parameters{
	real<lower=1e-10, upper=1> r0;	// R0
	real<lower=0> vmr_minus_one;	// variance to mean ratio of NegBin offspring distribution minus one, 1+R0/kappa-1= R0/kappa
	real<lower=0, upper=1> rho;		// sampling probability					
}

model{
	// priors
	target+= beta_lpdf(r0 | 2, 2);
	target+= exponential_lpdf(vmr_minus_one | 1);
	target+= beta_lpdf(rho | sampling_k+0.5, sampling_n-sampling_k+0.5); 	
}
"
	stan.model <- stan_model(model_name= 'Blumberg2013', model_code = gsub('\t',' ',stan.code))
	#stan.model.prior <- stan_model(model_name= 'Blumberg2013_prior', model_code = gsub('\t',' ',stan.code.samplefromprior ))
	
	
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
	stan.data$cs_obs <- tmp$N
	stan.data$N_cs_obs <- length(stan.data$cs_obs)
	stan.data$N_cs_actual <- ceiling(stan.data$N_cs_obs / (args$size.inf.pop.sampled/args$size.inf.pop) * args$upper.bound.multiplier)	#	set upper bound for infinite sum approximation 
	stan.data$sampling_n <- args$size.inf.pop				# replace with actual number infected HSX Amsterdam
	stan.data$sampling_k <- args$size.inf.pop.sampled 		# replace with actual number sequenced infected HSX Amsterdam	
	fit <- rstan::sampling(stan.model, data=stan.data, iter=1e3, warmup=5e2, 
			chains=3, control = list(max_treedepth= 15, adapt_delta= 0.999),
			init= list(list(r0=0.5, vmr_minus_one=0.05, rho=0.5), list(r0=0.25, vmr_minus_one=0.05, rho=0.5), list(r0=0.75, vmr_minus_one=0.05, rho=0.5)))
	
	save(fit, file=file.path(outdir, "200528_Blumberg2013_stanfit.rda"))
	#	runs in 2 minutes
	

	#	examine neff and rhat
	fit.target.pars <- c('r0','vmr_minus_one','rho','kappa')
	rstan::monitor(rstan::extract(fit, pars=fit.target.pars, permuted = FALSE, inc_warmup = TRUE))
	#	neff too low, should ideally be 2e3
	
	#	traces	
	color_scheme_set("mix-blue-red")
	p <- rstan::traceplot(fit, pars=fit.target.pars, inc_warmup=TRUE, ncol = 1)
	pdf(file=file.path(outdir, "200528_Blumberg2013_mcmc_traces.pdf"), w=10, h=8)
	print(p)
	dev.off()
	
	#	keep only what we need to avoid memory meltdown
	fit.po <- rstan::extract(fit, permuted=TRUE, inc_warmup=FALSE)	
	for(x in c('log_rho','bin_lpmf','log_one_minus_rho','log_r0_div_kappa','log_vmr','tmp','ones'))
		fit.po[[x]] <- NULL
	gc()
	
	fit.plot <- matrix(NA, nrow= length(fit.po[[1]]), ncol= length(fit.target.pars))
	colnames(fit.plot) <- fit.target.pars
	for(x in fit.target.pars)
	{
		fit.plot[,x] <- fit.po[[x]]
	}
	
	
	#	pair plots	
	p <- mcmc_pairs(rstan::extract(fit, pars=fit.target.pars, permuted=FALSE, inc_warmup=FALSE), diag_fun = "dens", off_diag_fun = "hex")
	pdf(file=file.path(outdir, "200528_Blumberg2013_mcmc_pairs.pdf"), w=10, h=10)
	print(p)
	dev.off()
	
	# generate posterior predictive samples of actual chain sizes until we reach sampling_n individuals
	cs_actual_cdf <- exp(fit.po$cs_actual_lpmf)
	cs_actual_cdf <- t(apply( cs_actual_cdf, 1, cumsum))	
	samples.n <- 5e3
	cs_actual_postpred <- lapply( seq_len(nrow(cs_actual_cdf)), function(i){
				tmp <- runif(samples.n)				
				cs_actual_postpred <- sapply(tmp, function(x){  head( which( x < cs_actual_cdf[i,] ), 1 )  })
				tmp2 <- which( args$size.inf.pop <= cumsum(cs_actual_postpred) )[1] 
				cs_actual_postpred[1:tmp2]				
			})
	
	
	# posterior estimate of transmissions originating from outside Amsterdam
	sources_external <- sapply(cs_actual_postpred, function(x) length(x)/sum(x) )
	quantile(sources_external, p=c(0.5, 0.025, 0.975)) 
	#	50% 	  2.5%      97.5% 
	#	0.2601837 0.2322827 0.2869928 
}

roadmap.200203.phyloscanner.get.dated.subgraphs<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	require(reshape)
	require(phyloscannerR)
	
	# working directory with phyloscanner output
	home <- '/rds/general/project/ratmann_roadmap_data_analysis/live/analysis_200407'
	home <- '~/Box/OR_Work/AIDSFonds/analysis_200407'
	indir.phsc  <- file.path(home,'phyloscanner_dated')
	outdir <- file.path(home,'subgraphs_dated')
	infiles  <- data.table(F=list.files(indir.phsc, pattern='_annotated_dated_tree.rda$', full.names=TRUE, recursive=TRUE))
	infiles[, SELECT:= gsub('^.*_rerooted_([A-Za-z0-9]+)_.*$','\\1',basename(F))]
	infiles <- subset(infiles, grepl('subtype_B',F))
	
	# extract subgraphs of dated trees
	for(i in seq_len(nrow(infiles)))
	{
		#i <- 3
		cat('process', i,'\n')
		infile <- infiles[i, F]
		host <- infiles[i,SELECT]
		load(infile)
		mrcas <- which( attr(ph, 'SUBGRAPH_MRCA') )
		# some of the tips states are "unknown" and this clashes internally with the NA state, so we need to take some extra care
		# this is not a problem because the "unknown" and NA state mean the same thing
		attr(ph, 'INDIVIDUAL') <- as.character(attr(ph, 'INDIVIDUAL'))
		attr(ph, 'INDIVIDUAL')[is.na(attr(ph, 'INDIVIDUAL'))] <- 'Unknown'
		mrcas <- mrcas[ attr(ph, 'INDIVIDUAL')[mrcas]==host ]
		
		stopifnot( !any(is.na(mrcas)) )
		# check that tree is of class simmap
		stopifnot( any(attr(ph,'class')=='simmap') )
		
		# extract subgraphs
		subgraphs <- lapply(mrcas, function(mrca) extract.subgraph(ph, mrca))
		# save
		outfile <- gsub('_annotated_dated_tree',paste0('_datedsubgraphs_',host),basename(infile))
		outfile <- file.path(outdir,outfile)
		save(subgraphs, file=outfile)
	} 
}


Amsterdam.200317.phyloscanner.B <- function()
{
  require(dplyr)
  require(data.table)
  require(ape)
  require(adephylo)
  require(phytools)
  require(phangorn)
  
  
  plot.phylogenies <- 1
  max.Ntip <- 5e3 
  if(1)
  {
    prog.phyloscanner_analyse_trees <- '/Users/alexb/Documents/software/phyloscanner/phyloscanner/phyloscanner_analyse_trees.R'
    home <- '/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223'    
  }
  if(0)
  {
    prog.phyloscanner_analyse_trees <- '/rdsgpfs/general/user/ablenkin/home/phyloscanner_analyse_trees.R'
    home <- '/rdsgpfs/general/project/ratmann_roadmap_data_analysis/live'   
  }
  
  
  indir.trees <- file.path(home,'fasttree')
  infile.labels <- data.table(FLABEL= c(  file.path(home,'misc','200319_sequence_labels_Ams.csv'),
                                          file.path(home,'misc','200319_sequence_labels_AmsMSM.csv'),
                                          file.path(home,'misc','200319_sequence_labels_AmsHSX.csv'))
  )
  outdir <- file.path(home,'phyloscanner')
  infiles <- data.table(FIN=list.files(indir.trees, pattern='\\.newick$',full.names=TRUE))
  infiles[, DUMMY:= 1L]
  infile.labels[, DUMMY:= 1L]
  infiles <- merge(infiles, infile.labels, by='DUMMY', allow.cartesian=TRUE)
  infiles[, DUMMY:= NULL]
  infiles[, ST:= gsub('.*_subtype_([A-Za-z0-9]+)_.*','\\1',basename(FIN))]    
  infiles <- subset(infiles, ST=='B')
  set.seed(42L)
  for(i in seq_len(nrow(infiles)))
  {
    #i<- 1
    cat('\nprocess ',infiles[i,FIN],'\nprocess ',infiles[i,FLABEL])
    infile <- infiles[i,FIN]
    ph <- read.tree(infile)
    ph$node.label <- NULL
    #   update tip labels: add world region to start of label
    local.world.region <- gsub('^.*_labels_([A-Za-z]+)\\.csv$','\\1', basename(infiles[i,FLABEL]))
    dl <- as.data.table(read.csv(infiles[i,FLABEL], stringsAsFactors=FALSE))
    dl[, X:=NULL]   
    stopifnot(!any(ph$tip.label==''))
    dp <- data.table(IDX=1:Ntip(ph), TAXA=ph$tip.label)
    dp <- merge(dp,dl,by='TAXA',all.x=TRUE)
    dp$TAXA_NEW[grepl('K03455',dp$TAXA)] <- 'HXB2'
    dp$SUBTYPE[grepl('K03455',dp$TAXA)] <- 'B'
    dp$WRLD[grepl('K03455',dp$TAXA)] <- 'Europe'
    dp$SUBTYPE[grepl('Outgroup',dp$TAXA)] <-  gsub('Outgroup.([0-9_A-Z]+).*','\\1',basename(dp$TAXA)[grepl('Outgroup',dp$TAXA)])   
    dp$TAXA_NEW[grepl('Outgroup',dp$TAXA)] <- dp$TAXA[grepl('Outgroup',dp$TAXA)]
    #dp$SUBTYPE[grepl('Outgroup',dp$TAXA)] <- '28_BF'
    stopifnot( nrow(dp[is.na(TAXA_NEW),])==0 )
    stopifnot( !any(duplicated(ph$tip.label)) )
    ph$tip.label <- dp[order(IDX),][, TAXA_NEW]
    #   re-root
    tmp <- dp[, list(NST=length(TAXA)), by='SUBTYPE']
    tmp <- tmp[order(-NST),][2,SUBTYPE]
    tmp <- subset(dp, SUBTYPE==tmp,)[,TAXA_NEW]
    root <- getMRCA(ph, tmp)
    ph <- reroot(ph, root, ph$edge.length[which(ph$edge[,2]==root)]/2)
    #   drop other subtypes
    tmp <- dp[, list(NST=length(TAXA)), by='SUBTYPE']
    # Drop all subtypes which are not the most common
    tmp <- tmp[order(-NST),][-1,SUBTYPE]
    tmp <- subset(dp, SUBTYPE%in%tmp,)[,TAXA_NEW]
    ph <- drop.tip(ph, tmp)
    #   split very large tree into small trees if necessary
    if(is.finite(max.Ntip))
    {
      split.phs <- list()
      repeat({
        local.tips <- which(grepl(paste0('^',local.world.region),ph$tip.label))
        local.tip.ancestors <- Ancestors(ph, sample(local.tips,1))
        cat('\nNumer of local tips left to divide into smaller trees, n=', length(local.tips))
        subtree.mrca <- Ntip(ph)+1L
        for(k in seq_along(local.tip.ancestors))
        {               
          tmp <- length(Descendants(ph, local.tip.ancestors[k], type='tips')[[1]])                        
          if(tmp>max.Ntip)
          {
            subtree.mrca <- local.tip.ancestors[max(1,k-1)]
            break
          }
        }                   
        stopifnot( subtree.mrca>Ntip(ph) )
        tmp <- extract.clade(ph, subtree.mrca)
        if(Ntip(ph)==Ntip(tmp))
        {
          split.phs[[length(split.phs)+1L]] <- tmp
          break
        }
        if(Ntip(ph)-Ntip(tmp)>50)                   
        {                       
          split.phs[[length(split.phs)+1L]] <- tmp
          ph <- drop.tip(ph, split.phs[[length(split.phs)]][['tip.label']])                       
        }
        if(Ntip(ph)-Ntip(tmp)<=50)
        {
          cat('\nFound very small final tree, retry sampling a local.tip')    
        }                   
      })          
    }
    if(!is.finite(max.Ntip))
      split.phs[[1]] <- ph
    #   write to files
    if(length(split.phs)>1)
    {
      for(k in seq_along(split.phs))
      {
        intree.phsc <- gsub('_(subtype_[A-Za-z0-9]+)_',paste0('_\\1c',k,'_'),gsub('\\.newick',paste0('_rerooted_',local.world.region,'.newick'),basename(infile)))
        intree.phsc <- file.path(outdir,intree.phsc)
        stopifnot(is.binary(split.phs[[k]]))
        write.tree(split.phs[[k]], file=intree.phsc)
        if(plot.phylogenies)
        {
          pdf(file=gsub('newick','pdf',intree.phsc), w=20, h=10+Ntip(split.phs[[k]])/10)
          plot(split.phs[[k]], show.node.label=TRUE, cex=0.3)
          dev.off()           
        }
      }                           
    }
    if(length(split.phs)==1)
    {
      intree.phsc <- file.path(outdir,gsub('\\.newick',paste0('_rerooted_',local.world.region,'.newick'),basename(infile)))
      write.tree(split.phs[[1]], file=intree.phsc)
      if(plot.phylogenies)
      {
        pdf(file=gsub('newick','pdf',intree.phsc), w=20, h=10+Ntip(split.phs[[1]])/10)
        plot(split.phs[[1]], show.node.label=TRUE, cex=0.3)
        dev.off()           
      }
    }
  }
  
  infiles <- data.table(FIN=list.files(outdir, pattern='\\.newick$',full.names=TRUE))
  infiles[, ST:= gsub('^.*_Subtype([A-Za-z0-9]+)_.*$','\\1',basename(FIN))]
  infiles[, DUMMY:= gsub('\\.newick$','',basename(FIN))]
  tmp <- data.table(FOUT=list.files(outdir, pattern='workspace.rda$',full.names=TRUE))
  tmp[, DUMMY:= gsub('__workspace.rda$','',basename(FOUT))]
  infiles <- merge(infiles, tmp, by='DUMMY', all.x=TRUE)
  infiles <- subset(infiles, is.na(FOUT), c(ST, FIN))
  infiles <- subset(infiles, grepl('Bc[0-9]+',ST))
  cmds <- vector('list',nrow(infiles))    
  for(i in seq_len(nrow(infiles)))
  {               
    #   make phyloscanner UNIX command
    infile <- infiles[i,FIN]
    outputString <- paste0(gsub('\\.newick','_',infile))
    tip.regex <- "^([A-Za-z]+)___.*$"
    cmd <- paste("CWD=$(pwd)\n",sep='')
    cmd <- paste(cmd,"echo $CWD\n",sep='')  
    tmpdir.prefix <- paste('phsc_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
    tmpdir <- paste("$CWD/",tmpdir.prefix,sep='')
    tmp.in <- file.path(tmpdir, basename(infile))
    cmd <- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
    cmd <- paste(cmd,'cp "',infile,'" ',tmp.in,'\n', sep='')
    cmd <- paste(cmd,'cd ', tmpdir,'\n', sep='')
    cmd <- paste0(cmd,'Rscript ',prog.phyloscanner_analyse_trees,' ',basename(infile),' ',basename(outputString))
    #   don t use -m multifurcation threshold to ensure that output tree is still binary
    cmd <- paste0(cmd,' s,0 -x "',tip.regex,'" -v 1 -ow -rda\n')
    cmd <- paste0(cmd,'mv ',basename(outputString),'* ','"',dirname(outputString),'"','\n') 
    cmd <- paste(cmd, "cd $CWD\n",sep='')
    cmd <- paste(cmd, "rm ", tmpdir,'\n',sep='')
    cmds[[i]] <- cmd
  }
  infiles[, CMD:= unlist(cmds)]   
  
  #   submit jobs like this one:
  cat(infiles[1,CMD])
  
  #   run on HPC as array job
  df <- copy(infiles)
  df[, CASE_ID:= 1:nrow(df)]  
  #   make PBS header
  hpc.load    <- "module load anaconda3/personal"
  hpc.select  <- 1                        # number of nodes
  hpc.nproc   <- 1                        # number of processors on node
  hpc.walltime<- 23                       # walltime
  hpc.q       <- NA #"pqeelab"                        # PBS queue
  hpc.mem     <- "4gb"                    # RAM
  hpc.array   <- length(unique(df$CASE_ID))   # number of runs for job array
  pbshead     <- "#!/bin/sh"
  tmp         <- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  tmp         <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  pbshead     <- paste(pbshead, "#PBS -j oe", sep = "\n") 
  if(!is.na(hpc.array))
    pbshead <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')        
  if(!is.na(hpc.q)) 
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead     <- paste(pbshead, hpc.load, sep = "\n")         
  #   make array job
  cmd     <- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
  cmd     <- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]         
  cmd     <- paste(pbshead,cmd,sep='\n')  
  #   submit job
  outfile     <- gsub(':','',paste("phs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
  outfile     <- file.path(outdir, outfile)
  cat(cmd, file=outfile)
  cmd         <- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))  
}


amsterdam.200306.alignment.make.bootstraps()
amsterdam.200306.fastree()
