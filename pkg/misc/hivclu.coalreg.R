cr.various<- function()
{
	#cr.various.master()
	cr.various.pangea()
}

cr.hpc.submit<- function()
{
	if(1)
	{
		par.maxNodeDepth	<- Inf
		par.maxHeight		<- 10	
		par.lasso			<- 5
		par.bias			<- 1
		par.scale			<- 0
		par.climb			<- 'BFGS'		
		par.noise			<- 0
		for(i in 44:51)
		{
			par.base.pattern	<- paste0('PANGEA-AcuteHigh-InterventionNone-cov11.8-seed',i)
			cmd					<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -par.base.pattern=',par.base.pattern,' -par.maxNodeDepth=',par.maxNodeDepth,' -par.maxHeight=',par.maxHeight, ' -par.lasso=',par.lasso, ' -par.noise=', par.noise, ' -par.bias=', par.bias, ' -par.climb=', par.climb, ' -par.scale=', par.scale)			
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=171, hpc.mem="3600mb")
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)		
		}
		quit("no")	 
	}
	if(0)
	{
		par.maxNodeDepth	<- Inf
		par.maxHeight		<- 10	
		par.lasso			<- 5
		par.scale			<- 0
		par.climb			<- 'BFGS'
		for(par.bias in c(0.5,1,2))
			for(par.noise in c(0.5,1,2))
				for(i in 44:51)
				{
					par.base.pattern	<- paste0('PANGEA-AcuteHigh-InterventionNone-cov11.8-seed',i)
					cmd					<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -par.base.pattern=',par.base.pattern,' -par.maxNodeDepth=',par.maxNodeDepth,' -par.maxHeight=',par.maxHeight, ' -par.lasso=',par.lasso, ' -par.noise=', par.noise, ' -par.bias=', par.bias, ' -par.climb=', par.climb, ' -par.scale=', par.scale)			
					cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeph', hpc.walltime=171, hpc.mem="3600mb")
					cat(cmd)	
					outdir		<- paste(DATA,"tmp",sep='/')
					outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
					hivc.cmd.hpccaller(outdir, outfile, cmd)		
				}
		quit("no")	
	}
	if(0)
	{
		par.maxNodeDepth	<- Inf
		par.maxHeight		<- 10	
		par.lasso			<- 5
		par.noise			<- 0
		for(par.bias in c(0.5,1,2))
			for(par.climb in c('BFGS','Nelder-Mead'))
				for(par.scale in c(0,1))
					for(i in 44:51)
					{
						par.base.pattern	<- paste0('PANGEA-AcuteHigh-InterventionNone-cov11.8-seed',i)
						cmd					<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -par.base.pattern=',par.base.pattern,' -par.maxNodeDepth=',par.maxNodeDepth,' -par.maxHeight=',par.maxHeight, ' -par.lasso=',par.lasso, ' -par.noise=', par.noise, ' -par.bias=', par.bias, ' -par.climb=', par.climb, ' -par.scale=', par.scale)			
						cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=171, hpc.mem="3600mb")
						cat(cmd)	
						outdir		<- paste(DATA,"tmp",sep='/')
						outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
						hivc.cmd.hpccaller(outdir, outfile, cmd)		
					}	
		quit("no")	
	}
}

cr.various.master<- function()
{
	
	indir				<- '/work/or105/ATHENA_2016/master_examples'
	par.base.pattern	<- 'm3.RR5.n1250_seed123'	
	if(0)
	{
		par.s				<- 0.5
		cr.master.ex3.runcoalreg.using.TYPE.BFGS2(indir, par.base.pattern, par.s)	
	}	
	if(0)
	{
		par.s				<- 0.5
		cr.master.ex3.runcoalreg.using.TYPE.ETFI.vanilla.BFGS2(indir, par.base.pattern, par.s)	
	}
	if(1)
	{
		par.s				<- 0.5
		par.maxNodeDepth	<- 15
		par.maxHeight		<- Inf
		cr.master.ex3.runcoalreg.using.TYPE.ETFI.vanilla.BFGS3(indir, par.base.pattern, par.s, par.maxNodeDepth, par.maxHeight=par.maxHeight)	
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

cr.various.pangea<- function()
{
	indir				<- '/work/or105/coalreg/pangea_examples'
	par.base.pattern	<- 'PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43'
	par.maxHeight		<- Inf
	par.maxNodeDepth	<- Inf
	par.climb			<- 'BFGS'
	par.lasso			<- 5
	par.noise			<- 0
	par.bias			<- 1
	par.scale			<- 0
	#
	#	read args
	#
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,17), par.base.pattern= return(substr(arg,19,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.base.pattern<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,17), par.maxNodeDepth= return(substr(arg,19,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.maxNodeDepth<- as.numeric(tmp[1])		
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,14), par.maxHeight= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.maxHeight<- as.numeric(tmp[1])
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,10), par.lasso= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.lasso<- as.numeric(tmp[1])		
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,10), par.noise= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.noise<- as.numeric(tmp[1])
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,9), par.bias= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.bias<- as.numeric(tmp[1])		
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,10), par.scale= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.scale<- as.numeric(tmp[1])
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,10), par.climb= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.climb<- tmp[1]
	}	
	cat('input args\n',par.base.pattern,'\n',par.maxNodeDepth,'\n',par.maxHeight,'\n', par.noise,'\n', par.bias, '\n', par.climb, '\n', par.scale)	
	if(0)
	{		
		cr.png.runcoalreg.using.TRSTAGE.ETSI.vanilla.BFGS3(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso)		
	}	
	if(0)
	{
		cr.png.runcoalreg.using.TRSTAGE.TRRISK.ETSI.vanilla.BFGS3(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso)
	}
	if(0)
	{
		cr.png.runcoalreg.using.TRGENDER.ETSI.vanilla.BFGS3(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso)
	}
	if(0)
	{
		cr.png.runcoalreg.using.TRRISK.ETSI.vanilla.BFGS3(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso)
	}	
	if(0)
	{
		cr.png.runcoalreg.using.TRSTAGE.TRRISK.TRGENDER.ETSI.vanilla.BFGS3(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso)
	}
	if(0)
	{
		cr.png.runcoalreg.using.TRSTAGE.TRRISK.TRGENDER.TRAGEDIAG.ETSI.vanilla.BFGS3(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso)
	}
	if(0)
	{
		cr.png.runcoalreg.using.TRSTAGE.TRRISK.ETSI.noise.BFGS3(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale)
	}
	if(0)
	{
		cr.png.runcoalreg.using.TRSTAGE.TRRISK.TRGENDER.TRAGEDIAG.ETSI.noise.BFGS3(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale)
	}
	if(0)
	{
		cr.png.runcoalreg.using.TRSTAGE.ETSI.unmodelled.het.BFGS3(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale)
	}
	if(0)
	{
		cr.png.runcoalreg.using.TRRISK.ETSI.unmodelled.het.BFGS3(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale)
	}
	if(1)
	{
		#cr.png.runcoalreg.using.TRvariable.ETSInoise.unmodelled.het(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_IN_RECENT'))
		cr.png.runcoalreg.using.TRvariable.ETSInoise.aoi.model1(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_IN_RECENT'))
	}
	if(0)
	{
		#cr.png.runcoalreg.using.TRvariable.ETSInoise.unmodelled.het(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_IN_SIXM'))
		cr.png.runcoalreg.using.TRvariable.ETSInoise.aoi.model1(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_IN_SIXM'))
	}
	if(0)
	{
		#cr.png.runcoalreg.using.TRvariable.ETSInoise.unmodelled.het(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_IN_THREEM'))
		cr.png.runcoalreg.using.TRvariable.ETSInoise.aoi.model1(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_IN_THREEM'))
	}
	if(0)
	{
		#cr.png.runcoalreg.using.TRvariable.ETSInoise.unmodelled.het(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_CD4_STAGE2'))
		cr.png.runcoalreg.using.TRvariable.ETSInoise.aoi.model1(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_CD4_STAGE2'))
	}	
	if(0)
	{
		#cr.png.runcoalreg.using.TRvariable.ETSInoise.unmodelled.het(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_IN_RECENT','MALE','AGE_AT_DIAG'))
		cr.png.runcoalreg.using.TRvariable.ETSInoise.aoi.model1(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_IN_RECENT','MALE','AGE_AT_DIAG'))
	}
	if(0)
	{
		#cr.png.runcoalreg.using.TRvariable.ETSInoise.unmodelled.het(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_IN_SIXM','MALE','AGE_AT_DIAG'))
		cr.png.runcoalreg.using.TRvariable.ETSInoise.aoi.model1(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_IN_SIXM','MALE','AGE_AT_DIAG'))
	}
	if(0)
	{
		#cr.png.runcoalreg.using.TRvariable.ETSInoise.unmodelled.het(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_IN_THREEM','MALE','AGE_AT_DIAG'))
		cr.png.runcoalreg.using.TRvariable.ETSInoise.aoi.model1(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_IN_THREEM','MALE','AGE_AT_DIAG'))
	}
	if(0)
	{
		#cr.png.runcoalreg.using.TRvariable.ETSInoise.unmodelled.het(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_CD4_STAGE2','MALE','AGE_AT_DIAG'))
		cr.png.runcoalreg.using.TRvariable.ETSInoise.aoi.model1(indir, par.base.pattern, par.maxNodeDepth=par.maxNodeDepth, par.maxHeight=par.maxHeight, par.lasso=par.lasso, par.noise=par.noise, par.bias=par.bias, par.climb=par.climb, par.scale=par.scale, trm.factors= c('DIAG_CD4_STAGE2','MALE','AGE_AT_DIAG'))		
	}
	
}

cr.png.generate.data<- function()
{
	require(data.table)
	require(ape)
	require(phangorn)
	indir	<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations'
	infiles	<- data.table(FR=list.files(indir, pattern='_INTERNAL.R', full.names=TRUE, recursive=TRUE))
	infiles[, SC:= basename(dirname(dirname(FR)))]
	infiles[, SEQCOV:= as.numeric(gsub('.*cov([0-9]+\\.[0-9]+).*','\\1',SC))/100]
	infiles[, SEED:= as.numeric(gsub('.*seed([0-9]+).*','\\1',SC))]	
	tmp		<- data.table(FT=list.files(indir, pattern='_DATEDTREE.newick', full.names=TRUE, recursive=TRUE))
	tmp[, SC:= basename(dirname(dirname(FT)))]
	infiles	<- merge(infiles, tmp, by='SC', all=1)
	
	infiles[, table(SEQCOV)]
	
	outfile	<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-'
	t.start	<- 0
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
	
	FR		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed50/150129_HPTN071_scHN_INTERNAL/150129_HPTN071_scHN_SIMULATED_INTERNAL.R'
	#FR		<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov41.3-seed46/150129_HPTN071_scHN_INTERNAL/150129_HPTN071_scHN_SIMULATED_INTERNAL.R'
	#FR		<- '~/Downloads/150129_HPTN071_scHN_SIMULATED_INTERNAL.R'
	cat('\n',FR)
	load(FR)
	dfo		<- data.table(TS=c(0, 2000, 2005, 2005, 2010), TE=c(2020, 2010, 2015,2020,2020))
	dfr		<- dfo[, {
				t.start	<- TS
				t.end	<- TE
				#
				#	denominator population: all adults by 2020 that were infected before 2020
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & TIME_TR>t.start & HIV=='Y' &  TIME_TR<t.end)
				#
				#	transmission rate given risk category: 
				#	what I can calculate is the #transmissions per year when in stage X
				#	so I need the duration that individuals spend in stage X and the #transmissions while in stage X
				#	this is easy for risk group, because individuals stay in the risk group forever
				#
				#	calculate: time between infection until death or end of observation 2020
				dfd[, DUR_RISK:= DOD]
				#tmp		<- dfd[, which(!is.na(ART1_T))]
				#set(dfd, tmp, 'DUR_RISK', dfd[tmp, ART1_T])				
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
				set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-pmax(t.start,TIME_TR)])				
				#	calculate: number of infections until 2020
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by='IDTR']
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				dfr		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='RISK'], measure.vars='RISK', variable.name='GROUP', value.name='FACTOR')
				
				#	for individuals that started ART or did not start ART
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR<t.end & TIME_TR>t.start)
				dfd[, EVER_ART:= factor(is.na(ART1_T), levels=c(TRUE,FALSE), labels=c('N','Y'))]
				dfd[, DUR_RISK:= DOD]
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
				set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-pmax(t.start,TIME_TR)])
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by='IDTR']
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				tmp		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='EVER_ART'], measure.vars='EVER_ART', variable.name='GROUP', value.name='FACTOR') 
				dfr		<- rbind(dfr, tmp)
				#tmp		<- melt(subset(dfd, RISK=='M')[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='ART_STARTED'], measure.vars='ART_STARTED', variable.name='GROUP', value.name='FACTOR')
				#tmp[, GROUP:='ART_STARTED_RISKM']
				#dfr		<- rbind(dfr, tmp)
				#tmp		<- melt(subset(dfd, RISK=='M' & GENDER=='F')[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='ART_STARTED'], measure.vars='ART_STARTED', variable.name='GROUP', value.name='FACTOR')
				#tmp[, GROUP:='ART_STARTED_RISKM_GENDERF']
				#dfr		<- rbind(dfr, tmp)
				
				
				#	for males
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR<t.end & TIME_TR>t.start)
				dfd	<- subset(dfd, is.na(ART1_T))
				dfd[, DUR_RISK:= DOD]
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
				set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-pmax(t.start,TIME_TR)])
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by='IDTR']
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				tmp		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='GENDER'], measure.vars='GENDER', variable.name='GROUP', value.name='FACTOR') 
				dfr		<- rbind(dfr, tmp)
				#	conditional trm risk of males for RISK group= M
				#tmp		<- melt(subset(dfd, DUR_RISK>0 & RISK=='M')[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='GENDER'], measure.vars='GENDER', variable.name='GROUP', value.name='FACTOR')
				#tmp[, GROUP:='GENDER_RISKM']
				#dfr		<- rbind(dfr, tmp)	
				
				#	for period before / after ART
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR<t.end & TIME_TR>t.start & !is.na(DIAG_T))						
				tmp		<- as.data.table(expand.grid(ART_STAGE=c('before start','after start'),IDPOP=unique(dfd$IDPOP)))
				dfd		<- merge(dfd,tmp,by='IDPOP')
				dfd[, DUR_RISK:= DOD]
				dfd[, DUMMY:= dfd[, pmax(t.start,TIME_TR)]]
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
				tmp		<- dfd[, which(ART_STAGE=='before start' & !is.na(ART1_T) & ART1_T<t.end)]
				set(dfd, tmp, 'DUR_RISK', dfd[tmp, ART1_T])
				tmp		<- dfd[, which(ART_STAGE=='after start' & !is.na(ART1_T) & ART1_T<t.end)]
				set(dfd, tmp, 'DUMMY', dfd[tmp, ART1_T])
				tmp		<- dfd[, which(ART_STAGE=='after start' & (is.na(ART1_T) | !is.na(ART1_T) & ART1_T>=t.end))]
				set(dfd, tmp, 'DUMMY', dfd[tmp, DUR_RISK])
				set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-DUMMY])
				tmp		<- copy(df.trms)
				setnames(tmp, 'IDTR', 'IDPOP')
				tmp		<- merge(tmp, subset(dfd, select=c(IDPOP,ART1_T)), by='IDPOP')				
				tmp[, ART_STAGE:= as.character(factor(TIME_TR<=ART1_T,levels=c(TRUE,FALSE),labels=c('before start','after start')))]
				set(tmp, tmp[, which(is.na(ART_STAGE))], 'ART_STAGE', 'before start')				
				tmp		<- subset(tmp, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by=c('IDPOP','ART_STAGE')]
				dfd		<- merge(dfd, tmp, by=c('IDPOP','ART_STAGE'), all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				#	marginal trm 
				tmp		<- melt(subset(dfd, DUR_RISK>0)[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='ART_STAGE'], measure.vars='ART_STAGE', variable.name='GROUP', value.name='FACTOR')
				dfr		<- rbind(dfr, tmp)
				#	conditional trm risk for RISK group= M
				#tmp		<- melt(subset(dfd, DUR_RISK>0 & RISK=='M')[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='ART_STAGE'], measure.vars='ART_STAGE', variable.name='GROUP', value.name='FACTOR')
				#tmp[, GROUP:='ART_STAGE_RISKM']
				#dfr		<- rbind(dfr, tmp)
				#	conditional trm risk for RISK group= M, GENDER=F
				#tmp		<- melt(subset(dfd, DUR_RISK>0 & RISK=='M' & GENDER=='F')[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='ART_STAGE'], measure.vars='ART_STAGE', variable.name='GROUP', value.name='FACTOR')
				#tmp[, GROUP:='ART_STAGE_RISKM_GENDERF']
				#dfr		<- rbind(dfr, tmp)
				
				
				#	for ART start-DIAG_T < 12m or >12m
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR<t.end & TIME_TR>t.start & !is.na(DIAG_T))
				dfd[, TIME_TO_ART:= cut(ART1_T-DIAG_T, breaks=c(-1,.5,1e4),labels=c('<6m to ART','>6m to ART'))]
				set(dfd, dfd[, which(is.na(TIME_TO_ART))], 'TIME_TO_ART', 'no ART')  				
				dfd[, DUR_RISK:= DOD]
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
				set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-pmax(t.start,TIME_TR)])
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by='IDTR']
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				tmp		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='TIME_TO_ART'], measure.vars='TIME_TO_ART', variable.name='GROUP', value.name='FACTOR') 
				dfr		<- rbind(dfr, tmp)
				#tmp		<- melt(subset(dfd, RISK=='M')[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='TIME_TO_ART'], measure.vars='TIME_TO_ART', variable.name='GROUP', value.name='FACTOR')
				#tmp[, GROUP:='TIME_TO_ART_RISKM']
				#dfr		<- rbind(dfr, tmp)
				#tmp		<- melt(subset(dfd, RISK=='M' & GENDER=='F')[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='TIME_TO_ART'], measure.vars='TIME_TO_ART', variable.name='GROUP', value.name='FACTOR')
				#tmp[, GROUP:='TIME_TO_ART_RISKM_GENDERF']
				#dfr		<- rbind(dfr, tmp)
				
				
				
				#	for transmission during early infection < 1 year vs late infection > 1 year
				#	(of course this will never be known, but include as ideal case)
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR<t.end & DOD>t.start)
				tmp		<- as.data.table(expand.grid(TR_STAGE=c('yes','no'),IDPOP=unique(dfd$IDPOP)))
				dfd		<- merge(dfd,tmp,by='IDPOP')			
				dfd[, DUR_RISK:= DOD]
				tmp		<- 12/12
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)				 
				set(dfd, dfd[, which(TR_STAGE=='yes' & (TIME_TR+tmp)<t.start)], 'DUR_RISK', 0)	# early stage not in observation period
				set(dfd, dfd[, which(TR_STAGE=='yes' & (TIME_TR+tmp)>=t.start)], 'DUR_RISK', dfd[TR_STAGE=='yes' & (TIME_TR+tmp)>=t.start, pmin(tmp, pmin(TIME_TR+tmp,DUR_RISK)-t.start)]) # some of the early stage in observation period				
				set(dfd, dfd[, which(TR_STAGE=='no' & (TIME_TR+tmp)<t.start)], 'DUR_RISK', dfd[TR_STAGE=='no' & (TIME_TR+tmp)<t.start, DUR_RISK-pmax(t.start,TIME_TR)])	# only late stage in observation period
				set(dfd, dfd[, which(TR_STAGE=='no' & (TIME_TR+tmp)>=t.start)], 'DUR_RISK', dfd[TR_STAGE=='no' & (TIME_TR+tmp)>=t.start, (DUR_RISK-pmax(t.start,TIME_TR)) - (pmin(tmp, pmin(TIME_TR+tmp,DUR_RISK)-t.start))]) # some of the early stage in observation period
				df.trms[, TR_STAGE:= as.character(factor(TIME_TR-IDTR_TIME_INFECTED<=tmp,levels=c(TRUE,FALSE),labels=c('yes','no')))]
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by=c('IDTR','TR_STAGE')]
				df.trms[, TR_STAGE:=NULL]
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by=c('IDPOP','TR_STAGE'), all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				#	marginal trm 
				tmp		<- melt(subset(dfd, DUR_RISK>0)[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='TR_STAGE'], measure.vars='TR_STAGE', variable.name='GROUP', value.name='FACTOR')
				tmp[, GROUP:='TRM_IN_TWM']
				dfr		<- rbind(dfr, tmp)
				#	conditional trm risk for RISK group= M
				#tmp		<- melt(subset(dfd, DUR_RISK>0 & RISK=='M')[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='TR_STAGE'], measure.vars='TR_STAGE', variable.name='GROUP', value.name='FACTOR')
				#tmp[, GROUP:='TR_STAGE_RISKM']
				#dfr		<- rbind(dfr, tmp)
				#	conditional trm risk for RISK group= M, GENDER=F
				#tmp		<- melt(subset(dfd, DUR_RISK>0 & RISK=='M' & GENDER=='F')[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='TR_STAGE'], measure.vars='TR_STAGE', variable.name='GROUP', value.name='FACTOR')
				#tmp[, GROUP:='TR_STAGE_RISKM_GENDERF']
				#dfr		<- rbind(dfr, tmp)
				
				#	for transmission during first 6 months vs after 6 months of infection
				#	(of course this will never be known, but include as ideal case)
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR<t.end & TIME_TR>t.start)
				tmp		<- as.data.table(expand.grid(TR_STAGE=c('yes','no'),IDPOP=unique(dfd$IDPOP)))
				dfd		<- merge(dfd,tmp,by='IDPOP')			
				dfd[, DUR_RISK:= DOD]
				tmp		<- 6/12
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)				 
				set(dfd, dfd[, which(TR_STAGE=='yes' & (TIME_TR+tmp)<t.start)], 'DUR_RISK', 0)	# early stage not in observation period
				set(dfd, dfd[, which(TR_STAGE=='yes' & (TIME_TR+tmp)>=t.start)], 'DUR_RISK', dfd[TR_STAGE=='yes' & (TIME_TR+tmp)>=t.start, pmin(tmp, pmin(TIME_TR+tmp,DUR_RISK)-t.start)]) # some of the early stage in observation period				
				set(dfd, dfd[, which(TR_STAGE=='no' & (TIME_TR+tmp)<t.start)], 'DUR_RISK', dfd[TR_STAGE=='no' & (TIME_TR+tmp)<t.start, DUR_RISK-pmax(t.start,TIME_TR)])	# only late stage in observation period
				set(dfd, dfd[, which(TR_STAGE=='no' & (TIME_TR+tmp)>=t.start)], 'DUR_RISK', dfd[TR_STAGE=='no' & (TIME_TR+tmp)>=t.start, (DUR_RISK-pmax(t.start,TIME_TR)) - (pmin(tmp, pmin(TIME_TR+tmp,DUR_RISK)-t.start))]) # some of the early stage in observation period
				df.trms[, TR_STAGE:= as.character(factor(TIME_TR-IDTR_TIME_INFECTED<=tmp,levels=c(TRUE,FALSE),labels=c('yes','no')))]
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by=c('IDTR','TR_STAGE')]
				df.trms[, TR_STAGE:=NULL]
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by=c('IDPOP','TR_STAGE'), all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				#	marginal trm 
				tmp		<- melt(subset(dfd, DUR_RISK>0)[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='TR_STAGE'], measure.vars='TR_STAGE', variable.name='GROUP', value.name='FACTOR')
				tmp[, GROUP:='TRM_IN_SIXM']
				dfr		<- rbind(dfr, tmp)
				
				
				#	for transmission during first 3 months vs after 3 months of infection
				#	(of course this will never be known, but include as ideal case)
				dfd		<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR<t.end & TIME_TR>t.start)
				tmp		<- as.data.table(expand.grid(TR_STAGE=c('yes','no'),IDPOP=unique(dfd$IDPOP)))
				dfd		<- merge(dfd,tmp,by='IDPOP')			
				dfd[, DUR_RISK:= DOD]
				tmp		<- 3/12
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)				 
				set(dfd, dfd[, which(TR_STAGE=='yes' & (TIME_TR+tmp)<t.start)], 'DUR_RISK', 0)	# early stage not in observation period
				set(dfd, dfd[, which(TR_STAGE=='yes' & (TIME_TR+tmp)>=t.start)], 'DUR_RISK', dfd[TR_STAGE=='yes' & (TIME_TR+tmp)>=t.start, pmin(tmp, pmin(TIME_TR+tmp,DUR_RISK)-t.start)]) # some of the early stage in observation period				
				set(dfd, dfd[, which(TR_STAGE=='no' & (TIME_TR+tmp)<t.start)], 'DUR_RISK', dfd[TR_STAGE=='no' & (TIME_TR+tmp)<t.start, DUR_RISK-pmax(t.start,TIME_TR)])	# only late stage in observation period
				set(dfd, dfd[, which(TR_STAGE=='no' & (TIME_TR+tmp)>=t.start)], 'DUR_RISK', dfd[TR_STAGE=='no' & (TIME_TR+tmp)>=t.start, (DUR_RISK-pmax(t.start,TIME_TR)) - (pmin(tmp, pmin(TIME_TR+tmp,DUR_RISK)-t.start))]) # some of the early stage in observation period
				df.trms[, TR_STAGE:= as.character(factor(TIME_TR-IDTR_TIME_INFECTED<=tmp,levels=c(TRUE,FALSE),labels=c('yes','no')))]
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by=c('IDTR','TR_STAGE')]
				df.trms[, TR_STAGE:=NULL]
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by=c('IDPOP','TR_STAGE'), all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 	# if individuals have early stage before observation period, TR_N is set to 0. This is OK because in the next step we only consider individuals with DUR_RISK>0, and this excludes the individuals in question
				#	marginal trm 
				tmp		<- melt(subset(dfd, DUR_RISK>0)[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='TR_STAGE'], measure.vars='TR_STAGE', variable.name='GROUP', value.name='FACTOR')
				tmp[, GROUP:='TRM_IN_THREEM']
				dfr		<- rbind(dfr, tmp)


				#	for diagnosed w CD4<500 vs diagnosed w CD4>500
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR<t.end & TIME_TR>t.start & !is.na(DIAG_T))
				dfd[, DIAG_CD4_STAGE2:= cut(DIAG_CD4, breaks=c(-1,500,1e4),labels=c('l500','g500'))]			
				dfd[, DUR_RISK:= DOD]
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
				set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-pmax(t.start,TIME_TR)])
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by='IDTR']
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				tmp		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='DIAG_CD4_STAGE2'], measure.vars='DIAG_CD4_STAGE2', variable.name='GROUP', value.name='FACTOR') 
				dfr		<- rbind(dfr, tmp)
				
				
				#	for diagnosed w CD4<350 vs diagnosed w CD4>350
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR<t.end & TIME_TR>t.start & !is.na(DIAG_T))
				dfd[, DIAG_CD4_STAGE:= cut(DIAG_CD4, breaks=c(-1,350,1e4),labels=c('l350','g350'))]			
				dfd[, DUR_RISK:= DOD]
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
				set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-pmax(t.start,TIME_TR)])
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by='IDTR']
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				tmp		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='DIAG_CD4_STAGE'], measure.vars='DIAG_CD4_STAGE', variable.name='GROUP', value.name='FACTOR') 
				dfr		<- rbind(dfr, tmp)
				
				
				
				#	for diagnosed within 12 months of infection vs later
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR>t.start & TIME_TR<t.end & !is.na(DIAG_T))
				dfd[, DIAG_IN_RECENT:= as.character(factor((DIAG_T-TIME_TR)>=1, levels=c(TRUE,FALSE), labels=c('N','Y')))]
				dfd[, DUR_RISK:= DOD]
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
				set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-pmax(t.start,TIME_TR)])
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by='IDTR']
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				tmp		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='DIAG_IN_RECENT'], measure.vars='DIAG_IN_RECENT', variable.name='GROUP', value.name='FACTOR') 
				dfr		<- rbind(dfr, tmp)
				
				
				#	for diagnosed within three months of infection vs later 
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR>t.start & TIME_TR<t.end & !is.na(DIAG_T))
				dfd[, DIAG_IN_THREEM:= as.character(factor((DIAG_T-TIME_TR)>=3/12, levels=c(TRUE,FALSE), labels=c('N','Y')))]
				dfd[, DUR_RISK:= DOD]
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
				set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-pmax(t.start,TIME_TR)])
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by='IDTR']
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				tmp		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='DIAG_IN_THREEM'], measure.vars='DIAG_IN_THREEM', variable.name='GROUP', value.name='FACTOR') 
				dfr		<- rbind(dfr, tmp)				
				#melt(subset(dfd, TIME_TR<2017 & TIME_TR>2016)[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='DIAG_IN_THREEM'], measure.vars='DIAG_IN_THREEM', variable.name='GROUP', value.name='FACTOR')
				#subset(dfd, TIME_TR>2010 & ART1_T<2015)[, table(DIAG_IN_THREEM)]
				#melt(subset(dfd, TIME_TR>2010 & ART1_T<2015)[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='DIAG_IN_THREEM'], measure.vars='DIAG_IN_THREEM', variable.name='GROUP', value.name='FACTOR')
				
				#	for diagnosed within six months of infection vs later 
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR>t.start & TIME_TR<t.end & !is.na(DIAG_T))
				dfd[, DIAG_IN_SIXM:= as.character(factor((DIAG_T-TIME_TR)>=.5, levels=c(TRUE,FALSE), labels=c('N','Y')))]
				dfd[, DUR_RISK:= DOD]
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
				set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-pmax(t.start,TIME_TR)])
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by='IDTR']
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				tmp		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='DIAG_IN_SIXM'], measure.vars='DIAG_IN_SIXM', variable.name='GROUP', value.name='FACTOR') 
				dfr		<- rbind(dfr, tmp)
								
				
				#	age at diagnosis	
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR>t.start & TIME_TR<t.end )
				#dfd	<- subset(dfd, is.na(ART1_T))
				dfd[, AGE_AT_DIAG:= DIAG_T-DOB]
				set(dfd, NULL, 'AGE_AT_DIAG', dfd[, cut(AGE_AT_DIAG, breaks=c(-1,25,30,1e4),labels=c('less25','25to29','30'))])
				set(dfd, dfd[, which(is.na(AGE_AT_DIAG))], 'AGE_AT_DIAG', 'unknown')
				dfd[, DUR_RISK:= DOD]
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
				set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-pmax(t.start,TIME_TR)])
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by='IDTR']
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				tmp		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='AGE_AT_DIAG'], measure.vars='AGE_AT_DIAG', variable.name='GROUP', value.name='FACTOR')
				dfr		<- rbind(dfr, tmp)
				#ggplot(dfd, aes(x= floor(TIME_TR), fill=AGE_AT_DIAG)) + geom_bar(position='fill')
				#ggplot(dfd, aes(x= floor(TIME_TR), fill=cut(TIME_TR-DOB, breaks=c(-1,25,30,1e4)))) + geom_bar(position='fill')
			
				#	age at infection	
				dfd	<- subset(df.inds, IDPOP>0 & DOB<(t.end-t.adult) & HIV=='Y' & TIME_TR>t.start & TIME_TR<t.end )
				dfd	<- subset(dfd, is.na(ART1_T))
				dfd[, AGE_AT_INFECTION:= TIME_TR-DOB]
				set(dfd, NULL, 'AGE_AT_INFECTION', dfd[, cut(AGE_AT_INFECTION, breaks=c(-1,25,30,1e4),labels=c('less25','25to29','30'))])
				set(dfd, dfd[, which(is.na(AGE_AT_INFECTION))], 'AGE_AT_INFECTION', 'unknown')
				dfd[, DUR_RISK:= DOD]
				set(dfd, dfd[, which(is.na(DUR_RISK) | DUR_RISK>t.end)], 'DUR_RISK', t.end)
				set(dfd, NULL, 'DUR_RISK', dfd[, DUR_RISK-pmax(t.start,TIME_TR)])
				tmp		<- subset(df.trms, TIME_TR<t.end & TIME_TR>t.start)[, list(TR_N=length(IDREC)), by='IDTR']
				setnames(tmp, 'IDTR', 'IDPOP')
				dfd		<- merge(dfd, tmp, by='IDPOP', all.x=1)
				set(dfd, dfd[, which(is.na(TR_N))], 'TR_N', 0L) 
				tmp		<- melt(dfd[, list(TR_RATE= mean(TR_N/DUR_RISK)), by='AGE_AT_INFECTION'], measure.vars='AGE_AT_INFECTION', variable.name='GROUP', value.name='FACTOR')				
				dfr		<- rbind(dfr, tmp)									
				#	
				
				dfr				
			}, by=c('TS','TE')]
	#	generate transmission risk ratios
	dtrr	<- dfr[, {
				tmp		<- c(	'b_TRM_IN_TWM_yes'=			log( TR_RATE[which(GROUP=='TRM_IN_TWM' & FACTOR=='yes')] 		/ TR_RATE[which(GROUP=='TRM_IN_TWM' & FACTOR=='no')]),
								'b_TRM_IN_SIXM_yes'=		log( TR_RATE[which(GROUP=='TRM_IN_SIXM' & FACTOR=='yes')] 		/ TR_RATE[which(GROUP=='TRM_IN_SIXM' & FACTOR=='no')]),
								'b_TRM_IN_THREEM_yes'=		log( TR_RATE[which(GROUP=='TRM_IN_THREEM' & FACTOR=='yes')] 	/ TR_RATE[which(GROUP=='TRM_IN_THREEM' & FACTOR=='no')]),
								'b_DIAG_IN_RECENT_yes'=		log( TR_RATE[which(GROUP=='DIAG_IN_RECENT' & FACTOR=='Y')] 		/ TR_RATE[which(GROUP=='DIAG_IN_RECENT' & FACTOR=='N')]),
								'b_DIAG_IN_SIXM_yes'=		log( TR_RATE[which(GROUP=='DIAG_IN_SIXM' & FACTOR=='Y')] 		/ TR_RATE[which(GROUP=='DIAG_IN_SIXM' & FACTOR=='N')]),
								'b_DIAG_IN_THREEM_yes'=		log( TR_RATE[which(GROUP=='DIAG_IN_THREEM' & FACTOR=='Y')] 		/ TR_RATE[which(GROUP=='DIAG_IN_THREEM' & FACTOR=='N')]),								
								'b_DIAG_CD4_STAGE2_g500'=	log( TR_RATE[which(GROUP=='DIAG_CD4_STAGE2' & FACTOR=='g500')] 	/ TR_RATE[which(GROUP=='DIAG_CD4_STAGE2' & FACTOR=='l500')]),
								'b_DIAG_CD4_STAGE2_unknown'=log( TR_RATE[which(GROUP=='DIAG_CD4_STAGE2' & is.na(FACTOR))] 	/ TR_RATE[which(GROUP=='DIAG_CD4_STAGE2' & FACTOR=='l500')]),								
								'b_RISK_L'=					log( TR_RATE[which(GROUP=='RISK' & FACTOR=='L')] 				/ TR_RATE[which(GROUP=='RISK' & FACTOR=='M')]),
								'b_RISK_H'=					log( TR_RATE[which(GROUP=='RISK' & FACTOR=='H')] 				/ TR_RATE[which(GROUP=='RISK' & FACTOR=='M')]),
								'b_MALE'=					log( TR_RATE[which(GROUP=='GENDER' & FACTOR=='M')] 				/ TR_RATE[which(GROUP=='GENDER' & FACTOR=='F')]),
								'b_MALE_RISKM'=				log( TR_RATE[which(GROUP=='GENDER_RISKM' & FACTOR=='M')] 		/ TR_RATE[which(GROUP=='GENDER_RISKM' & FACTOR=='F')]),
								'b_AGE_AT_DIAG_unknown'=	log( TR_RATE[which(GROUP=='AGE_AT_DIAG' & FACTOR=='unknown')] 	/ TR_RATE[which(GROUP=='AGE_AT_DIAG' & FACTOR=='30')]),
								'b_AGE_AT_DIAG_less25'=		log( TR_RATE[which(GROUP=='AGE_AT_DIAG' & FACTOR=='less25')] 	/ TR_RATE[which(GROUP=='AGE_AT_DIAG' & FACTOR=='30')]),
								'b_AGE_AT_DIAG_25to29'=		log( TR_RATE[which(GROUP=='AGE_AT_DIAG' & FACTOR=='25to29')] 	/ TR_RATE[which(GROUP=='AGE_AT_DIAG' & FACTOR=='30')])
								)		
				list(TRUE_TR_RATE_RATIO=unname(tmp), FACTOR=names(tmp))
			}, by=c('TS','TE')]	
	
	dcast.data.table(dfr,  GROUP+FACTOR~TS+TE, value.var='TR_RATE')	
	dcast.data.table(subset(dtrr, grepl('b_DIAG_IN|b_TRM_IN',FACTOR)),  FACTOR~TS+TE, value.var='TRUE_TR_RATE_RATIO')
	dcast.data.table(subset(dtrr, grepl('b_AGE_AT_DIAG',FACTOR)),  FACTOR~TS+TE, value.var='TRUE_TR_RATE_RATIO')
	
	
	save(dfr, dtrr, file=paste0(outfile,'empirical_transmission_rates.rda'))
	
	ggplot(dtrr, aes(x=paste(FACTOR, TS, TE), y=TRUE_TR_RATE_RATIO)) + 
			geom_point() + coord_flip()
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

cr.png.compare<- function()
{
	require(coalreg)
	require(viridis)
	indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_results'
	infiles	<- data.table(F=list.files(indir, pattern='rda$',full.names=TRUE))
	#	update file names if needed
	if(0)
	{
		infiles	<- subset(infiles, !grepl('scale', F))
		#invisible(infiles[,file.rename(F, gsub('\\.rda$','_lasso5.rda',F)),by='F'])
		#invisible(infiles[,file.rename(F, gsub('\\.rda$','_ETSIbias1_ETSInoise0.rda',F)),by='F'])
		invisible(infiles[,file.rename(F, gsub('\\.rda$','_scale0.rda',F)),by='F'])
	}
	#	load true values
	load("~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-empirical_transmission_rates.rda")
	dtrr	<- subset(dtrr, TS==2005 & TE==2020)
	set(dtrr, NULL, c('TS','TE'), NULL)
	#	parse results	
	infiles	<- subset(infiles, !grepl('results\\.rda',F))
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_results/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43_coalreg_using_TRM-FROM-RECENTtrf_ETSIaoi_BFGSargs3_maxNodeDepth3.rda'
				load(F)
				list(	FACTOR= names(fit$bestfit$par),
						MLE= unname(fit$bestfit$par),
						MLE_CONVERGED= fit$bestfit$convergence==0
						#PROF_MEDIAN=unname(apply(pci$sample, 2, median))
						)				
			}, by='F']
	res[, SEED:= as.integer(gsub('.*seed([0-9]+)_.*','\\1',F))]
	res[, COV:= as.numeric(gsub('.*cov([0-9]+\\.[0-9]+)-.*','\\1',F))/100]
	res[, MAXNODE:=as.numeric(gsub('.*_maxNodeDepth([0-9A-Za-z]+).*','\\1',F))]
	res[, MAXHEIGHT:=as.numeric(gsub('.*_maxHeight([0-9A-Za-z]+).*','\\1',F))]
	res[, LASSO:=as.numeric(gsub('.*_lasso([0-9]+).*','\\1',F))]
	res[, SCALE:=as.numeric(gsub('.*_scale([0-9]).*','\\1',F))]
	res[, CLIMB:= gsub('.*_(BFGS|Nelder-Mead).*','\\1',F)]
	res[, CR_VERSION:='standard']
	set(res, res[, which(grepl('coalreguh',F))], 'CR_VERSION','unmodelled heterogeneity')
	set(res, res[, which(grepl('coalregaoimodel1',F))], 'CR_VERSION','aoi model 1')
	res[, TRF_RECENT:= as.numeric(grepl('TRM-FROM-RECENTtrf',F))]
	res[, TRF_RISK:= as.numeric(grepl('RISKtrf',F))]
	res[, TRF_GENDER:= as.numeric(grepl('MALEtrf',F))]
	res[, TRF_AGEATDIAG:= as.numeric(grepl('AGEATDIAGtrf|AGE_AT_DIAGtrf',F))]	
	res[, TRF_DIAG3M:= as.numeric(grepl('DIAG_IN_THREEMtrf',F))]
	res[, TRF_DIAG6M:= as.numeric(grepl('DIAG_IN_SIXMtrf',F))]
	res[, TRF_DIAG12M:= as.numeric(grepl('DIAG_IN_RECENTtrf',F))]
	res[, TRF_DIAGCD4g500:= as.numeric(grepl('DIAG_CD4_STAGE2trf',F))]
	res[, AOI_ETSI:= as.numeric(grepl('ETSIaoi',F))]
	res[, AOI_ETSINOISE:=as.numeric(gsub('.*_ETSInoise([0-9]+\\.?[0-9]*).*','\\1',F))]
	res[, AOI_ETSIBIAS:=as.numeric(gsub('.*_ETSIbias([0-9]+\\.?[0-9]*).*','\\1',F))]	
	#
	res[, TRF:=NA_character_]
	set(res, res[, which(TRF_RECENT==1)],'TRF','recent transmission')
	set(res, res[, which(TRF_RISK==1)],'TRF','risk low, high')
	set(res, res[, which(TRF_GENDER==1)],'TRF','male')	
	set(res, res[, which(TRF_DIAG3M==1)],'TRF','diagnosed within 3 month of infection')
	set(res, res[, which(TRF_DIAG6M==1)],'TRF','diagnosed within 6 month of infection')
	set(res, res[, which(TRF_DIAG12M==1)],'TRF','diagnosed within 12 month of infection')
	set(res, res[, which(TRF_DIAGCD4g500==1)],'TRF','diagnosed with CD4>500')
	set(res, res[, which(TRF_RECENT==1 & TRF_RISK==1)],'TRF','recent transmission\nrisk low, high')
	set(res, res[, which(TRF_RECENT==1 & TRF_RISK==1 & TRF_GENDER==1)],'TRF','recent transmission\nrisk low, high\nmale')
	set(res, res[, which(TRF_RECENT==1 & TRF_RISK==1 & TRF_GENDER==1 & TRF_AGEATDIAG==1)],'TRF','recent transmission\nrisk low, high\nmale\nage at diagnosis')	
	set(res, res[, which(TRF_DIAG3M==1 & TRF_GENDER==1 & TRF_AGEATDIAG==1)],'TRF','diagnosed within 3 month of infection\nmale\nage at diagnosis')
	set(res, res[, which(TRF_DIAG6M==1 & TRF_GENDER==1 & TRF_AGEATDIAG==1)],'TRF','diagnosed within 6 month of infection\nmale\nage at diagnosis')
	set(res, res[, which(TRF_DIAG12M==1 & TRF_GENDER==1 & TRF_AGEATDIAG==1)],'TRF','diagnosed within 12 month of infection\nmale\nage at diagnosis')
	set(res, res[, which(TRF_DIAGCD4g500==1 & TRF_GENDER==1 & TRF_AGEATDIAG==1)],'TRF','diagnosed with CD4>500\nmale\nage at diagnosis')
	stopifnot(res[, !any(is.na(TRF))])
	res[, AOI:=NA_character_]
	set(res, res[, which(AOI_ETSI==1 & AOI_ETSINOISE==0 & AOI_ETSIBIAS==1)],'AOI','exact time since infection')
	tmp		<- res[, which(AOI_ETSI==1 & AOI_ETSINOISE==0 & AOI_ETSIBIAS!=1)]
	set(res, tmp,'AOI',res[tmp, paste0('time since infection\nwith mult bias ',AOI_ETSIBIAS,'\nno noise')])
	tmp		<- res[, which(AOI_ETSI==1 & AOI_ETSINOISE!=0 & AOI_ETSIBIAS!=1)]
	set(res, tmp,'AOI',res[tmp, paste0('time since infection\nwith mult bias ',AOI_ETSIBIAS,'\nwith logn noise sigma ',AOI_ETSINOISE)])
	tmp		<- res[, which(AOI_ETSI==1 & AOI_ETSINOISE!=0 & AOI_ETSIBIAS==1)]
	set(res, tmp,'AOI',res[tmp, paste0('time since infection\nno bias\nwith logn noise sigma ',AOI_ETSINOISE)])
	stopifnot(res[, !any(is.na(AOI))])
	res[, R_RANGE:= NA_character_]	
	set(res, res[, which(grepl('args3',F))], 'R_RANGE', 'r in -4,2')	
	#
	tmp		<- res[, which(TRF=='recent transmission\nrisk low, high' & FACTOR=='b_TRM_FROM_RECENT')]
	set(res, tmp, 'FACTOR', 'b_TRM_FROM_RECENT_RISKM')	
	tmp		<- res[, which(TRF=='recent transmission\nrisk low, high\nmale' & FACTOR=='b_TRM_FROM_RECENT')]
	set(res, tmp, 'FACTOR', 'b_TRM_FROM_RECENT_RISKM')
	tmp		<- res[, which(TRF=='recent transmission\nrisk low, high\nmale' & FACTOR=='b_MALE')]
	set(res, tmp, 'FACTOR', 'b_MALE_RISKM')	
	tmp		<- res[, which(TRF=='recent transmission\nrisk low, high\nmale\nage at diagnosis' & FACTOR=='b_TRM_FROM_RECENT')]
	set(res, tmp, 'FACTOR', 'b_TRM_FROM_RECENT_RISKM')	
	tmp		<- res[, which(TRF=='recent transmission\nrisk low, high\nmale\nage at diagnosis' & FACTOR=='b_MALE')]
	set(res, tmp, 'FACTOR', 'b_MALE_RISKM')
	tmp		<- res[, which(TRF=='recent transmission\nrisk low, high\nmale\nage at diagnosis' & FACTOR=='b_AGE_AT_DIAG_unknown')]
	set(res, tmp, 'FACTOR', 'b_AGE_AT_DIAG_unknown_RISKM_GENDERF')	
	tmp		<- res[, which(TRF=='recent transmission\nrisk low, high\nmale\nage at diagnosis' & FACTOR=='b_AGE_AT_DIAG_less25')]
	set(res, tmp, 'FACTOR', 'b_AGE_AT_DIAG_less25_RISKM_GENDERF')
	tmp		<- res[, which(TRF=='recent transmission\nrisk low, high\nmale\nage at diagnosis' & FACTOR=='b_AGE_AT_DIAG_25to29')]
	set(res, tmp, 'FACTOR', 'b_AGE_AT_DIAG_25to29_RISKM_GENDERF')	
	#		
	res		<- merge(res, dtrr, by='FACTOR', all.x=1)
	#
	set(res, res[, which(FACTOR=='b_DIAG_IN_RECENT_yes')], 'FACTOR', 'b_DIAG_IN_TWELVEM_yes')
	set(res, res[, which(FACTOR=='b_DIAG_IN_RECENT_unknown')], 'FACTOR', 'b_DIAG_IN_TWELVEM_unknown')	
	save(res, file=file.path(indir, 'PANGEA-AcuteHigh-InterventionNone-results.rda'))
	#
	#
	resp	<- melt(res, measure.vars=c('MLE','TRUE_TR_RATE_RATIO'), variable.name='STAT', value.name='V')
	resp	<- subset(resp, !is.na(V))
	#
	#	 compare MLEs: given exact time since infection, consider trm risk from recent infection, vary max node depth 
	#
	tmp	<- subset(resp, CR_VERSION=='standard' & CLIMB=='BFGS' & SCALE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & TRF=='recent transmission' & AOI=='exact time since infection')
	subset(tmp,FACTOR=='b_TRM_FROM_RECENT' & STAT=='MLE')[, table(MAXNODE, useNA='if')]	
	tmp[, LABEL:= factor(MAXNODE, levels= sort(unique(tmp$MAXNODE)), labels=paste0('exact time since infection\nmaxHeight 10\nmaxDepth',sort(unique(tmp$MAXNODE))))]					
	ggplot(tmp, aes(x=FACTOR, y=V, colour=STAT)) + 
			#geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			geom_point() +
			theme_bw() + 
			labs(x='', y='log risk ratio\n') +	
			coord_cartesian(ylim=c(-5,5)) +
			scale_colour_manual(values=c('MLE'='black', 'TRUE_TR_RATE_RATIO'='red')) +
			scale_y_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			facet_grid(~LABEL) +
			theme(legend.position='bottom')
	ggsave(file= file.path(indir,'compare_png_MLEs_BFGSargs3_trfRECENT_aoiETSI.pdf'), w=7,h=4)
	#
	#	compare MLEs: given exact time since infection
	#		vary trm risk from recent infection AND PERHAPS risk  
	#		fix lasso=5 maxHeight=10
	#	result: adding RISK_L and RISK_H really improves inference
	tmp	<- subset(resp, CR_VERSION=='standard' & CLIMB=='BFGS' & SCALE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & MAXNODE==Inf & AOI=='exact time since infection')
	subset(tmp,FACTOR=='b_TRM_FROM_RECENT' & STAT=='MLE')[, table(MAXNODE, TRF, useNA='if')]	
	tmp[, LABEL:= paste0(TRF,'\n\nmaxHeight ',MAXHEIGHT,'\nmaxDepth ',MAXNODE,'\nlasso ',LASSO)]
	ggplot(tmp, aes(x=V, y=FACTOR, colour=STAT)) + 
			geom_vline(xintercept=0, colour='grey50', size=1) +
			geom_point() +
			theme_bw() + 
			labs(y='', x='\nlog risk ratio') +	
			coord_cartesian(xlim=c(-7,7)) +
			scale_x_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			scale_colour_manual(values=c('MLE'='black', 'TRUE_TR_RATE_RATIO'='red')) +
			facet_grid(~LABEL) +
			theme(legend.position='bottom')
	ggsave(file= file.path(indir,'compare_png_MLEs_BFGSargs3_varyTRF_fixETSI.pdf'), w=25,h=6)
	#
	#	vary time since infection  
	#		fix lasso=5 maxHeight=10
	#	result: this is pretty good! noise + bias in the time since infection model do not have much of an impact 
	tmp	<- subset(resp, CR_VERSION=='standard' & CLIMB=='BFGS' & SCALE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & MAXNODE==Inf & MAXHEIGHT%in%c(5,10) & TRF=='recent transmission\nrisk low, high\nmale\nage at diagnosis')
	subset(tmp,FACTOR=='b_TRM_FROM_RECENT_RISKM' & STAT=='MLE')[, table(MAXHEIGHT, AOI, useNA='if')]	
	#tmp[, LABEL:= paste0(AOI,'\nmaxHeight ',MAXHEIGHT,'\nmaxDepth ',MAXNODE,'\nlasso ',LASSO)]
	tmp[, LABEL:= AOI]
	ggplot(tmp, aes(x=V, y=LABEL, colour=STAT)) + 
			geom_vline(xintercept=0, colour='grey50', size=1) +
			geom_point() +
			theme_bw() + 
			labs(y='', x='\nlog risk ratio') +	
			coord_cartesian(xlim=c(-7,7)) +
			scale_x_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			scale_colour_manual(values=c('MLE'='black', 'TRUE_TR_RATE_RATIO'='red')) +
			facet_grid(~FACTOR, scales='free_x') +
			theme(legend.position='bottom')
	ggsave(file= file.path(indir,'compare_png_MLEs_BFGSargs3_fixTRF_comapre_noisebias_in_timesinceinfection.pdf'), w=25,h=6)
	#
	#	check unmodelled heterogeneity
	#
	tmp	<- subset(resp, MLE_CONVERGED & CLIMB=='BFGS' & SCALE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & MAXNODE==Inf & MAXHEIGHT==10 & AOI=='exact time since infection' & !grepl('diagnosed',TRF))
	tmp	<- subset(tmp, TRF%in%c('recent transmission','risk low, high'))
	subset(tmp,FACTOR=='b_TRM_FROM_RECENT' & STAT=='MLE')[, table(TRF, CR_VERSION, useNA='if')]
	tmp[, LABEL:= paste0(TRF,'\n',CR_VERSION)]
	ggplot(tmp, aes(x=V, y=FACTOR, colour=STAT)) + 
			geom_vline(xintercept=0, colour='grey50', size=1) +
			geom_point() +
			theme_bw() + 
			labs(y='', x='\nlog risk ratio') +	
			coord_cartesian(xlim=c(-3,3)) +
			scale_x_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			scale_colour_manual(values=c('MLE'='black', 'TRUE_TR_RATE_RATIO'='red')) +
			facet_grid(~LABEL) +
			theme(legend.position='bottom')
	ggsave(file= file.path(indir,'compare_png_MLEs_BFGSargs3_varyTRF_fixETSI_compare_unmodelledheterogeneity.pdf'), w=15,h=6)
	#
	#	check unmodelled heterogeneity
	#
	tmp	<- subset(resp, MLE_CONVERGED & CLIMB=='BFGS' & SCALE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & MAXNODE==Inf & MAXHEIGHT==10 & AOI=='exact time since infection' & grepl('diagnosed',TRF))
	set(tmp, NULL, 'LABEL', tmp[, factor(TRF, levels=c(	'diagnosed within 3 month of infection','diagnosed within 6 month of infection','diagnosed within 12 month of infection','diagnosed with CD4>500',
														'diagnosed within 3 month of infection\nmale\nage at diagnosis','diagnosed within 6 month of infection\nmale\nage at diagnosis','diagnosed within 12 month of infection\nmale\nage at diagnosis','diagnosed with CD4>500\nmale\nage at diagnosis'
														))])
	subset(tmp, STAT=='MLE')[, table(TRF, CR_VERSION, useNA='if')]	
	ggplot(tmp, aes(x=V, y=FACTOR, colour=STAT)) + 
			geom_vline(xintercept=0, colour='grey50', size=1) +
			geom_point() +
			theme_bw() + 
			labs(y='', x='\nlog risk ratio') +	
			coord_cartesian(xlim=c(-3,3)) +
			scale_x_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			scale_colour_manual(values=c('MLE'='black', 'TRUE_TR_RATE_RATIO'='red')) +
			facet_grid(CR_VERSION~LABEL) +
			theme(legend.position='bottom')
	ggsave(file= file.path(indir,'compare_png_MLEs_BFGSargs3_varyObservedTRF_fixETSI_compare_aoimodel1_unmodelledheterogeneity.pdf'), w=20,h=12)
	
	
	
	#
	#	compare MLEs: given exact time since infection, given trm recent infection and risk 
	#		vary lasso
	#	result: lasso has large impact, lasso=5 much better than lasso=10
	#			maxDepth=3 does not lead to significantly worse performace
	#			maxHeight=5 does not lead to significantly worse performace
	tmp	<- subset(resp, CR_VERSION=='standard' & CLIMB=='BFGS' & SCALE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & AOI=='exact time since infection')
	subset(tmp,FACTOR=='b_TRM_FROM_RECENT_RISKM' & STAT=='MLE')[, table(MAXNODE, MAXHEIGHT, LASSO, useNA='if')]	
	tmp[, LABEL:= paste0(TRF,'\nmaxHeight ',MAXHEIGHT,'\nmaxDepth ',MAXNODE,'\nlasso ',LASSO)]
	ggplot(tmp, aes(x=V, y=FACTOR, colour=STAT)) + 
			geom_vline(xintercept=0, colour='grey50', size=1) +
			geom_point() +
			theme_bw() + 
			labs(y='', x='\nlog risk ratio') +	
			coord_cartesian(xlim=c(-10,10)) +
			scale_x_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			scale_colour_manual(values=c('MLE'='black', 'TRUE_TR_RATE_RATIO'='red')) +
			facet_grid(~LABEL) +
			theme(legend.position='bottom')
	ggsave(file= file.path(indir,'compare_png_MLEs_BFGSargs3_trfRECENT_trfRISK_aoiETSI_vary_LASSO_MAXHEIGHT_MAXDEPTH.pdf'), w=30,h=6)
	#
	#	compare CLIMB / SCALE 
	#
	tmp	<- subset(resp, CR_VERSION=='standard' & AOI_ETSINOISE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & MAXNODE==Inf & MAXHEIGHT==10 & TRF=='recent transmission\nrisk low, high' & MLE_CONVERGED)
	subset(tmp,FACTOR=='b_TRM_FROM_RECENT_RISKM' & STAT=='MLE')[, table(CLIMB, SCALE, useNA='if')]	
	tmp[, LABEL:= paste0(AOI,'\nclimb ',CLIMB,'\nscale ',SCALE,'\nlasso ',LASSO,'\nmaxHeight ',MAXHEIGHT)]
	ggplot(tmp, aes(x=V, y=LABEL, colour=STAT)) + 
			geom_vline(xintercept=0, colour='grey50', size=1) +
			geom_point() +
			theme_bw() + 
			labs(y='', x='\nlog risk ratio') +	
			coord_cartesian(xlim=c(-7,7)) +
			scale_x_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			scale_colour_manual(values=c('MLE'='black', 'TRUE_TR_RATE_RATIO'='red')) +
			facet_grid(AOI_ETSIBIAS~FACTOR, scales='free') +
			theme(legend.position='bottom')
	ggsave(file= file.path(indir,'compare_png_MLEs_BFGSargs3_modelWRISKRECENT_aoiETSINOISE_checkCLIMB.pdf'), w=15,h=15)	
	
	tmp	<- subset(resp, AOI_ETSINOISE==0 & AOI_ETSIBIAS==1 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & MAXNODE==Inf & MAXHEIGHT==10 & TRF=='recent transmission\nrisk low, high' & AOI=='exact time since infection' & SCALE==0)
	load('/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_results/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed51_coalreg_using_TRM-FROM-RECENTtrf_RISKtrf_ETSIaoi_BFGSargs3_maxNodeDepthInf_maxHeight10_lasso5_ETSIbias1_ETSInoise0_scale0.rda')
	fit.bfgs	<- fit
	load('~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/tmp/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed51_coalreg_using_TRM-FROM-RECENTtrf_RISKtrf_ETSIaoi_BFGSargs3_maxNodeDepthInf_maxHeight10_lasso5_ETSIbias1_ETSInoise0_scale0.rda')
	fit.bfgs.n	<- fit
	load('/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_results/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed51_coalreg_using_TRM-FROM-RECENTtrf_RISKtrf_ETSIaoi_Nelder-Meadargs3_maxNodeDepthInf_maxHeight10_lasso5_ETSIbias1_ETSInoise0_scale0.rda')
	fit.nm		<- fit
	lapply(list(fit.bfgs, fit.bfgs.n, fit.nm), function(x) x[['theta0']])
	lapply(list(fit.bfgs, fit.bfgs.n, fit.nm), function(x) x[['theta_starts']])
	lapply(list(fit.bfgs, fit.bfgs.n, fit.nm), function(x) x[['method']])
	lapply(list(fit.bfgs, fit.bfgs.n, fit.nm), function(x) x[['s0']])
	lapply(list(fit.bfgs, fit.bfgs.n, fit.nm), function(x) x[['b0']])
	lapply(list(fit.bfgs, fit.bfgs.n, fit.nm), function(x) x[['X']][c("ID_51400","ID_51453","ID_88472","ID_47654","ID_35266","ID_23203","ID_19449","ID_34289","ID_48652","ID_33527"),])
	lapply(list(fit.bfgs, fit.bfgs.n, fit.nm), function(x) x[['bdts']])
	lapply(list(fit.bfgs, fit.bfgs.n, fit.nm), function(x) x[['bestfit']])
	lapply(list(fit.bfgs, fit.bfgs.n, fit.nm), function(x) x[['fits']])
	
	lapply(list(fit.bfgs, fit.bfgs.n, fit.nm), function(x) x[['bdts']])
	
	fit.bfgs[['bdts']]
	#
	#	MLE correlations
	#
	tmp	<- subset(resp, CLIMB=='BFGS' & SCALE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & MAXNODE==Inf & TRF=='recent transmission\nrisk low, high\nmale\nage at diagnosis' & AOI=='exact time since infection')
	tmp	<- dcast.data.table(subset(tmp, STAT=='MLE'), SEED~FACTOR, value.var='V')
	setnames(tmp, colnames(tmp), gsub('_','\n',gsub('b_','',colnames(tmp))))
	pdf(file=file.path(indir,'corr_trfRECENT_trfRISK_trfMALE_trfAGEATDIAG.pdf'), w=10, h=10)
	pairs(tmp[, 2:10], col='black', pch=16)
	dev.off()
	
	tmp	<- subset(resp, CLIMB=='BFGS' & SCALE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & MAXHEIGHT==10 & MAXNODE==Inf & TRF=='recent transmission\nrisk low, high\nmale' & AOI=='exact time since infection')
	tmp	<- dcast.data.table(subset(tmp, STAT=='MLE'), SEED~FACTOR, value.var='V')
	setnames(tmp, colnames(tmp), gsub('_','\n',gsub('b_','',colnames(tmp))))
	pdf(file=file.path(indir,'corr_trfRECENT_trfRISK_trfMALE.pdf'), w=7, h=7)
	pairs(tmp[, 2:7], col='black', pch=16)
	dev.off()
	
	tmp	<- subset(resp, CLIMB=='BFGS' & SCALE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & MAXHEIGHT==10 & MAXNODE==Inf & TRF=='recent transmission\nrisk low, high' & AOI=='exact time since infection')
	tmp	<- dcast.data.table(subset(tmp, STAT=='MLE'), SEED~FACTOR, value.var='V')
	setnames(tmp, colnames(tmp), gsub('_','\n',gsub('b_','',colnames(tmp))))
	pdf(file=file.path(indir,'corr_trfRECENT_trfRISK.pdf'), w=5, h=5)
	pairs(tmp[, 2:6], col='black', pch=16)
	dev.off()
	
	tmp	<- subset(resp, CLIMB=='BFGS' & SCALE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & MAXHEIGHT==10 & MAXNODE==Inf & TRF=='recent transmission' & AOI=='exact time since infection')
	tmp	<- dcast.data.table(subset(tmp, STAT=='MLE'), SEED~FACTOR, value.var='V')
	setnames(tmp, colnames(tmp), gsub('_','\n',gsub('b_','',colnames(tmp))))
	pdf(file=file.path(indir,'corr_trfRECENT.pdf'), w=4, h=4)
	pairs(tmp[, 2:4], col='black', pch=16)
	dev.off()
	
	tmp	<- subset(resp, CLIMB=='BFGS' & SCALE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & MAXHEIGHT==10 & MAXNODE==Inf & TRF=='male' & AOI=='exact time since infection')
	tmp	<- dcast.data.table(subset(tmp, STAT=='MLE'), SEED~FACTOR, value.var='V')
	setnames(tmp, colnames(tmp), gsub('_','\n',gsub('b_','',colnames(tmp))))
	pdf(file=file.path(indir,'corr_trfMALE.pdf'), w=4, h=4)
	pairs(tmp[, 2:4], col='black', pch=16)
	dev.off()
	
	tmp	<- subset(resp, CLIMB=='BFGS' & SCALE==0 & R_RANGE=='r in -4,2' & STAT%in%c('MLE','TRUE_TR_RATE_RATIO') & LASSO==5 & MAXHEIGHT==10 & MAXNODE==Inf & TRF=='risk low, high' & AOI=='exact time since infection')
	tmp	<- dcast.data.table(subset(tmp, STAT=='MLE'), SEED~FACTOR, value.var='V')
	setnames(tmp, colnames(tmp), gsub('_','\n',gsub('b_','',colnames(tmp))))
	pdf(file=file.path(indir,'corr_trfRISK.pdf'), w=5, h=5)
	pairs(tmp[, 2:5], col='black', pch=16)
	dev.off()
	
}
	
cr.master.ex3.compare<- function()
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
	set(res, res[, which(grepl('BFGSargs3',F) & !grepl('maxNodeDepthInf',F))], 'MAXNODE', 2)
	tmp	<- res[, which(grepl('maxNodeDepth[0-9]+',F))]
	set(res, tmp, 'MAXNODE', res[tmp,as.numeric(gsub('.*maxNodeDepth([0-9]+)_.*','\\1',F))])
	res[, MAXHEIGHT:=Inf]
	tmp	<- res[, which(grepl('maxHeight[0-9]+',F))]
	set(res, tmp, 'MAXHEIGHT', res[tmp,as.numeric(gsub('.*maxHeight([0-9]+)_.*','\\1',F))])	
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
	tmp	<- subset(resp, BFGSargs=='r in -4,2' & STAT=='MLE' & TAXA_SAMPLED==.5 & NOISE_MODEL%in%c('lognormal','none') & MAXNODE==Inf & MAXHEIGHT==Inf)
	subset(tmp,THETA=='b_TYPE')[, table(N, TSI_COEFF, useNA='if')]
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
	tmp	<- subset(resp, BFGSargs=='r in -4,2' & grepl('BFGSargs3',F) & STAT=='MLE' & TAXA_SAMPLED==.5 & NOISE_MODEL%in%c('lognormal','none') & MAXHEIGHT==Inf)
	subset(tmp,THETA=='b_TYPE')[, table(N, MAXNODE, useNA='if')]	
	tmp[, LABEL:= factor(MAXNODE, levels= sort(unique(tmp$MAXNODE)), labels=paste0('exact time since infection\nmaxNodeDepth\n',sort(unique(tmp$MAXNODE))))]					
	ggplot(tmp, aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +	
			coord_cartesian(ylim=c(-5,5)) +
			scale_y_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			facet_grid(N~LABEL) 
	ggsave(file= file.path(indir,'compare_MLEs_BFGSargs3_maxNodeDepth.pdf'), w=10,h=7)	
	#
	#	compare maxHeight results
	#
	tmp	<- subset(resp, BFGSargs=='r in -4,2' & grepl('BFGSargs3',F) & STAT=='MLE' & TAXA_SAMPLED==.5 & NOISE_MODEL%in%c('lognormal','none') & MAXNODE==Inf)
	subset(tmp,THETA=='b_TYPE')[, table(N, MAXHEIGHT, useNA='if')]	
	tmp[, LABEL:= factor(MAXHEIGHT, levels= sort(unique(tmp$MAXHEIGHT)), labels=paste0('exact time since infection\nmaxHeight\n',sort(unique(tmp$MAXHEIGHT))))]					
	ggplot(tmp, aes(x=THETA, y=V)) + 
			geom_violin(trim=TRUE, scale='width') + geom_boxplot(outlier.shape=NA, fill="#482576FF", width=0.2, size=0.3, colour="#FCA50AFF") +
			theme_bw() + 
			labs(x='', y='log risk ratio I1 vs baseline I0\n') +	
			coord_cartesian(ylim=c(-5,5)) +
			scale_y_continuous(breaks=seq(-10,10,1), expand=c(0,0)) +
			facet_grid(N~LABEL) 
	ggsave(file= file.path(indir,'compare_MLEs_BFGSargs3_maxHeight.pdf'), w=10,h=4)	
	

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

cr.png.evalcoalreg.using.TRSTAGE.ETSI.vanilla.BFGS3<- function()
{
	require(coalreg)
	load('~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-empirical_transmission_rates.rda')
	subset(dfr, GROUP=='TR_STAGE' & SC=='PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43')
	
	infiles	<- data.table(F=list.files('/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_results', pattern='*rda', full.names=TRUE))
	
	load(infiles[1,F])
	f.maxD3<- fit$bestfit
	load(infiles[2,F])
	f.maxDInf<- fit$bestfit
}

cr.png.runcoalreg.using.TRSTAGE.ETSI.vanilla.BFGS3<- function(indir, par.base.pattern, par.maxNodeDepth=3, par.maxHeight=10, par.lasso=5)
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
		par.maxNodeDepth	<- 3
		par.maxHeight		<- 10
		par.lasso			<- 5
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#		 
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
										IDPOP=paste0('ID_',sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
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
				ph 				<- DatedTree( ph, setNames( phi$TIME_SEQ, ph$tip.label) ) 
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
											maxNodeDepth=par.maxNodeDepth,
											maxHeight=par.maxHeight,
											lasso_threshold=par.lasso, method = 'BFGS', lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )
				tmp		<- file.path(dirname(F), gsub('\\.newick',paste0('_coalreg_using_TRM-FROM-RECENTtrf_ETSIaoi_BFGSargs3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_lasso',par.lasso,'.rda'),basename(F)))				
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.png.runcoalreg.using.TRSTAGE.ETSI.unmodelled.het.BFGS3<- function(indir, par.base.pattern, par.maxNodeDepth=3, par.maxHeight=10, par.lasso=5, par.climb='BFGS', par.bias=1, par.noise=0, par.scale=0)
{
	require(coalreguh)
	require(data.table)
	require(ape)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_results'
		par.base.pattern	<- 'PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43'
		par.maxNodeDepth	<- 3
		par.maxHeight		<- 10
		par.lasso			<- 5
		par.climb			<- 'BFGS' 
		par.bias			<- 1 
		par.noise			<- 0 
		par.scale			<- 0
	}
	stopifnot(par.bias==1, par.noise==0)
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#		 
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
						IDPOP=paste0('ID_',sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
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
				ph 				<- DatedTree( ph, setNames( phi$TIME_SEQ, ph$tip.label) ) 
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
						maxNodeDepth=par.maxNodeDepth,
						maxHeight=par.maxHeight,
						lasso_threshold=par.lasso, 
						method=unname(par.climb), 
						lnr0 = -2, lnrLimits = c(-4, 2), 
						scale=as.logical(par.scale))							
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )
				tmp		<- file.path(dirname(F), gsub('\\.newick',paste0('_coalreguh_using_TRM-FROM-RECENTtrf_ETSIaoi_',par.climb,'args3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_lasso',par.lasso,'_ETSIbias',par.bias,'_ETSInoise',par.noise,'_scale',par.scale,'.rda'),basename(F)))								
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.png.runcoalreg.using.TRvariable.ETSInoise.unmodelled.het<- function(indir, par.base.pattern, par.maxNodeDepth=3, par.maxHeight=10, par.lasso=5, par.climb='BFGS', par.bias=1, par.noise=0, par.scale=0, trm.factors= c('DIAG_IN_RECENT','MALE','AGE_AT_DIAG'))
{
	require(coalreguh)
	require(data.table)
	require(ape)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_results'
		par.base.pattern	<- 'PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43'
		par.maxNodeDepth	<- 3
		par.maxHeight		<- 10
		par.lasso			<- 5
		par.climb			<- 'BFGS' 
		par.bias			<- 1 
		par.noise			<- 0 
		par.scale			<- 0
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#		 
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
										TAXA_ID=seq_along(ph$tip.label),
										IDPOP=paste0('ID_',sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
										GENDER=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',2),
										DOB=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',3)),
										RISK=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',4),
										TIME_TR=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',5)),
										DIAG_T=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',6)),
										DIAG_CD4=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',7)),																				
										DIAG_IN_ACUTE=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',9),
										TIME_SEQ=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',10)),
										ETSI=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',11)),
										SAMPLED_TR=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',12),
										TRM_FROM_RECENT=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',13),
										TRM_FROM_ACUTE=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',14),
										SEQ_COV_2020=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',15))	
										)				
				set(phi, NULL, 'AGE_AT_DIAG', phi[, cut(DIAG_T-DOB, breaks=c(-1,25,30,1e4),labels=c('less25','25to29','30'))])								
				set(phi, NULL, 'AGE_AT_DIAG_25to29', phi[,factor(!(!is.na(AGE_AT_DIAG) & AGE_AT_DIAG=='25to29'),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'AGE_AT_DIAG_less25', phi[,factor(!(!is.na(AGE_AT_DIAG) & AGE_AT_DIAG=='less25'),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'AGE_AT_DIAG_unknown', phi[,factor(!is.na(AGE_AT_DIAG),levels=c(TRUE,FALSE),labels=c('N','Y'))])				
				set(phi, NULL, 'MALE', phi[,factor(GENDER,levels=c('F','M'))])
				set(phi, NULL, 'DIAG_CD4_STAGE2', phi[, cut(DIAG_CD4, breaks=c(-1,500,1e4),labels=c('l500','g500'))])			
				set(phi, NULL, 'DIAG_CD4_STAGE2_unknown', phi[,factor(!is.na(DIAG_CD4_STAGE2),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'DIAG_CD4_STAGE2_g500', phi[,factor(!(!is.na(DIAG_CD4_STAGE2) & DIAG_CD4_STAGE2=='g500') ,levels=c(TRUE,FALSE),labels=c('N','Y'))])								
				set(phi, NULL, 'DIAG_IN_RECENT', phi[, as.character(factor((DIAG_T-TIME_TR)<1, levels=c(TRUE,FALSE),labels=c('Y','N')))])
				set(phi, NULL, 'DIAG_IN_RECENT_unknown', phi[,factor(!is.na(DIAG_IN_RECENT),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'DIAG_IN_RECENT_yes', phi[,factor(!(!is.na(DIAG_IN_RECENT) & DIAG_IN_RECENT=='Y') ,levels=c(TRUE,FALSE),labels=c('N','Y'))])				
				set(phi, NULL, 'DIAG_IN_SIXM', phi[, as.character(factor((DIAG_T-TIME_TR)<.5, levels=c(TRUE,FALSE),labels=c('Y','N')))])
				set(phi, NULL, 'DIAG_IN_SIXM_unknown', phi[,factor(!is.na(DIAG_IN_SIXM),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'DIAG_IN_SIXM_yes', phi[,factor(!(!is.na(DIAG_IN_SIXM) & DIAG_IN_SIXM=='Y') ,levels=c(TRUE,FALSE),labels=c('N','Y'))])				
				set(phi, NULL, 'DIAG_IN_THREEM', phi[, as.character(factor((DIAG_T-TIME_TR)<3/12, levels=c(TRUE,FALSE),labels=c('Y','N')))])
				set(phi, NULL, 'DIAG_IN_THREEM_unknown', phi[,factor(!is.na(DIAG_IN_THREEM),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'DIAG_IN_THREEM_yes', phi[,factor(!(!is.na(DIAG_IN_THREEM) & DIAG_IN_THREEM=='Y') ,levels=c(TRUE,FALSE),labels=c('N','Y'))])
				#	add noise and bias
				tmp		<- phi[, list(ETSI_NOISE=rlnorm(1, meanlog= log(ETSI)-0.5*log(par.noise*par.noise+1), sdlog= sqrt(log(par.noise*par.noise+1)))), by='TAXA']				
				phi		<- merge(phi, tmp, by='TAXA')					
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE*par.bias])
				#	zero mean ETSI_NOISE so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE-mean(ETSI_NOISE)])
				setkey(phi, TAXA_ID)
				#	prepare coalreg input
				ph$tip.label	<- phi[, IDPOP]
				#	select columns
				tmp		<- 'ETSI_NOISE'
				if(any(trm.factors=='DIAG_IN_THREEM'))
					tmp	<- c(tmp, 'DIAG_IN_THREEM_unknown', 'DIAG_IN_THREEM_yes')
				if(any(trm.factors=='DIAG_IN_SIXM'))
					tmp	<- c(tmp, 'DIAG_IN_SIXM_unknown', 'DIAG_IN_SIXM_yes')
				if(any(trm.factors=='DIAG_IN_RECENT'))
					tmp	<- c(tmp, 'DIAG_IN_RECENT_unknown', 'DIAG_IN_RECENT_yes')
				if(any(trm.factors=='DIAG_CD4_STAGE2'))
					tmp	<- c(tmp, 'DIAG_CD4_STAGE2_unknown', 'DIAG_CD4_STAGE2_g500')
				if(any(trm.factors=='MALE'))
					tmp	<- c(tmp, 'MALE')
				if(any(trm.factors=='AGE_AT_DIAG'))
					tmp	<- c(tmp, 'AGE_AT_DIAG_unknown', 'AGE_AT_DIAG_less25', 'AGE_AT_DIAG_25to29')
				tmp				<- data.matrix(subset(phi, select=tmp))
				rownames(tmp)	<- phi[, IDPOP]
				stopifnot(!any(is.na(tmp)))
				ph 				<- DatedTree( ph, setNames( phi$TIME_SEQ, ph$tip.label) ) 
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
				fit 	<- trf.lasso(ph, tmp, trf_names=colnames(tmp)[!grepl('ETSI',colnames(tmp))], aoi_names=colnames(tmp)[grepl('ETSI',colnames(tmp))], 
						maxNodeDepth=par.maxNodeDepth,
						maxHeight=par.maxHeight,
						lasso_threshold=par.lasso, 
						method=unname(par.climb), 
						lnr0 = -2, lnrLimits = c(-4, 2), 
						scale=as.logical(par.scale))							
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )
				tmp		<- file.path(dirname(F), gsub('\\.newick',paste0('_coalreguh_using_',paste(trm.factors,'trf',collapse='_',sep=''),'_ETSIaoi_',par.climb,'args3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_lasso',par.lasso,'_ETSIbias',par.bias,'_ETSInoise',par.noise,'_scale',par.scale,'.rda'),basename(F)))								
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.png.runcoalreg.using.TRvariable.ETSInoise.aoi.model1<- function(indir, par.base.pattern, par.maxNodeDepth=3, par.maxHeight=10, par.lasso=5, par.climb='BFGS', par.bias=1, par.noise=0, par.scale=0, trm.factors= c('DIAG_IN_RECENT','MALE','AGE_AT_DIAG'))
{
	require(coalregaoiModel1)
	require(data.table)
	require(ape)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_results'
		par.base.pattern	<- 'PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43'
		par.maxNodeDepth	<- 3
		par.maxHeight		<- 10
		par.lasso			<- 5
		par.climb			<- 'BFGS' 
		par.bias			<- 1 
		par.noise			<- 0 
		par.scale			<- 0
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#		 
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
						TAXA_ID=seq_along(ph$tip.label),
						IDPOP=paste0('ID_',sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
						GENDER=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',2),
						DOB=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',3)),
						RISK=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',4),
						TIME_TR=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',5)),
						DIAG_T=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',6)),
						DIAG_CD4=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',7)),																				
						DIAG_IN_ACUTE=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',9),
						TIME_SEQ=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',10)),
						ETSI=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',11)),
						SAMPLED_TR=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',12),
						TRM_FROM_RECENT=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',13),
						TRM_FROM_ACUTE=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',14),
						SEQ_COV_2020=as.numeric(sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',15))	
				)				
				set(phi, NULL, 'AGE_AT_DIAG', phi[, cut(DIAG_T-DOB, breaks=c(-1,25,30,1e4),labels=c('less25','25to29','30'))])								
				set(phi, NULL, 'AGE_AT_DIAG_25to29', phi[,factor(!(!is.na(AGE_AT_DIAG) & AGE_AT_DIAG=='25to29'),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'AGE_AT_DIAG_less25', phi[,factor(!(!is.na(AGE_AT_DIAG) & AGE_AT_DIAG=='less25'),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'AGE_AT_DIAG_unknown', phi[,factor(!is.na(AGE_AT_DIAG),levels=c(TRUE,FALSE),labels=c('N','Y'))])				
				set(phi, NULL, 'MALE', phi[,factor(GENDER,levels=c('F','M'))])
				set(phi, NULL, 'DIAG_CD4_STAGE2', phi[, cut(DIAG_CD4, breaks=c(-1,500,1e4),labels=c('l500','g500'))])			
				set(phi, NULL, 'DIAG_CD4_STAGE2_unknown', phi[,factor(!is.na(DIAG_CD4_STAGE2),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'DIAG_CD4_STAGE2_g500', phi[,factor(!(!is.na(DIAG_CD4_STAGE2) & DIAG_CD4_STAGE2=='g500') ,levels=c(TRUE,FALSE),labels=c('N','Y'))])								
				set(phi, NULL, 'DIAG_IN_RECENT', phi[, as.character(factor((DIAG_T-TIME_TR)<1, levels=c(TRUE,FALSE),labels=c('Y','N')))])
				set(phi, NULL, 'DIAG_IN_RECENT_unknown', phi[,factor(!is.na(DIAG_IN_RECENT),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'DIAG_IN_RECENT_yes', phi[,factor(!(!is.na(DIAG_IN_RECENT) & DIAG_IN_RECENT=='Y') ,levels=c(TRUE,FALSE),labels=c('N','Y'))])				
				set(phi, NULL, 'DIAG_IN_SIXM', phi[, as.character(factor((DIAG_T-TIME_TR)<.5, levels=c(TRUE,FALSE),labels=c('Y','N')))])
				set(phi, NULL, 'DIAG_IN_SIXM_unknown', phi[,factor(!is.na(DIAG_IN_SIXM),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'DIAG_IN_SIXM_yes', phi[,factor(!(!is.na(DIAG_IN_SIXM) & DIAG_IN_SIXM=='Y') ,levels=c(TRUE,FALSE),labels=c('N','Y'))])				
				set(phi, NULL, 'DIAG_IN_THREEM', phi[, as.character(factor((DIAG_T-TIME_TR)<3/12, levels=c(TRUE,FALSE),labels=c('Y','N')))])
				set(phi, NULL, 'DIAG_IN_THREEM_unknown', phi[,factor(!is.na(DIAG_IN_THREEM),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'DIAG_IN_THREEM_yes', phi[,factor(!(!is.na(DIAG_IN_THREEM) & DIAG_IN_THREEM=='Y') ,levels=c(TRUE,FALSE),labels=c('N','Y'))])
				#	add noise and bias
				tmp		<- phi[, list(ETSI_NOISE=rlnorm(1, meanlog= log(ETSI)-0.5*log(par.noise*par.noise+1), sdlog= sqrt(log(par.noise*par.noise+1)))), by='TAXA']				
				phi		<- merge(phi, tmp, by='TAXA')					
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE*par.bias])
				#	zero mean ETSI_NOISE so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE-mean(ETSI_NOISE)])
				setkey(phi, TAXA_ID)
				#	prepare coalreg input
				ph$tip.label	<- phi[, IDPOP]
				#	select columns
				tmp		<- 'ETSI_NOISE'
				if(any(trm.factors=='DIAG_IN_THREEM'))
					tmp	<- c(tmp, 'DIAG_IN_THREEM_unknown', 'DIAG_IN_THREEM_yes')
				if(any(trm.factors=='DIAG_IN_SIXM'))
					tmp	<- c(tmp, 'DIAG_IN_SIXM_unknown', 'DIAG_IN_SIXM_yes')
				if(any(trm.factors=='DIAG_IN_RECENT'))
					tmp	<- c(tmp, 'DIAG_IN_RECENT_unknown', 'DIAG_IN_RECENT_yes')
				if(any(trm.factors=='DIAG_CD4_STAGE2'))
					tmp	<- c(tmp, 'DIAG_CD4_STAGE2_unknown', 'DIAG_CD4_STAGE2_g500')
				if(any(trm.factors=='MALE'))
					tmp	<- c(tmp, 'MALE')
				if(any(trm.factors=='AGE_AT_DIAG'))
					tmp	<- c(tmp, 'AGE_AT_DIAG_unknown', 'AGE_AT_DIAG_less25', 'AGE_AT_DIAG_25to29')
				tmp				<- data.matrix(subset(phi, select=tmp))
				rownames(tmp)	<- phi[, IDPOP]
				stopifnot(!any(is.na(tmp)))
				ph 				<- DatedTree( ph, setNames( phi$TIME_SEQ, ph$tip.label) ) 
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
				fit 	<- trf.lasso(ph, tmp, trf_names=colnames(tmp)[!grepl('ETSI',colnames(tmp))], aoi_names=colnames(tmp)[grepl('ETSI',colnames(tmp))], 
						maxNodeDepth=par.maxNodeDepth,
						maxHeight=par.maxHeight,
						lasso_threshold=par.lasso, 
						method=unname(par.climb), 
						lnr0 = -2, lnrLimits = c(-4, 2), 
						scale=as.logical(par.scale))							
				#fci 	<- fisher.ci(fit)	 
				#pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )
				tmp		<- file.path(dirname(F), gsub('\\.newick',paste0('_coalregaoimodel1_using_',paste(trm.factors,'trf',collapse='_',sep=''),'_ETSIaoi_',par.climb,'args3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_lasso',par.lasso,'_ETSIbias',par.bias,'_ETSInoise',par.noise,'_scale',par.scale,'.rda'),basename(F)))								
				save( fit, file=tmp)
			}, by='F']				
}

cr.png.runcoalreg.using.TRSTAGE.TRRISK.ETSI.vanilla.BFGS3<- function(indir, par.base.pattern, par.maxNodeDepth=3, par.maxHeight=10, par.lasso=5)
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
		par.maxNodeDepth	<- 3
		par.maxHeight		<- 10
		par.lasso			<- 5
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#		 
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
						IDPOP=paste0('ID_',sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
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
				set(phi, NULL, 'RISK_L', phi[,factor(RISK!='L',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'RISK_H', phi[,factor(RISK!='H',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				#	zero mean ETSI so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI', phi[, ETSI-mean(ETSI)])
				#	prepare coalreg input
				ph$tip.label	<- phi[, IDPOP]
				tmp				<- data.matrix(subset(phi, select=c(TRM_FROM_RECENT,RISK_L,RISK_H,ETSI)))
				rownames(tmp)	<- phi[, IDPOP]
				stopifnot(!any(is.na(tmp)))
				ph 				<- DatedTree( ph, setNames( phi$TIME_SEQ, ph$tip.label) ) 
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
				fit 	<- trf.lasso(ph, tmp, trf_names=c('TRM_FROM_RECENT','RISK_L','RISK_H'), aoi_names=c( 'ETSI' ), 
											maxNodeDepth=par.maxNodeDepth,
											maxHeight=par.maxHeight,
											lasso_threshold=par.lasso, 
											method = 'BFGS', lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )
				tmp		<- file.path(dirname(F), gsub('\\.newick',paste0('_coalreg_using_TRM-FROM-RECENTtrf_RISKtrf_ETSIaoi_BFGSargs3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_lasso',par.lasso,'.rda'),basename(F)))				
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.png.runcoalreg.using.TRSTAGE.TRRISK.ETSI.noise.BFGS3<- function(indir, par.base.pattern, par.maxNodeDepth=3, par.maxHeight=10, par.lasso=5, par.noise=1, par.bias=1, par.climb='BFGS', par.scale=0)
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
		par.maxNodeDepth	<- 3
		par.maxHeight		<- 10
		par.lasso			<- 5
		par.climb			<- 'BFGS'
		par.scale			<- 0
		par.noise			<- 1
		par.bias			<- 1
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#		 
	set.seed(42)	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'.*.newick'),full.names=TRUE))
	infiles[, {
				#F		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed51.newick'
				ph		<- read.tree( F )
				ph		<- multi2di(ladderize(ph),random=FALSE)
				#	create data.table with infection type				
				#	paste(	IDREC,GENDER,DOB,RISK,round(TIME_TR,d=3),round(DIAG_T,d=3),DIAG_CD4,
				#			DIAG_IN_RECENT,DIAG_IN_ACUTE,round(TIME_SEQ,d=3),round(ETSI,d=3),SAMPLED_TR,TRM_FROM_RECENT,TRM_FROM_ACUTE,SEQ_COV_2020,sep='|')				
				phi		<- data.table(	TAXA=ph$tip.label,
						TAXA_ID=seq_along(ph$tip.label),
						IDPOP=paste0('ID_',sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
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
				set(phi, NULL, 'RISK_L', phi[,factor(RISK!='L',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'RISK_H', phi[,factor(RISK!='H',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				#	add noise and bias
				tmp		<- phi[, list(ETSI_NOISE=rlnorm(1, meanlog= log(ETSI)-0.5*log(par.noise*par.noise+1), sdlog= sqrt(log(par.noise*par.noise+1)))), by='TAXA']				
				phi		<- merge(phi, tmp, by='TAXA')					
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE*par.bias])
				#	zero mean ETSI_NOISE so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE-mean(ETSI_NOISE)])
				setkey(phi, TAXA_ID)
				#	prepare coalreg input
				ph$tip.label	<- phi[, IDPOP]
				tmp				<- data.matrix(subset(phi, select=c(TRM_FROM_RECENT,RISK_L,RISK_H,ETSI_NOISE)))
				rownames(tmp)	<- phi[, IDPOP]
				stopifnot(!any(is.na(tmp)))
				ph 				<- DatedTree( ph, setNames( phi$TIME_SEQ, ph$tip.label) ) 
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
				fit 	<- trf.lasso(ph, tmp, trf_names=c('TRM_FROM_RECENT','RISK_L','RISK_H'), aoi_names=c( 'ETSI_NOISE' ), 
						maxNodeDepth=par.maxNodeDepth,
						maxHeight=par.maxHeight,
						lasso_threshold=par.lasso, 
						method=par.climb, 
						lnr0= -2, lnrLimits= c(-4, 2), 
						scale=as.logical(par.scale))	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )
				tmp		<- file.path(dirname(F), gsub('\\.newick',paste0('_coalreg_using_TRM-FROM-RECENTtrf_RISKtrf_ETSIaoi_',par.climb,'args3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_lasso',par.lasso,'_ETSIbias',par.bias,'_ETSInoise',par.noise,'_scale',par.scale,'.rda'),basename(F)))				
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.png.runcoalreg.using.TRSTAGE.TRRISK.TRGENDER.ETSI.vanilla.BFGS3<- function(indir, par.base.pattern, par.maxNodeDepth=3, par.maxHeight=10, par.lasso=5)
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
		par.maxNodeDepth	<- Inf
		par.maxHeight		<- 10
		par.lasso			<- 5
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#		 
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
						IDPOP=paste0('ID_',sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
						GENDER=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',2),
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
				set(phi, NULL, 'MALE', phi[,factor(GENDER,levels=c('F','M'))])
				set(phi, NULL, 'RISK_L', phi[,factor(RISK!='L',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'RISK_H', phi[,factor(RISK!='H',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				#	zero mean ETSI so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI', phi[, ETSI-mean(ETSI)])
				#	prepare coalreg input
				ph$tip.label	<- phi[, IDPOP]
				tmp				<- data.matrix(subset(phi, select=c(TRM_FROM_RECENT,RISK_L,RISK_H,MALE,ETSI)))
				rownames(tmp)	<- phi[, IDPOP]
				stopifnot(!any(is.na(tmp)))
				ph 				<- DatedTree( ph, setNames( phi$TIME_SEQ, ph$tip.label) ) 
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
				fit 	<- trf.lasso(ph, tmp, trf_names=c('TRM_FROM_RECENT','RISK_L','RISK_H','MALE'), aoi_names=c( 'ETSI' ), 
						maxNodeDepth=par.maxNodeDepth,
						maxHeight=par.maxHeight,
						lasso_threshold=par.lasso, 
						method = 'BFGS', lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )
				tmp		<- file.path(dirname(F), gsub('\\.newick',paste0('_coalreg_using_TRM-FROM-RECENTtrf_RISKtrf_MALEtrf_ETSIaoi_BFGSargs3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_lasso',par.lasso,'.rda'),basename(F)))				
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.png.runcoalreg.using.TRGENDER.ETSI.vanilla.BFGS3<- function(indir, par.base.pattern, par.maxNodeDepth=3, par.maxHeight=10, par.lasso=5)
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
		par.maxNodeDepth	<- Inf
		par.maxHeight		<- 10
		par.lasso			<- 5
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#		 
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
						IDPOP=paste0('ID_',sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
						GENDER=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',2),
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
				set(phi, NULL, 'MALE', phi[,factor(GENDER,levels=c('F','M'))])
				#	zero mean ETSI so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI', phi[, ETSI-mean(ETSI)])
				#	prepare coalreg input
				ph$tip.label	<- phi[, IDPOP]
				tmp				<- data.matrix(subset(phi, select=c(MALE,ETSI)))
				rownames(tmp)	<- phi[, IDPOP]
				stopifnot(!any(is.na(tmp)))
				ph 				<- DatedTree( ph, setNames( phi$TIME_SEQ, ph$tip.label) ) 
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
				fit 	<- trf.lasso(ph, tmp, trf_names=c('MALE'), aoi_names=c( 'ETSI' ), 
						maxNodeDepth=par.maxNodeDepth,
						maxHeight=par.maxHeight,
						lasso_threshold=par.lasso, 
						method = 'BFGS', lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )
				tmp		<- file.path(dirname(F), gsub('\\.newick',paste0('_coalreg_using_MALEtrf_ETSIaoi_BFGSargs3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_lasso',par.lasso,'.rda'),basename(F)))				
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.png.runcoalreg.using.TRRISK.ETSI.vanilla.BFGS3<- function(indir, par.base.pattern, par.maxNodeDepth=3, par.maxHeight=10, par.lasso=5)
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
		par.maxNodeDepth	<- Inf
		par.maxHeight		<- 10
		par.lasso			<- 5
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#		 
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
						IDPOP=paste0('ID_',sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
						GENDER=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',2),
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
				set(phi, NULL, 'RISK_L', phi[,factor(RISK!='L',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'RISK_H', phi[,factor(RISK!='H',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				#	zero mean ETSI so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI', phi[, ETSI-mean(ETSI)])
				#	prepare coalreg input
				ph$tip.label	<- phi[, IDPOP]
				tmp				<- data.matrix(subset(phi, select=c(RISK_L,RISK_H,ETSI)))
				rownames(tmp)	<- phi[, IDPOP]
				stopifnot(!any(is.na(tmp)))
				ph 				<- DatedTree( ph, setNames( phi$TIME_SEQ, ph$tip.label) ) 
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
				fit 	<- trf.lasso(ph, tmp, trf_names=c('RISK_L','RISK_H'), aoi_names=c( 'ETSI' ), 
						maxNodeDepth=par.maxNodeDepth,
						maxHeight=par.maxHeight,
						lasso_threshold=par.lasso, 
						method = 'BFGS', lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )
				tmp		<- file.path(dirname(F), gsub('\\.newick',paste0('_coalreg_using_RISKtrf_ETSIaoi_BFGSargs3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_lasso',par.lasso,'.rda'),basename(F)))				
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.png.runcoalreg.using.TRRISK.ETSI.unmodelled.het.BFGS3<- function(indir, par.base.pattern, par.maxNodeDepth=3, par.maxHeight=10, par.lasso=5, par.climb='BFGS', par.bias=1, par.noise=0, par.scale=0)
{
	require(coalreguh)
	require(data.table)
	require(ape)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/png_results'
		par.base.pattern	<- 'PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43'
		par.maxNodeDepth	<- Inf
		par.maxHeight		<- 10
		par.lasso			<- 5
	}
	stopifnot(par.bias==1, par.noise==0)
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#		 
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
						IDPOP=paste0('ID_',sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
						GENDER=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',2),
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
				set(phi, NULL, 'RISK_L', phi[,factor(RISK!='L',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'RISK_H', phi[,factor(RISK!='H',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				#	zero mean ETSI so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI', phi[, ETSI-mean(ETSI)])
				#	prepare coalreg input
				ph$tip.label	<- phi[, IDPOP]
				tmp				<- data.matrix(subset(phi, select=c(RISK_L,RISK_H,ETSI)))
				rownames(tmp)	<- phi[, IDPOP]
				stopifnot(!any(is.na(tmp)))
				ph 				<- DatedTree( ph, setNames( phi$TIME_SEQ, ph$tip.label) ) 
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
				fit 	<- trf.lasso(ph, tmp, trf_names=c('RISK_L','RISK_H'), aoi_names=c( 'ETSI' ), 
						maxNodeDepth=par.maxNodeDepth,
						maxHeight=par.maxHeight, 
						lasso_threshold=par.lasso, 
						method=unname(par.climb), 
						lnr0 = -2, lnrLimits = c(-4, 2), 
						scale=as.logical(par.scale))	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )				
				tmp		<- file.path(dirname(F), gsub('\\.newick',paste0('_coalreguh_using_RISKtrf_ETSIaoi_',par.climb,'args3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_lasso',par.lasso,'_ETSIbias',par.bias,'_ETSInoise',par.noise,'_scale',par.scale,'.rda'),basename(F)))				
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.png.runcoalreg.using.TRSTAGE.TRRISK.TRGENDER.TRAGEDIAG.ETSI.vanilla.BFGS3<- function(indir, par.base.pattern, par.maxNodeDepth=3, par.maxHeight=10, par.lasso=5)
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
		par.maxNodeDepth	<- Inf
		par.maxHeight		<- 10
		par.lasso			<- 5
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#		 
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
						IDPOP=paste0('ID_',sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
						GENDER=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',2),
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
				set(phi, NULL, 'AGE_AT_DIAG', phi[, cut(DIAG_T-DOB, breaks=c(-1,25,30,1e4),labels=c('less25','25to29','30'))])								
				set(phi, NULL, 'AGE_AT_DIAG_25to29', phi[,factor(is.na(AGE_AT_DIAG) | AGE_AT_DIAG!='25to29',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'AGE_AT_DIAG_less25', phi[,factor(is.na(AGE_AT_DIAG) | AGE_AT_DIAG!='less25',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'AGE_AT_DIAG_unknown', phi[,factor(!is.na(AGE_AT_DIAG),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'TRM_FROM_RECENT', phi[,factor(TRM_FROM_RECENT,levels=c('N','Y'))])
				set(phi, NULL, 'MALE', phi[,factor(GENDER,levels=c('F','M'))])
				set(phi, NULL, 'RISK_L', phi[,factor(RISK!='L',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'RISK_H', phi[,factor(RISK!='H',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				#	zero mean ETSI so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI', phi[, ETSI-mean(ETSI)])
				#	prepare coalreg input
				ph$tip.label	<- phi[, IDPOP]
				tmp				<- data.matrix(subset(phi, select=c(TRM_FROM_RECENT,RISK_L,RISK_H,MALE,AGE_AT_DIAG_25to29,AGE_AT_DIAG_less25,AGE_AT_DIAG_unknown,ETSI)))
				rownames(tmp)	<- phi[, IDPOP]
				stopifnot(!any(is.na(tmp)))
				ph 				<- DatedTree( ph, setNames( phi$TIME_SEQ, ph$tip.label) ) 
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
				fit 	<- trf.lasso(ph, tmp, trf_names=c('TRM_FROM_RECENT','RISK_L','RISK_H','MALE','AGE_AT_DIAG_25to29', 'AGE_AT_DIAG_less25', 'AGE_AT_DIAG_unknown'), aoi_names=c( 'ETSI' ), 
						maxNodeDepth=par.maxNodeDepth,
						maxHeight=par.maxHeight,
						lasso_threshold=par.lasso, 
						method = 'BFGS', lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )
				tmp		<- file.path(dirname(F), gsub('\\.newick',paste0('_coalreg_using_TRM-FROM-RECENTtrf_RISKtrf_MALEtrf_AGEATDIAGtrf_ETSIaoi_BFGSargs3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_lasso',par.lasso,'.rda'),basename(F)))				
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.png.runcoalreg.using.TRSTAGE.TRRISK.TRGENDER.TRAGEDIAG.ETSI.noise.BFGS3<- function(indir, par.base.pattern, par.maxNodeDepth=3, par.maxHeight=10, par.lasso=5, par.noise=1, par.bias=1, par.climb='BFGS', par.scale=0)
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
		par.maxNodeDepth	<- Inf
		par.maxHeight		<- 10
		par.lasso			<- 5
		par.noise			<- 1
		par.scale			<- 0
		par.climb			<- 'BFGS'
		par.bias			<- 1
	}
	#
	#	run coalreg	run using exact time to infection
	#	with extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#		 
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
						TAXA_ID=seq_along(ph$tip.label),
						IDPOP=paste0('ID_',sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',1)),
						GENDER=sapply(strsplit(ph$tip.label,'|',fixed=TRUE),'[[',2),
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
				set(phi, NULL, 'AGE_AT_DIAG', phi[, cut(DIAG_T-DOB, breaks=c(-1,25,30,1e4),labels=c('less25','25to29','30'))])								
				set(phi, NULL, 'AGE_AT_DIAG_25to29', phi[,factor(is.na(AGE_AT_DIAG) | AGE_AT_DIAG!='25to29',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'AGE_AT_DIAG_less25', phi[,factor(is.na(AGE_AT_DIAG) | AGE_AT_DIAG!='less25',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'AGE_AT_DIAG_unknown', phi[,factor(!is.na(AGE_AT_DIAG),levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'TRM_FROM_RECENT', phi[,factor(TRM_FROM_RECENT,levels=c('N','Y'))])
				set(phi, NULL, 'MALE', phi[,factor(GENDER,levels=c('F','M'))])
				set(phi, NULL, 'RISK_L', phi[,factor(RISK!='L',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				set(phi, NULL, 'RISK_H', phi[,factor(RISK!='H',levels=c(TRUE,FALSE),labels=c('N','Y'))])
				#	add noise and bias
				tmp		<- phi[, list(ETSI_NOISE=rlnorm(1, meanlog= log(ETSI)-0.5*log(par.noise*par.noise+1), sdlog= sqrt(log(par.noise*par.noise+1)))), by='TAXA']				
				phi		<- merge(phi, tmp, by='TAXA')					
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE*par.bias])
				#	zero mean ETSI_NOISE so the coefficients can be interpreted as RR
				set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE-mean(ETSI_NOISE)])								
				#	prepare coalreg input
				setkey(phi, TAXA_ID)
				ph$tip.label	<- phi[, IDPOP]
				tmp				<- data.matrix(subset(phi, select=c(TRM_FROM_RECENT,RISK_L,RISK_H,MALE,AGE_AT_DIAG_25to29,AGE_AT_DIAG_less25,AGE_AT_DIAG_unknown,ETSI_NOISE)))
				rownames(tmp)	<- phi[, IDPOP]
				stopifnot(!any(is.na(tmp)))
				ph 				<- DatedTree( ph, setNames( phi$TIME_SEQ, ph$tip.label) ) 
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
				fit 	<- trf.lasso(ph, tmp, trf_names=c('TRM_FROM_RECENT','RISK_L','RISK_H','MALE','AGE_AT_DIAG_25to29', 'AGE_AT_DIAG_less25', 'AGE_AT_DIAG_unknown'), aoi_names=c( 'ETSI_NOISE' ), 
						maxNodeDepth=par.maxNodeDepth,
						maxHeight=par.maxHeight,
						lasso_threshold=par.lasso, 
						method=unname(par.climb), 
						lnr0= -2, lnrLimits= c(-4, 2), 
						scale=as.logical(par.scale))								
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )				
				tmp		<- file.path(dirname(F), gsub('\\.newick',paste0('_coalreg_using_TRM-FROM-RECENTtrf_RISKtrf_MALEtrf_AGEATDIAGtrf_ETSIaoi_',par.climb,'args3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_lasso',par.lasso,'_ETSIbias',par.bias,'_ETSInoise',par.noise,'_scale',par.scale,'.rda'),basename(F)))				
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.master.ex3.runcoalreg.using.TYPE.ETFI.vanilla.BFGS3<- function(indir, par.base.pattern, par.s, par.maxNodeDepth, par.maxHeight)
{
	require(coalreg)
	require(viridis)
	require(data.table)
	if(0)
	{
		indir				<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_coalregression/master_examples'			
		par.base.pattern	<- 'm3.RR5.n150_seed123'	
		par.maxNodeDepth	<- 3
		par.maxHeight		<- 10
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
				fit 	<- trf.lasso(	dph, 	tmp, trf_names = c( 'TYPE'), aoi_names = c( 'ETSI' ), 
										maxNodeDepth=par.maxNodeDepth, maxHeight=par.maxHeight, 
										lasso_threshold=5, method = 'BFGS', lnr0 = -2, lnrLimits = c(-4, 2), scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  ) 
				#print(fit$bestfit$par )
				#print( fci$ci )	
				tmp		<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_s',par.s*100,'.rda'),basename(F)))
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
