cr.various<- function()
{
	cr.various.master()
	#cr.various.pangea()
}

cr.hpc.submit.170803<- function()
{
	if(0)
	{
		indir						<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
		indir						<- '/work/or105/ATHENA_2016/master_examples'	
		par.base.pattern			<- 'm3.RR5.n1250_seed123'
		infiles				<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))	
		infiles[, REP:= as.numeric(gsub('.*rep([0-9]+)\\.nwk','\\1',basename(F)))]
		infiles				<- subset(infiles, REP<=30)
		setkey(infiles, REP)
		
		formula.tr					<- '~TYPE'
		formula.inf					<- '~TYPE'		
		par.maxNodeDepth			<- Inf			
		par.hetInflation_logprior	<- 0
		par.noise					<- 0
		par.bias					<- 1		
		par.mincladesize			<- 100		
		for(i in seq_len(nrow(infiles)))
		{
			infile						<- infiles[i,F]
			par.s						<- 0.5
			par.maxHeight				<- 10
			extra						<- '_170803'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q='pqeelab', hpc.walltime=71, hpc.mem="5600mb", hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)
			
			par.s						<- 1
			par.maxHeight				<- 10
			extra						<- '_170803'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q='pqeelab', hpc.walltime=71, hpc.mem="5600mb", hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)
			
			par.s						<- 0.5
			par.maxHeight				<- 25
			extra						<- '_170713'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q='pqeelab', hpc.walltime=71, hpc.mem="5600mb", hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
			
			par.s						<- 1
			par.maxHeight				<- 25
			extra						<- '_170713'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q='pqeelab', hpc.walltime=71, hpc.mem="5600mb", hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
		}
	}
}

cr.hpc.submit.170917<- function()
{
	if(1)
	{		
		hpc.q						<- 'pqeelab'
		hpc.mem						<- '5600mb'
		hpc.walltime				<- 71
				
		cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS')
		cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
		cat(cmd)	
		outdir		<- paste(DATA,"tmp",sep='/')
		outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		hivc.cmd.hpccaller(outdir, outfile, cmd)			
	}
}

cr.hpc.submit.170914<- function()
{	
	if(1)
	{
		indir						<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
		indir						<- '/work/or105/ATHENA_2016/master_examples'	
		par.base.pattern			<- 'm3.RR5.n1250_seed123'
		infiles				<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))	
		infiles[, REP:= as.numeric(gsub('.*rep([0-9]+)\\.nwk','\\1',basename(F)))]
		infiles				<- subset(infiles, REP<=30)
		setkey(infiles, REP)
		
		hpc.q						<- 'pqeelab'
		hpc.mem						<- '5600mb'
		hpc.walltime				<- 71
		formula.tr					<- '~TYPE'
		formula.inf					<- '~TYPE'
		formula.inf					<- '~ETSI'
		par.maxNodeDepth			<- Inf		
		par.maxHeight				<- 10
		par.hetInflation_logprior	<- 0
		par.noise					<- 0
		par.bias					<- 1					
		for(i in seq_len(nrow(infiles)))
		{
			infile						<- infiles[i,F]
			par.s						<- 1
			par.mincladesize			<- 50
			extra						<- '170913'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
			
			
			infile						<- infiles[i,F]
			par.s						<- 1
			par.mincladesize			<- 75
			extra						<- '170913'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
			
			
			par.s						<- 1
			par.mincladesize			<- 100
			extra						<- '170913'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)
			
			par.s						<- 1
			par.mincladesize			<- 150
			extra						<- '170913'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)
			
			par.s						<- 1
			par.mincladesize			<- 200
			extra						<- '170913'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)							
		}
	}
}

cr.hpc.submit<- function()
{	
	if(1)
	{
		indir						<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
		indir						<- '/work/or105/ATHENA_2016/master_examples'	
		par.base.pattern			<- 'm3.RR5.n1250_seed123'
		infiles				<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))	
		infiles[, REP:= as.numeric(gsub('.*rep([0-9]+)\\.nwk','\\1',basename(F)))]
		infiles				<- subset(infiles, REP<=30)
		setkey(infiles, REP)
		
		hpc.q						<- 'pqeelab'
		hpc.mem						<- '5600mb'
		hpc.walltime				<- 71
		formula.tr					<- '~TYPE'	
		formula.inf					<- '~TYPE'
		#formula.inf					<- '~ETSI'
		#formula.inf					<- '~ETSI+TYPE'
		extra						<- '170919'
		par.maxNodeDepth			<- Inf				
		par.hetInflation_logprior	<- 1 	# cannot pass function in cmd, define within
		par.hetInflation_logprior	<- 0 	
		par.noise					<- 0
		par.bias					<- 1	
		#par.infprior.mean			<- -2
		par.infprior.mean			<- log(5)
		for(i in seq_len(nrow(infiles)))
		{
			infile						<- infiles[i,F]
			par.s						<- 1
			par.mincladesize			<- 75	
			par.maxHeight				<- 10
			par.infprior.sd				<- 1/5
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior,
					' -par.infprior.mean=',par.infprior.mean,
					' -par.infprior.sd=',par.infprior.sd,
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
			
			infile						<- infiles[i,F]
			par.s						<- 0.5
			par.mincladesize			<- 75	
			par.maxHeight				<- 10
			par.infprior.sd				<- 1/5
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior,
					' -par.infprior.mean=',par.infprior.mean,
					' -par.infprior.sd=',par.infprior.sd,					
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
			
			infile						<- infiles[i,F]
			par.s						<- 1
			par.mincladesize			<- 75	
			par.maxHeight				<- 10
			par.infprior.sd				<- 1/2
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior,
					' -par.infprior.mean=',par.infprior.mean,
					' -par.infprior.sd=',par.infprior.sd,
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
			
			infile						<- infiles[i,F]
			par.s						<- 0.5
			par.mincladesize			<- 75	
			par.maxHeight				<- 10
			par.infprior.sd				<- 1/2
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior,
					' -par.infprior.mean=',par.infprior.mean,
					' -par.infprior.sd=',par.infprior.sd,					
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
		
			infile						<- infiles[i,F]
			par.s						<- 1
			par.mincladesize			<- 75	
			par.maxHeight				<- 10
			par.infprior.sd				<- 1
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior,
					' -par.infprior.mean=',par.infprior.mean,
					' -par.infprior.sd=',par.infprior.sd,
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
			
			infile						<- infiles[i,F]
			par.s						<- 0.5
			par.mincladesize			<- 75	
			par.maxHeight				<- 10
			par.infprior.sd				<- 1
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior,
					' -par.infprior.mean=',par.infprior.mean,
					' -par.infprior.sd=',par.infprior.sd,					
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
			
		}
	}
}

cr.hpc.submit.170919<- function()
{	
	if(1)
	{
		indir						<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
		indir						<- '/work/or105/ATHENA_2016/master_examples'	
		par.base.pattern			<- 'm3.RR5.n1250_seed123'
		infiles				<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))	
		infiles[, REP:= as.numeric(gsub('.*rep([0-9]+)\\.nwk','\\1',basename(F)))]
		infiles				<- subset(infiles, REP<=30)
		setkey(infiles, REP)
		
		hpc.q						<- 'pqeelab'
		hpc.mem						<- '5600mb'
		hpc.walltime				<- 71
		formula.tr					<- '~TYPE'
		formula.inf					<- '~TYPE'
		#formula.inf					<- '~ETSI'
		extra						<- '170919'
		par.maxNodeDepth			<- Inf				
		par.hetInflation_logprior	<- 1 	# cannot pass function in cmd, define within
		par.noise					<- 0
		par.bias					<- 1					
		for(i in seq_len(nrow(infiles)))
		{
			infile						<- infiles[i,F]
			par.s						<- 1
			par.mincladesize			<- 100	
			par.maxHeight				<- 25
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
			
			infile						<- infiles[i,F]
			par.s						<- 1
			par.mincladesize			<- 50	
			par.maxHeight				<- 10
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
			
			
			infile						<- infiles[i,F]
			par.s						<- 1
			par.mincladesize			<- 75
			par.maxHeight				<- 10
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)	
			
			
			par.s						<- 1
			par.mincladesize			<- 100	
			par.maxHeight				<- 10
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)
			
			par.s						<- 1
			par.mincladesize			<- 150	
			par.maxHeight				<- 10
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)
			
			par.s						<- 1
			par.mincladesize			<- 200	
			par.maxHeight				<- 10
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=hpc.q, hpc.walltime=hpc.walltime, hpc.mem=hpc.mem, hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)							
		}
	}
}

cr.hpc.submit.170810<- function()
{	
	if(1)
	{
		indir						<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
		indir						<- '/work/or105/ATHENA_2016/master_examples'	
		par.base.pattern			<- 'm3.RR5.n1250_seed123'
		infiles				<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))	
		infiles[, REP:= as.numeric(gsub('.*rep([0-9]+)\\.nwk','\\1',basename(F)))]
		infiles				<- subset(infiles, REP<=30)
		setkey(infiles, REP)
		
		formula.tr					<- '~TYPE'
		formula.inf					<- '~ETSI'		
		par.maxNodeDepth			<- Inf	
		par.maxHeight				<- 25
		par.hetInflation_logprior	<- 0
		par.noise					<- 0
		par.bias					<- 1		
		par.mincladesize			<- 100		
		for(i in seq_len(nrow(infiles)))
		{
			infile						<- infiles[i,F]
			par.s						<- 0.5			
			extra						<- '_170810'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=NA, hpc.walltime=23, hpc.mem="1800mb", hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)
			
			par.s						<- 1
			extra						<- '_170810'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q=NA, hpc.walltime=23, hpc.mem="1800mb", hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)							
		}
	}
	if(0)
	{
		indir						<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
		indir						<- '/work/or105/ATHENA_2016/master_examples'	
		par.base.pattern			<- 'm3.RR5.n1250_seed123'
		infiles				<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))	
		infiles[, REP:= as.numeric(gsub('.*rep([0-9]+)\\.nwk','\\1',basename(F)))]
		infiles				<- subset(infiles, REP<=30)
		setkey(infiles, REP)
		
		par.maxNodeDepth			<- Inf
		par.maxHeight				<- 10	
		par.hetInflation_logprior	<- 0
		par.noise					<- 0
		par.bias					<- 1
		par.s						<- 0.5
		par.mincladesize			<- 100		
		for(i in seq_len(nrow(infiles)))
		{
			infile						<- infiles[i,F]
			formula.tr					<- '~TYPE'
			formula.inf					<- '~ETSI'
			extra						<- '_170803'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
											' -formula.tr=',formula.tr,
											' -formula.inf=',formula.inf,
											' -extra=',extra,
											' -par.maxNodeDepth=',par.maxNodeDepth,
											' -par.maxHeight=',par.maxHeight, 
											' -par.hetInflation_logprior=',par.hetInflation_logprior, 
											' -par.noise=', par.noise, 
											' -par.bias=', par.bias, 
											' -par.s=', par.s, 
											' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q='pqeelab', hpc.walltime=71, hpc.mem="5600mb", hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)
			stop()
			formula.tr					<- '~TYPE'
			formula.inf					<- '~TYPE'
			extra						<- '_170803'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q='pqeelab', hpc.walltime=71, hpc.mem="5600mb", hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)
		
			formula.tr					<- '~TYPE'
			formula.inf					<- '~ETSIC'
			extra						<- '_170713'			
			cmd		<- paste0(CODE.HOME, '/misc/hivclu.startme.R -exe=VARIOUS -infile=',infile,
					' -formula.tr=',formula.tr,
					' -formula.inf=',formula.inf,
					' -extra=',extra,
					' -par.maxNodeDepth=',par.maxNodeDepth,
					' -par.maxHeight=',par.maxHeight, 
					' -par.hetInflation_logprior=',par.hetInflation_logprior, 
					' -par.noise=', par.noise, 
					' -par.bias=', par.bias, 
					' -par.s=', par.s, 
					' -par.mincladesize=', par.mincladesize)
			cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.nproc=1, hpc.q='pqeelab', hpc.walltime=71, hpc.mem="5600mb", hpc.load='module load intel-suite R/3.3.3')
			cat(cmd)	
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("cr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			hivc.cmd.hpccaller(outdir, outfile, cmd)			
		}
	}
	if(0)
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
	#cr.various.master.MLE()
	cr.various.master.Bayes()
}

cr.various.master.Bayes<- function()
{
	infile						<- "/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep2.nwk"
	par.maxNodeDepth			<- Inf
	par.maxHeight				<- 10	
	par.hetInflation_logprior	<- NA
	par.s						<- 0.5
	par.mincladesize			<- 100
	formula.tr					<- ~TYPE
	formula.inf					<- ~ETSI
	extra						<- ''
	par.infprior.mean			<- NA
	par.infprior.sd				<- NA
	#
	#	read args
	#
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,7), infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,11), formula.tr= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) formula.tr<- as.formula(tmp[1])	
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,12), formula.inf= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) formula.inf<- as.formula(tmp[1])
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,6), extra= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) extra<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,17), par.maxNodeDepth= return(substr(arg,19,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.maxNodeDepth<- as.numeric(tmp[1])		
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,14), par.maxHeight= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.maxHeight<- as.numeric(tmp[1])
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,10), par.noise= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.noise<- as.numeric(tmp[1])
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,9), par.bias= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.bias<- as.numeric(tmp[1])		
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,26), par.hetInflation_logprior= return(substr(arg,28,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.hetInflation_logprior<- as.numeric(tmp[1])
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,6), par.s= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.s<- as.numeric(tmp[1])
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,17), par.mincladesize= return(substr(arg,19,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.mincladesize<- as.numeric(tmp[1])		
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,18), par.infprior.mean= return(substr(arg,20,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.infprior.mean<- as.numeric(tmp[1])
		tmp<- na.omit(sapply(argv,function(arg){	switch(substr(arg,2,16), par.infprior.sd= return(substr(arg,18,nchar(arg))),NA)	}))
		if(length(tmp)>0) par.infprior.sd<- as.numeric(tmp[1])		
	}	
	#	process args
	par.tsimn	<- par.noise
	par.tsimb	<- par.bias		
	#	
	cr.master.ex3.adMCMC(infile, formula.tr, formula.inf, par.s, par.maxNodeDepth, par.maxHeight, par.mincladesize, par.tsimb, par.tsimn, par.hetInflation_logprior, par.infprior.mean, par.infprior.sd, extra)		
}


cr.various.master.MLE<- function()
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
		par.s				<- 1	
		par.maxNodeDepth	<- Inf 
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
	if(0)
	{
		par.s				<- 0.5
		par.maxNodeDepth	<- Inf
		par.maxHeight		<- 10
		cr.master.ex3.runcoalreg.using.TYPE.ETFI.vanilla.MCMC1(indir, par.base.pattern, par.s, par.maxNodeDepth, par.maxHeight=par.maxHeight)	
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
	indir	<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
	infiles	<- data.table(FR=list.files(indir, pattern='_INTERNAL.R', full.names=TRUE, recursive=TRUE))
	infiles[, SC:= basename(dirname(dirname(FR)))]
	infiles[, SEQCOV:= as.numeric(gsub('.*cov([0-9]+\\.[0-9]+).*','\\1',SC))/100]
	infiles[, SEED:= as.numeric(gsub('.*seed([0-9]+).*','\\1',SC))]	
	tmp		<- data.table(FT=list.files(indir, pattern='_DATEDTREE.newick', full.names=TRUE, recursive=TRUE))
	tmp[, SC:= basename(dirname(dirname(FT)))]
	infiles	<- merge(infiles, tmp, by='SC', all=1)
	
	infiles[, table(SEQCOV)]
	
	outfile	<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-'
	t.start	<- 0
	t.end	<- 2020
	t.adult	<- 14
	t.early	<- 1
	t.acute	<- 3/12	
	#
	#	prepare tree files
	#
	infiles[, {
				#FR		<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43/150129_HPTN071_scHN_INTERNAL/150129_HPTN071_scHN_SIMULATED_INTERNAL.R'
				#FT		<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43/150129_HPTN071_scHN_SIMULATED_TREE/150129_HPTN071_scHN_DATEDTREE.newick'
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
	
	FR		<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed50/150129_HPTN071_scHN_INTERNAL/150129_HPTN071_scHN_SIMULATED_INTERNAL.R'
	#FR		<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov41.3-seed46/150129_HPTN071_scHN_INTERNAL/150129_HPTN071_scHN_SIMULATED_INTERNAL.R'
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

cr.master.ex3.runcoalreg.using.TYPE.BFGS2<- function(indir, par.base.pattern, par.s)
{
	require(coalreg)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
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
				#F		<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep1.nwk'
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

cr.png.compare.170914<- function()
{
	require(coalreg)
	require(viridis)
	indir				<- '~/Box Sync/OR_Work/2017/2017_coalregression/master_results_5'
	infiles	<- data.table(F=list.files(indir, pattern='rda$',full.names=TRUE))
	infiles	<- subset(infiles, grepl('maxNodeDepth',F))
	infiles	<- subset(infiles, !grepl('results\\.rda',F))	
	#	parse results	
	
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43_coalreg_using_TRM-FROM-RECENTtrf_ETSIaoi_BFGSargs3_maxNodeDepth3.rda'
				load(F)
				list(	FACTOR= names(fit$bestfit$par),
						MLE= unname(fit$bestfit$par),
						MLE_CONVERGED= fit$bestfit$convergence==0 )				
			}, by='F']
	res[, SEED:= as.integer(gsub('.*seed([0-9]+)_.*','\\1',F))]
	res[, REP:= as.numeric(gsub('.*_rep([0-9]+)_.*','\\1',basename(F)))]
	res[, MAXNODE:=as.numeric(gsub('.*_maxNodeDepth([0-9A-Za-z]+).*','\\1',F))]
	res[, MAXHEIGHT:=as.numeric(gsub('.*_maxHeight([0-9A-Za-z]+).*','\\1',F))]
	res[, AOI_ETSI:= as.numeric(grepl('ETFIaoi',F))]
	res[, SAMPLING:= as.numeric(gsub('.*_s([0-9]+).rda','\\1',basename(F)))]
	#
	res	<- subset(res, MLE_CONVERGED)
	#
	#	load(file.path(indir,'results.rda'))
	dfm	<- dcast.data.table(res, SAMPLING+MAXNODE+MAXHEIGHT+REP~FACTOR, value.var='MLE')
	set(res, NULL, 'SAMPLING', res[, paste0('sampling=',SAMPLING,'%')])
	set(res, NULL, 'MAXNODE', res[, paste0('max node=',MAXNODE)])
	set(res, NULL, 'MAXHEIGHT', res[, paste0('max height=',MAXHEIGHT)])
	tmp	<- unique(subset(res, FACTOR=='b_TYPE', select=c('SAMPLING','MAXNODE','MAXHEIGHT','FACTOR')))
	tmp[, TRUTH:=log(5)]	
	ggplot(res, aes(x=MAXHEIGHT)) + 
			geom_boxplot(aes(y=MLE)) +
			geom_hline(data=tmp, aes(yintercept=TRUTH), colour='red') +
			facet_grid(FACTOR~SAMPLING+MAXNODE, scales='free')
	ggsave(file=file.path(indir,'boxplots_trTYPE_infETSI_sampling100_maxheight10_vary_maxcladesize.pdf'),w=10,h=10)	
}

cr.png.compare<- function()
{
	require(coalreg)
	require(viridis)
	indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
	load("~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-empirical_transmission_rates.rda")
	dtrr	<- subset(dtrr, TS==2005 & TE==2020)
	set(dtrr, NULL, c('TS','TE'), NULL)
	#	parse results	
	infiles	<- subset(infiles, !grepl('results\\.rda',F))
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43_coalreg_using_TRM-FROM-RECENTtrf_ETSIaoi_BFGSargs3_maxNodeDepth3.rda'
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
	load('/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed51_coalreg_using_TRM-FROM-RECENTtrf_RISKtrf_ETSIaoi_BFGSargs3_maxNodeDepthInf_maxHeight10_lasso5_ETSIbias1_ETSInoise0_scale0.rda')
	fit.bfgs	<- fit
	load('~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/tmp/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed51_coalreg_using_TRM-FROM-RECENTtrf_RISKtrf_ETSIaoi_BFGSargs3_maxNodeDepthInf_maxHeight10_lasso5_ETSIbias1_ETSInoise0_scale0.rda')
	fit.bfgs.n	<- fit
	load('/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed51_coalreg_using_TRM-FROM-RECENTtrf_RISKtrf_ETSIaoi_Nelder-Meadargs3_maxNodeDepthInf_maxHeight10_lasso5_ETSIbias1_ETSInoise0_scale0.rda')
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
	
cr.master.ex3.adMCMC.evaluate.170803<- function()
{
	require(data.table)
	require(coda)
	
	indir	<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_results_2'
	infiles	<- data.table(F=list.files(indir, pattern='rda$',full.names=TRUE))
	infiles[, REP:= as.numeric(gsub('.*_rep([0-9]+)_.*','\\1',basename(F)))]
	infiles[, FORMULA_INF:= gsub('.*_inf([A-Z]+)_.*','\\1',basename(F))]
	infiles[, FORMULA_TR:= gsub('.*_tr([A-Z]+)_.*','\\1',basename(F))]
	infiles[, SAMPLING:= 100*as.numeric(gsub('.*_s([0-1]\\.?[0-9]*).*','\\1',basename(F)))]
	infiles[, MAXHEIGHT:= as.numeric(gsub('.*_mh([0-9]+)_.*','\\1',basename(F)))]
	
	tmp		<- lapply(seq_len(nrow(infiles)), function(i)
			{
				infile	<- infiles[i,F]
				cat(basename(infile),'\n')
				#infile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_results_2/m3.RR5.n1250_seed123_rep1_aMCMC170803_170713_trTYPE_infTYPE_mndInf_mh25_hetinf0_mcs100_s0.5.rda'
				load(infile)	
				fit.mcmc			<- cbind(fit$trace, fit$trace_tr, fit$trace_inf)				
				fit.mcmc			<- fit.mcmc[,-which(colnames(fit.mcmc)=='logHetInflation')]
				colnames(fit.mcmc)	<- paste0(c('','tr','inf'),colnames(fit.mcmc))				
				fit.mcmc			<- mcmc(fit.mcmc)
				#	raw plots
				pdf(file.path(dirname(infile), gsub('\\.rda',paste0('_coeff_trace.pdf'),basename(infile))), w=10, h=7)
				plot(fit.mcmc)
				dev.off()
				fit.mcmc.a	<- mcmc(fit$trace_qsd)
				pdf(file.path(dirname(infile), gsub('\\.rda',paste0('_sdadapt.pdf'),basename(infile))), w=5, h=10)
				plot(fit.mcmc.a)
				dev.off()
				#	rm burn in and thin
				tmp <- seq.int(1e3,1e4,4*5)
				dfm	<- melt(cbind(as.data.table(fit.mcmc[tmp,]), data.table(IT=tmp)), id.vars='IT')
				dfm	<- dfm[, list(STAT=paste0('q',c(0.025,0.25,0.5,0.75,0.975)), V=quantile(value,p=c(0.025,0.25,0.5,0.75,0.975))), by=c('variable')]
				#	add acceptance
				tmp	<- 1-rejectionRate(window(fit.mcmc, start=1e3, end=1e4, thin=4))
				tmp	<- data.table(variable=names(tmp), STAT='acceptance',V=tmp)
				dfm	<- rbind(dfm, tmp)
				#	add ESS
				tmp	<- effectiveSize(window(fit.mcmc, start=1e3, end=1e4, thin=4))
				tmp	<- data.table(variable=names(tmp), STAT='ESS',V=tmp)
				dfm	<- rbind(dfm, tmp)
				#
				#fit.mcmc	<- as.data.table(fit.mcmc)
				#fit.mcmc[, IT:= seq_len(nrow(fit.mcmc))]
				#fit.mcmc	<- subset(fit.mcmc, IT%%20==1)
				#ggplot(fit.mcmc, aes(y=trTYPE, x=infETSI)) + geom_point() + geom_hline(yintercept=log(5), colour='red')				
				#
				dfm[, REP:=infiles[i,REP]]
				dfm[, FORMULA_INF:=infiles[i,FORMULA_INF]]
				dfm[, FORMULA_TR:=infiles[i,FORMULA_TR]]
				dfm[, SAMPLING:=infiles[i,SAMPLING]]
				dfm[, MAXHEIGHT:=infiles[i,MAXHEIGHT]]
				dfm
			})
	dfm	<- do.call('rbind',tmp)
	save(dfm, file=file.path(indir,'results.rda'))
	
	dfm	<- dcast.data.table(dfm, SAMPLING+MAXHEIGHT+REP+variable~STAT, value.var='V')
	set(dfm, NULL, 'SAMPLING', dfm[, paste0('sampling=',SAMPLING,'%')])
	set(dfm, NULL, 'MAXHEIGHT', dfm[, paste0('max height=',MAXHEIGHT)])
	tmp	<- unique(subset(dfm, variable=='trTYPE', select=c('SAMPLING','MAXHEIGHT','REP','variable')))
	tmp[, TRUTH:=log(5)]	
	ggplot(dfm, aes(x=REP)) + 
			geom_boxplot(aes(middle=q0.5, lower=q0.25, upper=q0.75, ymin=q0.025, ymax=q0.975),stat = "identity") +
			geom_line(data=tmp, aes(x=REP, y=TRUTH), colour='red') +
			facet_grid(variable~SAMPLING+MAXHEIGHT, scales='free')
	ggsave(file=file.path(indir,'boxplots.pdf'),w=10,h=10)
}

cr.master.ex3.adMCMC.evaluate.170913<- function()
{
	require(data.table)
	require(coda)
	require(ggplot2)
	
	indir	<- '~/Box Sync/OR_Work/2017/2017_coalregression/master_results_3'
	infiles	<- data.table(F=list.files(indir, pattern='^m3.*rda$',full.names=TRUE))
	infiles[, REP:= as.numeric(gsub('.*_rep([0-9]+)_.*','\\1',basename(F)))]
	infiles[, FORMULA_INF:= gsub('.*_inf([A-Z]+)_.*','\\1',basename(F))]
	infiles[, FORMULA_TR:= gsub('.*_tr([A-Z]+)_.*','\\1',basename(F))]
	infiles[, SAMPLING:= 100*as.numeric(gsub('.*_s([0-1]\\.?[0-9]*).*','\\1',basename(F)))]
	infiles[, MAXHEIGHT:= as.numeric(gsub('.*_mh([0-9]+)_.*','\\1',basename(F)))]
	
	tmp		<- lapply(seq_len(nrow(infiles)), function(i)
			{
				infile	<- infiles[i,F]
				cat(basename(infile),'\n')
				#infile	<- '/Users/Oliver/Box Sync/OR_Work/2017/2017_coalregression/master_results_3/m3.RR5.n1250_seed123_rep9_aMCMC170803_170810_trTYPE_infETSI_mndInf_mh25_hetinf0_mcs100_s1.rda'
				load(infile)	
				fit.mcmc			<- cbind(fit$trace, fit$trace_tr, fit$trace_inf)				
				fit.mcmc			<- fit.mcmc[,-which(colnames(fit.mcmc)=='logHetInflation')]
				colnames(fit.mcmc)	<- paste0(c('','tr','inf'),colnames(fit.mcmc))				
				fit.mcmc			<- mcmc(fit.mcmc)
				#	raw plots
				pdf(file.path(dirname(infile), gsub('\\.rda',paste0('_coeff_trace.pdf'),basename(infile))), w=10, h=7)
				plot(fit.mcmc)
				dev.off()
				fit.mcmc.a	<- mcmc(fit$trace_qsd)
				pdf(file.path(dirname(infile), gsub('\\.rda',paste0('_sdadapt.pdf'),basename(infile))), w=5, h=10)
				plot(fit.mcmc.a)
				dev.off()
				#	rm burn in and thin
				tmp <- seq.int(1e3,1e4,4*5)
				dfm	<- melt(cbind(as.data.table(fit.mcmc[tmp,]), data.table(IT=tmp)), id.vars='IT')
				dfm	<- dfm[, list(STAT=paste0('q',c(0.025,0.25,0.5,0.75,0.975)), V=quantile(value,p=c(0.025,0.25,0.5,0.75,0.975))), by=c('variable')]
				#	add acceptance
				tmp	<- 1-rejectionRate(window(fit.mcmc, start=1e3, end=1e4, thin=4))
				tmp	<- data.table(variable=names(tmp), STAT='acceptance',V=tmp)
				dfm	<- rbind(dfm, tmp)
				#	add ESS
				tmp	<- effectiveSize(window(fit.mcmc, start=1e3, end=1e4, thin=4))
				tmp	<- data.table(variable=names(tmp), STAT='ESS',V=tmp)
				dfm	<- rbind(dfm, tmp)
				#
				dfm[, REP:=infiles[i,REP]]
				dfm[, FORMULA_INF:=infiles[i,FORMULA_INF]]
				dfm[, FORMULA_TR:=infiles[i,FORMULA_TR]]
				dfm[, SAMPLING:=infiles[i,SAMPLING]]
				dfm[, MAXHEIGHT:=infiles[i,MAXHEIGHT]]
				dfm
			})
	dfm	<- do.call('rbind',tmp)
	save(dfm, file=file.path(indir,'results.rda'))
	
	dfm	<- dcast.data.table(dfm, SAMPLING+MAXHEIGHT+REP+variable~STAT, value.var='V')
	set(dfm, NULL, 'SAMPLING', dfm[, paste0('sampling=',SAMPLING,'%')])
	set(dfm, NULL, 'MAXHEIGHT', dfm[, paste0('max height=',MAXHEIGHT)])
	tmp	<- unique(subset(dfm, variable=='trTYPE', select=c('SAMPLING','MAXHEIGHT','REP','variable')))
	tmp[, TRUTH:=log(5)]	
	ggplot(dfm, aes(x=REP)) + 
			geom_boxplot(aes(middle=q0.5, lower=q0.25, upper=q0.75, ymin=q0.025, ymax=q0.975),stat = "identity") +
			geom_line(data=tmp, aes(x=REP, y=TRUTH), colour='red') +
			facet_grid(variable~SAMPLING+MAXHEIGHT, scales='free')
	ggsave(file=file.path(indir,'boxplots_trTYPE_infETSI.pdf'),w=10,h=10)
	
	
	tmp		<- lapply(seq_len(nrow(infiles)), function(i)
			{
				infile	<- infiles[i,F]				
				cat(basename(infile),'\n')
				#infile	<- '/Users/Oliver/Box Sync/OR_Work/2017/2017_coalregression/master_results_3/m3.RR5.n1250_seed123_rep9_aMCMC170803_170810_trTYPE_infETSI_mndInf_mh25_hetinf0_mcs100_s1.rda'
				load(infile)	
				fit.mcmc			<- cbind(IT=seq_along(fit$loglik), fit$trace, fit$trace_tr, fit$trace_inf, LLKL=fit$loglik)				
				fit.mcmc			<- fit.mcmc[,-which(colnames(fit.mcmc)=='logHetInflation')]
				colnames(fit.mcmc)	<- paste0(c('','','tr','inf',''),colnames(fit.mcmc))
				fit.mcmc			<- subset(as.data.table(fit.mcmc), IT>=1e3)				
				fit.mcmc			<- unique(fit.mcmc, by=setdiff(colnames(fit.mcmc),'IT'))
				fit.mcmc			<- melt(fit.mcmc, id.vars=c('IT','logscale','trTYPE','LLKL'))
				fit.mcmc[, REP:=infiles[i,REP]]
				fit.mcmc[, SAMPLING:=infiles[i,SAMPLING]]
				fit.mcmc[, MAXHEIGHT:=infiles[i,MAXHEIGHT]]
				fit.mcmc
			})
	dfl	<- do.call('rbind',tmp)	
	save(dfl, file=file.path(indir,'results_lkl.rda'))
	for(x in c('infETSI'))
	{	
		tmp	<- unique(subset(dfl, SAMPLING==100 & variable==x, c(REP, trTYPE, LLKL, variable, value)), by=c('REP','trTYPE','value'))	
		ggplot(tmp, aes(x=trTYPE, y=value, z=LLKL)) + 
				stat_density2d(colour='black') +
				labs(x='H coeff for I0', y='S coeff for ETSI', colour='log likelihood') +
				geom_vline(xintercept=log(5), colour='red') + 
				theme_bw() +
				facet_wrap(~REP, ncol=6, scales='free_y')
		ggsave(file=file.path(indir,paste0('loglklsurface_with_formula.inf_',x,'_s100pc.pdf')),w=15,h=10)
	}
}

cr.master.ex3.adMCMC.evaluate.170914<- function()
{
	require(data.table)
	require(coda)
	require(ggplot2)
	
	indir	<- '~/Box Sync/OR_Work/2017/2017_coalregression/master_results_4'
	infiles	<- data.table(F=list.files(indir, pattern='^m3.*rda$',full.names=TRUE))
	infiles[, REP:= as.numeric(gsub('.*_rep([0-9]+)_.*','\\1',basename(F)))]
	infiles[, FORMULA_INF:= gsub('.*_inf([A-Z]+)_.*','\\1',basename(F))]
	infiles[, FORMULA_TR:= gsub('.*_tr([A-Z]+)_.*','\\1',basename(F))]
	infiles[, SAMPLING:= 100*as.numeric(gsub('.*_s([0-1]\\.?[0-9]*).*','\\1',basename(F)))]
	infiles[, MAXHEIGHT:= as.numeric(gsub('.*_mh([0-9]+)_.*','\\1',basename(F)))]
	infiles[, MAXCLADE:= as.numeric(gsub('.*_mcs([0-9]+)_.*','\\1',basename(F)))]
	
	tmp		<- lapply(seq_len(nrow(infiles)), function(i)
			{
				infile	<- infiles[i,F]
				cat(basename(infile),'\n')
				#infile	<- '/Users/Oliver/Box Sync/OR_Work/2017/2017_coalregression/master_results_3/m3.RR5.n1250_seed123_rep9_aMCMC170803_170810_trTYPE_infETSI_mndInf_mh25_hetinf0_mcs100_s1.rda'
				load(infile)	
				fit.mcmc			<- cbind(fit$trace, fit$trace_tr, fit$trace_inf)				
				fit.mcmc			<- fit.mcmc[,-which(colnames(fit.mcmc)=='logHetInflation')]
				colnames(fit.mcmc)	<- paste0(c('','tr','inf'),colnames(fit.mcmc))				
				fit.mcmc			<- mcmc(fit.mcmc)
				#	raw plots
				pdf(file.path(dirname(infile), gsub('\\.rda',paste0('_coeff_trace.pdf'),basename(infile))), w=10, h=7)
				plot(fit.mcmc)
				dev.off()
				fit.mcmc.a	<- mcmc(fit$trace_qsd)
				pdf(file.path(dirname(infile), gsub('\\.rda',paste0('_sdadapt.pdf'),basename(infile))), w=5, h=10)
				plot(fit.mcmc.a)
				dev.off()
				#	rm burn in and thin
				tmp <- seq.int(1e3,1e4,4*5)
				dfm	<- melt(cbind(as.data.table(fit.mcmc[tmp,]), data.table(IT=tmp)), id.vars='IT')
				dfm	<- dfm[, list(STAT=paste0('q',c(0.025,0.25,0.5,0.75,0.975)), V=quantile(value,p=c(0.025,0.25,0.5,0.75,0.975))), by=c('variable')]
				#	add acceptance
				tmp	<- 1-rejectionRate(window(fit.mcmc, start=1e3, end=1e4, thin=4))
				tmp	<- data.table(variable=names(tmp), STAT='acceptance',V=tmp)
				dfm	<- rbind(dfm, tmp)
				#	add ESS
				tmp	<- effectiveSize(window(fit.mcmc, start=1e3, end=1e4, thin=4))
				tmp	<- data.table(variable=names(tmp), STAT='ESS',V=tmp)
				dfm	<- rbind(dfm, tmp)
				#
				dfm[, REP:=infiles[i,REP]]
				dfm[, FORMULA_INF:=infiles[i,FORMULA_INF]]
				dfm[, FORMULA_TR:=infiles[i,FORMULA_TR]]
				dfm[, SAMPLING:=infiles[i,SAMPLING]]
				dfm[, MAXHEIGHT:=infiles[i,MAXHEIGHT]]
				dfm[, MAXCLADE:=infiles[i,MAXCLADE]]				
				dfm
			})
	dfm	<- do.call('rbind',tmp)
	save(dfm, file=file.path(indir,'results.rda'))
	#	load(file.path(indir,'results.rda'))
	dfm	<- dcast.data.table(dfm, SAMPLING+MAXCLADE+MAXHEIGHT+REP+variable+FORMULA_INF~STAT, value.var='V')
	set(dfm, NULL, 'SAMPLING', dfm[, paste0('sampling=',SAMPLING,'%')])
	set(dfm, NULL, 'MAXCLADE', dfm[, paste0('max clade=',MAXCLADE)])
	set(dfm, NULL, 'MAXHEIGHT', dfm[, paste0('max height=',MAXHEIGHT)])
	tmp	<- unique(subset(dfm, FORMULA_INF=='ETSI' & variable=='trTYPE', select=c('SAMPLING','MAXCLADE','MAXHEIGHT','REP','variable')))
	tmp[, TRUTH:=log(5)]	
	ggplot(subset(dfm, FORMULA_INF=='ETSI' & MAXHEIGHT=='max height=10'), aes(x=REP)) + 
			geom_boxplot(aes(middle=q0.5, lower=q0.25, upper=q0.75, ymin=q0.025, ymax=q0.975),stat = "identity") +
			geom_line(data=tmp, aes(x=REP, y=TRUTH), colour='red') +
			facet_grid(variable~SAMPLING+MAXCLADE+MAXHEIGHT, scales='free')
	ggsave(file=file.path(indir,'boxplots_trTYPE_infETSI_sampling100_maxheight10_vary_maxcladesize.pdf'),w=10,h=10)
	
	tmp	<- unique(subset(dfm, FORMULA_INF=='TYPE' & variable=='trTYPE' & MAXHEIGHT=='max height=10', select=c('SAMPLING','MAXCLADE','MAXHEIGHT','REP','variable')))
	tmp[, TRUTH:=log(5)]	
	ggplot(subset(dfm, FORMULA_INF=='TYPE' & MAXHEIGHT=='max height=10'), aes(x=REP)) + 
			geom_boxplot(aes(middle=q0.5, lower=q0.25, upper=q0.75, ymin=q0.025, ymax=q0.975),stat = "identity") +
			geom_line(data=tmp, aes(x=REP, y=TRUTH), colour='red') +
			facet_grid(variable~SAMPLING+MAXCLADE+MAXHEIGHT, scales='free')
	ggsave(file=file.path(indir,'boxplots_trTYPE_infTYPE_sampling100_maxheight10_vary_maxcladesize.pdf'),w=10,h=10)
}

cr.master.ex3.adMCMC.evaluate.170924<- function()
{
	require(data.table)
	require(coda)
	require(ggplot2)
	
	indir	<- '~/Box Sync/OR_Work/2017/2017_coalregression/master_results_6'
	infiles	<- data.table(F=list.files(indir, pattern='^m3.*rda$',full.names=TRUE))
	infiles[, REP:= as.numeric(gsub('.*_rep([0-9]+)_.*','\\1',basename(F)))]
	infiles[, FORMULA_INF:= gsub('.*_inf([A-Z]+)_.*','\\1',basename(F))]
	infiles[, FORMULA_TR:= gsub('.*_tr([A-Z]+)_.*','\\1',basename(F))]
	infiles[, SAMPLING:= 100*as.numeric(gsub('.*_s([0-1]\\.?[0-9]*).*','\\1',basename(F)))]
	infiles[, MAXHEIGHT:= as.numeric(gsub('.*_mh([0-9]+)_.*','\\1',basename(F)))]
	infiles[, MAXCLADE:= as.numeric(gsub('.*_mcs([0-9]+)_.*','\\1',basename(F)))]
	infiles[, HETINF:= as.numeric(gsub('.*_hetinf([0-9]+)_.*','\\1',basename(F)))]
	
	tmp		<- lapply(seq_len(nrow(infiles)), function(i)
			{
				infile	<- infiles[i,F]
				cat(basename(infile),'\n')
				#infile	<- '/Users/Oliver/Box Sync/OR_Work/2017/2017_coalregression/master_results_6/m3.RR5.n1250_seed123_rep1_aMCMC170919_trTYPE_infTYPE_mndInf_mh10_hetinf1_mcs100_s1.rda'
				load(infile)	
				fit.mcmc			<- cbind(fit$trace, fit$trace_tr, fit$trace_inf)				
				#fit.mcmc			<- fit.mcmc[,-which(colnames(fit.mcmc)=='logHetInflation')]
				colnames(fit.mcmc)	<- paste0(c('','','tr','inf'),colnames(fit.mcmc))				
				fit.mcmc			<- mcmc(fit.mcmc)
				#	raw plots
				pdf(file.path(dirname(infile), gsub('\\.rda',paste0('_coeff_trace.pdf'),basename(infile))), w=10, h=7)
				plot(fit.mcmc)
				dev.off()
				fit.mcmc.a	<- mcmc(fit$trace_qsd)
				pdf(file.path(dirname(infile), gsub('\\.rda',paste0('_sdadapt.pdf'),basename(infile))), w=5, h=10)
				plot(fit.mcmc.a)
				dev.off()
				#	rm burn in and thin
				tmp <- seq.int(1e3,1e4,4*5)
				dfm	<- melt(cbind(as.data.table(fit.mcmc[tmp,]), data.table(IT=tmp)), id.vars='IT')
				dfm	<- dfm[, list(STAT=paste0('q',c(0.025,0.25,0.5,0.75,0.975)), V=quantile(value,p=c(0.025,0.25,0.5,0.75,0.975))), by=c('variable')]
				#	add acceptance
				tmp	<- 1-rejectionRate(window(fit.mcmc, start=1e3, end=1e4, thin=4))
				tmp	<- data.table(variable=names(tmp), STAT='acceptance',V=tmp)
				dfm	<- rbind(dfm, tmp)
				#	add ESS
				tmp	<- effectiveSize(window(fit.mcmc, start=1e3, end=1e4, thin=4))
				tmp	<- data.table(variable=names(tmp), STAT='ESS',V=tmp)
				dfm	<- rbind(dfm, tmp)
				#
				dfm[, REP:=infiles[i,REP]]
				dfm[, FORMULA_INF:=infiles[i,FORMULA_INF]]
				dfm[, FORMULA_TR:=infiles[i,FORMULA_TR]]
				dfm[, SAMPLING:=infiles[i,SAMPLING]]
				dfm[, MAXHEIGHT:=infiles[i,MAXHEIGHT]]
				dfm[, MAXCLADE:=infiles[i,MAXCLADE]]
				dfm[, HETINF:=infiles[i,HETINF]]
				dfm
			})
	dfm	<- do.call('rbind',tmp)
	save(dfm, file=file.path(indir,'results.rda'))
	#	load(file.path(indir,'results.rda'))
	dfm	<- dcast.data.table(dfm, SAMPLING+MAXCLADE+MAXHEIGHT+REP+variable+FORMULA_INF+HETINF~STAT, value.var='V')
	set(dfm, NULL, 'SAMPLING', dfm[, paste0('sampling=',SAMPLING,'%')])
	set(dfm, NULL, 'MAXCLADE', dfm[, paste0('max clade=',MAXCLADE)])
	set(dfm, NULL, 'MAXHEIGHT', dfm[, paste0('max height=',MAXHEIGHT)])
	tmp	<- unique(subset(dfm, FORMULA_INF=='TYPE' & variable=='trTYPE', select=c('SAMPLING','MAXCLADE','MAXHEIGHT','REP','variable')))
	tmp[, TRUTH:=log(5)]	
	ggplot(subset(dfm, FORMULA_INF=='TYPE'), aes(x=REP)) + 
			geom_boxplot(aes(middle=q0.5, lower=q0.25, upper=q0.75, ymin=q0.025, ymax=q0.975),stat = "identity") +
			geom_line(data=tmp, aes(x=REP, y=TRUTH), colour='red') +
			facet_grid(variable~SAMPLING+MAXCLADE+MAXHEIGHT, scales='free')
	ggsave(file=file.path(indir,'boxplots_trTYPE_infTYPE_sampling100_maxheight10_hetinf1_vary_maxcladesize.pdf'),w=10,h=10)
	
	tmp	<- unique(subset(dfm, FORMULA_INF=='TYPE' & variable=='trTYPE' & MAXHEIGHT=='max height=10', select=c('SAMPLING','MAXCLADE','MAXHEIGHT','REP','variable')))
	tmp[, TRUTH:=log(5)]	
	ggplot(subset(dfm, FORMULA_INF=='TYPE' & MAXHEIGHT=='max height=10'), aes(x=REP)) + 
			geom_boxplot(aes(middle=q0.5, lower=q0.25, upper=q0.75, ymin=q0.025, ymax=q0.975),stat = "identity") +
			geom_line(data=tmp, aes(x=REP, y=TRUTH), colour='red') +
			facet_grid(variable~SAMPLING+MAXCLADE+MAXHEIGHT, scales='free')
	ggsave(file=file.path(indir,'boxplots_trTYPE_infTYPE_sampling100_maxheight10_vary_maxcladesize.pdf'),w=10,h=10)
}

cr.master.ex3.adMCMC.evaluate.170925<- function()
{
	require(data.table)
	require(coda)
	require(ggplot2)
	
	indir	<- '~/Box Sync/OR_Work/2017/2017_coalregression/master_results_7'
	infiles	<- data.table(F=list.files(indir, pattern='^m3.*rda$',full.names=TRUE))
	infiles[, REP:= as.numeric(gsub('.*_rep([0-9]+)_.*','\\1',basename(F)))]
	infiles[, FORMULA_INF:= gsub('.*_inf([A-Z\\+]+)_.*','\\1',basename(F))]
	infiles[, FORMULA_TR:= gsub('.*_tr([A-Z]+)_.*','\\1',basename(F))]
	infiles[, SAMPLING:= 100*as.numeric(gsub('.*_s([0-1]\\.?[0-9]*).*','\\1',basename(F)))]
	infiles[, MAXHEIGHT:= as.numeric(gsub('.*_mh([0-9]+)_.*','\\1',basename(F)))]
	infiles[, MAXCLADE:= as.numeric(gsub('.*_mcs([0-9]+)_.*','\\1',basename(F)))]
	infiles[, HETINF:= as.numeric(gsub('.*_hetinf([0-9]+)_.*','\\1',basename(F)))]	
	infiles[, INFPRIOR_S:= as.numeric(gsub('.*_infs([0-1]\\.?[0-9]*)_.*','\\1',basename(F)))]
	infiles[, INFPRIOR_M:= as.numeric(gsub('.*_infm(-?[0-1]?\\.?[0-9]*).*','\\1',basename(F)))]
	
	tmp		<- lapply(seq_len(nrow(infiles)), function(i)
			{
				infile	<- infiles[i,F]
				cat(basename(infile),'\n')
				#infile	<- '/Users/Oliver/Box Sync/OR_Work/2017/2017_coalregression/master_results_7/m3.RR5.n1250_seed123_rep9_aMCMC170919_trTYPE_infETSI+TYPE_mndInf_mh10_hetinf1_mcs75_infm-2_infs0.5_s0.5.rda'
				load(infile)					
				fit.mcmc			<- cbind(fit$trace, fit$trace_tr, fit$trace_inf)
				tmp					<- c(rep('',ncol(fit$trace)), rep('tr',ncol(fit$trace_tr)), rep('inf',ncol(fit$trace_inf)))
				colnames(fit.mcmc)	<- paste0(tmp,colnames(fit.mcmc))				
				fit.mcmc			<- mcmc(fit.mcmc)
				#	raw plots
				pdf(file.path(dirname(infile), gsub('\\.rda',paste0('_coeff_trace.pdf'),basename(infile))), w=10, h=7)
				plot(fit.mcmc)
				dev.off()
				fit.mcmc.a	<- mcmc(fit$trace_qsd)
				pdf(file.path(dirname(infile), gsub('\\.rda',paste0('_sdadapt.pdf'),basename(infile))), w=5, h=10)
				plot(fit.mcmc.a)
				dev.off()
				#	rm burn in and thin
				tmp <- seq.int(1e3,1e4,4*5)
				dfm	<- melt(cbind(as.data.table(fit.mcmc[tmp,]), data.table(IT=tmp)), id.vars='IT')
				dfm	<- dfm[, list(STAT=paste0('q',c(0.025,0.25,0.5,0.75,0.975)), V=quantile(value,p=c(0.025,0.25,0.5,0.75,0.975))), by=c('variable')]
				#	add acceptance
				tmp	<- 1-rejectionRate(window(fit.mcmc, start=1e3, end=1e4, thin=4))
				tmp	<- data.table(variable=names(tmp), STAT='acceptance',V=tmp)
				dfm	<- rbind(dfm, tmp)
				#	add ESS
				tmp	<- effectiveSize(window(fit.mcmc, start=1e3, end=1e4, thin=4))
				tmp	<- data.table(variable=names(tmp), STAT='ESS',V=tmp)
				dfm	<- rbind(dfm, tmp)
				#
				dfm[, REP:=infiles[i,REP]]
				dfm[, FORMULA_INF:=infiles[i,FORMULA_INF]]
				dfm[, FORMULA_TR:=infiles[i,FORMULA_TR]]
				dfm[, SAMPLING:=infiles[i,SAMPLING]]
				dfm[, MAXHEIGHT:=infiles[i,MAXHEIGHT]]
				dfm[, MAXCLADE:=infiles[i,MAXCLADE]]
				dfm[, HETINF:=infiles[i,HETINF]]				
				dfm[, INFPRIOR_M:=infiles[i,INFPRIOR_M]]
				dfm[, INFPRIOR_S:=infiles[i,INFPRIOR_S]]
				dfm
			})
	dfm	<- do.call('rbind',tmp)
	save(dfm, file=file.path(indir,'results.rda'))
	#	load(file.path(indir,'results.rda'))
	dfm	<- dcast.data.table(dfm, SAMPLING+MAXCLADE+MAXHEIGHT+REP+variable+FORMULA_INF+HETINF+INFPRIOR_M+INFPRIOR_S~STAT, value.var='V')
	set(dfm, NULL, 'SAMPLING', dfm[, paste0('sampling=',SAMPLING,'%')])
	set(dfm, NULL, 'MAXCLADE', dfm[, paste0('max clade=',MAXCLADE)])
	set(dfm, NULL, 'MAXHEIGHT', dfm[, paste0('max height=',MAXHEIGHT)])
	set(dfm, NULL, 'INFPRIOR_S', dfm[, paste0('inf prior sd=',INFPRIOR_S)])
	tmp	<- unique(subset(dfm, FORMULA_INF=='ETSI' & variable=='trTYPE', select=c('SAMPLING','MAXCLADE','MAXHEIGHT','REP','variable')))
	tmp[, TRUTH:=log(5)]	
	ggplot(subset(dfm, FORMULA_INF=='ETSI'), aes(x=REP)) + 
			geom_boxplot(aes(middle=q0.5, lower=q0.25, upper=q0.75, ymin=q0.025, ymax=q0.975),stat = "identity") +
			geom_line(data=tmp, aes(x=REP, y=TRUTH), colour='red') +
			facet_grid(variable~SAMPLING+MAXCLADE+MAXHEIGHT+INFPRIOR_S, scales='free')
	ggsave(file=file.path(indir,'boxplots_trTYPE_infETSI_maxheight10_hetinf1_vary_infpriorsd.pdf'),w=10,h=10)
	
	tmp	<- unique(subset(dfm, FORMULA_INF=='TYPE' & variable=='trTYPE' & MAXHEIGHT=='max height=10', select=c('SAMPLING','MAXCLADE','MAXHEIGHT','REP','variable')))
	tmp[, TRUTH:=log(5)]	
	ggplot(subset(dfm, FORMULA_INF=='TYPE'), aes(x=REP)) + 
			geom_boxplot(aes(middle=q0.5, lower=q0.25, upper=q0.75, ymin=q0.025, ymax=q0.975),stat = "identity") +
			geom_line(data=tmp, aes(x=REP, y=TRUTH), colour='red') +
			facet_grid(variable~SAMPLING+MAXCLADE+MAXHEIGHT+INFPRIOR_S, scales='free')
	ggsave(file=file.path(indir,'boxplots_trTYPE_infTYPE_maxheight10_hetinf1_vary_infpriorsd.pdf'),w=10,h=10)
	
	tmp	<- unique(subset(dfm, FORMULA_INF=='ETSI+TYPE' & variable=='trTYPE' & MAXHEIGHT=='max height=10', select=c('SAMPLING','MAXCLADE','MAXHEIGHT','REP','variable')))
	tmp[, TRUTH:=log(5)]	
	ggplot(subset(dfm, FORMULA_INF=='ETSI+TYPE'), aes(x=REP)) + 
			geom_boxplot(aes(middle=q0.5, lower=q0.25, upper=q0.75, ymin=q0.025, ymax=q0.975),stat = "identity") +
			geom_line(data=tmp, aes(x=REP, y=TRUTH), colour='red') +
			facet_grid(variable~SAMPLING+MAXCLADE+MAXHEIGHT+INFPRIOR_S, scales='free')
	ggsave(file=file.path(indir,'boxplots_trTYPE_infETSI+TYPE_maxheight10_hetinf1_vary_infpriorsd.pdf'),w=10,h=10)
	
}

cr.master.ex3.adMCMC.evaluate<- function()
{
	require(data.table)
	require(coda)
	
	indir	<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_results'
	infiles	<- data.table(F=list.files(indir, pattern='rda$',full.names=TRUE))
	infiles[, REP:= as.numeric(gsub('.*_rep([0-9]+)_.*','\\1',basename(F)))]
	infiles[, FORMULA_INF:= gsub('.*_inf([A-Z]+)_.*','\\1',basename(F))]
	tmp		<- lapply(seq_len(nrow(infiles)), function(i)
			{
				infile	<- infiles[i,F]
				cat(basename(infile),'\n')
				#infile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_results/m3.RR5.n1250_seed123_rep1_aMCMC170710_170713_trTYPE_infETSI_mndInf_mh10_hetinf0_mcs100_s0.5.rda'
				load(infile)	
				fit.mcmc			<- cbind(fit$trace, fit$trace_tr, fit$trace_inf)				
				fit.mcmc			<- fit.mcmc[,-which(colnames(fit.mcmc)=='logHetInflation')]
				colnames(fit.mcmc)	<- paste0(c('','tr','inf'),colnames(fit.mcmc))				
				fit.mcmc			<- mcmc(fit.mcmc)
				#	raw plots
				pdf(file.path(dirname(infile), gsub('\\.rda',paste0('_coeff_trace.pdf'),basename(infile))), w=10, h=7)
				plot(fit.mcmc)
				dev.off()
				fit.mcmc.a	<- mcmc(fit$trace_qsd)
				pdf(file.path(dirname(infile), gsub('\\.rda',paste0('_sdadapt.pdf'),basename(infile))), w=5, h=10)
				plot(fit.mcmc.a)
				dev.off()
				#	rm burn in and thin
				tmp <- seq.int(1e3,1e4,4*5)
				dfm	<- melt(cbind(as.data.table(fit.mcmc[tmp,]), data.table(IT=tmp)), id.vars='IT')
				dfm	<- dfm[, list(STAT=paste0('q',c(0.025,0.25,0.5,0.75,0.975)), V=quantile(value,p=c(0.025,0.25,0.5,0.75,0.975))), by=c('variable')]
				#	add acceptance
				tmp	<- 1-rejectionRate(window(fit.mcmc, start=1e3, end=1e4, thin=4))
				tmp	<- data.table(variable=names(tmp), STAT='acceptance',V=tmp)
				dfm	<- rbind(dfm, tmp)
				#	add ESS
				tmp	<- effectiveSize(window(fit.mcmc, start=1e3, end=1e4, thin=4))
				tmp	<- data.table(variable=names(tmp), STAT='ESS',V=tmp)
				dfm	<- rbind(dfm, tmp)
				#
				dfm[, REP:=infiles[i,REP]]
				dfm[, FORMULA_INF:=infiles[i,FORMULA_INF]]
				dfm
			})
	dfm	<- do.call('rbind',tmp)
	
	dfm	<- dcast.data.table(dfm, FORMULA_INF+REP+variable~STAT, value.var='V') 
	tmp	<- unique(subset(dfm, variable=='trTYPE', select=c('FORMULA_INF','REP','variable')))
	tmp[, TRUTH:=log(5)]	
	ggplot(dfm, aes(x=REP)) + 
			geom_boxplot(aes(middle=q0.5, lower=q0.25, upper=q0.75, ymin=q0.025, ymax=q0.975),stat = "identity") +
			geom_line(data=tmp, aes(x=REP, y=TRUTH), colour='red') +
			facet_grid(variable~FORMULA_INF, scales='free')
	ggsave(file=file.path(indir,'boxplots.pdf'),w=10,h=10)
	
	
	tmp		<- lapply(seq_len(nrow(infiles)), function(i)
			{
				infile	<- infiles[i,F]
				cat(basename(infile),'\n')
				#infile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_results/m3.RR5.n1250_seed123_rep1_aMCMC170710_170713_trTYPE_infTYPE_mndInf_mh10_hetinf0_mcs100_s0.5.rda'
				load(infile)	
				fit.mcmc			<- cbind(IT=seq_along(fit$loglik), fit$trace, fit$trace_tr, fit$trace_inf, LOGPO=fit$loglik)				
				fit.mcmc			<- fit.mcmc[,-which(colnames(fit.mcmc)=='logHetInflation')]
				colnames(fit.mcmc)	<- paste0(c('','','tr','inf',''),colnames(fit.mcmc))
				fit.mcmc			<- subset(as.data.table(fit.mcmc), IT>=1e3)				
				fit.mcmc			<- unique(fit.mcmc, by=setdiff(colnames(fit.mcmc),'IT'))
				fit.mcmc			<- melt(fit.mcmc, id.vars=c('IT','logscale','trTYPE','LOGPO'))
				fit.mcmc[, REP:=infiles[i,REP]]
				fit.mcmc
			})
	dfl	<- do.call('rbind',tmp)	
	for(x in c('infTYPE','infETSI','infETSIC'))
	{
		tmp	<- subset(dfl, variable==x)	
		ggplot(tmp, aes(x=trTYPE, y=LOGPO, colour=value)) + 
				geom_point() + 
				labs(x='H coeff for I0', y='log posterior', colour='S coeff') +
				geom_vline(xintercept=log(5), colour='red') + 
				theme_bw() +
				facet_wrap(~REP, ncol=6, scales='free_y')
		ggsave(file=file.path(indir,paste0('logpoprofile_with_formula.inf_',x,'.pdf')),w=15,h=10)
	}
	for(x in c('infTYPE','infETSI','infETSIC'))
	{	
		tmp	<- unique(subset(dfl, variable==x, c(REP, trTYPE, LOGPO, variable, value)), by=c('REP','trTYPE','value'))	
		ggplot(tmp, aes(x=trTYPE, y=value, z=LOGPO)) + 
				stat_density2d(colour='black') +
				labs(x='H coeff for I0', y='S coeff', colour='log posterior') +
				geom_vline(xintercept=log(5), colour='red') + 
				theme_bw() +
				facet_wrap(~REP, ncol=6, scales='free_y')
		ggsave(file=file.path(indir,paste0('logposurface_with_formula.inf_',x,'.pdf')),w=15,h=10)
	}
	
	tmp		<- lapply(seq_len(nrow(infiles)), function(i)
			{
				infile	<- infiles[i,F]
				cat(basename(infile),'\n')
				#infile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_results/m3.RR5.n1250_seed123_rep1_aMCMC170710_170713_trTYPE_infTYPE_mndInf_mh10_hetinf0_mcs100_s0.5.rda'
				load(infile)	
				fit.bdt	<- data.table(SLICE_ID=seq_along(fit$bdt), SLICE_NTIP=sapply(fit$bdt, Ntip))
				fit.bdt[, REP:=infiles[i,REP]]
				fit.bdt[, FORMULA_INF:=infiles[i,FORMULA_INF]]
				fit.bdt				
			})
	fit.bdt	<- do.call('rbind',tmp)
	fit.bdt	<- fit.bdt[, list(SLICE_N=length(SLICE_ID), SLICE_ALLTIPS=sum(SLICE_NTIP), SLICE_MTIPS=mean(SLICE_NTIP)), by=c('FORMULA_INF','REP')]
	fit.bdt	<- melt(fit.bdt, id.vars=c('FORMULA_INF','REP'))
	ggplot(fit.bdt, aes(x=REP, y=value)) + geom_point() + facet_grid(variable~FORMULA_INF, scales='free')
	ggsave(file=file.path(indir,'slice_info.pdf'),w=10,h=10)
}

cr.master.ex3.compare<- function()
{
	require(coalreg)
	require(viridis)
	indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
	infiles	<- data.table(F=list.files(indir, pattern='rda$',full.names=TRUE))
	res		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep99_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'			
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
				#F		<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep1.nwk'
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
				#F	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep99_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'			
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
				#F		<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep1.nwk'
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
				#F	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep99_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
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
	load('~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-empirical_transmission_rates.rda')
	subset(dfr, GROUP=='TR_STAGE' & SC=='PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43')
	
	infiles	<- data.table(F=list.files('/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results', pattern='*rda', full.names=TRUE))
	
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
				#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43.newick'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
				#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43.newick'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
				#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43.newick'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
				#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43.newick'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
				#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43.newick'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
				#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed51.newick'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
				#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43.newick'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
				#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43.newick'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
				#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43.newick'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
				#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43.newick'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
				#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43.newick'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations'
		outdir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_results'
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
				#F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/png_simulations/PANGEA-AcuteHigh-InterventionNone-cov11.8-seed43.newick'
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
	require(ape)
	require(lhs)
	require(akima)
	require(mvtnorm)
	require(coalreg)
	require(viridis)
	require(data.table)	
	if(0)
	{
		indir				<- '~/Box Sync/OR_Work/2017/2017_coalregression/master_examples'			
		par.base.pattern	<- 'm3.RR5.n1250_seed123'	
		par.maxNodeDepth	<- 15
		par.maxHeight		<- 10
		par.s				<- 0.5
		par.lnHetInflation0	<- NA
		outfile				<- paste0('~/Box Sync/OR_Work/2017/2017_coalregression/master_results_5/m3.RR5.n1250_rep1_trflasso_d',par.maxNodeDepth,'_h',par.maxHeight,'_s',par.s,'_h',par.lnHetInflation0,'.rda')
	}
	#
	#	run coalreg	run using exact time to infection
	#	with maxNodeDepth=2, extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#	with 50% sampling	 
	set.seed(42)	
	infiles	<- data.table(F=list.files(indir, pattern=paste0(par.base.pattern,'_rep[0-9]+.nwk'),full.names=TRUE))
	infiles[, {
				#F		<- '/Users/Oliver/Box Sync/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep1.nwk'
				ph		<- read.tree( F )
				if(par.s<1)
					ph 	<- drop.tip(ph, sample(ph$tip.label, replace=FALSE, size=length(ph$tip.label)*par.s))
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
				fit 	<- trf.lasso(	dph, 	tmp, 
										trf_names = c( 'TYPE'), 
										aoi_names = c( 'ETSI' ), 
										maxNodeDepth=par.maxNodeDepth, 
										maxHeight=par.maxHeight, 
										lasso_threshold=5, 
										method = 'BFGS', 
										lnHetInflation0= par.lnHetInflation0,
										lnr0 = -2, 
										lnrLimits = c(-4, 2), 
										scale=FALSE)	
				fci 	<- fisher.ci(fit)	 
				pci 	<- prof.ci(fit, fci  )
				# save(fit, file=outfile)
				#print(fit$bestfit$par )
				#print( fci$ci )	
				tmp		<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs3_maxNodeDepth',par.maxNodeDepth,'_maxHeight',par.maxHeight,'_s',par.s*100,'.rda'),basename(F)))
				save( fit, fci, pci, file=tmp)
			}, by='F']				
}

cr.master.ex3.dev.mcmc.with.hetinflation<- function()
{
	#	TODO try not do zero mean
	require(coalreg)
	require(viridis)
	require(data.table)
	require(ape)
	require(coda)
	infile						<- '~/Box Sync/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep1.nwk'	
	par.maxNodeDepth			<- Inf
	par.maxHeight				<- 25
	par.hetInflation_logprior	<- NA #function(x) dnorm(x, mean=0, sd=10, log=TRUE)
	par.mincladesize			<- 100
	par.s						<- 1
	formula.tr					<- ~TYPE
	formula.inf					<- ~TYPE		
	
	phylo	<- read.tree( infile )	
	if(par.s<1)
		phylo 	<- multi2di(drop.tip(phylo, sample(phylo$tip.label, replace=FALSE, size=length(phylo$tip.label)*par.s)), random=FALSE)
	#	create data.table with infection type
	phi		<- data.table(	TAXA=phylo$tip.label, 
							TYPE=sapply(strsplit(phylo$tip.label,'_',),'[[',2),
							ETSI=as.numeric(sapply(strsplit(phylo$tip.label,'_',),'[[',3)))
	set(phi, NULL, 'TYPE', phi[, as.numeric(as.character(factor(TYPE,levels=c('I1','I0'), labels=c('0','1'))))])
	#	zero mean	
	set(phi, NULL, 'ETSI', phi[,ETSI-mean(ETSI)])	
	#	print
	print(phylo)
	print(phi)
	#	init coreg
	tmp				<- unique(c(unlist(strsplit(gsub(' ','',as.character(formula.tr))[2],'+',fixed=TRUE)),unlist(strsplit(gsub(' ','',as.character(formula.inf))[2],'+',fixed=TRUE))))	
	X				<- data.matrix(phi[,tmp,with=FALSE])
	rownames(X)		<- phi[, TAXA]
	colnames(X)		<- tmp	
	#	run coreg
	fit		<- coreg.adaptiveMHwithinGibbs.or170919(
			X, phylo,
			transmission=formula.tr,
			infection=formula.inf,
			hetInflation_logprior=par.hetInflation_logprior,
			adapt.batch=function(accn){  ifelse(accn<1e3, 20, 50) },
			adapt.schedule=function(b){ 1+1*(1/b)^(1/3) }, 
			mhsteps=5e3,													
			maxHeight=par.maxHeight,
			maxNodeDepth=par.maxNodeDepth,
			mincladesize=par.mincladesize,
			coef_logprior_sd=10, 			
			verbose=TRUE
			)	
	fit.mcmc	<- mcmc(cbind(fit$trace, fit$trace_tr, fit$trace_inf))					
	#pdf(file.path(dirname(outfile), gsub('\\.rda',paste0('_coeff_trace.pdf'),basename(outfile))), w=10, h=7)
	plot(fit.mcmc)
	#dev.off()	
	fit.mcmc.a	<- mcmc(fit$trace_qsd)
	#pdf(file.path(dirname(outfile), gsub('\\.rda',paste0('_sdadapt.pdf'),basename(outfile))), w=5, h=10)
	plot(fit.mcmc.a)
	#dev.off()
	
}

cr.master.ex3.adMCMC<- function(infile, formula.tr, formula.inf, par.s, par.maxNodeDepth, par.maxHeight, par.mincladesize, par.tsimb, par.tsimn, par.hetInflation_logprior, par.infprior.mean, par.infprior.sd, extra='')
{
	require(coalreg)
	require(data.table)
	require(ape)
	require(coda)	
	if(0)
	{
		infile						<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep1.nwk'
		infile						<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep2.nwk'
		par.maxNodeDepth			<- Inf
		par.maxHeight				<- 10
		par.hetInflation_logprior	<- function(x) dnorm(x, mean=0, sd=10, log=TRUE)
		par.mincladesize			<- 100
		par.s						<- 0.5
		par.tsimb					<- 2
		par.tsimn					<- 1
		formula.tr					<- ~TYPE
		formula.inf					<- ~ETSI	
		par.infprior.mean			<- log(5)
		par.infprior.sd				<- 1/5	
	}
	#
	#	run coalreg	run using exact time to infection
	#	with maxNodeDepth=2, extra args lnr0 = -2, lnrLimits = c(-4, 2), scale=F, lasso_threshold=5, method = 'BFGS'
	#	with 50% sampling	 
	set.seed(42)	
	cat(	'\ninfile=',infile,
			'\nformula.tr=',as.character(formula.tr), 
			'\nformula.inf=',as.character(formula.inf), 
			'\npar.s=',par.s, 
			'\npar.maxNodeDepth=',par.maxNodeDepth, 
			'\npar.maxHeight=',par.maxHeight, 
			'\npar.mincladesize=',par.mincladesize, 
			'\npar.tsimb=',par.tsimb, 
			'\npar.tsimn=',par.tsimn, 
			'\npar.hetInflation_logprior=',as.character(par.hetInflation_logprior),
			'\npar.infprior.mean=',par.infprior.mean,
			'\npar.infprior.sd=',par.infprior.sd,			
			'\nextra=',extra)
	if(par.hetInflation_logprior==0)
		par.hetInflation_logprior	<- NA
	if(!is.na(par.hetInflation_logprior))
		par.hetInflation_logprior	<- function(x) dnorm(x, mean=0, sd=10, log=TRUE)	
	infection_logpriors	<- list()
	if(!is.na(par.infprior.mean) & !is.na(par.infprior.sd))
		infection_logpriors<- list( ETSI=function(x) dnorm(x, par.infprior.mean, par.infprior.sd, log=TRUE) )
	phylo	<- read.tree( infile )
	if(par.s<1)
		phylo 	<- multi2di(drop.tip(phylo, sample(phylo$tip.label, replace=FALSE, size=length(phylo$tip.label)*par.s)), random=FALSE)
	#	create data.table with infection type
	phi		<- data.table(	TAXA=phylo$tip.label, 
							TYPE=sapply(strsplit(phylo$tip.label,'_',),'[[',2),
							ETSI=as.numeric(sapply(strsplit(phylo$tip.label,'_',),'[[',3)))
	set(phi, NULL, 'TYPE', phi[, as.numeric(as.character(factor(TYPE,levels=c('I1','I0'), labels=c('0','1'))))])
	phi[, ETSIC:= as.numeric(as.character(factor(ETSI<0.25,levels=c(TRUE,FALSE),labels=c('0','1'))))]
	#	add noise: 	lognormal model with mean= ETSI and std dev= ETSI*par.tsimn
	tmp		<- phi[, list(ETSI_NOISE=rlnorm(1, meanlog= log(ETSI)-0.5*log(par.tsimn*par.tsimn+1), sdlog= sqrt(log(par.tsimn*par.tsimn+1)))), by='TAXA']				
	phi		<- merge(phi, tmp, by='TAXA')			
	#	add multiplicative bias	
	set(phi, NULL, 'ETSI_NOISE', phi[, ETSI_NOISE*par.tsimb])
	#	zero mean	
	set(phi, NULL, 'ETSI', phi[,ETSI-mean(ETSI)])
	set(phi, NULL, 'ETSI_NOISE', phi[,ETSI_NOISE-mean(ETSI_NOISE)])
	#	print
	print(phylo)
	print(phi)
	#	init coreg
	tmp				<- unique(c(unlist(strsplit(gsub(' ','',as.character(formula.tr))[2],'+',fixed=TRUE)),unlist(strsplit(gsub(' ','',as.character(formula.inf))[2],'+',fixed=TRUE))))	
	X				<- data.matrix(phi[,tmp,with=FALSE])
	rownames(X)		<- phi[, TAXA]
	colnames(X)		<- tmp	
	#	run coreg
	#fit		<- coreg.adaptiveMHwithinGibbs.or170710(
	fit		<- coreg.adaptiveMHwithinGibbs.or170919(
								X, phylo,
								transmission=formula.tr,
								infection=formula.inf,
								adapt.batch=function(accn){  ifelse(accn<1e3, 20, 50) },
								adapt.schedule=function(b){ 1+1*(1/b)^(1/3) },
								infection_logpriors=infection_logpriors,
								hetInflation_logprior=par.hetInflation_logprior,
								mhsteps=10e3,													
								maxHeight=par.maxHeight,
								maxNodeDepth=par.maxNodeDepth,
								mincladesize=par.mincladesize,
								coef_logprior_sd=10, 
								verbose=FALSE
								)	
	outfile	<- file.path(dirname(infile), gsub( '\\.nwk', 
												paste0(	'_aMCMC',
														extra,
														'_tr',gsub(' ','',as.character(formula.tr)[2]),
														'_inf',gsub(' ','',as.character(formula.inf)[2]),
														'_mnd',par.maxNodeDepth,
														'_mh',par.maxHeight,
														'_hetinf',as.numeric(is.function(par.hetInflation_logprior)),
														'_mcs',par.mincladesize,
														ifelse(!is.na(par.infprior.mean), paste0('_infm',par.infprior.mean), ''),
														ifelse(!is.na(par.infprior.sd), paste0('_infs',par.infprior.sd), ''),
														'_s',par.s,														
														'.rda'),basename(infile)))
	cat('\noutfile=',outfile)
	save(fit, file=outfile)	
	cat('\nDONE\n')
}

cr.master.ex3.dev.MCMC1<- function()
{
	require(coalreg)
	require(viridis)
	require(data.table)
	require(ape)
	require(coda)
	
	indir						<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'			
	par.base.pattern			<- 'm3.RR5.n150_seed123'	
	par.maxNodeDepth			<- Inf
	par.maxHeight				<- 10
	par.hetInflation_logprior	<- NA
	par.mincladesize			<- 100
	par.s						<- 0.5
	
	F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep1.nwk'
	F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep2.nwk'
	
	set.seed(42)
	phylo	<- read.tree( F )					
	phylo 	<- multi2di(drop.tip(phylo, sample(phylo$tip.label, replace=FALSE, size=length(phylo$tip.label)*par.s)), random=FALSE)
	#	create data.table with infection type
	phi		<- data.table(	TAXA=phylo$tip.label, 
			TYPE=sapply(strsplit(phylo$tip.label,'_',),'[[',2),
			ETSI=as.numeric(sapply(strsplit(phylo$tip.label,'_',),'[[',3)))
	set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
	set(phi, NULL, 'ETSI', phi[,ETSI-mean(ETSI)])
	#	coreg
	X				<- data.matrix(phi[,2:3,with=FALSE])
	rownames(X)	<- phi[, TAXA]
	colnames(X)	<- colnames(phi)[2:3]				
	fit 			<- coreg.or170611(	X, phylo,
			transmission=~TYPE,
			infection=~ETSI,
			transmission_coef0= list(),
			infection_coef0= list(ETSI=log(2)),
			mhsteps=1e2,
			mhscale=setNames(c(1/4,1,1/4),c('logscale','TYPE','ETSI')),
			maxHeight=par.maxHeight,
			maxNodeDepth=par.maxNodeDepth,
			mincladesize=par.mincladesize,
			infection_logpriors=list(ETSI= function(x) dnorm(x, mean=0, sd=0.1, log=TRUE)),
			hetInflation_logprior=par.hetInflation_logprior
	)
	
	fit 			<- coreg(	X, phylo,
			transmission=~TYPE,
			infection=~ETSI,
			mhsteps=3e3,
			mhscale=1/4,
			maxHeight=10,
			maxNodeDepth=Inf,
			mincladesize=100,
			infection_logpriors=list(ETSI= function(x) dnorm(x, mean=0, sd=0.1, log=TRUE)),
			hetInflation_logprior=NA
	)
	
	fit				<- coreg.adaptiveMHwithinGibbs.or170710(
			X, phylo,
			transmission=~TYPE,
			infection=~ETSI,
			transmission_coef0= list(),
			infection_coef0= list(),
			adapt.batch=function(accn){  ifelse(accn<1e3, 20, 50) },
			adapt.schedule=function(b){ 1+1*(1/b)^(1/3) }, 
			mhsteps=3e3,													
			maxHeight=10,
			maxNodeDepth=Inf,
			mincladesize=100,
			coef_logprior_sd=10,
			infection_logpriors=list(ETSI= function(x) dnorm(x, mean=0, sd=0.1, log=TRUE)),
			hetInflation_logprior=NA,
			debug.nolkl=FALSE
	)
	
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_20100.rda'),basename(F)))
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_50100_v2.rda'),basename(F)))
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_50100_v3.rda'),basename(F)))
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_100100.rda'),basename(F)))
	save(fit, file=outfile)	
	fit				<- coreg.adaptiveMHwithinGibbs.or170710(
			X, phylo,
			transmission=~TYPE,
			infection=~ETSI,
			transmission_coef0= list(),
			infection_coef0= list(),
			adapt.batch=20,
			adapt.schedule=function(b){ 1+1*(1/b)^(1/3) }, 
			mhsteps=2e5,													
			maxHeight=10,
			maxNodeDepth=Inf,
			mincladesize=100,
			coef_logprior_sd=5,
			infection_logpriors=list(ETSI= function(x) dnorm(x, mean=0, sd=0.1, log=TRUE)),
			hetInflation_logprior=NA,
			debug.nolkl=TRUE
	)	
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_2025.rda'),basename(F)))
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_2050.rda'),basename(F)))
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_20100.rda'),basename(F)))
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_20100_debug.rda'),basename(F)))
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_20100_logscalednorm_debug.rda'),basename(F)))
	save(fit, file=outfile)									
	
	
	fit.mcmc	<- cbind(fit$trace, fit$trace_tr, fit$trace_inf)				
	if(is.na(par.hetInflation_logprior))
		fit.mcmc<- fit.mcmc[,-which(colnames(fit.mcmc)=='logHetInflation')]
	fit.mcmc	<- mcmc(fit.mcmc)
	pdf(file.path(dirname(outfile), gsub('\\.rda',paste0('_coeff_trace.pdf'),basename(outfile))), w=10, h=7)
	plot(fit.mcmc)
	dev.off()
	
	#
	fit.mcmc.a	<- mcmc(fit$trace_qsd)
	pdf(file.path(dirname(outfile), gsub('\\.rda',paste0('_sdadapt.pdf'),basename(outfile))), w=5, h=10)
	plot(fit.mcmc.a)
	dev.off()
	#actual prior for TYPE and ETSI
	prior.logTYPE.sd	<- apply(X,2,sd)[1]*5
	dfp	<- rbind(	data.table(variable='logscale', x=seq(-15,15,1e-2), y=dnorm(seq(-15,15,1e-2), 0, 5)),
			data.table(variable='TYPE', x=seq(-15,15,1e-2), y=dnorm(seq(-15,15,1e-2), 0, prior.logTYPE.sd)),
			data.table(variable='ETSI', x=seq(-0.4,0.4,1e-3), y=dnorm(seq(-0.4,0.4,1e-3), 0, 0.1))
	)				
	prior.logTYPE.d		<- function(x){ dnorm(x,0, prior.logTYPE.sd) }
	prior.logETSI.d		<- function(x){ dnorm(x, 0, 0.1) }				
	#sampled parameters
	tmp <- seq.int(1e3,5e4,4*5)
	dfm	<- melt(cbind(as.data.table(fit.mcmc[tmp,]), data.table(IT=tmp)), id.vars='IT')
	ggplot(dfm) + geom_histogram(aes(x=value, y=..density..), bins=50) + facet_wrap(~variable, scale='free', ncol=3) +
			geom_line(data=dfp, aes(x=x, y=y, colour='red'))
	ggsave(file=file.path(dirname(outfile), gsub('\\.rda',paste0('_coeffhistograms_from_prior.pdf'),basename(outfile))), w=10, h=7)
	
	#
	da	<- as.data.table(t(matrix(data=fit$accept, nrow=4, dimnames=list(c('logscale','logHetInflation','ETSI','TYPE'),NULL))))				
	da[, BATCH:= ceiling(seq_len(nrow(da))/20 )]
	da	<- melt(da, id.vars='BATCH')
	da	<- da[, list(ACC=mean(value)), by=c('BATCH','variable')]
	ggplot(da, aes(x=BATCH, y=ACC, colour=variable)) + geom_line() + geom_hline(yintercept=0.44) + facet_grid(variable~.)
	ggsave(file=file.path(dirname(outfile), gsub('\\.rda',paste0('_accept.pdf'),basename(outfile))), w=10, h=10)
	
	#
	tr	<- as.data.table(fit$trace_qsd)
	tr[, IT:=seq.int(1,nrow(tr))]
	tr	<- subset(tr, IT%%(20*4)==0)
	
	
	1-rejectionRate(window(fit.mcmc, start=1e3, end=5e3, thin=4))
	
}

cr.master.ex3.dev.MCMC2<- function()
{
	require(coalreg)	
	require(data.table)
	require(ape)
	require(coda)
		
	par.s	<- 0.5		
	F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep2.nwk'
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_50100_v2.rda'),basename(F)))
	
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
	fit				<- coreg.adaptiveMHwithinGibbs.or170710(
			X, phylo,
			transmission=~TYPE,
			infection=~ETSI,
			adapt.batch=function(accn){  ifelse(accn<1e3, 20, 50) },
			adapt.schedule=function(b){ 1+1*(1/b)^(1/3) }, 
			mhsteps=3e3,													
			maxHeight=10,
			maxNodeDepth=Inf,
			mincladesize=10,
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
	
	#actual prior for TYPE and ETSI
	prior.logTYPE.sd	<- apply(X,2,sd)[1]*5
	dfp	<- rbind(	data.table(variable='logscale', x=seq(-15,15,1e-2), y=dnorm(seq(-15,15,1e-2), 0, 5)),
			data.table(variable='TYPE', x=seq(-15,15,1e-2), y=dnorm(seq(-15,15,1e-2), 0, prior.logTYPE.sd)),
			data.table(variable='ETSI', x=seq(-0.4,0.4,1e-3), y=dnorm(seq(-0.4,0.4,1e-3), 0, 0.1))
	)				
	prior.logTYPE.d		<- function(x){ dnorm(x,0, prior.logTYPE.sd) }
	prior.logETSI.d		<- function(x){ dnorm(x, 0, 0.1) }				
	#sampled parameters
	tmp <- seq.int(1e3,5e4,4*5)
	dfm	<- melt(cbind(as.data.table(fit.mcmc[tmp,]), data.table(IT=tmp)), id.vars='IT')
	ggplot(dfm) + geom_histogram(aes(x=value, y=..density..), bins=50) + facet_wrap(~variable, scale='free', ncol=3) +
			geom_line(data=dfp, aes(x=x, y=y, colour='red'))
	ggsave(file=file.path(dirname(outfile), gsub('\\.rda',paste0('_coeffhistograms_from_prior.pdf'),basename(outfile))), w=10, h=7)
	
	#
	da	<- as.data.table(t(matrix(data=fit$accept, nrow=4, dimnames=list(c('logscale','logHetInflation','ETSI','TYPE'),NULL))))				
	da[, BATCH:= ceiling(seq_len(nrow(da))/20 )]
	da	<- melt(da, id.vars='BATCH')
	da	<- da[, list(ACC=mean(value)), by=c('BATCH','variable')]
	ggplot(da, aes(x=BATCH, y=ACC, colour=variable)) + geom_line() + geom_hline(yintercept=0.44) + facet_grid(variable~.)
	ggsave(file=file.path(dirname(outfile), gsub('\\.rda',paste0('_accept.pdf'),basename(outfile))), w=10, h=10)
	
	#
	tr	<- as.data.table(fit$trace_qsd)
	tr[, IT:=seq.int(1,nrow(tr))]
	tr	<- subset(tr, IT%%(20*4)==0)
	
	
	1-rejectionRate(window(fit.mcmc, start=1e3, end=5e3, thin=4))
	
}

cr.master.ex3.dev.MCMC3<- function()
{
	require(coalreg)	
	require(data.table)
	require(ape)
	require(coda)
	
	par.s	<- 0.5		
	F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep3.nwk'
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_50100.rda'),basename(F)))
	
	set.seed(42)
	phylo	<- read.tree( F )					
	phylo 	<- multi2di(drop.tip(phylo, sample(phylo$tip.label, replace=FALSE, size=length(phylo$tip.label)*par.s)), random=FALSE)
	#	create data.table with infection type
	phi		<- data.table(	TAXA=phylo$tip.label, 
			TYPE=sapply(strsplit(phylo$tip.label,'_',),'[[',2),
			ETSI=as.numeric(sapply(strsplit(phylo$tip.label,'_',),'[[',3)))
	set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
	phi[, ETSI_C:= as.numeric(as.character(factor(ETSI<0.25,levels=c(TRUE,FALSE),labels=c('0','1'))))]
	
	#phi[, table(TYPE,ETSI_C)]	
	#     0   1
	#I1  53 232
	#I0 267  73
	
	#set(phi, NULL, 'ETSI', phi[,ETSI-mean(ETSI)])
	#	coreg
	tmp		<- c('TYPE','ETSI_C')
	X				<- data.matrix(phi[,tmp,with=FALSE])
	rownames(X)	<- phi[, TAXA]
	colnames(X)	<- tmp
	# 	make sure baseline is coded as x=0, and that acute stage is coded x=1 
	#	hence exp(beta*x)=exp(beta) for individuals with x=1 and 1 for individuals with x=0
	X[, 'TYPE']	<- X[, 'TYPE']-1 	
	#
	#	make 
	
	fit				<- coreg.adaptiveMHwithinGibbs.or170710(
							X, phylo,
							transmission=~TYPE,
							infection=~ETSI_C,
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
	
	#actual prior for TYPE and ETSI
	prior.logTYPE.sd	<- apply(X,2,sd)[1]*5
	dfp	<- rbind(	data.table(variable='logscale', x=seq(-15,15,1e-2), y=dnorm(seq(-15,15,1e-2), 0, 5)),
			data.table(variable='TYPE', x=seq(-15,15,1e-2), y=dnorm(seq(-15,15,1e-2), 0, prior.logTYPE.sd)),
			data.table(variable='ETSI', x=seq(-0.4,0.4,1e-3), y=dnorm(seq(-0.4,0.4,1e-3), 0, 0.1))
	)				
	prior.logTYPE.d		<- function(x){ dnorm(x,0, prior.logTYPE.sd) }
	prior.logETSI.d		<- function(x){ dnorm(x, 0, 0.1) }				
	#sampled parameters
	tmp <- seq.int(1e3,5e4,4*5)
	dfm	<- melt(cbind(as.data.table(fit.mcmc[tmp,]), data.table(IT=tmp)), id.vars='IT')
	ggplot(dfm) + geom_histogram(aes(x=value, y=..density..), bins=50) + facet_wrap(~variable, scale='free', ncol=3) +
			geom_line(data=dfp, aes(x=x, y=y, colour='red'))
	ggsave(file=file.path(dirname(outfile), gsub('\\.rda',paste0('_coeffhistograms_from_prior.pdf'),basename(outfile))), w=10, h=7)
	
	#
	da	<- as.data.table(t(matrix(data=fit$accept, nrow=4, dimnames=list(c('logscale','logHetInflation','ETSI','TYPE'),NULL))))				
	da[, BATCH:= ceiling(seq_len(nrow(da))/20 )]
	da	<- melt(da, id.vars='BATCH')
	da	<- da[, list(ACC=mean(value)), by=c('BATCH','variable')]
	ggplot(da, aes(x=BATCH, y=ACC, colour=variable)) + geom_line() + geom_hline(yintercept=0.44) + facet_grid(variable~.)
	ggsave(file=file.path(dirname(outfile), gsub('\\.rda',paste0('_accept.pdf'),basename(outfile))), w=10, h=10)
	
	#
	tr	<- as.data.table(fit$trace_qsd)
	tr[, IT:=seq.int(1,nrow(tr))]
	tr	<- subset(tr, IT%%(20*4)==0)
	
	
	1-rejectionRate(window(fit.mcmc, start=1e3, end=5e3, thin=4))
	
}

cr.master.ex3.dev.MCMC5<- function()
{
	require(coalreg)	
	require(data.table)
	require(ape)
	require(coda)
	
	par.s	<- 0.5	
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
	infile	<- file.path(indir, 'm3.RR5.n150_seed123_rep3.nwk')
	infile	<- file.path(indir, 'm3.RR5.n1250_seed123_rep3.nwk')
	
	set.seed(42)
	tr 			<- read.tree( infile )	
	tr 			<- multi2di(drop.tip(tr, sample(tr$tip.label, replace=FALSE, size=length(tr$tip.label)*par.s)), random=FALSE)
	X 			<- data.frame( TYPE=grepl('I0', tr$tip.label) )
	rownames(X) <- tr$tip.label	
	f0 			<- coreg ( X, tr
							, ~TYPE
							, ~TYPE
							, mhsteps = 5e3, ncpu=1)
	set.seed(42)
	tr 			<- read.tree( infile )	
	tr 			<- multi2di(drop.tip(tr, sample(tr$tip.label, replace=FALSE, size=length(tr$tip.label)*par.s)), random=FALSE)
	X 			<- data.frame( TYPE=grepl('I0', tr$tip.label) )
	rownames(X) <- tr$tip.label		
	f2 			<-  coreg.adaptiveMHwithinGibbs.or170710(X, tr
							, ~TYPE
							, ~TYPE
							, adapt.batch=function(accn){  ifelse(accn<1e3, 20, 50) },
							adapt.schedule=function(b){ 1+1*(1/b)^(1/3) }, 
							mhsteps=5e3,													
							maxHeight=10, #***
							maxNodeDepth=Inf,
							mincladesize=100,
							coef_logprior_sd=10,	
							verbose=TRUE)	
	outfile	<- file.path(dirname(infile), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_TYPEaoi_simpleversion.rda'),basename(infile)))
	save(f0, f2, file=outfile)
	
					
	
	fit	<- f0
	fit.mcmc	<- cbind(fit$trace, fit$trace_tr, fit$trace_inf)				
	fit.mcmc<- fit.mcmc[,-which(colnames(fit.mcmc)=='logHetInflation')]
	fit.mcmc	<- mcmc(fit.mcmc)
	pdf(file.path(dirname(outfile), gsub('\\.rda',paste0('_coeff_trace_f0.pdf'),basename(outfile))), w=10, h=7)
	plot(fit.mcmc)
	dev.off()	
	
	fit	<- f2
	fit.mcmc	<- cbind(fit$trace, fit$trace_tr, fit$trace_inf)					
	fit.mcmc<- fit.mcmc[,-which(colnames(fit.mcmc)=='logHetInflation')]
	colnames(fit.mcmc)	<- c('logscale','TR_I0','INF_I0')
	fit.mcmc	<- mcmc(fit.mcmc)
	pdf(file.path(dirname(outfile), gsub('\\.rda',paste0('_coeff_trace_f2.pdf'),basename(outfile))), w=10, h=7)
	plot(fit.mcmc)
	dev.off()
	fit.mcmc.a	<- mcmc(fit$trace_qsd)
	pdf(file.path(dirname(outfile), gsub('\\.rda',paste0('_sdadapt_f2.pdf'),basename(outfile))), w=5, h=10)
	plot(fit.mcmc.a)
	dev.off()
	tmp <- seq.int(1e3,5e3,4*5)
	dfm	<- cbind(as.data.table(fit.mcmc[tmp,]), data.table(IT=tmp))
	
	
	
	fs				<- data.table( IT=seq_along(f2$logprior), PR=f2$logprior, LKL=f2$loglik, PO=f2$logposterior, ACC=f2$accept )
	tmp				<- cbind(fit$trace, fit$trace_tr, fit$trace_inf)
	colnames(tmp)	<- c('logscale','logHetInflation','TR_I0','INF_I0')
	fs				<- cbind(fs, as.data.table(tmp))
	
	subset(fs, TR_I0>(log(5)-.1) & TR_I0<(log(5)+.1))[, mean(LKL)]
	#	-2943.644
	subset(fs, TR_I0>(5-.1) & TR_I0<(5+.1))[, mean(LKL)]
	#	-2943.338
	#	--> lkl is indeed larger for theta=5 than theta=log 5

	#	is this the same on the full tree?	
	coefs0		<- c(logscale=0.11, logHetInflation=0); coef_inf<- 1.033; 
	coef_trs	<- as.data.table(expand.grid(COEF_TR=c(1,log(5),2,2.5,3,3.5,4,5,6,7,8,10), SAMPLE=c(1,0.8,0.5)))
	infile		<- file.path(indir, 'm3.RR5.n1250_seed123_rep5.nwk')
	coef_trs	<- coef_trs[, {
				coef_tr		<- COEF_TR
				par.s		<- SAMPLE
				set.seed(42)
				cat('\nprocess',coef_tr, ', ',par.s)				
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
											exp(coefs0['logscale']), 
											maxHeight = maxHeight, 
											maxNodeDepth = maxNodeDepth, 
											hetInflation = exp(coefs0['logHetInflation']))
								}))
				list(LKL=loglkl)
			}, by=c('COEF_TR','SAMPLE')]
	outfile	<- file.path(dirname(infile), gsub('\\.nwk',paste0('_coalreg_likelihoods_for_mean_other_params.rda'),basename(infile)))
	save(infile, coefs0, coef_inf, coef_trs, file=outfile)			
}

cr.master.ex3.dev.MCMC4<- function()
{
	require(coalreg)	
	require(data.table)
	require(ape)
	require(coda)
	
	par.s	<- 0.5		
	F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep4.nwk'
	F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep2.nwk'
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_50100_v4.rda'),basename(F)))
	
	set.seed(42)
	phylo	<- read.tree( F )					
	phylo 	<- multi2di(drop.tip(phylo, sample(phylo$tip.label, replace=FALSE, size=length(phylo$tip.label)*par.s)), random=FALSE)
	#	create data.table with infection type
	phi		<- data.table(	TAXA=phylo$tip.label, 
			TYPE=sapply(strsplit(phylo$tip.label,'_',),'[[',2),
			ETSI=as.numeric(sapply(strsplit(phylo$tip.label,'_',),'[[',3)))
	set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
	phi[, ETSI_C:= as.numeric(as.character(factor(ETSI<0.25,levels=c(TRUE,FALSE),labels=c('0','1'))))]
	
	#set(phi, NULL, 'ETSI', phi[,ETSI-mean(ETSI)])
	#	coreg
	tmp		<- unique(c('TYPE','TYPE'))
	X				<- data.matrix(phi[,tmp,with=FALSE])
	rownames(X)	<- phi[, TAXA]
	colnames(X)	<- tmp
	# 	make sure baseline is coded as x=0, and that acute stage is coded x=1 
	#	hence exp(beta*x)=exp(beta) for individuals with x=1 and 1 for individuals with x=0
	X[, 'TYPE']	<- X[, 'TYPE']-1 	
	#
	#	make 
	
	fit				<- coreg.adaptiveMHwithinGibbs.or170710(
							X, phylo,
							transmission=~TYPE,
							infection=~TYPE,
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
	
	#actual prior for TYPE and ETSI
	prior.logTYPE.sd	<- apply(X,2,sd)[1]*5
	dfp	<- rbind(	data.table(variable='logscale', x=seq(-15,15,1e-2), y=dnorm(seq(-15,15,1e-2), 0, 5)),
			data.table(variable='TYPE', x=seq(-15,15,1e-2), y=dnorm(seq(-15,15,1e-2), 0, prior.logTYPE.sd)),
			data.table(variable='ETSI', x=seq(-0.4,0.4,1e-3), y=dnorm(seq(-0.4,0.4,1e-3), 0, 0.1))
	)				
	prior.logTYPE.d		<- function(x){ dnorm(x,0, prior.logTYPE.sd) }
	prior.logETSI.d		<- function(x){ dnorm(x, 0, 0.1) }				
	#sampled parameters
	tmp <- seq.int(1e3,5e4,4*5)
	dfm	<- melt(cbind(as.data.table(fit.mcmc[tmp,]), data.table(IT=tmp)), id.vars='IT')
	ggplot(dfm) + geom_histogram(aes(x=value, y=..density..), bins=50) + facet_wrap(~variable, scale='free', ncol=3) +
			geom_line(data=dfp, aes(x=x, y=y, colour='red'))
	ggsave(file=file.path(dirname(outfile), gsub('\\.rda',paste0('_coeffhistograms_from_prior.pdf'),basename(outfile))), w=10, h=7)
	
	#
	da	<- as.data.table(t(matrix(data=fit$accept, nrow=4, dimnames=list(c('logscale','logHetInflation','ETSI','TYPE'),NULL))))				
	da[, BATCH:= ceiling(seq_len(nrow(da))/20 )]
	da	<- melt(da, id.vars='BATCH')
	da	<- da[, list(ACC=mean(value)), by=c('BATCH','variable')]
	ggplot(da, aes(x=BATCH, y=ACC, colour=variable)) + geom_line() + geom_hline(yintercept=0.44) + facet_grid(variable~.)
	ggsave(file=file.path(dirname(outfile), gsub('\\.rda',paste0('_accept.pdf'),basename(outfile))), w=10, h=10)
	
	#
	tr	<- as.data.table(fit$trace_qsd)
	tr[, IT:=seq.int(1,nrow(tr))]
	tr	<- subset(tr, IT%%(20*4)==0)
	
	
	1-rejectionRate(window(fit.mcmc, start=1e3, end=5e3, thin=4))
	
}

cr.master.ex3.dev.MCMC6<- function()
{
	require(coalreg)	
	require(data.table)
	require(ape)
	require(coda)
	
	par.s	<- 1			
	F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep5.nwk'
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_ETFIaoi_adaptiveMCMC_50100_v6.rda'),basename(F)))
	
	set.seed(42)
	phylo	<- read.tree( F )
	if(par.s<1)
		phylo <- multi2di(drop.tip(phylo, sample(phylo$tip.label, replace=FALSE, size=length(phylo$tip.label)*par.s)), random=FALSE)
	#	create data.table with infection type
	phi		<- data.table(	TAXA=phylo$tip.label, 
			TYPE=sapply(strsplit(phylo$tip.label,'_',),'[[',2),
			ETSI=as.numeric(sapply(strsplit(phylo$tip.label,'_',),'[[',3)))
	set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
	phi[, ETSI_C:= as.numeric(as.character(factor(ETSI<0.25,levels=c(TRUE,FALSE),labels=c('0','1'))))]
	
	#set(phi, NULL, 'ETSI', phi[,ETSI-mean(ETSI)])
	#	coreg
	tmp		<- unique(c('TYPE','TYPE'))
	X				<- data.matrix(phi[,tmp,with=FALSE])
	rownames(X)	<- phi[, TAXA]
	colnames(X)	<- tmp
	# 	make sure baseline is coded as x=0, and that acute stage is coded x=1 
	#	hence exp(beta*x)=exp(beta) for individuals with x=1 and 1 for individuals with x=0
	X[, 'TYPE']	<- X[, 'TYPE']-1 	
	#
	#	make 
	
	fit				<- coreg.adaptiveMHwithinGibbs.or170710(
			X, phylo,
			transmission=~TYPE,
			infection=~TYPE,
			adapt.batch=function(accn){  ifelse(accn<1e3, 20, 50) },
			adapt.schedule=function(b){ 1+1*(1/b)^(1/3) }, 
			mhsteps=5e3,													
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

cr.master.ex3.dev.MCMC7<- function()
{
	require(coalreg)	
	require(data.table)
	require(ape)
	require(coda)
	
	par.s	<- 0.5			
	F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep3.nwk'
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_TYPEaoi_mincladesize200_v7.rda'),basename(F)))
	
	set.seed(42)
	phylo	<- read.tree( F )
	if(par.s<1)
		phylo <- drop.tip(phylo, sample(phylo$tip.label, replace=FALSE, size=length(phylo$tip.label)*par.s))
		#phylo <- multi2di(drop.tip(phylo, sample(phylo$tip.label, replace=FALSE, size=length(phylo$tip.label)*par.s)), random=FALSE)
	#	create data.table with infection type
	phi		<- data.table(	TAXA=phylo$tip.label, 
							TYPE=sapply(strsplit(phylo$tip.label,'_',),'[[',2),
							ETSI=as.numeric(sapply(strsplit(phylo$tip.label,'_',),'[[',3)))
	set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
	phi[, ETSI_C:= as.numeric(as.character(factor(ETSI<0.25,levels=c(TRUE,FALSE),labels=c('0','1'))))]
	
	#set(phi, NULL, 'ETSI', phi[,ETSI-mean(ETSI)])
	#	coreg
	tmp			<- unique(c('TYPE','TYPE'))
	X			<- data.matrix(phi[,tmp,with=FALSE])
	rownames(X)	<- phi[, TAXA]
	colnames(X)	<- tmp
	# 	make sure baseline is coded as x=0, and that acute stage is coded x=1 
	#	hence exp(beta*x)=exp(beta) for individuals with x=1 and 1 for individuals with x=0
	X[, 'TYPE']	<- X[, 'TYPE']-1 	
	#
	#	make 
	
	fit				<- coreg.adaptiveMHwithinGibbs.or170710(
			X, phylo,
			transmission=~TYPE,
			infection=~TYPE,
			adapt.batch=function(accn){  ifelse(accn<1e3, 20, 50) },
			adapt.schedule=function(b){ 1+1*(1/b)^(1/3) }, 
			mhsteps=5e3,													
			maxHeight=10,
			maxNodeDepth=Inf,
			mincladesize=200,
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

cr.master.ex3.dev.MCMC8<- function()
{
	require(coalreg)	
	require(data.table)
	require(ape)
	require(coda)
	
	par.s	<- 0.5			
	F		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n1250_seed123_rep3.nwk'
	outfile	<- file.path(dirname(F), gsub('\\.nwk',paste0('_coalreg_using_TYPEtrf_TYPEaoi_mincladesize200_maxHeightInf_v8.rda'),basename(F)))
	
	set.seed(42)
	phylo	<- read.tree( F )
	if(par.s<1)
		phylo <- drop.tip(phylo, sample(phylo$tip.label, replace=FALSE, size=length(phylo$tip.label)*par.s))
	#phylo <- multi2di(drop.tip(phylo, sample(phylo$tip.label, replace=FALSE, size=length(phylo$tip.label)*par.s)), random=FALSE)
	#	create data.table with infection type
	phi		<- data.table(	TAXA=phylo$tip.label, 
			TYPE=sapply(strsplit(phylo$tip.label,'_',),'[[',2),
			ETSI=as.numeric(sapply(strsplit(phylo$tip.label,'_',),'[[',3)))
	set(phi, NULL, 'TYPE', phi[,factor(TYPE,levels=c('I1','I0'))])
	phi[, ETSI_C:= as.numeric(as.character(factor(ETSI<0.25,levels=c(TRUE,FALSE),labels=c('0','1'))))]
	
	#set(phi, NULL, 'ETSI', phi[,ETSI-mean(ETSI)])
	#	coreg
	tmp			<- unique(c('TYPE','TYPE'))
	X			<- data.matrix(phi[,tmp,with=FALSE])
	rownames(X)	<- phi[, TAXA]
	colnames(X)	<- tmp
	# 	make sure baseline is coded as x=0, and that acute stage is coded x=1 
	#	hence exp(beta*x)=exp(beta) for individuals with x=1 and 1 for individuals with x=0
	X[, 'TYPE']	<- X[, 'TYPE']-1 	
	#
	#	make 
	
	fit				<- coreg.adaptiveMHwithinGibbs.or170710(
			X, phylo,
			transmission=~TYPE,
			infection=~TYPE,
			adapt.batch=function(accn){  ifelse(accn<1e3, 20, 50) },
			adapt.schedule=function(b){ 1+1*(1/b)^(1/3) }, 
			mhsteps=5e3,													
			maxHeight=Inf,
			maxNodeDepth=Inf,
			mincladesize=200,
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

cr.master.ex3.runcoalreg.using.TYPE.ETFI.mbias.BFGS1<- function(indir, par.base.pattern, par.s, par.tsimb)
{
	require(coalreg)
	require(viridis)
	if(0)
	{
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
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
				#F		<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11.nwk'
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
				#F	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
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
				#F		<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11.nwk'
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
				#F	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
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
				#F		<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11.nwk'
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
				#F	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
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
		indir				<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples'
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
				#F		<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11.nwk'
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
				#F	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_coalregression/master_examples/m3.RR5.n150_seed123_rep11_coalreg_using_TYPEtrf_ETFIaoi_BFGSargs_s50.rda'
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
