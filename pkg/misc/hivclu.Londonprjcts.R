######################################################################################
project.London.rmDRM.151130<- function()
{
	infile	<- '~/Dropbox (Infectious Disease)/2016_LondonMSM/subUKogC.fasta'
	seql	<- read.dna(infile, format='fasta')	
	stopifnot( any(grepl('HXB2',rownames(seql))) )
	
	require(big.phylo)
	drmout	<- seq.rm.drugresistance(seql)
	outfile	<- '~/Dropbox (Infectious Disease)/2016_LondonMSM/subUKogC_noDRM.fasta'	
	save(drmout, file=gsub('\\.fasta','_output.R',outfile))
	seq		<- drmout$nodr.seq
	save(seq, file=gsub('\\.fasta','.R',outfile))
	write.dna(seq, file=outfile, format='fasta', colsep='', nbcol=-1)
}
######################################################################################
project.London.FirstExaml.151202<- function()
{
	require(big.phylo)	
	indir		<- '~/Dropbox (Infectious Disease)/2016_LondonMSM'	#TODO: currently no whitespace or brackets in file name: escape all files with ""
	indir		<- '/Users/Oliver/duke/2015_various'
	infile		<- "subUKogC_noDRM"				#TODO: use fasta instead of R	
	#	ExaML bootstrap args
	bs.from		<- 0		# 0 is the actual data alignment
	bs.to		<- 50
	bs.n		<- 50		# optional
	outdir		<- indir
	
	# arguments for ExaML parser to create binaries
	args.parser	<- paste("-m DNA")
	#args.parser	<- paste("-m DNA -q",PARTITION)	# PARTITION is a text file that specifies gene partitions
	
	# arguments for ExaML
	#args.examl	<- "-f d -m GAMMA"		#	 -- ran for weeks
	#args.examl	<- "-f o -D -m GAMMA"	#	 -- not followed up until 'default' worked
	args.examl	<- "-f d -D -m GAMMA"	#	 -- this is the default that worked in 24 hours	
	cmd			<- cmd.examl.bootstrap(indir, infile, bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, opt.bootstrap.by="codon", args.parser=args.parser, args.examl=args.examl, tmpdir.prefix="examl")
	
	# create UNIX commands: list, each element for one bs tree	
	dummy		<- lapply(cmd, function(x)
			{				
				#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
				#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=79, hpc.q="pqeph", hpc.mem="1800mb", hpc.nproc=1)
				x		<- cmd.hpcwrapper(x, hpc.walltime=129, hpc.q="pqeelab", hpc.mem="5800mb", hpc.nproc=1)
				signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
				outfile	<- paste("exa",signat,sep='.')
				cat(x)
				cmd.hpccaller(outdir, outfile, x)
				Sys.sleep(1)
			})	
}
