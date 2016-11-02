######################################################################################
project.examl.Bezemer.141104<- function()
{
	require(ape)
	
	dir.name	<- DATA  	
	indir		<- paste(dir.name,'bezemer',sep='/')
	infile		<- "NLB10BLAST_cutRmu"
	insignat	<- "Tue_Nov_4_2014"
	
	if(0)
	{
		seq.PROT.RT	<- read.dna(paste(indir, '/', infile, '.fasta', sep=''), format='fasta')
		file		<- paste(indir, '/', infile, '_', insignat, '.R', sep='')
		save(seq.PROT.RT, file=file)
	}	
	#	ExaML bootstrap args
	bs.from		<- 0
	bs.to		<- 2
	bs.n		<- 200
	outdir		<- indir
	cmd			<- hivc.cmd.examl.bootstrap(indir, infile, insignat, insignat, bs.from=bs.from, bs.to=bs.to,bs.n=bs.n,outdir=outdir, resume=1, verbose=1)
		
	dummy		<- lapply(cmd, function(x)
			{				
				#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
				x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=34, hpc.q="pqeph", hpc.mem="450mb", hpc.nproc=1)
				signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
				outfile	<- paste("exa",signat,sep='.')
				#cat(x)
				hivc.cmd.hpccaller(outdir, outfile, x)
				Sys.sleep(1)
			})	
}
######################################################################################
project.examl.ATHENA1610.161102<- function()
{
	require(big.phylo)
	indir		<- '/work/or105/ATHENA_2016/data/examl'
	infile		<- "ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B.fasta"	
	#	ExaML bootstrap args
	bs.from		<- 0
	bs.to		<- 0
	bs.n		<- 500
	outdir		<- indir
	
	#args.parser	<- paste("-m DNA -q",PARTITION)
	args.parser	<- paste("-m DNA")
	args.examl	<- "-f d -D -m GAMMA"	#	 -- this is the default that worked in 24 hours	
	cmd			<- cmd.examl.bootstrap(indir, infile, bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, opt.bootstrap.by="codon", args.parser=args.parser, args.examl=args.examl, tmpdir.prefix="examl")					
	dummy		<- lapply(cmd, function(x)
			{				
				#x		<- cmd.hpcwrapper(x, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
				#x		<- cmd.hpcwrapper(x, hpc.walltime=79, hpc.q="pqeph", hpc.mem="1800mb", hpc.nproc=1)
				x		<- cmd.hpcwrapper(x, hpc.walltime=79, hpc.q="pqeelab", hpc.mem="5800mb", hpc.nproc=1)
				signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
				outfile	<- paste("exa",signat,sep='.')
				cat(x)
				cmd.hpccaller(outdir, outfile, x)
				Sys.sleep(1)
			})		
}