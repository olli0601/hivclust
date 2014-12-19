
######################################################################################
hivc.pipeline.recombination<- function()
{
	if(0)	#generate candidate recombinants	- 	this is computationally expensive
	{
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"
		#infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences100"
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		
		batch.n			<- 100		#for 1e4 sequences, does about 100 in 25hrs, so request 100 batches for walltime 25 expected + 10hrs grace
		file			<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
		load(file)
		batch.seq		<- round(seq.int(0,nrow(seq.PROT.RT),len=batch.n),d=0)
		batch.seq		<- rbind(batch.seq[-length(batch.seq)], batch.seq[-1]-1)
		#batch.seq		<- batch.seq[,1:10]	#test run
		file			<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".phylip",sep='')
		lapply(seq_len(ncol(batch.seq)),function(j)
				{					
					cmd			<- hivc.cmd.recombination.run.3seq(infile=file, outfile=paste(indir,'/',infile,'_',batch.seq[1,j],'-',batch.seq[2,j],'_',gsub('/',':',insignat),".3seq",sep=''), recomb.3seq.siglevel=0.1, nproc=1, recomb.3seq.testvsall.beginatseq=batch.seq[1,j], recomb.3seq.testvsall.endatseq=batch.seq[2,j], verbose=1)
					cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=35, hpc.q="pqeph", hpc.mem="3850mb",  hpc.nproc=1)
					cat(cmd)
					outdir		<- indir
					outfile		<- paste("r3seq",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')									
					hivc.cmd.hpccaller(outdir, outfile, cmd)			
				})		
		stop()
	}
	if(1)	#validate candidate recombinants
	{
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		resume		<- 0
		verbose		<- 1
		
		argv				<<-	hivc.cmd.recombination.process.3SEQ.output(indir, infile, insignat, resume=1, verbose=1) 
		argv				<<- unlist(strsplit(argv,' '))
		df.recomb			<- hivc.prog.recombination.process.3SEQ.output()	
		
		triplets			<- seq_len(nrow(df.recomb))
		triplets			<- 147:nrow(df.recomb)
		dummy	<- lapply(triplets, function(i)
				{					
					if(verbose)	cat(paste("\nprocess triplet number",i,"\n"))
					argv				<<- hivc.cmd.recombination.check.candidates(indir, infile, insignat, i, resume=resume, verbose=1)
					argv				<<- unlist(strsplit(argv,' '))
					hivc.prog.recombination.check.candidates()		#this starts ExaML for the ith triplet			
				})
		stop()
	}
	if(1)	#process validation step to produce plots
	{
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		resume		<- 0
		verbose		<- 1
		
		argv				<<- hivc.cmd.recombination.plot.incongruence(indir, infile, insignat, triplet.id=NA, opt.select="ng2", verbose=1)
		argv				<<- unlist(strsplit(argv,' '))
		hivc.prog.recombination.plot.incongruence()		
		
		argv				<<- hivc.cmd.recombination.plot.incongruence(indir, infile, insignat, triplet.id=NA, opt.select="g2", verbose=1)
		argv				<<- unlist(strsplit(argv,' '))
		hivc.prog.recombination.plot.incongruence()
		
	}
	if(0)	#	collect likely recombinants or those likely confounding the phylogeny		-- 	identified by eye
	{
		verbose				<- 1
		
		argv				<<-	hivc.cmd.recombination.process.3SEQ.output(indir, infile, insignat, resume=1, verbose=1) 
		argv				<<- unlist(strsplit(argv,' '))
		df.recomb			<- hivc.prog.recombination.process.3SEQ.output()
		setnames(df.recomb, "dummy", "triplet.id")
		setkey(df.recomb, triplet.id)
		#	collect likely recombinants or those likely confounding the phylogeny
		recombinants.ng2	<-	c(	subset( df.recomb, triplet.id==60 )[,child],		#likely confounding
									subset( df.recomb, triplet.id==52 )[,child],		
									subset( df.recomb, triplet.id==55 )[,parent2],"R09-26706","R11-12152","R12-15108",		#others cluster in addition
									subset( df.recomb, triplet.id==48 )[,parent2],"2007G319",								#others cluster in addition
									subset( df.recomb, triplet.id==112 )[,child],
									subset( df.recomb, triplet.id==129 )[,child],
									subset( df.recomb, triplet.id==148 )[,child],		#length only 50
									subset( df.recomb, triplet.id==135 )[,child],
									subset( df.recomb, triplet.id==102 )[,child],		#likely confounding
									subset( df.recomb, triplet.id==85 )[,child],
									subset( df.recomb, triplet.id==81 )[,child],		#likely confounding
									subset( df.recomb, triplet.id==125 )[,child],		#likely confounding but might be recombinant
									subset( df.recomb, triplet.id==120 )[,child]	)
		recombinants.g2		<-	c(	"2006G052", "2007G263", "M3621708072010", "M4048713072011", 
									"M4203226082011", "R03-14636", "TN_B.HT.2005.05HT_129336.EU439719", "TN_B.HT.2004.04HT_129732.EU439728", "TN_B.PH.2008.08R_01_361.AB587101", "TN_B.ZA.2011.patient_1720_seq_1746.KC423374", "TN_B.ZA.2012.patient_1720_seq_2734.KC423805", "TN_B.DO.2008.HIV_PRRT_PJ01967_48.JN713614",					#these cluster with  M4048713072011									
									"R11-15440", "R12-00343", "R10-09812",				#these are children of R08-20970
									"R12-07939" )
		recombinants		<- c( recombinants.ng2, recombinants.g2 )
		if(verbose)	cat(paste("\ncollected recombinants or likely confounding sequences, n=",length(recombinants)))
		#
		#	save non-recombinant sequence dataset
		#
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"	
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		outfile		<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences"	
		outsignat	<- "Fri_Nov_01_16/07/23_2013"
		
		file		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
		load(file)
		tmp			<- setdiff(rownames(seq.PROT.RT), recombinants)
		seq.PROT.RT	<- seq.PROT.RT[tmp,]
		if(verbose)	cat(paste("\nnumber of sequences without (likely) recombinants, n=",nrow(seq.PROT.RT)))
		file		<- paste(indir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
		if(verbose)	cat(paste("\nsave new file to",file))
		save(seq.PROT.RT, file=file)
	}
}
######################################################################################
hivc.pipeline.ExaML<- function()
{
	dir.name<- DATA
	if(0)	#compute one ExaML tree, no bootstrapping
	{		
		indir	<- paste(dir.name,"tmp",sep='/')
		infile	<- "ATHENA_2013_03_FirstAliSequences_PROTRT"
		outdir	<- paste(dir.name,"tmp",sep='/')
		cmd		<- paste(cmd,hivc.cmd.examl(indir,infile,gsub('/',':',signat.out),gsub('/',':',signat.out),outdir=outdir,resume=1,verbose=1),sep='')
		cmd		<- paste(cmd,hivc.cmd.examl.cleanup(outdir),sep='')
	}
	if(0)	#compute ExaML trees with bootstrap values. Bootstrap is over initial starting trees to start ML search.
	{
		bs.from	<- 0
		bs.to	<- 0
		bs.n	<- 100
		signat.in	<- "Sat_May_11_14/23/46_2013"
		signat.out	<- "Sat_May_11_14/23/46_2013"				
		indir	<- paste(dir.name,"tmp",sep='/')
		infile	<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
		infile	<- "ATHENA_2013_03_FirstCurSequences_PROTRTCD3"
		outdir	<- paste(dir.name,"tmp",sep='/')
		cmd		<- hivc.cmd.examl.bsstarttree(indir,infile,gsub('/',':',signat.out),gsub('/',':',signat.out),bs.from=bs.from,bs.to=bs.to,bs.n=bs.n,outdir=outdir, resume=1, verbose=1)
		#check if we have all bs.n files. if yes, combine and cleanup
		
		outdir	<- paste(dir.name,"tmp",sep='/')							
		lapply(cmd, function(x)
				{				
					x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=36, hpc.q="pqeph")
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("pipeline",signat,sep='.')
					cat(x)
					#hivc.cmd.hpccaller(outdir, outfile, x)
					#Sys.sleep(1)
				})
		stop()
	}
	if(1)	#compute ExaML trees with bootstrap values. Bootstrap is over codon in alignment and over initial starting trees to start ML search.
	{
		bs.from		<- 0
		bs.to		<- 100
		bs.n		<- 500
		
		indir		<- paste(dir.name,"tmp",sep='/')
		#infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences"
		#signat.in	<- "Sat_Jun_16_17:23:46_2013"								
		#infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"
		#signat.in	<- "Thu_Aug_01_17/05/23_2013"	
		#infile		<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences"	
		#signat.in	<- "Fri_Nov_01_16/07/23_2013"
		infile		<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		signat.in	<- "Wed_Dec_18_11/37/00_2013"
		
		#infile		<- "UKCA_2013_07_TNTPHIVnGTR"
		#signat.in	<- "Mon_Sep_22_17/23/46_2013"
		
		outdir		<- paste(dir.name,"tmp",sep='/')
		cmd			<- hivc.cmd.examl.bootstrap(indir,infile,gsub('/',':',signat.in),gsub('/',':',signat.in),bs.from=bs.from,bs.to=bs.to,bs.n=bs.n,outdir=outdir, resume=1, verbose=1)				
		outdir		<- paste(dir.name,"tmp",sep='/')							
		lapply(cmd, function(x)
				{				
					x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=24, hpc.q=NA, hpc.mem="3850mb", hpc.nproc=8)
					#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=24, hpc.q="pqeph", hpc.mem="3850mb", hpc.nproc=8)
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("exa",signat,sep='.')
					#cat(x)
					hivc.cmd.hpccaller(outdir, outfile, x)
					Sys.sleep(1)
				})
		stop()
	}	
}
######################################################################################
hivc.pipeline.clustering<- function()
{	
	if(1)	#clustering: precompute clustering objects, evaluate TPTN, get default clustering, refine to capture MSM transmission
	{	
		resume		<- 1
		verbose		<- 1
		#seq project
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
		insignat	<- "Sat_Jun_16_17/23/46_2013"
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"
		#infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs500"
		insignat	<- "Thu_Aug_01_17/05/23_2013"		
		infile		<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences_examlbs500"	
		insignat	<- "Fri_Nov_01_16/07/23_2013"
		infile		<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500"	
		insignat	<- "Wed_Dec_18_11/37/00_2013"
		
		
		#seq covariates
		indircov	<- paste(DATA,"derived",sep='/')
		infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"
		#default clustering tptn		
		patient.n	<- 15700
		#default clustering	
		opt.brl		<- "dist.brl.casc" 
		thresh.brl	<- 0.096
		thresh.bs	<- 0.8		
		
		cmd			<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov)		
		cmd			<- paste(cmd, hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume),sep='')
		cmd			<- paste(cmd, hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.max", patient.n=patient.n, resume=resume),sep='')
		cmd			<- paste(cmd, hivc.cmd.clustering(indir, infile, insignat, opt.brl, thresh.brl, thresh.bs, resume=resume),sep='')		
		cmd			<- paste(cmd, hivc.cmd.clustering.msm(indir, infile, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=resume),sep='')			
		cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=71, hpc.q="pqeph", hpc.mem="15850mb",  hpc.nproc=1)
		
		cat(cmd)
		stop()
		outdir		<- paste(DATA,"tmp",sep='/')
		outfile		<- paste("clust",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		hivc.cmd.hpccaller(outdir, outfile, cmd)
		quit("no")
	}
}
######################################################################################
hivc.pipeline.BEAST<- function()
{
	if(0)	#run BEAST 1.7.5 GMRF skyline
	{
		indir				<- paste(DATA,"tmp",sep='/')		
		indircov			<- paste(DATA,"derived",sep='/')
		infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
		#infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		#insignat			<- "Thu_Aug_01_17/05/23_2013"		
		#infiletree			<- paste(infile,"examlbs100",sep="_")
		#infilexml			<- paste(infile,'_',"beast",'_',"seroneg",sep='')		
		infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infilexml			<- paste(infile,'_',"beast",'_',"all",sep='')
		infiletree			<- paste(infile,"examlbs500",sep="_")
		insignat			<- "Wed_Dec_18_11:37:00_2013"
		
		#infilexml.template	<- "um22rhU2050"
		#infilexml.template	<- "um22rhG202018"
		#infilexml.template	<- "rhU65rho753"
		#infilexml.template	<- "rhU65rho903"
		#infilexml.template	<- "rhU65rho906"
		#infilexml.template	<- "rhU65rho909"	
		#infilexml.template	<- "um181rhU2045"
		#infilexml.template	<- "um182rhU2045"
		infilexml.template	<- "um192rhU2080"
		#infilexml.template	<- "unhum192rhU2080"
		#infilexml.template	<- "um183rhU2045"
		#infilexml.template	<- "um182us45"
		#infilexml.template	<- "um182us60"
		#infilexml.template	<- "um182rhU2045ay"
		#infilexml.template	<- "um232rhU2045"
		#infilexml.template	<- "um232rhU2045ay"
		#infilexml.opt		<- "txs4clu"
		#infilexml.opt		<- "txs4clufx03"
		#infilexml.opt		<- "mph4clu"
		#infilexml.opt		<- "mph4clutx4tip"
		#infilexml.opt		<- "mph4clufx03"
		#infilexml.opt		<- "mph4cluLdTd"
		#infilexml.opt		<- "mph4cluUmTd"
		#infilexml.opt		<- "mph4cluLsTd"
		#infilexml.opt		<- "mph4cluNoTd"
		#infilexml.opt		<- "mph4cluu4tipLdTd"
		#infilexml.opt		<- "mph4clutx4tipLdTd"
		#infilexml.opt		<- "mph4clutx4tipLsTd"
		infilexml.opt		<- "mph4clutx4tip"
		
		outdir				<- indir
		outsignat			<- insignat
		
		opt.brl				<- "dist.brl.casc" 
		thresh.brl			<- 0.096
		thresh.bs			<- 0.8
		#pool.ntip			<- 130		
		#pool.ntip			<- 150
		pool.ntip			<- 190
		#pool.ntip			<- 400
		resume				<- 1
		verbose				<- 1
		
		argv				<<- hivc.cmd.beast.poolrunxml(indir, infile, insignat, indircov, infilecov, infiletree, infilexml, outsignat, pool.ntip, infilexml.opt=infilexml.opt, infilexml.template=infilexml.template, opt.brl=opt.brl, thresh.brl=thresh.brl, thresh.bs=thresh.bs, resume=resume, verbose=1)
		argv				<<- unlist(strsplit(argv,' '))
		hivc.prog.BEAST.generate.xml()		
		quit("no")
	}
	if(1)		#generate BEAST2 BDSKYline xml file
	{
		indir				<- paste(DATA,"tmp",sep='/')
		indircov			<- paste(DATA,"derived",sep='/')
		outdir				<- indir
		infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"

		#infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		#insignat			<- "Thu_Aug_01_17/05/23_2013"
		#infiletree			<- paste(infile,"examlbs100",sep="_")
		#infilexml			<- paste(infile,'_',"seroneg",sep='')		
		#outsignat			<- "Tue_Aug_26_09/13/47_2013"		
		infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree			<- paste(infile,"examlbs500",sep="_")
		insignat			<- "Wed_Dec_18_11:37:00_2013"
		infilexml			<- paste(infile,'_',"all",sep='')
		outsignat			<- insignat
		
		opt.brl				<- "dist.brl.casc" 
		thresh.brl			<- 0.096
		thresh.bs			<- 0.8
		pool.ntip			<- NA		
		resume				<- 1
		verbose				<- 1
		
		infilexml.template	<- "bdsky_hky" 
		infilexml.template	<- "sasky_hky"
		infilexml.template	<- "sasky_sdr06"
		infilexml.template	<- "sasky_sdr06fr"
		infilexml.opt		<- "S4p"
		#infilexml.opt		<- "S5p"
		#infilexml.opt		<- "S8p"
		infilexml.opt		<- "s0106"
		#infilexml.opt		<- "s0108"
		#infilexml.opt		<- "s00106"
		infilexml.opt		<- "s124"
		infilexml.opt		<- "s024"
		infilexml.opt		<- "s424"
		infilexml.opt		<- "s184"
		infilexml.opt		<- "sartest"
		infilexml.opt		<- "d999"
		infilexml.opt		<- "d774"
		infilexml.opt		<- "d543"
		infilexml.opt		<- "dg543"
		infilexml.opt		<- "r54350"
		infilexml.opt		<- "r54399"
		infilexml.opt		<- "run210"
		infilexml.opt		<- "rbe910"
		infilexml.opt		<- "rbe415"
		infilexml.opt		<- "rbe420"
		infilexml.opt		<- "rbe425"
		infilexml.opt		<- "pmct15"
		infilexml.opt		<- "pmct25"
		infilexml.opt		<- "pmct35"
		infilexml.opt		<- "pfn635"
		#infilexml.opt		<- "pfn670"
		#infilexml.opt		<- "pfn6100"
		infilexml.opt		<- "piv470"
		infilexml.opt		<- "piv490"
		infilexml.opt		<- "rfn815"
		infilexml.opt		<- "rfn835"
		infilexml.opt		<- "rse815"
		infilexml.opt		<- "rse835"
		infilexml.opt		<- "rsu815"
		infilexml.opt		<- "alsu50"
		infilexml.opt		<- "alrh40"
		infilexml.opt		<- "alrh80"
		#infilexml.opt		<- "alrh160"
		infilexml.opt		<- "clrh80"
		#infilexml.opt		<- "rsu835"		
		#infilexml.opt		<- "ori40"
		#infilexml.opt		<- "ori50"
		#infilexml.opt		<- "ori60"
		#infilexml.opt		<- "ori70"
		argv				<<- hivc.cmd.beast.poolrunxml(indir, infile, insignat, indircov, infilecov, infiletree, infilexml, outsignat, pool.ntip, infilexml.opt=infilexml.opt, infilexml.template=infilexml.template, opt.brl=opt.brl, thresh.brl=thresh.brl, thresh.bs=thresh.bs, resume=resume, verbose=1)
		argv				<<- unlist(strsplit(argv,' '))		
		hivc.prog.BEAST2.generate.xml()
	}
	if(0)		#run BEAST2 BDSKYline
	{
		indir				<- paste(DATA,"tmp",sep='/')
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_standard"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxGpInvS4"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxGpInvS6"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxGpInvS8"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxGpInvS10"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxS4"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxS5"
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_bdsky_seroneg-130-1_standard_fxS8"
		insignat			<- "Tue_Aug_26_09/13/47_2013"
		hpc.ncpu			<- 1

		cmd					<- hivc.cmd.beast2.runxml(indir, infile, insignat, hpc.ncpu=hpc.ncpu, prog.opt.Xmx="1200m")
		cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=hpc.ncpu, hpc.walltime=72, hpc.mem="1200mb")
		
		cat(cmd)
		stop()
		signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
		outdir		<- paste(DATA,"tmp",sep='/')
		outfile		<- paste("b2",signat,sep='.')					
		hivc.cmd.hpccaller(outdir, outfile, cmd)
	}
	
}
######################################################################################
hivc.pipeline.BEASTout.get.cluster.trees<- function()
{
	require(ape)
	require(data.table)
	require(hivclust)
	#	program options
	opt.pool					<- NA
	opt.burnin					<- 5e6			
	resume						<- TRUE
	verbose						<- TRUE
	#	input files			
	indir					<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/beast/beast2_140201'
	infile					<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_seroneg-130"
	insignat				<- "Tue_Aug_26_09:13:47_2013"
	infilexml.opt			<- "rsu815"
	infilexml.template		<- "sasky_sdr06"
	#
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,14),
									infilexml.opt= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.opt<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									infilexml.template= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.template<- tmp[1]
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									burnin= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) opt.burnin<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,5),
									pool= return(as.numeric(substr(arg,7,nchar(arg)))),NA)	}))
		if(length(tmp)>0) opt.pool<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}
	#
	outdir					<- indir	
	outsignat				<- insignat
	#
	if(verbose)
	{
		print(opt.pool)
		print(opt.burnin)		
		print(resume)
		print(indir)		
		print(infile)
		print(insignat)
		print(infilexml.opt)
		print(infilexml.template)		
	}
	
	cmd			<- hivc.cmd.beast2.getclustertrees(indir, infile, insignat, infilexml.opt, infilexml.template, opt.burnin, outdir=indir, outsignat=insignat, prog= PR.BEAST2CLUTREES, opt.pool=opt.pool, verbose=verbose, resume=resume)
	#cat(cmd)
	cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=1, hpc.walltime=21, hpc.mem="3800mb")
	outdir		<- paste(DATA,"tmp",sep='/')
	outfile		<- paste("b2m.",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='')					
	hivc.cmd.hpccaller(outdir, outfile, cmd)
}
######################################################################################
hivc.pipeline.BEASTout<- function()
{
	if(1)
	{
		indircov			<- paste(DATA,"derived",sep='/')
		infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
		#
		indir				<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/beast/beast_131011"		
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		insignat			<- "Tue_Aug_26_09:13:47_2013"				
		infilexml.template	<- "um232rhU2045"
		infilexml.opt		<- "mph4clutx4tip"
		burnin				<- 2e7
		#
		indir				<- paste(DATA,"tmp",sep='/')
		indir				<- paste(DATA,"zip",sep='/')
		
		infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"		
		insignat			<- "Wed_Dec_18_11:37:00_2013"		
		infilexml.template	<- "sasky_sdr06fr"
		infilexml.opt		<- "alrh160"
		infilexml.opt		<- "clrh80"
		burnin				<- 5e6
		#
		#infilexml.template	<- "um192rhU2080"
		#infilexml.opt		<- "mph4clutx4tip"
		#burnin				<- 2e7
		#
		outdir				<- indir
		outsignat			<- insignat
		opt.pool			<- NA
		#
		files		<- list.files(indir)
		files		<- files[ sapply(files, function(x) grepl(infile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(paste('_',infilexml.opt,'_',sep=''), x, fixed=1) & grepl(paste('_',infilexml.template,'_',sep=''), x, fixed=1) & grepl('_pool_[0-9]+',x) & grepl('trees$',x) ) ]				
		if(!length(files))	stop('no input files matching criteria')
		tmp			<- regmatches( files, regexpr('_pool_[0-9]+',files)) 
		pool		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
		file.info	<- data.table(file=files, pool=pool)
		setkey(file.info, pool)
		#
		dummy		<- sapply( file.info[,unique(pool)], function(pool)
				{
					cmd				<- hivc.cmd.beast2.getclustertrees(indir, infile, insignat, infilexml.opt, infilexml.template, burnin=burnin, opt.pool=pool, verbose=1, resume=1)
					#cmd			<- paste(cmd, hivc.cmd.beast2.processclustertrees(indir, infile, insignat, infilexml.opt, infilexml.template, verbose=1, resume=1), sep='')
					#cmd			<- paste(cmd, hivc.cmd.beast2.plotclustertrees(indir, infile, insignat, indircov, infilecov, infilexml.opt, infilexml.template, resume=1, verbose=1), sep='')
					cat(cmd)
					stop()
					cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=1, hpc.walltime=21, hpc.mem="3800mb")
					outdir		<- paste(DATA,"tmp",sep='/')
					outfile		<- paste("b2p.",strsplit(date(),split=' ')[[1]],collapse='_',sep='')					
					hivc.cmd.hpccaller(outdir, outfile, cmd)			
				})
	}
	if(0)
	{
		indir				<- paste(DATA,"tmp",sep='/')
		indircov			<- paste(DATA,"derived",sep='/')
		outdir				<- indir
		infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"		
		infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"		
		insignat			<- "Wed_Dec_18_11:37:00_2013"		
		outsignat			<- insignat
		#infilexml.template	<- "sasky_sdr06"		
		#infilexml.opt		<- "alsu50"
		infilexml.template	<- "sasky_sdr06fr"
		infilexml.opt		<- "alrh160"
		infilexml.opt		<- "clrh80"
		#infilexml.template	<- "um192rhU2080"
		#infilexml.opt		<- "mph4clutx4tip"	
		
		#indir					<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/beast/beast2_140201'
		#infile					<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_seroneg-130"
		#insignat				<- "Tue_Aug_26_09:13:47_2013"
		#infilexml.opt			<- "rsu815"
		#infilexml.template		<- "sasky_sdr06"
		
		files		<- list.files(indir)
		files		<- files[ sapply(files, function(x) grepl(infile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(paste('_',infilexml.opt,'_',sep=''), x, fixed=1) & grepl(paste('_',infilexml.template,'_',sep=''), x, fixed=1) & grepl('_pool_[0-9]+',x) & grepl('_clu_[0-9]+',x) & grepl('R$',x) ) ]				
		if(!length(files))	stop('no input files matching criteria')
		tmp			<- regmatches( files, regexpr('_clu_[0-9]+',files)) 
		cluster		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
		file.info	<- data.table(file=files, cluster=cluster)
		setkey(file.info, cluster)
		
		#file.info	<- subset(file.info, cluster%in%c(23, 77, 126, 152, 315))
		#print(file.info)
		
		dummy		<- sapply( file.info[,unique(cluster)], function(clu)
				{
					cmd			<- hivc.cmd.beast2.processclustertrees(indir, infile, insignat, infilexml.opt, infilexml.template, cluster=clu, verbose=1, resume=1)					
					cat(cmd)
					cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=1, hpc.walltime=80, hpc.mem="1600mb")
					outdir		<- paste(DATA,"tmp",sep='/')
					outfile		<- paste("b2m.",strsplit(date(),split=' ')[[1]],collapse='_',sep='')					
					hivc.cmd.hpccaller(outdir, outfile, cmd)			
				})
		
	}	
}
######################################################################################
hivc.pipeline.props_univariate<- function()
{
	#stop()
	indir					<- paste(DATA,"fisheretal_data",sep='/')		
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	outdir					<- paste(DATA,"fisheretal",sep='/')
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"
	
	if(0)
	{
		method					<- '3d'
		method.recentctime		<- '2011-01-01'
		#method.recentctime		<- '2013-03-01'
		method.nodectime		<- 'any'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"					
		infilexml.opt			<- "mph4clutx4tip"
		infilexml.template		<- "um192rhU2080"	
		outfile					<- paste(infile,'Ac=MY_D=35_gmrf',sep='_')
	}
	if(0)
	{
		method					<- '3d'
		method.recentctime		<- '2013-03-01'
		method.nodectime		<- 'any'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		infilexml.opt			<- "clrh80"
		infilexml.template		<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
	}
	if(0)
	{
		method					<- '3d'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		infilexml.opt			<- "clrh80"
		infilexml.template		<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
	}
	if(0)
	{
		method					<- '3i'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		infilexml.opt			<- "clrh80"
		infilexml.template		<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
	}
	if(0)
	{
		method					<- '3j'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		infilexml.opt			<- "clrh80"
		infilexml.template		<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
	}
	if(0)
	{
		method					<- '3k'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		infilexml.opt			<- "clrh80"
		infilexml.template		<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
	}
	if(0)
	{
		method					<- '3l'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		infilexml.opt			<- "clrh80"
		infilexml.template		<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
	}
	if(0)
	{
		method					<- '3m'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		infilexml.opt			<- "clrh80"
		infilexml.template		<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
	}
	if(0)
	{
		method					<- '3n'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		infilexml.opt			<- "clrh80"
		infilexml.template		<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
	}
	if(0)
	{
		method					<- '3o'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		infilexml.opt			<- "clrh80"
		infilexml.template		<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
	}
	if(1)
	{
		method					<- '3p'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		infilexml.opt			<- "clrh80"
		infilexml.template		<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
	}
	method.lRNA.supp			<- 100
	method.use.AcuteSpec		<- 1
	method.minQLowerU			<- 0.135	#0.135, 0.18, 0.23
	method.thresh.pcoal			<- 0.3
	method.PDT					<- 'SEQ'
	#method.Acute				<- 'empirical'
	method.Acute				<- 'higher'
	#method.Acute				<- 'lower'
	#method.Acute				<- 'central'
	method.minLowerUWithNegT	<- 1
	#method.risk	<- c('m3.nicv','m3.nicv.clu','m3.tnicv','m3.tnicv.clu','m21st.cas.clu','m2wmx.cas.clu','m2t.cas.clu','m2wmx.tp.clu','m3.i.clu','m3.ni.clu','m3.nic.clu','m3.tni.clu','m3.tnic.clu','m3.tniv.clu','m3.tnicvNo.clu','m21st.cas','m2wmx.cas','m2t.cas','m2wmx.tp','m3.i','m3.ni','m3.nic','m3.tni','m3.tnic','m3.tniv','m3.tnicvNo')
	#method.risk	<- c('m3.nic.clu.adj','m3.tnic.clu.adj','m3.tnicvNo.clu.adj','m21st.cas.clu.adj','m2t.cas.clu.adj','m2wmx.cas.clu.adj','m2wmx.tp.clu.adj','m3.nicv.clu.adj','m3.tnicv.clu.adj')
	#method.risk	<- c('m21st.cas.adj','m2t.cas.adj','m2wmx.cas.adj','m2wmx.tp1.adj', 'm2wmx.tp2.adj', 'm2wmx.tp3.adj', 'm2wmx.tp4.adj','m3.nic.adj','m3.nicv.adj','m3.tnic.adj','m3.tnicv.adj','m3.tnicvNo.adj')
	method.risk	<- c(	'm2B1st.cas','m2Bwmx.cas','m2Bt.cas','m2Bwmx.tp1', 'm2Bwmx.tp2', 'm2Bwmx.tp3', 'm2Bwmx.tp4',
						'm2B1st.cas.adj','m2Bwmx.cas.adj','m2Bt.cas.adj','m2Bwmx.tp1.adj', 'm2Bwmx.tp2.adj', 'm2Bwmx.tp3.adj', 'm2Bwmx.tp4.adj',
						'm2B1st.cas.clu.adj','m2Bwmx.cas.clu.adj','m2Bt.cas.clu.adj','m2Bwmx.tp1.clu.adj', 'm2Bwmx.tp2.clu.adj', 'm2Bwmx.tp3.clu.adj', 'm2Bwmx.tp4.clu.adj')
	method.risk	<- c(	'm2Bwmx.tp1', 'm2Bwmx.tp2', 'm2Bwmx.tp3', 'm2Bwmx.tp4',
						'm2Bwmx.tp1.adj', 'm2Bwmx.tp2.adj', 'm2Bwmx.tp3.adj', 'm2Bwmx.tp4.adj',
						'm2Bwmx.tp1.clu.adj', 'm2Bwmx.tp2.clu.adj', 'm2Bwmx.tp3.clu.adj', 'm2Bwmx.tp4.clu.adj')
	method.risk	<- c(	'm2Bwmx.tp1.cens','m2Bwmx.tp2.cens','m2Bwmx.tp3.cens','m2Bwmx.tp4.cens','m2Bwmx.tp1.clu.cens','m2Bwmx.tp2.clu.cens','m2Bwmx.tp3.clu.cens','m2Bwmx.tp4.clu.cens')	
	method.risk	<- c(	'm3.nicv','m3.tnicv','m3.tnicvNo','m3.nicv.adj','m3.tnicv.adj','m3.tnicvNo.adj','m3.nicv.clu','m3.tnicv.clu','m3.tnicvNo.clu','m2Bwmx.tp1.cens')	
	method.risk	<- c(	'm2B1st.cas.cens','m2B1st.cas.clu.cens','m2B1st.cas.censp','m2B1st.cas.clu.censp','m2Bt.cas.cens','m2Bt.cas.clu.cens','m2Bt.cas.censp','m2Bt.cas.clu.censp'		)
	#	basic model 3 runs 	mem=1800
	#method.risk	<- c(	'm3.tnic','m3.tnic.clu','m3.tnic.censp','m3.tnic.clu.censp','m3.tnicv.censp','m3.tnicv.clu.censp' )
	#	basic model 3 runs, mem=3800 				
	#method.risk	<- c(	'm3.tnicMV.censp','m3.tnicMV.clu.censp','m3.tnicNo.censp','m3.tnicNo.clu.censp' )
	#method.risk	<- c(	'm3.tnicNo.censp','m3.tnicNo.clu.censp' )
	#method.risk		<- c(	'm3.tnicMV','m3.tnicMV.adj','m3.tnicMV.clu.adj'	)#,'m3.tnicMV.clu.censp','m3.tnicMv','m3.tnicMv.adj','m3.tnicMv.censp','m3.tnicMv.clu.censp' )
	#method.risk		<- c(	'm3.tnicNoMV','m3.tnicNoMV.adj','m3.tnicNoMV.clu.adj' )
	##method.risk		<- c(	'm3.tnicMV.adj','m3.tnicMV.clu.adj','m3.tnicNoMV.adj','m3.tnicNoMV.clu.adj' )
	##method.risk		<- c(	'm3.tnicMV','m3.tnicMV.clu','m3.tnicNoMV','m3.tnicNoMV.clu' )
	#method.risk		<- c(	'm3.atnicMV','m3.atnicMV.clu','m3.atnicNoMV','m3.atnicNoMV.clu' )
	#method.risk			<- c('m3.btnicMV.clu','m3.btnicNoMV.clu','m3.btnicMV.clu.wstar','m3.btnicNoMV.clu.wstar')
	#!##method.risk				<- c('m3.ind','m3.indNo','m3.indmx','m3.indmxNo','m3.n3mx','m3.indMV','m3.indNoMV','m3.indmxMV','m3.indmxNoMV','m3.n3mxMV',
	#!##							 'm3.ind.wstar','m3.indNo.wstar','m3.indmx.wstar','m3.indmxNo.wstar','m3.n3mx.wstar','m3.indMV.wstar','m3.indNoMV.wstar','m3.indmxMV.wstar','m3.indmxNoMV.wstar','m3.n3mxMV.wstar')	
	#	basic m2Bwmx runs	mem=1800
	method.risk		<- c(	'm2BwmxMv.tp1','m2BwmxMv.tp2','m2BwmxMv.tp3','m2BwmxMv.tp4')
	method.risk		<- c(	'm2CwmxMv.wtn.tp1','m2CwmxMv.wtn.tp2','m2CwmxMv.wtn.tp3','m2CwmxMv.wtn.tp4')
	#method.risk		<- c(	'm2CwmxMv.wtn.tp1','m2CwmxMv.wtn.tp2','m2CwmxMv.wtn.tp3','m2CwmxMv.wtn.tp4')
	#method.risk		<- c(	'm2BwmxMv.tp1' )
	#!##method.risk		<- c(	'm2BtMv.tp1','m2BtMv.tp2','m2BtMv.tp3','m2BtMv.tp4','m2BtMv.tp1.wstar','m2BtMv.tp2.wstar','m2BtMv.tp3.wstar','m2BtMv.tp4.wstar')
	#	basic m2Bwmx runs	mem=3800
	#method.risk	<- c(	'm2Bwmx.cas','m2Bwmx.cas.clu','m2Bwmx.cas.censp')
	#	basic m2BXXXMv runs	mem=7800
	#method.risk	<- c(	'm2BwmxMv.cas','m2BwmxMv.cas.adj','m2BwmxMv.cas.censp','m2BwmxMv.cas.clu.censp', 
	#		 			'm2B1stMv.cas', 'm2B1stMv.cas.adj', 'm2B1stMv.cas.censp', 'm2B1stMv.cas.clu.censp', 
	#					'm2BtMv.cas', 'm2BtMv.cas.adj', 'm2BtMv.cas.censp', 'm2BtMv.cas.clu.censp' )
	##method.risk	<- c(	'm2BwmxMv.cas.clu.censp', 'm2B1stMv.cas.clu.censp', 'm2BtMv.cas.clu.censp' )
	##method.risk	<- c(	'm2BwmxMv.cas.clu', 'm2B1stMv.cas.clu', 'm2BtMv.cas.clu','m2BwmxMv.cas', 'm2B1stMv.cas', 'm2BtMv.cas' )
	#	basic m2B1st runs	mem=3800
	#method.risk	<- c(	'm2B1st.cas','m2B1st.cas.clu','m2B1st.cas.censp')
	#	basic m2Bt runs	mem=3800
	#method.risk	<- c(	'm2Bt.cas','m2Bt.cas.clu','m2Bt.cas.censp')
	#	m5 runs
	#method.risk		<- c(	'm5.tA.tp1.clu','m5.tA.tp2.clu','m5.tA.tp3.clu','m5.tA.tp4.clu','m5.tA.tp1.clu.wstar','m5.tA.tp2.clu.wstar','m5.tA.tp3.clu.wstar','m5.tA.tp4.clu.wstar')
	#!##method.risk		<- c(	'm5.tAc.tp1','m5.tAc.tp2','m5.tAc.tp3','m5.tAc.tp4','m5.tAc.tp1.wstar','m5.tAc.tp2.wstar','m5.tAc.tp3.wstar','m5.tAc.tp4.wstar')
	#	Acute higher than VL, which we can check after diagnosis
	#method.risk	<- c( 'm4.Bwmxv','m4.Bwmxv.adj','m4.Bwmxv.censp','m4.Bwmxv.clu.censp','m4.BwmxvNo','m4.BwmxvNo.adj','m4.BwmxvNo.censp','m4.BwmxvNo.clu.censp','m4.BwmxvMv','m4.BwmxvMv.adj','m4.BwmxvMv.censp','m4.BwmxvMv.clu.censp' )
	#	NRTI+NNRTI puzzle
	#method.risk	<- c( 	'm3.tnicMv', 'm3.tnicMv.adj','m3.tnicMv.clu.adj', 'm3.tnicMv.censp','m3.tnicMv.clu.censp'	)
	#	all runs combined
	#method.risk		<- c(	'm2BwmxMv.tp1','m2BwmxMv.tp2','m2BwmxMv.tp3','m2BwmxMv.tp4','m2BtMv.tp1','m2BtMv.tp2','m2BtMv.tp3','m2BtMv.tp4','m3.n3mx','m3.n3mxMV','m5.tAc.tp1','m5.tAc.tp2','m5.tAc.tp3','m5.tAc.tp4')
	#method.risk			<- c('m5.tAc.tp1','m5.tAc.tp2','m5.tAc.tp3','m5.tAc.tp4')
	# use to pre-compute tables
	#method.risk		<- c( 	'm2B1st.cas.clu.adj','m2Bt.cas.clu.adj','m2Bwmx.cas.clu.adj','m2Bwmx.tp1.clu.adj', 'm2Bwmx.tp2.clu.adj', 'm2Bwmx.tp3.clu.adj', 'm2Bwmx.tp4.clu.adj', 'm4.Bwmxv.clu.adj')	
	#method.risk				<- c('m3.ind.clu.adj','m3.indNo.clu.adj','m3.indmx.clu.adj','m3.indmxNo.clu.adj','m3.n3mx.clu.adj')
	#method.risk			<- c( 	'm2Bwmx.tp1.clu.adj', 'm2Bwmx.tp2.clu.adj', 'm2Bwmx.tp3.clu.adj', 'm2Bwmx.tp4.clu.adj','m2Cwmx.tp1.clu.adj', 'm2Cwmx.tp2.clu.adj', 'm2Cwmx.tp3.clu.adj', 'm2Cwmx.tp4.clu.adj')
	#method.risk			<- c( 	'm2Bwmx.tp1.clu.adj', 'm2Bwmx.tp2.clu.adj', 'm2Bwmx.tp3.clu.adj', 'm2Bwmx.tp4.clu.adj','m2Bt.tp1.clu.adj', 'm2Bt.tp2.clu.adj', 'm2Bt.tp3.clu.adj', 'm2Bt.tp4.clu.adj','m3.n3mx.clu.adj','m3.indmx.clu.adj','m3.indmxNo.clu.adj','m5.tAc.tp1.clu.adj','m5.tAc.tp2.clu.adj','m5.tAc.tp3.clu.adj','m5.tAc.tp4.clu.adj'	)
	#method.risk			<- c( 	'm5.tAc.tp1.clu.adj','m5.tAc.tp2.clu.adj','m5.tAc.tp3.clu.adj','m5.tAc.tp4.clu.adj'	)
	#method.risk			<- c( 	'm5.tA.clu.adj','m5.tAb.clu.adj','m5.tAc.clu.adj','m5.tiA.clu.adj','m5.tiAb.clu.adj','m5.tiAc.clu.adj')
	#method.risk			<- c( 	'm5.tA.clu.adj','m5.tA.tp1.clu.adj','m5.tiA.tp1.clu.adj','m5.tA.tp2.clu.adj','m5.tiA.tp2.clu.adj','m5.tA.tp3.clu.adj','m5.tiA.tp3.clu.adj','m5.tA.tp4.clu.adj','m5.tiA.tp4.clu.adj')
	#method.risk			<- c( 	'm5.tAc.clu.adj','m5.tAc.tp1.clu.adj','m5.tAc.tp2.clu.adj','m5.tAc.tp3.clu.adj','m5.tAc.tp4.clu.adj')
	#method.risk			<- c( 	'm2Bwmx.tp1.clu.adj', 'm2Cwmx.tp1.clu.adj' )
	method.risk			<- c( 	'm2Cwmx.tp1.clu.adj','m2Cwmx.tp2.clu.adj','m2Cwmx.tp3.clu.adj','m2Cwmx.tp4.clu.adj','m2Cwmx.tp5.clu.adj','m2Cwmx.tp6.clu.adj' )
	#method.risk			<- c( 	'm2Cwmx.tp1.clu.adj' )
	dummy	<- sapply(method.risk, function(x)
			{
				cmd	<- hivc.cmd.props.estimate(indir, infile, insignat, indircov, infilecov, infiletree, infilexml.opt, infilexml.template, method, method.nodectime, x, method.recentctime, method.PDT, method.Acute, method.use.AcuteSpec, method.minQLowerU, method.lRNA.supp, method.thresh.pcoal, method.minLowerUWithNegT, outdir=outdir, outfile=outfile, resume=1, verbose=1)
				cat(cmd)
				#stop()
				#cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=3, hpc.mem="1800mb")
				#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=22, hpc.mem="1900mb")
				#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q='pqeph', hpc.nproc=1, hpc.walltime=71, hpc.mem="3800mb")
				#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="4000mb")
				#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q='pqeph', hpc.nproc=1, hpc.walltime=71, hpc.mem="7800mb")
				#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="64000mb")
				cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="120000mb")
				outdir		<- paste(DATA,"tmp",sep='/')
				outfile		<- paste("beta.",strsplit(date(),split=' ')[[1]],collapse='_',sep='')					
				hivc.cmd.hpccaller(outdir, outfile, cmd)			
			})	
}
######################################################################################
hivc.pipeline.various<- function()
{
	dir.name	<- DATA		 
	signat.in	<- "Wed_May__1_17/08/15_2013"
	signat.out	<- "Wed_May__1_17/08/15_2013"		
	

	if(0)	#align sequences in fasta file with Clustalo
	{
		indir	<- paste(dir.name,"tmp",sep='/')
		outdir	<- paste(dir.name,"tmp",sep='/')
		infile	<- "ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_Wed_May__1_17:08:15_2013.fasta"
		infile	<- "ATHENA_2013_03_All+LANL_Sequences_Sat_Jun_16_17:23:46_2013.fasta"
		cmd		<- hivc.cmd.clustalo(indir, infile, signat='', outdir=outdir)
		cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=1)
	}	
	if(0)	#extract first sequences for each patient as available from ATHENA data set
	{
		indir		<- paste(dir.name,"tmp",sep='/')
		infile		<- "ATHENA_2013_03_SeqMaster.R"		
		outdir		<- paste(dir.name,"tmp",sep='/')
		cmd			<- ''
		cmd			<- paste(cmd,hivc.cmd.get.firstseq(indir, infile, signat.in, signat.out, outdir=outdir),sep='')
	}
	if(0)	#compute branch length distances between tips in phylogenetic trees
	{
		indir		<- paste(DATA,"tmp/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlout_Wed_Dec_18_11:37:00_2013",sep='/')
		infiles		<- list.files(indir)
		infiles		<- infiles[ grepl('^ExaML_result.*finaltree\\.[0-9]{3}$',infiles)  ]	
		if(!length(infiles))	stop('cannot find files matching criteria')				
		cmd		<- sapply( seq_along(infiles), function(i)
			{
				infile		<- infiles[i]
				file		<- paste(indir, '/', infile, sep='')
				hivc.cmd.ph.dist.tips(file)					
			})	
		#put 50 into one job
		n		<- 50
		dummy<- lapply( seq_len(length(cmd)/n), function(i)
				{					
					tmp			<- paste( cmd[ seq.int((i-1)*n+1, min(i*n, length(cmd))) ], collapse='' )
					tmp			<- hivc.cmd.hpcwrapper(tmp, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="31400mb")
					signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outdir		<- indir
					outfile		<- paste("dtp",signat,sep='.')					
					cat(tmp)
					hivc.cmd.hpccaller(outdir, outfile, tmp)
					stop()			
				})
		stop()
	}
	if(0)	#compute raw pairwise genetic distances accounting correctly for ambiguous IUPAC nucleotides 
	{				
		gd.max		<- 0.045
		indir		<- paste(DATA,"tmp",sep='/')
		#infile		<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
		#insignat	<- "Sat_May_11_14/23/46_2013"				
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		cmd			<- hivc.cmd.get.geneticdist(indir, infile, insignat, gd.max, outdir=indir)
		
		outfile		<- paste("pipeline",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')		
		cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q="pqeph")
		cat(cmd)
		hivc.cmd.hpccaller(outdir, outfile, cmd)
		quit("no")
	}	
	if(0)	#run ExaML
	{
		project.hivc.examl()
		quit("no")
	}
	if(1)
	{
		hivc.pipeline.BEAST()
		quit("no")
	}
	if(0)	#run fisher similar analysis
	{
		indir		<- paste(dir.name,"tmp",sep='/')
		infile		<- "ATHENA_2013_03_SeqMaster.R"		
		outdir		<- paste(dir.name,"tmp",sep='/')
		cmd			<- ''
		cmd			<- paste(cmd,hivc.cmd.get.firstseq(indir, infile, signat.in, signat.out, outdir=outdir),sep='')
	}
	if(0)
	{
		outdir		<- dir.name		
		cmd			<- paste(CODE.HOME,"misc/hivclu.startme.R -exe=BETAREG.NUMBERS\n",sep='/')
		cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="149000mb")
	}
	if(0)
	{
		project.Gates()
		quit("no")
	}
	if(0)
	{
		project.hivc.examl.median.brl()
		quit("no")
	}						
}



