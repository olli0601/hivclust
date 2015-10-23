
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
		#thresh.brl			<- 0.096
		thresh.brl			<- 1000
		thresh.bs			<- 0.8
		#thresh.bs			<- 0.85
		#thresh.bs			<- 0.9
		#thresh.bs			<- 0.95
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
	if(0)
	{
		#indircov			<- paste(DATA,"derived",sep='/')
		#infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
		#
		#indir				<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/beast/beast_131011"		
		#infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		#insignat			<- "Tue_Aug_26_09:13:47_2013"				
		#infilexml.template	<- "um232rhU2045"
		#infilexml.opt		<- "mph4clutx4tip"
		#burnin				<- 2e7
		#
		indir				<- paste(DATA,"tmp70",sep='/')
		#indir				<- paste(DATA,"zip",sep='/')		
		infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"		
		insignat			<- "Wed_Dec_18_11:37:00_2013"		
		infilexml.template	<- "sasky_sdr06fr"
		#infilexml.opt		<- "alrh160"
		infilexml.opt		<- "clrh80"
		burnin				<- 1e7
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
		files		<- files[ grepl(infile, files, fixed=1) & grepl(gsub('/',':',insignat), files, fixed=1) & grepl(paste('_',infilexml.opt,'_',sep=''), files, fixed=1) & grepl(paste('_',infilexml.template,'_',sep=''), files, fixed=1) & grepl('_pool_[0-9]+',files) & grepl('trees$',files) ]				
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
					#argv			<<- unlist(strsplit(cmd,' '))
					#hivc.prog.BEAST2.get.cluster.trees()
					#stop()
					cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeelab", hpc.nproc=1, hpc.walltime=31, hpc.mem="6200mb")
					outdir		<- paste(DATA,"tmp",sep='/')
					outfile		<- paste("b2p.",strsplit(date(),split=' ')[[1]],collapse='_',sep='')					
					hivc.cmd.hpccaller(outdir, outfile, cmd)			
				})
		stop()
	}
	if(1)
	{
		indir				<- paste(DATA,"tmp70",sep='/')
		indircov			<- paste(DATA,"derived",sep='/')
		outdir				<- indir
		infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"		
		infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"		
		insignat			<- "Wed_Dec_18_11:37:00_2013"		
		outsignat			<- insignat
		#infilexml.template	<- "sasky_sdr06"		
		#infilexml.opt		<- "alsu50"
		infilexml.template	<- "sasky_sdr06fr"
		#infilexml.opt		<- "alrh160"
		infilexml.opt		<- "clrh80_bs0.7_brl1000"
		#infilexml.template	<- "um192rhU2080"
		#infilexml.opt		<- "mph4clutx4tip"	
		
		#indir					<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/beast/beast2_140201'
		#infile					<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_seroneg-130"
		#insignat				<- "Tue_Aug_26_09:13:47_2013"
		#infilexml.opt			<- "rsu815"
		#infilexml.template		<- "sasky_sdr06"
		
		files		<- list.files(indir)
		files		<- files[ grepl(infile, files, fixed=1) & grepl(gsub('/',':',insignat), files, fixed=1) & grepl(paste('_',infilexml.opt,'_',sep=''), files, fixed=1) & grepl(paste('_',infilexml.template,'_',sep=''), files, fixed=1) & grepl('_pool_[0-9]+',files) & grepl('_clu_[0-9]+',files) & grepl('R$',files) ]				
		if(!length(files))	stop('no input files matching criteria')
		tmp			<- regmatches( files, regexpr('_clu_[0-9]+',files)) 
		cluster		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
		file.info	<- data.table(file=files, cluster=cluster)
		setkey(file.info, cluster)		
		#80
		#file.info	<- subset(file.info, cluster%in%c(23, 77, 126, 152))
		#85
		#file.info	<- subset(file.info, cluster%in%c(462, 567, 1462))
		#90
		#file.info	<- subset(file.info, cluster%in%c(74,  221,  457, 1445))
		#95
		#file.info	<- subset(file.info, cluster%in%c(143,  210,  256,  330,  357,  365,  438,  444,  525,  721,  727,  734,  764,  768, 1146, 1190, 1194, 1211, 1371, 1393))
	
		#print(file.info)
		
		dummy		<- sapply( file.info[,unique(cluster)], function(clu)
				{
					cmd			<- hivc.cmd.beast2.processclustertrees(indir, infile, insignat, infilexml.opt, infilexml.template, cluster=clu, verbose=1, resume=1)					
					cat(cmd)
					cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeelab", hpc.nproc=1, hpc.walltime=300, hpc.mem="32000mb")
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
	if(0)	#	iterate over first large run
	{
		df.method				<- list(
				#	infection times
				as.data.table(expand.grid( DUMMY=1, method.minQLowerU=c(0.194, 0.148, 0.109), method.thresh.pcoal=0.2, method.thresh.bs=0.8, method.cut.brl=Inf )), 
				#	phylo exclusion criteria
				as.data.table(expand.grid( DUMMY=1, method.minQLowerU=c(0.148), method.thresh.pcoal=c(0.1,0.2,0.3), method.thresh.bs=c(0.7,0.8,0.85), method.cut.brl=Inf )),
				#	branch lengths
				as.data.table(expand.grid( DUMMY=1, method.minQLowerU=c(0.148), method.thresh.pcoal=c(0.2), method.thresh.bs=c(0.8), method.cut.brl=c(0.02, 0.04) )) 
		)
		df.method				<- do.call('rbind', df.method)
		df.method				<- df.method[-c(1,2,3,13,14,11,12),]
		tmp						<- data.table(DUMMY=1, method.lRNA.supp=100, method.use.AcuteSpec=1, method.PDT='SEQ', 	method.Acute='higher', method.minLowerUWithNegT=1)
		df.method				<- merge(df.method, tmp, by='DUMMY')
		df.method[, method.risk:='m2Awmx.tp1.clu.adj']
		df.method[, DUMMY:=seq_len(nrow(df.method))]
		
		df.method[,{
					cmd			<- hivc.cmd.props.estimate(indir, infile, insignat, indircov, infilecov, infiletree, infilexml.opt, infilexml.template, method, method.nodectime, method.risk, method.recentctime, method.PDT, method.Acute, method.use.AcuteSpec, method.minQLowerU, method.lRNA.supp, method.thresh.pcoal, method.minLowerUWithNegT, method.cut.brl, method.thresh.bs, outdir=outdir, outfile=outfile, resume=1, verbose=1)
					cat(cmd)	
					cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="130000mb")
					hivc.cmd.hpccaller(paste(DATA,"tmp",sep='/'), paste("beta.",paste(strsplit(date(),split=' ')[[1]],collapse='_'),sep=''), cmd)
				}, by='DUMMY']		
	}
	if(0)	#	iterate over second large run: bootstrap censoring
	{
		df.method				<- list(
				#	infection times
				as.data.table(expand.grid( DUMMY=1, method.minQLowerU=c(0.194, 0.109), method.thresh.pcoal=0.2, method.thresh.bs=0.8, method.cut.brl=Inf )), 
				#	phylo exclusion criteria
				as.data.table(expand.grid( DUMMY=1, method.minQLowerU=c(0.148), method.thresh.pcoal=c(0.1,0.2,0.3), method.thresh.bs=c(0.7,0.8,0.85), method.cut.brl=Inf )),
				#	branch lengths
				as.data.table(expand.grid( DUMMY=1, method.minQLowerU=c(0.148), method.thresh.pcoal=c(0.2), method.thresh.bs=c(0.8), method.cut.brl=c(0.02, 0.04) )) 
		)
		df.method				<- do.call('rbind', df.method)				
		tmp						<- data.table(method.risk=paste('m2Awmx.tp',1:6,'.clu.adj', sep=''), DUMMY=1, method.lRNA.supp=100, method.use.AcuteSpec=1, method.PDT='SEQ', 	method.Acute='higher', method.minLowerUWithNegT=1)
		df.method				<- merge(df.method, tmp, by='DUMMY', allow.cartesian=TRUE)
		df.method[, DUMMY:=seq_len(nrow(df.method))]		
		df.method[,{
					cmd			<- hivc.cmd.props.estimate(indir, infile, insignat, indircov, infilecov, infiletree, infilexml.opt, infilexml.template, method, method.nodectime, method.risk, method.recentctime, method.PDT, method.Acute, method.use.AcuteSpec, method.minQLowerU, method.lRNA.supp, method.thresh.pcoal, method.minLowerUWithNegT, method.cut.brl, method.thresh.bs, outdir=outdir, outfile=outfile, resume=1, verbose=1)
					cat(cmd)	
					cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="95000mb")
					hivc.cmd.hpccaller(paste(DATA,"tmp",sep='/'), paste("beta.",paste(strsplit(date(),split=' ')[[1]],collapse='_'),sep=''), cmd)
				}, by='DUMMY']		
	}
	if(0)	#	iterate over short runs: prop estimates
	{
		df.method				<- list(
				#	infection times
				as.data.table(expand.grid( DUMMY=1, method.minQLowerU=c(0.194, 0.109), method.thresh.pcoal=0.2, method.thresh.bs=0.8, method.cut.brl=Inf )), 
				#	phylo exclusion criteria
				as.data.table(expand.grid( DUMMY=1, method.minQLowerU=c(0.148), method.thresh.pcoal=c(0.1,0.2,0.3), method.thresh.bs=c(0.7,0.8,0.85), method.cut.brl=Inf )),
				#	branch lengths
				as.data.table(expand.grid( DUMMY=1, method.minQLowerU=c(0.148), method.thresh.pcoal=c(0.2), method.thresh.bs=c(0.8), method.cut.brl=c(0.02, 0.04) )) 
		)
		df.method				<- do.call('rbind', df.method)
		df.method				<- df.method[7,]
		df.method[, DUMMY2:=seq_len(nrow(df.method))]
		tmp						<- data.table(method.risk=paste('m2Awmx.wtn.tp',1:6, sep=''), DUMMY=1, method.lRNA.supp=100, method.use.AcuteSpec=1, method.PDT='SEQ', 	method.Acute='higher', method.minLowerUWithNegT=1)
		tmp						<- data.table(method.risk=c(paste('m2Awmx.tp',1:6, sep=''),paste('m2Awmx.nophyloscore.tp',1:6, sep=''),paste('m2Awmx.noscore.tp',1:6, sep='')), DUMMY=1, method.lRNA.supp=100, method.use.AcuteSpec=1, method.PDT='SEQ', 	method.Acute='higher', method.minLowerUWithNegT=1)
		#tmp					<- data.table(method.risk=paste('m2Cwmx.nophyloscore.tp',1:6, sep=''), DUMMY=1, method.lRNA.supp=100, method.use.AcuteSpec=1, method.PDT='SEQ', 	method.Acute='higher', method.minLowerUWithNegT=1)
		#tmp					<- data.table(method.risk=paste('m2Cwmx.noscore.tp',1:6, sep=''), DUMMY=1, method.lRNA.supp=100, method.use.AcuteSpec=1, method.PDT='SEQ', 	method.Acute='higher', method.minLowerUWithNegT=1)		
		df.method				<- merge(df.method, tmp, by='DUMMY', allow.cartesian=TRUE)		
		df.method				<- df.method[,{
					cmd			<- hivc.cmd.props.estimate(indir, infile, insignat, indircov, infilecov, infiletree, infilexml.opt, infilexml.template, method, method.nodectime, method.risk, method.recentctime, method.PDT, method.Acute, method.use.AcuteSpec, method.minQLowerU, method.lRNA.supp, method.thresh.pcoal, method.minLowerUWithNegT, method.cut.brl, method.thresh.bs, outdir=outdir, outfile=outfile, resume=1, verbose=1)
					list(CMD=cmd)
				}, by=c('DUMMY','DUMMY2')]	
		df.method[,{
					cmd			<- paste(CMD, collapse='', sep='')
					cat(cmd)					
					cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q='pqeelab', hpc.nproc=1, hpc.walltime=10, hpc.mem="4000mb")
					hivc.cmd.hpccaller(paste(DATA,"tmp",sep='/'), paste("beta.",paste(strsplit(date(),split=' ')[[1]],collapse='_'),sep=''), cmd)					
				}, by=c('DUMMY2')]	
	}
	if(1)	#	iterate over short runs: method.realloc
	{
		#testing.cov				<- seq(30,70,10)
		testing.cov				<- c(18,seq(30,70,20))
		#prep.cov				<- c(33,50,66)
		prep.cov				<- 50
		df.method				<- list(
				#	infection times
				as.data.table(expand.grid( DUMMY=1, method.minQLowerU=c(0.194, 0.109), method.thresh.pcoal=0.2, method.thresh.bs=0.8, method.cut.brl=Inf )), 
				#	phylo exclusion criteria
				as.data.table(expand.grid( DUMMY=1, method.minQLowerU=c(0.148), method.thresh.pcoal=c(0.1,0.2,0.3), method.thresh.bs=c(0.7,0.8,0.85), method.cut.brl=Inf )),
				#	branch lengths
				as.data.table(expand.grid( DUMMY=1, method.minQLowerU=c(0.148), method.thresh.pcoal=c(0.2), method.thresh.bs=c(0.8), method.cut.brl=c(0.02, 0.04) )) 
		)
		df.method				<- do.call('rbind', df.method)				
		#df.method				<- df.method[c(12,13),]
		#df.method				<- df.method[-7,]
		df.method				<- df.method[7,]
		if(1)
			df.method[, method.risk:='m2Awmx.wtn.tp4']
		if(0)
			df.method			<- merge(df.method, data.table(DUMMY=1, method.risk=c('m2Awmx.wtn.tp4','m2Awmx.noscore.tp4','m2Awmx.nophyloscore.tp4')), by='DUMMY')
		if(0)
			tmp					<- data.table(	DUMMY=1, method.lRNA.supp=100, method.use.AcuteSpec=1, method.PDT='SEQ', 	method.Acute='higher', method.minLowerUWithNegT=1,												
												method.realloc=  c(	'ImmediateART','ARTat500',	
							#	testing 12m
							paste('TestC12m',testing.cov,'pc',sep=''), paste('TestA12m',testing.cov,'pc',sep=''),
							#	testing 6m
							paste('TestC6m',testing.cov,'pc',sep=''), paste('TestA6m',testing.cov,'pc',sep=''),																	
							#	test + PrEP 
							paste('PrestIPrEXC12m',testing.cov,'pc',testing.cov,'pc',sep=''), paste('PrestPROUDC12m',testing.cov,'pc',testing.cov,'pc',sep=''),																	
							#	Test + PrEP + Immediate ART
							paste('PrestIPrEXC12m',testing.cov,'pc',testing.cov,'pc+ImmediateART',sep=''), paste('PrestPROUDC12m',testing.cov,'pc',testing.cov,'pc+ImmediateART',sep=''),																	
							#	Test + PrEP + ARTat500
							paste('PrestIPrEXC12m',testing.cov,'pc',testing.cov,'pc+ARTat500',sep=''), paste('PrestPROUDC12m',testing.cov,'pc',testing.cov,'pc+ARTat500',sep=''),
							#	test + treat
							paste('TestC12m',testing.cov,'pc+ImmediateART',sep=''), paste('TestC12m',testing.cov,'pc+ARTat500',sep=''),
							#	test (Acute)  + treat		
							paste('TestA12m',seq(30,70,10),'pc+ImmediateART',sep=''), paste('TestA12m',seq(30,70,10),'pc+ARTat500',sep='')
					))
		if(1)
		{
			tmp					<- as.vector(sapply(testing.cov, function(x)
										{
											c(	paste('TestC12m',x,'pc',sep=''), paste('TestA12m',x,'pc',sep=''), 
													#	test + PrEP 
													paste('PrestIPrEXC12m',x,'pc',prep.cov,'pc',sep=''), paste('PrestPROUDC12m',x,'pc',prep.cov,'pc',sep=''),																	
													#	Test + PrEP + Immediate ART
													paste('PrestIPrEXC12m',x,'pc',prep.cov,'pc+ImmediateART',sep=''), paste('PrestPROUDC12m',x,'pc',prep.cov,'pc+ImmediateART',sep=''),																	
													#	Test + PrEP + ARTat500
													paste('PrestIPrEXC12m',x,'pc',prep.cov,'pc+ARTat500',sep=''), paste('PrestPROUDC12m',x,'pc',prep.cov,'pc+ARTat500',sep=''),							
													#	test + Targeted PrEP 
													paste('TrestIPrEXC12m',x,'pc',prep.cov,'pc30y',sep=''), paste('TrestPROUDC12m',x,'pc',prep.cov,'pc30y',sep=''),																	
													#	Test + Targeted PrEP + Immediate ART
													paste('TrestIPrEXC12m',x,'pc',prep.cov,'pc30y+ImmediateART',sep=''), paste('TrestPROUDC12m',x,'pc',prep.cov,'pc30y+ImmediateART',sep=''),																	
													#	Test + Targeted PrEP + ARTat500
													paste('TrestIPrEXC12m',x,'pc',prep.cov,'pc30y+ARTat500',sep=''), paste('TrestPROUDC12m',x,'pc',prep.cov,'pc30y+ARTat500',sep=''),
													#	test + treat
													paste('TestC12m',x,'pc+ImmediateART',sep=''), paste('TestC12m',x,'pc+ARTat500',sep=''),
													#	test (Acute)  + treat		
													paste('TestA12m',x,'pc+ImmediateART',sep=''), paste('TestA12m',x,'pc+ARTat500',sep='')
											)
										}))
			tmp					<- data.table(	DUMMY=1, method.lRNA.supp=100, method.use.AcuteSpec=1, method.PDT='SEQ', 	method.Acute='higher', method.minLowerUWithNegT=1,	method.realloc=tmp)
		}					
		if(0)
			tmp					<- data.table(	DUMMY=1, method.lRNA.supp=100, method.use.AcuteSpec=1, method.PDT='SEQ', 	method.Acute='higher', method.minLowerUWithNegT=1,												
					method.realloc=  c(	#	test (Acute)  + treat		
							paste('TestA12m',testing.cov,'pc+ImmediateART',sep=''), paste('TestA12m',testing.cov,'pc+ARTat500',sep='')
					)
			)		
		df.method				<- merge(df.method, tmp, by='DUMMY', allow.cartesian=TRUE)
		df.method[, DUMMY:=seq_len(nrow(df.method))]
		#df.method				<- df.method[c(13:22,30:32,40:42,50:52),]
		df.method[, DUMMY2:=rep(1:(nrow(df.method)%/%3+1), each=3)[seq_len(nrow(df.method))]]
		#df.method[, DUMMY2:=seq_len(nrow(df.method))]
		df.method				<- df.method[,{
					cmd			<- hivc.cmd.props.estimate(indir, infile, insignat, indircov, infilecov, infiletree, infilexml.opt, infilexml.template, method, method.nodectime, method.risk, method.recentctime, method.PDT, method.Acute, method.use.AcuteSpec, method.minQLowerU, method.lRNA.supp, method.thresh.pcoal, method.minLowerUWithNegT, method.cut.brl, method.thresh.bs, method.realloc=method.realloc, outdir=outdir, outfile=outfile, resume=1, verbose=1)
					list(CMD=cmd)
				}, by=c('DUMMY','DUMMY2')]	
		df.method[,{
					cmd			<- paste(CMD, collapse='', sep='')
					cat(cmd)						
					cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q='pqeelab', hpc.nproc=1, hpc.walltime=10, hpc.mem="4000mb")
					hivc.cmd.hpccaller(paste(DATA,"tmp",sep='/'), paste("beta.",paste(strsplit(date(),split=' ')[[1]],collapse='_'),sep=''), cmd)
					#stop()
				}, by=c('DUMMY2')]	
	}
	if(0)	#custom
	{
		method.lRNA.supp			<- 100
		method.use.AcuteSpec		<- 1
		method.minQLowerU			<- 0.148	#0.194 0.148 0.109
		method.thresh.pcoal			<- 0.2
		method.PDT					<- 'SEQ'
		#method.Acute				<- 'empirical'
		method.Acute				<- 'higher'
		#method.Acute				<- 'lower'
		#method.Acute				<- 'central'
		method.minLowerUWithNegT	<- 1
		method.cut.brl				<- Inf
		method.thresh.bs			<- 0.8
		
		#	basic m2Bwmx runs	mem=1800 MB	
		#method.risk		<- c(	'm2Cwmx.wtn.tp1','m2Cwmx.wtn.tp2','m2Cwmx.wtn.tp3','m2Cwmx.wtn.tp4','m2Cwmx.wtn.tp5','m2Cwmx.wtn.tp6')
		method.risk		<- c(	'm2Awmx.wtn.tp1','m2Awmx.wtn.tp2','m2Awmx.wtn.tp3','m2Awmx.wtn.tp4','m2Awmx.wtn.tp5','m2Awmx.wtn.tp6',
				'm2Bwmx.wtn.tp1','m2Bwmx.wtn.tp2','m2Bwmx.wtn.tp3','m2Bwmx.wtn.tp4','m2Bwmx.wtn.tp5','m2Bwmx.wtn.tp6',
				'm2Cwmx.wtn.tp1','m2Cwmx.wtn.tp2','m2Cwmx.wtn.tp3','m2Cwmx.wtn.tp4','m2Cwmx.wtn.tp5','m2Cwmx.wtn.tp6')
		#method.risk		<- c(	'm2Cwmx.tp1','m2Cwmx.tp2','m2Cwmx.tp3','m2Cwmx.tp4','m2Cwmx.tp5','m2Cwmx.tp6',
		#		'm2Cwmx.nophyloscore.tp1','m2Cwmx.nophyloscore.tp2','m2Cwmx.nophyloscore.tp3','m2Cwmx.nophyloscore.tp4','m2Cwmx.nophyloscore.tp5','m2Cwmx.nophyloscore.tp6',
		#		'm2Cwmx.noscore.tp1','m2Cwmx.noscore.tp2','m2Cwmx.noscore.tp3','m2Cwmx.noscore.tp4','m2Cwmx.noscore.tp5','m2Cwmx.noscore.tp6')
		#	m5 runs
		#method.risk		<- c(	'm5.tA.tp1.clu','m5.tA.tp2.clu','m5.tA.tp3.clu','m5.tA.tp4.clu','m5.tA.tp1.clu.wstar','m5.tA.tp2.clu.wstar','m5.tA.tp3.clu.wstar','m5.tA.tp4.clu.wstar')
		#!##method.risk		<- c(	'm5.tAc.tp1','m5.tAc.tp2','m5.tAc.tp3','m5.tAc.tp4','m5.tAc.tp1.wstar','m5.tAc.tp2.wstar','m5.tAc.tp3.wstar','m5.tAc.tp4.wstar')
		# use to pre-compute tables mem 95 GB
		method.risk			<- c( 	'm2Awmx.tp1.clu.adj','m2Awmx.tp2.clu.adj','m2Awmx.tp3.clu.adj','m2Awmx.tp4.clu.adj','m2Awmx.tp5.clu.adj','m2Awmx.tp6.clu.adj',
				'm2Bwmx.tp1.clu.adj','m2Bwmx.tp2.clu.adj','m2Bwmx.tp3.clu.adj','m2Bwmx.tp4.clu.adj','m2Bwmx.tp5.clu.adj','m2Bwmx.tp6.clu.adj',
				'm2Cwmx.tp1.clu.adj','m2Cwmx.tp2.clu.adj','m2Cwmx.tp3.clu.adj','m2Cwmx.tp4.clu.adj','m2Cwmx.tp5.clu.adj','m2Cwmx.tp6.clu.adj')
		# use to pre-compute tables mem 130 GB
		method.risk				<- c( 	'm2Cwmx.tp1.clu.adj' )	
		#method.risk				<- 'm2Cwmx.wtn.tp4'
		dummy	<- sapply(method.risk, function(x)
				{
					cmd	<- hivc.cmd.props.estimate(indir, infile, insignat, indircov, infilecov, infiletree, infilexml.opt, infilexml.template, method, method.nodectime, x, method.recentctime, method.PDT, method.Acute, method.use.AcuteSpec, method.minQLowerU, method.lRNA.supp, method.thresh.pcoal, method.minLowerUWithNegT, method.cut.brl, method.thresh.bs, outdir=outdir, outfile=outfile, resume=1, verbose=1)
					cat(cmd)
					#stop()
					#cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=3, hpc.mem="1800mb")
					#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=20, hpc.mem="1900mb")
					cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q='pqeelab', hpc.nproc=1, hpc.walltime=71, hpc.mem="4000mb")
					#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="4000mb")
					#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q='pqeph', hpc.nproc=1, hpc.walltime=71, hpc.mem="7800mb")
					#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="95000mb")
					#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="130000mb")
					outdir		<- paste(DATA,"tmp",sep='/')
					outfile		<- paste("beta.",strsplit(date(),split=' ')[[1]],collapse='_',sep='')					
					hivc.cmd.hpccaller(outdir, outfile, cmd)			
				})	
	}
}
######################################################################################
hivc.pipeline.ages<- function()
{
	#stop()
	indir					<- paste(DATA,"fisheretal_data",sep='/')		
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	outdir					<- paste(DATA,"tpairs_age",sep='/')
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"
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
	if(1)	#custom
	{
		method.lRNA.supp			<- 100
		method.use.AcuteSpec		<- 1
		method.minQLowerU			<- 0.148	#0.194 0.148 0.109
		method.thresh.pcoal			<- 0.2
		method.PDT					<- 'SEQ'
		#method.Acute				<- 'empirical'
		method.Acute				<- 'higher'
		#method.Acute				<- 'lower'
		#method.Acute				<- 'central'
		method.minLowerUWithNegT	<- 1
		method.cut.brl				<- Inf
		method.thresh.bs			<- 0.8
		
		# use to pre-compute tables mem 130 GB		
		#method.risk				<- c( 	'm5A.tp1.clu.adj' )
		#method.risk				<- c( 	'm5B.tp1.clu.adj' )
		#method.risk				<- c( 	'm5C.tp1.clu.adj' )
		method.risk				<- c( 	'm5D.tp1.clu.adj' )
		#method.risk				<- c( 	'm5E.tp1.clu.adj' )
		#method.risk				<- c( 	'm5F.tp1.clu.adj' )
		dummy	<- sapply(method.risk, function(x)
				{
					cmd	<- hivc.cmd.age.estimate(indir, infile, insignat, indircov, infilecov, infiletree, infilexml.opt, infilexml.template, method, method.nodectime, x, method.recentctime, method.PDT, method.Acute, method.use.AcuteSpec, method.minQLowerU, method.lRNA.supp, method.thresh.pcoal, method.minLowerUWithNegT, method.cut.brl, method.thresh.bs, outdir=outdir, outfile=outfile, resume=1, verbose=1)
					cat(cmd)
					#stop()
					#cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=3, hpc.mem="1800mb")
					#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=20, hpc.mem="1900mb")
					#cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q='pqeelab', hpc.nproc=1, hpc.walltime=71, hpc.mem="4000mb")
					#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="4000mb")
					#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q='pqeph', hpc.nproc=1, hpc.walltime=71, hpc.mem="7800mb")
					#cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="95000mb")
					cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q=NA, hpc.nproc=1, hpc.walltime=71, hpc.mem="130000mb")
					outdir		<- paste(DATA,"tmp",sep='/')
					outfile		<- paste("ags.",gsub(' ','_',date()),sep='')					
					hivc.cmd.hpccaller(outdir, outfile, cmd)			
				})	
	}
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
		indir		<- paste(DATA,"ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlout",sep='/')
		infiles		<- list.files(indir)
		infiles		<- infiles[ grepl('^ExaML_result.*finaltree\\.[0-9]{3}$',infiles)  ]	
		if(!length(infiles))	stop('cannot find files matching criteria')				
		cmd		<- sapply( seq_along(infiles), function(i)
			{
				infile		<- infiles[i]
				file		<- paste(indir, '/', infile, sep='')
				hivc.cmd.ph.dist.tips(file)					
			})	
		#put 1 into one job
		n		<- 1
		dummy<- lapply( seq_len(length(cmd)/n), function(i)
				{					
					tmp			<- paste( cmd[ seq.int((i-1)*n+1, min(i*n, length(cmd))) ], collapse='' )
					tmp			<- hivc.cmd.hpcwrapper(tmp, hpc.q='pqeelab', hpc.nproc=1, hpc.walltime=71, hpc.mem="31400mb")
					signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outdir		<- indir
					outfile		<- paste("dtp",signat,sep='.')					
					cat(tmp)
					hivc.cmd.hpccaller(outdir, outfile, tmp)
					#stop()			
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
	if(0)
	{
		hivc.pipeline.BEAST()
		quit("no")
	}
	if(0)
	{
		hivc.pipeline.BEASTout()
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
	if(1)
	{
		#project.hivc.examl.median.brl()
		cmd			<- hivc.cmd.various()
		cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q=NA, hpc.walltime=71, hpc.mem="95000mb")
		cat(cmd)		
		outdir		<- paste(DATA,"tmp",sep='/')
		outfile		<- paste("vrs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		hivc.cmd.hpccaller(outdir, outfile, cmd)
		quit("no")		
	}	
	if(0)	#	run ExamML with partition for tree comparison
	{		
		indir.wgaps	<- '/Users/Oliver/git/HPTN071sim/treec150623/withgapstrees'
		indir.wgaps	<- '/work/or105/Gates_2014/tree_comparison'
		infiles		<- data.table(FILE=list.files(indir.wgaps, pattern='\\.R$'))
		infiles[,SIGNAT:=infiles[, regmatches(FILE,regexpr('BWC.*|UGC.*', FILE))]]
		tmp			<- infiles[, list(BASE=gsub(paste('_',SIGNAT,sep=''),'',FILE)), by='FILE']
		infiles		<- merge(infiles, tmp, by='FILE')
		set(infiles, NULL, 'SIGNAT', infiles[, gsub('\\.R','',SIGNAT)])
		infiles[, PARTITION:= gsub('\\.R','_codon.txt',FILE)]
		bs.from		<- 0
		bs.to		<- 19
		bs.n		<- 20
		outdir		<- indir.wgaps
		invisible(infiles[, {					
							infile		<- BASE
							signat.in	<- signat.out	<- SIGNAT
							args.parser	<- paste("-m DNA -q",PARTITION)
							cmd			<- hivc.cmd.examl.bootstrap(indir.wgaps, infile, signat.in, signat.out, bs.from=bs.from, bs.to=bs.to, prog.bscreate=PR.EXAML.BSCREATE, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, prog.examl=PR.EXAML.EXAML, opt.bootstrap.by="codon", args.parser=args.parser, args.examl="-f o -m GAMMA", prog.supportadder=PR.EXAML.BS, tmpdir.prefix="examl")					
							invisible(lapply(cmd, function(x)
											{												
												x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=21, hpc.q="pqeelab", hpc.mem="950mb", hpc.nproc=1)
												signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
												outfile	<- paste("ex",signat,sep='.')
												#cat(x)
												hivc.cmd.hpccaller(outdir, outfile, x)
												Sys.sleep(1)
											}))
							NULL					
						}, by='FILE'])				
	}
		
}



