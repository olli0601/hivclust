
######################################################################################
project.hivc.check.nt.table<- function(dir.name= DATA, verbose=1)
{
	#> subset(nt.table, Patient=='M14909')
    #Patient  risk      factor X.clu     X.msm  X.seq YX
 	#1:  M14909 stage ART.suA.N.1     0  8242  3855  0
 	#2:  M14909 stage ART.suA.Y.1     0     8    12  0
	lRNA.supp	<- log10(51)
	X.seq.p		<- subset(X.seq, t.Patient=='M14909')
	X.msm.p		<- subset(X.msm, t.Patient=='M14909')
	#	--> t.Patient is identical in X.msm and X.seq
	#	why is there a difference in the nt.table????
	setdiff( X.msm.p[, unique(Patient)], X.seq.p[, unique(Patient)] )
	setdiff( X.seq.p[, unique(Patient)], X.msm.p[, unique(Patient)] )
	
	subset(X.seq.p, stage=='ART.started' & lRNA<=lRNA.supp)
	
	subset(X.seq.p, stage=='ART.started')[,table(ART.F)]
	#ART.F
	#No  	Yes 
	#8627   46 
	subset(X.msm.p, stage=='ART.started')[,table(ART.F)]
	#ART.F
	#No  	Yes 
	#8627   46 
	X.seq.p2	<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.seq.p, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp, plot.file.varyvl=NA, plot.file.or=NA )
	X.msm.p2	<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.msm.p, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp, plot.file.varyvl=NA, plot.file.or=NA )
	#	stratify seems OK
	X.msm.p1	<- subset(X.msm, t.Patient=='M32608')
	X.seq.p1	<- subset(X.seq, t.Patient=='M32608')
	X.seq.p1	<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.seq.p1, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp, plot.file.varyvl=NA, plot.file.or=NA )
	X.msm.p1	<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.msm.p1, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp, plot.file.varyvl=NA, plot.file.or=NA )	
	X.seq.p1[, length(unique(Patient))]
	X.msm.p1[, length(unique(Patient))]
	X.seq.p1[, table(t.period)]
	X.msm.p1[, table(t.period)]	
	setdiff( subset(X.seq.p1, CD4b.tperiod=='ART.suA.Y.3')[, Patient], subset(X.msm.p1, CD4b.tperiod=='ART.suA.Y.3')[, Patient] )
	#	"M37323" "M37456" "M37459" "M37576" "M37609" "M38159" "M38870"
	subset(X.seq.p1, CD4b.tperiod=='ART.suA.Y.3' & Patient=="M37323")		#--> lRNA	
	subset(X.seq.p1, CD4b.tperiod=='ART.suA.Y.3' & Patient=="M37456")		#--> same lRNA
	X.seq.p1[, table(CD4b.tperiod)]
	X.msm.p1[, table(CD4b.tperiod)]
	#
	#
	X.msm.p1	<- subset(X.msm, t.Patient=='M36262')
	X.seq.p1	<- subset(X.seq, t.Patient=='M36262')
	X.seq.p1	<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.seq.p1, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp, plot.file.varyvl=NA, plot.file.or=NA )
	X.msm.p1	<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.msm.p1, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp, plot.file.varyvl=NA, plot.file.or=NA )	
	
	X.msm.p1	<- subset(X.msm, t.Patient=='M36567')
	X.seq.p1	<- subset(X.seq, t.Patient=='M36567')
	X.seq.p1	<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.seq.p1, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp, plot.file.varyvl=NA, plot.file.or=NA )
	X.msm.p1	<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.msm.p1, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp, plot.file.varyvl=NA, plot.file.or=NA )	
	setdiff( subset(X.seq.p1, CD4b.tperiod=='ART.suA.Y.4')[, Patient], subset(X.msm.p1, CD4b.tperiod=='ART.suA.Y.4')[, Patient] )
	# [1] "M38061" "M38146" "M38401" "M39012" "M39073" "M39096" "M39134" "M39143" "M39200" "M39228" "M39281" "M39309" "M39447"
	#[14] "M39707"
	subset(X.seq.p1, CD4b.tperiod=='ART.suA.Y.4' & Patient=="M38061")
	subset(X.msm.p1, CD4b.tperiod=='ART.suA.N.4' & Patient=="M38061")
}
######################################################################################
project.hivc.check.CD4interpolation<- function(dir.name= DATA, verbose=1)
{
	
	df.tpairs	<- subset(YX, !is.na(CD4), select=unique(t.Patient))
	X.incare	<- project.athena.Fisheretal.X.incare(df.tpairs, df.all, df.viro, df.immu, df.treatment, t.period=1/8, t.endctime=2013.)
	X.incare	<- subset(X.incare, !is.na(CD4), select=c(t.Patient, t, CD4))
	
	df.cd4		<- copy(df.immu)
	setnames(df.cd4, 'Patient','t.Patient')
	
	df.cd4	<- merge( data.table( t.Patient=X.incare[, unique(t.Patient)] ), df.cd4, by='t.Patient' )
	set(df.cd4, NULL, 'PosCD4', hivc.db.Date2numeric(df.cd4[,PosCD4]))
	tmp		<- subset(df.cd4[, list(n=length(CD4), cg500= all(CD4>500), cg350=all(CD4>350), csd=sd(CD4)), by='t.Patient'], n>=1, select=c(t.Patient, cg500, cg350, csd))
	tmp[, cgroup:='l350']
	set(tmp, tmp[, which(cg350)], 'cgroup', 'l500')
	set(tmp, tmp[, which(cg500)], 'cgroup', 'g500')
	setkey(tmp, cgroup, csd)
	tmp		<- tmp[,  list(t.Patient=t.Patient, csd= ceiling(seq_along(csd)/length(csd)*3)), by='cgroup']
	
	df.cd4		<- merge( tmp, df.cd4, by='t.Patient' )
	X.incare	<- merge( X.incare, tmp, by='t.Patient')
	
	set(X.incare, NULL, 't.Patient', X.incare[,factor(t.Patient)])
	set(X.incare, NULL, 'CD4', X.incare[,as.double(CD4)])
	set(df.cd4, NULL, 't.Patient', df.cd4[,factor(t.Patient)])
	#
	#	
	tmp			<- subset( X.incare, cgroup=='g500' )
	setnames(tmp, 't', 'PosCD4')
	tmp2		<- merge( df.cd4, unique(subset(tmp, select=t.Patient)), by='t.Patient' )	
	tmp3		<- lapply( tmp2[, unique(t.Patient)], function(x)
			{
				z		<- subset(tmp, t.Patient==x, PosCD4)
				z2		<- subset(tmp2, t.Patient==x)
				if(nrow(z2)>1)
				{
					cd4.d	<- ifelse(z2[, diff(range(PosCD4))<2], 1, min(15,ceiling(nrow(z2)/8))  )
					cd4.ml	<- gamlss(CD4 ~ PosCD4, data=z2, family='NO', trace = FALSE)
					cd4.m	<- gamlss(CD4 ~ bs(PosCD4, degree=cd4.d), data=z2, family='NO', trace = FALSE)
					suppressWarnings( cd4.sl	<- predict(cd4.ml, type='response', newdata=z, data=z2) )
					suppressWarnings( cd4.s	<- predict(cd4.m, type='response', newdata=z, data=z2) )											
					dev		<- deviance(cd4.ml)-deviance(cd4.m)
				}
				if(nrow(z2)==1)
				{
					cd4.s	<- cd4.sl	<- rep(z2[1,][,CD4],nrow(z))
					dev		<- NA_real_
				}
				data.table(t.Patient=x, t=z[,PosCD4], CD4sl=cd4.sl, CD4s=cd4.s, dev=dev)
			})
	setnames(tmp, 'PosCD4', 't')
	tmp			<- merge(tmp, do.call('rbind',tmp3), by=c('t.Patient','t'))	
	df.cd4.g500	<- melt(tmp, measure.var=c('CD4', 'CD4s', 'CD4sl'), variable.name='method', value.name='CD4')
	ggplot(df.cd4.g500, aes(x=t, y=CD4, group=interaction(t.Patient, method))) + geom_line(aes(colour=method)) + 
			geom_point(data=tmp2, aes(x=PosCD4)) +
			facet_grid(t.Patient ~ csd)
	file		<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_checkCD4interpolation_g500.pdf'
	ggsave(file=file, w=6, h=50, limitsize=FALSE)
	#
	#
	tmp			<- subset( X.incare, cgroup=='l500' )
	setnames(tmp, 't', 'PosCD4')
	tmp2		<- merge( df.cd4, unique(subset(tmp, select=t.Patient)), by='t.Patient' )	
	tmp3		<- lapply( tmp2[, unique(t.Patient)], function(x)
			{
				z		<- subset(tmp, t.Patient==x, PosCD4)
				z2		<- subset(tmp2, t.Patient==x)
				if(nrow(z2)>1)
				{
					cd4.d	<- ifelse(z2[, diff(range(PosCD4))<2], 1, min(15,ceiling(nrow(z2)/8))  )
					cd4.ml	<- gamlss(CD4 ~ PosCD4, data=z2, family='NO', trace = FALSE)
					cd4.m	<- gamlss(CD4 ~ bs(PosCD4, degree=cd4.d), data=z2, family='NO', trace = FALSE)
					suppressWarnings( cd4.sl	<- predict(cd4.ml, type='response', newdata=z, data=z2) )
					suppressWarnings( cd4.s	<- predict(cd4.m, type='response', newdata=z, data=z2) )											
					dev		<- deviance(cd4.ml)-deviance(cd4.m)
				}
				if(nrow(z2)==1)
				{
					cd4.s	<- cd4.sl	<- rep(z2[1,][,CD4],nrow(z))
					dev		<- NA_real_
				}
				data.table(t.Patient=x, t=z[,PosCD4], CD4sl=cd4.sl, CD4s=cd4.s, dev=dev)
			})
	setnames(tmp, 'PosCD4', 't')
	tmp			<- merge(tmp, do.call('rbind',tmp3), by=c('t.Patient','t'))	
	df.cd4.l500	<- melt(tmp, measure.var=c('CD4', 'CD4s', 'CD4sl'), variable.name='method', value.name='CD4')
	ggplot(df.cd4.l500, aes(x=t, y=CD4, group=interaction(t.Patient, method))) + geom_line(aes(colour=method)) + 
			geom_point(data=tmp2, aes(x=PosCD4)) +
			facet_grid(t.Patient ~ csd)
	file		<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_checkCD4interpolation_l500.pdf'
	ggsave(file=file, w=6, h=170, limitsize=FALSE)
	#
	#
	tmp			<- subset( X.incare, cgroup=='l350' )
	setnames(tmp, 't', 'PosCD4')
	tmp2		<- merge( df.cd4, unique(subset(tmp, select=t.Patient)), by='t.Patient' )	
	tmp3		<- lapply( tmp2[, unique(t.Patient)], function(x)
			{
				z		<- subset(tmp, t.Patient==x, PosCD4)
				z2		<- subset(tmp2, t.Patient==x)
				if(nrow(z2)>1)
				{
					cd4.d	<- ifelse(z2[, diff(range(PosCD4))<2], 1, min(15,ceiling(nrow(z2)/8))  )
					cd4.ml	<- gamlss(CD4 ~ PosCD4, data=z2, family='NO', trace = FALSE)
					cd4.m	<- gamlss(CD4 ~ bs(PosCD4, degree=cd4.d), data=z2, family='NO', trace = FALSE)
					suppressWarnings( cd4.sl	<- predict(cd4.ml, type='response', newdata=z, data=z2) )
					suppressWarnings( cd4.s	<- predict(cd4.m, type='response', newdata=z, data=z2) )											
					dev		<- deviance(cd4.ml)-deviance(cd4.m)
				}
				if(nrow(z2)==1)
				{
					cd4.s	<- cd4.sl	<- rep(z2[1,][,CD4],nrow(z))
					dev		<- NA_real_
				}
				data.table(t.Patient=x, t=z[,PosCD4], CD4sl=cd4.sl, CD4s=cd4.s, dev=dev)
			})
	setnames(tmp, 'PosCD4', 't')
	tmp			<- merge(tmp, do.call('rbind',tmp3), by=c('t.Patient','t'))	
	df.cd4.l350	<- melt(tmp, measure.var=c('CD4', 'CD4s', 'CD4sl'), variable.name='method', value.name='CD4')
	ggplot(df.cd4.l350, aes(x=t, y=CD4, group=interaction(t.Patient, method))) + geom_line(aes(colour=method)) + 
			geom_point(data=tmp2, aes(x=PosCD4)) +
			facet_grid(t.Patient ~ csd)
	file		<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_checkCD4interpolation_l350.pdf'
	ggsave(file=file, w=6, h=1100, limitsize=FALSE)
		
	#min deviance to select nonlinear prediction
	subset(df.cd4.l500, t.Patient%in%c('M11902','M19233','M25370','M27565','M3728'))[, min(dev)]
	#max deviance to select linear prediction
	subset(df.cd4.l500, t.Patient%in%c('M34587','M37014'))[, max(dev)]
	
	#
	#
	#
	tmp			<- subset( X.incare, cgroup=='l350' )		
	setnames(tmp, 't', 'PosCD4')
	tmp2		<- merge( df.cd4, unique(subset(tmp, select=t.Patient)), by='t.Patient' )	
	tmp3		<- lapply( tmp2[, unique(t.Patient)], function(x)
			{
				z		<- subset(tmp, t.Patient==x, PosCD4)
				z2		<- subset(tmp2, t.Patient==x)
				if(nrow(z2)>1)
				{
					cd4.d	<- ifelse(z2[, diff(range(PosCD4))<2], 1, min(15,ceiling(nrow(z2)/8))  )
					cd4.ml	<- gamlss(CD4 ~ PosCD4, data=z2, family='NO', trace = FALSE)
					cd4.m	<- gamlss(CD4 ~ bs(PosCD4, degree=cd4.d), data=z2, family='NO', trace = FALSE)
					if(deviance(cd4.ml)-deviance(cd4.m)>10)
						suppressWarnings( cd4.s	<- predict(cd4.m, type='response', newdata=z, data=z2) )
					if(deviance(cd4.ml)-deviance(cd4.m)<=10)
						suppressWarnings( cd4.s	<- predict(cd4.ml, type='response', newdata=z, data=z2) )					
				}
				if(nrow(z2)==1)
				{
					cd4.s	<- rep(z2[1,][,CD4],nrow(z))
				}
				data.table(t.Patient=x, t=z[,PosCD4], CD4s=cd4.s)
			})
	setnames(tmp, 'PosCD4', 't')
	tmp			<- merge(tmp, do.call('rbind',tmp3), by=c('t.Patient','t'))	
	tmp			<- melt(tmp, measure.var=c('CD4', 'CD4s'), variable.name='method', value.name='CD4')
	ggplot(tmp, aes(x=t, y=CD4, group=interaction(t.Patient, method))) + geom_line(aes(colour=method)) + 
			geom_point(data=tmp2, aes(x=PosCD4)) +
			facet_grid(t.Patient ~ csd)
	file		<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_checkCD4interpolation_l350used.pdf'
	ggsave(file=file, w=6, h=1100, limitsize=FALSE)
	
}

######################################################################################
project.Gates.RootSeqSim.getrootseq<- function()
{
	#dir.name	<- "/work/or105/Gates_2014"
	dir.name	<- '/Users/Oliver/duke/2014_Gates'  
	
	if(1)	#	devel
	{		
		indir				<- paste(dir.name,'methods_comparison_rootseqsim/140730',sep='/')
		outdir				<- indir
		infile.xml			<- 'working.xml'
		infile.beastparsed	<- 'working.R'
		outfile				<- 'working_ancseq.R'
		
		#	load BEAST PARSER output
		file		<- paste(indir, '/',infile.beastparsed, sep='')
		load(file)	#	expect tree, node.stat		
		#	get original sequences
		file		<- paste(indir, '/',infile.xml, sep='')
		bxml		<- xmlTreeParse(file, useInternalNodes=TRUE, addFinalizer = TRUE)
		bseq		<- hivc.beast.get.sequences(bxml, verbose=1)	
		bseq		<- merge(bseq, data.table(ALIGNMENT_ID=paste('alignment',1:3,sep=''), GENE=c('env','gag','pol')), by='ALIGNMENT_ID')				
		#	compute gag pol env		
		tmp			<- PANGEA.RootSeqSim.get.ancestral.seq(tree, node.stat, bseq, tree.id.sep='_', tree.id.idx.mcmcit=2, tree.id.burnin=1, label.sep='|', label.idx.ctime=2)
		ancseq.gag	<- tmp$GAG
		ancseq.env	<- tmp$ENV
		ancseq.pol	<- tmp$POL
		#	save as R
		file		<- paste(outdir, outfile, sep='/')
		save(ancseq.gag, ancseq.env, ancseq.pol, file=file)			
		#	save as FASTA
		file		<- paste(outdir, paste(substr(outfile,1,nchar(outfile)-1),'fasta',sep=''), sep='/')		
		write.dna(cbind(ancseq.gag, ancseq.env, ancseq.pol), file, format = "fasta")		
		#
		#	sample ancestral sequences between 1980-2000 and reconstruct tree with RAxML
		#
		ancseq				<- cbind(ancseq.gag, ancseq.env, ancseq.pol)
		label.sep			<- '|'
		label.idx.tree.id	<- 1
		label.idx.node.id	<- 2
		label.idx.ctime		<- 3
		ancseq.label		<- data.table(	TREE_ID= sapply(strsplit(rownames(ancseq), label.sep, fixed=1),'[[',label.idx.tree.id),
											NODE_ID= as.numeric(sapply(strsplit(rownames(ancseq), label.sep, fixed=1),'[[',label.idx.node.id)),
											CALENDAR_TIME= as.numeric(sapply(strsplit(rownames(ancseq), label.sep, fixed=1),'[[',label.idx.ctime)))
		hist( ancseq.label[, CALENDAR_TIME], breaks=100 )
	}	
}
######################################################################################
project.Gates.RootSeqSim.BEAST.SSAfg.checkancestralseq.runExaML<- function()
{
	DATA				<<- "/work/or105/Gates_2014"
	#DATA				<<- '/Users/Oliver/duke/2014_Gates'
	dir.name			<- DATA  	
	indir				<- paste(dir.name,'methods_comparison_rootseqsim/140813',sep='/')
	#	ExaML bootstrap args
	bs.from		<- 0
	bs.to		<- 100
	bs.n		<- 100
	
	#	search for 'checkdraw' files
	infiles		<- list.files(indir)
	infiles		<- infiles[ sapply(infiles, function(x) grepl('.*checkdraw[0-9]+.*R$',x) ) ]	
	if(!length(infiles))	stop('cannot find files matching criteria')
	
	outdir		<- indir
	for(i in seq_along(infiles))
	{
		infile		<- infiles[i]
		infile		<- substr(infile, 1, nchar(infile)-2)
		insignat	<- regmatches(infile, regexpr('checkdraw[0-9]+_.*', infile))
		insignat	<- regmatches(insignat,regexpr('_.*',insignat))
		insignat	<- substr(insignat,2,nchar(insignat))
		infile		<- regmatches(infile, regexpr('.*checkdraw[0-9]+', infile))
		
		
		cmd			<- hivc.cmd.examl.bootstrap(indir, infile, insignat, insignat, bs.from=bs.from, bs.to=bs.to,bs.n=bs.n,outdir=outdir, resume=1, verbose=1)
		dummy		<- lapply(cmd, function(x)
				{				
					x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=24, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
					#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=24, hpc.q="pqeph", hpc.mem="3850mb", hpc.nproc=8)
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("exa",signat,sep='.')
					#cat(x)
					hivc.cmd.hpccaller(outdir, outfile, x)
					Sys.sleep(1)
				})
		if(i==3)
			stop()
	}
	
}
######################################################################################
project.Gates.test.runxml<- function()
{
	DATA		<<- "/work/or105/Gates_2014"
	#DATA		<<- '/Users/Oliver/duke/2014_Gates'		
	if(1)
	{			
		indir		<- paste(DATA,'methods_comparison_pipeline/150130',sep='/')
		#search for XML files in indir
		infiles		<- list.files(indir, pattern=paste(".xml$",sep=''))
		insignat	<- ''	
		hpc.ncpu	<- 8		
		for(infile in infiles)
		{
			infile		<- substr(infile, 1, nchar(infile)-4) 		
			cmd			<- hivc.cmd.beast.runxml(indir, infile, insignat, prog.beast=PR.BEAST, prog.beast.opt=" -beagle -working", hpc.tmpdir.prefix="beast", hpc.ncpu=hpc.ncpu)
			cat(cmd)	
			cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q='pqeelab', hpc.nproc=hpc.ncpu, hpc.walltime=71, hpc.mem="6000mb")		
			outdir		<- indir
			outfile		<- paste("bpg.",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='')					
			hivc.cmd.hpccaller(outdir, outfile, cmd)		
		}
	}
}
######################################################################################
project.Gates.RootSeqSim.runxml<- function()
{
	DATA		<<- "/work/or105/Gates_2014"
	#DATA		<<- '/Users/Oliver/duke/2014_Gates'		
	if(0)	#all tasks combined
	{
		indir		<- paste(DATA,'methods_comparison_rootseqsim/140830',sep='/')
		#search for XML files in indir
		infiles		<- list.files(indir, pattern=paste(".xml$",sep=''))
		insignat	<- ''	
		hpc.ncpu	<- 8		
		for(infile in infiles)
		{
			infile		<- substr(infile, 1, nchar(infile)-4) 		
			cmd			<- hivc.cmd.beast.runxml(indir, infile, insignat, prog.beast=PR.BEAST, prog.beast.opt=" -beagle -working", hpc.tmpdir.prefix="beast", hpc.ncpu=hpc.ncpu)
			tmp			<- paste(infile,'.timetrees',sep='')	
			cmd			<- paste(cmd, hivc.cmd.beast.read.nexus(indir, tmp, indir, tree.id=NA, method.node.stat='any.node'), sep='\n')
			cmd			<- paste(cmd, hivc.cmd.beast.run.treeannotator(indir, infile, insignat, prog.beastmcc=PR.BEASTMCC, beastmcc.burnin=500, beastmcc.heights="median"), sep='\n')
			cat(cmd)	
			cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=hpc.ncpu, hpc.walltime=791, hpc.mem="7700mb")		
			outdir		<- indir
			outfile		<- paste("b2m.",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='')					
			hivc.cmd.hpccaller(outdir, outfile, cmd)		
		}
	}
	if(0)	#only parse existing output
	{
		indir		<- paste(DATA,'methods_comparison_rootseqsim/140813/save',sep='/')
		#search for XML files in indir
		infiles		<- list.files(indir, pattern=paste(".timetrees$",sep=''))							
		for(infile in infiles)
		{
			cmd			<- hivc.cmd.beast.read.nexus(indir, infile, indir, tree.id=NA, method.node.stat='any.node')
			cat(cmd)
			cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=1, hpc.walltime=348, hpc.mem="3700mb")
			outdir		<- indir
			outfile		<- paste("b2p.",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='')					
			hivc.cmd.hpccaller(outdir, outfile, cmd)
		}				
	}
	if(1)
	{			
		indir		<- paste(DATA,'methods_comparison_rootseqsim/140907',sep='/')
		#search for XML files in indir
		infiles		<- list.files(indir, pattern=paste(".xml$",sep=''))
		insignat	<- ''	
		hpc.ncpu	<- 8
		
		for(infile in infiles)
		{
			infile		<- substr(infile, 1, nchar(infile)-4) 		
			cmd			<- hivc.cmd.beast.runxml(indir, infile, insignat, prog.beast=PR.BEAST, prog.beast.opt=" -beagle -working", hpc.tmpdir.prefix="beast", hpc.ncpu=hpc.ncpu)
			tmp			<- paste(infile,'.timetrees',sep='')	
			cmd			<- paste(cmd, hivc.cmd.beast.read.nexus(indir, tmp, indir, tree.id=NA, method.node.stat='any.node'), sep='\n')
			cmd			<- paste(cmd, hivc.cmd.beast.run.treeannotator(indir, infile, insignat, prog.beastmcc=PR.BEASTMCC, beastmcc.burnin=500, beastmcc.heights="median"), sep='\n')
			cat(cmd)	
			cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=hpc.ncpu, hpc.walltime=791, hpc.mem="3700mb")		
			outdir		<- indir
			outfile		<- paste("bpg.",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='')					
			hivc.cmd.hpccaller(outdir, outfile, cmd)		
		}
	}
}
######################################################################################
project.Gates.dual.ExaMLrun<- function()
{
	require(phytools)
	require(hivclust)
	tree.id.labelsep		<- '|'
	tree.id.label.idx.ctime	<- 4 
	DATA		<<- "/work/or105/Gates_2014"
	#DATA		<<- "/Users/Oliver/duke/2014_Gates"
	indir		<- paste(DATA,'dual/141128',sep='/')	  
	outdir		<- indir
	infile		<- 'Contigs_141126_BLASTstrict750_OR_linsi.R'
	file		<- paste(indir, '/', infile, sep='')
	cat(paste('\nLoading file', file))
	load(file)		#expect "seqa"	
	
	infile.seq.sig	<- "Fri_Nov_28_12:59:06_2013"
	infile.seq		<- substr(infile,1,nchar(infile)-2)
	file			<- paste( outdir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
	seq				<- seqa
	save(seq, file=file)	
	#
	#	run ExaML 
	#
	if(1)
	{		
		bs.from		<- 0
		bs.to		<- 200
		bs.n		<- 200
		outdir		<- indir
		cmd			<- hivc.cmd.examl.bootstrap(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=bs.from, bs.to=bs.to,bs.n=bs.n,outdir=outdir, resume=1, verbose=1)
		dummy		<- lapply(cmd, function(x)
				{				
					#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
					x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=21, hpc.q="pqeph", hpc.mem="450mb", hpc.nproc=1)
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("du",signat,sep='.')
					#cat(x)
					hivc.cmd.hpccaller(outdir, outfile, x)
					Sys.sleep(1)
				})	
	}
	if(0)
	{
		#	
		cmd				<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=500, verbose=1)
		cmd				<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=71, hpc.q='pqeph', hpc.mem="3600mb", hpc.nproc=8)
		hivc.cmd.hpccaller(outdir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		Sys.sleep(1)	
	}		
}
######################################################################################
project.Gates.test.ExaMLrun<- function()
{
	require(phytools)
	require(hivclust)
	tree.id.labelsep		<- '|'
	tree.id.label.idx.ctime	<- 4 
	DATA		<<- "/work/or105/Gates_2014"
	#DATA		<<- "/Users/Oliver/duke/2014_Gates"
	indir		<- paste(DATA,'methods_comparison_pipeline/150130',sep='/')	  
	outdir		<- indir
	infiles		<- list.files(indir, '.*INTERNAL.R$', full.names=FALSE)
	#
	#	run ExaML 
	#
	for(i in seq_along(infiles))
	{
		infile		<- infiles[i]
		#	load simulated data
		file			<- paste(indir, '/', infile, sep='')
		cat(paste('\nLoading file', file))
		load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
		#	load outgroup sequences
		file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
		cat(paste('\nLoading outgroup seq from file', file))
		load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
		
		tmp				<- tolower(do.call('rbind',strsplit(df.seq[, GAG],'')))
		rownames(tmp)	<- df.seq[, LABEL]
		df.seq.gag		<- as.DNAbin(tmp)
		tmp				<- tolower(do.call('rbind',strsplit(df.seq[, POL],'')))
		rownames(tmp)	<- df.seq[, LABEL]
		df.seq.pol		<- as.DNAbin(tmp)	
		tmp				<- tolower(do.call('rbind',strsplit(df.seq[, ENV],'')))
		rownames(tmp)	<- df.seq[, LABEL]
		df.seq.env		<- as.DNAbin(tmp)
		if(0)
		{
			#
			#	run ExaML on gag
			#
			seq				<- df.seq.gag
			seq				<- rbind(seq, outgroup.seq.gag[, seq_len(ncol(seq))])
			infile.seq.sig	<- "Sun_Sep_14_12:59:06_2013"
			infile.seq		<- paste(substr(infile,1,nchar(infile)-20),'INFO_simu_gagseq',sep='')
			file			<- paste( outdir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
			save(seq, file=file)
			#	run ExaML
			cmd				<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
			cmd				<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= 'pqeelab', hpc.mem="450mb", hpc.nproc=1)
			hivc.cmd.hpccaller(outdir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
			Sys.sleep(1)						
		}
		if(1)
		{
			#
			#	run ExaML on pol
			#
			seq				<- df.seq.pol
			seq				<- rbind(seq, outgroup.seq.pol[, seq_len(ncol(seq))])
			infile.seq.sig	<- "Sun_Sep_14_12:59:06_2013"
			infile.seq		<- paste(substr(infile,1,nchar(infile)-20),'INFO_simu_polseq',sep='')
			file			<- paste( outdir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
			save(seq, file=file)
			#	run ExaML
			cmd				<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
			cmd				<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= 'pqeelab', hpc.mem="450mb", hpc.nproc=1)
			hivc.cmd.hpccaller(outdir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
			Sys.sleep(1)				
		}
		if(0)
		{
			#
			#	run ExaML on env
			#
			seq				<- df.seq.env
			seq				<- rbind(seq, outgroup.seq.env[, seq_len(ncol(seq))])
			infile.seq.sig	<- "Sun_Sep_14_12:59:06_2013"
			infile.seq		<- paste(substr(infile,1,nchar(infile)-20),'INFO_simu_envseq',sep='')
			file			<- paste( outdir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
			save(seq, file=file)
			#	run ExaML
			cmd				<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
			cmd				<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= 'pqeelab', hpc.mem="450mb", hpc.nproc=1)
			hivc.cmd.hpccaller(outdir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
			Sys.sleep(1)						
		}
		#
		#	run ExaML on concatenated
		#
		if(0)
		{
			seq				<- cbind(df.seq.gag,df.seq.pol,df.seq.env)
			tmp				<- cbind(outgroup.seq.gag[,1:ncol(df.seq.gag)], outgroup.seq.pol, outgroup.seq.env)
			seq				<- rbind(seq,tmp)
			infile.seq.sig	<- "Sun_Sep_14_12:59:06_2013"
			infile.seq		<- paste(substr(infile,1,nchar(infile)-20),'INFO_simu_concseq',sep='')
			file			<- paste( outdir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
			save(seq, file=file)
			#	run ExaML
			cmd				<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
			cmd				<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
			hivc.cmd.hpccaller(outdir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
			Sys.sleep(1)	
		}			
	}		
}
######################################################################################
project.Gates<- function()
{
	#project.Gates.RootSeqSim.runxml()
	#project.Gates.RootSeqSim.BEAST.SSAfg.checkancestralseq.runExaML()
	project.Gates.test.ExaMLrun()
	#project.Gates.test.runxml()
	#project.Gates.dual.ExaMLrun()
}
######################################################################################
project.hivc.check<- function()
{
	if(1) project.hivc.check.DateRes.after.HIVPosTest()
	if(0) project.hivc.check.DateRes.after.T0()	
}
######################################################################################
project.hivc.check.samplingmodel<- function()
{	
	ggplot(YXf, aes(x=score.Y, fill=stage)) + geom_histogram(binwidth=0.05) + facet_grid(stage ~ ., scales='free', margins=FALSE)
	
	z	<- gamlss(score.Y ~ stage-1, sigma.formula= ~ stage-1, family=BE, data=as.data.frame(subset(YX.m3, select=c(stage, score.Y))))
	#extract coefficients
	z	<- data.table(  risk='stage', factor=names(z[['mu.coefficients']]), mu= 1/(1+exp(-z[['mu.coefficients']])), sigma=1/(1+exp(-z[['sigma.coefficients']]))  )
	set(z, NULL, 'factor', z[, substr(factor, 6, nchar(factor))])
	z2	<- z[,  {
				score.Y<- seq(0.001,0.999,0.001)
				list(score.Y=score.Y, dBE=dBE(score.Y, mu=mu, sigma=sigma))
			}, by=c('risk','factor')]
	ggplot(z2, aes(x=score.Y, y=dBE, group=factor, colour=factor)) + geom_line() + facet_grid(factor ~ ., scales='free', margins=FALSE)
}
######################################################################################
project.hivc.check.Coal.threshold<- function()
{
	#Y.coal for linked --> usually high Y.coal but can go down to 0
	tmp		<- merge( subset(cluphy.info, select=c(FASTASampleCode, cluster)), subset(Y.rawbrl.linked, select=c(Patient, FASTASampleCode, t.FASTASampleCode)), by='FASTASampleCode' )		
	tmp		<- merge( data.table(t.FASTASampleCode= cluphy.info[, unique(FASTASampleCode)]), tmp, by='t.FASTASampleCode' )
	Y.U		<- project.athena.Fisheretal.Y.infectiontime(tmp, df.all, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, score.set.value=NA, method='for.infected')
	setnames(tmp, 'Patient', 't.Patient')
	setnames(Y.U, c('Patient','score.Inf'), c('t.Patient','U.score'))
	Y.coal.l	<- project.athena.Fisheretal.Y.coal(tmp, df.all, Y.U, cluphy, cluphy.info, cluphy.map.nodectime, coal.t.Uscore.min=0.01, coal.within.inf.grace= 0.25, t.period=t.period, save.file=NA, resume=0 )
	#Y.coal for unlinked --> 95% quantile is large, 
	tmp		<- subset(Y.rawbrl.unlinked, select=c(t.Patient, sc.FASTASampleCode, t.FASTASampleCode))
	setnames(tmp, 'sc.FASTASampleCode', 'FASTASampleCode')
	tmp		<- merge( subset(cluphy.info, select=c(FASTASampleCode, cluster)), tmp, by='FASTASampleCode' )		
	tmp		<- merge( data.table(t.FASTASampleCode= cluphy.info[, unique(FASTASampleCode)]), tmp, by='t.FASTASampleCode' )
	Y.U		<- project.athena.Fisheretal.Y.infectiontime(tmp, df.all, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, score.set.value=NA, method='for.transmitter')
	Y.coal.u	<- project.athena.Fisheretal.Y.coal(tmp, df.all, Y.U, cluphy, cluphy.info, cluphy.map.nodectime, coal.t.Uscore.min=0.01, coal.within.inf.grace= 0.25, t.period=t.period, save.file=NA, resume=0 )
	#Y.coal for female - female pairs
	tmp		<- merge( subset(cluphy.info, select=c(cluster, FASTASampleCode)), subset(df.all, select=c(FASTASampleCode, Patient, Sex)), by='FASTASampleCode' )
	tmp[, table(Sex)]
	#females were dropped from clusters. ugh.	
}
######################################################################################
project.hivc.check.Pjx<- function()
{
	#
	#
	#
	#
	z	<- subset(tmp, !is.nan(Pjx) & factor=='ART.suA.Y.4')
	hist( z[, Pjx], breaks=seq(0,1,0.01) )
	hist( z[, Pjx.e0cp], breaks=seq(0,1,0.01) )
	
	z	<- subset(tmp, !is.nan(Pjx) & factor=='ART.suA.Y.4' & Pjx.e0cp>0.2)
	setkey(z, Pjx.e0cp)
	#     risk Patient      factor       Pjx    Pjx.e0  Pjx.e0cp             coef
	 #1: stage  M38112 ART.suA.Y.4 0.2456374 0.2262577 0.2052394 stageART.suA.Y.4
	 #2: stage  M39212 ART.suA.Y.4 0.2638938 0.2374935 0.2148885 stageART.suA.Y.4
	 #3: stage  M39252 ART.suA.Y.4 0.2638938 0.2376439 0.2159022 stageART.suA.Y.4
	 #4: stage  M38097 ART.suA.Y.4 0.3029802 0.2660489 0.2312217 stageART.suA.Y.4
	 #5: stage  M38876 ART.suA.Y.4 0.3099701 0.2782368 0.2490324 stageART.suA.Y.4
	 #6: stage  M38750 ART.suA.Y.4 0.3109900 0.2815144 0.2539564 stageART.suA.Y.4
	 #7: stage  M39137 ART.suA.Y.4 0.3431932 0.3103555 0.2831176 stageART.suA.Y.4
	 #8: stage  M41876 ART.suA.Y.4 0.3976996 0.3476053 0.3031130 stageART.suA.Y.4
	 #9: stage  M39614 ART.suA.Y.4 0.7455575 0.4325278 0.3132440 stageART.suA.Y.4
	#10: stage  M37781 ART.suA.Y.4 0.4381747 0.3770513 0.3227846 stageART.suA.Y.4
	#11: stage  M38881 ART.suA.Y.4 0.5905905 0.4804145 0.4031972 stageART.suA.Y.4
	#12: stage  M38623 ART.suA.Y.4 0.6825126 0.5088493 0.4047791 stageART.suA.Y.4
	#13: stage  M38138 ART.suA.Y.4 0.6239798 0.5043751 0.4132018 stageART.suA.Y.4
	#14: stage  M37935 ART.suA.Y.4 1.0000000 0.7251842 0.5512027 stageART.suA.Y.4
	
	z2	<- merge( YX, unique(subset(z, select=Patient)), by='Patient')
	z2[,list(np=length(unique(t.Patient)), nt=length(unique(t))), by='Patient']
	
	z2	<- subset(z2, CD4b.tperiod=='ART.suA.Y.4', c(Patient, t.Patient, t, score.Y, FASTASampleCode, t.FASTASampleCode))
	setkey(z2, Patient, t.Patient)
	unique(z2)
	
	subset( missing, factor=='ART.suA.Y.4'  & yYX.sum>1)
	z	<- subset(merge( YX, unique(subset( missing, factor=='ART.suA.Y.4'  & yYX.sum>1, Patient)), by='Patient'), CD4b.tperiod=='ART.suA.Y.4', c(Patient, t.Patient, t, score.Y, FASTASampleCode, t.FASTASampleCode))
	
	#	try set COAL to include these
	subset(df.tpairs, Patient=='M37935' & t.Patient=='M38982')
	subset(df.tpairs, Patient=='M38138' & t.Patient=='M34329')
	subset(df.tpairs, Patient=='M37781' & t.Patient=='M38876')
	subset(df.tpairs, Patient=='M37781' & t.Patient=='M38506')
	subset(df.tpairs, Patient=='M39614' & t.Patient=='M35997')
	subset(df.tpairs, Patient=='M41876' & t.Patient=='M36082')
	
	
}
######################################################################################
project.hivc.get.geneticdist.from.sdc<- function(dir.name= DATA)
{	
	tmp<- hivc.clu.geneticdist.cutoff(dir.name=dir.name, plot=1, verbose=1, level.retain.unlinked=0.05)
	print(tmp)
}
######################################################################################
project.seq.dataset.mDR.mRC.mSH.pLANL<- function()	
{
	require(data.table)
	require(phangorn)
	
	file		<- paste(DATA,"/tmp/","ATHENA_2013_03_NoRCDRAll+LANL_Sequences","_",gsub('/',':',"Fri_Nov_01_16/07/23_2013"),".R",sep='')
	if(verbose) cat(paste("\nread",file))
	tmp			<- load(file)
	seq.len		<- seq.length(seq.PROT.RT)		
	hist(seq.len, breaks=20)
	seq.amb		<- seq.proportion.ambiguous(seq.PROT.RT)
	hist(seq.amb, breaks=40)
	
	
	#
	# get clusters for No Recombination + No Drug resistance mutations, single linkage criterion		
	#						
	infilecov		<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infile			<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences_examlbs500"			
	insignat		<- "Fri_Nov_01_16/07/23_2013"		
	argv			<<- hivc.cmd.preclustering(paste(DATA,"/tmp",sep=''), infile, insignat, paste(DATA,"/derived",sep=''), infilecov, resume=1)				 
	argv			<<- unlist(strsplit(argv,' '))
	nrc.clu.pre		<- hivc.prog.get.clustering.precompute()
	#
	ph				<- nrc.clu.pre$ph
	ph.node.bs		<- as.numeric(ph$node.label)
	
	seq.stat		<- data.table(FASTASampleCode=names(seq.len), tip:=match(seq.stat[,FASTASampleCode],ph$tip.label), length=seq.len, pamb=seq.amb)
	seq.stat		<- seq.stat[,	{
				tmp		<- Ancestors(ph, tip, type='all') - Ntip(ph)
				list(bs.mx=max(ph.node.bs[tmp]), length=length, pamb=pamb, tip=tip)			
			},by=FASTASampleCode]		
	seq.short		<- subset(seq.stat, length<400)
	seq.long		<- subset(seq.stat, length>=400)
	seq.short.nfrgn	<- subset(seq.short,substr(seq.short[,FASTASampleCode],1,2)!="TN")		
	seq.closefrgn	<- subset(seq.stat, substr(seq.stat[,FASTASampleCode],1,8)=="PROT+P51") 
	seq.tn			<- subset(seq.stat, substr(seq.stat[,FASTASampleCode],1,2)=="TN")	
	
	hist(seq.short[,bs.mx], breaks=20)			
	hist(seq.long[,bs.mx], breaks=20)
	hist(seq.short.nfrgn[,bs.mx], breaks=20)
	hist(seq.closefrgn[,bs.mx], breaks=20)
	#	-> longer sequences have larger bootstrap
	#	-> short ATHENA sequences have small bootstap
	#	-> almost all short sequences (500 / 548 ) are TN.	so TNs might not cluster simply because they are short. restrict TNs to foreign sequences > 600
	#	-> all PROT+P51 sequences are long, and do not cluster as well as the NL sequences
	seq.athena.exclude	<- subset(seq.short.nfrgn, bs.mx<0.6 )[, FASTASampleCode]		
	#			
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences"	
	insignat	<- "Fri_Nov_01_16/07/23_2013"
	outfile		<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"	
	outsignat	<- "Wed_Dec_18_11/37/00_2013"		
	
	file		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
	load(file)
	#	exclude short ATHENA sequences that have BS<0.6 in pilot run AND all TN sequences because they are <400 nt
	seq.keep	<- setdiff(rownames(seq.PROT.RT), seq.athena.exclude)
	seq.keep	<- seq.keep[ substr(seq.keep,1,2)!="TN" ]				
	seq.PROT.RT	<- seq.PROT.RT[seq.keep,]
	if(verbose)	cat(paste("\nnumber of long sequences, n=",nrow(seq.PROT.RT)))
	file		<- paste(indir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
	if(verbose)	cat(paste("\nsave new file to",file))
	save(seq.PROT.RT, file=file)
}

######################################################################################
project.WTprop.ATHENAmobility.160309<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	indir		<- "~/Dropbox (Infectious Disease)/2015_ATHENA_May_Update"
	infile.main	<- file.path(indir,"ATHENA_1502_All_PatientKeyCovariates.R")
	infile.ggd	<- file.path(indir,"ATHENA_1502_All_Region_GGD_v02.R")
	outdir		<- "/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2016/2016_GuidingTransmissionElimination"
	
	load(infile.main)
	df.all		<- subset(df.all, select=c(Patient, DateBorn, Sex, CountryBorn, CountryInfection, Region_RegT1, Region_now, RegionOrigin, Trm, DateLastContact, DateDied, isDead, AnyPos_T1, PosSeqT))
	set(df.all, NULL, "DateBorn", df.all[,hivc.db.Date2numeric( DateBorn )])
	set(df.all, NULL, "DateLastContact", df.all[,hivc.db.Date2numeric( DateLastContact )])	
	set(df.all, NULL, "DateDied", df.all[,hivc.db.Date2numeric( DateDied )])
	set(df.all, NULL, "AnyPos_T1", df.all[,hivc.db.Date2numeric( AnyPos_T1 )])
	set(df.all, df.all[, which(is.na(DateDied))], 'DateDied', 2015+5/12)
	#time since diagnosis
	set(df.all, NULL, 'TsD', df.all[, DateDied-AnyPos_T1])	
	#	age at diagnosis
	set(df.all, NULL, "Age_AnyPosT1", df.all[, AnyPos_T1-DateBorn])		
	#	redefine RegionOrigin
	set(df.all, df.all[, which(RegionOrigin%in%c("Central_EU"))], "RegionOrigin", "Eastern_EU_stans")
	set(df.all, df.all[, which(RegionOrigin%in%c("Austr_NZ","Sout_SouthEast_Asia","Oceania_Pacific","North_Africa_Middle_East"))], "RegionOrigin", "Other")	
	set(df.all, NULL, "RegionOrigin", df.all[, factor(RegionOrigin)])
	#	redefine Tansmission
	set(df.all, df.all[, which(Trm%in%c("MSM",'BI'))], "Trm", "MSM")
	set(df.all, df.all[, which(Trm%in%c("HET",'HETfa'))], "Trm", "HET")
	set(df.all, df.all[, which(Trm%in%c("IDU","BLOOD","NEEACC","PREG","BREAST","SXCH"))], "Trm", "OTH")
	set(df.all, df.all[,which(is.na(Trm))], "Trm", "Unknown")		
	set(df.all, df.all[, which(Trm=='HET' & Sex=='M')], 'Trm', 'HETM')
	set(df.all, df.all[, which(Trm=='HET' & Sex=='F')], 'Trm', 'HETF')
	set(df.all, NULL, "Trm", df.all[, factor(Trm)])	
	#	reduce to patients with one sequence
	df.all[, DUMMY:=seq_len(nrow(df.all))]
	tmp			<- df.all[, list(DUMMY=ifelse(all(is.na(PosSeqT)), DUMMY, DUMMY[which.min(PosSeqT)] )), by='Patient']
	df.all		<- merge(df.all, tmp, by=c('Patient','DUMMY'))
	df.all[, DUMMY:=NULL]
	df.all[, PosSeqT:=NULL]
	#	define Amsterdam
	df.all[, AMST:= NA_character_]
	set(df.all, df.all[, which( Region_RegT1=='Amst')], 'AMST', 'Y')
	set(df.all, df.all[, which( Region_RegT1!='Amst')], 'AMST', 'N')
	set(df.all, df.all[, which( is.na(Region_RegT1) & Region_now=='Amst')], 'AMST', 'Y')
	set(df.all, df.all[, which( is.na(Region_RegT1) & Region_now!='Amst')], 'AMST', 'N')
	set(df.all, df.all[, which( !is.na(CountryInfection) & CountryInfection!='NL')], 'AMST', 'N')	
	stopifnot( !nrow(subset(df.all, is.na(AMST) & AnyPos_T1>2010)) )
	#	define time period
	set(df.all, NULL, 'TP', df.all[, cut(AnyPos_T1, breaks=c(-Inf, 2010, 2011, 2012, 2013, 2014, 2015, 2016), labels=c('<2010','2010','2011','2012','2013','2014','2015'))])
	#	define young
	set(df.all, NULL, 'YOUNG', df.all[, cut(Age_AnyPosT1, breaks=c(0,28,100), labels=c('16-27','28-80'))])
	#	define migrant
	set(df.all, NULL, 'MIGRANT', df.all[, RegionOrigin])
	#
	#
	#
	load(infile.ggd)
	set(df, NULL, "GGD_Reg", df[,hivc.db.Date2numeric( GGD_Reg )])	
	tmp			<- subset(df, !is.na(GGD))[, list(GGD_N=length(GGD_Reg)), by='Patient']
	df			<- merge(df.all, tmp, by='Patient', all.x=1)
	#
	#	mobile with GGD_N>1 of those diagnosed since 2010 (numbers small, may not be all transmitters)
	#
	tmp			<- subset(df, !is.na(MIGRANT) & Trm%in%c('HETF','HETM','MSM') & AnyPos_T1>=2010)	
	tmp[, table(GGD_N, useNA='if')]
	tmp[, list(Patient_N=length(Patient), GGD_CH_ATLEAST_ONCE=mean(GGD_N>1, na.rm=TRUE) ), by='MIGRANT']
	tmp2		<- tmp[, list(Patient_N=length(Patient), GGD_CH_ATLEAST_ONCE=mean(GGD_N>1, na.rm=TRUE) ), by=c('MIGRANT','Trm')]
	setkey(tmp2, MIGRANT, Trm)
	tmp2
	#
	#	mobile with GGD_N>1 of those diagnosed since 2001 (perhaps including all transmitters..)
	#
	tmp			<- subset(df, !is.na(MIGRANT) & Trm%in%c('HETF','HETM','MSM') & AnyPos_T1>=2001 & AMST=='Y')	
	tmp[, table(GGD_N, useNA='if')]
	tmp[, list(Patient_N=length(Patient), GGD_CH_ATLEAST_ONCE=mean(GGD_N>1, na.rm=TRUE) ), by='MIGRANT']
	tmp2		<- tmp[, list(Patient_N=length(Patient), GGD_CH_ATLEAST_ONCE=mean(GGD_N>1, na.rm=TRUE) ), by=c('MIGRANT','Trm')]
	setkey(tmp2, MIGRANT, Trm)
	ggplot(tmp2, aes(x=Trm, fill=MIGRANT, y=GGD_CH_ATLEAST_ONCE)) +
		theme_bw() + theme(legend.position='bottom') +
		geom_bar(stat='identity', position=position_dodge(width=0.9)) +
		geom_text(aes(label=round(Patient_N*GGD_CH_ATLEAST_ONCE)), position=position_dodge(width=0.9), vjust=-0.25) +
		scale_y_continuous(labels=percent) +
		scale_fill_brewer(palette='Dark2') +				
		labs(x='\nDiagnosed in Amsterdam after 2000', y='At least 1 GGD change\n', fill='Region of Origin')
	ggsave(file=file.path(outdir, 'ATHENA0502_NewDiag0114AMST_GGDCH_by_regionorigin.pdf'), w=8, h=4)
	#
	#	mobile with GGD_N>1 of those diagnosed since 2001 (perhaps including all transmitters..)
	#
	tmp			<- subset(df, !is.na(MIGRANT) & !is.na(AMST) & Trm%in%c('HETF','HETM','MSM') & AnyPos_T1>=2001)
	set(tmp, NULL, 'MIGRANT', tmp[,factor(as.character(MIGRANT)=='NL', levels=c(TRUE,FALSE), labels=c('Dutch origin','Migrant'))])
	set(tmp, NULL, 'AMST', tmp[,factor(as.character(AMST), levels=c('Y','N'), labels=c('First registered in Amsterdam','First registered outside Amsterdam'))])
	tmp			<- tmp[, list(Patient_N=length(Patient), GGD_CH_ATLEAST_ONCE=mean(GGD_N>1, na.rm=TRUE) ), by=c('MIGRANT','Trm','AMST')]
	setkey(tmp, MIGRANT, Trm, AMST)
	ggplot(tmp, aes(x=Trm, fill=MIGRANT, y=GGD_CH_ATLEAST_ONCE)) +
			theme_bw() + theme(legend.position='bottom') +
			geom_bar(stat='identity', position=position_dodge(width=0.9)) +
			geom_text(aes(label=round(Patient_N*GGD_CH_ATLEAST_ONCE)), position=position_dodge(width=0.9), vjust=-0.25) +
			scale_y_continuous(labels=percent, limits=c(0,0.3)) +
			scale_fill_manual(values=c('Dutch origin'="#41B6C4", 'Migrant'="#FE9929")) +	
			facet_grid(~AMST) +
			labs(x='\nDiagnosed in the Netherlands since 2001', y='At least 1 GGD change\n', fill='Region of Origin')
	ggsave(file=file.path(outdir, 'ATHENA0502_NewDiag0114AMST_GGDCH_by_Amsterdam.pdf'), w=8, h=4)
	#
	#	mobile with GGD_N>1 of those diagnosed since 2010 (perhaps including all transmitters..)
	#
	tmp			<- subset(df, !is.na(MIGRANT) & !is.na(AMST) & Trm%in%c('HETF','HETM','MSM') & AnyPos_T1>=2010)
	set(tmp, NULL, 'MIGRANT', tmp[,factor(as.character(MIGRANT)=='NL', levels=c(TRUE,FALSE), labels=c('Dutch origin','Migrant'))])
	set(tmp, NULL, 'AMST', tmp[,factor(as.character(AMST), levels=c('Y','N'), labels=c('First registered in Amsterdam','First registered outside Amsterdam'))])
	tmp			<- tmp[, list(Patient_N=length(Patient), GGD_CH_ATLEAST_ONCE=mean(GGD_N>1, na.rm=TRUE) ), by=c('MIGRANT','Trm','AMST')]
	setkey(tmp, MIGRANT, Trm, AMST)
	ggplot(tmp, aes(x=Trm, fill=MIGRANT, y=GGD_CH_ATLEAST_ONCE)) +
			theme_bw() + theme(legend.position='bottom') +
			geom_bar(stat='identity', position=position_dodge(width=0.9)) +
			geom_text(aes(label=round(Patient_N*GGD_CH_ATLEAST_ONCE)), position=position_dodge(width=0.9), vjust=-0.25) +
			scale_y_continuous(labels=percent, limits=c(0,0.3)) +
			scale_fill_manual(values=c('Dutch origin'="#41B6C4", 'Migrant'="#FE9929")) +	
			facet_grid(~AMST) +
			labs(x='\nDiagnosed in the Netherlands since 2010', y='At least 1 GGD change\n', fill='Region of Origin')
	ggsave(file=file.path(outdir, 'ATHENA0502_NewDiag114AMST_GGDCH_by_Amsterdam.pdf'), w=8, h=4)
	
}
######################################################################################





######################################################################################

######################################################################################
project.hivc.clustering.get.linked.and.unlinked<- function(dir.name= DATA)
{
	require(data.table)
	require(ape)
	if(1)	#extract unlinked and linked pairs -- this is version Sat_Jun_16_17/23/46_2013
	{
		verbose				<- 1
		
		#load all+enriched sequences
		indir				<- paste(dir.name,"tmp",sep='/')
		infile				<- "ATHENA_2013_03_CurAll+LANL_Sequences"
		insignat			<- "Sat_Jun_16_17/23/46_2013"		
		file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
		if(verbose)	cat(paste("\nread file",file))
		load(file)
		print(seq.PROT.RT)
		#exctract geographically distant seqs that are assumed to be truly unlinked to NL seqs
		seq.PROT.RT.TN		<- seq.PROT.RT[substr(rownames(seq.PROT.RT),1,2)=="TN", ]
		#extract ATHENA seqs
		seq.PROT.RT.NL		<- seq.PROT.RT[substr(rownames(seq.PROT.RT),1,2)!="TN", ]
		seq.PROT.RT.NL		<- seq.PROT.RT.NL[substr(rownames(seq.PROT.RT.NL),1,8)!="PROT+P51", ]
		
		#need patient id and PosSeq for each sample code		
		indir				<- paste(dir.name,"derived",sep='/')
		infile.seqinfo		<- "ATHENA_2013_03_Sequences"
		file				<- paste(indir,'/',paste(infile.seqinfo,".R", sep=''), sep='')
		if(verbose) cat(paste("\nloading file",file))
		load(file)	
		df.PROT				<- as.data.table(df[["PROT"]][,c("Patient","SampleCode","DateRes")])
		setkey(df.PROT,"SampleCode")
		df.RT				<- as.data.table(df[["RT"]][,c("Patient","SampleCode","DateRes")])
		setkey(df.RT,"SampleCode")
		df					<- merge(df.PROT,df.RT,all=TRUE)
		if( !all( df[,Patient.x==Patient.y], na.rm=1) ) stop("Patient names per sampling code in PROT and RT don t equal")
		if( !all( df[,DateRes.x==DateRes.y], na.rm=1) ) stop("Sampling dates per sampling code in PROT and RT don t equal")
		tmp					<- df[,DateRes.x]
		tmp2				<- df[,is.na(DateRes.x)]
		tmp[tmp2]			<- df[tmp2,DateRes.y]	
		df[,PosSeqT:= tmp]
		tmp					<- df[,Patient.x]
		tmp2				<- df[,is.na(Patient.x)]
		tmp[tmp2]			<- df[tmp2,Patient.y]	
		df[,Patient:= tmp]
		if(verbose) cat(paste("\nfound ATHENA seqs info, n=",nrow(df)))
		df.seqinfo			<- subset(df, !is.na(PosSeqT), select=c(SampleCode, Patient, PosSeqT))
		tmp					<- df.seqinfo[,gsub(' ','',SampleCode,fixed=1)]
		df.seqinfo[,FASTASampleCode:=tmp]
		
		if(verbose) cat(paste("\nfound ATHENA seqs with known sampling date, n=",nrow(df.seqinfo)))
		if(verbose) cat(paste("\nfound ATHENA patients whose sequences have known sampling date, n=",length(unique(df.seqinfo[,Patient]))))
		#extract list of truly linked sample codes
		setkey(df.seqinfo, "Patient")
		tmp					<- subset(df.seqinfo[, length(SampleCode)>1, by=Patient], V1==T, select=Patient)
		linked.bypatient	<- merge(tmp, df.seqinfo, all.x=1)
		setkey(linked.bypatient, "FASTASampleCode")
		linked.bypatient	<- subset(linked.bypatient, select=c(Patient, FASTASampleCode, PosSeqT))
		
		#extract list of truly unlinked seqs -- temporal separation
		indir				<- paste(dir.name,"derived",sep='/')
		infile.patient		<- "ATHENA_2013_03_All1stPatientCovariates"
		file				<- paste(indir,'/',paste(infile.patient,".R", sep=''), sep='')
		if(verbose) cat(paste("\nloading file",file))
		load(file)	
		#extract seqs of dead invidividuals
		df.dead						<- subset(df.all, !is.na(Died), c(Patient,Died))
		df.dead						<- merge(df.dead, df.seqinfo, by="Patient")
		setkey(df.dead,Died)
		#extract seroconverters
		df.serocon.acc				<- subset(df.all, NegT_Acc=="Yes" & NegT>=df.dead[1,Died],)
		#add seroconverters with inaccurate info
		df.serocon.nacc				<- subset(df.all, NegT_Acc=="No" & !is.na(NegT) & NegT>=df.dead[1,Died], )
		df.serocon.nacc.dy			<- subset(df.serocon.nacc, as.POSIXlt(NegT)$mday==15, )						#for inaccurate days, we (conservatively) assume the patient was only seronegative at the start of the month
		tmp							<- as.POSIXlt(df.serocon.nacc.dy[,NegT] )
		tmp$mday					<- 1
		df.serocon.nacc.dy[,NegT:=as.Date(tmp)]
		df.serocon.nacc.mody		<- subset(df.serocon.nacc, as.POSIXlt(NegT)$mon==6 & as.POSIXlt(NegT)$mday==1, )		#for inaccurate months and days, we (conservatively) assume the patient was only seronegative at the start of the year
		tmp							<- as.POSIXlt(df.serocon.nacc.mody[,NegT] )
		tmp$mon						<- 0
		df.serocon.nacc.mody[,NegT:=as.Date(tmp)]
		df.serocon					<- rbind(df.serocon.acc, df.serocon.nacc.dy, df.serocon.nacc.mody)		
		if(verbose) cat(paste("\nnumber of seroconverters with at least 1 preceeding dead HIV+",nrow(df.serocon)))
		#for each accurate seroconverter, extract HIV+ seqs that are dead before seroconversion
		if( any(as.logical(df.serocon[,is.na(NegT)])) )	warning("Found accurate seroconverters with missing NegT")
		df.serocon					<- merge(df.serocon,df.seqinfo,all.x=1,by="Patient")		
		setkey(df.serocon,NegT)
		unlinked.bytime				<- lapply(seq_len(nrow(df.serocon)), function(i)
				{
					tmp				<- subset(df.dead, Died<=df.serocon[i,NegT],select=c(Patient,FASTASampleCode,Died))
					tmp2			<- rep(df.serocon[i,FASTASampleCode],nrow(tmp))
					tmp[,query.FASTASampleCode:= tmp2]
					setkey(tmp,"FASTASampleCode")
					tmp					
				})
		names(unlinked.bytime)	<- df.serocon[,FASTASampleCode]			
		
		#need also NegT for each seq			
		df.seqinfo				<- merge(df.seqinfo, df.serocon[,list(NegT=min(NegT)),by=Patient], all.x=1, by="Patient")
		
		#extract list of truly unlinked seqs -- geographical separation
		unlinked.byspace		<- data.table(FASTASampleCode=rownames(seq.PROT.RT.TN), key="FASTASampleCode")
						
		outdir					<- paste(dir.name,"tmp",sep='/')
		outfile					<- "ATHENA_2013_03_Unlinked_and_Linked"
		outsignat				<- "Sat_Jun_16_17/23/46_2013"
		file					<- paste(outdir,'/',outfile,'_',gsub('/',':',insignat),".R",sep='')
		if(verbose) cat(paste("\nwrite linked and unlinked to file",file))
		save(unlinked.byspace,unlinked.bytime,linked.bypatient,df.seqinfo,file=file)
		
		#
		#some quick statistics
		#
		if(0)
		{
			#number of unlinked HIV+ by date (this date is when someone else is still seroneg)
			y		<- sapply(unlinked.bytime, function(x) nrow(x) )
			y2		<- nrow(unlinked.byspace)
			x		<- df.serocon[,NegT]
			xlim	<- range(x)
			ylim	<- c(0, y2+max(y))
			tmp		<- as.POSIXlt(xlim[1])
			tmp$mday<- 1
			tmp$mon	<- 1
			xlim[1]	<- as.Date(tmp)
			plot(1,1,type='n',bty='n',xlab="time of seroconversion",ylab="number unlinked",xaxt='n',xlim=xlim,ylim=ylim)
			axis.Date(1,Year,at=seq(xlim[1], xlim[2],by="12 months"))			
			polygon(c(x[1],x[length(x)],x[length(x)],x[1]), c(y2,y2,0,0), border=NA, col="blue" )
			polygon(c(x,c(x[length(x)],x[1])), c(y+y2,y2,y2), border=NA, col="grey60" )
			legend("topleft",bty='n',fill=c("blue","grey60"),legend=c("by space","by time"),border=NA)						
		}
		quit("no")
	}
	if(0)	#extract unlinked pairs by temporal separation -- this is version Fri_May_24_12/59/06_2013
	{
		verbose				<- 1
		unlinked.closest.n	<- NA
		indir				<- paste(dir.name,"derived",sep='/')
		outdir				<- paste(dir.name,"derived",sep='/')
		infile				<- "ATHENA_2013_03_All1stPatientCovariates.R"
		outfile				<- "ATHENA_2013_03_Unlinked_SeroConv_Dead"
		outsignat			<- "Fri_May_24_12/59/06_2013"
		file				<- paste(indir,infile,sep='/')
		if(verbose)
		{
			cat(paste("\nunlinked.closest.n",unlinked.closest.n))			
		}
		
		if(verbose)	cat(paste("\nread file",file))
		load(file)
		if(verbose) str(df.all)
		
		df.dead						<- subset(df.all, !is.na(Died), c(Patient,Died))
		setkey(df.dead,Died)
		df.serocon.acc				<- subset(df.all, NegT_Acc=="Yes" & NegT>=df.dead[1,Died],)
		if(1)
		{
			df.serocon.nacc				<- subset(df.all, NegT_Acc=="No" & !is.na(NegT) & NegT>=df.dead[1,Died], )
			#for inaccurate days, we (conservatively) assume the patient was only seronegative at the start of the month
			df.serocon.nacc.dy			<- subset(df.serocon.nacc, as.POSIXlt(NegT)$mday==15, )
			tmp							<- as.POSIXlt(df.serocon.nacc.dy[,NegT] )
			tmp$mday					<- 1
			df.serocon.nacc.dy[,NegT:=as.Date(tmp)]
			#for inaccurate months and days, we (conservatively) assume the patient was only seronegative at the start of the year
			df.serocon.nacc.mody		<- subset(df.serocon.nacc, as.POSIXlt(NegT)$mon==6 & as.POSIXlt(NegT)$mday==1, )
			tmp							<- as.POSIXlt(df.serocon.nacc.mody[,NegT] )
			tmp$mon						<- 0
			df.serocon.nacc.mody[,NegT:=as.Date(tmp)]
			#merge all
			df.serocon					<- rbind(df.serocon.acc, df.serocon.nacc.dy, df.serocon.nacc.mody)
		}
		else
			df.serocon					<- df.serocon.acc
		
		if(verbose) cat(paste("\nnumber of seroconverters with at least 1 preceeding dead HIV+",nrow(df.serocon)))
		#for each accurate seroconverter, extract HIV+ that are dead before seroconversion
		if( any(as.logical(df.serocon[,is.na(NegT)])) )	warning("Found accurate seroconverters with missing NegT")		
		setkey(df.serocon,NegT)
		unlinked.bytime				<- lapply(seq_len(nrow(df.serocon)), function(i)
				{
					tmp<- subset(df.dead, Died<=df.serocon[i,NegT],)											
					#tmp<- tmp[,Patient]
					#if(length(tmp)<unlinked.closest.n)	
					#	tmp<- c(rep(NA,unlinked.closest.n-length(tmp)),tmp)
					#rev(tmp)[1:unlinked.closest.n]
					if(is.na(unlinked.closest.n))	return( tmp[,Patient] )
					else							return( rev(tmp)[1:min(length(tmp),unlinked.closest.n)] )
				})
		names(unlinked.bytime)	<- df.serocon[,Patient]			
		#
		#some quick statistics
		#
		if(1)
		{
			#number of unlinked HIV+ by date (this date is when someone else is still seroneg)
			y		<- sapply(unlinked.bytime, function(x) length(x) )
			x		<- df.serocon[,NegT]
			xlim	<- range(x)
			tmp		<- as.POSIXlt(xlim[1])
			tmp$mday<- 1
			tmp$mon	<- 1
			xlim[1]	<- as.Date(tmp)
			plot(1,1,type='n',bty='n',xlab="time of seroconversion",ylab="number unlinked",xaxt='n',xlim=xlim,ylim=range(y))
			axis.Date(1,Year,at=seq(xlim[1], xlim[2],by="12 months"))
			polygon(c(x,c(x[length(x)],x[1])), c(y,0,0), border=NA, col="grey60" )
			
			#average PosT for all unlinked patients by date (this date is when someone else is still seroneg)
			unlinked.PosT	<- sapply(seq_along(unlinked.bytime), function(i)
					{
						tmp<- data.table(Patient=unlinked.bytime[[i]], key="Patient")
						tmp<- subset(df.all[tmp], select=c(Patient, PosT))
						tmp<- tmp[,mean(as.numeric(difftime(df.serocon[i,NegT],PosT,units="weeks"))/52,na.rm=1)]
						tmp
					})
			y		<- unlinked.PosT				
			par(mar=c(5,6,1,1))
			xlim[1]	<- as.Date("2000-01-01")
			plot(1,1,type='n',bty='n',xlab="time of seroconversion",ylab="median time difference\nbetween unlinked sequences [yrs]",xaxt='n',xlim=xlim,ylim=range(y))
			axis.Date(1,Year,at=seq(xlim[1], xlim[2],by="12 months"))
			lines(x,y)
		}
		#
		#print(unlinked.bytime.n)
		#print(any(diff(unlinked.bytime.n)<0))
		#
		#save
		#
		if(is.na(unlinked.closest.n))			
			file						<- paste(outdir,paste(outfile,"_UnlinkedAll_",gsub('/',':',outsignat),".R",sep=''),sep='/')
		else
			file						<- paste(outdir,paste(outfile,"_Unlinked",unlinked.closest.n,'_',gsub('/',':',outsignat),".R",sep=''),sep='/')
		if(verbose) cat(paste("\nwrite unlinked pairs to file",file))
		save(unlinked.bytime, df.serocon, df.all, file=file)		
		quit("no")
	}
}
######################################################################################
project.hivc.tables.fixup.m2Bwmx<- function()
{
	tps	<- 1:4
	tmp	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3da_tablesCLU_m2Bwmx.tp'
	#tmp	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3ea_tablesCLU_m2Bwmx.tp'
	#tmp	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3fa_tablesCLU_m2Bwmx.tp'
	#tmp	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3ga_tablesCLU_m2Bwmx.tp'

	dummy	<- sapply(tps, function(tp)
			{
				file	<- paste(tmp, tp, '.R',sep='')
				load(file)
				set(ans$cens.table, NULL, 'factor', ans$cens.table[,gsub('SuA','ART.suA', factor)])
				set(ans$cens.table, NULL, 'factor2', ans$cens.table[,gsub('SuA','ART.suA', factor2)])
				set(ans$risk.table, NULL, 'factor', ans$risk.table[,gsub('SuA','ART.suA', factor)])
				set(ans$adj.clu, NULL, 'factor', ans$adj.clu[,gsub('SuA','ART.suA', factor)])
				set(ans$adj.seq, NULL, 'factor', ans$adj.seq[,gsub('SuA','ART.suA', factor)])
				cat(paste('\n save to file',file))
				save(ans, file=file)			
			})
	
}
project.hivc.clustering.NoRecombNoDR.to.NoShort<- function()
{	
	verbose		<- 1
	resume		<- 1
	#patient.n	<- 15700; 	
	thresh.brl	<- 1000
	opt.brl		<- "dist.brl.casc" 
	indircov	<- paste(DATA,"derived",sep='/')
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"
	indir		<- paste(DATA,"tmp",sep='/')
	#
	# get clusters for No Recombination + No Drug resistance mutations + No short sequences, single linkage criterion		
	#						
	infile			<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500"			
	insignat		<- "Wed_Dec_18_11/37/00_2013"
	#
	argv			<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	nsh.clu.pre		<- hivc.prog.get.clustering.precompute()
	#
	resume			<- 0
	#argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume)
	#argv			<<- unlist(strsplit(argv,' '))
	#nsh.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=nsh.clu.pre, with.plot=0)
	#			
	
	thresh.bs		<- 0.8	
	thresh.bss		<- c( 0.85, 0.9, 0.95 )
	dummy			<- lapply(thresh.bss, function(thresh.bs)
			{
				argv			<<- hivc.cmd.clustering(indir, infile, insignat, opt.brl="dist.brl.casc", thresh.brl, thresh.bs, resume=resume)				 
				argv			<<- unlist(strsplit(argv,' '))
				nsh.clu			<- hivc.prog.get.clustering()		
				#
				argv			<<- hivc.cmd.clustering.msm(indir, infile, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=resume)
				argv			<<- unlist(strsplit(argv,' '))		
				nsh.msm			<- hivc.prog.get.clustering.MSM()			
			})			
}
######################################################################################
project.hivc.clustering.forStephane.onUK.ExaML<- function()
{	
	verbose			<- 1
	resume			<- 1
	patient.n		<- 15700; 	thresh.brl<- 0.096; 	thresh.bs<- 0.8;	opt.brl<- "dist.brl.casc" 
	indir			<- '~/Dropbox (Infectious Disease)/2015_OptimalClusteringThresholds/data'	
	outdir			<- '~/Dropbox (Infectious Disease)/2015_OptimalClusteringThresholds/data'
	
	#	read tree
	if(0)
	{
		infile			<- "TNTP_PRRT_withoutDoublon_withoutDRM_examlbs500_151023.newick"
		ph				<- read.tree(paste(indir,infile,sep='/'))
		ph 				<- ladderize( ph )
		ph.node.bs						<- as.numeric( ph$node.label )		
		ph.node.bs[is.na(ph.node.bs)]	<- 0
		ph.node.bs[c(2,3,4,5,6)]		<- 0
		ph$node.label					<- ph.node.bs / 100	
		ph$tip.label	<- toupper(ph$tip.label)
		
		tips.unlinked	<- as.data.table(read.csv2(file=paste(indir,'TN_correct2.csv',sep='/'), stringsAsFactors=FALSE))
	}
	if(1)
	{
		infile			<- "TNTP_PRRT_FastTree.newick"
		ph				<- read.tree(paste(indir,infile,sep='/'))
		ph 				<- ladderize( ph )
		ph.node.bs						<- as.numeric( ph$node.label )		
		ph.node.bs[is.na(ph.node.bs)]	<- 0		
		ph$node.label					<- ph.node.bs 	
		ph$tip.label	<- toupper(ph$tip.label)	
		
		tips.unlinked	<- as.data.table(read.csv2(file=paste(indir,'TN_correct.csv',sep='/'), stringsAsFactors=FALSE))
	}
	#	plot tree
	pdf(file=paste(outdir,'/',gsub('\\.newick','\\.pdf',infile),sep=''), w=10, h=150)
	plot(ph, show.tip.label=FALSE, show.node.label=TRUE, cex=0.4)
	dev.off()
	#	read TP 	
	tips.linked		<- as.data.table(read.csv2(file=paste(indir,'TP.csv',sep='/'), stringsAsFactors=FALSE))
	setnames(tips.linked, c('seq1','seq2'), c('FASTASampleCode','t.FASTASampleCode'))
	#	read TN
	
	setnames(tips.unlinked, c('seq1','seq2'), c('FASTASampleCode','t.FASTASampleCode'))
	#	get unique tips
	tmp				<- tips.unlinked[, which(FASTASampleCode<t.FASTASampleCode)]
	z				<- tips.unlinked[tmp, t.FASTASampleCode]
	set(tips.unlinked,tmp,'t.FASTASampleCode',tips.unlinked[tmp, FASTASampleCode])
	set(tips.unlinked,tmp,'FASTASampleCode',z)
	setkey(tips.unlinked, FASTASampleCode, t.FASTASampleCode)
	tips.unlinked	<- unique(tips.unlinked)	#they were unique :-)	
	tmp				<- tips.linked[, which(FASTASampleCode<t.FASTASampleCode)]
	z				<- tips.linked[tmp, t.FASTASampleCode]
	set(tips.linked,tmp,'t.FASTASampleCode',tips.linked[tmp, FASTASampleCode])
	set(tips.linked,tmp,'FASTASampleCode',z)
	setkey(tips.linked, FASTASampleCode, t.FASTASampleCode)
	tips.linked		<- unique(tips.linked)	#they were also unique :-)
	#	some TP and TNs are not in the fasta file..
	tips.linked		<- subset(tips.linked, FASTASampleCode%in%ph$tip.label & t.FASTASampleCode%in%ph$tip.label)
	tips.unlinked	<- subset(tips.unlinked, FASTASampleCode%in%ph$tip.label & t.FASTASampleCode%in%ph$tip.label)
	
	#	get BRLs and MRCAS 
	ph.dist				<- cophenetic.phylo(ph)
	ph.mrca				<- mrca(ph)
	rownames(ph.dist)	<- toupper(rownames(ph.dist))
	colnames(ph.dist)	<- toupper(colnames(ph.dist))
	rownames(ph.mrca)	<- toupper(rownames(ph.mrca))
	colnames(ph.mrca)	<- toupper(colnames(ph.mrca))
	tmp					<- tips.linked[, list(BRL=ph.dist[FASTASampleCode,t.FASTASampleCode]) ,by=c('FASTASampleCode','t.FASTASampleCode')]
	tips.linked			<- merge(tmp, tips.linked, by=c('FASTASampleCode','t.FASTASampleCode'))
	tmp					<- tips.unlinked[, list(BRL=ph.dist[FASTASampleCode,t.FASTASampleCode]) ,by=c('FASTASampleCode','t.FASTASampleCode')]
	tips.unlinked		<- merge(tmp, tips.unlinked, by=c('FASTASampleCode','t.FASTASampleCode'))
	tmp					<- tips.linked[, list(MRCA= ph.mrca[FASTASampleCode,t.FASTASampleCode]), by=c('FASTASampleCode','t.FASTASampleCode')]
	tips.linked			<- merge(tips.linked, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))
	tmp					<- tips.unlinked[, list(MRCA= ph.mrca[FASTASampleCode,t.FASTASampleCode]) ,by=c('FASTASampleCode','t.FASTASampleCode')]
	tips.unlinked		<- merge(tmp, tips.unlinked, by=c('FASTASampleCode','t.FASTASampleCode'))
	tips.linked[, BSatMRCA:= ph$node.label[ tips.linked[, MRCA-Ntip(ph)] ]]
	tmp					<- tips.linked[, {													
									tmp			<- Ancestors(ph, MRCA)		
									anc.bs		<- ph$node.label[ tmp-Ntip(ph) ]
									list(	BSbelowMRCA= ifelse(length(tmp),max(anc.bs),0)		)													
								}, by=c('FASTASampleCode','t.FASTASampleCode')]
	tips.linked			<- merge(tips.linked, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))	
	tmp					<- tips.linked[,	list(BS=max(BSbelowMRCA,BSatMRCA)),		by=c('FASTASampleCode','t.FASTASampleCode')]
	tips.linked			<- merge(tips.linked, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))		
	tips.unlinked[, BSatMRCA:= ph$node.label[ tips.unlinked[, MRCA-Ntip(ph)] ]]
	tmp					<- tips.unlinked[, {													
				tmp			<- Ancestors(ph, MRCA)		
				anc.bs		<- ph$node.label[ tmp-Ntip(ph) ]
				list(	BSbelowMRCA= ifelse(length(tmp),max(anc.bs),0)		)													
			}, by=c('FASTASampleCode','t.FASTASampleCode')]
	tips.unlinked		<- merge(tips.unlinked, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))
	tmp					<- tips.unlinked[,	list(BS=max(BSbelowMRCA,BSatMRCA)),		by=c('FASTASampleCode','t.FASTASampleCode')]
	tips.unlinked		<- merge(tips.unlinked, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))	
	#	save	
	file			<- paste(outdir,'/',gsub('\\.newick','_TPTN.R',infile),sep='')
	cat('\nsave to file',file)
	save(tips.linked, tips.unlinked, ph, ph.dist, ph.mrca, file=file)
	
	set(tips.linked, NULL, 'BS', tips.linked[, 100*BS])
	set(tips.linked, NULL, 'BSatMRCA', tips.linked[, 100*BSatMRCA])
	set(tips.linked, NULL, 'BSbelowMRCA', tips.linked[, 100*BSbelowMRCA])
	set(tips.unlinked, NULL, 'BS', tips.unlinked[, 100*BS])
	set(tips.unlinked, NULL, 'BSatMRCA', tips.unlinked[, 100*BSatMRCA])
	set(tips.unlinked, NULL, 'BSbelowMRCA', tips.unlinked[, 100*BSbelowMRCA])
	
	
	tips.linked[, TYPE:='LINKED']
	tips.unlinked[, TYPE:='UNLINKED']
	tips.lul		<- rbind(tips.linked, tips.unlinked)
	ggplot(tips.lul, aes(x=BRL)) + geom_histogram(binwidth=0.002) + scale_x_continuous(lim=c(0,0.1)) + facet_grid(~TYPE)
	ggsave(file=paste(outdir,'/',gsub('\\.newick','_LinkedUnlinkedBRLHisto.pdf',infile),sep=''), w=7, h=4)
	ggplot(tips.lul, aes(x=BS)) + geom_histogram(binwidth=0.02) + facet_grid(~TYPE)
	ggsave(file=paste(outdir,'/',gsub('\\.newick','_LinkedUnlinkedBSHisto.pdf',infile),sep=''), w=7, h=4)
	
	heat			<- as.data.table(expand.grid(BS_CUT=seq(0,1,0.01), BRL_CUT=seq(0, 0.2, 0.005), DATA='Dataset Surrogate Linked', DECLARE='Declare Linked'))
	tmp				<- heat[, list( PC= nrow(subset(tips.linked, BS>=BS_CUT & BRL<=BRL_CUT)) / nrow(tips.linked) ), by=c('BS_CUT','BRL_CUT')]
	heat			<- merge(heat, tmp, by=c('BS_CUT','BRL_CUT'))
	
	ggplot(subset(heat, DATA=='Dataset Surrogate Linked' & DECLARE=='Declare Linked' & BS_CUT>=0.5), aes(x=100*BS_CUT, y=100*BRL_CUT, fill=cut(PC, breaks=c(-0.01, 0.01, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.01), labels=c('<1%','1-10%','10-30%','30-50%','50-60%','60-70%','70-80%','80-90%','90-95%','>95%')))) + geom_tile(colour = "white") +
			scale_fill_brewer(palette='PiYG') + theme_bw() +
			scale_x_continuous(expand=c(0,0), breaks=c(50,60,70,80,90,95,100)) + scale_y_continuous(expand=c(0,0), breaks=c(0,2,4,6,8,10,15,20)) +
			labs(x='bootstrap support threshold\n(%)', y='single linkage branch length threshold\n(%)', fill='Linked sequences\nthat are in same cluster') +
			#facet_grid(DECLARE~DATA) +
			theme_bw() + theme(legend.position='bottom') + guides(fill = guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/',infile,'_Bias.pdf',sep=''), w=5, h=6)
	
		
}
######################################################################################
project.hivc.clustering.forStephane.onUK.FastTree<- function()
{	
	verbose		<- 1
	resume		<- 1
	patient.n	<- 15700; 	thresh.brl<- 0.096; 	thresh.bs<- 0.8;	opt.brl<- "dist.brl.casc" 
	indircov	<- paste(DATA,"derived",sep='/')
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"
	indir		<- paste(DATA,"tmp",sep='/')	
	outdir		<- paste(DATA,"brl_surrlink",sep='/')
	
	infile		<- "~/duke/2015_OptimalClusteringThresholds/TNTP_PRRT_FastTree.nwk"
	ph			<- read.tree(infile)
	ph 			<- ladderize( ph )
	ph.node.bs						<- as.numeric( ph$node.label )		
	ph.node.bs[is.na(ph.node.bs)]	<- 0
	ph$node.label					<- ph.node.bs
	stat.fun						<- hivc.clu.min.transmission.cascade
	#stat.fun						<- max
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)	
	save(ph, dist.brl,ph.node.bs, file=gsub('\\.nwk','_PRECLU\\.R',infile))
	
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=0.8, thresh.brl=0.04, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	print(clustering)		
	
	
	#read bootstrap support values
	thresh.bs						<- 0.9
	ph.node.bs[c(13,15,27,41,43)]	<- thresh.bs-0.05*seq(0.01,length.out=5)
	ph.node.bs[ph.node.bs==1]		<- seq_along(which(ph.node.bs==1))*0.005 + 0.7
	ph$node.label					<- ph.node.bs
	#read patristic distances
	stat.fun						<- hivc.clu.min.transmission.cascade
	#stat.fun						<- max
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)
	print(dist.brl)
	thresh.brl						<- quantile(dist.brl,seq(0.1,1,by=0.05))["100%"]
	print(quantile(dist.brl,seq(0.1,0.5,by=0.05)))
	print(thresh.brl)
	#produce clustering 
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	print(clustering)		
	hivc.clu.plot(ph, clustering[["clu.mem"]], highlight.edge.of.tiplabel=c("TN_","TP_"), highlight.edge.of.tiplabel.col= c("red","blue") )
	#produce some tip states
	ph.tip.state<- rep(1:20, each=ceiling( length( ph$tip.label )/20 ))[seq_len(Ntip(ph))]
	states		<- data.table(state.text=1:20, state.col=rainbow(20))
	states		<- states[ph.tip.state,]
	hivc.clu.plot.tiplabels( seq_len(Ntip(ph)), matrix(states[,state.text], nrow=1,ncol=nrow(states)), matrix(states[,state.col], nrow=1,ncol=nrow(states)), cex=0.4, adj=c(-1,0.5) )
	
}
######################################################################################
project.hivc.clustering.forStephane.examl<- function()
{
	require(ape)
	
	dir.name	<- DATA  	
	indir		<- paste(dir.name,'hue',sep='/')
	infile		<- "TNTP_PRRT_withoutDoublon_withoutDRM"
	insignat	<- "151023"
	
	if(0)
	{
		seq.PROT.RT	<- read.dna(paste(indir, '/', infile, '.fasta', sep=''), format='fasta')
		file		<- paste(indir, '/', infile, '_', insignat, '.R', sep='')
		save(seq.PROT.RT, file=file)
	}	
	#	ExaML bootstrap args
	bs.from		<- 0
	bs.to		<- 0
	bs.n		<- 500
	outdir		<- indir
	args.parser	<- "-m DNA"
	args.examl	<- "-f d -D -m GAMMA"	#	 -- this is the default that worked in 24 hours	
	cmd			<- hivc.cmd.examl.bootstrap(indir, infile, insignat, insignat, bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, opt.bootstrap.by="codon", args.parser=args.parser, args.examl=args.examl, tmpdir.prefix="examl")					
	
	
	dummy		<- lapply(cmd, function(x)
			{				
				#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
				x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=50, hpc.q="pqeelab", hpc.mem="5500mb", hpc.nproc=1)
				signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
				outfile	<- paste("exa",signat,sep='.')
				#cat(x)
				hivc.cmd.hpccaller(outdir, outfile, x)
				Sys.sleep(1)
			})	
}
######################################################################################
project.hivc.clustering.forStephane.onNL<- function()
{	
	verbose		<- 1
	resume		<- 1
	patient.n	<- 15700; 	thresh.brl<- 0.096; 	thresh.bs<- 0.8;	opt.brl<- "dist.brl.casc" 
	indircov	<- paste(DATA,"derived",sep='/')
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"
	indir		<- paste(DATA,"tmp",sep='/')	
	outdir		<- paste(DATA,"brl_surrlink",sep='/')
	#
	# get clusters for No Recombination + No Drug resistance mutations + No short sequences, single linkage criterion		
	#						
	infile			<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500"			
	insignat		<- "Wed_Dec_18_11/37/00_2013"
	#
	argv			<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	nsh.clu.pre		<- hivc.prog.get.clustering.precompute()
	ph				<- nsh.clu.pre$ph
	#	read TP tips and TN tips
	tips.linked		<- subset(nsh.clu.pre$bs.linked.bypatient, select=c('tip1','tip2','bs','PosSeqT.diff'))
	tips.unlinked	<- do.call('rbind',nsh.clu.pre$unlinked.bytime)
	tips.unlinked	<- subset(tips.unlinked, select=c('query.FASTASampleCode','FASTASampleCode'))
	setnames(tips.linked, c('tip1','tip2'), c('FASTASampleCode','t.FASTASampleCode'))
	setnames(tips.unlinked, c('query.FASTASampleCode','FASTASampleCode'), c('FASTASampleCode','t.FASTASampleCode'))
	#	get unique tips
	tmp				<- tips.unlinked[, which(FASTASampleCode<t.FASTASampleCode)]
	z				<- tips.unlinked[tmp, t.FASTASampleCode]
	set(tips.unlinked,tmp,'t.FASTASampleCode',tips.unlinked[tmp, FASTASampleCode])
	set(tips.unlinked,tmp,'FASTASampleCode',z)
	setkey(tips.unlinked, FASTASampleCode, t.FASTASampleCode)
	#
	#	get branch lengths of TP 
	#	get max bootstrap from mrca to root
	#
	indir.dists		<- paste(DATA,"fisheretal_data",sep='/')
	if(0)
	{
		tmp				<- subset(tips.linked, select=c('FASTASampleCode','t.FASTASampleCode'))
		ph.dist			<- cophenetic.phylo(ph)
		tmp				<- tmp[, list(BRL=ph.dist[FASTASampleCode,t.FASTASampleCode]) ,by=c('FASTASampleCode','t.FASTASampleCode')]
		tips.linked		<- merge(tmp, tips.linked, by=c('FASTASampleCode','t.FASTASampleCode'))
		setnames(tips.linked, c('bs'), c('BS'))
		file			<- paste(outdir,'/',infile,'_TPbrl.R',sep='')
		cat('\nsave to file',file)
		save(tips.linked, file=file)		
		load(file)
	}
	#
	#	get branch lengths of TN
	#	get max bootstrap from mrca to root
	#
	if(1)
	{
		tips.unlinked.d	<- copy(tips.unlinked)
		tmp				<- subset(tips.unlinked.d, select=c('FASTASampleCode','t.FASTASampleCode'))
		tmp				<- tmp[, list(BRL=ph.dist[FASTASampleCode,t.FASTASampleCode]) ,by=c('FASTASampleCode','t.FASTASampleCode')]
		tips.unlinked.d	<- merge(tmp, tips.unlinked.d, by=c('FASTASampleCode','t.FASTASampleCode'))		
		
		ph					<- nsh.clu.pre$ph
		ph.mrca				<- mrca(ph)		
		tmp					<- tips.unlinked.d[, list(MRCA= ph.mrca[FASTASampleCode,t.FASTASampleCode]), by=c('FASTASampleCode','t.FASTASampleCode')]
		tips.unlinked.d		<- merge(tips.unlinked.d, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))
		tips.unlinked.d[, BSatMRCA:= ph$node.label[ tips.unlinked.d[, MRCA-Ntip(ph)] ]]
		tmp					<- tips.unlinked.d[, {													
					tmp			<- Ancestors(ph, MRCA)		
					anc.bs		<- ph$node.label[ tmp-Ntip(ph) ]
					list(	BSbelowMRCA= ifelse(length(tmp),max(anc.bs),0)		)													
				}, by=c('FASTASampleCode','t.FASTASampleCode')]
		tips.unlinked.d		<- merge(tips.unlinked.d, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))
		tmp					<- tips.unlinked.d[,	list(BS=max(BSbelowMRCA,BSatMRCA)),		by=c('FASTASampleCode','t.FASTASampleCode')]
		tips.unlinked.d		<- merge(tips.unlinked.d, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))	
		
		file			<- paste(outdir,'/',infile,'_TNbrl.R',sep='')
		cat('\nsave to file',file)
		save(tips.unlinked.d, file=file)		
	}	
	if(0)	
	{
		#get Shankarappa TP data set
		require(gtools)
		file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/brl_surrlink/set1_Shankarappa.fasta'
		shk		<- read.dna(file=file, format='fasta')
		shk.brl	<- dist.dna(shk, model='raw', as.matrix=T)
		shk.df	<- data.table(TAXON1= rownames(shk))
		shk.df[, PATIENT:=shk.df[, regmatches(TAXON1,regexpr('[0-9]+V',TAXON1))]] 
		set(shk.df, NULL, 'PATIENT', shk.df[ , substr(PATIENT,1,nchar(PATIENT)-1) ])
		shk.df	<- shk.df[, 
				{
					z<- combinations(length(TAXON1),2)
					list(TAXON1=TAXON1[z[,1]], TAXON2=TAXON1[z[,2]])
				}, by='PATIENT']
		shk.df	<- shk.df[, list(BRL= shk.brl[TAXON1,TAXON2]), by=c('TAXON1','TAXON2')]
		file	<- paste(outdir,'/Shankarappa_TPbrl.R',sep='')
		save(shk.df, file=file)
		#
		#		analysis of Shankarappa TP data
		#
		shk.df[, quantile(BRL, p=c(0.95,0.975,0.99))]
		#	95%      	97.5%        	99% 
		#	0.06250000 	0.06854839 		0.07459677 
		ggplot(shk.df, aes(x=BRL)) +
				scale_x_continuous(breaks=seq(0,0.3,0.01), limit=c(0,0.15)) +
				geom_histogram() + theme_bw() + labs(x='genetic distance among within-host Shankarappa sequences\n(substitutions/site)')
		ggsave(file= paste(outdir,'/Shankarappa_TPbrl.pdf',sep=''), w=8, h=5)		
	}
	if(0)
	{
		tips.linked[, TYPE:='LINKED']
		tips.unlinked[, TYPE:='UNLINKED']
		file			<- paste(outdir,'/',infile,'_TNTP.R',sep='')
		cat('\nsave to file',file)
		save(tips.unlinked, tips.linked, ph, file=file)
		#		
		tips.lul		<- rbind( subset(tips.linked, select=c(FASTASampleCode, t.FASTASampleCode, BRL, BS, TYPE)), subset(tips.unlinked, select=c(FASTASampleCode, t.FASTASampleCode, BRL, BS, TYPE)) )
		ggplot(tips.lul, aes(x=BRL)) + geom_histogram(binwidth=0.002) + scale_x_continuous(lim=c(0,0.1)) + facet_grid(~TYPE)
		ggsave(file=paste(outdir,'/',infile,'_LinkedUnlinkedBRLHisto.pdf',sep=''), w=7, h=4)
		ggplot(tips.lul, aes(x=BS)) + geom_histogram(binwidth=0.02) + facet_grid(~TYPE)
		ggsave(file=paste(outdir,'/',infile,'_LinkedUnlinkedBSHisto.pdf',sep=''), w=7, h=4)
		ggplot(tips.lul, aes(x=BS)) + geom_histogram(binwidth=0.02) + coord_cartesian(ylim=c(0,2.1e3)) + facet_grid(~TYPE) 
		ggsave(file=paste(outdir,'/',infile,'_LinkedUnlinkedBSHisto_Zoom.pdf',sep=''), w=7, h=4)		
		#
		#		analysis of ATHENA TP / TN data
		#
		subset(tips.linked, BS==0)[, quantile(BRL, p=c(0.95,0.975,0.99))]
		#       95%      	97.5%       99% 
		#		0.03780011 	0.08691466 	0.21071687
		tips.unlinked.bs0	<- subset(tips.unlinked, BS==0)
		tips.unlinked.bs0[, quantile(BRL, p=c(0.05,0.025,0.01, 0.001, 0.0001))]
		#	    5%       	2.5%         	1%       	0.1%      	0.01% 
		#		0.16053774 	0.15067962 		0.13792314 	0.10836011 	0.04658164
		#		this is **much** larger than SH calculated, but still the TN pairs are the majority  
	
		brl					<- rbind( tips.unlinked.bs0, subset(tips.linked, BS==0) )	
		
		ggplot(subset(tips.linked, BS==0), aes(x=BRL)) +
				scale_x_continuous(breaks=seq(0,0.3,0.01), limit=c(0,0.15)) +
				geom_histogram() + theme_bw() + labs(x='genetic distance among within-host sequences\n(substitutions/site)')
		ggsave(file= paste(outdir,'/',infile,'_TPbrl.pdf',sep=''), w=8, h=5)
		
		ggplot(brl, aes(x=BRL, group=TYPE, fill=TYPE)) +
				scale_x_continuous(breaks=seq(0,0.6,0.05)) +
				scale_y_continuous(breaks=c(10,100,1000,10000,20000), minor_breaks=NULL, trans ="log1p") +
				scale_fill_brewer(palette='Set1') +
				geom_histogram(alpha=0.4) + theme_bw() + labs(x='genetic distance among surrogate pairs\n(substitutions/site)')
		ggsave(file= paste(outdir,'/',infile,'_TPTNbrl.pdf',sep=''), w=8, h=6)
		
		ggplot(brl, aes(x=BRL, group=TYPE, fill=TYPE)) +
				scale_x_continuous(breaks=seq(0,0.3,0.01), limit=c(0,0.2)) +
				geom_histogram() + theme_bw() + labs(x='genetic distance among surrogate pairs\n(substitutions/site)')
	}	
	#
	#	get heatmaps
	#
	if(0)
	{
		tips.linked.d	<- copy(tips.linked)
		#tips.linked.d	<- subset(tips.linked, BRL_TREE_ID==0)
		heat			<- as.data.table(expand.grid(BS_CUT=seq(0,1,0.01), BRL_CUT=seq(0, 0.2, 0.005), DATA='Dataset Surrogate Linked', DECLARE='Declare Linked'))
		tmp				<- heat[, list( PC= nrow(subset(tips.linked, BS>=BS_CUT & BRL<=BRL_CUT)) / nrow(tips.linked) ), by=c('BS_CUT','BRL_CUT')]
		heat			<- merge(heat, tmp, by=c('BS_CUT','BRL_CUT'))
		heat.df			<- copy(heat)
		
		ggplot(subset(heat.df, DATA=='Dataset Surrogate Linked' & DECLARE=='Declare Linked' & BS_CUT>=0.5), aes(x=100*BS_CUT, y=100*BRL_CUT, fill=cut(PC, breaks=c(-0.01, 0.01, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.01), labels=c('<1%','1-10%','10-30%','30-50%','50-60%','60-70%','70-80%','80-90%','90-95%','>95%')))) + geom_tile(colour = "white") +
				scale_fill_brewer(palette='PiYG') + theme_bw() +
				scale_x_continuous(expand=c(0,0), breaks=c(50,60,70,80,90,95,100)) + scale_y_continuous(expand=c(0,0), breaks=c(0,2,4,6,8,10,15,20)) +
				labs(x='bootstrap support threshold\n(%)', y='single linkage branch length threshold\n(%)', fill='Linked sequences\nthat are in same cluster') +
				#facet_grid(DECLARE~DATA) +
				theme_bw() + theme(legend.position='bottom') + guides(fill = guide_legend(ncol=4))
		ggsave(file=paste(outdir,'/',infile,'_Bias.pdf',sep=''), w=5, h=6)
		
		
		heat			<- as.data.table(expand.grid(BS_CUT=seq(0,1,0.01), BRL_CUT=seq(0, 0.2, 0.005), DATA='Dataset Unlinked', DECLARE='Declare Linked'))
		tmp				<- heat[, list( PC= nrow(subset(tips.unlinked.d, BS>=BS_CUT & BRL<=BRL_CUT)) / nrow(tips.unlinked.d) ), by=c('BS_CUT','BRL_CUT')]
		heat			<- merge(heat, tmp, by=c('BS_CUT','BRL_CUT'))
		heat.df			<- rbind(heat.df, heat)
		
		heat			<- copy(heat.df)
		set(heat, NULL, 'DECLARE', 'Declare Unlinked')
		set(heat, NULL, 'PC', heat[, 1-PC])
		heat.df			<- rbind(heat.df, heat)
		
		ggplot(heat.df, aes(x=100*BS_CUT, y=BRL_CUT, fill=cut(PC, breaks=c(-0.01, 0.01, 0.05, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.01)))) + geom_tile(colour = "white") +
							scale_fill_brewer(palette='PiYG') + theme_bw() +
							scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
							labs(x='bootstrap support threshold\n(%)', y='single linkage branch length threshold\n(%)', fill='%') +
							facet_grid(DECLARE~DATA) +
							theme_bw() + theme(legend.position='bottom') + guides(fill = guide_legend(ncol=5))
		
		
		
		
		ggplot(subset(heat.df, DATA=='Dataset Unlinked' & DECLARE=='Declare Unlinked' & BS_CUT>=0.5), aes(x=100*BS_CUT, y=BRL_CUT, fill=cut(PC, include.lowest=0, breaks=c(-1,0.99, 0.999, 0.9999, 0.99995, 0.99999, 1.00, 1.01)))) + geom_tile(colour = "white") +
				scale_fill_brewer(palette='PiYG') + theme_bw() +
				scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
				labs(x='bootstrap support threshold\n(%)', y='single linkage branch length threshold\n(%)', fill='%') +
				facet_grid(DECLARE~DATA) +
				theme_bw() + theme(legend.position='bottom') + guides(fill = guide_legend(ncol=5))
		
		
		ggplot(subset(heat.df, DATA=='Dataset Unlinked' & DECLARE=='Declare Linked' & BS_CUT>=0.5), aes(x=100*BS_CUT, y=BRL_CUT, fill=cut(PC*nrow(tips.unlinked.d), breaks=c(100,80,60,40,20,10,5,1,-1), labels=c('<2','2-5','6-10','11-20','21-40','41-60','61-80','81-100')))) + 
				geom_tile(colour = "white") +
				scale_fill_brewer(palette='PiYG') + 
				theme_bw() +
				scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
				labs(x='bootstrap support threshold\n(%)', y='single linkage branch length threshold\n(%)', fill='%') +
				facet_grid(DECLARE~DATA) +
				theme_bw() + theme(legend.position='bottom') + guides(fill = guide_legend(ncol=5))
		
		
		
		#tips.linked.d	<- subset(tips.linked, BRL_TREE_ID==0)
		heat			<- as.data.table(expand.grid(BS_CUT=seq(0,1,0.05), BRL_CUT=seq(0, 0.2, 0.005)))
		tmp				<- heat[, {
										linked.declared.linked		<- nrow(subset(tips.linked, BS>=BS_CUT & BRL<=BRL_CUT))
										linked.declared.unlinked	<- nrow(tips.linked)-linked.declared.linked
										unlinked.declared.linked	<- nrow(subset(tips.unlinked.d, BS>=BS_CUT & BRL<=BRL_CUT))
										unlinked.declared.unlinked	<- nrow(tips.unlinked.d)-unlinked.declared.linked
										list( 	LdL=linked.declared.linked,
												LdU=linked.declared.unlinked,
												UdL=unlinked.declared.linked,
												UdU=unlinked.declared.unlinked,
												FDR_L= unlinked.declared.linked/(unlinked.declared.linked+linked.declared.linked),
											  	FDR_U= linked.declared.unlinked/(linked.declared.unlinked+unlinked.declared.unlinked)	)						
									}, by=c('BS_CUT','BRL_CUT')]
		heat			<- merge(heat, tmp, by=c('BS_CUT','BRL_CUT'))
		set(heat, heat[, which(is.nan(FDR_L))], 'FDR_L', 0)
		set(heat, heat[, which(is.nan(FDR_U))], 'FDR_U', 0)
		tmp				<- melt(heat, measure.vars=c('FDR_L','FDR_U'))
		ggplot(tmp, aes(x=100*BS_CUT, y=BRL_CUT, fill=cut(value, breaks=c(-0.01, 0.01, 0.05, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.01)))) + 
				geom_tile(colour = "white") +
				scale_fill_brewer(palette='PiYG') + theme_bw() +
				scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
				labs(x='bootstrap support threshold\n(%)', y='single linkage branch length threshold\n(%)', fill='%') +
				facet_grid(~variable) +
				theme_bw() + theme(legend.position='bottom') + guides(fill = guide_legend(ncol=5))
		
		
		
	}
	#
	#	get false discovery rate
	#
	if(0)
	{
		setkey(heat.df, BS_CUT, BRL_CUT)
		fdr.df			<- subset(unique(heat.df), select=c(BS_CUT, BRL_CUT))		
		tntp.df			<- fdr.df[, {
					clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=BS_CUT, thresh.brl=BRL_CUT, dist.brl=nsh.clu.pre$dist.brl.casc, nodesupport=ph$node.label, retval="all")
					tmp			<- data.table(CLU_N= clustering$size.tips)
					tmp[, TP:= CLU_N-1]
					tmp[, TN:= choose(CLU_N-1,2)]
					stopifnot( tmp[, !any(is.na(TP))], tmp[, !any(is.na(TN))])
					list( SUM_TP= tmp[, sum(TP)], SUM_TN= tmp[, sum(TN)])			
				}, by=c('BS_CUT','BRL_CUT')]
		fdr.df			<- merge(heat.df, tntp.df, by=c('BS_CUT','BRL_CUT'))		
		set( fdr.df, NULL, 'DATA', fdr.df[, gsub(' ','xx',DATA)]) 
		set( fdr.df, NULL, 'DECLARE', fdr.df[, gsub(' ','xx',DECLARE)]) 
		fdr.df				<- dcast.data.table(fdr.df, BS_CUT+BRL_CUT+SUM_TP+SUM_TN~DATA+DECLARE, value.var='PC')
		fdr.df[, ETL:= DatasetxxSurrogatexxLinked_DeclarexxLinked*SUM_TP]
		fdr.df[, EFL:= DatasetxxUnlinked_DeclarexxLinked*SUM_TN]		
		fdr.df[, ETU:= DatasetxxSurrogatexxLinked_DeclarexxUnlinked*SUM_TP]
		fdr.df[, EFU:= DatasetxxUnlinked_DeclarexxUnlinked*SUM_TN]
		fdr.df[, EFDRL:= EFL / (ETL+EFL)]
		fdr.df[, EFDRU:= EFU / (ETU+EFU)]
		fdr.df			<- melt(fdr.df, id.vars=c('BS_CUT','BRL_CUT','SUM_TP','SUM_TN'), measure.vars=c('EFDRL','EFDRU'), value.name='FDR', variable.name='DECLARE')
		set(fdr.df, fdr.df[, which(DECLARE=='EFDRL')], 'DECLARE', 'Declare Linked')
		set(fdr.df, fdr.df[, which(DECLARE=='EFDRU')], 'DECLARE', 'Declare Unlinked')
		
		ggplot(fdr.df, aes(x=100*BS_CUT, y=100*BRL_CUT, fill=cut(FDR, breaks=c(-0.01, 0.01, 0.05, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.01)))) + geom_tile(colour = "white") +
				scale_fill_brewer(palette='PiYG') + theme_bw() +
				scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
				labs(x='bootstrap support threshold\n(%)', y='single linkage branch length threshold\n(%)', fill='%') +
				facet_grid(DECLARE~.) +
				theme_bw() + theme(legend.position='bottom') + guides(fill = guide_legend(ncol=5))
		
	}
}
######################################################################################
project.hivc.clustering.compare.NoDR.to.NoRecombNoDR.to.NoShort<- function()
{	
	verbose		<- 1
	resume		<- 1
	patient.n	<- 15700; 	thresh.brl<- 0.096; 	thresh.bs<- 0.8;	opt.brl<- "dist.brl.casc" 
	indircov	<- paste(DATA,"derived",sep='/')
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"
	indir		<- paste(DATA,"tmp",sep='/')
	#
	# get clusters for No Drug resistance mutations, single linkage criterion		
	#					
	infile			<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs500"			
	insignat		<- "Thu_Aug_01_17/05/23_2013"
	#
	argv			<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	ndr.clu.pre		<- hivc.prog.get.clustering.precompute()
	#
	argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume)
	argv			<<- unlist(strsplit(argv,' '))
	ndr.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=ndr.clu.pre)
	#			
	argv			<<- hivc.cmd.clustering(indir, infile, insignat, opt.brl="dist.brl.casc", thresh.brl, thresh.bs, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	ndr.clu			<- hivc.prog.get.clustering()
	#
	# get clusters for No Recombination + No Drug resistance mutations, single linkage criterion		
	#						
	infile			<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences_examlbs500"			
	insignat		<- "Fri_Nov_01_16/07/23_2013"
	#
	argv			<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	nrc.clu.pre		<- hivc.prog.get.clustering.precompute()
	#
	argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume)
	argv			<<- unlist(strsplit(argv,' '))
	nrc.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=nrc.clu.pre)
	#			
	argv			<<- hivc.cmd.clustering(indir, infile, insignat, opt.brl="dist.brl.casc", thresh.brl, thresh.bs, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	nrc.clu			<- hivc.prog.get.clustering()
	#
	# get clusters for No Recombination + No Drug resistance mutations + No short sequences, single linkage criterion		
	#						
	infile			<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500"			
	insignat		<- "Wed_Dec_18_11/37/00_2013"
	#
	argv			<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	nsh.clu.pre		<- hivc.prog.get.clustering.precompute()
	#
	resume			<- 0
	argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume)
	argv			<<- unlist(strsplit(argv,' '))
	nsh.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=nsh.clu.pre, with.plot=0)
	#			
	argv			<<- hivc.cmd.clustering(indir, infile, insignat, opt.brl="dist.brl.casc", thresh.brl, thresh.bs, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	nsh.clu			<- hivc.prog.get.clustering()		
	#
	argv			<<- hivc.cmd.clustering.msm(indir, infile, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=resume)
	argv			<<- unlist(strsplit(argv,' '))		
	nsh.msm			<- hivc.prog.get.clustering.MSM()		
	
	seq.n	<- c( Ntip( ndr.clu.pre$ph ), Ntip( nrc.clu.pre$ph ), Ntip( nsh.clu.pre$ph ))
	
	bs.small		<- 0.5
	bs.ok			<- 0.7
	bs.large		<- 0.95
	bs.su			<- matrix( c(
								c( length( which( ndr.clu.pre$ph.node.bs<bs.small ) ), length( which( nrc.clu.pre$ph.node.bs<bs.small ) ), length( which( nsh.clu.pre$ph.node.bs<bs.small ) ) ),
								c( length( which( ndr.clu.pre$ph.node.bs>bs.ok ) ), length( which( nrc.clu.pre$ph.node.bs>bs.ok ) ), length( which( nsh.clu.pre$ph.node.bs>bs.ok ) ) ),
								c( length( which( ndr.clu.pre$ph.node.bs>bs.large ) ), length( which( nrc.clu.pre$ph.node.bs>bs.large ) ), length( which( nsh.clu.pre$ph.node.bs>bs.large ) ) ) 
								), byrow=1, nrow=3, ncol=3, dimnames=list(c(),c("nDR","nDRRC","nDRRCSH"))	)
	bs.su			<- cbind( as.data.table(bs.su), info=c('bs.small','bs.ok','bs.large'))
	bs.su[,nDRRC-nDRRCSH] 
	#
	#	compare distribution of bootstrap values
	#
	dir.name	<- paste( DATA,'tmp',sep='/' )
	file		<- paste( dir.name, 'ATHENA_2013_03_compare_-DR-RC-SH_seqdatasets.pdf', sep='/')
	cat(paste("\nplot to file",file))
	pdf(file, 5, 5)
	hist( ndr.clu.pre$ph.node.bs, main='', xlab='bootstrap score' )
	hist( nrc.clu.pre$ph.node.bs, border="blue", add=1 )
	hist( nsh.clu.pre$ph.node.bs, border="red", add=1 )
	legend( "topright", border=NA, bty='n', fill= c('black','blue','red'), legend=c("-DR","-DR -RC","-DR -RC -SH" ))
	dev.off()
	#	compare TP, TN on a common set, use the one from -DR -RC -SH
	#
	# 	need to reset ph.unlinked[[i]]	ph.unlinked.info	ph.linked		miss NodeIdx in ph.unlinked.info
	ref						<- nsh.clu.pre
	ref$ph.unlinked.info	<- ref$ph.unlinked.info[, -c(2,3), with=0]	
	ref$ph.linked			<- ref$ph.linked[,-3,with=0]
	#
	#	-DR
	#
	run						<- ndr.clu.pre
	run.tips				<- data.table(Node=seq_len(Ntip(run$ph)), FASTASampleCode=run$ph$tip.label, key="FASTASampleCode" )
	run$ph.unlinked.info	<- merge( run.tips, ref$ph.unlinked.info, by="FASTASampleCode" )	
	run$ph.linked			<- merge( run.tips, ref$ph.linked, by="FASTASampleCode" )	
	run$ph.unlinked			<- lapply(seq_along(ref$ph.unlinked),function(i)
		{
			ref$ph.unlinked[[i]]	<- ref$ph.unlinked[[i]][,-2,with=0]
			z						<- merge( run.tips, ref$ph.unlinked[[i]], by="FASTASampleCode")
			setkey(z, Node)
			z
		})
	names(run$ph.unlinked)	<- names(ref$ph.unlinked)
	tmp						<- data.table(FASTASampleCode=names(run$ph.unlinked), NodeIdx= seq_along(run$ph.unlinked) )
	run$ph.unlinked.info	<- merge(tmp, run$ph.unlinked.info, by="FASTASampleCode" )	
	#	re-run tptn	
	infile			<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs500"			
	insignat		<- "Thu_Aug_01_17/05/23_2013"
	argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=0)
	argv			<<- unlist(strsplit(argv,' '))
	ndr.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=run)
	#
	#	-DR -RC
	#
	run						<- nrc.clu.pre
	run.tips				<- data.table(Node=seq_len(Ntip(run$ph)), FASTASampleCode=run$ph$tip.label, key="FASTASampleCode" )
	run$ph.unlinked.info	<- merge( run.tips, ref$ph.unlinked.info, by="FASTASampleCode" )	
	run$ph.linked			<- merge( run.tips, ref$ph.linked, by="FASTASampleCode" )	
	run$ph.unlinked			<- lapply(seq_along(ref$ph.unlinked),function(i)
			{
				ref$ph.unlinked[[i]]	<- ref$ph.unlinked[[i]][,-2,with=0]
				z						<- merge( run.tips, ref$ph.unlinked[[i]], by="FASTASampleCode")
				setkey(z, Node)
				z
			})
	names(run$ph.unlinked)	<- names(ref$ph.unlinked)
	tmp						<- data.table(FASTASampleCode=names(run$ph.unlinked), NodeIdx= seq_along(run$ph.unlinked) )
	run$ph.unlinked.info	<- merge(tmp, run$ph.unlinked.info, by="FASTASampleCode" )
	#	re-run tptn	
	infile			<- "ATHENA_2013_03_NoRCDRAll+LANL_Sequences_examlbs500"			
	insignat		<- "Fri_Nov_01_16/07/23_2013"
	argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=0)
	argv			<<- unlist(strsplit(argv,' '))
	nrc.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=run)	
	#
	#	compare TP
	#
	tptn.cmp		<- lapply( c('0.8','0.95'), function(bs)
			{
				tmp	<- list( 	data.table( run='-DR',tp= as.numeric(ndr.clu.tptn$tp.by.all[bs, ]), bs=as.numeric(bs), brl=as.numeric(colnames(ndr.clu.tptn$tp.by.all)) ), 
						data.table( run='-DR -RC',tp= as.numeric(nrc.clu.tptn$tp.by.all[bs, ]), bs=as.numeric(bs), brl=as.numeric(colnames(nrc.clu.tptn$tp.by.all)) ),
						data.table( run='-DR -RC -SH',tp= as.numeric(nsh.clu.tptn$tp.by.all[bs, ]), bs=as.numeric(bs), brl=as.numeric(colnames(nsh.clu.tptn$tp.by.all)) )		)
				tmp	<- do.call('rbind',tmp)						
			})
	tptn.cmp		<- do.call('rbind',tptn.cmp)
	#
	file		<- paste( dir.name, 'ATHENA_2013_03_compare_-DR-RC-SH_tp.pdf', sep='/')
	cat(paste("\nplot to file",file))
	pdf(file, 5, 5)	
	ltys	<- c(3,1)
	cols	<- c('red','blue')
	runs	<- c('-DR -RC','-DR -RC -SH')
	bss		<- c(0.8, 0.95)
	plot(bty='n',type='n',1,1,xlab="Branch length cut-off",ylab="% of within patient seq in same cluster",xlim=c(0.02,0.125),ylim=range(c(1,tptn.cmp[,tp])))	
	dummy	<- sapply(seq_along(bss), function(j)
			{
				sapply(seq_along(runs), function(i)
						{
							lines( subset(tptn.cmp,run==runs[i] & bs==bss[j], )[, brl], subset(tptn.cmp,run==runs[i] & bs==bss[j], )[, tp], lty=ltys[i], col=cols[j] )			
						})				
			})
	legend("topright", bty='n', lty=ltys, legend=runs)
	legend("topleft", bty='n', border=NA, fill=cols, legend=paste('bootstrap= ',bss*100,'%',sep=''))
	dev.off()
	#
	#	compare FP
	#
	ndr.clu.tptn$fpn.by.sum
	nrc.clu.tptn$fpn.by.sum	
	nsh.clu.tptn$fpn.by.sum	
	#
}
######################################################################################
project.hivc.clustering.computeclusterstatistics.fordistbrl<- function(thresh, clusters, with.withinpatientclu=0, dist.brl="dist.brl.med")
{
	#select clusters for analysis
	clusters			<- clusters[[dist.brl]]		
	patient.nNL			<- 15700	#estimate from 2012 annual report	
	#length(df.seqinfo[,unique(Patient)])		
	patient.n			<- patient.nNL	 + length(which(substr(ph$tip.label,1,2)=="TN")) + length(which(substr(ph$tip.label,1,8)=="PROT+P51"))
	
	clusters.pairs		<- sapply(clusters,function(x){	table(x$size.tips)[1]  	})
	clusters.seqinclu	<- sapply(clusters,function(x){	sum(x$size.tips) / patient.n 	})
	clusters.NLseqinclu	<- sapply(clusters,function(x){	sum(x$size.tips) / patient.nNL 	})
	clusters.totalclu	<- sapply(clusters,function(x){	length(na.omit( unique(x$clu.mem) ))  	})	
	clusters.maxsize		<- sapply(clusters,function(x){	max(x$size.tips)  	})
	clusters.info		<- as.data.table( cbind(thresh,clusters.totalclu,clusters.pairs,clusters.maxsize,clusters.seqinclu,clusters.NLseqinclu) )
	clusters.info		<- subset(clusters.info, brl>0.025 & bs>0.6)	
	clusters.info.bs.sd	<- clusters.info[, list(sd.totalclu=sd(clusters.totalclu), sd.pairs=sd(clusters.pairs), sd.maxsize=sd(clusters.maxsize), sd.seqinclu=sd(clusters.seqinclu), sd.NLseqinclu=sd(clusters.seqinclu)), by=bs]
	clusters.info.brl.sd<- clusters.info[, list(sd.totalclu=sd(clusters.totalclu), sd.pairs=sd(clusters.pairs), sd.maxsize=sd(clusters.maxsize), sd.seqinclu=sd(clusters.seqinclu), sd.NLseqinclu=sd(clusters.seqinclu)), by=brl]
	print(clusters.info)
	print(clusters.info.bs.sd)
	print(clusters.info.brl.sd)	
}		
######################################################################################
project.hivc.beast<- function(dir.name= DATA)
{
	stop()
	require(ape)
	require(data.table)
	require(RColorBrewer)
	require(XML)
	if(1)
	{
		indir					<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/beast/beast2_140201'
		infile					<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_seroneg-130"
		insignat				<- "Tue_Aug_26_09:13:47_2013"
		infilexml.opt			<- "rsu815"
		infilexml.template		<- "sasky_sdr06"
		
		argv		<<- hivc.cmd.beast2.getclustertrees(indir, infile, insignat, infilexml.opt, infilexml.template, burnin=5e6, outdir=indir, outsignat=insignat, prog= PR.BEAST2CLUTREES, verbose=1, resume=1)
		argv		<<- unlist(strsplit(argv,' '))
		hivc.prog.BEAST2.get.cluster.trees()		

		argv		<<- hivc.cmd.beast2.processclustertrees(indir, infile, insignat, infilexml.opt, infilexml.template, outdir=indir, outsignat=insignat, prog= PR.BEAST2CLUPOSTERIOR, verbose=1, resume=0)
		argv		<<- unlist(strsplit(argv,' '))
		hivc.prog.BEAST2.process.cluster.trees()

		indircov	<- paste(DATA,"derived",sep='/')	
		infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"					
		argv		<<- hivc.cmd.beast2.plotclustertrees(indir, infile, insignat, indircov, infilecov, infilexml.opt, infilexml.template, outdir=indir, outsignat=insignat, prog=PR.BEAST2.PLOTCLUTREES, resume=1, verbose=1)
		argv		<<- unlist(strsplit(argv,' '))
		hivc.prog.BEAST2.plot.cluster.trees()
	}
	if(1)	#select single likely transmitter for every recent as in the Brighton study Fisher et al
	{
		project.athena.Fisheretal()		
	}
	if(0)	#get BEAST nexus file for seroconverters
	{		
		indir			<- paste(DATA,"tmp",sep='/')		
		infile			<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		insignat		<- "Thu_Aug_01_17/05/23_2013"
		indircov		<- paste(DATA,"derived",sep='/')
		infilecov		<- "ATHENA_2013_03_AllSeqPatientCovariates"
		infiletree		<- paste(infile,"examlbs100",sep="_")
		infilexml		<- paste(infile,'_',"beast",'_',"seroneg",sep='')
		
		outdir			<- indir
		outsignat		<- "Tue_Aug_26_09/13/47_2013"
		
		opt.brl			<- "dist.brl.casc" 
		thresh.brl		<- 0.096
		thresh.bs		<- 0.8
		pool.ntip		<- 130
		infilexml.opt	<- "txs4clu"
		resume			<- 1
		verbose			<- 1
		
		argv		<<- hivc.cmd.beast.poolrunxml(indir, infile, insignat, indircov, infilecov, infiletree, infilexml, outsignat, pool.ntip, infilexml.opt=infilexml.opt, opt.brl=opt.brl, thresh.brl=thresh.brl, thresh.bs=thresh.bs, resume=resume, verbose=1)
		cat(argv)
		argv		<<- unlist(strsplit(argv,' '))
		hivc.prog.BEAST.generate.xml()		
	}
	if(0)	#get distribution of TMRCAs of tip nodes for a particular BEAST run
	{
		indir				<- paste(DATA,"beast/beast_131011",sep='/')		
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_beast_seroneg"
		insignat			<- "Tue_Aug_26_09/13/47_2013"
		infilexml.opt		<- "mph4clutx4tipLdTd"
		infilexml.template	<- "um232rhU2045"
		verbose				<- 1
		
		
		tmp			<- list.files(indir, pattern=paste(".log$",sep=''))
		files		<- tmp[ grepl(paste('_',infilexml.template,'_',sep=''), tmp) & grepl(paste('_',infilexml.opt,'_',sep=''), tmp) & grepl(gsub('/',':',insignat), tmp) ]
		if(verbose)	cat(paste("\nFound files matching input args, n=", length(files)))
		
		#	read length of tip stems
		df.tstem	<- lapply( seq_along(files), function(i)
					{
						file.log	<- files[i]
						if(verbose)	cat(paste("\nprocess file,", file.log))
						file.log	<- paste(indir,file.log,sep='/')
						file.xml	<- paste(substr(file.log,1,nchar(file.log)-3), "xml", sep='')
						if(verbose)	cat(paste("\nand file,", file.xml))
						hivc.beast.read.log2tstem(file.log, file.xml, beastlabel.idx.samplecode=6, burn.in= 5e6, breaks.n= 30, verbose=1)					
					})
		df.tstem<- rbindlist( df.tstem )		
	}
	if(1)	#plot BEAST clusters
	{					
		indir				<- paste(DATA,"beast/beast_130908",sep='/')		
		infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
		insignat			<- "Tue_Aug_26_09/13/47_2013"
		indircov			<- paste(DATA,"derived",sep='/')
		infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
		infiletree			<- paste(infile,"examlbs100",sep="_")
		infilexml			<- paste(infile,'_',"beast",'_',"seroneg",sep='')
		#infilexml.template	<- "um22rhU2050"
		infilexml.template	<- "um22rhG202018"
		#infilexml.template	<- "rhU65rho753"
		#infilexml.template	<- "rhU65rho903"
		#infilexml.template	<- "rhU65rho906"
		#infilexml.template	<- "rhU65rho909"	
		#infilexml.template	<- "um181rhU2045"
		#infilexml.template	<- "um182rhU2045"
		#infilexml.template	<- "um183rhU2045"
		#infilexml.template	<- "um182us45"
		#infilexml.template	<- "um182us60"
		infilexml.opt		<- "txs4clu"
		#infilexml.opt		<- "txs4clufx03"
		#infilexml.opt		<- "mph4clu"
		#infilexml.opt		<- "mph4clumph4tu"
		#infilexml.opt		<- "mph4clufx03"
	
		argv				<<- hivc.cmd.beast.evalrun(indir, infilexml, insignat, infilexml.opt, infilexml.template, pool.n, verbose=verbose)
		argv				<<- unlist(strsplit(argv,' '))
		hivc.prog.BEAST.evalpoolrun()
	}
}
######################################################################################
project.ukca.demographics<- function()
{
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/tblBASIC.csv'
	df.b	<- read.csv(f, stringsAsFactors=FALSE)
	df.b	<- as.data.table(df.b)
	set(df.b, NULL, 'NUCSEQ_D', hivc.db.Date2numeric(df.b[, as.Date(NUCSEQ_D, "%d/%m/%Y")]))
	set(df.b, NULL, 'HIVPOS_D', hivc.db.Date2numeric(df.b[, as.Date(HIVPOS_D, "%d/%m/%Y")]))
	set(df.b, NULL, 'BIRTH_D', hivc.db.Date2numeric(df.b[, as.Date(BIRTH_D, "%d/%m/%Y")]))
	set(df.b, NULL, 'DEATH_D', hivc.db.Date2numeric(df.b[, as.Date(DEATH_D, "%d/%m/%Y")]))
	set(df.b, NULL, 'LTALIV_D', hivc.db.Date2numeric(df.b[, as.Date(LTALIV_D, "%d/%m/%Y")]))	
	set(df.b, NULL, 'SEROCO_D', hivc.db.Date2numeric(df.b[, as.Date(SEROCO_D, "%d/%m/%Y")]))
	set(df.b, NULL, 'MODE', df.b[, factor(MODE, levels=c(1,2,3,4,5,6,7,8,90,99), labels=c('MSM','IDU','MSM+IDU','HAEM','5','HET','HET+IDU','8','OTH','NA'))])
	
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/tblLAB_RES.csv'
	df.s	<- read.csv(f, stringsAsFactors=FALSE)
	df.s	<- as.data.table(df.s)
	set(df.s, NULL, 'NUCSEQ_D', hivc.db.Date2numeric(df.s[, as.Date(NUCSEQ_D, "%Y-%m-%d")]))
	set(df.s, NULL, 'SEQ_DT', hivc.db.Date2numeric(df.s[, as.Date(SEQ_DT, "%Y-%m-%d")]))
	
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/tblLAB_CD4.csv'
	df.c	<- read.csv(f, stringsAsFactors=FALSE)
	df.c	<- as.data.table(df.c)
	set(df.c, NULL, 'CD4_D', hivc.db.Date2numeric(df.c[, as.Date(CD4_D, "%Y-%m-%d")]))
	
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/tblLAB_RNA.csv'
	df.v	<- read.csv(f, stringsAsFactors=FALSE)
	df.v	<- as.data.table(df.v)
	set(df.v, NULL, 'RNA_D', hivc.db.Date2numeric(df.v[, as.Date(RNA_D, "%Y-%m-%d")]))
	set(df.v, NULL, 'RNA_Vl10', df.v[, log10(RNA_V)])
	
	#individuals with sequence
	ggplot(df.s, aes(x=NUCSEQ_D)) + geom_histogram(binwidth=1) +
			scale_x_continuous(breaks=seq(1985,2020,1)) + scale_y_continuous(breaks=seq(0,1600,200))
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/expl_seqtotal.pdf'
	ggsave(f, w=14,h=6)	
	
	#seroconverters with sequence
	df.a	<- merge(df.b, subset(df.s, select=c(PATIENT, NUCSEQ_ID)), by='PATIENT')
	setkey(df.a, PATIENT)
	df.a<- unique(df.a)	
	df.a[, table(MODE)]
	#MODE
    #	MSM     IDU 	MSM+IDU    HAEM       5     HET HET+IDU       8     OTH      NA 
   	#	6638    1313      90      46      75    3222     227       6     173     387 
	ggplot(subset(df.a, !is.na(SEROCO_D)), aes(x=SEROCO_D)) + geom_histogram(binwidth=1) +
			scale_x_continuous(breaks=seq(1980,2020,1)) + scale_y_continuous(breaks=seq(0,500,50))
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/expl_seroconvtotal.pdf'
	ggsave(f, w=14,h=6)
	#	last alive if not dead
	tmp		<- subset(df.a, DEATH_Y==0)	
	tmp[, NOTSEEN:=max(LTALIV_D, na.rm=1)-LTALIV_D]
	ggplot(tmp, aes(x=NOTSEEN)) + geom_histogram(binwidth=0.5)
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/expl_ltalive.pdf'
	ggsave(f, w=8,h=6)
	
	#seroconverters with CD4 within 1 year of HIV+
	tmp		<- merge(df.a, df.c,by='PATIENT')
	setkey(tmp, PATIENT, CD4_D)
	tmp		<- unique(tmp)
	tmp		<- subset(tmp, abs(CD4_D-HIVPOS_D)<1 )
	setkey(tmp, PATIENT)
	tmp		<- unique(tmp)
	#	3505
	#ggplot(tmp, aes(x=CD4_V)) + geom_histogram(binwidth=200) + scale_x_continuous(breaks=seq(0,4000,200)) 
	
	#seroconverters with VL
	tmp		<- merge(df.a, df.v,by='PATIENT')
	setkey(tmp, PATIENT, RNA_D)
	tmp		<- unique(tmp)
	tmp[, list(n=length(RNA_Vl10)), by='PATIENT'][, median(n)]
	subset(tmp, !is.na(SEROCO_D))[, length(unique(PATIENT))]
	#	3870		
	subset(tmp, !is.na(SEROCO_D) & !is.nan(RNA_Vl10) & RNA_D-HIVPOS_D>0.5 & RNA_D-HIVPOS_D<2)[, length(unique(PATIENT))]
	#	609
	ggplot(subset(tmp, !is.na(SEROCO_D) & !is.nan(RNA_Vl10) & RNA_D-HIVPOS_D>0.5 & RNA_D-HIVPOS_D<2), aes(x=RNA_Vl10)) + 
			scale_x_continuous(breaks=seq(0,8,1)) +
			geom_histogram(binwidth=0.25)
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/expl_seroconvSPVL.pdf'
	ggsave(f, w=8,h=8)
	
	

}
######################################################################################
project.ukca.TPTN.bootstrapvalues<- function(dir.name= DATA)
{
	require(ape)
	require(data.table)
	require(RColorBrewer)
	
	indir				<- paste(DATA,"tmp",sep='/')
	infile				<- "UKCA_2013_07_TNTPHIVnGTR_examlbs100"
	insignat			<- "Mon_Sep_22_17/23/46_2013"		
	file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".newick",sep='')
	plot.file			<- paste(indir,'/',infile,"_preclust_distbs_",gsub('/',':',insignat),".pdf", sep='')
	ph					<- read.tree(file)
	ph$node.label		<- as.numeric(ph$node.label)/100
	ph$node.label[1]	<- 0
	ph$node.label[3]	<- 0.01		#little phylo cleaning hack ;-)
	ph.mrca				<- mrca(ph)
	dist.brl.casc		<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
	
	infile				<- "UKCA_2013_07_B_true_pos"
	file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".csv",sep='')
	df.tp				<- as.data.table(read.csv(file, header=F, sep=' ', col.names=c("tip1","tip2")))		
	set(df.tp,NULL,"tip1", as.character(df.tp[,tip1]))
	set(df.tp,NULL,"tip2", as.character(df.tp[,tip2]))
	if(verbose)	cat(paste("\nNumber of TP pairs read, n=",nrow(df.tp)))
	df.tp				<- subset(df.tp, df.tp[,tip1]%in%rownames(ph.mrca)  &  df.tp[,tip2]%in%rownames(ph.mrca))
	if(verbose)	cat(paste("\nNumber of TP pairs with tips in ph, n=",nrow(df.tp)))
	
	infile				<- "UKCA_2013_07_B_true_neg"
	file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".csv",sep='')
	df.tn				<- as.data.table(read.csv(file, header=F, sep=' ', col.names=c("tip1","tip2")))		
	set(df.tn,NULL,"tip1", as.character(df.tn[,tip1]))
	set(df.tn,NULL,"tip2", as.character(df.tn[,tip2]))
	if(verbose)	cat(paste("\nNumber of TN pairs read, n=",nrow(df.tn)))		
	df.tn				<- subset(df.tn, df.tn[,tip1]%in%rownames(ph.mrca)  &  df.tn[,tip2]%in%rownames(ph.mrca))
	if(verbose)	cat(paste("\nNumber of TN pairs with tips in ph, n=",nrow(df.tn)))
	
	#set MRCAs
	df.tn[,dummy:=seq_len(nrow(df.tn))]
	df.tp[,dummy:=seq_len(nrow(df.tp))]
	df.tn				<- df.tn[, list(tip1=tip1, tip2=tip2, mrca= ph.mrca[tip1,tip2]),by="dummy"]
	df.tp				<- df.tp[, list(tip1=tip1, tip2=tip2, mrca= ph.mrca[tip1,tip2]),by="dummy"]
	
	bs.linked.bypatient	<- df.tp
	bs.unlinkedpairs	<- df.tn		
	
	hivc.phy.get.TP.and.TN.bootstrapvalues(ph, bs.linked.bypatient, ph.mrca=ph.mrca ,df.seqinfo=NULL, bs.unlinkedpairs=bs.unlinkedpairs, bs.unlinked.byspace=NULL, dist.brl=NULL, thresh.brl=0.096, plot.file=plot.file, verbose= 1)	
}

######################################################################################
project.athena.Niaetal.similar<- function()
{
	require(data.table)
	require(ape)
	indir						<- paste(DATA,"tmp",sep='/')		
	indircov					<- paste(DATA,"derived",sep='/')
	outdir						<- paste(DATA,"NiaMSC",sep='/')	
	infile.cov					<- paste(indircov,"ATHENA_2013_03_AllSeqPatientCovariates.R",sep='/')
	infile.viro					<- paste(indircov,"ATHENA_2013_03_Viro.R",sep='/')
	infile.immu					<- paste(indircov,"ATHENA_2013_03_Immu.R",sep='/')
	infile.treatment			<- paste(indircov,"ATHENA_2013_03_Regimens.R",sep='/')
	infile.seq					<- paste(indir,'/',"ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_",gsub('/',':',"Wed_Dec_18_11:37:00_2013"),".R",sep='')
	#
	#	load patient RNA
	#
	load(infile.viro)
	set(df, NULL, 'Patient', df[, as.character(Patient)])
	df.viro				<- df
	#
	#	load patient CD4
	#
	load(infile.immu)
	set(df, NULL, 'Patient', df[, as.character(Patient)])
	df.immu				<- df
	#
	#	load patient regimen
	#
	load(infile.treatment)
	set(df, NULL, 'Patient', df[, as.character(Patient)])
	df.treatment		<- df		
	#
	#	load all covariates
	#
	load(infile.cov)
	set(df.all, NULL, 'Patient', df.all[, as.character(Patient)])
	df.demo				<- df.all		
	#
	#	load sequences
	#
	load(infile.seq)	
	#
	#	anonymize
	#
	tmp			<- unique(subset( df.demo, select=Patient))
	tmp[, Patient.new:= paste('P',seq_len(nrow(tmp)),sep='')]
	df.demo	<- merge(df.demo, tmp, by='Patient')
	df.demo	<- subset( df.demo, select=c(Patient, Patient.new, FASTASampleCode, PosSeqT, DateBorn, Sex, CountryBorn, RegionOrigin, DateDied, Subtype, isAcute, NegT,
						NegT_Acc, PosT, PosT_Acc, CountryInfection, Trm, DateLastContact, RegionHospital, DateFirstEverCDCC, isDead, AnyPos_T1, PoslRNA_T1, lRNA_T1, PosCD4_T1, CD4_T1, AnyT_T1 ))
	
	tmp			<- unique(subset( df.demo, select=c(Patient.new, FASTASampleCode, PosSeqT)))
	setkey(tmp, Patient.new, PosSeqT)
	tmp			<- cbind( tmp, tmp[, list(FASTASampleCode.new= paste('S',substr(Patient.new,2,nchar(Patient.new)),'-',seq_along(FASTASampleCode), sep='')),by=Patient.new])
	df.demo		<- merge(df.demo, subset(tmp, select=c(FASTASampleCode, FASTASampleCode.new)), by='FASTASampleCode')
	
	setnames( df.demo, c('Patient', 'Patient.new', 'FASTASampleCode', 'FASTASampleCode.new'), c('Patient.old', 'Patient', 'FASTASampleCode.old', 'FASTASampleCode') )	
	setnames( df.treatment, 'Patient', 'Patient.old')
	setnames( df.viro, 'Patient', 'Patient.old')
	setnames( df.immu, 'Patient', 'Patient.old')
	df.treatment	<- merge(unique(subset(df.demo, select=c(Patient, Patient.old))), df.treatment, by='Patient.old')
	df.viro			<- merge(unique(subset(df.demo, select=c(Patient, Patient.old))), df.viro, by='Patient.old')
	df.immu			<- merge(unique(subset(df.demo, select=c(Patient, Patient.old))), df.immu, by='Patient.old')
	
	tmp				<- subset(df.demo, select=c(FASTASampleCode, FASTASampleCode.old) )
	setkey(tmp, FASTASampleCode.old)
	rownames(seq.PROT.RT)	<- tmp[rownames(seq.PROT.RT), ][, FASTASampleCode]
	#	save with keys
	file			<- paste(outdir, 'data_with_keys.R',sep='/')
	save(df.demo, df.treatment, df.viro, df.immu, seq.PROT.RT, file=file)
	
	df.treatment[, Patient.old:=NULL]
	df.viro[, Patient.old:=NULL]
	df.immu[, Patient.old:=NULL]
	df.demo[, Patient.old:=NULL]
	df.demo[, FASTASampleCode.old:=NULL]
	#	save without keys
	file			<- paste(outdir, 'data_without_keys.R',sep='/')
	save(df.demo, df.treatment, df.viro, df.immu, seq.PROT.RT, file=file)

}
######################################################################################
project.athena.Xavieretal.similar<- function()
{
	require(data.table)
	require(ape)
	stop()
	indir					<- paste(DATA,"tmp",sep='/')		
	indircov				<- paste(DATA,"derived",sep='/')
	outdir					<- paste(DATA,"fisheretal",sep='/')
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	
	if(0)
	{
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"				
		clu.indir				<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/sasky_sdr06_-DR-RC-SH+LANL_alrh160'
		clu.infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		clu.insignat			<- "Wed_Dec_18_11:37:00_2013"
		clu.infilexml.opt		<- "alrh160"
		clu.infilexml.template	<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=2_sasky',sep='_')
	}
	if(1)
	{
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"				
		clu.indir				<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/gmrf_sdr06fr_-DR-RC-SH+LANL_um192rhU2080'
		clu.infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		clu.insignat			<- "Wed_Dec_18_11:37:00_2013"
		clu.infilexml.opt		<- "mph4clutx4tip"
		clu.infilexml.template	<- "um192rhU2080"	
		outfile					<- paste(infile,'gmrf',sep='_')
		outdir					<- paste(DATA,'Xavieretal',sep='/')
	}
	#	fixed input args
	opt.brl			<- "dist.brl.casc" 
	thresh.brl		<- 0.096
	thresh.bs		<- 0.8	
	#
	# 	recent msm 
	#
	load(paste(indircov,'/',infilecov,'.R',sep=''))	
	#
	#	load clustering msm
	#
	argv			<<- hivc.cmd.clustering.msm(indir, infiletree, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=1)
	argv			<<- unlist(strsplit(argv,' '))		
	msm				<- hivc.prog.get.clustering.MSM()
	clumsm.info		<- msm$df.cluinfo
	#	update clumsm.info with latest values from df.all
	clumsm.info		<- merge( df.all, subset(clumsm.info, select=c(FASTASampleCode, cluster, Node, clu.npat, clu.ntip, clu.nFrgnInfection, clu.fPossAcute, clu.AnyPos_T1, clu.bwpat.medbrl)), by='FASTASampleCode')
	#
	set(clumsm.info, NULL, 'PosSeqT', hivc.db.Date2numeric(clumsm.info[,PosSeqT]))
	set(clumsm.info, NULL, 'DateBorn', hivc.db.Date2numeric(clumsm.info[,DateBorn]))
	set(clumsm.info, NULL, 'DateDied', hivc.db.Date2numeric(clumsm.info[,DateDied]))
	set(clumsm.info, NULL, 'NegT', hivc.db.Date2numeric(clumsm.info[,NegT]))
	set(clumsm.info, NULL, 'PosT', hivc.db.Date2numeric(clumsm.info[,PosT]))
	set(clumsm.info, NULL, 'DateLastContact', hivc.db.Date2numeric(clumsm.info[,DateLastContact]))
	set(clumsm.info, NULL, 'DateFirstEverCDCC', hivc.db.Date2numeric(clumsm.info[,DateFirstEverCDCC]))
	set(clumsm.info, NULL, 'AnyPos_T1', hivc.db.Date2numeric(clumsm.info[,AnyPos_T1]))
	set(clumsm.info, NULL, 'PoslRNA_T1', hivc.db.Date2numeric(clumsm.info[,PoslRNA_T1]))
	set(clumsm.info, NULL, 'PoslRNA_TS', hivc.db.Date2numeric(clumsm.info[,PoslRNA_TS]))
	set(clumsm.info, NULL, 'lRNA.hb4tr_LT', hivc.db.Date2numeric(clumsm.info[,lRNA.hb4tr_LT]))
	set(clumsm.info, NULL, 'PosCD4_T1', hivc.db.Date2numeric(clumsm.info[,PosCD4_T1]))
	set(clumsm.info, NULL, 'PosCD4_TS', hivc.db.Date2numeric(clumsm.info[,PosCD4_TS]))
	set(clumsm.info, NULL, 'AnyT_T1', hivc.db.Date2numeric(clumsm.info[,AnyT_T1]))
	set(clumsm.info, NULL, 'clu.AnyPos_T1', hivc.db.Date2numeric(clumsm.info[,clu.AnyPos_T1]))
	#	make selection of pilot clusters
	cluster.pilot	<- unique( subset(clumsm.info, Patient%in%c('M32179','M28593','M37538','M15184','M12774'), cluster) )	
	clumsm.info		<- merge(clumsm.info, cluster.pilot, by='cluster')
	clumsm.info		<- subset(clumsm.info, select=c(cluster, FASTASampleCode, Patient,  DateBorn, NegT, AnyPos_T1, PosSeqT, AnyT_T1, DateLastContact, DateDied, isAcute))
	setkey(clumsm.info, cluster, Patient, FASTASampleCode)
	#
	#	
	tmp						<- unique(subset(clumsm.info, select=Patient))	
	file.viro				<- paste(indircov,"/ATHENA_2013_03_Viro.R",sep='/')
	file.immu				<- paste(indircov,"/ATHENA_2013_03_Immu.R",sep='/')
	file.treatment			<- paste(indircov,"/ATHENA_2013_03_Regimens.R",sep='/')		
	#	load patient RNA
	load(file.viro)
	clumsm.viro				<- merge(tmp, subset(df, select=c(Patient, PosRNA, lRNA)), by='Patient')
	set(clumsm.viro, NULL, 'PosRNA', hivc.db.Date2numeric(clumsm.viro[,PosRNA]))
	
	#	load patient CD4				
	load(file.immu)
	clumsm.cd4				<- merge(tmp, subset(df, select=c(Patient, PosCD4, CD4)), by='Patient')
	set(clumsm.cd4, NULL, 'PosCD4', hivc.db.Date2numeric(clumsm.cd4[,PosCD4]))
	set(clumsm.cd4, NULL, 'CD4', clumsm.cd4[,round(CD4)])
	
	#	load patient regimen
	load(file.treatment)
	clumsm.art				<- merge(tmp, subset(df, select=c(Patient, AnyT_T1, StartTime, StopTime, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel)), by='Patient')
	set(clumsm.art, NULL, 'StartTime', hivc.db.Date2numeric(clumsm.art[,StartTime]))
	set(clumsm.art, NULL, 'StopTime', hivc.db.Date2numeric(clumsm.art[,StopTime]))
	set(clumsm.art, NULL, 'AnyT_T1', hivc.db.Date2numeric(clumsm.art[,AnyT_T1]))
	setnames(clumsm.art, c('StartTime','StopTime','TrI','TrCh.failure','TrCh.adherence','TrCh.patrel'), c('ART.TS','ART.TE','ART.Intrpt','ART.Fail','ART.Adh','ART.PatDec'))

	#	anonymize	
	tmp			<- unique(subset( clumsm.info, select=Patient))
	tmp[, Patient.new:= paste('P',seq_len(nrow(tmp)),sep='')]
	clumsm.info	<- merge(clumsm.info, tmp, by='Patient')
	tmp			<- unique(subset( clumsm.info, select=c(Patient.new, FASTASampleCode, PosSeqT)))
	setkey(tmp, Patient.new, PosSeqT)
	tmp			<- cbind( tmp, tmp[, list(FASTASampleCode.new= paste('S',substr(Patient.new,2,nchar(Patient.new)),'-',seq_along(FASTASampleCode), sep='')),by=Patient.new]) 
	clumsm.info	<- merge(clumsm.info, subset(tmp, select=c(FASTASampleCode, FASTASampleCode.new)), by='FASTASampleCode')
	setnames( clumsm.info, c('Patient', 'Patient.new', 'FASTASampleCode', 'FASTASampleCode.new'), c('Patient.old', 'Patient', 'FASTASampleCode.old', 'FASTASampleCode') )	
	setnames( clumsm.art, 'Patient', 'Patient.old')
	setnames( clumsm.viro, 'Patient', 'Patient.old')
	setnames( clumsm.cd4, 'Patient', 'Patient.old')
	clumsm.art	<- merge(unique(subset(clumsm.info, select=c(Patient, Patient.old))), clumsm.art, by='Patient.old')
	clumsm.viro	<- merge(unique(subset(clumsm.info, select=c(Patient, Patient.old))), clumsm.viro, by='Patient.old')
	clumsm.cd4	<- merge(unique(subset(clumsm.info, select=c(Patient, Patient.old))), clumsm.cd4, by='Patient.old')
	clumsm.art[, Patient.old:=NULL]
	clumsm.viro[, Patient.old:=NULL]
	clumsm.cd4[, Patient.old:=NULL]
	
	#	load MAP cluster topology with branch lengths
	files		<- list.files(clu.indir)
	if(!length(files))	stop('no input files matching criteria')
	files		<- files[ sapply(files, function(x) grepl(clu.infile, x, fixed=1) & grepl(gsub('/',':',clu.insignat), x, fixed=1) & grepl(clu.infilexml.opt, x, fixed=1) & grepl(clu.infilexml.template,x, fixed=1) & grepl('_cluposterior_[0-9]+',x) & grepl('R$',x) ) ]		
	if(!length(files))	stop('no input files matching criteria')
	tmp			<- regmatches( files, regexpr('_cluposterior_[0-9]+',files)) 
	cluster		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
	file.info	<- data.table(file=files, cluster=cluster)
	setkey(file.info, cluster)
	#	subset of clusters: those in df.tpairs.select	
	file.info	<- merge(file.info, unique(subset(clumsm.info, select=cluster)), by='cluster')
	clumsm.info	<- merge(clumsm.info, subset(file.info, select=cluster), by='cluster')				
	#	combine dated cluster phylogenies
	clumsm.ph	<- lapply(seq_len(nrow(file.info)), function(i)
		{				
			#	load dated cluster phylogenies
			file					<- paste(clu.indir, file.info[i,file], sep='/')
			cat(paste('\nload file=',file,'i=',i))
			tmp						<- load(file)
			topo.map				<- mph.clu.dtopo[which.max(freq),]					
			tmp						<- which( grepl( paste('mph.i=',topo.map[,mph.i],'_',sep=''), names(ph.consensus) ) )
			if(length(tmp)!=1)	stop('unexpected ph.consensus index')
			topo.map.ph				<- ph.consensus[[ tmp ]]
			topo.map.SA				<- subset(mph.SA.cnt, equal.to==topo.map[,mph.i])
			set(topo.map.SA, NULL, 'SA.freq', topo.map.SA[,SA.freq/n])													
			topo.map.nodectime		<- subset(mph.node.ctime, equal.to==topo.map[,mph.i])
			topo.map				<- merge(subset(topo.map, select=c(cluster, dens)), subset(topo.map.SA, select=c(cluster, tip, SA.freq)), by='cluster')
			topo.map.nodectime		<- subset(topo.map.nodectime, select=c(cluster, node, q, cdf, pdf))
			#
			tmp						<- merge( data.table( BEASTlabel=topo.map.ph$tip.label, FASTASampleCode.old= sapply( strsplit(topo.map.ph$tip.label, '_'), '[[', 6 ) ), subset(clumsm.info, select=c(FASTASampleCode, FASTASampleCode.old)), by='FASTASampleCode.old')
			setkey(tmp, BEASTlabel)
			topo.map.ph$tip.label	<- tmp[topo.map.ph$tip.label, ][,FASTASampleCode]
			topo.map.ph$node.dated	<- topo.map.nodectime
			topo.map.ph			
		})

	#	save KEY used for anonymization
	clumsm.info
	file		<- paste(outdir,'/',outfile,'_key.R',sep='')
	save(clumsm.info, file=file)
	clumsm.info[, Patient.old:=NULL]
	clumsm.info[, FASTASampleCode.old:=NULL]
	clumsm.info	<- subset(clumsm.info, select=c(cluster, Patient, FASTASampleCode, DateBorn, NegT, AnyPos_T1,  PosSeqT,  AnyT_T1, DateLastContact, DateDied, isAcute ))
	
	#	save info per cluster
	setkey(clumsm.info, FASTASampleCode)	
	setkey(clumsm.art, Patient)
	setkey(clumsm.cd4, Patient)
	setkey(clumsm.viro, Patient)
	tmp			<- lapply(seq_len(nrow(file.info)), function(i)
		{
			clu			<- list()
			clu$ph		<- clumsm.ph[[i]]
			clu$info	<- clumsm.info[ clu$ph$tip.label, ]
			clu$art		<- clumsm.art[ clu$info[, unique(Patient)], ]
			clu$viro	<- clumsm.viro[ clu$info[, unique(Patient)], ]
			clu$cd4		<- clumsm.cd4[ clu$info[, unique(Patient)], ]
			
			file		<- paste(outdir,'/',outfile,'_clu',clu$info[,unique(cluster)],'_raw.rda',sep='')
			cat(paste('\nsave to file', file))
			save(file=file, clu)			
		})
}
######################################################################################
project.athena.TPTN.bootstrapvalues<- function(dir.name= DATA)
{	
	require(ape)
	require(data.table)
	require(RColorBrewer)
	
	verbose		<- 1
	resume		<- 1
	
	# load sequences
	indir								<- paste(dir.name,"tmp",sep='/')
	infile								<- "ATHENA_2013_03_CurAll+LANL_Sequences"			
	insignat							<- "Sat_Jun_16_17/23/46_2013"
	file								<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')			
	load(file)		
	#
	# precompute clustering stuff		
	#
	patient.n	<- 15700
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
	insignat	<- "Sat_Jun_16_17/23/46_2013"
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	indircov	<- paste(DATA,"derived",sep='/')
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"							
	argv		<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv		<<- unlist(strsplit(argv,' '))
	clu.pre		<- hivc.prog.get.clustering.precompute()
	
	load(paste(indir,'/',"mrcas.R",sep=''))
	
	outfile				<- paste(infile,"preclust",sep='_')
	
	
	dist.brl			<- clu.pre$dist.brl.casc
	linked.bypatient	<- clu.pre$linked.bypatient
	unlinked			<- clu.pre$ph.unlinked
	unlinked.byspace	<- clu.pre$unlinked.byspace
	unlinked.bytime		<- clu.pre$unlinked.bytime
	unlinked.info		<- clu.pre$ph.unlinked.info
	ph					<- clu.pre$ph
	df.seqinfo			<- clu.pre$df.seqinfo
	#ph.mrca			<- mrca(ph)		#probably fastest
	thresh.brl			<- 0.096
	plot.file			<- paste(indir,'/',outfile,"_distbs_",gsub('/',':',insignat),".pdf", sep='')
	
	#prepare data.tables with mrca		
	bs.unlinkedpairs	<- lapply(unlinked.bytime, function(x)
			{						
				set(x, NULL, "query.FASTASampleCode", as.character(x[,query.FASTASampleCode]))
				set(x, NULL, "FASTASampleCode", as.character(x[,FASTASampleCode]))
				ans					<- merge(x, x[, list(mrca= ph.mrca[query.FASTASampleCode,FASTASampleCode]) ,by=FASTASampleCode],by="FASTASampleCode")
				setnames(ans, c("FASTASampleCode","query.FASTASampleCode"), c("tip2","tip1"))
				subset(ans, select=c(tip1, tip2, mrca))
			})
	bs.unlinkedpairs	<- rbindlist(bs.unlinkedpairs)
	#
	unlinked.byspace[,dummy:=seq_len(nrow(unlinked.byspace))]
	set(unlinked.byspace, NULL, "FASTASampleCode", as.character(unlinked.byspace[,FASTASampleCode]))	
	seq.indb			<- colnames(ph.mrca)[ which( substr(colnames(ph.mrca),1,2)!="TN" & substr(colnames(ph.mrca),1,8)!="PROT+P51" ) ]
	bs.unlinked.byspace	<- unlinked.byspace[,	list(tip1=FASTASampleCode,tip2=seq.indb, mrca= ph.mrca[seq.indb,FASTASampleCode]), by="dummy"]
	#
	setkey(linked.bypatient,Patient)
	bs.linked.bypatient	<- linked.bypatient[, {
				tmp					<- match(FASTASampleCode, ph$tip.label)
				ans					<- t(combn(tmp, 2 ))
				ans					<- cbind(ans, apply(ans, 1, function(z)  ph.mrca[z[1],z[2]]))
				data.table(tip1=ans[,1], tip2=ans[,2], mrca=ans[,3])
			}, by="Patient"]
	#	get bootstrap values
	ph.mrca				<- mrca(ph)
	tmp					<- hivc.phy.get.TP.and.TN.bootstrapvalues(ph,  bs.linked.bypatient, ph.mrca=ph.mrca, clu.pre$df.seqinfo, bs.unlinkedpairs=bs.unlinkedpairs, bs.unlinked.byspace=bs.unlinked.byspace, dist.brl=clu.pre$dist.brl.casc, thresh.brl=0.096, plot.file=paste(indir,'/',outfile,"_distbs_",gsub('/',':',insignat),".pdf"), verbose= 1)	
	#	further analysis of those with pairs with PosSeq.diff==0		
	if(0)
	{
		bs.linked.bypatient.eqPosSeqT		<- subset(bs.linked.bypatient, PosSeqT.diff==0)
		# compute raw genetic distance between sequences
		indir								<- paste(dir.name,"tmp",sep='/')
		infile								<- "ATHENA_2013_03_CurAll+LANL_Sequences"			
		insignat							<- "Sat_Jun_16_17/23/46_2013"
		file								<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')			
		load(file)
		#dist.dna( rbind( seq.PROT.RT[bs.linked.bypatient.eqPosSeqT[i,tip1], ], seq.PROT.RT[bs.linked.bypatient.eqPosSeqT[i,tip2], ]), model="raw", as.matrix=1)
		#tmp		<- seq.dist( rbind( seq.PROT.RT[bs.linked.bypatient.eqPosSeqT[i,tip1], ], seq.PROT.RT[bs.linked.bypatient.eqPosSeqT[i,tip2], ]) )
		dummy	<- 0
		tmp		<- sapply(seq_len(nrow(bs.linked.bypatient.eqPosSeqT)),function(i)
				{
					1 - .C("hivc_dist_ambiguous_dna", seq.PROT.RT[bs.linked.bypatient.eqPosSeqT[i,tip1], ], seq.PROT.RT[bs.linked.bypatient.eqPosSeqT[i,tip2], ], ncol(seq.PROT.RT), dummy )[[4]]
				})
		bs.linked.bypatient.eqPosSeqT[,dist:=tmp]
		#
		#	conclusion: despite identical sequences, BS can be suprisingly low! Unclear if this improves with more BS replications
		#
		save(bs.linked.bypatient.eqPosSeqT, file="/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/notes/20130917_ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100_poorBSvalues.R")
	}
}
######################################################################################
project.hivc.clustering<- function(dir.name= DATA)
{
	require(ape)
	require(data.table)
	require(RColorBrewer)
	if(1)
	{
		require(adephylo)
		indir			<- paste(DATA, 'tmp', sep='/')
		infiletree		<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500"
		insignat		<- "Wed_Dec_18_11:37:00_2013"								
		file			<- paste(indir, '/', infiletree, '_', gsub('/',':',insignat),".R",sep='')
		cat(paste('\nget dist tips for file=',file))
		load(file)	#loads ph		
		brl				<- distTips(ph , method='patristic')		
		file			<- paste(indir, '/', infiletree, '_', gsub('/',':',insignat),"_distTips.R",sep='')
		cat(paste('\nsave dist tips to file=',file))
		save(brl, file=file)
		stop()			
	}	
	if(0)	#plot composition of selected MSM clusters
	{
		hivc.prog.eval.clustering.bias()
	}	
	if(0)
	{
		hivc.prog.recombination.checkcandidates()
	}
	if(1)
	{
		project.hivc.clustering.compare.NoDR.to.NoRecombNoDR()
		project.seq.dataset.mDR.mRC.mSH.pLANL()
	}
	if(0)	#min brl to get a transmission cascade from brl matrix
	{
		brlmat	<- matrix( c(0,0.1,0.1, 	0.1,0,0.2,	0.1,0.2,0), ncol=3, nrow=3)
		brlmat	<- matrix( c(0,0.1,0.1,0.2, 	0.1,0,0.2,0.3,	0.1,0.2,0,0.1,	0.2,0.3,0.1,0), ncol=4, nrow=4)
		#brlmat	<- matrix( c(0,0.1,0.5,0.5, 	0.1,0,0.5,0.5,	0.5,0.5,0,0.1,	0.5,0.5,0.1,0), ncol=4, nrow=4)
		#brlmat	<- matrix( c(0,0.5,0.5, 	0.5,0,0.1,	0.5,0.1,0), ncol=3, nrow=3)
		#brlmat	<- matrix( c(0,0.5,0.1, 	0.5,0,0.5,	0.1,0.5,0), ncol=3, nrow=3)
		#brlmat	<- 0.5
		#brlmat	<- matrix(2, ncol=1, nrow=1)
		#brlmat	<- brlmat[upper.tri(brlmat)]
		brl		<- hivc.clu.min.transmission.cascade(brlmat)
		str(brl)
		stop()
	}
	if(0)	#test clustering on simple test tree
	{		
		ph	<- "(Wulfeniopsis:0.196108,(((TN_alpinus:0.459325,TN_grandiflora:0.259364)1.00:0.313204,uniflora:1.155678)1.00:0.160549,(((TP_angustibracteata:0.054609,(TN_brevituba:0.085086,TP_stolonifera:0.086001)0.76:0.035958)1.00:0.231339,(((axillare:0.017540,liukiuense:0.018503)0.96:0.038019,stenostachyum:0.049803)1.00:0.083104,virginicum:0.073686)1.00:0.103843)1.00:0.086965,(carinthiaca:0.018150,orientalis:0.019697)1.00:0.194784)1.00:0.077110)1.00:0.199516,(((((abyssinica:0.077714,glandulosa:0.063758)1.00:0.152861,((((allionii:0.067154,(morrisonicola:0.033595,officinalis:0.067266)1.00:0.055175)1.00:0.090694,(alpina:0.051894,baumgartenii:0.024152,(bellidioides:0.016996,nutans:0.063292)0.68:0.031661,urticifolia:0.032044)0.96:0.036973,aphylla:0.117223)0.67:0.033757,(((japonensis:0.018053,miqueliana:0.033676)1.00:0.160576,vandellioides:0.099761)0.69:0.036188,montana:0.050690)1.00:0.058380)1.00:0.115874,scutellata:0.232093)0.99:0.055014)1.00:0.209754,((((((acinifolia:0.112279,reuterana:0.108698)0.94:0.055829,pusilla:0.110550)1.00:0.230282,((davisii:0.053261,serpyllifolia:0.087290)0.89:0.036820,(gentianoides:0.035798,schistosa:0.038522)0.95:0.039292)1.00:0.092830)1.00:0.169662,(((anagalloides:0.018007,scardica:0.017167)1.00:0.135357,peregrina:0.120179)1.00:0.098045,beccabunga:0.069515)1.00:0.103473)1.00:0.287909,(((((((((((agrestis:0.017079,filiformis:0.018923)0.94:0.041802,ceratocarpa:0.111521)1.00:0.072991,amoena:0.229452,(((argute_serrata:0.017952,campylopoda:0.075210)0.64:0.034411,capillipes:0.022412)0.59:0.034547,biloba:0.037143)1.00:0.141513,intercedens:0.339760,((opaca:0.019779,persica:0.035744)0.94:0.038558,polita:0.036762)1.00:0.108620,rubrifolia:0.186799)1.00:0.144789,(((bombycina_11:0.033926,bombycina_bol:0.035290,cuneifolia:0.017300,jacquinii:0.054249,oltensis:0.045755,paederotae:0.051579,turrilliana:0.017117)0.85:0.049052,czerniakowskiana:0.089983)0.93:0.051111,farinosa:0.138075)1.00:0.080565)1.00:0.104525,((albiflora:0.017984,ciliata_Anna:0.032685,vandewateri:0.017610)0.97:0.045649,arguta:0.063057,(catarractae:0.022789,decora:0.049785)0.96:0.048220,((cheesemanii:0.040125,cupressoides:0.146538)1.00:0.067761,macrantha:0.038130)1.00:0.088158,(densifolia:0.090044,formosa:0.116180)0.71:0.046353,(elliptica:0.038650,(odora:0.019325,salicornioides:0.021228)0.94:0.042950,salicifolia:0.020829)0.92:0.043978,(nivea:0.070429,(papuana:0.035003,tubata:0.031140)0.98:0.064379)0.93:0.065336,raoulii:0.109101)0.97:0.076607)0.93:0.085835,chamaepithyoides:0.485601)0.57:0.072713,(ciliata_157:0.069943,lanuginosa:0.052833)1.00:0.098638,(densiflora:0.069429,macrostemon:0.118926)0.92:0.124911,(fruticulosa:0.086891,saturejoides:0.041181)0.94:0.086148,kellererii:0.083762,lanosa:0.263033,mampodrensis:0.103384,nummularia:0.191180,pontica:0.128944,thessalica:0.129197)0.65:0.031006,(arvensis:0.342138,(((((chamaedrys:0.043720,micans:0.032021,vindobonensis:0.033309)0.51:0.034053,micrantha:0.019084)0.64:0.037906,krumovii:0.020175)1.00:0.103875,verna:0.254017)0.81:0.057105,magna:0.112657)1.00:0.104070)1.00:0.101845)1.00:0.149208,(((aznavourii:0.664103,glauca:0.405588)0.85:0.209945,praecox:0.447238)1.00:0.185614,(donii:0.260827,triphyllos:0.176032)1.00:0.194928)1.00:0.611079)0.74:0.055152,((crista:0.591702,(((cymbalaria_Avlan:0.017401,panormitana:0.017609)1.00:0.229508,((cymbalaria_Istanbul:0.028379,trichadena_332:0.016891,trichadena_Mugla:0.019131)1.00:0.196417,lycica_333:0.146772)1.00:0.097646,lycica_192:0.154877)1.00:0.234748,(((hederifolia:0.018068,triloba:0.075784)1.00:0.084865,(sibthorpioides:0.122542,sublobata:0.136951)1.00:0.074683)0.89:0.043623,stewartii:0.040679)1.00:0.596859)1.00:0.237324)0.58:0.057120,javanica:0.133802)1.00:0.137214)1.00:0.269201,(missurica:0.016685,rubra:0.019696)1.00:0.351184)0.54:0.058275)0.52:0.062485,((dahurica:0.023542,longifolia:0.016484,spicata:0.018125)0.95:0.042294,(nakaiana:0.016270,schmidtiana:0.058451)0.88:0.037207)1.00:0.261643)0.55:0.056458)1.00:0.229509,kurrooa:0.100611)0.74:0.068198,(bonarota:0.040842,lutea:0.115316)1.00:0.241657)0.99:0.085772);"
		ph 	<- ladderize( read.tree(text = ph) )				
		#read bootstrap support values
		thresh.bs						<- 0.9
		ph.node.bs						<- as.numeric( ph$node.label )		
		ph.node.bs[is.na(ph.node.bs)]	<- 0
		ph.node.bs[c(13,15,27,41,43)]	<- thresh.bs-0.05*seq(0.01,length.out=5)
		ph.node.bs[ph.node.bs==1]		<- seq_along(which(ph.node.bs==1))*0.005 + 0.7
		ph$node.label					<- ph.node.bs
		#read patristic distances
		stat.fun						<- hivc.clu.min.transmission.cascade
		#stat.fun						<- max
		dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)
		print(dist.brl)
		thresh.brl						<- quantile(dist.brl,seq(0.1,1,by=0.05))["80%"]
		#produce clustering 
		clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
		print(clustering)		
		hivc.clu.plot(ph, clustering[["clu.mem"]], highlight.edge.of.tiplabel=c("TN_","TP_"), highlight.edge.of.tiplabel.col= c("red","blue") )
		#produce some tip states
		ph.tip.state<- rep(1:20, each=ceiling( length( ph$tip.label )/20 ))[seq_len(Ntip(ph))]
		states		<- data.table(state.text=1:20, state.col=rainbow(20))
		states		<- states[ph.tip.state,]
		hivc.clu.plot.tiplabels( seq_len(Ntip(ph)), matrix(states[,state.text], nrow=1,ncol=nrow(states)), matrix(states[,state.col], nrow=1,ncol=nrow(states)), cex=0.4, adj=c(-1,0.5) )
				
		#get nodes in each cluster
		#clu.nodes		<- sapply(unique(clustering[["clu.mem"]])[-1],function(clu){		which(clustering[["clu.mem"]]==clu)	})
		#get boostrap support of index cases
		clu.idx			<- clustering[["clu.idx"]]-Ntip(ph)
		print(clu.idx)
		print(ph.node.bs[clu.idx])
		stop()
	}
	if(0)
	{
		#plot properties of concentrated clusters by region
		verbose							<- 1
		thresh.bs						<- 0.85
		thresh.brl						<- 0.06				
		
		infile		<- paste("ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_withinpatientclusters_bs",thresh.bs*100,"_brlmed",thresh.brl*100,"_clusterinfo",sep='')
		insignat	<- "Sat_Jun_16_17/23/46_2013"
		file		<- paste(dir.name,"tmp",paste(infile,'_',gsub('/',':',insignat),".R",sep=''),sep='/')
		if(verbose)	cat(paste("read cluster info from file",file))
		load(file)
		if(verbose)	cat(paste("found seq, n=",nrow(df.seqinfo)))		
		print(df.seqinfo)
		df.seqinfo	<- subset(df.seqinfo, !is.na(cluster) )
		if(verbose)	cat(paste("found seq in clu, n=",nrow(df.seqinfo)))
		
		#remove within patient clusters before out-of-NL taken out		
		set(df.seqinfo, which( df.seqinfo[,is.na(Patient)] ), "Patient", "NA")
		tmp			<- df.seqinfo[, list(isWithinCluster= length(unique(Patient))==1), by=cluster]
		df.seqinfo	<- merge( subset(tmp,!isWithinCluster), df.seqinfo, by="cluster" )
		#from between patient clusters, remove within patient seqs and foreign seqs
		df.pat		<- df.seqinfo[,list(CountryInfection=CountryInfection[1], Trm=Trm[1], Sex=Sex[1], NegT=NegT[1], AnyPos_T1=AnyPos_T1[1], RegionHospital=RegionHospital[1]), by=Patient]			
		df.clu		<- df.seqinfo[, list(Patient=unique(Patient)), by=cluster]		
		df.clu		<- merge(df.clu, df.pat, all.x=1, by="Patient")
		
		
		#determine if cluster concentrated in region, and where
		setkey(df.clu, cluster)
		if(verbose)	cat(paste("found pat in clu, n=",nrow(df.clu)))
		clu.reg		<- table( df.clu[,cluster,RegionHospital] )
		clu.reg		<- clu.reg / matrix(apply(clu.reg, 2, sum), nrow=nrow(clu.reg), ncol=ncol(clu.reg), byrow=1)				
		clu.regconc	<- apply( clu.reg, 2, function(x)	any(x>0.4) && length(which(x!=0))>2 ) |
					   apply( clu.reg, 2, function(x)	any(x>0.5) && length(which(x!=0))<=2 )
		if(verbose)	cat(paste("number of clusters, n=",length(clu.regconc)))
		if(verbose)	cat(paste("number of spatially concentrated clusters, n=",length(which(clu.regconc))))
		tmp			<- rownames(clu.reg)
		clu.regconc2<- sapply( seq_len(ncol(clu.reg)), function(j)
						{
							ifelse(!clu.regconc[j], NA, tmp[ which.max(clu.reg[,j]) ] )						
						})
		df.clureg	<- data.table(cluster= as.numeric(colnames( clu.reg )), IsRegionConcentrated= clu.regconc, RegionConcentrated=clu.regconc2, key="cluster" )
		
		#determine median AnyPos_T1 time for cluster and add
		tmp			<- df.clu[,list(cluster.PosT=mean(AnyPos_T1, na.rm=T)),by=cluster]				
		df.clureg	<- merge(df.clureg, tmp, all.x=1, by="cluster")

		#merge all
		df.seqinfo	<- merge( df.seqinfo, df.clureg, all.x=1, by="cluster" )
						
		
		#plot Bezemer style
		file				<- paste(dir.name,"tmp",paste(infile,"_spatialvariation_",gsub('/',':',insignat),".pdf",sep=''),sep='/')
		df					<- df.seqinfo
		df[,time:=AnyPos_T1 ]		
		set(df, which(df[,!IsRegionConcentrated]), "RegionConcentrated", "Bridge")
		tmp					<- numeric(nrow(df))
		tmp2				<- c("N","E","S","W","Amst","Bridge")
		for( i in seq_along(tmp2))
			tmp[which(df[,RegionConcentrated==tmp2[i]])]	<- i
		df[,sort1:=factor(tmp, levels=1:6, labels=tmp2) ]
		df[,sort2:=cluster.PosT]		
		tmp					<- numeric(nrow(df))		
		tmp2				<- c("N","E","S","W","Amst")
		for( i in seq_along(tmp2))
			tmp[which(df[,RegionHospital==tmp2[i]])]	<- i
		df[,covariate:=factor(tmp, levels=1:5, labels=tmp2) ]
		xlab				<- "first HIV+ event"
		ylab				<- paste("clusters BS=",thresh.bs*100," BRL.med=",thresh.brl*100," version=",gsub('/',':',insignat),sep="")
		cex.points			<- 0.5
		
		#extract clusters one by one in desired order
		tmp					<- df[,list(sort1=sort1[1], sort2=sort2[1]),by=cluster]
		clusters.sortby1	<- as.numeric( tmp[,sort1] )
		clusters.sortby2	<- as.numeric( tmp[,sort2] )
		clusters.sortby		<- order(clusters.sortby1,clusters.sortby2, decreasing=F)
		clusters			<- tmp[clusters.sortby,cluster]								
		clusters.levels		<- levels(df[,covariate])					
		clusters			<- lapply(clusters,function(x)	subset(df,cluster==x,c(time,covariate))		)
		
		ylim				<- c(-2,length(clusters))
		cols				<- brewer.pal(length(clusters.levels),"Set1")	#[c(2,1,3,4,5)]		
#cols			<- c("green","red","blue","grey30")
		cols				<- sapply(cols, function(x) my.fade.col(x,0.7))
		pch					<- c(rep(19,length(clusters.levels)))
		xlim				<- range( df[,time], na.rm=1 )
		xlim[1]				<- xlim[1]-700
		
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=14)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab=xlab)
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		dummy<- sapply(seq_along(clusters),function(i)
				{							
					cluster.cex	<- cex.points	#(1 - nrow(subset(clusters[[i]],hregion=="unknown")) / nrow(clusters[[i]])) * cex.points
					cluster.ix	<- order(as.numeric(clusters[[i]][,time]))
					cluster.x	<- as.numeric(clusters[[i]][,time])[cluster.ix]										
					cluster.y	<- rep(i,nrow(clusters[[i]]))
					cluster.z	<- as.numeric(clusters[[i]][,covariate])[cluster.ix]
					#print(clusters[[i]][,covariate])
					#if(cluster.z[1]!=4)
					#lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
					points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=pch[1], cex=cluster.cex )										
				})	
		pch[pch==19]	<- 21
		legend("topleft",bty='n',pt.bg=cols,pch=pch,legend=levels(df[,covariate]), col=rep("transparent",length(clusters.levels)))				
		dev.off()	
		stop()
	}
	if(0)
	{
		#check BEEHIVE sequences
		indircov	<- paste(DATA,"derived",sep='/')
		infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"							
		load(paste(indircov,'/',infilecov,".R",sep=''))
		
		df.bee		<- as.data.table(read.csv("~/duke/2013_HIV_NL/ATHENA_2013/data/BEEHIVE1_data.csv", stringsAsFactors=0))
		setnames(df.bee,"mcode","Patient")
		df.bee		<- merge(df.all, subset(df.bee,select=Patient), by="Patient")
		save(df.bee, file="~/duke/2013_HIV_NL/ATHENA_2013/data/BEEHIVE1_data.R")
		
	}
	if(0)
	{		
		verbose		<- 1
		resume		<- 1
		#
		# precompute clustering stuff		
		#
		patient.n	<- 15700
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
		insignat	<- "Sat_Jun_16_17/23/46_2013"
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"			
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		indircov	<- paste(DATA,"derived",sep='/')
		infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"							
		argv		<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
		argv		<<- unlist(strsplit(argv,' '))
		clu.pre		<- hivc.prog.get.clustering.precompute()
		#
		# evaluate TPTN for various thresholds
		#
		if(verbose) cat(paste("compute TPTN for dist.brl.casc"))
		argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume)
		argv			<<- unlist(strsplit(argv,' '))
		clu.tptn.casc	<- hivc.prog.get.clustering.TPTN(clu.pre=clu.pre)
		if(verbose) cat(paste("compute TPTN for dist.brl.max"))
		argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.max", patient.n=patient.n, resume=resume)
		argv			<<- unlist(strsplit(argv,' '))
		clu.tptn.max	<- hivc.prog.get.clustering.TPTN(clu.pre=clu.pre)
		#
		# compare dist.brl vs dist.max in terms of total patients in clusters
		#
		file			<- paste(indir, paste(infile,"_clust_cascvsmax_patients_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
		pdf(file=file,width=5,height=5)
		par(mar=c(4,4,1,0))
		dbwpat.cmp		<- list(	max	=lapply(seq_along(select.bs), function(i)	clu.tptn.max[["clusters.dbwpat"]][[select.brl[i]]][[select.bs[i]]]		),
									casc=lapply(seq_along(select.bs), function(i)	clu.tptn.casc[["clusters.dbwpat"]][[select.brl[i]]][[select.bs[i]]]		)	)
		ylim			<- c(0,3100)
		#sapply(dbwpat.cmp, function(x)	sapply(x, function(z) sum(as.numeric(names(z))*z)  ) )
		xlim			<- c(1,max( sapply(dbwpat.cmp, function(x)	sapply(x, function(z) max(as.numeric(names(z)))) ) ))
		
		
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim, xlab="cluster size", ylab="patients in clusters of size <=x")
		lapply(seq_along(dbwpat.cmp), function(i)
				{
					lapply(seq_along(dbwpat.cmp[[i]]),function(j)
							{	
								z	<- dbwpat.cmp[[i]][[j]][-1]
								z2	<- as.numeric(names(z))								
								lines(z2,cumsum(z*z2), col=cols[i], lty=j)
							})
				})
		legend("topleft",bty='n', border=NA, fill= cols, legend= c("max","casc"))
		legend("topright",bty='n', border=NA, lty= 1:2, legend= c("BS=0.8, BRL=0.1","BS=0.95, BRL=0.04"))
		dev.off()
		#
		# compare dist.brl vs dist.max in terms of %cov 
		#
		file			<- paste(indir, paste(infile,"_clust_cascvsmax_covepi_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
		pdf(file=file,width=5,height=5)
		par(mar=c(4,4,1,0))
		ylim			<- c(0,0.2)				
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim, xlab="cluster size", ylab="%coverage (of epi) in clusters of size <=x")
		lapply(seq_along(dbwpat.cmp), function(i)
				{
					lapply(seq_along(dbwpat.cmp[[i]]),function(j)
							{	
								z	<- dbwpat.cmp[[i]][[j]][-1]
								z2	<- as.numeric(names(z))								
								lines(z2,cumsum(z*z2)/patient.n, col=cols[i], lty=j)
							})
				})		
		legend("topleft",bty='n', border=NA, fill= cols, legend= c("max","casc"))
		legend("topright",bty='n', border=NA, lty= 1:2, legend= c("BS=0.8, BRL=0.1","BS=0.95, BRL=0.04"))
		dev.off()
		#
		# compare dist.brl vs dist.max in terms of %cov vs %fp for BS=0.8/BRL=0.1/casc  vs BS=0.95/BRL=0.05/max 
		#
		covfp.cmp.x	<- rbind( clu.tptn.casc[["fp.by.all"]]["0.8",], clu.tptn.max[["fp.by.all"]]["0.95",]	)
		covfp.cmp.y	<- rbind( clu.tptn.casc[["clusters.cov.epidemic"]]["0.8",], clu.tptn.max[["clusters.cov.epidemic"]]["0.95",]	)		
		
		file			<- paste(indir, paste(infile,"_clust_cascvsmax_covepifp_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
		pdf(file=file,width=5,height=5)
		par(mar=c(4,4,1,0))		
		plot(1,1,type='n',bty='n',xlim=range(c(0.01,covfp.cmp.x)),ylim=range(c(0.2,covfp.cmp.y)),xlab="%FP (among all)",ylab="%coverage (of epi)")
		dummy	<- sapply(seq_len(nrow(covfp.cmp.x)),function(i)
				{					
					points(covfp.cmp.x[i,],covfp.cmp.y[i,],col=cols[i],type='b')
					text(covfp.cmp.x[i,],covfp.cmp.y[i,],col=cols[i],labels=as.numeric(colnames(covfp.cmp.y)),adj=c(-0.8,0.5),cex=0.5)
				})
		legend("bottomright",border=NA,bty='n',fill=cols,legend=c("max","casc"))
		dev.off()				
	}
	if(0)
	{
		#get branch lengths between F2F transmissions, where there must be a missed intermediary
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"hetsamegender",sep='')
		outsignat				<- insignat		
		plot.file				<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep='')		
		clu.female2female		<- hivc.clu.getplot.female2female( clu.pre$ph, clu$clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep='') )				
		clu.missedintermediary	<- names(cluphy.brl.bwpat) %in% as.character( unique( clu.female2female$cluphy.df[,cluster] ) )
		brl.missedintermediary	<- unlist( lapply( which(clu.missedintermediary), function(i)	cluphy.brl.bwpat[[i]] ) )
		brl.others				<- na.omit( unlist( lapply( which(!clu.missedintermediary), function(i)	cluphy.brl.bwpat[[i]] ) ) )
		brl.breaks				<- seq(0,max(brl.others,brl.missedintermediary)*1.1,by=0.02)
		brl.cols				<- sapply( brewer.pal(3, "Set1"), function(x) my.fade.col(x,0.5) )
		
		par(mfcol=c(1,2))
		hist( brl.others, breaks=brl.breaks, col=brl.cols[1], border=NA, freq=T, main='', xlab="branch lengths, others"  )
		#legend("topright", fill= brl.cols[1:2], legend= c("others","missed intermediary"), bty='n', border=NA)
		hist( brl.missedintermediary, breaks=brl.breaks, col=brl.cols[2], border=NA, freq=T, main='', xlab="branch lengths, missed intermediary" )		
		par(mfcol=c(1,1))
		
		
		#highlight seq of same patient not in cluster
		#select clusters with TN and P51

		
		stop()
		
	}
	if(0)	#count how many unlinked pairs in clustering
	{		
		verbose	<- 1
		file	<- paste(dir.name,"derived/ATHENA_2013_03_Unlinked_SeroConv_Dead_UnlinkedAll.R", sep='/')		
		if(verbose)	cat(paste("read file",file))
		load(file)
		#str(unlinked.bytime)
		
		infile		<- "ATHENA_2013_03_FirstCurSequences_PROTRT_examlbs100"
		signat.in	<- "Fri_May_24_12/59/06_2013"
		signat.out	<- "Fri_Jun_07_09/59/23_2013"
		file		<- paste(dir.name,"tmp",paste(infile,'_',gsub('/',':',signat.in),".newick",sep=''),sep='/')
		cat(paste("read file",file))
		
		#
		#set up tree, get boostrap values and patristic distances between leaves
		#
		ph								<- ladderize( read.tree(file) )		
		ph.node.bs						<- as.numeric( ph$node.label )
		ph.node.bs[is.na(ph.node.bs)]	<- 0
		ph.node.bs						<- ph.node.bs/100
		ph$node.label					<- ph.node.bs
		#print(quantile(ph.node.bs,seq(0.1,1,by=0.1)))		
		dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf")		#read patristic distances -- this is the expensive step but still does not take very long
		#print(quantile(dist.brl,seq(0.1,1,by=0.1)))
		
		#
		#convert truly unlinked pairs from Patient name to ph node index
		#
		df.tips							<- data.table(Node=seq_len(Ntip(ph)), Patient=ph$tip.label )
		setkey(df.tips, Patient)																#use data.table to speed up search
		df.tips							<- df.all[df.tips]
		ph.unlinked.dead				<- lapply(seq_along(unlinked.bytime), function(j)
											{
												as.numeric( sapply(seq_along(unlinked.bytime[[j]]), function(i)	subset(df.tips, unlinked.bytime[[j]][i]==Patient, Node) ))					
											})		
		ph.unlinked.seroneg				<- as.numeric( sapply(names(unlinked.bytime), function(x)	subset(df.tips, x==Patient, Node) ))		
		tmp								<- sort( ph.unlinked.seroneg, index.return=1)			#sort ph.unlinked.seroneg and return index, so we can also sort ph.unlinked.dead
		ph.unlinked.seroneg				<- tmp$x
		names(ph.unlinked.seroneg)		<- names(unlinked.bytime)[tmp$ix]
		ph.unlinked.dead				<- lapply(tmp$ix, function(i){		ph.unlinked.dead[[i]]	})	
		names(ph.unlinked.dead)			<- names(unlinked.bytime)[tmp$ix]		
		ph.unlinked.seroneg				<- data.table(PhNode=ph.unlinked.seroneg, PhNodeUnlinked= seq_along(ph.unlinked.seroneg), Patient= names(ph.unlinked.seroneg))
		setkey(ph.unlinked.seroneg, Patient)													#set key to take right outer join with df.serocon -- temporarily destroy correspondance with 'ph.unlinked.dead'
		setkey(df.serocon, Patient)
		ph.unlinked.seroneg				<- df.serocon[ph.unlinked.seroneg]		
		setkey(ph.unlinked.seroneg, PhNode)			#use data.table to speed up search -- restore correspondance with 'ph.unlinked.dead'		
		#print(ph.unlinked.seroneg); print(ph.unlinked.dead); stop()

		#
		#determine threshold for patristic distances based on truly unlinked pairs
		#
		thresh.brl	<- NULL
		thresh.bs	<- 0.9
		ans			<- hivc.clu.clusterbytypeIerror(ph, dist.brl, ph.node.bs, ph.unlinked.seroneg, ph.unlinked.dead, thresh.brl=thresh.brl, thresh.bs=thresh.bs, level= 0.001, verbose=1)		
		#
		#determine quick estimate of expected false pos for these thresholds if clustering was random
		#				
		check		<- hivc.clu.exp.typeIerror.randomclu(ph, dist.brl, ph.node.bs, ph.unlinked.seroneg, ph.unlinked.dead, ans[["thresh.brl"]], ans[["thresh.bs"]])
		cat(paste("\n#clusters with at least one FP=",check[["fp.n"]]," , %clusters with at least one FP=",check[["fp.rate"]],sep=''))
		#
		#max cluster size distribution
		barplot(cumsum(h$counts), names.arg=seq_along(h$counts)+1, axisnames=1, xlab="maximum cluster size")
		#
		#number of sequences in cluster, and %
		cat(paste("\n#seq in cluster=",sum(ans[["clu"]][["size.tips"]]),", %seq in cluster=",sum(ans[["clu"]][["size.tips"]])/Ntip(ph),sep='')) 		
		#
		#plot final clusters, save  clusters
		#								
		file		<- paste(dir.name,"tmp",paste(infile,"_clu_",gsub('/',':',signat.out),".pdf",sep=''),sep='/')
		cat(paste("plot clusters on tree to file",file))
		hivc.clu.plot(ph, ans[["clu"]][["clu.mem"]], file=file, pdf.scaley=25)
		file		<- paste(dir.name,"tmp",paste(infile,"_clu_",gsub('/',':',signat.out),".R",sep=''),sep='/')
		cat(paste("save analysis to file",file))
		save(ans,check,ph,dist.brl,ph.node.bs,ph.unlinked.seroneg,ph.unlinked.dead, file=file)				
	}
}
######################################################################################
project.hivc.clustering.selectparticularclusters<- function()
{	
	if(1)
	{
		verbose		<- 1
		resume		<- 1
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"			
		insignat	<- "Thu_Aug_01_17/05/23_2013"
		indircov	<- paste(DATA,"derived",sep='/')
		infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"
		opt.brl		<- "dist.brl.casc" 
		thresh.brl	<- 0.096
		thresh.bs	<- 0.8
		argv		<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
		argv		<<- unlist(strsplit(argv,' '))
		clu.pre		<- hivc.prog.get.clustering.precompute()
		argv		<<- hivc.cmd.clustering(indir, infile, insignat, opt.brl, thresh.brl, thresh.bs, resume=resume)				 
		argv		<<- unlist(strsplit(argv,' '))
		clu			<- hivc.prog.get.clustering()		
		#
		# remove singletons
		#
		if(verbose) cat(paste("\nnumber of seq in tree is n=", nrow(clu$df.cluinfo)))
		df.cluinfo	<- subset(clu$df.seqinfo, !is.na(cluster) )
		if(verbose) cat(paste("\nnumber of seq in clusters is n=", nrow(df.cluinfo)))
		if(verbose) cat(paste("\nnumber of clusters is n=", length(unique(df.cluinfo[,cluster]))))
		#
		# remove within patient clusters
		#
		tmp			<- subset(df.cluinfo[,list(clu.is.bwpat=length(unique(Patient))>1),by="cluster"], clu.is.bwpat, cluster )
		df.cluinfo	<- merge(tmp, df.cluinfo, by="cluster", all.x=1)
		#
		# plot merged clusters that have shared patients
		#
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"sharingpatclu2",sep='')
		outsignat		<- insignat									
		tmp				<- hivc.clu.getplot.potentialsuperinfections(clu.pre$ph, clu$clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep='') )		 		
		# identify potential superinfections by name
		cluphy.df		<- tmp$cluphy.
		cluphy.df		<- subset(cluphy.df, select=c(cluster, FASTASampleCode, Patient, PosSeqT))
		tmp				<- cluphy.df[, list(count= length(unique(cluster)), cluster=cluster, FASTASampleCode=FASTASampleCode, PosSeqT=PosSeqT),by="Patient"]
		tmp				<- subset(tmp,count>1)
		unique(tmp[,Patient])
	}
	if(1)
	{	
		msm<- hivc.prog.get.clustering.MSM()		
		#
		#1) plot clusters with 		small brl / npat 	--> explosive for targeted testing -- mostly acute -- serial or starlike or what ?
		#
		df.cluinfo				<- msm$df.cluinfo
		tmp						<- df.cluinfo[,list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], select=clu.bwpat.medbrl[1]/clu.npat[1]),by="cluster"]
		cumsum( table( tmp[,clu.npat] ) / nrow(tmp) )
		tmp						<- subset(tmp, clu.npat>4)
		tmp						<- subset( tmp, select<quantile( tmp[,select], probs=0.2 ))
		cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
		if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))		
		cluphy.subtrees			<- lapply( as.character(tmp[,cluster]), function(x)  msm$cluphy.subtrees[[x]]	)
		names(cluphy.subtrees)	<- as.character(tmp[,cluster])
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectexplosive",sep='')
		outsignat				<- insignat									
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=10, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.3, pdf.xlim=0.36)
		cluphy					<- tmp$cluphy		
		cluphy.df				<- subset(cluphy.df, select=c(cluster, FASTASampleCode, Patient,    PosSeqT,   DateBorn, Sex,  CountryBorn, CountryInfection, RegionHospital))			
		outfile					<- paste(DATA,'/',"tmp",'/',infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectexplosive.csv",sep='')
		write.table(cluphy.df, file=outfile, sep=",", row.names=F)
		outfile					<- paste(DATA,'/',"tmp",'/',infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectexplosive.R",sep='')
		save(cluphy.df, file=outfile)
		#
		#2) plot clusters with 		npat>4 	very small acute	small brl / npat 	--> explosive for targeted testing -- non - acute
		#
		df.cluinfo				<- msm$df.cluinfo
		tmp						<- df.cluinfo[,list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], select=clu.bwpat.medbrl[1]/clu.npat[1]/clu.npat[1]),by="cluster"]
		cumsum( table( tmp[,clu.npat] ) / nrow(tmp) )
		tmp						<- subset(tmp, clu.npat>4)
		tmp						<- subset( tmp,  clu.fPossAcute<quantile( tmp[,clu.fPossAcute], probs=0.2 ) )
		cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
		if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))		
		cluphy.subtrees			<- lapply( as.character(tmp[,cluster]), function(x)  msm$cluphy.subtrees[[x]]	)
		names(cluphy.subtrees)	<- as.character(tmp[,cluster])
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectlargenonacute",sep='')
		outsignat				<- insignat									
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=10, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.3, pdf.xlim=0.36)
		#
		#3) plot clusters with 		high VLI after treat			--> likely to infect -- check manually ?	large clusters are a subset of (2) so plot all
		#
		df.cluinfo				<- msm$df.cluinfo
		tmp						<- df.cluinfo[,	{
					z			<- which(PosSeqT>AnyT_T1)
					hlRNA.mx	<- ifelse(length(z), max(lRNA_aTS[z]), 0)
					hlRNA.sm	<- ifelse(length(z), sum(lRNA_aTS[z]), 0)
					list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], hlRNA.mx=hlRNA.mx, hlRNA.sm=hlRNA.sm)
				},by="cluster"]
		tmp						<- subset( tmp, hlRNA.sm>=quantile(tmp[,hlRNA.sm], probs=0.8) )
		cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
		if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))		
		cluphy.subtrees			<- lapply( as.character(tmp[,cluster]), function(x)  msm$cluphy.subtrees[[x]]	)
		names(cluphy.subtrees)	<- as.character(tmp[,cluster])
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selecthighVLduringtreat",sep='')
		outsignat				<- insignat									
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=15, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.25, pdf.xlim=0.36)
		#
		#4) plot clusters with 		long treatment interruptions			--> likely to infect -- check manually ?	large clusters are a subset of (2) so plot all
		#
		df.cluinfo				<- msm$df.cluinfo
		tmp						<- df.cluinfo[,	{
					z		<- which(!is.na(AnyT_T1))
					TrI.mx	<- ifelse(length(z),	max(TrImo_bTS[z] + TrImo_aTS[z]), 0)
					TrI.sm	<- ifelse(length(z),	sum(TrImo_bTS[z] + TrImo_aTS[z]), 0)
					list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], TrI.mx=TrI.mx, TrI.sm=TrI.sm)
				},by="cluster"]
		tmp						<- subset(tmp, TrI.mx>=quantile(tmp[,TrI.mx], probs=0.8)  & clu.fPossAcute<0.7 )
		cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
		if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))		
		cluphy.subtrees			<- lapply( as.character(tmp[,cluster]), function(x)  msm$cluphy.subtrees[[x]]	)
		names(cluphy.subtrees)	<- as.character(tmp[,cluster])
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectLongTRI",sep='')
		outsignat				<- insignat									
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=12, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.27, pdf.xlim=0.36)
		#
		#5) plot clusters with 		NegT			--> might help to better understand what is going on ?	
		#
		df.cluinfo				<- msm$df.cluinfo
		tmp						<- df.cluinfo[,	list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], fNegT=length(which(!is.na(NegT))) / clu.ntip[1]),by="cluster"]										
		tmp						<- subset(tmp, fNegT>=quantile(tmp[,fNegT], probs=0.8) )
		cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
		if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))		
		cluphy.subtrees			<- lapply( as.character(tmp[,cluster]), function(x)  msm$cluphy.subtrees[[x]]	)
		names(cluphy.subtrees)	<- as.character(tmp[,cluster])
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectAllWithNegT",sep='')
		outsignat				<- insignat									
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=10, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.3, pdf.xlim=0.36)
		#
		#6) plot clusters with 		long treatment interruptions			--> long TRI vs Acute ?	
		#		
		df.cluinfo				<- msm$df.cluinfo
		tmp						<- df.cluinfo[,	list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], fNegT=length(which(!is.na(NegT))) / clu.ntip[1]),by="cluster"]										
		tmp						<- subset(tmp, fNegT>=quantile(tmp[,fNegT], probs=0.8) )
		tmp						<- subset(tmp, clu.npat>3)
		cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
		if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))		
		cluphy.subtrees			<- lapply( as.character(tmp[,cluster]), function(x)  msm$cluphy.subtrees[[x]]	)
		names(cluphy.subtrees)	<- as.character(tmp[,cluster])
		outdir					<- indir
		outfile					<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_selectLargeLongTRI",sep='')
		outsignat				<- insignat									
		tmp						<- hivc.clu.polyphyletic.clusters(cluphy.df, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=3, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.3, pdf.xlim=0.36)
	}
}
######################################################################################
project.hivc.examl.median.brl<- function()
{
	require(reshape2)
	require(data.table)
	require(ape)	
	#
	#	get data relating to study population (subtype B sequ)
	#	
	indir					<- paste(DATA,"fisheretal_data",sep='/')		
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	outdir					<- paste(DATA,"fisheretal",sep='/')
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	infiletree				<- paste(infile,"examlbs500",sep="_")
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	infile.cov.study		<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infile.viro.study		<- paste(indircov,"ATHENA_2013_03_Viro.R",sep='/')
	infile.immu.study		<- paste(indircov,"ATHENA_2013_03_Immu.R",sep='/')
	infile.treatment.study	<- paste(indircov,"ATHENA_2013_03_Regimens.R",sep='/')	
	infile.cov.all			<- "ATHENA_2013_03_AllSeqPatientCovariates_AllMSM"
	infile.viro.all			<- paste(indircov,"ATHENA_2013_03_Viro_AllMSM.R",sep='/')
	infile.immu.all			<- paste(indircov,"ATHENA_2013_03_Immu_AllMSM.R",sep='/')
	infile.treatment.all	<- paste(indircov,"ATHENA_2013_03_Regimens_AllMSM.R",sep='/')	
	infile.trm.model		<- NA
	t.period				<- 1/8
	t.recent.startctime		<- hivc.db.Date2numeric(as.Date("1996-07-15"))
	t.recent.startctime		<- floor(t.recent.startctime) + floor( (t.recent.startctime%%1)*100 %/% (t.period*100) ) * t.period
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))	
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	method.recentctime		<- '2011-01-01'
	t.recent.endctime		<- hivc.db.Date2numeric(as.Date(method.recentctime))	
	t.recent.endctime		<- floor(t.recent.endctime) + floor( (t.recent.endctime%%1)*100 %/% (t.period*100) ) * t.period		
	adjust.AcuteByNegT		<- 1
	method.use.AcuteSpec	<- 1
	tmp				<- project.athena.Fisheretal.select.denominator(	indir, infile, insignat, indircov, infile.cov.study, infile.viro.study, infile.immu.study, infile.treatment.study, infiletree=infiletree, 
																		adjust.AcuteByNegT=adjust.AcuteByNegT, adjust.NegT4Acute=NA, adjust.NegTByDetectability=0.25, adjust.minSCwindow=0.25, adjust.AcuteSelect=c('Yes','Maybe'), use.AcuteSpec=method.use.AcuteSpec, t.recent.endctime=t.recent.endctime, t.recent.startctime=t.recent.startctime)																
	df.all			<- tmp$df.all							#this is all men only
	df.denom.CLU	<- tmp$df.select
	df.denom.SEQ	<- tmp$df.select.SEQ	
	ri.CLU			<- unique(subset(df.denom.CLU, select=Patient))
	ri.SEQ			<- unique(subset(df.denom.SEQ, select=Patient))	
	df.ris			<- merge(df.all, ri.SEQ, by='Patient')	#this is all recent
	
	
	indir.dtp		<- paste(DATA,"ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlout",sep='/')
	infiles.dtp		<- list.files(indir.dtp)
	infiles.dtp		<- infiles.dtp[ grepl('^ExaML_result.*finaltree.*distTips.R$',infiles.dtp)  ]	
	if(!length(infiles.dtp))	stop('cannot find files matching criteria')
	for(i in seq_along(infiles.dtp))
	{
		infile.dtp	<- infiles.dtp[i]
		cat(paste('\nprocess file',infile.dtp))
		file	<- paste(indir.dtp, '/', infile.dtp, sep='')
		load(file)	#expect brl	
		#	get labels 
		brl.info	<- data.table(FASTASampleCode=attr(brl, "Labels" ))
		brl.info[, IDX:=seq_along(FASTASampleCode)]
		brl.info	<- subset( brl.info, !grepl('PROT+P51',FASTASampleCode, fixed=T) )
		#	get pairs of recent with men
		tmp			<- merge(brl.info, subset(df.ris, select=FASTASampleCode), by='FASTASampleCode')				
		tmp2		<- merge(brl.info, subset(df.all, select=c(FASTASampleCode)), by='FASTASampleCode')
		set(tmp, NULL, 'FASTASampleCode', tmp[,factor(FASTASampleCode)])
		set(tmp2, NULL, 'FASTASampleCode', tmp2[,factor(FASTASampleCode)])
		setnames(tmp, c('FASTASampleCode','IDX'),c('FASTASampleCode_R','IDX_R'))
		setnames(tmp2, c('FASTASampleCode','IDX'),c('FASTASampleCode_T','IDX_T'))
		df.brl		<- as.data.table(expand.grid(IDX_R= tmp[,IDX_R], IDX_T=tmp2[,IDX_T]))
		df.brl		<- merge(df.brl, tmp, by='IDX_R')
		df.brl		<- merge(df.brl, tmp2, by='IDX_T')		
		#	get brl on this tree
		n			<- attr(brl, "Size")
		df.brl		<- subset(df.brl, IDX_R!=IDX_T)	
		#	key: need to make sure we have indices on lower triangle
		tmp			<- df.brl[, which(IDX_T<IDX_R)]
		tmp2		<- df.brl[tmp, IDX_T]
		set(df.brl, tmp, 'IDX_T', df.brl[tmp,IDX_R])
		set(df.brl, tmp, 'IDX_R', tmp2)
		set(df.brl, NULL, 'BRL', df.brl[, brl[my.lower.tri.index(n, IDX_T, IDX_R)]])
		set(df.brl, NULL, c('IDX_T','IDX_R'), NULL)
		file		<- gsub('distTips','distTipsToRec',file)
		cat(paste('\nsave to file=',file))
		save(df.brl,file=file,compress='xz')
		#	
	}
stop()	
	load("/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_data/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_Wed_Dec_18_11:37:00_2013_tpairs.R")
	for(i in seq_along(infiles.dtp))
	{
		i	<- 1
		infile.dtp	<- infiles.dtp[i]
		cat(paste('\nprocess file', infile.dtp))
		infile.exa	<- substr(infile.dtp, 1, nchar(infile.dtp)-11)
		
		#	read phylogeny for which distTips was generated
		file	<- paste( indir.dtp, '/', infile.exa, sep='' )	
		ph		<- read.tree(file)
		#	load brl
		file	<- paste(indir.dtp, '/', infile.dtp, sep='')
		load(file)	#expect brl	
		#
		#	compute branch lengths between infected and all potential transmitters
		#
		df.tpairs.brl		<- copy(df.tpairs)
		tmp					<- unique( subset(df.tpairs.brl, select=FASTASampleCode) )	
		tmp					<- tmp[, list(sc.i=match(FASTASampleCode, ph$tip.label)), by='FASTASampleCode']
		df.tpairs.brl		<- merge(df.tpairs.brl, tmp, by='FASTASampleCode')
		tmp					<- unique( subset(df.tpairs.brl, select=t.FASTASampleCode) )
		tmp					<- tmp[, list(sc.t=match(t.FASTASampleCode, ph$tip.label)), by='t.FASTASampleCode']
		df.tpairs.brl		<- merge(df.tpairs.brl, tmp, by='t.FASTASampleCode')
		tmp					<- df.tpairs.brl[, which(sc.t<sc.i)]
		tmp2				<- df.tpairs.brl[tmp, sc.t]
		set(df.tpairs.brl, tmp, 'sc.t', df.tpairs.brl[tmp,sc.i])
		set(df.tpairs.brl, tmp, 'sc.i', tmp2)
		brl.n				<- attr(brl,'Size')
		df.tpairs.brl		<- df.tpairs.brl[, 	list( brl= brl[ my.lower.tri.index(brl.n, sc.t, sc.i) ]) ,by=c('FASTASampleCode','t.FASTASampleCode')]
		df.tpairs.brl		<- subset(df.tpairs.brl, !is.na(brl))		#some requested brl may be missing because sequences have been deselected (too short or recombinants)
		#
		file				<- paste( indir.dtp, '/', infile.exa, '_distPairs.R',sep='' )
		cat(paste('\nsave to file',file))
		save(df.tpairs.brl, file= file)
		#	remove distTips
		file	<- paste(indir.dtp, '/', infile.dtp, sep='')
		file.remove(file)	
		stop()
	}
	#
	#	part 2 after computation on cluster
	#
	indir		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tmp/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlout_Wed_Dec_18_11:37:00_2013_distPairs'
	infiles		<- list.files(indir)
	infiles		<- infiles[ grepl('*distPairs.R$',infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')
	if(length(infiles))		cat(paste('\nfound files, n=', length(infiles)))
	
	brl.tpairs	<- lapply(seq_along(infiles), function(i)
			{
				infile	<- infiles[i]
				file	<- paste(indir, '/', infile, sep='')
				load(file)	#expect df.tpairs.brl
				tmp		<- regmatches(infile, regexpr('finaltree\\.[0-9]{3}', infile))
				df.tpairs.brl[, BS:= as.numeric(substr(tmp, 11, nchar(tmp)))]
				df.tpairs.brl
			})
	brl.tpairs	<- do.call('rbind', brl.tpairs)	
	setkey(brl.tpairs, FASTASampleCode, t.FASTASampleCode)
	file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tmp/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlout_Wed_Dec_18_11:37:00_2013_distPairs.R'
	save(brl.tpairs, file=file)
	
	file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_data/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlout_Wed_Dec_18_11:37:00_2013_distPairs.R'
	#
	#	plot random sample
	#
	tmp			<- unique(brl.tpairs)
	tmp			<- tmp[sample(nrow(tmp),5e2),]
	tmp			<- merge(subset(tmp, select=c(FASTASampleCode, t.FASTASampleCode)), brl.tpairs, by=c('FASTASampleCode','t.FASTASampleCode'))
	tmp[, interaction:= tmp[, paste(FASTASampleCode, t.FASTASampleCode, sep='::')]]
	setkey(tmp, interaction)
	tmp[, variable:='from bootstrap\nreplicate\nof data']
	tmp			<- rbind(tmp, tmp[, list(FASTASampleCode=FASTASampleCode[1], t.FASTASampleCode=t.FASTASampleCode[1], brl=median(brl), BS=NA_real_, variable='median\npatristic distance'), by='interaction'], use.names=TRUE)
	set(tmp, tmp[, which(BS==0)], 'variable', 'from data')
	#set(tmp, tmp[, which(BS==409)], 'variable', 'maximum\nlikelihood\tree')
	set(tmp, tmp[, which(BS==288)], 'variable', 'with largest\nmaximum likelihood\nacross bootstrap replicates')		#second best tree
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('from data','from bootstrap\nreplicate\nof data','with largest\nmaximum likelihood\nacross bootstrap replicates','median\npatristic distance'))])
	
	ggplot(tmp, aes(x=interaction, y=brl, colour=variable)) + geom_point(size=1) + 
			geom_point(data=subset(tmp, variable=='from data'), size=2) + geom_point(data=subset(tmp, variable=='with largest\nmaximum likelihood\nacross bootstrap replicates'), size=2) +
			geom_point(data=subset(tmp, variable=='median\npatristic distance'), size=2) +
			scale_y_continuous(expand = c(0.01, 0.001)) +
			scale_colour_manual(values=c("#7FC97F",'black', "#BEAED4", "#FDC086")) +
			labs(colour='maximum\nlikelihood\ntree', y='estimated patristic distance', x='sequence pairs,\nfirst sequence from recipient MSM and second sequence from potential transmitter') +
			theme_bw() +
			theme(axis.text.x=element_blank(), axis.ticks=element_blank(),legend.key.size=unit(11,'mm'))
	file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tmp/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlout_Wed_Dec_18_11:37:00_2013_distPairs.pdf'
	ggsave(file=file, h=8, w=16)
	
}
######################################################################################
project.PANGEA.data.preprocess<- function()
{
	#	this cleans and pre-processes the PANGEA metadata
	#	in PANGEA.HIV.sim/misc/treecomparison, I merge this with other study data etc
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_metadata/original_161128/MetaData_161219.csv'
	dfp		<- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	setnames(dfp, colnames(dfp), toupper(colnames(dfp)))
	#	clean up trailing/leading whitespace
	for(x in colnames(dfp))
		if(any(class(x)=='character'))
			set(dfp, NULL, x, gsub(' +$','',gsub('^ +','',dfp[[x]])))
	#	dates to numeric
	set(dfp, NULL, 'SAMPLE_DATE', dfp[, as.Date(SAMPLE_DATE, format="%d/%m/%Y")])
	set(dfp, NULL, 'VL_DATE', dfp[, as.Date(VL_DATE, format="%d/%m/%Y")])
	set(dfp, NULL, 'CD4_DATE', dfp[, as.Date(CD4_DATE, format="%d/%m/%Y")])
	for(x in colnames(dfp))
		if(class(dfp[[x]])=='Date')
			set(dfp, NULL, x, hivc.db.Date2numeric(dfp[[x]]))
	set(dfp, NULL, 'DOB_YEAR', dfp[, as.numeric(DOB_YEAR)])
	#	NA cohorts
	set(dfp, dfp[, which(COHORT_ID=='')], 'COHORT_ID', NA_character_)
	#	handle CD4 ranges VL ranges
	set(dfp, NULL, 'VL_RANGE', dfp[, gsub('to ','',gsub('> ','',VL_RANGE))])
	set(dfp, dfp[, which(VL_RANGE=='<50')], 'VL_RANGE', '0 49')
	set(dfp, dfp[, which(VL_RANGE=='>400 1000')], 'VL_RANGE', '401 1000')
	set(dfp, dfp[, which(VL_RANGE=='>1000 5000')], 'VL_RANGE', '1001 5000')
	set(dfp, dfp[, which(VL_RANGE=='>5000 10000')], 'VL_RANGE', '5001 10000')
	set(dfp, dfp[, which(VL_RANGE=='>50000 100000')], 'VL_RANGE', '50001 100000')
	set(dfp, dfp[, which(VL_RANGE=='>100000 500000')], 'VL_RANGE', '100001 500000')
	set(dfp, dfp[, which(VL_RANGE=='>10000 50000')], 'VL_RANGE', '10001 50000')
	set(dfp, dfp[, which(VL_RANGE=='>500000 1000000')], 'VL_RANGE', '500001 1000000')
	set(dfp, dfp[, which(VL_RANGE=='1000000')], 'VL_RANGE', '1000000 Inf')	
	dfp[, VL_U:=NA_character_]
	dfp[, VL_L:=NA_character_]
	tmp		<- dfp[, which(VL_RANGE!='')]
	set(dfp, tmp, 'VL_L', dfp[tmp, gsub(' [0-9]+$| Inf+$','',VL_RANGE)])
	set(dfp, tmp, 'VL_U', dfp[tmp, gsub('^[0-9]+ ','',VL_RANGE)])
	set(dfp, NULL, 'VL_L', dfp[, as.numeric(VL_L)])
	set(dfp, NULL, 'VL_U', dfp[, as.numeric(VL_U)])
	set(dfp, NULL, 'CD4_RANGE', dfp[, gsub('to ','',CD4_RANGE)])	
	set(dfp, dfp[, which(CD4_RANGE=='350 <500')], 'CD4_RANGE', '350 499')
	set(dfp, dfp[, which(CD4_RANGE=='200 <350')], 'CD4_RANGE', '200 349')
	set(dfp, dfp[, which(CD4_RANGE=='100 <200')], 'CD4_RANGE', '100 199')
	set(dfp, dfp[, which(CD4_RANGE=='<100')], 'CD4_RANGE', '0 99')
	set(dfp, dfp[, which(CD4_RANGE=='>1600')], 'CD4_RANGE', '1601 Inf')
	dfp[, CD4_U:=NA_character_]
	dfp[, CD4_L:=NA_character_]
	tmp		<- dfp[, which(CD4_RANGE!='')]
	set(dfp, tmp, 'CD4_L', dfp[tmp, gsub(' [0-9]+$| Inf+$','',CD4_RANGE)])
	set(dfp, tmp, 'CD4_U', dfp[tmp, gsub('^[0-9]+ ','',CD4_RANGE)])
	set(dfp, NULL, 'CD4_L', dfp[, as.numeric(CD4_L)])
	set(dfp, NULL, 'CD4_U', dfp[, as.numeric(CD4_U)])
	#
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEA_metadata/processed_metadata/PANGEA_meta_161128.rda'
	save(dfp, file=outfile)	
}
######################################################################################
project.hivc.examlclock<- function()
{
	#
	#	root to tip divergence
	#	
	require(ape)
	require(adephylo)
	require(phytools)
	indir				<- paste(DATA,"tmp",sep='/')
	infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	infiletree			<- paste(infile,"examlbs500",sep="_")
	indircov			<- paste(DATA,"derived",sep='/')
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	insignat			<- "Wed_Dec_18_11:37:00_2013"		
	file				<- paste(indir,'/',infiletree,'_',gsub('/',':',insignat),".R",sep='')
	load(file)
	#
	file	<- paste( indir, '/', infiletree, '_', insignat,'_ExaMLTree.pdf', sep='' )
	pdf(file=file, h=150, w=10)
	plot(ph, cex=0.7)
	dev.off()
	#	visual inspection: H0576-19 is 'outgoup' to Dutch seqs (and not in the Dutch database)
	tmp		<- which(ph$tip.label=="H0576-19")
	ph		<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
	ph		<- ladderize(ph)	
	tmp		<- node.depth.edgelength(ph)
	df		<- data.table(FASTASampleCode= ph$tip.label, height=tmp[seq_len(Ntip(ph))])
	#
	file				<- paste(indircov,'/',infilecov,".R",sep='')
	load(file)
	df			<- merge( subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyPos_T1, AnyT_T1)), df, by='FASTASampleCode' )
	set(df, NULL, 'PosSeqT', hivc.db.Date2numeric(df[,PosSeqT]))
	set(df, NULL, 'AnyPos_T1', hivc.db.Date2numeric(df[,AnyPos_T1]))
	set(df, NULL, 'AnyT_T1', hivc.db.Date2numeric(df[,AnyT_T1]))
	df			<- subset(df, !is.na(PosSeqT))
	ggplot(df, aes(x=PosSeqT, y=height)) + geom_point()
	
	df.clock	<- lm(height~PosSeqT, data=as.data.frame(df))
	summary(df.clock)
	#Adjusted R-squared: 0.1056
	#             Estimate Std. Error t value Pr(>|t|)    
	#(Intercept) -5.690e+00  1.865e-01  -30.51   <2e-16 ***
	#PosSeqT      2.917e-03  9.297e-05   31.37   <2e-16 ***
	dfbh		<- df[, {
						z<- which.min(PosSeqT)
						list(FASTASampleCode=FASTASampleCode[z], PosSeqT= PosSeqT[z], AnyT_T1=AnyT_T1[z], height=height[z])
					}, by='Patient']
	dfbh.clock	<- lm(height~PosSeqT, data=as.data.frame(dfbh))
	#dfbh.clock	<- lm(height~PosSeqT, data=as.data.frame(subset(dfbh, PosSeqT>1996.6)))
	summary(dfbh.clock)		
	#Adjusted R-squared: 0.1271	    
	#(Intercept) -6.3631412  0.2165435  -29.39   <2e-16 ***
	#PosSeqT      0.0032520  0.0001079   30.13   <2e-16 ***
	dfbh[, b4T:= factor(is.na(AnyT_T1) | PosSeqT<AnyT_T1, labels=c('No','Yes'))]
	dfb4T.clock	<- lm(height~PosSeqT, data=as.data.frame(subset(dfbh, b4T)))
	summary(dfb4T.clock)
	#Adjusted R-squared: 0.105
	#(Intercept) -6.702884   0.303098  -22.11   <2e-16 ***
	#PosSeqT      0.003422
	ggplot(dfbh, aes(x=PosSeqT, y=height, colour=b4T)) + geom_point(alpha=0.75) +
		scale_x_continuous(breaks=seq(1980,2020,2)) +
		scale_colour_brewer(name='Sampled before ART start', palette='Set1') +
		stat_smooth(method="lm", colour='black') + stat_smooth(method='lm', data=subset(dfbh, b4T=='Yes'), colour="black", linetype=2) +
		labs(x='Sequence sampling date', y='root-to-tip divergence') +
		theme(legend.position=c(0,1), legend.justification=c(0,1))
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_clock.pdf'
	
	ggplot(dfbh, aes(x=PosSeqT, y=height, colour=b4T, pch=b4T)) + geom_point(alpha=0.75) +
			scale_x_continuous(breaks=seq(1980,2020,2)) +
			scale_colour_manual(values=c('red','grey50')) +
			stat_smooth(method="lm", colour='black', aes(pch='Yes')) + 
			labs(x='Sequence sampling date', y='root-to-tip divergence', colour='Sampled before ART start', pch='Sampled before ART start') +
			theme_bw() +
			theme(legend.position=c(0,1), legend.justification=c(0,1))	
	file	<- '/Users/Oliver/duke/2014_HIVNL_monitoringreport/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_clock.pdf'
	ggsave(file=file, w=15, h=5)
	#
	#	repeat R2 for all bootstrap trees, just to check
	#
	indir			<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tmp/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlout_Wed_Dec_18_11:37:00_2013'
	infiles			<- list.files(indir)
	infiles			<- infiles[ grepl('^ExaML_result.*finaltree*',infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')	
	df.r2t			<- lapply(seq_along(infiles), function(j)
			{					
				infile	<- infiles[j]
				cat(paste('\nprocess file', infile))
				file	<- paste(indir, infile, sep='/')
				ph		<- read.tree(file)
				tmp		<- which(ph$tip.label=="H0576-19")
				ph		<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
				ph		<- ladderize(ph)
				tmp		<- node.depth.edgelength(ph)
				df		<- data.table(FASTASampleCode= ph$tip.label, height=tmp[seq_len(Ntip(ph))])
				df		<- merge( subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyPos_T1, AnyT_T1)), df, by='FASTASampleCode' )
				set(df, NULL, 'PosSeqT', hivc.db.Date2numeric(df[,PosSeqT]))
				set(df, NULL, 'AnyT_T1', hivc.db.Date2numeric(df[,AnyT_T1]))
				df		<- subset(df, !is.na(PosSeqT))
				df		<- df[, 	{
										z<- which.min(PosSeqT)
										list(FASTASampleCode=FASTASampleCode[z], PosSeqT= PosSeqT[z], AnyT_T1=AnyT_T1[z], height=height[z])
									}, by='Patient']
				set(df, NULL, 'BS', as.numeric(substr(infile, nchar(infile)-2, nchar(infile))) )
				df
			})
	df.r2t		<- do.call('rbind',df.r2t)
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_R2T.R'
	save(df.r2t, file=file)
	df.r2t.su	<- df.r2t[,	{
								tmp		<- lm(height ~ PosSeqT)		 					 	
								list( R2=round(summary(tmp)$r.squared,d=3) )				
							},by='BS']						
	ggplot( df.r2t.su, aes(x=R2)) + geom_histogram(binwidth=0.05)
	#> summary(df.r2t.su[, R2])
	#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
	#0.0570  0.1020  0.1140  0.1136  0.1280  0.1740

	#
	#	within host divergence
	#
	file				<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_nbrlraw3kaH'				
	tmp					<- project.athena.Fisheretal.Y.rawbrl(NULL, NULL, NULL, NULL, NULL, NULL, df.tpairs.tptn=NULL, save.file=paste(file, '.R', sep=''), resume=1, plot.file=NA, method.restrictTPtoRI=0)
	Y.rawbrl			<- tmp$tpairs
	Y.rawbrl.linked		<- tmp$linked
	Y.rawbrl.unlinked	<- tmp$unlinked
	#	compute sequence sampling times	to see if time between seq sampling times could be useful
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )	
	tmp					<- merge( data.table(FASTASampleCode=Y.rawbrl.linked[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, tmp, by='t.FASTASampleCode')
	#	order with respect to PosSeqT
	tmp					<- subset(Y.rawbrl.linked, t.PosSeqT>PosSeqT)
	setnames(tmp, c("t.FASTASampleCode","t.PosSeqT","t.AnyT_T1"), c("a.FASTASampleCode","a.PosSeqT","a.AnyT_T1"))
	setnames(tmp, c("FASTASampleCode","PosSeqT","AnyT_T1"), c("t.FASTASampleCode","t.PosSeqT","t.AnyT_T1"))
	setnames(tmp, c("a.FASTASampleCode","a.PosSeqT","a.AnyT_T1"), c("FASTASampleCode","PosSeqT","AnyT_T1"))
	Y.rawbrl.linked		<- rbind( subset(Y.rawbrl.linked, t.PosSeqT<=PosSeqT), tmp, use.names=TRUE )	
	set(Y.rawbrl.linked, NULL, 'PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,PosSeqT]))
	set(Y.rawbrl.linked, NULL, 't.PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,t.PosSeqT]))
	set(Y.rawbrl.linked, NULL, 'AnyT_T1', hivc.db.Date2numeric(Y.rawbrl.linked[,AnyT_T1]))
	set(Y.rawbrl.linked, NULL, 't.AnyT_T1', hivc.db.Date2numeric(Y.rawbrl.linked[,t.AnyT_T1]))	
	Y.rawbrl.linked[, b4T:= 'both']
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT<=AnyT_T1)], 'b4T', 'only.RI')
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT<=t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'only.T')
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'none')
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, data.table(b4T= Y.rawbrl.linked[, unique(b4T)], col=c('green','red','orange')), by='b4T' )	
	set(Y.rawbrl.linked, NULL, 'b4T', Y.rawbrl.linked[, factor(b4T)])
	Y.rawbrl.linked[, dt:= abs(PosSeqT-t.PosSeqT)]		
	Y.rawbrl.linked[, brlr:= brl/dt]
	#	add treatment complications: any of ART.P ART.F ART.I ART.A yes before PosSeqT
	if(1)
	{
		tmp				<- subset(df.treatment, select=c(Patient, StopTime, TrI, TrCh.failure, TrVL.failure))
		set(tmp, NULL, 'StopTime', hivc.db.Date2numeric(tmp[,StopTime]))
		tmp				<- merge(tmp, subset( Y.rawbrl.linked, select=c(Patient, PosSeqT) ), by='Patient', allow.cartesian=TRUE)
		tmp				<- subset(tmp, StopTime<=PosSeqT)	
		tmp				<- tmp[, list( ART.C=any( sapply(.SD, function(z) any(z=='Yes', na.rm=TRUE)) ) ), by=c('Patient','PosSeqT'), .SDcols=c('TrI','TrCh.failure','TrVL.failure')]
		Y.rawbrl.linked	<- merge( Y.rawbrl.linked, tmp, by=c('Patient','PosSeqT'), all.x=TRUE )		
		set(Y.rawbrl.linked, Y.rawbrl.linked[, which(is.na(ART.C))], 'ART.C', FALSE)		
		Y.rawbrl.linked[, table(b4T, ART.C)]
		Y.rawbrl.linked[, b4T2:=b4T]		
		tmp				<- Y.rawbrl.linked[, which(b4T!='both' & ART.C==FALSE)]
		set(Y.rawbrl.linked, tmp, 'b4T2', Y.rawbrl.linked[tmp, paste(b4T2,'ART.C.No',sep='.')])
		tmp				<- Y.rawbrl.linked[, which(b4T!='both' & ART.C==TRUE)]
		set(Y.rawbrl.linked, tmp, 'b4T2', Y.rawbrl.linked[tmp, paste(b4T2,'ART.C.yes',sep='.')])
		set(Y.rawbrl.linked, Y.rawbrl.linked[, which(b4T2=='none.ART.C.yes' | b4T2=='only.T.ART.C.yes')], 'b4T2', 'ART.C.yes')	
		#Y.rawbrl.linked	<- subset(Y.rawbrl.linked, b4T2!='none.ART.C.No' | (b4T2=='none.ART.C.No' & dt<6) )
		#	labels
		tmp				<- c('both sampling dates\nbefore cART initiation','sampling dates before and\nafter cART initiation,\nno interruption or failure','both sampling dates\nafter cART initiation,\nno interruption or failure','at least one sampling date\nafter cART initiation,\nwith interruption or failure')
		Y.rawbrl.linked	<- merge( Y.rawbrl.linked, data.table(b4T2= c('both','only.T.ART.C.No','none.ART.C.No','ART.C.yes'), b4T2.long=tmp), by='b4T2' )
		set(Y.rawbrl.linked, NULL, 'b4T2.long', Y.rawbrl.linked[, factor(b4T2.long, levels=tmp)])			
		#	explore clock
		Y.rawbrl.linked[, excluded:= factor( dt>brl.linked.max.dt | dt<brl.linked.min.dt | brlr>brl.linked.max.brlr, levels=c(FALSE,TRUE), labels=c('No','Yes'))]
		ggplot(Y.rawbrl.linked, aes(x=dt, y=brl, colour=b4T2.long)) + geom_point(size=1, aes(shape=excluded)) + facet_grid(. ~ b4T2.long) +			 
			scale_x_continuous(limits=c(0,10)) + scale_y_continuous(limits=c(0,0.35), breaks=seq(0,0.35,0.02)) + scale_colour_brewer(palette='Set1', guide=FALSE) + 					
			stat_smooth(aes(colour=NULL), colour='black', data=subset(Y.rawbrl.linked, excluded=='No')) +			
			labs(x="years between within-host sampling dates", y='within-host divergence')
	
	}
	if(1)
	{
		tmp				<- c('both sampling dates\nbefore ART start','sampling dates before and\nafter ART start','both sampling dates\nafter ART start') 
		Y.rawbrl.linked	<- merge( Y.rawbrl.linked, data.table(b4T= c('both','only.T','none'), b4T.long=factor( tmp, levels=tmp )), by='b4T' )
		set(Y.rawbrl.linked, NULL, 'b4T.long', Y.rawbrl.linked[, factor(b4T.long)])		
	}	
	
	setkey(Y.rawbrl.linked, b4T.long, dt)
	#
	#	explore clock
	#
	brl.linked.max.dt= 10; brl.linked.min.dt= 1; brl.linked.max.brlr=0.02
	
	tmp		<- subset(Y.rawbrl.linked, b4T=='both' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T.long))	
	tmpg.ZAGA.B	<- gamlss(as.formula('brl ~ bs(dt, degree=4)'), sigma.formula=as.formula('~ bs(dt, degree=4)'), nu.formula=as.formula('~ bs(dt, degree=4)'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
	tmp[, y.b:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)]
	tmp[, y.u:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.B, type='response', se.fit=TRUE)$se.fit]
	tmp[, y.l:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.B, type='response', se.fit=TRUE)$se.fit]
	ans		<- copy(tmp)
	tmp		<- subset(Y.rawbrl.linked, b4T=='only.T' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T.long))
	tmpg.ZAGA.O	<- gamlss(as.formula('brl ~ bs(dt, degree=4)'), sigma.formula=as.formula('~ bs(dt, degree=4)'), nu.formula=as.formula('~ bs(dt, degree=4)'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
	tmp[, y.b:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)]
	tmp[, y.u:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.O, type='response', se.fit=TRUE)$se.fit]
	tmp[, y.l:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.O, type='response', se.fit=TRUE)$se.fit]
	ans		<- rbind(ans, tmp)
	tmp		<- subset(Y.rawbrl.linked, b4T=='none' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T.long))
	tmpg.ZAGA.N	<- gamlss(as.formula('brl ~ bs(dt, degree=4)'), sigma.formula=as.formula('~ bs(dt, degree=4)'), nu.formula=as.formula('~ bs(dt, degree=4)'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
	tmp[, y.b:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)]
	tmp[, y.u:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.N, type='response', se.fit=TRUE)$se.fit]
	tmp[, y.l:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.N, type='response', se.fit=TRUE)$se.fit]
	ans		<- rbind(ans, tmp)
	Y.rawbrl.linked[, excluded:= factor( dt>brl.linked.max.dt | dt<brl.linked.min.dt | brlr>brl.linked.max.brlr, levels=c(FALSE,TRUE), labels=c('No','Yes'))]
	ggplot(Y.rawbrl.linked, aes(x=dt, y=brl, colour=b4T.long)) + 
			geom_point(size=1, aes(shape=excluded)) + scale_x_continuous(limits=c(0,10)) + scale_y_continuous(limits=c(0,0.35)) + scale_colour_brewer(palette='Set1', guide=FALSE) + facet_grid(. ~ b4T.long) +		
			#geom_line(aes(x=dt, y=y.b), colour='black', data=ans) +
			stat_smooth(aes(colour=NULL), colour='black', data=subset(Y.rawbrl.linked, excluded=='No')) +
			#geom_ribbon(aes(x=dt, ymin=y.l, ymax=y.u), alpha=0.2, data=ans, inherit.aes = FALSE) +
			labs(x="years between within-host sampling dates", y='within-host divergence')
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examl_clockwhexplore.pdf'
	ggsave(file=file, w=8, h=8)	
	#
	#	mean evol rate
	#
	tmp			<- subset(Y.rawbrl.linked, b4T=='both' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T.long))		
	tmpg.ZAGA.B	<- gamlss(as.formula('brl ~ dt-1'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
	tmp[, y.b:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)]
	tmp[, y.u:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.B, type='response', se.fit=TRUE)$se.fit]
	tmp[, y.l:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.B, type='response', se.fit=TRUE)$se.fit]
	ans			<- copy(tmp)
	tmp			<- subset(Y.rawbrl.linked, b4T=='only.T' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T.long))
	tmpg.ZAGA.O	<- gamlss(as.formula('brl ~ dt-1'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
	tmp[, y.b:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)]
	tmp[, y.u:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.O, type='response', se.fit=TRUE)$se.fit]
	tmp[, y.l:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.O, type='response', se.fit=TRUE)$se.fit]
	ans			<- rbind(ans, tmp)
	tmp			<- subset(Y.rawbrl.linked, b4T=='none' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T.long))
	tmpg.ZAGA.N	<- gamlss(as.formula('brl ~ dt-1'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
	tmp[, y.b:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)]
	tmp[, y.u:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.N, type='response', se.fit=TRUE)$se.fit]
	tmp[, y.l:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.N, type='response', se.fit=TRUE)$se.fit]
	ans			<- rbind(ans, tmp)
	#	
	Y.rawbrl.linked[, excluded:= factor( dt>brl.linked.max.dt | dt<brl.linked.min.dt | brlr>brl.linked.max.brlr, levels=c(FALSE,TRUE), labels=c('No','Yes'))]
	ggplot(Y.rawbrl.linked, aes(x=dt, y=brl, colour=b4T.long)) + 
			geom_point(size=1, aes(shape=excluded)) + scale_x_continuous(limits=c(0,10)) + scale_y_continuous(limits=c(0,0.35)) + scale_colour_brewer(palette='Set1', guide=FALSE) + facet_grid(. ~ b4T.long) +					
			stat_smooth(aes(colour=NULL), colour='black', data=subset(Y.rawbrl.linked, excluded=='No')) +
			geom_line(aes(x=dt, y=y.b), colour='black', linetype='dotdash', data=ans) +	
			#geom_ribbon(aes(x=dt, ymin=y.l, ymax=y.u), alpha=0.2, data=ans, inherit.aes = FALSE) +
			labs(x="years between within-host sampling dates", y='within-host divergence')  
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examl_clockwh.pdf'
	ggsave(file=file, w=8, h=8)
	if(1)
	{
		#
		#	mean evol rate
		#
		tmp			<- subset(Y.rawbrl.linked, b4T2=='both' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T2.long))		
		tmpg.ZAGA.B	<- gamlss(as.formula('brl ~ dt'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
		tmp[, y.b:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)]
		tmp[, y.u:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.B, type='response', se.fit=TRUE)$se.fit]
		tmp[, y.l:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.B, type='response', se.fit=TRUE)$se.fit]
		ans			<- copy(tmp)
		tmp			<- subset(Y.rawbrl.linked, b4T2=='only.T.ART.C.No' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T2.long))
		tmpg.ZAGA.O	<- gamlss(as.formula('brl ~ dt'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
		tmp[, y.b:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)]
		tmp[, y.u:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.O, type='response', se.fit=TRUE)$se.fit]
		tmp[, y.l:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.O, type='response', se.fit=TRUE)$se.fit]
		ans			<- rbind(ans, tmp)		
		tmp			<- subset(Y.rawbrl.linked, b4T2=='none.ART.C.No' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T2.long))
		tmpg.ZAGA.N	<- gamlss(as.formula('brl ~ dt'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
		tmp[, y.b:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)]
		tmp[, y.u:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.N, type='response', se.fit=TRUE)$se.fit]
		tmp[, y.l:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.N, type='response', se.fit=TRUE)$se.fit]
		ans			<- rbind(ans, tmp)		
		tmp			<- subset(Y.rawbrl.linked, b4T2=='ART.C.yes' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T2.long))
		tmpg.ZAGA.C	<- gamlss(as.formula('brl ~ dt'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
		tmp[, y.b:=predict(tmpg.ZAGA.C, type='response', se.fit=FALSE)]
		tmp[, y.u:=predict(tmpg.ZAGA.C, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.C, type='response', se.fit=TRUE)$se.fit]
		tmp[, y.l:=predict(tmpg.ZAGA.C, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.C, type='response', se.fit=TRUE)$se.fit]
		ans			<- rbind(ans, tmp)
		#	
		Y.rawbrl.linked[, excluded:= factor( dt>brl.linked.max.dt | dt<brl.linked.min.dt | brlr>brl.linked.max.brlr, levels=c(FALSE,TRUE), labels=c('No','Yes'))]
		ggplot(Y.rawbrl.linked, aes(x=dt, y=brl, colour=b4T2.long)) + 
				geom_point(size=1, aes(shape=excluded)) + facet_grid(. ~ b4T2.long) +
				scale_x_continuous(limits=c(0,10)) + scale_y_continuous(limits=c(0,0.35), breaks=seq(0,0.35,0.02)) + scale_colour_brewer(palette='Set2', guide=FALSE) +
				stat_smooth(aes(colour=NULL), colour='black', data=subset(Y.rawbrl.linked, excluded=='No')) +
				geom_line(aes(x=dt, y=y.b), colour='black', linetype='dotdash', data=ans) +					
				labs(x="years between within-host sampling dates", y='within-host divergence')  
		file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examl_clockwh_ARTC.pdf'
		ggsave(file=file, w=12, h=8)
		#
		sapply( list(tmpg.ZAGA.B, tmpg.ZAGA.O, tmpg.ZAGA.N, tmpg.ZAGA.C), Rsq )
		#0.23127742  0.04869653 0.10655731  0.14558462
		#	excluded	outside 1<dt<10 & brlr>0.02
		print(tmpg.ZAGA.B)	#dt		0.002663		#sigma	0.9496      -0.1229
		print(tmpg.ZAGA.O)	#dt		0.003733		#sigma	0.71293     -0.03297
		print(tmpg.ZAGA.N)	#dt  	0.005023		#sigma	0.50178      0.05858
		print(tmpg.ZAGA.C)	#dt  	0.00621			#sigma	0.56372     -0.06687 
	}
	
	
	sapply( list(tmpg.ZAGA.B, tmpg.ZAGA.O, tmpg.ZAGA.N), Rsq )
	#0.4411087 0.2959325 0.4299490
	#	excluded	only outside 1<dt<10
	print(tmpg.ZAGA.B)	#		0.001567
	print(tmpg.ZAGA.O)	#dt		0.002675
	print(tmpg.ZAGA.N)	#dt  	0.003985
	#	excluded	outside 1<dt<10 & brlr>0.01
	print(tmpg.ZAGA.B)	#		0.002177
	print(tmpg.ZAGA.O)	#dt		0.002849
	print(tmpg.ZAGA.N)	#dt  	0.004247 
	#	excluded	outside 1<dt<10 & brlr>0.02
	print(tmpg.ZAGA.B)	#		0.002663
	print(tmpg.ZAGA.O)	#dt		0.004007
	print(tmpg.ZAGA.N)	#dt  	0.006346 
	#	excluded	outside 1<dt<10 & brlr>0.025
	print(tmpg.ZAGA.B)	#		0.00275
	print(tmpg.ZAGA.O)	#dt		0.004131  
	print(tmpg.ZAGA.N)	#dt  	0.006926 
	
	#
	#	EXCLUDE
	#	
	Y.rawbrl.linked		<- subset( Y.rawbrl.linked, dt!=0 )
	cat(paste('\nY.rawbrl.linked all, n=',nrow(Y.rawbrl.linked)))
	Y.rawbrl.linked		<- subset(Y.rawbrl.linked,  brlr<=brl.linked.max.brlr)
	cat(paste('\nY.rawbrl.linked max brlr ',brl.linked.max.brlr,' excluded, n=',nrow(Y.rawbrl.linked)))
	Y.rawbrl.linked		<- subset(Y.rawbrl.linked,  dt<=brl.linked.max.dt)		
	cat(paste('\nY.rawbrl.linked max dt ',brl.linked.max.dt,' excluded, n=',nrow(Y.rawbrl.linked)))	
	Y.rawbrl.linked		<- subset(Y.rawbrl.linked,  dt>brl.linked.min.dt)		
	cat(paste('\nY.rawbrl.linked min dt ',brl.linked.min.dt,' excluded, n=',nrow(Y.rawbrl.linked)))
	Y.rawbrl.linked[, brlz:=brl]
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(brlz<1e-4)], 'brlz', 0)	
	#
	# wilcox test for same mean
	#
	test	<- list( 	both.only.T= 	wilcox.test( subset(Y.rawbrl.linked, b4T=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.T')[, log10(brl)] ),				
						both.none= 		wilcox.test( subset(Y.rawbrl.linked, b4T=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='none')[, log10(brl)] ),						
						none.only.T= 	wilcox.test( subset(Y.rawbrl.linked, b4T=='none')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.T')[, log10(brl)] )		)
	cat(paste('\nwilcox test for same mean'))
	print( sapply(test, '[[', 'p.value') )
	if(1)
	{
		test	<- list( 	both.only.T= 	wilcox.test( subset(Y.rawbrl.linked, b4T2=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T2=='only.T.ART.C.No')[, log10(brl)] ),
							both.none.CNo= 	wilcox.test( subset(Y.rawbrl.linked, b4T2=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T2=='none.ART.C.No')[, log10(brl)] ),
							both.CYes= 	wilcox.test( subset(Y.rawbrl.linked, b4T2=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T2=='ART.C.yes')[, log10(brl)] ),							
							both.none= 		wilcox.test( subset(Y.rawbrl.linked, b4T2=='only.T.ART.C.No')[, log10(brl)], subset(Y.rawbrl.linked, b4T2=='none.ART.C.No')[, log10(brl)] ),						
							none.only.T= 	wilcox.test( subset(Y.rawbrl.linked, b4T2=='none.ART.C.No')[, log10(brl)], subset(Y.rawbrl.linked, b4T2=='ART.C.yes')[, log10(brl)] )		)
		cat(paste('\nwilcox test for same mean'))
		print( sapply(test, '[[', 'p.value') )
		#	  both.only.T both.none.CNo     both.CYes     both.none   none.only.T 
 		#3.529053e-02  7.964509e-01  1.651674e-21  3.479906e-01  2.830199e-05 		
	}
	#wilcox test for same mean		excluded brl>0.02
	# 	 both.only.T    both.none  none.only.T 
	#2.071745e-04 1.279865e-22 1.322711e-17
	#								excluded brl>0.025
	# both.only.T    both.none  none.only.T 
	#2.238328e-04 6.624977e-24 4.412482e-20
	#
	# AIC between EXP GA ZAGA	-- AIC not comparable between brl and brlz
	#
	tmpd.B		<- subset(Y.rawbrl.linked, b4T=='both', select=c(brl, brlz, dt, b4T.long))
	tmpd.O		<- subset(Y.rawbrl.linked, b4T=='only.T', select=c(brl, brlz, dt, b4T.long))
	tmpd.N		<- subset(Y.rawbrl.linked, b4T=='none', select=c(brl, brlz, dt, b4T.long))
	nrow(tmpd.B)		#[1] 134	excluded <0.01			141		excluded <0.02		142		excluded <0.025
	nrow(tmpd.O)		#[1] 328							373							376
	nrow(tmpd.N)		#[1] 1576							2154						2257
	tmpg.E.B	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.B), family=EXP, n.cyc = 40)	
	tmpg.E.O	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.O), family=EXP, n.cyc = 40)
	tmpg.E.N	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.N), family=EXP, n.cyc = 40)	
	tmpg.GA.B	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.B), family=GA, n.cyc = 40)
	tmpg.GA.O	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.O), family=GA, n.cyc = 40)
	tmpg.GA.N	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.N), family=GA, n.cyc = 40)
	tmpg.ZAGA.B	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.B), family=ZAGA, n.cyc = 40)			
	tmpg.ZAGA.O	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.O), family=ZAGA, n.cyc = 40)		
	tmpg.ZAGA.N	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.N), family=ZAGA, n.cyc = 40)
	AIC(tmpg.E.B, tmpg.GA.B, tmpg.ZAGA.B, tmpg.E.O, tmpg.GA.O, tmpg.ZAGA.O, tmpg.E.N, tmpg.GA.N, tmpg.ZAGA.N)
	#            df         AIC
	#tmpg.GA.N    2 -13443.5388		tmpg.GA.O    2  -3006.6705		tmpg.GA.B    2  -1322.0507
	#tmpg.ZAGA.N  3 -13441.5387		tmpg.ZAGA.O  3  -3004.6703		tmpg.ZAGA.B  3  -1320.0505
	#tmpg.E.N     1 -12494.9765		tmpg.E.O     1  -2351.8905		tmpg.E.B     1   -925.7479		
	tmpp.ZAGA	<- sapply( list(tmpg.ZAGA.B, tmpg.ZAGA.O, tmpg.ZAGA.N), function(z){	c(mu=as.double(exp(z[['mu.coefficients']])), sigma=as.double(exp(z[['sigma.coefficients']])), nu=as.double(1/(1+exp(-z[['nu.coefficients']]))) ) 	})
	colnames(tmpp.ZAGA)	<- c('B','O','N')
	tmpp.GA	<- sapply( list(tmpg.GA.B, tmpg.GA.O, tmpg.GA.N), function(z){	c(mu=as.double(exp(z[['mu.coefficients']])), sigma=as.double(exp(z[['sigma.coefficients']])) ) 	})
	colnames(tmpp.GA)	<- c('B','O','N')	
	tmpp.E	<- as.matrix(t(sapply( list(tmpg.E.B, tmpg.E.O, tmpg.E.N), FUN=function(z){	c(mu=as.double(exp(z[['mu.coefficients']])) ) 	}, simplify=FALSE)))
	colnames(tmpp.E)	<- c('B','O','N')
	rownames(tmpp.E)	<- c('mu')	
	#
	# ks.test between EXP GA ZAGA	
	#	
	test		<- list(	E.B= ks.test(tmpd.B[,brl], pEXP, mu=unlist(tmpp.E['mu','B'])),
							E.O= ks.test(tmpd.O[,brl], pEXP, mu=unlist(tmpp.E['mu','O'])),
							E.N= ks.test(tmpd.N[,brl], pEXP, mu=unlist(tmpp.E['mu','N'])),
							GA.B= ks.test(tmpd.B[,brl], pGA, mu=unlist(tmpp.GA['mu','B']), sigma=unlist(tmpp.GA['sigma','B'])),
							GA.O= ks.test(tmpd.O[,brl], pGA, mu=unlist(tmpp.GA['mu','O']), sigma=unlist(tmpp.GA['sigma','O'])),
							GA.N= ks.test(tmpd.N[,brl], pGA, mu=unlist(tmpp.GA['mu','N']), sigma=unlist(tmpp.GA['sigma','N'])),
							ZAGA.B= ks.test(subset(tmpd.B, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','B']), sigma=unlist(tmpp.ZAGA['sigma','B']), nu=0),
							ZAGA.O= ks.test(subset(tmpd.B, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','O']), sigma=unlist(tmpp.ZAGA['sigma','O']), nu=0),
							ZAGA.N= ks.test(subset(tmpd.B, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','N']), sigma=unlist(tmpp.ZAGA['sigma','N']), nu=0)
							)	
	print( sapply(test, '[[', 'p.value') )
	#	excluded > 0.02
	#         E.B          E.O          E.N         GA.B         GA.O         GA.N       ZAGA.B       ZAGA.O       ZAGA.N 
	#0.000000e+00 0.000000e+00 0.000000e+00 5.000988e-08 1.911804e-13 0.000000e+00 4.456872e-01 1.358416e-06 4.620748e-13
	#	excluded > 0.025
	#         E.B          E.O          E.N         GA.B         GA.O         GA.N       ZAGA.B       ZAGA.O       ZAGA.N 
	#0.000000e+00 0.000000e+00 0.000000e+00 5.514628e-08 2.016165e-13 0.000000e+00 4.273192e-01 2.096416e-06 7.571721e-14 

	#0 is <2.2e-16
	#
	#	Parameters
	#
	tmpp		<- do.call('rbind',	list( 	data.table(d='ZAGA',b4T=c('B','O','N'),mu=unlist(tmpp.ZAGA['mu',]), sigma=tmpp.ZAGA['sigma',], nu=tmpp.ZAGA['nu',]),
											data.table(d='GA',b4T=c('B','O','N'),mu=unlist(tmpp.GA['mu',]), sigma=tmpp.GA['sigma',], nu=NA),
											data.table(d='EXP',b4T=c('B','O','N'),mu=unlist(tmpp.E['mu',]), sigma=NA, nu=NA)	)	)
	tmp			<- c('both sampling dates\nbefore ART start','sampling dates before and\nafter ART start','both sampling dates\nafter ART start')
	tmp			<- data.table(b4T= c('B','O','N'), b4T.long=factor( tmp, levels=tmp ))
	tmpp		<- merge( tmpp, tmp, by='b4T' )
	#	excluded > 0.02	
	#   b4T    d          mu     sigma        nu                                   b4T.long
	#1:   B ZAGA 0.013857413 0.7738736 0.3829787      both sampling dates\nbefore ART start
	#2:   B   GA 0.008553248 1.9267489        NA      both sampling dates\nbefore ART start
	#3:   B  EXP 0.008553248        NA        NA      both sampling dates\nbefore ART start
	#4:   N ZAGA 0.027019369 0.7546553 0.1429898       both sampling dates\nafter ART start
	#5:   N   GA 0.023157079 1.4154861        NA       both sampling dates\nafter ART start
	#6:   N  EXP 0.023157079        NA        NA       both sampling dates\nafter ART start
	#7:   O ZAGA 0.021351489 0.8010297 0.2975871 sampling dates before and\nafter ART start
	#8:   O   GA 0.014999959 1.8104345        NA sampling dates before and\nafter ART start
	#9:   O  EXP 0.014999959        NA        NA sampling dates before and\nafter ART start
	#	excluded > 0.025
	#b4T    d          mu     sigma        nu                                   b4T.long
	#1:   B ZAGA 0.014198191 0.7819489 0.3802817      both sampling dates\nbefore ART start
	#2:   B   GA 0.008801788 1.9263632        NA      both sampling dates\nbefore ART start
	#3:   B  EXP 0.008801788        NA        NA      both sampling dates\nbefore ART start
	#4:   N ZAGA 0.028178663 0.7533713 0.1364643       both sampling dates\nafter ART start
	#5:   N   GA 0.024334430 1.3967320        NA       both sampling dates\nafter ART start
	#6:   N  EXP 0.024334430        NA        NA       both sampling dates\nafter ART start
	#7:   O ZAGA 0.021496155 0.7984599 0.2952128 sampling dates before and\nafter ART start
	#8:   O   GA 0.015152594 1.8054715        NA sampling dates before and\nafter ART start
	#9:   O  EXP 0.015152594        NA        NA sampling dates before and\nafter ART start
	#	CDF plot
	#
	tmp			<- tmpp[, {
								if(b4T=='O')
									brl			<- tmpd.O[,brlz]
								if(b4T=='B')
									brl			<- tmpd.B[,brlz]
								if(b4T=='N')
									brl			<- tmpd.N[,brlz]	
								brl.cdf		<- data.table( brl=sort(brl), empirical=seq_along(brl)/length(brl) )[, list(empirical=tail(empirical, 1)), by='brl']
								if(d=='EXP')
									brl.cdf[, predicted:= pEXP( brl.cdf[,brl], mu=mu ) ]
								if(d=='GA')
									brl.cdf[, predicted:= pGA( brl.cdf[,brl], mu=mu, sigma=sigma ) ]
								if(d=='ZAGA')
									brl.cdf[, predicted:= pZAGA( brl.cdf[,brl], mu=mu, sigma=sigma, nu=nu ) ]
								brl.cdf				
							}, by=c('b4T.long', 'd')]
	tmp			<- melt(tmp, id=c('b4T.long','d', 'brl'))	
	ggplot( data=tmp, aes(x=brl, y=value, colour=variable)) + 
			labs(x='within-host divergence', y='c.d.f.') + 
			scale_colour_brewer(name='', palette='Paired') + 
			geom_point(data=subset(tmp, variable=='empirical'), show_guide=FALSE ) +
			geom_line() + facet_grid(d ~ b4T.long, scales='free_x', margins=FALSE) +
			theme(legend.key.size=unit(10,'mm'))
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examl_whcdf.pdf'
	ggsave(file=file, w=10,h=10)		
	#
	#
	#
	tmpd.B		<- subset(Y.rawbrl.linked, b4T2=='both', select=c(brl, brlz, dt, b4T2.long))
	tmpd.O		<- subset(Y.rawbrl.linked, b4T2=='only.T.ART.C.No', select=c(brl, brlz, dt, b4T2.long))
	tmpd.N		<- subset(Y.rawbrl.linked, b4T2=='none.ART.C.No', select=c(brl, brlz, dt, b4T2.long))
	tmpd.C		<- subset(Y.rawbrl.linked, b4T2=='ART.C.yes', select=c(brl, brlz, dt, b4T2.long))
	nrow(tmpd.B)		#141		excluded <0.02	
	nrow(tmpd.O)		#165						
	nrow(tmpd.N)		#28						
	nrow(tmpd.C)		#2334
	tmpg.E.B	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.B), family=EXP, n.cyc = 40)	
	tmpg.E.O	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.O), family=EXP, n.cyc = 40)
	tmpg.E.N	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.N), family=EXP, n.cyc = 40)
	tmpg.E.C	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.C), family=EXP, n.cyc = 40)
	tmpg.GA.B	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.B), family=GA, n.cyc = 40)
	tmpg.GA.O	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.O), family=GA, n.cyc = 40)
	tmpg.GA.N	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.N), family=GA, n.cyc = 40)
	tmpg.GA.C	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.C), family=GA, n.cyc = 40)
	tmpg.ZAGA.B	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.B), family=ZAGA, n.cyc = 40)			
	tmpg.ZAGA.O	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.O), family=ZAGA, n.cyc = 40)		
	tmpg.ZAGA.N	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.N), family=ZAGA, n.cyc = 40)
	tmpg.ZAGA.C	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.C), family=ZAGA, n.cyc = 40)
	#
	tmpp.ZAGA	<- sapply( list(tmpg.ZAGA.B, tmpg.ZAGA.O, tmpg.ZAGA.N, tmpg.ZAGA.C), function(z){	c(mu=as.double(exp(z[['mu.coefficients']])), sigma=as.double(exp(z[['sigma.coefficients']])), nu=as.double(1/(1+exp(-z[['nu.coefficients']]))) ) 	})
	colnames(tmpp.ZAGA)	<- c('B','O','N','C')
	tmpp.GA	<- sapply( list(tmpg.GA.B, tmpg.GA.O, tmpg.GA.N, tmpg.GA.C), function(z){	c(mu=as.double(exp(z[['mu.coefficients']])), sigma=as.double(exp(z[['sigma.coefficients']])) ) 	})
	colnames(tmpp.GA)	<- c('B','O','N','C')	
	tmpp.E	<- as.matrix(t(sapply( list(tmpg.E.B, tmpg.E.O, tmpg.E.N, tmpg.E.C), FUN=function(z){	c(mu=as.double(exp(z[['mu.coefficients']])) ) 	}, simplify=FALSE)))
	colnames(tmpp.E)	<- c('B','O','N','C')
	rownames(tmpp.E)	<- c('mu')	
	tmpp		<- do.call('rbind',	list( 	data.table(d='ZAGA',b4T2=c('B','O','N','C'),mu=unlist(tmpp.ZAGA['mu',]), sigma=tmpp.ZAGA['sigma',], nu=tmpp.ZAGA['nu',]),
											data.table(d='GA',b4T2=c('B','O','N','C'),mu=unlist(tmpp.GA['mu',]), sigma=tmpp.GA['sigma',], nu=NA),
											data.table(d='EXP',b4T2=c('B','O','N','C'),mu=unlist(tmpp.E['mu',]), sigma=NA, nu=NA)	)	)
	tmp			<- c('both sampling dates\nbefore cART initiation','sampling dates before and\nafter cART initiation,\nno interruption or failure','both sampling dates\nafter cART initiation,\nno interruption or failure','at least one sampling date\nafter cART initiation,\nwith interruption or failure')						
	tmp			<- data.table(b4T2= c('B','O','N','C'), b4T2.long=factor( tmp, levels=tmp, labels=tmp ))
	tmpp		<- merge( tmpp, tmp, by='b4T2' )
	#   b4T2    d          mu     sigma        nu
 	#1:    B ZAGA 0.013857413 0.7738736 0.3829787
 	#2:    B   GA 0.008553248 1.9267489        NA
 	#3:    B  EXP 0.008553248        NA        NA
 	#4:    C ZAGA 0.026909219 0.7605336 0.1516710
 	#5:    C   GA 0.022829153 1.4444899        NA
 	#6:    C  EXP 0.022829153        NA        NA
 	#7:    N ZAGA 0.013327744 0.5810444 0.3571429
 	#8:    N   GA 0.008570276 1.8545254        NA
 	#9:    N  EXP 0.008570276        NA        NA
	#10:    O ZAGA 0.017742656 0.7320898 0.3333333
	#11:    O   GA 0.011831054 1.8523769        NA
	#12:    O  EXP 0.011831054        NA        NA
	#
	# ks.test between EXP GA ZAGA	
	#	
	test		<- list(	E.B= ks.test(tmpd.B[,brl], pEXP, mu=unlist(tmpp.E['mu','B'])),
							E.O= ks.test(tmpd.O[,brl], pEXP, mu=unlist(tmpp.E['mu','O'])),
							E.N= ks.test(tmpd.N[,brl], pEXP, mu=unlist(tmpp.E['mu','N'])),
							E.C= ks.test(tmpd.C[,brl], pEXP, mu=unlist(tmpp.E['mu','C'])),
							GA.B= ks.test(tmpd.B[,brl], pGA, mu=unlist(tmpp.GA['mu','B']), sigma=unlist(tmpp.GA['sigma','B'])),
							GA.O= ks.test(tmpd.O[,brl], pGA, mu=unlist(tmpp.GA['mu','O']), sigma=unlist(tmpp.GA['sigma','O'])),
							GA.N= ks.test(tmpd.N[,brl], pGA, mu=unlist(tmpp.GA['mu','N']), sigma=unlist(tmpp.GA['sigma','N'])),
							GA.C= ks.test(tmpd.C[,brl], pGA, mu=unlist(tmpp.GA['mu','C']), sigma=unlist(tmpp.GA['sigma','C'])),
							ZAGA.B= ks.test(subset(tmpd.B, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','B']), sigma=unlist(tmpp.ZAGA['sigma','B']), nu=0),
							ZAGA.O= ks.test(subset(tmpd.O, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','O']), sigma=unlist(tmpp.ZAGA['sigma','O']), nu=0),
							ZAGA.N= ks.test(subset(tmpd.N, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','N']), sigma=unlist(tmpp.ZAGA['sigma','N']), nu=0),
							ZAGA.C= ks.test(subset(tmpd.C, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','C']), sigma=unlist(tmpp.ZAGA['sigma','C']), nu=0)		)	
	print( sapply(test, '[[', 'p.value') )
	#          E.B          E.O          E.N          E.C         GA.B         GA.O         GA.N         GA.C       ZAGA.B 
	#0.000000e+00 3.330669e-16 1.651107e-03 0.000000e+00 5.000988e-08 1.990942e-07 2.632056e-02 0.000000e+00 4.456872e-01 
    #  ZAGA.O       ZAGA.N       ZAGA.C 
	#6.337170e-01 3.623706e-01 8.584151e-03
	
	#	CDF plot
	#
	tmp			<- tmpp[, {
				if(b4T2=='O')
					brl			<- tmpd.O[,brlz]
				if(b4T2=='B')
					brl			<- tmpd.B[,brlz]
				if(b4T2=='N')
					brl			<- tmpd.N[,brlz]	
				if(b4T2=='C')
					brl			<- tmpd.C[,brlz]									
				brl.cdf		<- data.table( brl=sort(brl), empirical=seq_along(brl)/length(brl) )[, list(empirical=tail(empirical, 1)), by='brl']
				if(d=='EXP')
					brl.cdf[, predicted:= pEXP( brl.cdf[,brl], mu=mu ) ]
				if(d=='GA')
					brl.cdf[, predicted:= pGA( brl.cdf[,brl], mu=mu, sigma=sigma ) ]
				if(d=='ZAGA')
					brl.cdf[, predicted:= pZAGA( brl.cdf[,brl], mu=mu, sigma=sigma, nu=nu ) ]
				brl.cdf				
			}, by=c('b4T2.long', 'd')]
	tmp			<- melt(tmp, id=c('b4T2.long','d', 'brl'))	
	ggplot( data=tmp, aes(x=brl, y=value, colour=variable)) + 
			labs(x='within-host divergence', y='c.d.f.') + 
			scale_colour_brewer(name='', palette='Paired') + 
			geom_point(data=subset(tmp, variable=='empirical'), show_guide=FALSE ) +
			geom_line() + facet_grid(d ~ b4T2.long, scales='free_x', margins=FALSE) +
			theme(legend.key.size=unit(10,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5))
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examl_whcdf_ARTC.pdf'
	ggsave(file=file, w=12,h=10)		
	
	#
	#	QQ plot -- TODO can call stat_qq, need qZAGA etc and params only, and slope + int
	#
	tmp			<- tmpp[, {
				if(b4T=='O')
					brl			<- tmpd.O[,brlz]
				if(b4T=='B')
					brl			<- tmpd.B[,brlz]
				if(b4T=='N')
					brl			<- tmpd.N[,brlz]	
				if(b4T=='C')
					brl			<- tmpd.C[,brlz]					
				brl.cdf		<- data.table( brl=sort(brl), empirical=seq_along(brl)/length(brl) )[, list(empirical=tail(empirical, 1)), by='brl']
				if(d=='EXP')
				{
					brl.cdf[, predicted:= qEXP( brl.cdf[,brl], mu=mu ) ]
					qt 	<- qEXP(c(.25, .75), mu=mu )
				}
				if(d=='GA')
				{
					brl.cdf[, predicted:= qGA( brl.cdf[,brl], mu=mu, sigma=sigma ) ]
					qt 	<- qGA(c(.25, .75), mu=mu, sigma=sigma )
				}
				if(d=='ZAGA')
				{
					brl.cdf[, predicted:= qZAGA( brl.cdf[,brl], mu=mu, sigma=sigma, nu=nu ) ]										
					qt 	<- qZAGA(c(.25, .75), mu=mu, sigma=sigma, nu=nu )											
				}
				qv 		<- quantile(brl, c(.25, .75))
				slope 	<- diff(qv)/diff(qt)
				int 	<- qv[1] - slope * qt[1]
				brl.cdf				
			}, by=c('b4T.long', 'd')]
	#qplot(sample=brl, geom = "point", stat = "qq", pch='both sampling dates\nbefore cART initiation', distribution = qZAGA, dparams = list(mu=exp(ml.zaga.one[['mu.coefficients']]), sigma=exp(ml.zaga.one[['sigma.coefficients']]), nu=1/(1+exp(-ml.zaga.one[['nu.coefficients']])))) + 
	#		geom_abline(slope = slope, intercept = int) + labs(x='Zero-inflated Gamma', y='empirical') +
	#		theme(legend.justification=c(0,1), legend.position=c(0,1), legend.key.size=unit(13,'mm'), legend.title = element_blank())
	#ggsave(file=paste(substr(plot.file.one,1,nchar(plot.file.one)-4),'_zagaqq','.pdf',sep=''), w=4,h=6)		
}

######################################################################################
project.Bezemer.VLIntros.DataFile<- function()
{
	infile	<- '~/Dropbox (Infectious Disease)/2017_NL_Introductions/seq_info/NONB_flowinfo.csv'
	df		<- as.data.table(read.csv(infile))
}

######################################################################################
project.Bezemer.VLIntros.FastTrees<- function()
{	
	require(big.phylo)
	#	run FastTree on HPC
	if(0)
	{
		#indir			<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle'
		indir			<- '/work/or105/ATHENA_2016/vlintros'
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
			cmd				<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=998, hpc.q="pqeph", hpc.mem="5800mb",  hpc.nproc=1, hpc.load='module load R/3.3.2')
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
		indir			<- '/Users/Oliver/Dropbox (Infectious Disease)/2017_NL_Introductions/trees_ft'
		infiles			<- data.table(F=list.files(indir, pattern='newick$', recursive=TRUE, full.names=TRUE))		
		infiles[, {
					#F	<- '/Users/Oliver/Dropbox (Infectious Disease)/2017_NL_Introductions/trees_ft/06cpx_withD_bootstrap_trees/06cpx_withD_ft.000.newick'
					ph	<- read.newick(F)		
					tmp				<- which(grepl('subtype',ph$tip.label))[1]		
					ph				<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
					ph				<- ladderize(ph)	
					write.tree(ph, file=gsub('_ft\\.','_rr.',F))
					pdf(file=paste0(gsub('_ft\\.','_rr.',F),'.pdf'), w=20, h=10+Ntip(ph)/10)
					plot(ph, show.node.label=TRUE, cex=0.3)
					dev.off()
				}, by='F']
		
	}
}

######################################################################################
project.Bezemer.clusters<- function()
{
	#clustering
	dir.name	<- DATA  	
	indir		<- paste(dir.name,'bezemer',sep='/')
	infile		<- 'NLB10BLAST_cutRmu_examlbs200_Tue_Nov_4_2014'

	file		<- paste(indir, '/', infile, '.newick', sep='')
	ph			<- ladderize(read.tree(file=file))
	#	easy: extract bootstrap support
	#ph$node.label[2]				<- 0								#little hack so that clustering works
	ph.node.bs						<- as.numeric( ph$node.label )
	ph.node.bs[is.na(ph.node.bs)]	<- 0
	ph.node.bs						<- ph.node.bs/100
	ph$node.label					<- ph.node.bs
	#
	#	memory consuming: extract branch length statistic of subtree
	#
	dist.brl.casc					<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
	gc()		
	#	
	thresh.bs		<- 0.7
	thresh.brl		<- 0.1
	clustering.70.10<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl.casc, nodesupport=ph.node.bs,retval="all")
	
	thresh.bs		<- 0.8
	thresh.brl		<- 0.1
	clustering.80.10<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl.casc, nodesupport=ph.node.bs,retval="all")
	
	thresh.bs		<- 0.8
	thresh.brl		<- 0.08
	clustering.80.08<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl.casc, nodesupport=ph.node.bs,retval="all")
	
	thresh.bs		<- 0.8
	thresh.brl		<- 0.06
	clustering.80.06<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl.casc, nodesupport=ph.node.bs,retval="all")
	
	thresh.bs		<- 0.9
	thresh.brl		<- 0.1
	clustering.90.10<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl.casc, nodesupport=ph.node.bs,retval="all")
	
	thresh.bs		<- 0.9
	thresh.brl		<- 0.08
	clustering.90.08<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl.casc, nodesupport=ph.node.bs,retval="all")
	
	thresh.bs		<- 0.9
	thresh.brl		<- 0.06
	clustering.90.06<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl.casc, nodesupport=ph.node.bs,retval="all")
	
	thresh.bs		<- 0.95
	thresh.brl		<- 0.06
	clustering.95.06<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl.casc, nodesupport=ph.node.bs,retval="all")
	
	#	save
	file		<- paste(indir, '/', infile, '_CLUSTERING.R', sep='')
	save(file=file, ph, dist.brl.casc, clustering.70.10, clustering.80.10, clustering.80.08, clustering.80.06, clustering.90.10, clustering.90.08, clustering.90.06, clustering.95.06)
		
	tmp			<- subset( data.table( LABEL=ph$tip.label, CLUSTER=clustering.70.10$clu.mem[ seq_len(Ntip(ph))] ), !is.na(CLUSTER) )
	file		<- paste(indir, '/', infile, '_CLUSTERING_70_10.csv', sep='')
	write.csv(tmp, file=file)
	tmp			<- subset( data.table( LABEL=ph$tip.label, CLUSTER=clustering.80.10$clu.mem[ seq_len(Ntip(ph))] ), !is.na(CLUSTER) )
	file		<- paste(indir, '/', infile, '_CLUSTERING_80_10.csv', sep='')
	write.csv(tmp, file=file)
	tmp			<- subset( data.table( LABEL=ph$tip.label, CLUSTER=clustering.80.08$clu.mem[ seq_len(Ntip(ph))] ), !is.na(CLUSTER) )
	file		<- paste(indir, '/', infile, '_CLUSTERING_80_08.csv', sep='')
	write.csv(tmp, file=file)
	tmp			<- subset( data.table( LABEL=ph$tip.label, CLUSTER=clustering.80.06$clu.mem[ seq_len(Ntip(ph))] ), !is.na(CLUSTER) )
	file		<- paste(indir, '/', infile, '_CLUSTERING_80_06.csv', sep='')
	write.csv(tmp, file=file)
	tmp			<- subset( data.table( LABEL=ph$tip.label, CLUSTER=clustering.90.10$clu.mem[ seq_len(Ntip(ph))] ), !is.na(CLUSTER) )
	file		<- paste(indir, '/', infile, '_CLUSTERING_90_10.csv', sep='')
	write.csv(tmp, file=file)
	tmp			<- subset( data.table( LABEL=ph$tip.label, CLUSTER=clustering.90.08$clu.mem[ seq_len(Ntip(ph))] ), !is.na(CLUSTER) )
	file		<- paste(indir, '/', infile, '_CLUSTERING_90_08.csv', sep='')
	write.csv(tmp, file=file)
	tmp			<- subset( data.table( LABEL=ph$tip.label, CLUSTER=clustering.90.06$clu.mem[ seq_len(Ntip(ph))] ), !is.na(CLUSTER) )
	file		<- paste(indir, '/', infile, '_CLUSTERING_90_06.csv', sep='')
	write.csv(tmp, file=file)	
	tmp			<- subset( data.table( LABEL=ph$tip.label, CLUSTER=clustering.95.06$clu.mem[ seq_len(Ntip(ph))] ), !is.na(CLUSTER) )
	file		<- paste(indir, '/', infile, '_CLUSTERING_95_06.csv', sep='')
	write.csv(tmp, file=file)
}
######################################################################################
project.hivc.clustering.plottree<- function()
{
	#
	#	plot ExaML tree with 80% BS clusters
	#	
	require(ape)
	require(adephylo)
	require(phytools)
	indir				<- paste(DATA,"tmp",sep='/')
	infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	infiletree			<- paste(infile,"examlbs500",sep="_")
	indircov			<- paste(DATA,"derived",sep='/')
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	insignat			<- "Wed_Dec_18_11:37:00_2013"		
	file				<- paste(indir,'/',infiletree,'_',gsub('/',':',insignat),".R",sep='')
	load(file)
	#	visual inspection: H0576-19 is 'outgoup' to Dutch seqs (and not in the Dutch database)
	tmp					<- which(ph$tip.label=="H0576-19")
	ph					<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
	ph					<- ladderize(ph)	
	#	produce clustering
	ph.node.bs			<- as.numeric( ph$node.label )
	ph.node.bs[is.na(ph.node.bs)]	<- 0
	dist.brl			<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)	 
	thresh.bs	<- 0.8
	thresh.brl	<- 0.09
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.brl=thresh.brl, dist.brl=dist.brl, thresh.nodesupport=thresh.bs, nodesupport=ph.node.bs, retval="all")	
	#print(clustering)	
	#	produce tip labels
	df.seqinfo	<- merge(data.table(FASTASampleCode=ph$tip.label, TIPID=seq_along(ph$tip.label), FRGN='No'), subset(df.all, select=c(FASTASampleCode, Patient, isAcute)), by='FASTASampleCode', all.x=1)
	set(df.seqinfo, df.seqinfo[, which(grepl('PROT+P51',FASTASampleCode,fixed=1))], 'FRGN', 'Yes')
	set(df.seqinfo, df.seqinfo[, which(is.na(isAcute) & !is.na(Patient))], 'isAcute', 'Unconfirmed')
	df.seqinfo[, COL:= 'transparent']		
	set(df.seqinfo, df.seqinfo[, which(isAcute=='No')], 'COL', "#EF6548")
	set(df.seqinfo, df.seqinfo[, which(isAcute=='Yes')], 'COL', "#990000")
	set(df.seqinfo, df.seqinfo[, which(isAcute=='Unconfirmed')], 'COL', "#FDBB84")
	set(df.seqinfo, df.seqinfo[, which(FRGN=='Yes')], 'COL', 'grey50')
	df.seqinfo[, TXT:= '']
	set(df.seqinfo, df.seqinfo[, which(isAcute=='No')], 'TXT', 'Chr')
	set(df.seqinfo, df.seqinfo[, which(isAcute=='Yes')], 'TXT', 'Rec')
	set(df.seqinfo, df.seqinfo[, which(isAcute=='Unconfirmed')], 'TXT', '?')
	set(df.seqinfo, df.seqinfo[, which(FRGN=='Yes')], 'TXT', 'FOREIGN')
	setkey(df.seqinfo, TIPID)
	file	<- paste( indir, '/', infiletree, '_', insignat,'_ExaMLCluTree.pdf', sep='' )
	hivc.clu.plot(ph, clustering[["clu.mem"]], file=file, pdf.off=0, pdf.width=10, pdf.height=170, show.node.label= F,  cex.edge.incluster=1 )
	hivc.clu.plot.tiplabels( seq_len(Ntip(ph)), matrix(df.seqinfo[,TXT], nrow=1,ncol=nrow(df.seqinfo)), matrix(df.seqinfo[,COL], nrow=1,ncol=nrow(df.seqinfo)), cex=0.05, adj=c(-2,0.5) )
	dev.off()
	#	
	
	ph	<- "(Wulfeniopsis:0.196108,(((TN_alpinus:0.459325,TN_grandiflora:0.259364)1.00:0.313204,uniflora:1.155678)1.00:0.160549,(((TP_angustibracteata:0.054609,(TN_brevituba:0.085086,TP_stolonifera:0.086001)0.76:0.035958)1.00:0.231339,(((axillare:0.017540,liukiuense:0.018503)0.96:0.038019,stenostachyum:0.049803)1.00:0.083104,virginicum:0.073686)1.00:0.103843)1.00:0.086965,(carinthiaca:0.018150,orientalis:0.019697)1.00:0.194784)1.00:0.077110)1.00:0.199516,(((((abyssinica:0.077714,glandulosa:0.063758)1.00:0.152861,((((allionii:0.067154,(morrisonicola:0.033595,officinalis:0.067266)1.00:0.055175)1.00:0.090694,(alpina:0.051894,baumgartenii:0.024152,(bellidioides:0.016996,nutans:0.063292)0.68:0.031661,urticifolia:0.032044)0.96:0.036973,aphylla:0.117223)0.67:0.033757,(((japonensis:0.018053,miqueliana:0.033676)1.00:0.160576,vandellioides:0.099761)0.69:0.036188,montana:0.050690)1.00:0.058380)1.00:0.115874,scutellata:0.232093)0.99:0.055014)1.00:0.209754,((((((acinifolia:0.112279,reuterana:0.108698)0.94:0.055829,pusilla:0.110550)1.00:0.230282,((davisii:0.053261,serpyllifolia:0.087290)0.89:0.036820,(gentianoides:0.035798,schistosa:0.038522)0.95:0.039292)1.00:0.092830)1.00:0.169662,(((anagalloides:0.018007,scardica:0.017167)1.00:0.135357,peregrina:0.120179)1.00:0.098045,beccabunga:0.069515)1.00:0.103473)1.00:0.287909,(((((((((((agrestis:0.017079,filiformis:0.018923)0.94:0.041802,ceratocarpa:0.111521)1.00:0.072991,amoena:0.229452,(((argute_serrata:0.017952,campylopoda:0.075210)0.64:0.034411,capillipes:0.022412)0.59:0.034547,biloba:0.037143)1.00:0.141513,intercedens:0.339760,((opaca:0.019779,persica:0.035744)0.94:0.038558,polita:0.036762)1.00:0.108620,rubrifolia:0.186799)1.00:0.144789,(((bombycina_11:0.033926,bombycina_bol:0.035290,cuneifolia:0.017300,jacquinii:0.054249,oltensis:0.045755,paederotae:0.051579,turrilliana:0.017117)0.85:0.049052,czerniakowskiana:0.089983)0.93:0.051111,farinosa:0.138075)1.00:0.080565)1.00:0.104525,((albiflora:0.017984,ciliata_Anna:0.032685,vandewateri:0.017610)0.97:0.045649,arguta:0.063057,(catarractae:0.022789,decora:0.049785)0.96:0.048220,((cheesemanii:0.040125,cupressoides:0.146538)1.00:0.067761,macrantha:0.038130)1.00:0.088158,(densifolia:0.090044,formosa:0.116180)0.71:0.046353,(elliptica:0.038650,(odora:0.019325,salicornioides:0.021228)0.94:0.042950,salicifolia:0.020829)0.92:0.043978,(nivea:0.070429,(papuana:0.035003,tubata:0.031140)0.98:0.064379)0.93:0.065336,raoulii:0.109101)0.97:0.076607)0.93:0.085835,chamaepithyoides:0.485601)0.57:0.072713,(ciliata_157:0.069943,lanuginosa:0.052833)1.00:0.098638,(densiflora:0.069429,macrostemon:0.118926)0.92:0.124911,(fruticulosa:0.086891,saturejoides:0.041181)0.94:0.086148,kellererii:0.083762,lanosa:0.263033,mampodrensis:0.103384,nummularia:0.191180,pontica:0.128944,thessalica:0.129197)0.65:0.031006,(arvensis:0.342138,(((((chamaedrys:0.043720,micans:0.032021,vindobonensis:0.033309)0.51:0.034053,micrantha:0.019084)0.64:0.037906,krumovii:0.020175)1.00:0.103875,verna:0.254017)0.81:0.057105,magna:0.112657)1.00:0.104070)1.00:0.101845)1.00:0.149208,(((aznavourii:0.664103,glauca:0.405588)0.85:0.209945,praecox:0.447238)1.00:0.185614,(donii:0.260827,triphyllos:0.176032)1.00:0.194928)1.00:0.611079)0.74:0.055152,((crista:0.591702,(((cymbalaria_Avlan:0.017401,panormitana:0.017609)1.00:0.229508,((cymbalaria_Istanbul:0.028379,trichadena_332:0.016891,trichadena_Mugla:0.019131)1.00:0.196417,lycica_333:0.146772)1.00:0.097646,lycica_192:0.154877)1.00:0.234748,(((hederifolia:0.018068,triloba:0.075784)1.00:0.084865,(sibthorpioides:0.122542,sublobata:0.136951)1.00:0.074683)0.89:0.043623,stewartii:0.040679)1.00:0.596859)1.00:0.237324)0.58:0.057120,javanica:0.133802)1.00:0.137214)1.00:0.269201,(missurica:0.016685,rubra:0.019696)1.00:0.351184)0.54:0.058275)0.52:0.062485,((dahurica:0.023542,longifolia:0.016484,spicata:0.018125)0.95:0.042294,(nakaiana:0.016270,schmidtiana:0.058451)0.88:0.037207)1.00:0.261643)0.55:0.056458)1.00:0.229509,kurrooa:0.100611)0.74:0.068198,(bonarota:0.040842,lutea:0.115316)1.00:0.241657)0.99:0.085772);"
	ph <- ladderize( read.tree(text = ph) )				
	#read bootstrap support values
	thresh.bs						<- 0.9
	ph.node.bs						<- as.numeric( ph$node.label )		
	ph.node.bs[is.na(ph.node.bs)]	<- 0
	ph.node.bs[c(13,15,27,41,43)]	<- thresh.bs-0.05*seq(0.01,length.out=5)
	ph.node.bs[ph.node.bs==1]		<- seq_along(which(ph.node.bs==1))*0.005 + 0.7
	ph$node.label					<- ph.node.bs
	#read patristic distances
	stat.fun						<- hivc.clu.min.transmission.cascade
	#stat.fun						<- max
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)
	print(dist.brl)
	thresh.brl						<- quantile(dist.brl,seq(0.1,1,by=0.05))["100%"]
	print(quantile(dist.brl,seq(0.1,0.5,by=0.05)))
	print(thresh.brl)
	#produce clustering 
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	print(clustering)		
	hivc.clu.plot(ph, clustering[["clu.mem"]] )
	#produce some tip states
	ph.tip.state<- rep(1:20, each=ceiling( length( ph$tip.label )/20 ))[seq_len(Ntip(ph))]
	states		<- data.table(state.text=1:20, state.col=rainbow(20))
	states		<- states[ph.tip.state,]
	hivc.clu.plot.tiplabels( seq_len(Ntip(ph)), matrix(states[,state.text], nrow=1,ncol=nrow(states)), matrix(states[,state.col], nrow=1,ncol=nrow(states)), cex=0.4, adj=c(-1,0.5) )
	
}
######################################################################################
inla.tut1<- function()
{
	install.packages("INLA", repos="http://www.math.ntnu.no/inla/R/testing")
	library(INLA)
	data(Seeds)
	head(Seeds)
	formula	<- r ~ 1 + x1+x2+x1:x2
	#INLA f() random effect
	formula	<- r ~ 1 + x1+x2+x1:x2+f(plate, model='iid')
	result	<- inla(formula, family='binomial', data=Seeds, Ntrials=n)
	summary(result)
	plot(result)
	#	INLA parameters must be small 10-15 at max
	#	can give fixed covariate u an AR smoothing effect -> remains Gaussian
	#	AR parameter is an INLA parameter
	#	precision matrix shows correlation of response and covariates and hyperparameters
	#	ie which of the 'locations' are correlated to explain 'transmissions'	
}
######################################################################################
inla.tut.basic1<- function()
{
	library(spatstat)
	library(mvtnorm)
	library(lattice)
	library(mgcv)
	library(pixmap)
	library(numDeriv)
	library(fields)
	library(INLA)
	
	data(Seeds)
	Y		<- Seeds$r
	Z1		<- Seeds$x1
	Z2		<- Seeds$x2
	n		<- Seeds$n
	formula <- Y ~ 1 + Z1 * Z2
	data	<- data.frame(Y, Z1, Z2, n)
	result	<- inla(formula, data=data, family="binomial", Ntrials=n, verbose=TRUE, 
								control.compute=list(dic=TRUE))
	summary(result)
	plot(result)
	#	fixed effects: OK I get the posterior means etc. great.
	
	plate	<- Seeds$plate
	data	<- data.frame(Y, Z1, Z2, n, plate)
	formula <- Y ~ Z1*Z2 + f(plate, model="iid")
	result2	<- inla(formula, data=data, family="binomial", Ntrials=n, verbose=0, 
								control.compute=list(dic=TRUE))
	summary(result2)	
	
	data(trees)
	Y		<- log(trees$Volume)
	Z1		<-trees$Girth
	Z2		<-trees$Height
	prior.c	<- c(1,0.05)
	hyper	<- list(theta=list(param=prior.c))
	formula <- Y ~ 1 + f(Z1, model = "rw2", hyper=hyper) + f(Z2, model = "rw2", hyper=hyper)
	data	<- data.frame(Y, Z1, Z2)
	result3	<- inla(formula, data=data, family="normal", verbose=0, 
								control.compute=list(dic=TRUE))
	summary(result3)
	plot(result3)
	#	with random effects, the mean is always zero and the output reports marginal posterior precision
	#	Z1 linear, Z2 tips off at the end -- num issue?
	#	what is linear predictor? fitted values?
	
	prior.c	<- c(20,0.05)
	#	what are the equations behind the prior parameters?
	#	x ~ N(0, tau^-1 Q^-1) Q depends on knots
	#	with scale.model Q is re-scaled so tau is comparable across models
	#	important for prior spec, otherwise prior data may concentrate oddly
	#	scale.model ~ prior on precision * knots
	hyper	<- list(theta=list(param=prior.c))
	formula <- Y ~ 1 + f(Z1, model = "rw2", hyper=hyper, scale.model=TRUE) + f(Z2, model = "rw2", hyper=hyper, scale.model=TRUE)
	result4	<- inla(formula, data=data, family="normal", verbose=0, 
								control.compute=list(dic=TRUE))
	summary(result4)
	plot(result4)
	#	how can I 'see' the model, ie plot data and prediction on data with CIs?
	
	prior.c	<- c(1,0.05)
	hyper	<- list(theta=list(param=prior.c))
	formula <- Y ~ 1 + f(Z1, model = "rw1", hyper=hyper) + f(Z2, model = "rw1", hyper=hyper)
	data	<- data.frame(Y, Z1, Z2)
	result5	<- inla(formula, data=data, family="normal", verbose=0, 
			control.compute=list(dic=TRUE))
	summary(result5)
	plot(result5)
	#	for predict, need to: response<- c(Y, rep(n, NA))	data<- c(X, X)
	#	need to do it this way because bunch of results are thrown away

	#	inla.mcmc
	#	check demo: secret args
}
######################################################################################
inla.tut.rainfall<- function()
{
	library(lattice)
	library(INLA)
	indir	<- '~/duke/2015_tutorial_INLA/Practicals'
	file	<- paste(indir, 'rainfall.RData', sep='/')
	load(file) 
	head(rainfall)	
	loc = cbind(rainfall$lon, rainfall$lat)
		
	mesh	<- inla.mesh.create.helper(	points.domain=cbind(rainfall$lon,rainfall$lat),
									max.edge=c(2,100))
	plot(mesh)
	spde 	<- inla.spde2.matern(mesh, alpha=2)
	A 		<- inla.spde.make.A(mesh=mesh, loc =loc)
	formula <- z ~ intercept + f(nodes, model = spde) - 1
	stack 	<- inla.stack(	data= list(z=rainfall$z), 
							A= list(1,A),
							effects= list(intercept=rep(1,length(rainfall$z)),
							nodes= 1:spde$n.spde ))
	result 	<- inla(formula,
							data=inla.stack.data(stack),
							control.predictor=list(A=inla.stack.A(stack)))
	result$summary.fixed	
	result$summary.hyperpar
	#	transform marginal with delta method
	post.se <- inla.tmarginal(function(x) sqrt(1/x), result$marginals.hyperpar[[1]])
	inla.emarginal(function(x) x, post.se)
	inla.qmarginal(c(0.025, 0.5, 0.975), post.se)
	inla.hpdmarginal(0.95, post.se)
	inla.pmarginal(c(0.5, 0.7), post.se)
	
	result.field 	<- inla.spde2.result(result, 'nodes', spde, do.transf=TRUE)
	
	inla.emarginal(function(x) x, result.field$marginals.kappa[[1]]) 
	inla.emarginal(function(x) x, result.field$marginals.variance.nominal[[1]]) 
	inla.emarginal(function(x) x, result.field$marginals.range.nominal[[1]])
	
	
	require(maps) #install.packages(c("maps","maptools")) 
	require(maptools)
	## Calculate mapping between triangulation vertices and grid points:
	proj 			<- inla.mesh.projector(mesh, dims=c(200,200))
	## Construct greyscale palette function:
	my.grey.palette <- function (n,...) { return (grey.colors(n,0.05,0.95,...))}
	my.palette 		<- my.grey.palette
	## Construct map data appropriate for easy plotting:
	levelplotmap 	<- function(..., mm) {
											panel.levelplot(...)
											panel.lines(mm$x, mm$y, col="black")
										}
	#####################
	## Plot results:
	plotdata 		<- inla.mesh.project(proj, result.field$summary.values$mean)
	dev.new()
	bbb 			<- levelplot(	row.values=proj$x, column.values=proj$y, x=plotdata,								
									mm=map("usa"), panel=levelplotmap,
									col.regions=my.palette,
									xlim=range(proj$x), ylim=range(proj$y), aspect="iso",
									contour=TRUE, cuts=11, labels=FALSE, pretty=TRUE,
									xlab="Easting",ylab="Northing")

	print(bbb)
	#	geostatsp
	#	diseasemapping
}

######################################################################################
hivc.project.remove.resistancemut.save<- function()
{
	dr		<- as.data.table( read.csv( paste( CODE.HOME,"/data/IAS_primarydrugresistance_201303.csv",sep='' ), stringsAsFactors=F ) )	
	dr[,Alignment.nuc.pos:= (Gene.codon.number-1)*3+Gene.HXB2pos ]		
	dr		<- dr[,	{	tmp<- unique(Mutant); list(Mutant=tmp, Gene.codon.number=Gene.codon.number[1], Wild.type=Wild.type[1], DR.name=DR.name[1])	}, by=Alignment.nuc.pos]
	#select nucleotide codes that are consistent with drug resistance mutants
	nt2aa	<- as.data.table( read.csv( paste( CODE.HOME,"/data/standard_nt_code.csv",sep='' ), stringsAsFactors=F ) )
	setnames(nt2aa,c("AA","NTs"),c("Mutant","Mutant.NTs"))
	nt2aa	<- subset(nt2aa, select=c(Mutant,Mutant.NTs))
	dr		<- merge(dr, nt2aa, all.x=1, by="Mutant", allow.cartesian=TRUE)
	setkey(dr, "Alignment.nuc.pos")
	#print(dr, nrows=250)
	dr		<- subset(dr, select=c(Alignment.nuc.pos, Mutant.NTs, DR.name))
	set(dr, NULL, "Mutant.NTs", tolower(dr[,Mutant.NTs]))
	IAS_primarydrugresistance_201303	<- data.frame(dr)
	save(IAS_primarydrugresistance_201303, file=paste( CODE.HOME,"/data/IAS_primarydrugresistance_201303.rda",sep='' ))	
}
######################################################################################
project.Tchain.MTC.sensecheck	<- function()
{	
	#
	#	check if drug resistance mutations are masked	
	if(0)
	{
		indir	<- '/Users/Oliver/duke/2014_HIV_MTCNL'
		infile	<- 'Mother_Child_cohort.fasta'
		file	<- paste(indir, '/',infile, sep='')
		seq.pol	<- read.dna(file, format='fasta')		
		#
		file	<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
		cat(paste('\nLoading outgroup seq from file', file))
		load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
		#	align against HXB2
		seq.pol <- c(as.list(outgroup.seq.pol), as.list(seq.pol))		
		infile	<- '141105_MTCNL58_Drcheck.fasta'
		file	<- paste(indir,'/',infile, sep='')
		write.dna(file=file, seq.pol, format='fasta')
		#	hand-curate: just a bit of noise at the end
		infile	<- '141105_MTCNL58_Drcheck_HIValign2.fasta'
		file	<- paste(indir,'/',infile, sep='')
		seq.pol	<- read.dna(file=file, format='fa')		
		seq.pol	<- as.character(seq.pol	)
		seq.pol	<- as.DNAbin(seq.pol[, 169:1469])		
		# first site is 2253 in HXB2
		
		#load drug resistance mutations and select unique mutants by codon
		load( system.file(package="hivclust", "data",'IAS_primarydrugresistance_201303.rda') )		
		dr				<- as.data.table(IAS_primarydrugresistance_201303)		
		alignment.start	<- 2253	
		set(dr, NULL, "Alignment.nuc.pos", dr[,Alignment.nuc.pos]-alignment.start+1)
		#if alignment equals any of the drug resistance mutants, replace with NNN	
		seq.pol			<- seq.rm.drugresistance(as.character(seq.pol), dr, verbose=1, rtn.DNAbin=1 )	
		#rm HX1B and close gaps
		seq.pol			<- seq.pol[!grepl('HXB2',rownames(seq.pol)),]
		#	save
		infile			<- "141105_MTCNL58_nodr.fasta"
		file			<- paste(indir,'/',infile, sep='')
		write.dna(file=file, seq.pol, format='fasta')			
	}
	#	
	#	collect seq, ExaML phylos and patristic distances
	#
	if(0)
	{
		indir			<- '/Users/Oliver/duke/2014_HIV_MTCNL'		
		infile.pol		<- "141105_MTCNL58_nodr.fasta"
		file			<- paste(indir,'/',infile.pol, sep='')
		seq.pol			<- read.dna(file, format='fa')
		#
		#	run ExaML to get trees in same way
		#
		seq					<- seq.pol	
		infile				<- infile.pol
		infile.seq.sig		<- "Wed_Sep_24_12:59:06_2013"
		infile.seq			<- paste(substr(infile,1,nchar(infile)-6),'',sep='')
		file				<- paste( indir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
		save(seq, file=file)
		cmd					<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
		cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
		tmp					<- hivc.cmd.hpccaller(indir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		#
		#	get patristic distance
		#	
		file				<- list.files(indir, paste('^ExaML_result.141105_.*','finaltree.000$',sep=''), full.names=TRUE)
		ph.pol				<- read.tree(file)
		brl.pol				<- as.matrix( distTips(ph.pol , method='patristic') )
		tmp					<- as.data.table(expand.grid(TIP1=rownames(brl.pol), TIP2=colnames(brl.pol), stringsAsFactors = FALSE))
		tmp					<- tmp[, list(BRL= brl.pol[TIP1, TIP2]), by=c('TIP1','TIP2')]
		tmp					<- subset(tmp, TIP1!=TIP2)
		brl.pol				<- tmp
		#
		infile				<- '141105_MTCNL58_nodr_INFO.R'
		file				<- paste(indir, '/',infile, sep='')
		save(seq.pol, brl.pol, ph.pol, file=file)
	}	
	#	
	#	get time elapsed versus brl in transmission pairs
	#
	if(1)
	{
		# 	read metavariables
		indir	<- '/Users/Oliver/duke/2014_HIV_MTCNL'
		infile	<- 'Mother_Child_cohort_data.csv'
		file	<- paste(indir,'/',infile, sep='')
		df.ind	<- read.csv(file=file, stringsAsFactors=FALSE )
		df.ind	<- as.data.table(df.ind)
		df.ind[, PosSeqT:=df.ind[, as.Date(gsub('\\s','',sampling.date), format='%d-%m-%Y')]]
		set(df.ind, NULL, 'PosSeqT', df.ind[, hivc.db.Date2numeric(PosSeqT)])
		set(df.ind, NULL, 'patient.code', df.ind[, gsub('\\s|#','',patient.code)])
		df.ind[, PAIR:= df.ind[, sapply(strsplit(patient.code,'-'),'[[',2)]]
		df.ind[, PATIENT:= df.ind[, sapply(strsplit(patient.code,'-'),'[[',1)]]				
		df.ind[, TIME_TR:= paste('15-',month.of.birth,'-',year.of.birth,sep='')]
		set(df.ind, df.ind[, which(is.na(year.of.birth))], 'TIME_TR', NA_character_)
		set(df.ind, NULL, 'TIME_TR', df.ind[, as.Date(gsub('\\s','',TIME_TR), format='%d-%m-%Y')])
		set(df.ind, NULL, 'TIME_TR', df.ind[, hivc.db.Date2numeric(TIME_TR)])
		tmp		<- subset(df.ind, PATIENT=='mother' & (PAIR=='08' | PAIR=='17'))
		set(tmp, NULL, 'PAIR', tmp[, paste(PAIR,'b',sep='')]) 
		df.ind	<- rbind(df.ind, tmp)
		tmp		<- df.ind[, which(PATIENT=='mother' & (PAIR=='08' | PAIR=='17'))]
		set(df.ind, tmp, 'PAIR', df.ind[tmp, paste(PAIR,'a',sep='')])
		#	get transmissions
		tmp		<- subset(df.ind, PATIENT=='mother', select=c(patient.code, PosSeqT, PAIR))
		df.trm	<- subset(df.ind, PATIENT=='child', select=c(patient.code, PosSeqT, PAIR, TIME_TR))
		setnames(df.trm, c('patient.code','PosSeqT'), c('TO','TO_PosSeqT'))
		setnames(tmp, c('patient.code','PosSeqT'), c('FROM','FROM_PosSeqT'))
		df.trm	<- merge(tmp, df.trm, by='PAIR')		
		#
		#	
		#	merge branch lengths and transmissions
		infile	<- '141105_MTCNL58_nodr_INFO.R'
		file	<- paste(indir,'/',infile, sep='')
		load(file)
		setnames(brl.pol, c('TIP1','TIP2'), c('FROM','TO'))
		trm.pol.mtc	<- merge(df.trm, brl.pol, by=c('FROM','TO'))
		#	set times between sampling times etc
		trm.pol.mtc[, d_SeqT:=abs(TO_PosSeqT-FROM_PosSeqT)]
		trm.pol.mtc[, d_TSeqT:= abs(TO_PosSeqT-TIME_TR) + abs(FROM_PosSeqT-TIME_TR)]
		
		file	<- paste(indir,'/','141105_MTCNL58_nodr_INFO_TRM.R', sep='')
		save(file=file, trm.pol.mtc)
		#
		#	explore
		#
		if(1)
		{
			ggplot(trm.pol.mtc, aes(x=d_TSeqT, y=BRL)) + geom_point() + 
					labs(x='time elapsed\n(years)', y='evolutionary divergence\nbetween confirmed transmission pairs\n(nucleotide substitutions / site)') +
					#coord_cartesian(xlim=c(0,22), ylim=c(0,0.125)) +					
					scale_y_continuous(breaks=seq(0, 0.2, 0.01)) +
					scale_x_continuous(breaks=seq(0,30,2)) + theme_bw() + theme(panel.grid.minor=element_line(colour="grey75", size=0.2), panel.grid.major=element_line(colour="grey75")) 
			file	<- paste(indir, '/',"141105_MTCNL58_nodr_patristic_dTS_all.pdf", sep='')
			ggsave(file=file, w=5, h=5)	
		}
		#
		#	are there multi drug resistant sites for a particular person?
		#
		tmp		<- apply( as.character(seq.pol), 1, function(x) length(which(x=='n')))
		tmp		<- data.table(DR= tmp/3, PATIENT=names(tmp))
		setnames(tmp, c('PATIENT','DR'), c('FROM','FROM_DR'))
		trm.pol.mtc	<- merge(trm.pol.mtc, tmp, by='FROM')
		setnames(tmp, c('FROM','FROM_DR'), c('TO','TO_DR'))
		trm.pol.mtc	<- merge(trm.pol.mtc, tmp, by='TO')
		if(1)
		{
			ggplot(trm.pol.mtc, aes(x=d_TSeqT, y=BRL, colour=FROM_DR>0 | TO_DR>0)) + geom_point() + 
					labs(x='time elapsed\n(years)', y='evolutionary divergence\nbetween confirmed transmission pairs\n(nucleotide substitutions / site)') +
					#coord_cartesian(xlim=c(0,22), ylim=c(0,0.125)) +
					scale_colour_brewer(name='DR present\nin mother or child',palette='Set1') +
					scale_y_continuous(breaks=seq(0, 0.2, 0.01)) +
					scale_x_continuous(breaks=seq(0,30,2)) + theme_bw() + theme(panel.grid.minor=element_line(colour="grey75", size=0.2), panel.grid.major=element_line(colour="grey75")) 
			file	<- paste(indir, '/',"141105_MTCNL58_nodr_patristic_dTS_all_byDR.pdf", sep='')
			ggsave(file=file, w=7, h=5)	
		}
	}
}
######################################################################################
project.Tchain.Swedish.sensecheck	<- function()
{
	require(adephylo)	
	#	get readable fasta file	
	if(0)
	{
		indir	<- '/Users/Oliver/duke/2014_HIV_TChainSwedish'
		#
		#	gag
		#
		infile	<- 'set2_seqs_gag.fasta'
		file	<- paste(indir, '/',infile, sep='')				
		txt		<- scan(file = file, what = "", sep = "\n", quiet = TRUE)
		#	search for '>' names 
		tmp			<- c( which( grepl('^>',txt) ), length(txt)+1 )
		tmp			<- lapply(seq_len(length(tmp)-1), function(i)
				{				
					data.table( LABEL= gsub('>','',txt[tmp[i]]), SEQ= paste(txt[(tmp[i]+1):(tmp[i+1]-1)], sep='', collapse='') )
				})
		seq				<- do.call('rbind', tmp)	
		tmp				<- tolower(do.call('rbind',strsplit(seq[, SEQ],'')))
		rownames(tmp)	<- seq[, LABEL]
		seq				<- as.DNAbin(tmp)		
		infile			<- '141105_set2_seqs_gag.fasta'
		file			<- paste(indir,'/',infile, sep='')
		write.dna(file=file, seq, format='fasta')
		#
		#	env
		#
		infile	<- 'set2_seqs_env.fasta'
		file	<- paste(indir, '/',infile, sep='')				
		txt		<- scan(file = file, what = "", sep = "\n", quiet = TRUE)
		#	search for '>' names 
		tmp			<- c( which( grepl('^>',txt) ), length(txt)+1 )
		tmp			<- lapply(seq_len(length(tmp)-1), function(i)
				{				
					data.table( LABEL= gsub('>','',txt[tmp[i]]), SEQ= paste(txt[(tmp[i]+1):(tmp[i+1]-1)], sep='', collapse='') )
				})
		seq				<- do.call('rbind', tmp)	
		tmp				<- tolower(do.call('rbind',strsplit(seq[, SEQ],'')))
		rownames(tmp)	<- seq[, LABEL]
		seq				<- as.DNAbin(tmp)		
		infile		<- '141105_set2_seqs_env.fasta'
		file		<- paste(indir,'/',infile, sep='')
		write.dna(file=file, seq, format='fasta')		
	}
	#	
	#	collect seq, ExaML phylos and patristic distances
	#
	if(1)
	{
		#
		#	get ExaML tree
		#		
		indir				<- '/Users/Oliver/duke/2014_HIV_TChainSwedish'
		infile				<- '141105_set2_seqs_gag.fasta'
		file				<- paste(indir,'/',infile, sep='')
		seq					<- read.dna(file, format='fa')
		infile.seq.sig		<- "Wed_Sep_24_12:59:06_2013"
		infile.seq			<- paste(substr(infile,1,nchar(infile)-6),'',sep='')
		file				<- paste( indir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
		save(seq, file=file)
		cmd					<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
		cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
		tmp					<- hivc.cmd.hpccaller(indir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		
		infile				<- '141105_set2_seqs_env.fasta'
		file				<- paste(indir,'/',infile, sep='')
		seq					<- read.dna(file, format='fa')
		infile.seq.sig		<- "Wed_Sep_24_12:59:06_2013"
		infile.seq			<- paste(substr(infile,1,nchar(infile)-6),'',sep='')
		file				<- paste( indir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
		save(seq, file=file)
		cmd					<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
		cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
		tmp					<- hivc.cmd.hpccaller(indir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		
		infile				<- '141105_set2_seqs_gag.fasta'
		file				<- paste(indir,'/',infile, sep='')
		seq					<- read.dna(file, format='fa')
		seq.gag				<- seq
		infile				<- '141105_set2_seqs_env.fasta'
		file				<- paste(indir,'/',infile, sep='')
		seq.env				<- read.dna(file, format='fa')
		stopifnot( all(rownames(seq)==rownames(seq.env)) )
		seq					<- cbind(seq, seq.env)
		seq.conc			<- seq
		infile.seq.sig		<- "Wed_Sep_24_12:59:06_2013"
		infile				<- '141105_set2_seqs_conc.fasta'
		infile.seq			<- paste(substr(infile,1,nchar(infile)-6),'',sep='')
		file				<- paste( indir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
		save(seq, file=file)
		cmd					<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
		cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
		tmp					<- hivc.cmd.hpccaller(indir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		#
		#	get patristic distance
		#	
		file				<- list.files(indir, paste('^ExaML_result.141105_.*','gag','.*finaltree.000$',sep=''), full.names=TRUE)
		ph					<- read.tree(file)
		brl					<- as.matrix( distTips(ph , method='patristic') )
		tmp					<- as.data.table(expand.grid(TIP1=rownames(brl), TIP2=colnames(brl), stringsAsFactors = FALSE))
		tmp					<- tmp[, list(BRL= brl[TIP1, TIP2]), by=c('TIP1','TIP2')]
		tmp					<- subset(tmp, TIP1!=TIP2)
		tmp					<- tmp[, {
					if(TIP1<TIP2)	z<- list(TIP1n=TIP1, TIP2n=TIP2, BRL=BRL)
					if(TIP1>TIP2)	z<- list(TIP1n=TIP2, TIP2n=TIP1, BRL=BRL)
					z
				}, by=c('TIP1','TIP2')]
		set(tmp, NULL, c('TIP1','TIP2'), NULL)
		setnames(tmp, c('TIP1n','TIP2n'), c('TIP1','TIP2'))
		setkey(tmp, TIP1, TIP2)
		ph.gag				<- ph
		brl.gag				<- unique(tmp)		
		#							
		file				<- list.files(indir, paste('^ExaML_result.141105_.*','env','.*finaltree.000$',sep=''), full.names=TRUE)
		ph					<- read.tree(file)
		brl					<- as.matrix( distTips(ph , method='patristic') )
		tmp					<- as.data.table(expand.grid(TIP1=rownames(brl), TIP2=colnames(brl), stringsAsFactors = FALSE))
		tmp					<- tmp[, list(BRL= brl[TIP1, TIP2]), by=c('TIP1','TIP2')]
		tmp					<- subset(tmp, TIP1!=TIP2)
		tmp					<- tmp[, {
					if(TIP1<TIP2)	z<- list(TIP1n=TIP1, TIP2n=TIP2, BRL=BRL)
					if(TIP1>TIP2)	z<- list(TIP1n=TIP2, TIP2n=TIP1, BRL=BRL)
					z
				}, by=c('TIP1','TIP2')]
		set(tmp, NULL, c('TIP1','TIP2'), NULL)
		setnames(tmp, c('TIP1n','TIP2n'), c('TIP1','TIP2'))
		setkey(tmp, TIP1, TIP2)
		ph.env				<- ph
		brl.env				<- unique(tmp)		
		#
		file				<- list.files(indir, paste('^ExaML_result.141105_.*','conc','.*finaltree.000$',sep=''), full.names=TRUE)
		ph					<- read.tree(file)
		brl					<- as.matrix( distTips(ph , method='patristic') )
		tmp					<- as.data.table(expand.grid(TIP1=rownames(brl), TIP2=colnames(brl), stringsAsFactors = FALSE))
		tmp					<- tmp[, list(BRL= brl[TIP1, TIP2]), by=c('TIP1','TIP2')]
		tmp					<- subset(tmp, TIP1!=TIP2)
		tmp					<- tmp[, {
					if(TIP1<TIP2)	z<- list(TIP1n=TIP1, TIP2n=TIP2, BRL=BRL)
					if(TIP1>TIP2)	z<- list(TIP1n=TIP2, TIP2n=TIP1, BRL=BRL)
					z
				}, by=c('TIP1','TIP2')]
		set(tmp, NULL, c('TIP1','TIP2'), NULL)
		setnames(tmp, c('TIP1n','TIP2n'), c('TIP1','TIP2'))
		setkey(tmp, TIP1, TIP2)
		ph.conc				<- ph
		brl.conc			<- unique(tmp)		
		#
		#
		infile				<- '141105_set2_seqs_INFO.R'
		file				<- paste(indir, '/',infile, sep='')
		save(seq.env, seq.gag, seq.conc, brl.env, brl.gag, brl.conc, ph.env, ph.gag, ph.conc, file=file)
	}	
	#
	# 	get divergence between transmission pairs
	#
	if(1)
	{
		indir	<- '/Users/Oliver/duke/2014_HIV_TChainSwedish'
		infile	<- 'set2_anno.csv'
		df.s	<- as.data.table(read.csv( paste(indir, '/', infile, sep=''), stringsAsFactors=FALSE ))
		df.s	<- subset( df.s, select=c(Patient.code, PAT.id, Patient.sex, Risk.factor, Infection.country, Sampling.year, Accession, Name) )		
		df.s	<- subset( df.s, select=c(Patient.code, Sampling.year, Name) )
		setnames(df.s, c('Patient.code', 'Sampling.year', 'Name'), c('PATIENT','PosSeqT', 'LABEL'))
		set(df.s, NULL, 'PosSeqT', df.s[, PosSeqT+.5]) 
		setkey(df.s, LABEL)
		df.s	<- unique(df.s)
		set(df.s, NULL, 'LABEL', df.s[, paste(sapply(strsplit(LABEL, '_', fixed=1),'[[',1), '.',sapply(strsplit(LABEL, '_', fixed=1),'[[',2),sep='') ] )		
		#	from comments variable and from what T Leitner sent
		df.trm	<- list(data.table(FROM='P1', TO=c('P2','P5','P7','P8','P11'), TIME_TR_MIN=c('01-04-1983','01-08-1982','01-06-1982','01-10-1980','01-11-1981'), TIME_TR_MAX=c('01-05-1983','01-11-1982','01-9-1982','01-05-1981','01-01-1982'), MODE=c('HET','HET','HET','HET','HET')),
					data.table(FROM='P8', TO='P10', TIME_TR_MIN=c('01-05-1986'), TIME_TR_MAX=c('01-12-1986'), MODE= 'HET'),
					data.table(FROM='P8', TO='P9', TIME_TR_MIN=c('24-03-1984'), TIME_TR_MAX=c('24-03-1984'), MODE= 'MTC'),
					data.table(FROM='P5', TO='P6', TIME_TR_MIN=c('01-01-1983'), TIME_TR_MAX=c('01-05-1983'), MODE= 'HET'),
					data.table(FROM='P2', TO='P3', TIME_TR_MIN=c('04-09-1985'), TIME_TR_MAX=c('04-09-1985'), MODE= 'MTC')	)
		df.trm	<- do.call('rbind',df.trm)				
		set(df.trm, NULL, 'TIME_TR_MIN', df.trm[, hivc.db.Date2numeric(as.Date(TIME_TR_MIN, format='%d-%m-%Y'))])
		set(df.trm, NULL, 'TIME_TR_MAX', df.trm[, hivc.db.Date2numeric(as.Date(TIME_TR_MAX, format='%d-%m-%Y'))])
		df.trm[, TIME_TR:=df.trm[, (TIME_TR_MIN+TIME_TR_MAX)/2]]
		#
		#	load seq divergence
		#	
		infile				<- '141105_set2_seqs_INFO.R'
		file				<- paste(indir, '/',infile, sep='')
		load(file)
		#
		#	GAG	overall brl distribution in transmission pairs
		#		
		brl					<- brl.gag		
		brl[, PATIENT1:=brl[, sapply(strsplit(TIP1,'.',fixed=1),'[[',1)]]
		brl[, PATIENT2:=brl[, sapply(strsplit(TIP2,'.',fixed=1),'[[',1)]]		
		brl	<- subset(brl, PATIENT1!=PATIENT2)
		#	extract sequence sampling year
		setnames(df.s, c('LABEL','PosSeqT'),c('TIP1', 'PATIENT1_PosSeqT'))
		brl		<- merge(brl, subset(df.s, select=c(TIP1, PATIENT1_PosSeqT)), by='TIP1')
		setnames(df.s, c('TIP1', 'PATIENT1_PosSeqT'), c('TIP2', 'PATIENT2_PosSeqT'))
		brl		<- merge(brl, subset(df.s, select=c(TIP2, PATIENT2_PosSeqT)), by='TIP2')		
		#	get other way round in transmission pair
		tmp		<- copy(brl)
		setnames(tmp, c('PATIENT1','PATIENT2','PATIENT1_PosSeqT','PATIENT2_PosSeqT'),c('PATIENT2','PATIENT1','PATIENT2_PosSeqT','PATIENT1_PosSeqT'))
		brl		<- rbind(brl, tmp, use.names=TRUE)
		#	merge with transmission pairs
		setnames(brl, c('PATIENT1','PATIENT2','PATIENT1_PosSeqT','PATIENT2_PosSeqT'), c('FROM','TO','FROM_SeqT','TO_SeqT'))
		set(brl, NULL, 'FROM', brl[, gsub('p','P',FROM)])
		set(brl, NULL, 'TO', brl[, gsub('p','P',TO)])
		trm		<- merge(df.trm, brl, by=c('FROM','TO'))
		#
		#	set times between sampling times etc
		trm[, d_SeqT:=abs(TO_SeqT-FROM_SeqT)]
		trm[, d_TSeqT:= abs(TO_SeqT-TIME_TR) + abs(FROM_SeqT-TIME_TR)]
		trm.gag	<- copy(trm)
		brl.gag	<- brl		
		#
		#	ENV	overall brl distribution in transmission pairs
		#		
		brl					<- brl.env		
		brl[, PATIENT1:=brl[, sapply(strsplit(TIP1,'.',fixed=1),'[[',1)]]
		brl[, PATIENT2:=brl[, sapply(strsplit(TIP2,'.',fixed=1),'[[',1)]]		
		brl	<- subset(brl, PATIENT1!=PATIENT2)
		#	extract sequence sampling year
		brl		<- merge(brl, subset(df.s, select=c(TIP2, PATIENT2_PosSeqT)), by='TIP2')
		setnames(df.s, c('TIP2', 'PATIENT2_PosSeqT'), c('TIP1', 'PATIENT1_PosSeqT'))
		brl		<- merge(brl, subset(df.s, select=c(TIP1, PATIENT1_PosSeqT)), by='TIP1')
		#	get other way round in transmission pair
		tmp		<- copy(brl)
		setnames(tmp, c('PATIENT1','PATIENT2','PATIENT1_PosSeqT','PATIENT2_PosSeqT'),c('PATIENT2','PATIENT1','PATIENT2_PosSeqT','PATIENT1_PosSeqT'))
		brl		<- rbind(brl, tmp, use.names=TRUE)
		#	merge with transmission pairs
		setnames(brl, c('PATIENT1','PATIENT2','PATIENT1_PosSeqT','PATIENT2_PosSeqT'), c('FROM','TO','FROM_SeqT','TO_SeqT'))
		set(brl, NULL, 'FROM', brl[, gsub('p','P',FROM)])
		set(brl, NULL, 'TO', brl[, gsub('p','P',TO)])
		trm		<- merge(df.trm, brl, by=c('FROM','TO'))
		#
		#	set times between sampling times etc
		trm[, d_SeqT:=abs(TO_SeqT-FROM_SeqT)]
		trm[, d_TSeqT:= abs(TO_SeqT-TIME_TR) + abs(FROM_SeqT-TIME_TR)]
		trm.env	<- copy(trm)
		brl.env	<- brl
		#
		#	CONC	overall brl distribution in transmission pairs
		#		
		brl					<- brl.conc		
		brl[, PATIENT1:=brl[, sapply(strsplit(TIP1,'.',fixed=1),'[[',1)]]
		brl[, PATIENT2:=brl[, sapply(strsplit(TIP2,'.',fixed=1),'[[',1)]]		
		brl	<- subset(brl, PATIENT1!=PATIENT2)
		#	extract sequence sampling year
		brl		<- merge(brl, subset(df.s, select=c(TIP1, PATIENT1_PosSeqT)), by='TIP1')
		setnames(df.s, c('TIP1', 'PATIENT1_PosSeqT'), c('TIP2', 'PATIENT2_PosSeqT'))
		brl		<- merge(brl, subset(df.s, select=c(TIP2, PATIENT2_PosSeqT)), by='TIP2')
		#	get other way round in transmission pair
		tmp		<- copy(brl)
		setnames(tmp, c('PATIENT1','PATIENT2','PATIENT1_PosSeqT','PATIENT2_PosSeqT'),c('PATIENT2','PATIENT1','PATIENT2_PosSeqT','PATIENT1_PosSeqT'))
		brl		<- rbind(brl, tmp, use.names=TRUE)
		#	merge with transmission pairs
		setnames(brl, c('PATIENT1','PATIENT2','PATIENT1_PosSeqT','PATIENT2_PosSeqT'), c('FROM','TO','FROM_SeqT','TO_SeqT'))
		set(brl, NULL, 'FROM', brl[, gsub('p','P',FROM)])
		set(brl, NULL, 'TO', brl[, gsub('p','P',TO)])
		trm		<- merge(df.trm, brl, by=c('FROM','TO'))
		#
		#	set times between sampling times etc
		trm[, d_SeqT:=abs(TO_SeqT-FROM_SeqT)]
		trm[, d_TSeqT:= abs(TO_SeqT-TIME_TR) + abs(FROM_SeqT-TIME_TR)]
		trm.conc<- copy(trm)
		brl.conc<- brl
		
		#
		#	save
		#
		infile	<- '141105_set2_seqs_INFO_TRM.R'		
		file	<- paste(indir, '/',infile, sep='')
		save(file=file, seq.env, seq.gag, seq.conc, brl.env, brl.gag, brl.conc, ph.env, ph.gag, ph.conc, trm.env, trm.gag, trm.conc) 		
	}
	#	plot association 
	if(1)
	{
		indir	<- '/Users/Oliver/duke/2014_HIV_TChainSwedish'
		infile	<- '141105_set2_seqs_INFO_TRM.R'
		file	<- paste(indir, '/',infile, sep='')
		load(file)
		
		if(1)
		{
			#
			#	plot GAG	divergence
			#
			tmp				<- subset(trm.gag, select=c(d_SeqT, d_TSeqT, BRL))			
			trm.gag.swGA11	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(tmp), family=GA(mu.link='identity',sigma.link='identity'), n.cyc = 40)			
			
			trm.gag.p2		<- data.table(d_TSeqT=seq(0,15,0.1))	
			trm.gag.p2[, BRL_p:=predict(trm.gag.swGA11, data=tmp, newdata=as.data.frame(trm.gag.p2), type='response', se.fit=FALSE)]
			trm.gag.p3		<- gamlss.centiles.get(trm.gag.swGA11, tmp$d_TSeqT, cent = c(2.5, 97.5) )		
			setkey(trm.gag.p3, x)
			trm.gag.p3		<- unique(trm.gag.p3)			
						
			ggplot(trm.gag, aes(x=d_TSeqT)) + geom_point(aes(y=BRL, pch=MODE), size=2) + geom_line(data=trm.gag.p2, aes(y=BRL_p)) +
					geom_line(data=trm.gag.p3, aes(x=x, y=q2.5), lty='dashed') +
					geom_line(data=trm.gag.p3, aes(x=x, y=q97.5), lty='dashed') +					
					labs(x='time elapsed\n(years)', y='evolutionary divergence\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,13), ylim=c(0,0.115)) +
					scale_shape_manual(values=c(16,15), name='transmission mode') + 					
					scale_y_continuous(breaks=seq(0, 0.2, 0.02)) + scale_x_continuous(breaks=seq(0,30,2)) + 
					theme_bw() + theme(legend.key = element_blank(), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), legend.justification=c(0,1), legend.position=c(0,1), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.2)) +
					guides(pch=guide_legend(override.aes = list(size=4)))
			file	<- paste(indir, '/',"141105_set2_seqs_gag_patristic_dTS_data.pdf", sep='')
			ggsave(file=file, w=6, h=6)
			#
			#slope: 0.00396 for GAG
			#
			
			
			#
			#	plot ENV	divergence
			#
			tmp				<- subset(trm.env, select=c(d_SeqT, d_TSeqT, BRL))			
			trm.env.swGA11	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(tmp), family=GA(mu.link='identity',sigma.link='identity'), n.cyc = 40)			
			
			trm.env.p2		<- data.table(d_TSeqT=seq(0,15,0.1))	
			trm.env.p2[, BRL_p:=predict(trm.env.swGA11, data=tmp, newdata=as.data.frame(trm.env.p2), type='response', se.fit=FALSE)]
			trm.env.p3		<- gamlss.centiles.get(trm.env.swGA11, tmp$d_TSeqT, cent = c(2.5, 97.5) )		
			setkey(trm.env.p3, x)
			trm.env.p3		<- unique(trm.env.p3)			
			
			ggplot(trm.env, aes(x=d_TSeqT)) + geom_point(aes(y=BRL, pch=MODE), size=2) + geom_line(data=trm.env.p2, aes(y=BRL_p)) +
					geom_line(data=trm.env.p3, aes(x=x, y=q2.5), lty='dashed') +
					geom_line(data=trm.env.p3, aes(x=x, y=q97.5), lty='dashed') +					
					labs(x='time elapsed\n(years)', y='evolutionary divergence\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,13), ylim=c(0,0.2)) +
					scale_shape_manual(values=c(16,15), name='transmission mode') + 					
					scale_y_continuous(breaks=seq(0, 0.2, 0.02)) + scale_x_continuous(breaks=seq(0,30,2)) + 
					theme_bw() + theme(legend.key = element_blank(), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), legend.justification=c(0,1), legend.position=c(0,1), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.2)) +
					guides(pch=guide_legend(override.aes = list(size=4)))
			file	<- paste(indir, '/',"141105_set2_seqs_env_patristic_dTS_data.pdf", sep='')
			ggsave(file=file, w=6, h=6)
			#
			#	ENV slope	0.01137
			#			
		}
	}
}
######################################################################################
project.Eduan.cluster.results<- function()
{
	require(data.table)
	require(ggplot2)
	
	df	<- as.data.table(read.csv('~/Downloads/HighCov.csv'))
	set(df, 24L, 'FACTOR', 'Time.to.Diag>6M')
	setnames(df, colnames(df), gsub('COR','clustering odds ratio',gsub('TTR','true transmission rate ratio',gsub('TOD','true odds ratio being transmitter',colnames(df)))))
	df	<- suppressWarnings(melt(df, id.vars=c('Coverage','FACTOR')))
	df[, STAT:=gsub('.*_([A-Za-z]+)$','\\1', variable)]
	set(df, NULL, 'variable', df[, gsub('(.*)_[A-Za-z]+$','\\1', variable)])	
	df	<- dcast.data.table(df, FACTOR+variable~STAT, value.var='value')	
	df	<- subset(df, FACTOR%in%c(	'Age.at.Diag.25to29','Age.at.Diag.less25',
									'ART+','Male','HighRisk',
									'Recency<1Y','Recency<6M','Recency<3M',
									'Time.to.Diag<1Y','Time.to.Diag<6M',
									'Time.to.ART<6M'))
	set(df, NULL, 'FACTOR', df[, factor(FACTOR, levels=c(	'Male','HighRisk','ART+','Time.to.ART<6M','Age.at.Diag.less25','Age.at.Diag.25to29','Recency<3M','Recency<6M','Recency<1Y','Time.to.Diag<6M','Time.to.Diag<1Y'),
												labels=c(	'Male vs\nFemale',
															'high risk behaviour vs\nmedium risk behaviour',
															'ever ART started vs\nnever ART started',
															'ART started in first 6 monthsvs\nART started after first 6 months or never',
															'Age at diagnosis <25 yrs vs\n>30 yrs',
															'Age at diagnosis 25-29 yrs vs\n>30 yrs',
															'In the first 3 months of infection vs\nafter the first 3 months',
															'In the first 6 months of infection vs\nafter the first 6 months',
															'In the first 12 months of infection vs\nafter the first 12 months',															
															'Diagnosed in the first 6 months of infection vs\nafter the first 6 months',
															'Diagnosed in the first 12 months of infection vs\nafter the first 12 months'
															))])
	ggplot(df, aes(x=FACTOR, y=Mean, ymin=Lower, ymax=Upper, color=variable)) + 
			geom_point(position=position_dodge(width=0.5)) + 
			geom_errorbar(position=position_dodge(width=0.5)) + 
			coord_flip() + 
			scale_y_log10(breaks=c(0.1, 0.33, 0.5, 0.66, 1, 1.5, 2, 3, 10), limits=c(1/3,3)) +
			theme_bw() + theme(legend.position='bottom') +
			labs(x='', colour='', y='')
}
######################################################################################
project.Tchain.Belgium.sensecheck	<- function()
{
	require(adephylo)
	indir		<- '/Users/Oliver/duke/2014_HIV_TChainBelgium'
	#	get readable fasta file
	if(0)
	{
		#infile.online	<- 'set7_pol.fasta'
		infile.online	<- 'set7_env.fasta'
		#infile			<- '140919_set7_pol.fasta'
		infile			<- '140919_set7_env.fasta'
		
		file		<- paste(indir,'/',infile.online, sep='')
		txt			<- scan(file = file, what = "", sep = "\n", quiet = TRUE)
		#	search for '>' names 
		tmp			<- c( which( grepl('^>',txt) ), length(txt)+1 )
		tmp			<- lapply(seq_len(length(tmp)-1), function(i)
				{				
					data.table( LABEL= gsub('>','',txt[tmp[i]]), SEQ= paste(txt[(tmp[i]+1):(tmp[i+1]-1)], sep='', collapse='') )
				})
		seq				<- do.call('rbind', tmp)	
		tmp				<- tolower(do.call('rbind',strsplit(seq[, SEQ],'')))
		rownames(tmp)	<- seq[, LABEL]
		seq.pol			<- as.DNAbin(tmp)
		
		file			<- paste(indir,'/',infile, sep='')
		write.dna(file=file, seq.pol, format='fasta')
	}
	#
	#	check if there are zero distance seqs
	if(0)
	{
		indir	<- '/Users/Oliver/duke/2014_HIV_TChainBelgium'
		infile	<- '141105_Vrancken_pol_2.fasta'
		file	<- paste(indir,'/',infile, sep='')
		seq.pol	<- read.dna(file, format='fa')		
		seq.pol	<- as.character(seq.pol)
		seq.pol <- seq.pol[, 2097:3515]
		seq.ext	<- as.DNAbin(seq.pol)
		
		infile	<- '140921_set7_Drcheck_HIValign2.fasta'
		file	<- paste(indir,'/',infile, sep='')
		seq.pol	<- read.dna(file=file, format='fa')		
		seq.pol	<- as.character(seq.pol	)
		seq.pol	<- as.DNAbin(seq.pol[, 13:1431])
		rownames(seq.pol)	<- paste(rownames(seq.pol),'STUDY',sep='_')
		
		seq.ext	<- rbind(seq.pol, seq.ext)
		infile	<- '141105_Vrancken_pol_ext.fasta'
		file	<- paste(indir,'/',infile, sep='')
		write.dna(file=file, seq.ext, format='fasta')
		#	cleaning of alignment
		seq.ext			<- as.character(seq.ext)	
		query.yes		<- seq.find(seq.ext, 89, c("-","-","-",'g','a','g','g'))	
		seq.ext[query.yes,89:95]<- matrix( c('g','a','g','g',"-","-","-"), nrow=length(query.yes), ncol=7, byrow=1 )
		query.yes		<- seq.find(seq.ext, 89, c("-","-","-",'g','a','a','g'))	
		seq.ext[query.yes,89:95]<- matrix( c('g','a','a','g',"-","-","-"), nrow=length(query.yes), ncol=7, byrow=1 )
		query.yes		<- seq.find(seq.ext, 89, c("-","-","-",'g','g','g','g'))	
		seq.ext[query.yes,89:95]<- matrix( c('g','g','g','g',"-","-","-"), nrow=length(query.yes), ncol=7, byrow=1 )		
		query.yes		<- seq.find(seq.ext, 101, c("-","-","-"))
		seq.ext[query.yes,96:103]<- cbind( matrix( '-', nrow=length(query.yes), ncol=3 ), seq.ext[query.yes,96:100])		
		query.yes		<- seq.find(seq.ext, 143, c("-","-","-"))
		seq.ext[query.yes,140:145]<- cbind( matrix( '-', nrow=length(query.yes), ncol=3 ), seq.ext[query.yes,140:142])		
		infile	<- '141105_Vrancken_pol_ext2.fasta'
		file	<- paste(indir,'/',infile, sep='')
		write.dna(file=file, seq.ext, format='fasta')
		#	manual clean
		seq.ext	<- read.dna(file, format='fa')		
		#	distances		
		seq.d	<- dist.dna(seq.ext, model='N', as.matrix=1)		
		tmp		<- lapply(seq_len(nrow(seq.d)), function(i)
					{
						tmp	<- setdiff(which(seq.d[i,]==0),i)
						if(!length(tmp))	tmp<- NA_integer_
						data.table(FROM=i, EQU=tmp)
					})
		tmp		<- do.call('rbind',tmp)
		seq.eq	<- subset(tmp, !is.na(EQU))
		set(seq.eq, NULL, 'FROM', seq.eq[, rownames(seq.d)[FROM]])
		set(seq.eq, NULL, 'EQU', seq.eq[, colnames(seq.d)[EQU]])
		
		seq.eq	<- subset( seq.eq, grepl('STUDY',FROM) & !grepl('HXB2', FROM) )
		seq.eq[, EQU.ID:=NA_character_]
		tmp		<- seq.eq[, which(!grepl('STUDY',EQU))]
		set(seq.eq, tmp, 'EQU.ID', seq.eq[tmp, sapply(strsplit(EQU,'.',fixed=1),'[[',4) ])
		set(seq.eq, tmp, 'EQU.ID', seq.eq[tmp, paste(EQU.ID,'STUDY',sep='_') ])
		
		subset( seq.eq, FROM!=EQU.ID )		
		# 	none are equal between individuals ... 		
	}
	#
	#	check if drug resistance mutations are masked
	if(0)
	{
		infile	<- '140921_set7_INFO_withDR.R'
		file	<- paste(indir, '/',infile, sep='')
		load(file)	#expect seq.env, seq.pol, brl.env, brl.pol, ph.env, ph.pol
		#
		file	<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
		cat(paste('\nLoading outgroup seq from file', file))
		load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
		#	align against HXB2
		seq.pol <- c(as.list(outgroup.seq.pol), as.list(seq.pol))		
		infile	<- '140921_set7_Drcheck.fasta'
		file	<- paste(indir,'/',infile, sep='')
		write.dna(file=file, seq.pol, format='fasta')
		#	hand-curate: there are a few double runs!
		infile	<- '140921_set7_Drcheck_HIValign2.fasta'
		file	<- paste(indir,'/',infile, sep='')
		seq.pol	<- read.dna(file=file, format='fa')		
		seq.pol	<- as.character(seq.pol	)
		seq.pol	<- as.DNAbin(seq.pol[, 13:1431])
		# first site is 2097 in HXB2
	
		#load drug resistance mutations and select unique mutants by codon
		load( system.file(package="hivclust", "data",'IAS_primarydrugresistance_201303.rda') )		
		dr				<- as.data.table(IAS_primarydrugresistance_201303)		
		alignment.start	<- 2097	
		set(dr, NULL, "Alignment.nuc.pos", dr[,Alignment.nuc.pos]-alignment.start+1)
		#if alignment equals any of the drug resistance mutants, replace with NNN	
		seq.pol			<- seq.rm.drugresistance(as.character(seq.pol), dr, verbose=1, rtn.DNAbin=1 )	
		#rm HX1B and close gaps
		seq.pol			<- seq.pol[!grepl('HXB2',rownames(seq.pol)),]
		seq.pol			<- seq.rmallchar(seq.pol, rm.char='-', verbose=1)
		
		#	check visually
		infile			<- "140921_pol_norecombnodr.fasta"
		file			<- paste(indir,'/',infile, sep='')
		write.dna(file=file, seq.pol, format='fasta')
	}
	#	
	#	collect seq, ExaML phylos and patristic distances
	if(1)
	{
		indir			<- '/Users/Oliver/duke/2014_HIV_TChainBelgium'
		infile.pol		<- "140921_pol_norecombnodr.fasta"
		file			<- paste(indir,'/',infile.pol, sep='')
		seq.pol			<- read.dna(file, format='fa')
		infile.env		<- "140921_env_norecomb.fasta"
		file			<- paste(indir,'/',infile.env, sep='')
		seq.env			<- read.dna(file, format='fa')
		#
		#	run ExaML to get trees in same way
		#
		seq					<- seq.pol	
		infile				<- infile.pol
		infile.seq.sig		<- "Wed_Sep_24_12:59:06_2013"
		infile.seq			<- paste(substr(infile,1,nchar(infile)-6),'',sep='')
		file				<- paste( indir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
		save(seq, file=file)
		cmd					<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
		cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
		tmp					<- hivc.cmd.hpccaller(indir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		Sys.sleep(1)	
		seq					<- seq.env	
		infile				<- infile.env
		infile.seq.sig		<- "Wed_Sep_24_12:59:06_2013"
		infile.seq			<- paste(substr(infile,1,nchar(infile)-6),'',sep='')
		file				<- paste( indir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
		save(seq, file=file)
		#	run ExaML
		cmd					<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
		cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
		tmp					<- hivc.cmd.hpccaller(indir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		Sys.sleep(1)
		#
		#	get patristic distance
		#	
		file				<- list.files(indir, paste('^ExaML_result.140921_*','pol_norecombnodr','.*finaltree.000$',sep=''), full.names=TRUE)
		ph.pol				<- read.tree(file)
		brl.pol				<- as.matrix( distTips(ph.pol , method='patristic') )
		tmp					<- as.data.table(expand.grid(TIP1=rownames(brl.pol), TIP2=colnames(brl.pol), stringsAsFactors = FALSE))
		tmp					<- tmp[, list(BRL= brl.pol[TIP1, TIP2]), by=c('TIP1','TIP2')]
		tmp					<- subset(tmp, TIP1!=TIP2)
		tmp					<- tmp[, {
					if(TIP1<TIP2)	z<- list(TIP1n=TIP1, TIP2n=TIP2, BRL=BRL)
					if(TIP1>TIP2)	z<- list(TIP1n=TIP2, TIP2n=TIP1, BRL=BRL)
					z
				}, by=c('TIP1','TIP2')]
		set(tmp, NULL, c('TIP1','TIP2'), NULL)
		setnames(tmp, c('TIP1n','TIP2n'), c('TIP1','TIP2'))
		setkey(tmp, TIP1, TIP2)
		brl.pol				<- unique(tmp)
		#							
		file				<- list.files(indir, paste('^ExaML_result.140921_*','env','.*finaltree.000$',sep=''), full.names=TRUE)
		ph.env				<- read.tree(file)
		brl.env				<- as.matrix( distTips(ph.pol , method='patristic') )
		tmp					<- as.data.table(expand.grid(TIP1=rownames(brl.env), TIP2=colnames(brl.env), stringsAsFactors = FALSE))
		tmp					<- tmp[, list(BRL= brl.env[TIP1, TIP2]), by=c('TIP1','TIP2')]
		tmp					<- subset(tmp, TIP1!=TIP2)
		tmp					<- tmp[, {
					if(TIP1<TIP2)	z<- list(TIP1n=TIP1, TIP2n=TIP2, BRL=BRL)
					if(TIP1>TIP2)	z<- list(TIP1n=TIP2, TIP2n=TIP1, BRL=BRL)
					z
				}, by=c('TIP1','TIP2')]
		set(tmp, NULL, c('TIP1','TIP2'), NULL)
		setnames(tmp, c('TIP1n','TIP2n'), c('TIP1','TIP2'))
		setkey(tmp, TIP1, TIP2)
		brl.env				<- unique(tmp)
		#
		infile				<- '140921_set7_INFO.R'
		file				<- paste(indir, '/',infile, sep='')
		save(seq.env, seq.pol, brl.env, brl.pol, ph.env, ph.pol, file=file)
	}	
	#	collect seq, phylos and patristic distances
	if(0)
	{
		infile.pol			<- '140919_set7_pol.fasta'
		infile.env			<- '140919_set7_env.fasta'
		file				<- paste(indir,'/',infile.pol, sep='')
		seq.pol				<- read.dna(file, format='fasta')
		file				<- paste(indir,'/',infile.env, sep='')
		seq.env				<- read.dna(file, format='fasta')
		#
		#	run ExaML to get trees in same way
		#
		seq					<- seq.pol	
		infile				<- infile.pol
		infile.seq.sig		<- "Sun_Sep_14_12:59:06_2013"
		infile.seq			<- paste(substr(infile,1,nchar(infile)-6),'',sep='')
		file				<- paste( indir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
		save(seq, file=file)
		cmd					<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
		cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
		tmp					<- hivc.cmd.hpccaller(indir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		Sys.sleep(1)	
		seq					<- seq.env	
		infile				<- infile.env
		infile.seq.sig		<- "Sun_Sep_14_12:59:06_2013"
		infile.seq			<- paste(substr(infile,1,nchar(infile)-6),'',sep='')
		file				<- paste( indir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
		save(seq, file=file)
		#	run ExaML
		cmd					<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
		cmd					<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
		tmp					<- hivc.cmd.hpccaller(indir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		Sys.sleep(1)
		#
		#	get patristic distance
		#
		file				<- list.files(indir, paste('^ExaML_result.*','pol','.*finaltree.000$',sep=''), full.names=TRUE)
		ph.pol				<- read.tree(file)
		brl.pol				<- as.matrix( distTips(ph.pol , method='patristic') )
		tmp					<- as.data.table(expand.grid(TIP1=rownames(brl.pol), TIP2=colnames(brl.pol), stringsAsFactors = FALSE))
		tmp					<- tmp[, list(BRL= brl.pol[TIP1, TIP2]), by=c('TIP1','TIP2')]
		tmp					<- subset(tmp, TIP1!=TIP2)
		tmp					<- tmp[, {
					if(TIP1<TIP2)	z<- list(TIP1n=TIP1, TIP2n=TIP2, BRL=BRL)
					if(TIP1>TIP2)	z<- list(TIP1n=TIP2, TIP2n=TIP1, BRL=BRL)
					z
				}, by=c('TIP1','TIP2')]
		set(tmp, NULL, c('TIP1','TIP2'), NULL)
		setnames(tmp, c('TIP1n','TIP2n'), c('TIP1','TIP2'))
		setkey(tmp, TIP1, TIP2)
		brl.pol				<- unique(tmp)
		#							
		file				<- list.files(indir, paste('^ExaML_result.*','env','.*finaltree.000$',sep=''), full.names=TRUE)
		ph.env				<- read.tree(file)
		brl.env				<- as.matrix( distTips(ph.env , method='patristic') )
		tmp					<- as.data.table(expand.grid(TIP1=rownames(brl.env), TIP2=colnames(brl.env), stringsAsFactors = FALSE))
		tmp					<- tmp[, list(BRL= brl.env[TIP1, TIP2]), by=c('TIP1','TIP2')]
		tmp					<- subset(tmp, TIP1!=TIP2)
		tmp					<- tmp[, {
					if(TIP1<TIP2)	z<- list(TIP1n=TIP1, TIP2n=TIP2, BRL=BRL)
					if(TIP1>TIP2)	z<- list(TIP1n=TIP2, TIP2n=TIP1, BRL=BRL)
					z
				}, by=c('TIP1','TIP2')]
		set(tmp, NULL, c('TIP1','TIP2'), NULL)
		setnames(tmp, c('TIP1n','TIP2n'), c('TIP1','TIP2'))
		setkey(tmp, TIP1, TIP2)
		brl.env				<- unique(tmp)
		#
		infile				<- '140919_set7_INFO.R'
		file				<- paste(indir, '/',infile, sep='')
		save(seq.env, seq.pol, brl.env, brl.pol, ph.env, ph.pol, file=file)
	}
	# 	get divergence-pair matrix
	if(1)
	{
		require(splines)
		require(gamlss)
		#
		#	sampling dates
		#
		infile	<- '140921_set7_sample_dates.csv'
		df.s	<- read.csv( paste(indir, '/', infile, sep=''), stringsAsFactors=FALSE )
		df.s	<- as.data.table(df.s)[1:25,]
		set(df.s, NULL, 'sample', df.s[, gsub('\\s','',sample)])
		df.s[, PosSeqT:=df.s[, as.Date(gsub('\\s','',sample.date), format='%d-%b-%Y')]]
		tmp		<- df.s[, which(is.na(PosSeqT))]
		set(df.s, tmp, 'PosSeqT', df.s[tmp, as.Date(gsub('\\s','',sample.date), format='%d-%m-%Y')])
		set(df.s, NULL, 'PosSeqT', df.s[, hivc.db.Date2numeric(PosSeqT)])		
		df.s	<- subset(df.s, select=c(sample, PosSeqT))
		tmp		<- df.s[, which(grepl('_1',sample,fixed=1))]
		set(df.s,tmp,'sample',df.s[tmp, substr(sample, 1, 3)])
		#
		#	transmission dates
		#
		infile	<- "140921_set7_transmission_dates.csv"
		df.trm	<- read.csv( paste(indir, '/', infile, sep=''), stringsAsFactors=FALSE )
		df.trm	<- as.data.table(df.trm)[1:11,]
		df.trm[, FROM:=df.trm[, sapply(strsplit(patients, '->'),'[[',1)]]
		df.trm[, TO:=df.trm[, sapply(strsplit(patients, '->'),'[[',2)]]
		set(df.trm, NULL, 'TO', df.trm[, gsub('\\s|\\*','',TO)])
		set(df.trm, NULL, 'FROM', df.trm[, gsub('\\s|\\*','',FROM)])
		df.trm[, TIME_TR_MIN:=df.trm[, sapply(strsplit(window.of.transmission, '-'),'[[',1)]]
		df.trm[, TIME_TR_MAX:=df.trm[, sapply(strsplit(window.of.transmission, '-'),'[[',2)]]
		set(df.trm, NULL, 'TIME_TR_MIN', df.trm[, gsub('\\s|\\*|\\?','',TIME_TR_MIN)])
		set(df.trm, NULL, 'TIME_TR_MAX', df.trm[, gsub('\\s|\\*|\\?','',TIME_TR_MAX)])		
		set(df.trm, NULL, 'TIME_TR_MIN', df.trm[, hivc.db.Date2numeric(as.Date(TIME_TR_MIN, format='%d/%m/%y'))])
		set(df.trm, NULL, 'TIME_TR_MAX', df.trm[, hivc.db.Date2numeric(as.Date(TIME_TR_MAX, format='%d/%m/%y'))])
		tmp		<- df.trm[, which(is.na(TIME_TR_MIN))]
		set(df.trm, tmp, 'TIME_TR_MIN', df.trm[tmp,TIME_TR_MAX-1])
		df.trm[, TIME_TR:=df.trm[, (TIME_TR_MIN+TIME_TR_MAX)/2]]
		df.trm	<- subset( df.trm, select=c(FROM, TO, TIME_TR) )
		#	A	M	MB 1995-08								H 1996-01
		#	B	F	MB 1990-09 -- 1992-12, 1996-01->death
		#	C	M											H 2000-01
		#	D	F											H 1998-01
		#	E	F
		#	F	F											H 2002-05
		#	G	M											H 2002-05
		#	H	M	MB 1996-08 -- 1997-11					H 1997-12
		#	I	M	MB 1997-08 -- 1999-12					H 2000-01
		#	J
		#	K
		#	L
		df.ind	<- data.table( 	ID		= c('A','B','C','D','E','F','G','H','I','J','K','L'),
								AnyT_T1	= c(1995+7/12, 1990+8/12, 2000, 1998, Inf, 2002+4/12, 2002+4/12, 1996+7/12, 1997+7/12, NA, NA, NA))						
		indir	<- '/Users/Oliver/duke/2014_HIV_TChainBelgium'
		infile	<- '140919_set7_INFO.R'
		file	<- paste(indir, '/',infile, sep='')
		load(file)	#expect seq.env, seq.pol, brl.env, brl.pol, ph.env, ph.pol
		#
		#	POL	overall brl distribution in transmission pairs
		#
		brl.pol[, PATIENT1:=brl.pol[,substr(TIP1,1,1)]]
		brl.pol[, PATIENT2:=brl.pol[,substr(TIP2,1,1)]]
		brl.pol	<- subset(brl.pol, PATIENT1!=PATIENT2)
		#	extract sequence sampling year
		brl.pol[, sample:=brl.pol[,regmatches(TIP1, regexpr('.*c',TIP1))]]
		set(brl.pol, NULL, 'sample', brl.pol[, substr(sample,1,nchar(sample)-1)])
		stopifnot(!length(setdiff(brl.pol[, unique(sample)], df.s[, sample])))		
		brl.pol	<- merge(brl.pol, df.s, by='sample')
		setnames(brl.pol, 'PosSeqT', 'PATIENT1_SeqT')
		brl.pol[, sample:=brl.pol[,regmatches(TIP2, regexpr('.*c',TIP2))]]
		set(brl.pol, NULL, 'sample', brl.pol[, substr(sample,1,nchar(sample)-1)])
		stopifnot(!length(setdiff(brl.pol[, unique(sample)], df.s[, sample])))		
		brl.pol	<- merge(brl.pol, df.s, by='sample')
		setnames(brl.pol, 'PosSeqT', 'PATIENT2_SeqT')
		brl.pol[, sample:=NULL]
		#	add Tb4S manually		
		df.tb4s	<- data.table( TIP1=unique(c(subset( brl.pol, grepl('A', TIP1) )[, TIP1], subset( brl.pol, grepl('A', TIP2) )[, unique(TIP2)])), Tb4S='Yes' )
		tmp		<- unique(c(subset( brl.pol, grepl('B', TIP1) )[, unique(TIP1)], subset( brl.pol, grepl('B', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp[grepl('B90',tmp)], Tb4S='No' ), data.table( TIP1=tmp[!grepl('B90',tmp)], Tb4S='Yes' ))
		tmp		<- unique(c(subset( brl.pol, grepl('C', TIP1) )[, unique(TIP1)], subset( brl.pol, grepl('C', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp, Tb4S='Yes' ) )
		tmp		<- unique(c(subset( brl.pol, grepl('D', TIP1) )[, unique(TIP1)], subset( brl.pol, grepl('D', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp, Tb4S='Yes' ) )
		tmp		<- unique(c(subset( brl.pol, grepl('E', TIP1) )[, unique(TIP1)], subset( brl.pol, grepl('E', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp, Tb4S='Unknown' ) )
		tmp		<- unique(c(subset( brl.pol, grepl('F', TIP1) )[, unique(TIP1)], subset( brl.pol, grepl('F', TIP2) )[, unique(TIP2)]))
		#subset( brl.pol, grepl('G', TIP2) )[, unique(TIP2)]
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp[grepl('F02',tmp)], Tb4S='No' ), data.table( TIP1=tmp[!grepl('F02',tmp)], Tb4S='Yes' ))
		tmp		<- unique(c(subset( brl.pol, grepl('H', TIP1) )[, unique(TIP1)], subset( brl.pol, grepl('H', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp[grepl('H96',tmp)], Tb4S='No' ), data.table( TIP1=tmp[!grepl('H96',tmp)], Tb4S='Yes' ))
		tmp		<- unique(c(subset( brl.pol, grepl('I', TIP1) )[, unique(TIP1)], subset( brl.pol, grepl('I', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp, Tb4S='Yes' ) )
		#subset( brl.pol, grepl('J', TIP2) )[, unique(TIP2)]
		tmp		<- unique(c(subset( brl.pol, grepl('K', TIP1) )[, unique(TIP1)], subset( brl.pol, grepl('K', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp, Tb4S='Unknown' ) )
		tmp		<- unique(c(subset( brl.pol, grepl('L', TIP1) )[, unique(TIP1)], subset( brl.pol, grepl('L', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp, Tb4S='Unknown' ) )
		setnames(df.tb4s, 'Tb4S', 'PATIENT1_Tb4S')
		brl.pol	<- merge(brl.pol, df.tb4s, by='TIP1')
		setnames(df.tb4s, c('TIP1','PATIENT1_Tb4S'), c('TIP2','PATIENT2_Tb4S'))
		brl.pol	<- merge(brl.pol, df.tb4s, by='TIP2')		
		#	get other way round in transmission pair
		tmp		<- copy(brl.pol)
		setnames(tmp, c('PATIENT1','PATIENT2','PATIENT1_SeqT','PATIENT2_SeqT','PATIENT1_Tb4S','PATIENT2_Tb4S'),c('PATIENT2','PATIENT1','PATIENT2_SeqT','PATIENT1_SeqT','PATIENT2_Tb4S','PATIENT1_Tb4S'))
		brl.pol	<- rbind(brl.pol, tmp, use.names=TRUE)
		#	merge with transmission pairs
		setnames(brl.pol, c('PATIENT1','PATIENT2','PATIENT1_SeqT','PATIENT2_SeqT','PATIENT1_Tb4S','PATIENT2_Tb4S'), c('FROM','TO','FROM_SeqT','TO_SeqT','FROM_Tb4S','TO_Tb4S'))
		trm.pol	<- merge(df.trm, brl.pol, by=c('FROM','TO'))
		trm.pol[, table(FROM_Tb4S, TO_Tb4S)]
		#
		tmp		<- trm.pol[, which( substr(TIP1,1,1)!=FROM )]
		tmp2	<- trm.pol[tmp, TIP2]
		set(trm.pol, tmp, 'TIP2', trm.pol[tmp, TIP1])
		set(trm.pol, tmp, 'TIP1', tmp2)
		#	set times between sampling times etc
		trm.pol[, d_SeqT:=abs(TO_SeqT-FROM_SeqT)]
		trm.pol[, d_TSeqT:= abs(TO_SeqT-TIME_TR) + abs(FROM_SeqT-TIME_TR)]		
		#	set Tb4S	
		trm.pol[, Tb4S:=NA_character_]
		set(trm.pol, trm.pol[, which(FROM_Tb4S=='Yes' & TO_Tb4S=='Yes')], 'Tb4S', 'both')
		set(trm.pol, trm.pol[, which(FROM_Tb4S=='Yes' & TO_Tb4S%in%c('No','Unknown'))], 'Tb4S', 'one')
		set(trm.pol, trm.pol[, which(FROM_Tb4S%in%c('No','Unknown') & TO_Tb4S=='Yes')], 'Tb4S', 'one')
		set(trm.pol, trm.pol[, which(FROM_Tb4S=='No' & TO_Tb4S=='No')], 'Tb4S', 'none')
		set(trm.pol, trm.pol[, which(FROM_Tb4S=='Unknown' & TO_Tb4S=='Unknown')], 'Tb4S', 'Unknown')
		#
		trm.pol[, withA:= trm.pol[,grepl('A',TIP1) | grepl('A',TIP2)]]
		trm.pol[, withB:= trm.pol[,grepl('B',TIP1) | grepl('B',TIP2)]]
		trm.pol[, withC:= trm.pol[,grepl('C',TIP1) | grepl('C',TIP2)]]
		trm.pol[, withD:= trm.pol[,grepl('D',TIP1) | grepl('D',TIP2)]]
		trm.pol[, withE:= trm.pol[,grepl('E',TIP1) | grepl('E',TIP2)]]
		trm.pol[, withF:= trm.pol[,grepl('F',TIP1) | grepl('F',TIP2)]]
		#trm.pol[, withG:= trm.pol[,grepl('G',TIP1) | grepl('G',TIP2)]]
		trm.pol[, withH:= trm.pol[,grepl('H',TIP1) | grepl('H',TIP2)]]
		trm.pol[, withI:= trm.pol[,grepl('I',TIP1) | grepl('I',TIP2)]]
		#trm.pol[, withJ:= trm.pol[,grepl('J',TIP1) | grepl('J',TIP2)]]
		trm.pol[, withK:= trm.pol[,grepl('K',TIP1) | grepl('K',TIP2)]]
		trm.pol[, withL:= trm.pol[,grepl('L',TIP1) | grepl('L',TIP2)]]
		#
		#	ENV	overall brl distribution in transmission pairs
		#
		brl.env[, PATIENT1:=brl.env[,substr(TIP1,1,1)]]
		brl.env[, PATIENT2:=brl.env[,substr(TIP2,1,1)]]
		brl.env	<- subset(brl.env, PATIENT1!=PATIENT2)
		#	extract sequence sampling year
		brl.env[, sample:=brl.env[,regmatches(TIP1, regexpr('.*c',TIP1))]]
		set(brl.env, NULL, 'sample', brl.env[, substr(sample,1,nchar(sample)-1)])
		stopifnot(!length(setdiff(brl.env[, unique(sample)], df.s[, sample])))		
		brl.env	<- merge(brl.env, df.s, by='sample')
		setnames(brl.env, 'PosSeqT', 'PATIENT1_SeqT')
		brl.env[, sample:=brl.env[,regmatches(TIP2, regexpr('.*c',TIP2))]]
		set(brl.env, NULL, 'sample', brl.env[, substr(sample,1,nchar(sample)-1)])
		stopifnot(!length(setdiff(brl.env[, unique(sample)], df.s[, sample])))		
		brl.env	<- merge(brl.env, df.s, by='sample')
		setnames(brl.env, 'PosSeqT', 'PATIENT2_SeqT')
		brl.env[, sample:=NULL]
		#	add Tb4S manually		
		df.tb4s	<- data.table( TIP1=unique(c(subset( brl.env, grepl('A', TIP1) )[, TIP1], subset( brl.env, grepl('A', TIP2) )[, unique(TIP2)])), Tb4S='Yes' )
		tmp		<- unique(c(subset( brl.env, grepl('B', TIP1) )[, unique(TIP1)], subset( brl.env, grepl('B', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp[grepl('B90',tmp)], Tb4S='No' ), data.table( TIP1=tmp[!grepl('B90',tmp)], Tb4S='Yes' ))
		tmp		<- unique(c(subset( brl.env, grepl('C', TIP1) )[, unique(TIP1)], subset( brl.env, grepl('C', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp, Tb4S='Yes' ) )
		tmp		<- unique(c(subset( brl.env, grepl('D', TIP1) )[, unique(TIP1)], subset( brl.env, grepl('D', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp, Tb4S='Yes' ) )
		tmp		<- unique(c(subset( brl.env, grepl('E', TIP1) )[, unique(TIP1)], subset( brl.env, grepl('E', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp, Tb4S='Unknown' ) )
		tmp		<- unique(c(subset( brl.env, grepl('F', TIP1) )[, unique(TIP1)], subset( brl.env, grepl('F', TIP2) )[, unique(TIP2)]))
		#subset( brl.env, grepl('G', TIP2) )[, unique(TIP2)]
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp[grepl('F02',tmp)], Tb4S='No' ), data.table( TIP1=tmp[!grepl('F02',tmp)], Tb4S='Yes' ))
		tmp		<- unique(c(subset( brl.env, grepl('H', TIP1) )[, unique(TIP1)], subset( brl.env, grepl('H', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp[grepl('H96',tmp)], Tb4S='No' ), data.table( TIP1=tmp[!grepl('H96',tmp)], Tb4S='Yes' ))
		tmp		<- unique(c(subset( brl.env, grepl('I', TIP1) )[, unique(TIP1)], subset( brl.env, grepl('I', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp, Tb4S='Yes' ) )
		#subset( brl.env, grepl('J', TIP2) )[, unique(TIP2)]
		tmp		<- unique(c(subset( brl.env, grepl('K', TIP1) )[, unique(TIP1)], subset( brl.env, grepl('K', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp, Tb4S='Unknown' ) )
		tmp		<- unique(c(subset( brl.env, grepl('L', TIP1) )[, unique(TIP1)], subset( brl.env, grepl('L', TIP2) )[, unique(TIP2)]))
		df.tb4s	<- rbind( df.tb4s, data.table( TIP1=tmp, Tb4S='Unknown' ) )
		setnames(df.tb4s, 'Tb4S', 'PATIENT1_Tb4S')
		brl.env	<- merge(brl.env, df.tb4s, by='TIP1')
		setnames(df.tb4s, c('TIP1','PATIENT1_Tb4S'), c('TIP2','PATIENT2_Tb4S'))
		brl.env	<- merge(brl.env, df.tb4s, by='TIP2')		
		#	get other way round in transmission pair
		tmp		<- copy(brl.env)
		setnames(tmp, c('PATIENT1','PATIENT2','PATIENT1_SeqT','PATIENT2_SeqT','PATIENT1_Tb4S','PATIENT2_Tb4S'),c('PATIENT2','PATIENT1','PATIENT2_SeqT','PATIENT1_SeqT','PATIENT2_Tb4S','PATIENT1_Tb4S'))
		brl.env	<- rbind(brl.env, tmp, use.names=TRUE)
		#	merge with transmission pairs
		setnames(brl.env, c('PATIENT1','PATIENT2','PATIENT1_SeqT','PATIENT2_SeqT','PATIENT1_Tb4S','PATIENT2_Tb4S'), c('FROM','TO','FROM_SeqT','TO_SeqT','FROM_Tb4S','TO_Tb4S'))
		trm.env	<- merge(df.trm, brl.env, by=c('FROM','TO'))
		trm.env[, table(FROM_Tb4S, TO_Tb4S)]
		#
		tmp		<- trm.env[, which( substr(TIP1,1,1)!=FROM )]
		tmp2	<- trm.env[tmp, TIP2]
		set(trm.env, tmp, 'TIP2', trm.env[tmp, TIP1])
		set(trm.env, tmp, 'TIP1', tmp2)
		#	set times between sampling times etc
		trm.env[, d_SeqT:=abs(TO_SeqT-FROM_SeqT)]
		trm.env[, d_TSeqT:= abs(TO_SeqT-TIME_TR) + abs(FROM_SeqT-TIME_TR)]		
		#	set Tb4S	
		trm.env[, Tb4S:=NA_character_]
		set(trm.env, trm.env[, which(FROM_Tb4S=='Yes' & TO_Tb4S=='Yes')], 'Tb4S', 'both')
		set(trm.env, trm.env[, which(FROM_Tb4S=='Yes' & TO_Tb4S%in%c('No','Unknown'))], 'Tb4S', 'one')
		set(trm.env, trm.env[, which(FROM_Tb4S%in%c('No','Unknown') & TO_Tb4S=='Yes')], 'Tb4S', 'one')
		set(trm.env, trm.env[, which(FROM_Tb4S=='No' & TO_Tb4S=='No')], 'Tb4S', 'none')
		set(trm.env, trm.env[, which(FROM_Tb4S=='Unknown' & TO_Tb4S=='Unknown')], 'Tb4S', 'Unknown')
		#
		#
		trm.env[, withA:= trm.env[,grepl('A',TIP1) | grepl('A',TIP2)]]
		trm.env[, withB:= trm.env[,grepl('B',TIP1) | grepl('B',TIP2)]]
		trm.env[, withC:= trm.env[,grepl('C',TIP1) | grepl('C',TIP2)]]
		trm.env[, withD:= trm.env[,grepl('D',TIP1) | grepl('D',TIP2)]]
		trm.env[, withE:= trm.env[,grepl('E',TIP1) | grepl('E',TIP2)]]
		trm.env[, withF:= trm.env[,grepl('F',TIP1) | grepl('F',TIP2)]]
		#trm.env[, withG:= trm.env[,grepl('G',TIP1) | grepl('G',TIP2)]]
		trm.env[, withH:= trm.env[,grepl('H',TIP1) | grepl('H',TIP2)]]
		trm.env[, withI:= trm.env[,grepl('I',TIP1) | grepl('I',TIP2)]]
		#trm.env[, withJ:= trm.env[,grepl('J',TIP1) | grepl('J',TIP2)]]
		trm.env[, withK:= trm.env[,grepl('K',TIP1) | grepl('K',TIP2)]]
		trm.env[, withL:= trm.env[,grepl('L',TIP1) | grepl('L',TIP2)]]
		#
		#	save
		#
		infile	<- '140921_set7_INFO_TRM.R'
		file	<- paste(indir, '/',infile, sep='')
		save(file=file, seq.env, seq.pol, brl.env, brl.pol, ph.env, ph.pol, trm.env, trm.pol) 
	}	
}

project.Tchain.Belgium.sensecheck.explore<- function()
{
	#
	#	explore
	#
	if(0)
	{
		indir	<- '/Users/Oliver/duke/2014_HIV_TChainBelgium'
		infile	<- '140921_set7_INFO_TRM.R'
		file	<- paste(indir, '/',infile, sep='')
		load(file)
		#	exclude B->A, same as A->B
		trm.pol	<- subset( trm.pol, !(FROM=='B' & TO=='A') )		
		
		if(0)
		{
			length(union(brl.pol[, FROM], brl.pol[, TO]))
			unique(subset(trm.pol, select=c(FROM,TO)))
			#	numbers before A removed
			trm.pol.nA		<- subset(trm.pol, withA==FALSE)
			unique(subset(trm.pol.nA, select=c(FROM,TO)))
			length(union(trm.pol.nA[, TIP1], trm.pol.nA[, TIP2]))
		}
		if(0)
		{
			#
			#	plot all
			#	divergence
			ggplot(trm.pol, aes(x=BRL, fill=Tb4S)) + geom_histogram(binwidth=0.0005) +
					theme(legend.position='bottom') +
					labs(fill='Treatment start before\nsequence sampling date', x='patristic distance\n(estimated subst/site)', y='sequence pairs of transmission pairs\n(number)')
			file	<- paste(indir, '/',"140921_set7_pol_patristic_all.pdf", sep='')
			ggsave(file=file, w=12, h=6)
			#
			ggplot(trm.pol, aes(x=d_TSeqT, y=BRL, colour=Tb4S)) + geom_point() + stat_smooth(method = "lm", formula= y ~ ns(x,3)) + facet_grid(.~Tb4S)
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_all.pdf", sep='')
			ggsave(file=file, w=12, h=6)		
			ggplot(trm.pol, aes(x=d_SeqT, y=BRL, colour=Tb4S)) + geom_point() + stat_smooth(method = "lm", formula= y ~ ns(x,3)) + facet_grid(.~Tb4S)
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dSS_all.pdf", sep='')
			ggsave(file=file, w=12, h=6)				
		}
		#
		#	fit linear model through origin
		#
		if(0)
		{
			trm.pol.data.b	<- subset(trm.pol, Tb4S=='both')
			trm.pol.ZAGA.b	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), nu.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.b), family=ZAGA(mu.link='identity'), n.cyc = 40)
			trm.pol.GA.b	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.b), family=GA(mu.link='identity'), n.cyc = 40)
			trm.pol.EXP.b	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), data=as.data.frame(trm.pol.data.b), family=EXP, n.cyc = 40)
			trm.pol.NO.b	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.b), family=NO, n.cyc = 40)
			AIC(trm.pol.ZAGA.b, trm.pol.GA.b, trm.pol.EXP.b, trm.pol.NO.b)
			trm.pol.data.o	<- subset(trm.pol, Tb4S=='one')
			trm.pol.ZAGA.o	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), nu.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.o), family=ZAGA(mu.link='identity'), n.cyc = 40)
			trm.pol.GA.o	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.o), family=GA(mu.link='identity'), n.cyc = 40)
			trm.pol.EXP.o	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), data=as.data.frame(trm.pol.data.o), family=EXP, n.cyc = 40)
			trm.pol.NO.o	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.o), family=NO, n.cyc = 40)
			AIC(trm.pol.ZAGA.o, trm.pol.GA.o, trm.pol.EXP.o, trm.pol.NO.o)
			sapply( list(trm.pol.GA.b, trm.pol.GA.o, trm.pol.GA.n, trm.pol.GA.u), Rsq )
			#	normal model fits best, but we need a positive model -- just Gamma?
			trm.pol.data.n	<- subset(trm.pol, Tb4S=='none')
			trm.pol.GA.n	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.n), family=GA(mu.link='identity'), n.cyc = 40)
			trm.pol.data.u	<- subset(trm.pol, Tb4S=='Unknown')
			trm.pol.GA.u	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.u), family=GA(mu.link='identity'), n.cyc = 40)
			#	predict
			trm.pol.p		<- data.table(d_TSeqT=seq(0,25,0.1))
			trm.pol.p[, both:=predict(trm.pol.GA.b, data=trm.pol.data.b, newdata=as.data.frame(trm.pol.p), type='response', se.fit=FALSE)]
			trm.pol.p[, one:=predict(trm.pol.GA.o, data=trm.pol.data.o, newdata=as.data.frame(subset(trm.pol.p, select=d_TSeqT)), type='response', se.fit=FALSE)]
			trm.pol.p[, none:=predict(trm.pol.GA.n, data=trm.pol.data.n, newdata=as.data.frame(subset(trm.pol.p, select=d_TSeqT)), type='response', se.fit=FALSE)]
			trm.pol.p[, Unknown:=predict(trm.pol.GA.u, data=trm.pol.data.u, newdata=as.data.frame(subset(trm.pol.p, select=d_TSeqT)), type='response', se.fit=FALSE)]
			trm.pol.p		<- melt(trm.pol.p, id.vars=c('d_TSeqT'), variable.name='Tb4S', value='BRL_p')
			#	plot
			ggplot(trm.pol, aes(x=d_TSeqT, y=BRL, colour=Tb4S)) + geom_point() + geom_line(data=trm.pol.p, aes(y=BRL_p)) + facet_grid(.~Tb4S) + scale_y_continuous(limits=c(0,0.15))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_all_linear.pdf", sep='')
			ggsave(file=file, w=12, h=6)
		}
		if(1)
		{
			trm.pol.data.b	<- subset(trm.pol, !grepl('A',TIP1) & !grepl('A',TIP2) & Tb4S=='both')
			trm.pol.ZAGA.b	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), nu.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.b), family=ZAGA(mu.link='identity'), n.cyc = 40)
			trm.pol.GA.b	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.b), family=GA(mu.link='identity'), n.cyc = 40)
			trm.pol.EXP.b	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), data=as.data.frame(trm.pol.data.b), family=EXP, n.cyc = 40)
			trm.pol.NO.b	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.b), family=NO, n.cyc = 40)
			AIC(trm.pol.ZAGA.b, trm.pol.GA.b, trm.pol.EXP.b, trm.pol.NO.b)
			trm.pol.data.o	<- subset(trm.pol, !grepl('A',TIP1) & !grepl('A',TIP2) & Tb4S=='one')
			trm.pol.ZAGA.o	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), nu.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.o), family=ZAGA(mu.link='identity'), n.cyc = 40)
			trm.pol.GA.o	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.o), family=GA(mu.link='identity'), n.cyc = 40)
			trm.pol.EXP.o	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), data=as.data.frame(trm.pol.data.o), family=EXP, n.cyc = 40)
			trm.pol.NO.o	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.o), family=NO, n.cyc = 40)
			AIC(trm.pol.ZAGA.o, trm.pol.GA.o, trm.pol.EXP.o, trm.pol.NO.o)
			sapply( list(trm.pol.GA.b, trm.pol.GA.o, trm.pol.GA.n, trm.pol.GA.u), Rsq )
			#	normal model fits best, but we need a positive model -- just Gamma?
			trm.pol.data.n	<- subset(trm.pol, !grepl('A',TIP1) & !grepl('A',TIP2) & Tb4S=='none')
			trm.pol.GA.n	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.n), family=GA(mu.link='identity'), n.cyc = 40)
			trm.pol.data.u	<- subset(trm.pol, !grepl('A',TIP1) & !grepl('A',TIP2) & Tb4S=='Unknown')
			trm.pol.GA.u	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.data.u), family=GA(mu.link='identity'), n.cyc = 40)
			#	predict
			trm.pol.p		<- data.table(d_TSeqT=seq(0,25,0.1))
			trm.pol.p[, both:=predict(trm.pol.GA.b, data=trm.pol.data.b, newdata=as.data.frame(trm.pol.p), type='response', se.fit=FALSE)]
			trm.pol.p[, one:=predict(trm.pol.GA.o, data=trm.pol.data.o, newdata=as.data.frame(subset(trm.pol.p, select=d_TSeqT)), type='response', se.fit=FALSE)]
			trm.pol.p[, none:=predict(trm.pol.GA.n, data=trm.pol.data.n, newdata=as.data.frame(subset(trm.pol.p, select=d_TSeqT)), type='response', se.fit=FALSE)]
			trm.pol.p[, Unknown:=predict(trm.pol.GA.u, data=trm.pol.data.u, newdata=as.data.frame(subset(trm.pol.p, select=d_TSeqT)), type='response', se.fit=FALSE)]
			trm.pol.p		<- melt(trm.pol.p, id.vars=c('d_TSeqT'), variable.name='Tb4S', value='BRL_p')
			#	plot
			ggplot(subset(trm.pol, !grepl('A',TIP1) & !grepl('A',TIP2)), aes(x=d_TSeqT, y=BRL, colour=Tb4S)) + geom_point() + geom_line(data=trm.pol.p, aes(y=BRL_p)) + facet_grid(.~Tb4S) + scale_y_continuous(limits=c(0,0.15))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_all_linear_noA.pdf", sep='')
			ggsave(file=file, w=12, h=6)
		}
		#
		#	since fits not too dissimilar, lump all together
		#
		if(0)	
		{
			#	all patients
			trm.pol.GA		<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol), family=GA(mu.link='identity'), n.cyc = 40)
			trm.pol.p2		<- data.table(d_TSeqT=seq(0,25,0.1))
			trm.pol.p2[, BRL_p:=predict(trm.pol.GA, data=trm.pol, newdata=as.data.frame(trm.pol.p2), type='response', se.fit=FALSE)]
			trm.pol.p3		<- gamlss.centiles.get(trm.pol.GA, trm.pol$d_TSeqT, cent = c(2.5, 97.5) )		
			setkey(trm.pol.p3, x)
			trm.pol.p3		<- unique(trm.pol.p3)
			
			#tmp				<- gamlss.centiles.get(trm.pol.GA, trm.pol.p2$d_TSeqT, cent = c(2.5, 97.5), with.ordering=TRUE )
			#setnames(tmp, c('x','q2.5','q97.5'), c('d_TSeqT', 'BRLql', 'BRLqu'))
			#trm.pol.p2		<- merge(trm.pol.p2, tmp, by='d_TSeqT')
			trm.pol[, summary(BRL)]
			#0.006431 0.041570 0.055770 0.062380 0.086620 0.121200 
			ggplot(trm.pol, aes(x=d_TSeqT)) + geom_jitter(aes(y=BRL), size=1.2, alpha=0.5, position = position_jitter(width = .1)) + geom_line(data=trm.pol.p2, aes(y=BRL_p)) + 
					geom_ribbon(data=trm.pol.p3, aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) +
					#geom_vline(xintercept=0.5+3+2*4.1, color="#E41A1C", linetype=2) +
					labs(x='time elapsed\n(years)', y='evolutionary divergence\nbetween confirmed transmission pairs\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,22), ylim=c(0,0.125)) +
					scale_y_continuous(breaks=seq(0, 0.2, 0.01)) +
					scale_x_continuous(breaks=seq(0,30,2)) + theme_bw() + theme(panel.grid.minor=element_line(colour="grey75", size=0.2), panel.grid.major=element_line(colour="grey75"))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_all_GAfit.pdf", sep='')
			ggsave(file=file, w=5, h=5)		
			Rsq(trm.pol.GA)	#0.5884969
		}
		
		trm.pol[, summary(BRL)]
		#0.006431 0.041570 0.055770 0.062380 0.086620 0.121200
		
		if(0)
		{
			#	stratify by patients
			tmp		<- melt(trm.pol, measure.vars=c('withA','withB','withC','withD','withE','withF','withH','withI','withK','withL'))
			set(tmp, NULL, 'variable', tmp[, gsub('with','pairs including\npatient ',variable)])
			ggplot(tmp, aes(x=d_TSeqT)) +
					geom_ribbon(aes(x=c(10.17,30), ymin=0, ymax=1), fill='grey30', alpha=0.2) +
					geom_ribbon(aes(x=c(0,10.17), ymin=0, ymax=1), fill='red', alpha=0.2) +
					geom_jitter(aes(y=BRL, colour=value), size=0.7, alpha=0.5, position = position_jitter(width = .1)) + 
					labs(x='time elapsed\n(years)', y='evolutionary divergence\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,22), ylim=c(0,0.125)) +
					scale_colour_manual(values=c('black',"#F0027F"), name='pairs with sequences\nfrom selected patient', guide = FALSE) +
					scale_y_continuous(breaks=seq(0, 0.2, 0.02)) +
					scale_x_continuous(breaks=seq(0,30,4), minor_breaks=seq(0,30,1)) + theme_bw() + 
					theme(panel.grid.minor=element_line(colour="grey75", size=0.2), panel.grid.major=element_line(colour="grey75")) +
					facet_wrap( ~ variable, ncol=3, scales = "free_x")
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_all_GAfit_byPatient.pdf", sep='')
			ggsave(file=file, w=10, h=10)	
		}		
		if(1)
		{
			#	
			#	all together, not patient A, POL
			#
			trm.pol.nA		<- subset(trm.pol, withA==FALSE, select=c(d_SeqT, d_TSeqT, BRL)) 
			#trm.pol.nA		<- rbind(trm.pol.nA, data.table(BRL=rep(1e-3,100), d_SeqT=5e-4, d_TSeqT=5e-4), fill=TRUE)			
			#trm.pol.GA		<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity', sigma.link='identity'), n.cyc = 40)
			#trm.pol.GA		<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=2)'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity'), n.cyc = 40)			
			#trm.pol.GA		<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=3)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=2)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity'), n.cyc = 40)
			if(0)
			{
				trm.pol.GA1		<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=3)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity', sigma.link='identity'), n.cyc = 40)
				trm.pol.GA2		<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=2)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=3)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity', sigma.link='identity'), n.cyc = 40)
				trm.pol.GA3		<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=3)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=3)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity', sigma.link='identity'), n.cyc = 40)
				trm.pol.GA4		<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=4)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=3)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity', sigma.link='identity'), n.cyc = 40)
				trm.pol.GA5		<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=5)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=3)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity', sigma.link='identity'), n.cyc = 40)
				trm.pol.GA6		<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=6)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=3)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity', sigma.link='identity'), n.cyc = 40)			
				#	AIC is a joke ..
				AIC(trm.pol.GA1,trm.pol.GA2,trm.pol.GA3,trm.pol.GA4,trm.pol.GA5,trm.pol.GA6)				
			}
			trm.pol.GA42	<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=4)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=2)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity', sigma.link='identity'), n.cyc = 40)
			trm.pol.GA		<- trm.pol.GA42
			#trm.pol.GA		<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=4, Boundary.knots=c(0,15))'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=1, knots=c(5), Boundary.knots=c(0,15))'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity', sigma.link='identity'), n.cyc = 40)
			#trm.pol.GA		<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=4, Boundary.knots=c(0,15))'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=2)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity'), n.cyc = 40)
			#trm.pol.GA		<- gamlss(as.formula('BRL ~ ns(d_TSeqT, knots=c(2,4), Boundary.knots=c(0,15))'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=2)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity'), n.cyc = 40)			
			#trm.pol.GA		<- gamlss(as.formula('BRL ~ ns(d_TSeqT, knots=c(2,6,9))'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=2)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity'), n.cyc = 40)
			trm.pol.p2		<- data.table(d_TSeqT=seq(0,15,0.1))	
			trm.pol.p2[, BRL_p:=predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(trm.pol.p2), type='response', se.fit=FALSE)]
			trm.pol.p3		<- gamlss.centiles.get(trm.pol.GA, trm.pol.nA$d_TSeqT, cent = c(2.5, 97.5) )		
			setkey(trm.pol.p3, x)
			trm.pol.p3		<- unique(trm.pol.p3)			
			ggplot(trm.pol.nA, aes(x=d_TSeqT)) + geom_jitter(aes(y=BRL), size=1.2, alpha=0.5, position = position_jitter(width = .1)) + geom_line(data=trm.pol.p2, aes(y=BRL_p)) + 
					geom_ribbon(data=trm.pol.p3, aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) +
					#geom_vline(xintercept=0.5+3+2*4.1, color="#E41A1C", linetype=2) +
					labs(x='time elapsed\n(years)', y='evolutionary divergence\nbetween confirmed transmission pairs\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,15), ylim=c(0,0.08)) +
					scale_colour_brewer(name='pairs with sequences\nfrom selected patient', palette='Set1') +
					scale_y_continuous(breaks=seq(0, 0.2, 0.01)) +
					scale_x_continuous(breaks=seq(0,30,2)) + theme_bw() + theme(panel.grid.minor=element_line(colour="grey75", size=0.2), panel.grid.major=element_line(colour="grey75"))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_notA_GAfit.pdf", sep='')
			ggsave(file=file, w=5, h=5)	
			Rsq(trm.pol.GA42)	#0.895
		}
		if(0)
		{
			#	
			#	all together, not patient A, POL
			#
			#	NORMAL MODEL
			trm.pol.nA		<- subset(trm.pol, withA==FALSE, select=c(d_SeqT, d_TSeqT, BRL))
			trm.pol.nA		<- rbind(trm.pol.nA, data.table(BRL=seq(-0.005-0.005,0.001-0.005,len=20), d_SeqT=-0.025, d_TSeqT=-0.025), fill=TRUE)
			trm.pol.NO42	<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=4)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=2)'), data=as.data.frame(trm.pol.nA), family=NO(mu.link='identity', sigma.link='identity'), n.cyc = 40)
			trm.pol.NO		<- trm.pol.NO42
			trm.pol.p2		<- data.table(d_TSeqT=seq(0,15,0.1))	
			trm.pol.p2[, BRL_p:=predict(trm.pol.NO, data=trm.pol.nA, newdata=as.data.frame(trm.pol.p2), type='response', se.fit=FALSE)]
			trm.pol.p3		<- gamlss.centiles.get(trm.pol.NO, trm.pol.nA$d_TSeqT, cent = c(2.5, 97.5) )		
			setkey(trm.pol.p3, x)
			trm.pol.p3		<- unique(trm.pol.p3)			
			ggplot(trm.pol.nA, aes(x=d_TSeqT)) + geom_jitter(aes(y=BRL), size=1.2, alpha=0.5, position = position_jitter(width = .1)) + geom_line(data=trm.pol.p2, aes(y=BRL_p)) + 
					geom_ribbon(data=trm.pol.p3, aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) +
					#geom_vline(xintercept=0.5+3+2*4.1, color="#E41A1C", linetype=2) +
					labs(x='time elapsed\n(years)', y='evolutionary divergence\nbetween confirmed transmission pairs\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,15), ylim=c(0,0.08)) +
					scale_colour_brewer(name='pairs with sequences\nfrom selected patient', palette='Set1') +
					scale_y_continuous(breaks=seq(0, 0.2, 0.01)) +
					scale_x_continuous(breaks=seq(0,30,2)) + theme_bw() + theme(panel.grid.minor=element_line(colour="grey75", size=0.2), panel.grid.major=element_line(colour="grey75"))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_notA_NOfit.pdf", sep='')
			ggsave(file=file, w=5, h=5)	
			Rsq(trm.pol.NO42)	#0.783
		}
		if(1)
		{
			#	
			#	all together, not patient A, POL
			#
			#	GAMMA MODEL
			trm.pol.nA		<- subset(trm.pol, withA==FALSE, select=c(d_SeqT, d_TSeqT, BRL))			
			#trm.pol.nA		<- rbind(trm.pol.nA, data.table(BRL=0.005, d_SeqT=seq(0.5,2,len=40), d_TSeqT=seq(0.05,2,len=40)), fill=TRUE)
			trm.pol.nA		<- rbind(trm.pol.nA, data.table(BRL=0.001, d_SeqT=rep(0.1,40), d_TSeqT=rep(0.1,40)), fill=TRUE)			
			trm.pol.GA32	<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=4)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=2)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity',sigma.link='identity'), n.cyc = 40)			
			#trm.pol.NO42	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), nu.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.nA), family=ZAGA(mu.link='identity',sigma.link='identity'), n.cyc = 40)			
			trm.pol.GA		<- trm.pol.GA32
			trm.pol.p2		<- data.table(d_TSeqT=seq(0,15,0.1))	
			trm.pol.p2[, BRL_p:=predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(trm.pol.p2), type='response', se.fit=FALSE)]
			trm.pol.p3		<- gamlss.centiles.get(trm.pol.GA, trm.pol.nA$d_TSeqT, cent = c(2.5, 97.5) )		
			setkey(trm.pol.p3, x)
			trm.pol.p3		<- unique(trm.pol.p3)			
			ggplot(trm.pol.nA, aes(x=d_TSeqT)) + geom_jitter(aes(y=BRL), size=1.2, alpha=0.5, position = position_jitter(width = .1), data=subset(trm.pol.nA, BRL>0.002)) + geom_line(data=trm.pol.p2, aes(y=BRL_p)) + 
					geom_ribbon(data=trm.pol.p3, aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) +
					#geom_vline(xintercept=0.5+3+2*4.1, color="#E41A1C", linetype=2) +
					labs(x='time elapsed\n(years)', y='evolutionary divergence\nbetween confirmed transmission pairs\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,15), ylim=c(0,0.08)) +
					scale_colour_brewer(name='pairs with sequences\nfrom selected patient', palette='Set1') +
					scale_y_continuous(breaks=seq(0, 0.2, 0.01)) +
					scale_x_continuous(breaks=seq(0,30,2)) + theme_bw() + theme(panel.grid.minor=element_line(colour="grey75", size=0.2), panel.grid.major=element_line(colour="grey75"))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_notA_GAfit.pdf", sep='')
			ggsave(file=file, w=5, h=5)
			
			ggplot(trm.pol.nA, aes(x=d_TSeqT)) + geom_jitter(aes(y=BRL), size=1.2, alpha=0.5, position = position_jitter(width = .1), data=subset(trm.pol.nA, BRL>0.002)) +  
					labs(x='time elapsed\n(years)', y='evolutionary divergence\nbetween confirmed transmission pairs\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,15), ylim=c(0,0.08)) +
					scale_colour_brewer(name='pairs with sequences\nfrom selected patient', palette='Set1') +
					scale_y_continuous(breaks=seq(0, 0.2, 0.01)) +
					scale_x_continuous(breaks=seq(0,30,2)) + theme_bw() + theme(panel.grid.minor=element_line(colour="grey75", size=0.2), panel.grid.major=element_line(colour="grey75"))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_notA_data.pdf", sep='')
			ggsave(file=file, w=5, h=5)	
			
			Rsq(trm.pol.GA32)	#0.90
		}
		if(1)
		{
			#	
			#	EXPLORE all except patient A together, POL, prior for low brl
			#
			#	GAMMA MODEL
			trm.pol.nA		<- subset(trm.pol, withA==FALSE, select=c(d_SeqT, d_TSeqT, BRL))			
			trm.pol.nA		<- rbind(trm.pol.nA, data.table(BRL=0.0025, d_SeqT=seq(2,4,len=200), d_TSeqT=seq(2,8,len=200)), fill=TRUE)
			#trm.pol.nA		<- rbind(trm.pol.nA, data.table(BRL=0.001, d_SeqT=rep(0.1,40), d_TSeqT=rep(0.1,40)), fill=TRUE)			
			#trm.pol.GA32	<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=4)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=2)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity',sigma.link='identity'), n.cyc = 40)			
			trm.pol.GA11e	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity',sigma.link='identity'), n.cyc = 40)			
			trm.pol.GA		<- trm.pol.GA11e
			trm.pol.p2		<- data.table(d_TSeqT=seq(0,18,0.1))	
			trm.pol.p2[, BRL_p:=predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(trm.pol.p2), type='response', se.fit=FALSE)]
			trm.pol.p3		<- gamlss.centiles.get(trm.pol.GA, trm.pol.nA$d_TSeqT, cent = c(2.5, 10, 25, 75, 90, 97.5) )		
			setkey(trm.pol.p3, x)
			trm.pol.p3		<- unique(trm.pol.p3)	
			colors			<- c("#377EB8","#E41A1C")
			line.size		<- 0.8
			ggplot(trm.pol, aes(x=d_TSeqT)) +
					geom_ribbon(aes(x=c(10.17,20), ymin=0, ymax=1), fill='grey30', alpha=0.2) +
					geom_ribbon(aes(x=c(0,10.17), ymin=0, ymax=1), fill=colors[2], alpha=0.2) +
					geom_jitter(aes(y=BRL), size=1.2, alpha=0.5, position = position_jitter(width = .1), data=subset(trm.pol, withA==FALSE & BRL>0.003)) +									
					geom_line(data=trm.pol.p2, aes(y=BRL_p)) +
					geom_line(data=trm.pol.p3, aes(x=x, y=q2.5), lty='dotted', color=colors[1], size=line.size) + geom_line(data=trm.pol.p3, aes(x=x, y=q97.5), lty='dotted', color=colors[1], size=line.size) +
					#geom_line(data=trm.pol.p3, aes(x=x, y=q1), lty='twodash') + geom_line(data=trm.pol.p3, aes(x=x, y=q99), lty='twodash') +
					#geom_line(data=trm.pol.p3, aes(x=x, y=q25), lty='dashed') + geom_line(data=trm.pol.p3, aes(x=x, y=q75), lty='dashed') +
					geom_line(data=trm.pol.p3, aes(x=x, y=q10), lty='twodash', color=colors[1], size=line.size) + geom_line(data=trm.pol.p3, aes(x=x, y=q90), lty='twodash', color=colors[1], size=line.size) +
					#geom_vline(xintercept=10.17, color=colors[2] ) +					
					labs(x='time elapsed\n(years)', y='evolutionary divergence\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,13), ylim=c(0,0.1)) +
					scale_y_continuous(breaks=seq(0, 0.2, 0.02)) + scale_x_continuous(breaks=seq(0,30,2)) + 
					theme_bw() + theme(legend.key = element_blank(), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), legend.justification=c(0,1), legend.position=c(0,1), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.2))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_notA_GA11eLINfit.pdf", sep='')
			ggsave(file=file, w=5, h=5)		
			#
			#	WITH DRUG RESISTANCE
			#
			trm.pol[, withA.long:= trm.pol[, factor(as.numeric(withA), levels=c(0,1), labels=c('No','Yes'))]]
			ggplot(trm.pol, aes(x=d_TSeqT)) +					
					geom_jitter(aes(y=BRL, colour=withA.long), size=1.2, alpha=0.5, position = position_jitter(width = .1), data=subset(trm.pol, BRL>0.003)) +
					#geom_point(aes(y=BRL), size=1.6, colour="#4DAF4A", pch=15, data=trm.pol.mtc) +					
					geom_line(data=trm.pol.p2, aes(y=BRL_p)) +
					geom_line(data=trm.pol.p3, aes(x=x, y=q2.5), lty='dotted') + geom_line(data=trm.pol.p3, aes(x=x, y=q97.5), lty='dotted') +
					#geom_line(data=trm.pol.p3, aes(x=x, y=q1), lty='twodash') + geom_line(data=trm.pol.p3, aes(x=x, y=q99), lty='twodash') +
					#geom_line(data=trm.pol.p3, aes(x=x, y=q25), lty='dashed') + geom_line(data=trm.pol.p3, aes(x=x, y=q75), lty='dashed') +
					geom_line(data=trm.pol.p3, aes(x=x, y=q10), lty='twodash') + geom_line(data=trm.pol.p3, aes(x=x, y=q90), lty='twodash') +
					labs(x='time elapsed\n(years)', y='evolutionary divergence\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,13), ylim=c(0,0.1)) +
					scale_colour_brewer(name='patient with\nmulti-drug resistance', palette='Set1') +
					scale_y_continuous(breaks=seq(0, 0.2, 0.02)) + scale_x_continuous(breaks=seq(0,30,2)) + 
					theme_bw() + theme(legend.key = element_blank(), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), legend.justification=c(0,1), legend.position=c(0,1), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.2)) +
					guides(colour=guide_legend(override.aes = list(size=4)))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_wTDR_GA11eLINfit.pdf", sep='')
			ggsave(file=file, w=6, h=6)	
			#
			#	WITH MTC cohort
			#	WITH Swedish chain
			#
			file	<- "/Users/Oliver/duke/2014_HIV_MTCNL/141105_MTCNL58_nodr_INFO_TRM.R"
			load(file)	
			file	<- '/Users/Oliver/duke/2014_HIV_TChainSwedish/141105_set2_seqs_INFO_TRM.R'			
			load(file)
			
			ggplot(trm.pol, aes(x=d_TSeqT)) +
					#geom_point(aes(y=BRL, colour="MTC cohort,\nAmsterdam"), size=2, data=trm.pol.mtc) +
					geom_point(aes(y=BRL, colour="Confirmed pairs,\nSweden"), size=3, alpha=1, data=trm.gag) +					
					geom_jitter(aes(y=BRL, colour="Confirmed pairs,\nBelgium"), size=1.2, alpha=0.5, position = position_jitter(width = .1), data=subset(trm.pol, withA==FALSE & BRL>0.003)) +										
					geom_line(data=trm.pol.p2, aes(y=BRL_p)) +
					geom_line(data=trm.pol.p3, aes(x=x, y=q2.5), lty='dotted', size=0.8) + geom_line(data=trm.pol.p3, aes(x=x, y=q97.5), lty='dotted', size=0.8) +
					#geom_line(data=trm.pol.p3, aes(x=x, y=q1), lty='twodash') + geom_line(data=trm.pol.p3, aes(x=x, y=q99), lty='twodash') +
					#geom_line(data=trm.pol.p3, aes(x=x, y=q25), lty='dashed') + geom_line(data=trm.pol.p3, aes(x=x, y=q75), lty='dashed') +
					geom_line(data=trm.pol.p3, aes(x=x, y=q10), lty='twodash', size=0.8) + geom_line(data=trm.pol.p3, aes(x=x, y=q90), lty='twodash', size=0.8) +
					labs(x='time elapsed\n(years)', y='evolutionary divergence\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,13), ylim=c(0,0.1)) +
					scale_colour_manual(name='Data set', values=c('grey50',"#F0027F")) +
					scale_y_continuous(breaks=seq(0, 0.2, 0.02)) + scale_x_continuous(breaks=seq(0,30,2)) + 
					theme_bw() + theme(legend.key = element_blank(), legend.key.height=unit(2.5,"line"), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), legend.justification=c(0,1), legend.position=c(0,1), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.2)) +
					guides(colour=guide_legend(override.aes = list(size=4)))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_wSWwMTC_GA11eLINfit.pdf", sep='')
			ggsave(file=file, w=6, h=6)	
			
			
			if(1)
			{
				#
				#	plot GAG	divergence
				#
				tmp				<- subset(trm.gag, select=c(d_SeqT, d_TSeqT, BRL))			
				trm.gag.swGA11	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(tmp), family=GA(mu.link='identity',sigma.link='identity'), n.cyc = 40)			
				
				trm.gag.p2		<- data.table(d_TSeqT=seq(0,15,0.1))	
				trm.gag.p2[, BRL_p:=predict(trm.gag.swGA11, data=tmp, newdata=as.data.frame(trm.gag.p2), type='response', se.fit=FALSE)]
				trm.gag.p3		<- gamlss.centiles.get(trm.gag.swGA11, tmp$d_TSeqT, cent = c(2.5, 97.5) )		
				setkey(trm.gag.p3, x)
				trm.gag.p3		<- unique(trm.gag.p3)			
				
				ggplot(trm.gag, aes(x=d_TSeqT)) + geom_point(aes(y=BRL, pch=MODE), size=2) + geom_line(data=trm.gag.p2, aes(y=BRL_p)) +
						geom_line(data=trm.gag.p3, aes(x=x, y=q2.5), lty='dashed') +
						geom_line(data=trm.gag.p3, aes(x=x, y=q97.5), lty='dashed') +					
						labs(x='time elapsed\n(years)', y='evolutionary divergence\n(nucleotide substitutions / site)') +
						coord_cartesian(xlim=c(0,13), ylim=c(0,0.115)) +
						scale_shape_manual(values=c(16,15), name='transmission mode') + 					
						scale_y_continuous(breaks=seq(0, 0.2, 0.02)) + scale_x_continuous(breaks=seq(0,30,2)) + 
						theme_bw() + theme(legend.key = element_blank(), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), legend.justification=c(0,1), legend.position=c(0,1), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.2)) +
						guides(pch=guide_legend(override.aes = list(size=4)))
				file	<- paste(indir, '/',"141105_set2_seqs_gag_patristic_dTS_data.pdf", sep='')
				ggsave(file=file, w=6, h=6)
			}

		}
		if(1)
		{
			file	<- "/Users/Oliver/duke/2014_HIV_MTCNL/141105_MTCNL58_nodr_INFO_TRM.R"
			load(file)
			trm.pol[, withA.long:= trm.pol[, factor(as.numeric(withA), levels=c(0,1), labels=c('No','Yes'))]]
			ggplot(trm.pol, aes(x=d_TSeqT)) +
					geom_point(aes(y=BRL, pch="MTC cohort"), size=1.4, colour="black", data=trm.pol.mtc) +
					geom_jitter(aes(y=BRL, colour=withA.long, pch="Confirmed pairs\nBelgium"), size=1.2, alpha=0.5, position = position_jitter(width = .1), data=subset(trm.pol, BRL>0.003)) +
					#geom_point(aes(y=BRL), size=1.6, colour="#4DAF4A", pch=15, data=trm.pol.mtc) +					
					geom_line(data=trm.pol.p2, aes(y=BRL_p)) +
					geom_line(data=trm.pol.p3, aes(x=x, y=q2.5), lty='dotted') + geom_line(data=trm.pol.p3, aes(x=x, y=q97.5), lty='dotted') +
					#geom_line(data=trm.pol.p3, aes(x=x, y=q1), lty='twodash') + geom_line(data=trm.pol.p3, aes(x=x, y=q99), lty='twodash') +
					#geom_line(data=trm.pol.p3, aes(x=x, y=q25), lty='dashed') + geom_line(data=trm.pol.p3, aes(x=x, y=q75), lty='dashed') +
					geom_line(data=trm.pol.p3, aes(x=x, y=q10), lty='twodash') + geom_line(data=trm.pol.p3, aes(x=x, y=q90), lty='twodash') +
					labs(x='time elapsed\n(years)', y='evolutionary divergence\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,13), ylim=c(0,0.1)) +
					scale_shape_manual(values=c(16,15), name='data set') + 
					scale_colour_brewer(name='including patient with\nmulti-drug resistance', palette='Set1') +
					scale_y_continuous(breaks=seq(0, 0.2, 0.02)) + scale_x_continuous(breaks=seq(0,30,2)) + 
					theme_bw() + theme(legend.key = element_blank(), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), legend.justification=c(0,1), legend.position=c(0,1), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.2)) +
					guides(colour=guide_legend(override.aes = list(size=4)), pch=guide_legend(override.aes = list(size=4)))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_data.pdf", sep='')
			ggsave(file=file, w=6, h=6)	
			
			Rsq(trm.pol.GA11e)	#-0.02 , 1e-5; if the points at zero are included, this is 22%
			summary(trm.pol.GA11e)	#all components significant, n=2807
		}		
		if(1)
		{
			#	
			#	EXPLORE all together, not patient A, POL
			#
			#	GAMMA MODEL
			trm.pol.nA		<- subset(trm.pol, withA==FALSE, select=c(d_SeqT, d_TSeqT, BRL))			
			#trm.pol.nA		<- rbind(trm.pol.nA, data.table(BRL=0.005, d_SeqT=seq(0.5,2,len=40), d_TSeqT=seq(0.05,2,len=40)), fill=TRUE)
			#trm.pol.nA		<- rbind(trm.pol.nA, data.table(BRL=0.001, d_SeqT=rep(0.1,40), d_TSeqT=rep(0.1,40)), fill=TRUE)			
			#trm.pol.GA32	<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=4)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=2)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity',sigma.link='identity'), n.cyc = 40)			
			trm.pol.GA11	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity',sigma.link='identity'), n.cyc = 40)			
			trm.pol.GA		<- trm.pol.GA11
			trm.pol.p2		<- data.table(d_TSeqT=seq(0,15,0.1))	
			trm.pol.p2[, BRL_p:=predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(trm.pol.p2), type='response', se.fit=FALSE)]
			trm.pol.p3		<- gamlss.centiles.get(trm.pol.GA, trm.pol.nA$d_TSeqT, cent = c(2.5, 97.5) )		
			setkey(trm.pol.p3, x)
			trm.pol.p3		<- unique(trm.pol.p3)			
			ggplot(trm.pol.nA, aes(x=d_TSeqT)) + geom_jitter(aes(y=BRL), size=1.2, alpha=0.5, position = position_jitter(width = .1), data=subset(trm.pol.nA, BRL>0.002)) + geom_line(data=trm.pol.p2, aes(y=BRL_p)) + 
					geom_ribbon(data=trm.pol.p3, aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) +
					#geom_vline(xintercept=0.5+3+2*4.1, color="#E41A1C", linetype=2) +
					labs(x='time elapsed\n(years)', y='evolutionary divergence\nbetween confirmed transmission pairs\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,15), ylim=c(0,0.08)) +
					scale_colour_brewer(name='pairs with sequences\nfrom selected patient', palette='Set1') +
					scale_y_continuous(breaks=seq(0, 0.2, 0.01)) +
					scale_x_continuous(breaks=seq(0,30,2)) + theme_bw() + theme(panel.grid.minor=element_line(colour="grey75", size=0.2), panel.grid.major=element_line(colour="grey75"))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_notA_GALINfit.pdf", sep='')
			ggsave(file=file, w=5, h=5)
			
			ggplot(trm.pol.nA, aes(x=d_TSeqT)) + geom_jitter(aes(y=BRL), size=1.2, alpha=0.5, position = position_jitter(width = .1), data=subset(trm.pol.nA, BRL>0.002)) +  
					labs(x='time elapsed\n(years)', y='evolutionary divergence\nbetween confirmed transmission pairs\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,15), ylim=c(0,0.08)) +
					scale_colour_brewer(name='pairs with sequences\nfrom selected patient', palette='Set1') +
					scale_y_continuous(breaks=seq(0, 0.2, 0.01)) +
					scale_x_continuous(breaks=seq(0,30,2)) + theme_bw() + theme(panel.grid.minor=element_line(colour="grey75", size=0.2), panel.grid.major=element_line(colour="grey75"))
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_notA_data.pdf", sep='')
			ggsave(file=file, w=5, h=5)	
			
			Rsq(trm.pol.GA11)	#-0.02 , 1e-5; if the points at zero are included, this is 22%
			summary(trm.pol.GA11)	#all components significant, n=2807
		}
		if(1)
		{
			#	
			#	all together, not patient A, ENV
			#
			trm.env.nA		<- subset(trm.env, withA==FALSE, select=c(d_SeqT, d_TSeqT, BRL))
			trm.env.nA		<- rbind(trm.env.nA, data.table(BRL=seq(-0.0075-0.01,0.0075-0.01,len=20), d_SeqT=-0.01, d_TSeqT=-0.01), fill=TRUE)
			trm.env.NO42	<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=4)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=2)'), data=as.data.frame(trm.env.nA), family=NO(mu.link='identity', sigma.link='identity'), n.cyc = 40)
			trm.env.NO		<- trm.env.NO42
			trm.env.p2		<- data.table(d_TSeqT=seq(0,15,0.1))	
			trm.env.p2[, BRL_p:=predict(trm.env.NO, data=trm.env.nA, newdata=as.data.frame(trm.env.p2), type='response', se.fit=FALSE)]
			trm.env.p3		<- gamlss.centiles.get(trm.env.NO, trm.env.nA$d_TSeqT, cent = c(2.5, 97.5) )		
			setkey(trm.env.p3, x)
			trm.env.p3		<- unique(trm.env.p3)			
			ggplot(trm.env.nA, aes(x=d_TSeqT)) + geom_jitter(aes(y=BRL), size=1.2, alpha=0.5, position = position_jitter(width = .1)) + geom_line(data=trm.env.p2, aes(y=BRL_p)) + 
					geom_ribbon(data=trm.env.p3, aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) +
					#geom_vline(xintercept=0.5+3+2*4.1, color="#E41A1C", linetype=2) +
					labs(x='time elapsed\n(years)', y='evolutionary divergence\nbetween confirmed transmission pairs\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,15), ylim=c(0,0.18)) +
					scale_colour_brewer(name='pairs with sequences\nfrom selected patient', palette='Set1') +
					scale_y_continuous(breaks=seq(0, 0.2, 0.01)) +
					scale_x_continuous(breaks=seq(0,30,2)) + theme_bw() + theme(panel.grid.minor=element_line(colour="grey75", size=0.2), panel.grid.major=element_line(colour="grey75"))
			file	<- paste(indir, '/',"140921_set7_env_patristic_dTS_notA_GAfit.pdf", sep='')
			ggsave(file=file, w=5, h=5)	
			Rsq(trm.env.NO42)	#0.895
		}
		#
		if(0)
		{
			#	final model, pol, comparing to d_SeqT
			trm.pol.nAb		<- subset(trm.pol, withA==FALSE, select=c(d_SeqT, d_TSeqT, BRL))
			set(trm.pol.nAb, trm.pol.nAb[,which(d_SeqT==0)], 'd_SeqT', 0.5)			
			trm.pol.NO42b	<- gamlss(as.formula('BRL ~ bs(d_SeqT, degree=4)'), sigma.formula=as.formula('~ bs(d_SeqT, degree=2)'), data=as.data.frame(trm.pol.nAb), family=NO(mu.link='identity', sigma.link='identity'), n.cyc = 40)
			trm.pol.p2b		<- data.table(d_SeqT=seq(0,15,0.1))
			trm.pol.p2b[, BRL_p:=predict(trm.pol.NO42b, data=trm.pol.nAb, newdata=as.data.frame(trm.pol.p2b), type='response', se.fit=FALSE)]
			trm.pol.p3b		<- gamlss.centiles.get(trm.pol.NO42b, trm.pol.nAb$d_SeqT, cent = c(2.5, 97.5) )		
			setkey(trm.pol.p3b, x)
			trm.pol.p3b		<- unique(trm.pol.p3b)			
			ggplot(trm.pol.nAb, aes(x=d_SeqT)) + geom_jitter(aes(y=BRL), size=1.2, alpha=0.5, position = position_jitter(width = .1)) + geom_line(data=trm.pol.p2b, aes(y=BRL_p)) + 
					geom_ribbon(data=trm.pol.p3b, aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) +
					#geom_vline(xintercept=0.5+3+2*4.1, color="#E41A1C", linetype=2) +
					labs(x='time elapsed\n(years)', y='evolutionary divergence\nbetween confirmed transmission pairs\n(nucleotide substitutions / site)') +
					coord_cartesian(xlim=c(0,15), ylim=c(0,0.08)) +
					scale_colour_brewer(name='pairs with sequences\nfrom selected patient', palette='Set1') +
					scale_y_continuous(breaks=seq(0, 0.2, 0.01)) +
					scale_x_continuous(breaks=seq(0,30,2)) + theme_bw() + theme(panel.grid.minor=element_line(colour="grey75", size=0.2), panel.grid.major=element_line(colour="grey75"))
			
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dSS_notA_GAfit.pdf", sep='')
			ggsave(file=file, w=9, h=6)		
			Rsq(trm.pol.NO42b)	#0.243			
		}
		#
		#	plot likelihood
		#
		if(1)
		{
			#	NORMAL model
			trm.lkl			<- as.data.table(expand.grid(BRL=seq(0.0001, 0.1, 0.0001), d_TSeqT=c(0.25, 0.5, 1, 3, 10, 15)))
			trm.lkl[, mu:= predict(trm.pol.NO42, data=trm.pol.nA, newdata=as.data.frame(subset(trm.lkl, select=d_TSeqT)), what='mu', type='link')]
			trm.lkl[, sigma:= predict(trm.pol.NO42, data=trm.pol.nA, newdata=as.data.frame(subset(trm.lkl, select=d_TSeqT)), what='sigma', type='link') ]
			#trm.lkl[, nu:= as.double(1/(1+exp(-predict(trm.pol.GA32, data=trm.pol.nA, newdata=as.data.frame(subset(trm.lkl, select=d_TSeqT)), what='nu', type='link')))) ]
			set( trm.lkl, NULL, 'lkl', trm.lkl[, dNO(BRL, mu=mu, sigma=sigma)/pNO(0,mu=mu,sigma=sigma, lower.tail=FALSE)])					
			set( trm.lkl, NULL, 'd_TSeqT', trm.lkl[, factor(d_TSeqT, levels=c(0.25, 0.5, 1,3,10, 15), labels=c('3mo','6 mo','1 yr','3 yrs','10 yrs', '15 yrs'))])
			ggplot(trm.lkl, aes(x=BRL, y=lkl, group=d_TSeqT, colour=d_TSeqT)) + geom_line() + facet_grid(d_TSeqT~., scales='free_y') + 
					scale_x_continuous(breaks=seq(0, 0.2, 0.01)) + theme_bw() + theme(strip.background = element_blank(), strip.text = element_blank(), legend.justification=c(1,1), legend.position=c(1,1)) +
					labs(y='likelihood of direct HIV transmission', x='evolutionary divergence\nbetween sequences from confirmed transmission pairs\n(nucleotide substitutions / site)', 
							colour='cumulated time\nsince transmission\nin recipient and\nsource')
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_notA_NO42_likelihood.pdf", sep='')
			ggsave(file=file, w=5, h=5)	
			
			trm.pol.NO		<- trm.pol.NO42
			file	<- paste(indir, '/',"141105_set7_pol_GAmodel_nA_NO_INFO.R", sep='')
			save(file=file, trm.pol.NO, trm.pol.nA)			
		}
		if(1)
		{
			#	GAMMA 32 model
			trm.lkl			<- as.data.table(expand.grid(BRL=seq(0.0001, 0.1, 0.0001), d_TSeqT=c(0.25, 0.5, 1, 3, 10, 15)))
			trm.lkl[, mu:= predict(trm.pol.GA32, data=trm.pol.nA, newdata=as.data.frame(subset(trm.lkl, select=d_TSeqT)), what='mu', type='link')]
			trm.lkl[, sigma:= predict(trm.pol.GA32, data=trm.pol.nA, newdata=as.data.frame(subset(trm.lkl, select=d_TSeqT)), what='sigma', type='link') ]
			#trm.lkl[, nu:= as.double(1/(1+exp(-predict(trm.pol.GA32, data=trm.pol.nA, newdata=as.data.frame(subset(trm.lkl, select=d_TSeqT)), what='nu', type='link')))) ]
			set( trm.lkl, NULL, 'lkl', trm.lkl[, dGA(BRL, mu=mu, sigma=sigma)])
			#set( trm.lkl, NULL, 'lkl', trm.lkl[, dZAGA(BRL, mu=mu, sigma=sigma, nu=nu)])
			set( trm.lkl, NULL, 'd_TSeqT', trm.lkl[, factor(d_TSeqT, levels=c(0.25, 0.5, 1,3,10, 15), labels=c('3mo','6 mo','1 yr','3 yrs','10 yrs', '15 yrs'))])
			ggplot(trm.lkl, aes(x=BRL, y=lkl, group=d_TSeqT, colour=d_TSeqT)) + geom_line() + facet_grid(d_TSeqT~., scales='free_y') + 
					scale_x_continuous(breaks=seq(0, 0.2, 0.01)) + theme_bw() + theme(strip.background = element_blank(), strip.text = element_blank(), legend.justification=c(1,1), legend.position=c(1,1)) +
					labs(y='likelihood of direct HIV transmission', x='evolutionary divergence\nbetween sequences from confirmed transmission pairs\n(nucleotide substitutions / site)', 
							colour='cumulated time\nsince transmission\nin recipient and\nsource')
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_notA_GA32_likelihood.pdf", sep='')
			ggsave(file=file, w=5, h=5)	
			
			trm.pol.GA		<- trm.pol.GA32
			file	<- paste(indir, '/',"141105_set7_pol_GAmodel_nA_INFO.R", sep='')
			save(file=file, trm.pol.GA, trm.pol.nA)			
		}
		if(1)
		{
			#	GAMMA 11 model
			trm.lkl			<- as.data.table(expand.grid(BRL=seq(0.0001, 0.1, 0.0001), d_TSeqT=c(0.25, 0.5, 1, 3, 10, 15)))
			trm.lkl[, mu:= predict(trm.pol.GA11, data=trm.pol.nA, newdata=as.data.frame(subset(trm.lkl, select=d_TSeqT)), what='mu', type='link')]
			trm.lkl[, sigma:= predict(trm.pol.GA11, data=trm.pol.nA, newdata=as.data.frame(subset(trm.lkl, select=d_TSeqT)), what='sigma', type='link') ]			
			set( trm.lkl, NULL, 'lkl', trm.lkl[, dGA(BRL, mu=mu, sigma=sigma)])			
			set( trm.lkl, NULL, 'd_TSeqT', trm.lkl[, factor(d_TSeqT, levels=c(0.25, 0.5, 1,3,10, 15), labels=c('3mo','6 mo','1 yr','3 yrs','10 yrs', '15 yrs'))])
			ggplot(trm.lkl, aes(x=BRL, y=lkl, group=d_TSeqT, colour=d_TSeqT)) + geom_line() + facet_grid(d_TSeqT~., scales='free_y') + 
					scale_x_continuous(breaks=seq(0, 0.2, 0.01)) + theme_bw() + theme(strip.background = element_blank(), strip.text = element_blank(), legend.justification=c(1,1), legend.position=c(1,1)) +
					labs(y='likelihood of direct HIV transmission', x='evolutionary divergence\nbetween sequences from confirmed transmission pairs\n(nucleotide substitutions / site)', 
							colour='cumulated time\nsince transmission\nin recipient and\nsource')
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_notA_GA11_likelihood.pdf", sep='')
			ggsave(file=file, w=5, h=5)	
			
			trm.pol.GA		<- trm.pol.GA11
			file	<- paste(indir, '/',"141105_set7_pol_GA11model_nA_INFO.R", sep='')
			save(file=file, trm.pol.GA, trm.pol.nA)			
		}
		if(1)
		{
			#	GAMMA 11e model
			trm.lkl			<- as.data.table(expand.grid(BRL=seq(0.0001, 0.1, 0.0001), d_TSeqT=c(0.25, 0.5, 1, 3, 10, 15)))
			trm.lkl[, mu:= predict(trm.pol.GA11e, data=trm.pol.nA, newdata=as.data.frame(subset(trm.lkl, select=d_TSeqT)), what='mu', type='link')]
			trm.lkl[, sigma:= predict(trm.pol.GA11e, data=trm.pol.nA, newdata=as.data.frame(subset(trm.lkl, select=d_TSeqT)), what='sigma', type='link') ]			
			set( trm.lkl, NULL, 'lkl', trm.lkl[, dGA(BRL, mu=mu, sigma=sigma)])			
			set( trm.lkl, NULL, 'd_TSeqT', trm.lkl[, factor(d_TSeqT, levels=c(0.25, 0.5, 1,3,10, 15), labels=c('3mo','6 mo','1 yr','3 yrs','10 yrs', '15 yrs'))])
			ggplot(trm.lkl, aes(x=BRL, y=lkl, group=d_TSeqT, colour=d_TSeqT)) + geom_line() + facet_grid(d_TSeqT~., scales='free_y') + 
					scale_x_continuous(breaks=seq(0, 0.2, 0.01)) + theme_bw() + theme(strip.background = element_blank(), strip.text = element_blank(), legend.justification=c(1,1), legend.position=c(1,1)) +
					labs(y='likelihood of direct HIV transmission', x='evolutionary divergence\nbetween sequences from confirmed transmission pairs\n(nucleotide substitutions / site)', 
							colour='cumulated time\nsince transmission\nin recipient and\nsource')
			file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_notA_GA11e_likelihood.pdf", sep='')
			ggsave(file=file, w=5, h=5)	
			
			trm.pol.GA		<- trm.pol.GA11e
			file	<- paste(indir, '/',"141105_set7_pol_GA11emodel_nA_INFO.R", sep='')
			save(file=file, trm.pol.GA, trm.pol.nA)			
		}
		#
		#	reduce to sampling years closest to transmission
		#
		tmp		<- subset(trm.pol, select=c(FROM, TO, TIP1, TIP2, FROM_SeqT, TO_SeqT, TIME_TR))
		tmp		<- tmp[, 	{
					z1	<- FROM_SeqT[ which.min(abs(FROM_SeqT-TIME_TR)) ]
					z2	<- TO_SeqT[ which.min(abs(TO_SeqT-TIME_TR)) ]
					z	<- FROM_SeqT==z1 & TO_SeqT==z2
					list(TIP1=TIP1[z], TIP2=TIP2[z])
				}, by=c('FROM','TO')]
		#
		#	plot close
		#	divergence		
		tclose.pol	<- merge(trm.pol, subset(tmp, select=c(TIP1, TIP2)), c('TIP1','TIP2'))
		ggplot(tclose.pol, aes(x=BRL, fill=Tb4S)) + geom_histogram(binwidth=0.0005) +
				theme(legend.position='bottom') +
				labs(fill='Treatment start before\nsequence sampling date', x='patristic distance\n(estimated subst/site)', y='sequence pairs of transmission pairs\n(number)')
		tclose.pol[, d_SeqT:=abs(TO_SeqT-FROM_SeqT)]
		tclose.pol[, d_TSeqT:= abs(TO_SeqT-TIME_TR) + abs(FROM_SeqT-TIME_TR)]
		ggplot(tclose.pol, aes(x=d_TSeqT, y=BRL, colour=Tb4S)) + geom_point() + stat_smooth(method = "lm", formula= y ~ ns(x,3)) + facet_grid(.~Tb4S)
		file	<- paste(indir, '/',"140921_set7_pol_patristic_dTS_closestyr.pdf", sep='')
		ggsave(file=file, w=12, h=6)		
		
		ggplot(tclose.pol, aes(x=d_SeqT, y=BRL, colour=Tb4S)) + geom_point() + stat_smooth(method = "lm", formula= y ~ ns(x,3)) + facet_grid(.~Tb4S)
		file	<- paste(indir, '/',"140921_set7_pol_patristic_dSS_closestyr.pdf", sep='')
		ggsave(file=file, w=12, h=6)		
		
	}
}
######################################################################################
project.hivc.test<- function()
{
	require(ape)
	if(1)
	{
		x<- as.DNAbin( c("a","c","g","t","m","r","w","s","y","k","v","h","d","b","n") )
		x<- as.DNAbin( as.matrix(c("a","c","g","t","m","r","w","s","y","k","v","h","d","b","n"),1,15) )
		print(x)
		.C("hivc_printdna", x, length(x) ) 
	}
	if(0)
	{
		x<- as.DNAbin( c("a","c","g","t","m","r","w","s","y","k","v","h","d","b","n") )
		y<- as.DNAbin( c("s","c","g","t","m","r","w","s","y","k","v","h","d","b","n") )
		print(as.character(x))
		print(as.character(y))
		tmp<- 0
		gd<- .C("hivc_dist_ambiguous_dna", x, y, length(x), tmp )[[4]]
		print(gd)
	}
	quit("no")
}

