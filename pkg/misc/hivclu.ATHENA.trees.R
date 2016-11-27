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
project.examl.ATHENA1610.161102.B.isolate.clade.alignments<- function()
{
	#1- run FastTree on subtype B alignment
	#2- select two large clades in FigTree
	#3- now split alignment
	infile.c1	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessing_sequences/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B_clade1.nexus"
	infile.c2	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessing_sequences/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B_clade2.nexus"
	infile.B	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B.fasta"
	infile.D	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_D.fasta"
	
	seq.b		<- read.dna(infile, format='fa')
	ph.c1		<- read.nexus(infile.c1)
	ph.c2		<- read.nexus(infile.c2)
	
	seq.b1		<- seq.b[gsub("'",'',ph.c1$tip.label),]
	seq.b2		<- seq.b[gsub("'",'',ph.c2$tip.label),]
	seq.b3		<- seq.b[setdiff( rownames(seq.b), c(rownames(seq.b1),rownames(seq.b2)) ),]
	stopifnot( nrow(seq.b1)+nrow(seq.b2)+nrow(seq.b3)==nrow(seq.b) )
	
	#4- add D outgroup sequence
	seq.D		<- read.dna(infile.D, format='fa')
	seq.b1		<- rbind(seq.b1,seq.D["LANL.D.BR.1996.patient_96BRRJ100.DQ141203",])
	seq.b2		<- rbind(seq.b2,seq.D["LANL.D.BR.1996.patient_96BRRJ100.DQ141203",])
	seq.b3		<- rbind(seq.b3,seq.D["LANL.D.BR.1996.patient_96BRRJ100.DQ141203",])
	
	#	save
	write.dna(seq.b1, gsub('B','B_sub1',infile.B), format='fasta', colsep='', nbcol=-1)
	write.dna(seq.b2, gsub('B','B_sub2',infile.B), format='fasta', colsep='', nbcol=-1)
	write.dna(seq.b3, gsub('B','B_sub3',infile.B), format='fasta', colsep='', nbcol=-1)
}
######################################################################################
project.examl.ATHENA1610.LSD.prepare.dates.161110<- function()
{
	infile.info		<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_LANL.rda"	
	outfile.dates	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_trees/ATHENA_1610_Sequences_LANL_Dates.csv"
	#
	#	prepare csv file with sampling times for all taxa
	#
	load(infile.info)
	#	process ATHENA dates
	lsdi		<- subset(ds, select=c(FASTASampleCode, PosSeqT))
	setnames(lsdi, c('FASTASampleCode','PosSeqT'), c('TAXA','DATE'))
	set(lsdi, NULL, 'DATE', lsdi[, hivc.db.Date2numeric(DATE)])
	#	process LANL dates
	tmp			<- unique(subset(db, select=TAXA_L))
	tmp[, DATE:=tmp[, sapply( strsplit(TAXA_L,'.',fixed=TRUE), function(x) x[[3]])]]
	#tmp[, table(DATE, useNA='if')]
	set(tmp, tmp[, which(DATE=='-')], 'DATE', NA_character_)
	set(tmp, NULL, 'DATE', tmp[, as.numeric(DATE)])
	set(tmp, NULL, 'TAXA_L', tmp[, gsub('_(stripped)', '', TAXA_L, fixed=TRUE)])
	set(tmp, NULL, 'TAXA_L', tmp[, paste('LANL.',TAXA_L,sep='')])
	setnames(tmp, c('TAXA_L','DATE'), c('TAXA','DATE'))
	lsdi		<- rbind(lsdi, tmp)
	#	process OUTGROUP DATES
	infile.tree	<- "/Users/Oliver/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE2_subtype_B_wOutgroup.finaltree.000"
	ph			<- read.tree(infile.tree)
	tmp			<- subset(data.table(TAXA=ph$tip.label), grepl('OUTGROUP',TAXA))
	tmp[, DATE:=tmp[,as.numeric(sapply(strsplit(TAXA, '.', fixed=1),'[[',4))]]
	lsdi		<- rbind(lsdi, tmp)
	#	duplicate HXB2
	tmp			<- lsdi[which(grepl('HXB2',TAXA)),]
	set(tmp, NULL, 'TAXA', 'HXB2')
	lsdi		<- rbind(lsdi, tmp)
	write.csv(lsdi, file=outfile.dates, row.names=FALSE)
}

######################################################################################
project.examl.ATHENA1610.LSD.process.after.outgroup.tree<- function()
{
	#	aim is to remove taxa from alignment that are odd in tree
	require(big.phylo)
	require(phytools)	
	infile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE2_subtype_B_wOutgroup.finaltree.000.lsd.date.newick"
	
	ph			<- read.tree(infile.tree)	
	ph			<- ladderize(ph)	
	
	phd			<- data.table(TAXA= ph$tip.label)
	phd[, NL:= !grepl('LANL|OUTGROUP',TAXA)]
	phd[, TIP.CLR:= 'black']
	set(phd, phd[, which(NL)], 'TIP.CLR','red')
	
	pdf(file=paste(infile.tree,'.pdf',sep=''), width=40, height=800)
	plot(ph, tip.color= phd[, TIP.CLR], cex=0.5, adj=.05)
	axisPhylo(1)
	dev.off()	
	#
	#	remove some taxa from outgroup tree -- tree no 3
	#
	infile.tree		<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE2_subtype_B_wOutgroup.finaltree.000"
	ph				<- read.tree(infile.tree)
	infile.info		<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_LANL.rda"
	load(infile.info)
	outfile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE3_subtype_B_wOutgroup.finaltree.000"
	#	remove all 'singletons' defined by very long tip branches >5%
	#	rm also all non-B COMET sequences
	phd		<- data.table(FASTASampleCode= ph$tip.label)
	phd[, BRL:=sapply( seq_len(Ntip(ph)), function(i)	ph$edge.length[ which(ph$edge[,2]==i) ] )]
	phd		<- merge(phd, ds, by="FASTASampleCode",all.x=1)
	phd		<- subset(phd, grepl('OUTGROUP',FASTASampleCode) | ((is.na(SUBTYPE_C) | SUBTYPE_C=='B') & BRL<=0.05))
	#	rm short seq and recombinant and some LANL
	phd		<- subset(phd, !FASTASampleCode%in%c("131028040783PR","2006G020","LANL.D.UG.2008.CP076303105.HQ995281", "LANL.19B.CU.2010.CUPL1081.JN000053", "LANL.19_cpx.CU.2010.CUPL1073.JN000045", "LANL.19_cpx.CU.2009.AIDS-RP28.KP688103", "LANL.19_cpx.CU.2009.MCB1.HQ108364",
					"LANL.05_DF.TG.2006.06TG.HT156.FM955743", 'LANL.06_cpx.BF.2003.ORAV203.GU207150', 'LANL.06_cpx.BF.2003.ORAV210.GU207132', 'LANL.F1.RO.2011.F1_DU2837_2011.KJ194766', 'LANL.F1.RO.-.SR7.JX966533',
					'LANL.01_AE.TH.-.045133.AF191195','LANL.01_AE.TH.-.GEN5WK4.AF191204','LANL.01_AE.TH.-.GEN5WK1.AF191203',
					'LANL.F1.RO.2012.F1_ND1088_2012.KJ194679','LANL.F1.ES.2010.X2922_3s_nfl.KJ883142','LANL.F1.ES.2009.X2687_1.GU326171','LANL.F1.ES.2010.X3079_2i_nfl.KJ883146'))
	tmp		<- setdiff(ph$tip.label, phd[, FASTASampleCode])
	ph		<- drop.tip(ph, tmp)
	write.tree(ph, file=outfile.tree)
	#
	#	remove some taxa from outgroup tree -- tree no 4
	#
	infile.tree		<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE2_subtype_B_wOutgroup.finaltree.000"
	ph				<- read.tree(infile.tree)
	infile.info		<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_LANL.rda"
	load(infile.info)
	outfile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE4_subtype_B_wOutgroup.finaltree.000"	
	#	rm all non-B COMET sequences
	phd		<- data.table(FASTASampleCode= ph$tip.label)
	phd[, BRL:=sapply( seq_len(Ntip(ph)), function(i)	ph$edge.length[ which(ph$edge[,2]==i) ] )]
	phd		<- merge(phd, ds, by="FASTASampleCode",all.x=1)
	phd		<- subset(phd, grepl('OUTGROUP',FASTASampleCode) | grepl('LANL',FASTASampleCode) | SUBTYPE_C=='B')
	#	rm short seq and recombinant and some LANL
	phd		<- subset(phd, !FASTASampleCode%in%c("131028040783PR","2006G020","LANL.D.UG.2008.CP076303105.HQ995281", "LANL.19B.CU.2010.CUPL1081.JN000053", "LANL.19_cpx.CU.2010.CUPL1073.JN000045", "LANL.19_cpx.CU.2009.AIDS-RP28.KP688103", "LANL.19_cpx.CU.2009.MCB1.HQ108364",
					"LANL.05_DF.TG.2006.06TG.HT156.FM955743", 'LANL.06_cpx.BF.2003.ORAV203.GU207150', 'LANL.06_cpx.BF.2003.ORAV210.GU207132', 'LANL.F1.RO.2011.F1_DU2837_2011.KJ194766', 'LANL.F1.RO.-.SR7.JX966533',
					'LANL.01_AE.TH.-.045133.AF191195','LANL.01_AE.TH.-.GEN5WK4.AF191204','LANL.01_AE.TH.-.GEN5WK1.AF191203',
					'LANL.F1.RO.2012.F1_ND1088_2012.KJ194679','LANL.F1.ES.2010.X2922_3s_nfl.KJ883142','LANL.F1.ES.2009.X2687_1.GU326171','LANL.F1.ES.2010.X3079_2i_nfl.KJ883146'))
	tmp		<- setdiff(ph$tip.label, phd[, FASTASampleCode])
	ph		<- drop.tip(ph, tmp)
	write.tree(ph, file=outfile.tree)
	
}

######################################################################################
project.examl.ATHENA1610.LSD.process.after.first.tree<- function()
{
	#	aim is to remove taxa from alignment that are odd in tree
	require(big.phylo)
	require(phytools)	
	infile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE_subtype_B.finaltree.000.lsd.date.newick"
	infile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE2_subtype_B.finaltree.000.lsd_ra.date.newick"
	infile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE3_subtype_B.finaltree.000.lsd_ra.date.newick"
	infile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE4_subtype_B.finaltree.000.lsd_ra.date.newick"
	ph			<- read.tree(infile.tree)	
	ph			<- ladderize(ph)	
	
	phd			<- data.table(TAXA= ph$tip.label)
	phd[, NL:= !grepl('LANL',TAXA)]
	phd[, TIP.CLR:= 'black']
	set(phd, phd[, which(NL)], 'TIP.CLR','red')
	
	pdf(file=paste(infile.tree,'.pdf',sep=''), width=40, height=800)
	plot(ph, tip.color= phd[, TIP.CLR], cex=0.5, adj=.05)
	axisPhylo(1)
	dev.off()	
	#
	#
	#
	infile.info		<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_LANL.rda"
	load(infile.info)
	dodd	<- data.table(FASTASampleCode= c("M4546917032015", "M4579808092015",  
					"M2567209052006", "M2567228072009","M2567209012007", "M2567203072007", "21412569"))
	dodd	<- merge(ds, dodd, by='FASTASampleCode', all.y=1)	
	
	infile.tree		<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE_subtype_B.finaltree.000"
	outfile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE2_subtype_B.finaltree.000"
	ph				<- read.tree(infile.tree)
	ph				<- drop.tip(ph, dodd[, FASTASampleCode])
	write.tree(ph, file=outfile.tree)
	#
	#	remove all 'singletons' defined by very long tip branches >10%
	#
	infile.tree		<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE2_subtype_B.finaltree.000"
	ph				<- read.tree(infile.tree)	
	ph				<- ladderize(ph)	
	phd				<- data.table(TAXA= ph$tip.label)
	phd[, BRL:=sapply( seq_len(Ntip(ph)), function(i)	ph$edge.length[ which(ph$edge[,2]==i) ] )]
	ph				<- drop.tip(ph, subset(phd, BRL>0.1)[, TAXA])
	outfile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE3_subtype_B.finaltree.000"
	write.tree(ph, file=outfile.tree)
	#
	#	remove all 'singletons' defined by very long tip branches >7.5%
	#
	infile.tree		<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE3_subtype_B.finaltree.000"
	ph				<- read.tree(infile.tree)	
	ph				<- ladderize(ph)	
	phd				<- data.table(TAXA= ph$tip.label)
	phd[, BRL:=sapply( seq_len(Ntip(ph)), function(i)	ph$edge.length[ which(ph$edge[,2]==i) ] )]
	ph				<- drop.tip(ph, subset(phd, BRL>0.075)[, TAXA])
	outfile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE4_subtype_B.finaltree.000"
	write.tree(ph, file=outfile.tree)
	#
	#	remove all 'singletons' defined by very long tip branches >5%
	#
	infile.tree		<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE3_subtype_B.finaltree.000"
	ph				<- read.tree(infile.tree)	
	ph				<- ladderize(ph)	
	phd				<- data.table(TAXA= ph$tip.label)
	phd[, BRL:=sapply( seq_len(Ntip(ph)), function(i)	ph$edge.length[ which(ph$edge[,2]==i) ] )]
	ph				<- drop.tip(ph, subset(phd, BRL>0.075)[, TAXA])
	outfile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE4_subtype_B.finaltree.000"
	write.tree(ph, file=outfile.tree)
	
}
######################################################################################
project.examl.ATHENA1610.LSD.run.161110<- function()
{
	require(big.phylo)
	#	on MAC
	if(0)
	{
		infile.dates	<- "/Users/Oliver/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_trees/ATHENA_1610_Sequences_LANL_Dates.csv"	
		indir.tree		<- "/Users/Oliver/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees"
		outdir			<- indir.tree
	}
	#	on HPC
	if(1)
	{
		infile.dates	<- "/work/or105/ATHENA_2016/data/lsd/ATHENA_1610_Sequences_LANL_Dates.csv"
		indir.tree		<- "/work/or105/ATHENA_2016/data/examl"
		outdir			<- "/work/or105/ATHENA_2016/data/lsd"		
	}
	ali.len					<- 1289			
	#lsd.args				<- '-v 2 -c -b 10 -r as'	# extremely slow ..
	lsd.args				<- '-v 2 -c -b 10'			# root at HXB2 and keep the root there	
	#lsd.args				<- '-v 2 -c -b 10 -r a'	
	#
	infile.tree			<- data.table(FT=list.files(indir.tree, pattern='^ExaML_result.*finaltree\\.[0-9]+$', full.names=TRUE))
	infile.tree[, FD:= file.path(outdir,basename(paste(FT, '.lsd.dates', sep='')))]
	infile.tree[, FT_PRUNED:= file.path(outdir,basename(gsub('\\.finaltree', '_OnlyDates.finaltree', FT)))]
	#infile.tree[, FL:= file.path(outdir,basename(paste(FT, '.lsd_ra', sep='')))]
	infile.tree[, FL:= file.path(outdir,basename(paste(FT, '.lsd', sep='')))]
	infile.tree			<- subset(infile.tree, grepl('noROGUE4', FT) & grepl('wOutgroup', FT))
	#	
	dlsd	<- infile.tree[, {
				#cmd		<- cmd.lsd.dates(infile.dates, FT, FD, run.lsd=FALSE, root=root, exclude.missing.dates=exclude.missing.dates, outfile.tree=FT_PRUNED)
				#cmd		<- paste(cmd, cmd.lsd(FT_PRUNED, FD, ali.len, outfile=FL, pr.args=lsd.args), sep='\n')
				#cmd		<- cmd.lsd.dates(infile.dates, FT, FD, run.lsd=FALSE, exclude.missing.dates=TRUE, outfile.tree=FT_PRUNED)
				#cmd		<- paste(cmd, cmd.lsd(FT_PRUNED, FD, ali.len, outfile=FL, pr.args=lsd.args), sep='\n')
				cmd			<- cmd.lsd.dates(infile.dates, FT, FD, run.lsd=FALSE, exclude.missing.dates=FALSE)
				cmd			<- paste(cmd, cmd.lsd(FT, FD, ali.len, outfile=FL, pr.args=lsd.args), sep='\n')
				list(CMD=cmd)
			}, by='FT']	
	#dlsd[1, cat(CMD)]
	dlsd[, {
				x		<- cmd.hpcwrapper(CMD, hpc.walltime=500, hpc.q="pqeelab", hpc.mem="5800mb", hpc.nproc=1)
				signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
				outfile	<- paste("lsd",signat,sep='.')
				cat(x)
				cmd.hpccaller(outdir, outfile, x)
				Sys.sleep(1)
				stop()
			}, by='FT']
}

######################################################################################
project.examl.ATHENA1610.examl.process.bstree<- function()
{
	#	aim is to remove taxa from alignment that are odd in tree
	require(big.phylo)
	require(phytools)
	infile.info	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_LANL.rda"
	infile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_trees/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE_subtype_B_examl_dtbs500.newick"	
	ph			<- read.tree(infile.tree)
	#	re-root at HXB2	
	tmp			<- which(grepl('HXB2',ph$tip.label))
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
	ph 			<- ladderize( ph )
	#	process bootstrap support values	
	ph.node.bs						<- as.numeric( ph$node.label )		
	ph.node.bs[is.na(ph.node.bs)]	<- 0
	ph.node.bs						<- ph.node.bs/100
	ph$node.label					<- ph.node.bs
	
	#
	#read patristic distances
	stat.fun						<- hivc.clu.min.transmission.cascade	
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)	
	thresh.brl						<- 0.1
	thresh.bs						<- 0.7
	#produce clustering 
	clustering						<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	
	tmp								<- paste(infile.tree,'.pdf',sep='')
	hivc.clu.plot(ph, clustering[["clu.mem"]], file=tmp, pdf.off=0, pdf.width=10, pdf.height=200, show.node.label=TRUE,  cex.edge.incluster=1, highlight.edge.of.tiplabel=c("LANL"), highlight.edge.of.tiplabel.col= c("red") )
	#hivc.clu.plot.tiplabels( seq_len(Ntip(ph)), matrix(df.seqinfo[,TXT], nrow=1,ncol=nrow(df.seqinfo)), matrix(df.seqinfo[,COL], nrow=1,ncol=nrow(df.seqinfo)), cex=0.05, adj=c(-2,0.5) )
	dev.off()
	
	
	#print(clustering)		
	hivc.clu.plot(ph, clustering[["clu.mem"]], highlight.edge.of.tiplabel=c("LANL"), highlight.edge.of.tiplabel.col= c("red") )
	#produce some tip states
	ph.tip.state<- rep(1:20, each=ceiling( length( ph$tip.label )/20 ))[seq_len(Ntip(ph))]
	states		<- data.table(state.text=1:20, state.col=rainbow(20))
	states		<- states[ph.tip.state,]
	hivc.clu.plot.tiplabels( seq_len(Ntip(ph)), matrix(states[,state.text], nrow=1,ncol=nrow(states)), matrix(states[,state.col], nrow=1,ncol=nrow(states)), cex=0.4, adj=c(-1,0.5) )
	
}

######################################################################################
project.examl.ATHENA1610.examl.process.after.second.tree<- function()
{
	#	aim is to remove taxa from alignment that are odd in tree
	require(big.phylo)
	require(phytools)
	infile.info	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_LANL.rda"
	infile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_trees/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE_subtype_B_examl_dtbs500.newick"
	outfile.path<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE_subtype_B_examl_dtbs500.nex"
	ph			<- read.tree(infile.tree)
	#
	#	get dates
	#
	load(infile.info)
	#	process ATHENA dates
	lsdi		<- subset(ds, select=c(FASTASampleCode, PosSeqT))
	setnames(lsdi, c('FASTASampleCode','PosSeqT'), c('TAXA','DATE'))
	set(lsdi, NULL, 'DATE', lsdi[, hivc.db.Date2numeric(DATE)])
	#	process LANL dates
	tmp			<- unique(subset(db, select=TAXA_L))
	tmp[, DATE:=tmp[, sapply( strsplit(TAXA_L,'.',fixed=TRUE), function(x) x[[3]])]]
	#tmp[, table(DATE, useNA='if')]
	set(tmp, tmp[, which(DATE=='-')], 'DATE', NA_character_)
	set(tmp, NULL, 'DATE', tmp[, as.numeric(DATE)])
	set(tmp, NULL, 'TAXA_L', tmp[, gsub('_(stripped)', '', TAXA_L, fixed=TRUE)])
	set(tmp, NULL, 'TAXA_L', tmp[, paste('LANL.',TAXA_L,sep='')])
	setnames(tmp, c('TAXA_L','DATE'), c('TAXA','DATE'))	
	lsdi		<- rbind(lsdi, tmp)
	#
	#	add dates to taxon name
	#
	phd			<- data.table(TAXA= ph$tip.label)
	phd[ ,DUMMY:= seq_len(nrow(phd))]
	phd			<- merge(phd, lsdi, by='TAXA', all.x=1)
	phd[, NEW_TAXA:= as.character(DATE)]	
	set(phd, NULL, 'NEW_TAXA', phd[, paste(TAXA,'_',NEW_TAXA,sep='')])
	setkey(phd, DUMMY)
	ph$tip.label<- phd[, NEW_TAXA]
	ph			<- drop.tip(ph, subset(phd, is.na(DATE))[, NEW_TAXA])
	#	write in nexus for TempEst
	require(phytools)
	writeNexus(ph, file=outfile.path)
}

######################################################################################
project.examl.ATHENA1610.examl.process.after.first.tree<- function()
{
	#	aim is to remove taxa from alignment that are odd in tree
	require(big.phylo)
	require(phytools)
	infile.info	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_Sequences_LANL.rda"
	infile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessed/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B.finaltree.000"
	ph			<- read.tree(infile.tree)	
	tmp			<- which(grepl('HXB2',ph$tip.label))
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])	
	ph			<- ladderize(ph)
	
	phd			<- data.table(TAXA= ph$tip.label)
	phd[, NL:= !grepl('LANL',TAXA)]
	phd[, TIP.CLR:= 'black']
	set(phd, phd[, which(NL)], 'TIP.CLR','red')
	
	pdf(file=paste(infile.tree,'.pdf',sep=''), width=40, height=800)
	plot(ph, tip.color= phd[, TIP.CLR], cex=0.5, adj=.05)
	dev.off()	
	
	load(infile.info)
	#	odd taxa
	dodd	<- data.table(FASTASampleCode= c("21444311", #CPX??
		"M4330616072013", #long branch 
		"R07-03480", "M4213909122011","M4532929042015", #BF?	
		"M4536216042015", "M4579808092015", "M4546917032015", #long branch
		"131028040783PR", # very long branch
		"2006G020", #F1? --> there are a bunch of other LANL sequences around that which we may want to remove too
		"M4334603092013", #D ?
		"M4416814022014","21531193","150518055302","150515000202"))
	dodd	<- merge(ds, dodd, by='FASTASampleCode', all.y=1)	
	dodd	<- subset(dodd, SUBTYPE!='B'|SUBTYPE_C!='B'|SEQ_L<500)
	#	taxa to remove from alignment
	tmp		<- c( dodd[, FASTASampleCode] , c("LANL.D.UG.2008.CP076303105.HQ995281", "LANL.19B.CU.2010.CUPL1081.JN000053", "LANL.19_cpx.CU.2010.CUPL1073.JN000045", "LANL.19_cpx.CU.2009.AIDS-RP28.KP688103", "LANL.19_cpx.CU.2009.MCB1.HQ108364",
		"LANL.05_DF.TG.2006.06TG.HT156.FM955743", 'LANL.06_cpx.BF.2003.ORAV203.GU207150', 'LANL.06_cpx.BF.2003.ORAV210.GU207132', 'LANL.F1.RO.2011.F1_DU2837_2011.KJ194766', 'LANL.F1.RO.-.SR7.JX966533',
		'LANL.01_AE.TH.-.045133.AF191195','LANL.01_AE.TH.-.GEN5WK4.AF191204','LANL.01_AE.TH.-.GEN5WK1.AF191203',
		'LANL.F1.RO.2012.F1_ND1088_2012.KJ194679','LANL.F1.ES.2010.X2922_3s_nfl.KJ883142','LANL.F1.ES.2009.X2687_1.GU326171','LANL.F1.ES.2010.X3079_2i_nfl.KJ883146'))
	#
	infile		<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B.fasta'
	seq.b		<- read.dna(infile, format='fa')
	stopifnot( length(setdiff( tmp, rownames(seq.b)))==0 )	
	seq.b		<- seq.b[setdiff( rownames(seq.b), tmp ),]
	#
	write.dna(seq.b, gsub('subtype','noROGUE_subtype',infile), format='fasta', colsep='', nbcol=-1)
	
	#
	#	generate pdf for tree number 2
	#
	require(big.phylo)
	require(phytools)
	infile.info	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_Sequences_LANL.rda"
	infile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessed/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE_subtype_B.finaltree.000"
	
	ph			<- read.tree(infile.tree)		
	#ph$tip.label[ which(grepl('LANL',ph$tip.label)) ]	
	tmp			<- which(grepl('HXB2',ph$tip.label))
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])	
	ph			<- ladderize(ph)
	phd			<- data.table(TAXA= ph$tip.label)
	phd[, NL:= !grepl('LANL',TAXA)]
	phd[, TIP.CLR:= 'black']
	set(phd, phd[, which(NL)], 'TIP.CLR','red')
	
	pdf(file=paste(infile.tree,'.pdf',sep=''), width=40, height=800)
	plot(ph, tip.color= phd[, TIP.CLR], cex=0.5, adj=.05)
	dev.off()	
	
	load(infile.info)
	dodd	<- data.table(FASTASampleCode= c("21532360", "M3656124112008", "M4414807042014", #long branch
				"M4595123112015",
				"21552782","R10-24464","21537053","150515000202",
				"LANL.56_cpx.FR.2010.patient_B.KC852172","LANL.56_cpx.FR.2010.patient_C.KC852174","LANL.56_cpx.FR.2010.URF5_patient_A.JN882655",#long branch
				"R08-05623","R08-03756","2004G030","R12-15425","R08-21316","2010G216","2001G047","2008G208",	#recombinants?
				"M4546917032015","M4579808092015"))
	dodd	<- merge(ds, dodd, by='FASTASampleCode', all.y=1)	
	
	dodd	<- subset(dodd, SUBTYPE!='B'|SUBTYPE_C!='B'|is.na(Patient))
	tmp		<- dodd[, FASTASampleCode]
	#
	infile		<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B.fasta'
	seq.b		<- read.dna(infile, format='fa')
	stopifnot( length(setdiff( tmp, rownames(seq.b)))==0 )
	seq.b		<- seq.b[setdiff( ph$tip.label, tmp ),]
	#
	write.dna(seq.b, gsub('subtype','noROGUE_subtype',infile), format='fasta', colsep='', nbcol=-1)
	
	#
	#	remove same taxa from alignment that were odd in previous tree
	#
	infile.info	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_LANL.rda"	
	load(infile.info)
	#	odd taxa
	dodd	<- data.table(FASTASampleCode= c("21444311", #CPX??
					"M4330616072013", #long branch 
					"R07-03480", "M4213909122011","M4532929042015", #BF?	
					"M4536216042015", "M4579808092015", "M4546917032015", #long branch
					"131028040783PR", # very long branch
					"2006G020", #F1? --> there are a bunch of other LANL sequences around that which we may want to remove too
					"M4334603092013", #D ?
					"M4416814022014","21531193","150518055302","150515000202"))
	dodd	<- merge(ds, dodd, by='FASTASampleCode', all.y=1)	
	dodd	<- subset(dodd, SUBTYPE!='B'|SUBTYPE_C!='B'|SEQ_L<500)
	#	taxa to remove from alignment
	tmp		<- c( dodd[, FASTASampleCode] , c("LANL.D.UG.2008.CP076303105.HQ995281", "LANL.19B.CU.2010.CUPL1081.JN000053", "LANL.19_cpx.CU.2010.CUPL1073.JN000045", "LANL.19_cpx.CU.2009.AIDS-RP28.KP688103", "LANL.19_cpx.CU.2009.MCB1.HQ108364",
					"LANL.05_DF.TG.2006.06TG.HT156.FM955743", 'LANL.06_cpx.BF.2003.ORAV203.GU207150', 'LANL.06_cpx.BF.2003.ORAV210.GU207132', 'LANL.F1.RO.2011.F1_DU2837_2011.KJ194766', 'LANL.F1.RO.-.SR7.JX966533',
					'LANL.01_AE.TH.-.045133.AF191195','LANL.01_AE.TH.-.GEN5WK4.AF191204','LANL.01_AE.TH.-.GEN5WK1.AF191203',
					'LANL.F1.RO.2012.F1_ND1088_2012.KJ194679','LANL.F1.ES.2010.X2922_3s_nfl.KJ883142','LANL.F1.ES.2009.X2687_1.GU326171','LANL.F1.ES.2010.X3079_2i_nfl.KJ883146'))
	#
	infile		<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B_wOutgroup.fasta'
	seq.b		<- read.dna(infile, format='fa')
	stopifnot( length(setdiff( tmp, rownames(seq.b)))==0 )	
	seq.b		<- seq.b[setdiff( rownames(seq.b), tmp ),]
	#
	write.dna(seq.b, gsub('subtype','noROGUE_subtype',infile), format='fasta', colsep='', nbcol=-1)
	
	#
	#	read tree number 3 with outgroup
	#
	infile.info	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_Sequences_LANL.rda"
	infile.tree	<- "/Users/Oliver/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE_subtype_B_wOutgroup.finaltree.000"
	
	ph			<- read.tree(infile.tree)				
	tmp			<- which(grepl('OUTGROUP\\.D',ph$tip.label))
	tmp			<- getMRCA(ph,tmp)
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])	
	ph			<- ladderize(ph)
	phd			<- data.table(TAXA= ph$tip.label)
	phd[, NL:= !grepl('LANL|OUTGROUP',TAXA)]
	phd[, TIP.CLR:= 'black']
	set(phd, phd[, which(NL)], 'TIP.CLR','red')
	
	pdf(file=paste(infile.tree,'.pdf',sep=''), width=40, height=800)
	plot(ph, tip.color= phd[, TIP.CLR], cex=0.5, adj=.05)
	dev.off()	
	
	write.tree(ph, file=gsub('noROGUE','noROGUE2',infile.tree))
	
	load(infile.info)
	dodd	<- data.table(FASTASampleCode= c("21532360", "M3656124112008", "M4414807042014", #long branch
					"M4595123112015",
					"21552782","R10-24464","21537053","150515000202",
					"LANL.56_cpx.FR.2010.patient_B.KC852172","LANL.56_cpx.FR.2010.patient_C.KC852174","LANL.56_cpx.FR.2010.URF5_patient_A.JN882655",#long branch
					"R08-05623","R08-03756","2004G030","R12-15425","R08-21316","2010G216","2001G047","2008G208",	#recombinants?
					"M4546917032015","M4579808092015"))
	dodd	<- merge(ds, dodd, by='FASTASampleCode', all.y=1)	
	
	dodd	<- subset(dodd, SUBTYPE!='B'|SUBTYPE_C!='B'|is.na(Patient))
	tmp		<- dodd[, FASTASampleCode]
	#
	infile		<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B.fasta'
	seq.b		<- read.dna(infile, format='fa')
	stopifnot( length(setdiff( tmp, rownames(seq.b)))==0 )
	seq.b		<- seq.b[setdiff( ph$tip.label, tmp ),]
	#
	write.dna(seq.b, gsub('subtype','noROGUE_subtype',infile), format='fasta', colsep='', nbcol=-1)
	#
	#	read tree 4 (using latest COMET subtype assignments)
	#
	#infile.info	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_Sequences_LANL.rda"
	infile.tree	<- "~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processing_trees/ExaML_result.ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B_wOutgroup.finaltree.000"
	ph			<- read.tree(infile.tree)
	tmp			<- which(grepl('OUTGROUP\\.D',ph$tip.label))
	tmp			<- getMRCA(ph,tmp)
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])	
	ph			<- ladderize(ph)
	
	phd			<- data.table(TAXA= ph$tip.label)
	phd[, NL:= !grepl('LANL',TAXA)]
	phd[, TIP.CLR:= 'black']
	set(phd, phd[, which(NL)], 'TIP.CLR','red')
	
	pdf(file=paste(infile.tree,'.pdf',sep=''), width=40, height=800)
	plot(ph, tip.color= phd[, TIP.CLR], cex=0.5, adj=.05)
	dev.off()	

}
######################################################################################
project.examl.ATHENA1610.161102<- function()
{
	require(big.phylo)
	project.examl.ATHENA1610.examl.run.161102()
	#project.examl.ATHENA1610.LSD.run.161110()	
}
######################################################################################
project.examl.ATHENA1610.examl.inspect.trees.161102<- function()
{
	indir	<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_trees'
	tmp		<- data.table(FT=list.files(indir, pattern='ExaML_result.*finaltree\\.[0-9]+$', full.names=TRUE))
	tmp		<- tmp[, {
				ph<- read.tree(FT)
				list(NT=Ntip(ph))
			}, by='FT']
}
######################################################################################
project.examl.ATHENA1610.examl.run.161102<- function()
{
	require(big.phylo)
	if(0)	#	MAC
	{
		indir		<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_sequences_latest'
		infile		<- "ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B_wOutgroup.fasta"
		outdir		<- '~/Dropbox (Infectious Disease)/2016_ATHENA_Oct_Update/processed_trees'
	}
	if(1)	#	HPC
	{
		indir		<- '/work/or105/ATHENA_2016/data/examl'
		infile		<- "ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B_wOutgroup.fasta"
		infile		<- "ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B_wOutgroup_p.fasta"
		#infile		<- "ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE_subtype_B.fasta"
		#infile		<- "ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B_sub2.fasta"
		#infile		<- "ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B_sub3.fasta"
		outdir		<- indir
	}
	
	args.parser	<- paste("-m DNA")
	if(1)
	{
		inpartition	<- gsub('\\.fasta','_partition.txt',infile)
		cat('DNA, p12=1-1287\\3,2-1287\\3\nDNA, p3=3-1287\\3\n', file=file.path(indir,inpartition))
		args.parser	<- paste("-m DNA -q",inpartition)
	}
	#	ExaML bootstrap args
	bs.from		<- 0
	bs.to		<- 0
	bs.n		<- 500
	
	args.examl	<- "-f d -D -m GAMMA"	#	 -- this is the default that worked in 24 hours	
	cmd			<- cmd.examl.bootstrap(indir, infile, bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, opt.bootstrap.by="codon", args.parser=args.parser, args.examl=args.examl, tmpdir.prefix="examl")					
	dummy		<- lapply(cmd, function(x)
			{				
				#x		<- cmd.hpcwrapper(x, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
				#x		<- cmd.hpcwrapper(x, hpc.walltime=79, hpc.q="pqeph", hpc.mem="1800mb", hpc.nproc=1)
				x		<- cmd.hpcwrapper(x, hpc.walltime=79, hpc.q="pqeelab", hpc.mem="5800mb", hpc.nproc=8)
				signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
				outfile	<- paste("exa",signat,sep='.')
				cat(x)
				cmd.hpccaller(outdir, outfile, x)
				Sys.sleep(1)
			})		
}