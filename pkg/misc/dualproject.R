project.dual<- function()
{
	HOME		<<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA'
	#HOME		<<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA"	
	#project.dual.distances.231015()
	#project.dual.examl.231015()
	project.dualinfecions.phylotypes.pipeline.fasta.160110()
	#project.dualinfecions.phylotypes.pipeline.examl.160110()
}

project.dual.distances.231015<- function()
{
	indir		<- paste(HOME,"alignments_160110",sep='/')
	#indir		<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_151023'
	infiles		<- list.files(indir, pattern='R$')
	
	for(i in seq_along(infiles))
	{
		load( paste(indir,'/',infiles[i],sep='') )
		d		<- seq.dist(seq)
		save(d, seq, file= gsub('\\.R','_dist\\.R',paste(indir,'/',infiles[i],sep='')))
		gc()
	}		
}

project.dual.examl.231015<- function() 
{
	require(big.phylo)
	#indir		<- paste(HOME,"alignments_151023",sep='/')
	indir		<- paste(HOME,"alignments_160110",sep='/')	
	#indir		<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_151023'
	#indir		<- '~/duke/2015_various'
	#infiles		<- list.files(indir, pattern='PANGEA_HIV_n5003_Imperial_.*\\.R')	
	infiles		<- 'PANGEA_HIV_n5003_Imperial_v160110_BW.R'
	infiles		<- 'PANGEA_HIV_n5003_Imperial_v160110_UG.R'
	infiles		<- 'PANGEA_HIV_n5003_Imperial_v160110_ZA.R'
	
	for(i in seq_along(infiles))
	{
		bs.from		<- 0
		bs.to		<- 0
		bs.n		<- 500
		outdir		<- indir
		infile		<- gsub(paste('\\.R',sep=''),'',infiles[i])
		args.parser	<- "-m DNA"
		args.examl	<- "-f d -D -m GAMMA"	#	 -- this is the default that worked in 24 hours				
		cmd			<- cmd.examl.bootstrap(indir, infile, bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, opt.bootstrap.by="nucleotide", args.parser=args.parser, args.examl=args.examl, tmpdir.prefix="examl")					
		dummy		<- lapply(cmd, function(x)
				{				
					#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=33, hpc.q= NA, hpc.mem="1800mb", hpc.nproc=1)
					x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=50, hpc.q="pqeelab", hpc.mem="5500mb", hpc.nproc=1)
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("exa",signat,sep='.')
					#cat(x)
					hivc.cmd.hpccaller(outdir, outfile, x)
					Sys.sleep(1)
				})			
	}		
}

project.dual.alignments.160110<- function()
{
	outdir	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110'
	#	read info
	#file	<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	#si		<- as.data.table(read.csv(file, stringsAsFactors=FALSE))
	#setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	#set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	#setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	
	#	read global PANGEA alignment w SA seqs and split by site	
	file			<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/AfricaCentreSeqs/GlobalAln_PlusAllPANGEA.fasta'
	sq				<- read.dna(file, format='fasta')
	sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
	tmp				<- sqi[, which(duplicated(TAXA))]
	set(sqi, tmp, 'TAXA', sqi[tmp, paste(TAXA,'-R2',sep='')])
	set(sqi, NULL, 'TAXA', sqi[, gsub('_BaseFreqs','',TAXA)])
	set(sqi, NULL, 'TAXA', sqi[, gsub('-','_',TAXA)])
	setkey(sqi, DUMMY)
	rownames(sq)	<- sqi[,TAXA]	
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	sqi[, PNG:= sqi[, factor(grepl('^PG|^R[0-9]',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]		
	sqi[, SITE:= NA_character_]
	tmp				<- sqi[, which(PNG=='Y' & grepl('^PG',TAXA))]
	set(sqi, tmp, 'SITE', sqi[tmp, substring(sapply(strsplit(TAXA,'_'),'[[',2),1,2)])
	tmp				<- sqi[, which(PNG=='Y' & grepl('^R[0-9]',TAXA))]
	set(sqi, tmp, 'SITE','ZA')
	sqi[, SEQLOC:= 'LosAlamos']
	set(sqi, sqi[, which(grepl('^PG',TAXA))], 'SEQLOC','Sanger')
	set(sqi, sqi[, which(grepl('^R[0-9]',TAXA))], 'SEQLOC','AfricaCentre')	
	sqi[, PANGEA_ID:= gsub('_R[0-9]+$','',TAXA)]
	save(sqi, sq, file=paste(outdir, '/PANGEA_HIV_n5003_Imperial_v160110_GlobalAlignment.rda',sep=''))
	
	sqi				<- subset(sqi, COV>0)
	seq				<- sq[ subset(sqi, SITE=='UG' | PNG=='N')[, TAXA], ]	
	write.dna( seq, file=paste(outdir,'/PANGEA_HIV_n5003_Imperial_v160110_UG.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEA_HIV_n5003_Imperial_v160110_UG.R',sep=''))
	seq		<- sq[ subset(sqi, SITE=='BW' | PNG=='N')[, TAXA], ]
	write.dna( seq, file=paste(outdir,'/PANGEA_HIV_n5003_Imperial_v160110_BW.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEA_HIV_n5003_Imperial_v160110_BW.R',sep=''))	
	seq		<- sq[ subset(sqi, SITE=='ZA' | PNG=='N')[, TAXA], ]
	write.dna( seq, file=paste(outdir,'/PANGEA_HIV_n5003_Imperial_v160110_ZA.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEA_HIV_n5003_Imperial_v160110_ZA.R',sep=''))
}

project.dual.alignments.151023<- function()
{
	outdir	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_151113'
	#	read info
	file	<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	si		<- as.data.table(read.csv(file, stringsAsFactors=FALSE))
	setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	
	#	read global PANGEA alignment and split by site	
	file			<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_GlobalAlignment.fasta"
	file			<- '~/Dropbox (Infectious Disease)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.fasta'
	sq				<- read.dna(file, format='fasta')
	sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
	tmp				<- sqi[, which(duplicated(TAXA))]
	set(sqi, tmp, 'TAXA', sqi[tmp, paste(TAXA,'-R2',sep='')])
	setkey(sqi, DUMMY)
	rownames(sq)	<- sqi[,TAXA]	
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]		
	sqi[, SITE:= NA_character_]
	tmp				<- sqi[, which(PNG=='Y')]
	set(sqi, tmp, 'SITE', sqi[tmp, substring(sapply(strsplit(TAXA,'-'),'[[',2),1,2)])
	sqi[, PANGEA_ID:= gsub('-R[0-9]+','',TAXA)]
	save(sqi, sq, si, file=paste(outdir, '/', gsub('fasta','rda',basename(file)),sep=''))
	
	sqi				<- subset(sqi, COV>0)
	sqi				<- merge(sqi, si, by=c('PANGEA_ID','COV'), all.x=1)
	seq				<- sq[ subset(sqi, SITE=='UG' | PNG=='N')[, TAXA], ]	
	write.dna( seq, file=paste(outdir,'/PANGEAconsensuses_2015-09_Imperial_UG_151113.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEAconsensuses_2015-09_Imperial_UG_151113.R',sep=''))
	seq		<- sq[ subset(sqi, SITE=='BW' | PNG=='N')[, TAXA], ]
	write.dna( seq, file=paste(outdir,'/PANGEAconsensuses_2015-09_Imperial_BW_151113.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEAconsensuses_2015-09_Imperial_BW_151113.R',sep=''))
	
	
	#	read contig alignment and split by site
	file	<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAcontigs_2015-09_Imperial/contigs_cnsalign_PNGIDn3366_CNTGSn6120_stripped99.fasta"
	cr		<- read.dna(file, format='fasta')
	cri		<- data.table(TAXA=rownames(cr))
	cri[, PNG:= cri[, factor(grepl('^[0-9]+_[0-9]+_[0-9]+.*',gsub('^\\.', '', TAXA)),levels=c(TRUE,FALSE),labels=c('Y','N'))]]
	cri[, SANGER_ID:=NA_character_]
	tmp		<- cri[, which(PNG=='Y')]
	set(cri, tmp, 'SANGER_ID', cri[tmp, sapply(strsplit( gsub('^\\.', '', TAXA), '.', fixed=1),'[[',1)] )	
	tmp		<- subset(cri, PNG=='Y')[, list(TAXA=TAXA, CNTG_ID_NEW=seq_along(TAXA), CONTG_ID= gsub(SANGER_ID,'',gsub('^\\.', '', TAXA))), by='SANGER_ID']
	cri		<- merge(cri, tmp, all.x=1, by=c('SANGER_ID','TAXA'))
	setnames(sqi, 'TAXA', 'PANGEA_ID_WDUP')
	cri[, PNG:=NULL]
	cri		<- merge(cri, subset(sqi, !is.na(SANGER_ID)), by='SANGER_ID', all.x=1)
	stopifnot( nrow(subset(cri, is.na(SANGER_ID) & PNG=='Y'))==0 )
	cri[, SITE:= NA_character_]
	tmp		<- cri[, which(PNG=='Y')]
	set(cri, tmp, 'SITE', cri[tmp, substring(sapply(strsplit(PANGEA_ID,'-'),'[[',2),1,2)])		
	tmp		<- cri[, list(TAXA_NEW= ifelse( is.na(CNTG_ID_NEW), paste('Ref.',TAXA,sep=''), paste(PANGEA_ID_WDUP,'-C',CNTG_ID_NEW,sep='') )), by='TAXA']
	cri		<- merge(cri, tmp, by='TAXA')
	setkey(cri, TAXA)
	rownames(cr)	<- cri[rownames(cr),][, TAXA_NEW]
	seq		<- cr[ subset(cri, SITE=='UG' | PNG=='N')[, TAXA_NEW], ]
	write.dna( seq, file=paste(outdir,'/PANGEAcontigs_2015-09_Imperial_UG_151023.fasta',sep=''), format='fasta', colsep='', nbcol=-1)
	save( seq, file=paste(outdir,'/PANGEAcontigs_2015-09_Imperial_UG_151023.R',sep=''))
	seq		<- cr[ subset(cri, SITE=='BW' | PNG=='N')[, TAXA_NEW], ]
	write.dna( seq, file=paste(outdir,'/PANGEAcontigs_2015-09_Imperial_BW_151023.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEAcontigs_2015-09_Imperial_BW_151023.R',sep=''))
	
	#	save info on consensus and contigs
	save( sqi, cri, file=paste(outdir,'/PANGEAinfo_2015-09_Imperial.R',sep=''))
	
	#	next: distances
}

project.dualinfections.build.blastdb<- function()
{
	#	make BLAST database
	#	rm all gaps in subtype alignment
	file	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/R_code/HIV1_REF_2010_genome_DNA.fasta'
	ref		<- read.dna(file, format='fasta')
	ref		<- seq.rmgaps(ref, rm.only.col.gaps=0)
	file	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/R_code/HIV1_REF_2010_genome_DNA_nogap.fasta'
	write.dna(ref, file=file, format="fasta", append=0, colsep='', colw=80, blocksep=0)
	cmd		<- hivc.cmd.blast.makedb('/Users/Oliver/duke/2014_Gates1125_dualinfection/R_code', 'HIV1_REF_2010_genome_DNA_nogap', with.mask=0, verbose=0)
	cat ( cmd )	
}

project.dualinfections.ExaML.AllTrees<- function()
{
	#	load R
	indir	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis'
	infile	<- 'Contigs_141126_BLASTstrict_OR_linsi.fasta'
	file	<- paste(indir, '/', substr(infile, 1, nchar(infile)-5),'R',sep='')
	tmp		<- load(file)
	print(tmp)
	#	for each tree, see if contigs monophyletic. If no, flag as potential dual infection
	
}

project.dualinfections.ExaML.DataTree<- function()
{
	if(0)
	{
		#prepare ExaML
		indir	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis'
		infile	<- 'Contigs_141126_BLASTstrict_OR_linsi.fasta'
		#	load R
		file	<- paste(indir, '/', substr(infile, 1, nchar(infile)-5),'R',sep='')
		tmp		<- load(file)
		print(tmp)
		#	get BLASTstrict selection + phylo requirement of 750nt
		tmp		<- subset(seqn, BLAST_STRICT & CONTIG_LEN>750)[, ID]
		tmp		<- seqa[tmp,]
		#	add refseqs		
		file	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/R_code/HIV1_REF_2010_genome_DNA.fasta'
		ref		<- read.dna(file, format='fasta')	
		tmp		<- rbind(tmp, ref)
		#	
		seqa	<- seq.rmgaps(tmp, rm.only.col.gaps=1, rm.char='-', verbose=1)
		indir	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis'
		infile	<- 'Contigs_141126_BLASTstrict750_OR_linsi.R'
		file	<- paste(indir, '/', infile, sep='')
		save(file=file, seqa)		
	}
	if(1)
	{
		#look at ExaML data tree with bootstrap values
		indir	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis'
		#infile	<- 'Contigs_141126_BLASTstrict750_OR_linsi_bs168.newick'
		infile	<- gsub('/',':','Contigs_141126_BLASTstrict750_OR_linsi_examlbs200_Fri_Nov_28_12/59/06_2013.newick')
		
		ph			<- read.tree(paste(indir,'/',infile,sep=''))
		#	re-root at SIV outgroup	http://www.ncbi.nlm.nih.gov/nuccore/AF103818
		tmp			<- which(ph$tip.label=='Ref.CPZ.US.85.US_Marilyn.AF103818')
		ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
		#	drop irrelevant tips
		tmp			<- c('Ref.CPZ.CM.05.SIVcpzMT145.DQ373066','Ref.CPZ.CD.90.ANT.U42720','Ref.P.FR.09.RBF168.GU111555','Ref.P.CM.06.U14788.HQ179987','Ref.O.CM.91.MVP5180.L20571','Ref.O.BE.87.ANT70.L20587','Ref.O.CM.98.98CMU2901.AY169812','Ref.O.SN.99.99SE_MP1300.AJ302647','Ref.N.CM.95.YBF30.AJ006022','Ref.N.CM.02.DJO0131.AY532635','Ref.N.CM.97.YBF106.AJ271370')
		ph			<- drop.tip(ph, tmp)
		#	*** ONLY for visualization purposes ***
		ph$edge.length[1]	<- 0.31
		#	bootstrap support		
		tmp			<- as.numeric( ph$node.label ) / 100
		tmp[which(is.na(tmp))]	<- 0
		ph$node.label	<- tmp		
		ph				<- ladderize(ph)		
		#	clustering
		dist.brl	<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)		
		thresh.brl	<- quantile(dist.brl,seq(0.1,1,by=0.05))["100%"]		
		thresh.bs	<- 0.8
		thresh.brl	<- 0.3
		clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph$node.label, retval="all")
		#	plot
		file		<- paste(indir,'/',gsub('.','_',infile,fixed=1),'.pdf',sep='')
		tmp			<- hivc.clu.plot(ph, clustering[["clu.mem"]], edge.col.basic="grey70", highlight.cluster=unique(na.omit(clustering[["clu.mem"]])), highlight.cluster.col='black', cex.edge.incluster=2, show.tip.label=TRUE, file=file, pdf.width=100 )
		#
		#	get patristic distance to find contigs with large distance 
		#		
		patd		<- as.matrix( distTips(ph , method='patristic') )
		tmp			<- as.data.table(expand.grid(TIP1=rownames(patd), TIP2=colnames(patd), stringsAsFactors = FALSE))
		tmp			<- tmp[, list(BRL= patd[TIP1, TIP2]), by=c('TIP1','TIP2')]
		tmp			<- subset(tmp, TIP1!=TIP2)
		tmp			<- tmp[, {
								if(TIP1<TIP2)	z<- list(TIP1n=TIP1, TIP2n=TIP2, BRL=BRL)
								if(TIP1>TIP2)	z<- list(TIP1n=TIP2, TIP2n=TIP1, BRL=BRL)
								z
							}, by=c('TIP1','TIP2')]
		set(tmp, NULL, c('TIP1','TIP2'), NULL)
		setnames(tmp, c('TIP1n','TIP2n'), c('TIP1','TIP2'))
		setkey(tmp, TIP1, TIP2)
		patd		<- unique(tmp)
		#	add approx subtype
		tmp			<- unique( subset(seqn, BLAST_STRICT, select=c(ID, BLAST_SUBTYPE)) )
		setnames(tmp, c('ID','BLAST_SUBTYPE'), c('TIP1','ST1'))
		patd		<- merge(patd, tmp, all.x=1, by='TIP1')
		setnames(tmp, c('TIP1','ST1'), c('TIP2','ST2'))
		patd		<- merge(patd, tmp, all.x=1, by='TIP2')		
		#	get branch lengths among contigs from same PATIENT
		set(patd, NULL, 'PATIENT1', patd[, sapply(strsplit(TIP1,'|',fixed=1),'[[',1)])
		set(patd, NULL, 'PATIENT2', patd[, sapply(strsplit(TIP2,'|',fixed=1),'[[',1)])		
		patdp		<- subset(patd, PATIENT1==PATIENT2, , select=c(TIP1, TIP2, BRL, ST1, ST2) )
		setnames(patdp, c('TIP1','TIP2'), c('ID1','ID2'))
		patdp		<- merge(patdp, subset(seqd.str, BS==0), by=c('ID1','ID2'), all.x=1)
		#
		ggplot(patdp, aes(x=BRL, fill= ST1==ST2)) + geom_histogram(binwidth=0.05) + theme_bw() +labs(x='patristic distance among contigs from same patient\n(subst/site)') + theme(legend.position='bottom',panel.grid.minor=element_line(colour="grey70", size=0.4), panel.grid.major=element_line(colour="grey70", size=0.4)) +
			scale_fill_brewer(palette='Set1', name='same BLAST subtype') +
			scale_x_continuous(breaks=seq(0,2,0.2))
		file		<- paste(indir,'/',gsub('.','_',infile,fixed=1),'_ContigSamePat_PatDist.pdf',sep='')
		ggsave(file=file, w=8, h=5)
		#	patdp[, summary(BRL)]
		#
		#	patristic distance super large: median 0.21 !!
		#
		#	for now consider only those with rel large overlap
		patdp[, summary(ID12_N)]
		patdpo		<- subset( patdp, ID12_N>750 )		#based on Hamming figure
		#	add clustering
		patdpo		<- merge(patdpo, data.table(ID1= ph$tip.label, CLU1= clustering$clu.mem[seq_len(Ntip(ph))]), by='ID1')
		patdpo		<- merge(patdpo, data.table(ID2= ph$tip.label, CLU2= clustering$clu.mem[seq_len(Ntip(ph))]), by='ID2')
		print( subset(patdpo, CLU1!=CLU2 | !is.na(CLU1) & is.na(CLU2) | is.na(CLU1) & !is.na(CLU2) | is.na(CLU1) & is.na(CLU2) & BRL>0.2 ) )
		select		<- unique( subset(patdpo, CLU1!=CLU2 | !is.na(CLU1) & is.na(CLU2) | is.na(CLU1) & !is.na(CLU2) | is.na(CLU1) & is.na(CLU2) & BRL>0.2, PATIENT) )
		#	plot selected in colours for every patient		
		require(colorspace)
		set(select, NULL, 'COL', rainbow_hcl(nrow(select), start = 90, end = -30) ) 				
		select[, {
					cat(paste('\nplot patient',PATIENT))
					select.lab	<- data.table(TIP= seq_len(Ntip(ph)), ID= ph$tip.label, COL2='transparent')
					tmp			<- select.lab[, which(sapply(strsplit(ID,'|',fixed=1),'[[',1)==PATIENT)]
					set(select.lab, tmp, 'COL2', COL[1])
					setkey(select.lab, TIP)		
					#print(select.lab)
					#	plot phylogeny
					file		<- paste(indir,'/',gsub('.','_',infile,fixed=1),'_dualinfection_',PATIENT,'.pdf',sep='')
					tmp			<- hivc.clu.plot(ph, clustering[["clu.mem"]], edge.col.basic="grey70", highlight.cluster=unique(na.omit(clustering[["clu.mem"]])), highlight.cluster.col='black', cex.edge.incluster=2, show.tip.label=FALSE,										
													file=file, pdf.width=100, pdf.off=0 )		
					hivc.clu.plot.tiplabels( select.lab[, TIP], matrix(select.lab[, ID], nrow=1,ncol=nrow(select.lab)), matrix(select.lab[, COL2], nrow=1,ncol=nrow(select.lab)), cex=0.4, adj=c(-0.05,0.5) )
					dev.off()
					#	save alignment					
					tmp			<- seqa[select.lab[COL2!='transparent', ID],] 
					file		<- paste(indir,'/',gsub('.','_',infile,fixed=1),'_dualinfection_',PATIENT,'.fasta',sep='')
					write.dna(file=file, tmp, format='fasta')
					#
				}, by='PATIENT']		
	}	
}

project.dualinfections.summary<- function()
{
	indir	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis'
	infile	<- 'Contigs_141126_BLASTstrict_OR_linsi.fasta'
	#	load R
	file	<- paste(indir, '/', substr(infile, 1, nchar(infile)-5),'R',sep='')
	tmp		<- load(file)
	print(tmp)
	
	tmp		<- contigs.summarize(seq, seqn, seqa, seqd.str, outdir)
}

project.dualinfecions.copy.bam<- function()
{
	tmp		<- subset(si, grepl('^PG14-ZA', PANGEA_ID), c(PANGEA_ID, SANGER_ID))[, SANGER_ID]
	tmp		<- subset(pty.clu, !is.na(SANGER_ID))[, SANGER_ID]
	paste('zip SangerAC.zip ',paste(tmp,'_ref* ',collapse='',sep=''),paste(tmp,'.bam ',collapse='',sep=''),sep='')
}

project.dualinfecions.phylotypes.setup.ZA.160110<- function()
{
	#
	#	input args
	#
	pty.gd		<- 0.2
	pty.sel.n	<- 20
	
	infile		<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_GlobalAlignment.rda'
	load(infile)	#loads sqi, sq
	indir		<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data"
	outdir		<- indir
	infile		<-  "PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500.rda"
	load(paste(indir,infile,sep='/'))	#loads "ph" "dist.brl" "ph.gdtr"  "ph.mrca"
	#	add SANGER_ID
	infile.s	<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	si			<- as.data.table(read.csv(infile.s, stringsAsFactors=FALSE))
	setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	tmp			<- subset(si, grepl('^PG14-ZA', PANGEA_ID), c(PANGEA_ID, SANGER_ID))
	set(tmp, NULL, 'PANGEA_ID', tmp[,gsub('-','_',PANGEA_ID)])
	sqi			<- merge(sqi, tmp, by='PANGEA_ID', all.x=1)	
	#	delete duplicates identified by Tulio from tree
	infile.dup	<- '/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_ZA_Duplicates.csv'
	dupl		<- as.data.table(read.csv(infile.dup, stringsAsFactors=FALSE))
	dupl		<- subset(dupl, Duplicated==1, select=c(strains, Duplicated))
	dupl[, DUP_ID:= seq_len(nrow(dupl))]
	set(dupl, NULL, 'strains', dupl[, gsub(' +$','',gsub('^ +','',strains))])	
	dupl		<- dupl[, list(TAXA= strsplit(strains,' ')[[1]]), by='DUP_ID']
	dupl		<- merge(dupl, sqi, by='TAXA', all.x=1)
	setkey(dupl, DUP_ID)
	#write.dna( sq[dupl[, TAXA],], format='fasta', file=gsub('.csv','.fasta',infile.dup))
	dupl		<- dupl[, list(TAXA= TAXA[which.min(COV)]),by='DUP_ID']
	ph$tip.label<- gsub('-','_',ph$tip.label) 
	ph			<- drop.tip(ph, dupl[,TAXA])
	#	need to recalculate stats..
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)	
	ph.gdtr							<- cophenetic.phylo(ph)
	ph.mrca							<- mrca(ph)	
	#
	#	determine large clusters
	#
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.brl=pty.gd, dist.brl=dist.brl, retval="all")	
	pty.clu		<- subset(data.table(TX_IDX=seq_len(Ntip(ph)), TAXA= ph$tip.label, CLU_ID=clustering$clu.mem[seq_len(Ntip(ph))]), !is.na(CLU_ID))	
	#	reduce to clusters containing at least one ZA sequence
	pty.clu		<- merge(pty.clu, sqi,by='TAXA')
	tmp			<- pty.clu[, list(ANY_NOT_INCOUNTRY= any(is.na(SITE) | SITE!='ZA'), ALL_NOT_INCOUNTRY= all(is.na(SITE) | SITE!='ZA')), by='CLU_ID']
	cat('\nInspecting clusters if not in-country')
	print(tmp[, table(ANY_NOT_INCOUNTRY, ALL_NOT_INCOUNTRY)])
	tmp			<- subset(tmp, !ALL_NOT_INCOUNTRY)
	pty.clu		<- merge(pty.clu, subset(tmp, select=c(CLU_ID)), by='CLU_ID')
	#	reduce to clusters of ZA sequences
	pty.clu		<- subset(pty.clu, PNG=='Y')
	pty.clu		<- merge(pty.clu, pty.clu[, list(CLU_N= length(TAXA)), by='CLU_ID'], by='CLU_ID')
	pty.clu		<- subset(pty.clu, CLU_N>1)
	cat('\nFound in-country clusters of size:')	
	print( unique(subset(pty.clu, select=c(CLU_ID, CLU_N)))[, table(CLU_N)] )
	#
	#	get all combinations within each cluster
	#
	pty.fill	<- which(grepl('^R[0-9]+|^PG',ph$tip.label))	
	pty.runs	<- pty.clu[,{
				#CLU_ID	<- 1; TX_IDX	<- subset(pty.clu, CLU_ID==1)[, TX_IDX]
				ans	<- pty.get.taxa.combinations(length(TX_IDX), pty.sel.n)	# get all combinations within cluster
				setnames(ans, 'TX_IDX', 'IDX')
				tx	<- TX_IDX
				tmp	<- pty.sel.n-length(TX_IDX)				
				if(tmp>0)
					tx	<- c(tx, sample(setdiff(pty.fill,TX_IDX), tmp))
				ans	<- merge(ans, data.table(TX_IDX=tx, IDX=seq_along(tx), FILL= c(rep(0,length(TX_IDX)),rep(1,tmp))), by='IDX')
				subset(ans, select=c(RUN, TX_IDX, FILL))
			}, by='CLU_ID']
	setnames(pty.runs, 'RUN', 'RUN_OF_CLU')
	tmp			<- unique(subset(pty.runs, select=c(CLU_ID,RUN_OF_CLU)))
	setkey(tmp, CLU_ID, RUN_OF_CLU)
	tmp[, PTY_RUN:= seq_len(nrow(tmp))]
	pty.runs	<- merge(tmp, pty.runs,by=c('CLU_ID','RUN_OF_CLU'))
	pty.runs[, RUN_OF_CLU:=NULL]
	pty.runs[, TAXA:= ph$tip.label[ TX_IDX] ]
	#	
	cat('\nNumber of clusters=', pty.runs[, length(unique(CLU_ID))])
	cat('\nNumber of scheduled phylotype runs=', pty.runs[, max(PTY_RUN)])
	cat('\nNumber of selected taxa=', subset(pty.runs, !FILL)[, length(unique(TAXA))])
	cat('\nLargest cluster size=', unique(subset(pty.clu, select=c(CLU_ID, CLU_N)))[,  max(as.numeric(names(table(CLU_N))))] )
	tmp			<- paste(indir, '/', gsub('\\.rda','_ptyrunsinput\\.rda',infile), sep='')
	save(pty.runs, pty.clu, ph, dist.brl, ph.gdtr, ph.mrca, clustering, sqi, sq, file= tmp)
	
	
	
}

pty.cmd<- function(file.bam, file.ref, window.coord, prog=PROG.PTY, prog.raxml='raxmlHPC-AVX', prog.mafft='mafft', merge.threshold=1, min.read.count=2, quality.trim.ends=30, min.internal.quality=2, merge.paired.reads='-P',no.trees='',keep.overhangs='',out.dir='.')
{	
	stopifnot(is.character(file.bam),is.character(file.ref),is.numeric(window.coord), !length(window.coord)%%2)
	#	create local tmp dir
	cmd		<- paste("CWD=$(pwd)\n",sep='\n')
	cmd		<- paste(cmd,"echo $CWD\n",sep='')
	tmpdir	<- paste('pty','_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
	tmpdir	<- paste("$CWD/",tmpdir,sep='')
	cmd		<- paste(cmd,'mkdir -p "',tmpdir,'"\n',sep='')
	#	copy files to local tmp dir
	cmd		<- paste(cmd,'cp "',file.bam,'" "',tmpdir,'"\n',sep='')
	cmd		<- paste(cmd,'cp "',file.ref,'" "',tmpdir,'"\n',sep='')
	#	cd to tmp dir
	cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
	cmd		<- paste(cmd, prog,' ',merge.threshold,' ',min.read.count,' "',basename(file.bam),'" "',basename(file.ref),'" ',paste(as.character(window.coord), collapse=' '),sep='')	
	cmd		<- paste(cmd, '-Q1', quality.trim.ends, '-Q2',min.internal.quality, merge.paired.reads)
	if(nchar(no.trees))		
		cmd	<- paste(cmd, no.trees)	
	if(nchar(keep.overhangs))
		cmd	<- paste(cmd, keep.overhangs)
	cmd		<- paste(cmd, '--x-raxml',prog.raxml,'--x-mafft',prog.mafft,'\n')
	tmp		<- gsub('_bam.txt','',basename(file.bam))	
	cmd		<- paste(cmd, 'for file in RAxML_bestTree\\.*.tree; do\n\tmv "$file" "${file//RAxML_bestTree\\./',tmp,'_}"\ndone\n',sep='')
	cmd		<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tmv "$file" "${file//AlignedReads/',tmp,'_}"\ndone\n',sep='')
	cmd		<- paste(cmd, 'mv ',tmp,'*tree "',out.dir,'"\n',sep='')
	cmd		<- paste(cmd, 'mv ',tmp,'*fasta "',out.dir,'"\n',sep='')
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')
	cmd
}

pty.cmdwrap <- function(pty.runs, si, pty.args) 		
{
	#
	#	associate BAM and REF files with each scheduled phylotype run
	#	
	#	determine FILE_ID: pty.data stored in terms of SANGER_ID	
	tmp			<- subset(si, select=c(SANGER_ID, PANGEA_ID))
	set(tmp, NULL, 'PANGEA_ID', tmp[, gsub('-','_',PANGEA_ID)])
	setnames(tmp, c('PANGEA_ID','SANGER_ID'), c('TAXA','FILE_ID'))
	pty.runs	<- merge(pty.runs, tmp, by='TAXA', all.x=1)
	tmp			<- pty.runs[, which(is.na(FILE_ID))]
	set(pty.runs, tmp,'FILE_ID', pty.runs[tmp, TAXA])	
	#	get available Bam files
	ptyd		<- data.table(FILE=list.files(pty.args[['data.dir']]))
	ptyd[, TYPE:=NA_character_]
	set(ptyd, ptyd[, which(grepl('_ref.fasta$',FILE))], 'TYPE', 'REF')
	set(ptyd, ptyd[, which(grepl('.bam$',FILE))], 'TYPE', 'BAM')
	ptyd		<- subset(ptyd, !is.na(TYPE))
	ptyd[, FILE_ID:= gsub('\\.bam|_ref\\.fasta','',FILE)]
	ptyd		<- dcast.data.table(ptyd, FILE_ID~TYPE, value.var='FILE')
	#	merge
	pty.runs	<- merge(pty.runs, ptyd, by='FILE_ID', all.x=1)
	tmp			<- subset(pty.runs, is.na(BAM) | is.na(REF))
	if(nrow(tmp))
	{
		print(tmp)
		stop()	#check we have all BAM files		
	}
	pty.runs	<- subset(pty.runs, !is.na(BAM) & !is.na(REF)) 	
	#	determine length of each ref file
	tmp			<- unique(subset(pty.runs, select=REF))
	tmp			<- tmp[, list(REF_LEN=ncol(read.dna(file.path(pty.args[['data.dir']],REF), format='fasta'))), by='REF']
	pty.runs	<- merge(pty.runs, tmp, by='REF')
	setkey(pty.runs, PTY_RUN)
	#
	#	write pty.run files and get pty command lines
	#
	pty.win		<- pty.args[['win']]
	pty.cmd		<- pty.runs[, {
				if(0)
				{
					PTY_RUN		<- 2
					BAM			<- subset(pty.runs, PTY_RUN==2)[, BAM]
					REF			<- subset(pty.runs, PTY_RUN==2)[, REF]
					REF_LEN		<- subset(pty.runs, PTY_RUN==2)[, REF_LEN]
				}
				file.bam	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_bam.txt',sep=''))
				file.ref	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_ref.txt',sep=''))
				cat( paste(file.path(pty.args[['data.dir']],BAM),collapse='\n'), file= file.bam	)
				cat( paste(file.path(pty.args[['data.dir']],REF),collapse='\n'), file= file.ref	)				
				windows		<- seq(1,by=pty.win,len=ceiling(max(REF_LEN)/pty.win))
				windows		<- as.vector(rbind( windows,windows-1+pty.win ))
				cmd			<- pty.cmd(	file.bam, file.ref, window.coord=windows, 
										prog=pty.args[['prog']], prog.raxml=pty.args[['raxml']], prog.mafft=pty.args[['mafft']], 
										merge.threshold=pty.args[['merge.threshold']], min.read.count=pty.args[['min.read.count']], quality.trim.ends=pty.args[['quality.trim.ends']], min.internal.quality=pty.args[['min.internal.quality']], 
										merge.paired.reads=pty.args[['merge.paired.reads']], no.trees=pty.args[['no.trees']], keep.overhangs=pty.args[['keep.overhangs']],
										out.dir=pty.args[['out.dir']])
				#cat(cmd)
				list(CMD= cmd)				
			},by='PTY_RUN']
	pty.cmd
}	

project.dualinfecions.phylotypes.test<- function()
{
	tfd	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/phylotypes'
	tfd	<- '~/duke/2016_PANGEAphylotypes/phylotypes'
	tf	<- data.table(FILE=list.files(tfd, pattern='^ptyr1_*'))
	tf	<- subset(tf, !grepl('fasta2',FILE))
	tf[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]
	tf[, W_FROM:=gsub('InWindow_','',regmatches(FILE,regexpr('InWindow_[0-9]+',FILE)))] 
	tf[, W_TO:=gsub('to_','',regmatches(FILE,regexpr('to_[0-9]+',FILE)))]
	tf[, Q1:=2]
	set(tf, tf[, which(grepl('Q1-30',FILE))], 'Q1', 30)
	invisible(	tf[,{
						system(paste("sed 's/<unknown description>//' ",file.path(tfd,FILE)," > ",file.path(tfd,paste(FILE,'2',sep='')),sep=''))
					}	, by='FILE']	)
	
	
	tmp	<- tf[, list(BAM= rownames(read.dna(file.path(tfd,paste(FILE,'2',sep='')),format='fasta'))), by='FILE']
	tf	<- merge(tf, tmp, by='FILE')	
	tf[, FILE_ID:= gsub('_read.*','',BAM)]	
	tf	<- merge(unique(subset(pty.runs, select=c( TAXA ,FILE_ID))), tf, by='FILE_ID')
	tmp	<- tf[, list(BAM_N=length(BAM)), by=c('W_FROM','Q1','FILE_ID')]
	tmp	<- dcast.data.table(tmp, W_FROM+FILE_ID~Q1, value.var='BAM_N')
	setnames(tmp, c('2','30'), c('Q1_2','Q1_30'))
	subset(tmp, Q1_2<Q1_30)
	
	ggplot(tmp, aes(x=as.numeric(W_FROM),y=BAM_N,colour=FILE_ID,group=FILE_ID)) + geom_line() + facet_grid(~Q1)
	#no R8_RES486_S8_L001 at all!

	load("/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda")
}

project.dualinfecions.phylotypes.mltrees.160115<- function() 
{
	require(ggtree)
	require(phytools)
	
	load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda") )
	tmp				<- subset(si, select=c(SANGER_ID,PANGEA_ID))
	set(tmp, NULL, 'SANGER_ID', tmp[, gsub('-','_',SANGER_ID)])
	set(tmp, NULL, 'PANGEA_ID', tmp[, gsub('-','_',PANGEA_ID)])
	setnames(tmp, c('SANGER_ID','PANGEA_ID'),c('FILE_ID','TAXA'))
	pty.runs		<- merge(pty.runs, tmp, by='TAXA',all.x=1)
	tmp				<- pty.runs[, which(is.na(FILE_ID))]
	set(pty.runs, tmp, 'FILE_ID', pty.runs[tmp, TAXA])
	
	indir.tr		<- file.path(HOME,"phylotypes_160115")	
	#	collect ML tree files
	infiles		<- data.table(FILE=list.files(indir.tr, 'newick$'))
	infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]
	infiles[, W_FROM:= as.numeric(gsub('InWindow_','',regmatches(FILE,regexpr('InWindow_[0-9]+',FILE))))] 
	infiles[, W_TO:= as.numeric(gsub('to_','',regmatches(FILE,regexpr('to_[0-9]+',FILE))))]
	infiles[, IDX:=seq_len(nrow(infiles))]	
	infiles		<- subset(infiles, W_FROM<9000)
	#
	infiles[, table(PTY_RUN)]
	#	raw trees w/o any attributes
	pty.ph		<- lapply( seq_len(nrow(infiles)), function(i)
			{				
				ph			<- read.tree(file.path(indir.tr,infiles[i, FILE]))
				#	node labels
				tmp				<- ph$node.label
				tmp[which(tmp=='Root'|tmp=='')]	<- '0'
				ph$node.label	<- as.numeric(tmp)
				ph
			})
	names(pty.ph)<- infiles[, FILE]	
	#
	#	determine root for each run: find taxon with largest distance from BAM of select individuals
	#
	pty.root	<- lapply(infiles[, unique(PTY_RUN)], function(ptyr){	
			#print(ptyr)	
				#ptyr<- 14
				tmp			<- subset(infiles, PTY_RUN==ptyr)[,FILE]
				phs			<- lapply(tmp, function(x) pty.ph[[x]])	
				names(phs)	<- tmp
				#get patristic distances between target and fill
				#if no target, select target_ind with lowest taxon index (just makes an unambiguous selection when windows are considered one by one)
				phpd		<- do.call('rbind',lapply(seq_along(phs),function(i){
							#i			<- 1
							ph			<- phs[[i]]							
							phb			<- data.table(IDX=seq_along(ph$tip.label), BAM=ph$tip.label, FILE_ID= gsub('_read.*','',ph$tip.label))				
							phb			<- merge(phb, pty.runs[pty.runs$PTY_RUN==ptyr,], by='FILE_ID', all=1)
							phgd		<- cophenetic.phylo(ph)
							tmp			<- subset(phb, !FILL & !is.na(BAM))[, BAM]
							if(length(tmp)==0)
							{
								tmp		<- subset(phb, !is.na(BAM))[, FILE_ID[which.min(TX_IDX)]]
								tmp		<- subset(phb, !is.na(BAM) & FILE_ID==tmp)[, BAM]
							}								
							tmp			<- phgd[ tmp, setdiff(rownames(phgd), tmp), drop=FALSE]			
							ans			<- as.data.table(melt(tmp, value.name='PATR'))
							setnames(ans, c('Var1','Var2'), c('TARGET','FILL'))
							ans[, FILE:= names(phs)[i]]	
							ans[, IDX:=i]							
							ans							
						}))
		#print(phpd[, unique(IDX)])				
				phpd[, FILL_IND:= gsub('_read.*','',FILL)]				
		#print(phpd)		
				#	try consider as root only individual present across all BAM files
				tmp		<- phpd[, list(FILL_IND_N=length(FILL)), by=c('FILL_IND','FILE')]
				tmp		<- dcast.data.table(tmp, FILE~FILL_IND, value.var='FILL_IND_N')
				tmp		<- apply(tmp[,-1, with=FALSE], 2, function(x) !any(is.na(x))	)
				tmp		<- data.table(FILL_IND=names(tmp)[tmp])
				if(nrow(tmp)==0)	# 	may be empty
				{
					print(PTY_RUN)
					stop()
				}				
				tmp		<- merge(phpd, tmp, by='FILL_IND')
				#	select individual with average largest distance
				tmp		<- tmp[, list(PATR= median(PATR)), by='FILL_IND']	
				root	<- tmp[which.max(PATR),][,FILL_IND]
				ans		<- subset(phpd, FILL_IND==root)[, list(ROOT=FILL[which.max(PATR)]), by=c('IDX','FILE')]
				tmp		<- ans[, list(CHECK= ROOT%in%phs[[IDX]]$tip.label) , by='IDX']
				stopifnot( tmp[, all(CHECK)] )		
				ans[, IDX:=NULL]
				ans[, PTY_RUN:= ptyr]
				ans
			})
	pty.root	<- do.call('rbind',pty.root)
	infiles		<- merge(infiles,pty.root, by=c('FILE','PTY_RUN'))
	#
	#	root and ladderize trees
	#
	for(i in seq_along(pty.ph))
	{
		#print(i)
		#i<- 44
		#i<- 1
		root			<- subset(infiles, FILE==names(pty.ph)[i])[, ROOT]
		ph				<- pty.ph[[i]]
		tmp										<- which(ph$tip.label==root)
		ph										<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
		ph$node.label[ph$node.label=='Root']	<- 0
		ph$node.label							<- as.numeric(ph$node.label)			
		ph				<- ladderize(ph)
		phb				<- data.table(IDX=seq_along(ph$tip.label), BAM=ph$tip.label, FILE_ID= gsub('_read.*','',ph$tip.label))
		#	group edges by individual				
		tmp				<- lapply( phb[, unique(FILE_ID)], function(x)	subset(phb, FILE_ID==x)[, BAM]	)
		names(tmp)		<- phb[, unique(FILE_ID)]
		ph				<- groupOTU(ph, tmp, group='INDIVIDUAL')		
		z	<- merge(data.table(FROM=ph$edge[,1],IDX=ph$edge[,2]), phb, by='IDX', all=1)
		z[, GROUP:= attr(ph,'INDIVIDUAL')[1:nrow(ph$edge)]]
		z	<- unique(subset(z, !is.na(FILE_ID), select=c(FILE_ID, GROUP)))
		attr(ph,'INDIVIDUAL')	<- factor(attr(ph,'INDIVIDUAL'), levels=c(0,z[,as.character(GROUP)]), labels=c('not characterized',z[,FILE_ID]))
		#	group edges by FILL
		tmp				<- as.numeric(gsub('ptyr','',regmatches(names(pty.ph)[i], regexpr('ptyr[0-9]+',names(pty.ph)[i]))))
		phb				<- merge(phb, subset(pty.runs, PTY_RUN==tmp), by='FILE_ID')
		set(phb, NULL, 'FILL', phb[, factor(FILL, levels=c(0,1), labels=c('target','filler'))])
		tmp				<- lapply( phb[, unique(FILL)], function(x)	subset(phb, FILL==x)[, BAM]	)
		names(tmp)		<- phb[, unique(FILL)]
		ph				<- groupOTU(ph, tmp, group='TYPE')				
		z	<- merge(data.table(FROM=ph$edge[,1],IDX=ph$edge[,2]), phb, by='IDX', all=1)
		z[, GROUP:= attr(ph,'TYPE')[1:nrow(ph$edge)]]
		z	<- unique(subset(z, !is.na(FILE_ID), select=c(FILL, GROUP)))
		attr(ph,'TYPE')	<- factor(attr(ph,'TYPE'), levels=c(0,z[,as.character(GROUP)]), labels=c('not characterized',z[,as.character(FILL)]))
		pty.ph[[i]]		<- ph
	}		
	#
	#	get statistics
	#
	#	node heights
	tmp					<- infiles[,{										
										ph		<- pty.ph[[FILE]]
										tmp		<- node.depth.edgelength(ph)[1:Ntip(ph)]
										list(BAM=ph$tip.label, HEIGHT=tmp)
									}, by='FILE']
	pty.stat			<- merge(infiles, tmp, by='FILE')
	#	extract individual from read name
	pty.stat[, IND:= gsub('_read.*','',BAM)]
	#	extract number of identical reads from read name 
	pty.stat[, READ_N:= as.integer(gsub('count_','',regmatches(BAM, regexpr('count_[0-9]+',BAM))))]
	#	check if reads from same individual are monophyletic
	infiles[,{										
				ph		<- pty.ph[[FILE]]
				
				tmp		<- node.depth.edgelength(ph)[1:Ntip(ph)]
				list(BAM=ph$tip.label, HEIGHT=tmp)
			}, by='FILE']

	#	TODO change bootstrap code to return data tree too
	
	
	infiles				<- merge(infiles, pty.stat[, list(HMX= max(HEIGHT)), by='FILE'], by='FILE')
	
	
	# For each patient: record whether his/her tips are monophyletic, find the
	# pairwise patristic distances between the tips - the 'cophenetic distances' -
	# and characterise those distances. 
	dummy.p.value<-0
	
	for (i in 1:num.ids) {
		id <- ids[i]
		num.leaves <- length(patient.tips[[id]])
		num.reads<-0
		for (tip in patient.tips[[id]]) num.reads <- num.reads + as.numeric(unlist(strsplit(tip,"count_"))[2])
		if (num.leaves>0) {
			monophyletic <- as.numeric(is.monophyletic(tree, patient.tips[[id]]))
			if (num.leaves > 1) {
				subtree <- drop.tip(phy=tree,
						tip=tree$tip.label[!(tree$tip.label %in% patient.tips[[id]])])
				subtree.dist.matrix <- cophenetic(subtree)
				subtree.dist <- subtree.dist.matrix[lower.tri(subtree.dist.matrix)]
				mean.size <- mean(subtree.dist)
				variance <- var(subtree.dist)
				coeff.of.var.size <- ifelse(num.leaves > 2, sqrt(variance)/mean.size, 0)
				## Distances are not weighted for the number of reads for each tip. 
				## Corrections are included because mean and variance of distance matrices
				## include zeros on diagonal, but shouldn't.
				#subtree.dist.cf <- as.vector(subtree.dist.matrix)
				#mean.size.cf <- mean(subtree.dist.cf)/(1-1/num.leaves)
				#variance.cf <- var(subtree.dist.cf)*(1+1/num.leaves)-mean.size.cf^2/num.leaves
				#cat("CF mean & var: ", mean.size.cf, variance.cf, '\n')
				#cat("CW mean & var: ", mean.size.cw, variance.cw, '\n')
				#cat('\n')
			} else {
				mean.size <- NA
				coeff.of.var.size <- NA
			}
			root.to.tip<-0
			for (i in 1:length(subtree$tip.label)) root.to.tip<-root.to.tip+nodeheight(subtree,i)
			root.to.tip<-root.to.tip/length(subtree$tip.label)
		} else {
			monophyletic<-NA
			mean.size<-NA
			coeff.of.var.size<-NA
			root.to.tip<-NA
		}
		pat.stats <- rbind(pat.stats, c(id, window, num.leaves, num.reads, monophyletic,
						mean.size, coeff.of.var.size,root.to.tip))
	}
	
	
	
	#
	#	plot trees
	#
	require(gridExtra)
	require(colorspace)
	setkey(infiles, PTY_RUN, W_FROM)	
	for(ptyr in infiles[, unique(PTY_RUN)])
	{
		tmp			<- subset(infiles, PTY_RUN==ptyr)
		#	title
		tmp[, TITLE:=paste('run',PTY_RUN,', window [',W_FROM,',',W_TO,']',sep='')]	
		phs			<- lapply(tmp[, FILE], function(x) pty.ph[[x]])
		names(phs)	<- tmp[, TITLE] 
		#	colours	
		tmp2		<- unique(unlist(lapply(seq_along(phs), function(i)	levels(attr(phs[[i]],'INDIVIDUAL'))	)))
		col			<- c('black',rainbow_hcl(length(tmp2)-1, start = 270, end = -30, c=100, l=50))
		names(col)	<- tmp2	
		phps		<- lapply(seq_along(phs), function(i){
					max.node.height	<- tmp[i,][, HMX]
					p				<- ggtree(phs[[i]], aes(color=INDIVIDUAL, linetype=TYPE)) + 
							geom_nodepoint(size=phs[[i]]$node.label/100*3) +
							geom_text(aes(label=label), size=2,  hjust=-.25) + 
							scale_color_manual(values=col, guide = FALSE) +											 
							scale_linetype_manual(values=c('target'='solid','filler'='dotted'),guide = FALSE) +
							theme_tree2() +
							theme(legend.position="bottom") + ggplot2::xlim(0, max.node.height*1.3) +
							labs(x='subst/site', title=names(phs)[i])
					p
				})	
		file	<- file.path( indir.tr, tmp[1,gsub('newick','pdf',gsub('_InWindow_[0-9]+_to_[0-9]+','',FILE))] )
		pdf(file=file, w=120, h=40)
		tmp	 	<- matrix(seq(1, 2*ceiling(length(phps)/2)), ncol=ceiling(length(phps)/2), nrow=2)
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(tmp), ncol(tmp))))
		for(i in seq_along(phps))
			print(phps[[i]], vp = viewport(layout.pos.row=which(tmp==i, arr.ind=TRUE)[1,'row'], layout.pos.col=which(tmp==i, arr.ind=TRUE)[1,'col']))	
		dev.off()	
	}
		
}

project.dualinfecions.phylotypes.pipeline.examl.160110<- function() 
{
	require(big.phylo)
	#
	#	input args
	#	(used function project.dualinfecions.phylotypes.pipeline.fasta.160110 to create all fasta files)
	#
	indir		<- file.path(HOME,"phylotypes")
	infiles		<- data.table(FILE=list.files(indir, pattern='fasta$'))
	infiles[, PTY_RUN:= sapply(strsplit(FILE,'_'),'[[',1)]
	
	args.examl	<- "-f d -D -m GAMMA"
	bs.to		<- 49
	bs.n		<- 50
	
	exa.cmd		<- infiles[,{
				infile		<- sub("\\.[^.]*$", "", FILE) 
				cmd			<- cmd.examl.bootstrap.on.one.machine(indir, infile, bs.from=0, bs.to=bs.to, bs.n=bs.n, outdir=indir, opt.bootstrap.by="nucleotide", args.examl=args.examl, verbose=1)
				list(CMD=cmd)
			}, by=c('PTY_RUN','FILE')]
	exa.cmd		<- exa.cmd[, list(CMD=paste(CMD,collapse='\n\n',sep='')), by='PTY_RUN']
	# submit
	invisible(exa.cmd[,	{
						cmd			<- hivc.cmd.hpcwrapper(CMD, hpc.walltime=10, hpc.q="pqeelab", hpc.mem="5000mb",  hpc.nproc=1)
						#cmd		<- hivc.cmd.hpcwrapper(CMD, hpc.walltime=10, hpc.q="pqeph", hpc.mem="1800mb",  hpc.nproc=1)
						outdir		<- file.path(HOME,"ptyruns")
						outfile		<- paste("ptp",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
						hivc.cmd.hpccaller(outdir, outfile, cmd)
					}, by='PTY_RUN'])
}

project.dualinfecions.phylotypes.pipeline.fasta.160110<- function() 
{
	#
	#	input args
	#	(used function project.dualinfecions.phylotypes.setup.ZA.160110 to select bam files for one run)
	#
	load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda") )	
	pty.data.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'
	work.dir			<- '/Users/Oliver/duke/2016_PANGEAphylotypes/ptyruns'
	pty.prog			<- '/Users/Oliver/Dropbox\\ \\(Infectious\\ Disease\\)/PhylotypesCode_localcopy_160118/phylotypes.py'
	raxml				<- 'raxmlHPC-AVX'
	no.trees			<- '-T'
	#	on HPC
	if(1)
	{
		raxml			<- 'raxml'
		pty.data.dir	<- '/work/or105/PANGEA_mapout/data'
		work.dir		<- file.path(HOME,"ptyruns")
		pty.prog		<- '/work/or105/libs/phylotypes/phylotypes.py'
		no.trees		<- '-T'
		HPC.LOAD		<<- "module load intel-suite/2015.1 mpi R/3.2.0 raxml/8.2.4 mafft/7 anaconda/2.3.0 samtools"
	}
	#
	#	set up all temporary files and create bash commands
	#
	#	run 160115	window length 300
	if(0)
	{
		pty.args			<- list(	prog=pty.prog, mafft='mafft', raxml=raxml,
										data.dir=pty.data.dir, work.dir=work.dir, out.dir=file.path(HOME,"phylotypes"),
										merge.threshold=1, min.read.count=2, quality.trim.ends=30, min.internal.quality=2, merge.paired.reads='-P',no.trees=no.trees, win=300, keep.overhangs='' )		
		pty.c				<- pty.cmdwrap(pty.runs, si, pty.args)		
	}
	#	run 160118	window length 60 & Q1 18 & keep overhangs
	if(1)
	{		
		pty.args			<- list(	prog=pty.prog, mafft='mafft', raxml=raxml,
										data.dir=pty.data.dir, work.dir=work.dir, out.dir=file.path(HOME,"phylotypes"),
										merge.threshold=1, min.read.count=2, quality.trim.ends=18, min.internal.quality=2, merge.paired.reads='-P',no.trees=no.trees, win=60, keep.overhangs='--keep-overhangs' )		
		pty.c				<- pty.cmdwrap(pty.runs, si, pty.args)		
	}
	if(no.trees=='-T')
	{
		invisible(pty.c[,	{					
					cmd			<- hivc.cmd.hpcwrapper(CMD, hpc.walltime=1, hpc.q="pqeelab", hpc.mem="4800mb",  hpc.nproc=1)
					#cat(cmd)
					outdir		<- file.path(HOME,"ptyruns")
					outfile		<- paste("pty",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
					hivc.cmd.hpccaller(outdir, outfile, cmd)
				}, by='PTY_RUN'])
	}
	if(0)
	{
		#
		#	add HPC header and submit
		#
		invisible(pty.c[,	{
							#cmd		<- hivc.cmd.hpcwrapper(CMD, hpc.walltime=5, hpc.q="pqeelab", hpc.mem="5000mb",  hpc.nproc=1)
							cmd			<- hivc.cmd.hpcwrapper(CMD, hpc.walltime=1, hpc.q="pqeph", hpc.mem="1800mb",  hpc.nproc=2)
							outdir		<- file.path(HOME,"ptyruns")
							outfile		<- paste("pty",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							hivc.cmd.hpccaller(outdir, outfile, cmd)
							stop()
						}, by='PTY_RUN'])
	}	
}

pty.get.taxa.combinations<- function(pty.taxan, pty.sel.n)
{
	pty.maxn	<- ceiling( (pty.taxan-pty.sel.n) / (pty.sel.n-1) )
	pty.maxn	<- pty.maxn*(pty.sel.n-1)
	pty.idx		<- matrix( rev(seq_len(pty.maxn))+pty.sel.n, nrow=pty.sel.n-1 )
	#	for each col with lowest number x, repeat x-1 times
	tmp			<- lapply( seq_len(ncol(pty.idx)), function(j)
			{
				z	<- seq_len(min(pty.idx[,j])-1)
				z	<- t(matrix(data=z, nrow=length(z), ncol=nrow(pty.idx)+1))
				z[seq_len(nrow(pty.idx))+1,]	<- NA				
				z[seq_len(nrow(pty.idx))+1,]	<- pty.idx[,j]
				t(z)
			})
	pty.idx		<- do.call('rbind',tmp)
	pty.idx		<- rbind(matrix(seq_len(pty.sel.n),nrow=1), pty.idx)
	pty.idx		<- do.call('rbind',lapply( seq_len(nrow(pty.idx)), function(i) data.table(RUN=i, TX_IDX=pty.idx[i,])))
	pty.idx
}

project.dualinfecions.NJExaTrees.160110<- function()
{
	require(phytools)
	verbose		<- 1
	resume		<- 1
	thresh.brl	<- 0.07; 	thresh.bs<- 0.8;	 
	indir		<- paste(HOME,"data",sep='/')
	outdir		<- paste(HOME,"data",sep='/')
	#
	#	BW contigs
	#
	load('~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_151023/PANGEAinfo_2015-09_Imperial.R')
	
	infile		<- "PANGEAcontigs_2015-09_Imperial_BW_151023_dist.R"
	load( paste(indir,'/',infile,sep='') )
	tmp		<- as.character(seq)
	tx		<- data.table(	TAXA1= rownames(seq), 
							FIRST= apply( tmp, 1, function(x) which(x!='-')[1] ),
							LAST= ncol(seq)-apply( tmp, 1, function(x) which(rev(x)!='-')[1] )+1L		)
	tx[, LEN1:= LAST-FIRST+1L]
	d		<- merge(d, subset(tx, select=c(TAXA1,LEN1)), by='TAXA1')
	setnames(tx, c('TAXA1','LEN1'), c('TAXA2','LEN2'))
	d		<- merge(d, subset(tx, select=c(TAXA2,LEN2)), by='TAXA2')
	#
	#	BW ExaML consensus
	#	
	infile		<- "PANGEAconsensuses_2015-09_Imperial_BW_examlbs100_151023.newick"
	ph			<- read.tree( paste(indir,'/',infile,sep='') )	
	#	drop a few bizarre tips
	ph			<- drop.tip(ph, c('PG14-BW000365-S06313','PG14-BW000374-S06322','PG14-BW000057-S01150-R2','PG14-BW000064-S01157-R2','PG14-BW000059-S01152-R2','PG14-BW000058-S01151-R2','PG14-BW000062-S01155-R2','PG14-BW000063-S01156-R2'))
	#	reroot at SIV	
	tmp			<- which(ph$tip.label=='CPZ.CM.05.SIVcpzMT145.DQ373066')
	tmp			<- which(ph$tip.label=='CPZ.TZ.06.SIVcpzTAN13.JQ768416')
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])	
	ph 			<- ladderize( ph )	
	ph.node.bs						<- ph$node.label
	ph.node.bs[ph.node.bs=='Root']	<- NA
	ph.node.bs						<- as.numeric(ph.node.bs)
	ph.node.bs[is.na(ph.node.bs)]	<- 0
	ph$node.label					<- ph.node.bs / 100
	ph$node.label[ph$node.label>1]	<- 1	
	stat.fun						<- hivc.clu.min.transmission.cascade
	#stat.fun						<- max
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)	
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	print(clustering)		
	tip.color	<- rep('black',Ntip(ph))
	tip.color[ grepl('PG14',ph$tip.label) ]	<- 'DarkRed'
	file		<- paste(outdir,"/", gsub('\\.newick',paste('_gd',100*thresh.brl,'bs',thresh.bs*100,'\\.pdf',sep=''), infile), sep='')
	hivc.clu.plot(ph, clustering[["clu.mem"]], cex.edge.incluster=3, tip.color=tip.color, file=file, pdf.scaley=25, show.tip.label=TRUE, pdf.width=10)	
	save(ph, dist.brl, file=paste(indir,'/',gsub('\\.newick','\\.R',infile),sep=''))
	#
	#	UG consensus data tree
	#
	indir		<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data"
	infile		<- "ExaML_result.PANGEAconsensuses_2015-09_Imperial_UG_151023.finaltree.040"
	ph			<- read.tree( paste(indir,'/',infile,sep='') )	
	#	no branch lengths??
	
	#
	#	UG consensus BEST tree
	#
	infile		<- "PANGEAconsensuses_2015-09_Imperial_UG_examlbs100_151023.newick"
	ph			<- read.tree( paste(indir,'/',infile,sep='') )	
	#	drop a few bizarre tips
	#ph			<- drop.tip(ph, c('PG14-BW000365-S06313','PG14-BW000374-S06322','PG14-BW000057-S01150-R2','PG14-BW000064-S01157-R2','PG14-BW000059-S01152-R2','PG14-BW000058-S01151-R2','PG14-BW000062-S01155-R2','PG14-BW000063-S01156-R2'))
	#	reroot at SIV		
	tmp			<- which(ph$tip.label=='CPZ.TZ.06.SIVcpzTAN13.JQ768416')
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])	
	ph 			<- ladderize( ph )
	ph.node.bs						<- ph$node.label
	ph.node.bs[ph.node.bs=='Root']	<- NA
	ph.node.bs						<- as.numeric(ph.node.bs)
	ph.node.bs[is.na(ph.node.bs)]	<- 0
	ph$node.label					<- ph.node.bs / 100
	ph$node.label[ph$node.label>1]	<- 1	
	stat.fun						<- hivc.clu.min.transmission.cascade
	#stat.fun						<- max
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)	
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")	
	tip.color	<- rep('black',Ntip(ph))
	tip.color[ grepl('PG14',ph$tip.label) ]	<- 'DarkRed'
	file		<- paste(outdir,"/", gsub('\\.newick',paste('_gd',100*thresh.brl,'bs',thresh.bs*100,'\\.pdf',sep=''), infile), sep='')
	hivc.clu.plot(ph, clustering[["clu.mem"]], cex.edge.incluster=3, tip.color=tip.color, file=file, pdf.scaley=50, show.tip.label=TRUE, pdf.width=30)	
	save(ph, dist.brl, file=paste(indir,'/',gsub('\\.newick','\\.R',infile),sep=''))
	
	#
	#	ZA consensus BEST tree
	#	
	indir		<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data"
	outdir		<- indir
	infile		<-  "PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500.newick"
	ph			<- read.tree( paste(indir,'/',infile,sep='') )
	#	reroot at SIV		
	tmp			<- which(ph$tip.label=='CPZ.TZ.06.SIVcpzTAN13.JQ768416')
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])	
	ph 			<- ladderize( ph )
	ph.node.bs						<- ph$node.label
	ph.node.bs[ph.node.bs=='Root']	<- NA
	ph.node.bs						<- as.numeric(ph.node.bs)
	ph.node.bs[is.na(ph.node.bs)]	<- 0
	ph$node.label					<- ph.node.bs / 100
	ph$node.label[ph$node.label>1]	<- 1	
	stat.fun						<- hivc.clu.min.transmission.cascade
	#stat.fun						<- max
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)
	thresh.brl						<- 0.1
	thresh.bs						<- 0.7
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")	
	tip.color	<- rep('black',Ntip(ph))
	tip.color[ grepl('^PG14|^R[0-9]',ph$tip.label) ]	<- 'DarkRed'
	tmp			<- paste(outdir,"/", gsub('\\.newick',paste('_gd',100*thresh.brl,'bs',thresh.bs*100,'\\.pdf',sep=''), infile), sep='')
	invisible(hivc.clu.plot(ph, clustering[["clu.mem"]], cex.edge.incluster=3, tip.color=tip.color, file=tmp, pdf.scaley=25, show.tip.label=TRUE, pdf.width=10))
	
	ph.gdtr							<- cophenetic.phylo(ph)
	ph.mrca							<- mrca(ph)
	save(ph, dist.brl, ph.gdtr, ph.mrca, file=paste(indir,'/',gsub('\\.newick','\\.rda',infile),sep=''))
	
	#	
	#pdf(file=paste(indir,'/',gsub('\\.newick','\\.pdf',infile),sep=''), w=10, h=250)
	#plot(ph, show.node.label=TRUE, cex=0.4)
	#dev.off()
}

project.dualinfections.basic<- function()
{
	indir	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126'
	outdir	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis'

	tmp		<- contigs.read.v141126(indir)
	seq		<- tmp$seq
	seqn	<- copy(tmp$seq.names)
	set(seqn, NULL, 'CONTIG_LEN', sapply(seq, length))
	#	save contigs in single fasta file
	file	<- paste(outdir, 'Contigs_141126.fasta', sep='/')
	write.dna(seq, format='fasta', file=file )
	
	#### stop and do the blast in command line here ####	
	#   Match the contigs against a database of sequences that should be similar (in our case, HIV sequences)
	#   in order to filter out other nucleic acid sequences that were included in the amplification by mistake.
	#   We use the program blast and the Los Alamos HIV-1/SIVcpz 'subtype reference' set from 2010, consisting of 170 sequences.
	#   We download the latter to the file path/to/MyReferenceDatabase.fasta, say. The command we use is:
	#   makeblastdb -dbtype nucl -input_type fasta -in path/to/MyReferenceDatabase.fasta -out path/to/MyReferenceDatabase.db
	#   converts the reference sequences to a blast database. Then the command
	#   blastn -query path/to/AllMyContigs.fasta -db path/to/MyReferenceDatabase.db -out path/to/MyBlastMatches.csv -max_target_seqs 1 -outfmt "10 qacc sacc sseqid evalue pident qstart qend sstart send"
	#   performs the matching; -outmft are the output formatting options.
	#   If each sequence has multiple matches in the reference database, it will have multiple lines (comma separated, containing the fields just described)
	#   in the output file path/to/MyBlastMatches.csv; if it has no matches, it will have no lines in that file.	
	#### see other notes on how to do it in Francois_notebook/notes on sequence data
	indir	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis'
	infile	<- 'Contigs_141126'
	cmd		<- hivc.cmd.blast(indir, infile, '', '/Users/Oliver/duke/2014_Gates1125_dualinfection/R_code', 'HIV1_REF_2010_genome_DNA_nogap', dbsignat='', blast.task="blastn", blast.max_target_seqs=1, blast.evalue=10, blast.wordsize=11, blast.gapopen=5, blast.gapextend=2, blast.penalty=-3, blast.reward= 2, blast.dust= "no", nproc=1, verbose=1)
	cat ( cmd )	
	seqb	<- contigs.getblast(seqn, indir, infile)	
	seqb	<- subset(seqb, select=c(ID, BLAST_SUBTYPE, BLAST_ALIGNMENT_LENGTH, BLAST_IDENTITY, BLAST_EVALUE, BLAST_STRICT, BLAST_RELAXED))
		
	#	consensus subtype per patient	
	seqn	<- merge(seqn, seqb, by='ID', all.x=1)
	tmp		<- subset(seqn, !is.na(BLAST_SUBTYPE))[, list(BLAST_ALIGNMENT_LENGTH=max(BLAST_ALIGNMENT_LENGTH)), by=c('PATIENT','BLAST_SUBTYPE')]	
	tmp		<- tmp[, list( BLAST_CSUBTYPE=BLAST_SUBTYPE[ which.max(BLAST_ALIGNMENT_LENGTH) ] ), by='PATIENT']
	seqn	<- merge(seqn, tmp, by='PATIENT', all.x=1)

	#	run mafft
	contigs.mafft.run.vFrancois(seq, seqn, outdir)
	#	run mafft-linsi for each contig
	contigs.mafft.run.v141127(seq, seqn,outdir)
	
	if(0)
	{
		#   process mafft	FRANCOIS
		indir	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis'
		infile	<- 'Contigs_141126_BLASTstrict_Francois_mafft.fasta'
		seqa	<- contigs.mafft.get.vFrancois(indir, infile)
		#
		#	just count distances among contigs from same patient
		#	can also do bootstrap to get uncertainty
		#
		seqn.str<- subset( seqn, BLAST_STRICT )	
		seqn.str<- merge(seqn.str, seqn.str[, list(CONTIG_N=length(ID)), by='PATIENT'], by='PATIENT')
		seqd.str<- contigs.rawdist.bs(seqn.str, seqa, 100)
		#
		#	add prob of >= DN mutations occurring under Poisson model  
		#
		seqd.str<- contigs.probofdist.poissonmodelglobal(seqd.str, indir, infile, mu=0.01, t.elapsed= 7*2)
		
		subset(seqn, BLAST_STRICT)[, length(unique(PATIENT))]	#205
		seqd.str[, length(unique(PATIENT))]	#72
		subset(seqd.str, BS==0 & POIS_GE_DN<0.1)[, length(unique(PATIENT))]	#10
		subset(seqd.str, BS==0 & POIS_GE_DN>0.9)[, length(unique(PATIENT))]	#35
		
		#	save all to R
		file	<- paste(indir, '/', substr(infile, 1, nchar(infile)-5),'R',sep='')
		save(file=file, seq, seqn, seqb, seqa, seqd.str)		
	}
	if(1)
	{
		#   process mafft	OR
		indir	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis'
		infile	<- 'Contigs_141126_BLASTstrict_OR'		
		#	note: this does not run out of the box, go into function, line by line
		seqa	<- contigs.mafft.get.v141127(indir, infile)		
		#
		#	just count distances among contigs from same patient
		#	can also do bootstrap to get uncertainty
		#
		seqn.str<- subset( seqn, BLAST_STRICT )	
		seqn.str<- merge(seqn.str, seqn.str[, list(CONTIG_N=length(ID)), by='PATIENT'], by='PATIENT')
		seqd.str<- contigs.rawdist.bs(seqn.str, seqa, 100)
		#
		#	add prob of >= DN mutations occurring under Poisson model  
		#
		infile	<- 'Contigs_141126_BLASTstrict_OR_linsi.fasta'
		seqd.str<- contigs.probofdist.poissonmodelglobal(seqd.str, indir, infile, mu=0.01, t.elapsed= 7*2)
		
		subset(seqn, BLAST_STRICT)[, length(unique(PATIENT))]	#205
		seqd.str[, length(unique(PATIENT))]	#72
		subset(seqd.str, BS==0 & POIS_GE_DN<0.1)[, length(unique(PATIENT))]	#8
		subset(seqd.str, BS==0 & POIS_GE_DN>0.9)[, length(unique(PATIENT))]	#31
		
		#	save all to R
		file	<- paste(indir, '/', substr(infile, 1, nchar(infile)-5),'R',sep='')
		save(file=file, seq, seqn, seqb, seqa, seqd.str)		
	}	
}




contigs.summarize<- function(seq, seqn, seqb, seqa, seqd.str, outdir)
{
	require(ggplot2)
	require(plyr)
	
	seqn	<- merge(seqn, seqn[, list(CONTIG_LEN_MEDIAN=median(as.double(CONTIG_LEN))), by='COUNTRY'], by='COUNTRY')	
	seqn	<- merge(seqn, seqn[, list(CONTIG_N=length(CONTIG)), by='PATIENT'], by='PATIENT')	
	
	ggplot(seqn, aes(x=CONTIG_LEN)) + geom_histogram() + geom_vline(aes(xintercept=CONTIG_LEN_MEDIAN)) +  theme_bw() + facet_grid(COUNTRY~., margins=FALSE) + labs(x='contig length')
	file	<- paste(outdir, '/INFO_contig_length_RAW.pdf',sep='')
	ggsave(file=file, w=6, h=8)	
	ggplot(subset(seqn, BLAST_STRICT), aes(x=CONTIG_LEN)) + geom_histogram() +  theme_bw() + facet_grid(COUNTRY~., margins=FALSE) + labs(x='contig length')
	file	<- paste(outdir, '/INFO_contig_length_BLASTstrict.pdf',sep='')
	ggsave(file=file, w=6, h=8)	
	
	tmp		<- merge( subset(seqn, CONTIG_LEN>350 & BLAST_ALIGNMENT_LENGTH>100, select=c(ID, COUNTRY, BLAST_CSUBTYPE)), seqb, by='ID' )
	ggplot(tmp, aes(x=BLAST_ALIGNMENT_LENGTH)) + geom_histogram()  +  theme_bw() + facet_grid(COUNTRY~., margins=FALSE) + labs(x='match length to reference sequence', y='contigs with len>350 and match len>100\n(#)')
	file	<- paste(outdir, '/INFO_match_length.pdf',sep='')
	ggsave(file=file, w=6, h=8)	
	
	tmp		<- merge( subset(seqn, select=c(ID, COUNTRY, BLAST_CSUBTYPE, CONTIG_LEN)), seqb, by='ID' )
	ggplot(tmp, aes(x=BLAST_ALIGNMENT_LENGTH, y=BLAST_ALIGNMENT_LENGTH/CONTIG_LEN)) + geom_point()  +  theme_bw() + facet_grid(COUNTRY~., margins=FALSE) + labs(x='match length', y='proportion of contig matched') +
			coord_trans(xtrans="log10", limx=c(10, 12000)) +
			scale_x_continuous(breaks=c(20, 50, 100, 200, 1000, 5000, 1000), minor_breaks=c(60, 70, 80, 90)) +
			scale_y_continuous(breaks=seq(0,1,0.1), minor_breaks=c(0.11,0.12, 0.13, 0.14, 0.15, 0.2, 0.21, 0.22, 0.23, 0.24)) 
	file	<- paste(outdir, '/INFO_match_length_vs_propmatchted.pdf',sep='')
	ggsave(file=file, w=6, h=8)	
			
	tmp		<- unique(subset(seqn, select=c(PATIENT, CONTIG_N, COUNTRY)))
	ggplot(tmp, aes(x=CONTIG_N)) + geom_histogram() +  theme_bw() + facet_grid(COUNTRY~., margins=FALSE) + labs(x='number contigs per individual')
	file	<- paste(outdir, '/INFO_contig_number_per_individual_RAW.pdf',sep='')
	ggsave(file=file, w=6, h=8)	
	tmp		<- subset(seqn, BLAST_STRICT)[, list(CONTIG_N=length(CONTIG), COUNTRY=COUNTRY[1]), by='PATIENT']
	ggplot(tmp, aes(x=CONTIG_N)) + geom_histogram() +  theme_bw() + facet_grid(COUNTRY~., margins=FALSE) + labs(x='number contigs per individual')
	file	<- paste(outdir, '/INFO_contig_number_per_individual_BLASTstrict.pdf',sep='')
	ggsave(file=file, w=6, h=8)	
	
	
	tmp		<- seqs[, list(CONTIG_LEN_MAX=max(CONTIG_LEN), COUNTRY=COUNTRY[1]), by='PATIENT']
	ggplot(tmp, aes(x=CONTIG_LEN_MAX)) + geom_histogram() +  theme_bw() + facet_grid(COUNTRY~., margins=FALSE) + labs(x='max contig length per individual')
	file	<- paste(outdir, '/INFO_contig_maxnumber_per_individual.pdf',sep='')
	ggsave(file=file, w=6, h=8)	
	
	tmp		<- subset(seqn, BLAST_STRICT, select=c(PATIENT, BLAST_CSUBTYPE, COUNTRY))
	ggplot(unique(tmp), aes(x=BLAST_CSUBTYPE, fill=COUNTRY)) + geom_bar(position='dodge') + theme_bw() +
			scale_fill_brewer(palette='Set1') + labs(x='consensus subtype from Blast match')
	file	<- paste(outdir, '/INFO_BLAST_CSUBTYPE_BLASTstrict.pdf',sep='')
	ggsave(file=file, w=8, h=6)	
}

contigs.probofdist.poissonmodelglobal<- function(seqd.str, indir, infile, mu=0.01, t.elapsed= 7*2)
{
	seqd.str[, POIS_GE_DN:= seqd.str[, ppois(DN, mu*t.elapsed*ID12_N, lower.tail=FALSE)]]
	set(seqd.str, seqd.str[, which(ID12_N==0)], 'POIS_GE_DN', NA_real_)
	
	#	plot 	
	ggplot(seqd.str, aes(x=ID12_N, y=DR, colour='bootstrap')) + geom_point(alpha=0.5, size=1.2) + geom_point(data=subset(seqd.str, BS==0), aes(colour='data')) +
			scale_colour_manual(name='Contig alignment', values=c('black','red')) +
			theme_bw() + labs(x='overlap between contigs\n(# sites)', y='evol distance between contigs of same patient\n(Hamming)' ) +
			theme(legend.position='bottom')
	file	<- paste(indir, '/', substr(infile, 1, nchar(infile)-6),'_HammingDistance_ContigOverlap.pdf',sep='')
	ggsave(file=file, w=6, h=6)
	
	ggplot(subset(seqd.str, !is.nan(DR)), aes(x=interaction(ID1, ID2), y=DR)) + geom_point(alpha=0.5, size=1.2, aes(colour='bootstrap'), data=subset(seqd.str, !is.nan(DR) & BS>0))  +			
			geom_point(aes(colour='data'), data=subset(seqd.str, !is.nan(DR) & BS==0))  +
			scale_colour_manual(name='Contig alignment', values=c('black','red')) +
			scale_y_continuous(breaks=seq(0,1,0.1)) + scale_x_discrete(breaks=NULL) +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey70", size=0.4), panel.grid.major=element_line(colour="grey70", size=0.4), axis.text.x=element_blank(), legend.position='bottom') + labs(x='contig pair of same patient', y='evol distance between contigs of same patient\n(Hamming)')
	file	<- paste(indir, '/', substr(infile, 1, nchar(infile)-6),'_HammingDistance_ContigPair.pdf',sep='')
	ggsave(file=file, w=10, h=6)
	
	ggplot(subset(seqd.str, !is.nan(DR)), aes(x=interaction(ID1, ID2), y=POIS_GE_DN)) + geom_boxplot(outlier.shape=NA, fill='grey70',aes(colour='bootstrap'), data=subset(seqd.str, !is.nan(DR) & BS>0))  +			
			geom_point(aes(colour='data'), data=subset(seqd.str, !is.nan(DR) & BS==0))  +
			scale_colour_manual(name='Contig alignment', values=c('black','red')) +
			scale_y_continuous(breaks=seq(0,1,0.1)) + scale_x_discrete(breaks=NULL) + 
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey70", size=0.4), panel.grid.major=element_line(colour="grey70", size=0.4), axis.text.x=element_blank(), legend.position='bottom') + labs(x='contig pair of same patient', y='prob evol distance accumulated within patient\n(Poisson model)')
	file	<- paste(indir, '/', substr(infile, 1, nchar(infile)-6),'_ProbPoisson_ContigPair.pdf',sep='')
	ggsave(file=file, w=10, h=6)
	
	ggplot(subset(seqd.str, !is.nan(DR)), aes(x=ID12_N, y=POIS_GE_DN)) + geom_point(alpha=0.5, size=1.2,aes(colour='bootstrap'), data=subset(seqd.str, !is.nan(DR) & BS>0))  +			
			geom_point(aes(colour='data'), data=subset(seqd.str, !is.nan(DR) & BS==0))  +
			scale_colour_manual(name='Contig alignment', values=c('black','red')) +
			scale_y_continuous(breaks=seq(0,1,0.1)) +  
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey70", size=0.4), panel.grid.major=element_line(colour="grey70", size=0.4), legend.position='bottom') + labs(x='overlap between contigs\n(# sites)', y='prob evol distance accumulated within patient\n(Poisson model)')
	file	<- paste(indir, '/', substr(infile, 1, nchar(infile)-6),'_ProbPoisson_ContigOverlap.pdf',sep='')
	ggsave(file=file, w=6, h=6)
	
	
	#	fasta of contigs dual infection NO
	tmp		<- subset(seqd.str, BS==0 & POIS_GE_DN>0.9)
	tmp[, CONTIG_PAIR:=seq_len(nrow(tmp))]
	tmp		<- tmp[ , list(ID=c(ID1, ID2), LABEL= paste(c(ID1, ID2),'_pair',CONTIG_PAIR, '_p',round(POIS_GE_DN*100), sep='')), by='CONTIG_PAIR']
	z		<- seqa[tmp[,ID],]
	rownames(z)<- tmp[, LABEL]
	file	<- paste(indir, '/', substr(infile, 1, nchar(infile)-6),'_DualInfectionNo.fasta',sep='')
	write.dna(z, file=file, format='fasta')
	#	fasta of contigs dual infection YES
	tmp		<- subset(seqd.str, BS==0 & POIS_GE_DN<0.1)
	tmp[, CONTIG_PAIR:=seq_len(nrow(tmp))]
	tmp		<- tmp[ , list(ID=c(ID1, ID2), LABEL= paste(c(ID1, ID2),'_pair',CONTIG_PAIR, '_p',round(POIS_GE_DN*100), sep='')), by='CONTIG_PAIR']
	z		<- seqa[tmp[,ID],]
	rownames(z)<- tmp[, LABEL]
	file	<- paste(indir, '/', substr(infile, 1, nchar(infile)-6),'_DualInfectionYes.fasta',sep='')
	write.dna(z, file=file, format='fasta')
	seqd.str	
}

contigs.rawdist.bs<- function(seqn.str, seqa, bs.n)
{
	#	distance on actual alignment
	seqd		<- contigs.rawdist(seqn.str, seqa)
	seqd[, BS:=0]
	#	distance on bootstrapped alignments
	tmp			<- lapply(seq_len(bs.n), function(bs.i)
			{
				cat(paste('\ncomputing distance on bootstrap alignment',bs.i))
				tmp			<- floor( ncol(seqa )/3)
				tmp			<- sample(seq_len(tmp), tmp, replace=T)-1
				tmp			<- as.numeric( sapply(tmp,function(x)		3*x+c(1,2,3)		) )
				seqa.bs		<- seqa[, tmp]
				tmp			<- contigs.rawdist(seqn.str, seqa.bs)
				tmp[, BS:=bs.i]
				tmp
			})
	seqd		<- rbind(seqd, do.call('rbind', tmp))
	seqd
}	

contigs.rawdist<- function(seqn.str, seq.tmp)
{
	seqd.str<- subset(seqn.str, CONTIG_N>1)[, {
				tmp	<- seq.tmp[ID,]
				dn	<- dist.dna(tmp, model='N',as.matrix=1)			#number of differences, ignoring gaps in one seq and no gap in other seq 				
				di	<- dist.dna(tmp, model='indel',as.matrix=1)		#gaps in one seq and no gap in other seq
				dr	<- dist.dna(tmp, model='raw',as.matrix=1)		#N where none is gapped
				tmp	<- combn(length(ID),2)
				tmp	<- data.table(ID1= ID[tmp[1,]], ID2=ID[tmp[2,]], DN=NA_real_, DI=NA_real_, DR=NA_real_, ID1_N=NA_real_, ID2_N=NA_real_, ID12_N=NA_real_)
				set(tmp, NULL, 'DN', sapply(seq_len(nrow(tmp)), function(i)  dn[tmp[i,ID1], tmp[i,ID2]])   )
				set(tmp, NULL, 'DI', sapply(seq_len(nrow(tmp)), function(i)  di[tmp[i,ID1], tmp[i,ID2]])   )
				set(tmp, NULL, 'DR', sapply(seq_len(nrow(tmp)), function(i)  dr[tmp[i,ID1], tmp[i,ID2]])   )
				set(tmp, NULL, 'ID1_N', ncol(seq.tmp)-sapply( seq_len(nrow(tmp)), function(i)		base.freq( seq.tmp[ tmp[i, ID1], ], freq=1, all=1 )['-']		))
				set(tmp, NULL, 'ID2_N', ncol(seq.tmp)-sapply( seq_len(nrow(tmp)), function(i)		base.freq( seq.tmp[ tmp[i, ID2], ], freq=1, all=1 )['-']		))
				set(tmp, NULL, 'ID12_N', tmp[, DN/DR])
				tmp
			}, by='PATIENT']
	#if DN==0 and DR==0, then there is complete agreement betw seqs
	tmp		<- seqd.str[, which(DN==0 & DR==0)]
	set(seqd.str, tmp, 'ID12_N', seqd.str[tmp, abs(ID1_N-ID2_N)])	
	#if DN==0 and DR==NaN, then there is no overlap betw seqs
	set(seqd.str, seqd.str[, which(DN==0 & is.nan(DR))], 'ID12_N', 0)
	seqd.str
}

contigs.getblast<- function(seqn, indir, infile)
{
	#### parse BLAST	
	file	<- paste(indir, '/', infile, '.blast', sep='')
	seqb	<- seq.blast.read(file, sep = "\t")
	setnames(seqb, colnames(seqb), toupper(colnames(seqb)))
	# evalue means Expect value
	# pident means Percentage of identical matches
	# qstart means Start of alignment in query
	# qend means End of alignment in query
	# sstart means Start of alignment in subject
	# send means End of alignment in subject
	
	#	for each contig, retain best reference	
	tmp		<- seqb[, 	{
				tmp	<- which( EVALUE==min(EVALUE) )
				if(length(tmp)>1)
					tmp	<- tmp[ ALIGNMENT.LENGTH[tmp]==max(ALIGNMENT.LENGTH[tmp]) ]
				if(length(tmp)>1)
				{
					cat(paste('\nFound same EVALUE and ALIGNMENT.LENGTH for', QUERY.ID, 'ALIGNMENT.LENGTH=',ALIGNMENT.LENGTH[1]))
					tmp	<- tmp[1]
				}
				lapply(.SD, '[', tmp)
				#list(SUBJECT.ID= SUBJECT.ID[tmp])							
			} , by='QUERY.ID']
	stopifnot( length(setdiff( tmp[,QUERY.ID],  unique(subset(seqb, select=QUERY.ID))[, QUERY.ID] ))==0 )	#must keep all QUERY.IDs in seqb
	seqb	<- tmp
	# 	contigs with match in blast
	#seqn[, BLAST:= seqn[, ID%in%seqb[, QUERY.ID] ]]
	#cat(paste('\nContigs with no match in BLAST, n=', seqn[, length(which(BLAST==FALSE))]))
	#	to determine first cutoff, use SIV and everything not in group M O N
	setnames(seqb, 'QUERY.ID', 'ID')
	seqb	<- merge( seqb, subset(seqn, select=c(ID, CONTIG_LEN, COUNTRY)), by='ID' )
	seqb[, SUBTYPE:= seqb[, sapply( strsplit(SUBJECT.ID, '.', fixed=1), '[[', 2 ) ]]
	seqb[, SIV:= grepl('SIV', SUBJECT.ID)]
	#	criteria		
	ggplot(seqb, aes(x=ALIGNMENT.LENGTH, y=ALIGNMENT.LENGTH/CONTIG_LEN), log='x' ) +
			geom_point(alpha=0.5) + labs(x='alignment length', y='proportion of contig aligned') +
			geom_point(data=subset(seqb, SUBTYPE=='CPZ'), aes(colour='CPZ')) +
			geom_point(data=subset(seqb, SIV==TRUE), aes(colour='SIV')) + 
			coord_trans(xtrans="log10", limx=c(10, 12000)) + theme_bw() +
			scale_colour_brewer(name='not M N O group', palette='Set1') +
			scale_x_continuous(breaks=c(20, 50, 100, 200, 1000, 5000, 1000), minor_breaks=c(60, 70, 80, 90)) +
			scale_y_continuous(breaks=seq(0,1,0.1), minor_breaks=c(0.11,0.12, 0.13, 0.14, 0.15, 0.2, 0.21, 0.22, 0.23, 0.24)) +
			theme( panel.grid.minor=element_line(colour="grey70", size=0.4), panel.grid.major=element_line(colour="grey70", size=0.4) )
	file	<- paste(indir, '/INFO_BLAST_selection_criteria_SIVCPZ.pdf', sep='')
	ggsave(file=file, w=6, h=6)
	
	ggplot(seqb, aes(x=ALIGNMENT.LENGTH, y=ALIGNMENT.LENGTH/CONTIG_LEN), log='x' ) +
			geom_point(alpha=0.5) + labs(x='alignment length', y='proportion of contig aligned') +
			geom_point(data=subset(seqb, SUBTYPE%in%c('A1','A2','B','C','D','F1','F2','G','H','J','K','N','O','P')), aes(colour=SUBTYPE)) +
			coord_trans(xtrans="log10", limx=c(10, 12000)) + theme_bw() +
			scale_colour_discrete(name='subtype') +
			scale_x_continuous(breaks=c(20, 50, 100, 200, 1000, 5000, 1000), minor_breaks=c(60, 70, 80, 90)) +
			scale_y_continuous(breaks=seq(0,1,0.1), minor_breaks=c(0.11,0.12, 0.13, 0.14, 0.15, 0.2, 0.21, 0.22, 0.23, 0.24)) +
			theme( panel.grid.minor=element_line(colour="grey70", size=0.4), panel.grid.major=element_line(colour="grey70", size=0.4) )
	file	<- paste(indir, '/INFO_BLAST_selection_criteria_allsubtypes.pdf', sep='')
	ggsave(file=file, w=6, h=6)
	
	ggplot(seqb, aes(x=ALIGNMENT.LENGTH, y=ALIGNMENT.LENGTH/CONTIG_LEN), log='x' ) +
			geom_point(alpha=0.5) + labs(x='alignment length', y='proportion of contig aligned') +
			geom_point(data=subset(seqb, SUBTYPE%in%c('A1','A2','D','G')), aes(colour=SUBTYPE)) +
			coord_trans(xtrans="log10", limx=c(10, 12000)) + theme_bw() +
			scale_colour_discrete(name='subtype') +
			scale_x_continuous(breaks=c(20, 50, 100, 200, 1000, 5000, 1000), minor_breaks=c(60, 70, 80, 90)) +
			scale_y_continuous(breaks=seq(0,1,0.1), minor_breaks=c(0.11,0.12, 0.13, 0.14, 0.15, 0.2, 0.21, 0.22, 0.23, 0.24)) +
			theme( panel.grid.minor=element_line(colour="grey70", size=0.4), panel.grid.major=element_line(colour="grey70", size=0.4) )
	file	<- paste(indir, '/INFO_BLAST_selection_criteria_ADG.pdf', sep='')
	ggsave(file=file, w=6, h=6)
	
	#seqb[, table(SUBTYPE, SIV)]
	seqb[, BLAST_STRICT:= ALIGNMENT.LENGTH>350 & ALIGNMENT.LENGTH>0.9*CONTIG_LEN]
	seqb[, BLAST_RELAXED:= ALIGNMENT.LENGTH>350 & ALIGNMENT.LENGTH>0.1*CONTIG_LEN]
	print( seqb[, table(BLAST_STRICT, BLAST_RELAXED)] )
	#				BLAST_RELAXED
	#BLAST_STRICT 	FALSE  TRUE
	#FALSE 			10046   56
	#TRUE      		0   	296
	setnames(seqb, c('SUBTYPE', 'ALIGNMENT.LENGTH', 'IDENTITY', 'EVALUE'), c('BLAST_SUBTYPE', 'BLAST_ALIGNMENT_LENGTH', 'BLAST_IDENTITY', 'BLAST_EVALUE'))
	seqb
}

contigs.read.v141126<- function(indir)
{
	require(data.table)
		
	files	<- list.files(path=indir)
	#	get Extract_ID should be unique 
	files	<- data.table(FNAME=files)
	set(files, NULL, 'EXTRACT_ID', files[, paste( sapply(strsplit(FNAME,'-'), '[[', 1), '-', sapply(strsplit(FNAME,'-'), '[[', 2) , sep='') ] )
	#	stop if files not unique
	stopifnot( files[, length(unique(EXTRACT_ID))] == nrow(files) )
	#	get seqs 
	seq		<- lapply( seq_len(nrow(files)), function(i)
			{				
				tmp			<- paste(indir, files[i, FNAME], sep='/')
				cat(paste('\nreading', files[i, FNAME]))
				tmp			<- read.dna(tmp, format='fasta')
				cat(paste('\nread contigs, n=', length(tmp)))
				if(!is.list(tmp))
					tmp		<- as.list(tmp)
				names(tmp)	<- paste(files[i, EXTRACT_ID],names(tmp),sep='.')
				tmp
			})
	#	get DNAbin	
	seq		<- do.call('c', seq)
	#	deal with seq names
	seqn	<- data.table(SNAME=names(seq), INDEX=seq_along(seq))
	set(seqn, NULL, 'EXTRACT_ID', seqn[, sapply(strsplit(SNAME,'.',fixed=1), '[[', 1) ] )
	set(seqn, NULL, 'CONTIG', seqn[, sapply(strsplit(SNAME,'.',fixed=1), '[[', 4) ] )
	set(seqn, NULL, 'PATIENT', seqn[, sapply(strsplit(EXTRACT_ID,'-',fixed=1), '[[', 2) ] )
	set(seqn, NULL, 'COUNTRY', seqn[, substr(PATIENT,1,2) ] )	
	set(seqn, NULL, 'ID', seqn[, paste(PATIENT,CONTIG,COUNTRY,sep='|')])	
	#tmp		<- seqn[, list(ID= paste(COUNTRY, seq_along(COUNTRY), CONTIG, sep='|'), SNAME=SNAME ), by='COUNTRY']
	#seqn	<- merge(seqn, subset(tmp, select=c(SNAME, ID)), by='SNAME')
	setkey(seqn, INDEX)
	names(seq)	<- seqn[, ID]
	list(seq=seq, seq.names=seqn)
}

contigs.mafft.run.v141127<- function(seq, seqn, outdir)
{
	#	save selected contigs in single fasta file
	outfile	<- 'Contigs_141126_BLASTstrict_141127'
	file	<- paste(outdir, '/', outfile, '.fasta', sep='')
	write.dna(seq[ subset(seqn, BLAST_STRICT)[, ID] ], format='fasta', file=file )
	#	BLAST
	cmd		<- hivc.cmd.blast(outdir, outfile, '', '/Users/Oliver/duke/2014_Gates1125_dualinfection/R_code', 'HIV1_REF_2010_genome_DNA_nogap', dbsignat='', blast.task="blastn", blast.max_target_seqs=6, blast.evalue=10, blast.wordsize=11, blast.gapopen=5, blast.gapextend=2, blast.penalty=-3, blast.reward= 2, blast.dust= "no", nproc=1, verbose=1)
	cat ( cmd )	
	
	file	<- paste(outdir, '/', outfile, '.blast', sep='')
	seqb	<- seq.blast.read(file, sep = "\t")
	setnames(seqb, colnames(seqb), toupper(colnames(seqb)))
	setnames(seqb, 'QUERY.ID', 'ID')
	seqb	<- subset(seqb, ALIGNMENT.LENGTH>350)
	#	number of ref seqs
	print( seqb[, list(REF_N=length(SUBJECT.ID)), by='ID'][, table(REF_N)] )
	#	load refseq	
	file	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/R_code/HIV1_REF_2010_genome_DNA.fasta'
	ref		<- read.dna(file, format='fasta')	
	#	write fasta files
	df		<- seqb[,{
				tmp	<- as.list( ref[SUBJECT.ID,] )
				tmp	<- c(seq[ID],tmp)
				z	<- paste( outdir, '/', outfile, '_', ID, '.fasta', sep='' )
				cat(paste('\nwrite to file', z))
				write.dna( tmp, file=z, format='fasta', append=0, colsep='', colw=80, blocksep=0)
				list(FNAME= z)
			}, by='ID']
	#	accurate mafft for every fasta file
	#	linsi --thread -1 /Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis/Contigs_141126_BLASTstrict.fasta > /Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis/Contigs_141126_BLASTstrict_MAFFTlinsi.fasta
	tmp		<- df[,	{	
				cat(paste('\nprocess',ID))
				cmd	<- paste( 'linsi --quiet --thread -4 \"', FNAME, '\" > \"', substr(FNAME, 1, nchar(FNAME)-6), '_linsi.fasta\"', sep='' )
				system(cmd)
			},by='ID']
		
}

contigs.mafft.get.v141127<- function(indir, infile)
{
	require(data.table)
	infile	<- 'Contigs_141126_BLASTstrict_OR'
	#files	<- data.table(FROM= list.files(path=indir, 'ATHENA'))
	#files[, TO:= files[, gsub('.*Sequences','Contigs_141126_BLASTstrict_OR',FROM)]]
	#files[, file.rename(paste(indir,'/',FROM,sep=''), paste(indir,'/',TO,sep='')), by='FROM']
	files	<- data.table(FNAME=list.files(path=indir, infile))
	files[, ID:=files[, sapply(strsplit(FNAME,'_'), '[[', 5) ]]
	set(files, NULL, 'ID', files[, gsub('.fasta','',ID)])
	files[, TYPE:= 'FNAME_IN']
	set(files, files[, which(grepl('linsi', FNAME))], 'TYPE', 'FNAME_OUT')
	set(files, files[, which(grepl('realigned', FNAME))], 'TYPE', 'REALIGNED')
	
	files	<- dcast.data.table(files, ID~TYPE, value.var='FNAME')
	#	for not realigned OUT files:
	#	read fasta IN and OUT and re-align OUT so that it aligns with IN --> we **ASSUME** that no gaps were introduced into the reference alignment, otherwise this may fail 
	files	<- subset(files, is.na(REALIGNED))	
	files[, {
				seqi	<- read.dna( paste(indir, '/', FNAME_IN, sep=''), format='fasta' )
				seqa	<- as.character( read.dna( paste(indir, '/', FNAME_OUT, sep=''), format='fasta' ) )
				print(dim(seqa))
				stopifnot( setequal(names(seqi), rownames(seqa)) )
				
				ref		<- names(seqi)[ grepl('^Ref',names(seqi)) ][1]
				cat(paste('\nrefseq is', ref))
				seqri	<- as.character( seqi[ref] )[[1]]
				seqri	<- c(seqri, rep('-', max( length(seqri), ncol(seqa) )-length(seqri)))
				seqa	<- cbind(seqa, matrix('-', nrow=nrow(seqa), ncol=max( length(seqri), ncol(seqa) )-ncol(seqa)))
				print(dim(seqa))
				seqra	<- seqa[ref, ]
				tmp		<- which( seqri=='-' & seqra!='-' )		#gap in ref that is not in alignment -> expand
				tmpa	<- which( seqri!='-' & seqra=='-' )		#gap in alignment that is not in ref -> swap
				#print(tail(seqri, 200))
				#print(tail(seqra,1500))				
				print( rbind(seqri,seqra)[,1:50] )
				print( rbind(seqri,seqra)[,1:50] )
				#stop()
				#for loop
				while(length(tmp) | length(tmpa))
				{
					if((length(tmp) & length(tmpa) & tmp[1]<tmpa[1]) | (length(tmp) & !length(tmpa)))
					{
						cat(paste('\nalignment: adding gap at pos', tmp[1]))
						if(tmp[1]==1)
							seqa	<- cbind( matrix('-',nrow=nrow(seqa), ncol=1), seqa[, -ncol(seqa)] )
						if(tmp[1]>1)
							seqa	<- cbind( seqa[, 1:(tmp[1]-1), drop=0], matrix('-',nrow=nrow(seqa), ncol=1), seqa[, tmp[1]:(ncol(seqa)-1)] )
					}						
					if((length(tmpa) & length(tmp) & tmpa[1]<tmp[1]) | (length(tmpa) & !length(tmp)))
					{
						z			<- which( seqra[tmpa[1]:length(seqra)]!='-' )	#swap with first non-gap
						if(length(z)==0)
						{
							print(tmpa[1])
							print(length(seqra))
							#print( rbind(seqri,seqra) )							
						}
						#print( rbind(seqri,seqra)[,9850:9870] )						
						#print(tail(seqri, 200))
						#print(tail(seqra,1500))						
						stopifnot(length(z)>0)
						z			<- tmpa[1]	+ z[1] -1
						cat(paste('\nalignment: swap pos', tmpa[1], 'with pos',z))						
						if(z==(tmpa[1]+1))
							seqa	<- seqa[, c(1:(tmpa[1]-1), tmpa[1]+1, tmpa[1], (tmpa[1]+2):ncol(seqa))]
						if(z>(tmpa[1]+1))
							seqa	<- seqa[, c(1:(tmpa[1]-1), z, (tmpa[1]+1):(z-1), tmpa[1], (z+1):ncol(seqa))]
						#if((tmpa[1]+1)==ncol(seqa))
						#	seqa	<- seqa[, c(1:(tmpa[1]-1), tmpa[1]+1, tmpa[1])]
					}
					#print(tail(seqra[seqra!='-'],5))
					#if(!all( tail(seqra[seqra!='-'],5)==c("a", "g", "t", "g", "c")	))
					#		print(seqra)
					#stopifnot( all( tail(seqra[seqra!='-'],5)==c("a", "g", "t", "g", "c")	) )
					#print(dim(seqa))
					seqra	<- seqa[ref, ] 
					#seqra	<- c(seqra, rep('-', max( length(seqri), length(seqra) )-length(seqra)))
					tmp		<- which( seqri=='-' & seqra!='-' )
					tmpa	<- which( seqri!='-' & seqra=='-' )
					#print( tmp )
					#print( tmpa )
					
					#print( rbind(seqri,seqra)[,410:420] )
					#print( rbind(seqri,seqra)[,830:840] )
					#if(tmp[1]==10576 & tmp[1]<tmpa[1])
					#{
					#	print( rbind(seqri,seqra) )
					#	stop()
					#} 
				}											
				tmp	<-  which(seqri!=seqra)
				if(length(tmp))
					print(FNAME_OUT)
				stopifnot(length(tmp)==0)
				tmp	<- paste(indir,'/',gsub('linsi','realigned',FNAME_OUT),sep='')
				cat(paste('\nWrite to file',tmp))
				write.dna(as.DNAbin(seqa), file=tmp, format='fasta')
				
			}, by='ID']
	# read realigned files
	files	<- data.table(FNAME=list.files(path=indir, infile))
	files[, ID:=files[, sapply(strsplit(FNAME,'_'), '[[', 5) ]]
	set(files, NULL, 'ID', files[, gsub('.fasta','',ID)])
	files[, TYPE:= 'FNAME_IN']
	set(files, files[, which(grepl('linsi', FNAME))], 'TYPE', 'FNAME_OUT')
	set(files, files[, which(grepl('realigned', FNAME))], 'TYPE', 'REALIGNED')	
	files	<- dcast.data.table(files, ID~TYPE, value.var='FNAME')

	seqa	<- lapply( files[, REALIGNED], function(x)	read.dna(file=paste(indir,'/',x,sep=''),format='fasta') )
	stopifnot( all( ncol(seqa[[1]]) == sapply(seqa, ncol) ) ) 
	seqa	<- do.call('rbind',seqa)
	seqa	<- seqa[!grepl('^Ref',rownames(seqa)),]
	
	file	<- paste(indir, '/', infile, '_linsi.fasta',sep='')
	write.dna(seqa, file=file, format='fasta', append=0, colsep='', colw=80, blocksep=0)
	seqa
}


contigs.mafft.run.vFrancois<- function(seq, seqn, outdir)
{
	#	save seqs that satisfy BLAST criteria along with refseq	
	file	<- '/Users/Oliver/duke/2014_Gates1125_dualinfection/R_code/HIV1_REF_2010_genome_DNA.fasta'
	ref		<- read.dna(file, format='fasta')
	ref		<- seq.rmgaps(ref, rm.only.col.gaps=0)
	tmp		<- seq[ subset(seqn, BLAST_STRICT)[, ID] ]
	tmp		<- c(tmp, ref)
	file	<- paste(outdir, 'Contigs_141126_BLASTstrict.fasta', sep='/')
	write.dna(tmp, file=file, format="fasta", append=0, colsep='', colw=80, blocksep=0)
	tmp		<- seq[ subset(seqn, BLAST_RELAXED)[, ID] ]
	tmp		<- c(tmp, ref)
	file	<- paste(outdir, 'Contigs_141126_BLASTrelaxed.fasta', sep='/')
	write.dna(tmp, file=file, format="fasta", append=0, colsep='', colw=80, blocksep=0)
	
	#### ... here align the contigs using MAFFT (done outside of R)... ####
	# We perform an alignment of HIV contigs + 170 reference genomes
	# using the best MAFFT settings (i.e opening penalty 2, extension penalty 0.5)
	# cf file alignment_quality.csv
	#	fast mafft for now:
	#	mafft --op 2 --ep 0.5 --thread -4 /Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis/Contigs_141126_BLASTstrict.fasta > /Users/Oliver/duke/2014_Gates1125_dualinfection/Contigs_141126_analysis/Contigs_141126_BLASTstrict_Francois_mafft.fasta		
}

contigs.mafft.get.vFrancois<- function(indir, infile)
{
	file	<- paste(indir, '/', infile, sep='')
	seqa	<- read.dna(file, format='fasta')
	#	cut for now at 941 and 13250
	seqa	<- as.character(seqa)
	#seqa[460:466, 935:945]	seqa[260:266, 13245:13255]
	seqa	<- as.DNAbin(seqa[,941:13255])		
}
