project.dual<- function()
{
	HOME		<<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA'
	project.dual.distances.231015()
}

project.dual.distances.231015<- function()
{
	indir		<- paste(HOME,"alignments_151023",sep='/')
	indir		<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_151023'
	infiles		<- list.files(indir, pattern='R$')
	
	for(i in seq_along(infiles))
	{
		load( paste(indir,'/',infiles[i],sep='') )
		d		<- seq.dist(seq)
		save(d, seq, file= gsub('\\.R','_dist\\.R',paste(indir,'/',infiles[i],sep='')))
		gc()
	}		
}

project.dual.alignments.231015<- function()
{
	outdir<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_151023'
	#	read info
	file	<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	si		<- as.data.table(read.csv(file, stringsAsFactors=FALSE))
	setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	setkey(si, PANGEA_ID)
	
	#	read global PANGEA alignment and split by site
	file	<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_GlobalAlignment.fasta"
	sq		<- read.dna(file, format='fasta')
	sqi		<- data.table(TAXA=rownames(sq))
	sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]		
	sqi[, SITE:= NA_character_]
	tmp		<- sqi[, which(PNG=='Y')]
	set(sqi, tmp, 'SITE', sqi[tmp, substring(sapply(strsplit(TAXA,'-'),'[[',2),1,2)])
	setnames(sqi, 'TAXA', 'PANGEA_ID')
	sqi		<- merge(sqi, unique(si), by='PANGEA_ID', all.x=1)
	sqi		<- subset(sqi, is.na(CLINICAL_GENOME_COVERAGE) | CLINICAL_GENOME_COVERAGE>0)
	seq		<- sq[ subset(sqi, SITE=='UG' | PNG=='N')[, PANGEA_ID], ]
	write.dna( seq, file=paste(outdir,'/PANGEAconsensuses_2015-09_Imperial_UG.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEAconsensuses_2015-09_Imperial_UG.R',sep=''))
	seq		<- sq[ subset(sqi, SITE=='BW' | PNG=='N')[, PANGEA_ID], ]
	write.dna( seq, file=paste(outdir,'/PANGEAconsensuses_2015-09_Imperial_BW.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEAconsensuses_2015-09_Imperial_BW.R',sep=''))
	
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
	cri		<- merge(cri, si, by='SANGER_ID', all.x=1)
	stopifnot( nrow(subset(cri, is.na(SANGER_ID) & PNG=='Y'))==0 )
	cri[, SITE:= NA_character_]
	tmp		<- cri[, which(PNG=='Y')]
	set(cri, tmp, 'SITE', cri[tmp, substring(sapply(strsplit(PANGEA_ID,'-'),'[[',2),1,2)])		
	tmp		<- cri[, list(TAXA_NEW= ifelse( is.na(CNTG_ID_NEW), paste('Ref.',TAXA,sep=''), paste(PANGEA_ID,'-C',CNTG_ID_NEW,sep='') )), by='TAXA']
	cri		<- merge(cri, tmp, by='TAXA')
	setkey(cri, TAXA)
	rownames(cr)	<- cri[rownames(cr),][, TAXA_NEW]
	seq		<- cr[ subset(cri, SITE=='UG' | PNG=='N')[, TAXA_NEW], ]
	write.dna( seq, file=paste(outdir,'/PANGEAcontigs_2015-09_Imperial_UG.fasta',sep=''), format='fasta', colsep='', nbcol=-1)
	save( seq, file=paste(outdir,'/PANGEAcontigs_2015-09_Imperial_UG.R',sep=''))
	seq		<- cr[ subset(cri, SITE=='BW' | PNG=='N')[, TAXA_NEW], ]
	write.dna( seq, file=paste(outdir,'/PANGEAcontigs_2015-09_Imperial_BW.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEAcontigs_2015-09_Imperial_BW.R',sep=''))
	
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
