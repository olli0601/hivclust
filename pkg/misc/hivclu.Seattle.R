seattle.wrapper<- function()
{
	seattle.170621.fastree()
}

seattle.190723.subgraph.empirics<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	require(reshape)
	
	#	working directory with phyloscanner output		
	indir.phsc	<- '~/Box Sync/OR_Work/Seattle/phyloscanner_out_190723/KCHSX'
	outdir <- '~/Box Sync/OR_Work/Seattle/analysis_190723/KCHSX'
	
	#
	#	extract subgraph taxa	
	#
	infiles		<- data.table(F=list.files(indir.phsc, pattern='subgraphs_KCHSX.rda$', full.names=TRUE, recursive=TRUE))
	infiles[, SELECT:= gsub('^(.*)_subgraphs_([A-Z]+)\\.rda','\\2',basename(F))]	
	infiles[, ST:= gsub('^Subtype_(.*)_([0-9]+)_subgraphs_([A-Z]+)\\.rda','\\1',basename(F))]
	infiles[, REP:= gsub('^Subtype_(.*)_([0-9]+)_subgraphs_([A-Z]+)\\.rda','\\2',basename(F))]
	dsubgraphtaxa <- infiles[, {
				#i<- 1
				#infile <- infiles[i, F]			
				infile <- F
				cat('Process',infile,'\n')
				load(infile)			
				subgraph.names <- unlist(sapply( subgraphs, function(subgraph)  rep(subgraph$subgraph.name, length(subgraph$tip.label)) ))
				subgraph.taxa <- unlist(sapply( subgraphs, function(subgraph)  subgraph$tip.label))				
				list(	NAME=subgraph.names, 
						TAXA= subgraph.taxa 
				)				
			}, by=c('ST','REP','SELECT')]	
	#	add meta data from taxa names
	dsubgraphtaxa[, ID:= as.numeric(gsub('^([A-Z]+)_+PR/RT-([0-9]+)-([0-9]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z_]+)-([0-9]+)$','\\3',TAXA))]
	dsubgraphtaxa[, LOC:= gsub('^([A-Z]+)_+PR/RT-([0-9]+-[0-9]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z_]+)-([0-9]+)$','\\3',TAXA)]
	dsubgraphtaxa[, ETH:= gsub('^([A-Z]+)_+PR/RT-([0-9]+-[0-9]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z_]+)-([0-9]+)$','\\4',TAXA)]
	dsubgraphtaxa[, BORN:= gsub('^([A-Z]+)_+PR/RT-([0-9]+-[0-9]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z_]+)-([0-9]+)$','\\5',TAXA)]
	dsubgraphtaxa[, SEX:= gsub('^([A-Z]+)_+PR/RT-([0-9]+-[0-9]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z_]+)-([0-9]+)$','\\6',TAXA)]
	dsubgraphtaxa[, TRM:= gsub('_[M|F]','',gsub('^([A-Z]+)_+PR/RT-([0-9]+-[0-9]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z_]+)-([0-9]+)$','\\7',TAXA))]
	dsubgraphtaxa[, POSDATE:= as.numeric(gsub('^([A-Z]+)_+PR/RT-([0-9]+-[0-9]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z]+)-([a-zA-Z_]+)-([0-9]+)$','\\8',TAXA))]
	dsubgraphtaxa[, FULL_NAME:= paste0(ST,'_',REP,'_',SELECT,'_',NAME)]
	#	add meta data from persons file
	infile.meta <- '~/Box Sync/OR_Work/Seattle/PHSKC-2018-07-09/person.rds'
	dmeta <- as.data.table(readRDS(infile.meta))
	setnames(dmeta, c('newnum','b_yr','race','sex','gender'), c('ID','BIRTH_YEAR2','RACE2','SEX2','GENDER2'))
	tmp <- subset(dmeta, select=c(ID,BIRTH_YEAR2,RACE2,SEX2,GENDER2))
	set(tmp, NULL, 'BIRTH_YEAR2', tmp[, as.numeric(BIRTH_YEAR2)])
	dsubgraphtaxa <- merge(dsubgraphtaxa,tmp,by='ID')
	dsubgraphtaxa[, AGE2020:= cut(2020-BIRTH_YEAR2, breaks=c(-Inf,20,25,30,35,40,45,50,55,60,Inf))]
	
	# ethnicity
	dsubgraphtaxa[, DUMMY:= ETH]
	dcnt <- dsubgraphtaxa[, list(N= length(ID)), by=c('FULL_NAME','DUMMY')]
	ggplot(dcnt, aes(x=FULL_NAME, y=N, fill=DUMMY)) +
			geom_bar(stat='identity') + 
			coord_flip() +
			labs(y='#individuals',x='phylo subgraphs',fill='ethnicity')
	ggsave(file=file.path(outdir, 'empirical_subgraph_composition_ethnicity.pdf'), h=15,w=6)
	# birth place
	dsubgraphtaxa[, DUMMY:= BORN]
	dcnt <- dsubgraphtaxa[, list(N= length(ID)), by=c('FULL_NAME','DUMMY')]
	ggplot(dcnt, aes(x=FULL_NAME, y=N, fill=DUMMY)) +
			geom_bar(stat='identity') + 
			coord_flip() +
			labs(y='#individuals',x='phylo subgraphs',fill='birth place')
	ggsave(file=file.path(outdir, 'empirical_subgraph_composition_birthplace.pdf'), h=15,w=6)
	# sex
	dsubgraphtaxa[, DUMMY:= SEX]
	dcnt <- dsubgraphtaxa[, list(N= length(ID)), by=c('FULL_NAME','DUMMY')]
	ggplot(dcnt, aes(x=FULL_NAME, y=N, fill=DUMMY)) +
			geom_bar(stat='identity') + 
			coord_flip() +
			labs(y='#individuals',x='phylo subgraphs',fill='sex')
	ggsave(file=file.path(outdir, 'empirical_subgraph_composition_sex.pdf'), h=15,w=6)
	# age2020
	dsubgraphtaxa[, DUMMY:= AGE2020]
	dcnt <- dsubgraphtaxa[, list(N= length(ID)), by=c('FULL_NAME','DUMMY')]
	ggplot(dcnt, aes(x=FULL_NAME, y=N, fill=DUMMY)) +
			geom_bar(stat='identity') + 
			coord_flip() +
			labs(y='#individuals',x='phylo subgraphs',fill='age in 2020')
	ggsave(file=file.path(outdir, 'empirical_subgraph_composition_age2020.pdf'), h=15,w=6)
	
}

seattle.190723.test.simmap.casting <- function()
{
	require(data.table)
	require(tidyverse)
	require(ape)
	require(phytools)
	require(phangorn)
	require(phyloscannerR)
	
	infile <- "/Users/Oliver/Box Sync/OR_Work/Seattle/analysis_190723/phyloscanner_HSX/Subtype_01_AE_0_workspace.rda"
	load(infile)
	ptree <- phyloscanner.trees[[1]]
	tree <- ptree[['tree']]
	ph <- phyloscanner.to.simmap(tree)	
	attr(ph, "SPLIT") <- NULL
	attr(ph, "INDIVIDUAL") <- NULL
	attr(ph, "BRANCH_COLOURS") <- NULL
	attr(ph, "SUBGRAPH_MRCA") <- NULL
	ph <- simmap.to.phyloscanner(ph)
	ph[["maps"]] <- ph[["mapped.edge"]] <- ph[["node.states"]] <- NULL
	attr(ph, 'map.order') <- NULL
	attr(ph, 'class') <- 'phylo'
	
	stopifnot(all( tree$tip.label == ph$tip.label ))
	stopifnot(all( tree$edge == ph$edge ))
	stopifnot(all( tree$edge.length == ph$edge.length ))
	stopifnot(all( tree$Nnode == ph$Nnode ))
	
	tmp <- as.character(attr(tree,'SPLIT'))
	tmp2 <- as.character(attr(ph,'SPLIT'))
	tmp[is.na(tmp)] <- 'Unknown'
	tmp2[is.na(tmp2)] <- 'Unknown'
	stopifnot( all(tmp==tmp2) )
	
	tmp <- as.character(attr(tree,'INDIVIDUAL'))
	tmp2 <- as.character(attr(ph,'INDIVIDUAL'))
	tmp[is.na(tmp)] <- 'Unknown'
	tmp2[is.na(tmp2)] <- 'Unknown'
	stopifnot( all(tmp==tmp2) )
	
	tmp <- as.character(attr(tree,'BRANCH_COLOURS'))
	tmp2 <- as.character(attr(ph,'BRANCH_COLOURS'))
	tmp[is.na(tmp)] <- 'Unknown'
	tmp2[is.na(tmp2)] <- 'Unknown'
	stopifnot( all(tmp==tmp2) )
		
	tmp <- attr(tree,'SUBGRAPH_MRCA')
	tmp2 <- attr(ph,'SUBGRAPH_MRCA')
	stopifnot( all(tmp==tmp2) )
	
}

seattle.190723.treedater <- function()
{
	require(data.table)
	require(tidyverse)
	require(ape)
	require(phytools)
	require(treedater)
	
	#infile.seqinfo <- '~/Box Sync/OR_Work/Seattle/PHSKC-2019_07_24/sequences_meta.rds'
	#new file above not compatible with old data: cannot find sequence ids.
	infile.seqinfo <- '~/Box Sync/OR_Work/Seattle/PHSKC-2018-07-09/sequences_meta.rds'	
	indir.trees	<- '~/Box Sync/OR_Work/Seattle/analysis_190723/trees'	
	infiles.trees <- data.table(F=list.files(indir.trees, pattern='rerooted\\.newick$', recursive=TRUE, full.names=TRUE))
	indir.phsc	<- '~/Box Sync/OR_Work/Seattle/analysis_190723/phyloscanner_HSX'
	infiles.phsc <- data.table(F=list.files(indir.phsc, pattern='workspace.rda$', full.names=TRUE, recursive=TRUE))	
	alignment.length <- 1100
	
	#
	# read Seattle sampling data 
	#
	
	dseq <- readRDS(infile.seqinfo) 
	dseq <- dseq %>% 
			select(seqID, newnum, type, seqy, seqm) %>%
			mutate( SEQ_DATE:= seqy + (seqm-1)/12 + 15/365,
					SEQ_DATE_LOWER:= seqy + (seqm-1)/12,
					SEQ_DATE_UPPER:= seqy + (seqm-1)/12 + 30/365
					) %>%
			select(-seqy, -seqm) %>%
			rename(TAXA= seqID, PID= newnum, TYPE=type)
	dseq <- as.data.table(dseq)
	
	#
	#	for each tree: 
	#	make data.table of sequence sampling times 
	#	remove taxa without data on sampling times
	#
	for(i in seq_len(nrow(infiles.phsc)))
	{
		infile <- infiles.phsc[i,F]	
		load(infile)
		ph <- phyloscanner.trees[[1]][['tree']]
		stopifnot( !any( ph$tip.label=='' ) )
		
		#
		#	drop tips without sequence dates, 
		#	while conserving the ancestral state reconstructions
		#
		
		#	extract taxa names
		dph <- data.table(	TAXA_LABEL=ph$tip.label,
				TAXA= gsub('^(PR/RT-[0-9]+).*','\\1',gsub('^[^_]+___(.*)','\\1',ph$tip.label)),
				TAXA_ID= seq_along(ph$tip.label))
		#	add Seattle sequence dates
		dph <- merge(dph, dseq, by='TAXA', all.x=TRUE)
		#	extract GenBank sequence dates from taxa names where possible 
		dph[, GENBANK_ID:= gsub('^([^\\.]+)\\.([^\\.]+)\\.([^\\.]+)\\.([^\\.]+)\\.([^\\.]+)\\.([^\\.]+)$','\\5',TAXA)]
		dph[, GENBANK_SEQDATE:= gsub('^([^\\.]+)\\.([^\\.]+)\\.([^\\.]+)\\.([^\\.]+)\\.([^\\.]+)\\.([^\\.]+)$','\\3',TAXA)]
		set(dph, dph[, which(GENBANK_SEQDATE=='-' | !is.na(SEQ_DATE))],'GENBANK_SEQDATE',NA_character_)	
		set(dph, NULL, 'GENBANK_SEQDATE_LOWER', dph[, as.numeric(GENBANK_SEQDATE) ])
		set(dph, NULL, 'GENBANK_SEQDATE_UPPER', dph[, as.numeric(GENBANK_SEQDATE) + 364/365 ])	
		set(dph, NULL, 'GENBANK_SEQDATE', dph[, as.numeric(GENBANK_SEQDATE) + 1/2 ])
		tmp <- dph[, which(is.na(SEQ_DATE))]
		set(dph, tmp, 'SEQ_DATE', dph[tmp, GENBANK_SEQDATE])
		set(dph, tmp, 'SEQ_DATE_LOWER', dph[tmp, GENBANK_SEQDATE_LOWER])
		set(dph, tmp, 'SEQ_DATE_UPPER', dph[tmp, GENBANK_SEQDATE_UPPER])
		set(dph, NULL, c('GENBANK_SEQDATE','GENBANK_SEQDATE_LOWER','GENBANK_SEQDATE_UPPER'), NULL)	
		#	drop tips 
		dph.old <- subset(dph, select=c(TAXA, SEQ_DATE, SEQ_DATE_LOWER, SEQ_DATE_UPPER))	
		tmp <- subset(dph, is.na(SEQ_DATE))[, TAXA_ID]
		cat('Dropping tips without sampling date from ', infile,' n=', length(tmp), 'of Ntips=', Ntip(ph), '\n')
		ph <- phyloscanner.to.simmap(ph)	
		ph <- phytools:::drop.tip.simmap(ph, ph$tip.label[tmp])
		ph <- simmap.to.phyloscanner(ph)
		
		
		#
		#	date tree
		#	
		
		#	make data.table of sequence sampling times
		dph <- data.table(	TAXA_LABEL=ph$tip.label,
				TAXA= gsub('^(PR/RT-[0-9]+).*','\\1',gsub('^[^_]+___(.*)','\\1',ph$tip.label)),
				TAXA_ID= seq_along(ph$tip.label))
		dph <- merge(dph, dph.old, by= 'TAXA')
		stopifnot( !any(is.na(dph$SEQ_DATE)) )
		dph <- dph[order(TAXA_ID),]
		#	get into format needed for tree.dater	
		#sampling.times <- dph$SEQ_DATE
		#names(sampling.times) <- dph$TAXA_LABEL	
		sampling.times.init <- dph$SEQ_DATE
		names(sampling.times.init) <- dph$TAXA_LABEL
		sampling.times.bounds <- as.data.frame(subset(dph, select=c(SEQ_DATE_LOWER, SEQ_DATE_UPPER)))
		rownames(sampling.times.bounds) <- dph$TAXA_LABEL		
		colnames(sampling.times.bounds) <- c('lower','upper')
		#	date tree
		#ph.dated <- dater(ph, sampling.times, alignment.length, numStartConditions=1)
		ph.dated <- dater(ph, sampling.times.init, alignment.length, numStartConditions=1, estimateSampleTimes=sampling.times.bounds )
		stopifnot( all(  ph.dated$tip.label == ph$tip.label ) )
		stopifnot( all( ph.dated$edge == ph$edge ) )
		#	since the tree topology is unchanged, we can copy
		#	the branch lenghts in units of time onto the original tree
		#	that has the ancestral state reconstructions
		ph$edge.length <- ph.dated$edge.length
		
		#
		#	plot dated tree to spot obvious errors
		#
		outfile <- gsub('workspace\\.rda','annotated_dated_tree.pdf',infile)
		tmp <- vector('list')
		tmp[['tree']] <- ph
		tmp[['tree']][['node.states']] <- tmp[['tree']][['mapped.edge']] <- tmp[['tree']][['maps']] <- NULL
		attr(tmp[['tree']],'map.order') <- NULL
		attr(tmp[['tree']],'class') <- 'phylo'
		tmp[['read.counts']] <- rep(1, Ntip(ph))	
		write.annotated.tree(tmp, outfile, format="pdf", pdf.scale.bar.width = 0.01, pdf.w = 40, pdf.hm = 0.2, verbose = FALSE)
		
		#
		#	save phyloscanner.tree
		#
		outfile <- gsub('workspace\\.rda','annotated_dated_tree.rda',infile)
		save(ph, file=outfile)
	}	
}

seattle.190723.get.newick.from.phsc<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	require(reshape)
	require(phyloscannerR)
	
	#	working directory with phyloscanner output		
	indir.phsc	<- '~/Box Sync/OR_Work/Seattle/analysis_190723/phyloscanner_HSX'
	infiles		<- data.table(F=list.files(indir.phsc, pattern='workspace.rda$', full.names=TRUE, recursive=TRUE))
	for(i in seq_len(nrow(infiles)))
	{
		#i<- 1
		cat('process', i,'\n')
		infile <- infiles[i, F]	
		outfile <- gsub('_workspace\\.rda$', '.newick', infile)
		load(infile)
		ph <- phyloscanner.trees[[1]][['tree']]
		write.tree(ph, file=outfile)
	}
}

seattle.190723.get.subgraphs<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	require(reshape)
	require(phyloscannerR)
	#	working directory with phyloscanner output		
	indir.phsc	<- '~/Box Sync/OR_Work/Seattle/analysis_190723/phyloscanner_HSX'
	infiles		<- data.table(F=list.files(indir.phsc, pattern='workspace.rda$', full.names=TRUE, recursive=TRUE))
	infiles[, SELECT:= 'KCHSX']	
	
	#
	#	make pdf file
	#
	for(i in seq_len(nrow(infiles)))
	{
		#i<- 1
		cat('process', i,'\n')
		infile <- infiles[i, F]
		load(infile)
		outfile <- gsub('_workspace\\.rda','_annotated_tree.pdf',infile)
		ptree <- phyloscanner.trees[[1]]
		write.annotated.tree(ptree, outfile, format="pdf", pdf.scale.bar.width = 0.01, pdf.w = 40, pdf.hm = 0.2, verbose = FALSE)
	}
	
	#
	#	extract subgraphs of undated trees	
	#
	for(i in seq_len(nrow(infiles)))
	{
		#i<- 1
		cat('process', i,'\n')
		infile <- infiles[i, F]
		host <- infiles[i,SELECT]
		load(infile)	
		ph <- phyloscanner.trees[[1]][['tree']]		
		mrcas <- which( attr(ph, 'SUBGRAPH_MRCA') )
		mrcas <- mrcas[ attr(ph, 'INDIVIDUAL')[mrcas]==host ]	
		# convert tree to class simmap
		ph <- phyloscanner.to.simmap(ph)
		# extract subgraphs
		subgraphs <- lapply(mrcas, function(mrca) extract.subgraph(ph, mrca))
		# save
		outfile <- gsub('_workspace',paste0('_subgraphs_',host),infile)
		save(subgraphs, file=outfile)
	}
	
	#
	#	extract subgraphs of dated trees
	#
	infiles		<- data.table(F=list.files(indir.phsc, pattern='_annotated_dated_tree.rda$', full.names=TRUE, recursive=TRUE))
	infiles[, SELECT:= 'KCHSX']	
	for(i in seq_len(nrow(infiles)))
	{
		#i<- 1
		cat('process', i,'\n')
		infile <- infiles[i, F]
		host <- infiles[i,SELECT]
		load(infile)					
		mrcas <- which( attr(ph, 'SUBGRAPH_MRCA') )
		mrcas <- mrcas[ attr(ph, 'INDIVIDUAL')[mrcas]==host ]	
		# check that tree is of class simmap
		stopifnot( any(attr(ph,'class')=='simmap') )
		# extract subgraphs
		subgraphs <- lapply(mrcas, function(mrca) extract.subgraph(ph, mrca))
		# save
		outfile <- gsub('_annotated_dated_tree',paste0('_datedsubgraphs_',host),infile)
		save(subgraphs, file=outfile)
	}
	
}

seattle.levu.180613<- function()
{
	#	https://github.com/slevu/tenbrit
	require(tenbrit)
	
	# 	clean data and save to file
	#	this has to be modified
	#	TODO get UK data just for replication purposes
	import_rawdata()
	#	-- remove all UK seqs from LANL
	
	
	#	https://github.com/slevu/garel
	require(garel)
	
	#	- make BLAST database by subtype	
	#	-- split LANL alignment into group of 200
	#	-- for each group, remove DRMs	
	pipeline_blastdb
		
	#	- find closest LANL sequences for target alignment
	#	-- remove DRMs from target alignment	
	#	-- BLAST against LANL seqs, save BLAST output
	#	-- keep best hit for each target sequence
	#	-- MAFFT to align target alignment with best hit LANL sequences (not codon aligned)
	#	-- refer to this alignment as input alignment
	pipeline_msa
	
	#	- create tree on input alignment with RaXML
	#	-- to create UNIX scripts, use funr package https://cran.r-project.org/web/packages/funr/index.html
	#	-- date with LSD 
	pipeline_trees
	
	#	- creates bootstrap trees 
	#	-- create bootstrap alignment from input alignment
	#	-- run RAxML to create bootstrap trees
	#	-- run LSD
	pipeline_bstrees
	
	#	since subtype B trees are very large, we need to split into interesting 
	#	subtrees
	#	- create UK specific subtrees
	#	-- option1: ancestral state reconstruction (Fitch algorithm)
	#	-- option2: slice by date to obtain subtrees with 1000 tips and then prune non-UK tips
	#	-- these subtrees can be small. TODO: Stephane to write avg, 50% IQR, min max.
	extract_clades
	
	#	run Source Attribution for each subtype separately (for B: each clade)
	#	requires phylodynR
	#	
	#	-- dated tree in ape format
	#	-- MH max height 20 years
	#	-- requires point estimates of PLWHIV and INCIDENCE. assumes constant numbers
	#	-- needs PLWHIV; p is proportion of people with subtype specific infection
	#	   p can be anything because the method uses the ratio NEWINF*p/ PLWIHIV*p
	#	-- runs in seconds
	#	-- output: list of donor, recip (tip index), ip (probability)
	#	   for this bootstrap replicate
	sa_by_clades

	#	run on simulated data
	#	https://github.com/mrc-ide/londonMSM_tree_simulator/blob/master/model1-simBaseline0.R
	
	#	process this output in many ways
	#	
}

seattle.180423.make.subtype.alignments<- function()
{	
	require(ape)
	
	infile.fasta		<- '~/Box Sync/OR_Work/2017/2017_Seattle/20180423/PHSKC.HXB2.March2018.fasta'
	sq					<- read.dna(infile.fasta, format='fa')
	
	#	get subtyping
	infile.st			<- '~/Box Sync/OR_Work/2017/2017_Seattle/20180423/PHSKC.HXB2.March2018.subtype.txt'	
	st					<- as.data.table(read.table(infile.st, skip=1, fill=TRUE, sep='\t'))
	setnames(st, c('V1','V2','V3','V4','V5'), c('TAXA','RESULT','MAJ_ST','BS_MAJ_ST','ORDERED_LIST'))
	st[, table(RESULT)]
	st[, table(MAJ_ST)]
	#	MAJ_ST
  	#01_AE   02_AG  06_cpx   07_BC  09_cpx  11_cpx   12_BF  13_cpx  15_01B  16_A2D 22_01A1   24_BG   29_BF   44_BF  45_cpx      A1      A2       B       C       D      F1       G 
    # 56      63       7       4       1       6       1       2       7       2       1       1       6       1       1      97       1    5636     230      17       7       5 
    #  H 
    #  2 
	
	outfile.base		<- '~/Box Sync/OR_Work/2017/2017_Seattle/20180423/PHSKC.HXB2.March2018'
	#	B + H as outgroup
	tmp					<- subset(st, MAJ_ST=='B' & !grepl('unassigned', RESULT))[, as.character(TAXA)]
	tmp					<- c(tmp, subset(st, MAJ_ST=='H')[, as.character(TAXA)])
	sqst				<- sq[tmp,]
	write.dna(sqst, file=paste0(outfile.base,'_B.fasta'))
	
	#	A1 + H as outgroup + HXB2 for referencing
	tmp					<- subset(st, MAJ_ST=='A1' & !grepl('unassigned', RESULT))[, as.character(TAXA)]
	tmp					<- c(tmp, "HXB2", subset(st, MAJ_ST=='H')[, as.character(TAXA)])
	sqst				<- sq[tmp,]
	write.dna(sqst, file=paste0(outfile.base,'_A1.fasta'))
	
	#	C + H as outgroup + HXB2 for referencing
	tmp					<- subset(st, MAJ_ST=='C' & !grepl('unassigned', RESULT))[, as.character(TAXA)]
	tmp					<- c(tmp, "HXB2", subset(st, MAJ_ST=='H')[, as.character(TAXA)])
	sqst				<- sq[tmp,]
	write.dna(sqst, file=paste0(outfile.base,'_C.fasta'))
	
	#	02_AG + H as outgroup + HXB2 for referencing
	tmp					<- subset(st, MAJ_ST=='02_AG')[, as.character(TAXA)]
	tmp					<- c(tmp, "HXB2", subset(st, MAJ_ST=='H')[, as.character(TAXA)])
	sqst				<- sq[tmp,]
	write.dna(sqst, file=paste0(outfile.base,'_02AG.fasta'))
	
	#	01_AE + H as outgroup + HXB2 for referencing
	tmp					<- subset(st, MAJ_ST=='01_AE')[, as.character(TAXA)]
	tmp					<- c(tmp, "HXB2", subset(st, MAJ_ST=='H')[, as.character(TAXA)])
	sqst				<- sq[tmp,]
	write.dna(sqst, file=paste0(outfile.base,'_01AE.fasta'))
}


seattle.170601.rm.drug.resistance.mutations<- function()
{
	require(big.phylo)
	#	handle DRMs
	infile.fasta		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.fasta'
	outfile				<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm.fasta'
	sq					<- read.dna(infile.fasta, format='fa')	
	tmp					<- which(grepl("HXB2",rownames(sq)))
	rownames(sq)[tmp]	<- 'HXB2'
	tmp					<- big.phylo:::seq.rm.drugresistance(sq)
	sq					<- tmp$nodr.seq
	write.dna(sq, file= outfile, format='fasta')
}

seattle.170607.rm.drug.resistance.mutations<- function()
{
	require(big.phylo)
	#	handle DRMs
	infile.fasta		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.plus.LANL.USA.refse.fasta'
	outfile				<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.plus.LANL.USA.refse.ndrm.fasta'	
	infile.fasta		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.June2017.PHSKCLANLe.refse.fasta'
	outfile				<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.June2017.PHSKCLANLe.refse.ndrm.fasta'
	sq					<- read.dna(infile.fasta, format='fa')	
	tmp					<- which(grepl("HXB2",rownames(sq)))
	rownames(sq)[tmp]	<- 'HXB2'
	tmp					<- big.phylo:::seq.rm.drugresistance(sq)
	sq					<- tmp$nodr.seq
	write.dna(sq, file= outfile, format='fasta')
}

seattle.170601.fastree<- function()
{	
	require(big.phylo)
	
	#	first tree
	infile.fasta	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.June2017.alignment.fasta'
	outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
	tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
	#	second tree, DRMs removed + including references for rooting
	infile.fasta	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm.fasta'
	outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
	tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
	#	run on command line
	cat(tmp)		
	
	#
	# 	bootstrap trees
	
	#	second tree, DRMs removed + including references for rooting
	infile.fasta	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm.fasta'
	bs.id			<- 2
	outfile.ft		<- gsub('\\.fasta',paste0('_ft.',sprintf("%03d",bs.id),'.newick'),infile.fasta)
	tmp				<- cmd.fasttree.one.bootstrap(infile.fasta, bs.id, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
	#	run on command line
	cat(tmp)		
	
	infile.fasta	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm.fasta'
	bs.n			<- 3
	bs.dir			<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm_bootstrap_trees'
	outfile			<- paste0('/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm_ft_bs',bs.n,'.newick')
	tmp				<- cmd.fasttree.many.bootstraps(infile.fasta, bs.dir, bs.n, outfile)
	cat(paste0('#!/bin/sh\n',tmp), file=paste0(outfile,'.sh'))
	Sys.chmod(paste0(outfile,'.sh'), mode = "777")	
}

seattle.170607.fastree<- function()
{	
	require(big.phylo)
	
	infile.fasta	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.plus.LANL.USA.refse.ndrm.fasta'
	outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
	tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma', check.binary=TRUE)
	#	run on command line
	cat(tmp)	
}

seattle.170621.fastree<- function()
{	
	require(big.phylo)
		
	#indir			<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle'
	indir			<- '/work/or105/SEATTLE'
	infile.fasta	<- file.path(indir,'Seattle.June2017.PHSKCLANLe.refse.ndrm.fasta')
	
	bs.n			<- 100	
	bs.dir			<- file.path(indir,'Seattle.June2017.PHSKCLANLe.refse.ndrm_bootstrap_trees')	
	outfile.ft		<- gsub('\\.fasta',paste0('_ft_bs',bs.n,'.newick'),infile.fasta)
	tmp				<- cmd.fasttree.many.bootstraps(infile.fasta, bs.dir, bs.n, outfile.ft, pr.args='-nt -gtr -gamma', opt.bootstrap.by='nucleotide')
	
	#	run on HPC
	cmd				<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=998, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=1, hpc.load='module load R/3.2.0')
	cmd				<- paste(cmd,tmp,sep='\n')
	cat(cmd)					
	outfile.cmd		<- paste("sea",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
	cmd.hpccaller(indir, outfile.cmd, cmd)	
}

seattle.170814.metadata.regionorigin<- function()
{
	infile	<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.2017_08_04.metadata/sequences_meta.RData'
	outfile	<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.2017_08_04.metadata/sequences_first_OR.rda'
	load(infile)
	ds	<- as.data.table(sequences_meta)
	ds	<- subset(ds, select=c(seqID, newnum, seqy, seqm))
	setnames(ds, colnames(ds), toupper(colnames(ds)))
	set(ds, NULL, 'NEWNUM', ds[, as.integer(NEWNUM)])
	set(ds, NULL, 'SEQID', ds[, as.character(SEQID)])
	set(ds, NULL, 'SEQY', ds[, as.integer(SEQY)])
	set(ds, NULL, 'SEQM', ds[, as.integer(SEQM)])
	set(ds, NULL, 'SEQ_DATE', ds[, as.numeric(SEQY)+(as.numeric(SEQM)-1+.5)/12])
	
	dsf	<- ds[, {
					z<- which(SEQ_DATE==min(SEQ_DATE))
					list(SEQID=SEQID[z], SEQ_DATE1=min(SEQ_DATE))
				}, by='NEWNUM']
	dsf[, GENE:= dsf[,gsub('(.*)-[0-9]+$','\\1',SEQID)]]	
	save(ds, dsf, file=outfile)
	
	infile	<- '~/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.2017_08_04.metadata/person.RData'
	outfile	<- gsub('\\.RData','_OR.rda',infile)
	load(infile)
	dp	<- as.data.table(person)
	dp	<- subset(dp, select=c(newnum, race, b_yr, gender, transm, birthCountry, rsh_county_name, rsa_county_name, cur_county_name, hivyr, hivmo))
	setnames(dp, colnames(dp), toupper(colnames(dp)))
	
	set(dp, NULL, 'NEWNUM', dp[, as.integer(NEWNUM)])
	set(dp, NULL, 'RACE', dp[, as.character(RACE)])
	set(dp, NULL, 'B_YR', dp[, as.integer(B_YR)])
	set(dp, NULL, 'GENDER', dp[, as.character(GENDER)])
	set(dp, NULL, 'TRANSM', dp[, as.character(TRANSM)])
	set(dp, NULL, 'BIRTHCOUNTRY', dp[, as.character(BIRTHCOUNTRY)])
	set(dp, NULL, 'RSH_COUNTY_NAME', dp[, as.character(RSH_COUNTY_NAME)])
	set(dp, NULL, 'RSA_COUNTY_NAME', dp[, as.character(RSA_COUNTY_NAME)])
	set(dp, NULL, 'CUR_COUNTY_NAME', dp[, as.character(CUR_COUNTY_NAME)])
	set(dp, NULL, 'HIVYR', dp[, as.character(HIVYR)])
	set(dp, NULL, 'HIVMO', dp[, as.character(HIVMO)])
	
	set(dp, NULL, 'RACE', dp[, gsub('\\([0-9]+\\)(.*)','\\1',RACE)])
	set(dp, NULL, 'RACE', dp[, gsub('-','_',gsub('/','_',gsub(' ','_',gsub('Not Hispanic, |, All races','',RACE))))])
	set(dp, NULL, 'TRANSM', dp[, gsub('[0-9]+\\.(.*)','\\1',TRANSM)])	
	#set(dp, NULL, 'HIV_DIAG', dp[, as.numeric(HIVYR)+(as.numeric(HIVMO)-1+.5)/12])
	
	# use line below to create the following vector
	# dp[, cat(paste(unique(BIRTHCOUNTRY), collapse='\", \"'))]
	tmp	<- c( "(USA)United States"="USA", 
			"(ZMB)Zambia"="AFR", 
			"(SLV)El Salvador"="CAR", 
			"(X99)Not Specified/Unknown"="UNKNOWN", 
			"(MEX)Mexico"="MEX", 
			"(VNM)Viet Nam"="SASIA", 
			"(TZA)Tanzania, United Republic of"="AFR", 
			"(MWI)Malawi"="AFR", 
			"(PHL)Philippines"="SASIA", 
			"(IRL)Ireland"="EU", 
			"(PER)Peru"="SAM", 
			"(ETH)Ethiopia"="AFR", 
			"(KEN)Kenya"="AFR", 
			"(PRI)Puerto Rico"="CAR", 
			"(GRC)Greece"="EU", 
			"(GTM)Guatemala"="CAR", 
			"(CUB)Cuba"="CAR", 
			"(X98)Other"="UNKNOWN", 
			"(NIC)Nicaragua"="CAR", 
			"(CIV)Cote d'Ivoire (Ivory Coast)"="AFR", 
			"(ARG)Argentina"="SAM", 
			"(ESP)Spain"="EU", 
			"(ERI)Eritrea"="AFR", 
			"(HKG)Hong Kong"="SASIA", 
			"(VEN)Venezuela"="SAM", 
			"(LBR)Liberia"="AFR", 
			"(SOM)Somalia"="AFR", 
			"(DOM)Dominican Republic"="CAR", 
			"(ITA)Italy"="EU", 
			"(THA)Thailand"="SASIA", 
			"(UGA)Uganda"="AFR", 
			"(GUF)French Guiana"="SAM", 
			"(JAM)Jamaica"="CAR", 
			"(CAN)Canada"="CAN", 
			"(DMA)Dominica"="CAR", 
			"(MRT)Mauritania"="AFR", 
			"(EGY)Egypt"="AFR", 
			"(CHN)China"="NASIA", 
			"(JPN)Japan"="NASIA", 
			"(DEU)Germany"="EU", 
			"(ZWE)Zimbabwe"="AFR", 
			"(PAN)Panama"="CAR", 
			"(TWN)Taiwan, Province of China"="NASIA", 
			"(IND)India"="NASIA", 
			"(GUM)Guam"="PAC", 
			"(KOR)Korea, Republic of (South)"="NASIA", 
			"(LAO)Lao People's Democratic Republic"="SASIA", 
			"(KHM)Cambodia"="SASIA", 
			"(GBR)United Kingdom"="EU", 
			"(ASM)American Samoa"="PAC", 
			"(FJI)Fiji"="PAC", 
			"(RUS)Russian Federation"="NASIA", 
			"(SAU)Saudi Arabia"="MIDEAST", 
			"(TGO)Togo"="AFR", 
			"(COD)Congo, Democratic Republic of Zaire"="AFR", 
			"(SLE)Sierra Leone"="AFR", 
			"(MMR)Myanmar (Burma)"="SASIA", 
			"(NGA)Nigeria"="AFR", 
			"(COL)Colombia"="SAM", 
			"(COG)Congo"="AFR", 
			"(NOR)Norway"="EU", 
			"(SDN)Sudan"="AFR", 
			"(SVK)Slovakia"="EU", 
			"(TON)Tonga"="PAC", 
			"(BRA)Brazil"="SAM", 
			"(POL)Poland"="EU", 
			"(LUX)Luxembourg"="EU", 
			"(NZL)New Zealand"="PAC", 
			"(AFG)Afghanistan"="NASIA", 
			"(MLI)Mali"="AFR", 
			"(RWA)Rwanda"="AFR", 
			"(CHL)Chile"="SAM", 
			"(BDI)Burundi"="AFR", 
			"(ECU)Ecuador"="CAR", 
			"(GHA)Ghana"="AFR", 
			"(FRA)France"="EU", 
			"(AUS)Australia"="PAC", 
			"(HND)Honduras"="CAR", 
			"(SYR)Syrian Arab Republic"="MIDEAST", 
			"(HTI)Haiti"="CAR", 
			"(MNP)Northern Mariana Islands"="PAC", 
			"(IDN)Indonesia"="SASIA", 
			"(ZAF)South Africa"="AFR", 
			"(CMR)Cameroon"="AFR", 
			"(LKA)Sri Lanka"="NASIA", 
			"(KAZ)Kazakhstan"="NASIA", 
			"(IRN)Iran, Islamic Republic of"="NASIA", 
			"(VIR)Virgin Islands, U.S."="CAR", 
			"(MHL)Marshall Islands"="PAC", 
			"(GIN)Guinea"="AFR", 
			"(MAR)Morocco"="AFR", 
			"(ROU)Romania"="EU", 
			"(YUG)Yugoslavia"="EU", 
			"(CRI)Costa Rica"="SAM", 
			"(MYS)Malaysia"="SASIA", 
			"(GMB)Gambia"="AFR", 
			"(YEM)Yemen"="MIDEAST", 
			"(SWZ)Swaziland"="AFR", 
			"(UKR)Ukraine"="EU", 
			"(AUT)Austria"="EU", 
			"(TCD)Chad"="AFR", 
			"(NLD)Netherlands"="EU", 
			"(BFA)Burkina Faso"="AFR", 
			"(URY)Uruguay"="SAM", 
			"(NAM)Namibia"="AFR", 
			"(NER)Niger"="AFR", 
			"(SGP)Singapore"="SASIA", 
			"(LBN)Lebanon"="MIDEAST", 
			"(PRT)Portugal"="EU", 
			"(BIH)Bosnia and Herzegovina"="EU",
			"(BEN)Benin"="AFR", 
			"(ISR)Israel"="MIDEAST", 
			"(HUN)Hungary"="EU", 
			"(MDA)Moldova, Republic of"="EU", 
			"(BGR)Bulgaria"="EU", 
			"(SEN)Senegal"="AFR", 
			"(DNK)Denmark"="EU", 
			"(LSO)Lesotho"="AFR", 
			"(BHR)Bahrain"="MIDEAST", 
			"(ALB)Albania"="EU", 
			"(UZB)Uzbekistan"="NASIA", 
			"(TTO)Trinidad and Tobago"="CAR", 
			"(BWA)Botswana"="AFR", 
			"(PAK)Pakistan"="NASIA", 
			"(IRQ)Iraq"="NASIA", 
			"(DJI)Djibouti"="AFR", 
			"(BOL)Bolivia"="SAM", 
			"(ATG)Antigua and Barbuda"="CAR", 
			"(MOZ)Mozambique"="AFR", 
			"(GAB)Gabon"="AFR", 
			"(BEL)Belgium"="EU", 
			"(BMU)Bermuda"="CAR", 
			"(LVA)Latvia"="EU", 
			"(LBY)Libyan Arab Jamahiriya"="AFR", 
			"(FSM)Micronesia, Federated States of"="PAC", 
			"(MNG)Mongolia"="NASIA", 
			"(SWE)Sweden"="EU", 
			"(TUR)Turkey"="MIDEAST", 
			"(GEO)Georgia"="EU", 
			"(WSM)Samoa"="PAC", 
			"(AGO)Angola"="AFR", 
			"(PNG)Papua New Guinea"="SASIA", 
			"(GUY)Guyana"="SAM", 
			"(AZE)Azerbaijan"="NASIA", 
			"(NPL)Nepal"="NASIA", 
			"(CZE)Czech Republic"="EU", 
			"(DZA)Algeria"="AFR", 
			"(PRY)Paraguay"="SAM", 
			"(QAT)Qatar"="MIDEAST", 
			"(BLZ)Belize"="SAM", 
			"(MTQ)Martinique"="CAR", 
			"(FIN)Finland"="EU", 
			"(BRB)Barbados"="CAR", 
			"(ISL)Iceland"="EU", 
			"(PRK)Korea, Dem People's Rep of (North)"="NASIA", 
			"(CHE)Switzerland"="EU", 
			"(UMI)U.S. Minor Outlying Areas"="PAC")
	dp	<- merge(dp, data.table(BIRTHCOUNTRY=names(tmp), BIRTHLOC=tmp), all.x=1, by='BIRTHCOUNTRY')
	set(dp, dp[, which(is.na(BIRTHLOC))], 'BIRTHLOC', 'UNKNOWN')
	save(dp, file=outfile)
}

seattle.170814.LSD<- function()
{	
	require(data.table)
	require(big.phylo)
	#
	#	write the overall dates file
	#
	if(0)
	{
		infile.seq	<- '~/Box Sync/OR_Work/2017/2017_Seattle/Seattle.2017_08_04.metadata/sequences_first_OR.rda'
		infile.pers	<- "~/Box Sync/OR_Work/2017/2017_Seattle/Seattle.2017_08_04.metadata/person_OR.rda"
		infile.tree	<- '~/Box Sync/OR_Work/2017/2017_Seattle/Seattle.plus.LANL.USA/Seattle.plus.LANL.USA.refse.ndrm_ft_reroot.newick'
		load(infile.seq)
		load(infile.pers)
		
		ph			<- read.tree(infile.tree)
		dph			<- data.table(TAXA=ph$tip.label)
		dph[, SEA:= as.numeric(grepl('^SEA',TAXA))]
		tmp			<- dph[, which(SEA==1)]
		set(dph, tmp, 'NEWNUM', as.integer(dph[tmp, gsub('SEATTLE\\.([0-9]+)\\|.*','\\1',TAXA)]))
		#	there is one sequence per individual
		dph			<- merge(dph, unique(subset(dsf, select=c(NEWNUM, SEQ_DATE1))), by='NEWNUM', all.x=TRUE)
		
		write.csv(subset(dph, !is.na(NEWNUM) & is.na(SEQ_DATE1)), file='~/Box Sync/OR_Work/2017/2017_Seattle/Seattle.plus.LANL.USA/Seattle.plus.LANL.USA.refse.ndrm_ft_reroot_missingsamplingdates.csv')
		
		
		infile.dates	<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo_OR.csv'
		df				<- as.data.table(read.csv(infile.dates))
		df				<- subset(df, select=c(ID, SUBTYPE, SAMPLING_DATE))
		setnames(df, c('ID','SAMPLING_DATE'), c('TAXA','DATE'))
		tmp				<- copy(df)
		set(tmp, NULL, 'TAXA', tmp[,gsub('[0-9]$','',paste0(TAXA,'_subtype',SUBTYPE))])
		df				<- rbind(df, tmp)	
		infile.dates	<- gsub('\\.csv','_lsddates.csv',infile.dates)
		write.csv(df, row.names=FALSE, infile.dates)
		#
		#	write subtype specific dates files
		#
		indir.ft		<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/trees_ft'
		infiles			<- data.table(F=list.files(indir.ft, pattern='*newick$', full.names=TRUE))
		infiles[, DATES_FILE:= gsub('_ft_bs100\\.newick','_lsddates.csv',F)]
		invisible(infiles[, {
							cmd		<- cmd.lsd.dates(infile.dates, F, DATES_FILE, run.lsd=FALSE)
							system(cmd)
							NULL
						}, by='F'])
	}	
	#
	#	run LSD 
	#
	#indir.ft		<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/trees_ft'
	#indir.dates		<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/trees_ft'
	#outdir			<- '~/Dropbox (SPH Imperial College)/2017_NL_Introductions/trees_ft'
	
	indir.ft		<- '/work/or105/ATHENA_2016/vlintros/trees_ft'
	indir.dates		<- '/work/or105/ATHENA_2016/vlintros/trees_lsd'
	outdir			<- '/work/or105/ATHENA_2016/vlintros/trees_lsd'
	
	ali.len			<- 1000				
	lsd.args		<- '-v 1 -c -b 10'			# no rooting, no re-estimation of rates				
	#	get files	
	infiles	<- data.table(F=list.files(indir.ft, pattern='*newick$', full.names=TRUE, recursive=TRUE))
	infiles	<- subset(infiles, !grepl('RAxML',F))
	infiles[, SUBTYPE:= gsub('_.*','',basename(F))]
	tmp		<- data.table(FD=list.files(indir.dates, pattern='*_lsddates.csv$', full.names=TRUE, recursive=TRUE))
	tmp[, SUBTYPE:= gsub('_.*','',basename(FD))]
	infiles	<- merge(infiles, tmp, by='SUBTYPE')
	infiles[, FL:= file.path(outdir,gsub('_ft\\.|_ft_','_lsd_',gsub('\\.newick','',basename(F))))]	
	#	build LSD commands
	dlsd	<- infiles[, {
				cmd		<- cmd.lsd(F, FD, ali.len, outfile=FL, pr.args=lsd.args)
				list(CMD=cmd)
			}, by='F']
	tmp		<- dlsd[, paste(CMD, collapse='\n')]
	#	qsub
	cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=50, hpc.q="pqeelab", hpc.mem="5800mb",  hpc.nproc=1, hpc.load='module load R/3.3.3')							
	cmd		<- paste(cmd,tmp,sep='\n')
	#cat(cmd)					
	outfile	<- paste("scRAr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
	cmd.hpccaller(indir.ft, outfile, cmd)	
}


seattle.170601.clustering<- function()
{	
	require(big.phylo)
	#	first tree
	infile.ft		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm_ft.newick'
	ph				<- read.tree(infile.ft)
	#	reroot at AY371169
	tmp				<- which(grepl('AY371169',ph$tip.label))
	ph				<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
	ph 				<- ladderize( ph )
	#	process bootstrap support values	
	tmp				<- as.numeric( ph$node.label )		
	tmp[is.na(tmp)]	<- 0
	ph$node.label	<- tmp
	#	
	write.tree(ph, gsub('\\.newick','_reroot.newick',infile.ft))
	pdf(gsub('newick','pdf',infile.ft), w=15, h=350)
	plot(ph, cex=0.4, show.node.label=TRUE)
	dev.off()
	#
	#read patristic distances
	require(hivclust)
	stat.fun						<- hivc.clu.min.transmission.cascade	
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)	
	thresh.brl						<- 0.045
	thresh.bs						<- 0.8
	#produce clustering 
	clustering						<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph$node.label, retval="all")	
	outfile							<- gsub('\\.newick',paste0('_reroot.clu',thresh.brl*1e2,'-',thresh.bs*1e2,'.pdf'),infile.ft)	
	hivc.clu.plot(ph, 
					clustering[["clu.mem"]], 
					file=outfile, 
					pdf.off=0, pdf.width=15, pdf.height=1000, 
					show.node.label=TRUE, show.edge.label=TRUE,  
					cex.edge.incluster=2)
	dev.off()
	
	dfc		<- data.table(TAXA=ph$tip.label, CLU=clustering$clu.mem[1:Ntip(ph)])
	outfile	<- gsub('\\.newick',paste0('_reroot.clu',thresh.brl*1e2,'-',thresh.bs*1e2,'.rda'),infile.ft)
	save(ph, clustering, dfc, file, file=outfile)	
}

seattle.170607.clustering<- function()
{	
	require(big.phylo)
	#	first tree
	infile.ft		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.plus.LANL.USA.refse.ndrm_ft.newick'
	ph				<- read.tree(infile.ft)
	#	reroot at AY371169
	tmp				<- which(grepl('AY371169',ph$tip.label))
	ph				<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
	ph 				<- ladderize( ph )
	#	process bootstrap support values	
	tmp				<- as.numeric( ph$node.label )		
	tmp[is.na(tmp)]	<- 0
	ph$node.label	<- tmp
	#	
	write.tree(ph, gsub('\\.newick','_reroot.newick',infile.ft))
	pdf(gsub('\\.newick','_reroot.pdf',infile.ft), w=15, h=2e3)
	plot(ph, cex=0.4, show.node.label=TRUE)
	dev.off()
	#
	#read patristic distances
	require(hivclust)
	stat.fun						<- hivc.clu.min.transmission.cascade	
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)	
	thresh.brl						<- 0.045
	thresh.bs						<- 0.8
	#produce clustering 
	clustering						<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph$node.label, retval="all")	
	outfile							<- gsub('\\.newick',paste0('_reroot.clu',thresh.brl*1e2,'-',thresh.bs*1e2,'.pdf'),infile.ft)	
	hivc.clu.plot(ph, 
			clustering[["clu.mem"]], 
			file=outfile, 
			pdf.off=0, pdf.width=15, pdf.height=1000, 
			show.node.label=TRUE, show.edge.label=TRUE,  
			cex.edge.incluster=2)
	dev.off()
	
	dfc		<- data.table(TAXA=ph$tip.label, CLU=clustering$clu.mem[1:Ntip(ph)])
	outfile	<- gsub('\\.newick',paste0('_reroot.clu',thresh.brl*1e2,'-',thresh.bs*1e2,'.rda'),infile.ft)
	save(ph, clustering, dfc, file, file=outfile)	
}

seattle.170601.clustering.plot<- function()
{
	
	require(data.table)
	infile	<- "/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm_ft_reroot.clu4.5-80.rda"
	load(infile)
	#	exclude reference sequences non-Seattle
	dfc		<- subset(dfc, !grepl('AY371169|HXB2',TAXA))
	#	read meta-data 
	dfc[, ID:= sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',1)]
	dfc[, SXO:= sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',2)]
	dfc[, RACE:= sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',3)]
	dfc[, SEQ_YR:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',4))]
	dfc[, SEQ_Q:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',5))]
	dfc[, HIV_YR:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',6))]
	dfc[, BIRTH_YR:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',7))]
	dfc[, AGE_AT_SEQ:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',8))]
	dfc[, AGE_AT_HIV:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',9))]
	#	fixup
	set(dfc, dfc[, which(SXO=='UNKNOWN')], 'SXO', NA_character_)
	set(dfc, dfc[, which(SXO=='')], 'SXO', NA_character_)
	set(dfc, dfc[, which(SXO=='BLOOD')], 'SXO', 'Other')
	set(dfc, dfc[, which(SXO=='PERINAT')], 'SXO', 'Other')
	set(dfc, dfc[, which(RACE=='FALSE')], 'RACE', NA_character_)
	#	give each singleton a new cluster ID (cluster of size 1)
	tmp	<- dfc[, which(is.na(CLU))]
	set(dfc, tmp, 'CLU', max(dfc$CLU, na.rm=TRUE)+seq_along(tmp))
	
	#	compute cluster statistics 
	tmp	<- dfc[, list(	CLU_TS= min(HIV_YR, na.rm=TRUE),
				CLU_TE= max(HIV_YR, na.rm=TRUE),
				CLU_TM=mean(HIV_YR, na.rm=TRUE),
				CLU_N=length(ID),
				CLU_N_MSM=length(which(grepl('MSM',SXO))),
				CLU_N_IDU=length(which(grepl('IDU',SXO))),
				CLU_N_HET=length(which(grepl('HET',SXO)))
				), by='CLU']
	dfc	<- merge(dfc, tmp, by='CLU') 
	#	determine max lkl SXO of cluster
	tmp	<- melt(dfc, id.vars='CLU',measure.vars=c('CLU_N_MSM','CLU_N_IDU','CLU_N_HET'))
	tmp	<- tmp[, list(CLU_MX_SXO=variable[which.max(value)]), by='CLU']
	set(tmp, NULL, 'CLU_MX_SXO', tmp[, gsub('CLU_N_','',CLU_MX_SXO)])
	dfc	<- merge(dfc, tmp, by='CLU')
	#	set factors 
	set(dfc, NULL, 'CLU_MX_SXO', dfc[, factor(CLU_MX_SXO, levels=c('MSM','HET','IDU'), labels=c('MSM','heterosexual','injecting\ndrug\nuser'))])
	set(dfc, dfc[, which(is.na(SXO))], 'SXO', 'Unknown')
	set(dfc, dfc[, which(SXO=='HETEROp')], 'SXO', 'HETERO')
	set(dfc, NULL, 'SXO', dfc[, factor(SXO,  levels=c('MSM','MSM_IDU','IDU','HETERO','Other','Unknown'), 
											 labels=c('MSM','MSM/injecting drug user','injecting drug user','heterosexual','other','unknown'))])
	#
	#	plot clusters only
	dfp	<- subset(dfc, CLU_N>1 & !is.na(CLU_TS))
	#	order by CLU_MX_SXO and start date of cluster
	tmp	<- unique(subset(dfp, select=c(CLU, CLU_MX_SXO, CLU_TS)))
	tmp	<- tmp[order(CLU_MX_SXO, CLU_TS), ]	
	set(tmp, NULL, 'CLU_ORDERED', tmp[, factor(CLU, levels=tmp$CLU)])
	dfp	<- merge(dfp, subset(tmp, select=c(CLU, CLU_ORDERED)), by='CLU')
	tmp	<- melt(unique(subset(dfp, select=c(CLU_ORDERED, CLU_TS, CLU_TE, CLU_MX_SXO))), id.vars=c('CLU_ORDERED','CLU_MX_SXO'))
	ggplot(dfp, aes(y=CLU_ORDERED)) +
			scale_colour_brewer(palette='Set1') +
			scale_x_continuous(breaks=seq(1950,2030,5), expand=c(0,0)) +
			geom_line(data=tmp, aes(x=value, y=CLU_ORDERED), colour='black') +
			geom_point(aes(x=HIV_YR, colour=SXO)) +
			facet_grid(CLU_MX_SXO~., scales='free_y', space='free_y') +
			labs(x='\nyear of diagnosis', y='phylogenetic clusters\n', colour='sexual orientation') +
			theme_bw() +
			theme(	strip.text.y=element_text(angle=0),
					panel.grid.minor.x=element_line(colour="grey70", size=0.25), 
					panel.grid.major.x=element_line(colour="grey70", size=0.5),
					panel.grid.minor.y=element_blank(), 
					panel.grid.major.y=element_blank(),					
					axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					legend.position='bottom')
	ggsave(gsub('\\.rda','_cluevolution_by_sxo.pdf',infile), w=12,h=18)	
}