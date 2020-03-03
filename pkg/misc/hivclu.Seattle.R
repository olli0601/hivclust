## ---- rmd.chunk.seattle.wrapper ----
seattle.start.HPC<- function()
{
	CODE.HOME	<<- "/rds/general/user/or105/home/libs/hivclust/pkg"
	HOME		<<- '/rds/general/user/or105/home' 
	
	#	various   
	if(0) 
	{				
		#hpc.load	<- "module load anaconda3/personal"		
		hpc.load	<- "module load R/3.3.3"
		hpc.select	<- 1						# number of nodes
		hpc.nproc	<- 1						# number of processors on node
		hpc.walltime<- 123						# walltime 
		hpc.q		<- "pqeelab"				# PBS queue
		hpc.mem		<- "6gb" 					# RAM		
		pbshead		<- "#!/bin/sh"
		tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
		pbshead		<- paste(pbshead, tmp, sep = "\n")
		tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
		pbshead 	<- paste(pbshead, tmp, sep = "\n")
		pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
		if(!is.na(hpc.q)) 
			pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
		pbshead 	<- paste(pbshead, hpc.load, sep = "\n")			
		tmp			<- paste('Rscript ',file.path(CODE.HOME, "misc/hivclu.startme.R"), ' -exe=VARIOUS', '\n', sep='')
		cmd			<- paste(pbshead,tmp,sep='\n')
		cat(cmd)	 								
		outfile		<- gsub(':','',paste("pv",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
		outfile		<- file.path(HOME, outfile)
		cat(cmd, file=outfile)
		cmd 		<- paste("qsub", outfile)		
		cat(system(cmd, intern= TRUE))	
		quit("no")	
	}	
	#	various job array
	if(1) 
	{		
		cmds		<- paste0('Rscript ',file.path(CODE.HOME, "misc/hivclu.startme.R"), ' -exe=VARIOUS', ' -input=', 1:4999, '\n')
						
		#	make PBS header
		hpc.load	<- "module load anaconda3/personal\nsource activate base_backup2"
		#hpc.load	<- "module load R/3.3.3"
		hpc.select	<- 1						# number of nodes
		hpc.nproc	<- 1						# number of processors on node
		hpc.walltime<- 23						# walltime
		hpc.q		<- NA #"pqeelab"				# PBS queue
		hpc.mem		<- "2gb" 					# RAM
		hpc.array	<- length(cmds)	# number of runs for job array
		pbshead		<- "#!/bin/sh"
		tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
		pbshead		<- paste(pbshead, tmp, sep = "\n")
		tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
		pbshead 	<- paste(pbshead, tmp, sep = "\n")
		pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
		if(!is.na(hpc.array))
			pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')		
		if(!is.na(hpc.q)) 
			pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
		pbshead 	<- paste(pbshead, hpc.load, sep = "\n")			
		#	make array job
		cmd <- sapply(seq_along(cmds), function(i) paste0(i,')\n',cmds[i],';;\n'))
		cmd <- paste0('case $PBS_ARRAY_INDEX in\n',paste0(cmd, collapse=''),'esac')
		cmd		<- paste(pbshead,cmd,sep='\n')	
		#	submit job
		outfile		<- gsub(':','',paste("phy",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
		outfile		<- file.path(HOME, outfile)
		cat(cmd, file=outfile)
		cmd 		<- paste("qsub", outfile)
		cat(cmd)
		cat(system(cmd, intern= TRUE))
	} 
}

## ---- rmd.chunk.seattle.various ----
seattle.various<- function()
{
	#seattle.191017.phydyn.volz.msmUK.mle()
	#seattle.191017.phydyn.olli.SIT01.sim()
	#seattle.191017.phydyn.olli.SITmf01.sim()
	#seattle.191017.phydyn.olli.SITmf01.longsim()
	#seattle.191017.phydyn.olli.SITmf01.add.true.ll()
	#seattle.191017.phydyn.olli.SITmf01.mle()	
	#seattle.191017.phydyn.olli.SITmf01yrsXX.sim()
	seattle.191017.phydyn.olli.SITmf01yrsXX.mle()
}


## ---- rmd.chunk.seattle.191017.phyloscanner.subgraph.empirics ----
seattle.191017.phyloscanner.subgraph.empirics<- function()
{
	require(data.table) 
	require(phangorn)
	require(ggplot2)
	require(reshape)
	
	#	working directory with phyloscanner output		
	home <- '~/Box Sync/OR_Work/Seattle/analysis_191017'
	indir.phsc	<- file.path(home,'phyloscanner_dated')
	infile.meta <- file.path(home,'misc','180709_sequence_labels.rda')	
	
	#
	#	extract subgraph taxa	
	#
	infiles <- data.table(F=list.files(indir.phsc, pattern='rda$', full.names=TRUE, recursive=TRUE))	
	infiles <- subset(infiles, grepl('datedsubgraphs',F))
	#	SELECT defines if KC, KCHSX, KCMSM
	infiles[, SELECT:= gsub('^(.*)subgraphs_([A-Z]+)\\.rda','\\2',basename(F))]
	#	for subtype B, I ran separate analyses based on very large subtrees for comp efficiency
	#	ST stores the subtype, and ST_CLADE the large subtree
	infiles[, ST:= gsub('^.*Subtype([^_]+)_.*\\.rda','\\1',basename(F))]
	infiles[, ST_CLADE:= as.integer(gsub('[^c]*c([0-9]+)','\\1',ST))]
	infiles[, ST:= gsub('([^c]*)c([0-9]+)','\\1',ST)]
	#	ID of bootstrap replicate, 000 for analysis on real alignment
	infiles[, REP:= gsub('^.*ndrm_([0-9]+)_.*\\.rda','\\1',basename(F))]
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
			}, by=c('ST','ST_CLADE','REP','SELECT')]	
	
	#	add meta data from taxa names
	regex.tip.label <- '^([A-Z]+)_+PR/RT-([0-9]+)_([0-9]+)_([a-zA-Z]+)_([a-zA-Z]+)_([a-zA-Z]+)_([a-zA-Z]+)_([a-zA-Z_]+)_([0-9]+)$'
	dsubgraphtaxa[, ID:= as.numeric(gsub(regex.tip.label,'\\3',TAXA))]	
	dsubgraphtaxa[, SEQ_YEAR:= as.numeric(gsub(regex.tip.label,'\\9',TAXA))]
	dsubgraphtaxa[, FULL_NAME:= paste0(ST,'_',REP,'_',SELECT,'_',NAME)]
	
	#	add meta data from pre-processed persons file	
	load(infile.meta)
	dind <- as.data.table(dind)	
	setnames(dind, c('newnum','b_yr','birthCountry2','county','race2','transm2','Gender2'), c('ID','BIRTH_YEAR2','BIRTH_COUNTRY2','COUNTY','RACE2','TRANSM2','GENDER2'))
	tmp <- subset(dind, select=c(ID,BIRTH_YEAR2,BIRTH_COUNTRY2,COUNTY,RACE2,TRANSM2,GENDER2))
	set(tmp, NULL, 'BIRTH_YEAR2', tmp[, as.numeric(BIRTH_YEAR2)])
	dsubgraphtaxa <- merge(dsubgraphtaxa,tmp,by='ID')
	dsubgraphtaxa[, AGE2020:= cut(2020-BIRTH_YEAR2, breaks=c(-Inf,20,25,30,35,40,45,50,55,60,Inf))]
		
	#	subgraph size 
	dsubgraphsize <- dsubgraphtaxa[, list(SIZE=length(ID)), by=c('ST','REP','SELECT','NAME')]
	dsubgraphsize[, SINGLETON:= as.character(factor(SIZE==1, levels=c(TRUE,FALSE), labels=c('subgraph size 1','subgraph size >1')))]
	ggplot(subset(dsubgraphsize, SELECT!='KC'), aes(x=SIZE, fill=ST)) + 
			geom_bar() + 
			facet_grid(SELECT~.) +
			labs(x='\nsize of phylogenetic subgraphs attributed to King County/Seattle')
	
	#	number of individuals in subgraphs of size == 1 or size > 1
	dsubgraphind <-  dsubgraphsize[, list(N_INDIVIDUALS=sum(SIZE)), by=c('ST','REP','SELECT','SINGLETON')]
	dsubgraphind <- dcast.data.table(dsubgraphind, ST+REP+SINGLETON~SELECT, value.var='N_INDIVIDUALS')
	
	
	tmp <- subset(dsubgraphsize, SIZE>1, c(ST, REP, SELECT, NAME))
	dsubgraphtaxa2 <- merge(dsubgraphtaxa, tmp, by=c('ST','REP','SELECT','NAME'))
	
	# ethnicity
	dsubgraphtaxa2[, DUMMY:= RACE2]
	dcnt <- dsubgraphtaxa2[, list(N= length(ID)), by=c('FULL_NAME','SELECT','ST','DUMMY')]
	tmp <- dsubgraphtaxa2[, list(TOTAL= length(ID)), by=c('FULL_NAME','SELECT','ST')]
	dcnt <- merge(dcnt, tmp, by=c('FULL_NAME','SELECT','ST'))
	ggplot(subset(dcnt, SELECT!='KC'), aes(x=FULL_NAME, y=N, fill=DUMMY)) +
			geom_bar(stat='identity') + 
			labs(y='#individuals',x='phylo subgraphs',fill='ethnicity') +
			theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
			facet_grid(SELECT~ST, scale='free_x', space='free_x')
	tmp <- dcnt[, list(R= sum( (N/TOTAL)^2), PMAX= max(N/TOTAL) ), by=c('FULL_NAME','SELECT','ST')]
	tmp[, list(R_AVG=mean(R), PMAX_AVG=mean(PMAX)), by=c('SELECT')]
	
	#ggsave(file=file.path(outdir, 'empirical_subgraph_composition_ethnicity.pdf'), h=15,w=6)
	
	# birth place
	dsubgraphtaxa2[, DUMMY:= as.character(factor(BIRTH_COUNTRY2=='USA',levels=c(TRUE,FALSE),labels=c('US-born','other')))]
	dcnt <- dsubgraphtaxa2[, list(N= length(ID)), by=c('FULL_NAME','SELECT','ST','DUMMY')]
	tmp <- dsubgraphtaxa2[, list(TOTAL= length(ID)), by=c('FULL_NAME','SELECT','ST')]
	dcnt <- merge(dcnt, tmp, by=c('FULL_NAME','SELECT','ST'))
	ggplot(subset(dcnt, SELECT!='KC'), aes(x=FULL_NAME, y=N, fill=DUMMY)) +
			geom_bar(stat='identity') + 
			labs(y='#individuals',x='phylo subgraphs',fill='birth place') +
			theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
			facet_grid(SELECT~ST, scale='free_x', space='free_x')
	tmp <- dcnt[, list(R= sum( (N/TOTAL)^2), PMAX= max(N/TOTAL) ), by=c('FULL_NAME','SELECT','ST')]
	tmp[, list(R_AVG=mean(R), PMAX_AVG=mean(PMAX)), by=c('SELECT')]
	
	#ggsave(file=file.path(outdir, 'empirical_subgraph_composition_birthplace.pdf'), h=15,w=6)

	
	# age2020
	dsubgraphtaxa2[, DUMMY:= AGE2020]
	dcnt <- dsubgraphtaxa2[, list(N= length(ID)), by=c('FULL_NAME','SELECT','ST','DUMMY')]
	tmp <- dsubgraphtaxa2[, list(TOTAL= length(ID)), by=c('FULL_NAME','SELECT','ST')]
	dcnt <- merge(dcnt, tmp, by=c('FULL_NAME','SELECT','ST'))
	ggplot(subset(dcnt, SELECT!='KC'), aes(x=FULL_NAME, y=N, fill=DUMMY)) +
			geom_bar(stat='identity') + 
			labs(y='#individuals',x='phylo subgraphs',fill='age in 2020') +
			theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
			facet_grid(SELECT~ST, scale='free_x', space='free_x')
	tmp <- dcnt[, list(R= sum( (N/TOTAL)^2), PMAX= max(N/TOTAL) ), by=c('FULL_NAME','SELECT','ST')]
	tmp[, list(R_AVG=mean(R), PMAX_AVG=mean(PMAX)), by=c('SELECT')]
	
}

## ---- rmd.chunk.seattle.190723.test.simmap.casting ----
seattle.190723.test.simmap.casting <- function()
{
	require(data.table)
	require(tidyverse)
	require(ape)
	require(phytools)
	require(phangorn)
	require(phyloscannerR)
	
	infile <- "/Users/or105/Box/OR_Work/Seattle/analysis_190723/phyloscanner_HSX/Subtype_01_AE_0_workspace.rda"
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

## ---- rmd.chunk.seattle.191017.phydy.select.taxa ----
seattle.191017.phydy.select.taxa <- function()
{
	require(data.table)
	require(tidyverse)
	require(ape)
	require(phytools)
	require(phangorn)
	require(phyloscannerR)
	
	home <- '~/Box Sync/OR_Work/Seattle/analysis_191017'
	indir.phsc	<- file.path(home,'phyloscanner_dated')
	outdir <- file.path(home,'phydy')
	infiles.phsc <- data.table(F=list.files(indir.phsc, pattern='_annotated_dated_tree.rda$', full.names=TRUE, recursive=TRUE))
	infiles.phsc[, REP:= as.numeric(gsub('^.*_ndrm_([0-9]+)_ft_.*$','\\1',basename(F)))]
	infiles.phsc[, ST:= gsub('^.*_Subtype([0-9A-Za-z]+)_.*$','\\1',basename(F))]
	infiles.phsc[, LOCAL_STATE:= gsub('^.*_rerooted_([A-Za-z]+)_.*$','\\1',basename(F))]
	
	control	<- list( 	localstate.regex='^KCMSM.*',
						#localstate.regex='^KCHSX.*',
						bs.replicate=0,
						keep.n.nonlocal.tips=2,
						remove.localstate.subtrees.of.size.1= 1
						)
	verbose <- 1
	
	localstate.regex <- control$localstate.regex
	localstate <- gsub('[^A-Za-z]','',localstate.regex)
	keep.n.nonlocal.tips <- control$keep.n.nonlocal.tips
	remove.localstate.subtrees.of.size.1 <- control$remove.localstate.subtrees.of.size.1
	bs.replicate <- control$bs.replicate
	infiles.phsc.idx <- infiles.phsc[, which(REP==bs.replicate & LOCAL_STATE==localstate)]
	#	process all subtypes
	phs <- list()
	while(length(infiles.phsc.idx))
	{
		#i<- 197
		i <- infiles.phsc.idx[1]
		#	read dated tree
		if(verbose)
			cat('\nProcessing ',infiles.phsc[i,F])		
		infile <- infiles.phsc[i,F]
		load(infile)
		if(verbose)
			cat('\nIs binary ',is.binary(ph), ' has zero edge len ', any(ph$edge.length<1e-12))
		
		#	for each tip in local state, determine closest tips that are not in local state
		#	option: either keep all of them, or those in subtrees of size >1
		if(remove.localstate.subtrees.of.size.1==0)
		{
			localstate.tips.idx <- which(grepl(localstate.regex,ph$tip.label))	
		}
		if(remove.localstate.subtrees.of.size.1>0)
		{
			subtrees.of.tips <- data.table(	TIP=1:Ntip(ph), 
					ST=as.character(attr(ph, "SPLIT")[1:Ntip(ph)])
			)
			subtrees.of.tips <- subset( subtrees.of.tips, grepl(localstate.regex,ST) )								
			subtrees.of.tips <- merge(subtrees.of.tips, subtrees.of.tips[, list(N=length(TIP)),by='ST'], by='ST')
			localstate.tips.idx <- subset(subtrees.of.tips, N>1)[, TIP]		
		}
		if(verbose) 
			cat('\nFound tips in local state, n=',length(localstate.tips.idx))
		
		if(length(localstate.tips.idx)>0)
		{
			#	get subtrees of grand-parent of MRCA. these should contain at least 
			#	one tip in non-local state, often a few more
			#	keep the n closest, which should contain "embedded" non-local tips
			#	that may indicate implausible phylogeographic reconstructions
			mrcas <- which(attr(ph,'SUBGRAPH_MRCA'))
			subtrees.of.tips <- unique(as.character(attr(ph,'SPLIT')[localstate.tips.idx]))	
			mrcas <- mrcas[ as.character(attr(ph, 'SPLIT')[mrcas])%in%subtrees.of.tips ]
			mrcas2 <- sapply(mrcas, function(mrca){
						z <- Ancestors(ph, mrca)
						z[min(2,length(z))]
					})
			subtrees <- lapply(mrcas2, function(x) extract.clade(ph, x))
			nonlocal.tips <- lapply(subtrees, function(x){
						subtree.dist <- as.data.table(melt(cophenetic.phylo(x)))
						setnames(subtree.dist, c('Var1','Var2','value'), c('TAXA1','TAXA2','D'))
						set(subtree.dist, NULL, 'TAXA1', as.character(subtree.dist$TAXA1))
						set(subtree.dist, NULL, 'TAXA2', as.character(subtree.dist$TAXA2))
						subtree.dist <- subset(subtree.dist, TAXA1!=TAXA2 & grepl(localstate.regex, TAXA1) & !grepl(localstate.regex, TAXA2))
						subtree.dist <- subtree.dist[order(D),]
						nonlocal.tips <- subtree.dist[1:min(nrow(subtree.dist),keep.n.nonlocal.tips),TAXA2]
						nonlocal.tips
					})
			nonlocal.tips <- unique(unlist(nonlocal.tips))
			#	ensure that "encompassing" non-local tips are also included
			siblings <- Siblings(ph, mrcas, include.self=FALSE)
			nonlocal.tips2 <- lapply(siblings, function(x){
						if(x<=Ntip(ph))
						{
							nonlocal.tips2 <- ph$tip.label[x]
						}
						if(x>Ntip(ph))
						{
							tmp <- extract.clade(ph,x)
							tmp <- data.table(TAXA=tmp$tip.label, DEPTH=node.depth.edgelength(tmp)[1:Ntip(tmp)])
							tmp <- subset(tmp, !grepl(localstate.regex, TAXA))
							tmp <- tmp[order(DEPTH),]
							nonlocal.tips2 <- tmp[1:min(nrow(tmp), keep.n.nonlocal.tips), TAXA]		
						}
						nonlocal.tips2
					})
			nonlocal.tips <- unique(c(nonlocal.tips, unlist(nonlocal.tips2)))
			
			#	drop all other tips  	
			tmp <- unique(c(ph$tip.label[localstate.tips.idx],nonlocal.tips))
			drop.nonlocal.tips <- setdiff(ph$tip.label,tmp)
			if(verbose)
				cat('\nDropping non local tips from phylogeny, n=', length(drop.nonlocal.tips))
			ph <- phyloscanner.to.simmap(ph)	
			ph <- phytools:::drop.tip.simmap(ph, drop.nonlocal.tips)
			ph <- simmap.to.phyloscanner(ph)
			if(verbose)
				cat('\nRetained tips from phylogeny, n=', Ntip(ph))
			
			if(verbose)
				cat('\nIs binary ',is.binary(ph), ' has zero edge len ', any(ph$edge.length<1e-12))
						
			# 	plot and save
			outfile <- gsub('annotated_dated_tree\\.rda',paste0('annotated_dated_collapsed_tree_',gsub('[^A-Za-z0-9]','',paste(control, collapse='')),'.pdf'),basename(infile))
			outfile <- file.path(outdir,outfile)
			tmp <- vector('list')
			tmp[['tree']] <- ph
			tmp[['tree']][['node.label']] <- NULL 
			tmp[['tree']][['node.states']] <- tmp[['tree']][['mapped.edge']] <- tmp[['tree']][['maps']] <- NULL
			attr(tmp[['tree']],'map.order') <- NULL
			attr(tmp[['tree']],'class') <- 'phylo'
			tmp[['read.counts']] <- rep(1, Ntip(ph))	
			write.annotated.tree(tmp, outfile, format="pdf", pdf.scale.bar.width = 5, pdf.w = 40, pdf.hm = 0.4, verbose = FALSE)
			outfile <- gsub('annotated_dated_tree\\.rda',paste0('annotated_dated_collapsed_tree_',gsub('[^A-Za-z0-9]','',paste(control, collapse='')),'.rda'),basename(infile))
			outfile <- file.path(outdir,outfile)
			save(ph, file=outfile)
			
			phs[[length(phs)+1L]] <- ph
		}
		infiles.phsc.idx <- infiles.phsc.idx[-1]		
	}
	
	#	collect node depths to root
	node.depths <- sapply(seq_along(phs), function(i) max(node.depth.edgelength(phs[[i]])[1:Ntip(phs[[i]])]) )
	node.depths.max <- max(node.depths)
		
	#	strip attributes and add root edge to preserve dating 
	for(i in seq_along(phs))
	{
		ph <- phs[[i]]
		ph$maps <- ph$mapped.edge <- ph$node.states<- ph$node.label <- NULL
		attr(ph,'SPLIT') <- attr(ph,'INDIVIDUAL') <-  attr(ph,'SUBGRAPH_MRCA') <-  attr(ph,'map.order') <- NULL
		attr(ph, "class") <- 'phylo'
		ph$root.edge <- node.depths.max - node.depths[i] + 1
		phs[[i]] <- ph
	}
	
	# concatenate trees as binary tree
	#options(expressions=1e4)
	ph <- eval(parse(text=paste('phs[[',seq_along(phs),']]', sep='',collapse='+')))
	#options(expressions=5e3)
	ph <- multi2di(ph, random=TRUE)
	outfile <- file.path( outdir, paste0('SubtypeALL_annotated_dated_collapsed_tree_',gsub('[^A-Za-z0-9]','',paste(control, collapse='')),'.newick') )
	write.tree(ph, file=outfile)	
}

## ---- rmd.chunk.seattle.190723.phydy.kchsx ----
seattle.191017.phydy.kchsx <- function()
{
	require(data.table)
	require(ape)
	
	home <- '~/Box Sync/OR_Work/Seattle/analysis_191017'
	indir <- file.path(home,'phydy')	
	outdir <- file.path(home,'phydy')
	infiles <- data.table(F=list.files(indir, pattern='^SubtypeALL_', full.names=TRUE, recursive=TRUE))
	infiles[, OPT:= gsub('^.*_([A-Za-z0-9]+).newick$','\\1',basename(F))]
	localstate.regex <- '^KCHSX.*'
	
	#	select KCHSX analysis
	localstate <- gsub('[^A-Za-z]','',localstate.regex)
	infile <- subset(infiles, grepl(localstate,OPT))[,F]
	ph <- read.tree(infile)
	
	#	process all subtypes
	df <- data.table(TAXA=ph$tip.label)
	df[, LOCAL_TIP:= grepl(localstate.regex, TAXA)]
	df[, table(LOCAL_TIP)]
	
	#	plot phylogeny
	tip.col <- rep('blue',Ntip(ph))	
	tip.col[ grepl(localstate.regex, ph$tip.label) ] <- 'red'
	outfile <- gsub('\\.newick','.pdf',infile)
	pdf(file=outfile, w=6, h=12)
	plot(ph, show.tip.label=TRUE, tip.col=tip.col, cex=0.25)
	axisPhylo()
	dev.off()
}

## ---- rmd.chunk.seattle.191017.phydyn.volz.acuteHIV.mle.post ----
seattle.191017.phydyn.volz.msmUK.mle.post <- function()
{
	require(data.table)
	require(ggplot2)	
	
	home <- '/Users/or105/Box/OR_Work/Seattle'
	simdir <- file.path(home,'phydyn_vignettes','volz_acuteHIV_mle_sim')
	simfiles <- list.files(simdir, pattern='parms.csv', full.names=TRUE)
	outfile.base<- file.path(home,'phydyn_vignettes','volz_acuteHIV_')
	#	read true parameters
	dp <- lapply(seq_along(simfiles), function(i)
			{
				dp <- read.csv(simfiles[i], stringsAsFactors=FALSE)
				dp$F <- basename(simfiles[i])
				dp
			})
	dp <- as.data.table(do.call('rbind',dp))
	dp[, F_ID:= gsub('-parms\\.csv','',F)]
	setnames(dp, c('S0','tfin','wacute','wrisk2'), c('S0_true','tfin_true','wacute_true','wrisk2_true'))
	dp[, F:= NULL]
	
	#	read results
	infiles <- data.table(FR= list.files(file.path(home,'phydyn_vignettes'), , pattern='rds$', full.names=TRUE, recursive=TRUE))
	infiles[, LKL_APPROX:= gsub('^(.*)_([A-Z0-9]+)\\.rds','\\2',basename(FR))]
	infiles[, F_ID:= gsub('^(.*)_([A-Z0-9]+)\\.rds','\\1',basename(FR))]	
	dr <- infiles[,{
				z<- readRDS(FR)
				z<- c(z$par, lkl=z$value, fcounts=z$counts['function'], convergence=z$convergence)
				as.list(z)
			},by=c('LKL_APPROX','F_ID')]
	
	#	evaluate
	dr <- merge(dr, dp, by='F_ID')
	ggplot(dr, aes(x=wacute_true, y=exp(wacute), colour=LKL_APPROX, shape=as.factor(convergence))) + 
			geom_point() +
			geom_abline(slope=1, intercept=0) +
			coord_cartesian(ylim=c(-2,10))
	ggsave(file=paste0(outfile.base,'_wacute_accuracy.pdf'), w=8, h=6)		
	ggplot(dr, aes(x=wrisk2_true, y=exp(wrisk2), colour=LKL_APPROX, shape=as.factor(convergence))) + 
			geom_point() +
			geom_abline(slope=1, intercept=0) +
			coord_cartesian(ylim=c(-2,10))
	ggsave(file=paste0(outfile.base,'_wrisk2_accuracy.pdf'), w=8, h=6)		
	
	
}

## ---- rmd.chunk.seattle.191017.phydyn.olli.SITmf01.sim ----
seattle.191017.phydyn.olli.SITmf01.sim <- function()
{
	require(methods)
	require(inline)
	require(phydynR)
	require(data.table)
	require(ggplot2)
	
	home <- '/Users/or105/Box/OR_Work/Seattle'
	home <- '/rds/general/project/ratmann_seattle_data_analysis/live'
	simdir <- file.path(home,'phydyn_olli','olli_SITmf01_sim')		
	
	#	check for files already processed
	simfiles <- data.table(FSIM=list.files(simdir, pattern='rda$', full.name=TRUE))
	simfiles <- subset(simfiles, grepl('tree',FSIM))
	simfiles[, SIM:= as.integer(gsub('^sim([0-9]+)_[a-z]+_sample([0-9]+)_([0-9]+)\\.rda$', '\\1', basename(FSIM)))]
	simfiles[, SAMPLE:= as.integer(gsub('^sim([0-9]+)_[a-z]+_sample([0-9]+)_([0-9]+)\\.rda$', '\\2', basename(FSIM)))]
	simfiles[, REP:= as.integer(gsub('^sim([0-9]+)_[a-z]+_sample([0-9]+)_([0-9]+)\\.rda$', '\\3', basename(FSIM)))]
	
	
	#		
	#	setup model equations
	demes <- c('If0','If1','Im0','Im1')
	m <- length(demes)
	bir <- matrix( '0.', nrow=m, ncol=m, dimnames=list(demes,demes))
	bir['Im0', 'If0'] <- 'beta*beta00*Sf0*Im0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['Im0', 'If1'] <- 'beta*beta01*Sf1*Im0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['Im1', 'If0'] <- 'beta*beta10*Sf0*Im1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'	
	bir['Im1', 'If1'] <- 'beta*beta11*Sf1*Im1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If0', 'Im0'] <- 'beta*beta00*Sm0*If0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If0', 'Im1'] <- 'beta*beta01*Sm1*If0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If1', 'Im0'] <- 'beta*beta10*Sm0*If1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If1', 'Im1'] <- 'beta*beta11*Sm1*If1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	#	non deme dynamics
	nondemes <- c('Sf0','Sf1','Sm0','Sm1','Tf0','Tf1','Tm0','Tm1')
	mm <- length(nondemes)
	ndd <- setNames(rep('0.', mm), nondemes) 
	ndd['Sf0'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sf0/(Sf0+Sf1+Sm0+Sm1) -mu*Sf0 - beta*(beta00*Sf0*Im0+beta10*Sf0*Im1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Sf1'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sf1/(Sf0+Sf1+Sm0+Sm1) -mu*Sf1 - beta*(beta01*Sf1*Im0+beta11*Sf1*Im1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'	
	ndd['Sm0'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sm0/(Sf0+Sf1+Sm0+Sm1) -mu*Sm0 - beta*(beta00*Sm0*If0+beta10*Sm0*If1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Sm1'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sm1/(Sf0+Sf1+Sm0+Sm1) -mu*Sm1 - beta*(beta01*Sm1*If0+beta11*Sm1*If1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Tf0'] <- 'gamma*If0 - mu*Tf0'
	ndd['Tf1'] <- 'gamma*If1 - mu*Tf1'
	ndd['Tm0'] <- 'gamma*Im0 - mu*Tm0'
	ndd['Tm1'] <- 'gamma*Im1 - mu*Tm1'
	#	no changes in ethnicity / foreign-born status
	mig <- matrix('0.', nrow=m, ncol=m, dimnames=list(demes,demes))	
	#	deaths
	death <- setNames(rep('0.', m), demes)
	death['If0'] <- '(gamma+mu)*If0'
	death['If1'] <- '(gamma+mu)*If1'
	death['Im0'] <- '(gamma+mu)*Im0'
	death['Im1'] <- '(gamma+mu)*Im1'
	#
	model.par.names <- c('beta','beta00','beta01','beta10','beta11','gamma','mu')
	#	build stochastic model	
	dmd <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=model.par.names , rcpp=TRUE, sde=FALSE)
	dms <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=model.par.names , rcpp=TRUE, sde=TRUE)
	
	#
	#	set up distinct simulations that we want to tell apart
	waifms <- list()
	#	no spread between 0 and 1
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 1
	waifm['0','1'] <- 0
	waifm['1','0'] <- 0
	waifm['1','1'] <- 1			
	waifms[[1]] <- waifm
	#	homogeneous spread 
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.55
	waifm['0','1'] <- 0.45
	waifm['1','0'] <- 0.45
	waifm['1','1'] <- 0.55
	waifms[[2]] <- waifm
	#	symmetric 0->1 25%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.75
	waifm['0','1'] <- 0.25
	waifm['1','0'] <- 0.25
	waifm['1','1'] <- 0.75			
	waifms[[3]] <- waifm
	#	symmetric 0->1 15%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.85
	waifm['0','1'] <- 0.15
	waifm['1','0'] <- 0.15
	waifm['1','1'] <- 0.85	
	waifms[[4]] <- waifm
	#	asymmetric 0->1 25% 1->0 50%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.75
	waifm['0','1'] <- 0.25
	waifm['1','0'] <- 0.45
	waifm['1','1'] <- 0.55	
	waifms[[5]] <- waifm
	#	asymmetric 0->1 50% 1->0 25%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.55
	waifm['0','1'] <- 0.45
	waifm['1','0'] <- 0.25
	waifm['1','1'] <- 0.75			
	waifms[[6]] <- waifm	
	#	asymmetric 0->1 25% 1->0 15%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.75
	waifm['0','1'] <- 0.25
	waifm['1','0'] <- 0.15
	waifm['1','1'] <- 0.85	
	waifms[[7]] <- waifm
	#	asymmetric 0->1 15% 1->0 25%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.85
	waifm['0','1'] <- 0.15
	waifm['1','0'] <- 0.25
	waifm['1','1'] <- 0.75			
	waifms[[8]] <- waifm
	
	#for(kk in seq_along(waifms))
	for(kk in 1:2)
	{
		#	find equilibrium parameters
		propf <- 0.5
		pS <- 0.99; pI<- (1-pS)*0.2
		pSm <- (1-propf)*pS; pSf <- propf*pS
		pIm<- (1-propf)*pI; pIf<- propf*pI  		
		mu<- 1/40
		popN <- 1e6
		
		waifm <- waifms[[kk]]		
		beta00<- waifm['0','0']; beta01<- waifm['0','1']; beta10<-  waifm['1','0']; beta11 <-  waifm['1','1']			
		pI0 <- (beta11-beta10)/(beta00+beta11-beta01-beta10)
		pIf0 <- pI0*pIf; pIf1 <- (1-pI0)*pIf
		pIm0 <- pI0*pIm; pIm1 <- (1-pI0)*pIm
		beta <- mu*2*(1-pS)/(pS*pI) * ( beta00+beta11-beta01-beta10 )/( beta00*(beta11-beta10)+beta10*(beta00-beta01) )
		gamma <- mu*(1-pS-pI)/pI 
		pS0 <- pI0
		pSf0 <- pS0*pSf; pSf1 <- (1-pS0)*pSf
		pSm0 <- pS0*pSm; pSm1 <- (1-pS0)*pSm
		
		all.pars <- c(beta=beta, beta00=beta00,beta01=beta01,beta10=beta10,beta11=beta11,
				gamma=gamma, mu=mu, 
				popN=popN,
				Sf0_init= round(pSf0*popN), Sm0_init= round(pSm0*popN),
				Sf1_init= round(pSf1*popN), Sm1_init= round(pSm1*popN),
				If0_init= round(pIf0*popN), Im0_init= round(pIm0*popN),
				If1_init= round(pIf1*popN), Im1_init= round(pIm1*popN),
				Tf0_init= round(pI0*(1-pS-pI)*popN/2), Tm0_init= round(pI0*(1-pS-pI)*popN/2),
				Tf1_init= round((1-pI0)*(1-pS-pI)*popN/2), Tm1_init= round((1-pI0)*(1-pS-pI)*popN/2)
		)
		all.pars['N_init'] <- sum(all.pars[c('Sf0_init','Sf1_init','Sm0_init','Sm1_init','If0_init','If1_init','Im0_init','Im1_init','Tf0_init','Tf1_init','Tm0_init','Tm1_init')])
		
		
		#	simulate trajectories from stochastic model
		model.pars <- all.pars[c('beta','beta00','beta01','beta10','beta11','gamma','mu')]
		t0 <- 0; t1 <- 50; 
		x0 <- all.pars[c('Sf0_init','Sf1_init','Sm0_init','Sm1_init','If0_init','If1_init','Im0_init','Im1_init','Tf0_init','Tf1_init','Tm0_init','Tm1_init')]
		names(x0) <- gsub('_init','',names(x0))
		tfgy <- dmd(model.pars, x0, t0, t1, res = 1000, integrationMethod='adams')
		dbir <- as.data.table( t( sapply(tfgy[['births']], function(x){ as.numeric(x) }) ) )
		tmp <- as.vector(t(sapply( demes, function(x) paste(x,'->',demes))))
		setnames(dbir, colnames(dbir), tmp)			
		dbir[, time:= tfgy[['times']] ]		
		dsim <- as.data.table( tfgy[[5]] )
		#	plot trajectories
		tmp <- melt(dsim, id.vars=c('time'))
		ggplot(tmp, aes(x=time, colour=variable, y=value)) + 
				geom_line() + 
				theme_bw() +
				scale_y_log10()
		outfile <- file.path(simdir,paste0('sim',kk,'_trajectories_det.pdf'))
		ggsave(file=outfile, w=8, h=6)		
		
		#	simulate trajectories from stochastic model
		dsim <- list()
		dbir <- list()
		for(i in 1:5)
		{
			tfgy <- dms(model.pars, x0, t0, t1, res = 1000, integrationMethod='adams')
			dbir[[i]] <- as.data.table( t( sapply(tfgy[['births']], function(x){ as.numeric(x) }) ) )
			tmp <- as.vector(t(sapply( demes, function(x) paste(x,'->',demes))))
			setnames(dbir[[i]], colnames(dbir[[i]]), tmp)			
			dbir[[i]][, time:= tfgy[['times']] ]
			dbir[[i]][, RUN:= i]
			dsim[[i]] <- as.data.table( tfgy[[5]] )
			dsim[[i]][, RUN:= i]
		}
		dsim <- do.call('rbind',dsim)
		dbir <- do.call('rbind',dbir)
		dsim <- subset(dsim, RUN==1)
		dbir <- subset(dbir, RUN==1)
		
		#	plot trajectories
		tmp <- melt(dsim, id.vars=c('RUN','time'))
		ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
				geom_line() + 
				theme_bw() +
				scale_y_log10()
		outfile <- file.path(simdir,paste0('sim',kk,'_trajectories.pdf'))
		ggsave(file=outfile, w=8, h=6)
		
		
		
		#	plot births
		tmp <- dbir[, lapply(.SD, function(x) any(x!=0) ), .SDcols=setdiff(colnames(dbir),c('RUN','time'))]
		cols.nnzero <- subset(melt(tmp, id.vars=NULL, measure.vars=names(tmp)), value)[, as.character(variable)]
		tmp <- melt(dbir, id.vars=c('RUN','time'), measure.vars=cols.nnzero)
		ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
				geom_line() + 
				theme_bw() +
				scale_y_log10()
		outfile <- file.path(simdir,paste0('sim',kk,'_births.pdf'))
		ggsave(file=outfile, w=8, h=6)
		
		#	plot flows
		dprop <- melt(dbir, id.vars=c('RUN','time'))
		dprop <- dprop[, list(variable=variable, value=value/sum(value)), by=c('RUN','time')]
		dprop <- dcast.data.table(dprop, RUN+time~variable)
		tmp <- melt(dprop, id.vars=c('RUN','time'), measure.vars=cols.nnzero)
		ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
				geom_line() + 
				theme_bw() +
				scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1))
		outfile <- file.path(simdir,paste0('sim',kk,'_flows.pdf'))
		ggsave(file=outfile, w=8, h=6)
		
		
		#	plot onward transmissions
		donw <- melt(dbir, id.vars=c('RUN','time'))
		donw[, recipient:= gsub('^([A-Za-z0-9]+) -> ([A-Za-z0-9]+)','\\1',variable)]
		donw[, source:= gsub('^([A-Za-z0-9]+) -> ([A-Za-z0-9]+)','\\2',variable)]
		donw <- donw[, list(recipient=recipient, value=value/sum(value)),  by=c('RUN','time','source')]
		donw[, variable:= paste0(source,' -> ',recipient)]
		set(donw, NULL, c('source','recipient'), NULL)
		donw <- dcast.data.table(donw, RUN+time~variable)
		tmp <- melt(donw, id.vars=c('RUN','time'), measure.vars=cols.nnzero)
		ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
				geom_line() + 
				theme_bw() +
				scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1))		
		outfile <- file.path(simdir,paste0('sim',kk,'_onwardtransmissions.pdf'))
		ggsave(file=outfile, w=8, h=6)
		
		
		#	plot simulated waifm
		dwaifm <- melt(dbir, id.vars=c('RUN','time'))		
		dwaifm[, source:= gsub('^([A-Za-z0-9]+) -> ([A-Za-z0-9]+)','\\1',variable)]
		dwaifm[, recipient:= gsub('^([A-Za-z0-9]+) -> ([A-Za-z0-9]+)','\\2',variable)]
		tmp <- melt(dsim, id.vars=c('RUN','time'))		
		tmp <- subset(tmp, substr(variable,1,1)=='S')
		tmp[, gender:= gsub('S([m|f])[0-9]','\\1',variable)]
		tmp <- tmp[, list(variable=variable, value=value/sum(value)), by=c('RUN','time','gender')]
		set(tmp, NULL, 'variable', tmp[, gsub('^S','I', variable)])
		setnames(tmp, c('variable','value'), c('recipient','susceptible.prop'))
		dwaifm <- merge(dwaifm, tmp, by=c('RUN','time','recipient'))
		set(dwaifm, NULL, 'value', dwaifm[, value/susceptible.prop])		
		dwaifm <- dwaifm[, list(recipient=recipient, value=value/sum(value)),  by=c('RUN','time','source')]
		dwaifm[, variable:= paste0(source,' -> ',recipient)]
		set(dwaifm, NULL, c('source','recipient'), NULL)
		dwaifm <- dcast.data.table(dwaifm, RUN+time~variable)
		tmp <- melt(dwaifm, id.vars=c('RUN','time'), measure.vars=cols.nnzero)
		ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
				geom_line() + 
				theme_bw() +
				scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1))		
		outfile <- file.path(simdir,paste0('sim',kk,'_waifm.pdf'))
		ggsave(file=outfile, w=8, h=6)
		
		
		outfile <- file.path(simdir,paste0('sim',kk,'.rda'))
		save(dsim, dbir, dprop, donw, dwaifm, dms, dmd, all.pars, file=outfile)
		
		#
		#	simulate dated trees from deterministic model 
		#sampleNs <- c(1e2,5e2,1e3)
		sampleNs <- c(2e3,3e3)
		simR <- 100
		for(sampleN in sampleNs)
		{
			dprev <- melt(subset(dsim, time==t1), id.vars=c('RUN','time'))		
			dprev <- subset(dprev, substr(variable,1,1)=='I')
			dprev <- dprev[, list(variable=variable, value=value/sum(value)), by=c('RUN')]
			dprev <- dprev[, list(value=mean(value)), by='variable']
			state.prob <- setNames(vector('double', m), demes)
			state.prob[dprev$variable] <- dprev$value
			tmp <- subset(simfiles, SIM==kk & SAMPLE==sampleN)[, sort(REP)]
			sim.todo <- setdiff(1:simR, tmp)			
			for(i in sim.todo)
			{
				sampleTimes <- seq( t1-10, t1, length.out=sampleN)
				sampleStates <- t(rmultinom(sampleN, size = 1, prob=state.prob ))
				colnames(sampleStates) <- demes
				tree <- sim.co.tree(model.pars, dmd, x0, t0, sampleTimes, sampleStates, res=1e3)
				tree$all.pars <- all.pars								
				tree$ll <- phydynR:::colik( tree, 
						theta=model.pars, 
						dmd, 
						x0=x0, 
						t0=0, 
						res=200, 
						forgiveAgtY = 0, 
						AgtY_penalty = Inf, 
						step_size_res =10,
						likelihood='PL2',
						maxHeight= floor( tree$maxHeight-1 )
				)								
				save(tree, file=file.path(simdir,paste0('sim',kk,'_dettree_sample',sampleN,'_',i,'.rda')))				
				pdf(file=file.path(simdir,paste0('sim',kk,'_dettree_sample',sampleN,'_',i,'.pdf')), w=8, h=0.15*sampleN)
				plot.phylo(tree)
				dev.off()
				#ltt.plot(tree)					
			}				
		}	
		#
		#	simulate dated trees from stochastic model 
		#sampleNs <- c(1e2,5e2,1e3)
		sampleNs <- c()
		simR <- 100
		for(sampleN in sampleNs)
		{
			dprev <- melt(subset(dsim, time==t1), id.vars=c('RUN','time'))		
			dprev <- subset(dprev, substr(variable,1,1)=='I')
			dprev <- dprev[, list(variable=variable, value=value/sum(value)), by=c('RUN')]
			dprev <- dprev[, list(value=mean(value)), by='variable']
			state.prob <- setNames(vector('double', m), demes)
			state.prob[dprev$variable] <- dprev$value
			tmp <- subset(simfiles, SIM==kk & SAMPLE==sampleN)[, sort(REP)]
			sim.todo <- setdiff(1:simR, tmp)						
			for(i in sim.todo)
			{
				sampleTimes <- seq( t1-10, t1, length.out=sampleN)
				sampleStates <- t(rmultinom(sampleN, size = 1, prob=state.prob ))
				colnames(sampleStates) <- demes
				tree <- sim.co.tree(model.pars, dms, x0, t0, sampleTimes, sampleStates, res=1e3)
				tree$all.pars <- all.pars
				tree$ll <- phydynR:::colik( tree, 
						theta=model.pars, 
						dms, 
						x0=x0, 
						t0=0, 
						res=200, 
						forgiveAgtY = 0, 
						AgtY_penalty = Inf, 
						step_size_res =10,
						likelihood='PL2',
						maxHeight= floor( tree$maxHeight-1 )
				)
				save(tree, file=file.path(simdir,paste0('sim',kk,'_stotree_sample',sampleN,'_',i,'.rda')))				
				pdf(file=file.path(simdir,paste0('sim',kk,'_stotree_sample',sampleN,'_',i,'.pdf')), w=8, h=0.15*sampleN)
				plot.phylo(tree)
				dev.off()
				#ltt.plot(tree)					
			}				
		}	
	}
}

## ---- rmd.chunk.seattle.191017.phydyn.olli.SITmf01yrsXX.mle ----
seattle.191017.phydyn.olli.SITmf01yrsXX.mle <- function()
{
	require(methods)
	require(inline)
	require(phydynR)
	require(data.table)
	require(ggplot2)
	
	home <- '/Users/or105/Box/OR_Work/Seattle'
	home <- '/rds/general/project/ratmann_seattle_data_analysis/live'	
	simdir <- file.path(home,'phydyn_olli','olli_SITmf01yrsXX_sim')
	outdir <- file.path(home,'phydyn_olli','olli_SITmf01yrsXX_mle')
	simfiles <- data.table(FIN=list.files(simdir, pattern='rda$', full.names=TRUE))
	simfiles <- subset(simfiles, grepl('sim1|sim2',basename(FIN)) & grepl('_dettree_|_stotree_',basename(FIN)))
	tmp <- 'sim([0-9])_t1([0-9]+)_yrs([0-9]+)_([a-z]+)_sample([0-9]+)_([0-9]+).(rda|rds)'
	simfiles[, SIM_ID:= as.integer(gsub(tmp,'\\1', basename(FIN)))]
	simfiles[, SIM_T1:= as.integer(gsub(tmp,'\\2', basename(FIN)))]
	simfiles[, SIM_YRS:= as.integer(gsub(tmp,'\\3', basename(FIN)))]
	simfiles[, SIM_TREE:= gsub(tmp,'\\4', basename(FIN))]
	simfiles[, SIM_SAMPLE:= as.integer(gsub(tmp,'\\5', basename(FIN)))]
	simfiles[, SIM_REP:= as.integer(gsub(tmp,'\\6', basename(FIN)))]
	#	setup simulations
	tmp <- as.data.table(expand.grid(DUMMY=1,SIM_RES=c(200,400,600,1000), SIM_STEP_SIZE_RES=c(10,100)))
	simfiles[, DUMMY:= 1]
	simfiles <- merge(simfiles, tmp, by='DUMMY',allow.cartesian=TRUE)
	#	check what is done and remove
	tmp <- 'sim([0-9])_t1([0-9]+)_yrs([0-9]+)_([a-z]+)_sample([0-9]+)_res([0-9]+)_stepsizeres([0-9]+)_([0-9]+)_result.rds'
	outfiles <- data.table(FO=list.files(outdir, pattern='rds$', full.names=TRUE))
	outfiles[, SIM_ID:= as.integer(gsub(tmp,'\\1', basename(FO)))]
	outfiles[, SIM_T1:= as.integer(gsub(tmp,'\\2', basename(FO)))]
	outfiles[, SIM_YRS:= as.integer(gsub(tmp,'\\3', basename(FO)))]
	outfiles[, SIM_TREE:= gsub(tmp,'\\4', basename(FO))]
	outfiles[, SIM_SAMPLE:= as.integer(gsub(tmp,'\\5', basename(FO)))]
	outfiles[, SIM_RES:= as.integer(gsub(tmp,'\\6', basename(FO)))]
	outfiles[, SIM_STEP_SIZE_RES:= as.integer(gsub(tmp,'\\7', basename(FO)))]
	outfiles[, SIM_REP:= as.integer(gsub(tmp,'\\8', basename(FO)))]
	simfiles <- merge(simfiles, outfiles, by=c('SIM_ID','SIM_T1','SIM_YRS','SIM_TREE','SIM_SAMPLE','SIM_RES','SIM_STEP_SIZE_RES','SIM_REP'), all=TRUE)
	simfiles <- subset(simfiles, is.na(FO))
	
	#	read tree id and select subset of input files to process
	tree.id <- NA_integer_
	tmp <- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,6),
								input= return(substr(arg,8,nchar(arg))),
								NA)
					}))
	if(length(tmp)!=0)	
		tree.id<- as.integer(tmp[1])
	if(!is.na(tree.id) & tree.id>nrow(simfiles))
		stop('Found tree.id outside number of simfiles', tree.id)
	if(!is.na(tree.id))
		simfiles <- simfiles[tree.id,]
	
	
	#		
	#	setup model equations
	demes <- c('If0','If1','Im0','Im1')
	m <- length(demes)
	bir <- matrix( '0.', nrow=m, ncol=m, dimnames=list(demes,demes))
	bir['Im0', 'If0'] <- 'beta*beta00*Sf0*Im0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['Im0', 'If1'] <- 'beta*beta01*Sf1*Im0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['Im1', 'If0'] <- 'beta*beta10*Sf0*Im1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'	
	bir['Im1', 'If1'] <- 'beta*beta11*Sf1*Im1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If0', 'Im0'] <- 'beta*beta00*Sm0*If0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If0', 'Im1'] <- 'beta*beta01*Sm1*If0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If1', 'Im0'] <- 'beta*beta10*Sm0*If1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If1', 'Im1'] <- 'beta*beta11*Sm1*If1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	#	non deme dynamics
	nondemes <- c('Sf0','Sf1','Sm0','Sm1','Tf0','Tf1','Tm0','Tm1')
	mm <- length(nondemes)
	ndd <- setNames(rep('0.', mm), nondemes) 
	ndd['Sf0'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sf0/(Sf0+Sf1+Sm0+Sm1) -mu*Sf0 - beta*(beta00*Sf0*Im0+beta10*Sf0*Im1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Sf1'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sf1/(Sf0+Sf1+Sm0+Sm1) -mu*Sf1 - beta*(beta01*Sf1*Im0+beta11*Sf1*Im1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'	
	ndd['Sm0'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sm0/(Sf0+Sf1+Sm0+Sm1) -mu*Sm0 - beta*(beta00*Sm0*If0+beta10*Sm0*If1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Sm1'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sm1/(Sf0+Sf1+Sm0+Sm1) -mu*Sm1 - beta*(beta01*Sm1*If0+beta11*Sm1*If1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Tf0'] <- 'gamma*If0 - mu*Tf0'
	ndd['Tf1'] <- 'gamma*If1 - mu*Tf1'
	ndd['Tm0'] <- 'gamma*Im0 - mu*Tm0'
	ndd['Tm1'] <- 'gamma*Im1 - mu*Tm1'
	#	no changes in ethnicity / foreign-born status
	mig <- matrix('0.', nrow=m, ncol=m, dimnames=list(demes,demes))	
	#	deaths
	death <- setNames(rep('0.', m), demes)
	death['If0'] <- '(gamma+mu)*If0'
	death['If1'] <- '(gamma+mu)*If1'
	death['Im0'] <- '(gamma+mu)*Im0'
	death['Im1'] <- '(gamma+mu)*Im1'
	#
	model.par.names <- c('beta','beta00','beta01','beta10','beta11','gamma','mu')
	#	build deterministic model	
	dm <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=model.par.names , rcpp=TRUE, sde=FALSE)
	
	
	#	
	#	define objective function	
	obj.fun <- function(est.pars, tree, dm, fixed.pars)
	{				
		x0 <- fixed.pars[c('Sf0_init','Sf1_init','Sm0_init','Sm1_init',
						'If0_init','If1_init','Im0_init','Im1_init',
						'Tf0_init','Tf1_init','Tm0_init','Tm1_init')]
		names(x0) <- gsub('_init','',names(x0))		 
		gamma <- unname(fixed.pars['gamma'])
		mu <- unname(fixed.pars['mu'])
		t0 <- unname(fixed.pars['t0'])
		res <- unname(fixed.pars['res'])
		step_size_res <- unname(fixed.pars['step_size_res'])
		maxheight <- unname(fixed.pars['colik.maxheight'])
		r0 <- unname(exp(est.pars['log_r0']))
		beta00 <- unname(exp(fixed.pars['log_beta00']))
		beta01 <- unname(exp(est.pars['log_beta01']))
		beta10 <- unname(exp(est.pars['log_beta10']))
		beta11 <- unname(exp(est.pars['log_beta11']))
		beta <- r0*(gamma+mu)
		model.pars <- c(beta=beta, beta00=beta00, beta01=beta01, beta10=beta10, beta11=beta11, gamma=gamma, mu=mu)
		mll <- phydynR:::colik( tree, 
				theta=model.pars, 
				dm, 
				x0=x0, 
				t0=t0, 
				res=res, 
				forgiveAgtY = 0, 
				AgtY_penalty = Inf, 
				step_size_res =step_size_res,
				likelihood='PL2',
				maxHeight= maxheight
		)
		# track progress:
		print(c(mll=mll, r0=beta/(gamma+mu), beta=beta, beta01=beta01, beta10=beta10, beta11=beta11) )
		mll
	}
	
	
	for(i in seq_len(nrow(simfiles)))
	{
		#	load dated tree
		#i <- 1
		simfile <- simfiles[i,FIN]
		load(simfile)	
		all.pars.truth <- tree$all.pars			
		fixed.pars <- all.pars.truth[c('beta00','gamma','mu',
						'Sf0_init','Sf1_init','Sm0_init','Sm1_init',
						'If0_init','If1_init','Im0_init','Im1_init',
						'Tf0_init','Tf1_init','Tm0_init','Tm1_init')]
		fixed.pars['log_beta00'] <- log(fixed.pars['beta00'])
		fixed.pars['t0'] <- 0
		fixed.pars['t1'] <- simfiles[i,SIM_T1]
		fixed.pars['colik.maxheight'] <- floor( tree$maxHeight-1 )
		fixed.pars['res'] <- simfiles[i,SIM_RES]
		fixed.pars['step_size_res'] <- simfiles[i,SIM_STEP_SIZE_RES]
		#	find init values with non-zero likelihood
		for( r0 in seq(1,8,0.25))
		{
			est.pars.init <- c(log_r0=log(r0), log_beta01=log(0.5), log_beta10=log(0.5), log_beta11=log(0.5))
			ll <- obj.fun(est.pars.init, tree, dm, fixed.pars)
			print(ll)
			if(is.finite(ll))
				break
		}
		if(!is.finite(ll))
			stop('Could not find initial values')
		fit <- optim(par=est.pars.init, fn=obj.fun,  
				tree=tree, dm=dm, fixed.pars=fixed.pars, 
				control=list(fnscale=-1, reltol=1e-5, trace=6),
				method= 'Nelder-Mead'				
				)
		
		fit$tree <- tree
		fit$fixed.pars <- fixed.pars
		fit$est.pars.init <- est.pars.init
		
		outfile <- simfiles[i, paste0('sim',SIM_ID,'_t1',SIM_T1,'_yrs',SIM_YRS,'_',SIM_TREE,'_sample',SIM_SAMPLE,'_res',SIM_RES,'_stepsizeres',SIM_STEP_SIZE_RES,'_',SIM_REP,'_result.rds')]
		outfile <- file.path(outdir, outfile)
		saveRDS(fit, file=outfile)	
	}	
}

## ---- rmd.chunk.seattle.191017.phydyn.olli.SITmf01yrsXX.sim ----
seattle.191017.phydyn.olli.SITmf01yrsXX.sim <- function()
{
	require(methods)
	require(inline)
	require(phydynR)
	require(data.table)
	require(ggplot2)
	
	home <- '/Users/or105/Box/OR_Work/Seattle'
	home <- '/rds/general/project/ratmann_seattle_data_analysis/live'
	simdir <- file.path(home,'phydyn_olli','olli_SITmf01yrsXX_sim')		
	
	#	check for files already processed
	simfiles <- data.table(FSIM=list.files(simdir, pattern='rda$', full.name=TRUE))
	simfiles <- subset(simfiles, grepl('tree',FSIM))
	simfiles[, SIM:= as.integer(gsub('^sim([0-9]+)_[a-z]+_sample([0-9]+)_([0-9]+)\\.rda$', '\\1', basename(FSIM)))]
	simfiles[, SAMPLE:= as.integer(gsub('^sim([0-9]+)_[a-z]+_sample([0-9]+)_([0-9]+)\\.rda$', '\\2', basename(FSIM)))]
	simfiles[, REP:= as.integer(gsub('^sim([0-9]+)_[a-z]+_sample([0-9]+)_([0-9]+)\\.rda$', '\\3', basename(FSIM)))]
	
	
	#		
	#	setup model equations
	demes <- c('If0','If1','Im0','Im1')
	m <- length(demes)
	bir <- matrix( '0.', nrow=m, ncol=m, dimnames=list(demes,demes))
	bir['Im0', 'If0'] <- 'beta*beta00*Sf0*Im0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['Im0', 'If1'] <- 'beta*beta01*Sf1*Im0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['Im1', 'If0'] <- 'beta*beta10*Sf0*Im1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'	
	bir['Im1', 'If1'] <- 'beta*beta11*Sf1*Im1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If0', 'Im0'] <- 'beta*beta00*Sm0*If0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If0', 'Im1'] <- 'beta*beta01*Sm1*If0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If1', 'Im0'] <- 'beta*beta10*Sm0*If1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If1', 'Im1'] <- 'beta*beta11*Sm1*If1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	#	non deme dynamics
	nondemes <- c('Sf0','Sf1','Sm0','Sm1','Tf0','Tf1','Tm0','Tm1')
	mm <- length(nondemes)
	ndd <- setNames(rep('0.', mm), nondemes) 
	ndd['Sf0'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sf0/(Sf0+Sf1+Sm0+Sm1) -mu*Sf0 - beta*(beta00*Sf0*Im0+beta10*Sf0*Im1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Sf1'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sf1/(Sf0+Sf1+Sm0+Sm1) -mu*Sf1 - beta*(beta01*Sf1*Im0+beta11*Sf1*Im1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'	
	ndd['Sm0'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sm0/(Sf0+Sf1+Sm0+Sm1) -mu*Sm0 - beta*(beta00*Sm0*If0+beta10*Sm0*If1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Sm1'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sm1/(Sf0+Sf1+Sm0+Sm1) -mu*Sm1 - beta*(beta01*Sm1*If0+beta11*Sm1*If1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Tf0'] <- 'gamma*If0 - mu*Tf0'
	ndd['Tf1'] <- 'gamma*If1 - mu*Tf1'
	ndd['Tm0'] <- 'gamma*Im0 - mu*Tm0'
	ndd['Tm1'] <- 'gamma*Im1 - mu*Tm1'
	#	no changes in ethnicity / foreign-born status
	mig <- matrix('0.', nrow=m, ncol=m, dimnames=list(demes,demes))	
	#	deaths
	death <- setNames(rep('0.', m), demes)
	death['If0'] <- '(gamma+mu)*If0'
	death['If1'] <- '(gamma+mu)*If1'
	death['Im0'] <- '(gamma+mu)*Im0'
	death['Im1'] <- '(gamma+mu)*Im1'
	#
	model.par.names <- c('beta','beta00','beta01','beta10','beta11','gamma','mu')
	#	build stochastic model	
	dmd <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=model.par.names , rcpp=TRUE, sde=FALSE)
	dms <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=model.par.names , rcpp=TRUE, sde=TRUE)
	
	#
	#	set up distinct simulations that we want to tell apart
	waifms <- list()
	#	no spread between 0 and 1
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 1
	waifm['0','1'] <- 0
	waifm['1','0'] <- 0
	waifm['1','1'] <- 1			
	waifms[[1]] <- waifm
	#	homogeneous spread 
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.55
	waifm['0','1'] <- 0.45
	waifm['1','0'] <- 0.45
	waifm['1','1'] <- 0.55
	waifms[[2]] <- waifm
	
	t1s <- c(100,200,400)
	for(t1 in t1s)
	{
		for(kk in seq_along(waifms))	
		{
			t0 <- 0; 
			#t1 <- 200; 
			
			#	find equilibrium parameters
			propf <- 0.5
			pS <- 0.99; pI<- (1-pS)*0.2
			pSm <- (1-propf)*pS; pSf <- propf*pS
			pIm<- (1-propf)*pI; pIf<- propf*pI  		
			mu<- 1/40
			popN <- 1e6
			
			waifm <- waifms[[kk]]		
			beta00<- waifm['0','0']; beta01<- waifm['0','1']; beta10<-  waifm['1','0']; beta11 <-  waifm['1','1']			
			pI0 <- (beta11-beta10)/(beta00+beta11-beta01-beta10)
			pIf0 <- pI0*pIf; pIf1 <- (1-pI0)*pIf
			pIm0 <- pI0*pIm; pIm1 <- (1-pI0)*pIm
			beta <- mu*2*(1-pS)/(pS*pI) * ( beta00+beta11-beta01-beta10 )/( beta00*(beta11-beta10)+beta10*(beta00-beta01) )
			gamma <- mu*(1-pS-pI)/pI 
			pS0 <- pI0
			pSf0 <- pS0*pSf; pSf1 <- (1-pS0)*pSf
			pSm0 <- pS0*pSm; pSm1 <- (1-pS0)*pSm
			
			all.pars <- c(beta=beta, beta00=beta00,beta01=beta01,beta10=beta10,beta11=beta11,
					gamma=gamma, mu=mu, 
					popN=popN,
					t1=t1,
					Sf0_init= round(pSf0*popN), Sm0_init= round(pSm0*popN),
					Sf1_init= round(pSf1*popN), Sm1_init= round(pSm1*popN),
					If0_init= round(pIf0*popN), Im0_init= round(pIm0*popN),
					If1_init= round(pIf1*popN), Im1_init= round(pIm1*popN),
					Tf0_init= round(pI0*(1-pS-pI)*popN/2), Tm0_init= round(pI0*(1-pS-pI)*popN/2),
					Tf1_init= round((1-pI0)*(1-pS-pI)*popN/2), Tm1_init= round((1-pI0)*(1-pS-pI)*popN/2)
			)
			all.pars['N_init'] <- sum(all.pars[c('Sf0_init','Sf1_init','Sm0_init','Sm1_init','If0_init','If1_init','Im0_init','Im1_init','Tf0_init','Tf1_init','Tm0_init','Tm1_init')])
			
			
			#	simulate trajectories from stochastic model
			model.pars <- all.pars[c('beta','beta00','beta01','beta10','beta11','gamma','mu')]		
			x0 <- all.pars[c('Sf0_init','Sf1_init','Sm0_init','Sm1_init','If0_init','If1_init','Im0_init','Im1_init','Tf0_init','Tf1_init','Tm0_init','Tm1_init')]
			names(x0) <- gsub('_init','',names(x0))
			tfgy <- dmd(model.pars, x0, t0, t1, res = 1000, integrationMethod='adams')
			dbir <- as.data.table( t( sapply(tfgy[['births']], function(x){ as.numeric(x) }) ) )
			tmp <- as.vector(t(sapply( demes, function(x) paste(x,'->',demes))))
			setnames(dbir, colnames(dbir), tmp)			
			dbir[, time:= tfgy[['times']] ]		
			dsim <- as.data.table( tfgy[[5]] )
			#	plot trajectories
			tmp <- melt(dsim, id.vars=c('time'))
			ggplot(tmp, aes(x=time, colour=variable, y=value)) + 
					geom_line() + 
					theme_bw() +
					scale_y_log10()
			outfile <- file.path(simdir,paste0('sim',kk,'_t1',t1,'_trajectories_det.pdf'))
			ggsave(file=outfile, w=8, h=6)		
			
			#	simulate trajectories from stochastic model
			dsim <- list()
			dbir <- list()
			for(i in 1:5)
			{
				tfgy <- dms(model.pars, x0, t0, t1, res = 1000, integrationMethod='adams')
				dbir[[i]] <- as.data.table( t( sapply(tfgy[['births']], function(x){ as.numeric(x) }) ) )
				tmp <- as.vector(t(sapply( demes, function(x) paste(x,'->',demes))))
				setnames(dbir[[i]], colnames(dbir[[i]]), tmp)			
				dbir[[i]][, time:= tfgy[['times']] ]
				dbir[[i]][, RUN:= i]
				dsim[[i]] <- as.data.table( tfgy[[5]] )
				dsim[[i]][, RUN:= i]
			}
			dsim <- do.call('rbind',dsim)
			dbir <- do.call('rbind',dbir)
			dsim <- subset(dsim, RUN==1)
			dbir <- subset(dbir, RUN==1)
			
			#	plot trajectories
			tmp <- melt(dsim, id.vars=c('RUN','time'))
			ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
					geom_line() + 
					theme_bw() +
					scale_y_log10()
			outfile <- file.path(simdir,paste0('sim',kk,'_t1',t1,'_trajectories.pdf'))
			ggsave(file=outfile, w=8, h=6)
			
			
			
			#	plot births
			tmp <- dbir[, lapply(.SD, function(x) any(x!=0) ), .SDcols=setdiff(colnames(dbir),c('RUN','time'))]
			cols.nnzero <- subset(melt(tmp, id.vars=NULL, measure.vars=names(tmp)), value)[, as.character(variable)]
			tmp <- melt(dbir, id.vars=c('RUN','time'), measure.vars=cols.nnzero)
			ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
					geom_line() + 
					theme_bw() +
					scale_y_log10()
			outfile <- file.path(simdir,paste0('sim',kk,'_t1',t1,'_births.pdf'))
			ggsave(file=outfile, w=8, h=6)
			
			#	plot flows
			dprop <- melt(dbir, id.vars=c('RUN','time'))
			dprop <- dprop[, list(variable=variable, value=value/sum(value)), by=c('RUN','time')]
			dprop <- dcast.data.table(dprop, RUN+time~variable)
			tmp <- melt(dprop, id.vars=c('RUN','time'), measure.vars=cols.nnzero)
			ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
					geom_line() + 
					theme_bw() +
					scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1))
			outfile <- file.path(simdir,paste0('sim',kk,'_t1',t1,'_flows.pdf'))
			ggsave(file=outfile, w=8, h=6)
			
			
			#	plot onward transmissions
			donw <- melt(dbir, id.vars=c('RUN','time'))
			donw[, recipient:= gsub('^([A-Za-z0-9]+) -> ([A-Za-z0-9]+)','\\1',variable)]
			donw[, source:= gsub('^([A-Za-z0-9]+) -> ([A-Za-z0-9]+)','\\2',variable)]
			donw <- donw[, list(recipient=recipient, value=value/sum(value)),  by=c('RUN','time','source')]
			donw[, variable:= paste0(source,' -> ',recipient)]
			set(donw, NULL, c('source','recipient'), NULL)
			donw <- dcast.data.table(donw, RUN+time~variable)
			tmp <- melt(donw, id.vars=c('RUN','time'), measure.vars=cols.nnzero)
			ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
					geom_line() + 
					theme_bw() +
					scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1))		
			outfile <- file.path(simdir,paste0('sim',kk,'_t1',t1,'_onwardtransmissions.pdf'))
			ggsave(file=outfile, w=8, h=6)
			
			
			#	plot simulated waifm
			dwaifm <- melt(dbir, id.vars=c('RUN','time'))		
			dwaifm[, source:= gsub('^([A-Za-z0-9]+) -> ([A-Za-z0-9]+)','\\1',variable)]
			dwaifm[, recipient:= gsub('^([A-Za-z0-9]+) -> ([A-Za-z0-9]+)','\\2',variable)]
			tmp <- melt(dsim, id.vars=c('RUN','time'))		
			tmp <- subset(tmp, substr(variable,1,1)=='S')
			tmp[, gender:= gsub('S([m|f])[0-9]','\\1',variable)]
			tmp <- tmp[, list(variable=variable, value=value/sum(value)), by=c('RUN','time','gender')]
			set(tmp, NULL, 'variable', tmp[, gsub('^S','I', variable)])
			setnames(tmp, c('variable','value'), c('recipient','susceptible.prop'))
			dwaifm <- merge(dwaifm, tmp, by=c('RUN','time','recipient'))
			set(dwaifm, NULL, 'value', dwaifm[, value/susceptible.prop])		
			dwaifm <- dwaifm[, list(recipient=recipient, value=value/sum(value)),  by=c('RUN','time','source')]
			dwaifm[, variable:= paste0(source,' -> ',recipient)]
			set(dwaifm, NULL, c('source','recipient'), NULL)
			dwaifm <- dcast.data.table(dwaifm, RUN+time~variable)
			tmp <- melt(dwaifm, id.vars=c('RUN','time'), measure.vars=cols.nnzero)
			ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
					geom_line() + 
					theme_bw() +
					scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1))		
			outfile <- file.path(simdir,paste0('sim',kk,'_t1',t1,'_waifm.pdf'))
			ggsave(file=outfile, w=8, h=6)
			
			
			outfile <- file.path(simdir,paste0('sim',kk,'_t1',t1,'.rda'))
			save(dsim, dbir, dprop, donw, dwaifm, dms, dmd, all.pars, file=outfile)
			
			#
			#	simulate dated trees from deterministic model 
			sampleNs <- c(5e2,1e3,2e3)	
			sampleNs <- c(1e3)	
			yrss <- c(20,40,60)
			simR <- 100
			for(yrs in yrss)
			{
				for(sampleN in sampleNs)
				{
					dprev <- melt(subset(dsim, time==t1), id.vars=c('RUN','time'))		
					dprev <- subset(dprev, substr(variable,1,1)=='I')
					dprev <- dprev[, list(variable=variable, value=value/sum(value)), by=c('RUN')]
					dprev <- dprev[, list(value=mean(value)), by='variable']
					state.prob <- setNames(vector('double', m), demes)
					state.prob[dprev$variable] <- dprev$value
					tmp <- subset(simfiles, SIM==kk & SAMPLE==sampleN)[, sort(REP)]
					sim.todo <- setdiff(1:simR, tmp)			
					for(i in sim.todo)
					{
						sampleTimes <- seq( t1-yrs/100*t1, t1, length.out=sampleN)
						sampleStates <- t(rmultinom(sampleN, size = 1, prob=state.prob ))
						colnames(sampleStates) <- demes
						tree <- sim.co.tree(model.pars, dmd, x0, t0, sampleTimes, sampleStates, res=1e3)
						tree$all.pars <- all.pars								
						tree$ll <- phydynR:::colik( tree, 
								theta=model.pars, 
								dmd, 
								x0=x0, 
								t0=0, 
								res=1e3, 
								forgiveAgtY = 0, 
								AgtY_penalty = Inf, 
								step_size_res =10,
								likelihood='PL2',
								maxHeight= floor( tree$maxHeight-1 )
						)								
						save(tree, file=file.path(simdir,paste0('sim',kk,'_t1',t1,'_yrs',yrs,'_dettree_sample',sampleN,'_',i,'.rda')))				
						pdf(file=file.path(simdir,paste0('sim',kk,'_t1',t1,'_yrs',yrs,'_dettree_sample',sampleN,'_',i,'.pdf')), w=8, h=0.15*sampleN)
						plot.phylo(tree)
						dev.off()
						#ltt.plot(tree)					
					}				
				}		
			}			
		}
	}
}


## ---- rmd.chunk.seattle.191017.phydyn.olli.SITmf01.longsim ----
seattle.191017.phydyn.olli.SITmf01long.sim <- function()
{
	require(methods)
	require(inline)
	require(phydynR)
	require(data.table)
	require(ggplot2)
	
	home <- '/Users/or105/Box/OR_Work/Seattle'
	home <- '/rds/general/project/ratmann_seattle_data_analysis/live'
	simdir <- file.path(home,'phydyn_olli','olli_SITmf01long_sim')
	
	#	check for files already processed
	simfiles <- data.table(FSIM=list.files(simdir, pattern='rda$', full.name=TRUE))
	simfiles <- subset(simfiles, grepl('tree',FSIM))
	simfiles[, SIM:= as.integer(gsub('^sim([0-9]+)_[a-z]+_sample([0-9]+)_([0-9]+)\\.rda$', '\\1', basename(FSIM)))]
	simfiles[, SAMPLE:= as.integer(gsub('^sim([0-9]+)_[a-z]+_sample([0-9]+)_([0-9]+)\\.rda$', '\\2', basename(FSIM)))]
	simfiles[, REP:= as.integer(gsub('^sim([0-9]+)_[a-z]+_sample([0-9]+)_([0-9]+)\\.rda$', '\\3', basename(FSIM)))]
	
	#		
	#	setup model equations
	demes <- c('If0','If1','Im0','Im1')
	m <- length(demes)
	bir <- matrix( '0.', nrow=m, ncol=m, dimnames=list(demes,demes))
	bir['Im0', 'If0'] <- 'beta*beta00*Sf0*Im0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['Im0', 'If1'] <- 'beta*beta01*Sf1*Im0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['Im1', 'If0'] <- 'beta*beta10*Sf0*Im1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'	
	bir['Im1', 'If1'] <- 'beta*beta11*Sf1*Im1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If0', 'Im0'] <- 'beta*beta00*Sm0*If0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If0', 'Im1'] <- 'beta*beta01*Sm1*If0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If1', 'Im0'] <- 'beta*beta10*Sm0*If1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If1', 'Im1'] <- 'beta*beta11*Sm1*If1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	#	non deme dynamics
	nondemes <- c('Sf0','Sf1','Sm0','Sm1','Tf0','Tf1','Tm0','Tm1')
	mm <- length(nondemes)
	ndd <- setNames(rep('0.', mm), nondemes) 
	ndd['Sf0'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sf0/(Sf0+Sf1+Sm0+Sm1) -mu*Sf0 - beta*(beta00*Sf0*Im0+beta10*Sf0*Im1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Sf1'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sf1/(Sf0+Sf1+Sm0+Sm1) -mu*Sf1 - beta*(beta01*Sf1*Im0+beta11*Sf1*Im1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'	
	ndd['Sm0'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sm0/(Sf0+Sf1+Sm0+Sm1) -mu*Sm0 - beta*(beta00*Sm0*If0+beta10*Sm0*If1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Sm1'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sm1/(Sf0+Sf1+Sm0+Sm1) -mu*Sm1 - beta*(beta01*Sm1*If0+beta11*Sm1*If1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Tf0'] <- 'gamma*If0 - mu*Tf0'
	ndd['Tf1'] <- 'gamma*If1 - mu*Tf1'
	ndd['Tm0'] <- 'gamma*Im0 - mu*Tm0'
	ndd['Tm1'] <- 'gamma*Im1 - mu*Tm1'
	#	no changes in ethnicity / foreign-born status
	mig <- matrix('0.', nrow=m, ncol=m, dimnames=list(demes,demes))	
	#	deaths
	death <- setNames(rep('0.', m), demes)
	death['If0'] <- '(gamma+mu)*If0'
	death['If1'] <- '(gamma+mu)*If1'
	death['Im0'] <- '(gamma+mu)*Im0'
	death['Im1'] <- '(gamma+mu)*Im1'
	#
	model.par.names <- c('beta','beta00','beta01','beta10','beta11','gamma','mu')
	#	build stochastic model	
	dmd <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=model.par.names , rcpp=TRUE, sde=FALSE)	
	
	#
	#	set up distinct simulations that we want to tell apart
	waifms <- list()
	#	no spread between 0 and 1
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 1
	waifm['0','1'] <- 0
	waifm['1','0'] <- 0
	waifm['1','1'] <- 1			
	waifms[[1]] <- waifm
	#	homogeneous spread 
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.55
	waifm['0','1'] <- 0.45
	waifm['1','0'] <- 0.45
	waifm['1','1'] <- 0.55
	waifms[[2]] <- waifm
	#	symmetric 0->1 25%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.75
	waifm['0','1'] <- 0.25
	waifm['1','0'] <- 0.25
	waifm['1','1'] <- 0.75			
	waifms[[3]] <- waifm
	#	symmetric 0->1 15%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.85
	waifm['0','1'] <- 0.15
	waifm['1','0'] <- 0.15
	waifm['1','1'] <- 0.85	
	waifms[[4]] <- waifm
	#	asymmetric 0->1 25% 1->0 50%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.75
	waifm['0','1'] <- 0.25
	waifm['1','0'] <- 0.45
	waifm['1','1'] <- 0.55	
	waifms[[5]] <- waifm
	#	asymmetric 0->1 50% 1->0 25%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.55
	waifm['0','1'] <- 0.45
	waifm['1','0'] <- 0.25
	waifm['1','1'] <- 0.75			
	waifms[[6]] <- waifm	
	#	asymmetric 0->1 25% 1->0 15%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.75
	waifm['0','1'] <- 0.25
	waifm['1','0'] <- 0.15
	waifm['1','1'] <- 0.85	
	waifms[[7]] <- waifm
	#	asymmetric 0->1 15% 1->0 25%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.85
	waifm['0','1'] <- 0.15
	waifm['1','0'] <- 0.25
	waifm['1','1'] <- 0.75			
	waifms[[8]] <- waifm
	
	for(kk in seq_along(waifms))
	{
		#	find equilibrium parameters
		propf <- 0.5
		pS <- 0.99; pI<- (1-pS)*0.2
		pSm <- (1-propf)*pS; pSf <- propf*pS
		pIm<- (1-propf)*pI; pIf<- propf*pI  		
		mu<- 1/40
		popN <- 1e6
		
		waifm <- waifms[[kk]]		
		beta00<- waifm['0','0']; beta01<- waifm['0','1']; beta10<-  waifm['1','0']; beta11 <-  waifm['1','1']			
		pI0 <- (beta11-beta10)/(beta00+beta11-beta01-beta10)
		pIf0 <- pI0*pIf; pIf1 <- (1-pI0)*pIf
		pIm0 <- pI0*pIm; pIm1 <- (1-pI0)*pIm
		beta <- mu*2*(1-pS)/(pS*pI) * ( beta00+beta11-beta01-beta10 )/( beta00*(beta11-beta10)+beta10*(beta00-beta01) )
		gamma <- mu*(1-pS-pI)/pI 
		
		pS0 <- pI0
		pSf0 <- pS0*pSf; pSf1 <- (1-pS0)*pSf
		pSm0 <- pS0*pSm; pSm1 <- (1-pS0)*pSm
		
		all.pars <- c(beta=beta, beta00=beta00,beta01=beta01,beta10=beta10,beta11=beta11,
				gamma=gamma, mu=mu, 
				popN=popN,
				Sf0_init= round(propf*pS0*popN), Sm0_init= round((1-propf)*pS0*popN),
				Sf1_init= round(propf*(1-pS0)*popN), Sm1_init= round((1-propf)*(1-pS0)*popN)-1,
				If0_init= 0, Im0_init= 0,
				If1_init= 0, Im1_init= 1,
				Tf0_init= 0, Tm0_init= 0,
				Tf1_init= 0, Tm1_init= 0
				)
		all.pars['N_init'] <- sum(all.pars[c('Sf0_init','Sf1_init','Sm0_init','Sm1_init','If0_init','If1_init','Im0_init','Im1_init','Tf0_init','Tf1_init','Tm0_init','Tm1_init')])
		
		
		#	simulate trajectories from stochastic model
		model.pars <- all.pars[c('beta','beta00','beta01','beta10','beta11','gamma','mu')]
		t0 <- 0; t1 <- 2e4; 
		x0 <- all.pars[c('Sf0_init','Sf1_init','Sm0_init','Sm1_init','If0_init','If1_init','Im0_init','Im1_init','Tf0_init','Tf1_init','Tm0_init','Tm1_init')]
		names(x0) <- gsub('_init','',names(x0))
		tfgy <- dmd(model.pars, x0, t0, t1, res = 1000, integrationMethod='adams')
		dbir <- as.data.table( t( sapply(tfgy[['births']], function(x){ as.numeric(x) }) ) )
		tmp <- as.vector(t(sapply( demes, function(x) paste(x,'->',demes))))
		setnames(dbir, colnames(dbir), tmp)			
		dbir[, time:= tfgy[['times']] ]		
		dsim <- as.data.table( tfgy[[5]] )
		dsim[, RUN:= 1]
		#	plot trajectories
		tmp <- melt(dsim, id.vars=c('time'))
		ggplot(tmp, aes(x=time, colour=variable, y=value)) + 
				geom_line() + 
				theme_bw() +
				scale_y_log10()
		outfile <- file.path(simdir,paste0('sim',kk,'_trajectories_det.pdf'))
		ggsave(file=outfile, w=8, h=6)		
		
		outfile <- file.path(simdir,paste0('sim',kk,'.rda'))
		save(dsim, dbir, dmd, all.pars, file=outfile)
		
		#
		#	simulate dated trees from deterministic model 
		sampleNs <- c(1e2,5e2,1e3)
		simR <- 100
		for(sampleN in sampleNs)
		{
			dprev <- melt(subset(dsim, time==t1), id.vars=c('RUN','time'))		
			dprev <- subset(dprev, substr(variable,1,1)=='I')
			dprev <- dprev[, list(variable=variable, value=value/sum(value)), by=c('RUN')]
			dprev <- dprev[, list(value=mean(value)), by='variable']
			state.prob <- setNames(vector('double', m), demes)
			state.prob[dprev$variable] <- dprev$value		
			tmp <- subset(simfiles, SIM==kk & SAMPLE==sampleN)[, sort(REP)]
			sim.todo <- setdiff(1:simR, tmp)			
			for(i in sim.todo)
			{
				sampleTimes <- seq( t1-10, t1, length.out=sampleN)
				sampleStates <- t(rmultinom(sampleN, size = 1, prob=state.prob ))
				colnames(sampleStates) <- demes
				tree <- sim.co.tree(model.pars, dmd, x0, t0, sampleTimes, sampleStates, res=1e3)
				tree$all.pars <- all.pars								
				tree$ll <- phydynR:::colik( tree, 
						theta=model.pars, 
						dmd, 
						x0=x0, 
						t0=0, 
						res=200, 
						forgiveAgtY = 0, 
						AgtY_penalty = Inf, 
						step_size_res =10,
						likelihood='PL2',
						maxHeight= floor( tree$maxHeight-1 )
				)								
				save(tree, file=file.path(simdir,paste0('sim',kk,'_dettree_sample',sampleN,'_',i,'.rda')))				
				pdf(file=file.path(simdir,paste0('sim',kk,'_dettree_sample',sampleN,'_',i,'.pdf')), w=8, h=0.15*sampleN)
				plot.phylo(tree)
				dev.off()
				#ltt.plot(tree)					
			}				
		}		
	}
}

## ---- rmd.chunk.seattle.191017.phydyn.olli.SIT01.sim ----
seattle.191017.phydyn.olli.SIT01.sim <- function()
{
	require(methods)
	require(inline)
	require(phydynR)
	require(data.table)
	require(ggplot2)
	
	home <- '/Users/or105/Box/OR_Work/Seattle'
	home <- '/rds/general/project/ratmann_seattle_data_analysis/live'
	simdir <- file.path(home,'phydyn_olli','olli_SIT01_sim')		
	
	
	#		
	#	setup model equations
	demes <- c('I0','I1')
	m <- length(demes)
	bir <- matrix( '0.', nrow=m, ncol=m, dimnames=list(demes,demes))
	bir['I0', 'I0'] <- 'beta*beta00*S0*I0/(S0+I0+T0+S1+I1+T1)'	
	bir['I1', 'I0'] <- 'beta*beta10*S0*I1/(S0+I0+T0+S1+I1+T1)'			
	bir['I0', 'I1'] <- 'beta*beta01*S1*I0/(S0+I0+T0+S1+I1+T1)'	
	bir['I1', 'I1'] <- 'beta*beta11*S1*I1/(S0+I0+T0+S1+I1+T1)'
	#	non deme dynamics
	nondemes <- c('S0','S1','T0','T1')
	mm <- length(nondemes)
	ndd <- setNames(rep('0.', mm), nondemes) 
	ndd['S0'] <- 'mu*(S0+I0+T0+S1+I1+T1)*S0/(S0+S1) -mu*S0 - beta*(beta00*I0+beta10*I1)*S0/(S0+I0+T0+S1+I1+T1)'
	ndd['S1'] <- 'mu*(S0+I0+T0+S1+I1+T1)*S1/(S0+S1) -mu*S1 - beta*(beta01*I0+beta11*I1)*S1/(S0+I0+T0+S1+I1+T1)'	
	ndd['T0'] <- 'gamma*I0 - mu*T0'
	ndd['T1'] <- 'gamma*I1 - mu*T1'
	#	no changes in ethnicity / foreign-born status
	mig <- matrix('0.', nrow=m, ncol=m, dimnames=list(demes,demes))	
	#	deaths
	death <- setNames(rep('0.', m), demes)
	death['I0'] <- '(gamma+mu)*I0'
	death['I1'] <- '(gamma+mu)*I1'
	#
	model.par.names <- c('beta','beta00','beta01','beta10','beta11','gamma','mu')
	#	build stochastic model	
	dmd <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=model.par.names , rcpp=TRUE, sde=FALSE)
	dms <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=model.par.names , rcpp=TRUE, sde=TRUE)
	
	#
	#	set up distinct simulations that we want to tell apart
	waifms <- list()
	#	no spread between 0 and 1
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 1
	waifm['0','1'] <- 0
	waifm['1','0'] <- 0
	waifm['1','1'] <- 1			
	waifms[[1]] <- waifm
	#	homogeneous spread 
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.55
	waifm['0','1'] <- 0.45
	waifm['1','0'] <- 0.45
	waifm['1','1'] <- 0.55
	waifms[[2]] <- waifm
	#	symmetric 0->1 25%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.75
	waifm['0','1'] <- 0.25
	waifm['1','0'] <- 0.25
	waifm['1','1'] <- 0.75			
	waifms[[3]] <- waifm
	#	symmetric 0->1 15%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.85
	waifm['0','1'] <- 0.15
	waifm['1','0'] <- 0.15
	waifm['1','1'] <- 0.85	
	waifms[[4]] <- waifm
	#	asymmetric 0->1 25% 1->0 50%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.75
	waifm['0','1'] <- 0.25
	waifm['1','0'] <- 0.45
	waifm['1','1'] <- 0.55	
	waifms[[5]] <- waifm
	#	asymmetric 0->1 50% 1->0 25%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.55
	waifm['0','1'] <- 0.45
	waifm['1','0'] <- 0.25
	waifm['1','1'] <- 0.75			
	waifms[[6]] <- waifm	
	#	asymmetric 0->1 25% 1->0 15%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.75
	waifm['0','1'] <- 0.25
	waifm['1','0'] <- 0.15
	waifm['1','1'] <- 0.85	
	waifms[[7]] <- waifm
	#	asymmetric 0->1 15% 1->0 25%
	waifm <- matrix(0, 2, 2, dimnames=list(c(0,1),c(0,1)))
	waifm['0','0'] <- 0.85
	waifm['0','1'] <- 0.15
	waifm['1','0'] <- 0.25
	waifm['1','1'] <- 0.75			
	waifms[[8]] <- waifm
	
	for(kk in seq_along(waifms))
	{
		#	find equilibrium parameters
		pS <- 0.99; pI<- (1-pS)*0.2; pT<- (1-pS)*0.8; mu<- 1/40
		waifm <- waifms[[kk]]
		beta00<- waifm['0','0']; beta01<- waifm['0','1']; beta10<-  waifm['1','0']; beta11 <-  waifm['1','1']
		pI0 <- (beta11-beta10)/(beta00+beta11-beta01-beta10)
		beta <- mu*(1-pS)/(pS*pI) * ( beta00+beta11-beta01-beta10 )/( beta00*(beta11-beta10)+beta10*(beta00-beta01) )
		gamma <- mu*(1-pS-pI)/pI
		pS0 <- (beta11-beta10)/(beta00+beta11-beta01-beta10)		
		popN <- 1e6	
		all.pars <- c(beta=beta, beta00=beta00,beta01=beta01,beta10=beta10,beta11=beta11,
				gamma=gamma, mu=mu, 
				popN=popN,
				S0_init= round(pS0*pS*popN),
				S1_init= round((1-pS0)*pS*popN),
				I0_init= round(pI0*pI*popN),
				I1_init= round((1-pI0)*pI*popN),
				T0_init= round(pI0*(1-pS-pI)*popN),
				T1_init= round((1-pI0)*(1-pS-pI)*popN)
				)
		all.pars['N_init'] <- sum(all.pars[c('S0_init','S1_init','I0_init','I1_init','T0_init','T1_init')])
		
		
		#	simulate trajectories from deterministic model
		model.pars <- all.pars[c('beta','beta00','beta01','beta10','beta11','gamma','mu')]		
		t0 <- 0; t1 <- 50; 
		x0 <- all.pars[c('S0_init','S1_init','I0_init','I1_init','T0_init','T1_init')]
		names(x0) <- gsub('_init','',names(x0))
		
		tfgy <- dmd(model.pars, x0, t0, t1, res = 1000, integrationMethod='adams')
		dbir <- as.data.table( t( sapply(tfgy[['births']], function(x){ as.numeric(x) }) ) )
		tmp <- as.vector(t(sapply( demes, function(x) paste(x,'->',demes))))
		setnames(dbir, colnames(dbir), tmp)			
		dbir[, time:= tfgy[['times']] ]		
		dsim <- as.data.table( tfgy[[5]] )
			
		#	plot trajectories
		tmp <- melt(dsim, id.vars=c('time'))
		ggplot(tmp, aes(x=time, colour=variable, y=value)) + 
				geom_line() + 
				theme_bw() +
				scale_y_log10()
		outfile <- file.path(simdir,paste0('sim',kk,'_trajectories_det.pdf'))
		ggsave(file=outfile, w=8, h=6)
		
		#	simulate trajectories from stochastic model
		dsim <- list()
		dbir <- list()
		for(i in 1:5)
		{
			tfgy <- dms(model.pars, x0, t0, t1, res = 1000, integrationMethod='adams')
			dbir[[i]] <- as.data.table( t( sapply(tfgy[['births']], function(x){ as.numeric(x) }) ) )
			tmp <- as.vector(t(sapply( demes, function(x) paste(x,'->',demes))))
			setnames(dbir[[i]], colnames(dbir[[i]]), tmp)			
			dbir[[i]][, time:= tfgy[['times']] ]
			dbir[[i]][, RUN:= i]
			dsim[[i]] <- as.data.table( tfgy[[5]] )
			dsim[[i]][, RUN:= i]
		}
		dsim <- do.call('rbind',dsim)
		dbir <- do.call('rbind',dbir)
		#	plot trajectories
		tmp <- melt(dsim, id.vars=c('RUN','time'))
		ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
				geom_line() + 
				theme_bw() +
				scale_y_log10()
		outfile <- file.path(simdir,paste0('sim',kk,'_trajectories.pdf'))
		ggsave(file=outfile, w=8, h=6)
		
		#	plot births
		tmp <- dbir[, lapply(.SD, function(x) any(x!=0) ), .SDcols=setdiff(colnames(dbir),c('RUN','time'))]
		cols.nnzero <- subset(melt(tmp, id.vars=NULL, measure.vars=names(tmp)), value)[, as.character(variable)]
		tmp <- melt(dbir, id.vars=c('RUN','time'), measure.vars=cols.nnzero)
		ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
				geom_line() + 
				theme_bw() +
				scale_y_log10()
		outfile <- file.path(simdir,paste0('sim',kk,'_births.pdf'))
		ggsave(file=outfile, w=8, h=6)
		
		#	plot flows
		dprop <- melt(dbir, id.vars=c('RUN','time'))
		dprop <- dprop[, list(variable=variable, value=value/sum(value)), by=c('RUN','time')]
		dprop <- dcast.data.table(dprop, RUN+time~variable)
		tmp <- melt(dprop, id.vars=c('RUN','time'), measure.vars=cols.nnzero)
		ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
				geom_line() + 
				theme_bw() +
				scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1))
		outfile <- file.path(simdir,paste0('sim',kk,'_flows.pdf'))
		ggsave(file=outfile, w=8, h=6)
		
		
		#	plot onward transmissions
		donw <- melt(dbir, id.vars=c('RUN','time'))
		donw[, recipient:= gsub('^([A-Za-z0-9]+) -> ([A-Za-z0-9]+)','\\1',variable)]
		donw[, source:= gsub('^([A-Za-z0-9]+) -> ([A-Za-z0-9]+)','\\2',variable)]
		donw <- donw[, list(recipient=recipient, value=value/sum(value)),  by=c('RUN','time','source')]
		donw[, variable:= paste0(source,' -> ',recipient)]
		set(donw, NULL, c('source','recipient'), NULL)
		donw <- dcast.data.table(donw, RUN+time~variable)
		tmp <- melt(donw, id.vars=c('RUN','time'), measure.vars=cols.nnzero)
		ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
				geom_line() + 
				theme_bw() +
				scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1))		
		outfile <- file.path(simdir,paste0('sim',kk,'_onwardtransmissions.pdf'))
		ggsave(file=outfile, w=8, h=6)
		
		
		#	plot simulated waifm
		dwaifm <- melt(dbir, id.vars=c('RUN','time'))		
		dwaifm[, source:= gsub('^([A-Za-z0-9]+) -> ([A-Za-z0-9]+)','\\1',variable)]
		dwaifm[, recipient:= gsub('^([A-Za-z0-9]+) -> ([A-Za-z0-9]+)','\\2',variable)]
		tmp <- melt(dsim, id.vars=c('RUN','time'))		
		tmp <- subset(tmp, substr(variable,1,1)=='S')
		tmp[, gender:= gsub('S([m|f])[0-9]','\\1',variable)]
		tmp <- tmp[, list(variable=variable, value=value/sum(value)), by=c('RUN','time','gender')]
		set(tmp, NULL, 'variable', tmp[, gsub('^S','I', variable)])
		setnames(tmp, c('variable','value'), c('recipient','susceptible.prop'))
		dwaifm <- merge(dwaifm, tmp, by=c('RUN','time','recipient'))
		set(dwaifm, NULL, 'value', dwaifm[, value/susceptible.prop])		
		dwaifm <- dwaifm[, list(recipient=recipient, value=value/sum(value)),  by=c('RUN','time','source')]
		dwaifm[, variable:= paste0(source,' -> ',recipient)]
		set(dwaifm, NULL, c('source','recipient'), NULL)
		dwaifm <- dcast.data.table(dwaifm, RUN+time~variable)
		tmp <- melt(dwaifm, id.vars=c('RUN','time'), measure.vars=cols.nnzero)
		ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
				geom_line() + 
				theme_bw() +
				scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1))		
		outfile <- file.path(simdir,paste0('sim',kk,'_waifm.pdf'))
		ggsave(file=outfile, w=8, h=6)
		
		
		outfile <- file.path(simdir,paste0('sim',kk,'.rda'))
		save(dsim, dbir, dprop, donw, dwaifm, dms, dmd, all.pars, file=outfile)
		
		#
		#	simulate dated trees from stochastic model 
		sampleNs <- c(1e2,5e2,1e3)
		simR <- 100
		for(sampleN in sampleNs)
		{
			dprev <- melt(subset(dsim, time==t1), id.vars=c('RUN','time'))		
			dprev <- subset(dprev, substr(variable,1,1)=='I')
			dprev <- dprev[, list(variable=variable, value=value/sum(value)), by=c('RUN')]
			dprev <- dprev[, list(value=mean(value)), by='variable']
			state.prob <- setNames(vector('double', m), demes)
			state.prob[dprev$variable] <- dprev$value		
			for(i in 1:simR)
			{
				sampleTimes <- seq( t1-10, t1, length.out=sampleN)
				sampleStates <- t(rmultinom(sampleN, size = 1, prob=state.prob ))
				colnames(sampleStates) <- demes
				tree <- sim.co.tree(model.pars, dms, x0, t0, sampleTimes, sampleStates, res=1e3)
				tree$all.pars <- all.pars
				save(tree, file=file.path(simdir,paste0('sim',kk,'_tree_sample',sampleN,'_',i,'.rda')))
				
				pdf(file=file.path(simdir,paste0('sim',kk,'_tree_sample',sampleN,'_',i,'.pdf')), w=8, h=0.15*sampleN)
				plot.phylo(tree)
				dev.off()
				#ltt.plot(tree)					
			}				
		}		
	}
}

## ---- rmd.chunk.seattle.191017.phydyn.olli.SITmf01.add.true.ll ----
seattle.191017.phydyn.olli.SITmf01.add.true.ll <- function()
{
	require(methods)
	require(inline)
	require(phydynR)
	require(data.table)
	require(ggplot2)
	
	#home <- '/Users/or105/Box/OR_Work/Seattle'
	home <- '/rds/general/project/ratmann_seattle_data_analysis/live'
	simdir <- file.path(home,'phydyn_olli','olli_SITmf01_sim')
	simfiles <- data.table(FSIMTREE=list.files(simdir, pattern='rda$', full.names=TRUE))	
	simfiles <- subset(simfiles, grepl('tree',basename(FSIMTREE)))
	
	#	create true likelihood values
	#	(need to re-build model)	
	#	setup model equations
	demes <- c('If0','If1','Im0','Im1')
	m <- length(demes)
	bir <- matrix( '0.', nrow=m, ncol=m, dimnames=list(demes,demes))
	bir['Im0', 'If0'] <- 'beta*beta00*Sf0*Im0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['Im0', 'If1'] <- 'beta*beta01*Sf1*Im0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['Im1', 'If0'] <- 'beta*beta10*Sf0*Im1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'	
	bir['Im1', 'If1'] <- 'beta*beta11*Sf1*Im1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If0', 'Im0'] <- 'beta*beta00*Sm0*If0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If0', 'Im1'] <- 'beta*beta01*Sm1*If0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If1', 'Im0'] <- 'beta*beta10*Sm0*If1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If1', 'Im1'] <- 'beta*beta11*Sm1*If1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	#	non deme dynamics
	nondemes <- c('Sf0','Sf1','Sm0','Sm1','Tf0','Tf1','Tm0','Tm1')
	mm <- length(nondemes)
	ndd <- setNames(rep('0.', mm), nondemes) 
	ndd['Sf0'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sf0/(Sf0+Sf1+Sm0+Sm1) -mu*Sf0 - beta*(beta00*Sf0*Im0+beta10*Sf0*Im1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Sf1'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sf1/(Sf0+Sf1+Sm0+Sm1) -mu*Sf1 - beta*(beta01*Sf1*Im0+beta11*Sf1*Im1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'	
	ndd['Sm0'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sm0/(Sf0+Sf1+Sm0+Sm1) -mu*Sm0 - beta*(beta00*Sm0*If0+beta10*Sm0*If1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Sm1'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sm1/(Sf0+Sf1+Sm0+Sm1) -mu*Sm1 - beta*(beta01*Sm1*If0+beta11*Sm1*If1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Tf0'] <- 'gamma*If0 - mu*Tf0'
	ndd['Tf1'] <- 'gamma*If1 - mu*Tf1'
	ndd['Tm0'] <- 'gamma*Im0 - mu*Tm0'
	ndd['Tm1'] <- 'gamma*Im1 - mu*Tm1'
	#	no changes in ethnicity / foreign-born status
	mig <- matrix('0.', nrow=m, ncol=m, dimnames=list(demes,demes))	
	#	deaths
	death <- setNames(rep('0.', m), demes)
	death['If0'] <- '(gamma+mu)*If0'
	death['If1'] <- '(gamma+mu)*If1'
	death['Im0'] <- '(gamma+mu)*Im0'
	death['Im1'] <- '(gamma+mu)*Im1'
	#
	model.par.names <- c('beta','beta00','beta01','beta10','beta11','gamma','mu')
	#	build deterministic model	
	dm <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=model.par.names , rcpp=TRUE, sde=FALSE)
	
	#	calculate true ll and attach to simulated tree
	for(i in seq_len(nrow(simfiles)))
	{
		simfile <- simfiles[i,FSIMTREE]
		load(simfile)
		if(is.null(tree$ll))
		{
			x0 <- tree$all.pars[c('Sf0_init','Sf1_init','Sm0_init','Sm1_init',
							'If0_init','If1_init','Im0_init','Im1_init',
							'Tf0_init','Tf1_init','Tm0_init','Tm1_init')]
			names(x0) <- gsub('_init','',names(x0))		 
			model.pars <- tree$all.pars[c('beta','beta00','beta01','beta10','beta11','gamma','mu')]
			maxheight <- floor( tree$maxHeight-1 )
			ll <- phydynR:::colik( tree, 
					theta=model.pars, 
					dm, 
					x0=x0, 
					t0=0, 
					res=200, 
					forgiveAgtY = 0, 
					AgtY_penalty = Inf, 
					step_size_res =10,
					likelihood='PL2',
					maxHeight= maxheight
			)
			tree$ll <- ll
			save(tree, file=simfile)	
		}				
	}
}

## ---- rmd.chunk.seattle.191017.phydyn.olli.SITmf01.assess ----
seattle.191017.phydyn.olli.SITmf01.assess <- function()
{
	library(RcppArmadillo)
	require(methods)
	require(inline)
	require(phydynR)
	require(data.table)
	require(ggplot2)
	
	home <- '/Users/or105/Box/OR_Work/Seattle'
	simdir <- file.path(home,'phydyn_olli','olli_SITmf01_sim')
	simfiles <- data.table(FSIM=list.files(simdir, pattern='rda$', full.names=TRUE))
	simfiles[, SIM_ID:= as.integer(gsub('^sim([0-9]+).*.rda$','\\1',basename(FSIM)))]
	tmp <- subset(simfiles, !grepl('tree',basename(FSIM)))
	simfiles <- subset(simfiles, grepl('tree',basename(FSIM)))
	setnames(simfiles, 'FSIM', 'FSIMTREE')
	simfiles <- merge(tmp, simfiles, by='SIM_ID')
	simfiles[, SAMPLE_SIZE:= as.integer(gsub('^sim([0-9]+)_[a-z]*tree_sample([0-9]+)_([0-9]+).rda$','\\2',basename(FSIMTREE)))]
	simfiles[, REP_ID:= as.integer(gsub('^sim([0-9]+)_[a-z]*tree_sample([0-9]+)_([0-9]+).rda$','\\3',basename(FSIMTREE)))]
	simfiles[, SIM_ARGS:= gsub('^sim([0-9]+)_([a-z]*tree)_sample([0-9]+)_([0-9]+).rda$','\\2',basename(FSIMTREE))]
	indir <- file.path(home,'phydyn_olli','olli_SITmf01_mle')	
	infiles <- data.table(FIN=list.files(indir, pattern='rds$', full.names=TRUE))
	infiles[, SIM_ID:= as.integer(gsub('^sim([0-9]+)_[a-z]*tree_sample([0-9]+)_([0-9]+)_res.rds$','\\1',basename(FIN)))]
	infiles[, SIM_ARGS:= gsub('^sim([0-9]+)_([a-z]*tree)_sample([0-9]+)_([0-9]+)_res.rds$','\\2',basename(FIN))]
	infiles[, SAMPLE_SIZE:= as.integer(gsub('^sim([0-9]+)_[a-z]*tree_sample([0-9]+)_([0-9]+)_res.rds$','\\2',basename(FIN)))]
	infiles[, REP_ID:= as.integer(gsub('^sim([0-9]+)_[a-z]*tree_sample([0-9]+)_([0-9]+)_res.rds$','\\3',basename(FIN)))]
	infiles <- merge(infiles, simfiles, by=c('SIM_ID','SIM_ARGS','SAMPLE_SIZE','REP_ID'))
	#	read estimates
	df <- infiles[, {
				#infile <- infiles[1,FIN]
				#print(FIN)
				infile <- FIN
				fit <- readRDS(infile)
				ans <- fit$tree$all.pars
				names(ans) <- paste0(gsub('_','.',names(ans)),'_true')
				tmp <- fit$par
				names(tmp) <- paste0(names(tmp),'_est')
				as.list( c(ans, tmp, 
							convergence=fit$convergence, 
							ll_fit= fit$value,
							ll_sim= fit$tree$ll) )
			}, by=c('SIM_ID','SIM_ARGS','SAMPLE_SIZE','REP_ID')]
	
	#	rescale all beta values so that beta=1
	setnames(df, colnames(df), gsub('\\.init','',colnames(df)))
	df[, beta00_est:= 1]
	df[, beta_est:= exp(log_r0_est)*(gamma_true+mu_true)]
	df[, beta00_est:= beta00_true/beta_est]
	df[, beta10_est:= exp(log_beta10_est)/beta_est]
	df[, beta01_est:= exp(log_beta01_est)/beta_est]
	df[, beta11_est:= exp(log_beta11_est)/beta_est]
	df[, beta_est:= 1]
	df[, beta00_true:= beta00_true/beta_true]
	df[, beta10_true:= beta10_true/beta_true]
	df[, beta01_true:= beta01_true/beta_true]
	df[, beta11_true:= beta11_true/beta_true]
	df[, beta_true:= 1]
	
	#	keep only some columns for analysis
	df <- subset(df, select=c(SIM_ID, SIM_ARGS, SAMPLE_SIZE, REP_ID, 
					beta00_true, beta01_true, beta10_true, beta11_true,
					beta00_est, beta01_est, beta10_est, beta11_est, 
					convergence, ll_fit, ll_sim
					))
	df[, SAMPLE_SIZE:= factor(SAMPLE_SIZE)]
	df[, SIM_ID:= paste0('simulation ',SIM_ID)]
	df[, convergence:= factor(convergence, levels=c(1,0), labels=c('yes','no'))]
	df[, MAE_PARS:= (abs(beta00_true-beta00_est) + 
					abs(beta01_true-beta01_est) +
					abs(beta10_true-beta10_est) +
					abs(beta11_true-beta11_est))/4]	
	df[, MRAE_LL:= abs((ll_sim-ll_fit)/ll_sim)]
	#	how many estimations have converged?
	ggplot(df, aes(x=SAMPLE_SIZE, fill=convergence)) + geom_bar() +
			theme_bw() +
			facet_grid(~SIM_ID, scales='free_y') 
	ggsave(file= file.path(indir,'runs_convergence.pdf'), w=10, h=5)
	
	#	is the log likelihood under the MLE estimates close to the value for the simulated tree?
	ggplot(df, aes(x=SAMPLE_SIZE, fill=convergence, y=abs((ll_fit-ll_sim)/ll_sim) )) +
			geom_boxplot() +
			theme_bw() +
			facet_grid(SIM_ARGS~SIM_ID, scales='free_y')
	ggsave(file= file.path(indir,'runs_convergence_lldifference.pdf'), w=9, h=9)
	#	while most runs have not converged, there is actually
	#	no difference between the runs that have/have not converged
	#	in how close the best fit MLEs are to the true MLEs 
	
	
	ggplot(df, aes(x=MRAE_LL, y=MAE_PARS, colour=SIM_ARGS)) + geom_point() +
		theme_bw() +
		scale_y_log10() +
		scale_x_continuous(label=scales:::percent) +
		facet_grid(SAMPLE_SIZE~SIM_ID, scales='free_y') +
		labs(x='(ll_sim-ll_MLE_fit)/ll_sim', y='MAE of parameters')
	ggsave(file= file.path(indir,'runs_MAE.pdf'), w=9, h=9)
	#	for SAMPLE_SIZE=500, we obtain very good results for sim1
	#	however for sim2, the LL under the MLEs are all close to the true MLEs, but fit can be very different

	ggplot(df, aes(x=ll_fit, y=ll_sim, colour= as.factor(MAE_PARS<1))) +
			geom_abline(intercept=0, slope=1) +
			geom_point() +
			facet_wrap(~SIM_ID+SIM_ARGS+SAMPLE_SIZE, scales='free', ncol=3) 
	ggsave(file= file.path(indir,'runs_llsim_vs_llfit.pdf'), w=9, h=12)
	
	#	overall results
	df2 <- melt(df, id.vars=c('SIM_ID','SIM_ARGS','SAMPLE_SIZE','REP_ID','convergence','ll_sim','ll_fit'))
	df2[, TYPE:= gsub('^([a-z0-9]+)_([a-z0-9]+)$','\\2',variable)]
	df2[, VAR:= gsub('^([a-z0-9]+)_([a-z0-9]+)$','\\1',variable)]
	
	dft <- df2[, list(REP_ID=min(REP_ID)), by=c('SIM_ID','SIM_ARGS','SAMPLE_SIZE','TYPE')]
	dft <- merge(dft, df2, by=c('SIM_ID','SIM_ARGS','SAMPLE_SIZE','TYPE','REP_ID'))	
	dft <- subset(dft, TYPE=='true' & SAMPLE_SIZE==100)
	#dfe <- subset(df2, TYPE=='est' & convergence=='yes')
	dfe <- subset(df2, TYPE=='est')
	ggplot() + 			
			geom_jitter(data=dfe, aes(x=SAMPLE_SIZE, y=value, colour=SAMPLE_SIZE, shape=convergence), width=0.25, height=0) +
			geom_boxplot(data=dfe, aes(x=SAMPLE_SIZE, y=value, group=SAMPLE_SIZE), fill='transparent', outlier.shape = NA) +
			facet_grid(SIM_ID+SIM_ARGS~VAR, scales='free_y') +
			coord_cartesian(ylim=c(-10,10)) +
			theme_bw() +
			geom_hline(data=dft, aes(yintercept=value)) +
			labs(x='', y='parameter value\n', colour='number tips in phylo', shape='ML optim converged') +
			theme(legend.position='top')
	ggsave(file= file.path(indir,'accuracy_marginals.pdf'), w=10, h=15)
}

## ---- rmd.chunk.seattle.191017.phydyn.olli.SITmf01.longassess ----
seattle.191017.phydyn.olli.SITmf01long.assess <- function()
{
	library(RcppArmadillo)
	require(methods)
	require(inline)
	require(phydynR)
	require(data.table)
	require(ggplot2)
	
	home <- '/Users/or105/Box/OR_Work/Seattle'
	indir <- file.path(home,'phydyn_olli','olli_SITmf01long_mle')	
	infiles <- data.table(FIN=list.files(indir, pattern='rds$', full.names=TRUE))
	infiles[, SIM_ID:= as.integer(gsub('^sim([0-9]+)_[a-z]*tree_sample([0-9]+)_([0-9]+)_res.rds$','\\1',basename(FIN)))]
	infiles[, SIM_ARGS:= gsub('^sim([0-9]+)_([a-z]*tree)_sample([0-9]+)_([0-9]+)_res.rds$','\\2',basename(FIN))]
	infiles[, SAMPLE_SIZE:= as.integer(gsub('^sim([0-9]+)_[a-z]*tree_sample([0-9]+)_([0-9]+)_res.rds$','\\2',basename(FIN)))]
	infiles[, REP_ID:= as.integer(gsub('^sim([0-9]+)_[a-z]*tree_sample([0-9]+)_([0-9]+)_res.rds$','\\3',basename(FIN)))]
	
	#	read estimates
	df <- infiles[, {
				#infile <- infiles[1,FIN]
				#print(FIN)
				infile <- FIN
				fit <- readRDS(infile)
				ans <- fit$tree$all.pars
				names(ans) <- paste0(gsub('_','.',names(ans)),'_true')
				tmp <- fit$par
				names(tmp) <- paste0(names(tmp),'_est')
				as.list( c(ans, tmp, 
								convergence=fit$convergence, 
								ll_fit= fit$value,
								ll_sim= fit$tree$ll) )
			}, by=c('SIM_ID','SIM_ARGS','SAMPLE_SIZE','REP_ID')]
	
	#	rescale all beta values so that beta=1
	setnames(df, colnames(df), gsub('\\.init','',colnames(df)))
	df[, beta00_est:= 1]
	df[, beta_est:= exp(log_r0_est)*(gamma_true+mu_true)]
	df[, beta00_est:= beta00_true/beta_est]
	df[, beta10_est:= exp(log_beta10_est)/beta_est]
	df[, beta01_est:= exp(log_beta01_est)/beta_est]
	df[, beta11_est:= exp(log_beta11_est)/beta_est]
	df[, beta_est:= 1]
	df[, beta00_true:= beta00_true/beta_true]
	df[, beta10_true:= beta10_true/beta_true]
	df[, beta01_true:= beta01_true/beta_true]
	df[, beta11_true:= beta11_true/beta_true]
	df[, beta_true:= 1]
	
	#	keep only some columns for analysis
	df <- subset(df, select=c(SIM_ID, SAMPLE_SIZE, REP_ID, 
					beta00_true, beta01_true, beta10_true, beta11_true,
					beta00_est, beta01_est, beta10_est, beta11_est, 
					convergence, ll_fit, ll_sim
			))
	df[, SAMPLE_SIZE:= factor(SAMPLE_SIZE)]
	df[, SIM_ID:= paste0('simulation ',SIM_ID)]
	df[, convergence:= factor(convergence, levels=c(1,0), labels=c('yes','no'))]
	df[, MAE_PARS:= (abs(beta00_true-beta00_est) + 
						abs(beta01_true-beta01_est) +
						abs(beta10_true-beta10_est) +
						abs(beta11_true-beta11_est))/4]	
	df[, MRAE_LL:= abs((ll_sim-ll_fit)/ll_sim)]
	#	how many estimations have converged?
	ggplot(df, aes(x=SAMPLE_SIZE, fill=convergence)) + geom_bar() +
			theme_bw() +
			facet_grid(~SIM_ID, scales='free_y') 
	ggsave(file= file.path(indir,'runs_convergence.pdf'), w=10, h=5)
	
	#	is the log likelihood under the MLE estimates close to the value for the simulated tree?
	ggplot(df, aes(x=SAMPLE_SIZE, fill=convergence, y=abs((ll_fit-ll_sim)/ll_sim) )) +
			geom_boxplot() +
			theme_bw() +
			facet_grid(~SIM_ID, scales='free_y')
	ggsave(file= file.path(indir,'runs_convergence_lldifference.pdf'), w=10, h=5)
	#	while most runs have not converged, there is actually
	#	no difference between the runs that have/have not converged
	#	in how close the best fit MLEs are to the true MLEs 
	
	
	ggplot(df, aes(x=MRAE_LL, y=MAE_PARS)) + geom_point() +
			theme_bw() +
			scale_y_log10() +
			scale_x_continuous(label=scales:::percent) +
			facet_grid(SAMPLE_SIZE~SIM_ID, scales='free_y') +
			labs(x='(ll_sim-ll_MLE_fit)/ll_sim', y='MAE of parameters')
	ggsave(file= file.path(indir,'runs_MAE.pdf'), w=9, h=9)
	#	for SAMPLE_SIZE=500, we obtain very good results for sim1
	#	however for sim2, the LL under the MLEs are all close to the true MLEs, but fit can be very different
	
	ggplot(df, aes(x=ll_fit, y=ll_sim, colour= as.factor(MAE_PARS<1))) +
			geom_abline(intercept=0, slope=1) +
			geom_point() +
			facet_wrap(SAMPLE_SIZE~SIM_ID, scales='free', ncol=3) 
	ggsave(file= file.path(indir,'runs_llsim_vs_llfit.pdf'), w=9, h=9)
	
	#	overall results
	df2 <- melt(df, id.vars=c('SIM_ID','SAMPLE_SIZE','REP_ID','convergence','ll_sim','ll_fit'))
	df2[, TYPE:= gsub('^([a-z0-9]+)_([a-z0-9]+)$','\\2',variable)]
	df2[, VAR:= gsub('^([a-z0-9]+)_([a-z0-9]+)$','\\1',variable)]
	
	dft <- df2[, list(REP_ID=min(REP_ID)), by=c('SIM_ID','SAMPLE_SIZE','TYPE')]
	dft <- merge(dft, df2, by=c('SIM_ID','SAMPLE_SIZE','TYPE','REP_ID'))	
	dft <- subset(dft, TYPE=='true' & SAMPLE_SIZE==100)
	dfe <- subset(df2, TYPE=='est' & convergence=='yes')
	dfe <- subset(df2, TYPE=='est')
	ggplot() + 			
			geom_jitter(data=dfe, aes(x=SAMPLE_SIZE, y=value, colour=SAMPLE_SIZE, shape=convergence), width=0.25, height=0) +
			geom_boxplot(data=dfe, aes(x=SAMPLE_SIZE, y=value, group=SAMPLE_SIZE), fill='transparent', outlier.shape = NA) +
			facet_grid(SIM_ID~VAR, scales='free_y') +
			coord_cartesian(ylim=c(-10,10)) +
			theme_bw() +
			geom_hline(data=dft, aes(yintercept=value)) +
			labs(x='', y='parameter value\n', colour='number tips in phylo', shape='ML optim converged') +
			theme(legend.position='top')
	ggsave(file= file.path(indir,'accuracy_marginals.pdf'), w=10, h=15)
}

## ---- rmd.chunk.seattle.191017.phydyn.olli.SITmf01.mle ----
seattle.191017.phydyn.olli.SITmf01.mle <- function()
{
	require(methods)
	require(inline)
	require(phydynR)
	require(data.table)
	require(ggplot2)
	
	home <- '/Users/or105/Box/OR_Work/Seattle'
	home <- '/rds/general/project/ratmann_seattle_data_analysis/live'
	if(0)
	{
		simdir <- file.path(home,'phydyn_olli','olli_SITmf01_sim')
		outdir <- file.path(home,'phydyn_olli','olli_SITmf01_mle')
		t1 <- 50		
	}
	if(1)
	{
		simdir <- file.path(home,'phydyn_olli','olli_SITmf01long_sim')
		outdir <- file.path(home,'phydyn_olli','olli_SITmf01long_mle')
		t1 <- 2e4		
	}
	 
	simfiles <- data.table(FIN=list.files(simdir, pattern='rda$', full.names=TRUE))
	simfiles <- subset(simfiles, grepl('sim1|sim2',basename(FIN)) & grepl('_dettree_|_stotree_',basename(FIN)))
	
	#	read tree id and select subset of input files to process
	tree.id <- NA_integer_
	tmp <- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,6),
								input= return(substr(arg,8,nchar(arg))),
								NA)
					}))
	if(length(tmp)!=0)	
		tree.id<- as.integer(tmp[1])
	if(!is.na(tree.id))
		simfiles <- simfiles[tree.id,]
				
	
	#		
	#	setup model equations
	demes <- c('If0','If1','Im0','Im1')
	m <- length(demes)
	bir <- matrix( '0.', nrow=m, ncol=m, dimnames=list(demes,demes))
	bir['Im0', 'If0'] <- 'beta*beta00*Sf0*Im0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['Im0', 'If1'] <- 'beta*beta01*Sf1*Im0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['Im1', 'If0'] <- 'beta*beta10*Sf0*Im1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'	
	bir['Im1', 'If1'] <- 'beta*beta11*Sf1*Im1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If0', 'Im0'] <- 'beta*beta00*Sm0*If0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If0', 'Im1'] <- 'beta*beta01*Sm1*If0/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If1', 'Im0'] <- 'beta*beta10*Sm0*If1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	bir['If1', 'Im1'] <- 'beta*beta11*Sm1*If1/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	#	non deme dynamics
	nondemes <- c('Sf0','Sf1','Sm0','Sm1','Tf0','Tf1','Tm0','Tm1')
	mm <- length(nondemes)
	ndd <- setNames(rep('0.', mm), nondemes) 
	ndd['Sf0'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sf0/(Sf0+Sf1+Sm0+Sm1) -mu*Sf0 - beta*(beta00*Sf0*Im0+beta10*Sf0*Im1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Sf1'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sf1/(Sf0+Sf1+Sm0+Sm1) -mu*Sf1 - beta*(beta01*Sf1*Im0+beta11*Sf1*Im1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'	
	ndd['Sm0'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sm0/(Sf0+Sf1+Sm0+Sm1) -mu*Sm0 - beta*(beta00*Sm0*If0+beta10*Sm0*If1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Sm1'] <- 'mu*(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)*Sm1/(Sf0+Sf1+Sm0+Sm1) -mu*Sm1 - beta*(beta01*Sm1*If0+beta11*Sm1*If1)/(Sm0+Im0+Tm0+Sf0+If0+Tf0+Sm1+Im1+Tm1+Sf1+If1+Tf1)'
	ndd['Tf0'] <- 'gamma*If0 - mu*Tf0'
	ndd['Tf1'] <- 'gamma*If1 - mu*Tf1'
	ndd['Tm0'] <- 'gamma*Im0 - mu*Tm0'
	ndd['Tm1'] <- 'gamma*Im1 - mu*Tm1'
	#	no changes in ethnicity / foreign-born status
	mig <- matrix('0.', nrow=m, ncol=m, dimnames=list(demes,demes))	
	#	deaths
	death <- setNames(rep('0.', m), demes)
	death['If0'] <- '(gamma+mu)*If0'
	death['If1'] <- '(gamma+mu)*If1'
	death['Im0'] <- '(gamma+mu)*Im0'
	death['Im1'] <- '(gamma+mu)*Im1'
	#
	model.par.names <- c('beta','beta00','beta01','beta10','beta11','gamma','mu')
	#	build deterministic model	
	dm <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=model.par.names , rcpp=TRUE, sde=FALSE)
	
	
	#	
	#	define objective function	
	obj.fun <- function(est.pars, tree, dm, fixed.pars)
	{				
		x0 <- fixed.pars[c('Sf0_init','Sf1_init','Sm0_init','Sm1_init',
							'If0_init','If1_init','Im0_init','Im1_init',
							'Tf0_init','Tf1_init','Tm0_init','Tm1_init')]
		names(x0) <- gsub('_init','',names(x0))		 
		gamma <- unname(fixed.pars['gamma'])
		mu <- unname(fixed.pars['mu'])
		t0 <- unname(fixed.pars['t0'])
		maxheight <- unname(fixed.pars['colik.maxheight'])
		r0 <- unname(exp(est.pars['log_r0']))
		beta00 <- unname(exp(fixed.pars['log_beta00']))
		beta01 <- unname(exp(est.pars['log_beta01']))
		beta10 <- unname(exp(est.pars['log_beta10']))
		beta11 <- unname(exp(est.pars['log_beta11']))
		beta <- r0*(gamma+mu)
		model.pars <- c(beta=beta, beta00=beta00, beta01=beta01, beta10=beta10, beta11=beta11, gamma=gamma, mu=mu)
		mll <- phydynR:::colik( tree, 
				theta=model.pars, 
				dm, 
				x0=x0, 
				t0=t0, 
				res=200, 
				forgiveAgtY = 0, 
				AgtY_penalty = Inf, 
				step_size_res =10,
				likelihood='PL2',
				maxHeight= maxheight
				)
		# track progress:
		print(c(mll=mll, r0=beta/(gamma+mu), beta=beta, beta01=beta01, beta10=beta10, beta11=beta11) )
		mll
	}
	
	
	for(i in seq_len(nrow(simfiles)))
	{
		#	load dated tree
		#i <- 1
		simfile <- simfiles[i,FIN]
		load(simfile)	
		all.pars.truth <- tree$all.pars			
		fixed.pars <- all.pars.truth[c('beta00','gamma','mu',
						'Sf0_init','Sf1_init','Sm0_init','Sm1_init',
						'If0_init','If1_init','Im0_init','Im1_init',
						'Tf0_init','Tf1_init','Tm0_init','Tm1_init')]
		fixed.pars['log_beta00'] <- log(fixed.pars['beta00'])
		fixed.pars['t0'] <- 0
		fixed.pars['t1'] <- t1	
		fixed.pars['colik.maxheight'] <- floor( tree$maxHeight-1 ) 	
		#	find init values with non-zero likelihood
		for( r0 in seq(1,8,0.25))
		{
			est.pars.init <- c(log_r0=log(r0), log_beta01=log(0.5), log_beta10=log(0.5), log_beta11=log(0.5))
			ll <- obj.fun(est.pars.init, tree, dm, fixed.pars)
			print(ll)
			if(is.finite(ll))
				break
		}
		if(!is.finite(ll))
			stop('Could not find initial values')
		fit <- optim(par=est.pars.init, fn=obj.fun,  
				tree=tree, dm=dm, fixed.pars=fixed.pars, 
				control=list(fnscale=-1, reltol=1e-8, trace=6),
				method= 'Nelder-Mead'				
			)
			
		fit$tree <- tree
		fit$fixed.pars <- fixed.pars
		fit$est.pars.init <- est.pars.init
		outfile <- file.path(outdir, gsub('\\.rda','_res.rds',basename(simfile)))
		saveRDS(fit, file=outfile)	
	}
	
}

## ---- rmd.chunk.seattle.191017.phydyn.olli.SITmf.mle ----
seattle.191017.phydyn.olli.SITmf.mle <- function()
{
	require(methods)
	require(inline)
	require(phydynR)
	require(data.table)
	require(ggplot2)
	
	home <- '/Users/or105/Box/OR_Work/Seattle'
	#home <- '/rds/general/project/ratmann_seattle_data_analysis/live'
	simdir <- file.path(home,'phydyn_olli','olli_SITmf_sim')
	outdir <- file.path(home,'phydyn_olli','olli_SITmf_mle')
	simfiles <- list.files(simdir, pattern='rda$', full.names=TRUE)
	
		
	#	setup model equations
	demes <- c('Im','If')
	m <- length(demes)
	bir <- matrix( '0.', nrow=m, ncol=m, dimnames=list(demes,demes))
	bir['Im', 'If'] <- 'beta*Sf*Im/(Sm+Im+Tm+Sf+If+Tf)'
	bir['If', 'Im'] <- 'beta*Sm*If/(Sm+Im+Tm+Sf+If+Tf)'
	nondemes <- c('Sm','Sf','Tm','Tf')
	mm <- length(nondemes)
	ndd <- setNames(rep('0.', mm), nondemes) 
	ndd['Sm'] <- 'mu*(Sm+Im+Tm+Sf+If+Tf)/2 -mu*Sm - beta*Sm*If/(Sm+Im+Tm+Sf+If+Tf)'
	ndd['Sf'] <- 'mu*(Sm+Im+Tm+Sf+If+Tf)/2 -mu*Sf - beta*Sf*Im/(Sm+Im+Tm+Sf+If+Tf)'
	ndd['Tm'] <- 'gamma*Im - mu*Tm'
	ndd['Tf'] <- 'gamma*If - mu*Tf'
	mig <- matrix('0.', nrow=m, ncol=m, dimnames=list(demes,demes))	
	death <- setNames(rep('0.', m), demes)
	death['Im'] <- '(gamma+mu)*Im'
	death['If'] <- '(gamma+mu)*If'
	#	build deterministic model
	dm <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=c('beta','gamma','mu') , rcpp=TRUE, sde=FALSE)
				
	#	define objective function
	obj.fun.191127 <- function(r0, tree, dm, fixed.pars)
	{				
		# converged to R0=10 without maxHeight that excludes the polytomy at the root
		x0 <- setNames(fixed.pars[c("Sf0","Sm0","If0","Im0","Tf0","Tm0")], c('Sf','Sm','If','Im','Tf','Tm'))
		gamma <- unname(fixed.pars['gamma'])
		mu <- unname(fixed.pars['mu'])
		t0 <- unname(fixed.pars['t0'])
		beta <- r0*(gamma+mu)
		model.pars <- c(beta=beta, gamma=gamma, mu=mu)
		mll <- phydynR:::colik( tree, 
				theta=model.pars, 
				dm, 
				x0=x0, 
				t0=t0, 
				res=200, 
				forgiveAgtY = 1, 
				AgtY_penalty = 1, 
				step_size_res=10,
				likelihood='PL2' 	
		)
		# track progress:
		print(c(mll=mll, r0=beta/(gamma+mu)) )
		mll
	}
	obj.fun <- function(r0, tree, dm, fixed.pars)
	{				
		x0 <- setNames(fixed.pars[c("Sf0","Sm0","If0","Im0","Tf0","Tm0")], c('Sf','Sm','If','Im','Tf','Tm'))
		gamma <- unname(fixed.pars['gamma'])
		mu <- unname(fixed.pars['mu'])
		t0 <- unname(fixed.pars['t0'])
		maxheight <- unname(fixed.pars['colik.maxheight'])
		beta <- r0*(gamma+mu)
		model.pars <- c(beta=beta, gamma=gamma, mu=mu)
		mll <- phydynR:::colik( tree, 
				theta=model.pars, 
				dm, 
				x0=x0, 
				t0=t0, 
				res=200, 
				forgiveAgtY = 0, 
				AgtY_penalty = Inf, 
				step_size_res =10,
				likelihood='PL2',
				maxHeight= maxheight
		)
		# track progress:
		print(c(mll=mll, r0=beta/(gamma+mu)) )
		mll
	}
	
	
	for(i in seq_along(simfiles))
	{
		#	load dated tree
		simfile <- simfiles[i]
		load( simfile )
		
		#	setup parameters and objective function
		all.pars.truth <- tree$all.pars
		cat('\ntrue par to estimate: r0=', all.pars.truth['beta']/(all.pars.truth['gamma']+all.pars.truth['mu']))		
		fixed.pars<- all.pars.truth[c('gamma','mu',"Sf0","Sm0","If0","Im0","Tf0","Tm0")]
		fixed.pars['t0'] <- 0
		fixed.pars['t1'] <- 50	
		fixed.pars['colik.maxheight'] <- floor( tree$maxHeight-1 ) 	
				
		#	optimise
		r0.init <- 1
		fit <- optim(par=r0.init, fn=obj.fun,  
				tree=tree, dm=dm, fixed.pars=fixed.pars, 
				control=list(fnscale=-1, abstol=1e-4, trace=6),
				method= 'Brent',
				hessian= FALSE,
				lower=0,
				upper=10
		)
		
		outfile <- file.path(outdir, gsub('\\.rda','_res.rds',basename(simfile)))
		saveRDS(fit, file=outfile)	
	}
		
	df <- data.table(FO=list.files(outdir, pattern='rds$', full.names=TRUE))
	df[, SIM:= gsub('_res\\.rds','',basename(FO))]
	df <- df[, {
				z<- readRDS(FO)
				list(r0_est=z$par, ll= z$value, convergence=z$convergence)
			}, by='SIM']
	tmp <- data.table(FI=simfiles)
	tmp[, SIM:= gsub('\\.rda','',basename(FI))]
	tmp <- tmp[, {
				load(FI)
				as.list(tree$all.pars)
			},by='SIM']
	df <- merge(df, tmp, by='SIM')
	df[, r0_true:= beta/(gamma+mu)]
	
	ggplot(df, aes(x=r0_true, y=r0_est, shape=as.factor(convergence))) + 
			geom_point() +
			geom_abline(slope=1, intercept=0) 
	ggsave(file=file.path(outdir,'accuracy_r0.pdf'), w=8, h=6)		
	
	#	questions:
	#	1. tree sim: problems associated with tree with many singletons, what to do?
	#	2. estimate start time of process
	#	3. estimate initial conditions
}

## ---- rmd.chunk.seattle.191017.phydyn.olli.SITmf.sim ----
seattle.191017.phydyn.olli.SITmf.sim <- function()
{
	require(methods)
	require(inline)
	require(phydynR)
	require(data.table)
	require(ggplot2)
	
	home <- '/Users/or105/Box/OR_Work/Seattle'
	home <- '/rds/general/project/ratmann_seattle_data_analysis/live'
	simdir <- file.path(home,'phydyn_olli','olli_SITmf_sim')		
	
	#	
	#	setup model equations
	demes <- c('Im','If')
	m <- length(demes)
	bir <- matrix( '0.', nrow=m, ncol=m, dimnames=list(demes,demes))
	bir['Im', 'If'] <- 'beta*Sf*Im/(Sm+Im+Tm+Sf+If+Tf)'
	bir['If', 'Im'] <- 'beta*Sm*If/(Sm+Im+Tm+Sf+If+Tf)'
	nondemes <- c('Sm','Sf','Tm','Tf')
	mm <- length(nondemes)
	ndd <- setNames(rep('0.', mm), nondemes) 
	ndd['Sm'] <- 'mu*(Sm+Im+Tm+Sf+If+Tf)/2 -mu*Sm - beta*Sm*If/(Sm+Im+Tm+Sf+If+Tf)'
	ndd['Sf'] <- 'mu*(Sm+Im+Tm+Sf+If+Tf)/2 -mu*Sf - beta*Sf*Im/(Sm+Im+Tm+Sf+If+Tf)'
	ndd['Tm'] <- 'gamma*Im - mu*Tm'
	ndd['Tf'] <- 'gamma*If - mu*Tf'
	mig <- matrix('0.', nrow=m, ncol=m, dimnames=list(demes,demes))	
	death <- setNames(rep('0.', m), demes)
	death['Im'] <- '(gamma+mu)*Im'
	death['If'] <- '(gamma+mu)*If'				
	#	build stochastic model
	dm <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=c('beta','gamma','mu') , rcpp=TRUE, sde=TRUE)
	
	#
	#	simulate, change proportion susceptibles at equilibrium
	pSs <- seq(0.90,0.99,0.005)
	for( kk in seq_along(pSs) )
	{
		#	find SIT equilibrium parameters
		pS <- pSs[kk]
		pSm <- 0.5*pS; pSf <- 0.5*pS; 
		pIm <- (1-pSm-pSf)*0.2*0.5; pIf<- (1-pSm-pSf)*0.2*0.5; mu<- 1/40; 
		gamma <- mu*(1-(pSf+pIf+pSm+pIm))/pIm
		beta <- (gamma+mu)*pIm/pIf/pSm
		
		#	other parameters
		popN <- 1e5
		sampleN <- 1e2
		simR <- 1	
		t0 <- 0; t1 <- 50
		
		#	collect all parameters
		all.pars <- c(beta=beta,gamma=gamma,mu=mu, pS=pS,
						N=popN, Sf0= round(pSf*popN), Sm0= round(pSm*popN), 
						If0= round(pIf*popN), Im0= round(pIm*popN), 
						Tf0= ( popN-round(pIf*popN)-round(pIm*popN)-round(pSf*popN)-round(pSm*popN) )/2,
						Tm0= ( popN-round(pIf*popN)-round(pIm*popN)-round(pSf*popN)-round(pSm*popN) )/2,
						t0=t0, t1=t1,
						sampleN=sampleN, simR=simR
						)
				
		#	simulate trajectories from stochastic model
		model.pars <- all.pars[c('beta','gamma','mu')]					
		x0 <- setNames(all.pars[c('Sm0','Sf0','Im0','If0','Tm0','Tf0')], c('Sm','Sf','Im','If','Tm','Tf'))
		dsim <- list()
		for(i in 1:5)
		{
			tfgy <- dm(model.pars, x0, t0, t1, res = 1000, integrationMethod='adams')
			dsim[[i]] <- as.data.table( tfgy[[5]] )
			dsim[[i]][, RUN:= i]
		}
		dsim <- do.call('rbind',dsim)
		
		#	plot
		tmp <- melt(dsim, id.vars=c('RUN','time'))
		ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
				geom_line() + 
				theme_bw() +
				scale_y_log10()
		ggsave(file=file.path(simdir,paste0('sim_traj_',kk,'.pdf')), w=8, h=6)
		
		#	simulate dated trees from model 
		for(i in 1:simR)
		{
			sampleTimes <- seq( 15, 25, length.out=sampleN)
			sampleStates <- t(rmultinom(sampleN, size = 1, prob = c(.5, 0.5) ))
			colnames(sampleStates) <- demes
			tree <- sim.co.tree(model.pars, dm, x0, t0, sampleTimes, sampleStates, res=1e3)
			tree$all.pars <- all.pars
			
			pdf(file=file.path(simdir,paste0('sim_tree_',kk,'_',i,'.pdf')), w=8,h=8)
			plot.phylo(tree)
			dev.off()
			
			save(tree, file=file.path(simdir,paste0('sim_tree_',kk,'_',i,'.rda')))		
		}	
	}
}

## ---- rmd.chunk.seattle.191017.phydyn.olli.SIT.sim ----
seattle.191017.phydyn.olli.SIT.sim <- function()
{
	require(methods)
	require(inline)
	require(phydynR)
	require(data.table)
	require(ggplot2)
	
	home <- '/Users/or105/Box/OR_Work/Seattle'
	#home <- '/rds/general/project/ratmann_seattle_data_analysis/live'
	simdir <- file.path(home,'phydyn_olli','olli_SIT_sim')		
	
	#	find SIT equilibrium parameters
	pS <- 0.99; pI<- (1-pS)*0.2; mu<- 1/40; 
	gamma <- mu*(1-(pS+pI))/pI
	beta <- (gamma+mu)/pS
	r0 <- beta/(gamma+mu)
	
	popN <- 1e5	
	all.pars <- c(beta=beta,gamma=gamma,R0=r0,mu=mu, N=popN, S0= round(pS*popN), I0= round(pI*popN), T0=popN-round(pI*popN)-round(pS*popN) )
	
	#	setup model equations
	demes <- c('I')
	m <- length(demes)
	bir <- matrix( '0.', nrow=m, ncol=m, dimnames=list(demes,demes))
	bir['I', 'I'] <- 'beta*S*I/(S+I+T)'	
	nondemes <- c('S','T')
	mm <- length(nondemes)
	ndd <- setNames(rep('0.', mm), nondemes) 
	ndd['S'] <- 'mu*(I+T) - beta*S*I/(S+I+T)'
	ndd['T'] <- 'gamma*I - mu*T'	
	mig <- matrix('0.', nrow=m, ncol=m, dimnames=list(demes,demes))	
	death <- setNames(rep('0.', m), demes)
	death['I'] <- '(gamma+mu)*I' 
	model.pars <- all.pars[c('beta','gamma','mu')]
	
	#	build stochastic model
	dm <- phydynR:::build.demographic.process(bir, migrations=mig, death=death, nonDeme=ndd, parameter=names(model.pars) , rcpp=TRUE, sde=TRUE)
	
	#	simulate trajectories from stochastic model
	theta <- model.pars
	t0 <- 0; t1 <- 50; 
	x0 <- setNames(all.pars[c('S0','I0','T0')], c('S','I','T'))
	dsim <- list()
	for(i in 1:5)
	{
		tfgy <- dm(theta, x0, t0, t1, res = 1000, integrationMethod='adams')
		dsim[[i]] <- as.data.table( tfgy[[5]] )
		dsim[[i]][, RUN:= i]
	}
	dsim <- do.call('rbind',dsim)
	tmp <- melt(dsim, id.vars=c('RUN','time'))
	ggplot(tmp, aes(x=time, colour=variable, y=value, linetype=as.factor(RUN))) + 
			geom_line() + 
			theme_bw() +
			scale_y_log10()
	ggsave(file=file.path(simdir,'simulation_trajectories.pdf'), w=8, h=6)	
}

## ---- rmd.chunk.seattle.191017.phydyn.volz.acuteHIV.mle ----
seattle.191017.phydyn.volz.msmUK.mle <- function()
{
	require(methods)
	require(inline)
	require(phydynR)
			
	home <- '/Users/or105/Box/OR_Work/Seattle'
	home <- '/rds/general/project/ratmann_seattle_data_analysis/live'
	#likelihood.approx <- 'QL'
	likelihood.approx <- 'PL2'	
	simdir <- file.path(home,'phydyn_vignettes','volz_acuteHIV_mle_sim')
	outdir <- file.path(home,'phydyn_vignettes',paste0('volz_acuteHIV_mle_',likelihood.approx))	
	tree.id <- 1
	
	#	read tree id
	tmp <- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,6),
								input= return(substr(arg,8,nchar(arg))),
								NA)
					}))
	if(length(tmp)!=0)	
		tree.id<- as.integer(tmp[1])
	
	## parms
	## to estimate: wrisk2, wacute, beta, init_acute, assortprob
	parms0 <- list( 
			mu = 1/40		# natural mortality rate
			, gamma0 = 1	# rate per year of progressing from 1st stage to 2nd stage of infection
			, gamma1 = 1/9	# rate per year of death in 2nd stage
			, beta = 1/4	# transmission rate
			, wrisk2 = 5	# transmission risk ratio for high risk group
			, wacute = 5	# transmission risk ratio for early stage of infection
			, init_acute = 1	# initial conditions
			, S0 = 5e3			# initial susceptible size 
			, prisk2 = .30		# expected proportion in high risk group
			, assortprob = 0.80 # proportion of transmissions reserved for within-group
			, agerate = 1/5.	# rate of 'aging' out of high risk group to low risk group
	)
	demes <- c('acute', 'chron', 'acute2', 'chron2')
	m <- 4
	
	## the model 
	#~ f = (beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2))
	#~ W = (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)
	B <- matrix( '0.', nrow=m, ncol=m)
	rownames(B) = colnames(B) <- demes
	B['acute', 'acute'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * 
			(wacute * acute / (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) * 
			(assortprob + (1-assortprob) * (1-prisk2))'
	B['acute', 'acute2'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (wacute * acute / (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) * 
			( (1-assortprob) * (prisk2))'
	B['chron', 'acute'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (chron / (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) *
			(assortprob + (1-assortprob) * (1-prisk2))'
	B['chron', 'acute2'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (chron / (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) * 
			( (1-assortprob) * (prisk2))'
	
	B['acute2', 'acute'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (wacute *wrisk2 * acute2 / (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) * 
			( (1-assortprob) * (1-prisk2))'
	B['acute2', 'acute2'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (wacute*wrisk2 * acute2 / (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) * 
			(assortprob + (1-assortprob) * prisk2)'
	B['chron2', 'acute'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (chron2 *wrisk2/ (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) *
			( (1-assortprob) * (1-prisk2))'
	B['chron2', 'acute2'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (chron2 *wrisk2/ (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) *
			(assortprob + (1-assortprob) * prisk2)'	
	M <- matrix( '0.', nrow=m, ncol=m)
	rownames(M) = colnames(M) <- demes
	M['acute', 'chron'] <- 'gamma0 * acute'
	M['acute2', 'chron2'] <- 'gamma0 * acute2'
	M['acute2','acute'] <- 'agerate * acute2'
	M['chron2', 'chron'] <- 'agerate * chron2'	
	D <- setNames(rep('0.', m), demes)
	D['acute'] <- 'mu * acute'
	D['chron'] <- 'mu * chron + gamma1*chron'
	D['acute2'] <- 'mu * acute2'
	D['chron2'] <- 'mu * chron2 + gamma1 * chron2'	
	NDD <- c( S = '-(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) + mu * S0 - mu *S' )	
	dm <- build.demographic.process ( B, migrations=M, death=D, nonDeme=NDD, parameter=names(parms0) , rcpp=TRUE)
		
	
	## initial paramaters
	x0 <- c( acute = 4, chron = 0.001, acute2= 0.001, chron2=.001, S=parms0$S0)
	theta0 <- log( c( wrisk2 = 5, wacute = 5, beta = .25, init_acute=3, assortprob = .75/(1-.75) ))
	
	#use QL approximation or PL2 approximation
	obj.fun <- function(theta, bdt=NULL, dm=NULL, likelihood.approx=NULL, parms0=NULL, x0=NULL){
		parms1 <- parms0
		parms1$wrisk2 = unname( exp( theta['wrisk2'] ))
		parms1$wacute = unname( exp( theta['wacute'] ))
		parms1$beta = unname( exp( theta['beta'] ))
		ap <- theta['assortprob'] 
		parms1$assortprob = unname( exp(ap) / (1 + exp(ap))  )
		x1 <- x0 
		x1['acute'] <- unname( exp( theta['init_acute'] ))
		print( unlist( parms1 ))
		-phydynR:::colik( bdt, theta=parms1, dm, 
				x0=x1, 
				t0=0, 
				res=200, 
				forgiveAgtY = 1, 
				AgtY_penalty = 1, 
				step_size_res=10,
				likelihood=likelihood.approx
		)
	}
	#system.time( print( obj.fun( theta0, bdt )))
	#system.time( {o <- dm( parms0, x0, 0, 250)} )
	
	# not used
	sim.dm <- function( theta )
	{
		parms1 <- parms0
		parms1$wrisk2 = unname( exp( theta['wrisk2'] ))
		parms1$wacute = unname( exp( theta['wacute'] ))
		parms1$beta = unname( exp( theta['beta'] ))
		ap <- theta['assortprob'] 
		parms1$assortprob = unname( exp(ap) / (1 + exp(ap))  )
		x1 <- x0 
		x1['acute'] <- unname( exp( theta['init_acute'] ))
		print( unlist( parms1 ))
		dm(  parms1,  x0 = x1, t0=0, t1 = 100, res =50)
	}
	#~ x <- sim.dm( theta0 )
	
	trefns <- list.files( path=simdir, pattern = '[0-9]+.nwk' , full.names=TRUE)
	treefn <- trefns[tree.id]
	#for(i in seq_len(length(trefns)))
	#{
	#	treefn <- trefns[1]
		cat('\nprocess ',treefn)
		cat('\nusing ',likelihood.approx)
		#paste0(simdir, cargs[1])
		tr <- read.tree( treefn)
		n <- length(tr$tip.label)
		ssts <- matrix(0, nrow = n, ncol = m)
		colnames(ssts) <- demes
		rownames(ssts) <- tr$tip.label
		annots <- sapply( strsplit( tr$tip.label, '_' ), '[[', 2 )
		ssts[ cbind( tr$tip.label, annots) ] <- 1
		sts <- setNames( node.depth.edgelength( tr )[1:n], tr$tip.label)
		bdt <- DatedTree( tr,sts,  sampleStates = ssts )		
		fit <- optim(par= theta0, fn=obj.fun,  control=list(abstol=1e-4, trace=6), hessian=FALSE, 
				bdt=bdt, dm=dm, likelihood.approx=likelihood.approx, parms0=parms0, x0=x0 )
		outfile <- file.path(outdir, gsub('\\.nwk',paste0('_',likelihood.approx, '.rds'),basename(treefn)))
		saveRDS(fit, file=outfile)		
	#}	
}

## ---- rmd.chunk.seattle.191017.phydyn.volz.acuteHIV.mle ----
seattle.191017.phydyn.volz.acuteHIV.mle <- function()
{
	#	install.packages("rcolgem", repos="http://R-Forge.R-project.org")
	library(ape) 
	library(akima) 	 	
	library(treedater)	
	library(rcolgem)
	library(bbmle)
	#library(phydynR)
	
	parms_truth <- list(gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, 
			beta0 = 12./10, beta1=3./100, beta2=9./100, 
			S_0=3000, I0_0=1, I1_0=0.01, I2_0=0.01, 
			m=3, mm=1)
	INFECTEDNAMES <- c('I0','I1','I2')
	births <- rbind(c('parms$beta0 * S * I0/(S+I0+I1+I2)','0','0'),
			c('parms$beta1 * S * I1/(S+I0+I1+I2)','0','0'),
			c('parms$beta2 * S * I2/(S+I0+I1+I2)','0','0')
	)
	rownames(births) <- colnames(births) <- INFECTEDNAMES
	migrations <- rbind(c('0', 'parms$gamma0 * I0', '0'),
			c('0', '0', 'parms$gamma1 * I1'),
			c('0', '0', '0')
	)
	rownames(migrations) <- colnames(migrations) <- INFECTEDNAMES
	deaths <- c('parms$mu*I0', 'parms$mu*I1', 'parms$mu*I2 + parms$gamma2 * I2')
	names(deaths) <- INFECTEDNAMES
	nonDemeDynamics <- c( S = '-parms$mu*S + parms$mu*(S + I0 + I1 + I2) - S * (parms$beta0*I0+parms$beta1*I1+parms$beta2*I2) / (S + I0 + I1 + I2)')
	
	# read the tree
	tree <- read.tree(system.file('extdata/hivSimulation.nwk', package='rcolgem'))
	# ~ the sample times are the same, because it is a homochronous sample at 50 years
	sampleTimes <- rep(50, length(tree$tip.label))
	names(sampleTimes) <- tree$tip.label
	# create a tree with dated tips and internal nodes,
	# will infer the sample states from tip labels
	bdt <- binaryDatedTree(tree, sampleTimes, sampleStatesAnnotations=INFECTEDNAMES)
	print(system.time(print(
			coalescent.log.likelihood( bdt, births,  deaths, nonDemeDynamics, 
							t0=0, 
							x0=c(I0=1, I1=.01, I2=.01, S = parms_truth$S_0), 
							migrations = migrations, 
							parms=parms_truth, 
							fgyResolution=1000, 
							integrationMethod='euler'
							)
			)))

	obj_fun <- function(lnbeta0, lnbeta1, lnbeta2, t0)
	{
		parms <- parms_truth
		parms$beta0 <- exp(lnbeta0)
		parms$beta1 <- exp(lnbeta1)
		parms$beta2 <- exp(lnbeta2)
		mll <- -coalescent.log.likelihood( bdt, births, deaths, nonDemeDynamics, 
											t0 = t0, 
											x0=c(I0=1, I1=.01, I2=.01, S = parms$S_0), 
											migrations = migrations, 
											parms=parms, 
											fgyResolution = 1000, 
											integrationMethod ='rk4')
		# track progress:
		print(c(mll, exp(c(lnbeta0, lnbeta1, lnbeta2) ), t0) )
		mll
	}
	fit <- mle2(obj_fun, start=list(lnbeta0=log(.6), lnbeta1=log(.2), lnbeta2=log(.05), t0=0), 
			method='Nelder-Mead', 
			optimizer='optim',  
			control=list(trace=6, reltol=1e-8))
	AIC(fit)
	logLik(fit)
	coef(fit)
	exp(coef(fit))
	profbeta <- profile(fit, which=1, alpha=.05, std.err=.5, trace=TRUE, tol.newmin=1 )
	c( exp( confint( profbeta ) ), TrueVal=parms_truth$beta0 )
	plot(profbeta)
	abline( v = log( parms_truth$beta0) , col='red')	
}

## ---- rmd.chunk.seattle.191017.phydyn.volz.acuteHIV.simulate.tree ----
seattle.191017.phydyn.volz.acuteHIV.simulate.tree <- function()
{
	#	install.packages("rcolgem", repos="http://R-Forge.R-project.org")
	library(ape) 
	library(akima) 	 	
	library(treedater)	
	library(rcolgem)
	#library(phydynR)
	
	#	phylody model structure
	INFECTEDNAMES <- c('I0', 'I1', 'I2')
	births <- rbind(c('beta0 * S * I0/(S+I0+I1+I2)','0','0'),
		c('beta1 * S * I1/(S+I0+I1+I2)','0','0'),
		c('beta2 * S * I2/(S+I0+I1+I2)','0','0')
		)
	rownames(births) <- colnames(births) <- INFECTEDNAMES
	migrations <- rbind(c('0', 'gamma0 * I0', '0'),
		c('0', '0', 'gamma1 * I1'),
		c('0', '0', '0')
		)
	rownames(migrations) <- colnames(migrations) <- INFECTEDNAMES
	deaths <- c('mu*I0', 'mu*I1', 'mu*I2 + gamma2 * I2')
	names(deaths) <- INFECTEDNAMES
	nonDemeDynamics <- c( S = '-mu*S + mu*(S + I0 + I1 + I2) - S * (beta0*I0+beta1*I1+beta2*I2) / (S + I0 + I1 + I2)')
	
	#	build demographic process
	demo.model <- build.demographic.process(
			births
			, nonDemeDynamics
			, migrations=migrations
			, deaths=deaths
			, parameterNames = c('beta0'
				, 'beta1'
				, 'beta2'
				, 'gamma0'
				, 'gamma1'
				, 'gamma2'
				, 'mu')
			, rcpp = TRUE
			, sde=TRUE
			)

	#	parameters			
	theta <- c( gamma0 = 1, gamma1=1/7, gamma2=1/2, mu=1/30, beta0 = 12./10, beta1=3./100, beta2=9./100 )
	t0 <- 0; t1 <- 50; x0 <- c(S = 999, I0 = 1, I1 =.1, I2 = .1)	
	show.demographic.process(demo.model, theta, x0, t0, t1)
	
	#	simulate dated tree from model 
	n <- 100
	sampleTimes <- seq( 15, 25, length.out = n)
	sampleStates <- t(rmultinom( n, size = 1, prob = c(.025, .9, .075) ))
	head(sampleStates)
	tree <- sim.co.tree(theta, demo.model, x0, t0, sampleTimes, sampleStates, res = 1e3)
	tree
	
	plot.phylo(tree)
	ltt.plot(tree)
}

## ---- rmd.chunk.seattle.191017.phydyn.nascimento.senegal ----
seattle.191017.phydyn.nascimento.senegal <- function()
{
	#install.packages("remotes")
	#remotes::install_github("thednainus/senegalHIVmodel")
	library(ape) 
	library(akima) 
	library(BayesianTools) 
	library(phydynR)
	library(treedater)
	
	
	#	demes and non-demes
	demes <- c("gpm", "gpf", "msm", "src")
	nondemes <- c()
	m <- length(demes)
	mm <- length(nondemes)
	#	basic structure of objects needed for likelihood	
	births <- matrix("0.", nrow = m, ncol = m)
	migs <- matrix("0.", nrow = m, ncol = m)
	rownames(births) = rownames(migs) = colnames(births) = colnames(migs) <- demes
	deaths <- setNames(rep("0.", m), demes)
	nonDemeDynamics <- setNames(rep("0.", mm), nondemes)
	# 	model parameters
	T0 <- 1978		
	T1 <- 2014		
	GAMMA <- 1/10
	theta.gpspline <-  function(t, parms)
		{
			if (t < T0) 
				return(parms$gpsp0)
			if (t > T1) 
				return (parms$gpsp2) 
			with(parms, aspline(x = c(T0, gpsploc, T1), y=c(gpsp0, gpsp1, gpsp2), xout = t)$y)		
		}
	theta.msmspline <- function(t, parms)
		{
			if (t < T0) 
				return(parms$msmsp0)
			if (t > T1) 
				return (parms$msmsp2) 
			with(parms, aspline(x = c(T0, msmsploc, T1), y=c(msmsp0, msmsp1, msmsp2), xout = t)$y)
		}	
	THETA <- list(	gpsp0 = 6/10, # par incidence spline GP
			gpsp1 = 4/10, 	# par incidence spline GP
			gpsp2 = 1/10,	# par incidence spline GP
			gpsploc = 1987,
			msmsp0 = 4/10,	# par incidence spline MSM
			msmsp1 = 4/10,	# par incidence spline MSM
			msmsp2 = 2/10,	# par incidence spline MSM
			msmsploc = 1995,
			maleX = 2.0,
			import = 1/20,
			srcNe = 1/10,
			gpspline = theta.gpspline,
			msmspline = theta.msmspline,
			pmsm2msm = 0.85,
			pgpf2gpm = 0.85,
			initmsm = 1,
			initgp = 1)
	SRCSIZE <<- 1e5
	X0 <- c(gpm = unname(THETA$initgp/2),
			gpf = unname(THETA$initgp/2), 
			msm = unname(THETA$initmsm), 
			src = SRCSIZE)
	#	model specfication: birth matrix
	births['msm', 'msm'] <- "parms$msmspline(t, parms) * msm * parms$pmsm2msm"
	births['msm', 'gpf'] <- "parms$msmspline(t, parms) * msm * (1-parms$pmsm2msm)"
	births['gpm', 'gpf'] <- "parms$gpspline(t, parms) * gpm * parms$maleX"
	births['gpf', 'gpm'] <- "parms$gpspline(t, parms) * gpf * parms$pgpf2gpm"
	births['gpf', 'msm'] <- "parms$gpspline(t, parms) * gpf * (1-parms$pgpf2gpm)"	
	births['src', 'src'] <- "0.5 * SRCSIZE^2 / parms$srcNe"
	#	model specfication: migration matrix
	migs['src', 'gpm'] <- "parms$import * gpm"
	migs['src', 'gpf'] <- "parms$import * gpf"
	migs['src', 'msm'] <- "parms$import * msm"
	migs['gpm', 'src'] <- "parms$import * gpm"
	migs['gpf', 'src'] <- "parms$import * gpf"
	migs['msm', 'src'] <- "parms$import * msm"
	#	model specfication: death vector
	deaths['msm'] <- "GAMMA * msm"
	deaths['gpf'] <- "GAMMA * gpf"
	deaths['gpm'] <- "GAMMA * gpm"
	deaths['src'] <- "0.5 * SRCSIZE^2 / parms$srcNe"
	
	#	build demographic model
	dm <- phydynR:::build.demographic.process(births = births, 
			deaths = deaths,
			migrations = migs, 
			parameterNames = names(THETA), 
			rcpp = FALSE,
			sde = FALSE)

	#	load data, 512 taxa in total	
	tree.all <- read.tree( system.file("data/bindTree_CGR_GTR+Gp12+3_droppedTip.tre", package = "senegalHIVmodel"))
	#	data= subtype, risk group, location, sampling year
	all.data.cgr <- read.csv(system.file("data/HIV_subtypes_summary_CGR.csv", package = "senegalHIVmodel"))
	all.data.SN <- read.csv( system.file("data/HIV_subtypes_summary_SENEGAL_noDups.csv", package = "senegalHIVmodel"))
	all_data <- senegalHIVmodel::organize_metadata(all.data.cgr, all.data.SN)
	#	read estimated sampling times (real value, from treedater)
	times <- readRDS( system.file("data/bindTree_CGR_GTR+Gp12+3_droppedTip_sts.RDS", package = "senegalHIVmodel"))
	#	make tip states matrix
	gpm <- gpf <- msm <- src <- rep(0, length(tree.all$tip.label))	
	gpm[all_data$States == "gpm"] <- 1 
	gpf[all_data$States == "gpf"] <- 1 
	msm[all_data$States == "msm"] <- 1 
	src[all_data$States == "src"] <- 1
	sampleStates <- cbind(gpm, gpf, msm, src) 
	rownames(sampleStates) <- all_data$tip.name
	#	make DatedTree object
	dated.tree <- phydynR::DatedTree(phylo = tree.all, 
			sampleTimes = times,
			sampleStates = sampleStates, 
			minEdgeLength = 2/52,
			tol = 0.1)
	
	obj_fun <- function(parameters){
		# we use unname here because "parameters" can be a vector or matrix, and 
		# sometimes it comes with column names, which I chose to remove these nam 
		# in here.
		parameters <- unname(parameters)
		# add the values of THETA to a new variable named THETA.new
		THETA.new <- THETA
		# change the values in THETA.new to the new proposals that will be evalua
		THETA.new$gpsp0 <- parameters[1] 
		THETA.new$gpsp1 <- parameters[2] 
		THETA.new$gpsp2 <- parameters[3] 
		THETA.new$gpsploc <- parameters[4] 
		THETA.new$msmsp0 <- parameters[5] 
		THETA.new$msmsp1 <- parameters[6] 
		THETA.new$msmsp2 <- parameters[7] 
		THETA.new$msmsploc <- parameters[8] 
		THETA.new$import <- parameters[9] 
		THETA.new$srcNe <- parameters[10] 
		THETA.new$pmsm2msm <- parameters[11] 
		THETA.new$pgpf2gpm <- parameters[12]
		mll <- phydynR::colik(tree = dated.tree, 
				theta = THETA.new,
				demographic.process.model = dm,
				x0 = X0,
				t0 = 1978,
				res = 1e3,
				timeOfOriginBoundaryCondition = FALSE,
				AgtY_penalty = 1,
				maxHeight = 41)
		return(mll) 
	}
	densities <- function(par)
	{
		# d1 to d3 and d5 to d7 I am using a gamma distribution with
		#    mean = R0 = 1.1 and sigma = 1
		# d4, d8 uniform distribution between the start time and
		# the final time of our simualtions (t0 = 1978, and t1 = 2014)
		# d9 exponential distribution with mean around 1/30
		# d10 exponential distribution with mean around 1/20
		d1 = dgamma(par[1], shape = 3, rate = 3/1.1, log= TRUE) #gsp0
		d2 = dgamma(par[2], shape = 3, rate = 3/1.1, log= TRUE) #gsp1
		d3 = dgamma(par[3], shape = 3, rate = 3/1.1, log= TRUE) #gsp2
		d4 = dunif(par[4], min = 1978, max = 2014, log= TRUE)	#gsploc
		d5 = dgamma(par[5], shape = 3, rate = 3/1.1, log= TRUE) #msm0
		d6 = dgamma(par[6], shape = 3, rate = 3/1.1, log= TRUE) #msm1
		d7 = dgamma(par[7], shape = 3, rate = 3/1.1, log= TRUE) #msm2
		d8 = dunif(par[8], min = 1978, max = 2014, log= TRUE)	#msmloc
		d9 = dexp(par[9], rate = 30, log = TRUE) #import
		d10 = dexp(par[10], rate = 20, log = TRUE) #srcNe
		d11 = dbeta(par[11], shape1 = 16, shape2 = 4, log = TRUE) #pmsm2msm
		d12 = dbeta(par[12], shape1 = 16, shape2 = 4, log = TRUE) #pgpf2gpm
		return(d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9 + d10 + d11 + d12)
	}
	# Create sampling, this is optional.
	#The MCMCs can automatically generate starting # conditions if a sampler is provided
	sampler <- function(n=1)
	{
		d1 = rgamma(n, shape = 3, rate = 3/1.1) #gpsp0 
		d2 = rgamma(n, shape = 3, rate = 3/1.1) #gpsp1 
		d3 = rgamma(n, shape = 3, rate = 3/1.1) #gpsp2 
		d4 = runif(n, min = 1978, max = 2014) #gpsploc 
		d5 = rgamma(n, shape = 3, rate = 3/1.1) #msmsp0 
		d6 = rgamma(n, shape = 3, rate = 3/1.1) #msmsp1 
		d7 = rgamma(n, shape = 3, rate = 3/1.1) #msmsp2 
		d8 = runif(n, min = 1978, max = 2014) #msmsploc 
		d9 = rexp(n, rate = 30) #import
		d10 = rexp(n, rate = 20) #srcNe
		d11 = rbeta(n, shape1 = 16, shape2 = 4) #pmsm2msm 
		d12 = rbeta(n, shape1 = 16, shape2 = 4) #pgpf2gpm
		return(cbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12)) 
	}
	# Create prior (necessary for the BayesianTools package)
	prior <- createPrior(
			density = densities, 
			sampler = sampler,
			lower = c(0.01, 0.01, 0.01, 1978, 0.01, 0.01, 0.01, 1978, 0, 0.0001, 0.3, 0.3), 
			upper = c(3, 3, 3,2014, 3, 3, 3, 2014, 0.15, 0.15, 1, 1)
			)
	bayesianSetup <- createBayesianSetup(likelihood=obj_fun , prior = prior)
	settings = list(iterations = 18000, nrChains = 1, thin = 1) 
	out <- runMCMC(bayesianSetup = bayesianSetup,
			sampler = "DEzs",
			settings = settings)
}

## ---- rmd.chunk.seattle.191017.treedater.run ----
seattle.191017.treedater.run <- function()
{
	require(big.phylo)
	require(data.table)
	require(tidyverse)
	require(ape)
	require(phytools)
	require(treedater)
	require(phyloscannerR)
	
	home <- '~/Box Sync/OR_Work/Seattle/analysis_191017'
	home <- '/rds/general/project/ratmann_seattle_data_analysis/live/olli_191017'
	infile.seqinfo <- file.path(dirname(home),'PHSKC-2018-07-09','sequences_meta.rds')	
	indir.phsc	<- file.path(home,'phyloscanner')
	outdir <- file.path(home,'phyloscanner_dated')
	infiles.phsc <- data.table(F=list.files(indir.phsc, pattern='workspace.rda$', full.names=TRUE, recursive=TRUE))
	#infiles.phsc <- subset(infiles.phsc, grepl('000',F))
	alignment.length <- 1000
	
	#infiles.phsc <- subset(infiles.phsc, grepl('180709_LANL_SubtypeBc11_mafft_ndrm_008_ft_rerooted_KCHSX', F))
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
	cmds <- vector('list',nrow(infiles.phsc))
	for(i in seq_len(nrow(infiles.phsc)))
	{
		#	i<- 4
		cat('\nProcess',i)
		infile <- infiles.phsc[i,F]	
		load(infile)
		ph <- phyloscanner.trees[[1]][['tree']]
		stopifnot( !any( ph$tip.label=='' ) )		
		stopifnot( is.binary(ph) )
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
		cat('\nDropping tips without sampling date from ', infile,' n=', length(tmp), 'of Ntips=', Ntip(ph), '\n')
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
		sampling.times.init <- dph$SEQ_DATE
		names(sampling.times.init) <- dph$TAXA_LABEL
		sampling.times.bounds <- as.data.frame(subset(dph, select=c(SEQ_DATE_LOWER, SEQ_DATE_UPPER)))
		rownames(sampling.times.bounds) <- dph$TAXA_LABEL		
		colnames(sampling.times.bounds) <- c('lower','upper')
		#	save files for dating
		outfile.collapsed.phsc <- file.path(outdir, gsub('workspace\\.rda','collapsed_workspace\\.rda',basename(infile)))
		outfile.tree <- file.path(outdir, gsub('workspace\\.rda','collapsed\\.newick',basename(infile)))
		outfile.sampling.times.bounds <- file.path(outdir, gsub('workspace\\.rda','sampling_times_bounds\\.csv',basename(infile)))
		outfile.sampling.times.init <- file.path(outdir, gsub('workspace\\.rda','sampling_times_init\\.csv',basename(infile)))
		save(ph, file=outfile.collapsed.phsc)		
		write.tree(ph, file=outfile.tree)
		write.csv(sampling.times.bounds, file=outfile.sampling.times.bounds, row.names=TRUE)		
		write.csv(data.table(TAXA=names(sampling.times.init), SEQ_DATE=sampling.times.init), file=outfile.sampling.times.init, row.names=FALSE)
		control <- list( outfile=gsub('.newick$','_dated.newick',outfile.tree), 
						 ali.len=alignment.length, 
						 root=NA, 
						 omega0=NA, 
						 temporalConstraints=TRUE, 
						 strictClock=FALSE, 
						 estimateSampleTimes=outfile.sampling.times.bounds 
						 )
		cmd <- cmd.treedater.script(outfile.tree, outfile.sampling.times.init, control=control)
		cmds[[i]] <- cmd
	}	
	infiles.phsc[, CMD:= unlist(cmds)]	
	
	#	submit jobs like this one:
	cat(infiles.phsc[1,CMD])
	
	#	run on HPC as array job
	df <- copy(infiles.phsc)
	df[, CASE_ID:= 1:nrow(df)]	
	#	make PBS header
	hpc.load	<- "module load anaconda3/personal"
	hpc.select	<- 1						# number of nodes
	hpc.nproc	<- 1						# number of processors on node
	hpc.walltime<- 23						# walltime
	hpc.q		<- NA #"pqeelab"						# PBS queue
	hpc.mem		<- "4gb" 					# RAM
	hpc.array	<- length(unique(df$CASE_ID))	# number of runs for job array
	pbshead		<- "#!/bin/sh"
	tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
	pbshead		<- paste(pbshead, tmp, sep = "\n")
	tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
	pbshead 	<- paste(pbshead, tmp, sep = "\n")
	pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
	if(!is.na(hpc.array))
		pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')		
	if(!is.na(hpc.q)) 
		pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
	pbshead 	<- paste(pbshead, hpc.load, sep = "\n")			
	#	make array job
	cmd		<- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
	cmd		<- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]			
	cmd		<- paste(pbshead,cmd,sep='\n')	
	#	submit job
	outfile		<- gsub(':','',paste("phs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
	outfile		<- file.path(outdir, outfile)
	cat(cmd, file=outfile)
	cmd 		<- paste("qsub", outfile)
	cat(cmd)
	cat(system(cmd, intern= TRUE))	
}	

## ---- rmd.chunk.seattle.191017.treedater.rerun ----
seattle.191017.treedater.rerun <- function()
{
	require(big.phylo)
	require(data.table)
	require(tidyverse)
	require(ape)
	require(phytools)
	require(treedater)
	require(phyloscannerR)
	
	
	home <- '~/Box Sync/OR_Work/Seattle/analysis_191017'
	home <- '/rds/general/project/ratmann_seattle_data_analysis/live/olli_191017'
	
	#	determine failed runs
	indir <- file.path(home,'phyloscanner_dated')
	infiles <- data.table(F_PHSC=list.files(indir, pattern='collapsed_workspace.rda$', full.names=TRUE, recursive=TRUE))
	infiles[, BASENAME:= gsub('collapsed_workspace.rda$','',basename(F_PHSC))]
	tmp <- data.table(F_TREE=list.files(indir, pattern='collapsed_dated.newick$', full.names=TRUE, recursive=TRUE))
	tmp[, BASENAME:= gsub('collapsed_dated.newick$','',basename(F_TREE))]
	infiles <- merge(infiles, tmp, by='BASENAME', all.x=TRUE)	
	infiles <- subset(infiles, is.na(F_TREE))
	outdir <- indir
	alignment.length <- 1000
		
	#
	#	for each tree: 
	#	make data.table of sequence sampling times 
	#	remove taxa without data on sampling times
	#
	cmds <- vector('list',nrow(infiles))
	for(i in seq_len(nrow(infiles)))
	{
		#	i<- 1
		cat('\nProcess',i)
		
		outfile.collapsed.phsc <- infiles[i,F_PHSC]
		outfile.tree <- file.path(outdir, gsub('collapsed_workspace\\.rda','collapsed\\.newick',basename(outfile.collapsed.phsc)) )
		outfile.sampling.times.bounds <- file.path(outdir, gsub('collapsed_workspace\\.rda','sampling_times_bounds\\.csv',basename(outfile.collapsed.phsc)))
		outfile.sampling.times.init <- file.path(outdir, gsub('collapsed_workspace\\.rda','sampling_times_init\\.csv',basename(outfile.collapsed.phsc)))
		control <- list( outfile=gsub('.newick$','_dated.newick',outfile.tree), 
				ali.len=alignment.length, 
				root=NA, 
				omega0=NA, 
				temporalConstraints=TRUE, 
				strictClock=FALSE, 
				estimateSampleTimes=outfile.sampling.times.bounds 
		)
		cmd <- cmd.treedater.script(outfile.tree, outfile.sampling.times.init, control=control)
		cmds[[i]] <- cmd
		
		cat('\n',Ntip(read.tree(outfile.tree)))
	}	
	infiles[, CMD:= unlist(cmds)]	
	
	#	submit jobs like this one:
	cat(infiles[1,CMD])
	
	#	run on HPC as array job
	df <- copy(infiles)
	df[, CASE_ID:= 1:nrow(df)]	
	#	make PBS header
	hpc.load	<- "module load R/3.3.3" # "module load anaconda3/personal"
	hpc.select	<- 1						# number of nodes
	hpc.nproc	<- 1						# number of processors on node
	hpc.walltime<- 123						# walltime
	hpc.q		<- "pqeelab"						# PBS queue
	hpc.mem		<- "12gb" 					# RAM
	hpc.array	<- length(unique(df$CASE_ID))	# number of runs for job array
	pbshead		<- "#!/bin/sh"
	tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
	pbshead		<- paste(pbshead, tmp, sep = "\n")
	tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
	pbshead 	<- paste(pbshead, tmp, sep = "\n")
	pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
	if(!is.na(hpc.array))
		pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')		
	if(!is.na(hpc.q)) 
		pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
	pbshead 	<- paste(pbshead, hpc.load, sep = "\n")			
	#	make array job
	cmd		<- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
	cmd		<- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]			
	cmd		<- paste(pbshead,cmd,sep='\n')	
	#	submit job
	outfile		<- gsub(':','',paste("phs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
	outfile		<- file.path(outdir, outfile)
	cat(cmd, file=outfile)
	cmd 		<- paste("qsub", outfile)
	cat(cmd)
	cat(system(cmd, intern= TRUE))	
}	


## ---- rmd.chunk.seattle.191017.treedater.postprocess ----
seattle.191017.treedater.postprocess <- function()
{
	require(ape)
	require(data.table)	
	require(phyloscannerR)
	
	home <- '~/Box Sync/OR_Work/Seattle/analysis_191017'
	#home <- '/rds/general/project/ratmann_seattle_data_analysis/live/olli_191017'		
	indir <- file.path(home,'phyloscanner_dated')
	infiles <- data.table(F_PHSC=list.files(indir, pattern='collapsed_workspace.rda$', full.names=TRUE, recursive=TRUE))
	infiles[, BASENAME:= gsub('collapsed_workspace.rda$','',basename(F_PHSC))]
	tmp <- data.table(F_TREE=list.files(indir, pattern='collapsed_dated.newick$', full.names=TRUE, recursive=TRUE))
	tmp[, BASENAME:= gsub('collapsed_dated.newick$','',basename(F_TREE))]
	infiles <- merge(infiles, tmp, by='BASENAME', all.x=TRUE)	
	infiles <- subset(infiles, !is.na(F_TREE))	
	stopifnot( nrow(subset(infiles, is.na(F_TREE)))==0 )
	
	for(i in seq_len(nrow(infiles)))
	{
		#i<- 1
		cat('\nProcess ',infiles[i,F_TREE])
		#	load data
		load(infiles[i,F_PHSC])
		ph.dated <- read.tree(infiles[i,F_TREE])		
		stopifnot( all(  ph.dated$tip.label == ph$tip.label ) )
		stopifnot( all( ph.dated$edge == ph$edge ) )
		
		#	since the tree topology is unchanged, we can copy
		#	the branch lenghts in units of time onto the original tree
		#	that has the ancestral state reconstructions
		ph$edge.length <- ph.dated$edge.length
		
		#	plot dated tree to spot obvious errors
		outfile <- gsub('collapsed_workspace\\.rda','annotated_dated_tree.pdf',infiles[i,F_PHSC])
		tmp <- vector('list')
		tmp[['tree']] <- ph
		tmp[['tree']][['node.states']] <- tmp[['tree']][['mapped.edge']] <- tmp[['tree']][['maps']] <- NULL
		attr(tmp[['tree']],'map.order') <- NULL
		attr(tmp[['tree']],'class') <- 'phylo'
		tmp[['read.counts']] <- rep(1, Ntip(ph))	
		write.annotated.tree(tmp, outfile, format="pdf", pdf.scale.bar.width = 0.01, pdf.w = 40, pdf.hm = 0.2, verbose = FALSE)
		
		#	save phyloscanner.tree
		outfile <- gsub('collapsed_workspace\\.rda','annotated_dated_tree.rda',infiles[i,F_PHSC])
		save(ph, file=outfile)
	}	
}

## ---- rmd.chunk.seattle.191017.sequence.labels ----	
seattle.191017.sequence.labels<- function()
{
	require(data.table)
	require(ape)
	require(tidyverse)
	
	home <- '~/Box Sync/OR_Work/Seattle'
	infile.indinfo <- file.path(home,'PHSKC-2018-07-09','person.rds')
	infile.seqinfo <- file.path(home,'PHSKC-2018-07-09','sequences_meta.rds')	
	infile.countryinfo <- file.path(home,'analysis_191017','misc/Country_db.rds')
	infiles.lanl <- file.path(home,'analysis_191017','alignments',c('180709_LANL_Subtype01AE_mafft.fasta','180709_LANL_SubtypeA1_mafft.fasta','180709_LANL_SubtypeB_mafft.fasta','180709_LANL_SubtypeC_mafft.fasta'))
	infile.subtype <- file.path(home,'analysis_191017','misc','180709_Subtype.csv')
	outfile.base <- file.path(home,'analysis_191017','misc','180709_')
	
	#
	# collect country codes
	#
	dco <- readRDS(infile.countryinfo)
	dco <- lapply(seq_len(length(dco)), function(i) tibble(WRLD=names(dco)[i], CNTRY=dco[[i]]))
	dco <- do.call('rbind',dco)
	dco <- dco %>% mutate(WRLD=gsub('Unkno','Unknown',gsub('\\.|\\_','',WRLD)))	
	tmp <- tibble(	CNTRY= c("BO", "SD", "EC", "AG", "DM", "GD", "GY", "VC", "SR", "LC", "UZ", "AZ", "GH", "HT", "LV", "BZ", "MN", "GF", "TJ"),
					NAME= c("Bolivia","Sudan","Ecuador","Antigua","Dominica","Grenada","Guyana","Saint Vincent","Suriname","Saint Lucia","Uzbekistan","Azerbaijan","Ghana","Haiti","Latvia","Belize","Mongolia","French Guiana","Tajikistan"),
					WRLD= c("SMAm","NorAf","SMAm","SMAm","SMAm","SMAm","SMAm","SMAm","SMAm","SMAm","Asia","Asia","SSA","SMAm","Europe","SMAm","Asia","SMAm","Asia")
					)
	dco <- tmp %>%	
		select(WRLD, CNTRY) %>%			
		rbind(dco) %>% 
		mutate(WRLD:= case_when(CNTRY=='MX'~"SMAm",CNTRY=='CA'~'Canada',CNTRY=='GL'~'Canada',CNTRY!='GL'&CNTRY!='CA'&CNTRY!='MX'~WRLD)) %>%		
		arrange(WRLD, CNTRY)
	table(dco$WRLD, dco$CNTRY)	

	#
	# read Seattle indidivual data 
	#	
	dind <- readRDS(infile.indinfo) 
	dind <- dind %>% 
			select(newnum, race, b_yr, birthCountry, gender, sex, sex_with_male_and_female, transm, rsh_county_name, rsa_county_name, cur_county_name) %>%
			mutate( birthCountry2:= gsub('^\\(([A-Z0-9]+)\\).*$','\\1', birthCountry),
					county:= case_when(	rsh_county_name!="KING CO."&rsa_county_name!="KING CO."&cur_county_name!="KING CO."~'other',
										rsh_county_name=="KING CO."|rsa_county_name=="KING CO."|cur_county_name=="KING CO."~'king'),
					Gender2:= case_when(gender=='FM'~'Trns',
										gender=='MF'~'Trns',
										gender=='F'&sex=='F'~'F',
										gender=='F'&sex=='M'~'Trns',
										gender=='M'&sex=='M'~'M',
										gender=='M'&sex=='F'~'Trns'),
					transm2:= case_when(grepl('MSM',transm)~'MSM',
										grepl('HETERO',transm)~'HSX',
										grepl('BLOOD|PERINAT|OTHER|IDU',transm)~'OTH',
										grepl('UNKNOWN',transm)~'UNKNOWN'),
					race2:= case_when( grepl('Black',race)~'Black',
									grepl('White',race)~'White',
									grepl('\\(1\\)Hispanic',race)~'Hispanic',
									grepl('Hawaiian|Indian|Asian|Multi-race',race)~'Other',
									grepl('Unknown',race)~'Unknown'
									)) %>%
			mutate( birthCountry2:= case_when(	birthCountry2==''~'Unknown',
												birthCountry2=='X98'~'Unknown',
												birthCountry2=='X99'~'Unknown',
												birthCountry2!=''&birthCountry2!='X98'&birthCountry2!='X99'~birthCountry2
												))
	#	setdiff( dind$birthCountry2, dco$CNTRY )	
	#	let s not resolve country of birth to world region for now..
	for(x in colnames(dind))
	 attr(dind[[x]], "label") <- NULL			
	#			
	#table(dind$Gender2, dind$transm2, dind$sex_with_male_and_female) # sex_with_male_and_female does not look very reliable

	#
	# read Seattle sequence data 
	#	
	dseq <- readRDS(infile.seqinfo)
	dseq <- dseq %>% filter(type=='PR/RT') %>% select(seqID, newnum, seqy)
	for(x in colnames(dseq))
		attr(dseq[[x]], "label") <- NULL	
	dseq <- dind %>% inner_join(dseq, by='newnum')	
	dseq <- dseq %>% select(newnum, seqID, seqy, county, Gender2, transm2, race2, birthCountry2)	
	tmp <- as_tibble(read.csv(infile.subtype, stringsAsFactors=FALSE)) %>%
			rename(seqID=name, ST=subtype) %>%
			select(seqID, ST) %>%
			mutate(seqID= gsub('PRRT','PR/RT',seqID))
	dseq <- dseq %>% left_join(tmp, by='seqID')
	
	#
	#	collect LANL labels
	#
	dl <- sapply(seq_along(infiles.lanl), function(i) rownames(read.dna(infiles.lanl[i],format='fa')))
	dl <- tibble(TAXA=unlist(dl)) %>% 
			filter(!grepl('PR\\/RT',TAXA) & TAXA!='HXB2') %>% 
			mutate(	ST:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',1),
					CNTRY:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',2),
					YEAR:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',3),
					GENBANK:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',5)) %>% 
			distinct()
	dl <- dl %>%
			filter(grepl('HXB2',TAXA)) %>%
			mutate(TAXA:='HXB2') %>%
			rbind(dl)	
	#
	outfile <- paste0(outfile.base,'sequence_labels.rda')
	save(dseq, dind, dl, dco, file= outfile)
	
	#
	#	make sequence labels for KCMSM
	#
	tmp <- dseq %>% rename(TAXA:= seqID) %>%
			mutate(WRLD:= case_when(county=='king'&transm2=='MSM'~'KCMSM',
										county=='king'&transm2!='MSM'~'KCnonMSM',
										county!='king'~'WAState')) %>%
			mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA,'_',newnum,'_',county,'_',race2,'_',birthCountry2,'_',Gender2,'_',transm2,'_',seqy)) %>%
			select(ST,WRLD,TAXA,TAXA_NEW)
	tmp <- dl %>% left_join(dco, by='CNTRY') %>%
			mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA)) %>%
			select(ST,WRLD,TAXA,TAXA_NEW) %>%
			rbind(tmp) %>%
			arrange(WRLD)
	outfile <- paste0(outfile.base,'sequence_labels_KCMSM.csv')
	write.csv(tmp, file=outfile)
	#
	#	make sequence labels for KCHSX
	#
	tmp <- dseq %>% rename(TAXA:= seqID) %>%
			mutate(WRLD:= case_when(county=='king'&transm2=='HSX'~'KCHSX',
							county=='king'&transm2!='HSX'~'KCnonHSX',
							county!='king'~'WAState')) %>%
			mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA,'_',newnum,'_',county,'_',race2,'_',birthCountry2,'_',Gender2,'_',transm2,'_',seqy)) %>%
			select(ST,WRLD,TAXA,TAXA_NEW)
	tmp <- dl %>% left_join(dco, by='CNTRY') %>%
			mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA)) %>%
			select(ST,WRLD,TAXA,TAXA_NEW) %>%
			rbind(tmp) %>%
			arrange(WRLD)
	outfile <- paste0(outfile.base,'sequence_labels_KCHSX.csv')
	write.csv(tmp, file=outfile)
	#
	#	make sequence labels for KC overall
	#
	tmp <- dseq %>% rename(TAXA:= seqID) %>%
			mutate(WRLD:= case_when(county=='king'~'KC',							
							county!='king'~'WAState')) %>%
			mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA,'_',newnum,'_',county,'_',race2,'_',birthCountry2,'_',Gender2,'_',transm2,'_',seqy)) %>%
			select(ST,WRLD,TAXA,TAXA_NEW)
	tmp <- dl %>% left_join(dco, by='CNTRY') %>%
			mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA)) %>%
			select(ST,WRLD,TAXA,TAXA_NEW) %>%
			rbind(tmp) %>%
			arrange(WRLD)
	outfile <- paste0(outfile.base,'sequence_labels_KC.csv')
	write.csv(tmp, file=outfile)
}

## ---- rmd.chunk.seattle.191017.phyloscanner.nonB ----
seattle.191017.phyloscanner.nonB <- function()
{
	require(data.table)
	require(ape)
	require(adephylo)
	require(phytools)
	
	plot.phylogenies <- 0
	if(0)
	{
		prog.phyloscanner_analyse_trees <- '/rds/general/user/or105/home/libs_sandbox/phyloscanner/phyloscanner_analyse_trees.R'
		home <- '/rds/general/project/ratmann_seattle_data_analysis/live/olli_191017'
	}
	if(1)
	{
		prog.phyloscanner_analyse_trees <- '/Users/Oliver/git/phyloscanner/phyloscanner_analyse_trees.R'
		home <- '/Users/or105/Box/OR_Work/Seattle/analysis_191017'	
	}
	
	
	
	#
	
	indir.trees <- file.path(home,'fasttree')
	infile.labels <- data.table(FLABEL= c( 	file.path(home,'misc','180709_sequence_labels_KC.csv'),
										 	file.path(home,'misc','180709_sequence_labels_KCMSM.csv'),
											file.path(home,'misc','180709_sequence_labels_KCHSX.csv'))
								)
	outdir <- file.path(home,'phyloscanner')
	
	
	infiles <- data.table(FIN=list.files(indir.trees, pattern='\\.newick$',full.names=TRUE))
	infiles[, DUMMY:= 1L]
	infile.labels[, DUMMY:= 1L]
	infiles <- merge(infiles, infile.labels, by='DUMMY', allow.cartesian=TRUE)
	infiles[, DUMMY:= NULL]
	infiles[, ST:= gsub('.*_Subtype([A-Z0-9]+)_.*','\\1',basename(FIN))]	
	infiles[, SELECT:= gsub('^.*_labels_([A-Z]+)\\.csv$','\\1', basename(FLABEL))]	
	infiles[, DUMMY:= paste0(gsub('\\.newick$','_rerooted_',basename(infiles$FIN)), infiles$SELECT)]
	tmp <- data.table(FOUT=list.files(outdir, pattern='workspace.rda$',full.names=TRUE))
	tmp[, DUMMY:= gsub('__workspace.rda$','',basename(FOUT))]
	infiles <- merge(infiles, tmp, by='DUMMY', all.x=TRUE)
	infiles <- subset(infiles, is.na(FOUT) & ST!='B', c(ST, FIN, FLABEL))	
	cmds <- vector('list',nrow(infiles))
	for(i in seq_len(nrow(infiles)))
	{
		#i<- 1
		cat('\nprocess ',infiles[i,FIN],'\nprocess ',infiles[i,FLABEL])
		infile <- infiles[i,FIN]
		ph <- read.tree(infile)
		ph$node.label <- NULL
		#	update tip labels: add world region to start of label
		local.world.region <- gsub('^.*_labels_([A-Z]+)\\.csv$','\\1', basename(infiles[i,FLABEL]))
		dl <- as.data.table(read.csv(infiles[i,FLABEL], stringsAsFactors=FALSE))
		dl[, X:=NULL]	
		stopifnot(!any(ph$tip.label==''))
		dp <- data.table(IDX=1:Ntip(ph), TAXA=ph$tip.label)
		dp <- merge(dp,dl,by='TAXA',all.x=TRUE)
		stopifnot( nrow(dp[is.na(TAXA_NEW),])==0 )
		stopifnot( !any(duplicated(ph$tip.label)) )
		ph$tip.label <- dp[order(IDX),][, TAXA_NEW]
		#	re-root
		tmp <- dp[, list(NST=length(TAXA)), by='ST']
		tmp <- tmp[order(-NST),][2,ST]
		tmp <- subset(dp, ST==tmp,)[,TAXA_NEW]
		root <- getMRCA(ph, tmp)
		ph <- reroot(ph, root, ph$edge.length[which(ph$edge[,2]==root)]/2)
		#	drop other subtypes
		tmp <- dp[, list(NST=length(TAXA)), by='ST']
		tmp <- tmp[order(-NST),][-1,ST]
		tmp <- subset(dp, ST%in%tmp,)[,TAXA_NEW]
		ph <- drop.tip(ph, tmp)
		stopifnot(is.binary(ph))
		#	write to file
		intree.phsc <- file.path(outdir,gsub('\\.newick',paste0('_rerooted_',local.world.region,'.newick'),basename(infile)))
		write.tree(ph, file=intree.phsc)
		if(plot.phylogenies)
		{
			pdf(file=gsub('newick','pdf',intree.phsc), w=20, h=10+Ntip(ph)/10)
			plot(ph, show.node.label=TRUE, cex=0.3)
			dev.off()			
		}
		# 	make phyloscanner UNIX command
		infile <- intree.phsc
		outputString <- paste0(gsub('\\.newick','_',intree.phsc))
		tip.regex <- "^([A-Za-z]+)___.*$"
		cmd <- paste("CWD=$(pwd)\n",sep='')
		cmd <- paste(cmd,"echo $CWD\n",sep='')	
		tmpdir.prefix <- paste('phsc_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
		tmpdir <- paste("$CWD/",tmpdir.prefix,sep='')
		tmp.in <- file.path(tmpdir, basename(infile))
		cmd <- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
		cmd <- paste(cmd,'cp "',infile,'" ',tmp.in,'\n', sep='')
		cmd <- paste(cmd,'cd ', tmpdir,'\n', sep='')
		cmd <- paste0(cmd,'Rscript ',prog.phyloscanner_analyse_trees,' ',basename(infile),' ',basename(outputString))
		#	don t use -m multifurcation threshold to ensure that output tree is still binary
		cmd <- paste0(cmd,' s,0 -x "',tip.regex,'" -v 1 -ow -rda\n')
		cmd <- paste0(cmd,'mv ',basename(outputString),'* ','"',dirname(outputString),'"','\n')	
		cmd <- paste(cmd, "cd $CWD\n",sep='')
		cmd <- paste(cmd, "rm ", tmpdir,'\n',sep='')
		cmds[[i]] <- cmd
	}
	infiles[, CMD:= unlist(cmds)]	
	
	#	submit jobs like this one:
	cat(infiles[1,CMD])
	
	#	run on HPC as array job
	df <- copy(infiles)
	df[, CASE_ID:= 1:nrow(df)]	
	#	make PBS header
	hpc.load	<- "module load anaconda3/personal"
	hpc.select	<- 1						# number of nodes
	hpc.nproc	<- 1						# number of processors on node
	hpc.walltime<- 23						# walltime
	hpc.q		<- NA #"pqeelab"						# PBS queue
	hpc.mem		<- "4gb" 					# RAM
	hpc.array	<- length(unique(df$CASE_ID))	# number of runs for job array
	pbshead		<- "#!/bin/sh"
	tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
	pbshead		<- paste(pbshead, tmp, sep = "\n")
	tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
	pbshead 	<- paste(pbshead, tmp, sep = "\n")
	pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
	if(!is.na(hpc.array))
		pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')		
	if(!is.na(hpc.q)) 
		pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
	pbshead 	<- paste(pbshead, hpc.load, sep = "\n")			
	#	make array job
	cmd		<- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
	cmd		<- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]			
	cmd		<- paste(pbshead,cmd,sep='\n')	
	#	submit job
	outfile		<- gsub(':','',paste("phs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
	outfile		<- file.path(outdir, outfile)
	cat(cmd, file=outfile)
	cmd 		<- paste("qsub", outfile)
	cat(cmd)
	cat(system(cmd, intern= TRUE))	
}

## ---- rmd.chunk.seattle.191017.phyloscanner.B ----
seattle.191017.phyloscanner.B <- function()
{
	require(data.table)
	require(ape)
	require(adephylo)
	require(phytools)
	require(phangorn)
	
	
	plot.phylogenies <- 1
	max.Ntip <- 5e3	
	if(1)
	{
		prog.phyloscanner_analyse_trees <- '/Users/Oliver/git/phyloscanner/phyloscanner_analyse_trees.R'
		home <- '/Users/or105/Box/OR_Work/Seattle/analysis_191017'	
	}
	if(0)
	{
		prog.phyloscanner_analyse_trees <- '/rds/general/user/or105/home/libs_sandbox/phyloscanner/phyloscanner_analyse_trees.R'
		home <- '/rds/general/project/ratmann_seattle_data_analysis/live/olli_191017'	
	}
	
	
	indir.trees <- file.path(home,'fasttree')
	infile.labels <- data.table(FLABEL= c( 	file.path(home,'misc','180709_sequence_labels_KC.csv'),
					file.path(home,'misc','180709_sequence_labels_KCMSM.csv'),
					file.path(home,'misc','180709_sequence_labels_KCHSX.csv'))
					)
	outdir <- file.path(home,'phyloscanner')
	infiles <- data.table(FIN=list.files(indir.trees, pattern='\\.newick$',full.names=TRUE))
	infiles[, DUMMY:= 1L]
	infile.labels[, DUMMY:= 1L]
	infiles <- merge(infiles, infile.labels, by='DUMMY', allow.cartesian=TRUE)
	infiles[, DUMMY:= NULL]
	infiles[, ST:= gsub('.*_Subtype([A-Z0-9]+)_.*','\\1',basename(FIN))]	
	infiles <- subset(infiles, ST=='B')
	set.seed(42L)
	for(i in seq_len(nrow(infiles)))
	{
		#i<- 1
		cat('\nprocess ',infiles[i,FIN],'\nprocess ',infiles[i,FLABEL])
		infile <- infiles[i,FIN]
		ph <- read.tree(infile)
		ph$node.label <- NULL
		#	update tip labels: add world region to start of label
		local.world.region <- gsub('^.*_labels_([A-Z]+)\\.csv$','\\1', basename(infiles[i,FLABEL]))
		dl <- as.data.table(read.csv(infiles[i,FLABEL], stringsAsFactors=FALSE))
		dl[, X:=NULL]	
		stopifnot(!any(ph$tip.label==''))
		dp <- data.table(IDX=1:Ntip(ph), TAXA=ph$tip.label)
		dp <- merge(dp,dl,by='TAXA',all.x=TRUE)
		stopifnot( nrow(dp[is.na(TAXA_NEW),])==0 )
		stopifnot( !any(duplicated(ph$tip.label)) )
		ph$tip.label <- dp[order(IDX),][, TAXA_NEW]
		#	re-root
		tmp <- dp[, list(NST=length(TAXA)), by='ST']
		tmp <- tmp[order(-NST),][2,ST]
		tmp <- subset(dp, ST==tmp,)[,TAXA_NEW]
		root <- getMRCA(ph, tmp)
		ph <- reroot(ph, root, ph$edge.length[which(ph$edge[,2]==root)]/2)
		#	drop other subtypes
		tmp <- dp[, list(NST=length(TAXA)), by='ST']
		tmp <- tmp[order(-NST),][-1,ST]
		tmp <- subset(dp, ST%in%tmp,)[,TAXA_NEW]
		ph <- drop.tip(ph, tmp)
		#	split very large tree into small trees if necessary
		if(is.finite(max.Ntip))
		{
			split.phs <- list()
			repeat({
					local.tips <- which(grepl(paste0('^',local.world.region),ph$tip.label))
					local.tip.ancestors <- Ancestors(ph, sample(local.tips,1))
					cat('\nNumer of local tips left to divide into smaller trees, n=', length(local.tips))
					subtree.mrca <- Ntip(ph)+1L
					for(k in seq_along(local.tip.ancestors))
					{				
						tmp <- length(Descendants(ph, local.tip.ancestors[k], type='tips')[[1]])						
						if(tmp>max.Ntip)
						{
							subtree.mrca <- local.tip.ancestors[max(1,k-1)]
							break
						}
					}					
					stopifnot( subtree.mrca>Ntip(ph) )
					tmp <- extract.clade(ph, subtree.mrca)
					if(Ntip(ph)==Ntip(tmp))
					{
						split.phs[[length(split.phs)+1L]] <- tmp
						break
					}
					if(Ntip(ph)-Ntip(tmp)>50)					
					{						
						split.phs[[length(split.phs)+1L]] <- tmp
						ph <- drop.tip(ph, split.phs[[length(split.phs)]][['tip.label']])						
					}
					if(Ntip(ph)-Ntip(tmp)<=50)
					{
						cat('\nFound very small final tree, retry sampling a local.tip')	
					}					
				})			
		}
		if(!is.finite(max.Ntip))
			split.phs[[1]] <- ph
		#	write to files
		if(length(split.phs)>1)
		{
			for(k in seq_along(split.phs))
			{
				intree.phsc <- gsub('_(Subtype[A-Za-z0-9]+)_',paste0('_\\1c',k,'_'),gsub('\\.newick',paste0('_rerooted_',local.world.region,'.newick'),basename(infile)))
				intree.phsc <- file.path(outdir,intree.phsc)
				stopifnot(is.binary(split.phs[[k]]))
				write.tree(split.phs[[k]], file=intree.phsc)
				if(plot.phylogenies)
				{
					pdf(file=gsub('newick','pdf',intree.phsc), w=20, h=10+Ntip(split.phs[[k]])/10)
					plot(split.phs[[k]], show.node.label=TRUE, cex=0.3)
					dev.off()			
				}
			}							
		}
		if(length(split.phs)==1)
		{
			intree.phsc <- file.path(outdir,gsub('\\.newick',paste0('_rerooted_',local.world.region,'.newick'),basename(infile)))
			write.tree(split.phs[[1]], file=intree.phsc)
			if(plot.phylogenies)
			{
				pdf(file=gsub('newick','pdf',intree.phsc), w=20, h=10+Ntip(split.phs[[1]])/10)
				plot(split.phs[[1]], show.node.label=TRUE, cex=0.3)
				dev.off()			
			}
		}
	}
	
	infiles <- data.table(FIN=list.files(outdir, pattern='\\.newick$',full.names=TRUE))
	infiles[, ST:= gsub('^.*_Subtype([A-Za-z0-9]+)_.*$','\\1',basename(FIN))]
	infiles[, DUMMY:= gsub('\\.newick$','',basename(FIN))]
	tmp <- data.table(FOUT=list.files(outdir, pattern='workspace.rda$',full.names=TRUE))
	tmp[, DUMMY:= gsub('__workspace.rda$','',basename(FOUT))]
	infiles <- merge(infiles, tmp, by='DUMMY', all.x=TRUE)
	infiles <- subset(infiles, is.na(FOUT), c(ST, FIN))
	infiles <- subset(infiles, grepl('Bc[0-9]+',ST))
	cmds <- vector('list',nrow(infiles))	
	for(i in seq_len(nrow(infiles)))
	{				
		# 	make phyloscanner UNIX command
		infile <- infiles[i,FIN]
		outputString <- paste0(gsub('\\.newick','_',infile))
		tip.regex <- "^([A-Za-z]+)___.*$"
		cmd <- paste("CWD=$(pwd)\n",sep='')
		cmd <- paste(cmd,"echo $CWD\n",sep='')	
		tmpdir.prefix <- paste('phsc_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
		tmpdir <- paste("$CWD/",tmpdir.prefix,sep='')
		tmp.in <- file.path(tmpdir, basename(infile))
		cmd <- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
		cmd <- paste(cmd,'cp "',infile,'" ',tmp.in,'\n', sep='')
		cmd <- paste(cmd,'cd ', tmpdir,'\n', sep='')
		cmd <- paste0(cmd,'Rscript ',prog.phyloscanner_analyse_trees,' ',basename(infile),' ',basename(outputString))
		#	don t use -m multifurcation threshold to ensure that output tree is still binary
		cmd <- paste0(cmd,' s,0 -x "',tip.regex,'" -v 1 -ow -rda\n')
		cmd <- paste0(cmd,'mv ',basename(outputString),'* ','"',dirname(outputString),'"','\n')	
		cmd <- paste(cmd, "cd $CWD\n",sep='')
		cmd <- paste(cmd, "rm ", tmpdir,'\n',sep='')
		cmds[[i]] <- cmd
	}
	infiles[, CMD:= unlist(cmds)]	
	
	#	submit jobs like this one:
	cat(infiles[1,CMD])
	
	#	run on HPC as array job
	df <- copy(infiles)
	df[, CASE_ID:= 1:nrow(df)]	
	#	make PBS header
	hpc.load	<- "module load anaconda3/personal"
	hpc.select	<- 1						# number of nodes
	hpc.nproc	<- 1						# number of processors on node
	hpc.walltime<- 23						# walltime
	hpc.q		<- NA #"pqeelab"						# PBS queue
	hpc.mem		<- "4gb" 					# RAM
	hpc.array	<- length(unique(df$CASE_ID))	# number of runs for job array
	pbshead		<- "#!/bin/sh"
	tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
	pbshead		<- paste(pbshead, tmp, sep = "\n")
	tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
	pbshead 	<- paste(pbshead, tmp, sep = "\n")
	pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
	if(!is.na(hpc.array))
		pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')		
	if(!is.na(hpc.q)) 
		pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
	pbshead 	<- paste(pbshead, hpc.load, sep = "\n")			
	#	make array job
	cmd		<- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
	cmd		<- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]			
	cmd		<- paste(pbshead,cmd,sep='\n')	
	#	submit job
	outfile		<- gsub(':','',paste("phs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
	outfile		<- file.path(outdir, outfile)
	cat(cmd, file=outfile)
	cmd 		<- paste("qsub", outfile)
	cat(cmd)
	cat(system(cmd, intern= TRUE))	
}

## ---- rmd.chunk.seattle.190723.get.newick.from.phsc ----
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

## ---- rmd.chunk.seattle.190723.phyloscanner.make.pdf ----
seattle.191017.phyloscanner.make.pdf<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	require(reshape)
	require(phyloscannerR)
	#	working directory with phyloscanner output		
	indir.phsc	<- '/Users/or105/Box/OR_Work/Seattle/analysis_191017/phyloscanner'
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
}

## ---- rmd.chunk.seattle.190723.get.undated.subgraphs ----
seattle.191017.phyloscanner.get.undated.subgraphs<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	require(reshape)
	require(phyloscannerR)
	#	working directory with phyloscanner output		
	indir.phsc	<- '/Users/or105/Box/OR_Work/Seattle/analysis_191017/phyloscanner'
	infiles		<- data.table(F=list.files(indir.phsc, pattern='workspace.rda$', full.names=TRUE, recursive=TRUE))
	infiles[, SELECT:= 'KCHSX']	
	
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
}

## ---- rmd.chunk.seattle.191017.get.dated.subgraphs ----
seattle.191017.phyloscanner.get.dated.subgraphs<- function()
{
	require(data.table)
	require(phangorn)
	require(ggplot2)
	require(reshape)
	require(phyloscannerR)
	
	#	working directory with phyloscanner output
	home <- '/Users/or105/Box/OR_Work/Seattle/analysis_191017'
	indir.phsc	<- file.path(home,'phyloscanner_dated')
	outdir <- file.path(home,'subgraphs_dated')
	infiles		<- data.table(F=list.files(indir.phsc, pattern='_annotated_dated_tree.rda$', full.names=TRUE, recursive=TRUE))	
	infiles[, SELECT:= gsub('^.*_rerooted_([A-Za-z0-9]+)_.*$','\\1',basename(F))]	
	
	#	extract subgraphs of dated trees
	for(i in seq_len(nrow(infiles)))
	{
		#i<- 21
		cat('process', i,'\n')
		infile <- infiles[i, F]
		host <- infiles[i,SELECT]
		load(infile)					
		mrcas <- which( attr(ph, 'SUBGRAPH_MRCA') )
		#	some of the tips states are "unknown" and this clashes internally with the NA state, so we need to take some extra care
		#	this is not a problem because the "unknown" and NA state mean the same thing
		attr(ph, 'INDIVIDUAL') <- as.character(attr(ph, 'INDIVIDUAL'))
		attr(ph, 'INDIVIDUAL')[is.na(attr(ph, 'INDIVIDUAL'))] <- 'Unknown'
		mrcas <- mrcas[ attr(ph, 'INDIVIDUAL')[mrcas]==host ]
		stopifnot( !any(is.na(mrcas)) )		
		# check that tree is of class simmap
		stopifnot( any(attr(ph,'class')=='simmap') )
		# extract subgraphs
		subgraphs <- lapply(mrcas, function(mrca) extract.subgraph(ph, mrca))
		# save
		outfile <- gsub('_annotated_dated_tree',paste0('_datedsubgraphs_',host),basename(infile))
		outfile <- file.path(outdir,outfile)
		save(subgraphs, file=outfile)
	}	
}

## ---- rmd.chunk.mueller.tree.height.test ----
mueller.tree.height.test<- function()
{
	require(data.table)
	#	esco takes about 8 hours, masco about 1.5 hours per 1M samples and 1M needed
	#	need to run on HPC
	
	cmd.beast <- '"/Applications/BEAST 2.6.0/lib/launcher.jar"'
	cmd.beast <- '"/Users/Oliver/git/The-Structured-Coalescent/jar/esco.jar"'
	cmd.launcher <- paste0('java -Xmx16g -jar ',cmd.beast,' -threads 4')
	
	
	indir <- '~/Box Sync/OR_Work/Seattle/Mueller/TreeHeightTest/xmls'
	outdir <- '~/Box Sync/OR_Work/Seattle/Mueller/TreeHeightTest/out_191015'
	infiles <- data.table(FIN=list.files(indir, pattern='xml$', full.names=TRUE))
	
	i<- 4
	cmd <- paste0(cmd.launcher, ' "', infiles[i, FIN], '"')
	cat(cmd)	
}

## ---- rmd.chunk.mueller.aiv ----
mueller.aiv<- function()
{
	require(data.table)
	#	masco takes about 50 hours for 1M iterations and 10M needed
	#	need to run on HPC
	
	cmd.beast <- '"/Applications/BEAST 2.6.0/lib/launcher.jar"'
	cmd.beast <- '"/Users/Oliver/git/The-Structured-Coalescent/jar/esco.jar"'
	cmd.launcher <- paste0('java -Xmx16g -jar ',cmd.beast,' -threads 4')
	
	
	indir <- '~/Box Sync/OR_Work/Seattle/Mueller/AIV/xmls'
	outdir <- '~/Box Sync/OR_Work/Seattle/Mueller/AIV/out_191015'
	infiles <- data.table(FIN=list.files(indir, pattern='xml$', full.names=TRUE))
	
	i<- 2
	cmd <- paste0(cmd.launcher, ' "', infiles[i, FIN], '"')
	cat(cmd)	
}

## ---- rmd.chunk.cmd.beast2 ----
cmd.beast2<- function(infile.xml, control=list(beast.launcher='/Applications/BEAST 2.6.0/lib/launcher.jar', out.dir=dirname(infile.xml)))
{			
	stopifnot(is.character(control$beast.launcher))
	stopifnot(is.character(control$out.dir))
	stopifnot(grepl('\\.xml$',infile.xml))
	if(!any(names(control)=='Xmx'))
		control$Xmx <- '2g'
	if(!any(names(control)=='threads'))
		control$threads <- 1L
	cmd				<- paste("CWD=$(pwd)\n",sep='')
	cmd				<- paste(cmd,"echo $CWD\n",sep='')	
	tmpdir.prefix	<- paste('ft_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	tmpdir			<- paste("$CWD/",tmpdir.prefix,sep='')
	tmp.in			<- file.path(tmpdir, basename(infile.xml))
	tmp.out			<- gsub('\\.xml','',basename(infile.xml))	
	cmd				<- paste0(cmd,"mkdir -p ",tmpdir,'\n')
	cmd				<- paste0(cmd,'cp "',infile.xml,'" ',tmp.in,'\n')
	cmd				<- paste0(cmd,'cd ', tmpdir,'\n')
	cmd				<- paste0(cmd, 'java -Xmx',control$Xmx,' -jar "', control$beast.launcher,'" -threads ',control$threads)
	cmd				<- paste0(cmd,' ',basename(infile.xml),'\n')
	cmd				<- paste0(cmd, 'mv ', tmp.out,'* "',control$out.dir,'"\n')
	cmd				<- paste0(cmd,'cd $CWD\n')
	cmd				<- paste(cmd, "rm -r ", tmpdir,'\n',sep='')	
	cmd
}

## ---- rmd.chunk.siveroni.flu ----
siveroni.flu<- function()
{
	require(data.table)
	#	masco takes about 50 hours for 1M iterations and 10M needed
	#	need to run on HPC
	
	if(0)
	{
		beast.launcher <- '/Applications/BEAST 2.6.0/lib/launcher.jar'
		#beast.launcher <- '/Users/Oliver/git/The-Structured-Coalescent/jar/esco.jar'	
		in.dir <- '~/Box Sync/OR_Work/Seattle/Siveroni/Flu/xmls'
		out.dir <- '~/Box Sync/OR_Work/Seattle/Siveroni/Flu/out_191015'		
	}
	if(1)
	{
		beast.launcher <- '/rds/general/user/or105/home/libs/beast/lib/launcher.jar'
		in.dir <- '/rds/general/project/ratmann_seattle_data_analysis/live/siveroni/flu/xmls'
		out.dir <- '/rds/general/project/ratmann_seattle_data_analysis/live/siveroni/flu/out'		
	}
	
	#	create BEAST2 run
	infiles <- data.table(FIN=list.files(in.dir, pattern='xml$', full.names=TRUE))
	control <- list(beast.launcher=beast.launcher,
			threads= 1L,
			Xmx= '6g',
			out.dir= out.dir)	
	i<- 1
	cmd <- cmd.beast2(infiles[i, FIN], control=control)
	cat(cmd)	
	
	#	submit to HPC
	hpc.load	<- "module load java/jdk-8u66"
	hpc.select	<- 1						# number of nodes
	hpc.nproc	<- 1						# number of processors on node
	hpc.walltime<- 923						# walltime
	hpc.q		<- "pqeelab"				# PBS queue
	hpc.mem		<- "6gb" 					# RAM		
	pbshead		<- "#!/bin/sh"
	tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
	pbshead		<- paste(pbshead, tmp, sep = "\n")
	tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
	pbshead 	<- paste(pbshead, tmp, sep = "\n")
	pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
	if(!is.na(hpc.q)) 
		pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
	pbshead 	<- paste(pbshead, hpc.load, sep = "\n")			
	cmd			<- paste(pbshead,cmd,sep='\n')
	cat(cmd)	 								
	outfile		<- gsub(':','',paste("pydy",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
	outfile		<- file.path(out.dir, outfile)
	cat(cmd, file=outfile)
	cmd 		<- paste("qsub", outfile)		
	cat(system(cmd, intern= TRUE))	
	quit("no")		
}

## ---- rmd.chunk.seattle.191017.alignment.rm.drug.resistance.mutations ----
seattle.191017.alignment.rm.drug.resistance.mutations <- function()
{	
	require(big.phylo)
	home <- '~/Box Sync/OR_Work/Seattle/analysis_191017'
	indir <- file.path('alignments')	
	df <- data.table(F=list.files(indir, pattern='fasta$',full.names=TRUE))
	df <- subset(df, !grepl('nodrm',F))
	for(i in seq_len(nrow(df)))
	{
		infile.fasta <- df[i,F]
		outfile <- gsub('\\.fasta','_ndrm.fasta',infile.fasta)
		
		cat('\nprocess ',infile.fasta)
		sq					<- read.dna(infile.fasta, format='fa')	
		tmp					<- which(grepl("HXB2",rownames(sq)))
		rownames(sq)[tmp]	<- 'HXB2'
		tmp					<- big.phylo:::seq.rm.drugresistance(sq)
		sq					<- tmp$nodr.seq
		cat('\nwrite to ',outfile)
		write.dna(sq, file= outfile, format='fasta')
	}	
}

## ---- rmd.chunk.seattle.191017.alignment.update ----
seattle.191017.alignment.update <- function()
{	
	require(tidyverse)
	
	home <- '~/Box Sync/OR_Work/Seattle'
	infile.seqinfo.current <- file.path(home,'PHSKC-2018-07-09/sequences_meta.rds')
	infile.seq.current <- file.path(home,'PHSKC-2018-07-09/Sequences.fas')
	infile.seqinfo.update <- file.path(home,'Seattle/PHSKC-2019_10_11/sequences_meta.rds')
	infile.seq.update <- file.path(home,'PHSKC-2019_10_11/Sequences.fas')
	
	sc <- readRDS(infile.seqinfo.current) %>%
			filter(type!='IN') %>%
			select(seqID, newnum, seqy, seqm) %>%	
			rename(cur.newnum=newnum, cur.seqy=seqy, cur.seqm=seqm)			
	su <- readRDS(infile.seqinfo.update) %>%
			filter(type!='IN') %>%
			select(seqID, newnum, type, seqy, seqm) %>%
			rename(upd.newnum=newnum, upd.seqy=seqy, upd.seqm=seqm)
	ss <- sc %>% full_join(su, by='seqID')
	which(is.na(ss$upd.newnum))
	
	#	! FAIL !
	ss %>% filter(is.na(upd.newnum))
	
	#	ok
	ss %>% 
			filter(!is.na(upd.newnum) & !is.na(cur.newnum)) %>%			
			filter(cur.newnum!=upd.newnum | cur.seqy!=upd.seqy | cur.seqm!=upd.seqm ) 

	seq.ids <- ss %>% 
			filter(!is.na(upd.newnum) & !is.na(cur.newnum)) %>%
			pull(seqID)
	seqc <- read.dna(infile.seq.current, format='fa')
	sequ <- read.dna(infile.seq.update, format='fa')
	test <- sapply(seq_along(seq.ids), function(i){
				cur <- paste(unlist(as.character(seqc[ seq.ids[i] ])), collapse='')
				upd <- paste(unlist(as.character(sequ[ seq.ids[i] ])), collapse='')
				cur==upd
			})
	test <- tibble(seqID=seq.ids, identical=test)
	test.fail <- test %>% filter(!identical)
	test.fail
	
	seq.name <- test.fail[1, ] %>% pull(seqID )
	paste(unlist(as.character(seqc[ seq.name ])), collapse='')
	paste(unlist(as.character(sequ[ seq.name ])), collapse='')
	
		
}

## ---- rmd.chunk.seattle.191017.alignment.make.bootstraps ----
seattle.191017.alignment.make.bootstraps <- function()
{	
	home <- '~/Box Sync/OR_Work/Seattle/analysis_191017'
	indir <- file.path(home,'alignments')
	outdir <- file.path(home,'alignments_bs')
	bsn <- 10
	seed <- 42L
	#
	df <- data.table(F=list.files(indir, pattern='fasta$',full.names=TRUE))
	df <- subset(df, grepl('ndrm',F))
	for(i in seq_len(nrow(df)))
	{
		infile.fasta <- df[i,F]
		cat('\nprocess ',infile.fasta)
		sq					<- read.dna(infile.fasta, format='fa')	
		set.seed(seed)
		for(b in 0:bsn)
		{
			outfile <- gsub('\\.fasta',paste0('_',sprintf("%03d",b),'.fasta'),infile.fasta)
			outfile <- file.path(outdir,basename(outfile))
			bs <- 1:ncol(sq)
			if(b>0)
			{
				bs<- sample(1:ncol(sq), ncol(sq), replace=TRUE)	
			}			
			write.dna(sq[,bs], file=outfile, format='fa', colsep='', nbcol=-1)
		}		
	}	
}

## ---- rmd.chunk.seattle.191017.fastree ----
seattle.191017.fastree<- function()
{	
	require(big.phylo)
	require(data.table)
	
	home <- '~/Box Sync/OR_Work/Seattle/analysis_191017'
	home <- '/rds/general/project/ratmann_seattle_data_analysis/live/olli_191017'
	
	indir <- file.path(home,'alignments_bs')
	outdir <- file.path(home,'fasttree')
	
	#	get UNIX commands for all trees to be processed
	df <- data.table(FIN=list.files(indir, pattern='fasta$',full.names=TRUE))
	df <- subset(df, grepl('ndrm',FIN))
	df[, ST:= gsub('.*_Subtype([A-Z0-9]+)_.*','\\1',basename(FIN))]
	df[, FOUT:= file.path(outdir,gsub('\\.fasta','_ft.newick',basename(FIN)))]
	cmds <- vector('list', nrow(df))
	for(i in seq_len(nrow(df)))
	{
		infile.fasta <- df[i,FIN]
		outfile <- df[i,FOUT]			
		tmp <- cmd.fasttree(infile.fasta, outfile=outfile, pr.args='-nt -gtr -gamma', check.binary=TRUE)
		cmds[[i]] <- tmp
	}
	df[, CMD:= unlist(cmds)]	
	
	#	submit jobs like this one:
	cat(df[1,CMD])
	
	#	run on HPC as array job
	df <- subset(df, ST=='B')
	df[, CASE_ID:= 1:nrow(df)]
	
	#	make PBS header
	hpc.load	<- "module load R/3.3.3"
	hpc.select	<- 1						# number of nodes
	hpc.nproc	<- 1						# number of processors on node
	hpc.walltime<- 923						# walltime
	hpc.q		<- "pqeelab"						# PBS queue
	hpc.mem		<- "6gb" 					# RAM
	hpc.array	<- length(unique(df$CASE_ID))	# number of runs for job array
	pbshead		<- "#!/bin/sh"
	tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
	pbshead		<- paste(pbshead, tmp, sep = "\n")
	tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
	pbshead 	<- paste(pbshead, tmp, sep = "\n")
	pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
	if(!is.na(hpc.array))
		pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')		
	if(!is.na(hpc.q)) 
		pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
	pbshead 	<- paste(pbshead, hpc.load, sep = "\n")			
	#	make array job
	cmd		<- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
	cmd		<- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]			
	cmd		<- paste(pbshead,cmd,sep='\n')	
	#	submit job
	outfile		<- gsub(':','',paste("trs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
	outfile		<- file.path(outdir, outfile)
	cat(cmd, file=outfile)
	cmd 		<- paste("qsub", outfile)
	cat(cmd)
	cat(system(cmd, intern= TRUE))
}

## ---- rmd.chunk.seattle.levu.180613 ----
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

## ---- rmd.chunk.seattle.180423.make.subtype.alignments ----
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

## ---- rmd.chunk.seattle.170601.fastree ----
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

## ---- rmd.chunk.seattle.170621.fastree ----
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

## ---- rmd.chunk.seattle.170814.metadata.regionorigin ----
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

## ---- rmd.chunk.seattle.170814.LSD ----
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

## ---- rmd.chunk.seattle.170601.clustering ----
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

## ---- rmd.chunk.seattle.170607.clustering ----
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

## ---- rmd.chunk.seattle.170601.clustering.plot ----
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