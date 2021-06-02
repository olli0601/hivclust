## ---- rmd.chunk.athena.Sequences.LANL.preprocess.210115 ----
Sequences.LANL.preprocess.210115 <- function(dir.name= DATA)
{
	require(ape)
	require(gdata)
	require(data.table)
	
	dir.name		<- '/Users/alexb/Box Sync/Roadmap/Data/data_210114'
	infile <- file.path(dir.name,'LANL210218.txt')
	outfile <- file.path(dir.name,'data_210218_LANL_raw.fasta')
	
	ds <- read.dna(infile, format='fasta')

	# Takes a while to save
	write.dna(ds, file=outfile, format='fasta', colsep='', nbcol=-1)

	outfileNL			<- file.path(dir.name,"data_210218_LANL_noNL.fasta")
	
	# Remove dutch sequences
	ids <- rownames(ds)
	dutch <- grepl('.NL.',ids,fixed=1)
	df <- ds[!dutch,]
	df <- as.list(df)
	
	# Remove any gaps from end of taxa labels
	rownames(df) <- gsub('([ ]+)$','',rownames(df))

	# Takes a while to save
	write.dna(df, file=outfileNL, format='fasta', colsep='', nbcol=-1)
	
	# Remove DRMs
	lfile <- file.path(dir.name,'data_210218_LANL_noNL.fasta')
	rfile  <- file.path(dir.name,'K03455.fasta')
	loutfile <- file.path(dir.name,'data_210218_LANL_alignment_noDRM.fasta')
	seq					<- read.dna(lfile, format='fa')
	hxb2					<- read.dna(rfile, format='fa')
	# replace first sequence with hxb2 ref sequence
	seq <- rbind(hxb2,seq[-1,])
	tmp					<- which(grepl("K03455",rownames(seq)))
	rownames(seq)[tmp]	<- 'HXB2'
	tmp					<- seq.rm.drugresistance(seq)
	nodr.info		<- tmp$nodr.info
	seq					<- tmp$nodr.seq
	write.dna(seq, file= loutfile, format='fasta')
	save(seq, nodr.info, file= gsub('fasta','rda',loutfile))
	
}
###################

Sequences.CW.removeDRMs.210115<- function(dir.name= DATA)
{
	require(ape)
	require(big.phylo)
	dir.name		<- '~/Box Sync/Roadmap/Data/data_200821/'
	dir.name <- '/rds/general/project/ratmann_roadmap_data_analysis/live/alignments'
	file <- file.path(dir.name,'CW_seqs_aligned_trimmed.fasta')
	outfile	<- file.path(dir.name,'CW_seqs_aligned_noDRMs.fasta')
	
	seq	<- read.dna(file, format='fa')
	tmp	<- which(grepl("HXB2",rownames(seq)))
	rownames(seq)[tmp]	<- 'HXB2'
	tmp	<- big.phylo:::seq.rm.drugresistance(seq)
	nodr.info	<- tmp$nodr.info
	seq	<- tmp$nodr.seq
	write.dna(seq, file= outfile, format='fasta')
	save(seq, nodr.info, file= gsub('fasta','rda',outfile))
}

## ---- rmd.chunk.athena.Sequences.addLANL.210115 ----
Sequences.addLANL.210115<- function(dir.name= DATA)
{
	require(plyr)
	require(dplyr)
	require(tidyr)
	require(ape)
	require(data.table)
	
	home <- '~/Box Sync/Roadmap'
	dir.name		<- '~/Box Sync/Roadmap/Data/data_200821'
	file.c			<- file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_alignment_200908_noDRM_COMET.csv")
	file <- file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_alignment_200908_noDRM.fasta")
	
	#virulign '/Users/alexb/Box Sync/Roadmap/Data/K03455-pol.fasta' '/Users/alexb/Box Sync/Roadmap/Data/data_200821/CW_seqs.fasta'   --exportKind GlobalAlignment   --exportAlphabet Nucleotides --exportReferenceSequence yes --nt-debug Failed > '/Users/alexb/Box Sync/Roadmap/Data/data_200821/CW_seqs_aligned.fasta'

	dir.name		<- '~/Box Sync/Roadmap/Data/data_200821'
	bfile <- file.path(dir.name,'SHM_1902_ROADMAP_200821_tblLAB_alignment_200821_LANL_BLAST_20.txt')
	lfile <- file.path(dir.name,'data_210218_LANL_alignment_noDRM.fasta')
	afile <- file.path(dir.name,'SHM_1902_ROADMAP_200821_tblLAB_alignment_200821_noDRM.fasta')
	#afile <- file.path(dir.name,'ATHENA_210114_NL_Sequences.fasta')
	sfile <- file.path(dir.name,'Complete_subtype_classifications.rds')
	cfile <- file.path(dir.name,'CUseq.txt')
	infile.cw <- file.path(dir.name,'CW_seqs_aligned_noDRMs.fasta')
	infile.seqinfo			<- file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_seq.rda")
	infile.indinfo <- file.path(dir.name,'SHM_1902_ROADMAP_191223_tblBAS.csv')
	outfile.cw <- file.path(dir.name,'CUseq.fasta')
	outfile.cwst <- file.path(dir.name,'CW_seqs.fasta')
	
	dir.name <- '~/Box Sync/Roadmap/Data/data_210114'
	lfile <- file.path(dir.name,'data_210218_LANL_alignment_noDRM.fasta')
	sfile <- file.path(dir.name,'Complete_subtype_classifications.rds')
	outfile.matches      <- file.path(dir.name,'ATHENA_210114_All_Taxa_withsubtype.rda')
	outfile.lanl <- file.path(dir.name,'ATHENA_210218_LANL_alignment_noDRM_withCW.fasta')
	outfile.nl.seqs <- file.path(dir.name,'ATHENA_210114_NL_Sequences.fasta')
	outfile.seqs <- file.path(dir.name,'ATHENA_210218_All_Sequences.fasta')
	outfile.trimseqs <- file.path(dir.name,'ATHENA_210218_All_Sequences_LANL_aligned.fasta')
	
	## add Curacao seqs
	cw.dna <- read.delim(cfile, header = TRUE, sep = "\t", dec = ".")
	cw.dna <- as.data.table(cw.dna)
	cw.dna$DATE <- as.Date(cw.dna$DateRes,format="%d/%m/%Y")
	cw.dna$YEAR <- as.numeric(format(cw.dna$DATE, "%Y"))
	cw.dna$PID <- seq(1,nrow(cw.dna),1)
	cw.dna$lab <- paste0("CW.",cw.dna$YEAR,".",cw.dna$PID,".-")
	cw.dna$DateRes <- NULL
	cw.dna$DATE <- NULL
	cw.dna$YEAR <- NULL
	cw.dna$PID <- NULL
	cw.dna <- cw.dna[, c(2,1)]
	cw.dna$sequence <- as.character(cw.dna$sequence)
	cw <- as.DNAbin(strsplit(cw.dna$sequence,''))
	names(cw)	<- cw.dna[,lab]
	write.dna(cw, file=outfile.cw, format='fasta', colsep='', nbcol=-1)
	
	csfile <- file.path(dir.name,'CUseq_COMET.csv')
	cw.st <- read.csv(csfile, header = TRUE)
	cw.dna <- merge(cw.dna,cw.st,by.x='lab',by.y='name')
	cw.dna <- cw.dna[!grepl("unassigned",cw.dna$lab),]
	cw.dna$lab <- paste0(cw.dna$subtype,".",cw.dna$lab)
	cw.dna$subtype <- NULL
	cw.dna$bootstrap.support <- NULL
		
	cw <- as.DNAbin(strsplit(cw.dna$sequence,''))
	names(cw)	<- cw.dna[,lab]
	write.dna(cw,file=outfile.cwst, format='fasta', colsep='', nbcol=-1) 
	# add hxb2
	# mafft --reorder --add new_sequences --auto input 
	
	st <- readRDS(sfile)
	dl <- read.dna(lfile,format='fa')
	db <- data.table(TAXA_L=rownames(dl))
	
	# Get subtype of each LANL match
	db <- as.data.table(db)
	db[, SUBTYPE_L := substr(TAXA_L,1,regexpr("\\.",TAXA_L)-1)]
	db[TAXA_L=="HXB2", SUBTYPE_L:='B']
	
	# Read in Dutch seqs
	dbs			<- read.dna(lfile, format='fa')
	aseq <- read.dna(afile,format='fa')
	
	rownames(dbs) <- gsub('([ ]+)$','',rownames(dbs))
	hxb2lab <- rownames(seqs.all)[grepl('K03455',rownames(seqs.all))]
	# Rename HXB2 with Taxa label for matching
	rownames(dbs)[rownames(dbs)=='HXB2']	<- hxb2lab

	# Make sure HXB2 is in the file
	ids <- rownames(dbs)
	hxb2 <- grepl('K03455',ids,fixed=1)
	hxb2lab <- ids[hxb2]
	
	# Add Dutch sequences to LANL ones
	ids <- rownames(aseq)
	hxb2 <- grepl('HXB2',ids,fixed=1)
	rownames(aseq)[hxb2] <- hxb2lab
	
	# Only keep seqs longer than 750 nucleotides, first sequence per patient
	load(infile.seqinfo)
	dseq <- ds
	dseq <- as.data.table(dseq)
	setkey(dseq, PATIENT, SEQ_D)
	dseq <- subset(dseq,SEQ_L>=750)
	first	<- dseq[, list(SEQUENCE=SEQ_LABEL[1]), by='PATIENT']
	aseq <- aseq[labels(aseq) %in% c(hxb2lab,first$SEQUENCE),]

	# remove sequences with no HIV pos date and diagnosed after 31st Dec (lookup from metadata)
	infile.indinfo <- file.path(home,'Data','data_191223_Amsterdam','SHM_1902_ROADMAP_191223_tblBAS.csv')
	infile.nlinfo <- file.path(home,'Data','data_200316_Netherlands','SHM_1902_ROADMAP_200316_tblBAS.csv')
	dind <- read.csv(infile.indinfo,header=T)
	dnl <- read.csv(infile.nlinfo,header=T)
	dind[,'CITY'] <- 'Amsterdam'
	dnl[,'CITY'] <- 'Non-Amsterdam'
	dind <- merge(dind,dnl,all=T)
	dind[,'HIV1_POS_D'] <- as.Date(dind[,'HIV1_POS_D'],format="%d/%m/%Y")
	unknown_posd <- subset(dind, is.na(HIV1_POS_D) | HIV1_POS_D=="1911-11-11" | HIV1_POS_D>='2018-12-31')
	
	`%notin%` <- Negate(`%in%`)
	keep <- dseq[c(dseq$PATIENT %notin% unknown_posd$PATIENT),]
	aseq <- aseq[labels(aseq) %in%  c(hxb2lab,keep$SEQ_LABEL),]
	outfile.NL.final <- file.path(dir.name,'ATHENA_210114_NL_Sequences_final.fasta')
	write.dna(aseq,file=outfile.NL.final, format='fasta', colsep='', nbcol=-1) 
	
	# Add Curacao seqs to LANL
	cw.seqs <- read.dna(infile.cw,format='fa')
	ids <- rownames(cw.seqs)
	hxb2 <- grepl('HXB2',ids,fixed=1)
	rownames(cw.seqs)[hxb2] <- hxb2lab
	seqs.all <- seq.align.based.on.common.reference(dbs, cw.seqs, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
	write.dna(seqs.all,file=outfile.lanl, format='fasta', colsep='', nbcol=-1) 
	
	# exclude LANL sequences with no origin or sequence date
	seqs.all <- read.dna(outfile.lanl,format='fa')
	labs <- tibble(TAXA=unlist(rownames(seqs.all))) %>%
		mutate( SUBTYPE:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',1),
						Alpha.2.code:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',2),
						YEAR:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',3),
						GENBANK:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',5)) %>%
		distinct()
	labs <- subset(labs,YEAR!='-' & Alpha.2.code!='-')
	seqs.all <- seqs.all[labels(seqs.all) %in%  labs$TAXA,]
	
	# Add LANL seqs
	seqs.all <- seq.align.based.on.common.reference(aseq, seqs.all, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
	
	write.dna(seqs.all,file=outfile.seqs, format='fasta', colsep='', nbcol=-1) 
	seqs.all <- seqs.all[,1:1302]
	write.dna(seqs.all,file=outfile.trimseqs, format='fasta', colsep='', nbcol=-1) 
	
	# Read back in
	seqs.all <- read.dna(outfile.trimseqs,format='fa')

	# Read in outgroup reference sequences
	rfile <- file.path(dir.name,'HIV1_SUBTYPE_REF_2010_2253-3870_DNA.fasta')
	rs <- read.dna(rfile,format='fa')
	rownames(rs) <- gsub('Ref','Outgroup',rownames(rs))
	ids <- labels(rs)
	
	#	select closest sequences for subtype B
	outfile.b <- file.path(dir.name,'alignments','ROADMAP_210218_Sequences_LANL_aligned_noDRM_subtype_B_wOutgroup.fasta')
	db[1,TAXA_L:=hxb2lab]
	db.l			<- unique(subset(db, SUBTYPE_L=='B', TAXA_L))
	db.a			<- unique(subset(st, SUBTYPE=='B', FASTASampleCode))
	db.cw			<- labels(cw.seqs[grepl('B.CW.',labels(cw.seqs),0),])
	db.b <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],db.cw),]
	st.b <- ids[c(grep('K03455',ids,fixed=1),grep('.D.',ids,fixed=1))] 
	sub.b <- rs[labels(rs) %in% st.b,]
	scn				<- seq.align.based.on.common.reference(db.b, sub.b, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
	scn <- scn[,1:1302]
	# replace label names with NL__
	rownames(scn) <- gsub('Amst_','NL_',rownames(scn))
	write.dna(scn, file=outfile.b, format='fa', colsep='', nbcol=-1)
	
	# Subtype 02_AG
	outfile.02ag <- file.path(dir.name,'alignments','ROADMAP_210218_Sequences_LANL_aligned_noDRM_subtype_02AG_wOutgroup.fasta')
	db.l			<- unique(subset(db, SUBTYPE_L=='02_AG', TAXA_L))
	db.a			<- unique(subset(st, SUBTYPE=='02_AG', FASTASampleCode))
	db.cw			<- labels(cw.seqs[grepl('02_AG.CW.',labels(cw.seqs),0),])
	db.02ag <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],db.cw,hxb2lab),]
	st.02ag <- ids[c(grep('K03455',ids,fixed=1),grep('.D.',ids,fixed=1))] 
	sub.02ag <- rs[labels(rs) %in% st.02ag,]
	scn				<- seq.align.based.on.common.reference(db.02ag, sub.02ag, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
	scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
	scn <- scn[,1:1302]
	# replace label names with NL__
	rownames(scn) <- gsub('Amst_','NL_',rownames(scn))
	write.dna(scn, file=outfile.02ag, format='fa', colsep='', nbcol=-1)
	
	# Subtype C
	outfile.c <- file.path(dir.name,'alignments','ROADMAP_210218_Sequences_LANL_aligned_noDRM_subtype_C_wOutgroup.fasta')
	db.l			<- unique(subset(db, SUBTYPE_L=='C', TAXA_L))
	db.a			<- unique(subset(st, SUBTYPE=='C', FASTASampleCode))
	db.cw			<- labels(cw.seqs[grepl('C.CW.',labels(cw.seqs),0),])
	db.c <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],db.cw,hxb2lab),]
	st.c <- ids[c(grep('K03455',ids,fixed=1),grep('.D.',ids,fixed=1))] 
	sub.c <- rs[labels(rs) %in% st.c,]
	scn				<- seq.align.based.on.common.reference(db.c, sub.c, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
	scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
	scn <- scn[,1:1302]
	# replace label names with NL__
	rownames(scn) <- gsub('Amst_','NL_',rownames(scn))
	write.dna(scn, file=outfile.c, format='fa', colsep='', nbcol=-1)
	
	# Subtype A1
	outfile.a1 <- file.path(dir.name,'alignments','ROADMAP_210218_Sequences_LANL_aligned_noDRM_subtype_A1_wOutgroup.fasta')
	db.l			<- unique(subset(db, SUBTYPE_L %in% c('A1','A'), TAXA_L))
	db.a			<- unique(subset(st, SUBTYPE=='A1', FASTASampleCode))
	db.a1 <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
	st.a1 <- ids[c(grep('K03455',ids,fixed=1),grep('.C.',ids,fixed=1))] 
	sub.a1 <- rs[labels(rs) %in% st.a1,]
	scn				<- seq.align.based.on.common.reference(db.a1, sub.a1, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
	scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
	scn <- scn[,1:1302]
	# replace label names with NL__
	rownames(scn) <- gsub('Amst_','NL_',rownames(scn))
	write.dna(scn, file=outfile.a1, format='fa', colsep='', nbcol=-1)
	
	# Subtype 01_AE
	outfile.01ae <- file.path(dir.name,'alignments','ROADMAP_210218_Sequences_LANL_aligned_noDRM_subtype_01AE_wOutgroup.fasta')
	db.l			<- unique(subset(db, SUBTYPE_L=='01_AE', TAXA_L))
	db.a			<- unique(subset(st, SUBTYPE=='01_AE', FASTASampleCode))
	db.01ae <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
	st.01ae <- ids[c(grep('K03455',ids,fixed=1),grep('.D.',ids,fixed=1))] 
	sub.01ae <- rs[labels(rs) %in% st.a1,]
	scn				<- seq.align.based.on.common.reference(db.01ae, sub.01ae, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
	scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
	scn <- scn[,1:1302]
	# replace label names with NL__
	rownames(scn) <- gsub('Amst_','NL_',rownames(scn))
	write.dna(scn, file=outfile.01ae, format='fa', colsep='', nbcol=-1)
	
	# Subtype G
	outfile.g <- file.path(dir.name,'alignments','ROADMAP_210218_Sequences_LANL_aligned_noDRM_subtype_G_wOutgroup.fasta')
	db.l			<- unique(subset(db, SUBTYPE_L=='G', TAXA_L))
	db.a			<- unique(subset(st, SUBTYPE=='G', FASTASampleCode))
	db.g <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
	st.g <- ids[c(grep('K03455',ids,fixed=1),grep('.D.',ids,fixed=1))] 
	sub.g <- rs[labels(rs) %in% st.g,]
	scn				<- seq.align.based.on.common.reference(db.g, sub.g, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
	scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
	scn <- scn[,1:1302]
	# replace label names with NL__
	rownames(scn) <- gsub('Amst_','NL_',rownames(scn))
	write.dna(scn, file=outfile.g, format='fa', colsep='', nbcol=-1)
	
	# Subtype D
	outfile.d <- file.path(dir.name,'alignments','ROADMAP_210218_Sequences_LANL_aligned_noDRM_subtype_D_wOutgroup.fasta')
	db.l			<- unique(subset(db, SUBTYPE_L=='D', TAXA_L))
	db.a			<- unique(subset(st, SUBTYPE=='D' , FASTASampleCode))
	db.d <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
	st.d <- ids[c(grep('K03455',ids,fixed=1),grep('.C.',ids,fixed=1))] 
	sub.d <- rs[labels(rs) %in% st.d,]
	scn				<- seq.align.based.on.common.reference(db.d, sub.d, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
	scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
	scn <- scn[,1:1302]
	# replace label names with NL__
	rownames(scn) <- gsub('Amst_','NL_',rownames(scn))
	write.dna(scn, file=outfile.d, format='fa', colsep='', nbcol=-1)
	
	# Subtype 06_cpx (check if count>50)
	outfile.06cpx <- file.path(dir.name,'alignments','ROADMAP_210218_Sequences_LANL_aligned_noDRM_subtype_06cpx_wOutgroup.fasta')
	db.l			<- unique(subset(db, SUBTYPE_L=='06_cpx', TAXA_L))
	db.a			<- unique(subset(st, SUBTYPE=='06_cpx', FASTASampleCode))
	db.06cpx <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
	st.06cpx <- ids[c(grep('K03455',ids,fixed=1),grep('.D.',ids,fixed=1))]
	sub.06cpx <- rs[labels(rs) %in% st.06cpx,]
	scn				<- seq.align.based.on.common.reference(db.06cpx, sub.06cpx, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
	scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
	scn <- scn[,1:1302]
	# replace label names with NL__
	rownames(scn) <- gsub('Amst_','NL_',rownames(scn))
	write.dna(scn, file=outfile.06cpx, format='fa', colsep='', nbcol=-1)
	
	# Subtype F
	outfile.F <- file.path(dir.name,'alignments','ROADMAP_210218_Sequences_LANL_aligned_noDRM_subtype_F1_wOutgroup.fasta')
	db.l			<- unique(subset(db, SUBTYPE_L %in% c('F1','F'), TAXA_L))
	db.a			<- unique(subset(st, SUBTYPE=='F1', FASTASampleCode))
	db.F <- seqs.all[labels(seqs.all) %in% c(db.l[, TAXA_L],db.a[,FASTASampleCode],hxb2lab),]
	st.F <- ids[c(grep('K03455',ids,fixed=1),grep('.A1.',ids,fixed=1))]
	sub.F <- rs[labels(rs) %in% st.F,]
	scn				<- seq.align.based.on.common.reference(db.F, sub.F, return.common.sites=FALSE, regexpr.reference='K03455', regexpr.nomatch='-|\\?')
	scn       <- scn[!(grepl('K03455',labels(scn),fixed=1)),]
	scn <- scn[,1:1302]
	# replace label names with NL__
	rownames(scn) <- gsub('Amst_','NL_',rownames(scn))
	write.dna(scn, file=outfile.F, format='fa', colsep='', nbcol=-1)
	
}


## ---- rmd.chunk.athena.220115.alignment.make.bootstraps ----
athena.220115.alignment.make.bootstraps <- function(analysis)
{
  require(ape)
  require(data.table)
  
	analysis <- 'analysis_210504_ATHENA'
	
  home <- '/Users/alexb/Box Sync/Roadmap'
  #home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
  indir <- file.path(home,analysis,'alignments_subtype')
  outdir <- file.path(home,analysis,'alignments_bs')
  bsn <- 100
  seed <- 42L
  #
  df <- data.table(F=list.files(indir, pattern='fasta$',full.names=TRUE))
  df <- subset(df, grepl('noDRM',F))
  for(i in seq_len(nrow(df)))
  {
    infile.fasta <- df[i,F]
    cat('\nprocess ',infile.fasta)
    sq                  <- read.dna(infile.fasta, format='fa')
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

## ---- rmd.chunk.athena.220115.splitbalignment ----
athena.220115.splitbalignment<- function(analysis)
{
	# after running phyloscanner split B alignment into clades and re-run fasttree
	require(big.phylo)
	require(data.table)
	
	analysis <- 'analysis_210504_ATHENA'
	
	home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	
	indir.trees <- file.path(home,analysis,'phyloscanner_B')
	indir.aln <- file.path(home,'analysis_210226_ATHENA','alignments_bs')
	
	infile.labels <- data.table(FLABEL= c(  file.path(home,analysis,'misc','200917_sequence_labels_NLMSM.csv'),
																					file.path(home,analysis,'misc','200917_sequence_labels_NLHSX.csv'))
	)
	outdir <- file.path(home,analysis,'phyloscanner')
	infiles <- data.table(FIN=list.files(indir.trees, pattern='\\.newick$',full.names=TRUE))
	infile.labels[, SELECT:= gsub('.*_labels_([A-Z]+).*','\\1',basename(FLABEL))]
	infiles[, ST:= gsub('.*_subtype_([A-Za-z0-9]+)_.*','\\1',basename(FIN))]
	infiles[, REP:= gsub('.*_wOutgroup_([A-Za-z0-9]+)_.*','\\1',basename(FIN))]
	infiles[, SELECT:= gsub('.*_rerooted_([A-Z]+).*','\\1',basename(FIN))]
	infiles <- merge(infiles, infile.labels, by='SELECT', allow.cartesian=TRUE)

	indir <- file.path(home,analysis,'alignments_bs')
	outdir <- file.path(home,analysis,'fasttree_B')
	
	#   get UNIX commands for all trees to be processed
	df <- data.table(FIN=list.files(indir.aln, pattern='fasta$',full.names=TRUE))
	df <- subset(df, grepl('noDRM',FIN))
	df[, ST:= gsub('.*_subtype_([A-Z0-9a-z]+)_.*','\\1',basename(FIN))]
	df[, FOUT:= file.path(outdir,gsub('\\.fasta','_ft.newick',basename(FIN)))]
	df[, REP:= gsub('.*_wOutgroup_([A-Za-z0-9]+).*','\\1',basename(FIN))]
	df <- subset(df, ST=='B')
	
	set.seed(42L)
	for(i in seq_len(nrow(infiles)))
		{
			#i<- 1
			cat('\nprocess ',infiles[i,FIN],'\nprocess ',infiles[i,FLABEL])
			infile <- infiles[i,FIN]
			rep <- infiles[i,REP]
			ph <- read.tree(infile)
			ph$node.label <- NULL
			#   update tip labels: add world region to start of label
			local.world.region <- gsub('^.*_labels_([A-Za-z]+)\\.csv$','\\1', basename(infiles[i,FLABEL]))
			dl <- as.data.table(read.csv(infiles[i,FLABEL], stringsAsFactors=FALSE))
			dl[, X:=NULL]
			dl <- dl[dl$TAXA_NEW %in% ph$tip.label,]
			infile.fasta <- df[REP==rep,FIN]
			sq <- read.dna(infile.fasta, format='fa')
			sq <- sq[labels(sq) %in% c(dl[, TAXA],labels(sq[grep('Outgroup',labels(sq),fixed=1),])),]
			outfile <- gsub('_B_',paste0('_',infiles[i,ST],'_'),infile.fasta)
			outfile <- gsub('_ft',paste0('_ft_',infiles[i,SELECT]),outfile)
			outfile <- file.path(outdir,basename(outfile))
			write.dna(sq, file=outfile, format='fa', colsep='', nbcol=-1)
		}
}

## ---- rmd.chunk.athena.220115.fastree ----
athena.220115.fastree<- function(analysis)
{
  require(big.phylo)
  require(data.table)

	analysis <- 'analysis_210504_ATHENA'
	
  home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'

  indir <- file.path(home,analysis,'alignments_bs')
  outdir <- file.path(home,analysis,'fasttree')
  
  #   get UNIX commands for all trees to be processed
  df <- data.table(FIN=list.files(indir, pattern='fasta$',full.names=TRUE))
  df <- subset(df, grepl('noDRM',FIN))
  df[, ST:= gsub('.*_subtype_([A-Z0-9a-z]+)_.*','\\1',basename(FIN))]
  df[, FOUT:= file.path(outdir,gsub('\\.fasta','_ft.newick',basename(FIN)))]

  dt <- data.table(FT=list.files(outdir, pattern='newick$',full.names=TRUE))
  dt[, IN:=1]
  
  df <- merge(df,dt,by.x=c('FOUT'),by.y=('FT'),all=T)
  
  df <- subset(df, is.na(IN))
  
  cmds <- vector('list', nrow(df))
  for(i in seq_len(nrow(df)))
  {
    infile.fasta <- df[i,FIN]
    outfile <- df[i,FOUT]
    tmp <- cmd.fasttree(infile.fasta, outfile=outfile, pr.args='-nt -gtr -gamma', check.binary=TRUE)
    cmds[[i]] <- tmp
  }
  df[, CMD:= unlist(cmds)]
  
  #   submit jobs like this one:
  cat(df[1,CMD])
  
  #   run on HPC as array job
  df[, CASE_ID:= 1:nrow(df)]
  
  #   make PBS header
  #hpc.load    <- "module load R/3.4.0"
  hpc.load    <- "module load anaconda3/personal"
  r.activate  <- "source activate Renv3"
  hpc.select  <- 1                        # number of nodes
  hpc.nproc   <- 1                        # number of processors on node
  hpc.walltime<- 72                      # walltime
  hpc.q       <- NA                       # PBS queue
  hpc.mem     <- "6gb"                    # RAM
  hpc.array   <- length(unique(df$CASE_ID))   # number of runs for job array
  pbshead     <- "#!/bin/sh"
  tmp         <- paste("#PBS -l walltime=", hpc.walltime, ":00:00,pcput=", hpc.walltime, ":00:00", sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  tmp         <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  pbshead     <- paste(pbshead, "#PBS -j oe", sep = "\n")
  if(!is.na(hpc.array))
    pbshead <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
  if(!is.na(hpc.q))
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead     <- paste(pbshead, hpc.load, r.activate, sep = "\n")
  #   make array job
  cmd     <- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
  cmd     <- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]
  cmd     <- paste(pbshead,cmd,sep='\n')
  #   submit job
  outfile     <- gsub(':','',paste("trs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
  outfile     <- file.path(outdir, outfile)
  cat(cmd, file=outfile)
  cmd         <- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))
}


## ---- rmd.chunk.athena.220115.sequence.labels ----
athena.220115.sequence.labels<- function(analysis)
{
	require(data.table)
  require(ape)
  require(tidyverse)
	require(big.phylo)
  
	analysis <- 'analysis_210504_ATHENA'
	source('/rds/general/project/ratmann_roadmap_data_analysis/live/R/bigphylo functions.R')
  home <- '/Users/alexb/Box Sync/Roadmap'
  home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'

    infile.indinfo <- file.path(home,'Data','data_191223_Amsterdam','SHM_1902_ROADMAP_191223_tblBAS.csv')
    infile.nlinfo <- file.path(home,'Data','data_200316_Netherlands','SHM_1902_ROADMAP_200316_tblBAS.csv')
  	infile.subtypes <- file.path(home,'Data','data_200821','ROADMAP_200821_All_Taxa_withsubtype.rda')
  	infile.subtypes.NL <- file.path(home,analysis,'misc','ATHENA_210114_withsubtype.rda')
  	infile.subtypes.LANL <- file.path(home,analysis,'misc','ATHENA_210114_LANL_withsubtype.rda')
    infile.seqinfo <- file.path(home,'Data','data_200821','SHM_1902_ROADMAP_200821_tblLAB_seq.rda')
    infile.georeg <- file.path(home,'misc','UN_geographic_locations.csv')
    infile.georegnew <- file.path(home,'misc','NEWGEO.csv')
    infile.countrycodes <- file.path(home,'misc','ISO_country_codes.csv')
    infiles.lanl <- file.path(home,analysis,'alignments','ATHENA_210218_LANL_alignment_noDRM_withCW.fasta')
    outfile.noorigin <- file.path(home,analysis,'misc','No_birth_country.csv')
    outfile.demodata <- file.path(home,analysis,'misc','Baseline_data_Ams_NL.csv')
    outfile.base <- file.path(home,analysis,'misc','200917_')
    
    
  dgeo <- read.csv(infile.georegnew,header=T)
  dgeo <- dgeo %>% select(c('Alpha_2_code', 'Alpha_3_code','CNTRY','WRLD')) %>%
  					rename('Alpha.2.code'='Alpha_2_code') %>% rename('Alpha.3.code'='Alpha_3_code') %>%
  	mutate(Alpha.2.code:=as.character(Alpha.2.code),
  				 WRLD:= as.character(WRLD),
  				 CNTRY:= as.character(CNTRY))
    
  #
  # read Amsterdam individual data and merge with Dutch data
  #
  dind <- read.csv(infile.indinfo,header=T)
  dnl <- read.csv(infile.nlinfo,header=T)
  dind[,'CITY'] <- 'Amsterdam'
  dnl[,'CITY'] <- 'Non-Amsterdam'
  dind <- merge(dind,dnl,all=T)
  
# Format dates to numeric and add ranges for ones which may be inaccurate
  dind[,'BIRTH_D'] <- as.Date(dind[,'BIRTH_D'])
  dind[,'MIG_D'] <- as.Date(dind[,'MIG_D'],format="%d/%m/%Y")
  dind[,'MIG_D_pq'] <- as.Date(dind[,'MIG_D_pq'],format="%d/%m/%Y")
  dind[,'MIG_D_aq'] <- as.Date(dind[,'MIG_D_aq'],format="%d/%m/%Y")
  dind$MIG_D_lower <- hivc.db.Date.lower(dind$MIG_D)
  dind$MIG_D_upper <- hivc.db.Date.upper(dind$MIG_D)
  dind$MIG_D_lower[!is.na(as.POSIXlt(dind$MIG_D_pq))] <- dind$MIG_D_pq[!is.na(as.POSIXlt(dind$MIG_D_pq))]
  dind$MIG_D_upper[!is.na(as.POSIXlt(dind$MIG_D_aq))] <- dind$MIG_D_aq[!is.na(as.POSIXlt(dind$MIG_D_aq))]
  dind$MIG_D_lower <- as.Date(dind$MIG_D_lower)
  dind$MIG_D_upper <- as.Date(dind$MIG_D_upper)
  
  dind[,'HIV1_NEG_D'] <- as.Date(dind[,'HIV1_NEG_D'],format="%d/%m/%Y")
  dind[,'HIV1_NEG_D_pq'] <- as.Date(dind[,'HIV1_NEG_D_pq'],format="%d/%m/%Y")
  dind[,'HIV1_NEG_D_aq'] <- as.Date(dind[,'HIV1_NEG_D_aq'],format="%d/%m/%Y")
  dind$HIV1_NEG_D_lower <- hivc.db.Date.lower(dind$HIV1_NEG_D)
  dind$HIV1_NEG_D_upper <- hivc.db.Date.upper(dind$HIV1_NEG_D)
  dind$HIV1_NEG_D_lower[!is.na(as.POSIXlt(dind$HIV1_NEG_D_pq))] <- dind$HIV1_NEG_D_pq[!is.na(as.POSIXlt(dind$HIV1_NEG_D_pq))]
  dind$HIV1_NEG_D_upper[!is.na(as.POSIXlt(dind$HIV1_NEG_D_aq))] <- dind$HIV1_NEG_D_aq[!is.na(as.POSIXlt(dind$HIV1_NEG_D_aq))]
  dind$HIV1_NEG_D_lower <- as.Date(dind$HIV1_NEG_D_lower)
  dind$HIV1_NEG_D_upper <- as.Date(dind$HIV1_NEG_D_upper)
  
  dind[,'HIV1_POS_D'] <- as.Date(dind[,'HIV1_POS_D'],format="%d/%m/%Y")
  dind[,'HIV1_POS_D_pq'] <- as.Date(dind[,'HIV1_POS_D_pq'],format="%d/%m/%Y")
  dind[,'HIV1_POS_D_aq'] <- as.Date(dind[,'HIV1_POS_D_aq'],format="%d/%m/%Y")
  dind$HIV1_POS_D_lower <- hivc.db.Date.lower(dind$HIV1_POS_D)
  dind$HIV1_POS_D_upper <- hivc.db.Date.upper(dind$HIV1_POS_D)
  dind$HIV1_POS_D_lower[!is.na(as.POSIXlt(dind$HIV1_POS_D_pq))] <- dind$HIV1_POS_D_pq[!is.na(as.POSIXlt(dind$HIV1_POS_D_pq))]
  dind$HIV1_POS_D_upper[!is.na(dind$HIV1_POS_D_aq)] <- dind$HIV1_POS_D_aq[!is.na(dind$HIV1_POS_D_aq)]
  dind$HIV1_POS_D_lower <- as.Date(dind$HIV1_POS_D_lower)
  dind$HIV1_POS_D_upper <- as.Date(dind$HIV1_POS_D_upper)
  
  dind[,'T0'] <- as.Date(dind[,'T0'],format="%d/%m/%Y")
  dind$T0_lower <- as.Date(hivc.db.Date.lower(dind$T0))
  dind$T0_upper <- as.Date(hivc.db.Date.upper(dind$T0))
  
  dind[,'RECART_D'] <- as.Date(dind[,'RECART_D'],format="%Y-%m-%d")
  dind$RECART_D_lower <- as.Date(hivc.db.Date.lower(dind$RECART_D))
  dind$RECART_D_upper <- as.Date(hivc.db.Date.upper(dind$RECART_D))
  
  dind[,'CARE_FRS_D'] <- as.Date(dind[,'CARE_FRS_D'],format="%d/%m/%Y")
  dind$CARE_FRS_D_lower <- as.Date(hivc.db.Date.lower(dind$CARE_FRS_D))
  dind$CARE_FRS_D_upper <- as.Date(hivc.db.Date.upper(dind$CARE_FRS_D))
  dind[,'DROP_D'] <- as.Date(dind[,'DROP_D'],format="%Y-%m-%d")
  dind$DROP_D_lower <- as.Date(hivc.db.Date.lower(dind$DROP_D))
  dind$DROP_D_upper <- as.Date(hivc.db.Date.upper(dind$DROP_D))
  
  dind[,'BIRTH_D'] <- hivc.db.Date2numeric(dind[,'BIRTH_D'])
  dind[,'BIRTH_Y'] <- floor(dind[,'BIRTH_D'])
  dind[,'MIG_D'] <- hivc.db.Date2numeric(dind$MIG_D)
  dind[,'MIG_D_lower'] <- hivc.db.Date2numeric(dind$MIG_D_lower)
  dind[,'MIG_D_upper'] <- hivc.db.Date2numeric(dind$MIG_D_upper)
  dind[,'HIV1_NEG_D'] <- hivc.db.Date2numeric(dind[,'HIV1_NEG_D'])
  dind[,'HIV1_NEG_D_lower'] <- hivc.db.Date2numeric(dind$HIV1_NEG_D_lower)
  dind[,'HIV1_NEG_D_upper'] <- hivc.db.Date2numeric(dind$HIV1_NEG_D_upper)
  dind[,'HIV1_POS_D'] <- hivc.db.Date2numeric(dind[,'HIV1_POS_D'])
  dind[,'HIV1_POS_D_lower'] <- hivc.db.Date2numeric(dind$HIV1_POS_D_lower)
  dind[,'HIV1_POS_D_upper'] <- hivc.db.Date2numeric(dind$HIV1_POS_D_upper)
  dind[,'T0'] <- hivc.db.Date2numeric(dind[,'T0'])
  dind[,'T0_lower'] <- hivc.db.Date2numeric(dind[,'T0_lower'])
  dind[,'T0_upper'] <- hivc.db.Date2numeric(dind[,'T0_upper'])
  dind[,'RECART_D'] <- hivc.db.Date2numeric(dind[,'RECART_D'])
  dind[,'RECART_D_lower'] <- hivc.db.Date2numeric(dind[,'RECART_D_lower'])
  dind[,'RECART_D_upper'] <- hivc.db.Date2numeric(dind[,'RECART_D_upper'])
  dind[,'CARE_FRS_D_lower'] <- hivc.db.Date2numeric(dind[,'CARE_FRS_D_lower'])
  dind[,'CARE_FRS_D_upper'] <- hivc.db.Date2numeric(dind[,'CARE_FRS_D_upper'])
  dind[,'DROP_D'] <-  hivc.db.Date2numeric(dind[,'DROP_D'])
  dind[,'DROP_D_lower'] <-  hivc.db.Date2numeric(dind[,'DROP_D_lower'])
  dind[,'DROP_D_upper'] <-  hivc.db.Date2numeric(dind[,'DROP_D_upper'])
  
  # Recode baseline variables
    dind <- as_tibble(dind) %>%
      select(PATIENT,BIRTH_Y,GENDER,MODE,MODE_OTH,ORIGIN,MIG_D,MIG_D_lower,MIG_D_upper,HIV1_NEG_D,HIV1_NEG_D_lower,HIV1_NEG_D_upper,
             HIV1_NEG_D_A,HIV1_POS_D,HIV1_POS_D_lower,HIV1_POS_D_upper,HIV1_POS_D_A,INF_NL,INF_COUNTRY_1,
             REG_FRS_GGD,REG_LST_GGD,CITY,RECART_Y,T0,T0_lower,T0_upper,RECART_D,RECART_D_lower,RECART_D_upper,NAIVE_Y,
             CARE_FRS_D,CARE_FRS_D_lower,CARE_FRS_D_upper,DROP_D,DROP_D_lower,DROP_D_upper,DROP_RS) %>%
      mutate( GENDER:= case_when(GENDER==1~'Male',
                                 GENDER==2~'Female'),
              TRANSM:= case_when(MODE==1~'MSM',
                                MODE==2~'IDU',
                                MODE==4~'Other',
                                MODE==5~'Other',
                                MODE==6~'HSX',
                                MODE==8~'Other',
                                MODE==9~'Other',
                                MODE==90~'Other',
                                MODE==99~'Unknown'),
              DUTCH_BORN:= case_when(ORIGIN=='Unknown'~'Unknown',
                                     ORIGIN==''~'Unknown',
                                     ORIGIN=="NL"~'Dutch',
                                     TRUE~'Other'),
              INF_NL:= case_when(INF_COUNTRY_1=='NL'~'Dutch',
                                 INF_COUNTRY_1=='-1'~'Unknown',
                                 INF_COUNTRY_1==""~'Unknown',
                                 INF_COUNTRY_1!='NL'&INF_COUNTRY_1!='-1'&INF_COUNTRY_1!=""~'Non-NL'),
              GGD_FIRST:= case_when(REG_FRS_GGD==111~'Groningen',
                                    REG_FRS_GGD==706~'Drenthe',
                                    REG_FRS_GGD==1009~'IJsselland',
                                    REG_FRS_GGD==1106~'Twente',
                                    REG_FRS_GGD==1406~'Gelre_IJssel',
                                    REG_FRS_GGD==1906~'Hulpverlening_Gelderland_Midden',
                                    REG_FRS_GGD==2006~'Rivierenland',
                                    REG_FRS_GGD==2106~'Nijmegen',
                                    REG_FRS_GGD==2209~'Flevoland',
                                    REG_FRS_GGD==2406~'Utrecht',
                                    REG_FRS_GGD==2506~'Midden_Nederland',
                                    REG_FRS_GGD==2707~'Hollands_Noorden',
                                    REG_FRS_GGD==3109~'Kennemerland',
                                    REG_FRS_GGD==3406~'Amsterdam',
                                    REG_FRS_GGD==3606~'Gooi_Vechtstreek',
                                    REG_FRS_GGD==3906~'Den_Haag',
                                    REG_FRS_GGD==4106~'Zuid_Holland_West',
                                    REG_FRS_GGD==4506~'Hollands_Midden',
                                    REG_FRS_GGD==4607~'Rotterdam_Rijnmond',
                                    REG_FRS_GGD==4810~'Zuid_Holland_Zuid',
                                    REG_FRS_GGD==5006~'Zeeland',
                                    REG_FRS_GGD==5206~'West_Brabant',
                                    REG_FRS_GGD==5406~'Hart_voor_Brabant',
                                    REG_FRS_GGD==5608~'Brabant_Zuidoost',
                                    REG_FRS_GGD==6011~'Limburg-Noord',
                                    REG_FRS_GGD==6106~'Zuid_Limburg',
                                    REG_FRS_GGD==7206~'Fryslan',
                                    REG_FRS_GGD==7306~'Zaanstreek_Waterland'),
              GGD_LAST:= case_when(REG_LST_GGD==111~'Groningen',
                                    REG_LST_GGD==706~'Drenthe',
                                    REG_LST_GGD==1009~'IJsselland',
                                    REG_LST_GGD==1106~'Twente',
                                    REG_LST_GGD==1406~'Gelre_IJssel',
                                    REG_LST_GGD==1906~'Hulpverlening_Gelderland_Midden',
                                    REG_LST_GGD==2006~'Rivierenland',
                                    REG_LST_GGD==2106~'Nijmegen',
                                    REG_LST_GGD==2209~'Flevoland',
                                    REG_LST_GGD==2406~'Utrecht',
                                    REG_LST_GGD==2506~'Midden_Nederland',
                                    REG_LST_GGD==2707~'Hollands_Noorden',
                                    REG_LST_GGD==3109~'Kennemerland',
                                    REG_LST_GGD==3406~'Amsterdam',
                                    REG_LST_GGD==3606~'Gooi_Vechtstreek',
                                    REG_LST_GGD==3906~'Den_Haag',
                                    REG_LST_GGD==4106~'Zuid_Holland_West',
                                    REG_LST_GGD==4506~'Hollands_Midden',
                                    REG_LST_GGD==4607~'Rotterdam_Rijnmond',
                                    REG_LST_GGD==4810~'Zuid_Holland_Zuid',
                                    REG_LST_GGD==5006~'Zeeland',
                                    REG_LST_GGD==5206~'West_Brabant',
                                    REG_LST_GGD==5406~'Hart_voor_Brabant',
                                    REG_LST_GGD==5608~'Brabant_Zuidoost',
                                    REG_LST_GGD==6011~'Limburg-Noord',
                                    REG_LST_GGD==6106~'Zuid_Limburg',
                                    REG_LST_GGD==7206~'Fryslan',
                                    REG_LST_GGD==7306~'Zaanstreek_Waterland'
                                 )) %>%
      mutate(REGION_FIRST:= case_when(GGD_FIRST=='Amsterdam'~'Amsterdam',
                                      GGD_FIRST=='Rotterdam'~'Rotterdam',
                                      GGD_FIRST=='The Hague'~'The Hague',
                                      GGD_FIRST=='Utrecht'~'Utrecht',
                                      TRUE~'Other'),
            REGION_LAST:= case_when(GGD_LAST=='Amsterdam'~'Amsterdam',
                                      GGD_LAST=='Rotterdam'~'Rotterdam',
                                      GGD_LAST=='The Hague'~'The Hague',
                                      GGD_LAST=='Utrecht'~'Utrecht',
                                      TRUE~'Other'
                                      )) %>%
      mutate( Alpha.2.code:= as.character(ORIGIN),
              ORIGIN:= as.character(ORIGIN)) %>%
      mutate(ORIGIN:= case_when(ORIGIN==""~'Unknown',
                                is.na(ORIGIN)~'Unknown',
                         TRUE ~ ORIGIN))
    
    # Obtain world region for country of origin

    dind <- dind %>% left_join(dgeo,by='Alpha.2.code')   %>%
            mutate(BIRTH_CNTRY:= case_when(CNTRY=='Unknown'~'Unknown',
                                           is.na(CNTRY)~'Unknown',
                                           TRUE ~ CNTRY),
                   LOC_BIRTH:= case_when(CNTRY=='Unknown'~'Unknown',
                                         is.na(CNTRY)~'Unknown',
                                         TRUE ~ WRLD))   %>%
            select(-WRLD,-CNTRY,-Alpha.2.code,-Alpha.3.code) %>%
            mutate(Alpha.2.code:= as.character(INF_COUNTRY_1)) %>%
            left_join(dgeo,by='Alpha.2.code') %>%
            mutate(INF_CNTRY:= case_when(CNTRY=='Unknown'~'Unknown',
                                     is.na(CNTRY)~'Unknown',
                                     TRUE ~ CNTRY),
                   LOC_INF:= case_when(CNTRY=='Unknown'~'Unknown',
                                   is.na(CNTRY)~'Unknown',
                                   TRUE ~ WRLD))

  # Export individuals without a birth country recorded for investigation
  noorigin <- dind[is.na(dind$ORIGIN) | dind$ORIGIN=="" | dind$ORIGIN=="Unknown",]
  write.csv(noorigin,file=outfile.noorigin)
  
  # Save demographic data
  write.csv(dind,file=outfile.demodata)
  
  # Read in subtype file
  load(infile.subtypes.NL)
  load(infile.subtypes.LANL)
  st.l <- as_tibble(db) %>% select(TAXA_L,SUBTYPE_L) %>% distinct()
  st.nl <- as_tibble(st) %>% select(FASTASampleCode,SUBTYPE) %>% distinct()
  colnames(st.nl)[1] <- 'SEQ_LABEL'
  
  #
  # read sequence labels for Amsterdam and NL and add in subtype and geographical data
  #
  load(infile.seqinfo)
  dseq <- ds
  dseq <- as_tibble(dseq) %>% inner_join(st.nl, by='SEQ_LABEL')
  
  dseq <- dind %>% inner_join(dseq, by='PATIENT')
  dseq <- dseq %>% select(PATIENT, SEQ_LABEL, SEQ_ID, SEQ_D, SEQ_L, SUBTYPE, ORIGIN, GENDER, TRANSM, CITY, BIRTH_CNTRY, LOC_BIRTH, INF_CNTRY, LOC_INF) %>%
          rename(SEQ_DATE=SEQ_D, SEQ_LENGTH=SEQ_L)

  #
  # collect LANL labels
  #
  dl <- read.dna(infiles.lanl,format='fa')
  dl <- tibble(TAXA=unlist(rownames(dl))) %>%
    filter(!grepl('-.-.-.',TAXA))  %>%
    mutate( SUBTYPE:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',1),
            Alpha.2.code:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',2),
            YEAR:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',3),
            GENBANK:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',5)) %>%
    distinct()
  dl$TAXA[grepl('K03455',dl$TAXA)] <- 'HXB2'
  
  outfile <- paste0(outfile.base,'sequence_labels.rda')
  save(dseq, dind, dl, dgeo, file= outfile)
  
  
 	# create labels for NL 
  #
  #   make sequence labels for NLMSM
  #
  dseqnl <- copy(dseq)
  dseqnl$SEQ_LABEL <- gsub('Amst_','NL_',dseqnl$SEQ_LABEL)
  tmp <- dseqnl %>% rename(TAXA:= SEQ_LABEL) %>%
  	mutate(GRP:= case_when(TRANSM=='MSM'~'NLMSM',
  												 TRANSM=='HSX'~'NLHSX',
  												 TRANSM=='IDU'~'NLIDU',
  												 TRUE~'NLOTH')) %>%
  	mutate(TAXA_NEW:= paste0(GRP,'___',SEQ_ID,'_',PATIENT,'_',CITY,'_',ORIGIN,'_',LOC_BIRTH,'-',GENDER,'_',TRANSM,'_',SEQ_DATE,'_',SEQ_LENGTH)) %>%
  	select(SUBTYPE,GRP,TAXA,TAXA_NEW)
  tmp <- dl %>% left_join(dgeo, by='Alpha.2.code') %>%
  	mutate(WRLD:= case_when(is.na(WRLD)~'Unknown',
  													TRUE~WRLD)) %>%
  	mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA)) %>%
  	select(SUBTYPE,WRLD,TAXA,TAXA_NEW)  %>%
  	rename(GRP=WRLD) %>%
  	rbind(tmp)  %>%
  	arrange(GRP) %>%
  	distinct()
  outfile <- paste0(outfile.base,'sequence_labels_NLMSM.csv')
  write.csv(tmp, file=outfile)
  #
  #   make sequence labels for NLHSX
  #
  tmp <- dseqnl %>% rename(TAXA:= SEQ_LABEL) %>%
  	mutate(GRP:= case_when(TRANSM=='HSX'~'NLHSX',
  												 TRANSM=='MSM'~'NLMSM',
  												 TRANSM=='IDU'~'NLIDU',
  												 TRUE~'NLOTH')) %>%
  	mutate(TAXA_NEW:= paste0(GRP,'___',SEQ_ID,'_',PATIENT,'_',CITY,'_',ORIGIN,'_',LOC_BIRTH,'-',GENDER,'_',TRANSM,'_',SEQ_DATE,'_',SEQ_LENGTH)) %>%
  	select(SUBTYPE,GRP,TAXA,TAXA_NEW)
  tmp <- dl %>% left_join(dgeo, by='Alpha.2.code') %>%
  	mutate(WRLD:= case_when(is.na(WRLD)~'Unknown',
  													TRUE~WRLD)) %>%
  	mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA)) %>%
  	select(SUBTYPE,WRLD,TAXA,TAXA_NEW)  %>%
  	rename(GRP=WRLD) %>%
  	rbind(tmp)  %>%
  	arrange(GRP) %>%
  	distinct()
  outfile <- paste0(outfile.base,'sequence_labels_NLHSX.csv')
  write.csv(tmp, file=outfile)
  
  #
  #   make sequence labels for NLIDU
  #
  tmp <- dseqnl %>% rename(TAXA:= SEQ_LABEL) %>%
  	mutate(GRP:= case_when(TRANSM=='IDU'~'NLIDU',
  												 TRANSM=='HSX'~'NLHSX',
  												 TRANSM=='MSM'~'NLMSM',
  												 TRUE~'NLOTH')) %>%
  	mutate(TAXA_NEW:= paste0(GRP,'___',SEQ_ID,'_',PATIENT,'_',CITY,'_',ORIGIN,'_',LOC_BIRTH,'-',GENDER,'_',TRANSM,'_',SEQ_DATE,'_',SEQ_LENGTH)) %>%
  	select(SUBTYPE,GRP,TAXA,TAXA_NEW)
  tmp <- dl %>% left_join(dgeo, by='Alpha.2.code') %>%
  	mutate(WRLD:= case_when(is.na(WRLD)~'Unknown',
  													TRUE~WRLD)) %>%
  	mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA)) %>%
  	select(SUBTYPE,WRLD,TAXA,TAXA_NEW)  %>%
  	rename(GRP=WRLD) %>%
  	rbind(tmp)  %>%
  	arrange(GRP) %>%
  	distinct()
  outfile <- paste0(outfile.base,'sequence_labels_NLIDU.csv')
  write.csv(tmp, file=outfile)
}


## ---- rmd.chunk.athena.220115.phyloscanner.nonB ----
athena.220115.phyloscanner.nonB <- function(analysis)
{
  require(dplyr)
  require(data.table)
  require(ape)
  require(adephylo)
  require(phytools)
  require(phangorn)
	
	analysis <- 'analysis_210504_ATHENA'
	
  
  plot.phylogenies <- 1
  max.Ntip <- 5e3
  if(0)
  {
    prog.phyloscanner_analyse_trees <- '/Users/alexb/Documents/software/phyloscanner/phyloscanner/phyloscanner_analyse_trees.R'
    home <- '/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223'
  }
  if(1)
  {
    prog.phyloscanner_analyse_trees <- '/rds/general/project/ratmann_roadmap_data_analysis/live/R/phyloscanner_analyse_trees.R'
    home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
  }
  
  
  indir.trees <- file.path(home,analysis,'fasttree')
  infile.labels <- data.table(FLABEL= c(  file.path(home,analysis,'misc','200917_sequence_labels_NLMSM.csv'),
  																				file.path(home,analysis,'misc','200917_sequence_labels_NLHSX.csv'),
  																				file.path(home,analysis,'misc','200917_sequence_labels_NLIDU.csv'))
  )
  outdir <- file.path(home,analysis,'phyloscanner')
  infiles <- data.table(FIN=list.files(indir.trees, pattern='\\.newick$',full.names=TRUE))
  infiles[, DUMMY:= 1L]
  infile.labels[, DUMMY:= 1L]
  infiles <- merge(infiles, infile.labels, by='DUMMY', allow.cartesian=TRUE)
  infiles[, DUMMY:= NULL]
  infiles[, ST:= gsub('.*_subtype_([A-Za-z0-9]+)_.*','\\1',basename(FIN))]
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
    #   update tip labels: add world region to start of label
    local.world.region <- gsub('^.*_labels_([A-Za-z]+)\\.csv$','\\1', basename(infiles[i,FLABEL]))
    dl <- as.data.table(read.csv(infiles[i,FLABEL], stringsAsFactors=FALSE))
    dl[, X:=NULL]
    stopifnot(!any(ph$tip.label==''))
    dp <- data.table(IDX=1:Ntip(ph), TAXA=ph$tip.label)
    dp <- merge(dp,dl,by='TAXA',all.x=TRUE)
    dp$TAXA_NEW[grepl('K03455',dp$TAXA)] <- 'HXB2'
    dp$GRP[grepl('K03455',dp$TAXA)] <- 'WEurope'
    dp$SUBTYPE[grepl('K03455',dp$TAXA)] <- 'B'
    dp$SUBTYPE[grepl('Outgroup',dp$TAXA)] <-  gsub('Outgroup.([0-9_A-Z]+).*','\\1',basename(dp$TAXA)[grepl('Outgroup',dp$TAXA)])
    dp$TAXA_NEW[grepl('Outgroup',dp$TAXA)] <- dp$TAXA[grepl('Outgroup',dp$TAXA)]
    stopifnot( nrow(dp[is.na(TAXA_NEW),])==0 )
    stopifnot( !any(duplicated(ph$tip.label)) )
    ph$tip.label <- dp[order(IDX),][, TAXA_NEW]
    
    # include subtypes F/A with F1/A1
    dp$SUBTYPE[dp$SUBTYPE=='F'] <- 'F1'  
    dp$SUBTYPE[dp$SUBTYPE=='A'] <- 'A1'  
		#   re-root
    tmp <- dp[, list(NST=length(TAXA)), by='SUBTYPE']
    tmp <- tmp[order(-NST),][2,SUBTYPE]
    tmp <- subset(dp, SUBTYPE==tmp,)[,TAXA_NEW]
    root <- getMRCA(ph, tmp)
    ph <- reroot(ph, root, ph$edge.length[which(ph$edge[,2]==root)]/2)
    stopifnot(is.binary(ph))
    #   write to file
    intree.phsc <- file.path(outdir,gsub('\\.newick',paste0('_rerooted_',local.world.region,'.newick'),basename(infile)))
    write.tree(ph, file=intree.phsc)
    if(plot.phylogenies)
    {
      pdf(file=gsub('newick','pdf',intree.phsc), w=20, h=10+Ntip(ph)/10)
      plot(ph, show.node.label=TRUE, cex=0.3)
      dev.off()
    }
    #   make phyloscanner UNIX command
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
    #   don t use -m multifurcation threshold to ensure that output tree is still binary
    #cmd <- paste0(cmd,' s,0 -x "',tip.regex,'" -v 1 -ow -rda -m 1e-5\n')
    cmd <- paste0(cmd,' s,0 -x "',tip.regex,'" -v 1 -ow -rda\n')
    cmd <- paste0(cmd,'mv ',basename(outputString),'* ','"',dirname(outputString),'"','\n')
    cmd <- paste(cmd, "cd $CWD\n",sep='')
    cmd <- paste(cmd, "rm ", "-r ", tmpdir,'\n',sep='')
    cmds[[i]] <- cmd
  }
  infiles[, CMD:= unlist(cmds)]
  
  #   submit jobs like this one:
  cat(infiles[1,CMD])
  
  #   run on HPC as array job
  df <- copy(infiles)
  df[, CASE_ID:= 1:nrow(df)]
  #   make PBS header
  hpc.load    <- "module load anaconda3/personal"
  r.activate  <- "source activate phylo"
  export.path <- 'export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH'
  hpc.select  <- 1                        # number of nodes
  hpc.nproc   <- 1                        # number of processors on node
  hpc.walltime<- 23                       # walltime
  hpc.q       <- NA #"pqeelab"                        # PBS queue
  hpc.mem     <- "30gb"                    # RAM
  hpc.array   <- length(unique(df$CASE_ID))   # number of runs for job array
  pbshead     <- "#!/bin/sh"
  tmp         <- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  tmp         <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  pbshead     <- paste(pbshead, "#PBS -j oe", sep = "\n")
  if(!is.na(hpc.array))
    pbshead <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
  if(!is.na(hpc.q))
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead     <- paste(pbshead, hpc.load, r.activate, export.path, sep = "\n")
  #   make array job
  cmd     <- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
  cmd     <- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]
  cmd     <- paste(pbshead,cmd,sep='\n')
  #   submit job
  outfile     <- gsub(':','',paste("phs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
  outfile     <- file.path(outdir, outfile)
  cat(cmd, file=outfile)
  cmd         <- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))
}


## ---- rmd.chunk.athena.220115.phyloscanner.B ----
athena.220115.phyloscanner.B <- function(analysis)
{
  require(dplyr)
  require(data.table)
  require(ape)
  require(adephylo)
  require(phytools)
  require(phangorn)
 
	analysis <- 'analysis_210504_ATHENA'
	
  plot.phylogenies <- 1
  max.Ntip <- 5e3
  if(0)
  {
    prog.phyloscanner_analyse_trees <- '/Users/alexb/Documents/software/phyloscanner/phyloscanner/phyloscanner_analyse_trees.R'
    home <- '/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223'
  }
  if(1)
  {
    prog.phyloscanner_analyse_trees <- '/rds/general/project/ratmann_roadmap_data_analysis/live/R/phyloscanner_analyse_trees.R'
    home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
  }
  
  
  indir.trees <- file.path(home,analysis,'fasttree')
  infile.labels <- data.table(FLABEL= c(  file.path(home,analysis,'misc','200917_sequence_labels_NLMSM.csv'),
  																				file.path(home,analysis,'misc','200917_sequence_labels_NLHSX.csv'),
  																				file.path(home,analysis,'misc','200917_sequence_labels_NLIDU.csv'))
  )
  outdir <- file.path(home,analysis,'phyloscanner')
  infiles <- data.table(FIN=list.files(indir.trees, pattern='\\.newick$',full.names=TRUE))
  infiles[, DUMMY:= 1L]
  infile.labels[, DUMMY:= 1L]
  infiles <- merge(infiles, infile.labels, by='DUMMY', allow.cartesian=TRUE)
  infiles[, DUMMY:= NULL]
  infiles[, ST:= gsub('.*_subtype_([A-Za-z0-9]+)_.*','\\1',basename(FIN))]
  infiles <- subset(infiles, ST=='B')
  set.seed(42L)
  for(i in seq_len(nrow(infiles)))
  	{
    #i<- 1
    cat('\nprocess ',infiles[i,FIN],'\nprocess ',infiles[i,FLABEL])
    infile <- infiles[i,FIN]
    ph <- read.tree(infile)
    ph$node.label <- NULL
    #   update tip labels: add world region to start of label
    local.world.region <- gsub('^.*_labels_([A-Za-z]+)\\.csv$','\\1', basename(infiles[i,FLABEL]))
    dl <- as.data.table(read.csv(infiles[i,FLABEL], stringsAsFactors=FALSE))
    dl[, X:=NULL]
    stopifnot(!any(ph$tip.label==''))
    dp <- data.table(IDX=1:Ntip(ph), TAXA=ph$tip.label)
    dp <- merge(dp,dl,by='TAXA',all.x=TRUE)
    dp$TAXA_NEW[grepl('K03455',dp$TAXA)] <- 'HXB2'
    dp$SUBTYPE[grepl('K03455',dp$TAXA)] <- 'B'
    dp$SUBTYPE[dp$SUBTYPE=='Other'] <- 'B'
    dp$GRP[grepl('K03455',dp$TAXA)] <- 'WEurope'
    dp$SUBTYPE[grepl('Outgroup',dp$TAXA)] <-  gsub('Outgroup.([0-9_A-Z]+).*','\\1',basename(dp$TAXA)[grepl('Outgroup',dp$TAXA)])
    dp$TAXA_NEW[grepl('Outgroup',dp$TAXA)] <- dp$TAXA[grepl('Outgroup',dp$TAXA)]
    stopifnot( nrow(dp[is.na(TAXA_NEW),])==0 )
    stopifnot( !any(duplicated(ph$tip.label)) )
    ph$tip.label <- dp[order(IDX),][, TAXA_NEW]
    #   re-root
    tmp <- dp[, list(NST=length(TAXA)), by='SUBTYPE']
    tmp <- tmp[order(-NST),][2,SUBTYPE]
    tmp <- subset(dp, SUBTYPE==tmp,)[,TAXA_NEW]
    root <- getMRCA(ph, tmp)
    ph <- reroot(ph, root, ph$edge.length[which(ph$edge[,2]==root)]/2)
    #   split very large tree into small trees if necessary
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
    #   write to files
    if(length(split.phs)>1)
    {
      for(k in seq_along(split.phs))
      {
        intree.phsc <- gsub('_(subtype_[A-Za-z0-9]+)_',paste0('_\\1c',k,'_'),gsub('\\.newick',paste0('_rerooted_',local.world.region,'.newick'),basename(infile)))
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
  infiles[, ST:= gsub('^.*_subtype_([A-Za-z0-9]+)_.*$','\\1',basename(FIN))]
  infiles[, DUMMY:= gsub('\\.newick$','',basename(FIN))]
  tmp <- data.table(FOUT=list.files(outdir, pattern='workspace.rda$',full.names=TRUE))
  tmp[, DUMMY:= gsub('__workspace.rda$','',basename(FOUT))]
  infiles <- merge(infiles, tmp, by='DUMMY', all.x=TRUE)
  infiles <- subset(infiles, is.na(FOUT), c(ST, FIN))
  infiles <- subset(infiles, grepl('Bc[0-9]+',ST))
  cmds <- vector('list',nrow(infiles))
  for(i in seq_len(nrow(infiles)))
  {
    #   make phyloscanner UNIX command
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
    #   don t use -m multifurcation threshold to ensure that output tree is still binary
    #cmd <- paste0(cmd,' s,0 -x "',tip.regex,'" -v 1 -ow -rda -m 1e-5\n')
    cmd <- paste0(cmd,' s,0 -x "',tip.regex,'" -v 1 -ow -rda\n')
    cmd <- paste0(cmd,'mv ',basename(outputString),'* ','"',dirname(outputString),'"','\n')
    cmd <- paste(cmd, "cd $CWD\n",sep='')
    cmd <- paste(cmd, "rm ", "-r ", tmpdir,'\n',sep='')
    cmds[[i]] <- cmd
  }
  infiles[, CMD:= unlist(cmds)]
  
  #   submit jobs like this one:
  cat(infiles[1,CMD])
  
  #   run on HPC as array job
  df <- copy(infiles)
  df[, CASE_ID:= 1:nrow(df)]
  #   make PBS header
  hpc.load    <- "module load anaconda3/personal"
  r.activate  <- "source activate phylo"
  hpc.select  <- 1                        # number of nodes
  hpc.nproc   <- 1                        # number of processors on node
  hpc.walltime<- 23                       # walltime
  hpc.q       <- NA #"pqeelab"                        # PBS queue
  hpc.mem     <- "30gb"                    # RAM
  hpc.array   <- length(unique(df$CASE_ID))   # number of runs for job array
  pbshead     <- "#!/bin/sh"
  tmp         <- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  tmp         <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  pbshead     <- paste(pbshead, "#PBS -j oe", sep = "\n")
  if(!is.na(hpc.array))
    pbshead <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
  if(!is.na(hpc.q))
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead     <- paste(pbshead, hpc.load, r.activate, sep = "\n")
  #   make array job
  cmd     <- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
  cmd     <- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]
  cmd     <- paste(pbshead,cmd,sep='\n')
  #   submit job
  outfile     <- gsub(':','',paste("phs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
  outfile     <- file.path(outdir, outfile)
  cat(cmd, file=outfile)
  cmd         <- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))
}

## ---- rmd.chunk.athena.220115.fastree.B ----
athena.220115.fastree.B<- function(analysis)
{
	require(big.phylo)
	require(data.table)
	
	analysis <- 'analysis_210504_ATHENA'
	
	home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	
	indir.trees <- file.path(home,analysis,'phyloscanner_B')
	indir.aln <- file.path(home,'analysis_210226_ATHENA','alignments_bs')
	indir.trees <- '/rds/general/project/ratmann_roadmap_data_analysis/live/analysis_210226_ATHENA/phyloscanner_B'
	
	infile.labels <- data.table(FLABEL= c(  file.path(home,analysis,'misc','200917_sequence_labels_NLMSM.csv'),
																					file.path(home,analysis,'misc','200917_sequence_labels_NLHSX.csv'))
	)
	outdir <- file.path(home,analysis,'phyloscanner')
	infiles <- data.table(FIN=list.files(indir.trees, pattern='\\.newick$',full.names=TRUE))
	infiles[, DUMMY:= 1L]
	infile.labels[, DUMMY:= 1L]
	infiles <- merge(infiles, infile.labels, by='DUMMY', allow.cartesian=TRUE)
	infiles[, DUMMY:= NULL]
	infiles[, ST:= gsub('.*_subtype_([A-Za-z0-9]+)_.*','\\1',basename(FIN))]
	infiles[, REP:= gsub('.*_wOutgroup_([A-Za-z0-9]+)_.*','\\1',basename(FIN))]
	#infiles <- subset(infiles, ST=='B')
	
	indir <- file.path(home,analysis,'alignments_bs')
	outdir <- file.path(home,analysis,'fasttree_B')
	
	#   get UNIX commands for all trees to be processed
	df <- data.table(FIN=list.files(indir.aln, pattern='fasta$',full.names=TRUE))
	df <- subset(df, grepl('noDRM',FIN))
	df[, ST:= gsub('.*_subtype_([A-Z0-9a-z]+)_.*','\\1',basename(FIN))]
	df[, FOUT:= file.path(outdir,gsub('\\.fasta','_ft.newick',basename(FIN)))]
	df[, REP:= gsub('.*_wOutgroup_([A-Za-z0-9]+).*','\\1',basename(FIN))]
	df <- subset(df, ST=='B')
	
	set.seed(42L)
	for(i in seq_len(nrow(infiles)))
	{
		#i<- 1
		cat('\nprocess ',infiles[i,FIN],'\nprocess ',infiles[i,FLABEL])
		infile <- infiles[i,FIN]
		rep <- infiles[i,REP]
		ph <- read.tree(infile)
		ph$node.label <- NULL
		#   update tip labels: add world region to start of label
		local.world.region <- gsub('^.*_labels_([A-Za-z]+)\\.csv$','\\1', basename(infiles[i,FLABEL]))
		dl <- as.data.table(read.csv(infiles[i,FLABEL], stringsAsFactors=FALSE))
		dl[, X:=NULL]
		dl <- dl[dl$TAXA_NEW %in% ph$tip.label,]
		infile.fasta <- df[REP==rep,FIN]
		sq <- read.dna(infile.fasta, format='fa')
		sq <- sq[labels(sq) %in% c(dl[, TAXA],grep('Outgroup',labels(sq),fixed=1)),]
		outfile <- gsub('_B_',paste0('_',infiles[i,ST],'_'),infile.fasta)
		outfile <- file.path(outdir,basename(outfile))
		write.dna(sq, file=outfile, format='fa', colsep='', nbcol=-1)
	}
	
	#df <- subset(df, grepl('subtype_B',FIN))
	cmds <- vector('list', nrow(df))
	for(i in seq_len(nrow(df)))
	{
		infile.fasta <- df[REP==rep,FIN]
		sq <- read.dna(infile.fasta, format='fa')
		sq <- sq[labels(sq) %in% c(dl[, TAXA],grep('Outgroup',labels(sq),fixed=1)),]

		outfile <- gsub('\\.fasta',paste0('_',sprintf("%03d",b),'.fasta'),infile.fasta)
			outfile <- file.path(outdir,basename(outfile))
			write.dna(sq[,bs], file=outfile, format='fa', colsep='', nbcol=-1)

		outfile <- df[i,FOUT]
		tmp <- cmd.fasttree(infile.fasta, outfile=outfile, pr.args='-nt -gtr -gamma', check.binary=TRUE)
		cmds[[i]] <- tmp
	}
	df[, CMD:= unlist(cmds)]
	
	#   submit jobs like this one:
	cat(df[1,CMD])
	
	#   run on HPC as array job
	df[, CASE_ID:= 1:nrow(df)]
	
	#   make PBS header
	#hpc.load    <- "module load R/3.4.0"
	hpc.load    <- "module load anaconda3/personal"
	r.activate  <- "source activate Renv3"
	hpc.select  <- 1                        # number of nodes
	hpc.nproc   <- 1                        # number of processors on node
	hpc.walltime<- 72                      # walltime
	hpc.q       <- NA                       # PBS queue
	hpc.mem     <- "6gb"                    # RAM
	hpc.array   <- length(unique(df$CASE_ID))   # number of runs for job array
	pbshead     <- "#!/bin/sh"
	tmp         <- paste("#PBS -l walltime=", hpc.walltime, ":00:00,pcput=", hpc.walltime, ":00:00", sep = "")
	pbshead     <- paste(pbshead, tmp, sep = "\n")
	tmp         <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
	pbshead     <- paste(pbshead, tmp, sep = "\n")
	pbshead     <- paste(pbshead, "#PBS -j oe", sep = "\n")
	if(!is.na(hpc.array))
		pbshead <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
	if(!is.na(hpc.q))
		pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
	pbshead     <- paste(pbshead, hpc.load, r.activate, sep = "\n")
	#   make array job
	cmd     <- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
	cmd     <- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]
	cmd     <- paste(pbshead,cmd,sep='\n')
	#   submit job
	outfile     <- gsub(':','',paste("trs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
	outfile     <- file.path(outdir, outfile)
	cat(cmd, file=outfile)
	cmd         <- paste("qsub", outfile)
	cat(cmd)
	cat(system(cmd, intern= TRUE))
}


## ---- rmd.chunk.athena.220115.get.subgraphs ----
athena.220115.phyloscanner.get.subgraphs<- function(analysis)
{
  require(data.table)
  require(phangorn)
  require(ggplot2)
  require(reshape)
  require(phyloscannerR)
	require(big.phylo)
  
  #    working directory with phyloscanner output
	analysis <- 'analysis_210226_ATHENA'
	
  home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
  indir.phsc    <- file.path(home,analysis,'phyloscanner')
  outdir <- file.path(home,analysis,'subgraphs')
  infiles        <- data.table(F=list.files(indir.phsc, pattern='_workspace.rda$', full.names=TRUE, recursive=TRUE))
  infiles[, SELECT:= gsub('^.*_rerooted_([A-Za-z0-9]+)_.*$','\\1',basename(F))]

  #    extract subgraphs of dated trees
  for(i in seq_len(nrow(infiles)))
    {
    cat('process', i,'\n')
    infile <- infiles[i, F]
    host <- infiles[i,SELECT]
    load(infile)
    ph <- phyloscanner.trees[[1]][['tree']]
    
    mrcas <- which( attr(ph, 'SUBGRAPH_MRCA') )
    #    some of the tips states are "unknown" and this clashes internally with the NA state, so we need to take some extra care
    #    this is not a problem because the "unknown" and NA state mean the same thing
    attr(ph, 'INDIVIDUAL') <- as.character(attr(ph, 'INDIVIDUAL'))
    attr(ph, 'INDIVIDUAL')[is.na(attr(ph, 'INDIVIDUAL'))] <- 'Unknown'
    mrcas <- mrcas[ attr(ph, 'INDIVIDUAL')[mrcas]==host ]
    stopifnot( !any(is.na(mrcas)) )
    ph <- phyloscanner.to.simmap(ph)
    # check that tree is of class simmap
    stopifnot( any(attr(ph,'class')=='simmap') )
    # extract subgraphs
    subgraphs <- lapply(mrcas, function(mrca) extract.subgraph(ph, mrca))
    # save
    outfile <- gsub('_workspace',paste0('_subgraphs_',host),basename(infile))
    outfile <- file.path(outdir,outfile)
    save(subgraphs, file=outfile)
  }
}

### Read phylogenetic subgraphs
require(data.table)
require(phangorn)
require(ggplot2)
require(reshape)

#    working directory with phyloscanner output
home <- '~/Box Sync/Roadmap'
analysis <- 'analysis_210226_ATHENA'
home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
indir.phsc    <- file.path(home,analysis,'subgraphs_dated')
infile.meta <- file.path(home,analysis,'misc','200917_sequence_labels.rda')
outdir <- file.path(home,analysis,'subgraph_analysis')

#    extract subgraph taxa
infiles <- data.table(F=list.files(indir.phsc, pattern='rda$', full.names=TRUE, recursive=TRUE))
infiles <- subset(infiles, grepl('subgraphs',F))
#    SELECT defines if Ams, AmsHSX, AmsMSM
infiles[, SELECT:= gsub('^(.*)subgraphs_([A-Za-z]+).rda','\\2',basename(F))]
#    for subtype B, I ran separate analyses based on very large subtrees for comp efficiency
#    ST stores the subtype, and ST_CLADE the large subtree
infiles[, ST:= gsub('^.*subtype_([^_]+)_.*\\.rda','\\1',basename(F))]
# For B only
infiles[, ST_CLADE:= as.integer(gsub('[^c]*c([0-9]+)','\\1',ST))]
infiles[, ST:= gsub('([^c]*)c([0-9]+)','\\1',ST)]
#    ID of bootstrap replicate, 000 for analysis on real alignment
infiles[, REP:= gsub('^.*wOutgroup_([0-9]+)_.*\\.rda','\\1',basename(F))]
infiles[,REPSEL:=paste0(ST,"_",SELECT,"_",REP)]
dsubgraphtaxa <- infiles[, {
  infile <- F
  cat('Process',infile,'\n')
  load(infile)
    if(length(subgraphs)==1){
        subgraph.names <- rep(subgraphs[[1]]$subgraph.name, length(subgraphs[[1]]$tip.label))
        subgraph.taxa <- subgraphs[[1]]$tip.label
        subgraph.parent.state <- subgraphs[[1]]$subgraph.parent.state
  }  else{
    subgraph.names <- unlist(sapply( subgraphs, function(subgraph)  rep(subgraph$subgraph.name, length(subgraph$tip.label)) ))
    subgraph.taxa <- unlist(sapply( subgraphs, function(subgraph)  subgraph$tip.label))
    subgraph.parent.state <- unlist(sapply( subgraphs, function(subgraph)  rep(subgraph$subgraph.parent.state, length(subgraph$tip.label))))
  }
  list(NAME=subgraph.names,
      TAXA= subgraph.taxa,
      ORIGINHOST= subgraph.parent.state
  )
}, by=c('ST','ST_CLADE','REP','SELECT')]
#    add meta data from taxa names
regex.tip.label <- '^([A-Za-z]+)_+(T[0-9]+)_([0-9]+)_([a-zA-Z-]+)_([a-zA-Z]+)_([a-zA-Z]+)-([a-zA-Z]+)_([a-zA-Z_]+)_([0-9]+)-([0-9-]+)-([0-9-]+)_([0-9]+)$'
dsubgraphtaxa[, ID:= as.numeric(gsub(regex.tip.label,'\\3',TAXA))]
dsubgraphtaxa[, SEQ_YEAR:= as.numeric(gsub(regex.tip.label,'\\9',TAXA))]
dsubgraphtaxa[, FULL_NAME:= paste0(ST,'_',REP,'_',SELECT,'_',NAME)]

# Number of distinct chains
nsubgraph <- dsubgraphtaxa[, list(N_SUBGRAPHS=length(NAME)), by=c('REP','SELECT')]

#   add meta data from pre-processed persons file
load(infile.meta)
dind <- as.data.table(dind)
setnames(dind, c('PATIENT','BIRTH_Y','BIRTH_CNTRY'), c('ID','BIRTH_YEAR','BIRTH_COUNTRY'))
tmp <- subset(dind, select=c(ID,BIRTH_YEAR,BIRTH_COUNTRY,LOC_BIRTH,CITY,TRANSM,GENDER))
set(tmp, NULL, 'BIRTH_YEAR', tmp[, as.numeric(BIRTH_YEAR)])
dsubgraphtaxa <- merge(dsubgraphtaxa,tmp,by='ID')

dsubgraphtaxa <- unique(dsubgraphtaxa)
outfile <- file.path(outdir,'subgraphs_withmetadata.rda')
save(dsubgraphtaxa, file=outfile)

# summarise subgraph data and individuals in subgraphs
regex.tip.label <- '^([A-Za-z]+)_+(T[0-9]+)_([0-9]+)_([a-zA-Z-]+)_([a-zA-Z]+)_([a-zA-Z]+)-([a-zA-Z]+)_([a-zA-Z_]+)_([0-9-]+)_([0-9]+)$' # seq update has T prefix to seqID

dsubgraphtaxa[, SEQ_DATE:= as.Date(gsub(regex.tip.label,'\\9',TAXA))]
sgs <- dsubgraphtaxa[,list(MIN_SAMPLING_DATE=min(SEQ_DATE),MAX_SAMPLING_DATE=max(SEQ_DATE),
													 SIZE=length(TAXA)),by=c('ST','ST_CLADE','REP','SELECT','NAME','ORIGINHOST')]

setnames(sgs,c('SELECT','ORIGINHOST'),c('TRM_GROUP','PARENT_STATE'))
setnames(dsubgraphtaxa,c('SELECT','ORIGINHOST'),c('TRM_GROUP','PARENT_STATE'))

outfile <- file.path(outdir,'subgraph_data_210511.csv')
write.csv(sgs,file=outfile)

dsubgraphtaxa[, SAMPLING_COUNTRY:='Netherlands']
dsubgraphtaxa[is.na(PARENT_STATE),PARENT_STATE:='Unknown']

outfile <- file.path(outdir,'individuals_in_subgraph_data_210511.csv')
write.csv(dsubgraphtaxa,file=outfile)