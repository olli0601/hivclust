######################################################################################
Sequences.200205<- function(dir.name= DATA)
{
	require(data.table)
	require(gdata)
	require(ape)
	
	#	read SEQUENCE csv data file and preprocess
	dir.name		<- '~/Box/OR_Work/AIDSFonds/data_191223'
	file			<- file.path(dir.name,"SHM_1902_ROADMAP_191223_tblLAB_RES.csv")
	outfile			<- file.path(dir.name,"SHM_1902_ROADMAP_191223_tblLAB_seq.rda")
	#					
	verbose			<- 1	
	NA.DateRes		<- c("1911-11-11","")
	NA.Seq			<- c('')
	min.seq.len		<- 21
	proc.GeneCode	<- c(1,2)
	
	
	
	
	df <- read.csv(file, stringsAsFactors=FALSE)
	df <- as.data.table(df)
	setnames(df, colnames(df), toupper(colnames(df)))
	setnames(df, 'TEST_ID', 'SEQ_ID')
	#set(df, NULL, 'GENECODE', df[, as.character(factor(GENECODE, levels=c('1','2'), labels=c('PROT','RT')))])
	#	process sequencing dates
	tmp <- df[, which(SAMPLE_D%in%NA.DateRes)]
	if(verbose) cat(paste("\nentries with missing sequence dates, n=", length(tmp)))
	#	exclude seqs with missing dates
	df <- subset(df, !SAMPLE_D%in%NA.DateRes)
	#	exclude integrase
	df <- subset(df, SEQTYPE!='IN')
			
	set(df, NULL, 'SAMPLE_D', df[,as.Date(SAMPLE_D, format="%Y-%m-%d")])
	if(verbose) cat(paste("\nrange of sequence dates is",paste(range(df[,SAMPLE_D], na.rm=1),collapse=', ')))	
	
	#	define partial pol
	df				<- df[, {
				if(length(SEQ_NUC)>2) stop('more than two sequence fragments for same ID?', SEQ_ID[1])
				z	<- ''
				if(any(SEQTYPE=='PRO'))
					z<- paste0(z,SEQ_NUC[SEQTYPE=='PRO'])
				if(any(SEQTYPE=='RT'))
					z<- paste0(z,SEQ_NUC[SEQTYPE=='RT'])				
				list(	SEQ_D=SAMPLE_D[1], 						
						SEQUENCE=z)
			}, by=c('PATIENT','SEQ_ID')]
	
	#	define gene length
	tmp			<- df[, list(SEQ_L=nchar(SEQUENCE)), by=c('PATIENT','SEQ_ID')]	
	df			<- merge(df, tmp, by=c('PATIENT','SEQ_ID'))
	#	warning if no sequence
	if( nrow(subset(df, SEQ_L==0)) ) 
		cat('\nwarning: no sequence found for', subset(df, SEQ_L==0)[, paste(SEQ_ID, collapse=', ')])
	df			<- subset(df, SEQ_L>0)
	#	fixup seq IDs
	set(df, NULL, 'SEQ_ID', df[, gsub(' ','-',SEQ_ID)])
	#	add AMST to sequence labels
	set(df, NULL, 'SEQ_LABEL', df[, paste0('Amst_',PATIENT,'_',SEQ_D,'_',SEQ_ID)])
	
		
	#
	#	done reading. transform data to two objects: sequence info (data.table) and sequences (DNAbin)
	#
	ds			<- copy(df)
	set(ds, NULL, 'SEQUENCE', ds[, gsub('?','n',tolower(SEQUENCE),fixed=1)])
	
	seq			<- as.DNAbin(strsplit(ds[, SEQUENCE],''))
	names(seq)	<- ds[,SEQ_LABEL]
	write.dna(seq, file=gsub('rda','fasta',outfile), format='fasta', colsep='', nbcol=-1)
	
	ds[, SEQUENCE:=NULL]		
	#
	#	save Dutch sequences
	#	
	save(ds, file=outfile)		
	
	quit("no")	
}