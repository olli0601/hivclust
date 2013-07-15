project.gccontent.collect.data<- function(dir.name, select.MSM=0, select.max.months.after.PosEarliest= 30, select.min.months.after.PosEarliest=6)
{
	#read ATHENA sequences 		
	indir				<- paste(dir.name,"tmp",sep='/')
	infile.seq			<- "ATHENA_2013_03_CurAll+LANL_Sequences"	
	insignat			<- "Sat_Jun_16_17/23/46_2013"
	file				<- paste(indir,'/',paste(infile.seq,'_',gsub('/',':',insignat),".R", sep=''), sep='')
	if(verbose) cat(paste("\nloading file",file))
	load(file)	
	seq.PROT.RT			<- seq.PROT.RT[substr(rownames(seq.PROT.RT),1,2)!="TN", ]
	seq.PROT.RT			<- seq.PROT.RT[substr(rownames(seq.PROT.RT),1,8)!="PROT+P51", ]
	if(verbose) cat(paste("\nfound ATHENA seqs, n=",nrow(seq.PROT.RT)))
	
	#read patient data
	indir				<- paste(dir.name,"derived",sep='/')
	infile.patient		<- "ATHENA_2013_03_All1stPatientCovariates"
	file				<- paste(indir,'/',paste(infile.patient,".R", sep=''), sep='')
	if(verbose) cat(paste("\nloading file",file))
	load(file)
	df.patient			<- df.all	
	df.all				<- NULL
	#define first pos event
	tmp					<- df.patient[, AnyTherT1]
	tmp2				<- df.patient[, is.na(AnyTherT1) | (!is.na(PosRNA1) &  PosRNA1<AnyTherT1) ]
	tmp[tmp2]			<- df.patient[tmp2, PosRNA1]
	tmp2				<- df.patient[, !is.na(PosCD41) & !is.na(tmp) & PosCD41<tmp ]
	tmp[tmp2]			<- df.patient[tmp2, PosCD41]
	tmp2				<- df.patient[, !is.na(PosSeq) & !is.na(tmp) & PosSeq<tmp ]
	tmp[tmp2]			<- df.patient[tmp2, PosSeq]
	tmp2				<- df.patient[, !is.na(PosT) & !is.na(tmp) & PosT<tmp ]
	if(verbose) cat(paste("\nfound PosT is earliest pos event, n=",length(which(tmp2))," out of n=",nrow(df.patient)))
	tmp[tmp2]			<- df.patient[tmp2, PosT]
	df.patient[, PosEarliest:=tmp]
	if(	length(which(df.patient[, PosEarliest>PosT])) 		||
		length(which(df.patient[, PosEarliest>PosSeq])) 	|| 
		length(which(df.patient[, PosEarliest>PosRNA1]))	||
		length(which(df.patient[, PosEarliest>PosCD41]))	||
		length(which(df.patient[, PosEarliest>AnyTherT1]))		) 	stop("something is wrong with PosEarliest")			
	
	#read RNA data and merge with first PosEarliest, NegT, NegT_Acc, isAcute 
	indir				<- paste(dir.name,"derived",sep='/')
	infile.viro			<- "ATHENA_2013_03_Viro"
	file				<- paste(indir,'/',paste(infile.viro,".R", sep=''), sep='')
	if(verbose) cat(paste("\nloading file",file))
	load(file)
	df.viro				<- as.data.table(df)
	setkey(df.viro, "Patient")
	df.viro				<- subset(df.viro, !is.na(DateRNA) & Undetectable=="No", c(Patient, DateRNA, RNA))
	df					<- NULL
	df.viro				<- merge(df.viro, subset(df.patient, select=c(Patient,PosEarliest, NegT, NegT_Acc, isAcute)))
	
	#read sampling date of sequences
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
	#sanity checks
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
	gc.seqinfo			<- subset(df, !is.na(PosSeqT), select=c(SampleCode, Patient, PosSeqT))
	if(verbose) cat(paste("\nfound ATHENA seqs with known sampling date, n=",nrow(gc.seqinfo)))
	if(verbose) cat(paste("\nfound ATHENA patients whose sequences have known sampling date, n=",length(unique(gc.seqinfo[,Patient]))))
		
	#keep only ATHENA sequences with known sampling date	and 	sampling date < treatment start	
	df					<- subset(df.patient,!is.na(AnyTherT1),select=c(Patient,AnyTherT1))
	if(verbose) cat(paste("\nfound ATHENA patients with known therapy date, n=",nrow(df)))
	setkey(gc.seqinfo,"Patient")	
	gc.seqinfo			<- merge(gc.seqinfo,df,all=TRUE)
	gc.seqinfo			<- subset(gc.seqinfo, !is.na(PosSeqT) & !is.na(AnyTherT1))
	if(verbose) cat(paste("\nfound ATHENA seqs with known sampling date and known therapy date, n=",nrow(gc.seqinfo)))
	gc.seqinfo			<- subset(gc.seqinfo, PosSeqT<=AnyTherT1)
	if(verbose) cat(paste("\nfound ATHENA seqs with sampling date < therapy date, n=",nrow(gc.seqinfo)))	
	tmp					<- gsub(' ','',gc.seqinfo[,SampleCode])
	gc.seqinfo[,SampleCodeFASTA:=tmp]	
	
	#keep only ATHENA sequences with known MSM transmission route
	if(select.MSM)
	{
		df					<- subset(df.patient,!is.na(Trm) & Trm=="MSM",Patient)
		if(verbose) cat(paste("\nfound ATHENA patients with known MSM, n=",nrow(df)))
		gc.seqinfo			<- merge(gc.seqinfo, df)
		if(verbose) cat(paste("\nfound ATHENA patients with sampling date < therapy date and known MSM, n=",nrow(df)))
	}
	
	#keep those without AIDS defining event before sampling date -- would like to but don t have it
	
	#consider RNA for these patients and keep those before therapy 			 
	df.viro				<- merge(gc.seqinfo, df.viro)
	if(verbose) cat(paste("\nfound viral loads for candidate sequences, n=",nrow(df.viro)))
	df.viro				<- subset(df.viro, DateRNA<=AnyTherT1)
	if(verbose) cat(paste("\nfound viral loads for candidate sequences with RNA<=AnyTherapyT1, n=",nrow(df.viro)))
	#keep RNA for which there is at least one RNA before 'select.max.months.after.PosEarliest' mo after PosEarliest		
	df.viro				<- subset(df.viro, DateRNA<=(PosEarliest+30*select.max.months.after.PosEarliest))
	if(verbose) cat(paste("\nfound viral loads for candidate sequences with RNA<=max months after PosEarliest, n=",nrow(df.viro)))
	#for those for which there is no reliable seroNeg date, 
	#keep RNA for which there is at least one RNA after 'select.min.months.after.PosEarliest' mo after PosEarliest	
	if(0)
	{
		#TODO treat seroneg and seropos separately
		#TODO need not lump NegT_Acc="No" separately into seropos
		df.viro.seropos		<- subset(df.viro, is.na(NegT) | is.na(NegT_Acc) | NegT_Acc=="No" )							
		if(verbose) cat(paste("\nfound viral loads for candidate sequences with RNA<=max months after PosEarliest for seropos, n=",nrow(df.viro.seropos)))
		df.viro.seropos		<- subset(df.viro.seropos, DateRNA>=(PosEarliest+30*select.min.months.after.PosEarliest))
		if(verbose) cat(paste("\nfound viral loads for candidate sequences with min months after PosEarliest<=RNA<=max months after PosEarliest for seropos, n=",nrow(df.viro.seropos)))
	}
	else
	{
		df.viro				<- subset(df.viro, DateRNA>=(PosEarliest+30*select.min.months.after.PosEarliest))
		if(verbose) cat(paste("\nfound viral loads for candidate sequences with min months after PosEarliest<=RNA<=max months after PosEarliest, n=",nrow(df.viro)))
	}
	if(1)	#check for DateRNA duplicates
	{
		DateRNA.duplicate	<- df.viro[, length(unique(Patient))<length(Patient), by=DateRNA]
		DateRNA.duplicate	<- subset(DateRNA.duplicate, V1, DateRNA)
		DateRNA.duplicate	<- merge(DateRNA.duplicate, df.viro, all.x=1)
		tmp					<- DateRNA.duplicate[, length(unique(DateRNA))<length(DateRNA), by=Patient]
		tmp					<- subset(tmp, V1, Patient )
		DateRNA.duplicate	<- merge(tmp, DateRNA.duplicate, all.x=1)		
		file				<- paste(indir,'/',paste("DateRNAduplicates.csv", sep=''), sep='')
		write.table(subset(DateRNA.duplicate,select=c(Patient, DateRNA, RNA, SampleCode)), file=file)
	}
	gc.seq				<- seq.PROT.RT[ rownames(seq.PROT.RT)%in%gc.seqinfo[,SampleCodeFASTA], ]
	
	file				<- paste(outdir,'/',paste(outfile,'_',gsub('/',':',outsignat),".R", sep=''), sep='')
	if(verbose) cat(paste("\nwrite file",file))
	save(gc.seq, gc.seqinfo, file=file)
	list(gc.seq=gc.seq, gc.seqinfo=gc.seqinfo)
}

project.gccontent<- function(dir.name=DATA)
{	
	require(data.table)
	
	verbose				<- 1
	resume				<- 0
	outdir				<- paste(dir.name,"gcprjct",sep='/')
	outfile				<- "ATHENA_2013_03_GC_Master"	
	outsignat			<- "Thu_Jun_19_12/02/34_2013"
	
	#read master data 
	if(resume)		
	{
		options(show.error.messages = FALSE)
		file				<- paste(outdir,'/',paste(outfile,'_',gsub('/',':',outsignat),".R", sep=''), sep='')
		readAttempt<-try(suppressWarnings(load(file)))
		if(verbose && !inherits(readAttempt, "try-error"))			cat(paste("\nresumed file",file))
		options(show.error.messages = TRUE)
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{
		if(verbose) cat(paste("\ncall project.gccontent.collect.data",file))
		tmp			<- project.gccontent.collect.data(dir.name)
		gc.seq		<- tmp[["gc.seq"]]
		gc.seqinfo	<- tmp[["gc.seqinfo"]]
	}
	
	#start from PROT15 ATAGGGGGGCAA after the slip site; this is pos 43 in gc.seq; allow extension past p51
	gc.seq			<- gc.seq[,43:ncol(gc.seq)]
	
	gc.seq.codon1	<- gc.seq[,seq.int(1,ncol(gc.seq),3)]
	gc.seq.codon2	<- gc.seq[,seq.int(2,ncol(gc.seq),3)]
	gc.seq.codon3	<- gc.seq[,seq.int(3,ncol(gc.seq),3)]	
	gc.codon1		<- hiv.seq.gc.content(gc.seq.codon1)
	gc.codon2		<- hiv.seq.gc.content(gc.seq.codon2)
	gc.codon3		<- hiv.seq.gc.content(gc.seq.codon3)
	#TODO https://en.wikipedia.org/wiki/Genetic_code
	#depending on the amino acids, the nucleotides may be differentially variable ie for F, only the 3rd codon is variable, but for L also the first one is variable etc
	#different sites are also differentially conserved

	#TODO 
	#compute viral load measure - set point viral load ?
}
