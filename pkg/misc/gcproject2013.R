project.gccontent.collect.data<- function(dir.name, select.MSM=0, select.max.months.after.AnyPos_T1= 30, select.min.months.after.AnyPos_T1=6)
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
	infile.patient		<- "ATHENA_2013_03_AllSeqPatientCovariates"
	file				<- paste(indir,'/',paste(infile.patient,".R", sep=''), sep='')
	if(verbose) cat(paste("\nloading file",file))
	load(file)
	df.patient			<- df.all
	df.all				<- NULL
	setnames(df.patient, c("PosRNA","lRNA"), c("PosRNA.closestToPosSeqT","lRNA.closestToPosSeqT"))
	if(verbose) cat(paste("\nnumber of ATHENA seqs info, n=",nrow(df.patient)))
	df.patient			<- subset(df.patient, !is.na(PosSeqT))
	if(verbose) cat(paste("\nnumber of ATHENA seqs with known sampling date, n=",nrow(df.patient)))
	if(verbose) cat(paste("\nnumber of ATHENA patients whose sequences have known sampling date, n=",length(unique(df.patient[,Patient]))))
	#
	# keep only ATHENA sequences with known sampling date	and 	sampling date < treatment start
	#
	tmp					<- subset(df.patient, is.na(AnyT_T1))
	if(verbose) cat(paste("\nnumber of ATHENA seqs with is.na(AnyT_T1), n=",nrow(tmp),"IGNORE FOR NOW"))	
	df.patient			<- subset(df.patient, !is.na(AnyT_T1) & PosSeqT<=AnyT_T1)
	if(verbose) cat(paste("\nnumber of ATHENA seqs with PosSeqT<=AnyT_T1, n=",nrow(df.patient)))
	#
	# keep only ATHENA sequences with known MSM transmission route
	#
	if(select.MSM)
	{
		df.patient		<- subset(df.patient,Trm=="MSM")
		if(verbose) cat(paste("\nfound ATHENA seqs with known MSM, n=",nrow(df.patient)))		
	}
	#
	# keep only ATHENA sequences with sampling date < CDC-C event
	#
	tmp					<- subset(df.patient, !is.na(DateFirstEverCDCC) & PosSeqT>DateFirstEverCDCC)
	if(verbose) cat(paste("\nnumber of ATHENA seqs with CDC-C event and PosSeqT>DateFirstEverCDCC, n=",nrow(tmp),"-- ignoring these"))
	df.patient			<- subset(df.patient, is.na(DateFirstEverCDCC) | (!is.na(DateFirstEverCDCC) & PosSeqT<=DateFirstEverCDCC))
	if(verbose) cat(paste("\nnumber of ATHENA seqs with PosSeqT<=DateFirstEverCDCC, n=",nrow(df.patient)))
	#
	# reset seroconverters with inaccurate info, so they are definitely Neg at NegT
	#
	df.patient			<- hivc.db.reset.inaccurateNegT(df.patient, nacc.dy.dy= 1, nacc.mody.mo= 0, nacc.mody.dy= 1, verbosee=1)	
	#
	# read RNA data and merge with first AnyPos_T1, NegT, NegT_Acc, isAcute
	#
	indir				<- paste(dir.name,"derived",sep='/')
	infile.viro			<- "ATHENA_2013_03_Viro"
	file				<- paste(indir,'/',paste(infile.viro,".R", sep=''), sep='')
	if(verbose) cat(paste("\nloading file",file))
	load(file)
	df.viro				<- df
	df					<- NULL
	df.viro				<- subset(df.viro, !is.na(PosRNA) & Undetectable=="No", c(Patient, PosRNA, lRNA))
	if(verbose) cat(paste("\nnumber of viral detectable and non-missing VL measurements, n=",nrow(df.viro)))	
	setkey(df.viro, Patient)
	if(verbose) cat(paste("\nnumber of patients with viral detectable and non-missing VL measurments, n=",nrow(unique(df.viro))))		
	setkey(df.patient, Patient)
	df.viro				<- merge(df.viro, df.patient, by="Patient", allow.cartesian=T)
	if(verbose) cat(paste("\nnumber of VL measurements crossed with ATHENA seq, n=",nrow(df.viro)))
	if(verbose) cat(paste("\nnumber of selected patients with VL measurements and ATHENA seq, n=",nrow(unique(df.viro))))
	#
	# keep measurements before therapy
	#
	df.viro				<- subset(df.viro, PosRNA<AnyT_T1)
	if(verbose) cat(paste("\nnumber of VL measurements with PosRNA<AnyT_T1, n=",nrow(df.viro)))
	#
	df.viro[, NegT2PosRNA:=df.viro[,PosRNA-NegT]]
	#
	#  keep measurements after 6mo of first pos event
	#
	select.max.months.AnyPos_T1.after.NegT<- 12

	df.viro.before6mo		<- subset(df.viro, PosRNA<(AnyPos_T1+30*select.min.months.after.AnyPos_T1))
	table( df.viro.before6mo[,list(nRNA= length(lRNA)),by="Patient"][,nRNA] )
	
	
	
	df.viro.before6mo.sn	<- subset(df.viro.before6mo, !is.na(NegT) & AnyPos_T1<(NegT+30*select.max.months.AnyPos_T1.after.NegT) )
	table( df.viro.before6mo.sn[,list(nRNA= length(lRNA)),by="Patient"][,nRNA] )
	
	subset(df.viro.before6mo.sn, lRNA<3, c(Patient, PosRNA, lRNA, isAcute, NegT, NegT_Acc, Trm))
	
	df.viro.spvl			<- subset(df.viro, PosRNA>=(AnyPos_T1+30*select.min.months.after.AnyPos_T1) & PosRNA<=(AnyPos_T1+30*select.max.months.after.AnyPos_T1))
	table( df.viro.spvl[,list(nRNA= length(lRNA)),by="Patient"][,nRNA] )
	
	df.viro.spvl.sn			<- subset(df.viro.spvl, !is.na(NegT) & AnyPos_T1<(NegT+30*select.max.months.AnyPos_T1.after.NegT) )
	table( df.viro.spvl.sn[,list(nRNA= length(lRNA)),by="Patient"][,nRNA] )
	
	
	plot( df.viro.before6mo.sn[,NegT2PosRNA], df.viro.before6mo.sn[,lRNA], pch=19 )
	
	
	
	hist(df.viro.before6mo[,lRNA])
	hist(df.viro.before6mo.sn[,lRNA])	
	hist(df.viro.spvl[,lRNA])
	hist(df.viro.spvl.sn[,lRNA])
	
stop()



	
	#consider RNA for these patients and keep those before therapy 			 
	df.viro				<- merge(gc.seqinfo, df.viro)
	if(verbose) cat(paste("\nfound viral loads for candidate sequences, n=",nrow(df.viro)))
	df.viro				<- subset(df.viro, DateRNA<=AnyTherT1)
	if(verbose) cat(paste("\nfound viral loads for candidate sequences with RNA<=AnyTherapyT1, n=",nrow(df.viro)))
	#keep RNA for which there is at least one RNA before 'select.max.months.after.AnyPos_T1' mo after AnyPos_T1		
	df.viro				<- subset(df.viro, DateRNA<=(AnyPos_T1+30*select.max.months.after.AnyPos_T1))
	if(verbose) cat(paste("\nfound viral loads for candidate sequences with RNA<=max months after AnyPos_T1, n=",nrow(df.viro)))
	#for those for which there is no reliable seroNeg date, 
	#keep RNA for which there is at least one RNA after 'select.min.months.after.AnyPos_T1' mo after AnyPos_T1	
	if(0)
	{
		#TODO treat seroneg and seropos separately
		#TODO need not lump NegT_Acc="No" separately into seropos
		df.viro.seropos		<- subset(df.viro, is.na(NegT) | is.na(NegT_Acc) | NegT_Acc=="No" )							
		if(verbose) cat(paste("\nfound viral loads for candidate sequences with RNA<=max months after AnyPos_T1 for seropos, n=",nrow(df.viro.seropos)))
		df.viro.seropos		<- subset(df.viro.seropos, DateRNA>=(AnyPos_T1+30*select.min.months.after.AnyPos_T1))
		if(verbose) cat(paste("\nfound viral loads for candidate sequences with min months after AnyPos_T1<=RNA<=max months after AnyPos_T1 for seropos, n=",nrow(df.viro.seropos)))
	}
	else
	{
		df.viro				<- subset(df.viro, DateRNA>=(AnyPos_T1+30*select.min.months.after.AnyPos_T1))
		if(verbose) cat(paste("\nfound viral loads for candidate sequences with min months after AnyPos_T1<=RNA<=max months after AnyPos_T1, n=",nrow(df.viro)))
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
