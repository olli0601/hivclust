######################################################################################
project.hivc.check.DateRes.after.HIVPosTest<- function(dir.name= DATA, verbose=1)
{
	require(data.table)
	require(ggplot2)
	require(plyr)
	
	NL.HIV.phases<- as.Date( c("1980-01-01","1984-01-01","1996-01-01","2000-01-01","2004-01-01") )
	NL.possibly.Acute<- c(1,2)
	
	file	<- paste(dir.name,"tmp/ATHENA_2013_03_Sequences.R",sep='/')
	load(file)
	df.seq	<- df	
	file	<- paste(dir.name,"tmp/ATHENA_2013_03_Patients.R",sep='/')
	load(file)
	df.pat	<- df
	
	out		<- lapply( seq_along(df.seq),function(i)
			{
				cat(paste("\nprocess", names(df.seq)[i]))
				#no missing DateRes
				nok.idx<- which( is.na(df.seq[[i]][,"DateRes"]) )
				if(verbose)		cat(paste("\nentries in Sequences", nrow(df.seq[[i]])))
				if(verbose)		cat(paste("\nentries with missing DateRes", length(nok.idx)))
				seq.ok.idx<- which( !is.na(df.seq[[i]][,"DateRes"]) )
				df.seq[[i]]<- df.seq[[i]][seq.ok.idx,]								
				
				#extract min DateRes per patient
				df<- data.table( df.seq[[i]][,c("Patient","DateRes")], key="Patient" )
				df.minDateRes<- df[,min(DateRes),by=Patient]
				setnames(df.minDateRes,"V1","DateRes")
				#print(df.minDateRes)				
				if(verbose)		cat(paste("\npatients with DateRes", nrow(df.minDateRes)))
				
				#no missing MyDatePos1
				if(verbose)		cat(paste("\nentries in Patients", nrow(df.pat)))
				nok.idx<- which( is.na(df.pat[,"MyDatePos1"]) )				
				if(verbose)		cat(paste("\nentries with missing MyDatePos1", length(nok.idx)))				
				ok.idx<- which( !is.na(df.pat[,"MyDatePos1"]) )
				df.pos<- data.table( df.pat[ok.idx,c("Patient","MyDatePos1","MyDatePos1_Acc","AcuteInfection")], key="Patient" )				
				if(verbose)		cat(paste("\npatients with MyDatePos1", nrow(df.pos)))
				
				print(df.pos)
				#exclude Acute for checking HIVPos against DateRes
				df.pos<- subset(df.pos, !AcuteInfection%in%NL.possibly.Acute)
				if(verbose)		cat(paste("\npatients with MyDatePos1 and not known to acute or possibly acute", nrow(df.pos)))				
				
				#merge
				df<- merge(df.minDateRes,df.pos)
				if(verbose)		cat(paste("\npatients with HIVPos and DateRes", nrow(df)))
				#print(df)
				
				df.cmpHIVPos<- df[,difftime(MyDatePos1,DateRes,units="days")]
				df.cmpHIVPos<- as.numeric(df.cmpHIVPos)/12
				tmp<- which(df.cmpHIVPos>0)
				if(verbose)		cat(paste("\npatients with HIVPos>DateRes, n=", length(tmp)))	
				patient.laterHIVPos<- df[tmp,Patient]
				
				#in cases where MyDatePos1>DateRes, the guess of MyDatePos1 if MyDatePos1_Acc=0 may not be appropriate, improve this guess.								
				tmp<- which( df[,df.cmpHIVPos>0 & MyDatePos1_Acc==0 & as.POSIXlt(MyDatePos1)$mday==15] )
				if(verbose)		cat(paste("\npatients with unclear MyDatePos1 and day=15 and HIVPos>DateRes, n=", length(tmp)))				
				patient.betterHIVPos<- subset(df, df.cmpHIVPos>0 & MyDatePos1_Acc==0 & as.POSIXlt(MyDatePos1)$mday==15 & as.POSIXlt(MyDatePos1)$mon==as.POSIXlt(DateRes)$mon )
				if(verbose)		cat(paste("\npatients above for which HIVPos is in same month as DateRes (ie HIVPos can be fixed), n=", nrow(patient.betterHIVPos)))
				tmp<- subset( df, df.cmpHIVPos>0 & MyDatePos1_Acc==0 &  as.POSIXlt(MyDatePos1)$mday==1 & as.POSIXlt(MyDatePos1)$mon==6 )
				if(verbose)		cat(paste("\npatients with unclear MyDatePos1 and day=1, month=7 and HIVPos>DateRes, n=", nrow(tmp)))
				tmp<- subset( df, df.cmpHIVPos>0 & MyDatePos1_Acc==0 &  as.POSIXlt(MyDatePos1)$mday==1 & as.POSIXlt(MyDatePos1)$mon==6 & as.POSIXlt(MyDatePos1)$year==as.POSIXlt(DateRes)$year )
				if(verbose)		cat(paste("\npatients above for which HIVPos is in same year as DateRes (ie HIVPos can be fixed), n=", nrow(tmp)))				
				patient.betterHIVPos<- rbind(patient.betterHIVPos, tmp)
				setkey(patient.betterHIVPos, Patient)						#important to sort for next step, which assumes sorted
				tmp<- df.pos[,MyDatePos1]
				tmp[ which( df.pos[,Patient%in%patient.betterHIVPos$Patient]) ]<- patient.betterHIVPos$DateRes
				df.pos[,orDatePos1:=tmp]
				
				#re-merge
				df<- merge(df.minDateRes,df.pos)
				df.cmpHIVPos<- df[,difftime(orDatePos1,DateRes,units="days")]
				df.cmpHIVPos<- as.numeric(df.cmpHIVPos)/12
				tmp<- which(df.cmpHIVPos>0)
				if(verbose)		cat(paste("\npatients with HIVPos>DateRes after FIX, n=", length(tmp)))
				hist(df.cmpHIVPos[tmp], breaks=100)
				print( df[tmp,] )
				
				#tmp<- which(df.cmpHIVPos>1)
				#if(verbose)		cat(paste("\npatients with HIVPos>DateRes after FIX and 1 month grace, n=", length(tmp)))
				#print( df[tmp,] )
				
				patient.laterHIVPos<- df[tmp,Patient]				
				#if(verbose)		cat(paste("\npatients with HIVPos>DateRes, Patient", paste(df[tmp,Patient],collapse=', ') ))
				#print(df[tmp,])
				#hist(df.cmpHIVPos, breaks=50)
				
				if(0)
				{
					df.earlierHIVPos<- -df.cmpHIVPos[ which(df.cmpHIVPos<=0) ]
					df.earlierHIVPos[ df.earlierHIVPos < 1 ]<- 1		#if sequenced within 1mo after diagnosis, ignore
					df.earlierHIVPos<- sort(df.earlierHIVPos)
					col<- "grey70"
					x<- df.earlierHIVPos
					y<- seq_along(df.earlierHIVPos) / length(df.earlierHIVPos)				
					plot(1,1,type='n',bty='n',xlab="time [months]", ylab="%sequenced after diagnosis", log='x', xlim=range(x),ylim=range(y))				
					apply( matrix(seq(1,max(x),by=6),2) ,2,function(x)
							{
								polygon(c(x,rev(x)),c(0,0,1,1),border=NA,col=col)
							})
					lines(x,y,type='s')
				}
				if(1)
				{											
					tmp		<- numeric(nrow(df))
					tmp[ which( df[,orDatePos1]<NL.HIV.phases[2] ) ]<- 1
					tmp[ which( df[,orDatePos1]>=NL.HIV.phases[2] & df[,orDatePos1]<NL.HIV.phases[3]) ]<- 2
					tmp[ which( df[,orDatePos1]>=NL.HIV.phases[3] & df[,orDatePos1]<NL.HIV.phases[4]) ]<- 3
					tmp[ which( df[,orDatePos1]>=NL.HIV.phases[4] & df[,orDatePos1]<NL.HIV.phases[5]) ]<- 4
					tmp[ which( df[,orDatePos1]>NL.HIV.phases[5] ) ]<- 5					
					tmp		<- data.frame( HIVPos=df[,orDatePos1], TimeToSeq= -df.cmpHIVPos, HIVPhase=tmp )
					tmp		<- tmp[ tmp$TimeToSeq>=0, ]
					tmp$TimeToSeq[tmp$TimeToSeq<1]<- 1
					#cumsum(phase) / nrow
					#print( cumsum( as.numeric( tmp$HIVPhase==1 ) ) )
					#quit("no")
					tmp		<- cbind( tmp[ with(tmp, order(TimeToSeq)), ], seq_len(nrow(tmp)) )
					counts	<- t(ddply(tmp, .(HIVPhase), function(x)
							{								
								x<- x[with(x,order(TimeToSeq)),]
								z<- numeric(nrow(tmp))
								for(j in seq_len(nrow(x)))
								{
									z[x[j,4]:length(z)]<- 1+z[x[j,4]:length(z)]
								}
								z
							}))
					counts<- counts[-1,]		
					tmp<- as.data.frame( cbind(tmp$TimeToSeq,counts) )
					colnames(tmp)<- c("TimeToSeq",paste("ph",1:5,sep=''))
					
					cols	<- rainbow(5)
					z		<- numeric(nrow(tmp))
					xlim	<- c(1,max(tmp$TimeToSeq))
					plot(1,1,type='n',xlab="Time from Diag To Seq [months]",ylab="number seq",xlim=xlim,ylim=c(0,nrow(df)), log='x')					
					for(j in 2:ncol(tmp))
							{
								polygon( c(tmp$TimeToSeq,rev(tmp$TimeToSeq)), c(z,rev(z+tmp[,j])), border=NA, col=cols[j-1] )
								z<- z+tmp[,j]
							}
					legend("topleft",fill=cols,legend=colnames(tmp)[2:6],bty='n', border=NA)		
					#need different data frame see  http://stackoverflow.com/questions/5030389/getting-a-stacked-area-plot-in-r 
					#p<- 	ggplot(data=tmp,aes(x=TimeToSeq,y=val, group=HIVPhase, colour= HIVPhase))
					#p<- 	p+geom_line(aes(col=HIVPhase))
					#p<-		p+geom_area(position = "fill")	
					#print(p)					
					quit("no")
				}												
				list(patient.laterHIVPos=patient.laterHIVPos)
			})
	patient.laterHIVPos	<- unique( unlist( lapply(out,function(x) x[["patient.laterHIVPos"]]) ) )
	if(verbose)		cat(paste("\npatients with HIVPos>DateRes RT or PROT, n=", length(patient.laterHIVPos)))
	if(verbose)		cat(paste("\npatients with HIVPos>DateRes RT or PROT, Patient", paste(patient.laterHIVPos,collapse=', ') ))		
}
######################################################################################
project.hivc.check.DateRes.after.T0<- function(dir.name= DATA, verbose=1)
{
	require(data.table)
	
	file<- paste(dir.name,"tmp/ATHENA_2013_03_Sequences.R",sep='/')
	load(file)
	df.seq<- df

	file<- paste(dir.name,"tmp/ATHENA_2013_03_Regimens.R",sep='/')
	load(file)
	df.reg<- df
	
	lapply( seq_along(df.seq),function(i)
			{
				cat(paste("\nprocess", names(df.seq)[i]))
				#no missing DateRes
				nok.idx<- which( is.na(df.seq[[i]][,"DateRes"]) )
				if(verbose)		cat(paste("\nentries in Sequences", nrow(df.seq[[i]])))
				if(verbose)		cat(paste("\nentries with missing DateRes", length(nok.idx)))
				seq.ok.idx<- which( !is.na(df.seq[[i]][,"DateRes"]) )
				df.seq[[i]]<- df.seq[[i]][seq.ok.idx,]
				#print( range(df.seq[[i]][,"DateRes"]) )
				
				
				df<- data.table( df.seq[[i]][,c("Patient","DateRes")], key="Patient" )
				df.minDateRes<- df[,min(DateRes),by=Patient]
				setnames(df.minDateRes,"V1","DateRes")
				print(df.minDateRes)
				#print( range(df.minDateRes[,DateRes]) )
				if(verbose)		cat(paste("\npatients with DateRes", nrow(df.minDateRes)))
				
				#no missing T0				
				nok.idx<- which( is.na(df.reg[,"T0"]) )
				if(verbose)		cat(paste("\nentries in Regimens", nrow(df.reg)))
				if(verbose)		cat(paste("\nentries with missing T0", length(nok.idx)))
				reg.ok.idx<- which( !is.na(df.reg[,"T0"]) )
				df.reg<- df.reg[reg.ok.idx,]
				
				#extract min T0 per patient
				df<- data.table( df.reg[,c("Patient","T0")], key="Patient" )
				df.minT0<- df[,min(T0),by=Patient]
				setnames(df.minT0,"V1","T0")
				if(verbose)		cat(paste("\npatients with T0", nrow(df.minT0)))
				print(df.minT0)
				
				#merge
				df<- merge(df.minDateRes,df.minT0)
				if(verbose)		cat(paste("\npatients with T0 and DateRes", nrow(df)))
				
				#print(df[6258:6262,])
				df.laterT0<- df[,difftime(T0,DateRes,units="days")]
				df.laterT0<- as.numeric(df.laterT0)/365.25
				tmp<- which(df.laterT0>0)
				if(verbose)		cat(paste("\npatients with T0>DateRes, n=", length(tmp)))				
				#if(verbose)		cat(paste("\npatients with T0>DateRes, Patient", paste(df[tmp,Patient],collapse=', ') ))
				#print(df[tmp,])
				hist(df.laterT0, breaks=50)
				#print(range(df.laterT0))				
			})
}
######################################################################################
project.hivc.check<- function()
{
	if(1) project.hivc.check.DateRes.after.HIVPosTest()
	if(0) project.hivc.check.DateRes.after.T0()	
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
project.hivc.collectpatientdata<- function(dir.name= DATA, verbose=1, resume=0)
{	
	require(data.table)		
	
	#input files generated with "project.hivc.Excel2dataframe"
	resume			<- 0
	verbose			<- 1
	file.seq		<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.R",sep='/')
	file.patient	<- paste(dir.name,"derived/ATHENA_2013_03_Patients.R",sep='/')
	file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	file.immu		<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
	file.treatment	<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	
	#compute file AllSeqPatientCovariates
	file.out.name	<- "ATHENA_2013_03_AllSeqPatientCovariates"
	file.out		<- paste(dir.name,"/derived/",file.out.name,".R",sep='')
	if(resume)												#//load if there is R Master data.table
	{
		options(show.error.messages = FALSE)		
		readAttempt<-try(suppressWarnings(load(file.out)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresume file",file.out))
		if(!inherits(readAttempt, "try-error"))	str(df.all)
		options(show.error.messages = TRUE)		
	}
	if(!resume || inherits(readAttempt, "try-error"))		#else generate R Master data.table
	{
		#get data.table of seqs		
		load(file.seq)			
		df		<- lapply( df,function(x)	data.table(x[,c("Patient","DateRes","SampleCode")], key="SampleCode")	)
		df		<- merge(df[[1]],df[[2]],all.x=1,all.y=1)
		setkey(df, Patient.x)
		tmp		<- which(df[,is.na(Patient.x)])
		set(df, tmp, "Patient.x", df[tmp, Patient.y])
		set(df, tmp, "DateRes.x", df[tmp, DateRes.y])
		set(df, NULL, "SampleCode", gsub(' ','', df[, SampleCode]))
		setnames(df, c("SampleCode","Patient.x","DateRes.x"), c("FASTASampleCode","Patient","PosSeqT"))
		set(df, NULL, "FASTASampleCode", factor(df[,FASTASampleCode]))
		set(df, NULL, "Patient", factor(df[,Patient]))		
		df.all	<- subset(df, select=c(FASTASampleCode,Patient,PosSeqT))
		if(verbose)		cat(paste("\nnumber of sequences found, n=", nrow(df.all)))
		#
		# add Patient data
		#
		if(verbose)		cat(paste("\nadding patient data"))
		load(file.patient)		
		df.all	<- merge(df.all, df, all.x=1, by="Patient")
		setnames(df.all, c("MyDateNeg1","MyDatePos1"), c("NegT","PosT"))
		#	reset PosT_Acc=='No' conservatively to end of year / month
		df.all	<- hivc.db.reset.inaccuratePosT(df.all, nacc.dy.dy= 30, nacc.mody.mo= 11, nacc.mody.dy= 31, verbose=1)
		#	reset NegT_Acc=='No' conservatively to start of year / month
		df.all	<- hivc.db.reset.inaccurateNegT(df.all)
		#	check for clashes in NegT and PosT
		tmp		<- subset(df.all, !is.na(PosT) & !is.na(NegT) & PosT<NegT)
		if(verbose)		cat(paste("\nchecking for clashes in NegT and PosT"))
		if(verbose)		print(tmp)
		if(verbose)		cat(paste("\nmanually resolving clashes -- M41654 has wrong PosT since PosRNA much later -- M12982 has wrong NegT as PosRNA before that time"))
		#	resolving M41654
		tmp		<- which(df.all[,FASTASampleCode=="M4165429052012"])
		set(df.all, tmp, "PosT", df.all[tmp, PosSeqT]) 		
		# 	add preliminary AnyPos_T1	-- since PosT conservative we are good to set irrespective of PosT_Acc	
		df.all[, AnyPos_T1:=PosSeqT]
		tmp		<- which( df.all[, !is.na(PosT) & ( is.na(PosSeqT) | PosT<PosSeqT ) ] )
		if(verbose)		cat(paste("\nbuilding prel AnyPos_T1. Number of seq with !is.na(PosT) & ( is.na(PosSeqT) | PosT<PosSeqT ), n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, PosT])		
		if(nrow(subset(df.all, is.na(AnyPos_T1))))	stop("unexpected NA in AnyPos_T1")
		#	check for invalid NegT and set to NA	-- we would only know that NegT is before AnyPosT and this is not helpful
		df.all	<- hivc.db.resetNegTbyAnyPosT(df.all)				
		#
		#	add first RNA Virology date
		#
		if(verbose)		cat(paste("\nadding virology data"))
		load(file.viro)
		df.cross	<- merge( subset(df.all, select=c(FASTASampleCode,Patient,AnyPos_T1,PosT,PosSeqT,NegT,NegT_Acc)), df, allow.cartesian=T, by="Patient" )
		#	checking manually AnyPos_T1>PosRNA & lRNA<3
		if(verbose)		cat(paste("\ncheck manually AnyPos_T1>PosRNA & lRNA<3 -- THIS IS ASSUMED OK including 2004G180"))
		tmp	<- subset(df.cross, AnyPos_T1>PosRNA & lRNA<3)
		tmp[,diff:=tmp[, difftime(PosRNA,AnyPos_T1, units="weeks")]]
		print( subset(tmp, diff< -10) )		
		#subset(df.cross, FASTASampleCode=="2004G180")
		#		
		#	checking manually NegT>PosRNA
		#
		if(verbose)		cat(paste("\ncheck manually NegT_Acc=='Yes' & NegT>PosRNA -- THIS IS ASSUMED OK"))
		print( subset(df.cross, NegT_Acc=="Yes" & NegT>PosRNA) )
		if(verbose)		cat(paste("\ncheck manually NegT_Acc=='No' & NegT>PosRNA -- THIS IS ASSUMED OK"))
		print( subset(df.cross, NegT_Acc=="No" & NegT>PosRNA) )
		# 	compute lRNA_T1 and lRNA_TS
		tmp			<- hivc.db.getlRNA.T1andTS(df.cross, lRNA.bTS.quantile= 0.75, lRNA.aTS.quantile= 0.25, lRNAi.min= log10(1e4), verbose=1)
		df.all		<- merge(df.all, tmp, all.x=1, by="FASTASampleCode")
		#	reset NegT by lRNA_T1
		df.all<- hivc.db.resetNegTbyPoslRNA_T1(df.all)
		# 	reset preliminary AnyPos_T1
		tmp		<- df.all[, list(AnyPos_T1=min(AnyPos_T1,PoslRNA_T1, na.rm=1)), by="FASTASampleCode"]
		if(verbose)	cat(paste("\nnumber of seq with PoslRNA_T1<AnyPos_T1, n=",length(tmp),"RESETTING"))
		set(df.all,NULL,"AnyPos_T1",tmp[,AnyPos_T1])		
		#
		#	add CD4 count data
		#
		if(verbose)		cat(paste("\nadding CD4 data"))
		load(file.immu)
		if(length(which(df[,is.na(PosCD4_T1)])))	stop("unexpected NA in PosCD4_T1")
		if(length(which(df[, is.na(PosCD4)]))) stop("unexpected NA in PosCD4")
		if(length(which(df[, is.na(CD4)]))) stop("unexpected NA in CD4")
		df.cross	<- merge( subset(df.all, select=c(FASTASampleCode,Patient,AnyPos_T1,PosT,PosSeqT,PoslRNA_T1,NegT,NegT_Acc)), subset(df, select=c(Patient,PosCD4,CD4)), allow.cartesian=T, by="Patient" )
		# delete entries where NegT>PosCD4
		df.cross	<- hivc.db.resetCD4byNegT(df.cross, with.NegT_Acc.No=1, verbose=1)
		# compute CD4_T1 and CD4_TS -- do not use on df[,CD4_T1] because some CD4 measurements might correspond to healthy patients
		tmp		<- hivc.db.getCD4.T1andTS(df.cross, CD4.HIVNeg.min= 500)
		df.all	<- merge(df.all, tmp, by="FASTASampleCode", all.x=1)
		#
		# manually checked remaining PosCD4_T1 < AnyPos_T1 -- overall plausible
		#
		#	tmp		<- subset(df.all,!is.na(PosCD4_T1) & difftime(PosCD4_T1,AnyPos_T1, units="weeks")< 0)
		#	tmp[, diff:=as.numeric(tmp[,difftime(PosCD4_T1,AnyPos_T1, units="weeks")])]
		# 	subset(tmp,diff< -10)
		#
		if(verbose)		cat(paste("\ncheck manually PosCD4_T1 < AnyPos_T1 -- THIS IS ASSUMED OK"))
		tmp	<- subset(df.all, PosCD4_T1 < AnyPos_T1, c(FASTASampleCode, Patient, AnyPos_T1, PosSeqT, PosT, PosT_Acc, PoslRNA_T1, lRNA_T1,  PosCD4_T1, CD4_T1))
		print( tmp , nrow=400)
		tmp		<- which(df.all[,!is.na(PosCD4_T1) & PosCD4_T1<AnyPos_T1])
		if(verbose)		cat(paste("\nnumber of seq with !is.na(PosCD4_T1) & PosCD4_T1<AnyPos_T1, n=",length(tmp)))
		if(verbose)		cat(paste("\nnumber of patients with !is.na(PosCD4_T1) & PosCD4_T1<AnyPos_T1, n=",length(unique(df.all[tmp,Patient]))))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp,PosCD4_T1])
		#
		#	add Treatment dates
		#
		if(verbose)		cat(paste("\nadding treatment data"))
		load(file.treatment)
		df			<- subset(df, select=c(Patient, StartTime, StopTime, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel, TrI.n, TrI.mo, TrI.p, AnyT_T1, AnyT_T1_Acc, HAART_T1))
		tmp			<- subset(df, select=c(Patient, TrI.n, AnyT_T1, AnyT_T1_Acc ))
		setkey(tmp, Patient)
		df.all		<- merge(df.all, unique(tmp), all.x=1, by="Patient")		
		#	compare treatment history relative to PosSeqT		
		df.cross	<- merge( subset(df.all, select=c(FASTASampleCode,Patient,PosSeqT)), df, allow.cartesian=T, by="Patient" )
		tmp			<- hivc.db.getTrIMo(df.cross)
		df.all		<- merge(df.all, tmp, all.x=1, by="FASTASampleCode")
		#		
		setkey(df.all, PosSeqT)
		setkey(df.all, Patient)
		#		
		#	manually checked remaining AnyT_T1<AnyPos_T1 -- overall plausible
		#
		if(verbose)		cat(paste("\ncheck manually AnyT_T1<AnyPos_T1 -- THIS IS ASSUMED OK"))
		tmp	<- subset(df.all, !is.na(AnyT_T1) & AnyT_T1<AnyPos_T1, c(Patient, FASTASampleCode, PosCD4_T1, CD4_T1, TrI.n, PosSeqT, AnyPos_T1,  AnyT_T1, AnyT_T1_Acc))
		tmp[, diff:=as.numeric(tmp[,difftime(AnyT_T1,AnyPos_T1, units="weeks")])]
		subset(tmp,diff< -10)
		tmp		<- which( df.all[,!is.na(AnyT_T1) & AnyT_T1<AnyPos_T1] )
		if(verbose)		cat(paste("\nnumber of seq with !is.na(AnyT_T1) & AnyT_T1<AnyPos_T1, n=",length(tmp),"SET AnyPos_T1 to lower value"))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp,AnyT_T1])
		df.all[, AnyT_T1_Acc:=NULL]
		#	round numbers
		#set(df.all, NULL, "TrI.p", round(df.all[,TrI.p],d=2))
		#set(df.all, NULL, "TrI.mo", round(df.all[,TrI.mo],d=1))
		if(verbose)	cat(paste("\nsave to file",file.out))
		#setkey(df.all, FASTASampleCode)
		save(df.all,file=file.out)
		str(df.all)		
	}	
	#
	#	new diagnoses by CD4
	#
	df.newdiag		<- copy(subset(df.all, select=c(Patient,AnyPos_T1, CD4_T1)))
	setkey(df.newdiag,Patient)
	df.newdiag		<- unique(df.newdiag)
	df.newdiagCD4 	<- hivc.db.getplot.newdiagnosesbyCD4(df.newdiag, plot.file= paste(dir.name,"/derived/",file.out.name,"_NewDiagByCD4.pdf",sep=''), plot.ylab= "New diagnoses HIV-1 subtype B,\n seq available")
	#
	#	seem in care by risk group
	#
	df.living		<- copy(subset(df.all, select=c(Patient, AnyPos_T1, Trm, DateDied, DateLastContact)))
	setkey(df.living,Patient)
	df.living		<- unique(df.living)
	tmp				<- hivc.db.getplot.livingbyexposure(df.living, plot.file=paste(dir.name,"/derived/",file.out.name,"_Seen4CareByExpGroup.pdf",sep=''), plot.ylab="Seen for care with HIV-1 subtype B,\n seq available", db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)
	
	quit("no")
	#
	#compute coverage of all data covariates
	#
	#coverage by seq
	df.covbyseq	<- hivc.db.getcoverage(df)	
	#coverage by patient
	df			<- subset(df.all, select= c(	Patient, DateBorn, Sex, CountryBorn, RegionOrigin, DateDied, Subtype, isAcute, 
												NegT, NegT_Acc, PosT, PosT_Acc, CountryInfection, Trm,  DateLastContact, RegionHospital, DateFirstEverCDCC,
												isDead, PoslRNA_T1, lRNA_T1, PosCD4_T1, CD4_T1, AnyT_T1))
	setkey(df, Patient)
	df			<- unique(df)
	df.covbypat	<- hivc.db.getcoverage(df)							
	
	quit("no")
	#
	#compute file All1stPatientCovariates
	#
	file.out		<- paste(dir.name,"derived/ATHENA_2013_03_All1stPatientCovariates.R",sep='/')
	if(resume)												#//load if there is R Master data.table
	{
		options(show.error.messages = FALSE)		
		readAttempt<-try(suppressWarnings(load(file.out)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresume file",file.out))
		if(!inherits(readAttempt, "try-error"))	str(df.all)
		options(show.error.messages = TRUE)		
	}
	if(!resume || inherits(readAttempt, "try-error"))		#else generate R Master data.table
	{
		#get PosSeqPROT   PosSeqRT		
		load(file.seq)
		df.seq<- df	
		df.minDateRes	<- lapply( seq_along(df.seq),function(i)
				{
					cat(paste("\nprocess", names(df.seq)[i]))
					#extract entries without missing DateRes
					df.miss<- data.table(subset(df.seq[[i]], is.na(DateRes), c(Patient,DateRes,SampleCode)), key="Patient")
					if( nrow(df.miss)!=length(unique(df.miss[,Patient])) ) stop("handling of missing DateRes not appropriate")
					if(verbose)		cat(paste("\npatients without DateRes", nrow(df.miss)))
					
					df<- subset(df.seq[[i]], !is.na(DateRes), c(Patient,DateRes,SampleCode))
					df<- data.table(df, key="Patient")								
					if(verbose)		cat(paste("\nentries without missing DateRes", nrow(df)))							
					#extract min DateRes per patient	
					df<- df[,.SD[which.min(DateRes)],by=Patient]
					if(verbose)		cat(paste("\npatients with DateRes", nrow(df)))					
					df<- rbindlist(list(df,subset( df.miss, !Patient%in%df$Patient )))
					setkey(df, Patient)
					df[, "SampleCode"]	<- factor(df[,SampleCode])
					setnames(df,"SampleCode",paste("Seq",names(df.seq)[i],sep=''))
					setnames(df,"DateRes",paste("Pos",names(df.seq)[i],sep=''))							
					df
				})
		df.minDateRes<- merge(df.minDateRes[[1]],df.minDateRes[[2]],all=1)
		if(verbose)		cat(paste("\npatients RT or PROT", nrow(df.minDateRes)))
		if(nrow(df.minDateRes)!=length(unique( df.minDateRes[,Patient] )))	stop("non-unique patients at this point")		
		#str(df.minDateRes)
	
		#add Patient data		
		load(file.patient)		
		df.all<- df[df.minDateRes]
		if(verbose)		cat(paste("\npatients in combined data table", nrow(df.all)))
		
		#add Treatment dates
		load(file.treatment)
		df		<- data.table(df, key="Patient")
		setnames(df, "T0","HAART_T1")
		df		<- subset(df, select= c(Patient, HAART_T1, StartTime))		
		if(verbose)		cat(paste("\n\nadding treatment data\nentries in regimens data table", nrow(df)))
		df		<- subset(df, !is.na(StartTime) )				
		df		<- df[, { tmp<- which.min(StartTime); list(HAART_T1=HAART_T1[tmp], AnyT_T1=StartTime[tmp])}, by=Patient]		
		if( nrow(subset(df,!is.na(HAART_T1) & HAART_T1<AnyT_T1)) )	stop("found AnyT_T1 that is older than HAART_T1")
		if(verbose)		cat(paste("\npatients with at least one non-missing treatment date", nrow(df)))
		df		<- subset(df,select=c(Patient,AnyT_T1))
		setnames(df, "AnyT_T1","AnyTherT1")		
		df.all	<- df[df.all]
		#setnames(df, "T0","PREHAART_T1")
		
		#add first RNA Virology date		
		load(file.viro)
		df		<- data.table(df, key="Patient")		
		#str(df)		
		setnames(df, "DateRNA","PosRNA")
		if(verbose)		cat(paste("\n\nadding virology\nentries in viro data table", nrow(df)))
		df		<- subset(df, !is.na(PosRNA) & Undetectable!="Yes" )
		if(verbose)		cat(paste("\nentries in viro data table without missing or non-detectable", nrow(df)))
		df		<- df[,{tmp<- which.min(PosRNA); list(PosRNA_T1= PosRNA[tmp], RNA_T1= RNA[tmp]) }, by=Patient]		
		if(verbose)		cat(paste("\npatients with at least one non-missing and detectable RNA date", nrow(df)))
		if(0)
		{
			print(range(df[,PosRNA], na.rm=1))		
		}
		df.all<- df[df.all]
		
		#add first CD4 count date		
		load(file.immu)
		df		<- data.table(df, key="Patient")
		setnames(df, "DateImm","PosCD4")
		if(verbose)		cat(paste("\nadding immunology\nentries in immu data table", nrow(df)))
		df<- subset(df, !is.na(PosCD4) )
		if(verbose)		cat(paste("\nentries in immu data table, non-missing", nrow(df)))
		df		<- df[,{tmp<- which.min(PosCD4); list(PosCD4_T1= PosCD4[tmp], CD4_T1= CD4A[tmp]) }, by=Patient]
		if(verbose)		cat(paste("\npatients with at least one non-missing CD4 date", nrow(df)))
		if(0)
		{
			print(range(df[,PosCD4_T1], na.rm=1))		
		}
		df.all<- df[df.all]
		df.all<- subset(df.all,select=c(Patient,Trm,SeqPROT,SeqRT,PosPROT,PosRT,PosT,PosT_Acc,PosCD4_T1,PosRNA_T1, NegT,NegT_Acc, AnyTherT1, isAcute,isDead,Died, RNA_T1, CD4_T1))
		
		
		#if 	PosPROT!=PosRT		consider only earliest sequence and set other to missing
		tmp<- subset(df.all, is.na(PosPROT) & is.na(PosRT))
		if(verbose)		cat(paste("\n\nkeep only earliest sequence PROT or RT if PosPROT!=PosRT\npatients with is.na(PosPROT) & is.na(PosRT):", nrow(tmp)))				
		tmp<- subset(df.all, !is.na(PosPROT) & !is.na(PosRT) & PosPROT!=PosRT)
		if(verbose)		cat(paste("\npatients with PosPROT!=PosRT & !is.na(PosPROT) & !is.na(PosRT):", nrow(tmp)))
		tmp<- subset(df.all, !is.na(PosPROT) & !is.na(PosRT) & PosPROT<PosRT, Patient)
		df.all[tmp,c("PosRT","SeqRT")]<- as.factor(NA)
		tmp<- subset(df.all, !is.na(PosPROT) & !is.na(PosRT) & PosPROT>PosRT, Patient)
		df.all[tmp,c("PosPROT","SeqPROT")]<- as.factor(NA)	
		tmp<- subset(df.all, PosPROT!=PosRT & !is.na(PosPROT) & !is.na(PosRT))
		if(verbose)		cat(paste("\npatients with PosPROT!=PosRT & !is.na(PosPROT) & !is.na(PosRT)	AFTER FIX:", nrow(tmp)))		
		#replace PosPROT, PosRT with PosSeq
		tmp					<- df.all[,PosPROT]
		tmp[ is.na(tmp) ]	<- df.all[is.na(tmp),PosRT]	
		df.all[,PosSeq:=tmp]		
		df.all<- subset(df.all, select=c(Patient,Trm,SeqPROT,SeqRT,PosSeq,PosT,PosT_Acc,PosCD4_T1,PosRNA_T1,NegT,NegT_Acc, AnyTherT1, isAcute,isDead,Died,RNA_T1, CD4_T1))
		
		
		if(verbose)		cat(paste("\nsave df.all to", file.out))				
		save(df.all,file=file.out)
		str(df.all)
	}
}
######################################################################################
project.hivc.Excel2dataframe.Regimen<- function(dir.name= DATA, verbose=1)
{
	NA.time			<- c("01/01/1911","01/11/1911","11/11/1911")	
	MAX.time		<- c("")
	TR.notyet		<- "30/03/2013"
	TR.failure 		<- c(21,31,32)
	TR.adherence	<- c(47)
	TR.patrel		<- c(23, 24, 33, 34, 36)	#either patient s decision, toxicity or interaction with other medication
	#read REGIMEN csv data file and preprocess
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.csv",sep='/')
	df				<- read.csv(file, stringsAsFactors=FALSE)							
	
	date.var		<- c("T0","StartTime","StopTime")		
	for(x in date.var)
	{
		cat(paste("\nprocess Time", x))
		nok.idx			<- which( df[,x]==NA.time[1] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[1],", n=", length(nok.idx), "set to NA"))
		if(length(nok.idx))	
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[2] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[2],", n=", length(nok.idx), "set to NA"))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[3] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[3],", n=", length(nok.idx), "set to NA"))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA			
		nok.idx			<- which( df[,x]==MAX.time[1] )
		if(verbose)	cat(paste("\nentries with format ",MAX.time[1],", n=", length(nok.idx), "set to max time 01/01/2999"))
		if(length(nok.idx))
			df[nok.idx,x]	<- TR.notyet
		df[,x]			<- as.Date(df[,x], format="%d/%m/%Y")	
	}
	TR.notyet				<- as.Date(TR.notyet, format="%d/%m/%Y")
	
	df						<- data.table(df, key="Patient")
	setnames(df, "T0","HAART_T1")
	set(df,NULL,"Patient",factor(df[,Patient]))
	if(verbose)	cat(paste("\nnumber of entries, n=",nrow(df)))
	#	fix entry 	M17493  1996-01-01 manually
	tmp						<- which(df[, Patient=="M17493" & StartTime=="1996-01-01"])
	if(verbose)	cat(paste("\nfix entry 	M17493  1996-01-01 manually to somewhere in 1996"))
	set(df, tmp, "StartTime", as.Date("1996-07-01"))		
	#	fix entry 	M41688  2011-01-01 manually
	tmp						<- which(df[, Patient=="M41688" & StartTime=="2011-01-01"])
	if(verbose)	cat(paste("\nfix entry 	M41688  2011-01-01 manually to 2011-11-01 (Ard)"))
	set(df, tmp, "StartTime", as.Date("2011-11-01"))	
	set(df, which(df[, Patient=="M41688"]), "HAART_T1", as.Date("2011-11-01"))		
	#	fix entry 	M41688  2011-08-23 manually
	tmp						<- which(df[, Patient=="M42092" & StartTime=="2011-08-23"])
	if(verbose)	cat(paste("\nfix entry 	M42092  2011-08-23 to 2012-08-23"))
	set(df, tmp, "StartTime", as.Date("2012-08-23"))
	set(df, tmp, "HAART_T1", as.Date("2012-08-23"))
	#	fix entry 	M41688  2010-10-23 manually
	tmp						<- which(df[, Patient=="M42186" & StartTime=="2010-10-23"])
	if(verbose)	cat(paste("\nfix entry 	M42186  2010-10-23 to 2012-10-23"))
	set(df, tmp, "StartTime", as.Date("2012-10-23"))
	set(df, tmp, "HAART_T1", as.Date("2012-10-23"))		
	#	fix entry 	M37531  2006-08-15 manually
	tmp						<- which(df[, Patient=="M37531" & StartTime=="2006-08-15"])
	if(verbose)	cat(paste("\nfix entry 	M37531  2006-08-15 to 2008-08-15"))
	set(df, tmp, "StartTime", as.Date("2008-08-15"))
	tmp						<- which(df[, Patient=="M37531"])
	set(df, tmp, "HAART_T1", as.Date("2008-08-15"))				
	#	set potentially inaccurate StartTime	
	tmp																										<- rep(0, nrow(df))
	tmp[is.na(df[,StartTime])]																				<- NA
	tmp[which( df[, !is.na(StartTime) & as.POSIXlt(StartTime)$mday==15] )]									<- 1
	tmp[which( df[, !is.na(StartTime) & as.POSIXlt(StartTime)$mon==6 & as.POSIXlt(StartTime)$mday==1 ] )]	<- 2
	if(verbose) cat(paste("\nnumber of inaccurate StartTime entries, n=",length(tmp[na.omit(tmp>0)])))
	df[,StartTime_Acc:= factor(tmp,levels=c(0,1,2),labels=c("Acc","NAccD","NAccMD"))]
	#	set potentially inaccurate StopTime
	tmp																										<- rep(0, nrow(df))
	tmp[is.na(df[,StopTime])]																				<- NA
	tmp[which( df[, !is.na(StopTime) & as.POSIXlt(StopTime)$mday==15] )]									<- 1
	tmp[which( df[, !is.na(StopTime) & as.POSIXlt(StopTime)$mon==6 & as.POSIXlt(StopTime)$mday==1 ] )]		<- 2
	if(verbose) cat(paste("\nnumber of inaccurate StopTime entries, n=",length(tmp[na.omit(tmp>0)])))
	df[,StopTime_Acc:= factor(tmp,levels=c(0,1,2),labels=c("Acc","NAccD","NAccMD"))]
	#
	#	removing Patients not on treatment
	#
	tmp						<- which(df[, is.na(StartTime) & StopTime==TR.notyet & NoDrug==0])
	if(verbose)	cat(paste("\nnumber of entries with is.na(StartTime) & StopTime==TR.notyet & NoDrug==0",length(tmp),"SETTING TO NA"))
	set(df, tmp, "StopTime", NA)
	tmp						<- which(df[, StartTime==TR.notyet & StopTime==TR.notyet])
	if(verbose)	cat(paste("\nnumber of entries with StartTime==TR.notyet & StopTime==TR.notyet",length(tmp),"SETTING TO NA"))
	set(df, tmp, "StartTime", NA)
	set(df, tmp, "StopTime", NA)
	tmp						<- which(df[, StartTime==TR.notyet])
	if(verbose)	cat(paste("\nnumber of entries with StartTime==TR.notyet & StopTime!=TR.notyet",length(tmp),"MISCLASSIFIED StartTime - setting to NA"))
	set(df, tmp, "StartTime", NA)		
	df						<- subset(df, !is.na(StopTime))
	if(verbose)	cat(paste("\nnumber of entries with !is.na(StartTime) & !is.na(StopTime), n=",nrow(df)))		
	#
	#	TR.interrupted
	#
	if(nrow(df[which(is.na(df[,NoDrug]) & StartTime!=TR.notyet),])) stop("unexpected NA in NoDrug when on treatment")
	tmp								<- rep(0, nrow(df))
	tmp[ which( df[,NoDrug==0] ) ]	<- 1
	df[, TrI:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#
	#	process treatment change reasons
	#
	if(any(!is.na(df[, Reason7])))	stop("unexpected !NA after Reason7")
	df.TrCh.noreason		<- which( df[, is.na(Reason1)&is.na(Reason2)&is.na(Reason3)&is.na(Reason4)&is.na(Reason5)&is.na(Reason6)] )		
	#	TR.failure
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.failure | Reason2%in%TR.failure | Reason3%in%TR.failure | Reason4%in%TR.failure | Reason5%in%TR.failure | Reason6%in%TR.failure ]) ]<- 1
	df[, TrCh.failure:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#	TR.adherence
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.adherence | Reason2%in%TR.adherence | Reason3%in%TR.adherence | Reason4%in%TR.adherence | Reason5%in%TR.adherence | Reason6%in%TR.adherence ]) ]<- 1
	df[, TrCh.adherence:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#	TR.patient related
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.patrel | Reason2%in%TR.patrel | Reason3%in%TR.patrel | Reason4%in%TR.patrel | Reason5%in%TR.patrel | Reason6%in%TR.patrel ]) ]<- 1
	df[, TrCh.patrel:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#
	#	removing Patient entries before treatment started
	#
	tmp						<- which( df[,is.na(StartTime) & TrI=="Yes"]  )
	tmp						<- df[tmp,][, df.idx:=tmp] 
	tmp						<- merge(df, subset(tmp,select=c(Patient,df.idx)), by="Patient")		
	tmp						<- tmp[, 	{
				x<- matrix( as.numeric(c(StartTime, StopTime)), ncol=2 )
				x<- x[order(x[,2]),]
				list(NAStartTime.canberemoved=which(is.na(x[,1]))==1, df.idx=df.idx[1])					
			}, by=Patient]
	tmp						<- subset(tmp, NAStartTime.canberemoved)
	if(verbose)	cat(paste("\nnumber of entries with is.na(StartTime) & TrI=='Yes' before start of treatment , n=",nrow(tmp),"REMOVE"))
	set(df, tmp[,df.idx], "StopTime", NA)		
	df						<- subset(df, !is.na(StopTime))
	if(verbose)	cat(paste("\nnumber of entries with !is.na(StartTime) & !is.na(StopTime), n=",nrow(df)))
	#
	#	fix inconsistent timings
	#
	tmp		<- which(df[,StartTime>StopTime])		
	tmp		<- cbind( tmp, sapply(tmp,function(x)		which(df[, Patient==df[x,Patient] & StopTime==df[x,StartTime]])	) )		#second col contains StopTime
	for(i in seq_len(nrow(tmp)))
	{
		if(verbose)	cat(paste("\nprocess",df[tmp[i,1],Patient]))
		if(df[tmp[i,1],StartTime_Acc=="Acc"])			#set StartTime_Acc=="NAccD"
		{
			set(df,tmp[i,1],"StartTime_Acc","NAccD")			
			set(df,tmp[i,2],"StopTime_Acc","NAccD")			
		}
		if(df[tmp[i,1],StartTime_Acc=="NAccD"])			#reset StopTime		StartTime_Acc=="NAccD"
		{
			z		<- as.POSIXlt(df[tmp[i,2],StopTime])
			z$mday	<- 1
			z		<- as.Date(z)
			z		<- max( df[tmp[i,2],StartTime],z )
			set(df,tmp[i,1],"StartTime",z)
			set(df,tmp[i,2],"StopTime",z)
		}
		if(df[tmp[i,1],StartTime]>df[tmp[i,1],StopTime])
		{
			set(df,tmp[i,1],"StartTime_Acc","NAccMD")			
			set(df,tmp[i,2],"StopTime_Acc","NAccMD")
		}
		if(df[tmp[i,1],StartTime_Acc=="NAccMD"])		#set to midpoint
		{
			z<- df[tmp[i,2],StartTime] + difftime(df[tmp[i,1],StopTime],df[tmp[i,2],StartTime],units="days") / 2
			set(df,tmp[i,1],"StartTime",z)
			set(df,tmp[i,2],"StopTime",z)
		}		
	}
	#
	#	get 1st sort by Patient, 2nd sort by StopTime -> treatment history by Patient is now in order
	#
	setkey(df, StopTime)
	setkey(df, Patient)	
	#subset( df,StartTime_Acc=="NAccMD" & StopTime_Acc!="NAccMD")
	#subset( df, Patiet=="M10027")		
	#
	#	simple statistics of Patient history: number of treatment interruptions, total length of interruption in months and proportion of treatment interruption in treatment history 		
	#
	tmp		<- df[, 	{ 
				x		<- data.table( StartTime, StopTime, TrI, StartTime_Acc )
				x		<- subset(x, !is.na(StartTime))		#discard first time period if unknown StartTime
				z		<- subset(x, TrI=="Yes")
				list( 	TrI.n		= nrow(z), 
						TrI.mo		= as.numeric(sum(z[,difftime(StopTime,StartTime,units="days")/30])), 
						TrI.p		= as.numeric(sum(z[,difftime(StopTime,StartTime,units="days")/30])) / as.numeric(sum(x[,difftime(StopTime,StartTime,units="days")/30])), 
						HAART_T1	= HAART_T1[1], 
						AnyT_T1		= subset(x, TrI=="No")[,StartTime][1], 
						AnyT_T1_Acc	= subset(x, TrI=="No")[,StartTime_Acc][1]		)																					
			}, by=Patient]
	if( nrow(subset(tmp,!is.na(HAART_T1) & HAART_T1<AnyT_T1)) )	stop("found AnyT_T1 that is older than HAART_T1")		
	df		<- merge(subset(tmp,select=c(Patient,AnyT_T1,AnyT_T1_Acc,TrI.n, TrI.mo, TrI.p)), df, all.y=1,by="Patient")
	#
	#	reset AnyT_T1 conservatively
	#
	nacc				<- which(df[,AnyT_T1_Acc=="NAccD"])		
	tmp					<- as.POSIXlt(df[nacc,AnyT_T1] )
	tmp$mday			<- 30
	set(df, nacc, "AnyT_T1", as.Date(tmp))
	nacc				<- which(df[,AnyT_T1_Acc=="NAccMD"])
	tmp					<- as.POSIXlt(df[nacc,AnyT_T1] )
	tmp$mday			<- 31
	tmp$mon				<- 11
	set(df, nacc, "AnyT_T1", as.Date(tmp))		
	#
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)
}
######################################################################################
project.hivc.Excel2dataframe.CD4<- function(dir.name= DATA, verbose=1)
{
	NA.time			<- c("","01/01/1911","11/11/1911","24/06/1923")		
	verbose			<- 1
	#read CD4 csv data file and preprocess
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Immu.csv",sep='/')
	df				<- read.csv(file, stringsAsFactors=FALSE)											
	date.var		<- c("DateImm")		
	for(x in date.var)
	{
		cat(paste("\nprocess Time", x))
		nok.idx			<- which( df[,x]==NA.time[1] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[1],", n=", length(nok.idx)))
		#if(verbose && length(nok.idx))	cat(paste("\nentries with format ",NA.time[1],", Patient", paste(df[nok.idx,"Patient"],collapse=', ')))			
		if(length(nok.idx))	
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[2] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[2],", n=", length(nok.idx)))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[3] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[3],", n=", length(nok.idx)))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[4] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[4],", n=", length(nok.idx)))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		df[,x]			<- as.Date(df[,x], format="%d/%m/%Y")	
	}		
	df		<- data.table(df, key="Patient")
	setnames(df, "DateImm","PosCD4")
	set(df, NULL, "Patient", factor(df[,Patient]))
	if(verbose) cat(paste("\nnumber of entries read, n=",nrow(df)))
	tmp		<- which(df[,is.na(CD4A)])
	if(verbose) cat(paste("\nnumber of entries with is.na(CD4A), n=",length(tmp),"SETTING PosCD4 to NA"))
	set(df, tmp, "PosCD4", NA)
	df		<- subset(df,!is.na(PosCD4))
	if(verbose) cat(paste("\nnumber of entries with !is.na(PosCD4), n=",nrow(df)))
	#
	#	data corrections from Ard
	#
	tmp		<- which(df[, Patient=="M11392" & PosCD4=="2000-10-30" & CD4A==3001])
	if(verbose) cat(paste("\nnumber of entries with corrected CD4 units Patient=='M11392' & CD4A>3001, n=",length(tmp),"double entry, remove"))
	set(df, tmp, "CD4A", NA)	
	tmp		<- which(df[, Patient=="M26334" & PosCD4=="2008-09-11" & CD4A==2230])
	if(verbose) cat(paste("\nnumber of entries with corrected CD4 units Patient==M26334 & PosCD4==2008-09-11 & CD4A==2230, n=",length(tmp),"double entry, remove"))
	set(df, tmp, "CD4A", NA)		
	tmp		<- which(df[, 	Patient=='M26293' & PosCD4=='2003-09-12' & CD4A==800])
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M26293' & PosCD4=='2003-09-12' & CD4A==800, n=",length(tmp),"set to 16/9/2003"))
	set(df, tmp, "PosCD4", as.Date('2003-09-16'))
	tmp		<- which(df[, 	Patient=='M11233' & PosCD4=='2001-01-10' & CD4A==24])
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M11233' & PosCD4=='2001-01-10' & CD4A==24, n=",length(tmp),"set to 25/4/1996"))
	set(df, tmp, "PosCD4", as.Date('1996-04-25'))
	tmp		<- which(df[, 	Patient=='M13124' & PosCD4=='1998-08-07'])
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M13124' & PosCD4=='1998-08-07', n=",length(tmp),"set to 7/9/1998"))
	set(df, tmp, "PosCD4", as.Date('1998-09-07'))	
	tmp		<- which(df[, 	Patient=='M13124' & PosCD4=='2001-09-14'])
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M13124' & PosCD4=='2001-09-14', n=",length(tmp),"set to 14/9/2000"))
	set(df, tmp, "PosCD4", as.Date('2000-09-14'))
	tmp		<- which(df[, 	Patient=="M10544" & PosCD4=="2003-02-17" & CD4A==850	|
							Patient=="M11099" & PosCD4=="1997-12-30" & CD4A==1240	|
							Patient=="M11133" & PosCD4=="2003-06-16" & CD4A==170	|
							Patient=="M11137" & PosCD4=="2003-06-25" & CD4A==460	|
							Patient=="M11167" & PosCD4=="2006-09-04" & CD4A==400	|	
							Patient=="M11351" & PosCD4=="1996-10-08" & CD4A==150	|
							Patient=="M11713" & PosCD4=="1996-07-03" & CD4A==0.37	|
							Patient=="M12577" & PosCD4=="2000-09-25" & CD4A==210	|
							Patient=="M12884" & PosCD4=="1997-11-04" & CD4A==350	|
							Patient=="M13051" & PosCD4=="1998-06-08" & CD4A==460	|
							Patient=="M13124" & PosCD4=="2001-10-09" & CD4A==1.17	|
							Patient=="M13124" & PosCD4=="2003-02-12" & CD4A==0.74	|
							Patient=="M13124" & PosCD4=="2003-03-21" & CD4A==0.59	|
							Patient=="M13124" & PosCD4=="2003-06-17" & CD4A==0.61	|
							Patient=="M13126" & PosCD4=="2001-01-08" & CD4A==0.5	|
							Patient=="M13126" & PosCD4=="2001-03-05" & CD4A==0.43	|
							Patient=="M13126" & PosCD4=="2003-01-17" & CD4A==0.48	|
							Patient=="M13126" & PosCD4=="2003-04-28" & CD4A==0.48	|
							Patient=="M13126" & PosCD4=="2003-07-24" & CD4A==0.46	|	#no CD4 anymore at 24/7/2003, delete both
							Patient=="M13126" & PosCD4=="2003-07-24" & CD4A==460	|
							Patient=="M13298" & PosCD4=="1997-12-11" & CD4A==760	|
							Patient=="M14834" & PosCD4=="1998-12-10" & CD4A==990	|
							Patient=="M14945" & PosCD4=="2000-09-20" & CD4A==1620	|
							Patient=="M14986" & PosCD4=="2001-03-29" & CD4A==640	|
							Patient=="M14995" & PosCD4=="1999-01-12" & CD4A==370	|
							Patient=="M15071" & PosCD4=="1999-10-12" & CD4A==100	|
							Patient=="M15234" & PosCD4=="1998-09-02" & CD4A==25		|
							Patient=="M15234" & PosCD4=="1998-11-04" & CD4A==32		|
							Patient=="M15234" & PosCD4=="1998-12-16" & CD4A==35		|
							Patient=="M16018" & PosCD4=="1998-09-09" & CD4A==1400	|
							Patient=="M16570" & PosCD4=="2003-06-17" & CD4A==280	|
							Patient=="M16622" & PosCD4=="2000-09-27" & CD4A==100	|
							Patient=="M17154" & PosCD4=="2011-01-05" & CD4A==495	|
							Patient=="M17819" & PosCD4=="2000-03-15" & CD4A==1010	|
							Patient=="M17951" & PosCD4=="1999-07-22" & CD4A==530	|
							Patient=="M18712" & PosCD4=="2000-07-31" & CD4A==600	|
							Patient=="M25530" & PosCD4=="2002-04-09" & CD4A==59		|
							Patient=="M25530" & PosCD4=="2002-04-19" & CD4A==66		|	
							Patient=="M28189" & PosCD4=="2007-09-04" & CD4A==31		|
							Patient=="M30605" & PosCD4=="2011-08-22" & CD4A==83		|
							Patient=="M31573" & PosCD4=="2011-03-24" & CD4A==0		|	#import from hospital db, remove unlikely entry
							Patient=="M33353" & PosCD4=="2006-03-22" & CD4A==101	|
							Patient=="M33924" & PosCD4=="2007-11-01" & CD4A==0		|
							Patient=="M37294" & PosCD4=="2011-07-18" & CD4A==820	|
							Patient=="M39055" & PosCD4=="2012-02-06" & CD4A==6850							
							])
	if(verbose) cat(paste("\nnumber of entries with incorrect CD4  should be 45, n=",length(tmp),"set to NA"))	
	set(df, tmp, "CD4A", NA)	
	#	adjust likely wrong units > 20000
	tmp		<- which(df[, CD4A>20000])
	if(verbose) cat(paste("\nnumber of entries with likely wrong CD4 units > 20000, n=",length(tmp),"DIVIDE BY 1e3"))
	if(verbose) print(df[tmp,])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-3,d=0))
	#	adjust likely wrong units > 10000
	if(verbose) cat(paste("\npatient M11368"))
	tmp		<- which(df[, CD4A>10000 & Patient=="M11368"])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-3,d=0))
	if(verbose) cat(paste("\npatient M32958"))
	tmp		<- which(df[, Patient=="M32958" & PosCD4<="2010-05-03"])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	if(verbose) cat(paste("\npatient M20633"))
	tmp		<- which(df[, Patient=="M20633" & CD4A>5000])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-3,d=0))
	tmp		<- which(df[, CD4A>3000])
	if(verbose) cat(paste("\nnumber of entries with likely wrong CD4 units > 3000, n=",length(tmp),"DIVIDE BY 10"))
	if(verbose) print(df[tmp,])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	#
	#	check above 1700 manually
	#
	tmp<- merge(df, subset(df, CD4A>1700, Patient), by="Patient")
	tmp<- tmp[,	list(CD4.med= median(CD4A), CD4.max=max(CD4A)),by="Patient"]
	print( subset(tmp, CD4.med*1.5<CD4.max) )
	#	divide by 10
	tmp		<- which(df[, Patient%in%c("M10212","M14927","M15431","M15519","M20720","M26334","M27643","M27571") & CD4A>1700])
	if(verbose) cat(paste("\nnumber of entries with likely wrong CD4 units > 1700  M10212 M27571 M14927  M15431  M15519  M20720  M26334  M27643, n=",length(tmp),"DIVIDE BY 10"))
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	if(verbose) cat(paste("\npatient M17554"))
	tmp		<- which(df[, Patient=="M17554" & CD4A%in%c(2760,2170)])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	if(verbose) cat(paste("\npatient M17554"))
	tmp		<- which(df[, Patient=="1500" & CD4A%in%c(1500)])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	#	set to NA
	tmp		<- which(df[, Patient%in%c("M12953","M13340","M26537","M26669","M35668") & CD4A>1900])
	if(verbose) cat(paste("\nnumber of entries with likely wrong CD4 units > 1900  M12953   M13340  M26537  M26669  M35668, n=",length(tmp),"SET to NA"))
	set(df, tmp, "CD4A", NA)
	tmp		<- which(df[, Patient%in%c("M15743") & CD4A>2500])
	if(verbose) cat(paste("\nnumber of entries with likely wrong CD4 units > 2500  M15743, n=",length(tmp),"SET to NA"))
	set(df, tmp, "CD4A", NA)
	tmp		<- which(df[, Patient%in%c("M30155") & CD4A>1000])
	if(verbose) cat(paste("\nnumber of entries with likely wrong CD4 units > 1000  M30155, n=",length(tmp),"SET to NA"))
	set(df, tmp, "CD4A", NA)
	#
	df		<- subset(df, !is.na(CD4A))
	if(verbose) cat(paste("\nnumber of entries with !is.na(PosCD4), n=",nrow(df)))
	#
	#	check for too small entries
	#
	tmp		<- which(df[,CD4A<1 & CD4A>0])
	if(verbose) cat(paste("\nnumber of entries with too small CD4 CHECK MANUALLY, n=",length(tmp),"SET TO *1e3"))
	if(verbose)	print(df[tmp,])
	set(df, tmp, "CD4A", df[tmp,CD4A]*1e3)
	#
	#	check for double entries
	#
	tmp		<- df[,	{
						z	<- which(as.numeric(difftime(PosCD4[-1],PosCD4[-length(PosCD4)],units="days"))==0)
						z	<- rbind(z,z+1)
						z2	<- rep( apply( z, 2, function(z2)	abs(CD4A[z2[1]]-CD4A[z2[2]])	), each=2 )
						z3	<- rep( apply( z, 2, function(z2)	mean(c(CD4A[z2[1]],CD4A[z2[2]]))	), each=2 )
						z	<- as.numeric(z)
						list( PosCD4=PosCD4[z], CD4A=CD4A[z], CD4d=z2, CD4mean=z3 ) 					
					}, by="Patient"]
	if(verbose) cat(paste("\nfound double entries, n=",nrow(tmp),"SET TO MEAN -- NOT ALWAYS OK"))
	if(verbose)	print(tmp)	
	for( i in seq.int(1,nrow(tmp),2))
	{
		z<- which(df[, Patient==tmp[i,Patient] & PosCD4==tmp[i,PosCD4]])
		set(df,z[1],"CD4A", mean(df[z,CD4A]))
		set(df,z[-1],"CD4A", NA)		
	}
	df		<- subset(df, !is.na(CD4A))
	if(verbose) cat(paste("\nnumber of entries with !is.na(PosCD4), n=",nrow(df)))
	#	
	setkey(df, PosCD4)
	setkey(df, Patient)
	#
	df		<- df[, 	{
				z<- which.min(PosCD4)
				list(PosCD4=PosCD4, CD4=CD4A, PosCD4_T1=PosCD4[z], CD4_T1=CD4A[z] ) 	
			},by=Patient]
	
	
	
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)		
}
######################################################################################
project.hivc.Excel2dataframe.Viro<- function()		
{
	verbose				<- 1
	dir.name			<- DATA
	DB.locktime			<- HIVC.db.locktime
	
	#need for checking of VL data
	file.treatment		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	load(file.treatment)
	df.treat			<- subset(df, select=c(Patient, StartTime, StopTime, AnyT_T1, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel))
	
	NA.time				<- c("","01/01/1911","11/11/1911","24/06/1923")		
	RNA.min				<- 400	#seems to be standard value
	RNA.max				<- 5e6
	lRNA.min.infectious	<- log10(1e4)
	lRNA.min.early		<- log10(1e5)
	lRNA.max.b4early	<- log10(2e4)
	
	#read VIROLOGY csv data file and preprocess
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Viro.csv",sep='/')
	df				<- read.csv(file, stringsAsFactors=FALSE)									
	df$Undetectable	<- factor(df$Undetectable, levels=c(0,1,2),labels=c("No","Yes","LargerThan"))
	date.var		<- c("DateRNA")		
	for(x in date.var)
	{
		cat(paste("\nprocess Time", x))
		nok.idx			<- which( df[,x]==NA.time[1] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[1],", n=", length(nok.idx)))
		#if(verbose && length(nok.idx))	cat(paste("\nentries with format ",NA.time[1],", Patient", paste(df[nok.idx,"Patient"],collapse=', ')))			
		if(length(nok.idx))	
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[2] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[2],", n=", length(nok.idx)))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[3] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[3],", n=", length(nok.idx)))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		nok.idx			<- which( df[,x]==NA.time[4] )
		if(verbose)	cat(paste("\nentries with format ",NA.time[4],", n=", length(nok.idx)))
		if(length(nok.idx))
			df[nok.idx,x]	<- NA
		df[,x]			<- as.Date(df[,x], format="%d/%m/%Y")	
	}		
	
	df		<- data.table(df, key="Patient")
	if(verbose) cat(paste("\nnumber of entries read, n=",nrow(df)))
	setnames(df, "DateRNA","PosRNA")
	set(df, NULL, "Patient", factor(df[,Patient]))
	#
	#	checking manually NegT>PosRNA
	#
	if(verbose)		cat(paste("\nset entry Patient=='M38400' & as.character(PosRNA)=='2005-11-24' to NA -- seems unlikely"))
	if(verbose)		print( subset(df, Patient=="M38400") )
	tmp		<- which( df[, Patient=="M38400" & as.character(PosRNA)=="2005-11-24"] )
	if(verbose)		cat(paste("\nsetting number of entries to NA, n=",length(tmp)))
	set(df,tmp,"PosRNA",NA)		
	if(verbose)		cat(paste("\nset entry Patient=='M36146' & as.character(PosRNA)=='2005-11-19' to NA -- seems unlikely"))
	if(verbose)		print( subset(df, Patient=="M36146") )
	tmp		<- which( df[, Patient=="M36146" & as.character(PosRNA)=="2005-11-19"] )
	if(verbose)		cat(paste("\nsetting number of entries to NA, n=",length(tmp)))
	set(df,tmp,"PosRNA",NA)
	# remove is.na(PosRNA) and !is.na(RNA)
	df		<- subset(df, !is.na(PosRNA) & !is.na(RNA))
	if(verbose)		cat(paste("\nnumber of entries with !is.na(PosRNA) & !is.na(RNA), n=",nrow(df)))		
	#
	#	#checking for duplicates -- do not matter that much - leave if any
	#
	#subset(df[,	{
	#				z<- diff(RNA); tmp<- (RNA>1e4 & !is.na(Undetectable))[-1]; length(which(z[tmp]==0))
	#			}				
	#			,by="Patient"], V1>0 )
	#
	#	combine Undetectable=="LargerThan" with Undetectable=="No"
	#
	tmp		<- which(df[,Undetectable=="LargerThan"])
	if(verbose)		cat(paste("\nsetting Undetectable=='LargerThan' to Undetectable=='No', n=",length(tmp)))
	set(df,tmp,"Undetectable","No")
	if(any(df[,is.na(Undetectable)]))	stop("unexpected is.na(Undetectable)")		
	set(df, NULL, "Undetectable",factor(as.numeric(df[,Undetectable]), levels=c(1,2),labels=c("No","Yes")))		
	#
	#	set Undetectable=="Yes" and RNA<RNA.min to Undetectable=="No" and RNA.min
	#
	tmp<- which( df[, Undetectable=="Yes" & RNA<1e4] )
	if(verbose)		cat(paste("\nsetting Undetectable=='Yes' and RNA<RNA.min to Undetectable=='No' and RNA.min, n=",length(tmp)))
	set(df,tmp,"Undetectable","No")
	set(df,tmp,"RNA",RNA.min)
	#
	#	set RNA<RNA.min to RNA.min
	#		
	tmp<- which( df[, RNA<RNA.min] )
	if(verbose)		cat(paste("\nsetting RNA<RNA.min to RNA.min, n=",length(tmp)))		
	set(df,tmp,"RNA",RNA.min)
	#
	#	wrong units ? adjust manually
	#				
	tmp		<- which(df[, Patient%in%c("M11314","M11331","M40782","M14759") & RNA>5e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 5e6  -- M11314  M11331  M40782  M14759, n=",length(tmp),"SET to 5e6"))
	set(df, tmp, "RNA", 5e6)
	#
	tmp		<- which(df[, Patient=="M27377" & PosRNA>"2007-08-11"])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units -- patient M27377 after 2007-08-11, n=",length(tmp),"DIV by 1e3"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-3)
	#
	tmp		<- which(df[, Patient%in%c("M13134","M18385","M18767","M20308","M35814","M35852","M36515","M41877") & RNA>1e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M13134  M18385  M18767  M20308  M35814  M35852  M36515  M41877, n=",length(tmp),"DIV by 10"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M38031") & RNA>1e5])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e5  -- M38031, n=",length(tmp),"DIV by 10"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M38036","M33131","M33668","M34200","M34302","M20350") & RNA>1e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M38036  M33131  M33668  M34200  M34302 M20350, n=",length(tmp),"DIV by 100"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-2)
	#
	tmp		<- which(df[, Patient%in%c("M11995","M13266","M14486","M15621","M16588") & RNA>5e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 5e6  -- M11995 M13266 M14486 M15621 M16588, n=",length(tmp),"DIV by 10"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)		
	#
	tmp		<- which(df[, Patient%in%c("M16570") & RNA>2e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 2e6  -- M16570, n=",length(tmp),"DIV by 10"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M17655","M37746","M30733","M27377") & RNA>2e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 2e6  -- M17655 M37746 M30733 M27377, n=",length(tmp),"NA"))
	set(df, tmp, "RNA", NA)
	#
	tmp		<- which(df[, Patient%in%c("M17554") & RNA>1e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M17554, n=",length(tmp),"DIV by 100"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-2)
	#
	tmp		<- which(df[, Patient%in%c("M10607","M10969","M11428","M28707","M31388","M31455","M32401","M32877","M33406","M33839","M33918","M34062","M35280","M30788") & RNA>=1e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M10607 M10969 M11428 M28707 M31388 M31455 M32401 M32877 M33406 M33839 M33918 M34062 M35280  M30788, n=",length(tmp),"DIV by 10"))
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	# double check entries manually in range >5e6
	#		
	tmp		<- merge(subset(df, RNA>5e6, c(Patient,PosRNA)), df.treat, by="Patient")
	tmp		<- tmp[, 	{
				z<- which( difftime(StopTime,PosRNA,units="days")>0 )
				z<- z[1]
				list(PosRNA=PosRNA[z], StartTime=StartTime[z], StopTime=StopTime[z], TrI=TrI[z], TrCh.failure=TrCh.failure[z], TrCh.adherence=TrCh.adherence[z], TrCh.patrel=TrCh.patrel[z])				
			}, by="Patient"]
	setnames(tmp, "PosRNA","PosRNAh")
	tmp		<- merge(tmp, df,by="Patient")		
	tmp[, print(data.table(Patient,RNA,PosRNA,PosRNAh,StartTime,StopTime,TrI,TrCh.failure,TrCh.adherence,TrCh.patrel)), by="Patient"]
	#
	# double check entries manually in range >RNA<5e6 & RNA>2e6
	#		
	tmp		<- merge(subset(df, RNA<5e6 & RNA>2e6, c(Patient,PosRNA)), df.treat, by="Patient")
	tmp		<- tmp[, 	{
				z<- which( difftime(StopTime,PosRNA,units="days")>0 )
				z<- z[1]
				list(PosRNA=PosRNA[z], StartTime=StartTime[z], StopTime=StopTime[z], TrI=TrI[z], TrCh.failure=TrCh.failure[z], TrCh.adherence=TrCh.adherence[z], TrCh.patrel=TrCh.patrel[z])				
			}, by="Patient"]
	setnames(tmp, "PosRNA","PosRNAh")
	tmp		<- merge(tmp, df,by="Patient")		
	tmp[, print(data.table(Patient,RNA,PosRNA,PosRNAh,StartTime,StopTime,TrI,TrCh.failure,TrCh.adherence,TrCh.patrel)), by="Patient"]
	#
	tmp		<- which(df[, Patient=="M12612" & PosRNA=="1996-06-13"])
	if(verbose) cat(paste("\nset  M12612 1996-06-13 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	tmp		<- which(df[, Patient=="M17044" & PosRNA=="2004-07-19"])
	if(verbose) cat(paste("\nset  M17044 2004-07-19 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	#
	# check entries manually with Undetectable=="Yes" & RNA>1e4
	#		
	tmp		<- merge(subset(df, Undetectable=="Yes" & RNA>1e4, c(Patient,PosRNA)), df.treat, by="Patient")
	tmp		<- tmp[, 	{
				z<- which( difftime(StopTime,PosRNA,units="days")>0 )
				z<- z[1]
				list(PosRNA=PosRNA[z], StartTime=StartTime[z], StopTime=StopTime[z], TrI=TrI[z], TrCh.failure=TrCh.failure[z], TrCh.adherence=TrCh.adherence[z], TrCh.patrel=TrCh.patrel[z])				
			}, by="Patient"]
	setnames(tmp, "PosRNA","PosRNAh")
	tmp		<- merge(tmp, df,by="Patient")		
	tmp[, print(data.table(Patient,RNA,PosRNA,Undetectable,StartTime,StopTime,TrI,TrCh.failure,TrCh.adherence,TrCh.patrel)), by="Patient"]
	#
	tmp		<- which(df[, Patient=="M30584" & PosRNA=="2004-10-14"])
	if(verbose) cat(paste("\nset  M30584 2004-10-14 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	tmp		<- which(df[, Patient=="M26369" & PosRNA=="2003-10-13"])[1]
	if(verbose) cat(paste("\nset  M26369 2003-10-13 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	tmp		<- which(df[, Patient=="M16566" & PosRNA=="2000-11-16"])
	if(verbose) cat(paste("\nset  M16566 2000-11-16 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	tmp		<- which(df[, Patient=="M12818" & PosRNA=="1996-09-18" & RNA>=1e4])
	if(verbose) cat(paste("\nset  M12818 1996-09-18 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)		
	#
	#	set Undetectable=='Yes' & RNA>1e4
	#
	tmp		<- which(df[, Undetectable=="Yes" & RNA>=1e4])
	if(verbose) cat(paste("\nsetting Undetectable=='Yes' & RNA>1e4 to Undetectable=='No', n=",length(tmp)))
	set(df, tmp, "Undetectable", "No")		
	#
	#	set RNA>RNA.max to RNA.max
	#		
	tmp		<- which(df[, RNA>RNA.max])
	if(verbose) cat(paste("\nsetting RNA>RNA.max to RNA.max, n=",length(tmp)))
	set(df, tmp, "RNA", RNA.max)		
	#
	df		<- subset(df,!is.na(RNA), select=c(Patient, PosRNA, RNA))
	if(verbose) cat(paste("\nentries with !is.na(RNA), n=",nrow(df)))
	#
	#	average duplicate entries out
	#						
	if(verbose) cat(paste("\nremoving duplicate entries"))
	df		<- df[,	{
				z		<- as.numeric(difftime(PosRNA[-1], PosRNA[-length(PosRNA)],units="days"))					
				if(length(which(z==0)))
				{
					RNA.d				<- sapply(which(z==0), function(i)  mean( RNA[seq.int(i,i+1)] ) )					
					PosRNA.d			<- PosRNA[z!=0]		#keep the last of the duplicates
					z2					<- RNA				
					z2[which(z==0)+1]	<- RNA.d			#overwrite the last of the duplicates
					RNA.d				<- z2[z!=0]			#keep the last of the duplicates
					#print(z); print(which(z==0)); print(RNA.d); print(list( PosRNA= PosRNA.d, RNA= RNA.d))
				}
				else
					ans<- list( PosRNA= PosRNA, RNA= RNA)
				ans					
			},by="Patient"]
	if(verbose) cat(paste("\nnumber of entries, n=",nrow(df)))
	setkey(df, PosRNA)
	setkey(df, Patient)
	#
	df[,"lRNA":=round(log10( df[,RNA] ), d=3)]		
	#	
	#	compute several statistics on lRNA life history
	#	PoslRNA_T1		time of first lRNA
	#	lRNA_T1			first lRNA 
	#	lRNA.i 			proportion of time spent above 'lRNA.min.infectious'
	#	lRNA.hb4tr_LT 	last time of lRNA above 'lRNA.min.early'
	#	lRNA.early		is there increasing lRNA before treatment 
	#	
	df					<- df[, {		
									z<- which.min(PosRNA)
									list(PosRNA=PosRNA, RNA=RNA, lRNA=lRNA, PoslRNA_T1=PosRNA[z], lRNA_T1=lRNA[z])
								},by=Patient]
	tmp					<- subset(df.treat,select=c(Patient,AnyT_T1))
	setkey(tmp,Patient)		 
	df					<- merge(df, unique(tmp), all.x=1, by="Patient")
	set(df,which(df[,is.na(AnyT_T1)]),"AnyT_T1",DB.locktime)
	df[, lRNA.infectious:=lRNA>=lRNA.min.infectious]
	df[, lRNA.high.b4tr	:=lRNA>=lRNA.min.early & PosRNA<AnyT_T1]
	
	tmp	<- df[,		{
						z				<- data.table(PosRNA, StopRNA=c(PosRNA[-1],DB.locktime), lRNA, AnyT_T1, lRNA.infectious, lRNA.high.b4tr)
						p.infectious	<- z[,as.numeric(difftime(StopRNA, PosRNA,units="days")/30)]							#difftime between subsequent PosRNA
						p.infectious	<- sum(p.infectious[ which(z[,lRNA.infectious]) ]) / sum(p.infectious)	#prop of time in lRNA.infectious
						lt.highb4tr		<- subset(z,lRNA.high.b4tr)
						if(nrow(lt.highb4tr))
							lt.highb4tr	<- lt.highb4tr[,StopRNA][nrow(lt.highb4tr)]	#last time before tr that VL was high
						else
							lt.highb4tr	<- as.Date(NA)
						if(!is.na(lt.highb4tr) && lt.highb4tr>AnyT_T1[1])
							lt.highb4tr	<- AnyT_T1[1]
						early			<- ifelse(	any(lRNA.high.b4tr)  && any( lRNA[seq_len( which(lRNA.high.b4tr)[1] )]<lRNA.max.b4early )	,TRUE,FALSE)
						list(lRNA.i= p.infectious, lRNA.hb4tr_LT=lt.highb4tr, lRNA.early= early )										
					}, by="Patient"]
	df<- merge(subset(df,select=c(Patient, PosRNA, RNA, lRNA, PoslRNA_T1, lRNA_T1)), subset(tmp,select=c(Patient,lRNA.i,lRNA.hb4tr_LT,lRNA.early)), all.x=1, by="Patient")
	
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)
}	
######################################################################################
project.hivc.Excel2dataframe<- function(dir.name= DATA, min.seq.len=21, verbose=1)
{
	if(0)
	{
		#read SEQUENCE csv data file and preprocess				
		names.GeneCode	<- c("PROT","RT")
		NA.DateRes		<- as.Date("1911-11-11") 
		
		file			<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.csv",sep='/')
		df				<- read.csv(file, stringsAsFactors=FALSE)				
		proc.GeneCode	<- c(1,2)
		df				<- lapply(proc.GeneCode,function(gene)
				{
					cat(paste("\nprocess GeneCode", gene))
					tmp				<- df[ df[,"GeneCode"]==gene, c("Patient","SampleCode","DateRes","Sequence"), drop=0 ]
					tmp[,"DateRes"]	<- as.Date(tmp[,"DateRes"], format="%d/%m/%Y")
					
					nok.idx<- which( tmp[,"DateRes"]==NA.DateRes )					
					if(verbose) cat(paste("\nentries with missing DateRes, n=", length(nok.idx)))					
					if(verbose) cat(paste("\nentries with missing DateRes, SampleCode", paste(tmp[nok.idx,"SampleCode"],collapse=', ')))
					tmp[nok.idx,"DateRes"]<- NA
					if(verbose) cat(paste("\nrange of DateRes is",paste(range(tmp[,"DateRes"], na.rm=1),collapse=', ')))
					cat(paste("\nfound n=", nrow(tmp)))
					seq.len			<- nchar(tmp[,"Sequence"])
					nok.idx		<- which(seq.len<min.seq.len)
					seq.ok.idx		<- which(seq.len>=min.seq.len)
					if(verbose)	cat(paste("\ndiscarding sequences with insufficient length, n=",length(nok.idx),'\n'))
					if(verbose)	cat(paste("\ndiscarding sequence with insufficient length, SampleCode",paste(tmp[nok.idx,"SampleCode"],collapse=', ')))
					tmp				<- tmp[seq.ok.idx,]
					if(verbose) cat(paste("\nfinal n=", nrow(tmp)))
					tmp
				})		
		names(df)	<- names.GeneCode
		file		<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.R",sep='/')
		cat(paste("\nsave to", file))
		save(df, file=file)
		quit("no")
	}
	if(0)
	{
		project.hivc.Excel2dataframe.Regimen()						
	}
	if(0)
	{
		NA.Acute			<- c(NA,9)
		NA.CountryInfection	<- c(NA,"")
		NA.time				<- c("","01/01/1911","11/11/1911")		
		NA.transmission		<- 900
		#read PATIENT csv data file and preprocess
		file				<- paste(dir.name,"derived/ATHENA_2013_03_Patient.csv",sep='/')
		df					<- read.csv(file, stringsAsFactors=FALSE)									
		df$isDead			<- as.numeric( df[,"DateDied"]!="")
		df[which(df[,"Transmission"]==NA.transmission),"Transmission"]<- NA

		date.var			<- c("DateBorn","MyDateNeg1","MyDatePos1","DateDied","DateLastContact","DateFirstEverCDCC")		
		for(x in date.var)
		{
			cat(paste("\nprocess Time", x))
			nok.idx			<- which( df[,x]==NA.time[1] )
			if(verbose)	cat(paste("\nentries with format ",NA.time[1],", n=", length(nok.idx)))
			#if(verbose && length(nok.idx))	cat(paste("\nentries with format ",NA.time[1],", Patient", paste(df[nok.idx,"Patient"],collapse=', ')))
			
			if(length(nok.idx))	
				df[nok.idx,x]	<- NA
			nok.idx			<- which( df[,x]==NA.time[2] )
			if(verbose)	cat(paste("\nentries with format ",NA.time[2],", n=", length(nok.idx)))
			if(length(nok.idx))
				df[nok.idx,x]	<- NA
			nok.idx			<- which( df[,x]==NA.time[3] )
			if(verbose)	cat(paste("\nentries with format ",NA.time[3],", n=", length(nok.idx)))
			if(length(nok.idx))
				df[nok.idx,x]	<- NA
			df[,x]			<- as.Date(df[,x], format="%d/%m/%Y")	
		}
		file		<- paste(dir.name,"derived/ATHENA_2013_03_Patients.R",sep='/')
		if(verbose) cat(paste("\nsave to", file))
		
		df<- data.table(df)
		setnames(df, c("MyDateNeg1_Acc","MyDatePos1_Acc","AcuteInfection","Transmission","HospitalRegion"), c("NegT_Acc","PosT_Acc","isAcute","Trm","RegionHospital"))
		set(df, NULL, "Patient", factor(df[,Patient]))
		set(df, which( df[,isAcute%in%NA.Acute] ), "isAcute", NA )		
		set(df, NULL, "isAcute", factor(df[,isAcute], levels=c(0,1,2), labels=c("No","Yes","Maybe")) )
		set(df, NULL, "Sex", factor(df[,Sex], levels=c(1,2), labels=c("M","F")))
		set(df, NULL, "Subtype", factor(df[,Subtype]))
		set(df, NULL, "CountryBorn", factor(df[,CountryBorn]))
		set(df, which( df[,CountryInfection%in%NA.CountryInfection] ), "CountryInfection", NA_character_ )		
		set(df, NULL, "CountryInfection", factor(df[,CountryInfection]))
		set(df, NULL, "RegionOrigin", factor(df[,RegionOrigin]))
		set(df, NULL, "NegT_Acc", factor(df[,NegT_Acc], levels=c(0,1), labels=c("No","Yes")))
		set(df, NULL, "PosT_Acc", factor(df[,PosT_Acc], levels=c(0,1), labels=c("No","Yes")))		
		set(df, NULL, "isDead", factor(df[,isDead], levels=c(0,1), labels=c("No","Yes")))
		set(df, NULL, "Trm", factor(df[, Trm], levels=c(100, 101,  102,  202, 103,  104,  105,  106,  107, 108,  110), labels= c("MSM","BI","HET","HETfa","IDU","BLOOD","NEEACC", "PREG", "BREAST", "OTH", "SXCH")) )
		set(df, NULL, "RegionHospital", factor(df[,RegionHospital], levels=c(1,2,3,4,5,6), labels=c("Amst","N","E","S","W","Curu")))
		setkey(df,Patient)
		str(df)
		save(df, file=file)		
	}
	if(0)
	{
		project.hivc.Excel2dataframe.Viro()	
	}
	if(0)
	{
		project.hivc.Excel2dataframe.CD4()		
	}
}
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
	nrc.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=clu.pre)
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
	argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume)
	argv			<<- unlist(strsplit(argv,' '))
	nsh.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=nsh.clu.pre)
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
project.athena.Fisheretal.X.Trm.Region<- function(df.tpairs, clumsm.info)
{
	t.group	<- merge( data.table(Patient=df.tpairs[, unique(t.Patient)]), unique(subset( clumsm.info, select=c(Patient, RegionHospital, Trm))), by='Patient' )
	#	Exposure group - MSM or Other or NA
	set(t.group, t.group[, which(Trm!='MSM')], 'Trm', 'OTH' )
	set(t.group, NULL, 'Trm', as.factor(as.character(t.group[,Trm])) )
	set(t.group, NULL, 'RegionHospital', as.factor(as.character(t.group[,RegionHospital])) )
	tmp	<- table(t.group[,Trm])
	cat(paste('\nExposure group categories', paste(names(tmp), tmp, collapse=', ', sep='=')))
	tmp	<- table(t.group[,RegionHospital])
	cat(paste('\nRegion group categories', paste(names(tmp), tmp, collapse=', ', sep='=')))
	setnames(t.group, 'Patient','t.Patient')
	t.group		
}
######################################################################################
project.athena.Fisheretal.X.calendarperiod<- function(df.tpairs, clumsm.info, t.period= 0.25, c.nperiod= 4)
{
	i.diag		<- merge( unique(subset(df.tpairs, select=Patient)), unique(subset(clumsm.info, select=c(Patient, AnyPos_T1))), by='Patient' )				
	set(i.diag, NULL, 'AnyPos_T1', i.diag[, floor(AnyPos_T1) + floor( (AnyPos_T1%%1)*100 %/% (t.period*100) ) * t.period] )
	tmp			<- table(i.diag[,AnyPos_T1])
	i.cgroup	<- data.table(AnyPos_T1= as.numeric(names(tmp)), n= tmp, nsum= cumsum(tmp))	
	i.cgroup.n	<- ceiling( i.cgroup[ nrow(i.cgroup), nsum] / c.nperiod )
	cat(paste('\nTarget group size is',i.cgroup.n))
	i.cgroup[, t.period:= cut(i.cgroup[,nsum], breaks= seq(from= 0, by=i.cgroup.n, length.out=c.nperiod+1), labels= seq_len(c.nperiod))]
	tmp			<- i.cgroup[, list(n=sum(n), min=min(AnyPos_T1), max=max(AnyPos_T1)),by='t.period']
	cat(paste('\nGroup sizes are=',tmp[, paste(n, collapse=', ',sep='')]))
	cat(paste('\nGroup min diag times are=',tmp[, paste(min, collapse=', ',sep='')]))
	cat(paste('\nGroup max diag times are=',tmp[, paste(max, collapse=', ',sep='')]))	
	i.diag		<- merge(i.diag, subset(i.cgroup, select=c(AnyPos_T1, t.period)), by='AnyPos_T1')
	subset(i.diag, select=c(Patient, t.period)) 		
}
######################################################################################
project.athena.Fisheretal.X.followup<- function(df.tpairs, df.immu, t.period= 0.25, t.endctime=2013.0)
{
	follow		<- subset(df.immu, select=c(Patient, PosCD4))
	setkey(follow, Patient)
	#	get follow.up time periods for every potential transmitter
	#	follow.up timeline starts at first PosCD4 and ends at endctime
	follow		<- merge( data.table(Patient= df.tpairs[, unique(t.Patient)]),  follow, by='Patient' )
	set(follow, NULL, 'PosCD4', hivc.db.Date2numeric(follow[,PosCD4]))
	follow.t	<- follow[, list(PosCD4_T1= min(PosCD4)),by='Patient']
	set(follow.t, NULL, 'PosCD4_T1', follow.t[, floor(PosCD4_T1) + floor( (PosCD4_T1%%1)*100 %/% (t.period*100) ) * t.period] )
	follow.t	<- follow.t[, list(t= seq(PosCD4_T1, t.endctime, by=t.period)),by='Patient']
	#	get follow up quantiles that do not depend on the time periods t
	follow.q	<- follow[, {
				if(length(PosCD4)>1)
					tmp	<- diff(PosCD4)
				else
					tmp	<- NA_real_
				list(fw.up.mean= mean(tmp), fw.up.mx=max(tmp))
			},by='Patient']
	#plot(follow[, fw.up.mean], follow[, log(fw.up.mx)], pch=18)
	#hist(follow[, log(fw.up.mx)])
	follow.q		<- round( follow.q[, quantile(fw.up.mx, p=c(0, 0.75, 0.95, 1), na.rm=TRUE)], d=3)
	cat(paste('\nCD4 follow up quantiles are=',paste(follow.q, collapse=' ', sep=''), sep=''))
	#	determine fw.up.mx per patient up to time period t
	follow.t		<- merge(follow, follow.t, by='Patient', allow.cartesian=TRUE)
	cat(paste('\nmerge gives nrows=',nrow(follow.t)))
	follow.t		<- follow.t[,	list(fw.up.mx= ifelse(length(which(PosCD4<=t))>1, max(diff(PosCD4[PosCD4<=t])), NA_real_) ),		by=c('Patient','t')]
	set(follow.t, NULL, 'fw.up.mx', cut( follow.t[, fw.up.mx], breaks=follow.q, labels=c('<75pc','75-95pc','>95pc') )		)
	setnames(follow.t, c('Patient','fw.up.mx'), c('t.Patient','fw.up'))
	follow.t	
}
######################################################################################
project.athena.Fisheretal.X.followup.compareCD4toVL<- function(df.tpairs, df.immu, clumsm.info)
{
	follow		<- subset(df.immu, select=c(Patient, PosCD4))
	follow		<- merge( data.table(Patient= df.tpairs[, unique(t.Patient)]),  follow, by='Patient' )
	follow.cd4	<- follow[, list(nCD4=length(unique(PosCD4))), by='Patient']	
	follow		<- subset(df.viro, select=c(Patient, PosRNA))
	follow		<- merge( data.table(Patient= df.tpairs[, unique(t.Patient)]),  follow, by='Patient' )
	follow.vl	<- follow[, list(nVL=length(unique(PosRNA))), by='Patient']
	follow		<- merge(follow.cd4, follow.vl, by='Patient')
	#plot(follow[,nCD4], follow[,nVL], pch=18)
	#abline(a=0,b=1, col='red')
	#cat(paste('\n number of CD4 and VL tests=',follow[,sum(nCD4)], follow[,sum(nVL)]))
	
	follow		<- subset(df.immu, select=c(Patient, PosCD4))
	follow		<- merge( data.table(Patient= df.tpairs[, unique(t.Patient)]),  follow, by='Patient' )
	set(follow, NULL, 'PosCD4', hivc.db.Date2numeric(follow[,PosCD4]))
	follow.cd4	<- follow[, list( fw.up.mx.cd4= ifelse(length(PosCD4)>1, max(diff(PosCD4)), NA_real_) ),by='Patient']
	
	follow		<- subset(df.viro, select=c(Patient, PosRNA))
	follow		<- merge( data.table(Patient= df.tpairs[, unique(t.Patient)]),  follow, by='Patient' )
	set(follow, NULL, 'PosRNA', hivc.db.Date2numeric(follow[,PosRNA]))
	follow.vl	<- follow[, list( fw.up.mx.vl= ifelse(length(PosRNA)>1, max(diff(PosRNA)), NA_real_) ),by='Patient']
	follow		<- merge(follow.cd4, follow.vl, by='Patient')
	
	#plot(follow[,fw.up.mx.cd4], follow[,fw.up.mx.vl], pch=18)
	#abline(a=0,b=1, col='red')
	#subset( follow, 2*fw.up.mx.cd4<= fw.up.mx.vl )
	
	follow.select	<- subset( follow, fw.up.mx.cd4>1 & 2*fw.up.mx.cd4<= fw.up.mx.vl )
	follow.select	<- merge( follow.select, subset( clumsm.info, select=c(Patient, cluster) ), by='Patient' )
	setkey(follow.select, cluster)
	follow.select
}
######################################################################################
project.athena.Fisheretal.Y.brlweight<- function(Y.rawbrl, Y.rawbrl.linked	, Y.rawbrl.unlinked, linked.x= 0, plot.file=NA)
{
	#	check if Exp model would be reasonable
	require(MASS)
	unlinked.x	<- quantile( Y.rawbrl.unlinked[,brl], p=1e-3 )	 
	rawbrl.l	<- subset( Y.rawbrl.linked, brl>=0 )[, sort(brl)]		#allow for one mutation
	rawbrl.exp	<- fitdistr(rawbrl.l, "exponential") 
	#plot(log(rawbrl.l[,y]))
	#rawbrl.h	<- hist(Y.rawbrl.linked[,brl], breaks=1e3, freq=0, xlim=c(0,0.05))
	#tmp			<- seq(min(Y.rawbrl.linked[,brl]), max(Y.rawbrl.linked[,brl]), length.out=1024) 
	#lines( tmp, dexp(tmp, rate = rawbrl.exp$estimate), col='red')
	Y.brl		<- Y.rawbrl
	Y.brl[, wbrl:= brl]
	tmp			<- pexp(unlinked.x, rate = rawbrl.exp$estimate)
	set(Y.brl, NULL, 'wbrl', pexp(Y.brl[,brl], rate = rawbrl.exp$estimate, lower.tail=FALSE) / tmp)		
	set(Y.brl, Y.brl[,which(brl>unlinked.x)], 'wbrl', 0 )
	Y.brl[, w2brl:= pgamma(Y.brl[,brl], shape=2, rate = rawbrl.exp$estimate, lower.tail=FALSE)]		
	#tmp			<- seq(min(Y.rawbrl.linked[,brl]), max(Y.rawbrl.linked[,brl]), length.out=1024)
	#plot(tmp, pgamma(tmp, shape=2, rate = rawbrl.exp$estimate, lower.tail=FALSE), type='l')
	#lines(tmp, pgamma(tmp, shape=1, rate = rawbrl.exp$estimate, lower.tail=FALSE), col='red')
	if(!is.na(plot.file))
	{
		pdf(plot.file, w=5, h=5)
		cat(paste('\nplot to file',plot.file))
		require(RColorBrewer)
		par(mar=c(4,4,0.5,0.5))
		cols		<- brewer.pal(3, 'Set1')[c(1,3,2)]
		legend.txt	<- c('same host', 'potential transmission pairs', 'branch length weight')	
		xlim		<- c(0, max(Y.rawbrl[, max(brl)],  Y.rawbrl.linked[, max(brl)])*1.1 )
		tmp			<- seq(from=xlim[1], to=xlim[2], by=diff(xlim)/200)
		hist( Y.rawbrl.linked[, brl], breaks=tmp , col=my.fade.col(cols[1],0.5), add=0, freq=0, xlab='branch length', main='', ylim=c(0,100) )
		hist( Y.rawbrl[, brl], breaks=tmp, col=my.fade.col(cols[2],0.5), freq=0, add=1 )
		tmp			<- seq(min(Y.rawbrl.linked[,brl]), max(Y.rawbrl.linked[,brl]), length.out=1024)
		tmp2		<- pgamma(tmp, shape=2, rate = rawbrl.exp$estimate, lower.tail=FALSE) 
		lines(tmp, tmp2/(sum(tmp2)*diff(tmp)[1]), col=cols[3], lwd=2)		
		legend('topright', bty='n', fill= cols, border=NA, legend=legend.txt)		
		dev.off()
	}	
	Y.brl		<- subset( Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, w2brl) )
	setnames(Y.brl, 'w2brl', 'brl')
	Y.brl		
}
######################################################################################
project.athena.Fisheretal.Y.rawbrl<- function(df.tpairs, indir, insignat, indircov, infilecov, infiletree, save.file=NA, plot.file=NA)
{
	if(!is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(is.na(save.file) || inherits(readAttempt, "try-error"))
	{
		require(adephylo)
		#
		#	load precomputed true pos pairs and true neg pairs for full population
		#	
		argv				<<- hivc.cmd.preclustering(indir, infiletree, insignat, indircov, infilecov, resume=1)				 
		argv				<<- unlist(strsplit(argv,' '))
		msm.pre				<- hivc.prog.get.clustering.precompute()
		msm.linked			<- msm.pre$ph.linked
		msm.unlinked.bytime	<- do.call('rbind',msm.pre$unlinked.bytime)
		#
		#	select any true negative transmitter for every seroconverter in the denominator population
		#
		tmp					<- unique( subset( df.tpairs, select='FASTASampleCode') ) 
		setnames(tmp, 'FASTASampleCode', 'sc.FASTASampleCode')	
		setnames(msm.unlinked.bytime, 'query.FASTASampleCode','sc.FASTASampleCode')
		msm.unlinked.bytime	<- merge(msm.unlinked.bytime, tmp, by='sc.FASTASampleCode')
		setnames(msm.unlinked.bytime, c('Patient','FASTASampleCode','DateDied'), c('t.Patient','t.FASTASampleCode','t.DateDied'))	
		cat(paste('\nfound viral sequences that are not linked to sequences of denominator group based on Died<NegT, n=',msm.unlinked.bytime[, length(unique(t.FASTASampleCode))]))
		cat(paste('\nfound unlinked pairs based on Died<NegT, n=',nrow(msm.unlinked.bytime)))
		#
		#	select any other within host sequence for every individual in the denominator population as 'true positive'
		#
		tmp					<- unique( subset( df.tpairs, select=Patient ) )
		msm.linked			<- merge( msm.linked, tmp, by='Patient' )
		cat(paste('\nfound patients in denominator group with multiple sequences that can be taken as truly linked, n=',msm.linked[,length(unique(Patient))]))
		cat(paste('\nfound seq in denominator group with multiple sequences that can be taken as truly linked, n=',msm.linked[,length(unique(t.FASTASampleCode))]))
		msm.linked			<- msm.linked[, {
					tmp<- combn(length(FASTASampleCode),2)
					list( t.FASTASampleCode= FASTASampleCode[tmp[1,]], FASTASampleCode=FASTASampleCode[tmp[2,]] )					
				},by='Patient']
		cat(paste('\nfound viral sequences pairs of patients in denominator group that can be considered linked, n=',nrow(msm.linked)))
		#
		#	load full MLE phylogeny and patristic distances of full MLE phylogeny
		#
		file	<- paste(indir, '/', infiletree, '_', gsub('/',':',insignat),".R",sep='')
		load(file)	#loads ph	
		file	<- paste(indir, '/', infiletree, '_', gsub('/',':',insignat),"_distTips.R",sep='')
		load(file)	#loads brl	
		#require(adephylo)
		#brl		<- distTips(ph , method='patristic')
		#save(brl, file=file)
		#
		#	compute branch lengths between truly linked and truly unlinked
		#
		brl.n				<- attr(brl,'Size')
		#	get tip indices in ph$tip.label to determine index of brl between the two tips
		tmp					<- unique( subset(msm.unlinked.bytime, select=sc.FASTASampleCode) )	
		tmp					<- tmp[, list(sc.i=match(sc.FASTASampleCode, ph$tip.label)), by='sc.FASTASampleCode']
		msm.unlinked.bytime	<- merge(msm.unlinked.bytime, tmp, by='sc.FASTASampleCode')
		tmp					<- unique( subset(msm.unlinked.bytime, select=t.FASTASampleCode) )
		tmp					<- tmp[, list(sc.t=match(t.FASTASampleCode, ph$tip.label)), by='t.FASTASampleCode']
		msm.unlinked.bytime	<- merge(msm.unlinked.bytime, tmp, by='t.FASTASampleCode')
		#	make sure that sc.i is always smaller than sc.t
		tmp					<- msm.unlinked.bytime[, which(sc.t<sc.i)]
		tmp2				<- msm.unlinked.bytime[tmp, sc.t]
		set(msm.unlinked.bytime, tmp, 'sc.t', msm.unlinked.bytime[tmp,sc.i])
		set(msm.unlinked.bytime, tmp, 'sc.i', tmp2)
		#	get brl between true negatives
		msm.unlinked.bytime	<- msm.unlinked.bytime[, 	list( brl= brl[ my.lower.tri.index(brl.n, sc.t, sc.i) ], t.Patient=t.Patient) ,by=c('sc.FASTASampleCode','t.FASTASampleCode')]
		#	get tip indices in ph$tip.label to determine index of brl between the two tips
		tmp					<- unique( subset(msm.linked, select=FASTASampleCode) )	
		tmp					<- tmp[, list(sc.i=match(FASTASampleCode, ph$tip.label)), by='FASTASampleCode']
		msm.linked			<- merge(msm.linked, tmp, by='FASTASampleCode')
		tmp					<- unique( subset(msm.linked, select=t.FASTASampleCode) )	
		tmp					<- tmp[, list(sc.t=match(t.FASTASampleCode, ph$tip.label)), by='t.FASTASampleCode']
		msm.linked			<- merge(msm.linked, tmp, by='t.FASTASampleCode')
		#	make sure that sc.i is always smaller than sc.t
		tmp					<- msm.linked[, which(sc.t<sc.i)]
		tmp2				<- msm.linked[tmp, sc.t]
		set(msm.linked, tmp, 'sc.t', msm.linked[tmp,sc.i])
		set(msm.linked, tmp, 'sc.i', tmp2)
		#	get brl between true negatives
		msm.linked			<- msm.linked[, 	list( brl= brl[ my.lower.tri.index(brl.n, sc.t, sc.i) ], Patient=Patient) ,by=c('FASTASampleCode','t.FASTASampleCode')]
		#
		#	compute branch lengths between infected and all potential transmitters
		#
		df.tpairs.brl		<- df.tpairs
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
		df.tpairs.brl		<- df.tpairs.brl[, 	list( brl= brl[ my.lower.tri.index(brl.n, sc.t, sc.i) ]) ,by=c('FASTASampleCode','t.FASTASampleCode')]
		
		if(!is.na(save.file))
		{
			cat(paste('\nsave to file',save.file))
			save(file= save.file, df.tpairs.brl=df.tpairs.brl, msm.linked=msm.linked, msm.unlinked.bytime=msm.unlinked.bytime)
		}					
	}
	#
	#	plot brl for comparison
	#
	if(!is.na(plot.file))
	{
		pdf(plot.file, w=5, h=5)
		cat(paste('\nplot to file',plot.file))
		require(RColorBrewer)
		par(mar=c(4,4,0.5,0.5))
		cols		<- brewer.pal(3, 'Set1')
		legend.txt	<- c('same host', 'Died < last HIV- test', 'potential transmission pairs')	
		xlim		<- c(0, max( msm.linked[, max(brl)], msm.unlinked.bytime[, max(brl)], df.tpairs.brl[, max(brl)])*1.1 )
		tmp			<- seq(from=xlim[1], to=xlim[2], by=diff(xlim)/200)
		hist( msm.linked[, brl], breaks=tmp , col=my.fade.col(cols[1],0.5), add=0, freq=0, xlab='branch length', main='' )
		hist( msm.unlinked.bytime[, brl], breaks=tmp, col=my.fade.col(cols[2],0.5), freq=0, add=1 )
		hist( df.tpairs.brl[, brl], breaks=tmp, col=my.fade.col(cols[3],0.5), freq=0, add=1 )
		legend('topright', bty='n', fill= cols, border=NA, legend=legend.txt)
		
		tmp			<- msm.unlinked.bytime[, quantile(brl, p=c(1e-5, 1e-4, 0.001, 0.01))]
		ltys		<- seq_along(tmp)+1
		abline(v=tmp, col=cols[2], lty=ltys)
		legend('bottomright', bty='n', lty= c(ltys,NA), col=cols[2], border=NA, legend=c('1e-5', '1e-4', '1e-3', '1e-2',''))
		dev.off()
	}
	
	list(tpairs=df.tpairs.brl, linked=msm.linked, unlinked=msm.unlinked.bytime)	
}
######################################################################################
project.athena.Fisheretal.select.denominator<- function(indir, infile, insignat, indircov, infilecov, infiletree, adjust.AcuteByNegT=NA, adjust.NegT4Acute=NA, adjust.AcuteSelect=c('Yes'))
{
	#	fixed input args
	opt.brl			<- "dist.brl.casc" 
	thresh.brl		<- 0.096
	thresh.bs		<- 0.8					
	#
	# 	recent msm 
	#
	load(paste(indircov,'/',infilecov,'.R',sep=''))
	#	adjust Acute=='Maybe' by NegT 
	if(!is.na(adjust.AcuteByNegT))
	{
		tmp		<- which( df.all[, (is.na(isAcute) | isAcute=='No') & !is.na(NegT) & AnyPos_T1<=NegT+adjust.AcuteByNegT*365])
		cat(paste('\nmsm.recent: set Acute==Maybe when NegT is close to AnyPos_T1, n=',length(tmp)))
		set(df.all, tmp, 'isAcute', 'Maybe')
	}		
	msm.recent		<- subset( df.all, Sex=='M' & !Trm%in%c('OTH','IDU','HET','BLOOD','PREG','HETfa','NEEACC','SXCH') )
	setkey(msm.recent, isAcute)
	msm.recent		<- msm.recent[adjust.AcuteSelect,]
	cat(paste('\nmsm isAcute: #seq=',nrow(msm.recent),'#patient=',length(msm.recent[,unique(Patient)])))
	# 	recent msm seroconverters
	msm.recentsc	<- subset( msm.recent, !is.na(NegT))
	cat(paste('\nmsm isAcute & !is.na(NegT): #seq=',nrow(msm.recentsc),'#patient=',length(msm.recentsc[,unique(Patient)])))
	#
	#	load clustering msm
	#
	argv			<<- hivc.cmd.clustering.msm(indir, infiletree, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=1)
	argv			<<- unlist(strsplit(argv,' '))		
	msm				<- hivc.prog.get.clustering.MSM()
	clumsm.info		<- msm$df.cluinfo
	set(clumsm.info, NULL, 'PosSeqT', hivc.db.Date2numeric(clumsm.info[,PosSeqT]))
	#set(clumsm.info, NULL, 'DateBorn', hivc.db.Date2numeric(clumsm.info[,DateBorn]))
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
	#	fixup
	setkey(clumsm.info, Patient)
	tmp		<- clumsm.info[, list(ok= all(AnyPos_T1==min(AnyPos_T1))) ,by='Patient']
	for(x in subset(tmp, !ok)[, unique(Patient)])
	{
		z	<- clumsm.info[, which(Patient==x)]		
		set( clumsm.info, z, 'AnyPos_T1', clumsm.info[z, min(AnyPos_T1)])
	}
	#	adjust Acute=='Maybe' by NegT 
	if(!is.na(adjust.AcuteByNegT))
	{
		tmp		<- which( clumsm.info[, (is.na(isAcute) | isAcute=='No') & !is.na(NegT) & AnyPos_T1<=NegT+adjust.AcuteByNegT])
		cat(paste('\nclumsm.info: set Acute==Maybe when NegT is close to AnyPos_T1, n=',length(tmp)))
		set(clumsm.info, tmp, 'isAcute', 'Maybe')
	}
	#	adjust missing NegT when Acute=='Yes'
	if(!is.na(adjust.NegT4Acute))
	{		
		tmp		<- which( clumsm.info[, !is.na(isAcute) & isAcute=='Yes' & is.na(NegT)] )	
		cat(paste('\nset NegT for Acute==Yes and NegT missing to one year before AnyPos_T1, n=',length(tmp)))
		set(clumsm.info, tmp,'NegT', subset(clumsm.info, !is.na(isAcute) & isAcute=='Yes' & is.na(NegT) )[, AnyPos_T1-adjust.NegT4Acute] )					
	}	
	# 	recent individuals in cluster
	setkey(clumsm.info, isAcute)	
	clumsm.recent	<- clumsm.info[adjust.AcuteSelect,] 
	# 	recent seroconcerverters in cluster
	clumsm.recentsc	<- subset( clumsm.recent, !is.na(NegT))
	cat(paste('\nmsm clustering: #seq=',nrow(clumsm.info),'#patient=',length(clumsm.info[,unique(Patient)]),'#cluster=',length(clumsm.info[,unique(cluster)])))		
	cat(paste('\nmsm clustering isAcute: #seq=',nrow(clumsm.recent),'#patient=',length(clumsm.recent[,unique(Patient)]),'#cluster=',length(clumsm.recent[,unique(cluster)])))
	cat(paste('\nmsm clustering isAcute & !is.na(NegT): #seq=',nrow(clumsm.recentsc),'#patient=',length(clumsm.recentsc[,unique(Patient)]),'#cluster=',length(clumsm.recentsc[,unique(cluster)])))
	#
	#	make selection
	#
	list(clumsm.info=clumsm.info, df.select=clumsm.recent, clumsm.subtrees=msm$cluphy.subtrees, clumsm.ph=msm$cluphy)
}
######################################################################################
project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos<- function(clumsm.info, df.denom, any.pos.grace.yr= 0.5, select.if.transmitter.seq.unique=TRUE)
{
	df.transmitters	<- clumsm.info[J(df.denom[, unique(cluster)]),]
	df.transmitters	<- subset(df.transmitters, select=c(cluster, Patient, AnyPos_T1, FASTASampleCode))		
	setnames(df.transmitters, colnames(df.transmitters), paste('t.',colnames(df.transmitters),sep=''))
	# 	all possible sequence pairs with diagnosis criterium
	df.tpairs			<- df.denom[,	{
				tmp<- subset(df.transmitters, t.cluster==cluster & t.Patient!=Patient[1] & t.AnyPos_T1<=any.pos.grace.yr+AnyPos_T1[1])
				as.list(tmp)
			}, by=c('cluster','FASTASampleCode')]								
	df.tpairs			<- merge( subset(clumsm.info, select=c(FASTASampleCode, Patient)), df.tpairs, by='FASTASampleCode')
	cat(paste('\nnumber possible transmission sequence pairs, n=',nrow(df.tpairs)))
	# 	select unique Infected->pot Transmitter based on this criterium
	df.tpairs.stat		<- df.tpairs[,list(t.potpat.n=length(unique(t.Patient)), t.potpatseq.n=length(unique(t.FASTASampleCode)) ),by=c('cluster','Patient')]
	cat(paste('\ntable of #possible transmitter seqs for every infected patient'))
	tmp					<- table( df.tpairs.stat[,t.potpatseq.n] )
	print( tmp )
	print( cumsum( tmp / sum(tmp ) ) )	
	#	allows for multiple infected seq pointing to same transmitter seq
	if(select.if.transmitter.seq.unique)		
		df.tpairs.anypos	<- merge(df.tpairs, subset(df.tpairs.stat, t.potpatseq.n==1), by=c('cluster','Patient'))
	else
		df.tpairs.anypos	<- df.tpairs
	
	cat(paste('\nselected tpairs: total', nrow(df.tpairs.anypos),' recent #seq=',length(df.tpairs.anypos[,unique(FASTASampleCode)]),'#patient=',length(df.tpairs.anypos[,unique(Patient)]),'#cluster=',length(df.tpairs.anypos[,unique(cluster)])))
	subset( df.tpairs.anypos, select=c(cluster, Patient, FASTASampleCode, t.Patient, t.FASTASampleCode) )	
}
######################################################################################
project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos.ClosestSeqSim<- function(clumsm.info, df.denom, any.pos.grace.yr= 0.5)
{
	indir			<- paste(DATA,"tmp",sep='/')		
	infile			<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
	insignat		<- "Thu_Aug_01_17/05/23_2013"		
	#
	#	load sequences
	#
	file		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
	cat(paste("\nload sequences from file",file))
	load( file )	# loads seq.PROT.RT
	
	#	get potential transmitters
	df.transmitters				<- clumsm.info[J(df.denom[, unique(cluster)]),]
	df.transmitters				<- subset(df.transmitters, select=c(cluster, Patient, AnyPos_T1, FASTASampleCode))		
	setnames(df.transmitters, colnames(df.transmitters), paste('t.',colnames(df.transmitters),sep=''))
	# 	all possible sequence pairs with diagnosis criterium
	df.tpairs					<- df.denom[,	{
				tmp<- subset(df.transmitters, t.cluster==cluster & t.Patient!=Patient[1] & t.AnyPos_T1<=any.pos.grace.yr+AnyPos_T1[1])
				as.list(tmp)
			}, by=c('cluster','FASTASampleCode')]								
	df.tpairs					<- merge( subset(clumsm.info, select=c(FASTASampleCode, Patient, AnyPos_T1)), df.tpairs, by='FASTASampleCode')
	cat(paste('\nnumber possible transmission sequence pairs, n=',nrow(df.tpairs)))	
	#	compute genetic distances
	dummy						<- 0.
	tmp							<- sapply(seq_len(nrow(df.tpairs)),function(i)		.C("hivc_dist_ambiguous_dna", seq.PROT.RT[df.tpairs[i,FASTASampleCode], ], seq.PROT.RT[df.tpairs[i,t.FASTASampleCode], ], ncol(seq.PROT.RT), dummy )[[4]]		)
	df.tpairs[, seq.similar:= tmp]
	#	select closest transmitter seqs
	df.tpairs					<- df.tpairs[, .SD[seq.similar==max(seq.similar)], by=c('cluster','FASTASampleCode')]
	df.tpairs.stat				<- df.tpairs[,list(t.potpat.n=length(unique(t.Patient)), t.potpatseq.n=length(unique(t.FASTASampleCode)) ),by=c('cluster','Patient')]
	cat(paste('\ntable of #possible transmitter seqs for every infected patient'))
	tmp							<- table( df.tpairs.stat[,t.potpatseq.n] )
	print( tmp )
	print( cumsum( tmp / sum(tmp ) ) )	
	#	return uniquely identified pot transmitters
	df.tpairs.anypos.seqsim	<- merge(df.tpairs, subset(df.tpairs.stat, t.potpatseq.n==1), by=c('cluster','Patient'))		
	cat(paste('\nselected Patients with one potential transmitter', length(df.tpairs.anypos.seqsim[, unique(Patient)])))		
	subset( df.tpairs.anypos.seqsim, select=c(cluster, Patient, FASTASampleCode, t.Patient, t.FASTASampleCode) )		
}
######################################################################################
project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos.MinBrlMLETree<- function(df.denom, clumsm.subtrees, any.pos.grace.yr= 0.5)
{
	require(adephylo)
	#	get potential transmitters
	df.transmitters			<- clumsm.info[J(df.denom[, unique(cluster)]),]
	df.transmitters			<- subset(df.transmitters, select=c(cluster, Patient, AnyPos_T1, FASTASampleCode))		
	setnames(df.transmitters, colnames(df.transmitters), paste('t.',colnames(df.transmitters),sep=''))
	# 	all possible sequence pairs with diagnosis criterium
	df.tpairs				<- df.denom[,	{
				tmp<- subset(df.transmitters, t.cluster==cluster & t.Patient!=Patient[1] & t.AnyPos_T1<=any.pos.grace.yr+AnyPos_T1[1])
				as.list(tmp)
			}, by=c('cluster','FASTASampleCode')]								
	df.tpairs				<- merge( subset(clumsm.info, select=c(FASTASampleCode, Patient, AnyPos_T1)), df.tpairs, by='FASTASampleCode')
	cat(paste('\nnumber possible transmission sequence pairs, n=',nrow(df.tpairs)))	
	#
	df.tpairs.stat			<- df.tpairs[,list(t.potpat.n=length(unique(t.Patient)), t.potpatseq.n=length(unique(t.FASTASampleCode)) ),by=c('cluster','Patient')]
	cat(paste('\ntable of #possible transmitter seqs for every infected patient that satisfy AnyPos'))
	tmp						<- table( df.tpairs.stat[,t.potpatseq.n] )
	print( tmp )
	print( cumsum( tmp / sum(tmp ) ) )		
	#	compute branch length distances
	clumsm.subtrees.clu		<- data.table(cluster=as.numeric(names(clumsm.subtrees)))	
	tmp						<- setdiff( df.tpairs[,unique(cluster)], clumsm.subtrees.clu[, unique(cluster)] )
	cat(paste('\nnumber of clusters missing for analysis of pot transmitters, n=',length(tmp)))
	cat(paste('\nnumber of denominator individuals for which clusters are available, n=',subset( df.tpairs, !cluster%in%tmp )[, length(unique(Patient))]))
	if(length(tmp))	cat(paste('\nnumber of denominator individuals  for which clusters are missing, n=',subset( df.tpairs, cluster%in%tmp )[, length(unique(Patient))]))
	#
	df.tpairs				<- merge(df.tpairs, clumsm.subtrees.clu, by='cluster')
	cat(paste('\nnumber possible transmission sequence pairs for which dated clusters are available, n=',nrow(df.tpairs)))
	tmp						<- lapply( df.tpairs[,unique(cluster)], function(clu)
			{
				cluphy.subtree				<- clumsm.subtrees[[ which(clumsm.subtrees.clu[,cluster==clu]) ]]								 
				tmp							<- subset(df.tpairs, cluster==clu)
				dist						<- as.matrix( distTips(cluphy.subtree , method='patristic') ) 
				dist						<- sapply(seq_len(nrow(tmp)), function(i)		dist[ tmp[i,FASTASampleCode], tmp[i, t.FASTASampleCode]]		)
				data.table(FASTASampleCode= tmp[,FASTASampleCode,], t.FASTASampleCode=tmp[, t.FASTASampleCode], patristic=dist)					
			})
	tmp						<- do.call('rbind', tmp)
	df.tpairs				<- merge(df.tpairs, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))
	#	select closest transmitter seqs
	df.tpairs				<- df.tpairs[, .SD[patristic==min(patristic)], by=c('cluster','FASTASampleCode')]
	df.tpairs.stat			<- df.tpairs[,list(t.potpat.n=length(unique(t.Patient)), t.potpatseq.n=length(unique(t.FASTASampleCode)) ),by=c('cluster','Patient')]
	cat(paste('\ntable of #possible transmitter seqs for every infected patient that satisfy AnyPos and MLE.MinBrl'))
	tmp						<- table( df.tpairs.stat[,t.potpatseq.n] )
	print( tmp )
	print( cumsum( tmp / sum(tmp ) ) )	
	#	return uniquely identified pot transmitters
	ans						<- merge(df.tpairs, subset(df.tpairs.stat, t.potpatseq.n==1), by=c('cluster','Patient'))		
	cat(paste('\nselected Patients with one potential transmitter', length(ans[, unique(Patient)])))		
	subset( ans, select=c(cluster, Patient, FASTASampleCode, t.Patient, t.FASTASampleCode) )		
}	
######################################################################################
project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos.MinBrl<- function(clumsm.info, df.denom, cluphy.subtrees, cluphy.info, any.pos.grace.yr= 0.5)
{
	require(adephylo)
	#	get potential transmitters
	df.transmitters				<- clumsm.info[J(df.denom[, unique(cluster)]),]
	df.transmitters				<- subset(df.transmitters, select=c(cluster, Patient, AnyPos_T1, FASTASampleCode))		
	setnames(df.transmitters, colnames(df.transmitters), paste('t.',colnames(df.transmitters),sep=''))
	# 	all possible sequence pairs with diagnosis criterium
	df.tpairs					<- df.denom[,	{
				tmp<- subset(df.transmitters, t.cluster==cluster & t.Patient!=Patient[1] & t.AnyPos_T1<=any.pos.grace.yr+AnyPos_T1[1])
				as.list(tmp)
			}, by=c('cluster','FASTASampleCode')]								
	df.tpairs					<- merge( subset(clumsm.info, select=c(FASTASampleCode, Patient, AnyPos_T1)), df.tpairs, by='FASTASampleCode')
	cat(paste('\nnumber possible transmission sequence pairs, n=',nrow(df.tpairs)))	
	#	compute branch length distances
	setkey(cluphy.info, BEASTlabel)
	cluphy.subtrees.clu		<- data.table(cluster=as.numeric(names(cluphy.subtrees)))	
	tmp						<- setdiff( df.tpairs[,unique(cluster)], cluphy.subtrees.clu[, unique(cluster)] )
	cat(paste('\nnumber of clusters missing for analysis of pot transmitters, n=',length(tmp)))
	if(length(tmp))
		cat(paste('\nnumber of pot transmitters for which clusters are missing, n=',subset( df.tpairs, cluster%in%tmp )[, length(unique(t.Patient))]))		
	df.tpairs				<- merge(df.tpairs, cluphy.subtrees.clu, by='cluster')
	cat(paste('\nnumber possible transmission sequence pairs for which dated clusters are available, n=',nrow(df.tpairs)))
	tmp			<- lapply( df.tpairs[,unique(cluster)], function(clu)
			{
				cluphy.subtree				<- cluphy.subtrees[[ which(cluphy.subtrees.clu[,cluster==clu]) ]]		
				cluphy.subtree$tip.label	<- cluphy.info[ cluphy.subtree$tip.label, ][,FASTASampleCode]				 
				tmp							<- subset(df.tpairs, cluster==clu)
				dist						<- as.matrix( distTips(cluphy.subtree , method='patristic') ) 
				dist						<- sapply(seq_len(nrow(tmp)), function(i)		dist[ tmp[i,FASTASampleCode], tmp[i, t.FASTASampleCode]]		)
				data.table(FASTASampleCode= tmp[,FASTASampleCode,], t.FASTASampleCode=tmp[, t.FASTASampleCode], patristic=dist)					
			})
	tmp			<- do.call('rbind', tmp)
	df.tpairs	<- merge(df.tpairs, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))
	#	select closest transmitter seqs
	df.tpairs		<- df.tpairs[, .SD[patristic==max(patristic)], by=c('cluster','FASTASampleCode')]
	df.tpairs.stat	<- df.tpairs[,list(t.potpat.n=length(unique(t.Patient)), t.potpatseq.n=length(unique(t.FASTASampleCode)) ),by=c('cluster','Patient')]
	cat(paste('\ntable of #possible transmitter seqs for every infected patient'))
	tmp				<- table( df.tpairs.stat[,t.potpatseq.n] )
	print( tmp )
	print( cumsum( tmp / sum(tmp ) ) )	
	#	return uniquely identified pot transmitters
	ans				<- merge(df.tpairs, subset(df.tpairs.stat, t.potpatseq.n==1), by=c('cluster','Patient'))		
	cat(paste('\nselected Patients with one potential transmitter', length(ans[, unique(Patient)])))		
	subset( ans, select=c(cluster, Patient, FASTASampleCode, t.Patient, t.FASTASampleCode) )	
}
######################################################################################
project.athena.Fisheretal.get.data.for.selection<- function(df.tpairs, clu.indir, clu.infile, clu.insignat, clu.infilexml.opt, clu.infilexml.template, indircov, infilecov)
{
	file.cov				<- paste(indircov,"/",infilecov,".R",sep='')
	file.viro				<- paste(indircov,"/ATHENA_2013_03_Viro.R",sep='/')
	file.immu				<- paste(indircov,"/ATHENA_2013_03_Immu.R",sep='/')
	file.treatment			<- paste(indircov,"/ATHENA_2013_03_Regimens.R",sep='/')		
	#	load patient demographic data	(loads df.all)
	load(file.cov)				
	#	load patient RNA
	load(file.viro)
	df.viro				<- df				
	#	load patient CD4				
	load(file.immu)
	df.immu				<- df
	#	load patient regimen
	load(file.treatment)
	df.treatment		<- df		
	#
	#	collect files containing dated cluster phylogenies
	files		<- list.files(clu.indir)
	if(!length(files))	stop('no input files matching criteria')
	files		<- files[ sapply(files, function(x) grepl(clu.infile, x, fixed=1) & grepl(gsub('/',':',clu.insignat), x, fixed=1) & grepl(clu.infilexml.opt, x, fixed=1) & grepl(clu.infilexml.template,x, fixed=1) & grepl('_cluposterior_[0-9]+',x) & grepl('R$',x) ) ]		
	if(!length(files))	stop('no input files matching criteria')
	tmp			<- regmatches( files, regexpr('_cluposterior_[0-9]+',files)) 
	cluster		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
	file.info	<- data.table(file=files, cluster=cluster)
	setkey(file.info, cluster)
	#	subset of clusters: those in df.tpairs.select	
	if(!is.null(df.tpairs))
	{
		tmp					<- setdiff( df.tpairs[,unique(cluster)], file.info[, unique(cluster)] )
		cat(paste('\nnumber of clusters missing for analysis of pot transmitters, n=',length(tmp)))
		file.info	<- merge(file.info, unique(subset(df.tpairs, select=cluster)), by='cluster')
		df.tpairs.reduced	<- merge(df.tpairs, subset(file.info, select=cluster), by='cluster')
		cat(paste('\nnumber of remaining Patients with one potential transmitter', length(df.tpairs.reduced[, unique(Patient)])))
	}		
	#	combine dated cluster phylogenies
	clu			<- hivc.beast2out.combine.clu.trees(clu.indir, file.info)		
	
	list(clu=clu, df.all=df.all, df.viro=df.viro, df.immu=df.immu, df.treatment=df.treatment, df.tpairs=df.tpairs.reduced)
}
######################################################################################
project.athena.Fisheretal.plot.selected.transmitters<- function(df.all, df.immu, df.viro, df.treatment, df.tpairs, cluphy, cluphy.info, cluphy.subtrees, cluphy.map.nodectime, outfile)
{
	require(RColorBrewer)
	beastlabel.idx.clu			<- 1
	beastlabel.idx.hivn			<- 2
	beastlabel.idx.hivd			<- 3
	beastlabel.idx.hivs			<- 4
	beastlabel.idx.samplecode	<- 6
	beastlabel.idx.rate			<- NA			
	#
	#	prepare df.tpairs.select for plotting
	#
	if(!is.null(df.tpairs))
	{
		df.tpairs.plot		<- subset( df.tpairs, select=c(cluster, FASTASampleCode, t.FASTASampleCode))
		setkey(df.tpairs.plot, cluster)
		tmp					<- df.tpairs.plot[, list(col=rev(brewer.pal( max(3,length(FASTASampleCode)), 'Dark2' )[seq_along(FASTASampleCode)])), by='cluster']
		df.tpairs.plot[, col:=tmp[,col]]
		df.tpairs.plot[, pch:=0]
		df.tpairs.plot[, t.pch:=5]
		df.tpairs.plot[, cex:=1]
		df.tpairs.plot[, t.cex:=1.4]
		tmp					<- subset(df.tpairs.plot, select=c(cluster, t.FASTASampleCode, col, t.pch, t.cex))
		setnames(tmp, c('t.FASTASampleCode', 't.pch', 't.cex'), c('FASTASampleCode', 'pch', 'cex'))
		df.tpairs.plot		<- rbind( subset(df.tpairs.plot, select=c(cluster, FASTASampleCode, col, pch, cex)), tmp )
		setkey(df.tpairs.plot, 'cluster','FASTASampleCode')			
	}		
	tmp							<- setdiff( df.tpairs.plot[, unique(FASTASampleCode)], cluphy$tip.label )
	cat(paste('cannot find tip labels for df.tpairs, n=',length(tmp)))
	df.tpairs.plot				<- merge( df.tpairs.plot, data.table(FASTASampleCode=cluphy$tip.label), by='FASTASampleCode' )
	#	determine pdf.xlim: subset of timeline that contains the MRCA of all clusters
	cluphy.root.ctime			<- cluphy.info[1,TipT]-node.depth.edgelength(cluphy)[1]
	cluphy.tip.ctime			<- cluphy.info[, TipT]
	cluphy.subtrees.root.ctime	<- sapply(seq_along(cluphy.subtrees), function(i)
			{
				tmp		<- hivc.treeannotator.tiplabel2df(cluphy.subtrees[[i]], beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
				tmp[1,TipT]-node.depth.edgelength(cluphy.subtrees[[i]])[1]
			})		
	cluphy.subtrees.root.ctime	<- min(cluphy.subtrees.root.ctime)		
	end.ctime					<- 2013.3
	pdf.xlim					<- c( cluphy.subtrees.root.ctime-cluphy.root.ctime-1, ceiling(end.ctime)-cluphy.root.ctime )
	pdf.xlim					<- pdf.xlim+c(0,diff(pdf.xlim)*0.35)
	#	plot	
	cat(paste('\nplot to file=',outfile))
	dummy				<- hivc.beast2out.plot.cluster.trees(df.all, df.immu, df.viro, df.treatment, cluphy, cluphy.root.ctime, cluphy.tip.ctime, ph.prob=NA, df.node.ctime=cluphy.map.nodectime, df.rates=NULL, df.tips=df.tpairs.plot, end.ctime=end.ctime,  cex.nodelabel=0.5,  cex.tiplabel=0.5,  file=outfile,  pdf.width=7, pdf.height=140, pdf.xlim=pdf.xlim)	
}
######################################################################################
project.athena.Fisheretal.X.viro<- function(df.tpairs, df.viro, t.period=0.25, lRNA.cutoff= log10(1e3))
{
	viro	<- subset( df.viro, select=c(Patient, PosRNA, lRNA) )
	setnames(viro, 'Patient','t.Patient')
	viro	<- merge(viro, unique(subset(df.tpairs, select=t.Patient)), by='t.Patient')	
	set(viro, NULL, 'PosRNA', hivc.db.Date2numeric(viro[,PosRNA]))
	tmp		<- subset(viro, select=c(t.Patient, PosRNA))[, list(ts=min(PosRNA), te=max(PosRNA)), by='t.Patient']
	set(tmp, NULL, 'ts', tmp[, floor(ts) + floor( (ts%%1)*100 %/% (t.period*100) ) * t.period] )
	set(tmp, NULL, 'te', tmp[, floor(te) + floor( (te%%1)*100 %/% (t.period*100) ) * t.period] )
	tmp		<- tmp[, list(t= seq(ts, te, by=t.period)),by='t.Patient']
	
	setnames(tmp, 't.Patient', 'tt.Patient')
	setkey(tmp, tt.Patient)
	viro	<- viro[, {
				z		<- subset(tmp, tt.Patient==t.Patient[1])	
				if(length(PosRNA)<2)
				{
					y	<- rep(NA, nrow(z))
					y[z[,which( PosRNA<=t+t.period )[1]]]<- lRNA
				}
				else
					y	<- approx(PosRNA , lRNA, xout=z[,t]+t.period/2, yleft=NA_real_, yright=NA_real_, rule=2)$y
				list(t=z[,t], lRNA=y)
			},by='t.Patient']
	if(!is.na(lRNA.cutoff))
	{
		set(viro, viro[, which(lRNA<lRNA.cutoff)], 'lRNA', 0.)
		set(viro, viro[, which(lRNA>=lRNA.cutoff)], 'lRNA', 1.)
		set(viro, NULL, 'lRNA', as.integer(viro[,lRNA]))		
	}
	viro
}
######################################################################################
project.athena.Fisheretal.Y.coal.and.inf<- function(df.tpairs, cluphy, cluphy.info, cluphy.map.nodectime, t.period= 0.25 )
{
	setkey(df.tpairs, cluster)	
	#	select tpairs for which dated phylogenies are available
	tmp					<- setdiff( df.tpairs[,unique(cluster)], cluphy.info[, unique(cluster)] )
	cat(paste('\nnumber of clusters missing for analysis of pot transmitters, n=',length(tmp)))
	df.tpairs			<- df.tpairs[J(cluphy.info[, unique(cluster)]),]	
	cat(paste('\nnumber of pot transmitters for which dated phylogenies are available, n=',df.tpairs[,length(unique(Patient))]))
	#	compute mrcas of tpairs 
	tmp					<- df.tpairs[,	list(node=getMRCA(cluphy, c(FASTASampleCode, t.FASTASampleCode))), by=c('FASTASampleCode','t.FASTASampleCode')]
	if(nrow(tmp)!=nrow(df.tpairs))	stop('unexpected length of tmp')
	df.tpairs			<- merge(df.tpairs, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))
	#	compute posterior prob that coalescence is after last NegT and before TipT of individual in denominator population
	tmp					<- unique( subset( clumsm.info, select=c(Patient, NegT)) )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	df.tpairs			<- merge(df.tpairs, tmp, by='t.Patient')
	tmp					<- subset( df.tpairs, select=c(cluster, FASTASampleCode, t.FASTASampleCode, node,  t.NegT) )
	tmp					<- merge(tmp, subset(cluphy.info, select=c(FASTASampleCode, AnyPos_T1)), by='FASTASampleCode')
	tmp					<- merge( cluphy.map.nodectime, tmp, by=c('cluster','node'))	
	coal				<- tmp[,  {
									z		<- 1-approx(q ,cdf, xout=c(t.NegT[1], AnyPos_T1[1]), yleft=0., yright=1., rule=2)$y									
									list(coal.after.t.NegT= z[1], coal.after.i.AnyPos_T1=z[2])
									} , by=c('cluster','FASTASampleCode','t.FASTASampleCode')]
	#	earliest time of infection: either lowest coal or NegT(infected)
	tmp					<- merge( subset( df.tpairs, select=c(cluster,  node, FASTASampleCode) ), subset( cluphy.info, select=c(FASTASampleCode, NegT, AnyPos_T1)), by='FASTASampleCode' )
	inf					<- merge(tmp, cluphy.map.nodectime[, list(min.coalT= min(q)), by=c('cluster','node')], by=c('cluster','node'))
	tmp					<- inf[,which(is.na(NegT))]
	set(inf, tmp, 'NegT', inf[tmp, min.coalT])
	setnames(inf, 'NegT', 'InfT')	
	#	get list of time periods for every infected
	set(inf, NULL, 'InfT', inf[, floor(InfT) + floor( (InfT%%1)*100 %/% (t.period*100) ) * t.period] )
	set(inf, NULL, 'AnyPos_T1', inf[, floor(AnyPos_T1) + floor( (AnyPos_T1%%1)*100 %/% (t.period*100) ) * t.period] )
	inf					<- inf[, list(InfT= seq(InfT, AnyPos_T1, by=t.period), cluster=cluster, node=node),by='FASTASampleCode']
	#	evaluated coal at InfT
	setkey(cluphy.map.nodectime, cluster, node)
	setnames(inf, c('cluster','node'), c('i.cluster','i.node'))
	inf					<- inf[, {
										tmp  <- subset(cluphy.map.nodectime, cluster==i.cluster[1] & node==i.node[1])
										InfV <- approx(tmp[,q] ,tmp[,cdf], xout=InfT+t.period/2, yleft=0., yright=1., rule=2)$y
										list(t=InfT, c.inf=InfV)
									},by='FASTASampleCode']
	#
	df.tpairs			<- merge(df.tpairs, coal, by=c('cluster','FASTASampleCode','t.FASTASampleCode'))			
	if(nrow(subset( df.tpairs, !is.na(t.NegT) & is.na(coal.after.t.NegT) )))	stop('unexpected NA')
	set(df.tpairs, which(df.tpairs[,is.na(coal.after.t.NegT)]), 'prob.coal.after.t.NegT', 1.) 
	#	filter df.tpairs based on coalescence
	df.tpairs.excluded	<- subset(df.tpairs, coal.after.t.NegT<=EPS)
	df.tpairs			<- subset(df.tpairs, coal.after.t.NegT>EPS)	
	cat(paste('\nnumber of tpairs for which t.prob.coal.after.NegT==0, n=',length(df.tpairs.excluded[,unique(Patient)])))
	cat(paste('\nnumber of tpairs for which t.prob.coal.after.NegT>0, n=',length(df.tpairs[,unique(Patient)])))
	coal				<- subset(df.tpairs, select=c(FASTASampleCode, t.FASTASampleCode, coal.after.t.NegT, coal.after.i.AnyPos_T1))
	df.tpairs			<- subset(df.tpairs, select=c(cluster, Patient, FASTASampleCode, t.Patient, t.FASTASampleCode))	
	inf					<- merge( subset(df.tpairs, select=c(FASTASampleCode, t.FASTASampleCode)), inf, by='FASTASampleCode')	
	list(df.tpairs=df.tpairs, inf=inf, coal=coal, exluded=df.tpairs.excluded)		
}
######################################################################################
project.athena.Fisheretal.exact.repro<- function()
{
	require(data.table)
	require(ape)
	stop()
	indir					<- paste(DATA,"tmp",sep='/')		
	indircov				<- paste(DATA,"derived",sep='/')
	outdir					<- paste(DATA,"fisheretal",sep='/')
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"
	if(1)
	{
		infile					<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs100",sep="_")
		insignat				<- "Thu_Aug_01_17/05/23_2013"	
		clu.indir				<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/beast/beast2_140201'
		clu.infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_seroneg-130"
		clu.insignat			<- "Tue_Aug_26_09:13:47_2013"
		clu.infilexml.template	<- "sasky_sdr06"
		clu.infilexml.opt		<- "rsu815"				
	}
	if(0)
	{
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"				
		clu.indir				<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/sasky_sdr06_-DR-RC-SH+LANL_alsu50'
		clu.infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		clu.insignat			<- "Wed_Dec_18_11:37:00_2013"
		clu.infilexml.opt		<- "alsu50"
		clu.infilexml.template	<- "sasky_sdr06"		
	}	
	#
	#	select infected individuals and return in df.select
	#
	tmp				<- project.athena.Fisheretal.select.denominator(indir, infile, insignat, indircov, infilecov, infiletree, adjust.AcuteByNegT=0.75, adjust.NegT4Acute=1, adjust.AcuteSelect=c('Yes','Maybe'))
	df.denom		<- tmp$df.select
	clumsm.subtrees	<- tmp$clumsm.subtrees
	clumsm.info		<- tmp$clumsm.info
	setkey(clumsm.info, cluster)
	#
	#	select potential transmitters on MLE tree
	#
	df.tpairs		<- project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos.MinBrlMLETree(df.denom, clumsm.subtrees, any.pos.grace.yr= 0.5)
	tmp				<- merge( subset(df.tpairs, select=Patient), subset(clumsm.info, select=c(Patient, AnyPos_T1)), by='Patient' )
	outfile			<- paste(outdir,'/',infile, '_', clu.infilexml.template, '_', clu.infilexml.opt, '_', gsub('/',':',insignat), '_', 'nrecentlyinfected', '.pdf',sep='')
	pdf(file=outfile, w=5, h=5)
	par(mar=c(3,5,0.5,0.5))
	barplot( table( tmp[, round(AnyPos_T1)] ), ylab="# recently infected\n with unique potential transmitter" )
	dev.off()
	#
	#	get clinical data and dated phylogenies for potential transmitters
	#	
	tmp						<- project.athena.Fisheretal.get.data.for.selection(df.tpairs, clu.indir, clu.infile, clu.insignat, clu.infilexml.opt, clu.infilexml.template, indircov, infilecov)	
	df.tpairs				<- tmp$df.tpairs
	df.all					<- tmp$df.all
	df.viro					<- tmp$df.viro
	df.immu					<- tmp$df.immu
	df.treatment			<- tmp$df.treatment
	cluphy					<- tmp$clu$cluphy
	cluphy.subtrees			<- tmp$clu$cluphy.subtrees
	cluphy.info				<- tmp$clu$cluphy.info
	cluphy.map.nodectime	<- tmp$clu$cluphy.map.nodectime	
	#	
	#	plot MLE tree
	#
	tmp								<- merge( data.table( cluster=as.numeric( names(clumsm.subtrees) ), clu.i=seq_along(clumsm.subtrees) ), df.tpairs, by='cluster' )
	df.tpairs.subtrees 				<- lapply( tmp[, unique(clu.i)], function(i)	clumsm.subtrees[[i]] )
	df.tpairs.mleph					<- hivc.clu.polyphyletic.clusters(cluphy.subtrees=df.tpairs.subtrees)$cluphy
	
	df.tpairs.info	<- merge( data.table(FASTASampleCode= df.tpairs.mleph$tip.label), subset(clumsm.info, select=c(FASTASampleCode, Patient, cluster, AnyPos_T1)), by='FASTASampleCode' )
	df.tpairs.info[, tiplabel:='']
	tmp				<- which( df.tpairs.info[, FASTASampleCode]%in%df.tpairs[, FASTASampleCode] )
	set(df.tpairs.info, tmp, 'tiplabel', paste('I',df.tpairs.info[tmp,tiplabel],sep='') )
	tmp				<- which( df.tpairs.info[, FASTASampleCode]%in%df.tpairs[, t.FASTASampleCode] )
	set(df.tpairs.info, tmp, 'tiplabel', paste('T',df.tpairs.info[tmp,tiplabel],sep='') )	
	set(df.tpairs.info, NULL, 'tiplabel', paste(df.tpairs.info[,tiplabel],'_',df.tpairs.info[,Patient],'_clu=',df.tpairs.info[,cluster],'_d=',df.tpairs.info[,AnyPos_T1],sep='') )
	setkey(df.tpairs.info, FASTASampleCode)
	df.tpairs.mleph$tip.label		<- df.tpairs.info[df.tpairs.mleph$tip.label, ][, tiplabel]
	outfile							<- paste(outdir,'/',infile, '_', clu.infilexml.template, '_', clu.infilexml.opt, '_', gsub('/',':',insignat), '_', 'pt_anypos_0.5_minbrl_MLE', '.pdf',sep='')
	pdf(file=outfile, w=8, h=30)
	plot.phylo(df.tpairs.mleph, show.node.label=1, cex=0.4, no.margin=1, label.offset=0.005, edge.width=0.5)
	dev.off()
	#			
	#	plot clinical data + dated phylogenies	
	#	
	df.all					<- merge(df.all, subset(clumsm.info, select=c(FASTASampleCode, cluster)), by='FASTASampleCode', all.x=1)
	outfile					<- paste(outdir,'/',infile, '_', clu.infilexml.template, '_', clu.infilexml.opt, '_', gsub('/',':',insignat), '_', 'pt_anypos_0.5_minbrl', '.pdf',sep='')
	project.athena.Fisheretal.plot.selected.transmitters(df.all, df.immu, df.viro, df.treatment, df.tpairs, cluphy, cluphy.info, cluphy.subtrees, cluphy.map.nodectime, outfile)
	
	
}
######################################################################################
project.athena.Fisheretal.similar<- function()
{
	require(data.table)
	require(ape)
	stop()
	indir					<- paste(DATA,"tmp",sep='/')		
	indircov				<- paste(DATA,"derived",sep='/')
	outdir					<- paste(DATA,"fisheretal",sep='/')
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"
	if(0)
	{
		infile					<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs100",sep="_")
		insignat				<- "Thu_Aug_01_17/05/23_2013"	
		clu.indir				<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/beast/beast2_140201'
		clu.infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_seroneg-130"
		clu.insignat			<- "Tue_Aug_26_09:13:47_2013"
		clu.infilexml.template	<- "sasky_sdr06"
		clu.infilexml.opt		<- "rsu815"				
	}
	if(1)
	{
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"				
		clu.indir				<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/sasky_sdr06_-DR-RC-SH+LANL_alrh160'
		clu.infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		clu.insignat			<- "Wed_Dec_18_11:37:00_2013"
		clu.infilexml.opt		<- "alrh160"
		clu.infilexml.template	<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=2',sep='_')
	}	
	#
	#	select infected individuals and return in df.select
	#
	tmp				<- project.athena.Fisheretal.select.denominator(indir, infile, insignat, indircov, infilecov, infiletree, adjust.AcuteByNegT=0.75, adjust.NegT4Acute=1, adjust.AcuteSelect=c('Yes','Maybe'))	
	df.denom		<- tmp$df.select
	clumsm.subtrees	<- tmp$clumsm.subtrees
	clumsm.info		<- tmp$clumsm.info
	clumsm.ph		<- tmp$clumsm.ph
	setkey(clumsm.info, cluster)
	if(0)
	{
		files		<- list.files(clu.indir)
		if(!length(files))	stop('no input files matching criteria')
		files		<- files[ sapply(files, function(x) grepl(clu.infile, x, fixed=1) & grepl(gsub('/',':',clu.insignat), x, fixed=1) & grepl(clu.infilexml.opt, x, fixed=1) & grepl(clu.infilexml.template,x, fixed=1) & grepl('_cluposterior_[0-9]+',x) & grepl('R$',x) ) ]		
		if(!length(files))	stop('no input files matching criteria')
		tmp			<- regmatches( files, regexpr('_cluposterior_[0-9]+',files)) 
		cluster		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
		file.info	<- data.table(file=files, cluster=cluster)
		setkey(file.info, cluster)
		
		clu.missing	<- sort( setdiff( clumsm.info[,unique(cluster)], file.info[,cluster] ) )
		clu.missing	<- unique( subset( clumsm.info, cluster%in%clu.missing, c(cluster, clu.npat, clu.ntip) ) )
		setkey(clu.missing, cluster)		
	}
	#
	#	select potential transmitters
	#
	df.tpairs		<- project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos(clumsm.info, df.denom, any.pos.grace.yr= 2, select.if.transmitter.seq.unique=FALSE)	
	#df.tpairs		<- project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos.MinBrlMLETree(df.denom, clumsm.subtrees, any.pos.grace.yr= 2)
	#	plot number of potential transmitters
	if(0)
	{
		tmp				<- merge( subset(df.tpairs, select=Patient), subset(df.denom, select=c(Patient, AnyPos_T1)), by='Patient' )
		file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'nrecentlyinfected', '.pdf',sep='')
		pdf(file=file, w=5, h=5)
		par(mar=c(3,5,0.5,0.5))
		barplot( table( tmp[, round(AnyPos_T1)] ), ylab="# recently infected\n with unique potential transmitter" )
		dev.off()	
	}	
	#
	#	compute Y score: raw branch length between pot transmitter and infected
	#
	file				<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'nbrlraw',sep='')	
	tmp					<- project.athena.Fisheretal.Y.rawbrl(df.tpairs, indir, insignat, indircov, infilecov, infiletree, save.file=paste(file, '.R', sep=''), plot.file=paste(file, '.pdf', sep=''))
	Y.rawbrl			<- tmp$tpairs
	Y.rawbrl.linked		<- tmp$linked
	Y.rawbrl.unlinked	<- tmp$unlinked
	#
	#	compute Y score: [0,1]: branch length weight between pot transmitter and infected
	#	
	file				<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl.pdf',sep='')
	Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked	, Y.rawbrl.unlinked, linked.x= 0, plot.file=file)
	#
	#	compute Y score: [0,1]: prob that viral lineage coalesces within seroconversion window and [0,1]: infection at time t happens after coalescence within potential transmitter
	#
	tmp						<- project.athena.Fisheretal.get.data.for.selection(df.tpairs, clu.indir, clu.infile, clu.insignat, clu.infilexml.opt, clu.infilexml.template, indircov, infilecov)	
	cluphy.map.nodectime	<- tmp$clu$cluphy.map.nodectime
	cluphy.subtrees			<- tmp$clu$cluphy.subtrees
	cluphy.info				<- tmp$clu$cluphy.info
	cluphy					<- tmp$clu$cluphy	
	df.viro					<- tmp$df.viro
	df.immu					<- tmp$df.immu
	df.treatment			<- tmp$df.treatment
	#
	#	compute X: follow up
	X.fwup					<- project.athena.Fisheretal.X.followup(df.tpairs, df.immu, t.period= 0.25, t.endctime=2013.0)
	#tmp					<- project.athena.Fisheretal.X.followup.compareCD4toVL(df.tpairs, df.immu, clumsm.info)
	#	compute X: calendar time period
	X.tperiod				<- project.athena.Fisheretal.X.calendarperiod(df.tpairs, clumsm.info, t.period= 0.25, c.nperiod= 4)
	#	compute X: RegionHospital and Exposure group
	X.Trm					<- project.athena.Fisheretal.X.Trm.Region(df.tpairs, clumsm.info)
	#	compute X: viral load of potential transmitter  	
	X.viro					<- project.athena.Fisheretal.X.viro(df.tpairs, df.viro, t.period=0.25, lRNA.cutoff=NA)
	
	
	t.period=0.25
	t.endctime= 2013.
	#	prepare incare timeline for potential transmitters
	incare		<- merge( data.table(Patient=df.tpairs[, unique(t.Patient)]), unique( subset(clumsm.info, select=c(Patient, AnyPos_T1, DateDied)) ), by='Patient' )
	set(incare, NULL, 'AnyPos_T1', incare[, floor(AnyPos_T1) + floor( (AnyPos_T1%%1)*100 %/% (t.period*100) ) * t.period] )
	set(incare, incare[,which(is.na(DateDied))], 'DateDied', t.endctime)
	set(incare, NULL, 'DateDied', incare[, floor(DateDied) + floor( (DateDied%%1)*100 %/% (t.period*100) ) * t.period] )	
	incare.t	<- incare[, list(t= seq(AnyPos_T1, DateDied, by=t.period)),by='Patient']	
	#	prepare treatment variables for potential transmitters
	treat		<- subset(df.treatment, select=c(Patient, AnyT_T1, StartTime, StopTime, TrI)) 
	treat		<- merge( data.table(Patient=df.tpairs[, unique(t.Patient)]), treat, by='Patient' )
	set(treat, NULL, 'AnyT_T1', hivc.db.Date2numeric(treat[,AnyT_T1]))
	set(treat, NULL, 'StopTime', hivc.db.Date2numeric(treat[,StopTime]))
	set(treat, NULL, 'StartTime', hivc.db.Date2numeric(treat[,StartTime]))
	#	set ever on ART per period t
	incare.t	<- merge(unique(subset(treat, select=c(Patient, AnyT_T1))), incare.t, by='Patient' )
	incare.t[, stage:='Diag']
	set(incare.t, incare.t[, which(AnyT_T1<=t+t.period/2)], 'stage', 'ART.started')
	#	identify periods t when treatment interrupted
	treat		<- subset(treat, TrI=='Yes', select=c(Patient, StartTime, StopTime))
	set(treat, NULL, 'StartTime', treat[, floor(StartTime) + floor( (StartTime%%1)*100 %/% (t.period*100) ) * t.period] )
	set(treat, NULL, 'StopTime', treat[, floor(StopTime) + floor( (StopTime%%1)*100 %/% (t.period*100) ) * t.period] )
	#	identify distinct periods t when treatment interrupted
	treat		<- treat[, {
								if(length(StartTime)>1)
								{
									tmp			<- StopTime[-length(StopTime)]+t.period < StartTime[-1]
									StartTime	<- StartTime[c(TRUE,tmp)]
									StopTime	<- StopTime[c(tmp,TRUE)]					
								}
								else
								{
									StartTime	<- StartTime
									StopTime	<- StopTime
								}
								list(StartTime=StartTime, StopTime=StopTime)
							} ,by='Patient']
	treat		<- treat[, list(t= seq(from=StartTime, to=StopTime, by=t.period), ART.interrupted='Yes'), by=c('Patient', 'StartTime')]
	#	add ART.interrupted time periods to incare.t
	incare.t	<- merge( subset(incare.t, select=c(Patient, t, stage)), subset(treat, select=c(Patient, t, ART.interrupted)), by=c('Patient','t'), all.x=1)
	set(incare.t, incare.t[,which(is.na(ART.interrupted))], 'ART.interrupted', 'No')
	#	add viro time periods to incare.t
	setnames(incare.t, 'Patient','t.Patient')

	incare.t	<- merge(incare.t, X.viro, by=c('t.Patient','t'), all.x=1)
	#	TODO lost M42459 ? --> should be because clusters are missing
	#	TODO add CD4 by time period for completeness, so that we can easily make any changes as needed later



	
	treat[, length(StartTime), by='Patient']
	tmp						<- project.athena.Fisheretal.Y.coal.and.inf(df.tpairs, cluphy, cluphy.info, cluphy.map.nodectime, t.period= 0.25 )
	df.tpairs				<- tmp$df.tpairs
	Y.inf					<- tmp$inf
	Y.coal					<- tmp$coal	
	df.tpairs.excluded		<- tmp$excluded
		
	
	
	
	
	#
	#	get data for selection
	#
	df.tpairs<- subset( df.tpairs, cluster%in%c(1205, 1510) )
	tmp						<- project.athena.Fisheretal.get.data.for.selection(df.tpairs, clu.indir, clu.infile, clu.insignat, clu.infilexml.opt, clu.infilexml.template, indircov, infilecov)	
	df.all					<- tmp$df.all
	df.viro					<- tmp$df.viro
	df.immu					<- tmp$df.immu
	df.treatment			<- tmp$df.treatment
	cluphy					<- tmp$clu$cluphy
	cluphy.subtrees			<- tmp$clu$cluphy.subtrees
	cluphy.info				<- tmp$clu$cluphy.info
	cluphy.map.nodectime	<- tmp$clu$cluphy.map.nodectime	
	#			
	#	plot	
	#	
	outfile					<- paste(indir,'/',infile, '_', clu.infilexml.template, '_', clu.infilexml.opt, '_', gsub('/',':',insignat), '_', 'pt_anypos_0.5_minbrl_coalb4', '.pdf',sep='')
	project.athena.Fisheretal.plot.selected.transmitters(df.all, df.immu, df.viro, df.treatment, df.tpairs, cluphy, cluphy.info, cluphy.subtrees, cluphy.map.nodectime, outfile)

	
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
	
project.hivc.examl<- function(dir.name= DATA)
{
	require(ape)
	
	indir		<- paste(DATA,"derived",sep='/')
	outdir		<- paste(DATA,"tmp",sep='/')
	signat.out	<- signat.in	<- "Wed_May__1_17/08/15_2013"
	verbose		<- resume		<- 1
	infile		<- "ATHENA_2013_03_FirstAliSequences_PROTRT"
	file		<- paste(indir,'/',infile,"_",gsub('/',':',signat.in),".R",sep='')
	if(verbose) cat(paste("\nload ",file))
	load(file)
	str(seq.PROT.RT)
	
	#debug	
	
	
	cat(cmd)
	stop()
	#create ExaML binary file from phylip
	#create Parsimonator starting tree
	#run ExaML starting tree
	#delete phylip file
	
}

project.hivc.clustalo<- function(dir.name= DATA, min.seq.len=21)
{			
	if(0)
	{
		#generate fasta files per gene
		verbose<- 1
		analyze<- 1
		
		file<- paste(dir.name,"tmp/ATHENA_2013_03_Sequences.R",sep='/')
		load(file)		
		lapply(seq_along(df),function(i)
				{
					df.gene<- df[[i]]
					#print(length(df.gene[,"SampleCode"]))
					#print(length(unique(df.gene[,"SampleCode"])))
					#stop()
					df.gene[,"SampleCode"]	<- paste('>',df.gene[,"SampleCode"],sep='')
					seq.len		<- nchar(df.gene[,"Sequence"])	
					if(analyze)
					{
						cat(paste("\nprocess",names(df)[i]))
						cat(paste("\nnumber of sequences:",nrow(df.gene),"\nnumber of unique SampleCodes:", length(unique(df.gene[,"SampleCode"])),"\n") )
						hist(seq.len, breaks=100)	
						if(nrow(df.gene)!=length(unique(df.gene[,"SampleCode"])))
							o<- sapply(unique(df.gene[,"SampleCode"]),function(x)
									{
										tmp<- which(df.gene[,"SampleCode"]==x)
										if(length(tmp)>1)
										{ 
											cat("\n found duplicates\n")
											print(df.gene[tmp,])	
										}
									})
						nok.idx	<- which(seq.len<min.seq.len)
						if(verbose)	cat(paste("discarding sequences with insufficient length, n=",length(nok.idx),'\n'))
						stop("something is odd here - code acc deleted ?")
					}
					tmp			<- as.vector(t(df.gene))
					tmp			<- paste(tmp,collapse='\n')					
					file		<- paste(dir.name,"/derived/ATHENA_2013_03_Sequences_",names(df)[i],".fa",sep='')
					cat(paste("\nwrite file to",file))
					cat(tmp,file=file)	
				})			
	}
	if(0)	#edit multiple sequence alignment: remove gaps
	{
		verbose<- 1
		require(ape)
		signat.in	<- "Wed_Apr_24_09/02/03_2013"
		signat.out	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
		indir		<- paste(dir.name,"tmp",sep='/')
		pattern 	<- gsub('/',':',paste(signat.in,".clustalo$",sep=''))
		files		<- list.files(path=indir, pattern=pattern, full.names=0)
		lapply(files,function(x)
				{
					cat(paste("\nprocess",x))
					tmp		<- read.dna( paste(indir,x,sep='/'), format="fa", as.matrix=1 )
					tmp		<- seq.rmgaps(tmp, verbose=verbose)
					file	<- paste(dir.name, "tmp", gsub(pattern, paste(signat.out,"clustalo",sep='.'),x),sep='/')
					cat(paste("\nwrite to",file))
					write.dna(tmp, file, format= "fasta" )
				})		
	}	
	if(0)
	{
		#generate initial clustalo command		 
		indir	<- paste(dir.name,"derived",sep='/')
		outdir	<- paste(dir.name,"tmp",sep='/')
		infiles	<- list.files(path=indir, pattern=".fa$", full.names=0)		
		signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
		cmd		<- hiv.cmd.clustalo(indir, infiles, signat=signat, outdir=outdir)
		cmd		<- hiv.cmd.hpcwrapper(cmd, hpc.q="pqeph")
		
		cat(cmd[[1]])
		#lapply(cmd, cat)
		stop()
		lapply(cmd,function(x)
				{
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("clustalo",signat,sep='.')					
					hiv.cmd.hpccaller(outdir, outfile, x)			
				})			
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

