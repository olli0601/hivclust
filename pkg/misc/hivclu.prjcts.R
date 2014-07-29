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
project.Gates.RootSeqSim<- function()
{
	DATA		<<- "/work/or105/Gates_2014"
	indir		<- paste(DATA,'methods_comparison_rootseqsim/140727',sep='/')
	infile		<- 'ALLv02.n100.rlx.gmrf' 
	insignat	<- 'Sun_Jul_27_09-00-00_2014'
	cmd			<- hivc.cmd.beast.runxml(indir, infile, insignat, prog.beast=PR.BEAST, prog.beastmcc=PR.BEASTMCC, beastmcc.burnin=500, beastmcc.heights="median", hpc.tmpdir.prefix="beast", hpc.ncpu=1)
	cat(cmd)
	
	cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=1, hpc.walltime=91, hpc.mem="1800mb")
	outdir		<- indir
	outfile		<- paste("b2m.",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='')					
	hivc.cmd.hpccaller(outdir, outfile, cmd)
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
project.hivc.Excel2dataframe.AllPatientCovariates<- function(dir.name= DATA, verbose=1, resume=0)
{	
	require(data.table)		
	
	#input files generated with "project.hivc.Excel2dataframe"
	resume			<- 0
	verbose			<- 1
	if(1)
	{
		file.seq		<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.R",sep='/')
		file.patient	<- paste(dir.name,"derived/ATHENA_2013_03_Patient.R",sep='/')
		file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
		file.immu		<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
		file.treatment	<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
		file.out		<- paste(dir.name,"derived/ATHENA_2013_03_AllSeqPatientCovariates.R",sep='/')
	}
	if(0)
	{
		file.seq		<- paste(dir.name,"derived/ATHENA_2013_03_Sequences_AllMSM.R",sep='/')
		file.patient	<- paste(dir.name,"derived/ATHENA_2013_03_Patient_AllMSM.R",sep='/')
		file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro_AllMSM.R",sep='/')
		file.immu		<- paste(dir.name,"derived/ATHENA_2013_03_Immu_AllMSM.R",sep='/')
		file.treatment	<- paste(dir.name,"derived/ATHENA_2013_03_Regimens_AllMSM.R",sep='/')
		file.out		<- paste(dir.name,"derived/ATHENA_2013_03_AllSeqPatientCovariates_AllMSM.R",sep='/')
	}
	
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
		#set(df, NULL, "FASTASampleCode", factor(df[,FASTASampleCode]))
		#set(df, NULL, "Patient", factor(df[,Patient]))		
		df.all	<- subset(df, select=c(FASTASampleCode,Patient,PosSeqT))
		if(verbose)		cat(paste("\nnumber of sequences found, n=", nrow(df.all)))
		#
		# add Patient data
		#
		if(verbose)		cat(paste("\nadding patient data"))
		load(file.patient)		
		df.all	<- merge(df.all, df, all.x=1, all.y=1, by="Patient")
		setnames(df.all, c("MyDateNeg1","MyDatePos1"), c("NegT","PosT"))
		#	reset PosT_Acc=='No' conservatively to end of year / month
		df.all	<- hivc.db.reset.inaccuratePosT(df.all, nacc.dy.dy= 30, nacc.mody.mo= 11, nacc.mody.dy= 31, verbose=1)
		#	reset NegT_Acc=='No' conservatively to start of year / month
		df.all	<- hivc.db.reset.inaccurateNegT(df.all)
		#
		#	check for clashes in NegT and PosT
		# 	add viro data to check PosT and NegT
		#
		load(file.viro)	
		tmp		<- subset(df, select=c(Patient, PoslRNA_T1, lRNA_T1))
		setkey(tmp, Patient)
		df.all	<- merge(df.all, unique(tmp), by='Patient', all.x=1)
		load(file.treatment)
		tmp		<- subset(df, select=c(Patient, AnyT_T1))
		setkey(tmp, Patient)
		df.all	<- merge(df.all, unique(tmp), by='Patient', all.x=1)
		tmp		<- subset(df.all, !is.na(PosT) & !is.na(NegT) & PosT<NegT)
		if(verbose)		cat(paste("\nchecking for clashes in NegT and PosT"))
		if(verbose)		print(tmp)		
		#	wrong NegT	M12982 M31263 M37982 M40070 M41206	M11399 M12702	
		#	wrong PosT M41654
		tmp		<- which(df.all[,FASTASampleCode=="M4165429052012"])
		set(df.all, tmp, "PosT", df.all[tmp, PosSeqT]) 		
		tmp		<- which(df.all[,Patient%in%c("M12982","M31263","M37982","M40070","M41206","M11399","M12702")])
		set(df.all, tmp, "NegT", NA_integer_) 	
		#	set missing PosT		
		tmp		<- df.all[, which(is.na(PosT) & !is.na(AnyT_T1))]
		if(verbose)		cat(paste("\nSet missing PosT to  AnyT_T1, n=", length(tmp)))
		set(df.all,tmp,'PosT',df.all[tmp, AnyT_T1])
		tmp		<- df.all[, which(is.na(PosT) & !is.na(PoslRNA_T1))]
		if(verbose)		cat(paste("\nSet missing PosT to  PoslRNA_T1, n=", length(tmp)))
		set(df.all,tmp,'PosT',df.all[tmp, PoslRNA_T1])
		#	remove rest
		if(verbose)		cat(paste("\nREMOVE"))
		print(subset(df.all, is.na(PosT)))
		df.all	<- subset(df.all, !is.na(PosT))
		df.all	<- df.all[, setdiff( colnames(df.all), c("PoslRNA_T1","lRNA_T1")), with=FALSE]
		# 	add preliminary AnyPos_T1	-- since PosT conservative we are good to set irrespective of PosT_Acc	
		df.all[, AnyPos_T1:=PosSeqT]
		tmp		<- which( df.all[, !is.na(PosT) & ( is.na(PosSeqT) | PosT<PosSeqT ) ] )
		if(verbose)		cat(paste("\nbuilding prel AnyPos_T1. Number of seq with !is.na(PosT) & ( is.na(PosSeqT) | PosT<PosSeqT ), n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, PosT])	
		if(nrow(subset(df.all, is.na(AnyPos_T1))))	stop("unexpected NA in AnyPos_T1")
		#	check for invalid NegT and set to NA	-- we would only know that NegT is before AnyPosT and this is not helpful
		df.all	<- hivc.db.resetNegTbyAnyPosT(df.all)				
		#
		#	checking for SeqT before AnyPos_T1
		#
		df.all	<- merge(df.all, df.all[, list(PosSeqT.min=min(PosSeqT)),by='Patient'], by='Patient')
		tmp		<- df.all[, which(AnyPos_T1>PosSeqT.min)]
		if(verbose)		cat(paste("\nFound AnyPos_T1>PosSeqT.min, reset. n=", length(tmp)))
		set(df.all, tmp, 'AnyPos_T1', df.all[tmp, PosSeqT.min])
		#
		#	add first RNA Virology date
		#
		if(verbose)		cat(paste("\nadding virology data"))
		load(file.viro)
		df.all[, idx:=seq_len(nrow(df.all))]
		df.cross	<- merge( subset(df.all, select=c(idx,FASTASampleCode,Patient,AnyPos_T1,PosT,PosSeqT,AnyT_T1,NegT,NegT_Acc)), df, allow.cartesian=T, by="Patient" )
		#	checking manually AnyPos_T1>PosRNA & lRNA<3
		if(verbose)		cat(paste("\ncheck manually AnyPos_T1>PosRNA & lRNA<3 -- THIS IS ASSUMED OK"))
		tmp	<- subset(df.cross, AnyPos_T1>PosRNA & AnyT_T1>PosRNA  & lRNA<3)
		tmp[,diff:=tmp[, difftime(PosRNA,AnyPos_T1, units="weeks")]]
		print( subset(tmp, diff< -10) )		
		#		
		#	checking manually NegT>PosRNA
		#
		if(verbose)		cat(paste("\ncheck manually NegT_Acc=='Yes' & NegT>PosRNA -- THIS IS ASSUMED OK"))
		print( subset(df.cross, NegT_Acc=="Yes" & NegT>PosRNA) )
		if(verbose)		cat(paste("\ncheck manually NegT_Acc=='No' & NegT>PosRNA -- THIS IS ASSUMED OK"))
		print( subset(df.cross, NegT_Acc=="No" & NegT>PosRNA) )
		# 	compute lRNA_T1 and lRNA_TS
		tmp			<- hivc.db.getlRNA.T1andTS(df.cross, lRNA.bTS.quantile= 0.75, lRNA.aTS.quantile= 0.25, lRNAi.min= log10(1e4), verbose=1)
		df.all		<- merge(df.all, tmp, all.x=1, by="idx")
		#	reset NegT by lRNA_T1 -- there is an error in here
		df.all		<- hivc.db.resetNegTbyPoslRNA_T1(df.all)
		# 	reset preliminary AnyPos_T1
		tmp			<- df.all[, which(AnyPos_T1>PoslRNA_T1)]
		if(verbose)	cat(paste("\nnumber AnyPos_T1>PoslRNA_T1 --> RESETTING",length(tmp)))
		set(df.all,tmp,"AnyPos_T1",df.all[tmp,PoslRNA_T1])	
		if(verbose)	cat(paste("\nnumber of seq with NegT==AnyPos_T1 --> RESETTING NegT"))
		print(subset(df.all, NegT==AnyPos_T1, c(idx, Patient, NegT, NegT_Acc, PosT, PosT_Acc, AnyPos_T1, AnyT_T1, PoslRNA_T1, lRNA_T1)))
		tmp		<- df.all[, which(NegT==AnyPos_T1)]
		set(df.all, tmp, 'NegT', df.all[tmp, AnyPos_T1-3*30])
		set(df.all, tmp, 'NegT_Acc', 'No')
		df.all[,AnyT_T1:=NULL]
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
		df.cross	<- merge( subset(df.all, select=c(idx,FASTASampleCode,Patient,PosSeqT)), df, allow.cartesian=T, by="Patient" )
		tmp			<- hivc.db.getTrIMo(df.cross)
		df.all		<- merge(df.all, tmp, all.x=1, by="idx")		
		tmp			<- subset(df.all, AnyT_T1<AnyPos_T1)
		if(verbose)		cat(paste("\ncheck manually AnyT_T1<AnyPos_T1 -- THIS IS ASSUMED OK"))
		print(tmp, n=300)
		tmp			<- df.all[, which(AnyT_T1<AnyPos_T1)]
		set(df.all, tmp, 'AnyPos_T1', df.all[tmp, AnyT_T1])		
		#
		#	add CD4 count data
		#
		if(verbose)		cat(paste("\nadding CD4 data"))
		load(file.immu)
		if(length(which(df[,is.na(PosCD4_T1)])))	stop("unexpected NA in PosCD4_T1")
		if(length(which(df[, is.na(PosCD4)]))) stop("unexpected NA in PosCD4")
		if(length(which(df[, is.na(CD4)]))) stop("unexpected NA in CD4")
		df.cross	<- merge( subset(df.all, select=c(idx,FASTASampleCode,Patient,AnyPos_T1,PosT,PosSeqT,PoslRNA_T1,NegT,NegT_Acc)), subset(df, select=c(Patient,PosCD4,CD4)), allow.cartesian=T, by="Patient" )
		# delete entries where NegT>PosCD4
		df.cross	<- hivc.db.resetCD4byNegT(df.cross, with.NegT_Acc.No=1, verbose=1)
		# compute CD4_T1 and CD4_TS -- do not use on df[,CD4_T1] because some CD4 measurements might correspond to healthy patients
		tmp		<- hivc.db.getCD4.T1andTS(df.cross, CD4.HIVNeg.min= 500)
		df.all	<- merge(df.all, tmp, by="idx", all.x=1)
		# manually checked remaining PosCD4_T1 < AnyPos_T1
		df.all	<- hivc.db.reset.PosT1byCD4T1(df.all, verbose=1)
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
		
		df.all[, idx:=NULL]
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
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.csv",sep='/')
	file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Regimens_AllMSM.csv",sep='/')	
	file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro_AllMSM.R",sep='/')
	
	load(file.viro)
	df.viro		<- subset(df, select=c(Patient, PosRNA, lRNA))
	
	NA.time			<- c("01/01/1911","01/11/1911","11/11/1911")	
	MAX.time		<- c("")
	TR.notyet		<- "30/03/2013"
	TR.failure 		<- c(21, 31 ) 							#viro failure
	TR.immufailure	<- c(32, 35)							#either immu failure or new CDC-B/C event
	TR.toxicity  	<- c(24, 34) 							#toxicity
	TR.adherence	<- c(47)
	TR.patrel		<- c(23, 33, 42, 43)					#either patient s decision, desired pregnancy or pregnancy
	#read REGIMEN csv data file and preprocess	
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
	set(df, NULL, 'Patient', df[, as.character(Patient)])
	setnames(df, "T0","HAART_T1")
	#set(df,NULL,"Patient",factor(df[,Patient]))
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
	#	check NoDrug==0 entries
	if(verbose)	cat(paste("\nnumber of patients with is.na(NoDrug), n=",nrow(subset( df[, list(select=any(is.na(NoDrug))), by='Patient'], select))))
	#	remove patients for which all periods have NoDrug==0
	tmp						<- df[,  list(select=!all(NoDrug==0)) , by='Patient']
	if(verbose)	cat(paste("\nnumber of patients with all(NoDrug==0), n=",nrow(subset(tmp, !select))))
	df						<- merge( df, subset(tmp, select, Patient), by='Patient')
	#	sort by Patient and StopTime 
	setkey(df, Patient, StopTime)
	#	for each patient find first periods with NoDrug>0
	tmp						<- subset( df[, list( select= which(NoDrug>0)[1]>1 ), by='Patient'], select)[, Patient]
	if(verbose)	cat(paste("\nnumber of patients with ! NoDrug>0 as first StartTime for any patient, n=",length(tmp)))	
	tmp						<- df[, list( StopTime= StopTime[ seq.int( which(NoDrug>0)[1], length(NoDrug)) ] ), by='Patient']
	if(verbose)	cat(paste("\nnumber of entries with NoDrug>0 as first StartTime for any patient, n=",nrow(tmp)))
	df						<- merge(df, tmp, by=c('Patient','StopTime'))	
	tmp						<- merge( df, subset( df[,list(check= any( is.na(StartTime) & NoDrug==0) ), by='Patient'], check, Patient), by='Patient') 
	if(verbose)	cat(paste("\nnumber of entries with NoDrug==0 and is.na StartTime, n=",nrow(tmp)))
	#
	#	remove duplicate entries
	#
	setkey(df, Patient, StartTime, StopTime)
	tmp		<- which(duplicated(df))
	if(verbose)	cat(paste("\nnumber of duplicate entries with key  Patient, StartTime, StopTime, n=",nrow(tmp)))
	df		<- unique(df)
	#M15198 1998-05-28 1996-08-22 1998-02-24
	#M15198 1998-05-28 1996-08-22 1998-07-01
	#M17774 1997-03-30 2013-03-30 1996-08-26
	#M17774 1997-03-30 2013-03-30 1998-06-21
	tmp						<- which(df[, Patient=="M15198" & StartTime=="1998-07-01"])
	if(verbose)	cat(paste("\nrm entry 	M15198 1998-05-28 1996-08-22 1998-07-01"))
	df						<- subset(df, !(Patient=="M15198" & StartTime=="1998-07-01"))
	#
	#	fix StartTime>StopTime
	#	
	tmp		<- which(df[,StartTime>StopTime])	
	if(verbose)	cat(paste("\nnumber of entries with StartTime>StopTime, n=",length(tmp)))
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
	#	see if we can merge consecutive NoDrug==0 periods	
	df.nodrug		<- subset(df, NoDrug==0)
	if(verbose)	cat(paste("\nnumber of entries with NoDrug==0, n=",nrow(tmp)))
	tmp				<- df.nodrug[, {
										if(length(StartTime)>1)
										{
											tmp			<- StopTime[-length(StopTime)] != StartTime[-1]
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
	if(verbose)	cat(paste("\nnumber of merged entries with NoDrug==0, n=",nrow(tmp)))						
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
	df.TrCh.noreason		<- which( df[, is.na(Reason1)&is.na(Reason2)&is.na(Reason3)&is.na(Reason4)&is.na(Reason5)&is.na(Reason6)&(NoDrug!=0) ] )
	#	TR.addition
	tmp	<- na.omit(c( 	df[df.TrCh.noreason, AZT]<=df[df.TrCh.noreason+1, AZT] &
				df[df.TrCh.noreason, ddI]<=df[df.TrCh.noreason+1, ddI] &
				df[df.TrCh.noreason, ddC]<=df[df.TrCh.noreason+1, ddC] &
				df[df.TrCh.noreason, d4T]<=df[df.TrCh.noreason+1, d4T] &
				df[df.TrCh.noreason, dTC]<=df[df.TrCh.noreason+1, dTC] &
				df[df.TrCh.noreason, IDV]<=df[df.TrCh.noreason+1, IDV] & 
				df[df.TrCh.noreason, SAQ]<=df[df.TrCh.noreason+1, SAQ] &
				df[df.TrCh.noreason, RTV]<=df[df.TrCh.noreason+1, RTV] &
				df[df.TrCh.noreason, NVP]<=df[df.TrCh.noreason+1, NVP] &
				df[df.TrCh.noreason, EFV]<=df[df.TrCh.noreason+1, EFV] &
				df[df.TrCh.noreason, ETR]<=df[df.TrCh.noreason+1, ETR] &
				df[df.TrCh.noreason, RPV]<=df[df.TrCh.noreason+1, RPV] &
				df[df.TrCh.noreason, NFV]<=df[df.TrCh.noreason+1, NFV] &
				df[df.TrCh.noreason, ABC]<=df[df.TrCh.noreason+1, ABC] &
				df[df.TrCh.noreason, LPV]<=df[df.TrCh.noreason+1, LPV] &
				df[df.TrCh.noreason, ATV]<=df[df.TrCh.noreason+1, ATV] &
				df[df.TrCh.noreason, ENF]<=df[df.TrCh.noreason+1, ENF] &
				df[df.TrCh.noreason, TDF]<=df[df.TrCh.noreason+1, TDF] &
				df[df.TrCh.noreason, TPV]<=df[df.TrCh.noreason+1, TPV] &
				df[df.TrCh.noreason, FPV]<=df[df.TrCh.noreason+1, FPV] &
				df[df.TrCh.noreason, FTC]<=df[df.TrCh.noreason+1, FTC] &
				df[df.TrCh.noreason, DRV]<=df[df.TrCh.noreason+1, DRV] &
				df[df.TrCh.noreason, MVC]<=df[df.TrCh.noreason+1, MVC] &
				df[df.TrCh.noreason, RAL]<=df[df.TrCh.noreason+1, RAL]		))	
	#tmp						<- df[df.TrCh.noreason, NoDrug+1]==df[df.TrCh.noreason+1, NoDrug]	
	df[, TrCh.addition:=0]
	set(df, df.TrCh.noreason[tmp], 'TrCh.addition', 1)
	df[, TrCh.addition:= factor(TrCh.addition,levels=c(0,1),labels=c("No","Yes"))]
	df.TrCh.noreason		<- df.TrCh.noreason[!tmp]
	#	TR.failure
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.failure | Reason2%in%TR.failure | Reason3%in%TR.failure | Reason4%in%TR.failure | Reason5%in%TR.failure | Reason6%in%TR.failure ]) ]<- 1
	df[, TrCh.failure:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#	TR.immufailure
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.immufailure | Reason2%in%TR.immufailure | Reason3%in%TR.immufailure | Reason4%in%TR.immufailure | Reason5%in%TR.immufailure | Reason6%in%TR.immufailure ]) ]<- 1
	df[, TrCh.immufailure:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]	
	#	TR.toxicity
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.toxicity | Reason2%in%TR.toxicity | Reason3%in%TR.toxicity | Reason4%in%TR.toxicity | Reason5%in%TR.toxicity | Reason6%in%TR.toxicity ]) ]<- 1
	df[, TrCh.toxicity:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]		
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
	#	add largest lRNA at least 3 months after therapy start.
	tmp					<- subset(df, select=c(Patient, AnyT_T1)) 
	setkey(tmp, Patient)
	tmp					<- merge(unique(tmp), df.viro, by='Patient')	
	tmp					<- merge(subset(df, select=c(Patient, AnyT_T1, StartTime, StopTime)), subset(tmp, difftime(PosRNA, AnyT_T1, units='days')>90,select=c(Patient, PosRNA, lRNA)), by='Patient', allow.cartesian=TRUE, all.x=TRUE)
	tmp					<- subset(tmp, is.na(PosRNA) | (StartTime<=PosRNA & PosRNA<StopTime))
	tmp					<- tmp[, {
										z<- which.max(lRNA)
										list(AnyT_T1=AnyT_T1[z], PosRNA.mx= PosRNA[z], lRNA.mx=lRNA[z])
									}, by=c('Patient','StartTime','StopTime')]
	df					<- merge( df, subset(tmp, select=c(Patient, StartTime, StopTime, PosRNA.mx, lRNA.mx)), by=c('Patient','StartTime','StopTime'), all.x=TRUE)
	#	set failure by viral load during treatment period
	df[, TrVL.failure:=factor( df[, lRNA.mx>3], levels=c(FALSE, TRUE), labels=c('No','Yes'))]	
	#subset(df, TrVL.failure=='No' & TrCh.failure=='Yes')	viro failure can be indicated for VL < 1e3
	#subset(df, TrVL.failure=='No' & TrCh.failure=='Yes' & lRNA.mx<log10(200), select=c(Patient, StartTime, StopTime, AnyT_T1, PosRNA.mx, lRNA.mx))	#some are really suspicious!
	#
	file		<- paste(substr(file, 1, nchar(file)-3),'R',sep='')
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)
}
######################################################################################
project.hivc.Excel2dataframe.Regimen.CheckStartVL<- function(dir.name= DATA, verbose=1)
{
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	load(file)
	df.treatment<- df
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	load(file)
	df.v		<- df
	
	check<- subset(df.v, Patient%in% c("M35938","M36809","M34392","M31657","M32854","M32608","M31261","M31145","M31104","M31101","M30994","M28495","M27763","M17749","M17487","M17218","M16622") )	
	check	<- merge(check, unique(subset(df.treatment, select=c(Patient, AnyT_T1))), by='Patient')
	print(check,n=550)
	#M36809
	
	#	check if VL decreases just before treatment start
	tmp			<- subset(df.treatment, !is.na(AnyT_T1), select=c(Patient, AnyT_T1, StartTime, StopTime, AnyT_T1_Acc, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel))
	
	
	if(verbose)	cat(paste("\ncheck VL during ART for patients, n=",nrow(tmp)))
	check		<- merge( subset(df.v, select=c(Patient, PosRNA, RNA, lRNA)), tmp, by='Patient', allow.cartesian=TRUE )
	check		<- check[, {
								z<- StartTime<=PosRNA & PosRNA<=StopTime
								lapply(.SD, '[', z)
							},  by=c('Patient','PosRNA')]	
	setkey(check, Patient, PosRNA)
	check		<- subset(check,  difftime(PosRNA, AnyT_T1, units='days')<= 0.25*365)
	
	
}
######################################################################################
project.hivc.Excel2dataframe.Regimen.CheckOnVL<- function(dir.name= DATA, verbose=1)
{
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	load(file)
	df.treatment<- df
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	load(file)
	df.v		<- df
	#	check if VL decreases just before treatment start
	tmp			<- subset(df.treatment, !is.na(AnyT_T1), select=c(Patient, AnyT_T1, StartTime, StopTime, AnyT_T1_Acc, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel))
	
	
	if(verbose)	cat(paste("\ncheck VL during ART for patients, n=",nrow(tmp)))
	check		<- merge( subset(df.v, select=c(Patient, PosRNA, RNA, lRNA)), tmp, by='Patient', allow.cartesian=TRUE )
	check		<- check[, {
								z<- StartTime<=PosRNA & PosRNA<=StopTime
								lapply(.SD, '[', z)
							},  by=c('Patient','PosRNA')]	
	setkey(check, Patient, PosRNA)	
	check		<- subset(check,  difftime(PosRNA, AnyT_T1, units='days')>= 0.5)
	if(verbose)	cat(paste("\nPatients with VL after ART start, n=",check[,length(unique(Patient))]))	
	#	select patients for which lRNA stays above 4.5 
	check		<- merge(check, subset( check[, list(select= any(lRNA>4.5)), by='Patient'], select ), by='Patient' )
	if(verbose)	cat(paste("\nPatients with high VL after ART start, n=",check[,length(unique(Patient))]))
	
	check.run	<- check[, list( VLhighrun=max(table(cumsum(c(0,as.numeric(diff(which(lRNA>4.5))!=1))))) ), by='Patient']
	check1		<- merge(check, subset(check.run, VLhighrun>20, Patient), by='Patient')
	check1[, length(unique(Patient))]
	print(check1, n=400)
	
	check2		<- merge(check, subset(check.run, VLhighrun>10, Patient), by='Patient')
	check2[, length(unique(Patient))]	
	check2		<- check2[,  {
				tmp<- lRNA>4.5
				list( 	TrI.na=length(which(is.na(TrI))), TrI.y=length(which(TrI=='Yes')), TrF.na=length(which(is.na(TrCh.failure))), TrF.y=length(which(TrCh.failure=='Yes')),
						TrA.na=length(which(is.na(TrCh.adherence))), TrA.y=length(which(TrCh.adherence=='Yes')),
						TrP.na=length(which(is.na(TrCh.patrel))), TrP.y=length(which(TrCh.patrel=='Yes')), Run=length(tmp) )
			},by='Patient']
	check2		<- check2[, lapply(.SD,'/',Run), by='Patient']
	check2[, lapply(.SD,mean), .SDcols=c('TrI.na','TrI.y','TrF.na','TrF.y','TrA.na','TrA.y','TrP.na','TrP.y')]
	
	check3		<- merge(check, subset(check.run, VLhighrun<10 & VLhighrun>4, Patient), by='Patient')
	check3[, length(unique(Patient))]	
	check3		<- check3[,  {
				tmp<- lRNA>4.5
				list( 	TrI.na=length(which(is.na(TrI))), TrI.y=length(which(TrI=='Yes')), TrF.na=length(which(is.na(TrCh.failure))), TrF.y=length(which(TrCh.failure=='Yes')),
						TrA.na=length(which(is.na(TrCh.adherence))), TrA.y=length(which(TrCh.adherence=='Yes')),
						TrP.na=length(which(is.na(TrCh.patrel))), TrP.y=length(which(TrCh.patrel=='Yes')), Run=length(tmp) )
			},by='Patient']
	check3		<- check3[, lapply(.SD,'/',Run), by='Patient']
	check3[, lapply(.SD,mean), .SDcols=c('TrI.na','TrI.y','TrF.na','TrF.y','TrA.na','TrA.y','TrP.na','TrP.y')]
	

	hist( subset( check, !is.na(TrCh.failure) & TrCh.failure=='Yes' )[, lRNA], xlab= 'log10 VL for ART failure', breaks=c(0,3,4,5,7) )
	hist( subset( check, !is.na(TrI) & TrI=='Yes' )[, lRNA], xlab= 'log10 VL for ART interruption', breaks=c(0,3,4,5,7) )	
	hist( subset( check, !is.na(TrCh.adherence) & TrCh.adherence=='Yes' )[, lRNA], xlab= 'log10 VL for ART adherence', breaks=c(0,3,4,5,7) )	
	hist( subset( check, !is.na(TrCh.patrel) & TrCh.patrel=='Yes' )[, lRNA], xlab= 'log10 VL for ART patient request', breaks=c(0,3,4,5,7) )
	
	hist( subset( check, !is.na(TrCh.failure) & TrCh.failure=='Yes' & lRNA>2.8)[, lRNA], xlab= 'log10 VL for ART failure', breaks=100 )
	tmp			<- c( subset( check, !is.na(TrI) & TrI=='Yes' & lRNA>2.8)[, lRNA], subset( check, !is.na(TrCh.adherence) & TrCh.adherence=='Yes' & lRNA>2.8)[, lRNA], subset( check, !is.na(TrCh.patrel) & TrCh.patrel=='Yes' & lRNA>2.8 )[, lRNA] )
	hist( tmp, xlab= 'log10 VL for ART I/A/P', breaks=100 )
	
	check[, lRNAc:= cut(check[,lRNA], breaks=c(0,3,4,5,20), labels=c('S','3-4','4-5','>5'), right=FALSE)]
	table( subset( check, !is.na(TrI) & TrI=='Yes' )[, lRNAc], xlab= 'log10 VL for ART interruption' )
}
######################################################################################
project.hivc.Excel2dataframe.Regimen.CheckARTStartDate<- function(dir.name= DATA, verbose=1)
{
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	load(file)
	df.treatment<- df
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	load(file)
	df.v		<- df
	#	check if VL decreases just before treatment start
	tmp			<- unique(subset(df.treatment, !is.na(AnyT_T1), select=c(Patient, AnyT_T1, AnyT_T1_Acc)))
	if(verbose)	cat(paste("\ncheck ART start date for patients, n=",nrow(tmp)))
	check		<- merge( subset(df.v, select=c(Patient, PosRNA, RNA, lRNA)), tmp, by='Patient' )
	setkey(check, Patient, PosRNA)
	check		<- subset(check,  difftime(AnyT_T1, PosRNA, units='days')<= 1.25*365 & difftime(PosRNA, AnyT_T1, units='days')<= 0.25*365)	
	check		<- merge(check, subset( check[, list(select=difftime(AnyT_T1[1], PosRNA[1], units='days')>0), by='Patient'], select, Patient ), by='Patient')	
	if(verbose)	cat(paste("\nPatients with PosRNA preceeding AnyT_T1, n=",check[,length(unique(Patient))]))
	tmp			<- check[,  {
								z<- tail(which( difftime(AnyT_T1[1], PosRNA, units='days')>=0 ),1)
								list(lRNA.drop.before.ART=lRNA[1] - lRNA[z], lRNA.before.ART=lRNA[z], lRNA1=lRNA[1], drop.day=PosRNA[z], drop.days2ART= difftime(AnyT_T1[1], PosRNA[z], units='days'))					
							}  ,by='Patient']
	check		<- merge(check, tmp, by='Patient')
					
	check1		<- subset(check, lRNA.drop.before.ART>0.5 & lRNA.before.ART<4 & AnyT_T1_Acc=='NAccD')
	if(verbose)	cat(paste("\nPatients lRNA.drop.before.ART>0.5 & AnyT_T1_Acc=='NAccD', n=",check1[,length(unique(Patient))]))
	print(check1[,unique(Patient)])
	#print(check1,n=500)
	#	reset AnyT_T1 to drop.day	
	reset		<- unique( subset(check1, select=c(Patient, AnyT_T1, drop.day)) )	
	set(reset, reset[,which(Patient=='M10607')], 'drop.day', as.Date('1992-09-15'))
	set(reset, reset[,which(Patient=='M13163')], 'drop.day', as.Date('1996-12-13'))
	set(reset, reset[,which(Patient=='M34362')], 'drop.day', as.Date('2007-11-29'))
	
	check2		<- subset(check, lRNA.drop.before.ART>0.5 & lRNA.before.ART<4 & AnyT_T1_Acc=='NAccMD')
	if(verbose)	cat(paste("\nPatients lRNA.drop.before.ART>0.5 & AnyT_T1_Acc=='NAccMD', n=",check2[,length(unique(Patient))]))
	print(check2[,unique(Patient)])
	#print(check2,n=500)
	#	checked manually: reset AnyT_T1 to first PosRNA with VL < 3.44
	tmp			<- check2[, list(AnyT_T1=AnyT_T1[1], drop.day=min(AnyT_T1[1], PosRNA[which(lRNA<3.44)[1]])), by='Patient']
	reset		<- rbind(reset, tmp)
	
	check3		<- subset(check, lRNA.drop.before.ART>0.5 & lRNA.before.ART<4 & AnyT_T1_Acc=='Acc' & drop.days2ART>0 & drop.days2ART<=30)
	if(verbose)	cat(paste("\nPatients lRNA.drop.before.ART>0.5 & AnyT_T1_Acc=='Acc', n=",check3[,length(unique(Patient))]))
	print(check3[,unique(Patient)])
	#print(check3,n=500)
	#	M13105 1996-10-14     400 2.602		(don t set min to >3.1)
	tmp			<- check3[, list(AnyT_T1=AnyT_T1[1], drop.day=min(AnyT_T1[1], PosRNA[which(lRNA<3.09)[1]])), by='Patient']		
	set(tmp, tmp[,which(Patient=='M11331')], 'drop.day', as.Date('1996-01-18'))
	set(tmp, tmp[,which(Patient=='M13105')], 'drop.day', as.Date('1996-10-14'))
	set(tmp, tmp[,which(Patient=='M12882')], 'drop.day', as.Date('1997-01-14'))
	set(tmp, tmp[,which(Patient=='M15067')], 'drop.day', as.Date('1998-02-24'))
	set(tmp, tmp[,which(Patient=='M15166')], 'drop.day', as.Date('1998-03-20'))
	set(tmp, tmp[,which(Patient=='M15178')], 'drop.day', as.Date('2000-04-04'))
	set(tmp, tmp[,which(Patient=='M15692')], 'drop.day', as.Date('1997-05-02'))
	set(tmp, tmp[,which(Patient=='M17154')], 'drop.day', as.Date('1999-10-20'))
	set(tmp, tmp[,which(Patient=='M19100')], 'drop.day', as.Date('2000-10-17'))
	set(tmp, tmp[,which(Patient=='M19105')], 'drop.day', as.Date('2000-10-17'))
	set(tmp, tmp[,which(Patient=='M19208')], 'drop.day', as.Date('2000-11-22'))
	set(tmp, tmp[,which(Patient=='M19233')], 'drop.day', as.Date('2000-11-22'))
	set(tmp, tmp[,which(Patient=='M20373')], 'drop.day', as.Date('2001-07-19'))
	set(tmp, tmp[,which(Patient=='M26537')], 'drop.day', as.Date('2002-08-20'))
	set(tmp, tmp[,which(Patient=='M28831')], 'drop.day', as.Date('2003-09-01'))
	set(tmp, tmp[,which(Patient=='M29510')], 'drop.day', as.Date('2004-03-25'))
	set(tmp, tmp[,which(Patient=='M29641')], 'drop.day', as.Date('2006-09-20'))
	set(tmp, tmp[,which(Patient=='M32779')], 'drop.day', as.Date('2006-05-17'))
	set(tmp, tmp[,which(Patient=='M33857')], 'drop.day', as.Date('2007-04-05'))
	set(tmp, tmp[,which(Patient=='M34098')], 'drop.day', as.Date('2011-05-04'))
	set(tmp, tmp[,which(Patient=='M34416')], 'drop.day', as.Date('2008-08-28'))
	set(tmp, tmp[,which(Patient=='M34453')], 'drop.day', as.Date('2008-01-14'))
	set(tmp, tmp[,which(Patient=='M34587')], 'drop.day', as.Date('2007-10-15'))	
	set(tmp, tmp[,which(Patient=='M34708')], 'drop.day', as.Date('2008-09-26'))
	set(tmp, tmp[,which(Patient=='M34926')], 'drop.day', as.Date('2009-05-22'))
	set(tmp, tmp[,which(Patient=='M35398')], 'drop.day', as.Date('2007-11-27'))
	set(tmp, tmp[,which(Patient=='M36082')], 'drop.day', as.Date('2011-01-07'))
	set(tmp, tmp[,which(Patient=='M36159')], 'drop.day', as.Date('2008-09-11'))
	set(tmp, tmp[,which(Patient=='M36619')], 'drop.day', as.Date('2010-09-28'))
	set(tmp, tmp[,which(Patient=='M37307')], 'drop.day', as.Date('2009-05-14'))
	set(tmp, tmp[,which(Patient=='M39134')], 'drop.day', as.Date('2011-05-31'))
	set(tmp, tmp[,which(Patient=='M39188')], 'drop.day', as.Date('2010-09-09'))
	set(tmp, tmp[,which(Patient=='M40166')], 'drop.day', as.Date('2011-12-15'))	
	set(tmp, tmp[,which(Patient=='M15184')], 'drop.day', as.Date('1998-03-20'))
	reset		<- rbind(reset, subset(tmp, !is.na(drop.day))) 
	
	check4		<- subset(check, lRNA.drop.before.ART>0.5 & lRNA.before.ART<4 & AnyT_T1_Acc=='Acc' & drop.days2ART>0 & drop.days2ART>30)
	if(verbose)	cat(paste("\nPatients lRNA.drop.before.ART>0.5 & AnyT_T1_Acc=='Acc' >30, n=",check4[,length(unique(Patient))]))
	print(check4[,unique(Patient)])
	#print(check4,n=500)		
	tmp			<- check4[, list(AnyT_T1=AnyT_T1[1], drop.day=min(AnyT_T1[1], PosRNA[which(lRNA<3.09)[1]])), by='Patient']
	set(tmp, tmp[,which(Patient=='M14785')], 'drop.day', as.Date('1997-05-14'))
	set(tmp, tmp[,which(Patient=='M14796')], 'drop.day', as.Date('1997-05-14'))
	set(tmp, tmp[,which(Patient=='M15151')], 'drop.day', as.Date('1998-03-20'))
	set(tmp, tmp[,which(Patient=='M15715')], 'drop.day', as.Date('1997-07-08'))
	set(tmp, tmp[,which(Patient=='M16342')], 'drop.day', as.Date('1998-09-08'))
	set(tmp, tmp[,which(Patient=='M16378')], 'drop.day', as.Date('1997-10-07'))
	set(tmp, tmp[,which(Patient=='M16701')], 'drop.day', as.Date('1998-05-07'))
	set(tmp, tmp[,which(Patient=='M18948')], 'drop.day', as.Date('2000-08-15'))
	set(tmp, tmp[,which(Patient=='M20546')], 'drop.day', as.Date('2009-03-03'))
	set(tmp, tmp[,which(Patient=='M20732')], 'drop.day', as.Date('2002-07-01'))
	set(tmp, tmp[,which(Patient=='M28930')], 'drop.day', as.Date('2004-01-13'))
	set(tmp, tmp[,which(Patient=='M32484')], 'drop.day', as.Date('2009-06-24'))
	set(tmp, tmp[,which(Patient=='M35687')], 'drop.day', as.Date('2008-09-02'))
	set(tmp, tmp[,which(Patient=='M37131')], 'drop.day', as.Date('2009-08-19'))
	set(tmp, tmp[,which(Patient=='M10623')], 'drop.day', as.Date('1995-02-23'))
	set(tmp, tmp[,which(Patient=='M17719')], 'drop.day', as.Date('1999-05-20'))	 
	reset		<- rbind(reset, subset(tmp, !is.na(drop.day)))
	#
	reset		<- subset( reset, AnyT_T1!=drop.day )
	setnames(reset, c('Patient','drop.day'), c('xPatient','new.AnyT_T1'))	
	#file<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/derived/ATHENA_2013_03_changed_ART_startdate.csv"
	#tmp			<- merge( subset(check, select=c(Patient, PosRNA, lRNA , AnyT_T1_Acc)), reset, by='Patient' )
	#write.csv(tmp, file=file)		
	#
	for( i in seq_len(nrow(reset)))
	{
		cat(paste('\nprocess',i))
		tmp	<- df.treatment[, which(Patient==reset[i,xPatient])]
		if(!length(tmp)) stop('could not find Patient')
		set(df.treatment, tmp, 'AnyT_T1', reset[i,new.AnyT_T1])		
		set(df.treatment, tmp[1], 'StartTime', reset[i,new.AnyT_T1])
	}
	
	set(df.treatment, df.treatment[,which(Patient=='M38411')][1], 'AnyT_T1', as.Date('2011-05-24 '))
	set(df.treatment, df.treatment[,which(Patient=='M38411')][1], 'StartTime', as.Date('2011-05-24 '))
	
	df<- df.treatment
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	save(df, file=file)
}
######################################################################################
project.hivc.Excel2dataframe.CD4<- function(dir.name= DATA, verbose=1)
{
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Immu.csv",sep='/')
	#file			<- paste(dir.name,"derived/ATHENA_2013_03_Immu_AllMSM.csv",sep='/')
	NA.time			<- c("","01/01/1911","11/11/1911","24/06/1923")		
	verbose			<- 1
	#read CD4 csv data file and preprocess
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
	#set(df, NULL, "Patient", factor(df[,Patient]))
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
							Patient=="M39055" & PosCD4=="2012-02-06" & CD4A==6850	|
							Patient=="M36408" & PosCD4=="2011-12-09" & CD4A==3701		#unclear
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
	#	remove too small entries and digit entries
	#	
	df		<- subset(df, CD4A>10)	
	#subset(df, (CD4A%%1)!=0 & ((CD4A*10)%%1)==0 & CD4A<170)
	df		<- subset(df, (CD4A%%1)==0 )
	#
	#	remove duplicate entries 
	#
	setkey(df, Patient, PosCD4, CD4A)
	tmp		<- which(duplicated(df))
	if(verbose) cat(paste("\nnumber of duplicated entries by Patient, PosCD4, CD4A, DELETE n=",length(tmp)))
	df		<- unique(df)
	setkey(df, Patient, PosCD4)
	tmp		<- which(duplicated(df))	
	if('Source'%in%colnames(df))
		tmp		<- sapply(tmp, function(i) c(i-1,i)[ which.min(df[c(i-1,i), Source]) ] )
	if(!'Source'%in%colnames(df))
		tmp		<- sapply(tmp, function(i) c(i-1,i)[ which.max(df[c(i-1,i), CD4A]) ] )
	if(verbose) cat(paste("\nnumber of duplicated entries by Patient, PosCD4, DELETE WITH LOWER SOURCE n=",length(tmp)))
	df		<- df[-tmp,]
	#
	#	remove consecutive jumps >1e3 
	#
	setkey(df, Patient, PosCD4)
	tmp		<- subset( df[, {  tmp<- which( length(CD4A)>1 & diff(CD4A)>1000 ); list(PosCD4=ifelse(length(tmp), as.character(PosCD4[tmp+1]), NA_character_), CD4A=ifelse(length(tmp), CD4A[tmp+1], NA_real_)) }, by='Patient'], !is.na(PosCD4))
	if(verbose) cat(paste("\nnumber of consecutive jumps >1e3 DELETE n=",nrow(tmp)))
	tmp		<- sapply(seq_len(nrow(tmp)), function(i)		df[,which(Patient==tmp[i,Patient] & PosCD4==tmp[i,PosCD4] & CD4A==tmp[i,CD4A])])
	df		<- df[-tmp,]
	#
	#	remove consecutive jumps >7e2 within 60days 
	#
	tmp		<- subset( df[, {  tmp<- which( length(CD4A)>1 & diff(CD4A)>700 & difftime(PosCD4[-1],PosCD4[-length(PosCD4)],units='days')<60 ); list(PosCD4=ifelse(length(tmp), as.character(PosCD4[tmp+1]), NA_character_), CD4A=ifelse(length(tmp), CD4A[tmp+1], NA_real_)) }, by='Patient'], !is.na(PosCD4))
	if(verbose) cat(paste("\nnumber of consecutive jumps >7e2 within 60days DELETE n=",nrow(tmp)))
	tmp		<- sapply(seq_len(nrow(tmp)), function(i)		df[,which(Patient==tmp[i,Patient] & PosCD4==tmp[i,PosCD4] & CD4A==tmp[i,CD4A])])
	df		<- df[-tmp,]
	#
	#	check above 1700 manually
	#
	tmp<- merge(df, subset(df, CD4A>1700, Patient), by="Patient")
	setkey(tmp, Patient, PosCD4)
	tmp<- tmp[,	{
					z<- which(CD4A>1700)
					list(PosCD4_T1= min(PosCD4), CD4_T1= CD4A[1], CD4.med= median(CD4A), CD4.q= quantile(CD4A, p=0.9), CD4.h=CD4A[z], PosCD4=PosCD4[z])
			},by="Patient"]	
	#leave those with high CD4 at start for now
	tmp		<- subset(tmp, CD4_T1!=CD4.h)	 	
	#leave those with last non-zero digit for now
	tmp		<- subset(tmp, ((CD4.h/10)%%1)==0 )		
	#	divide by 10
	#subset(tmp, CD4.q*3<CD4.h)
	z		<- df[, which( CD4A>1700 & Patient%in%subset(tmp, CD4.q*3<CD4.h)[,Patient])]
	set(df, z, 'CD4A', df[z,CD4A]/10)
	tmp		<- subset(tmp, CD4.q*3>=CD4.h)
	#	remove
	z		<- df[, which( CD4A>1700 & Patient%in%subset(tmp, CD4.q*2<CD4.h)[,Patient])]
	df		<- df[-z]
	tmp		<- subset(tmp, CD4.q*2>=CD4.h)
	#	remove
	z		<- df[, which( CD4A>1700 & Patient%in%subset(tmp, CD4.q<CD4.h & CD4.q<1200)[,Patient])]
	df		<- df[-z]
	tmp		<- subset(tmp, !(CD4.q<CD4.h & CD4.q<1200))
	if(verbose) cat(paste("\nstill high entries NOT DELETED n=",nrow(tmp)))
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
	df		<- df[, 	{
							z<- which.min(PosCD4)
							list(PosCD4=PosCD4, CD4=CD4A, PosCD4_T1=PosCD4[z], CD4_T1=CD4A[z] ) 	
						},by=Patient]
	file		<- paste(substr(file, 1, nchar(file)-3),'R',sep='')	
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)		
}
######################################################################################
project.hivc.Excel2dataframe.Viro<- function()		
{
	file.treatment		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	file				<- paste(dir.name,"derived/ATHENA_2013_03_Viro.csv",sep='/')
	file.treatment		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens_AllMSM.R",sep='/')
	file				<- paste(dir.name,"derived/ATHENA_2013_03_Viro_AllMSM.csv",sep='/')
	
	verbose				<- 1
	dir.name			<- DATA
	DB.locktime			<- HIVC.db.locktime
	
	#need for checking of VL data
	
	load(file.treatment)
	df.treat			<- subset(df, select=c(Patient, StartTime, StopTime, AnyT_T1, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel))
	
	NA.time					<- c("","01/01/1911","11/11/1911","24/06/1923")		
	RNA.min					<- 50	#400	#seems to be standard value
	RNA.min.b4T				<- 400
	RNA.stdvl.udetect.aTS	<- 1e3
	RNA.max					<- 5e6
	lRNA.min.infectious		<- log10(1e4)
	lRNA.min.early			<- log10(1e5)
	lRNA.max.b4early		<- log10(2e4)
	
	#read VIROLOGY csv data file and preprocess	
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
	#set(df, NULL, "Patient", factor(df[,Patient]))
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
	tmp		<- which( df[, Patient=="M12735" & as.character(PosRNA)=="1986-04-29"] )
	set(df,tmp,"PosRNA",NA)		
	tmp		<- which( df[, Patient=="M29967" & as.character(PosRNA)=="1998-07-09"] )
	set(df,tmp,"PosRNA",NA)		
	# remove is.na(PosRNA) and !is.na(RNA)
	df		<- subset(df, !is.na(PosRNA) & !is.na(RNA))
	if(verbose)		cat(paste("\nnumber of entries with !is.na(PosRNA) & !is.na(RNA), n=",nrow(df)))	
	#
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
	#	remove RNA values that are probably real but have wrong units
	#
	df	<- merge( df, unique(subset(df.treat, select=c(Patient, AnyT_T1))), by='Patient', all.x=1 )
	tmp	<- df[, which((is.na(AnyT_T1) | difftime(AnyT_T1, PosRNA, units='days')>=0) & ((RNA%%1)!=0) ) ]
	if(verbose)		cat(paste("\nremove RNA values with .XXX before ART start, n=",length(tmp)))	
	set(df, tmp, 'RNA', NA_real_)
	df	<- subset(df, !is.na(RNA))
	#
	#	remove RNA values < 400 if undetectable and before ART start 
	#	
	tmp	<- df[, which(Undetectable=='Yes' & RNA<RNA.min.b4T & (is.na(AnyT_T1) | difftime(AnyT_T1, PosRNA, units='days')>0) )]
	if(verbose)		cat(paste("\nremove RNA values < ",RNA.min.b4T," if undetectable and before ART start, n=",length(tmp)))
	set(df, tmp, 'RNA', NA_real_)
	df	<- subset(df, !is.na(RNA))
	#
	#	remove RNA values < 400 if detectable and before ART start 
	#	
	tmp	<- df[, which(Undetectable=='No' & RNA<RNA.min.b4T & (is.na(AnyT_T1) | difftime(AnyT_T1, PosRNA, units='days')>0) )]
	if(verbose)		cat(paste("\nremove RNA values < 400 if detectable and before ART start, n=",length(tmp)))
	set(df, tmp, 'RNA', NA_real_)
	df	<- subset(df, !is.na(RNA))
	#
	#	set undetectable & 1e3<RNA<1e4 after treatment start to RNA
	#	
	tmp		<- which( df[, Undetectable=="Yes" & RNA>RNA.stdvl.udetect.aTS & RNA<1e4 & difftime(PosRNA, AnyT_T1, units='days')>=0] )
	if(verbose)		cat(paste("\nsetting Undetectable=='Yes' and RNA<1e4 and PosRNA>AnyT_T1 to Undetectable=='No', n=",length(tmp)))
	set(df,tmp,"Undetectable","No")
	#
	#	set undetectable & RNA<1e4 after treatment start to 1e3
	#	
	tmp		<- which( df[, Undetectable=="Yes" & RNA<=RNA.stdvl.udetect.aTS & difftime(PosRNA, AnyT_T1, units='days')>=0] )
	if(verbose)		cat(paste("\nsetting Undetectable=='Yes' and RNA<1e4 and PosRNA>AnyT_T1 to Undetectable=='No' and RNA.min, n=",length(tmp)))
	set(df,tmp,"Undetectable","No")
	set(df,tmp,"RNA",RNA.stdvl.udetect.aTS)
	#
	#	remove undetectable RNA (before treatment start) that are before the first detectable RNA
	#
	tmp		<- df[, list(select= !all(Undetectable=='Yes')), by='Patient']
	if(verbose)		cat(paste("\nPatients with all Undetectable RNA, DELETE n=",nrow(subset(tmp, !select))))
	df		<- merge( df, subset(tmp, select, Patient), by='Patient' )
	df		<- df[,  {
							tmp<- seq.int(which(Undetectable=='No')[1], length(Undetectable)) 
							lapply(.SD,'[',tmp)
						},by='Patient']
	if(verbose)		cat(paste("\nnumber of entries after removing undetectable RNA before any detectable RNA, n=",nrow(df)))
	if(verbose)		cat(paste("\nnumber of remaining undetectable RNA, n=",df[,length(which(Undetectable=='Yes'))]))
	#M38913, M38722, M38411, M37285, M37211, M35710, M33521, M31924, M28495, M15900
	# 	remove further suspicious entries
	tmp		<- df[, which(Patient=='M38913' & PosRNA==as.Date("2010-07-01"))]	
	tmp		<- c(tmp, df[, which(Patient=='M38722' & PosRNA==as.Date("2010-05-11"))])
	tmp		<- c(tmp, df[, which(Patient=='M37211' & PosRNA==as.Date("2009-08-18"))])
	tmp		<- c(tmp, df[, which(Patient=='M33521' & PosRNA==as.Date("2006-12-18"))])
	tmp		<- c(tmp, df[, which(Patient=='M15900' & PosRNA==as.Date("1997-01-13"))])
	set(df, tmp, 'RNA', NA_real_)
	df	<- subset(df, !is.na(RNA))
	#
	#	remove undetectable if RNA<1e4 before ART
	#
	tmp		<- df[, which(Undetectable=='Yes' & RNA<=1e4 & (is.na(AnyT_T1) | difftime(AnyT_T1,PosRNA,units='days')>0 ) ) ]
	if(verbose)		cat(paste("\nremove undetectable RNA with RNA<1e4 before ART start, n=",length(tmp)))
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	if(verbose)		cat(paste("\nnumber of remaining undetectable RNA, n=",df[,length(which(Undetectable=='Yes'))]))
	print(subset(df, Undetectable=='Yes'))
	if(verbose)		cat(paste("\nremove these unclear RNA too"))
	df		<- subset(df, Undetectable=='No')
	#
	#	set RNA<RNA.min to RNA.min
	#		
	tmp<- which( df[, RNA<RNA.min] )
	if(verbose)		cat(paste("\nsetting RNA<RNA.min to RNA.min=",RNA.min,", n=",length(tmp)))		
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
	#	remove duplicate entries
	#
	setkey(df, Patient, PosRNA, RNA)
	df		<- unique(df)
	setkey(df, Patient, PosRNA)
	tmp		<- which(duplicated(df))
	tmp		<- sapply(tmp, function(i)		c(i-1,i)[which.min(df[c(i-1,i),Source])]	)
	df		<- df[-tmp,]			
	#			
	#
	tmp		<- merge(unique(subset(df, RNA>5e6 & PosRNA>AnyT_T1, Patient)), df, by="Patient")
	tmp		<- which(df[, Patient%in%c("M12736","M27885","M14799","M12612") & RNA>5e6])
	set(df, tmp, "RNA", NA_real_)	
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
	#tmp[, print(data.table(Patient,RNA,PosRNA,Undetectable,StartTime,StopTime,TrI,TrCh.failure,TrCh.adherence,TrCh.patrel)), by="Patient"]
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
	#if(verbose) cat(paste("\nremoving duplicate entries"))
	#df		<- df[, list(RNA=mean(RNA)), by=c('Patient','PosRNA')]
	#if(verbose) cat(paste("\nnumber of entries after duplicates removed, n=",nrow(df)))
	setkey(df, Patient, PosRNA)
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
	
	file		<- paste(substr(file, 1, nchar(file)-3),'R',sep='')	
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)	
}
######################################################################################
project.hivc.Excel2dataframe.Patients<- function(dir.name= DATA, min.seq.len=21, verbose=1)
{
	file				<- paste(dir.name,"derived/ATHENA_2013_03_Patient.csv",sep='/')
	file				<- paste(dir.name,"derived/ATHENA_2013_03_Patient_AllMSM.csv",sep='/')
	
	NA.Acute			<- c(NA,9)
	NA.CountryInfection	<- c(NA,"")
	NA.CountryBorn		<- c(NA,"",'XX')
	NA.RegionOrigin		<- c(NA,"",'XX')
	NA.Subtype			<- c(NA,"")
	NA.time				<- c("","01/01/1911","11/11/1911")		
	NA.transmission		<- 900
	#read PATIENT csv data file and preprocess	
	df					<- read.csv(file, stringsAsFactors=FALSE)									
	df$isDead			<- as.numeric( df[,"DateDied"]!="")
	tmp					<- which(df[,"Transmission"]==NA.transmission)
	if(length(tmp))
		df[tmp,"Transmission"]<- NA
	
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
	
	df<- data.table(df)
	setnames(df, c("MyDateNeg1_Acc","MyDatePos1_Acc","AcuteInfection","Transmission","HospitalRegion"), c("NegT_Acc","PosT_Acc","isAcute","Trm","RegionHospital"))
	#set(df, NULL, "Patient", factor(df[,Patient]))
	set(df, which( df[,isAcute%in%NA.Acute] ), "isAcute", NA )		
	set(df, NULL, "isAcute", factor(df[,isAcute], levels=c(0,1,2), labels=c("No","Yes","Maybe")) )
	set(df, NULL, "Sex", factor(df[,Sex], levels=c(1,2), labels=c("M","F")))
	set(df, which( df[,Subtype%in%NA.Subtype] ), "Subtype", NA_character_ )
	set(df, NULL, "Subtype", factor(df[,Subtype]))
	set(df, which( df[,CountryBorn%in%NA.CountryBorn] ), "CountryBorn", NA_character_ )	
	set(df, NULL, "CountryBorn", factor(df[,CountryBorn]))
	set(df, which( df[,CountryInfection%in%NA.CountryInfection] ), "CountryInfection", NA_character_ )		
	set(df, NULL, "CountryInfection", factor(df[,CountryInfection]))
	set(df, which( df[,RegionOrigin%in%NA.RegionOrigin] ), "RegionOrigin", NA_character_ )	
	set(df, NULL, "RegionOrigin", factor(df[,RegionOrigin]))
	
	set(df, NULL, "NegT_Acc", factor(df[,NegT_Acc], levels=c(0,1), labels=c("No","Yes")))
	set(df, NULL, "PosT_Acc", factor(df[,PosT_Acc], levels=c(0,1), labels=c("No","Yes")))		
	set(df, NULL, "isDead", factor(df[,isDead], levels=c(0,1), labels=c("No","Yes")))
	set(df, NULL, "Trm", factor(df[, Trm], levels=c(100, 101,  102,  202, 103,  104,  105,  106,  107, 108,  110), labels= c("MSM","BI","HET","HETfa","IDU","BLOOD","NEEACC", "PREG", "BREAST", "OTH", "SXCH")) )
	set(df, NULL, "RegionHospital", factor(df[,RegionHospital], levels=c(1,2,3,4,5,6), labels=c("Amst","N","E","S","W","Curu")))
	setkey(df,Patient)
	str(df)
	
	#reset NegT date
	tmp	<- df[, which(Patient=='M40895' & MyDateNeg1=='2011-04-06')]	
 	set(df, tmp, 'MyDateNeg1', as.Date('2010-04-06'))
	set(df, tmp, 'NegT_Acc', 'No')
	tmp	<- df[, which(Patient=='M29536' & MyDateNeg1=='1994-07-07')]	
	set(df, tmp, 'MyDateNeg1', as.Date('1994-07-02'))	
	tmp	<- df[, which(Patient=='M14759' & MyDateNeg1=='1991-09-15')]
	set(df, tmp, 'MyDateNeg1', as.Date('1991-08-015'))	
	tmp	<- df[, which(Patient=='M29967' & MyDateNeg1=='2000-07-15')]
	set(df, tmp, 'MyDateNeg1', NA)
	set(df, tmp, 'NegT_Acc', NA_character_)	
	tmp	<- df[, which(Patient=='M30654' & MyDateNeg1=='2004-05-15')]
	set(df, tmp, 'MyDateNeg1', NA)
	set(df, tmp, 'NegT_Acc', NA_character_)	
	tmp	<- df[, which(Patient=='M35339' & MyDateNeg1=='2007-08-15')]
	set(df, tmp, 'MyDateNeg1', NA)
	set(df, tmp, 'NegT_Acc', NA_character_)
	tmp	<- df[, which(Patient=='M35513' & MyDateNeg1=='2007-10-15')]
	set(df, tmp, 'MyDateNeg1', NA)
	set(df, tmp, 'NegT_Acc', NA_character_)
	
	
	file		<- paste(substr(file, 1, nchar(file)-3),'R',sep='')	
	if(verbose) cat(paste("\nsave to", file))	
	save(df, file=file)	
}
######################################################################################
project.hivc.Excel2dataframe.Sequences<- function(dir.name= DATA)
{
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.csv",sep='/')
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Sequences_AllMSM.csv",sep='/')
	#read SEQUENCE csv data file and preprocess				
	verbose			<- 1
	names.GeneCode	<- c("PROT","RT")
	NA.DateRes		<- as.Date("1911-11-11")
	NA.Seq			<- c('')
	min.seq.len		<- 21
	proc.GeneCode	<- c(1,2)
	
	df				<- read.csv(file, stringsAsFactors=FALSE)
	df[,"DateRes"]	<- as.Date(df[,"DateRes"], format="%d/%m/%Y")	
	nok.idx<- which( df[,"DateRes"]==NA.DateRes )	
	if(verbose) cat(paste("\nrange of DateRes is",paste(range(df[,"DateRes"], na.rm=1),collapse=', ')))
	if(verbose) cat(paste("\nentries with missing DateRes, n=", length(nok.idx)))					
	if(verbose) cat(paste("\nentries with missing DateRes, SampleCode", paste(df[nok.idx,"SampleCode"],collapse=', ')))
	df[nok.idx,"DateRes"]					<- NA	
	df[ which(df[,"Sequence"]%in%NA.Seq),]	<- NA 
	tmp										<- nchar(df[,"Sequence"])
	tmp										<- which(tmp<min.seq.len)
	if(verbose)	cat(paste("\ndiscarding sequences with insufficient length, n=",length(tmp),'\n'))
	df										<- df[-tmp,]
	
	
	df				<- lapply(proc.GeneCode,function(gene)
			{
				cat(paste("\nprocess GeneCode", gene))
				tmp				<- df[ df[,"GeneCode"]==gene, c("Patient","SampleCode","DateRes","Sequence"), drop=0 ]								
				tmp
			})		
	names(df)	<- names.GeneCode
	
	file		<- paste(substr(file, 1, nchar(file)-3),'R',sep='')
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)
	quit("no")	
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
project.hivc.tables.fixup.m2Bwmx<- function()
{
	tps	<- 1:4
	tmp	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3da_tablesCLU_m2Bwmx.tp'
	#tmp	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3ea_tablesCLU_m2Bwmx.tp'
	#tmp	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3fa_tablesCLU_m2Bwmx.tp'
	#tmp	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3ga_tablesCLU_m2Bwmx.tp'

	dummy	<- sapply(tps, function(tp)
			{
				file	<- paste(tmp, tp, '.R',sep='')
				load(file)
				set(ans$cens.table, NULL, 'factor', ans$cens.table[,gsub('SuA','ART.suA', factor)])
				set(ans$cens.table, NULL, 'factor2', ans$cens.table[,gsub('SuA','ART.suA', factor2)])
				set(ans$risk.table, NULL, 'factor', ans$risk.table[,gsub('SuA','ART.suA', factor)])
				set(ans$adj.clu, NULL, 'factor', ans$adj.clu[,gsub('SuA','ART.suA', factor)])
				set(ans$adj.seq, NULL, 'factor', ans$adj.seq[,gsub('SuA','ART.suA', factor)])
				cat(paste('\n save to file',file))
				save(ans, file=file)			
			})
	
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
	nrc.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=nrc.clu.pre)
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
	resume			<- 0
	argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume)
	argv			<<- unlist(strsplit(argv,' '))
	nsh.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=nsh.clu.pre, with.plot=0)
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
project.ukca.demographics<- function()
{
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/tblBASIC.csv'
	df.b	<- read.csv(f, stringsAsFactors=FALSE)
	df.b	<- as.data.table(df.b)
	set(df.b, NULL, 'NUCSEQ_D', hivc.db.Date2numeric(df.b[, as.Date(NUCSEQ_D, "%d/%m/%Y")]))
	set(df.b, NULL, 'HIVPOS_D', hivc.db.Date2numeric(df.b[, as.Date(HIVPOS_D, "%d/%m/%Y")]))
	set(df.b, NULL, 'BIRTH_D', hivc.db.Date2numeric(df.b[, as.Date(BIRTH_D, "%d/%m/%Y")]))
	set(df.b, NULL, 'DEATH_D', hivc.db.Date2numeric(df.b[, as.Date(DEATH_D, "%d/%m/%Y")]))
	set(df.b, NULL, 'LTALIV_D', hivc.db.Date2numeric(df.b[, as.Date(LTALIV_D, "%d/%m/%Y")]))	
	set(df.b, NULL, 'SEROCO_D', hivc.db.Date2numeric(df.b[, as.Date(SEROCO_D, "%d/%m/%Y")]))
	set(df.b, NULL, 'MODE', df.b[, factor(MODE, levels=c(1,2,3,4,5,6,7,8,90,99), labels=c('MSM','IDU','MSM+IDU','HAEM','5','HET','HET+IDU','8','OTH','NA'))])
	
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/tblLAB_RES.csv'
	df.s	<- read.csv(f, stringsAsFactors=FALSE)
	df.s	<- as.data.table(df.s)
	set(df.s, NULL, 'NUCSEQ_D', hivc.db.Date2numeric(df.s[, as.Date(NUCSEQ_D, "%Y-%m-%d")]))
	set(df.s, NULL, 'SEQ_DT', hivc.db.Date2numeric(df.s[, as.Date(SEQ_DT, "%Y-%m-%d")]))
	
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/tblLAB_CD4.csv'
	df.c	<- read.csv(f, stringsAsFactors=FALSE)
	df.c	<- as.data.table(df.c)
	set(df.c, NULL, 'CD4_D', hivc.db.Date2numeric(df.c[, as.Date(CD4_D, "%Y-%m-%d")]))
	
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/tblLAB_RNA.csv'
	df.v	<- read.csv(f, stringsAsFactors=FALSE)
	df.v	<- as.data.table(df.v)
	set(df.v, NULL, 'RNA_D', hivc.db.Date2numeric(df.v[, as.Date(RNA_D, "%Y-%m-%d")]))
	set(df.v, NULL, 'RNA_Vl10', df.v[, log10(RNA_V)])
	
	#individuals with sequence
	ggplot(df.s, aes(x=NUCSEQ_D)) + geom_histogram(binwidth=1) +
			scale_x_continuous(breaks=seq(1985,2020,1)) + scale_y_continuous(breaks=seq(0,1600,200))
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/expl_seqtotal.pdf'
	ggsave(f, w=14,h=6)	
	
	#seroconverters with sequence
	df.a	<- merge(df.b, subset(df.s, select=c(PATIENT, NUCSEQ_ID)), by='PATIENT')
	setkey(df.a, PATIENT)
	df.a<- unique(df.a)	
	df.a[, table(MODE)]
	#MODE
    #	MSM     IDU 	MSM+IDU    HAEM       5     HET HET+IDU       8     OTH      NA 
   	#	6638    1313      90      46      75    3222     227       6     173     387 
	ggplot(subset(df.a, !is.na(SEROCO_D)), aes(x=SEROCO_D)) + geom_histogram(binwidth=1) +
			scale_x_continuous(breaks=seq(1980,2020,1)) + scale_y_continuous(breaks=seq(0,500,50))
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/expl_seroconvtotal.pdf'
	ggsave(f, w=14,h=6)
	#	last alive if not dead
	tmp		<- subset(df.a, DEATH_Y==0)	
	tmp[, NOTSEEN:=max(LTALIV_D, na.rm=1)-LTALIV_D]
	ggplot(tmp, aes(x=NOTSEEN)) + geom_histogram(binwidth=0.5)
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/expl_ltalive.pdf'
	ggsave(f, w=8,h=6)
	
	#seroconverters with CD4 within 1 year of HIV+
	tmp		<- merge(df.a, df.c,by='PATIENT')
	setkey(tmp, PATIENT, CD4_D)
	tmp		<- unique(tmp)
	tmp		<- subset(tmp, abs(CD4_D-HIVPOS_D)<1 )
	setkey(tmp, PATIENT)
	tmp		<- unique(tmp)
	#	3505
	#ggplot(tmp, aes(x=CD4_V)) + geom_histogram(binwidth=200) + scale_x_continuous(breaks=seq(0,4000,200)) 
	
	#seroconverters with VL
	tmp		<- merge(df.a, df.v,by='PATIENT')
	setkey(tmp, PATIENT, RNA_D)
	tmp		<- unique(tmp)
	tmp[, list(n=length(RNA_Vl10)), by='PATIENT'][, median(n)]
	subset(tmp, !is.na(SEROCO_D))[, length(unique(PATIENT))]
	#	3870		
	subset(tmp, !is.na(SEROCO_D) & !is.nan(RNA_Vl10) & RNA_D-HIVPOS_D>0.5 & RNA_D-HIVPOS_D<2)[, length(unique(PATIENT))]
	#	609
	ggplot(subset(tmp, !is.na(SEROCO_D) & !is.nan(RNA_Vl10) & RNA_D-HIVPOS_D>0.5 & RNA_D-HIVPOS_D<2), aes(x=RNA_Vl10)) + 
			scale_x_continuous(breaks=seq(0,8,1)) +
			geom_histogram(binwidth=0.25)
	f		<- '~/duke/2013_HIV_Hue/UKCA_1309/data/original_Cascade_1307/expl_seroconvSPVL.pdf'
	ggsave(f, w=8,h=8)
	
	

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
project.athena.Niaetal.similar<- function()
{
	require(data.table)
	require(ape)
	indir						<- paste(DATA,"tmp",sep='/')		
	indircov					<- paste(DATA,"derived",sep='/')
	outdir						<- paste(DATA,"NiaMSC",sep='/')	
	infile.cov					<- paste(indircov,"ATHENA_2013_03_AllSeqPatientCovariates.R",sep='/')
	infile.viro					<- paste(indircov,"ATHENA_2013_03_Viro.R",sep='/')
	infile.immu					<- paste(indircov,"ATHENA_2013_03_Immu.R",sep='/')
	infile.treatment			<- paste(indircov,"ATHENA_2013_03_Regimens.R",sep='/')
	infile.seq					<- paste(indir,'/',"ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_",gsub('/',':',"Wed_Dec_18_11:37:00_2013"),".R",sep='')
	#
	#	load patient RNA
	#
	load(infile.viro)
	set(df, NULL, 'Patient', df[, as.character(Patient)])
	df.viro				<- df
	#
	#	load patient CD4
	#
	load(infile.immu)
	set(df, NULL, 'Patient', df[, as.character(Patient)])
	df.immu				<- df
	#
	#	load patient regimen
	#
	load(infile.treatment)
	set(df, NULL, 'Patient', df[, as.character(Patient)])
	df.treatment		<- df		
	#
	#	load all covariates
	#
	load(infile.cov)
	set(df.all, NULL, 'Patient', df.all[, as.character(Patient)])
	df.demo				<- df.all		
	#
	#	load sequences
	#
	load(infile.seq)	
	#
	#	anonymize
	#
	tmp			<- unique(subset( df.demo, select=Patient))
	tmp[, Patient.new:= paste('P',seq_len(nrow(tmp)),sep='')]
	df.demo	<- merge(df.demo, tmp, by='Patient')
	df.demo	<- subset( df.demo, select=c(Patient, Patient.new, FASTASampleCode, PosSeqT, DateBorn, Sex, CountryBorn, RegionOrigin, DateDied, Subtype, isAcute, NegT,
						NegT_Acc, PosT, PosT_Acc, CountryInfection, Trm, DateLastContact, RegionHospital, DateFirstEverCDCC, isDead, AnyPos_T1, PoslRNA_T1, lRNA_T1, PosCD4_T1, CD4_T1, AnyT_T1 ))
	
	tmp			<- unique(subset( df.demo, select=c(Patient.new, FASTASampleCode, PosSeqT)))
	setkey(tmp, Patient.new, PosSeqT)
	tmp			<- cbind( tmp, tmp[, list(FASTASampleCode.new= paste('S',substr(Patient.new,2,nchar(Patient.new)),'-',seq_along(FASTASampleCode), sep='')),by=Patient.new])
	df.demo		<- merge(df.demo, subset(tmp, select=c(FASTASampleCode, FASTASampleCode.new)), by='FASTASampleCode')
	
	setnames( df.demo, c('Patient', 'Patient.new', 'FASTASampleCode', 'FASTASampleCode.new'), c('Patient.old', 'Patient', 'FASTASampleCode.old', 'FASTASampleCode') )	
	setnames( df.treatment, 'Patient', 'Patient.old')
	setnames( df.viro, 'Patient', 'Patient.old')
	setnames( df.immu, 'Patient', 'Patient.old')
	df.treatment	<- merge(unique(subset(df.demo, select=c(Patient, Patient.old))), df.treatment, by='Patient.old')
	df.viro			<- merge(unique(subset(df.demo, select=c(Patient, Patient.old))), df.viro, by='Patient.old')
	df.immu			<- merge(unique(subset(df.demo, select=c(Patient, Patient.old))), df.immu, by='Patient.old')
	
	tmp				<- subset(df.demo, select=c(FASTASampleCode, FASTASampleCode.old) )
	setkey(tmp, FASTASampleCode.old)
	rownames(seq.PROT.RT)	<- tmp[rownames(seq.PROT.RT), ][, FASTASampleCode]
	#	save with keys
	file			<- paste(outdir, 'data_with_keys.R',sep='/')
	save(df.demo, df.treatment, df.viro, df.immu, seq.PROT.RT, file=file)
	
	df.treatment[, Patient.old:=NULL]
	df.viro[, Patient.old:=NULL]
	df.immu[, Patient.old:=NULL]
	df.demo[, Patient.old:=NULL]
	df.demo[, FASTASampleCode.old:=NULL]
	#	save without keys
	file			<- paste(outdir, 'data_without_keys.R',sep='/')
	save(df.demo, df.treatment, df.viro, df.immu, seq.PROT.RT, file=file)

}
######################################################################################
project.athena.Xavieretal.similar<- function()
{
	require(data.table)
	require(ape)
	stop()
	indir					<- paste(DATA,"tmp",sep='/')		
	indircov				<- paste(DATA,"derived",sep='/')
	outdir					<- paste(DATA,"fisheretal",sep='/')
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	
	if(0)
	{
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"				
		clu.indir				<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/sasky_sdr06_-DR-RC-SH+LANL_alrh160'
		clu.infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		clu.insignat			<- "Wed_Dec_18_11:37:00_2013"
		clu.infilexml.opt		<- "alrh160"
		clu.infilexml.template	<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=2_sasky',sep='_')
	}
	if(1)
	{
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"				
		clu.indir				<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/gmrf_sdr06fr_-DR-RC-SH+LANL_um192rhU2080'
		clu.infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		clu.insignat			<- "Wed_Dec_18_11:37:00_2013"
		clu.infilexml.opt		<- "mph4clutx4tip"
		clu.infilexml.template	<- "um192rhU2080"	
		outfile					<- paste(infile,'gmrf',sep='_')
		outdir					<- paste(DATA,'Xavieretal',sep='/')
	}
	#	fixed input args
	opt.brl			<- "dist.brl.casc" 
	thresh.brl		<- 0.096
	thresh.bs		<- 0.8	
	#
	# 	recent msm 
	#
	load(paste(indircov,'/',infilecov,'.R',sep=''))	
	#
	#	load clustering msm
	#
	argv			<<- hivc.cmd.clustering.msm(indir, infiletree, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=1)
	argv			<<- unlist(strsplit(argv,' '))		
	msm				<- hivc.prog.get.clustering.MSM()
	clumsm.info		<- msm$df.cluinfo
	#	update clumsm.info with latest values from df.all
	clumsm.info		<- merge( df.all, subset(clumsm.info, select=c(FASTASampleCode, cluster, Node, clu.npat, clu.ntip, clu.nFrgnInfection, clu.fPossAcute, clu.AnyPos_T1, clu.bwpat.medbrl)), by='FASTASampleCode')
	#
	set(clumsm.info, NULL, 'PosSeqT', hivc.db.Date2numeric(clumsm.info[,PosSeqT]))
	set(clumsm.info, NULL, 'DateBorn', hivc.db.Date2numeric(clumsm.info[,DateBorn]))
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
	#	make selection of pilot clusters
	cluster.pilot	<- unique( subset(clumsm.info, Patient%in%c('M32179','M28593','M37538','M15184','M12774'), cluster) )	
	clumsm.info		<- merge(clumsm.info, cluster.pilot, by='cluster')
	clumsm.info		<- subset(clumsm.info, select=c(cluster, FASTASampleCode, Patient,  DateBorn, NegT, AnyPos_T1, PosSeqT, AnyT_T1, DateLastContact, DateDied, isAcute))
	setkey(clumsm.info, cluster, Patient, FASTASampleCode)
	#
	#	
	tmp						<- unique(subset(clumsm.info, select=Patient))	
	file.viro				<- paste(indircov,"/ATHENA_2013_03_Viro.R",sep='/')
	file.immu				<- paste(indircov,"/ATHENA_2013_03_Immu.R",sep='/')
	file.treatment			<- paste(indircov,"/ATHENA_2013_03_Regimens.R",sep='/')		
	#	load patient RNA
	load(file.viro)
	clumsm.viro				<- merge(tmp, subset(df, select=c(Patient, PosRNA, lRNA)), by='Patient')
	set(clumsm.viro, NULL, 'PosRNA', hivc.db.Date2numeric(clumsm.viro[,PosRNA]))
	
	#	load patient CD4				
	load(file.immu)
	clumsm.cd4				<- merge(tmp, subset(df, select=c(Patient, PosCD4, CD4)), by='Patient')
	set(clumsm.cd4, NULL, 'PosCD4', hivc.db.Date2numeric(clumsm.cd4[,PosCD4]))
	set(clumsm.cd4, NULL, 'CD4', clumsm.cd4[,round(CD4)])
	
	#	load patient regimen
	load(file.treatment)
	clumsm.art				<- merge(tmp, subset(df, select=c(Patient, AnyT_T1, StartTime, StopTime, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel)), by='Patient')
	set(clumsm.art, NULL, 'StartTime', hivc.db.Date2numeric(clumsm.art[,StartTime]))
	set(clumsm.art, NULL, 'StopTime', hivc.db.Date2numeric(clumsm.art[,StopTime]))
	set(clumsm.art, NULL, 'AnyT_T1', hivc.db.Date2numeric(clumsm.art[,AnyT_T1]))
	setnames(clumsm.art, c('StartTime','StopTime','TrI','TrCh.failure','TrCh.adherence','TrCh.patrel'), c('ART.TS','ART.TE','ART.Intrpt','ART.Fail','ART.Adh','ART.PatDec'))

	#	anonymize	
	tmp			<- unique(subset( clumsm.info, select=Patient))
	tmp[, Patient.new:= paste('P',seq_len(nrow(tmp)),sep='')]
	clumsm.info	<- merge(clumsm.info, tmp, by='Patient')
	tmp			<- unique(subset( clumsm.info, select=c(Patient.new, FASTASampleCode, PosSeqT)))
	setkey(tmp, Patient.new, PosSeqT)
	tmp			<- cbind( tmp, tmp[, list(FASTASampleCode.new= paste('S',substr(Patient.new,2,nchar(Patient.new)),'-',seq_along(FASTASampleCode), sep='')),by=Patient.new]) 
	clumsm.info	<- merge(clumsm.info, subset(tmp, select=c(FASTASampleCode, FASTASampleCode.new)), by='FASTASampleCode')
	setnames( clumsm.info, c('Patient', 'Patient.new', 'FASTASampleCode', 'FASTASampleCode.new'), c('Patient.old', 'Patient', 'FASTASampleCode.old', 'FASTASampleCode') )	
	setnames( clumsm.art, 'Patient', 'Patient.old')
	setnames( clumsm.viro, 'Patient', 'Patient.old')
	setnames( clumsm.cd4, 'Patient', 'Patient.old')
	clumsm.art	<- merge(unique(subset(clumsm.info, select=c(Patient, Patient.old))), clumsm.art, by='Patient.old')
	clumsm.viro	<- merge(unique(subset(clumsm.info, select=c(Patient, Patient.old))), clumsm.viro, by='Patient.old')
	clumsm.cd4	<- merge(unique(subset(clumsm.info, select=c(Patient, Patient.old))), clumsm.cd4, by='Patient.old')
	clumsm.art[, Patient.old:=NULL]
	clumsm.viro[, Patient.old:=NULL]
	clumsm.cd4[, Patient.old:=NULL]
	
	#	load MAP cluster topology with branch lengths
	files		<- list.files(clu.indir)
	if(!length(files))	stop('no input files matching criteria')
	files		<- files[ sapply(files, function(x) grepl(clu.infile, x, fixed=1) & grepl(gsub('/',':',clu.insignat), x, fixed=1) & grepl(clu.infilexml.opt, x, fixed=1) & grepl(clu.infilexml.template,x, fixed=1) & grepl('_cluposterior_[0-9]+',x) & grepl('R$',x) ) ]		
	if(!length(files))	stop('no input files matching criteria')
	tmp			<- regmatches( files, regexpr('_cluposterior_[0-9]+',files)) 
	cluster		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
	file.info	<- data.table(file=files, cluster=cluster)
	setkey(file.info, cluster)
	#	subset of clusters: those in df.tpairs.select	
	file.info	<- merge(file.info, unique(subset(clumsm.info, select=cluster)), by='cluster')
	clumsm.info	<- merge(clumsm.info, subset(file.info, select=cluster), by='cluster')				
	#	combine dated cluster phylogenies
	clumsm.ph	<- lapply(seq_len(nrow(file.info)), function(i)
		{				
			#	load dated cluster phylogenies
			file					<- paste(clu.indir, file.info[i,file], sep='/')
			cat(paste('\nload file=',file,'i=',i))
			tmp						<- load(file)
			topo.map				<- mph.clu.dtopo[which.max(freq),]					
			tmp						<- which( grepl( paste('mph.i=',topo.map[,mph.i],'_',sep=''), names(ph.consensus) ) )
			if(length(tmp)!=1)	stop('unexpected ph.consensus index')
			topo.map.ph				<- ph.consensus[[ tmp ]]
			topo.map.SA				<- subset(mph.SA.cnt, equal.to==topo.map[,mph.i])
			set(topo.map.SA, NULL, 'SA.freq', topo.map.SA[,SA.freq/n])													
			topo.map.nodectime		<- subset(mph.node.ctime, equal.to==topo.map[,mph.i])
			topo.map				<- merge(subset(topo.map, select=c(cluster, dens)), subset(topo.map.SA, select=c(cluster, tip, SA.freq)), by='cluster')
			topo.map.nodectime		<- subset(topo.map.nodectime, select=c(cluster, node, q, cdf, pdf))
			#
			tmp						<- merge( data.table( BEASTlabel=topo.map.ph$tip.label, FASTASampleCode.old= sapply( strsplit(topo.map.ph$tip.label, '_'), '[[', 6 ) ), subset(clumsm.info, select=c(FASTASampleCode, FASTASampleCode.old)), by='FASTASampleCode.old')
			setkey(tmp, BEASTlabel)
			topo.map.ph$tip.label	<- tmp[topo.map.ph$tip.label, ][,FASTASampleCode]
			topo.map.ph$node.dated	<- topo.map.nodectime
			topo.map.ph			
		})

	#	save KEY used for anonymization
	clumsm.info
	file		<- paste(outdir,'/',outfile,'_key.R',sep='')
	save(clumsm.info, file=file)
	clumsm.info[, Patient.old:=NULL]
	clumsm.info[, FASTASampleCode.old:=NULL]
	clumsm.info	<- subset(clumsm.info, select=c(cluster, Patient, FASTASampleCode, DateBorn, NegT, AnyPos_T1,  PosSeqT,  AnyT_T1, DateLastContact, DateDied, isAcute ))
	
	#	save info per cluster
	setkey(clumsm.info, FASTASampleCode)	
	setkey(clumsm.art, Patient)
	setkey(clumsm.cd4, Patient)
	setkey(clumsm.viro, Patient)
	tmp			<- lapply(seq_len(nrow(file.info)), function(i)
		{
			clu			<- list()
			clu$ph		<- clumsm.ph[[i]]
			clu$info	<- clumsm.info[ clu$ph$tip.label, ]
			clu$art		<- clumsm.art[ clu$info[, unique(Patient)], ]
			clu$viro	<- clumsm.viro[ clu$info[, unique(Patient)], ]
			clu$cd4		<- clumsm.cd4[ clu$info[, unique(Patient)], ]
			
			file		<- paste(outdir,'/',outfile,'_clu',clu$info[,unique(cluster)],'_raw.rda',sep='')
			cat(paste('\nsave to file', file))
			save(file=file, clu)			
		})
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

project.hivc.examlclock<- function()
{
	#
	#	root to tip divergence
	#	
	require(ape)
	require(adephylo)
	indir				<- paste(DATA,"tmp",sep='/')
	infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	infiletree			<- paste(infile,"examlbs500",sep="_")
	indircov			<- paste(DATA,"derived",sep='/')
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	insignat			<- "Wed_Dec_18_11:37:00_2013"		
	file				<- paste(indir,'/',infiletree,'_',gsub('/',':',insignat),".R",sep='')
	load(file)
	tmp					<- node.depth.edgelength(ph)
	df					<- data.table(FASTASampleCode= ph$tip.label, height=tmp[seq_len(Ntip(ph))])
	
	file				<- paste(indircov,'/',infilecov,".R",sep='')
	load(file)
	df			<- merge( subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyPos_T1, AnyT_T1)), df, by='FASTASampleCode' )
	set(df, NULL, 'PosSeqT', hivc.db.Date2numeric(df[,PosSeqT]))
	set(df, NULL, 'AnyPos_T1', hivc.db.Date2numeric(df[,AnyPos_T1]))
	set(df, NULL, 'AnyT_T1', hivc.db.Date2numeric(df[,AnyT_T1]))
	df			<- subset(df, !is.na(PosSeqT))
	ggplot(df, aes(x=PosSeqT, y=height)) + geom_point()
	
	df.clock	<- lm(height~PosSeqT, data=as.data.frame(df))
	summary(df.clock)
	#Adjusted R-squared: 0.1056
	#             Estimate Std. Error t value Pr(>|t|)    
	#(Intercept) -5.690e+00  1.865e-01  -30.51   <2e-16 ***
	#PosSeqT      2.917e-03  9.297e-05   31.37   <2e-16 ***
	dfbh		<- df[, {
						z<- which.min(PosSeqT)
						list(FASTASampleCode=FASTASampleCode[z], PosSeqT= PosSeqT[z], AnyT_T1=AnyT_T1[z], height=height[z])
					}, by='Patient']
	dfbh.clock	<- lm(height~PosSeqT, data=as.data.frame(dfbh))
	#dfbh.clock	<- lm(height~PosSeqT, data=as.data.frame(subset(dfbh, PosSeqT>1996.6)))
	summary(dfbh.clock)		
	#Adjusted R-squared: 0.1271	    
	#(Intercept) -6.3631412  0.2165435  -29.39   <2e-16 ***
	#PosSeqT      0.0032520  0.0001079   30.13   <2e-16 ***
	dfbh[, b4T:= factor(is.na(AnyT_T1) | PosSeqT<AnyT_T1, labels=c('No','Yes'))]
	dfb4T.clock	<- lm(height~PosSeqT, data=as.data.frame(subset(dfbh, b4T)))
	summary(dfb4T.clock)
	#Adjusted R-squared: 0.105
	#(Intercept) -6.702884   0.303098  -22.11   <2e-16 ***
	#PosSeqT      0.003422
	ggplot(dfbh, aes(x=PosSeqT, y=height, colour=b4T)) + geom_point(alpha=0.75) +
		scale_x_continuous(breaks=seq(1980,2020,2)) +
		scale_colour_brewer(name='Sampled before ART start', palette='Set1') +
		stat_smooth(method="lm", colour='black') + stat_smooth(method='lm', data=subset(dfbh, b4T=='Yes'), colour="black", linetype=2) +
		labs(x='Sequence sampling date', y='root-to-tip divergence') +
		theme(legend.position=c(0,1), legend.justification=c(0,1))
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_clock.pdf'
	ggsave(file=file, w=15, h=5)
	#
	#	within host divergence
	#
	file				<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_nbrlraw3kaH'				
	tmp					<- project.athena.Fisheretal.Y.rawbrl(NULL, NULL, NULL, NULL, NULL, NULL, df.tpairs.tptn=NULL, save.file=paste(file, '.R', sep=''), resume=1, plot.file=NA, method.restrictTPtoRI=0)
	Y.rawbrl			<- tmp$tpairs
	Y.rawbrl.linked		<- tmp$linked
	Y.rawbrl.unlinked	<- tmp$unlinked
	#	compute sequence sampling times	to see if time between seq sampling times could be useful
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )	
	tmp					<- merge( data.table(FASTASampleCode=Y.rawbrl.linked[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, tmp, by='t.FASTASampleCode')
	#	order with respect to PosSeqT
	tmp					<- subset(Y.rawbrl.linked, t.PosSeqT>PosSeqT)
	setnames(tmp, c("t.FASTASampleCode","t.PosSeqT","t.AnyT_T1"), c("a.FASTASampleCode","a.PosSeqT","a.AnyT_T1"))
	setnames(tmp, c("FASTASampleCode","PosSeqT","AnyT_T1"), c("t.FASTASampleCode","t.PosSeqT","t.AnyT_T1"))
	setnames(tmp, c("a.FASTASampleCode","a.PosSeqT","a.AnyT_T1"), c("FASTASampleCode","PosSeqT","AnyT_T1"))
	Y.rawbrl.linked		<- rbind( subset(Y.rawbrl.linked, t.PosSeqT<=PosSeqT), tmp, use.names=TRUE )	
	set(Y.rawbrl.linked, NULL, 'PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,PosSeqT]))
	set(Y.rawbrl.linked, NULL, 't.PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,t.PosSeqT]))
	set(Y.rawbrl.linked, NULL, 'AnyT_T1', hivc.db.Date2numeric(Y.rawbrl.linked[,AnyT_T1]))
	set(Y.rawbrl.linked, NULL, 't.AnyT_T1', hivc.db.Date2numeric(Y.rawbrl.linked[,t.AnyT_T1]))	
	Y.rawbrl.linked[, b4T:= 'both']
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT<=AnyT_T1)], 'b4T', 'only.RI')
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT<=t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'only.T')
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'none')
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, data.table(b4T= Y.rawbrl.linked[, unique(b4T)], col=c('green','red','orange')), by='b4T' )	
	set(Y.rawbrl.linked, NULL, 'b4T', Y.rawbrl.linked[, factor(b4T)])
	Y.rawbrl.linked[, dt:= abs(PosSeqT-t.PosSeqT)]		
	Y.rawbrl.linked[, brlr:= brl/dt]
	#	add treatment complications: any of ART.P ART.F ART.I ART.A yes before PosSeqT
	if(1)
	{
		tmp				<- subset(df.treatment, select=c(Patient, StopTime, TrI, TrCh.failure, TrVL.failure))
		set(tmp, NULL, 'StopTime', hivc.db.Date2numeric(tmp[,StopTime]))
		tmp				<- merge(tmp, subset( Y.rawbrl.linked, select=c(Patient, PosSeqT) ), by='Patient', allow.cartesian=TRUE)
		tmp				<- subset(tmp, StopTime<=PosSeqT)	
		tmp				<- tmp[, list( ART.C=any( sapply(.SD, function(z) any(z=='Yes', na.rm=TRUE)) ) ), by=c('Patient','PosSeqT'), .SDcols=c('TrI','TrCh.failure','TrVL.failure')]
		Y.rawbrl.linked	<- merge( Y.rawbrl.linked, tmp, by=c('Patient','PosSeqT'), all.x=TRUE )		
		set(Y.rawbrl.linked, Y.rawbrl.linked[, which(is.na(ART.C))], 'ART.C', FALSE)		
		Y.rawbrl.linked[, table(b4T, ART.C)]
		Y.rawbrl.linked[, b4T2:=b4T]		
		tmp				<- Y.rawbrl.linked[, which(b4T!='both' & ART.C==FALSE)]
		set(Y.rawbrl.linked, tmp, 'b4T2', Y.rawbrl.linked[tmp, paste(b4T2,'ART.C.No',sep='.')])
		tmp				<- Y.rawbrl.linked[, which(b4T!='both' & ART.C==TRUE)]
		set(Y.rawbrl.linked, tmp, 'b4T2', Y.rawbrl.linked[tmp, paste(b4T2,'ART.C.yes',sep='.')])
		set(Y.rawbrl.linked, Y.rawbrl.linked[, which(b4T2=='none.ART.C.yes' | b4T2=='only.T.ART.C.yes')], 'b4T2', 'ART.C.yes')	
		#Y.rawbrl.linked	<- subset(Y.rawbrl.linked, b4T2!='none.ART.C.No' | (b4T2=='none.ART.C.No' & dt<6) )
		#	labels
		tmp				<- c('both sampling dates\nbefore cART initiation','sampling dates before and\nafter cART initiation,\nno interruption or failure','both sampling dates\nafter cART initiation,\nno interruption or failure','at least one sampling date\nafter cART initiation,\nwith interruption or failure')
		Y.rawbrl.linked	<- merge( Y.rawbrl.linked, data.table(b4T2= c('both','only.T.ART.C.No','none.ART.C.No','ART.C.yes'), b4T2.long=tmp), by='b4T2' )
		set(Y.rawbrl.linked, NULL, 'b4T2.long', Y.rawbrl.linked[, factor(b4T2.long, levels=tmp)])			
		#	explore clock
		Y.rawbrl.linked[, excluded:= factor( dt>brl.linked.max.dt | dt<brl.linked.min.dt | brlr>brl.linked.max.brlr, levels=c(FALSE,TRUE), labels=c('No','Yes'))]
		ggplot(Y.rawbrl.linked, aes(x=dt, y=brl, colour=b4T2.long)) + geom_point(size=1, aes(shape=excluded)) + facet_grid(. ~ b4T2.long) +			 
			scale_x_continuous(limits=c(0,10)) + scale_y_continuous(limits=c(0,0.35), breaks=seq(0,0.35,0.02)) + scale_colour_brewer(palette='Set1', guide=FALSE) + 					
			stat_smooth(aes(colour=NULL), colour='black', data=subset(Y.rawbrl.linked, excluded=='No')) +			
			labs(x="years between within-host sampling dates", y='within-host divergence')
	
	}
	if(1)
	{
		tmp				<- c('both sampling dates\nbefore ART start','sampling dates before and\nafter ART start','both sampling dates\nafter ART start') 
		Y.rawbrl.linked	<- merge( Y.rawbrl.linked, data.table(b4T= c('both','only.T','none'), b4T.long=factor( tmp, levels=tmp )), by='b4T' )
		set(Y.rawbrl.linked, NULL, 'b4T.long', Y.rawbrl.linked[, factor(b4T.long)])		
	}	
	
	setkey(Y.rawbrl.linked, b4T.long, dt)
	#
	#	explore clock
	#
	brl.linked.max.dt= 10; brl.linked.min.dt= 1; brl.linked.max.brlr=0.02
	
	tmp		<- subset(Y.rawbrl.linked, b4T=='both' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T.long))	
	tmpg.ZAGA.B	<- gamlss(as.formula('brl ~ bs(dt, degree=4)'), sigma.formula=as.formula('~ bs(dt, degree=4)'), nu.formula=as.formula('~ bs(dt, degree=4)'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
	tmp[, y.b:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)]
	tmp[, y.u:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.B, type='response', se.fit=TRUE)$se.fit]
	tmp[, y.l:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.B, type='response', se.fit=TRUE)$se.fit]
	ans		<- copy(tmp)
	tmp		<- subset(Y.rawbrl.linked, b4T=='only.T' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T.long))
	tmpg.ZAGA.O	<- gamlss(as.formula('brl ~ bs(dt, degree=4)'), sigma.formula=as.formula('~ bs(dt, degree=4)'), nu.formula=as.formula('~ bs(dt, degree=4)'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
	tmp[, y.b:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)]
	tmp[, y.u:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.O, type='response', se.fit=TRUE)$se.fit]
	tmp[, y.l:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.O, type='response', se.fit=TRUE)$se.fit]
	ans		<- rbind(ans, tmp)
	tmp		<- subset(Y.rawbrl.linked, b4T=='none' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T.long))
	tmpg.ZAGA.N	<- gamlss(as.formula('brl ~ bs(dt, degree=4)'), sigma.formula=as.formula('~ bs(dt, degree=4)'), nu.formula=as.formula('~ bs(dt, degree=4)'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
	tmp[, y.b:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)]
	tmp[, y.u:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.N, type='response', se.fit=TRUE)$se.fit]
	tmp[, y.l:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.N, type='response', se.fit=TRUE)$se.fit]
	ans		<- rbind(ans, tmp)
	Y.rawbrl.linked[, excluded:= factor( dt>brl.linked.max.dt | dt<brl.linked.min.dt | brlr>brl.linked.max.brlr, levels=c(FALSE,TRUE), labels=c('No','Yes'))]
	ggplot(Y.rawbrl.linked, aes(x=dt, y=brl, colour=b4T.long)) + 
			geom_point(size=1, aes(shape=excluded)) + scale_x_continuous(limits=c(0,10)) + scale_y_continuous(limits=c(0,0.35)) + scale_colour_brewer(palette='Set1', guide=FALSE) + facet_grid(. ~ b4T.long) +		
			#geom_line(aes(x=dt, y=y.b), colour='black', data=ans) +
			stat_smooth(aes(colour=NULL), colour='black', data=subset(Y.rawbrl.linked, excluded=='No')) +
			#geom_ribbon(aes(x=dt, ymin=y.l, ymax=y.u), alpha=0.2, data=ans, inherit.aes = FALSE) +
			labs(x="years between within-host sampling dates", y='within-host divergence')
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examl_clockwhexplore.pdf'
	ggsave(file=file, w=8, h=8)	
	#
	#	mean evol rate
	#
	tmp			<- subset(Y.rawbrl.linked, b4T=='both' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T.long))		
	tmpg.ZAGA.B	<- gamlss(as.formula('brl ~ dt-1'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
	tmp[, y.b:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)]
	tmp[, y.u:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.B, type='response', se.fit=TRUE)$se.fit]
	tmp[, y.l:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.B, type='response', se.fit=TRUE)$se.fit]
	ans			<- copy(tmp)
	tmp			<- subset(Y.rawbrl.linked, b4T=='only.T' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T.long))
	tmpg.ZAGA.O	<- gamlss(as.formula('brl ~ dt-1'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
	tmp[, y.b:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)]
	tmp[, y.u:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.O, type='response', se.fit=TRUE)$se.fit]
	tmp[, y.l:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.O, type='response', se.fit=TRUE)$se.fit]
	ans			<- rbind(ans, tmp)
	tmp			<- subset(Y.rawbrl.linked, b4T=='none' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T.long))
	tmpg.ZAGA.N	<- gamlss(as.formula('brl ~ dt-1'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
	tmp[, y.b:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)]
	tmp[, y.u:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.N, type='response', se.fit=TRUE)$se.fit]
	tmp[, y.l:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.N, type='response', se.fit=TRUE)$se.fit]
	ans			<- rbind(ans, tmp)
	#	
	Y.rawbrl.linked[, excluded:= factor( dt>brl.linked.max.dt | dt<brl.linked.min.dt | brlr>brl.linked.max.brlr, levels=c(FALSE,TRUE), labels=c('No','Yes'))]
	ggplot(Y.rawbrl.linked, aes(x=dt, y=brl, colour=b4T.long)) + 
			geom_point(size=1, aes(shape=excluded)) + scale_x_continuous(limits=c(0,10)) + scale_y_continuous(limits=c(0,0.35)) + scale_colour_brewer(palette='Set1', guide=FALSE) + facet_grid(. ~ b4T.long) +					
			stat_smooth(aes(colour=NULL), colour='black', data=subset(Y.rawbrl.linked, excluded=='No')) +
			geom_line(aes(x=dt, y=y.b), colour='black', linetype='dotdash', data=ans) +	
			#geom_ribbon(aes(x=dt, ymin=y.l, ymax=y.u), alpha=0.2, data=ans, inherit.aes = FALSE) +
			labs(x="years between within-host sampling dates", y='within-host divergence')  
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examl_clockwh.pdf'
	ggsave(file=file, w=8, h=8)
	if(1)
	{
		#
		#	mean evol rate
		#
		tmp			<- subset(Y.rawbrl.linked, b4T2=='both' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T2.long))		
		tmpg.ZAGA.B	<- gamlss(as.formula('brl ~ dt'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
		tmp[, y.b:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)]
		tmp[, y.u:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.B, type='response', se.fit=TRUE)$se.fit]
		tmp[, y.l:=predict(tmpg.ZAGA.B, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.B, type='response', se.fit=TRUE)$se.fit]
		ans			<- copy(tmp)
		tmp			<- subset(Y.rawbrl.linked, b4T2=='only.T.ART.C.No' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T2.long))
		tmpg.ZAGA.O	<- gamlss(as.formula('brl ~ dt'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
		tmp[, y.b:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)]
		tmp[, y.u:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.O, type='response', se.fit=TRUE)$se.fit]
		tmp[, y.l:=predict(tmpg.ZAGA.O, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.O, type='response', se.fit=TRUE)$se.fit]
		ans			<- rbind(ans, tmp)		
		tmp			<- subset(Y.rawbrl.linked, b4T2=='none.ART.C.No' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T2.long))
		tmpg.ZAGA.N	<- gamlss(as.formula('brl ~ dt'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
		tmp[, y.b:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)]
		tmp[, y.u:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.N, type='response', se.fit=TRUE)$se.fit]
		tmp[, y.l:=predict(tmpg.ZAGA.N, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.N, type='response', se.fit=TRUE)$se.fit]
		ans			<- rbind(ans, tmp)		
		tmp			<- subset(Y.rawbrl.linked, b4T2=='ART.C.yes' & dt>0 & dt>brl.linked.min.dt & dt<=brl.linked.max.dt & brlr<brl.linked.max.brlr, select=c(brl, dt, b4T2.long))
		tmpg.ZAGA.C	<- gamlss(as.formula('brl ~ dt'), sigma.formula=as.formula('~ dt'), nu.formula=as.formula('~ dt'), data=as.data.frame(tmp), family=ZAGA(mu.link='identity'), n.cyc = 40)
		tmp[, y.b:=predict(tmpg.ZAGA.C, type='response', se.fit=FALSE)]
		tmp[, y.u:=predict(tmpg.ZAGA.C, type='response', se.fit=FALSE)+2*predict(tmpg.ZAGA.C, type='response', se.fit=TRUE)$se.fit]
		tmp[, y.l:=predict(tmpg.ZAGA.C, type='response', se.fit=FALSE)-2*predict(tmpg.ZAGA.C, type='response', se.fit=TRUE)$se.fit]
		ans			<- rbind(ans, tmp)
		#	
		Y.rawbrl.linked[, excluded:= factor( dt>brl.linked.max.dt | dt<brl.linked.min.dt | brlr>brl.linked.max.brlr, levels=c(FALSE,TRUE), labels=c('No','Yes'))]
		ggplot(Y.rawbrl.linked, aes(x=dt, y=brl, colour=b4T2.long)) + 
				geom_point(size=1, aes(shape=excluded)) + facet_grid(. ~ b4T2.long) +
				scale_x_continuous(limits=c(0,10)) + scale_y_continuous(limits=c(0,0.35), breaks=seq(0,0.35,0.02)) + scale_colour_brewer(palette='Set2', guide=FALSE) +
				stat_smooth(aes(colour=NULL), colour='black', data=subset(Y.rawbrl.linked, excluded=='No')) +
				geom_line(aes(x=dt, y=y.b), colour='black', linetype='dotdash', data=ans) +					
				labs(x="years between within-host sampling dates", y='within-host divergence')  
		file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examl_clockwh_ARTC.pdf'
		ggsave(file=file, w=12, h=8)
		#
		sapply( list(tmpg.ZAGA.B, tmpg.ZAGA.O, tmpg.ZAGA.N, tmpg.ZAGA.C), Rsq )
		#0.23127742  0.04869653 0.10655731  0.14558462
		#	excluded	outside 1<dt<10 & brlr>0.02
		print(tmpg.ZAGA.B)	#dt		0.002663		#sigma	0.9496      -0.1229
		print(tmpg.ZAGA.O)	#dt		0.003733		#sigma	0.71293     -0.03297
		print(tmpg.ZAGA.N)	#dt  	0.005023		#sigma	0.50178      0.05858
		print(tmpg.ZAGA.C)	#dt  	0.00621			#sigma	0.56372     -0.06687 
	}
	
	
	sapply( list(tmpg.ZAGA.B, tmpg.ZAGA.O, tmpg.ZAGA.N), Rsq )
	#0.4411087 0.2959325 0.4299490
	#	excluded	only outside 1<dt<10
	print(tmpg.ZAGA.B)	#		0.001567
	print(tmpg.ZAGA.O)	#dt		0.002675
	print(tmpg.ZAGA.N)	#dt  	0.003985
	#	excluded	outside 1<dt<10 & brlr>0.01
	print(tmpg.ZAGA.B)	#		0.002177
	print(tmpg.ZAGA.O)	#dt		0.002849
	print(tmpg.ZAGA.N)	#dt  	0.004247 
	#	excluded	outside 1<dt<10 & brlr>0.02
	print(tmpg.ZAGA.B)	#		0.002663
	print(tmpg.ZAGA.O)	#dt		0.004007
	print(tmpg.ZAGA.N)	#dt  	0.006346 
	#	excluded	outside 1<dt<10 & brlr>0.025
	print(tmpg.ZAGA.B)	#		0.00275
	print(tmpg.ZAGA.O)	#dt		0.004131  
	print(tmpg.ZAGA.N)	#dt  	0.006926 
	
	#
	#	EXCLUDE
	#	
	Y.rawbrl.linked		<- subset( Y.rawbrl.linked, dt!=0 )
	cat(paste('\nY.rawbrl.linked all, n=',nrow(Y.rawbrl.linked)))
	Y.rawbrl.linked		<- subset(Y.rawbrl.linked,  brlr<=brl.linked.max.brlr)
	cat(paste('\nY.rawbrl.linked max brlr ',brl.linked.max.brlr,' excluded, n=',nrow(Y.rawbrl.linked)))
	Y.rawbrl.linked		<- subset(Y.rawbrl.linked,  dt<=brl.linked.max.dt)		
	cat(paste('\nY.rawbrl.linked max dt ',brl.linked.max.dt,' excluded, n=',nrow(Y.rawbrl.linked)))	
	Y.rawbrl.linked		<- subset(Y.rawbrl.linked,  dt>brl.linked.min.dt)		
	cat(paste('\nY.rawbrl.linked min dt ',brl.linked.min.dt,' excluded, n=',nrow(Y.rawbrl.linked)))
	Y.rawbrl.linked[, brlz:=brl]
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(brlz<1e-4)], 'brlz', 0)	
	#
	# wilcox test for same mean
	#
	test	<- list( 	both.only.T= 	wilcox.test( subset(Y.rawbrl.linked, b4T=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.T')[, log10(brl)] ),				
						both.none= 		wilcox.test( subset(Y.rawbrl.linked, b4T=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='none')[, log10(brl)] ),						
						none.only.T= 	wilcox.test( subset(Y.rawbrl.linked, b4T=='none')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.T')[, log10(brl)] )		)
	cat(paste('\nwilcox test for same mean'))
	print( sapply(test, '[[', 'p.value') )
	if(1)
	{
		test	<- list( 	both.only.T= 	wilcox.test( subset(Y.rawbrl.linked, b4T2=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T2=='only.T.ART.C.No')[, log10(brl)] ),
							both.none.CNo= 	wilcox.test( subset(Y.rawbrl.linked, b4T2=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T2=='none.ART.C.No')[, log10(brl)] ),
							both.CYes= 	wilcox.test( subset(Y.rawbrl.linked, b4T2=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T2=='ART.C.yes')[, log10(brl)] ),							
							both.none= 		wilcox.test( subset(Y.rawbrl.linked, b4T2=='only.T.ART.C.No')[, log10(brl)], subset(Y.rawbrl.linked, b4T2=='none.ART.C.No')[, log10(brl)] ),						
							none.only.T= 	wilcox.test( subset(Y.rawbrl.linked, b4T2=='none.ART.C.No')[, log10(brl)], subset(Y.rawbrl.linked, b4T2=='ART.C.yes')[, log10(brl)] )		)
		cat(paste('\nwilcox test for same mean'))
		print( sapply(test, '[[', 'p.value') )
		#	  both.only.T both.none.CNo     both.CYes     both.none   none.only.T 
 		#3.529053e-02  7.964509e-01  1.651674e-21  3.479906e-01  2.830199e-05 		
	}
	#wilcox test for same mean		excluded brl>0.02
	# 	 both.only.T    both.none  none.only.T 
	#2.071745e-04 1.279865e-22 1.322711e-17
	#								excluded brl>0.025
	# both.only.T    both.none  none.only.T 
	#2.238328e-04 6.624977e-24 4.412482e-20
	#
	# AIC between EXP GA ZAGA	-- AIC not comparable between brl and brlz
	#
	tmpd.B		<- subset(Y.rawbrl.linked, b4T=='both', select=c(brl, brlz, dt, b4T.long))
	tmpd.O		<- subset(Y.rawbrl.linked, b4T=='only.T', select=c(brl, brlz, dt, b4T.long))
	tmpd.N		<- subset(Y.rawbrl.linked, b4T=='none', select=c(brl, brlz, dt, b4T.long))
	nrow(tmpd.B)		#[1] 134	excluded <0.01			141		excluded <0.02		142		excluded <0.025
	nrow(tmpd.O)		#[1] 328							373							376
	nrow(tmpd.N)		#[1] 1576							2154						2257
	tmpg.E.B	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.B), family=EXP, n.cyc = 40)	
	tmpg.E.O	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.O), family=EXP, n.cyc = 40)
	tmpg.E.N	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.N), family=EXP, n.cyc = 40)	
	tmpg.GA.B	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.B), family=GA, n.cyc = 40)
	tmpg.GA.O	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.O), family=GA, n.cyc = 40)
	tmpg.GA.N	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.N), family=GA, n.cyc = 40)
	tmpg.ZAGA.B	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.B), family=ZAGA, n.cyc = 40)			
	tmpg.ZAGA.O	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.O), family=ZAGA, n.cyc = 40)		
	tmpg.ZAGA.N	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.N), family=ZAGA, n.cyc = 40)
	AIC(tmpg.E.B, tmpg.GA.B, tmpg.ZAGA.B, tmpg.E.O, tmpg.GA.O, tmpg.ZAGA.O, tmpg.E.N, tmpg.GA.N, tmpg.ZAGA.N)
	#            df         AIC
	#tmpg.GA.N    2 -13443.5388		tmpg.GA.O    2  -3006.6705		tmpg.GA.B    2  -1322.0507
	#tmpg.ZAGA.N  3 -13441.5387		tmpg.ZAGA.O  3  -3004.6703		tmpg.ZAGA.B  3  -1320.0505
	#tmpg.E.N     1 -12494.9765		tmpg.E.O     1  -2351.8905		tmpg.E.B     1   -925.7479		
	tmpp.ZAGA	<- sapply( list(tmpg.ZAGA.B, tmpg.ZAGA.O, tmpg.ZAGA.N), function(z){	c(mu=as.double(exp(z[['mu.coefficients']])), sigma=as.double(exp(z[['sigma.coefficients']])), nu=as.double(1/(1+exp(-z[['nu.coefficients']]))) ) 	})
	colnames(tmpp.ZAGA)	<- c('B','O','N')
	tmpp.GA	<- sapply( list(tmpg.GA.B, tmpg.GA.O, tmpg.GA.N), function(z){	c(mu=as.double(exp(z[['mu.coefficients']])), sigma=as.double(exp(z[['sigma.coefficients']])) ) 	})
	colnames(tmpp.GA)	<- c('B','O','N')	
	tmpp.E	<- as.matrix(t(sapply( list(tmpg.E.B, tmpg.E.O, tmpg.E.N), FUN=function(z){	c(mu=as.double(exp(z[['mu.coefficients']])) ) 	}, simplify=FALSE)))
	colnames(tmpp.E)	<- c('B','O','N')
	rownames(tmpp.E)	<- c('mu')	
	#
	# ks.test between EXP GA ZAGA	
	#	
	test		<- list(	E.B= ks.test(tmpd.B[,brl], pEXP, mu=unlist(tmpp.E['mu','B'])),
							E.O= ks.test(tmpd.O[,brl], pEXP, mu=unlist(tmpp.E['mu','O'])),
							E.N= ks.test(tmpd.N[,brl], pEXP, mu=unlist(tmpp.E['mu','N'])),
							GA.B= ks.test(tmpd.B[,brl], pGA, mu=unlist(tmpp.GA['mu','B']), sigma=unlist(tmpp.GA['sigma','B'])),
							GA.O= ks.test(tmpd.O[,brl], pGA, mu=unlist(tmpp.GA['mu','O']), sigma=unlist(tmpp.GA['sigma','O'])),
							GA.N= ks.test(tmpd.N[,brl], pGA, mu=unlist(tmpp.GA['mu','N']), sigma=unlist(tmpp.GA['sigma','N'])),
							ZAGA.B= ks.test(subset(tmpd.B, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','B']), sigma=unlist(tmpp.ZAGA['sigma','B']), nu=0),
							ZAGA.O= ks.test(subset(tmpd.B, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','O']), sigma=unlist(tmpp.ZAGA['sigma','O']), nu=0),
							ZAGA.N= ks.test(subset(tmpd.B, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','N']), sigma=unlist(tmpp.ZAGA['sigma','N']), nu=0)
							)	
	print( sapply(test, '[[', 'p.value') )
	#	excluded > 0.02
	#         E.B          E.O          E.N         GA.B         GA.O         GA.N       ZAGA.B       ZAGA.O       ZAGA.N 
	#0.000000e+00 0.000000e+00 0.000000e+00 5.000988e-08 1.911804e-13 0.000000e+00 4.456872e-01 1.358416e-06 4.620748e-13
	#	excluded > 0.025
	#         E.B          E.O          E.N         GA.B         GA.O         GA.N       ZAGA.B       ZAGA.O       ZAGA.N 
	#0.000000e+00 0.000000e+00 0.000000e+00 5.514628e-08 2.016165e-13 0.000000e+00 4.273192e-01 2.096416e-06 7.571721e-14 

	#0 is <2.2e-16
	#
	#	Parameters
	#
	tmpp		<- do.call('rbind',	list( 	data.table(d='ZAGA',b4T=c('B','O','N'),mu=unlist(tmpp.ZAGA['mu',]), sigma=tmpp.ZAGA['sigma',], nu=tmpp.ZAGA['nu',]),
											data.table(d='GA',b4T=c('B','O','N'),mu=unlist(tmpp.GA['mu',]), sigma=tmpp.GA['sigma',], nu=NA),
											data.table(d='EXP',b4T=c('B','O','N'),mu=unlist(tmpp.E['mu',]), sigma=NA, nu=NA)	)	)
	tmp			<- c('both sampling dates\nbefore ART start','sampling dates before and\nafter ART start','both sampling dates\nafter ART start')
	tmp			<- data.table(b4T= c('B','O','N'), b4T.long=factor( tmp, levels=tmp ))
	tmpp		<- merge( tmpp, tmp, by='b4T' )
	#	excluded > 0.02	
	#   b4T    d          mu     sigma        nu                                   b4T.long
	#1:   B ZAGA 0.013857413 0.7738736 0.3829787      both sampling dates\nbefore ART start
	#2:   B   GA 0.008553248 1.9267489        NA      both sampling dates\nbefore ART start
	#3:   B  EXP 0.008553248        NA        NA      both sampling dates\nbefore ART start
	#4:   N ZAGA 0.027019369 0.7546553 0.1429898       both sampling dates\nafter ART start
	#5:   N   GA 0.023157079 1.4154861        NA       both sampling dates\nafter ART start
	#6:   N  EXP 0.023157079        NA        NA       both sampling dates\nafter ART start
	#7:   O ZAGA 0.021351489 0.8010297 0.2975871 sampling dates before and\nafter ART start
	#8:   O   GA 0.014999959 1.8104345        NA sampling dates before and\nafter ART start
	#9:   O  EXP 0.014999959        NA        NA sampling dates before and\nafter ART start
	#	excluded > 0.025
	#b4T    d          mu     sigma        nu                                   b4T.long
	#1:   B ZAGA 0.014198191 0.7819489 0.3802817      both sampling dates\nbefore ART start
	#2:   B   GA 0.008801788 1.9263632        NA      both sampling dates\nbefore ART start
	#3:   B  EXP 0.008801788        NA        NA      both sampling dates\nbefore ART start
	#4:   N ZAGA 0.028178663 0.7533713 0.1364643       both sampling dates\nafter ART start
	#5:   N   GA 0.024334430 1.3967320        NA       both sampling dates\nafter ART start
	#6:   N  EXP 0.024334430        NA        NA       both sampling dates\nafter ART start
	#7:   O ZAGA 0.021496155 0.7984599 0.2952128 sampling dates before and\nafter ART start
	#8:   O   GA 0.015152594 1.8054715        NA sampling dates before and\nafter ART start
	#9:   O  EXP 0.015152594        NA        NA sampling dates before and\nafter ART start
	#	CDF plot
	#
	tmp			<- tmpp[, {
								if(b4T=='O')
									brl			<- tmpd.O[,brlz]
								if(b4T=='B')
									brl			<- tmpd.B[,brlz]
								if(b4T=='N')
									brl			<- tmpd.N[,brlz]	
								brl.cdf		<- data.table( brl=sort(brl), empirical=seq_along(brl)/length(brl) )[, list(empirical=tail(empirical, 1)), by='brl']
								if(d=='EXP')
									brl.cdf[, predicted:= pEXP( brl.cdf[,brl], mu=mu ) ]
								if(d=='GA')
									brl.cdf[, predicted:= pGA( brl.cdf[,brl], mu=mu, sigma=sigma ) ]
								if(d=='ZAGA')
									brl.cdf[, predicted:= pZAGA( brl.cdf[,brl], mu=mu, sigma=sigma, nu=nu ) ]
								brl.cdf				
							}, by=c('b4T.long', 'd')]
	tmp			<- melt(tmp, id=c('b4T.long','d', 'brl'))	
	ggplot( data=tmp, aes(x=brl, y=value, colour=variable)) + 
			labs(x='within-host divergence', y='c.d.f.') + 
			scale_colour_brewer(name='', palette='Paired') + 
			geom_point(data=subset(tmp, variable=='empirical'), show_guide=FALSE ) +
			geom_line() + facet_grid(d ~ b4T.long, scales='free_x', margins=FALSE) +
			theme(legend.key.size=unit(10,'mm'))
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examl_whcdf.pdf'
	ggsave(file=file, w=10,h=10)		
	#
	#
	#
	tmpd.B		<- subset(Y.rawbrl.linked, b4T2=='both', select=c(brl, brlz, dt, b4T2.long))
	tmpd.O		<- subset(Y.rawbrl.linked, b4T2=='only.T.ART.C.No', select=c(brl, brlz, dt, b4T2.long))
	tmpd.N		<- subset(Y.rawbrl.linked, b4T2=='none.ART.C.No', select=c(brl, brlz, dt, b4T2.long))
	tmpd.C		<- subset(Y.rawbrl.linked, b4T2=='ART.C.yes', select=c(brl, brlz, dt, b4T2.long))
	nrow(tmpd.B)		#141		excluded <0.02	
	nrow(tmpd.O)		#165						
	nrow(tmpd.N)		#28						
	nrow(tmpd.C)		#2334
	tmpg.E.B	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.B), family=EXP, n.cyc = 40)	
	tmpg.E.O	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.O), family=EXP, n.cyc = 40)
	tmpg.E.N	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.N), family=EXP, n.cyc = 40)
	tmpg.E.C	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.C), family=EXP, n.cyc = 40)
	tmpg.GA.B	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.B), family=GA, n.cyc = 40)
	tmpg.GA.O	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.O), family=GA, n.cyc = 40)
	tmpg.GA.N	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.N), family=GA, n.cyc = 40)
	tmpg.GA.C	<- gamlss(as.formula('brl ~ 1'), data=as.data.frame(tmpd.C), family=GA, n.cyc = 40)
	tmpg.ZAGA.B	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.B), family=ZAGA, n.cyc = 40)			
	tmpg.ZAGA.O	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.O), family=ZAGA, n.cyc = 40)		
	tmpg.ZAGA.N	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.N), family=ZAGA, n.cyc = 40)
	tmpg.ZAGA.C	<- gamlss(as.formula('brlz ~ 1'), data=as.data.frame(tmpd.C), family=ZAGA, n.cyc = 40)
	#
	tmpp.ZAGA	<- sapply( list(tmpg.ZAGA.B, tmpg.ZAGA.O, tmpg.ZAGA.N, tmpg.ZAGA.C), function(z){	c(mu=as.double(exp(z[['mu.coefficients']])), sigma=as.double(exp(z[['sigma.coefficients']])), nu=as.double(1/(1+exp(-z[['nu.coefficients']]))) ) 	})
	colnames(tmpp.ZAGA)	<- c('B','O','N','C')
	tmpp.GA	<- sapply( list(tmpg.GA.B, tmpg.GA.O, tmpg.GA.N, tmpg.GA.C), function(z){	c(mu=as.double(exp(z[['mu.coefficients']])), sigma=as.double(exp(z[['sigma.coefficients']])) ) 	})
	colnames(tmpp.GA)	<- c('B','O','N','C')	
	tmpp.E	<- as.matrix(t(sapply( list(tmpg.E.B, tmpg.E.O, tmpg.E.N, tmpg.E.C), FUN=function(z){	c(mu=as.double(exp(z[['mu.coefficients']])) ) 	}, simplify=FALSE)))
	colnames(tmpp.E)	<- c('B','O','N','C')
	rownames(tmpp.E)	<- c('mu')	
	tmpp		<- do.call('rbind',	list( 	data.table(d='ZAGA',b4T2=c('B','O','N','C'),mu=unlist(tmpp.ZAGA['mu',]), sigma=tmpp.ZAGA['sigma',], nu=tmpp.ZAGA['nu',]),
											data.table(d='GA',b4T2=c('B','O','N','C'),mu=unlist(tmpp.GA['mu',]), sigma=tmpp.GA['sigma',], nu=NA),
											data.table(d='EXP',b4T2=c('B','O','N','C'),mu=unlist(tmpp.E['mu',]), sigma=NA, nu=NA)	)	)
	tmp			<- c('both sampling dates\nbefore cART initiation','sampling dates before and\nafter cART initiation,\nno interruption or failure','both sampling dates\nafter cART initiation,\nno interruption or failure','at least one sampling date\nafter cART initiation,\nwith interruption or failure')						
	tmp			<- data.table(b4T2= c('B','O','N','C'), b4T2.long=factor( tmp, levels=tmp, labels=tmp ))
	tmpp		<- merge( tmpp, tmp, by='b4T2' )
	#   b4T2    d          mu     sigma        nu
 	#1:    B ZAGA 0.013857413 0.7738736 0.3829787
 	#2:    B   GA 0.008553248 1.9267489        NA
 	#3:    B  EXP 0.008553248        NA        NA
 	#4:    C ZAGA 0.026909219 0.7605336 0.1516710
 	#5:    C   GA 0.022829153 1.4444899        NA
 	#6:    C  EXP 0.022829153        NA        NA
 	#7:    N ZAGA 0.013327744 0.5810444 0.3571429
 	#8:    N   GA 0.008570276 1.8545254        NA
 	#9:    N  EXP 0.008570276        NA        NA
	#10:    O ZAGA 0.017742656 0.7320898 0.3333333
	#11:    O   GA 0.011831054 1.8523769        NA
	#12:    O  EXP 0.011831054        NA        NA
	#
	# ks.test between EXP GA ZAGA	
	#	
	test		<- list(	E.B= ks.test(tmpd.B[,brl], pEXP, mu=unlist(tmpp.E['mu','B'])),
							E.O= ks.test(tmpd.O[,brl], pEXP, mu=unlist(tmpp.E['mu','O'])),
							E.N= ks.test(tmpd.N[,brl], pEXP, mu=unlist(tmpp.E['mu','N'])),
							E.C= ks.test(tmpd.C[,brl], pEXP, mu=unlist(tmpp.E['mu','C'])),
							GA.B= ks.test(tmpd.B[,brl], pGA, mu=unlist(tmpp.GA['mu','B']), sigma=unlist(tmpp.GA['sigma','B'])),
							GA.O= ks.test(tmpd.O[,brl], pGA, mu=unlist(tmpp.GA['mu','O']), sigma=unlist(tmpp.GA['sigma','O'])),
							GA.N= ks.test(tmpd.N[,brl], pGA, mu=unlist(tmpp.GA['mu','N']), sigma=unlist(tmpp.GA['sigma','N'])),
							GA.C= ks.test(tmpd.C[,brl], pGA, mu=unlist(tmpp.GA['mu','C']), sigma=unlist(tmpp.GA['sigma','C'])),
							ZAGA.B= ks.test(subset(tmpd.B, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','B']), sigma=unlist(tmpp.ZAGA['sigma','B']), nu=0),
							ZAGA.O= ks.test(subset(tmpd.O, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','O']), sigma=unlist(tmpp.ZAGA['sigma','O']), nu=0),
							ZAGA.N= ks.test(subset(tmpd.N, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','N']), sigma=unlist(tmpp.ZAGA['sigma','N']), nu=0),
							ZAGA.C= ks.test(subset(tmpd.C, brlz>0)[,brlz], pZAGA, mu=unlist(tmpp.ZAGA['mu','C']), sigma=unlist(tmpp.ZAGA['sigma','C']), nu=0)		)	
	print( sapply(test, '[[', 'p.value') )
	#          E.B          E.O          E.N          E.C         GA.B         GA.O         GA.N         GA.C       ZAGA.B 
	#0.000000e+00 3.330669e-16 1.651107e-03 0.000000e+00 5.000988e-08 1.990942e-07 2.632056e-02 0.000000e+00 4.456872e-01 
    #  ZAGA.O       ZAGA.N       ZAGA.C 
	#6.337170e-01 3.623706e-01 8.584151e-03
	
	#	CDF plot
	#
	tmp			<- tmpp[, {
				if(b4T2=='O')
					brl			<- tmpd.O[,brlz]
				if(b4T2=='B')
					brl			<- tmpd.B[,brlz]
				if(b4T2=='N')
					brl			<- tmpd.N[,brlz]	
				if(b4T2=='C')
					brl			<- tmpd.C[,brlz]									
				brl.cdf		<- data.table( brl=sort(brl), empirical=seq_along(brl)/length(brl) )[, list(empirical=tail(empirical, 1)), by='brl']
				if(d=='EXP')
					brl.cdf[, predicted:= pEXP( brl.cdf[,brl], mu=mu ) ]
				if(d=='GA')
					brl.cdf[, predicted:= pGA( brl.cdf[,brl], mu=mu, sigma=sigma ) ]
				if(d=='ZAGA')
					brl.cdf[, predicted:= pZAGA( brl.cdf[,brl], mu=mu, sigma=sigma, nu=nu ) ]
				brl.cdf				
			}, by=c('b4T2.long', 'd')]
	tmp			<- melt(tmp, id=c('b4T2.long','d', 'brl'))	
	ggplot( data=tmp, aes(x=brl, y=value, colour=variable)) + 
			labs(x='within-host divergence', y='c.d.f.') + 
			scale_colour_brewer(name='', palette='Paired') + 
			geom_point(data=subset(tmp, variable=='empirical'), show_guide=FALSE ) +
			geom_line() + facet_grid(d ~ b4T2.long, scales='free_x', margins=FALSE) +
			theme(legend.key.size=unit(10,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5))
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examl_whcdf_ARTC.pdf'
	ggsave(file=file, w=12,h=10)		
	
	#
	#	QQ plot -- TODO can call stat_qq, need qZAGA etc and params only, and slope + int
	#
	tmp			<- tmpp[, {
				if(b4T=='O')
					brl			<- tmpd.O[,brlz]
				if(b4T=='B')
					brl			<- tmpd.B[,brlz]
				if(b4T=='N')
					brl			<- tmpd.N[,brlz]	
				if(b4T=='C')
					brl			<- tmpd.C[,brlz]					
				brl.cdf		<- data.table( brl=sort(brl), empirical=seq_along(brl)/length(brl) )[, list(empirical=tail(empirical, 1)), by='brl']
				if(d=='EXP')
				{
					brl.cdf[, predicted:= qEXP( brl.cdf[,brl], mu=mu ) ]
					qt 	<- qEXP(c(.25, .75), mu=mu )
				}
				if(d=='GA')
				{
					brl.cdf[, predicted:= qGA( brl.cdf[,brl], mu=mu, sigma=sigma ) ]
					qt 	<- qGA(c(.25, .75), mu=mu, sigma=sigma )
				}
				if(d=='ZAGA')
				{
					brl.cdf[, predicted:= qZAGA( brl.cdf[,brl], mu=mu, sigma=sigma, nu=nu ) ]										
					qt 	<- qZAGA(c(.25, .75), mu=mu, sigma=sigma, nu=nu )											
				}
				qv 		<- quantile(brl, c(.25, .75))
				slope 	<- diff(qv)/diff(qt)
				int 	<- qv[1] - slope * qt[1]
				brl.cdf				
			}, by=c('b4T.long', 'd')]
	#qplot(sample=brl, geom = "point", stat = "qq", pch='both sampling dates\nbefore cART initiation', distribution = qZAGA, dparams = list(mu=exp(ml.zaga.one[['mu.coefficients']]), sigma=exp(ml.zaga.one[['sigma.coefficients']]), nu=1/(1+exp(-ml.zaga.one[['nu.coefficients']])))) + 
	#		geom_abline(slope = slope, intercept = int) + labs(x='Zero-inflated Gamma', y='empirical') +
	#		theme(legend.justification=c(0,1), legend.position=c(0,1), legend.key.size=unit(13,'mm'), legend.title = element_blank())
	#ggsave(file=paste(substr(plot.file.one,1,nchar(plot.file.one)-4),'_zagaqq','.pdf',sep=''), w=4,h=6)		
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

hivc.project.remove.resistancemut.save<- function()
{
	dr		<- as.data.table( read.csv( paste( CODE.HOME,"/data/IAS_primarydrugresistance_201303.csv",sep='' ), stringsAsFactors=F ) )	
	dr[,Alignment.nuc.pos:= (Gene.codon.number-1)*3+Gene.HXB2pos ]		
	dr		<- dr[,	{	tmp<- unique(Mutant); list(Mutant=tmp, Gene.codon.number=Gene.codon.number[1], Wild.type=Wild.type[1], DR.name=DR.name[1])	}, by=Alignment.nuc.pos]
	#select nucleotide codes that are consistent with drug resistance mutants
	nt2aa	<- as.data.table( read.csv( paste( CODE.HOME,"/data/standard_nt_code.csv",sep='' ), stringsAsFactors=F ) )
	setnames(nt2aa,c("AA","NTs"),c("Mutant","Mutant.NTs"))
	nt2aa	<- subset(nt2aa, select=c(Mutant,Mutant.NTs))
	dr		<- merge(dr, nt2aa, all.x=1, by="Mutant", allow.cartesian=TRUE)
	setkey(dr, "Alignment.nuc.pos")
	#print(dr, nrows=250)
	dr		<- subset(dr, select=c(Alignment.nuc.pos, Mutant.NTs, DR.name))
	set(dr, NULL, "Mutant.NTs", tolower(dr[,Mutant.NTs]))
	IAS_primarydrugresistance_201303	<- data.frame(dr)
	save(IAS_primarydrugresistance_201303, file=paste( CODE.HOME,"/data/IAS_primarydrugresistance_201303.rda",sep='' ))	
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

