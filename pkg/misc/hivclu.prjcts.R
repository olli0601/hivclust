
#' @export
HIVC.COUNTRY.TABLE<- data.table(matrix(c(	"AG","ANTIGUA AND BARBUDA",		"AL","ALBANIA",		"AO","ANGOLA",	"AR","ARGENTINA",	"AT","AUSTRIA",		"AU","AUSTRALIA",				
						"BB","BARBADOS",	"BE","BELGIUM",	"BG","BULGARIA",	"BR","BRAZIL",	"CA","CANADA",	"CH","SWITZERLAND",		"CI","COTE D'IVOIRE",	"CL","CHILE",	
						"CN","CHINA",	"CO","COLOMBIA",	"CU","CUBA",	"CV","CAPE VERDE",	"CY","CYPRUS",	"CZ","CZECH REPUBLIC",	"DE","GERMANY",	"DK","DENMARK",	
						"DO","DOMINICAN REPUBLIC",	"DZ","ALGERIA",		"EC","ECUADOR",	"ES","SPAIN",	"FI","FINLAND",	"FR","FRANCE",	"GA","GABON",	"GB","UNITED KINGDOM",	
						"GD","GRENADA",	"GE","GIORGIA",	"GH","GHANA",	"GR","GREECE",	"GW","GUINEA-BISSAU",	"GY","GUYANA",	"HN","HONDURAS", "HT","HAITI",	
						"HU","HUNGARY",	"IL","ISRAEL",	"IN","INDIA",	"IT","ITALY",	"JM","JAMAICA",	"JP","JAPAN",	"KE","KENYA",	"KR","SOUTH KOREA",	"LB","LEBANON",				
						"LU","LUXEMBOURG",	"LV","LATVIA",	"MA","MOROCCO",	"ME","MONTENEGRO",	"MG","MADAGASCAR",	"ML","MALI",	"MN","MONGOLIA",	"MX","MEXICO",	"MY","MALAYSIA",
						"NG","NIGERIA",	"NL","NETHERLANDS",	"NO","NORWAY",	"PA","PANAMA",	"PE","PERU",	"PH","PHILIPPINES",	"PL","POLAND",	"PR","PUERTO RICO",	"PT","PORTUGAL",			
						"PY","PARAGUAY",	"RO","ROMANIA",	"RS","SERBIA",	"RU","RUSSIAN FEDERATION",	"SC","SEYCHELLES",	"SD","SUDAN",	"SE","SWEDEN",	"SG","SINGAPORE",	
						"SI","SLOVENIA",	"SK","SLOVAKIA",	"SN","SENEGAL",	"SR","SURINAME",	"SV","EL SALVADOR",	"TH","THAILAND",	"TT","TRINIDAD AND TOBAGO",	"TW","TAIWAN",			
						"UA","UKRAINE",	"UG","UGANDA",	"US","UNITED STATES",	"UY","URUGUAY",	"VE","VENEZUELA",	"VN","VIETNAM",	"YE","YEMEN","ZA","SOUTH AFRICA"),ncol=2,byrow=1,dimnames=list(c(),c("key","country"))))

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
					#stop()
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
					stop()
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
project.hivc.collectpatientdata<- function(dir.name= DATA, verbose=1, resume=0)
{	
	require(data.table)		
	
	#input files generated with "project.hivc.Excel2dataframe"
	file.seq		<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.R",sep='/')
	file.patient	<- paste(dir.name,"derived/ATHENA_2013_03_Patients.R",sep='/')
	file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	file.immu		<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
	file.treatment	<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	
	#compute file AllSeqPatientCovariates
	file.out		<- paste(dir.name,"derived/ATHENA_2013_03_AllSeqPatientCovariates.R",sep='/')
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
		
		#add Patient data
		if(verbose)		cat(paste("\nadding patient data"))
		load(file.patient)		
		df.all	<- merge(df.all, df, all.x=1, by="Patient")
		setnames(df.all, c("MyDateNeg1","MyDatePos1"), c("NegT","PosT"))
		#add Treatment dates
		if(verbose)		cat(paste("\nadding treatment data"))
		load(file.treatment)
		df		<- subset(df, select=c(Patient,AnyT_T1))
		df		<- df[,list(AnyT_T1= min(AnyT_T1)),by=Patient]
		df.all	<- merge(df.all, df, all.x=1, by="Patient")
		#add first RNA Virology date
		if(verbose)		cat(paste("\nadding virology data"))
		load(file.viro)
		tmp		<- subset(df, !is.na(PosRNA1), select=c(Patient,PosRNA1))
		tmp		<- tmp[,list(PosRNA1=min(PosRNA1)),by=Patient]
		df.all	<- merge(df.all, tmp, all.x=1, by="Patient")
		#add closest lRNA		
		df		<- subset(df,!is.na(lRNA), select=c(Patient,PosRNA,lRNA))
		df		<- merge( subset(df.all, select=c(FASTASampleCode,Patient,PosSeqT)), df, allow.cartesian=T, by="Patient" )
		setkey(df, FASTASampleCode)
		df		<- df[, { tmp<- which.min(difftime(PosSeqT, PosRNA, units="weeks")); list(PosRNA=PosRNA[tmp], lRNA=lRNA[tmp]) } , by=FASTASampleCode]		
		df.all	<- merge(df.all, df, all.x=1, by="FASTASampleCode")		
		#add first CD4 count date
		if(verbose)		cat(paste("\nadding CD4 data"))
		load(file.immu)
		tmp		<- subset(df, !is.na(PosCD41), select=c(Patient,PosCD41))
		tmp		<- tmp[,list(PosCD41=min(PosCD41)),by=Patient]
		df.all	<- merge(df.all, tmp, all.x=1, by="Patient")
		
		#get first positive event time
		if(verbose)		cat(paste("\ncompute first pos event"))
		#define first pos event
		tmp					<- df.all[, PosT]
		df.all[,AnyPos_T1:=tmp]
		tmp2				<- which( df.all[, is.na(AnyPos_T1) | (!is.na(AnyT_T1) &  AnyPos_T1>AnyT_T1) ] )
		if(verbose)	cat(paste("\nfound PosT>AnyT_T1, n=",length(tmp2)))
		set(df.all, tmp2, "AnyPos_T1", df.all[tmp2, AnyT_T1])		
		tmp2				<- which( df.all[, is.na(AnyPos_T1) | (!is.na(PosRNA1) &  AnyPos_T1>PosRNA1) ] )
		if(verbose)	cat(paste("\nfound AnyPos_T1>PosRNA1, n=",length(tmp2)))				
		set(df.all, tmp2, "AnyPos_T1", df.all[tmp2, PosRNA1])		
		tmp2				<- which( df.all[, is.na(AnyPos_T1) | (!is.na(PosCD41) &  AnyPos_T1>PosCD41) ] )
		if(verbose)	cat(paste("\nfound AnyPos_T1>PosCD41, n=",length(tmp2)))				
		set(df.all, tmp2, "AnyPos_T1", df.all[tmp2, PosCD41])
		
		tmp2				<- which( df.all[, is.na(AnyPos_T1) | (!is.na(PosCD41) &  AnyPos_T1>PosSeqT) ] )
		if(verbose)	cat(paste("\nfound AnyPos_T1>PosCD41, n=",length(tmp2)))				
		set(df.all, tmp2, "AnyPos_T1", df.all[tmp2, PosSeqT])
		
		
		if(	length(which(df.all[, AnyPos_T1>PosT])) 		||
				length(which(df.all[, AnyPos_T1>PosSeqT])) 	|| 
				length(which(df.all[, AnyPos_T1>PosRNA1]))	||
				length(which(df.all[, AnyPos_T1>PosCD41]))	||
				length(which(df.all[, AnyPos_T1>AnyT_T1]))		) 	stop("something is wrong with PosEarliest")			
		
		if(verbose)	cat(paste("\nsave to file",file.out))
		setkey(df.all, FASTASampleCode)
		save(df.all,file=file.out)
		str(df.all)		
	}
	stop()
	#compute file All1stPatientCovariates
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
		df		<- df[,{tmp<- which.min(PosRNA); list(PosRNA1= PosRNA[tmp], RNA1= RNA[tmp]) }, by=Patient]		
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
		df		<- df[,{tmp<- which.min(PosCD4); list(PosCD41= PosCD4[tmp], CD41= CD4A[tmp]) }, by=Patient]
		if(verbose)		cat(paste("\npatients with at least one non-missing CD4 date", nrow(df)))
		if(0)
		{
			print(range(df[,PosCD41], na.rm=1))		
		}
		df.all<- df[df.all]
		df.all<- subset(df.all,select=c(Patient,Trm,SeqPROT,SeqRT,PosPROT,PosRT,PosT,PosT_Acc,PosCD41,PosRNA1, NegT,NegT_Acc, AnyTherT1, isAcute,isDead,Died, RNA1, CD41))
		
		
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
		df.all<- subset(df.all, select=c(Patient,Trm,SeqPROT,SeqRT,PosSeq,PosT,PosT_Acc,PosCD41,PosRNA1,NegT,NegT_Acc, AnyTherT1, isAcute,isDead,Died,RNA1, CD41))
		
		
		if(verbose)		cat(paste("\nsave df.all to", file.out))				
		save(df.all,file=file.out)
		str(df.all)
	}
	
	
	
	
	stop()
	#compute PosMeta
	tmp<- df.all[, min(PosT, PosCD4, PosRNA), by=Patient]$V1
	df.all[,PosMeta:=tmp]
	print(df.all)
	stop()

	tmp<- subset(df.all, !is.na(PosSeq) & PosT>PosSeq)
	print(tmp)
	stop()
	print(df.all)	
	stop()	
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
		stop()
	}
	if(0)
	{
		NA.time			<- c("01/01/1911","01/11/1911","11/11/1911")	
		MAX.time		<- c("")
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
				df[nok.idx,x]	<- "01/01/2999"
			df[,x]			<- as.Date(df[,x], format="%d/%m/%Y")	
		}
		
		df		<- data.table(df, key="Patient")
		setnames(df, "T0","HAART_T1")
		set(df,NULL,"Patient",factor(df[,Patient]))
		tmp		<- df[,StartTime]
		df[,"AnyT_T1":=tmp]
		tmp		<- which(is.na(df[,AnyT_T1]))
		set(df, tmp, "AnyT_T1", df[tmp,HAART_T1])		
		tmp		<- df[, { tmp<- which.min(StartTime); list(HAART_T1=HAART_T1[tmp], AnyT_T1=StartTime[tmp])}, by=Patient]
		if( nrow(subset(tmp,!is.na(HAART_T1) & HAART_T1<AnyT_T1)) )	stop("found AnyT_T1 that is older than HAART_T1")
		tmp		<- subset(tmp,select=c(Patient,AnyT_T1))
		df		<- merge(df,tmp,all.x=1,by="Patient")
		setnames(df, "AnyT_T1.y","AnyT_T1")
		df		<- subset(df, select= 1-ncol(df) )

		file		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
		if(verbose) cat(paste("\nsave to", file))
		save(df, file=file)		
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
		NA.time			<- c("","01/01/1911","11/11/1911","24/06/1923")		
		#read VIROLOGY csv data file and preprocess
		file			<- paste(dir.name,"derived/ATHENA_2013_03_Viro.csv",sep='/')
		df				<- read.csv(file, stringsAsFactors=FALSE)									
		df$Undetectable	<- factor(df$Undetectable, levels=c(0,1),labels=c("No","Yes"))
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
		setnames(df, "DateRNA","PosRNA")
		set(df, NULL, "Patient", factor(df[,Patient]))				
		tmp		<- round(log10( df[,RNA] ), d=3) 
		df[,"lRNA":=tmp]
		tmp		<- which(df[, Undetectable=="No"])
		set(df, tmp, "lRNA", NA)		
		tmp		<- df[,list(PosRNA1=min(PosRNA)),by=Patient]
		df		<- merge(df, tmp, all.x=1, by="Patient")
		
		file		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
		if(verbose) cat(paste("\nsave to", file))
		save(df, file=file)		
	}
	if(0)
	{
		NA.time			<- c("","01/01/1911","11/11/1911","24/06/1923")		
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
		tmp		<- df[,list(PosCD41=min(PosCD4)),by=Patient]
		df		<- merge(df, tmp, all.x=1, by="Patient")
		
		file		<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
		if(verbose) cat(paste("\nsave to", file))
		save(df, file=file)		
	}
}
######################################################################################
#create PROT+RT data set of first sequences from all patients
hivc.prog.get.bootstrapseq<- function(check.any.bs.identical=0)
{	
	library(ape)
	library(data.table)
	library(hivclust)
	
	indir		<- outdir		<- paste(DATA,"tmp",sep='/')
	infile		<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
	signat.out	<- signat.in	<- "Sat_May_11_14/23/46_2013"
	verbose		<- resume		<- 1
	bs			<- 0
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.in<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.out<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									bootstrap= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) bs<- tmp[1]
	}
	if(1)
	{
		print(indir)
		print(outdir)
		print(infile)
		print(signat.in)
		print(signat.out)
		print(verbose)
		print(resume)
		print(bs)
		print(signat.in)
		print(signat.out)
	}
	pattern 	<- paste(infile,"_",gsub('/',':',signat.out),".phylip.",sprintf("%03d",bs),sep='')
	file		<- list.files(path=outdir, pattern=pattern, full.names=1)
	if(!resume || !length(file))	
	{					
		file		<- paste(outdir,"/",infile,"_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nread",file))
		load(file)
		print(seq.PROT.RT)
		#print(bs)
		if(bs)		#keep bs0 intact
		{
			dummy			<- 0
			any.eq			<- 1
			j				<- 0
			while(any.eq)
			{
				j			<- j+1
				bs.blocks.n	<- floor( ncol(seq.PROT.RT )/3)
				bs.blocks.s	<- sample(seq_len(bs.blocks.n),bs.blocks.n,replace=1)-1
				bs.seq.s	<- as.numeric( sapply(bs.blocks.s,function(x)		3*x+c(1,2,3)		) )
				seq.BS		<- seq.PROT.RT[,bs.seq.s]
				if(check.any.bs.identical)
				{
					if(verbose) cat(paste("\ncheck for identity proposed boostrap seq alignment no",j))
					#check no seqs identical								
					for(i1 in seq_len(nrow(seq.BS)-1))
					{
						seq1		<- seq.BS[i1,]
						tmp			<- 1-sapply(seq.int(i1+1,nrow(seq.BS)),function(i2)
													{		
														.C("hivc_dist_ambiguous_dna", seq1, seq.BS[i2,], ncol(seq1), dummy )[[4]]			
													})
						#print(tmp)
						if(any(tmp==0))
						{
							print(tmp)
							break
						}									
						if(i1==nrow(seq.BS)-1)
							any.eq	<- 0
					}
					if(verbose) cat(paste("\nchecked for identity proposed boostrap seq alignment no",j,"is any identical",any.eq))
				}
				else
					any.eq	<- 0
			}					
		}
		else
		{
			cat(paste("\nkeep boostrap seq alignment no",bs,"as original"))
			seq.BS	<- seq.PROT.RT
		}
		file		<- paste(outdir,"/",infile,"_",gsub('/',':',signat.out),".phylip.",sprintf("%03d",bs),sep='')
		cat(paste("\nsave boostrap seq alignment to",file))
		hivc.seq.write.dna.phylip(seq.BS, file=file)
	}
	else
		cat("\nfound boostrap sequence alignment")
}
######################################################################################
#create PROT+RT data set of first sequences from all patients
hivc.prog.get.allseq<- function()
{	
	library(ape)
	library(data.table)
	library(RFLPtools)
	library(hivclust)
	
	indir		<- paste(DATA,"derived",sep='/')
	outdir		<- paste(DATA,"tmp",sep='/')
	infile		<- "ATHENA_2013_03_All1stPatientCovariates.R"
	signat.out	<- "Sat_Jun_16_17/23/46_2013"
	signat.in	<- "Wed_May__1_17/08/15_2013"
	verbose		<- 1 
	resume		<- 0
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.in<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.out<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
	}
	if(1)
	{
		print(indir)
		print(outdir)
		print(infile)
		print(signat.in)
		print(signat.out)
		print(verbose)
		print(resume)
	}
	pattern 	<- gsub('/',':',paste("ATHENA_2013_03_All+LANL_Sequences_",signat.out,".R$",sep=''))
	file		<- list.files(path=outdir, pattern=pattern, full.names=1)
	if(!resume || !length(file))	
	{						
		#create file of sequences that are in preliminary cluster. use these as seed to enrich data set to get more reliable boostrap
		dir.name		<- DATA
		infile.tree		<- "ATHENA_2013_03_FirstCurSequences_PROTRT_examlbs100"
		infile.seq		<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
		signat.in		<- "Sat_May_11_14/23/46_2013"
		
		
		infile.enrich	<- "LosAlamos_HIV1B_Prot_P51_n4243_gapsyes"
		outfile.enrich	<- "LosAlamos_HIV1B_Prot_P51_n4243_gapsno"
		infile.enrich	<- "LosAlamos_HIV1B_Prot_n85928_gapsyes"
		outfile.enrich	<- "LosAlamos_HIV1B_Prot_n85928_gapsno"

		if(0) 
		{
			insignat				<- "Sat_Jun_15_18/23/46_2013"
			#build BLAST database
			#remove all gaps from FASTA file
			file<- paste(paste(dir.name,"original",sep='/'),'/',paste(infile.enrich,'_',gsub('/',':',insignat),".fasta", sep=''), sep='')						
			seq.enrich				<- read.dna( file, format="fa" )
			seq.enrich				<- hivc.seq.rmgaps(seq.enrich, rm.only.col.gaps=0)
			names(seq.enrich)		<- gsub('-','NA',names(seq.enrich))					
			file<- paste(paste(dir.name,"derived",sep='/'),'/',paste(outfile.enrich,'_',gsub('/',':',insignat),".fasta", sep=''), sep='')			
			write.dna(seq.enrich, file=file, format="fasta", append=0, colsep='', colw=80, blocksep=0)
			#build BLAST database
			cmd						<- hivc.cmd.blast.makedb(paste(dir.name,"derived",sep='/'), outfile.enrich, signat=gsub('/',':',insignat), with.mask=0, verbose=0)
			cat ( cmd )
			#cat( system(cmd, intern=TRUE) )
		}
		if(0)
		{	
			plot							<- 0
			#identify query sequences for BLAST search from large trial clusters
			signat.in						<- "Fri_May_24_12/59/06_2013"						
			#read tree to construct large trial clusters
			file							<- paste(dir.name,"tmp",paste(infile.tree,'_',gsub('/',':',signat.in),".newick",sep=''),sep='/')
			cat(paste("read file",file))
			ph								<- ladderize( read.tree(file) )
			#read bootstrap support values		
			ph.node.bs						<- as.numeric( ph$node.label )
			ph.node.bs[is.na(ph.node.bs)]	<- 0
			ph.node.bs						<- ph.node.bs/100
			ph$node.label					<- ph.node.bs
			dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf")			
			#produce trial clustering
			thresh.bs						<- 0.9
			thresh.brl						<- 0.105		#subst rate 3.4 10^-3  so for 15 years either way have 10.5 * 10^-2
			clustering						<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
			#get tip names in trial cluster
			ph.tips.in.cluster				<- which( !is.na(clustering[["clu.mem"]][seq_len(Ntip(ph))]) )
			ph.tips.in.cluster.names		<- ph$tip.label[ph.tips.in.cluster]
			#get sequences in trial cluster			
			file							<- paste(dir.name,"tmp",paste(infile.seq,'_',gsub('/',':',signat.in),".R",sep=''),sep='/')
			load(file)
			if(verbose) cat(paste("\nload trial sequences from",file,sep=''))
			seq.Trial.PROT.RT				<- seq.PROT.RT[ph.tips.in.cluster.names,]
			seq.Trial.PROT.RT				<- hivc.seq.replace(seq.Trial.PROT.RT, code.from='?', code.to='n')			
			if(verbose) cat(paste("\nfound trial sequences, n=",nrow(seq.Trial.PROT.RT),sep=''))			
			outsignat						<- "Sat_Jun_15_18/23/46_2013"
			file							<- paste(paste(dir.name,"tmp",sep='/'),'/',paste(infile.tree,"_intrialcluster_",gsub('/',':',outsignat),".fasta", sep=''), sep='')
			if(verbose) cat(paste("\nwrite trial sequences to",file,sep=''))
			write.dna(seq.Trial.PROT.RT, file=file, format="fasta", append=0, colsep='', colw=ncol(seq.Trial.PROT.RT), blocksep=0)
			
			#BLAST search against PROT+P51 database
			indir							<- paste(dir.name,"tmp",sep='/')
			infile							<- paste(infile.tree,"_intrialcluster",sep='')
			insignat						<- "Sat_Jun_15_18/23/46_2013"
			dbdir							<- paste(dir.name,"derived",sep='/')
			dbfile							<- "LosAlamos_HIV1B_Prot_P51_n4243_gapsno"			
			outsignat						<- "Sat_Jun_15_18/23/46_2013"
			cmd<- hivc.cmd.blast(indir, infile, gsub('/',':',insignat), dbdir, dbfile, gsub('/',':',insignat), outdir=indir, outsignat=gsub('/',':',outsignat), blast.max_target_seqs=20)			
			cat(cmd)
			system(cmd, intern=TRUE)	#wait until finished	
			#read BLAST file against PROT+P51 database and take unique best 10 hits
			file							<- paste(paste(dir.name,"tmp",sep='/'),'/',paste(infile.tree,"_intrialcluster_",gsub('/',':',signat.out),".blast", sep=''), sep='')
			seq.enrich						<- hivc.seq.blast.read(file=file)
			if(plot)
			{
				tmp							<- 1-seq.enrich[,identity]/100
				summary(tmp)
				hist(tmp)
			}
			seq.enrich.unique				<- data.table(subject.id=unique(seq.enrich[,subject.id]), key="subject.id")			
			tmp								<- seq.enrich.unique[, substr(subject.id,3,4)]
			tmp[tmp=="NA"]					<- NA
			seq.enrich.unique[,subject.location:=tmp]
			tmp								<- sapply(seq.enrich.unique[, strsplit(subject.id,'.',fixed=1) ], function(x) x[3])
			tmp[tmp=="NA"]					<- NA
			seq.enrich.unique[,subject.PosSeq:=as.numeric(tmp)]
			if(plot)
			{
				hist(seq.enrich.unique[,subject.PosSeq])
				seq.enrich.unique.loc		<- sort(table(seq.enrich.unique[,subject.location]),decreasing=TRUE)
				barplot(seq.enrich.unique.loc, cex.names=.75)
			}			
			#best PROT+p51 hits to enrich NL data set
			indir				<- paste(dir.name,"derived",sep='/')
			infile				<- "LosAlamos_HIV1B_Prot_P51_n4243_gapsno"
			insignat			<- "Sat_Jun_15_18/23/46_2013"
			file				<- paste(indir,'/',paste(infile,'_',gsub('/',':',insignat),".fasta", sep=''), sep='')
			seq.enrich.PROT.P51	<- read.dna( file, format="fa" )
			tmp					<- which( names(seq.enrich.PROT.P51)%in%seq.enrich.unique[,subject.id] )
			tmp2				<- lapply( tmp, function(i)	seq.enrich.PROT.P51[[i]] )
			names(tmp2)			<- paste("PROT+P51",names(seq.enrich.PROT.P51)[tmp],sep="_")
			class(tmp2)			<- "DNAbin"
			seq.enrich.PROT.P51	<- tmp2
			
			
			#BLAST search against PROT database
			indir							<- paste(dir.name,"tmp",sep='/')
			infile							<- paste(infile.tree,"_intrialcluster",sep='')			
			dbdir							<- paste(dir.name,"derived",sep='/')			
			dbfile							<- "LosAlamos_HIV1B_Prot_n85928_gapsno"
			outsignat						<- "Sat_Jun_16_17/23/46_2013"
			cmd<- hivc.cmd.blast(indir, infile, gsub('/',':',insignat), dbdir, dbfile, gsub('/',':',insignat), outdir=indir, outsignat=gsub('/',':',outsignat), blast.max_target_seqs=20)			
			cat(cmd)
			#cat( system(cmd, intern=TRUE) )
			
			#read BLAST file against PROT database and take unique best 10 hits that are geographically distant as TRUE NEGATIVES
			file							<- paste(paste(dir.name,"tmp",sep='/'),'/',paste(infile.tree,"_intrialcluster_",gsub('/',':',outsignat),".blast", sep=''), sep='')
			seq.enrich						<- hivc.seq.blast.read(file=file)
			if(plot)
			{
				tmp							<- 1-seq.enrich[,identity]/100
				summary(tmp)
				hist(tmp)
			}
			seq.enrich.unique				<- data.table(subject.id=unique(seq.enrich[,subject.id]), key="subject.id")			
			tmp								<- seq.enrich.unique[, substr(subject.id,3,4)]
			tmp[tmp=="NA"]					<- NA
			seq.enrich.unique[,subject.location:=tmp]
			tmp								<- sapply(seq.enrich.unique[, strsplit(subject.id,'.',fixed=1) ], function(x) x[3])
			tmp[tmp=="NA"]					<- NA
			seq.enrich.unique[,subject.PosSeq:=as.numeric(tmp)]
			if(plot)
			{
				hist(seq.enrich.unique[,subject.PosSeq])
				seq.enrich.unique.loc			<- sort(table(seq.enrich.unique[,subject.location]),decreasing=TRUE)
				barplot(seq.enrich.unique.loc, cex.names=.75)
			}
			
			#to build HIVC.COUNTRY.TABLE, use this to search LANL
			#tmp								<- seq.enrich.unique[,min(subject.id),by=subject.location]
			#tmp								<- sapply(tmp[, strsplit(V1,'.',fixed=1) ], function(x) x[5])			
			#extract seq.enrich.unique that are very unlikely to transmit to NL		seq.enrich.unique.loc.excl taken from Bezemer2013
			seq.enrich.unique.loc.level		<- "0.01"			
			seq.enrich.unique.loc.excl		<- structure(list(	`0.01` = c("GREENLAND", "CHILE", "ALBANIA", "Suriname", "SINGAPORE", "PANAMA", "SENEGAL", "FINLAND", "TAIWAN", "JAPAN", 
																			"SLOVAKIA", "LUXEMBOURG", "AUSTRALIA", "THAILAND", "SERBIA", "ROMANIA", "NORWAY", "VENEZUELA", "SLOVENIA", "RUSSIAN FEDERATION", 
																			"AUSTRIA", "SWEDEN", "FRANCE", "SOUTH KOREA", "CYPRUS", "MONTENEGRO", "CUBA", "HONDURAS", "DENMARK", "POLAND", "PORTUGAL", "CHINA", 
																			"SWITZERLAND", "BRAZIL", "ARGENTINA", "GERMANY", "BELGIUM", "CANADA","CZECH REPUBLIC", "UNITED STATES", "SPAIN", "ITALY", "UNITED KINGDOM","NETHERLANDS"), 
																`0.05` = c("SERBIA", "ROMANIA", "NORWAY", "VENEZUELA", "SLOVENIA","RUSSIAN FEDERATION", "AUSTRIA", "SWEDEN", "FRANCE", "SOUTH KOREA", 
																			"CYPRUS", "MONTENEGRO", "CUBA", "HONDURAS", "DENMARK", "POLAND","PORTUGAL", "CHINA", "SWITZERLAND", "BRAZIL", "ARGENTINA", "GERMANY",
																			"BELGIUM", "CANADA", "CZECH REPUBLIC", "UNITED STATES", "SPAIN","ITALY", "UNITED KINGDOM","NETHERLANDS")), .Names = c("0.01", "0.05"))																										
			seq.enrich.unique.loc.incl		<- subset(HIVC.COUNTRY.TABLE, !country%in%seq.enrich.unique.loc.excl[[seq.enrich.unique.loc.level]] )
			seq.enrich.unique.loc.incl		<- subset(seq.enrich.unique.loc.incl, !country%in%c("BULGARIA","GRENADA","GIORGIA","GREECE","HUNGARY","UKRAINE","LATVIA") )						
			seq.enrich.unlikelytransmission	<- subset(seq.enrich.unique,subject.location%in%as.character(seq.enrich.unique.loc.incl[,key]) )
			if(plot)
			{				
				hist(seq.enrich.unlikelytransmission[,subject.PosSeq])
				seq.enrich.unique.loc			<- sort(table(seq.enrich.unlikelytransmission[,subject.location]),decreasing=TRUE)
				barplot(seq.enrich.unique.loc, cex.names=.75)
			}
			#best geographically distant PROT hits to enrich NL data set
			indir				<- paste(dir.name,"derived",sep='/')
			infile				<- "LosAlamos_HIV1B_Prot_n85928_gapsno"
			insignat			<- "Sat_Jun_15_18/23/46_2013"
			file				<- paste(indir,'/',paste(infile,'_',gsub('/',':',insignat),".fasta", sep=''), sep='')
			seq.enrich.TN		<- read.dna( file, format="fa" )
			tmp					<- which( names(seq.enrich.TN)%in%seq.enrich.unlikelytransmission[,subject.id] )
			tmp2				<- lapply( tmp, function(i)	seq.enrich.TN[[i]] )
			names(tmp2)			<- paste("TN",names(seq.enrich.TN)[tmp],sep="_")
			class(tmp2)			<- "DNAbin"
			seq.enrich.TN		<- tmp2

			#combine all sequence data sets and save
			seq.enrich			<- c(seq.enrich.PROT.P51,seq.enrich.TN)
			outdir				<- paste(dir.name,"derived",sep='/')
			outsignat			<- "Sat_Jun_16_17/23/46_2013"			
			file				<- paste(outdir,"/LosAlamos_EnrichSequences_For_ATHENA201303_",gsub('/',':',outsignat),".R",sep='')
			if(verbose) cat(paste("\nwrite to",file))
			save(seq.enrich, file=file)						
		}
		if(0)
		{
			#get all ATHENA sequences into PROT+RT format, and add foreign sequences
			indir		<- paste(dir.name,"tmp",sep='/')
			outdir		<- paste(dir.name,"tmp",sep='/')
			insignat	<- "Wed_May__1_17/08/15_2013"
			outsignat	<- "Sat_Jun_16_17/23/46_2013"
			outfile		<- "ATHENA_2013_03_All+LANL_Sequences"
			
			if(verbose)	cat(paste("\ncreate LosAlamos_EnrichSequences_For_ATHENA201303 file"))						
			#create matrix of all PROT+RT sequences				
			pattern 	<- gsub('/',':',paste(insignat,".clustalo$",sep=''))
			files		<- list.files(path=indir, pattern=pattern, full.names=1)
			#read all sequences and add a missing one with name "NA" if not in union of all possible sequences take at a sampling date		
			seq.PROT	<- read.dna( files[ grep("PROT",files) ], format="fa", as.matrix=1 )
			if(verbose) cat(paste("\nfound PROT sequences, n=",nrow(seq.PROT),"\n"))
			seq.RT		<- read.dna( files[ grep("RT",files) ], format="fa", as.matrix=1 )
			if(verbose) cat(paste("\nfound RT sequences, n=",nrow(seq.RT),"\n"))
			seq.nam.all	<- union( rownames(seq.PROT), rownames(seq.RT) )
			if(verbose) cat(paste("\nfound unique sampling IDs, n=",length(seq.nam.all),"\n"))
			seq.nam.PROT<- seq.nam.all
			seq.nam.RT	<- seq.nam.all
			seq.nam.PROT[ !seq.nam.all%in%rownames(seq.PROT) ]	<- "NA"		#those sampling dates missing among PROT get name NA
			seq.nam.RT[ !seq.nam.all%in%rownames(seq.RT) ]		<- "NA"
			#prepare PROT and RT sequences: add a missing one with name "NA" 		
			seq.PROT	<- read.dna( files[ grep("PROT",files) ], format="fa", as.matrix=1 )
			tmp			<- as.DNAbin( matrix(rep('-',ncol(seq.PROT)),1,ncol(seq.PROT), dimnames=list(c("NA"),c())) )
			seq.PROT	<- rbind(seq.PROT,tmp)
			seq.RT		<- read.dna( files[ grep("RT",files) ], format="fa", as.matrix=1 )						 				
			tmp			<- as.DNAbin( matrix(rep('-',ncol(seq.RT)),1,ncol(seq.RT), dimnames=list(c("NA"),c())) )
			seq.RT		<- rbind(seq.RT,tmp)
			#add missing sequence where needed to combine PROT and RT 
			seq.PROT	<- seq.PROT[seq.nam.PROT,]
			seq.RT		<- seq.RT[seq.nam.RT,]
			rownames(seq.PROT)	<- seq.nam.all
			rownames(seq.RT)	<- seq.nam.all
			#combine PROT and RT
			seq.PROT.RT	<- cbind(seq.PROT,seq.RT)
			if(verbose) cat(paste("\ncombined PROT and RT sequences, n=",nrow(seq.PROT.RT),"\n"))
			print(seq.PROT.RT)
			
			#load EnrichSequences
			indir				<- paste(dir.name,"derived",sep='/')
			insignat			<- "Sat_Jun_16_17/23/46_2013"
			file				<- paste(indir,"/LosAlamos_EnrichSequences_For_ATHENA201303_",gsub('/',':',insignat),".R",sep='')
			if(verbose) cat(paste("\nloading file",file))
			load(file)	
			print(seq.enrich)
			#load HXB2 reference sequence
			data( refseq_hiv1_hxb2 )
			hxb2				<- as.character( data.table( hxb2 )[, HXB2.K03455 ] )
			hxb2				<- hxb2[seq_len(length(hxb2)-2)]
			
			#write all sequences to fasta file
			file		<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".fasta",sep='')			
			if(verbose) cat(paste("\nwrite reference sequence to",file,"\n"))
			cat( paste(">HXB2\n",paste( hxb2, collapse='',sep='' ),"\n",sep=''), file=file, append=0)			
			if(verbose) cat(paste("\nappend ATHENA combined PROT and RT sequences to",file,"\n"))
			write.dna(seq.PROT.RT, file=file, format="fasta", append=1, colsep='', colw=length(hxb2), blocksep=0)
			if(verbose) cat(paste("\nappend LosAlamos_EnrichSequences to",file,"\n"))
			write.dna(seq.enrich, file=file, format="fasta", append=1, colsep='', colw=length(hxb2), blocksep=0)				 
		}		
		if(1)	#curate alignment with reference 
		{
			indir								<- paste(dir.name,"tmp",sep='/')
			infile								<- "ATHENA_2013_03_All+LANL_Sequences"
			outdir								<- paste(dir.name,"tmp",sep='/')
			outfile								<- "ATHENA_2013_03_CurAll+LANL_Sequences"
			insignat							<- "Sat_Jun_16_17/23/46_2013"
			outsignat							<- "Sat_Jun_16_17/23/46_2013"
			file								<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo",sep='')			
			if(verbose) cat(paste("\nread ",file))
			seq.PROT.RT							<- read.dna(file, format="fasta", as.character=1)
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2273, c("g","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2273:2275]<- matrix( c("-","-","g"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2348, c("?","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2348:2349]<- matrix( c("g","?"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2350, c("-","-","-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2350:2353]<- matrix( c("t","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			#always delete drug resistance from pos 2356 because this is really cryptic; leave g at 2355 always OK
			seq.PROT.RT[,2356:2363]				<- matrix( c("-","-","-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=8, byrow=1 )			
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo1",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)			
			#manual edits seq 3151 3243 9049
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2624, c("t","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2624:2625]<- matrix( c("-","t"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2634, c("-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2634:2635]<- matrix( c("g","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			seq.PROT.RT.sort.by					<- apply(seq.PROT.RT,1,function(x)		which(rev(x)!="-" )[1]  )
			seq.PROT.RT							<- seq.PROT.RT[sort(seq.PROT.RT.sort.by, index.return=1)$ix,]
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo2",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2759, c("-","g","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("g","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2759, c("-","a","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("a","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2760, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2760, c("-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("g","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2760, c("c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("-","c"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2761, c("-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2761:2762]<- matrix( c("c","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			#always delete drug resistance from pos 2752 because this is really cryptic; leave a at 2751 always OK
			seq.PROT.RT[,2752:2759]				<- matrix( c("-","-","-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=8, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo3",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			#double fixup needed
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2759, c("-","g","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("g","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2759, c("-","g","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("g","g","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2759, c("-","a","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("a","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2759, c("-","r","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("r","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			seq.PROT.RT[,2752:2760]				<- matrix( c("-","-","-","-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=9, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2760, c("-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("g","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 2760, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			seq.PROT.RT[,2752:2760]				<- matrix( c("-","-","-","-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=9, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo4",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3149, c("-"))
			if(length(query.yes))
			{	
				seq.PROT.RT[query.yes,3149:3156]<- seq.PROT.RT[query.yes,3150:3157]
				seq.PROT.RT[query.yes,3157]		<- "-"
			}	
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo5",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)				
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3156, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3156:3157]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3157, c("a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3157:3158]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
			#rm col 3213 - not in HB2 and we keep the standard numbering
			seq.PROT.RT							<- seq.PROT.RT[,-3213]
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo6",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3243, c("-"))
			if(length(query.yes))
			{	
				seq.PROT.RT[query.yes,3241:3243]<- seq.PROT.RT[query.yes,3240:3242]
				seq.PROT.RT[query.yes,3240]		<- "-"
			}	
			#fixup: rm col 2671 to 2676 - not in HB2 and we return to the standard numbering
			seq.PROT.RT							<- seq.PROT.RT[,c(1:2670,2677:ncol(seq.PROT.RT))]
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo7",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)						
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3433, c("-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3433:3434]<- matrix( c("t","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3433, c("-","-","-","-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3433:3437]<- matrix( c("t","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3434, c("-","-","-","-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3434:3438]<- matrix( c("a","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3434, c("g","-","-","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3434:3438]<- matrix( c("-","-","-","-","g"), nrow=length(query.yes), ncol=5, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo8",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3434, c("-","-","-","-","g","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3434:3439]<- matrix( c("g","-","-","-","-","g"), nrow=length(query.yes), ncol=6, byrow=1 )						
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3435, c("-","-","-","-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3435:3439]<- matrix( c("a","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3435, c("-","-","-","-","r"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3435:3439]<- matrix( c("r","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3436, c("-","-","-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3439]<- matrix( c("c","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3436, c("-","-","-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3439]<- matrix( c("t","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3436, c("-","-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3438]<- matrix( c("c","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3436, c("-","-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3438]<- matrix( c("t","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3436, c("-","-","y"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3438]<- matrix( c("y","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3437, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3437:3438]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3438, c("-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3438:3439]<- matrix( c("g","-"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3438, c("a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3437:3438]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
			
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo9",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3491, c("c","a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3491:3493]<- matrix( c("-","c","a"), nrow=length(query.yes), ncol=3, byrow=1 )															
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3501, c("-","t","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3501:3503]<- matrix( c("t","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3501, c("a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3501:3502]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3503, c("-"))
			for(i in query.yes)
				seq.PROT.RT[i,3501:3503]		<- c("-",seq.PROT.RT[i,3501:3502])			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3504, c("-","-","-","t"))	
			for(i in query.yes)
				seq.PROT.RT[i,3504:3524]		<- c(seq.PROT.RT[i,3507:3524],"-","-","-")
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3506, c("-","-","-","t","g"))
			for(i in query.yes)
				seq.PROT.RT[i,3506:3526]		<- c(seq.PROT.RT[i,3509:3526],"-","-","-")
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3507, c("-","-","-","g","a"))
			for(i in query.yes)
				seq.PROT.RT[i,3507:3537]		<-	c(seq.PROT.RT[i,3510:3537],"-","-","-")			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3543, c("-","c"))
			if(length(query.yes))
			{	
				seq.PROT.RT[query.yes,3543:3550]<- seq.PROT.RT[query.yes,3544:3551]
				seq.PROT.RT[query.yes,3551]		<- "-"
			}
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("-","c","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("-","g","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("g","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("-","g","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("g","t","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","g","c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3551:3554]<- matrix( c("-","g","g","c"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3495, c("g","-","-","-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3495:3499]<- matrix( c("g","g","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo10",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3497, c("-","-","-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3497:3500]<- matrix( c("a","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )						
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3497, c("-","-","-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3497:3500]<- matrix( c("g","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )									
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3499, c("t","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3499:3501]<- matrix( c("-","-","t"), nrow=length(query.yes), ncol=3, byrow=1 )			
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3578, c("g","a","g","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3578:3581]<- matrix( c("-","g","a","g"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3578, c("g","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3578:3579]<- matrix( c("-","g"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","g"))
			for(i in query.yes)
				seq.PROT.RT[i,3551:3600]		<-	c("-",seq.PROT.RT[i,3551:3599])
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3552, c("g","c"))
			for(i in query.yes)
				seq.PROT.RT[i,3552:3600]		<-	c("-",seq.PROT.RT[i,3552:3599])
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3498, c("c","g","t"))
			for(i in query.yes)
				seq.PROT.RT[i,3498:3520]		<-	c("-",seq.PROT.RT[i,3498:3519])
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("a","a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","a","c"))
			for(i in query.yes)
				seq.PROT.RT[i,3551:3560]		<-	c("-",seq.PROT.RT[i,3551:3559])
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","a","y","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3551:3554]<- matrix( c("-","g","a","y"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("-","c","c","g","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","-","c"), nrow=length(query.yes), ncol=3, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo11",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			
			#manual curation on the end
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3422, c("t","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3422:3424]<- matrix( c("-","-","t"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3494, c("g","g","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3494:3496]<- matrix( c("-","g","g"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3498, c("-","-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3498:3500]<- matrix( c("g","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3578, c("g","a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3578:3580]<- matrix( c("-","g","a"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3618, c("g","-","-","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3618:3622]<- matrix( c("-","-","-","-","g"), nrow=length(query.yes), ncol=5, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3499, c("-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3499:3500]<- matrix( c("c","-"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3504, c("a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3504:3505]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3437, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3437:3438]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3438, c("c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3438:3439]<- matrix( c("-","c"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 3444, c("c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3444:3445]<- matrix( c("-","c"), nrow=length(query.yes), ncol=2, byrow=1 )
			#cut before PROT pol start and at max len in database
			seq.PROT.RT							<- seq.PROT.RT[,2253:ncol(seq.PROT.RT)] 
			seq.PROT.RT							<- seq.PROT.RT[,1:1624]
			seq.PROT.RT.sort.by					<- apply(seq.PROT.RT,1,function(x)		which(rev(x)!="-" )[1]  )
			seq.PROT.RT							<- seq.PROT.RT[sort(seq.PROT.RT.sort.by, index.return=1)$ix,]
			#always delete drug resistance at G333D/E because this is at the end of the alignment and alignment unreliable - could have picked up other ends
			seq.PROT.RT[,1294:1299]				<- matrix( c("-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=6, byrow=1 )			
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo12",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","g","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1302]<- matrix( c("g","g","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","g","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1302]<- matrix( c("g","g","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","g","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1302]<- matrix( c("g","g","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","a","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1302]<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","g","-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("g","g","g","c"), nrow=length(query.yes), ncol=4, byrow=1 )			
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo14",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)						
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("a","a","c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("-","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","a","t","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("g","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("t","g","t","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("t","g","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("a","t","c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("-","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("-","g","t","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("-","g","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- hivc.seq.find(seq.PROT.RT, 1300, c("g","a","y","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("g","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			
			seq.PROT.RT							<- as.DNAbin( seq.PROT.RT )
			seq.PROT.RT							<- hivc.seq.replace(seq.PROT.RT, code.from='?', code.to='n')
			seq.PROT.RT							<- seq.PROT.RT[-1,]												#remove HXB2
			rownames(seq.PROT.RT)				<- gsub(' ','',rownames(seq.PROT.RT))
			
			file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".phylip",sep='')
			if(verbose)	cat(paste("\nwrite phylip file to",file))
			hivc.seq.write.dna.phylip(seq.PROT.RT, file=file)			
			file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
			if(verbose)	cat(paste("\nwrite R file to",file))
			save(seq.PROT.RT, file=file)
			file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".fasta",sep='')			
			if(verbose)	cat(paste("\nwrite fasta file to",file))
			write.dna(seq.PROT.RT, file=file, format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)						
		}
		stop()				
	}
	else
	{
		if(verbose)	cat(paste("\nload",file))
		load(file)
	}	
}
######################################################################################
#create PROT+RT data set of first sequences from all patients
hivc.prog.get.firstseq<- function()
{	
	library(ape)
	library(data.table)
	library(hivclust)
	
	indir		<- outdir		<- paste(DATA,"tmp",sep='/')
	infile		<- "ATHENA_2013_03_SeqMaster.R"
	signat.out	<- signat.in	<- "Wed_May__1_17/08/15_2013"
	verbose		<- resume		<- 1
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.in<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.out<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
	}
	if(1)
	{
		print(indir)
		print(outdir)
		print(infile)
		print(signat.in)
		print(signat.out)
		print(verbose)
		print(resume)
	}
	pattern 	<- gsub('/',':',paste("FirstAliSequences_PROTRT_",signat.in,".R$",sep=''))
	file		<- list.files(path=outdir, pattern=pattern, full.names=1)
	if(!resume || !length(file))	
	{		
		if(verbose)	cat(paste("\ncreate FirstAliSequences_PROTRT file"))			
		file		<- paste(indir,infile,sep='/')
		if(verbose)	cat(paste("\nload",file,"\n"))
		load(file)
		#str(df.all)
		
		#get correct order of sequence SampleCodes corresponding to first seq of Patient
		seq.PROT.nam<- as.character( df.all[,SeqPROT] )
		seq.PROT.nam[ which(is.na(seq.PROT.nam)) ]	<- "NA"
		seq.RT.nam	<- as.character( df.all[,SeqRT] )
		seq.RT.nam[ which(is.na(seq.RT.nam)) ]		<- "NA"
					
		pattern 	<- gsub('/',':',paste(signat.in,".clustalo$",sep=''))
		files		<- list.files(path=indir, pattern=pattern, full.names=1)
		#read all sequences and add a missing one with name "NA" 		
		seq.PROT	<- read.dna( files[ grep("PROT",files) ], format="fa", as.matrix=1 )
		tmp			<- as.DNAbin( matrix(rep('-',ncol(seq.PROT)),1,ncol(seq.PROT), dimnames=list(c("NA"),c())) )
		seq.PROT	<- rbind(seq.PROT,tmp)
		seq.RT		<- read.dna( files[ grep("RT",files) ], format="fa", as.matrix=1 )						 				
		tmp			<- as.DNAbin( matrix(rep('-',ncol(seq.RT)),1,ncol(seq.RT), dimnames=list(c("NA"),c())) )
		seq.RT		<- rbind(seq.RT,tmp)
		
		seq.PROT	<- seq.PROT[seq.PROT.nam,]
		seq.RT		<- seq.RT[seq.RT.nam,]
		rownames(seq.PROT)<- as.character( df.all[,Patient] )
		rownames(seq.RT)<- as.character( df.all[,Patient] )
		
		seq.PROT.RT	<- cbind(seq.PROT,seq.RT)
		if(verbose) print(seq.PROT.RT)
		file		<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_PROTRT_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nwrite to",file))
		save(seq.PROT.RT, file=file)	
	}
	else
	{
		if(verbose)	cat(paste("\nload",file))
		load(file)
	}		
	if(0){	print("DEBUG FirstAliSequences"); seq.PROT.RT<- seq.PROT.RT[1:10,]	}
	if(0)	#create full fasta file with reference sequence and align once more
	{
		file		<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.out),".fasta",sep='')
		if(verbose) cat(paste("\nwrite to ",file))		
		data( refseq_hiv1_hxb2 )
		hxb2		<- data.table( hxb2 )
		hxb2		<- as.character( hxb2[, HXB2.K03455 ] )
		cat( paste(">HXB2\n",paste( hxb2[ seq.int(1,length(hxb2)-2) ], collapse='',sep='' ),"\n",sep=''), file=file, append=1)
		write.dna(seq.PROT.RT, file=file, format="fasta", append=1, colsep='', colw=length(hxb2)-2, blocksep=0)
				
		file<- paste("ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.out),".fasta",sep='')
		cmd<- hivc.cmd.clustalo(outdir, file, signat='')
		system(cmd)		
	}
	if(0)	#curate alignment with reference 
	{
		file								<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.out),".fasta.clustalo",sep='')
		if(verbose) cat(paste("\nread ",file))
		seq.PROT.RT							<- read.dna(file, format="fasta", as.character=1)
		query.yes							<- hivc.seq.find(seq.PROT.RT, 2550, c("-","-","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,2550:2552]<- matrix( c("c","c","y"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 2550, c("c","c","y","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,2550:2553]<- matrix( c("-","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 2751, c("-","-","-","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,2751:2754]<- matrix( c("a","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3143, c("-"))
		for(i in query.yes)	
			seq.PROT.RT[i,3143:3151]		<-	c(seq.PROT.RT[i,3144:3151],"-") 	#align pos 3143 and move gap into 3rd codon
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3151, c("a","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3151:3152]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
		#always delete drug resistance insertion at pos 69 because this is really cryptic
		seq.PROT.RT[,2752:2754]				<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3237, c("-"))
		for(i in query.yes)	
			seq.PROT.RT[i,3233:3237]<- c("-",seq.PROT.RT[i,3233:3236]) 		#align pos 3237 and move gap into 3rd codon
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3433, c("-","t"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3433:3434]<- matrix( c("t","-"), nrow=length(query.yes), ncol=2, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3433, c("-","-","-","-","t"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3433:3437]<- matrix( c("t","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3434, c("-","-","-","-","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3434:3438]<- matrix( c("a","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3503, c("-"))
		for(i in query.yes)
			seq.PROT.RT[i,3501:3503]		<- c("-",seq.PROT.RT[i,3501:3502])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3504, c("-","-","-","t"))	
		for(i in query.yes)
			seq.PROT.RT[i,3504:3524]		<- c(seq.PROT.RT[i,3507:3524],"-","-","-")
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3506, c("-","-","-","t","g"))
		for(i in query.yes)
			seq.PROT.RT[i,3506:3526]		<- c(seq.PROT.RT[i,3509:3526],"-","-","-")
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3507, c("-","-","-","g","a"))
		for(i in query.yes)
			seq.PROT.RT[i,3507:3537]		<-	c(seq.PROT.RT[i,3510:3537],"-","-","-")
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3504, c("g","a","c"))
		for(i in query.yes)
			seq.PROT.RT[i,3504:3530]		<-	c("-","-","-",seq.PROT.RT[i,3504:3527])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3495, c("g","-","-","-","g"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3495:3499]<- matrix( c("g","g","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3501, c("-","t","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3501:3503]<- matrix( c("t","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3501, c("a","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3501:3502]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","g"))
		for(i in query.yes)
			seq.PROT.RT[i,3551:3600]		<-	c("-",seq.PROT.RT[i,3551:3599])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3552, c("g","c"))
		for(i in query.yes)
			seq.PROT.RT[i,3552:3600]		<-	c("-",seq.PROT.RT[i,3552:3599])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3501, c("-","t","t"))
		for(i in query.yes)
			seq.PROT.RT[i,3501:3600]		<-	c("-",seq.PROT.RT[i,3501:3599])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3498, c("c","g","t"))
		for(i in query.yes)
			seq.PROT.RT[i,3498:3520]		<-	c("-",seq.PROT.RT[i,3498:3519])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("-","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3543, c("-","c","a","g","g","g","g","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","g","g","g","g","c","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3543, c("-","c","a","g","g","a","g","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","g","g","a","g","c","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3543, c("-","c","a","k","g","g","a","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","k","g","g","a","c","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3543, c("-","c","a","g","g","g","g","g","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","g","g","g","g","g","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("a","a","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3549:3551]<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","a","c"))
		for(i in query.yes)
			seq.PROT.RT[i,3551:3560]		<-	c("-",seq.PROT.RT[i,3551:3559])
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3551, c("g","a","y","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3551:3554]<- matrix( c("-","g","a","y"), nrow=length(query.yes), ncol=4, byrow=1 )
		query.yes							<- hivc.seq.find(seq.PROT.RT, 3549, c("-","c","c","g","g"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","-","c"), nrow=length(query.yes), ncol=3, byrow=1 )
		#some more manual editing at the end
		#cut before PROT pol start and at max len in database
		seq.PROT.RT							<- seq.PROT.RT[,2253:ncol(seq.PROT.RT)] 
		seq.PROT.RT							<- seq.PROT.RT[,1:1624]
		seq.PROT.RT.sort.by					<- apply(seq.PROT.RT,1,function(x)		which(rev(x)!="-" )[1]  )
		seq.PROT.RT							<- seq.PROT.RT[sort(seq.PROT.RT.sort.by, index.return=1)$ix,]
		#some more manual editing at all ends after sorting
		#always delete drug resistance at G333D/E because this is at the end of the alignment and alignment unreliable - could have picked up other ends
		seq.PROT.RT[,1294:1299]				<- matrix( c("-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=6, byrow=1 )

		seq.PROT.RT							<- as.DNAbin( seq.PROT.RT )
		file								<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.out),".fasta.clustalo2",sep='')		
		write.dna(seq.PROT.RT, file=file, format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
	}
	if(0)	#create final HXB2PROTRT R and phylip files; need phylip for ExaML
	{
		signat.in	<- "Wed_May__1_17/08/15_2013"
		signat.out	<- "Sat_May_11_14:23:46_2013"
		file<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.in),".fasta.clustalo4",sep='')
		seq.PROT.RT	<- read.dna(file, format="fasta")
		print(seq.PROT.RT)
						
		seq.PROT.RT	<- seq.PROT.RT[ rownames(seq.PROT.RT)!="HXB2", ]
		
		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRT_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nwrite to",file))
		save(seq.PROT.RT, file=file)
		
		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRT_",gsub('/',':',signat.out),".phylip",sep='')
		if(verbose) cat(paste("\nwrite to ",file))
		hivc.seq.write.dna.phylip(seq.PROT.RT, file=file)				
	}
	if(1)	#retain only third codon positions
	{
		signat.in	<- "Sat_May_11_14:23:46_2013"
		signat.out	<- "Sat_May_11_14:23:46_2013"

		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRT_",gsub('/',':',signat.in),".R",sep='')
		if(verbose) cat(paste("\nload ",file))			
		load(file)
		
		codon3.idx	<- seq.int(3,ncol(seq.PROT.RT),3)
		seq.PROT.RT3<- seq.PROT.RT[, codon3.idx]
		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRTCD3_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nwrite to ",file))
		save(seq.PROT.RT3, file=file)
		
		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRTCD3_",gsub('/',':',signat.out),".phylip",sep='')
		if(verbose) cat(paste("\nwrite to ",file))
		hivc.seq.write.dna.phylip(seq.PROT.RT3, file=file)
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
		stop()
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
		stop()
	}
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
	
	#compute cluster size distributions
	clusters.distr		<- lapply(clusters,function(x){	tmp<- table(x$size.tips); tmp  	})
	clusters.distr.x	<- lapply(clusters.distr, function(x){	as.numeric(names(x))	})
	clusters.distr.plot	<- merge( thresh, data.table(bs= unique(thresh[,"bs"]), bs.pch= 20+seq_along(unique(thresh[,"bs"]))), all.x=1, by="bs")	
	clusters.distr.plot	<- merge( clusters.distr.plot, data.table(brl= unique(thresh[,"brl"]), brl.col= seq_along(unique(thresh[,"brl"]))), all.x=1, by="brl")
	clusters.distr.plot	<- as.data.table(clusters.distr.plot)
	setkey(clusters.distr.plot, bs)	
	cols				<- brewer.pal( length(unique(clusters.distr.plot[,brl.col])), "RdYlBu")
	
	if(with.withinpatientclu)
		outfile				<- paste("ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_clusterdistr_",dist.brl,sep='')
	else
		outfile				<- paste("ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_transmclusterdistr_",dist.brl,sep='')
	
	outsignat			<- "Sat_Jun_16_17/23/46_2013"
	file				<- paste(dir.name,"tmp",paste(outfile,"_BS=0,6_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	pdf(file=file, width=5,height=5)
	select				<- 1:7		#BS=0.6
	plot(1,1,bty='n',type='n',xlim=range(unlist(clusters.distr.x)), ylim=range((unlist(clusters.distr))), xlab="size",ylab="frequency")
	dummy<- lapply(seq_along(clusters.distr)[select],function(i)
			{
				points(clusters.distr.x[[i]],(clusters.distr[[i]]),type='b',col=cols[clusters.distr.plot[i,brl.col]],pch=clusters.distr.plot[i,bs.pch])
			})
	legend("topright", fill=cols[clusters.distr.plot[1:7,brl.col]], legend= thresh[1:7,"brl"], bty='n', border=NA)
	legend("bottomright", pch=unique(clusters.distr.plot[,bs.pch]), legend= thresh[seq.int(1,len=4,by=7),"bs"], bty='n', border=NA)
	dev.off()
	
	file				<- paste(dir.name,"tmp",paste(outfile,"_BS=0,7_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	pdf(file=file, width=5,height=5)
	select				<- 8:14		#BS=0.7
	plot(1,1,bty='n',type='n',xlim=range(unlist(clusters.distr.x)), ylim=range((unlist(clusters.distr))), xlab="size",ylab="frequency")
	dummy<- lapply(seq_along(clusters.distr)[select],function(i)
			{
				points(clusters.distr.x[[i]],(clusters.distr[[i]]),type='b',col=cols[clusters.distr.plot[i,brl.col]],pch=clusters.distr.plot[i,bs.pch])
			})
	legend("topright", fill=cols[clusters.distr.plot[1:7,brl.col]], legend= thresh[1:7,"brl"], bty='n', border=NA)
	legend("bottomright", pch=unique(clusters.distr.plot[,bs.pch]), legend= thresh[seq.int(1,len=4,by=7),"bs"], bty='n', border=NA)
	dev.off()
	
	file				<- paste(dir.name,"tmp",paste(outfile,"_BS=0,8_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	pdf(file=file, width=5,height=5)
	select				<- 15:21		#BS=0.8
	plot(1,1,bty='n',type='n',xlim=range(unlist(clusters.distr.x)), ylim=range((unlist(clusters.distr))), xlab="size",ylab="frequency")
	dummy<- lapply(seq_along(clusters.distr)[select],function(i)
			{
				points(clusters.distr.x[[i]],(clusters.distr[[i]]),type='b',col=cols[clusters.distr.plot[i,brl.col]],pch=clusters.distr.plot[i,bs.pch])
			})
	legend("topright", fill=cols[clusters.distr.plot[1:7,brl.col]], legend= thresh[1:7,"brl"], bty='n', border=NA)
	legend("bottomright", pch=unique(clusters.distr.plot[,bs.pch]), legend= thresh[seq.int(1,len=4,by=7),"bs"], bty='n', border=NA)
	dev.off()
	
	file				<- paste(dir.name,"tmp",paste(outfile,"_BS=0,9_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	pdf(file=file, width=5,height=5)
	select				<- 22:28		#BS=0.9
	plot(1,1,bty='n',type='n',xlim=range(unlist(clusters.distr.x)), ylim=range((unlist(clusters.distr))), xlab="size",ylab="frequency")
	dummy<- lapply(seq_along(clusters.distr)[select],function(i)
			{
				points(clusters.distr.x[[i]],(clusters.distr[[i]]),type='b',col=cols[clusters.distr.plot[i,brl.col]],pch=clusters.distr.plot[i,bs.pch])
			})
	legend("topright", fill=cols[clusters.distr.plot[1:7,brl.col]], legend= thresh[1:7,"brl"], bty='n', border=NA)
	legend("bottomright", pch=unique(clusters.distr.plot[,bs.pch]), legend= thresh[seq.int(1,len=4,by=7),"bs"], bty='n', border=NA)
	dev.off()
	
		
	file				<- paste(dir.name,"tmp",paste(outfile,"_BRL=0,12_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	pdf(file=file, width=5,height=5)
	select				<- seq.int(7,len=4,by=7)		#BRL=0.12
	plot(1,1,bty='n',type='n',xlim=range(unlist(clusters.distr.x)), ylim=range((unlist(clusters.distr))), xlab="size",ylab="frequency")
	dummy<- lapply(seq_along(clusters.distr)[select],function(i)
			{
				points(clusters.distr.x[[i]],(clusters.distr[[i]]),type='b',col=cols[clusters.distr.plot[i,brl.col]],pch=clusters.distr.plot[i,bs.pch])
			})
	legend("topright", fill=cols[clusters.distr.plot[1:7,brl.col]], legend= thresh[1:7,"brl"], bty='n', border=NA)
	legend("bottomright", pch=unique(clusters.distr.plot[,bs.pch]), legend= thresh[seq.int(1,len=4,by=7),"bs"], bty='n', border=NA)
	dev.off()
}
######################################################################################
project.hivc.clustering.computeclusterstatistics<- function(with.withinpatientclu=0)	
{	
	#load precomputed R objects for clustering 
	infile							<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_preclustlink"		
	insignat						<- "Sat_Jun_16_17/23/46_2013"
	file							<- paste(dir.name,"tmp",paste(infile,'_',gsub('/',':',insignat),".R",sep=''),sep='/')
	load(file)
	
	#see if we can resume computations
	infile							<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_clusterdistrpre"		
	insignat						<- "Sat_Jun_16_17/23/46_2013"
	file							<- paste(dir.name,"tmp",paste(infile,'_',gsub('/',':',insignat),".R",sep=''),sep='/')		
	options(show.error.messages = FALSE)		
	readAttempt<-try(suppressWarnings(load(file)))	
	options(show.error.messages = TRUE)
	#if not, compute	
	if(inherits(readAttempt, "try-error"))
	{
		thresh			<- expand.grid( brl=c(0.015,seq(0.02,0.06,0.01),0.12), bs=seq(0.6,0.9,0.1))
		clusters		<- lapply(list(dist.brl.med,dist.brl.max,dist.brl.casc),function(dist.brl)
							{
								tmp	<- lapply(seq_len(nrow(thresh)), function(i)
										{
											print(thresh[i,])
											hivc.clu.clusterbythresh(ph,thresh.brl=thresh[i,"brl"],dist.brl=dist.brl,thresh.nodesupport=thresh[i,"bs"],nodesupport=ph.node.bs,retval="all")
										})
								names(tmp)<- sapply(seq_len(nrow(thresh)), function(i) paste(thresh[i,],collapse='',sep='') )
								tmp
							})
		names(clusters)	<- 	c("dist.brl.med","dist.brl.max","dist.brl.casc")		
		save(thresh,clusters,file=file)		
	}	
	
	#remove within patient clusters
	if(!with.withinpatientclu)
	{
		#see if we can resume computations
		infile							<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_clusterdistrprewithin"		
		insignat						<- "Sat_Jun_16_17/23/46_2013"
		file							<- paste(dir.name,"tmp",paste(infile,'_',gsub('/',':',insignat),".R",sep=''),sep='/')		
		options(show.error.messages = FALSE)		
		readAttempt<-try(suppressWarnings(load(file)))	
		options(show.error.messages = TRUE)
		#if not, compute	
		if(inherits(readAttempt, "try-error"))
		{		
			tmp			<- lapply(seq_along(clusters),function(i)
							{
								tmp<- lapply(clusters[[i]],function(x){
											clu.onlytp			<- hivc.clu.truepos(x, ph.linked, Ntip(ph))$clu.onlytp
											tmp					<- which( x[["clu.mem"]]%in%clu.onlytp[,clu] )
											x[["clu.mem"]][tmp]	<- NA
											x[["size.all"]]		<- table(x[["clu.mem"]])
											x[["size.tips"]]	<- table(x[["clu.mem"]][1:Ntip(ph)])
											x
										})
								names(tmp)<- names(clusters[[i]])
								tmp
							})
			names(tmp)	<- names(clusters)
			clusters.w	<- tmp
			save(thresh,clusters.w,file=file)	
		}	
	}
	#now have clusters and clusters.w

	#compute cluster distributions etc for clusters
	project.hivc.clustering.computeclusterstatistics.fordistbrl(thresh, clusters.w, with.withinpatientclu=0, dist.brl="dist.brl.med")
	project.hivc.clustering.computeclusterstatistics.fordistbrl(thresh, clusters.w, with.withinpatientclu=0, dist.brl="dist.brl.max")
	project.hivc.clustering.computeclusterstatistics.fordistbrl(thresh, clusters.w, with.withinpatientclu=0, dist.brl="dist.brl.casc")
	
	#compare BS=0.8 BRL=0.06 across dist.brl methods
	outfile				<- paste("ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_transmclusterdistr_dist.brl.cmp",sep='')
	outsignat			<- "Sat_Jun_16_17/23/46_2013"
	file				<- paste(dir.name,"tmp",paste(outfile,"_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')

	select.thresh			<- c(20,22)
	select.brl				<- c("dist.brl.max","dist.brl.casc")
	cluster.d				<- lapply(select.thresh, function(i)
				{
					tmp			<- lapply( select.brl, function(x)	cumsum(table( clusters.w[[x]][[i]][["size.tips"]])	)	)
					names(tmp)	<- select.brl
					tmp
				})
	names(cluster.d)		<- 	apply( thresh[select.thresh,], 1, function(x) paste(x,collapse='_',sep='')	)
	xlim					<- c(2,45)
	ylim					<- c(0, 800)
	pdf(file=file, width=5,height=5)
	plot(1,1,type='n',bty='n',xlim=xlim, ylim=ylim, xlab="cluster size",ylab="number clusters")
	cols					<- brewer.pal(3, "Set1")
	dummy					<- sapply(seq_along(cluster.d), function(i)
			{
				sapply(seq_along(cluster.d[[i]]), function(j)
						{
							lines( as.numeric(names(cluster.d[[i]][[j]])), cluster.d[[i]][[j]], lty=i, col=cols[j] )
						})
			})
	legend("bottomright", legend= select.brl, fill= cols[1:2], bty='n', border=NA)
	legend("bottomleft", legend= names(cluster.d), lty= 1:2, bty='n')
	dev.off()
}		
######################################################################################
project.hivc.clustering.computeTPstatistics<- function(ph, ph.node.bs, dist.brl, ph.linked, ph.unlinked.info, ph.unlinked, outdir=NA, outfile=NA, outsignat=NA, thresh.brl= c(seq(0.02,0.06,0.01),seq(0.08,0.26,0.04)), thresh.bs= c(0.8,0.85,0.9,0.95), resume=1)
{
	file		<- paste(outdir, paste(outfile,"_pre_",gsub('/',':',outsignat),".R",sep=''),sep='/')
	if(resume)												#//load if there is R Master data.table
	{
		options(show.error.messages = FALSE)		
		readAttempt<-try(suppressWarnings(load(file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresume file",file))			
		options(show.error.messages = TRUE)		
	}
	if(!resume || inherits(readAttempt, "try-error"))		#else generate R Master data.table
	{
		cat(paste("\ncreate file",file))
		clusters		<- lapply(thresh.bs, function(bs)
							{
								clusters	<- lapply(thresh.brl, function(brl)
										{
											hivc.clu.clusterbythresh(ph,thresh.brl=brl,dist.brl=dist.brl,thresh.nodesupport=bs,nodesupport=ph.node.bs,retval="all")
										})
								clusters.tp	<- lapply(clusters,function(x)
										{
											hivc.clu.truepos(x, ph.linked, Ntip(ph), verbose=0)
										})		
								clusters.fp	<- lapply(clusters,function(x)
										{
											hivc.clu.trueneg(x, ph.unlinked.info, ph.unlinked, Ntip(ph), verbose=0)
										})
								list(clu= clusters, tp= clusters.tp, fp= clusters.fp)							
							})	
		names(clusters)	<- thresh.bs
		file			<- paste(outdir, paste(outfile,"_pre_",gsub('/',':',outsignat),".R",sep=''),sep='/')
		save(clusters,file=file)
	}
	
	patient.nNL						<- 15700
	clusters.NLseqinclu				<- sapply(seq_along(clusters),function(i)
										{				
											sapply(seq_along(clusters[[i]][["clu"]]), function(j)
													{
														#print(clusters[[i]][["tp"]])
														clustering						<- clusters[[i]][["clu"]][[j]]
														tmp								<- clusters[[i]][["tp"]][[j]][["clu.onlytp"]]
														tmp								<- tmp[,clu]
														tmp								<- which( clustering[["clu.mem"]] %in% tmp )
														clustering[["clu.mem"]][tmp]	<- NA
														clustering[["size.tips"]]		<- table(clustering[["clu.mem"]][1:Ntip(ph)])
														sum(clustering[["size.tips"]]) / patient.nNL
													})
										})
	rownames(clusters.NLseqinclu)	<- thresh.brl
	colnames(clusters.NLseqinclu)	<- thresh.bs
	clusters.NLseqinclu				<- t(clusters.NLseqinclu)
	
	method.tp				<- "tp.rate.tpclu" #"tp.rate.tpavg"
	tp.stat					<- sapply(seq_along(clusters),function(i)
								{				
									sapply(clusters[[i]][["tp"]], function(x) x[method.tp] )
								})
	rownames(tp.stat)		<- thresh.brl
	colnames(tp.stat)		<- thresh.bs
	tp.stat					<- t(tp.stat)
	
	method.fp				<- "fp.n" #"fp.rate"
	fp.stat					<- sapply(seq_along(clusters),function(i)
								{				
									sapply(clusters[[i]][["fp"]], function(x) x[method.fp] )
								})
	rownames(fp.stat)		<- thresh.brl
	colnames(fp.stat)		<- thresh.bs
	fp.stat					<- t(fp.stat)		
	#print(fp.stat); print(tp.stat); print(clusters.NLseqinclu)
		
	cols		<- diverge_hcl(nrow(fp.stat), h = c(246, 40), c = 96, l = c(65, 90))
	names(cols)	<- thresh.bs
	
	file		<- paste(outdir, paste(outfile,"_TPvsFP_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	xlim		<- c(0,16)	#range(fp.stat)
	xlim[2]		<- xlim[2] + 2
	cat(paste("\nplot file",file))
	pdf(file=file,width=5,height=5)
	plot(1,1,type='n',bty='n',xlim=xlim,ylim=range(tp.stat),xlab="#FP",ylab="%TP")
	dummy	<- sapply(seq_len(nrow(fp.stat)),function(i)
			{
				points(fp.stat[i,],tp.stat[i,],col=cols[i],type='b')
				text(fp.stat[i,],tp.stat[i,],col=cols[i],labels=thresh.brl,adj=c(-0.8,0.5),cex=0.5)
			})
	legend("bottomright",border=NA,bty='n',fill=cols,legend=names(cols))
	dev.off()
	
	file		<- paste(outdir, paste(outfile,"_COvsFP_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	cat(paste("\nplot file",file))
	pdf(file=file,width=5,height=5)
	xlim		<- c(0,16)	#range(fp.stat)
	xlim[2]		<- xlim[2] + 2	
	plot(1,1,type='n',bty='n',xlim=xlim,ylim=range(clusters.NLseqinclu),xlab="#FP",ylab="%NL seq in cluster")
	dummy	<- sapply(seq_len(nrow(fp.stat)),function(i)
			{
				points(fp.stat[i,],clusters.NLseqinclu[i,],col=cols[i],type='b')
				text(fp.stat[i,],clusters.NLseqinclu[i,],col=cols[i],labels=thresh.brl,adj=c(-0.8,0.5),cex=0.5)
			})
	legend("bottomright",border=NA,bty='n',fill=cols,legend=names(cols))
	dev.off()
	list(clusters.NLseqinclu=clusters.NLseqinclu, fp.stat=fp.stat, tp.stat=tp.stat)
}
######################################################################################
project.hivc.clustering<- function(dir.name= DATA)
{
	require(ape)
	require(data.table)
	require(RColorBrewer)
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
	if(1)	#test clustering on simple test tree
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
	if(0)	#clustering on ATHENA tree
	{
		opt.dist.brl	<- "dist.brl.max"
		opt.dist.brl	<- "dist.brl.casc"
		opt.dist.brl	<- "dist.brl.med"
		if(1)
		{
			#read tree, get boostrap values and patristic distances between leaves
			infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_preclustlink"
			signat.in	<- "Sat_Jun_16_17/23/46_2013"			
			file		<- paste(dir.name,"tmp",paste(infile,'_',gsub('/',':',signat.in),".R",sep=''),sep='/')
			cat(paste("load file",file))		
			load(file)
		}
		else
		{
			stop(paste("cannot find file",infile))
		}
		dist.brl			<- switch(	opt.dist.brl, 
										"dist.brl.max"		= dist.brl.max,
										"dist.brl.med"		= dist.brl.med,
										"dist.brl.casc"		= dist.brl.casc,
										NA)		
		thresh.bs			<- 0.85
		thresh.brl			<- 0.06				
		clustering			<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
		
		#print(clustering)		
	
		#for each tip, collect 'Patient' and 'TN or PROT+P51 or ATHENA'
		indir				<- paste(dir.name,"tmp",sep='/')
		infile				<- "ATHENA_2013_03_Unlinked_and_Linked"
		insignat			<- "Sat_Jun_16_17/23/46_2013"
		file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
		if(verbose) cat(paste("\nread linked and unlinked from file",file))
		load(file)								
		df.seqinfo			<- merge( data.table(FASTASampleCode= ph$tip.label, Node=seq_along(ph$tip.label), key="FASTASampleCode"), subset(df.seqinfo,select=c(FASTASampleCode, Patient)), all.x=1, by="FASTASampleCode")
		tmp					<- rep("NL", Ntip(ph))
		tmp[ df.seqinfo[,substr(FASTASampleCode,1,2)=="TN"] ]		<- "FRGNTN"
		tmp[ df.seqinfo[,substr(FASTASampleCode,1,8)=="PROT+P51"] ]	<- "FRGN"
		df.seqinfo[,InAthena:=tmp]
				
		#for each tip, collect seq data 
		indir				<- paste(dir.name,"derived",sep='/')
		infile				<- "ATHENA_2013_03_AllSeqPatientCovariates"		
		file				<- paste(indir,'/',infile,".R",sep='')
		if(verbose) cat(paste("\nread patient data from file",file))
		load(file)
		df.seqinfo			<- merge(df.seqinfo, df.all, all.x=1, by="FASTASampleCode")
		setnames(df.seqinfo,"Patient.x","Patient")
		df.seqinfo			<- subset( df.seqinfo,select=-7 )
				
		#focus now only on some seq data entries:
		df.seqinfo			<- subset(df.seqinfo, select=c(Node, FASTASampleCode, InAthena, CountryInfection, Patient, Trm, Sex, PosSeqT , NegT, AnyPos_T1, RegionHospital))	
		#merge InAthena CountryInfection
		tmp					<- which(df.seqinfo[,InAthena=="NL"])
		set(df.seqinfo, tmp, "InAthena", df.seqinfo[tmp,CountryInfection])
		df.seqinfo			<- subset( df.seqinfo,select=-4 )
		setnames(df.seqinfo,"InAthena","CountryInfection")
		
		#add clustering to data table
		setkey(df.seqinfo, Node)
		tmp					<- clustering[["clu.mem"]][seq_len(Ntip(ph))]
		df.seqinfo[,"cluster":=tmp]
		
		outfile				<- paste("ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_clustPHY_",opt.dist.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,sep='')
		outsignat			<- "Sat_Jun_16_17/23/46_2013"
		file				<- paste(dir.name,"tmp",paste(outfile,"_clusterinfo_",gsub('/',':',outsignat),".R",sep=''),sep='/')
		cat(paste("write cluster info to file",file))
		save(df.seqinfo, file=file)
		
		
		#set colors CountryInfection
		tmp					<- rep("transparent",nrow(df.seqinfo))
		df.seqinfo[,CountryInfection.col:=tmp]
		set(df.seqinfo, which(df.seqinfo[,CountryInfection=="NL"]), "CountryInfection.col", "#EF9708")
		#set colors Patient
		tmp					<- unique( df.seqinfo[,Patient] )
		tmp2				<- diverge_hcl(length(tmp), h = c(246, 40), c = 96, l = c(85, 90))		
		tmp					<- data.table(Patient=tmp, Patient.col=tmp2, key="Patient")
		df.seqinfo			<- merge(  df.seqinfo, tmp, all.x=1, by="Patient" )
		#set colors Sex
		tmp					<- rep("transparent",nrow(df.seqinfo))
		df.seqinfo[,Sex.col:=tmp]						
		set(df.seqinfo, which(df.seqinfo[,Sex=="M"]), "Sex.col", "#EF9708")
		#set colors MSM or BI
		tmp					<- rep("transparent",nrow(df.seqinfo))
		df.seqinfo[,Trm.col:=tmp]						
		set(df.seqinfo, which(df.seqinfo[,Trm%in%c("MSM","BI")]), "Trm.col", "#EF9708")
		#set colors RegionHospital
		tmp					<- unique( df.seqinfo[,RegionHospital] )		
		tmp2				<- c("transparent",brewer.pal(length(tmp)-1, "Dark2"))
		tmp					<- data.table(RegionHospital=tmp, RegionHospital.col=tmp2, key="RegionHospital")
		df.seqinfo			<- merge(  df.seqinfo, tmp, all.x=1, by="RegionHospital" )		
		#set colors time		
		tmp					<- range( range(df.seqinfo[, NegT],na.rm=1), range(df.seqinfo[, AnyPos_T1],na.rm=1) )
		tmp					<- as.POSIXlt( seq.Date(tmp[1],tmp[2]+365,by="years") )$year
		tmp2				<- heat_hcl(length(tmp), h = c(0, -100), l = c(75, 40), c = c(40, 80), power = 1)
		yearcols			<- data.table(Year=tmp, Year.col=tmp2)
		#set colors PosSeqT
		tmp					<- data.table(PosSeqT= unique( df.seqinfo[,PosSeqT] ), key="PosSeqT" )
		tmp2				<- tmp[, as.POSIXlt(PosSeqT)$year]
		tmp[,"Year":=tmp2]
		tmp					<- subset( merge(tmp, yearcols, all.x=1, by="Year"), select=c(PosSeqT,Year.col))
		setnames(tmp,"Year.col","PosSeqT.col")
		df.seqinfo			<- merge(  df.seqinfo, tmp, all.x=1, by="PosSeqT" )
		#set colors NegT
		tmp					<- data.table(NegT= unique( df.seqinfo[,NegT] ), key="NegT" )
		tmp2				<- tmp[, as.POSIXlt(NegT)$year]
		tmp[,"Year":=tmp2]
		tmp					<- subset( merge(tmp, yearcols, all.x=1, by="Year"), select=c(NegT,Year.col))
		setnames(tmp,"Year.col","NegT.col")
		df.seqinfo			<- merge(  df.seqinfo, tmp, all.x=1, by="NegT" )
		#set colors AnyPos_T1
		tmp					<- data.table(AnyPos_T1= unique( df.seqinfo[,AnyPos_T1] ), key="AnyPos_T1" )
		tmp2				<- tmp[, as.POSIXlt(AnyPos_T1)$year]
		tmp[,"Year":=tmp2]
		tmp					<- subset( merge(tmp, yearcols, all.x=1, by="Year"), select=c(AnyPos_T1,Year.col))
		setnames(tmp,"Year.col","AnyPos_T1.col")
		df.seqinfo			<- merge(  df.seqinfo, tmp, all.x=1, by="AnyPos_T1" )									
		
		#convert time to string
		set(df.seqinfo,NULL,"AnyPos_T1",substr(as.character( df.seqinfo[,AnyPos_T1] ),1,7))
		set(df.seqinfo,NULL,"NegT",substr(as.character( df.seqinfo[,NegT] ),1,7))
		set(df.seqinfo,NULL,"PosSeqT",substr(as.character( df.seqinfo[,PosSeqT] ),1,7))
				
		#handle missing entries
		setkey(df.seqinfo, Patient)
		tmp					<- which(is.na(df.seqinfo[,Patient]))
		set(df.seqinfo, tmp, "Patient", '')
		set(df.seqinfo, tmp, "Patient.col", "transparent")
		tmp					<- which(is.na(df.seqinfo[,RegionHospital]))
		set(df.seqinfo, tmp, "RegionHospital", '-')
		set(df.seqinfo, tmp, "RegionHospital.col", "transparent")		
		tmp					<- which(is.na(df.seqinfo[,CountryInfection]))
		set(df.seqinfo, tmp, "CountryInfection", "--")
		set(df.seqinfo, tmp, "CountryInfection.col", "transparent")				
		tmp					<- which(is.na(df.seqinfo[,Trm]))
		set(df.seqinfo, tmp, "Trm", '')
		set(df.seqinfo, tmp, "Trm.col", "transparent")
		tmp					<- which(is.na(df.seqinfo[,Sex]))
		set(df.seqinfo, tmp, "Sex", '')
		set(df.seqinfo, tmp, "Sex.col", "transparent")
		tmp					<- which(is.na(df.seqinfo[,AnyPos_T1]))
		set(df.seqinfo, tmp, "AnyPos_T1", "-------")
		set(df.seqinfo, NULL, "AnyPos_T1", paste("HIV+:",df.seqinfo[,AnyPos_T1],sep=''))
		set(df.seqinfo, tmp, "AnyPos_T1.col", "transparent")
		tmp					<- which(is.na(df.seqinfo[,NegT]))
		set(df.seqinfo, tmp, "NegT", "-------")
		set(df.seqinfo, NULL, "NegT", paste("HIV-:",df.seqinfo[,NegT],sep=''))
		set(df.seqinfo, tmp, "NegT.col", "transparent")
		tmp					<- which(is.na(df.seqinfo[,PosSeqT]))
		set(df.seqinfo, tmp, "PosSeqT", "-------")
		set(df.seqinfo, NULL, "PosSeqT", paste("HIVS:",df.seqinfo[,PosSeqT],sep=''))
		set(df.seqinfo, tmp, "PosSeqT.col", "transparent")
		
		setkey(df.seqinfo, Node)				
		#select text and col matrix 
		text				<- t( as.matrix( subset(df.seqinfo,select=c(CountryInfection, Trm, Sex, NegT, AnyPos_T1, PosSeqT, Patient, RegionHospital)) ) )
		col					<- t( as.matrix( subset(df.seqinfo,select=c(CountryInfection.col, Trm.col, Sex.col, NegT.col, AnyPos_T1.col, PosSeqT.col, Patient.col, RegionHospital.col)) ) )
						
		clu.onlytp			<- hivc.clu.truepos(clustering, ph.linked, Ntip(ph))$clu.onlytp
		#double check 'clu.onlytp' 
		if(0)
		{
			dummy<- sapply(clu.onlytp[,clu], function(i)
				{
					tmp<- which( clustering[["clu.mem"]]==i)
					tmp<- tmp[tmp<Ntip(ph)]
					tmp<- merge( data.table(FASTASampleCode= ph$tip.label[ tmp ]), df.seqinfo, by="FASTASampleCode" )
					if( length(unique(tmp[,Patient]))!=1 )
						print(tmp) 					
				})
		}
		
		#plot tree to file
		file		<- paste(dir.name,"tmp",paste(outfile,'_',gsub('/',':',outsignat),".pdf",sep=''),sep='/')
		cat(paste("write tree to file",file))
		hivc.clu.plot(ph, clustering[["clu.mem"]], file=file, pdf.scaley=25, pdf.off=0, highlight.cluster= list( clu.onlytp[,clu] ), highlight.cluster.col="grey50", highlight.edge.of.tiplabel=c("TN","PROT+P51"), highlight.edge.of.tiplabel.col= c("red","blue"), cex.nodelabel=0.1 )
		hivc.clu.plot.tiplabels( seq_len(Ntip(ph)), text, col, cex=0.12, adj=c(-0.15,0.5), add.xinch=0, add.yinch=0 )
		dev.off()
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
	if(1)
	{		
		verbose		<- 1
		resume		<- 1
		#precompute clustering stuff		#KEY1
		indir		<- paste(DATA,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
		insignat	<- "Sat_Jun_16_17/23/46_2013"		
		indircov	<- paste(DATA,"derived",sep='/')
		infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"							
		argv		<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
		argv		<<- unlist(strsplit(argv,' '))
		clu.pre		<- hivc.prog.precompute.clustering()
		#precompute clustering stuff for particular thresholds etc
		opt.brl		<- "dist.brl.casc" 
		thresh.brl	<- 0.096
		thresh.bs	<- 0.8		
		argv		<<- hivc.cmd.clustering(indir, infile, insignat, opt.brl, thresh.brl, thresh.bs, resume=resume)				 
		argv		<<- unlist(strsplit(argv,' '))
		clu			<- hivc.prog.get.clustering()
				
		#remove singletons
		if(verbose) cat(paste("\nnumber of seq in tree is n=", nrow(df.cluinfo)))
		df.cluinfo	<- subset(clu$df.seqinfo, !is.na(cluster) )
		if(verbose) cat(paste("\nnumber of seq in clusters is n=", nrow(df.cluinfo)))
		if(verbose) cat(paste("\nnumber of clusters is n=", length(unique(df.cluinfo[,cluster]))))
		#remove within patient clusters		
		tmp			<- subset(df.cluinfo[,list(clu.is.bwpat=length(unique(Patient))>1),by="cluster"], clu.is.bwpat, cluster )
		df.cluinfo	<- merge(tmp, df.cluinfo, by="cluster", all.x=1)
		if(verbose) cat(paste("\nnumber of seq in clusters between patients is n=", nrow(df.cluinfo)))
		if(verbose) cat(paste("\nnumber of clusters between patients is n=", length(unique(df.cluinfo[,cluster]))))		
		
		#get in-country clusters. this splits clusters with a foreign sequence
		#char.frgn  	='CountryInfection=="FRGN"'; char.frgntn	='CountryInfection=="FRGNTN"'; ph			<- clu.pre$ph; clustering	<- clu$clustering		
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"incountry",sep='')
		outsignat		<- insignat							
		plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep='')
		tmp				<- hivc.clu.getplot.incountry(clu.pre$ph, clu$clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''))
		clu$clustering	<- tmp$clustering
		df.cluinfo		<- tmp$df.cluinfo			
		length(unique(df.cluinfo[,cluster]))
		
		#get msm exposure group clusters. this splits clusters with HET-F
		set(df.cluinfo, which( df.cluinfo[,Trm%in%c("BLOOD","BREAST","PREG","NEEACC")] ), "Trm", "OTH" )
		set(df.cluinfo, which( df.cluinfo[,Trm=="HETfa"] ), "Trm", "HET" )		
		set(df.cluinfo, NULL, "Trm", factor(df.cluinfo[,Trm]) )		
		clustering	<- clu$clustering; 		ph<- clu.pre$ph; 		plot.file	<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''); levels.msm=c("BI","MSM","IDU","NA"); levels.het=c("BI","HET","IDU","NA"); levels.mixed=c("BI","MSM","HET","IDU","NA")
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr",sep='')
		outsignat		<- insignat							
		tmp				<- hivc.clu.getplot.msmexposuregroup(clu.pre$ph, clu$clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''))
		clu$clustering	<- tmp$clustering
		df.cluinfo		<- tmp$df.cluinfo
		
		
		#get branch lengths between all clusters
		tmp					<- hivc.clu.polyphyletic.clusters(clu.pre$ph, clu$clustering, df.cluinfo )
		cluphy.subtrees		<- tmp$cluphy.subtrees 	
		cluphy.brl.bwpat	<- hivc.clu.brl.bwpat(cluphy.subtrees, df.cluinfo)				
		
		#plot clusters with mixed exposure group
		outdir		<- indir
		outfile		<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"mixedexpgr",sep='')
		outsignat	<- insignat		
		hivc.clu.getplot.mixedexposuregroup( clu.pre$ph, clu$clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep='') )
				
		
		
		
		
		#rm NA from clustering[["clu.idx"]], leave in clustering[["clu.mem"]]
		#rebuild subtrees
		
		print(tmp)
		#if more than one cluster, need to introduce new cluster in cluphy.hetF

		
		
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
	if(0)	#investigate various cluster size distributions by size etc
	{	
		#project.hivc.clustering.computeclusterstatistics(with.withinpatientclu=1)
		project.hivc.clustering.computeclusterstatistics(with.withinpatientclu=0)		
		stop()
	}
	if(1)
	{
		#analyze clusters 
		infile							<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_preclustlink"		
		insignat						<- "Sat_Jun_16_17/23/46_2013"
		file							<- paste(dir.name,"tmp",paste(infile,'_',gsub('/',':',insignat),".R",sep=''),sep='/')
		load(file)
		
		#illustrate saturation effect -- level cannot be higher than 90%
		opt.dist.brl					<- "dist.brl.med"
		outdir							<- paste(dir.name,"tmp",sep='/')
		outfile							<- paste("ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_clustTP",opt.dist.brl,sep='_')
		outsignat						<- "Sat_Jun_16_17/23/46_2013"
		stat.med						<- project.hivc.clustering.computeTPstatistics(ph, ph.node.bs, dist.brl.med, ph.linked, ph.unlinked.info, ph.unlinked, outdir=outdir, outfile=outfile, outsignat=outsignat)
		
		opt.dist.brl					<- "dist.brl.max"
		outfile							<- paste("ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_clustTP",opt.dist.brl,sep='_')
		stat.max						<- project.hivc.clustering.computeTPstatistics(ph, ph.node.bs, dist.brl.max, ph.linked, ph.unlinked.info, ph.unlinked, outdir=outdir, outfile=outfile, outsignat=outsignat)
		
		opt.dist.brl					<- "dist.brl.casc"
		outfile							<- paste("ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_clustTP",opt.dist.brl,sep='_')
		stat.casc						<- project.hivc.clustering.computeTPstatistics(ph, ph.node.bs, dist.brl.casc, ph.linked, ph.unlinked.info, ph.unlinked, outdir=outdir, outfile=outfile, outsignat=outsignat)
		#get BRL so that %TP ~ 90%  					 
		#ans80							<- hivc.clu.clusterbytruepos(ph, dist.brl, ph.node.bs, ph.linked, thresh.bs=0.8, level= 0.9, verbose=1)
		#ans90							<- hivc.clu.clusterbytruepos(ph, dist.brl, ph.node.bs, ph.linked, thresh.bs=0.9, level= 0.9, verbose=1)
		#for BS=0.8 BRL is 0.058	-- corresponds to BRL.max 0.106
		#for BS=0.9 BRL is 0.067	-- corresponds to BRL.max 0.096
	
		outfile			<- paste("ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100_clustTP_dist.brl.cmp",sep='_')
		
		select.bs		<- as.character(c(0.8,0.95))
		stat.cmp		<- lapply(list("clusters.NLseqinclu","fp.stat","tp.stat"), function(x)
				{
					tmp	<- lapply(select.bs,function(bs)
							{
								tmp				<- rbind(stat.max[[x]][bs,],stat.casc[[x]][bs,])
								rownames(tmp)	<- c("max","casc")
								tmp
							})
					names(tmp)	<- select.bs
					
					tmp
				})
		names(stat.cmp)	<- c("clusters.NLseqinclu","fp.stat","tp.stat")
		
		thresh.brl		<- as.numeric(colnames(stat.max[[1]]))
		xlim			<- c(0,17)
		ylim			<- c(0.75,0.95)
		cols			<- diverge_hcl(4, h = c(246, 40), c = 96, l = c(65, 90))[c(1,4)]
		names(cols)		<- c("max","casc")
		file			<- paste(outdir, paste(outfile,"_TPvsFP_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
		pdf(file=file,width=5,height=5)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xlab="#FP",ylab="%TP")
		sapply(list("max","casc"),function(x)
				{
					sapply(seq_along(stat.cmp[[1]]),function(j)
							{								
								tmp.x<- as.numeric(stat.cmp[["fp.stat"]][[j]][x,])
								tmp.y<- as.numeric(stat.cmp[["tp.stat"]][[j]][x,])
								points(tmp.x, tmp.y, type='b',col=cols[x], lty=j, pch=15+j)
								text(tmp.x, tmp.y,col=cols[x],labels=thresh.brl,adj=c(-0.8,0.5),cex=0.5)
							})
				})
		legend("bottomright",bty='n',border=NA,fill= cols, legend= names(cols))
		legend("topleft",bty='n',legend= select.bs, pch= 15+1:2)
		dev.off()
		
		xlim			<- c(0,17)
		ylim			<- c(0.05,0.25)
		cols			<- diverge_hcl(4, h = c(246, 40), c = 96, l = c(65, 90))[c(1,4)]
		names(cols)		<- c("max","casc")
		file			<- paste(outdir, paste(outfile,"_COvsFP_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
		pdf(file=file,width=5,height=5)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xlab="#FP",ylab="%NL in cluster")
		sapply(list("max","casc"),function(x)
				{
					sapply(seq_along(stat.cmp[[1]]),function(j)
							{								
								tmp.x<- as.numeric(stat.cmp[["fp.stat"]][[j]][x,])
								tmp.y<- as.numeric(stat.cmp[["clusters.NLseqinclu"]][[j]][x,])
								points(tmp.x, tmp.y, type='b',col=cols[x], lty=j, pch=15+j)
								text(tmp.x, tmp.y,col=cols[x],labels=thresh.brl,adj=c(-0.8,0.5),cex=0.5)
							})
				})
		legend("bottomright",bty='n',border=NA,fill= cols, legend= names(cols))
		legend("topleft",bty='n',legend= select.bs, pch= 15+1:2)
		dev.off()
		
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
					tmp		<- hivc.seq.rmgaps(tmp, verbose=verbose)
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
					outfile	<- paste("clustalo",signat,"qsub",sep='.')					
					hiv.cmd.hpccaller(outdir, outfile, x)			
				})			
	}
}
######################################################################################
hivc.prog.get.clustering<- function()
{
	library(ape)
	library(data.table)	
			
	opt.dist.brl	<- "dist.brl.casc"
	thresh.bs		<- 0.8
	thresh.brl		<- 0.096	#with mut rate 6e-3/yr expect 2*0.048 brl for 8 yrs apart	
	verbose			<- 1
	resume			<- 0
	indir			<- paste(DATA,"tmp",sep='/')
	infile			<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
	insignat		<- "Sat_Jun_16_17/23/46_2013"

	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									thresh.bs= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.bs<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									thresh.brl= return(as.numeric(substr(arg,13,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.brl<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									opt.brl= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.dist.brl<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}	
	
	if(verbose)
	{
		cat(paste("\nopt.dist.brl",opt.dist.brl))
		cat(paste("\nthresh.bs",thresh.bs))
		cat(paste("\nthresh.brl",thresh.brl))
		cat(paste("\nindir",indir))
		cat(paste("\ninfile",infile))
		cat(paste("\ninsignat",insignat))				
	}
	outdir			<- indir
	outfile			<- paste(infile,"_clust_",opt.dist.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,sep='')
	outsignat		<- insignat	
	
	if(resume)
	{
		file		<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')		
		options(show.error.messages = FALSE)		
		if(verbose)
			cat(paste("\ntry to resume file ",file))
		readAttempt<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)
			cat(paste("\nresumed file ",file))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{		
		#
		#	load preclustlink files
		#
		file		<- paste(indir,paste(infile,"_preclust_",gsub('/',':',insignat),".R",sep=''),sep='/')
		if(verbose) cat(paste("load file",file))
		options(show.error.messages = FALSE)
		readAttempt<-	try(suppressWarnings(load(file)))
		if(inherits(readAttempt, "try-error"))	stop(paste("cannot load required file", file))
		options(show.error.messages = TRUE)
		#
		#	generate clustering
		#
		dist.brl		<- switch(	opt.dist.brl, 
									"dist.brl.max"		= dist.brl.max,
									"dist.brl.med"		= dist.brl.med,
									"dist.brl.casc"		= dist.brl.casc,
									NA)
		if(any(is.na(dist.brl)))	stop("unexpected NA in dist.brl")					
		clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
		#
		#	add cluster membership to df.seqinfo
		#
		setkey(df.seqinfo, Node)
		df.seqinfo[,"cluster":= clustering[["clu.mem"]][seq_len(Ntip(ph))]]
		#	set CountryInfection for non-ATHENA seqs
		tmp				<- which( df.seqinfo[,substr(FASTASampleCode,1,2)=="TN"] )
		set(df.seqinfo, tmp, "CountryInfection", "FRGNTN")
		tmp				<- which( df.seqinfo[,substr(FASTASampleCode,1,8)=="PROT+P51"] )
		set(df.seqinfo, tmp, "CountryInfection", "FRGN")
		#
		#	compute TP and TN clusters
		#
		clusters.tp		<- hivc.clu.truepos(clustering, ph.linked, Ntip(ph), verbose= 0)
		clusters.tn		<- hivc.clu.trueneg(clustering, ph.unlinked.info, ph.unlinked, Ntip(ph), verbose=0)
		#
		#	save
		#
		file			<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
		cat(paste("\nwrite cluster info to file",file))
		save(df.seqinfo, clustering, clusters.tp, clusters.tn, file=file)
	}
	
	ans	<- list(thresh.bs=thresh.bs, thresh.brl=thresh.brl, opt.dist.brl=opt.dist.brl, df.seqinfo=df.seqinfo, clustering=clustering, clusters.tp=clusters.tp, clusters.tn=clusters.tn)
	ans
}
######################################################################################
hivc.prog.precompute.clustering<- function()
{
	library(ape)
	#library(adephylo)
	library(data.table)	
	
	verbose		<- 1
	resume		<- 0
	indir		<- paste(DATA,"tmp",sep='/')	
	infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
	insignat	<- "Sat_Jun_16_17/23/46_2013"
	indircov	<- paste(DATA,"derived",sep='/')	
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indircov= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indircov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infilecov= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilecov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}	
	outdir			<- indir
	outsignat		<- insignat
	outfile			<- paste(infile,"preclust",sep='_')
	outfile.dtips	<- paste(infile,"predtips",sep='_')
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infilecov)
		print(resume)
	}
	if(resume)
	{
		file		<- paste(outdir,paste(outfile,'_',gsub('/',':',outsignat),".R",sep=''),sep='/')
		options(show.error.messages = FALSE)		
		if(verbose)
			cat(paste("\ntry to resume file ",file))
		readAttempt<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)
			cat(paste("\nresumed file ",file))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{	
		file							<- paste(indir,paste(infile,'_',gsub('/',':',insignat),".newick",sep=''),sep='/')
		if(verbose) cat(paste("read phylo from file",file))		
		ph								<- ladderize( read.tree(file) )
		gc()
		#
		#	easy: extract bootstrap support
		#
		ph$node.label[2]				<- 0								#little hack so that clustering works
		ph.node.bs						<- as.numeric( ph$node.label )
		ph.node.bs[is.na(ph.node.bs)]	<- 0
		ph.node.bs						<- ph.node.bs/100
		ph$node.label					<- ph.node.bs
		#
		#	memory consuming: extract branch length statistic of subtree
		#
		if(verbose) cat("\ncompute dist.brl.med")
		dist.brl.med					<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=median)
		gc()
		if(verbose) cat("\ncompute dist.brl.max")
		dist.brl.max					<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=max)		#read patristic distances -- this is the expensive step but still does not take very long
		gc()
		if(verbose) cat("\ncompute dist.brl.casc")
		dist.brl.casc					<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
		gc()		
		#
		#	easy: extract tree specific TP and FN data sets
		#		
		file							<- paste(indircov,"/",infilecov,".R",sep='')
		load(file)
		if(verbose) cat(paste("\nfound covariates for patients, n=",nrow(df.all)))
		df.seqinfo						<- subset(df.all, !is.na(PosSeqT) )
		if(verbose) cat(paste("\nfound covariates for patients with non-NA PosSeqT, n=",nrow(df.seqinfo)))
		if(verbose) cat(paste("\nstart: compute TP and TN data tables for phylogeny"))
		tmp								<- hivc.phy.get.TP.and.TN(ph, df.seqinfo, verbose=verbose)
		if(verbose) cat(paste("\nend: compute TP and TN data tables for phylogeny"))
		unlinked.byspace				<- tmp[["unlinked.byspace"]]
		unlinked.bytime					<- tmp[["unlinked.bytime"]]
		linked.bypatient				<- tmp[["linked.bypatient"]]	
		ph.linked						<- tmp[["ph.linked"]]
		ph.unlinked.info				<- tmp[["ph.unlinked.info"]]
		ph.unlinked						<- tmp[["ph.unlinked"]]
		#
		#	add node number to df.seqinfo
		#
		df.seqinfo						<- merge( df.all, data.table(Node=seq_len(Ntip(ph)), FASTASampleCode=ph$tip.label, key="FASTASampleCode" ), all.y=1)
		#
		#	memory consuming: compute branch length matrix between tips
		#
		#if(verbose) cat("\ncompute dist.root")
		#dist.root						<-  distRoot(ph, method= "patristic")
		#gc()
		#if(verbose) cat("\ncompute dist.tips.mat")
		#dist.tips.mat					<-  distTips(ph, method= "patristic")
		#gc()		
		#
		#	save output
		#
		file							<- paste(outdir,paste(outfile,'_',gsub('/',':',outsignat),".R",sep=''),sep='/')	
		if(verbose)	cat(paste("write to file",file))
		save(ph, dist.brl.max, dist.brl.med, dist.brl.casc, ph.node.bs, ph.linked, ph.unlinked.info, ph.unlinked, df.seqinfo, unlinked.byspace, unlinked.bytime, linked.bypatient,  file=file )
		#save(ph, dist.brl.max, dist.brl.med, dist.brl.casc, ph.node.bs, ph.linked, ph.unlinked.info, ph.unlinked, df.seqinfo, unlinked.byspace, unlinked.bytime, linked.bypatient,  file=file )		
	}
	
	ans	<- list(	ph=ph, #dist.tips.mat=dist.tips.mat, #dist.root=dist.root,  
					dist.brl.max=dist.brl.max, dist.brl.med=dist.brl.med, dist.brl.casc=dist.brl.casc, 
					ph.node.bs=ph.node.bs, ph.linked=ph.linked, ph.unlinked.info=ph.unlinked.info, ph.unlinked=ph.unlinked, 
					df.seqinfo=df.seqinfo, unlinked.byspace=unlinked.byspace, unlinked.bytime=unlinked.bytime, linked.bypatient=linked.bypatient
					)
	ans				
}

hivc.prog.remove.resistancemut<- function()
{
	library(ape)
	library(data.table)
	
	#load drug resistance mutations and select unique mutants by codon
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

	indir			<- paste(DATA,"tmp",sep='/')
	infile			<- "ATHENA_2013_03_CurAll+LANL_Sequences"
	insignat		<- "Sat_Jun_16_17/23/46_2013"
	outdir			<- paste(DATA,"tmp",sep='/')
	outfile			<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"
	outsignat		<- "Thu_Aug_01_17/05/23_2013"
	alignment.start	<- 2253	
	verbose			<- 1
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) insignat<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									outfile= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) outsignat<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,16),
									alignment.start= return(substr(arg,18,nchar(arg))),NA)	}))
		if(length(tmp)>0) alignment.start<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
	}
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(outdir)
		print(outfile)		
		print(outsignat)
		print(alignment.start)
		print(verbose)
	}	
	#load alignment
	file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
	load(file)	
	#modify dr table for particular alignment	
	set(dr, NULL, "Alignment.nuc.pos", dr[,Alignment.nuc.pos]-alignment.start+1)
	
	#remove	likely.nonB.outliers
	cat("\nchange infile: remove	likely.nonB.outliers")
	likely.nonB.outliers	<- c("R03-07193","2006G206","PROT+P51_B.AU.1995.C92.AF538307","2008G084")
	likely.nonB.outliers	<- which(rownames(seq.PROT.RT) %in% likely.nonB.outliers)
	seq.PROT.RT				<- seq.PROT.RT[-likely.nonB.outliers,]
	
	#if alignment equals any of the drug resistance mutants, replace with NNN	
	seq.PROT.RT			<- hivc.seq.rm.drugresistance(as.character(seq.PROT.RT), dr, verbose=verbose, rtn.DNAbin=1 )	
	
	#save No Drug resistance alignment to file
	file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
	if(verbose)	cat(paste("\nwrite R file to",file))
	save(seq.PROT.RT, file=file)	
	file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".phylip",sep='')
	if(verbose)	cat(paste("\nwrite phylip file to",file))
	hivc.seq.write.dna.phylip(seq.PROT.RT, file=file)					
	file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".fasta",sep='')			
	if(verbose)	cat(paste("\nwrite fasta file to",file))
	write.dna(seq.PROT.RT, file=file, format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
	
	seq.PROT.RT
}
######################################################################################
hivc.prog.get.geneticdist<- function()
{
	library(bigmemory)
	library(ape)
	
	indir	<- outdir<- paste(DATA,"tmp",sep='/')
	infile	<- "ATHENA_2013_03_FirstAliSequences_PROTRT"
	resume	<- verbose <- 1	
	signat 	<- "Wed_May__1_17/08/15_2013"
	gd.max	<- 0.045 	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))		
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									signat= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									maxgd= return(as.numeric(substr(arg,8,nchar(arg)))),NA)	}))
		if(length(tmp)>0) gd.max<- tmp[1]
	}
	if(1)
	{
		print(indir)
		print(infile)
		print(outdir)
		print(signat)
		print(resume)
		print(gd.max)
	}	
	
	pattern 	<- paste(infile,"_Gd",gd.max*1000,"_",gsub('/',':',signat),".R",sep='')
	file		<- list.files(path=outdir, pattern=pattern, full.names=1)
	if(!resume || !length(file))
	{
		#extract sequences and store as phylip
		#see if gdm file is already there
		pattern 	<- gsub('/',':',paste(signat,".gdm$",sep=''))		
		file		<- list.files(path=outdir, pattern=pattern, full.names=1)		
		if(!resume || !length(file))	
		{		
			if(verbose)	cat(paste("\ncreate gdm file"))
			file				<- paste(indir,"/",infile,"_",gsub('/',':',signat),".R",sep='')
			if(verbose)	cat(paste("\nload",file))
			load(file)
			str(seq.PROT.RT)		
			#tmp				<- tmp[1:10,]
			gd.bigmat			<- hivc.seq.dist(  seq.PROT.RT )
			file				<- paste(outdir,"/",infile,"_",gsub('/',':',signat),".gdm",sep='')
			if(verbose) cat(paste("\nwrite to",file))
			write.big.matrix(gd.bigmat, file, row.names= 1, col.names=0, sep=',')		
			#gd.bigmat.d 		<- describe( gd.bigmat )
			#file				<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_PROTRT_",gsub('/',':',signat),".gdm",sep='')		
			#dput(gd.bigmat.d, file=file)
		}
		else 
		{
			if(verbose)	cat(paste("\nloading",file))
			gd.bigmat			<- read.big.matrix(file, has.row.names=1, sep=',', type="char")
		}
		stop()
		#now have pairwise genetic distances in gd.bigmat
		gd.bigmat.min	<- sapply(seq.int(1,nrow(gd.bigmat)-1),function(i)
								{
									tmp<- which.min(gd.bigmat[i,])
									c(i,tmp, gd.bigmat[i,tmp])
								})
		rownames(gd.bigmat.min)<- c("seq.idx1","seq.idx2","gd")
		gd.bigmat.seq.i	<- which( gd.bigmat.min["gd",] < (gd.max * 1e3) )
		gd.seqs			<- gd.bigmat.min[c("seq.idx1","seq.idx2"),gd.bigmat.seq.i]
		gd.seqs			<- sort(unique(as.vector(gd.seqs)))
		if(verbose) cat(paste("\ncomputed sequences sth at least one pair with less than ",gd.max*100,"% genetic distance, n=",length(gd.seqs)))		
		file			<- paste(indir,"/",infile,"_",gsub('/',':',signat),".R",sep='')
		if(verbose)	cat(paste("\nload",file))
		load(file)
		seq.PROT.RT.gd	<- seq.PROT.RT[ rownames(gd.bigmat)[gd.seqs], ]
		
		file 			<- paste(outdir,"/",infile,"_Gd",gd.max*1000,"_",gsub('/',':',signat),".R",sep='')		
		save(seq.PROT.RT.gd, file=file)
	}
	else
	{
		if(verbose)	cat(paste("\nloading",file))
		load(file)
	}
	#now have seq.PROT.RT.gd
	file			<- paste(outdir,"/ATHENA_2013_03_Gd",gd.max*1000,"Sequences_PROTRT_",gsub('/',':',signat),".phylip",sep='')
	if(verbose) cat(paste("\nwrite to ",file))
	hivc.seq.write.dna.phylip(seq.PROT.RT.gd, file=file)								
}

hivc.proj.pipeline<- function()
{
	dir.name<- DATA		 
	signat.in	<- "Wed_May__1_17/08/15_2013"
	signat.out	<- "Wed_May__1_17/08/15_2013"		
	cmd			<- ''

	if(0)	#align sequences in fasta file with Clustalo
	{
		indir	<- paste(dir.name,"tmp",sep='/')
		outdir	<- paste(dir.name,"tmp",sep='/')
		infile	<- "ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_Wed_May__1_17:08:15_2013.fasta"
		infile	<- "ATHENA_2013_03_All+LANL_Sequences_Sat_Jun_16_17:23:46_2013.fasta"
		cmd		<- hivc.cmd.clustalo(indir, infile, signat='', outdir=outdir)
		cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=1)
	}
	if(0)	#extract first sequences for each patient as available from ATHENA data set
	{
		indir		<- paste(dir.name,"tmp",sep='/')
		infile		<- "ATHENA_2013_03_SeqMaster.R"		
		outdir		<- paste(dir.name,"tmp",sep='/')
		cmd			<- paste(cmd,hivc.cmd.get.firstseq(indir, infile, signat.in, signat.out, outdir=outdir),sep='')
	}
	if(0)	#compute raw pairwise genetic distances accounting correctly for ambiguous IUPAC nucleotides 
	{				
		gd.max		<- 0.045
		signat.in	<- "Sat_May_11_14/23/46_2013"
		signat.out	<- "Sat_May_11_14/23/46_2013"				
		indir		<- paste(dir.name,"tmp",sep='/')
		infile		<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
		outdir		<- paste(dir.name,"tmp",sep='/')
		cmd			<- paste(cmd,hivc.cmd.get.geneticdist(indir, infile, signat.out, gd.max, outdir=outdir),sep='', resume=0)
		
		outfile		<- paste("pipeline",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),"qsub",sep='.')					
		cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q="pqeph")
		cat(cmd)
		hivc.cmd.hpccaller(outdir, outfile, cmd)
		stop()
	}	
	if(1)	#compute bootstrap support and branch length metrics for clustering
	{								
		indir		<- paste(dir.name,"tmp",sep='/')		
		infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
		insignat	<- "Sat_Jun_16_17/23/46_2013"
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"
		insignat	<- "Thu_Aug_01_17/05/23_2013"

		indircov	<- paste(dir.name,"derived",sep='/')
		infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"
		cmd			<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov)
		cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=71, hpc.q="pqeph", hpc.mem="15850mb",  hpc.nproc=1)
		cat(cmd)		
		outdir		<- paste(dir.name,"tmp",sep='/')
		outfile		<- paste("preclust",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),"qsub",sep='.')
		hivc.cmd.hpccaller(outdir, outfile, cmd)
		stop()
	}
	if(0)	#compute one ExaML tree, no bootstrapping
	{		
		indir	<- paste(dir.name,"tmp",sep='/')
		infile	<- "ATHENA_2013_03_FirstAliSequences_PROTRT"
		outdir	<- paste(dir.name,"tmp",sep='/')
		cmd		<- paste(cmd,hivc.cmd.examl(indir,infile,gsub('/',':',signat.out),gsub('/',':',signat.out),outdir=outdir,resume=1,verbose=1),sep='')
		cmd		<- paste(cmd,hivc.cmd.examl.cleanup(outdir),sep='')
	}
	if(0)	#compute ExaML trees with bootstrap values. Bootstrap is over nucleotides in alignment and over initial starting trees to start ML search.
	{
		bs.from		<- 0
		bs.to		<- 2
		bs.n		<- 100
		signat.in	<- "Sat_Jun_16_17:23:46_2013"
		signat.out	<- "Sat_Jun_16_17:23:46_2013"				
		indir		<- paste(dir.name,"tmp",sep='/')
		infile		<- "ATHENA_2013_03_CurAll+LANL_Sequences"

		signat.in	<- "Thu_Aug_01_17/05/23_2013"
		signat.out	<- "Thu_Aug_01_17/05/23_2013"						
		infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"

		outdir		<- paste(dir.name,"tmp",sep='/')
		cmd			<- hivc.cmd.examl.bsalignment(indir,infile,gsub('/',':',signat.out),gsub('/',':',signat.out),bs.from=bs.from,bs.to=bs.to,bs.n=bs.n,outdir=outdir, resume=1, verbose=1)				
		outdir		<- paste(dir.name,"tmp",sep='/')							
		lapply(cmd, function(x)
				{				
					x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=24, hpc.q=NA, hpc.mem="3850mb", hpc.nproc=8)
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("pipeline",signat,"qsub",sep='.')
					#cat(x)					
					hivc.cmd.hpccaller(outdir, outfile, x)
					Sys.sleep(1)
				})
		stop()
	}
	if(0)	#compute ExaML trees with bootstrap values. Bootstrap is over initial starting trees to start ML search.
	{
		bs.from	<- 0
		bs.to	<- 0
		bs.n	<- 100
		signat.in	<- "Sat_May_11_14/23/46_2013"
		signat.out	<- "Sat_May_11_14/23/46_2013"				
		indir	<- paste(dir.name,"tmp",sep='/')
		infile	<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
		infile	<- "ATHENA_2013_03_FirstCurSequences_PROTRTCD3"
		outdir	<- paste(dir.name,"tmp",sep='/')
		cmd		<- hivc.cmd.examl.bsstarttree(indir,infile,gsub('/',':',signat.out),gsub('/',':',signat.out),bs.from=bs.from,bs.to=bs.to,bs.n=bs.n,outdir=outdir, resume=1, verbose=1)
		#check if we have all bs.n files. if yes, combine and cleanup
			
		outdir	<- paste(dir.name,"tmp",sep='/')							
		lapply(cmd, function(x)
				{				
					x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=36, hpc.q="pqeph")
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("pipeline",signat,"qsub",sep='.')
					cat(x)
					#hivc.cmd.hpccaller(outdir, outfile, x)
					#Sys.sleep(1)
				})
		stop()
	}
	
	signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
	outdir	<- paste(dir.name,"tmp",sep='/')
	outfile	<- paste("pipeline",signat,"qsub",sep='.')					
	lapply(cmd, function(x)
			{				
				#x<- hivc.cmd.hpcwrapper(x, hpc.q="pqeph")
				cat(x)
				hivc.cmd.hpccaller(outdir, outfile, x)
			})							
}


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
	stop()
}

