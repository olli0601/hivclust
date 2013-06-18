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

project.hivc.check<- function()
{
	if(1) project.hivc.check.DateRes.after.HIVPosTest()
	if(0) project.hivc.check.DateRes.after.T0()	
}

project.hivc.get.geneticdist.from.sdc<- function(dir.name= DATA)
{	
	tmp<- hivc.clu.geneticdist.cutoff(dir.name=dir.name, plot=1, verbose=1, level.retain.unlinked=0.05)
	print(tmp)
}

project.hivc.collectpatientdata<- function(dir.name= DATA, verbose=1, resume=0)
{	
	require(data.table)
	MISSING.Acute<- c(NA,9)
	
	file.out		<- paste(dir.name,"derived/ATHENA_2013_03_All1stPatientCovariates.R",sep='/')
	#input files generated with "project.hivc.Excel2dataframe"
	file.seq		<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.R",sep='/')
	file.patient	<- paste(dir.name,"derived/ATHENA_2013_03_Patients.R",sep='/')
	file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	file.immu		<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
	
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
		df.pat	<- df
		cat(paste("\n\nadding patient data"))	
		#preprocess df.pat
		df<- data.table( df.pat[,c("Patient","AcuteInfection","MyDateNeg1","MyDateNeg1_Acc","MyDatePos1","MyDatePos1_Acc","isDead","DateDied","Transmission")], key="Patient" )
		tmp<- df[, AcuteInfection]
		tmp[ which( df[,AcuteInfection%in%MISSING.Acute] ) ]<- NA
		df[, isAcute:=tmp]			
	
		df[, "NegT_Acc"]	<- factor(df[,MyDateNeg1_Acc], labels=c("No","Yes"))
		df[, "PosT_Acc"]	<- factor(df[,MyDatePos1_Acc], labels=c("No","Yes"))
		df[, "isAcute"]		<- factor(df[,isAcute], labels=c("No","Yes","Maybe"))
		df[, "isDead"]		<- factor(df[,isDead], labels=c("No","Yes"))	
		df[, "Trm"]			<- factor(df[, Transmission], levels=c(100, 101,  102,  202, 103,  104,  105,  106,  107, 108,  110), labels= c("MSM","BI","HET","HETfa","IDU","BLOOD","NEEACC", "PREG", "BREAST", "OTH", "SXCH"))
		df<- subset(df, select=c("Patient","Trm","isAcute","MyDateNeg1","NegT_Acc","MyDatePos1","PosT_Acc","isDead","DateDied"))
		setnames(df, "MyDateNeg1","NegT")	
		setnames(df, "MyDatePos1","PosT")	
		setnames(df, "DateDied","Died")
		if(0)
		{
			print(table(df[, NegT_Acc]))
			print(table(df[, PosT_Acc]))
			print(table(df[, Trm]))
			print(range(df[, NegT], na.rm=1))
			print(range(df[, PosT], na.rm=1))
			print(range(df[, Died], na.rm=1))
		}				
		df.all<- df[df.minDateRes]
		if(verbose)		cat(paste("\npatients in combined data table", nrow(df.all)))
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
		df.all<- subset(df.all,select=c(Patient,Trm,SeqPROT,SeqRT,PosPROT,PosRT,PosT,PosT_Acc,PosCD41,PosRNA1,NegT,NegT_Acc,isAcute,isDead,Died,RNA1, CD41))
		
		
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
		df.all<- subset(df.all, select=c(Patient,Trm,SeqPROT,SeqRT,PosSeq,PosT,PosT_Acc,PosCD41,PosRNA1,NegT,NegT_Acc,isAcute,isDead,Died,RNA1, CD41))
		
		
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

project.hivc.Excel2dataframe<- function(dir.name= DATA, min.seq.len=21, verbose=1)
{
	if(0)
	{
		#read SEQUENCE csv data file and preprocess		
		analyze			<- 1
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
	}
	if(0)
	{
		NA.time			<- c("","01/01/1911","11/11/1911")		
		#read REGIMEN csv data file and preprocess
		file			<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.csv",sep='/')
		df				<- read.csv(file, stringsAsFactors=FALSE)							
				
		date.var		<- c("T0","StartTime","StopTime")		
		for(x in date.var)
		{
			cat(paste("\nprocess Time", x))
			nok.idx			<- which( df[,x]==NA.time[1] )
			if(verbose)	cat(paste("\nentries with format ",NA.time[1],", n=", length(nok.idx)))
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
		file		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
		if(verbose) cat(paste("\nsave to", file))
		save(df, file=file)		
	}
	if(0)
	{
		NA.time			<- c("","01/01/1911","11/11/1911")		
		NA.transmission	<- 900
		#read PATIENT csv data file and preprocess
		file			<- paste(dir.name,"derived/ATHENA_2013_03_Patient.csv",sep='/')
		df				<- read.csv(file, stringsAsFactors=FALSE)									
		df$isDead		<- as.numeric( df[,"DateDied"]!="")
		df[which(df[,"Transmission"]==NA.transmission),"Transmission"]<- NA
		
		date.var		<- c("DateBorn","MyDateNeg1","MyDatePos1","DateDied","DateLastContact","DateFirstEverCDCC")		
		for(x in date.var)
		{
			cat(paste("\nprocess Time", x))
			nok.idx			<- which( df[,x]==NA.time[1] )
			if(verbose)	cat(paste("\nentries with format ",NA.time[1],", n=", length(nok.idx)))
			if(verbose && length(nok.idx))	cat(paste("\nentries with format ",NA.time[1],", Patient", paste(df[nok.idx,"Patient"],collapse=', ')))
			
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
		file		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
		if(verbose) cat(paste("\nsave to", file))
		save(df, file=file)		
	}
	if(1)
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
		file		<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
		if(verbose) cat(paste("\nsave to", file))
		save(df, file=file)		
	}
}

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
		if(1)
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
		stop()				
	}
	else
	{
		if(verbose)	cat(paste("\nload",file))
		load(file)
	}	
}

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
		query.yes							<- hivx.seq.find(seq.PROT.RT, 2550, c("-","-","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,2550:2552]<- matrix( c("c","c","y"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 2550, c("c","c","y","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,2550:2553]<- matrix( c("-","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 2751, c("-","-","-","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,2751:2754]<- matrix( c("a","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3143, c("-"))
		for(i in query.yes)	
			seq.PROT.RT[i,3143:3151]		<-	c(seq.PROT.RT[i,3144:3151],"-") 	#align pos 3143 and move gap into 3rd codon
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3151, c("a","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3151:3152]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
		#always delete drug resistance insertion at pos 69 because this is really cryptic
		seq.PROT.RT[,2752:2754]				<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3237, c("-"))
		for(i in query.yes)	
			seq.PROT.RT[i,3233:3237]<- c("-",seq.PROT.RT[i,3233:3236]) 		#align pos 3237 and move gap into 3rd codon
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3433, c("-","t"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3433:3434]<- matrix( c("t","-"), nrow=length(query.yes), ncol=2, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3433, c("-","-","-","-","t"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3433:3437]<- matrix( c("t","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3434, c("-","-","-","-","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3434:3438]<- matrix( c("a","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3503, c("-"))
		for(i in query.yes)
			seq.PROT.RT[i,3501:3503]		<- c("-",seq.PROT.RT[i,3501:3502])
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3504, c("-","-","-","t"))	
		for(i in query.yes)
			seq.PROT.RT[i,3504:3524]		<- c(seq.PROT.RT[i,3507:3524],"-","-","-")
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3506, c("-","-","-","t","g"))
		for(i in query.yes)
			seq.PROT.RT[i,3506:3526]		<- c(seq.PROT.RT[i,3509:3526],"-","-","-")
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3507, c("-","-","-","g","a"))
		for(i in query.yes)
			seq.PROT.RT[i,3507:3537]		<-	c(seq.PROT.RT[i,3510:3537],"-","-","-")
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3504, c("g","a","c"))
		for(i in query.yes)
			seq.PROT.RT[i,3504:3530]		<-	c("-","-","-",seq.PROT.RT[i,3504:3527])
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3495, c("g","-","-","-","g"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3495:3499]<- matrix( c("g","g","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3501, c("-","t","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3501:3503]<- matrix( c("t","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3501, c("a","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3501:3502]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3551, c("g","g"))
		for(i in query.yes)
			seq.PROT.RT[i,3551:3600]		<-	c("-",seq.PROT.RT[i,3551:3599])
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3552, c("g","c"))
		for(i in query.yes)
			seq.PROT.RT[i,3552:3600]		<-	c("-",seq.PROT.RT[i,3552:3599])
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3501, c("-","t","t"))
		for(i in query.yes)
			seq.PROT.RT[i,3501:3600]		<-	c("-",seq.PROT.RT[i,3501:3599])
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3498, c("c","g","t"))
		for(i in query.yes)
			seq.PROT.RT[i,3498:3520]		<-	c("-",seq.PROT.RT[i,3498:3519])
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3549, c("-","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3543, c("-","c","a","g","g","g","g","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","g","g","g","g","c","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3543, c("-","c","a","g","g","a","g","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","g","g","a","g","c","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3543, c("-","c","a","k","g","g","a","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","k","g","g","a","c","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3543, c("-","c","a","g","g","g","g","g","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","g","g","g","g","g","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3549, c("a","a","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3549:3551]<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3551, c("g","a","c"))
		for(i in query.yes)
			seq.PROT.RT[i,3551:3560]		<-	c("-",seq.PROT.RT[i,3551:3559])
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3551, c("g","a","y","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3551:3554]<- matrix( c("-","g","a","y"), nrow=length(query.yes), ncol=4, byrow=1 )
		query.yes							<- hivx.seq.find(seq.PROT.RT, 3549, c("-","c","c","g","g"))
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

project.hivc.clustering<- function(dir.name= DATA)
{
	require(data.table)
	require(ape)
	if(0)	#test clustering on simple test tree
	{		
		ph	<- "(Wulfeniopsis:0.196108,(((alpinus:0.459325,grandiflora:0.259364)1.00:0.313204,uniflora:1.155678)1.00:0.160549,(((angustibracteata:0.054609,(brevituba:0.085086,stolonifera:0.086001)0.76:0.035958)1.00:0.231339,(((axillare:0.017540,liukiuense:0.018503)0.96:0.038019,stenostachyum:0.049803)1.00:0.083104,virginicum:0.073686)1.00:0.103843)1.00:0.086965,(carinthiaca:0.018150,orientalis:0.019697)1.00:0.194784)1.00:0.077110)1.00:0.199516,(((((abyssinica:0.077714,glandulosa:0.063758)1.00:0.152861,((((allionii:0.067154,(morrisonicola:0.033595,officinalis:0.067266)1.00:0.055175)1.00:0.090694,(alpina:0.051894,baumgartenii:0.024152,(bellidioides:0.016996,nutans:0.063292)0.68:0.031661,urticifolia:0.032044)0.96:0.036973,aphylla:0.117223)0.67:0.033757,(((japonensis:0.018053,miqueliana:0.033676)1.00:0.160576,vandellioides:0.099761)0.69:0.036188,montana:0.050690)1.00:0.058380)1.00:0.115874,scutellata:0.232093)0.99:0.055014)1.00:0.209754,((((((acinifolia:0.112279,reuterana:0.108698)0.94:0.055829,pusilla:0.110550)1.00:0.230282,((davisii:0.053261,serpyllifolia:0.087290)0.89:0.036820,(gentianoides:0.035798,schistosa:0.038522)0.95:0.039292)1.00:0.092830)1.00:0.169662,(((anagalloides:0.018007,scardica:0.017167)1.00:0.135357,peregrina:0.120179)1.00:0.098045,beccabunga:0.069515)1.00:0.103473)1.00:0.287909,(((((((((((agrestis:0.017079,filiformis:0.018923)0.94:0.041802,ceratocarpa:0.111521)1.00:0.072991,amoena:0.229452,(((argute_serrata:0.017952,campylopoda:0.075210)0.64:0.034411,capillipes:0.022412)0.59:0.034547,biloba:0.037143)1.00:0.141513,intercedens:0.339760,((opaca:0.019779,persica:0.035744)0.94:0.038558,polita:0.036762)1.00:0.108620,rubrifolia:0.186799)1.00:0.144789,(((bombycina_11:0.033926,bombycina_bol:0.035290,cuneifolia:0.017300,jacquinii:0.054249,oltensis:0.045755,paederotae:0.051579,turrilliana:0.017117)0.85:0.049052,czerniakowskiana:0.089983)0.93:0.051111,farinosa:0.138075)1.00:0.080565)1.00:0.104525,((albiflora:0.017984,ciliata_Anna:0.032685,vandewateri:0.017610)0.97:0.045649,arguta:0.063057,(catarractae:0.022789,decora:0.049785)0.96:0.048220,((cheesemanii:0.040125,cupressoides:0.146538)1.00:0.067761,macrantha:0.038130)1.00:0.088158,(densifolia:0.090044,formosa:0.116180)0.71:0.046353,(elliptica:0.038650,(odora:0.019325,salicornioides:0.021228)0.94:0.042950,salicifolia:0.020829)0.92:0.043978,(nivea:0.070429,(papuana:0.035003,tubata:0.031140)0.98:0.064379)0.93:0.065336,raoulii:0.109101)0.97:0.076607)0.93:0.085835,chamaepithyoides:0.485601)0.57:0.072713,(ciliata_157:0.069943,lanuginosa:0.052833)1.00:0.098638,(densiflora:0.069429,macrostemon:0.118926)0.92:0.124911,(fruticulosa:0.086891,saturejoides:0.041181)0.94:0.086148,kellererii:0.083762,lanosa:0.263033,mampodrensis:0.103384,nummularia:0.191180,pontica:0.128944,thessalica:0.129197)0.65:0.031006,(arvensis:0.342138,(((((chamaedrys:0.043720,micans:0.032021,vindobonensis:0.033309)0.51:0.034053,micrantha:0.019084)0.64:0.037906,krumovii:0.020175)1.00:0.103875,verna:0.254017)0.81:0.057105,magna:0.112657)1.00:0.104070)1.00:0.101845)1.00:0.149208,(((aznavourii:0.664103,glauca:0.405588)0.85:0.209945,praecox:0.447238)1.00:0.185614,(donii:0.260827,triphyllos:0.176032)1.00:0.194928)1.00:0.611079)0.74:0.055152,((crista:0.591702,(((cymbalaria_Avlan:0.017401,panormitana:0.017609)1.00:0.229508,((cymbalaria_Istanbul:0.028379,trichadena_332:0.016891,trichadena_Mugla:0.019131)1.00:0.196417,lycica_333:0.146772)1.00:0.097646,lycica_192:0.154877)1.00:0.234748,(((hederifolia:0.018068,triloba:0.075784)1.00:0.084865,(sibthorpioides:0.122542,sublobata:0.136951)1.00:0.074683)0.89:0.043623,stewartii:0.040679)1.00:0.596859)1.00:0.237324)0.58:0.057120,javanica:0.133802)1.00:0.137214)1.00:0.269201,(missurica:0.016685,rubra:0.019696)1.00:0.351184)0.54:0.058275)0.52:0.062485,((dahurica:0.023542,longifolia:0.016484,spicata:0.018125)0.95:0.042294,(nakaiana:0.016270,schmidtiana:0.058451)0.88:0.037207)1.00:0.261643)0.55:0.056458)1.00:0.229509,kurrooa:0.100611)0.74:0.068198,(bonarota:0.040842,lutea:0.115316)1.00:0.241657)0.99:0.085772);"
		ph <- ladderize( read.tree(text = ph) )
		#read bootstrap support values
		thresh.bs						<- 0.9
		ph.node.bs						<- as.numeric( ph$node.label )		
		ph.node.bs[is.na(ph.node.bs)]	<- 0
		ph.node.bs[c(13,15,27,41,43)]	<- thresh.bs-0.05*seq(0.01,length.out=5)
		ph.node.bs[ph.node.bs==1]		<- seq_along(which(ph.node.bs==1))*0.005 + 0.7
		ph$node.label					<- ph.node.bs
		#read patristic distances
		dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf")
		thresh.brl						<- quantile(dist.brl,seq(0.1,1,by=0.05))["100%"]
		print(quantile(dist.brl,seq(0.1,0.5,by=0.05)))
		print(thresh.brl)
		#produce clustering and plot
		clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
		print(clustering)		
		hivc.clu.plot(ph, clustering[["clu.mem"]])
		
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
		infile		<- "ATHENA_2013_03_FirstCurSequences_PROTRT_examlbs100"
		signat.in	<- "Sat_May_11_14/23/46_2013"
		signat.out	<- "Sat_May_11_14/23/46_2013"
		file		<- paste(dir.name,"tmp",paste(infile,'_',gsub('/',':',signat.in),".newick",sep=''),sep='/')
		cat(paste("read file",file))
		ph			<- ladderize( read.tree(file) )
		#read bootstrap support values
		thresh.bs						<- 0.9
		ph.node.bs						<- as.numeric( ph$node.label )
		ph.node.bs[is.na(ph.node.bs)]	<- 0
		ph.node.bs						<- ph.node.bs/100
		ph$node.label					<- ph.node.bs
		print(quantile(ph.node.bs,seq(0.1,1,by=0.1)))
		#read patristic distances -- this is the expensive step but still does not take very long
		dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf")
		thresh.brl						<- quantile(dist.brl,seq(0.1,0.5,by=0.05))["50%"]
		print(quantile(dist.brl,seq(0.1,1,by=0.1)))
		print(thresh.brl)
		#produce clustering and plot
		clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
		print(clustering)
		
		file		<- paste(dir.name,"tmp",paste(infile,'_',gsub('/',':',signat.out),".pdf",sep=''),sep='/')
		cat(paste("write tree to file",file))
		hivc.clu.plot(ph, clustering[["clu.mem"]], file=file, pdf.scaley=25)
		
		#get boostrap support of index cases
		clu.idx			<- clustering[["clu.idx"]]-Ntip(ph)
		print(any(ph.node.bs[clu.idx]<0.9))
		stop()	
	}
	if(1)	#extract unlinked pairs by temporal separation
	{
		verbose				<- 1
		unlinked.closest.n	<- NA
		indir				<- paste(dir.name,"derived",sep='/')
		outdir				<- paste(dir.name,"derived",sep='/')
		infile				<- "ATHENA_2013_03_SeqMaster.R"
		outfile				<- "ATHENA_2013_03_Unlinked_SeroConv_Dead"
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
			file						<- paste(outdir,paste(outfile,"_UnlinkedAll.R",sep=''),sep='/')
		else
			file						<- paste(outdir,paste(outfile,"_Unlinked",unlinked.closest.n,".R",sep=''),sep='/')
		if(verbose) cat(paste("\nwrite unlinked pairs to file",file))
		save(unlinked.bytime, df.serocon, df.all, file=file)		
		stop()
	}
	if(1)	#count how many unlinked pairs in clustering
	{
		require(data.table)
		require(ape)
		
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

	if(1)	#align sequences in fasta file
	{
		indir	<- paste(dir.name,"tmp",sep='/')
		outdir	<- paste(dir.name,"tmp",sep='/')
		infile	<- "ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_Wed_May__1_17:08:15_2013.fasta"
		infile	<- "ATHENA_2013_03_All+LANL_Sequences_Sat_Jun_16_17:23:46_2013.fasta"
		cmd		<- hivc.cmd.clustalo(indir, infile, signat='', outdir=outdir)
		cmd		<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=1)
	}
	if(0)	#extract first sequences for each patient as available
	{
		indir		<- paste(dir.name,"tmp",sep='/')
		infile		<- "ATHENA_2013_03_SeqMaster.R"		
		outdir		<- paste(dir.name,"tmp",sep='/')
		cmd			<- paste(cmd,hivc.cmd.get.firstseq(indir, infile, signat.in, signat.out, outdir=outdir),sep='')
	}
	if(0)	#compute genetic distances
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
	if(0)	#compute ExaML tree
	{		
		indir	<- paste(dir.name,"tmp",sep='/')
		infile	<- "ATHENA_2013_03_FirstAliSequences_PROTRT"
		outdir	<- paste(dir.name,"tmp",sep='/')
		cmd		<- paste(cmd,hivc.cmd.examl(indir,infile,gsub('/',':',signat.out),gsub('/',':',signat.out),outdir=outdir,resume=1,verbose=1),sep='')
		cmd		<- paste(cmd,hivc.cmd.examl.cleanup(outdir),sep='')
	}
	if(0)	#compute ExaML trees with bootstrap values, bootstrap over alignment
	{
		bs.from	<- 0
		bs.to	<- 0
		bs.n	<- 100
		signat.in	<- "Fri_May_24_12/59/06_2013"
		signat.out	<- "Fri_May_24_12/59/06_2013"				
		indir	<- paste(dir.name,"tmp",sep='/')
		infile	<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
		#infile	<- "ATHENA_2013_03_FirstCurSequences_PROTRTCD3"
		outdir	<- paste(dir.name,"tmp",sep='/')
		cmd		<- hivc.cmd.examl.bsalignment(indir,infile,gsub('/',':',signat.out),gsub('/',':',signat.out),bs.from=bs.from,bs.to=bs.to,bs.n=bs.n,outdir=outdir, resume=1, verbose=1)				
		outdir	<- paste(dir.name,"tmp",sep='/')							
		lapply(cmd, function(x)
				{				
					x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=36, hpc.q="pqeph")
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("pipeline",signat,"qsub",sep='.')
					cat(x)
					#stop()
					#hivc.cmd.hpccaller(outdir, outfile, x)
					#Sys.sleep(1)
				})
		stop()
	}
	if(0)	#compute ExaML trees with bootstrap values, bootstrap over starttree
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

