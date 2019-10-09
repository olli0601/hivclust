######################################################################################
project.hivc.Excel2dataframe.AllPatientCovariates.checkPosSeqT<- function()
{
	dir.name	<- DATA
	file		<- paste(dir.name,"derived/ATHENA_2013_03_AllSeqPatientCovariates.R",sep='/')
	load(file)
	
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
	load(file)
	df.treatment<- df
	file		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	load(file)
	df.viro		<- df
	
	tmp			<- subset(df.all, PosSeqT>AnyT_T1, select=c(Patient, FASTASampleCode, PosSeqT))
	tmp			<- merge(tmp, subset(df.treatment, TrVL.failure=='No', select=c(Patient, StartTime, StopTime, lRNA.mx)), by='Patient')
	tmp			<- subset(tmp, PosSeqT>StartTime & PosSeqT<=StopTime)
	tmp2		<- subset(tmp, lRNA.mx<2)
	tmp2		<- merge(tmp2, df.viro,by='Patient')
}
Patients.161027<- function()
{
	#	CONSTANTS
	NA.Acute			<- c(NA,9)
	NA.CountryInfection	<- c(NA,"")
	NA.CountryBorn		<- c(NA,"",'XX')
	NA.RegionOrigin		<- c(NA,"",'XX')
	NA.Subtype			<- c(NA,"")
	NA.time				<- c("","1911-01-01","1911-11-11")		
	date.format			<- "%Y-%m-%d"
	#	INPUT	
	file				<- "~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/original_latest/SHM_1602_161102_OR_ALL_Patient.xlsx"
	file.ggd			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_Region_GGD.rda'
	outfile				<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_Patient.rda'
	#
	#	read data files	
	#
	load(file.ggd)
	dg			<- copy(df)
	require(gdata)
	df			<- read.xls(file, stringsAsFactors=FALSE)
	#
	#	reset Dates
	#
	df			<- hivc.db.reset.Date(df, 	col=c("BIRTH_D","DateNeg1","DatePos1","FirstMed","DateLastContact","DEATH_D","DateAIDS","DateInCare","T0"), 
			NA.time=c("","1911-01-01","1911-11-11"), 
			date.format="%Y-%m-%d")
	#
	#	rename columns
	#
	df			<- as.data.table(df)	
	setnames(df, 	c("DateNeg1","DatePos1","BIRTH_D","DEATH_D","GENDER","ORIGIN","ORIGIN_Region","REG_First_Region","REG_Last_Region","REG_First_GGD","REG_Last_GGD","DateNeg1_Acc","DatePos1_Acc","Acute","MODE_TRANS","REG_Last_Reason","ACS_ID"), 
			c("MyDateNeg1","MyDatePos1","DateBorn","DateDied","Sex","CountryBorn","RegionOrigin","Region_first",'Region_now','GGD_first','GGD_now',"NegT_Acc","PosT_Acc","isAcute","Trm","GGD_LastCh","ACS"))
	#
	#	set factors
	#
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
	set(df, NULL, 'GGD_first', 	df[, factor(GGD_first, 	levels=c(111, 		 			706, 				1009, 		 		1106,			1406, 				1906, 
									2006, 		 			2106, 				2209, 				2406, 			2506, 				2707, 
									3109, 					3406, 				3606, 				3906, 			4106, 				4506, 
									4607, 					4810, 				5006, 				5206, 			5406, 				5608, 
									6011, 					6106, 				7206, 				7306), 
							labels=c('Groningen',			'Drenthe',			'IJsselland',		'Twente',		'Gelre_IJssel',		'Hulpverlening_Gelderland_Midden',
									'Rivierenland',		'Nijmegen', 		'Flevoland',		'Utrecht',		'Midden_Nederland',	'Hollands_Noorden',
									'Kennemerland',		'Amsterdam',		'Gooi_Vechtstreek',	'Den_Haag',		'Zuid_Holland_West','Hollands_Midden',
									'Rotterdam_Rijnmond',	'Zuid_Holland_Zuid','Zeeland',			'West_Brabant',	'Hart_voor_Brabant', 'Brabant_Zuidoost',
									'Limburg-Noord',		'Zuid_Limburg',		'Fryslan',			'Zaanstreek_Waterland'))])
	set(df, NULL, 'GGD_now', 	df[, factor(GGD_now, 	levels=c(111, 		 			706, 				1009, 		 		1106,			1406, 				1906, 
									2006, 		 			2106, 				2209, 				2406, 			2506, 				2707, 
									3109, 					3406, 				3606, 				3906, 			4106, 				4506, 
									4607, 					4810, 				5006, 				5206, 			5406, 				5608, 
									6011, 					6106, 				7206, 				7306), 
							labels=c('Groningen',			'Drenthe',			'IJsselland',		'Twente',		'Gelre_IJssel',		'Hulpverlening_Gelderland_Midden',
									'Rivierenland',		'Nijmegen', 		'Flevoland',		'Utrecht',		'Midden_Nederland',	'Hollands_Noorden',
									'Kennemerland',		'Amsterdam',		'Gooi_Vechtstreek',	'Den_Haag',		'Zuid_Holland_West','Hollands_Midden',
									'Rotterdam_Rijnmond',	'Zuid_Holland_Zuid','Zeeland',			'West_Brabant',	'Hart_voor_Brabant', 'Brabant_Zuidoost',
									'Limburg-Noord',		'Zuid_Limburg',		'Fryslan',			'Zaanstreek_Waterland'))])	
	#
	#	try complete GGD_first etc
	#
	tmp		<- merge(subset(df, is.na(GGD_now), c(Patient)), dg, by='Patient')
	tmp		<- subset(tmp, !is.na(GGD))
	tmp		<- tmp[, list(GGD_now_NEW=GGD[which.max(GGD_Reg)]), by='Patient']
	cat('\nRecovered missing GGD_now for n=', nrow(tmp))
	df		<- merge(df, tmp, all.x=1, by='Patient')
	tmp		<- merge(subset(df, is.na(GGD_first), c(Patient, DateInCare)), dg, by='Patient')
	tmp		<- subset(tmp, !is.na(GGD))	
	tmp		<- tmp[, 	{
							z<- which.min(GGD_Reg)
							list(GGD_first_NEW=GGD[z], GGD_first_Date=GGD_Reg[z])
						}, by='Patient']
	cat('\nRecovered missing GGD_first for n=', nrow(tmp))
	df		<- merge(df, tmp, all.x=1, by='Patient')
	tmp		<- df[, which(!is.na(GGD_now_NEW))]
	set(df, tmp, 'GGD_now', df[tmp,GGD_now_NEW])
	tmp		<- df[, which(!is.na(GGD_first_NEW))]
	set(df, tmp, 'GGD_first', df[tmp,GGD_first_NEW])
	set(df, NULL, c('GGD_now_NEW','GGD_first_NEW'), NULL)
	#
	#	keep biggest cities separate for now
	#
	stopifnot( !nrow(subset(df, Region_first!='' & is.na(GGD_first))) )
	stopifnot( !nrow(subset(df, Region_now!='' & is.na(GGD_now))) )
	set(df, NULL, c('Region_first','Region_now'), NA_character_)
	set(df, df[,which(GGD_first=='Den_Haag')], 'Region_first', 'Den_Haag')
	set(df, df[,which(GGD_now=='Den_Haag')], 'Region_now', 'Den_Haag')	
	set(df, df[,which(GGD_first=='Amsterdam')], 'Region_first', 'Amsterdam')
	set(df, df[,which(GGD_now=='Amsterdam')], 'Region_now', 'Amsterdam')		
	set(df, df[,which(GGD_first=='Rotterdam_Rijnmond')], 'Region_first', 'Rotterdam_Rijnmond')
	set(df, df[,which(GGD_now=='Rotterdam_Rijnmond')], 'Region_now', 'Rotterdam_Rijnmond')
	set(df, df[,which(GGD_first=='Utrecht')], 'Region_first', 'Utrecht')
	set(df, df[,which(GGD_now=='Utrecht')], 'Region_now', 'Utrecht')		
	#	update rest of Region_first etc
	set(df, df[, which(GGD_now%in%c('Groningen','Drenthe','Fryslan'))], 'Region_now', 'North')
	set(df, df[, which(GGD_first%in%c('Groningen','Drenthe','Fryslan'))], 'Region_first', 'North')
	set(df, df[, which(GGD_now%in%c('IJsselland','Twente','Gelre_IJssel','Hulpverlening_Gelderland_Midden','Flevoland',"Rivierenland","Nijmegen","Midden_Nederland"))], 'Region_now', 'East')
	set(df, df[, which(GGD_first%in%c('IJsselland','Twente','Gelre_IJssel','Hulpverlening_Gelderland_Midden','Flevoland',"Rivierenland","Nijmegen","Midden_Nederland"))], 'Region_first', 'East')
	set(df, df[, which(GGD_now%in%c('Gooi_Vechtstreek','Zuid_Holland_Zuid','Zeeland','Hollands_Midden','Kennemerland','Zaanstreek_Waterland','Hollands_Noorden',"Zuid_Holland_West"))], 'Region_now', 'West')
	set(df, df[, which(GGD_first%in%c('Gooi_Vechtstreek','Zuid_Holland_Zuid','Zeeland','Hollands_Midden','Kennemerland','Zaanstreek_Waterland','Hollands_Noorden',"Zuid_Holland_West"))], 'Region_first', 'West')
	set(df, df[, which(GGD_now%in%c('West_Brabant','Hart_voor_Brabant','Brabant_Zuidoost','Limburg-Noord','Zuid_Limburg'))], 'Region_now', 'South')
	set(df, df[, which(GGD_first%in%c('West_Brabant','Hart_voor_Brabant','Brabant_Zuidoost','Limburg-Noord','Zuid_Limburg'))], 'Region_first', 'South')	
	set(df, df[, which(is.na(GGD_first))], 'Region_first', 'Unknown')
	set(df, df[, which(is.na(GGD_now))], 'Region_now', 'Unknown')
	#
	#	
	#
	set(df, df[, which(Trm==900)], 'Trm', NA_integer_)
	set(df, NULL, "Trm", factor(df[, Trm], 	levels= c(	100, 	110, 	150, 	200,  202, 		300,  	400,  		450,  		600,  	620, 		800),
					labels= c(	'MSM',	'BI',	'SXCH',	'HET','HETfa',	'IDU',	'BLOOD',	'NEEACC',	'PREG',	'BREAST',	'OTH')))
	set(df, NULL, 'RegionOrigin', df[, factor(RegionOrigin, levels=c('AUS',		'Car',		'EUC',			'EUO',				'EUW',		'Lat',					'NAM',						'NL',	'OAP',				'SSA',					'ZAz'),
							labels=c('Austr_NZ','Caribbean','Central_EU',	'Eastern_EU_stans',	'Western_EU','Latin_South_America',	'North_Africa_Middle_East',	'NL',	'Oceania_Pacific',	'Sub_Saharan_Africa',	'Sout_SouthEast_Asia'))])
	set(df, NULL, 'GGD_LastCh', df[, factor(GGD_LastCh, levels=c(-1,0,1,2,3,4,5,6,7,8),
							labels=c('unknown','unknown2','other_hospital','objection_patient','objection_parents','death','moved_abroad','lost','objection','protocol'))])
	set(df, df[, which(grepl('unknown',GGD_LastCh))],'GGD_LastCh',NA_integer_)	
	set(df, NULL, 'GGD_LastCh', df[, factor(as.character(GGD_LastCh))])
	set(df, NULL, 'isDead', df[, factor(is.na(DateDied), levels=c(TRUE,FALSE), labels=c("No","Yes"))] )
	set(df, NULL, "NegT_Acc", factor(df[,NegT_Acc], levels=c(0,1), labels=c("No","Yes")))
	set(df, NULL, "PosT_Acc", factor(df[,PosT_Acc], levels=c(0,1), labels=c("No","Yes")))
	#
	#	process Acute fields
	#
	tmp	<- df[, {
				tmp		<- NA_character_
				z		<- Acute_Spec_1%in%c(1L,2L) | Acute_Spec_2%in%c(1L,2L) | Acute_Spec_3%in%c(1L,2L) | Acute_Spec_4%in%c(1L,2L)
				if(!is.na(z) & z)
					tmp	<- 'LAB'
				z		<- Acute_Spec_1%in%c(3L,4L,8L,9L) | Acute_Spec_2%in%c(3L,4L,8L,9L) | Acute_Spec_3%in%c(3L,4L,8L,9L) | Acute_Spec_4%in%c(3L,4L,8L,9L)
				if(is.na(tmp) & !is.na(z) & z)
					tmp	<- 'SYM'
				list(Acute_Spec=tmp)
			}, by='Patient']
	df	<- merge(df, tmp, by='Patient')
	set(df, NULL, 'Acute_Spec', df[, factor(Acute_Spec)])
	setkey(df,Patient)
	str(df)
	#
	#	reset NegT date (checked w Ard)
	#
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
	#df[, table(is.na(isAcute), Acute_Spec)]; subset(df, is.na(isAcute) & !is.na(Acute_Spec))
	#df[, table(isAcute, Acute_Spec, useNA='if')]
	#
	#	reset isAcute
	#
	tmp	<- df[, which(Acute_Spec=='SYM' & is.na(isAcute))]
	cat(paste('\nFound entries with AcuteSpec and is.na(isAcute), n=', length(tmp)))
	set(df, tmp, 'isAcute', 'Maybe') 
	tmp	<- df[, which(Acute_Spec=='LAB' & is.na(isAcute))]
	cat(paste('\nFound entries with AcuteSpec and is.na(isAcute), n=', length(tmp)))
	set(df, tmp, 'isAcute', 'Yes') 	
	tmp	<- df[, which(Acute_Spec=='LAB' & isAcute=='No')]
	cat(paste('\nFound entries with AcuteSpec and isAcute==No, n=', length(tmp)))
	set(df, tmp, 'isAcute', 'Yes') 
	tmp	<- df[, which(Acute_Spec=='SYM' & isAcute=='No')]
	cat(paste('\nFound entries with AcuteSpec and isAcute==No, n=', length(tmp)))
	set(df, tmp, 'isAcute', 'Maybe') 
	#	reset AcuteSpec
	tmp	<- df[, which(is.na(Acute_Spec) & isAcute=='Yes')]
	cat(paste('\nFound entries with is.na(Acute_Spec) and isAcute==Yes, n=', length(tmp)))
	set(df, tmp, 'Acute_Spec', 'CLIN')
	tmp	<- df[, which(Acute_Spec=='SYM' & isAcute=='Yes')]
	cat(paste('\nFound entries with Acute_Spec==SYM and isAcute==Yes, n=', length(tmp)))
	set(df, tmp, 'Acute_Spec', 'CLIN') 	
	tmp	<- df[, which(is.na(Acute_Spec) & isAcute=='Maybe')]
	cat(paste('\nFound entries with is.na(Acute_Spec) and isAcute==Maybe, n=', length(tmp)))
	set(df, tmp, 'Acute_Spec', 'SYM') 
	#
	#
	cat("\nsave to", outfile)	
	save(df, file=outfile)	
}
######################################################################################
Patients.130301<- function(dir.name= DATA, min.seq.len=21, verbose=1)
{
	if(1)
	{
		file				<- paste(dir.name,"derived/ATHENA_2013_03_Patient.csv",sep='/')			#all with sequence
		file.update			<- paste(dir.name,"derived/ATHENA_2014_06_Patient_All.csv",sep='/')
		date.format			<- "%d/%m/%Y"
	}
	if(1)
	{
		file				<- paste(dir.name,"derived/ATHENA_2013_03_Patient_AllMSM.csv",sep='/')	#all MSM
		file.update			<- paste(dir.name,"derived/ATHENA_2014_06_Patient_AllMSM.csv",sep='/')	#all MSM update
		date.format			<- "%d/%m/%Y"
	}
	if(1)
	{
		dir.name			<- "~/Dropbox (SPH Imperial College)/2015_ATHENA_May_Update"
		file				<- NA
		file.update			<- paste(dir.name,"ATHENA_1502_All_Patient.csv",sep='/')
		date.format			<- "%d/%m/%y"
	}
	
	
	NA.Acute			<- c(NA,9)
	NA.CountryInfection	<- c(NA,"")
	NA.CountryBorn		<- c(NA,"",'XX')
	NA.RegionOrigin		<- c(NA,"",'XX')
	NA.Subtype			<- c(NA,"")
	NA.time				<- c("","01/01/1911","11/11/1911")		
	NA.transmission		<- 900
	#	read PATIENT csv data file	
	if(!is.na(file))
	{
		df					<- read.csv(file, stringsAsFactors=FALSE)	
		#
		#	UPDATE	
		#	
		df.update			<- read.csv(file.update, stringsAsFactors=FALSE)
		tmp					<- as.data.table( merge( subset(df, select=c('Patient','Transmission')), subset(df.update, select=c('Patient','Transmission')), by='Patient' ) )
		setkey(tmp, Transmission.x)
		#	Transmission codes changed	!!!
		#	Transmission.x Transmission.y
		#   100            100
		#	101            110
		#	102            200
		#	103            300
		#	104            400
		#	105            450
		#	106            600
		#	108            800
		#	110            150
		#	202            202
		#	900            900
		df.update$Transmission[df.update$Transmission==110]	<- 101
		df.update$Transmission[df.update$Transmission==200]	<- 102
		df.update$Transmission[df.update$Transmission==300]	<- 103
		df.update$Transmission[df.update$Transmission==400]	<- 104
		df.update$Transmission[df.update$Transmission==450]	<- 105
		df.update$Transmission[df.update$Transmission==600]	<- 106
		df.update$Transmission[df.update$Transmission==800]	<- 108
		df.update$Transmission[df.update$Transmission==150]	<- 110	
		cat(paste('\ndifferent variables that are new in update=',paste(setdiff(colnames(df.update), colnames(df)), collapse=' ')))
		#	no update for patients in df that are not in df.update
		tmp					<- merge( data.frame(Patient=setdiff( df[, 'Patient'], df.update[, 'Patient'])), df, by='Patient' )
		tmp$DateAIDS		<- ''
		tmp$Acute			<- tmp$AcuteInfection
		cat(paste('\nFound patients in df that are not in update, n=',nrow(tmp)))
		#	update for all other patients in df
		df					<- merge( subset(df, select=c(Patient, DateFirstEverCDCC)), df.update, by='Patient', all.y=1 )
		for(x in setdiff(colnames(df), colnames(tmp)))
		{
			cat(paste('\nadding new row to those entries with no update',x))
			tmp[[x]]= NA
		}		
		df					<- rbind( df, subset(tmp, select=colnames(df)) )		
	}
	if(is.na(file))
	{
		df.update			<- read.csv(file.update, stringsAsFactors=FALSE)
		df.update$Transmission[df.update$Transmission==100]	<- 101
		df.update$Transmission[df.update$Transmission==110]	<- 101
		df.update$Transmission[df.update$Transmission==150]	<- 110
		df.update$Transmission[df.update$Transmission==200]	<- 102
		df.update$Transmission[df.update$Transmission==202]	<- 202
		df.update$Transmission[df.update$Transmission==300]	<- 103
		df.update$Transmission[df.update$Transmission==400]	<- 104
		df.update$Transmission[df.update$Transmission==450]	<- 105
		df.update$Transmission[df.update$Transmission==600]	<- 106
		df.update$Transmission[df.update$Transmission==800]	<- 108
		df.update$Transmission[df.update$Transmission==900]	<- NA_real_
		df					<- copy(df.update)
	}
	#	
	#
	df$isDead			<- as.numeric( df[,"DateDied"]!="")
	tmp					<- which(df[,"Transmission"]==NA.transmission)
	if(length(tmp))
		df[tmp,"Transmission"]<- NA
	
	#date.var			<- c("DateBorn","MyDateNeg1","MyDatePos1","DateDied","DateLastContact","DateFirstEverCDCC","DateAIDS")
	date.var			<- c("DateBorn","MyDateNeg1","MyDatePos1","DateDied","DateLastContact","DateAIDS","FirstMed","DateInCare","T0")
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
		df[,x]			<- as.Date(df[,x], format=date.format)	
	}
	
	df<- data.table(df)
	#setnames(df, c("MyDateNeg1_Acc","MyDatePos1_Acc","Acute","Transmission","HospitalRegion"), c("NegT_Acc","PosT_Acc","isAcute","Trm","RegionHospital"))
	setnames(df, c("MyDateNeg1_Acc","MyDatePos1_Acc","Acute","Transmission"), c("NegT_Acc","PosT_Acc","isAcute","Trm"))
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
	#set(df, NULL, "RegionHospital", factor(df[,RegionHospital], levels=c(1,2,3,4,5,6), labels=c("Amst","N","E","S","W","Curu")))
	#set(df, NULL, "ReasonStopRegistration", factor(df[,ReasonStopRegistration], levels=c(0,1,2,3,4), labels=c("OptOut","Died","Moved","Lost","PrntsOptOut")))	
	tmp	<- df[, {
				tmp		<- NA_character_
				z		<- Acute_Spec_1%in%c(1L,2L) | Acute_Spec_2%in%c(1L,2L) | Acute_Spec_3%in%c(1L,2L) | Acute_Spec_4%in%c(1L,2L)
				if(!is.na(z) & z)
					tmp	<- 'LAB'
				z		<- Acute_Spec_1%in%c(3L,4L,8L,9L) | Acute_Spec_2%in%c(3L,4L,8L,9L) | Acute_Spec_3%in%c(3L,4L,8L,9L) | Acute_Spec_4%in%c(3L,4L,8L,9L)
				if(is.na(tmp) & !is.na(z) & z)
					tmp	<- 'SYM'
				list(Acute_Spec=tmp)
			}, by='Patient']
	df	<- merge(df, tmp, by='Patient')
	set(df, NULL, 'Acute_Spec', df[, factor(Acute_Spec)])
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
	
	#df[, table(is.na(isAcute), Acute_Spec)]; subset(df, is.na(isAcute) & !is.na(Acute_Spec))
	#df[, table(isAcute, Acute_Spec, useNA='if')]
	
	#reset isAcute
	tmp	<- df[, which(Acute_Spec=='SYM' & is.na(isAcute))]
	cat(paste('\nFound entries with AcuteSpec and is.na(isAcute), n=', length(tmp)))
	set(df, tmp, 'isAcute', 'Maybe') 
	tmp	<- df[, which(Acute_Spec=='LAB' & is.na(isAcute))]
	cat(paste('\nFound entries with AcuteSpec and is.na(isAcute), n=', length(tmp)))
	set(df, tmp, 'isAcute', 'Yes') 	
	tmp	<- df[, which(Acute_Spec=='LAB' & isAcute=='No')]
	cat(paste('\nFound entries with AcuteSpec and isAcute==No, n=', length(tmp)))
	set(df, tmp, 'isAcute', 'Yes') 
	tmp	<- df[, which(Acute_Spec=='SYM' & isAcute=='No')]
	cat(paste('\nFound entries with AcuteSpec and isAcute==No, n=', length(tmp)))
	set(df, tmp, 'isAcute', 'Maybe') 
	#reset AcuteSpec
	tmp	<- df[, which(is.na(Acute_Spec) & isAcute=='Yes')]
	cat(paste('\nFound entries with is.na(Acute_Spec) and isAcute==Yes, n=', length(tmp)))
	set(df, tmp, 'Acute_Spec', 'CLIN')
	tmp	<- df[, which(Acute_Spec=='SYM' & isAcute=='Yes')]
	cat(paste('\nFound entries with Acute_Spec==SYM and isAcute==Yes, n=', length(tmp)))
	set(df, tmp, 'Acute_Spec', 'CLIN') 	
	tmp	<- df[, which(is.na(Acute_Spec) & isAcute=='Maybe')]
	cat(paste('\nFound entries with is.na(Acute_Spec) and isAcute==Maybe, n=', length(tmp)))
	set(df, tmp, 'Acute_Spec', 'SYM') 
	
	if(verbose) cat(paste("\nsave to", paste(substr(file.update, 1, nchar(file.update)-3),'R',sep='')))	
	save(df, file=paste(substr(file.update, 1, nchar(file.update)-3),'R',sep=''))	
}
######################################################################################
Sequences.160227<- function(dir.name= DATA)
{
	dir.name		<- '~/Dropbox (SPH Imperial College)/2015_ATHENA_May_Update'
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.csv",sep='/')
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Sequences_AllMSM.csv",sep='/')
	date.format		<- "%d/%m/%Y"
	file			<- paste(dir.name,"ATHENA_1502_All_Sequences.csv",sep='/')
	date.format		<- "%d/%m/%y"
	#read SEQUENCE csv data file and preprocess				
	verbose			<- 1
	names.GeneCode	<- c("PROT","RT")
	NA.DateRes		<- as.Date("1911-11-11")
	NA.Seq			<- c('')
	min.seq.len		<- 21
	proc.GeneCode	<- c(1,2)
	
	df				<- read.csv(file, stringsAsFactors=FALSE)
	df[,"DateRes"]	<- as.Date(df[,"DateRes"], format=date.format)	
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
				tmp				<- df[ df[,"GeneCode"]==gene, c("Patient","SampleCode","DateRes","Subtype","Sequence"), drop=0 ]								
				tmp
			})		
	names(df)	<- names.GeneCode
	
	file		<- paste(substr(file, 1, nchar(file)-3),'R',sep='')
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)
	quit("no")	
}
######################################################################################
Sequences.addCOMET.161027<- function(dir.name= DATA)
{
	#	run
	#	./comet_osx -bs 1000 < /Users/Oliver/Downloads/ATHENA_1610_All_Sequences.fasta > /Users/Oliver/Downloads/ATHENA_1610_All_Sequences_COMETv2.1.csv
	
	file	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences.rda'
	outfile <- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_withCOMET.rda'
	file.c	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_COMETv2.1.csv'
	
	tmp		<- as.data.table(read.csv(file=file.c, sep='\t', stringsAsFactors=FALSE))
	setnames(tmp, c('name','subtype','bootstrap.support'), c('FASTASampleCode','SUBTYPE_C','SUBTYPE_CBS'))
	set(tmp, tmp[, which(grepl('unassigned',SUBTYPE_C))],'SUBTYPE_C','pot_recombinant')
	tmp2	<- tmp[, which(grepl('check',SUBTYPE_C))]
	set(tmp, tmp2, 'SUBTYPE_C', tmp[tmp2, gsub(' \\(check.*','-check',SUBTYPE_C)])	
	set(tmp, tmp[, which(grepl('_',SUBTYPE_C))],'SUBTYPE_C','pot_recombinant')
	
	load(file)
	ds		<- merge(ds, tmp, by='FASTASampleCode', all.x=1)	
	#ds[, table(SUBTYPE, SUBTYPE_C)]	
	set(ds, ds[, which(SUBTYPE=='02_AG' & SUBTYPE_C=='A1-check')], 'SUBTYPE_C', 'pot_recombinant')
	set(ds, ds[, which(SUBTYPE=='A' & SUBTYPE_C=='A1-check')], 'SUBTYPE_C', 'A1')
	set(ds, ds[, which(SUBTYPE=='A1' & SUBTYPE_C=='A1-check')], 'SUBTYPE_C', 'A1')
	set(ds, ds[, which(SUBTYPE=='REC' & SUBTYPE_C=='A1-check')], 'SUBTYPE_C', 'pot_recombinant')
	set(ds, ds[, which(SUBTYPE=='CRF37CPX' & SUBTYPE_C=='A1-check')], 'SUBTYPE_C', 'pot_recombinant')
	set(ds, ds[, which(SUBTYPE=='A' & SUBTYPE_C=='A1')], 'SUBTYPE', 'A1')
	set(ds, ds[, which(SUBTYPE=='B' & SUBTYPE_C=='B-check')], 'SUBTYPE_C', 'B')
	set(ds, ds[, which(SUBTYPE=='CRF04CPX' & SUBTYPE_C=='B-check')], 'SUBTYPE_C', 'pot_recombinant')	
	set(ds, ds[, which(SUBTYPE_C=='C-check')], 'SUBTYPE_C', 'pot_recombinant')	
	set(ds, ds[, which(SUBTYPE=='REC' & SUBTYPE_C=='D-check')], 'SUBTYPE_C', 'pot_recombinant')
	set(ds, ds[, which(SUBTYPE_C=='F1-check')], 'SUBTYPE_C', 'pot_recombinant')
	set(ds, ds[, which(SUBTYPE=='G' & SUBTYPE_C=='G-check')], 'SUBTYPE_C', 'G')
	set(ds, ds[, which(SUBTYPE=='02_AG' & SUBTYPE_C=='G-check')], 'SUBTYPE_C', 'pot_recombinant')
	set(ds, ds[, which(SUBTYPE=='06_cpx' & SUBTYPE_C=='G-check')], 'SUBTYPE_C', 'pot_recombinant')
	set(ds, ds[, which(SUBTYPE=='CRF11CPX' & SUBTYPE_C=='G-check')], 'SUBTYPE_C', 'pot_recombinant')
	set(ds, ds[, which(SUBTYPE=='REC' & SUBTYPE_C=='G-check')], 'SUBTYPE_C', 'pot_recombinant')	
	set(ds, ds[, which(SUBTYPE_C=='J-check')], 'SUBTYPE_C', 'J')
	set(ds, ds[, which(SUBTYPE_C=='K-check')], 'SUBTYPE_C', 'K')
	set(ds, ds[, which(SUBTYPE_C=='pot_recombinant')], 'SUBTYPE_C', 'REC')
	#	get B and OTHER subtype data sets	
	seq.b		<- seq[subset(ds, SUBTYPE_C=='B'|SUBTYPE=='B')[, FASTASampleCode]]
	seq.other	<- seq[subset(ds, !(SUBTYPE_C=='B'|SUBTYPE=='B'))[, FASTASampleCode]]
	#	save
	save(seq, seq.b, seq.other, ds, file=outfile)
}

######################################################################################
Sequences.addOutgroupToSubtypeB.161121<- function(dir.name= DATA)
{
	require(big.phylo)
	if(0)	#done once with some manual editing 
	{
		#	add outgroup to subtype B sequences
		file.new		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processing_sequences/LANL_outgroup_For_Subtype_B.fasta'
		file.previous	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_sequences_161102/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B.fasta'
		outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processing_sequences/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B_wOutgroup.fasta'
		scp				<- read.dna(file.previous, format='fa')	
		sce				<- read.dna(file.new, format='fa')
		tmp				<- which(grepl('HXB2',rownames(sce)))
		rownames(sce)[tmp]	<- 'HXB2'
		
		scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
		write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
		#
		#	cut ends manually and mask DRMs manually (easy on 8 seqs)
		#	
	}
	if(1)	#next iteration: simply copy from previous file
	{
		file.new		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE_subtype_B.fasta'
		file.previous	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_sequences_161102/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B_wOutgroup.fasta'
		outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_LANL_codonaligned_noDRM_noROGUE_subtype_B_wOutgroup.fasta'
		scp				<- read.dna(file.previous, format='fa')
		scp				<- scp[grepl('OUTGROUP',rownames(scp)),]
		sce				<- read.dna(file.new, format='fa')
		sce				<- sce[,1:1288]
		sce				<- rbind(sce, scp)
		write.dna(sce, file=outfile, format='fa', colsep='', nbcol=-1)
	}
}

######################################################################################
Sequences.align.161027<- function(dir.name= DATA)
{
	file.previous.align	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_data/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Wed_Dec_18_11:37:00_2013.R'
	file				<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_withCOMET_withLANL.rda'
	outfile.ali			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_aligned.fasta'
	outfile.nali		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_notaligned.fasta'
	load(file)
	#
	#	write to file only those sequences that have not yet been aligned
	#
	load(file.previous.align)
	seq.PROT.RT		<- seq.PROT.RT[intersect( rownames(seq.PROT.RT), names(seq) ),]	# exclude previous LANL and seqs w different labels
	tmp				<- seq.lanl[setdiff(names(seq.lanl), rownames(seq.PROT.RT))]
	write.dna(seq.PROT.RT, file=outfile.ali, format='fasta', colsep='', nbcol=-1)
	write.dna(tmp, file=outfile.nali, format='fasta', colsep='', nbcol=-1)	
	tmp2			<- data.table(FROM=seq(1,8501,500), TO=seq(1,8501,500)+499)
	
	#	checking how long alignments can be to get good codon alignment at HIValign
	#write.dna(tmp[1:200], file=gsub('\\.fasta','_200.fasta',outfile.nali), format='fasta', colsep='', nbcol=-1)	
	#write.dna(tmp[1:1000], file=gsub('\\.fasta','_1000.fasta',outfile.nali), format='fasta', colsep='', nbcol=-1)
	
	tmp2[, {
				write.dna(tmp[FROM:min(TO, nrow(tmp))], file=gsub('\\.fasta',paste('_',FROM,'_',min(TO, nrow(tmp)),'.fasta',sep=''),outfile.nali), format='fasta', colsep='', nbcol=-1)				
			}, by='FROM']
	#	two taxa deleted because not recognized as nucleotide at HIValign: 'M2633604022003' '07-916950'
	
	# run these batches manually at HIValign (LANL). Codon-align. Use the curated HMMER pol model from LANL.
	# the curated LANL model produces much better alignments than de-novo models that are learned at runtime with mafft clustalx   	
	# 
	# then edit manually to check for alignment errors (should be very few)
	# then edit manually to remove the pol DRMs insert/del regions
	
	# add HXB2 to prev alignment
	outfile.hxb2	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/hxb2.fasta'
	load('~/git/hivclust/pkg/data/refseq_hiv1_hxb2.rda')
	tmp				<- as.DNAbin(as.matrix(t(as.character(hxb2[2000:4000, HXB2.K03455]))))
	rownames(tmp)	<- 'HXB2'
	write.dna(tmp, file=outfile.hxb2, format='fa', colsep='', nbcol=-1)
	# edited manually. change some gaps to N
	seq.PROT.RT		<- read.dna('~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_HXB2_e.fasta',format='fa')
	seq.PROT.RT		<- as.character(seq.PROT.RT)	
	tmp				<- seq.find(seq.PROT.RT, 504, "-")	
	seq.PROT.RT[tmp,504]<- matrix( c('n'), nrow=length(tmp), ncol=1, byrow=1 )
	tmp				<- seq.find(seq.PROT.RT, 505, "-")	
	seq.PROT.RT[tmp,505]<- matrix( c('n'), nrow=length(tmp), ncol=1, byrow=1 )
	write.dna(seq.PROT.RT, file='~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_HXB2_e2.fasta', format='fa', colsep='', nbcol=-1)
	
	# merge alignments
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_HXB2_e3.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_1_500.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_500.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')	
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_500.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_501_1000.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_1000.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')	
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_1000.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_1001_1500.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_1500.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')	
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_1500.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_1501_2000.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_2000.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')	
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_2000.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_2001_2500.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_2500.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')	
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_2500.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_2501_3000.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_3000.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')	
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_3000.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_3001_3500.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_3500.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')	
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_3500.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_3501_4000.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_4000.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')	
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_4000.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_4001_4500.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_4500.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')	
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_4500.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_4501_5000.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_5000.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')	
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_5000.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_5001_5500.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_5500.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')
	tmp				<- which(grepl('HXB2',rownames(sce)))
	sce				<- sce[-tmp[2],]
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_5500.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_5501_6000.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_6000.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_6000.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_6001_6500.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_6500.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_6500.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_6501_7000.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_7000.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_7000.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_7001_7500.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_7500.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_7500.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_7501_8000.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_8000.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_8000.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_8001_8500.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_8500.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_8500.fasta'
	file.codonedited<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignededited_8501_9000.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmerged_9000.fasta'
	scp				<- read.dna(file, format='fa')	
	sce				<- read.dna(file.codonedited, format='fa')
	scn				<- seq.align.based.on.common.reference(sce, scp, return.common.sites=FALSE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	write.dna(scn, file=outfile, format='fa', colsep='', nbcol=-1)
	#
	#	mask major DRM sites
	#
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_9000_e1.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_9000_e2.fasta'
	s				<- read.dna(file, format='fa')
	s				<- as.character(s)
	s[which(!grepl('HXB2',rownames(s))),105:112]	<- 'n'
	s[which(!grepl('HXB2',rownames(s))),494:502]	<- 'n'
	s				<- as.DNAbin(s)
	write.dna(s, file=outfile, format='fa', colsep='', nbcol=-1)
	#
	#	curate end bit
	#
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_9000_e3.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_9000_e4.fasta'
	s				<- read.dna(file, format='fa')
	s				<- as.character(s)
	s[ which(s[,1184]=='-' & s[,1185]!='-' & s[,1186]=='-'), 1185]<- '-'
	s[ which(s[,1185]=='-' & s[,1186]!='-' & s[,1187]=='-'), 1186]<- '-'
	s[ which(s[,1186]=='-' & s[,1187]!='-' & s[,1188]=='-'), 1187]<- '-'
	s[ which(s[,1187]=='-' & s[,1188]!='-' & s[,1189]=='-'), 1188]<- '-'
	#	remove end bit that does not agree with consensus
	ds				<- data.table(	TAXA=rownames(s),
									LAST=apply(s,1,function(x) max(which(!x%in%c('-','?','n'))))
									)						
	ref.idx			<- which(grepl('HXB2', rownames(s)))		
	ds				<- ds[, {
				z<- rev( s[TAXA,seq.int(1,LAST)]==s[ref.idx,seq.int(1,LAST)] )
				list(LAST=LAST, LAST_AGREE=LAST+1L-which(z)[1])
			}, by='TAXA']
	ds				<- subset(ds, LAST>LAST_AGREE)
	for(i in seq_len(nrow(ds)))
	{
		s[ds[i,TAXA], seq.int(ds[i,LAST_AGREE+1],ds[i,LAST])]<- '-'
	}	
	s				<- as.DNAbin(s)	
	write.dna(s, file=outfile, format='fa', colsep='', nbcol=-1)
	#
	#	curate dangling end bit
	#
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_9000_e7.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_9000_e8.fasta'
	s				<- read.dna(file, format='fa')
	s				<- as.character(s)
	s[ which(s[,1296]=='-' & s[,1297]!='-' & s[,1298]=='-'), 1297]<- '-'
	s				<- as.DNAbin(s)	
	write.dna(s, file=outfile, format='fa', colsep='', nbcol=-1)
	#
	#	cut for our purposes
	#
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_9000_e10.fasta'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_9000_e11.fasta'
	s				<- read.dna(file, format='fa')
	s				<- s[, seq.int(1,1289)]	
	write.dna(s, file=outfile, format='fa', colsep='', nbcol=-1)	
}
######################################################################################
Sequences.select.by.subtype.161027<- function(dir.name= DATA)
{	
	file.seq.noDRM		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_LANL_codonaligned_noDRM.fasta'
	file.info 			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_withCOMET_withLANL.rda'
	s					<- read.dna(file.seq.noDRM, format='fa')
	load(file.info)
	#	141028040728 is pol but starting at pos 2550. Exclude.
	ds					<- subset(ds,FASTASampleCode!='141028040728')
	
	#	B
	tmp					<- subset(ds, SUBTYPE_C=='B')
	tmp2				<- setdiff(tmp[, FASTASampleCode], rownames(s))
	if(length(tmp2))	cat('\ncannot find sequence',tmp2)	
	seq					<- s[tmp[, FASTASampleCode],]
	tmp2				<- merge(subset(tmp, select=FASTASampleCode), db, by='FASTASampleCode')
	tmp2				<- intersect(rownames(s), tmp2[, paste('LANL.',gsub('_(stripped)','',TAXA_L,fixed=TRUE),sep='')])
	tmp2				<- c(tmp2, 'HXB2')
	seq					<- rbind(seq, s[tmp2,])
	write.dna(seq, file=gsub('\\.fasta','_subtype_B.fasta',gsub('All_','',file.seq.noDRM)), format='fa', colsep='', nbcol=-1)
	
	#	A
	tmp					<- subset(ds, SUBTYPE_C=='A1'|SUBTYPE_C=='A2')
	tmp2				<- setdiff(tmp[, FASTASampleCode], rownames(s))
	if(length(tmp2))	cat('\ncannot find sequence',tmp2)	
	seq					<- s[tmp[, FASTASampleCode],]
	tmp2				<- merge(subset(tmp, select=FASTASampleCode), db, by='FASTASampleCode')
	tmp2				<- intersect(rownames(s), tmp2[, paste('LANL.',gsub('_(stripped)','',TAXA_L,fixed=TRUE),sep='')])
	tmp2				<- c(tmp2, 'HXB2')
	seq					<- rbind(seq, s[tmp2,])
	write.dna(seq, file=gsub('\\.fasta','_subtype_A.fasta',gsub('All_','',file.seq.noDRM)), format='fa', colsep='', nbcol=-1)
	
	#	C
	tmp					<- subset(ds, SUBTYPE_C=='C')
	tmp2				<- setdiff(tmp[, FASTASampleCode], rownames(s))
	if(length(tmp2))	cat('\ncannot find sequence',tmp2)	
	tmp					<- subset(tmp, !FASTASampleCode%in%tmp2)
	seq					<- s[tmp[, FASTASampleCode],]
	tmp2				<- merge(subset(tmp, select=FASTASampleCode), db, by='FASTASampleCode')
	tmp2				<- intersect(rownames(s), tmp2[, paste('LANL.',gsub('_(stripped)','',TAXA_L,fixed=TRUE),sep='')])
	tmp2				<- c(tmp2, 'HXB2')
	seq					<- rbind(seq, s[tmp2,])
	write.dna(seq, file=gsub('\\.fasta','_subtype_C.fasta',gsub('All_','',file.seq.noDRM)), format='fa', colsep='', nbcol=-1)
	
	#	D
	tmp					<- subset(ds, SUBTYPE_C=='D')
	tmp2				<- setdiff(tmp[, FASTASampleCode], rownames(s))
	if(length(tmp2))	cat('\ncannot find sequence',tmp2)	
	tmp					<- subset(tmp, !FASTASampleCode%in%tmp2)
	seq					<- s[tmp[, FASTASampleCode],]
	tmp2				<- merge(subset(tmp, select=FASTASampleCode), db, by='FASTASampleCode')
	tmp2				<- intersect(rownames(s), tmp2[, paste('LANL.',gsub('_(stripped)','',TAXA_L,fixed=TRUE),sep='')])
	tmp2				<- c(tmp2, 'HXB2')
	seq					<- rbind(seq, s[tmp2,])
	write.dna(seq, file=gsub('\\.fasta','_subtype_D.fasta',gsub('All_','',file.seq.noDRM)), format='fa', colsep='', nbcol=-1)
	
	#	F
	tmp					<- subset(ds, SUBTYPE_C=='F1'|SUBTYPE_C=='F2')
	tmp2				<- setdiff(tmp[, FASTASampleCode], rownames(s))
	if(length(tmp2))	cat('\ncannot find sequence',tmp2)	
	tmp					<- subset(tmp, !FASTASampleCode%in%tmp2)
	seq					<- s[tmp[, FASTASampleCode],]
	tmp2				<- merge(subset(tmp, select=FASTASampleCode), db, by='FASTASampleCode')
	tmp2				<- intersect(rownames(s), tmp2[, paste('LANL.',gsub('_(stripped)','',TAXA_L,fixed=TRUE),sep='')])
	tmp2				<- c(tmp2, 'HXB2')
	seq					<- rbind(seq, s[tmp2,])
	write.dna(seq, file=gsub('\\.fasta','_subtype_F.fasta',gsub('All_','',file.seq.noDRM)), format='fa', colsep='', nbcol=-1)
	
	#	G
	tmp					<- subset(ds, SUBTYPE_C=='G')
	tmp2				<- setdiff(tmp[, FASTASampleCode], rownames(s))
	if(length(tmp2))	cat('\ncannot find sequence',tmp2)	
	tmp					<- subset(tmp, !FASTASampleCode%in%tmp2)
	seq					<- s[tmp[, FASTASampleCode],]
	tmp2				<- merge(subset(tmp, select=FASTASampleCode), db, by='FASTASampleCode')
	tmp2				<- intersect(rownames(s), tmp2[, paste('LANL.',gsub('_(stripped)','',TAXA_L,fixed=TRUE),sep='')])
	tmp2				<- c(tmp2, 'HXB2')
	seq					<- rbind(seq, s[tmp2,])
	write.dna(seq, file=gsub('\\.fasta','_subtype_G.fasta',gsub('All_','',file.seq.noDRM)), format='fa', colsep='', nbcol=-1)
	
	#
	#	dont do recombinants for now
	#
	if(0)
	{
		#	01_AE
		tmp					<- subset(ds, SUBTYPE=='01_AE')
		tmp2				<- setdiff(tmp[, FASTASampleCode], rownames(s))
		if(length(tmp2))	cat('\ncannot find sequence',tmp2)	
		tmp					<- subset(tmp, !FASTASampleCode%in%tmp2)
		seq					<- s[tmp[, FASTASampleCode],]
		tmp2				<- merge(subset(tmp, select=FASTASampleCode), db, by='FASTASampleCode')
		tmp2				<- intersect(rownames(s), tmp2[, paste('LANL.',gsub('_(stripped)','',TAXA_L,fixed=TRUE),sep='')])
		tmp2				<- c(tmp2, 'HXB2')
		seq					<- rbind(seq, s[tmp2,])
		write.dna(seq, file=gsub('\\.fasta','_subtype_01AE.fasta',gsub('All_','',file.seq.noDRM)), format='fa', colsep='', nbcol=-1)
		
		#	02_AG
		tmp					<- subset(ds, SUBTYPE=='02_AG')
		tmp2				<- setdiff(tmp[, FASTASampleCode], rownames(s))
		if(length(tmp2))	cat('\ncannot find sequence',tmp2)	
		tmp					<- subset(tmp, !FASTASampleCode%in%tmp2)
		seq					<- s[tmp[, FASTASampleCode],]
		tmp2				<- merge(subset(tmp, select=FASTASampleCode), db, by='FASTASampleCode')
		tmp2				<- intersect(rownames(s), tmp2[, paste('LANL.',gsub('_(stripped)','',TAXA_L,fixed=TRUE),sep='')])
		tmp2				<- c(tmp2, 'HXB2')
		seq					<- rbind(seq, s[tmp2,])
		write.dna(seq, file=gsub('\\.fasta','_subtype_02AG.fasta',gsub('All_','',file.seq.noDRM)), format='fa', colsep='', nbcol=-1)
		
		#	06_cpx
		tmp					<- subset(ds, SUBTYPE=='06_cpx')
		tmp2				<- setdiff(tmp[, FASTASampleCode], rownames(s))
		if(length(tmp2))	cat('\ncannot find sequence',tmp2)	
		tmp					<- subset(tmp, !FASTASampleCode%in%tmp2)
		seq					<- s[tmp[, FASTASampleCode],]
		tmp2				<- merge(subset(tmp, select=FASTASampleCode), db, by='FASTASampleCode')
		tmp2				<- intersect(rownames(s), tmp2[, paste('LANL.',gsub('_(stripped)','',TAXA_L,fixed=TRUE),sep='')])
		tmp2				<- c(tmp2, 'HXB2')
		seq					<- rbind(seq, s[tmp2,])
		write.dna(seq, file=gsub('\\.fasta','_subtype_06CPX.fasta',gsub('All_','',file.seq.noDRM)), format='fa', colsep='', nbcol=-1)
	}
}
######################################################################################
Sequences.removeDRMs.161027<- function(dir.name= DATA)
{
	require(big.phylo)
	file				<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonalignedmergededited_9000_e11.fasta'
	outfile				<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_codonaligned_noDRM.fasta'
	s					<- read.dna(file, format='fa')
	tmp					<- which(grepl("HXB2",rownames(s)))
	rownames(s)[tmp]	<- 'HXB2'
	tmp					<- big.phylo:::seq.rm.drugresistance(s)
	nodr.info			<- tmp$nodr.info
	seq					<- tmp$nodr.seq
	write.dna(seq, file= outfile, format='fasta')
	save(seq, nodr.info, file= gsub('fasta','rda',outfile))
}
######################################################################################
Sequences.addcloseLANL.161027<- function(dir.name= DATA)
{
	#
	#	the blast db does not accept gaps '-' 
	lfile		<- '~/Downloads/LANL_background_pol_n64000_selected.fasta'
	s			<- read.dna(lfile, format='fa')
	s2			<- seq.strip.gap(s, strip.pc=0.99998)
	s2			<- seq.replace(s2, code.from='-', code.to='n', verbose=0)
	write.dna(s2,file='~/Downloads/LANL_background_pol_n64000_strippedn.fasta', format='fasta', colsep='', nbcol=-1)
	# run
	#	makeblastdb -in /Users/Oliver/Downloads/LANL_background_polany_n64000_strippedn.fasta -dbtype nucl -title LANL_background_polany_n64000_strippedn -out /Users/Oliver/Downloads/LANL_background_polany_n64000_strippedn.db
	lfile		<- '~/Downloads/LANL_background_protrt300complete_n25455_stripped.fasta'
	s2			<- read.dna(lfile, format='fa')
	s2			<- seq.replace(s2, code.from='-', code.to='n', verbose=0)
	write.dna(s2,file='~/Downloads/LANL_background_protrt300complete_n25455_strippedn.fasta', format='fasta', colsep='', nbcol=-1)
	# run
	#	makeblastdb -in /Users/Oliver/Downloads/LANL_background_protrt300complete_n25455_strippedn.fasta -dbtype nucl -title LANL_background_protrt300complete_n25000_strippedn -out /Users/Oliver/Downloads/LANL_background_protrt300complete_n25000_strippedn.db
	
	#run
	#	blastn -query /Users/Oliver/Downloads/ATHENA_1610_All_Sequences.fasta -db /Users/Oliver/Downloads/LANL_background_protrt300complete_n25000_strippedn.db -out /Users/Oliver/Downloads/ATHENA_1610_All_Sequences.txt -max_target_seqs 5 -outfmt 6
	
	file		<- "~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_withCOMET.rda"	
	bfile		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_Closest5LANL_protrt300complete.txt'
	lfile		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/LANL_background_protrt300complete_n25455_strippedn.fasta'
	outfile		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_withCOMET_withLANL.rda'
	outfile.B	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_B_LANL.fasta'
	outfile.O	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_Other_LANL.fasta'	
	outfile.nali<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Sequences_LANL_notaligned.fasta'
	load(file)
	db			<- as.data.table(read.csv(bfile, sep='\t',stringsAsFactors=FALSE, header=FALSE))		
	setnames(db, paste('V',1:12,sep=''),c('FASTASampleCode','TAXA_L','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'))
	#	get the 3 closest LANL sequences
	setkey(db, FASTASampleCode, bitscore)
	tmp			<- db[, list(TAXA_L=TAXA_L[(length(TAXA_L)-2):length(TAXA_L)] ), by='FASTASampleCode']
	db			<- merge(db, tmp, by=c('FASTASampleCode','TAXA_L'))
	db			<- merge(db, subset(ds, select=c(FASTASampleCode, SUBTYPE, SUBTYPE_C)), by='FASTASampleCode', all=1)
	
	#	select closest sequences for all
	dbc			<- unique(subset(db, select=TAXA_L))	
	dbs			<- read.dna(lfile, format='fa')
	dbs			<- dbs[dbc[, TAXA_L],]
	rownames(dbs)	<- paste('LANL.',gsub('_(stripped)','',rownames(dbs),fixed=1),sep='')
	seq.lanl	<- c(seq, as.list(dbs))
	
	#	select closest sequences for subtype B
	dbc			<- unique(subset(db, SUBTYPE=='B'|SUBTYPE_C=='B', TAXA_L))	
	dbs			<- read.dna(lfile, format='fa')
	dbs			<- dbs[dbc[, TAXA_L],]
	rownames(dbs)	<- paste('LANL.',gsub('_(stripped)','',rownames(dbs),fixed=1),sep='')
	seq.b.lanl	<- c(seq.b, as.list(dbs))
	
	#	select closest sequences for other subtypes
	dbc			<- unique(subset(db, !(SUBTYPE=='B'|SUBTYPE_C=='B'), TAXA_L))	
	dbs			<- read.dna(lfile, format='fa')
	dbs			<- dbs[dbc[, TAXA_L],]
	rownames(dbs)	<- paste('LANL.',gsub('_(stripped)','',rownames(dbs),fixed=1),sep='')
	seq.other.lanl	<- c(seq.other, as.list(dbs))
	
	#
	#	write to file
	#
	save(ds, db, seq, seq.b, seq.other, seq.lanl, seq.b.lanl, seq.other.lanl, file=outfile)
	write.dna(seq.b.lanl, file=outfile.B, format='fasta', colsep='', nbcol=-1)
	write.dna(seq.other.lanl, file=outfile.O, format='fasta', colsep='', nbcol=-1)
}
######################################################################################
Sequences.161027<- function(dir.name= DATA)
{
	#	read SEQUENCE csv data file and preprocess
	dir.name		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update'
	file			<- paste(dir.name,"original","SHM_1602_161013_OR_ALL_Sequences.xlsx",sep='/')
	outfile			<- paste(dir.name,"preprocessed","ATHENA_1610_All_Sequences.rda",sep='/')
	#					
	verbose			<- 1	
	NA.DateRes		<- "1911-11-11"
	NA.Seq			<- c('')
	min.seq.len		<- 21
	proc.GeneCode	<- c(1,2)
	
	require(gdata)
	require(data.table)
	
	df				<- read.xls(file, stringsAsFactors=FALSE)
	df				<- as.data.table(df)
	setnames(df, colnames(df), toupper(colnames(df)))
	set(df, NULL, 'GENECODE', df[, as.character(factor(GENECODE, levels=c('1','2'), labels=c('PROT','RT')))])
	#	process sequencing dates
	tmp				<- df[, which(SEQ_D%in%NA.DateRes)]
	if(verbose) cat(paste("\nentries with missing sequence dates, n=", length(tmp)))
	set(df, tmp, 'SEQ_D', NA_character_)		
	set(df, NULL, 'SEQ_D', df[,as.Date(SEQ_D, format="%Y-%m-%d")])
	if(verbose) cat(paste("\nrange of sequence dates is",paste(range(df[,SEQ_D], na.rm=1),collapse=', ')))	
	if(verbose) cat(paste("\nentries with missing sequence dates, SampleCode", paste(df[is.na(SEQ_D), SEQ_ID],collapse=', ')))
	#	define partial pol
	tmp				<- df[, {
				if(length(SEQUENCE)>2) stop('more than two sequence fragments for same ID?', SEQ_ID[1])
				z	<- ''
				if(any(GENECODE=='PROT'))
					z<- paste(z,SEQUENCE[GENECODE=='PROT'],sep='')
				if(any(GENECODE=='RT'))
					z<- paste(z,SEQUENCE[GENECODE=='RT'],sep='')				
				list(	GENECODE='PARTIALPOL', 
						SEQ_D=SEQ_D[1], 
						SUBTYPE=ifelse(length(unique(SUBTYPE))==1, SUBTYPE[1], 'discordant'), 
						SEQUENCE=z)
			}, by=c('PATIENT','SEQ_ID')]
	df			<- rbind(df, tmp)
	#	define gene length
	tmp			<- df[, list(SEQ_L=nchar(SEQUENCE)), by=c('PATIENT','SEQ_ID','GENECODE')]	
	df			<- merge(df, tmp, by=c('PATIENT','SEQ_ID','GENECODE'))
	#	warning if no sequence
	if(verbose) cat('\nwarning: no sequence found for', subset(df, GENECODE=='PARTIALPOL' & SEQ_L==0)[, paste(SEQ_ID, collapse=', ')])
	df			<- subset(df, SEQ_L>0)
	#	fixup seq IDs
	set(df, NULL, 'SEQ_ID', df[, gsub(' ','-',SEQ_ID)])
	
	#
	#	done reading. transform data to two objects: sequence info (data.table) and sequences (DNAbin)
	#
	ds			<- subset(df, GENECODE=='PARTIALPOL')
	ds[, GENECODE:=NULL]	
	set(ds, NULL, 'SEQUENCE', ds[, gsub('?','n',tolower(SEQUENCE),fixed=1)])
	
	seq			<- as.DNAbin(strsplit(ds[, SEQUENCE],''))
	names(seq)	<- ds[,SEQ_ID]
	ds[, SEQUENCE:=NULL]
	setnames(ds, c('PATIENT','SEQ_ID','SEQ_D'), c('Patient','FASTASampleCode','PosSeqT'))	
	#
	#	save Dutch sequences
	#	
	save(ds, seq, file=outfile)		
	write.dna(seq, file=gsub('rda','fasta',outfile), format='fasta', colsep='', nbcol=-1)
	quit("no")	
}
######################################################################################
CD4.130301<- function(dir.name= DATA, verbose=1)
{
	file			<- paste(dir.name,"derived/ATHENA_2013_03_Immu.csv",sep='/')
	#file			<- paste(dir.name,"derived/ATHENA_2013_03_Immu_AllMSM.csv",sep='/')
	
	verbose			<- 1
	#read CD4 csv data file and preprocess
	df				<- read.csv(file, stringsAsFactors=FALSE)											
	df				<- project.hivc.Excel2dataframe.CD4.process.NA( df, NA.time=c("","01/01/1911","11/11/1911","24/06/1923"), date.var=c("DateImm") )
	df		<- data.table(df, key="Patient")
	if(!grepl('AllMSM',file))
		set(df, NULL, 'Source', NA_integer_)	
	setnames(df, "DateImm","PosCD4")
	df		<- project.hivc.Excel2dataframe.CD4.process.QC( df )
	df		<- subset(df,!is.na(PosCD4) & !is.na(CD4A))
	if(verbose) cat(paste("\nnumber of entries read, n=",nrow(df)))	
	#
	#	update CD4 measurements	
	#
	tmp				<- ifelse(	!grepl('AllMSM',file),	
			paste(dir.name,"derived/ATHENA_2013_03_Immu_AllMSM.csv",sep='/'),
			paste(dir.name,"derived/ATHENA_2013_03_Immu.csv",sep='/')
	)
	df.up		<- read.csv(tmp, stringsAsFactors=FALSE)
	df.up		<- project.hivc.Excel2dataframe.CD4.process.NA( df.up, NA.time=c("","01/01/1911","11/11/1911","24/06/1923"), date.var=c("DateImm") )
	df.up		<- data.table(df.up, key="Patient")
	setnames(df.up, "DateImm","PosCD4")	
	df.up		<- project.hivc.Excel2dataframe.CD4.process.QC(df.up)
	df.up		<- subset(df.up, !is.na(PosCD4) & !is.na(CD4A) )	
	#	all distinct viro measurements combined		
	df.up		<- merge(df, df.up, by=c('Patient','PosCD4','CD4A'), all.x=1, all.y=1, allow.cartesian=1)		
	tmp			<- nrow(merge( unique(subset(df, select=Patient)), df.up, by='Patient' ))
	if(verbose) cat(paste("\nnumber of entries updated [should be larger or equal than read], n=",tmp))		
	#
	#	if duplicate per day, keep those with higher Source variable
	#	
	if(grepl('AllMSM',file))	#just to make sure we can use same code in both cases
		setnames(df.up, c('Source'), c('Source.y'))	
	df.up		<- df.up[, {
				tmp	<- which(Source.y==max(Source.y))
				list(CD4A=CD4A[tmp], Source=Source.y[tmp])	
			}, by=c('Patient','PosCD4')]
	df			<- merge( unique(subset(df, select=Patient)), df.up, by='Patient' )			
	#	remove identical entries
	setkey(df, Patient, PosCD4, CD4A)
	df		<- unique(df)
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
	setkey(df, Patient, PosCD4)
	df		<- df[, list(CD4A=round(mean(CD4A)), Source=Source[1] ), by=c('Patient','PosCD4')]	
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
Viro.checkVersions<- function()		
{
	file				<- paste(dir.name,"derived/ATHENA_2013_03_Viro_AllMSM.R",sep='/')
	load(file)
	df2					<- df
	file				<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	load(file)
	tmp					<- merge(subset(df, select=c(Patient, PosRNA, RNA)), subset(df2, select=c(Patient, PosRNA, RNA)), by=c('Patient','PosRNA'), all.x=1, all.y=1)
	tmp					<- subset( tmp, RNA.x!=RNA.y )
	tmp[, d:=abs(log(RNA.x/RNA.y))]
	setkey(tmp, d)
	subset( tmp, abs(log(RNA.x/RNA.y))>log(2) )	
}
######################################################################################
CD4.checkVersions<- function()		
{
	file				<- paste(dir.name,"derived/ATHENA_2013_03_Immu_AllMSM.R",sep='/')
	load(file)
	df2					<- df
	file				<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
	load(file)
	tmp					<- merge(subset(df, select=c(Patient, PosCD4, CD4)), subset(df2, select=c(Patient, PosCD4, CD4)), by=c('Patient','PosCD4'), all.x=1, all.y=1)
	tmp					<- subset( tmp, CD4.x!=CD4.y )	
}
######################################################################################
#	QC corrections from Ard
Viro.process.QC.Corrections<- function(df)
{
	tmp		<- which( df[, Patient=="M32210" & as.character(PosRNA)=="2010-05-10"] )
	set(df,tmp,"PosRNA", 50.)
	set(df,tmp,"Undetectable", 'No')
	tmp		<- which( df[, Patient=="M36914" & as.character(PosRNA)=="2010-04-12"] )
	set(df,tmp,"PosRNA", NA_real_)	
	df
}
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
Viro.130301<- function()		
{
	if(0)
	{
		file.treatment		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
		file				<- paste(dir.name,"derived/ATHENA_2013_03_Viro.csv",sep='/')			
	}
	if(0)
	{
		file.treatment		<- paste(dir.name,"derived/ATHENA_2013_03_Regimens_AllMSM.R",sep='/')
		file				<- paste(dir.name,"derived/ATHENA_2013_03_Viro_AllMSM.csv",sep='/')		
	}	
	
	verbose				<- 1
	dir.name			<- DATA
	DB.locktime			<- HIVC.db.locktime
	
	#need for checking of VL data
	
	load(file.treatment)
	df.treat			<- subset(df, select=c(Patient, StartTime, StopTime, AnyT_T1, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel))
	
	
	RNA.min					<- 50	#400	#seems to be standard value
	RNA.min.b4T				<- 400
	RNA.stdvl.udetect.aTS	<- 1e3
	RNA.max					<- 5e6
	lRNA.min.infectious		<- log10(1e4)
	lRNA.min.early			<- log10(1e5)
	lRNA.max.b4early		<- log10(2e4)
	
	#	read VIROLOGY csv data file and preprocess	
	df				<- read.csv(file, stringsAsFactors=FALSE)										
	df				<- project.hivc.Excel2dataframe.Viro.process.NA(df, NA.time=c("","01/01/1911","11/11/1911","24/06/1923"), date.var=c("DateRNA")	)	
	df				<- data.table(df, key="Patient")
	set(df, NULL, 'Undetectable', df[, factor(Undetectable, levels=c(0,1,2),labels=c("No","Yes","LargerThan"))])
	if(!grepl('AllMSM',file))
		set(df, NULL, 'Source', NA_integer_)
	setnames(df, "DateRNA","PosRNA")	
	df				<- project.hivc.Excel2dataframe.Viro.process.QC.Corrections(df)
	df				<- subset(df, !is.na(RNA) & !is.na(PosRNA) & !is.na(RNA))
	if(verbose) cat(paste("\nnumber of entries read, n=",nrow(df)))
	#
	#	update Viro measurements	
	#
	tmp				<- ifelse(	!grepl('AllMSM',file),	
			paste(dir.name,"derived/ATHENA_2013_03_Viro_AllMSM.csv",sep='/'),
			paste(dir.name,"derived/ATHENA_2013_03_Viro.csv",sep='/')
	)
	df.up		<- read.csv(tmp, stringsAsFactors=FALSE)
	df.up		<- project.hivc.Excel2dataframe.Viro.process.NA(df.up, NA.time=c("","01/01/1911","11/11/1911","24/06/1923"), date.var=c("DateRNA")	)
	df.up		<- data.table(df.up, key="Patient")
	set(df.up, NULL, 'Undetectable', df.up[, factor(Undetectable, levels=c(0,1,2),labels=c("No","Yes","LargerThan"))])			
	setnames(df.up, "DateRNA","PosRNA")		
	df.up		<- project.hivc.Excel2dataframe.Viro.process.QC.Corrections(df.up)
	df.up		<- subset(df.up, !is.na(RNA) & !is.na(PosRNA) & !is.na(RNA))
	#	all distinct viro measurements combined		
	df.up		<- merge(df, df.up, by=c('Patient','PosRNA','RNA'), all.x=1, all.y=1, allow.cartesian=1)		
	tmp			<- nrow(merge( unique(subset(df, select=Patient)), df.up, by='Patient' ))
	if(verbose) cat(paste("\nnumber of entries updated [should be larger or equal than read], n=",tmp))		
	#
	#	corrections
	#
	set(df, df[, which(Patient=='M12818' & PosRNA=='1996-09-18')], 'RNA', 400)
	set(df, df[, which(Patient=='M32210' & PosRNA=='2010-10-05')], 'RNA', 50)
	set(df, df[, which(Patient=='M36914' & PosRNA=='2010-12-04')], 'RNA', 0)
	#
	#	if duplicate per day, keep those with higher Source variable
	#	
	if(grepl('AllMSM',file))	#just to make sure we can use same code in both cases
		setnames(df.up, c('Source','Undetectable.x','Undetectable.y'), c('Source.y','Undetectable.y','Undetectable.x'))
	
	df.up		<- df.up[, {
				tmp	<- which(Source.y==max(Source.y))
				list(RNA=RNA[tmp], Source=Source.y[tmp],  Undetectable.x=Undetectable.x[tmp], Undetectable.y=Undetectable.y[tmp])	
			}, by=c('Patient','PosRNA')]
	#	inconsistent Undetectable --> set to Yes				
	tmp			<- df.up[, which( !is.na(Undetectable.x) & !is.na(Undetectable.y) & Undetectable.x!=Undetectable.y) ]
	cat(paste('\ndf.up: Inconsistent Undetectable for data versions. Set to Undetectable==Yes, n=', length(tmp)))
	set(df.up, tmp, c('Undetectable.x','Undetectable.y'), 'Yes')
	
	tmp			<- df.up[, which( !is.na(Undetectable.x) & is.na(Undetectable.y)) ]
	cat(paste('\ndf.up: !is.na(Undetectable.x) & is.na(Undetectable.y), [should be zero] n=', length(tmp)))
	setnames(df.up, 'Undetectable.y', 'Undetectable')
	df.up[, Undetectable.x:=NULL]
	df			<- merge( unique(subset(df, select=Patient)), df.up, by='Patient' )			
	#
	#	Handle duplicates with Undetectable Yes/No: rm undetectables='yes'
	#
	tmp		<- subset( df[, list(n=length(RNA), nD=length(unique(Undetectable)) ), by=c('Patient','PosRNA')], n>1 & nD>1)	
	if(verbose) cat(paste("\nFound duplicate entries, one of which is undetectable. Remove Undetectable one, n=",nrow(tmp)))
	setnames(tmp, c('Patient','PosRNA'), c('xPatient','xPosRNA'))
	for(i in seq_len(nrow(tmp)))
		set(df, df[, which(Patient==tmp[i,xPatient] & PosRNA==tmp[i,xPosRNA] & Undetectable=='Yes')], 'RNA', NA_real_ )
	df		<- subset(df, !is.na(RNA))
	#	remove identical entries
	setkey(df, Patient, PosRNA, RNA)
	df		<- unique(df)
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
	if(verbose)		print( subset(df, Patient=="M12735") )
	tmp		<- which( df[, Patient=="M12735" & as.character(PosRNA)=="1986-04-29"] )
	set(df,tmp,"PosRNA",NA)		
	if(verbose)		print( subset(df, Patient=="M29967") )
	tmp		<- which( df[, Patient=="M29967" & as.character(PosRNA)=="1998-07-09"] )
	set(df,tmp,"PosRNA",NA)		
	# remove is.na(PosRNA) and !is.na(RNA)
	df		<- subset(df, !is.na(PosRNA) & !is.na(RNA))
	if(verbose)		cat(paste("\nnumber of entries with !is.na(PosRNA) & !is.na(RNA), n=",nrow(df)))	
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
	if(verbose)		print(df[tmp,])
	set(df, tmp, 'RNA', NA_real_)
	df	<- subset(df, !is.na(RNA))
	#
	#	remove RNA values < 400 if undetectable and before ART start 
	#	
	tmp	<- df[, which(Undetectable=='Yes' & RNA<RNA.min.b4T & (is.na(AnyT_T1) | difftime(AnyT_T1, PosRNA, units='days')>0) )]
	if(verbose)		cat(paste("\nremove RNA values < ",RNA.min.b4T," if undetectable and before ART start, n=",length(tmp)))
	if(verbose)		print(df[tmp,][, Patient])
	set(df, tmp, 'RNA', NA_real_)
	df	<- subset(df, !is.na(RNA))
	#
	#	remove RNA values < 400 if detectable and before ART start 
	#	
	tmp	<- df[, which(Undetectable=='No' & RNA<RNA.min.b4T & (is.na(AnyT_T1) | difftime(AnyT_T1, PosRNA, units='days')>0) )]
	if(verbose)		cat(paste("\nremove RNA values < 400 if detectable and before ART start, n=",length(tmp)))
	if(verbose)		print(df[tmp,][, Patient])
	set(df, tmp, 'RNA', NA_real_)
	df	<- subset(df, !is.na(RNA))
	#
	#	remove undetectable & 1e3<RNA<1e4 after treatment start to RNA
	#	
	tmp		<- which( df[, Undetectable=="Yes" & RNA>RNA.stdvl.udetect.aTS & RNA<1e4 & difftime(PosRNA, AnyT_T1, units='days')>=0] )
	if(verbose)		cat(paste("\nremove Undetectable=='Yes' and 1e3<RNA<1e4 and PosRNA>AnyT_T1', n=",length(tmp)))
	if(verbose)		print(df[tmp,])
	set(df, tmp, 'RNA', NA_real_)
	df	<- subset(df, !is.na(RNA))	
	#
	#	set undetectable & RNA<1e3 after treatment start to RNA
	#	
	tmp		<- which( df[, Undetectable=="Yes" & RNA<=RNA.stdvl.udetect.aTS & difftime(PosRNA, AnyT_T1, units='days')>=0] )
	if(verbose)		cat(paste("\nsetting Undetectable=='Yes' and RNA<1e3 and PosRNA>AnyT_T1 to Undetectable=='No' and RNA, n=",length(tmp)))
	set(df,tmp,"Undetectable","No")
	#
	#	set undetectable & RNA>1e4 after treatment start to RNA
	#	
	tmp		<- which( df[, Undetectable=="Yes" & RNA>=1e4 & difftime(PosRNA, AnyT_T1, units='days')>=0] )
	if(verbose)		cat(paste("\nsetting Undetectable=='Yes' and RNA>1e4 and PosRNA>AnyT_T1 to Undetectable=='No' and RNA, n=",length(tmp)))
	if(verbose)		print(df[tmp,])
	set(df,tmp,"Undetectable","No")	
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
	#	remove undetectable if RNA<1e3 before ART
	#
	tmp		<- df[, which(Undetectable=='Yes' & RNA<=1e3 & (is.na(AnyT_T1) | difftime(AnyT_T1,PosRNA,units='days')>0 ) ) ]
	if(verbose)		cat(paste("\nremove undetectable RNA with RNA<1e3 before ART start, n=",length(tmp)))
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	#
	#	remove undetectable if RNA<1e3 before ART
	#	
	tmp		<- df[, which(Undetectable=='Yes' & RNA<=1e4 & (is.na(AnyT_T1) | difftime(AnyT_T1,PosRNA,units='days')>0 ) ) ]
	if(verbose)		cat(paste("\nremove undetectable RNA with 1e3<RNA<1e4 before ART start, n=",length(tmp)))
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	if(verbose)		cat(paste("\nnumber of remaining undetectable RNA, n=",df[,length(which(Undetectable=='Yes'))]))
	#
	#	check remaining undetectable manually
	#	
	if(verbose)		cat(paste("\nCheck these measurements manually -- keeping all"))
	print(subset(df, Undetectable=='Yes'))
	tmp		<- df[, which(Undetectable=='Yes')]
	set(df,tmp,"Undetectable","No")	
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
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, "RNA", 5e6)
	#
	tmp		<- which(df[, Patient=="M27377" & PosRNA>"2007-08-11"])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units -- patient M27377 after 2007-08-11, n=",length(tmp),"DIV by 1e3"))
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-3)
	#
	tmp		<- which(df[, Patient%in%c("M13134","M18385","M18767","M20308","M35814","M35852","M36515","M41877") & RNA>1e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M13134  M18385  M18767  M20308  M35814  M35852  M36515  M41877, n=",length(tmp),"DIV by 10"))
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M38031") & RNA>1e5])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e5  -- M38031, n=",length(tmp),"DIV by 10"))
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M38036","M33131","M33668","M34200","M34302","M20350") & RNA>1e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M38036  M33131  M33668  M34200  M34302 M20350, n=",length(tmp),"DIV by 100"))
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-2)
	#
	tmp		<- which(df[, Patient%in%c("M11995","M13266","M14486","M15621","M16588") & RNA>5e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 5e6  -- M11995 M13266 M14486 M15621 M16588, n=",length(tmp),"DIV by 10"))
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)		
	#
	tmp		<- which(df[, Patient%in%c("M16570") & RNA>2e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 2e6  -- M16570, n=",length(tmp),"DIV by 10"))
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M17655","M37746","M30733","M27377") & RNA>2e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 2e6  -- M17655 M37746 M30733 M27377, n=",length(tmp),"NA"))
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, "RNA", NA)
	#
	tmp		<- which(df[, Patient%in%c("M17554") & RNA>1e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M17554, n=",length(tmp),"DIV by 100"))
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-2)
	#
	tmp		<- which(df[, Patient%in%c("M10607","M10969","M11428","M28707","M31388","M31455","M32401","M32877","M33406","M33839","M33918","M34062","M35280","M30788") & RNA>=1e6])
	if(verbose) cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M10607 M10969 M11428 M28707 M31388 M31455 M32401 M32877 M33406 M33839 M33918 M34062 M35280  M30788, n=",length(tmp),"DIV by 10"))
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	#	take geom mean on remaining duplicate entries
	#
	setkey(df, Patient, PosRNA)
	df		<- df[, list(RNA=as.integer(round(exp(mean(log(RNA))))), Source=Source[1], Undetectable=Undetectable[1], AnyT_T1=AnyT_T1[1] ), by=c('Patient','PosRNA')]	
	# double check entries manually in range >5e6		
	tmp		<- merge(unique(subset(df, RNA>5e6 & PosRNA>AnyT_T1, Patient)), df, by="Patient")
	tmp		<- which(df[, Patient%in%c("M12736","M27885","M14799","M12612") & RNA>5e6])
	if(verbose)		print(df[tmp,], n=300)
	set(df, tmp, "RNA", NA_real_)	
	# double check entries manually in range >5e6		
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
	#tmp		<- merge(subset(df, Undetectable=="Yes" & RNA>1e4, c(Patient,PosRNA)), df.treat, by="Patient")
	#tmp		<- tmp[, 	{
	#			z<- which( difftime(StopTime,PosRNA,units="days")>0 )
	#			z<- z[1]
	#			list(PosRNA=PosRNA[z], StartTime=StartTime[z], StopTime=StopTime[z], TrI=TrI[z], TrCh.failure=TrCh.failure[z], TrCh.adherence=TrCh.adherence[z], TrCh.patrel=TrCh.patrel[z])				
	#		}, by="Patient"]
	#setnames(tmp, "PosRNA","PosRNAh")
	#tmp		<- merge(tmp, df,by="Patient")		
	#tmp[, print(data.table(Patient,RNA,PosRNA,Undetectable,StartTime,StopTime,TrI,TrCh.failure,TrCh.adherence,TrCh.patrel)), by="Patient"]
	#
	#tmp		<- which(df[, Patient=="M30584" & PosRNA=="2004-10-14"])
	#if(verbose) cat(paste("\nset  M30584 2004-10-14 to NA, n=",length(tmp)))
	#set(df, tmp, "RNA", NA)
	#tmp		<- which(df[, Patient=="M26369" & PosRNA=="2003-10-13"])[1]
	#if(verbose) cat(paste("\nset  M26369 2003-10-13 to NA, n=",length(tmp)))
	#set(df, tmp, "RNA", NA)
	#tmp		<- which(df[, Patient=="M16566" & PosRNA=="2000-11-16"])
	#if(verbose) cat(paste("\nset  M16566 2000-11-16 to NA, n=",length(tmp)))
	#set(df, tmp, "RNA", NA)
	#tmp		<- which(df[, Patient=="M12818" & PosRNA=="1996-09-18" & RNA>=1e4])
	#if(verbose) cat(paste("\nset  M12818 1996-09-18 to NA, n=",length(tmp)))
	#set(df, tmp, "RNA", NA)		
	#
	#	set Undetectable=='Yes' & RNA>1e4
	#
	#tmp		<- which(df[, Undetectable=="Yes" & RNA>=1e4])
	#if(verbose) cat(paste("\nsetting Undetectable=='Yes' & RNA>1e4 to Undetectable=='No', n=",length(tmp)))
	#set(df, tmp, "Undetectable", "No")		
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
Viro.160227<- function()		
{
	#	CONSTANTS
	RNA.min					<- 50	#400	#seems to be standard value
	RNA.min.b4T				<- 400
	RNA.stdvl.udetect.aTS	<- 1e3
	RNA.max					<- 5e6
	lRNA.min.infectious		<- log10(1e4)
	lRNA.min.early			<- log10(1e5)
	lRNA.max.b4early		<- log10(2e4)
	verbose					<- 1	
	DB.locktime				<- as.Date("2015-05-01")
	
	#	INPUT FILE NAMES
	dir.name		<- '~/Dropbox (SPH Imperial College)/2015_ATHENA_May_Update'
	file.treatment	<- paste(dir.name,"ATHENA_1502_All_Regimens.R",sep='/')
	file			<- file.path(dir.name, 'ATHENA_1502_All_Viro.csv')		
	
	#	load ART data
	load(file.treatment)
	df.treat			<- subset(df, select=c(Patient, StartTime, StopTime, AnyT_T1, AnyT_T1_Acc, StartTime_Acc, StopTime_Acc, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel))
	#
	#	read VIROLOGY csv data file and preprocess
	#
	df				<- read.csv(file, stringsAsFactors=FALSE)
	df				<- hivc.db.reset.Date(df, col="DateRNA", NA.time=c("","1911-01-01","1911-11-11","1923-06-24"), date.format="%Y-%m-%d")
	df				<- data.table(df, key="Patient")
	set(df, NULL, 'Undetectable', df[, factor(Undetectable, levels=c(0,1,2),labels=c("No","Yes","LargerThan"))])
	setnames(df, "DateRNA","PosRNA")	
	#
	#	Manual Corrections checked with Ard
	#
	tmp		<-  df[, which(Patient=="M32210" & PosRNA=="2010-05-10")]
	if(length(tmp))
	{
		cat('\nFix M32210 2010-05-10')
		set(df,tmp,"PosRNA", 50.)
		set(df,tmp,"Undetectable", 'No')		
	}
	tmp		<- df[, which( Patient=="M36914" & PosRNA=="2010-04-12")]
	{
		cat('\nDelete M36914 2010-04-12')		
		set(df,tmp,"PosRNA", NA_real_)	
	}
	set(df, df[, which(Patient=='M12818' & PosRNA=='1996-09-18')], 'RNA', 400)
	set(df, df[, which(Patient=='M32210' & PosRNA=='2010-10-05')], 'RNA', 50)
	set(df, df[, which(Patient=='M36914' & PosRNA=='2010-12-04')], 'RNA', 0)	
	df				<- subset(df, !is.na(RNA) & !is.na(PosRNA) & !is.na(RNA))
	cat("\nnumber of entries read, n=",nrow(df))	
	#
	#	if duplicate per day, keep those with higher Source variable
	#	
	df[, DUMMY:= seq_len(nrow(df))]
	tmp		<- df[, list(N=list(length(RNA))), by=c('Patient','PosRNA')]	
	tmp		<- merge(subset(tmp, N>1), df, by=c('Patient','PosRNA'))
	tmp		<- tmp[, list(IDX=DUMMY[which(Source!=max(Source))]), by=c('Patient','PosRNA')]
	cat('\nremoving day duplicates with lower Source, n=', nrow(tmp))
	set(df, tmp[,IDX], 'PosRNA', NA_real_)
	set(df, NULL, 'DUMMY', NULL)
	df		<- subset(df, !is.na(PosRNA))	
	#
	#	Handle duplicates with Undetectable Yes/No: rm undetectables='yes'
	#
	df[, DUMMY:= seq_len(nrow(df))]
	tmp		<- df[, list(N=length(RNA), ND=length(unique(Undetectable)) ), by=c('Patient','PosRNA')]
	tmp		<- merge(subset(tmp, N>1 & ND>1), df, by=c('Patient','PosRNA'))
	tmp		<- subset(tmp, Undetectable=='Yes')
	cat("\nFound duplicate entries, one of which is undetectable. Remove Undetectable one, n=",nrow(tmp))
	set(df, tmp[, DUMMY], 'PosRNA', NA_real_) 
	set(df, NULL, 'DUMMY', NULL)
	df		<- subset(df, !is.na(PosRNA))	
	#	
	#	remove identical entries
	#
	setkey(df, Patient, PosRNA, RNA)
	df		<- unique(df)
	#
	#	checking manually NegT>PosRNA
	#
	cat(paste("\nset entry Patient=='M38400' & as.character(PosRNA)=='2005-11-24' to NA -- seems unlikely"))
	#print( subset(df, Patient=="M38400") )
	tmp		<- df[, which( Patient=="M38400" & PosRNA=="2005-11-24")]
	cat(paste("\nsetting number of entries to NA, n=",length(tmp)))
	set(df,tmp,"PosRNA",NA)		
	cat(paste("\nset entry Patient=='M36146' & as.character(PosRNA)=='2005-11-19' to NA -- seems unlikely"))
	#print( subset(df, Patient=="M36146") )
	tmp		<- df[, which( Patient=="M36146" & PosRNA=="2005-11-19")] 
	cat(paste("\nsetting number of entries to NA, n=",length(tmp)))
	set(df,tmp,"PosRNA",NA)
	#print( subset(df, Patient=="M12735") )
	tmp		<- df[, which( Patient=="M12735" & PosRNA=="1986-04-29")] 
	set(df,tmp,"PosRNA",NA)		
	#print( subset(df, Patient=="M29967") )
	tmp		<- df[, which( Patient=="M29967" & PosRNA=="1998-07-09")]
	set(df,tmp,"PosRNA",NA)		
	df		<- subset(df, !is.na(PosRNA) & !is.na(RNA))		
	#
	#	combine Undetectable=="LargerThan" with Undetectable=="No"
	#
	tmp		<- which(df[, Undetectable=="LargerThan"])
	cat(paste("\nsetting Undetectable=='LargerThan' to Undetectable=='No', n=",length(tmp)))
	set(df, tmp,"Undetectable", "No")
	stopifnot(!any(df[,is.na(Undetectable)]))		
	set(df, NULL, "Undetectable",factor(as.numeric(df[,Undetectable]), levels=c(1,2),labels=c("No","Yes")))
	#
	#	remove RNA values that are probably real but have wrong units
	#
	df		<- merge( df, unique(subset(df.treat, select=c(Patient, AnyT_T1))), by='Patient', all.x=1 )
	tmp		<- df[, which((is.na(AnyT_T1) | difftime(AnyT_T1, PosRNA, units='days')>=0) & ((RNA%%1)!=0) ) ]
	cat(paste("\nremove RNA values with .XXX before ART start, n=",length(tmp)))
	#print(df[tmp,])
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	#
	#	remove RNA values < 400 if undetectable and before ART start 
	#	
	tmp		<- df[, which(Undetectable=='Yes' & RNA<RNA.min.b4T & (is.na(AnyT_T1) | difftime(AnyT_T1, PosRNA, units='days')>0) )]
	cat(paste("\nremove RNA values < ",RNA.min.b4T," if undetectable and before ART start, n=",length(tmp)))
	#print(df[tmp,][, Patient])
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	#
	#	remove RNA values < 400 if detectable and before ART start 
	#	
	tmp		<- df[, which(Undetectable=='No' & RNA<RNA.min.b4T & (is.na(AnyT_T1) | difftime(AnyT_T1, PosRNA, units='days')>0) )]
	cat(paste("\nremove RNA values < 400 if detectable and before ART start, n=",length(tmp)))
	#print(df[tmp,][, Patient])
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	#	
	#	remove undetectable & 1e3<RNA<1e4 after treatment start to RNA
	#	
	tmp		<- which( df[, Undetectable=="Yes" & RNA>RNA.stdvl.udetect.aTS & RNA<1e4 & difftime(PosRNA, AnyT_T1, units='days')>=0] )
	cat(paste("\nremove Undetectable=='Yes' and 1e3<RNA<1e4 and PosRNA>AnyT_T1', n=",length(tmp)))
	#print(df[tmp,])
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))	
	#
	#	set undetectable & RNA<1e3 after treatment start to RNA
	#	
	tmp		<- which( df[, Undetectable=="Yes" & RNA<=RNA.stdvl.udetect.aTS & difftime(PosRNA, AnyT_T1, units='days')>=0] )
	cat(paste("\nsetting Undetectable=='Yes' and RNA<1e3 and PosRNA>AnyT_T1 to Undetectable=='No' and RNA, n=",length(tmp)))
	set(df,tmp,"Undetectable","No")
	#
	#	set undetectable & RNA>1e4 after treatment start to RNA
	#	
	tmp		<- which( df[, Undetectable=="Yes" & RNA>=1e4 & difftime(PosRNA, AnyT_T1, units='days')>=0] )
	cat(paste("\nsetting Undetectable=='Yes' and RNA>1e4 and PosRNA>AnyT_T1 to Undetectable=='No' and RNA, n=",length(tmp)))
	#print(df[tmp,])
	set(df,tmp,"Undetectable","No")	
	#
	#	remove undetectable RNA (before treatment start) that are before the first detectable RNA
	#
	tmp		<- df[, list(select= !all(Undetectable=='Yes')), by='Patient']
	cat(paste("\nPatients with all Undetectable RNA, DELETE n=",nrow(subset(tmp, !select))))
	df		<- merge( df, subset(tmp, select, Patient), by='Patient' )
	df		<- df[,  {
				tmp<- seq.int(which(Undetectable=='No')[1], length(Undetectable)) 
				lapply(.SD,'[',tmp)
			},by='Patient']
	cat(paste("\nnumber of entries after removing undetectable RNA before any detectable RNA, n=",nrow(df)))
	cat(paste("\nnumber of remaining undetectable RNA, n=",df[,length(which(Undetectable=='Yes'))]))
	#M38913, M38722, M38411, M37285, M37211, M35710, M33521, M31924, M28495, M15900
	# 	remove further suspicious entries
	tmp		<- df[, which(Patient=='M38913' & PosRNA==as.Date("2010-07-01"))]	
	tmp		<- c(tmp, df[, which(Patient=='M38722' & PosRNA==as.Date("2010-05-11"))])
	tmp		<- c(tmp, df[, which(Patient=='M37211' & PosRNA==as.Date("2009-08-18"))])
	tmp		<- c(tmp, df[, which(Patient=='M33521' & PosRNA==as.Date("2006-12-18"))])
	tmp		<- c(tmp, df[, which(Patient=='M15900' & PosRNA==as.Date("1997-01-13"))])
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	#
	#	remove undetectable if RNA<1e3 before ART
	#
	tmp		<- df[, which(Undetectable=='Yes' & RNA<=1e3 & (is.na(AnyT_T1) | difftime(AnyT_T1,PosRNA,units='days')>0 ) ) ]
	cat(paste("\nremove undetectable RNA with RNA<1e3 before ART start, n=",length(tmp)))
	#print(df[tmp,], n=300)
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	#
	#	remove undetectable if RNA<1e3 before ART
	#	
	tmp		<- df[, which(Undetectable=='Yes' & RNA<=1e4 & (is.na(AnyT_T1) | difftime(AnyT_T1,PosRNA,units='days')>0 ) ) ]
	cat(paste("\nremove undetectable RNA with 1e3<RNA<1e4 before ART start, n=",length(tmp)))
	#print(df[tmp,], n=300)
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	cat(paste("\nnumber of remaining undetectable RNA, n=",df[,length(which(Undetectable=='Yes'))]))
	#
	#	check remaining undetectable manually
	#	
	cat(paste("\nCheck these measurements manually -- keeping all"))
	print(subset(df, Undetectable=='Yes'))
	tmp		<- df[, which(Undetectable=='Yes')]
	set(df,tmp,"Undetectable","No")	
	#
	#	set RNA<RNA.min to RNA.min
	#		
	tmp<- which( df[, RNA<RNA.min] )
	cat(paste("\nsetting RNA<RNA.min to RNA.min=",RNA.min,", n=",length(tmp)))		
	set(df,tmp,"RNA",RNA.min)
	#
	#	wrong units ? adjust manually
	#				
	tmp		<- which(df[, Patient%in%c("M11314","M11331","M40782","M14759") & RNA>5e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 5e6  -- M11314  M11331  M40782  M14759, n=",length(tmp),"SET to 5e6"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", 5e6)
	#
	tmp		<- which(df[, Patient=="M27377" & PosRNA>"2007-08-11"])
	cat(paste("\nnumber of entries with likely wrong RNA units -- patient M27377 after 2007-08-11, n=",length(tmp),"DIV by 1e3"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-3)
	#
	tmp		<- which(df[, Patient%in%c("M13134","M18385","M18767","M20308","M35814","M35852","M36515","M41877") & RNA>1e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M13134  M18385  M18767  M20308  M35814  M35852  M36515  M41877, n=",length(tmp),"DIV by 10"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M38031") & RNA>1e5])
	cat(paste("\nnumber of entries with likely wrong RNA units > 1e5  -- M38031, n=",length(tmp),"DIV by 10"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M38036","M33131","M33668","M34200","M34302","M20350") & RNA>1e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M38036  M33131  M33668  M34200  M34302 M20350, n=",length(tmp),"DIV by 100"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-2)
	#
	tmp		<- which(df[, Patient%in%c("M11995","M13266","M14486","M15621","M16588") & RNA>5e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 5e6  -- M11995 M13266 M14486 M15621 M16588, n=",length(tmp),"DIV by 10"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)		
	#
	tmp		<- which(df[, Patient%in%c("M16570") & RNA>2e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 2e6  -- M16570, n=",length(tmp),"DIV by 10"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M17655","M37746","M30733","M27377") & RNA>2e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 2e6  -- M17655 M37746 M30733 M27377, n=",length(tmp),"NA"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", NA)
	#
	tmp		<- which(df[, Patient%in%c("M17554") & RNA>1e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M17554, n=",length(tmp),"DIV by 100"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-2)
	#
	tmp		<- which(df[, Patient%in%c("M10607","M10969","M11428","M28707","M31388","M31455","M32401","M32877","M33406","M33839","M33918","M34062","M35280","M30788") & RNA>=1e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M10607 M10969 M11428 M28707 M31388 M31455 M32401 M32877 M33406 M33839 M33918 M34062 M35280  M30788, n=",length(tmp),"DIV by 10"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	#	take geom mean on remaining duplicate entries
	#
	setkey(df, Patient, PosRNA)
	df		<- df[, list(RNA=as.integer(round(exp(mean(log(RNA))))), Source=Source[1], Undetectable=Undetectable[1], AnyT_T1=AnyT_T1[1] ), by=c('Patient','PosRNA')]
	#
	# 	double check entries manually in range >5e6
	#
	tmp		<- merge(unique(subset(df, RNA>5e6 & PosRNA>AnyT_T1, Patient)), df, by="Patient")
	tmp		<- df[, which(Patient%in%c("M12736","M27885","M14799","M12612") & RNA>5e6)]
	print(df[tmp,], n=300)
	set(df, tmp, "RNA", NA_integer_)	
	# double check entries manually in range >5e6		
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
	cat(paste("\nset  M12612 1996-06-13 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	tmp		<- which(df[, Patient=="M17044" & PosRNA=="2004-07-19"])
	cat(paste("\nset  M17044 2004-07-19 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	#
	#	set RNA>RNA.max to RNA.max
	#		
	tmp		<- which(df[, RNA>RNA.max])
	cat(paste("\nsetting RNA>RNA.max to RNA.max, n=",length(tmp)))
	set(df, tmp, "RNA", RNA.max)		
	#
	df		<- subset(df,!is.na(RNA), select=c(Patient, PosRNA, RNA))
	cat(paste("\nentries with !is.na(RNA), n=",nrow(df)))	
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
			}, by=Patient]
	tmp					<- subset(df.treat,select=c(Patient,AnyT_T1))
	setkey(tmp,Patient)		 
	df					<- merge(df, unique(tmp), all.x=1, by="Patient")
	set(df, which(df[,is.na(AnyT_T1)]),"AnyT_T1", DB.locktime)		#only need AnyT_T1 for internal calculations
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
	df	<- merge(subset(df,select=c(Patient, PosRNA, RNA, lRNA, PoslRNA_T1, lRNA_T1)), subset(tmp,select=c(Patient,lRNA.i,lRNA.early)), all.x=1, by="Patient")
	
	file		<- paste(substr(file, 1, nchar(file)-3),'R',sep='')	
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)	
}
######################################################################################
Viro.161027<- function()		
{
	#	CONSTANTS
	RNA.min					<- 50	#400	#seems to be standard value
	RNA.min.b4T				<- 400
	RNA.stdvl.udetect.aTS	<- 1e3
	RNA.max					<- 5e6
	lRNA.min.infectious		<- log10(1e4)
	lRNA.min.early			<- log10(1e5)
	lRNA.max.b4early		<- log10(2e4)
	verbose					<- 1	
	DB.locktime				<- as.Date("2016-06-01")
	
	#	INPUT FILE NAMES
	file.treatment	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_ART.rda'
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/original_latest/SHM_1602_161102_OR_ALL_RNA.csv'		
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_RNA.rda'
	
	#	load ART data
	load(file.treatment)
	df.treat			<- subset(df, select=c(Patient, StartTime, StopTime, AnyT_T1, AnyT_T1_Acc, StartTime_Acc, StopTime_Acc, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel))
	#
	#	read VIROLOGY csv data file and preprocess
	#
	df				<- read.csv(file, stringsAsFactors=FALSE)	
	df				<- hivc.db.reset.Date(df, col="RNA_D", NA.time=c("","01/01/1911","11/11/1911","24/06/1923"), date.format="%d/%m/%Y")
	df				<- data.table(df, key="Patient")
	set(df, NULL, 'RNA_UDT', df[, factor(RNA_UDT, levels=c(0,1,2),labels=c("No","Yes","LargerThan"))])
	setnames(df, c("RNA_D","RNA_UDT"),c("PosRNA","Undetectable"))	
	#
	#	Manual Corrections checked with Ard
	#
	tmp		<-  df[, which(Patient=="M32210" & PosRNA=="2010-05-10")]
	if(length(tmp))
	{
		cat('\nFix M32210 2010-05-10')
		set(df,tmp,"PosRNA", 50.)
		set(df,tmp,"Undetectable", 'No')		
	}
	tmp		<- df[, which( Patient=="M36914" & PosRNA=="2010-04-12")]
	{
		cat('\nDelete M36914 2010-04-12')		
		set(df,tmp,"PosRNA", NA_real_)	
	}
	set(df, df[, which(Patient=='M12818' & PosRNA=='1996-09-18')], 'RNA', 400)
	set(df, df[, which(Patient=='M32210' & PosRNA=='2010-10-05')], 'RNA', 50)
	set(df, df[, which(Patient=='M36914' & PosRNA=='2010-12-04')], 'RNA', 0)	
	df				<- subset(df, !is.na(RNA) & !is.na(PosRNA) & !is.na(RNA))
	cat("\nnumber of entries read, n=",nrow(df))	
	#
	#	if duplicate per day, keep those with higher Source variable
	#	
	df[, DUMMY:= seq_len(nrow(df))]
	tmp		<- df[, list(N=list(length(RNA))), by=c('Patient','PosRNA')]	
	tmp		<- merge(subset(tmp, N>1), df, by=c('Patient','PosRNA'))
	tmp		<- tmp[, list(IDX=DUMMY[which(Source!=max(Source))]), by=c('Patient','PosRNA')]
	cat('\nremoving day duplicates with lower Source, n=', nrow(tmp))
	set(df, tmp[,IDX], 'PosRNA', NA_real_)
	set(df, NULL, 'DUMMY', NULL)
	df		<- subset(df, !is.na(PosRNA))	
	#
	#	Handle duplicates with Undetectable Yes/No: rm undetectables='yes'
	#
	df[, DUMMY:= seq_len(nrow(df))]
	tmp		<- df[, list(N=length(RNA), ND=length(unique(Undetectable)) ), by=c('Patient','PosRNA')]
	tmp		<- merge(subset(tmp, N>1 & ND>1), df, by=c('Patient','PosRNA'))
	tmp		<- subset(tmp, Undetectable=='Yes')
	cat("\nFound duplicate entries, one of which is undetectable. Remove Undetectable one, n=",nrow(tmp))
	set(df, tmp[, DUMMY], 'PosRNA', NA_real_) 
	set(df, NULL, 'DUMMY', NULL)
	df		<- subset(df, !is.na(PosRNA))	
	#	
	#	remove identical entries
	#
	setkey(df, Patient, PosRNA, RNA)
	df		<- unique(df)
	#
	#	checking manually NegT>PosRNA
	#
	cat(paste("\nset entry Patient=='M38400' & as.character(PosRNA)=='2005-11-24' to NA -- seems unlikely"))
	#print( subset(df, Patient=="M38400") )
	tmp		<- df[, which( Patient=="M38400" & PosRNA=="2005-11-24")]
	cat(paste("\nsetting number of entries to NA, n=",length(tmp)))
	set(df,tmp,"PosRNA",NA)		
	cat(paste("\nset entry Patient=='M36146' & as.character(PosRNA)=='2005-11-19' to NA -- seems unlikely"))
	#print( subset(df, Patient=="M36146") )
	tmp		<- df[, which( Patient=="M36146" & PosRNA=="2005-11-19")] 
	cat(paste("\nsetting number of entries to NA, n=",length(tmp)))
	set(df,tmp,"PosRNA",NA)
	#print( subset(df, Patient=="M12735") )
	tmp		<- df[, which( Patient=="M12735" & PosRNA=="1986-04-29")] 
	set(df,tmp,"PosRNA",NA)		
	#print( subset(df, Patient=="M29967") )
	tmp		<- df[, which( Patient=="M29967" & PosRNA=="1998-07-09")]
	set(df,tmp,"PosRNA",NA)		
	df		<- subset(df, !is.na(PosRNA) & !is.na(RNA))		
	#
	#	combine Undetectable=="LargerThan" with Undetectable=="No"
	#
	tmp		<- which(df[, Undetectable=="LargerThan"])
	cat(paste("\nsetting Undetectable=='LargerThan' to Undetectable=='No', n=",length(tmp)))
	set(df, tmp,"Undetectable", "No")
	stopifnot(!any(df[,is.na(Undetectable)]))		
	set(df, NULL, "Undetectable",factor(as.numeric(df[,Undetectable]), levels=c(1,2),labels=c("No","Yes")))
	#
	#	remove RNA values that are probably real but have wrong units
	#
	df		<- merge( df, unique(subset(df.treat, select=c(Patient, AnyT_T1))), by='Patient', all.x=1 )
	tmp		<- df[, which((is.na(AnyT_T1) | difftime(AnyT_T1, PosRNA, units='days')>=0) & ((RNA%%1)!=0) ) ]
	cat(paste("\nremove RNA values with .XXX before ART start, n=",length(tmp)))
	#print(df[tmp,])
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	#
	#	remove RNA values < 400 if undetectable and before ART start 
	#	
	tmp		<- df[, which(Undetectable=='Yes' & RNA<RNA.min.b4T & (is.na(AnyT_T1) | difftime(AnyT_T1, PosRNA, units='days')>0) )]
	cat(paste("\nremove RNA values < ",RNA.min.b4T," if undetectable and before ART start, n=",length(tmp)))
	#print(df[tmp,][, Patient])
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	#
	#	remove RNA values < 400 if detectable and before ART start 
	#	
	tmp		<- df[, which(Undetectable=='No' & RNA<RNA.min.b4T & (is.na(AnyT_T1) | difftime(AnyT_T1, PosRNA, units='days')>0) )]
	cat(paste("\nremove RNA values < 400 if detectable and before ART start, n=",length(tmp)))
	#print(df[tmp,][, Patient])
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	#	
	#	remove undetectable & 1e3<RNA<1e4 after treatment start to RNA
	#	
	tmp		<- which( df[, Undetectable=="Yes" & RNA>RNA.stdvl.udetect.aTS & RNA<1e4 & difftime(PosRNA, AnyT_T1, units='days')>=0] )
	cat(paste("\nremove Undetectable=='Yes' and 1e3<RNA<1e4 and PosRNA>AnyT_T1', n=",length(tmp)))
	#print(df[tmp,])
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))	
	#
	#	set undetectable & RNA<1e3 after treatment start to RNA
	#	
	tmp		<- which( df[, Undetectable=="Yes" & RNA<=RNA.stdvl.udetect.aTS & difftime(PosRNA, AnyT_T1, units='days')>=0] )
	cat(paste("\nsetting Undetectable=='Yes' and RNA<1e3 and PosRNA>AnyT_T1 to Undetectable=='No' and RNA, n=",length(tmp)))
	set(df,tmp,"Undetectable","No")
	#
	#	set undetectable & RNA>1e4 after treatment start to RNA
	#	
	tmp		<- which( df[, Undetectable=="Yes" & RNA>=1e4 & difftime(PosRNA, AnyT_T1, units='days')>=0] )
	cat(paste("\nsetting Undetectable=='Yes' and RNA>1e4 and PosRNA>AnyT_T1 to Undetectable=='No' and RNA, n=",length(tmp)))
	#print(df[tmp,])
	set(df,tmp,"Undetectable","No")	
	#
	#	remove undetectable RNA (before treatment start) that are before the first detectable RNA
	#
	tmp		<- df[, list(select= !all(Undetectable=='Yes')), by='Patient']
	cat(paste("\nPatients with all Undetectable RNA, DELETE n=",nrow(subset(tmp, !select))))
	df		<- merge( df, subset(tmp, select, Patient), by='Patient' )
	df		<- df[,  {
				tmp<- seq.int(which(Undetectable=='No')[1], length(Undetectable)) 
				lapply(.SD,'[',tmp)
			},by='Patient']
	cat(paste("\nnumber of entries after removing undetectable RNA before any detectable RNA, n=",nrow(df)))
	cat(paste("\nnumber of remaining undetectable RNA, n=",df[,length(which(Undetectable=='Yes'))]))
	#M38913, M38722, M38411, M37285, M37211, M35710, M33521, M31924, M28495, M15900
	# 	remove further suspicious entries
	tmp		<- df[, which(Patient=='M38913' & PosRNA==as.Date("2010-07-01"))]	
	tmp		<- c(tmp, df[, which(Patient=='M38722' & PosRNA==as.Date("2010-05-11"))])
	tmp		<- c(tmp, df[, which(Patient=='M37211' & PosRNA==as.Date("2009-08-18"))])
	tmp		<- c(tmp, df[, which(Patient=='M33521' & PosRNA==as.Date("2006-12-18"))])
	tmp		<- c(tmp, df[, which(Patient=='M15900' & PosRNA==as.Date("1997-01-13"))])
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	#
	#	remove undetectable if RNA<1e3 before ART
	#
	tmp		<- df[, which(Undetectable=='Yes' & RNA<=1e3 & (is.na(AnyT_T1) | difftime(AnyT_T1,PosRNA,units='days')>0 ) ) ]
	cat(paste("\nremove undetectable RNA with RNA<1e3 before ART start, n=",length(tmp)))
	#print(df[tmp,], n=300)
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	#
	#	remove undetectable if RNA<1e3 before ART
	#	
	tmp		<- df[, which(Undetectable=='Yes' & RNA<=1e4 & (is.na(AnyT_T1) | difftime(AnyT_T1,PosRNA,units='days')>0 ) ) ]
	cat(paste("\nremove undetectable RNA with 1e3<RNA<1e4 before ART start, n=",length(tmp)))
	#print(df[tmp,], n=300)
	set(df, tmp, 'RNA', NA_real_)
	df		<- subset(df, !is.na(RNA))
	cat(paste("\nnumber of remaining undetectable RNA, n=",df[,length(which(Undetectable=='Yes'))]))
	#
	#	check remaining undetectable manually
	#	
	cat(paste("\nCheck these measurements manually -- keeping all"))
	print(subset(df, Undetectable=='Yes'))
	tmp		<- df[, which(Undetectable=='Yes')]
	set(df,tmp,"Undetectable","No")	
	#
	#	set RNA<RNA.min to RNA.min
	#		
	tmp<- which( df[, RNA<RNA.min] )
	cat(paste("\nsetting RNA<RNA.min to RNA.min=",RNA.min,", n=",length(tmp)))		
	set(df,tmp,"RNA",RNA.min)
	#
	#	wrong units ? adjust manually
	#				
	tmp		<- which(df[, Patient%in%c("M11314","M11331","M40782","M14759") & RNA>5e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 5e6  -- M11314  M11331  M40782  M14759, n=",length(tmp),"SET to 5e6"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", 5e6)
	#
	tmp		<- which(df[, Patient=="M27377" & PosRNA>"2007-08-11"])
	cat(paste("\nnumber of entries with likely wrong RNA units -- patient M27377 after 2007-08-11, n=",length(tmp),"DIV by 1e3"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-3)
	#
	tmp		<- which(df[, Patient%in%c("M13134","M18385","M18767","M20308","M35814","M35852","M36515","M41877") & RNA>1e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M13134  M18385  M18767  M20308  M35814  M35852  M36515  M41877, n=",length(tmp),"DIV by 10"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M38031") & RNA>1e5])
	cat(paste("\nnumber of entries with likely wrong RNA units > 1e5  -- M38031, n=",length(tmp),"DIV by 10"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M38036","M33131","M33668","M34200","M34302","M20350") & RNA>1e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M38036  M33131  M33668  M34200  M34302 M20350, n=",length(tmp),"DIV by 100"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-2)
	#
	tmp		<- which(df[, Patient%in%c("M11995","M13266","M14486","M15621","M16588") & RNA>5e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 5e6  -- M11995 M13266 M14486 M15621 M16588, n=",length(tmp),"DIV by 10"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)		
	#
	tmp		<- which(df[, Patient%in%c("M16570") & RNA>2e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 2e6  -- M16570, n=",length(tmp),"DIV by 10"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	tmp		<- which(df[, Patient%in%c("M17655","M37746","M30733","M27377") & RNA>2e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 2e6  -- M17655 M37746 M30733 M27377, n=",length(tmp),"NA"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", NA)
	#
	tmp		<- which(df[, Patient%in%c("M17554") & RNA>1e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M17554, n=",length(tmp),"DIV by 100"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-2)
	#
	tmp		<- which(df[, Patient%in%c("M10607","M10969","M11428","M28707","M31388","M31455","M32401","M32877","M33406","M33839","M33918","M34062","M35280","M30788") & RNA>=1e6])
	cat(paste("\nnumber of entries with likely wrong RNA units > 1e6  -- M10607 M10969 M11428 M28707 M31388 M31455 M32401 M32877 M33406 M33839 M33918 M34062 M35280  M30788, n=",length(tmp),"DIV by 10"))
	#print(df[tmp,], n=300)
	set(df, tmp, "RNA", df[tmp,RNA]*1e-1)
	#
	#	take geom mean on remaining duplicate entries
	#
	setkey(df, Patient, PosRNA)
	df		<- df[, list(RNA=as.integer(round(exp(mean(log(RNA))))), Source=Source[1], Undetectable=Undetectable[1], AnyT_T1=AnyT_T1[1] ), by=c('Patient','PosRNA')]
	#
	# 	double check entries manually in range >5e6
	#
	tmp		<- merge(unique(subset(df, RNA>5e6 & PosRNA>AnyT_T1, Patient)), df, by="Patient")
	tmp		<- df[, which(Patient%in%c("M12736","M27885","M14799","M12612") & RNA>5e6)]
	print(df[tmp,], n=300)
	set(df, tmp, "RNA", NA_integer_)	
	# double check entries manually in range >5e6		
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
	cat(paste("\nset  M12612 1996-06-13 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	tmp		<- which(df[, Patient=="M17044" & PosRNA=="2004-07-19"])
	cat(paste("\nset  M17044 2004-07-19 to NA, n=",length(tmp)))
	set(df, tmp, "RNA", NA)
	#
	#	set RNA>RNA.max to RNA.max
	#		
	tmp		<- which(df[, RNA>RNA.max])
	cat(paste("\nsetting RNA>RNA.max to RNA.max, n=",length(tmp)))
	set(df, tmp, "RNA", RNA.max)		
	#
	df		<- subset(df,!is.na(RNA), select=c(Patient, PosRNA, RNA))
	cat(paste("\nentries with !is.na(RNA), n=",nrow(df)))	
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
			}, by=Patient]
	tmp					<- subset(df.treat,select=c(Patient,AnyT_T1))
	setkey(tmp,Patient)		 
	df					<- merge(df, unique(tmp), all.x=1, by="Patient")
	set(df, which(df[,is.na(AnyT_T1)]),"AnyT_T1", DB.locktime)		#only need AnyT_T1 for internal calculations
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
	df	<- merge(subset(df,select=c(Patient, PosRNA, RNA, lRNA, PoslRNA_T1, lRNA_T1)), subset(tmp,select=c(Patient,lRNA.i,lRNA.early)), all.x=1, by="Patient")
	
	
	if(verbose) cat(paste("\nsave to", outfile))
	save(df, file=outfile)	
}
######################################################################################
Patients.160227<- function()
{
	#	CONSTANTS
	NA.Acute			<- c(NA,9)
	NA.CountryInfection	<- c(NA,"")
	NA.CountryBorn		<- c(NA,"",'XX')
	NA.RegionOrigin		<- c(NA,"",'XX')
	NA.Subtype			<- c(NA,"")
	NA.time				<- c("","1911-01-01","1911-11-11")		
	date.format			<- "%Y-%m-%d"
	#	INPUT
	dir.name			<- "~/Dropbox (SPH Imperial College)/2015_ATHENA_May_Update"
	file				<- file.path(dir.name, "ATHENA_1502_All_Patient.csv")
	outfile				<- gsub('csv','R',file)
	#
	#	read PATIENT csv data file	
	#
	df			<- read.csv(file, stringsAsFactors=FALSE)
	#
	#	reset Dates
	#
	df			<- hivc.db.reset.Date(df, 	col=c("DateBorn","MyDateNeg1","MyDatePos1","DateDied","DateLastContact","DateAIDS","FirstMed","DateInCare","T0"), 
			NA.time=c("","1911-01-01","1911-11-11"), 
			date.format="%Y-%m-%d")
	#
	#	rename columns
	#
	df			<- as.data.table(df)	
	setnames(df, c("MyDateNeg1_Acc","MyDatePos1_Acc","Acute","Transmission","Reg_Last_Reason","ACS_ID"), c("NegT_Acc","PosT_Acc","isAcute","Trm","GGD_LastCh","ACS"))
	#
	#	set factors
	#
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
	set(df, NULL, 'Region_first', df[, factor(Region_first, levels=c('AMSTERDAM','ROTTERDAM','REST'), labels=c('Amst','Rott','Other'))])
	set(df, NULL, 'Region_now', df[, factor(Region_now, levels=c('AMSTERDAM','ROTTERDAM','REST'), labels=c('Amst','Rott','Other'))])
	set(df, NULL, 'GGD_first', 	df[, factor(GGD_first, 	levels=c(111, 		 			706, 				1009, 		 		1106,			1406, 				1906, 
									2006, 		 			2106, 				2209, 				2406, 			2506, 				2707, 
									3109, 					3406, 				3606, 				3906, 			4106, 				4506, 
									4607, 					4810, 				5006, 				5206, 			5406, 				5608, 
									6011, 					6106, 				7206, 				7306), 
							labels=c('Groningen',			'Drenthe',			'IJsselland',		'Twente',		'Gelre_IJssel',		'Hulpverlening_Gelderland_Midden',
									'Rivierenland',		'Nijmegen', 		'Flevoland',		'Utrecht',		'Midden_Nederland',	'Hollands_Noorden',
									'Kennemerland',		'Amsterdam',		'Gooi_Vechtstreek',	'Den_Haag',		'Zuid_Holland_West','Hollands_Midden',
									'Rotterdam_Rijnmond',	'Zuid_Holland_Zuid','Zeeland',			'West_Brabant',	'Hart_voor_Brabant', 'Brabant_Zuidoost',
									'Limburg-Noord',		'Zuid_Limburg',		'Fryslan',			'Zaanstreek_Waterland'))])
	set(df, NULL, 'GGD_now', 	df[, factor(GGD_now, 	levels=c(111, 		 			706, 				1009, 		 		1106,			1406, 				1906, 
									2006, 		 			2106, 				2209, 				2406, 			2506, 				2707, 
									3109, 					3406, 				3606, 				3906, 			4106, 				4506, 
									4607, 					4810, 				5006, 				5206, 			5406, 				5608, 
									6011, 					6106, 				7206, 				7306), 
							labels=c('Groningen',			'Drenthe',			'IJsselland',		'Twente',		'Gelre_IJssel',		'Hulpverlening_Gelderland_Midden',
									'Rivierenland',		'Nijmegen', 		'Flevoland',		'Utrecht',		'Midden_Nederland',	'Hollands_Noorden',
									'Kennemerland',		'Amsterdam',		'Gooi_Vechtstreek',	'Den_Haag',		'Zuid_Holland_West','Hollands_Midden',
									'Rotterdam_Rijnmond',	'Zuid_Holland_Zuid','Zeeland',			'West_Brabant',	'Hart_voor_Brabant', 'Brabant_Zuidoost',
									'Limburg-Noord',		'Zuid_Limburg',		'Fryslan',			'Zaanstreek_Waterland'))])	
	set(df, df[, which(Trm==900)], 'Trm', NA_integer_)
	set(df, NULL, "Trm", factor(df[, Trm], 	levels= c(	100, 	110, 	150, 	200,  202, 		300,  	400,  		450,  		600,  	620, 		800),
					labels= c(	'MSM',	'BI',	'SXCH',	'HET','HETfa',	'IDU',	'BLOOD',	'NEEACC',	'PREG',	'BREAST',	'OTH')))
	set(df, NULL, 'RegionOrigin', df[, factor(RegionOrigin, levels=c('AUS',		'Car',		'EUC',			'EUO',				'EUW',		'Lat',					'NAM',						'NL',	'OAP',				'SSA',					'ZAz'),
							labels=c('Austr_NZ','Caribbean','Central_EU',	'Eastern_EU_stans',	'Western_EU','Latin_South_America',	'North_Africa_Middle_East',	'NL',	'Oceania_Pacific',	'Sub_Saharan_Africa',	'Sout_SouthEast_Asia'))])
	set(df, NULL, 'isDead', df[, factor(is.na(DateDied), levels=c(TRUE,FALSE), labels=c("No","Yes"))] )
	set(df, NULL, "NegT_Acc", factor(df[,NegT_Acc], levels=c(0,1), labels=c("No","Yes")))
	set(df, NULL, "PosT_Acc", factor(df[,PosT_Acc], levels=c(0,1), labels=c("No","Yes")))
	
	#
	#	process Acute fields
	#
	tmp	<- df[, {
				tmp		<- NA_character_
				z		<- Acute_Spec_1%in%c(1L,2L) | Acute_Spec_2%in%c(1L,2L) | Acute_Spec_3%in%c(1L,2L) | Acute_Spec_4%in%c(1L,2L)
				if(!is.na(z) & z)
					tmp	<- 'LAB'
				z		<- Acute_Spec_1%in%c(3L,4L,8L,9L) | Acute_Spec_2%in%c(3L,4L,8L,9L) | Acute_Spec_3%in%c(3L,4L,8L,9L) | Acute_Spec_4%in%c(3L,4L,8L,9L)
				if(is.na(tmp) & !is.na(z) & z)
					tmp	<- 'SYM'
				list(Acute_Spec=tmp)
			}, by='Patient']
	df	<- merge(df, tmp, by='Patient')
	set(df, NULL, 'Acute_Spec', df[, factor(Acute_Spec)])
	setkey(df,Patient)
	str(df)
	#
	#	reset NegT date (checked w Ard)
	#
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
	#df[, table(is.na(isAcute), Acute_Spec)]; subset(df, is.na(isAcute) & !is.na(Acute_Spec))
	#df[, table(isAcute, Acute_Spec, useNA='if')]
	#
	#	reset isAcute
	#
	tmp	<- df[, which(Acute_Spec=='SYM' & is.na(isAcute))]
	cat(paste('\nFound entries with AcuteSpec and is.na(isAcute), n=', length(tmp)))
	set(df, tmp, 'isAcute', 'Maybe') 
	tmp	<- df[, which(Acute_Spec=='LAB' & is.na(isAcute))]
	cat(paste('\nFound entries with AcuteSpec and is.na(isAcute), n=', length(tmp)))
	set(df, tmp, 'isAcute', 'Yes') 	
	tmp	<- df[, which(Acute_Spec=='LAB' & isAcute=='No')]
	cat(paste('\nFound entries with AcuteSpec and isAcute==No, n=', length(tmp)))
	set(df, tmp, 'isAcute', 'Yes') 
	tmp	<- df[, which(Acute_Spec=='SYM' & isAcute=='No')]
	cat(paste('\nFound entries with AcuteSpec and isAcute==No, n=', length(tmp)))
	set(df, tmp, 'isAcute', 'Maybe') 
	#	reset AcuteSpec
	tmp	<- df[, which(is.na(Acute_Spec) & isAcute=='Yes')]
	cat(paste('\nFound entries with is.na(Acute_Spec) and isAcute==Yes, n=', length(tmp)))
	set(df, tmp, 'Acute_Spec', 'CLIN')
	tmp	<- df[, which(Acute_Spec=='SYM' & isAcute=='Yes')]
	cat(paste('\nFound entries with Acute_Spec==SYM and isAcute==Yes, n=', length(tmp)))
	set(df, tmp, 'Acute_Spec', 'CLIN') 	
	tmp	<- df[, which(is.na(Acute_Spec) & isAcute=='Maybe')]
	cat(paste('\nFound entries with is.na(Acute_Spec) and isAcute==Maybe, n=', length(tmp)))
	set(df, tmp, 'Acute_Spec', 'SYM') 
	#
	#
	cat("\nsave to", outfile)	
	save(df, file=outfile)	
}
######################################################################################
CD4.160227<- function()
{
	#	TODO check with Ard that max source should be kept?
	verbose			<- 1
	
	#	INPUT FILE NAMES
	dir.name		<- '~/Dropbox (SPH Imperial College)/2015_ATHENA_May_Update'	
	file			<- paste(dir.name,"ATHENA_1502_All_Immu.csv",sep='/')
	#
	#	read CD4 csv data file and preprocess
	#
	df				<- read.csv(file, stringsAsFactors=FALSE)		
	df				<- hivc.db.reset.Date(df, col="DateImm", NA.time=c("","1911-01-01","1911-11-11","1923-06-24"), date.format="%Y-%m-%d")
	df				<- data.table(df, key="Patient")
	setnames(df, "DateImm","PosCD4")
	#
	#	data corrections from Ard
	#
	tmp		<- df[, which( Patient=="M11392" & PosCD4=="2000-10-30" & CD4A==3001)]
	if(verbose) cat(paste("\nnumber of entries with corrected CD4 units Patient=='M11392' & CD4A>3001, n=",length(tmp),"double entry, remove"))
	set(df, tmp, "CD4A", NA)	
	tmp		<- df[, which( Patient=="M26334" & PosCD4=="2008-09-11" & CD4A==2230)]
	if(verbose) cat(paste("\nnumber of entries with corrected CD4 units Patient==M26334 & PosCD4==2008-09-11 & CD4A==2230, n=",length(tmp),"double entry, remove"))
	set(df, tmp, "CD4A", NA)		
	tmp		<- df[, which( 	Patient=='M26293' & PosCD4=='2003-09-12' & CD4A==800)]
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M26293' & PosCD4=='2003-09-12' & CD4A==800, n=",length(tmp),"set to 16/9/2003"))
	set(df, tmp, "PosCD4", as.Date('2003-09-16'))
	tmp		<- df[, which( 	Patient=='M11233' & PosCD4=='2001-01-10' & CD4A==24)]
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M11233' & PosCD4=='2001-01-10' & CD4A==24, n=",length(tmp),"set to 25/4/1996"))
	set(df, tmp, "PosCD4", as.Date('1996-04-25'))
	tmp		<- df[, which( 	Patient=='M13124' & PosCD4=='1998-08-07')]
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M13124' & PosCD4=='1998-08-07', n=",length(tmp),"set to 7/9/1998"))
	set(df, tmp, "PosCD4", as.Date('1998-09-07'))	
	tmp		<- df[, which( 	Patient=='M13124' & PosCD4=='2001-09-14')]
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M13124' & PosCD4=='2001-09-14', n=",length(tmp),"set to 14/9/2000"))
	set(df, tmp, "PosCD4", as.Date('2000-09-14'))
	tmp		<- df[, which(	Patient=="M10544" & PosCD4=="2003-02-17" & CD4A==850	|
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
			)]
	cat(paste("\nnumber of entries with incorrect CD4  should be 45, n=",length(tmp),"set to NA"))	
	set(df, tmp, "CD4A", NA)
	df		<- subset(df, !is.na(CD4A))
	cat(paste("\nnumber of entries read, n=",nrow(df)))	
	#
	#	if duplicate per day, keep those with higher Source variable
	#	
	df[, DUMMY:= seq_len(nrow(df))]
	tmp		<- df[, list(N=list(length(CD4A))), by=c('Patient','PosCD4')]	
	tmp		<- merge(subset(tmp, N>1), df, by=c('Patient','PosCD4'))
	tmp		<- tmp[, list(IDX=DUMMY[which(Source!=max(Source))]), by=c('Patient','PosCD4')]
	cat('\nremoving day duplicates with lower Source, n=', nrow(tmp))
	set(df, tmp[,IDX], 'PosCD4', NA_real_)
	set(df, NULL, 'DUMMY', NULL)
	df		<- subset(df, !is.na(PosCD4))	
	#
	#	remove identical entries
	#
	setkey(df, Patient, PosCD4, CD4A)
	df		<- unique(df)
	#
	#	adjust likely wrong units > 20000
	#
	tmp		<- df[, which( CD4A>20000)]
	cat(paste("\nnumber of entries with likely wrong CD4 units > 20000, n=",length(tmp),"DIVIDE BY 1e3"))
	#print(df[tmp,])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-3,d=0))
	#	adjust likely wrong units > 10000
	cat(paste("\npatient M11368"))
	tmp		<- df[, which( CD4A>10000 & Patient=="M11368")]
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-3,d=0))
	cat(paste("\npatient M32958"))
	tmp		<- df[, which( Patient=="M32958" & PosCD4<="2010-05-03")]
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	cat(paste("\npatient M20633"))
	tmp		<- df[, which( Patient=="M20633" & CD4A>5000)]
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-3,d=0))
	tmp		<- df[, which( CD4A>3000)]
	cat(paste("\nnumber of entries with likely wrong CD4 units > 3000, n=",length(tmp),"DIVIDE BY 10"))
	#print(df[tmp,])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	#
	#	remove too small entries and digit entries
	#	
	df		<- subset(df, CD4A>10)		
	df		<- subset(df, (CD4A%%1)==0 )
	#
	#	remove duplicate entries 
	#
	setkey(df, Patient, PosCD4)
	df		<- df[, list(CD4A=round(mean(CD4A)), Source=Source[1] ), by=c('Patient','PosCD4')]	
	#
	#	remove consecutive jumps >1e3 
	#
	setkey(df, Patient, PosCD4)
	tmp		<- subset( df[, {  tmp<- which( length(CD4A)>1 & diff(CD4A)>1000 ); list(PosCD4=ifelse(length(tmp), as.character(PosCD4[tmp+1]), NA_character_), CD4A=ifelse(length(tmp), CD4A[tmp+1], NA_real_)) }, by='Patient'], !is.na(PosCD4))
	cat(paste("\nnumber of consecutive jumps >1e3 DELETE n=",nrow(tmp)))
	tmp		<- sapply(seq_len(nrow(tmp)), function(i)		df[,which(Patient==tmp[i,Patient] & PosCD4==tmp[i,PosCD4] & CD4A==tmp[i,CD4A])])
	df		<- df[-tmp,]
	#
	#	remove consecutive jumps >7e2 within 60days 
	#
	tmp		<- subset( df[, {  tmp<- which( length(CD4A)>1 & diff(CD4A)>700 & difftime(PosCD4[-1],PosCD4[-length(PosCD4)],units='days')<60 ); list(PosCD4=ifelse(length(tmp), as.character(PosCD4[tmp+1]), NA_character_), CD4A=ifelse(length(tmp), CD4A[tmp+1], NA_real_)) }, by='Patient'], !is.na(PosCD4))
	cat(paste("\nnumber of consecutive jumps >7e2 within 60days DELETE n=",nrow(tmp)))
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
	cat(paste("\nstill high entries NOT DELETED n=",nrow(tmp)))
	#	divide by 10
	tmp		<- which(df[, Patient%in%c("M10212","M14927","M15431","M15519","M20720","M26334","M27643","M27571") & CD4A>1700])
	cat(paste("\nnumber of entries with likely wrong CD4 units > 1700  M10212 M27571 M14927  M15431  M15519  M20720  M26334  M27643, n=",length(tmp),"DIVIDE BY 10"))
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	cat(paste("\npatient M17554"))
	tmp		<- which(df[, Patient=="M17554" & CD4A%in%c(2760,2170)])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	cat(paste("\npatient M17554"))
	tmp		<- which(df[, Patient=="1500" & CD4A%in%c(1500)])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	#	set to NA
	tmp		<- which(df[, Patient%in%c("M12953","M13340","M26537","M26669","M35668") & CD4A>1900])
	cat(paste("\nnumber of entries with likely wrong CD4 units > 1900  M12953   M13340  M26537  M26669  M35668, n=",length(tmp),"SET to NA"))
	set(df, tmp, "CD4A", NA)
	tmp		<- which(df[, Patient%in%c("M15743") & CD4A>2500])
	cat(paste("\nnumber of entries with likely wrong CD4 units > 2500  M15743, n=",length(tmp),"SET to NA"))
	set(df, tmp, "CD4A", NA)
	tmp		<- which(df[, Patient%in%c("M30155") & CD4A>1000])
	cat(paste("\nnumber of entries with likely wrong CD4 units > 1000  M30155, n=",length(tmp),"SET to NA"))
	set(df, tmp, "CD4A", NA)
	#
	df		<- subset(df, !is.na(CD4A))
	cat(paste("\nnumber of entries with !is.na(PosCD4), n=",nrow(df)))
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
CD4.161027<- function()
{
	#	TODO check with Ard that max source should be kept?
	verbose			<- 1
	
	#	INPUT FILE NAMES	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/original_latest/SHM_1602_161102_OR_ALL_Immu.csv'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_CD4.rda'
	#
	#	read CD4 csv data file and preprocess
	#
	df				<- read.csv(file, stringsAsFactors=FALSE)		
	df				<- hivc.db.reset.Date(df, col="CD4_D", NA.time=c("","01/01/1911","11/11/1911","24/06/1923"), date.format="%d/%m/%Y")
	df				<- data.table(df, key="Patient")
	setnames(df, "CD4_D","PosCD4")
	#
	#	data corrections from Ard
	#
	tmp		<- df[, which( Patient=="M11392" & PosCD4=="2000-10-30" & CD4A==3001)]
	if(verbose) cat(paste("\nnumber of entries with corrected CD4 units Patient=='M11392' & CD4A>3001, n=",length(tmp),"double entry, remove"))
	set(df, tmp, "CD4A", NA)	
	tmp		<- df[, which( Patient=="M26334" & PosCD4=="2008-09-11" & CD4A==2230)]
	if(verbose) cat(paste("\nnumber of entries with corrected CD4 units Patient==M26334 & PosCD4==2008-09-11 & CD4A==2230, n=",length(tmp),"double entry, remove"))
	set(df, tmp, "CD4A", NA)		
	tmp		<- df[, which( 	Patient=='M26293' & PosCD4=='2003-09-12' & CD4A==800)]
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M26293' & PosCD4=='2003-09-12' & CD4A==800, n=",length(tmp),"set to 16/9/2003"))
	set(df, tmp, "PosCD4", as.Date('2003-09-16'))
	tmp		<- df[, which( 	Patient=='M11233' & PosCD4=='2001-01-10' & CD4A==24)]
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M11233' & PosCD4=='2001-01-10' & CD4A==24, n=",length(tmp),"set to 25/4/1996"))
	set(df, tmp, "PosCD4", as.Date('1996-04-25'))
	tmp		<- df[, which( 	Patient=='M13124' & PosCD4=='1998-08-07')]
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M13124' & PosCD4=='1998-08-07', n=",length(tmp),"set to 7/9/1998"))
	set(df, tmp, "PosCD4", as.Date('1998-09-07'))	
	tmp		<- df[, which( 	Patient=='M13124' & PosCD4=='2001-09-14')]
	if(verbose) cat(paste("\nnumber of entries with corrected PosCD4  Patient=='M13124' & PosCD4=='2001-09-14', n=",length(tmp),"set to 14/9/2000"))
	set(df, tmp, "PosCD4", as.Date('2000-09-14'))
	tmp		<- df[, which(	Patient=="M10544" & PosCD4=="2003-02-17" & CD4A==850	|
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
			)]
	cat(paste("\nnumber of entries with incorrect CD4  should be 45, n=",length(tmp),"set to NA"))	
	set(df, tmp, "CD4A", NA)
	df		<- subset(df, !is.na(CD4A))
	cat(paste("\nnumber of entries read, n=",nrow(df)))	
	#
	#	if duplicate per day, keep those with higher Source variable
	#	
	df[, DUMMY:= seq_len(nrow(df))]
	tmp		<- df[, list(N=list(length(CD4A))), by=c('Patient','PosCD4')]	
	tmp		<- merge(subset(tmp, N>1), df, by=c('Patient','PosCD4'))
	tmp		<- tmp[, list(IDX=DUMMY[which(Source!=max(Source))]), by=c('Patient','PosCD4')]
	cat('\nremoving day duplicates with lower Source, n=', nrow(tmp))
	set(df, tmp[,IDX], 'PosCD4', NA_real_)
	set(df, NULL, 'DUMMY', NULL)
	df		<- subset(df, !is.na(PosCD4))	
	#
	#	remove identical entries
	#
	setkey(df, Patient, PosCD4, CD4A)
	df		<- unique(df)
	#
	#	adjust likely wrong units > 20000
	#
	tmp		<- df[, which( CD4A>20000)]
	cat(paste("\nnumber of entries with likely wrong CD4 units > 20000, n=",length(tmp),"DIVIDE BY 1e3"))
	#print(df[tmp,])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-3,d=0))
	#	adjust likely wrong units > 10000
	cat(paste("\npatient M11368"))
	tmp		<- df[, which( CD4A>10000 & Patient=="M11368")]
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-3,d=0))
	cat(paste("\npatient M32958"))
	tmp		<- df[, which( Patient=="M32958" & PosCD4<="2010-05-03")]
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	cat(paste("\npatient M20633"))
	tmp		<- df[, which( Patient=="M20633" & CD4A>5000)]
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-3,d=0))
	tmp		<- df[, which( CD4A>3000)]
	cat(paste("\nnumber of entries with likely wrong CD4 units > 3000, n=",length(tmp),"DIVIDE BY 10"))
	#print(df[tmp,])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	#
	#	remove too small entries and digit entries
	#	
	df		<- subset(df, CD4A>10)		
	df		<- subset(df, (CD4A%%1)==0 )
	#
	#	remove duplicate entries 
	#
	setkey(df, Patient, PosCD4)
	df		<- df[, list(CD4A=round(mean(CD4A)), Source=Source[1] ), by=c('Patient','PosCD4')]	
	#
	#	remove consecutive jumps >1e3 
	#
	setkey(df, Patient, PosCD4)
	tmp		<- subset( df[, {  tmp<- which( length(CD4A)>1 & diff(CD4A)>1000 ); list(PosCD4=ifelse(length(tmp), as.character(PosCD4[tmp+1]), NA_character_), CD4A=ifelse(length(tmp), CD4A[tmp+1], NA_real_)) }, by='Patient'], !is.na(PosCD4))
	cat(paste("\nnumber of consecutive jumps >1e3 DELETE n=",nrow(tmp)))
	tmp		<- sapply(seq_len(nrow(tmp)), function(i)		df[,which(Patient==tmp[i,Patient] & PosCD4==tmp[i,PosCD4] & CD4A==tmp[i,CD4A])])
	df		<- df[-tmp,]
	#
	#	remove consecutive jumps >7e2 within 60days 
	#
	tmp		<- subset( df[, {  tmp<- which( length(CD4A)>1 & diff(CD4A)>700 & difftime(PosCD4[-1],PosCD4[-length(PosCD4)],units='days')<60 ); list(PosCD4=ifelse(length(tmp), as.character(PosCD4[tmp+1]), NA_character_), CD4A=ifelse(length(tmp), CD4A[tmp+1], NA_real_)) }, by='Patient'], !is.na(PosCD4))
	cat(paste("\nnumber of consecutive jumps >7e2 within 60days DELETE n=",nrow(tmp)))
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
	cat(paste("\nstill high entries NOT DELETED n=",nrow(tmp)))
	#	divide by 10
	tmp		<- which(df[, Patient%in%c("M10212","M14927","M15431","M15519","M20720","M26334","M27643","M27571") & CD4A>1700])
	cat(paste("\nnumber of entries with likely wrong CD4 units > 1700  M10212 M27571 M14927  M15431  M15519  M20720  M26334  M27643, n=",length(tmp),"DIVIDE BY 10"))
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	cat(paste("\npatient M17554"))
	tmp		<- which(df[, Patient=="M17554" & CD4A%in%c(2760,2170)])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	cat(paste("\npatient M17554"))
	tmp		<- which(df[, Patient=="1500" & CD4A%in%c(1500)])
	set(df, tmp, "CD4A", round(df[tmp, CD4A]*1e-1,d=0))
	#	set to NA
	tmp		<- which(df[, Patient%in%c("M12953","M13340","M26537","M26669","M35668") & CD4A>1900])
	cat(paste("\nnumber of entries with likely wrong CD4 units > 1900  M12953   M13340  M26537  M26669  M35668, n=",length(tmp),"SET to NA"))
	set(df, tmp, "CD4A", NA)
	tmp		<- which(df[, Patient%in%c("M15743") & CD4A>2500])
	cat(paste("\nnumber of entries with likely wrong CD4 units > 2500  M15743, n=",length(tmp),"SET to NA"))
	set(df, tmp, "CD4A", NA)
	tmp		<- which(df[, Patient%in%c("M30155") & CD4A>1000])
	cat(paste("\nnumber of entries with likely wrong CD4 units > 1000  M30155, n=",length(tmp),"SET to NA"))
	set(df, tmp, "CD4A", NA)
	#
	df		<- subset(df, !is.na(CD4A))
	cat(paste("\nnumber of entries with !is.na(PosCD4), n=",nrow(df)))
	#
	df		<- df[, 	{
				z<- which.min(PosCD4)
				list(PosCD4=PosCD4, CD4=CD4A, PosCD4_T1=PosCD4[z], CD4_T1=CD4A[z] ) 	
			},by=Patient]		
	if(verbose) cat(paste("\nsave to", outfile))
	save(df, file=outfile)		
}
######################################################################################
CombinedTable.130301<- function(dir.name= DATA, verbose=1, resume=0)
{	
	require(data.table)		
	
	#input files generated with "project.hivc.Excel2dataframe"
	resume			<- 0
	verbose			<- 1
	if(1)
	{
		file.seq		<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.R",sep='/')
		file.primoSHM	<- paste(dir.name,"original/ATHENA_1501_primoSHM.csv",sep='/')
		file.patient	<- paste(dir.name,"derived/ATHENA_2013_03_Patient.R",sep='/')
		file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
		file.immu		<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
		file.treatment	<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.R",sep='/')
		file.out		<- paste(dir.name,"derived/ATHENA_2013_03_AllSeqPatientCovariates.R",sep='/')
	}
	if(0)
	{
		file.seq		<- paste(dir.name,"derived/ATHENA_2013_03_Sequences_AllMSM.R",sep='/')
		file.primoSHM	<- paste(dir.name,"original/ATHENA_1501_primoSHM.csv",sep='/')
		file.patient	<- paste(dir.name,"derived/ATHENA_2013_03_Patient_AllMSM.R",sep='/')
		file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro_AllMSM.R",sep='/')
		file.immu		<- paste(dir.name,"derived/ATHENA_2013_03_Immu_AllMSM.R",sep='/')
		file.treatment	<- paste(dir.name,"derived/ATHENA_2013_03_Regimens_AllMSM.R",sep='/')
		file.out		<- paste(dir.name,"derived/ATHENA_2013_03_AllSeqPatientCovariates_AllMSM.R",sep='/')
	}
	if(1)
	{
		dir.name		<- "~/Dropbox (SPH Imperial College)/2015_ATHENA_May_Update"
		file.seq		<- paste(dir.name,"ATHENA_1502_All_Sequences.R",sep='/')
		file.patient	<- paste(dir.name,"ATHENA_1502_All_Patient.R",sep='/')
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
		df		<- do.call('rbind',lapply( seq_along(df),	function(i)
						{
							z	<- data.table(df[[i]][,	c("Patient","DateRes","Subtype","SampleCode")], key="SampleCode")
							z[, GENE:=names(df)[i]]
							z
						}))
		df		<- dcast.data.table(df, Patient+DateRes+Subtype+SampleCode~GENE, value.var='GENE')
		set(df, NULL, "SampleCode", gsub(' ','', df[, SampleCode]))
		setnames(df, c("SampleCode","DateRes"), c("FASTASampleCode","PosSeqT"))
		#set(df, NULL, "FASTASampleCode", factor(df[,FASTASampleCode]))
		#set(df, NULL, "Patient", factor(df[,Patient]))		
		df.all	<- subset(df, select=c(FASTASampleCode,Patient,Subtype,PosSeqT))
		if(verbose)		cat(paste("\nnumber of sequences found, n=", nrow(df.all)))
		#
		# add Patient data
		#
		if(verbose)		cat(paste("\nadding patient data"))
		load(file.patient)
		df[, Subtype:=NULL]
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
		set(tmp, NULL, 'Patient', tmp[, as.character(Patient)])
		df.all	<- merge(df.all, unique(tmp), by='Patient', all.x=1)
		load(file.treatment)
		tmp		<- subset(df, select=c(Patient, AnyT_T1))		
		set(tmp, NULL, 'Patient', tmp[, as.character(Patient)])
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
		set(df, NULL, 'Patient', df[, as.character(Patient)])
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
		set(df, NULL, 'Patient', df[, as.character(Patient)])
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
		set(df, NULL, 'Patient', df[, as.character(Patient)])
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
		
		#	reset Acute field		
		df.all[, DUMMY:= df.all[, as.numeric(difftime(AnyPos_T1, NegT, units="days"))/365]]		
		#tmp	<- df.all[, which(!is.na(DUMMY) & DUMMY<=1 & (Acute_Spec=='CLIN' | Acute_Spec=='SYM' | is.na(Acute_Spec)))]
		tmp	<- df.all[, which(!is.na(DUMMY) & DUMMY<=1)]
		cat(paste('\nFound patients with NegT<1yr, n=', length(tmp)))
		set( df.all, tmp, 'Acute_Spec', 'NEGTEST' )
		# df.all[, table(!is.na(DUMMY) & DUMMY<=1, Acute_Spec, useNA='if')]
		df.all[, isAcuteNew:= NA_character_]
		set(df.all, df.all[, which(Acute_Spec%in%c('CLIN','LAB','NEGTEST'))], 'isAcuteNew', 'Yes')		
		set(df.all, df.all[, which(isAcute=='No' & is.na(Acute_Spec))], 'isAcuteNew', 'No')			
		#df.all[, table(isAcuteNew, Acute_Spec, useNA='if')]
		#df.all[, table(isAcuteNew, isAcute, useNA='if')]
		set(df.all, NULL, 'isAcuteNew', df.all[, factor(isAcuteNew)])
		#
		df.all[, DUMMY:=NULL]
		df.all[, Acute_Spec_1:=NULL]
		df.all[, Acute_Spec_2:=NULL]
		df.all[, Acute_Spec_3:=NULL]
		df.all[, Acute_Spec_4:=NULL]
		#
		#	add participation in primo-SHM
		#
		df.pr	<- as.data.table(read.csv(file.primoSHM, stringsAsFactors=FALSE, header=FALSE, col.names='Patient'))
		df.pr[, WithPrimoSHMID:='Yes']
		df.all	<- merge(df.all, df.pr, by='Patient', all.x=1)
		set(df.all, df.all[, which(is.na(WithPrimoSHMID))], 'WithPrimoSHMID', 'No' )
		set(df.all, NULL, 'WithPrimoSHMID', df.all[, factor(WithPrimoSHMID, levels=c('No','Yes'), labels=c('No','Yes'))])
		#setkey(df.all, Patient); unique(df.all)[, table(WithPrimoSHMID)]		
		if(verbose)	cat(paste("\nsave to file",file.out))
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
Regimen.160227<- function()
{
	#	BASIC VARIABLE DEFINITIONS		
	TR.failure 		<- c(21, 31 ) 								#viro failure
	TR.immufailure	<- c(32, 35)								#either immu failure or new CDC-B/C event
	TR.toxicity  	<- c(24, 34) 								#toxicity
	TR.adherence	<- c(47)
	TR.patrel		<- c(23, 33, 42, 43)						#either patient s decision, desired pregnancy or pregnancy
	date.var		<- c("StartTime","StopTime")
	date.format		<- '%Y-%m-%d'
	TR.notyet		<- as.Date("2016-01-01")	#set stop time for ongoing ART episodes
	verbose			<- 1
	#
	#	INPUT FILE
	#
	dir.name		<- '~/Dropbox (SPH Imperial College)/2015_ATHENA_May_Update'
	file			<- file.path(dir.name, 'ATHENA_1502_All_Regimens.csv')
	#file			<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.csv",sep='/')
	#file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	#file			<- paste(dir.name,"derived/ATHENA_2013_03_Regimens_AllMSM.csv",sep='/')	
	#file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro_AllMSM.R",sep='/')	
	#load(file.viro)
	#df.viro		<- subset(df, select=c(Patient, PosRNA, lRNA))
	#
	#	read REGIMEN csv data file and preprocess
	#
	df		<- as.data.table(read.csv(file, stringsAsFactors=FALSE))	
	set(df, NULL,'StartTime', df[, as.Date(StartTime, format=date.format)])	
	set(df, NULL,'StopTime', df[, as.Date(StopTime, format=date.format)])
	#	remove patients with empty fields
	tmp		<- df[, which(is.na(StopTime) & is.na(StartTime))]
	cat('\nRemove patients with no StartTime and no StopTime, n=', length(tmp))
	df		<- subset(df, !is.na(StopTime) | !is.na(StartTime))
	#	set stop times for ongoing ART episodes
	tmp		<- df[, which(is.na(StopTime))]
	cat('\nPatients with no StopTime: ART ongoing. Set to StopTime ',as.character(TR.notyet),', n=', length(tmp))
	set(df, tmp, 'StopTime', TR.notyet)
	cat(paste("\nnumber of entries, n=",nrow(df)))
	#
	#	Manual curations discussed with Ard (I think these are fixed)
	# 
	tmp						<- which(df[, Patient=="M18204" & StartTime=="1997-02-15"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M18204  1997-02-15"))
		set(df, tmp, "StopTime", as.Date("1999-03-15"))				
	}
	tmp						<- which(df[, Patient=="M29191" & StartTime=="2009-07-01"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M29191  2009-07-01"))
		set(df, tmp, "StopTime", as.Date("2010-06-04"))				
	}
	tmp						<- which(df[, Patient=="M33034" & StartTime=="2006-11-18"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M33034  2008-07-01"))
		set(df, tmp, "StopTime", as.Date("2008-06-15"))				
	}
	tmp						<- which(df[, Patient=="M35936" & StartTime=="2011-02-15"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M35936  2011-07-01"))
		set(df, tmp, "StopTime", as.Date("2011-04-29"))				
	}
	tmp						<- which(df[, Patient=="M37778" & StartTime=="2012-05-02"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M37778  2012-07-01"))
		set(df, tmp, "StopTime", as.Date("2012-06-20"))				
	}
	tmp						<- which(df[, Patient=="M41695" & StartTime=="2012-06-28"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M41695  2013-07-01"))
		set(df, tmp, "StopTime", as.Date("2013-05-21"))				
	}
	tmp						<- which(df[, Patient=="M29531" & StartTime=="1997-07-01"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M29531  1997-07-01"))
		set(df, tmp, "StartTime", as.Date("1997-05-30"))				
	}
	tmp						<- which(df[, Patient=="M33706" & StartTime=="2008-12-23"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M33706  2010-07-01"))
		set(df, tmp, "StopTime", as.Date("2010-11-15"))				
	}
	tmp						<- which(df[, Patient=="M33706" & StartTime=="2010-07-01" & StopTime=="2011-07-01"])
	if(length(tmp))
	{
		cat(paste("\nrm double entry 	M33706  2010-07-01"))
		set(df, tmp, c("StartTime",'StopTime'), NA_real_)				
	}
	tmp						<- which(df[, Patient=="M33706" & StartTime=="2010-11-15" & StopTime=="2011-07-01"])
	if(length(tmp))
	{
		cat(paste("\nrm double entry part 1 	M33706  2010-11-15"))
		set(df, tmp, "StopTime", as.Date("2011-01-17"))								
	}
	tmp						<- which(df[, Patient=="M33706" & StartTime=="2011-07-01" & StopTime=="2011-01-17"])
	if(length(tmp))
	{
		cat(paste("\nrm double entry part 2 	M33706  2011-07-01"))	
		set(df, tmp, c("StartTime",'StopTime'), NA_real_)				
	}	
	tmp						<- which(df[, Patient=="M15834" & StartTime=="1998-06-26"])
	if(length(tmp))
	{
		cat(paste("\nrm double entry 	M15834  1998-06-26"))
		set(df, tmp, c("StartTime",'StopTime'), NA_real_)				
	}
	tmp						<- which(df[, Patient=="M17493" & StartTime=="1996-01-01"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M17493  1996-01-01 manually to somewhere in 1996"))
		set(df, tmp, "StartTime", as.Date("1996-07-01"))				
	}	
	tmp						<- which(df[, Patient=="M41688" & StartTime=="2011-01-01"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M41688  2011-01-01 manually to 2011-11-01 (Ard)"))
		set(df, tmp, "StartTime", as.Date("2011-11-01"))
	}				
	tmp						<- which(df[, Patient=="M42092" & StartTime=="2011-08-23"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M42092  2011-08-23 to 2012-08-23"))
		set(df, tmp, "StartTime", as.Date("2012-08-23"))
	}	
	tmp						<- which(df[, Patient=="M42186" & StartTime=="2010-10-23"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M42186  2010-10-23 to 2012-10-23"))
		set(df, tmp, "StartTime", as.Date("2012-10-23"))
	}				
	tmp						<- which(df[, Patient=="M37531" & StartTime=="2006-08-15"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M37531  2006-08-15 to 2008-08-15"))
		set(df, tmp, "StartTime", as.Date("2008-08-15"))
	}
	tmp						<- which(df[, Patient=="M15834" & StartTime=="1998-06-26"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M15834  1911-11-11 to 1998-05-18"))
		set(df, tmp, "StopTime", as.Date("1997-12-08")+difftime( as.Date("1998-05-18"),as.Date("1997-12-08"), units='days' ))
	}
	tmp						<- which(df[, Patient=="M11814" & StartTime=="1996-12-15"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M11814  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1996-12-15")+difftime( as.Date("1998-10-30"),as.Date("1996-12-15"), units='days' ) )
	}
	df		<- subset(df, !is.na(StopTime) | !is.na(StartTime))
	#
	#	set potentially inaccurate StartTime
	#
	df[, StartTime_Acc:='Acc'] 
	set(df, df[, which(is.na(StartTime))], 'StartTime_Acc', NA_character_)
	set(df, df[, which(!is.na(StartTime) & as.POSIXlt(StartTime)$mday==15)], 'StartTime_Acc', "NAccD")
	set(df, df[, which(!is.na(StartTime) & as.POSIXlt(StartTime)$mon==6 & as.POSIXlt(StartTime)$mday==1)], 'StartTime_Acc', "NAccMD")
	set(df, df[, which(StartTime%in%as.Date(c("1911-01-01","1911-11-01","1911-11-11")))], 'StartTime_Acc', "NAccYMD")
	cat('\nNumber inaccurate StartTimes:')
	df[, table(StartTime_Acc)]	
	#
	#	set potentially inaccurate StopTime
	#
	df[, StopTime_Acc:='Acc'] 
	set(df, df[, which(is.na(StopTime))], 'StopTime_Acc', NA_character_)
	set(df, df[, which(!is.na(StopTime) & as.POSIXlt(StopTime)$mday==15)], 'StopTime_Acc', "NAccD")
	set(df, df[, which(!is.na(StopTime) & as.POSIXlt(StopTime)$mon==6 & as.POSIXlt(StopTime)$mday==1)], 'StopTime_Acc', "NAccMD")
	set(df, df[, which(StopTime%in%as.Date(c("1911-01-01","1911-11-01","1911-11-11")))], 'StopTime_Acc', "NAccYMD")
	cat('\nNumber inaccurate StopTimes:')
	df[, table(StopTime_Acc)]
	#
	#	handle Ard's bogus values	#0: manual curation
	#		
	tmp		<- which(df[, Patient=="M26087" & StartTime=="2005-05-09" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M26087  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2008-12-09"))
	}
	tmp		<- which(df[, Patient=="M17324" & StartTime=="1993-02-10" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M17324  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1995-03-15"))
	}
	tmp		<- which(df[, Patient=="M17324" & StartTime=="1995-03-15" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M17324  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1996-04-17"))
	}
	tmp		<- which(df[, Patient=="M17324" & StartTime=="1996-04-17" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M17324  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1996-10-14"))
	}
	tmp		<- which(df[, Patient=="M25962" & StartTime=="1996-07-01" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M25962  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1999-04-15"))
	}
	tmp		<- which(df[, Patient=="M25962" & StartTime=="1999-04-15" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M25962  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2003-03-20"))
	}
	tmp		<- which(df[, Patient=="M26044" & StartTime=="2008-04-01" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M26044  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2008-09-25"))
	}
	tmp		<- which(df[, Patient=="M26044" & StartTime=="2008-09-25" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M26044  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2009-05-07"))
	}
	tmp		<- which(df[, Patient=="M31241" & StartTime=="1989-09-15" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M31241  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1994-01-15"))
	}
	tmp		<- which(df[, Patient=="M31241" & StartTime=="1994-01-15" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M31241  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1996-02-16"))
	}
	tmp		<- which(df[, Patient=="M33356" & StartTime=="1991-07-01" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M33356  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1994-07-01"))
	}
	tmp		<- which(df[, Patient=="M33356" & StartTime=="1994-07-01" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M33356  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2005-07-01"))
	}
	tmp		<- which(df[, Patient=="M39099" & StartTime=="2002-07-01" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M39099  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2004-12-15"))
	}
	tmp		<- which(df[, Patient=="M39099" & StartTime=="2004-12-15" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M39099  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2010-09-08"))
	}	
	#
	#	handle Ard's bogus values	#1
	#
	setkey(df, Patient, StopTime)
	tmp		<- c( 	df[, which(is.na(StartTime) & StopTime%in%as.Date(c("1911-01-01","1911-07-01","1911-11-01","1911-11-11")))],
			df[, which(StartTime%in%as.Date(c("1911-01-01","1911-07-01","1911-11-01","1911-11-11")) & StopTime%in%as.Date(c("1911-01-01","1911-11-01","1911-11-11")))]	)
	if(length(tmp))
	{
		cat('\nremove bogus entries with NA StartTime and "1911-01-01","1911-11-01","1911-07-01","1911-11-11" StopTime, n=', length(tmp))
		set(df, tmp, 'StopTime', NA_real_)
	}
	df		<- subset(df, !is.na(StopTime))
	#
	#	set all bogus values to 1911-11-11 for simplicity
	#
	set(df, df[,which(StartTime%in%as.Date(c("1911-01-01","1911-07-01","1911-11-01")))], 'StartTime',as.Date("1911-11-11"))
	set(df, df[,which(StopTime%in%as.Date(c("1911-01-01","1911-07-01","1911-11-01")))], 'StopTime',as.Date("1911-11-11")) 
	#
	#	handle Ard's bogus values	#2: set bogus values for last episodes to StartTime+1
	#
	df[, DUMMY:=seq_len(nrow(df))]
	tmp		<- df[, list(N=length(which(StopTime=="1911-11-11"))), by='Patient']	
	tmp		<- merge(subset(tmp, N>=1), df, by='Patient')
	tmp2	<- tmp[, {
				z		<- as.character(StartTime[ which(StopTime=="1911-11-11") ])
				z		<- z[ z==as.character(max( max(StartTime), max(StopTime) )) ]
				if(length(z)==0)
					z	<- NA_character_
				list(StartTime=z)
			}, by='Patient']
	tmp2	<- subset(tmp2, !is.na(StartTime))
	for(i in seq_len(nrow(tmp2)))
	{
		cat('\nSet 1911-11-11 of',tmp2[i,Patient],'to',as.character(tmp2[i,as.Date(StartTime)+1]))
		z	<- df[, which(Patient==tmp2[i,Patient] & StartTime==tmp2[i,StartTime])]
		stopifnot(length(z)==1)
		set(df, z, 'StopTime', tmp2[i,as.Date(StartTime)+1])		
	}
	#
	#	handle Ard's bogus values	#3: remove 1911-11-11 pairs by midpoint 
	#
	tmp		<- df[, list(Start_N=length(which(StartTime=="1911-11-11")), Stop_N=length(which(StopTime=="1911-11-11"))), by='Patient']		
	tmp		<- merge(subset(tmp, Start_N==1 & Stop_N%in%c(1,2)), df, by='Patient')
	tmp		<- tmp[, {
				z		<- StartTime[ which(StopTime=="1911-11-11") ]
				zz		<- which(StartTime=="1911-11-11")
				if(length(zz)>0)
				{
					z	<- z[ which.min(difftime(StopTime[zz], z, units='days')) ]					
				}
				if(length(zz)==0)
				{
					z	<- NA_character_
					zz	<- NA_integer_
				}
				list(StartTime=z, StopTime=StopTime[zz], IDX1=DUMMY[z==StartTime], IDX2=DUMMY[zz])
			}, by='Patient']
	for(i in seq_len(nrow(tmp)))
	{
		z	<- difftime(tmp[i, as.Date(StopTime)], tmp[i, as.Date(StartTime)], units='days')/2
		stopifnot(z>0)
		z	<- tmp[i, as.Date(StartTime)]+z
		cat('\nSet 1911-11-11 of',tmp[i,Patient],'to',as.character(z))
		set(df, tmp[i,IDX1], 'StopTime',z)
		set(df, tmp[i,IDX2], 'StartTime',z)				
	}
	#
	#	handle Ard's bogus values	#4: set 1911-11-11 to next StartTime (only after #2)
	#	
	tmp		<- df[, list(Start_N=length(which(StartTime=="1911-11-11")), Stop_N=length(which(StopTime=="1911-11-11"))), by='Patient']
	stopifnot( !nrow(subset(tmp, Stop_N>1)) )
	tmp		<- merge(subset(tmp, Start_N==0 & Stop_N==1), df, by='Patient')
	tmp		<- tmp[, {
				z	<- which(StopTime=="1911-11-11")
				zz	<- sort(StartTime)
				list(IDX1=DUMMY[z], StopTime=zz[ which(zz==StartTime[z])+1L ])		
			}, by='Patient']
	for(i in seq_len(nrow(tmp)))
	{		
		cat('\nSet 1911-11-11 of',tmp[i,Patient],'to',tmp[i,as.character(StopTime)])
		stopifnot( df[tmp[i,IDX1],StartTime]<=tmp[i,StopTime])
		set(df, tmp[i,IDX1], 'StopTime',tmp[i,StopTime])			
	}
	tmp		<- df[, list(Start_N=length(which(StartTime=="1911-11-11")), Stop_N=length(which(StopTime=="1911-11-11"))), by='Patient']
	stopifnot( !nrow(subset(tmp, Stop_N>0)) )
	df[, DUMMY:=NULL]
	#
	#	removing Patients not on treatment (this should be legacy)
	#
	tmp						<- which(df[, is.na(StartTime) & StopTime==TR.notyet & NoDrug==0])
	if(length(tmp))
	{
		cat("\nnumber of entries with is.na(StartTime) & StopTime==TR.notyet & NoDrug==0",length(tmp),"Set StopTime to NA")
		set(df, tmp, "StopTime", NA)		
	}
	tmp						<- which(df[, StartTime==TR.notyet & StopTime==TR.notyet])
	if(length(tmp))
	{		
		cat("\nnumber of entries with StartTime==TR.notyet & StopTime==TR.notyet",length(tmp),"Set StartTime and StopTime to NA")
		set(df, tmp, "StartTime", NA)
		set(df, tmp, "StopTime", NA)
	}
	tmp						<- which(df[, StartTime==TR.notyet])
	if(length(tmp))
	{		
		cat("\nnumber of entries with StartTime==TR.notyet & StopTime!=TR.notyet",length(tmp),"Misclassified StartTime. Set to NA")
		set(df, tmp, "StartTime", NA)
	}	
	df						<- subset(df, !is.na(StopTime))
	#
	#	check patients with all NoDrug==NA (legacy)
	#
	tmp			<- nrow(subset( df[, list(select=any(is.na(NoDrug))), by='Patient'], select))
	cat("\nnumber of patients with is.na(NoDrug), n=", tmp)
	stopifnot( !length(tmp)	)
	#
	#	remove patients for which all periods have NoDrug==0
	#
	tmp			<- df[,  list(select=!all(NoDrug==0)) , by='Patient']	
	cat("\nPatients with all(NoDrug==0). Remove. n=",nrow(subset(tmp, !select)))
	df			<- merge( df, subset(tmp, select, Patient), by='Patient')
	#
	#	for each patient find first ART period with NoDrug>0 and keep only these
	#
	setkey(df, Patient, StopTime)	
	tmp			<- subset( df[, list( select= which(NoDrug>0)[1]>1 ), by='Patient'], select)[, Patient]
	cat(paste("\nPatients with NoDrug==0 in first episode(s). Remove these episodes. n=",length(tmp)))		
	df[, DUMMY:= seq_len(nrow(df))]
	tmp			<- df[, {					#need some akward code since StopTimes are not yet unique
				z<- seq.int( which(NoDrug>0)[1], length(NoDrug))
				list( DUMMY= DUMMY[z] )
			}, by='Patient']
	df			<- df[tmp[, DUMMY],]
	df[, DUMMY:=NULL]
	#	check that no NoDrug==0 at start
	tmp			<- subset( df[,list(check= any( is.na(StartTime) & NoDrug==0) ), by='Patient'], check, Patient)
	stopifnot( !nrow(tmp)	)
	tmp			<- subset( df[,list(check= NoDrug[1]==0 ), by='Patient'], check, Patient)
	stopifnot( !nrow(tmp)	)	
	#
	#	check Ards bogus values #5: set dangling StartTime to StopTime-1
	#
	tmp		<- df[, which(StartTime=="1911-11-11")]
	for( i in tmp)
		set(df, i, 'StartTime', df[i,StopTime-1])	
	#
	#	fix StartTime > StopTime: check if StopTime is StartTime
	#	
	tmp		<- data.table(IDX=which(df[, StartTime>StopTime]))	
	cat(paste("\nnumber of entries with StartTime>StopTime, n=",nrow(tmp)))
	df[, DUMMY:=seq_len(nrow(df))]
	tmp[, Patient:= df[tmp[,IDX],Patient]]
	tmp[, StartTime:= df[tmp[,IDX],StartTime]]
	tmp[, StopTime:= df[tmp[,IDX],StopTime]]
	for(i in seq_len(nrow(tmp)))
		if(nrow(subset(df, Patient==tmp[i,Patient] & StartTime==tmp[i,StopTime] & StopTime==tmp[i,StartTime])))
		{
			cat('\nStopTime is StartTime for patient',tmp[i,Patient])
			set(df, tmp[i,IDX], c('StartTime','StopTime'), NA_real_)
		}
	df		<- subset(df, !is.na(StartTime) | !is.na(StopTime))
	df[, DUMMY:=seq_len(nrow(df))]
	#
	#	fix StartTime > StopTime: fix inaccurate StopTime
	#		
	tmp		<- data.table(IDX=which(df[, StartTime>StopTime]))	
	cat(paste("\nnumber of entries with StartTime>StopTime, n=",nrow(tmp)))	
	tmp[, Patient:= df[tmp[,IDX],Patient]]
	tmp[, StartTime_Acc:=df[tmp[,IDX],StartTime_Acc]]
	tmp[, StopTime_Acc:=df[tmp[,IDX],StopTime_Acc]]	
	stopifnot( !nrow(subset(tmp,  StartTime_Acc=='Acc' & StopTime_Acc=='Acc')))	
	tmp2	<- merge(subset(tmp,  StartTime_Acc!='NAccMD'), df, by='Patient')	
	tmp2	<- tmp2[, list(Start_N=length(which(StartTime==StopTime[IDX==DUMMY])), Stop_N=length(which(StopTime==StopTime[IDX==DUMMY]))), by=c('Patient','IDX')]
	stopifnot( !nrow(subset(tmp2, Start_N>1))	)
	tmp2	<- merge(subset(tmp2, Start_N>0), df, by='Patient')
	tmp2	<- tmp2[, {
				z		<- StopTime[ IDX==DUMMY ]				# that s the inaccurate query time to be changed, eg 2003-07-01
				zz		<- which(StartTime==z)
				stopifnot(length(zz)==1)						# below we only compare just one start time 
				z		<- StartTime[ which(StopTime==z) ]		# that s the start times of all stop times that equal the query time
				stopifnot(length(z)>0)	
				z		<- z[ z<StopTime[zz] ]					
				stopifnot(length(z)>0)
				z		<- as.character( z[ which.min(difftime(StopTime[zz], z, units='days')) ] )									
				list(StartTime=z, StopTime=StopTime[zz], IDX1=DUMMY[z==StartTime], IDX2=DUMMY[zz])
			}, by=c('Patient','IDX')]
	for(i in seq_len(nrow(tmp2)))
	{
		cat('\nSet',as.character( df[tmp2[i,IDX1],StopTime] ), as.character( df[tmp2[i,IDX2],StartTime] ),'of',tmp2[i,Patient])
		z	<- difftime(tmp2[i, as.Date(StopTime)], tmp2[i, as.Date(StartTime)], units='days')/2
		stopifnot(z>0)
		z	<- tmp2[i, as.Date(StartTime)]+z
		stopifnot( df[tmp2[i,IDX1],StopTime]==df[tmp2[i,IDX2],StartTime], df[tmp2[i,IDX1],StopTime_Acc!='Acc'], df[tmp2[i,IDX2],StartTime_Acc!='Acc'] )		
		set(df, tmp2[i,IDX1], 'StopTime',z)
		set(df, tmp2[i,IDX2], 'StartTime',z)				
	}
	#
	#	fix StartTime > StopTime: fix dangling inaccurate StopTimes of StartTime>StopTime 
	#			
	tmp		<- data.table(IDX=which(df[, StartTime>StopTime]))	
	cat(paste("\nnumber of entries with StartTime>StopTime, n=",nrow(tmp)))	
	tmp[, Patient:= df[tmp[,IDX],Patient]]
	tmp[, StartTime_Acc:=df[tmp[,IDX],StartTime_Acc]]
	tmp[, StopTime_Acc:=df[tmp[,IDX],StopTime_Acc]]	
	tmp2	<- merge(subset(tmp,  StartTime_Acc!='NAccMD'), df, by=c('Patient','StartTime_Acc','StopTime_Acc'))
	tmp2	<- tmp2[, list(Start_N=length(which(StartTime==StopTime[IDX==DUMMY])), Stop_N=length(which(StopTime==StopTime[IDX==DUMMY]))), by=c('Patient','IDX')]
	stopifnot( !nrow(subset(tmp2, Start_N>0))	)	
	tmp2	<- merge(tmp2, df, by='Patient')
	tmp2	<- tmp2[, {
				z		<- StopTime[ IDX==DUMMY ]
				z		<- which(StopTime==z)
				zz		<- sort(StartTime)
				if( which(zz==StartTime[z])<length(zz))
					ans	<- zz[ which(zz==StartTime[z])+1L ]		#set to next start time
				if( which(zz==StartTime[z])==length(zz))
					ans	<- zz[length(zz)]+1
				list(IDX1=DUMMY[z], StopTime=ans)		
			}, by='Patient']
	for(i in seq_len(nrow(tmp2)))
	{		
		cat('\nSet ',as.character( df[tmp2[i,IDX1],StopTime] ),'of',tmp2[i,Patient],'to',tmp2[i,as.character(StopTime)])
		stopifnot( df[tmp2[i,IDX1],StartTime]<=tmp2[i,StopTime])
		set(df, tmp2[i,IDX1], 'StopTime',tmp2[i,StopTime])			
	}
	#
	#	fix StartTime > StopTime: fix inaccurate StartTime
	#		
	tmp		<- data.table(IDX=which(df[, StartTime>StopTime]))	
	cat(paste("\nnumber of entries with StartTime>StopTime, n=",nrow(tmp)))	
	tmp[, Patient:= df[tmp[,IDX],Patient]]
	tmp[, StartTime_Acc:=df[tmp[,IDX],StartTime_Acc]]
	tmp[, StopTime_Acc:=df[tmp[,IDX],StopTime_Acc]]	
	stopifnot( !nrow(subset(tmp,  StartTime_Acc=='Acc' & StopTime_Acc=='Acc')))	
	tmp2	<- merge(subset(tmp,  StopTime_Acc!='NAccMD'), df, by='Patient')	
	tmp2	<- tmp2[, list(Start_N=length(which(StartTime==StartTime[IDX==DUMMY])), Stop_N=length(which(StopTime==StartTime[IDX==DUMMY]))), by=c('Patient','IDX')]
	stopifnot( !nrow(subset(tmp2, Start_N>1)), !nrow(subset(tmp2, Start_N==0))	)
	tmp2	<- merge(subset(tmp2, Stop_N>0), df, by='Patient')
	stopifnot( !nrow(subset(tmp2, Start_N>1)), !nrow(subset(tmp2, Start_N==0))	)
	tmp2	<- tmp2[, {				
				z		<- StartTime[ IDX==DUMMY ]
				zz		<- which(StartTime==z)
				stopifnot(length(zz)==1)	
				z		<- StartTime[ which(StopTime==z) ]
				stopifnot(length(z)>0)	
				z		<- z[ z<StopTime[zz] ]					
				stopifnot(length(z)>0)				
				z		<- as.character( z[ which.min(difftime(StopTime[zz], z, units='days')) ] )									
				list(StartTime=z, StopTime=StopTime[zz], IDX1=DUMMY[z==StartTime], IDX2=DUMMY[zz])
			}, by=c('Patient','IDX')]
	for(i in seq_len(nrow(tmp2)))
	{
		cat('\nSet',as.character( df[tmp2[i,IDX1],StopTime] ), as.character( df[tmp2[i,IDX2],StartTime] ),'of',tmp2[i,Patient])
		z	<- difftime(tmp2[i, as.Date(StopTime)], tmp2[i, as.Date(StartTime)], units='days')/2
		stopifnot(z>0)
		z	<- tmp2[i, as.Date(StartTime)]+z
		stopifnot( df[tmp2[i,IDX1],StopTime]==df[tmp2[i,IDX2],StartTime], df[tmp2[i,IDX1],StopTime_Acc!='Acc'], df[tmp2[i,IDX2],StartTime_Acc!='Acc'] )		
		set(df, tmp2[i,IDX1], 'StopTime',z)
		set(df, tmp2[i,IDX2], 'StartTime',z)				
	}
	#
	#	fix StartTime > StopTime: fix inaccurate StartTime
	#		
	tmp		<- data.table(IDX=which(df[, StartTime>StopTime]))
	stopifnot( !nrow(tmp))	
	df[, DUMMY:=NULL]
	#
	#	remove duplicate entries
	#
	setkey(df, Patient, StartTime, StopTime)
	tmp		<- which(duplicated(df))
	if(length(tmp))
	{
		cat("\nnumber of duplicate entries with key  Patient, StartTime, StopTime, n=",nrow(tmp))
		df		<- unique(df)		
	}
	#
	#	remove equal entries
	#	
	tmp			<- df[, which(StartTime==StopTime)]
	if(length(tmp))
	{
		cat("\nrm entries with StartTime==StopTime, n=",length(tmp))
		df		<- df[-tmp,]		
	}
	#	
	#	should have manually fixed double entries (these were just 6)
	#
	setkey(df, Patient, StopTime)
	tmp		<- data.table(IDX=which(duplicated(df)))	
	tmp[, Patient:= df[tmp[,IDX],Patient]]
	stopifnot( !nrow(tmp) )
	#
	#	see if we can merge consecutive NoDrug==0 periods
	#
	df[, DUMMY:= seq_len(nrow(df))]
	tmp		<- unique(subset(df, NoDrug==0, Patient))
	tmp		<- merge(tmp, df, by='Patient')	
	tmp		<- tmp[,	{
				z		<- which(NoDrug==0)
				ans		<- NA_integer_
				if(length(z)>1)
					ans	<- DUMMY[ z[ which(diff(z)==1)+1L ] ]
				list(IDX=ans)
			} ,by='Patient']
	tmp		<- subset(tmp, !is.na(IDX))	
	cat("\nnumber of no drug episodes with follow up no drug episodes, n=",nrow(tmp))
	for(i in seq_len(nrow(tmp)))
	{
		cat('\nmerge no drug episodes for patient', tmp[i, Patient])
		stopifnot(df[tmp[i,IDX]-1L, Patient]==df[tmp[i,IDX], Patient],  df[tmp[i,IDX], NoDrug]==0,  df[tmp[i,IDX]-1L, NoDrug]==0)
		set(df, tmp[i,IDX], 'StartTime', df[tmp[i,IDX]-1L, StartTime])
		set(df, tmp[i,IDX], 'StartTime_Acc', df[tmp[i,IDX]-1L, StartTime_Acc])
		set(df, tmp[i,IDX]-1L, c('StartTime','StopTime'), NA_real_)
	}
	df[, DUMMY:=NULL]
	df		<- subset(df, !is.na(StartTime) | !is.na(StopTime))
	#
	#	TR.interrupted
	#
	if(nrow(df[which(is.na(df[,NoDrug]) & StartTime!=TR.notyet),])) 
		stop("unexpected NA in NoDrug when on treatment")
	tmp								<- rep(0, nrow(df))
	tmp[ which( df[,NoDrug==0] ) ]	<- 1
	df[, TrI:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#
	#	process treatment change reasons
	#
	if(any(!is.na(df[, Reason9])))	stop("unexpected !NA after Reason7")
	df.TrCh.noreason	<-  df[, which(is.na(Reason1)&is.na(Reason2)&is.na(Reason3)&is.na(Reason4)&is.na(Reason5)&is.na(Reason6)&is.na(Reason7)&is.na(Reason8)&(NoDrug!=0)) ] 
	#	TR.addition
	tmp	<- na.omit(c( 	df[df.TrCh.noreason, AZT]<=df[df.TrCh.noreason+1, AZT] &
							df[df.TrCh.noreason, DDI]<=df[df.TrCh.noreason+1, DDI] &
							df[df.TrCh.noreason, DDC]<=df[df.TrCh.noreason+1, DDC] &
							df[df.TrCh.noreason, D4T]<=df[df.TrCh.noreason+1, D4T] &
							df[df.TrCh.noreason, DTC]<=df[df.TrCh.noreason+1, DTC] &
							df[df.TrCh.noreason, ABC]<=df[df.TrCh.noreason+1, ABC] &
							df[df.TrCh.noreason, TDF]<=df[df.TrCh.noreason+1, TDF] &
							df[df.TrCh.noreason, FTC]<=df[df.TrCh.noreason+1, FTC] &
							df[df.TrCh.noreason, SQV]<=df[df.TrCh.noreason+1, SQV] &
							df[df.TrCh.noreason, IDV]<=df[df.TrCh.noreason+1, IDV] & 
							df[df.TrCh.noreason, RTV]<=df[df.TrCh.noreason+1, RTV] &
							df[df.TrCh.noreason, NFV]<=df[df.TrCh.noreason+1, NFV] &
							df[df.TrCh.noreason, LPV]<=df[df.TrCh.noreason+1, LPV] &
							df[df.TrCh.noreason, FPV]<=df[df.TrCh.noreason+1, FPV] &
							df[df.TrCh.noreason, ATV]<=df[df.TrCh.noreason+1, ATV] &
							df[df.TrCh.noreason, TPV]<=df[df.TrCh.noreason+1, TPV] &
							df[df.TrCh.noreason, DRV]<=df[df.TrCh.noreason+1, DRV] &
							df[df.TrCh.noreason, NVP]<=df[df.TrCh.noreason+1, NVP] &
							df[df.TrCh.noreason, EFV]<=df[df.TrCh.noreason+1, EFV] &
							df[df.TrCh.noreason, ETR]<=df[df.TrCh.noreason+1, ETR] &
							df[df.TrCh.noreason, RPV]<=df[df.TrCh.noreason+1, RPV] &				
							df[df.TrCh.noreason, DOR]<=df[df.TrCh.noreason+1, DOR] &
							df[df.TrCh.noreason, ENF]<=df[df.TrCh.noreason+1, ENF] &
							df[df.TrCh.noreason, RAL]<=df[df.TrCh.noreason+1, RAL] &
							df[df.TrCh.noreason, MVC]<=df[df.TrCh.noreason+1, MVC] &
							df[df.TrCh.noreason, EVG]<=df[df.TrCh.noreason+1, EVG] &
							df[df.TrCh.noreason, DTG]<=df[df.TrCh.noreason+1, DTG] &
							df[df.TrCh.noreason, COBI]<=df[df.TrCh.noreason+1, COBI] 
			))	
	#tmp						<- df[df.TrCh.noreason, NoDrug+1]==df[df.TrCh.noreason+1, NoDrug]	
	df[, TrCh.addition:=0]
	set(df, df.TrCh.noreason[tmp], 'TrCh.addition', 1)
	df[, TrCh.addition:= factor(TrCh.addition,levels=c(0,1),labels=c("No","Yes"))]
	df.TrCh.noreason		<- df.TrCh.noreason[!tmp]
	#	TR.failure
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.failure | Reason2%in%TR.failure | Reason3%in%TR.failure | Reason4%in%TR.failure | Reason5%in%TR.failure | Reason6%in%TR.failure | Reason7%in%TR.failure | Reason8%in%TR.failure ]) ]<- 1
	df[, TrCh.failure:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#	TR.immufailure
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.immufailure | Reason2%in%TR.immufailure | Reason3%in%TR.immufailure | Reason4%in%TR.immufailure | Reason5%in%TR.immufailure | Reason6%in%TR.immufailure | Reason7%in%TR.immufailure | Reason8%in%TR.immufailure ]) ]<- 1
	df[, TrCh.immufailure:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]	
	#	TR.toxicity
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.toxicity | Reason2%in%TR.toxicity | Reason3%in%TR.toxicity | Reason4%in%TR.toxicity | Reason5%in%TR.toxicity | Reason6%in%TR.toxicity | Reason7%in%TR.toxicity | Reason8%in%TR.toxicity]) ]<- 1
	df[, TrCh.toxicity:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]		
	#	TR.adherence
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.adherence | Reason2%in%TR.adherence | Reason3%in%TR.adherence | Reason4%in%TR.adherence | Reason5%in%TR.adherence | Reason6%in%TR.adherence | Reason7%in%TR.adherence | Reason8%in%TR.adherence]) ]<- 1
	df[, TrCh.adherence:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#	TR.patient related
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.patrel | Reason2%in%TR.patrel | Reason3%in%TR.patrel | Reason4%in%TR.patrel | Reason5%in%TR.patrel | Reason6%in%TR.patrel | Reason7%in%TR.patrel | Reason8%in%TR.patrel]) ]<- 1
	df[, TrCh.patrel:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#
	#	set last stop time to start time - 1
	#
	df[, DUMMY:=seq_len(nrow(df))]
	tmp	<- df[, {
				ans		<- NA_integer_
				if(length(StartTime)>1)
					ans	<- DUMMY[ 1L+which(StartTime[-1]!=StopTime[-length(StopTime)]+1) ]
				list(IDX=ans)
			} , by='Patient']
	set(df, tmp[,IDX]-1L, 'StopTime', df[tmp[,IDX],StartTime]-1)	
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
						AnyT_T1		= subset(x, TrI=="No")[,StartTime][1], 
						AnyT_T1_Acc	= subset(x, TrI=="No")[,StartTime_Acc][1]		)																					
			}, by='Patient']			
	df		<- merge(subset(tmp,select=c(Patient,AnyT_T1,AnyT_T1_Acc,TrI.n, TrI.mo, TrI.p)), df, all.y=1, by="Patient")
	#
	#	reset AnyT_T1 conservatively
	#
	df[, AnyT_T1_Crude:=AnyT_T1]
	df	<- hivc.db.reset.ARTStartFromAccurate(df, col='AnyT_T1')			
	#	add largest lRNA at least 3 months after therapy start.
	if(0)
	{
		tmp					<- subset(df, select=c(Patient, AnyT_T1)) 
		setkey(tmp, Patient)
		tmp					<- merge(unique(tmp), df.viro, by='Patient')	
		tmp					<- merge(subset(df, select=c(Patient, AnyT_T1, StartTime, StopTime)), subset(tmp, difftime(PosRNA, AnyT_T1, units='days')>90,select=c(Patient, PosRNA, lRNA)), by='Patient', allow.cartesian=TRUE)
		tmp					<- subset(tmp, (StartTime<=PosRNA & PosRNA<=StopTime))
		tmp					<- tmp[, {
					z<- which.max(lRNA)
					list(AnyT_T1=AnyT_T1[z], PosRNA.mx= PosRNA[z], lRNA.mx=lRNA[z])
				}, by=c('Patient','StartTime','StopTime')]
		df					<- merge( df, subset(tmp, select=c(Patient, StartTime, StopTime, PosRNA.mx, lRNA.mx)), by=c('Patient','StartTime','StopTime'), all.x=TRUE)
		#	set failure by viral load during treatment period
		df[, TrVL.failure:=factor( df[, lRNA.mx>3], levels=c(FALSE, TRUE), labels=c('No','Yes'))]	
		#subset(df, TrVL.failure=='No' & TrCh.failure=='Yes')	viro failure can be indicated for VL < 1e3
		#subset(df, TrVL.failure=='No' & TrCh.failure=='Yes' & lRNA.mx<log10(200), select=c(Patient, StartTime, StopTime, AnyT_T1, PosRNA.mx, lRNA.mx))	#some are really suspicious!		
	}
	#
	df[, DUMMY:=NULL]
	file		<- paste(substr(file, 1, nchar(file)-3),'R',sep='')
	if(verbose) cat(paste("\nsave to", file))
	save(df, file=file)
}
######################################################################################
Regimen.161027<- function()
{
	require(gdata)
	#	BASIC VARIABLE DEFINITIONS		
	TR.failure 		<- c(21, 31 ) 								#viro failure
	TR.immufailure	<- c(32, 35)								#either immu failure or new CDC-B/C event
	TR.toxicity  	<- c(24, 34) 								#toxicity
	TR.adherence	<- c(47)
	TR.patrel		<- c(23, 33, 42, 43)						#either patient s decision, desired pregnancy or pregnancy
	date.var		<- c("StartTime","StopTime")
	date.format		<- '%Y-%m-%d'
	TR.notyet		<- as.Date("2016-06-01")	#set stop time for ongoing ART episodes
	verbose			<- 1
	#
	#	INPUT FILE
	#	
	file			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/original_latest/SHM_1602_161102_OR_ALL_Regimens.xlsx'
	outfile			<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_ART.rda'
	#file			<- paste(dir.name,"derived/ATHENA_2013_03_Regimens.csv",sep='/')
	#file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro.R",sep='/')
	#file			<- paste(dir.name,"derived/ATHENA_2013_03_Regimens_AllMSM.csv",sep='/')	
	#file.viro		<- paste(dir.name,"derived/ATHENA_2013_03_Viro_AllMSM.R",sep='/')	
	#load(file.viro)
	#df.viro		<- subset(df, select=c(Patient, PosRNA, lRNA))
	#
	#	read REGIMEN csv data file and preprocess
	#
	#dfo		<- copy(df)
	df		<- as.data.table(read.xls(file, stringsAsFactors=FALSE))
	#df		<- copy(dfo)
	set(df, NULL,'StartTime', df[, as.Date(StartTime, format=date.format)])	
	set(df, NULL,'StopTime', df[, as.Date(StopTime, format=date.format)])
	#	remove patients with empty fields
	tmp		<- df[, which(is.na(StopTime) & is.na(StartTime))]
	cat('\nRemove patients with no StartTime and no StopTime, n=', length(tmp))
	df		<- subset(df, !is.na(StopTime) | !is.na(StartTime))
	#	set stop times for ongoing ART episodes
	tmp		<- df[, which(is.na(StopTime))]
	cat('\nPatients with no StopTime: ART ongoing. Set to StopTime ',as.character(TR.notyet),', n=', length(tmp))
	set(df, tmp, 'StopTime', TR.notyet)
	cat(paste("\nnumber of entries, n=",nrow(df)))
	#
	#	Manual curations discussed with Ard (I think these are fixed)
	# 
	tmp		<- which(df[, Patient=="M46130" & StopTime=="1985-07-01"])
	set(df, tmp, "StartTime", NA_real_)	
	tmp						<- which(df[, Patient=="M35950" & StartTime=="2011-07-01"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M35950  2011-07-01"))
		set(df, tmp, "StopTime", TR.notyet)				
	}
	tmp						<- which(df[, Patient=="M29191" & StartTime=="2009-07-01"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M29191  2009-07-01"))
		set(df, tmp, "StopTime", as.Date("2010-06-04"))				
	}
	tmp						<- which(df[, Patient=="M33034" & StartTime=="2006-11-18"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M33034  2008-07-01"))
		set(df, tmp, "StopTime", as.Date("2008-06-15"))				
	}
	tmp						<- which(df[, Patient=="M37778" & StartTime=="2012-05-02"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M37778  2012-07-01"))
		set(df, tmp, "StopTime", as.Date("2012-06-20"))				
	}
	tmp						<- which(df[, Patient=="M41695" & StartTime=="2012-06-28"])
	if(length(tmp))
	{
		cat(paste("\nfix entry 	M41695  2013-07-01"))
		set(df, tmp, "StopTime", as.Date("2013-05-21"))				
	}
	tmp						<- which(df[, Patient=="M15834" & StartTime=="1998-06-26"])
	if(length(tmp))
	{
		cat(paste("\nrm double entry 	M15834  1998-06-26"))
		set(df, tmp, c("StartTime",'StopTime'), NA_real_)				
	}
	tmp						<- which(df[, Patient=="M15834" & StartTime=="1998-06-26"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M15834  1911-11-11 to 1998-05-18"))
		set(df, tmp, "StopTime", as.Date("1997-12-08")+difftime( as.Date("1998-05-18"),as.Date("1997-12-08"), units='days' ))
	}
	tmp						<- which(df[, Patient=="M11814" & StartTime=="1996-12-15"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M11814  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1996-12-15")+difftime( as.Date("1998-10-30"),as.Date("1996-12-15"), units='days' ) )
	}
	df		<- subset(df, !is.na(StopTime) | !is.na(StartTime))
	#
	#	set potentially inaccurate StartTime
	#
	df[, StartTime_Acc:='Acc'] 
	set(df, df[, which(is.na(StartTime))], 'StartTime_Acc', NA_character_)
	set(df, df[, which(!is.na(StartTime) & as.POSIXlt(StartTime)$mday==15)], 'StartTime_Acc', "NAccD")
	set(df, df[, which(!is.na(StartTime) & as.POSIXlt(StartTime)$mon==6 & as.POSIXlt(StartTime)$mday==1)], 'StartTime_Acc', "NAccMD")
	set(df, df[, which(StartTime%in%as.Date(c("1911-01-01","1911-11-01","1911-11-11")))], 'StartTime_Acc', "NAccYMD")
	cat('\nNumber inaccurate StartTimes:')
	df[, table(StartTime_Acc)]	
	#
	#	set potentially inaccurate StopTime
	#
	df[, StopTime_Acc:='Acc'] 
	set(df, df[, which(is.na(StopTime))], 'StopTime_Acc', NA_character_)
	set(df, df[, which(!is.na(StopTime) & as.POSIXlt(StopTime)$mday==15)], 'StopTime_Acc', "NAccD")
	set(df, df[, which(!is.na(StopTime) & as.POSIXlt(StopTime)$mon==6 & as.POSIXlt(StopTime)$mday==1)], 'StopTime_Acc', "NAccMD")
	set(df, df[, which(StopTime%in%as.Date(c("1911-01-01","1911-11-01","1911-11-11")))], 'StopTime_Acc', "NAccYMD")
	cat('\nNumber inaccurate StopTimes:')
	df[, table(StopTime_Acc)]
	#
	#	handle Ard's bogus values	#0: manual curation
	#	
	tmp		<- which(df[, Patient=="M26087" & StartTime=="2005-05-09" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M26087  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2008-12-09"))
	}
	tmp		<- which(df[, Patient=="M17324" & StartTime=="1993-02-10" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M17324  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1995-03-15"))
	}
	tmp		<- which(df[, Patient=="M17324" & StartTime=="1995-03-15" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M17324  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1996-04-17"))
	}
	tmp		<- which(df[, Patient=="M17324" & StartTime=="1996-04-17" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M17324  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1996-10-14"))
	}
	tmp		<- which(df[, Patient=="M25962" & StartTime=="1996-07-01" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M25962  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1999-04-15"))
	}
	tmp		<- which(df[, Patient=="M25962" & StartTime=="1999-04-15" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M25962  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2003-03-20"))
	}
	tmp		<- which(df[, Patient=="M26044" & StartTime=="2008-04-01" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M26044  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2008-09-25"))
	}
	tmp		<- which(df[, Patient=="M26044" & StartTime=="2008-09-25" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M26044  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2009-05-07"))
	}
	tmp		<- which(df[, Patient=="M31241" & StartTime=="1989-09-15" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M31241  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1994-01-15"))
	}
	tmp		<- which(df[, Patient=="M31241" & StartTime=="1994-01-15" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M31241  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1996-02-16"))
	}
	tmp		<- which(df[, Patient=="M33356" & StartTime=="1991-07-01" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M33356  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("1994-07-01"))
	}
	tmp		<- which(df[, Patient=="M33356" & StartTime=="1994-07-01" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M33356  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2005-07-01"))
	}
	tmp		<- which(df[, Patient=="M39099" & StartTime=="2002-07-01" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M39099  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2004-12-15"))
	}
	tmp		<- which(df[, Patient=="M39099" & StartTime=="2004-12-15" & StopTime=="1911-11-11"])
	if(length(tmp))
	{		 
		cat(paste("\nfix entry 	M39099  1911-11-11"))
		set(df, tmp, "StopTime", as.Date("2010-09-08"))
	}	
	#
	#	handle Ard's bogus values	#1
	#
	setkey(df, Patient, StopTime)
	tmp		<- c( 	df[, which(is.na(StartTime) & StopTime%in%as.Date(c("1911-01-01","1911-07-01","1911-11-01","1911-11-11")))],
			df[, which(StartTime%in%as.Date(c("1911-01-01","1911-07-01","1911-11-01","1911-11-11")) & StopTime%in%as.Date(c("1911-01-01","1911-11-01","1911-11-11")))]	)
	if(length(tmp))
	{
		cat('\nremove bogus entries with NA StartTime and "1911-01-01","1911-11-01","1911-07-01","1911-11-11" StopTime, n=', length(tmp))
		set(df, tmp, 'StopTime', NA_real_)
	}
	df		<- subset(df, !is.na(StopTime))
	#
	#	set all bogus values to 1911-11-11 for simplicity
	#
	set(df, df[,which(StartTime%in%as.Date(c("1911-01-01","1911-07-01","1911-11-01")))], 'StartTime',as.Date("1911-11-11"))
	set(df, df[,which(StopTime%in%as.Date(c("1911-01-01","1911-07-01","1911-11-01")))], 'StopTime',as.Date("1911-11-11")) 
	#
	#	handle Ard's bogus values	#2: set bogus values for last episodes to StartTime+1
	#
	df[, DUMMY:=seq_len(nrow(df))]
	tmp		<- df[, list(N=length(which(StopTime=="1911-11-11"))), by='Patient']	
	tmp		<- merge(subset(tmp, N>=1), df, by='Patient')
	tmp2	<- tmp[, {
				z		<- as.character(StartTime[ which(StopTime=="1911-11-11") ])
				z		<- z[ z==as.character(max( max(StartTime), max(StopTime) )) ]
				if(length(z)==0)
					z	<- NA_character_
				list(StartTime=z)
			}, by='Patient']
	tmp2	<- subset(tmp2, !is.na(StartTime))
	for(i in seq_len(nrow(tmp2)))
	{
		cat('\nSet 1911-11-11 of',tmp2[i,Patient],'to',as.character(tmp2[i,as.Date(StartTime)+1]))
		z	<- df[, which(Patient==tmp2[i,Patient] & StartTime==tmp2[i,StartTime])]
		stopifnot(length(z)==1)
		set(df, z, 'StopTime', tmp2[i,as.Date(StartTime)+1])		
	}
	#
	#	handle Ard's bogus values	#3: remove 1911-11-11 pairs by midpoint 
	#
	tmp		<- df[, list(Start_N=length(which(StartTime=="1911-11-11")), Stop_N=length(which(StopTime=="1911-11-11"))), by='Patient']		
	tmp		<- merge(subset(tmp, Start_N==1 & Stop_N%in%c(1,2)), df, by='Patient')
	tmp		<- tmp[, {
				z		<- StartTime[ which(StopTime=="1911-11-11") ]
				zz		<- which(StartTime=="1911-11-11")
				if(length(zz)>0)
				{
					z	<- z[ which.min(difftime(StopTime[zz], z, units='days')) ]					
				}
				if(length(zz)==0)
				{
					z	<- NA_character_
					zz	<- NA_integer_
				}
				list(StartTime=z, StopTime=StopTime[zz], IDX1=DUMMY[z==StartTime], IDX2=DUMMY[zz])
			}, by='Patient']
	for(i in seq_len(nrow(tmp)))
	{
		z	<- difftime(tmp[i, as.Date(StopTime)], tmp[i, as.Date(StartTime)], units='days')/2
		stopifnot(z>0)
		z	<- tmp[i, as.Date(StartTime)]+z
		cat('\nSet 1911-11-11 of',tmp[i,Patient],'to',as.character(z))
		set(df, tmp[i,IDX1], 'StopTime',z)
		set(df, tmp[i,IDX2], 'StartTime',z)				
	}
	#
	#	handle Ard's bogus values	#4: set 1911-11-11 to next StartTime (only after #2)
	#	
	tmp		<- df[, list(Start_N=length(which(StartTime=="1911-11-11")), Stop_N=length(which(StopTime=="1911-11-11"))), by='Patient']
	stopifnot( !nrow(subset(tmp, Stop_N>1)) )
	tmp		<- merge(subset(tmp, Start_N==0 & Stop_N==1), df, by='Patient')
	tmp		<- tmp[, {
				z	<- which(StopTime=="1911-11-11")
				zz	<- sort(StartTime)
				list(IDX1=DUMMY[z], StopTime=zz[ which(zz==StartTime[z])+1L ])		
			}, by='Patient']
	for(i in seq_len(nrow(tmp)))
	{		
		cat('\nSet 1911-11-11 of',tmp[i,Patient],'to',tmp[i,as.character(StopTime)])
		stopifnot( df[tmp[i,IDX1],StartTime]<=tmp[i,StopTime])
		set(df, tmp[i,IDX1], 'StopTime',tmp[i,StopTime])			
	}
	tmp		<- df[, list(Start_N=length(which(StartTime=="1911-11-11")), Stop_N=length(which(StopTime=="1911-11-11"))), by='Patient']
	stopifnot( !nrow(subset(tmp, Stop_N>0)) )
	df[, DUMMY:=NULL]
	#
	#	removing Patients not on treatment (this should be legacy)
	#
	tmp						<- which(df[, is.na(StartTime) & StopTime==TR.notyet & NoDrug==0])
	if(length(tmp))
	{
		cat("\nnumber of entries with is.na(StartTime) & StopTime==TR.notyet & NoDrug==0",length(tmp),"Set StopTime to NA")
		set(df, tmp, "StopTime", NA)		
	}
	tmp						<- which(df[, StartTime==TR.notyet & StopTime==TR.notyet])
	if(length(tmp))
	{		
		cat("\nnumber of entries with StartTime==TR.notyet & StopTime==TR.notyet",length(tmp),"Set StartTime and StopTime to NA")
		set(df, tmp, "StartTime", NA)
		set(df, tmp, "StopTime", NA)
	}
	tmp						<- which(df[, StartTime==TR.notyet])
	if(length(tmp))
	{		
		cat("\nnumber of entries with StartTime==TR.notyet & StopTime!=TR.notyet",length(tmp),"Misclassified StartTime. Set to NA")
		set(df, tmp, "StartTime", NA)
	}	
	df						<- subset(df, !is.na(StopTime))
	#
	#	check patients with all NoDrug==NA (legacy)
	#
	tmp			<- nrow(subset( df[, list(select=any(is.na(NoDrug))), by='Patient'], select))
	cat("\nnumber of patients with is.na(NoDrug), n=", tmp)
	stopifnot( !tmp	)
	#
	#	remove patients for which all periods have NoDrug==0
	#
	tmp			<- df[,  list(select=!all(NoDrug==0)) , by='Patient']	
	cat("\nPatients with all(NoDrug==0). Remove. n=",nrow(subset(tmp, !select)))
	df			<- merge( df, subset(tmp, select, Patient), by='Patient')
	#
	#	for each patient find first ART period with NoDrug>0 and keep only these
	#
	setkey(df, Patient, StopTime)	
	tmp			<- subset( df[, list( select= which(NoDrug>0)[1]>1 ), by='Patient'], select)[, Patient]
	cat(paste("\nPatients with NoDrug==0 in first episode(s). Remove these episodes. n=",length(tmp)))		
	df[, DUMMY:= seq_len(nrow(df))]
	tmp			<- df[, {					#need some akward code since StopTimes are not yet unique
				z<- seq.int( which(NoDrug>0)[1], length(NoDrug))
				list( DUMMY= DUMMY[z] )
			}, by='Patient']
	df			<- df[tmp[, DUMMY],]
	df[, DUMMY:=NULL]
	#	check that no NoDrug==0 at start
	tmp			<- subset( df[,list(check= any( is.na(StartTime) & NoDrug==0) ), by='Patient'], check, Patient)
	stopifnot( !nrow(tmp)	)
	tmp			<- subset( df[,list(check= NoDrug[1]==0 ), by='Patient'], check, Patient)
	stopifnot( !nrow(tmp)	)	
	#
	#	check Ards bogus values #5: set dangling StartTime to StopTime-1
	#
	tmp		<- df[, which(StartTime=="1911-11-11")]
	for( i in tmp)
		set(df, i, 'StartTime', df[i,StopTime-1])	
	#
	#	fix StartTime > StopTime: check if StopTime is StartTime
	#	
	tmp		<- data.table(IDX=which(df[, StartTime>StopTime]))	
	cat(paste("\nnumber of entries with StartTime>StopTime, n=",nrow(tmp)))
	df[, DUMMY:=seq_len(nrow(df))]
	tmp[, Patient:= df[tmp[,IDX],Patient]]
	tmp[, StartTime:= df[tmp[,IDX],StartTime]]
	tmp[, StopTime:= df[tmp[,IDX],StopTime]]
	for(i in seq_len(nrow(tmp)))
		if(nrow(subset(df, Patient==tmp[i,Patient] & StartTime==tmp[i,StopTime] & StopTime==tmp[i,StartTime])))
		{
			cat('\nStopTime is StartTime for patient',tmp[i,Patient])
			set(df, tmp[i,IDX], c('StartTime','StopTime'), NA_real_)
		}
	df		<- subset(df, !is.na(StartTime) | !is.na(StopTime))
	df[, DUMMY:=seq_len(nrow(df))]
	#
	#	fix suspicious StopTimes: patients still on ART? (TODO need Ard to confirm)
	#
	tmp		<- subset(df, StartTime==StopTime & (StartTime_Acc!='Acc' | StopTime_Acc!='Acc'), select=c(Patient, StopTime, DUMMY) )
	setnames(tmp, 'StopTime', 'OddStopTime')
	tmp		<- merge(subset(df, select=c(Patient, StopTime)), tmp, by='Patient')
	tmp		<- tmp[, list(OddStopTime=OddStopTime[1], DUMMY=DUMMY[1], LAST=max(StopTime)<=OddStopTime[1]), by='Patient']
	tmp		<- subset(tmp, LAST)
	#write.csv(tmp, file='~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/check/check_inaccurateStopTime_equal_StartTime_and_Last_2.csv')
	for(i in seq_len(nrow(tmp)))
	{
		set(df, tmp[i,DUMMY], 'StopTime', TR.notyet)
		set(df, tmp[i,DUMMY], 'StopTime_Acc', 'Acc')
	}
	#
	#	fix StartTime > StopTime: fix inaccurate StopTime
	#		
	tmp		<- data.table(IDX=which(df[, StartTime>StopTime]))	
	cat(paste("\nnumber of entries with StartTime>StopTime, n=",nrow(tmp)))	
	tmp[, Patient:= df[tmp[,IDX],Patient]]
	tmp[, StartTime_Acc:=df[tmp[,IDX],StartTime_Acc]]
	tmp[, StopTime_Acc:=df[tmp[,IDX],StopTime_Acc]]	
	stopifnot( !nrow(subset(tmp,  StartTime_Acc=='Acc' & StopTime_Acc=='Acc')))	
	tmp2	<- merge(subset(tmp,  StartTime_Acc!='NAccMD'), df, by='Patient')	
	tmp2	<- tmp2[, list(Start_N=length(which(StartTime==StopTime[IDX==DUMMY])), Stop_N=length(which(StopTime==StopTime[IDX==DUMMY]))), by=c('Patient','IDX')]
	stopifnot( !nrow(subset(tmp2, Start_N>1))	)
	tmp2	<- merge(subset(tmp2, Start_N>0), df, by='Patient')
	tmp2	<- tmp2[, {
				cat('\n',Patient)
				z		<- StopTime[ IDX==DUMMY ]				# that s the inaccurate query time to be changed, eg 2003-07-01
				zz		<- which(StartTime==z)
				stopifnot(length(zz)==1)						# below we only compare just one start time 
				z		<- StartTime[ which(StopTime==z) ]		# that s the start times of all stop times that equal the query time
				stopifnot(length(z)>0)	
				z		<- z[ z<StopTime[zz] ]					
				stopifnot(length(z)>0)
				z		<- as.character( z[ which.min(difftime(StopTime[zz], z, units='days')) ] )									
				list(StartTime=z, StopTime=StopTime[zz], IDX1=DUMMY[z==StartTime], IDX2=DUMMY[zz])
			}, by=c('Patient','IDX')]
	for(i in seq_len(nrow(tmp2)))
	{
		cat('\nSet',as.character( df[tmp2[i,IDX1],StopTime] ), as.character( df[tmp2[i,IDX2],StartTime] ),'of',tmp2[i,Patient])
		z	<- difftime(tmp2[i, as.Date(StopTime)], tmp2[i, as.Date(StartTime)], units='days')/2
		stopifnot(z>0)
		z	<- tmp2[i, as.Date(StartTime)]+z
		stopifnot( df[tmp2[i,IDX1],StopTime]==df[tmp2[i,IDX2],StartTime], df[tmp2[i,IDX1],StopTime_Acc!='Acc'], df[tmp2[i,IDX2],StartTime_Acc!='Acc'] )		
		set(df, tmp2[i,IDX1], 'StopTime',z)
		set(df, tmp2[i,IDX2], 'StartTime',z)				
	}
	#
	#	fix StartTime > StopTime: fix dangling inaccurate StopTimes of StartTime>StopTime 
	#			
	tmp		<- data.table(IDX=which(df[, StartTime>StopTime]))	
	cat(paste("\nnumber of entries with StartTime>StopTime, n=",nrow(tmp)))	
	tmp[, Patient:= df[tmp[,IDX],Patient]]
	tmp[, StartTime_Acc:=df[tmp[,IDX],StartTime_Acc]]
	tmp[, StopTime_Acc:=df[tmp[,IDX],StopTime_Acc]]	
	tmp2	<- merge(subset(tmp,  StartTime_Acc!='NAccMD'), df, by=c('Patient','StartTime_Acc','StopTime_Acc'))
	tmp2	<- tmp2[, list(Start_N=length(which(StartTime==StopTime[IDX==DUMMY])), Stop_N=length(which(StopTime==StopTime[IDX==DUMMY]))), by=c('Patient','IDX')]
	stopifnot( !nrow(subset(tmp2, Start_N>0))	)	
	tmp2	<- merge(tmp2, df, by='Patient')
	tmp2	<- tmp2[, {
				z		<- StopTime[ IDX==DUMMY ]
				z		<- which(StopTime==z)
				zz		<- sort(StartTime)
				if( which(zz==StartTime[z])<length(zz))
					ans	<- zz[ which(zz==StartTime[z])+1L ]		#set to next start time
				if( which(zz==StartTime[z])==length(zz))
					ans	<- zz[length(zz)]+1
				list(IDX1=DUMMY[z], StopTime=ans)		
			}, by='Patient']
	for(i in seq_len(nrow(tmp2)))
	{		
		cat('\nSet ',as.character( df[tmp2[i,IDX1],StopTime] ),'of',tmp2[i,Patient],'to',tmp2[i,as.character(StopTime)])
		stopifnot( df[tmp2[i,IDX1],StartTime]<=tmp2[i,StopTime])
		set(df, tmp2[i,IDX1], 'StopTime',tmp2[i,StopTime])			
	}
	#
	#	fix StartTime > StopTime: fix inaccurate StartTime
	#		
	tmp		<- data.table(IDX=which(df[, StartTime>StopTime]))	
	cat(paste("\nnumber of entries with StartTime>StopTime, n=",nrow(tmp)))	
	tmp[, Patient:= df[tmp[,IDX],Patient]]
	tmp[, StartTime_Acc:=df[tmp[,IDX],StartTime_Acc]]
	tmp[, StopTime_Acc:=df[tmp[,IDX],StopTime_Acc]]	
	stopifnot( !nrow(subset(tmp,  StartTime_Acc=='Acc' & StopTime_Acc=='Acc')))	
	tmp2	<- merge(subset(tmp,  StopTime_Acc!='NAccMD'), df, by='Patient')	
	tmp2	<- tmp2[, list(Start_N=length(which(StartTime==StartTime[IDX==DUMMY])), Stop_N=length(which(StopTime==StartTime[IDX==DUMMY]))), by=c('Patient','IDX')]
	stopifnot( !nrow(subset(tmp2, Start_N>1)), !nrow(subset(tmp2, Start_N==0))	)
	tmp2	<- merge(subset(tmp2, Stop_N>0), df, by='Patient')
	stopifnot( !nrow(subset(tmp2, Start_N>1)), !nrow(subset(tmp2, Start_N==0))	)
	tmp2	<- tmp2[, {				
				z		<- StartTime[ IDX==DUMMY ]
				zz		<- which(StartTime==z)
				stopifnot(length(zz)==1)	
				z		<- StartTime[ which(StopTime==z) ]
				stopifnot(length(z)>0)	
				z		<- z[ z<StopTime[zz] ]					
				stopifnot(length(z)>0)				
				z		<- as.character( z[ which.min(difftime(StopTime[zz], z, units='days')) ] )									
				list(StartTime=z, StopTime=StopTime[zz], IDX1=DUMMY[z==StartTime], IDX2=DUMMY[zz])
			}, by=c('Patient','IDX')]
	for(i in seq_len(nrow(tmp2)))
	{
		cat('\nSet',as.character( df[tmp2[i,IDX1],StopTime] ), as.character( df[tmp2[i,IDX2],StartTime] ),'of',tmp2[i,Patient])
		z	<- difftime(tmp2[i, as.Date(StopTime)], tmp2[i, as.Date(StartTime)], units='days')/2
		stopifnot(z>0)
		z	<- tmp2[i, as.Date(StartTime)]+z
		stopifnot( df[tmp2[i,IDX1],StopTime]==df[tmp2[i,IDX2],StartTime], df[tmp2[i,IDX1],StopTime_Acc!='Acc'], df[tmp2[i,IDX2],StartTime_Acc!='Acc'] )		
		set(df, tmp2[i,IDX1], 'StopTime',z)
		set(df, tmp2[i,IDX2], 'StartTime',z)				
	}
	#
	#	fix StartTime > StopTime: fix inaccurate StartTime
	#		
	tmp		<- data.table(IDX=which(df[, StartTime>StopTime]))
	stopifnot( !nrow(tmp))	
	df[, DUMMY:=NULL]
	#
	#	remove duplicate entries
	#
	setkey(df, Patient, StartTime, StopTime)
	tmp		<- which(duplicated(df))
	if(length(tmp))
	{
		cat("\nnumber of duplicate entries with key  Patient, StartTime, StopTime, n=",nrow(tmp))
		df		<- unique(df)		
	}
	#
	#	remove equal entries
	#	
	tmp			<- df[, which(StartTime==StopTime)]
	if(length(tmp))
	{
		cat("\nrm entries with StartTime==StopTime, n=",length(tmp))
		df		<- df[-tmp,]		
	}
	#	
	#	should have manually fixed double entries (these were just 6)
	#
	setkey(df, Patient, StopTime)
	tmp		<- data.table(IDX=which(duplicated(df)))	
	tmp[, Patient:= df[tmp[,IDX],Patient]]
	stopifnot( !nrow(tmp) )
	#
	#	see if we can merge consecutive NoDrug==0 periods
	#
	df[, DUMMY:= seq_len(nrow(df))]
	tmp		<- unique(subset(df, NoDrug==0, Patient))
	tmp		<- merge(tmp, df, by='Patient')	
	tmp		<- tmp[,	{
				z		<- which(NoDrug==0)
				ans		<- NA_integer_
				if(length(z)>1)
					ans	<- DUMMY[ z[ which(diff(z)==1)+1L ] ]
				list(IDX=ans)
			} ,by='Patient']
	tmp		<- subset(tmp, !is.na(IDX))	
	cat("\nnumber of no drug episodes with follow up no drug episodes, n=",nrow(tmp))
	for(i in seq_len(nrow(tmp)))
	{
		cat('\nmerge no drug episodes for patient', tmp[i, Patient])
		stopifnot(df[tmp[i,IDX]-1L, Patient]==df[tmp[i,IDX], Patient],  df[tmp[i,IDX], NoDrug]==0,  df[tmp[i,IDX]-1L, NoDrug]==0)
		set(df, tmp[i,IDX], 'StartTime', df[tmp[i,IDX]-1L, StartTime])
		set(df, tmp[i,IDX], 'StartTime_Acc', df[tmp[i,IDX]-1L, StartTime_Acc])
		set(df, tmp[i,IDX]-1L, c('StartTime','StopTime'), NA_real_)
	}
	df[, DUMMY:=NULL]
	df		<- subset(df, !is.na(StartTime) | !is.na(StopTime))
	#
	#	TR.interrupted
	#
	if(nrow(df[which(is.na(df[,NoDrug]) & StartTime!=TR.notyet),])) 
		stop("unexpected NA in NoDrug when on treatment")
	tmp								<- rep(0, nrow(df))
	tmp[ which( df[,NoDrug==0] ) ]	<- 1
	df[, TrI:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#
	#	process treatment change reasons
	#
	if(any(!is.na(df[, Reason9])))	stop("unexpected !NA after Reason7")
	df.TrCh.noreason	<-  df[, which(is.na(Reason1)&is.na(Reason2)&is.na(Reason3)&is.na(Reason4)&is.na(Reason5)&is.na(Reason6)&is.na(Reason7)&is.na(Reason8)&(NoDrug!=0)) ] 
	#	TR.addition
	tmp	<- na.omit(c( 	df[df.TrCh.noreason, AZT]<=df[df.TrCh.noreason+1, AZT] &
							df[df.TrCh.noreason, DDI]<=df[df.TrCh.noreason+1, DDI] &
							df[df.TrCh.noreason, DDC]<=df[df.TrCh.noreason+1, DDC] &
							df[df.TrCh.noreason, D4T]<=df[df.TrCh.noreason+1, D4T] &
							df[df.TrCh.noreason, DTC]<=df[df.TrCh.noreason+1, DTC] &
							df[df.TrCh.noreason, ABC]<=df[df.TrCh.noreason+1, ABC] &
							df[df.TrCh.noreason, TDF]<=df[df.TrCh.noreason+1, TDF] &
							df[df.TrCh.noreason, FTC]<=df[df.TrCh.noreason+1, FTC] &
							df[df.TrCh.noreason, SQV]<=df[df.TrCh.noreason+1, SQV] &
							df[df.TrCh.noreason, IDV]<=df[df.TrCh.noreason+1, IDV] & 
							df[df.TrCh.noreason, RTV]<=df[df.TrCh.noreason+1, RTV] &
							df[df.TrCh.noreason, NFV]<=df[df.TrCh.noreason+1, NFV] &
							df[df.TrCh.noreason, LPV]<=df[df.TrCh.noreason+1, LPV] &
							df[df.TrCh.noreason, FPV]<=df[df.TrCh.noreason+1, FPV] &
							df[df.TrCh.noreason, ATV]<=df[df.TrCh.noreason+1, ATV] &
							df[df.TrCh.noreason, TPV]<=df[df.TrCh.noreason+1, TPV] &
							df[df.TrCh.noreason, DRV]<=df[df.TrCh.noreason+1, DRV] &
							df[df.TrCh.noreason, NVP]<=df[df.TrCh.noreason+1, NVP] &
							df[df.TrCh.noreason, EFV]<=df[df.TrCh.noreason+1, EFV] &
							df[df.TrCh.noreason, ETR]<=df[df.TrCh.noreason+1, ETR] &
							df[df.TrCh.noreason, RPV]<=df[df.TrCh.noreason+1, RPV] &				
							df[df.TrCh.noreason, DOR]<=df[df.TrCh.noreason+1, DOR] &
							df[df.TrCh.noreason, ENF]<=df[df.TrCh.noreason+1, ENF] &
							df[df.TrCh.noreason, RAL]<=df[df.TrCh.noreason+1, RAL] &
							df[df.TrCh.noreason, MVC]<=df[df.TrCh.noreason+1, MVC] &
							df[df.TrCh.noreason, EVG]<=df[df.TrCh.noreason+1, EVG] &
							df[df.TrCh.noreason, DTG]<=df[df.TrCh.noreason+1, DTG] &
							df[df.TrCh.noreason, COBI]<=df[df.TrCh.noreason+1, COBI] 
			))	
	#tmp						<- df[df.TrCh.noreason, NoDrug+1]==df[df.TrCh.noreason+1, NoDrug]	
	df[, TrCh.addition:=0]
	set(df, df.TrCh.noreason[tmp], 'TrCh.addition', 1)
	df[, TrCh.addition:= factor(TrCh.addition,levels=c(0,1),labels=c("No","Yes"))]
	df.TrCh.noreason		<- df.TrCh.noreason[!tmp]
	#	TR.failure
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.failure | Reason2%in%TR.failure | Reason3%in%TR.failure | Reason4%in%TR.failure | Reason5%in%TR.failure | Reason6%in%TR.failure | Reason7%in%TR.failure | Reason8%in%TR.failure ]) ]<- 1
	df[, TrCh.failure:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#	TR.immufailure
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.immufailure | Reason2%in%TR.immufailure | Reason3%in%TR.immufailure | Reason4%in%TR.immufailure | Reason5%in%TR.immufailure | Reason6%in%TR.immufailure | Reason7%in%TR.immufailure | Reason8%in%TR.immufailure ]) ]<- 1
	df[, TrCh.immufailure:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]	
	#	TR.toxicity
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.toxicity | Reason2%in%TR.toxicity | Reason3%in%TR.toxicity | Reason4%in%TR.toxicity | Reason5%in%TR.toxicity | Reason6%in%TR.toxicity | Reason7%in%TR.toxicity | Reason8%in%TR.toxicity]) ]<- 1
	df[, TrCh.toxicity:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]		
	#	TR.adherence
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.adherence | Reason2%in%TR.adherence | Reason3%in%TR.adherence | Reason4%in%TR.adherence | Reason5%in%TR.adherence | Reason6%in%TR.adherence | Reason7%in%TR.adherence | Reason8%in%TR.adherence]) ]<- 1
	df[, TrCh.adherence:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#	TR.patient related
	tmp						<-	rep(0, nrow(df))
	tmp[ df.TrCh.noreason ]	<- NA
	tmp[ which(df[, Reason1%in%TR.patrel | Reason2%in%TR.patrel | Reason3%in%TR.patrel | Reason4%in%TR.patrel | Reason5%in%TR.patrel | Reason6%in%TR.patrel | Reason7%in%TR.patrel | Reason8%in%TR.patrel]) ]<- 1
	df[, TrCh.patrel:= factor(tmp,levels=c(0,1),labels=c("No","Yes"))]
	#
	#	set last stop time to start time - 1
	#
	df[, DUMMY:=seq_len(nrow(df))]
	tmp	<- df[, {
				ans		<- NA_integer_
				if(length(StartTime)>1)
					ans	<- DUMMY[ 1L+which(StartTime[-1]!=StopTime[-length(StopTime)]+1) ]
				list(IDX=ans)
			} , by='Patient']
	set(df, tmp[,IDX]-1L, 'StopTime', df[tmp[,IDX],StartTime]-1)	
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
						AnyT_T1		= subset(x, TrI=="No")[,StartTime][1], 
						AnyT_T1_Acc	= subset(x, TrI=="No")[,StartTime_Acc][1]		)																					
			}, by='Patient']			
	df		<- merge(subset(tmp,select=c(Patient,AnyT_T1,AnyT_T1_Acc,TrI.n, TrI.mo, TrI.p)), df, all.y=1, by="Patient")
	#
	#	reset AnyT_T1 conservatively
	#
	df[, AnyT_T1_Crude:=AnyT_T1]
	df	<- hivc.db.reset.ARTStartFromAccurate(df, col='AnyT_T1')			
	#	add largest lRNA at least 3 months after therapy start.
	if(0)
	{
		tmp					<- subset(df, select=c(Patient, AnyT_T1)) 
		setkey(tmp, Patient)
		tmp					<- merge(unique(tmp), df.viro, by='Patient')	
		tmp					<- merge(subset(df, select=c(Patient, AnyT_T1, StartTime, StopTime)), subset(tmp, difftime(PosRNA, AnyT_T1, units='days')>90,select=c(Patient, PosRNA, lRNA)), by='Patient', allow.cartesian=TRUE)
		tmp					<- subset(tmp, (StartTime<=PosRNA & PosRNA<=StopTime))
		tmp					<- tmp[, {
					z<- which.max(lRNA)
					list(AnyT_T1=AnyT_T1[z], PosRNA.mx= PosRNA[z], lRNA.mx=lRNA[z])
				}, by=c('Patient','StartTime','StopTime')]
		df					<- merge( df, subset(tmp, select=c(Patient, StartTime, StopTime, PosRNA.mx, lRNA.mx)), by=c('Patient','StartTime','StopTime'), all.x=TRUE)
		#	set failure by viral load during treatment period
		df[, TrVL.failure:=factor( df[, lRNA.mx>3], levels=c(FALSE, TRUE), labels=c('No','Yes'))]	
		#subset(df, TrVL.failure=='No' & TrCh.failure=='Yes')	viro failure can be indicated for VL < 1e3
		#subset(df, TrVL.failure=='No' & TrCh.failure=='Yes' & lRNA.mx<log10(200), select=c(Patient, StartTime, StopTime, AnyT_T1, PosRNA.mx, lRNA.mx))	#some are really suspicious!		
	}
	#
	df[, DUMMY:=NULL]	
	if(verbose) cat(paste("\nsave to", outfile))
	save(df, file=outfile)
}
######################################################################################
Regimen.CheckStartVL<- function(dir.name= DATA, verbose=1)
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
Regimen.CheckOnVL<- function(dir.name= DATA, verbose=1)
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
Regimen.CheckARTStartDate<- function()
{
	verbose		<- 1
	file		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_ART.rda'
	load(file)
	df.treatment<- df
	file		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_RNA.rda'
	load(file)
	df.v		<- df
	outfile		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_democlin/ATHENA_1610_All_ART_checked_DropVLBeforeART.rda'
	#
	#	check if VL decreases just before treatment start
	#
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
	
	check1		<- subset(check, lRNA.drop.before.ART>0.5 & lRNA.before.ART<4 & AnyT_T1_Acc=='NAccD' & drop.days2ART>20)
	if(verbose)	cat(paste("\nPatients lRNA.drop.before.ART>0.5 & AnyT_T1_Acc=='NAccD', n=",check1[,length(unique(Patient))]))
	print(check1[,unique(Patient)])
	#print(check1,n=500)
	#	reset AnyT_T1 to drop.day	
	reset		<- unique( subset(check1, select=c(Patient, AnyT_T1, drop.day)) )	
	set(reset, reset[,which(Patient=='M10607')], 'drop.day', as.Date('1992-09-15'))
	set(reset, reset[,which(Patient=='M13163')], 'drop.day', as.Date('1996-12-13'))
	set(reset, reset[,which(Patient=='M34362')], 'drop.day', as.Date('2007-11-29'))	
	set(reset, reset[,which(Patient=='M16397')], 'drop.day', as.Date('1998-07-21'))	
	set(reset, reset[,which(Patient=='M16420')], 'drop.day', as.Date('1997-08-04'))
	set(reset, reset[,which(Patient=='M16533')], 'drop.day', as.Date('1996-07-16'))
	set(reset, reset[,which(Patient=='M17326')], 'drop.day', as.Date('1996-07-16'))
	set(reset, reset[,which(Patient=='M17590')], 'drop.day', as.Date('1998-01-01'))	
	set(reset, reset[,which(Patient=='M19679')], 'drop.day', as.Date('2000-08-14'))
	set(reset, reset[,which(Patient=='M20670')], 'drop.day', as.Date('2002-03-01'))
	set(reset, reset[,which(Patient=='M27026')], 'drop.day', as.Date('2002-10-01'))
	set(reset, reset[,which(Patient=='M28553')], 'drop.day', as.Date('2003-06-02'))
	set(reset, reset[,which(Patient=='M29601')], 'drop.day', as.Date('2011-09-02'))
	set(reset, reset[,which(Patient=='M34816')], 'drop.day', as.Date('2007-09-25'))
	
	check2		<- subset(check, lRNA.drop.before.ART>0.5 & lRNA.before.ART<4 & AnyT_T1_Acc=='NAccMD' & drop.days2ART>20)
	if(verbose)	cat(paste("\nPatients lRNA.drop.before.ART>0.5 & AnyT_T1_Acc=='NAccMD', n=",check2[,length(unique(Patient))]))
	print(check2[,unique(Patient)])
	#print(check2,n=500)
	tmp			<- unique( subset(check2, select=c(Patient, AnyT_T1, drop.day)) )	
	set(tmp, tmp[,which(Patient=='M15637')], 'drop.day', as.Date('1998-08-24'))
	set(tmp, tmp[,which(Patient=='M15637')], 'drop.day', as.Date('1998-08-24'))
	set(tmp, tmp[,which(Patient=='M16547')], 'drop.day', as.Date('1998-08-08'))
	set(tmp, tmp[,which(Patient=='M17105')], 'drop.day', as.Date('1996-07-02'))
	set(tmp, tmp[,which(Patient=='M17999')], 'drop.day', as.Date('1999-07-18'))
	set(tmp, tmp[,which(Patient=='M18215')], 'drop.day', as.Date('1999-06-21'))
	set(tmp, tmp[,which(Patient=='M19512')], 'drop.day', as.Date('1999-09-16'))
	set(tmp, tmp[,which(Patient=='M26502')], 'drop.day', as.Date('1999-05-15'))
	set(tmp, tmp[,which(Patient=='M26798')], 'drop.day', as.Date('1999-08-14'))
	set(tmp, tmp[,which(Patient=='M30905')], 'drop.day', as.Date('2007-06-10'))
	set(tmp, tmp[,which(Patient=='M31720')], 'drop.day', as.Date('2008-06-16'))
	set(tmp, tmp[,which(Patient=='M32174')], 'drop.day', as.Date('2007-07-02'))
	set(tmp, tmp[,which(Patient=='M33446')], 'drop.day', as.Date('2003-04-03'))
	set(tmp, tmp[,which(Patient=='M34089')], 'drop.day', as.Date('2008-07-02'))
	set(tmp, tmp[,which(Patient=='M34413')], 'drop.day', as.Date('2011-07-10'))
	set(tmp, tmp[,which(Patient=='M37608')], 'drop.day', as.Date('2009-06-21'))
	set(tmp, tmp[,which(Patient=='M38858')], 'drop.day', as.Date('2010-07-16'))
	set(tmp, tmp[,which(Patient=='M40331')], 'drop.day', as.Date('2011-07-08'))	
	set(tmp, tmp[,which(Patient=='M42775')], 'drop.day', as.Date('2011-04-25'))
	set(tmp, tmp[,which(Patient=='M45852')], 'drop.day', as.Date('2015-08-08'))
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
		set(df.treatment, tmp, 'AnyT_T1_Acc', 'NAccD')
		set(df.treatment, tmp[1], 'StartTime_Acc', 'NAccD')		
	}
	df			<- df.treatment	
	save(df, file=outfile)
}
######################################################################################
CD4.process.QC<- function( df )
{
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
	subset(df, !is.na(CD4A))
}
######################################################################################
CD4.process.NA<- function( df, NA.time=c("","01/01/1911","11/11/1911","24/06/1923"), date.var=c("DateImm") )
{
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
	df
}
######################################################################################
GGD.160227<- function()
{
	#	INPUT FILE NAMES
	dir.name	<- '~/Dropbox (SPH Imperial College)/2015_ATHENA_May_Update'	
	file		<- paste(dir.name,"ATHENA_1502_All_Region_GGD_v02.csv",sep='/')
	df			<- read.csv(file, stringsAsFactors=FALSE)
	df			<- hivc.db.reset.Date(df, col="DateReg", NA.time=c("","01/01/1911","11/11/1911","24/06/1923"), date.format="%Y-%m-%d")
	df			<- data.table(df, key="Patient")
	setnames(df, c("DateReg","HospitalRegion"),c("GGD_Reg","Region_RegT1"))
	set(df, NULL, 'Region', df[, factor(Region, levels=c('AMSTERDAM','ROTTERDAM','REST'), labels=c('Amst','Rott','Other'))])
	set(df, NULL, 'GGD', 	df[, factor(GGD, 	levels=c(111, 		 			706, 				1009, 		 		1106,			1406, 				1906, 
									2006, 		 			2106, 				2209, 				2406, 			2506, 				2707, 
									3109, 					3406, 				3606, 				3906, 			4106, 				4506, 
									4607, 					4810, 				5006, 				5206, 			5406, 				5608, 
									6011, 					6106, 				7206, 				7306), 
							labels=c('Groningen',			'Drenthe',			'IJsselland',		'Twente',		'Gelre_IJssel',		'Hulpverlening_Gelderland_Midden',
									'Rivierenland',		'Nijmegen', 		'Flevoland',		'Utrecht',		'Midden_Nederland',	'Hollands_Noorden',
									'Kennemerland',		'Amsterdam',		'Gooi_Vechtstreek',	'Den_Haag',		'Zuid_Holland_West','Hollands_Midden',
									'Rotterdam_Rijnmond',	'Zuid_Holland_Zuid','Zeeland',			'West_Brabant',	'Hart_voor_Brabant', 'Brabant_Zuidoost',
									'Limburg-Noord',		'Zuid_Limburg',		'Fryslan',			'Zaanstreek_Waterland'))])
	set(df, NULL, 'Region_RegT1', df[, factor(Region_RegT1, levels=c(1,		2,		3,		4,		 5,		  6), 
							labels=c('Amst','North','East',	'South', 'West', 'Curuacao'))])					 
	#	exclude patients with all missing entries
	tmp			<- subset(df, !(!is.na(Region_RegT1) | !is.na(GGD_Reg) | !is.na(GGD) | !is.na(Region)))
	cat('\npatient with no info on GGD/Region', tmp[, Patient])	
	df			<- subset(df, !is.na(Region_RegT1) | !is.na(GGD_Reg) | !is.na(GGD) | !is.na(Region))
	
	#	remove duplicates
	set(df, df[, which(Patient=='M30073' & GGD_Reg=='2005-04-08' & is.na(Region))], 'GGD_Reg', NA_real_)
	df			<- subset(df, !is.na(GGD_Reg))
	setkey(df, Patient, GGD_Reg)
	tmp		<- data.table(IDX=which(duplicated(df)))
	tmp[, Patient:=df[tmp[,IDX],Patient]]
	print( merge(df, tmp, by='Patient') )
	set(df, tmp[, unique(IDX)], 'Patient', NA_character_)
	df			<- subset(df, !is.na(Patient))
	#	
	file		<- paste(substr(file, 1, nchar(file)-3),'R',sep='')	
	cat(paste("\nsave to", file))
	save(df, file=file)
}
######################################################################################
Population.161209<- function()
{
	indir	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/original_latest'
	chdir	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/check'
	inf.ggd	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/original_latest/CBS_Gemeenden_GGD.csv'
	outfile	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_democlin/CBS_1612_Population_GGD_Age.rda'
	#
	#	read GGD table
	#	
	fis		<- list.files(indir, pattern='CBS_Gemeenden_GGD_[0-9]+.csv', full.names=TRUE)
	dg		<- lapply(fis, function(fi){
						tmp		<- as.data.table(read.csv(fi, skip=6 ,header=FALSE, sep=';', quote='', col.names=c('COM','COM_ID','COM2','GGD_ID','GGD')))
						tmp		<- tmp[-nrow(tmp),]						
						tmp
					})
	dg		<- do.call('rbind', dg)	
	tmp		<- as.data.table(read.csv(file.path(indir, 'CBS_Gemeenden_GGD.csv'), header=FALSE, sep=';', quote='', col.names=c('COM','COM_ID','COM2','GGD_ID','GGD')))
	dg		<- rbind(dg,tmp)
	for(x in colnames(dg))
		set(dg, NULL, x, gsub(',','',gsub('"','',dg[[x]])))
	set(dg, NULL, 'GGD_ID', dg[, gsub('GG','',GGD_ID)])
	set(dg, NULL, 'COM_ID', dg[, gsub('GM','',COM_ID)])
	set(dg, NULL, 'COM2', NULL)
	set(dg, NULL, 'COM', dg[,gsub('[^[:alnum:]]','',COM)])
	set(dg, NULL, 'GGD', dg[,gsub('GGenGD','',gsub('Hulpverlening','',gsub('Hulpverleningsdienst','',gsub('GGD','',gsub('[^[:alnum:]]','',GGD)))))])
	#
	#	GGD table fixup corrupted entries
	#
	# dg[, sort(unique(GGD_ID))]
	tmp		<- dg[, which(grepl('^D|^ en',GGD_ID))]
	set(dg, tmp, 'GGD', dg[tmp,GGD_ID])
	set(dg, tmp, 'GGD_ID', NA_character_)
	set(dg, NULL, 'GGD', dg[,gsub('DG','',gsub('enG','',gsub('[^[:alnum:]]','',gsub('Regio','',gsub('DGD','',gsub(' en GD ','',gsub('D ','',GGD)))))))])
	#	fix spelling
	set(dg, NULL, 'GGD', dg[,	gsub('NoordLimburg|NoordenMiddenLimburg','LimburgNoord',
								gsub('GezondheidsdienstStadsgewestBreda','WestBrabant',
								gsub('DienstGezondheidJeugdZHZ','ZuidHollandZuid',
								gsub('Westfriesland','WestFriesland',
								gsub('ZuidoostBrabant','BrabantZuidoost',
								gsub('Delftland','ZuidHollandWest',
								gsub('VeiligheidsezondheidsregioGelderlandMidden|VeilighezGelderlandMidden','GelderlandMidden',
								gsub('MiddenNederland[0-9]','MiddenNederland',
								gsub('IJsselVecht','IJsselland',
								gsub('Delfland|Delftland',"ZuidHollandWest",
								gsub('MiddenHolland','HollandsMidden',
								gsub('KopvanNoordHolland','KopVanNoordHolland',
								gsub('WestFriesland|Frysln|Friesland','Fryslan',
								gsub('Achterhoek','GelreIJssel',
								gsub('AmstellanddeMeerlanden|AmstellandDeMeerlanden','Amsterdam',
								gsub('Rotterdameo','RotterdamRijnmond',GGD))))))))))))))))])
	#	update GGD assignments to current GGD
	dg.cur	<- c(106, 1, 11, 14, 44, 36, 56, 20,28,46,31,3107,59,50,26,19,111,706,1009,1106,1406,1906,2006,2106,2209,2406,2506,2707,3109,3406,3606,3906,4106,4506,4607,4810,5006,5206,5406,5608,6011,6106,7206,7306)
	dg		<- dg[, {
				z		<- sapply(as.numeric(GGD_ID), function(x) x%in%dg.cur)
				ggd		<- GGD
				ggd_id	<- GGD_ID
				if(any(z))
				{
					ggd_id[!z]	<- ggd[!z]	<- NA_character_
				}
				list(COM=COM, GGD_ID=ggd_id, GGD=ggd, ANY_CURRENT_GGD=any(z))
			}, by='COM_ID']
	dg		<- unique(subset(dg, !is.na(GGD)))
	#	check GGD manually by closest larger city / actual municipality
	#	subset(dg, !ANY_CURRENT_GGD)[1:50,]
	set(dg, dg[, which(COM%in%c('AmbtDelden'))], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='TerAar')], 'GGD', 'HollandsMidden')
	set(dg, dg[, which(COM=='TerAar')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM%in%c('Limmen','Akersloot'))], 'GGD', 'HollandsNoorden')
	set(dg, dg[, which(COM%in%c('Limmen','Akersloot'))], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Zelhem')], 'GGD', 'GelderlandMidden')
	set(dg, dg[, which(COM=='Zelhem')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Zwartsluis')], 'GGD', 'IJsselland')
	set(dg, dg[, which(COM=='Zwartsluis')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='BergenDal')], 'GGD', 'Nijmegen')
	set(dg, dg[, which(COM=='BergenDal')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Wognum')], 'GGD', 'HollandsNoorden')
	set(dg, dg[, which(COM=='Wognum')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='AmbtMontfort')], 'GGD', 'LimburgNoord')
	set(dg, dg[, which(COM=='AmbtMontfort')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Wisch')], 'GGD', 'IJsselland')
	set(dg, dg[, which(COM=='Wisch')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Zwolle')], 'GGD', 'IJsselland')
	set(dg, dg[, which(COM=='Zwolle')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Wehl')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Wehl')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Avereest')], 'GGD', 'IJsselland')
	set(dg, dg[, which(COM=='Avereest')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Bathmen')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Bathmen')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Bergh')], 'GGD', 'GelderlandMidden')
	set(dg, dg[, which(COM=='Bergh')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Borculo')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Borculo')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Eibergen')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Eibergen')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Neede')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Neede')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Ruurlo')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Ruurlo')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Born')], 'GGD', 'ZuidLimburg')
	set(dg, dg[, which(COM=='Born')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Geleen')], 'GGD', 'ZuidLimburg')
	set(dg, dg[, which(COM=='Geleen')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Brederwiede')], 'GGD', 'IJsselland')
	set(dg, dg[, which(COM=='Brederwiede')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Dinxperlo')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Dinxperlo')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Echt')], 'GGD', 'LimburgNoord')
	set(dg, dg[, which(COM=='Echt')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Susteren')], 'GGD', 'LimburgNoord')
	set(dg, dg[, which(COM=='Susteren')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Egmond')], 'GGD', 'LimburgNoord')
	set(dg, dg[, which(COM=='Egmond')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Schoorl')], 'GGD', 'LimburgNoord')
	set(dg, dg[, which(COM=='Schoorl')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Genemuiden')], 'GGD', 'IJsselland')
	set(dg, dg[, which(COM=='Genemuiden')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Gramsbergen')], 'GGD', 'IJsselland')
	set(dg, dg[, which(COM=='Gramsbergen')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='sGraveland')], 'GGD', 'GooienVechtstreek')
	set(dg, dg[, which(COM=='sGraveland')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='sGravendeel')], 'GGD', 'ZuidHollandZuid')
	set(dg, dg[, which(COM=='sGravendeel')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Harmelen')], 'GGD', 'MiddenNederland')
	set(dg, dg[, which(COM=='Harmelen')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Hasselt')], 'GGD', 'IJsselland')
	set(dg, dg[, which(COM=='Hasselt')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Leidschendam')], 'GGD', 'ZuidHollandWest')
	set(dg, dg[, which(COM=='Leidschendam')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Voorburg')], 'GGD', 'ZuidHollandWest')
	set(dg, dg[, which(COM=='Voorburg')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Liemeer')], 'GGD', 'HollandsMidden')
	set(dg, dg[, which(COM=='Liemeer')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Maartensdijk')], 'GGD', 'MiddenNederland')
	set(dg, dg[, which(COM=='Maartensdijk')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Ravenstein')], 'GGD', 'HartvoorBrabant')
	set(dg, dg[, which(COM=='Ravenstein')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Rijnsburg')], 'GGD', 'HollandsMidden')
	set(dg, dg[, which(COM=='Rijnsburg')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Sassenheim')], 'GGD', 'HollandsMidden')
	set(dg, dg[, which(COM=='Sassenheim')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Sittard')], 'GGD', 'ZuidLimburg')
	set(dg, dg[, which(COM=='Sittard')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Steenderen')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Steenderen')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Valkenburg')], 'GGD', 'HollandsMidden')
	set(dg, dg[, which(COM=='Valkenburg')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='VleutenDeMeern')], 'GGD', 'Utrecht')
	set(dg, dg[, which(COM=='VleutenDeMeern')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Voorhout')], 'GGD', 'HollandsMidden')
	set(dg, dg[, which(COM=='Voorhout')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Warmond')], 'GGD', 'HollandsMidden')
	set(dg, dg[, which(COM=='Warmond')], 'ANY_CURRENT_GGD', TRUE)		
	set(dg, dg[, which(COM=='Deventer')], 'GGD', 'IJsselland')
	set(dg, dg[, which(COM=='Deventer')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Heerde')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Heerde')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Andijk')], 'GGD', 'HollandsNoorden')
	set(dg, dg[, which(COM=='Andijk')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Hoorn')], 'GGD', 'HollandsNoorden')
	set(dg, dg[, which(COM=='Hoorn')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Medemblik')], 'GGD', 'HollandsNoorden')
	set(dg, dg[, which(COM=='Medemblik')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Castricum')], 'GGD', 'HollandsNoorden')
	set(dg, dg[, which(COM=='Castricum')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Castricum')], 'GGD_ID', NA_character_)
	set(dg, dg[, which(COM=='Enkhuizen')], 'GGD', 'HollandsNoorden')
	set(dg, dg[, which(COM=='Enkhuizen')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Opmeer')], 'GGD', 'HollandsNoorden')
	set(dg, dg[, which(COM=='Opmeer')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Wervershoof')], 'GGD', 'HollandsNoorden')
	set(dg, dg[, which(COM=='Wervershoof')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Drechterland')], 'GGD', 'HollandsNoorden')
	set(dg, dg[, which(COM=='Drechterland')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='StedeBroec')], 'GGD', 'HollandsNoorden')
	set(dg, dg[, which(COM=='StedeBroec')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='MookenMiddelaar')], 'GGD', 'LimburgNoord')
	set(dg, dg[, which(COM=='MookenMiddelaar')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='OlstWijhe')], 'GGD', 'IJsselland')
	set(dg, dg[, which(COM=='OlstWijhe')], 'ANY_CURRENT_GGD', TRUE)
	set(dg, dg[, which(COM=='Olst')], 'GGD', 'IJsselland')
	set(dg, dg[, which(COM=='Olst')], 'ANY_CURRENT_GGD', TRUE)		
	set(dg, dg[, which(COM=='Bergschenhoek')], 'GGD', 'ZuidHollandWest')
	set(dg, dg[, which(COM=='Bergschenhoek')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='BerkelenRodenrijs')], 'GGD', 'ZuidHollandWest')
	set(dg, dg[, which(COM=='BerkelenRodenrijs')], 'ANY_CURRENT_GGD', TRUE)		
	set(dg, dg[, which(COM=='Apeldoorn')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Apeldoorn')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Brummen')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Brummen')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Gorssel')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Gorssel')], 'ANY_CURRENT_GGD', TRUE)		
	set(dg, dg[, which(COM=='Epe')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Epe')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Vorden')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Vorden')], 'ANY_CURRENT_GGD', TRUE)		
	set(dg, dg[, which(COM=='Warnsveld')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Warnsveld')], 'ANY_CURRENT_GGD', TRUE)			
	set(dg, dg[, which(COM=='Lochem')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Lochem')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Voorst')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Voorst')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Zutphen')], 'GGD', 'GelreIJssel')
	set(dg, dg[, which(COM=='Zutphen')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(COM=='Bleiswijk')], 'GGD', 'RotterdamRijnmond')
	set(dg, dg[, which(COM=='Bleiswijk')], 'ANY_CURRENT_GGD', TRUE)	
	set(dg, dg[, which(GGD_ID=='31')], 'GGD', 'Kennemerland')
	set(dg, dg[, which(GGD_ID=='26')], 'GGD', 'MiddenNederland')
	set(dg, dg[, which(GGD_ID=='28')], 'GGD', 'HollandsNoorden')
	set(dg, dg[, which(GGD_ID=='19')], 'GGD', 'GelderlandMidden')
	set(dg, dg[, which(GGD_ID=='36')], 'GGD', 'GooienVechtstreek')		
	set(dg, dg[, which(GGD_ID=='29')], 'GGD', 'HollandsNoorden')	
	set(dg, dg[, which(GGD=='Groningen')], 'GGD_ID', '0111')
	set(dg, dg[, which(GGD=='Drenthe')], 'GGD_ID', '0706')
	set(dg, dg[, which(GGD=='IJsselland')], 'GGD_ID', '1009')
	set(dg, dg[, which(GGD=='Twente')], 'GGD_ID', '1106')
	set(dg, dg[, which(GGD=='GelreIJssel')], 'GGD_ID', '1406')
	set(dg, dg[, which(GGD=='GelderlandMidden')], 'GGD_ID', '1906')
	set(dg, dg[, which(GGD=='Rivierenland')], 'GGD_ID', '2006')
	set(dg, dg[, which(GGD=='Nijmegen')], 'GGD_ID', '2106')
	set(dg, dg[, which(GGD=='Flevoland')], 'GGD_ID', '2209')
	set(dg, dg[, which(GGD=='Utrecht')], 'GGD_ID', '2406')
	set(dg, dg[, which(GGD=='MiddenNederland')], 'GGD_ID', '2506')
	set(dg, dg[, which(GGD=='HollandsNoorden')], 'GGD_ID', '2707')
	set(dg, dg[, which(GGD=='Kennemerland')], 'GGD_ID', '3107')
	set(dg, dg[, which(GGD=='Amsterdam')], 'GGD_ID', '3406')
	set(dg, dg[, which(GGD=='GooienVechtstreek')], 'GGD_ID', '3606')
	set(dg, dg[, which(GGD=='DenHaag')], 'GGD_ID', '3906')
	set(dg, dg[, which(GGD=='ZuidHollandWest')], 'GGD_ID', '4106')
	set(dg, dg[, which(GGD=='HollandsMidden')], 'GGD_ID', '4506')
	set(dg, dg[, which(GGD=='RotterdamRijnmond')], 'GGD_ID', '4607')
	set(dg, dg[, which(GGD=='ZuidHollandZuid')], 'GGD_ID', '4810')
	set(dg, dg[, which(GGD=='Zeeland')], 'GGD_ID', '5006')
	set(dg, dg[, which(GGD=='WestBrabant')], 'GGD_ID', '5206')
	set(dg, dg[, which(GGD=='HartvoorBrabant')], 'GGD_ID', '5406')
	set(dg, dg[, which(GGD=='BrabantZuidoost')], 'GGD_ID', '5608')
	set(dg, dg[, which(GGD=='LimburgNoord')], 'GGD_ID', '6011')
	set(dg, dg[, which(GGD=='ZuidLimburg')], 'GGD_ID', '6106')
	set(dg, dg[, which(GGD=='Fryslan')], 'GGD_ID', '7206')
	set(dg, dg[, which(GGD=='ZaanstreekWaterland')], 'GGD_ID', '7306')
	set(dg, NULL, 'ANY_CURRENT_GGD', NULL)
	setkey(dg, COM_ID, COM, GGD_ID, GGD)
	dg		<- unique(dg)
	# 	find non-unique COM_ID
	dg		<- merge(dg, dg[, list(COM_N=length(COM)), by='COM_ID'], by='COM_ID')
	#	fixup
	set(dg, dg[, which(COM_ID=='0014')], 'COM', 'Groningen')
	set(dg, dg[, which(COM_ID=='0164')], 'COM', 'Hengelo')
	set(dg, dg[, which(COM_ID=='0248')], 'COM', 'HengeloGld')
	set(dg, dg[, which(COM_ID=='0344')], 'COM', 'Utrecht')
	set(dg, dg[, which(COM_ID=='0393')], 'COM', 'HaarlemmerliedeenSpaarnwoude')
	set(dg, dg[, which(COM_ID=='0417')], 'COM', 'Laren')
	set(dg, dg[, which(COM_ID=='0518')], 'COM', 'sGravenhage')
	set(dg, dg[, which(COM_ID=='0603')], 'COM', 'Rijswijk')
	set(dg, dg[, which(COM_ID=='0687')], 'COM', 'Middelburg')
	set(dg, dg[, which(COM_ID=='0820')], 'COM', 'NuenenGerwenenNederwetten')
	set(dg, dg[, which(COM_ID=='0888')], 'COM', 'Beek')
	set(dg, dg[, which(COM_ID=='0893')], 'COM', 'Bergen')
	set(dg, dg[, which(COM_ID=='0971')], 'COM', 'Stein')
	dg[, COM_N:=NULL]
	dg		<- merge(dg, dg[, list(COM_N=length(unique(GGD))), by='COM_ID'], by='COM_ID')
	stopifnot(nrow(subset(dg, COM_N>1))==0)
	dg[, COM_N:=NULL]
	setkey(dg, COM_ID, COM, GGD_ID, GGD)
	dg		<- unique(dg)
	#	
	#	AGE COUNTS	
	#
	fis		<- list.files(indir, pattern='CBS_Bevolking_[0-9]+_[0-9]+.csv', full.names=TRUE)
	df		<- lapply(fis, function(fi){
				z	<- as.data.table(read.csv(fi, header=FALSE, sep=ifelse(grepl('2004',fi),',',';'), quote=ifelse(grepl('2004',fi),'"',''), col.names=c('COM','YEAR','AGE','M','F','COM_ID')))
				for(x in colnames(z))
					set(z, NULL, x, gsub('"','',z[[x]]))
				z
			})
	df		<- do.call('rbind', df)	
	set(df, NULL, 'YEAR', df[, as.integer(YEAR)])
	set(df, NULL, 'M', df[, as.integer(M)])
	set(df, NULL, 'F', df[, as.integer(F)])	
	set(df, NULL, 'AGE', df[, gsub(' of ouder',' 100',gsub('Jonger dan ','0 ',gsub(' jaar','',gsub(' tot ',' ',AGE))))])
	set(df, NULL, 'AGE_END', df[, as.integer(gsub('[0-9]+ ','',AGE))-1])
	set(df, NULL, 'AGE_START', df[, as.integer(gsub('([0-9]+) [0-9]+','\\1',AGE))])
	set(df, NULL, 'AGE', df[, paste(AGE_START,'-',AGE_END,sep='')])		
	#	AGE COUNTS fixup corrupted entries
	df		<- subset(df, !is.na(M) | !is.na(F))
	set(df, NULL, 'COM', df[,gsub('[^[:alnum:]]','',COM)])	
	tmp		<- unique(subset(df, COM_ID=='' | is.na(COM_ID), COM))
	stopifnot(nrow(tmp)==0)
	tmp		<- df[, which(nchar(COM_ID)==1)]
	set(df, tmp, 'COM_ID', df[tmp, paste('000',COM_ID,sep='')])	
	tmp		<- df[, which(nchar(COM_ID)==2)]
	set(df, tmp, 'COM_ID', df[tmp, paste('00',COM_ID,sep='')])
	tmp		<- df[, which(nchar(COM_ID)==3)]
	set(df, tmp, 'COM_ID', df[tmp, paste('0',COM_ID,sep='')])	
	#	fixup COM names
	set(df, df[, which(COM_ID=='0014')], 'COM', 'Groningen')
	set(df, df[, which(COM_ID=='0164')], 'COM', 'Hengelo')	
	set(df, df[, which(COM_ID=='0344')], 'COM', 'Utrecht')	
	set(df, df[, which(COM_ID=='0417')], 'COM', 'Laren')
	set(df, df[, which(COM_ID=='0518')], 'COM', 'sGravenhage')
	set(df, df[, which(COM_ID=='0603')], 'COM', 'Rijswijk')
	set(df, df[, which(COM_ID=='0687')], 'COM', 'Middelburg')	
	set(df, df[, which(COM_ID=='0888')], 'COM', 'Beek')
	set(df, df[, which(COM_ID=='0893')], 'COM', 'Bergen')
	set(df, df[, which(COM_ID=='0971')], 'COM', 'Stein')
	set(df, df[, which(COM_ID=='0619')], 'COM', 'Valkenburg')
	stopifnot( !nrow(unique(subset(df, COM%in%setdiff(df[, COM], dg[,COM]), c(COM, COM_ID)))) )	
	stopifnot( length(setdiff(df[, COM_ID], dg[,COM_ID]))==0 )
	#
	#	MERGE
	#
	dpop	<- merge(dg, df, by=c('COM_ID','COM'))
	#
	#	further re-assignments of former municipalities to GGD of current municipality
	#
	set(dpop, dpop[, which(COM%in%c("Bergh","Didam","Zelhem","Wisch"))], 'GGD', 'GelreIJssel' )
	set(dpop, dpop[, which(COM%in%c("Bergh","Didam","Zelhem","Wisch"))], 'GGD_ID', '1406' )
	set(dpop, dpop[, which(COM%in%c("BerkelenRodenrijs","Bergschenhoek"))], 'GGD', 'RotterdamRijnmond' )
	set(dpop, dpop[, which(COM%in%c("BerkelenRodenrijs","Bergschenhoek"))], 'GGD_ID', '4607' )	
	save(dpop, file=outfile)
	write.table(dpop, file=gsub('rda','csv',outfile))
	
	#
	#	further checks
	#
	z<- subset(dpop, GGD=='GelderlandMidden' & YEAR==2004)[, sort(unique(COM))]
	zz<- subset(dpop, GGD=='GelderlandMidden' & YEAR==2005)[, sort(unique(COM))]
	setdiff(z, zz)	#	"Angerlo" "Bergh"   "Didam"   "Zelhem"
	#	Angerlo now part of Zevenaar
	#	Bergh Didam merged into Montferland in 2005
	#	Zelhem now part of Bronckhorst
	subset(dpop, grepl('Zevenaar',COM))[, unique(YEAR)]			#	Zevenaar has counts since 2000 but jump in 2005
	subset(dpop, grepl('Zevenaar',COM))[, list(M=sum(M), F=sum(F)),by='YEAR']
	subset(dpop, grepl('Montferland',COM))[, unique(YEAR)]		#	Montferland has counts since 2005
	subset(dpop, grepl('Bronckhorst',COM))[, unique(YEAR)]		#	Bronckhorst has counts since 2005
	# 	Montferland Bronckhorst in GelreIJssel
	
	z<- subset(dpop, GGD=='ZuidHollandZuid' & YEAR==2006)[, sort(unique(COM))]
	zz<- subset(dpop, GGD=='ZuidHollandZuid' & YEAR==2007)[, sort(unique(COM))]
	setdiff(z, zz)	#	"sGravendeel"
	#	sGravendeel now part of Binnenmaas also ZuidHollandZuid
	subset(dpop, grepl('Binnenmaas',COM))[, list(M=sum(M), F=sum(F)),by='YEAR']
	
	z<- subset(dpop, GGD=='RotterdamRijnmond' & YEAR==2006)[, sort(unique(COM))]
	zz<- subset(dpop, GGD=='RotterdamRijnmond' & YEAR==2007)[, sort(unique(COM))]
	setdiff(zz, z)	#	"Lansingerland"	
	#	BerkelenRodenrijs, Bleiswijk and Bergschenhoek
	#	BerkelenRodenrijs and Bergschenhoek are in ZuidHollandWest

	z<- subset(dpop, GGD=='IJsselland' & YEAR==2004)[, sort(unique(COM))]
	zz<- subset(dpop, GGD=='IJsselland' & YEAR==2005)[, sort(unique(COM))]
	setdiff(z, zz)	#	"Wisch"
	#  Wisch merged with Oude IJsselstreek

}
######################################################################################
GGD.161027<- function()
{
	#	INPUT FILE NAMES
	file		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/original_161102/SHM_1602_161102_OR_ALL_Region_GGD.csv'		
	outfile		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/preprocessed/ATHENA_1610_All_Region_GGD.rda'
	
	df			<- read.csv(file, stringsAsFactors=FALSE)
	df			<- hivc.db.reset.Date(df, col="REG_D", NA.time=c("","01/01/1911","11/11/1911","24/06/1923"), date.format="%d/%m/%Y")
	df			<- data.table(df, key="Patient")
	setnames(df, c("REG_D"),c("GGD_Reg"))
	set(df, NULL, 'Region', df[, factor(Region, levels=c('AMSTERDAM','ROTTERDAM','REST'), labels=c('Amst','Rott','Other'))])
	set(df, NULL, 'GGD', 	df[, factor(GGD, 	levels=c(111, 		 			706, 				1009, 		 		1106,			1406, 				1906, 
									2006, 		 			2106, 				2209, 				2406, 			2506, 				2707, 
									3109, 					3406, 				3606, 				3906, 			4106, 				4506, 
									4607, 					4810, 				5006, 				5206, 			5406, 				5608, 
									6011, 					6106, 				7206, 				7306), 
							labels=c('Groningen',			'Drenthe',			'IJsselland',		'Twente',		'Gelre_IJssel',		'Hulpverlening_Gelderland_Midden',
									'Rivierenland',		'Nijmegen', 		'Flevoland',		'Utrecht',		'Midden_Nederland',	'Hollands_Noorden',
									'Kennemerland',		'Amsterdam',		'Gooi_Vechtstreek',	'Den_Haag',		'Zuid_Holland_West','Hollands_Midden',
									'Rotterdam_Rijnmond',	'Zuid_Holland_Zuid','Zeeland',			'West_Brabant',	'Hart_voor_Brabant', 'Brabant_Zuidoost',
									'Limburg-Noord',		'Zuid_Limburg',		'Fryslan',			'Zaanstreek_Waterland'))])
	#	exclude patients with all missing entries
	tmp			<- subset(df, !(!is.na(GGD_Reg) | !is.na(GGD) | !is.na(Region)))
	cat('\npatient with no info on GGD/Region', tmp[, Patient])	
	df			<- subset(df, !is.na(GGD_Reg) | !is.na(GGD) | !is.na(Region))
	
	#	remove duplicates
	set(df, df[, which(Patient=='M30073' & GGD_Reg=='2005-04-08' & is.na(Region))], 'GGD_Reg', NA_real_)
	set(df, df[, which(Patient=='M39550' & GGD_Reg=='2010-12-15')], 'GGD_Reg', NA_real_)	
	df			<- subset(df, !is.na(GGD_Reg))
	setkey(df, Patient, GGD_Reg)
	tmp		<- data.table(IDX=which(duplicated(df)))
	tmp[, Patient:=df[tmp[,IDX],Patient]]
	print( merge(df, tmp, by='Patient') )
	set(df, tmp[, unique(IDX)], 'Patient', NA_character_)
	df			<- subset(df, !is.na(Patient))
	#			
	cat(paste("\nsave to", outfile))
	save(df, file=outfile)
}
######################################################################################
CombinedTable.160227<- function()
{	
	require(data.table)			
	#	CONSTANTS
	resume			<- 0
	verbose			<- 1	
	#
	#	INPUT FILES generated with "project.hivc.Excel2dataframe.XXX"
	#	
	dir.name		<- "~/Dropbox (SPH Imperial College)/2015_ATHENA_May_Update"
	file.seq		<- paste(dir.name,"ATHENA_1502_All_Sequences.R",sep='/')
	file.patient	<- paste(dir.name,"ATHENA_1502_All_Patient.R",sep='/')		
	file.primoSHM	<- paste(dir.name,"ATHENA_1501_PrimoSHM.csv",sep='/')
	file.viro		<- paste(dir.name,"ATHENA_1502_All_Viro.R",sep='/')
	file.immu		<- paste(dir.name,"ATHENA_1502_All_Immu.R",sep='/')
	file.treatment	<- paste(dir.name,"ATHENA_1502_All_Regimens.R",sep='/')
	file.ggd		<- paste(dir.name,"ATHENA_1502_All_Region_GGD_v02.R",sep='/')
	file.out		<- paste(dir.name,"ATHENA_1502_All_PatientKeyCovariates.R",sep='/')
	
	if(resume)												#//load if there is R Master data.table
	{
		options(show.error.messages = FALSE)		
		readAttempt	<-try(suppressWarnings(load(file.out)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresume file",file.out))
		if(!inherits(readAttempt, "try-error"))	str(df.all)
		options(show.error.messages = TRUE)		
	}
	if(!resume || inherits(readAttempt, "try-error"))		#else generate R Master data.table
	{
		#
		#	start with sequences
		#
		load(file.seq)			
		df		<- do.call('rbind',lapply( seq_along(df),	function(i)
						{
							z	<- data.table(df[[i]][,	c("Patient","DateRes","Subtype","SampleCode")], key="SampleCode")
							z[, GENE:=names(df)[i]]
							z
						}))
		df		<- dcast.data.table(df, Patient+DateRes+Subtype+SampleCode~GENE, value.var='GENE')
		set(df, NULL, "SampleCode", gsub(' ','', df[, SampleCode]))
		setnames(df, c("SampleCode","DateRes"), c("FASTASampleCode","PosSeqT"))
		#set(df, NULL, "FASTASampleCode", factor(df[,FASTASampleCode]))
		#set(df, NULL, "Patient", factor(df[,Patient]))		
		df.all	<- subset(df, select=c(FASTASampleCode,Patient,Subtype,PosSeqT))
		cat(paste("\nnumber of sequences found, n=", nrow(df.all)))
		#
		# add Patient data
		#
		cat(paste("\nadding patient data"))
		load(file.patient)		
		stopifnot( !length(setdiff( df.all[, unique(Patient)], df[, unique(Patient)] ))	)
		df[, Subtype:=NULL]		
		df.all	<- merge(df, df.all, all.x=1, by="Patient", allow.cartesian=TRUE)
		setnames(df.all, c("MyDateNeg1","MyDatePos1"), c("NegT","PosT"))
		set(df.all, NULL, c('Acute_Spec_1','Acute_Spec_2','Acute_Spec_3','Acute_Spec_4'), NULL)
		#	reset PosT_Acc=='No' conservatively to end of year / month
		df.all[, PosT_Crude:= PosT]	
		df.all	<- hivc.db.reset.inaccuratePosT(df.all, nacc.dy.dy= 30, nacc.mody.mo= 11, nacc.mody.dy= 31, verbose=1)
		#	reset NegT_Acc=='No' conservatively to start of year / month
		df.all[, NegT_Crude:= NegT]
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
		cat("\nchecking manually clashes in NegT and PosT, n=", nrow(tmp))		
		#	wrong NegT	M28679 M30745 M31091 M31091	
		set(df.all, df.all[, which(Patient%in%c('M28679','M30745','M31091','M31091'))],'NegT',NA_integer_)
		set(df.all, df.all[, which(Patient%in%c('M28679','M30745','M31091','M31091'))],'NegT_Acc',NA_character_)
		df.all	<- df.all[, setdiff( colnames(df.all), c("PoslRNA_T1","lRNA_T1")), with=FALSE]
		#
		#	set AnyPos_T1
		#
		# 	add preliminary AnyPos_T1	-- since PosT conservative we are good to set irrespective of PosT_Acc
		tmp		<- df.all[, {
					z		<- NA_real_
					if( any(!is.na(PosSeqT)) )
						z	<- min(PosSeqT, na.rm=TRUE)
					list(AnyPos_T1=z)	
				},by='Patient']
		df.all	<- merge(df.all, tmp, by='Patient')
		tmp		<- which( df.all[, !is.na(PosT) & ( is.na(PosSeqT) | PosT<AnyPos_T1 ) ] )
		cat(paste("\nbuilding prel AnyPos_T1. Number of entries with !is.na(PosT) & ( is.na(PosSeqT) | PosT<AnyPos_T1 ), n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, PosT])
		tmp		<- which( df.all[, !is.na(AnyT_T1) & !is.na(AnyPos_T1) & AnyT_T1<AnyPos_T1 ] )
		cat(paste("\nbuilding prel AnyPos_T1. Number of entries with !is.na(AnyT_T1) & !is.na(AnyPos_T1) & AnyT_T1<AnyPos_T1 ), n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, AnyT_T1])		
		tmp		<- df.all[, which(!is.na(AnyT_T1) & is.na(AnyPos_T1)) ] 
		cat(paste("\nbuilding prel AnyPos_T1. Number of entries with !is.na(AnyT_T1) & is.na(AnyPos_T1), n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, AnyT_T1])						
		#	check for invalid NegT and set to NA 
		#	-- we would only know that NegT is before AnyPosT and this is not helpful
		df.all	<- hivc.db.resetNegTbyAnyPosT(df.all)				
		#
		#	add first RNA Virology date
		#
		cat(paste("\nadding virology data"))
		load(file.viro)
		#	calculate first RNA above log10(500)
		tmp		<- df[, {
					ans	<- NA_real_
					z	<- which(RNA>500)
					if(length(z))
						ans	<- PosRNA[min(z)]
					list(PoslRNAg500_T1=ans)
				}, by='Patient']
		setkey(df, Patient)		
		tmp		<- merge(unique(df), tmp, by='Patient')
		df.all	<- merge(df.all, subset(tmp, select=c(Patient, PoslRNA_T1, lRNA_T1,PoslRNAg500_T1)), by='Patient')		
		tmp		<-  df.all[, which(!is.na(PoslRNA_T1) & lRNA_T1>log10(500) &  ( is.na(AnyPos_T1) | PoslRNA_T1<AnyPos_T1 )) ]
		
		cat(paste("\nbuilding prel AnyPos_T1. Number of entries with !is.na(PoslRNA_T1) & lRNA_T1>log10(500) &  PoslRNA_T1<AnyPos_T1, n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, PoslRNA_T1])
		tmp		<-  df.all[, which(!is.na(PoslRNA_T1) & lRNA_T1<=log10(500) &  !is.na(PoslRNAg500_T1) & ( is.na(AnyPos_T1) | PoslRNAg500_T1<AnyPos_T1 )) ]
		cat(paste("\nbuilding prel AnyPos_T1. Number of entries with !is.na(PoslRNA_T1) & lRNA_T1<=log10(500) &  !is.na(PoslRNAg500_T1) & PoslRNAg500_T1<AnyPos_T1, n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, PoslRNAg500_T1])		
		#		
		#	checking manually NegT>PosRNA
		#
		tmp		<- subset(df.all, NegT>AnyPos_T1)[, unique(Patient)]
		set(df.all, df.all[, which(Patient%in%tmp)], 'NegT', NA_integer_)
		set(df.all, df.all[, which(Patient%in%tmp)], 'NegT_Acc', NA_character_)
		#
		#	add Treatment dates
		#
		cat(paste("\nadding treatment data"))
		load(file.treatment)				
		tmp			<- subset(df, select=c(Patient, AnyT_T1, AnyT_T1_Acc, AnyT_T1_Crude, TrI.n, TrI.p ))
		setkey(tmp, Patient)
		df.all[, AnyT_T1:=NULL]
		df.all		<- merge(df.all, unique(tmp), all.x=1, by="Patient")
		#
		#	add CD4 count data
		#
		cat(paste("\nadding CD4 data"))
		load(file.immu)		
		if(length(which(df[,is.na(PosCD4_T1)])))	stop("unexpected NA in PosCD4_T1")
		if(length(which(df[, is.na(PosCD4)]))) stop("unexpected NA in PosCD4")
		if(length(which(df[, is.na(CD4)]))) stop("unexpected NA in CD4")
		df.all[, idx:= seq_len(nrow(df.all))]
		df.cross	<- merge( subset(df.all, select=c(idx,FASTASampleCode,Patient,AnyPos_T1,PosT,PosSeqT,PoslRNA_T1,NegT,NegT_Acc)), subset(df, select=c(Patient,PosCD4,CD4)), allow.cartesian=T, by="Patient" )
		# 	ignore CD4 entries where NegT>PosCD4
		df.cross	<- hivc.db.resetCD4byNegT(df.cross, with.NegT_Acc.No=1, verbose=1)
		#	compute CD4_T1 and CD4_TS -- do not use on df[,CD4_T1] because some CD4 measurements might correspond to healthy patients
		#	ignore CD4 entries with CD4>500 and PosCD4<AnyPos_T1
		tmp		<- hivc.db.getCD4.T1andTS(df.cross, CD4.HIVNeg.min= 500)
		tmp		<- subset(tmp, select=c(idx, PosCD4_T1, CD4_T1))
		df.all	<- merge(df.all, tmp, by="idx", all.x=1)
		# manually checked remaining PosCD4_T1 < AnyPos_T1
		tmp		<-  df.all[, which(!is.na(PosCD4_T1) & (is.na(AnyPos_T1) | PosCD4_T1<AnyPos_T1) ) ]
		cat(paste("\nbuilding prel AnyPos_T1. Number of patients with !is.na(PosCD4_T1) & PosCD4_T1<AnyPos_T1, n=",df.all[tmp,][, length(unique(Patient))]))
		print( subset(df.all, PosCD4_T1<AnyPos_T1) )
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, PosCD4_T1])
		#				
		#	reset Acute fields
		#
		df.all[, DUMMY:= df.all[, as.numeric(difftime(AnyPos_T1, NegT, units="days"))/365]]				
		tmp		<- df.all[, which(!is.na(DUMMY) & DUMMY<=1)]		
		set( df.all, tmp, 'Acute_Spec', 'NEGTEST_le12m' )
		set(df.all, df.all[, which(Acute_Spec%in%c('CLIN','LAB','NEGTEST_le12m'))], 'isAcute', 'Yes')				
		df.all[, DUMMY:=NULL]
		df.all[, idx:=NULL]
		#
		#	add participation in primo-SHM
		#
		df.pr	<- as.data.table(read.csv(file.primoSHM, stringsAsFactors=FALSE, header=FALSE, col.names='Patient'))
		df.pr[, PrimoSHM:='Yes']
		df.all	<- merge(df.all, df.pr, by='Patient', all.x=1)
		set(df.all, df.all[, which(is.na(PrimoSHM))], 'PrimoSHM', 'No' )
		set(df.all, NULL, 'PrimoSHM', df.all[, factor(PrimoSHM, levels=c('No','Yes'), labels=c('No','Yes'))])
		#
		#	get GGD at time of diagnosis
		#
		load(file.ggd)		
		tmp		<- subset(df, select=c(Patient, GGD_Reg))
		setkey(tmp, Patient, GGD_Reg)
		tmp		<- df[, {
					list(GGD_Reg=GGD_Reg, GGD_Last=c(GGD_Reg[-1]-1, as.Date('2015-05-01')))
				}, by='Patient']
		df		<- merge(df, tmp, by=c('Patient','GGD_Reg'))		
		tmp		<- subset(df.all, select=c(Patient, AnyPos_T1))
		tmp		<- merge(unique(tmp), df, by='Patient', allow.cartesian=TRUE)
		setkey(tmp, Patient, GGD_Reg)
		tmp		<- tmp[,	{
					z		<- which(GGD_Reg<=AnyPos_T1 & AnyPos_T1<GGD_Last)
					z2		<- 'At_Diag'
					if(length(z)==0)
					{
						z	<- which.min(difftime(GGD_Reg, AnyPos_T1))
						z2	<- 'Earliest_After_Diag'
					}
					list(Region_RegT1=Region_RegT1, RegionGGD_AnyPosT1_Type=z2, GGD_RecordTime=GGD_Reg[z[1]], Region_AnyPosT1=Region[z[1]], GGD_AnyPosT1=GGD[z[1]])
				} , by='Patient']	
		#	
		stopifnot( !nrow(subset(tmp, RegionGGD_AnyPosT1_Type=='At_Diag' & (Region_RegT1=='Amst' & Region_RegT1!='Amst') | (Region_RegT1!='Amst' & Region_RegT1=='Amst'))) )
		z		<- tmp[, which(is.na(Region_RegT1))]
		set(tmp, z, 'Region_RegT1', tmp[z, GGD_AnyPosT1])
		set(tmp, tmp[, which(Region_RegT1%in%c('Groningen','Drenthe','Fryslan'))], 'Region_RegT1', 'North')
		set(tmp, tmp[, which(Region_RegT1%in%c('IJsselland','Twente','Gelre_IJssel','Hulpverlening_Gelderland_Midden','Flevoland',"Rivierenland","Nijmegen","Midden_Nederland"))], 'Region_RegT1', 'East')		
		set(tmp, tmp[, which(Region_RegT1%in%c('Amsterdam'))], 'Region_RegT1', 'Amst')
		set(tmp, tmp[, which(Region_RegT1%in%c('Utrecht','Gooi_Vechtstreek','Zuid_Holland_Zuid','Zeeland','Rotterdam_Rijnmond','Den_Haag','Hollands_Midden','Kennemerland','Zaanstreek_Waterland','Hollands_Noorden',"Zuid_Holland_West"))], 'Region_RegT1', 'West')
		set(tmp, tmp[, which(Region_RegT1%in%c('West_Brabant','Hart_voor_Brabant','Brabant_Zuidoost','Limburg-Noord','Zuid_Limburg'))], 'Region_RegT1', 'South')
		set(tmp, NULL, 'Region_RegT1', tmp[, factor(as.character(Region_RegT1))])		
		cat('\nPatients with no entry in GGD file', setdiff( df.all[, Patient], tmp[, Patient] ))		
		df.all	<- merge(df.all, tmp, by='Patient', all.x=1)
		#
		#	re-arrange cols
		#
		df.all	<- df.all[, c(	"Patient", "DateBorn", "Sex", "CountryBorn", "RegionOrigin", "DateLastContact", "DateDied", "isDead","Region_now","GGD_now",      
						"NegT","NegT_Acc","NegT_Crude","ACS",
						"AnyPos_T1","CountryInfection","Trm","isAcute","Acute_Spec","RegionGGD_AnyPosT1_Type","GGD_RecordTime","Region_AnyPosT1","GGD_AnyPosT1","Region_RegT1",
						"DateInCare", "FirstMed","PoslRNA_T1", "lRNA_T1","PoslRNAg500_T1","PosCD4_T1","CD4_T1",
						"AnyT_T1","AnyT_T1_Acc","AnyT_T1_Crude","PrimoSHM","TrI.n","TrI.p","DateAIDS",					
						"PosSeqT","FASTASampleCode","Subtype"), with=FALSE]          
		if(verbose)	cat(paste("\nsave to file",file.out))
		stopifnot( !nrow(subset(df.all, is.na(AnyPos_T1))) )
		save(df.all,file=file.out)
	}	
}
######################################################################################
CombinedTable.161027<- function()
{	
	require(data.table)			
	#	CONSTANTS
	resume			<- 0
	verbose			<- 1	
	#
	#	INPUT FILES generated with "project.hivc.Excel2dataframe.XXX"
	#	
	dir.name		<- "~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_democlin"
	file.seq		<- "~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_sequences_latest/ATHENA_1610_Sequences_withCOMET_withLANL.rda"
	file.patient	<- paste(dir.name,"ATHENA_1610_All_Patient.rda",sep='/')		
	file.primoSHM	<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/original_latest/ATHENA_1501_PrimoSHM.csv'
	file.viro		<- paste(dir.name,"ATHENA_1610_All_RNA.rda",sep='/')
	file.immu		<- paste(dir.name,"ATHENA_1610_All_CD4.rda",sep='/')
	file.treatment	<- paste(dir.name,"ATHENA_1610_All_ART_checked_DropVLBeforeART.rda",sep='/')
	file.ggd		<- paste(dir.name,"ATHENA_1610_All_Region_GGD.rda",sep='/')
	file.out		<- paste(dir.name,"ATHENA_1610_All_PatientKeyCovariates.rda",sep='/')
	file.out.numeric<- paste(dir.name,"ATHENA_1610_All_PatientKeyCovariates_Numeric.rda",sep='/')
	if(resume)												#//load if there is R Master data.table
	{
		options(show.error.messages = FALSE)		
		readAttempt	<-try(suppressWarnings(load(file.out)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresume file",file.out))
		if(!inherits(readAttempt, "try-error"))	str(df.all)
		options(show.error.messages = TRUE)		
	}
	if(!resume || inherits(readAttempt, "try-error"))		#else generate R Master data.table
	{
		#
		#	start with sequences
		#
		load(file.seq)		
		df.all	<- subset(ds, select=c(FASTASampleCode,Patient,SEQ_L,SUBTYPE_C,PosSeqT))
		cat(paste("\nnumber of sequences found, n=", nrow(df.all)))
		#
		# add Patient data
		#
		cat(paste("\nadding patient data"))
		load(file.patient)		
		stopifnot( !length(setdiff( df.all[, unique(Patient)], df[, unique(Patient)] ))	)
		df[, Subtype:=NULL]		
		df.all	<- merge(df, df.all, all.x=1, by="Patient", allow.cartesian=TRUE)
		setnames(df.all, c("MyDateNeg1","MyDatePos1"), c("NegT","PosT"))
		set(df.all, NULL, c('Acute_Spec_1','Acute_Spec_2','Acute_Spec_3','Acute_Spec_4'), NULL)
		#	reset PosT_Acc=='No' conservatively to end of year / month
		df.all[, PosT_Crude:= PosT]	
		df.all	<- hivc.db.reset.inaccuratePosT(df.all, nacc.dy.dy= 30, nacc.mody.mo= 11, nacc.mody.dy= 31, verbose=1)
		#	reset NegT_Acc=='No' conservatively to start of year / month
		df.all[, NegT_Crude:= NegT]
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
		cat("\nchecking manually clashes in NegT and PosT, n=", nrow(tmp))		
		#	wrong NegT	M28679 M30745 M31091 M38032	M45869
		set(df.all, df.all[, which(Patient%in%c('M28679','M30745','M31091','M38032','M45869'))],'NegT',NA_integer_)
		set(df.all, df.all[, which(Patient%in%c('M28679','M30745','M31091','M38032','M45869'))],'NegT_Acc',NA_character_)		
		df.all	<- df.all[, setdiff( colnames(df.all), c("PoslRNA_T1","lRNA_T1")), with=FALSE]		
		#
		#	set AnyPos_T1
		#
		# 	add preliminary AnyPos_T1	-- since PosT conservative we are good to set irrespective of PosT_Acc
		tmp		<- df.all[, {
					z		<- NA_real_
					if( any(!is.na(PosSeqT)) )
						z	<- min(PosSeqT, na.rm=TRUE)
					list(AnyPos_T1=z)	
				},by='Patient']
		df.all	<- merge(df.all, tmp, by='Patient')
		#	for trm through birth, use date born
		tmp		<- which( df.all[, Trm=='PREG' & !is.na(NegT) ] )
		set(df.all, tmp, 'NegT', NA_integer_)
		set(df.all, tmp, 'NegT_Acc', NA_character_)
		tmp		<- which( df.all[, Trm=='PREG' ] )
		set(df.all, tmp, 'AnyPos_T1', df.all[tmp, DateBorn])
		#	update if PosSeqT earlier
		tmp		<- which( df.all[, !is.na(PosT) & ( is.na(PosSeqT) | PosT<AnyPos_T1 ) ] )
		cat(paste("\nbuilding prel AnyPos_T1. Number of entries with !is.na(PosT) & ( is.na(PosSeqT) | PosT<AnyPos_T1 ), n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, PosT])
		#	update if AnyT_T1 earlier
		tmp		<- which( df.all[, !is.na(AnyT_T1) & !is.na(AnyPos_T1) & AnyT_T1<AnyPos_T1 ] )
		cat(paste("\nbuilding prel AnyPos_T1. Number of entries with !is.na(AnyT_T1) & !is.na(AnyPos_T1) & AnyT_T1<AnyPos_T1 ), n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, AnyT_T1])		
		tmp		<- df.all[, which(!is.na(AnyT_T1) & is.na(AnyPos_T1)) ] 
		cat(paste("\nbuilding prel AnyPos_T1. Number of entries with !is.na(AnyT_T1) & is.na(AnyPos_T1), n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, AnyT_T1])						
		#	check for invalid NegT and set to NA 
		#	-- we would only know that NegT is before AnyPosT and this is not helpful
		df.all	<- hivc.db.resetNegTbyAnyPosT(df.all)				
		#
		#	add first RNA Virology date
		#
		cat(paste("\nadding virology data"))
		load(file.viro)
		#	calculate first RNA above log10(500)
		tmp		<- df[, {
					ans	<- NA_real_
					z	<- which(RNA>500)
					if(length(z))
						ans	<- PosRNA[min(z)]
					list(PoslRNAg500_T1=ans)
				}, by='Patient']
		setkey(df, Patient)		
		tmp		<- merge(unique(df), tmp, by='Patient')
		df.all	<- merge(df.all, subset(tmp, select=c(Patient, PoslRNA_T1, lRNA_T1,PoslRNAg500_T1)), by='Patient')		
		tmp		<-  df.all[, which(!is.na(PoslRNA_T1) & lRNA_T1>log10(500) &  ( is.na(AnyPos_T1) | PoslRNA_T1<AnyPos_T1 )) ]
		
		cat(paste("\nbuilding prel AnyPos_T1. Number of entries with !is.na(PoslRNA_T1) & lRNA_T1>log10(500) &  PoslRNA_T1<AnyPos_T1, n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, PoslRNA_T1])
		tmp		<-  df.all[, which(!is.na(PoslRNA_T1) & lRNA_T1<=log10(500) &  !is.na(PoslRNAg500_T1) & ( is.na(AnyPos_T1) | PoslRNAg500_T1<AnyPos_T1 )) ]
		cat(paste("\nbuilding prel AnyPos_T1. Number of entries with !is.na(PoslRNA_T1) & lRNA_T1<=log10(500) &  !is.na(PoslRNAg500_T1) & PoslRNAg500_T1<AnyPos_T1, n=",length(tmp)))
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, PoslRNAg500_T1])		
		#		
		#	checking manually NegT>PosRNA
		#
		tmp		<- subset(df.all, NegT>AnyPos_T1)[, unique(Patient)]
		set(df.all, df.all[, which(Patient%in%tmp)], 'NegT', NA_integer_)
		set(df.all, df.all[, which(Patient%in%tmp)], 'NegT_Acc', NA_character_)
		#
		#	add Treatment dates
		#
		cat(paste("\nadding treatment data"))
		load(file.treatment)				
		tmp			<- subset(df, select=c(Patient, AnyT_T1, AnyT_T1_Acc, AnyT_T1_Crude, TrI.n, TrI.p ))
		setkey(tmp, Patient)
		df.all[, AnyT_T1:=NULL]
		df.all		<- merge(df.all, unique(tmp), all.x=1, by="Patient")
		#
		#	add CD4 count data
		#
		cat(paste("\nadding CD4 data"))
		load(file.immu)		
		if(length(which(df[,is.na(PosCD4_T1)])))	stop("unexpected NA in PosCD4_T1")
		if(length(which(df[, is.na(PosCD4)]))) stop("unexpected NA in PosCD4")
		if(length(which(df[, is.na(CD4)]))) stop("unexpected NA in CD4")
		df.all[, idx:= seq_len(nrow(df.all))]
		df.cross	<- merge( subset(df.all, select=c(idx,FASTASampleCode,Patient,AnyPos_T1,PosT,PosSeqT,PoslRNA_T1,NegT,NegT_Acc)), subset(df, select=c(Patient,PosCD4,CD4)), allow.cartesian=T, by="Patient" )
		# 	ignore CD4 entries where NegT>PosCD4
		df.cross	<- hivc.db.resetCD4byNegT(df.cross, with.NegT_Acc.No=1, verbose=1)
		#	compute CD4_T1 and CD4_TS -- do not use on df[,CD4_T1] because some CD4 measurements might correspond to healthy patients
		#	ignore CD4 entries with CD4>500 and PosCD4<AnyPos_T1
		tmp		<- hivc.db.getCD4.T1andTS(df.cross, CD4.HIVNeg.min= 500)
		tmp		<- subset(tmp, select=c(idx, PosCD4_T1, CD4_T1))
		df.all	<- merge(df.all, tmp, by="idx", all.x=1)
		# manually checked remaining PosCD4_T1 < AnyPos_T1
		tmp		<-  df.all[, which(!is.na(PosCD4_T1) & (is.na(AnyPos_T1) | PosCD4_T1<AnyPos_T1) ) ]
		cat(paste("\nbuilding prel AnyPos_T1. Number of patients with !is.na(PosCD4_T1) & PosCD4_T1<AnyPos_T1, n=",df.all[tmp,][, length(unique(Patient))]))
		print( subset(df.all, PosCD4_T1<AnyPos_T1) )
		set(df.all, tmp, "AnyPos_T1", df.all[tmp, PosCD4_T1])
		#				
		#	reset Acute fields
		#
		df.all[, DUMMY:= df.all[, as.numeric(difftime(AnyPos_T1, NegT, units="days"))/365]]				
		tmp		<- df.all[, which(!is.na(DUMMY) & DUMMY<=1)]		
		set( df.all, tmp, 'Acute_Spec', 'NEGTEST_le12m' )
		set(df.all, df.all[, which(Acute_Spec%in%c('CLIN','LAB','NEGTEST_le12m'))], 'isAcute', 'Yes')				
		df.all[, DUMMY:=NULL]
		df.all[, idx:=NULL]
		#
		#	add participation in primo-SHM
		#
		df.pr	<- as.data.table(read.csv(file.primoSHM, stringsAsFactors=FALSE, header=FALSE, col.names='Patient'))
		df.pr[, PrimoSHM:='Yes']
		df.all	<- merge(df.all, df.pr, by='Patient', all.x=1)
		set(df.all, df.all[, which(is.na(PrimoSHM))], 'PrimoSHM', 'No' )
		set(df.all, NULL, 'PrimoSHM', df.all[, factor(PrimoSHM, levels=c('No','Yes'), labels=c('No','Yes'))])
		#
		#	get GGD at time of diagnosis
		#
		#load(file.ggd)		
		#tmp		<- subset(df, select=c(Patient, GGD_Reg))
		#setkey(tmp, Patient, GGD_Reg)
		#tmp		<- df[, {
		#			list(GGD_Reg=GGD_Reg, GGD_Last=c(GGD_Reg[-1]-1, as.Date('2015-05-01')))
		#		}, by='Patient']
		#df		<- merge(df, tmp, by=c('Patient','GGD_Reg'))		
		#tmp		<- subset(df.all, select=c(Patient, AnyPos_T1))
		#tmp		<- merge(unique(tmp), df, by='Patient', allow.cartesian=TRUE)
		#setkey(tmp, Patient, GGD_Reg)
		#tmp		<- tmp[,	{
		#			z		<- which(GGD_Reg<=AnyPos_T1 & AnyPos_T1<GGD_Last)
		#			z2		<- 'At_Diag'
		#			if(length(z)==0)
		#			{
		#				z	<- which.min(difftime(GGD_Reg, AnyPos_T1))
		#				z2	<- 'Earliest_After_Diag'
		#			}
		#			list(RegionGGD_AnyPosT1_Type=z2, GGD_RecordTime=GGD_Reg[z[1]], Region_AnyPosT1=Region[z[1]], GGD_AnyPosT1=GGD[z[1]])
		#		} , by='Patient']	
		#			
		#set(tmp, tmp[, which(Region_AnyPosT1!='Rott' & GGD_AnyPosT1%in%c('Groningen','Drenthe','Fryslan'))], 'Region_AnyPosT1', 'North')
		#set(tmp, tmp[, which(Region_AnyPosT1!='Rott' & GGD_AnyPosT1%in%c('IJsselland','Twente','Gelre_IJssel','Hulpverlening_Gelderland_Midden','Flevoland',"Rivierenland","Nijmegen","Midden_Nederland"))], 'Region_AnyPosT1', 'East')		
		#set(tmp, tmp[, which(Region_AnyPosT1!='Rott' & GGD_AnyPosT1%in%c('Amsterdam'))], 'Region_AnyPosT1', 'Amst')
		#set(tmp, tmp[, which(Region_AnyPosT1!='Rott' & GGD_AnyPosT1%in%c('Utrecht','Gooi_Vechtstreek','Zuid_Holland_Zuid','Zeeland','Rotterdam_Rijnmond','Den_Haag','Hollands_Midden','Kennemerland','Zaanstreek_Waterland','Hollands_Noorden',"Zuid_Holland_West"))], 'Region_AnyPosT1', 'West')
		#set(tmp, tmp[, which(Region_AnyPosT1!='Rott' & GGD_AnyPosT1%in%c('West_Brabant','Hart_voor_Brabant','Brabant_Zuidoost','Limburg-Noord','Zuid_Limburg'))], 'Region_AnyPosT1', 'South')
		#set(tmp, NULL, 'Region_AnyPosT1', tmp[, factor(as.character(Region_AnyPosT1))])		
		#cat('\nPatients with no entry in GGD file', setdiff( df.all[, Patient], tmp[, Patient] ))		
		#df.all	<- merge(df.all, tmp, by='Patient', all.x=1)
		#
		#	re-arrange cols
		#
		df.all	<- df.all[, c(	"Patient", "DateBorn", "Sex", "CountryBorn", "RegionOrigin", "DateLastContact", "DateDied", "isDead",      
						"NegT","NegT_Acc","NegT_Crude","ACS",
						"AnyPos_T1","CountryInfection","Trm","isAcute","Acute_Spec","Region_first","GGD_first","GGD_first_Date","Region_now","GGD_now",
						"DateInCare", "FirstMed","PoslRNA_T1", "lRNA_T1","PoslRNAg500_T1","PosCD4_T1","CD4_T1",
						"AnyT_T1","AnyT_T1_Acc","AnyT_T1_Crude","PrimoSHM","TrI.n","TrI.p","DateAIDS",					
						"PosSeqT","FASTASampleCode","SUBTYPE_C","SEQ_L"), with=FALSE]          
		if(verbose)	cat(paste("\nsave to file",file.out))
		#stopifnot( !nrow(subset(df.all, is.na(AnyPos_T1))) )
		save(df.all,file=file.out)
		#
		#	save with Dates in numeric format
		#
		for(x in colnames(df.all))
			if(class(df.all[[x]])=='Date')
				set(df.all, NULL, x, hivc.db.Date2numeric(df.all[[x]]))		
		save(df.all,file=file.out.numeric)
	}
}
