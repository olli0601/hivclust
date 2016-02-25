######################################################################################
project.ACdata.size.proposal.160223<- function()
{
	#	get seroconverters M/F; last neg test; first pos test; second available DBS; bounded ID		
	indir	<- '/Users/Oliver/duke/20160_AC/original_160212'
	load(file.path(indir,'RD02-01_ACDIS_Demography_Curated.rda'))
		
	#	IND with a neg test
	#dneg	<- unique(subset(adem, HIVNegative=='Y' | !is.na(LatestHIVNegative) | !is.na(EarliestHIVNegative), IIntId))
	dneg	<- unique(subset(adem, !is.na(LastHIVNegative), IIntId))
	#	IND with a pos test
	#dpos	<- unique(subset(adem, HIVPositive=='Y' | !is.na(EarliestHIVPositive), IIntId))
	dpos	<- unique(subset(adem, !is.na(FirstHIVPositive), IIntId))
	#	seroconverters
	dsn		<- merge(dneg, dpos, by='IIntId')
	#	seroconverters data
	dsnd	<- merge(dsn, adem, by='IIntId')
	setkey(dsnd, IIntId, Episode)
	tmp		<- unique(subset(dsnd, LastHIVNegative>FirstHIVPositive, IIntId))
	if(nrow(tmp))	
		cat('\nWarning: Seroconverters with EarliestHIVPositive<=LatestHIVNegative', tmp[, IIntId])	
	#	seroconverters data: add high/low incidence in bounded structure
	load(file.path(indir,'HS01-01_HotSpot_Dataset.rda'))
	dsnd	<- merge(dsnd, ahs, by='BSIntId',all.x=1)
	setnames(dsnd, 'Cluster_2','BSIncidence')
	tmp		<- subset(dsnd, !is.na(BSIntId) & is.na(BSIncidence))[, unique(BSIntId)]
	cat('\nWarning: BSIntId without Low/high incidence category', tmp)
	#
	#	seroconverters data: add visits with HIV+ result
	#
	load(file.path(indir,'RD05-99_ACDIS_HIV_All_Curated.rda'))
	ahiv	<- subset(ahiv, select=c(IIntId, VisitDate, HIVSpecimenType, HIVResult, KnowsHIVStatus))
	ahiv	<- subset(ahiv, HIVResult=='Pos')	
	stopifnot( !length(setdiff( ahiv[, IIntId], dpos[,IIntId]))	)	#check that there are no HIV+ in HIV Individuals that are not in Demography
	#	Note: Some people seem to forget their status
	df		<- subset(dsnd, select=c(IIntId, ObservationStart, ObservationEnd ))
	df		<- merge(df, df[, list(ObservationStart_Mn=min(ObservationStart), ObservationEnd_Mx= max(ObservationEnd) ), by='IIntId'], by='IIntId')	
	df		<- merge(df, ahiv, by='IIntId', allow.cartesian=1)
	df		<- subset(df, ObservationEnd_Mx<VisitDate | VisitDate<ObservationStart_Mn | (ObservationStart<=VisitDate & VisitDate<=ObservationEnd))	
	stopifnot( !nrow(subset(df, VisitDate<ObservationStart_Mn))	)
	tmp		<- subset(df, ObservationEnd_Mx<VisitDate & ObservationEnd==ObservationEnd_Mx)
	df		<- rbind( subset(df, ObservationStart<=VisitDate & VisitDate<=ObservationEnd), tmp )
	#	Note: df now contains all Episodes with any HIV+ result
	cat('\nFound Episodes with an HIV+ visit, n=', nrow(df))	
	#	Note: there are multiple HIV+ visits per episode; we are only interested in the latest one
	df		<- merge(df, df[, list(VisitDate=max(VisitDate)), by=c('IIntId','ObservationStart')], by=c('IIntId','ObservationStart','VisitDate'))
	df		<- subset(df, select=c(IIntId, ObservationStart, ObservationEnd, VisitDate, HIVSpecimenType))
	setnames(df, 'VisitDate', 'VisitDate_HIVPos')
	dsnd	<- merge(dsnd, df, by=c('IIntId','ObservationStart','ObservationEnd'), all.x=1)
	#
	# 	done: seroconverters data: add visits with HIV+ result
	#	
	#
	#	seroconverters data: add birth date
	#
	load(file.path(indir,'RD01-01_ACDIS_Individuals_OR.rda'))	
	dsnd	<- merge(dsnd, subset(aind, select=c(IIntId, DateOfBirth)), by='IIntId', all.x=1)
	stopifnot( dsnd[, !any(is.na(DateOfBirth))]	)
	#	seroconverters data: add age at midpoint of episode
	dsnd[, Age_Mid:= (ObservationStart+ObservationEnd)/2-DateOfBirth]
	dsnd[, AgeGrp:=NULL]
	#	now: get seroconverter positive visits
	#		 start with first visits
	dsn		<- unique(subset(dsnd, select=c(IIntId, FirstHIVPositive)))
	setnames(dsn, 'FirstHIVPositive', 'VisitDate_HIVPos')
	dsn[, Visit_ID:= 1L]
	#		add HIV+ follow up visits to dsn
	df		<- subset(dsnd, !is.na(VisitDate_HIVPos) & FirstHIVPositive<VisitDate_HIVPos, c(IIntId, VisitDate_HIVPos))	
	df		<- merge(df, df[, list(VisitDate_HIVPos=sort(VisitDate_HIVPos), Visit_N=1+length(VisitDate_HIVPos), Visit_ID= 1+seq_len(length(VisitDate_HIVPos))), by='IIntId'], by=c('IIntId','VisitDate_HIVPos'))
	dsn		<- rbind(dsn, df, use.names=TRUE, fill=TRUE)
	set(dsn, dsn[, which(is.na(Visit_N))], 'Visit_N', 1L)
	#		add info for each seroconverter
	dsn		<- merge(dsn, dsnd[, list(EverResident=all(EpisodeType!='resident'), LastHIVNegative=LastHIVNegative[1], FirstHIVPositive=FirstHIVPositive[1], Sex=Sex[1]), by='IIntId'], by='IIntId') 
	set(dsn, NULL, 'EverResident', dsn[, factor(EverResident, levels=c(FALSE,TRUE), labels=c('Y','N'))])
	dsn[, SC_mid:= (LastHIVNegative+FirstHIVPositive)/2]
	dsn[, SC_mid_C:= floor(SC_mid)]
	#		add info for each visit
	tmp		<- subset(dsnd, !is.na(VisitDate_HIVPos), select=c(IIntId, VisitDate_HIVPos, BSIntId, BSIncidence, Area, Age_Mid))
	dsn		<- merge(dsn, tmp, by=c('IIntId', 'VisitDate_HIVPos'), all.x=1)
	#		add info for each visit: in some cases EarliestPositive is outside episode
	tmp		<- unique(subset(dsn, is.na(Age_Mid), IIntId))
	tmp		<- merge(dsnd, tmp, by='IIntId')
	tmp		<- tmp[, {
				z<- which.min( abs((ObservationStart+ObservationEnd)/2-FirstHIVPositive) )
				list(VisitDate_HIVPos= FirstHIVPositive[1], BSIntId=BSIntId[z], BSIncidence=BSIncidence[z], Area=Area[z], Age_Mid=Age_Mid[z])
			}, by='IIntId']
	for(i in seq_len(nrow(tmp)))
	{
		tmp2	<- dsn[, which(IIntId==tmp$IIntId[i] & VisitDate_HIVPos==tmp$VisitDate_HIVPos[i])]
		stopifnot( length(tmp2)==1 )
		set(dsn, tmp2, 'BSIntId', tmp$BSIntId[i])
		set(dsn, tmp2, 'BSIncidence', tmp$BSIncidence[i])
		set(dsn, tmp2, 'Area', tmp$Area[i])
		set(dsn, tmp2, 'Age_Mid', tmp$Age_Mid[i])						       
	}	
	stopifnot( !nrow(subset(dsn, is.na(Age_Mid)))	)
	#
	#	
	#
	df		<- as.data.table(expand.grid(select.mininf=seq(0, 5, 0.5), select.max.scw=c(2,3,4,5,6,7,8,Inf)))	
	df		<- df[, {
						z		<- subset(dsn, VisitDate_HIVPos-SC_mid>=select.mininf & FirstHIVPositive-LastHIVNegative<select.max.scw )
						z		<- merge(z, z[, list(VisitDate_HIVPos=min(VisitDate_HIVPos)), by='IIntId'], by=c('IIntId','VisitDate_HIVPos'))
						area	<- z[, table(Area)]
						sex		<- z[, table(Sex)]
						hl		<- z[, table(BSIncidence)]
						c(list(SC_n= nrow(z)), as.list(sex), as.list(area), as.list(hl))
			}, by=c('select.max.scw','select.mininf')]
	
	subset(df, select.max.scw%in%c(5,8) & select.mininf>=2)
	subset(df, select.max.scw%in%c(5,8) & select.mininf>=2)*0.1
	
	dsn[, length(unique(IIntId))]
	
	subset(dsn, VisitDate_HIVPos-SC_mid<0 )

	setkey(dsn, SC_mid_C, Sex, Area)
	dsn[, list(SC_N=length(IIntId)), by=c('SC_mid_C','Sex','Area')]

	
	#	for most seroconverters, the seroconversion time falls within an observation period
	#	get age, BSId, LastNeg test, Area IF it did not change during the seroconversion window
	setkey(dsnd, IIntId, Episode)
	tmp		<- subset(dsnd, ObservationStart<=EarliestHIVPositive & EarliestHIVPositive<=ObservationEnd, c(IIntId, LatestHIVNegative, EarliestHIVPositive, BSIntId, BSIncidence, Area))	
	tmp2	<- tmp[, {
				#IIntId<- 236; LatestHIVNegative<- 2005.288; EarliestHIVPositive<- 2005.658
				z<- which(dsnd$IIntId==IIntId & LatestHIVNegative<=dsnd$ObservationEnd & dsnd$ObservationStart<EarliestHIVPositive)
				z<- dsnd[z,]
				list(Area_N=z[,length(unique(Area))], BS_N=z[,length(unique(BSIntId))], BSIncidence_N=z[,length(unique(BSIncidence))])				
			}, by='IIntId']
	tmp		<- merge(tmp, tmp2, by='IIntId')	
	dsn		<- merge(dsn, subset(tmp, select=c(IIntId, BSIntId, BSIncidence, Area, BS_N, BSIncidence_N, Area_N)), by='IIntId', all.x=1)
	set(dsn, dsn[, which(Area_N>1)],'Area', NA_character_)
	set(dsn, dsn[, which(BS_N>1)],'BSIntId', NA_integer_)
	set(dsn, dsn[, which(BSIncidence_N>1)],'BSIncidence', NA_character_)
	dsn[, SC_mid:= (LatestHIVNegative+EarliestHIVPositive)/2]
	#
	#	add HIV+ follow up visits to dsn
	#
	df		<- subset(dsnd, !is.na(VisitDate_HIVPos) & EarliestHIVPositive<VisitDate_HIVPos, c(IIntId, VisitDate_HIVPos, EarliestHIVPositive))	
	df		<- merge(df, df[, list(VisitDate_HIVPos=sort(VisitDate_HIVPos), FollowUpVisit_N=length(VisitDate_HIVPos), FollowUpVisit_ID= seq_len(length(VisitDate_HIVPos))), by='IIntId'], by=c('IIntId','VisitDate_HIVPos'))
	df		<- merge(df, subset(dsn, select=c(IIntId, SC_mid)), by='IIntId')
	
	
	df[, ]
	
	
	
}
######################################################################################
project.ACdata.rda.basic.HIV_Individuals.160224<- function()
{
	#	HIV INDIVIDUALS
	indir	<- '/Users/Oliver/duke/20160_AC/original_160212'
	load(file.path(indir,'RD05-99_ACDIS_HIV_All.rda'))
	#	dates into numerical format
	set(ahiv, NULL, 'VisitDate', ahiv[, gsub(' 00:00:00.000','',VisitDate)])
	set(ahiv, NULL, 'VisitDate', ahiv[, hivc.db.Date2numeric(VisitDate)])
	#	Sex
	set(ahiv, ahiv[, which(Sex==9)],'Sex',NA_integer_)
	set(ahiv, NULL, 'Sex', ahiv[, factor(Sex, levels=c(1,2),labels=c('M','F'))])
	#	HIV refused
	set(ahiv, ahiv[, which(HIVRefused%in%c(98,99))], 'HIVRefused', NA_integer_)
	#	TODO: what is default? Not applicable?
	set(ahiv, NULL, 'HIVRefused', ahiv[, factor(HIVRefused, levels=c(1,2,95,97),labels=c('Y','N','Default','Not_applicable'))])
	#	HIV refused by
	#	TODO: what is Not applicable?
	set(ahiv, ahiv[, which(HIVRefusedBy==99)], 'HIVRefusedBy', NA_integer_)
	set(ahiv, NULL, 'HIVRefusedBy', ahiv[, factor(HIVRefusedBy, levels=c(1,2,3,4,5,97),labels=c('Self','Partner','HH_head','BS_owner','Other','Not_applicable'))])
	#	HIV sample available
	set(ahiv, NULL, 'HIVSampleAvailable', ahiv[, factor(HIVSampleAvailable, levels=c(1,2),labels=c('Y','N'))])
	#	HIV specimen type
	set(ahiv, ahiv[, which(HIVSpecimenType==97)], 'HIVSpecimenType', NA_integer_)
	set(ahiv, NULL, 'HIVSpecimenType', ahiv[, factor(HIVSpecimenType, levels=c(1,2,3,4),labels=c('DBS','Cap-Ind','Cap-Pool','Oral'))])
	#	HIV result
	set(ahiv, ahiv[, which(HIVResult%in%c(97,99))], 'HIVResult', NA_integer_)
	set(ahiv, NULL, 'HIVResult', ahiv[, factor(HIVResult, levels=c(0,1,2,3,4),labels=c('Neg','Pos','Indeterminate','Insufficient Volume','Spurious'))])
	#	KnowsHIVStatus
	#	TODO what is not applicable?
	set(ahiv, ahiv[, which(KnowsHIVStatus%in%c(99))], 'KnowsHIVStatus', NA_integer_)
	set(ahiv, NULL, 'KnowsHIVStatus', ahiv[, factor(KnowsHIVStatus, levels=c(1,2,93,97),labels=c('Y','N','Refused','Not_applicable'))])
	#	save curated
	save(ahiv, file=file.path(indir,'RD05-99_ACDIS_HIV_All_OR.rda'))
}
######################################################################################
project.ACdata.rda.basic.Individuals.160225<- function()
{
	indir	<- '/Users/Oliver/duke/20160_AC/original_160212'
	load(file.path(indir,'RD01-01_ACDIS_Individuals.rda'))
	#	Sex
	set(aind, aind[, which(Sex==9)],'Sex',NA_integer_)
	set(aind, NULL, 'Sex', aind[, factor(Sex, levels=c(1,2),labels=c('M','F'))])
	#	DateOfBirth
	set(aind, NULL, 'DateOfBirth', aind[, hivc.db.Date2numeric(DateOfBirth)])
	#	IEndDate
	set(aind, NULL, 'IEndDate', aind[, hivc.db.Date2numeric(IEndDate)])
	#	IEndType
	set(aind, NULL, 'IEndType', aind[, factor(IEndType, levels=c(10, 52, 60, 61, 64, 90, 102, 110, 111, 114),labels=c('Visit','HH_dissolution','HH_migration','HH_migration_internal','HH_migration_out','Death','HH_member_end','IND_migration','IND_migration_internal','IND_migration_out'))])
	#	save curated
	save(ahiv, file=file.path(indir,'RD01-01_ACDIS_Individuals_OR.rda'))	
}
######################################################################################
project.ACdata.rda.curate<- function()
{
	#
	#	HIV INDIVIDUAL
	#	curate latest negative
	#
	load(file.path(indir,'RD05-99_ACDIS_HIV_All_OR.rda'))
	setkey(ahiv, IIntId, VisitDate)
	tmp		<- subset(ahiv, HIVResult=='Neg', select=c(IIntId, VisitDate, HIVSampleAvailable, HIVSpecimenType, HIVResult))
	df		<- tmp[, list(LatestHIVNegative=max(VisitDate)), by='IIntId']
	tmp		<- subset(ahiv, HIVResult=='Pos')[, list(EarliestHIVPositive= min(VisitDate)), by='IIntId']
	df		<- merge(df, tmp, by='IIntId', all=1)
	tmp		<- subset(df, LatestHIVNegative>EarliestHIVPositive)
	cat('\nFound HIV- after HIV+, for IIntId', tmp[, unique(IIntId)])
	#	get evidence for HIV+ or HIV-
	df		<- merge(ahiv, unique(subset(tmp, select=IIntId)), by='IIntId')
	df		<- merge(df, tmp, by='IIntId')
	tmp2	<- subset(df, HIVResult=='Neg' & VisitDate>=EarliestHIVPositive)[, list(Neg_N=length(HIVResult), Neg_Cap_Any= any(HIVSpecimenType=='Cap-Ind')), by='IIntId']
	tmp		<- merge(tmp, tmp2, by='IIntId')	
	tmp2	<- subset(df, HIVResult=='Pos' & VisitDate<=LatestHIVNegative)[, list(Pos_N=length(HIVResult), Pos_Cap_Any= any(HIVSpecimenType=='Cap-Ind')), by='IIntId']
	tmp		<- merge(tmp, tmp2, by='IIntId')
	tmp[, HIVResult_Curated:= 'Pos']
	#	keep HIV- if there is a capillary neg test after the pos test and the pos test was NA or DBS
	stopifnot( !nrow(subset(tmp, Neg_Cap_Any & Pos_Cap_Any)) )
	set(tmp, tmp[, which(Neg_Cap_Any)], 'HIVResult_Curated', 'Neg')
	#	keep HIV- if there is more neg tests than pos tests
	set(tmp, tmp[, which(Neg_N>Pos_N)], 'HIVResult_Curated', 'Neg')
	#	set Indeterminate if no capillary HIV+ and number pos==number neg 
	set(tmp, tmp[, which(Neg_N==Pos_N & (is.na(Pos_Cap_Any) | !Pos_Cap_Any))], 'HIVResult_Curated', 'Indeterminate')
	df		<- copy(tmp)
	#	Note: the following only works because there is at most only one visit per individual to reset	
	stopifnot( subset(df, HIVResult_Curated=='Neg')[, all(Pos_N==1)] )
	stopifnot( subset(df, HIVResult_Curated=='Indeterminate')[, all(Pos_N==1|Neg_N==1)] )
	#	curate to neg
	tmp		<- subset(df, HIVResult_Curated=='Neg')
	setnames(tmp, 'EarliestHIVPositive', 'VisitDate')
	ahiv	<- merge(ahiv, subset(tmp, select=c(IIntId, VisitDate, HIVResult_Curated)), by=c('IIntId','VisitDate'), all.x=1)
	tmp		<- ahiv[, which(!is.na(HIVResult_Curated))]
	set(ahiv, tmp, 'HIVResult', ahiv[tmp, HIVResult_Curated])
	set(ahiv, NULL, 'HIVResult_Curated', NULL)
	#	curate to indeterminate
	tmp		<- subset(df, HIVResult_Curated=='Indeterminate')
	tmp		<- data.table:::melt.data.table(tmp, measure.vars=c('LatestHIVNegative','EarliestHIVPositive'), value.name='VisitDate')
	ahiv	<- merge(ahiv, subset(tmp, select=c(IIntId, VisitDate, HIVResult_Curated)), by=c('IIntId','VisitDate'), all.x=1)
	tmp2	<- ahiv[, which(!is.na(HIVResult_Curated))]
	set(ahiv, tmp2, 'HIVResult', ahiv[tmp2, HIVResult_Curated])
	set(ahiv, NULL, 'HIVResult_Curated', NULL)	
	#	curate to pos
	df		<- subset(df, HIVResult_Curated=='Pos')
	tmp		<- subset(ahiv, select=c(IIntId, VisitDate, HIVResult))
	tmp		<- merge(tmp, df, by='IIntId')
	tmp		<- subset(tmp, VisitDate<=LatestHIVNegative & EarliestHIVPositive<=VisitDate & HIVResult=='Neg')
	ahiv	<- merge(ahiv, subset(tmp, select=c(IIntId, VisitDate, HIVResult_Curated)), by=c('IIntId','VisitDate'), all.x=1)
	tmp		<- ahiv[, which(!is.na(HIVResult_Curated))]
	set(ahiv, tmp, 'HIVResult', ahiv[tmp, HIVResult_Curated])
	set(ahiv, NULL, 'HIVResult_Curated', NULL)
	#
	setkey(ahiv, IIntId, VisitDate)
	save(ahiv, file=file.path(indir,'RD05-99_ACDIS_HIV_All_Curated.rda'))	
	#	HIV INDIVIDUAL
	#	end: curate latest negative
	#
	#
	#	HIV DEMOGRAPHY
	#	reset EarliestHIVPositive, LatestHIVNegative
	#	add LastHIVNegative FirstHIVPositive
	tmp		<- subset(ahiv, HIVResult=='Neg', select=c(IIntId, VisitDate, HIVSampleAvailable, HIVSpecimenType, HIVResult))
	df		<- tmp[, list(LatestHIVNegative=max(VisitDate)), by='IIntId']
	tmp		<- subset(ahiv, HIVResult=='Pos')[, list(EarliestHIVPositive= min(VisitDate)), by='IIntId']
	df		<- merge(df, tmp, by='IIntId', all=1)	
	stopifnot( !nrow(subset(df, LatestHIVNegative>EarliestHIVPositive))	)	
	load(file.path(indir,'RD02-01_ACDIS_Demography_OR.rda'))
	tmp		<- setdiff( subset(adem, !is.na(LatestHIVNegative))[, unique(IIntId)], subset(df, !is.na(LatestHIVNegative))[, unique(IIntId)]	)	
	setnames(df, c('LatestHIVNegative','EarliestHIVPositive'), c('LatestHIVNegative_Cur','EarliestHIVPositive_Cur'))
	adem	<- merge(adem, df, by='IIntId', all.x=1)
	tmp		<- adem[, which(!is.na(EarliestHIVPositive_Cur))]
	set(adem, tmp, 'EarliestHIVPositive', adem[tmp, EarliestHIVPositive_Cur])
	tmp		<- adem[, which(!is.na(EarliestHIVPositive_Cur & LatestHIVNegative>=EarliestHIVPositive_Cur))]
	set(adem, tmp, 'LatestHIVNegative', adem[tmp, EarliestHIVPositive_Cur]-1/365)	
	tmp		<- adem[, which(!is.na(LatestHIVNegative_Cur))]
	set(adem, tmp, 'LatestHIVNegative', adem[tmp, LatestHIVNegative_Cur])
	tmp		<- adem[, which(!is.na(LatestHIVNegative_Cur & EarliestHIVPositive<=LatestHIVNegative_Cur))]
	set(adem, tmp, 'EarliestHIVPositive', adem[tmp, LatestHIVNegative_Cur]+1/365)
	setnames(adem, c('LatestHIVNegative_Cur','EarliestHIVPositive_Cur'), c('LastHIVNegative','FirstHIVPositive'))
	#
	tmp		<- adem[, which(LatestHIVNegative<EarliestHIVNegative)]
	set(adem, tmp, 'EarliestHIVNegative', adem[tmp, LatestHIVNegative])
	#
	tmp		<- subset(adem, !is.na(LatestHIVNegative) & is.na(LastHIVNegative))[, unique(IIntId)]
	cat('\nWarning: Individuals with no Neg test but LatestHIVNegative in Demography', tmp)
	save(adem, file=file.path(indir,'RD02-01_ACDIS_Demography_Curated.rda'))
	#	HIV DEMOGRAPHY
	#	end: reset EarliestHIVPositive, LatestHIVNegative
}
######################################################################################
project.ACdata.rda.basic.Demography.160223<- function()
{
	#	DEMOGRAPHY
	indir	<- '/Users/Oliver/duke/20160_AC/original_160212'
	load(file.path(indir,'RD02-01_ACDIS_Demography.rda'))
	#	dates into numerical format
	set(adem, adem[, which(EarliestHIVPositive=='')], 'EarliestHIVPositive', NA_character_)	
	set(adem, NULL, 'EarliestHIVPositive', adem[, hivc.db.Date2numeric(EarliestHIVPositive)])
	set(adem, adem[, which(EarliestHIVNegative=='')], 'EarliestHIVNegative', NA_character_)
	set(adem, NULL, 'EarliestHIVNegative', adem[, hivc.db.Date2numeric(EarliestHIVNegative)])
	set(adem, adem[, which(LatestHIVNegative=='')], 'LatestHIVNegative', NA_character_)
	set(adem, NULL, 'LatestHIVNegative', adem[, hivc.db.Date2numeric(LatestHIVNegative)])
	set(adem, adem[, which(ObservationStart=='')], 'ObservationStart', NA_character_)
	set(adem, NULL, 'ObservationStart', adem[, hivc.db.Date2numeric(ObservationStart)])
	set(adem, NULL, 'ObservationEnd', adem[, hivc.db.Date2numeric(ObservationEnd)])
	#	Sex
	set(adem, adem[, which(Sex==9)],'Sex',NA_integer_)
	set(adem, NULL, 'Sex', adem[, factor(Sex, levels=c(1,2),labels=c('M','F'))])
	#	Residency
	set(adem, NULL, 'EpisodeType', adem[, factor(EpisodeType, levels=c(1,2,3),labels=c('resident','not_resident','not_resident_before_resident'))])
	#	Semester
	set(adem, NULL, 'ExpSem', adem[, factor(ExpSem, levels=c(1,2),labels=c('Jan-Jun','Jul-Dec'))])
	#	HIV+ and HIV- in episode
	set(adem, NULL, 'HIVPositive', adem[, factor(HIVPositive, levels=c(0,1),labels=c('N','Y'))])
	set(adem, NULL, 'HIVNegative', adem[, factor(HIVNegative, levels=c(0,1),labels=c('N','Y'))])
	#	Area
	set(adem, NULL, 'Area', adem[, factor(Area, levels=c(0,1,2,3),labels=c('Rural','PeriUrban','Urban','Outside'))])
	#	save curated
	save(adem, file=file.path(indir,'RD02-01_ACDIS_Demography_OR.rda'))	
}
######################################################################################
project.ACdata.csv.to.rda.160223<- function()
{
	indir	<- '/Users/Oliver/duke/20160_AC/original_160212'
	#
	#	save data sets as rda
	#	
	#	DEMOGRAPHY
	#	RD02-01 ACDIS Demography.csv
	#	De-anonymize IIntID: 	RD02-001 Demography Anonymised Individuals.csv
	#	De-anonymize BSIntID:	RD02-002 Demography Anonymised BoundedStructures.csv 	
	#	TODO: 		RD02-02 ACDIS DemographyYear.dta has no csv
	#	[	RecID Sex Episodes Episode EpisodeType ExpYear ExpSem AgeGrp ExpDays PeriodStart 
	#		AgedIn MemberStart InMigr InMigrEx Died MemberEnd OutMigr OutMigrEx Lost
	#		PeriodEnd ObservationStart ObservationEnd PositiveExp NegativeExp UnknownExp HIVPositive 
	#		HIVNegative EarliestHIVPositive LatestHIVNegative EarliestHIVNegative
	#		Area C_Unknown C_CMPN C_AIDS_TB C_NonComm C_Injuries LBCnt Parity DeliveryDate 
	#		IIntId BSIntId		]
	file	<- file.path(indir, 'RD02-01 ACDIS Demography.csv')
	adem	<- as.data.table(read.delim(file, stringsAsFactors=FALSE))	
	setnames(adem,'IIntID','id_anonymised')
	tmp		<- as.data.table(read.csv(file.path(indir,'RD02-001 Demography Anonymised Individuals.csv'))) 	#	merge with anonymised code
	adem	<- merge(adem, tmp, by='id_anonymised', all.x=1)
	stopifnot( adem[, !any(is.na(IIntId))] )	#check for error in anonymisation
	set(adem, NULL, 'id_anonymised', NULL)
	tmp		<- as.data.table(read.csv(file.path(indir,'RD02-002 Demography Anonymised BoundedStructures.csv'))) 	#	merge with anonymised code
	setnames(adem,'BSIntID','id_anonymised')
	tmp2	<- nrow(subset(adem, is.na(id_anonymised)))	#1085544	
	adem	<- merge(adem, tmp, by='id_anonymised', all.x=1)
	stopifnot( nrow(subset(adem, is.na(BSIntId)))==tmp2 )	
	set(adem, NULL, 'id_anonymised', NULL)	
	setkey(adem, IIntId)
	save( adem, file=gsub('csv','rda',gsub(' ','_',file)))
	
	
	#	CORE INDIVIDUALS
	#	RD01-01 ACDIS Individuals.csv
	#	De-anonymize IIntID: RD01-001 Core Anonymised Individuals.csv
	#	[Sex DateOfBirth IEndDate IEndType MotherIntID MotherPlaceAtBirth FatherIntID IIntId]
	file	<- file.path(indir,'RD01-01 ACDIS Individuals.csv')
	aind	<- as.data.table(read.csv(file, stringsAsFactors=FALSE))
	tmp		<- as.data.table(read.csv(file.path(indir,'RD01-001 Core Anonymised Individuals.csv')))	
	setnames(aind,'IIntID','id_anonymised')
	aind	<- merge(aind, tmp, by='id_anonymised', all.x=1)
	stopifnot( aind[, !any(is.na(IIntId))] )
	set(aind, NULL, 'id_anonymised', NULL)
	setkey(aind, IIntId)
	save( aind, file=gsub('csv','rda',gsub(' ','_',file)))
	
	
	#	BOUNDED STRUCTURES
	#	RD01-03 ACDIS BoundedStructures.csv
	#	De-anonymize BSIntId: RD01-003 Core Anonymised BoundedStructures.csv
	#	TODO: the name is BSIntId. Is this anonymized? (Usually the anonymized names are ending in ID)
	#	[ 	BSStartDate BSStartType  BSEndDate BSEndType Isigodi LocalArea IsUrbanOrRural NearestClinic 
	#		KmToNearestClinic NearestSecondarySchool KmToNearestSecondarySchool NearestPrimarySchool KmToNearestPrimarySchool 
	#		KmToNearestLevel1Road KmToNearestLevel2Road BSIntId	]
	file	<- file.path(indir,'RD01-03 ACDIS BoundedStructures.csv')
	abs		<- as.data.table(read.csv(file, stringsAsFactors=FALSE))
	tmp		<- as.data.table(read.csv(file.path(indir,'RD01-003 Core Anonymised BoundedStructures.csv')))	
	setnames(abs,'BSIntId','id_anonymised')
	abs	<- merge(abs, tmp, by='id_anonymised', all.x=1)
	stopifnot( abs[, !any(is.na(BSIntId))] )
	set(abs, NULL, 'id_anonymised', NULL)
	setkey(abs, BSIntId)	
	save( abs, file=gsub('csv','rda',gsub(' ','_',file)))
	
	
	#	HOTSPOTS
	#	HS01-01 HotSpot Dataset.csv
	#	Not anonymized	
	file	<- file.path(indir, 'HS01-01 HotSpot Dataset.csv')
	ahs		<- as.data.table(read.delim(file, stringsAsFactors=FALSE))	
	set(ahs, NULL, 'Cluster_2', ahs[, factor(Cluster_2,levels=c(0,1),labels=c('LowInc','HighInc'))])
	save( ahs, file=gsub('csv','rda',gsub(' ','_',file)))
	
	
	
	#	HIV ALL
	#	De-anonymize IIntId: RD05-001 HIV Anonymised Individuals.csv
	#	De-anonymize ResidencyBSIntId:	RD05-003 HIV Anonymised BoundedStructures.csv
	file	<- file.path(indir,'RD05-99 ACDIS HIV All.csv')
	ahiv	<- as.data.table(read.csv(file, stringsAsFactors=FALSE))
	tmp		<- as.data.table(read.csv(file.path(indir,'RD05-001 HIV Anonymised Individuals.csv')))	
	setnames(ahiv,'IIntId','id_anonymised')
	ahiv	<- merge(ahiv, tmp, by='id_anonymised', all.x=1)
	stopifnot( ahiv[, !any(is.na(IIntId))] )
	set(ahiv, NULL, 'id_anonymised', NULL)	
	tmp		<- as.data.table(read.csv(file.path(indir,'RD05-003 HIV Anonymised BoundedStructures.csv')))
	tmp2	<- nrow(subset(ahiv, is.na(ResidencyBSIntId)))	#1085544
	setnames(ahiv,'ResidencyBSIntId','id_anonymised')
	ahiv	<- merge(ahiv, tmp, by='id_anonymised', all.x=1)
	stopifnot( nrow(subset(ahiv, is.na(BSIntId)))==tmp2 )
	set(ahiv, NULL, 'id_anonymised', NULL)
	save( ahiv, file=gsub('csv','rda',gsub(' ','_',file)))
	
	
	
	#	RELATIONSHIPS
	#	RD01-09 ACDIS Conjugal Relationships.csv
	#	De-anonymize CRId: RD01-012 Core Anonymised ConjugalRelationships.csv
	#	[	CRId FemalePartnerId MalePartnerId CRStartDate CRStartType  CREndDate CREndType MarriageDate 
	#		MarriageType CRFirstRecorded MarriageEventFirstRecorded]
	file	<- file.path(indir,'RD01-09 ACDIS Conjugal Relationships.csv')
	arls	<- as.data.table(read.csv(file, stringsAsFactors=FALSE))
	tmp		<- as.data.table(read.csv(file.path(indir,'RD01-012 Core Anonymised ConjugalRelationships.csv')))	
	setnames(arls,'CRId','id_anonymised')
	arls	<- merge(arls, tmp, by='id_anonymised', all.x=1)
	stopifnot( arls[, !any(is.na(CRId))] )
	set(arls, NULL, 'id_anonymised', NULL)
	save( arls, file=gsub('csv','rda',gsub(' ','_',file)))
	
	
	#
	#	LAB RESULTS
	#	RD10-03 Lab Results.csv
	#	Confirmed: ACDIS_IIntID is not anonymized
	file	<- file.path(indir, 'RD10-03 Lab Results.csv')
	alb		<- as.data.table(read.delim(file, stringsAsFactors=FALSE))	
	save( alb, file=gsub('csv','rda',gsub(' ','_',file)))
}
######################################################################################
project.AC.cleanalignment<- function()
{
	sx				<- read.dna("~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data/AC_Geneious1012Seq.gag_nef.fasta", format='fasta')
	sxd				<- data.table(TAXA= rownames(sx))
	set(sxd, NULL, 'TAXA', sxd[, gsub('_$','',TAXA)])
	set(sxd, NULL, 'TAXA', sxd[, gsub(' ','R',TAXA)])
	tmp				<- sxd[, which(duplicated(TAXA))]
	set(sxd, tmp, 'TAXA', sxd[tmp, paste(TAXA,'-R',sep='')])
	rownames(sx)	<- sxd[, TAXA]
	
	
	tmp	<- sxd[, which(grepl('AF411967',TAXA))]
	sx	<- sx[ -tmp[-1], ]
	write.dna(sx, file="~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data/AC_Geneious1012Seq.gag_nef.OR.fasta", format='fasta')
	seq.write.dna.phylip(sx, file="~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data/AC_Geneious1012Seq.gag_nef.OR.phylip")
}
######################################################################################
project.ACpolext.rmDRM.150913<- function()
{
	indir	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_SA/ACpolext150831'
	infile	<- 'SATURN150831.csv'
	dfstr	<- data.table(read.csv(file=paste(indir, '/', infile, sep=''), stringsAsFactors=FALSE))
	infile	<- 'Eduan_DRT_170815.csv'
	dfdrt	<- data.table(read.csv(file=paste(indir, '/', infile, sep=''), stringsAsFactors=FALSE))
	infile	<- 'Comet150831.csv'
	dfcmt	<- data.table(read.csv(file=paste(indir, '/', infile, sep=''), stringsAsFactors=FALSE))	
	#	what subtypes do we have? 
	dfcmt[, table(subtype)]
	dfdrt[, table(assignment)]
	#	OK separate all C and the non-C sequences.  
	dfc		<- subset(dfdrt, assignment=='HIV-1 Subtype C', name) 
	dfnc	<- subset(dfdrt, assignment!='HIV-1 Subtype C', name)
	#	clarify: not all C? do we have subtype assignment for all others that are not in SATURN?
	
	
	indir				<- '/Users/Oliver/Dropbox\ (Infectious Disease)/2015_AC_Origin_Cascade'
	infile				<- 'ZA.SubC.12294.13.09.15.aln.fasta'
	acxs				<- read.dna(file=paste(indir, '/', infile,sep=''), format='fasta')
	tmp					<- which(grepl("B.FR.K03455.1983",rownames(acxs)))
	rownames(acxs)[tmp]	<- 'HXB2'
	tmp					<- big.phylo:::seq.rm.drugresistance(acxs)
	nodr.info			<- tmp$nodr.info
	seq					<- tmp$nodr.seq
	write.dna(seq, file= paste(indir, '/', gsub('13.09.15.aln','alnnoDRM_150913',infile), sep=''), format='fasta', colsep='', nbcol=-1)
	save(seq, nodr.info, file= paste(indir, '/', gsub('13.09.15.aln\\.fasta','alnnoDRM_150913\\.R',infile), sep=''))
}
######################################################################################
project.ACpolext.rmDRM.160120<- function()
{
	require(big.phylo)
	infile				<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_SA/ACpolext160120/ZA_aln4_160120.fasta'
	acxs				<- read.dna(file=infile, format='fasta')
	tmp					<- which(grepl("B.FR.K03455.1983",rownames(acxs)))
	rownames(acxs)[tmp]	<- 'HXB2'
	tmp					<- big.phylo:::seq.rm.drugresistance(acxs)
	nodr.info			<- tmp$nodr.info
	seq					<- tmp$nodr.seq
	write.dna(seq, file= gsub('aln4','aln4noDRM',infile), format='fasta')
	save(seq, nodr.info, file= gsub('aln4','aln4noDRM',infile))
}
######################################################################################
project.ACpolext.rmDRM.160209<- function()
{
	require(big.phylo)
	infile				<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_SA/ACpolext160209/ViroIntro_16442_20160208_aln.fasta'
	acxs				<- read.dna(file=infile, format='fasta')
	tmp					<- which(grepl("B.FR.K03455.1983",rownames(acxs)))
	rownames(acxs)[tmp]	<- 'HXB2'
	tmp					<- big.phylo:::seq.rm.drugresistance(acxs)
	nodr.info			<- tmp$nodr.info
	seq					<- tmp$nodr.seq
	write.dna(seq, file= gsub('ViroIntro_','ViroIntro_noDRM_',infile), format='fasta')
	save(seq, nodr.info, file= gsub('fasta','rda',gsub('ViroIntro_','ViroIntro_noDRM_',infile)))
}
######################################################################################
project.ACpolext.rmDRM.150907<- function()
{
	indir	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_SA/ACpolext150831'
	infile	<- 'SATURN150831.csv'
	dfstr	<- data.table(read.csv(file=paste(indir, '/', infile, sep=''), stringsAsFactors=FALSE))
	infile	<- 'Eduan_DRT_170815.csv'
	dfdrt	<- data.table(read.csv(file=paste(indir, '/', infile, sep=''), stringsAsFactors=FALSE))
	infile	<- 'Comet150831.csv'
	dfcmt	<- data.table(read.csv(file=paste(indir, '/', infile, sep=''), stringsAsFactors=FALSE))	
	#	what subtypes do we have? 
	dfcmt[, table(subtype)]
	dfdrt[, table(assignment)]
	#	OK separate all C and the non-C sequences.  
	dfc		<- subset(dfdrt, assignment=='HIV-1 Subtype C', name) 
	dfnc	<- subset(dfdrt, assignment!='HIV-1 Subtype C', name)
	#	clarify: not all C? do we have subtype assignment for all others that are not in SATURN?
	seq				<- seq[c("B.FR.83.HXB2_LAI_IIIB_BRU.K03455",arvdat[, PNG_ID_FULL]),]
	rownames(seq)[1]	<- 'HXB2'
	outfile			<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_PotentialDRMs_OR_150910.R'
	tmp				<- big.phylo:::seq.rm.drugresistance(seq, outfile=outfile)
	nodr.info		<- tmp$nodr.info
	nodr.seq		<- tmp$nodr.seq
	
	#	determine alignment start relative to HXB2
	#	pos 60 is pos 2300 in HXB2 ie pos 1 is pos 2241 in HXB2
	paste(as.character(acxs[ which(grepl("B.FR.K03455.1983",rownames(acxs))), ]), collapse='')
	acxs.st	<- 2241
	#	rm primary DRMs
	load( system.file(package="hivclust", "data",'IAS_primarydrugresistance_201303.rda') )		
	dr				<- as.data.table(IAS_primarydrugresistance_201303)				
	set(dr, NULL, "Alignment.nuc.pos", dr[,Alignment.nuc.pos]-acxs.st+1)		
	seq				<- seq.rm.drugresistance(as.character(acxs), dr, verbose=1, rtn.DNAbin=1 )
	#	removed DR mutations, n= 19340 n per patient= 1.55566280566281
	#	rm any additional DMRs that are indicated in the Excel sheet
	dfdrm			<- subset(dfstr, select=c(Seq_Code, PI_MAJOR, PI_MINOR, NRTI, NNRTI))
	dfdrm			<- subset(melt(dfdrm, id.vars=c('Seq_Code'), variable.name='TYPE', value.name='DRM'), !DRM%in%c('None','No PR seq','No RT seq'))
	dfdrm			<- dfdrm[, list(DRM=unlist(strsplit(DRM, ','))), by=c('Seq_Code','TYPE')]
	dfdrm[, Gene.codon.number:=dfdrm[, regmatches(DRM,regexpr('[0-9]+', DRM))]]
	dfdrm[, Wild.type:=dfdrm[, regmatches(DRM,regexpr('^[^0-9]*', DRM))]]
	dfdrm[, Mutant:=dfdrm[, regmatches(DRM,regexpr('[^0-9]*$', DRM))]]
	dfdrm			<- subset(dfdrm, !grepl('insertion|deletion',Mutant))
	set(dfdrm, NULL, 'Mutant', dfdrm[, gsub('*','',Mutant,fixed=1)])
	dfdrm			<- dfdrm[, list(Mutant= unlist(strsplit(Mutant,''))), by=c('Seq_Code','TYPE','DRM','Gene.codon.number','Wild.type')]	
	nt2aa			<- as.data.table( read.csv( system.file(package="hivclust", "data",'standard_nt_code.csv.gz'), stringsAsFactors=F ) )
	setnames(nt2aa,c("AA","NTs"),c("Mutant","Mutant.NTs"))
	nt2aa			<- subset(nt2aa, select=c(Mutant,Mutant.NTs))
	dfdrm			<- merge(dfdrm, nt2aa, all.x=1, by="Mutant", allow.cartesian=TRUE)
	dfdrm			<- subset(dfdrm, Mutant!=Wild.type)
	set(dfdrm, NULL, "Mutant.NTs", tolower(dfdrm[,Mutant.NTs]))
	#	get pos of wild type in alignment: align HXB2s
	load( system.file(package="hivclust", "data",'refseq_hiv1_hxb2.rda') )
	tmp				<- hxb2[acxs.st:(acxs.st+ncol(acxs)+100-1L), as.list(as.DNAbin(as.character(HXB2.K03455)))]
	names(tmp)	<- 'HXB2'
	tmp				<- c(as.list(acxs[ which(grepl("B.FR.K03455.1983",rownames(acxs))), ]), tmp)
	write.dna(tmp, file=paste(indir, 'ZA.SubC.HXB2s.fasta', sep='/'), format='fasta')
	#	OK HXB2 is aligned in one stretch :-)
	#	HXB2 pos for PI mutations
	set(dfdrm, NULL,'HXB2.pos', dfdrm[,(as.numeric(Gene.codon.number)-15)*3+2295])
	#	HXB2 pos for RT mutations
	tmp			<- dfdrm[, which(grepl('RT', TYPE))]
	set(dfdrm, tmp, 'HXB2.pos', dfdrm[tmp,(as.numeric(Gene.codon.number)-1)*3+2550])
	#	save all DRM locations for future use
	dr			<- unique(subset(dfdrm, select=c(DRM, Gene.codon.number, Wild.type, Mutant.NTs, HXB2.pos)))
	setnames(dr, c('DRM'), c('DR.name'))
	save(dr, file=paste( CODE.HOME,"/data/AC_drugresistance_201508.rda",sep=''))
	#,'Alignment.nuc.pos'
	set(dr, NULL,'Alignment.nuc.pos', dr[, HXB2.pos-acxs.st+1L])
	tmp			<- subset(dr, select=c(DR.name, Mutant.NTs, Alignment.nuc.pos))
	seq			<- seq.rm.drugresistance(as.character(acxs), tmp, verbose=1, rtn.DNAbin=1 )
	outdir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_SA/ACextpol150831'
	outfile		<- 'ZA_SubC_12432_nDRM.fasta'
	write.dna(seq, paste(outdir, '/', outfile, sep=''), format='fasta', colsep='', nbcol=-1)
	save(seq, file=paste(outdir, '/', gsub('\\.fasta','\\.R',outfile), sep=''))	
}
######################################################################################
project.ACpolext.intros<- function()
{	  	
	indir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_SA/ACpolext160209'
	#	David's FastTree2 + LSD tree
	infile		<- 'ViroIntro_noDRM_16442_20160208_F2_LSD.nex'
	ph			<- read.nexus(file.path(indir, infile))
	phd			<- data.table(TAXA=ph$tip.label)
	phd[, LOC:= sapply(strsplit(TAXA,'_'),'[[',1)]
	phd[, LOCL:= factor(LOC, 	levels=c('AC','AO','BW','CAP','CD','CM','MW','MZ','SZ','TZ','ZA','ZM','ZW'),
								labels=c('Africa Centre','Angola','Botswana','CAPRISA','DR Congo','Cameroon','Malawi','Mozambique','Swaziland','Tanzania','South Africa','Zambia','Zimbabwe'))]
	phd[, SEQYR:= as.numeric(regmatches(TAXA,regexpr('[0-9]+$',TAXA)))]		
	
}
######################################################################################
project.ACpolext.Fasttree.160209<- function()
{	  	
	indir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_SA/ACpolext160209'
	infiles		<- data.table(FILE=list.files(indir, pattern='_ft2.newick$', full.names=1))
	require(phytools)
	i			<- 1
	ph			<- read.tree(infiles[i,FILE])	
	tmp			<- which(grepl('HXB2',ph$tip.label))
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])	
	ph			<- ladderize(ph)
	
	phd			<- data.table(TAXA= ph$tip.label)
	phd[, ZA:= grepl('^AC_|^ZA|^CAPR',TAXA)]
	phd[, TIP.CLR:= 'black']
	set(phd, phd[, which(ZA)], 'TIP.CLR','red')
	
	pdf(file= infiles[i,paste(FILE,'.pdf',sep='')], width=40, height=700)
	plot(ph, tip.color= phd[, TIP.CLR], cex=0.5, adj=1)
	dev.off()	
	write.tree(ph, file=infiles[i,paste(gsub('.newick','_reroot.newick',FILE),sep='')])
	
	save(ph, phd, file=infiles[i,paste(FILE,'.rda',sep='')])	
}
######################################################################################
project.ACpolext.trees.inspect<- function()
{	  	
	indir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_SA/ACpolext150913'
	infiles		<- data.table(FILE=list.files(indir, pattern='^ExaML_result', full.names=1))
	require(phytools)
	i			<- 1
	ph			<- read.tree(infiles[i,FILE])	
	tmp			<- which(grepl('HXB2',ph$tip.label))
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
	
	ph			<- ladderize(ph)
	phd			<- data.table(TAXA= ph$tip.label)
	phd[, AC:= grepl('AC_',TAXA)]
	phd[, TIP.CLR:= 'black']
	set(phd, phd[, which(AC)], 'TIP.CLR','red')
	
	pdf(file= infiles[i,paste(FILE,'.pdf',sep='')], width=40, height=250)
	plot(ph, tip.color= phd[, TIP.CLR], cex=0.5, adj=1)
	dev.off()	
	save(ph, phd, file=infiles[i,paste(FILE,'.R',sep='')])
	
	infile		<- "ZA_SubC_12432_nDRM"
	signat.in	<- '150831' 
	signat.out	<- '150831'
	#	ExaML bootstrap args
	bs.from		<- 1
	bs.to		<- 1
	bs.n		<- 500
	outdir		<- indir
	
	#args.parser	<- paste("-m DNA -q",PARTITION)
	args.parser	<- paste("-m DNA")
	args.examl	<- "-f d -m GAMMA"	
	cmd			<- hivc.cmd.examl.bootstrap(indir, infile, signat.in, signat.out, bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, opt.bootstrap.by="codon", args.parser=args.parser, args.examl=args.examl, tmpdir.prefix="examl")					
	dummy		<- lapply(cmd, function(x)
			{				
				#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
				x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=34, hpc.q="pqeph", hpc.mem="450mb", hpc.nproc=1)
				signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
				outfile	<- paste("exa",signat,sep='.')
				#cat(x)
				hivc.cmd.hpccaller(outdir, outfile, x)
				Sys.sleep(1)
			})	
}
######################################################################################
project.ACpolext.examl<- function()
{	  	
	require(big.phylo)
	indir		<- DATA
	infile		<- "ZA.SubC.12294.alnnoDRM"	
	infile		<- "ZA_aln4noDRM_160120"
	#signat.in   <- signat.out   <- '150913'
	#	ExaML bootstrap args
	bs.from		<- 0
	bs.to		<- 2
	bs.n		<- 500
	outdir		<- indir
	
	#args.parser	<- paste("-m DNA -q",PARTITION)
	args.parser	<- paste("-m DNA")
	#args.examl	<- "-f d -m GAMMA"		#	 -- ran for weeks
	#args.examl	<- "-f o -D -m GAMMA"	#	 -- not followed up until 'default' worked
	args.examl	<- "-f d -D -m GAMMA"	#	 -- this is the default that worked in 24 hours	
	cmd			<- cmd.examl.bootstrap(indir, infile, bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, opt.bootstrap.by="codon", args.parser=args.parser, args.examl=args.examl, tmpdir.prefix="examl")					
	dummy		<- lapply(cmd, function(x)
			{				
				#x		<- cmd.hpcwrapper(x, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
				#x		<- cmd.hpcwrapper(x, hpc.walltime=79, hpc.q="pqeph", hpc.mem="1800mb", hpc.nproc=1)
				x		<- cmd.hpcwrapper(x, hpc.walltime=79, hpc.q="pqeelab", hpc.mem="5800mb", hpc.nproc=1)
				signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
				outfile	<- paste("exa",signat,sep='.')
				cat(x)
				cmd.hpccaller(outdir, outfile, x)
				Sys.sleep(1)
			})	
}