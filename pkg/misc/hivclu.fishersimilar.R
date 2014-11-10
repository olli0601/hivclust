######################################################################################
project.athena.Fisheretal.YX.Trm.Region<- function(YX, clumsm.info)
{
	t.group	<- merge( data.table(Patient=YX[, unique(t.Patient)]), unique(subset( clumsm.info, select=c(Patient, RegionHospital, Trm))), by='Patient' )
	#	Exposure group for transmitter - MSM or Other or NA
	set(t.group, t.group[, which(Trm!='MSM')], 'Trm', 'OTH' )
	set(t.group, NULL, 'Trm', as.factor(as.character(t.group[,Trm])) )
	set(t.group, NULL, 'RegionHospital', as.factor(as.character(t.group[,RegionHospital])) )
	tmp	<- table(t.group[,Trm])
	cat(paste('\nExposure group categories for transmitter', paste(names(tmp), tmp, collapse=', ', sep='=')))
	#	Region code for transmitter
	tmp	<- table(t.group[,RegionHospital])
	cat(paste('\nRegion group categories for transmitter', paste(names(tmp), tmp, collapse=', ', sep='=')))
	setnames(t.group, colnames(t.group), paste('t.',colnames(t.group),sep=''))
	#	
	YX	<- merge(YX, t.group, by='t.Patient')
	#
	t.group	<- merge( data.table(Patient=YX[, unique(Patient)]), unique(subset( clumsm.info, select=c(Patient, RegionHospital, Trm))), by='Patient' )
	#	Exposure group for infected - MSM or Other or NA
	set(t.group, t.group[, which(Trm!='MSM')], 'Trm', 'OTH' )
	set(t.group, NULL, 'Trm', as.factor(as.character(t.group[,Trm])) )
	set(t.group, NULL, 'RegionHospital', as.factor(as.character(t.group[,RegionHospital])) )
	tmp	<- table(t.group[,Trm])
	cat(paste('\nExposure group categories for infected', paste(names(tmp), tmp, collapse=', ', sep='=')))
	#	Region code for infected
	tmp	<- table(t.group[,RegionHospital])
	cat(paste('\nRegion group categories for infected', paste(names(tmp), tmp, collapse=', ', sep='=')))
	#
	YX	<- merge(YX, t.group, by='Patient')
	YX
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
	i.diag		<- merge( subset(i.diag, select=c(Patient, t.period)), unique(subset(df.tpairs, select=c(Patient,t.Patient))), by='Patient')
	cat(paste('\nReturn X for #t.Patients=',i.diag[, length(unique(t.Patient))]))	
	i.diag
}
######################################################################################
project.athena.Fisheretal.YX.calendarperiod.byweight<- function(YX, clumsm.info, t.period= 0.25, c.nperiod= 4)
{
	i.diag		<- YX[, list(w=sum(w)), by='Patient']
	i.diag		<- merge( i.diag, unique(subset(clumsm.info, select=c(Patient, AnyPos_T1))), by='Patient' )					
	setkey(i.diag, AnyPos_T1)
	i.diag[, nsum:= cumsum(w)]
	
	i.cgroup.n	<- ceiling( i.diag[ nrow(i.diag), nsum] / c.nperiod )	
	cat(paste('\nTarget group size is',i.cgroup.n))
	i.diag[, t.period:= cut(i.diag[,nsum], breaks= seq(from= 0, by=i.cgroup.n, length.out=c.nperiod+1), labels= seq_len(c.nperiod))]
	tmp			<- i.diag[, list(n=sum(w), min=min(AnyPos_T1), max=max(AnyPos_T1)),by='t.period']
	cat(paste('\nGroup sizes are=',tmp[, paste(n, collapse=', ',sep='')]))
	cat(paste('\nGroup min diag times are=',tmp[, paste(min, collapse=', ',sep='')]))
	cat(paste('\nGroup max diag times are=',tmp[, paste(max, collapse=', ',sep='')]))
	
	YX		<- merge( YX, subset(i.diag, select=c(Patient, t.period)), by='Patient')
	cat(paste('\nReturn X for #t.Patients=',YX[, length(unique(t.Patient))]))	
	YX
}
######################################################################################
project.athena.Fisheretal.X.followup<- function(X.incare, clumsm.info, df.immu, t.period= 0.25, t.endctime=2013.0)
{
	#	get follow.up time periods for every potential transmitter
	#	follow.up timeline starts at AnyPos_T1 and ends at DateDied or t.endctime	
	follow		<- subset(df.immu, select=c(Patient, PosCD4))
	setkey(follow, Patient)
	follow		<- merge( data.table(Patient= X.incare[, unique(t.Patient)]),  follow, by='Patient' )	
	tmp			<- subset( clumsm.info, select=c(Patient, AnyPos_T1, DateDied) )
	setkey(tmp, Patient)
	follow		<- merge(follow, unique(tmp), by='Patient')	
	set(follow, NULL, 'PosCD4', hivc.db.Date2numeric(follow[,PosCD4]))
	set(follow, follow[, which(is.na(DateDied))], 'DateDied', t.endctime)
	follow		<- merge(follow, follow[, list(PosCD4_T1= min(PosCD4)),by='Patient'], by='Patient')
	#	add fake first CD4 times at AnyPos_T1 so that follow up counts also the time to first CD4
	follow		<- rbind( follow, subset(follow, PosCD4_T1>AnyPos_T1)[, list(PosCD4=AnyPos_T1[1], AnyPos_T1=AnyPos_T1[1], DateDied=DateDied[1], PosCD4_T1=PosCD4_T1[1]), by='Patient'] )
	setkey(follow, Patient, PosCD4)
	#
	#	get follow up quantiles that do not depend on the time periods t
	#
	follow.q	<- follow[, {
				if(length(PosCD4)>1)
					tmp	<- diff(PosCD4)
				else
					tmp	<- NA_real_
				list(fw.up.mean= mean(tmp), fw.up.med= median(tmp), fw.up.mx=max(tmp))
			},by='Patient']
	#plot(follow.q[, fw.up.med], follow.q[, fw.up.mx], pch=18)
	#hist(follow[, log(fw.up.mx)])
	follow.q.mx		<- round( follow.q[, quantile(fw.up.mx, p=c(0, 0.75, 0.95, 1), na.rm=TRUE)], d=3)
	cat(paste('\nCD4 max follow up quantiles are=',paste(follow.q.mx, collapse=' ', sep=''), sep=''))
	follow.q.med	<- round( follow.q[, quantile(fw.up.med, p=c(0, 0.75, 0.95, 1), na.rm=TRUE)], d=3)
	cat(paste('\nCD4 median follow up quantiles are=',paste(follow.q.med, collapse=' ', sep=''), sep=''))
	follow.q.mx[1]	<- follow.q.med[1]<- -1
	follow.q.mx[4]	<- follow.q.med[4]<- 100
	#
	#	set up follow up timeline per patient
	#
	follow.t	<- subset(follow, select=c(Patient, AnyPos_T1, DateDied))
	setkey(follow.t, Patient)
	follow.t	<- unique(follow.t)
	set(follow.t, NULL, 'AnyPos_T1', follow.t[, floor(AnyPos_T1) + round( (AnyPos_T1%%1)*100 %/% (t.period*100) ) * t.period] )
	set(follow.t, NULL, 'DateDied', follow.t[, floor(DateDied) + round( (DateDied%%1)*100 %/% (t.period*100) ) * t.period] )	
	follow.t	<- follow.t[, list(t= seq(AnyPos_T1, DateDied-t.period, by=t.period)),by='Patient']
	#
	#	determine fw.up.mx per patient up to time period t
	#
	follow.t		<- merge(follow, follow.t, by='Patient', allow.cartesian=TRUE)
	cat(paste('\nmerge gives nrows=',nrow(follow.t)))
	follow.t		<- follow.t[,	{
										tmp<- PosCD4[which(PosCD4<=t)]
										list(	fw.up.mx= ifelse(length(tmp)>1, max(diff(tmp)), NA_real_),
												fw.up.med= ifelse(length(tmp)>1, median(diff(tmp)), NA_real_)
												)
									},		by=c('Patient','t')]
	follow.t		<- follow.t[, {
										tmp<- which(!is.na(fw.up.mx))
										if(!length(tmp))
											ans<- list(t=t, fw.up.mx=NA_real_, fw.up.med=NA_real_)
										else if(tmp[1]==1)
											ans<- list(t=t, fw.up.mx=fw.up.mx, fw.up.med=fw.up.med)
										else
											ans<- list(t=t, fw.up.mx=c(rep(fw.up.mx[tmp[1]], tmp[1]-1), fw.up.mx[tmp]), fw.up.med=c(rep(fw.up.med[tmp[1]], tmp[1]-1), fw.up.med[tmp]))
										ans						
									}, by='Patient']	
							
	#set(follow.t, NULL, 'fw.up.mx', cut( follow.t[, fw.up.mx], breaks=follow.q.mx, labels=c('<=75pc','<=95pc','>95pc') )		)
	#set(follow.t, NULL, 'fw.up.med', cut( follow.t[, fw.up.med], breaks=follow.q.med, labels=c('<=75pc','<=95pc','>95pc') )		)
	setnames(follow.t, 'Patient', 't.Patient')
	cat(paste('\nreturn X for #t.Patient=',follow.t[, length(unique(t.Patient))]))	
	X.incare		<- merge(X.incare, follow.t, by=c('t','t.Patient'), all.x=1)	
	X.incare
}
######################################################################################
project.athena.Fisheretal.X.nocontact<- function(X.incare, df.viro, df.immu, df.tpairs, clumsm.info, contact.grace=0.5, t.period=0.25, t.endctime= 2013.)
{	
	contact		<- merge( data.table(Patient=df.tpairs[, unique(t.Patient)]), unique(subset(clumsm.info,select=c(Patient, DateLastContact, DateDied))), by='Patient')
	set(contact, contact[, which(is.na(DateDied) | DateDied>t.endctime)], 'DateDied', t.endctime)
	tmp			<- contact[, which(DateLastContact>=DateDied)]
	set(contact, tmp, 'DateLastContact', contact[tmp, DateDied])
	tmp			<- merge(data.table(Patient=df.tpairs[, unique(t.Patient)]), subset(df.viro, select=c(Patient, PosRNA)), by='Patient')
	set(tmp, NULL, 'PosRNA', hivc.db.Date2numeric(tmp[,PosRNA]))
	tmp			<- tmp[, list(PosRNA_TL=max(PosRNA)), by='Patient']	
	contact		<- merge(contact, tmp, by='Patient', all.x=1)
	tmp			<- merge(data.table(Patient=df.tpairs[, unique(t.Patient)]), subset(df.immu, select=c(Patient, PosCD4)), by='Patient')
	set(tmp, NULL, 'PosCD4', hivc.db.Date2numeric(tmp[,PosCD4]))
	tmp			<- tmp[, list(PosCD4_TL=max(PosCD4)), by='Patient']	
	contact		<- merge(contact, tmp, by='Patient', all.x=1)
	#contact[, summary(abs(PosRNA_TL-PosCD4_TL))]
	set(contact, NULL, 'PosRNA_TL', contact[,PosRNA_TL+contact.grace])
	set(contact, NULL, 'PosCD4_TL', contact[,PosCD4_TL+contact.grace])
	#	check if next step OK
	tmp			<- contact[, list(allNA= is.na(DateLastContact) & is.na(DateDied) & is.na(PosRNA_TL) & is.na(PosCD4_TL)), by='Patient']
	if(tmp[,any(allNA)])	stop('unexpected NAs')	
	contact		<- contact[, 	list(	DateDied=DateDied, DateLastContact.old=DateLastContact, 
										PosRNA_TL=PosRNA_TL, PosCD4_TL=PosCD4_TL,
										DateLastContact=ifelse(all(is.na(c(DateLastContact, PosRNA_TL, PosCD4_TL))), DateDied, min( c( max(c(DateLastContact, PosRNA_TL, PosCD4_TL), na.rm=TRUE), DateDied), na.rm=TRUE))	), by='Patient']
	cat(paste('\nsetting stricter DateLastContact for n=',contact[, length(which(DateLastContact.old>DateLastContact))]))
	cat(paste('\nsetting relaxed DateLastContact for n=',contact[, length(which(DateLastContact.old<DateLastContact))]))
	if(nrow(subset( contact, DateDied<DateLastContact )))	stop('unexpected DateDied<DateLastContact')
	set(contact, NULL, 'DateLastContact', contact[, floor(DateLastContact) + ceiling( (DateLastContact%%1)*100 %/% (t.period*100) ) * t.period] )
	set(contact, NULL, 'DateDied', contact[, floor(DateDied) + round( (DateDied%%1)*100 %/% (t.period*100) ) * t.period] )
	contact		<- subset(contact, DateLastContact<DateDied)	
	contact		<- contact[, list(t= seq(DateLastContact, DateDied-t.period, by=t.period), contact='No'),by='Patient']	
	setnames(contact, 'Patient','t.Patient')
	X.incare	<- merge(X.incare, contact, by=c('t.Patient','t'), all.x=1)
	set(X.incare, X.incare[,which(is.na(contact))],'contact','Yes')
	if(X.incare[, length(which(is.na(stage)))])	stop('unexpected NA stage')
	X.incare
}
######################################################################################
project.athena.Fisheretal.X.time.diag2firstVLandCD4<- function(df.tpairs, clumsm.info, df.viro, df.immu, t2.care.t1.q=c(0.25,0.5))
{
	#	recent follow up: Max time of follow up either by CD4 or VL within the first 12 months		
	follow	<- merge( data.table(Patient=df.tpairs[, unique(t.Patient)]),unique(subset(clumsm.info, select=c(Patient, AnyPos_T1))),by='Patient' )
	viro	<- subset( df.viro, select=c(Patient, PosRNA) )
	viro	<- merge(viro, follow, by='Patient')
	set(viro, NULL, 'PosRNA', hivc.db.Date2numeric(viro[,PosRNA]))
	tmp		<- viro[, list(t2.vl.t1= min(PosRNA)-AnyPos_T1[1]) ,by='Patient']	
	immu	<- subset( df.immu, select=c(Patient, PosCD4) )
	immu	<- merge(immu, follow, by='Patient')
	set(immu, NULL, 'PosCD4', hivc.db.Date2numeric(immu[,PosCD4]))
	immu	<- subset(immu, PosCD4>=AnyPos_T1)
	follow	<- immu[, list(t2.cd4.t1= min(PosCD4)-AnyPos_T1[1]) ,by='Patient']
	follow	<- merge(follow, tmp, by='Patient')
	follow	<- merge(follow, follow[, list(t2.care.t1= max(t2.cd4.t1,t2.vl.t1)), by='Patient'], by='Patient')	
	tmp		<- cut(follow[, t2.care.t1], breaks=c(0, t2.care.t1.q, follow[, max(t2.care.t1)+1]), right=FALSE, label= c('other',paste('>=',t2.care.t1.q,'y',sep='')))
	follow[, t2.care.t1.c:= tmp]
	set(follow, NULL, 't2.care.t1', tmp)	
	tmp		<- round( cumsum(rev(table(follow[,t2.care.t1]))) / nrow(follow), d=3 )
	cat(paste('\ncumulative right tail probability for recent follow up quantiles ',paste( names(tmp), ' = ',tmp*100, '%',sep='', collapse=', ')))
	setnames(follow, 'Patient', 't.Patient')
	follow	<- subset(follow, select=c(t.Patient, t2.care.t1))
	cat(paste('\nreturn X for #t.Patient=',follow[, length(unique(t.Patient))]))
	follow
}
######################################################################################
project.athena.Fisheretal.X.time.diag2suppressed<- function(df.tpairs, clumsm.info, df.viro, lRNA.suppressed= log10(1e3), t2.vl.supp.p=c(0.1, 0.25))
{
	#	do category because it is not so easy to deal with those not suppressed otherwise
	tmp		<- merge( data.table( Patient=df.tpairs[, unique(t.Patient)] ), unique(subset(clumsm.info, select=c(Patient, AnyT_T1, AnyPos_T1))), by='Patient' )
	viro	<- subset( df.viro, select=c(Patient, PosRNA, lRNA) )
	tmp2	<- setdiff( tmp[,Patient], viro[, unique(Patient)])
	cat(paste('\nPotential transmitters for which we have not a single VL',length(tmp2)))
	viro	<- merge(viro, tmp, by='Patient')	
	set(viro, NULL, 'PosRNA', hivc.db.Date2numeric(viro[,PosRNA]))
	setnames(viro, 'Patient', 't.Patient')
	#	assume we have at least one lRNA for every potential transmitter
	#	so everyone not in the following subset is either not yet on ART or has not yet suppressed VL
	viro.supp		<- subset(viro, PosRNA>AnyT_T1 & lRNA<=lRNA.suppressed)	
	viro.supp		<- viro.supp[, list(t2.vl.supp= min(PosRNA)-AnyPos_T1[1]),by='t.Patient']	
	#hist(viro[,t2.vl.supp], breaks=50)
	#t2.vl.supp.q	<- round(quantile(viro.supp[,t2.vl.supp], prob=c(0,t2.vl.supp.p,1)), d=3)
	#cat(paste('\nquantiles for time to viral suppression after diagnosis',paste(names(t2.vl.supp.q), t2.vl.supp.q, collapse=', ',sep='=')))
	#tmp				<- cut(viro.supp[,t2.vl.supp], breaks= t2.vl.supp.q, right=FALSE, label= c( paste( '<',t2.vl.supp.p*100,'pc',sep='' ), 'other' ) )
	#set(viro.supp,NULL,'t2.vl.supp',tmp)
	viro			<- merge( unique( subset(viro, select=t.Patient) ), viro.supp, by='t.Patient', all.x=1 )
	#set(viro,viro[,which(is.na(t2.vl.supp))],'t2.vl.supp','other')
	cat(paste('\nreturn X for #t.Patient=',viro[, length(unique(t.Patient))]))
	viro
}
######################################################################################
project.athena.Fisheretal.median.followup<- function(df.tpairs, df.viro, df.immu, X.tperiod)
{
	follow		<- subset(df.viro, select=c(Patient, PosRNA))
	setkey(follow, Patient)
	follow		<- merge( data.table(Patient= df.tpairs[, unique(t.Patient)]),  follow, by='Patient' )
	set(follow, NULL, 'PosRNA', hivc.db.Date2numeric(follow[,PosRNA]))
	tmp			<- follow[, list(nRNA=length(PosRNA)) , by='Patient']
	follow.vl	<- merge(follow, subset(tmp, nRNA>1, Patient), by='Patient')
	follow.vl	<- follow.vl[, list(dRNA=median(diff(PosRNA))), by='Patient']
	
	follow		<- subset(df.immu, select=c(Patient, PosCD4))
	setkey(follow, Patient)
	follow		<- merge( data.table(Patient= df.tpairs[, unique(t.Patient)]),  follow, by='Patient' )
	set(follow, NULL, 'PosCD4', hivc.db.Date2numeric(follow[,PosCD4]))
	tmp			<- follow[, list(nCD4=length(PosCD4)) , by='Patient']
	follow.cd4	<- merge(follow, subset(tmp, nCD4>1, Patient), by='Patient')
	follow.cd4	<- follow.cd4[, list(dCD4=median(diff(PosCD4))), by='Patient']
	
	follow		<- merge(follow.vl, follow.cd4, by='Patient')
	tmp			<- follow[, list(d=min(dRNA,dCD4)),by='Patient']
	setnames(tmp, 'Patient', 't.Patient')	
	tmp			<- merge( unique(subset(X.tperiod, select=c(t.Patient, t.period))), tmp, by='t.Patient' )
	setkey(tmp, t.period)
	tmp
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
project.athena.Fisheretal.tripletweights<- function(Y.score, save.file=NA)
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
		triplet.weight				<- subset(Y.score, select=c(t, Patient, t.Patient))	
		setkey(triplet.weight, t, Patient, t.Patient)
		tmp						<- copy(triplet.weight)
		setnames(tmp, colnames(tmp), paste('q.',colnames(tmp),sep='') )
		triplet.weight				<- triplet.weight[,	{ 	
															z	<- tmp[, which(q.t==t & q.Patient==t.Patient & q.t.Patient==Patient)]													
															list(equal.triplet.idx= ifelse(!length(z), NA_integer_, z))						
														}	, by=c('t','Patient','t.Patient')]
		triplet.weight				<- subset(triplet.weight, !is.na(equal.triplet.idx))										
		if(!is.na(save.file))
		{
			cat(paste('\nsave triplet.weight to file=',save.file))
			save(triplet.weight, file=save.file)
		}
	}
	triplet.weight
}
######################################################################################
project.athena.Fisheretal.YX.part2<- function(YX.part1, df.all, df.treatment, df.viro, predict.t2inf, t2inf.args, indir, insignat, indircov, infilecov, infiletree, infile.trm.model, outdir, outfile, cluphy=NULL, cluphy.info=NULL, cluphy.map.nodectime=NULL, df.tpairs.4.rawbrl=NULL, dur.Acute=NULL, rm.zero.score=FALSE, any.pos.grace.yr=3.5, thresh.pcoal=0.5, brl.bwhost.multiplier=1, method.minLowerUWithNegT=TRUE, lRNA.supp=log10(51), t.period=0.25, save.file=NA, save.all=FALSE, resume=1, method='3aa')
{
	if(resume && !is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(!resume || is.na(save.file) || inherits(readAttempt, "try-error"))
	{		
		YX.tpairs	<- subset(YX.part1, select=c(cluster, Patient, t.Patient, FASTASampleCode, t.FASTASampleCode, class))
		setkey(YX.tpairs, Patient, t.Patient, FASTASampleCode, t.FASTASampleCode)
		YX.tpairs	<- unique(YX.tpairs)
		#
		if(is.null(df.tpairs.4.rawbrl))
			df.tpairs.4.rawbrl	<- YX.tpairs	
		#
		#	compute Y scores for potential transmitters in ATHENA.clu for which dated phylogenies are available
		#
		#	BRL	[0,1]: raw branch length between pot transmitter and infected
		method.restrictTPtoRI	<- ifelse(substr(method,1,2)%in%c('3a','3b'),1,0)
		#file					<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'nbrlraw',method,sep='')
		file					<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'nbrlrawmed',method,sep='')
		tmp						<- project.athena.Fisheretal.Y.rawbrl(YX.tpairs, indir, insignat, indircov, infilecov, infiletree, df.tpairs.tptn=df.tpairs.4.rawbrl, save.file=paste(file, '.R', sep=''), resume=resume, plot.file=paste(file, '.pdf', sep=''), method.restrictTPtoRI=method.restrictTPtoRI)
		Y.rawbrl				<- tmp$tpairs
		Y.rawbrl.linked			<- tmp$linked
		Y.rawbrl.unlinked		<- tmp$unlinked
		#	BRL [0,1]: branch length weight between pot transmitter and infected
		plot.file.score			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl',method,'_score','.pdf',sep='')
		plot.file.both			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl',method,'_both','.pdf',sep='')
		plot.file.one			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl',method,'_one','.pdf',sep='')
		#brl.linked.max.brlr= 0.02; brl.linked.min.brl= 1e-4; brl.linked.max.dt= 10; brl.linked.min.dt= 1; 
		if(substr(method,1,2)=='3c')					
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.treatment, brl.linked.max.brlr= 0.04, brl.linked.min.brl= 1e-4, brl.linked.max.dt= 10, brl.linked.min.dt= 1, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one)
		if(substr(method,1,2)=='3d')	#expect central brl	
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.treatment, brl.linked.max.brlr= 0.01, brl.linked.min.brl= 1e-12, brl.linked.max.dt= 10, brl.linked.min.dt= 1, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one)		
		if(substr(method,1,2)=='3j')	#expect central brl	
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.treatment, brl.linked.max.brlr= 0.02, brl.linked.min.brl= 1e-12, brl.linked.max.dt= 10, brl.linked.min.dt= 1, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one)				
		if(substr(method,1,2)=='3i')	#expect central brl	
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.treatment, brl.linked.max.brlr= 0.02, brl.linked.min.brl= 1e-12, brl.linked.max.dt= 10, brl.linked.min.dt= 1, brl.bwhost.multiplier=brl.bwhost.multiplier, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one)
		if(substr(method,1,2)=='3k')	#expect central brl	
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.treatment, brl.linked.max.brlr= 0.02, brl.linked.min.brl= 1e-12, brl.linked.max.dt= 10, brl.linked.min.dt= 1, brl.bwhost.multiplier=brl.bwhost.multiplier, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one)
		if(substr(method,1,2)=='3l')	#expect central brl
			Y.brl				<- project.athena.Fisheretal.Y.brlweight.3l(infile.trm.model, Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.viro, dur.Acute=dur.Acute, t.period=t.period, lRNA.supp=lRNA.supp, plot.file.score=plot.file.score)					
		if(substr(method,1,2)=='3m')	#expect central brl
			Y.brl				<- project.athena.Fisheretal.Y.brlweight.3m(infile.trm.model, Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.viro, dur.Acute=dur.Acute, t.period=t.period, lRNA.supp=lRNA.supp, plot.file.score=plot.file.score)							
		if(substr(method,1,2)=='3n')	#expect central brl
			Y.brl				<- project.athena.Fisheretal.Y.brlweight.3n(infile.trm.model, Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.viro, dur.Acute=dur.Acute, t.period=t.period, lRNA.supp=lRNA.supp, plot.file.score=plot.file.score)									
		if(substr(method,1,2)=='3o')	#expect central brl
			Y.brl				<- project.athena.Fisheretal.Y.brlweight.3o(infile.trm.model, Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.viro, dur.Acute=dur.Acute, t.period=t.period, lRNA.supp=lRNA.supp, plot.file.score=plot.file.score)											
		if(substr(method,1,2)=='3p')	#expect central brl
			Y.brl				<- project.athena.Fisheretal.Y.brlweight.3p(infile.trm.model, Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.viro, dur.Acute=dur.Acute, t.period=t.period, lRNA.supp=lRNA.supp, plot.file.score=plot.file.score)													
		if(substr(method,1,2)=='3e')	#expect large brl				
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.treatment, brl.linked.max.brlr= 0.0175, brl.linked.min.brl= 1e-4, brl.linked.max.dt= 10, brl.linked.min.dt= 1, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one)
		if(substr(method,1,2)=='3f')	#expect small brl				
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.treatment, brl.linked.max.brlr= 0.005, brl.linked.min.brl= 1e-12, brl.linked.max.dt= 10, brl.linked.min.dt= -1, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one)		
		if(substr(method,1,2)=='3g')	#similar to central brl, use dependence on time between sampling dates				
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.treatment, brl.linked.max.brlr= 0.01, brl.linked.min.brl= 1e-12, brl.linked.max.dt= 10, brl.linked.min.dt= -1, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one)		
		if(substr(method,1,2)=='3h')	#similar to central brl, use dependence on time between sampling dates
		{
			save.file.3ha		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl',method,'_convolution','.R',sep='')
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.treatment, brl.linked.max.brlr= 0.01, brl.linked.min.brl= 1e-12, brl.linked.max.dt= 10, brl.linked.min.dt= -1, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one, save.file.3ha=save.file.3ha)		
		}
		if(all(substr(method,1,2)!=c('3c','3d','3e','3f','3g','3i','3j','3k','3l','3m','3n','3o','3p')))	stop('brlweight: method not supported')		
		#	COAL [0,1]: prob that coalescence is within the transmitter
		Y.coal					<- NULL
		if(!is.null(cluphy) & !is.null(cluphy.info) & !is.null(cluphy.map.nodectime)  )
		{
			#	U [0,1]: prob that pot transmitter is still infected at time t. Needed to determine time of infection for transmitter (as quantile of the surival distribution)
			Y.U					<- project.athena.Fisheretal.Y.infectiontime(YX.tpairs, df.all, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, score.set.value=NA, method='for.transmitter', method.minLowerUWithNegT=method.minLowerUWithNegT)			
			file				<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'coalraw',method,'.R',sep='')		
			Y.coal				<- project.athena.Fisheretal.Y.coal(YX.tpairs, df.all, Y.U, cluphy, cluphy.info, cluphy.map.nodectime, coal.t.Uscore.min=0.01, coal.within.inf.grace= 0.25, t.period=t.period, save.file=file, resume=resume, method.minLowerUWithNegT=method.minLowerUWithNegT )			
		}
		#	U [0,1]: prob that pot transmitter is still infected during infection window
		Y.U						<- project.athena.Fisheretal.Y.transmitterinfected(YX.part1)		
		#	screen for likely missed intermediates/sources
		plot.file				<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'missed_',method,'.pdf',sep='')		
		YX.tpairs.select		<- project.athena.Fisheretal.Y.rm.missedtransmitter(YX.tpairs, df.all, Y.brl, Y.U, Y.coal=Y.coal, thresh.pcoal=thresh.pcoal, cut.brl=1e-3, any.pos.grace.yr= any.pos.grace.yr, rm.zero.score=rm.zero.score, plot.file=plot.file, t.period=t.period)
		#	take branch length weight ONLY to derive "score.Y"
		Y.score					<- merge( YX.tpairs.select, subset(Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, telapsed, brl, score.brl.TPp, score.brl.TPd)), by=c('FASTASampleCode','t.FASTASampleCode'))
		#	keep only one (Patient, t.Patient) pair if there are multiple sequences per individual. Take the one with largest brl score
		tmp						<- Y.score[,	{	
													z<- which.min(brl)
													list(FASTASampleCode=FASTASampleCode[z], t.FASTASampleCode=t.FASTASampleCode[z])
												},by=c('Patient', 't.Patient')]	
		Y.score					<- merge( Y.score, subset(tmp, select=c(FASTASampleCode=FASTASampleCode, t.FASTASampleCode)), by=c('FASTASampleCode', 't.FASTASampleCode') )
		set(Y.score, NULL, 'score.Y', Y.score[, score.Y*score.brl.TPp])
		Y.score[, score.brl.TPp:=NULL]
		#
		tmp						<- Y.score[, list(n= length(FASTASampleCode)), by=c('Patient', 't.Patient')]
		if(tmp[,any(n!=1)])	stop('unexpected multiple entries per Patient, t.Patient')
		#	merge over time periods t
		YX						<- merge( YX.part1, subset(Y.score, select=c(FASTASampleCode, t.FASTASampleCode, score.Y, telapsed, brl, score.brl.TPd)), by=c('FASTASampleCode', 't.FASTASampleCode'))
		#	only triplets currently lost because of missing BEAST2 dated clusters
		YX.part1				<- NULL
		print(YX[, table(score.Y>0, class)])
		gc()
		#	add weights 		
		YX						<- project.athena.Fisheretal.YX.weight(YX)
		#	add balanced calendar periods 		
		YX						<- project.athena.Fisheretal.YX.calendarperiod.byweight(YX, df.all, t.period=t.period, c.nperiod= 4)
		#	re-arrange a little
		YX						<- subset(YX, select=	c(	t, t.Patient, Patient, cluster, score.Y, telapsed, brl, 																	#triplets identifiers and Y score
															t.period, stage, U.score, contact, CDCC, lRNA, nlRNA.supp, nlRNA.nsupp, lRNA_T1.supp, CD4, 														#main covariates							
															t.Age, Age, t.RegionHospital, RegionHospital, ART.I, ART.F, ART.A, ART.P, ART.T, ART.pulse, ART.nDrug, ART.nNRT, ART.nNNRT, ART.nPI, ART.nBoost, ART.TDF.EFV.FTC, fw.up.mx, fw.up.med, t2.care.t1, t2.vl.supp, 		#secondary covariates
															t.AnyPos_T1,  AnyPos_T1, t.AnyT_T1, StartTime, StopTime, lRNAc, t.isAcute, t.Trm, Trm,												#other 
															FASTASampleCode, t.FASTASampleCode, w, w.i, w.in, w.t, w.tn, class										#other							
															))
		setkey(YX, t, t.Patient, Patient)		
		#
		if(!is.na(save.file))
		{
			cat(paste('\nsave YX to file', save.file))
			if(!save.all)
				save(YX, file=save.file)
			if(save.all)
				save(YX, YX.tpairs, df.all, Y.brl, Y.U, Y.coal, file=save.file)
		}
	}
	if(!save.all)
		return( YX )
	if(save.all)
		return( list(YX=YX, YX.tpairs=YX.tpairs, df.all=df.all, Y.brl=Y.brl, Y.U=Y.U, Y.coal=Y.coal) )
}
######################################################################################
project.athena.Fisheretal.Y.TSupp<- function(df.TS, df.viro, t.period=0.25, lRNA.supp=3)
{
	#	get time in [meanInfT,  PosSeqT] under which individual is on ART 
	tmp			<- df.TS[, {
								tmp	<-  range(PosSeqT, meanInfT)
								if(is.na(AnyT_T1) | tmp[2]<AnyT_T1)						#no time potentially suppressed	
									tmp[2]<- tmp[1]
								if(!is.na(AnyT_T1) & tmp[1]<AnyT_T1 & tmp[2]>=AnyT_T1)	#only time between ART start and tmp[2] potentially suppressed
									tmp[1]<- AnyT_T1
								list( ts= tmp[1], te=tmp[2]  )
							}, by='DUMMY']
	df.TS		<- merge(subset(df.TS, select=c(DUMMY,Patient)), tmp, by='DUMMY')
	df.TS		<- subset(df.TS, te>ts)
	#	adjust ts, te to t.period intervals and expand
	set(df.TS, NULL, 'ts', df.TS[, floor(ts) + floor( (ts%%1)*100 %/% (t.period*100) ) * t.period] )
	set(df.TS, NULL, 'te', df.TS[, floor(te) + ceiling( (te%%1)*100 %/% (t.period*100) ) * t.period] )
	df.TS		<- df.TS[, list(t= seq(ts, te, by=t.period), Patient=Patient), by='DUMMY']
	#	get viral loads for time [FromT, ToT] under which individual is on ART
	tmp			<- data.table(t.Patient=df.TS[, unique(Patient)])	
	tmp			<- project.athena.Fisheretal.X.viro(tmp, df.viro, t.period=t.period, lRNA.supp=lRNA.supp)
	setnames(tmp, 't.Patient','Patient')	
	df.TS		<- merge( df.TS, subset(tmp, !is.na(lRNA)), by=c('Patient','t'), all.x=TRUE )	
	tmp			<- df.TS[, {
								tmp			<- which(!is.na(lRNA))
								lRNAext		<- lRNA
								if( length(tmp)>0 && tail(tmp,1)!=length(lRNA))	#if last non-NA is not final entry
								{
									lRNAext	<- c(lRNA[1:tail(tmp,1)], rep( lRNA[tail(tmp,1)], length(lRNA)-tail(tmp,1) ))
								}
								list(t=t, lRNAext=lRNAext)
							}, by='DUMMY']
	df.TS		<- merge( df.TS, tmp, by=c('DUMMY','t'))				
	#	count time suppressed
	df.TS[, list( tsupp= t.period*length(which(!is.na(lRNA) & lRNAext<=lRNA.supp))), by='DUMMY']		
}
######################################################################################
project.athena.Fisheretal.Y.brlweight.3m<- function(infile.trm.model, Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.viro, dur.Acute=NULL, t.period=0.25, lRNA.supp=log10(51), plot.file.score=NA)
{
	require(reshape2)
	require(ggplot2)
	require(gamlss)		
	#
	#	compute cdf for pairs given they are unlinked
	#
	setkey(Y.rawbrl.unlinked, brl)
	p.rawbrl.unlinked	<- Y.rawbrl.unlinked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	setkey(Y.rawbrl.linked, brl)
	p.rawbrl.linked		<- Y.rawbrl.linked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	#
	#	add sequence sampling times	
	#
	Y.brl				<- copy(Y.rawbrl)
	Y.brl				<- merge( Y.brl, unique(subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyT_T1, isAcute))), by='FASTASampleCode' )	
	tmp					<- merge( data.table(FASTASampleCode=Y.brl[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	Y.brl				<- merge( Y.brl, tmp, by='t.FASTASampleCode')
	stopifnot( Y.brl[, all(isAcute%in%c('Yes','Maybe'))] )
	#
	#	calculate d_TSeqT
	#
	Y.brl				<- merge( Y.brl, data.table(isAcute=c('Yes', 'Maybe'), meanInfWindow= dur.Acute[c('Yes','Maybe')] ), by='isAcute' )		
	Y.brl[, meanInfT:= PosSeqT-meanInfWindow/365.25]
	Y.brl[, DUMMY:=seq_len(nrow(Y.brl))]
	#	calculate time spent suppressed for transmitter
	df.TS		<- subset(Y.brl, select=c(DUMMY, t.Patient, t.AnyT_T1, meanInfT, t.PosSeqT))
	setnames(df.TS, colnames(df.TS), gsub('t.', '', colnames(df.TS), fixed=TRUE))	
	df.TS		<- project.athena.Fisheretal.Y.TSupp(df.TS, df.viro, t.period=t.period, lRNA.supp=lRNA.supp)
	setnames(df.TS, 'tsupp', 't.tsupp')
	Y.brl		<- merge(Y.brl, df.TS, by='DUMMY', all.x=TRUE)
	set(Y.brl, Y.brl[, which(is.na(t.tsupp))], 't.tsupp', 0.)
	#	calculate time spent suppressed for infected
	df.TS		<- subset(Y.brl, select=c(DUMMY, Patient, AnyT_T1, meanInfT, PosSeqT))	
	df.TS		<- project.athena.Fisheretal.Y.TSupp(df.TS, df.viro, t.period=t.period, lRNA.supp=lRNA.supp)
	Y.brl		<- merge(Y.brl, df.TS, by='DUMMY', all.x=TRUE)
	set(Y.brl, Y.brl[, which(is.na(tsupp))], 'tsupp', 0.)
	cat(paste('\nFound person-years of suppressed evolution, ', Y.brl[, sum(t.tsupp)+sum(tsupp)]))
	#	calculate time evolving in transmitter and  infected
	Y.brl[, t.d_TSeqT:=  abs(t.PosSeqT-meanInfT)-t.tsupp ]
	set(Y.brl, Y.brl[, which(t.d_TSeqT<0)], 't.d_TSeqT', 0.)	
	Y.brl[, i.d_TSeqT:=  abs(PosSeqT-meanInfT)-tsupp ]
	set(Y.brl, Y.brl[, which(i.d_TSeqT<0)], 'i.d_TSeqT', 0.)	
	Y.brl[, d_TSeqT:=t.d_TSeqT+i.d_TSeqT]
	cat(paste('\nFound person-years of unsuppressed evolution, ', Y.brl[, sum(d_TSeqT)]))
	#
	#	load Gamma model and compute likelihood
	#
	load(infile.trm.model)		#expect "trm.pol.GA" "trm.pol.nA"   
	Y.brl[, mu:= predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='mu', type='link')]
	Y.brl[, sigma:= predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='sigma', type='link')]	
	Y.brl[, score.brl.TPd:= Y.brl[, dGA(brl, mu=mu, sigma=sigma)]]
	Y.brl[, score.brl.TPp:= Y.brl[, score.brl.TPd]]
	#
	#	set branch length quantiles
	#	
	Y.brl[, score.brl.TN:= p.rawbrl.unlinked(brl)]
	Y.brl[, score.brl.TP:= p.rawbrl.linked(brl)]
	#	
	if(!is.na(plot.file.score))
	{
		#ggplot(Y.brl, aes(x=score.brl.TPd)) + geom_histogram()		
	}
	#
	Y.brl		<- subset( Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, score.brl.TPd, score.brl.TPp, score.brl.TN, score.brl.TP, d_TSeqT, brl) )
	setnames(Y.brl, 'd_TSeqT','telapsed')
	cat(paste('\nReturn brl score for #t.Patient seqs=',Y.brl[, length(unique(t.FASTASampleCode))]))
	Y.brl		
}
######################################################################################
project.athena.Fisheretal.Y.brlweight.3o<- function(infile.trm.model, Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.viro, dur.Acute=NULL, t.period=0.25, lRNA.supp=log10(51), plot.file.score=NA)
{
	require(reshape2)
	require(ggplot2)
	require(gamlss)		
	#
	#	compute cdf for pairs given they are unlinked
	#
	setkey(Y.rawbrl.unlinked, brl)
	p.rawbrl.unlinked	<- Y.rawbrl.unlinked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	setkey(Y.rawbrl.linked, brl)
	p.rawbrl.linked		<- Y.rawbrl.linked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	#
	#	add sequence sampling times	
	#
	Y.brl				<- copy(Y.rawbrl)
	Y.brl				<- merge( Y.brl, unique(subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyT_T1, isAcute))), by='FASTASampleCode' )	
	tmp					<- merge( data.table(FASTASampleCode=Y.brl[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	Y.brl				<- merge( Y.brl, tmp, by='t.FASTASampleCode')
	stopifnot( Y.brl[, all(isAcute%in%c('Yes','Maybe'))] )
	#
	#	calculate d_TSeqT
	#
	Y.brl				<- merge( Y.brl, data.table(isAcute=c('Yes', 'Maybe'), meanInfWindow= dur.Acute[c('Yes','Maybe')] ), by='isAcute' )		
	Y.brl[, meanInfT:= PosSeqT-meanInfWindow/365.25]
	Y.brl[, DUMMY:=seq_len(nrow(Y.brl))]
	#	calculate time spent suppressed for transmitter
	df.TS		<- subset(Y.brl, select=c(DUMMY, t.Patient, t.AnyT_T1, meanInfT, t.PosSeqT))
	setnames(df.TS, colnames(df.TS), gsub('t.', '', colnames(df.TS), fixed=TRUE))	
	df.TS		<- project.athena.Fisheretal.Y.TSupp(df.TS, df.viro, t.period=t.period, lRNA.supp=lRNA.supp)
	setnames(df.TS, 'tsupp', 't.tsupp')
	Y.brl		<- merge(Y.brl, df.TS, by='DUMMY', all.x=TRUE)
	set(Y.brl, Y.brl[, which(is.na(t.tsupp))], 't.tsupp', 0.)
	#	calculate time spent suppressed for infected
	df.TS		<- subset(Y.brl, select=c(DUMMY, Patient, AnyT_T1, meanInfT, PosSeqT))	
	df.TS		<- project.athena.Fisheretal.Y.TSupp(df.TS, df.viro, t.period=t.period, lRNA.supp=lRNA.supp)
	Y.brl		<- merge(Y.brl, df.TS, by='DUMMY', all.x=TRUE)
	set(Y.brl, Y.brl[, which(is.na(tsupp))], 'tsupp', 0.)
	cat(paste('\nFound person-years of suppressed evolution, ', Y.brl[, sum(t.tsupp)+sum(tsupp)]))
	#	calculate time evolving in transmitter and  infected
	Y.brl[, t.d_TSeqT:=  abs(t.PosSeqT-meanInfT)-t.tsupp ]
	set(Y.brl, Y.brl[, which(t.d_TSeqT<0)], 't.d_TSeqT', 0.)	
	Y.brl[, i.d_TSeqT:=  abs(PosSeqT-meanInfT)-tsupp ]
	set(Y.brl, Y.brl[, which(i.d_TSeqT<0)], 'i.d_TSeqT', 0.)	
	Y.brl[, d_TSeqT:=t.d_TSeqT+i.d_TSeqT]
	cat(paste('\nFound person-years of unsuppressed evolution, ', Y.brl[, sum(d_TSeqT)]))
	#
	#	load Gamma model and compute likelihood
	#
	load(infile.trm.model)		#expect "trm.pol.GA" "trm.pol.nA"   
	Y.brl[, mu:= predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='mu', type='link')]
	Y.brl[, sigma:= predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='sigma', type='link')]	
	Y.brl[, score.brl.TPd:= Y.brl[, dGA(brl, mu=mu, sigma=sigma)]]
	Y.brl[, score.brl.TPp:= Y.brl[, score.brl.TPd]]
	#
	#	set branch length quantiles
	#	
	Y.brl[, score.brl.TN:= p.rawbrl.unlinked(brl)]
	Y.brl[, score.brl.TP:= p.rawbrl.linked(brl)]
	#	
	if(!is.na(plot.file.score))
	{
		#ggplot(Y.brl, aes(x=score.brl.TPd)) + geom_histogram()		
	}
	#
	Y.brl		<- subset( Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, score.brl.TPd, score.brl.TPp, score.brl.TN, score.brl.TP, d_TSeqT, brl) )
	setnames(Y.brl, 'd_TSeqT','telapsed')
	cat(paste('\nReturn brl score for #t.Patient seqs=',Y.brl[, length(unique(t.FASTASampleCode))]))
	Y.brl		
}
######################################################################################
project.athena.Fisheretal.Y.brlweight.3p<- function(infile.trm.model, Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.viro, dur.Acute=NULL, t.period=0.25, lRNA.supp=log10(51), plot.file.score=NA)
{
	require(reshape2)
	require(ggplot2)
	require(gamlss)		
	#
	#	compute cdf for pairs given they are unlinked
	#
	setkey(Y.rawbrl.unlinked, brl)
	p.rawbrl.unlinked	<- Y.rawbrl.unlinked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	setkey(Y.rawbrl.linked, brl)
	p.rawbrl.linked		<- Y.rawbrl.linked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	#
	#	add sequence sampling times	
	#
	Y.brl				<- copy(Y.rawbrl)
	Y.brl				<- merge( Y.brl, unique(subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyT_T1, isAcute))), by='FASTASampleCode' )	
	tmp					<- merge( data.table(FASTASampleCode=Y.brl[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	Y.brl				<- merge( Y.brl, tmp, by='t.FASTASampleCode')
	stopifnot( Y.brl[, all(isAcute%in%c('Yes','Maybe'))] )
	#
	#	calculate d_TSeqT
	#
	Y.brl				<- merge( Y.brl, data.table(isAcute=c('Yes', 'Maybe'), meanInfWindow= dur.Acute[c('Yes','Maybe')] ), by='isAcute' )		
	Y.brl[, meanInfT:= PosSeqT-meanInfWindow/365.25]
	Y.brl[, DUMMY:=seq_len(nrow(Y.brl))]
	#	calculate time spent suppressed for transmitter
	df.TS		<- subset(Y.brl, select=c(DUMMY, t.Patient, t.AnyT_T1, meanInfT, t.PosSeqT))
	setnames(df.TS, colnames(df.TS), gsub('t.', '', colnames(df.TS), fixed=TRUE))	
	df.TS		<- project.athena.Fisheretal.Y.TSupp(df.TS, df.viro, t.period=t.period, lRNA.supp=lRNA.supp)
	setnames(df.TS, 'tsupp', 't.tsupp')
	Y.brl		<- merge(Y.brl, df.TS, by='DUMMY', all.x=TRUE)
	set(Y.brl, Y.brl[, which(is.na(t.tsupp))], 't.tsupp', 0.)
	#	calculate time spent suppressed for infected
	df.TS		<- subset(Y.brl, select=c(DUMMY, Patient, AnyT_T1, meanInfT, PosSeqT))	
	df.TS		<- project.athena.Fisheretal.Y.TSupp(df.TS, df.viro, t.period=t.period, lRNA.supp=lRNA.supp)
	Y.brl		<- merge(Y.brl, df.TS, by='DUMMY', all.x=TRUE)
	set(Y.brl, Y.brl[, which(is.na(tsupp))], 'tsupp', 0.)
	cat(paste('\nFound person-years of suppressed evolution, ', Y.brl[, sum(t.tsupp)+sum(tsupp)]))
	#	calculate time evolving in transmitter and  infected
	Y.brl[, t.d_TSeqT:=  abs(t.PosSeqT-meanInfT)-t.tsupp ]
	set(Y.brl, Y.brl[, which(t.d_TSeqT<0)], 't.d_TSeqT', 0.)	
	Y.brl[, i.d_TSeqT:=  abs(PosSeqT-meanInfT)-tsupp ]
	set(Y.brl, Y.brl[, which(i.d_TSeqT<0)], 'i.d_TSeqT', 0.)	
	Y.brl[, d_TSeqT:=t.d_TSeqT+i.d_TSeqT]
	cat(paste('\nFound person-years of unsuppressed evolution, ', Y.brl[, sum(d_TSeqT)]))
	#
	#	load Gamma model and compute likelihood
	#
	load(infile.trm.model)		#expect "trm.pol.GA" "trm.pol.nA"   
	Y.brl[, mu:= predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='mu', type='link')]
	Y.brl[, sigma:= predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='sigma', type='link')]	
	Y.brl[, score.brl.TPd:= Y.brl[, dGA(brl, mu=mu, sigma=sigma)]]
	Y.brl[, score.brl.TPp:= Y.brl[, score.brl.TPd]]
	#
	#	set branch length quantiles
	#	
	Y.brl[, score.brl.TN:= p.rawbrl.unlinked(brl)]
	Y.brl[, score.brl.TP:= p.rawbrl.linked(brl)]
	#	
	if(!is.na(plot.file.score))
	{
		#ggplot(Y.brl, aes(x=score.brl.TPd)) + geom_histogram()		
	}
	#
	Y.brl		<- subset( Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, score.brl.TPd, score.brl.TPp, score.brl.TN, score.brl.TP, d_TSeqT, brl) )
	setnames(Y.brl, 'd_TSeqT','telapsed')
	cat(paste('\nReturn brl score for #t.Patient seqs=',Y.brl[, length(unique(t.FASTASampleCode))]))
	Y.brl		
}
######################################################################################
project.athena.Fisheretal.Y.brlweight.3n<- function(infile.trm.model, Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.viro, dur.Acute=NULL, t.period=0.25, lRNA.supp=log10(51), plot.file.score=NA)
{
	require(reshape2)
	require(ggplot2)
	require(gamlss)		
	#
	#	compute cdf for pairs given they are unlinked
	#
	setkey(Y.rawbrl.unlinked, brl)
	p.rawbrl.unlinked	<- Y.rawbrl.unlinked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	setkey(Y.rawbrl.linked, brl)
	p.rawbrl.linked		<- Y.rawbrl.linked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	#
	#	add sequence sampling times	
	#
	Y.brl				<- copy(Y.rawbrl)
	Y.brl				<- merge( Y.brl, unique(subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyT_T1, isAcute))), by='FASTASampleCode' )	
	tmp					<- merge( data.table(FASTASampleCode=Y.brl[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	Y.brl				<- merge( Y.brl, tmp, by='t.FASTASampleCode')
	stopifnot( Y.brl[, all(isAcute%in%c('Yes','Maybe'))] )
	#
	#	calculate d_TSeqT
	#
	Y.brl				<- merge( Y.brl, data.table(isAcute=c('Yes', 'Maybe'), meanInfWindow= dur.Acute[c('Yes','Maybe')] ), by='isAcute' )		
	Y.brl[, meanInfT:= PosSeqT-meanInfWindow/365.25]
	Y.brl[, DUMMY:=seq_len(nrow(Y.brl))]
	#	calculate time spent suppressed for transmitter
	df.TS		<- subset(Y.brl, select=c(DUMMY, t.Patient, t.AnyT_T1, meanInfT, t.PosSeqT))
	setnames(df.TS, colnames(df.TS), gsub('t.', '', colnames(df.TS), fixed=TRUE))	
	df.TS		<- project.athena.Fisheretal.Y.TSupp(df.TS, df.viro, t.period=t.period, lRNA.supp=lRNA.supp)
	setnames(df.TS, 'tsupp', 't.tsupp')
	Y.brl		<- merge(Y.brl, df.TS, by='DUMMY', all.x=TRUE)
	set(Y.brl, Y.brl[, which(is.na(t.tsupp))], 't.tsupp', 0.)
	#	calculate time spent suppressed for infected
	df.TS		<- subset(Y.brl, select=c(DUMMY, Patient, AnyT_T1, meanInfT, PosSeqT))	
	df.TS		<- project.athena.Fisheretal.Y.TSupp(df.TS, df.viro, t.period=t.period, lRNA.supp=lRNA.supp)
	Y.brl		<- merge(Y.brl, df.TS, by='DUMMY', all.x=TRUE)
	set(Y.brl, Y.brl[, which(is.na(tsupp))], 'tsupp', 0.)
	cat(paste('\nFound person-years of suppressed evolution, ', Y.brl[, sum(t.tsupp)+sum(tsupp)]))
	#	calculate time evolving in transmitter and  infected
	Y.brl[, t.d_TSeqT:=  abs(t.PosSeqT-meanInfT)-t.tsupp ]
	set(Y.brl, Y.brl[, which(t.d_TSeqT<0)], 't.d_TSeqT', 0.)	
	Y.brl[, i.d_TSeqT:=  abs(PosSeqT-meanInfT)-tsupp ]
	set(Y.brl, Y.brl[, which(i.d_TSeqT<0)], 'i.d_TSeqT', 0.)	
	Y.brl[, d_TSeqT:=t.d_TSeqT+i.d_TSeqT]
	cat(paste('\nFound person-years of unsuppressed evolution, ', Y.brl[, sum(d_TSeqT)]))
	#
	#	load Gamma model and compute likelihood
	#
	load(infile.trm.model)		#expect "trm.pol.NO" "trm.pol.nA"   
	Y.brl[, mu:= predict(trm.pol.NO, data=trm.pol.nA, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='mu', type='link')]
	Y.brl[, sigma:= predict(trm.pol.NO, data=trm.pol.nA, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='sigma', type='link')]	
	Y.brl[, score.brl.TPd:= Y.brl[, dNO(brl, mu=mu, sigma=sigma)/pNO(0,mu=mu,sigma=sigma, lower.tail=FALSE)]]
	Y.brl[, score.brl.TPp:= Y.brl[, score.brl.TPd]]
	#
	#	set branch length quantiles
	#	
	Y.brl[, score.brl.TN:= p.rawbrl.unlinked(brl)]
	Y.brl[, score.brl.TP:= p.rawbrl.linked(brl)]
	#	
	if(!is.na(plot.file.score))
	{
		#ggplot(Y.brl, aes(x=score.brl.TPd)) + geom_histogram()		
	}
	#
	Y.brl		<- subset( Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, score.brl.TPd, score.brl.TPp, score.brl.TN, score.brl.TP, d_TSeqT, brl) )
	setnames(Y.brl, 'd_TSeqT','telapsed')
	cat(paste('\nReturn brl score for #t.Patient seqs=',Y.brl[, length(unique(t.FASTASampleCode))]))
	Y.brl		
}
######################################################################################
project.athena.Fisheretal.Y.brlweight.3l<- function(infile.trm.model, Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.viro, dur.Acute=NULL, t.period=0.25, lRNA.supp=log10(51), plot.file.score=NA)
{
	require(reshape2)
	require(ggplot2)
	require(gamlss)		
	#
	#	compute cdf for pairs given they are unlinked
	#
	setkey(Y.rawbrl.unlinked, brl)
	p.rawbrl.unlinked	<- Y.rawbrl.unlinked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	setkey(Y.rawbrl.linked, brl)
	p.rawbrl.linked		<- Y.rawbrl.linked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	#
	#	add sequence sampling times	
	#
	Y.brl				<- copy(Y.rawbrl)
	Y.brl				<- merge( Y.brl, unique(subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyT_T1, isAcute))), by='FASTASampleCode' )	
	tmp					<- merge( data.table(FASTASampleCode=Y.brl[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	Y.brl				<- merge( Y.brl, tmp, by='t.FASTASampleCode')
	stopifnot( Y.brl[, all(isAcute%in%c('Yes','Maybe'))] )
	#
	#	calculate d_TSeqT
	#
	Y.brl				<- merge( Y.brl, data.table(isAcute=c('Yes', 'Maybe'), meanInfWindow= dur.Acute[c('Yes','Maybe')] ), by='isAcute' )		
	Y.brl[, meanInfT:= PosSeqT-meanInfWindow/365.25]
	Y.brl[, DUMMY:=seq_len(nrow(Y.brl))]
	#	calculate time spent suppressed for transmitter
	df.TS		<- subset(Y.brl, select=c(DUMMY, t.Patient, t.AnyT_T1, meanInfT, t.PosSeqT))
	setnames(df.TS, colnames(df.TS), gsub('t.', '', colnames(df.TS), fixed=TRUE))	
	df.TS		<- project.athena.Fisheretal.Y.TSupp(df.TS, df.viro, t.period=t.period, lRNA.supp=lRNA.supp)
	setnames(df.TS, 'tsupp', 't.tsupp')
	Y.brl		<- merge(Y.brl, df.TS, by='DUMMY', all.x=TRUE)
	set(Y.brl, Y.brl[, which(is.na(t.tsupp))], 't.tsupp', 0.)
	#	calculate time spent suppressed for infected
	df.TS		<- subset(Y.brl, select=c(DUMMY, Patient, AnyT_T1, meanInfT, PosSeqT))	
	df.TS		<- project.athena.Fisheretal.Y.TSupp(df.TS, df.viro, t.period=t.period, lRNA.supp=lRNA.supp)
	Y.brl		<- merge(Y.brl, df.TS, by='DUMMY', all.x=TRUE)
	set(Y.brl, Y.brl[, which(is.na(tsupp))], 'tsupp', 0.)
	cat(paste('\nFound person-years of suppressed evolution, ', Y.brl[, sum(t.tsupp)+sum(tsupp)]))
	#	calculate time evolving in transmitter and  infected
	Y.brl[, t.d_TSeqT:=  abs(t.PosSeqT-meanInfT)-t.tsupp ]
	set(Y.brl, Y.brl[, which(t.d_TSeqT<0)], 't.d_TSeqT', 0.)	
	Y.brl[, i.d_TSeqT:=  abs(PosSeqT-meanInfT)-tsupp ]
	set(Y.brl, Y.brl[, which(i.d_TSeqT<0)], 'i.d_TSeqT', 0.)	
	Y.brl[, d_TSeqT:=t.d_TSeqT+i.d_TSeqT]
	cat(paste('\nFound person-years of unsuppressed evolution, ', Y.brl[, sum(d_TSeqT)]))
	#
	#	load Gamma model and compute likelihood
	#
	load(infile.trm.model)		#expect "trm.pol.GA" "trm.pol"   
	Y.brl[, mu:= predict(trm.pol.GA, data=trm.pol, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='mu', type='link')]
	Y.brl[, sigma:= exp( predict(trm.pol.GA, data=trm.pol, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='sigma', type='link') )]	
	Y.brl[, score.brl.TPd:= Y.brl[, dGA(brl, mu=mu, sigma=sigma)]]
	Y.brl[, score.brl.TPp:= Y.brl[, score.brl.TPd]]
	#
	#	set branch length quantiles
	#	
	Y.brl[, score.brl.TN:= p.rawbrl.unlinked(brl)]
	Y.brl[, score.brl.TP:= p.rawbrl.linked(brl)]
	#	
	if(!is.na(plot.file.score))
	{
		#ggplot(Y.brl, aes(x=score.brl.TPd)) + geom_histogram()		
	}
	#
	Y.brl		<- subset( Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, score.brl.TPd, score.brl.TPp, score.brl.TN, score.brl.TP, brl) )
	cat(paste('\nReturn brl score for #t.Patient seqs=',Y.brl[, length(unique(t.FASTASampleCode))]))
	Y.brl		
}
######################################################################################
project.athena.Fisheretal.Y.brlweight.3l.v1<- function(infile.trm.model, Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, dur.Acute=NULL, plot.file.score=NA)
{
	require(reshape2)
	require(ggplot2)
	require(gamlss)		
	#
	#	compute cdf for pairs given they are unlinked
	#
	setkey(Y.rawbrl.unlinked, brl)
	p.rawbrl.unlinked	<- Y.rawbrl.unlinked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	setkey(Y.rawbrl.linked, brl)
	p.rawbrl.linked		<- Y.rawbrl.linked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	#
	#	add sequence sampling times	
	#
	Y.brl				<- copy(Y.rawbrl)
	Y.brl				<- merge( Y.brl, unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1, isAcute))), by='FASTASampleCode' )	
	tmp					<- merge( data.table(FASTASampleCode=Y.brl[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	Y.brl				<- merge( Y.brl, tmp, by='t.FASTASampleCode')
	stopifnot( Y.brl[, all(isAcute%in%c('Yes','Maybe'))] )
	#
	#	calculate d_TSeqT
	#
	Y.brl				<- merge( Y.brl, data.table(isAcute=c('Yes', 'Maybe'), meanInfWindow= dur.Acute[c('Yes','Maybe')] ), by='isAcute' )		
	Y.brl[, d_TSeqT:= meanInfWindow/365.25 + abs(t.PosSeqT-(PosSeqT-meanInfWindow/365.25))]
	#
	#	load Gamma model and compute likelihood
	#
	load(infile.trm.model)		#expect "trm.pol.GA" "trm.pol"   
	Y.brl[, mu:= predict(trm.pol.GA, data=trm.pol, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='mu', type='link')]
	Y.brl[, sigma:= exp( predict(trm.pol.GA, data=trm.pol, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='sigma', type='link') )]	
	Y.brl[, score.brl.TPd:= Y.brl[, dGA(brl, mu=mu, sigma=sigma)]]
	Y.brl[, score.brl.TPp:= Y.brl[, score.brl.TPd]]
	#
	#	set branch length quantiles
	#	
	Y.brl[, score.brl.TN:= p.rawbrl.unlinked(brl)]
	Y.brl[, score.brl.TP:= p.rawbrl.linked(brl)]
	#	
	if(!is.na(plot.file.score))
	{
		#ggplot(Y.brl, aes(x=score.brl.TPd)) + geom_histogram()		
	}
	#
	Y.brl		<- subset( Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, score.brl.TPd, score.brl.TPp, score.brl.TN, score.brl.TP, brl) )
	cat(paste('\nReturn brl score for #t.Patient seqs=',Y.brl[, length(unique(t.FASTASampleCode))]))
	Y.brl		
}
######################################################################################
project.athena.Fisheretal.Y.brlweight<- function(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, df.treatment, brl.linked.max.brlr= 0.04, brl.linked.min.brl= 1e-4, brl.linked.max.dt= 10, brl.linked.min.dt= 1.5, brl.bwhost.multiplier=1, plot.file.score=NA, method='3c', plot.file.both=NA, plot.file.one=NA, save.file.3ha=NA)
{
	#brl.linked.max.brlr= 0.01; brl.linked.min.brl= 1e-12; brl.linked.max.dt= 10; brl.linked.min.dt= 1; plot.file.score=NA; method='3aa'; plot.file.both=NA; plot.file.one=NA
	#	check if Exp model would be reasonable	
	require(MASS)	
	require(reshape2)
	require(ggplot2)
	require(gamlss)	
	require(fitdistrplus)
	require(ADGofTest)
	require(gamlss.mx)
	library(distr)
	#	compute cdf for pairs given they are unlinked
	setkey(Y.rawbrl.unlinked, brl)
	p.rawbrl.unlinked	<- Y.rawbrl.unlinked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	setkey(Y.rawbrl.linked, brl)
	p.rawbrl.linked		<- Y.rawbrl.linked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]
	
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
	if(class(df.all$PosSeqT)=='Date')
		set(Y.rawbrl.linked, NULL, 'PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,PosSeqT]))
	if(class(df.all$t.PosSeqT)=='Date')
		set(Y.rawbrl.linked, NULL, 't.PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,t.PosSeqT]))
	if(class(df.all$AnyT_T1)=='Date')
		set(Y.rawbrl.linked, NULL, 'AnyT_T1', hivc.db.Date2numeric(Y.rawbrl.linked[,AnyT_T1]))
	if(class(df.all$t.AnyT_T1)=='Date')
		set(Y.rawbrl.linked, NULL, 't.AnyT_T1', hivc.db.Date2numeric(Y.rawbrl.linked[,t.AnyT_T1]))
	#	compute b4T
	Y.rawbrl.linked[, b4T:= 'both']
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT<=AnyT_T1)], 'b4T', 'only.RI')
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT<=t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'only.T')
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'none')
	set(Y.rawbrl.linked, NULL, 'b4T', Y.rawbrl.linked[, factor(b4T)])
	tmp				<- c('both sampling dates\nafter cART initiation','both sampling dates\nbefore cART initiation','sampling dates before and\nafter cART initiation')
	Y.rawbrl.linked	<- merge( Y.rawbrl.linked, data.table(b4T= c('none','both','only.T'), b4T.long=tmp), by='b4T' )
	set(Y.rawbrl.linked, NULL, 'b4T.long', Y.rawbrl.linked[, factor(b4T.long)])	
	#	compute ART.C
	tmp				<- subset(df.treatment, select=c(Patient, StopTime, TrI, TrCh.failure, TrVL.failure))
	set(tmp, NULL, 'StopTime', hivc.db.Date2numeric(tmp[,StopTime]))
	tmp				<- merge(tmp, subset( Y.rawbrl.linked, select=c(Patient, PosSeqT) ), by='Patient', allow.cartesian=TRUE)
	tmp				<- subset(tmp, StopTime<=PosSeqT)	
	tmp				<- tmp[, list( ART.C=any( sapply(.SD, function(z) any(z=='Yes', na.rm=TRUE)) ) ), by=c('Patient','PosSeqT'), .SDcols=c('TrI','TrCh.failure','TrVL.failure')]
	Y.rawbrl.linked	<- merge( Y.rawbrl.linked, tmp, by=c('Patient','PosSeqT'), all.x=TRUE )		
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(is.na(ART.C))], 'ART.C', FALSE)
	set(Y.rawbrl.linked, NULL, 'ART.C', Y.rawbrl.linked[, factor(ART.C, levels=c(FALSE, TRUE), labels=c('No','Yes'))] )
	#	t.ART.C does not add anything of course
	#	because we consider at least one of the two sampling time after complication, and time are ordered
	#tmp				<- subset(df.treatment, select=c(Patient, StopTime, TrI, TrCh.failure, TrVL.failure))
	#set(tmp, NULL, 'StopTime', hivc.db.Date2numeric(tmp[,StopTime]))	
	#tmp				<- merge(tmp, subset( Y.rawbrl.linked, select=c(Patient, t.PosSeqT) ), by='Patient', allow.cartesian=TRUE)
	#tmp				<- subset(tmp, StopTime<=t.PosSeqT)		
	#tmp				<- tmp[, list( t.ART.C=any( sapply(.SD, function(z) any(z=='Yes', na.rm=TRUE)) ) ), by=c('Patient','t.PosSeqT'), .SDcols=c('TrI','TrCh.failure','TrVL.failure')]
	#Y.rawbrl.linked	<- merge( Y.rawbrl.linked, tmp, by=c('Patient','t.PosSeqT'), all.x=TRUE )		
	#set(Y.rawbrl.linked, Y.rawbrl.linked[, which(is.na(t.ART.C))], 't.ART.C', FALSE)
	#set(Y.rawbrl.linked, NULL, 't.ART.C', Y.rawbrl.linked[, factor(t.ART.C, levels=c(FALSE, TRUE), labels=c('No','Yes'))] )
	tmp				<- c('Both sampling dates\nbefore interruption or failure','At least one sampling date\nafter interruption or failure')
	Y.rawbrl.linked	<- merge( Y.rawbrl.linked, data.table(ART.C= c('No','Yes'), ART.C.long=tmp), by='ART.C' )
	set(Y.rawbrl.linked, NULL, 'ART.C.long', Y.rawbrl.linked[, factor(ART.C.long, levels=tmp)])	
	#
	Y.rawbrl.linked[, dt:= abs(PosSeqT-t.PosSeqT)]		
	Y.rawbrl.linked[, brlr:= brl/dt]
	#
	#	EXCLUDE options
	#
	Y.rawbrl.linked		<- subset( Y.rawbrl.linked, dt!=0 )
	#	brl.linked.max.brlr<- 0.01; brl.linked.min.brl<- 1e-12; brl.linked.min.dt<- -1
	cat(paste('\nY.rawbrl.linked all, n=',nrow(Y.rawbrl.linked)))
	Y.rawbrl.linked		<- subset(Y.rawbrl.linked,  brlr<=brl.linked.max.brlr)
	cat(paste('\nY.rawbrl.linked max brlr ',brl.linked.max.brlr,' excluded, n=',nrow(Y.rawbrl.linked)))
	Y.rawbrl.linked		<- subset(Y.rawbrl.linked,  dt<=brl.linked.max.dt)		
	cat(paste('\nY.rawbrl.linked max dt ',brl.linked.max.dt,' excluded, n=',nrow(Y.rawbrl.linked)))	
	Y.rawbrl.linked		<- subset(Y.rawbrl.linked,  brlr>brl.linked.min.brl)
	cat(paste('\nY.rawbrl.linked min brl ',brl.linked.min.brl,' excluded, n=',nrow(Y.rawbrl.linked)))
	Y.rawbrl.linked		<- subset(Y.rawbrl.linked,  dt>brl.linked.min.dt)		
	cat(paste('\nY.rawbrl.linked min dt ',brl.linked.min.dt,' excluded, n=',nrow(Y.rawbrl.linked)))
	#
	#	fit distributions to data 
	#
	#
	#	zero adjusted Gamma 
	#
	Y.rawbrl.linked[, brlz:=brl]
	tmp					<- Y.rawbrl.linked[, which(brlz<1e-4)]
	set(Y.rawbrl.linked, tmp, 'brlz', 0)
	brl					<- Y.rawbrl.linked[,brlz]
	ml.zaga.all			<- gamlss(brl~1, family=ZAGA)	
	#	ART.C
	brl					<- subset(Y.rawbrl.linked, ART.C=='No'  )[,brlz]
	ml.zaga.CN			<- gamlss(brl~1, family=ZAGA)
	brl					<- subset(Y.rawbrl.linked, ART.C=='Yes'  )[,brlz]
	ml.zaga.CY			<- gamlss(brl~1, family=ZAGA)		
	#	B4T
	brl					<- subset(Y.rawbrl.linked, b4T=='both'  )[,brlz]
	ml.zaga				<- gamlss(brl~1, family=ZAGA)
	ml.zaga.p			<- ad.test(brl, pZAGA, mu=exp(ml.zaga[['mu.coefficients']]), sigma=exp(ml.zaga[['sigma.coefficients']]), nu=1/(1+exp(-ml.zaga[['nu.coefficients']])) )$p.value
	cat(paste('\nbest ZAGA pvalue for both b4T=', ml.zaga.p))
	brl					<- subset(Y.rawbrl.linked, b4T!='both' & b4T!='none')[,brlz]
	ml.zaga.one			<- gamlss(brl~1, family=ZAGA)
	ml.zaga.one.p		<- ad.test(brl, pZAGA, mu=exp(ml.zaga.one[['mu.coefficients']]), sigma=exp(ml.zaga.one[['sigma.coefficients']]), nu=1/(1+exp(-ml.zaga.one[['nu.coefficients']])) )$p.value
	cat(paste('\nbest ZAGA pvalue for one b4T=', ml.zaga.one.p))	
	ml.zaga.dt.data		<- as.data.frame(subset(Y.rawbrl.linked, b4T=='both', select=c(brlz, dt)  ))		
	#ml.zaga.dt			<- gamlss(brlz~dt, family=ZAGA, data=ml.zaga.dt.data)
	ml.zaga.dt			<- gamlss(brlz~dt, sigma.formula= ~dt, nu.formula= ~dt, family=ZAGA, data=ml.zaga.dt.data)		
	ml.zaga.dt.data.one	<- as.data.frame(subset(Y.rawbrl.linked, b4T=='only.T', select=c(brlz, dt)  ))		
	ml.zaga.dt.one		<- gamlss(brlz~dt, sigma.formula= ~dt, nu.formula= ~dt, family=ZAGA, data=ml.zaga.dt.data.one)			
	#
	#	fit Gamma for b4T==both and b4T one in ART
	#
	Y.brl			<- copy(Y.rawbrl)	
	Y.brl[, score.brl.TN:= p.rawbrl.unlinked(brl)]
	Y.brl[, score.brl.TP:= p.rawbrl.linked(brl)]
	if(substr(method,1,2)%in%c('3c','3e'))	#ignore zeros use Gamma
	{
		cat(paste('\nBoth b4T: MLE param:',paste( mle.ga$estimate, collapse=' ')))
		cat(paste('\nOne b4T: MLE param:',paste( mle.ga.one$estimate, collapse=' ')))
		
		Y.brl			<- merge( Y.brl, unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )	
		tmp				<- merge( data.table(FASTASampleCode=Y.brl[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		Y.brl			<- merge( Y.brl, tmp, by='t.FASTASampleCode')
		Y.brl[, b4T:= 'both']
		set(Y.brl, Y.brl[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT<=AnyT_T1)], 'b4T', 'one')
		set(Y.brl, Y.brl[, which(t.PosSeqT<=t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'one')
		set(Y.brl, Y.brl[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'one')		
		print(Y.brl[,table(b4T)])
		
		Y.brl[, score.brl.TPp:= pgamma(Y.brl[,brl*brl.bwhost.multiplier], shape=2*mle.ga$estimate['shape'], rate=mle.ga$estimate['rate'], lower.tail=FALSE)]
		tmp<- Y.brl[, which(b4T=='one')]
		set(Y.brl, tmp, 'score.brl.TPp', pgamma(Y.brl[tmp,brl*brl.bwhost.multiplier], shape=2*mle.ga.one$estimate['shape'], rate=mle.ga.one$estimate['rate'], lower.tail=FALSE) )
		Y.brl[, score.brl.TPd:= score.brl.TPp]		
	}
	if(substr(method,1,2)%in%c('3d','3f'))	#dont ignore zeros - use zero adjusted Gamma
	{
		tmp			<- c( exp(ml.zaga[['mu.coefficients']]), exp(ml.zaga[['sigma.coefficients']]), 1/(1+exp(-ml.zaga[['nu.coefficients']])) )
		cat(paste('\nZAGA Both b4T: MLE param:',paste( tmp, collapse=' ')))
		tmp			<- c( exp(ml.zaga.one[['mu.coefficients']]), exp(ml.zaga.one[['sigma.coefficients']]), 1/(1+exp(-ml.zaga.one[['nu.coefficients']])) )
		cat(paste('\nZAGA One b4T: MLE param:',paste( tmp, collapse=' ')))
		
		Y.brl[, brlz:=brl]		 
		set(Y.brl, Y.brl[, which(brlz<1e-4)], 'brlz', 0)
		Y.brl			<- merge( Y.brl, unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )	
		tmp				<- merge( data.table(FASTASampleCode=Y.brl[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		Y.brl			<- merge( Y.brl, tmp, by='t.FASTASampleCode')
		Y.brl[, b4T:= 'both']
		set(Y.brl, Y.brl[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT<=AnyT_T1)], 'b4T', 'one')
		set(Y.brl, Y.brl[, which(t.PosSeqT<=t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'one')
		set(Y.brl, Y.brl[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'one')		
		print(Y.brl[,table(b4T)])
		
		#	convolution of zero adjusted Gamma	
		ml.zaga.pa		<- c(mu=as.double(exp(ml.zaga[['mu.coefficients']])), sigma=as.double(exp(ml.zaga[['sigma.coefficients']])), nu=as.double(1/(1+exp(-ml.zaga[['nu.coefficients']]))) )
		ml.zaga.pa		<- c(ml.zaga.pa, shape= as.double(1/(ml.zaga.pa['sigma']*ml.zaga.pa['sigma'])), scale= as.double(ml.zaga.pa['mu']*ml.zaga.pa['sigma']*ml.zaga.pa['sigma']) )
		ml.zaga.one.pa	<- c(mu=as.double(exp(ml.zaga.one[['mu.coefficients']])), sigma=as.double(exp(ml.zaga.one[['sigma.coefficients']])), nu=as.double(1/(1+exp(-ml.zaga.one[['nu.coefficients']]))) )
		ml.zaga.one.pa	<- c(ml.zaga.one.pa, shape= as.double(1/(ml.zaga.one.pa['sigma']*ml.zaga.one.pa['sigma'])), scale= as.double(ml.zaga.one.pa['mu']*ml.zaga.one.pa['sigma']*ml.zaga.one.pa['sigma']) )
		#	sense check -- should evaluate to 0, ~1
		tmp				<- c(0, 0.1)
		tmp				<- 2*ml.zaga.pa['nu']*(1-ml.zaga.pa['nu'])*pGA( tmp,  mu=ml.zaga.pa['mu'], sigma=ml.zaga.pa['sigma'] ) + (1-ml.zaga.pa['nu'])*(1-ml.zaga.pa['nu'])*pgamma( tmp,  shape=2*ml.zaga.pa['shape'], scale=ml.zaga.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.pa['nu']*ml.zaga.pa['nu'])				
		cat(paste('\nZAGA Both b4T convolution sense check: should be ~0 ~1',paste( tmp, collapse=' ')))
		#
		tmp				<- 2*ml.zaga.pa['nu']*(1-ml.zaga.pa['nu'])*pGA( Y.brl[,brlz*brl.bwhost.multiplier],  mu=ml.zaga.pa['mu'], sigma=ml.zaga.pa['sigma'] ) + (1-ml.zaga.pa['nu'])*(1-ml.zaga.pa['nu'])*pgamma( Y.brl[,brlz*brl.bwhost.multiplier],  shape=2*ml.zaga.pa['shape'], scale=ml.zaga.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.pa['nu']*ml.zaga.pa['nu'])				
		Y.brl[, score.brl.TPp:= 1-tmp]
		if(substr(method,1,2)=='3d')
		{
			tmp				<- Y.brl[, which(b4T=='one')]
			tmp				<- 2*ml.zaga.one.pa['nu']*(1-ml.zaga.one.pa['nu'])*pGA( Y.brl[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.one.pa['mu'], sigma=ml.zaga.one.pa['sigma'] ) + (1-ml.zaga.one.pa['nu'])*(1-ml.zaga.one.pa['nu'])*pgamma( Y.brl[tmp,brlz*brl.bwhost.multiplier],  shape=2*ml.zaga.one.pa['shape'], scale=ml.zaga.one.pa['scale'] ) 
			tmp				<- tmp / (1-ml.zaga.one.pa['nu']*ml.zaga.one.pa['nu'])				
			set(Y.brl, Y.brl[, which(b4T=='one')], 'score.brl.TPp', 1-tmp )			
		}		
		Y.brl[, score.brl.TPd:= score.brl.TPp]		
		Y.brl[, dt:= abs(PosSeqT-t.PosSeqT)]		
	}
	if(substr(method,1,2)%in%c('3i'))	#dont ignore zeros - use zero adjusted Gamma
	{
		tmp			<- c( exp(ml.zaga.all[['mu.coefficients']]), exp(ml.zaga.all[['sigma.coefficients']]), 1/(1+exp(-ml.zaga.all[['nu.coefficients']])) )
		cat(paste('\nZAGA.all: MLE param:',paste( tmp, collapse=' ')))
		#ZAGA.all: MLE param: 0.0258211573917653 0.771293214006904 0.177286356821589b4T
		Y.brl[, brlz:=brl]		 
		set(Y.brl, Y.brl[, which(brlz<1e-4)], 'brlz', 0)		
		Y.brl			<- merge( Y.brl, unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )	
		tmp				<- merge( data.table(FASTASampleCode=Y.brl[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		Y.brl			<- merge( Y.brl, tmp, by='t.FASTASampleCode')
		Y.brl[, b4T:= 'all']
		print(Y.brl[,table(b4T)])
		#	convolution of zero adjusted Gamma	
		ml.zaga.all.pa	<- c(mu=as.double(exp(ml.zaga.all[['mu.coefficients']])), sigma=as.double(exp(ml.zaga.all[['sigma.coefficients']])), nu=as.double(1/(1+exp(-ml.zaga.all[['nu.coefficients']]))) )
		ml.zaga.all.pa	<- c(ml.zaga.all.pa, shape= as.double(1/(ml.zaga.all.pa['sigma']*ml.zaga.all.pa['sigma'])), scale= as.double(ml.zaga.all.pa['mu']*ml.zaga.all.pa['sigma']*ml.zaga.all.pa['sigma']) )
		#	sense check -- should evaluate to 0, ~1
		tmp				<- c(0, 0.1)
		tmp				<- 2*ml.zaga.all.pa['nu']*(1-ml.zaga.all.pa['nu'])*pGA( tmp,  mu=ml.zaga.all.pa['mu'], sigma=ml.zaga.all.pa['sigma'] ) + (1-ml.zaga.all.pa['nu'])*(1-ml.zaga.all.pa['nu'])*pgamma( tmp,  shape=2*ml.zaga.all.pa['shape'], scale=ml.zaga.all.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.all.pa['nu']*ml.zaga.all.pa['nu'])				
		cat(paste('\nZAGA Both b4T convolution sense check: should be ~0 ~1',paste( tmp, collapse=' ')))		
		#	set score for cases 'both' and 'one'
		tmp				<- 2*ml.zaga.all.pa['nu']*(1-ml.zaga.all.pa['nu'])*pGA( Y.brl[,brlz*brl.bwhost.multiplier],  mu=ml.zaga.all.pa['mu'], sigma=ml.zaga.all.pa['sigma'] ) + (1-ml.zaga.all.pa['nu'])*(1-ml.zaga.all.pa['nu'])*pgamma( Y.brl[,brlz*brl.bwhost.multiplier],  shape=2*ml.zaga.all.pa['shape'], scale=ml.zaga.all.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.all.pa['nu']*ml.zaga.all.pa['nu'])
		set(Y.brl, NULL, 'score.brl.TPp', 1-tmp )		
		#
		Y.brl[, score.brl.TPd:= score.brl.TPp]		
		Y.brl[, dt:= abs(PosSeqT-t.PosSeqT)]		
	}
	if(substr(method,1,2)%in%c('3j'))	#dont ignore zeros - use zero adjusted Gamma
	{
		tmp			<- c( exp(ml.zaga[['mu.coefficients']]), exp(ml.zaga[['sigma.coefficients']]), 1/(1+exp(-ml.zaga[['nu.coefficients']])) )
		cat(paste('\nZAGA Both b4T: MLE param:',paste( tmp, collapse=' ')))
		tmp			<- c( exp(ml.zaga.one[['mu.coefficients']]), exp(ml.zaga.one[['sigma.coefficients']]), 1/(1+exp(-ml.zaga.one[['nu.coefficients']])) )
		cat(paste('\nZAGA One b4T: MLE param:',paste( tmp, collapse=' ')))
		
		Y.brl[, brlz:=brl]		 
		set(Y.brl, Y.brl[, which(brlz<1e-4)], 'brlz', 0)
		Y.brl			<- merge( Y.brl, unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )	
		tmp				<- merge( data.table(FASTASampleCode=Y.brl[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		Y.brl			<- merge( Y.brl, tmp, by='t.FASTASampleCode')
		Y.brl[, b4T:= 'both']
		set(Y.brl, Y.brl[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT<=AnyT_T1)], 'b4T', 'mix')
		set(Y.brl, Y.brl[, which(t.PosSeqT<=t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'mix')
		set(Y.brl, Y.brl[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'one')		
		print(Y.brl[,table(b4T)])
		#
		#	cases 'both' and 'one'
		#
		#	convolution of zero adjusted Gamma	
		ml.zaga.pa		<- c(mu=as.double(exp(ml.zaga[['mu.coefficients']])), sigma=as.double(exp(ml.zaga[['sigma.coefficients']])), nu=as.double(1/(1+exp(-ml.zaga[['nu.coefficients']]))) )
		ml.zaga.pa		<- c(ml.zaga.pa, shape= as.double(1/(ml.zaga.pa['sigma']*ml.zaga.pa['sigma'])), scale= as.double(ml.zaga.pa['mu']*ml.zaga.pa['sigma']*ml.zaga.pa['sigma']) )
		ml.zaga.one.pa	<- c(mu=as.double(exp(ml.zaga.one[['mu.coefficients']])), sigma=as.double(exp(ml.zaga.one[['sigma.coefficients']])), nu=as.double(1/(1+exp(-ml.zaga.one[['nu.coefficients']]))) )
		ml.zaga.one.pa	<- c(ml.zaga.one.pa, shape= as.double(1/(ml.zaga.one.pa['sigma']*ml.zaga.one.pa['sigma'])), scale= as.double(ml.zaga.one.pa['mu']*ml.zaga.one.pa['sigma']*ml.zaga.one.pa['sigma']) )
		#	sense check -- should evaluate to 0, ~1
		tmp				<- c(0, 0.1)
		tmp				<- 2*ml.zaga.pa['nu']*(1-ml.zaga.pa['nu'])*pGA( tmp,  mu=ml.zaga.pa['mu'], sigma=ml.zaga.pa['sigma'] ) + (1-ml.zaga.pa['nu'])*(1-ml.zaga.pa['nu'])*pgamma( tmp,  shape=2*ml.zaga.pa['shape'], scale=ml.zaga.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.pa['nu']*ml.zaga.pa['nu'])				
		cat(paste('\nZAGA Both b4T convolution sense check: should be ~0 ~1',paste( tmp, collapse=' ')))
		Y.brl[, score.brl.TPp:= NA_real_]
		#	set score for cases 'both' and 'one'
		tmp				<- Y.brl[, which(b4T=='both')]
		tmp				<- 2*ml.zaga.pa['nu']*(1-ml.zaga.pa['nu'])*pGA( Y.brl[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.pa['mu'], sigma=ml.zaga.pa['sigma'] ) + (1-ml.zaga.pa['nu'])*(1-ml.zaga.pa['nu'])*pgamma( Y.brl[tmp,brlz*brl.bwhost.multiplier],  shape=2*ml.zaga.pa['shape'], scale=ml.zaga.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.pa['nu']*ml.zaga.pa['nu'])
		set(Y.brl, Y.brl[, which(b4T=='both')], 'score.brl.TPp', 1-tmp )
		tmp				<- Y.brl[, which(b4T=='one')]
		tmp				<- 2*ml.zaga.one.pa['nu']*(1-ml.zaga.one.pa['nu'])*pGA( Y.brl[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.one.pa['mu'], sigma=ml.zaga.one.pa['sigma'] ) + (1-ml.zaga.one.pa['nu'])*(1-ml.zaga.one.pa['nu'])*pgamma( Y.brl[tmp,brlz*brl.bwhost.multiplier],  shape=2*ml.zaga.one.pa['shape'], scale=ml.zaga.one.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.one.pa['nu']*ml.zaga.one.pa['nu'])				
		set(Y.brl, Y.brl[, which(b4T=='one')], 'score.brl.TPp', 1-tmp )			
		#
		#	case 'mix'
		#		
		#	get cdf of convolution GA_noART with GA_ART
		print('h')
		ml.zaga.pmix	<- p(convpow( Gammad(scale=ml.zaga.pa['scale'], shape=ml.zaga.pa['shape']) + Gammad(scale=ml.zaga.one.pa['scale'], shape=ml.zaga.one.pa['shape']), 1))
		print('h2')
		tmp				<- Y.brl[, which(b4T=='mix')]
		print('h3')
		#	get cdf of convolution ZAGA_mix (not normalized)
		tmp				<- ml.zaga.pa['nu']*(1-ml.zaga.pa['nu'])*pGA( Y.brl[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.pa['mu'], sigma=ml.zaga.pa['sigma'] ) +
							ml.zaga.one.pa['nu']*(1-ml.zaga.one.pa['nu'])*pGA( Y.brl[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.one.pa['mu'], sigma=ml.zaga.one.pa['sigma'] ) +
							(1-ml.zaga.pa['nu'])*(1-ml.zaga.one.pa['nu'])*ml.zaga.pmix( Y.brl[tmp,brlz*brl.bwhost.multiplier] )
		print('h4')			
		#	divide by normalizing constant
		tmp				<- tmp / ( ml.zaga.pa['nu']*(1-ml.zaga.pa['nu']) + ml.zaga.one.pa['nu']*(1-ml.zaga.one.pa['nu']) + (1-ml.zaga.pa['nu'])*(1-ml.zaga.one.pa['nu']) )
		set(Y.brl, Y.brl[, which(b4T=='mix')], 'score.brl.TPp', 1-tmp )
		#
		Y.brl[, score.brl.TPd:= score.brl.TPp]		
		Y.brl[, dt:= abs(PosSeqT-t.PosSeqT)]		
	}
	if(substr(method,1,2)%in%c('3k'))	#dont ignore zeros - use zero adjusted Gamma
	{
		cat('\nin method 3k')
		ml.zaga.CN.pa	<- c(mu=as.double(exp(ml.zaga.CN[['mu.coefficients']])), sigma=as.double(exp(ml.zaga.CN[['sigma.coefficients']])), nu=as.double(1/(1+exp(-ml.zaga.CN[['nu.coefficients']]))) )
		ml.zaga.CN.pa	<- c(ml.zaga.CN.pa, shape= as.double(1/(ml.zaga.CN.pa['sigma']*ml.zaga.CN.pa['sigma'])), scale= as.double(ml.zaga.CN.pa['mu']*ml.zaga.CN.pa['sigma']*ml.zaga.CN.pa['sigma']) )
		ml.zaga.CY.pa	<- c(mu=as.double(exp(ml.zaga.CY[['mu.coefficients']])), sigma=as.double(exp(ml.zaga.CY[['sigma.coefficients']])), nu=as.double(1/(1+exp(-ml.zaga.CY[['nu.coefficients']]))) )
		ml.zaga.CY.pa	<- c(ml.zaga.CY.pa, shape= as.double(1/(ml.zaga.CY.pa['sigma']*ml.zaga.CY.pa['sigma'])), scale= as.double(ml.zaga.CY.pa['mu']*ml.zaga.CY.pa['sigma']*ml.zaga.CY.pa['sigma']) )
		cat(paste('\nZAGA.CN: MLE param:',paste( ml.zaga.CN.pa, collapse=' ')))
		cat(paste('\nZAGA.CY: MLE param:',paste( ml.zaga.CY.pa, collapse=' ')))
		#
		Y.brl[, brlz:=brl]		 
		set(Y.brl, Y.brl[, which(brlz<1e-4)], 'brlz', 0)
		Y.brl			<- merge( Y.brl, unique(subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyT_T1))), by='FASTASampleCode' )	
		tmp				<- merge( data.table(FASTASampleCode=Y.brl[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, Patient, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		Y.brl			<- merge( Y.brl, tmp, by='t.FASTASampleCode')
		#	ART.C for recipient
		tmp				<- subset(df.treatment, select=c(Patient, StopTime, TrI, TrCh.failure, TrVL.failure))
		set(tmp, NULL, 'StopTime', hivc.db.Date2numeric(tmp[,StopTime]))
		tmp				<- merge(tmp, subset( Y.brl, select=c(Patient, PosSeqT) ), by='Patient', allow.cartesian=TRUE)
		tmp				<- subset(tmp, StopTime<=PosSeqT)	
		tmp				<- tmp[, list( ART.C=any( sapply(.SD, function(z) any(z=='Yes', na.rm=TRUE)) ) ), by=c('Patient','PosSeqT'), .SDcols=c('TrI','TrCh.failure','TrVL.failure')]
		Y.brl			<- merge( Y.brl, tmp, by=c('Patient','PosSeqT'), all.x=TRUE )		
		#	ART.C for transmitter
		tmp				<- subset(df.treatment, select=c(Patient, StopTime, TrI, TrCh.failure, TrVL.failure))
		set(tmp, NULL, 'StopTime', hivc.db.Date2numeric(tmp[,StopTime]))
		setnames(tmp, 'Patient','t.Patient')
		tmp				<- merge(tmp, subset( Y.brl, select=c(t.Patient, t.PosSeqT) ), by='t.Patient', allow.cartesian=TRUE)
		tmp				<- subset(tmp, StopTime<=t.PosSeqT)	
		tmp				<- tmp[, list( t.ART.C=any( sapply(.SD, function(z) any(z=='Yes', na.rm=TRUE)) ) ), by=c('t.Patient','t.PosSeqT'), .SDcols=c('TrI','TrCh.failure','TrVL.failure')]
		Y.brl			<- merge( Y.brl, tmp, by=c('t.Patient','t.PosSeqT'), all.x=TRUE )		
		set(Y.brl, Y.brl[, which(is.na(ART.C))], 'ART.C', FALSE)
		set(Y.brl, Y.brl[, which(is.na(t.ART.C))], 't.ART.C', FALSE)
		set(Y.brl, NULL, 'ART.C', Y.brl[, factor(ART.C, levels=c(FALSE, TRUE), labels=c('No','Yes'))] )
		set(Y.brl, NULL, 't.ART.C', Y.brl[, factor(t.ART.C, levels=c(FALSE, TRUE), labels=c('No','Yes'))] )
		print(Y.brl[, table(ART.C, t.ART.C)])
		#       t.ART.C
		#ART.C   FALSE TRUE
  		#FALSE  7505  383
  		#TRUE     78   11
		#
		#	cases F/F and T/T
		#
		#	convolution of zero adjusted Gamma	
		#	sense check -- should evaluate to 0, ~1
		tmp				<- c(0, 0.1)
		tmp				<- 2*ml.zaga.CN.pa['nu']*(1-ml.zaga.CN.pa['nu'])*pGA( tmp,  mu=ml.zaga.CN.pa['mu'], sigma=ml.zaga.CN.pa['sigma'] ) + (1-ml.zaga.CN.pa['nu'])*(1-ml.zaga.CN.pa['nu'])*pgamma( tmp,  shape=2*ml.zaga.CN.pa['shape'], scale=ml.zaga.CN.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.CN.pa['nu']*ml.zaga.CN.pa['nu'])				
		cat(paste('\nZAGA Both b4T convolution sense check: should be ~0 ~1',paste( tmp, collapse=' ')))
		Y.brl[, score.brl.TPp:= NA_real_]
		#	set score for cases F/F
		tmp				<- Y.brl[, which(t.ART.C=='No' & ART.C=='No')]
		tmp				<- 2*ml.zaga.CN.pa['nu']*(1-ml.zaga.CN.pa['nu'])*pGA( Y.brl[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.CN.pa['mu'], sigma=ml.zaga.CN.pa['sigma'] ) + (1-ml.zaga.CN.pa['nu'])*(1-ml.zaga.CN.pa['nu'])*pgamma( Y.brl[tmp,brlz*brl.bwhost.multiplier],  shape=2*ml.zaga.CN.pa['shape'], scale=ml.zaga.CN.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.CN.pa['nu']*ml.zaga.CN.pa['nu'])
		set(Y.brl, Y.brl[, which(t.ART.C=='No' & ART.C=='No')], 'score.brl.TPp', 1-tmp )
		#	set score for cases T/T
		tmp				<- Y.brl[, which(t.ART.C=='Yes' & ART.C=='Yes')]
		tmp				<- 2*ml.zaga.CY.pa['nu']*(1-ml.zaga.CY.pa['nu'])*pGA( Y.brl[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.CY.pa['mu'], sigma=ml.zaga.CY.pa['sigma'] ) + (1-ml.zaga.CY.pa['nu'])*(1-ml.zaga.CY.pa['nu'])*pgamma( Y.brl[tmp,brlz*brl.bwhost.multiplier],  shape=2*ml.zaga.CY.pa['shape'], scale=ml.zaga.CY.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.CY.pa['nu']*ml.zaga.CY.pa['nu'])				
		set(Y.brl, Y.brl[, which(t.ART.C=='Yes' & ART.C=='Yes')], 'score.brl.TPp', 1-tmp )			
		#
		#	case 'mix'
		#		
		#	get cdf of convolution GA_noART with GA_ART
		ml.zaga.pmix	<- p(convpow( Gammad(scale=ml.zaga.CN.pa['scale'], shape=ml.zaga.CN.pa['shape']) + Gammad(scale=ml.zaga.CY.pa['scale'], shape=ml.zaga.CY.pa['shape']), 1))
		ml.zaga.dmix	<- d(convpow( Gammad(scale=ml.zaga.CN.pa['scale'], shape=ml.zaga.CN.pa['shape']) + Gammad(scale=ml.zaga.CY.pa['scale'], shape=ml.zaga.CY.pa['shape']), 1))
		tmp				<- Y.brl[, which(t.ART.C!=ART.C)]
		#	get cdf of convolution ZAGA_mix (not normalized)
		tmp				<- ml.zaga.CN.pa['nu']*(1-ml.zaga.CN.pa['nu'])*pGA( Y.brl[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.CN.pa['mu'], sigma=ml.zaga.CN.pa['sigma'] ) +
							ml.zaga.CY.pa['nu']*(1-ml.zaga.CY.pa['nu'])*pGA( Y.brl[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.CY.pa['mu'], sigma=ml.zaga.CY.pa['sigma'] ) +
							(1-ml.zaga.CN.pa['nu'])*(1-ml.zaga.CY.pa['nu'])*ml.zaga.pmix( Y.brl[tmp,brlz*brl.bwhost.multiplier] )	
		#	divide by normalizing constant
		tmp				<- tmp / ( ml.zaga.CN.pa['nu']*(1-ml.zaga.CN.pa['nu']) + ml.zaga.CY.pa['nu']*(1-ml.zaga.CY.pa['nu']) + (1-ml.zaga.CN.pa['nu'])*(1-ml.zaga.CY.pa['nu']) )
		set(Y.brl, Y.brl[, which(t.ART.C!=ART.C)], 'score.brl.TPp', 1-tmp )
		#
		Y.brl[, score.brl.TPd:= score.brl.TPp]		
		Y.brl[, dt:= abs(PosSeqT-t.PosSeqT)]		
	}
	if(substr(method,1,2)=='3g')	#dont ignore zeros - use zero adjusted Gamma. Ignore after cART initiation. Ignore convolution.
	{		
		Y.brl[, brlz:=brl]		 
		set(Y.brl, Y.brl[, which(brlz<1e-4)], 'brlz', 0)
		Y.brl			<- merge( Y.brl, unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )	
		tmp				<- merge( data.table(FASTASampleCode=Y.brl[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		Y.brl			<- merge( Y.brl, tmp, by='t.FASTASampleCode')
		Y.brl[, dt:= abs(PosSeqT-t.PosSeqT)]				
		Y.brl[, mu:=exp( predict( ml.zaga.dt, data=ml.zaga.dt.data, newdata=as.data.frame(subset(Y.brl, select=c(brlz, dt))), what='mu', type='link' ) )]
		Y.brl[, sigma:=exp( predict( ml.zaga.dt, data=ml.zaga.dt.data, newdata=as.data.frame(subset(Y.brl, select=c(brlz, dt))), what='sigma', type='link' ) )]
		Y.brl[, nu:=1/(1+exp( -predict( ml.zaga.dt, data=ml.zaga.dt.data, newdata=as.data.frame(subset(Y.brl, select=c(brlz, dt))), what='nu', type='link' ) ))]		
		Y.brl[, score.brl.TPp:= pZAGA(brlz, mu=mu, sigma=sigma, nu=nu, lower.tail=FALSE)/(1-nu)]
		Y.brl[, score.brl.TPd:= score.brl.TPp]
	}
	if(substr(method,1,2)=='3h')	#dont ignore zeros - use zero adjusted Gamma. Ignore after cART initiation. 
	{		
		resume<- 1
		if(resume && !is.na(save.file.3ha))
		{
			options(show.error.messages = FALSE)		
			readAttempt		<- try(suppressWarnings(load(save.file.3ha)))
			if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file.3ha))					
			options(show.error.messages = TRUE)		
		}
		if(!resume || is.na(save.file.3ha) || inherits(readAttempt, "try-error"))
		{
			brlz.dt	<- 0.0005
			brlz.dt	<- 0.01
			Y.brlc	<- data.table(dt= seq(0,20,0.0125), b4T='both')
			Y.brlc[, mu:=exp( predict( ml.zaga.dt, data=ml.zaga.dt.data, newdata=as.data.frame(subset(Y.brlc, select=dt)), what='mu', type='link' ) )]
			Y.brlc[, sigma:=exp( predict( ml.zaga.dt, data=ml.zaga.dt.data, newdata=as.data.frame(subset(Y.brlc, select=dt)), what='sigma', type='link' ) )]
			Y.brlc[, nu:=1/(1+exp( -predict( ml.zaga.dt, data=ml.zaga.dt.data, newdata=as.data.frame(subset(Y.brlc, select=dt)), what='nu', type='link' ) ))]			
			tmp		<- data.table(dt= seq(0,20,0.0125), b4T='one')
			tmp[, mu:=exp( predict( ml.zaga.dt.one, data=ml.zaga.dt.data.one, newdata=as.data.frame(subset(tmp, select=dt)), what='mu', type='link' ) )]
			tmp[, sigma:=exp( predict( ml.zaga.dt.one, data=ml.zaga.dt.data.one, newdata=as.data.frame(subset(tmp, select=dt)), what='sigma', type='link' ) )]
			tmp[, nu:=1/(1+exp( -predict( ml.zaga.dt.one, data=ml.zaga.dt.data.one, newdata=as.data.frame(subset(tmp, select=dt)), what='nu', type='link' ) ))]		
			Y.brlc	<- rbind(Y.brlc, tmp)			
			Y.brlc	<- merge(Y.brlc, as.data.table( expand.grid(dt=seq(0,20,0.0125), brlz=seq(0.0005,0.25,brlz.dt) )), by='dt', allow.cartesian=TRUE )
			Y.brlc	<- subset(Y.brlc, dt>15)	#subset(Y.brlc, dt<=3)
			Y.brlc[, dZAGA:=dZAGA(brlz, mu=mu, sigma=sigma, nu=nu)/(1-nu)]
			Y.brlc	<- Y.brlc[, {
						tmp			<- c(dt,brlz)
						tmp2		<- b4T						
						tmp			<- subset(Y.brlc, dt<=tmp[1] & brlz==tmp[2] & b4T==tmp2)
						if(nrow(tmp)>1)	
						{
							tmp		<- cbind( subset(tmp, select=c(dt, dZAGA)), data.table( dtms= tmp[,rev(dt)], dZAGA2= tmp[,rev(dZAGA)]) )
							dZAGA2	<- tmp[-nrow(tmp), sum(dZAGA*dZAGA2)*diff(dt[1:2])]
						}		
						else
							dZAGA2	<- tmp[1, dZAGA*dZAGA*0.0125]
						list( dZAGA2=dZAGA2, dZAGA=dZAGA, mu=mu, sigma=sigma, nu=nu )
					}, by=c('b4T','brlz', 'dt')]
			Y.brlc	<- merge( subset(Y.brlc, select=c(b4T, brlz, dt, dZAGA=dZAGA, mu=mu, sigma=sigma, nu=nu)), Y.brlc[, list(brlz=brlz, dZAGA2=dZAGA2/sum(dZAGA2)), by=c('b4T','dt')], by=c('b4T','brlz', 'dt') )			
			Y.brlc	<- merge( Y.brlc, Y.brlc[, list(brlz=brlz, sZAGA2=1-cumsum(dZAGA2)), by=c('b4T','dt')], by=c('b4T','brlz', 'dt') )
			Y.brlc[, sZAGA:= pZAGA(brlz, mu=mu, sigma=sigma, nu=nu, lower.tail=FALSE)]
			if(!is.na(save.file.3ha))
			{
				cat(paste('\nsave Y.brlc to file',save.file.3ha))
				save(Y.brlc, file=save.file.3ha)
			}	
			stop()
						
		}		
		#score by dt
		ggplot( data=Y.brlc, aes(x=brlz, y=yc, group=dt, colour=dt))	+
				scale_colour_brewer(name='years between\nsequence sampling dates', palette='Set1') +
				theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(13,'mm')) +
				labs(x="substitutions / site", y='probability of\ndirect HIV-1 transmission') +
				geom_line()

	}
	if(!is.na(plot.file.score) && substr(method,1,2)%in%c('3c','3e'))
	{
		plot.df			<- subset(Y.rawbrl.linked, b4T!='none' )
		plot.df[, score.brl.TPp:= pgamma(plot.df[,brl], shape=2*mle.ga$estimate['shape'], rate=mle.ga$estimate['rate'], lower.tail=FALSE)]
		tmp<- plot.df[, which(b4T=='one')]
		set(plot.df, tmp, 'score.brl.TPp', pgamma(plot.df[tmp,brl], shape=2*mle.ga.one$estimate['shape'], rate=mle.ga.one$estimate['rate'], lower.tail=FALSE) )
		setkey(plot.df, b4T.long, brl)
		tmp				<- plot.df[, list(brl=brl, empirical=seq_len(length(brl))/length(brl)), by='b4T.long']			
		tmp				<- tmp[, list(empirical=tail(empirical,1)), by=c('b4T.long','brl')]
		set(tmp, NULL, 'empirical', tmp[,1-empirical])
		plot.df			<- merge(plot.df, tmp, by=c('b4T.long','brl'))		
		plot.df[, shape:= mle.ga$estimate['shape']]
		plot.df[, rate:= mle.ga$estimate['rate']]		
		tmp				<- plot.df[,which(b4T!='both')]
		set(plot.df, tmp, 'shape', mle.ga.one$estimate['shape'])
		set(plot.df, tmp, 'rate', mle.ga.one$estimate['rate'])		
		tmp2			<- plot.df[, list( b4T=b4T, brl=brl, v=1-pgamma(brl, rate=rate, shape=shape), group='Gamma' ), by='b4T.long']		
		tmp				<- subset(plot.df, select=c(b4T.long, b4T, brl, score.brl.TPp))
		tmp[, group:='convoluted\nGamma']
		setnames(tmp, 'score.brl.TPp', 'v')
		plot.df			<- subset(plot.df, select=c(b4T.long, b4T, brl, empirical))
		plot.df[, group:='empirical']
		setnames(plot.df, 'empirical', 'v')		
		plot.df			<- do.call('rbind', list(plot.df, tmp, tmp2) )
		
		ggplot( plot.df , aes(x=brl, y=v, colour=group, group=group)) + scale_colour_brewer(palette='Set1') + labs(x="substitutions / site", y='survival') +
				theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(13,'mm')) +
				geom_line() + facet_grid(. ~ b4T.long, scales='free_x', margins=FALSE)
		ggsave(file=paste(substr(plot.file.score,1,nchar(plot.file.score)-4),'_cmpwithinhost','.pdf',sep=''), w=8, h=6)
		
		
		Y.brl[, dummy:= score.brl.TPp*30]	
		set(Y.brl, Y.brl[, which(b4T=='both')], 'b4T.long', 'both sampling dates\nbefore cART initiation')
		set(Y.brl, Y.brl[, which(b4T!='both')], 'b4T.long', 'sampling dates before\nand after cART initiation')
		ggplot(data=Y.brl, aes(x=brl, fill=b4T.long)) + 
				theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(13,'mm')) + scale_fill_brewer(name='treatment status', palette='Paired') +
				labs(x="substitutions / site", y='Clustering sequence pairs\nwith at least one sequence\nfrom a recently infected MSM') +
				geom_histogram(aes(y= ..density..), binwidth=0.003, show_guide=FALSE) + geom_line(aes(x=brl, y=dummy)) + facet_grid(. ~ b4T.long, scales='free_y', margins=FALSE)
		ggsave(file=plot.file.score, w=8, h=6)			
	}		
	if(!is.na(plot.file.score) && substr(method,1,2)=='3g')
	{
		#ZAGA parameters
		plot.df	<- data.table(dt=seq(0,10,0.1))
		plot.df[, mu:=exp( predict( ml.zaga.dt, data=ml.zaga.dt.data, newdata=as.data.frame(subset(plot.df, select=c(dt))), what='mu', type='link' ) )]
		plot.df[, sigma:=exp( predict( ml.zaga.dt, data=ml.zaga.dt.data, newdata=as.data.frame(subset(plot.df, select=c(dt))), what='sigma', type='link' ) )]
		plot.df[, nu:=1/(1+exp( -predict( ml.zaga.dt, data=ml.zaga.dt.data, newdata=as.data.frame(subset(plot.df, select=c(dt))), what='nu', type='link' ) ))]				
		ggplot( data=melt( plot.df, id.vars='dt', variable.name='parameter' ), aes(x=dt, y=value, group=parameter, colour=parameter)) +
				theme(legend.justification=c(1,0.4), legend.position=c(1,0.4), legend.key.size=unit(13,'mm')) + 
				scale_colour_brewer(name='parameters of\nzero-inflated Gamma', palette='Set1') +
				labs(x="years between\nsequence sampling dates", y='parameter value') +
				geom_line()
		ggsave(file=paste(substr(plot.file.score,1,nchar(plot.file.score)-4),'_ZAGAparam','.pdf',sep=''), w=4, h=6)
		#score by dt
		plot.df	<- data.table(dt=seq(1,5,1))
		plot.df[, mu:=exp( predict( ml.zaga.dt, data=ml.zaga.dt.data, newdata=as.data.frame(subset(plot.df, select=c(dt))), what='mu', type='link' ) )]
		plot.df[, sigma:=exp( predict( ml.zaga.dt, data=ml.zaga.dt.data, newdata=as.data.frame(subset(plot.df, select=c(dt))), what='sigma', type='link' ) )]
		plot.df[, nu:=1/(1+exp( -predict( ml.zaga.dt, data=ml.zaga.dt.data, newdata=as.data.frame(subset(plot.df, select=c(dt))), what='nu', type='link' ) ))]				
		plot.df	<- merge(plot.df, as.data.table( expand.grid(brl=seq(0,0.05,0.001), dt=seq(1,5,1)) ), by='dt')
		plot.df[, y:= pZAGA(brl, mu=mu, sigma=sigma, nu=nu, lower.tail=FALSE)]
		plot.df[, yc:= pZAGA(brl, mu=mu, sigma=sigma, nu=nu, lower.tail=FALSE)/(1-nu)]
		set(plot.df, NULL, 'dt', plot.df[, factor(dt)])
		ggplot( data=plot.df, aes(x=brl, y=yc, group=dt, colour=dt))	+
				scale_colour_brewer(name='years between\nsequence sampling dates', palette='Set1') +
				theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(13,'mm')) +
				labs(x="substitutions / site", y='probability of\ndirect HIV-1 transmission') +
				geom_line()
		ggsave(file=paste(substr(plot.file.score,1,nchar(plot.file.score)-4),'_ZAGAdt','.pdf',sep=''), w=4, h=6)			
	}	
	if(!is.na(plot.file.score) && substr(method,1,2)%in%c('3i'))
	{
		Y.brl[, b4T.long:=b4T]
		ggplot(Y.brl, aes(x=brl, y=score.brl.TPp, group=b4T.long, colour=b4T.long)) + geom_line() +
				coord_cartesian(xlim=c(0, 0.16)) + scale_x_continuous(breaks=seq(0,0.2,0.02), minor_breaks=seq(0, 0.2, 0.005)) + scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0, 1, 0.05)) +
				theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(11,'mm')) +
				scale_colour_brewer(palette='Set1', name='treatment status at\nsequence sampling time') +
				labs(x="between-host divergence", y=expression('Conditional probability score ('*y[ijt]^{C}*')'))
		ggsave(file=plot.file.score, w=10, h=6)		
	}	
	if(!is.na(plot.file.score) && substr(method,1,2)%in%c('3j'))
	{
		Y.brl[, b4T.long:=b4T]
		tmp	<- c('both sampling times\nbefore ART start','one sampling time\nafter ART start','both sampling times\nafter ART start')
		set(Y.brl, Y.brl[, which(b4T.long=='both')], 'b4T.long', tmp[1])		
		set(Y.brl, Y.brl[, which(b4T.long=='mix')], 'b4T.long', tmp[2])
		set(Y.brl, Y.brl[, which(b4T.long=='one')], 'b4T.long', tmp[3])
		set(Y.brl, NULL, 'b4T.long', Y.brl[, factor(b4T.long, levels=tmp, labels=tmp)])
		ggplot(Y.brl, aes(x=brl, y=score.brl.TPp, group=b4T.long, colour=b4T.long)) + geom_line() +
				coord_cartesian(xlim=c(0, 0.16)) + scale_x_continuous(breaks=seq(0,0.2,0.02), minor_breaks=seq(0, 0.2, 0.005)) + scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0, 1, 0.05)) +
				theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(11,'mm')) +
				scale_colour_brewer(palette='Set1', name='treatment status at\nsequence sampling time') +
				labs(x="between-host divergence", y=expression('Conditional probability score ('*y[ijt]^{C}*')'))
		ggsave(file=plot.file.score, w=10, h=6)		
	}
	if(!is.na(plot.file.score) && substr(method,1,2)%in%c('3k'))
	{
	
		plot.df	<- data.table(brlz=rep( seq(0,0.15,0.005), each=3), ART.C=c('Yes','Yes','No'), t.ART.C=c('Yes','No','No'), score.brl.TPp=NA_real_)		
		plot.df	<- merge(plot.df, data.table(brlz=rep( seq(0,0.15,0.005), each=3), brl.bwhost.multiplier=c(1.5, 2, 2.5)), by='brlz', allow.cartesian=TRUE)
		#	set score for cases F/F
		tmp				<- plot.df[, which(t.ART.C=='No' & ART.C=='No')]
		tmp				<- 2*ml.zaga.CN.pa['nu']*(1-ml.zaga.CN.pa['nu'])*pGA( plot.df[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.CN.pa['mu'], sigma=ml.zaga.CN.pa['sigma'] ) + (1-ml.zaga.CN.pa['nu'])*(1-ml.zaga.CN.pa['nu'])*pgamma( plot.df[tmp,brlz*brl.bwhost.multiplier],  shape=2*ml.zaga.CN.pa['shape'], scale=ml.zaga.CN.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.CN.pa['nu']*ml.zaga.CN.pa['nu'])
		set(plot.df, plot.df[, which(t.ART.C=='No' & ART.C=='No')], 'score.brl.TPp', 1-tmp )
		#	set score for cases T/T
		tmp				<- plot.df[, which(t.ART.C=='Yes' & ART.C=='Yes')]
		tmp				<- 2*ml.zaga.CY.pa['nu']*(1-ml.zaga.CY.pa['nu'])*pGA( plot.df[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.CY.pa['mu'], sigma=ml.zaga.CY.pa['sigma'] ) + (1-ml.zaga.CY.pa['nu'])*(1-ml.zaga.CY.pa['nu'])*pgamma( plot.df[tmp,brlz*brl.bwhost.multiplier],  shape=2*ml.zaga.CY.pa['shape'], scale=ml.zaga.CY.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.CY.pa['nu']*ml.zaga.CY.pa['nu'])				
		set(plot.df, plot.df[, which(t.ART.C=='Yes' & ART.C=='Yes')], 'score.brl.TPp', 1-tmp )			
		tmp				<- plot.df[, which(t.ART.C!=ART.C)]
		#	set score for cases T/F
		tmp				<- ml.zaga.CN.pa['nu']*(1-ml.zaga.CN.pa['nu'])*pGA( plot.df[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.CN.pa['mu'], sigma=ml.zaga.CN.pa['sigma'] ) +
				ml.zaga.CY.pa['nu']*(1-ml.zaga.CY.pa['nu'])*pGA( plot.df[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.CY.pa['mu'], sigma=ml.zaga.CY.pa['sigma'] ) +
				(1-ml.zaga.CN.pa['nu'])*(1-ml.zaga.CY.pa['nu'])*ml.zaga.pmix( plot.df[tmp,brlz*brl.bwhost.multiplier] )	
		tmp				<- tmp / ( ml.zaga.CN.pa['nu']*(1-ml.zaga.CN.pa['nu']) + ml.zaga.CY.pa['nu']*(1-ml.zaga.CY.pa['nu']) + (1-ml.zaga.CN.pa['nu'])*(1-ml.zaga.CY.pa['nu']) )
		set(plot.df, plot.df[, which(t.ART.C!=ART.C)], 'score.brl.TPp', 1-tmp )
		#		
		plot.df[, ART.C.long:=NA_character_]		
		tmp	<- c('both sampling times\nbefore treatment interruption or failure','one sampling time\nafter treatment interruption or failure','both sampling times\nafter treatment interruption or failure')
		set(plot.df, plot.df[, which(ART.C=='No' & t.ART.C=='No')], 'ART.C.long', tmp[1])				
		set(plot.df, plot.df[, which(ART.C=='Yes' & t.ART.C=='No')], 'ART.C.long', tmp[2])
		set(plot.df, plot.df[, which(ART.C=='Yes' & t.ART.C=='Yes')], 'ART.C.long', tmp[3])
		set(plot.df, NULL, 'ART.C.long', plot.df[, factor(ART.C.long, levels=tmp, labels=tmp)])
		set(plot.df, NULL, 'brl.bwhost.multiplier', plot.df[, factor(brl.bwhost.multiplier)])
		ggplot(plot.df, aes(x=brlz, y=score.brl.TPp, group=interaction(ART.C.long, brl.bwhost.multiplier), colour=ART.C.long, linetype=brl.bwhost.multiplier)) + geom_line() + 
				coord_cartesian(xlim=c(0, 0.1)) + scale_x_continuous(breaks=seq(0,0.2,0.02), minor_breaks=seq(0, 0.2, 0.005) ) + scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0, 1, 0.05)) +
				#theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(11,'mm')) +
				theme(legend.position='bottom', legend.key.size=unit(11,'mm'), strip.background = element_blank(), strip.text = element_blank(), panel.margin = unit(2, "lines")) +
				scale_colour_brewer(palette='Set1', name='treatment status at\nsequence sampling time') +
				scale_linetype_manual(values=c(2,1,4)) + 
				facet_grid(.~brl.bwhost.multiplier, margins=FALSE ) +
				labs(x=expression(atop("evolutionary distance "*d[ij],"(estimated substitutions/site)")), y=expression(y[ij]), linetype=expression(alpha))				
		ggsave(file=plot.file.score, w=11, h=8)	
		#
		#	PDF
		#
		plot.df<- data.table(brlz=rep( seq(0,0.15,0.001), each=3), ART.C=c('Yes','Yes','No'), t.ART.C=c('Yes','No','No'), score.brl.TPp=NA_real_)
		#	set score for cases F/F
		tmp				<- plot.df[, which(t.ART.C=='No' & ART.C=='No')]
		tmp				<- 2*ml.zaga.CN.pa['nu']*(1-ml.zaga.CN.pa['nu'])*dGA( plot.df[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.CN.pa['mu'], sigma=ml.zaga.CN.pa['sigma'] ) + (1-ml.zaga.CN.pa['nu'])*(1-ml.zaga.CN.pa['nu'])*dgamma( plot.df[tmp,brlz*brl.bwhost.multiplier],  shape=2*ml.zaga.CN.pa['shape'], scale=ml.zaga.CN.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.CN.pa['nu']*ml.zaga.CN.pa['nu'])
		set(plot.df, plot.df[, which(t.ART.C=='No' & ART.C=='No')], 'score.brl.TPp', tmp )
		#	set score for cases T/T
		tmp				<- plot.df[, which(t.ART.C=='Yes' & ART.C=='Yes')]
		tmp				<- 2*ml.zaga.CY.pa['nu']*(1-ml.zaga.CY.pa['nu'])*dGA( plot.df[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.CY.pa['mu'], sigma=ml.zaga.CY.pa['sigma'] ) + (1-ml.zaga.CY.pa['nu'])*(1-ml.zaga.CY.pa['nu'])*dgamma( plot.df[tmp,brlz*brl.bwhost.multiplier],  shape=2*ml.zaga.CY.pa['shape'], scale=ml.zaga.CY.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.CY.pa['nu']*ml.zaga.CY.pa['nu'])				
		set(plot.df, plot.df[, which(t.ART.C=='Yes' & ART.C=='Yes')], 'score.brl.TPp', tmp )			
		tmp				<- plot.df[, which(t.ART.C!=ART.C)]
		#	set score for cases T/F
		tmp				<- ml.zaga.CN.pa['nu']*(1-ml.zaga.CN.pa['nu'])*dGA( plot.df[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.CN.pa['mu'], sigma=ml.zaga.CN.pa['sigma'] ) +
								ml.zaga.CY.pa['nu']*(1-ml.zaga.CY.pa['nu'])*dGA( plot.df[tmp,brlz*brl.bwhost.multiplier],  mu=ml.zaga.CY.pa['mu'], sigma=ml.zaga.CY.pa['sigma'] ) +
								(1-ml.zaga.CN.pa['nu'])*(1-ml.zaga.CY.pa['nu'])*ml.zaga.dmix( plot.df[tmp,brlz*brl.bwhost.multiplier] )	
		tmp				<- tmp / ( ml.zaga.CN.pa['nu']*(1-ml.zaga.CN.pa['nu']) + ml.zaga.CY.pa['nu']*(1-ml.zaga.CY.pa['nu']) + (1-ml.zaga.CN.pa['nu'])*(1-ml.zaga.CY.pa['nu']) )
		set(plot.df, plot.df[, which(t.ART.C!=ART.C)], 'score.brl.TPp', tmp )
		#
		plot.df[, ART.C.long:=NA_character_]		
		tmp	<- c('both sampling times\nbefore treatment interruption or failure','one sampling time\nafter treatment interruption or failure','both sampling times\nafter treatment interruption or failure')
		set(plot.df, plot.df[, which(ART.C=='No' & t.ART.C=='No')], 'ART.C.long', tmp[1])				
		set(plot.df, plot.df[, which(ART.C=='Yes' & t.ART.C=='No')], 'ART.C.long', tmp[2])
		set(plot.df, plot.df[, which(ART.C=='Yes' & t.ART.C=='Yes')], 'ART.C.long', tmp[3])
		set(plot.df, NULL, 'ART.C.long', plot.df[, factor(ART.C.long, levels=tmp, labels=tmp)])
		ggplot(plot.df, aes(x=brlz, y=score.brl.TPp, group=ART.C.long, colour=ART.C.long)) + geom_line() + 
				#geom_point() +
				coord_cartesian(xlim=c(0, 0.16)) + scale_x_continuous(breaks=seq(0,0.2,0.02), minor_breaks=seq(0, 0.2, 0.005)) + #scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0, 1, 0.05)) +
				theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(11,'mm')) +
				scale_colour_brewer(palette='Set1', name='treatment status at\nsequence sampling time') +
				labs(x="between-host divergence", y='generation divergence distribution\n (p.d.f.)')
		ggsave(file=paste(substr(plot.file.score,1,nchar(plot.file.score)-4),'_p','.pdf',sep=''), w=4, h=6)
	}
	
	if(!is.na(plot.file.score) && substr(method,1,2)%in%c('3d','3f'))
	{
		plot.df			<- subset(Y.rawbrl.linked, b4T!='none' )
		tmp				<- 2*ml.zaga.pa['nu']*(1-ml.zaga.pa['nu'])*pGA( plot.df[,brlz],  mu=ml.zaga.pa['mu'], sigma=ml.zaga.pa['sigma'] ) + (1-ml.zaga.pa['nu'])*(1-ml.zaga.pa['nu'])*pgamma( plot.df[,brlz],  shape=2*ml.zaga.pa['shape'], scale=ml.zaga.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.pa['nu']*ml.zaga.pa['nu'])				
		plot.df[, score.brl.TPp:= 1-tmp]
		tmp				<- plot.df[, which(b4T!='both')]
		tmp				<- 2*ml.zaga.one.pa['nu']*(1-ml.zaga.one.pa['nu'])*pGA( plot.df[tmp,brlz],  mu=ml.zaga.one.pa['mu'], sigma=ml.zaga.one.pa['sigma'] ) + (1-ml.zaga.one.pa['nu'])*(1-ml.zaga.one.pa['nu'])*pgamma( plot.df[tmp,brlz],  shape=2*ml.zaga.one.pa['shape'], scale=ml.zaga.one.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.one.pa['nu']*ml.zaga.one.pa['nu'])				
		set(plot.df, plot.df[, which(b4T!='both')], 'score.brl.TPp', 1-tmp )
		setkey(plot.df, b4T.long, brlz)
		tmp				<- plot.df[, list(brl=brl, empirical=seq_len(length(brlz))/length(brlz)), by='b4T.long']			
		tmp				<- tmp[, list(empirical=tail(empirical,1)), by=c('b4T.long','brl')]
		set(tmp, NULL, 'empirical', tmp[,1-empirical])
		plot.df			<- merge(plot.df, tmp, by=c('b4T.long','brl'))
		plot.df[, mu:= ml.zaga.pa['mu']]
		plot.df[, nu:= ml.zaga.pa['nu']]
		plot.df[, sigma:= ml.zaga.pa['sigma']]
		tmp				<- plot.df[,which(b4T!='both')]
		set(plot.df, tmp, 'mu', ml.zaga.one.pa['mu'])
		set(plot.df, tmp, 'nu', ml.zaga.one.pa['nu'])
		set(plot.df, tmp, 'sigma', ml.zaga.one.pa['sigma'])
		
		tmp2			<- plot.df[, list( b4T=b4T, brl=brl, v=1-pZAGA(brl, mu=mu, sigma=sigma, nu=nu), group='Zero-inflated Gamma' ), by='b4T.long']		
		tmp				<- subset(plot.df, select=c(b4T.long, b4T, brl, score.brl.TPp))
		tmp[, group:='convoluted\nzero-inflated Gamma']
		setnames(tmp, 'score.brl.TPp', 'v')
		plot.df			<- subset(plot.df, select=c(b4T.long, b4T, brl, empirical))
		plot.df[, group:='empirical']
		setnames(plot.df, 'empirical', 'v')		
		plot.df			<- do.call('rbind', list(plot.df, tmp, tmp2) )
				
		ggplot( plot.df , aes(x=brl, y=v, colour=group, group=group)) + scale_colour_brewer(palette='Set1') + labs(x="substitutions / site", y='survival') +
				theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(13,'mm')) +
				geom_line() + facet_grid(. ~ b4T.long, scales='free_x', margins=FALSE)
		ggsave(file=paste(substr(plot.file.score,1,nchar(plot.file.score)-4),'_cmpwithinhost','.pdf',sep=''), w=8, h=6)
		
		
		Y.brl[, dummy:= score.brl.TPp*30]	
		set(Y.brl, Y.brl[, which(b4T=='both')], 'b4T.long', 'both sampling dates\nbefore cART initiation')
		set(Y.brl, Y.brl[, which(b4T!='both')], 'b4T.long', 'sampling dates before\nand after cART initiation')
		#tmp				<- Y.brl[, which(b4T=='both')]
		#set(Y.brl, tmp, 'dummy', Y.brl[tmp, score.brl.TPp*length(which(brl<=0.002 & b4T=='both'))])
		#tmp				<- Y.brl[, which(b4T!='both')]
		#set(Y.brl, tmp, 'dummy', Y.brl[tmp, score.brl.TPp*length(which(brl<=0.002 & b4T!='both'))])		
		ggplot(data=Y.brl, aes(x=brl, fill=b4T.long)) + 
				theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(13,'mm')) + scale_fill_brewer(name='treatment status', palette='Paired') +
				labs(x="substitutions / site", y='Clustering sequence pairs\nwith at least one sequence\nfrom a recently infected MSM') +
				geom_histogram(aes(y= ..density..), binwidth=0.003, show_guide=FALSE) + geom_line(aes(x=brl, y=dummy)) + facet_grid(. ~ b4T.long, scales='free_y', margins=FALSE)
		ggsave(file=plot.file.score, w=8, h=6)
		
		#pdf(plot.file.score, w=5, h=5)
		#cat(paste('\nplot to file',plot.file.score))
		#require(RColorBrewer)
		#par(mar=c(4,4,0.5,0.5))
		#cols		<- brewer.pal(4, 'Set1')[c(1,3,2)]
		#legend.txt	<- c('same host', 'potential transmission pairs', 'branch length weight')	
		#xlim		<- c(0, max(Y.brl[, max(brl)],  Y.rawbrl.linked[, max(brl)])*1.1 )
		#tmp			<- seq(from=xlim[1], to=xlim[2], by=diff(xlim)/200)
		#hist( Y.rawbrl.linked[, brl], breaks=tmp , col=my.fade.col(cols[1],0.5), add=0, freq=0, xlab='branch length', main='', ylim=c(0,50) )
		#tmp2		<- hist( Y.brl[, brl], breaks=tmp, col=my.fade.col(cols[2],0.5), freq=0, add=1 )
		#setkey(Y.brl, brl)
		#lines(subset(Y.brl, b4T=='both')[, brl], subset(Y.brl, b4T=='both')[, score.brl.TPp]*max(tmp2$density), col=cols[3], lwd=2)
		#lines(subset(Y.brl, b4T=='one')[, brl], subset(Y.brl, b4T=='one')[, score.brl.TPp]*max(tmp2$density), col=cols[3], lwd=2)
		#legend('topright', bty='n', fill= cols, border=NA, legend=legend.txt)		
		#dev.off()		
	}		
	Y.brl		<- subset( Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, score.brl.TPd, score.brl.TPp, score.brl.TN, score.brl.TP, brl) )
	cat(paste('\nReturn brl score for #t.Patient seqs=',Y.brl[, length(unique(t.FASTASampleCode))]))
	Y.brl		
}
######################################################################################
project.athena.Fisheretal.Y.rm.missedtransmitter<- function(YX.tpairs, df.all, Y.brl, Y.U, Y.coal=NULL, thresh.pcoal=0.5, cut.brl=1e-3, any.pos.grace.yr= 3.5, rm.zero.score=FALSE, plot.file=NA, pyiw=NULL, t.period=0.25)
{
	#thresh.pcoal=0.5; cut.brl=1e-3; any.pos.grace.yr= 3.5; rm.zero.score=FALSE
	missed					<- merge(YX.tpairs, Y.brl, by=c('FASTASampleCode','t.FASTASampleCode'), all.x=TRUE)
	if(!is.null(Y.coal))
		missed				<- merge(missed, Y.coal, by=c('FASTASampleCode','t.FASTASampleCode'), all.x=TRUE)
	missed					<- merge(missed, Y.U, by=c('Patient','t.Patient'), allow.cartesian=TRUE, all.x=TRUE)
	if(!is.null(Y.coal))
		missed				<- subset(missed, select=c(cluster, Patient, t.Patient, FASTASampleCode, t.FASTASampleCode, t, score.brl.TN, score.brl.TP, score.brl.TPp, coal.after.i.AnyPos_T1, coal.after.t.NegT, score.Inf, score.t.inf, brl, class))
	if(is.null(Y.coal))
		missed				<- subset(missed, select=c(cluster, Patient, t.Patient, FASTASampleCode, t.FASTASampleCode, t, score.brl.TN, score.brl.TP, score.brl.TPp, score.Inf, score.t.inf, brl, class))
	#
	set(missed, missed[, which(is.na(brl))], 'score.brl.TN', missed[, max(score.brl.TN,na.rm=TRUE)])		#deselected short or recombinant seqs
	set(missed, missed[, which(is.na(brl))], 'brl', missed[, max(brl,na.rm=TRUE)])							#deselected short or recombinant seqs
	#if(any(is.na(missed)))	stop('unexpected NA entries in missed')
	missed[, score.Y:=1.]
	missed[, score.t.inf:=NULL]
	#	exclude coalescence within infected
	if(!is.null(Y.coal))
	{
		set(missed, NULL, 'coal.after.t.NegT', missed[,coal.after.t.NegT-coal.after.i.AnyPos_T1])
		set(missed, missed[, which(coal.after.t.NegT<0)], 'coal.after.t.NegT', NA_real_)		#transmission could simply have other direction		
	}			
	#	exclude transmitters if date of diagnosis exceeds date of diagnosis of RI by more than 'any.pos.grace.yr'
	if(!is.na(any.pos.grace.yr) & !is.infinite(any.pos.grace.yr))
	{
		tmp		<- subset( df.all, select=c(Patient, AnyPos_T1) )
		setkey(tmp, Patient)
		missed	<- merge(missed, unique(tmp), by= 'Patient')
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		missed	<- merge(missed, unique(tmp), by= 't.Patient')
		cat(paste('\nnumber of pairs not meeting AnyPos_T1 cut=', missed[,length(which(t.AnyPos_T1>=any.pos.grace.yr+AnyPos_T1))] ))		
		set(missed, missed[,which(t.AnyPos_T1>=any.pos.grace.yr+AnyPos_T1)],'score.Y',0.)
	}
	#	simple rule for screening out missed sources, using genetics and dates
	if(1)
	{	
		set(missed, missed[, which(score.brl.TN>=cut.brl)], 'score.Y', NA_real_)
		if(!is.null(pyiw))
			pyiw	<- rbind(pyiw, data.table(pop='with.unlinkedbydeath', t.py= nrow(subset(missed, score.Y>0))*t.period, t.n= subset(missed, score.Y>0)[, length(unique(t.Patient))], n= subset(missed, score.Y>0)[, length(unique(Patient))] ))		
		if(!is.null(Y.coal))
		{	
			cat(paste('\nprop of PYIW that meet the BRL.TN & COAL cutoff:', nrow(subset( missed, score.brl.TN<cut.brl &  coal.after.t.NegT>thresh.pcoal))/nrow(missed)  ))
			set(missed, missed[, which(!is.na(coal.after.t.NegT) & coal.after.t.NegT<=thresh.pcoal)], 'score.Y', NA_real_)
			if(!is.null(pyiw))
				pyiw	<- rbind(pyiw, data.table(pop='with.coal', t.py= nrow(subset(missed, score.Y>0))*t.period, t.n= subset(missed, score.Y>0)[, length(unique(t.Patient))], n= subset(missed, score.Y>0)[, length(unique(Patient))] ))
			cat(paste('\nnumber of pn PYIW with score.Y>0:',missed[, length(which(is.na(coal.after.t.NegT) & class=='pn' & score.Y>0))]))			
		}												
	}
	#	simple rule for setting zero: brl score Y essentially zero. Define essentially zero by 1/non NA rows.
	#	use score instead of within host BRL because these are not 'convoluted'
	#	not using regression any longer, switch this off
	if(0)
	{		
		tmp	<- missed[, 1/length(which(!is.na(score.Y)))]
		
		cat(paste('\nnumber of rows with non NA score - these are not tractable in the regression part, n=',1/tmp))
		tmp	<- missed[, which(!is.na(score.Y) & score.brl.TPp<tmp)]
		cat(paste('\nnumber of rows with !NA score and Yijt < 1/nrow(!is.na(Y)). Setting to NA score, n=',length(tmp)))
		set(missed, tmp, 'score.Y', NA_real_)
		tmp	<- missed[, 1/length(which(!is.na(score.Y)))]
		cat(paste('\nnumber of rows with non NA score - these are not tractable in the regression part, n=',1/tmp))
		tmp	<- missed[, which(is.na(score.Y) & score.brl.TPp<tmp)]
		cat(paste('\nnumber of rows with NA score, n=', missed[, length(which(is.na(score.Y)))]  ))
		cat(paste('\nnumber of rows with NA score and Yijt < 1/nrow(!is.na(Y)). Setting to ZERO, n=',length(tmp)))
		set(missed, tmp, 'score.Y', 0.)
		#	add also small COAL
		tmp		<- missed[,  mean(!is.na(coal.after.t.NegT) & coal.after.t.NegT>thresh.pcoal)]		#proportion of rows that pass coal criterion for pos Yijt
		cat(paste('\nproportion of rows that satisfy COAL criterion (for Yijt>0), p=',tmp))
		tmp		<- missed[, quantile(coal.after.t.NegT, p=tmp)]			#cut-off for zero Yijt so that the same proportion is retained
		cat(paste('\ncorresponding cut off for Yijt=0 is ',tmp))		#roughly p(COAL within transmitter)<0.25	
		tmp		<- missed[,  which(is.na(score.Y) & !is.na(coal.after.t.NegT) & coal.after.t.NegT<tmp)]		
		cat(paste('\nnumber of rows with NA score and COAL satisfying zero cut-off. Setting to ZERO, n=',length(tmp)))
		set(missed, tmp, 'score.Y', 0.)		
	}
	cat(paste('\nnumber of rows with NA score (will be discarded, undeterminate), n=',nrow(subset(missed, is.na(score.Y)))  ))
	cat(paste('\nproportion of rows with NA score (will be discarded, undeterminate), p=',nrow(subset(missed, is.na(score.Y)))/nrow(missed)  ))
	#
	setkey(missed, FASTASampleCode, t.FASTASampleCode)
	missed		<- unique(missed)
			

	if(!is.na(plot.file))
	{
		missed[, c2:= coal.after.t.NegT]
		set(missed, missed[, which(c2>0.999 | c2<0.001)],'c2', NA_real_)
		require(betareg)		
		tmp			<- betareg('c2 ~ brl', data=as.data.frame(subset(missed, select=c(c2, brl))))
		summary(tmp)
		#Pseudo R-squared: 0.04968
		tmp$predict	<- predict(tmp, type='response')
		set(missed, as.integer(names(tmp$predict)), 'c.b', tmp$predict)
		tmp$predict	<- predict(tmp, type='quantile', at=0.025)
		set(missed, as.integer(names(tmp$predict)), 'c.sl', tmp$predict)
		tmp$predict	<- predict(tmp, type='quantile', at=0.975)
		set(missed, as.integer(names(tmp$predict)), 'c.su', tmp$predict)
		
		ggplot(missed, aes(y=coal.after.t.NegT, x=brl)) + geom_point(size=0.7) + 
				geom_line(aes(x=brl, y=c.b), data=subset(missed, !is.na(c2)), col='blue') +
				geom_ribbon(aes(x=brl, ymin= c.sl, ymax= c.su, linetype=NA), alpha=0.2, data=subset(missed, !is.na(c2))) +
				labs(y="                Posterior probability that\n                coalesence is compatible with\n                direct HIV transmission", x="       Patristic distance\n       between transmitter and recipient sequence") + 				 
				theme(axis.title.x=element_text(hjust=0), axis.title.y=element_text(hjust=0))
		ggsave(paste(substr(plot.file,1,nchar(plot.file)-4),'correlation.pdf',sep=''), w=5, h=5)
	}	
	if(!is.na(plot.file))
	{
		missed[, cutoff:= 'posterior probability\n > 0.5']
		set(missed, missed[, which(coal.after.t.NegT<=thresh.pcoal)], 'cutoff', 'posterior probability\n <= 0.5')
		#set(missed, missed[, which(score.brl.TN>=cut.brl)], 'cutoff', 'GD > GD(Unlinked-by-Death)')
		tmp		<- missed[, list( FASTASampleCode=FASTASampleCode[which.min(brl)], t.FASTASampleCode=t.FASTASampleCode[which.min(brl)] ), by=c('Patient','t.Patient')]
		tmp		<- merge(subset(tmp, select=c(FASTASampleCode, t.FASTASampleCode)), missed, by=c('FASTASampleCode','t.FASTASampleCode'))
		setkey(tmp, Patient, t.Patient)
		set(tmp, NULL, 'cutoff', tmp[, factor(cutoff, levels=c('posterior probability\n > 0.5', 'posterior probability\n <= 0.5'))])
		ggplot(tmp, aes(x= brl, fill=cutoff)) + labs(x="Patristic distance", y='Number putative transmission pairs') +
				scale_fill_brewer(name='coalesence compatible with\n direct HIV transmission', palette='Paired') + theme(legend.justification=c(1,1), legend.position=c(1,1)) +
				geom_histogram(binwidth=0.002)	+ theme(legend.key.size=unit(8,'mm'))			
		cat(paste('\nsave to file',plot.file))
		ggsave(file=plot.file, w=5,h=5)			
	}
	missed		<- subset(missed, !is.na(score.Y))				
	#	not any longer needed because zero score.brl.TPp is now always NA or ZERO
	#if(1)
	#{
	#		cat(paste('\nnumber of pairs with zero brl.score=', missed[,length(which(score.brl.TPp==0 & score.Y>0))] ))
	#	missed	<- subset(missed, score.brl.TPp>0)
	#}
	#	effect on distribution of branch lengths
	
	#
	missed		<- subset( missed, select=c(cluster, Patient, FASTASampleCode, t.Patient, t.FASTASampleCode, score.Y, class))
	cat(paste('\n#intervals with score.Y>0, n=',  nrow(subset(missed, score.Y>0))))
	cat(paste('\n#infecteds with score.Y>0, n=',  subset(missed, score.Y>0)[, length(unique(Patient))]))
	cat(paste('\n#potential transmitters with score.Y>0, n=',  subset(missed, score.Y>0)[, length(unique(t.Patient))]))
	cat(paste('\n#intervals with score.Y==0, n=',  nrow(subset(missed, score.Y==0))))
	cat(paste('\n#infecteds with score.Y==0, n=',  subset(missed, score.Y==0)[, length(unique(Patient))]))
	cat(paste('\n#potential transmitters with score.Y==0, n=',  subset(missed, score.Y==0)[, length(unique(t.Patient))]))	
	if(rm.zero.score)
		missed	<- subset(missed, score.Y>0)
	if(is.null(pyiw))
		return(missed)
	if(!is.null(pyiw))
		return(list(missed=missed, pyiw=pyiw))
}
######################################################################################
project.athena.Fisheretal.YX.weight<- function(YX)
{	
	#	w.tn	number of transmission intervals between i and j
	#	w.in	number of transmitters for j
	#	w.i		number of transmission intervals between i and j, weighted by likelihood
	#	w.t		number of (i,j) (j,i) pairs
	#
	#	add tpair weight: every transmitter can infect only in one time period
	YX		<- merge(YX, YX[, list(w.tn= 1/length(t)),by=c('t.FASTASampleCode','FASTASampleCode')], by=c('t.FASTASampleCode','FASTASampleCode'))
	#	check weight
	if( abs(nrow(unique( subset(YX, select=c(Patient, t.Patient)) ))-YX[, sum(w.tn)])>5*EPS )	stop('unexpected weight')	
	print( YX[, table(w.tn)] )
	stopifnot( nrow(subset(YX, w.tn<=0))==0 )
	#	add infected weight: every infected can only be infected by one transmitter
	tmp		<- subset(YX, score.Y>0., select=c(Patient, t.Patient, score.brl.TPd))
	setkey(tmp, Patient, t.Patient)
	tmp		<- unique(tmp)	
	tmp		<- tmp[,	list(w.i=score.brl.TPd/sum(score.brl.TPd), w.in=1/length(t.Patient), t.Patient=t.Patient), by='Patient']
	if( tmp[,sum(w.i)]!=tmp[, length(unique(Patient))] )	stop('unexpected weight')	
	YX		<- merge(YX, tmp, by=c('Patient','t.Patient'), all.x=TRUE)
	set(YX, YX[, which(is.na(w.i))], c('w.i','w.in'), 0.)
	set(YX, NULL, 'w', YX[, w.tn*w.i])	
	if( abs(YX[,sum(w)]-subset(YX, score.Y>0)[, length(unique(Patient))])>5*EPS )	stop('unexpected weight')
	#	weights are now normalised per 'Patient' but there are fewer infection EVENTS because Patient may also be a transmitter to t.Patient
	triplet.weight				<- subset(YX, score.Y>0, select=c(Patient, t.Patient))	
	setkey(triplet.weight, Patient, t.Patient)
	tmp							<- copy(triplet.weight)
	setnames(tmp, colnames(tmp), paste('q.',colnames(tmp),sep='') )
	triplet.weight				<- triplet.weight[,	{ 	
				z	<- tmp[, which(q.Patient==t.Patient & q.t.Patient==Patient)]													
				list(equal.triplet.idx= ifelse(!length(z), NA_integer_, z))						
			}	, by=c('Patient','t.Patient')]														
	triplet.weight[, w.t:=equal.triplet.idx]
	triplet.weight				<- merge(triplet.weight, triplet.weight[, list(n.t.Patient=length(t.Patient)), by='Patient'], by='Patient')
	#	correct only mututally exclusive A->B and B->A
	# barplot( table( subset(triplet.weight, !is.na(equal.triplet.idx))[, n.t.Patient] ) )
	set(triplet.weight, triplet.weight[, which(!is.na(w.t) & n.t.Patient<2)], 'w.t', 2L)
	set(triplet.weight, triplet.weight[, which(!is.na(w.t) & n.t.Patient>=2)], 'w.t', 1L)
	set(triplet.weight, triplet.weight[, which(is.na(w.t))], 'w.t', 1L)	
	cat(paste('\nnumber of (i,j) and (j,i) pairs of the same infection event', triplet.weight[,length(which(w.t==2))] ))
	set(triplet.weight, NULL, 'w.t', triplet.weight[, 1/as.double(w.t)])
	YX		<- merge(YX, subset(triplet.weight, select=c(Patient, t.Patient, w.t)), by=c('Patient','t.Patient'), all.x=TRUE)
	set(YX, YX[, which(is.na(w.t))], 'w.t', 1.)
	set(YX, NULL, 'w', YX[, w*w.t])
	YX[,score.brl.TPd:=NULL]
	#	
	YX	
}
######################################################################################
project.athena.Fisheretal.Y.transmitterinfected<- function(YX.part1)
{
	Y.U			<- subset(YX.part1, select=c(Patient, t.Patient, t, U.score, score.Inf) )
	setkey(Y.U, Patient, t.Patient, t)
	Y.U			<- unique(Y.U)
	set(Y.U, Y.U[, which(is.na(U.score))], 'U.score', 1)
	setnames(Y.U, 'U.score', 'score.t.inf')
	Y.U
}
######################################################################################
project.athena.Fisheretal.Y.rawbrl<- function(YX.tpairs, indir, insignat, indircov, infilecov, infiletree, df.tpairs.tptn=NULL, save.file=NA, plot.file=NA, method.restrictTPtoRI=FALSE, resume=0)
{
	if(resume & !is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(!resume || is.na(save.file) || inherits(readAttempt, "try-error"))
	{
		require(adephylo)
		if(is.null(df.tpairs.tptn))
			df.tpairs.tptn	<- YX.tpairs
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
		tmp					<- unique( subset( df.tpairs.tptn, select='FASTASampleCode') ) 
		setnames(tmp, 'FASTASampleCode', 'sc.FASTASampleCode')	
		setnames(msm.unlinked.bytime, 'query.FASTASampleCode','sc.FASTASampleCode')
		msm.unlinked.bytime	<- merge(msm.unlinked.bytime, tmp, by='sc.FASTASampleCode')
		setnames(msm.unlinked.bytime, c('Patient','FASTASampleCode','DateDied'), c('t.Patient','t.FASTASampleCode','t.DateDied'))	
		cat(paste('\nfound viral sequences that are not linked to sequences of denominator group based on Died<NegT, n=',msm.unlinked.bytime[, length(unique(t.FASTASampleCode))]))
		cat(paste('\nfound unlinked pairs based on Died<NegT, n=',nrow(msm.unlinked.bytime)))
		#
		#	select any other within host sequence for every individual in the denominator population as 'true positive'
		#
		if(method.restrictTPtoRI)
		{
			tmp					<- unique( subset( df.tpairs.tptn, select=Patient ) )
			msm.linked			<- merge( msm.linked, tmp, by='Patient' )			
		}
		setkey(msm.linked, Patient)
		cat(paste('\nfound patients in denominator group with multiple sequences that can be taken as truly linked, n=',msm.linked[,length(unique(Patient))]))
		cat(paste('\nfound seq in denominator group with multiple sequences that can be taken as truly linked, n=',msm.linked[,length(unique(FASTASampleCode))]))
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
		file	<- paste(indir, '/', infiletree, '_', gsub('/',':',insignat),"_distPairs.R",sep='')
		load(file)	#loads brl.tpairs	
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
		#	get brl between true positives
		msm.linked			<- msm.linked[, 	list( brl= brl[ my.lower.tri.index(brl.n, sc.t, sc.i) ], Patient=Patient) ,by=c('FASTASampleCode','t.FASTASampleCode')]
		#
		#	compute branch lengths between infected and all potential transmitters
		if(1)
		{
			#	use median brl.tpairs
			df.tpairs.brl		<- brl.tpairs[, list(brl=median(brl)), by=c('FASTASampleCode','t.FASTASampleCode')]
			df.tpairs.brl		<- merge( subset(YX.tpairs, select=c(FASTASampleCode, t.FASTASampleCode)), df.tpairs.brl, by=c('FASTASampleCode','t.FASTASampleCode'), all.x=TRUE	)
			stopifnot(df.tpairs.brl[, !any(is.na(FASTASampleCode))])
			stopifnot(df.tpairs.brl[, !any(is.na(t.FASTASampleCode))])
			#	if any missing, use brl of tree with largest max likelihood
			df.brlmiss			<- subset(df.tpairs.brl, is.na(brl), select=c(FASTASampleCode, t.FASTASampleCode))
			cat(paste('\nWARNING: no median brl for pairs, using ML branch length. n=', nrow(df.brlmiss)))
			tmp					<- unique( subset(df.brlmiss, select=FASTASampleCode) )	
			tmp					<- tmp[, list(sc.i=match(FASTASampleCode, ph$tip.label)), by='FASTASampleCode']
			df.brlmiss			<- merge(df.brlmiss, tmp, by='FASTASampleCode')
			tmp					<- unique( subset(df.brlmiss, select=t.FASTASampleCode) )
			tmp					<- tmp[, list(sc.t=match(t.FASTASampleCode, ph$tip.label)), by='t.FASTASampleCode']
			df.brlmiss			<- merge(df.brlmiss, tmp, by='t.FASTASampleCode')
			tmp					<- df.brlmiss[, which(sc.t<sc.i)]
			tmp2				<- df.brlmiss[tmp, sc.t]
			set(df.brlmiss, tmp, 'sc.t', df.brlmiss[tmp,sc.i])
			set(df.brlmiss, tmp, 'sc.i', tmp2)		
			df.brlmiss			<- df.brlmiss[, list( brl= brl[ my.lower.tri.index(brl.n, sc.t, sc.i) ]) , by=c('FASTASampleCode','t.FASTASampleCode')]			
			df.tpairs.brl		<- rbind( subset(df.tpairs.brl, !is.na(brl)), df.brlmiss )			
		}
		if(0)
		{
			#	use brl of tree with largest max likelihood
			df.tpairs.brl		<- YX.tpairs
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
			df.tpairs.brl		<- subset(df.tpairs.brl, !is.na(brl))		#some requested brl may be missing because sequences have been deselected (too short or recombinants)			
		}
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
		require(grid)
		require(ggplot2)
		tmp		<- do.call('rbind', list(	data.table(brl= msm.linked[,brl], group='Sequence pairs\nfrom same host'), data.table(brl= msm.unlinked.bytime[,brl], group='Sequence pairs\nfrom patients unlinked by death'),
											data.table(brl= df.tpairs.brl[,brl], group='Clustering sequence pairs\nwith at least one sequence\nfrom a recently infected MSM')		))
		set(tmp, NULL, 'group', tmp[, factor(group, levels=c('Clustering sequence pairs\nwith at least one sequence\nfrom a recently infected MSM', 'Sequence pairs\nfrom same host', 'Sequence pairs\nfrom patients unlinked by death'))])
				
		ggplot(tmp, aes(x= brl, fill=group)) + labs(x="Substitutions / site", y='Sequence pairs') +
				scale_fill_brewer(palette='Set1') + theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(13,'mm')) +
				geom_histogram(binwidth=0.002)		
		cat(paste('\nsave to file',plot.file))
		ggsave(file=plot.file, w=6,h=6)		
		#tmp2	<- data.table( x= msm.unlinked.bytime[, rep(quantile(brl, p=c(1e-3, 1e-2)),2,each=TRUE)], y=rep(c(0,2e3),2), quantile= factor(rep(c(0.001, 0.01),2))) 	
		ggplot() + labs(x="Substitutions / site", y='Sequence pairs') + xlim(0, 0.15) +
				scale_fill_brewer(palette='Set1') + theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(13,'mm')) +			 
				geom_histogram(data=tmp, aes(x= brl, fill=group), binwidth=0.002)	
		ggsave(file=paste( substr(plot.file,1,nchar(plot.file)-3),'inset.pdf',sep='' ), w=6,h=6)
		
		#ggplot(tmp, aes(x= brl, fill=group)) + labs(x="Genetic distance", y='Sequence pairs') + xlim(0, 0.15) +
		#		scale_fill_brewer(palette='Set1') + theme(legend.justification=c(1,1), legend.position=c(1,1), legend.key.size=unit(13,'mm')) +
		#		geom_histogram(aes(y=..density..), colour="black", binwidth=0.002, position="identity", alpha=0.7)				
		#pdf(plot.file, w=5, h=5)
		#cat(paste('\nplot to file',plot.file))
		#require(RColorBrewer)
		#par(mar=c(4,4,0.5,0.5))
		#cols		<- brewer.pal(3, 'Set1')
		#legend.txt	<- c('same host', 'Died < last HIV- test', 'potential transmission pairs')	
		#xlim		<- c(0, max( msm.linked[, max(brl)], msm.unlinked.bytime[, max(brl)], df.tpairs.brl[, max(brl)])*1.1 )
		#tmp			<- seq(from=xlim[1], to=xlim[2], by=diff(xlim)/200)
		#hist( msm.linked[, brl], breaks=tmp , col=my.fade.col(cols[1],0.5), add=0, freq=0, xlab='branch length', main='' )
		#hist( msm.unlinked.bytime[, brl], breaks=tmp, col=my.fade.col(cols[2],0.5), freq=0, add=1 )
		#hist( df.tpairs.brl[, brl], breaks=tmp, col=my.fade.col(cols[3],0.5), freq=0, add=1 )
		#legend('topright', bty='n', fill= cols, border=NA, legend=legend.txt)		
		#tmp			<- msm.unlinked.bytime[, quantile(brl, p=c(1e-5, 1e-4, 5e-4, 0.001, 0.01))]
		#ltys		<- seq_along(tmp)+1
		#abline(v=tmp, col=cols[2], lty=ltys)
		#legend('bottomright', bty='n', lty= c(ltys,NA), col=cols[2], border=NA, legend=c('1e-5', '1e-4', '5e-4', '1e-3', '1e-2',''))
		#dev.off()
	}	
	list(tpairs=df.tpairs.brl, linked=msm.linked, unlinked=msm.unlinked.bytime)	
}
######################################################################################
#infilecov<- infile.cov.study; infile.viro<- infile.viro.study; infile.immu<- infile.immu.study; infile.treatment<- infile.treatment.study; adjust.AcuteByNegT=adjust.AcuteByNegT; adjust.NegT4Acute=NA; adjust.NegTByDetectability=0.25; adjust.minSCwindow=0.25; adjust.AcuteSelect=c('Yes','Maybe'); use.Acute_Spec=0; t.recent.endctime=t.recent.endctime; t.recent.startctime=t.recent.startctime; df.viro.part=NULL; df.immu.part=NULL; df.treatment.part=NULL; df.all.part=NULL
project.athena.Fisheretal.select.denominator<- function(indir, infile, insignat, indircov, infilecov, infile.viro, infile.immu, infile.treatment, 
															infiletree=NULL, adjust.AcuteByNegT=NA, adjust.NegT4Acute=NA, adjust.NegTByDetectability=NA, adjust.minSCwindow=NA, adjust.AcuteSelect=c('Yes'), use.AcuteSpec=0, t.recent.startctime=1996., t.recent.endctime=2013.,
															df.viro.part=NULL, df.immu.part=NULL, df.treatment.part=NULL, df.all.part=NULL)
{
	#	adjust.AcuteByNegT<- 0.75
	#	fixed input args in !is.na(infiletree)
	opt.brl			<- "dist.brl.casc" 
	thresh.brl		<- 0.096
	thresh.bs		<- 0.8	
	#
	#	load patient RNA
	#
	load(infile.viro)
	set(df, NULL, 'Patient', df[, as.character(Patient)])
	df.viro				<- df
	if(!is.null(df.viro.part))
	{
		df.viro				<- rbind(df.viro, df.viro.part, use.names=TRUE)
		setkey(df.viro, Patient, PosRNA)
		df.viro				<- unique(df.viro)	
	}	
	#
	#	load patient CD4
	#
	load(infile.immu)
	set(df, NULL, 'Patient', df[, as.character(Patient)])
	df.immu				<- df
	if(!is.null(df.immu.part))
	{
		df.immu				<- rbind(df.immu, df.immu.part, use.names=TRUE)
		setkey(df.immu, Patient, PosCD4)
		df.immu				<- unique(df.immu)
	}
	#
	#	load patient regimen
	#
	load(infile.treatment)
	set(df, NULL, 'Patient', df[, as.character(Patient)])
	df.treatment		<- df
	if(!is.null(df.treatment.part))
	{
		df.treatment		<- rbind(df.treatment, df.treatment.part, use.names=TRUE)
		setkey(df.treatment, Patient, StartTime)
		df.treatment		<- unique(df.treatment)
	}
	#
	# 	recent msm 
	#
	load(paste(indircov,'/',infilecov,'.R',sep=''))
	set(df.all, NULL, 'Patient', df.all[, as.character(Patient)])
	set(df.all, NULL, 'PosSeqT', hivc.db.Date2numeric(df.all[,PosSeqT]))
	set(df.all, NULL, 'DateBorn', hivc.db.Date2numeric(df.all[,DateBorn]))	
	set(df.all, NULL, 'NegT', hivc.db.Date2numeric(df.all[,NegT]))
	set(df.all, NULL, 'PosT', hivc.db.Date2numeric(df.all[,PosT]))
	set(df.all, NULL, 'DateDied', hivc.db.Date2numeric(df.all[,DateDied]))
	set(df.all, NULL, 'DateLastContact', hivc.db.Date2numeric(df.all[,DateLastContact]))
	set(df.all, NULL, 'DateFirstEverCDCC', hivc.db.Date2numeric(df.all[,DateFirstEverCDCC]))
	set(df.all, NULL, 'AnyPos_T1', hivc.db.Date2numeric(df.all[,AnyPos_T1]))
	set(df.all, NULL, 'PoslRNA_T1', hivc.db.Date2numeric(df.all[,PoslRNA_T1]))
	set(df.all, NULL, 'PoslRNA_TS', hivc.db.Date2numeric(df.all[,PoslRNA_TS]))
	set(df.all, NULL, 'lRNA.hb4tr_LT', hivc.db.Date2numeric(df.all[,lRNA.hb4tr_LT]))
	set(df.all, NULL, 'PosCD4_T1', hivc.db.Date2numeric(df.all[,PosCD4_T1]))
	set(df.all, NULL, 'PosCD4_TS', hivc.db.Date2numeric(df.all[,PosCD4_TS]))
	set(df.all, NULL, 'AnyT_T1', hivc.db.Date2numeric(df.all[,AnyT_T1]))	
	#
	#	adjust Acute=='Maybe' by NegT 
	if(!is.na(adjust.AcuteByNegT))
	{
		tmp		<- which( df.all[, (is.na(isAcute) | isAcute=='No') & !is.na(NegT) & AnyPos_T1<=NegT+adjust.AcuteByNegT])
		cat(paste('\ndf.all: set Acute==Maybe when NegT is close to AnyPos_T1, n=',length(tmp)))
		set(df.all, tmp, 'isAcute', 'Maybe')
	}	
	if(use.AcuteSpec)
	{
		cat(paste('\nUsing use.Acute_Spec=',use.AcuteSpec))
		setnames(df.all, c('isAcute','isAcuteNew'), c('isAcuteOld','isAcute'))
	}
	#	date all NegT back by adjust.NegTByDetectability to allow for infection that the test does not pick up 
	if(!is.na(adjust.NegTByDetectability))
	{
		tmp		<- df.all[, which(!is.na(NegT))]
		cat(paste('\ndate back NegT values to allow for undetected infection for n=',length(tmp)))
		set(df.all, tmp, 'NegT', df.all[tmp, NegT-adjust.NegTByDetectability])
	}
	#	adjust missing NegT when Acute=='Yes'
	if(!is.na(adjust.NegT4Acute))
	{		
		tmp		<- which( df.all[, !is.na(isAcute) & isAcute=='Yes' & is.na(NegT)] )	
		cat(paste('\nset NegT for Acute==Yes and NegT missing to one year before AnyPos_T1, n=',length(tmp)))
		set(df.all, tmp,'NegT', subset(df.all, !is.na(isAcute) & isAcute=='Yes' & is.na(NegT) )[, AnyPos_T1-adjust.NegT4Acute] )					
	}	
	#	make sure the SC interval is at least adjust.NegTByDetectability years
	if(!is.na(adjust.minSCwindow))
	{
		tmp		<- df.all[, which(!is.na(NegT) & AnyPos_T1-NegT<adjust.minSCwindow)]
		cat(paste('\nensure seroconversion interval is at least adjust.NegTByDetectability years, dating back NegT values for n=',length(tmp)))
		set(df.all, tmp, 'NegT', df.all[tmp, AnyPos_T1-adjust.minSCwindow])
	}
	#	merge with previous data entries if any
	if(!is.null(df.all.part))
	{
		df.all		<- merge(data.table(Patient=setdiff( df.all[, Patient], df.all.part[, unique(Patient)])), df.all, by='Patient')
		df.all		<- rbind(df.all.part, df.all, use.names=TRUE, fill=TRUE)		
	}	
	#
	#
	msm.recent		<- subset( df.all, Sex=='M' & !Trm%in%c('OTH','IDU','HET','BLOOD','PREG','HETfa','NEEACC','SXCH') )	
	msm.recent		<- subset(msm.recent, t.recent.startctime<=AnyPos_T1 & AnyPos_T1<t.recent.endctime)
	cat(paste('\nmsm after startctime',t.recent.startctime,'before endctime',t.recent.endctime,': #seq=',nrow(msm.recent),'#patient=',length(msm.recent[,unique(Patient)])))
	setkey(msm.recent, isAcute)
	msm.recent		<- subset(msm.recent, isAcute%in%adjust.AcuteSelect)	
	set(msm.recent, NULL, 'Trm', msm.recent[,factor(as.character(Trm))])
	cat(paste('\nmsm isAcute: #seq=',nrow(msm.recent),'#patient=',length(msm.recent[,unique(Patient)])))
	# 	recent msm seroconverters
	msm.recentsc	<- subset( msm.recent, !is.na(NegT))
	cat(paste('\nmsm isAcute & !is.na(NegT): #seq=',nrow(msm.recentsc),'#patient=',length(msm.recentsc[,unique(Patient)])))
	#
	#	load clustering msm
	#
	if(!is.null(infiletree))
	{
		argv			<<- hivc.cmd.clustering.msm(indir, infiletree, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=1)
		argv			<<- unlist(strsplit(argv,' '))		
		msm				<- hivc.prog.get.clustering.MSM()
		clumsm.info		<- msm$df.cluinfo
		#	update clumsm.info with latest values from df.all
		clumsm.info		<- merge( df.all, subset(clumsm.info, select=c(FASTASampleCode, cluster, Node, clu.npat, clu.ntip, clu.nFrgnInfection, clu.fPossAcute, clu.AnyPos_T1, clu.bwpat.medbrl)), by='FASTASampleCode')
		set(clumsm.info, NULL, 'clu.AnyPos_T1', hivc.db.Date2numeric(clumsm.info[,clu.AnyPos_T1]))
		#	fixup
		setkey(clumsm.info, Patient)
		tmp		<- clumsm.info[, list(ok= all(AnyPos_T1==min(AnyPos_T1))) ,by='Patient']
		for(x in subset(tmp, !ok)[, unique(Patient)])
		{
			z	<- clumsm.info[, which(Patient==x)]		
			set( clumsm.info, z, 'AnyPos_T1', clumsm.info[z, min(AnyPos_T1)])
		}
		# 	recently infected MSM in cluster
		setkey(clumsm.info, isAcute)			
		clumsm.recent	<- subset(clumsm.info, isAcute%in%adjust.AcuteSelect) 
		clumsm.recent	<- subset( clumsm.recent, Trm%in%c('MSM','BI') )
		set(clumsm.recent, NULL, 'Trm', clumsm.recent[,factor(as.character(Trm))])
		clumsm.recent		<- subset(clumsm.recent, t.recent.startctime<=AnyPos_T1 & AnyPos_T1<t.recent.endctime)
		cat(paste('\nmsm clustering after startctime',t.recent.startctime,'before endctime',t.recent.endctime,': #seq=',nrow(clumsm.recent),'#patient=',length(clumsm.recent[,unique(Patient)])))		
		# 	recent seroconcerverters in cluster
		clumsm.recentsc	<- subset( clumsm.recent, !is.na(NegT))		
		cat(paste('\nmsm clustering: #seq=',nrow(clumsm.info),'#patient=',length(clumsm.info[,unique(Patient)]),'#cluster=',length(clumsm.info[,unique(cluster)])))		
		cat(paste('\nmsm clustering isAcute: #seq=',nrow(clumsm.recent),'#patient=',length(clumsm.recent[,unique(Patient)]),'#cluster=',length(clumsm.recent[,unique(cluster)])))
		cat(paste('\nmsm clustering isAcute & !is.na(NegT): #seq=',nrow(clumsm.recentsc),'#patient=',length(clumsm.recentsc[,unique(Patient)]),'#cluster=',length(clumsm.recentsc[,unique(cluster)])))
		
		ans<- list(df.all=df.all, df.viro=df.viro, df.immu=df.immu, df.treatment=df.treatment, clumsm.info=clumsm.info, df.select=clumsm.recent, df.select.SEQ=msm.recent, clumsm.subtrees=msm$cluphy.subtrees, clumsm.ph=msm$cluphy)
	}
	else
	{
		ans<- list(df.all=df.all, df.viro=df.viro, df.immu=df.immu, df.treatment=df.treatment, df.select.SEQ=msm.recent )
	}
	ans	
}
######################################################################################
project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos<- function(clumsm.info, df.denom, any.pos.grace.yr= 0.5, select.if.transmitter.seq.unique=TRUE)
{
	df.transmitters	<- clumsm.info[J(df.denom[, unique(cluster)]),]		#anyone co-clustering with the recently infected in df.denom
	df.transmitters	<- subset(df.transmitters, select=c(cluster, Patient, AnyPos_T1, FASTASampleCode))		
	setnames(df.transmitters, colnames(df.transmitters), paste('t.',colnames(df.transmitters),sep=''))
	# 	all possible sequence pairs between recipient and co-clustering that meet diagnosis criterium
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
gamlss.centiles.get<- function(obj, xvar, cent = c(2.5, 97.5), with.ordering=TRUE, mu=NULL ) 
{
	fname 	<- obj$family[1]
	qfun 	<- paste("q", fname, sep = "")	
	if(with.ordering)
		xvaro	<- order(xvar)
	if(!with.ordering)
		xvaro	<- seq_along(xvar)
	if(is.null(mu))
		mu		<- fitted(obj, "mu")[xvaro]
	oxvar 	<- xvar[xvaro]	
	lpar 	<- length(obj$parameters)
	ii 		<- 0
	ans	<- data.table(x=oxvar)
	for (var in cent) {
		if (lpar == 1) {
			newcall <- call(qfun, var/100, mu = mu)
		}
		else if (lpar == 2) {
			newcall <- call(qfun, var/100, mu = mu, sigma = fitted(obj, "sigma")[xvaro])
		}
		else if (lpar == 3) {
			newcall <- call(qfun, var/100, mu = mu, sigma = fitted(obj, "sigma")[xvaro], nu = fitted(obj, "nu")[xvaro])
		}
		else {
			newcall <- call(qfun, var/100, mu = mu, sigma = fitted(obj, "sigma")[xvaro], nu = fitted(obj, "nu")[xvaro], tau = fitted(obj, "tau")[xvaro])
		}
		ii <- ii + 1
		set(ans, NULL, as.character(paste('q',var,sep='')), eval(newcall))						
	}
	ans
}
######################################################################################
project.athena.Fisheretal.t2inf<- function(indircov, infilecov, method.Acute='empirical', method.minQLowerU=0.01, adjust.AcuteByNegT=0.75, adjust.dt.CD4=1, adjust.AnyPos_y=2003, adjust.NegT=2, dur.AcuteYes= 365/2, dur.AcuteMaybe=320, use.AcuteSpec=0, t.recent.endctime=2011, plot.file=NULL)
{
	require(MASS)
	require(grid)
	require(reshape2)
	require(ggplot2)
	require(gamlss)	
	
	load(paste(indircov,'/',infilecov,'.R',sep=''))	
	#	adjust Acute=='Maybe' by NegT 
	if(!is.na(adjust.AcuteByNegT))
	{
		tmp		<- which( df.all[, (is.na(isAcute) | isAcute=='No') & !is.na(NegT) & AnyPos_T1<=NegT+adjust.AcuteByNegT*365])
		cat(paste('\nset Acute==Maybe when NegT is close to AnyPos_T1, n=',length(tmp)))
		set(df.all, tmp, 'isAcute', 'Maybe')
	}
	if(use.AcuteSpec)
	{
		cat(paste('\nUsing use.Acute_Spec=',use.AcuteSpec))
		setnames(df.all, c('isAcute','isAcuteNew'), c('isAcuteOld','isAcute'))
	}	
	tmp				<- subset(df.all, Sex=='M' & !Trm%in%c('OTH','IDU','HET','BLOOD','PREG','HETfa','NEEACC','SXCH') )
	df.sc.negT		<- subset(tmp, !is.na(NegT))
	setkey(df.sc.negT, Patient)
	df.sc.negT		<- unique(df.sc.negT)
	cat(paste('\nnumber of seroconverters with !is.na(NegT), n=',nrow(df.sc.negT)))
	
	tmp				<- df.sc.negT[, list(	dt.NegT=as.numeric(difftime(AnyPos_T1, NegT, units='days'))/365,
											mp.NegT=round(as.numeric(difftime(AnyPos_T1, NegT, units='days'))/adjust.NegT),
											dt.CD4=as.numeric(difftime(PosCD4_T1, AnyPos_T1, units='days'))/365,
											AnyPos_a=as.numeric(difftime(AnyPos_T1, DateBorn, units='days'))/365,
											AnyPos_y=hivc.db.Date2numeric(AnyPos_T1)), by='Patient']
	#	seroconverters
	df.sc.negT		<- merge(df.sc.negT, tmp, by='Patient')
	#	seroconverters with acute infection at diagnosis
	df.scay.cd4		<- subset( df.sc.negT, isAcute=='Yes' & mp.NegT>0 & AnyPos_y>adjust.AnyPos_y & AnyPos_y<t.recent.endctime & dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)
	#	seroconverters with acute infection at diagnosis
	df.scaym.cd4	<- subset( df.sc.negT, ( isAcute=='Maybe' ) & mp.NegT>0 & AnyPos_y>adjust.AnyPos_y & AnyPos_y<t.recent.endctime & dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)	
	#	seroconverters with chronic infection at diagnosis
	df.scchr.negT	<- subset( df.sc.negT, isAcute=='No' & AnyPos_y>adjust.AnyPos_y & AnyPos_y<t.recent.endctime)
	df.scana.negT	<- subset( df.sc.negT, is.na(isAcute) & AnyPos_y>adjust.AnyPos_y & AnyPos_y<t.recent.endctime)
	#	seroconverters with chronic infection at diagnosis and 1st CD4 count within 1 year after diagnosis and before ART
	df.scchr.cd4	<- subset( df.scchr.negT, dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)
	df.scana.cd4	<- subset( df.scana.negT, dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)
	#	no evidence from midpoints that calendar time matters
	#ggplot(subset(df.scchr.cd4, CD4_T1<250), aes(x=AnyPos_y, y=mp.NegT)) + geom_point() + geom_smooth()
	#ggplot(subset(df.scchr.cd4, CD4_T1>250 & CD4_T1<850), aes(x=AnyPos_y, y=mp.NegT)) + geom_point() + geom_smooth()	
	m4a	<- glm.nb(mp.NegT ~  1, data=subset(df.scchr.cd4, CD4_T1>850 & dt.NegT<4), link=identity, trace= 1, maxit= 50)	
	m4d	<- glm.nb(mp.NegT ~  AnyPos_a, data=subset(df.scchr.negT, AnyPos_a<=45), link=identity, trace= 1, maxit= 50)
	m4e	<- glm.nb(mp.NegT ~  AnyPos_a, data=subset(df.scana.negT, AnyPos_a<=45), link=identity, trace= 1, maxit= 50)	
	m4b	<- glm.nb(mp.NegT ~  AnyPos_a+0+offset(rep(m4d$coefficients["(Intercept)"],length(mp.NegT))), data=subset(df.scchr.cd4, CD4_T1<=850 & CD4_T1>250 & AnyPos_a<=45), link=identity, trace= 1, maxit= 50)
	m4c	<- glm.nb(mp.NegT ~  AnyPos_a+0+offset(rep(m4d$coefficients["(Intercept)"],length(mp.NegT))), data=subset(df.scchr.cd4, CD4_T1<250 & AnyPos_a<=45), link=identity, trace= 1, maxit= 50)
	#plots
	if(0 & !is.na(plot.file))
	{
		tmp	<- subset(df.sc.negT, AnyPos_y>1996 & dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)
		set(tmp, tmp[, which(is.na(isAcute))],'isAcute','not reported')
		set(tmp, tmp[, which(isAcute=='No')],'isAcute','Chronic')
		set(tmp, tmp[, which(isAcute=='Maybe')],'isAcute','Recent')
		set(tmp, tmp[, which(isAcute=='Yes')],'isAcute','Recent')
		set(tmp, NULL, 'isAcute', tmp[, factor(isAcute, levels=c('Recent','Chronic','not reported'))] )
		set(tmp, NULL, 'AnyPos_y', factor(tmp[, round(AnyPos_y, d=0)]))
		ggplot(tmp, aes(x=AnyPos_y, fill=isAcute)) + geom_bar(binwidth=1.5, position='dodge') +
				scale_y_continuous(breaks=seq(0,500,20)) +
				labs(x="year of diagnosis", y='HIV seroconverters among MSM\nwith first CD4 count within 12 mo of diagnosis\nand before ART start\n(total)') +
				scale_fill_brewer(palette='Set2',name='Stage of HIV infection\n at diagnosis')
		ggsave(paste(plot.file,'sc_numbers.pdf'), w=15, h=6)		
				
	}
	if(0 & !is.na(plot.file))
	{
		tmp	<- rbind(df.scay.cd4, df.scaym.cd4)
		tmp[, mpy.NegT:= mp.NegT/365.25]
		set(tmp, NULL, 'CD4_T1c', tmp[, factor(cut(CD4_T1, breaks=c(-Inf, 250, 850, Inf), labels=c('1st CD4 count\nwithin 12 months after diagnosis\n< 250','1st CD4 count\nwithin 12 months after diagnosis\nin [250-850]','1st CD4 count\nwithin 12 months after diagnosis\n> 850')))])
		setkey(tmp, CD4_T1c, AnyPos_a)
		tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)+CD4_T1c'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)+CD4_T1c'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
		tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
		tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
		tmp		<- cbind(tmp, tmp3)
		ggplot(tmp, aes(x=AnyPos_a, y=mpy.NegT)) + 
				geom_point(size=0.85) + scale_y_continuous(limits=c(0,11.4), breaks=seq(0,14,2)) + xlim(18, 61) +			
				geom_line(aes(y=y.b), colour='blue') +
				geom_ribbon(aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) +
				labs(x="age at diagnosis", y='time between midpoint of\nseroconversion interval and diagnosis\n(years)') +
				facet_grid(. ~ CD4_T1c) + theme_bw()
		ggsave(paste(plot.file,'scay_cd4_age.pdf'), w=15, h=6)
	}
	if(0 & !is.na(plot.file))
	{
		df.scaym.cd4[, mpy.NegT:= mp.NegT/365.25]
		set(df.scaym.cd4, NULL, 'CD4_T1c', df.scaym.cd4[, factor(cut(CD4_T1, breaks=c(-Inf, 250, 850, Inf), labels=c('1st CD4 count\nwithin 12 months after diagnosis\n< 250','1st CD4 count\nwithin 12 months after diagnosis\nin [250-850]','1st CD4 count\nwithin 12 months after diagnosis\n> 850')))])
		setkey(df.scaym.cd4, CD4_T1c, AnyPos_a)
		tmp	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)+CD4_T1c'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)+CD4_T1c'), data=as.data.frame(subset(df.scaym.cd4, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
		df.scaym.cd4[, y.b:= predict(tmp, type='response', se.fit=FALSE)]	
		tmp	<- gamlss.centiles.get(tmp, df.scaym.cd4$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
		df.scaym.cd4	<- cbind(df.scaym.cd4, tmp)
		ggplot(df.scaym.cd4, aes(x=AnyPos_a, y=mpy.NegT)) + 
				geom_point() + scale_y_continuous(limits=c(0,11.4), breaks=seq(0,14,2)) + xlim(18, 61) +			
				geom_line(aes(y=y.b), colour='blue') +
				geom_ribbon(aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) +
				labs(x="age at diagnosis", y='Midpoint of seroconversion interval\n(years from diagnosis)') +
				facet_grid(. ~ CD4_T1c) 
		ggsave(paste(plot.file,'scam_cd4_age.pdf'), w=10, h=6)
	}	
	if(0 & !is.na(plot.file))
	{
		df.scchr.cd4[, mpy.NegT:= mp.NegT/365.25]
		tmp	<- gamlss(as.formula('mpy.NegT ~ bs(CD4_T1, degree=7)'), sigma.formula=as.formula('~ bs(CD4_T1, degree=3)'), data=as.data.frame(subset(df.scchr.cd4, select=c(mpy.NegT, CD4_T1))), family=GA)		
		df.scchr.cd4[, y.b:= predict(tmp, type='response', se.fit=FALSE)]	
		tmp	<- gamlss.centiles.get(tmp, df.scchr.cd4$CD4_T1, cent = c(2.5, 97.5) )
		df.scchr.cd4	<- cbind(df.scchr.cd4, tmp)
		ggplot(df.scchr.cd4, aes(x=CD4_T1, y=mpy.NegT)) + 
				geom_point() + scale_y_continuous(limits=c(0,11), breaks=seq(0,10,2)) +			
				geom_line(aes(y=y.b), colour='blue') + theme_bw() +
				#geom_ribbon(aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) + 
				labs(x="1st CD4 count after diagnosis", y='Midpoint of seroconversion interval\n(years from diagnosis)')		
		ggsave(paste(plot.file,'scchr_cd4.pdf'), w=5, h=5)
		df.scchr.cd4[, y.b:=NULL]
		df.scchr.cd4[, x:=NULL]
		df.scchr.cd4[, q2.5:=NULL]
		df.scchr.cd4[, q97.5:=NULL]
	}
	if(0)
	{
		df.scana.cd4[, mpy.NegT:= mp.NegT/365.25]
		set(df.scana.cd4, NULL, 'CD4_T1c', df.scana.cd4[, factor(cut(CD4_T1, breaks=c(-Inf, 250, 850, Inf), labels=c('1st CD4 count\nwithin 12 months after diagnosis\n< 250','1st CD4 count\nwithin 12 months after diagnosis\nin [250-850]','1st CD4 count\nwithin 12 months after diagnosis\n> 850')))])
		setkey(df.scana.cd4, CD4_T1c, AnyPos_a)		
		tmp	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=3)+CD4_T1c'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)+CD4_T1c'), data=as.data.frame(subset(df.scana.cd4, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)
		#Rsq(tmp)	20.17%
		df.scana.cd4[, y.b:= predict(tmp, type='response', se.fit=FALSE)]
		tmp	<- gamlss.centiles.get(tmp, df.scana.cd4$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
		df.scana.cd4	<- cbind(df.scana.cd4, tmp)
		ggplot(df.scana.cd4, aes(x=AnyPos_a, y=mpy.NegT)) + 
				geom_point(size=1) + scale_y_continuous(limits=c(0,11.4), breaks=seq(0,14,2)) + xlim(18, 61) +			
				geom_line(aes(y=y.b), colour='blue') + theme_bw() +
				#geom_ribbon(aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) + 
				labs(x="age at diagnosis", y='time between midpoint of\nseroconversion interval and diagnosis\n(years)') +
				facet_grid(. ~ CD4_T1c)
		ggsave(paste(plot.file,'scana_cd4_age.pdf'), w=15, h=6)		
		df.scana.cd4[, y.b:=NULL]
		df.scana.cd4[, x:=NULL]
		df.scana.cd4[, q2.5:=NULL]
		df.scana.cd4[, q97.5:=NULL]		
	}
	if(!is.na(plot.file))
	{	
		df.scchr.cd4[, mpy.NegT:= mp.NegT/365.25]
		set(df.scchr.cd4, NULL, 'CD4_T1c', df.scchr.cd4[, factor(cut(CD4_T1, breaks=c(-Inf, 250, 850, Inf), labels=c('1st CD4 count\nwithin 12 months after diagnosis\n< 250','1st CD4 count\nwithin 12 months after diagnosis\nin [250-850]','1st CD4 count\nwithin 12 months after diagnosis\n> 850')))])
		setkey(df.scchr.cd4, CD4_T1c, AnyPos_a)
		tmp	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)+CD4_T1c'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)+CD4_T1c'), data=as.data.frame(subset(df.scchr.cd4, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)
		#Rsq(tmp)	20.10%
		df.scchr.cd4[, y.b:= predict(tmp, type='response', se.fit=FALSE)]	
		tmp	<- gamlss.centiles.get(tmp, df.scchr.cd4$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )		
		df.scchr.cd4	<- cbind(df.scchr.cd4, tmp)
		ggplot(df.scchr.cd4, aes(x=AnyPos_a, y=mpy.NegT)) + 
				geom_point(size=1) + scale_y_continuous(limits=c(0,11.4), breaks=seq(0,14,2)) + xlim(18, 61) +			
				geom_line(aes(y=y.b), colour='blue') + theme_bw() +
				#geom_ribbon(aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) + 
				labs(x="age at diagnosis", y='time between midpoint of\nseroconversion interval and diagnosis\n(years)') +
				facet_grid(. ~ CD4_T1c)
		ggsave(paste(plot.file,'scchr_cd4_age.pdf'), w=15, h=6)
		df.scchr.cd4[, y.b:=NULL]
		df.scchr.cd4[, x:=NULL]
		df.scchr.cd4[, q2.5:=NULL]
		df.scchr.cd4[, q97.5:=NULL]	
	}
	if(!is.na(plot.file))
	{
		df.scchr.negT[, mpy.NegT:= mp.NegT/365.25]
		tmp	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=7)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(df.scchr.negT, select=c(mpy.NegT, AnyPos_a))), family=GA)		 
		df.scchr.negT[, y.b:= predict(tmp, type='response', se.fit=FALSE)]
		tmp	<- gamlss.centiles.get(tmp, df.scchr.negT$AnyPos_a, cent = c(2.5, 97.5) )
		df.scchr.negT	<- cbind(df.scchr.negT, tmp)
		ggplot(df.scchr.negT, aes(x=AnyPos_a, y=mpy.NegT)) + xlim(18, 61) + scale_y_continuous(limits=c(0,11), breaks=seq(0,10,2)) + 
				geom_point() + geom_line(aes(y=y.b), colour='blue') + geom_ribbon(aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2) +
				labs(x="Age at diagnosis", y='Midpoint of seroconversion interval\n(years from diagnosis)')
		ggsave(paste(plot.file,'scchr_age.pdf'), w=5, h=5)
		df.scchr.negT[, y.b:=NULL]
		df.scchr.negT[, x:=NULL]
		df.scchr.negT[, q2.5:=NULL]
		df.scchr.negT[, q97.5:=NULL]
	}
	#		
	cat(paste('\navg time of acute.Yes phase=',dur.AcuteYes))
	cat(paste('\navg time of acute.Maybe phase=',dur.AcuteMaybe))	
	#
	t2inf.args						<- list()
	t2inf.args$cd4g850.intercept	<- m4a$coefficients["(Intercept)"]
	t2inf.args$cd4g850.theta		<- m4a$theta
	t2inf.args$cd4l850.intercept	<- m4d$coefficients["(Intercept)"]
	t2inf.args$cd4l850.coef			<- m4b$coefficients["AnyPos_a"]
	t2inf.args$cd4l850.theta		<- m4b$theta
	t2inf.args$cd4l250.intercept	<- m4d$coefficients["(Intercept)"]
	t2inf.args$cd4l250.coef			<- m4c$coefficients["AnyPos_a"]
	t2inf.args$cd4l250.theta		<- m4c$theta
	t2inf.args$cd4other.intercept	<- m4d$coefficients["(Intercept)"]
	t2inf.args$cd4other.coef		<- m4d$coefficients["AnyPos_a"]
	t2inf.args$cd4other.theta		<- m4d$theta	
	t2inf.args$method.minQLowerU	<- method.minQLowerU
	t2inf.args$dur.AcuteYes			<- dur.AcuteYes
	t2inf.args$dur.AcuteMaybe		<- dur.AcuteMaybe
	#
	if(method.Acute%in%c('lower','central','higher') && t2inf.args$method.minQLowerU==0.)
	{
		#recent stages are calibrated to 90% quantile, so assign NA afterwards
		#assign NA based on quantile 'method.minQLowerU' only to non-recent
		#if method.minQLowerU is set to 0, cut undiagnosed window at MEAN
		predict.t2inf	<- function(q, df, t2inf.args)
		{				
			df[, {
						print(Patient)
						print(CD4_T1)
						print(AnyPos_a)
						print(isAcute)
						if(!is.na(isAcute) & isAcute=='Yes') 
						{
							ans				<- pexp(q, 1/t2inf.args$dur.AcuteYes, lower.tail=FALSE)
							ans[ans<0.1]	<- NA_real_		
						}
						else if(!is.na(isAcute) & isAcute=='Maybe')
						{
							ans				<- pexp(q, 1/t2inf.args$dur.AcuteMaybe, lower.tail=FALSE)
							ans[ans<0.1]	<- NA_real_	
						}
						else if(is.na(PosCD4_T1) | (!is.na(AnyT_T1) & (AnyT_T1<PosCD4_T1 | PosCD4_T1>AnyPos_T1+1)))		#chronic or missing isAcute, first CD4 after ART start or first CD4 too far from PosDiag
						{
							ans				<- pnbinom(q, mu= t2inf.args$cd4other.intercept+t2inf.args$cd4other.coef*min(AnyPos_a,45), size=t2inf.args$cd4other.theta, lower.tail=FALSE)
							ans[ q>t2inf.args$cd4other.intercept+t2inf.args$cd4other.coef*min(AnyPos_a,45) ]<- NA_real_
						}
						else if(CD4_T1<250)		#chronic or missing isAcute, first CD4 before ART start
						{
							ans				<- pnbinom(q, mu= t2inf.args$cd4l250.intercept+t2inf.args$cd4l250.coef*min(AnyPos_a,45), size=t2inf.args$cd4l250.theta, lower.tail=FALSE)
							print(CD4_T1)
							print(AnyPos_a)
							print( t2inf.args$cd4l250.intercept+t2inf.args$cd4l250.coef*min(AnyPos_a,45) )
							ans[ q>t2inf.args$cd4l250.intercept+t2inf.args$cd4l250.coef*min(AnyPos_a,45) ]<- NA_real_
						}
						else if(CD4_T1<850)		#chronic or missing isAcute, first CD4 before ART start
						{
							ans				<- pnbinom(q, mu= t2inf.args$cd4l850.intercept+t2inf.args$cd4l850.coef*min(AnyPos_a,45), size=t2inf.args$cd4l850.theta, lower.tail=FALSE)
							ans[ q>t2inf.args$cd4l850.intercept+t2inf.args$cd4l850.coef*min(AnyPos_a,45) ]<- NA_real_
						}
						else					#chronic or missing isAcute, first CD4 before ART start
						{
							ans				<- pnbinom(q, mu= t2inf.args$cd4g850.intercept, size=t2inf.args$cd4g850.theta, lower.tail=FALSE)
							ans[ q>t2inf.args$cd4g850.intercept ]<- NA_real_							
						}		
						list(q=q, score=ans)
					}, by='Patient']			
		}
	}
	if(method.Acute%in%c('lower','central','higher') && t2inf.args$method.minQLowerU>0.)
	{
		#recent stages are calibrated to 90% quantile, so assign NA afterwards
		#assign NA based on quantile 'method.minQLowerU' only to non-recent
		predict.t2inf	<- function(q, df, t2inf.args)
		{				
			df[, {
						if(!is.na(isAcute) & isAcute=='Yes') 
						{
							ans				<- pexp(q, 1/t2inf.args$dur.AcuteYes, lower.tail=FALSE)
							ans[ans<0.1]	<- NA_real_		
						}
						else if(!is.na(isAcute) & isAcute=='Maybe')
						{
							ans				<- pexp(q, 1/t2inf.args$dur.AcuteMaybe, lower.tail=FALSE)
							ans[ans<0.1]	<- NA_real_	
						}
						else if(is.na(PosCD4_T1) | (!is.na(AnyT_T1) & (AnyT_T1<PosCD4_T1 | PosCD4_T1>AnyPos_T1+1)))		#chronic or missing isAcute, first CD4 after ART start or first CD4 too far from PosDiag
						{
							ans				<- pnbinom(q, mu= t2inf.args$cd4other.intercept+t2inf.args$cd4other.coef*min(AnyPos_a,45), size=t2inf.args$cd4other.theta, lower.tail=FALSE)
							ans[ans<t2inf.args$method.minQLowerU]	<- NA_real_
						}
						else if(CD4_T1<250)		#chronic or missing isAcute, first CD4 before ART start
						{
							ans				<- pnbinom(q, mu= t2inf.args$cd4l250.intercept+t2inf.args$cd4l250.coef*min(AnyPos_a,45), size=t2inf.args$cd4l250.theta, lower.tail=FALSE)
							ans[ans<t2inf.args$method.minQLowerU]	<- NA_real_
						}
						else if(CD4_T1<850)		#chronic or missing isAcute, first CD4 before ART start
						{
							ans				<- pnbinom(q, mu= t2inf.args$cd4l850.intercept+t2inf.args$cd4l850.coef*min(AnyPos_a,45), size=t2inf.args$cd4l850.theta, lower.tail=FALSE)
							ans[ans<t2inf.args$method.minQLowerU]	<- NA_real_
						}
						else					#chronic or missing isAcute, first CD4 before ART start
						{
							ans				<- pnbinom(q, mu= t2inf.args$cd4g850.intercept, size=t2inf.args$cd4g850.theta, lower.tail=FALSE)
							ans[ans<t2inf.args$method.minQLowerU]	<- NA_real_
						}		
						list(q=q, score=ans)
					}, by='Patient']			
		}
	}
	if(!method.Acute%in%c('lower','central','higher') )
	{
		#assign NA based on quantile 'method.minQLowerU' to all
		predict.t2inf	<- function(q, df, t2inf.args)
		{				
			df[, {
						if(!is.na(isAcute) & isAcute=='Yes') 
						{
							ans	<- pexp(q, 1/t2inf.args$dur.AcuteYes, lower.tail=FALSE)
						}
						else if(!is.na(isAcute) & isAcute=='Maybe')
						{
							ans	<- pexp(q, 1/t2inf.args$dur.AcuteMaybe, lower.tail=FALSE)
						}
						else if(is.na(PosCD4_T1) | (!is.na(AnyT_T1) & (AnyT_T1<PosCD4_T1 | PosCD4_T1>AnyPos_T1+1)))		#chronic or missing isAcute, first CD4 after ART start or first CD4 too far from PosDiag
						{
							ans	<- pnbinom(q, mu= t2inf.args$cd4other.intercept+t2inf.args$cd4other.coef*min(AnyPos_a,45), size=t2inf.args$cd4other.theta, lower.tail=FALSE)
						}
						else if(CD4_T1<250)		#chronic or missing isAcute, first CD4 before ART start
						{
							ans	<- pnbinom(q, mu= t2inf.args$cd4l250.intercept+t2inf.args$cd4l250.coef*min(AnyPos_a,45), size=t2inf.args$cd4l250.theta, lower.tail=FALSE)
						}
						else if(CD4_T1<850)		#chronic or missing isAcute, first CD4 before ART start
						{
							ans	<- pnbinom(q, mu= t2inf.args$cd4l850.intercept+t2inf.args$cd4l850.coef*min(AnyPos_a,45), size=t2inf.args$cd4l850.theta, lower.tail=FALSE)
						}
						else					#chronic or missing isAcute, first CD4 before ART start
						{
							ans	<- pnbinom(q, mu= t2inf.args$cd4g850.intercept, size=t2inf.args$cd4g850.theta, lower.tail=FALSE)
						}		
						ans[ans<t2inf.args$method.minQLowerU]	<- NA_real_
						list(q=q, score=ans)
					}, by='Patient']			
		}
	}
	if(!is.na(plot.file))
	{		
		t.period	<- 1/64
		b4care		<- do.call('rbind', list(	subset(df.scchr.cd4, AnyPos_a>25.0 & AnyPos_a<25.99 & CD4_T1<250, select=c(Patient, DateBorn, AnyPos_T1, isAcute, NegT, PosCD4_T1, CD4_T1,  AnyT_T1))[1,],
												subset(df.scchr.cd4, AnyPos_a>35.0 & AnyPos_a<35.99 & CD4_T1<250, select=c(Patient, DateBorn, AnyPos_T1, isAcute, NegT, PosCD4_T1, CD4_T1,  AnyT_T1))[1,],
												subset(df.scchr.cd4, AnyPos_a>35.0 & AnyPos_a<35.99 & CD4_T1>250 & CD4_T1<850, select=c(Patient, DateBorn, AnyPos_T1, isAcute, NegT, PosCD4_T1, CD4_T1,  AnyT_T1))[1,],
												subset(df.scchr.cd4, AnyPos_a>34.5 & AnyPos_a<35.99 & CD4_T1>850, select=c(Patient, DateBorn, AnyPos_T1, isAcute, NegT, PosCD4_T1, CD4_T1,  AnyT_T1))[1,]												
												))
		set(b4care, NULL, 'DateBorn', hivc.db.Date2numeric(b4care[,DateBorn]))
		set(b4care, NULL, 'AnyPos_T1', hivc.db.Date2numeric(b4care[,AnyPos_T1]))
		set(b4care, NULL, 'NegT', NA_real_)
		set(b4care, NULL, 'PosCD4_T1', hivc.db.Date2numeric(b4care[,PosCD4_T1]))
		set(b4care, NULL, 'AnyT_T1', hivc.db.Date2numeric(b4care[,AnyT_T1]))
		tmp			<- b4care[, list(AnyPos_a=AnyPos_T1-DateBorn, ts=ifelse(is.na(NegT), AnyPos_T1-10, NegT), te=AnyPos_T1), by='Patient']
		b4care		<- merge(b4care, tmp, by='Patient')
		set(b4care, NULL, 'ts', b4care[, floor(ts) + round( (ts%%1)*100 %/% (t.period*100) ) * t.period ] )
		set(b4care, NULL, 'te', b4care[, floor(te) + round( (te%%1)*100 %/% (t.period*100) ) * t.period - t.period] )	#last time undiagnosed is just before first time period diagnosed		
		tmp			<- predict.t2inf(seq(0,10,t.period)*365, b4care, t2inf.args)	
		tmp			<- merge( subset( tmp, !is.na(score) ), subset(b4care, select=c(Patient,isAcute,ts,te)), by='Patient' )
		set(tmp, NULL, 'q', tmp[,q/365])
		tmp[, label:='']
		labels		<- c(	'chronic HIV infection at diagnosis,\n CD4 < 250,\n 25 years at diagnosis',
							'chronic HIV infection at diagnosis,\n CD4 < 250,\n 35 years at diagnosis',
							'chronic HIV infection at diagnosis,\n CD4 in [250-850],\n 35 years at diagnosis',
							'chronic HIV infection at diagnosis,\n CD4 > 850,\n 35 years at diagnosis')
		set(tmp, tmp[, which(Patient=='M31752')], 'label', labels[2])
		set(tmp, tmp[, which(Patient=='M32501')], 'label', labels[3])
		set(tmp, tmp[, which(Patient=='M34515')], 'label', labels[1])
		set(tmp, tmp[, which(Patient=='M39133')], 'label', labels[4])
		set(tmp, NULL, 'label', tmp[, factor(label, levels=labels)] )
		setkey(tmp, label)
		ggplot(tmp, aes(x=q, y=score, group=label, colour=label)) + geom_line() + theme_bw() + theme(legend.key.size=unit(13,'mm'))  +
				scale_x_continuous(breaks=seq(0,10,1)) +
				labs(x="t\n(years)", y='Probability that\nsurrogate time between HIV infection and diagnosis is > t') + 
				scale_colour_brewer(palette='Set2',name='Individual with')
		ggsave(paste(plot.file,'sc_predict.pdf'), w=15, h=7)
	}
	#
	list(predict.t2inf=predict.t2inf, t2inf.args=t2inf.args)
}
######################################################################################
project.athena.Fisheretal.t2inf.explore.acute<- function(indircov, infilecov)
{
	load(paste(indircov,'/',infilecov,'.R',sep=''))
	adjust.sc.max= 3
	adjust.AcuteByNegT=0.75
	
	if(!is.na(adjust.AcuteByNegT))
	{
		tmp		<- which( df.all[, (is.na(isAcute) | isAcute=='No') & !is.na(NegT) & AnyPos_T1<=NegT+adjust.AcuteByNegT*365])
		cat(paste('\nset Acute==Maybe when NegT is close to AnyPos_T1, n=',length(tmp)))
		set(df.all, tmp, 'isAcute', 'Maybe')
	}
	tmp				<- subset(df.all, Sex=='M' & !Trm%in%c('OTH','IDU','HET','BLOOD','PREG','HETfa','NEEACC','SXCH') )
	df.sc.negT		<- subset(tmp, !is.na(NegT))
	setkey(df.sc.negT, Patient)
	df.sc.negT		<- unique(df.sc.negT)
	cat(paste('\nnumber of seroconverters with !is.na(NegT), n=',nrow(df.sc.negT)))
	
	tmp				<- df.sc.negT[, list(	dt.NegT=as.numeric(difftime(AnyPos_T1, NegT, units='days'))/365,
					mp.NegT=round(as.numeric(difftime(AnyPos_T1, NegT, units='days'))/2),
					dt.CD4=as.numeric(difftime(PosCD4_T1, AnyPos_T1, units='days'))/365,
					AnyPos_a=as.numeric(difftime(AnyPos_T1, DateBorn, units='days'))/365,
					AnyPos_y=hivc.db.Date2numeric(AnyPos_T1)), by='Patient']
	
	df.sc.negT		<- merge(df.sc.negT, tmp, by='Patient')
	hist(df.sc.negT[,dt.NegT], breaks=100)
	tmp				<- subset(df.sc.negT, NegT_Acc=='Yes' )
	cat(paste('\nnumber of seroconverters with accurate NegT, n=',nrow(tmp)))
	hist(tmp[,dt.NegT], breaks=100)
	tmp				<- subset(df.sc.negT, isAcute=='No' )
	cat(paste('\nnumber of seroconverters with chronic, n=',nrow(tmp)))
	hist(tmp[,dt.NegT], breaks=100)
	tmp				<- subset(df.sc.negT, NegT_Acc=='Yes' & isAcute=='No' )
	cat(paste('\nnumber of seroconverters with chronic and accurate NegT, n=',nrow(tmp)))
	hist(tmp[,dt.NegT], breaks=100)
	
	#
	#	Acute==Yes
	#
	#	CALENDAR YEAR	
	df.sca.negT	<- subset( df.sc.negT, isAcute=='Yes')	
	hist(df.sca.negT[, NegT], breaks=20)
	plot(df.sca.negT[, AnyPos_y], df.sca.negT[, NegT], pch=18)
	plot(df.sca.negT[, AnyPos_y], df.sca.negT[, dt.NegT], pch=18)
	
	#	restrict to >2003
	df.sca.negT	<- subset( df.sc.negT, isAcute=='Yes' & AnyPos_y>2003)	
	plot(df.sca.negT[, AnyPos_y], df.sca.negT[, dt.NegT], pch=18)

	mean( subset(df.sca.negT, dt.NegT<10)[, dt.NegT] )
	#	compare to Exp(1)
	tmp	<- subset(df.sca.negT, dt.NegT>=0)
	hist( tmp[, dt.NegT], breaks=seq(-0.0002,30,0.2), xlim=c(0,10), freq=0)
	lines( seq(0.0002,30,0.02), dexp(seq(0.0002,30,0.02),rate=1), col='red')
	lines( seq(0.0002,30,0.02), dexp(seq(0.0002,30,0.02),rate=1/0.75), col='blue')
	lines( seq(0.0002,30,0.02), dexp(seq(0.0002,30,0.02),rate=1/0.5), col='green')
	
	tmp	<- subset(df.sca.negT, dt.NegT>=0)
	hist( tmp[, mp.NegT], breaks=seq(0,365*30,50), xlim=c(0,365*4), freq=0)
	lines( seq(0,365*30,50), dexp(seq(0,365*30,50),rate=1/365), col='red')
	lines( seq(0,365*30,50), dexp(seq(0,365*30,50),rate=1/(365/2)), col='blue')
	lines( seq(0,365*30,50), dexp(seq(0,365*30,50),rate=1/(365/3)), col='green')
	#
	#	--> use 6 mo for simplicity
	#
	hist(subset(df.sca.negT, dt.CD4<1 & PosCD4_T1<AnyT_T1)[, CD4_T1], breaks=50)
	
	df.sca.negT[, CD4.col:='black']
	set(df.sca.negT, df.sca.negT[, which(dt.CD4<1 & PosCD4_T1<AnyT_T1 & CD4_T1<250)], 'CD4.col', 'red')
	set(df.sca.negT, df.sca.negT[, which(dt.CD4<1 & PosCD4_T1<AnyT_T1 & CD4_T1>=250 & CD4_T1<850)], 'CD4.col', 'blue')
	set(df.sca.negT, df.sca.negT[, which(dt.CD4<1 & PosCD4_T1<AnyT_T1 & CD4_T1>850)], 'CD4.col', 'green')
	plot(df.sca.negT[, AnyPos_a], df.sca.negT[, dt.NegT], pch=18, col=df.scchr.negT[,CD4.col], ylim=c(0,5))
	plot(df.sca.negT[, CD4_T1], df.sca.negT[, dt.NegT], pch=18 )	
	#
	#	Acute==Maybe
	#
	#	CALENDAR YEAR	
	df.scam.negT	<- subset( df.sc.negT, isAcute=='Maybe')	
	hist(df.scam.negT[, NegT], breaks=20)
	plot(df.scam.negT[, AnyPos_y], df.scam.negT[, NegT], pch=18)
	plot(df.scam.negT[, AnyPos_y], df.scam.negT[, dt.NegT], pch=18)
	df.scam.negT	<- subset( df.sc.negT, isAcute=='Maybe' & AnyPos_y>2003)
	#
	tmp	<- subset(df.scam.negT, dt.NegT>=0)
	hist( tmp[, mp.NegT], breaks=seq(0,365*30,50), xlim=c(0,365*4), freq=0)
	lines( seq(0,365*30,50), dexp(seq(0,365*30,50),rate=1/365), col='red')
	lines( seq(0,365*30,50), dexp(seq(0,365*30,50),rate=1/(365/2)), col='blue')
	lines( seq(0,365*30,50), dexp(seq(0,365*30,50),rate=1/(365/3)), col='green')
	tmp[, mean(mp.NegT)]
	#	almost twice as long time to midpoint
	plot(df.scam.negT[, AnyPos_y], df.scam.negT[, dt.NegT], pch=18)

	df.scam.negT[, CD4.col:='black']
	set(df.scam.negT, df.scam.negT[, which(dt.CD4<1 & PosCD4_T1<AnyT_T1 & CD4_T1<250)], 'CD4.col', 'red')
	set(df.scam.negT, df.scam.negT[, which(dt.CD4<1 & PosCD4_T1<AnyT_T1 & CD4_T1>=250 & CD4_T1<850)], 'CD4.col', 'blue')
	set(df.scam.negT, df.scam.negT[, which(dt.CD4<1 & PosCD4_T1<AnyT_T1 & CD4_T1>850)], 'CD4.col', 'green')
	
	# 	no dependence on age
	plot(df.scam.negT[, AnyPos_a], df.scam.negT[, dt.NegT], pch=18, col=df.scchr.negT[,CD4.col], ylim=c(0,5))			
	tmp				<- loess(dt.NegT ~ AnyPos_a, df.scam.negT, span=0.5)
	tmp				<- predict(tmp, data.frame(AnyPos_a = seq(20, 70, 1)), se=TRUE)
	lines(seq(20, 70, 1), tmp$fit, col='red')
	lines(seq(20, 70, 1), tmp$fit+2*tmp$se.fit, col='red')
	lines(seq(20, 70, 1), tmp$fit-2*tmp$se.fit, col='red')
	
	#	CD4 count	-- const below 850 and lower const above 850 
	plot( df.scam.negT[, CD4_T1], df.scam.negT[, dt.NegT], pch=18)
	tmp				<- loess(dt.NegT ~ CD4_T1, df.scam.negT, span=0.5)
	tmp				<- predict(tmp, data.frame(CD4_T1 = seq(20, 1100, 10)), se=TRUE)
	lines(seq(20, 1100, 10), tmp$fit, col='red')
	lines(seq(20, 1100, 10), tmp$fit+2*tmp$se.fit, col='red')
	lines(seq(20, 1100, 10), tmp$fit-2*tmp$se.fit, col='red')

}
######################################################################################
project.athena.Fisheretal.t2inf.explore.chronic<- function(indircov, infilecov)
{
	load(paste(indircov,'/',infilecov,'.R',sep=''))
	adjust.sc.max= 3
	adjust.AcuteByNegT=0.75
	#	adjust Acute=='Maybe' by NegT 
	if(!is.na(adjust.AcuteByNegT))
	{
		tmp		<- which( df.all[, (is.na(isAcute) | isAcute=='No') & !is.na(NegT) & AnyPos_T1<=NegT+adjust.AcuteByNegT*365])
		cat(paste('\nset Acute==Maybe when NegT is close to AnyPos_T1, n=',length(tmp)))
		set(df.all, tmp, 'isAcute', 'Maybe')
	}
	tmp				<- subset(df.all, Sex=='M' & !Trm%in%c('OTH','IDU','HET','BLOOD','PREG','HETfa','NEEACC','SXCH') )
	df.sc.negT		<- subset(tmp, !is.na(NegT))
	setkey(df.sc.negT, Patient)
	df.sc.negT		<- unique(df.sc.negT)
	cat(paste('\nnumber of seroconverters with !is.na(NegT), n=',nrow(df.sc.negT)))
	
	tmp				<- df.sc.negT[, list(	dt.NegT=as.numeric(difftime(AnyPos_T1, NegT, units='days'))/365,
					mp.NegT=round(as.numeric(difftime(AnyPos_T1, NegT, units='days'))/2),
					dt.CD4=as.numeric(difftime(PosCD4_T1, AnyPos_T1, units='days'))/365,
					AnyPos_a=as.numeric(difftime(AnyPos_T1, DateBorn, units='days'))/365,
					AnyPos_y=hivc.db.Date2numeric(AnyPos_T1)), by='Patient']
	
	df.sc.negT		<- merge(df.sc.negT, tmp, by='Patient')
	hist(df.sc.negT[,dt.NegT], breaks=100)
	tmp				<- subset(df.sc.negT, NegT_Acc=='Yes' )
	cat(paste('\nnumber of seroconverters with accurate NegT, n=',nrow(tmp)))
	hist(tmp[,dt.NegT], breaks=100)
	tmp				<- subset(df.sc.negT, isAcute=='No' )
	cat(paste('\nnumber of seroconverters with chronic, n=',nrow(tmp)))
	hist(tmp[,dt.NegT], breaks=100)
	tmp				<- subset(df.sc.negT, NegT_Acc=='Yes' & isAcute=='No' )
	cat(paste('\nnumber of seroconverters with chronic and accurate NegT, n=',nrow(tmp)))
	hist(tmp[,dt.NegT], breaks=100)
	
	#
	#	FOCUS on CHRONIC		isAcute=='No'
	#
	#	CALENDAR YEAR	
	df.scchr.negT	<- subset( df.sc.negT, isAcute=='No')	
	hist(df.scchr.negT[, NegT], breaks=20)
	plot(df.scchr.negT[, AnyPos_y], df.scchr.negT[, NegT], pch=18)
	plot(df.scchr.negT[, AnyPos_y], df.scchr.negT[, dt.NegT], pch=18)
	plot(df.scchr.negT[, AnyPos_y], df.scchr.negT[, AnyPos_a], pch=18)
	
	#	make sure 	time to last HIV-1 can be up to 20 yrs
	df.scchr.negT	<- subset( df.sc.negT, isAcute=='No' & AnyPos_y>2003)
	df.scchr.negT[, CD4.col:='black']
	set(df.scchr.negT, df.scchr.negT[, which(dt.CD4<1 & PosCD4_T1<AnyT_T1 & CD4_T1<250)], 'CD4.col', 'red')
	set(df.scchr.negT, df.scchr.negT[, which(dt.CD4<1 & PosCD4_T1<AnyT_T1 & CD4_T1>=250 & CD4_T1<850)], 'CD4.col', 'blue')
	set(df.scchr.negT, df.scchr.negT[, which(dt.CD4<1 & PosCD4_T1<AnyT_T1 & CD4_T1>850)], 'CD4.col', 'green')
	
	#AGE	
	subset(df.scchr.negT, AnyPos_a>65)
	subset(df.scchr.negT, AnyPos_a<20)
	plot(df.scchr.negT[, AnyPos_a], df.scchr.negT[, dt.NegT], pch=18, col=df.scchr.negT[,CD4.col])
	tmp				<- loess(dt.NegT ~ AnyPos_a, df.scchr.negT, span=0.5)
	tmp				<- predict(tmp, data.frame(AnyPos_a = seq(20, 70, 1)), se=TRUE)
	lines(seq(20, 70, 1), tmp$fit, col='red')
	lines(seq(20, 70, 1), tmp$fit+2*tmp$se.fit, col='red')
	lines(seq(20, 70, 1), tmp$fit-2*tmp$se.fit, col='red')
	
	#those with missing CD4 or PosCD4 after ART start
	df.scchr.cd4m	<- subset(df.scchr.negT,  is.na(PosCD4_T1) | dt.CD4>=1 | PosCD4_T1>=AnyT_T1)
	plot( df.scchr.cd4m[, CD4_T1], df.scchr.cd4m[, dt.NegT], pch=18)
	#numbers too small to make sense. Use regression model on all data, just using age<=45
	
	
	#CD4
	df.scchr.cd4	<- subset( df.scchr.negT, dt.CD4<1 & (is.na(AnyT_T1) | PosCD4_T1<AnyT_T1))
	subset(df.scchr.cd4, CD4_T1>800)
	subset( df.scchr.negT, dt.CD4<1 & PosCD4_T1<AnyT_T1 & NegT_Acc=='Yes' )
	plot( df.scchr.cd4[, CD4_T1], df.scchr.cd4[, dt.NegT], pch=18)
	tmp				<- loess(dt.NegT ~ CD4_T1, df.scchr.negT, span=0.5)
	tmp				<- predict(tmp, data.frame(CD4_T1 = seq(20, 1100, 10)), se=TRUE)
	lines(seq(20, 1100, 10), tmp$fit, col='red')
	lines(seq(20, 1100, 10), tmp$fit+2*tmp$se.fit, col='red')
	lines(seq(20, 1100, 10), tmp$fit-2*tmp$se.fit, col='red')
	
	plot( df.scchr.cd4[, AnyPos_a], df.scchr.cd4[, CD4_T1], pch=18)
	tmp				<- loess(CD4_T1 ~ AnyPos_a, df.scchr.negT, span=0.5)
	tmp				<- predict(tmp, data.frame(AnyPos_a = seq(20, 70, 1)), se=TRUE)
	lines(seq(20, 70, 1), tmp$fit, col='red')
	lines(seq(20, 70, 1), tmp$fit+2*tmp$se.fit, col='red')
	lines(seq(20, 70, 1), tmp$fit-2*tmp$se.fit, col='red')
	
	#	regression model
	require(MASS)	
	m1	<- glm.nb(mp.NegT ~  AnyPos_a + CD4_T1, data=df.scchr.negT, link=identity, trace= 1, maxit= 50)
	m1b	<- glm.nb(mp.NegT ~  AnyPos_a:CD4.col, data=df.scchr.negT, link=identity, trace= 1, maxit= 50)
	
	m2a	<- glm.nb(mp.NegT ~  AnyPos_a, data=subset(df.scchr.cd4, CD4_T1>850), link=identity, trace= 1, maxit= 50)
	#AIC: 153.51
	m2b	<- glm.nb(mp.NegT ~  AnyPos_a, data=subset(df.scchr.cd4, CD4_T1<=850 & CD4_T1>250), link=identity, trace= 1, maxit= 50)
	#AIC: 4568.4
	m2c	<- glm.nb(mp.NegT ~  AnyPos_a, data=subset(df.scchr.cd4, CD4_T1<250), link=identity, trace= 1, maxit= 50)
	#AIC: 2078.8
	m2d	<- glm.nb(mp.NegT ~  AnyPos_a:CD4_T1, data=subset(df.scchr.cd4, CD4_T1<250), link=identity, trace= 1, maxit= 50)
	
	#	segmented regression models
	df.scchr.cd4[, AnyPos_af:= factor(df.scchr.cd4[, AnyPos_a>=45], labels=c('<45','>=45'))]
	m3a	<- glm.nb(mp.NegT ~  AnyPos_af/AnyPos_a, data=subset(df.scchr.cd4, CD4_T1>850), link=identity, trace= 1, maxit= 50)
	#AIC: 157.03	
	m3b	<- glm.nb(mp.NegT ~  AnyPos_af/AnyPos_a, data=subset(df.scchr.cd4, CD4_T1<=850 & CD4_T1>250), link=identity, trace= 1, maxit= 50)
	#AIC: 4572.1
	m3c	<- glm.nb(mp.NegT ~  AnyPos_af/AnyPos_a, data=subset(df.scchr.cd4, CD4_T1<250), link=identity, trace= 1, maxit= 50)
	#AIC: 2081
	#	plot segmented regression models
	plot( df.scchr.cd4[, AnyPos_a], df.scchr.cd4[, mp.NegT], pch=18, col=df.scchr.cd4[,CD4.col])
	#"(Intercept)" "AnyPos_af>=45" "AnyPos_af<45:AnyPos_a" "AnyPos_af>=45:AnyPos_a"
	age1	<- seq(20, 45, 1)
	age2	<- seq(45, 70, 1)
	lines(age1, m3a$coefficients["(Intercept)"] + m3a$coefficients["AnyPos_af<45:AnyPos_a"]*age1, col='green')	
	lines(age2, m3a$coefficients["(Intercept)"] + m3a$coefficients["AnyPos_af>=45"] + m3a$coefficients["AnyPos_af>=45:AnyPos_a"]*age2, col='green')	
	lines(age1, m3b$coefficients["(Intercept)"] + m3b$coefficients["AnyPos_af<45:AnyPos_a"]*age1, col='blue')	
	lines(age2, m3b$coefficients["(Intercept)"] + m3b$coefficients["AnyPos_af>=45"] + m3b$coefficients["AnyPos_af>=45:AnyPos_a"]*age2, col='blue')
	lines(age1, m3c$coefficients["(Intercept)"] + m3c$coefficients["AnyPos_af<45:AnyPos_a"]*age1, col='red')	
	lines(age2, m3c$coefficients["(Intercept)"] + m3c$coefficients["AnyPos_af>=45"] + m3c$coefficients["AnyPos_af>=45:AnyPos_a"]*age2, col='red')
	
	#	reasonable to fit a+bx to CD4<250, CD4>250 & CD4<850 [using only age<45], const to CD4>850
	#	const model for high CD4
	m4a	<- glm.nb(mp.NegT ~  1, data=subset(df.scchr.cd4, CD4_T1>850 & dt.NegT<4), link=identity, trace= 1, maxit= 50)
	#	all CD4 model when first CD4 after ART start
	m4d	<- glm.nb(mp.NegT ~  AnyPos_a, data=subset(df.scchr.negT, AnyPos_a<=45), link=identity, trace= 1, maxit= 50)	
	#	linear model for medium and low CD4, using same intercept
	m4b	<- glm.nb(mp.NegT ~  AnyPos_a+0+offset(rep(m4d$coefficients["(Intercept)"],length(mp.NegT))), data=subset(df.scchr.cd4, CD4_T1<=850 & CD4_T1>250 & AnyPos_a<=45), link=identity, trace= 1, maxit= 50)
	m4c	<- glm.nb(mp.NegT ~  AnyPos_a+0+offset(rep(m4d$coefficients["(Intercept)"],length(mp.NegT))), data=subset(df.scchr.cd4, CD4_T1<250 & AnyPos_a<=45), link=identity, trace= 1, maxit= 50)
		
	age1	<- seq(20, 45, 1)
	age2	<- seq(45, 70, 1)
	plot( df.scchr.cd4[, AnyPos_a], df.scchr.cd4[, mp.NegT], pch=18, col=df.scchr.cd4[,CD4.col])
	lines(age1, rep(m4a$coefficients["(Intercept)"], length(age1)), col='green')	
	lines(age2, rep(m4a$coefficients["(Intercept)"], length(age2)), col='green')	
	lines(age1, m4d$coefficients["(Intercept)"] + m4b$coefficients["AnyPos_a"]*age1, col='blue')		
	lines(age1, m4d$coefficients["(Intercept)"] + m4c$coefficients["AnyPos_a"]*age1, col='red')
	lines(age1, m4d$coefficients["(Intercept)"] + m4d$coefficients["AnyPos_a"]*age1, col='black')
}
######################################################################################
project.athena.Fisheretal.get.dated.phylo.for.selection<- function(df.tpairs, clu.indir, clu.infile, clu.insignat, clu.infilexml.opt, clu.infilexml.template, method.nodectime='any')
{
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
		tmp					<- sort(setdiff( df.tpairs[,unique(cluster)], file.info[, unique(cluster)] ))
		cat(paste('\nnumber of clusters missing for analysis of pot transmitters, n=',length(tmp)))
		print(tmp)
		file.info	<- merge(file.info, unique(subset(df.tpairs, select=cluster)), by='cluster')
		df.tpairs.reduced	<- merge(df.tpairs, subset(file.info, select=cluster), by='cluster')
		cat(paste('\nnumber of remaining Patients with one potential transmitter', length(df.tpairs.reduced[, unique(Patient)])))
		cat(paste('\nnumber of remaining potential transmitter', length(df.tpairs.reduced[, unique(t.Patient)])))
	}			
	#	combine dated cluster phylogenies
	clu			<- hivc.beast2out.combine.clu.trees(clu.indir, file.info, method.nodectime=method.nodectime)
	list(clu=clu, df.tpairs=df.tpairs.reduced)
}
######################################################################################
project.athena.Fisheretal.plot.selected.transmitters<- function(df.all, df.immu, df.viro, df.treatment, df.tpairs, cluphy, cluphy.info, cluphy.subtrees, cluphy.map.nodectime, file, pdf.height=600, label.select=c("cluster","PatientA","Trm","isAcute","RegionHospital"))
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
	df.tpairs.plot	<- NULL
	if(0 && !is.null(df.tpairs))
	{
		df.tpairs.plot		<- subset( df.tpairs, select=c(cluster, FASTASampleCode, t.FASTASampleCode))
		setkey(df.tpairs.plot, cluster)		
		#tmp					<- df.tpairs.plot[, list(col=rev(brewer.pal( max(3,length(FASTASampleCode)), 'Dark2' )[seq_along(FASTASampleCode)])), by='cluster']
		tmp					<- df.tpairs.plot[, list(col=rainbow(length(FASTASampleCode))), by='cluster']
		df.tpairs.plot[, col:=tmp[,col]]
		df.tpairs.plot[, pch:=0]
		df.tpairs.plot[, t.pch:=5]
		df.tpairs.plot[, cex:=1]
		df.tpairs.plot[, t.cex:=1.4]
		tmp					<- subset(df.tpairs.plot, select=c(cluster, t.FASTASampleCode, col, t.pch, t.cex))
		setnames(tmp, c('t.FASTASampleCode', 't.pch', 't.cex'), c('FASTASampleCode', 'pch', 'cex'))
		df.tpairs.plot		<- rbind( subset(df.tpairs.plot, select=c(cluster, FASTASampleCode, col, pch, cex)), tmp )
		setkey(df.tpairs.plot, 'cluster','FASTASampleCode')
		tmp							<- setdiff( df.tpairs.plot[, unique(FASTASampleCode)], cluphy$tip.label )
		cat(paste('cannot find tip labels for df.tpairs, n=',length(tmp)))
		df.tpairs.plot				<- merge( df.tpairs.plot, data.table(FASTASampleCode=cluphy$tip.label), by='FASTASampleCode' )		
	}		
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
	cat(paste('\nplot to file=',file))
	#ph=cluphy; ph.root.ctime=cluphy.root.ctime; ph.tip.ctime=cluphy.tip.ctime; ph.prob=NA; df.node.ctime=cluphy.map.nodectime; df.rates=NULL; df.tips=df.tpairs.plot; end.ctime=end.ctime;  cex.nodelabel=0.5;  cex.tiplabel=0.5;  file=NA;  pdf.width=7
	if(!'cluster'%in%colnames(df.all))
		df.all			<- merge(df.all, subset(cluphy.info, select=c(FASTASampleCode, cluster)), by='FASTASampleCode')
	dummy				<- hivc.beast2out.plot.cluster.trees(	df.all, df.immu, df.viro, df.treatment, cluphy, cluphy.root.ctime, cluphy.tip.ctime, ph.prob=NA, df.node.ctime=cluphy.map.nodectime, df.rates=NULL, df.tips=df.tpairs.plot, 
																end.ctime=end.ctime,  cex.nodelabel=0.5,  cex.tiplabel=0.5,  file=file,  pdf.width=7, pdf.height=pdf.height, pdf.xlim=pdf.xlim, label.select=label.select)	
}
######################################################################################
project.athena.Fisheretal.X.b4care<- function(df.tpairs, df.all, predict.t2inf, t2inf.args, t.period=0.25, ts.min=1980, method.minLowerUWithNegT=TRUE)
{
	b4care	<- merge(data.table(Patient=df.tpairs[, unique(t.Patient)]), unique(subset(df.all, select=c(Patient, DateBorn, AnyPos_T1, isAcute, NegT, PosCD4_T1, CD4_T1, AnyT_T1))), by='Patient')
	if(method.minLowerUWithNegT)
	{
		cat('\nUsing option method.minLowerUWithNegT=TRUE')
		tmp		<- b4care[, list(AnyPos_a=AnyPos_T1-DateBorn, ts=ifelse(is.na(NegT), AnyPos_T1-10, NegT), te=AnyPos_T1), by='Patient']	
	}	
	if(!method.minLowerUWithNegT)
	{
		cat('\nUsing option method.minLowerUWithNegT=FALSE')
		tmp		<- b4care[, list(AnyPos_a=AnyPos_T1-DateBorn, ts=AnyPos_T1-10, te=AnyPos_T1), by='Patient']	
	}	
	b4care	<- merge(b4care, tmp, by='Patient')
	#check if ts before 1980 and if so clip
	set(b4care, b4care[,which(ts<ts.min)], 'ts', ts.min )	
	set(b4care, NULL, 'ts', b4care[, floor(ts) + round( (ts%%1)*100 %/% (t.period*100) ) * t.period ] )
	set(b4care, NULL, 'te', b4care[, floor(te) + round( (te%%1)*100 %/% (t.period*100) ) * t.period - t.period] )	#last time undiagnosed is just before first time period diagnosed
	#apply predict.t2inf	
	t2inf	<- predict.t2inf(seq(0,10,t.period)*365, b4care, t2inf.args)
	t2inf	<- merge( subset( t2inf, !is.na(score) ), subset(b4care, select=c(Patient,ts,te)), by='Patient' )
	set(t2inf, NULL, 'q', t2inf[,te-q/365])
	t2inf	<- subset(t2inf, ts<=q, c(Patient, q, score))
	b4care	<- merge(t2inf, subset(b4care, select=c(Patient, AnyPos_T1, AnyPos_a, isAcute)), by='Patient')	
	set(b4care, NULL, 'score', b4care[, round(score, d=3)])
	setnames(b4care, c('Patient','q','score'), c('t.Patient','t','U.score'))	
	cat(paste('\nReturn X for #t.Patients=',b4care[, length(unique(t.Patient))]))
	b4care
}
######################################################################################
project.athena.Fisheretal.YX.age<- function(YX, clumsm.info, t.period=0.25)
{
	tmp		<- unique(subset(clumsm.info, select=c(Patient, DateBorn)))
	setnames(tmp, 'Patient', 't.Patient')
	age		<- merge(subset(YX, select=c(t.Patient, Patient, t)), tmp, by='t.Patient')
	setkey(age, t, Patient, t.Patient)
	age		<- unique(age)
	set(age, NULL, 'DateBorn', age[, t+t.period/2-DateBorn])
	setnames(age, 'DateBorn', 't.Age')
	age		<- merge(age, unique(subset(clumsm.info, select=c(Patient, DateBorn))), by='Patient')
	set(age, NULL, 'DateBorn', age[, t+t.period/2-DateBorn])
	setnames(age, 'DateBorn', 'Age')
	YX	<- merge(YX, age, by=c('t','t.Patient','Patient'))
	YX
}
######################################################################################
project.athena.Fisheretal.X.ART.pulsed<- function(df.tpairs, clumsm.info, df.treatment, t.pulse.ART=0.75, t.pulse.sc=NA, t.pulse.ART.I=1, t.pulse.ART.I.d=1)
{
	pulse		<- merge( data.table(Patient= df.tpairs[, unique(t.Patient)]), unique(subset(clumsm.info, select=c(Patient, NegT, AnyPos_T1, AnyT_T1))), by='Patient' )
	pulse		<- subset(pulse, !is.na(AnyT_T1) & AnyPos_T1<AnyT_T1 & AnyPos_T1+t.pulse.ART>=AnyT_T1)	
	if(!is.na(t.pulse.sc))
		pulse	<- subset(pulse, !is.na(NegT) & NegT+t.pulse.sc>=AnyPos_T1 )
	
	pulse		<- merge( subset(pulse, select=c(Patient, NegT, AnyPos_T1)), df.treatment, by='Patient')
	set(pulse, NULL, 'StartTime', hivc.db.Date2numeric(pulse[,StartTime]))
	set(pulse, NULL, 'StopTime', hivc.db.Date2numeric(pulse[,StopTime]))
	set(pulse, NULL, 'AnyT_T1', hivc.db.Date2numeric(pulse[,AnyT_T1]))
	pulse		<- subset(pulse, TrI=='Yes' & AnyT_T1+t.pulse.ART.I>=StartTime & StopTime-StartTime>t.pulse.ART.I.d, c(Patient, NegT, AnyPos_T1, AnyT_T1, StartTime, StopTime, TrI ))
	setkey(pulse, Patient)
	pulse		<- unique(pulse)
	cat(paste('\n number of potential transmitters with pulse ART shortly after seroconversion/diagnosis, n=',nrow(pulse)))
	pulse		<- subset(pulse, select=Patient)
	setnames(pulse, 'Patient', 't.Patient')
	pulse[, ART.pulse:='Yes']
	pulse		<- merge( unique(subset(df.tpairs, select=t.Patient)), pulse, all.x=1, by='t.Patient' )
	set(pulse, pulse[, which(is.na(ART.pulse))], 'ART.pulse', 'No')
	pulse	
}
######################################################################################
project.athena.Fisheretal.X.incare<- function(df.tpairs, df.all, df.viro, df.immu, df.treatment, t.period=0.25, t.endctime= 2013., lRNA.supp=3, ART.start.period.tmax=0.5)
{
	#	prepare incare timeline for potential transmitters
	incare		<- merge( data.table(Patient=df.tpairs[, unique(t.Patient)]), unique( subset(df.all, select=c(Patient, AnyPos_T1, DateDied)) ), by='Patient' )
	cat(paste('\nPot transmitters, n=',incare[, length(unique(Patient))]))
	set(incare, NULL, 'AnyPos_T1', incare[, floor(AnyPos_T1) + round( (AnyPos_T1%%1)*100 %/% (t.period*100) ) * t.period] )
	set(incare, incare[,which(is.na(DateDied))], 'DateDied', t.endctime)
	set(incare, NULL, 'DateDied', incare[, floor(DateDied) + round( (DateDied%%1)*100 %/% (t.period*100) ) * t.period] )
	#	
	tmp			<- subset( incare, AnyPos_T1>=DateDied )
	if(nrow(tmp))	cat(paste('\nremoving potential transmitters from incare that are dead at diagnosis, n=',nrow(tmp)))
	incare		<- subset( incare, AnyPos_T1<DateDied )
	incare.t	<- incare[, list(t= seq(AnyPos_T1, DateDied-t.period, by=t.period)),by='Patient']	
	incare.t[, stage:='Diag']
	setnames(incare.t, 'Patient','t.Patient')
	cat(paste('\nincare.t entries, n=',nrow(incare.t)))		
	#	compute X: viral load of potential transmitter  
	X.viro		<- project.athena.Fisheretal.X.viro(df.tpairs, df.viro, t.period=t.period, lRNA.supp=lRNA.supp)	
	cat(paste('\nVL info available for pot transmitters, n=',X.viro[, length(unique(t.Patient))]))
	#	compute X: CD4 of potential transmitter
	X.cd4		<- project.athena.Fisheretal.X.cd4(df.tpairs, df.immu, t.period=t.period)
	cat(paste('\nCD4 info available for pot transmitters, n=',X.cd4[, length(unique(t.Patient))]))
	#	add viro and cd4 time periods to incare.t		
	incare.t	<- merge(incare.t, X.viro, by=c('t.Patient','t'), all.x=1)
	incare.t	<- merge(incare.t, X.cd4, by=c('t.Patient','t'), all.x=1)
	#	prepare treatment variables for potential transmitters
	treat		<- subset(df.treatment, select=c(Patient, AnyT_T1, StartTime, StopTime, TrI, TrCh.failure, TrVL.failure, TrCh.toxicity, TrCh.adherence, TrCh.patrel, NoDrug, NoNRT, NoNNRT, NoPI, NoBoost, TDF, EFV, FTC))
	treat[, TDF.EFV.FTC:= as.numeric(TDF & EFV & FTC)]	#Atripla ingredients	
	#set(treat, treat[, which(is.na(TrVL.failure) & TrCh.failure=='No')], 'TrVL.failure', 'No')
	#set(treat, treat[, which(is.na(TrVL.failure) & TrCh.failure=='Yes')], 'TrVL.failure', 'Yes')	
	treat		<- merge( data.table(Patient=df.tpairs[, unique(t.Patient)]), treat, by='Patient' )
	set(treat, NULL, 'AnyT_T1', hivc.db.Date2numeric(treat[,AnyT_T1]))
	set(treat, NULL, 'StopTime', hivc.db.Date2numeric(treat[,StopTime]))
	set(treat, NULL, 'StartTime', hivc.db.Date2numeric(treat[,StartTime]))
	set(treat, NULL, 'StartTime', treat[, floor(StartTime) + round( (StartTime%%1)*100 %/% (t.period*100) ) * t.period] )
	set(treat, NULL, 'StopTime', treat[, floor(StopTime) + round( (StopTime%%1)*100 %/% (t.period*100) ) * t.period] )
	treat		<- subset(treat, StartTime!=StopTime)	#delete too small Start to Stop time intervals	-- we can have interrupted at start
	cat(paste('\nTreatment info available for pot transmitters, n=',treat[, length(unique(Patient))]))
	setnames(treat, 'Patient', 't.Patient')
	#	if there are treatment interruptions at start, then this is because treatment periods have been rounded away.
	#	keep as is.
	#	remove first ART periods that begin with an interruption
	#treat		<- merge( treat, subset( treat[, list(select=!all(TrI=='Yes')), by='Patient'], select, Patient), by='Patient' )
	#cat(paste('\nTreatment info for pot transmitters for which not all periods are ART interrupted, n=',treat[, length(unique(Patient))]))
	#treat		<- treat[, {
	#							tmp<- seq.int(which(TrI!='Yes')[1], length(TrI))
	#							lapply(.SD,'[',tmp)
	#						},by='Patient']	
	#	get treatment timeline
	treat.t		<- merge(subset(incare.t, select=c(t.Patient, t)), treat, by='t.Patient', allow.cartesian=TRUE )	
	treat.t		<- subset( treat.t, StartTime<=t+t.period/2 & t+t.period/2<StopTime )
	cat(paste('\ntreat.t entries, n=',nrow(treat.t)))
	setnames(treat.t, c('TrI','TrVL.failure','TrCh.toxicity','TrCh.adherence','TrCh.patrel','NoDrug','NoNRT','NoNNRT','NoPI','NoBoost','TDF.EFV.FTC'), c('ART.I','ART.F','ART.T','ART.A','ART.P','ART.nDrug','ART.nNRT','ART.nNNRT','ART.nPI','ART.nBoost','ART.TDF.EFV.FTC'))
	treat.t		<- subset(treat.t, select=c(t.Patient,t,AnyT_T1,StartTime,StopTime,ART.I,ART.F,ART.T,ART.A,ART.P,ART.nDrug,ART.nNRT,ART.nNNRT,ART.nPI,ART.nBoost,ART.TDF.EFV.FTC))
	#	merge incare and treatment timelines
	incare.t	<- merge(incare.t, treat.t, all.x=1, by=c('t.Patient','t'))	
	#	set ever on ART per period t		(avoid AnyT_T1 as it may give periods that would start with treatment interruption)
	set(incare.t, incare.t[, which(!is.na(StartTime))], 'stage', 'ART.started')
	#	set virological failure by smoothed VL, keep NA ART.F
	set(incare.t, incare.t[, which(stage!='ART.started')], 'ART.F', NA_character_)
	set(incare.t, incare.t[, which(stage=='ART.started' & !is.na(lRNA) & lRNA>lRNA.supp)], 'ART.F', 'Yes')
	set(incare.t, incare.t[, which(stage=='ART.started' & !is.na(lRNA) & lRNA<=lRNA.supp)], 'ART.F', 'No')
	#	set ART start period
	incare.t[, ART.startperiod:= NA_character_]
	set(incare.t, incare.t[, which(stage=='ART.started' & t-AnyT_T1<=ART.start.period.tmax)], 'ART.startperiod', 'Yes')
	set(incare.t, incare.t[, which(stage=='ART.started' & t-AnyT_T1>ART.start.period.tmax)], 'ART.startperiod', 'No')
	cat(paste("\nNumber entries with ART.I=='Yes' & stage=='Diag', n=",nrow(subset(incare.t, ART.I=='Yes' & stage=='Diag'))))
	cat(paste('\nReturn X for #t.Patients=',incare.t[, length(unique(t.Patient))]))	
	incare.t
}
######################################################################################
project.athena.Fisheretal.YX.model1.whichregression<- function(YX.m1.info, outdir, outfile, outsignat)
{
	require(MASS)
	require(fitdistrplus)
	#z			<- subset(YX.m1.info, score.Y<0.999 & stage=='U')					
	#score.Y		<- z[,score.Y]
	#logit.Y		<- z[,logit.Y]
	#mle			<- fitdist(score.Y, 'beta', start=list(shape1 = 2, shape2 = 2))
	#denscomp(mle, addlegend=FALSE, xlab='beta all')
	#qqcomp(mle, addlegend=FALSE, main='beta all', cex=0.25)				
	#mle2		<- fitdist(score.Y, 'logis')
	#denscomp(mle2, addlegend=FALSE, xlab='logit all')
	#qqcomp(mle2, addlegend=FALSE, main='logit all', cex=0.25)		
	#tmp			<- gofstat(list(mle, mle2), fitnames=c('beta','logit'))		
	plot.file	<- paste(outdir,'/',outfile, '_',gsub('/',':',outsignat),'_regressionmodel1_3.pdf',sep='')
	pdf(file=plot.file, w=15, h=5)
	par(mfrow=c(1, 4))
	#	plot all
	z			<- subset(YX.m1.info, score.Y<0.999 & stage=='U')					
	score.Y		<- z[,score.Y]
	logit.Y		<- z[,logit.Y]
	mle			<- fitdist(score.Y, 'beta', start=list(shape1 = 2, shape2 = 2))
	qqcomp(mle, addlegend=FALSE, main='beta all')
	denscomp(mle, addlegend=FALSE, xlab='beta all')
	#	plot by strata				
	mle2		<- fitdist(score.Y, 'logis')
	qqcomp(mle2, addlegend=FALSE, main='logit all')
	denscomp(mle2, addlegend=FALSE, xlab='logit all')				
	tmp			<- gofstat(list(mle, mle2), fitnames=c('beta','logit'))
	
	YX.m1.gof	<- subset(YX.m1.info, score.Y<0.999)[,	
			{
				mle			<- fitdist(score.Y, 'beta', start=list(shape1 = 2, shape2 = 2))
				qqcomp(mle, addlegend=FALSE, main=paste('beta',stage))
				denscomp(mle, addlegend=FALSE, xlab=paste('beta',stage))
				
				mle2		<- fitdist(score.Y, 'logis')
				qqcomp(mle2, addlegend=FALSE, main=paste('logit',stage))
				denscomp(mle2, addlegend=FALSE, xlab=paste('logit',stage))
				
				tmp			<- gofstat(list(mle, mle2), fitnames=c('beta','logit'))
				ans			<- c( as.list( tmp$kstest ), as.list( tmp$aic ) )
				names(ans)	<- c('beta.ks','logit.ks','beta.aic','logit.aic')
				ans
			},by='stage']							
	dev.off()			
	#	get goodness of fit summaries
	YX.m1.gof	<- rbind( data.table(stage='all', beta.ks=tmp$kstest['beta'] , logit.ks=tmp$kstest['logit'], beta.aic=tmp$aic['beta'] , logit.aic=tmp$aic['logit']), YX.m1.gof)
	YX.m1.gof
}
######################################################################################
project.athena.Fisheretal.YX.model5.stratify<- function(YX)
{
	YX.m5	<- copy(YX)
	YX.m5[, U.score:=NULL]
	if('score.Y'%in%colnames(YX.m5) && YX.m5[, !any(score.Y>1.1)])
	{
		tmp	<- YX.m5[, score.Y>0]
		set(YX.m5, which(tmp), 'score.Y', YX.m5[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )
	}			
	#	remove Acute and ART
	set(YX.m5, YX.m5[, which(t.isAcute=='Yes' | t.isAcute=='Maybe')], 'stage', NA_character_)
	set(YX.m5, YX.m5[, which(grepl('ART',stage))], 'stage', NA_character_)
	#	remove early calendar times?
	#set(YX.m5, YX.m5[, which(t<=2006.25)], 'stage', NA_character_)
	YX.m5	<- subset(YX.m5, !is.na(stage))
	#	set Age of Transmitter tA				
	age.cut		<- c(-1, 20, 25, 30, 35, 45, 100)
	YX.m5[, tA:=NA_character_]
	for(i in seq_along(age.cut)[-1])
	{
		tmp	<- paste('t<=',age.cut[i],sep='')
		set(YX.m5, YX.m5[, which(stage%in%c('U','Diag') & t.Age<=age.cut[i] & t.Age>age.cut[i-1])], 'tA', tmp)
	}		
	set(YX.m5, NULL, 'tA', YX.m5[, factor(as.character(tA))])
	#	set Age of Transmitter and Age of recipient
	setkey(YX.m5, t.Age, Age)
	YX.m5[, tiA:=NA_character_]
	for(i in seq_along(age.cut)[-1])
		for(j in seq.int(2, length(age.cut)))
		{
			tmp	<- paste('t<=',age.cut[i],'-i<=',age.cut[j],sep='')
			set(YX.m5, YX.m5[, which(stage%in%c('U','Diag') & t.Age<=age.cut[i] & t.Age>age.cut[i-1] & Age<=age.cut[j] & Age>age.cut[j-1])], 'tiA', tmp)	
		}
	set(YX.m5, NULL, 'tiA', YX.m5[, factor(as.character(tiA))])
	#	add tperiod
	if('t.period'%in%colnames(YX.m5))
	{
		cat(paste('\nadding tA.tperiod and tiA.tperiod\n'))	
		YX.m5[, tA.tperiod:= paste(tA, t.period,sep='.')]
		set(YX.m5, NULL, 'tA.tperiod', YX.m5[, factor(as.character(tA.tperiod))])
		YX.m5[, tiA.tperiod:= paste(tiA, t.period,sep='.')]
		set(YX.m5, NULL, 'tiA.tperiod', YX.m5[, factor(as.character(tiA.tperiod))])
	}	
	gc()
	#	set Age of Transmitter tA				
	age.cut		<- c(-1, 20, 28, 38, 48, 58, 100)
	YX.m5[, tAb:=NA_character_]
	for(i in seq_along(age.cut)[-1])
	{
		tmp	<- paste('t<=',age.cut[i],sep='')
		set(YX.m5, YX.m5[, which(stage%in%c('U','Diag') & t.Age<=age.cut[i] & t.Age>age.cut[i-1])], 'tAb', tmp)
	}		
	set(YX.m5, NULL, 'tAb', YX.m5[, factor(as.character(tAb))])
	#	set Age of Transmitter and Age of recipient
	setkey(YX.m5, t.Age, Age)
	YX.m5[, tiAb:=NA_character_]
	for(i in seq_along(age.cut)[-1])
		for(j in seq.int(2, length(age.cut)))
		{
			tmp	<- paste('t<=',age.cut[i],'-i<=',age.cut[j],sep='')
			set(YX.m5, YX.m5[, which(stage%in%c('U','Diag') & t.Age<=age.cut[i] & t.Age>age.cut[i-1] & Age<=age.cut[j] & Age>age.cut[j-1])], 'tiAb', tmp)	
		}
	set(YX.m5, NULL, 'tiAb', YX.m5[, factor(as.character(tiAb))])
	#	add tperiod
	if('t.period'%in%colnames(YX.m5))
	{
		cat(paste('\nadding tAb.tperiod and tiAb.tperiod\n'))	
		YX.m5[, tAb.tperiod:= paste(tAb, t.period,sep='.')]
		set(YX.m5, NULL, 'tAb.tperiod', YX.m5[, factor(as.character(tAb.tperiod))])
		YX.m5[, tiAb.tperiod:= paste(tiAb, t.period,sep='.')]
		set(YX.m5, NULL, 'tiAb.tperiod', YX.m5[, factor(as.character(tiAb.tperiod))])
	}	
	gc()
	#	set Age of Transmitter tA				
	age.cut		<- c(-1, 20, 25, 30, 35, 40, 45, 100)
	YX.m5[, tAc:=NA_character_]
	for(i in seq_along(age.cut)[-1])
	{
		tmp	<- paste('t<=',age.cut[i],sep='')
		set(YX.m5, YX.m5[, which(stage%in%c('U','Diag') & t.Age<=age.cut[i] & t.Age>age.cut[i-1])], 'tAc', tmp)
	}		
	set(YX.m5, NULL, 'tAc', YX.m5[, factor(as.character(tAc))])
	#	set Age of Transmitter and Age of recipient
	setkey(YX.m5, t.Age, Age)
	YX.m5[, tiAc:=NA_character_]
	for(i in seq_along(age.cut)[-1])
		for(j in seq.int(2, length(age.cut)))
		{
			tmp	<- paste('t<=',age.cut[i],'-i<=',age.cut[j],sep='')
			set(YX.m5, YX.m5[, which(stage%in%c('U','Diag') & t.Age<=age.cut[i] & t.Age>age.cut[i-1] & Age<=age.cut[j] & Age>age.cut[j-1])], 'tiAc', tmp)	
		}
	set(YX.m5, NULL, 'tiAc', YX.m5[, factor(as.character(tiAc))])
	#	add tperiod
	if('t.period'%in%colnames(YX.m5))
	{
		cat(paste('\nadding tAc.tperiod and tiAc.tperiod\n'))	
		YX.m5[, tAc.tperiod:= paste(tAc, t.period,sep='.')]
		set(YX.m5, NULL, 'tAc.tperiod', YX.m5[, factor(as.character(tAc.tperiod))])
		YX.m5[, tiAc.tperiod:= paste(tiAc, t.period,sep='.')]
		set(YX.m5, NULL, 'tiAc.tperiod', YX.m5[, factor(as.character(tiAc.tperiod))])
	}		
	gc()
	cat(paste('\nsubset\n'))
	if('score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, score.Y, stage, CDCC, lRNA, contact, fw.up.med, t.period, w, w.i, w.in, w.t, w.tn, AnyPos_T1, t.AnyPos_T1, t.AnyT_T1, tA, tiA, tA.tperiod, tiA.tperiod, tAb, tiAb, tAb.tperiod, tiAb.tperiod, tAc, tiAc, tAc.tperiod, tiAc.tperiod, t.Age, t.RegionHospital  ))	
	if(!'score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, stage, CDCC, lRNA, contact, fw.up.med, t.period, AnyPos_T1, t.AnyPos_T1, t.AnyT_T1, tA, tiA, tA.tperiod, tiA.tperiod, tAb, tiAb, tAb.tperiod, tiAb.tperiod, tAc, tiAc, tAc.tperiod, tiAc.tperiod  ))	
	gc()
	YX.m5
}
######################################################################################
project.athena.Fisheretal.YX.model5<- function(YX, clumsm.info, vl.suppressed=log10(1e3), acute.select=c('Yes','Maybe'), cd4.cut= c(-1, 350, 550, 5000) )
{
	require(betareg)
	#	Age groups
	#	
	#	PQ:	Are the younger ages associated with high transmission prob? across Diag+U ? 
	#		--> by time period at age at diagnosis?
	#		--> across all parts of the treatment cascade?
	#	PQ: Are late presenters associated with high transmission prob? (CD4<=220, CDCC at diag, Age at diag > 45)
	#		--> by time period at age at diagnosis?
	#		--> across Diag+U ? across all parts of the treatment cascade?
	#	PQ:	Age differences ?
	YX.m5	<- copy(YX)
	YX.m5[, U.score:=NULL]
	#	if transmitter acute, set undiagnosed to undiagnosed and acute
	set(YX.m5, YX.m5[, which(stage=='U' & t.isAcute=='Yes')], 'stage', 'UAy')
	set(YX.m5, YX.m5[, which(stage=='U' & t.isAcute=='Maybe')], 'stage', 'UAm')
	set(YX.m5, YX.m5[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.5)], 'stage', 'DAy' )
	set(YX.m5, YX.m5[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.5)], 'stage', 'DAm' )					
	#
	#	stratify ART into suppressed during infection window and not suppressed during infection window
	#	minimize NAs by using extrapolation for 6 months after last known VL
	#
	VL.cur	<- 3
	tmp		<- unique(subset( clumsm.info, select=c(Patient, PoslRNA_T1) ))
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	YX.m5	<- merge(YX.m5, tmp, by= 't.Patient')					
	YX.m5	<- merge(YX.m5, YX.m5[, {	
						tmp<- which(!is.na(lRNA))
						list(lRNA_TL= ifelse(length(tmp)>0, lRNA[tail(tmp,1)], NA_real_), 	PoslRNA_TL=ifelse(length(tmp)>0, t[tail(tmp,1)], NA_real_))
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
	YX.m5	<- merge(YX.m5, YX.m5[, {
						tmp<- which(!is.na(lRNA))
						if(!length(tmp) & any(t>t.PoslRNA_T1[1]) & all(contact=='Yes') & !is.na(lRNA_TL[1]) & lRNA_TL[1]<=VL.cur)
							lRNA.c	<- 'SuA.Y'						
						else if(!length(tmp) & any(t>t.PoslRNA_T1[1]) & all(contact=='Yes') & !is.na(lRNA_TL[1]) & lRNA_TL[1]>VL.cur)
							lRNA.c	<- 'SuA.N'						
						else if(!length(tmp))
							lRNA.c	<- 'SuA.NA'
						else if(max(lRNA[tmp])>VL.cur)
							lRNA.c	<- 'SuA.N'
						else
							lRNA.c	<- 'SuA.Y'
						list( lRNA.c= lRNA.c )
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
	set(YX.m5, YX.m5[,which(stage=='ART.started' & lRNA.c=='SuA.Y')],'stage', 'ART.suA.Y' )
	set(YX.m5, YX.m5[,which(stage=='ART.started' & lRNA.c=='SuA.N')],'stage', 'ART.suA.N' )
	set(YX.m5, YX.m5[,which(stage=='ART.started')],'stage', 'ART.vlNA' )
	set(YX.m5, NULL, 'stage', YX.m5[, factor(as.character(stage))])			
	YX.m5[, lRNA.c:= NULL]
	#	age at diagnosis
	tmp		<- subset( clumsm.info, select=c(Patient, DateBorn, AnyPos_T1) )
	setkey(tmp, Patient)
	tmp		<- unique(tmp)
	tmp[, AnyPos_A:= AnyPos_T1-DateBorn]
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	tmp		<- merge( unique( subset(YX.m5, select=t.Patient) ), tmp, by='t.Patient' )		
	YX.m5	<- merge(YX.m5, subset(tmp, select=c(t.Patient, t.AnyPos_A)), by='t.Patient')
	#
	YX.m5[, stage.orig:= stage]
	YX.m5[, table(stage)]
	#	young age groups by age at diagnosis
	if(1)
	{		
		set(YX.m5, NULL, 'stage', YX.m5[, stage.orig])
		
		age.cut		<- c(-1, 20, 22, 45, 100)
		age.label	<- c('T1<=20','T1<=22','T1<=45','T1>45')
		for(i in seq_along(age.cut)[-1])
			set(YX.m5, YX.m5[, which(stage%in%c('U','Diag') & t.AnyPos_A<=age.cut[i] & t.AnyPos_A>age.cut[i-1])], 'stage', age.label[i-1])
		
		set(YX.m5, NULL, 'stage', YX.m5[, factor(as.character(stage))])
		YX.m5[, table(stage)]
		#	ART.suA.N ART.suA.Y  ART.vlNA       DAm       DAy    T1<=20    T1<=22    T1<=45     T1>45       UAm       UAy 
     	#	3172      4526        70       		436       523       293       373     10741      3512      1275      1291 
		YX.m5.fit1 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m5)
		YX.m5.fit1.or	<- data.table( 	T1l20= my.or.from.logit(YX.m5.fit1, 'stageT1<=20', 'stageT1<=45', subset(YX.m5, stage=='T1<=20')[, sum(w)], subset(YX.m5, stage=='T1<=45')[, sum(w)], 1.962),
										T1l22= my.or.from.logit(YX.m5.fit1, 'stageT1<=22', 'stageT1<=45', subset(YX.m5, stage=='T1<=22')[, sum(w)], subset(YX.m5, stage=='T1<=45')[, sum(w)], 1.962),
										T1g45= my.or.from.logit(YX.m5.fit1, 'stageT1>45', 'stageT1<=45', subset(YX.m5, stage=='T1>45')[, sum(w)], subset(YX.m5, stage=='T1<=45')[, sum(w)], 1.962)
										)		
		#	T1l20     T1l22     T1g45
		#1: 1.728989 0.5585191 1.1386048
		#2: 1.430177 0.4574844 0.9498748
		#3: 2.090232 0.6818671 1.3648333
	}
	#	young age groups by age at time period
	if(1)
	{		
		set(YX.m5, NULL, 'stage', YX.m5[, stage.orig])		
		age.cut		<- c(-1, 20, 22, 45, 100)
		age.label	<- c('tT<=20','tT<=22','tT<=45','tT>45')
		for(i in seq_along(age.cut)[-1])
			set(YX.m5, YX.m5[, which(stage%in%c('U','Diag') & t.Age<=age.cut[i] & t.Age>age.cut[i-1])], 'stage', age.label[i-1])
		
		set(YX.m5, NULL, 'stage', YX.m5[, factor(as.character(stage))])
		YX.m5[, table(stage)]
		#	ART.suA.N ART.suA.Y  ART.vlNA       DAm       DAy    tT<=20    tT<=22    tT<=45     tT>45       UAm       UAy 
     	#	3172      4526        70       		436       523       295       388     10639      3597      1275      1291  
		YX.m5.fit2 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m5)
		YX.m5.fit2.or	<- data.table( 	tTl20= my.or.from.logit(YX.m5.fit2, 'stagetT<=20', 'stagetT<=45', subset(YX.m5, stage=='tT<=20')[, sum(w)], subset(YX.m5, stage=='tT<=45')[, sum(w)], 1.962),
										tTl22= my.or.from.logit(YX.m5.fit2, 'stagetT<=22', 'stagetT<=45', subset(YX.m5, stage=='tT<=22')[, sum(w)], subset(YX.m5, stage=='tT<=45')[, sum(w)], 1.962),
										tTg45= my.or.from.logit(YX.m5.fit2, 'stagetT>45', 'stagetT<=45', subset(YX.m5, stage=='tT>45')[, sum(w)], subset(YX.m5, stage=='tT<=45')[, sum(w)], 1.962)	)		
		#	      tTl20    tTl22     tTg45
		#1: 	1.645678 2.163212 0.9566155
		#2: 	1.357772 1.785391 0.7974689
		#3: 	1.994633 2.620986 1.1475223
		#
		set(YX.m5, NULL, 'stage', YX.m5[, stage.orig])		
		age.cut		<- c(-1, 20, 25, 45, 100)
		age.label	<- c('tT<=20','tT<=25','tT<=45','tT>45')
		for(i in seq_along(age.cut)[-1])
			set(YX.m5, YX.m5[, which(stage%in%c('U','Diag') & t.Age<=age.cut[i] & t.Age>age.cut[i-1])], 'stage', age.label[i-1])
		
		set(YX.m5, NULL, 'stage', YX.m5[, factor(as.character(stage))])
		YX.m5[, table(stage)]
		#	ART.suA.N ART.suA.Y  ART.vlNA       DAm       DAy    tT<=20    tT<=25    tT<=45     tT>45       UAm       UAy 
		#	3172      4526        70       		436       523       295       388     10639      3597      1275      1291  
		YX.m5.fit2 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m5)
		YX.m5.fit2.or	<- data.table( 	tTl20= my.or.from.logit(YX.m5.fit2, 'stagetT<=20', 'stagetT<=45', subset(YX.m5, stage=='tT<=20')[, sum(w)], subset(YX.m5, stage=='tT<=45')[, sum(w)], 1.962),
				tTl25= my.or.from.logit(YX.m5.fit2, 'stagetT<=25', 'stagetT<=45', subset(YX.m5, stage=='tT<=25')[, sum(w)], subset(YX.m5, stage=='tT<=45')[, sum(w)], 1.962),
				tTg45= my.or.from.logit(YX.m5.fit2, 'stagetT>45', 'stagetT<=45', subset(YX.m5, stage=='tT>45')[, sum(w)], subset(YX.m5, stage=='tT<=45')[, sum(w)], 1.962)	)		
		
		#		tTl20    tTl25     tTg45
		#1: 1.660208 1.552963 0.9653929
		#2: 1.363555 1.273374 0.8021024
		#3: 2.021401 1.893941 1.1619259
		#	young age groups are generally associated with higher TP. T1>45 probably too once young age adjusted for
		#
		set(YX.m5, NULL, 'stage', YX.m5[, stage.orig])		
		age.cut		<- c(-1, 20, 25, 45, 100)
		age.label	<- c('tT<=20','tT<=25','tT<=45','tT>45')
		for(i in seq_along(age.cut)[-1])
			set(YX.m5, YX.m5[, which(stage%in%c('U','Diag','ART.suA.N','ART.vlNA') & t.Age<=age.cut[i] & t.Age>age.cut[i-1])], 'stage', age.label[i-1])
		
		set(YX.m5, NULL, 'stage', YX.m5[, factor(as.character(stage))])
		YX.m5[, table(stage)]
		#	ART.suA.Y       DAm       DAy    tT<=20    tT<=25    tT<=45     tT>45       UAm       UAy 
     	#		4526       436       523       308      1310     12023      4520      1275      1291 
		YX.m5.fit2 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m5)
		YX.m5.fit2.or	<- data.table( 	tTl20= my.or.from.logit(YX.m5.fit2, 'stagetT<=20', 'stagetT<=45', subset(YX.m5, stage=='tT<=20')[, sum(w)], subset(YX.m5, stage=='tT<=45')[, sum(w)], 1.962),
										tTl25= my.or.from.logit(YX.m5.fit2, 'stagetT<=25', 'stagetT<=45', subset(YX.m5, stage=='tT<=25')[, sum(w)], subset(YX.m5, stage=='tT<=45')[, sum(w)], 1.962),
										tTg45= my.or.from.logit(YX.m5.fit2, 'stagetT>45', 'stagetT<=45', subset(YX.m5, stage=='tT>45')[, sum(w)], subset(YX.m5, stage=='tT<=45')[, sum(w)], 1.962)	)		
		#      tTl20    tTl25     tTg45
		#1: 1.742191 1.456322 0.9776405
		#2: 1.453716 1.212899 0.8255276
		#3: 2.087913 1.748599 1.1577819
		#	young age groups remain significant. T1>45 probably too once young age adjusted for
	}
	#	try age pairs -- full matrix
	if(1)
	{
		set(YX.m5, NULL, 'stage', YX.m5[, stage.orig])		
		age.cut		<- c(-1, 20, 25, 35, 45, 100)		
		setkey(YX.m5, t.Age, Age)
		for(i in seq_along(age.cut)[-1])
			for(j in seq.int(2, length(age.cut)))
			{
				tmp	<- paste('t',age.cut[i],'-i',age.cut[j],sep='')
				set(YX.m5, YX.m5[, which(stage%in%c('U','Diag') & t.Age<=age.cut[i] & t.Age>age.cut[i-1] & Age<=age.cut[j] & Age>age.cut[j-1])], 'stage', tmp)	
			}
		set(YX.m5, NULL, 'stage', YX.m5[, factor(as.character(stage))])
		YX.m5[, table(stage)]
		
		YX.m5.fit3 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m5)
		summary(YX.m5.fit3)
		#	is the assortativity matrix symmetric?
		YX.m5.fit3.or	<- data.table( 	a100.20= my.or.from.logit(YX.m5.fit3, 'staget100-i20', 'staget20-i100', subset(YX.m5, stage=='t100-i20')[, sum(w)], subset(YX.m5, stage=='t20-i100')[, sum(w)], 1.962),
										a100.25= my.or.from.logit(YX.m5.fit3, 'staget100-i25', 'staget25-i100', subset(YX.m5, stage=='t100-i25')[, sum(w)], subset(YX.m5, stage=='t25-i100')[, sum(w)], 1.962),
										a100.35= my.or.from.logit(YX.m5.fit3, 'staget100-i35', 'staget35-i100', subset(YX.m5, stage=='t100-i35')[, sum(w)], subset(YX.m5, stage=='t35-i100')[, sum(w)], 1.962),
										a100.45= my.or.from.logit(YX.m5.fit3, 'staget100-i45', 'staget45-i100', subset(YX.m5, stage=='t100-i45')[, sum(w)], subset(YX.m5, stage=='t45-i100')[, sum(w)], 1.962),				
										a45.20= my.or.from.logit(YX.m5.fit3, 'staget45-i20', 'staget20-i45', subset(YX.m5, stage=='t45-i20')[, sum(w)], subset(YX.m5, stage=='t20-i45')[, sum(w)], 1.962),
										a45.25= my.or.from.logit(YX.m5.fit3, 'staget45-i25', 'staget25-i45', subset(YX.m5, stage=='t45-i25')[, sum(w)], subset(YX.m5, stage=='t25-i45')[, sum(w)], 1.962),
										a45.35= my.or.from.logit(YX.m5.fit3, 'staget45-i35', 'staget35-i45', subset(YX.m5, stage=='t45-i35')[, sum(w)], subset(YX.m5, stage=='t35-i45')[, sum(w)], 1.962),				
										a35.20= my.or.from.logit(YX.m5.fit3, 'staget35-i20', 'staget20-i35', subset(YX.m5, stage=='t35-i20')[, sum(w)], subset(YX.m5, stage=='t20-i35')[, sum(w)], 1.962),
										a35.25= my.or.from.logit(YX.m5.fit3, 'staget35-i25', 'staget25-i35', subset(YX.m5, stage=='t35-i25')[, sum(w)], subset(YX.m5, stage=='t25-i35')[, sum(w)], 1.962),				
										a25.20= my.or.from.logit(YX.m5.fit3, 'staget25-i20', 'staget20-i25', subset(YX.m5, stage=='t25-i20')[, sum(w)], subset(YX.m5, stage=='t20-i25')[, sum(w)], 1.962)	)
		#        a100.20   a100.25   a100.35   a100.45    a45.20    a45.25    a45.35    a35.20    a35.25     a25.20
		#1: 9.465185e-01 0.4252564 0.8268187 0.9779962 1.5275007 0.4203820 0.8572967 0.4304842 0.7486459 0.52542055
		#2: 1.220221e-06 0.1306335 0.4697377 0.6289970 0.6392072 0.1961459 0.5792472 0.1503698 0.3936640 0.07912175
		#3: 7.342093e+05 1.3843546 1.4553424 1.5206376 3.6502378 0.9009669 1.2688153 1.2324060 1.4237283 3.48913869										
		#	no evidence that the assortativity matrix is asymmetric
	
		#	take a35.35 as reference group
		YX.m5.fit3.or	<- data.table( 	a100.20= my.or.from.logit(YX.m5.fit3, 'staget100-i20', 'staget35-i35', subset(YX.m5, stage=='t100-i20')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),			
										a100.25= my.or.from.logit(YX.m5.fit3, 'staget100-i25', 'staget35-i35', subset(YX.m5, stage=='t100-i25')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),
										a100.35= my.or.from.logit(YX.m5.fit3, 'staget100-i35', 'staget35-i35', subset(YX.m5, stage=='t100-i35')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),
										a100.45= my.or.from.logit(YX.m5.fit3, 'staget100-i45', 'staget35-i35', subset(YX.m5, stage=='t100-i45')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),				
										a100.100= my.or.from.logit(YX.m5.fit3, 'staget100-i100', 'staget35-i35', subset(YX.m5, stage=='t100-i100')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),								
										a45.20= my.or.from.logit(YX.m5.fit3, 'staget45-i20', 'staget35-i35', subset(YX.m5, stage=='t45-i20')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),
										a45.25= my.or.from.logit(YX.m5.fit3, 'staget45-i25', 'staget35-i35', subset(YX.m5, stage=='t45-i25')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),
										a45.35= my.or.from.logit(YX.m5.fit3, 'staget45-i35', 'staget35-i35', subset(YX.m5, stage=='t45-i35')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),				
										a45.45= my.or.from.logit(YX.m5.fit3, 'staget45-i45', 'staget35-i35', subset(YX.m5, stage=='t45-i45')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),			
										a35.20= my.or.from.logit(YX.m5.fit3, 'staget35-i20', 'staget35-i35', subset(YX.m5, stage=='t35-i20')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),
										a35.25= my.or.from.logit(YX.m5.fit3, 'staget35-i25', 'staget35-i35', subset(YX.m5, stage=='t35-i25')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),							
										a25.20= my.or.from.logit(YX.m5.fit3, 'staget25-i20', 'staget35-i35', subset(YX.m5, stage=='t25-i20')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),
										a25.25= my.or.from.logit(YX.m5.fit3, 'staget25-i25', 'staget35-i35', subset(YX.m5, stage=='t25-i25')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),
										a20.20= my.or.from.logit(YX.m5.fit3, 'staget20-i20', 'staget35-i35', subset(YX.m5, stage=='t20-i20')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962)	)		
		#		a100.20   a100.25   a100.35   a100.45  a100.100   	 	a45.20    	a45.25    a45.35    a45.45    a35.20    a35.25    a25.20    a25.25		a20.20
		#1: 	0.6016153 0.7483158 0.9039788 1.0055359 0.7210333 		0.9031965 	0.8730449 0.6598175 0.8890139 0.8659097 0.8966013 0.7861447 1.2735921	1.582723
		#2: 	0.4007082 0.4718837 0.5945772 0.6760679 0.4870311 		0.5725774 	0.5643833 0.4548822 0.6309828 0.5469808 0.5877713 0.5272560 0.8105412	1.035397
		#3:	 	0.9032530 1.1866835 1.3743844 1.4955636 1.0674658 		1.4247226 	1.3505139 0.9570809 1.2525630 1.3707967 1.3676985 1.1721509 2.0011776	2.419375
		#		L*		  L			E		  E			L		  		E			E		  L*		E		  E 		E		  E			E			H*
		
		#	take diagonal as reference group
		YX.m5.fit3.or	<- data.table( 	a100.20= my.or.from.logit(YX.m5.fit3, 'staget100-i20', 'staget100-i100', subset(YX.m5, stage=='t100-i20')[, sum(w)], subset(YX.m5, stage=='t100-i100')[, sum(w)], 1.962),			
										a100.25= my.or.from.logit(YX.m5.fit3, 'staget100-i25', 'staget100-i100', subset(YX.m5, stage=='t100-i25')[, sum(w)], subset(YX.m5, stage=='t100-i100')[, sum(w)], 1.962),
										a100.35= my.or.from.logit(YX.m5.fit3, 'staget100-i35', 'staget100-i100', subset(YX.m5, stage=='t100-i35')[, sum(w)], subset(YX.m5, stage=='t100-i100')[, sum(w)], 1.962),
										a100.45= my.or.from.logit(YX.m5.fit3, 'staget100-i45', 'staget100-i100', subset(YX.m5, stage=='t100-i45')[, sum(w)], subset(YX.m5, stage=='t100-i100')[, sum(w)], 1.962),																																
										a45.20= my.or.from.logit(YX.m5.fit3, 'staget45-i20', 'staget45-i45', subset(YX.m5, stage=='t45-i20')[, sum(w)], subset(YX.m5, stage=='t45-i45')[, sum(w)], 1.962),
										a45.25= my.or.from.logit(YX.m5.fit3, 'staget45-i25', 'staget45-i45', subset(YX.m5, stage=='t45-i25')[, sum(w)], subset(YX.m5, stage=='t45-i45')[, sum(w)], 1.962),
										a45.35= my.or.from.logit(YX.m5.fit3, 'staget45-i35', 'staget45-i45', subset(YX.m5, stage=='t45-i35')[, sum(w)], subset(YX.m5, stage=='t45-i45')[, sum(w)], 1.962),				
										a45.45= my.or.from.logit(YX.m5.fit3, 'staget45-i45', 'staget45-i45', subset(YX.m5, stage=='t45-i45')[, sum(w)], subset(YX.m5, stage=='t45-i45')[, sum(w)], 1.962),													
										a35.20= my.or.from.logit(YX.m5.fit3, 'staget35-i20', 'staget35-i35', subset(YX.m5, stage=='t35-i20')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),
										a35.25= my.or.from.logit(YX.m5.fit3, 'staget35-i25', 'staget35-i35', subset(YX.m5, stage=='t35-i25')[, sum(w)], subset(YX.m5, stage=='t35-i35')[, sum(w)], 1.962),																	
										a25.20= my.or.from.logit(YX.m5.fit3, 'staget25-i20', 'staget25-i25', subset(YX.m5, stage=='t25-i20')[, sum(w)], subset(YX.m5, stage=='t25-i25')[, sum(w)], 1.962)	)		
		#     a100.20   a100.25   a100.35   a100.45    a45.20    a45.25    a45.35    a45.45    a35.20    a35.25    a25.20
		#1: 0.8343793 1.0378381 1.2537268 1.3945763 1.0159532 0.9820375 0.7421903 1.0000000 0.8659097 0.8966013 0.6172657
		#2: 0.4997380 0.5883207 0.7684496 0.8841247 0.6603627 0.6475498 0.5177010 0.7162479 0.5469808 0.5877713 0.2094327
		#3: 1.3931077 1.8308175 2.0454573 2.1997383 1.5630213 1.4893026 1.0640243 1.3961647 1.3707967 1.3676985 1.8192810
		#								  *

		set(YX.m5, NULL, 'stage', YX.m5[, stage.orig])		
		age.cut		<- c(-1, 25, 35, 45, 100)		
		setkey(YX.m5, t.Age, Age)
		for(i in seq_along(age.cut)[-1])
			for(j in seq.int(2, length(age.cut)))
			{
				tmp	<- paste('t',age.cut[i],'-i',age.cut[j],sep='')
				set(YX.m5, YX.m5[, which(stage%in%c('U','Diag') & t.Age<=age.cut[i] & t.Age>age.cut[i-1] & Age<=age.cut[j] & Age>age.cut[j-1])], 'stage', tmp)	
			}
		set(YX.m5, NULL, 'stage', YX.m5[, factor(as.character(stage))])
		YX.m5[, table(stage)]
		
		YX.m5.fit3 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m5)
		summary(YX.m5.fit3)

		YX.m5s			<- subset(YX.m5, !stage%in%c('ART.suA.N','ART.suA.Y','ART.suA.vlNA','DAm','DAy','UAm','UAy'))
		YX.m5s.p		<- YX.m5s[, list(n=round(sum(w),d=1),prop= round(sum(w)/YX.m5s[,sum(w)],d=3) ),by='stage']
		#t25-i25   t25-i45   t25-i35   t25-i100  ART.vlNA  t35-i45   t35-i35   t35-i25   t35-i100  t45-i45   t45-i25   t45-i35   t45-i100 	t100-i25  t100-i35  t100-i45  t100-i100
		#12.1  		4.9 		11.2  	1.7  	2.7 		43.9 		52.3 	26.6 		18.7 	59.4 		19.6 	43.7 		36.5  		7.4 		21.8 	30.0	 33.6
	}
}
######################################################################################
#	wrapper to betareg to come up with better starting values due to scores close to 0 and 1
#	gamlss.BE.limit.l<- c(1e-3, 1e-4, 1e-5, 1e-6, 0)
project.athena.Fisheretal.betareg<- function(YX.m3, formula, include.colnames, gamlss.BE.limit.u=c( seq(0.95, 0.99, 0.01),seq(0.991, 1, 0.001)), gamlss.BE.limit.l=c(0, 0), sigma.formula='~1', verbose=0)
{
	require(gamlss)
	options(warn=2)
	gamlss.BE.limit.i	<- gamlss.BE.limit.u[1]
	YX.m3b				<- copy(YX.m3)
	set(YX.m3b, YX.m3b[, which(score.Y>gamlss.BE.limit.i )], 'score.Y', gamlss.BE.limit.i)
	set(YX.m3b, YX.m3b[, which(score.Y<gamlss.BE.limit.l[1])], 'score.Y', gamlss.BE.limit.l[1])	
	betafit.or			<- gamlss(as.formula(formula), sigma.formula=as.formula(sigma.formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0)
	betafit.rr			<- gamlss(as.formula(formula), sigma.formula=as.formula(sigma.formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0)		#, method=mixed(10,50), control=gamlss.control(c.crit=1e-8, gd.tol=6)			
	#	keep lower limit fixed and try increase upper limit until 1
	tryCatch({
				for(i in seq_along(gamlss.BE.limit.u)[-1])
				{
					#print(gamlss.BE.limit.u[i])			
					YX.m3b	<- copy(YX.m3)
					set(YX.m3b, YX.m3b[, which(score.Y>gamlss.BE.limit.u[i])], 'score.Y', gamlss.BE.limit.u[i])
					set(YX.m3b, YX.m3b[, which(score.Y<gamlss.BE.limit.l[1])], 'score.Y', gamlss.BE.limit.l[1])
					tmp.or				<- gamlss(as.formula(formula), sigma.formula=as.formula(sigma.formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0, start.from=betafit.or)
					tmp.rr				<- gamlss(as.formula(formula), sigma.formula=as.formula(sigma.formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0, start.from=betafit.rr)		#, method=mixed(10,50), control=gamlss.control(c.crit=1e-8, gd.tol=6)
					#print(tmp.rr)
					gamlss.BE.limit		<- gamlss.BE.limit.u[i]	#not run if error in gamlss
					betafit.rr			<- tmp.rr				#not run if error in gamlss
					betafit.or			<- tmp.or				#not run if error in gamlss
				}
				gamlss.init	<- list(start.from.rr=betafit.rr, start.from.or=betafit.or, BE.limit.l=gamlss.BE.limit.l[1], BE.limit.u=gamlss.BE.limit)
			}, error=function(e){ print(e$message); gamlss.init<<- list(start.from.rr=betafit.rr, start.from.or=betafit.or, BE.limit.l=gamlss.BE.limit.l[1], BE.limit.u=gamlss.BE.limit)})
	#	keep highest upper limit that works and try to decrease lower limit 
	if(gamlss.BE.limit.l[1]>0)
	{
		tryCatch({
					for(i in seq_along(gamlss.BE.limit.l)[-1])
					{
						#print(gamlss.BE.limit.l[i])			
						YX.m3b	<- copy(YX.m3)	
						set(YX.m3b, YX.m3b[, which(score.Y>gamlss.init$BE.limit.u)], 'score.Y', gamlss.init$BE.limit.u)					
						set(YX.m3b, YX.m3b[, which(score.Y<gamlss.BE.limit.l[i])], 'score.Y', gamlss.BE.limit.l[i])
						tmp.or				<- gamlss(as.formula(formula), sigma.formula=as.formula(sigma.formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0, start.from=betafit.or )
						tmp.rr				<- gamlss(as.formula(formula), sigma.formula=as.formula(sigma.formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0, start.from=betafit.rr )		#, method=mixed(10,50), control=gamlss.control(c.crit=1e-8, gd.tol=6)
						#print(tmp.rr)
						gamlss.BE.limit		<- gamlss.BE.limit.l[i]	#not run if error in gamlss
						betafit.rr			<- tmp.rr				#not run if error in gamlss
						betafit.or			<- tmp.or				#not run if error in gamlss
					}
					gamlss.init	<- list(start.from.rr=betafit.rr, start.from.or=betafit.or, BE.limit.l=gamlss.BE.limit, BE.limit.u=gamlss.init$BE.limit.u)
				}, error=function(e){ print(e$message); gamlss.init<<- list(start.from.rr=betafit.rr, start.from.or=betafit.or, BE.limit.l=gamlss.BE.limit, BE.limit.u=gamlss.init$BE.limit.u)})				
	}	
	if(verbose) cat(paste('\nsuccess: gamlss beta regression for score<=',gamlss.init$BE.limit.u,'and score>=',gamlss.init$BE.limit.l))
	set(YX.m3, YX.m3[, which(score.Y>gamlss.init$BE.limit.u)], 'score.Y', gamlss.init$BE.limit.u)
	set(YX.m3, YX.m3[, which(score.Y<gamlss.init$BE.limit.l)], 'score.Y', gamlss.init$BE.limit.l)
	betafit.or			<- gamlss(as.formula(formula), sigma.formula=as.formula(sigma.formula), weights=w, data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0, start.from=gamlss.init$start.from.or)
	betafit.rr			<- gamlss(as.formula(formula), sigma.formula=as.formula(sigma.formula), weights=w, data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0, start.from=gamlss.init$start.from.rr)
	options(warn=1)
	list(betafit.or=betafit.or, betafit.rr=betafit.rr, gamlss.BE.limit.u=gamlss.init$BE.limit.u)
}
######################################################################################
project.athena.Fisheretal.betareg.exp<- function(YX.m3, formula, include.colnames, gamlss.BE.limit.u=c( seq(0.95, 0.99, 0.01),seq(0.991, 1, 0.001)), gamlss.BE.limit.l=c(0), verbose=0)
{
	require(gamlss)
	options(warn=2)
	gamlss.BE.limit.i	<- gamlss.BE.limit.u[1]
	YX.m3b				<- copy(YX.m3)
	set(YX.m3b, YX.m3b[, which(score.Y>gamlss.BE.limit.i )], 'score.Y', gamlss.BE.limit.i)
	set(YX.m3b, YX.m3b[, which(score.Y<gamlss.BE.limit.l[1])], 'score.Y', gamlss.BE.limit.l[1])	
	betafit.or			<- gamlss(as.formula(formula), weights=w.orig, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0)
	betafit.rr			<- gamlss(as.formula(formula), weights=w.orig, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0)		#, method=mixed(10,50), control=gamlss.control(c.crit=1e-8, gd.tol=6)		
	
	
	betafit.rr			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0, start.from=betafit.rr)
	
	tmp					<- YX.m3b[, which(score.Y>gamlss.BE.limit.i )]
	set(YX.m3b, tmp, 'score.Y', runif(length(tmp), min=gamlss.BE.limit.i, max=1))
	tmp					<- YX.m3b[, which(score.Y<gamlss.BE.limit.l[1])]
	set(YX.m3b,  tmp, 'score.Y', runif(length(tmp), min=0, max=gamlss.BE.limit.l[1]))
	betafit.or			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0)
	betafit.rr			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0)		#, method=mixed(10,50), control=gamlss.control(c.crit=1e-8, gd.tol=6)
	
	set(YX.m3b, YX.m3b[, which(score.Y>gamlss.BE.limit.i )], 'score.Y', gamlss.BE.limit.i)
	set(YX.m3b, YX.m3b[, which(score.Y<gamlss.BE.limit.l[1])], 'score.Y', gamlss.BE.limit.l[1])
	betafit.or			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0, start.from=betafit.or)
	betafit.rr			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0, start.from=betafit.rr)		#, method=mixed(10,50), control=gamlss.control(c.crit=1e-8, gd.tol=6)		
	
	if(0)
	{
		set(YX.m3b, YX.m3b[, which(score.Y>gamlss.BE.limit.i )], 'score.Y', NA_real_)
		set(YX.m3b, YX.m3b[, which(score.Y<gamlss.BE.limit.l[1])], 'score.Y', NA_real_)		
		betafit.or			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0)
		betafit.rr			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0)		#, control=gamlss.control(c.crit=1e-8, gd.tol=5), i.control=glim.control(glm.trace=TRUE, bf.trace=TRUE)
		YX.m3c				<- copy(YX.m3)
		set(YX.m3c, YX.m3c[, which(score.Y>gamlss.BE.limit.i )], 'score.Y', gamlss.BE.limit.i)
		set(YX.m3c, YX.m3c[, which(score.Y<gamlss.BE.limit.l[1])], 'score.Y', gamlss.BE.limit.l[1])	
		YX.m3c[, mu:=predict( betafit.rr, newdata=na.omit(as.data.frame(YX.m3c[, include.colnames, with=FALSE])), what='mu',  type='response')]
		YX.m3c[, sigma:=predict( betafit.rr, newdata=na.omit(as.data.frame(YX.m3c[, include.colnames, with=FALSE])), what='sigma',  type='response')]
		betafit.rr			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3c[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0, mu.start=YX.m3c[, mu], sigma.start=YX.m3c[, sigma])		#, method=mixed(10,50), control=gamlss.control(c.crit=1e-8, gd.tol=6)		
	}
	#	keep lower limit fixed and try increase upper limit until 1
	tryCatch({
				for(i in seq_along(gamlss.BE.limit.u)[-1])
				{
					print(gamlss.BE.limit.u[i])			
					YX.m3b	<- copy(YX.m3)
					#tmp					<- YX.m3b[, which(score.Y>gamlss.BE.limit.u[i] )]
					#set(YX.m3b, tmp, 'score.Y', runif(length(tmp), min=gamlss.BE.limit.u[i], max=mean(c(gamlss.BE.limit.u[i],1))))
					#tmp					<- YX.m3b[, which(score.Y<gamlss.BE.limit.l[1])]
					#set(YX.m3b,  tmp, 'score.Y', runif(length(tmp), min=mean(c(0,gamlss.BE.limit.l[1])), max=gamlss.BE.limit.l[1]))
					#set(YX.m3b, YX.m3b[, which(score.Y>gamlss.BE.limit.u[i])], 'score.Y', NA_real_)
					#set(YX.m3b, YX.m3b[, which(score.Y<gamlss.BE.limit.l[1])], 'score.Y', NA_real_)
					set(YX.m3b, YX.m3b[, which(score.Y>gamlss.BE.limit.u[i])], 'score.Y', gamlss.BE.limit.u[i])
					set(YX.m3b, YX.m3b[, which(score.Y<gamlss.BE.limit.l[1])], 'score.Y', gamlss.BE.limit.l[1])
					tmp.or				<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0, start.from=betafit.or)
					tmp.rr				<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0, start.from=betafit.rr)		#, method=mixed(10,50), control=gamlss.control(c.crit=1e-8, gd.tol=6)
					print(tmp.rr)
					gamlss.BE.limit		<- gamlss.BE.limit.u[i]	#not run if error in gamlss
					betafit.rr			<- tmp.rr				#not run if error in gamlss
					betafit.or			<- tmp.or				#not run if error in gamlss
				}
				gamlss.init	<- list(start.from.rr=betafit.rr, start.from.or=betafit.or, BE.limit.u=gamlss.BE.limit)
			}, error=function(e){ print(e$message); gamlss.init<<- list(start.from.rr=betafit.rr, start.from.or=betafit.or, BE.limit.u=gamlss.BE.limit)})
	#	keep highest upper limit that works and try to decrease lower limit 
	if(gamlss.BE.limit.l[1]>0)
	{
		tryCatch({
					for(i in seq_along(gamlss.BE.limit.l)[-1])
					{
						#print(gamlss.BE.limit.l[i])			
						YX.m3b	<- copy(YX.m3)	
						tmp					<- YX.m3b[, which(score.Y>gamlss.init$BE.limit.u )]
						set(YX.m3b, tmp, 'score.Y', runif(length(tmp), min=gamlss.init$BE.limit.u, max=mean(c(gamlss.init$BE.limit.u,1))))
						tmp					<- YX.m3b[, which(score.Y<gamlss.BE.limit.l[i])]
						set(YX.m3b,  tmp, 'score.Y', runif(length(tmp), min=mean(c(0,gamlss.BE.limit.l[i])), max=gamlss.BE.limit.l[i]))
						#set(YX.m3b, YX.m3b[, which(score.Y>gamlss.init$BE.limit.u)], 'score.Y', gamlss.init$BE.limit.u)					
						#set(YX.m3b, YX.m3b[, which(score.Y<gamlss.BE.limit.l[i])], 'score.Y', gamlss.BE.limit.l[i])
						tmp.or				<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0, start.from=betafit.or )
						tmp.rr				<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0, start.from=betafit.rr )		#, method=mixed(10,50), control=gamlss.control(c.crit=1e-8, gd.tol=6)
						#print(tmp.rr)
						gamlss.BE.limit		<- gamlss.BE.limit.l[i]	#not run if error in gamlss
						betafit.rr			<- tmp.rr				#not run if error in gamlss
						betafit.or			<- tmp.or				#not run if error in gamlss
					}
					gamlss.init	<- list(start.from.rr=betafit.rr, start.from.or=betafit.or, BE.limit.l=gamlss.BE.limit, BE.limit.u=gamlss.init$BE.limit.u)
				}, error=function(e){ print(e$message); gamlss.init<<- list(start.from.rr=betafit.rr, start.from.or=betafit.or, BE.limit.l=gamlss.BE.limit, BE.limit.u=gamlss.init$BE.limit.u)})				
	}	
	if(verbose) cat(paste('\nsuccess: gamlss beta regression for score<=',gamlss.init$BE.limit.u,'and score>=',gamlss.init$BE.limit.l))
	set(YX.m3, YX.m3[, which(score.Y>gamlss.init$BE.limit.u)], 'score.Y', gamlss.init$BE.limit.u)
	set(YX.m3, YX.m3[, which(score.Y<gamlss.init$BE.limit.l)], 'score.Y', gamlss.init$BE.limit.l)
	betafit.or			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0, start.from=gamlss.init$start.from.or)
	betafit.rr			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0, start.from=gamlss.init$start.from.rr)
	options(warn=1)
	list(betafit.or=betafit.or, betafit.rr=betafit.rr, gamlss.BE.limit.u=gamlss.init$BE.limit.u)
}
######################################################################################
#n	number of random variables to generate
#r	number of observed samples
#p	sampling probability
rznbinom<- function(n, r, p)
{
	stopifnot( all(r>=0), all(p>0), all(p<1))
	ans			<- rep(NA, n)
	tmp			<- r>0
	tmp2		<- which(tmp)
	if(length(tmp2))
	{
		ans[tmp]	<- rnbinom(length(tmp2), r[tmp2], p[tmp2])
	}
	tmp2		<- which(!tmp)
	if(length(tmp2))
	{
		tmp3		<- runif(length(tmp2))
		ans[!tmp]	<- ceiling( log(1-tmp3)/log(1-p[tmp2]) - 1 )
	}
	ans
}
######################################################################################
project.athena.Fisheretal.estimate.risk.core.noWadj<- function(YX.m3, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=1e3, gamlss.BE.required.limit=0.99, gamlss.BE.limit.u=c(0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0, 0), sigma.formula='~1' )
{
	require(gamlss)
	options(warn=0)
	#
	cat(paste('\nbeta regression for formula=',formula))	 
	tmp				<- project.athena.Fisheretal.betareg(YX.m3, formula, include.colnames, gamlss.BE.limit.u=gamlss.BE.limit.u, gamlss.BE.limit.l=gamlss.BE.limit.l, sigma.formula=sigma.formula, verbose=1)
	betafit.rr		<- tmp$betafit.rr
	stopifnot(tmp$gamlss.BE.limit>gamlss.BE.required.limit)
	#
	#	pre-processing
	#
	nt.table	<- copy(X.tables$nt.table)
	set(nt.table, NULL, 'Patient', nt.table[, factor(Patient)])
	set(nt.table, NULL, 'risk', nt.table[, factor(risk)])
	set(nt.table, NULL, 'stat', nt.table[, factor(stat)])	
	nt.table	<- dcast.data.table(nt.table, Patient + risk + factor ~ stat, value.var="nt")
	##	sense check nt.table
	tmp			<- nt.table[, which(X.seq>X.msm)]
	if(length(tmp))	cat(paste('\nWARNING: X.seq>X.msm for entries n=',length(tmp)))
	stopifnot(length(tmp)==0)
	set(nt.table, tmp, 'X.seq', nt.table[tmp, X.msm])	
	tmp			<- nt.table[, which(X.clu>X.seq)]
	if(length(tmp))	cat(paste('\nWARNING: X.clu>X.seq for entries n=',length(tmp)))
	stopifnot(length(tmp)==0)
	set(nt.table, tmp, 'X.clu', nt.table[tmp, X.seq])
	tmp			<- nt.table[, which(YX>X.clu)]
	if(length(tmp))	cat(paste('\nWARNING: YX>X.clu for entries n=',length(tmp)))
	#stopifnot(length(tmp)==0)	
	set(nt.table, tmp, 'YX', nt.table[tmp, X.clu])
	#	make sure all risk factors are in nt.table for every patient (even if zero)
	tmp			<- merge( unique(subset(nt.table, select=c(risk, Patient))), unique(subset(risk.df, select=c(risk,factor))), by='risk', allow.cartesian=TRUE)
	nt.table	<- merge(tmp, nt.table, by=c('risk','factor','Patient'))
	set(nt.table, nt.table[, which(is.na(YX))], c('X.clu','X.msm','X.seq','YX'), 0)	
	#	X.msm not adjusted for censoring
	setnames(nt.table, 'X.msm', 'X.msm.e0')
	#	X.msm adjusted for censoring	
	tmp			<- copy(X.tables$cens.table)
	setkey(tmp, stat, t.period, risk, factor)
	tmp2		<- copy(X.tables$cens.table.bs)
	setkey(tmp2, stat, t.period, risk, factor)
	tmp			<- project.athena.Fisheretal.censoring.model(unique(tmp), unique(tmp2), plot.file=NA )
	tmp			<- merge( nt.table[, list(X.msm.e0=sum(X.msm.e0)), by=c('risk','factor')], subset(tmp, select=c(risk, factor, p.cens)), by=c('risk','factor'))
	#	total censored number of potential transmission intervals 
	tmp			<- tmp[, list( PYe0cpr=round(X.msm.e0/p.cens-X.msm.e0)), by=c('risk','factor')]		
	#	need to allocate the PYe0cpr potential transmission intervals among all Patients with given risk factor for whom yijt>0
	#	simply take mean
	tmp			<- merge(tmp, nt.table[, list(X.msm.e0cp= length(unique(Patient))), by=c('risk','factor') ], by=c('risk','factor')) 
	set(tmp, NULL, 'X.msm.e0cp', round( tmp[, PYe0cpr/X.msm.e0cp] ))
	nt.table	<- merge(nt.table, subset(tmp, select=c(risk, factor, X.msm.e0cp)), by=c('risk','factor'))
	set(nt.table, NULL, 'X.msm.e0cp', nt.table[, X.msm.e0+X.msm.e0cp])	
	#	compute extra variables from nt.table. 
	#	PYs: potential transmission intervals including zero yijt
	#	PTx: prob of non-zero yijt 
	#	X.msm.e0: potential transmission intervals in cohort
	#	X.msm.e0.cp: potential transmission intervals in cohort adjusted for clustering
	#	SX.e0: fraction of potential transmission intervals in cohort that is sampled
	#	SX.e0.cp: fraction of potential transmission intervals in cohort, adjusted for censoring, that is sampled
	adj			<- nt.table[, list(YX=sum(YX), X.clu=sum(X.clu), X.seq=sum(X.seq), X.msm.e0=sum(X.msm.e0), X.msm.e0cp=sum(X.msm.e0cp)), by=c('risk','factor')]
	tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')
	adj[, PYs:= adj[[tmp]]]
	adj[, PTx:= adj[, YX] / adj[[tmp]]]
	adj[, Sx.e0:= PYs/X.msm.e0]					
	adj[, Sx.e0cp:= PYs/X.msm.e0cp]
	#	prepare risk.df and nt.table as needed
	adj.s		<- copy(adj)
	risk.df		<- merge(risk.df, subset(adj.s, select=c(risk, factor, PYs, PTx, X.msm.e0, X.msm.e0cp)), by=c('risk','factor'))
	nt.table	<- merge(nt.table, subset(adj, select=c(risk, factor, Sx.e0, Sx.e0cp)), by=c('risk','factor'))
	#	add PTx.ref to risk.df
	setkey(risk.df, risk, factor)
	tmp			<- subset(unique(risk.df), select=c(risk, factor, PTx))
	setnames(tmp, colnames(tmp), paste(colnames(tmp),'.ref',sep=''))
	risk.df		<- merge(risk.df, tmp, by=c('risk.ref','factor.ref'))
	set(risk.df, NULL, 'risk.ref', risk.df[, as.character(risk.ref)])
	set(risk.df, NULL, 'factor.ref', risk.df[, as.character(factor.ref)])
	set(risk.df, NULL, 'risk', risk.df[, as.character(risk)])
	set(risk.df, NULL, 'factor', risk.df[, as.character(factor)])	
	#	construct bias and censoring adjustments
	adj	<- subset(adj.s, !grepl('U',factor), c(risk, factor, PYs, X.msm.e0cp))[, list(risk=risk, factor=factor, X.msm.e0cp=sum(X.msm.e0cp)*PYs/sum(PYs))]
	adj	<- rbind(adj, subset(adj.s, grepl('U',factor), c(risk, factor, X.msm.e0cp)))		
	adj	<- merge(subset(adj.s, select=c(risk, factor, PTx, PYs, X.msm.e0)), adj, by=c('risk','factor'))	
	adj	<- adj[, list(	risk=risk, factor=factor, cens.e0= 1, cens.e0cp= X.msm.e0cp/PYs	)]		
	tmp	<- subset(adj.s, select=c(risk, factor, PTx, PYs, X.msm.e0, X.msm.e0cp))
	tmp	<- tmp[, list(	risk=risk, factor=factor,
					bias.e0= X.msm.e0/PYs, bias.e0cp= X.msm.e0cp/PYs, 
					P= PYs/sum(PYs), P.raw= PYs/sum(PYs), 
					P.e0= X.msm.e0/sum(X.msm.e0), P.e0cp= X.msm.e0cp/sum(X.msm.e0cp), 
					P.raw.e0= X.msm.e0/sum(X.msm.e0), P.raw.e0cp= X.msm.e0cp/sum(X.msm.e0cp)	)]
	adj	<- merge(adj, tmp, by=c('risk','factor'))	
	#
	#	ready to go:
	#
	#	term-wise risk ratio
	cat(paste('\nterm wise risk ratios based on evidence for transmission ONLY'))
	setkey(risk.df, coef.ref, coef)
	risk.ans	<- subset(risk.df, coef!=coef.ref)[, 	{
							tmp				<- copy(predict.df)
							tmp2			<- copy(predict.df)
							set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
							set(tmp2, NULL, risk, factor(factor.ref, levels=levels( predict.df[[risk]] )))
							#	with corisks E(exp(beta*X)) != exp(beta*E(X)) so it s all a bit more complicated:
							tryCatch({
										tmp		<- predict(betafit.rr, newdata=as.data.frame(tmp), data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), terms=risk, type='terms')
										tmp2	<- predict(betafit.rr, newdata=as.data.frame(tmp2), data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), terms=risk, type='terms')								
										tmp		<- exp( tmp[,risk] - tmp2[,risk] )
									}, error=function(e){ print(e$message); tmp<<- NA_real_ }, warning=function(w){ print(w$message); tmp<<- NA_real_ })
							list(stat= 'RR.term', risk=risk, factor=factor, risk.ref=risk.ref, factor.ref=factor.ref, v=tmp )										
						},by=c('coef','coef.ref')]	
	#	term-wise risk ratio for PDT including zeros	
	#	(this is not really valid, because the variance term for the mixed regression model contains terms that come from the zeros)
	tmp			<- merge( subset(risk.ans, stat=='RR.term'), subset(risk.df, select=c(risk, factor, risk.ref, factor.ref, PTx, PTx.ref)), by=c('risk','risk.ref','factor','factor.ref'))
	set(tmp, NULL, 'v', tmp[, v*PTx/PTx.ref])
	set(tmp, NULL, 'stat', 'RR.term.ptx')
	risk.ans	<- rbind( risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)) )		
	#	potential transmission intervals in denominator pop (either X.seq or X.clu)
	cat(paste('\nintervals with evidence for or against transmission (denominator)'))
	setkey(risk.df, coef)
	tmp			<- unique(risk.df)
	tmp[, stat:='PYs']
	tmp			<- subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, PYs))
	set(tmp, NULL, c('coef.ref','risk.ref', 'factor.ref'), 'None')
	setnames(tmp,'PYs','v')
	risk.ans	<- rbind(risk.ans, tmp)
	#
	#	number of transmissions and proportions by raw count
	#
	cat(paste('\nnumber and proportion of transmissions using raw evidence for transmission (y_ijt*w_ijt)'))
	missing		<- merge(nt.table, unique( subset( risk.df, select=c(risk, factor, PTx) ) ), by=c('risk','factor'))	
	# compute the mean Y by stage for non-zero Y
	tmp			<- YX.m3[, list(yYX.mean= mean(score.Y), yYX.med= median(score.Y)), by='stage']
	setnames(tmp, 'stage','factor')
	tmp[, risk:='stage']	
	missing		<- merge(missing, tmp, by=c('risk','factor'))
	#compute the mean Y by stage including zero Y
	set(missing, NULL, 'yYX.mean', missing[, yYX.mean*PTx])
	set(missing, NULL, 'yYX.med', missing[, yYX.med*PTx])
	#	compute the sum of observed Y's by stage for each recipient
	tmp			<- YX.m3[, list(yYX.sum= sum(score.Y), YX.n=length(score.Y), YX.w=w.t[1]), by=c('stage','Patient')]
	setnames(tmp, 'stage','factor')
	tmp[, risk:='stage']	
	missing		<- merge(missing, tmp, by=c('risk','factor','Patient'), all.x=TRUE)
	set(missing, missing[, which(is.na(YX.n))], 'YX.w', 1.) 
	set(missing, missing[, which(is.na(YX.n))], c('yYX.sum','YX.n'), 0.)
	#	expected missing potential transmission intervals for recipient j and factor x
	tmp			<- melt( subset(missing, select=c(Patient, risk, factor, X.seq, Sx.e0, Sx.e0cp)), measure.vars=c('Sx.e0','Sx.e0cp'), variable.name='Sx.method', value.name='Sx' )
	set( tmp, tmp[, which(X.seq==0)], 'X.seq', 1)	#mean of sampling for no observed intervals coincides with sampling model for 1 observed interval 
	tmp[, YXm.e:= tmp[, as.integer(round( X.seq*(1-Sx)/Sx )) ]]
	set(tmp, NULL, 'Sx.method', tmp[, gsub('Sx','YXm.e',Sx.method)])
	missing		<- merge(missing, dcast.data.table(tmp, Patient + risk + factor ~ Sx.method, value.var="YXm.e"), by=c('Patient','risk','factor'))
	#	calculate prob Pj(x) that recipient j got infected from x  - with and without adjustment
	#tmp		<- missing[, 	list(	factor=factor,  Pjx= yYX.sum*YX.w, Pjx.e0= (yYX.sum+YXm.e.e0*yYX.mean)*YX.w, Pjx.e0cp= (yYX.sum+YXm.e.e0cp*yYX.mean)*YX.w, coef=paste(risk,as.character(factor), sep='')), by=c('risk','Patient')]
	tmp			<- missing[, 	list(	factor=factor, 	
										Pjx= yYX.sum*YX.w/sum(yYX.sum*YX.w), 
										#Pjx.e0= (yYX.sum+YXm.e.e0*yYX.mean)*YX.w/sum((yYX.sum+YXm.e.e0*yYX.mean)*YX.w),
										Pjx.e0= (yYX.sum+YXm.e.e0*yYX.med)*YX.w/sum((yYX.sum+YXm.e.e0*yYX.med)*YX.w),
										#Pjx.e0cp= (yYX.sum+YXm.e.e0cp*yYX.mean)*YX.w/sum((yYX.sum+YXm.e.e0cp*yYX.mean)*YX.w),
										Pjx.e0cp= (yYX.sum+YXm.e.e0cp*yYX.med)*YX.w/sum((yYX.sum+YXm.e.e0cp*yYX.med)*YX.w),
										coef=paste(risk,as.character(factor), sep='')), by=c('risk','Patient')]
	#	exclude recipients with no evidence for direct transmission
	tmp			<- subset(tmp, !is.nan(Pjx))					
	#	compute N.raw etc	
	tmp			<- tmp[, list(	N.raw= sum(Pjx), N.raw.e0= sum(Pjx.e0), N.raw.e0cp=sum(Pjx.e0cp), 
					risk.ref='None', factor.ref='None', coef.ref='None', coef=coef[1]), by=c('risk','factor')]		
	tmp			<- melt(tmp, id.vars=c('coef','risk','factor','coef.ref','risk.ref','factor.ref'), variable.name='stat', value.name = "v")
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#compute proportions from the various N.raw
	tmp			<- tmp[, list(v= v/sum(v, na.rm=TRUE), risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef=coef, coef.ref='None') , by='stat']
	set(tmp, NULL, 'stat', tmp[, gsub('N','P',stat)])
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#
	#	number of transmissions and proportion of transmissions from regression
	#
	cat(paste('\nnumber and proportion of transmissions from regression'))
	setkey(risk.df, risk, factor)
	tmp			<- unique(risk.df)[, 		{
				tmp					<- copy(predict.df)
				set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
				corisk				<- intersect(colnames(risk.df), colnames(predict.df))
				if(!length(corisk))	
					corisk<- NULL				
				tmp					<- na.omit(merge(tmp[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3[which( YX.m3[, risk, with=FALSE][[1]]==factor ), c(risk,corisk,'w'), with=FALSE], by=risk))
				setkey(tmp, NULL)
				if(nrow(tmp))	
				{
					
					if( nrow(unique(tmp[, setdiff(colnames(tmp),'w'), with=FALSE]))==1 )
						tmp			<- tmp[1,]		
					tryCatch({
								tmp[, predict:= predict(betafit.rr, newdata=as.data.frame(tmp), data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), type='link')]				
								tmp				<- weighted.mean(exp(tmp[,predict]), w=tmp[,w], na.rm=TRUE)
							}, error=function(e){ print(e$message); tmp<<- NA_real_ }, warning=function(w){ print(w$message); tmp<<- NA_real_ })
				}
				else
					tmp				<- NA_real_	
				tmp[2]				<- sum( YX.m3[ which(as.character(YX.m3[, risk, with=FALSE][[1]])==factor), w] )
				list(risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef.ref='None', stat= 'N', expbeta=ifelse(tmp[2]<2*EPS, 0., tmp[1]), n=tmp[2] )
			},by= 'coef']
	##	use expbeta instead of yijt in calculations of N and P
	missing		<- merge( missing, subset(tmp, select=c(risk, factor, expbeta)), by=c('risk','factor'))
	#	calculate prob Pj(x) that recipient j got infected from x  - with and without adjustment 
	tmp			<- missing[, 	list(	factor=factor, 	
										Pjx= YX*expbeta*YX.w/sum(YX*expbeta*YX.w), 
										Pjx.e0= (YX*expbeta+YXm.e.e0*expbeta*PTx)*YX.w/sum((YX*expbeta+YXm.e.e0*expbeta*PTx)*YX.w),
										Pjx.e0cp= (YX*expbeta+YXm.e.e0cp*expbeta*PTx)*YX.w/sum((YX*expbeta+YXm.e.e0cp*expbeta*PTx)*YX.w),
										coef=paste(risk,as.character(factor), sep='')), by=c('risk','Patient')]
	#	exclude recipients with no evidence for direct transmission
	tmp			<- subset(tmp, !is.nan(Pjx))										
	tmp			<- tmp[, list(	N= sum(Pjx), N.e0= sum(Pjx.e0), N.e0cp=sum(Pjx.e0cp), 
					risk.ref='None', factor.ref='None', coef.ref='None', coef=coef[1]), by=c('risk','factor')]		
	tmp			<- melt(tmp, id.vars=c('coef','risk','factor','coef.ref','risk.ref','factor.ref'), variable.name='stat', value.name = "v")
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#	compute proportions from the various N
	tmp			<- tmp[, list(v= v/sum(v, na.rm=TRUE), risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef=coef, coef.ref='None') , by='stat']	
	set(tmp, NULL, 'stat', tmp[, gsub('N','P',stat)])
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#	relative transmissibilities
	cat(paste('\nrelative transmissibilities'))	
	tmp			<- subset( melt(adj, id.vars=c('risk','factor'), variable.name='stat', value.name='p'), !grepl('cens',stat) & !grepl('^bias',stat)	)	
	tmp			<- merge( subset(risk.ans, grepl('P.',stat) | stat=='P' ), tmp, by=c('risk','factor','stat'))
	set(tmp, NULL, 'v', tmp[, v/p])
	set(tmp, tmp[,which(is.nan(v) | v==0)],'v',NA_real_)	
	set(tmp, NULL, 'stat', tmp[, gsub('P','RI',stat)])
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#	risk ratio from relative transmissibilities
	cat(paste('\n term-wise risk ratio from identity RR(x,z)=RT(x)/RT(z).'))	
	setkey(risk.df, coef.ref, coef)
	tmp			<- subset(risk.df, coef!=coef.ref, select=c(coef.ref, coef))
	tmp.ref		<- merge(data.table(coef=tmp[, unique(coef.ref)]), subset(risk.ans, grepl('RI',stat), select=c(coef, risk, factor, stat, v)), by='coef')
	setnames(tmp.ref, c('coef','risk','factor','v'), c('coef.ref','risk.ref','factor.ref','v.ref'))
	tmp			<- merge(tmp, subset(risk.ans, grepl('RI',stat), select=c(coef, risk, factor, stat, v)), by='coef', allow.cartesian=TRUE)
	tmp			<- merge(tmp.ref, tmp, by=c('stat','coef.ref'))
	tmp.ref		<- NULL
	set(tmp, NULL, 'v', tmp[, v/v.ref])
	set(tmp, NULL, 'stat', tmp[, gsub('RI','RR',stat)])
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#	Bootstraps	on ratio and on prob
	cat(paste('\nregression on bootstrap data sets bs.n=',bs.n))
	tmp			<- lapply(seq_len(bs.n), function(bs.i)
			{
				if(bs.i%%100==0)	cat(paste('\nregression on bootstrap data sets bs.i=',bs.i))
				#bootstrap over recently infected Patient
				tmp					<- unique(subset(YX.m3, select=Patient))					
				#
				bs.repeat			<- 1
				while(bs.repeat)
				{
					YX.m3.bs			<- tmp[ sample( seq_len(nrow(tmp)), nrow(tmp), replace=TRUE ), ]
					#recipient MSM are not necessarily unique any longer - need to create unique bs id
					YX.m3.bs[, Patient.bs:=paste(Patient, seq_len(nrow(tmp)),sep='_bs' )]
					YX.m3.bs			<- merge( YX.m3, YX.m3.bs, by='Patient', allow.cartesian=TRUE )
					tryCatch({ 
								tmp2			<- project.athena.Fisheretal.betareg(YX.m3.bs, formula, include.colnames, gamlss.BE.limit.u=gamlss.BE.limit.u, gamlss.BE.limit.l=gamlss.BE.limit.l, sigma.formula=sigma.formula, verbose=1 )								
								betafit.rr.bs	<- tmp2$betafit.rr
								bs.repeat		<- ifelse( tmp2$gamlss.BE.limit>gamlss.BE.required.limit, 0, 1 )
							}, error=function(e){ bs.repeat<<- 1})																
				}					
				#	odds ratio and risk ratio
				setkey(risk.df, coef.ref, coef)
				#term wise risk ratio
				risk.ans.bs					<- subset(risk.df, coef!=coef.ref)[, 	{
							tmp				<- copy(predict.df)
							tmp2			<- copy(predict.df)
							set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
							set(tmp2, NULL, risk, factor(factor.ref, levels=levels( predict.df[[risk]] )))
							#	with corisks E(exp(beta*X)) != exp(beta*E(X)) so it s all a bit more complicated:
							tryCatch({
										tmp		<- predict(betafit.rr.bs, newdata=as.data.frame(tmp), data=na.omit(as.data.frame(YX.m3.bs[, include.colnames, with=FALSE])), terms=risk, type='terms')
										tmp2	<- predict(betafit.rr.bs, newdata=as.data.frame(tmp2), data=na.omit(as.data.frame(YX.m3.bs[, include.colnames, with=FALSE])), terms=risk, type='terms')								
										tmp		<- exp( tmp[,risk] - tmp2[,risk] )
									}, error=function(e){ print(e$message); tmp<<- NA_real_ }, warning=function(w){ print(w$message); tmp<<- NA_real_ })
							list(stat= 'RR.term', risk=risk, factor=factor, risk.ref=risk.ref, factor.ref=factor.ref, v=tmp )										
						},by=c('coef','coef.ref')]
				#
				#	raw number of transmissions
				#				
				missing		<- merge(nt.table, unique( subset( risk.df, select=c(risk, factor, PTx) ) ), by=c('risk','factor'))
				#	reduce to bootstrap sampled recipients
				tmp			<- subset(YX.m3.bs, select=c(Patient, Patient.bs))		
				setkey(tmp, Patient, Patient.bs)
				missing		<- merge(unique(tmp), missing, by='Patient', allow.cartesian=TRUE)
				# compute the mean Y by stage for non-zero Y
				tmp			<- YX.m3[, list(yYX.mean= mean(score.Y)), by='stage']
				setnames(tmp, 'stage','factor')
				tmp[, risk:='stage']	
				missing		<- merge(missing, tmp, by=c('risk','factor'))
				#compute the mean Y by stage including zero Y
				set(missing, NULL, 'yYX.mean', missing[, yYX.mean*PTx])
				#	compute the sum of observed Y's by stage for each recipient
				tmp			<- YX.m3.bs[, list(Patient=Patient[1], yYX.sum= sum(score.Y), YX.bs=length(score.Y), YX.w=w.t[1]), by=c('stage','Patient.bs')]
				setnames(tmp, 'stage','factor')
				tmp[, risk:='stage']	
				missing		<- merge(missing, tmp, by=c('risk','factor','Patient','Patient.bs'), all.x=TRUE)
				set(missing, missing[, which(is.na(YX.bs))], 'YX.w', 1.) 
				set(missing, missing[, which(is.na(YX.bs))], c('yYX.sum','YX.bs'), 0.)
				#draw number missing potential transmission intervals from neg binomial
				tmp			<- melt( subset(missing, select=c(Patient.bs, risk, factor, X.seq, PTx, Sx.e0, Sx.e0cp)), measure.vars=c('Sx.e0','Sx.e0cp'), variable.name='Sx.method', value.name='Sx' )				
				tmp[, YXm.r:= rznbinom(nrow(tmp), tmp[, X.seq], tmp[, Sx])]
				#tmp[, YXm.r:= tmp[, round((1-Sx)*X.seq/Sx)]]
				set(tmp, NULL, 'Sx.method', tmp[, gsub('Sx','YXm.r',Sx.method)])
				#draw number of missing non-zero potential transmission intervals under PTx
				set(tmp, NULL, 'YXm.r', rbinom(nrow(tmp), tmp[,YXm.r], tmp[,PTx]))
				#set(tmp, NULL, 'YXm.r', tmp[,round(YXm.r*PTx)])
				missing		<- merge(missing, dcast.data.table(tmp, Patient.bs + risk + factor ~ Sx.method, value.var="YXm.r"), by=c('Patient.bs','risk','factor'), all.x=1)				
				#draw missing scores from all yijt in that stage	for number missing YXm.r.e0
				set(missing, NULL, 'risk', missing[, as.character(risk)])
				set(missing, NULL, 'factor', missing[, as.character(factor)])
				tmp			<- missing[, 	{
												z	<- median( YX.m3[ which( YX.m3[[risk]]==factor ), ][['score.Y']] )	
												#z	<- YX.m3[ which( YX.m3[[risk]]==factor ), ][['score.Y']]	#this is on purpose YX.m3 instead of YX.m3.bs to make sure that we have scores for every factor
												#z	<- c(z, rep(0, length(z)/PTx[1]-length(z) ))
												#list(Patient.bs=rep(Patient.bs, YXm.r.e0), yYXm.r.e0=sample(c(0, z), sum(YXm.r.e0), prob=c(length(z)/PTx[1]-length(z), rep(1, length(z))), replace=TRUE)  )
												list(Patient.bs=rep(Patient.bs, YXm.r.e0), yYXm.r.e0=sample(z, sum(YXm.r.e0), replace=TRUE)  )
											}, by=c('risk','factor')]
				#subset(tmp, Patient.bs=='M37593_bs131')
				#hist( subset(tmp, factor=='ART.suA.Y.4')[, yYXm.r.e0], breaks=100)									
				missing		<- merge(missing, tmp[, list(yYXm.sum.e0=sum(yYXm.r.e0)), by=c('Patient.bs','risk','factor')], by=c('Patient.bs','risk','factor'), all.x=TRUE)
				#draw missing scores from all yijt in that stage	for number missing YXm.r.e0cp
				tmp			<- missing[, 	{
												z	<- median( YX.m3[ which( YX.m3[[risk]]==factor ), ][['score.Y']] )
												#z	<- YX.m3[ which( YX.m3[[risk]]==factor ), ][['score.Y']]
												#z	<- c(z, rep(0, length(z)/PTx[1]-length(z) ))
												#list(Patient.bs=rep(Patient.bs, YXm.r.e0cp), yYXm.r.e0cp=sample(c(0,z), sum(YXm.r.e0cp), prob=c(length(z)/PTx[1]-length(z), rep(1, length(z))), replace=TRUE)  )
												list(Patient.bs=rep(Patient.bs, YXm.r.e0cp), yYXm.r.e0cp=sample(z, sum(YXm.r.e0cp), replace=TRUE)  )
											}, by=c('risk','factor')]
				missing		<- merge(missing, tmp[, list(yYXm.sum.e0cp=sum(yYXm.r.e0cp)), by=c('Patient.bs','risk','factor')], by=c('Patient.bs','risk','factor'), all.x=TRUE)
				set(missing, missing[, which(is.na(yYXm.sum.e0))], 'yYXm.sum.e0', 0.)
				set(missing, missing[, which(is.na(yYXm.sum.e0cp))], 'yYXm.sum.e0cp', 0.)
				#	calculate prob Pj(x) that recipient j got infected from x  - with and without adjustment				
				tmp			<- missing[, 	list(	factor=factor, 	
													Pjx= yYX.sum*YX.w/sum(yYX.sum*YX.w), 
													Pjx.e0= (yYX.sum+yYXm.sum.e0)*YX.w/sum((yYX.sum+yYXm.sum.e0)*YX.w),
													Pjx.e0cp= (yYX.sum+yYXm.sum.e0cp)*YX.w/sum((yYX.sum+yYXm.sum.e0cp)*YX.w),
													coef=paste(risk,as.character(factor), sep='')), by=c('risk','Patient.bs')]
				#	exclude recipients with no evidence for direct transmission
				tmp			<- subset(tmp, !is.nan(Pjx))					
				#various N.raw									
				tmp			<- tmp[, list(	N.raw= sum(Pjx), N.raw.e0= sum(Pjx.e0), N.raw.e0cp=sum(Pjx.e0cp), 
								risk.ref='None', factor.ref='None', coef.ref='None', coef=coef[1]), by=c('risk','factor')]		
				tmp			<- melt(tmp, id.vars=c('coef','risk','factor','coef.ref','risk.ref','factor.ref'), variable.name='stat', value.name = "v")				
				risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
				#
				#	number of transmissions from regression
				#						
				setkey(risk.df, coef)
				tmp			<- unique(risk.df)[, 		{
							tmp				<- copy(predict.df)
							set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
							corisk			<- intersect(colnames(risk.df), colnames(predict.df))
							if(!length(corisk))	
								corisk<- NULL
							tmp				<- na.omit(merge(tmp[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3.bs[which( YX.m3.bs[, risk, with=FALSE][[1]]==factor ), c(risk,corisk,'w'), with=FALSE], by=risk))
							setkey(tmp, NULL)
							if(nrow(tmp))	
							{
								if( nrow(unique(tmp[, setdiff(colnames(tmp),'w'), with=FALSE]))==1 )
									tmp		<- tmp[1,]																			
								tryCatch({
											tmp[, predict:= predict(betafit.rr.bs, newdata=as.data.frame(tmp), data=na.omit(as.data.frame(YX.m3.bs[, include.colnames, with=FALSE])), type='link')]							
											tmp			<- weighted.mean(exp(tmp[,predict]), w=tmp[,w], na.rm=TRUE)
										}, error=function(e){ print(e$message); tmp<<- NA_real_ }, warning=function(w){ print(w$message); tmp<<- NA_real_ })
							}
							else
								tmp			<- NA_real_
							tmp[2]			<- sum( YX.m3.bs[ which(as.character(YX.m3.bs[, risk, with=FALSE][[1]])==factor), w] )																												
							list(risk=risk, factor=factor, expbeta=ifelse(tmp[2]<2*EPS, 0., tmp[1]) )
						},by= 'coef']
				missing		<- merge( missing, subset(tmp, select=c(risk, factor, expbeta)), by=c('risk','factor'))
				#compute probability Pjx for bootstrap sampled recipient MSM j				
				tmp			<- missing[, 	list(	factor=factor, 	
													Pjx= YX.bs*expbeta*YX.w/sum(YX.bs*expbeta*YX.w), 
													Pjx.e0= (YX.bs+YXm.r.e0)*expbeta*YX.w/sum((YX.bs+YXm.r.e0)*expbeta*YX.w),
													Pjx.e0cp= (YX.bs+YXm.r.e0cp)*expbeta*YX.w/sum((YX.bs+YXm.r.e0cp)*expbeta*YX.w),
													coef=paste(risk,as.character(factor), sep='')), by=c('risk','Patient.bs')]
				tmp			<- subset(tmp, !is.nan(Pjx))						
				#various N									
				tmp			<- tmp[, list(	N= sum(Pjx), N.e0= sum(Pjx.e0), N.e0cp=sum(Pjx.e0cp), 
								risk.ref='None', factor.ref='None', coef.ref='None', coef=coef[1]), by=c('risk','factor')]		
				tmp			<- melt(tmp, id.vars=c('coef','risk','factor','coef.ref','risk.ref','factor.ref'), variable.name='stat', value.name = "v")				
				risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
				#
				#	denominator population: y_ijt=0 AND y_ijt>0
				#
				tmp			<- subset(risk.ans, stat=='PYs')
				set(tmp, NULL, 'v', rbinom(nrow(tmp), tmp[, sum(v)], tmp[, v/sum(v)]))
				risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
				#	sampling probabilities
				tmp			<- risk.df[,	{
												z	<- table( YX.m3.bs[, risk, with=FALSE])
												list(factor=rownames(z), n=as.numeric(unclass(z)))												
											}, by='risk']
				tmp			<- merge( subset(risk.ans.bs, stat=='PYs'), tmp, by=c('risk','factor'))		
				set(tmp, NULL, 'v', tmp[, n/v])
				set(tmp, NULL, 'stat', 'PTx')
				risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
				risk.ans.bs[, bs:=bs.i]
			})
	risk.ans.bs	<- do.call('rbind',tmp)
	#
	#	compute adjustments for each bootstrap sample: collect PYe0 etc to compute sampling proportion relative to bootstrap PYs
	#
	tmp			<- dcast.data.table( subset(risk.ans.bs, stat=='PYs' | stat=='PTx', select=c(risk,factor,v,bs,stat)), risk+factor+bs ~ stat, value.var='v' )
	tmp			<- merge( tmp, subset(unique(risk.df), select=c(risk, factor, X.msm.e0, X.msm.e0cp)), by=c('risk','factor'))
	tmp2		<- tmp[, list(  	factor=factor, 
									X.msm.e0=rbinom(length(X.msm.e0), sum(X.msm.e0), X.msm.e0/sum(X.msm.e0)),
									X.msm.e0cp=rbinom(length(X.msm.e0cp), sum(X.msm.e0cp), X.msm.e0cp/sum(X.msm.e0cp))	), by=c('bs','risk')]
	tmp			<- merge( subset(tmp, select=c(risk, factor, bs, PYs, PTx)), tmp2, by=c('bs','risk','factor') )	
	adj.bs		<- tmp[, list(	factor=factor,
								P= PYs/sum(PYs), P.raw= PYs/sum(PYs), 
								P.e0= X.msm.e0/sum(X.msm.e0), P.e0cp= X.msm.e0cp/sum(X.msm.e0cp),
								P.raw.e0= X.msm.e0/sum(X.msm.e0), P.raw.e0cp= X.msm.e0cp/sum(X.msm.e0cp)	), by=c('bs','risk')]
	#
	# 	add Proportions
	#
	tmp			<- subset(risk.ans.bs, grepl('N',stat))	
	tmp			<- tmp[, list(v= v/sum(v, na.rm=TRUE), risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef=coef, coef.ref='None') , by=c('stat','bs')]	
	set(tmp, NULL, 'stat', tmp[, gsub('N','P',stat)])
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v, bs)))
	#
	#	add RR.term
	#	
	tmp			<- subset(risk.ans.bs, stat=='RR.term')	
	tmp2		<- subset(risk.ans.bs, grepl('PTx',stat), select=c(risk, factor, stat, v, bs))		
	tmp2		<- dcast.data.table(tmp2, bs + risk + factor ~ stat, value.var="v")
	setnames(tmp2, c('risk','factor','PTx'), c('risk.ref','factor.ref','PTx.ref'))
	tmp			<- merge(tmp, tmp2, by=c('bs','risk.ref','factor.ref'))	
	tmp2		<- subset(risk.ans.bs, grepl('PTx',stat), select=c(risk, factor, stat, v, bs))		
	tmp2		<- dcast.data.table(tmp2, bs + risk + factor ~ stat, value.var="v")
	tmp			<- merge(tmp, tmp2, by=c('bs','risk','factor'))
	set(tmp, NULL, 'v', tmp[, v*PTx/PTx.ref])
	set(tmp, NULL, 'stat', 'RR.term.ptx')
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, v, bs)))
	#	set missing values in N P to zero
	set(risk.ans.bs, risk.ans.bs[, which(is.na(v) & (grepl('P',stat) | grepl('PYs',stat) | grepl('N',stat)))], 'v', 0.)
	#
	#	add RIs
	#	
	tmp			<- melt(adj.bs, id.vars=c('risk','factor','bs'), variable.name='stat', value.name='p')
	tmp2		<- subset(risk.ans.bs, grepl('P',stat) & !grepl('PT',stat) & !grepl('PY',stat) )
	tmp			<- merge(tmp2, tmp, by=c('risk','factor','stat','bs'))
	set(tmp, NULL, 'v', tmp[, v/p])
	set(tmp, tmp[,which(is.nan(v) | v==0)],'v',NA_real_)	
	set(tmp, NULL, 'stat', tmp[, gsub('P','RI',stat)])
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v, bs)))
	#
	#	Risk ratio from relative transmissibilities
	#
	tmp			<- subset(risk.df, coef!=coef.ref, select=c(coef.ref, coef))
	tmp2		<- merge(data.table(coef=tmp[, unique(coef.ref)]), subset(risk.ans.bs, grepl('RI',stat), select=c(bs, coef, risk, factor, stat, v)), by='coef')
	setnames(tmp2, c('coef','risk','factor','v'), c('coef.ref','risk.ref','factor.ref','v.ref'))
	tmp			<- merge(tmp, subset(risk.ans.bs, grepl('RI',stat), select=c(bs, coef, risk, factor, stat, v)), by='coef', allow.cartesian=TRUE)
	tmp			<- merge(tmp2, tmp, by=c('bs','stat','coef.ref'))
	tmp2		<- NULL
	set(tmp, NULL, 'v', tmp[, v/v.ref])
	set(tmp, NULL, 'stat', tmp[, gsub('RI','RR',stat)])
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v, bs)))
	#
	tmp			<- risk.ans.bs[,	list(l95.bs=quantile(v, prob=0.025, na.rm=TRUE), u95.bs=quantile(v, prob=0.975, na.rm=TRUE)), by=c('coef','coef.ref','stat')]
	risk.ans	<- merge(risk.ans, tmp, by=c('coef','coef.ref','stat'), all.x=TRUE)
	#	
	setkey(risk.ans, stat, coef.ref, coef)
	list(risk=risk.ans, fit.rr=betafit.rr, risk.bs=risk.ans.bs)
}
######################################################################################
project.athena.Fisheretal.estimate.risk.core.Wallinga<- function(YX.m3, YXf, X.tables, method.risk, risk.df, bs.n=1e3, use.YXf=TRUE )
{
	options(warn=0)
	#
	#	pre-processing
	#
	nt.table	<- copy(X.tables$nt.table)
	set(nt.table, NULL, 'Patient', nt.table[, factor(Patient)])
	set(nt.table, NULL, 'risk', nt.table[, factor(risk)])
	set(nt.table, NULL, 'stat', nt.table[, factor(stat)])	
	nt.table	<- dcast.data.table(nt.table, Patient + risk + factor ~ stat, value.var="nt")
	##	sense check nt.table
	tmp			<- nt.table[, which(X.seq>X.msm)]
	if(length(tmp))	cat(paste('\nWARNING: X.seq>X.msm for entries n=',length(tmp)))
	stopifnot(length(tmp)==0)
	set(nt.table, tmp, 'X.seq', nt.table[tmp, X.msm])	
	tmp			<- nt.table[, which(X.clu>X.seq)]
	if(length(tmp))	cat(paste('\nWARNING: X.clu>X.seq for entries n=',length(tmp)))
	stopifnot(length(tmp)==0)
	set(nt.table, tmp, 'X.clu', nt.table[tmp, X.seq])
	tmp			<- nt.table[, which(YX>X.clu)]
	if(length(tmp))	cat(paste('\nWARNING: YX>X.clu for entries n=',length(tmp)))
	#	check recipients
	tmp			<- nt.table[, list(S=any(YX>0)), by='Patient']
	tmp			<- setdiff(unique(subset(YX.m3, select=Patient))$Patient, subset(tmp, S, Patient)$Patient)
	if(length(tmp))	cat(paste('\nWARNING: recipients not in nt.table=',paste(tmp, collapse=' ')))				
	#stopifnot(length(tmp)==0)	
	set(nt.table, tmp, 'YX', nt.table[tmp, X.clu])
	#	make sure all risk factors are in nt.table for every patient (even if zero)
	tmp			<- merge( unique(subset(nt.table, select=c(risk, Patient))), unique(subset(risk.df, select=c(risk,factor))), by='risk', allow.cartesian=TRUE)
	nt.table	<- merge(tmp, nt.table, by=c('risk','factor','Patient'), all.x=TRUE)
	set(nt.table, nt.table[, which(is.na(YX))], c('X.clu','X.msm','X.seq','YX'), 0)	
	#	X.msm not adjusted for censoring
	setnames(nt.table, 'X.msm', 'X.msm.e0')
	#	X.msm adjusted for censoring	
	tmp			<- copy(X.tables$cens.table)
	setkey(tmp, stat, t.period, risk, factor)
	ct			<- unique(tmp)
	tmp2		<- copy(X.tables$cens.table.bs)
	setkey(tmp2, stat, t.period, risk, factor)
	ctb			<- unique(tmp2)
	tmp			<- project.athena.Fisheretal.censoring.model(ct, ctb, plot.file=NA )
	ct			<- copy(tmp$ctn)
	tmp			<- ct[, list(n.adj.med=round(median(n.adj))), by=c('risk','factor','t.period')]
	tmp			<- merge( nt.table[, list(X.msm.e0=sum(X.msm.e0)), by=c('risk','factor')], tmp, by=c('risk','factor'))
	stopifnot( tmp[, all(X.msm.e0<=n.adj.med)] )
	#	total number of right censored potential transmission intervals 
	tmp			<- tmp[, list( PYe0cpr=round(n.adj.med-X.msm.e0)), by=c('risk','factor')]		
	#	need to allocate the PYe0cpr potential transmission intervals among all Patients with given risk factor for whom yijt>0
	#	the total is for all recipient MSM, not just those with a likely transmitter
	#	compute the total per recipient MSM (irrespective of whether the recipient has a likely transmitter)
	tmp			<- merge(tmp, nt.table[, list(X.msm.e0cp= length(unique(Patient))), by=c('risk','factor') ], by=c('risk','factor')) 
	set(tmp, NULL, 'X.msm.e0cp', round( tmp[, PYe0cpr/X.msm.e0cp] ))
	nt.table	<- merge(nt.table, subset(tmp, select=c(risk, factor, X.msm.e0cp)), by=c('risk','factor'))
	set(nt.table, NULL, 'X.msm.e0cp', nt.table[, X.msm.e0+X.msm.e0cp])	
	#	compute extra variables from nt.table. 
	#	PYs: potential transmission intervals including zero yijt
	#	PTx: prob of non-zero yijt 
	#	X.msm.e0: potential transmission intervals in cohort
	#	X.msm.e0.cp: potential transmission intervals in cohort adjusted for clustering
	#	SX.e0: fraction of potential transmission intervals in cohort that is sampled
	#	SX.e0.cp: fraction of potential transmission intervals in cohort, adjusted for censoring, that is sampled
	adj			<- nt.table[, list(YX=sum(YX), X.clu=sum(X.clu), X.seq=sum(X.seq), X.msm.e0=sum(X.msm.e0), X.msm.e0cp=sum(X.msm.e0cp), nRec=length(unique(Patient))), by=c('risk','factor')]
	tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')
	adj[, PYs:= adj[[tmp]]]
	adj[, PTx:= adj[, YX] / adj[[tmp]]]
	adj[, Sx.e0:= PYs/X.msm.e0]					
	adj[, Sx.e0cp:= PYs/X.msm.e0cp]
	#	prepare risk.df and nt.table as needed
	adj.s		<- copy(adj)
	risk.df		<- merge(risk.df, subset(adj.s, select=c(risk, factor, PYs, PTx, X.msm.e0, X.msm.e0cp)), by=c('risk','factor'))
	nt.table	<- merge(nt.table, subset(adj, select=c(risk, factor, Sx.e0, Sx.e0cp)), by=c('risk','factor'))
	#	add PTx.ref to risk.df
	setkey(risk.df, risk, factor)
	tmp			<- subset(unique(risk.df), select=c(risk, factor, PTx))
	setnames(tmp, colnames(tmp), paste(colnames(tmp),'.ref',sep=''))
	risk.df		<- merge(risk.df, tmp, by=c('risk.ref','factor.ref'))
	set(risk.df, NULL, 'risk.ref', risk.df[, as.character(risk.ref)])
	set(risk.df, NULL, 'factor.ref', risk.df[, as.character(factor.ref)])
	set(risk.df, NULL, 'risk', risk.df[, as.character(risk)])
	set(risk.df, NULL, 'factor', risk.df[, as.character(factor)])	
	#	construct bias and censoring adjustments
	adj	<- subset(adj.s, !grepl('U',factor), c(risk, factor, PYs, X.msm.e0cp))[, list(risk=risk, factor=factor, X.msm.e0cp=sum(X.msm.e0cp)*PYs/sum(PYs))]
	adj	<- rbind(adj, subset(adj.s, grepl('U',factor), c(risk, factor, X.msm.e0cp)))		
	adj	<- merge(subset(adj.s, select=c(risk, factor, PTx, PYs, X.msm.e0)), adj, by=c('risk','factor'))	
	adj	<- adj[, list(	risk=risk, factor=factor, cens.e0= 1, cens.e0cp= X.msm.e0cp/PYs	)]		
	tmp	<- subset(adj.s, select=c(risk, factor, PTx, PYs, X.msm.e0, X.msm.e0cp))
	tmp	<- tmp[, list(	risk=risk, factor=factor,
					bias.e0= X.msm.e0/PYs, bias.e0cp= X.msm.e0cp/PYs, 
					P= PYs/sum(PYs), P.raw= PYs/sum(PYs), 
					P.e0= X.msm.e0/sum(X.msm.e0), P.e0cp= X.msm.e0cp/sum(X.msm.e0cp), 
					P.raw.e0= X.msm.e0/sum(X.msm.e0), P.raw.e0cp= X.msm.e0cp/sum(X.msm.e0cp)	)]
	adj	<- merge(adj, tmp, by=c('risk','factor'))	
	#
	set(nt.table, NULL, 'risk', nt.table[, as.character(risk)])
	set(nt.table, NULL, 'factor', nt.table[, as.character(factor)])
	set(nt.table, NULL, 'Patient', nt.table[, as.character(Patient)])
	#
	#	ready to go:
	#
	#	potential transmission intervals in denominator pop (either X.seq or X.clu)
	cat(paste('\nintervals with evidence for or against transmission (denominator)'))
	setkey(risk.df, coef)
	tmp			<- unique(risk.df)
	tmp[, stat:='PYs']
	tmp			<- subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, PYs))
	set(tmp, NULL, c('coef.ref','risk.ref', 'factor.ref'), 'None')
	setnames(tmp,'PYs','v')
	risk.ans	<- copy(tmp)
	#	number of adjusted potential transmission intervals
	cat(paste('\nnumber of adjusted potential transmission intervals'))
	tmp			<- subset(adj.s, select=c(risk, factor, X.msm.e0cp))
	tmp[, stat:='X.msm.e0cp']
	set(tmp, NULL, c('coef.ref','risk.ref', 'factor.ref'), 'None')
	set(tmp, NULL, 'coef', tmp[, paste(risk,factor,sep='')])
	setnames(tmp,'X.msm.e0cp','v')	
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#	number of expected missing transmission intervals
	cat(paste('\nnumber of expected missing transmission intervals'))
	tmp			<- melt( subset(adj.s, select=c(risk, factor, YX, X.seq, PTx, Sx.e0, Sx.e0cp)), measure.vars=c('Sx.e0','Sx.e0cp'), value.name='Sx', variable.name='stat' )
	set(tmp, NULL, 'v', tmp[, round(YX*(1-Sx)/Sx)])
	set(tmp, NULL, c('risk.ref','factor.ref','coef.ref'), 'None')
	set(tmp, NULL, 'coef', tmp[, paste(risk,factor,sep='')])
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#	number of recipients (including those with no likely transmitters)
	tmp			<- subset(adj.s, select=c(risk, factor, nRec))
	set(tmp, NULL, c('risk.ref','factor.ref','coef.ref'), 'None')
	set(tmp, NULL, 'coef', tmp[, paste(risk,factor,sep='')])
	set(tmp, NULL, 'stat', 'nRec')
	setnames(tmp, 'nRec', 'v')
	risk.ans	<- rbind( risk.ans, tmp, use.names=TRUE )	
	#	number of recipients with likely transmitters
	set(tmp, NULL, 'v', nrow( subset( nt.table[, list(S=all(YX==0)), by='Patient'], !S) ))
	set(tmp, NULL, 'stat', 'nRecLkl')
	risk.ans	<- rbind( risk.ans, tmp, use.names=TRUE )
	#
	#	number of transmissions and proportions by raw count
	#	*** among recipients with likely transmitter only ***
	#
	cat(paste('\nnumber and proportion of transmissions using raw evidence for transmission (y_ijt*w_ijt)'))
	missing		<- merge(nt.table, unique( subset( risk.df, select=c(risk, factor, PTx) ) ), by=c('risk','factor'))
	#	reduce to sampled recipients
	tmp			<- unique(subset(YX.m3, select=Patient))		
	missing		<- merge(tmp, missing, by='Patient', allow.cartesian=TRUE)	
	#	compute the mean Y by stage for non-zero Y
	if(!use.YXf)
	{
		tmp			<- YX.m3[, list(yYX.mean= mean(score.Y), yYX.med= median(score.Y)), by='stage']	
	}
	if(use.YXf)
	{
		tmp			<- YXf[, list(yYX.mean= mean(score.Y), yYX.med= median(score.Y)), by='stage']	
	}
	setnames(tmp, 'stage','factor')
	tmp[, risk:='stage']	
	set(tmp, NULL, 'factor', tmp[, as.character(factor)])
	missing		<- merge(missing, tmp, by=c('risk','factor'))
	#	compute the sum of observed Y's by stage for each recipient
	tmp			<- YX.m3[, list(yYX.sum= sum(score.Y), YX.n=length(score.Y), YX.w=w.t[1]), by=c('stage','Patient')]
	setnames(tmp, 'stage','factor')
	set(tmp, NULL, 'factor', tmp[, as.character(factor)])
	tmp[, risk:='stage']		
	missing		<- merge(missing, tmp, by=c('risk','factor','Patient'), all.x=TRUE)
	set(missing, missing[, which(is.na(YX.n))], 'YX.w', 1.) 
	set(missing, missing[, which(is.na(YX.n))], c('yYX.sum','YX.n'), 0.)
	#	calculate number expected missing potential transmission intervals for recipient j and factor x
	tmp			<- melt( subset(missing, select=c(Patient, risk, factor, X.seq, Sx.e0, Sx.e0cp, PTx)), measure.vars=c('Sx.e0','Sx.e0cp'), variable.name='Sx.method', value.name='Sx' )
	set( tmp, tmp[, which(X.seq==0)], 'X.seq', 1)	#mean of sampling for no observed intervals coincides with sampling model for 1 observed interval 
	#tmp[, YXm.e:= tmp[, as.integer(round( X.seq*(1-Sx)/Sx*PTx )) ]]
	tmp[, YXm.e:= tmp[, X.seq*(1-Sx)/Sx*PTx ]]
	#subset(tmp, factor=='UA.1')[, sum(X.seq+YXm.e)]
	#subset(tmp, factor=='UA.1')[, sum(YXm.e)]
	set(tmp, NULL, 'Sx.method', tmp[, gsub('Sx','YXm.e',Sx.method)])
	missing		<- merge(missing, dcast.data.table(tmp, Patient + risk + factor ~ Sx.method, value.var="YXm.e"), by=c('Patient','risk','factor'))	
	#missing[, YXm.sum.e0:=YXm.e.e0*yYX.mean]
	#missing[, YXm.sum.e0cp:=YXm.e.e0cp*yYX.mean]	
	missing[, YXm.sum.e0:=YXm.e.e0*yYX.med]
	missing[, YXm.sum.e0cp:=YXm.e.e0cp*yYX.med]
	#subset(missing, factor=='UA.1')[, sum(YXm.sum.e0cp)]
	#missing[, list(CHECK=sum(YXm.sum.e0cp)), by=c('risk','factor')]

	#	calculate prob Pj(x) that recipient j got infected from x  - with and without adjustment	
	tmp			<- missing[, 	list(	factor=factor, 	
										Pjx= yYX.sum*YX.w/sum(yYX.sum*YX.w), 
										Pjx.e0= (yYX.sum+YXm.sum.e0)*YX.w/sum((yYX.sum+YXm.sum.e0)*YX.w),					
										Pjx.e0cp= (yYX.sum+YXm.sum.e0cp)*YX.w/sum((yYX.sum+YXm.sum.e0cp)*YX.w),
										coef=paste(risk,as.character(factor), sep='')), by=c('risk','Patient')]
	#subset(tmp, factor=='UA.1')[, hist(Pjx.e0cp, breaks=100)]
	#subset(tmp, factor=='UA.1' & Pjx.e0cp<=0)
	#	compute N.raw etc	
	tmp			<- tmp[, list(	N.raw= sum(Pjx), N.raw.e0= sum(Pjx.e0), N.raw.e0cp=sum(Pjx.e0cp), 
								#n=length(Pjx), PJx.raw.mea= mean(Pjx), PJx.raw= median(Pjx), PJx.raw.e0cp.mea= mean(Pjx.e0cp), PJx.raw.e0cp= median(Pjx.e0cp),  
								risk.ref='None', factor.ref='None', coef.ref='None', coef=coef[1]), by=c('risk','factor')]		
	tmp			<- melt(tmp, id.vars=c('coef','risk','factor','coef.ref','risk.ref','factor.ref'), variable.name='stat', value.name = "v")
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#compute proportions from the various N.raw
	tmp			<- tmp[, list(v= v/sum(v, na.rm=TRUE), risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef=coef, coef.ref='None') , by='stat']
	set(tmp, NULL, 'stat', tmp[, gsub('N','P',stat)])
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))	
	#
	#	relative transmissibilities
	#	*** compare proportions transmitted among recipient with likely transmitter against population proportion among all recipient ***
	#
	cat(paste('\nrelative transmissibilities'))	
	tmp			<- subset( melt(adj, id.vars=c('risk','factor'), variable.name='stat', value.name='p'), !grepl('cens',stat) & !grepl('^bias',stat)	)	
	tmp			<- merge( subset(risk.ans, grepl('P.',stat) | stat=='P' ), tmp, by=c('risk','factor','stat'))
	set(tmp, NULL, 'v', tmp[, v/p])
	set(tmp, tmp[,which(is.nan(v) | v==0)],'v',NA_real_)	
	set(tmp, NULL, 'stat', tmp[, gsub('P','RI',stat)])
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#	risk ratio from relative transmissibilities
	cat(paste('\n term-wise risk ratio from identity RR(x,z)=RT(x)/RT(z).'))	
	setkey(risk.df, coef.ref, coef)
	tmp			<- subset(risk.df, coef!=coef.ref, select=c(coef.ref, coef))
	tmp.ref		<- merge(data.table(coef=tmp[, unique(coef.ref)]), subset(risk.ans, grepl('RI',stat), select=c(coef, risk, factor, stat, v)), by='coef')
	setnames(tmp.ref, c('coef','risk','factor','v'), c('coef.ref','risk.ref','factor.ref','v.ref'))
	tmp			<- merge(tmp, subset(risk.ans, grepl('RI',stat), select=c(coef, risk, factor, stat, v)), by='coef', allow.cartesian=TRUE)
	tmp			<- merge(tmp.ref, tmp, by=c('stat','coef.ref'))
	tmp.ref		<- NULL
	set(tmp, NULL, 'v', tmp[, v/v.ref])
	set(tmp, NULL, 'stat', tmp[, gsub('RI','RR',stat)])
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#	Bootstraps	on ratio and on prob
	setkey(ct, t.period, risk, factor)
	tmp			<- ct[, seq_len( length(n.adj)/length(unique(factor)) )]
	ct[, bs:=rep(tmp, nrow(ct)/length(tmp))]
	cat(paste('\nregression on bootstrap data sets bs.n=',bs.n))
	tmp			<- lapply(seq_len(bs.n), function(bs.i)
			{
				if(bs.i%%100==0)	cat(paste('\nprocess bootstrap data sets bs.i=',bs.i))
				#
				#	bootstrap over recently infected Patient
				#
				tmp					<- unique(subset(YX.m3, select=Patient))					
				YX.m3.bs			<- tmp[ sample( seq_len(nrow(tmp)), nrow(tmp), replace=TRUE ), ]
				#	debug
				#YX.m3.bs			<- unique(subset(YX.m3, select=Patient))				
				#	recipient MSM are not necessarily unique any longer - need to create unique bs id
				YX.m3.bs[, Patient.bs:=paste(Patient, seq_len(nrow(tmp)),sep='_bs' )]
				YX.m3.bs			<- merge( YX.m3, YX.m3.bs, by='Patient', allow.cartesian=TRUE )
				#
				#	bootstrap denominator population: y_ijt=0 AND y_ijt>0
				#
				tmp					<- subset(risk.ans, stat=='PYs')
				set(tmp, NULL, 'v', rbinom(nrow(tmp), tmp[, sum(v)], tmp[, v/sum(v)]))
				risk.ans.bs			<- subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v))
				#	bootstrap recipient with likely transmitter				
				set(tmp, NULL, 'v', YX.m3.bs[, length(unique(Patient.bs))])
				set(tmp, NULL, 'stat', 'nRecLkl')
				risk.ans.bs			<- rbind( risk.ans.bs, tmp, use.names=TRUE )			
				#	bootstrap non-zero scores
				tmp			<- risk.df[,	{
												z	<- table( YX.m3.bs[, risk, with=FALSE])
												list(factor=rownames(z), n=as.numeric(unclass(z)))												
											}, by='risk']
				tmp			<- merge( subset(risk.ans.bs, stat=='PYs'), tmp, by=c('risk','factor'))		
				set(tmp, NULL, 'v', tmp[, n/v])
				set(tmp, NULL, 'stat', 'PTx')				
				risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
				#	bootstrap number of likely transmission intervals
				tmp			<- YX.m3.bs[, list(v=length(t), stat='YX'), by='stage']
				set(tmp, NULL, 'risk', 'stage')
				setnames(tmp, 'stage', 'factor')
				set(tmp, NULL, c('risk.ref','factor.ref','coef.ref'), 'None')
				set(tmp, NULL, 'coef', tmp[, paste(risk,factor,sep='')])				
				risk.ans.bs	<- rbind(risk.ans.bs, tmp, use.names=TRUE)
				#
				#	bootstrap sample censoring
				#
				tmp					<- subset( ct, bs==sample(ct[, unique(bs)], 1) )
				set(tmp, NULL, 'n.adj', tmp[, as.integer(round(n.adj))])	
				setnames(tmp, 'n.adj', 'n.adj.bs' )
				tmp					<- merge( nt.table[, list(X.msm.e0=sum(X.msm.e0)), by=c('risk','factor')], tmp, by=c('risk','factor'))
				stopifnot( tmp[, all(X.msm.e0<=n.adj.bs)] )			 
				tmp					<- tmp[, list( PYe0cpr=round(n.adj.bs-X.msm.e0)), by=c('risk','factor')]		
				tmp					<- merge(tmp, nt.table[, list(X.msm.e0cp= length(unique(Patient))), by=c('risk','factor') ], by=c('risk','factor')) 
				set(tmp, NULL, 'X.msm.e0cp', round( tmp[, PYe0cpr/X.msm.e0cp] ))				
				nt.table.bs			<- merge(subset(nt.table, select=c(risk, factor, Patient, X.clu, X.msm.e0, X.seq, YX)), subset(tmp, select=c(risk, factor, X.msm.e0cp)), by=c('risk','factor'))
				set(nt.table.bs, NULL, 'X.msm.e0cp', nt.table.bs[, X.msm.e0+X.msm.e0cp])	
				#	combine bootstrap censoring with boostrap denom pop / non-zero
				tmp					<- nt.table.bs[, list(YX=sum(YX), X.clu=sum(X.clu), X.seq=sum(X.seq), X.msm.e0=sum(X.msm.e0), X.msm.e0cp=sum(X.msm.e0cp), nRec=length(unique(Patient))), by=c('risk','factor')]
				tmp					<- merge(tmp, subset(risk.ans.bs, stat=='PYs', c(risk, factor, v)), by=c('risk','factor'))
				setnames(tmp, 'v', 'PYs')
				tmp[, Sx.e0:= PYs/X.msm.e0]					
				tmp[, Sx.e0cp:= PYs/X.msm.e0cp]
				tmp					<- merge(tmp, subset(risk.ans.bs, stat=='PTx', c(risk, factor, v)), by=c('risk','factor'))
				setnames(tmp, 'v', 'PTx')
				set(tmp, NULL, 'risk', tmp[, as.character(risk)])
				set(tmp, NULL, 'factor', tmp[, as.character(factor)])								
				#	prepare risk.df and nt.table.bs as needed
				nt.table.bs			<- merge(nt.table.bs, subset(tmp, select=c(risk, factor, PTx, Sx.e0, Sx.e0cp)), by=c('risk','factor'))				
				#	bootstrap number of recipients				
				set(tmp, NULL, c('risk.ref','factor.ref','coef.ref'), 'None')
				set(tmp, NULL, 'coef', tmp[, paste(risk,factor,sep='')])
				set(tmp, NULL, 'stat', 'nRec')
				setnames(tmp, 'nRec', 'v')
				risk.ans.bs	<- rbind( risk.ans.bs, subset(tmp, select=c(risk, factor, risk.ref, factor.ref, coef.ref, coef, stat, v)), use.names=TRUE )
				#	bootstrap number of adjusted potential transmission intervals
				set(tmp, NULL, 'stat', 'X.msm.e0cp')
				setnames(tmp, c('v','X.msm.e0cp'),  c('nRec','v'))
				risk.ans.bs	<- rbind( risk.ans.bs, subset(tmp, select=c(risk, factor, risk.ref, factor.ref, coef.ref, coef, stat, v)), use.names=TRUE )
				#
				#	raw number of transmissions
				#
				#	reduce to bootstrap sampled recipients
				tmp			<- subset(YX.m3.bs, select=c(Patient, Patient.bs))		
				setkey(tmp, Patient, Patient.bs)
				missing		<- merge(unique(tmp), nt.table.bs, by='Patient', allow.cartesian=TRUE)
				#	compute the sum of observed Y's by risk factor for each recipient
				tmp			<- YX.m3.bs[, list(Patient=Patient[1], yYX.sum= sum(score.Y), YX.bs=length(score.Y), YX.w=w.t[1]), by=c('stage','Patient.bs')]
				setnames(tmp, 'stage','factor')
				set(tmp, NULL, 'factor', tmp[, as.character(factor)])	
				tmp[, risk:='stage']	
				missing		<- merge(missing, tmp, by=c('risk','factor','Patient','Patient.bs'), all.x=TRUE)
				set(missing, missing[, which(is.na(YX.bs))], 'YX.w', 1.) 
				set(missing, missing[, which(is.na(YX.bs))], c('yYX.sum','YX.bs'), 0.)
				#	to avoid rounding issues, sample total and then re-allocate to individual recipient MSM
				tmp			<- missing[, {
												#	total missed through sampling
												YXm.r.e0	<- rznbinom( 1, sum(X.seq), Sx.e0[1] )
												YXm.r.e0cp	<- rznbinom( 1, sum(X.seq), Sx.e0cp[1] )
												#YXm.r.e0	<- as.integer(round(  (1-Sx.e0[1])*sum(X.seq)/Sx.e0[1]  		))
												#YXm.r.e0cp	<- as.integer(round(  (1-Sx.e0cp[1])*sum(X.seq)/Sx.e0cp[1]  	))
												#	total missed with nonzero score
												YXm.r.e0	<- rbinom(1, YXm.r.e0, PTx[1])
												YXm.r.e0cp	<- rbinom(1, YXm.r.e0cp, PTx[1])
												#YXm.r.e0	<- as.integer(round(YXm.r.e0*PTx[1]))
												#YXm.r.e0cp	<- as.integer(round(YXm.r.e0cp*PTx[1]))
												#	distribute among all recipients as fraction
												#YXm.r.e0	<- as.vector(rmultinom(1, YXm.r.e0, rep( 1/length(X.seq), length(X.seq))))
												#YXm.r.e0cp	<- as.vector(rmultinom(1, YXm.r.e0cp, rep( 1/length(X.seq), length(X.seq))))
												YXm.r.e0	<- YXm.r.e0 / length(X.seq)
												YXm.r.e0cp	<- YXm.r.e0cp / length(X.seq)
												list(Patient.bs=Patient.bs, YXm.r.e0=YXm.r.e0, YXm.r.e0cp=YXm.r.e0cp)
										}, by=c('risk','factor')]
				#subset(tmp, factor=='UA.1')[, sum(YXm.r.e0cp)]
				missing		<- merge(missing, tmp, by=c('risk','factor','Patient.bs'))				
				#	draw missing scores from all yijt in that stage	for number missing YXm.r.e0
				tmp			<- missing[, {			
												if(!use.YXf)
												{
													z	<- median( YX.m3[ which( YX.m3[[risk]]==factor ), ][['score.Y']] )	
												}
												if(use.YXf)
												{
													z	<- median( YXf[ which( YXf[[risk]]==factor ), ][['score.Y']] )	
												}												
												#	TODO: if activate, double check sample -- these were integers !?
												#list(Patient.bs=rep(Patient.bs, YXm.r.e0), yYXm.r.e0=sample(z, sum(YXm.r.e0), replace=TRUE)  )
												list(Patient.bs=Patient.bs, yYXm.sum.e0=YXm.r.e0*z, yYXm.sum.e0cp=YXm.r.e0cp*z )
										}, by=c('risk','factor')]									
				#tmp[, list(CHECK=sum(yYXm.sum.e0cp)), by=c('risk','factor')]								
				missing		<- merge(missing, tmp, by=c('Patient.bs','risk','factor'), all.x=TRUE)
				#	calculate prob Pj(x) that recipient j got infected from x  - with and without adjustment				
				tmp			<- missing[, 	list(	factor=factor, 	
													Pjx= yYX.sum*YX.w/sum(yYX.sum*YX.w), 
													Pjx.e0= (yYX.sum+yYXm.sum.e0)*YX.w/sum((yYX.sum+yYXm.sum.e0)*YX.w),
													Pjx.e0cp= (yYX.sum+yYXm.sum.e0cp)*YX.w/sum((yYX.sum+yYXm.sum.e0cp)*YX.w),
													coef=paste(risk,as.character(factor), sep='')), by=c('risk','Patient.bs')]
				#subset(tmp, factor=='UA.1')[, hist(Pjx.e0cp, breaks=100)]					
				#	exclude recipients with no evidence for direct transmission
				tmp			<- subset(tmp, !is.nan(Pjx))					
				#various N.raw									
				tmp			<- tmp[, list(	N.raw= sum(Pjx), N.raw.e0= sum(Pjx.e0), N.raw.e0cp=sum(Pjx.e0cp), 
								risk.ref='None', factor.ref='None', coef.ref='None', coef=coef[1]), by=c('risk','factor')]		
				tmp			<- melt(tmp, id.vars=c('coef','risk','factor','coef.ref','risk.ref','factor.ref'), variable.name='stat', value.name = "v")
				risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))				
				#
				risk.ans.bs[, bs:=bs.i]
			})
	risk.ans.bs	<- do.call('rbind',tmp)
	#
	#	compute adjustments for each bootstrap sample: collect PYe0 etc to compute sampling proportion relative to bootstrap PYs
	#
	tmp			<- dcast.data.table( subset(risk.ans.bs, stat=='PYs' | stat=='PTx', select=c(risk,factor,v,bs,stat)), risk+factor+bs ~ stat, value.var='v' )
	tmp			<- merge( tmp, unique(subset(risk.df, select=c(risk, factor, X.msm.e0, X.msm.e0cp))), by=c('risk','factor'))
	tmp2		<- tmp[, list(  	factor=factor, 
					X.msm.e0=rbinom(length(X.msm.e0), sum(X.msm.e0), X.msm.e0/sum(X.msm.e0)),
					X.msm.e0cp=rbinom(length(X.msm.e0cp), sum(X.msm.e0cp), X.msm.e0cp/sum(X.msm.e0cp))	), by=c('bs','risk')]
	tmp			<- merge( subset(tmp, select=c(risk, factor, bs, PYs, PTx)), tmp2, by=c('bs','risk','factor') )	
	adj.bs		<- tmp[, list(	factor=factor,
								P= PYs/sum(PYs), P.raw= PYs/sum(PYs), 
								P.e0= X.msm.e0/sum(X.msm.e0), P.e0cp= X.msm.e0cp/sum(X.msm.e0cp),
								P.raw.e0= X.msm.e0/sum(X.msm.e0), P.raw.e0cp= X.msm.e0cp/sum(X.msm.e0cp),
								Sx.e0= PYs/X.msm.e0, Sx.e0cp= PYs/X.msm.e0cp), by=c('bs','risk')]
	#
	# 	add Proportions
	#
	tmp			<- melt( subset(adj.bs, select=c(bs, risk, factor, Sx.e0, Sx.e0cp)), measure.vars=c('Sx.e0','Sx.e0cp'), value.name='s', variable.name='stat' )
	tmp			<- merge(tmp, subset(risk.ans.bs, stat=='YX', select=c(bs, risk, factor, v)), by=c('bs','risk','factor'))
	set(tmp, NULL, 'v', tmp[, round(v*(1-s)/s)])
	set(tmp, NULL, c('risk.ref','factor.ref','coef.ref'), 'None')
	set(tmp, NULL, 'coef', tmp[, paste(risk,factor,sep='')])
	set(tmp, NULL, 's', NULL)
	risk.ans.bs	<- rbind( risk.ans.bs, tmp, use.names=TRUE )
	#
	# 	add Proportions
	#
	tmp			<- subset(risk.ans.bs, grepl('N',stat))	
	tmp			<- tmp[, list(v= v/sum(v, na.rm=TRUE), risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef=coef, coef.ref='None') , by=c('stat','bs')]	
	set(tmp, NULL, 'stat', tmp[, gsub('N','P',stat)])
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v, bs)))
	#	set missing values in N P to zero
	set(risk.ans.bs, risk.ans.bs[, which(is.na(v) & (grepl('P',stat) | grepl('PYs',stat) | grepl('N',stat)))], 'v', 0.)
	#
	#	add RIs
	#	
	tmp			<- melt(adj.bs, id.vars=c('risk','factor','bs'), variable.name='stat', value.name='p')
	tmp2		<- subset(risk.ans.bs, grepl('P',stat) & !grepl('PT',stat) & !grepl('PY',stat) )
	tmp			<- merge(tmp2, tmp, by=c('risk','factor','stat','bs'))
	set(tmp, NULL, 'v', tmp[, v/p])
	set(tmp, tmp[,which(is.nan(v) | v==0)],'v',NA_real_)	
	set(tmp, NULL, 'stat', tmp[, gsub('P','RI',stat)])
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v, bs)))
	#
	#	Risk ratio from relative transmissibilities
	#
	tmp			<- subset(risk.df, coef!=coef.ref, select=c(coef.ref, coef))
	tmp2		<- merge(data.table(coef=tmp[, unique(coef.ref)]), subset(risk.ans.bs, grepl('RI',stat), select=c(bs, coef, risk, factor, stat, v)), by='coef')
	setnames(tmp2, c('coef','risk','factor','v'), c('coef.ref','risk.ref','factor.ref','v.ref'))
	tmp			<- merge(tmp, subset(risk.ans.bs, grepl('RI',stat), select=c(bs, coef, risk, factor, stat, v)), by='coef', allow.cartesian=TRUE)
	tmp			<- merge(tmp2, tmp, by=c('bs','stat','coef.ref'))
	tmp2		<- NULL
	set(tmp, NULL, 'v', tmp[, v/v.ref])
	set(tmp, NULL, 'stat', tmp[, gsub('RI','RR',stat)])
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v, bs)))
	#
	tmp			<- risk.ans.bs[,	list(l95.bs=quantile(v, prob=0.025, na.rm=TRUE), u95.bs=quantile(v, prob=0.975, na.rm=TRUE), m50.bs=quantile(v, prob=0.5, na.rm=TRUE)), by=c('coef','coef.ref','stat')]
	risk.ans	<- merge(risk.ans, tmp, by=c('coef','coef.ref','stat'), all.x=TRUE)
	#	
	setkey(risk.ans, stat, coef.ref, coef)
	list(risk=risk.ans, risk.bs=risk.ans.bs)
}
######################################################################################
project.athena.Fisheretal.estimate.risk.core<- function(YX.m3, X.seq, formula, predict.df, risk.df, include.colnames, bs.n=1e3, gamlss.BE.required.limit=0.99, gamlss.BE.limit.u=c(0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0, 0) )
{
	require(gamlss)
	options(warn=0)
	#
	cat(paste('\nbeta regression for formula=',formula))	 
	tmp				<- project.athena.Fisheretal.betareg(YX.m3, formula, include.colnames, gamlss.BE.limit.u=gamlss.BE.limit.u, gamlss.BE.limit.l=gamlss.BE.limit.l, verbose=1)
	betafit.rr		<- tmp$betafit.rr
	betafit.or		<- tmp$betafit.or	
	stopifnot(tmp$gamlss.BE.limit>gamlss.BE.required.limit)
		
	#AIC(YX.m3.fit1, YX.m3.fit2, YX.m3.fit3, YX.m3.fit2b, YX.m3.fit4, YX.m3.fit5, YX.m3.fit6, YX.m3.fit7 )
	#AIC(YX.m3.fit1, YX.m3.fit2, YX.m3.fit3, YX.m3.fit2b, YX.m3.fit4, YX.m3.fit5, YX.m3.fit6, YX.m3.fit7, k=YX.m3[, log(round(sum(w)))])
	#sapply( list(YX.m3.fit1, YX.m3.fit2, YX.m3.fit3, YX.m3.fit2b, YX.m3.fit4, YX.m3.fit5, YX.m3.fit6, YX.m3.fit7), '[[', 'pseudo.r.squared')
	#str(X.seq)
	#str(YX.m3)	
	#	add PTx.ref to risk.df
	setkey(risk.df, risk, factor)
	tmp			<- subset(unique(risk.df), select=c(risk, factor, PTx))
	setnames(tmp, colnames(tmp), paste(colnames(tmp),'.ref',sep=''))
	risk.df		<- merge(risk.df, tmp, by=c('risk.ref','factor.ref'))	
	#	odds ratio and risk ratio
	cat(paste('\nodds ratios and risk ratios'))
	setkey(risk.df, coef.ref, coef)
	risk.ans	<- subset(risk.df, coef!=coef.ref)[, 	{										
				tmp				<- copy(predict.df)
				tmp2			<- copy(predict.df)
				set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
				set(tmp2, NULL, risk, factor(factor.ref, levels=levels( predict.df[[risk]] )))
				corisk			<- intersect(colnames(risk.df), colnames(predict.df))
				if(!length(corisk))	
					corisk<- NULL				
				tmp				<- na.omit(merge(tmp[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3[which( YX.m3[, risk, with=FALSE][[1]]==factor ), c(risk,corisk,'w'), with=FALSE], by=risk))
				setkey(tmp, NULL)
				tmp2			<- na.omit(merge(tmp2[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3[which( YX.m3[, risk.ref, with=FALSE][[1]]==factor.ref ), c(risk.ref,corisk,'w'), with=FALSE], by=risk.ref))
				setkey(tmp2, NULL)
				if(nrow(tmp))	#	with corisks E(exp(beta*X)) != exp(beta*E(X)) so it s all a bit more complicated:
				{
					if( nrow(unique(tmp[, setdiff(colnames(tmp),'w'), with=FALSE]))==1 )
						tmp		<- tmp[1,]
					if( nrow(unique(tmp2[, setdiff(colnames(tmp2),'w'), with=FALSE]))==1 )
						tmp2	<- tmp2[1,]
					tryCatch({
					tmp[, predict:= predict(betafit.or, newdata=as.data.frame(tmp), data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), type='link')]
					tmp2[, predict:= predict(betafit.or, newdata=as.data.frame(tmp2), data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), type='link')]
					tmp			<- weighted.mean(exp(tmp[,predict]), w=tmp[,w], na.rm=TRUE)/weighted.mean(exp(tmp2[,predict]), w=tmp2[,w], na.rm=TRUE)			
					}, error=function(e){ print(e$message); tmp<<- NA_real_ }, warning=function(w){ print(w$message); tmp<<- NA_real_ })
				}
				else
					tmp			<- NA_real_		
				list(stat= 'OR', risk=risk, factor=factor, risk.ref=risk.ref, factor.ref=factor.ref, v=tmp )
			},by=c('coef','coef.ref')]
	tmp			<- subset(risk.df, coef!=coef.ref)[, 	{
				tmp				<- copy(predict.df)
				tmp2			<- copy(predict.df)
				set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
				set(tmp2, NULL, risk, factor(factor.ref, levels=levels( predict.df[[risk]] )))
				corisk			<- intersect(colnames(risk.df), colnames(predict.df))
				if(!length(corisk))	
					corisk<- NULL				
				tmp				<- na.omit(merge(tmp[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3[which( YX.m3[, risk, with=FALSE][[1]]==factor ), c(risk,corisk,'w'), with=FALSE], by=risk))
				setkey(tmp, NULL)
				tmp2			<- na.omit(merge(tmp2[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3[which( YX.m3[, risk.ref, with=FALSE][[1]]==factor.ref ), c(risk.ref,corisk,'w'), with=FALSE], by=risk.ref))
				setkey(tmp2, NULL)
				#	with corisks E(exp(beta*X)) != exp(beta*E(X)) so it s all a bit more complicated:
				if(nrow(tmp))	
				{
					if( nrow(unique(tmp[, setdiff(colnames(tmp),'w'), with=FALSE]))==1 )
						tmp		<- tmp[1,]
					if( nrow(unique(tmp2[, setdiff(colnames(tmp2),'w'), with=FALSE]))==1 )
						tmp2	<- tmp2[1,]									
					tryCatch({
					tmp[, predict:= predict(betafit.rr, newdata=as.data.frame(tmp), data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), type='link')]
					tmp2[, predict:= predict(betafit.rr, newdata=as.data.frame(tmp2), data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), type='link')]
					tmp			<- weighted.mean(exp(tmp[,predict]), w=tmp[,w], na.rm=TRUE)/weighted.mean(exp(tmp2[,predict]), w=tmp2[,w], na.rm=TRUE)
					}, error=function(e){ print(e$message); tmp<<- NA_real_ }, warning=function(w){ print(w$message); tmp<<- NA_real_ })
				}
				else
					tmp			<- NA_real_					
				list(stat= 'RR', risk=risk, factor=factor, risk.ref=risk.ref, factor.ref=factor.ref, v=tmp )										
			},by=c('coef','coef.ref')]
	risk.ans	<- rbind(risk.ans, tmp)
	#	risk ratio for PDT assuming no transmissions from Not.PT
	tmp			<- merge( subset(risk.ans, stat=='RR'), subset(risk.df, select=c(risk, factor, risk.ref, factor.ref, PTx, PTx.ref)), by=c('risk','risk.ref','factor','factor.ref'))
	set(tmp, NULL, 'v', tmp[, v*PTx/PTx.ref])
	set(tmp, NULL, 'stat', 'RR.ptx')
	risk.ans	<- rbind( risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)) )	
	#	term-wise risk ratio
	cat(paste('\nterm wise risk ratios'))	
	tmp			<- subset(risk.df, coef!=coef.ref)[, 	{
				tmp				<- copy(predict.df)
				tmp2			<- copy(predict.df)
				set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
				set(tmp2, NULL, risk, factor(factor.ref, levels=levels( predict.df[[risk]] )))
				#	with corisks E(exp(beta*X)) != exp(beta*E(X)) so it s all a bit more complicated:
				tryCatch({
								tmp		<- predict(betafit.rr, newdata=as.data.frame(tmp), data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), terms=risk, type='terms')
								tmp2	<- predict(betafit.rr, newdata=as.data.frame(tmp2), data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), terms=risk, type='terms')								
								tmp		<- exp( tmp[,risk] - tmp2[,risk] )
							}, error=function(e){ print(e$message); tmp<<- NA_real_ }, warning=function(w){ print(w$message); tmp<<- NA_real_ })
				list(stat= 'RR.term', risk=risk, factor=factor, risk.ref=risk.ref, factor.ref=factor.ref, v=tmp )										
			},by=c('coef','coef.ref')]
	risk.ans	<- rbind(risk.ans, tmp)	
	#	term-wise risk ratio for PDT assuming no transmissions from Not.PT	(this is not really valid, anyhow)
	tmp			<- merge( subset(risk.ans, stat=='RR.term'), subset(risk.df, select=c(risk, factor, risk.ref, factor.ref, PTx, PTx.ref)), by=c('risk','risk.ref','factor','factor.ref'))
	set(tmp, NULL, 'v', tmp[, v*PTx/PTx.ref])
	set(tmp, NULL, 'stat', 'RR.term.ptx')
	risk.ans	<- rbind( risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)) )		
	#	person years in infection window
	cat(paste('\nperson years across infection windows'))
	setkey(risk.df, coef)
	if(is.null(X.seq))
	{
		tmp			<- unique(risk.df)
		tmp[, stat:='PYs']
		tmp			<- subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, PYs))
		set(tmp, NULL, c('coef.ref','risk.ref', 'factor.ref'), 'None')
		setnames(tmp,'PYs','v')
		risk.ans	<- rbind(risk.ans, tmp)		
	}
	if(!is.null(X.seq))
	{		
		tmp			<- unique(risk.df)[, list(risk=risk, factor=factor, risk.ref="None", factor.ref="None", coef.ref="None", stat="PYs", v=length(which(as.character(X.seq[, risk, with=FALSE][[1]])==factor))), by='coef']
		risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	}
	#	number of transmissions and proportion of transmissions
	cat(paste('\nnumber and proportion of transmissions'))
	tmp			<- unique(risk.df)[, 		{
				tmp					<- copy(predict.df)
				set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
				corisk				<- intersect(colnames(risk.df), colnames(predict.df))
				if(!length(corisk))	
					corisk<- NULL				
				tmp					<- na.omit(merge(tmp[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3[which( YX.m3[, risk, with=FALSE][[1]]==factor ), c(risk,corisk,'w'), with=FALSE], by=risk))
				setkey(tmp, NULL)
				if(nrow(tmp))	
				{
					
					if( nrow(unique(tmp[, setdiff(colnames(tmp),'w'), with=FALSE]))==1 )
						tmp			<- tmp[1,]		
					tryCatch({
					tmp[, predict:= predict(betafit.rr, newdata=as.data.frame(tmp), data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), type='link')]				
					tmp				<- weighted.mean(exp(tmp[,predict]), w=tmp[,w], na.rm=TRUE)
					}, error=function(e){ print(e$message); tmp<<- NA_real_ }, warning=function(w){ print(w$message); tmp<<- NA_real_ })
				}
				else
					tmp				<- NA_real_	
				tmp[2]				<- sum( YX.m3[ which(as.character(YX.m3[, risk, with=FALSE][[1]])==factor), w] )
				list(risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef.ref='None', stat= 'N', expbeta=ifelse(tmp[2]<2*EPS, 0., tmp[1]), n=tmp[2] )
			},by= 'coef']
	tmp[, v:= tmp[,expbeta*n]]	
	set(tmp, tmp[, which(v==0)],'v', NA_real_)
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	set(tmp, NULL, 'stat', 'P')
	set(tmp, NULL, 'v', tmp[, v/sum(expbeta*n)])
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#	proportions of transmissions assuming no transmissions from Not.PT
	setkey(risk.df, risk, factor)
	tmp			<- merge(tmp, subset(unique(risk.df), select=c(risk, factor, PTx)), by=c('risk','factor'))
	set(tmp, NULL, 'v', tmp[, expbeta*n*PTx/sum(expbeta*n*PTx)])
	set(tmp, NULL, 'stat', 'P.ptx')
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#	number of transmissions and proportions by raw count
	tmp			<- unique(risk.df)[, 	{
							tmp	<- sum( YX.m3[ which(as.character(YX.m3[, risk, with=FALSE][[1]])==factor), score.Y*w] )
							list(risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef.ref='None', stat= 'N.raw', v=tmp )
						}, by='coef']
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	tmp			<- tmp[, {
							list(risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef.ref='None', stat= 'P.raw', v=v/sum(tmp[,v]) )
						}, by='coef']
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#	relative infectiousness and 
	#	relative infectiousness assuming no transmissions from Not.PT and
	#	raw relative infectiousness
	cat(paste('\nrelative infectiousness'))
	tmp			<- subset(risk.ans, stat=='PYs' )[, sum(v)]
	tmp			<- subset(risk.ans, risk.ref=='None')[, list( 	stat=c('RI','RI.ptx','RI.raw'), risk=risk[1], factor=factor[1], risk.ref='None', factor.ref='None', coef.ref='None', 
																v= c( v[stat=='P'], v[stat=='P.ptx'], v[stat=='P.raw']) / ( v[stat=='PYs'] / tmp ) ), by='coef']
	set(tmp, tmp[,which(is.nan(v) | v==0)],'v',NA_real_)
	risk.ans	<- rbind(risk.ans, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
	#	Bootstraps	on ratio and on prob
	cat(paste('\nregression on bootstrap data sets bs.n=',bs.n))
	tmp			<- lapply(seq_len(bs.n), function(bs.i)
			{
				if(bs.i%%100==0)	cat(paste('\nregression on bootstrap data sets bs.i=',bs.i))
				#bootstrap over recently infected Patient
				tmp					<- unique(subset(YX.m3, select=Patient))					
				#
				bs.repeat			<- 1
				while(bs.repeat)
				{
					YX.m3.bs			<- merge( YX.m3, tmp[ sample( seq_len(nrow(tmp)), nrow(tmp), replace=TRUE ), ], by='Patient', allow.cartesian=TRUE )
					tryCatch({ 
						tmp2			<- project.athena.Fisheretal.betareg(YX.m3.bs, formula, include.colnames, gamlss.BE.limit.u=gamlss.BE.limit.u, gamlss.BE.limit.l=gamlss.BE.limit.l, verbose=1 )
						betafit.or.bs	<- tmp2$betafit.or
						betafit.rr.bs	<- tmp2$betafit.rr
						bs.repeat		<- ifelse( tmp2$gamlss.BE.limit>gamlss.BE.required.limit, 0, 1 )
					}, error=function(e){ bs.repeat<<- 1})																
				}					
				#	odds ratio and risk ratio
				setkey(risk.df, coef.ref, coef)
				risk.ans.bs			<- subset(risk.df, coef!=coef.ref)[, 	{
																				tmp				<- copy(predict.df)
																				tmp2			<- copy(predict.df)
																				set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
																				set(tmp2, NULL, risk, factor(factor.ref, levels=levels( predict.df[[risk]] )))
																				corisk			<- intersect(colnames(risk.df), colnames(predict.df))
																				if(!length(corisk))	
																					corisk<- NULL
																				tmp				<- na.omit(merge(tmp[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3.bs[which( YX.m3.bs[, risk, with=FALSE][[1]]==factor ), c(risk,corisk,'w'), with=FALSE], by=risk))
																				setkey(tmp, NULL)
																				tmp2			<- na.omit(merge(tmp2[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3.bs[which( YX.m3.bs[, risk.ref, with=FALSE][[1]]==factor.ref ), c(risk.ref,corisk,'w'), with=FALSE], by=risk.ref))
																				setkey(tmp2, NULL)
																				#	with corisks E(exp(beta*X)) != exp(beta*E(X)) so it s all a bit more complicated:
																				if(nrow(tmp))	
																				{
																					if( nrow(unique(tmp[, setdiff(colnames(tmp),'w'), with=FALSE]))==1 )
																						tmp		<- tmp[1,]
																					if( nrow(unique(tmp2[, setdiff(colnames(tmp2),'w'), with=FALSE]))==1 )
																						tmp2	<- tmp2[1,]															
																					tryCatch({
																					tmp[, predict:= predict(betafit.or.bs, newdata=as.data.frame(tmp), data=na.omit(as.data.frame(YX.m3.bs[, include.colnames, with=FALSE])), type='link')]
																					tmp2[, predict:= predict(betafit.or.bs, newdata=as.data.frame(tmp2), data=na.omit(as.data.frame(YX.m3.bs[, include.colnames, with=FALSE])), type='link')]
																					tmp			<- weighted.mean(exp(tmp[,predict]), w=tmp[,w], na.rm=TRUE)/weighted.mean(exp(tmp2[,predict]), w=tmp2[,w], na.rm=TRUE)
																					}, error=function(e){ print(e$message); tmp<<- NA_real_ }, warning=function(w){ print(w$message); tmp<<- NA_real_ })
																				}
																				else
																					tmp			<- NA_real_	
																				list(stat= 'OR', risk=risk, factor=factor, risk.ref=risk.ref, factor.ref=factor.ref, v=tmp )
						},by=c('coef','coef.ref')]
				tmp					<- subset(risk.df, coef!=coef.ref)[, 	{
																				tmp				<- copy(predict.df)
																				tmp2			<- copy(predict.df)
																				set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
																				set(tmp2, NULL, risk, factor(factor.ref, levels=levels( predict.df[[risk]] )))
																				corisk			<- intersect(colnames(risk.df), colnames(predict.df))
																				if(!length(corisk))	
																					corisk<- NULL
																				tmp				<- na.omit(merge(tmp[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3.bs[which( YX.m3.bs[, risk, with=FALSE][[1]]==factor ), c(risk,corisk,'w'), with=FALSE], by=risk))
																				setkey(tmp, NULL)
																				tmp2			<- na.omit(merge(tmp2[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3.bs[which( YX.m3.bs[, risk.ref, with=FALSE][[1]]==factor.ref ), c(risk.ref,corisk,'w'), with=FALSE], by=risk.ref))
																				setkey(tmp2, NULL)
																				if(nrow(tmp))	
																				{
																					if( nrow(unique(tmp[, setdiff(colnames(tmp),'w'), with=FALSE]))==1 )
																						tmp		<- tmp[1,]
																					if( nrow(unique(tmp2[, setdiff(colnames(tmp2),'w'), with=FALSE]))==1 )
																						tmp2	<- tmp2[1,]															
																					tryCatch({
																					tmp[, predict:= predict(betafit.rr.bs, newdata=as.data.frame(tmp), data=na.omit(as.data.frame(YX.m3.bs[, include.colnames, with=FALSE])), type='link')]
																					tmp2[, predict:= predict(betafit.rr.bs, newdata=as.data.frame(tmp2), data=na.omit(as.data.frame(YX.m3.bs[, include.colnames, with=FALSE])), type='link')]
																					tmp			<- weighted.mean(exp(tmp[,predict]), w=tmp[,w], na.rm=TRUE)/weighted.mean(exp(tmp2[,predict]), w=tmp2[,w], na.rm=TRUE)
																					}, error=function(e){ print(e$message); tmp<<- NA_real_ }, warning=function(w){ print(w$message); tmp<<- NA_real_ })
																				}
																				else
																					tmp			<- NA_real_
																				list(stat= 'RR', risk=risk, factor=factor, risk.ref=risk.ref, factor.ref=factor.ref, v=tmp )
																			},by=c('coef','coef.ref')]
				risk.ans.bs			<- rbind(risk.ans.bs, tmp)		
				#term wise risk ratio
				tmp					<- subset(risk.df, coef!=coef.ref)[, 	{
																				tmp				<- copy(predict.df)
																				tmp2			<- copy(predict.df)
																				set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
																				set(tmp2, NULL, risk, factor(factor.ref, levels=levels( predict.df[[risk]] )))
																				#	with corisks E(exp(beta*X)) != exp(beta*E(X)) so it s all a bit more complicated:
																				tryCatch({
																							tmp		<- predict(betafit.rr.bs, newdata=as.data.frame(tmp), data=na.omit(as.data.frame(YX.m3.bs[, include.colnames, with=FALSE])), terms=risk, type='terms')
																							tmp2	<- predict(betafit.rr.bs, newdata=as.data.frame(tmp2), data=na.omit(as.data.frame(YX.m3.bs[, include.colnames, with=FALSE])), terms=risk, type='terms')								
																							tmp		<- exp( tmp[,risk] - tmp2[,risk] )
																						}, error=function(e){ print(e$message); tmp<<- NA_real_ }, warning=function(w){ print(w$message); tmp<<- NA_real_ })
																				list(stat= 'RR.term', risk=risk, factor=factor, risk.ref=risk.ref, factor.ref=factor.ref, v=tmp )										
																			},by=c('coef','coef.ref')]
				risk.ans.bs			<- rbind(risk.ans.bs, tmp)
				#	raw number of transmissions
				setkey(risk.df, coef)
				tmp			<- unique(risk.df)[, 	{
														tmp	<- sum( YX.m3.bs[ which(as.character(YX.m3.bs[, risk, with=FALSE][[1]])==factor), score.Y*w] )
														list(risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef.ref='None', stat= 'N.raw', v=tmp )
													}, by='coef']
				risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
				#	for proportions		
				setkey(risk.df, coef)
				tmp			<- unique(risk.df)[, 		{
							tmp				<- copy(predict.df)
							set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
							corisk			<- intersect(colnames(risk.df), colnames(predict.df))
							if(!length(corisk))	
								corisk<- NULL
							tmp				<- na.omit(merge(tmp[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3.bs[which( YX.m3.bs[, risk, with=FALSE][[1]]==factor ), c(risk,corisk,'w'), with=FALSE], by=risk))
							setkey(tmp, NULL)
							if(nrow(tmp))	
							{
								if( nrow(unique(tmp[, setdiff(colnames(tmp),'w'), with=FALSE]))==1 )
									tmp		<- tmp[1,]																			
								tryCatch({
								tmp[, predict:= predict(betafit.rr.bs, newdata=as.data.frame(tmp), data=na.omit(as.data.frame(YX.m3.bs[, include.colnames, with=FALSE])), type='link')]							
								tmp			<- weighted.mean(exp(tmp[,predict]), w=tmp[,w], na.rm=TRUE)
								}, error=function(e){ print(e$message); tmp<<- NA_real_ }, warning=function(w){ print(w$message); tmp<<- NA_real_ })
							}
							else
								tmp			<- NA_real_
							tmp[2]			<- sum( YX.m3.bs[ which(as.character(YX.m3.bs[, risk, with=FALSE][[1]])==factor), w] )																												
							list(risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef.ref='None', stat='prob', v=ifelse(tmp[2]<2*EPS, 0., tmp[1]) )
						},by= 'coef']
				risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))																	
				tmp			<- unique(risk.df)[, list(risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef.ref='None', stat= 'n', v=sum( YX.m3.bs[ which(as.character(YX.m3.bs[, risk, with=FALSE][[1]])==factor), w] ) ),by= 'coef']
				risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v)))
				risk.ans.bs[, bs:=bs.i]
			})
	risk.ans.bs	<- do.call('rbind',tmp)
	#
	tmp			<- risk.ans.bs[, list(stat='N', risk=risk[1], factor=factor[1], risk.ref='None', factor.ref='None', coef.ref='None', v= v[stat=='prob']*v[stat=='n'] ), by=c('bs','coef')]
	set(tmp, tmp[, which(v==0)],'v', NA_real_)
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, v, bs)))
	tmp			<- tmp[,	list(stat='P', coef=coef, risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef.ref='None', v= v/sum(v, na.rm=TRUE)),	by='bs']		
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, v, bs)))	
	tmp			<- merge( subset(tmp, stat=='P'), subset(risk.ans, stat=='PYs', select=c(coef, v)), by='coef' )
	tmp[, v:= v.x/( v.y/subset(risk.ans, stat=='PYs')[, sum(v)] )]
	set(tmp, NULL, 'stat', 'RI')
	set(tmp, tmp[,which(is.nan(v) | v==0)],'v',NA_real_)
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, v, bs)))
	#	add RR assuming no transmissions from Not.PT
	tmp			<- merge( subset(risk.ans.bs, stat=='RR'), subset(risk.df, select=c(risk, factor, risk.ref, factor.ref, PTx, PTx.ref)), by=c('risk','risk.ref','factor','factor.ref'))
	set(tmp, NULL, 'v', tmp[, v*PTx/PTx.ref])
	set(tmp, NULL, 'stat', 'RR.ptx')
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, v, bs)))
	#	add RR.term assuming no transmissions from Not.PT
	tmp			<- merge( subset(risk.ans.bs, stat=='RR.term'), subset(risk.df, select=c(risk, factor, risk.ref, factor.ref, PTx, PTx.ref)), by=c('risk','risk.ref','factor','factor.ref'))
	set(tmp, NULL, 'v', tmp[, v*PTx/PTx.ref])
	set(tmp, NULL, 'stat', 'RR.term.ptx')
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, v, bs)))	
	#	add proportion assuming no transmissions from Not.PT
	tmp			<- risk.ans.bs[, list(stat='P.ptx', risk=risk[1], factor=factor[1], risk.ref='None', factor.ref='None', coef.ref='None', expbeta= v[stat=='prob'], n=v[stat=='n'] ), by=c('bs','coef')]
	setkey(risk.df, risk, factor)
	tmp			<- merge(tmp, subset(unique(risk.df), select=c(risk, factor, PTx)), by=c('risk','factor'))
	tmp			<- merge(subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, bs)), tmp[, list(coef=coef, v= expbeta*n*PTx / sum(expbeta*n*PTx)),by='bs'], by=c('bs','coef')) 
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, v, bs)))
	#	add RI assuming no transmissions from Not.PT
	tmp			<- subset(risk.ans, stat=='PYs' )
	set(tmp, NULL, 'pPYs', tmp[, v/sum(v, na.rm=TRUE)])
	tmp			<- merge( subset(risk.ans.bs, stat=='P.ptx'), subset(tmp, select=c(risk, factor, pPYs)), by=c('risk','factor'))
	set(tmp, NULL, 'v', tmp[, v/pPYs])
	set(tmp, tmp[,which(is.nan(v) | v==0)],'v',NA_real_)
	set(tmp, NULL, 'stat', 'RI.ptx')
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, v, bs)))
	#
	tmp			<- subset(risk.ans.bs, stat=='N.raw')	
	tmp			<- tmp[, {
				list(coef=coef, risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef.ref='None', stat= 'P.raw', v=v/sum(v) )
			}, by='bs']
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, v, bs)))
	#	RI.raw
	tmp			<- merge( subset(risk.ans.bs, stat=='P.raw'), subset(risk.ans, stat=='PYs', select=c(coef, v)), by='coef' )
	tmp[, v:= v.x/( v.y/subset(risk.ans, stat=='PYs')[, sum(v)] )]
	set(tmp, NULL, 'stat', 'RI.raw')
	set(tmp, tmp[,which(is.nan(v) | v==0)],'v',NA_real_)
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, v, bs)))
	#
	tmp			<- risk.ans.bs[,	list(l95.bs=quantile(v, prob=0.025, na.rm=TRUE), u95.bs=quantile(v, prob=0.975, na.rm=TRUE)), by=c('coef','coef.ref','stat')]
	risk.ans	<- merge(risk.ans, tmp, by=c('coef','coef.ref','stat'), all.x=TRUE)
	#	
	setkey(risk.ans, stat, coef.ref, coef)
	list(risk=risk.ans, fit.or=betafit.or, fit.rr=betafit.rr, risk.bs=risk.ans.bs)
}
######################################################################################
project.athena.Fisheretal.estimate.risk.wrap.add2riskdf<- function(method.risk, risk.df, X.tables)
{	
	tp				<- regmatches(method.risk, regexpr('tp[0-9]', method.risk))
	tp				<- ifelse(length(tp), paste('.',substr(tp, 3, 3),sep=''), '')
	#tmp				<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')
	tmp				<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')
	risk.df			<- merge(risk.df, subset(X.tables$risk.table, stat==tmp, c(risk, factor, n)), by=c('risk','factor'), all.x=1)
	setnames(risk.df, 'n', 'PYs')		
	#rbind( risk.df[,	{
	#				z	<- table( YX[, risk, with=FALSE], useNA='ifany')
	#				list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
	#			},by='risk'],
	#	risk.df[,	{
	#				z	<- table( YX.wz[, risk, with=FALSE], useNA='ifany')
	#				list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX.wz')												
	#			},by='risk'])[, list(PTx= n[stat=='YX']/n[stat=='YX.wz']), by=c('risk','factor')]
	tmp				<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')
	risk.df			<- merge(risk.df, X.tables$risk.table[, list(PTx= n[stat=='YX']/n[stat==tmp]), by=c('risk','factor')], by=c('risk','factor'), all.x=1)	#	add Prob( i identified as PT to j | risk factor)	
	#PYe	
	risk.df	<- merge(risk.df, subset(X.tables$risk.table, stat=='X.msm', c(risk, factor, n)), by=c('risk','factor'), all.x=1)
	setnames(risk.df, 'n', 'PYe0')
	#PYe5
	setkey(risk.df, risk, factor)
	risk.df[, PYe5:= PYe0] 	
	tmp		<- unique(risk.df)[, sum(PYe5)/0.95*0.05]
	set(risk.df, risk.df[, which(factor==paste('U',tp,sep=''))], 'PYe5', risk.df[factor==paste('U',tp,sep=''), PYe5+tmp] )
	#PYe7
	risk.df[, PYe7:= PYe0]
	tmp		<- unique(risk.df)[, sum(PYe7)/0.93*0.07]
	set(risk.df, risk.df[, which(factor==paste('U',tp,sep=''))], 'PYe7', risk.df[factor==paste('U',tp,sep=''), PYe7+tmp] )
	#PYe3
	risk.df[, PYe3:= PYe0]
	tmp		<- unique(risk.df)[, sum(PYe3)/0.97*0.03]
	set(risk.df, risk.df[, which(factor==paste('U',tp,sep=''))], 'PYe3', risk.df[factor==paste('U',tp,sep=''), PYe3+tmp] )
	#PYe censoring by constant proportion
	if(grepl('tp[0-9]', method.risk) || grepl('TP', method.risk) )
		tmp		<- subset(X.tables$cens.table, stat=='X.msm')[, list(risk=risk, PYe0cp= n.adjbyPU) , by='factor']
	if(!grepl('tp[0-9]', method.risk) & !grepl('TP', method.risk) )
	{
		tmp		<- subset(X.tables$cens.table, stat=='X.msm')[, list(risk=risk[1], PYe0cp= sum(n.adjbyPU)) , by='factor2']
		setnames(tmp, 'factor2', 'factor')
	}					
	risk.df	<- merge(risk.df, tmp, by=c('risk','factor'), all.x=1)
	#PYe5cp
	setkey(risk.df, risk, factor)
	risk.df[, PYe5cp:= PYe0cp] 	
	if(any(risk.df[, grepl('U',factor)]))
	{
		tmp		<- unique(risk.df)[, sum(PYe5cp)/0.95*0.05]
		set(risk.df, risk.df[, which(factor==paste('U',tp,sep=''))], 'PYe5cp', risk.df[factor==paste('U',tp,sep=''), PYe5cp+tmp] )
	}
	#PYe7cp
	risk.df[, PYe7cp:= PYe0cp]
	if(any(risk.df[, grepl('U',factor)]))
	{		
		tmp		<- unique(risk.df)[, sum(PYe7cp)/0.93*0.07]
		set(risk.df, risk.df[, which(factor==paste('U',tp,sep=''))], 'PYe7cp', risk.df[factor==paste('U',tp,sep=''), PYe7cp+tmp] )
	}
	#PYe3cp
	risk.df[, PYe3cp:= PYe0cp]
	if(any(risk.df[, grepl('U',factor)]))
	{		
		tmp		<- unique(risk.df)[, sum(PYe3cp)/0.97*0.03]
		set(risk.df, risk.df[, which(factor==paste('U',tp,sep=''))], 'PYe3cp', risk.df[factor==paste('U',tp,sep=''), PYe3cp+tmp] )
	}
	risk.df
}
######################################################################################
project.athena.Fisheretal.estimate.risk.wrap<- function(YX, X.tables, tperiod.info, plot.file.or=NA, bs.n=1e3, resume=TRUE, save.file=NA, method.risk=NA)
{	
	if(resume & !is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(!resume || is.na(save.file) || inherits(readAttempt, "try-error"))
	{
		cat(paste('\nregression on data set by method', method.risk))
		if(grepl('m21st.cas',method.risk) | grepl('m2B1st.cas',method.risk) | grepl('m21stMv.cas',method.risk) | grepl('m2B1stMv.cas',method.risk))
		{  
			#	use weights?
			if(grepl('now',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i*w.in])
			if(grepl('wstar',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i])	
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk))
				set(YX, NULL, 'w', 1.)
			#	censoring adjustment
			if(grepl('cens', method.risk))
			{
				tmp			<- ifelse(grepl('m2wmx',method.risk),"CD41st.tperiod","CD4b.tperiod")	
				set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
				tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')	
				if(grepl('censp', method.risk))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			tmp				<- ifelse(grepl('m2wmx',method.risk),"CD41st","CD4b")
			set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
			#	sequence adjustment (not censoring)
			if(grepl('adj', method.risk))
			{
				tmp			<- subset( X.tables$risk.table, factor%in%YX[, levels(stage)] )
				tmp			<- tmp[,  list(risk=risk, factor=factor, n=n, p=n/sum(n)), by='stat']
				if(!grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.seq']),by=c('risk','factor')]
				if(grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.clu']),by=c('risk','factor')]
				tmp			<- data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] )
				set(tmp, tmp[, which(w.b>8.)], 'w.b', 8.)
				YX			<- merge( YX, tmp, by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}		
			if(grepl('Mv', method.risk))	
			{
				include.colnames<- c('score.Y','w','stage','t','t.Age','t.RegionHospital')
				formula			<- 'score.Y ~ bs(t.Age, knots=c(30,45), degree=1)+bs(t, knots=c(2007,2010), degree=2)+stage+t.RegionHospital-1'
				predict.df		<- data.table(	stage=factor('ART1.su.Y', levels=YX[, levels(stage)]), 												
									t.Age=subset(YX, stage=='ART1.su.Y')[, mean(t.Age, na.rm=TRUE)], t=subset(YX, stage=='ART1.su.Y')[, mean(t, na.rm=TRUE)],
									t.RegionHospital=factor('Amst', levels=YX[, levels(t.RegionHospital)]),	w=1.)										
			}
			if(!grepl('Mv', method.risk))
			{
				formula			<- 'score.Y ~ stage-1'
				include.colnames<- c('score.Y','w','stage')						
				predict.df		<- data.table(stage=factor('ART1.su.Y', levels=YX[, levels(stage)]), w=1.)				
			}						
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)], risk.ref='stage', factor.ref='ART1.su.Y')
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			if(grepl('Mv', method.risk))
			{
				setkey(risk.df, risk, factor)
				risk.df		<- merge(risk.df, unique(risk.df)[, {
									t.Age		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t.Age, na.rm=TRUE )]
									t.Age		<- ifelse(is.nan(t.Age), YX[, mean( t.Age, na.rm=TRUE )], t.Age)
									t			<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t, na.rm=TRUE )]
									t			<- ifelse(is.nan(t), YX[, mean( t, na.rm=TRUE )], t)
									list(t.Age=t.Age, t=t, t.RegionHospital='Amst')
								}, by=c('risk','factor')], by=c('risk','factor'))						
			}	
			#risk.df	<- project.athena.Fisheretal.estimate.risk.wrap.add2riskdf(method.risk, risk.df, X.tables)
			#ans		<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )	
			ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
		}	
		if(grepl('m2t.cas',method.risk) | grepl('m2Bt.cas',method.risk) | grepl('m2tMv.cas',method.risk) | grepl('m2BtMv.cas',method.risk) )
		{  
			#	use cluster weights?
			if(grepl('now',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i*w.in])
			if(grepl('wstar',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i])	
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk))
				set(YX, NULL, 'w', 1.)
			#	censoring adjustment
			if(grepl('cens', method.risk))
			{
				tmp			<- ifelse(grepl('m2wmx',method.risk),"CD41st.tperiod","CD4b.tperiod")	
				set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
				tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')	
				if(grepl('censp', method.risk))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				if(!grepl('censp', method.risk))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			tmp				<- ifelse(grepl('m2wmx',method.risk),"CD41st","CD4b")
			set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
			#	sequence adjustment (not censoring)
			if(grepl('adj', method.risk))
			{
				tmp			<- subset( X.tables$risk.table, factor%in%YX[, levels(stage)] )
				tmp			<- tmp[,  list(risk=risk, factor=factor, n=n, p=n/sum(n)), by='stat']
				if(!grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.seq']),by=c('risk','factor')]
				if(grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.clu']),by=c('risk','factor')]
				tmp			<- data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] )
				set(tmp, tmp[, which(w.b>8.)], 'w.b', 8.)
				YX			<- merge( YX, tmp, by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			if(grepl('Mv', method.risk))	
			{
				include.colnames<- c('score.Y','w','stage','t','t.Age','t.RegionHospital')
				formula			<- 'score.Y ~ bs(t.Age, knots=c(30,45), degree=1)+bs(t, knots=c(2007,2010), degree=2)+stage+t.RegionHospital-1'
				predict.df		<- data.table(	stage=factor('ART1.su.Y', levels=YX[, levels(stage)]), 												
						t.Age=subset(YX, stage=='ART.su.Y')[, mean(t.Age, na.rm=TRUE)], t=subset(YX, stage=='ART.su.Y')[, mean(t, na.rm=TRUE)],
						t.RegionHospital=factor('Amst', levels=YX[, levels(t.RegionHospital)]),	w=1.)										
			}
			if(!grepl('Mv', method.risk))
			{
				formula			<- 'score.Y ~ stage-1'
				include.colnames<- c('score.Y','w','stage')						
				predict.df		<- data.table(stage=factor('ART.su.Y', levels=YX[, levels(stage)]), w=1.)				
			}						
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)], risk.ref='stage', factor.ref='ART.su.Y')
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]			
			if(grepl('Mv', method.risk))
			{
				setkey(risk.df, risk, factor)
				risk.df		<- merge(risk.df, unique(risk.df)[, {
									t.Age		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t.Age, na.rm=TRUE )]
									t.Age		<- ifelse(is.nan(t.Age), YX[, mean( t.Age, na.rm=TRUE )], t.Age)
									t			<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t, na.rm=TRUE )]
									t			<- ifelse(is.nan(t), YX[, mean( t, na.rm=TRUE )], t)
									list(t.Age=t.Age, t=t, t.RegionHospital='Amst')
								}, by=c('risk','factor')], by=c('risk','factor'))						
			}		
			#risk.df	<- project.athena.Fisheretal.estimate.risk.wrap.add2riskdf(method.risk, risk.df, X.tables)			
			#ans		<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
		}		
		if(grepl('m2Bwmx.cas', method.risk) | grepl('m2BwmxMv.cas', method.risk) | grepl('m2wmx.cas', method.risk) | grepl('m2wmxMv.cas', method.risk))
		{  
			#	use cluster weights?
			if(grepl('now',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i*w.in])
			if(grepl('wstar',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i])	
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk))
				set(YX, NULL, 'w', 1.)
			#	censoring adjustment
			if(grepl('cens', method.risk))
			{
				tmp			<- ifelse(grepl('m2wmx',method.risk),"CD41st.tperiod","CD4b.tperiod")	
				set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
				tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')	
				if(grepl('censp', method.risk))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				if(!grepl('censp', method.risk))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )	
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
				tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')			
				#adjust		<- X.tables$cens.table[, list(n.adjbyPU= sum(n.adjbyPU), factor=factor2), by=c('risk','factor2','stat')]				
				#tmp			<- tmp[, list(factor=factor[stat=='X.msm'], n.adjbyPU.Top= n.adjbyPU[stat=='X.msm'], n.adjbyPU.Bottom= n.adjbyPU[stat=='X.seq'], n.adjbyPU.BottomS= sum(n.adjbyPU[stat=='X.seq']), n.adjbyPU.TopS= sum(n.adjbyPU[stat=='X.msm'])  ), by= 'risk']				
				#tmp[, w.b:= (n.adjbyPU.Top / n.adjbyPU.TopS) / (n.adjbyPU.Bottom / n.adjbyPU.BottomS) ]		
				#adjust		<- adjust[, list(factor=factor[stat=='X.msm'], w.b= n.adjbyPU[stat=='X.msm'] / n.adjbyPU[stat==tmp] * sum(n.adjbyPU[stat==tmp]) / sum(n.adjbyPU[stat=='X.msm'])  ), by= 'risk']				
			}
			tmp				<- ifelse(grepl('m2wmx',method.risk),"CD41st","CD4b")
			set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
			#	sequence adjustment (not censoring)
			if(grepl('adj', method.risk))
			{
				tmp			<- subset( X.tables$risk.table, factor%in%YX[, levels(stage)] )
				tmp			<- tmp[,  list(risk=risk, factor=factor, n=n, p=n/sum(n)), by='stat']
				if(!grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.seq']),by=c('risk','factor')]
				if(grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.clu']),by=c('risk','factor')]
				tmp			<- data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] )
				set(tmp, tmp[, which(w.b>8.)], 'w.b', 8.)
				YX			<- merge( YX, tmp, by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}					
			if(grepl('Mv', method.risk))	
			{
				#ggplot(YX, aes(x = t, y = score.Y)) + geom_point(size=0.4) + stat_smooth(method.risk = "loess") + facet_grid(. ~ stage, margins=TRUE)	#d=2, 2007,2010
				#ggplot(YX, aes(x = t.Age, y = score.Y)) + geom_point(size=0.4) + stat_smooth(method.risk = "loess") + facet_grid(. ~ stage, margins=TRUE)#d=1, 30,45				
				include.colnames<- c('score.Y','w','stage','t','t.Age','t.RegionHospital')
				formula			<- 'score.Y ~ bs(t.Age, knots=c(30,45), degree=1)+bs(t, knots=c(2007,2010), degree=2)+stage+t.RegionHospital-1'
				predict.df		<- data.table(	stage=factor('ART.suA.Y', levels=YX[, levels(stage)]), 												
												t.Age=subset(YX, stage=='ART.suA.Y')[, mean(t.Age, na.rm=TRUE)], t=subset(YX, stage=='ART.suA.Y')[, mean(t, na.rm=TRUE)],
												t.RegionHospital=factor('Amst', levels=YX[, levels(t.RegionHospital)]),	w=1.)										
			}
			if(!grepl('Mv', method.risk))
			{
				formula			<- 'score.Y ~ stage-1'
				include.colnames<- c('score.Y','w','stage')						
				predict.df		<- data.table(stage=factor('ART.suA.Y', levels=YX[, levels(stage)]), w=1.)				
			}			
			#
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)], risk.ref='stage', factor.ref='ART.suA.Y')
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			if(grepl('Mv', method.risk))
			{
				setkey(risk.df, risk, factor)
				risk.df		<- merge(risk.df, unique(risk.df)[, {
									t.Age		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t.Age, na.rm=TRUE )]
									t.Age		<- ifelse(is.nan(t.Age), YX[, mean( t.Age, na.rm=TRUE )], t.Age)
									t			<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t, na.rm=TRUE )]
									t			<- ifelse(is.nan(t), YX[, mean( t, na.rm=TRUE )], t)
									list(t.Age=t.Age, t=t, t.RegionHospital='Amst')
								}, by=c('risk','factor')], by=c('risk','factor'))						
			}			
			#ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans				<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
		}	
		if(grepl('m2Bt.tp', method.risk) | grepl('m2BtMv.tp', method.risk))
		{
			YXf				<- copy(YX)
			tp				<- regmatches(method.risk, regexpr('tp[0-9]', method.risk))
			cat(paste('\nprocess time period',tp))
			tp				<- substr(tp, 3, 3)			
			YX				<- subset(YX, t.period==tp)						
			tmp				<- ifelse(grepl('m2wmx',method.risk),"CD41st.tperiod","CD4b.tperiod")
			set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
			tmp				<- ifelse(grepl('m2wmx',method.risk),"CD41st","CD4b")
			set(YXf, NULL, 'stage', factor(as.character(YXf[[tmp]])))
			#YX.wz			<- copy(YX)
			#YX				<- subset(YX, score.Y>0.)
			#	use cluster weights?
			if(grepl('now',method.risk))
			{
				set(YX, NULL, 'w', YX[, w/w.i*w.in])
				set(YXf, NULL, 'w', YXf[, w/w.i*w.in])
			}					
			if(grepl('wstar',method.risk))
			{
				set(YX, NULL, 'w', YX[, w/w.i])
				set(YXf, NULL, 'w', YXf[, w/w.i])
			}						
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk))
			{
				set(YX, NULL, 'w', 1.)
				set(YXf, NULL, 'w', 1.)
			}					
			if(grepl('Mv', method.risk))	
			{
				include.colnames<- c('score.Y','w','stage','t.Age')
				formula			<- 'score.Y ~ bs(t.Age, knots=c(30,45), degree=1)+stage-1'				
				predict.df		<- data.table(	stage=factor(paste('ART.su.Y',tp,sep='.'), levels=YX[, levels(stage)]), 												
												t.Age=subset(YX, stage==paste('ART.su.Y',tp,sep='.'))[, mean(t.Age, na.rm=TRUE)], w=1.)										
			}
			if(!grepl('Mv', method.risk))
			{				
				formula			<- 'score.Y ~ stage-1'
				include.colnames<- c('score.Y','w','stage')						
				predict.df		<- data.table(stage=factor(paste('ART.su.Y',tp,sep='.'), levels=YX[, levels(stage)]), w=1.)
			}
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.su.Y',tp,sep='.'))
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=YX[, levels(stage)], risk.ref='stage', factor.ref=paste('U',tp,sep='.')) )
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('l350', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			if(grepl('Mv', method.risk))
			{
				setkey(risk.df, risk, factor)
				risk.df		<- merge(risk.df, unique(risk.df)[, {
									t.Age		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t.Age, na.rm=TRUE )]
									t.Age		<- ifelse(is.nan(t.Age), YX[, mean( t.Age, na.rm=TRUE )], t.Age)
									list(t.Age=t.Age)
								}, by=c('risk','factor')], by=c('risk','factor'))						
			}			
			#risk.df		<- project.athena.Fisheretal.estimate.risk.wrap.add2riskdf(method.risk, risk.df, X.tables)
			#ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c(0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0, 0) )
			#ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c(0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0) )
			ans				<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c( 0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0)  )
		}
		if(	grepl('m2wmx.tp', method.risk) | grepl('m2wmxMv.tp', method.risk) |
			grepl('m2wmx.wtn.tp', method.risk) | grepl('m2wmxMv.wtn.tp', method.risk) |
			grepl('m2Bwmx.tp', method.risk) | grepl('m2BwmxMv.tp', method.risk) |
			grepl('m2Bwmx.wtn.tp', method.risk) | grepl('m2BwmxMv.wtn.tp', method.risk) |
			grepl('m2Cwmx.tp', method.risk) | grepl('m2CwmxMv.tp', method.risk) |
			grepl('m2Cwmx.wtn.tp', method.risk) | grepl('m2CwmxMv.wtn.tp', method.risk) )
		{			
			YXf				<- copy(YX)
			tp				<- regmatches(method.risk, regexpr('tp[0-9]', method.risk))
			cat(paste('\nprocess time period',tp))
			tp				<- substr(tp, 3, 3)
			YX				<- subset(YX, t.period==tp)						
			if( grepl('m2wmx',method.risk) )
				tmp			<- c("CD41st.tperiod","CD41st")
			if( grepl('m2Bwmx',method.risk) )
				tmp			<- c("CD4b.tperiod", "CD4b")	
			if( grepl('m2Cwmx',method.risk) )
				tmp			<- c("CD4c.tperiod","CD4c")	
			set(YX, NULL, 'stage', factor(as.character(YX[[tmp[1]]])))
			set(YXf, NULL, 'stage', factor(paste(as.character(YXf[[tmp[2]]]),'.',tp,sep='')))			
			#	use cluster weights?
			#	use cluster weights?
			if(grepl('wtn',method.risk))
			{
				cat(paste('\nsetting likelihood to likelihood of pair / number transmission intervals'))
				set(YX, NULL, 'score.Y', YX[, score.Y*w.tn])
				set(YX, NULL, c('w','w.tn'), 1.)
				set(YXf, NULL, 'score.Y', YXf[, score.Y*w.tn])
				set(YXf, NULL, c('w','w.tn'), 1.)
			}								
			if(grepl('now',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i*w.in])	
			if(grepl('wstar',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i])				
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk) & !grepl('wtn',method.risk))
				set(YX, NULL, 'w', 1.)				
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.suA.Y',tp,sep='.'))
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=YX[, levels(stage)], risk.ref='stage', factor.ref=paste('U',tp,sep='.')) )
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('l350', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			ans			<- project.athena.Fisheretal.estimate.risk.core.Wallinga(YX, YXf, X.tables, method.risk, risk.df, bs.n=bs.n, use.YXf=1 )
			#ans		<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c(0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0) )
			#ans		<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c( 0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0)  )
		}
		if(grepl('m3.tnic',method.risk) & !grepl('m3.tnicNo',method.risk) & !grepl('m3.tnicv',method.risk))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			#	use cluster weights?
			if(grepl('now',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i*w.in])	
			if(grepl('wstar',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i])	
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk))
				set(YX, NULL, 'w', 1.)
			#	censoring adjustment
			if(grepl('cens', method.risk))
			{				
				YX[, ART.ntstage.c.tperiod:= YX[, paste(ART.ntstage.c, t.period, sep='.')]]
				set(YX, YX[,which(is.na(ART.ntstage.c))], 'ART.ntstage.c.tperiod', NA_character_)
				set(YX, NULL, 'stage', YX[, factor(as.character(ART.ntstage.c.tperiod))])
				tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')				
				if(grepl('censp', method.risk))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				if(!grepl('censp', method.risk))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )
				tmp			<- data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] )
				set(tmp, tmp[, which(w.b>8.)], 'w.b', 8.)
				YX			<- merge( YX, tmp, by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			set(YX, NULL, 'stage', YX[, factor(as.character(ART.ntstage.c))])
			#	sequence adjustment (not censoring)
			if(grepl('adj', method.risk))
			{
				tmp			<- subset( X.tables$risk.table, factor%in%YX[, levels(stage)] )
				tmp			<- tmp[,  list(risk=risk, factor=factor, n=n, p=n/sum(n)), by='stat']
				if(!grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.seq']),by=c('risk','factor')]
				if(grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.clu']),by=c('risk','factor')]
				tmp			<- data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] )
				set(tmp, tmp[, which(w.b>8.)], 'w.b', 8.)
				YX			<- merge( YX, tmp, by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			if(!grepl('MV', method.risk))
			{
				formula			<- 'score.Y ~ stage-1'
				include.colnames<- c('score.Y','w','stage')
				predict.df		<- data.table(stage=factor('ART.3.NRT.PI', levels=YX[, levels(stage)]), w=1.)				
			}
			if(grepl('MV', method.risk))
			{
				formula			<- 'score.Y ~ bs(t, knots=c(2007,2010), degree=2)+bs(t.Age, knots=c(30,45), degree=1)+stage+t.RegionHospital-1'
				include.colnames<- c('score.Y','w','stage','t','t.Age','t.RegionHospital')
				predict.df		<- data.table(	stage=factor('ART.3.NRT.PI', levels=YX[, levels(stage)]), t.RegionHospital= factor('Amst', levels=YX[, levels(t.RegionHospital)]),
												t=subset(YX, stage=='ART.3.NRT.PI')[, mean(t, na.rm=TRUE)], t.Age=subset(YX, stage=='ART.3.NRT.PI')[, mean(t.Age, na.rm=TRUE)], w=1.)
			}			
			risk.df			<- data.table(risk='stage', factor=YX[, levels(stage)], risk.ref='stage', factor.ref= 'ART.3.NRT.PI')
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			risk.df[, coef.ref:=paste(risk.df[,risk.ref],risk.df[,factor.ref],sep='')]			
			if(grepl('MV', method.risk))
			{
				setkey(risk.df, risk, factor)
				risk.df		<- merge(risk.df, unique(risk.df)[, {
									t.Age		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t.Age, na.rm=TRUE )]
									t.Age		<- ifelse(is.nan(t.Age), YX[, mean( t.Age, na.rm=TRUE )], t.Age)
									t			<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t, na.rm=TRUE )]
									t			<- ifelse(is.nan(t), YX[, mean( t, na.rm=TRUE )], t)
									list(t.Age=t.Age, t=t, t.RegionHospital='Amst')
								}, by=c('risk','factor')], by=c('risk','factor'))						
			}			
			#risk.df	<- project.athena.Fisheretal.estimate.risk.wrap.add2riskdf(method.risk, risk.df, X.tables)				
			#ans		<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.l= c(0.3, 0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0) )
			ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.l= c(0.3, 0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0)  )
		}
		if(grepl('m3.tnicv',method.risk) & !grepl('m3.tnicvNo',method.risk))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			#	use cluster weights?
			if(grepl('now',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i*w.in])	
			if(grepl('wstar',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i])	
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk))
				set(YX, NULL, 'w', 1.)
			#	censoring adjustment
			if(grepl('cens', method.risk))
			{				
				YX[, ART.ntstage.c.tperiod:= YX[, paste(ART.ntstage.c, t.period, sep='.')]]
				set(YX, YX[,which(is.na(ART.ntstage.c))], 'ART.ntstage.c.tperiod', NA_character_)
				set(YX, NULL, 'stage', YX[, factor(as.character(ART.ntstage.c.tperiod))])
				tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')				
				if(grepl('censp', method.risk))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )	
				tmp			<- data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] )
				set(tmp, tmp[, which(w.b>8.)], 'w.b', 8.)
				YX			<- merge( YX, tmp, by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			set(YX, NULL, 'stage', YX[, factor(as.character(ART.ntstage.c))])
			#	sequence adjustment (not censoring)
			if(grepl('adj', method.risk))
			{
				tmp			<- subset( X.tables$risk.table, factor%in%YX[, levels(stage)] )
				tmp			<- tmp[,  list(risk=risk, factor=factor, n=n, p=n/sum(n)), by='stat']
				if(!grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.seq']),by=c('risk','factor')]
				if(grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.clu']),by=c('risk','factor')]
				tmp			<- data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] )
				set(tmp, tmp[, which(w.b>8.)], 'w.b', 8.)
				YX			<- merge( YX, tmp, by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
			include.colnames<- c('score.Y','w','stage','lRNA.mx')
			predict.df		<- data.table( 	stage=factor('ART.3.NRT.PI', levels=YX[, levels(stage)]), 
											lRNA.mx=subset(YX, stage=='ART.3.NRT.PI')[, mean(lRNA.mx, na.rm=TRUE)], w=1.	)
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)])
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]			
			risk.df			<- merge(risk.df, risk.df[, {
															tmp				<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
															list(lRNA.mx=ifelse(is.nan(tmp), YX[, mean( lRNA.mx, na.rm=TRUE )], tmp))
														}, by='coef'], by='coef')			
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])
			#risk.df		<- project.athena.Fisheretal.estimate.risk.wrap.add2riskdf(method.risk, risk.df, X.tables)
			#ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n)
			ans				<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
		}		
		if(grepl('m3.btnicNo',method.risk))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			#	use cluster weights?
			if(grepl('now',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i*w.in])	
			if(grepl('wstar',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i])	
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk))
				set(YX, NULL, 'w', 1.)
			#	censoring adjustment
			if(grepl('cens', method.risk))
			{				
				YX[, ART.ntbstage.no.c.tperiod:= YX[, paste(ART.ntbstage.no.c, t.period, sep='.')]]
				set(YX, YX[,which(is.na(ART.ntbstage.no.c))], 'ART.ntbstage.no.c.tperiod', NA_character_)
				set(YX, NULL, 'stage', YX[, factor(as.character(ART.ntbstage.no.c.tperiod))])
				YX			<- subset(YX, !is.na(stage) & !is.na(lRNA.mx))
				tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')				
				if(grepl('censp', method.risk))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				YX[, w.orig:=w]
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			set(YX, NULL, 'stage', YX[, factor(as.character(ART.ntbstage.no.c))])
			YX				<- subset(YX, !is.na(stage) & !is.na(lRNA.mx))
			#	sequence adjustment (not censoring)
			if(grepl('adj', method.risk))
			{
				tmp			<- subset( X.tables$risk.table, factor%in%YX[, levels(stage)] )
				tmp			<- tmp[,  list(risk=risk, factor=factor, n=n, p=n/sum(n)), by='stat']
				if(!grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.seq']),by=c('risk','factor')]
				if(grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.clu']),by=c('risk','factor')]
				tmp			<- data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] )
				set(tmp, tmp[, which(w.b>8.)], 'w.b', 8.)
				YX			<- merge( YX, tmp, by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			if(!grepl('atnicv', method.risk) & !grepl('MV', method.risk))
			{
				formula			<- 'score.Y ~ stage-1'
				include.colnames<- c('score.Y','w','stage')
				predict.df		<- data.table( 	stage=factor('ART.3.NRT.PIB', levels=YX[, levels(stage)]), w=1.	)				
			}										
			if(grepl('atnicv', method.risk))
			{
				formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
				include.colnames<- c('score.Y','w','stage','lRNA.mx')
				predict.df		<- data.table( 	stage=factor('ART.3.NRT.PIB', levels=YX[, levels(stage)]), 
						lRNA.mx=subset(YX, stage=='ART.3.NRT.PIB')[, mean(lRNA.mx, na.rm=TRUE)], w=1.	)				
			}			
			if(grepl('MV', method.risk))
			{
				formula			<- 'score.Y ~ bs(t, knots=c(2007,2010), degree=2)+bs(t.Age, knots=c(30,45), degree=1)+stage+t.RegionHospital-1'
				include.colnames<- c('score.Y','w','stage','t','t.Age','t.RegionHospital')
				predict.df		<- data.table(	stage=factor('ART.3.NRT.PIB', levels=YX[, levels(stage)]), t.RegionHospital= factor('Amst', levels=YX[, levels(t.RegionHospital)]),
						t=subset(YX, stage=='ART.3.NRT.PIB')[, mean(t, na.rm=TRUE)], t.Age=subset(YX, stage=='ART.3.NRT.PIB')[, mean(t.Age, na.rm=TRUE)], w=1.)
			}						
			risk.df			<- data.table( risk='stage', factor=YX[, levels(stage)], risk.ref='stage', factor.ref='ART.3.NRT.PIB' )
			risk.df			<- rbind(risk.df, data.table(risk='stage', factor=YX[, levels(stage)], risk.ref='stage', factor.ref= 'ART.3.ATRIPLALIKE'))
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			risk.df[, coef.ref:=paste(risk.df[,risk.ref],risk.df[,factor.ref],sep='')]
			if(grepl('MV', method.risk))
			{
				setkey(risk.df, risk, factor)
				risk.df		<- merge(risk.df, unique(risk.df)[, {
									t.Age		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t.Age, na.rm=TRUE )]
									t.Age		<- ifelse(is.nan(t.Age), YX[, mean( t.Age, na.rm=TRUE )], t.Age)
									t			<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t, na.rm=TRUE )]
									t			<- ifelse(is.nan(t), YX[, mean( t, na.rm=TRUE )], t)
									list(t.Age=t.Age, t=t, t.RegionHospital='Amst')
								}, by=c('risk','factor')], by=c('risk','factor'))						
			}						
			if(grepl('btnicv', method.risk))
			{				
				risk.df		<- merge(risk.df, risk.df[, {
									tmp				<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
									list(lRNA.mx=ifelse(is.nan(tmp), YX[, mean( lRNA.mx, na.rm=TRUE )], tmp))
								}, by='coef'], by='coef')
			}			
			#risk.df	<- project.athena.Fisheretal.estimate.risk.wrap.add2riskdf(method.risk, risk.df, X.tables)
			#ans		<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c(0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0.3, 0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0))
			ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c( 0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0.3, 0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0))
		}
		if(grepl('m3.ind',method.risk))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			#	use cluster weights?
			if(grepl('now',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i*w.in])	
			if(grepl('wstar',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i])
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk))
				set(YX, NULL, 'w', 1.)
			if(grepl('m3.ind',method.risk))
				tmp			<- 'ART.ind.c'
			if(grepl('m3.indNo',method.risk))
				tmp			<- 'ART.ind.no.c'
			if(grepl('m3.indmx',method.risk))
				tmp			<- 'ART.indmx.c'
			if(grepl('m3.indmxNo',method.risk))
				tmp			<- 'ART.indmx.no.c'
			set(YX, NULL, 'stage', factor(as.character( YX[[tmp]] )))
			YX				<- subset(YX, !is.na(stage))
			#
			if(!grepl('MV', method.risk))
			{		
				formula			<- 'score.Y ~ stage-1'
				include.colnames<- c('score.Y','w','stage')
				predict.df		<- data.table( 	stage=factor('ART.OK', levels=YX[, levels(stage)]), w=1.	)
			}
			if(grepl('MV', method.risk))
			{
				formula			<- 'score.Y ~ bs(t, degree=3)+bs(t.Age, knots=c(30,45), degree=1)+stage+t.RegionHospital-1'
				include.colnames<- c('score.Y','w','stage','t','t.Age','t.RegionHospital')
				predict.df		<- data.table(	stage=factor('ART.OK', levels=YX[, levels(stage)]), t.RegionHospital= factor('Amst', levels=YX[, levels(t.RegionHospital)]),												
												t=subset(YX, stage=='ART.OK')[, mean(t, na.rm=TRUE)], t.Age=subset(YX, stage=='ART.OK')[, mean(t.Age, na.rm=TRUE)], w=1.)
			}	
			risk.df			<- data.table( risk='stage', factor=YX[, levels(stage)], risk.ref='stage', factor.ref='ART.OK' )
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			risk.df[, coef.ref:=paste(risk.df[,risk.ref],risk.df[,factor.ref],sep='')]
			if(grepl('MV', method.risk))
			{
				setkey(risk.df, risk, factor)
				risk.df		<- merge(risk.df, unique(risk.df)[, {
																	t.Age		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t.Age, na.rm=TRUE )]
																	t.Age		<- ifelse(is.nan(t.Age), YX[, mean( t.Age, na.rm=TRUE )], t.Age)
																	t			<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t, na.rm=TRUE )]
																	t			<- ifelse(is.nan(t), YX[, mean( t, na.rm=TRUE )], t)
																	list(t.Age=t.Age, t=t, t.RegionHospital='Amst')
																}, by=c('risk','factor')], by=c('risk','factor'))						
			}						
			#risk.df	<- project.athena.Fisheretal.estimate.risk.wrap.add2riskdf(method.risk, risk.df, X.tables)
			#ans		<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c(0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0, 0))
			ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c( 0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0, 0)  )
		}
		if(grepl('m3.n3mx',method.risk))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			#	use cluster weights?
			if(grepl('now',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i*w.in])	
			if(grepl('wstar',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i])
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk))
				set(YX, NULL, 'w', 1.)
			set(YX, NULL, 'stage', factor(as.character( YX[['ART.n3mx.c']] )))
			YX				<- subset(YX, !is.na(stage))
			#
			sigma.formula	<- '~ 1'
			if(!grepl('MV', method.risk))
			{
				formula			<- 'score.Y ~ stage-1'				
				include.colnames<- c('score.Y','w','stage')
				predict.df		<- data.table( 	stage=factor('ART.3.2NRT.X', levels=YX[, levels(stage)]), w=1.	)				
			}										
			if(grepl('MV', method.risk))
			{			
				formula			<- 'score.Y ~ bs(t, degree=3)+bs(t.Age, knots=c(30,45), degree=1)+stage+t.RegionHospital-1'
				include.colnames<- c('score.Y','w','stage','t','t.Age','t.RegionHospital')
				predict.df		<- data.table(	stage=factor('ART.3.2NRT.X', levels=YX[, levels(stage)]), t.RegionHospital= factor('Amst', levels=YX[, levels(t.RegionHospital)]),												
												t=subset(YX, stage=='ART.3.2NRT.X')[, mean(t, na.rm=TRUE)], t.Age=subset(YX, stage=='ART.3.2NRT.X')[, mean(t.Age, na.rm=TRUE)], w=1.)
			}						
			risk.df			<- data.table( risk='stage', factor=YX[, levels(stage)], risk.ref='stage', factor.ref='ART.3.2NRT.X' )			
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			risk.df[, coef.ref:=paste(risk.df[,risk.ref],risk.df[,factor.ref],sep='')]
			if(!grepl('MV', method.risk))	
			{
				setkey(risk.df, risk, factor)
				risk.df		<- merge(risk.df, unique(risk.df)[, {
									t.Age		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t.Age, na.rm=TRUE )]
									t.Age		<- ifelse(is.nan(t.Age), YX[, mean( t.Age, na.rm=TRUE )], t.Age)
									t			<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t, na.rm=TRUE )]
									t			<- ifelse(is.nan(t), YX[, mean( t, na.rm=TRUE )], t)
									list(t.Age=t.Age, t=t, t.RegionHospital='Amst')
								}, by=c('risk','factor')], by=c('risk','factor'))						
			}									
			if(grepl('MV', method.risk))
			{
				setkey(risk.df, risk, factor)
				risk.df		<- merge(risk.df, unique(risk.df)[, {
									t.Age		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t.Age, na.rm=TRUE )]
									t.Age		<- ifelse(is.nan(t.Age), YX[, mean( t.Age, na.rm=TRUE )], t.Age)
									t			<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t, na.rm=TRUE )]
									t			<- ifelse(is.nan(t), YX[, mean( t, na.rm=TRUE )], t)
									list(t.Age=t.Age, t=t, t.RegionHospital='Amst')
								}, by=c('risk','factor')], by=c('risk','factor'))						
			}						
			#risk.df	<- project.athena.Fisheretal.estimate.risk.wrap.add2riskdf(method.risk, risk.df, X.tables)
			#ans		<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c(0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0.3, 0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0))
			ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c( 0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0.3, 0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0)  )
		}
		if(grepl('m3.nnrtpiNo',method.risk))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			#	use cluster weights?
			if(grepl('now',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i*w.in])	
			if(grepl('wstar',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i])	
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk))
				set(YX, NULL, 'w', 1.)
			#	censoring adjustment
			if(grepl('cens', method.risk))
			{				
				YX[, ART.nnrtpi.no.c.tperiod:= YX[, paste(ART.nnrtpi.no.c, t.period, sep='.')]]
				set(YX, YX[,which(is.na(ART.nnrtpi.no.c))], 'ART.nnrtpi.no.c.tperiod', NA_character_)
				set(YX, NULL, 'stage', YX[, factor(as.character(ART.nnrtpi.no.c.tperiod))])
				YX			<- subset(YX, !is.na(stage) & !is.na(lRNA.mx))
				tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')				
				if(grepl('censp', method.risk))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				YX[, w.orig:=w]
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			set(YX, NULL, 'stage', YX[, factor(as.character(ART.nnrtpi.no.c))])
			YX				<- subset(YX, !is.na(stage))
			#	sequence adjustment (not censoring)
			if(grepl('adj', method.risk))
			{
				tmp			<- subset( X.tables$risk.table, factor%in%YX[, levels(stage)] )
				tmp			<- tmp[,  list(risk=risk, factor=factor, n=n, p=n/sum(n)), by='stat']
				if(!grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.seq']),by=c('risk','factor')]
				if(grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.clu']),by=c('risk','factor')]
				tmp			<- data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] )
				set(tmp, tmp[, which(w.b>8.)], 'w.b', 8.)
				YX			<- merge( YX, tmp, by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			if(!grepl('MV', method.risk))
			{
				#formula			<- 'score.Y ~ bs(t, degree=1)+bs(t.Age, knots=c(30,45), degree=1)+stage+t.RegionHospital-1'
				formula			<- 'score.Y ~ bs(t, degree=1)+bs(t.Age, degree=1)+stage+t.RegionHospital-1'
				include.colnames<- c('score.Y','w','stage','t','t.Age','t.RegionHospital')
				predict.df		<- data.table(	stage=factor('ART.3.NRT.PIB', levels=YX[, levels(stage)]), t.RegionHospital= factor('Amst', levels=YX[, levels(t.RegionHospital)]),												
												t=subset(YX, stage=='ART.3.NRT.PIB')[, mean(t, na.rm=TRUE)], t.Age=subset(YX, stage=='ART.3.NRT.PIB')[, mean(t.Age, na.rm=TRUE)], w=1.)				
			}										
			if(grepl('MV', method.risk))
			{
				YX				<- subset(YX, ART.F=='No' & ART.T=='No' )
				#formula			<- 'score.Y ~ bs(t, knots=c(2007,2010), degree=2)+bs(t.Age, knots=c(30,45), degree=1)+stage+t.RegionHospital-1'
				formula			<- 'score.Y ~ bs(t, degree=1)+bs(t.Age, degree=1)+stage+t.RegionHospital-1'
				include.colnames<- c('score.Y','w','stage','t','t.Age','t.RegionHospital')
				predict.df		<- data.table(	stage=factor('ART.3.NRT.PIB', levels=YX[, levels(stage)]), t.RegionHospital= factor('Amst', levels=YX[, levels(t.RegionHospital)]),						
												t=subset(YX, stage=='ART.3.NRT.PIB')[, mean(t, na.rm=TRUE)], t.Age=subset(YX, stage=='ART.3.NRT.PIB')[, mean(t.Age, na.rm=TRUE)], w=1.)
			}						
			risk.df			<- data.table( risk='stage', factor=YX[, levels(stage)], risk.ref='stage', factor.ref='ART.3.NRT.PIB' )
			risk.df			<- rbind(risk.df, data.table(risk='stage', factor=YX[, levels(stage)], risk.ref='stage', factor.ref= 'ART.3.ATRIPLALIKE'))
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			risk.df[, coef.ref:=paste(risk.df[,risk.ref],risk.df[,factor.ref],sep='')]
			if(!grepl('MV', method.risk))
			{
				setkey(risk.df, risk, factor)
				risk.df		<- merge(risk.df, unique(risk.df)[, {
									t.Age		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t.Age, na.rm=TRUE )]
									t.Age		<- ifelse(is.nan(t.Age), YX[, mean( t.Age, na.rm=TRUE )], t.Age)
									t			<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t, na.rm=TRUE )]
									t			<- ifelse(is.nan(t), YX[, mean( t, na.rm=TRUE )], t)
									list(t.Age=t.Age, t=t, t.RegionHospital='Amst')
								}, by=c('risk','factor')], by=c('risk','factor'))						
			}									
			if(grepl('MV', method.risk))
			{
				setkey(risk.df, risk, factor)
				risk.df		<- merge(risk.df, unique(risk.df)[, {
									t.Age		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t.Age, na.rm=TRUE )]
									t.Age		<- ifelse(is.nan(t.Age), YX[, mean( t.Age, na.rm=TRUE )], t.Age)
									t			<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t, na.rm=TRUE )]
									t			<- ifelse(is.nan(t), YX[, mean( t, na.rm=TRUE )], t)
									list(t.Age=t.Age, t=t, t.RegionHospital='Amst')
								}, by=c('risk','factor')], by=c('risk','factor'))						
			}						
			#risk.df	<- project.athena.Fisheretal.estimate.risk.wrap.add2riskdf(method.risk, risk.df, X.tables)
			#ans		<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c(0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0.25, 0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0))
			ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c( 0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0.25, 0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0)  )
		}
		if(grepl('m4',method.risk))				
		{  
			#	use cluster weights?
			if(grepl('now',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i*w.in])
			if(grepl('wstar',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i])
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk))
				set(YX, NULL, 'w', 1.)
			if(!grepl('No', method.risk))
				YX			<- subset(YX, !is.na(lRNA.mx) & !is.na(Acute) & CD4t!='ART.started' & substr(CD4t,1,1)!='U')
			if(grepl('No', method.risk))
				YX			<- subset(YX, !is.na(lRNA.mx) & !is.na(AcuteNo) & CD4t!='ART.started' & substr(CD4t,1,1)!='U')
			#	censoring adjustment of CD4t
			if(grepl('cens', method.risk))
			{	
				set(YX, NULL, 'stage', factor(as.character(YX[,CD4t.tperiod])))
				tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')	
				if(grepl('censp', method.risk))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				if(!grepl('censp', method.risk))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				set(tmp, tmp[, which(w.b>8.)], 'w.b', 8.)
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			set(YX, NULL, 'stage', factor(as.character(YX[,CD4t])))
			#	sequence adjustment (not censoring) of CD4t
			if(grepl('adj', method.risk))
			{
				tmp			<- subset( X.tables$risk.table, factor%in%YX[, levels(stage)] )
				tmp			<- tmp[,  list(risk=risk, factor=factor, n=n, p=n/sum(n)), by='stat']
				if(!grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.seq']),by=c('risk','factor')]
				if(grepl('clu', method.risk))
					tmp		<- tmp[, list(w.b= p[stat=='X.msm']/p[stat=='X.clu']),by=c('risk','factor')]
				tmp			<- data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] )
				set(tmp, tmp[, which(w.b>8.)], 'w.b', 8.)
				YX			<- merge( YX, tmp, by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			if(!grepl('No', method.risk) & !grepl('Mv', method.risk))
			{
				formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.5,4.5), degree=1)+stage+Acute-1'
				include.colnames<- c('score.Y','w','stage','lRNA.mx','Acute')			
				predict.df		<- data.table(	stage=factor('Dtl350', levels=YX[, levels(stage)]), 
						lRNA.mx=subset(YX, stage=='Dtl350')[, mean(lRNA.mx, na.rm=TRUE)],
						Acute=factor('No',levels=c('No','Maybe','Yes')), w=1.)						
				risk.df			<- data.table( risk='Acute', factor=c('Yes','Maybe','No'), risk.ref= 'Acute', factor.ref='No') 							
			}
			if(!grepl('No', method.risk) & grepl('Mv', method.risk))
			{
				formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+bs(t.Age, knots=c(30,45), degree=1)+bs(t, knots=c(2011), degree=2)+stage+Acute+t.RegionHospital+t2.care.t1+CDCC-1'
				include.colnames<- c('score.Y','w','stage','lRNA.mx','Acute','t.Age','t','t.RegionHospital','t2.care.t1','CDCC')			
				predict.df		<- data.table(	stage=factor('Dtl350', levels=YX[, levels(stage)]), 
						lRNA.mx=subset(YX, stage=='Dtl350')[, mean(lRNA.mx, na.rm=TRUE)],
						t.Age=subset(YX, stage=='Dtl350')[, mean(t.Age, na.rm=TRUE)],
						t=subset(YX, stage=='Dtl350')[, mean(t, na.rm=TRUE)],
						t.RegionHospital=factor('Amst', levels=YX[, levels(t.RegionHospital)]),
						t2.care.t1=factor('other', levels=YX[, levels(t2.care.t1)]),
						CDCC=factor('No', levels=YX[, levels(CDCC)]),
						Acute=factor('No',levels=c('No','Maybe','Yes')), w=1.)						
				risk.df			<- data.table( risk='Acute', factor=c('Yes','Maybe','No'), risk.ref= 'Acute', factor.ref='No') 							
			}
			if(grepl('No', method.risk))
			{
				formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.5,4.5), degree=1)+stage+AcuteNo-1'
				include.colnames<- c('score.Y','w','stage','lRNA.mx','AcuteNo')			
				predict.df		<- data.table(	stage=factor('Dtl350', levels=YX[, levels(stage)]), 
						lRNA.mx=subset(YX, stage=='Dtl350')[, mean(lRNA.mx, na.rm=TRUE)],
						AcuteNo=factor('No',levels=c('No','Maybe','Yes')), w=1.)						
				risk.df			<- data.table( risk='AcuteNo', factor=c('Yes','Maybe','No'), risk.ref= 'AcuteNo', factor.ref='No') 							
			}
			
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			risk.df[, coef.ref:=paste(risk.df[,risk.ref],risk.df[,factor.ref],sep='')]
			risk.df			<- merge(risk.df, risk.df[, {
								lRNA.mx		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
								lRNA.mx		<- ifelse(is.nan(lRNA.mx), YX[, mean( lRNA.mx, na.rm=TRUE )], lRNA.mx)
								t.Age		<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t.Age, na.rm=TRUE )]
								t.Age		<- ifelse(is.nan(t.Age), YX[, mean( t.Age, na.rm=TRUE )], t.Age)
								t			<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( t, na.rm=TRUE )]
								t			<- ifelse(is.nan(t), YX[, mean( t, na.rm=TRUE )], t)
								list(lRNA.mx=lRNA.mx, t.Age=t.Age, t=t)
							}, by='coef'], by='coef')		
			#risk.df			<- project.athena.Fisheretal.estimate.risk.wrap.add2riskdf(method.risk, risk.df, X.tables)
			#ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n)
			ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
		}
		if(grepl('m5',method.risk))				
		{  
			tp				<- regmatches(method.risk, regexpr('tp[0-9]', method.risk))
			if(length(tp))
			{
				cat(paste('\nprocess time period',tp))
				tp			<- substr(tp, 3, 3)
				YX			<- subset(YX, t.period==tp)		
				tp			<- paste('.', tp, sep='')				
			}
			if(!length(tp))
				tp			<- ''
			if(grepl('now',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i*w.in])	
			if(grepl('wstar',method.risk))
				set(YX, NULL, 'w', YX[, w/w.i])	
			if(!grepl('wstar',method.risk) & !grepl('now',method.risk))
				set(YX, NULL, 'w', 1.)
			if(!grepl('Ab',method.risk) & !grepl('Ac',method.risk))
				risk.col		<- ifelse(nchar(tp), 'tA.tperiod', 'tA')										
			if(grepl('Ac',method.risk))
				risk.col		<- ifelse(nchar(tp), 'tAc.tperiod', 'tAc')				
			factor.ref		<- paste('t<=25', tp, sep='')	
			set(YX, NULL, 'stage', factor(as.character(YX[[risk.col]])))
			#
			if(!grepl('Mv', method.risk))
			{
				formula			<- 'score.Y ~ stage-1'
				include.colnames<- c('score.Y','w','stage')
				predict.df		<- data.table( 	stage=factor(factor.ref, levels=YX[, levels(stage)]), w=1.	)				
			}										
			if(grepl('Mv', method.risk))
			{
				formula			<- 'score.Y ~ stage+t.RegionHospital-1'
				include.colnames<- c('score.Y','w','stage','t.RegionHospital')
				predict.df		<- data.table(	stage=factor(factor.ref, levels=YX[, levels(stage)]), t.RegionHospital= factor('Amst', levels=YX[, levels(t.RegionHospital)]),
												w=1.)
			}						
			risk.df			<- data.table( risk='stage', factor=YX[, levels(stage)], risk.ref='stage', factor.ref=factor.ref )			
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			risk.df[, coef.ref:=paste(risk.df[,risk.ref],risk.df[,factor.ref],sep='')]
			if(grepl('Mv', method.risk))
			{
				setkey(risk.df, risk, factor)
				risk.df		<- merge(risk.df, unique(risk.df)[, {	list( 	t.RegionHospital='Amst')		}, by=c('risk','factor')], by=c('risk','factor'))						
			}						
			#risk.df	<- project.athena.Fisheretal.estimate.risk.wrap.add2riskdf(method.risk, risk.df, X.tables)
			#ans		<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c(0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0.3, 0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0))
			ans			<- project.athena.Fisheretal.estimate.risk.core.noWadj(YX, X.tables, tperiod.info, method.risk, formula, predict.df, risk.df, include.colnames, bs.n=bs.n, gamlss.BE.limit.u=c( 0.8, 0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1), gamlss.BE.limit.l= c(0, 0)  )
		}
		if(!is.na(save.file))
		{
			ans$YX			<- YX
			ans$risk.table	<- X.tables$risk.table
			ans$cens.table	<- X.tables$cens.table
			ans$cens.table.all	<- X.tables$cens.table.all
			ans$adj.seq		<- X.tables$adj.seq
			ans$adj.clu		<- X.tables$adj.clu
			ans$adj.cens.seq<- X.tables$adj.cens.seq
			ans$adj.cens.clu<- X.tables$adj.cens.clu
			cat(paste('\nsave to file', save.file))
			save(ans, file=save.file)
		}
	}
	ans
}
######################################################################################
project.athena.Fisheretal.estimate.risk.table<- function(YX=NULL, X.den=NULL, X.msm=NULL, X.clu=NULL, tperiod.info=NULL, resume=TRUE, save.file=NA, method=NA, bs.n=100, bs.cdelta.min=2, bs.cdelta.max=3)
{	
	if(resume & !is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(!resume || is.na(save.file) || inherits(readAttempt, "try-error"))
	{
		ans	<- NULL		
		if(is.null(YX) || is.null(X.den) || is.null(X.msm) || is.null(X.clu))
		{
			cat(paste('\nreturn NULL table', method))
			return(ans)	
		}			
		cat(paste('\ntables by method', method))		
		if(substr(method,1,2)=='m2' || substr(method,1,2)=='m4' || substr(method,1,2)=='m5' || (substr(method,1,2)=='m3' && grepl('nic',method)) || grepl('m3.n3mx',method) || grepl('m3.ind',method) || grepl('m3.nnrtpiNo',method))
		{
			
			tp				<- regmatches(method, regexpr('tp[0-9]', method))
			tp				<- ifelse(length(tp), paste('.',substr(tp, 3, 3),sep=''), '')
			if(grepl('m21st',method))
			{
				factor.ref.v	<- paste('ART1.su.Y',tp,sep='')
				risktp.col		<- 'CD41st.tperiod'
				risk.col		<- 'CD41st'
			}
			if(grepl('m2B1st',method))
			{
				factor.ref.v	<- paste('ART1.su.Y',tp,sep='')
				#risktp.col		<- 'CD4t.tperiod'
				#risk.col		<- 'CD4t'
				risktp.col		<- 'CD4b.tperiod'
				risk.col		<- 'CD4b'
			}
			if(grepl('m2t',method))
			{
				factor.ref.v	<- paste('ART.su.Y',tp,sep='')
				risktp.col		<- 'CD41st.tperiod'
				risk.col		<- 'CD41st'
			}
			if(grepl('m2Bt',method))
			{
				factor.ref.v	<- paste('ART.su.Y',tp,sep='')
				risktp.col		<- 'CD4b.tperiod'
				risk.col		<- 'CD4b'
			}
			if(grepl('m2wmx',method))
			{
				factor.ref.v	<- paste('ART.suA.Y',tp,sep='')
				risktp.col		<- 'CD41st.tperiod'
				risk.col		<- 'CD41st'
			}
			if(grepl('m2Bwmx',method))
			{
				factor.ref.v	<- paste('ART.suA.Y',tp,sep='')
				risktp.col		<- 'CD4b.tperiod'
				risk.col		<- 'CD4b'
			}			
			if(grepl('m2Cwmx',method))
			{
				factor.ref.v	<- paste('ART.suA.Y2',tp,sep='')
				risktp.col		<- 'CD4c.tperiod'
				risk.col		<- 'CD4c'
			}		
			if(grepl('m3.ind',method))
			{
				factor.ref.v	<- paste('ART.OK',tp,sep='')
				if(grepl('m3.ind',method))
					risk.col		<- 'ART.ind.c'					
				if(grepl('m3.indNo',method))
					risk.col		<- 'ART.ind.no.c'									
				if(grepl('m3.indmx',method))
					risk.col		<- 'ART.indmx.c'					
				if(grepl('m3.indmxNo',method))
					risk.col		<- 'ART.indmx.no.c'									
				risktp.col		<- paste(risk.col,'.tperiod',sep='')
				YX				<- YX[ !is.na( YX[[risk.col]] ),]
				X.den			<- X.den[ !is.na( X.den[[risk.col]] ),]
				X.clu			<- X.clu[ !is.na( X.clu[[risk.col]] ),]
				X.msm			<- X.msm[ !is.na( X.msm[[risk.col]] ),]
				set(YX, NULL, risktp.col, paste( YX[, risk.col, with=FALSE][[1]], YX[, t.period], sep='.'))				
				set(X.den, NULL, risktp.col, paste( X.den[, risk.col, with=FALSE][[1]], X.den[, t.period], sep='.'))
				set(X.clu, NULL, risktp.col, paste( X.clu[, risk.col, with=FALSE][[1]], X.clu[, t.period], sep='.'))
				set(X.msm, NULL, risktp.col, paste( X.msm[, risk.col, with=FALSE][[1]], X.msm[, t.period], sep='.'))								
			}				
			if(grepl('m3.n3mx',method) & !grepl('No',method))
			{
				factor.ref.v	<- paste('ART.3.2NRT.X',tp,sep='')
				risktp.col		<- 'ART.n3mx.c.tperiod'
				risk.col		<- 'ART.n3mx.c'
				YX[, ART.n3mx.c.tperiod:= YX[, paste(ART.n3mx.c, t.period, sep='.')]]
				YX				<- subset(YX, !is.na(ART.n3mx.c))
				X.den[, ART.n3mx.c.tperiod:= X.den[, paste(ART.n3mx.c, t.period, sep='.')]]
				X.den			<- subset(X.den, !is.na(ART.n3mx.c))
				X.clu[, ART.n3mx.c.tperiod:= X.clu[, paste(ART.n3mx.c, t.period, sep='.')]]
				X.clu			<- subset(X.clu, !is.na(ART.n3mx.c))
				X.msm[, ART.n3mx.c.tperiod:= X.msm[, paste(ART.n3mx.c, t.period, sep='.')]]
				X.msm			<- subset(X.msm, !is.na(ART.n3mx.c))
			}	
			if(grepl('m3.nnrtpiNo',method))
			{
				factor.ref.v	<- paste('ART.3.NRT.PIB',tp,sep='')
				risktp.col		<- 'ART.nnrtpi.no.c.tperiod'
				risk.col		<- 'ART.nnrtpi.no.c'
				YX[, ART.nnrtpi.no.c.tperiod:= YX[, paste(ART.nnrtpi.no.c, t.period, sep='.')]]
				YX				<- subset(YX, !is.na(ART.nnrtpi.no.c))
				X.den[, ART.nnrtpi.no.c.tperiod:= X.den[, paste(ART.nnrtpi.no.c, t.period, sep='.')]]
				X.den			<- subset(X.den, !is.na(ART.nnrtpi.no.c))
				X.clu[, ART.nnrtpi.no.c.tperiod:= X.clu[, paste(ART.nnrtpi.no.c, t.period, sep='.')]]
				X.clu			<- subset(X.clu, !is.na(ART.nnrtpi.no.c))
				X.msm[, ART.nnrtpi.no.c.tperiod:= X.msm[, paste(ART.nnrtpi.no.c, t.period, sep='.')]]
				X.msm			<- subset(X.msm, !is.na(ART.nnrtpi.no.c))
			}
			if(grepl('m4.Bwmx',method))
			{
				#	adjust CD4t -- cannot adjust VL or the independent acute categories
				factor.ref.v	<- paste('Dtg500',tp,sep='')
				risktp.col		<- 'CD4t.tperiod'
				risk.col		<- 'CD4t'
			}			
			if(grepl('m5.tA',method))
			{				
				factor.ref.v	<- paste('t<=45',tp,sep='')
				risktp.col		<- 'tA.tperiod'
				risk.col		<- 'tA'
			}			
			if(grepl('m5.tAb',method))
			{				
				factor.ref.v	<- paste('t<=45',tp,sep='')
				risktp.col		<- 'tAb.tperiod'
				risk.col		<- 'tAb'
			}			
			if(grepl('m5.tAc',method))
			{				
				factor.ref.v	<- paste('t<=50',tp,sep='')
				risktp.col		<- 'tAc.tperiod'
				risk.col		<- 'tAc'
			}						
			if(grepl('m5.tiA',method))
			{				
				factor.ref.v	<- paste('t<=45-i<=45',tp,sep='')
				risktp.col		<- 'tiA.tperiod'
				risk.col		<- 'tiA'
			}			
			if(grepl('m5.tiAb',method))
			{				
				factor.ref.v	<- paste('t<=45-i<=45',tp,sep='')
				risktp.col		<- 'tiAb.tperiod'
				risk.col		<- 'tiAb'
			}			
			if(grepl('m5.tiAc',method))
			{				
				factor.ref.v	<- paste('t<=50-i<=50',tp,sep='')
				risktp.col		<- 'tiAc.tperiod'
				risk.col		<- 'tiAc'
			}						
			gc()
			#YX<- copy(YX.s); X.clu<- copy(X.clu.s); X.den<- copy(X.seq.s); X.msm<- copy(X.msm.s)
			#	get cens.table
			set(YX, NULL, 'stage', YX[[risktp.col]])				
			set(X.clu, NULL, 'stage', X.clu[[risktp.col]])
			set(X.den, NULL, 'stage', X.den[[risktp.col]])
			set(X.msm, NULL, 'stage', X.msm[[risktp.col]])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])									
			set(X.den, NULL, 'stage', X.den[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.clu, NULL, 'stage', X.clu[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(YX,    NULL, 'stage', YX[,    factor(as.character(stage), levels=X.msm[, levels(stage)])])
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)])	
			#	get number of recipients per time period			
			ans$cens.Patient.n	<- do.call('rbind',list( 	
											YX[, list(stat='YX', Patient.n=length(unique(Patient))), by='t.period'],
											X.clu[, list(stat='X.clu', Patient.n=length(unique(Patient))), by='t.period'],
											X.den[, list(stat='X.seq', Patient.n=length(unique(Patient))), by='t.period'],
											X.msm[, list(stat='X.msm', Patient.n=length(unique(Patient))), by='t.period']	))
			#	get diagnosis times of transmitters for each stage
			ans$cens.AnyPos_T1	<- risk.df[,	{
												tmp	<- subset( X.msm[ which(X.msm[[risk]]==factor), ], select=c(t.Patient, t.AnyPos_T1))
												setkey(tmp, t.Patient)			
												tmp	<- unique(tmp)
												tmp[[risk]]	<- factor[1]
												tmp
											}, by=c('risk','factor')]
			#						
			#	cens.table for all potential transmitters
			#	these will be treated as uncensored in bootstrap estimation of censoring
			#
			cat(paste('\ncompute cens.table'))
			cens.table		<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.clu[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.clu')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.seq')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')
									},by='risk']))
			cens.table[, t.period:=cens.table[, substr(factor, nchar(factor), nchar(factor))]]
			cens.table[, factor2:=cens.table[, substr(factor, 1, nchar(factor)-2)]]
			cens.table		<- merge(cens.table, cens.table[, list(factor=factor, sum=sum(n, na.rm=TRUE), p= n/sum(n, na.rm=TRUE)), by=c('stat','t.period')], by=c('stat','t.period','factor'))
			gc()
			#
			#	bootstrap approach to estimata right censoring
			#
			#	need temporarily factor without tp for X.msm
			set(X.msm, NULL, 'stage', X.msm[[risk.col]])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])
			#
			#	start boostrapping estimation of censoring
			#			
			bs.cdelta		<- runif(bs.n, bs.cdelta.min, bs.cdelta.max)
			#	table with bootstrap versions of calendar time periods for censoring
			tperiod.bs		<- copy(tperiod.info)			
			tmp				<- tperiod.bs[, max(t.period.max)]
			tperiod.bs		<- lapply(seq_along(bs.cdelta), function(b)
					{
						tmp	<- tperiod.bs[, list(t.period.min.bs= t.period.min-bs.cdelta[b], t.period.max.bs= t.period.max-bs.cdelta[b], cens.t=tmp-bs.cdelta[b], cens.delta=bs.cdelta[b]), by='t.period']
						tmp[, BS:=b]
						tmp
					})
			tperiod.bs		<- do.call('rbind', tperiod.bs)		
			#	expand tperiod.bs to include risk and factor
			tmp				<- tperiod.info[, unique(t.period)] 
			tmp				<- data.table(risk='stage',factor=X.msm[, levels(stage)])[, list(t.period=tmp), by=c('risk','factor')]
			tperiod.bs		<- merge(tperiod.bs, tmp, by='t.period',allow.cartesian=TRUE)
			# compute cens.table.bs	
			# for every bootstrap cdelta.bs
			#	compute tperiod.bs by	tperiod.start/end - cdelta.bs
			#	get table
			cat(paste('\ncompute cens.table.bs'))
			cens.table.bs	<- tperiod.bs[, {
													bs.breaks	<- c( t.period.min.bs, max(t.period.max.bs), Inf)
													tmp			<- which( X.msm[[risk]]==factor )							
													n.all		<- cut( X.msm[tmp, AnyPos_T1], breaks=bs.breaks, labels=paste(factor,seq.int(1,length(bs.breaks)-1),sep='.'), right=FALSE)
													n.all		<- table(n.all, useNA='ifany'  )
													n.pseudocens<- cut( subset(X.msm[tmp,], t.AnyPos_T1<cens.t[1])[, AnyPos_T1], breaks=bs.breaks, labels=paste(factor,seq.int(1,length(bs.breaks)-1),sep='.'), right=FALSE)
													n.pseudocens<- table(n.pseudocens, useNA='ifany'  )
													tmp			<- merge( data.table( factor.tp=rownames(n.all), n=as.numeric(unclass(n.all)) ), data.table( factor.tp=rownames(n.pseudocens), nc=as.numeric(unclass(n.pseudocens)) ), by='factor.tp', all.x=TRUE, all.y=TRUE )																										
													list(factor.tp=tmp[, factor.tp], n=tmp[, n], nc=tmp[, nc], t.period.min.bs=c(t.period.min.bs,NA_real_), t.period.max.bs=c(t.period.max.bs,NA_real_), cens.delta=cens.delta[1], cens.t=cens.t[1])
											}, by=c('BS','risk','factor')]	
			cens.table.bs	<- subset(cens.table.bs, !is.na(t.period.min.bs))
			setnames(cens.table.bs, c('factor.tp','factor'), c('factor', 'factor2'))			
			set(cens.table.bs, NULL, 'stat',  cens.table.bs[, paste('X.msm.cbs',BS,sep='')])
			cens.table.bs[, t.period:= cens.table.bs[, substr(factor, nchar(factor), nchar(factor))]]
			# end boostrapping
			#
			# revert factor to tp for X.msm
			set(X.msm, NULL, 'stage', X.msm[[risktp.col]])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])									
			#
			gc()
			ans$cens.table		<- cens.table
			ans$cens.table.bs	<- cens.table.bs
			#
			#	if 'tp' is not '', reduce data to t.period
			#
			if(tp!='')
			{
				tp			<- substr(tp, 2, 2)
				YX			<- subset(YX, t.period==tp)
				X.clu		<- subset(X.clu, t.period==tp)
				X.den		<- subset(X.den, t.period==tp)
				X.msm		<- subset(X.msm, t.period==tp)
				cens.table	<- subset(cens.table, t.period==tp)
				gc()
			}
			#	if 'tp' is '', reset 'stage' so we use the global risk.col instead of the one stratified by t.period
			if(tp=='')
			{
				set(YX, NULL, 'stage', YX[[risk.col]])				
				set(X.clu, NULL, 'stage', X.clu[[risk.col]])
				set(X.den, NULL, 'stage', X.den[[risk.col]])
				set(X.msm, NULL, 'stage', X.msm[[risk.col]])				
			}
			#	due to subsetting or resetting: need to readjust levels
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])									
			set(X.den, NULL, 'stage', X.den[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.clu, NULL, 'stage', X.clu[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(YX,    NULL, 'stage', YX[,    factor(as.character(stage), levels=X.msm[, levels(stage)])])
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref=factor.ref.v)
			if(grepl('m4.Bwmx',method))
			{
				risk.df		<- rbind(risk.df, data.table( risk='Acute', factor=c('No','Maybe','Yes'), risk.ref='Acute', factor.ref='No'))
				risk.df		<- rbind(risk.df, data.table( risk='AcuteNo', factor=c('No','Maybe','Yes'), risk.ref='AcuteNo', factor.ref='No'))
			}
			gc()
			cat(paste('\ncompute nt.table'))
			#	potential transmission intervals by stage and transmitter
			nt.table.pt		<- do.call('rbind', list( 	YX[, list(nt= length(t), stat='YX'), by=c('stage','t.Patient')],
														X.clu[, list(nt= length(t), stat='X.clu'), by=c('stage','t.Patient')],
														X.den[, list(nt= length(t), stat='X.seq'), by=c('stage','t.Patient')],
														X.msm[, list(nt= length(t), stat='X.msm'), by=c('stage','t.Patient')] 	)	)
			setnames(nt.table.pt, 'stage', 'factor')
			tmp				<- data.table( expand.grid(t.Patient=nt.table.pt[, unique(t.Patient)], factor=nt.table.pt[, unique(as.character(factor))], stat=nt.table.pt[, unique(as.character(stat))], stringsAsFactors=FALSE) )			
			nt.table.pt		<- merge(tmp, nt.table.pt, by=c('t.Patient','factor','stat'), all.x=1)
			set(nt.table.pt, nt.table.pt[, which(is.na(nt))], 'nt', 0L)
			nt.table.pt[, risk:='stage']					
			nt.table.pt		<- subset( nt.table.pt, select=c(risk, factor, t.Patient, nt, stat) )
			ans$nt.table.pt	<- nt.table.pt			
			#	potential transmission intervals by stage and recipient
			nt.table		<- do.call('rbind', list( 	YX[, list(nt= length(t), stat='YX'), by=c('stage','Patient')],
							X.clu[, list(nt= length(t), stat='X.clu'), by=c('stage','Patient')],
							X.den[, list(nt= length(t), stat='X.seq'), by=c('stage','Patient')],
							X.msm[, list(nt= length(t), stat='X.msm'), by=c('stage','Patient')] 	)	)
			setnames(nt.table, 'stage', 'factor')
			tmp				<- data.table( expand.grid(Patient=nt.table[, unique(Patient)], factor=nt.table[, unique(as.character(factor))], stat=nt.table[, unique(as.character(stat))], stringsAsFactors=FALSE) )			
			nt.table		<- merge(tmp, nt.table, by=c('Patient','factor','stat'), all.x=1)
			set(nt.table, nt.table[, which(is.na(nt))], 'nt', 0L)
			nt.table[, risk:='stage']					
			nt.table		<- subset( nt.table, select=c(risk, factor, Patient, nt, stat) )
			ans$nt.table	<- nt.table			
			#	risk tables
			cat(paste('\ncompute risk.table'))
			risk.table		<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE])
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.clu[, risk, with=FALSE])
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.clu')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE])
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.seq')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE])
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')
									},by='risk']))
			risk.table		<- merge(risk.table, risk.table[, list(factor=factor, p= n/sum(n, na.rm=TRUE)), by=c('stat','risk')], by=c('stat','risk','factor'))
			ans$risk.table	<- risk.table						
		}		
		if(method%in%c('m3.i','m3.i.clu'))
		{
			#	ART indicators only						
			#formula			<- 'score.Y ~ ART.pulse+ART.I+ART.F+ART.P+ART.A-1'
			set(YX, NULL, 'ART.pulse', YX[,factor(as.character(ART.pulse), levels=X.msm[, levels(ART.pulse)])])
			set(YX, NULL, 'ART.I', YX[,factor(as.character(ART.I), levels=X.msm[, levels(ART.I)])])
			set(YX, NULL, 'ART.F', YX[,factor(as.character(ART.F), levels=X.msm[, levels(ART.F)])])
			set(YX, NULL, 'ART.P', YX[,factor(as.character(ART.P), levels=X.msm[, levels(ART.P)])])
			set(YX, NULL, 'ART.A', YX[,factor(as.character(ART.A), levels=X.msm[, levels(ART.A)])])
			
			set(X.den, NULL, 'ART.pulse', X.den[,factor(as.character(ART.pulse), levels=X.msm[, levels(ART.pulse)])])
			set(X.den, NULL, 'ART.I', X.den[,factor(as.character(ART.I), levels=X.msm[, levels(ART.I)])])
			set(X.den, NULL, 'ART.F', X.den[,factor(as.character(ART.F), levels=X.msm[, levels(ART.F)])])
			set(X.den, NULL, 'ART.P', X.den[,factor(as.character(ART.P), levels=X.msm[, levels(ART.P)])])
			set(X.den, NULL, 'ART.A', X.den[,factor(as.character(ART.A), levels=X.msm[, levels(ART.A)])])
			
			set(X.clu, NULL, 'ART.pulse', X.clu[,factor(as.character(ART.pulse), levels=X.msm[, levels(ART.pulse)])])
			set(X.clu, NULL, 'ART.I', X.clu[,factor(as.character(ART.I), levels=X.msm[, levels(ART.I)])])
			set(X.clu, NULL, 'ART.F', X.clu[,factor(as.character(ART.F), levels=X.msm[, levels(ART.F)])])
			set(X.clu, NULL, 'ART.P', X.clu[,factor(as.character(ART.P), levels=X.msm[, levels(ART.P)])])
			set(X.clu, NULL, 'ART.A', X.clu[,factor(as.character(ART.A), levels=X.msm[, levels(ART.A)])])
							
			risk.df			<- data.table( risk= c('ART.pulse','ART.I','ART.F','ART.P','ART.A'), factor=rep('Yes',5) )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.clu[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.clu')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.seq')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))			
		}
		if(method%in%c('m3.ni','m3.ni.clu'))
		{
			#	ART indicators and number of drugs
			#	formula			<- 'score.Y ~ stage+ART.pulse+ART.I+ART.F+ART.P+ART.A-1'			
			set(YX, NULL, 'stage', YX[, ART.nDrug.c])
			set(X.clu, NULL, 'stage', X.clu[, ART.nDrug.c])
			set(X.den, NULL, 'stage', X.den[, ART.nDrug.c])
			set(X.msm, NULL, 'stage', X.msm[, ART.nDrug.c])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])
			
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(YX, NULL, 'ART.pulse', YX[,factor(as.character(ART.pulse), levels=X.msm[, levels(ART.pulse)])])
			set(YX, NULL, 'ART.I', YX[,factor(as.character(ART.I), levels=X.msm[, levels(ART.I)])])
			set(YX, NULL, 'ART.F', YX[,factor(as.character(ART.F), levels=X.msm[, levels(ART.F)])])
			set(YX, NULL, 'ART.P', YX[,factor(as.character(ART.P), levels=X.msm[, levels(ART.P)])])
			set(YX, NULL, 'ART.A', YX[,factor(as.character(ART.A), levels=X.msm[, levels(ART.A)])])
			
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'ART.pulse', X.den[,factor(as.character(ART.pulse), levels=X.msm[, levels(ART.pulse)])])
			set(X.den, NULL, 'ART.I', X.den[,factor(as.character(ART.I), levels=X.msm[, levels(ART.I)])])
			set(X.den, NULL, 'ART.F', X.den[,factor(as.character(ART.F), levels=X.msm[, levels(ART.F)])])
			set(X.den, NULL, 'ART.P', X.den[,factor(as.character(ART.P), levels=X.msm[, levels(ART.P)])])
			set(X.den, NULL, 'ART.A', X.den[,factor(as.character(ART.A), levels=X.msm[, levels(ART.A)])])
			
			set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.clu, NULL, 'ART.pulse', X.clu[,factor(as.character(ART.pulse), levels=X.msm[, levels(ART.pulse)])])
			set(X.clu, NULL, 'ART.I', X.clu[,factor(as.character(ART.I), levels=X.msm[, levels(ART.I)])])
			set(X.clu, NULL, 'ART.F', X.clu[,factor(as.character(ART.F), levels=X.msm[, levels(ART.F)])])
			set(X.clu, NULL, 'ART.P', X.clu[,factor(as.character(ART.P), levels=X.msm[, levels(ART.P)])])
			set(X.clu, NULL, 'ART.A', X.clu[,factor(as.character(ART.A), levels=X.msm[, levels(ART.A)])])
			
			risk.df			<- data.table( risk= c(rep('stage',3), 'ART.pulse','ART.I','ART.F','ART.P','ART.A'), factor=c('ART.3','ART.l3','ART.g3', rep('Yes',5)))			
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.clu[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.clu')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.seq')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))				
		}
		if(method%in%c('m3.tni','m3.tni.clu','m3.tniv','m3.tniv.clu'))
		{
			#	'score.Y ~ stage+ART.pulse+ART.I+ART.F+ART.P+ART.A-1'
			set(YX, NULL, 'stage', YX[, ART.tnDrug.c])
			set(X.clu, NULL, 'stage', X.clu[, ART.tnDrug.c])
			set(X.den, NULL, 'stage', X.den[, ART.tnDrug.c])
			set(X.msm, NULL, 'stage', X.msm[, ART.tnDrug.c])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])
			
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(YX, NULL, 'ART.pulse', YX[,factor(as.character(ART.pulse), levels=X.msm[, levels(ART.pulse)])])
			set(YX, NULL, 'ART.I', YX[,factor(as.character(ART.I), levels=X.msm[, levels(ART.I)])])
			set(YX, NULL, 'ART.F', YX[,factor(as.character(ART.F), levels=X.msm[, levels(ART.F)])])
			set(YX, NULL, 'ART.P', YX[,factor(as.character(ART.P), levels=X.msm[, levels(ART.P)])])
			set(YX, NULL, 'ART.A', YX[,factor(as.character(ART.A), levels=X.msm[, levels(ART.A)])])
			
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'ART.pulse', X.den[,factor(as.character(ART.pulse), levels=X.msm[, levels(ART.pulse)])])
			set(X.den, NULL, 'ART.I', X.den[,factor(as.character(ART.I), levels=X.msm[, levels(ART.I)])])
			set(X.den, NULL, 'ART.F', X.den[,factor(as.character(ART.F), levels=X.msm[, levels(ART.F)])])
			set(X.den, NULL, 'ART.P', X.den[,factor(as.character(ART.P), levels=X.msm[, levels(ART.P)])])
			set(X.den, NULL, 'ART.A', X.den[,factor(as.character(ART.A), levels=X.msm[, levels(ART.A)])])
			
			set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.clu, NULL, 'ART.pulse', X.clu[,factor(as.character(ART.pulse), levels=X.msm[, levels(ART.pulse)])])
			set(X.clu, NULL, 'ART.I', X.clu[,factor(as.character(ART.I), levels=X.msm[, levels(ART.I)])])
			set(X.clu, NULL, 'ART.F', X.clu[,factor(as.character(ART.F), levels=X.msm[, levels(ART.F)])])
			set(X.clu, NULL, 'ART.P', X.clu[,factor(as.character(ART.P), levels=X.msm[, levels(ART.P)])])
			set(X.clu, NULL, 'ART.A', X.clu[,factor(as.character(ART.A), levels=X.msm[, levels(ART.A)])])
			
			risk.df			<- data.table( risk= c(rep('stage',9), 'ART.pulse','ART.I','ART.F','ART.P','ART.A'), factor=c('ART.3.NRT.PI','ART.l3','ART.g3','ART.3.NRT','ART.3.PI','ART.3.NNRT','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NNRT.PI', rep('Yes',5)))
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.clu[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.clu')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.seq')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))
		}
		if(!is.na(save.file))
		{			
			cat(paste('\nsave to file', save.file))
			save(ans, file=save.file)
		}
	}
	ans
}
######################################################################################
project.athena.Fisheretal.YX.model3.stratify.ARTriskgroups<- function(YX.m3, df.all, cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), return.only.ART=TRUE, plot.file.cascade=NA, score.Y.cut=1e-2)
{	
	#	Treatment cascade & ART risk groups
	#	prop and odds by traditional categories 	undiagnosed, 	undiagnosed&evidence acute, 	diagnosed by CD4 at diagnosis
	#	
	#	PQ: are there ART risk groups?
	#	PQ: is VL the causal pathway?		--> currently problematic with uniform infection window approach
	#	SQ: is ART pulse associated with higher transmissibility?
	#YX.m3	<- copy(YX)
	if(is.null(YX.m3))	return(NULL)
	YX.m3[, U.score:=NULL]
	if('score.Y'%in%colnames(YX.m3) && YX.m3[, !any(score.Y>1.1)])
	{
		tmp		<- YX.m3[, score.Y>0] 
		set(YX.m3, which(tmp), 'score.Y', YX.m3[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )	
	}							
	#	transmitter acute
	if(0)
	{
		set(YX.m3, YX.m3[, which(stage=='U' & t.isAcute=='Yes')], 'stage', 'UAy')
		set(YX.m3, YX.m3[, which(stage=='U' & t.isAcute=='Maybe')], 'stage', 'UAm')				
		set(YX.m3, YX.m3[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
		set(YX.m3, YX.m3[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )				
		#	stratify Diag by first CD4 after diagnosis
		tmp		<- subset(df.all, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
		setkey(tmp, Patient)
		tmp		<- unique(tmp)	
		tmp[, CD4.c:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
		cat(paste('\ntransmitters with no CD4 after diagnosis\n'))
		print( subset( tmp, is.na(CD4.c) ) )
		setnames(tmp, 'Patient', 't.Patient')
		set(tmp, tmp[, which(is.na(CD4.c))], 'CD4.c', 'D1.NA')
		YX.m3	<- merge(YX.m3, subset(tmp, select=c(t.Patient, CD4.c)), by='t.Patient')
		tmp		<- YX.m3[, which(stage=='Diag')]
		set(YX.m3, tmp, 'stage', YX.m3[tmp, CD4.c])
		YX.m3[, CD4.c:=NULL]						
	}
	if(1)
	{
		set(YX.m3, YX.m3[, which(stage=='U' & t.isAcute=='Yes')], 'stage', 'UAy')
		set(YX.m3, YX.m3[, which(stage=='U' & t.isAcute=='Maybe')], 'stage', 'UAm')		
		YX.m3	<- merge(YX.m3, YX.m3[, {
											tmp<- which(stage=='ART.started')
											list( lRNA.mx= ifelse(length(tmp), ifelse(any(!is.na(lRNA[tmp])), max(lRNA[tmp], na.rm=TRUE), NA_real_), NA_real_) )							
										}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))		
		set(YX.m3, YX.m3[, which(ART.F=='No' & stage=='ART.started' & lRNA>3)], 'ART.F', 'Yes' )
		set(YX.m3, YX.m3[, which(ART.F=='Yes' & stage=='ART.started' & lRNA<=3)], 'ART.F', 'No' )
	}
	YX.m3[, stage.orig:= stage]
	if(return.only.ART)
		YX.m3	<- subset( YX.m3,stage=='ART.started' )
	gc()
	#
	#	ART indicators
	#	do this before we collapse the missing ART indicators into 'No'
	#
	#YX.m3[, ART.ind.no.c:=NULL]
	YX.m3[, ART.ind.no.c:=NA_integer_]
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'ART.ind.no.c', 4L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'ART.ind.no.c', 5L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.T=='Yes')], 'ART.ind.no.c', 7L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'ART.ind.no.c', 6L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'ART.ind.no.c', 3L)	
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='Yes')], 'ART.ind.no.c', 2L )			#ART.I excludes pulsed
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' &  ART.I=='No' &  ART.P=='No' & ART.A=='No' &  ART.F=='No' & ART.T=='No' & is.na(ART.ind.no.c))], 'ART.ind.no.c', 1L)		
	if(!return.only.ART)
	{
		set(YX.m3, YX.m3[, which(stage=='U')], 'ART.ind.no.c', 8L)
		set(YX.m3, YX.m3[, which(stage=='UAy')], 'ART.ind.no.c', 9L)
		set(YX.m3, YX.m3[, which(stage=='UAm')], 'ART.ind.no.c', 10L)
		set(YX.m3, YX.m3[, which(stage=='Diag')], 'ART.ind.no.c', 11L)
		set(YX.m3, NULL, 'ART.ind.no.c', YX.m3[, factor(ART.ind.no.c, levels=1:11, labels=c('ART.OK','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A','ART.F','ART.T','U','UAy','UAm','Diag'))])
	}
	if(return.only.ART)
		set(YX.m3, NULL, 'ART.ind.no.c', YX.m3[, factor(ART.ind.no.c, levels=1:7, labels=c('ART.OK','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A','ART.F','ART.T'))])
	gc()
	#
	#	ART indicators in any overlapping interval; add missing indicators to ART.OK
	#	
	YX.m3[, ART.ind.c:=ART.ind.no.c]
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.ind.no.c))], 'ART.ind.c', 'ART.OK')	
	#
	#	ART indicators in any overlapping interval
	#	
	tmp	<- YX.m3[,	{
						tmp<- which(stage=='ART.started')
						list(	ART.F.mx= ifelse(any(!is.na(ART.F[tmp])), any(ART.F[tmp]=='Yes', na.rm=TRUE), NA), 
								ART.T.mx= ifelse(any(!is.na(ART.T[tmp])), any(ART.T[tmp]=='Yes', na.rm=TRUE), NA), 
								ART.I.mx= ifelse(any(!is.na(ART.I[tmp])), any(ART.I[tmp]=='Yes', na.rm=TRUE), NA),
								ART.P.mx= ifelse(any(!is.na(ART.P[tmp])), any(ART.P[tmp]=='Yes', na.rm=TRUE), NA),
								ART.A.mx= ifelse(any(!is.na(ART.A[tmp])), any(ART.A[tmp]=='Yes', na.rm=TRUE), NA) 	)
					}, by=c('t.Patient','Patient')]
	set(tmp, NULL, 'ART.P.mx', tmp[, factor(as.numeric(ART.P.mx), levels=c(1,0), labels=c('Yes','No'))])
	set(tmp, NULL, 'ART.A.mx', tmp[, factor(as.numeric(ART.A.mx), levels=c(1,0), labels=c('Yes','No'))])
	set(tmp, NULL, 'ART.T.mx', tmp[, factor(as.numeric(ART.T.mx), levels=c(1,0), labels=c('Yes','No'))])
	set(tmp, NULL, 'ART.F.mx', tmp[, factor(as.numeric(ART.F.mx), levels=c(1,0), labels=c('Yes','No'))])	
	set(tmp, NULL, 'ART.I.mx', tmp[, factor(as.numeric(ART.I.mx), levels=c(1,0), labels=c('Yes','No'))])
	YX.m3	<- merge(YX.m3, tmp, by=c('t.Patient','Patient'))
	YX.m3[, ART.indmx.no.c:=NA_integer_]
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P.mx=='Yes')], 'ART.indmx.no.c', 4L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A.mx=='Yes')], 'ART.indmx.no.c', 5L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.T.mx=='Yes')], 'ART.indmx.no.c', 7L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F.mx=='Yes')], 'ART.indmx.no.c', 6L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I.mx=='Yes')], 'ART.indmx.no.c', 3L)	
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='Yes')], 'ART.indmx.no.c', 2L )			#ART.I excludes pulsed
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' &  ART.I.mx=='No' &  ART.P.mx=='No' & ART.A.mx=='No' &  ART.F.mx=='No' & ART.T.mx=='No' & is.na(ART.indmx.no.c))], 'ART.indmx.no.c', 1L)		
	if(!return.only.ART)
	{
		set(YX.m3, YX.m3[, which(stage=='U')], 'ART.indmx.no.c', 8L)
		set(YX.m3, YX.m3[, which(stage=='UAy')], 'ART.indmx.no.c', 9L)
		set(YX.m3, YX.m3[, which(stage=='UAm')], 'ART.indmx.no.c', 10L)
		set(YX.m3, YX.m3[, which(stage=='Diag')], 'ART.indmx.no.c', 11L)
		set(YX.m3, NULL, 'ART.indmx.no.c', YX.m3[, factor(ART.indmx.no.c, levels=1:11, labels=c('ART.OK','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A','ART.F','ART.T','U','UAy','UAm','Diag'))])
	}
	if(return.only.ART)
		set(YX.m3, NULL, 'ART.indmx.no.c', YX.m3[, factor(ART.indmx.no.c, levels=1:7, labels=c('ART.OK','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A','ART.F','ART.T'))])
	gc()
	#
	#	ART indicators in any overlapping interval; add missing indicators to ART.OK
	#	
	YX.m3[, ART.indmx.c:=ART.indmx.no.c]
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.indmx.no.c))], 'ART.indmx.c', 'ART.OK')
	#
	#	Basic treatment options - see if this lowers ART.OK additionally a bit
	#
	YX.m3[, ART.n3mx.c:=as.numeric(ART.indmx.c)]
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.n3mx.c==1L & ART.nDrug<3)], 'ART.n3mx.c', 8L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.n3mx.c==1L & ART.nDrug>3)], 'ART.n3mx.c', 9L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.n3mx.c==1L & ART.nDrug==3 & !(ART.nNRT==2 & (ART.nPI>0 | ART.nNNRT>0)))], 'ART.n3mx.c', 10L)
	if(!return.only.ART)
	{
		set(YX.m3, YX.m3[, which(stage=='U')], 'ART.n3mx.c', 11L)
		set(YX.m3, YX.m3[, which(stage=='UAy')], 'ART.n3mx.c', 12L)
		set(YX.m3, YX.m3[, which(stage=='UAm')], 'ART.n3mx.c', 13L)
		set(YX.m3, YX.m3[, which(stage=='Diag')], 'ART.n3mx.c', 14L)
		set(YX.m3, NULL, 'ART.n3mx.c', YX.m3[, factor(ART.n3mx.c, levels=1:14, labels=c('ART.3.2NRT.X','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A','ART.F','ART.T','ART.l3','ART.g3','ART.3.OTH','U','UAy','UAm','Diag'))])
	}
	if(return.only.ART)
		set(YX.m3, NULL, 'ART.n3mx.c', YX.m3[, factor(ART.n3mx.c, levels=1:10, labels=c('ART.3.2NRT.X','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A','ART.F','ART.T','ART.l3','ART.g3','ART.3.OTH'))])
	gc()
	#
	#	ART indication risk factors
	#
	set(YX.m3, YX.m3[,which(stage!='ART.started' | is.na(ART.pulse))], 'ART.pulse', 'No')
	set(YX.m3, NULL, 'ART.pulse', YX.m3[, factor(ART.pulse)])	
	if( YX.m3[, any(is.na(ART.pulse))] ) stop('unexpected NA in ART.pulse')
	set(YX.m3, YX.m3[,which(is.na(ART.I) | stage!='ART.started' | ART.pulse=='Yes')], 'ART.I', 'No')
	set(YX.m3, NULL, 'ART.I', YX.m3[, factor(as.character(ART.I))])
	set(YX.m3, YX.m3[,which(is.na(ART.F) | stage!='ART.started' | ART.pulse=='Yes')], 'ART.F', 'No')
	set(YX.m3, NULL, 'ART.F', YX.m3[, factor(as.character(ART.F))])	
	set(YX.m3, YX.m3[,which(is.na(ART.T) | stage!='ART.started' | ART.pulse=='Yes')], 'ART.T', 'No')
	set(YX.m3, NULL, 'ART.T', YX.m3[, factor(as.character(ART.T))])	
	set(YX.m3, YX.m3[,which(is.na(ART.A) | stage!='ART.started' | ART.pulse=='Yes')], 'ART.A', 'No')
	set(YX.m3, NULL, 'ART.A', YX.m3[, factor(as.character(ART.A))])
	set(YX.m3, YX.m3[,which(is.na(ART.P) | stage!='ART.started' | ART.pulse=='Yes')], 'ART.P', 'No')
	set(YX.m3, NULL, 'ART.P', YX.m3[, factor(as.character(ART.P))])		
	#	number of drugs -- independent. I can set nDrug==0 to NA but then I get 'contrasts can be applied only to factors with 2' because this is perf corr with ART.I
	#	set reference group for dummy coding to ART.3  
	#YX.m3[, ART.nDrug.c:=NULL]
	YX.m3[, ART.nDrug.c:=NA_integer_]
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug<3)], 'ART.nDrug.c', 2L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3)], 'ART.nDrug.c', 1L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug>3)], 'ART.nDrug.c', 3L)	
	if(!return.only.ART)
	{
		set(YX.m3, YX.m3[, which(stage=='U')], 'ART.nDrug.c', 4L)
		set(YX.m3, YX.m3[, which(stage=='UAy')], 'ART.nDrug.c', 5L)
		set(YX.m3, YX.m3[, which(stage=='UAm')], 'ART.nDrug.c', 6L)
		set(YX.m3, YX.m3[, which(stage=='Diag')], 'ART.nDrug.c', 7L)
		set(YX.m3, NULL, 'ART.nDrug.c', YX.m3[, factor(ART.nDrug.c, levels=1:7, labels=c('ART.3','ART.l3','ART.g3','U','UAy','UAm','Diag'))])
	}
	if(return.only.ART)
		set(YX.m3, NULL, 'ART.nDrug.c', YX.m3[, factor(ART.nDrug.c, levels=1:3, labels=c('ART.3','ART.l3','ART.g3'))])
	if( YX.m3[, any(is.na(ART.nDrug.c))] ) stop('unexpected NA in ART.nDrug.c')
	gc()
	#	focus after ART initiation
	#YX.m3	<- subset( YX.m3, !stage%in%c('UAm','UAy','DAm','DAy','D1<=350','D1<=550','D1>550','D1.NA','U'))
	set(YX.m3, NULL, 'stage', YX.m3[, factor(as.character(stage))])
	#
	YX.m3
}
######################################################################################
project.athena.Fisheretal.YX.model3.ARTriskgroups<- function(YX, clumsm.info, cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), plot.file.cascade=NA, score.Y.cut=1e-2)
{
	require(betareg)
	#	Treatment cascade & ART risk groups
	#	prop and odds by traditional categories 	undiagnosed, 	undiagnosed&evidence acute, 	diagnosed by CD4 at diagnosis
	#	
	#	PQ: are there ART risk groups?
	#	PQ: is VL the causal pathway?		--> currently problematic with uniform infection window approach
	#	SQ: is ART pulse associated with higher transmissibility?
	YX.m3	<- copy(YX)
	YX.m3[, U.score:=NULL]
	set(YX.m3, YX.m3[,which(score.Y>1-score.Y.cut)], 'score.Y', 1-score.Y.cut)
	set(YX.m3, YX.m3[,which(score.Y<score.Y.cut)], 'score.Y', score.Y.cut)
	set(YX.m3, NULL, 'CDCC', YX.m3[, factor(CDCC)])
	#	transmitter acute
	set(YX.m3, YX.m3[, which(stage=='U' & t.isAcute=='Yes')], 'stage', 'UAy')
	set(YX.m3, YX.m3[, which(stage=='U' & t.isAcute=='Maybe')], 'stage', 'UAm')	
	set(YX.m3, YX.m3[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
	set(YX.m3, YX.m3[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )				
	#	stratify Diag by first CD4 after diagnosis
	tmp		<- subset(clumsm.info, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
	setkey(tmp, Patient)
	tmp		<- unique(tmp)	
	tmp[, CD4.c:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
	cat(paste('\ntransmitters with no CD4 after diagnosis\n'))
	print( subset( tmp, is.na(CD4.c) ) )
	setnames(tmp, 'Patient', 't.Patient')
	set(tmp, tmp[, which(is.na(CD4.c))], 'CD4.c', 'D1.NA')
	YX.m3	<- merge(YX.m3, subset(tmp, select=c(t.Patient, CD4.c)), by='t.Patient')
	tmp		<- YX.m3[, which(stage=='Diag')]
	set(YX.m3, tmp, 'stage', YX.m3[tmp, CD4.c])
	YX.m3[, CD4.c:=NULL]		
	YX.m3[, stage.orig:= stage]
	#	ART indication risk factors
	set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )		
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')		
	#	ART no indication: nDrug
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug!=3)], 'stage', 'ART.3n')
	#	split weight for PI and NNRT
	tmp		<- YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT>0)]
	tmp2	<- YX.m3[tmp,]
	set(tmp2, NULL, 'w', tmp2[, w/2])
	set(tmp2, NULL, 'stage', 'ART.3.NRT.PI')
	YX.m3	<- rbind(YX.m3, tmp2)
	set(YX.m3, tmp, 'w', YX.m3[tmp, w/2])
	set(YX.m3, tmp, 'stage', 'ART.3.NRT.NNRT')
	#			
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 )], 'stage', 'ART.3.NRT.PI')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nNNRT>0)], 'stage', 'ART.3.NRT.NNRT')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 )], 'stage', 'ART.3.NRT')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT==0 & ART.nNNRT>0 & ART.nPI>0)], 'stage', 'ART.3.NNRT.PI')
	#	add viral load suppressed during infection window
	#VL.cur	<- 3
	#YX.m3	<- merge(YX.m3, YX.m3[, {
	#									tmp<- which(!is.na(lRNA))
	#									list( lRNA.c= ifelse(length(tmp), ifelse( max(lRNA[tmp])>VL.cur, 'SuA.N', 'SuA.Y'), 'SuA.NA') )
	#								}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
	#	problem: one factor disappears because X'X is ow not invertible
	YX.m3	<- merge(YX.m3, YX.m3[, {
						tmp<- which(!is.na(lRNA))
						list( lRNA.c= ifelse(length(tmp), max(lRNA[tmp]), NA_real_) )
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))								
	#
	YX.m3s	<- subset( YX.m3, !stage%in%c('UAm','UAy','DAm','DAy','D1<=350','D1<=550','D1>550','U'))
	set(YX.m3s, NULL, 'stage', YX.m3s[, factor(as.character(stage))])
	#set(YX.m3s, NULL, 'lRNA.c', YX.m3s[, factor(as.character(lRNA.c))])
	#YX.m3s[, table(stage)]
	#YX.m3s[, table(lRNA.c)]
	YX.m3s.fit10 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m3s)	#dummy variables act as baseline for lRNA, so that s fine
	#YX.m3s.fit10 		<- betareg(score.Y ~ poly(lRNA.c,3)+stage-1, link='logit', weights=w, data = YX.m3s)	#poly does not like NAs
	#require(splines)
	#YX.m3s.fit10 		<- betareg(score.Y ~ ns(lRNA.c,knot=median(lRNA.c))+stage-1, link='logit', weights=w, data = YX.m3s)	#exactly same r2 as linear model in x
	#YX.m3s.fit10 		<- betareg(score.Y ~ ns(lRNA.c,knot=3)+stage-1, link='logit', weights=w, data = YX.m3s)					#unclear what ns does
	
	summary(YX.m3s.fit10)
	YX.m3s.fit10.or	<- data.table( 	F= my.or.from.logit(YX.m3s.fit10, 'stageART.F', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.F')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
									I= my.or.from.logit(YX.m3s.fit10, 'stageART.I', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.I')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
									P= my.or.from.logit(YX.m3s.fit10, 'stageART.P', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.P')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
									Pulse= my.or.from.logit(YX.m3s.fit10, 'stageART.pulse.Y', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.pulse.Y')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
									A= my.or.from.logit(YX.m3s.fit10, 'stageART.A', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.A')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
									n3= my.or.from.logit(YX.m3s.fit10, 'stageART.3n', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.3n')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
									NRT= my.or.from.logit(YX.m3s.fit10, 'stageART.3.NRT', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.3.NRT')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
									NRT.NNRT= my.or.from.logit(YX.m3s.fit10, 'stageART.3.NRT.NNRT', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.3.NRT.NNRT')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
									NNRT.PI= my.or.from.logit(YX.m3s.fit10, 'stageART.3.NNRT.PI', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.3.NNRT.PI')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962) )
	
	YX.m3.fit10.p	<- YX.m3[, list(n=round(sum(w),d=1),prop= round(sum(w)/YX.m3[,sum(w)],d=3) ),by='stage']	
	if(!is.na(plot.file.cascade))
	{
		tmp			<- data.table(stage= YX.m3s[, levels(stage)], col=sapply( rainbow(YX.m3s[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m3s		<- merge(YX.m3s, tmp, by='stage')				
		pdf(file=plot.file.cascade, w=5, h=5)
		plot( seq_len(nrow(YX.m3s)), YX.m3s[, score.Y], pch=18, col= YX.m3s[, col], cex=YX.m3s[, w^0.4])
		tmp				<- subset(YX.m3s, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m3s.fit10, tmp, type='response'))	
		start			<- c( YX.m3s.fit10$coef$mean + head( 2*sqrt(diag(vcov(YX.m3s.fit10))), length(YX.m3s.fit10$coef$mean)), YX.m3s.fit10$coef$precision)
		YX.m3s.fit10.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3s, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m3s.fit10$coef$mean - head( 2*sqrt(diag(vcov(YX.m3s.fit10))), length(YX.m3s.fit10$coef$mean)), YX.m3s.fit10$coef$precision)
		YX.m3s.fit10.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3s, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
				c( predict(YX.m3s.fit10.sup, tmp, type='response'), rev(predict(YX.m3s.fit10.slw, tmp, type='response'))), 
				border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m3s[, levels(stage)], fill=YX.m3s[, unique(col)])
		YX.m3s[, col:=NULL]
		dev.off()
		#		
	}
	list(m3=YX.m3, m3.p=YX.m3.fit10.p, m3.fit.art=YX.m3s.fit10, m3.or.art=YX.m3s.fit10.or)
}
######################################################################################
project.athena.Fisheretal.YX.model3.viralload<- function(YX, df.all, vl.suppressed=log10(1e3), acute.select=c('Yes','Maybe'), cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), score.Y.cut=1e-2)
{
	require(betareg)
	#	Treatment cascade & ART risk groups
	#	prop and odds by traditional categories 	undiagnosed, 	undiagnosed&evidence acute, 	diagnosed by CD4 at diagnosis
	#	
	#	PQ: are there ART risk groups?
	#	PQ: is VL the causal pathway?		--> currently problematic with uniform infection window approach
	#	SQ: is ART pulse associated with higher transmissibility?
	YX.m3	<- copy(YX)
	YX.m3[, U.score:=NULL]
	set(YX.m3, YX.m3[,which(score.Y>1-score.Y.cut)], 'score.Y', 1-score.Y.cut)
	set(YX.m3, YX.m3[,which(score.Y<score.Y.cut)], 'score.Y', score.Y.cut)
	#	stratify Diag by first CD4 after diagnosis
	tmp		<- subset(df.all, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
	setkey(tmp, Patient)
	tmp		<- unique(tmp)	
	tmp[, CD4.c:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
	cat(paste('\ntransmitters with no CD4 after diagnosis\n'))
	print( subset( tmp, is.na(CD4.c) ) )
	setnames(tmp, 'Patient', 't.Patient')
	set(tmp, tmp[, which(is.na(CD4.c))], 'CD4.c', 'D1.NA')
	YX.m3	<- merge(YX.m3, subset(tmp, select=c(t.Patient, CD4.c)), by='t.Patient')
	tmp		<- YX.m3[, which(stage=='Diag')]
	set(YX.m3, tmp, 'stage', YX.m3[tmp, CD4.c])
	YX.m3[, CD4.c:=NULL]	
	#	if transmitter acute, set undiagnosed to undiagnosed and acute
	set(YX.m3, YX.m3[, which(t.isAcute=='Yes')], 'stage', 'Ay')
	set(YX.m3, YX.m3[, which(t.isAcute=='Maybe')], 'stage', 'Am')
	YX.m3[, stage.orig:= stage]
	#	number of drugs
}
######################################################################################
project.athena.Fisheretal.YX.model3<- function(YX, clumsm.info, vl.suppressed=log10(1e3), acute.select=c('Yes','Maybe'), cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'))
{
	require(betareg)
	#	Treatment cascade & ART risk groups
	#	prop and odds by traditional categories 	undiagnosed, 	undiagnosed&evidence acute, 	diagnosed by CD4 at diagnosis
	#	
	#	PQ: are there ART risk groups?
	#	PQ: is VL the causal pathway?		--> currently problematic with uniform infection window approach
	#	SQ: is ART pulse associated with higher transmissibility?
	YX.m3	<- copy(YX)
	YX.m3[, U.score:=NULL]
	#	stratify Diag by first CD4 after diagnosis
	tmp		<- subset(clumsm.info, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
	setkey(tmp, Patient)
	tmp		<- unique(tmp)	
	tmp[, CD4.c:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
	cat(paste('\ntransmitters with no CD4 after diagnosis\n'))
	print( subset( tmp, is.na(CD4.c) ) )
	setnames(tmp, 'Patient', 't.Patient')
	set(tmp, tmp[, which(is.na(CD4.c))], 'CD4.c', 'D1.NA')
	YX.m3	<- merge(YX.m3, subset(tmp, select=c(t.Patient, CD4.c)), by='t.Patient')
	tmp		<- YX.m3[, which(stage=='Diag')]
	set(YX.m3, tmp, 'stage', YX.m3[tmp, CD4.c])
	YX.m3[, CD4.c:=NULL]
	YX.m3[, table(stage)]
	#	if transmitter acute, set undiagnosed to undiagnosed and acute
	set(YX.m3, YX.m3[, which(stage=='U' & t.isAcute%in%acute.select)], 'stage', 'UA')
	YX.m3[, stage.orig:= stage]
	#
	#	stratify ART by 	ART.pulse	not ART.pulse
	#
	if(0)
	{
		YX.m3[, stage.orig:= stage]	
		set( YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )
		set( YX.m3, YX.m3[, which(stage=='ART.started')], 'stage', 'ART.pulse.N' )
		set( YX.m3, NULL, 'stage', YX.m3[, factor(as.character(stage))])
		tmp			<- data.table(stage= YX.m3[, levels(stage)], col=sapply( brewer.pal(YX.m3[, nlevels(stage)], 'Set1'), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m3		<- merge(YX.m3, tmp, by='stage')
		
		YX.m3.fit1 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m3)
		summary(YX.m3.fit1)
		#
		plot( seq_len(nrow(YX.m3)), YX.m3[, score.Y], pch=18, col= YX.m3[, col], cex=YX.m3[, w^0.4])
		tmp				<- subset(YX.m3, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m3.fit1, tmp, type='response'))	
		start			<- c( YX.m3.fit1$coef$mean + head( 2*sqrt(diag(vcov(YX.m3.fit1))), length(YX.m3.fit1$coef$mean)), YX.m3.fit1$coef$precision)
		YX.m3.fit1.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m3.fit1$coef$mean - head( 2*sqrt(diag(vcov(YX.m3.fit1))), length(YX.m3.fit1$coef$mean)), YX.m3.fit1$coef$precision)
		YX.m3.fit1.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
				c( predict(YX.m3.fit1.sup, tmp, type='response'), rev(predict(YX.m3.fit1.slw, tmp, type='response'))), 
				border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m3[, levels(stage)], fill=YX.m3[, unique(col)])
		#		
		YX.m3.fit1.or	<- data.table( 	ART1.pulse.YN= my.or.from.logit(YX.m3.fit1, 'stageART.pulse.Y', 'stageART.pulse.N', subset(YX.m2, stage=='ART.pulse.Y')[, sum(w)], subset(YX.m2, stage=='ART.pulse.N')[, sum(w)], 1.962),
				ART1.pulse.Y= my.or.from.logit(YX.m3.fit1, 'stageART.pulse.Y', 'stageU', subset(YX.m2, stage=='ART.pulse.Y')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				ART1.pulse.N= my.or.from.logit(YX.m3.fit1, 'stageART.pulse.N', 'stageU', subset(YX.m2, stage=='ART.pulse.N')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962)		)
		#  	 	ART1.pulse.YN 	ART1.pulse.Y 	ART1.pulse.N
		#1:     1.489437     	0.620418    	0.4165452
		#2:     0.684947     	0.527689    	0.3487689
		#3:     3.238825     	0.729442    	0.4974925
		YX.m3.p			<- YX.m3[, list(n=round(sum(w),d=1),prop= round(sum(w)/YX.m3[,sum(w)],d=3) ),by='stage']	
		#ART.pulse.N 	ART.pulse.Y 	D1<=350     D1<=550     D1>550      U           UA
		#176.3   		6.4  			41.9  		95.9 		124.7 		248.9  		73.3
		#0.230 			0.008 			0.055 		0.125 		0.162 		0.324 		0.095
		YX.m3[, col:=NULL]		
	}
	#
	#	stratify ART by 	ART.I	ART.F	ART.A	ART.P	ART.NAi 	ART.ni
	#	
	if(0)
	{	
		set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])
		#	only one t.Patient with ART.F & ART.I
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
		#	no overlap between ART.P & ART.A
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
		#	there is overlap between ART.F and ART.P, ART.A: since ART.F is largest goup, keep fewest in ART.F
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')
		#	ART and no other indications:
		set(YX.m3, YX.m3[, which(stage=='ART.started' & !is.na(ART.F) & !is.na(ART.I) & !is.na(ART.A) & !is.na(ART.P))], 'stage', 'ART.ni')	
		set(YX.m3, YX.m3[, which(stage=='ART.started' & (is.na(ART.F) | is.na(ART.I) | is.na(ART.A) | is.na(ART.P)))], 'stage', 'ART.NAi')	
		#
		YX.m3[, table(stage)]
		set(YX.m3, NULL, 'stage', YX.m3[, factor(stage)])
		tmp			<- data.table(stage= YX.m3[, levels(stage)], col=sapply( rainbow(YX.m3[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m3		<- merge(YX.m3, tmp, by='stage')
		YX.m3.fit2 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m3)				
		summary(YX.m3.fit2)
		#
		plot( seq_len(nrow(YX.m3)), YX.m3[, score.Y], pch=18, col= YX.m3[, col], cex=YX.m3[, w^0.4])
		tmp				<- subset(YX.m3, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m3.fit2, tmp, type='response'))	
		start			<- c( YX.m3.fit2$coef$mean + head( 2*sqrt(diag(vcov(YX.m3.fit2))), length(YX.m3.fit2$coef$mean)), YX.m3.fit2$coef$precision)
		YX.m3.fit2.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m3.fit2$coef$mean - head( 2*sqrt(diag(vcov(YX.m3.fit2))), length(YX.m3.fit2$coef$mean)), YX.m3.fit2$coef$precision)
		YX.m3.fit2.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
				c( predict(YX.m3.fit2.sup, tmp, type='response'), rev(predict(YX.m3.fit2.slw, tmp, type='response'))), 
				border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m3[, levels(stage)], fill=YX.m3[, unique(col)])
		YX.m3[, col:=NULL]
	}
	#
	#	stratify ART by 	ART.I	ART.F	ART.A	ART.P	-- group: ART.NAi 	ART.ni
	#	
	if(0)
	{	
		set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])
		#	only one t.Patient with ART.F & ART.I
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
		#	no overlap between ART.P & ART.A
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
		#	there is overlap between ART.F and ART.P, ART.A: since ART.F is largest goup, keep fewest in ART.F
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')
		#	ART and no other indications:
		set(YX.m3, YX.m3[, which(stage=='ART.started')], 'stage', 'ART.ni')				
		#
		set(YX.m3, NULL, 'stage', YX.m3[, factor(stage)])
		YX.m3[, table(stage)]
		#   ART.A   ART.F   ART.I  ART.ni   ART.P 	D1<=350 	D1<=550  D1>550     U      UA 
		#	7    	1609     639    5263     250    1050    	3242     3980    	7606   2566 		
		YX.m3.fit3 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m3)				
		summary(YX.m3.fit3)		
		YX.m3.fit3.or	<- data.table( 	F= my.or.from.logit(YX.m3.fit3, 'stageART.F', 'stageART.ni', subset(YX.m3, stage=='ART.F')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
				I= my.or.from.logit(YX.m3.fit3, 'stageART.I', 'stageART.ni', subset(YX.m3, stage=='ART.I')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
				P= my.or.from.logit(YX.m3.fit3, 'stageART.P', 'stageART.ni', subset(YX.m3, stage=='ART.P')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
				A= my.or.from.logit(YX.m3.fit3, 'stageART.A', 'stageART.ni', subset(YX.m3, stage=='ART.A')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
				ni= my.or.from.logit(YX.m3.fit3, 'stageART.ni', 'stageU', subset(YX.m3, stage=='ART.ni')[, sum(w)], subset(YX.m3, stage=='U')[, sum(w)], 1.962) )
		#      F        	I        P        A        ni
		#	1: 1.1514737 	1.560395 2.104169 1.348143 0.3669521
		#	2: 0.8368176 	1.132237 1.503246 1.058437 0.2985505
		#	3: 1.5844453 	2.150461 2.945311 1.717146 0.4510253	
		tmp			<- data.table(stage= YX.m3[, levels(stage)], col=sapply( rainbow(YX.m3[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m3		<- merge(YX.m3, tmp, by='stage')
		
		plot( seq_len(nrow(YX.m3)), YX.m3[, score.Y], pch=18, col= YX.m3[, col], cex=YX.m3[, w^0.4])
		tmp				<- subset(YX.m3, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m3.fit3, tmp, type='response'))	
		start			<- c( YX.m3.fit3$coef$mean + head( 2*sqrt(diag(vcov(YX.m3.fit3))), length(YX.m3.fit3$coef$mean)), YX.m3.fit3$coef$precision)
		YX.m3.fit3.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m3.fit3$coef$mean - head( 2*sqrt(diag(vcov(YX.m3.fit3))), length(YX.m3.fit3$coef$mean)), YX.m3.fit3$coef$precision)
		YX.m3.fit3.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
				c( predict(YX.m3.fit3.sup, tmp, type='response'), rev(predict(YX.m3.fit3.slw, tmp, type='response'))), 
				border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m3[, levels(stage)], fill=YX.m3[, unique(col)])
		YX.m3[, col:=NULL]
		#
		YX.m3.fit3.p	<- YX.m3[, list(n=round(sum(w),d=1),prop= round(sum(w)/YX.m3[,sum(w)],d=3) ),by='stage']	
		#ART.F   	ART.I   ART.ni  D1<=350 D1<=550 	D1>550  U       UA
		#30.7  		30.6 	110.5  	41.9  	95.9 		124.7 	248.9  	73.3
		#0.041 		0.040 	0.146 	0.055 	0.127 		0.165 	0.329 	0.097
		
	}
	#
	#	stratify ART by 	ART.pulse	ART.I	ART.F 	ART.A	ART.P	ART.ni
	#	
	if(1)
	{
		set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')
		set(YX.m3, YX.m3[, which(stage=='ART.started')], 'stage', 'ART.ni')				
		YX.m3[, table(stage)]
		# 	stage
		#	ART.A       ART.F       ART.I      ART.ni       ART.P 		ART.pulse.Y     D1<=350     D1<=550      D1>550           U          UA 
		#  	7        	1590         534        5210         250         177        	1050        3242        3980        	7606        2566
		set(YX.m3, NULL, 'stage', YX.m3[, factor(stage)])
		tmp			<- data.table(stage= YX.m3[, levels(stage)], col=sapply( rainbow(YX.m3[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m3		<- merge(YX.m3, tmp, by='stage')
		YX.m3.fit4 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m3)				
		summary(YX.m3.fit4)
		#
		YX.m3.fit4.or	<- data.table( 	F= my.or.from.logit(YX.m3.fit4, 'stageART.F', 'stageART.ni', subset(YX.m3, stage=='ART.F')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
				I= my.or.from.logit(YX.m3.fit4, 'stageART.I', 'stageART.ni', subset(YX.m3, stage=='ART.I')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
				P= my.or.from.logit(YX.m3.fit4, 'stageART.P', 'stageART.ni', subset(YX.m3, stage=='ART.P')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
				A= my.or.from.logit(YX.m3.fit4, 'stageART.A', 'stageART.ni', subset(YX.m3, stage=='ART.A')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
				pulse= my.or.from.logit(YX.m3.fit4, 'stageART.pulse.Y', 'stageART.ni', subset(YX.m3, stage=='ART.pulse.Y')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
				ni= my.or.from.logit(YX.m3.fit4, 'stageART.ni', 'stageU', subset(YX.m3, stage=='ART.ni')[, sum(w)], subset(YX.m3, stage=='U')[, sum(w)], 1.962) )	
		#	F        I        P        A    	pulse    ni
		#1: 1.160429 1.589336 2.130490 1.364599 1.707879 0.3623668
		#2: 0.841474 1.146290 1.519117 1.069798 1.215947 0.2947070
		#3: 1.600283 2.203623 2.987914 1.740636 2.398831 0.4455601
		YX.m3.fit4.p	<- YX.m3[, list(n=round(sum(w),d=1),prop= round(sum(w)/YX.m3[,sum(w)],d=3) ),by='stage']
		#	ART.A       ART.F       ART.I       ART.ni      ART.P       ART.pulse.Y D1<=350     D1<=550     D1>550      U           UA		
		#	0.9  		27.8  		22.9 		102.5   	9.2   		5.8  		38.3  		89.0 		115.6 		212.8  		66.8
		#	0.001		0.040 		0.033 		0.148 		0.013 		0.008 		0.055 		0.129 		0.167 		0.308 		0.097
		plot( seq_len(nrow(YX.m3)), YX.m3[, score.Y], pch=18, col= YX.m3[, col], cex=YX.m3[, w^0.4])
		tmp				<- subset(YX.m3, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m3.fit4, tmp, type='response'))	
		start			<- c( YX.m3.fit4$coef$mean + head( 2*sqrt(diag(vcov(YX.m3.fit4))), length(YX.m3.fit4$coef$mean)), YX.m3.fit4$coef$precision)
		YX.m3.fit4.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m3.fit4$coef$mean - head( 2*sqrt(diag(vcov(YX.m3.fit4))), length(YX.m3.fit4$coef$mean)), YX.m3.fit4$coef$precision)
		YX.m3.fit4.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
				c( predict(YX.m3.fit4.sup, tmp, type='response'), rev(predict(YX.m3.fit4.slw, tmp, type='response'))), 
				border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m3[, levels(stage)], fill=YX.m3[, unique(col)])
		YX.m3[, col:=NULL]
		#		
	}
	#	interaction of ART risk groups with VL 
	if(1)
	{
		VL.cur	<- 3
		set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])	
		YX.m3[, t.PoslRNA_T1:=NULL]
		YX.m3[, lRNA_T1:=NULL]
		tmp		<- unique(subset( clumsm.info, select=c(Patient, PoslRNA_T1) ))
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		YX.m3	<- merge(YX.m3, tmp, by= 't.Patient')				
		YX.m3	<- merge(YX.m3, YX.m3[, {	
							tmp<- which(!is.na(lRNA))
							list(lRNA_TL= ifelse(length(tmp)>0, lRNA[tail(tmp,1)], NA_real_))
						}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
		#		
		YX.m3[, lRNA.c:='su.NA']
		set(YX.m3, YX.m3[,which(stage=='ART.started' & !is.na(lRNA) & lRNA<=VL.cur)],'lRNA.c', 'ART.su.Y' )
		set(YX.m3, YX.m3[,which(stage=='ART.started' & !is.na(lRNA) & lRNA>VL.cur)],'lRNA.c', 'ART.su.N' )		
		set(YX.m3, YX.m3[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & !is.na(lRNA_TL) & lRNA_TL<=VL.cur)],'lRNA.c', 'ART.su.Y' )	
		set(YX.m3, YX.m3[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & !is.na(lRNA_TL) & lRNA_TL>VL.cur)],'lRNA.c', 'ART.su.N' )
		YX.m3[, table(lRNA.c)]
		#
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')
		set(YX.m3, YX.m3[, which(stage=='ART.started')], 'stage', 'ART.ni')				
		#
		table( YX.m3[, lRNA.c, stage] )
		#             lRNA.c
		#stage         ART.su.N ART.su.Y su.NA
		#ART.A              0        7     0
		#ART.F            181     1409     0
		#ART.I            510       19     5
		#ART.ni           327     4821    62
		#ART.P            114      133     3
		#ART.pulse.Y      108       69     0
		#D1<=350            0        0  1050
		#D1<=550            0        0  3242
		#D1>550             0        0  3980
		#U                  0        0  7606
		#UA                 0        0  2566	
		#
		set(YX.m3, YX.m3[, which(stage=='ART.pulse.Y' & lRNA.c=='ART.su.Y')], 'stage', 'ART.pulse.Y.su.Y' )
		set(YX.m3, YX.m3[, which(stage=='ART.A')], 'stage', 'DEL' )
		set(YX.m3, YX.m3[, which(stage=='ART.F' & lRNA.c=='ART.su.Y')], 'stage', 'ART.F.su.Y' )
		set(YX.m3, YX.m3[, which(stage=='ART.F' & lRNA.c=='ART.su.N')], 'stage', 'ART.F.su.N' )
		set(YX.m3, YX.m3[, which(stage=='ART.F' & lRNA.c=='su.NA')], 'stage', 'DEL' )
		set(YX.m3, YX.m3[, which(stage=='ART.I' & lRNA.c=='ART.su.Y')], 'stage', 'ART.I.su.Y' )
		set(YX.m3, YX.m3[, which(stage=='ART.I' & lRNA.c=='ART.su.N')], 'stage', 'ART.I.su.N' )
		set(YX.m3, YX.m3[, which(stage=='ART.I' & lRNA.c=='su.NA')], 'stage', 'DEL' )
		set(YX.m3, YX.m3[, which(stage=='ART.ni' & lRNA.c=='ART.su.Y')], 'stage', 'ART.ni.su.Y' )
		set(YX.m3, YX.m3[, which(stage=='ART.ni' & lRNA.c=='ART.su.N')], 'stage', 'ART.ni.su.N' )
		set(YX.m3, YX.m3[, which(stage=='ART.ni' & lRNA.c=='su.NA')], 'stage', 'DEL' )
		set(YX.m3, YX.m3[, which(stage=='ART.P' & lRNA.c=='ART.su.Y')], 'stage', 'ART.P.su.Y' )
		set(YX.m3, YX.m3[, which(stage=='ART.P' & lRNA.c=='ART.su.N')], 'stage', 'ART.P.su.N' )
		set(YX.m3, YX.m3[, which(stage=='ART.P' & lRNA.c=='su.NA')], 'stage', 'DEL' )
		set(YX.m3, YX.m3[, which(stage=='ART.pulse.Y' & lRNA.c=='ART.su.Y')], 'stage', 'ART.pulse.su.Y' )
		set(YX.m3, YX.m3[, which(stage=='ART.pulse.Y' & lRNA.c=='ART.su.N')], 'stage', 'ART.pulse.su.N' )
		set(YX.m3, YX.m3[, which(stage=='ART.pulse.Y' & lRNA.c=='su.NA')], 'stage', 'DEL' )
		YX.m3[, table(stage)]
		#
		YX.m3.s	<- subset(YX.m3, stage!='DEL')
		set(YX.m3.s, NULL, 'stage', YX.m3.s[, factor(stage)])
		tmp			<- data.table(stage= YX.m3.s[, levels(stage)], col=sapply( rainbow(YX.m3.s[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m3.s		<- merge(YX.m3.s, tmp, by='stage')
		YX.m3.fit5 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m3.s)				
		summary(YX.m3.fit5)
		#
		YX.m3.fit5.p	<- YX.m3.s[, list(n=round(sum(w),d=1),prop= round(sum(w)/YX.m3.s[,sum(w)],d=3) ),by='stage']
		#	ART.F.su.N      ART.F.su.Y      ART.I.su.N      ART.I.su.Y      ART.ni.su.N      ART.ni.su.Y      ART.P.su.N       ART.P.su.Y    ART.pulse.su.N   ART.pulse.Y.su.Y  D1<=350         D1<=550         D1>550          U       UA  		
		#	5.8  			22.0  			21.9   			0.8  			11.4  			 88.7   		  4.7   		   4.3   		 4.2   			  1.6  				38.3  			89.0 			115.6 			212.8  	66.8				
		#
		#	these numbers are too small to yield robust coefficients
		#
		plot( seq_len(nrow(YX.m3.s)), YX.m3.s[, score.Y], pch=18, col= YX.m3.s[, col], cex=YX.m3.s[, w^0.4])
		tmp				<- subset(YX.m3.s, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m3.fit5, tmp, type='response'))	
		start			<- c( YX.m3.fit5$coef$mean + head( 2*sqrt(diag(vcov(YX.m3.fit5))), length(YX.m3.fit5$coef$mean)), YX.m3.fit5$coef$precision)
		YX.m3.fit5.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3.s, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m3.fit5$coef$mean - head( 2*sqrt(diag(vcov(YX.m3.fit5))), length(YX.m3.fit5$coef$mean)), YX.m3.fit5$coef$precision)
		YX.m3.fit5.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3.s, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
				c( predict(YX.m3.fit5.sup, tmp, type='response'), rev(predict(YX.m3.fit5.slw, tmp, type='response'))), 
				border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m3.s[, levels(stage)], fill=YX.m3.s[, unique(col)])
		YX.m3.s[, col:=NULL]
		#
	}
	#	interaction of ART risk groups with VL 
	if(1)
	{
		VL.cur	<- log10(1e3)
		VL.cur	<- log10(4750)
		set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])	
		YX.m3[, t.PoslRNA_T1:=NULL]
		YX.m3[, lRNA_T1:=NULL]
		YX.m3[, lRNA.c:=NULL]
		tmp		<- unique(subset( clumsm.info, select=c(Patient, PoslRNA_T1) ))
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		YX.m3	<- merge(YX.m3, tmp, by= 't.Patient')				
		YX.m3	<- merge(YX.m3, YX.m3[, {	
							tmp<- which(!is.na(lRNA))
							list(lRNA_TL= ifelse(length(tmp)>0, lRNA[tail(tmp,1)], NA_real_))
						}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
		#		
		YX.m3[, lRNA.c:='su.NA']
		set(YX.m3, YX.m3[,which(stage=='ART.started' & !is.na(lRNA) & lRNA<=VL.cur)],'lRNA.c', 'ART.su.Y' )
		set(YX.m3, YX.m3[,which(stage=='ART.started' & !is.na(lRNA) & lRNA>VL.cur)],'lRNA.c', 'ART.su.N' )		
		set(YX.m3, YX.m3[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & !is.na(lRNA_TL) & lRNA_TL<=VL.cur)],'lRNA.c', 'ART.su.Y' )	
		set(YX.m3, YX.m3[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & !is.na(lRNA_TL) & lRNA_TL>VL.cur)],'lRNA.c', 'ART.su.N' )
		#
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')
		set(YX.m3, YX.m3[, which(stage=='ART.started')], 'stage', 'ART.ni')				
		#
		table( YX.m3[, lRNA.c, stage] )
		#             lRNA.c
		#stage         ART.su.N ART.su.Y su.NA
		#ART.A              0        7     0
		#ART.F            181     1409     0
		#ART.I            510       19     5
		#ART.ni           327     4821    62
		#ART.P            114      133     3
		#ART.pulse.Y      108       69     0
		#D1<=350            0        0  1050
		#D1<=550            0        0  3242
		#D1>550             0        0  3980
		#U                  0        0  7606
		#UA                 0        0  2566	
		#
		YX.m3.s			<- subset(YX.m3, stage%in%c('ART.A','ART.F','ART.I','ART.P','ART.ni','ART.pulse.Y'))
		set(YX.m3.s, NULL, 'stage', YX.m3.s[, factor(stage)])
		set(YX.m3.s, NULL, 'lRNA.c', YX.m3.s[, factor(lRNA.c)])
		
		YX.m3.s.fit6 	<- betareg(score.Y ~ lRNA.c + stage - 1 , link='logit', weights=w, data = YX.m3.s)
		#
		#	risk groups explain more of the variation than VL --> in line with  YX.m3.s.fit5. Likely an artifact of the method!
		#	BUT shows that most infection window contain periods of high and low VL
		#
		tmp			<- data.table(stage= YX.m3.s[, levels(stage)], col=sapply( rainbow(YX.m3.s[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m3.s		<- merge(YX.m3.s, tmp, by='stage')
		plot( seq_len(nrow(YX.m3.s)), YX.m3.s[, score.Y], pch=18, col= YX.m3.s[, col], cex=YX.m3.s[, w^0.4])
		tmp				<- subset(YX.m3.s, select=c(stage, lRNA.c))
		set(tmp, tmp[, which(lRNA.c!='ART.su.Y')], 'lRNA.c', 'ART.su.Y')		
		lines(seq_len(nrow(tmp)), predict(YX.m3.s.fit6, tmp, type='response'))	
		tmp				<- subset(YX.m3.s, select=c(stage, lRNA.c))
		set(tmp, tmp[, which(lRNA.c!='ART.su.N')], 'lRNA.c', 'ART.su.N')		
		lines(seq_len(nrow(tmp)), predict(YX.m3.s.fit6, tmp, type='response'), col='red')	
		legend('bottomright', bty='n', border=NA, legend= YX.m3.s[, levels(stage)], fill=YX.m3.s[, unique(col)])
		YX.m3.s[, col:=NULL]		
	}
	#	interaction of ART risk groups with VL -- see if linear model on VL makes any difference to the conclusion above
	if(1)
	{
		set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])	
		YX.m3[, t.PoslRNA_T1:=NULL]
		YX.m3[, lRNA_TL:=NULL]
		YX.m3[, lRNA.c:=NULL]
		tmp		<- unique(subset( clumsm.info, select=c(Patient, PoslRNA_T1) ))
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		YX.m3	<- merge(YX.m3, tmp, by= 't.Patient')				
		YX.m3	<- merge(YX.m3, YX.m3[, {	
							tmp<- which(!is.na(lRNA))
							list(lRNA_TL= ifelse(length(tmp)>0, lRNA[tail(tmp,1)], NA_real_))
						}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
		#
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')
		set(YX.m3, YX.m3[, which(stage=='ART.started')], 'stage', 'ART.ni')				
		#
		YX.m3.s			<- subset(YX.m3, stage%in%c('ART.A','ART.F','ART.I','ART.P','ART.ni','ART.pulse.Y'))
		set(YX.m3.s, NULL, 'stage', YX.m3.s[, factor(stage)])
		
		YX.m3.s.fit7 	<- betareg(score.Y ~ lRNA + stage, link='logit', weights=w, data = YX.m3.s)
		#
		#
		tmp			<- data.table(stage= YX.m3.s[, levels(stage)], col=sapply( rainbow(YX.m3.s[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m3.s		<- merge(YX.m3.s, tmp, by='stage')
		plot( seq_len(nrow(YX.m3.s)), YX.m3.s[, score.Y], pch=18, col= YX.m3.s[, col], cex=YX.m3.s[, w^0.4])
		tmp				<- subset(YX.m3.s, select=c(stage, lRNA))
		set(tmp, NULL, 'lRNA', 3)		
		lines(seq_len(nrow(tmp)), predict(YX.m3.s.fit7, tmp, type='response'))	
		tmp				<- subset(YX.m3.s, select=c(stage, lRNA))
		set(tmp, NULL, 'lRNA', 5)				
		lines(seq_len(nrow(tmp)), predict(YX.m3.s.fit7, tmp, type='response'), col='red')	
		legend('bottomright', bty='n', border=NA, legend= YX.m3.s[, levels(stage)], fill=YX.m3.s[, unique(col)])
		YX.m3.s[, col:=NULL]		
	}
	#	interaction of ART risk groups with VL -- ever high VL
	if(1)
	{
		VL.cur	<- log10(1e3)
		#VL.cur	<- log10(4750)
		set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])	
		YX.m3[, t.PoslRNA_T1:=NULL]
		YX.m3[, lRNA_TL:=NULL]
		YX.m3[, lRNA.c:=NULL]
		
		tmp		<- YX.m3[, {
					tmp<- which(!is.na(lRNA))
					list( lRNA.c= ifelse(length(tmp), ifelse( max(lRNA[tmp])>VL.cur, 'SuA.N', 'SuA.Y'), 'SuA.NA') )
				}, by=c('Patient','t.Patient')]
		YX.m3	<- merge(YX.m3, tmp, by=c('Patient','t.Patient'))		
		#
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')
		set(YX.m3, YX.m3[, which(stage=='ART.started')], 'stage', 'ART.ni')				
		table( YX.m3[, lRNA.c, stage] )
		#	
		YX.m3.s			<- subset(YX.m3, stage%in%c('ART.A','ART.F','ART.I','ART.P','ART.ni','ART.pulse.Y'))
		set(YX.m3.s, NULL, 'stage', YX.m3.s[, factor(stage)])
		set(YX.m3.s, NULL, 'lRNA.c', YX.m3.s[, factor(lRNA.c)])
		table( YX.m3.s[, lRNA.c, stage] )
		
		YX.m3.s.fit8 	<- betareg(score.Y ~ lRNA.c + stage - 1, link='logit', weights=w, data = YX.m3.s)
		#
		tmp				<- data.table(stage= YX.m3.s[, levels(stage)], col=sapply( rainbow(YX.m3.s[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m3.s			<- merge(YX.m3.s, tmp, by='stage')
		plot( seq_len(nrow(YX.m3.s)), YX.m3.s[, score.Y], pch=18, col= YX.m3.s[, col], cex=YX.m3.s[, w^0.4])
		tmp				<- subset(YX.m3.s, select=c(stage, lRNA.c))
		set(tmp, tmp[, which(lRNA.c!='SuA.N')], 'lRNA.c', 'SuA.N')
		lines(seq_len(nrow(tmp)), predict(YX.m3.s.fit8, tmp, type='response'))	
		tmp				<- subset(YX.m3.s, select=c(stage, lRNA.c))
		set(tmp, tmp[, which(lRNA.c!='SuA.Y')], 'lRNA.c', 'SuA.Y')				
		lines(seq_len(nrow(tmp)), predict(YX.m3.s.fit8, tmp, type='response'), col='red')	
		legend('bottomright', bty='n', border=NA, legend= YX.m3.s[, levels(stage)], fill=YX.m3.s[, unique(col)])
		YX.m3.s[, col:=NULL]
		#
		#	comparison of the 3 models with VL
		#
		YX.m3.vlcmp	<- c( YX.m3.s.fit6$pseudo.r.squared, YX.m3.s.fit7$pseudo.r.squared, YX.m3.s.fit8$pseudo.r.squared)
		#	0.006771914 	0.005305756 	0.013678193
		AIC( YX.m3.s.fit6, YX.m3.s.fit7, YX.m3.s.fit8 )	
		#	             df       AIC
		#YX.m3.s.fit6  9 -217.0417
		#YX.m3.s.fit7  8 -211.3867
		#YX.m3.s.fit8  9 -219.4119
		YX.m3.fit8.or	<- data.table( 	F= my.or.from.logit(YX.m3.s.fit8, 'stageART.F', 'stageART.ni', subset(YX.m3.s, stage=='ART.F')[, sum(w)], subset(YX.m3.s, stage=='ART.ni')[, sum(w)], 1.962),
				I= my.or.from.logit(YX.m3.s.fit8, 'stageART.I', 'stageART.ni', subset(YX.m3.s, stage=='ART.I')[, sum(w)], subset(YX.m3.s, stage=='ART.ni')[, sum(w)], 1.962),
				P= my.or.from.logit(YX.m3.s.fit8, 'stageART.P', 'stageART.ni', subset(YX.m3.s, stage=='ART.P')[, sum(w)], subset(YX.m3.s, stage=='ART.ni')[, sum(w)], 1.962),
				A= my.or.from.logit(YX.m3.s.fit8, 'stageART.A', 'stageART.ni', subset(YX.m3.s, stage=='ART.A')[, sum(w)], subset(YX.m3.s, stage=='ART.ni')[, sum(w)], 1.962),
				pulse= my.or.from.logit(YX.m3.s.fit8, 'stageART.pulse.Y', 'stageART.ni', subset(YX.m3.s, stage=='ART.pulse.Y')[, sum(w)], subset(YX.m3.s, stage=='ART.ni')[, sum(w)], 1.962) )
		#YX.m3.fit8.or
		#			F        	I         	P  			A     	pulse
		#1: 		1.0631908 	1.250178 	1.8097541 	NA 		1.4108150
		#2: 		0.5908315 	0.661459 	0.7035367 	NA 		0.4358457
		#3: 		1.9131931 	2.362876 	4.6553505 	NA 		4.5667518						
	}
}
######################################################################################
project.athena.Fisheretal.YX.model2.time<- function(YX, clumsm.info, df.viro, cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), plot.file=NA, score.Y.cut=1e-2)
{
	require(betareg)
	#	Treatment cascade by time period
	#	prop and odds by updated treatment cascade
	YX.m2	<- copy(YX)
	YX.m2[, U.score:=NULL]
	set(YX.m2, YX.m2[,which(score.Y>1-score.Y.cut)], 'score.Y', 1-score.Y.cut)
	set(YX.m2, YX.m2[,which(score.Y<score.Y.cut)], 'score.Y', score.Y.cut)
	#	pool acute given small numbers
	set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'Ay' )
	set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'Am' )
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Yes')], 'stage', 'Ay')
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Maybe')], 'stage', 'Am')
	#	stratify Diag by first CD4 after diagnosis
	tmp			<- subset(clumsm.info, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
	setkey(tmp, Patient)
	tmp			<- unique(tmp)	
	tmp[, CD4.c:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
	cat(paste('\ntransmitters with no CD4 after diagnosis\n'))
	print( subset( tmp, is.na(CD4.c) ) )
	setnames(tmp, 'Patient', 't.Patient')
	set(tmp, tmp[, which(is.na(CD4.c))], 'CD4.c', 'D1.NA')
	YX.m2	<- merge(YX.m2, subset(tmp, select=c(t.Patient, CD4.c)), by='t.Patient')
	tmp		<- YX.m2[, which(stage=='Diag')]
	set(YX.m2, tmp, 'stage', YX.m2[tmp, CD4.c])
	YX.m2[, CD4.c:=NULL]
	YX.m2[, table(stage)]
	
	#set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute%in%acute.select)], 'stage', 'UA')
	YX.m2[, stage.orig:= stage]
	if(0)
	{
		# collect PoslRNA_T1
		tmp		<- unique(subset( clumsm.info, select=c(Patient, PoslRNA_T1) ))
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		YX.m2	<- merge(YX.m2, tmp, by= 't.Patient')	
		# collect PoslRNA_TL and lRNA_TL
		tmp		<- merge( unique(subset(clumsm.info, select=Patient)), subset(df.viro, select=c(Patient, PosRNA, lRNA)), by='Patient' )
		set( tmp, NULL, 'PosRNA', hivc.db.Date2numeric(tmp[,PosRNA]) )
		tmp		<- subset( tmp[, 	{
							tmp<- which(!is.na(lRNA))
							list(t.Patient=Patient[1], lRNA_TL=lRNA[tail(tmp,1)], PoslRNA_TL=PosRNA[tail(tmp,1)])
						}, by='Patient'], select=c(t.Patient, lRNA_TL, PoslRNA_TL) )
		YX.m2	<- merge(YX.m2, tmp, by= 't.Patient', all.x=TRUE)
		set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
		tmp		<- YX.m2[, {
					tmp<- which(!is.na(lRNA))
					if(!length(tmp) & any(t>t.PoslRNA_T1[1]) & any(contact=='Yes') & !is.na(lRNA_TL[1]) & lRNA_TL[1]<=vl.suppressed)
						lRNA.c	<- 'SuA.Y'						
					else if(!length(tmp) & any(t>t.PoslRNA_T1[1]) & any(contact=='Yes') & !is.na(lRNA_TL[1]) & lRNA_TL[1]>vl.suppressed)
						lRNA.c	<- 'SuA.N'						
					else if(!length(tmp))
						lRNA.c	<- 'SuA.NA'
					else if(max(lRNA[tmp])>vl.suppressed)
						lRNA.c	<- 'SuA.N'
					else
						lRNA.c	<- 'SuA.Y'
					list( lRNA.c= lRNA.c )
				}, by=c('Patient','t.Patient')]	
		YX.m2	<- merge(YX.m2, tmp, by=c('Patient','t.Patient'))	
		set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.c=='SuA.Y')],'stage', 'ART.suA.Y' )
		set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.c=='SuA.N')],'stage', 'ART.suA.N' )
		set(YX.m2, YX.m2[,which(stage=='ART.started')],'stage', 'ART.vlNA' )
		YX.m2		<- subset(YX.m2, select=c(t, t.Patient, Patient, cluster, score.Y, stage, t.period, CDCC, lRNA, t.Age, t.isAcute, t.PoslRNA_T1, lRNA_TL, PoslRNA_TL, t.AnyT_T1, contact, fw.up.med, w, w.t, stage.orig  ))
		
		#stage	(  using any(contact==='yes')  )
		#   Am 		ART.suA.N 	ART.suA.Y  	ART.vlNA        Ay   		D1<=350   	D1<=550    	D1>550         U 
		# 	2522     3187      	4532        49      		2111       	865      	2836      	3612     	10528		
	}
	if(1)
	{
		YX.m2	<- merge(YX.m2, YX.m2[, {
											tmp		<- which(!is.na(lRNA))
											list( lRNA.c= ifelse(!length(tmp), 'SuA.NA', ifelse(max(lRNA[tmp])>vl.suppressed, 'SuA.N', 'SuA.Y')) )
										}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))				
		set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.c=='SuA.Y')],'stage', 'ART.suA.Y' )
		set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.c=='SuA.N')],'stage', 'ART.suA.N' )
		set(YX.m2, YX.m2[,which(stage=='ART.started')],'stage', 'ART.vlNA' )
		YX.m2		<- subset(YX.m2, select=c(t, t.Patient, Patient, cluster, score.Y, stage, t.period, CDCC, lRNA, t.Age, t.isAcute, t.AnyT_T1, contact, fw.up.med, w, w.t, stage.orig  ))		
		#stage
		#   Am ART.suA.N ART.suA.Y  ART.vlNA        Ay   D1<=350   D1<=550    D1>550         U 
		# 2522      3172      4526        70      2111       865      2836      3612     10528 		
	}
	if(0)
	{		
		ggplot(subset(YX.m2, stage=='ART.suA.Y'), aes(t.period, score.Y)) + geom_boxplot() + geom_jitter(size=0.4)
		ggplot(subset(YX.m2, stage=='ART.suA.Y'), aes(t.period, w)) + geom_boxplot() + geom_jitter(size=0.4) + scale_y_continuous(breaks = seq(0,1,0.05))
		#ggplot(YX.m2, aes(t.period, w)) + geom_boxplot()  + scale_y_continuous(breaks = seq(0,1,0.05))
		ggplot(YX.m2, aes(t.period, score.Y)) + geom_boxplot() + geom_jitter(size=0.4)
	}	
	#subset(YX.m2, stage=='ART.suA.Y')[, list(n=length(w), w=sum(w), wpn= mean(w)),by='t.period']	
	#YX.m2[, table(stage)]
	#subset( YX.m2, stage=='ART.vlNA' ) 	
	YX.m2.p		<- YX.m2[, 	list(n=round(sum(w),d=1)) , by= c('stage','t.period')]
	YX.m2.p		<- merge( YX.m2.p, YX.m2[, list(N= sum(w)), by='t.period'], by='t.period' )
	
	tmp			<- unique(subset(clumsm.info, select=c(Patient, AnyPos_T1)))
	tmp			<- merge( unique(subset(YX.m2, select=c(Patient, t.period))), tmp, by='Patient' )
	YX.m2.p		<- merge( YX.m2.p, tmp[, list(t.min=min(AnyPos_T1), t.max=max(AnyPos_T1)), by='t.period'], by='t.period' )
	YX.m2.p[, p:= n/N]
	setkey(YX.m2.p, stage)
	#	get Binomial confidence intervals
	require(Hmisc)
	YX.m2.p		<- merge( YX.m2.p, YX.m2.p[ , 	{
													tmp<- binconf( n, N, alpha=0.05, method= "wilson", include.x=FALSE, include.n=FALSE, return.df=TRUE)
													list(p.l95=tmp$Lower, p.u95=tmp$Upper)
												}, by=c('stage','t.period')], by=c('stage','t.period'))
	#	plot
	if(!is.na(plot.file))
	{
		pdf(file=plot.file, w=5, h=5)
		tmp			<- data.table( stage=YX.m2.p[, unique(stage)], h=seq(1, YX.m2.p[,length(unique(stage))]), col= rainbow(YX.m2.p[,length(unique(stage))] ) )
		set(tmp, NULL, 'h', tmp[,(h-mean(h))/50] )
		YX.m2.p		<- merge(YX.m2.p, tmp, by='stage')
		set(YX.m2.p, NULL, 't.period', YX.m2.p[,as.numeric(t.period)])
		plot(1,1,type='n',xlim=range(YX.m2.p[,t.period]), ylim=range(unlist( subset( YX.m2.p, select=c(p, p.l95, p.u95) ) )),xlab='time period',ylab='prop')
		YX.m2.p[, lines(t.period+h, p, col=col, type='b') ,by='stage']
		YX.m2.p[, lines(rep(t.period+h,2), c(p.l95,p.u95), col=col ) ,by=c('stage','t.period')]	
		legend('topright', bty='n',border=NA,fill=YX.m2.p[, unique(col)], legend=YX.m2.p[, unique(stage)] )	
		YX.m2.p[, col:=NULL]	
		YX.m2.p[, h:=NULL]	
		dev.off()
	}
	#	
	YX.m2.p
}
######################################################################################
project.athena.Fisheretal.YX.model2<- function(YX, clumsm.info, df.viro, vl.suppressed=log10(1e3), acute.select=c('Yes','Maybe'), cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), score.Y.cut=1e-2)
{
	require(betareg)
	#	Treatment cascade:
	#	prop and odds by traditional categories 	undiagnosed, 	undiagnosed&evidence acute, 	diagnosed by CD4 at diagnosis, 	 ART by first suppressed
	#
	#	PQ: what is the preventative effect of reaching the endpoint of the treatment cascade?
	#	PQ: single VL that maximizes difference between ART/firstsuppressed and ART/notyetsuppressed
	#	SQ: does acute matter for diagnosed?
	#	SQ: does CDC matter?
	#	SQ: does VL at diagnosis or age at diagnosis matter?	
	YX.m2	<- copy(YX)
	YX.m2[, U.score:=NULL]
	set(YX.m2, YX.m2[,which(score.Y>1-score.Y.cut)], 'score.Y', 1-score.Y.cut)
	set(YX.m2, YX.m2[,which(score.Y<score.Y.cut)], 'score.Y', score.Y.cut)
	
	#	diagnosed and acute
	set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.5)], 'stage', 'DAy' )
	set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.5)], 'stage', 'DAm' )		
	#	stratify Diag by first CD4 after diagnosis
	cd4.label	<- c('D1<=350','D1<=550','D1>550')
	tmp			<- subset(clumsm.info, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
	setkey(tmp, Patient)
	tmp			<- unique(tmp)	
	tmp[, CD4.c:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
	cat(paste('\ntransmitters with no CD4 after diagnosis\n'))
	print( subset( tmp, is.na(CD4.c) ) )
	setnames(tmp, 'Patient', 't.Patient')
	set(tmp, tmp[, which(is.na(CD4.c))], 'CD4.c', 'D1.NA')
	YX.m2	<- merge(YX.m2, subset(tmp, select=c(t.Patient, CD4.c)), by='t.Patient')
	tmp		<- YX.m2[, which(stage=='Diag')]
	set(YX.m2, tmp, 'stage', YX.m2[tmp, CD4.c])
	YX.m2[, CD4.c:=NULL]
	YX.m2[, table(stage)]
	#	if transmitter acute, set undiagnosed to undiagnosed and acute
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Yes')], 'stage', 'UAy')
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Maybe')], 'stage', 'UAm')
	#set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute%in%acute.select)], 'stage', 'UA')
	YX.m2[, stage.orig:= stage]
	#
	tmp		<- unique(subset( clumsm.info, select=c(Patient, PoslRNA_T1) ))
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	YX.m2	<- merge(YX.m2, tmp, by= 't.Patient')				
	YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, cluster, score.Y, stage, CDCC, lRNA, t.Age, t.isAcute, t.PoslRNA_T1, t.AnyT_T1, contact, fw.up.med, w, stage.orig  ))	
	#
	#	VL first suppressed (treatment cascade endpoint)
	#
	if(0)
	{
		#	VL first suppressed threshold: trial and error	
		set(YX.m2, NULL, 'stage', YX.m2[, factor(stage)])
		#		consider different VL thresholds - no effect
		VL.win		<- c( seq(400, 1000, by=100), seq(1250, 5000, by=250), 1e4 ) 
		YX.m2.VL	<- sapply(VL.win, function(VL.cur)
				{
					cat(paste('\nprocess VL=', VL.cur))
					VL.cur	<- log10(VL.cur)
					set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
					
					tmp		<- subset(df.viro, select=c(Patient, PosRNA, lRNA))
					setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
					set(tmp, NULL, 't.PosRNA', hivc.db.Date2numeric(tmp[,t.PosRNA]))
					tmp		<- merge( unique(subset(YX.m2, select=c(t.Patient, t.AnyT_T1))), tmp, by='t.Patient' )
					tmp		<- subset(tmp, t.AnyT_T1<=t.PosRNA)[, {
								z<- which(t.lRNA<VL.cur)
								list(t.ARTSu_T1= ifelse(length(z), t.PosRNA[z[1]], NA_real_) )	
							}, by='t.Patient']
					YX.m2	<- merge(YX.m2, tmp, by= 't.Patient')
					set(YX.m2, YX.m2[, which(stage=='ART.started' & !is.na(t.ARTSu_T1) & t>=t.ARTSu_T1)], 'stage', 'ART1.su.Y')
					set(YX.m2, YX.m2[, which(stage=='ART.started')], 'stage', 'ART1.su.N')
					set(YX.m2, NULL, 'stage', YX.m2[, factor(as.character(stage))])
					YX.m2[, t.ARTSu_T1:=NULL]
					
					YX.m2.fit1 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m2)				
					ans	<- c(		coef(YX.m2.fit1)['stageART1.su.Y'], coef(YX.m2.fit1)['stageART1.su.N'], coef(YX.m2.fit1)['stageU'],	
							sqrt(diag(vcov(YX.m2.fit1)))[c('stageART1.su.Y','stageART1.su.N')],
							my.or.from.logit(YX.m2.fit1, 'stageART1.su.Y', 'stageART1.su.N', subset(YX.m2, stage=='ART1.su.Y')[, sum(w)], subset(YX.m2, stage=='ART1.su.N')[, sum(w)], 1.962),						
							my.or.from.logit(YX.m2.fit1, 'stageART1.su.Y', 'stageU', subset(YX.m2, stage=='ART1.su.Y')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
							my.or.from.logit(YX.m2.fit1, 'stageART1.su.N', 'stageU', subset(YX.m2, stage=='ART1.su.N')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
							logLik(YX.m2.fit1), YX.m2.fit1$pseudo.r.squared, 10^VL.cur	)
					names(ans)	<- c('coef.su.Y','coef.su.N','coef.U','coef.su.Y.sd','coef.su.N.sd','or.YN','or.YN.l95','or.YN.u95','or.YU','or.YU.l95','or.YU.u95','or.YN','or.YN.l95','or.YN.u95','lkl','r2','VL.thresh')
					ans
				})
		YX.m2.VL	<- as.data.table(t(YX.m2.VL))
		plot(YX.m2.VL[, VL.thresh], YX.m2.VL[, or.YN], type='p', pch=18, ylim= c( min(YX.m2.VL[, or.YN.l95], na.rm=TRUE),max(YX.m2.VL[, or.YN.u95], na.rm=TRUE))  )
		dummy		<- YX.m2.VL[,	lines( rep(VL.thresh,2), c(or.YN.l95,or.YN.u95))	, by='VL.thresh']
		#
		#		plot VL thresh 1e3
		#
		VL.cur		<- log10(1e3)
		set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
		tmp			<- subset(df.viro, select=c(Patient, PosRNA, lRNA))
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		set(tmp, NULL, 't.PosRNA', hivc.db.Date2numeric(tmp[,t.PosRNA]))
		tmp			<- merge( unique(subset(YX.m2, select=c(t.Patient, t.AnyT_T1))), tmp, by='t.Patient' )
		tmp			<- subset(tmp, t.AnyT_T1<=t.PosRNA)[, 	{
					z<- which(t.lRNA<VL.cur)
					list(t.ARTSu_T1= ifelse(length(z), t.PosRNA[z[1]], NA_real_) )	
				}, by='t.Patient']
		YX.m2		<- merge(YX.m2, tmp, by= 't.Patient')
		set(YX.m2, YX.m2[, which(stage=='ART.started' & !is.na(t.ARTSu_T1) & t>=t.ARTSu_T1)], 'stage', 'ART1.su.Y')
		set(YX.m2, YX.m2[, which(stage=='ART.started')], 'stage', 'ART1.su.N')		
		set(YX.m2, NULL, 'stage', YX.m2[, factor(as.character(stage))])
		YX.m2[, t.ARTSu_T1:=NULL]
		YX.m2.fit1 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m2)				
		#	plot
		tmp			<- data.table(stage= YX.m2[, levels(stage)], col=sapply( rainbow(YX.m2[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m2		<- merge(YX.m2, tmp, by='stage')
		plot( seq_len(nrow(YX.m2)), YX.m2[, score.Y], pch=18, col= YX.m2[, col], cex=YX.m2[, w^0.4])
		tmp				<- subset(YX.m2, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m2.fit1, tmp, type='response'))	
		start			<- c( YX.m2.fit1$coef$mean + head( 2*sqrt(diag(vcov(YX.m2.fit1))), length(YX.m2.fit1$coef$mean)), YX.m2.fit1$coef$precision)
		YX.m2.fit1.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m2, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m2.fit1$coef$mean - head( 2*sqrt(diag(vcov(YX.m2.fit1))), length(YX.m2.fit1$coef$mean)), YX.m2.fit1$coef$precision)
		YX.m2.fit1.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m2, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
				c( predict(YX.m2.fit1.sup, tmp, type='response'), rev(predict(YX.m2.fit1.slw, tmp, type='response'))), 
				border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m2[, levels(stage)], fill=YX.m2[, unique(col)])
		YX.m2[, col:=NULL]
		#
		YX.m2.fit1.or	<- data.table( 	ART1.su.NY= my.or.from.logit(YX.m2.fit1, 'stageART1.su.N', 'stageART1.su.Y', subset(YX.m2, stage=='ART1.su.N')[, sum(w)], subset(YX.m2, stage=='ART1.su.Y')[, sum(w)], 1.962),
				ART1.su.N= my.or.from.logit(YX.m2.fit1, 'stageART1.su.N', 'stageU', subset(YX.m2, stage=='ART1.su.N')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				ART1.su.Y= my.or.from.logit(YX.m2.fit1, 'stageART1.su.Y', 'stageU', subset(YX.m2, stage=='ART1.su.Y')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				DAy= my.or.from.logit(YX.m2.fit1, 'stageDAy', 'stageU', subset(YX.m2, stage=='DAy')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				DAm= my.or.from.logit(YX.m2.fit1, 'stageDAm', 'stageU', subset(YX.m2, stage=='DAm')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				D1l350= my.or.from.logit(YX.m2.fit1, 'stageD1<=350', 'stageU', subset(YX.m2, stage=='D1<=350')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				D1l550= my.or.from.logit(YX.m2.fit1, 'stageD1<=550', 'stageU', subset(YX.m2, stage=='D1<=550')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				D1g550= my.or.from.logit(YX.m2.fit1, 'stageD1>550', 'stageU', subset(YX.m2, stage=='D1>550')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				UAm= my.or.from.logit(YX.m2.fit1, 'stageUAm', 'stageU', subset(YX.m2, stage=='UAm')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				UAy= my.or.from.logit(YX.m2.fit1, 'stageUAy', 'stageU', subset(YX.m2, stage=='UAy')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962) 	)
		# YX.m2.fit1.or
		#   ART1.su.NY ART1.su.N ART1.su.Y      DAy      DAm    D1l350    D1l550    D1g550      UAm      UAy
		#1:  1.0819401 0.4959514 0.4583908 1.478592 2.080791 0.7559947 0.5594259 0.6672697 1.977886 2.746282
		#2:  0.8151495 0.3866189 0.3692749 1.139423 1.610585 0.5824594 0.4373290 0.5251420 1.538850 2.142883
		#3:  1.4360489 0.6362022 0.5690127 1.918720 2.688273 0.9812323 0.7156107 0.8478636 2.542179 3.519587
		YX.m2.p		<- YX.m2[, list(n=round(sum(w),d=1),prop= round(sum(w)/YX.m2[,sum(w)],d=3) ),by='stage']
		#	ART1.su.N 	ART1.su.Y 	D1<=350 D1<=550 D1>550 	DAm 	DAy 	U 		UAm 	UAy
		#	48.5 		120.5  		26.1  	53.7  	65.9  	13.0  	12.2 	168.6  	25.8  24.5
		#	0.087 		0.216 		0.047 	0.096 	0.118 	0.023 	0.022 	0.302 	0.046 0.044
		tmp			<- cooks.distance(YX.m2.fit1)
		plot(seq_along(tmp), tmp, type='h', col=YX.m2[, col], xlab='index', ylab='Cooks D')
		legend('topright', bty='n', border=NA, legend= YX.m2[, levels(stage)], fill=YX.m2[, unique(col)])
		tmp			<- residuals(YX.m2.fit1)
		plot(seq_along(tmp), tmp, type='p', pch=18, col=YX.m2[, col], xlab='index', ylab='std residuals')
		legend('topright', bty='n', border=NA, legend= YX.m2[, levels(stage)], fill=YX.m2[, unique(col)])		
	}
	#
	#	VL suppressed (more rigorous endpoint under continued monitoring)
	#
	if(1)
	{		
		#	split missing VL into contact Yes/No and before first contact
		#VL.cur	<- log10(VL.cur)
		set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
		set(YX.m2, YX.m2[,which(stage=='ART.started' & !is.na(lRNA) & lRNA<=VL.cur)],'stage', 'ART.su.Y' )
		set(YX.m2, YX.m2[,which(stage=='ART.started' & !is.na(lRNA) & lRNA>VL.cur)],'stage', 'ART.su.N' )
		set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='No')],'stage', 'ART.vlNA.c.N' )
		set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes')],'stage', 'ART.vlNA.c.Y' )
		set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t<t.PoslRNA_T1 )],'stage', 'ART.vlNA.c.b4T1' )					
		set(YX.m2, NULL, 'stage', YX.m2[, factor(as.character(stage))])
		YX.m2[, table(stage)]
		
		YX.m2.fit2 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m2)
		#	plot
		tmp			<- data.table(stage= YX.m2[, levels(stage)], col=sapply( rainbow(YX.m2[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m2		<- merge(YX.m2, tmp, by='stage')
		plot( seq_len(nrow(YX.m2)), YX.m2[, score.Y], pch=18, col= YX.m2[, col], cex=YX.m2[, w^0.4])
		tmp				<- subset(YX.m2, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m2.fit2, tmp, type='response'))	
		start			<- c( YX.m2.fit2$coef$mean + head( 2*sqrt(diag(vcov(YX.m2.fit2))), length(YX.m2.fit2$coef$mean)), YX.m2.fit2$coef$precision)
		YX.m2.fit2.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m2, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m2.fit2$coef$mean - head( 2*sqrt(diag(vcov(YX.m2.fit2))), length(YX.m2.fit2$coef$mean)), YX.m2.fit2$coef$precision)
		YX.m2.fit2.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m2, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
				c( predict(YX.m2.fit2.sup, tmp, type='response'), rev(predict(YX.m2.fit2.slw, tmp, type='response'))), 
				border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m2[, levels(stage)], fill=YX.m2[, unique(col)])
		YX.m2[, col:=NULL]
	}
	#	split missing VL into contact Yes/No and before first contact
	#	and consider every single time period 
	if(1)
	{			
		YX.m2	<- merge(YX.m2, YX.m2[, {	
							tmp<- which(!is.na(lRNA))
							list(lRNA_TL= ifelse(length(tmp)>0, lRNA[tail(tmp,1)], NA_real_))
						}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
		
		VL.win		<- c( 50, seq(100, 1000, by=100), seq(1250, 1e4, by=250) ) 
		YX.m2.VL	<- sapply(VL.win, function(VL.cur)
				{
					cat(paste('\nprocess VL=', VL.cur))									
					VL.cur	<- log10(VL.cur)
					set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
					set(YX.m2, YX.m2[,which(stage=='ART.started' & !is.na(lRNA) & lRNA<=VL.cur)],'stage', 'ART.su.Y' )
					set(YX.m2, YX.m2[,which(stage=='ART.started' & !is.na(lRNA) & lRNA>VL.cur)],'stage', 'ART.su.N' )		
					set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & !is.na(lRNA_TL) & lRNA_TL<=VL.cur)],'stage', 'ART.su.Y' )	
					set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & !is.na(lRNA_TL) & lRNA_TL>VL.cur)],'stage', 'ART.su.N' )
					set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & lRNA>VL.cur)],'stage', 'ART.su.N' )
					set(YX.m2, YX.m2[,which(stage=='ART.started')],'stage', 'ART.vlNA' )
					YX.m2[, table(stage)]
					YX.m2s		<- subset(YX.m2, stage!='ART.vlNA') 		
					set(YX.m2s, NULL, 'stage', YX.m2s[, factor(as.character(stage))])																		
					YX.m2.fit3 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m2s)
					
					ans			<- c(		YX.m2s[, length(which(stage=='ART.su.Y'))], YX.m2s[, length(which(stage=='ART.su.N'))],
							subset(YX.m2s, stage=='ART.su.Y')[, sum(w)], subset(YX.m2s, stage=='ART.su.N')[, sum(w)],
							coef(YX.m2.fit3)['stageART.su.Y'], coef(YX.m2.fit3)['stageART.su.N'], coef(YX.m2.fit3)['stageU'],	
							sqrt(diag(vcov(YX.m2.fit3)))[c('stageART.su.Y','stageART.su.N')],
							my.or.from.logit(YX.m2.fit3, 'stageART.su.Y', 'stageART.su.N', subset(YX.m2s, stage=='ART.su.Y')[, sum(w)], subset(YX.m2s, stage=='ART.su.N')[, sum(w)], 1.962),						
							my.or.from.logit(YX.m2.fit3, 'stageART.su.Y', 'stageU', subset(YX.m2s, stage=='ART.su.Y')[, sum(w)], subset(YX.m2s, stage=='U')[, sum(w)], 1.962),
							logLik(YX.m2.fit3), YX.m2.fit3$pseudo.r.squared, 10^VL.cur	)
					names(ans)	<- c('n.su.Y','n.su.N','w.su.Y','w.su.N','coef.su.Y','coef.su.N','coef.U','coef.su.Y.sd','coef.su.N.sd','or.YN','or.YN.l95','or.YN.u95','or.YU','or.YU.l95','or.YU.u95','lkl','r2','VL.thresh')
					ans
				})
		YX.m2.VL	<- as.data.table(t(YX.m2.VL))
		
		plot(YX.m2.VL[, VL.thresh], YX.m2.VL[, or.YN], type='p', pch=18, ylim= c( min(YX.m2.VL[, or.YN.l95], na.rm=TRUE),max(YX.m2.VL[, or.YN.u95], na.rm=TRUE))  )
		dummy	<- YX.m2.VL[,	lines( rep(VL.thresh,2), c(or.YN.l95,or.YN.u95))	, by='VL.thresh']
		
		VL.cur	<- 3
		set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
		set(YX.m2, YX.m2[,which(stage=='ART.started' & !is.na(lRNA) & lRNA<=VL.cur)],'stage', 'ART.su.Y' )
		set(YX.m2, YX.m2[,which(stage=='ART.started' & !is.na(lRNA) & lRNA>VL.cur)],'stage', 'ART.su.N' )		
		set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & !is.na(lRNA_TL) & lRNA_TL<=VL.cur)],'stage', 'ART.su.Y' )	
		set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & !is.na(lRNA_TL) & lRNA_TL>VL.cur)],'stage', 'ART.su.N' )
		set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & lRNA>VL.cur)],'stage', 'ART.su.N' )		
		YX.m2s		<- subset(YX.m2, stage!='ART.vlNA') 		
		set(YX.m2s, NULL, 'stage', YX.m2s[, factor(as.character(stage))])																		
		YX.m2.fit3 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m2s)
		tmp			<- c(		YX.m2s[, length(which(stage=='ART.su.Y'))], YX.m2s[, length(which(stage=='ART.su.N'))],
				subset(YX.m2s, stage=='ART.su.Y')[, sum(w)], subset(YX.m2s, stage=='ART.su.N')[, sum(w)],
				coef(YX.m2.fit3)['stageART.su.Y'], coef(YX.m2.fit3)['stageART.su.N'], coef(YX.m2.fit3)['stageU'],	
				sqrt(diag(vcov(YX.m2.fit3)))[c('stageART.su.Y','stageART.su.N')],
				my.or.from.logit(YX.m2.fit3, 'stageART.su.Y', 'stageART.su.N', subset(YX.m2s, stage=='ART.su.Y')[, sum(w)], subset(YX.m2s, stage=='ART.su.N')[, sum(w)], 1.962),						
				my.or.from.logit(YX.m2.fit3, 'stageART.su.Y', 'stageU', subset(YX.m2s, stage=='ART.su.Y')[, sum(w)], subset(YX.m2s, stage=='U')[, sum(w)], 1.962),
				logLik(YX.m2.fit3), YX.m2.fit3$pseudo.r.squared, 10^VL.cur	)
		names(tmp)	<- c('n.su.Y','n.su.N','w.su.Y','w.su.N','coef.su.Y','coef.su.N','coef.U','coef.su.Y.sd','coef.su.N.sd','or.YN','or.YN.l95','or.YN.u95','or.YU','or.YU.l95','or.YU.u95','lkl','r2','VL.thresh')
		tmp
		
	}
	#	split missing VL into contact Yes/No and before first contact
	#	and take max VL 
	if(1)
	{	
		VL.cur	<- 3
		set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
		YX.m2[, lRNA_TL:=NULL]
		YX.m2[, lRNA.c:=NULL]
		YX.m2[, PoslRNA_TL:=NULL]
		YX.m2	<- merge(YX.m2, YX.m2[, {	
							tmp<- which(!is.na(lRNA))
							list(lRNA_TL= ifelse(length(tmp)>0, lRNA[tail(tmp,1)], NA_real_), 	PoslRNA_TL=ifelse(length(tmp)>0, t[tail(tmp,1)], NA_real_))
						}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
		
		VL.win		<- c( seq(100, 1000, by=100), seq(1250, 1e4, by=250) ) 
		YX.m2.VL	<- sapply(VL.win, function(VL.cur)
				{
					cat(paste('\nprocess VL=', VL.cur))									
					VL.cur	<- log10(VL.cur)
					set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
					tmp		<- YX.m2[, {
								tmp<- which(!is.na(lRNA))
								if(!length(tmp) & any(t>t.PoslRNA_T1[1]) & all(contact=='Yes') & !is.na(lRNA_TL[1]) & lRNA_TL[1]<=VL.cur)
									lRNA.c	<- 'SuA.Y'						
								else if(!length(tmp) & any(t>t.PoslRNA_T1[1]) & all(contact=='Yes') & !is.na(lRNA_TL[1]) & lRNA_TL[1]>VL.cur)
									lRNA.c	<- 'SuA.N'						
								else if(!length(tmp))
									lRNA.c	<- 'SuA.NA'
								else if(max(lRNA[tmp])>VL.cur)
									lRNA.c	<- 'SuA.N'
								else
									lRNA.c	<- 'SuA.Y'
								list( lRNA.c= lRNA.c )
							}, by=c('Patient','t.Patient')]
					YX.m2	<- merge(YX.m2, tmp, by=c('Patient','t.Patient'))
					
					set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.c=='SuA.Y')],'stage', 'ART.suA.Y' )
					set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.c=='SuA.N')],'stage', 'ART.suA.N' )
					set(YX.m2, YX.m2[,which(stage=='ART.started')],'stage', 'ART.vlNA' )
					set(YX.m2, NULL, 'stage', YX.m2[, factor(as.character(stage))])		
					YX.m2.fit4 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m2)
					YX.m2[, lRNA.c:=NULL]
					
					ans			<- c(		YX.m2[, length(which(stage=='ART.suA.Y'))], YX.m2[, length(which(stage=='ART.suA.N'))],
							subset(YX.m2, stage=='ART.suA.Y')[, sum(w)], subset(YX.m2, stage=='ART.suA.N')[, sum(w)],
							coef(YX.m2.fit4)['stageART.suA.Y'], coef(YX.m2.fit4)['stageART.suA.N'], coef(YX.m2.fit4)['stageU'],	
							sqrt(diag(vcov(YX.m2.fit4)))[c('stageART.suA.Y','stageART.suA.N')],
							my.or.from.logit(YX.m2.fit4, 'stageART.suA.Y', 'stageART.suA.N', subset(YX.m2, stage=='ART.suA.Y')[, sum(w)], subset(YX.m2, stage=='ART.suA.N')[, sum(w)], 1.962),						
							my.or.from.logit(YX.m2.fit4, 'stageART.suA.Y', 'stageU', subset(YX.m2, stage=='ART.suA.Y')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
							logLik(YX.m2.fit4), YX.m2.fit4$pseudo.r.squared, 10^VL.cur	)
					names(ans)	<- c('n.su.Y','n.su.N','w.su.Y','w.su.N','coef.su.Y','coef.su.N','coef.U','coef.su.Y.sd','coef.su.N.sd','or.YN','or.YN.l95','or.YN.u95','or.YU','or.YU.l95','or.YU.u95','lkl','r2','VL.thresh')
					ans
				})
		YX.m2.VL	<- as.data.table(t(YX.m2.VL))
		plot(YX.m2.VL[, VL.thresh], YX.m2.VL[, or.YN], type='p', pch=18, ylim= c( min(YX.m2.VL[, or.YN.l95], na.rm=TRUE),max(YX.m2.VL[, or.YN.u95], na.rm=TRUE))  )
		dummy	<- YX.m2.VL[,	lines( rep(VL.thresh,2), c(or.YN.l95,or.YN.u95))	, by='VL.thresh']
		#
		#	using VL threshold 1e3
		#
		VL.cur	<- 3
		set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
		tmp		<- YX.m2[, {
					tmp<- which(!is.na(lRNA))
					if(!length(tmp) & any(t>t.PoslRNA_T1[1]) & all(contact=='Yes') & !is.na(lRNA_TL[1]) & lRNA_TL[1]<=VL.cur)
						lRNA.c	<- 'SuA.Y'						
					else if(!length(tmp) & any(t>t.PoslRNA_T1[1]) & all(contact=='Yes') & !is.na(lRNA_TL[1]) & lRNA_TL[1]>VL.cur)
						lRNA.c	<- 'SuA.N'						
					else if(!length(tmp))
						lRNA.c	<- 'SuA.NA'
					else if(max(lRNA[tmp])>VL.cur)
						lRNA.c	<- 'SuA.N'
					else
						lRNA.c	<- 'SuA.Y'
					list( lRNA.c= lRNA.c )
				}, by=c('Patient','t.Patient')]
		# 	SuA.N 	SuA.NA  SuA.Y 
		#	1477    947    	583 
		YX.m2	<- merge(YX.m2, tmp, by=c('Patient','t.Patient'))	
		set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.c=='SuA.Y')],'stage', 'ART.suA.Y' )
		set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.c=='SuA.N')],'stage', 'ART.suA.N' )
		set(YX.m2, YX.m2[,which(stage=='ART.started')],'stage', 'ART.vlNA' )
		YX.m2[, table(stage)]
		#stage
		#ART.suA.N ART.suA.Y  ART.vlNA   D1<=350   D1<=550    D1>550       DAm       DAy         U       UAm       UAy 
		#3172      4526        69       818      2064      2440       379       383      5732      1012       933
		set(YX.m2, NULL, 'stage', YX.m2[, factor(as.character(stage))])		
		YX.m2.fit4 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m2)
		
		YX.m2.fit4.or	<- data.table( 	ART1.su.NY= my.or.from.logit(YX.m2.fit4, 'stageART.suA.N', 'stageART.suA.Y', subset(YX.m2, stage=='ART.suA.N')[, sum(w)], subset(YX.m2, stage=='ART1.su.Y')[, sum(w)], 1.962),
				ART1.su.NY= my.or.from.logit(YX.m2.fit4, 'stageART.suA.Y', 'stageART.suA.N', subset(YX.m2, stage=='ART.suA.N')[, sum(w)], subset(YX.m2, stage=='ART1.su.Y')[, sum(w)], 1.962),
				ART1.su.N= my.or.from.logit(YX.m2.fit4, 'stageART.suA.N', 'stageU', subset(YX.m2, stage=='ART.suA.N')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				ART1.su.NA= my.or.from.logit(YX.m2.fit4, 'stageART.vlNA', 'stageU', subset(YX.m2, stage=='ART.vlNA')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				ART1.su.Y= my.or.from.logit(YX.m2.fit4, 'stageART.suA.Y', 'stageU', subset(YX.m2, stage=='ART.suA.Y')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				DAy= my.or.from.logit(YX.m2.fit4, 'stageDAy', 'stageU', subset(YX.m2, stage=='DAy')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				DAm= my.or.from.logit(YX.m2.fit4, 'stageDAm', 'stageU', subset(YX.m2, stage=='DAm')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),										
				D1l350= my.or.from.logit(YX.m2.fit4, 'stageD1<=350', 'stageU', subset(YX.m2, stage=='D1<=350')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				D1l550= my.or.from.logit(YX.m2.fit4, 'stageD1<=550', 'stageU', subset(YX.m2, stage=='D1<=550')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				D1g550= my.or.from.logit(YX.m2.fit4, 'stageD1>550', 'stageU', subset(YX.m2, stage=='D1>550')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				UAm= my.or.from.logit(YX.m2.fit4, 'stageUAm', 'stageU', subset(YX.m2, stage=='UAm')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				UAy= my.or.from.logit(YX.m2.fit4, 'stageUAy', 'stageU', subset(YX.m2, stage=='UAy')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
				DAyUAy= my.or.from.logit(YX.m2.fit4, 'stageDAy', 'stageUAy', subset(YX.m2, stage=='DAy')[, sum(w)], subset(YX.m2, stage=='UAy')[, sum(w)], 1.962),
				DAmUAm= my.or.from.logit(YX.m2.fit4, 'stageDAm', 'stageUAm', subset(YX.m2, stage=='DAm')[, sum(w)], subset(YX.m2, stage=='UAm')[, sum(w)], 1.962)										
		)
		#	   ART1.su.NY ART1.su.NY ART1.su.N ART1.su.NA ART1.su.Y      DAy      DAm     D1l350    D1l550    D1g550      UAm      UAy
		#1:   1.624383  	0.6156184 0.5987443  0.2684546 0.3685980 	1.482267 2.088981 0.7542298 0.5565142 0.6649380 1.985317 2.759259
		#2:   1.223428  	0.4625772 0.4753364  0.2089221 0.2925060 	1.142601 1.617384 0.5812614 0.4351638 0.5234374 1.545052 2.153543
		#3:   2.156743  	0.8192925 0.7541916  0.3449509 0.4644844 	1.922909 2.698085 0.9786689 0.7117045 0.8446903 2.551037 3.535343
		YX.m2.p		<- YX.m2[, list(n=round(sum(w),d=1),prop= round(sum(w)/YX.m2[,sum(w)],d=3) ),by='stage']
		#	D1<=350   DAm       UAm       U         ART.suA.N 	D1>550    D1<=550   ART.suA.Y ART.vlNA  DAy       UAy
		#	26.1  		13.0  	25.8 	168.6  		84.3 		 65.9  		53.7  	82.1   		2.7  	12.2  		24.5
		#	0.047 		0.023 	0.046 	0.302 		0.151 		0.118 		0.096 	0.147 		0.005 	0.022 		0.044
		#	plot
		tmp			<- data.table(stage= YX.m2[, levels(stage)], col=sapply( rainbow(YX.m2[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m2		<- merge(YX.m2, tmp, by='stage')
		plot( seq_len(nrow(YX.m2)), YX.m2[, score.Y], pch=18, col= YX.m2[, col], cex=YX.m2[, w^0.4])
		tmp				<- subset(YX.m2, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m2.fit4, tmp, type='response'))	
		start			<- c( YX.m2.fit4$coef$mean + head( 2*sqrt(diag(vcov(YX.m2.fit4))), length(YX.m2.fit4$coef$mean)), YX.m2.fit4$coef$precision)
		YX.m2.fit4.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m2, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m2.fit4$coef$mean - head( 2*sqrt(diag(vcov(YX.m2.fit4))), length(YX.m2.fit4$coef$mean)), YX.m2.fit4$coef$precision)
		YX.m2.fit4.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m2, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
				c( predict(YX.m2.fit4.sup, tmp, type='response'), rev(predict(YX.m2.fit4.slw, tmp, type='response'))), 
				border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m2[, levels(stage)], fill=YX.m2[, unique(col)])
		YX.m2[, col:=NULL]		
	}	
}
######################################################################################
project.athena.Fisheretal.YX.model2.stratify.VL1stsu<- function(YX.m2, df.all, df.viro, df.immu, lRNA.supp=log10(1e3))
{
	#YX.m2	<- copy(YX)
	gc()
	if(is.null(YX.m2))	return(NULL)
	YX.m2[, U.score:=NULL]
	#score.Y.cut<- 1e-8
	#set(YX.m2, YX.m2[,which(score.Y<score.Y.cut)], 'score.Y', score.Y.cut)
	#score.Y.cut<- 1e-5
	#set(YX.m2, YX.m2[,which(score.Y>1-score.Y.cut)], 'score.Y', 1-score.Y.cut)
	if('score.Y'%in%colnames(YX.m2) && YX.m2[, !any(score.Y>1.1)])
	{
		tmp	<- YX.m2[, score.Y>0]
		set(YX.m2, which(tmp), 'score.Y', YX.m2[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )
	}		
	#	
	#	set CD41st
	#
	cat(paste('\nsetting CD4 1st\n'))
	cd4.label	<- c('D1l350','D1l500','D1g500')
	cd4.cut		<- c(-1, 350, 500, 5000)
	tmp			<- subset(df.all, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
	setkey(tmp, Patient)
	tmp			<- unique(tmp)	
	tmp[, CD41st:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
	setnames(tmp, 'Patient', 't.Patient')
	set(tmp, tmp[, which(is.na(CD41st))], 'CD41st', 'D1.NA')
	gc()
	YX.m2	<- merge(YX.m2, subset(tmp, select=c(t.Patient, CD41st)), by='t.Patient', all.x=TRUE)	
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Yes')], 'CD41st', 'UAy')
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Maybe')], 'CD41st', 'UAm')
	set(YX.m2, YX.m2[, which(stage=='U' & (t.isAcute=='No' | is.na(t.isAcute)))], 'CD41st', 'U')
	set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD41st', 'DAy' )
	set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD41st', 'DAm' )
	set(YX.m2, YX.m2[, which(stage=='ART.started')], 'CD41st', 'ART.started' )	
	#	
	#	set CD4t
	#	
	cat(paste('\nsetting CD4t\n'))
	cd4.label	<- c('Dtl350','Dtl500','Dtg500')
	cd4.cut		<- c(350, 500, 5000)
	tmp			<- copy(df.immu)
	set(tmp, NULL, 'PosCD4', hivc.db.Date2numeric(tmp[, PosCD4]) )
	tmp			<- tmp[, {
				z			<- lapply(cd4.cut, function(x)	ifelse(any(CD4<x), PosCD4[ which(CD4<x)[1] ], NA_real_)	)
				names(z)	<- cd4.label
				z
			}, by='Patient']
	setnames(tmp, 'Patient', 't.Patient')	
	gc()
	YX.m2		<- merge(YX.m2, tmp, by='t.Patient', all.x=TRUE)
	YX.m2[, CD4t:=NA_character_]
	set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtg500)], 'CD4t', 'Dtg500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtl500)], 'CD4t', 'Dtl500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtl350)], 'CD4t', 'Dtl350')
	#	missing data: all first CD4>500 go to Dtg500
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & !is.na(Dtg500) & is.na(Dtl500) & is.na(Dtl350))], 'CD4t', 'Dtg500')
	#	missing data: allow grace of one year
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtl350)], 'CD4t', 'Dtl350')
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtl500)], 'CD4t', 'Dtl500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtg500)], 'CD4t', 'Dtg500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) )], 'CD4t', 'Dt.NA')
	#	undiagnosed
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Yes')], 'CD4t', 'UAy')
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Maybe')], 'CD4t', 'UAm')
	set(YX.m2, YX.m2[, which(stage=='U' & (t.isAcute=='No' | is.na(t.isAcute)))], 'CD4t', 'U')
	#	separate acute
	set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD4t', 'DAy' )
	set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD4t', 'DAm' )		
	set(YX.m2, YX.m2[, which(stage=='ART.started')], 'CD4t', 'ART.started' )
	#	combine acute
	YX.m2[, CD4a:=CD4t]
	set(YX.m2, YX.m2[, which(CD4t=='UAy' | CD4t=='UAm')], 'CD4a', 'UA')
	set(YX.m2, YX.m2[, which(CD4t=='DAy' | CD4t=='DAm')], 'CD4a', 'DA')
	#	treat missing acute separately
	YX.m2[, CD4b:=CD4a]
	set(YX.m2, YX.m2[, which(stage=='U' & is.na(t.isAcute))], 'CD4b', 'UAna')		
	#	merge PoslRNA_T1
	cat(paste('\nsetting PoslRNA_T1\n'))
	tmp			<- unique(subset( df.all, select=c(Patient, PoslRNA_T1) ))
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	YX.m2		<- merge(YX.m2, tmp, by= 't.Patient', all.x=TRUE)	
	cat(paste('\nsubset\n'))
	if('score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, score.Y, stage, CDCC, lRNA, t.isAcute, t.PoslRNA_T1, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, contact, fw.up.med, t.period, w, w.i, w.in, w.t, w.tn, CD41st, CD4t, CD4a, CD4b, t.Age, t.RegionHospital  ))	
	if(!'score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, stage, CDCC, lRNA, t.isAcute, t.PoslRNA_T1, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, contact, fw.up.med, t.period, CD41st, CD4t, CD4a, CD4b  ))	
	gc()
	#
	#	add suppressed to 'stage'
	#
	YX.m2[, stage.orig:=stage]
	#
	#	using given VL threshold 
	#		
	cat(paste('\nsettting PosRNA, lRNA\n'))
	set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
	tmp			<- subset(df.viro, select=c(Patient, PosRNA, lRNA))
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	set(tmp, NULL, 't.PosRNA', hivc.db.Date2numeric(tmp[,t.PosRNA]))
	cat(paste('\nsettting t.ARTSu_T1\n'))
	tmp			<- merge( unique(subset(YX.m2, select=c(t.Patient, t.AnyT_T1))), tmp, by='t.Patient' )
	tmp			<- subset(tmp, t.AnyT_T1<=t.PosRNA)[, 	{
				z<- which(t.lRNA<lRNA.supp)
				list(t.ARTSu_T1= ifelse(length(z), t.PosRNA[z[1]], NA_real_) )	
			}, by='t.Patient']
	YX.m2		<- merge(YX.m2, tmp, by= 't.Patient', all.x=TRUE)
	set(YX.m2, YX.m2[, which(stage=='ART.started' & !is.na(t.ARTSu_T1) & t>=t.ARTSu_T1)], 'stage', 'ART1.su.Y')
	set(YX.m2, YX.m2[, which(stage=='ART.started')], 'stage', 'ART1.su.N')
	#	set CD41st and CD4t
	tmp			<- YX.m2[, which(stage.orig=='ART.started')] 
	set(YX.m2, tmp, 'CD4t', YX.m2[tmp, stage])	
	set(YX.m2, NULL, 'CD4t', YX.m2[, factor(as.character(CD4t))])
	set(YX.m2, tmp, 'CD4b', YX.m2[tmp, stage])	
	set(YX.m2, NULL, 'CD4b', YX.m2[, factor(as.character(CD4b))])	
	set(YX.m2, tmp, 'CD41st', YX.m2[tmp, stage])
	set(YX.m2, NULL, 'CD41st', YX.m2[, factor(as.character(CD41st))])
	#	set CD4a
	YX.m2		<- merge(YX.m2, YX.m2[, {
						ans	<- list(t=t, lRNA.c2=stage)
						if( !all(stage=='ART1.su.Y') )
							ans$lRNA.c2[ stage=='ART1.su.Y' ]	<- 'ART1.su.M'
						ans
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient','t'), all.x=TRUE)
	tmp			<- YX.m2[, which(stage.orig=='ART.started')]
	set(YX.m2, tmp, 'CD4a', YX.m2[tmp, lRNA.c2])
	set(YX.m2, NULL, 'CD4a', YX.m2[, factor(as.character(CD4a))])	
	YX.m2[, t.ARTSu_T1:=NULL]
	YX.m2[, stage.orig:=NULL]
	YX.m2[, lRNA.c2:=NULL]
	gc()
	#
	#	add tperiod
	#
	if('t.period'%in%colnames(YX.m2))
	{
		cat(paste('\nadding CD4t.tperiod\n'))	
		YX.m2[, CD41st.tperiod:= paste(CD41st, t.period,sep='.')]
		set(YX.m2, NULL, 'CD41st.tperiod', YX.m2[, factor(as.character(CD41st.tperiod))])
		YX.m2[, CD4t.tperiod:= paste(CD4t, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4t.tperiod', YX.m2[, factor(as.character(CD4t.tperiod))])
		YX.m2[, CD4a.tperiod:= paste(CD4a, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4a.tperiod', YX.m2[, factor(as.character(CD4a.tperiod))])
		YX.m2[, CD4b.tperiod:= paste(CD4b, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4b.tperiod', YX.m2[, factor(as.character(CD4b.tperiod))])
	}
	gc()	
	YX.m2
}
######################################################################################
project.athena.Fisheretal.YX.model2.stratify.VLt<- function(YX.m2, df.all, df.viro, df.immu, lRNA.supp=log10(1e3), plot.file.or=NA )
{
	#YX.m2	<- copy(YX)
	gc()
	if(is.null(YX.m2))	return(NULL)
	YX.m2[, U.score:=NULL]
	#score.Y.cut<- 1e-8
	#set(YX.m2, YX.m2[,which(score.Y<score.Y.cut)], 'score.Y', score.Y.cut)
	#score.Y.cut<- 1e-5
	#set(YX.m2, YX.m2[,which(score.Y>1-score.Y.cut)], 'score.Y', 1-score.Y.cut)
	if('score.Y'%in%colnames(YX.m2) && YX.m2[, !any(score.Y>1.1)])
	{
		tmp	<- YX.m2[, score.Y>0]
		set(YX.m2, which(tmp), 'score.Y', YX.m2[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )
	}		
	#	
	#	set CD41st
	#
	cat(paste('\nsetting CD4 1st\n'))
	cd4.label	<- c('D1l350','D1l500','D1g500')
	cd4.cut		<- c(-1, 350, 500, 5000)
	tmp			<- subset(df.all, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
	setkey(tmp, Patient)
	tmp			<- unique(tmp)	
	tmp[, CD41st:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
	setnames(tmp, 'Patient', 't.Patient')
	set(tmp, tmp[, which(is.na(CD41st))], 'CD41st', 'D1.NA')
	gc()
	YX.m2	<- merge(YX.m2, subset(tmp, select=c(t.Patient, CD41st)), by='t.Patient', all.x=TRUE)	
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Yes')], 'CD41st', 'UAy')
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Maybe')], 'CD41st', 'UAm')
	set(YX.m2, YX.m2[, which(stage=='U' & (t.isAcute=='No' | is.na(t.isAcute)))], 'CD41st', 'U')
	set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD41st', 'DAy' )
	set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD41st', 'DAm' )
	set(YX.m2, YX.m2[, which(stage=='ART.started')], 'CD41st', 'ART.started' )	
	#	
	#	set CD4t
	#	
	cat(paste('\nsetting CD4t\n'))
	cd4.label	<- c('Dtl350','Dtl500','Dtg500')
	cd4.cut		<- c(350, 500, 5000)
	tmp			<- copy(df.immu)
	set(tmp, NULL, 'PosCD4', hivc.db.Date2numeric(tmp[, PosCD4]) )
	tmp			<- tmp[, {
				z			<- lapply(cd4.cut, function(x)	ifelse(any(CD4<x), PosCD4[ which(CD4<x)[1] ], NA_real_)	)
				names(z)	<- cd4.label
				z
			}, by='Patient']
	setnames(tmp, 'Patient', 't.Patient')		
	gc()
	YX.m2		<- merge(YX.m2, tmp, by='t.Patient', all.x=TRUE)
	YX.m2[, CD4t:=NA_character_]
	set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtg500)], 'CD4t', 'Dtg500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtl500)], 'CD4t', 'Dtl500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtl350)], 'CD4t', 'Dtl350')
	#	missing data: all first CD4>500 go to Dtg500
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & !is.na(Dtg500) & is.na(Dtl500) & is.na(Dtl350))], 'CD4t', 'Dtg500')
	#	missing data: allow grace of one year
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtl350)], 'CD4t', 'Dtl350')
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtl500)], 'CD4t', 'Dtl500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtg500)], 'CD4t', 'Dtg500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) )], 'CD4t', 'Dt.NA')
	#	undiagnosed
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Yes')], 'CD4t', 'UAy')
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Maybe')], 'CD4t', 'UAm')
	set(YX.m2, YX.m2[, which(stage=='U' & (t.isAcute=='No' | is.na(t.isAcute)))], 'CD4t', 'U')
	#	separate acute
	set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD4t', 'DAy' )
	set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD4t', 'DAm' )		
	set(YX.m2, YX.m2[, which(stage=='ART.started')], 'CD4t', 'ART.started' )
	#	combine acute
	YX.m2[, CD4a:=CD4t]
	set(YX.m2, YX.m2[, which(CD4t=='UAy' | CD4t=='UAm')], 'CD4a', 'UA')
	set(YX.m2, YX.m2[, which(CD4t=='DAy' | CD4t=='DAm')], 'CD4a', 'DA')		
	#	treat missing acute separately
	YX.m2[, CD4b:=CD4a]
	set(YX.m2, YX.m2[, which(stage=='U' & is.na(t.isAcute))], 'CD4b', 'UAna')
	#	merge PoslRNA_T1
	cat(paste('\nsetting PoslRNA_T1\n'))
	tmp			<- unique(subset( df.all, select=c(Patient, PoslRNA_T1) ))
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	YX.m2		<- merge(YX.m2, tmp, by= 't.Patient', all.x=TRUE)
	# 	merge PoslRNA_TL and lRNA_TL
	tmp			<- merge( unique(subset(df.all, select=Patient)), subset(df.viro, select=c(Patient, PosRNA, lRNA)), by='Patient' )
	set( tmp, NULL, 'PosRNA', hivc.db.Date2numeric(tmp[,PosRNA]) )
	tmp			<- subset( tmp[, 	{
										tmp<- which(!is.na(lRNA))
										list(t.Patient=Patient[1], lRNA_TL=lRNA[tail(tmp,1)], PoslRNA_TL=PosRNA[tail(tmp,1)])
									}, by='Patient'], select=c(t.Patient, lRNA_TL, PoslRNA_TL) )
	YX.m2		<- merge(YX.m2, tmp, by= 't.Patient', all.x=TRUE)				
	#
	cat(paste('\nsubset\n'))
	if('score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, score.Y, stage, CDCC, lRNA, t.isAcute, t.PoslRNA_T1, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, contact, fw.up.med, t.period, w, w.i, w.in, w.t, w.tn, lRNA_TL, PoslRNA_TL, CD41st, CD4t, CD4a, CD4b, t.Age, t.RegionHospital  ))	
	if(!'score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, stage, CDCC, lRNA, t.isAcute, t.PoslRNA_T1, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, contact, fw.up.med, t.period, lRNA_TL, PoslRNA_TL, CD41st, CD4t, CD4a, CD4b  ))	
	#
	#	add suppressed to 'stage'
	#
	YX.m2[, stage.orig:=stage]
	gc()
	#
	#	using given VL threshold 
	#		
	set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
	set(YX.m2, YX.m2[,which(stage=='ART.started' & !is.na(lRNA) & lRNA<=lRNA.supp)],'stage', 'ART.su.Y' )
	set(YX.m2, YX.m2[,which(stage=='ART.started' & !is.na(lRNA) & lRNA>lRNA.supp)],'stage', 'ART.su.N' )		
	set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & !is.na(lRNA_TL) & lRNA_TL<=lRNA.supp)],'stage', 'ART.su.Y' )	
	set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & !is.na(lRNA_TL) & lRNA_TL>lRNA.supp)],'stage', 'ART.su.N' )
	set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & lRNA>lRNA.supp)],'stage', 'ART.su.N' )		
	set(YX.m2, YX.m2[,which(stage=='ART.started')],'stage', 'ART.vlNA' )
	#	set CD41st and CD4t
	tmp			<- YX.m2[, which(stage.orig=='ART.started')] 
	set(YX.m2, tmp, 'CD4t', YX.m2[tmp, stage])
	set(YX.m2, tmp, 'CD41st', YX.m2[tmp, stage])
	set(YX.m2, tmp, 'CD4b', YX.m2[tmp, stage])
	set(YX.m2, NULL, 'CD4t', YX.m2[, factor(as.character(CD4t))])
	set(YX.m2, NULL, 'CD41st', YX.m2[, factor(as.character(CD41st))])
	set(YX.m2, NULL, 'CD4b', YX.m2[, factor(as.character(CD4b))])
	#	set CD4a
	YX.m2		<- merge(YX.m2, YX.m2[, {
					ans	<- list(t=t, lRNA.c2=stage)
					if( !all(stage=='ART.su.Y') )
						ans$lRNA.c2[ stage=='ART.su.Y' ]	<- 'ART.su.M'
					ans
				}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient','t'), all.x=TRUE)
	tmp			<- YX.m2[, which(stage.orig=='ART.started')]
	set(YX.m2, tmp, 'CD4a', YX.m2[tmp, lRNA.c2])
	set(YX.m2, NULL, 'CD4a', YX.m2[, factor(as.character(CD4a))])	
	#	clean up
	YX.m2[, stage.orig:=NULL]
	YX.m2[, lRNA.c2:=NULL]
	YX.m2[, PoslRNA_TL:=NULL]
	YX.m2[, lRNA_TL:=NULL]
	gc()
	#
	#	add tperiod
	#
	if('t.period'%in%colnames(YX.m2))
	{
		cat(paste('\nadding CD4t.tperiod\n'))	
		YX.m2[, CD41st.tperiod:= paste(CD41st, t.period,sep='.')]
		set(YX.m2, NULL, 'CD41st.tperiod', YX.m2[, factor(as.character(CD41st.tperiod))])
		YX.m2[, CD4t.tperiod:= paste(CD4t, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4t.tperiod', YX.m2[, factor(as.character(CD4t.tperiod))])
		YX.m2[, CD4a.tperiod:= paste(CD4a, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4a.tperiod', YX.m2[, factor(as.character(CD4a.tperiod))])
		YX.m2[, CD4b.tperiod:= paste(CD4b, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4b.tperiod', YX.m2[, factor(as.character(CD4b.tperiod))])
	}
	gc()	
	YX.m2
}
######################################################################################
project.athena.Fisheretal.YX.model2.stratify.VLgm<- function(YX.m2, df.all, df.viro, df.immu, lRNA.supp=log10(1e3), plot.file.or=NA )
{
	#YX.m2	<- copy(YX)
	gc()
	if(is.null(YX.m2))	return(NULL)
	YX.m2[, U.score:=NULL]
	#score.Y.cut<- 1e-8
	#set(YX.m2, YX.m2[,which(score.Y<score.Y.cut)], 'score.Y', score.Y.cut)
	#score.Y.cut<- 1e-5
	#set(YX.m2, YX.m2[,which(score.Y>1-score.Y.cut)], 'score.Y', 1-score.Y.cut)
	if('score.Y'%in%colnames(YX.m2) && YX.m2[, !any(score.Y>1.1)])
	{
		tmp	<- YX.m2[, score.Y>0]
		set(YX.m2, which(tmp), 'score.Y', YX.m2[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )
	}		
	#	
	#	set CD41st
	#
	cat(paste('\nsetting CD4 1st\n'))
	cd4.label	<- c('D1l350','D1l500','D1g500')
	cd4.cut		<- c(-1, 350, 500, 5000)
	tmp			<- subset(df.all, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
	setkey(tmp, Patient)
	tmp			<- unique(tmp)	
	tmp[, CD41st:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
	setnames(tmp, 'Patient', 't.Patient')
	set(tmp, tmp[, which(is.na(CD41st))], 'CD41st', 'D1.NA')
	gc()
	YX.m2	<- merge(YX.m2, subset(tmp, select=c(t.Patient, CD41st)), by='t.Patient', all.x=TRUE)	
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Yes')], 'CD41st', 'UAy')
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Maybe')], 'CD41st', 'UAm')
	set(YX.m2, YX.m2[, which(stage=='U' & (t.isAcute=='No' | is.na(t.isAcute)))], 'CD41st', 'U')
	set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD41st', 'DAy' )
	set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD41st', 'DAm' )
	set(YX.m2, YX.m2[, which(stage=='ART.started')], 'CD41st', 'ART.started' )	
	#	
	#	set CD4t
	#	
	cat(paste('\nsetting CD4t\n'))
	cd4.label	<- c('Dtl350','Dtl500','Dtg500')
	cd4.cut		<- c(350, 500, 5000)
	tmp			<- copy(df.immu)
	set(tmp, NULL, 'PosCD4', hivc.db.Date2numeric(tmp[, PosCD4]) )
	tmp			<- tmp[, {
				z			<- lapply(cd4.cut, function(x)	ifelse(any(CD4<x), PosCD4[ which(CD4<x)[1] ], NA_real_)	)
				names(z)	<- cd4.label
				z
			}, by='Patient']
	setnames(tmp, 'Patient', 't.Patient')	
	gc()
	YX.m2		<- merge(YX.m2, tmp, by='t.Patient', all.x=TRUE)
	YX.m2[, CD4t:=NA_character_]
	set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtg500)], 'CD4t', 'Dtg500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtl500)], 'CD4t', 'Dtl500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtl350)], 'CD4t', 'Dtl350')
	#	missing data: all first CD4>500 go to Dtg500
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & !is.na(Dtg500) & is.na(Dtl500) & is.na(Dtl350))], 'CD4t', 'Dtg500')
	#	missing data: allow grace of one year
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtl350)], 'CD4t', 'Dtl350')
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtl500)], 'CD4t', 'Dtl500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtg500)], 'CD4t', 'Dtg500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) )], 'CD4t', 'Dt.NA')
	#	undiagnosed
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Yes')], 'CD4t', 'UAy')
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Maybe')], 'CD4t', 'UAm')
	set(YX.m2, YX.m2[, which(stage=='U' & (t.isAcute=='No' | is.na(t.isAcute)))], 'CD4t', 'U')
	#	separate acute
	set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD4t', 'DAy' )
	set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD4t', 'DAm' )		
	set(YX.m2, YX.m2[, which(stage=='ART.started')], 'CD4t', 'ART.started' )
	#	combine acute
	YX.m2[, CD4a:=CD4t]
	set(YX.m2, YX.m2[, which(CD4t=='UAy' | CD4t=='UAm')], 'CD4a', 'UA')
	set(YX.m2, YX.m2[, which(CD4t=='DAy' | CD4t=='DAm')], 'CD4a', 'DA')	
	#	treat missing acute separately
	YX.m2[, CD4b:=CD4a]
	set(YX.m2, YX.m2[, which(stage=='U' & is.na(t.isAcute))], 'CD4b', 'UAna')	
	# collect PoslRNA_TL and lRNA_TL
	cat(paste('\nmerge df.viro\n'))
	tmp		<- merge( unique(subset(clumsm.info, select=Patient)), subset(df.viro, select=c(Patient, PosRNA, lRNA)), by='Patient' )
	set( tmp, NULL, 'PosRNA', hivc.db.Date2numeric(tmp[,PosRNA]) )
	tmp		<- subset( tmp[, 	{
						tmp<- which(!is.na(lRNA))
						list(t.Patient=Patient[1], lRNA_TL=lRNA[tail(tmp,1)], PoslRNA_TL=PosRNA[tail(tmp,1)])
					}, by='Patient'], select=c(t.Patient, lRNA_TL, PoslRNA_TL) )
	YX.m2	<- merge(YX.m2, tmp, by= 't.Patient', all.x=TRUE)
	#define lRNA.mxw and lRNA.gm
	YX.m2	<- merge(YX.m2, YX.m2[, {
						tmp<- which(!is.na(lRNA))					
						if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1])
							lRNA.gm		<- lRNA_TL[1]						
						else if(!length(tmp))
							lRNA.gm		<- NA_real_
						else
							lRNA.gm		<- mean(lRNA[tmp])
						list( lRNA.gm= lRNA.gm )
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))		
	cat(paste('\nsubset\n'))
	if('score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, score.Y, stage, CDCC, lRNA, t.isAcute, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, contact, fw.up.med, t.period, w, w.i, w.in, w.t, w.tn, CD41st, CD4t, CD4a, CD4b, lRNA.gm, t.Age, t.RegionHospital  ))	
	if(!'score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, stage, CDCC, lRNA, t.isAcute, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, contact, fw.up.med, t.period, CD41st, CD4t, CD4a, CD4b, lRNA.gm  ))		
	gc()
	#
	#	add suppressed to 'stage'
	#
	YX.m2[, stage.orig:=stage]	
	#
	#	using given VL threshold 
	#	
	gc()
	set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.gm<=lRNA.supp)],'stage', 'ART.suA.Y' )
	set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.gm>lRNA.supp)],'stage', 'ART.suA.N' )
	set(YX.m2, YX.m2[,which(stage=='ART.started')],'stage', 'ART.vlNA' )
	#	set CD41st and CD4t
	tmp			<- YX.m2[, which(stage.orig=='ART.started')] 
	set(YX.m2, tmp, 'CD4t', YX.m2[tmp, stage])
	set(YX.m2, NULL, 'CD4t', YX.m2[, factor(as.character(CD4t))])
	set(YX.m2, tmp, 'CD4b', YX.m2[tmp, stage])
	set(YX.m2, NULL, 'CD4b', YX.m2[, factor(as.character(CD4b))])	
	set(YX.m2, tmp, 'CD41st', YX.m2[tmp, stage])	
	set(YX.m2, NULL, 'CD41st', YX.m2[, factor(as.character(CD41st))])	
	YX.m2[, stage.orig:=NULL]
	YX.m2[, lRNA.c:=NULL]
	#
	#	add tperiod
	#
	if('t.period'%in%colnames(YX.m2))
	{
		cat(paste('\nadding CD4t.tperiod\n'))	
		YX.m2[, CD41st.tperiod:= paste(CD41st, t.period,sep='.')]
		set(YX.m2, NULL, 'CD41st.tperiod', YX.m2[, factor(as.character(CD41st.tperiod))])
		YX.m2[, CD4t.tperiod:= paste(CD4t, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4t.tperiod', YX.m2[, factor(as.character(CD4t.tperiod))])
		YX.m2[, CD4b.tperiod:= paste(CD4b, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4b.tperiod', YX.m2[, factor(as.character(CD4b.tperiod))])		
	}
	gc()
	YX.m2
}
######################################################################################
project.athena.Fisheretal.YX.model4.stratify.Diagnosed<- function(YX.m2, df.immu, return.only.Diag=TRUE )
{
	#YX.m2	<- copy(YX)
	gc()
	if(is.null(YX.m2))	return(NULL)
	YX.m2[, U.score:=NULL]
	if('score.Y'%in%colnames(YX.m2) && YX.m2[, !any(score.Y>1.1)])
	{
		tmp	<- YX.m2[, score.Y>0]
		set(YX.m2, which(tmp), 'score.Y', YX.m2[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )
	}		
	#	
	#	set CD4t
	#	
	cat(paste('\nsetting CD4t\n'))
	cd4.label	<- c('Dtl350','Dtl500','Dtg500')
	cd4.cut		<- c(350, 500, 5000)
	tmp			<- copy(df.immu)
	set(tmp, NULL, 'PosCD4', hivc.db.Date2numeric(tmp[, PosCD4]) )
	tmp			<- tmp[, {
				z			<- lapply(cd4.cut, function(x)	ifelse(any(CD4<x), PosCD4[ which(CD4<x)[1] ], NA_real_)	)
				names(z)	<- cd4.label
				z
			}, by='Patient']
	setnames(tmp, 'Patient', 't.Patient')	
	gc()
	YX.m2		<- merge(YX.m2, tmp, by='t.Patient', all.x=TRUE)
	YX.m2[, CD4t:=NA_character_]
	set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtg500)], 'CD4t', 'Dtg500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtl500)], 'CD4t', 'Dtl500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtl350)], 'CD4t', 'Dtl350')
	#	missing data: all first CD4>500 go to Dtg500
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & !is.na(Dtg500) & is.na(Dtl500) & is.na(Dtl350))], 'CD4t', 'Dtg500')
	#	missing data: allow grace of one year
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtl350)], 'CD4t', 'Dtl350')
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtl500)], 'CD4t', 'Dtl500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtg500)], 'CD4t', 'Dtg500')
	set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) )], 'CD4t', 'Dt.NA')
	#	undiagnosed
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Yes')], 'CD4t', 'UAy')
	set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Maybe')], 'CD4t', 'UAm')
	set(YX.m2, YX.m2[, which(stage=='U' & (t.isAcute=='No' | is.na(t.isAcute)))], 'CD4t', 'U')
	#	treated, so all should be different from NA now
	set(YX.m2, YX.m2[, which(stage=='ART.started')], 'CD4t', 'ART.started' )
	stopifnot( YX.m2[, length(which(is.na(CD4t)))]==0 )
	#
	#	set acute as independent additive effect
	#
	YX.m2[, Acute:='No']
	set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'Acute', 'Yes' )
	set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'Acute', 'Maybe' )
	set(YX.m2, NULL, 'Acute', YX.m2[, factor(Acute, levels=c('No','Maybe','Yes'))])
	#	CDCC
	set(YX.m2, NULL, 'CDCC', YX.m2[, factor(CDCC)])
	#
	#	set acute as independent additive effect and ignore missing Acute
	#
	YX.m2[, AcuteNo:=NA_character_]	
	set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'AcuteNo', 'Yes' )
	set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 >= 0.25)], 'AcuteNo', 'No' )
	set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'AcuteNo', 'Maybe' )
	set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 >= 0.25)], 'AcuteNo', 'No' )
	set(YX.m2, YX.m2[, which(t.isAcute=='No' & stage=='Diag')], 'AcuteNo', 'No' )	
	set(YX.m2, NULL, 'AcuteNo', YX.m2[, factor(AcuteNo, levels=c('No','Maybe','Yes'))])	
	#
	#	add max viral load
	#
	gc()
	cat(paste('\nadding lRNA.mx\n'))
	YX.m2	<- merge(YX.m2, YX.m2[, {
						tmp<- which(!is.na(lRNA))
						list( lRNA.mx= ifelse(length(tmp), max(lRNA[tmp]), NA_real_) )
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)
	gc()
	#
	#	add tperiod
	#
	if('t.period'%in%colnames(YX.m2))
	{
		cat(paste('\nadding CD4t.tperiod\n'))	
		YX.m2[, CD4t.tperiod:= paste(CD4t, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4t.tperiod', YX.m2[, factor(as.character(CD4t.tperiod))])		
	}	
	#	subset
	cat(paste('\nsubset\n'))
	if('score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, score.Y, stage, CDCC, lRNA, lRNA.mx, t.isAcute, t.AnyT_T1, contact, fw.up.med, t.period, w, w.i, w.in, w.t, w.tn, CD4t, CD4t.tperiod, Acute, AcuteNo, t.Age, t.RegionHospital, t2.care.t1  ))	
	if(!'score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, stage, CDCC, lRNA, lRNA.mx, t.isAcute, t.AnyT_T1, contact, fw.up.med, t.period, CD4t, CD4t.tperiod, Acute, AcuteNo, t2.care.t1  ))
	if(return.only.Diag)
		YX.m2	<- subset(YX.m2, stage=='Diag')
	gc()
	YX.m2
}
######################################################################################
project.athena.Fisheretal.YX.model2.stratify.VLmxwindow<- function(YX.m2, df.all, df.viro, df.immu, lRNA.supp=3, plot.file.or=NA, resume=1, save.file=NA )
{
	if(resume && !is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(!resume || is.na(save.file) || inherits(readAttempt, "try-error"))
	{
		#YX.m2	<- copy(YX)
		gc()
		if(is.null(YX.m2))	return(NULL)
		YX.m2[, U.score:=NULL]
		#score.Y.cut<- 1e-8
		#set(YX.m2, YX.m2[,which(score.Y<score.Y.cut)], 'score.Y', score.Y.cut)
		#score.Y.cut<- 1e-5
		#set(YX.m2, YX.m2[,which(score.Y>1-score.Y.cut)], 'score.Y', 1-score.Y.cut)
		if('score.Y'%in%colnames(YX.m2) && YX.m2[, !any(score.Y>1.1)] )
		{
			tmp	<- YX.m2[, score.Y>0]
			set(YX.m2, which(tmp), 'score.Y', YX.m2[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )
		}		
		#	
		#	set CD41st
		#
		cat(paste('\nsetting CD4 1st\n'))
		cd4.label	<- c('D1l350','D1l500','D1g500')
		cd4.cut		<- c(-1, 350, 500, 5000)
		tmp			<- subset(df.all, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
		setkey(tmp, Patient)
		tmp			<- unique(tmp)	
		tmp[, CD41st:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
		setnames(tmp, 'Patient', 't.Patient')
		set(tmp, tmp[, which(is.na(CD41st))], 'CD41st', 'D1.NA')
		gc()
		YX.m2	<- merge(YX.m2, subset(tmp, select=c(t.Patient, CD41st)), by='t.Patient', all.x=TRUE)	
		set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Yes')], 'CD41st', 'UAy')
		set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Maybe')], 'CD41st', 'UAm')
		set(YX.m2, YX.m2[, which(stage=='U' & (t.isAcute=='No' | is.na(t.isAcute)))], 'CD41st', 'U')
		set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD41st', 'DAy' )
		set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD41st', 'DAm' )
		set(YX.m2, YX.m2[, which(stage=='ART.started')], 'CD41st', 'ART.started' )	
		#	
		#	set CD4t
		#	
		cat(paste('\nsetting CD4t\n'))
		cd4.label	<- c('Dtl350','Dtl500','Dtg500')
		cd4.cut		<- c(350, 500, 5000)
		tmp			<- copy(df.immu)
		set(tmp, NULL, 'PosCD4', hivc.db.Date2numeric(tmp[, PosCD4]) )
		tmp			<- tmp[, {
					z			<- lapply(cd4.cut, function(x)	ifelse(any(CD4<x), PosCD4[ which(CD4<x)[1] ], NA_real_)	)
					names(z)	<- cd4.label
					z
				}, by='Patient']
		setnames(tmp, 'Patient', 't.Patient')	
		gc()
		YX.m2		<- merge(YX.m2, tmp, by='t.Patient', all.x=TRUE)
		YX.m2[, CD4t:=NA_character_]
		set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtg500)], 'CD4t', 'Dtg500')
		set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtl500)], 'CD4t', 'Dtl500')
		set(YX.m2, YX.m2[, which(stage=='Diag' & t>Dtl350)], 'CD4t', 'Dtl350')
		#	missing data: all first CD4>500 go to Dtg500
		set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & !is.na(Dtg500) & is.na(Dtl500) & is.na(Dtl350))], 'CD4t', 'Dtg500')
		#	missing data: allow grace of one year
		if(0)
		{
			set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtl350)], 'CD4t', 'Dtl350')
			set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtl500)], 'CD4t', 'Dtl500')
			set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) & t+1>Dtg500)], 'CD4t', 'Dtg500')			
		}
		set(YX.m2, YX.m2[, which(stage=='Diag' & is.na(CD4t) )], 'CD4t', 'Dt.NA')
		#	undiagnosed
		set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Yes')], 'CD4t', 'UAy')
		set(YX.m2, YX.m2[, which(stage=='U' & t.isAcute=='Maybe')], 'CD4t', 'UAm')
		set(YX.m2, YX.m2[, which(stage=='U' & (t.isAcute=='No' | is.na(t.isAcute)))], 'CD4t', 'U')
		#	separate acute
		set(YX.m2, YX.m2[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD4t', 'DAy' )
		set(YX.m2, YX.m2[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'CD4t', 'DAm' )		
		set(YX.m2, YX.m2[, which(stage=='ART.started')], 'CD4t', 'ART.started' )	
		#	combine acute
		YX.m2[, CD4a:=CD4t]
		set(YX.m2, YX.m2[, which(CD4t=='UAy' | CD4t=='UAm')], 'CD4a', 'UA')
		set(YX.m2, YX.m2[, which(CD4t=='DAy' | CD4t=='DAm')], 'CD4a', 'DA')
		#	treat missing acute separately
		YX.m2[, CD4b:=CD4a]
		set(YX.m2, YX.m2[, which(stage=='U' & is.na(t.isAcute))], 'CD4b', 'UAna')	
		#	set UAna for CD4t
		set(YX.m2, YX.m2[, which(stage=='U' & is.na(t.isAcute))], 'CD4t', 'UAna')	
		#
		cat(paste('\nsubset\n'))
		if('score.Y'%in%colnames(YX.m2))
			YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, FASTASampleCode, t.FASTASampleCode, score.Y, telapsed, brl, stage, CDCC, lRNA, t.isAcute, nlRNA.supp, nlRNA.nsupp, t.AnyT_T1, lRNA_T1.supp, AnyPos_T1, t.AnyPos_T1, contact, fw.up.med, t.period, w, w.i, w.in, w.t, w.tn, CD41st, CD4t, CD4a, CD4b, t.Age, t.RegionHospital  ))	
		if(!'score.Y'%in%colnames(YX.m2))
			YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, stage, CDCC, lRNA, nlRNA.supp, nlRNA.nsupp, t.isAcute, t.AnyT_T1, lRNA_T1.supp, AnyPos_T1, t.AnyPos_T1, contact, fw.up.med, t.period, CD41st, CD4t, CD4a, CD4b  ))
		YX.m2[, CD4c:=CD4b]
		gc()
		#
		#	add suppressed to 'stage'
		#
		YX.m2[, stage.orig:=stage]	
		#
		#	using given VL threshold 
		#	
		gc()
		cat(paste('\nadding lRNA.mx\n'))
		set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
		YX.m2[,lRNA.c:=NULL]
		set(df.viro, NULL, 'PosRNA', df.viro[, hivc.db.Date2numeric(PosRNA)] )
		#	stratify suppressed by max VL at all times during infection window
		if(0)
		{
			YX.m2	<- merge(YX.m2, YX.m2[, {
								lRNA.c		<- 'ART.vlNA'
								tmp			<- which(!is.na(lRNA))
								if( length(tmp) && ( max(lRNA[tmp])>lRNA.supp || length(t)>length(which(stage=='ART.started')) ) )
								{									
									lRNA.c	<- 'ART.suA.N'
								}
								if( length(tmp) &&  max(lRNA[tmp])<=lRNA.supp && length(t)<=length(which(stage=='ART.started'))  )
								{
									lRNA.c	<- 'ART.suA.Y'
								} 
								list( lRNA.c= lRNA.c )
							}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)
		}	
		if(1)
		{
			YX.m2	<- merge(YX.m2, YX.m2[, {
								lRNA.c		<- 'ART.vlNA'
								tmp			<- which(!is.na(lRNA))
								if( length(tmp) && max(lRNA[tmp])>lRNA.supp )
								{									
									tmp2		<- sum(nlRNA.nsupp, na.rm=TRUE)
									if(tmp2<=0)
										lRNA.c	<- 'ART.vlNA'
									if(tmp2>0)
										lRNA.c	<- 'ART.suA.N'
									if(any(t<lRNA_T1.supp))
										lRNA.c	<- 'ART.NotYetFirstSu'
								}
								if( length(tmp) && max(lRNA[tmp])<=lRNA.supp && (length(t)>length(which(stage=='ART.started')) || any(t<lRNA_T1.supp) ) )
								{
									lRNA.c		<- 'ART.NotYetFirstSu'
								}								
								if( length(tmp) &&  max(lRNA[tmp])<=lRNA.supp && length(t)<=length(which(stage=='ART.started')) && !any(t<lRNA_T1.supp) )
								{
									tmp2		<- sum(nlRNA.supp, na.rm=TRUE)
									if(tmp2<=0)
										lRNA.c	<- 'ART.vlNA'
									if(tmp2>0)
										lRNA.c	<- 'ART.suA.Y'
								} 
								list( lRNA.c= lRNA.c )
							}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)			
		}	
		gc()
		#	stratify suppressed by current VL
		YX.m2	<- merge(YX.m2, YX.m2[, {
											ans	<- list(t=t, lRNA.c2=lRNA.c)
											if( !all(lRNA.c=='ART.suA.Y' & stage.orig=='ART.started') )
												ans$lRNA.c2[ lRNA.c=='ART.suA.Y' ]	<- 'ART.suA.M'
											ans
										}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient','t'), all.x=TRUE)
		#	stratify suppressed by max VL with no other possibility (ie no other stage in window, and at least two measurements)		
		if(0)
		{
			YX.m2	<- merge(YX.m2, YX.m2[, {
								lRNA.c3			<- 'ART.vlNA'
								tmp				<- which(!is.na(lRNA))
								if( length(tmp) && ( max(lRNA[tmp])>lRNA.supp || length(t)>length(which(stage=='ART.started')) ) )
									lRNA.c3		<- 'ART.suA.N'
								if( length(tmp) &&  max(lRNA[tmp])<=lRNA.supp && length(t)<=length(which(stage=='ART.started'))  )
								{
									tmp2		<- which( df.viro[['Patient']]==t.Patient[1] & df.viro[['PosRNA']]>=min(t[tmp]) & df.viro[['PosRNA']]<max(t[tmp])+0.125 )
									if(length(tmp2)<=1)
										lRNA.c3	<- 'ART.suA.Y1'
									if(length(tmp2)>1)
										lRNA.c3	<- 'ART.suA.Y2'
								} 
								list( lRNA.c3=lRNA.c3 )
							}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)		
		}
		if(0)
		{
			YX.m2	<- merge(YX.m2, YX.m2[, {
								lRNA.c3			<- 'ART.vlNA'
								tmp				<- which(!is.na(lRNA))
								if( length(tmp) && ( max(lRNA[tmp])>lRNA.supp || length(t)>length(which(stage=='ART.started')) ) )
								{
									tmp2		<- sum(nlRNA.nsupp, na.rm=TRUE)
									if(tmp2<=0)
										lRNA.c3	<- 'ART.vlNA'
									if(tmp2>0)
										lRNA.c3	<- 'ART.suA.N'									
								}							
								if( length(tmp) &&  max(lRNA[tmp])<=lRNA.supp && length(t)<=length(which(stage=='ART.started'))  )
								{
									tmp2		<- sum(nlRNA.supp, na.rm=TRUE)
									if(tmp2<=0)
										lRNA.c3	<- 'ART.vlNA'
									if(tmp2==1)
										lRNA.c3	<- 'ART.suA.Y1'
									if(tmp2>1)
										lRNA.c3	<- 'ART.suA.Y2'															
								} 
								if( length(tmp) &&  max(lRNA[tmp])<=lRNA.supp && length(t)>length(which(stage=='ART.started'))  )
								{
									lRNA.c3	<- 'ART.mixed'
								}
								list( lRNA.c3=lRNA.c3 )
							}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)	
		}
		if(1)
		{
			YX.m2	<- merge(YX.m2, YX.m2[, {
								lRNA.c3			<- 'ART.vlNA'
								tmp				<- which(!is.na(lRNA))
								if( length(tmp) && max(lRNA[tmp])>lRNA.supp  )
								{
									tmp2		<- sum(nlRNA.nsupp, na.rm=TRUE)
									if(tmp2<=0)
										lRNA.c3	<- 'ART.vlNA'
									if(tmp2>0)
										lRNA.c3	<- 'ART.suA.N'	
									if(any(t<lRNA_T1.supp))
										lRNA.c3	<- 'ART.NotYetFirstSu'
								}
								if( length(tmp) && max(lRNA[tmp])<=lRNA.supp && (length(t)>length(which(stage=='ART.started')) || any(t<lRNA_T1.supp) ) )
								{
									lRNA.c3		<- 'ART.NotYetFirstSu'
								}
								if( length(tmp) &&  max(lRNA[tmp])<=lRNA.supp && length(t)<=length(which(stage=='ART.started')) && !any(t<lRNA_T1.supp) )
								{
									tmp2		<- sum(nlRNA.supp, na.rm=TRUE)
									if(tmp2<=0)
										lRNA.c3	<- 'ART.vlNA'
									if(tmp2==1)
										lRNA.c3	<- 'ART.suA.Y1'
									if(tmp2>1)
										lRNA.c3	<- 'ART.suA.Y2'															
								} 						
								list( lRNA.c3=lRNA.c3 )
							}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)	
		}
		if(0)
		{
			YX.m2	<- merge(YX.m2, YX.m2[, {
								lRNA.c3			<- 'ART.vlNA'
								tmp				<- which(!is.na(lRNA))
								if( length(tmp) && ( max(lRNA[tmp])>lRNA.supp || length(t)>length(which(stage=='ART.started')) ) )
								{
									tmp2		<- which( df.viro[['Patient']]==t.Patient[1] & df.viro[['PosRNA']]>=min(t[tmp]) & df.viro[['PosRNA']]<max(t[tmp])+0.125 )
									tmp2		<- length(tmp2) / (diff(range(t))+0.125)
									if(tmp2<=0)
										lRNA.c3	<- 'ART.vlNA'
									if(tmp2>0)
										lRNA.c3	<- 'ART.suA.N'
								}							
								if( length(tmp) &&  max(lRNA[tmp])<=lRNA.supp && length(t)<=length(which(stage=='ART.started'))  )
								{
									tmp2		<- which( df.viro[['Patient']]==t.Patient[1] & df.viro[['PosRNA']]>=min(t[tmp]) & df.viro[['PosRNA']]<max(t[tmp])+0.125 )
									tmp2		<- length(tmp2) / (diff(range(t))+0.125)
									if(tmp2<=0)
										lRNA.c3	<- 'ART.vlNA'
									if(tmp2>0 & tmp2<2)
										lRNA.c3	<- 'ART.suA.Y1'
									if(tmp2>=2 & tmp2<5)
										lRNA.c3	<- 'ART.suA.Y24'
									if(tmp2>=5)
										lRNA.c3	<- 'ART.suA.Y5'							
								} 
								list( lRNA.c3=lRNA.c3 )
							}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)	
		}
		if(0)
		{
			YX.m2	<- merge(YX.m2, YX.m2[, {
								lRNA.c3			<- 'ART.vlNA'
								tmp				<- which(!is.na(lRNA))
								if( length(tmp) && ( max(lRNA[tmp])>lRNA.supp || length(t)>length(which(stage=='ART.started')) ) )
								{
									tmp2		<- length(which( df.viro[['Patient']]==t.Patient[1] & df.viro[['PosRNA']]>=(AnyPos_T1-1) & df.viro[['PosRNA']]<AnyPos_T1 ))
									if(tmp2<=0)
										lRNA.c3	<- 'ART.vlNA'
									if(tmp2>0)
										lRNA.c3	<- 'ART.suA.N'
								}							
								if( length(tmp) &&  max(lRNA[tmp])<=lRNA.supp && length(t)<=length(which(stage=='ART.started'))  )
								{
									tmp2		<- length(which( df.viro[['Patient']]==t.Patient[1] & df.viro[['PosRNA']]>=(AnyPos_T1-1) & df.viro[['PosRNA']]<AnyPos_T1 ))								
									if(tmp2<=0)
										lRNA.c3	<- 'ART.vlNA'
									if(tmp2>0 & tmp2<2)
										lRNA.c3	<- 'ART.suA.Y1'
									if(tmp2>=2 & tmp2<3)
										lRNA.c3	<- 'ART.suA.Y2'
									if(tmp2>=3)
										lRNA.c3	<- 'ART.suA.Y3'
								} 
								list( lRNA.c3=lRNA.c3 )
							}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)
		}
		gc()	
		#	set CD41st and CD4t
		tmp		<- YX.m2[, which(stage.orig=='ART.started')]
		set(YX.m2, tmp, 'CD41st', YX.m2[tmp, lRNA.c])		
		set(YX.m2, tmp, 'CD4t', YX.m2[tmp, lRNA.c])		
		set(YX.m2, tmp, 'CD4a', YX.m2[tmp, lRNA.c2])				
		set(YX.m2, tmp, 'CD4b', YX.m2[tmp, lRNA.c])		
		set(YX.m2, tmp, 'CD4c', YX.m2[tmp, lRNA.c3])
		YX.m2[, lRNA.c:=NULL]
		YX.m2[, lRNA.c2:=NULL]
		YX.m2[, lRNA.c3:=NULL]
		#	set Lost
		#stopifnot( nrow(subset(YX.m2, stage.orig!='U' & is.na(contact)))==0 )
		#tmp		<- YX.m2[, which(stage.orig!='U' & contact!='Yes')]
		#set(YX.m2, tmp, c('CD41st','CD4t','CD4a','CD4b','CD4c'), 'Lost')
		set(YX.m2, NULL, 'CD41st', YX.m2[, factor(as.character(CD41st))])
		set(YX.m2, NULL, 'CD4t', YX.m2[, factor(as.character(CD4t))])
		set(YX.m2, NULL, 'CD4a', YX.m2[, factor(as.character(CD4a))])
		set(YX.m2, NULL, 'CD4b', YX.m2[, factor(as.character(CD4b))])
		set(YX.m2, NULL, 'CD4c', YX.m2[, factor(as.character(CD4c))])		
		YX.m2[, stage.orig:=NULL]		
		#
		#	add tperiod
		#
		if('t.period'%in%colnames(YX.m2))
		{
			cat(paste('\nadding CD4t.tperiod\n'))	
			YX.m2[, CD41st.tperiod:= paste(CD41st, t.period,sep='.')]
			set(YX.m2, NULL, 'CD41st.tperiod', YX.m2[, factor(as.character(CD41st.tperiod))])
			YX.m2[, CD4t.tperiod:= paste(CD4t, t.period,sep='.')]
			set(YX.m2, NULL, 'CD4t.tperiod', YX.m2[, factor(as.character(CD4t.tperiod))])
			YX.m2[, CD4a.tperiod:= paste(CD4a, t.period,sep='.')]
			set(YX.m2, NULL, 'CD4a.tperiod', YX.m2[, factor(as.character(CD4a.tperiod))])
			YX.m2[, CD4b.tperiod:= paste(CD4b, t.period,sep='.')]
			set(YX.m2, NULL, 'CD4b.tperiod', YX.m2[, factor(as.character(CD4b.tperiod))])
			YX.m2[, CD4c.tperiod:= paste(CD4c, t.period,sep='.')]
			set(YX.m2, NULL, 'CD4c.tperiod', YX.m2[, factor(as.character(CD4c.tperiod))])
		}
		gc()
		#YX.m2[, list(LkL=median(score.Y*w.tn), n=length(score.Y)), by='CD4c']
		#YX.m2[, list(LkL=median(score.Y*w.tn), n=length(score.Y)), by='CD4c.tperiod']				
		if(!is.na(save.file))
		{
			cat(paste('\nsave YX.m2 to file', save.file))
			save(YX.m2, file=save.file)			
		}
	}
	YX.m2
}
######################################################################################
project.athena.Fisheretal.YX.model1<- function(YX, vl.suppressed=log10(1e3), cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('Diag<=350','Diag<=550','Diag>550'))
{
	YX.m1	<- copy(YX)
	#	set Y.score to min of Y.score and score.U
	#setkey(YX.m1, t, t.Patient, Patient)
	#tmp	<- YX.m1[, list(score.Y= min(score.Y, U.score, na.rm=TRUE)), by=c('t','t.Patient','Patient')]
	#set(YX.m1, NULL, 'score.Y', tmp[, score.Y])
	YX.m1[, U.score:=NULL]
	#	get Diag<350, Diag<550, Diag>550 Diag.NA
	tmp	<- YX.m1[, which(stage=='Diag' & is.na(CD4))]
	cat(paste('\nnumber entries with Diag and is.na(CD4), n=', length(tmp)))
	set(YX.m1, tmp, 'stage', 'Diag.NA')
	for(i in seq_along(cd4.cut)[-1])
	{
		tmp	<- YX.m1[, which(stage=='Diag' & !is.na(CD4) & CD4<=cd4.cut[i] & CD4>cd4.cut[i-1])]
		cat(paste('\nnumber entries with Diag and CD4 ',cd4.cut[i-1], cd4.cut[i],' , n=', length(tmp)))
		set(YX.m1, tmp, 'stage', cd4.label[i-1])
	}
	cat(paste('\nnumber entries with Diag (should be zero), n=',nrow(subset(YX.m1, stage=='Diag'))))
	#	get ART first time suppressed if !is.na(vl.suppressed)
	if(!is.na(vl.suppressed))
	{
		tmp	<- YX.m1[, which(stage=='ART.started' & is.na(lRNA))]
		cat(paste('\nnumber entries with ART.started and is.na(lRNA), n=', length(tmp)))
		set(YX.m1, tmp, 'stage', 'ART.su.NA')
		
		tmp	<- YX.m1[, which(stage=='ART.started' & !is.na(lRNA) & lRNA<=vl.suppressed)]
		cat(paste('\nnumber entries with ART.started and lRNA<=vl.suppressed, n=', length(tmp)))
		set(YX.m1, tmp, 'stage', 'ART.su.Y')
		
		tmp	<- YX.m1[, which(stage=='ART.started' & !is.na(lRNA) & lRNA>vl.suppressed)]
		cat(paste('\nnumber entries with ART.started and lRNA<vl.suppressed, n=', length(tmp)))
		set(YX.m1, tmp, 'stage', 'ART.su.N')
		
	}
	#
	YX.m1	<- subset(YX.m1, select=c(t, t.Patient, Patient, t.FASTASampleCode, FASTASampleCode, cluster, score.Y, stage, lRNA, t.period, RegionHospital, w ))
	set(YX.m1, NULL, 'stage', YX.m1[, as.factor(stage)])
	YX.m1[, score.Y.orig:= score.Y]
	set(YX.m2, YX.m2[,which(score.Y>0.99999)], 'score.Y', subset(YX, score.Y<1)[, max(c(score.Y,0.99))])
	set(YX.m2, YX.m2[,which(score.Y<1e-5)], 'score.Y', subset(YX, score.Y>0)[, min(c(score.Y,0.01))])
	YX.m1[, logit.Y:= YX.m1[, log( score.Y/(1-score.Y) )]]	
	#	add info
	print(YX.m1[, table(stage)])	
	tmp			<- YX.m1[, length(unique(Patient))]
	YX.m1.info	<- YX.m1[, list(N.notadj= length(w), N.eff=sum(w), Prop=sum(w)/tmp, E.logit.Y=weighted.mean(logit.Y, w=w), E.score.Y=weighted.mean(score.Y, w=w)), by='stage']
	print(YX.m1.info)	
	YX.m1.info	<- merge(YX.m1, YX.m1.info, by='stage')

	#	try basic beta regression
	if(0)
	{
		require(betareg)		
				
		data("GasolineYield", package = "betareg")
		GasolineYield	<- as.data.table(GasolineYield)
		GasolineYield[, yield.logit:= log(yield/(1-yield))]
		gy_logit 		<- betareg(yield ~ batch + temp - 1, data = GasolineYield)
		summary(gy_logit)
		plot(GasolineYield$temp, GasolineYield$yield)				
		tmp				<- subset(GasolineYield,select=c(temp, batch))
		set(tmp, NULL, 'batch', '6')
		setkey(tmp, temp)
		lines(tmp$temp, predict(gy_logit, tmp, type='response'), col='red', pch=18)
		lines(tmp$temp, predict(gy_logit, tmp, type='quantile', at=0.025), col='blue', pch=18)
		lines(tmp$temp, predict(gy_logit, tmp, type='quantile', at=0.975), col='blue', pch=18)
		#	betareg does not have predict, type='confidence'
		#plot(GasolineYield$temp, GasolineYield$yield)
		#lines(tmp$temp, predict(gy_logit, tmp, type='response'), col='red', pch=18)
		lines(tmp$temp, predict(gy_logit, tmp, type='response') + 2*sqrt( predict(gy_logit, tmp, type='variance') ), col='green' )
		lines(tmp$temp, predict(gy_logit, tmp, type='response') - 2*sqrt( predict(gy_logit, tmp, type='variance') ), col='green' )
		#	+- 2*std deviation in predicted response is similar to quantiles
			
		start		<- c( gy_logit$coef$mean + head( 2*sqrt(diag(vcov(gy_logit))), length(gy_logit$coef$mean)), gy_logit$coef$precision)
		gy_logit.s	<- betareg(yield ~ batch + temp - 1, data = GasolineYield, start=start, hessian=TRUE, maxit=0)
		lines(tmp$temp, predict(gy_logit.s, tmp, type='response'), col='pink' )
		start		<- c( gy_logit$coef$mean - head( 2*sqrt(diag(vcov(gy_logit))), length(gy_logit$coef$mean)), gy_logit$coef$precision)
		gy_logit.s	<- betareg(yield ~ batch + temp - 1, data = GasolineYield, start=start, hessian=TRUE, maxit=0)
		lines(tmp$temp, predict(gy_logit.s, tmp, type='response'), col='pink' )
		#	+- 2*std deviation in coefficients is broader to quantiles			
	}
	#	try inflated beta regression at one
	if(1)
	{
		require(gamlss)
		require(gamlss.dist)
		#	prepare data
		YX.m1[, score.Y:= score.Y.orig]
		YX.m1[, logit.Y:= NULL]
		set(YX.m1, YX.m1[, which(score.Y>0.999)], 'score.Y', 1)		
		tmp			<- data.table(stage= YX.m1[, levels(stage)], col=sapply( brewer.pal(YX.m1[, nlevels(stage)], 'Set1'), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m1.plot	<- merge(YX.m1, tmp, by='stage')
		
		show.link( BEZI )
		#	need log link on nu since nu is Prob(score==1)/ ( 1-Prob(score==1)  ) 
		YX.m1.fit3b	<- gamlss(formula= score.Y ~ stage-1, nu.formula=~stage-1, family=BEOI(mu.link='logit', nu.link='log'), weights=w, data=subset(YX.m1.plot, select=c(score.Y, stage, w)))
		#	this estimates the coefficients for mu and nu separately
		summary(YX.m1.fit3b)
		#	suggests for mu term coeff stageART.su.N, stageART.su.NA, stageDiag.NA, stageDiag<=350  not significantly different from 0
		#	suggests for nu term coeff stageART.su.NA  not significantly different from 0
		gamlss(formula= prop ~ cov.factor-1, nu.formula=~cov.factor-1, family=BEOI(mu.link='logit', nu.link='log'), data=X)
		#	I would like to have the SAME coefficients in mu and nu
	}
	if(1)
	{
		require(lmtest)
		require(betareg)
		require(gamlss)			
		tmp			<- data.table(stage= YX.m1[, levels(stage)], col=sapply( brewer.pal(YX.m1[, nlevels(stage)], 'Set1'), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m1.plot	<- merge(YX.m1, tmp, by='stage')
		
		YX.m1.fit0		<- lm(logit.Y ~ stage-1, weights=w, data = YX.m1.plot)
		summary(YX.m1.fit0)
		plot(YX.m1.fit0)
		bptest(YX.m1.fit0)
		#	suggests coeff stageART.su.NA  stageART.su.Y 	not significantly different from 0
		#	residuals clearly not normally distributed
		#	studentized Breusch-Pagan test clearly shows the scores are heteroscedastic
		#
		#	move to BETA regression
		#
		YX.m1.fit1 		<- betareg(score.Y ~ stage, link='logit', data = YX.m1.plot)
		summary(YX.m1.fit1)
		#	one factor is dropped because intercept is present
		YX.m1.fit2 		<- betareg(score.Y ~ stage-1, link='logit', data = YX.m1.plot)
		summary(YX.m1.fit2)
		#	suggests coeff stageART.su.N not significantly different from 0
		#	pseudo R2 = 0.1
		YX.m1.fit3 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m1.plot)
		summary(YX.m1.fit3)
		#	suggests coeff stageART.su.N, stageART.su.NA, stageART.su.Y not significantly different from 0
		#	pseudo R2 = 0.1				
		plot( YX.m1.fit3 )
		#	Cook's D suggests that the largest unexplained variation is in those Diag / not yet on treatment + those on ART
		plot( seq_len(nrow(YX.m1.plot)), YX.m1.plot[, score.Y], pch=18, col= YX.m1.plot[, col], cex=YX.m1.plot[, w^0.4])
		tmp				<- subset(YX.m1.plot, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m1.fit3, tmp, type='response'))
		lines(seq_len(nrow(tmp)), predict(YX.m1.fit3, tmp, type='quantile', at=0.5), lty=2, col='red')
		lines(seq_len(nrow(tmp)), predict(YX.m1.fit3, tmp, type='quantile', at=0.025), lty=2)
		lines(seq_len(nrow(tmp)), predict(YX.m1.fit3, tmp, type='quantile', at=0.975), lty=2)
		#lines(seq_len(nrow(tmp)), predict(YX.m1.fit3, tmp, type='response') + 2*sqrt( predict(YX.m1.fit3, tmp, type='variance') ), lty=3)
		#lines(seq_len(nrow(tmp)), predict(YX.m1.fit3, tmp, type='response') - 2*sqrt( predict(YX.m1.fit3, tmp, type='variance') ), lty=3)		
		start			<- c( YX.m1.fit3$coef$mean + head( 2*sqrt(diag(vcov(YX.m1.fit3))), length(YX.m1.fit3$coef$mean)), YX.m1.fit3$coef$precision)
		YX.m1.fit3.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m1.plot, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m1.fit3$coef$mean - head( 2*sqrt(diag(vcov(YX.m1.fit3))), length(YX.m1.fit3$coef$mean)), YX.m1.fit3$coef$precision)
		YX.m1.fit3.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m1.plot, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
					c( predict(YX.m1.fit3.sup, tmp, type='response'), rev(predict(YX.m1.fit3.slw, tmp, type='response'))), 
					border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m1[, levels(stage)], fill=YX.m1.plot[, unique(col)])
		#
		YX.m1.fit3.or	<- unique(subset(YX.m1.plot, select=stage))
		YX.m1.fit3.or	<- cbind(YX.m1.fit3.or, data.table(p=predict(YX.m1.fit3, YX.m1.fit3.or, type='response')))
		tmp				<- subset(YX.m1.fit3.or, stage=='U')[, p/(1-p) ] 		
		YX.m1.fit3.or[, or:=YX.m1.fit3.or[, p/(1-p)/tmp]]		
		#	
		YX.m1.fit4 		<- betareg(score.Y ~ stage-1 | stage, link='logit', weights=w, data = YX.m1.plot)
		summary(YX.m1.fit4)
		#	suggests that for the precision term, none of the coeff are signif non-zero 
			
	}
	
	#	plot log odds
	if(1)
	{
		require(RColorBrewer)
		#YX.m1		<- subset(YX.m1, score.Y<0.99)
		tmp			<- data.table(stage= YX.m1[, levels(stage)], col=sapply( brewer.pal(YX.m1[, nlevels(stage)], 'Set1'), my.fade.col, alpha=0.5) )
		YX.m1.plot	<- merge(YX.m1.info, tmp, by='stage')
		
		#	across time
		plot( YX.m1.plot[, t], YX.m1.plot[, logit.Y], pch=18, col= YX.m1.plot[, col], cex=YX.m1.plot[, w^0.4])
		legend('topleft', bty='n', border=NA, legend= tmp[,stage], fill=tmp[,col])
		#	by strata
		plot( seq_len(nrow(YX.m1.plot)), YX.m1.plot[, logit.Y], pch=18, col= YX.m1.plot[, col], cex=YX.m1.plot[, w^0.4])
		lines( seq_len(nrow(YX.m1.plot)), YX.m1.plot[, E.logit.Y], type='s' )
		abline(h=YX.m1.plot[,weighted.mean(logit.Y, w)], lty=2)
		legend('topleft', bty='n', border=NA, legend= tmp[,stage], fill=tmp[,col])
	}
	#	debug
	if(0)
	{
		tmp	<- subset( YX.m1.plot, stage=='ART.su.Y' & score.Y>0.9, select=c(stage, t, t.Patient, Patient, t.FASTASampleCode, FASTASampleCode, cluster, score.Y, w) )
		setkey(tmp, cluster)
		print(tmp, n=300)
	}
	YX.m1.info
}
######################################################################################
project.athena.Fisheretal.X.incare.check<- function(X.incare)
{
	#set stage: ART.interrupted
	cat(paste('\ncheck\tany(is.na(ART.interrupted))=',X.incare[, any(is.na(ART.I))]))
	cat(paste('\ncheck\tany(is.na(stage))=',X.incare[, any(is.na(stage))]))
	cat(paste('\ncheck\tany(Diag and Interrupted)=',nrow(subset(X.incare, stage=='Diag' & ART.I=='Yes'))))
	check	<- subset(X.incare, stage=='Diag' & !is.na(lRNAc) & lRNAc<3, select=c(t.Patient, t, lRNAc, stage))		#may include those in acute phase
	cat(paste('\ncheck\tpot transmitters with VL<1e3 in diagnosis=',check[,length(unique(t.Patient))]))
	tmp		<- subset(X.incare, stage=='ART.started')[, list(ART.t1=min(t)) , by='t.Patient']
	tmp		<- merge(tmp, X.incare, by='t.Patient')
	check	<- merge(unique(subset(check,select=t.Patient)), tmp, by='t.Patient')	
	check	<- subset(check,  t+1>=ART.t1 & t-0.25<=ART.t1 & (stage=='ART.started' | (stage=='Diag' & !is.na(lRNA))))
	print(check, n=3000)		
}
######################################################################################
project.athena.Fisheretal.X.CDCC<- function(X.incare, df.tpairs, clumsm.info, t.period=0.25, t.endctime=2013.)
{
	tmp		<- subset(clumsm.info,!is.na(DateFirstEverCDCC), select=c(Patient, AnyPos_T1, DateFirstEverCDCC, DateDied))
	set(tmp, tmp[,which(is.na(DateDied))], 'DateDied', t.endctime)
	setkey(tmp, Patient)		
	tmp		<- merge( data.table(Patient=df.tpairs[, unique(t.Patient)]), unique(tmp), by='Patient')	
	set(tmp, tmp[, which(DateFirstEverCDCC<AnyPos_T1)], 'DateFirstEverCDCC', tmp[which(DateFirstEverCDCC<AnyPos_T1),AnyPos_T1])
	set(tmp, NULL, 'DateDied', tmp[, floor(DateDied) + ceiling( (DateDied%%1)*100 %/% (t.period*100) ) * t.period] )
	set(tmp, NULL, 'DateFirstEverCDCC', tmp[, floor(DateFirstEverCDCC) + round( (DateFirstEverCDCC%%1)*100 %/% (t.period*100) ) * t.period] )
	cat(paste('\nnumber of potential transmitters with CDC-C, n=',nrow(tmp)))	
	tmp		<- subset(tmp, DateFirstEverCDCC<DateDied)[, list(t= seq(DateFirstEverCDCC, DateDied-t.period, by=t.period), CDCC='Yes'),by='Patient']	
	setnames(tmp, 'Patient','t.Patient')	
	X.incare	<- merge(X.incare, tmp, by=c('t.Patient','t'), all.x=1)
	set(X.incare, X.incare[,which(is.na(CDCC))],'CDCC','No')
	X.incare		
}
######################################################################################
project.athena.Fisheretal.X.cd4<- function(df.tpairs, df.immu, t.period=0.25)
{
	require(gamlss)
	stopifnot(class(df.immu$Patient)=='character',class(df.tpairs$t.Patient)=='character')
	immu	<- subset( df.immu, select=c(Patient, PosCD4, CD4) )
	setnames(immu, 'Patient','t.Patient')
	immu	<- merge(immu, unique(subset(df.tpairs, select=t.Patient)), by='t.Patient')	
	set(immu, NULL, 'PosCD4', hivc.db.Date2numeric(immu[,PosCD4]))
	tmp		<- subset(immu, select=c(t.Patient, PosCD4))[, list(ts=min(PosCD4), te=max(PosCD4)), by='t.Patient']
	set(tmp, NULL, 'ts', tmp[, floor(ts) + round( (ts%%1)*100 %/% (t.period*100) ) * t.period] )
	set(tmp, NULL, 'te', tmp[, floor(te) + round( (te%%1)*100 %/% (t.period*100) ) * t.period] )
	tmp		<- tmp[, list(PosCD4= seq(ts, te, by=t.period)),by='t.Patient']	
	setkey(tmp, t.Patient)
	immu	<- lapply( immu[, unique(t.Patient)], function(x)
			{								
				z		<- subset(tmp, t.Patient==x, PosCD4)
				z2		<- subset(immu, t.Patient==x)
				z3		<- subset(z, PosCD4>= z2[,min(PosCD4)] & PosCD4<= z2[,max(PosCD4)])
				if(!nrow(z3))
					cd4.s	<- rep(z2[1,][,mean(CD4)],nrow(z))
				else
				{
					z		<- z3 
					if(nrow(z2)==1)
						cd4.s	<- rep(z2[1,][,CD4],nrow(z))
					if(nrow(z2)>1 && z2[, all(CD4==CD4[1])])
						cd4.s	<- rep(z2[1,][,CD4],nrow(z))
					if(nrow(z2)>1 && z2[, any(CD4!=CD4[1])])
					{
						cd4.d	<- ifelse(z2[, diff(range(PosCD4))<2], 1, min(15,ceiling(nrow(z2)/8))  )
						cd4.ml	<- gamlss(CD4 ~ PosCD4, data=z2, family='NO', trace = FALSE)
						cd4.m	<- gamlss(CD4 ~ bs(PosCD4, degree=cd4.d), data=z2, family='NO', trace = FALSE)
						tryCatch({
									#	gamlss fit may go wild occasionally, in this case fall back to linear interpolation		
									if(deviance(cd4.ml)-deviance(cd4.m)>10)
										cd4.s	<- predict(cd4.m, type='response', newdata=z, data=z2) 
									if(deviance(cd4.ml)-deviance(cd4.m)<=10)
										cd4.s	<- predict(cd4.ml, type='response', newdata=z, data=z2) 	
								}, 
								warning=function(w)
								{ 
									cat(paste('\nfall back to approx for patient',x))
									cd4.s	<<- approx(z2[,PosCD4], z2[,CD4], xout=z[,PosCD4]+t.period/2, rule=2)$y	
								})						
					}
				}				
				data.table(t.Patient=x, t=z[,PosCD4], CD4=cd4.s)
			})
	immu	<- do.call('rbind',immu)
	#immu	<- immu[, {
	#			z		<- subset(tmp, tt.Patient==t.Patient[1])	
	#			if(length(PosCD4)<2)
	#			{
	#				y	<- rep(NA, nrow(z))
	#				y[z[,which( PosCD4<=t+t.period )[1]]]<- CD4
	#			}
	#			else
	#				y	<- approx(PosCD4 , CD4, xout=z[,t]+t.period/2, yleft=NA_real_, yright=NA_real_, rule=2)$y
	#			list(t=z[,t], CD4=y)
	#		},by='t.Patient']
	set(immu, NULL, 'CD4', immu[, as.integer(round(CD4))])
	tmp		<- setdiff(subset(immu, CD4>1500)[,t.Patient], subset(df.immu, CD4>1500)[, Patient])
	if(length(tmp))
		cat(paste('\nWARNING: Found patients with CD4>1500 in smoothed CD4 that do not have a CD4>1500 measurement=',paste(tmp,collapse=' ')))
	immu		
}
######################################################################################
project.athena.Fisheretal.X.viro<- function(df.tpairs, df.viro, t.period=0.25, lRNA.supp=log10(51))
{
	stopifnot(class(df.viro$Patient)=='character',class(df.tpairs$t.Patient)=='character') 
	viro	<- subset( df.viro, select=c(Patient, PosRNA, lRNA) )
	setnames(viro, 'Patient','t.Patient')
	viro	<- merge(viro, unique(subset(df.tpairs, select=t.Patient)), by='t.Patient')	
	set(viro, NULL, 'PosRNA', hivc.db.Date2numeric(viro[,PosRNA]))
	tmp		<- subset(viro, select=c(t.Patient, PosRNA))[, list(ts=min(PosRNA), te=max(PosRNA)), by='t.Patient']
	#set(tmp, NULL, 'ts', tmp[, floor(ts) + round( (ts%%1)*100 %/% (t.period*100) ) * t.period] )
	set(tmp, NULL, 'ts', tmp[, floor(ts) + floor( (ts%%1)*100 %/% (t.period*100) ) * t.period] )
	set(tmp, NULL, 'te', tmp[, floor(te) + round( (te%%1)*100 %/% (t.period*100) ) * t.period] )
	tmp		<- tmp[, list(t= seq(ts, te, by=t.period)),by='t.Patient']
	
	setnames(tmp, 't.Patient', 'tt.Patient')
	setkey(tmp, tt.Patient)
	#merge(data.table(t.Patient=tmp[, unique(tt.Patient)]), viro, by='t.Patient')	
	viro	<- viro[, {
						z		<- subset(tmp, tt.Patient==t.Patient[1])						#fails if t.Patient is factor because t.Patient[1] gets converted to integer
						yt		<- PosRNA[ which(lRNA<lRNA.supp) ]
						if(length(PosRNA)<2)
						{
							y	<- rep(NA, nrow(z))
							y[z[,which( PosRNA<=t+t.period )[1]]]<- lRNA
							y2	<- y
							ynl	<- ynu <- rep(0, nrow(z))
						}
						else
						{
							y	<- approx(PosRNA , lRNA, xout=z[,t]+t.period/2, rule=2, method='linear')$y
							y2	<- approx(PosRNA , lRNA, xout=z[,t]+t.period/2, rule=2, method='constant')$y
							ynl	<- table(cut(PosRNA[which(lRNA<lRNA.supp)], breaks=c(z[,t], z[,max(t)]+t.period), right=FALSE, labels=z[,t]))
							ynl	<- as.numeric(ynl)							
							ynu	<- table(cut(PosRNA[which(lRNA>=lRNA.supp)], breaks=c(z[,t], z[,max(t)]+t.period), right=FALSE, labels=z[,t]))							
							ynu	<- as.numeric(ynu)	
														
						}						
						list(t=z[,t], lRNA=y, lRNAc=y2, nlRNA.supp=ynl, nlRNA.nsupp=ynu, lRNA_T1.supp=ifelse(length(yt), yt[1], NA_real_))
					},by='t.Patient']		
	viro
}
######################################################################################
project.athena.Fisheretal.Y.infectiontime<- function(YX.tpairs, df.all, predict.t2inf, t2inf.args, t.period=0.25, ts.min=1980, score.set.value=1, method.infectiontime='for.infected', method.minLowerUWithNegT=TRUE)
{
	#ts.min=1980; score.min=0.01; score.set.value=1
	stopifnot(method.infectiontime%in%c('for.infected','for.transmitter'))
	cat(paste('\nAt Y.infectiontime for ',method.infectiontime,'. Using option method.minLowerUWithNegT=',method.minLowerUWithNegT))
	if(method.infectiontime=='for.infected')
	{
		b4care	<- merge(unique(subset( YX.tpairs, select=Patient )), unique(subset(df.all, select=c(Patient, DateBorn, AnyPos_T1, isAcute, NegT, PosCD4_T1, CD4_T1, AnyT_T1))), by='Patient')
	}		
	if(method.infectiontime=='for.transmitter')
	{		
		b4care	<- merge(data.table(Patient=YX.tpairs[, unique(t.Patient)]), unique(subset(df.all, select=c(Patient, DateBorn, AnyPos_T1, isAcute, NegT, PosCD4_T1, CD4_T1, AnyT_T1))), by='Patient')				
	}		
	if(method.minLowerUWithNegT)
	{		
		tmp		<- b4care[, list(AnyPos_a=AnyPos_T1-DateBorn, ts=ifelse(is.na(NegT), AnyPos_T1-10, NegT), te=AnyPos_T1), by='Patient']	#the cut off below is for NA below the minQLowerU quantile, which is calibrated to 1 year. So it is either NegT or the quantile, whichever is shorter
	}
	if(!method.minLowerUWithNegT)
	{
		tmp		<- b4care[, list(AnyPos_a=AnyPos_T1-DateBorn, ts=AnyPos_T1-10, te=AnyPos_T1), by='Patient']	#the cut off below is for NA below the minQLowerU quantile, which is calibrated to 1 year. So it is either NegT or the quantile, whichever is shorter
	}
	#	
	b4care	<- merge(b4care, tmp, by='Patient')
	#check if ts before 1980 and if so clip
	set(b4care, b4care[,which(ts<ts.min)], 'ts', ts.min )	
	set(b4care, NULL, 'ts', b4care[, floor(ts) + round( (ts%%1)*100 %/% (t.period*100) ) * t.period ] )
	set(b4care, NULL, 'te', b4care[, floor(te) + round( (te%%1)*100 %/% (t.period*100) ) * t.period - t.period] )	#last time undiagnosed is just before first time period diagnosed
	#apply predict.t2inf	
	t2inf	<- predict.t2inf(seq(0,10,t.period)*365, b4care, t2inf.args)
	t2inf	<- merge( subset( t2inf, !is.na(score)), subset(b4care, select=c(Patient,ts,te)), by='Patient' )
	set(t2inf, NULL, 'q', t2inf[,te-q/365])
	t2inf	<- subset(t2inf, ts<=q, c(Patient, q, score))
	b4care	<- merge(t2inf, subset(b4care, select=c(Patient, AnyPos_T1, AnyPos_a, isAcute)), by='Patient')
	if(!is.na(score.set.value))
		set(b4care, NULL, 'score', score.set.value)
	set(b4care, NULL, 'score', b4care[, round(score, d=3)])
	#
	if(method.infectiontime=='for.infected')
	{
		setnames(b4care, c('q','score'), c('t','score.Inf'))
		cat(paste('\nReturn infection times for #infected patients=',b4care[, length(unique(Patient))]))
	}			
	if(method.infectiontime=='for.transmitter')
	{
		setnames(b4care, c('Patient','q','score'), c('t.Patient','t','U.score'))
		cat(paste('\nReturn infection times for #transmitting patients=',b4care[, length(unique(t.Patient))]))
	}
	b4care			
}
######################################################################################
project.athena.Fisheretal.Y.coal<- function(YX.tpairs, df.all, Y.U, cluphy, cluphy.info, cluphy.map.nodectime, coal.t.Uscore.min=0.01, coal.within.inf.grace= 0.25, method.minLowerUWithNegT=TRUE, t.period= 0.25, save.file=NA, resume=0 )
{
	#coal.t.Uscore.min=0.01; coal.within.inf.grace= 0.25
	#YX.tpairs	<- subset(X.pt, select=c(t.Patient, Patient, cluster, FASTASampleCode, t.FASTASampleCode))
	#setkey(YX.tpairs, FASTASampleCode, t.FASTASampleCode)
	#YX.tpairs	<- unique(YX.tpairs)
	#
	if(resume & !is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(!resume || is.na(save.file) || inherits(readAttempt, "try-error"))
	{
		setkey(YX.tpairs, cluster)	
		#	select tpairs for which dated phylogenies are available
		tmp					<- na.omit( setdiff( YX.tpairs[,unique(cluster)], cluphy.info[, unique(cluster)] ) )
		cat(paste('\nnumber of clusters missing for analysis of pot transmitters, n=',length(tmp)))		
		df.tpairs.mrca		<- merge( YX.tpairs, unique(subset(cluphy.info, select=cluster)), by='cluster' )
		cat(paste('\nnumber of pot transmitters for which dated phylogenies are available, n=',df.tpairs.mrca[,length(unique(t.Patient))]))		
		#	compute mrcas of tpairs 
		tmp					<- df.tpairs.mrca[,	list(node=hivc.clu.mrca(cluphy, c(FASTASampleCode, t.FASTASampleCode))$mrca), by=c('FASTASampleCode','t.FASTASampleCode')]
		#if(nrow(tmp)!=nrow(df.tpairs.mrca))	stop('unexpected length of tmp')
		df.tpairs.mrca		<- merge(df.tpairs.mrca, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))		
		#	compute posterior prob that coalescence is after last NegT of transmitter or after 80% cutoff from U.score 
		#			posterior prob that coalescence is before AnyPos_T1 of individual in denominator population
		tmp					<- unique( subset( df.all, select=c(Patient, NegT)) )
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		df.tpairs.mrca		<- merge(df.tpairs.mrca, tmp, by='t.Patient')		
		
		stopifnot(setequal( df.tpairs.mrca[, unique(t.Patient)], Y.U[, unique(t.Patient)])==TRUE)
		tmp	<- subset(Y.U, U.score>coal.t.Uscore.min, c(t.Patient, t, U.score))		
		stopifnot(tmp[, any(is.na(U.score))]==FALSE)
		setkey(tmp, t.Patient, t)		
		tmp					<- tmp[, list(t.UT= t[1]), by='t.Patient']		
		df.tpairs.mrca		<- merge(df.tpairs.mrca, tmp, by='t.Patient')
		stopifnot(df.tpairs.mrca[, any(is.na(t.UT))]==FALSE)
		if(method.minLowerUWithNegT)
			tmp				<- df.tpairs.mrca[, list(t.queryT= max(t.NegT, t.UT, na.rm=TRUE)), by=c('FASTASampleCode','t.FASTASampleCode')]
		if(!method.minLowerUWithNegT)
			tmp				<- df.tpairs.mrca[, list(t.queryT= t.UT), by=c('FASTASampleCode','t.FASTASampleCode')]
		if(tmp[, any(is.na(t.queryT))])
		{
			print(subset(df.tpairs.mrca, is.na(t.queryT)))
			stop('XX')
		}
		df.tpairs.mrca		<- merge(df.tpairs.mrca, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))
		tmp					<- subset( df.tpairs.mrca, select=c(cluster, node, FASTASampleCode, t.FASTASampleCode, t.queryT) )
		tmp					<- merge(tmp, subset(df.all, select=c(FASTASampleCode, AnyPos_T1)), by='FASTASampleCode')
		set(tmp, NULL, 'AnyPos_T1', tmp[, AnyPos_T1+coal.within.inf.grace])
		tmp					<- merge( cluphy.map.nodectime, tmp, by=c('cluster','node'), allow.cartesian=TRUE)	
		coal				<- tmp[,  {
											if(length(which(!is.na(cdf)))<2)
											{
													print(c(FASTASampleCode,t.FASTASampleCode))
													print(q)
													print(cdf)
													print(t.queryT[1])
											}	
											z		<- 1-approx(q ,cdf, xout=c(t.queryT[1], AnyPos_T1[1]), yleft=0., yright=1., rule=2)$y
											#z[1] is prob survive in transmitter, z[2] is prob survive in infected
											list(coal.after.t.NegT= z[1], coal.after.i.AnyPos_T1=z[2], node=node[1])
										} , by=c('FASTASampleCode','t.FASTASampleCode')]						
		if(!is.na(save.file))
		{
			cat(paste('\nsave to file',save.file))
			save(file= save.file, coal)
		}									
	}
	cat(paste('\n number of transmitter-infected pairs, n=',nrow(coal)))
	cat(paste('\n number of potential transmitter seqs, n=',coal[, length(unique(t.FASTASampleCode))]))
	cat(paste('\n number of transmitter-infected pairs with coal !NA score (should be all), n=',nrow(subset(coal, !is.na(coal.after.t.NegT)))))
	cat(paste('\n number of transmitter-infected pairs with zero coal score, n=',nrow(subset(coal, !is.na(coal.after.t.NegT) & coal.after.t.NegT<=EPS))))	
	coal		
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
	infile.viro.study		<- paste(indircov,"/ATHENA_2013_03_Viro.R",sep='/')
	infile.immu.study		<- paste(indircov,"/ATHENA_2013_03_Immu.R",sep='/')
	infile.treatment.study	<- paste(indircov,"/ATHENA_2013_03_Regimens.R",sep='/')		
	
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
		outfile					<- paste(infile,'Ac=Y_D=0.5_sasky',sep='_')
	}
	if(0)
	{
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"				
		clu.indir				<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/gmrf_sdr06fr_-DR-RC-SH+LANL_um192rhU2080'
		clu.infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		clu.insignat			<- "Wed_Dec_18_11:37:00_2013"
		clu.infilexml.opt		<- "mph4clutx4tip"
		clu.infilexml.template	<- "um192rhU2080"	
		outfile					<- paste(infile,'Ac=Y_D=0.5_gmrf',sep='_')
	}
	#
	#	select infected individuals and return in df.select
	#
	tmp				<- project.athena.Fisheretal.select.denominator(indir, infile, insignat, indircov, infilecov, file.viro, infile.immu.study, infile.treatment.study, infiletree=infiletree, adjust.AcuteByNegT=NA, adjust.NegT4Acute=1, adjust.AcuteSelect='Yes')
	df.denom		<- tmp$df.select
	df.viro			<- tmp$df.viro
	df.immu			<- tmp$df.immu
	df.treatment	<- tmp$df.treatment
	clumsm.subtrees	<- tmp$clumsm.subtrees
	clumsm.info		<- tmp$clumsm.info
	setkey(clumsm.info, cluster)
	#
	#	select potential transmitters on MLE tree
	#
	df.tpairs		<- project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos.MinBrlMLETree(df.denom, clumsm.subtrees, any.pos.grace.yr= 0.5)	
	if(0)
	{
		tmp				<- merge( subset(df.tpairs, select=Patient), subset(clumsm.info, select=c(Patient, AnyPos_T1)), by='Patient' )
		outfile			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'nrecentlyinfected', '.pdf',sep='')
		pdf(file=outfile, w=5, h=5)
		par(mar=c(3,5,0.5,0.5))
		barplot( table( tmp[, round(AnyPos_T1)] ), ylab="# recently infected\n with unique potential transmitter" )
		dev.off()		
	}
	#
	#	get clinical data and dated phylogenies for potential transmitters
	#	
	#df.tpairs	<- subset(df.tpairs, cluster%in%c(1502,1508,1510))
	#df.tpairs.save<- copy(df.tpairs)
	#df.tpairs<- df.tpairs.save[1:100,]
	tmp						<- project.athena.Fisheretal.get.dated.phylo.for.selection(df.tpairs, clu.indir, clu.infile, clu.insignat, clu.infilexml.opt, clu.infilexml.template)	
	df.tpairs				<- tmp$df.tpairs
	cluphy					<- tmp$clu$cluphy
	cluphy.subtrees			<- tmp$clu$cluphy.subtrees
	cluphy.info				<- tmp$clu$cluphy.info
	cluphy.map.nodectime	<- tmp$clu$cluphy.map.nodectime	
	#	
	#	plot MLE tree
	#
	if(0)
	{
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
	}			
	#			
	#	plot clinical data + dated phylogenies	
	#	
	df.all					<- merge(df.all, subset(clumsm.info, select=c(FASTASampleCode, cluster)), by='FASTASampleCode', all.x=1)
	file					<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'map_clusters', '.pdf',sep='')
	project.athena.Fisheretal.plot.selected.transmitters(df.all, df.immu, df.viro, df.treatment, df.tpairs, cluphy, cluphy.info, cluphy.subtrees, cluphy.map.nodectime, file, pdf.height=800)
	#
	#	compute median height for nodes + tips 
	#
	file					<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'median_nodeheights', '.pdf',sep='')
	cluphy.nh	<- 	cluphy.map.nodectime[, {
				list(median=approx(cdf , max(q)-q, xout=0.5, yleft=NA_real_, yright=NA_real_, rule=2, method='linear')$y)				
			},by='node']
	pdf(file=file, h=3, w=5)
	par(mar=c(4,4,0.5,0.5))
	hist(cluphy.nh[,median], breaks=100, freq=0, xlab='median node height',main='',col='grey50', xlim=c(0,8.5))
	dev.off()
	quantile(cluphy.nh[,median], prob=seq(0,1,0.1))

}
######################################################################################
my.or.from.logit<- function(object, f1, f2, n1, n2, u=1.962)
{
	n1			<- round(n1)
	n2			<- round(n2)	
	coef		<- coef(object)[c(f1, f2)]
	sd.raw		<- sqrt(diag(vcov(object)))[c(f1,f2)]			
	sd.pooled	<- ifelse(n1<2 || n2<2, NA, sqrt( ((n1-1)*sd.raw[1]*sd.raw[1] 	+ (n2-1)*sd.raw[2]*sd.raw[2] )/(n1-1 + n2-1)	) )
	lor			<- -diff( coef(object)[c(f1, f2)] )
	or.ci		<- lor + c(-1,1)*u*sd.pooled
	as.double(exp(c( lor, or.ci )))
}
######################################################################################
my.rr.from.log<- function(object, f1, f2, n1, n2, u=1.962)
{
	n1			<- round(n1)
	n2			<- round(n2)
	coef		<- coef(object)[c(f1, f2)]
	sd.raw		<- sqrt(diag(vcov(object)))[c(f1,f2)]			
	#sd.pooled	<- sqrt(sum(1/c(n1, n2))) * sqrt( ((n1-1)*sd.raw[1]*sd.raw[1] 	+ (n2-1)*sd.raw[2]*sd.raw[2] )/(n1-1 + n2-1)	)
	sd.pooled	<- ifelse(n1<2 || n2<2, NA, sqrt( ((n1-1)*sd.raw[1]*sd.raw[1] 	+ (n2-1)*sd.raw[2]*sd.raw[2] )/(n1-1 + n2-1)	) )
	lrr			<- -diff( coef(object)[c(f1, f2)] )
	rr.ci		<- lrr + c(-1,1)*u*sd.pooled
	as.double(exp(c( lrr, rr.ci )))
}
######################################################################################
my.prop.from.log<- function(object, f1, u=1.962)
{
	coef		<- coef(object)[f1]
	sd			<- sqrt(diag(vcov(object)))[f1]						
	ci			<- coef + c(-1,1)*u*sd
	as.double(exp(c( coef, ci )))
}
######################################################################################
project.athena.Fisheretal.check<- function()
{
	load('~/duke/2013_HIV_NL/ATHENA_2013/data/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_Wed_Dec_18_11:37:00_2013_RICT_3ca_tATHENAclu.R')
	YX.part1hpc	<- YX.part1
	YX.part1o	<- YX.part1.save
	
	stopifnot( setequal( YX.part1o[, unique(Patient)], YX.part1hpc[, unique(Patient)] )==TRUE )
	stopifnot( setequal( YX.part1o[, unique(t.Patient)], YX.part1hpc[, unique(t.Patient)] )==TRUE )
	stopifnot( setequal( YX.part1o[, unique(FASTASampleCode)], YX.part1hpc[, unique(FASTASampleCode)] )==TRUE )
	stopifnot( setequal( YX.part1o[, unique(t.FASTASampleCode)], YX.part1hpc[, unique(t.FASTASampleCode)] )==TRUE )
	setkey(YX.part1o, FASTASampleCode, t.FASTASampleCode, t)
	setkey(YX.part1hpc, FASTASampleCode, t.FASTASampleCode, t)
	
	stopifnot( all(YX.part1o[, t]==YX.part1hpc[, t])==TRUE )
	tmp			<- colnames(YX.part1hpc)[ !colnames(YX.part1hpc)%in%c("FASTASampleCode","t.FASTASampleCode","t","Patient","t.Patient") ]
		
	dummy<- sapply(tmp, function(col)
			{
				cat(paste('\nprocess',col))
				z<- which( YX.part1o[, col, with=FALSE]!=YX.part1hpc[, col, with=FALSE] )				
				stopifnot( length(z)==0 )
			})
	
}
######################################################################################
project.athena.Fisheretal.characteristics<- function()
{
	require(data.table)
	require(ape)
	indir					<- paste(DATA,"fisheretal_data",sep='/')		
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	outdir					<- paste(DATA,"fisheretal",sep='/')
	infile.cov.study		<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infile.viro.study		<- paste(indircov,"ATHENA_2013_03_Viro.R",sep='/')
	infile.immu.study		<- paste(indircov,"ATHENA_2013_03_Immu.R",sep='/')
	infile.treatment.study	<- paste(indircov,"ATHENA_2013_03_Regimens.R",sep='/')
	infile.cov.all			<- "ATHENA_2013_03_AllSeqPatientCovariates_AllMSM"
	infile.viro.all			<- paste(indircov,"ATHENA_2013_03_Viro_AllMSM.R",sep='/')
	infile.immu.all			<- paste(indircov,"ATHENA_2013_03_Immu_AllMSM.R",sep='/')
	infile.treatment.all	<- paste(indircov,"ATHENA_2013_03_Regimens_AllMSM.R",sep='/')		
	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	resume					<- 1
	verbose					<- 1
	if(0)
	{
		method					<- '3c'
		method.nodectime		<- 'any'
		method.risk				<- 'm21st.cas'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"					
		clu.infilexml.opt		<- "mph4clutx4tip"
		clu.infilexml.template	<- "um192rhU2080"	
		outfile					<- paste(infile,'Ac=MY_D=35_gmrf',sep='_')
	}
	if(0)
	{
		method					<- '3d'
		method.nodectime		<- 'any'
		method.risk				<- 'm21st.cas'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"					
		clu.infilexml.opt		<- "mph4clutx4tip"
		clu.infilexml.template	<- "um192rhU2080"	
		outfile					<- paste(infile,'Ac=MY_D=35_gmrf',sep='_')
	}	
	if(0)
	{		
		method					<- '3d'
		method.nodectime		<- 'any'
		method.risk				<- 'm21st.cas'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		clu.infilexml.opt		<- "clrh80"
		clu.infilexml.template	<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
	}
	if(1)
	{		
		method					<- '3c'
		method.nodectime		<- 'any'
		method.risk				<- 'm21st.cas'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		clu.infilexml.opt		<- "clrh80"
		clu.infilexml.template	<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
	}
	
	
	if(method.nodectime=='any')
		method		<- paste(method,'a',sep='')
	if(method.nodectime=='map')
		method		<- paste(method,'m',sep='')	
	
	file	<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RIPDT_RIMSM_',method,'_info','.R',sep='')
	cat(paste('\nload',file))
	tmp		<- load(file)
	save(ri.PT.su, ri.allmsm.su, ri.allmsm.resolution, ri.allmsm.l2c, file=save.file)
	
	
	ri.PT.su$rin
	ri.allmsm.su$rin
	
	subset(ri.PT.su$rip, stat=='p')
	subset(ri.allmsm.su$rip, stat=='p')
	
	subset(ri.PT.su$ris, stat=='quantile_0.5')
	subset(ri.allmsm.su$ris, stat=='quantile_0.5')
	
	# Figure 1A
	tmp	<- subset(YX.save, select=c(Patient, t.Patient, score.Y))
	setkey(tmp, Patient, t.Patient)	
	tmp	<- unique(tmp)
	hist(tmp[, score.Y], breaks=100)
}
######################################################################################
project.athena.Fisheretal.compositionage<- function()
{
	tperiod.info	<- structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.653, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.408, 2008.057, 2009.512, 
							2011)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max"))
	tmp				<- 1:4
	method.risk		<- "m5.tA.tp1"
	save.file		<-  "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3da_tablesCLU_m5.tA.tp"
	tmp				<- lapply(tmp, function(tp)
						{				
							X.tables		<- project.athena.Fisheretal.estimate.risk.table(YX=NULL, X.den=NULL, X.msm=NULL, X.clu=NULL, resume=TRUE, save.file=paste(save.file, tp, '.R', sep=''), method=NULL)
							X.tables$risk.table
						})
	risk.table		<- do.call('rbind',tmp)	
	risk.table[, t.period:= substr(factor, nchar(factor), nchar(factor))]
	risk.table		<- merge(risk.table, tperiod.info, by='t.period')
	risk.table[, group:= substr(factor, 2, nchar(factor)-2) ]
	set(risk.table, risk.table[, which(group=='<=100')], 'group', '>45')
	
	# plot number potential transmission intervals
	ggplot(risk.table, aes(x=t.period, y=n, colour=group, group=group)) + labs(x="calendar time periods", y="# potential transmission intervals") +
			scale_colour_manual(name='Age group', values = colorRampPalette(brewer.pal(9, "Set1"),interpolate='spline',space="Lab")(risk.table[,length(unique(group))])) +
			geom_line() + 
			facet_grid(stat ~ ., scales='free', margins=FALSE)
	# plot proportion potential transmission intervals
	ggplot(risk.table, aes(x=t.period, y=p, colour=group, group=group)) + labs(x="calendar time periods", y="# potential transmission intervals") +
			scale_colour_manual(name='Age group', values = colorRampPalette(brewer.pal(9, "Set1"),interpolate='spline',space="Lab")(risk.table[,length(unique(group))])) +
			geom_line() + 
			facet_grid(stat ~ ., scales='free', margins=FALSE)
	# look at bias compared to X.msm
	ggplot(risk.table, aes(x=t.period, y=p, fill=stat, colour=stat)) + geom_bar(stat='identity', position='dodge')	+
			facet_grid(. ~ group, scales='free', margins=FALSE)	
}
######################################################################################
project.athena.Fisheretal.compositionclock<- function()
{	
	#
	#	load ML tree branch lengths
	#
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_data/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500_Wed_Dec_18_11:37:00_2013.R'
	load(file)	#loads ph		
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_data/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500_Wed_Dec_18_11:37:00_2013_distTips.R'
	load(file)	#loads brl	
	brl.n	<- attr(brl,'Size')
	tmp		<- unique( subset(coal, select=FASTASampleCode) )
	set(tmp, NULL, 'FASTASampleCode', tmp[, as.character(FASTASampleCode)])
	tmp		<- tmp[, list(sc.i=match(FASTASampleCode, ph$tip.label)), by='FASTASampleCode']
	coal	<- merge(coal, tmp, by='FASTASampleCode')	
	tmp		<- unique( subset(coal, select=t.FASTASampleCode) )
	set(tmp, NULL, 't.FASTASampleCode', tmp[, as.character(t.FASTASampleCode)])
	tmp		<- tmp[, list(sc.t=match(t.FASTASampleCode, ph$tip.label)), by='t.FASTASampleCode']
	coal	<- merge(coal, tmp, by='t.FASTASampleCode')	
	tmp		<- coal[, which(sc.t<sc.i)]
	tmp2	<- coal[tmp, sc.t]
	set(coal, tmp, 'sc.t', coal[tmp,sc.i])
	set(coal, tmp, 'sc.i', tmp2)	
	tmp		<- coal[, 	list( brl= brl[ my.lower.tri.index(brl.n, sc.t, sc.i) ]) , by=c('FASTASampleCode','t.FASTASampleCode')]
	coal	<- merge(coal, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))
	coal[, sc.i:=NULL]
	coal[, sc.t:=NULL]
	
	coal[, cor(scoal.after.t.NegT, brl)]
	#	-0.3377944
	coal[, cor(gcoal.after.t.NegT, brl)]
	#	-0.3793354		
}
######################################################################################
project.athena.Fisheretal.compositioncoalvsbrl<- function()
{	
	require(reshape2)
	require(data.table)
	require(ape)
	require(adephylo)
	require(caTools)
	#stop()
	indir					<- paste(DATA,"fisheretal_data",sep='/')		
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	outdir					<- paste(DATA,"fisheretal",sep='/')
	infile.cov.study		<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infile.viro.study		<- paste(indircov,"ATHENA_2013_03_Viro.R",sep='/')
	infile.immu.study		<- paste(indircov,"ATHENA_2013_03_Immu.R",sep='/')
	infile.treatment.study	<- paste(indircov,"ATHENA_2013_03_Regimens.R",sep='/')
	infile.cov.all			<- "ATHENA_2013_03_AllSeqPatientCovariates_AllMSM"
	infile.viro.all			<- paste(indircov,"ATHENA_2013_03_Viro_AllMSM.R",sep='/')
	infile.immu.all			<- paste(indircov,"ATHENA_2013_03_Immu_AllMSM.R",sep='/')
	infile.treatment.all	<- paste(indircov,"ATHENA_2013_03_Regimens_AllMSM.R",sep='/')		
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	infiletree				<- paste(infile,"examlbs500",sep="_")
	insignat				<- "Wed_Dec_18_11:37:00_2013"							
	clu.infilexml.opt		<- "clrh80"
	clu.infilexml.template	<- "sasky_sdr06fr"	
	method.nodectime		<- 'any'
	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))	
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	resume					<- 1
	verbose					<- 1
	#
	#	process SASKY
	#
	#	collect files containing dated cluster phylogenies
	files		<- list.files(indir)
	if(!length(files))	stop('no input files matching criteria')
	files		<- files[ sapply(files, function(x) grepl(infile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(clu.infilexml.opt, x, fixed=1) & grepl(clu.infilexml.template,x, fixed=1) & grepl('_cluposterior_[0-9]+',x) & grepl('R$',x) ) ]		
	if(!length(files))	stop('no input files matching criteria')
	tmp			<- regmatches( files, regexpr('_cluposterior_[0-9]+',files)) 
	cluster		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
	file.info	<- data.table(file=files, cluster=cluster)
	setkey(file.info, cluster)
	#	get dated cluster MAP phylogenies
	tmp			<- lapply(seq_len(nrow(file.info)), function(i)
			{				
				#	load dated cluster phylogenies
				file					<- paste(indir, file.info[i,file], sep='/')
				if(verbose) cat(paste('\nload file=',file,'i=',i))
				tmp						<- load(file)
				topo.map				<- mph.clu.dtopo[which.max(freq),]					
				tmp						<- which( grepl( paste('mph.i=',topo.map[,mph.i],'_',sep=''), names(ph.consensus) ) )
				if(length(tmp)!=1)	stop('unexpected ph.consensus index')
				topo.map.ph				<- ph.consensus[[ tmp ]]
				topo.map.ph$tip.label	<- sapply(strsplit(topo.map.ph$tip.label, '_'), function(x) x[[length(x)]] )
				list(ph=topo.map.ph)			
			})
	#	extract calendar time brl between all tips in clusters
	tmp			<- lapply(seq_along(tmp), function(i)
			{
				z	<- as.matrix(distTips(tmp[[i]]$ph))
				df	<- combs( sort(tmp[[i]]$ph$tip.label), 2 )
				df	<- data.table( FASTASampleCode.1=df[,1], FASTASampleCode.2=df[,2] )	
				df[, list(tl=z[FASTASampleCode.1,FASTASampleCode.2]), by=c('FASTASampleCode.1','FASTASampleCode.2')]		
			})
	df			<- do.call('rbind', tmp) 
	setnames(df, 'tl','tl.sasky')
	#
	#	process GMRF
	#					
	clu.infilexml.opt		<- "mph4clutx4tip"
	clu.infilexml.template	<- "um192rhU2080"	
	#	collect files containing dated cluster phylogenies
	files		<- list.files(indir)
	if(!length(files))	stop('no input files matching criteria')
	files		<- files[ sapply(files, function(x) grepl(infile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(clu.infilexml.opt, x, fixed=1) & grepl(clu.infilexml.template,x, fixed=1) & grepl('_cluposterior_[0-9]+',x) & grepl('R$',x) ) ]		
	if(!length(files))	stop('no input files matching criteria')
	tmp			<- regmatches( files, regexpr('_cluposterior_[0-9]+',files)) 
	cluster		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
	file.info	<- data.table(file=files, cluster=cluster)
	setkey(file.info, cluster)
	#	get dated cluster MAP phylogenies
	tmp			<- lapply(seq_len(nrow(file.info)), function(i)
			{				
				#	load dated cluster phylogenies
				file					<- paste(indir, file.info[i,file], sep='/')
				if(verbose) cat(paste('\nload file=',file,'i=',i))
				tmp						<- load(file)
				topo.map				<- mph.clu.dtopo[which.max(freq),]					
				tmp						<- which( grepl( paste('mph.i=',topo.map[,mph.i],'_',sep=''), names(ph.consensus) ) )
				if(length(tmp)!=1)	stop('unexpected ph.consensus index')
				topo.map.ph				<- ph.consensus[[ tmp ]]
				topo.map.ph$tip.label	<- sapply(strsplit(topo.map.ph$tip.label, '_'), function(x) x[[length(x)]] )
				list(ph=topo.map.ph)			
			})
	#	extract calendar time brl between all tips in clusters
	tmp			<- lapply(seq_along(tmp), function(i)
			{
				z	<- as.matrix(distTips(tmp[[i]]$ph))
				df	<- combs( sort(tmp[[i]]$ph$tip.label), 2 )
				df	<- data.table( FASTASampleCode.1=df[,1], FASTASampleCode.2=df[,2] )	
				df[, list(tl=z[FASTASampleCode.1,FASTASampleCode.2]), by=c('FASTASampleCode.1','FASTASampleCode.2')]		
			})
	tmp			<- do.call('rbind', tmp)
	setnames(tmp, 'tl','tl.gmrf')
	df			<- merge(df, tmp, by=c('FASTASampleCode.1','FASTASampleCode.2'))
	#
	#	get ML branch lengths
	#
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_data/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500_Wed_Dec_18_11:37:00_2013.R'
	load(file)	#loads ph		
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_data/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500_Wed_Dec_18_11:37:00_2013_distTips.R'
	load(file)	#loads brl	
	brl.n	<- attr(brl,'Size')
	tmp		<- unique( subset(df, select=FASTASampleCode.1) )
	tmp		<- tmp[, list(sc.i=match(FASTASampleCode.1, ph$tip.label)), by='FASTASampleCode.1']
	df		<- merge(df, tmp, by='FASTASampleCode.1')	
	tmp		<- unique( subset(df, select=FASTASampleCode.2) )
	tmp		<- tmp[, list(sc.t=match(FASTASampleCode.2, ph$tip.label)), by='FASTASampleCode.2']
	df		<- merge(df, tmp, by='FASTASampleCode.2')
	tmp		<- df[, which(sc.t<sc.i)]
	tmp2	<- df[tmp, sc.t]
	set(df, tmp, 'sc.t', df[tmp,sc.i])
	set(df, tmp, 'sc.i', tmp2)	
	tmp		<- df[, 	list( brl= brl[ my.lower.tri.index(brl.n, sc.t, sc.i) ]) , by=c('FASTASampleCode.1','FASTASampleCode.2')]
	df		<- merge(df, tmp, by=c('FASTASampleCode.1','FASTASampleCode.2'))
	df[, sc.i:=NULL]
	df[, sc.t:=NULL]
	#
	#	compare and plot
	#
	setnames(df, c('tl.sasky','tl.gmrf'),c('SA','GMRF'))
	df		<- melt(df, id.vars=c('FASTASampleCode.1','FASTASampleCode.2','brl'), value.name='ctl', variable.name='treeprior')
	set(df, df[, which(treeprior=='SA')], 'treeprior', 'serial SA birth-death model')
	set(df, df[, which(treeprior=='GMRF')], 'treeprior', 'GMRF coalescent model')
	#	compare to ML brl
	#	SAsky does not explain brl better ..
	ggplot(df, aes(x=ctl, y=brl)) + geom_point() + facet_grid(treeprior ~ ., margins=FALSE) + stat_smooth()
	g.sa	<- gamlss(brl~ctl-1, sigma.formula=~ctl, data=subset(df, treeprior=='SA'), family=GA(mu.link='identity'))	
	g.g		<- gamlss(brl~ctl-1, sigma.formula=~ctl, data=subset(df, treeprior=='GMRF'), family=GA(mu.link='identity'))
	Rsq(g.sa)
	#	0.6035012
	Rsq(g.g)
	#	0.623532
	GAIC(g.sa)
	#	-82574.98
	GAIC(g.g)
	#	-83385.81
	
	#
	#
	ggplot(df, aes(x=ctl, colour=treeprior)) + stat_ecdf() + coord_cartesian(xlim = c(0, 20)) +
			labs(x='calendar time between clustering tip pairs\n(years)', y='cumulative distribution function') +
			scale_colour_brewer(palette='Set1', name='Tree prior') +
			theme(legend.position=c(0,1), legend.justification=c(0,1) )
	file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_2011_Wed_Dec_18_11:37:00_2013_gmrf_sasky_CLEN.pdf'
	ggsave(file=file, w=4, h=6)
}
######################################################################################
project.athena.Fisheretal.compositioncoal<- function()
{	
	load( '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_gmrf_2011_Wed_Dec_18_11:37:00_2013_YXCLU3da_all.R' )
	g.coal	<- Y.coal
	set(g.coal, NULL, 'coal.after.t.NegT', g.coal[,coal.after.t.NegT-coal.after.i.AnyPos_T1])
	load( '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXCLU3da_all.R' )
	s.coal	<- Y.coal
	set(s.coal, NULL, 'coal.after.t.NegT', s.coal[,coal.after.t.NegT-coal.after.i.AnyPos_T1])
	setnames(s.coal, 'coal.after.t.NegT', 'scoal.after.t.NegT')
	setnames(g.coal, 'coal.after.t.NegT', 'gcoal.after.t.NegT')
	coal	<- merge( subset(s.coal, select=c(FASTASampleCode, t.FASTASampleCode, scoal.after.t.NegT)), subset(g.coal, select=c(FASTASampleCode, t.FASTASampleCode, gcoal.after.t.NegT)), by=c('FASTASampleCode','t.FASTASampleCode'))
	label	<- "                  Posterior probability that\n                  coalesence is compatible with\n                  direct HIV transmission"
	ggplot(coal, aes(x= scoal.after.t.NegT, y=gcoal.after.t.NegT)) + geom_point() +
		labs(x=paste(label,'\n                  under Sampled Ancestor tree prior',sep=''), y=paste(label,'\n                  under Coalescent tree prior',sep='')) +
		theme(axis.title.x=element_text(hjust=0), axis.title.y=element_text(hjust=0))
	file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_2011_Wed_Dec_18_11:37:00_2013_gmrf_sasky_COAL.pdf'
	ggsave(file=file, w=6, h=6)
	#
	#
	coal[, cor(scoal.after.t.NegT, gcoal.after.t.NegT)]		
	#	0.9618471
	coal[, table(scoal.after.t.NegT<0.5 , gcoal.after.t.NegT<0.5) ]
	
	coal[, d:= scoal.after.t.NegT-gcoal.after.t.NegT]
	coal[, d:= abs(scoal.after.t.NegT-gcoal.after.t.NegT)]
	#	anonmymize
	setkey(clumsm.info, cluster, AnyPos_T1)
	tmp						<- unique(subset(clumsm.info, select=Patient))
	set(tmp, NULL, 'PatientA', paste('P',seq_len(nrow(tmp)),sep=''))
	clumsm.info				<- merge(clumsm.info, tmp, by='Patient')
	tmp						<- clumsm.info[, list(FASTASampleCode=FASTASampleCode, FASTASampleCodeA=paste(PatientA, seq_along(FASTASampleCode), sep='-')) , by='PatientA']
	set(tmp, NULL, 'FASTASampleCodeA', tmp[, gsub('P','S',FASTASampleCodeA)])
	clumsm.info				<- merge(clumsm.info, subset(tmp, select=c(FASTASampleCode, FASTASampleCodeA)), by='FASTASampleCode')	
	tmp						<- subset(clumsm.info, select=c(cluster, Patient, FASTASampleCode, PatientA, FASTASampleCodeA))
	coal					<- merge(coal, tmp, by='FASTASampleCode')
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	coal					<- merge(coal, tmp, by='t.FASTASampleCode')
	setkey(coal, d, cluster)
	coal					<- subset(coal, select=c(t.PatientA, t.FASTASampleCodeA, PatientA, FASTASampleCodeA, scoal.after.t.NegT, gcoal.after.t.NegT, d))
	file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_2011_Wed_Dec_18_11:37:00_2013_gmrf_sasky_COAL.csv'
	write.csv(coal, file=file, eol='\r\n')
	
}
######################################################################################
project.athena.Fisheretal.compositionacute<- function()
{
	tmp<- df.all
	tmp<- df.all.allmsm
	setkey(tmp, Patient)
	tmp<- unique(tmp)
	tmp<- subset(tmp, AnyPos_T1>1996.6 & AnyPos_T1<2011, select=c(Patient, isAcute, AnyPos_T1))
	set(tmp, tmp[, which(is.na(isAcute))],'isAcute', 'Missing')
	
	ggplot(tmp, aes(x=AnyPos_T1, fill=isAcute, colour=isAcute)) + geom_bar(aes(y=..count../sum(..count..)), binwidth=1, position='dodge')
	ggplot(tmp, aes(x=AnyPos_T1, fill=isAcute, colour=isAcute)) + geom_bar(binwidth=0.5, position='dodge')
	ggplot(tmp, aes(x=AnyPos_T1, fill=isAcute, colour=isAcute)) + geom_bar(binwidth=1, position='fill')
	
	tmp<- subset(YX, select=c(Patient, t.Patient, w.i, score.Y))
	setkey(tmp, Patient, t.Patient)
	tmp<- unique(tmp)
}
######################################################################################
project.athena.Fisheretal.sensitivity.tables.m2<- function(df, factors, levels, file)
{
	#	time period	
	set(df, NULL, 'factor', df[, as.character(factor)])
	setkey(df, factor)
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	#	cascade group
	set(df, NULL, 'group', df[, paste(substr(factor, 1, 1),'total',sep='_')])
	#	reset factor order
	set(df, NULL, 'factor', df[, factor(factor, levels=levels, labels=levels)])
	setkey(df, stat, t.period, factor)
	df		<- unique(df)
	#	first col is potential transmission intervals w seq per 100 'recipient w seq' in time period
	tmp		<- subset(df, stat%in%c('X.seq','nRec'))
	tmp2	<- subset(df, stat=='X.seq', select=c(factor, t.period, p))
	tmp		<- dcast.data.table(tmp, t.period+group+factor+method.risk~stat, value.var='n')
	tmp		<- merge(tmp, tmp2, by=c('t.period','factor'))
	tmp[, n:= X.seq/nRec*100]	
	tmp2	<- tmp[, list(n=sum(n)), by=c('group','t.period')]
	tmp2	<- merge(tmp2, tmp2[, list(group=group, p=n/sum(n)), by='t.period'], by=c('t.period','group'))
	tmp2[, factor:=group]
	tmp		<- rbind(tmp, tmp2, fill=TRUE)
	set(tmp, NULL, 'p', tmp[, round(p, d=3)]*100)
	set(tmp, NULL, 'n', tmp[, round(n, d=0)])	
	set(tmp, NULL, 'Xseq.PY', tmp[, paste('=\"',n,' (',p,')\"',sep='')])
	ans		<- subset(tmp, select=c(t.period, factor,  Xseq.PY))
	
	#	second col is proportion of pt sequenced
	#	X.msm is computed against RI.SEQ, so can just divide X.seq/X.msm
	tmp		<- subset(df, stat%in%c('X.seq','X.msm'), select=c(stat, method.risk, group, factor, n, t.period))
	tmp		<- dcast.data.table(tmp, t.period+group+factor+method.risk~stat, value.var='n')	
	tmp2	<- tmp[, list(X.msm=sum(X.msm), X.seq=sum(X.seq)), by=c('group','t.period')]
	tmp2[, factor:=group]
	tmp		<- rbind(tmp, tmp2, fill=TRUE)
	set(tmp, NULL, 'p', tmp[, round(X.seq/X.msm, d=3)]*100)
	set(tmp, NULL, 'Xseq.P', tmp[, paste('=\"',p,'\"',sep='')])
	ans		<- merge(ans, subset(tmp, select=c(factor, t.period, Xseq.P)), by=c('t.period','factor'))
	
	#	third col is likely transmission intervals per 100 'recipient with lkl transmitter' in time period
	tmp		<- subset(df, stat%in%c('YX','nRecLkl'))
	tmp2	<- subset(df, stat=='YX', select=c(factor, t.period, p))
	tmp		<- dcast.data.table(tmp, t.period+group+factor+method.risk~stat, value.var='n')
	tmp		<- merge(tmp, tmp2, by=c('factor','t.period'))
	tmp[, n:= YX/nRecLkl*100]	
	tmp2	<- tmp[, list(n=sum(n)), by=c('t.period','group')]
	tmp2	<- merge(tmp2, tmp2[, list(group=group, p=n/sum(n)), by='t.period'], by=c('t.period','group'))
	tmp2[, factor:=group]	
	tmp		<- rbind(tmp, tmp2, fill=TRUE)
	set(tmp, NULL, 'p', tmp[, round(p, d=3)]*100)
	set(tmp, NULL, 'n', tmp[, round(n, d=0)])	
	set(tmp, NULL, 'YX.PY', tmp[, paste('=\"',n,' (',p,')\"',sep='')])
	ans		<- merge(ans, subset(tmp, select=c(factor, t.period, YX.PY)), by=c('t.period','factor'))
	
	#	extra col: missing only per 100 recipient w lkl transmitter
	tmp		<- subset(df, stat%in%c('Sx.e0cp','nRecLkl'))
	tmp2	<- subset(df, stat=='Sx.e0cp', select=c(factor, t.period, l95.bs, u95.bs))
	tmp		<- dcast.data.table(tmp, t.period+group+factor+method.risk~stat, value.var='n')
	tmp		<- merge(tmp, tmp2, by=c('factor','t.period'))
	tmp[, n:= Sx.e0cp/8/nRecLkl*100]
	tmp[, l95.bs:= l95.bs/8/nRecLkl*100]
	tmp[, u95.bs:= u95.bs/8/nRecLkl*100]	
	tmp2	<- tmp[, list(n=sum(n), l95.bs=sum(l95.bs), u95.bs=sum(u95.bs)), by=c('t.period','group')]
	tmp2	<- merge(tmp2, tmp2[, list(group=group), by='t.period'], by=c('t.period','group'))
	tmp2[, factor:=group]	
	tmp		<- rbind(tmp, tmp2, fill=TRUE)	
	set(tmp, NULL, 'YXm.n', tmp[, paste('=\"',round(n, d=1),' (',round(l95.bs, d=1),'-',round(u95.bs, d=1),')\"',sep='')])
	ans		<- merge(ans, subset(tmp, select=c(factor, t.period, YXm.n)), by=c('t.period','factor'))
	
	#	fourth col is likely transmission intervals + expected missing intervals per 100 'recipient with lkl transmitter' in time period
	tmp		<- subset(df, stat%in%c('YX','Sx.e0cp','nRecLkl'))
	tmp2	<- subset(df, stat=='Sx.e0cp', select=c(factor, t.period, l95.bs, u95.bs))
	tmp		<- dcast.data.table(tmp, t.period+group+factor+method.risk~stat, value.var='n')
	tmp		<- merge(tmp, tmp2, by=c('factor','t.period'))
	tmp[, n:= (YX+Sx.e0cp/8)/nRecLkl*100]
	tmp[, l95.bs:= (YX+l95.bs/8)/nRecLkl*100]
	tmp[, u95.bs:= (YX+u95.bs/8)/nRecLkl*100]	
	tmp		<- merge(tmp, tmp[, list(factor=factor, p=n/sum(n)), by='t.period'], by=c('t.period','factor'))	
	tmp2	<- tmp[, list(n=sum(n), l95.bs=sum(l95.bs), u95.bs=sum(u95.bs)), by=c('t.period','group')]
	tmp2	<- merge(tmp2, tmp2[, list(group=group, p=n/sum(n)), by='t.period'], by=c('t.period','group'))	
	tmp2[, factor:=group]
	tmp		<- rbind(tmp, tmp2, fill=TRUE)
	set(tmp, NULL, 'YXm.PY', tmp[, paste('=\"',round(n, d=1),' (',round(p, d=3)*100,')\"',sep='')])
	ans		<- merge(ans, subset(tmp, select=c(factor, t.period, YXm.PY)), by=c('t.period','factor'))
	
	#	last col is likelihood
	tmp		<- subset(df, stat%in%c('LkL'), select=c(group, factor, t.period, p))
	set(tmp, NULL, 'YX.LKL', tmp[, paste('=\"',round(p, d=1),'\"',sep='')])
	ans		<- merge(ans, subset(tmp, select=c(factor, t.period, YX.LKL)), by=c('t.period','factor'), all.x=TRUE)
	set(ans, ans[, which(is.na(YX.LKL))], 'YX.LKL', '=""')
	
	#	add legend
	ans		<- merge(ans, subset(factors, select=c(factor, factor.legend)), by='factor')
	ans		<- subset(ans, select=c(t.period, factor, Xseq.PY, Xseq.P, YX.PY, YXm.PY, YX.LKL, factor.legend))
	setkey(ans, t.period, factor)
	set(ans, NULL, 'factor', ans[,factor.legend])
	set(ans, NULL, 'factor.legend', NULL)
	setnames(ans, 	c('t.period','factor'), c('=\"Observation\"&CHAR(13)&\"period\"','=\"Treatment cascade stage\"'))
	setnames(ans, 	'Xseq.PY', '=\"Potential\"&CHAR(13)&\"transmission intervals\"&CHAR(13)&\"with a sequence per\"&CHAR(13)&\"100 recipient MSM\"&CHAR(13)&\"(person-years (%))\"')
	setnames(ans, 	'Xseq.P', '=\"Sequence\"&CHAR(13)&\"sampling\"&CHAR(13)&\"probability\"&CHAR(13)&\"(%)\"')
	setnames(ans, 	'YX.PY', '=\"Likely\"&CHAR(13)&\"transmission intervals\"&CHAR(13)&\"per 100 recipient MSM;\"&CHAR(13)&\"observed\"&CHAR(13)&\"(person-years (%))\"')
	setnames(ans, 	'YXm.PY', '=\"Likely\"&CHAR(13)&\"transmission intervals\"&CHAR(13)&\"per 100 recipient MSM;\"&CHAR(13)&\"observed + expected missing\"&CHAR(13)&\"(person-years (%))\"')
	setnames(ans, 	'YX.LKL', '=\"Likelihood of\"&CHAR(13)&\"direct HIV transmission\"&CHAR(13)&\"(median)\"')
	write.csv(ans, file=file, eol="\r\n", row.names=FALSE)
	#in Excel replace COMMA with ,
	#right click on table or Apple-1 and set wrap text
}
######################################################################################
project.athena.Fisheretal.sensitivity.tables.m2.v1<- function(df, levels, file)
{
	setkey(df, factor)
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	
	tmp		<- subset(df, stat=='P.rawbias.e0cp')
	set(tmp, NULL, 'v', tmp[, round(v, d=3)]*100)
	set(tmp, NULL, 'l95.bs', tmp[, round(l95.bs, d=3)]*100)
	set(tmp, NULL, 'u95.bs', tmp[, round(u95.bs, d=3)]*100)
	set(tmp, NULL, 'P.rawbias.e0cp', tmp[, paste('=\"',v,'\"&CHAR(13)&\"(',l95.bs,'COMMA ',u95.bs,')\"',sep='')])
	ans		<- subset(tmp, select=c(factor, t.period,P.rawbias.e0cp))
	tmp		<- subset(df, stat=='N.raw')
	set(tmp, NULL, 'v', tmp[, round(v, d=1)])
	set(tmp, NULL, 'l95.bs', tmp[, round(l95.bs, d=1)])
	set(tmp, NULL, 'u95.bs', tmp[, round(u95.bs, d=1)])
	set(tmp, NULL, 'N.raw', tmp[, paste('=\"',v,'\"&CHAR(13)&\"(',l95.bs,'COMMA ',u95.bs,')\"',sep='')])
	print(tmp[, list(v=sum(v), l95.bs=sum(l95.bs), u95.bs=sum(u95.bs)), by='t.period'])
	ans		<- merge(subset(tmp, select=c(factor, t.period,N.raw)), ans, by=c('factor','t.period'))
	tmp		<- subset(df, stat=='RI.rawbias.e0cp')
	set(tmp, NULL, 'v', tmp[, round(v, d=2)])
	set(tmp, NULL, 'l95.bs', tmp[, round(l95.bs, d=2)])
	set(tmp, NULL, 'u95.bs', tmp[, round(u95.bs, d=2)])
	set(tmp, NULL, 'RI.rawbias.e0cp', tmp[, paste('=\"',v,'\"&CHAR(13)&\"(',l95.bs,'COMMA ',u95.bs,')\"',sep='')])
	ans		<- merge(ans, subset(tmp, select=c(factor, t.period, RI.rawbias.e0cp)), by=c('factor','t.period'))
	tmp		<- subset(df, grepl('U',factor.ref) & stat=='RR.raw')
	set(tmp, NULL, 'v', tmp[, round(v, d=2)])
	set(tmp, NULL, 'l95.bs', tmp[, round(l95.bs, d=2)])
	set(tmp, NULL, 'u95.bs', tmp[, round(u95.bs, d=2)])
	set(tmp, NULL, 'RR.raw', tmp[, paste('=\"',v,'\"&CHAR(13)&\"(',l95.bs,'COMMA ',u95.bs,')\"',sep='')])
	ans		<- merge(ans, subset(tmp, select=c(factor, t.period, RR.raw)), by=c('factor','t.period'), all.x=TRUE)
	
	set(ans, NULL, 'factor', ans[, factor(factor, levels=levels)])
	setkey(ans, factor, t.period)
	write.csv(ans, file=file, eol="\r\n", row.names=FALSE)
	#in Excel replace COMMA with ,
	#right click on table or Apple-1 and set wrap text
}
######################################################################################
project.athena.Fisheretal.sensitivity.tables.m5<- function(df, levels, file)
{
	setkey(df, factor)
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	
	tmp		<- subset(df, stat=='P.rawbias.e0')
	set(tmp, NULL, 'v', tmp[, round(v, d=3)]*100)
	set(tmp, NULL, 'l95.bs', tmp[, round(l95.bs, d=3)]*100)
	set(tmp, NULL, 'u95.bs', tmp[, round(u95.bs, d=3)]*100)
	set(tmp, NULL, 'P.rawbias.e0', tmp[, paste('=\"',v,'\"&CHAR(13)&\"(',l95.bs,'COMMA ',u95.bs,')\"',sep='')])
	ans		<- subset(tmp, select=c(factor, t.period,P.rawbias.e0))
	tmp		<- subset(df, stat=='N.raw')
	set(tmp, NULL, 'v', tmp[, round(v, d=1)])
	set(tmp, NULL, 'l95.bs', tmp[, round(l95.bs, d=1)])
	set(tmp, NULL, 'u95.bs', tmp[, round(u95.bs, d=1)])
	set(tmp, NULL, 'N.raw', tmp[, paste('=\"',v,'\"&CHAR(13)&\"(',l95.bs,'COMMA ',u95.bs,')\"',sep='')])
	print(tmp[, list(v=sum(v), l95.bs=sum(l95.bs), u95.bs=sum(u95.bs)), by='t.period'])
	ans		<- merge(subset(tmp, select=c(factor, t.period,N.raw)), ans, by=c('factor','t.period'))
	tmp		<- subset(df, stat=='RI.rawbias.e0')
	set(tmp, NULL, 'v', tmp[, round(v, d=2)])
	set(tmp, NULL, 'l95.bs', tmp[, round(l95.bs, d=2)])
	set(tmp, NULL, 'u95.bs', tmp[, round(u95.bs, d=2)])
	set(tmp, NULL, 'RI.rawbias.e0', tmp[, paste('=\"',v,'\"&CHAR(13)&\"(',l95.bs,'COMMA ',u95.bs,')\"',sep='')])
	ans		<- merge(ans, subset(tmp, select=c(factor, t.period, RI.rawbias.e0)), by=c('factor','t.period'))
	tmp		<- subset(df, stat=='RR.raw')
	set(tmp, NULL, 'v', tmp[, round(v, d=2)])
	set(tmp, NULL, 'l95.bs', tmp[, round(l95.bs, d=2)])
	set(tmp, NULL, 'u95.bs', tmp[, round(u95.bs, d=2)])
	set(tmp, NULL, 'RR.raw', tmp[, paste('=\"',v,'\"&CHAR(13)&\"(',l95.bs,'COMMA ',u95.bs,')\"',sep='')])
	ans		<- merge(ans, subset(tmp, select=c(factor, t.period, RR.raw)), by=c('factor','t.period'), all.x=TRUE)
	
	set(ans, NULL, 'factor', ans[, factor(factor, levels=levels)])
	setkey(ans, factor, t.period)
	write.csv(ans, file=file, eol="\r\n", row.names=FALSE)
}

######################################################################################
project.athena.Fisheretal.sensitivity.getfigures<- function()
{	
	require(data.table)
	require(ape)
	require(grid)
	require(reshape2)
	require(ggplot2)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_140714",sep='/')
	outdir					<- paste(DATA,"fisheretal_140714",sep='/')
	indir					<- paste(DATA,"fisheretal_140722",sep='/')
	outdir					<- paste(DATA,"fisheretal_140722",sep='/')
	indir					<- paste(DATA,"fisheretal_140729",sep='/')
	outdir					<- paste(DATA,"fisheretal_140729",sep='/')
	indir					<- paste(DATA,"fisheretal_140814",sep='/')
	outdir					<- paste(DATA,"fisheretal_140814",sep='/')		
	indir					<- paste(DATA,"fisheretal_140817",sep='/')
	outdir					<- paste(DATA,"fisheretal_140817",sep='/')		
	indir					<- paste(DATA,"fisheretal_140828",sep='/')
	outdir					<- paste(DATA,"fisheretal_140828",sep='/')		
	indir					<- paste(DATA,"fisheretal_140905",sep='/')
	outdir					<- paste(DATA,"fisheretal_140905",sep='/')		
	indir					<- paste(DATA,"fisheretal_140929",sep='/')
	outdir					<- paste(DATA,"fisheretal_140929",sep='/')		
	indir					<- paste(DATA,"fisheretal_141022",sep='/')
	outdir					<- paste(DATA,"fisheretal_141022",sep='/')		
	indir					<- paste(DATA,"fisheretal_141108",sep='/')
	outdir					<- paste(DATA,"fisheretal_141108",sep='/')		
	
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))	
	#
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )	
	
	#	set up factor legends for m2Bwmx
	factor.color	<- c("#990000","#EF6548","#FDBB84","#0570B0","#74A9CF","#7A0177","#F768A1","#FCC5C0","#005824","#41AB5D","#ADDD8E")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
			'Undiagnosed,\n Chronic infection\n at diagnosis',
			'Undiagnosed,\n Unknown if recent',
			'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
			'Diagnosed,\n CD4 progression to >500',
			'Diagnosed,\n CD4 progression to [350-500]',
			'Diagnosed,\n CD4 progression to <350',			
			'Diagnosed,\n No CD4 measured',
			'ART initiated,\n No viral suppression',
			'ART initiated,\n No viral load measured',
			'ART initiated,\n Viral suppression')
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.suA.N","ART.vlNA","ART.suA.Y")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m2Bwmx')	
	#	m2Bwmx
	#	3la2H1V51
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3la2H1V100
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1V100'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3la2H1V200
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1V200'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3la2H1C3V51
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1C3V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3la2H1V100
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1C3V100'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3la2H1C3V200
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1C3V200'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#
	#	WTN
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1C3V51'
	method.RISK		<- 'm2BwmxMv.wtn.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1C3V100'
	method.RISK		<- 'm2BwmxMv.wtn.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1C3V200'
	method.RISK		<- 'm2BwmxMv.wtn.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)
	
	#	set up factor legends
	factor.color	<- c(	"#990000","#EF6548","#FDBB84",
			"#0C2C84","#0570B0","#74A9CF","#41B6C4","#35978F",  
			"#7A0177","#F768A1","#FCC5C0",
			"#1A9850","#A6D96A")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
			'Undiagnosed,\n Chronic infection\n at diagnosis',
			'Undiagnosed,\n Unknown if recent',
			'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
			'Diagnosed,\n CD4 progression to >500',
			'Diagnosed,\n CD4 progression to [350-500]',
			'Diagnosed,\n CD4 progression to <350',			
			'Diagnosed,\n No CD4 measured',
			'ART initiated,\n No viral suppression, >1 observations',
			'ART initiated,\n No viral suppression, 1 observation',
			'ART initiated,\n No viral load measured',			
			'ART initiated,\n Viral suppression, 1 observation',
			'ART initiated,\n Viral suppression, >1 observations')
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.suA.N2","ART.suA.N1","ART.vlNA","ART.suA.Y1","ART.suA.Y2")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m2Cwmx')
	#	updated stages
	#	set up factor legends
	factor.color	<- c(	"#990000","#EF6548","#FDBB84",
							"#0C2C84","#0570B0","#74A9CF","#41B6C4","#35978F",  
							"#FCC5C0","#F768A1","#7A0177",
							"#1A9850","#A6D96A")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
							'Undiagnosed,\n Chronic infection\n at diagnosis',
							'Undiagnosed,\n Unknown if recent',
							'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
							'Diagnosed,\n CD4 progression to >500',
							'Diagnosed,\n CD4 progression to [350-500]',
							'Diagnosed,\n CD4 progression to <350',			
							'Diagnosed,\n No CD4 measured',
							'ART initiated,\n Before first viral suppression',													
							'ART initiated,\n After first viral suppression\nNo viral load measured',		
							'ART initiated,\n After first viral suppression\nNo viral suppression',	
							'ART initiated,\n After first viral suppression\nViral suppression, 1 observation',
							'ART initiated,\n After first viral suppression\nViral suppression, >1 observations')
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.NotYetFirstSu","ART.vlNA","ART.suA.N","ART.suA.Y1","ART.suA.Y2")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m2Cwmx')	
	
	#
	#	WTN 3ma 	m2Cwmx
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ma2H1C3V100'
	method.RISK		<- 'm2CwmxMv.wtn.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Cwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ma1H1C3V100'
	method.RISK		<- 'm2CwmxMv.wtn.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Cwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#
	#	WTN 3na 	m2Cwmx
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3na2H1C3V100'
	method.RISK		<- 'm2CwmxMv.wtn.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Cwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3na1H1C3V100'
	method.RISK		<- 'm2CwmxMv.wtn.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Cwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#
	#	WTN 3oa 	m2Cwmx
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3oa2H1C3V100'
	method.RISK		<- 'm2CwmxMv.wtn.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Cwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3oa1H1C3V100'
	method.RISK		<- 'm2CwmxMv.wtn.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Cwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	
	#
	#	WTN m2Cwmx
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1C3V51'
	method.RISK		<- 'm2CwmxMv.wtn.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Cwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1C3V100'
	method.RISK		<- 'm2CwmxMv.wtn.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Cwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info, with.guide=FALSE)			
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1C3V200'
	method.RISK		<- 'm2CwmxMv.wtn.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Cwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	
	
	#	m2Bwmx
	#	3da
	method.DENOM	<- 'CLU'
	method.BRL		<- '3da'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0cp','P.bias.e0cp','P.raw','P.raw.e0cp','P.rawbias.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING, tmp, stat.select, outfile)
	#	m2Bwmx
	#	3ia			(zaga.all)
	method.DENOM	<- 'CLU'
	method.BRL		<- '3ia'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0cp','P.bias.e0cp','P.raw','P.raw.e0cp','P.rawbias.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile)
	#	m2Bwmx
	#	3ja			(corrected 3da)
	method.DENOM	<- 'CLU'
	method.BRL		<- '3ja'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0cp','P.bias.e0cp','P.raw','P.raw.e0cp','P.rawbias.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile)
	#	m2Bwmx
	#	3ka			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile)
	#	m2Bwmx
	#	3kaH			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3kaH'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)	
	#	m2Bwmx
	#	3kaH5			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3kaH5'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)	
	#	m2Bwmx
	#	3kaH2			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3kaH2'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)	
	#	m2Bwmx
	#	3kaH1			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3kaH1'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)		
	#	m2Bwmx
	#	3ka1.5H1			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka1.5H1'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3ka1.5H2			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka1.5H2'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3ka1.5H1C3V51			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka1.5H1C3V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)				
	#	m2Bwmx
	#	3ka2H1			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H1'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3ka2H1V51			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H1V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)				
	#	m2Bwmx
	#	3ka2H2			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H2'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3ka2H2V400			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H2V400'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3ka2H2V200			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H2V200'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3ka2H1C3V200			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H1C3V200'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3ka2H2V51			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H2V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3ka2H2C3V51			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H2C3V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING, tmp, stat.select, outfile, tperiod.info=tperiod.info)
	#	m2Bwmx
	#	3ka2H2C3V51			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H2C3V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'gmrf'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING, tmp, stat.select, outfile, tperiod.info=tperiod.info)					
	#	m2Bwmx
	#	3ka2H4C3V51			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H4C3V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3ka2H4C3V51			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H4C3V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'gmrf'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)				
	#	m2Bwmx
	#	3ka2H1C3V51			(ART.C)		**
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H1C3V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)					
	#	m2Bwmx
	#	3ka2H4C3V51			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H4C3V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''	
	method.DATING	<- 'gmrf'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)				
	#	m2Bwmx
	#	3ka2H5C3V51			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H5C3V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3ka2H0.5C3V51			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H0.5C3V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3ka2H0.3C3V51			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H0.3C3V51'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	
	#	m2Bwmx
	#	3ka2H2C3V51N0			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H2C3V51N0'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	#	m2Bwmx
	#	3ka2H4C3V51N0			(ART.C)
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3ka2H4C3V51N0'
	method.RISK		<- 'm2BwmxMv.tp'
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'
	tmp				<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	stat.select		<- c(	'P','P.e0','P.e0cp','P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  tmp, stat.select, outfile, tperiod.info=tperiod.info)			
	
	
	
	#RRterm
	run.tp			<- subset(runs.risk, method.denom=='CLU' & method.nodectime=='any' & method.brl=='3ka' & factor.ref=='U.4' & method.dating=='sasky' & 'm2BwmxMv.tp4'==method.risk  & stat=='RR.term' )
	run.tp			<- subset(runs.risk, method.denom=='CLU' & method.nodectime=='any' & method.brl=='3ia' & factor.ref=='U.4' & method.dating=='sasky' & 'm2BwmxMv.tp4'==method.risk  & stat=='RR.term' )
	run.tp			<- subset(runs.risk, method.denom=='CLU' & method.nodectime=='any' & method.brl=='3ja' & factor.ref=='U.4' & method.dating=='sasky' & 'm2BwmxMv.tp4'==method.risk  & stat=='RR.term' )
	#RRtermptx
	run.tp			<- subset(runs.risk, method.denom=='CLU' & method.nodectime=='any' & method.brl=='3ja' & factor.ref=='U.4' & method.dating=='sasky' & 'm2BwmxMv.tp4'==method.risk  & stat=='RR.term.ptx' )
	#
	subset(runs.risk, method.denom=='CLU' & method.nodectime=='any' & method.brl=='3ja' & method.dating=='sasky' & 'm2BwmxMv.tp4'==method.risk  & stat=='N.raw' )
	
	
	
	
	
	
	#tmp	<- subset(runs.risk, method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating=='sasky' & grepl('m2Bwmx.tp',method.risk) & grepl(method.clu,method.risk) & !grepl('now',method.risk) & (grepl('RR.',stat,fixed=1) | stat=='RR') )	
	#tmp	<- subset(runs.risk, method.denom==method.deno & method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m2BwmxMv.tp',method.risk) & grepl(method.clu,method.risk) & !grepl('now',method.risk) & stat=='RR.term.ptx' )
	#set(tmp, NULL, 'v', tmp[, round(v, d=3)])
	#set(tmp, NULL, 'l95.bs', tmp[, round(l95.bs, d=3)])
	#set(tmp, NULL, 'u95.bs', tmp[, round(u95.bs, d=3)])
	#subset(tmp, select=c(method.brl, factor, factor.ref, v, l95.bs, u95.bs))	
	#subset(tmp, grepl('ART.suA.Y',factor), select=c(method.brl, factor, factor.ref, v, l95.bs, u95.bs))
	
	#
	#	MODEL 2B time trends
	#
}
######################################################################################
project.athena.Fisheretal.sensitivity.getfigures.counts.m5<- function(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info)
{		
	df		<- subset(runs.table, method.denom==method.DENOM & method.BRL==method.brl & grepl(method.RISK,method.risk), c(stat, method.risk, factor, n, p))
	if(method.WEIGHT=='')
		df		<- subset(df, !grepl('wstar',method.risk) & !grepl('now',method.risk))
	if(method.WEIGHT!='')
		df		<- subset(df, grepl(method.WEIGHT,method.risk) )
	
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	set(df, NULL, 'factor', df[, factor(factor, levels=factors[, levels(factor)])])
	tmp		<- df[, 	list( 	n.coh=sum(n[stat=='X.msm'])/8, 
					p.seq= sum(n[stat=='X.seq']) / sum(n[stat=='X.msm']), 
					p.clu= sum(n[stat=='X.clu']) / sum(n[stat=='X.seq']), 
					p.pos= sum(n[stat=='YX']) / sum(n[stat=='X.clu']),
					p.pos.seq= sum(n[stat=='YX']) / sum(n[stat=='X.seq']),
					n.pos=sum(n[stat=='YX'])/8 ), by=c('factor','t.period')]	
	tmp		<- merge(tmp, df[, list( factor=factor[stat=='X.msm'], pp.coh=n[stat=='X.msm']/sum(n[stat=='X.msm']), pp.pos=n[stat=='YX']/sum(n[stat=='YX']) ), by='t.period'], by=c('factor','t.period'))
	ans		<- tmp
	set(ans, NULL, 'p.seq', ans[, round(p.seq, d=3)]*100)
	set(ans, NULL, 'p.clu', ans[, round(p.clu, d=3)]*100)
	set(ans, NULL, 'p.pos', ans[, round(p.pos, d=4)]*100)
	set(ans, NULL, 'p.pos.seq', ans[, round(p.pos.seq, d=5)]*100)
	set(ans, NULL, 'n.pos', ans[, round(n.pos, d=1)])
	set(ans, NULL, 'n.coh', ans[, round(n.coh, d=1)])	
	set(ans, NULL, 'pp.pos', ans[, round(pp.pos, d=3)]*100)
	set(ans, NULL, 'pp.coh', ans[, round(pp.coh, d=3)]*100)	
	#tperiod.long
	ans	<- merge(ans, tperiod.info, by='t.period')
	ans[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]		
	#factor.long	
	ans	<- merge(ans, factors, by='factor')		
	#group
	ans[, group:=NA_character_ ]
	set(ans, ans[,which(factor%in%c('t<=20','t<=25'))], 'group', '<30')
	set(ans, ans[,which(factor%in%c('t<=30','t<=35'))], 'group', '30-39')
	set(ans, ans[,which(factor%in%c('t<=40','t<=45','t<=100'))], 'group', '40-')
	set(ans, NULL, 'group', ans[, factor(group, levels=c('<30','30-39','40-'))])	
	#	plot numbers in cohort to illustrate censoring
	ggplot(ans, aes(x=t.period.long, y=n.coh, group=factor.legend, colour=factor.legend))  + labs(x='', y='potential transmission intervals in cohort\n(person-years)') +
			theme(legend.key.size=unit(11,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) +
			scale_colour_manual(name='from age group', values=ans[, unique(factor.color)]) +
			geom_line() + geom_point() + 
			facet_grid(. ~ group, scales='free', margins=FALSE)
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','censoring',".pdf", sep='')	
	ggsave(file=file, w=10, h=6)
	#	plot proportion in cohort to illustrate censoring
	ggplot(ans, aes(x=factor.legend, y=pp.coh, fill=factor.legend))  + labs(x='', y='potential transmission intervals in cohort\n(%)') +
			scale_y_continuous(breaks=seq(0,100,2)) +
			theme(legend.key.size=unit(11,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) +
			scale_fill_manual(name='from age group', values=ans[, unique(factor.color)]) +
			geom_bar(stat='identity',binwidth=1, position='dodge') + 
			facet_grid(. ~ t.period.long, scales='free', margins=FALSE)
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','censoringp',".pdf", sep='')	
	ggsave(file=file, w=18, h=5)		
	#	plot p.seq	
	ggplot(ans, aes(x=t.period.long, y=p.seq, fill=factor.legend)) +
			theme(legend.key.size=unit(11,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) +
			labs(x='', y=expression(frac('potential transmission intervals with a sequence','potential transmission intervals in cohort')*' * 100')) +
			scale_fill_manual(name='from age group', values=df[, unique(factor.color)]) + 
			geom_bar(stat='identity',binwidth=1, position='dodge', show_guide=FALSE) +
			theme(legend.key.size=unit(10.5,'mm'), plot.margin=unit(c(0,5,0,0),"mm")) + #coord_flip() +
			facet_grid(. ~ factor.legend, scales='free', margins=FALSE)
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','seqsampling',".pdf", sep='')	
	ggsave(file=file, w=18, h=5)
	#	plot p.clu	
	ggplot(ans, aes(x=t.period.long, y=p.clu, fill=factor.legend)) +
			theme(legend.key.size=unit(11,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) +
			labs(x='', y=expression(frac('potential transmission intervals in cluster','potential transmission intervals with a sequence' )*' * 100')) +
			scale_fill_manual(name='from age group', values=df[, unique(factor.color)]) + 
			geom_bar(stat='identity',binwidth=1, position='dodge', show_guide=FALSE) +
			theme(legend.key.size=unit(10.5,'mm'), plot.margin=unit(c(0,5,0,0),"mm")) + #coord_flip() +
			facet_grid(. ~ factor.legend, scales='free', margins=FALSE)
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','seqclustering',".pdf", sep='')	
	ggsave(file=file, w=18, h=5)
	#	plot p.pos.seq	
	ggplot(ans, aes(x=t.period.long, y=p.pos.seq, fill=factor.legend)) +
			theme(legend.key.size=unit(11,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) +
			labs(x='', y=expression(frac('potential transmission intervals with '*y[ijt]>0,'potential transmission intervals with a sequence')*' * 100')) +
			scale_fill_manual(name='from age group', values=df[, unique(factor.color)]) + 
			geom_bar(stat='identity',binwidth=1, position='dodge', show_guide=FALSE) +
			theme(legend.key.size=unit(10.5,'mm'), plot.margin=unit(c(0,5,0,0),"mm")) + #coord_flip() +
			facet_grid(. ~ factor.legend, scales='free', margins=FALSE)
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','seqYgreaterZero',".pdf", sep='')	
	ggsave(file=file, w=18, h=5)
	#	plot p.pos	
	ggplot(ans, aes(x=t.period.long, y=p.pos, fill=factor.legend)) +
			theme(legend.key.size=unit(11,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) +
			labs(x='', y=expression(frac('potential transmission intervals with '*y[ijt]>0,'potential transmission intervals in a cluster')*' * 100')) +
			scale_fill_manual(name='from age group', values=df[, unique(factor.color)]) + 
			geom_bar(stat='identity',binwidth=1, position='dodge', show_guide=FALSE) +
			theme(legend.key.size=unit(10.5,'mm'), plot.margin=unit(c(0,5,0,0),"mm")) + #coord_flip() +
			facet_grid(. ~ factor.legend, scales='free', margins=FALSE)
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','seqcluYgreaterZero',".pdf", sep='')	
	ggsave(file=file, w=18, h=5)
}	
######################################################################################
project.athena.Fisheretal.sensitivity.getfigures.comparetransmissionintervals<- function(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, group)
{			
	df		<- subset(runs.table, method.denom==method.DENOM & method.BRL==method.brl & grepl(method.RISK,method.risk), c(stat, method.risk, factor, n, p, l95.bs, u95.bs))
	if(method.WEIGHT=='')
		df		<- subset(df, !grepl('wstar',method.risk) & !grepl('now',method.risk))
	if(method.WEIGHT!='')
		df		<- subset(df, grepl(method.WEIGHT,method.risk) )
	set(df, NULL, 'factor', df[, as.character(factor)])	
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	set(df, NULL, 'factor', df[, factor(factor, levels=factors[, levels(factor)])])
		
	#	prop potential transmission intervals 
	tmp		<- subset(df, stat%in%c('X.msm.e0cp','nRec'))	
	tmp		<- dcast.data.table(tmp, t.period+factor+method.risk~stat, value.var='n')
	tmp[, n:= X.msm.e0cp/8/nRec*100]	
	tmp		<- merge(tmp, tmp[, list(factor=factor, p.coh=n/sum(n)), by='t.period'], by=c('t.period','factor'))
	#set(tmp, NULL, 'p', tmp[, round(p, d=3)]*100)
	set(tmp, NULL, 'n', tmp[, round(n, d=0)])			
	ans		<- subset(tmp, select=c(t.period, factor,  p.coh))
	#	prop likely transmission intervals 
	tmp		<- subset(df, stat%in%c('YX','nRecLkl'))
	tmp2	<- subset(df, stat=='YX', select=c(factor, t.period, p))
	tmp		<- dcast.data.table(tmp, t.period+factor+method.risk~stat, value.var='n')
	tmp		<- merge(tmp, tmp2, by=c('factor','t.period'))
	tmp[, n:= YX/nRecLkl*100]	
	#set(tmp, NULL, 'p', tmp[, round(p, d=3)]*100)
	set(tmp, NULL, 'n', tmp[, round(n, d=0)])	
	setnames(tmp, 'p', 'p.lkl')
	ans		<- merge(ans, subset(tmp, select=c(factor, t.period, p.lkl)), by=c('t.period','factor'))
	#	prop likely transmission intervals + expected missing intervals per 100 'recipient with lkl transmitter' in time period
	tmp		<- subset(df, stat%in%c('YX','Sx.e0cp','nRecLkl'))
	tmp2	<- subset(df, stat=='Sx.e0cp', select=c(factor, t.period, l95.bs, u95.bs))
	tmp		<- dcast.data.table(tmp, t.period+factor+method.risk~stat, value.var='n')
	tmp		<- merge(tmp, tmp2, by=c('factor','t.period'))
	tmp[, n:= (YX+Sx.e0cp/8)/nRecLkl*100]
	tmp[, l95.bs:= (YX+l95.bs/8)/nRecLkl*100]
	tmp[, u95.bs:= (YX+u95.bs/8)/nRecLkl*100]	
	tmp		<- merge(tmp, tmp[, list(factor=factor, p=n/sum(n)), by='t.period'], by=c('t.period','factor'))	
	#set(tmp, NULL, 'p', tmp[, round(p, d=3)]*100)
	setnames(tmp, 'p', 'p.lklm')
	ans		<- merge(ans, subset(tmp, select=c(factor, t.period, p.lklm)), by=c('t.period','factor'))
	#	ratio of proportions
	ans[, r.lkl:= p.lkl/p.coh]
	ans[, r.lklm:= p.lklm/p.coh]
	ans		<- melt(ans, measure.vars=c('r.lkl','r.lklm'))	
	# group
	if(group=='cascade')
	{
		scale.name	<- 'in treatment\ncascade stage'
		ans[, group:=ans[, substr(factor, 1, 1)]]
		set(ans, ans[,which(group=='A')], 'group', 'ART started')
		set(ans, ans[,which(group=='U')], 'group', 'Undiagnosed')
		set(ans, ans[,which(group=='D')], 'group', 'Diagnosed')
		set(ans, NULL, 'group', ans[, factor(group, levels=c('Undiagnosed','Diagnosed','ART started'))])		
	}
	if(group=='age')
	{		
		scale.name	<- 'from age group'
		ans[, group:=NA_character_ ]
		set(ans, ans[,which(factor%in%c('t<=20','t<=25'))], 'group', '<30')
		set(ans, ans[,which(factor%in%c('t<=30','t<=35'))], 'group', '30-39')
		set(ans, ans[,which(factor%in%c('t<=40','t<=45','t<=100'))], 'group', '40-')
		set(ans, NULL, 'group', ans[, factor(group, levels=c('<30','30-39','40-'))])
		theme.age <- function (base_size = 12, base_family = "") 
		{
			theme_grey(base_size = base_size, base_family = base_family) %+replace% 
					theme(	legend.key.size=unit(11,'mm'), 
							axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5), 
							panel.background = element_rect(fill = 'grey60'),
							panel.grid.minor.y = element_line(colour='grey70'),
							panel.grid.major = element_line(colour='grey70'),
							legend.key=element_rect(fill='grey60'))
		}
		theme_set(theme.age())
	}
	#tperiod.long
	ans	<- merge(ans, tperiod.info, by='t.period')
	ans[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]		
	#factor.long	
	ans	<- merge(ans, factors, by='factor')	
	ans	<- subset(ans, variable=='r.lklm')
	ggplot(ans, aes(x=t.period.long, y=value, group=factor.legend, colour=factor.legend, pch=variable))  + 				
			scale_colour_manual(name=scale.name, values=ans[, unique(factor.color)]) +
			scale_y_continuous(breaks=seq(0,10,1)) +			
			geom_line() + geom_point(shape=17) +
			labs(x='', y=expression(frac('% likely transmission intervals','% potential transmission intervals'))) +
			theme_bw() + theme(legend.key.size=unit(11,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) +
			facet_grid(. ~ group, scales='free_y', margins=FALSE)		
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','ratio',".pdf", sep='')	
	ggsave(file=file, w=9, h=7)
}
######################################################################################
project.athena.Fisheretal.sensitivity.getfigures.likelihood<- function()
{
	require(data.table)
	require(ape)
	require(grid)
	require(reshape2)
	require(ggplot2)
	#stop()
	#
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )	
	#
	#	CD4b
	#	
	indir			<- paste(DATA,"fisheretal",sep='/')
	outdir			<- paste(DATA,"fisheretal_140929",sep='/')
	infile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXm2_SEQ_3la2H1C3V100.R'
	outfile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXm2_SEQ_3la2H1C3V100_lkl.pdf'	
	file			<- paste(indir, '/', infile, sep='')
	load(file)
	YX				<- copy(YX.m2)
	set(YX, NULL, 'score.Y', YX[, score.Y*w.tn])
	YX				<- subset(YX, select=c(t.period, t, t.Patient, Patient, score.Y, CD4b, CD4c ))	
	setnames( YX, 'CD4b', 'factor')	
	#	set up factor legends
	factor.color	<- c("#990000","#EF6548","#FDBB84","#0570B0","#74A9CF","#7A0177","#F768A1","#FCC5C0","#005824","#41AB5D","#ADDD8E")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
							'Undiagnosed,\n Chronic infection\n at diagnosis',
							'Undiagnosed,\n Unknown if recent',
							'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
							'Diagnosed,\n CD4 progression to >500',
							'Diagnosed,\n CD4 progression to [350-500]',
							'Diagnosed,\n CD4 progression to <350',			
							'Diagnosed,\n No CD4 measured',
							'ART initiated,\n No viral suppression',
							'ART initiated,\n No viral load measured',
							'ART initiated,\n Viral suppression')
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.suA.N","ART.vlNA","ART.suA.Y")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='CD4b')	
	set(YX, NULL, 'factor', YX[, factor(as.character(factor), levels=tmp)])
	ans				<- merge(subset(factors, method.risk=='CD4b'), YX, by='factor')
	#
	scale.name	<- 'in treatment\ncascade stage'
	ans[, group:=ans[, substr(factor, 1, 1)]]
	set(ans, ans[,which(group=='A')], 'group', 'ART started')
	set(ans, ans[,which(group=='U')], 'group', 'Undiagnosed')
	set(ans, ans[,which(group=='D')], 'group', 'Diagnosed')
	set(ans, NULL, 'group', ans[, factor(group, levels=c('Undiagnosed','Diagnosed','ART started'))])
	#	tperiod.long
	ans				<- merge(ans, tperiod.info, by='t.period')
	ans[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]
	setkey(ans, factor)
	#	boxplot
	ggplot(ans, aes(x=factor.legend, y=score.Y, colour=factor.legend, fill=factor.legend), log='y')  + 
			labs(x='', y='likelihood of direct HIV transmission') +				
		scale_colour_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
		scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
		scale_y_continuous(breaks=c(0.1, 0.5, 1, 2, 5, 10, 20, 100), minor_breaks=c(0.1)) +
		coord_trans(ytrans="log10", limy=c(0.01, 1000)) +		
		geom_boxplot(outlier.shape = NA) +
		stat_summary(geom= "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +		
		theme_bw() + theme(legend.key.size=unit(11,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) 		
	file			<- paste(outdir, '/', outfile, sep='')
	ggsave(file=file, w=8, h=6)		
	
	#ggplot(YXf, aes(x=stage, y=brl, colour=stage), log='y') + geom_boxplot()
	#ggplot(YXf, aes(x=stage, y=telapsed, colour=stage), log='y') + geom_boxplot()
	#
	#	CD4c
	#
	infile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXm2_SEQ_3la2H1C3V100.R'
	outfile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXm2Cwmx.wtn_SEQ_3la2H1C3V100_lkl.pdf'
	file			<- paste(indir, '/', infile, sep='')
	load(file)
	YX				<- copy(YX.m2)
	set(YX, NULL, 'score.Y', YX[, score.Y*w.tn])
	YX				<- subset(YX, select=c(t.period, t, t.Patient, Patient, score.Y, CD4b, CD4c ))
	setnames( YX, 'CD4c', 'factor')	
	#	set up factor legends
	factor.color	<- c(	"#990000","#EF6548","#FDBB84",
							"#0C2C84","#0570B0","#74A9CF","#41B6C4","#35978F",  
							"#7A0177","#F768A1","#FCC5C0",
							"#1A9850","#A6D96A")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
			'Undiagnosed,\n Chronic infection\n at diagnosis',
			'Undiagnosed,\n Unknown if recent',
			'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
			'Diagnosed,\n CD4 progression to >500',
			'Diagnosed,\n CD4 progression to [350-500]',
			'Diagnosed,\n CD4 progression to <350',			
			'Diagnosed,\n No CD4 measured',
			'ART initiated,\n No viral suppression, >1 observations',
			'ART initiated,\n No viral suppression, 1 observation',
			'ART initiated,\n No viral load measured',			
			'ART initiated,\n Viral suppression, 1 observation',
			'ART initiated,\n Viral suppression, >1 observations')
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.suA.N2","ART.suA.N1","ART.vlNA","ART.suA.Y1","ART.suA.Y2")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='CD4c')
	
	set(YX, NULL, 'factor', YX[, factor(as.character(factor), levels=tmp)])
	ans				<- merge(subset(factors, method.risk=='CD4c'), YX, by='factor')
	#
	scale.name	<- 'in treatment\ncascade stage'
	ans[, group:=ans[, substr(factor, 1, 1)]]
	set(ans, ans[,which(group=='A')], 'group', 'ART started')
	set(ans, ans[,which(group=='U')], 'group', 'Undiagnosed')
	set(ans, ans[,which(group=='D')], 'group', 'Diagnosed')
	set(ans, NULL, 'group', ans[, factor(group, levels=c('Undiagnosed','Diagnosed','ART started'))])
	#	tperiod.long
	ans				<- merge(ans, tperiod.info, by='t.period')
	ans[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]
	setkey(ans, factor)
	#	boxplot
	ggplot(ans, aes(x=factor.legend, y=score.Y, colour=factor.legend, fill=factor.legend), log='y')  + 
			labs(x='', y='likelihood of direct HIV transmission') +				
			scale_colour_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
			scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)]) +
			scale_y_continuous(breaks=c(0.1, 0.5, 1, 2, 5, 10, 20, 100), minor_breaks=c(0.1)) +
			coord_trans(ytrans="log10", limy=c(0.01, 1000)) +		
			geom_boxplot(outlier.shape = NA) +
			stat_summary(geom= "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +		
			theme_bw() + theme(legend.key.size=unit(11,'mm'), legend.position = "bottom", legend.box = "vertical", legend.title.align = 0, axis.text.x=element_blank(), axis.text.x = element_blank(), axis.ticks.margin = unit(c(0,0,0,0), "lines")) +
			guides(fill=guide_legend(nrow=4))
	file			<- paste(outdir, '/', outfile, sep='')
	ggsave(file=file, w=8, h=6)		
	#
	#	CD4c updated
	#
	indir			<- paste(DATA,"fisheretal_141022",sep='/')
	outdir			<- paste(DATA,"fisheretal_141022",sep='/')
	infile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXm2_SEQ_3la2H1C3V100.R'
	outfile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXm2_SEQ_3la2H1C3V100_lkl.pdf'
	file			<- paste(indir, '/', infile, sep='')
	load(file)
	YX				<- copy(YX.m2)
	set(YX, NULL, 'score.Y', YX[, score.Y*w.tn])
	YX				<- subset(YX, select=c(t.period, t, t.Patient, Patient, score.Y, CD4b, CD4c ))
	setnames( YX, 'CD4c', 'factor')	
	factor.color	<- c(	"#990000","#EF6548","#FDBB84",
			"#0C2C84","#0570B0","#74A9CF","#41B6C4","#35978F",  
			"#FCC5C0","#F768A1","#7A0177",
			"#1A9850","#A6D96A")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
			'Undiagnosed,\n Chronic infection\n at diagnosis',
			'Undiagnosed,\n Unknown if recent',
			'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
			'Diagnosed,\n CD4 progression to >500',
			'Diagnosed,\n CD4 progression to [350-500]',
			'Diagnosed,\n CD4 progression to <350',			
			'Diagnosed,\n No CD4 measured',
			'ART initiated,\n Before first viral suppression',													
			'ART initiated,\n After first viral suppression\nNo viral load measured',		
			'ART initiated,\n After first viral suppression\nNo viral suppression',	
			'ART initiated,\n After first viral suppression\nViral suppression, 1 observation',
			'ART initiated,\n After first viral suppression\nViral suppression, >1 observations')
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.NotYetFirstSu","ART.vlNA","ART.suA.N","ART.suA.Y1","ART.suA.Y2")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='CD4c')	
	set(YX, NULL, 'factor', YX[, factor(as.character(factor), levels=tmp)])
	ans				<- merge(subset(factors, method.risk=='CD4c'), YX, by='factor')
	#
	scale.name	<- 'in treatment\ncascade stage'
	ans[, group:=ans[, substr(factor, 1, 1)]]
	set(ans, ans[,which(group=='A')], 'group', 'ART started')
	set(ans, ans[,which(group=='U')], 'group', 'Undiagnosed')
	set(ans, ans[,which(group=='D')], 'group', 'Diagnosed')
	set(ans, NULL, 'group', ans[, factor(group, levels=c('Undiagnosed','Diagnosed','ART started'))])
	#	tperiod.long
	ans				<- merge(ans, tperiod.info, by='t.period')
	ans[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]
	setkey(ans, factor)
	#	boxplot
	ggplot(ans, aes(x=factor.legend, y=score.Y, colour=factor.legend, fill=factor.legend), log='y')  + 
			labs(x='', y='likelihood of direct HIV transmission') +				
			scale_colour_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
			scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)]) +
			scale_y_continuous(breaks=c(0.1, 0.5, 1, 2, 5, 10, 20, 100), minor_breaks=c(0.1)) +
			coord_trans(ytrans="log10", limy=c(0.01, 1000)) +		
			geom_boxplot(outlier.shape = NA) +
			stat_summary(geom= "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +		
			theme_bw() + theme(legend.key.size=unit(11,'mm'), legend.position = "bottom", legend.box = "vertical", legend.title.align = 0, axis.text.x=element_blank(), axis.text.x = element_blank(), axis.ticks.margin = unit(c(0,0,0,0), "lines")) +
			guides(fill=guide_legend(nrow=4))
	file			<- paste(outdir, '/', outfile, sep='')
	ggsave(file=file, w=8, h=6)		
	
}
######################################################################################
project.athena.Fisheretal.sensitivity.getfigures.npotentialtransmissionintervals<- function(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, group)
{			
	df		<- subset(runs.table, method.denom==method.DENOM & method.BRL==method.brl & grepl(method.RISK,method.risk), c(stat, method.risk, factor, n, p))
	if(method.WEIGHT=='')
		df		<- subset(df, !grepl('wstar',method.risk) & !grepl('now',method.risk))
	if(method.WEIGHT!='')
		df		<- subset(df, grepl(method.WEIGHT,method.risk) )
	set(df, NULL, 'factor', df[, as.character(factor)])	
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	set(df, NULL, 'factor', df[, factor(factor, levels=factors[, levels(factor)])])
	
	tmp		<- subset(df, stat%in%c('X.seq','nRec'))
	tmp2	<- subset(df, stat=='X.seq', select=c(factor, t.period, p))
	tmp		<- dcast.data.table(tmp, t.period+factor+method.risk~stat, value.var='n')
	tmp		<- merge(tmp, tmp2, by=c('t.period','factor'))
	tmp[, n:= X.seq/nRec]	
	set(tmp, NULL, 'p', tmp[, round(p, d=3)]*100)
	set(tmp, NULL, 'n', tmp[, round(n, d=0)])	
	ans		<- tmp
	# group
	if(group=='cascade')
	{
		scale.name	<- 'in treatment\ncascade stage'
		ans[, group:=ans[, substr(factor, 1, 1)]]
		set(ans, ans[,which(group=='A')], 'group', 'ART started')
		set(ans, ans[,which(group=='U')], 'group', 'Undiagnosed')
		set(ans, ans[,which(group=='D')], 'group', 'Diagnosed')
		set(ans, NULL, 'group', ans[, factor(group, levels=c('Undiagnosed','Diagnosed','ART started'))])		
	}
	if(group=='age')
	{		
		scale.name	<- 'from age group'
		ans[, group:=NA_character_ ]
		set(ans, ans[,which(factor%in%c('t<=20','t<=25'))], 'group', '<30')
		set(ans, ans[,which(factor%in%c('t<=30','t<=35'))], 'group', '30-39')
		set(ans, ans[,which(factor%in%c('t<=40','t<=45','t<=100'))], 'group', '40-')
		set(ans, NULL, 'group', ans[, factor(group, levels=c('<30','30-39','40-'))])
		theme.age <- function (base_size = 12, base_family = "") 
		{
			theme_grey(base_size = base_size, base_family = base_family) %+replace% 
					theme(	legend.key.size=unit(11,'mm'), 
							axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5), 
							panel.background = element_rect(fill = 'grey60'),
							panel.grid.minor.y = element_line(colour='grey70'),
							panel.grid.major = element_line(colour='grey70'),
							legend.key=element_rect(fill='grey60'))
		}
		theme_set(theme.age())
	}
	#tperiod.long
	ans	<- merge(ans, tperiod.info, by='t.period')
	ans[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]		
	#factor.long	
	ans	<- merge(ans, factors, by='factor')	
	
	ggplot(ans, aes(x=t.period.long, y=n, group=factor.legend, colour=factor.legend))  + labs(x='', y='potential transmission intervals\nper recipient MSM\n(person-years)') +				
				scale_colour_manual(name=scale.name, values=ans[, unique(factor.color)]) +
				scale_y_continuous(breaks=seq(0,2000,200)) +
				geom_line() + geom_point() + 				
				theme_bw() + theme(legend.key.size=unit(11,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) +
				facet_grid(. ~ group, scales='free_y', margins=FALSE)		
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','censoring',".pdf", sep='')	
	ggsave(file=file, w=9, h=7)			
}
######################################################################################
project.athena.Fisheretal.sensitivity.getfigures.nlikelytransmissionintervals<- function(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, group)
{			
	df		<- subset(runs.table, method.denom==method.DENOM & method.BRL==method.brl & grepl(method.RISK,method.risk), c(stat, method.risk, factor, n, p, l95.bs, u95.bs))
	if(method.WEIGHT=='')
		df		<- subset(df, !grepl('wstar',method.risk) & !grepl('now',method.risk))
	if(method.WEIGHT!='')
		df		<- subset(df, grepl(method.WEIGHT,method.risk) )
	set(df, NULL, 'factor', df[, as.character(factor)])	
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	set(df, NULL, 'factor', df[, factor(factor, levels=factors[, levels(factor)])])
	
	tmp		<- subset(df, stat%in%c('YX','nRecLkl'))		
	tmp		<- dcast.data.table(tmp, t.period+factor+method.risk~stat, value.var='n')		
	tmp[, n:= YX/nRecLkl*100]
	ans		<- tmp
	ans[, stat:='observed']
	
	tmp		<- subset(df, stat%in%c('Sx.e0cp','nRecLkl'))	
	tmp2	<- subset(tmp, stat=='Sx.e0cp', c(t.period, factor, l95.bs,  u95.bs))
	tmp		<- dcast.data.table(tmp, t.period+factor+method.risk~stat, value.var='n')
	tmp		<- merge(tmp, tmp2, by=c('t.period','factor'))
	tmp[, Sx.e0cp:= Sx.e0cp/8/nRecLkl*100]
	tmp[, l95.bs:= l95.bs/8/nRecLkl*100]
	tmp[, u95.bs:= u95.bs/8/nRecLkl*100]
	setnames(tmp, 'Sx.e0cp', 'n')
	tmp[, stat:='expected missing']
	ans		<- rbind( subset(ans, select=c(t.period, factor, n, stat)), subset(tmp, select=c(t.period, factor, n, l95.bs, u95.bs, stat)), fill=TRUE)
	
	scale.name	<- 'in treatment\ncascade stage'
	ans[, group:=ans[, substr(factor, 1, 1)]]
	set(ans, ans[,which(group=='A')], 'group', 'ART started')
	set(ans, ans[,which(group=='U')], 'group', 'Undiagnosed')
	set(ans, ans[,which(group=='D')], 'group', 'Diagnosed')	
	
	tmp		<- subset(ans, stat=='expected missing', select=c(t.period, factor, stat, group, l95.bs, u95.bs))
	tmp		<- tmp[, list(stat=stat[1], factor=factor[1], l95.bs=sum(l95.bs), u95.bs=sum(u95.bs)), by=c('t.period','group')]	
	ans		<- merge( 	subset(ans, select=c(t.period, factor, group, stat, n)), 
						subset(tmp, select=c(t.period, factor, stat, l95.bs, u95.bs)), by=c('t.period','factor','stat'), all.x=1) 
	
	set(ans, NULL, 'group', ans[, paste(group,'\n(',stat,')',sep='')]) 
	set(ans, NULL, 'stat', NULL)		
	set(ans, NULL, 'group',	ans[, factor(group, levels=c('Undiagnosed\n(observed)','Undiagnosed\n(expected missing)','Diagnosed\n(observed)','Diagnosed\n(expected missing)','ART started\n(observed)','ART started\n(expected missing)'))])				
				
	#tperiod.long
	ans	<- merge(ans, tperiod.info, by='t.period')
	ans[, t.period.long:= paste(t.period.min, '-', t.period.max,sep='')]
	set(ans, NULL, 't.period.long', ans[,factor(t.period.long, levels= tperiod.info[, paste(t.period.min, '-', t.period.max,sep='')])])	
	#factor.long	
	ans	<- merge(ans, factors, by='factor')	
	
	ggplot(ans, aes(x=t.period.long, y=n, ymin=l95.bs, ymax=u95.bs, group=group, fill=factor.legend))  + 
			labs(x='', y='probable transmission intervals\nper 100 recipient MSM (person-years)') +				
			scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
			scale_y_continuous(breaks=seq(0,1000,50), minor_breaks=seq(0,1000,10)) +
			geom_bar(stat='identity') + geom_errorbar(stat='identity',width=.2) +  				
			theme_bw() + theme(legend.key.size=unit(12,'mm'),panel.grid.minor=element_line(colour="grey70", size=0.2), panel.grid.major=element_line(colour="grey70")) +
			facet_grid(. ~ group, scales='free_y', margins=FALSE)	
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','lkladj',".pdf", sep='')	
	ggsave(file=file, w=8/4*6, h=5)	
}
######################################################################################
project.athena.Fisheretal.sensitivity.getfigures.npotentialtransmissionintervals.adjusted<- function(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, group)
{			
	df		<- subset(runs.table, method.denom==method.DENOM & method.BRL==method.brl & grepl(method.RISK,method.risk), c(stat, method.risk, factor, n, p, l95.bs, u95.bs))
	if(method.WEIGHT=='')
		df		<- subset(df, !grepl('wstar',method.risk) & !grepl('now',method.risk))
	if(method.WEIGHT!='')
		df		<- subset(df, grepl(method.WEIGHT,method.risk) )
	set(df, NULL, 'factor', df[, as.character(factor)])	
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	set(df, NULL, 'factor', df[, factor(factor, levels=factors[, levels(factor)])])
	
	tmp		<- subset(df, stat%in%c('X.msm','X.msm.e0cp','nRec'))	
	tmp2	<- subset(df, stat=='X.msm.e0cp', select=c(factor, t.period, l95.bs, u95.bs))
	tmp		<- dcast.data.table(tmp, t.period+factor+method.risk~stat, value.var='n')
	tmp		<- merge(tmp, tmp2, by=c('t.period','factor'))	
	tmp[, n:= round(X.msm/nRec)]
	tmp[, X.msm.e0cp:= round(X.msm.e0cp/8/nRec)]
	tmp[, l95.bs:= round(l95.bs/8/nRec)]
	tmp[, u95.bs:= round(u95.bs/8/nRec)]
	ans		<- tmp
	
	scale.name	<- 'in treatment\ncascade stage'
	ans[, group:=ans[, substr(factor, 1, 1)]]
	set(ans, ans[,which(group=='A')], 'group', 'ART started')
	set(ans, ans[,which(group=='U')], 'group', 'Undiagnosed\n(observed)')
	set(ans, ans[,which(group=='D')], 'group', 'Diagnosed')
		
	
	tmp		<- subset(ans, group=='Undiagnosed\n(observed)', select=c(t.period, factor, X.msm.e0cp, l95.bs, u95.bs))
	setnames(tmp, 'X.msm.e0cp', 'n')
	tmp[, group:='Undiagnosed\n(adjusted for\nright censoring)']
	tmp2	<- tmp[, list(factor='U', l95.bs=sum(l95.bs), u95.bs=sum(u95.bs)), by=c('t.period','group')]
	tmp		<- merge(subset(tmp, select=c(t.period, factor, n, group)), subset(tmp2, select=c(t.period, factor, l95.bs, u95.bs)), by=c('t.period','factor'), all.x=1) 
	
	ans	<- subset(ans, select=c(t.period, factor, n, l95.bs, u95.bs, group))
	set(ans, NULL, c('l95.bs', 'u95.bs'), NA_real_)
	ans	<- rbind(ans, tmp, use.names=TRUE)
	set(ans, NULL, 'group', ans[, factor(group, levels=c('Undiagnosed\n(observed)','Undiagnosed\n(adjusted for\nright censoring)','Diagnosed','ART started'))])
	
	#tperiod.long
	ans	<- merge(ans, tperiod.info, by='t.period')	
	ans[, t.period.long:= paste(t.period.min, '-', t.period.max,sep='')]
	set(ans, NULL, 't.period.long', ans[,factor(t.period.long, levels= tperiod.info[, paste(t.period.min, '-', t.period.max,sep='')])])
	#factor.long	
	ans	<- merge(ans, factors, by='factor')	
	
	ggplot(ans, aes(x=t.period.long, y=n, ymin=l95.bs, ymax=u95.bs, group=group, fill=factor.legend))  + 
			labs(x='', y='potential transmission intervals\nper recipient MSM\n(person-years)') +				
			scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)]) +
			scale_y_continuous(breaks=seq(0,20000,500), minor_breaks=seq(0,20000,100)) +
			geom_bar(stat='identity') + geom_errorbar(stat='identity',width=.2) +  				
			theme_bw() + theme(legend.key.size=unit(12,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) +
			facet_grid(. ~ group, scales='free_y', margins=FALSE)	
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','ptadjwl',".pdf", sep='')	
	ggsave(file=file, w=9, h=8)
	if(1)
	{
		ggplot(ans, aes(x=t.period.long, y=n, ymin=l95.bs, ymax=u95.bs, group=group, fill=factor.legend))  + 
				labs(x='', y='potential transmission intervals\nper recipient MSM\n(person-years)') +				
				scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
				scale_y_continuous(breaks=seq(0,20000,1000), minor_breaks=seq(0,20000,200)) +
				geom_bar(stat='identity') + geom_errorbar(stat='identity',width=.2) +  				
				theme_bw() + theme(legend.key.size=unit(12,'mm'), axis.title.x = element_text(vjust=-5), panel.grid.minor=element_line(colour="grey70", size=0.2), panel.grid.major=element_line(colour="grey70")) +				
				facet_grid(. ~ group, scales='free_y', margins=FALSE)	
		file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','ptadj',".pdf", sep='')	
		ggsave(file=file, w=9, h=5)
	}
	if(0)
	{
		ggplot(ans, aes(x=t.period, y=n, ymin=l95.bs, ymax=u95.bs, group=group, fill=factor.legend, colour=t.period))  + 
				labs(x='', y='potential transmission intervals\nper recipient MSM (person-years)') +				
				scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
				scale_color_manual(name='time of diagnosis\nof recipient MSM', values=rep('transparent',nrow(tperiod.info)), label= paste('period ',tperiod.info[, t.period],':\n',tperiod.info[, t.period.min],' - ',tperiod.info[, t.period.max],sep='')) +
				scale_x_discrete(labels=paste('period',tperiod.info[, t.period])) +
				scale_y_continuous(breaks=seq(0,20000,1000), minor_breaks=seq(0,20000,200)) +
				geom_bar(stat='identity') + geom_errorbar(stat='identity',width=.2,colour='black') +  				
				theme_bw() + theme(legend.key.height=unit(12,'mm'), legend.key.width=unit(0,'mm'), legend.key=element_rect(colour='transparent', fill='transparent'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) +
				facet_grid(. ~ group, scales='free_y', margins=FALSE)	
		file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','ptadj',".pdf", sep='')	
		ggsave(file=file, w=9, h=5)	
	}
		
	
}
######################################################################################
project.athena.Fisheretal.sensitivity.getfigures.pseq<- function(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, group)
{	
	require(Hmisc)
	
	df		<- subset(runs.table, method.denom==method.DENOM & method.BRL==method.brl & grepl(method.RISK,method.risk), c(stat, method.risk, factor, n, p))
	if(method.WEIGHT=='')
		df		<- subset(df, !grepl('wstar',method.risk) & !grepl('now',method.risk))
	if(method.WEIGHT!='')
		df		<- subset(df, grepl(method.WEIGHT,method.risk) )
	set(df, NULL, 'factor', df[, as.character(factor)])	
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	set(df, NULL, 'factor', df[, factor(factor, levels=factors[, levels(factor)])])
	tmp		<- df[, 	{
				#BINCONF not appropriate because potential transmission intervals are not independent
				p.seq.cis	<- as.double(binconf( sum(n[stat=='X.seq']), sum(n[stat=='X.msm']), alpha=0.05, method= "wilson", include.x=FALSE, include.n=FALSE, return.df=FALSE))
				list( 	p.seq= sum(n[stat=='X.seq']) / sum(n[stat=='X.msm']),
						p.seq.l= p.seq.cis[2],
						p.seq.u= p.seq.cis[3] )
			}, by=c('factor','t.period')]
	tmp		<- merge(tmp, df[, list( factor=factor[stat=='X.msm'], pp.coh=n[stat=='X.msm']/sum(n[stat=='X.msm']), pp.pos=n[stat=='YX']/sum(n[stat=='YX']) ), by='t.period'], by=c('factor','t.period'))
	ans		<- tmp
	set(ans, NULL, 'p.seq', ans[, round(p.seq, d=3)]*100)		
	# group
	if(group=='cascade')
	{
		scale.name	<- 'from cascade stage'
		ans[, group:=ans[, substr(factor, 1, 1)]]
		set(ans, ans[,which(group=='A')], 'group', 'cART initiated')
		set(ans, ans[,which(group=='U')], 'group', 'Undiagnosed')
		set(ans, ans[,which(group=='D')], 'group', 'Diagnosed')
		set(ans, NULL, 'group', ans[, factor(group, levels=c('Undiagnosed','Diagnosed','cART initiated'))])
		theme.cascade <- function (base_size = 12, base_family = "") 
		{
			theme_grey(base_size = base_size, base_family = base_family) %+replace% 
					theme(	legend.key.size=unit(11,'mm'), 
							axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5))
		}		
		theme_set(theme.cascade())	
	}
	if(group=='age')
	{		
		scale.name	<- 'from age group'
		ans[, group:=NA_character_ ]
		set(ans, ans[,which(factor%in%c('t<=20','t<=25'))], 'group', '<30')
		set(ans, ans[,which(factor%in%c('t<=30','t<=35'))], 'group', '30-39')
		set(ans, ans[,which(factor%in%c('t<=40','t<=45','t<=100'))], 'group', '40-')
		set(ans, NULL, 'group', ans[, factor(group, levels=c('<30','30-39','40-'))])
		theme.age <- function (base_size = 12, base_family = "") 
		{
			theme_grey(base_size = base_size, base_family = base_family) %+replace% 
					theme(	legend.key.size=unit(11,'mm'), 
							axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5), 
							panel.background = element_rect(fill = 'grey60'),
							panel.grid.minor.y = element_line(colour='grey70'),
							panel.grid.major = element_line(colour='grey70'),
							legend.key=element_rect(fill='grey60'))
		}
		theme_set(theme.age())
	}
	#tperiod.long
	ans	<- merge(ans, tperiod.info, by='t.period')
	ans[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]		
	#factor.long	
	ans	<- merge(ans, factors, by='factor')	
	#	plot p.seq	
	if(0)
	{
		ggplot(ans, aes(x=factor.legend, y=p.seq, fill=factor.legend)) +			
				labs(x='', y=expression(frac(atop('potential transmission intervals','with a sequence'),atop('potential transmission intervals','in cohort'))*' * 100')) +			
				scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)]) + 
				geom_bar(stat='identity',binwidth=1, position='dodge', show_guide=FALSE) +
				theme(plot.margin=unit(c(0,5,0,0),"mm")) + #coord_flip() +
				facet_grid(. ~ t.period.long, scales='free', margins=FALSE)
		file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','seqsampling',".pdf", sep='')	
		ggsave(file=file, w=25, h=8)		
	}
	if(1)
	{
		ggplot(ans, aes(x=factor.legend, y=p.seq, pch=t.period.long)) + geom_point() +
				theme_bw() + theme(legend.position='bottom', axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) + 
				scale_y_continuous(breaks=seq(0,100,10)) +
				labs(x='', y='proportion of\npotentialtransmission intervals\nwith a sequence\n(%)',pch='observation\nperiod')
		file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','seqsampling',".pdf", sep='')
		cat(paste('plot to',file))
		ggsave(file=file, w=9, h=6)
	}	
}
######################################################################################
project.athena.Fisheretal.sensitivity.getfigures.counts<- function(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, group)
{	
	require(Hmisc)
	
	df		<- subset(runs.table, method.denom==method.DENOM & method.BRL==method.brl & grepl(method.RISK,method.risk), c(stat, method.risk, factor, n, p))
	if(method.WEIGHT=='')
		df		<- subset(df, !grepl('wstar',method.risk) & !grepl('now',method.risk))
	if(method.WEIGHT!='')
		df		<- subset(df, grepl(method.WEIGHT,method.risk) )
	set(df, NULL, 'factor', df[, as.character(factor)])	
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	set(df, NULL, 'factor', df[, factor(factor, levels=factors[, levels(factor)])])
	tmp		<- df[, 	{
							#BINCONF not appropriate because potential transmission intervals are not independent
							p.seq.cis	<- as.double(binconf( sum(n[stat=='X.seq']), sum(n[stat=='X.msm']), alpha=0.05, method= "wilson", include.x=FALSE, include.n=FALSE, return.df=FALSE))
							list( 	n.coh=sum(n[stat=='X.msm'])/8, 
									p.seq= sum(n[stat=='X.seq']) / sum(n[stat=='X.msm']),
									p.seq.l= p.seq.cis[2],
									p.seq.u= p.seq.cis[3],
									p.clu= sum(n[stat=='X.clu']) / sum(n[stat=='X.seq']), 
									p.pos= sum(n[stat=='YX']) / sum(n[stat=='X.clu']), 
									p.pos.seq= sum(n[stat=='YX']) / sum(n[stat=='X.seq']),
									n.pos=sum(n[stat=='YX'])/8 )
						}, by=c('factor','t.period')]
	tmp		<- merge(tmp, df[, list( factor=factor[stat=='X.msm'], pp.coh=n[stat=='X.msm']/sum(n[stat=='X.msm']), pp.pos=n[stat=='YX']/sum(n[stat=='YX']) ), by='t.period'], by=c('factor','t.period'))
	ans		<- tmp
	set(ans, NULL, 'p.seq', ans[, round(p.seq, d=3)]*100)
	set(ans, NULL, 'p.clu', ans[, round(p.clu, d=3)]*100)
	set(ans, NULL, 'p.pos', ans[, round(p.pos, d=4)]*100)
	set(ans, NULL, 'p.pos.seq', ans[, round(p.pos.seq, d=5)]*100)
	set(ans, NULL, 'n.pos', ans[, round(n.pos, d=1)])
	set(ans, NULL, 'n.coh', ans[, round(n.coh, d=1)])	
	set(ans, NULL, 'pp.pos', ans[, round(pp.pos, d=3)]*100)
	set(ans, NULL, 'pp.coh', ans[, round(pp.coh, d=3)]*100)	
	# group
	if(group=='cascade')
	{
		scale.name	<- 'from cascade stage'
		ans[, group:=ans[, substr(factor, 1, 1)]]
		set(ans, ans[,which(group=='A')], 'group', 'cART initiated')
		set(ans, ans[,which(group=='U')], 'group', 'Undiagnosed')
		set(ans, ans[,which(group=='D')], 'group', 'Diagnosed')
		set(ans, NULL, 'group', ans[, factor(group, levels=c('Undiagnosed','Diagnosed','cART initiated'))])
		theme.cascade <- function (base_size = 12, base_family = "") 
			{
				theme_grey(base_size = base_size, base_family = base_family) %+replace% 
				theme(	legend.key.size=unit(11,'mm'), 
						axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5))
			}		
		theme_set(theme.cascade())	
	}
	if(group=='age')
	{		
		scale.name	<- 'from age group'
		ans[, group:=NA_character_ ]
		set(ans, ans[,which(factor%in%c('t<=20','t<=25'))], 'group', '<30')
		set(ans, ans[,which(factor%in%c('t<=30','t<=35'))], 'group', '30-39')
		set(ans, ans[,which(factor%in%c('t<=40','t<=45','t<=100'))], 'group', '40-')
		set(ans, NULL, 'group', ans[, factor(group, levels=c('<30','30-39','40-'))])
		theme.age <- function (base_size = 12, base_family = "") 
		{
			theme_grey(base_size = base_size, base_family = base_family) %+replace% 
					theme(	legend.key.size=unit(11,'mm'), 
							axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5), 
							panel.background = element_rect(fill = 'grey60'),
							panel.grid.minor.y = element_line(colour='grey70'),
							panel.grid.major = element_line(colour='grey70'),
							legend.key=element_rect(fill='grey60'))
		}
		theme_set(theme.age())
	}
	#tperiod.long
	ans	<- merge(ans, tperiod.info, by='t.period')
	ans[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]		
	#factor.long	
	ans	<- merge(ans, factors, by='factor')	
	#	plot numbers in cohort to illustrate censoring
	ggplot(ans, aes(x=t.period.long, y=n.coh, group=factor.legend, colour=factor.legend))  + labs(x='', y='potential transmission intervals in cohort\n(person-years)') +				
				scale_colour_manual(name=scale.name, values=ans[, unique(factor.color)]) +
				geom_line() + geom_point() + 
				facet_grid(. ~ group, scales='free', margins=FALSE)		
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','censoring',".pdf", sep='')	
	ggsave(file=file, w=10, h=6)
	#	plot proportion in cohort to illustrate censoring
	ggplot(ans, aes(x=factor.legend, y=pp.coh, fill=factor.legend))  + labs(x='', y='potential transmission intervals in cohort\n(%)') +
			scale_y_continuous(breaks=seq(0,100,2)) +			
			scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)]) +
			geom_bar(stat='identity',binwidth=1, position='dodge', show_guide=FALSE) + 
			facet_grid(. ~ t.period.long, scales='free', margins=FALSE)
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','censoringp',".pdf", sep='')	
	ggsave(file=file, w=25, h=8)			
	#	plot p.clu	
	ggplot(ans, aes(x=factor.legend, y=p.clu, fill=factor.legend)) +			
			labs(x='', y=expression(frac(atop('potential transmission intervals','in a cluster'),atop('potential transmission intervals','with a sequence'))*' * 100')) +			
			scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)]) + 
			geom_bar(stat='identity',binwidth=1, position='dodge', show_guide=FALSE) +
			theme(plot.margin=unit(c(0,5,0,0),"mm")) + #coord_flip() +
			facet_grid(. ~ t.period.long, scales='free', margins=FALSE)
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','seqclustering',".pdf", sep='')	
	ggsave(file=file, w=25, h=8)
	#	plot p.pos.seq	
	ggplot(ans, aes(x=factor.legend, y=p.pos.seq, fill=factor.legend)) +
			labs(x='', y=expression(frac(atop('potential transmission intervals','with evidence for direct HIV transmission'),atop('potential transmission intervals','with a sequence'))*' * 100')) +
			scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)]) + 
			geom_bar(stat='identity',binwidth=1, position='dodge', show_guide=FALSE) +
			theme(plot.margin=unit(c(0,5,0,0),"mm")) + #coord_flip() +
			facet_grid(. ~ t.period.long, scales='free', margins=FALSE)
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','seqYgreaterZero',".pdf", sep='')	
	ggsave(file=file, w=25, h=8)
	#	plot p.pos	
	ggplot(ans, aes(x=factor.legend, y=p.pos, fill=factor.legend)) +
			labs(x='', y=expression(frac(atop('potential transmission intervals','with evidence for direct HIV transmission'),atop('potential transmission intervals','in a cluster'))*' * 100')) +
			scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)]) + 
			geom_bar(stat='identity',binwidth=1, position='dodge', show_guide=FALSE) +
			theme(plot.margin=unit(c(0,5,0,0),"mm")) + #coord_flip() +
			facet_grid(. ~ t.period.long, scales='free', margins=FALSE)
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','seqcluYgreaterZero',".pdf", sep='')	
	ggsave(file=file, w=25, h=8)
}
######################################################################################
project.athena.Fisheretal.sensitivity.getfigures.m2<- function(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING, factors, stat.select, outfile, tperiod.info=NULL, with.guide=FALSE)
{
	run.tp			<- subset(runs.risk, method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating==method.DATING & grepl(method.RISK,method.risk)  & (grepl('P.',stat,fixed=1) | stat=='P') )
	if(method.WEIGHT=='')
		run.tp		<- subset(run.tp, !grepl('wstar',method.risk) & !grepl('now',method.risk))
	if(method.WEIGHT!='')
		run.tp		<- subset(run.tp, grepl(method.WEIGHT,method.risk) )
	if(is.null(tperiod.info))
	{
		tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
		set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
		set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )		
	}		
	run.tp[, m50.bs:=NULL]
	ylab		<- "Proportion of transmissions"	
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])
	run.tp[, cascade.stage:=run.tp[, substr(factor, 1, 1)]]	
	set(run.tp, run.tp[,which(cascade.stage=='A')], 'cascade.stage', 'cART initiated')
	set(run.tp, run.tp[,which(cascade.stage=='U')], 'cascade.stage', 'Undiagnosed')
	set(run.tp, run.tp[,which(cascade.stage=='D' & factor%in%c("Dtl500","Dtl350"))], 'cascade.stage', 'Diagnosed\n CD4<=500')
	set(run.tp, run.tp[,which(cascade.stage=='D')], 'cascade.stage', 'Diagnosed')	
	set(run.tp, NULL, 'cascade.stage', run.tp[, factor(cascade.stage, levels=c('Undiagnosed','Diagnosed','Diagnosed\n CD4<=500','cART initiated'))])
	set(run.tp, NULL, 'factor', run.tp[, factor(factor, levels=factors[, levels(factor)])])
	run.tp	<- merge(run.tp, factors, by='factor')	
	run.tp	<- merge(run.tp, tperiod.info, by='t.period')
	run.tp[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]		
	dummy	<- lapply(seq_along(stat.select), function(i)
			{
				cat(paste('\nprocess', stat.select[i]))
				ggplot(subset(run.tp, !is.na(v) & stat==stat.select[i]), aes(x=factor, y=v, fill=factor.legend, colour=factor.legend)) + labs(x="", y=ylab) + 
						scale_y_continuous(breaks=seq(0,0.9,0.1), minor_breaks=seq(0,0.9,0.01), labels=paste(seq(0,0.9,0.1)*100,'%',sep='')) + 
						scale_x_discrete(breaks=NULL, limits=run.tp[, intersect(levels(factor), as.character(factor))]) +
						scale_fill_manual(name='from cascade stage', values=run.tp[, unique(factor.color)], guide=with.guide) + 
						scale_colour_manual(name='from cascade stage', values = rep('black',run.tp[, length(unique(factor))]), guide=with.guide) +
						#scale_fill_brewer(palette='PRGn',name='from cascade stage') + scale_colour_manual(name='from cascade stage', values = rep('black',11)) +					
						#guides(colour=FALSE, fill = guide_legend(override.aes = list(size=5))) +
						theme_bw() + theme(legend.key.size=unit(10.5,'mm'), plot.margin=unit(c(0,0,-5,0),"mm")) + #coord_flip() +
						geom_bar(stat='identity',binwidth=1, position='dodge') + geom_errorbar(aes(ymin=l95.bs, ymax=u95.bs), width=0.3, position=position_dodge(width=0.9))	+ 
						facet_grid(. ~ t.period.long, margins=FALSE)
				file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_',run.tp[1, method.recentctime],'_',run.tp[1, method.denom], '_',run.tp[1, method.dating], '_',run.tp[1, method.brl],'_',run.tp[1, method.risk],'_',stat.select[i],".pdf", sep='')
				cat(paste('\nsave to file',file))							
				ggsave(file=file, w=12,h=4)							
			})
	#
	ylab		<- "Relative transmissibility"
	run.tp		<- subset(runs.risk, method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating==method.DATING & grepl(method.RISK,method.risk)  & (grepl('RI.',stat,fixed=1) | stat=='RI') )
	stat.select	<- gsub('P','RI', stat.select)
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])
	run.tp[, group:=run.tp[, substr(factor, 1, 1)]]	
	set(run.tp, run.tp[,which(group=='A' & grepl('suA\\.Y',factor))], 'group', 'ART started and\nviral load < 100 / ml')
	set(run.tp, run.tp[,which(group=='A' & !grepl('suA\\.Y',factor))], 'group', 'ART started and\nviral load > 100 / ml\nor not measured')
	set(run.tp, run.tp[,which(group=='A')], 'group', 'ART initiated')
	set(run.tp, run.tp[,which(group=='U')], 'group', 'Undiagnosed')	
	set(run.tp, run.tp[,which(group=='D')], 'group', 'Diagnosed')	
	set(run.tp, NULL, 'group', run.tp[, factor(group, levels=c('Undiagnosed','Diagnosed','ART started and\nviral load > 100 / ml\nor not measured','ART started and\nviral load < 100 / ml'))])	
	set(run.tp, NULL, 'factor', run.tp[, factor(factor, levels=factors[, levels(factor)])])
	run.tp	<- merge(run.tp, factors, by='factor')		
	run.tp	<- merge(run.tp, tperiod.info, by='t.period')
	run.tp[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]
	color	<- run.tp[, unique(factor.color)]
	dummy	<- lapply(seq_along(stat.select), function(i)
			{
				cat(paste('\nprocess', stat.select[i]))
				file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_',run.tp[1, method.recentctime],'_',run.tp[1, method.denom], '_',run.tp[1, method.dating], '_',run.tp[1, method.brl],'_',run.tp[1, method.risk],'_',stat.select[i],".pdf", sep='')
				cat(paste('\nsave to file',file))
				if(0)
				{
					ggplot(subset(run.tp, !is.na(v) & stat==stat.select[i]), aes(x=t.period.long, y=v, fill=factor.legend, colour=factor.legend)) + 
							labs(x="", y=ylab) + 						
							#scale_y_continuous(breaks=seq(0,0.3,0.1), labels=paste(seq(0,0.3,0.1)*100,'%',sep='')) + scale_x_discrete(breaks=NULL, limits=rev(c("<20","20-24","25-29","30-34","35-39","40-44","45-"))) +
							scale_fill_manual(name='from cascade stage', values = color, guide=with.guide) + 
							scale_colour_manual(name='from cascade stage', values = rep('black',run.tp[, length(unique(factor))]), guide=with.guide) +
							#guides(colour=FALSE, fill = guide_legend(override.aes = list(size=5))) +
							theme_bw() + theme(legend.key.size=unit(10.5,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) + #coord_flip() +
							geom_bar(stat='identity',binwidth=1, position='dodge', show_guide=FALSE)	+ geom_errorbar(aes(ymin=l95.bs, ymax=u95.bs), width=0.3, position=position_dodge(width=0.9))	+ 
							facet_grid(. ~ factor.legend, margins=FALSE)
					ggsave(file=file, w=16,h=5.5)
				}
				if(1)
				{
					ggplot(subset(run.tp, !is.na(v) & stat==stat.select[i]), aes(x=t.period.long, y=v, ymin=l95.bs, ymax=u95.bs, fill=factor.legend, colour=factor.legend, group=factor.legend)) + 
							labs(x="", y=ylab) + #coord_trans(ytrans=my_trans, limy=c(0.01, 15)) +	
							#scale_y_continuous(breaks=c(seq(0,1,0.2), seq(1,4,1), seq(4,20,2))) +			
							scale_fill_manual(name='from cascade stage', values = color, guide=with.guide) + 
							scale_colour_manual(name='from cascade stage', values = color, guide=with.guide) +
							geom_hline(yintercept=1) + geom_errorbar(position=position_dodge(.5), width=0.4) + geom_point(position=position_dodge(.5) , shape=18, size=3) +  
							theme_bw() + theme(legend.key.size=unit(11,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) +
							facet_wrap( ~ group, scales='free_y', ncol=4)		
					ggsave(file=file, w=12,h=4)
				}									
			})
	#
	ylab		<- "Number of transmissions"
	run.tp		<- subset(runs.risk, method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating==method.DATING & grepl(method.RISK,method.risk)  & (grepl('N.',stat,fixed=1) | stat=='N' ) )
	stat.select	<- gsub('RI','N', stat.select)	
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])
	run.tp[, cascade.stage:=run.tp[, substr(factor, 1, 1)]]	
	set(run.tp, run.tp[,which(cascade.stage=='A')], 'cascade.stage', 'cART initiated')
	set(run.tp, run.tp[,which(cascade.stage=='U')], 'cascade.stage', 'Undiagnosed')
	set(run.tp, run.tp[,which(cascade.stage=='D' & factor%in%c("Dtl500","Dtl350"))], 'cascade.stage', 'Diagnosed\n CD4<=500')
	set(run.tp, run.tp[,which(cascade.stage=='D')], 'cascade.stage', 'Diagnosed')	
	set(run.tp, NULL, 'cascade.stage', run.tp[, factor(cascade.stage, levels=c('Undiagnosed','Diagnosed','Diagnosed\n CD4<=500','cART initiated'))])
	set(run.tp, NULL, 'factor', run.tp[, factor(factor, levels=factors[, levels(factor)])])
	run.tp	<- merge(run.tp, factors, by='factor')	
	run.tp	<- merge(run.tp, tperiod.info, by='t.period')
	run.tp[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]		
	dummy	<- lapply(seq_along(stat.select), function(i)
			{
				cat(paste('\nprocess', stat.select[i]))
				ggplot(subset(run.tp, !is.na(v) & stat%in%c(stat.select[i])), aes(x=factor, y=v, fill=factor.legend, colour=factor.legend)) + labs(x="", y=ylab) + 
						scale_y_continuous(breaks=seq(0,300,10), minor_breaks=seq(0,300,1)) + scale_x_discrete(breaks=NULL, limits=run.tp[, levels(factor)]) +
						scale_fill_manual(name='from cascade stage', values=run.tp[, unique(factor.color)], guide=with.guide) + 
						scale_colour_manual(name='from cascade stage', values = rep('black',run.tp[, length(unique(factor))]), guide=with.guide) +
						#scale_fill_brewer(palette='PRGn',name='from cascade stage') + scale_colour_manual(name='from cascade stage', values = rep('black',11)) +					
						#guides(colour=FALSE, fill = guide_legend(override.aes = list(size=5))) +
						theme_bw() + theme(legend.key.size=unit(10.5,'mm'), plot.margin=unit(c(0,0,-5,0),"mm")) + #coord_flip() +
						geom_bar(stat='identity',binwidth=1, position='dodge')	+ geom_errorbar(aes(ymin=l95.bs, ymax=u95.bs), width=0.3, position=position_dodge(width=0.9))	+ 
						facet_grid(. ~ t.period.long, margins=FALSE)
				file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_',run.tp[1, method.recentctime],'_',run.tp[1, method.denom], '_',run.tp[1, method.dating], '_',run.tp[1, method.brl],'_',run.tp[1, method.risk],'_',stat.select[i],".pdf", sep='')
				cat(paste('\nsave to file',file))							
				ggsave(file=file, w=16,h=5.5)							
			})
	
}
######################################################################################
project.athena.Fisheretal.censoring.model<- function(ct, ctb, ctn=NULL, plot.file=NA, factors=NULL)
{
	method.group	<- ifelse( ct[, any(grepl('<=100',factor2))], 'age', ifelse(ct[, any(grepl('UA',factor2))], 'cascade', 'art'))
	if(!is.null(factors))
		setnames(factors, 'factor', 'factor2')		 
	if(!is.na(plot.file))
	{		
		theme.cascade <- function (base_size = 12, base_family = "") 
		{
			theme_grey(base_size = base_size, base_family = base_family) %+replace% 
					theme(	legend.key.size=unit(11,'mm'), 
							axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5))
		}		
		theme.age <- function (base_size = 12, base_family = "") 
		{
			theme_grey(base_size = base_size, base_family = base_family) %+replace% 
					theme(	legend.key.size=unit(11,'mm'), 
							axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5), 
							panel.background = element_rect(fill = 'grey60'),
							panel.grid.minor.y = element_line(colour='grey70'),
							panel.grid.major = element_line(colour='grey70'),
							legend.key=element_rect(fill='grey60'))
		}
		if(method.group=='age')
			theme_set(theme.age())
		if(method.group!='age')
			theme_set(theme.cascade())	
	}
	#	ct.bs contains counts at times [t1bs, t2bs]
	#	divide these by the counts at times [t1,t2] --> bootstrap fraction
	set(ctb, NULL, 'c', ctb[, nc/n])
	set(ctb, ctb[, which(n==0)], 'c', 1.)
	set(ctb, NULL, 't2c', ctb[, cens.t-t.period.max.bs])
	if(!is.na(plot.file))
	{
		ct.plot		<- subset(ctb, select=c(risk, factor, factor2, c, t2c))	
		ct.plot		<- merge( ct.plot, factors, by='factor2')
		ct.plot		<- merge(ct.plot, ct.plot[, list(mc=median(c)), by='factor'], by='factor')
		if(method.group=='cascade')
			ct.plot	<- subset(ct.plot, grepl('U',factor2))
		set(ct.plot, NULL, 't2cf', ct.plot[, factor(round(t2c,d=2))])
		ggplot(ct.plot, aes(x=t2c, y=c, colour=factor.legend)) + geom_point(aes(y=mc)) + geom_boxplot(aes(position=factor(t2c))) +
				labs(x=expression(atop('time to censoring '*t[C]-t[2]^k,'( years )')), y=expression(atop('fraction of non-censored potential transmission intervals','( bootstrap samples '*hat(c)[b]*' )'))) +
				scale_colour_manual(values=ct.plot[, unique(factor.color)], guide = FALSE) +		
				theme(panel.grid.minor.x = element_blank(), legend.position = "bottom") +
				scale_y_continuous(breaks=seq(0,1,0.2), lim=c(0,1) ) +
				facet_grid(. ~ factor.legend, margins=FALSE)
		file		<- paste(plot.file,'_censoringfraction.pdf',sep='')
		if(method.group=='cascade')
			ggsave(file=file, w=10, h=6)		
	}
	#	compute adjusted numbers for each bootstrap run
	ctn		<- merge( subset(ct, stat=='X.msm', c(risk, factor, factor2, t.period, n)), subset(ctb, select=c(risk, factor, factor2, t.period, BS, c)), by=c('t.period','risk','factor','factor2'))
	ctn[, n.adj:= ctn[, n/c]]	
	#	compute sample mean and median of bootstrap fraction -- mean is OK
	ctb		<- ctb[, list(p.cens.med=median(c), p.cens=mean(c), p.cens.95l=quantile(c, p=0.025), p.cens.95u=quantile(c, p=0.975)), by=c('t.period','risk','factor')]
	ctb		<- merge(subset(ct, stat=='X.msm', c(risk, factor, factor2, t.period, n)), ctb, by=c('t.period','risk','factor'))	
	#	adjust number of potential transmission intervals
	ctb[, n.adj:= ctb[, n/p.cens]]	
	#
	if(!is.na(plot.file))
	{
		ct.plot		<- copy(ctb)
		ct.plot		<- merge(ct.plot, subset(ctn, select=c('t.period','Patient.n')), by='t.period')
		ct.plot		<- melt(ct.plot, measure.vars=c('n','n.adj'), variable.name='group', value.name='n')
		set(ct.plot, NULL, 'n', ct.plot[, n/Patient.n])
		ct.plot		<- merge( ct.plot, factors, by='factor2')
		#	add time period
		tmp			<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))		
		set(tmp, NULL, 't.period.min', tmp[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
		set(tmp, NULL, 't.period.max', tmp[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )		
		ct.plot		<- merge( ct.plot, tmp, by='t.period')
		ct.plot[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]
		#	set pretty labels
		set(ct.plot, ct.plot[,which(group=='n')], 'group', 'observed counts')
		set(ct.plot, ct.plot[,which(group=='n.adj')], 'group', 'estimate adjusted\nfor right censoring')		
		set(ct.plot, NULL, 'group', ct.plot[, factor(group, levels=c('observed counts','estimate adjusted\nfor right censoring'))])
		#	plot
		ggplot(ct.plot, aes(x=t.period.long, y=n, colour=factor.legend, shape=group)) + geom_point(size=3) +
				labs(x='', y='potential transmission intervals in cohort\n(person-years per recipient MSM)',shape='') +
				scale_colour_manual(values=ct.plot[, unique(factor.color)], guide = FALSE) +	
				scale_shape_manual(values=c(8,2)) +
				facet_grid(. ~ factor.legend, margins=FALSE) + theme(legend.position = "bottom")
		file		<- paste(plot.file,'_censoringmodel.pdf',sep='')
		ggsave(file=file, w=17, h=6)	
		#
		ct.plot		<- copy(ctb)			
		ct.plot		<- merge(ct.plot, ct.plot[, list(risk=risk, factor2=factor2, prop.adj=round(n.adj/sum(n.adj),d=4), prop=round(n/sum(n),d=4)), by='t.period'], by=c('t.period','risk','factor2'))
		ct.plot		<- melt(ct.plot, measure.vars=c('prop','prop.adj'), id.vars=c('t.period','risk','factor2','factor'), variable.name='group', value.name='prop')
		ct.plot		<- merge( ct.plot, factors, by='factor2')
		ct.plot		<- merge( ct.plot, tmp, by='t.period')
		ct.plot[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]
		set(ct.plot, ct.plot[,which(group=='prop')], 'group', 'observed proportion')
		set(ct.plot, ct.plot[,which(group=='prop.adj')], 'group', 'estimated proportion,\nadjusted for right censoring')
		set(ct.plot, NULL, 'group', ct.plot[, factor(group, levels=c('observed proportion','estimated proportion,\nadjusted for right censoring'))])
		ggplot(ct.plot, aes(x=t.period.long, y=prop, colour=factor.legend, shape=group)) + geom_point(size=3) +
				labs(x='', y='potential transmission intervals in cohort\n(% within each observation period)',shape='') +
				scale_y_continuous(breaks=seq(0,1,0.05)) +
				scale_colour_manual(values=ct.plot[, unique(factor.color)], guide = FALSE) +	
				scale_shape_manual(values=c(8,2)) +
				facet_grid(. ~ factor.legend, margins=FALSE) + theme(legend.position = "bottom")
		file		<- paste(plot.file,'_censoringmodelprop.pdf',sep='')
		ggsave(file=file, w=17, h=6)	
	}
	#
	ctb		<- subset(ctb, select=c(t.period, risk, factor, factor2, p.cens, p.cens.95l, p.cens.95u))
	ctn		<- subset(ctn, select=c(t.period, risk, factor, factor2, n.adj))
	list( ctb=ctb, ctn=ctn )
}	
######################################################################################
project.athena.Fisheretal.censoring.explore<- function()
{
	require(data.table)
	require(ape)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal",sep='/')
	outdir					<- paste(DATA,"fisheretal_140905",sep='/')		
	
	
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	outfile					<- infile
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period	
	select					<- 'Yscore3ka2H[0-9\\.]+C3V51'
	select					<- 'Yscore3ka2H[0\\.5|1|2]C3V51'
	select					<- 'Yscore3ka2H0.5C3V51|Yscore3ka2H1C3V51|Yscore3ka2H2C3V51'
	
	files					<- list.files(indir)
	files					<- files[ sapply(files, function(x) grepl('*tables.*R$',x) & grepl(select, x)) ]	
	if(!length(files))	stop('cannot find files matching criteria')
	runs.opt	<- lapply( files, function(z)
			{
				print(z)
				method.brl			<- regmatches(z, regexpr('Yscore[^_]*',z))
				method.brl			<- substr(method.brl, 7, nchar(method.brl))
				method.risk			<- regmatches(z, regexpr('[^_]*.R$',z))
				method.risk			<- substr(method.risk,1,nchar(method.risk)-2)
				method.dating		<- ifelse(grepl('sasky',z),'sasky','gmrf')
				method.recentctime	<- ifelse(grepl('2011',z),'2011','2013-03-01')
				method.nodectime	<- ifelse(grepl('a',method.brl),'any','map')
				data.table(file=z, method.brl=method.brl, method.nodectime=method.nodectime, method.dating=method.dating, method.risk=method.risk, method.recentctime=method.recentctime)				
			})
	runs.opt	<- do.call('rbind', runs.opt)
	setkey(runs.opt, method.dating, method.brl)	
	runs.opt	<- subset(runs.opt, !is.na(file))
	print(runs.opt)	
	#
	#	set up cens tables
	#
	cens.tables.bs	<- runs.opt[,	{
										tmp	<- paste(indir,'/',file,sep='')
										cat(paste('\nprocess file=',tmp))
										tmp	<- load(tmp)
										ans	<- ans$cens.table.bs
										ans[, method.risk:=method.risk]
										ans[, method.brl:=method.brl ]																		
										ans
									},by='file']
	cens.tables.bs[, file:=NULL]
	set(cens.tables.bs, NULL, 'n', cens.tables.bs[, n/8])
	set(cens.tables.bs, NULL, 'nc', cens.tables.bs[, nc/8])
	#
	cens.tables	<- runs.opt[,	{
										tmp	<- paste(indir,'/',file,sep='')
										cat(paste('\nprocess file=',tmp))
										tmp	<- load(tmp)
										ans	<- subset(ans$cens.table, select=c(stat, t.period, factor, risk, n, factor2, sum, p))										
										ans[, method.risk:=method.risk]
										ans[, method.brl:=method.brl ]																		
										ans
									},by='file']
	cens.tables[, file:=NULL]
	set(cens.tables, NULL, 'n', cens.tables[, n/8])		
	#
	cens.Patient.n	<- runs.opt[,	{
										tmp	<- paste(indir,'/',file,sep='')
										cat(paste('\nprocess file=',tmp))
										tmp	<- load(tmp)
										ans	<- ans$cens.Patient.n
										ans[, method.risk:=method.risk]
										ans[, method.brl:=method.brl ]								
										ans
									},by='file']
	cens.Patient.n[, file:=NULL]
	
	file	<- paste(outdir, '/', "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_tables.R", sep='')
	save(cens.tables.bs, cens.Patient.n, file=file)
	#"/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_140828/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_tables.R"
	#
	#	set up factor legends
	#
	factor.color	<- c("#990000","#EF6548","#FDBB84","#0570B0","#74A9CF","#7A0177","#F768A1","#FCC5C0","#005824","#41AB5D","#ADDD8E")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
			'Undiagnosed,\n Chronic infection\n at diagnosis',
			'Undiagnosed,\n Unknown if recent',
			'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
			'Diagnosed,\n CD4 progression to\n >500',
			'Diagnosed,\n CD4 progression to\n [351-500]',
			'Diagnosed,\n CD4 progression to\n <=350',			
			'Diagnosed,\n Unknown CD4',
			'cART initiated,\n no viral suppression',
			'cART initiated,\n viral suppression',
			'cART initiated,\n Unknown viral load'				)
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.suA.N","ART.suA.Y","ART.vlNA")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m2Bwmx')
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
			'Undiagnosed,\n Chronic infection\n at diagnosis',
			'Undiagnosed,\n Unknown if recent',
			'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
			'Diagnosed,\n CD4 progression to\n >500',
			'Diagnosed,\n CD4 progression to\n [351-500]',
			'Diagnosed,\n CD4 progression to\n <=350',			
			'Diagnosed,\n Unknown CD4',
			'cART initiated,\n no viral suppression',
			'cART initiated,\n viral suppression',
			'cART initiated,\n Unknown viral load'				)	
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.su.N","ART.su.Y","ART.vlNA")
	factors			<- rbind(factors, data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m2Bt'))	
	factor.long	<- c("<20","20-24","25-29","30-34","35-39","40-44","45-")
	factor.color<- brewer.pal(length(factor.long), 'RdBu')	
	tmp			<- c("t<=20","t<=25","t<=30","t<=35","t<=40","t<=45","t<=100")
	factors		<- rbind(factors, data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m5.tAc'))
	#
	#	set up t.period
	#
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	#
	#	CASCADE
	#		
	ctb			<- subset(cens.tables.bs, method.brl=='3ka2H1C3V51')
	setkey(ctb, stat, t.period, risk, factor)
	ctb			<- unique(ctb)
	ct			<- subset(cens.tables, method.brl=='3ka2H1C3V51')
	setkey(ct, stat, t.period, risk, factor)
	ct			<- unique(ct)	
	ctn			<- subset(cens.Patient.n, grepl('X.msm',stat) & method.brl=='3ka2H1C3V51')
	setkey(ctn, stat, t.period)
	ctn			<- unique(ctn)
	plot.file	<- paste(outdir, '/', outfile, '_', 'SEQ', '_',ct[1,method.brl],'_',ct[1, method.risk], sep='')
	tmp			<- subset(factors, method.risk=='m2Bwmx')
	ct1			<- project.athena.Fisheretal.censoring.model(ct, ctb, ctn, plot.file=plot.file, factors=tmp)
	
	ctb			<- subset(cens.tables.bs, method.brl=='3ka2H0.5C3V51')
	setkey(ctb, stat, t.period, risk, factor)
	ctb			<- unique(ctb)
	ct			<- subset(cens.tables, method.brl=='3ka2H0.5C3V51')
	setkey(ct, stat, t.period, risk, factor)
	ct			<- unique(ct)	
	ctn			<- subset(cens.Patient.n, grepl('X.msm',stat) & method.brl=='3ka2H0.5C3V51')
	setkey(ctn, stat, t.period)
	ctn			<- unique(ctn)
	plot.file	<- paste(outdir, '/', outfile, '_', 'SEQ', '_',ct[1,method.brl],'_',ct[1, method.risk], sep='')
	tmp			<- subset(factors, method.risk=='m2Bwmx')
	ct05		<- project.athena.Fisheretal.censoring.model(ct, ctb, ctn, plot.file=plot.file, factors=tmp)
	
	ctb			<- subset(cens.tables.bs, method.brl=='3ka2H2C3V51')
	setkey(ctb, stat, t.period, risk, factor)
	ctb			<- unique(ctb)
	ct			<- subset(cens.tables, method.brl=='3ka2H2C3V51')
	setkey(ct, stat, t.period, risk, factor)
	ct			<- unique(ct)	
	ctn			<- subset(cens.Patient.n, grepl('X.msm',stat) & method.brl=='3ka2H2C3V51')
	setkey(ctn, stat, t.period)
	ctn			<- unique(ctn)
	plot.file	<- paste(outdir, '/', outfile, '_', 'SEQ', '_',ct[1,method.brl],'_',ct[1, method.risk], sep='')
	tmp			<- subset(factors, method.risk=='m2Bwmx')
	ct2			<- project.athena.Fisheretal.censoring.model(ct, ctb, ctn, plot.file=plot.file, factors=tmp)	
	
	#
	#	AGE
	#
	file	<- '~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3kaH_tablesSEQ_m5.tAc.tp1.R'
	load(file)
	ct		<- ans$cens.table.all
	ctn		<- ans$cens.Patient.n
	setkey(ctn, stat, t.period)
	ctn		<- unique(ctn)
	tmp		<- subset(factors, method.risk=='m5.tAc')
	plot.file	<- paste(outdir, '/', outfile, '_', 'SEQ', '_','3kaH','_','m5.tAc.tp1', sep='')
	#ct.p.H2		<- project.athena.Fisheretal.censoring.model(ct, ctn, plot.file=plot.file, factors=tmp)
	#
	# illustrate extent of censoring
	#
	file			<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3ka2H1C3V51_tablesSEQ_m2Bwmx.tp1.R'
	load(file)
	cens.AnyPos_T1	<- ans$cens.AnyPos_T1
	tmp				<- subset(cens.AnyPos_T1, stage%in%c('U.2','UA.2'))
	set(tmp, tmp[, which(stage=='UA.2')], 'stage', 'Undiagnosed,\n Recent infection\n at diagnosis')
	set(tmp, tmp[, which(stage=='U.2')], 'stage', 'Undiagnosed,\n Chronic infection\n at diagnosis')
	ggplot( tmp, aes(x=t.AnyPos_T1, fill=stage)) + geom_histogram(binwidth=0.25) +
			labs(y='number of potential transmitters\nto recipient MSM diagnosed between 2006-5 and 2007-12', x='time of diagnosis of potential transmitter', fill='', colour='') +
			scale_fill_manual(values=c("#EF6548","#990000"), guide = FALSE) + geom_vline(xintercept = c(2006.408, 2008.)) +
			scale_y_continuous(breaks=seq(0,2000,20)) + scale_x_continuous(breaks=seq(2005, 2015, 0.5)) +
			facet_grid(. ~ stage, margins=FALSE)
			#facet_grid(. ~ stage, scales='free_x', space='free_x', margins=FALSE)
	plot.file	<- paste(outdir, '/', outfile, '_', 'SEQ', '_','3kaH','_','censoringUandUA.pdf', sep='')
	ggsave(file=plot.file, w=12, h=6)
	#
	#
	#
	file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_RISEQ_3kaH2_tATHENAmsm.R'
	tmp		<- load(file)
	print(tmp)
	X.msm	<- YX.part1
	#check tperiod.info
	X.msm	<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.msm, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, vl.suppressed=log10(1e3), plot.file.varyvl=NA, plot.file.or=NA )
	
	subset(X.msm, AnyPos_T1<2008.265)[, table(CD4b.tperiod)/8]
	subset(X.msm, AnyPos_T1<2009.655)[, table(CD4b.tperiod)/8]
	hist(subset(X.msm, AnyPos_T1<2008.265 & CD4b.tperiod=='U.2')[, t.AnyPos_T1])
	hist(subset(X.msm, AnyPos_T1<2009.655 & CD4b.tperiod=='U.3')[, t.AnyPos_T1])
	#
	#
	#
	ct	<- subset(cens.tables, factor2=='ART.suA.Y')
	setkey(ct, t.period, risk, factor, method.brl, stat)
	ct	<- unique(ct)
	setkey(ct, method.brl, t.period)
	subset(ct, stat=='X.clu')
	subset(ct, stat=='X.seq')	
	#
	#	debug
	#
	file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_RISEQ_3kaH_tATHENAseq.R'
	tmp		<- load(file)
	X.seq.H	<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(YX.part1, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA, plot.file.or=NA )
	file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_RISEQ_3kaH1_tATHENAseq.R'
	tmp		<- load(file)
	X.seq.H1<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(YX.part1, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA, plot.file.or=NA )
	
	
	X.seq.H[, table(CD4b.tperiod)]
	X.seq.H1[, table(CD4b.tperiod)]	
	X.seq.suy.H		<- subset(X.seq.H, CD4b.tperiod=='ART.suA.Y.1')
	X.seq.suy.H1	<- subset(X.seq.H1, CD4b.tperiod=='ART.suA.Y.1')	
	tmp	<- setdiff(X.seq.suy.H1[, unique(Patient)], X.seq.suy.H[, unique(Patient)])
	# [1] "M32793" "M32820" "M32838" "M32871" "M32876" "M32894" "M32901" "M32912" "M32921" "M32925" "M32926" "M32957" "M32958" "M32986" "M33001"
	#[16] "M33021" "M33028" "M33043" "M33046" "M33058" "M33066" "M33073" "M33080" "M33081" "M33083" "M33093" "M33095" "M33106" "M33113" "M33115"
	#[31] "M33123" "M33142" "M33147" "M33158" "M33183" "M33187" "M33227" "M33288" "M33298" "M33390" "M33425" "M33447" "M33458" "M33649" "M33736"
	#[46] "M33967" "M35966"
	nrow(subset(X.seq.suy.H1, !Patient%in%tmp))
	
	subset(X.seq.H, Patient=="M35966")[, length(unique(t.Patient))]
	subset(X.seq.H1, Patient=="M35966")[, length(unique(t.Patient))]
	
	subset(X.seq.H1, Patient=="M35966")[, table(CD4b.tperiod)]
	subset(X.seq.H, Patient=="M35966")[, table(CD4b.tperiod)]
	
	subset(X.seq.H, t.period==1)[, range(AnyPos_T1)]
	subset(X.seq.H1, t.period==1)[, range(AnyPos_T1)]
	
	X.seq.H[, length(unique(Patient))]
	X.seq.H1[, length(unique(Patient))]
	X.seq.H[, length(unique(t.Patient))]
	X.seq.H1[, length(unique(t.Patient))]
}
######################################################################################
project.athena.Fisheretal.sensitivity.gettables<- function()
{
	require(data.table)
	require(ape)
	require(RColorBrewer)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_140714",sep='/')	
	outdir					<- paste(DATA,"fisheretal_140714",sep='/')
	indir					<- paste(DATA,"fisheretal_140616",sep='/')	
	outdir					<- paste(DATA,"fisheretal_140616",sep='/')	
	indir					<- paste(DATA,"fisheretal_140929",sep='/')	
	outdir					<- paste(DATA,"fisheretal_140929",sep='/')	
	indir					<- paste(DATA,"fisheretal_141022",sep='/')	
	outdir					<- paste(DATA,"fisheretal_141022",sep='/')	
	
	
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
		
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.table.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))	
	
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )		
	#set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-5','2008-1','2009-7'), labels=c('Jul 1996','May 2006','Jan 2008','Jul 2009'))])
	#set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-4','2007-12','2009-6','2010-12'), labels=c('Apr 2006','Dec 2007','Jun 2009','Dec 2010'))])
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-5','2008-1','2009-7'), labels=c('96/07','\n\n06/05','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-4','2007-12','2009-6','2010-12'), labels=c('06/04','07/12','09/06','10/12'))])
	#
	#
	#	set up factor legends
	factor.color	<- c(	"#990000","#EF6548","#FDBB84",
							"#0C2C84","#0570B0","#74A9CF","#41B6C4","#35978F",  
							"#7A0177","#F768A1","#FCC5C0",
							"#1A9850","#A6D96A")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
							'Undiagnosed,\n Chronic infection\n at diagnosis',
							'Undiagnosed,\n Unknown if recent',
							'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
							'Diagnosed,\n CD4 progression to >500',
							'Diagnosed,\n CD4 progression to [350-500]',
							'Diagnosed,\n CD4 progression to <350',			
							'Diagnosed,\n No CD4 measured',
							'ART initiated,\n No viral suppression, >1 observations',
							'ART initiated,\n No viral suppression, 1 observation',
							'ART initiated,\n No viral load measured',			
							'ART initiated,\n Viral suppression, 1 observation',
							'ART initiated,\n Viral suppression, >1 observations')
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.suA.N2","ART.suA.N1","ART.vlNA","ART.suA.Y1","ART.suA.Y2")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m2Cwmx')	
	#factor.long	<- c("<20","20-24","25-29","30-34","35-39","40-44","45-")
	#factor.color<- brewer.pal(length(factor.long), 'RdBu')	
	#tmp			<- c("t<=20","t<=25","t<=30","t<=35","t<=40","t<=45","t<=100")
	#factors		<- rbind(factors, data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m5.tAc'))
	#
	#	get n.coh p.seq etc by time period: table
	#
	file	<- paste(outdir,'/','table_m2BwmxMv.wtn.tp','.csv',sep='')	
	df		<- subset(runs.table, grepl('m2BwmxMv.wtn.tp',method.risk) & method.brl=='3la2H1C3V51', select=c(stat, method.risk, factor, n, p, l95.bs, u95.bs))
	project.athena.Fisheretal.sensitivity.tables.m2(df, factors, levels, file)
	#
	#	get n.coh p.seq etc by time period: figure
	#
	#	set up factor legends
	factor.color	<- c("#990000","#EF6548","#FDBB84","#0570B0","#74A9CF","#7A0177","#F768A1","#FCC5C0","#005824","#41AB5D","#ADDD8E")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
							'Undiagnosed,\n Chronic infection\n at diagnosis',
							'Undiagnosed,\n Unknown if recent',
							'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
							'Diagnosed,\n CD4 progression to\n >500',
							'Diagnosed,\n CD4 progression to\n [351-500]',
							'Diagnosed,\n CD4 progression to\n <=350',			
							'Diagnosed,\n Unknown CD4',
							'cART initiated,\n no viral suppression',
							'cART initiated,\n viral suppression',
							'cART initiated,\n Unknown viral load'				)
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.suA.N","ART.suA.Y","ART.vlNA")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m2Bwmx')		
	factor.color	<- c("#990000","#EF6548","#FDBB84","#0570B0","#74A9CF","#7A0177","#F768A1","#FCC5C0","#005824","#41AB5D","#F7FCB9","#ADDD8E")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
							'Undiagnosed,\n Chronic infection\n at diagnosis',
							'Undiagnosed,\n Unknown if recent',
							'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
							'Diagnosed,\n CD4 progression to\n >500',
							'Diagnosed,\n CD4 progression to\n [351-500]',
							'Diagnosed,\n CD4 progression to\n <=350',			
							'Diagnosed,\n Unknown CD4',
							'cART initiated,\n no viral suppression',
							'cART initiated,\n viral suppression <2 obs',
							'cART initiated,\n viral suppression >=2 obs',
							'cART initiated,\n Unknown viral load'				)
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.suA.N","ART.suA.Y1","ART.suA.Y2","ART.vlNA")
	factors			<- rbind(factors, data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m2Cwmx'))
	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1C3V51'
	method.RISK		<- 'm2BwmxMv.wtn.tp'
	method.WEIGHT	<- ''	
	factors			<- subset(factors, grepl('m2Bwmx',method.risk), select=c(factor, factor.legend, factor.color))
	outfile			<- infile	
	project.athena.Fisheretal.sensitivity.getfigures.counts(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, 'cascade')
	#
	#		m2Cwmx
	#
	factor.color	<- c(	"#990000","#EF6548","#FDBB84",
			"#0C2C84","#0570B0","#74A9CF","#41B6C4","#35978F",  
			"#7A0177","#F768A1","#FCC5C0",
			"#1A9850","#A6D96A")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
			'Undiagnosed,\n Chronic infection\n at diagnosis',
			'Undiagnosed,\n Unknown if recent',
			'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
			'Diagnosed,\n CD4 progression to >500',
			'Diagnosed,\n CD4 progression to [350-500]',
			'Diagnosed,\n CD4 progression to <350',			
			'Diagnosed,\n No CD4 measured',
			'ART initiated,\n No viral suppression, >1 observations',
			'ART initiated,\n No viral suppression, 1 observation',
			'ART initiated,\n No viral load measured',			
			'ART initiated,\n Viral suppression, 1 observation',
			'ART initiated,\n Viral suppression, >1 observations')
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.suA.N2","ART.suA.N1","ART.vlNA","ART.suA.Y1","ART.suA.Y2")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m2Cwmx')
	#	
	#	updated stages
	#	set up factor legends
	factor.color	<- c(	"#990000","#EF6548","#FDBB84",
			"#0C2C84","#0570B0","#74A9CF","#41B6C4","#35978F",  
			"#FCC5C0","#F768A1","#7A0177",
			"#1A9850","#A6D96A")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
			'Undiagnosed,\n Chronic infection\n at diagnosis',
			'Undiagnosed,\n Unknown if recent',
			'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
			'Diagnosed,\n CD4 progression to >500',
			'Diagnosed,\n CD4 progression to [350-500]',
			'Diagnosed,\n CD4 progression to <350',			
			'Diagnosed,\n No CD4 measured',
			'ART initiated,\n Before first viral suppression',													
			'ART initiated,\n After first viral suppression\n No viral load measured',		
			'ART initiated,\n After first viral suppression\n No viral suppression',	
			'ART initiated,\n After first viral suppression\n Viral suppression, 1 observation',
			'ART initiated,\n After first viral suppression\n Viral suppression, >1 observations')
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.NotYetFirstSu","ART.vlNA","ART.suA.N","ART.suA.Y1","ART.suA.Y2")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m2Cwmx')	

	#
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3la2H1C3V100'
	method.RISK		<- 'm2CwmxMv.wtn.tp'
	method.WEIGHT	<- ''	
	factors			<- subset(factors, grepl('m2Cwmx',method.risk), select=c(factor, factor.legend, factor.color))
	outfile			<- infile	
	project.athena.Fisheretal.sensitivity.getfigures.npotentialtransmissionintervals.adjusted(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, 'cascade')
	project.athena.Fisheretal.sensitivity.getfigures.nlikelytransmissionintervals(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, 'cascade')
	project.athena.Fisheretal.sensitivity.getfigures.comparetransmissionintervals(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, 'cascade')
	
	#	set up factor legends
	factor.color	<- c(	"#990000","#EF6548","#FDBB84",
							"#0C2C84","#0570B0","#74A9CF","#41B6C4","#35978F",  
							"#7A0177","#F768A1","#FCC5C0",
							"#1A9850","#A6D96A")
	factor.long		<- c(	'Undiagnosed,\n Recent infection\n at diagnosis',	
							'Undiagnosed,\n Chronic infection\n at diagnosis',
							'Undiagnosed,\n Unknown if recent',
							'Diagnosed < 3mo,\n Recent infection\n at diagnosis',					
							'Diagnosed,\n CD4 progression\n to >500',
							'Diagnosed,\n CD4 progression\n to [350-500]',
							'Diagnosed,\n CD4 progression\n to <350',			
							'Diagnosed,\n No CD4 measured',
							'ART initiated,\n No viral suppression,\n >1 observations',
							'ART initiated,\n No viral suppression,\n 1 observation',
							'ART initiated,\n No viral load measured',			
							'ART initiated,\n Viral suppression,\n 1 observation',
							'ART initiated,\n Viral suppression,\n >1 observations')
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.suA.N2","ART.suA.N1","ART.vlNA","ART.suA.Y1","ART.suA.Y2")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m2Cwmx')		
	project.athena.Fisheretal.sensitivity.getfigures.pseq(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, 'cascade')
	
	
	project.athena.Fisheretal.sensitivity.getfigures.counts(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, 'cascade')
	
	
	
	
	

	#
	#	get table summing over time
	#
	levels	<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.suA.N","ART.suA.Y","ART.vlNA")
	df		<- subset(runs.table, grepl('m2BwmxMv.tp',method.risk) & !grepl('wstar',method.risk), c(stat, method.risk, factor, n, p))	
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	set(df, NULL, 'factor', df[, factor(factor, levels=levels)])
	tmp		<- df[, list( n.coh=sum(n[stat=='X.msm'])/8, p.seq= sum(n[stat=='X.seq']) / sum(n[stat=='X.msm']), p.clu= sum(n[stat=='X.clu']) / sum(n[stat=='X.seq']), p.pos= sum(n[stat=='YX']) / sum(n[stat=='X.clu']), n.pos=sum(n[stat=='YX'])/8 ), by=c('factor')]
	set(tmp, NULL, 'pp.coh', tmp[, n.coh/sum(n.coh)])
	set(tmp, NULL, 'pp.pos', tmp[, n.pos/sum(n.pos)])	
	setkey(tmp, factor)
	ans		<- tmp
	print(ans[, sum(n.coh)])	
	levels	<- c('t<=20','t<=25','t<=30','t<=35','t<=40','t<=45','t<=100')
	df		<- subset(runs.table, grepl('m5.tAc.tp',method.risk) & !grepl('wstar',method.risk), c(stat, method.risk, factor, n, p))	
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	set(df, NULL, 'factor', df[, factor(factor, levels=levels)])	
	tmp		<- df[, list( n.coh=sum(n[stat=='X.msm'])/8, p.seq= sum(n[stat=='X.seq']) / sum(n[stat=='X.msm']), p.clu= sum(n[stat=='X.clu']) / sum(n[stat=='X.seq']), p.pos= sum(n[stat=='YX']) / sum(n[stat=='X.clu']), n.pos=sum(n[stat=='YX'])/8 ), by=c('factor')]
	set(tmp, NULL, 'pp.coh', tmp[, n.coh/sum(n.coh)])
	set(tmp, NULL, 'pp.pos', tmp[, n.pos/sum(n.pos)])	
	setkey(tmp, factor)
	ans		<- rbind(ans, tmp)	
	levels	<- c("ART.I", "ART.F","ART.T","ART.P","ART.pulse.Y","ART.l3","ART.g3","ART.3.2NRT.X","ART.3.OTH")                                          	              
	df		<- subset(runs.table, 'm3.n3mx'==method.risk & !grepl('wstar',method.risk), c(stat, method.risk, factor, n, p))
	df		<- subset(df, factor%in%levels)
	set(df, NULL, 'factor', df[, factor(factor, levels=levels)])
	tmp		<- df[, list( n.coh=sum(n[stat=='X.msm'])/8, p.seq= sum(n[stat=='X.seq']) / sum(n[stat=='X.msm']), p.clu= sum(n[stat=='X.clu']) / sum(n[stat=='X.seq']), p.pos= sum(n[stat=='YX']) / sum(n[stat=='X.clu']), n.pos=sum(n[stat=='YX'])/8 ), by=c('factor')]
	set(tmp, NULL, 'pp.coh', tmp[, n.coh/sum(n.coh)])
	set(tmp, NULL, 'pp.pos', tmp[, n.pos/sum(n.pos)])
	setkey(tmp, factor)
	ans		<- rbind(ans, tmp)
	
	set(ans, NULL, 'p.seq', ans[, round(p.seq, d=3)]*100)
	set(ans, NULL, 'p.clu', ans[, round(p.clu, d=3)]*100)
	set(ans, NULL, 'p.pos', ans[, round(p.pos, d=4)]*100)
	set(ans, NULL, 'n.pos', ans[, round(n.pos, d=1)])
	set(ans, NULL, 'n.coh', ans[, round(n.coh, d=1)])	
	set(ans, NULL, 'pp.pos', ans[, round(pp.pos, d=3)]*100)
	set(ans, NULL, 'pp.coh', ans[, round(pp.coh, d=3)]*100)
	
	set(ans, NULL, 'pp.coh', ans[, paste('=\"(',pp.coh,')\"',sep='')])
	set(ans, NULL, 'pp.pos', ans[, paste('=\"(',pp.pos,')\"',sep='')])
	write.csv(subset(ans, select=c(factor, n.coh, pp.coh, p.seq, p.clu, p.pos, n.pos, pp.pos)), file=file, eol="\r\n", row.names=FALSE)
	
	tmp		<- df[, list( n.coh=n[stat=='X.msm']/8, p.seq= n[stat=='X.seq'] / n[stat=='X.msm'], p.clu= n[stat=='X.clu'] / n[stat=='X.seq'], p.pos= n[stat=='YX'] / n[stat=='X.clu'], n.pos=n[stat=='YX']/8 ), by=c('factor','t.period')]
	setkey(tmp, factor, t.period)
	
}
######################################################################################
project.athena.Fisheretal.sensitivity<- function()
{
	require(data.table)
	require(ape)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_140714",sep='/')
	outdir					<- paste(DATA,"fisheretal_140714",sep='/')
	indir					<- paste(DATA,"fisheretal_140729",sep='/')
	outdir					<- paste(DATA,"fisheretal_140729",sep='/')	
	#indir					<- paste(DATA,"fisheretal_140616",sep='/')
	#outdir					<- paste(DATA,"fisheretal_140616",sep='/')	
	indir					<- paste(DATA,"fisheretal_140801",sep='/')
	outdir					<- paste(DATA,"fisheretal_140801",sep='/')	
	indir					<- paste(DATA,"fisheretal_140804",sep='/')
	outdir					<- paste(DATA,"fisheretal_140804",sep='/')	
	indir					<- paste(DATA,"fisheretal_140814",sep='/')
	outdir					<- paste(DATA,"fisheretal_140814",sep='/')	
	indir					<- paste(DATA,"fisheretal_140817",sep='/')
	outdir					<- paste(DATA,"fisheretal_140817",sep='/')	
	indir					<- paste(DATA,"fisheretal_140828",sep='/')
	outdir					<- paste(DATA,"fisheretal_140828",sep='/')		
	indir					<- paste(DATA,"fisheretal_140905",sep='/')
	outdir					<- paste(DATA,"fisheretal_140905",sep='/')		
	indir					<- paste(DATA,"fisheretal_140929",sep='/')
	outdir					<- paste(DATA,"fisheretal_140929",sep='/')		
	indir					<- paste(DATA,"fisheretal_141022",sep='/')
	outdir					<- paste(DATA,"fisheretal_141022",sep='/')		
	indir					<- paste(DATA,"fisheretal_141108",sep='/')
	outdir					<- paste(DATA,"fisheretal_141108",sep='/')		
	
	
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	outfile					<- infile
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
									
	if(resume)
	{
		options(show.error.messages = FALSE)		
		file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.R", sep='')		
		readAttempt			<- try(suppressWarnings(load(file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))		
		file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.fit.R", sep='')
		readAttempt			<- try(suppressWarnings(load(file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))		
		file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.table.R", sep='')
		readAttempt			<- try(suppressWarnings(load(file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))				
		options(show.error.messages = TRUE)		
	}
	if(!resume)
	{
		files					<- list.files(indir)
		files					<- files[ sapply(files, function(x) grepl('.R$',x) ) ]	
		if(!length(files))	stop('cannot find files matching criteria')
		runs.opt	<- lapply( files, function(z)
				{
					print(z)
					method.brl			<- regmatches(z, regexpr('Yscore[^_]*',z))
					method.brl			<- substr(method.brl, 7, nchar(method.brl))
					method.denom		<- regmatches(z, regexpr('denom[[:alnum:]]*',z))
					method.denom		<- substr(method.denom, 6, nchar(method.denom))
					method.risk			<- regmatches(z, regexpr('[^_]*.R$',z))
					method.risk			<- substr(method.risk,1,nchar(method.risk)-2)
					method.dating		<- ifelse(grepl('sasky',z),'sasky','gmrf')
					method.recentctime	<- ifelse(grepl('2011',z),'2011','2013-03-01')
					method.nodectime	<- ifelse(grepl('a',method.brl),'any','map')
					data.table(file=z, method.brl=method.brl, method.nodectime=method.nodectime, method.dating=method.dating, method.risk=method.risk, method.denom=method.denom, method.recentctime=method.recentctime)				
				})
		runs.opt	<- do.call('rbind', runs.opt)
		setkey(runs.opt, method.dating, method.brl)	
		runs.opt	<- subset(runs.opt, !is.na(file))
		print(runs.opt)
		#	load risk estimates
		tmp			<- lapply(seq_len(nrow(runs.opt)), function(i)
				{
					tmp	<- paste(indir, runs.opt[i,file], sep='/')
					cat(paste('\nprocess file=',runs.opt[i,file]))
					tmp	<- load(tmp)
					ans	<- ans$risk
					if(!any(colnames(ans)=='t.period'))
						ans[, t.period:= NA_character_]										
					set(ans, NULL, 'factor', ans[, as.character(factor)])											
					ans[, method.risk:=runs.opt[i,method.risk]]
					ans[, method.dating:=runs.opt[i,method.dating]]
					ans[, method.nodectime:=runs.opt[i,method.nodectime]]
					ans[, method.brl:=runs.opt[i,method.brl ]]
					ans[, method.denom:=runs.opt[i,method.denom]]
					ans[, method.recentctime:=runs.opt[i,method.recentctime ]]
					ans
				})
		runs.risk	<- do.call('rbind', tmp)
		file			<- paste(indir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.Rdata", sep='')
		save(runs.risk, file=file)
						
		#	load risk tables
		tmp			<- lapply(seq_len(nrow(runs.opt)), function(i)
				{
					tmp	<- paste(indir, runs.opt[i,file], sep='/')
					cat(paste('\nprocess file=',runs.opt[i,file]))
					tmp	<- load(tmp)					
					ans	<- ans$risk.table
					if(!any(colnames(ans)=='p'))
						ans	<- merge(ans, ans[, list(risk=risk, factor=factor, p= n/sum(n)  ), by='stat'], by=c('risk','factor','stat'))
					ans[, method.risk:=runs.opt[i,method.risk]]
					ans[, method.dating:=runs.opt[i,method.dating]]
					ans[, method.nodectime:=runs.opt[i,method.nodectime]]
					ans[, method.brl:=runs.opt[i,method.brl ]]
					ans[, method.denom:=runs.opt[i,method.denom]]
					ans[, method.recentctime:=runs.opt[i,method.recentctime ]]
					ans
				})
		runs.table	<- do.call('rbind', tmp)
		set(runs.table, NULL, 'n', runs.table[, n*t.period])
		#	add number recipient
		tmp			<- subset( runs.risk, stat%in%c('nRec','nRecLkl'), select=c(stat, risk, factor, v, method.risk, method.dating, method.nodectime,   method.brl, method.denom, method.recentctime) )
		setnames(tmp, 'v','n')
		runs.table	<- rbind(runs.table, tmp, fill=TRUE)
		#	add expected missing
		tmp			<- subset( runs.risk, stat%in%c('Sx.e0cp'), select=c(stat, risk, factor, v, l95.bs,  u95.bs, method.risk, method.dating, method.nodectime,   method.brl, method.denom, method.recentctime) )
		setnames(tmp, 'v','n')
		runs.table	<- rbind(runs.table, tmp, fill=TRUE)
		#	add adjusted potential transmission intervals
		tmp			<- subset( runs.risk, stat%in%c('X.msm.e0cp'), select=c(stat, risk, factor, v, l95.bs,  u95.bs, method.risk, method.dating, method.nodectime,   method.brl, method.denom, method.recentctime) )
		setnames(tmp, 'v','n')
		runs.table	<- rbind(runs.table, tmp, fill=TRUE)
		#	add median likelihood of direct HIV transmission
		tmp			<- lapply(seq_len(nrow(runs.opt)), function(i)
				{
					tmp	<- paste(indir, runs.opt[i,file], sep='/')
					cat(paste('\nprocess file=',runs.opt[i,file]))
					load(tmp)
					ans	<- ans$YX[, list(p= median(score.Y)), by='stage']
					setnames(ans, 'stage', 'factor')
					ans[, risk:='stage']
					ans[, stat:='LkL']
					ans[, method.risk:=runs.opt[i,method.risk]]
					ans[, method.dating:=runs.opt[i,method.dating]]
					ans[, method.nodectime:=runs.opt[i,method.nodectime]]
					ans[, method.brl:=runs.opt[i,method.brl ]]
					ans[, method.denom:=runs.opt[i,method.denom]]
					ans[, method.recentctime:=runs.opt[i,method.recentctime ]]
					ans
				})
		tmp			<- do.call('rbind', tmp)
		runs.table	<- rbind(runs.table, tmp, fill=TRUE)
		file			<- paste(indir, '/', infile, '_', gsub('/',':',insignat), '_', "method.table.Rdata", sep='')
		save(runs.table, file=file)				
	}
	
	file<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_Wed_Dec_18_11:37:00_2013_Yscore3da_model2_m2Bwmx.tp1.cens.R"
	load(paste(indir, file, sep='/'))
	#
	#	MODEL 2
	#
	#	check cascade RR
	tmp	<- subset(runs.risk, grepl('m2',method.risk,fixed=1) & stat=='RR' & factor.ref!='None' & factor%in%c('ART.su.N','ART1.su.N','ART.suA.N'))
	setkey(tmp, method.risk, method.brl, method.dating)
	tmp	
	#	check cascade OR
	tmp	<- subset(runs.risk, method.risk%in%c("m21st.cas","m2wmx.cas","m2t.cas") & stat=='OR' & factor.ref!='None' & factor%in%c('ART.su.N','ART1.su.N','ART.suA.N'))
	setkey(tmp, method.risk, method.brl, method.dating)
	tmp	
	#	check cascade P
	#	no huge variation in sensitivity
	tmp	<- subset(runs.risk, grepl('m2',method.risk,fixed=1) & !grepl('tp',method.risk,fixed=1) & stat=='P' )
	tmp[,  list(cv= sd(v, na.rm=TRUE)/mean(v, na.rm=TRUE), mean=mean(v, na.rm=TRUE) ) , by='factor']
	#	check cascade m2mww.cas
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2wmx.cas' & stat=='RR')
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2wmx.cas' & stat=='P')
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2wmx.cas' & stat=='RI')	
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2wmx.cas' & stat=='PY')
	tmp[, PYpc:= v/sum(v, na.rm=TRUE)]
	#	check cascade m2mww.cas.adj
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2wmx.cas.adj' & stat=='RR')
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2wmx.cas.adj' & stat=='P')
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2wmx.cas.adj' & stat=='RI')	
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2wmx.cas.adj' & stat=='PY')
	tmp[, PYpc:= v/sum(v, na.rm=TRUE)]
	#
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & grepl('m2wmx',method.risk,fixed=1) & stat=='RR')
	tmp	<- subset(runs.risk, grepl('m2',method.risk,fixed=1) & stat=='RR' & factor.ref!='None' & factor%in%c('ART.su.N','ART1.su.N','ART.suA.N'))
	setkey(tmp, method.risk, method.brl, method.dating)
	tmp		
	pPWYI.UAy.X.seq<- subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk=='m2wmx.cas' & stat=='PY' & factor=='UAy')[,v] / subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk=='m2wmx.cas' & stat=='PY' )[, sum(v)]
	pPWYI.UAy.X.clu<- subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk=='m2wmx.cas.clu' & stat=='PY' & factor=='UAy')[,v] / subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk=='m2wmx.cas.clu' & stat=='PY' )[, sum(v)]
	pPWYI.UAy.X.clu.adj<- subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk=='m2wmx.cas.clu.adj' & stat=='PY' & factor=='UAy')[,v] / subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk=='m2wmx.cas.clu.adj' & stat=='PY' )[, sum(v)]	
	pPWYI.UAy.X.seq<- subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk=='m2wmx.cas' & stat=='PY' & factor=='UAy')[,v] / subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk=='m2wmx.cas' & stat=='PY' )[, sum(v)]
	pPWYI.UAy.X.clu<- subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk=='m2wmx.cas.clu' & stat=='PY' & factor=='UAy')[,v] / subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk=='m2wmx.cas.clu' & stat=='PY' )[, sum(v)]	
	pPWYI.SuAY.X.seq<- subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk=='m2wmx.cas' & stat=='PY' & factor=='ART.suA.Y')[,v] / subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk=='m2wmx.cas' & stat=='PY' )[, sum(v)]
	#	compare fit of cascade models with aic		 															 															 
	tmp	<- subset(runs.fit, 	method.risk%in%c("m21st.cas","m2wmx.cas","m2t.cas") )
	tmp	<- subset(runs.fit, 	method.risk%in%c("m21st.cas","m2wmx.cas","m2t.cas") )
	tmp	<- subset(runs.fit, 	method.risk%in%c("m2B1stMv.cas.censp","m2BwmxMv.cas.censp","m2BtMv.cas.censp") & method.recentctime=='2011')
	tmp	<- subset(runs.risk, 	method.risk%in%c("m2B1stMv.cas.censp","m2BwmxMv.cas.censp","m2BtMv.cas.censp") & method.recentctime=='2011' & stat=='RR.term' & grepl('ART',factor) & grepl('ART',factor.ref))
	
	setkey(tmp, method.brl, method.dating, method.risk)
	#	compare fit of models with reduced RR through cascade
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2B1st.cas' & stat=='RR')
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2Bt.cas' & stat=='RR')
	tmp
	#	cascade with m2wm.cas
	subset(runs.risk, method.risk=="m2wmx.cas" & stat=='RR' & factor.ref=='ART.suA.Y' & method.brl=='3da')
	subset(runs.risk, method.risk=="m2wmx.cas" & stat=='P' & factor.ref=='None' & method.brl=='3da')
	subset(runs.risk, method.risk=="m2wmx.cas" & stat=='RI' & factor.ref=='None' & method.brl=='3da')
	#
	#	MODEL 2B
	#
	#	check cascade RR
	tmp	<- subset(runs.risk, grepl('m2B',method.risk,fixed=1) & stat=='RR' & factor.ref!='None' & factor%in%c('ART.su.N','ART1.su.N','ART.suA.N'))
	setkey(tmp, method.risk, method.brl, method.dating)
	tmp	
	#	compare fit of cascade models with aic	
	tmp	<- subset(runs.fit, 	grepl('m2B',method.risk,fixed=1) & grepl('cas',method.risk,fixed=1) & method.recentctime=='2011' )
	tmp	<- subset(runs.fit, 	method.risk%in%c('m2BwmxMv.cas.censp','m2B1stMv.cas.censp','m2BtMv.cas.censp') & method.recentctime=='2011' )	
	tmp[, dAIC:= AIC-min(AIC)]
	setkey(tmp, method.brl, method.dating, dAIC)
	tmp
	#AIC   logLik        method.risk method.dating method.nodectime method.brl method.recentctime      dAIC
	#1: -646.9299 347.4649 m2BwmxMv.cas.censp         sasky              any        3da               2011  0.000000
	#2: -644.7719 346.3860   m2BtMv.cas.censp         sasky              any        3da               2011  2.157926
	#3: -631.2340 338.6170 m2B1stMv.cas.censp         sasky              any        3da               2011 15.695878	
	tmp	<- subset(runs.fit, 	method.risk%in%c('m2Bwmx.cas.censp','m2B1st.cas.censp','m2Bt.cas.censp') & method.recentctime=='2011' )	
	tmp[, dAIC:= AIC-min(AIC)]
	setkey(tmp, method.brl, method.dating, dAIC)
	tmp
	#AIC   logLik      method.risk method.dating method.nodectime method.brl method.recentctime       dAIC
	#1: -1278.502 652.2511 m2Bwmx.cas.censp         sasky              any        3da               2011  0.0000000
	#2: -1277.646 651.8229   m2Bt.cas.censp         sasky              any        3da               2011  0.8563968
	#3: -1252.392 638.1961 m2B1st.cas.censp         sasky              any        3da               2011 26.1100020
	
	#	check P cascade m2Bwmx.cas vs m2Bwmx.cas.adj for 3ca
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk%in%c('m2Bwmx.cas','m2Bwmx.cas.adj') & stat=='P')
	#	check P cascade m2Bwmx.cas across 3ca 3da gmrf sasky
	#	3ca vs 3da 	more sensitive than gmrf vs sasky
	tmp	<- subset(runs.risk, method.nodectime=='any'  & method.risk=='m2Bwmx.cas' & stat=='P')
		
	
	#
	#	risk ratio compared to undiagnosed
	#	
	run.tp		<- subset(runs.risk, method.denom==method.deno & method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m2BwmxMv.tp',method.risk) & !grepl(method.clu,method.risk) & !grepl('wstar',method.risk) & !grepl('now',method.risk) & (grepl('RR.',stat,fixed=1) | stat=='RR') )
	run.tp		<- subset(run.tp, stat%in%c("RR.term","RR.term.ptx","RR.bias.e0cp","RR.rawbias.e0cp"))
	run.tp		<- subset(run.tp, grepl('U',factor.ref), select=c('stat','risk','factor','risk.ref','factor.ref','v','l95.bs','u95.bs'))
	setkey(run.tp, factor.ref, factor)
	set(run.tp, NULL, 'v', run.tp[, round(v, d=3)])
	set(run.tp, NULL, 'l95.bs', run.tp[, round(l95.bs, d=3)])
	set(run.tp, NULL, 'u95.bs', run.tp[, round(u95.bs, d=3)])	
	
	
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])
	run.tp[, cascade.stage:=run.tp[, substr(factor, 1, 1)]]	
	set(run.tp, run.tp[,which(cascade.stage=='A')], 'cascade.stage', 'cART initiated')
	set(run.tp, run.tp[,which(cascade.stage=='U')], 'cascade.stage', 'Undiagnosed')
	set(run.tp, run.tp[,which(cascade.stage=='D' & factor%in%c("Dtl500","Dtl350"))], 'cascade.stage', 'Diagnosed\n CD4<=500')
	set(run.tp, run.tp[,which(cascade.stage=='D')], 'cascade.stage', 'Diagnosed')	
	set(run.tp, NULL, 'cascade.stage', run.tp[, factor(cascade.stage, levels=c('Undiagnosed','Diagnosed','Diagnosed\n CD4<=500','cART initiated'))])	
	run.tp	<- merge(run.tp, tmp, by='factor')	
	run.tp	<- merge(run.tp, tperiod.info, by='t.period')
	run.tp[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]
	
	
	subset( runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m2BwmxMv.tp',method.risk) & !grepl('wstar',method.risk) & grepl('clu',method.risk) & stat=='RR.term', c(factor.ref, factor, v, l95.bs, u95.bs, method.brl, method.risk) )	
	#
	#
	#
	tmp			<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m2Bwmx.tp',method.risk) & grepl('P',stat), c(method.risk,method.recentctime))
	setkey(tmp,method.risk,method.recentctime)
	tmp			<- unique(tmp)
	tmp[, postfix:=substr(method.risk, regexpr('tp',method.risk)+4, nchar(method.risk))]
	setkey(tmp,postfix,method.recentctime)
	run.tp.opt	<- unique(tmp)
	#	debug
	#file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_140502/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3da_model2_m2Bwmx.cas.censp.R'
	file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_140502/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3da_model2_m2Bwmx.tp1.censp.R'
	load(file)
	#	
	#	PYIW trends
	#
	run.tp	<- subset(runs.table, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m2Bwmx.tp',method.risk) & grepl('adj',method.risk) & !grepl('clu',method.risk) & !grepl('cens',method.risk), c(risk, factor, stat, n))	
	setkey(run.tp, factor)		
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	run.tp	<- merge( run.tp, run.tp[, list(risk=risk, factor=factor, prop= n/sum(n, na.rm=TRUE)), by=c('stat','t.period')], by=c('risk','factor','stat','t.period') )
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])
	run.tp[, cascade.stage:=run.tp[, substr(factor, 1, 1)]]
	set(run.tp, run.tp[,which(cascade.stage=='A')], 'cascade.stage', 'cART initiated')
	set(run.tp, run.tp[,which(cascade.stage=='U')], 'cascade.stage', 'Undiagnosed')
	set(run.tp, run.tp[,which(cascade.stage=='D')], 'cascade.stage', 'Diagnosed')
	set(run.tp, run.tp[,which(stat=='X')], 'stat', 'study population')
	set(run.tp, run.tp[,which(stat=='X.msm')], 'stat', 'cohort')
	set(run.tp, run.tp[,which(stat=='YX')], 'stat', 'potential transmitter')
	set(run.tp, NULL, 'stat', run.tp[, factor(stat, levels=c('potential transmitter','study population','cohort'))])
	set(run.tp, NULL, 'cascade.stage', run.tp[, factor(cascade.stage, levels=c('Undiagnosed','Diagnosed','cART initiated'))])
	tmp		<- c(	'Undiagnosed,\nEvidence for acute infection\nat diagnosis',
			'Undiagnosed,\nEvidence for recent infection\nat diagnosis',
			'Undiagnosed,\nNo recent infection',
			'Diagnosed < 3mo,\nEvidence for acute infection\nat diagnosis',
			'Diagnosed < 3mo,\nEvidence for recent infection\nat diagnosis',
			'Diagnosed,\nlowest CD4 to date\n> 500',
			'Diagnosed,\nlowest CD4 to date\n> 350',
			'Diagnosed,\nlowest CD4 to date\n<= 350',
			'Diagnosed,\nmissing 1st CD4',
			'ART initiated,\nVL missing',
			'ART initiated,\nVL continually\nsuppressed', 
			'ART initiated,\nVL not continually\nsuppressed')
	tmp		<- data.table( factor.legend= factor(tmp), factor=c("UAy","UAm","U","DAy","DAm","Dtg500","Dtl500","Dtl350","Dt.NA","ART.vlNA","ART.suA.Y","ART.suA.N"))
	run.tp	<- merge(run.tp, tmp, by='factor')	
	run.tp	<- merge( run.tp, run.tp[, list(sum=sum(n)),by=c('stat','t.period')],by=c('stat','t.period') )	
	#	PYIW adjusted for censoring
	run.tp[, na:= n] 
	run.tp[, propa:= prop] 
	for(z in run.tp[, unique(t.period)])
	{
		tmp		<- rbind(	unadjusted=sapply(c('U','UAy','UAm'), function(f2)	tmp2	<- subset(run.tp, factor==f2 & stat=='cohort' & t.period==z)[, propa ]),
							adjusted=sapply(c('U','UAy','UAm'), function(f2)	tmp2	<- subset(run.tp, factor==f2 & stat=='cohort' & t.period%in%c('2',z))[, max(propa) ])		)						
		tmp		<- run.tp[which(stat=='cohort' & t.period==z), sum[1] * tmp['adjusted',] * (1-sum(tmp['unadjusted',]))/(1-sum(tmp['adjusted',])) ]
		for(f in c('U','UAy','UAm'))
			set(run.tp, run.tp[, which(factor==f & stat=='cohort' & t.period==z)], 'na',  tmp[f])				
	}
	run.tp[, propa:= NULL]
	#for(z in run.tp[, unique(t.period)])
	#{		
	#	set(run.tp, run.tp[, which(factor=='U' & stat=='cohort' & t.period==z)], 'na',  run.tp[which(factor=='U' & stat=='cohort' & t.period%in%c('2',z)), max(na)])
	#	set(run.tp, run.tp[, which(factor=='UAy' & stat=='cohort' & t.period==z)], 'na',  run.tp[which(factor=='UAy' & stat=='cohort' & t.period%in%c('2',z)), max(na)])
	#	set(run.tp, run.tp[, which(factor=='UAm' & stat=='cohort' & t.period==z)], 'na',  run.tp[which(factor=='UAm' & stat=='cohort' & t.period%in%c('2',z)), max(na)])		
	#}
	#	
	#	Figures adjust for censoring	PYIW trends in %
	#
	run.tp	<- merge( run.tp, run.tp[, list(factor=factor, propa= na/sum(na, na.rm=TRUE)), by=c('t.period','stat')], by=c('t.period','stat','factor'))
	require(RColorBrewer)
	require(grid)		
	#	PYIW trends in # with adjustment for censoring
	ggplot(subset(run.tp, !is.na(n)), aes(x=t.period, y=n, colour=factor.legend, group=factor.legend)) + labs(x="calendar time periods", y="PYIW") + 
			scale_colour_manual(name='Transmission times', values = colorRampPalette(brewer.pal(9, "Set1"),interpolate='spline',space="Lab")(run.tp[,length(unique(factor.legend))])) +			
			theme(legend.key.size=unit(13,'mm')) + guides(colour = guide_legend(override.aes = list(size=5))) + 
			geom_point(data=subset(run.tp, !is.na(n) & n!=na), aes(y=na), shape=8, show_guide= FALSE) + geom_point() +  geom_line(aes(linetype=factor.legend), show_guide= FALSE) + facet_grid(stat ~ cascade.stage, scales='free', margins=FALSE) 
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "cens_m2wmx.tp_3da_PYIWtotal_PU_adjcens.pdf", sep='')
	ggsave(file=file, w=8,h=8)
	#	PYIW trends in % with adjustment for censoring
	ggplot(subset(run.tp, !is.na(n)), aes(x=t.period, y=prop, colour=factor.legend, group=factor.legend)) + labs(x="calendar time periods", y="% PYIW") + 
			scale_colour_manual(name='Transmission times', values = colorRampPalette(brewer.pal(9, "Set1"),interpolate='spline',space="Lab")(run.tp[,length(unique(factor.legend))])) +			
			theme(legend.key.size=unit(13,'mm')) + guides(colour = guide_legend(override.aes = list(size=5))) + 			
			geom_point(data=subset(run.tp, !is.na(n) & prop!=propa), aes(y=propa), shape=8, show_guide= FALSE) + geom_point()  + geom_line(aes(linetype=factor.legend), show_guide= FALSE) + facet_grid(stat ~ cascade.stage, margins=FALSE)
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "cens_m2wmx.tp_3da_PYIWprop_PU_adjcens.pdf", sep='')
	ggsave(file=file, w=8,h=8)	
	#	PYIW Margin in # before adjustment	
	ggplot(subset(run.tp, !is.na(n)), aes(x=t.period, y=n, fill=factor.legend)) + labs(x="calendar time periods", y="PYIW") +
			scale_fill_manual(name='Transmission times', values = colorRampPalette(brewer.pal(9, "Set1"),interpolate='spline',space="Lab")(run.tp[,length(unique(factor.legend))])) +			
			theme(legend.key.size=unit(13,'mm')) + guides(colour = guide_legend(override.aes = list(size=5))) + 						
			geom_bar(stat="identity") + facet_grid(stat ~ ., scales='free', margins=FALSE)
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "m2wmx.tp_sasky_3da_PYIWtotalmargin.pdf", sep='')
	ggsave(file=file, w=8,h=8)
	#	PYIW Margin in % before adjustment	
	ggplot(subset(run.tp, !is.na(n)), aes(x=t.period, y=prop, fill=factor.legend)) + labs(x="calendar time periods", y="% PYIW") +
			scale_fill_manual(name='Transmission times', values = colorRampPalette(brewer.pal(9, "Set1"),interpolate='spline',space="Lab")(run.tp[,length(unique(factor.legend))])) +			
			theme(legend.key.size=unit(13,'mm')) + guides(colour = guide_legend(override.aes = list(size=5))) + 						
			geom_bar(stat="identity") + facet_grid(stat ~ ., scales='free', margins=FALSE)
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "m2wmx.tp_sasky_3da_PYIWpropmargin.pdf", sep='')
	ggsave(file=file, w=8,h=8)
	#
	#	sequence selection bias
	#
	file		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_Wed_Dec_18_11:37:00_2013_Yscore3da_tables_m2Bwmx.R'
	load(file)
	run.bias	<- ans$risk.table
	run.bias[, stat.orig:=stat]
	set(run.bias, run.bias[,which(stat=='X.clu')], 'stat', 'clustering')
	set(run.bias, run.bias[,which(stat=='X.seq')], 'stat', 'sequenced')	
	set(run.bias, run.bias[,which(stat=='X.msm')], 'stat', 'cohort')
	set(run.bias, run.bias[,which(stat=='YX')], 'stat', 'potential transmitter')
	set(run.bias, NULL, 'stat', run.bias[, factor(stat, levels=c('cohort','sequenced','clustering','potential transmitter'))])
	tmp		<- c(	'Undiagnosed,\nEvidence for acute infection\nat diagnosis',
			'Undiagnosed,\nEvidence for recent infection\nat diagnosis',
			'Undiagnosed,\nNo recent infection',
			'Diagnosed < 3mo,\nEvidence for acute infection\nat diagnosis',
			'Diagnosed < 3mo,\nEvidence for recent infection\nat diagnosis',
			'Diagnosed,\nlowest CD4 to date\n> 500',
			'Diagnosed,\nlowest CD4 to date\n> 350',
			'Diagnosed,\nlowest CD4 to date\n<= 350',
			'Diagnosed,\nmissing 1st CD4',
			'ART initiated,\nVL missing',
			'ART initiated,\nVL continually\nsuppressed', 
			'ART initiated,\nVL not continually\nsuppressed')
	tmp		<- data.table( factor.legend= factor(tmp), factor=c("UAy","UAm","U","DAy","DAm","Dtg500","Dtl500","Dtl350","Dt.NA","ART.vlNA","ART.suA.Y","ART.suA.N"))
	run.bias<- merge(run.bias, tmp, by='factor')	
	#	barplot of proportions, x=stat, facet by factor
	ggplot(subset(run.bias, stat.orig!='YX' & stat.orig!='X.clu'), aes(x=stat, y=p, fill=factor.legend))	+ labs(x="transmission times", y="% PYIW") +
			scale_fill_manual(name='Transmission times', values = colorRampPalette(brewer.pal(9, "Set1"),interpolate='spline',space="Lab")(run.bias[,length(unique(factor.legend))])) +
			theme(strip.text.x = element_text(size = 8), legend.key.size=unit(13,'mm')) +
			geom_bar(stat="identity") + facet_grid(. ~ factor.legend, margins=FALSE)
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "m2wmx_sasky_3da_BIAS_sequence.pdf", sep='')
	ggsave(file=file, w=24,h=8)
	#	barplot of proportions, x=stat, facet by factor
	ggplot(subset(run.bias, stat.orig!='YX' & stat.orig!='X.msm'), aes(x=stat, y=p, fill=factor.legend))	+ labs(x="transmission times", y="% PYIW") +
			scale_fill_manual(name='Transmission times', values = colorRampPalette(brewer.pal(9, "Set1"),interpolate='spline',space="Lab")(run.bias[,length(unique(factor.legend))])) +
			theme(strip.text.x = element_text(size = 8), legend.key.size=unit(13,'mm')) +
			geom_bar(stat="identity") + facet_grid(. ~ factor.legend, margins=FALSE)
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "m2wmx_sasky_3da_BIAS_clustering.pdf", sep='')
	ggsave(file=file, w=24,h=8)
	

	#	check ART risk factor model type/number drugs" ART NRT.NNRTI vs NRT.PI
	subset(runs.risk, method.risk%in%c("m3.tnic","m3.tni",'m3.tnicvNo') & stat=='RR'	& factor.ref=='ART.3.NRT.PI' & factor%in%c('ART.3.NRT.NNRT'))
	
	#	compare fit of ART risk factor models 
	tmp	<- subset(runs.fit, 	method.risk%in%c(	"m3.i", "m3.ni",  "m3.nic",  "m3.tni", "m3.tnic", "m3.tniv", "m3.tnicv", "m3.tnicvNo",        
													"m3.i.clu", "m3.ni.clu", "m3.nic.clu", "m3.tni.clu", "m3.tnic.clu", "m3.tniv.clu", "m3.tnicvNo.clu",
													'm3.nic.adj','m3.nicv.adj','m3.tnic.adj','m3.tnicv.adj','m3.tnicvNo.adj' ))	
	tmp[, dAIC:= AIC-min(AIC)]
	setkey(tmp, dAIC, method.brl, method.dating, method.risk)
	tmp
	#	best sasky 3da		m3.nic.adj m3.nicv.adj m3.tnic.adj m3.tnicv.adj
	tmp	<- subset(runs.risk, grepl('m3.tnic',method.risk,fixed=1) & !grepl('clu',method.risk,fixed=1) & stat=='RR' & (factor=='ART.I' | risk=='ART.I'))
	tmp[, list(v=mean(v), l95.bs=mean(l95.bs), u95.bs=mean(u95.bs))]
	#	ART.I always significant among 		"m3.tnic"        "m3.tnicvNo"     "m3.tnic.adj"    "m3.tnicv.adj"   "m3.tnicvNo.adj"
	#	1.564469 1.273714 1.920072
	tmp	<- subset(runs.risk, grepl('m3.tnic',method.risk,fixed=1) & !grepl('clu',method.risk,fixed=1) & stat=='RR' & (factor=='ART.pulse.Y' | risk=='ART.pulse'))
	tmp[, list(v=mean(v), l95.bs=mean(l95.bs), u95.bs=mean(u95.bs))]
	#	ART.pulse always significant among 		"m3.tnic"        "m3.tnicvNo"     "m3.tnic.adj"    "m3.tnicv.adj"   "m3.tnicvNo.adj"
	#	1.564469 1.273714 1.920072
	tmp	<- subset(runs.risk, grepl('m3.tnic',method.risk,fixed=1) & !grepl('clu',method.risk,fixed=1) & stat=='RR' & (factor=='ART.F' | risk=='ART.F'))
	tmp[, list(v=mean(v), l95.bs=mean(l95.bs), u95.bs=mean(u95.bs))]
	#	ART.F NOT significant for	"m3.tnicvNo"   "m3.tnicvNo.adj"
	#	1.34313 1.066642 1.673151
	#	check ART risk factor model type/number drugs"
	tmp	<- subset(runs.risk, grepl('m3.tnic',method.risk,fixed=1) & !grepl('clu',method.risk,fixed=1) & stat=='RR' & (factor=='ART.P' | risk=='ART.P'))
	tmp[, list(v=mean(v), l95.bs=mean(l95.bs), u95.bs=mean(u95.bs))]
	#	ART.P NEVER significant	
	#	1.182005 0.8575869 1.581313
	tmp	<- subset(runs.risk, grepl('m3.tnic',method.risk,fixed=1) & !grepl('clu',method.risk,fixed=1) & stat=='RR' & (factor=='ART.A' | risk=='ART.A'))
	tmp[, list(v=mean(v), l95.bs=mean(l95.bs), u95.bs=mean(u95.bs))]
	#	ART.A typically NOT significant
	#	1.466913 0.9311993 1.757025
	
	#
	#	Table 3		2011
	#
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnic') & stat=='RR' & method.recentctime=='2011', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnicv.clu.censp') & method.recentctime=='2011' & stat=='RR.term', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )	
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnic.censp') & stat=='RR' & method.recentctime=='2011', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnicNo.censp') & stat=='RR' & method.recentctime=='2011', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnic.censp') & stat=='P' & method.recentctime=='2011', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnic.censp') & stat=='RI' & method.recentctime=='2011', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnicNo.censp') & stat=='RI' & method.recentctime=='2011', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnicNoMV') & stat=='RR.term' & method.recentctime=='2011', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnicNoMV') & stat=='P' & method.recentctime=='2011', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	tmp	<- subset(runs.table,  method.risk%in%c('m3.tnicNoMV') & method.recentctime=='2011' )
	
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnicMV.adj') & method.recentctime=='2011' & stat=='RR.term', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnicMV.adj') & method.recentctime=='2011' & stat=='P', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnicMV.adj') & method.recentctime=='2011' & stat=='RI', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnic.censp') & stat=='RR' & method.recentctime=='2011', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	
	subset(runs.risk,  method.risk%in%c('m3.indmx') & stat=='RR.raw' & method.brl=='3da' & !grepl('wstar',method.risk) & method.recentctime=='2011', c(factor.ref, factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	#factor.ref      factor         v     l95.bs    u95.bs method.risk method.dating method.nodectime method.brl
	#1:     ART.OK       ART.A 1.5809488 1.28377351 5.9106179    m3.indmx         sasky              any        3da
	#2:     ART.OK       ART.F 1.5807045 1.02548412 2.4239759    m3.indmx         sasky              any        3da
	#3:     ART.OK       ART.I 3.1522717 1.91880233 4.7885494    m3.indmx         sasky              any        3da
	#4:     ART.OK       ART.P 0.3567823 0.01803316 0.9181399    m3.indmx         sasky              any        3da
	#5:     ART.OK       ART.T 0.9163154 0.54822956 1.4262255    m3.indmx         sasky              any        3da
	#6:     ART.OK ART.pulse.Y 1.7401583 0.73739314 3.0250908    m3.indmx         sasky              any        3da
	subset(runs.risk,  method.risk%in%c('m3.indmx') & stat=='RR' & method.brl=='3da' & !grepl('wstar',method.risk) & method.recentctime=='2011', c(factor.ref, factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	#factor.ref      factor         v    l95.bs   u95.bs method.risk method.dating method.nodectime method.brl
	#1:     ART.OK       ART.A 1.1416174 0.9259452 4.262792    m3.indmx         sasky              any        3da
	#2:     ART.OK       ART.F 1.2895961 0.8692488 1.941736    m3.indmx         sasky              any        3da
	#3:     ART.OK       ART.I 2.7976588 1.7092032 4.248215    m3.indmx         sasky              any        3da
	#4:     ART.OK       ART.P 0.4625983 0.0692564 1.126041    m3.indmx         sasky              any        3da
	#5:     ART.OK       ART.T 0.8903901 0.5585588 1.366896    m3.indmx         sasky              any        3da
	#6:     ART.OK ART.pulse.Y 1.5200023 0.7179900 2.571243    m3.indmx         sasky              any        3da	
	subset(runs.risk,  method.risk%in%c('m3.indmxNo') & stat=='RR.raw' & method.brl=='3da' & !grepl('wstar',method.risk) & method.recentctime=='2011', c(factor.ref, factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	#factor.ref      factor         v     l95.bs    u95.bs method.risk method.dating method.nodectime method.brl
	#1:     ART.OK       ART.A 1.3430146 1.02334495 5.0424648  m3.indmxNo         sasky              any        3da
	#2:     ART.OK       ART.F 1.3428070 0.83103794 2.1739126  m3.indmxNo         sasky              any        3da
	#3:     ART.OK       ART.I 2.6778519 1.66006690 4.4094826  m3.indmxNo         sasky              any        3da
	#4:     ART.OK       ART.P 0.3030862 0.02019137 0.7817204  m3.indmxNo         sasky              any        3da
	#5:     ART.OK       ART.T 0.7784091 0.45062607 1.2936543  m3.indmxNo         sasky              any        3da
	#6:     ART.OK ART.pulse.Y 1.4782629 0.63785011 2.8238263  m3.indmxNo         sasky              any        3da
	subset(runs.risk,  method.risk%in%c('m3.indmxMV') & stat=='RR' & method.brl=='3da' & !grepl('wstar',method.risk) & method.recentctime=='2011', c(factor.ref, factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	#factor.ref      factor         v    l95.bs   u95.bs method.risk method.dating method.nodectime method.brl
	#1:     ART.OK       ART.A 1.1554814 0.9359808 4.498911  m3.indmxMV         sasky              any        3da
	#2:     ART.OK       ART.F 1.3049506 0.8697021 1.934548  m3.indmxMV         sasky              any        3da
	#3:     ART.OK       ART.I 2.8301160 1.7652080 4.324338  m3.indmxMV         sasky              any        3da
	#4:     ART.OK       ART.P 0.4453794 0.0922852 1.103289  m3.indmxMV         sasky              any        3da
	#5:     ART.OK       ART.T 0.9008525 0.5341716 1.416596  m3.indmxMV         sasky              any        3da
	#6:     ART.OK ART.pulse.Y 1.3855881 0.6196826 2.602033  m3.indmxMV         sasky              any        3da

	subset(runs.risk,  method.risk%in%c('m3.n3mx') & stat=='RR' & method.brl=='3da' & !grepl('wstar',method.risk) & method.recentctime=='2011', c(factor.ref, factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	#factor.ref      factor         v     l95.bs   u95.bs method.risk method.dating method.nodectime method.brl
	#1: ART.3.2NRT.X   ART.3.OTH 0.4500962 0.07505392 1.143388     m3.n3mx         sasky              any        3da
	#2: ART.3.2NRT.X       ART.A 1.1262118 0.91330526 4.243472     m3.n3mx         sasky              any        3da
	#3: ART.3.2NRT.X       ART.F 1.2687952 0.84479613 1.876250     m3.n3mx         sasky              any        3da
	#4: ART.3.2NRT.X       ART.I 2.7654109 1.75653222 4.225697     m3.n3mx         sasky              any        3da
	#5: ART.3.2NRT.X       ART.P 0.4530566 0.08225300 1.080912     m3.n3mx         sasky              any        3da
	#6: ART.3.2NRT.X       ART.T 0.8765734 0.53770979 1.398853     m3.n3mx         sasky              any        3da
	#7: ART.3.2NRT.X      ART.g3 0.9308021 0.21832473 1.935641     m3.n3mx         sasky              any        3da
	#8: ART.3.2NRT.X      ART.l3 2.9636815 1.38470273 5.205899     m3.n3mx         sasky              any        3da
	#9: ART.3.2NRT.X ART.pulse.Y 1.4981287 0.69358328 2.537469     m3.n3mx         sasky              any        3da
	subset(runs.risk,  method.risk%in%c('m3.n3mx') & stat=='RR.raw' & method.brl=='3da' & !grepl('wstar',method.risk) & method.recentctime=='2011', c(factor.ref, factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	#factor.ref      factor         v     l95.bs    u95.bs method.risk method.dating method.nodectime method.brl
	#1: ART.3.2NRT.X   ART.3.OTH 0.4438304 0.04521115 1.1082716     m3.n3mx         sasky              any        3da
	#2: ART.3.2NRT.X       ART.A 1.5562210 1.24193228 6.0115599     m3.n3mx         sasky              any        3da
	#3: ART.3.2NRT.X       ART.F 1.5559805 0.96662769 2.3868770     m3.n3mx         sasky              any        3da
	#4: ART.3.2NRT.X       ART.I 3.1029667 1.89455034 4.8118123     m3.n3mx         sasky              any        3da
	#5: ART.3.2NRT.X       ART.P 0.3512018 0.02539939 0.8589455     m3.n3mx         sasky              any        3da
	#6: ART.3.2NRT.X       ART.T 0.9019832 0.53588384 1.4316282     m3.n3mx         sasky              any        3da
	#7: ART.3.2NRT.X      ART.g3 0.8665310 0.13560173 1.9205182     m3.n3mx         sasky              any        3da
	#8: ART.3.2NRT.X      ART.l3 3.3574769 1.57825530 6.0175172     m3.n3mx         sasky              any        3da
	#9: ART.3.2NRT.X ART.pulse.Y 1.7129403 0.74646530 3.0757344     m3.n3mx         sasky              any        3da	
	subset(runs.risk,  method.risk%in%c('m3.n3mxMV') & stat=='RR' & method.brl=='3da' & !grepl('wstar',method.risk) & method.recentctime=='2011', c(factor.ref, factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	#factor.ref      factor         v     l95.bs   u95.bs method.risk method.dating method.nodectime method.brl
	#1: ART.3.2NRT.X   ART.3.OTH 0.4413032 0.08216935 1.121655   m3.n3mxMV         sasky              any        3da
	#2: ART.3.2NRT.X       ART.A 1.1454337 0.95690890 4.360743   m3.n3mxMV         sasky              any        3da
	#3: ART.3.2NRT.X       ART.F 1.2889582 0.85464553 1.943280   m3.n3mxMV         sasky              any        3da
	#4: ART.3.2NRT.X       ART.I 2.8008671 1.81053747 4.303314   m3.n3mxMV         sasky              any        3da
	#5: ART.3.2NRT.X       ART.P 0.4381825 0.07523197 1.151721   m3.n3mxMV         sasky              any        3da
	#6: ART.3.2NRT.X       ART.T 0.8913099 0.54372083 1.394639   m3.n3mxMV         sasky              any        3da
	#7: ART.3.2NRT.X      ART.g3 0.9663522 0.21397830 2.062098   m3.n3mxMV         sasky              any        3da
	#8: ART.3.2NRT.X      ART.l3 3.0231784 1.40781727 5.397656   m3.n3mxMV         sasky              any        3da
	#9: ART.3.2NRT.X ART.pulse.Y 1.3443834 0.66303983 2.496730   m3.n3mxMV         sasky              any        3da
	subset(runs.risk,  method.risk%in%c('m3.n3mx') & stat=='P.rawbias.e0' & method.brl=='3da' & !grepl('wstar',method.risk) & method.recentctime=='2011', c(factor.ref, factor, v, l95.bs, u95.bs,method.risk ) )
	#factor.ref       factor           v       l95.bs     u95.bs method.risk
	#1:       None ART.3.2NRT.X 0.358805848 0.2831631504 0.44485251     m3.n3mx
	#2:       None    ART.3.OTH 0.017591054 0.0017525908 0.04156421     m3.n3mx
	#3:       None        ART.A 0.004515470 0.0000000000 0.01462649     m3.n3mx
	#4:       None        ART.F 0.157042362 0.1118320361 0.20844342     m3.n3mx
	#5:       None        ART.I 0.219970392 0.1576249623 0.29026857     m3.n3mx
	#6:       None        ART.P 0.006861843 0.0004894542 0.01616562     m3.n3mx
	#7:       None        ART.T 0.145749218 0.0977239042 0.20333479     m3.n3mx
	#8:       None       ART.g3 0.038668205 0.0066903371 0.07920521     m3.n3mx
	#9:       None       ART.l3 0.029665692 0.0140151459 0.05086179     m3.n3mx
	#10:       None  ART.pulse.Y 0.021129916 0.0092062212 0.03526553     m3.n3mx
	subset(runs.risk,  method.risk%in%c('m3.n3mx') & stat=='RI.rawbias.e0' & method.brl=='3da' & !grepl('wstar',method.risk) & method.recentctime=='2011', c(factor.ref, factor, v, l95.bs, u95.bs,method.risk ) )
	#factor.ref       factor         v     l95.bs    u95.bs method.risk
	#1:       None ART.3.2NRT.X 0.8200816 0.64719374 1.0167487     m3.n3mx
	#2:       None    ART.3.OTH 0.3639772 0.03626292 0.8600066     m3.n3mx
	#3:       None        ART.A 1.2762283 1.15268180 4.5701982     m3.n3mx
	#4:       None        ART.F 1.2760310 0.90867933 1.6936848     m3.n3mx
	#5:       None        ART.I 2.5446860 1.82345463 3.3579172     m3.n3mx
	#6:       None        ART.P 0.2880141 0.02054400 0.6785242     m3.n3mx
	#7:       None        ART.T 0.7396998 0.49596394 1.0319555     m3.n3mx
	#8:       None       ART.g3 0.7106262 0.12295188 1.4555963     m3.n3mx
	#9:       None       ART.l3 2.7534051 1.30080815 4.7207090     m3.n3mx
	#10:       None  ART.pulse.Y 1.4047509 0.61204444 2.3445091     m3.n3mx
	subset(runs.risk,  method.risk%in%c('m3.n3mx') & stat=='N.raw' & method.brl=='3da' & !grepl('wstar',method.risk) & method.recentctime=='2011', c(factor.ref, factor, v, l95.bs, u95.bs,method.risk ) )
	#factor.ref       factor          v      l95.bs    u95.bs method.risk
	#1:       None ART.3.2NRT.X 12.2191253  9.24489421 15.611420     m3.n3mx
	#2:       None    ART.3.OTH  0.5794369  0.06367325  1.440342     m3.n3mx
	#3:       None        ART.A  0.4794567  0.00000000  1.438370     m3.n3mx
	#4:       None        ART.F 12.9526215  8.91916901 17.448664     m3.n3mx
	#5:       None        ART.I 17.5537989 11.96303877 23.153660     m3.n3mx
	#6:       None        ART.P  0.4214211  0.02915683  1.001240     m3.n3mx
	#7:       None        ART.T  6.1430989  4.00832333  8.688211     m3.n3mx
	#8:       None       ART.g3  2.4332843  0.41869599  5.234366     m3.n3mx
	#9:       None       ART.l3  1.1961089  0.58104452  2.034773     m3.n3mx
	#10:       None  ART.pulse.Y  1.7324934  0.77997668  2.900998     m3.n3mx	

	load('/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_140603/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3ga_denomCLU_model3_m3.btnicNoMV.clu.R')
	subset(ans$risk, stat=='P.e0cp')
	
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.btnicMV') & method.recentctime=='2011' & method.denom=='CLU' & method.dating=='sasky', c(stat, factor, v, l95.bs, u95.bs, method.risk, method.brl ) )
	set(tmp, NULL, 'v', tmp[, round(v, d=3)])
	set(tmp, NULL, 'l95.bs', tmp[, round(l95.bs, d=3)])
	set(tmp, NULL, 'u95.bs', tmp[, round(u95.bs, d=3)])	
	subset(tmp, stat=='P.raw.e0')
	subset(tmp, stat=='RI.raw.e0' & method.brl=='3da')
	
	tmp	<- subset(runs.risk,  stat=='RR.term.ptx' & method.risk%in%c('m3.atnicMV') & method.recentctime=='2011' & method.denom=='CLU' & method.dating=='sasky', c(stat, factor, factor.ref, v, l95.bs, u95.bs, method.risk, method.brl ) )
	subset(runs.risk,  stat=='RR.term' & method.risk%in%c('m3.atnicMV') & method.recentctime=='2011' & method.denom=='CLU' & method.dating=='sasky', c(stat, factor, factor.ref, v, l95.bs, u95.bs, method.risk, method.brl ) )
	
	risk.table	<- subset(runs.table,  method.risk%in%c('m3.tnicNoMV') & method.recentctime=='2011' )
	subset( risk.table, grepl('ART',factor) )[, sum(p), by='stat']
	#
	#	completeness of ATHENA data
	#
	if(1)		
		file	<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_Wed_Dec_18_11:37:00_2013_RIPDT_RIMSM_3da_info.R"
	else
	{
		file	<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_2011_sasky_Wed_Dec_18_11:37:00_2013_RIPDT_RIMSM_3da_info.R"	
		set(ri.allmsm.su$rin, ri.allmsm.su$rin[,which(t.period.max==2013.125)],'t.period.max',2011.)		
		set(ri.allmsm.su$ris, ri.allmsm.su$ris[,which(t.period.max==2013.125)],'t.period.max',2011.)
		set(ri.PT.su$rin, ri.PT.su$rin[,which(t.period.max==2013.125)],'t.period.max',2011.)		
		set(ri.PT.su$ris, ri.PT.su$ris[,which(t.period.max==2013.125)],'t.period.max',2011.)
		
	}	
	tmp		<- load(file)
	#layer: col1 number MSM in care		col2 number CD4 & VL per person		col3	%at least one per year
	ggplot(subset(ri.allmsm.resolution,t<2011), aes(x=t, y=nPatient)) + labs(x="year", y="HIV infected MSM in cohort") + geom_bar(stat="identity", fill="#E41A1C") 
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "RIMSM_totalMSM.pdf", sep='')
	ggsave(file=file, w=3,h=8)	
	#	frequency
	tmp		<- copy(ri.allmsm.resolution)
	setnames(tmp, c('mean.nRNA.pP','mean.nCD4.pP'), c('viral load','CD4'))
	ggplot( melt(subset(tmp,t<2011), measure=c('viral load','CD4'), variable.name='measured'), aes(x=t, y=value, colour=measured)) + labs(x="year", y='average number of samples per person per year') + scale_colour_brewer(palette='Set1') + 
		geom_point() + geom_line()
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "RIMSM_resolution.pdf", sep='')
	ggsave(file=file, w=4,h=8)
	#	completeness
	tmp		<- copy(ri.allmsm.resolution)
	setnames(tmp, c('one.RNA.pP','one.CD4.pP'), c('viral load','CD4'))
	ggplot( melt(subset(tmp,t<2011), measure=c('viral load','CD4'), variable.name='measured'), aes(x=t, y=value, colour=measured)) + labs(x="year", y='proportion of individuals with at least one sample per year') + scale_colour_brewer(palette='Set1') +			
			geom_point() + geom_line()
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "RIMSM_completeness.pdf", sep='')
	ggsave(file=file, w=4,h=8)
	#
	#	composition recently infected
	#
	
	#	could do paired colours and then have lines + ribbon	x=time period	pch=Study-RI vs Cohort-RI	colour=factor	
	#	col1	y=total number	
	tmp		<- subset( ri.allmsm.su$rin, select=c(t.period, t.period.min, t.period.max, Patient) )	
	tmp[, group:='Cohort-RI']	
	tmp2	<- subset( ri.PT.su$rin, select=c(t.period, t.period.min, t.period.max, Patient) )	
	tmp2[, group:='Study-RI']
	tmp		<- rbind(tmp, tmp2)
	tmp[, t.period.long:= paste(round(t.period.min,d=1), '-\n', round(t.period.max,d=1),sep='')]
	ggplot(tmp, aes(x=t.period.long, y=Patient, pch=group, linetype=group, group=group)) + scale_y_continuous(limits=c(0,1100)) + labs(x="time period", y="recently infected MSM") + geom_point() + geom_line()
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "RIPDT_RIMSM_Total.pdf", sep='')
	ggsave(file=file, w=5,h=8)
	
	#	col2	y=%Evidence Yes/Potential
	tmp		<- merge( subset( ri.allmsm.su$rip, select=c(t.period, isAcute.Yes, isAcute.Maybe, stat) ), subset( ri.allmsm.su$rin, select=c(t.period, t.period.min, t.period.max) ), by='t.period' )
	tmp[, group:='Cohort-RI']	
	tmp2	<- merge( subset( ri.PT.su$rip, select=c(t.period, isAcute.Yes, isAcute.Maybe, stat) ), subset( ri.allmsm.su$rin, select=c(t.period, t.period.min, t.period.max) ), by='t.period' )
	tmp2[, group:='Study-RI']
	tmp		<- rbind(tmp, tmp2)
	tmp[, t.period.long:= paste(round(t.period.min,d=1), '-\n', round(t.period.max,d=1),sep='')]	
	setnames(tmp, c('isAcute.Yes','isAcute.Maybe'), c('Yes','Potentially'))
	tmp		<- melt( tmp, measure=c('Yes','Potentially') )
	tmp		<- cbind( subset(tmp, stat=='p', select=c(group, t.period.long, variable, value)), subset(tmp, stat=='p.l95', select=value), subset(tmp, stat=='p.u95', select=value) )
	setnames(tmp, c(4,5,6), c('mean','l95','u95'))
	tmp[, dummy:=paste(group,variable,sep='')]	
	ggplot(tmp, aes(x=t.period.long, y=mean, pch=group, linetype=group, group=dummy, colour=variable, fill=variable)) + scale_y_continuous(limits=c(0,1)) +
				labs(x="time period", y="Proportion") + 
				scale_colour_brewer(name='Evidence for\nacute infection', palette='Set1') + scale_fill_brewer(name='Evidence for\nacute infection', palette='Set1') +
				geom_ribbon(aes(ymin=l95, ymax=u95, linetype=NA), alpha=0.3) + facet_grid(. ~ variable, margins=FALSE) + geom_point() + geom_line()
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "RIPDT_RIMSM_Acute.pdf", sep='')
	ggsave(file=file, w=7,h=8)
		
	#	col3	y=%Region
	tmp		<- merge( subset( ri.allmsm.su$rip, select=c(t.period, Amst, N, E, S, W, stat) ), subset( ri.allmsm.su$rin, select=c(t.period, t.period.min, t.period.max) ), by='t.period' )
	tmp[, group:='Cohort-RI']
	tmp2	<- merge( subset( ri.PT.su$rip, select=c(t.period, Amst, N, E, S, W, stat) ), subset( ri.allmsm.su$rin, select=c(t.period, t.period.min, t.period.max) ), by='t.period' )	
	tmp2[, group:='Study-RI']
	tmp		<- rbind(tmp, tmp2)
	tmp[, t.period.long:= paste(round(t.period.min,d=1), '-\n', round(t.period.max,d=1),sep='')]	
	tmp		<- melt( tmp, measure=c('Amst','N','E','S','W') )	
	tmp		<- cbind( subset(tmp, stat=='p', select=c(group, t.period.long, variable, value)), subset(tmp, stat=='p.l95', select=value), subset(tmp, stat=='p.u95', select=value) )
	setnames(tmp, c(4,5,6), c('mean','l95','u95'))
	tmp[, dummy:=paste(group,variable,sep='')]	
	ggplot(tmp, aes(x=t.period.long, y=mean, pch=group, linetype=group, group=dummy, colour=variable, fill=variable)) + scale_y_continuous(limits=c(0,1)) +
			labs(x="time period", y="Proportion") + 
			scale_colour_brewer(name='Region in care', palette='Set1') + scale_fill_brewer(name='Region in care', palette='Set1') +
			geom_ribbon(aes(ymin=l95, ymax=u95, linetype=NA), alpha=0.3) + facet_grid(. ~ variable, margins=FALSE) + geom_point() + geom_line()
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "RIPDT_RIMSM_Region.pdf", sep='')
	ggsave(file=file, w=15,h=8)	
	#	col4	y=Age
	tmp		<- subset( ri.allmsm.su$ris, select=c(t.period, t.period.min, t.period.max, AnyPos_A, stat) )
	tmp[, group:='Cohort-RI']
	tmp2	<- subset( ri.PT.su$ris, select=c(t.period, t.period.min, t.period.max, AnyPos_A, stat) )
	tmp2[, group:='Study-RI']
	tmp		<- rbind(tmp, tmp2)
	tmp[, t.period.long:= paste(round(t.period.min,d=1), '-\n', round(t.period.max,d=1),sep='')]	
	tmp		<- cbind( subset(tmp, stat=='quantile_0.5', select=c(group, t.period.long, AnyPos_A)), subset(tmp, stat=='quantile_0.05', select=AnyPos_A), subset(tmp, stat=='quantile_0.95', select=AnyPos_A) )
	setnames(tmp, c(3,4,5), c('mean','l95','u95'))
	ggplot(tmp, aes(x=t.period.long, y=mean, pch=group, linetype=group, group=group, fill=group, colour=group)) + 
			labs(x="time period", y="Age at diagnosis") +
			scale_colour_brewer(name='group', palette='Set1') + scale_fill_brewer(name='group', palette='Set1') +
			geom_ribbon(aes(ymin=l95, ymax=u95, linetype=NA), alpha=0.3) + geom_point() + geom_line() 
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "RIPDT_RIMSM_AnyPos_A.pdf", sep='')
	ggsave(file=file, w=5,h=8)
	#	col5	y=1st VL
	tmp		<- subset( ri.allmsm.su$ris, select=c(t.period, t.period.min, t.period.max, lRNA_T1, stat) )
	tmp[, group:='Cohort-RI']
	tmp2	<- subset( ri.PT.su$ris, select=c(t.period, t.period.min, t.period.max, lRNA_T1, stat) )
	tmp2[, group:='Study-RI']
	tmp		<- rbind(tmp, tmp2)
	tmp[, t.period.long:= paste(round(t.period.min,d=1), '-\n', round(t.period.max,d=1),sep='')]	
	tmp		<- cbind( subset(tmp, stat=='quantile_0.5', select=c(group, t.period.long, lRNA_T1)), subset(tmp, stat=='quantile_0.05', select=lRNA_T1), subset(tmp, stat=='quantile_0.95', select=lRNA_T1) )
	setnames(tmp, c(3,4,5), c('mean','l95','u95'))
	ggplot(tmp, aes(x=t.period.long, y=mean, pch=group, linetype=group, group=group, fill=group, colour=group)) + 
			labs(x="time period", y="First log10 viral load") +
			scale_colour_brewer(name='group', palette='Set1') + scale_fill_brewer(name='group', palette='Set1') +
			geom_ribbon(aes(ymin=l95, ymax=u95, linetype=NA), alpha=0.3) + geom_point() + geom_line() 
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "RIPDT_RIMSM_VL.pdf", sep='')
	ggsave(file=file, w=5,h=8)
	#	col6	y=1st CD4
	tmp		<- subset( ri.allmsm.su$ris, select=c(t.period, t.period.min, t.period.max, CD4_T1, stat) )
	tmp[, group:='Cohort-RI']
	tmp2	<- subset( ri.PT.su$ris, select=c(t.period, t.period.min, t.period.max, CD4_T1, stat) )
	tmp2[, group:='Study-RI']
	tmp		<- rbind(tmp, tmp2)
	tmp[, t.period.long:= paste(round(t.period.min,d=1), '-\n', round(t.period.max,d=1),sep='')]	
	tmp		<- cbind( subset(tmp, stat=='quantile_0.5', select=c(group, t.period.long, CD4_T1)), subset(tmp, stat=='quantile_0.05', select=CD4_T1), subset(tmp, stat=='quantile_0.95', select=CD4_T1) )
	setnames(tmp, c(3,4,5), c('mean','l95','u95'))
	ggplot(tmp, aes(x=t.period.long, y=mean, pch=group, linetype=group, group=group, fill=group, colour=group)) + 
			labs(x="time period", y="First CD4 count") +
			scale_colour_brewer(name='group', palette='Set1') + scale_fill_brewer(name='group', palette='Set1') +
			geom_ribbon(aes(ymin=l95, ymax=u95, linetype=NA), alpha=0.3) + geom_point() + geom_line() 
	file	<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', "RIPDT_RIMSM_CD4.pdf", sep='')
	ggsave(file=file, w=5,h=8)
	#
	#
	#
	subset(runs.risk, method.brl=='3da' & method.dating=='sasky' & grepl('m5.tAc.tp',method.risk) &  !grepl('wstar',method.risk) & stat=='P.rawbias.e0' , select=c(stat, factor.ref, factor, v, l95.bs, u95.bs))
	subset(runs.risk, method.brl=='3da' & method.dating=='sasky' & grepl('m5.tAc.tp',method.risk) &  !grepl('wstar',method.risk) & stat=='N.raw' , select=c(stat, factor.ref, factor, v, l95.bs, u95.bs))
	subset(runs.risk, method.brl=='3da' & method.dating=='sasky' & grepl('m5.tAc.tp',method.risk) &  !grepl('wstar',method.risk) & stat=='RI.rawbias.e0' , select=c(stat, factor.ref, factor, v, l95.bs, u95.bs))
	subset(runs.risk, method.brl=='3da' & method.dating=='sasky' & grepl('m5.tAc.tp',method.risk) &  !grepl('wstar',method.risk) & stat=='RR.raw' , select=c(stat, factor.ref, factor, v, l95.bs, u95.bs))	
	#
	#	MODEL 5 time trends
	#
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.412, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )

	ylab		<- "Proportion of transmissions"
	stat.select	<- c(	'P','P.e0','P.bias.e0','P.raw','P.raw.e0','P.rawbias.e0'	)
	#ylab		<- "Relative transmissibility"
	#stat.select	<- c(	'RI','RI.e0','RI.ptx','RI.ptx.e0','RI.raw','RI.raw.e0')
	method.clu	<- 'clu'; method.deno	<- 'CLU'	
	run.tp		<- subset(runs.risk, method.denom==method.deno & method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m5.tAc.tp',method.risk) & !grepl(method.clu,method.risk) & !grepl('wstar',method.risk) & !grepl('now',method.risk) & (grepl('P.',stat,fixed=1) | stat=='P') )
	#run.tp		<- subset(runs.risk, method.denom==method.deno & method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m5.tAc.tp',method.risk) & grepl(method.clu,method.risk) & grepl('wstar',method.risk) & !grepl('now',method.risk) & (grepl('P.',stat,fixed=1) | stat=='P') )
	#run.tp		<- subset(runs.risk, method.denom==method.deno & method.nodectime=='any' & method.brl=='3ga' & method.dating=='sasky' & grepl('m5.tA.tp',method.risk) & grepl(method.clu,method.risk) & !grepl('now',method.risk) & (grepl('P.',stat,fixed=1) | stat=='P') )
	#run.tp		<- subset(runs.risk, method.denom==method.deno & method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m5.tA.tp',method.risk) & grepl(method.clu,method.risk) & !grepl('now',method.risk) & (grepl('RI.',stat,fixed=1) | stat=='RI') )
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 2, nchar(factor)-2)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=100','45-',factor)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=20','<20',factor)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=25','20-24',factor)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=30','25-29',factor)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=35','30-34',factor)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=40','35-39',factor)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=45','40-44',factor)])
	run.tp[, group:=NA_character_]
	set(run.tp, run.tp[,which(grepl('3',factor) | grepl('4',factor))], 'group', '30 or older')
	set(run.tp, run.tp[,which(grepl('2',factor))], 'group', 'below 30')	
	set(run.tp, NULL, 'group', run.tp[, factor(group, levels=c('below 30','30 or older'))])
	run.tp		<- merge(run.tp, tperiod.info,by='t.period')
	run.tp[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]
	#
	dummy	<- lapply(seq_along(stat.select), function(i)
		{
			cat(paste('\nprocess', stat.select[i]))
			ggplot(subset(run.tp, !is.na(v) & stat==stat.select[i]), aes(x=factor, y=v, fill=factor, colour=factor)) + labs(x="", y=ylab) + 
					scale_y_continuous(breaks=seq(0,0.3,0.1), labels=paste(seq(0,0.3,0.1)*100,'%',sep='')) + scale_x_discrete(breaks=NULL) +
					scale_fill_brewer(palette='RdBu',name='from age group') + 
					scale_colour_manual(name='from age group', values = rep('black',7)) +					
					guides(colour=FALSE, fill = guide_legend(override.aes = list(size=5))) +
					theme(legend.key.size=unit(11,'mm'), plot.margin=unit(c(0,0,-10,0),"mm")) + #coord_flip() +
					geom_bar(stat='identity',binwidth=1, position='dodge')	+ geom_errorbar(aes(ymin=l95.bs, ymax=u95.bs), width=0.3, position=position_dodge(width=0.9))	+ 
					facet_grid(. ~ t.period.long, margins=FALSE)
			file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', stat.select[i],'_',subset(run.tp, !is.na(v) & stat==stat.select[i])[1, method.risk],"_prop_",method.clu,'_',subset(run.tp, !is.na(v) & stat==stat.select[i])[1, method.brl],'_denom',method.deno,'_', subset(run.tp, !is.na(v) & stat==stat.select[i])[1, method.recentctime],".pdf", sep='')
			cat(paste('\nsave to file',file))			
			ggsave(file=file, w=12,h=4)			
		})
	#
	#
	ylab		<- "Relative transmissibility"
	stat.select	<- c(	'RI','RI.e0','RI.bias.e0','RI.raw','RI.raw.e0','RI.rawbias.e0'	)
	method.clu	<- 'clu'; method.deno	<- 'CLU'	
	run.tp		<- subset(runs.risk, method.denom==method.deno & method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m5.tAc.tp',method.risk) & !grepl(method.clu,method.risk) & !grepl('wstar',method.risk) & !grepl('now',method.risk) & (grepl('RI.',stat,fixed=1) | stat=='RI') )	
	#run.tp		<- subset(runs.risk, method.denom==method.deno & method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m5.tAc.tp',method.risk) & grepl(method.clu,method.risk) & grepl('wstar',method.risk) & !grepl('now',method.risk) & (grepl('RI.',stat,fixed=1) | stat=='RI') )
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 2, nchar(factor)-2)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=100','45-',factor)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=20','<20',factor)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=25','20-24',factor)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=30','25-29',factor)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=35','30-34',factor)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=40','35-39',factor)])
	set(run.tp, NULL, 'factor', run.tp[, gsub('<=45','40-44',factor)])
	run.tp[, group:=NA_character_]
	set(run.tp, run.tp[,which(grepl('3',factor) | grepl('4',factor))], 'group', '30 or older')
	set(run.tp, run.tp[,which(grepl('2',factor))], 'group', 'below 30')	
	set(run.tp, NULL, 'group', run.tp[, factor(group, levels=c('below 30','30 or older'))])
	run.tp		<- merge(run.tp, tperiod.info,by='t.period')
	run.tp[, t.period.long:= paste(t.period.min, ' to\n ', t.period.max,sep='')]

	dummy	<- lapply(seq_along(stat.select), function(i)
			{
				cat(paste('\nprocess', stat.select[i]))
				ggplot(subset(run.tp, !is.na(v) & stat==stat.select[i]), aes(x=t.period.long, y=v, fill=factor, colour=factor)) + labs(x="", y=ylab) + 
						#scale_y_continuous(breaks=seq(0,0.3,0.1), labels=paste(seq(0,0.3,0.1)*100,'%',sep='')) + scale_x_discrete(breaks=NULL, limits=rev(c("<20","20-24","25-29","30-34","35-39","40-44","45-"))) +
						scale_fill_brewer(palette='RdBu',name='from age group') + 
						scale_colour_manual(name='from age group', values = rep('black',7)) +					
						guides(colour=FALSE, fill = guide_legend(override.aes = list(size=5))) +
						theme(legend.key.size=unit(11,'mm'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) + #coord_flip() +
						geom_bar(stat='identity',binwidth=1, position='dodge', show_guide=FALSE)	+ geom_errorbar(aes(ymin=l95.bs, ymax=u95.bs), width=0.3, position=position_dodge(width=0.9))	+ 
						facet_grid(. ~ factor, margins=FALSE)
				file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat), '_', stat.select[i],'_',subset(run.tp, !is.na(v) & stat==stat.select[i])[1, method.risk],"_prop_",method.clu,'_',subset(run.tp, !is.na(v) & stat==stat.select[i])[1, method.brl],'_denom',method.deno,'_', subset(run.tp, !is.na(v) & stat==stat.select[i])[1, method.recentctime],".pdf", sep='')
				cat(paste('\nsave to file',file))
				ggsave(file=file, w=12,h=5)				
			})
	set(run.tp, NULL, 'v', run.tp[, round(v, d=3)])
	set(run.tp, NULL, 'l95.bs', run.tp[, round(l95.bs, d=3)])
	set(run.tp, NULL, 'u95.bs', run.tp[, round(u95.bs, d=3)])
	subset(run.tp, stat=='RI.rawbias.e0', c(t.period.min, t.period.max, stat, method.brl, factor, v, l95.bs, u95.bs))
	subset(run.tp, stat=='RI.raw.e0', c(t.period.min, t.period.max, stat, method.brl, factor, v, l95.bs, u95.bs))
	
	
}
######################################################################################
project.athena.Fisheretal.RI.candidatetransmitters.table<- function(rit, df.all, tperiod.info, save.file=save.file, resume=TRUE, vl.suppressed=3, cd4.cut=c(-1,350,550,1e4), cd4.label=c('D1<=350','D1<=550','D1>550'), t.period=0.25)
{
	#vl.suppressed=3; cd4.cut=c(-1,350,550,1e4); cd4.label=c('D1<=350','D1<=550','D1>550')
	if(resume && !is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(!resume || is.na(save.file) || inherits(readAttempt, "try-error"))
	{								
		#	diagnosed and acute
		rit[, score.Inf:=NULL]
		rit[, U.score:=NULL]
		set(rit, rit[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
		set(rit, rit[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )		
		#	if transmitter acute, set undiagnosed to undiagnosed and acute
		set(rit, rit[, which(stage=='U' & t.isAcute=='Yes')], 'stage', 'UAy')
		set(rit, rit[, which(stage=='U' & t.isAcute=='Maybe')], 'stage', 'UAm')
		#	stratify Diag by first CD4 after diagnosis	
		tmp			<- subset(df.all, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
		setkey(tmp, Patient)
		tmp			<- unique(tmp)	
		tmp[, CD4.c:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
		cat(paste('\ntransmitters with no CD4 after diagnosis\n'))
		print( subset( tmp, is.na(CD4.c) ) )
		setnames(tmp, 'Patient', 't.Patient')
		set(tmp, tmp[, which(is.na(CD4.c))], 'CD4.c', 'D1.NA')
		rit			<- merge(rit, subset(tmp, select=c(t.Patient, CD4.c)), by='t.Patient')
		tmp			<- rit[, which(stage=='Diag')]
		set(rit, tmp, 'stage', rit[tmp, CD4.c])
		rit[, CD4.c:=NULL]
		rit[, table(stage)]
		
		rit[, stage.orig:= stage]
		# 	collect PoslRNA_TL and lRNA_TL etc, then get VL suppressed
		tmp		<- unique(subset( df.all, select=c(Patient, PoslRNA_T1) ))
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		rit		<- merge(rit, tmp, by= 't.Patient')	
		tmp		<- merge( unique(subset(df.all, select=Patient)), subset(df.viro, select=c(Patient, PosRNA, lRNA)), by='Patient' )
		set( tmp, NULL, 'PosRNA', hivc.db.Date2numeric(tmp[,PosRNA]) )	
		tmp		<- subset( tmp[, 	{
							tmp<- which(!is.na(lRNA))
							list(t.Patient=Patient[1], lRNA_TL=lRNA[tail(tmp,1)], PoslRNA_TL=PosRNA[tail(tmp,1)])
						}, by='Patient'], select=c(t.Patient, lRNA_TL, PoslRNA_TL) )
		rit		<- merge(rit, tmp, by= 't.Patient', all.x=TRUE)
		rit		<- merge(rit, rit[, {
							tmp	<- which(!is.na(lRNA))
							if(!length(tmp) & all(!is.na(contact)) & all(contact=='Yes') & !is.na(lRNA_TL[1]) & any(t>t.PoslRNA_T1) & lRNA_TL[1]<=vl.suppressed)
								lRNA.c	<- 'SuA.Y'
							else if(!length(tmp) & all(!is.na(contact)) & all(contact=='Yes') & !is.na(lRNA_TL[1]) & any(t>t.PoslRNA_T1) & lRNA_TL[1]>vl.suppressed)
								lRNA.c	<- 'SuA.N'
							else if(!length(tmp))
								lRNA.c	<- 'SuA.NA'
							list( lRNA.c= ifelse(length(tmp), ifelse( max(lRNA[tmp])>vl.suppressed, 'SuA.N', 'SuA.Y'), lRNA.c) )
						}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
		set(rit, rit[,which(stage=='ART.started' & lRNA.c=='SuA.Y')],'stage', 'ART.suA.Y' )
		set(rit, rit[,which(stage=='ART.started' & lRNA.c=='SuA.N')],'stage', 'ART.suA.N' )
		set(rit, rit[,which(stage=='ART.started')],'stage', 'ART.vlNA' )
		set(rit, NULL, 'stage', rit[, factor(as.character(stage))])
		#	compute totals per stage
		tmp				<- rit[ , table( t.period, stage, useNA='ifany') ]	
		ritn			<- cbind( data.table(t.period=factor(rownames(tmp))), as.data.table(unclass(tmp)) )
		ritn			<- merge(ritn, rit[, list(t.Patient=length(Patient)), by='t.period'], by='t.period')
		#	get proportions and confidence intervals	
		require(Hmisc)
		tmp				<- ritn[ , 	{
					tmp				<- lapply(.SD, function(x) as.double(binconf( x, t.Patient, alpha=0.05, method= "wilson", include.x=FALSE, include.n=FALSE, return.df=FALSE)) )								
					tmp[['stat']]	<- c('p','p.l95','p.u95')
					tmp
				}, by='t.period']
		ritn			<- merge(ritn, data.table(t.period=tperiod.info[,t.period],stat='n'), by='t.period')
		ritn			<- rbind(ritn, tmp)
		#	get person years during infection window
		tmp		<- t.period
		tmp		<- subset(ritn, stat=='n')[, lapply(.SD, '*', tmp ), by='t.period', .SDcol=which(!colnames(ritn)%in%c('stat','t.period'))]	
		ritn	<- rbind(ritn, merge(tmp, data.table(t.period=tperiod.info[,t.period],stat='pyiw'), by='t.period'))
		ritn	<- merge(tperiod.info, ritn, by='t.period')
		setkey(ritn, stat)
		if(!is.na(save.file))
		{
			cat(paste('\nsave ritn to',save.file))
			save(ritn, file=save.file)
		}
		rit		<- NULL
		gc()
	}
	ritn
}	
######################################################################################
project.athena.Fisheretal.CT.link2care<- function(info, link2care.cut=c(-Inf,15/365,30/365,3/12,Inf), link2care.label=c('<=15d','<=30d','<=3m','>3m'), t.startctime=1994)
{
	setkey(info, Patient)
	info	<- unique(info)
	#	link to care
	info[, PosCD4_dt:= PosCD4_T1-AnyPos_T1]
	info[, PoslRNA_dt:= PoslRNA_T1-AnyPos_T1]	
	info	<- merge(info, info[, list( PoslRNAorCD4_dt=ifelse(all(is.na(c(PoslRNA_dt,PosCD4_dt))),NA_real_,min(PoslRNA_dt,PosCD4_dt,na.rm=TRUE)) ), by='Patient'], by='Patient')	
	info[, Link2Care_dt:=cut(info[, PoslRNAorCD4_dt], breaks=link2care.cut, labels=link2care.label, right=TRUE)]
	#	link to care per year
	tmp			<- subset(info, select=c(Patient, AnyPos_T1, Link2Care_dt))
	tmp[, t:=floor(tmp[,AnyPos_T1])]
	tmp[, list(nLink2Care_dt=tabulate(Link2Care_dt, nbin=nlevels(tmp[, Link2Care_dt])), Link2Care_dt=levels(tmp[, Link2Care_dt])), by='t' ]
	tmp			<- subset(tmp, t>=t.startctime)
	info.t		<- tmp[, 	{
				ans				<- as.list( tabulate(Link2Care_dt, nbin=nlevels(tmp[, Link2Care_dt])) )
				names(ans)		<- paste('Link2Care.',levels(tmp[, Link2Care_dt]),sep='')
				ans[['nDiag']]	<- sum(as.numeric(ans))
				ans					
			}, by='t' ]
	info.t
}	
######################################################################################
project.athena.Fisheretal.CT.resolution<- function(info, df.viro, df.immu, t.startctime=1994, t.endctime=2013.)
{
	setkey(info, Patient)
	info	<- unique(info)
	#	number of VL per year after year of diagnosis
	tmp		<- merge(subset(info, select=c(Patient,AnyPos_T1,DateDied)),subset(df.viro, select=c(Patient,PosRNA)),by='Patient',all.x=TRUE,allow.cartesian=TRUE)
	set(tmp, NULL, 'PosRNA', floor(hivc.db.Date2numeric(tmp[,PosRNA])))
	set(tmp, NULL, 'AnyPos_T1', floor(tmp[,AnyPos_T1]))
	set(tmp, NULL, 'DateDied', ceiling(tmp[,DateDied]))
	set(tmp, tmp[, which(is.na(DateDied))], 'DateDied', ceiling(t.endctime))
	tmp		<- tmp[, 	{
				if( all(is.na(PosRNA)) )
					ans	<- list(t=seq.int(AnyPos_T1[1], DateDied[1]-1), nRNA=rep(0,DateDied[1]-AnyPos_T1[1]))
				else
					ans	<- list(t=seq.int(AnyPos_T1[1], DateDied[1]-1), nRNA=c(rep(0, min(PosRNA)-AnyPos_T1[1]), tabulate(PosRNA-min(PosRNA)+1, nbins=DateDied[1]-min(PosRNA))))
				ans
			}, by='Patient']
	
	tmp		<- subset(tmp, t>=t.startctime)
	info	<- merge( info, tmp[, list(mean.nRNA.py=mean(nRNA), sum.nRNA=sum(nRNA)) ,by='Patient'], by='Patient', all.x=TRUE )	
	set(info, info[, which(is.na(mean.nRNA.py))], 'mean.nRNA.py', 0.)
	set(info, info[, which(is.na(sum.nRNA))], 'sum.nRNA', 0L)	
	#	number of CD4 per year after year of diagnosis	
	tmp2		<- merge(subset(info, select=c(Patient,AnyPos_T1,DateDied)),subset(df.immu, select=c(Patient,PosCD4)),by='Patient',all.x=TRUE,allow.cartesian=TRUE)
	set(tmp2, NULL, 'PosCD4', floor(hivc.db.Date2numeric(tmp2[,PosCD4])))
	set(tmp2, NULL, 'AnyPos_T1', floor(tmp2[,AnyPos_T1]))
	set(tmp2, NULL, 'DateDied', ceiling(tmp2[,DateDied]))
	set(tmp2, tmp2[, which(is.na(DateDied))], 'DateDied', ceiling(t.endctime))
	tmp2		<- subset( tmp2, is.na(PosCD4) | PosCD4>=AnyPos_T1 )[, 	{
				if( all(is.na(PosCD4)) )
					ans	<- list(t=seq.int(AnyPos_T1[1], DateDied[1]-1), nCD4=rep(0,DateDied[1]-AnyPos_T1[1]))
				else
					ans	<- list(t=seq.int(AnyPos_T1[1], DateDied[1]-1), nCD4=c(rep(0, min(PosCD4)-AnyPos_T1[1]), tabulate(PosCD4-min(PosCD4)+1, nbins=DateDied[1]-min(PosCD4))))
				ans
			}, by='Patient']
	tmp2		<- subset(tmp2, t>=t.startctime)
	tmp			<- merge(tmp, tmp2, by=c('Patient','t'))
	info		<- merge( info, tmp[, list(mean.nCD4.py=mean(nCD4), sum.nCD4=sum(nCD4)) ,by='Patient'], by='Patient', all.x=TRUE )	
	set(info, info[, which(is.na(mean.nCD4.py))], 'mean.nCD4.py', 0.)
	set(info, info[, which(is.na(sum.nCD4))], 'sum.nCD4', 0L)
	#
	info.t		<- tmp[, list(	nPatient=length(Patient), 
					mean.nRNA.pP=mean(nRNA), one.RNA.pP=mean(nRNA>0), 
					mean.nCD4.pP=mean(nCD4), one.CD4.pP=mean(nCD4>0),
					one.contact.pP=mean((nRNA+nCD4)>0)), by='t']		
	info.t
}	
######################################################################################
#
#	three purposes
#	mode1	if df.tpairs is not null, return the transmitter-infected-time period triplets for all **sequence** pairs in df.tpairs
#	mode2	if ri is not null and nsample is NA, return the triplets for all **patient** pairs of infected in ri and any transmitter in df.all
#	mode3	if ri is not null and nsample is not NA, return a **sample** of the triplets for all **sequence pairs** of the sequences belonging to ri and to any transmitter in df.all
#lRNA.supp=method.lRNA.supp; method.minLowerUWithNegT=method.minLowerUWithNegT; t.period=t.period; t.endctime=t.endctime;
project.athena.Fisheretal.YX.part1<- function(df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=NULL, df.tpairs=NULL, tperiod.info=NULL, lRNA.supp=3, t.period=0.25, t.endctime=2013., sample.n=NA, sample.exclude=NULL, method.minLowerUWithNegT=TRUE, save.file=NA, resume=FALSE)
{
	#df.all<- df.all.allmsm; df.immu<-df.immu.allmsm; df.viro<-df.viro.allmsm; df.treatment<-df.treatment.allmsm; ri<-ri.SEQ; df.tpairs<- NULL	
	if(resume && !is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(!resume || is.na(save.file) || inherits(readAttempt, "try-error"))
	{					
		if(is.null(df.tpairs) & is.null(ri))	stop('expect either ri or df.tpairs as input')
		if(is.null(df.tpairs) & !is.na(sample.n))
		{						
			df.tpairs	<- subset(df.all, select=c(Patient, FASTASampleCode))				
			setnames(df.tpairs, colnames(df.tpairs), paste('t.',colnames(df.tpairs),sep=''))
		}		
		if(is.null(df.tpairs) & is.na(sample.n))
		{						
			df.tpairs	<- unique(subset(df.all, select=Patient))				
			setnames(df.tpairs, colnames(df.tpairs), paste('t.',colnames(df.tpairs),sep=''))
		}		
		if(is.null(ri))
			ri			<- unique(subset(df.tpairs, select=Patient))
		X.incare				<- project.athena.Fisheretal.X.incare(df.tpairs, df.all, df.viro, df.immu, df.treatment, lRNA.supp=lRNA.supp, t.period=t.period, t.endctime=t.endctime)
		set(X.incare, X.incare[, which(is.na(lRNA_T1.supp))], 'lRNA_T1.supp', t.endctime)		
		X.incare				<- project.athena.Fisheretal.X.nocontact(X.incare, df.viro, df.immu, df.tpairs, df.all, contact.grace=0.5, t.period=t.period, t.endctime= t.endctime)		
		X.incare				<- project.athena.Fisheretal.X.CDCC(X.incare, df.tpairs, df.all, t.period=t.period, t.endctime=t.endctime)
		X.incare				<- project.athena.Fisheretal.X.followup(X.incare, df.all, df.immu, t.period=t.period, t.endctime=t.endctime)
		X.b4care				<- project.athena.Fisheretal.X.b4care(df.tpairs, df.all, predict.t2inf, t2inf.args, t.period=t.period, method.minLowerUWithNegT=method.minLowerUWithNegT)
		tmp						<- merge(X.incare, X.b4care, by=c('t.Patient','t'))
		cat(paste('\nnumber entries (Patient,t) that overlap in and before care [should be zero], n=',nrow(tmp)))		
		X.pt					<- merge(X.incare, X.b4care, by=c('t.Patient','t'), all.x=1, all.y=1)
		set(X.pt, X.pt[,which(is.na(stage))], 'stage', 'U')							
		X.pt[, AnyT_T1:=NULL]
		X.pt[, AnyPos_T1:=NULL]
		X.pt[, AnyPos_a:=NULL]
		X.pt[, isAcute:=NULL]
		tmp						<- unique(subset( df.all, select=c(Patient, DateBorn, AnyPos_T1, AnyT_T1, isAcute) ))
		set(tmp, NULL, 'DateBorn', tmp[, AnyPos_T1-DateBorn])
		setnames( tmp, 'DateBorn', 'AnyPos_a')
		setnames( tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		X.pt					<- merge(X.pt, tmp, by='t.Patient', allow.cartesian=TRUE)		
		tmp						<- project.athena.Fisheretal.X.ART.pulsed(df.tpairs, df.all, df.treatment, t.pulse.ART=0.75, t.pulse.sc=NA, t.pulse.ART.I=1, t.pulse.ART.I.d=1)
		X.pt					<- merge( X.pt, tmp, by='t.Patient', all.x=1 )
		tmp						<- project.athena.Fisheretal.X.time.diag2suppressed(df.tpairs, df.all, df.viro, lRNA.suppressed= log10(1e3), t2.vl.supp.p=c(0.1, 0.25))
		X.pt					<- merge( X.pt, tmp, by='t.Patient', all.x=1 )
		tmp						<- project.athena.Fisheretal.X.time.diag2firstVLandCD4(df.tpairs, df.all, df.viro, df.immu, t2.care.t1.q=c(0.25,0.5))	
		X.pt					<- merge( X.pt, tmp, by='t.Patient', all.x=1 )		
		#	compute infection window of recipient for direct potential transmitters ri	
		Y.infwindow				<- project.athena.Fisheretal.Y.infectiontime( ri, df.all, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, score.set.value=1, method='for.infected', method.minLowerUWithNegT=TRUE)
		Y.infwindow[, AnyPos_a:=NULL]
		#	some plots
		if(0)
		{
			tmp	<- X.b4care[,  {
						z<- sapply(c(0.03, 0.05, 0.1, 0.2), function(x) max(which(U.score>=x)) )								
						list(U=t[z], Q=c(0.03, 0.05, 0.1, 0.2) , AnyPos_T1=AnyPos_T1[1], AnyPos_a=AnyPos_a[1], isAcute=isAcute[1])	
					}, by='t.Patient']
			tmp[, t2D:= tmp[, AnyPos_T1-U]]
			set(tmp,NULL,'Q',tmp[,factor(Q)]) 
			tmp	<- subset(tmp, U>1996.5 & U<2011)
			ggplot(tmp, aes(x=U, y=t2D, colour=Q)) + geom_point(data=subset(tmp, Q==0.1), position=position_jitter(w = 0.1, h=0), alpha=0.65) + geom_smooth() +
					scale_y_continuous(breaks=seq(0,20,1)) + scale_x_continuous(breaks=seq(1996,2020,1)) + labs(y='time to diagnosis\n(years)', x='time of HIV infection', colour='quantile\nparameter') +
					scale_colour_brewer(palette='Set2') + theme_bw()
			file	<- '/Users/Oliver/duke/2014_HIVMSMbyProp/fig/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Wed_Dec_18_11:37:00_2013_2011_Time2Diag.pdf'
			ggsave(file=file, w=8, h=6)			
		}
		if('FASTASampleCode'%in%colnames(df.tpairs))		#mode 1
		{			
			tmp						<- subset(df.tpairs, select=c(t.Patient, t.FASTASampleCode, FASTASampleCode))
			X.pt					<- merge(X.pt, tmp, by='t.Patient', allow.cartesian=TRUE)									#merge all seq pairs
			tmp						<- subset(df.tpairs, select=c(Patient, t.FASTASampleCode, FASTASampleCode))
			Y.infwindow				<- merge(Y.infwindow, tmp, by='Patient', allow.cartesian=TRUE)							#merge all seq pairs
			YX.part1				<- merge(Y.infwindow, X.pt, by=c('FASTASampleCode','t.FASTASampleCode','t'), allow.cartesian=TRUE)
			set(YX.part1,NULL,'FASTASampleCode', YX.part1[, as.character(FASTASampleCode)])
			set(YX.part1,NULL,'t.FASTASampleCode', YX.part1[, as.character(t.FASTASampleCode)])
			YX.part1				<- project.athena.Fisheretal.YX.age(YX.part1, df.all, t.period=t.period)
			gc()
			YX.part1				<- project.athena.Fisheretal.YX.Trm.Region(YX.part1, df.all)
			gc()
			setkey(YX.part1, FASTASampleCode, t.FASTASampleCode, t)
			YX.part1				<- unique(YX.part1)
		}			
		else if(is.na(sample.n))							#mode 2
		{			
			cat(paste('\nstart big merge of patient pairs'))
			#	the BIG merge. this is easier if we first reduce to the relevant time periods
			X.b4care				<- X.incare	<- X.ARTpulsed<- X.t2.vlsupp<- X.t2.care<- tmp<- NULL
			gc()
			X.pt					<- merge(X.pt, unique(subset( Y.infwindow, select=t )), by='t')		
			YX.part1				<- merge(Y.infwindow, X.pt, by='t', allow.cartesian=TRUE)
			gc()
			cat(paste('\ncompleted big merge of patient pairs, nrows=',nrow(YX.part1)))
			YX.part1				<- project.athena.Fisheretal.YX.age(YX.part1, df.all, t.period=t.period)
			gc()
			YX.part1				<- project.athena.Fisheretal.YX.Trm.Region(YX.part1, df.all)					
			setkey(YX.part1, Patient, t.Patient, t)
			YX.part1				<- unique(YX.part1)
		}	
		else if(!is.na(sample.n))							#mode 3
		{
			cat(paste('\nstart sampled merge of sequence pairs'))
			tmp		<- data.table(t.Patient=setdiff( df.tpairs[, unique(t.Patient)], sample.exclude[, t.Patient]))		#t.Patients not in sample.exclude
			X.pt	<- merge(X.pt, tmp, by='t.Patient')
			X.pt	<- merge(X.pt, unique(subset( Y.infwindow, select=t )), by='t')
			cat(paste('\nt.Patient timelines not in sample.exclude has nrows=',nrow(X.pt)))
			#	split in batches to make searching and merging easier
			df.ntpairs		<- subset(X.pt, select=c(t,t.Patient))
			ri.batch		<- c( seq.int(1, nrow(ri), by=ceiling(nrow(ri)/10)), ifelse( nrow(ri)%%10, nrow(ri), NULL) )
			YX.part1.batch	<- lapply( seq_along(ri.batch)[-1],function(i)
					{
						YX.part1.batch			<- merge(X.pt, merge(Y.infwindow, subset( ri[seq.int(ri.batch[i-1],ri.batch[i]),], select='Patient'), by='Patient'), by='t', allow.cartesian=TRUE )
						df.ntpairs				<- subset(YX.part1.batch, select=c(Patient, t.Patient))
						setkey(df.ntpairs, Patient, t.Patient)
						df.ntpairs				<- unique(df.ntpairs)			
						df.ntpairs				<- df.ntpairs[sample(seq_len(nrow(df.ntpairs)), min(nrow(df.ntpairs),ceiling(sample.n/10))),]		#subsample Patient, t.Patient pairs so that we get intact t.Patient timelines
						YX.part1.batch			<- merge(YX.part1.batch, df.ntpairs, by=c('Patient','t.Patient'))
						cat(paste('\nbatch=',i-1,'YX.batch has nrows=',nrow(YX.part1.batch)))
						YX.part1.batch
					})
			YX.part1		<- do.call('rbind',YX.part1.batch)
			#	need to add FASTASampleCodes			
			YX.part1		<- merge(YX.part1, merge(ri, subset(df.all, select=c(Patient, FASTASampleCode)), by='Patient'), by='Patient', allow.cartesian=TRUE)	
			tmp						<- subset(df.tpairs, select=c(t.Patient, t.FASTASampleCode))
			setkey(tmp, t.Patient, t.FASTASampleCode)
			YX.part1		<- merge(YX.part1, unique(tmp), by='t.Patient', allow.cartesian=TRUE)
			YX.part1.batch	<- NULL
			gc()
			cat(paste('\ncompleted sampled merge of sequence pairs, nrows=',nrow(YX.part1)))
			set(YX.part1,NULL,'FASTASampleCode', YX.part1[, as.character(FASTASampleCode)])
			set(YX.part1,NULL,'t.FASTASampleCode', YX.part1[, as.character(t.FASTASampleCode)])
			YX.part1		<- project.athena.Fisheretal.YX.age(YX.part1, df.all, t.period=t.period)		
			gc()
			YX.part1		<- project.athena.Fisheretal.YX.Trm.Region(YX.part1, df.all)
			gc()
			setkey(YX.part1, FASTASampleCode, t.FASTASampleCode, t)
			YX.part1		<- unique(YX.part1)
		}	
		Y.infwindow<- X.b4care<- X.incare<- X.pt<- X.ARTpulsed<- X.t2.vlsupp<- X.t2.care<- tmp<- NULL
		gc()
		if(!is.null(tperiod.info))
		{
			#set t.period based on *recipient*			
			YX.part1[, t.period:= cut(YX.part1[,AnyPos_T1], breaks=c(tperiod.info[,t.period.min],tperiod.info[nrow(tperiod.info),t.period.max]), labels=seq.int(1,nrow(tperiod.info)), right=FALSE)]
			YX.part1		<- merge(YX.part1, tperiod.info, by='t.period') 			
		}		
		if(!is.na(save.file))
		{						
			#
			gc()
			cat(paste('\nsave YX.part1 to file=',save.file))
			save(YX.part1, file=save.file)
			cat('\nsave completed')
		}			
	}
	YX.part1
}
######################################################################################
project.athena.Fisheretal.RI.summarize<- function(ri, df.all, tperiod.info, info=NULL, with.seq=TRUE, with.cluster=TRUE, cols=c('lRNA_T1','CD4_T1','AnyPos_A','SC_dt','PosCD4_dt','PoslRNA_dt'))
{
	# with.seq=FALSE; with.cluster=FALSE; cols=c('lRNA_T1','CD4_T1','AnyPos_A','SC_dt','PosCD4_dt','PoslRNA_dt'); adjust.AcuteByNegT=0.75
	#	clinical & demographic 
	if(!is.null(info))
	{
		tmp	<- subset( info, select=c(Patient, DateBorn, Sex, AnyPos_T1,  RegionHospital, isAcute, Trm, DateFirstEverCDCC, PoslRNA_T1, lRNA_T1, PosCD4_T1, CD4_T1 ))
		setkey(tmp, Patient)	
		ri	<- merge( ri, unique(tmp), by='Patient' )
		tmp	<- subset( df.all, select=c(Patient, NegT))
		setkey(tmp, Patient)	
		ri	<- merge( ri, unique(tmp), by='Patient' )	
		if(class(ri$NegT)=='Date')
			set(ri, NULL, 'NegT', hivc.db.Date2numeric(ri[,NegT]))			
	}
	if(is.null(info))
	{
		tmp	<- subset( df.all, select=c(Patient, DateBorn, Sex, AnyPos_T1,  RegionHospital, isAcute, Trm, NegT, DateFirstEverCDCC, PoslRNA_T1, lRNA_T1, PosCD4_T1, CD4_T1 ))
		setkey(tmp, Patient)	
		ri	<- merge( ri, unique(tmp), by='Patient' )	
		ri	<- subset( ri, isAcute=='Yes' | isAcute=='Maybe' )			
	}
	#	number of Sequences 
	if(with.seq)
	{
		tmp			<- merge( subset(ri, select=Patient), subset( info, select=c(Patient, PosSeqT)), by='Patient' )
		ri			<- merge( ri, tmp[, list(n.seq= length(PosSeqT), PosSeq_T1= min(PosSeqT)) , by='Patient'], by='Patient')
		ri[, PosSeq_dt:=PosSeq_T1-AnyPos_T1]
		if(!'PosSeq_dt'%in%cols)
			cols	<- c( cols,'PosSeq_dt' )
	}
	#	cluster
	if(with.cluster)
	{
		tmp	<- subset( info, select=c(Patient, cluster, clu.fPossAcute, clu.npat ))
		setkey(tmp, Patient)	
		ri	<- merge( ri, unique(tmp), by='Patient' )			
	}
	#	age at diagnosis
	ri[, AnyPos_A:= AnyPos_T1-DateBorn]
	#	time period
	ri[, t.period:= factor(cut( ri[, AnyPos_T1], breaks=c( tperiod.info[, t.period.min], tail( tperiod.info[, t.period.max], 1 ) ), right=FALSE, label= seq_len(nrow(tperiod.info)) ))]
	ri	<- merge(ri, tperiod.info, by='t.period')	
	#	seroconversion interval
	ri[, SC_dt:=AnyPos_T1-NegT]		
	#	time differences
	ri[, PosCD4_dt:=PosCD4_T1-AnyPos_T1]
	ri[, PoslRNA_dt:=PoslRNA_T1-AnyPos_T1]		
	#	number with sequence
	cat(paste('\n ri with sequence=',nrow( merge(ri, unique(subset(df.all, !is.na(PosSeqT.min), select=c(Patient, PosSeqT.min))), by='Patient') )))
	#	summarize by time period: quantiles 
	ris	<- ri[, {
				tmp	<- lapply(.SD, quantile, prob=c(0.05, 0.5, 0.95), na.rm=TRUE)
				tmp[['stat']]<- paste('quantile',c(0.05, 0.5, 0.95),sep='_')
				tmp
			}, by='t.period',.SDcols=cols ]
	#	summarize by time period: sd
	ris	<- rbind(ris, ri[, {
						tmp	<- lapply(.SD, sd, na.rm=TRUE)
						tmp[['stat']]<- 'sd'
						tmp
					}, by='t.period',.SDcols=cols ])
	#	summarize by time period: total not NA
	ris	<- rbind(ris, ri[, {
						tmp	<- lapply(.SD, function(x) length(which(!is.na(x))))
						tmp[['stat']]<- 'not.NA'
						tmp
					}, by='t.period',.SDcols=cols ])
		
	#	summarize by time period: tables
	tmp				<- ri[ , table( t.period, RegionHospital, useNA='ifany') ]	
	rin				<- cbind( data.table(t.period=factor(rownames(tmp))), as.data.table(unclass(tmp)) )
	tmp				<- ri[ , table( t.period, isAcute, useNA='ifany') ]
	colnames(tmp)	<- paste('isAcute',colnames(tmp),sep='.')
	rin				<- merge( rin, cbind( data.table(t.period=factor(rownames(tmp))), as.data.table(unclass(tmp)) ), by='t.period' )
	tmp				<- ri[ , table( t.period, Trm, useNA='ifany') ]
	colnames(tmp)	<- paste('Trm',colnames(tmp),sep='.')
	rin				<- merge( rin, cbind( data.table(t.period=factor(rownames(tmp))), as.data.table(unclass(tmp)) ), by='t.period' )	
	tmp				<- ri[ , table( t.period, Sex, useNA='ifany') ]
	colnames(tmp)	<- paste('Sex',colnames(tmp),sep='.')
	rin				<- merge( rin, cbind( data.table(t.period=factor(rownames(tmp))), as.data.table(unclass(tmp)) ), by='t.period' )	
	#	
	tmp		<- merge(tperiod.info, ri[, list(Patient=length(Patient)), by='t.period'], by='t.period')
	rin		<- merge(tmp, rin, by='t.period')
	ris		<- merge(tmp, ris, by='t.period')
	#	get proportions and confidence intervals	
	require(Hmisc)
	cols	<- setdiff( colnames(rin), c('t.period','t.period.min','t.period.max') )
	rip		<- rin[ , 	{
							tmp				<- lapply(.SD, function(x) as.double(binconf( x, Patient, alpha=0.05, method= "wilson", include.x=FALSE, include.n=FALSE, return.df=FALSE)) )								
							tmp[['stat']]	<- c('p','p.l95','p.u95')
							tmp
						}, by='t.period', .SDcols=cols]
	list(rin=rin, ris=ris, rip=rip)
}
######################################################################################
project.athena.Fisheretal.numbers<- function()
{
	require(data.table)
	require(ape)
	#stop()
	indir					<- paste(DATA,"fisheretal_data",sep='/')		
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	outdir					<- paste(DATA,"fisheretal",sep='/')
	infile.cov.study		<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infile.viro.study		<- paste(indircov,"ATHENA_2013_03_Viro.R",sep='/')
	infile.immu.study		<- paste(indircov,"ATHENA_2013_03_Immu.R",sep='/')
	infile.treatment.study	<- paste(indircov,"ATHENA_2013_03_Regimens.R",sep='/')
	infile.cov.all			<- "ATHENA_2013_03_AllSeqPatientCovariates_AllMSM"
	infile.viro.all			<- paste(indircov,"ATHENA_2013_03_Viro_AllMSM.R",sep='/')
	infile.immu.all			<- paste(indircov,"ATHENA_2013_03_Immu_AllMSM.R",sep='/')
	infile.treatment.all	<- paste(indircov,"ATHENA_2013_03_Regimens_AllMSM.R",sep='/')		
	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))	
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	resume					<- 1
	verbose					<- 1
	if(0)
	{
		method					<- '3c'
		method.recentctime		<- '2013-03-01'
		method.nodectime		<- 'any'
		method.risk				<- 'm21st.cas'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"					
		clu.infilexml.opt		<- "mph4clutx4tip"
		clu.infilexml.template	<- "um192rhU2080"	
		outfile					<- paste(infile,'_Ac=MY_D=35_gmrf',sep='')
	}
	if(0)
	{
		method					<- '3d'
		method.recentctime		<- '2013-03-01'
		method.nodectime		<- 'any'
		method.risk				<- 'm21st.cas'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"					
		clu.infilexml.opt		<- "mph4clutx4tip"
		clu.infilexml.template	<- "um192rhU2080"	
		outfile					<- paste(infile,'_Ac=MY_D=35_gmrf',sep='')
	}	
	if(0)
	{		
		method					<- '3c'
		method.recentctime		<- '2013-03-01'
		method.nodectime		<- 'any'
		method.risk				<- 'm21st.cas'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		clu.infilexml.opt		<- "clrh80"
		clu.infilexml.template	<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'_Ac=MY_D=35_sasky',sep='')
	}
	if(0)
	{		
		method					<- '3d'
		method.recentctime		<- '2013-03-01'
		method.nodectime		<- 'any'
		method.risk				<- 'm2B1st.cas'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		clu.infilexml.opt		<- "clrh80"
		clu.infilexml.template	<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'_Ac=MY_D=35_sasky',sep='')
	}	
	if(1)
	{		
		method					<- '3d'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		method.risk				<- 'm2B1st.cas'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		clu.infilexml.opt		<- "clrh80"
		clu.infilexml.template	<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'_Ac=MY_D=35_sasky',sep='')
	}
	
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
		if(length(tmp)>0) infile.cov.study<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infiletree= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infiletree<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,14),
									infilexml.opt= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) clu.infilexml.opt<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									infilexml.template= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) clu.infilexml.template<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									outfile= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile<- tmp[1]						
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									method= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) method<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									method.risk= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) method.risk<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,17),
									method.nodectime= return(substr(arg,19,nchar(arg))),NA)	}))
		if(length(tmp)>0) method.nodectime<- tmp[1]			
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									method.recentctime= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) method.recentctime<- tmp[1]	
		
		
	}	
	clu.infile			<- infile
	clu.indir			<- indir
	clu.insignat		<- insignat	
	t.recent.endctime	<- hivc.db.Date2numeric(as.Date(method.recentctime))	
	t.recent.endctime	<- floor(t.recent.endctime) + floor( (t.recent.endctime%%1)*100 %/% (t.period*100) ) * t.period	
	outfile				<- paste( outfile, ifelse(t.recent.endctime==t.endctime,'',paste('_',t.recent.endctime,sep='')), sep='')
	
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infile.cov.study)
		print(infiletree)
		print(clu.infilexml.opt)
		print(clu.infilexml.template)
		print(outdir)
		print(outfile)
		print(method)
		print(method.risk)
		print(method.nodectime)
	}
	if(method.nodectime=='any')
		method		<- paste(method,'a',sep='')
	if(method.nodectime=='map')
		method		<- paste(method,'m',sep='')	
	adjust.AcuteByNegT=0.75
	if(resume)
	{		
		files		<- list.files(outdir)		
		files		<- files[ sapply(files, function(x) grepl(outfile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(paste('Yscore',method,sep=''), x, fixed=1) & !grepl('tables', x, fixed=1) & grepl(paste(method.risk,'.R',sep=''),x, fixed=1)  ) ]		
		stopifnot(length(files)==0)		
	}
	ans				<- list()
	#
	#	get rough idea about (backward) time to infection from time to diagnosis, taking midpoint of SC interval as 'training data'
	#
	tmp				<- project.athena.Fisheretal.t2inf(indircov, infile.cov.study, adjust.AcuteByNegT=0.75, adjust.dt.CD4=1, adjust.AnyPos_y=2003, adjust.NegT=2)
	predict.t2inf	<- tmp$predict.t2inf
	t2inf.args		<- tmp$t2inf.args
	#
	#	get data relating to full population (MSM including those without seq)
	#
	tmp					<- project.athena.Fisheretal.select.denominator(indir, infile, insignat, indircov, infile.cov.all, infile.viro.all, infile.immu.all, infile.treatment.all, infiletree=NULL, adjust.AcuteByNegT=adjust.AcuteByNegT, adjust.NegT4Acute=NA, adjust.NegTByDetectability=0.25, adjust.minSCwindow=0.25, adjust.AcuteSelect=c('Yes','Maybe'), t.recent.endctime=t.recent.endctime)	
	df.all.allmsm		<- tmp$df.all	
	df.viro.allmsm		<- tmp$df.viro
	df.immu.allmsm		<- tmp$df.immu
	df.treatment.allmsm	<- tmp$df.treatment
	tmp					<- tmp$df.select
	setkey(tmp, Patient)
	ri.allmsm			<- unique(tmp)	
	ri.allmsm			<- subset(ri.allmsm, AnyPos_T1>=1994.416 & AnyPos_T1<2011)
	tmp					<- ri.allmsm
	ans$ri				<- data.table(ri='allmsm', recent.n=nrow(tmp), acute.n= tmp[, length(which(isAcute=='Yes'))], MSM.n=tmp[, length(which(Trm=='MSM'))], BI.n=tmp[, length(which(Trm=='BI'))] )
	save.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RIMSM_',method,'_tATHENAmsm','.R',sep='')
	tmp					<- project.athena.Fisheretal.YX.part1(df.all.allmsm, df.immu.allmsm, df.viro.allmsm, df.treatment.allmsm, predict.t2inf, t2inf.args, ri=ri.allmsm, df.tpairs=NULL, tperiod.info=NULL, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)
	ans$pyiw			<- data.table(pop='allmsm', t.py=nrow(tmp)*t.period, t.n=tmp[, length(unique(t.Patient))], n=tmp[, length(unique(Patient))])
	tmp					<- NULL
	gc()	
	#
	#	get data relating to study population (subtype B sequ)
	#
	tmp				<- project.athena.Fisheretal.select.denominator(indir, infile, insignat, indircov, infile.cov.study, infile.viro.study, infile.immu.study, infile.treatment.study, infiletree=infiletree, adjust.AcuteByNegT=adjust.AcuteByNegT, adjust.NegT4Acute=NA, adjust.NegTByDetectability=0.25, adjust.minSCwindow=0.25, adjust.AcuteSelect=c('Yes','Maybe'), t.recent.endctime=t.recent.endctime)	
	df.all			<- tmp$df.all
	df.denom		<- tmp$df.select
	df.viro			<- tmp$df.viro
	df.immu			<- tmp$df.immu
	df.treatment	<- tmp$df.treatment	
	clumsm.subtrees	<- tmp$clumsm.subtrees
	clumsm.info		<- tmp$clumsm.info
	clumsm.ph		<- tmp$clumsm.ph
	setkey(clumsm.info, cluster)
	#
	df.tpairs		<- project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos(clumsm.info, df.denom, any.pos.grace.yr= 3.5, select.if.transmitter.seq.unique=FALSE)
	#
	#	get time stamped data (if clusters missing, confine df.tpairs to available clusters)
	#
	tmp						<- project.athena.Fisheretal.get.dated.phylo.for.selection(df.tpairs, clu.indir, clu.infile, clu.insignat, clu.infilexml.opt, clu.infilexml.template, method.nodectime=method.nodectime)
	cluphy.map.nodectime	<- tmp$clu$cluphy.map.nodectime
	cluphy.subtrees			<- tmp$clu$cluphy.subtrees
	cluphy.info				<- tmp$clu$cluphy.info
	cluphy					<- tmp$clu$cluphy		
	#	
	#	RI with at least one sequence
	#	
	ri.allseq			<- subset(merge(subset(df.all,select=c(Patient)), ri.allmsm, by='Patient'), !is.na(PosSeqT))
	setkey(ri.allseq, Patient)
	ri.allseq			<- unique(ri.allseq)
	tmp					<- ri.allseq
	ans$ri				<- rbind(ans$ri, data.table(ri='allseq', recent.n=nrow(tmp), acute.n= tmp[, length(which(isAcute=='Yes'))], MSM.n=tmp[, length(which(Trm=='MSM'))], BI.n=tmp[, length(which(Trm=='BI'))] ))
	save.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RISEQ_',method,'_tATHENAseq','.R',sep='')
	tmp					<- project.athena.Fisheretal.YX.part1(df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=ri.allseq, df.tpairs=NULL, tperiod.info=NULL, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)
	ans$pyiw			<- rbind(ans$pyiw, data.table(pop='allseq', t.py=nrow(tmp)*t.period, t.n=tmp[, length(unique(t.Patient))], n=tmp[, length(unique(Patient))]))
	tmp					<- NULL
	gc()
	#	
	#	RI with at least one sequence > 500 & no evidence for recombination
	#	
	load( paste(indir,'/', infile, '_', gsub('/',':',insignat), '.R', sep='') )
	ri.allph			<- merge( data.table(FASTASampleCode=rownames(seq.PROT.RT)), ri.allseq, by='FASTASampleCode' )
	tmp					<- ri.allph
	ans$ri				<- rbind(ans$ri, data.table(ri='allph', recent.n=nrow(tmp), acute.n= tmp[, length(which(isAcute=='Yes'))], MSM.n=tmp[, length(which(Trm=='MSM'))], BI.n=tmp[, length(which(Trm=='BI'))] ))
	save.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RIPH_',method,'_tATHENAseq','.R',sep='')
	tmp					<- project.athena.Fisheretal.YX.part1(df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=ri.allph, df.tpairs=NULL, tperiod.info=NULL, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)
	ans$pyiw			<- rbind(ans$pyiw, data.table(pop='allph', t.py=nrow(tmp)*t.period, t.n=tmp[, length(unique(t.Patient))], n=tmp[, length(unique(Patient))]))
	tmp					<- NULL
	gc()		
	#	
	#	RI in MSM cluster
	#	
	ri.allclu			<- df.denom
	setkey(ri.allclu, Patient)
	ri.allclu			<- unique(ri.allclu)
	tmp					<- ri.allclu
	ans$ri				<- rbind(ans$ri, data.table(ri='allclu', recent.n=nrow(tmp), acute.n= tmp[, length(which(isAcute=='Yes'))], MSM.n=tmp[, length(which(Trm=='MSM'))], BI.n=tmp[, length(which(Trm=='BI'))] ))
	save.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RICLU_',method,'_tATHENAseq','.R',sep='')
	tmp					<- project.athena.Fisheretal.YX.part1(df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=ri.allclu, df.tpairs=NULL, tperiod.info=NULL, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)	
	ans$pyiw			<- rbind(ans$pyiw, data.table(pop='allclu', t.py=nrow(tmp)*t.period, t.n=tmp[, length(unique(t.Patient))], n=tmp[, length(unique(Patient))]))
	tmp					<- NULL
	gc()		
	#	
	#	RI in MSM cluster that meet AnyPos_T1 < 3.5 yr
	#	
	df.tpairs			<- project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos(clumsm.info, df.denom, any.pos.grace.yr= 3.5, select.if.transmitter.seq.unique=FALSE)
	save.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RICT_',method,'_tATHENAclu','.R',sep='')	
	tmp					<- project.athena.Fisheretal.YX.part1(df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=NULL, df.tpairs=df.tpairs, tperiod.info=NULL, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)
	ans$pyiw			<- rbind(ans$pyiw, data.table(pop='with.Diag', t.py=nrow(tmp)*t.period, t.n=tmp[, length(unique(t.Patient))], n=tmp[, length(unique(Patient))]))	
	#	
	#	RI in MSM cluster with AnyPos_T1 < 3.5 yr + genetic distance smaller than 0.1% quantile of unlinked by death
	#	
	save.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'YX',method,'_all.R',sep='')
	tmp					<- merge( tmp, subset( df.tpairs, select=c(FASTASampleCode, t.FASTASampleCode, cluster) ), by=c('FASTASampleCode','t.FASTASampleCode'), all.x=1)
	tmp[, class:='pt']
	tmp					<- project.athena.Fisheretal.YX.part2(tmp, df.all, predict.t2inf, t2inf.args, indir, insignat, indircov, infile.cov.study, infiletree, outdir, outfile, cluphy=cluphy, cluphy.info=cluphy.info, cluphy.map.nodectime=cluphy.map.nodectime, df.tpairs.4.rawbrl=df.tpairs, rm.zero.score=TRUE, t.period=t.period, save.file=save.file, save.all=TRUE, resume=resume, method=method)
	YX.tpairs			<- tmp$YX.tpairs
	df.all				<- tmp$df.all
	Y.brl				<- tmp$Y.brl
	Y.U					<- tmp$Y.U
	Y.coal				<- tmp$Y.coal
	plot.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'missed_',method,'.pdf',sep='')		
	tmp					<- project.athena.Fisheretal.Y.rm.missedtransmitter(YX.tpairs, df.all, Y.brl, Y.U, Y.coal=Y.coal, cut.date=0.5, cut.brl=1e-3, any.pos.grace.yr= 3.5, rm.zero.score=TRUE, plot.file=plot.file, pyiw=ans$pyiw)			
	ans$pyiw			<- tmp$pyiw
	
	save.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'numbers_',method,'.R',sep='')
	cat(paste('save numbers to ',save.file))
	save(ans, file=save.file)	
}
######################################################################################
hivc.prog.betareg.estimaterisks<- function()
{
	require(reshape2)
	require(data.table)
	require(ape)
	#stop()
	indir					<- paste(DATA,"fisheretal_data",sep='/')		
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	outdir					<- paste(DATA,"fisheretal",sep='/')
	infile.cov.study		<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infile.viro.study		<- paste(indircov,"ATHENA_2013_03_Viro.R",sep='/')
	infile.immu.study		<- paste(indircov,"ATHENA_2013_03_Immu.R",sep='/')
	infile.treatment.study	<- paste(indircov,"ATHENA_2013_03_Regimens.R",sep='/')	
	infile.cov.all			<- "ATHENA_2013_03_AllSeqPatientCovariates_AllMSM"
	infile.viro.all			<- paste(indircov,"ATHENA_2013_03_Viro_AllMSM.R",sep='/')
	infile.immu.all			<- paste(indircov,"ATHENA_2013_03_Immu_AllMSM.R",sep='/')
	infile.treatment.all	<- paste(indircov,"ATHENA_2013_03_Regimens_AllMSM.R",sep='/')		
	infile.trm.model		<- NA
	t.period				<- 1/8
	t.recent.startctime		<- hivc.db.Date2numeric(as.Date("1996-07-15"))
	t.recent.startctime		<- floor(t.recent.startctime) + floor( (t.recent.startctime%%1)*100 %/% (t.period*100) ) * t.period
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))	
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	resume					<- 1
	verbose					<- 1
	if(1)
	{		
		method					<- '3n'
		method					<- '3o'
		#method					<- '3m'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		method.risk				<- 'm2Cwmx.wtn.tp4'
		method.Acute			<- 'higher'	#'central'#'empirical'
		method.minQLowerU		<- 0.1
		method.use.AcuteSpec	<- 1
		method.brl.bwhost		<- 2
		method.lRNA.supp		<- 100
		method.thresh.pcoal		<- 0.3
		method.minLowerUWithNegT<- 1
		method.PDT				<- 'SEQ'	# 'PDT'		
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		clu.infilexml.opt		<- "clrh80"
		clu.infilexml.template	<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'_Ac=MY_D=35_sasky',sep='')
	}
	if(0)
	{		
		method					<- '3l'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		method.risk				<- 'm2Bwmx.wtn.tp4'
		method.Acute			<- 'higher'	#'central'#'empirical'
		method.use.AcuteSpec	<- 1
		method.minQLowerU		<- 0.1
		method.brl.bwhost		<- 2
		method.lRNA.supp		<- 100
		method.thresh.pcoal		<- 0.3
		method.minLowerUWithNegT<- 1
		method.PDT				<- 'SEQ'	# 'PDT'		
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		clu.infilexml.opt		<- "clrh80"
		clu.infilexml.template	<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'_Ac=MY_D=35_sasky',sep='')
	}
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
		if(length(tmp)>0) infile.cov.study<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infiletree= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infiletree<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,14),
									infilexml.opt= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) clu.infilexml.opt<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									infilexml.template= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) clu.infilexml.template<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									outfile= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile<- tmp[1]						
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									method= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) method<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									method.risk= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) method.risk<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,13),
									method.Acute= return(substr(arg,15,nchar(arg))),NA)	}))
		if(length(tmp)>0) method.Acute<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,21),
									method.use.AcuteSpec= return(as.numeric(substr(arg,23,nchar(arg)))),NA)	}))
		if(length(tmp)>0) method.use.AcuteSpec<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,18),
									method.minQLowerU= return(as.numeric(substr(arg,20,nchar(arg)))),NA)	}))
		if(length(tmp)>0) method.minQLowerU<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,17),
									method.nodectime= return(substr(arg,19,nchar(arg))),NA)	}))
		if(length(tmp)>0) method.nodectime<- tmp[1]			
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									method.recentctime= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) method.recentctime<- tmp[1]			
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,20),
									method.thresh.pcoal= return(as.numeric(substr(arg,22,nchar(arg)))),NA)	}))
		if(length(tmp)>0) method.thresh.pcoal<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									method.PDT= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) method.PDT<- tmp[1]			
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,17),
									method.lRNA.supp= return(as.numeric(substr(arg,19,nchar(arg)))),NA)	}))
		if(length(tmp)>0) method.lRNA.supp<- tmp[1]		
		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,25),
									method.minLowerUWithNegT= return(as.numeric(substr(arg,27,nchar(arg)))),NA)	}))
		if(length(tmp)>0) method.minLowerUWithNegT<- tmp[1]		
		
		
	}	
	clu.infile			<- infile
	clu.indir			<- indir
	clu.insignat		<- insignat	
	t.recent.endctime	<- hivc.db.Date2numeric(as.Date(method.recentctime))	
	t.recent.endctime	<- floor(t.recent.endctime) + floor( (t.recent.endctime%%1)*100 %/% (t.period*100) ) * t.period	
	outfile				<- paste( outfile, ifelse(t.recent.endctime==t.endctime,'',paste('_',t.recent.endctime,sep='')), sep='')
	
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infile.cov.study)
		print(infiletree)
		print(clu.infilexml.opt)
		print(clu.infilexml.template)
		print(outdir)
		print(outfile)
		print(method)
		print(method.risk)
		print(method.nodectime)
		print(method.PDT)
		print(method.Acute)
		print(method.minQLowerU)
		print(method.use.AcuteSpec)
		print(method.lRNA.supp)
		print(method.thresh.pcoal)
		print(method.minLowerUWithNegT)
	}	
	if(method=='3l')
	{
		infile.trm.model	<- paste(indircov,"TchainBelgium_set7_pol_GAmodel_INFO.R",sep='/')
		cat(paste('\nusing file', infile.trm.model))
	}		
	if(method=='3m')
	{
		infile.trm.model	<- paste(indircov,"TchainBelgium_set7_pol_GAmodel_nA_INFO.R",sep='/')
		cat(paste('\nusing file', infile.trm.model))
	}
	if(method=='3n')
	{
		infile.trm.model	<- paste(indircov,"TchainBelgium_set7_pol_GAmodel_nA_NO_INFO.R",sep='/')
		cat(paste('\nusing file', infile.trm.model))
	}
	if(method=='3o')
	{
		infile.trm.model	<- paste(indircov,"TchainBelgium_set7_pol_GA11model_nA_INFO.R",sep='/')
		cat(paste('\nusing file', infile.trm.model))
	}
	if(method=='3p')
	{
		infile.trm.model	<- paste(indircov,"TchainBelgium_set7_pol_GA11emodel_nA_INFO.R",sep='/')
		cat(paste('\nusing file', infile.trm.model))
	}	
	if(method.nodectime=='any')
		method				<- paste(method,'a',sep='')
	if(method.nodectime=='map')
		method				<- paste(method,'m',sep='')	
	if(method.use.AcuteSpec==0)
		method				<- paste(method, 2, sep='')		#replaces variable method.brl.bwhost since this is now fixed
	if(method.use.AcuteSpec==1)
		method				<- paste(method, method.use.AcuteSpec, sep='')
	if(method.Acute=='empirical')
	{
		dur.Acute			<- c(Yes= 365/2, Maybe=320)	
	}	
	if(method.Acute=='central')
	{
		dur.Acute			<- c(Yes= 2.9*30, Maybe=2.9*30)
		method				<- paste(method,'C',sep='')
	}
	if(method.Acute=='lower')
	{
		dur.Acute			<- c(Yes= 1.23*30, Maybe=1.23*30)
		method				<- paste(method,'L',sep='')
	}
	if(method.Acute=='higher')
	{
		dur.Acute			<- c(Yes= 5.28*30, Maybe=5.28*30)
		method				<- paste(method,'H',sep='')
	}		
	if(method.minQLowerU!=0.01)
		method				<- paste(method, method.minQLowerU*10,sep='')	
	if(method.minQLowerU==0)
		method				<- paste(method,'m',sep='')
	if(method.thresh.pcoal!=0.5)
		method				<- paste(method,'C',method.thresh.pcoal*10,sep='')	
	if(method.lRNA.supp<1e3)
		method				<- paste(method,'V',method.lRNA.supp,sep='')
	if(!method.minLowerUWithNegT)
		method				<- paste(method,'N',method.minLowerUWithNegT,sep='')
	
	adjust.AcuteByNegT		<- 1
	any.pos.grace.yr		<- Inf	
	method.lRNA.supp		<- log10(method.lRNA.supp)	
	#
	#	check if we have precomputed tables
	#
	X.tables			<- NULL
	if(1)
	{
		save.file		<- NA
		if(grepl('m21st',method.risk))		save.file	<- 'm21st'
		if(grepl('m2B1st',method.risk))		save.file	<- 'm2B1st'
		if(grepl('m2t',method.risk))		save.file	<- 'm2t'
		if(grepl('m2Bt',method.risk))		save.file	<- 'm2Bt'
		if(grepl('m2wmx',method.risk))		save.file	<- 'm2wmx'
		if(grepl('m2Bwmx',method.risk))		save.file	<- 'm2Bwmx'
		if(grepl('m2Cwmx',method.risk))		save.file	<- 'm2Cwmx'
		if(grepl('m3.n3mx',method.risk) & !grepl('No',method.risk))								save.file	<- 'm3.n3mx'
		if(grepl('m3.ind',method.risk) & !grepl('No',method.risk))								save.file	<- 'm3.ind'
		if(grepl('m3.ind',method.risk) & grepl('No',method.risk))								save.file	<- 'm3.indNo'
		if(grepl('m3.indmx',method.risk) & !grepl('No',method.risk))							save.file	<- 'm3.indmx'
		if(grepl('m3.indmx',method.risk) & grepl('No',method.risk))								save.file	<- 'm3.indmxNo'
		if(grepl('m3.nnrtpiNo',method.risk))													save.file	<- 'm3.nnrtpiNo'
		if(grepl('m4.Bwmx',method.risk))	save.file	<- 'm4.Bwmx'
		if(grepl('m5.tA',method.risk))		save.file	<- 'm5.tA'
		if(grepl('m5.tAb',method.risk))		save.file	<- 'm5.tAb'
		if(grepl('m5.tAc',method.risk))		save.file	<- 'm5.tAc'
		if(grepl('m5.tiA',method.risk))		save.file	<- 'm5.tiA'
		if(grepl('m5.tiAb',method.risk))	save.file	<- 'm5.tiAb'
		if(grepl('m5.tiAc',method.risk))	save.file	<- 'm5.tiAc'
		if(is.na(save.file))	stop('unknown method.risk')				
		tmp				<- regmatches(method.risk, regexpr('tp[0-9]', method.risk))		
		save.file		<- paste(save.file, ifelse(length(tmp), paste('.',tmp,sep=''), ''), sep='')
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_tables',method.PDT,'_',save.file,'.R',sep='')
		X.tables		<- project.athena.Fisheretal.estimate.risk.table(YX=NULL, X.den=NULL, X.msm=NULL, X.clu=NULL, resume=TRUE, save.file=save.file, method=method.risk)
		if(!is.null(X.tables))
		{
			cat('\nloaded X.tables')
			##	sense check that risk factors have been correctly computed
			nt.table	<- copy(X.tables$nt.table.pt)
			nt.table	<- dcast.data.table(nt.table, t.Patient + risk + factor ~ stat, value.var="nt")		
			tmp			<- nt.table[, which(X.seq>X.msm)]
			if(length(tmp))	cat(paste('\nWARNING: X.seq>X.msm for entries n=',length(tmp)))
			#stopifnot(length(tmp)==0)			
			tmp			<- nt.table[, which(X.clu>X.seq)]
			if(length(tmp))	cat(paste('\nWARNING: X.clu>X.seq for entries n=',length(tmp)))
			#stopifnot(length(tmp)==0)		
			tmp			<- nt.table[, which(YX>X.clu)]
			if(length(tmp))	cat(paste('\nWARNING: YX>X.clu for entries n=',length(tmp)))
			#stopifnot(length(tmp)==0)	#there s one recipient that is just on the boundary for m2Bwmx.tp1 - let pass		
			nt.table	<- NULL
		}
	}
	if(resume)
	{		
		files		<- list.files(outdir)			
		if(length(files))
			files	<- files[ sapply(files, function(x) grepl(outfile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(paste('Yscore',method,sep=''), x, fixed=1) & grepl(paste('denom',method.PDT,sep=''), x) & !grepl('tables', x, fixed=1) & grepl(paste(method.risk,'.R',sep=''),x, fixed=1)  ) ]		
		stopifnot(length(files)==0)		
	}	
	#
	#	get rough idea about (backward) time to infection from time to diagnosis, taking midpoint of SC interval as 'training data'
	#
	plot.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 't2inf',method.PDT,method,sep='')
	tmp				<- project.athena.Fisheretal.t2inf(	indircov, infile.cov.study,
														method.Acute=method.Acute, method.minQLowerU=method.minQLowerU,
														adjust.AcuteByNegT=0.75, adjust.dt.CD4=1, adjust.AnyPos_y=2003, adjust.NegT=2, dur.AcuteYes=dur.Acute['Yes'], dur.AcuteMaybe=dur.Acute['Maybe'], use.AcuteSpec=method.use.AcuteSpec, t.recent.endctime=t.recent.endctime, 
														plot.file=plot.file)
	predict.t2inf	<- tmp$predict.t2inf
	t2inf.args		<- tmp$t2inf.args
	#
	#	get data relating to study population (subtype B sequ)
	#
	tmp				<- project.athena.Fisheretal.select.denominator(indir, infile, insignat, indircov, infile.cov.study, infile.viro.study, infile.immu.study, infile.treatment.study, infiletree=infiletree, adjust.AcuteByNegT=adjust.AcuteByNegT, adjust.NegT4Acute=NA, adjust.NegTByDetectability=0.25, adjust.minSCwindow=0.25, adjust.AcuteSelect=c('Yes','Maybe'), use.AcuteSpec=method.use.AcuteSpec, t.recent.endctime=t.recent.endctime, t.recent.startctime=t.recent.startctime)	
	df.all			<- tmp$df.all	
	df.denom.CLU	<- tmp$df.select
	df.denom.SEQ	<- tmp$df.select.SEQ
	ri.CLU			<- unique(subset(df.denom.CLU, select=Patient))
	ri.SEQ			<- unique(subset(df.denom.SEQ, select=Patient))
	df.viro			<- tmp$df.viro
	df.immu			<- tmp$df.immu
	df.treatment	<- tmp$df.treatment	
	clumsm.subtrees	<- tmp$clumsm.subtrees
	clumsm.info		<- tmp$clumsm.info
	clumsm.ph		<- tmp$clumsm.ph
	setkey(clumsm.info, cluster)	
	#
	#	get data relating to full population (MSM including those without seq)
	#	this merges the patients with HIV 1 B sequences and the MSM patients without a sequence 
	tmp					<- project.athena.Fisheretal.select.denominator(	indir, infile, insignat, indircov, infile.cov.all, infile.viro.all, infile.immu.all, infile.treatment.all, 
																			infiletree=NULL, adjust.AcuteByNegT=adjust.AcuteByNegT, adjust.NegT4Acute=NA, adjust.NegTByDetectability=0.25, adjust.minSCwindow=0.25, adjust.AcuteSelect=c('Yes','Maybe'), use.AcuteSpec=method.use.AcuteSpec,
																			t.recent.endctime=t.recent.endctime, t.recent.startctime=t.recent.startctime,
																			df.viro.part=df.viro, df.immu.part=df.immu, df.treatment.part=df.treatment, df.all.part=df.all)	
	df.all.allmsm		<- tmp$df.all
	df.viro.allmsm		<- tmp$df.viro
	df.immu.allmsm		<- tmp$df.immu
	df.treatment.allmsm	<- tmp$df.treatment
	tmp					<- tmp$df.select.SEQ
	setkey(tmp, Patient)
	ri.ALLMSM			<- unique(tmp)	
	#
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
	#	select potential transmitters within same cluster
	#
	df.tpairs		<- project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos(clumsm.info, df.denom.CLU, any.pos.grace.yr= any.pos.grace.yr, select.if.transmitter.seq.unique=FALSE)	
	#df.tpairs		<- project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos.MinBrlMLETree(df.denom, clumsm.subtrees, any.pos.grace.yr= 2)
	#	
	#	plot MLE tree
	#
	if(0)
	{
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
		file							<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'MLEtree.pdf',sep='')
		pdf(file=file, w=8, h=200)
		plot.phylo(df.tpairs.mleph, show.node.label=1, cex=0.4, no.margin=1, label.offset=0.005, edge.width=0.5)
		dev.off()
	}			
	#	plot number of potential transmitters
	if(0)
	{
		tmp				<- merge( subset(df.tpairs, select=Patient), subset(df.denom.CLU, select=c(Patient, AnyPos_T1)), by='Patient' )
		file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'nrecentlyinfected', '.pdf',sep='')
		pdf(file=file, w=5, h=5)
		par(mar=c(3,5,0.5,0.5))
		barplot( table( tmp[, round(AnyPos_T1)] ), ylab="# recently infected\n with unique potential transmitter" )
		dev.off()	
	}	
	#
	#	get time stamped data (if clusters missing, confine df.tpairs to available clusters)
	#
	tmp						<- project.athena.Fisheretal.get.dated.phylo.for.selection(df.tpairs, clu.indir, clu.infile, clu.insignat, clu.infilexml.opt, clu.infilexml.template, method.nodectime=method.nodectime)
	cluphy.map.nodectime	<- tmp$clu$cluphy.map.nodectime
	cluphy.subtrees			<- tmp$clu$cluphy.subtrees
	cluphy.info				<- tmp$clu$cluphy.info
	cluphy					<- tmp$clu$cluphy	
	if(0)
	{
		#	anonymize
		setkey(clumsm.info, cluster, AnyPos_T1)
		tmp						<- unique(subset(clumsm.info, select=Patient))
		set(tmp, NULL, 'PatientA', paste('P',seq_len(nrow(tmp)),sep=''))
		clumsm.info				<- merge(clumsm.info, tmp, by='Patient')
		outfile					<- paste(indir,'/',infile, '_', clu.infilexml.template, '_', clu.infilexml.opt, '_', gsub('/',':',insignat), '_', 'pt_anypos_3.5_anynodectime', '.pdf',sep='')	
		project.athena.Fisheretal.plot.selected.transmitters(clumsm.info, df.immu, df.viro, df.treatment, df.tpairs, cluphy, cluphy.info, cluphy.subtrees, cluphy.map.nodectime, outfile, pdf.height=900)		
	}
	#
	#	get timelines for the candidate transmitters in ATHENA.clu to the recently infected RI.PT; remove zero scores
	#
	resume			<- 1	
	save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RICT',method.PDT,'_',method,'_tATHENAclu','.R',sep='')	
	YX.part1		<- project.athena.Fisheretal.YX.part1(df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=NULL, df.tpairs=df.tpairs, tperiod.info=NULL, lRNA.supp=method.lRNA.supp, method.minLowerUWithNegT=method.minLowerUWithNegT, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)
	YX.part1		<- merge( YX.part1, subset( df.tpairs, select=c(FASTASampleCode, t.FASTASampleCode, cluster) ), by=c('FASTASampleCode','t.FASTASampleCode'), all.x=1)
	YX.part1[, class:='pt']
	gc()	
	rm.zero.score	<- TRUE	
	save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'YX',method.PDT,method,'.R',sep='')
	save.all		<- FALSE
	#save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'YX',method.PDT,method,'_all.R',sep='')
	#save.all		<- TRUE
	# df.tpairs.4.rawbrl=df.tpairs; thresh.pcoal=method.thresh.pcoal; brl.bwhost.multiplier=method.brl.bwhost; method.minLowerUWithNegT=method.minLowerUWithNegT; lRNA.supp=method.lRNA.supp; infilecov= infile.cov.study
	YX				<- project.athena.Fisheretal.YX.part2(	YX.part1, df.all, df.treatment, df.viro, predict.t2inf, t2inf.args, indir, insignat, indircov, infile.cov.study, infiletree, infile.trm.model, outdir, outfile, cluphy=cluphy, cluphy.info=cluphy.info, cluphy.map.nodectime=cluphy.map.nodectime, df.tpairs.4.rawbrl=df.tpairs, dur.Acute=dur.Acute,
															rm.zero.score=rm.zero.score, any.pos.grace.yr=any.pos.grace.yr, thresh.pcoal=method.thresh.pcoal, brl.bwhost.multiplier=method.brl.bwhost, method.minLowerUWithNegT=method.minLowerUWithNegT, lRNA.supp=method.lRNA.supp,
															t.period=t.period, save.file=save.file, resume=resume, method=method, save.all=save.all)
	gc()
	tperiod.info	<- merge(df.all, unique( subset(YX, select=c(Patient, t.period)) ), by='Patient')
	tperiod.info	<- tperiod.info[, list(t.period.min=min(AnyPos_T1)), by='t.period']
	set(tperiod.info, 1L,'t.period.min',t.recent.startctime-0.1)	
	tperiod.info[, t.period.max:=c(tperiod.info[-1, t.period.min], t.recent.endctime)]
	ri.PT			<- subset(YX, score.Y>0, select=c(Patient, t))[, list(n.t.infw= length(unique(t))), by='Patient']	
	#
	#	Sensitivity analysis:
	#	get timelines for all patients in ATHENA.seq to the recently infected RI.PT; keep zero scores
	#	randomly sample for every RI a number of t.Patients -- otherwise this results in 60e6 rows which we just cannot handle
	#	
	if(0)
	{
		rm.zero.score	<- FALSE
		tmp				<- unique(subset(df.tpairs, select=c(Patient, FASTASampleCode)))
		sample.n		<- subset(YX, select=c(Patient,t.Patient))
		setkey(sample.n, Patient, t.Patient)	
		sample.n		<- nrow( unique(sample.n) )
		sample.exclude	<- unique(subset(YX, select=t.Patient))
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RICT_',method,'_tATHENAseq','.R',sep='')
		YXS.part1		<- project.athena.Fisheretal.YX.part1(df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=tmp, df.tpairs=NULL, tperiod.info=NULL, lRNA.supp=method.lRNA.supp, sample.n=sample.n, sample.exclude=sample.exclude, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)
		if( !all(YX.part1[,FASTASampleCode]%in%df.all[,FASTASampleCode]) | !all(YX.part1[,t.FASTASampleCode]%in%df.all[,FASTASampleCode]) | !all(YXS.part1[,FASTASampleCode]%in%df.all[,FASTASampleCode]) | !all(YXS.part1[,t.FASTASampleCode]%in%df.all[,FASTASampleCode]) ) stop('unexpected FASTASampleCode')		
		YX.part1[, class:='pt']
		YXS.part1[, class:='pn']
		YXS.part1		<- do.call('rbind',list(YX.part1, subset(YXS.part1, select=colnames(YX.part1))))		#put potential transmitters and potential non-transmitters together		
		YXS.part1		<- merge( YXS.part1, subset( df.tpairs, select=c(FASTASampleCode, t.FASTASampleCode, cluster) ), by=c('FASTASampleCode','t.FASTASampleCode'), all.x=1)	
		rm.zero.score	<- FALSE
		thresh.pcoal	<- 0.75
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'YX',method,'_RICT_tATHENAseq','.R',sep='')
		YXS				<- project.athena.Fisheretal.YX.part2(	YXS.part1, df.all, df.treatment, df.viro, predict.t2inf, t2inf.args, indir, insignat, indircov, infile.cov.study, infiletree, infile.trm.model, outdir, paste(outfile, 'RICT_tATHENAseq', sep='_'), cluphy=cluphy, cluphy.info=cluphy.info, cluphy.map.nodectime=cluphy.map.nodectime, df.tpairs.4.rawbrl=df.tpairs, dur.Acute=dur.Acute, 
																rm.zero.score=rm.zero.score, any.pos.grace.yr=any.pos.grace.yr, thresh.pcoal=thresh.pcoal, lRNA.supp=method.lRNA.supp, t.period=t.period, save.file=save.file, resume=resume, method=method)		
	}
	#
	#	get timelines for all clustering candidate transmitters to the recently infected RI.PT
	#	
	X.clu				<- NULL
	if(is.null(X.tables) && grepl('.clu',method.risk,fixed=1))
	{		
		tmp				<- copy(clumsm.info)
		setkey(tmp, Patient)
		tmp				<- unique(tmp)	
		tmp[, FASTASampleCode:=NULL]
		tmp[, cluster:=NULL]
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RICLU_',method,'_tATHENAclu','.R',sep='')
		X.clu			<- project.athena.Fisheretal.YX.part1(	tmp, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=ri.CLU, df.tpairs=NULL, 
																tperiod.info=tperiod.info, lRNA.supp=method.lRNA.supp, t.period=t.period, t.endctime=t.endctime, method.minLowerUWithNegT=method.minLowerUWithNegT, save.file=save.file, resume=resume)			
		gc()		
	}
	#
	#	get timelines for all candidate transmitters in df.all (anyone with subtype B sequ) to the recently infected RI.PT
	#	
	X.seq			<- NULL
	if(is.null(X.tables))
	{
		if(method.PDT=='CLU')
		{			
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RICLU_',method,'_tATHENAseq','.R',sep='')
			X.seq			<- project.athena.Fisheretal.YX.part1(	df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=ri.CLU, df.tpairs=NULL, 
																	tperiod.info=tperiod.info, lRNA.supp=method.lRNA.supp, t.period=t.period, t.endctime=t.endctime, method.minLowerUWithNegT=method.minLowerUWithNegT, save.file=save.file, resume=resume)
		}
		if(method.PDT=='SEQ')
		{						
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RISEQ_',method,'_tATHENAseq','.R',sep='')
			X.seq			<- project.athena.Fisheretal.YX.part1(	df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=ri.SEQ, df.tpairs=NULL, 
																	tperiod.info=tperiod.info, lRNA.supp=method.lRNA.supp, t.period=t.period, t.endctime=t.endctime, method.minLowerUWithNegT=method.minLowerUWithNegT, save.file=save.file, resume=resume)
		}
		gc()				
	}
	#
	#	get timelines for all candidate transmitters in ATHENA.MSM (anyone with MSM exposure group irrespective of sequ available or not) to the recently infected RI.PT
	#	
	X.msm				<- NULL	
	if(is.null(X.tables) & ( grepl('adj',method.risk) | grepl('cens',method.risk)))
	{
		if(method.PDT=='CLU')
		{						
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RICLU_',method,'_tATHENAmsm','.R',sep='')
			X.msm			<- project.athena.Fisheretal.YX.part1(	df.all.allmsm, df.immu.allmsm, df.viro.allmsm, df.treatment.allmsm, predict.t2inf, t2inf.args, ri=ri.CLU, df.tpairs=NULL, 
																	tperiod.info=tperiod.info, lRNA.supp=method.lRNA.supp, t.period=t.period, t.endctime=t.endctime, method.minLowerUWithNegT=method.minLowerUWithNegT, save.file=save.file, resume=resume)
		}
		if(method.PDT=='SEQ')
		{												
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RISEQ_',method,'_tATHENAmsm','.R',sep='')
			X.msm			<- project.athena.Fisheretal.YX.part1(	df.all.allmsm, df.immu.allmsm, df.viro.allmsm, df.treatment.allmsm, predict.t2inf, t2inf.args, ri=ri.SEQ, df.tpairs=NULL, 
																	tperiod.info=tperiod.info, lRNA.supp=method.lRNA.supp, t.period=t.period, t.endctime=t.endctime, method.minLowerUWithNegT=method.minLowerUWithNegT, save.file=save.file, resume=resume)
		}
		gc()		
	}
	#
	#X.clu	<- X.clu[sample(seq_len(nrow(X.clu)), 1e6),]
	#X.seq	<- X.seq[sample(seq_len(nrow(X.seq)), 2e6),]
	#X.msm	<- X.msm[sample(seq_len(nrow(X.msm)), 3e6),]
	#
	#	characteristics of RI(with potential direct transmitter) and RI(in all.msm)
	#	
	if(0)
	{
		ri.PT.su				<- project.athena.Fisheretal.RI.summarize(ri.PT, df.all, tperiod.info, info=clumsm.info, with.seq=TRUE, with.cluster=TRUE, cols=c('lRNA_T1','CD4_T1','AnyPos_A','SC_dt','PosCD4_dt','PoslRNA_dt'))
		#unique(subset(ri.ALLMSM, select=Patient)) 
		ri.allmsm.su			<- project.athena.Fisheretal.RI.summarize(unique(subset(ri.ALLMSM, select=Patient)), df.all.allmsm, tperiod.info, info=NULL, with.seq=FALSE, with.cluster=FALSE, cols=c('lRNA_T1','CD4_T1','AnyPos_A','SC_dt','PosCD4_dt','PoslRNA_dt'))		
		ri.allmsm.resolution	<- project.athena.Fisheretal.CT.resolution(df.all.allmsm, df.viro.allmsm, df.immu.allmsm, t.startctime=1994, t.endctime=t.endctime)
		ri.allmsm.l2c			<- project.athena.Fisheretal.CT.link2care(df.all.allmsm, link2care.cut=c(-Inf,15/365,30/365,3/12,Inf), link2care.label=c('<=15d','<=30d','<=3m','>3m'), t.startctime=1994)
		save.file				<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RIPDT_RIMSM_',method,'_info','.R',sep='')
		save(ri.PT.su, ri.allmsm.su, ri.allmsm.resolution, ri.allmsm.l2c, file=save.file)	
	}	
	#
	#	stratify
	#	
	if(!is.null(X.tables))
	{
		if(any(sapply(c('m2wmx','m2Bwmx','m2Cwmx'), grepl, x=method.risk)))
		{
			save.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'YXm2_',method.PDT,'_',method,'.R',sep='')
			YX					<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(YX, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp, plot.file.or=NA, resume=1, save.file=save.file ) 
		}						
		if(any(sapply(c('m21st','m2B1st'), grepl, x=method.risk)))
			YX					<- project.athena.Fisheretal.YX.model2.stratify.VL1stsu(YX, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp)
		if(any(sapply(c('m2t','m2Bt'), grepl, x=method.risk)))
			YX					<- project.athena.Fisheretal.YX.model2.stratify.VLt(YX, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp)
		if(any(sapply(c('m2gm','m2Bgm'), grepl, x=method.risk)))
			YX					<- project.athena.Fisheretal.YX.model2.stratify.VLgm(YX, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp)
		if(substr(method.risk,1,2)=='m3')			
			YX					<- project.athena.Fisheretal.YX.model3.stratify.ARTriskgroups(YX, df.all, return.only.ART=TRUE)
		if(substr(method.risk,1,2)=='m4')			
			YX					<- project.athena.Fisheretal.YX.model4.stratify.Diagnosed(YX, df.immu, return.only.Diag=TRUE )
		if(substr(method.risk,1,2)=='m5')			
			YX					<- project.athena.Fisheretal.YX.model5.stratify(YX)
	}
	if(is.null(X.tables))
	{		
		cat(paste('\nstratify by', method))
		if(substr(method.risk,1,2)=='m2')
		{
			if(any(sapply(c('m2wmx','m2Bwmx','m2Cwmx'), grepl, x=method.risk)))
			{
				save.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'YXm2_',method.PDT,'_',method,'.R',sep='')
				YX					<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(YX, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp, plot.file.or=NA, resume=1, save.file=save.file )
				YX					<- subset(YX, select=c(t, Patient, t.Patient, t.period, AnyPos_T1, t.AnyPos_T1, CD4b, CD4c, CD4b.tperiod, CD4c.tperiod))
				gc()
				save.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Xclum2_',method.PDT,'_',method,'.R',sep='')
				X.clu				<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.clu, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp, plot.file.or=NA, resume=1, save.file=save.file )
				X.clu				<- subset(X.clu, select=c(t, Patient, t.Patient, t.period, AnyPos_T1, t.AnyPos_T1, CD4b, CD4c, CD4b.tperiod, CD4c.tperiod))
				gc()
				save.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Xseqm2_',method.PDT,'_',method,'.R',sep='')
				X.seq				<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.seq, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp, plot.file.or=NA, resume=1, save.file=save.file )
				X.seq				<- subset(X.seq, select=c(t, Patient, t.Patient, t.period, AnyPos_T1, t.AnyPos_T1, CD4b, CD4c, CD4b.tperiod, CD4c.tperiod))
				gc()
				save.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Xmsmm2_',method.PDT,'_',method,'.R',sep='')
				X.msm				<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.msm, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, lRNA.supp=method.lRNA.supp, plot.file.or=NA, resume=1, save.file=save.file )
				X.msm				<- subset(X.msm, select=c(t, Patient, t.Patient, t.period, AnyPos_T1, t.AnyPos_T1, CD4b, CD4c, CD4b.tperiod, CD4c.tperiod))
				gc()
			}
			if(any(sapply(c('m21st','m2B1st'), grepl, x=method.risk)))
			{
				YX					<- project.athena.Fisheretal.YX.model2.stratify.VL1stsu(YX, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp)
				X.clu				<- project.athena.Fisheretal.YX.model2.stratify.VL1stsu(X.clu, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp)
				X.seq				<- project.athena.Fisheretal.YX.model2.stratify.VL1stsu(X.seq, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp)
				X.msm				<- project.athena.Fisheretal.YX.model2.stratify.VL1stsu(X.msm, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, lRNA.supp=method.lRNA.supp)
			}
			if(any(sapply(c('m2t','m2Bt'), grepl, x=method.risk)))
			{
				YX					<- project.athena.Fisheretal.YX.model2.stratify.VLt(YX, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp)
				X.clu				<- project.athena.Fisheretal.YX.model2.stratify.VLt(X.clu, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp)
				X.seq				<- project.athena.Fisheretal.YX.model2.stratify.VLt(X.seq, df.all, df.viro, df.immu, lRNA.supp=method.lRNA.supp)
				X.msm				<- project.athena.Fisheretal.YX.model2.stratify.VLt(X.msm, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, lRNA.supp=method.lRNA.supp)
			}
		}				
		if(substr(method.risk,1,2)=='m3')
		{						
			return.only.ART	<- FALSE		#need all including undiagnosed for cens.tables
			#return.only.ART	<- TRUE			
			YX				<- project.athena.Fisheretal.YX.model3.stratify.ARTriskgroups(YX, df.all, return.only.ART=return.only.ART)	
			X.clu			<- project.athena.Fisheretal.YX.model3.stratify.ARTriskgroups(X.clu, df.all, return.only.ART=return.only.ART)
			X.seq			<- project.athena.Fisheretal.YX.model3.stratify.ARTriskgroups(X.seq, df.all, return.only.ART=return.only.ART)
			X.msm			<- project.athena.Fisheretal.YX.model3.stratify.ARTriskgroups(X.msm, df.all.allmsm, return.only.ART=return.only.ART)
		}
		if(substr(method.risk,1,2)=='m4')			
		{
			YX				<- project.athena.Fisheretal.YX.model4.stratify.Diagnosed(YX, df.immu, return.only.Diag=FALSE )
			X.clu			<- project.athena.Fisheretal.YX.model4.stratify.Diagnosed(X.clu, df.immu, return.only.Diag=FALSE )
			X.seq			<- project.athena.Fisheretal.YX.model4.stratify.Diagnosed(X.seq, df.immu, return.only.Diag=FALSE )
			X.msm			<- project.athena.Fisheretal.YX.model4.stratify.Diagnosed(X.msm, df.immu, return.only.Diag=FALSE )
		}
		if(substr(method.risk,1,2)=='m5')			
		{
			YX				<- project.athena.Fisheretal.YX.model5.stratify(YX)
			X.clu			<- project.athena.Fisheretal.YX.model5.stratify(X.clu)
			X.seq			<- project.athena.Fisheretal.YX.model5.stratify(X.seq)
			X.msm			<- project.athena.Fisheretal.YX.model5.stratify(X.msm)
		}					
		stop()
		#	compute tables
		if(grepl('adj',method.risk) & grepl('clu',method.risk))
		{
			save.file		<- NA
			resume			<- 0
			if(grepl('m21st',method.risk))		save.file	<- 'm21st'
			if(grepl('m2B1st',method.risk))		save.file	<- 'm2B1st'
			if(grepl('m2t',method.risk))		save.file	<- 'm2t'
			if(grepl('m2Bt',method.risk))		save.file	<- 'm2Bt'
			if(grepl('m2wmx',method.risk))		save.file	<- 'm2wmx'
			if(grepl('m2Bwmx',method.risk))		save.file	<- 'm2Bwmx'
			if(grepl('m2Cwmx',method.risk))		save.file	<- 'm2Cwmx'
			if(grepl('m3.n3mx',method.risk) & !grepl('No',method.risk))								save.file	<- 'm3.n3mx'
			if(grepl('m3.ind',method.risk) & !grepl('No',method.risk))								save.file	<- 'm3.ind'
			if(grepl('m3.ind',method.risk) & grepl('No',method.risk))								save.file	<- 'm3.indNo'
			if(grepl('m3.indmx',method.risk) & !grepl('No',method.risk))							save.file	<- 'm3.indmx'
			if(grepl('m3.indmx',method.risk) & grepl('No',method.risk))								save.file	<- 'm3.indmxNo'
			if(grepl('m3.nnrtpiNo',method.risk))													save.file	<- 'm3.nnrtpiNo'
			if(grepl('m4.Bwmx',method.risk))	save.file	<- 'm4.Bwmx'
			if(grepl('m5.tA',method.risk))		save.file	<- 'm5.tA'
			if(grepl('m5.tAb',method.risk))		save.file	<- 'm5.tAb'
			if(grepl('m5.tAc',method.risk))		save.file	<- 'm5.tAc'
			if(grepl('m5.tiA',method.risk))		save.file	<- 'm5.tiA'
			if(grepl('m5.tiAb',method.risk))	save.file	<- 'm5.tiAb'
			if(grepl('m5.tiAc',method.risk))	save.file	<- 'm5.tiAc'
			if(is.na(save.file))	stop('unknown method.risk')				
			tmp				<- regmatches(method.risk, regexpr('tp[0-9]', method.risk))			
			save.file		<- paste(save.file, ifelse(length(tmp), paste('.',tmp,sep=''), ''), sep='')
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_tables',method.PDT,'_',save.file,'.R',sep='')
			X.tables		<- project.athena.Fisheretal.estimate.risk.table(YX, X.seq, X.msm, X.clu, tperiod.info=tperiod.info, resume=resume, save.file=save.file, method=method.risk)
			stop()
		}
		X.clu<- X.seq<- X.msm<- NULL
		gc()
	}
	stopifnot(is.null(X.tables)==FALSE)
	#
	#	endpoint: first VL suppressed
	#
	#X.seq<- X.seq[sample(1:nrow(X.seq),2e6),]
	resume				<- 1
	bs.n				<- 1e3
	if(method.PDT=='')
		method.PDT<- 'PDT'
	if(grepl('m2',method.risk))
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_denom',method.PDT,'_model2_',method.risk,'.R',sep='')
	if(grepl('m3',method.risk))
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_denom',method.PDT,'_model3_',method.risk,'.R',sep='')	
	if(grepl('m4',method.risk))
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_denom',method.PDT,'_model4_',method.risk,'.R',sep='')		
	if(grepl('m5',method.risk))
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_denom',method.PDT,'_model5_',method.risk,'.R',sep='')			
	tmp					<- project.athena.Fisheretal.estimate.risk.wrap(YX, X.tables, tperiod.info, plot.file.or=NA, bs.n=1e3, resume=resume, save.file=save.file, method.risk=method.risk)					
	
	stop()
		
		
	tmp		<- subset(YX, select=c(Patient, t.Patient, score.Y))	
	setkey(tmp, Patient, t.Patient)
	tmp		<- unique(tmp)
	YX.rawp	<- tmp[ , list(t.Patient=t.Patient, P=score.Y/sum(score.Y), Y=score.Y), by='Patient']
	#
	#	select ART + suppressed transmitters in last time period; based on observed counts alone
	#
	set(YX, NULL, 'score.Y', YX[, score.Y*w.tn])
	set(YX, NULL, c('w','w.tn'), 1.)
	YXtp					<- subset(YX, t.period==4)
	YXtp					<- merge(YXtp, YXtp[, list(t.Patient=t.Patient, t=t, y=score.Y/sum(score.Y)), by='Patient'], by=c('Patient','t.Patient','t'))
	tmp						<- YXtp[, list(yf= sum(y)), by=c('Patient','CD4c.tperiod')]
	tmp						<- subset(tmp, CD4c.tperiod=='ART.suA.Y2.4' & yf>0.5)
	df.tpairs				<- merge( df.tpairs, tmp, by='Patient' )
	tmp						<- project.athena.Fisheretal.get.dated.phylo.for.selection(df.tpairs, clu.indir, clu.infile, clu.insignat, clu.infilexml.opt, clu.infilexml.template, method.nodectime=method.nodectime)
	df.tpairs				<- tmp$df.tpairs
	cluphy					<- tmp$clu$cluphy
	cluphy.subtrees			<- tmp$clu$cluphy.subtrees
	cluphy.info				<- tmp$clu$cluphy.info
	cluphy.map.nodectime	<- tmp$clu$cluphy.map.nodectime
	outfile					<- paste(indir,'/',infile, '_', clu.infilexml.template, '_', clu.infilexml.opt, '_', gsub('/',':',insignat), '_', 'SELECT_ARTSUPPOBS', '.pdf',sep='')	
	project.athena.Fisheretal.plot.selected.transmitters(clumsm.info, df.immu, df.viro, df.treatment, df.tpairs, cluphy, cluphy.info, cluphy.subtrees, cluphy.map.nodectime, outfile, pdf.height=50, label.select=c("cluster","Patient","Trm","isAcute","RegionHospital"))
	#
	tmp	<- subset(df.tpairs, select=c(Patient, t.Patient, yf, cluster))
	setkey(tmp, Patient, t.Patient)
	tmp	<- unique(tmp)
	tmp	<- merge(YXtp, tmp, by=c('t.Patient','Patient'))
	#
	#	get data for selection
	#		
	df.tpairs				<- subset(df.tpairs, cluster%in%c(1512, 1510, 1508))
	tmp						<- project.athena.Fisheretal.get.dated.phylo.for.selection(df.tpairs, clu.indir, clu.infile, clu.insignat, clu.infilexml.opt, clu.infilexml.template, method.nodectime=method.nodectime)
	#df.all2					<- subset(clumsm.info, select=c(Patient, cluster))
	#setkey(df.all2, Patient)
	#df.all2					<- merge( tmp$df.all, unique(df.all2), by='Patient', all.x=TRUE )
	df.tpairs				<- tmp$df.tpairs
	cluphy					<- tmp$clu$cluphy
	cluphy.subtrees			<- tmp$clu$cluphy.subtrees
	cluphy.info				<- tmp$clu$cluphy.info
	cluphy.map.nodectime	<- tmp$clu$cluphy.map.nodectime
	#	anonymize
	setkey(clumsm.info, cluster, AnyPos_T1)
	tmp						<- unique(subset(clumsm.info, select=Patient))
	#set(tmp, NULL, 'PatientA', paste('P',seq_len(nrow(tmp)),sep=''))
	set(tmp, NULL, 'PatientA', tmp[, Patient])
	clumsm.info				<- merge(clumsm.info, tmp, by='Patient')
	#			
	#	plot	
	#	
	outfile					<- paste(indir,'/',infile, '_', clu.infilexml.template, '_', clu.infilexml.opt, '_', gsub('/',':',insignat), '_', 'pt_anypos_3.5_anynodectime', '.pdf',sep='')	
	project.athena.Fisheretal.plot.selected.transmitters(clumsm.info, df.immu, df.viro, df.treatment, df.tpairs, cluphy, cluphy.info, cluphy.subtrees, cluphy.map.nodectime, outfile, pdf.height=600)

	
}

