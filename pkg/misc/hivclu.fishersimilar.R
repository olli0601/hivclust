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
	set(contact, contact[,which(is.na(DateDied))], 'DateDied', t.endctime)
	set(contact, contact[, which(DateLastContact==DateDied)], 'DateLastContact', NA_real_)
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
project.athena.Fisheretal.v1.YX<- function(df.all, clumsm.info, df.tpairs, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, cluphy, cluphy.info, cluphy.map.nodectime, insignat, indircov, infilecov, infiletree, outdir, outfile, t.period=0.25, t.endctime=2013., save.file=NA)
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
		#	get an idea how many time periods per year we could have:
		#	tmp<- project.athena.Fisheretal.median.followup(df.tpairs, df.viro, df.immu, X.tperiod)
		#	compute X: follow up
		X.fwup					<- project.athena.Fisheretal.X.followup(df.tpairs, df.immu, t.period=t.period, t.endctime=t.endctime)
		#	compute X: InCare
		X.incare				<- project.athena.Fisheretal.X.incare(df.tpairs, clumsm.info, df.viro, df.immu, df.treatment, t.period=t.period, t.endctime=t.endctime)
		X.incare				<- project.athena.Fisheretal.X.nocontact(X.incare, df.tpairs, clumsm.info, t.period=t.period, t.endctime=t.endctime)
		X.incare				<- project.athena.Fisheretal.X.CDCC(X.incare, df.tpairs, clumsm.info, t.period=t.period, t.endctime=t.endctime)	
		#	compute X: before care
		X.b4care				<- project.athena.Fisheretal.X.b4care(df.tpairs, clumsm.info, predict.t2inf, t2inf.args, t.period=t.period)
		tmp						<- merge(X.incare, X.b4care, by=c('t.Patient','t'))
		cat(paste('\nnumber entries (Patient,t) that overlap in and before care [should be zero], n=',nrow(tmp)))
		X.pt					<- merge(X.incare, X.b4care, by=c('t.Patient','t'), all.x=1, all.y=1)
		set(X.pt, X.pt[,which(is.na(stage))], 'stage', 'U')
		X.pt					<- project.athena.Fisheretal.X.age(X.pt, clumsm.info, t.period=t.period)	
		X.pt					<- merge(X.pt, X.fwup, by=c('t.Patient', 't'), all.x=1)
		#	complete of transmitter: AnyPos_T1, AnyT_T1, AnyPos_a, isAcute
		X.pt[, AnyPos_T1:=NULL]
		X.pt[, AnyT_T1:=NULL]
		X.pt[, AnyPos_a:=NULL]
		X.pt[, isAcute:=NULL]
		tmp						<- unique(subset( clumsm.info, select=c(Patient, DateBorn, AnyPos_T1, AnyT_T1, isAcute) ))
		set(tmp, NULL, 'DateBorn', tmp[, AnyPos_T1-DateBorn])
		setnames( tmp, 'DateBorn', 'AnyPos_a')
		setnames( tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		X.pt					<- merge(X.pt, tmp, by='t.Patient', allow.cartesian=TRUE)
		
		#tmp					<- project.athena.Fisheretal.X.followup.compareCD4toVL(df.tpairs, df.immu, clumsm.info)
		#	compute X: potential transmitters with pulsed ART shortly after seroconversion / diagnosis
		X.ARTpulsed				<- project.athena.Fisheretal.X.ART.pulsed(df.tpairs, clumsm.info, t.pulse.ART=0.75, t.pulse.sc=NA, t.pulse.ART.I=1, t.pulse.ART.I.d=1)
		#	compute X: calendar time period
		X.tperiod				<- project.athena.Fisheretal.X.calendarperiod(df.tpairs, clumsm.info, t.period=t.period, c.nperiod= 4)
		#	compute X: RegionHospital and Exposure group
		X.Trm					<- project.athena.Fisheretal.X.Trm.Region(df.tpairs, clumsm.info)
		#	project.athena.Fisheretal.X.incare.check(X.incare)
		#	compute X: time to first viral suppression from diagnosis	
		X.t2.vlsupp				<- project.athena.Fisheretal.X.time.diag2suppressed(df.tpairs, clumsm.info, df.viro, lRNA.suppressed= log10(1e3), t2.vl.supp.p=c(0.1, 0.25))
		#	compute X: time to first VL and first CD4 measurement
		X.t2.care				<- project.athena.Fisheretal.X.time.diag2firstVLandCD4(df.tpairs, clumsm.info, df.viro, df.immu, t2.care.t1.q=c(0.25,0.5))
		#
		#	merge all static X that are independent of time periods
		#
		tmp						<- merge( df.tpairs, X.tperiod, by=c('Patient','t.Patient'), all.x=1 )
		tmp						<- merge( tmp, X.Trm, by='t.Patient', all.x=1 )
		tmp						<- merge( tmp, X.t2.care, by='t.Patient', all.x=1 )
		tmp						<- merge( tmp, X.t2.vlsupp, by='t.Patient', all.x=1 )
		tmp						<- merge( tmp, X.ARTpulsed, by='t.Patient', all.x=1 )	
		X.pt					<- merge(tmp, X.pt, by='t.Patient', allow.cartesian=TRUE)	
		#
		#	compute Y scores
		#
		#	BRL	[0,1]: raw branch length between pot transmitter and infected
		file					<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'nbrlraw',sep='')	
		tmp						<- project.athena.Fisheretal.Y.rawbrl(df.tpairs, indir, insignat, indircov, infilecov, infiletree, save.file=paste(file, '.R', sep=''), plot.file=paste(file, '.pdf', sep=''))
		Y.rawbrl				<- tmp$tpairs
		Y.rawbrl.linked			<- tmp$linked
		Y.rawbrl.unlinked		<- tmp$unlinked
		#	BRL [0,1]: branch length weight between pot transmitter and infected
		file					<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl.pdf',sep='')
		Y.brl					<- project.athena.Fisheretal.v1.Y.brlweight(Y.rawbrl, Y.rawbrl.linked	, Y.rawbrl.unlinked, df.all, linked.brl.min= -4, plot.file=file)
		cat(paste('\n number of transmitter-infected pairs with brl score, n=',nrow(Y.brl)))
		cat(paste('\n number of transmitter-infected pairs with zero brl score, n=',nrow(subset(Y.brl, score.brl==0))))
		#
		file					<- NA #paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'coalraw.R',sep='')	
		tmp						<- project.athena.Fisheretal.v1.Y.coal.and.inf(df.tpairs, clumsm.info, cluphy, cluphy.info, cluphy.map.nodectime, t.period=t.period, save.file=file )
		Y.inf					<- tmp$inf
		Y.coal					<- tmp$coal	
		cat(paste('\n number of transmitter-infected pairs, n=',nrow(Y.coal)))
		cat(paste('\n number of transmitter-infected pairs with coal !NA score, n=',nrow(subset(Y.coal, !is.na(score.coal)))))
		cat(paste('\n number of transmitter-infected pairs with zero coal score, n=',nrow(subset(Y.coal, !is.na(score.coal) & score.coal<=EPS))))
		cat(paste('\n number of transmitter-infected-t triplets with inf !NA score, n=',nrow(Y.inf)))
		cat(paste('\n number of transmitter-infected-t triplets with zero inf score, n=',nrow(subset(Y.inf, score.inf<=EPS))))
		#	merge and combine scores
		Y.score					<- merge( df.tpairs, subset( Y.coal, is.na(score.coal) | score.coal>0 ), by=c('FASTASampleCode', 't.FASTASampleCode') )
		Y.score					<- merge( Y.score, subset(Y.brl, is.na(score.brl) | score.brl>0), by=c('FASTASampleCode', 't.FASTASampleCode') )			
		#	before we merge over time periods t, keep only one (Patient, t.Patient) pair if there are multiple sequences per individual. Take the one with largest brl score
		tmp						<- Y.score[,	{	
						z<- which.max(score.brl)
						list(FASTASampleCode=FASTASampleCode[z], t.FASTASampleCode=t.FASTASampleCode[z])
					},by=c('Patient', 't.Patient')]	
		Y.score					<- merge( Y.score, subset(tmp, select=c(FASTASampleCode=FASTASampleCode, t.FASTASampleCode)), by=c('FASTASampleCode', 't.FASTASampleCode') )
		Y.score					<- merge( subset(Y.inf, is.na(score.inf) | score.inf>EPS), Y.score, by=c('FASTASampleCode', 't.FASTASampleCode') )
		#	merge over time periods t		
		#	check that each (t, Patient, t.Patient) triplet is only once in the data.table
		tmp	<- Y.score[, list(n= length(FASTASampleCode)), by=c('t','Patient', 't.Patient')]
		if(tmp[,any(n!=1)])	stop('unexpected multiple entries per t, Patient, t.Patient')
		#
		tmp						<- Y.score[, mean(score.coal,na.rm=TRUE) ]
		cat(paste('\nmean coalescence score is= ',tmp) )
		if( Y.score[, any(is.na(score.brl))] )	stop('unexpected missing score.brl')
		if( Y.score[, any(is.na(score.inf))] )	stop('unexpected missing score.inf')	
		tmp						<- Y.score[, list(score.Y= prod(c(score.inf, ifelse(is.na(score.coal), tmp, score.coal), score.brl))), by=c('FASTASampleCode', 't.FASTASampleCode', 't')]
		Y.score					<- merge(Y.score, tmp, by=c('FASTASampleCode', 't.FASTASampleCode', 't'))
		set(Y.score, NULL, 'score.Y', Y.score[, round(score.Y, d=5)])
		cat(paste('\nnumber of transmitter-infected-t triplets with score <=1e-3, n=',nrow(subset(Y.score, score.Y<=1e-3))))
		Y.score					<- subset(Y.score, score.Y>1e-3)
		cat(paste('\nnumber of transmitter-infected-t triplets with score >1e-3, n=',nrow(Y.score)))
		cat(paste('\npercentage of transmitter-infected-t triplets with score >0.5, n=',Y.score[, mean(score.Y>0.5)]))
		cat(paste('\npercentage of recently infected for which there is at least one potential transmitter with score >0.5, n=',Y.score[,  list(selected= any(score.Y>0.5)) ,by='Patient'][, mean(selected)]))
		#
		#	plot the correlation in the score. 
		#
		hist( Y.score[, score.Y], breaks=100 )
		plot( Y.score[, score.brl], Y.score[, score.inf], pch=18)
		plot( Y.score[, score.brl], Y.score[, score.coal], pch=18)
		plot( Y.score[, score.inf], Y.score[, score.coal], pch=18)
		#
		#	checks
		#
		check	<-	subset(Y.score, select=c(t.Patient, Patient))
		setkey(check, t.Patient, Patient)
		check	<- unique(check)
		tmp		<- subset(X.pt, select=c(t.Patient, Patient))
		setkey(tmp, t.Patient, Patient)
		tmp		<- unique(tmp)[check]
		if(nrow(tmp)!=nrow(check))	stop('unexpected potential transmission pairs in Y.score that are not in X.pt')
		cat(paste('\nnumber of transmitter-infected pairs with non-zero score, n=',nrow(check)))
		cat(paste('\nnumber of recently infected patients in transmitter-infected pairs with non-zero score, n=',length(check[, unique(Patient)])))
		cat(paste('\nnumber of transmitters in transmitter-infected pairs with non-zero score, n=',length(check[, unique(t.Patient)])))
		#	
		YX		<- merge( subset(Y.score, select=c(t, t.FASTASampleCode, FASTASampleCode, score.inf, score.coal, score.brl, score.Y)), X.pt, by= c('FASTASampleCode', 't.FASTASampleCode', 't'))
		#	re-arrange a little
		YX		<- subset(YX, select=	c(	t, t.Patient, Patient, cluster, score.Y, 																	#triplets identifiers and Y score
						t.period, stage, U.score, contact, CDCC, lRNA, CD4, 														#main covariates							
						t.Age, t.AnyPos_a, RegionHospital, ART.I, ART.F, ART.A, ART.P, ART.pulse, fw.up, t2.care.t1, t2.vl.supp, 		#secondary covariates
						t.AnyPos_T1,  t.AnyT_T1, StartTime, StopTime, lRNAc, t.isAcute, Trm,												#other 
						score.inf, score.coal, score.brl, FASTASampleCode, t.FASTASampleCode										#other							
				))
		setkey(YX, t, t.Patient, Patient)
		#
		if(!is.na(save.file))
		{
			cat(paste('\nsave YX to file', save.file))
			save(YX, file=save.file)
		}
	}
	YX	
}
######################################################################################
project.athena.Fisheretal.v2.YX<- function(df.all, clumsm.info, df.tpairs, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, cluphy, cluphy.info, cluphy.map.nodectime, insignat, indircov, infilecov, infiletree, outdir, outfile, t.period=0.25, t.endctime=2013., save.file=NA)
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
		#	get an idea how many time periods per year we could have:
		#	tmp<- project.athena.Fisheretal.median.followup(df.tpairs, df.viro, df.immu, X.tperiod)
		#	compute X: follow up
		X.fwup					<- project.athena.Fisheretal.X.followup(df.tpairs, df.immu, t.period=t.period, t.endctime=t.endctime)
		#	compute X: InCare
		X.incare				<- project.athena.Fisheretal.X.incare(df.tpairs, clumsm.info, df.viro, df.immu, df.treatment, t.period=t.period, t.endctime=t.endctime)
		X.incare				<- project.athena.Fisheretal.X.nocontact(X.incare, df.tpairs, clumsm.info, t.period=t.period, t.endctime=t.endctime)
		X.incare				<- project.athena.Fisheretal.X.CDCC(X.incare, df.tpairs, clumsm.info, t.period=t.period, t.endctime=t.endctime)	
		#	compute X: before care
		X.b4care				<- project.athena.Fisheretal.X.b4care(df.tpairs, clumsm.info, predict.t2inf, t2inf.args, t.period=t.period)
		tmp						<- merge(X.incare, X.b4care, by=c('t.Patient','t'))
		cat(paste('\nnumber entries (Patient,t) that overlap in and before care [should be zero], n=',nrow(tmp)))
		X.pt					<- merge(X.incare, X.b4care, by=c('t.Patient','t'), all.x=1, all.y=1)
		set(X.pt, X.pt[,which(is.na(stage))], 'stage', 'U')
		X.pt					<- project.athena.Fisheretal.X.age(X.pt, clumsm.info, t.period=t.period)	
		X.pt					<- merge(X.pt, X.fwup, by=c('t.Patient', 't'), all.x=1)
		#	complete of transmitter: AnyPos_T1, AnyT_T1, AnyPos_a, isAcute
		X.pt[, AnyPos_T1:=NULL]
		X.pt[, AnyT_T1:=NULL]
		X.pt[, AnyPos_a:=NULL]
		X.pt[, isAcute:=NULL]
		tmp						<- unique(subset( clumsm.info, select=c(Patient, DateBorn, AnyPos_T1, AnyT_T1, isAcute) ))
		set(tmp, NULL, 'DateBorn', tmp[, AnyPos_T1-DateBorn])
		setnames( tmp, 'DateBorn', 'AnyPos_a')
		setnames( tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		X.pt					<- merge(X.pt, tmp, by='t.Patient', allow.cartesian=TRUE)
		
		#tmp					<- project.athena.Fisheretal.X.followup.compareCD4toVL(df.tpairs, df.immu, clumsm.info)
		#	compute X: potential transmitters with pulsed ART shortly after seroconversion / diagnosis
		X.ARTpulsed				<- project.athena.Fisheretal.X.ART.pulsed(df.tpairs, clumsm.info, t.pulse.ART=0.75, t.pulse.sc=NA, t.pulse.ART.I=1, t.pulse.ART.I.d=1)
		#	compute X: calendar time period
		X.tperiod				<- project.athena.Fisheretal.X.calendarperiod(df.tpairs, clumsm.info, t.period=t.period, c.nperiod= 4)
		#	compute X: RegionHospital and Exposure group
		X.Trm					<- project.athena.Fisheretal.X.Trm.Region(df.tpairs, clumsm.info)
		#	project.athena.Fisheretal.X.incare.check(X.incare)
		#	compute X: time to first viral suppression from diagnosis	
		X.t2.vlsupp				<- project.athena.Fisheretal.X.time.diag2suppressed(df.tpairs, clumsm.info, df.viro, lRNA.suppressed= log10(1e3), t2.vl.supp.p=c(0.1, 0.25))
		#	compute X: time to first VL and first CD4 measurement
		X.t2.care				<- project.athena.Fisheretal.X.time.diag2firstVLandCD4(df.tpairs, clumsm.info, df.viro, df.immu, t2.care.t1.q=c(0.25,0.5))
		#
		#	merge all static X that are independent of time periods
		#
		tmp						<- merge( df.tpairs, X.tperiod, by=c('Patient','t.Patient'), all.x=1 )
		tmp						<- merge( tmp, X.Trm, by='t.Patient', all.x=1 )
		tmp						<- merge( tmp, X.t2.care, by='t.Patient', all.x=1 )
		tmp						<- merge( tmp, X.t2.vlsupp, by='t.Patient', all.x=1 )
		tmp						<- merge( tmp, X.ARTpulsed, by='t.Patient', all.x=1 )	
		X.pt					<- merge(tmp, X.pt, by='t.Patient', allow.cartesian=TRUE)
		#
		#	compute infection window of recipient
		#
		Y.infwindow				<- project.athena.Fisheretal.Y.infectiontime(df.tpairs, clumsm.info, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, score.min=0.1, score.set.value=1, method='for.infected')		
		#
		#	compute Y scores
		#
		#	BRL	[0,1]: raw branch length between pot transmitter and infected
		file					<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'nbrlraw2e',sep='')	
		tmp						<- project.athena.Fisheretal.Y.rawbrl(df.tpairs, indir, insignat, indircov, infilecov, infiletree, save.file=paste(file, '.R', sep=''), plot.file=paste(file, '.pdf', sep=''))
		Y.rawbrl				<- tmp$tpairs
		Y.rawbrl.linked			<- tmp$linked
		Y.rawbrl.unlinked		<- tmp$unlinked
		#	BRL [0,1]: branch length weight between pot transmitter and infected
		file					<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl.pdf',sep='')
		Y.brl					<- project.athena.Fisheretal.v2.Y.brlweight(Y.rawbrl, Y.rawbrl.linked	, Y.rawbrl.unlinked, df.all, brl.linked.min=-4, brl.linked.min.dt=1.5, plot.file=file)	
		#	COAL [0,1]: prob that coalescence is within the transmitter
		file					<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'coalraw2e.R',sep='')		
		Y.coal					<- project.athena.Fisheretal.Y.coal(df.tpairs, clumsm.info, X.pt, cluphy, cluphy.info, cluphy.map.nodectime, coal.t.Uscore.min=0.01, coal.within.inf.grace= 0.25, t.period=t.period, save.file=file )
		#	U [0,1]: prob that pot transmitter is still infected during infection window
		Y.U						<- project.athena.Fisheretal.Y.transmitterinfected(Y.infwindow, X.pt)
		#
		#	screen for likely missed intermediates/sources
		#
		df.tpairs				<- project.athena.Fisheretal.Y.rm.missedtransmitter(df.tpairs, Y.brl, Y.coal, Y.U, cut.date=0.5, cut.brl=1e-3)
		#
		#	take branch length weight ONLY to derive "score.Y"
		#		
		Y.score					<- merge( df.tpairs, subset(Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, brl, score.brl.TP)), by=c('FASTASampleCode','t.FASTASampleCode'))
		#	keep only one (Patient, t.Patient) pair if there are multiple sequences per individual. Take the one with largest brl score
		tmp						<- Y.score[,	{	
					z<- which.min(brl)
					list(FASTASampleCode=FASTASampleCode[z], t.FASTASampleCode=t.FASTASampleCode[z])
				},by=c('Patient', 't.Patient')]	
		Y.score					<- merge( Y.score, subset(tmp, select=c(FASTASampleCode=FASTASampleCode, t.FASTASampleCode)), by=c('FASTASampleCode', 't.FASTASampleCode') )
		tmp						<- Y.score[,	list(score.Y=score.brl.TP/sum(score.brl.TP), t.Patient=t.Patient), by='Patient']
		Y.score					<- merge( subset(Y.score, select=c(cluster, Patient, t.Patient, FASTASampleCode, t.FASTASampleCode )), tmp, by=c('Patient', 't.Patient') )		
		#	check that each (t, Patient, t.Patient) triplet is only once in the data.table
		tmp						<- Y.score[, list(n= length(FASTASampleCode)), by=c('Patient', 't.Patient')]
		if(tmp[,any(n!=1)])	stop('unexpected multiple entries per Patient, t.Patient')
		#	merge over time periods t
		Y.score					<- merge( Y.score, subset(Y.infwindow, select=c(Patient, t)), by='Patient', allow.cartesian=TRUE )
		YX						<- merge( subset(Y.score, select=c(FASTASampleCode, t.FASTASampleCode, score.Y, t)), X.pt, by= c('FASTASampleCode', 't.FASTASampleCode', 't'))
		#	add weights -- TODO same tpair can be counted multiple times
		YX						<- project.athena.Fisheretal.v2.Y.weight(YX)
		#	re-arrange a little
		YX		<- subset(YX, select=	c(	t, t.Patient, Patient, cluster, score.Y, 																	#triplets identifiers and Y score
						t.period, stage, U.score, contact, CDCC, lRNA, CD4, 														#main covariates							
						t.Age, t.AnyPos_a, RegionHospital, ART.I, ART.F, ART.A, ART.P, ART.pulse, fw.up, t2.care.t1, t2.vl.supp, 		#secondary covariates
						t.AnyPos_T1,  t.AnyT_T1, StartTime, StopTime, lRNAc, t.isAcute, Trm,												#other 
						FASTASampleCode, t.FASTASampleCode, w										#other							
				))
		setkey(YX, t, t.Patient, Patient)
		
		#
		if(!is.na(save.file))
		{
			cat(paste('\nsave YX to file', save.file))
			save(YX, file=save.file)
		}
	}
	YX	
}
######################################################################################
project.athena.Fisheretal.YX.part2<- function(YX.part1, df.all, predict.t2inf, t2inf.args, indir, insignat, indircov, infilecov, infiletree, outdir, outfile, cluphy=NULL, cluphy.info=NULL, cluphy.map.nodectime=NULL, df.tpairs.4.rawbrl=NULL, rm.zero.score=FALSE, t.period=0.25, save.file=NA, resume=1, method='3aa')
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
		file					<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'nbrlraw',method,sep='')	
		tmp						<- project.athena.Fisheretal.Y.rawbrl(YX.tpairs, indir, insignat, indircov, infilecov, infiletree, df.tpairs.tptn=df.tpairs.4.rawbrl, save.file=paste(file, '.R', sep=''), resume=resume, plot.file=paste(file, '.pdf', sep=''), method.restrictTPtoRI=method.restrictTPtoRI)
		Y.rawbrl				<- tmp$tpairs
		Y.rawbrl.linked			<- tmp$linked
		Y.rawbrl.unlinked		<- tmp$unlinked
		#	BRL [0,1]: branch length weight between pot transmitter and infected
		plot.file.score			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl',method,'_score','.pdf',sep='')
		plot.file.both			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl',method,'_both','.pdf',sep='')
		plot.file.one			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl',method,'_one','.pdf',sep='')
		if(substr(method,1,2)=='3c')		
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, brl.linked.max.brlr= 0.04, brl.linked.min.brl= 1e-4, brl.linked.max.dt= 10, brl.linked.min.dt= 1, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one)
		if(substr(method,1,2)=='3d')		
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, brl.linked.max.brlr= 0.01, brl.linked.min.brl= 1e-12, brl.linked.max.dt= 10, brl.linked.min.dt= 1, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one)
		if(all(substr(method,1,2)!=c('3c','3d')))	stop('brlweight: method not supported')		
		#	U [0,1]: prob that pot transmitter is still infected at time t. Needed to determine time of infection for transmitter (as quantile of the surival distribution)
		Y.U						<- project.athena.Fisheretal.Y.infectiontime(YX.tpairs, df.all, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, score.min=0.01, score.set.value=NA, method='for.transmitter')
		#	COAL [0,1]: prob that coalescence is within the transmitter
		Y.coal					<- NULL
		if(!is.null(cluphy) & !is.null(cluphy.info) & !is.null(cluphy.map.nodectime)  )
		{
			file				<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'coalraw',method,'.R',sep='')		
			Y.coal				<- project.athena.Fisheretal.Y.coal(YX.tpairs, df.all, Y.U, cluphy, cluphy.info, cluphy.map.nodectime, coal.t.Uscore.min=0.01, coal.within.inf.grace= 0.25, t.period=t.period, save.file=file, resume=resume )			
		}
		#	U [0,1]: prob that pot transmitter is still infected during infection window
		Y.U						<- project.athena.Fisheretal.Y.transmitterinfected(YX.part1)		
		#	screen for likely missed intermediates/sources
		plot.file				<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'missed_',method,'.pdf',sep='')		
		YX.tpairs				<- project.athena.Fisheretal.Y.rm.missedtransmitter(YX.tpairs, df.all, Y.brl, Y.U, Y.coal=Y.coal, cut.date=0.5, cut.brl=1e-3, any.pos.grace.yr= 3.5, rm.zero.score=rm.zero.score, plot.file=plot.file)
		#	take branch length weight ONLY to derive "score.Y"
		Y.score					<- merge( YX.tpairs, subset(Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, brl, score.brl.TPp, score.brl.TPd)), by=c('FASTASampleCode','t.FASTASampleCode'))
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
		YX						<- merge( YX.part1, subset(Y.score, select=c(FASTASampleCode, t.FASTASampleCode, score.Y, score.brl.TPd)), by=c('FASTASampleCode', 't.FASTASampleCode'))
		#	only triplets currently lost because of missing BEAST2 dated clusters
		YX.part1				<- NULL
		print(YX[, table(score.Y>0, class)])
		gc()
		#	add weights 		
		YX						<- project.athena.Fisheretal.YX.weight(YX)
		#	add balanced calendar periods 		
		YX						<- project.athena.Fisheretal.YX.calendarperiod.byweight(YX, df.all, t.period=t.period, c.nperiod= 4)
		#	add age of transmitter and recipient
		YX						<- project.athena.Fisheretal.YX.age(YX, df.all, t.period=t.period)
		#	add RegionHospital and Exposure group of transmitter and recipient
		YX						<- project.athena.Fisheretal.YX.Trm.Region(YX, df.all)		
		#	re-arrange a little
		YX						<- subset(YX, select=	c(	t, t.Patient, Patient, cluster, score.Y, 																	#triplets identifiers and Y score
										t.period, stage, U.score, contact, CDCC, lRNA, CD4, 														#main covariates							
										t.Age, Age, t.RegionHospital, RegionHospital, ART.I, ART.F, ART.A, ART.P, ART.pulse, ART.nDrug, ART.nNRT, ART.nNNRT, ART.nPI, fw.up.mx, fw.up.med, t2.care.t1, t2.vl.supp, 		#secondary covariates
										t.AnyPos_T1,  t.AnyT_T1, StartTime, StopTime, lRNAc, t.isAcute, t.Trm, Trm,												#other 
										FASTASampleCode, t.FASTASampleCode, w, w.i, w.t, class										#other							
								))
		setkey(YX, t, t.Patient, Patient)		
		#
		if(!is.na(save.file))
		{
			cat(paste('\nsave YX to file', save.file))
			save(YX, file=save.file)
		}
	}
	YX	
}
######################################################################################
project.athena.Fisheretal.v3.YX<- function(df.all, clumsm.info, df.tpairs, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, cluphy, cluphy.info, cluphy.map.nodectime, insignat, indircov, infilecov, infiletree, outdir, outfile, t.period=0.25, t.endctime=2013., save.file=NA, resume=1, method='3aa')
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
		#	get an idea how many time periods per year we could have:
		#	tmp<- project.athena.Fisheretal.median.followup(df.tpairs, df.viro, df.immu, X.tperiod)
		#	compute X: InCare
		X.incare				<- project.athena.Fisheretal.X.incare(df.tpairs, clumsm.info, df.viro, df.immu, df.treatment, t.period=t.period, t.endctime=t.endctime)
		X.incare				<- project.athena.Fisheretal.X.nocontact(X.incare, df.viro, df.immu, df.tpairs, clumsm.info, contact.grace=0.5, t.period=t.period, t.endctime= t.endctime)		
		X.incare				<- project.athena.Fisheretal.X.CDCC(X.incare, df.tpairs, clumsm.info, t.period=t.period, t.endctime=t.endctime)
		X.incare				<- project.athena.Fisheretal.X.followup(X.incare, clumsm.info, df.immu, t.period=t.period, t.endctime=t.endctime)
		#	compute X: before care
		X.b4care				<- project.athena.Fisheretal.X.b4care(df.tpairs, clumsm.info, predict.t2inf, t2inf.args, t.period=t.period)
		tmp						<- merge(X.incare, X.b4care, by=c('t.Patient','t'))
		cat(paste('\nnumber entries (Patient,t) that overlap in and before care [should be zero], n=',nrow(tmp)))
		X.pt					<- merge(X.incare, X.b4care, by=c('t.Patient','t'), all.x=1, all.y=1)
		set(X.pt, X.pt[,which(is.na(stage))], 'stage', 'U')					
		#	complete of transmitter: AnyPos_T1, AnyT_T1, AnyPos_a, isAcute
		X.pt[, AnyPos_T1:=NULL]
		X.pt[, AnyT_T1:=NULL]
		X.pt[, AnyPos_a:=NULL]
		X.pt[, isAcute:=NULL]
		tmp						<- unique(subset( clumsm.info, select=c(Patient, DateBorn, AnyPos_T1, AnyT_T1, isAcute) ))
		set(tmp, NULL, 'DateBorn', tmp[, AnyPos_T1-DateBorn])
		setnames( tmp, 'DateBorn', 'AnyPos_a')
		setnames( tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		X.pt					<- merge(X.pt, tmp, by='t.Patient', allow.cartesian=TRUE)
		
		#tmp					<- project.athena.Fisheretal.X.followup.compareCD4toVL(df.tpairs, df.immu, clumsm.info)
		#	compute X: potential transmitters with pulsed ART shortly after seroconversion / diagnosis
		X.ARTpulsed				<- project.athena.Fisheretal.X.ART.pulsed(df.tpairs, clumsm.info, t.pulse.ART=0.75, t.pulse.sc=NA, t.pulse.ART.I=1, t.pulse.ART.I.d=1)
		#	project.athena.Fisheretal.X.incare.check(X.incare)
		#	compute X: time to first viral suppression from diagnosis	
		X.t2.vlsupp				<- project.athena.Fisheretal.X.time.diag2suppressed(df.tpairs, clumsm.info, df.viro, lRNA.suppressed= log10(1e3), t2.vl.supp.p=c(0.1, 0.25))
		#	compute X: time to first VL and first CD4 measurement
		X.t2.care				<- project.athena.Fisheretal.X.time.diag2firstVLandCD4(df.tpairs, clumsm.info, df.viro, df.immu, t2.care.t1.q=c(0.25,0.5))
		#
		#	merge all static X that are independent of time periods
		#		
		tmp						<- merge( df.tpairs, X.t2.care, by='t.Patient', all.x=1 )
		tmp						<- merge( tmp, X.t2.vlsupp, by='t.Patient', all.x=1 )
		tmp						<- merge( tmp, X.ARTpulsed, by='t.Patient', all.x=1 )	
		X.pt					<- merge(tmp, X.pt, by='t.Patient', allow.cartesian=TRUE)
		#
		#	compute infection window of recipient
		#
		Y.infwindow				<- project.athena.Fisheretal.Y.infectiontime(df.tpairs, clumsm.info, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, score.min=0.1, score.set.value=1, method='for.infected')		
		#
		#	compute Y scores
		#
		#	BRL	[0,1]: raw branch length between pot transmitter and infected
		method.restrictTPtoRI	<- ifelse(substr(method,1,2)%in%c('3a','3b'),1,0)
		file					<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'nbrlraw',method,sep='')	
		tmp						<- project.athena.Fisheretal.Y.rawbrl(df.tpairs, indir, insignat, indircov, infilecov, infiletree, save.file=paste(file, '.R', sep=''), resume=resume, plot.file=paste(file, '.pdf', sep=''), method.restrictTPtoRI=method.restrictTPtoRI)
		Y.rawbrl				<- tmp$tpairs
		Y.rawbrl.linked			<- tmp$linked
		Y.rawbrl.unlinked		<- tmp$unlinked
		#	BRL [0,1]: branch length weight between pot transmitter and infected
		plot.file.score			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl',method,'_score','.pdf',sep='')
		plot.file.both			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl',method,'_both','.pdf',sep='')
		plot.file.one			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'wbrl',method,'_one','.pdf',sep='')
		if(substr(method,1,2)=='3c')		
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, brl.linked.max.brlr= 0.04, brl.linked.min.brl= 1e-4, brl.linked.max.dt= 10, brl.linked.min.dt= 1, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one)
		if(substr(method,1,2)=='3d')		
			Y.brl				<- project.athena.Fisheretal.Y.brlweight(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, brl.linked.max.brlr= 0.01, brl.linked.min.brl= 1e-12, brl.linked.max.dt= 10, brl.linked.min.dt= 1, plot.file.score=plot.file.score, method=method, plot.file.both=plot.file.both, plot.file.one=plot.file.one)
		if(all(substr(method,1,2)!=c('3c','3d')))	stop('brlweight: method not supported')
		#	COAL [0,1]: prob that coalescence is within the transmitter
		file					<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'coalraw',method,'.R',sep='')		
		Y.coal					<- project.athena.Fisheretal.Y.coal(df.tpairs, clumsm.info, X.pt, cluphy, cluphy.info, cluphy.map.nodectime, coal.t.Uscore.min=0.01, coal.within.inf.grace= 0.25, t.period=t.period, save.file=file, resume=resume )
		#	U [0,1]: prob that pot transmitter is still infected during infection window
		Y.U						<- project.athena.Fisheretal.Y.transmitterinfected(Y.infwindow, X.pt)
		#	screen for likely missed intermediates/sources
		df.tpairs				<- project.athena.Fisheretal.Y.rm.missedtransmitter(df.tpairs, Y.brl, Y.U, Y.coal=Y.coal, cut.date=0.5, cut.brl=1e-3)		
		#	take branch length weight ONLY to derive "score.Y"
		Y.score					<- merge( df.tpairs, subset(Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, brl, score.brl.TPp, score.brl.TPd)), by=c('FASTASampleCode','t.FASTASampleCode'))
		#	keep only one (Patient, t.Patient) pair if there are multiple sequences per individual. Take the one with largest brl score
		tmp						<- Y.score[,	{	
					z<- which.min(brl)
					list(FASTASampleCode=FASTASampleCode[z], t.FASTASampleCode=t.FASTASampleCode[z])
				},by=c('Patient', 't.Patient')]	
		Y.score					<- merge( Y.score, subset(tmp, select=c(FASTASampleCode=FASTASampleCode, t.FASTASampleCode)), by=c('FASTASampleCode', 't.FASTASampleCode') )
		setnames(Y.score, 'score.brl.TPp', 'score.Y')
		#
		tmp						<- Y.score[, list(n= length(FASTASampleCode)), by=c('Patient', 't.Patient')]
		if(tmp[,any(n!=1)])	stop('unexpected multiple entries per Patient, t.Patient')
		#	merge over time periods t
		Y.score					<- merge( Y.score, subset(Y.infwindow, select=c(Patient, t)), by='Patient', allow.cartesian=TRUE )
		YX						<- merge( subset(Y.score, select=c(FASTASampleCode, t.FASTASampleCode, score.Y, score.brl.TPd, t)), X.pt, by= c('FASTASampleCode', 't.FASTASampleCode', 't'))		
		#	add weights 		
		YX						<- project.athena.Fisheretal.YX.weight(YX)
		#	add balanced calendar periods 		
		YX						<- project.athena.Fisheretal.YX.calendarperiod.byweight(YX, clumsm.info, t.period=t.period, c.nperiod= 4)
		#	add age of transmitter and recipient
		YX						<- project.athena.Fisheretal.YX.age(YX, clumsm.info, t.period=t.period)
		#	add RegionHospital and Exposure group of transmitter and recipient
		YX						<- project.athena.Fisheretal.YX.Trm.Region(YX, clumsm.info)		
		#	re-arrange a little
		YX		<- subset(YX, select=	c(	t, t.Patient, Patient, cluster, score.Y, 																	#triplets identifiers and Y score
						t.period, stage, U.score, contact, CDCC, lRNA, CD4, 														#main covariates							
						t.Age, Age, t.RegionHospital, RegionHospital, ART.I, ART.F, ART.A, ART.P, ART.pulse, ART.nDrug, ART.nNRT, ART.nNNRT, ART.nPI, fw.up.mx, fw.up.med, t2.care.t1, t2.vl.supp, 		#secondary covariates
						t.AnyPos_T1,  t.AnyT_T1, StartTime, StopTime, lRNAc, t.isAcute, t.Trm, Trm,												#other 
						FASTASampleCode, t.FASTASampleCode, w, w.t										#other							
				))
		setkey(YX, t, t.Patient, Patient)		
		#
		if(!is.na(save.file))
		{
			cat(paste('\nsave YX to file', save.file))
			save(YX, file=save.file)
		}
	}
	YX	
}
######################################################################################
project.athena.Fisheretal.v1.Y.brlweight<- function(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, linked.brl.min= -4, plot.file=NA)
{
	#	check if Exp model would be reasonable
	require(MASS)
	unlinked.x			<- quantile( Y.rawbrl.unlinked[,brl], p=1e-3 )	 
	#	compute sequence sampling times	
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, unique(subset(df.all, select=c(FASTASampleCode, PosSeqT))), by='FASTASampleCode' )
	
	tmp					<- merge( data.table(FASTASampleCode=Y.rawbrl.linked[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, PosSeqT))), by='FASTASampleCode' )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, tmp, by='t.FASTASampleCode')
	set(Y.rawbrl.linked, NULL, 'PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,PosSeqT]))
	set(Y.rawbrl.linked, NULL, 't.PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,t.PosSeqT]))
	Y.rawbrl.linked[, dt:= Y.rawbrl.linked[,abs(t.PosSeqT-PosSeqT)]]
	
	#plot(Y.rawbrl.linked[,dt], Y.rawbrl.linked[,log(brl)], pch=18)		#no exponential clock within hosts 
	if(0)
	{
		rawbrl.exp			<- sapply(c(0,1,2,3,4), function(dt.min)
				{
					tmp					<- subset(Y.rawbrl.linked, dt>=dt.min)				
					rawbrl.l	<- subset( tmp, brl>=0 )[, sort(brl)]		#allow for one mutation
					#plot(log10(rawbrl.l))
					rawbrl.exp	<- fitdistr(rawbrl.l, "exponential")$estimate 				
				})
		#rate estimates very similar irrespective of whether the first 1 2 3 years excluded:
		#1/rawbrl.exp	:	0.009735224 0.011242407 0.009690996 0.010396757 0.012383081
	}	
	Y.rawbrl.linked	<- subset(Y.rawbrl.linked, log10(brl)>linked.brl.min & brl<0.2)
	rawbrl.l		<- Y.rawbrl.linked[, sort(brl)]		#suppose evolution not 'static'				
	#plot(log10(rawbrl.l))
	rawbrl.exp		<- fitdistr(rawbrl.l, "exponential")$estimate
	#rawbrl.h	<- hist(Y.rawbrl.linked[,brl], breaks=1e3, freq=0, xlim=c(0,0.05))
	#tmp			<- seq(min(Y.rawbrl.linked[,brl]), max(Y.rawbrl.linked[,brl]), length.out=1024) 
	#lines( tmp, dexp(tmp, rate = rawbrl.exp$estimate), col='red')
	Y.brl			<- Y.rawbrl
	Y.brl[, wbrl:= brl]
	tmp				<- pexp(unlinked.x, rate = rawbrl.exp)
	set(Y.brl, NULL, 'wbrl', pexp(Y.brl[,brl], rate = rawbrl.exp, lower.tail=FALSE) / tmp)		
	set(Y.brl, Y.brl[,which(brl>unlinked.x)], 'wbrl', 0 )
	Y.brl[, w2brl:= pgamma(Y.brl[,brl], shape=2, rate = rawbrl.exp, lower.tail=FALSE)]		
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
		tmp			<- seq(min(Y.rawbrl[,brl]), max(Y.rawbrl[,brl]), length.out=1024)
		tmp2		<- pgamma(tmp, shape=2, rate = rawbrl.exp, lower.tail=FALSE) 
		lines(tmp, tmp2/(sum(tmp2)*diff(tmp)[1]), col=cols[3], lwd=2)		
		legend('topright', bty='n', fill= cols, border=NA, legend=legend.txt)		
		dev.off()
	}	
	Y.brl		<- subset( Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, w2brl) )
	setnames(Y.brl, 'w2brl', 'score.brl')
	set(Y.brl, NULL, 'score.brl', Y.brl[, round(score.brl, d=3)])
	Y.brl		
}
######################################################################################
project.athena.Fisheretal.v2.Y.brlweight<- function(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, brl.linked.min.brl=-4, brl.linked.min.dt=1, plot.file=NA)
{
	#	check if Exp model would be reasonable
	require(MASS)
	setkey(Y.rawbrl.unlinked, brl)
	#	compute cdf for pairs given they are unlinked
	p.rawbrl.unlinked	<- Y.rawbrl.unlinked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]		 
	#	compute sequence sampling times	to see if time between seq sampling times could be useful
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, unique(subset(df.all, select=c(FASTASampleCode, PosSeqT))), by='FASTASampleCode' )	
	tmp					<- merge( data.table(FASTASampleCode=Y.rawbrl.linked[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, PosSeqT))), by='FASTASampleCode' )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, tmp, by='t.FASTASampleCode')
	set(Y.rawbrl.linked, NULL, 'PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,PosSeqT]))
	set(Y.rawbrl.linked, NULL, 't.PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,t.PosSeqT]))
	Y.rawbrl.linked[, dt:= Y.rawbrl.linked[,abs(t.PosSeqT-PosSeqT)]]	
	#plot(Y.rawbrl.linked[,dt], Y.rawbrl.linked[,log(brl)], pch=18)		#no exponential clock within hosts 
	if(0)
	{
		rawbrl.exp			<- sapply(c(0,1,2,3,4), function(dt.min)
				{
					tmp					<- subset(Y.rawbrl.linked, dt>=dt.min & log10(brl)>brl.linked.min.brl & brl<0.2)				
					rawbrl.l	<- subset( tmp, brl>=0 )[, sort(brl)]		#allow for one mutation
					#plot(log10(rawbrl.l))
					rawbrl.exp	<- fitdistr(rawbrl.l, "exponential")$estimate 				
				})
		#rate estimates very similar irrespective of whether the first 1 2 3 years excluded:
		#1/rawbrl.exp	:	0.009735224 0.011242407 0.009690996 0.010396757 0.012383081
	}	
	Y.rawbrl.linked	<- subset(Y.rawbrl.linked, log10(brl)>brl.linked.min.brl & brl<0.2 & dt>=brl.linked.min.dt)
	rawbrl.l		<- Y.rawbrl.linked[, sort(brl)]		#suppose evolution not 'static'				
	#plot(log10(rawbrl.l))
	rawbrl.exp		<- fitdistr(rawbrl.l, "exponential")$estimate
	cat(paste('\nestimated mean=',round(1/rawbrl.exp,d=4)))
	#rawbrl.h	<- hist(Y.rawbrl.linked[,brl], breaks=1e3, freq=0, xlim=c(0,0.05))
	#tmp			<- seq(min(Y.rawbrl.linked[,brl]), max(Y.rawbrl.linked[,brl]), length.out=1024) 
	#lines( tmp, dexp(tmp, rate = rawbrl.exp$estimate), col='red')
	Y.brl			<- Y.rawbrl
	Y.brl[, score.brl.TP:= dexp(brl, rate=rawbrl.exp/2) ]
	Y.brl[, score.brl.TN:= p.rawbrl.unlinked(brl)]	
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
		tmp			<- seq(min(Y.rawbrl[,brl]), max(Y.rawbrl[,brl]), length.out=1024)
		tmp2		<- dexp(tmp, rate = rawbrl.exp/2) 
		lines(tmp, tmp2, col=cols[3], lwd=2)		
		legend('topright', bty='n', fill= cols, border=NA, legend=legend.txt)		
		dev.off()
	}	
	Y.brl		<- subset( Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, score.brl.TP, score.brl.TN, brl) )
	cat(paste('\nReturn brl score for #t.Patient seqs=',Y.brl[, length(unique(t.FASTASampleCode))]))
	Y.brl		
}
######################################################################################
project.athena.Fisheretal.Y.brlweight<- function(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, brl.linked.max.brlr= 0.04, brl.linked.min.brl= 1e-4, brl.linked.max.dt= 10, brl.linked.min.dt= 1.5, plot.file.score=NA, method='3c', plot.file.both=NA, plot.file.one=NA)
{
	#brl.linked.max.brlr= 0.01; brl.linked.min.brl= 1e-12; brl.linked.max.dt= 10; brl.linked.min.dt= 1; plot.file.score=NA; method='3aa'; plot.file.both=NA; plot.file.one=NA
	#	check if Exp model would be reasonable
	require(MASS)	
	require(reshape2)
	require(ggplot2)
	require(gamlss)
	
	setkey(Y.rawbrl.unlinked, brl)
	#	compute cdf for pairs given they are unlinked
	p.rawbrl.unlinked	<- Y.rawbrl.unlinked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]		 
	#	compute sequence sampling times	to see if time between seq sampling times could be useful
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )	
	tmp					<- merge( data.table(FASTASampleCode=Y.rawbrl.linked[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, tmp, by='t.FASTASampleCode')
	if(class(df.all$PosSeqT)=='Date')
		set(Y.rawbrl.linked, NULL, 'PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,PosSeqT]))
	if(class(df.all$t.PosSeqT)=='Date')
		set(Y.rawbrl.linked, NULL, 't.PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,t.PosSeqT]))
	if(class(df.all$AnyT_T1)=='Date')
		set(Y.rawbrl.linked, NULL, 'AnyT_T1', hivc.db.Date2numeric(Y.rawbrl.linked[,AnyT_T1]))
	if(class(df.all$t.AnyT_T1)=='Date')
		set(Y.rawbrl.linked, NULL, 't.AnyT_T1', hivc.db.Date2numeric(Y.rawbrl.linked[,t.AnyT_T1]))
	Y.rawbrl.linked[, b4T:= 'both']
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT<=AnyT_T1)], 'b4T', 'only.RI')
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT<=t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'only.T')
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'none')
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, data.table(b4T= Y.rawbrl.linked[, unique(b4T)], col=c('green','red','orange','yellow')), by='b4T' )	
	set(Y.rawbrl.linked, NULL, 'b4T', Y.rawbrl.linked[, factor(b4T)])
	Y.rawbrl.linked[, dt:= abs(PosSeqT-t.PosSeqT)]		
	Y.rawbrl.linked[, brlr:= brl/dt]
	#
	#	EXCLUDE options
	#
	Y.rawbrl.linked		<- subset( Y.rawbrl.linked, dt!=0 )
	
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
	# wilcox test for same mean
	#
	test	<- list( 	both.only.T= 	wilcox.test( subset(Y.rawbrl.linked, b4T=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.T')[, log10(brl)] ),
						both.only.RI= 	wilcox.test( subset(Y.rawbrl.linked, b4T=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.RI')[, log10(brl)] ),
						both.none= 		wilcox.test( subset(Y.rawbrl.linked, b4T=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='none')[, log10(brl)] ),
						only.RI.T= 		wilcox.test( subset(Y.rawbrl.linked, b4T=='only.RI')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.T')[, log10(brl)] ),
						none.only.T= 	wilcox.test( subset(Y.rawbrl.linked, b4T=='none')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.T')[, log10(brl)] ),
						none.only.RI= 	wilcox.test( subset(Y.rawbrl.linked, b4T=='none')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.RI')[, log10(brl)] )	)
	cat(paste('\nwilcox test for same mean'))
	print( sapply(test, '[[', 'p.value') )
	# KS test for same distribution
	test	<- list( 	both.only.T= 	ks.test( subset(Y.rawbrl.linked, b4T=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.T')[, log10(brl)]),
						both.only.RI= 	ks.test( subset(Y.rawbrl.linked, b4T=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.RI')[, log10(brl)] ),
						both.none= 		ks.test( subset(Y.rawbrl.linked, b4T=='both')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='none')[, log10(brl)] ),
						only.RI.T= 		ks.test( subset(Y.rawbrl.linked, b4T=='only.RI')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.T')[, log10(brl)] ),
						none.only.T= 	ks.test( subset(Y.rawbrl.linked, b4T=='none')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.T')[, log10(brl)] ),
						none.only.RI= 	ks.test( subset(Y.rawbrl.linked, b4T=='none')[, log10(brl)], subset(Y.rawbrl.linked, b4T=='only.RI')[, log10(brl)] )	)
	cat(paste('\nKS test for same distribution'))
	print( sapply(test, '[[', 'p.value') )	
	#
	#	fit distributions to data 
	#
	require(fitdistrplus)
	require(ADGofTest)
	require(gamlss.mx)
	#
	#	Gamma shape mixture 
	#
	if(0)
	{
		require(GSM)
		tmp			<- estim.gsm_theta(subset(Y.rawbrl.linked, b4T=='both')[,brl], 2, 300, burnin + mcmcsim, 1, 1, 1/2)
		tmp3		<- estim.gsm_theta(subset(Y.rawbrl.linked, b4T=='both')[,brl], 3, 300, burnin + mcmcsim, 1, 1, 1/3)
		tmp4		<- estim.gsm_theta(subset(Y.rawbrl.linked, b4T=='both')[,brl], 4, 300, burnin + mcmcsim, 1, 1, 1/4)
		plot(tmp, ndens = 0, nbin = 50, histogram = TRUE, start = (burnin + 1))
		ad.test(subset(Y.rawbrl.linked, b4T=='both')[,brl], pgsm, weight=apply(tmp@weight[seq.int(burnin+1,burnin+mcmcsim),], 2, median), rateparam=median(tmp@theta[seq.int(burnin+1,burnin+mcmcsim)]))
		#does not work: pvalue is 1e-6 almost Exponential		
	}
	#
	#	Gamma mixture 
	#
	if(0)
	{
		brl				<- subset(Y.rawbrl.linked, b4T=='both')[,brl]
		ml.gam			<- lapply(1:5, function(gam.n) gamlssMX(brl~1, family=GA, K=gam.n) )
		ml.gam.n		<- which.min(sapply(ml.gam, AIC))				
		cat(paste('\nbest AIC Gamma mixture has components n=', ml.gam.n))
		ml.gam.p		<- ad.test(brl, pMX, 	mu=as.list(exp(sapply(ml.gam[[ml.gam.n]]$models,'[[','mu.coefficients'))), 
				sigma=as.list(exp(sapply(ml.gam[[ml.gam.n]]$models,'[[','sigma.coefficients'))), 
				pi=as.list(ml.gam[[ml.gam.n]]$prob), family=ml.gam[[ml.gam.n]]$family 	)$p.value
		cat(paste('\nbest AIC Gamma mixture has pval=', ml.gam.p))		
		
		brl				<- subset(Y.rawbrl.linked, b4T!='both' & b4T!='none')[,brl]
		ml.gam.one		<- lapply(1:5, function(gam.n) gamlssMX(brl~1, family=GA, K=gam.n) )
		ml.gam.one.n	<- which.min(sapply(ml.gam.one, AIC))				
		cat(paste('\nbest AIC Gamma mixture has components n=', ml.gam.one.n))
		ml.gam.one.p	<- ad.test(brl, pMX, 	mu=as.list(exp(sapply(ml.gam.one[[ml.gam.one.n]]$models,'[[','mu.coefficients'))), 
				sigma=as.list(exp(sapply(ml.gam.one[[ml.gam.one.n]]$models,'[[','sigma.coefficients'))), 
				pi=as.list(ml.gam.one[[ml.gam.one.n]]$prob), family=ml.gam.one[[ml.gam.one.n]]$family 	)$p.value
		cat(paste('\nbest AIC Gamma mixture has pval=', ml.gam.one.p))									
	}
	#
	#	zero adjusted Gamma 
	#
	if(1)
	{
		Y.rawbrl.linked[, brlz:=brl]
		tmp				<- Y.rawbrl.linked[, which(brlz<1e-4)]
		set(Y.rawbrl.linked, tmp, 'brlz', 0)
		brl				<- subset(Y.rawbrl.linked, b4T=='both'  )[,brlz]
		ml.zaga			<- gamlss(brl~1, family=ZAGA)
		ml.zaga.p		<- ad.test(brl, pZAGA, mu=exp(ml.zaga[['mu.coefficients']]), sigma=exp(ml.zaga[['sigma.coefficients']]), nu=1/(1+exp(-ml.zaga[['nu.coefficients']])) )$p.value
		cat(paste('\nbest ZAGA pvalue for both b4T=', ml.zaga.p))
		brl				<- subset(Y.rawbrl.linked, b4T!='both' & b4T!='none'  & brlr<=0.01)[,brlz]
		ml.zaga.one		<- gamlss(brl~1, family=ZAGA)
		ml.zaga.one.p	<- ad.test(brl, pZAGA, mu=exp(ml.zaga.one[['mu.coefficients']]), sigma=exp(ml.zaga.one[['sigma.coefficients']]), nu=1/(1+exp(-ml.zaga.one[['nu.coefficients']])) )$p.value
		cat(paste('\nbest ZAGA pvalue for one b4T=', ml.zaga.one.p))			
	}
	#
	#	others
	#
	mle.exp		<- fitdist(subset(Y.rawbrl.linked, b4T=='both')[,brl], 'exp')
	mle.ga		<- fitdist(subset(Y.rawbrl.linked, b4T=='both')[,brl], 'gamma')
	mle.ln		<- fitdist(subset(Y.rawbrl.linked, b4T=='both')[,brl], 'lnorm')
	mle.w		<- fitdist(subset(Y.rawbrl.linked, b4T=='both')[,brl], 'weibull')
	mle.exp.one	<- fitdist(subset(Y.rawbrl.linked, b4T!='both' & b4T!='none')[,brl], 'exp')
	mle.ga.one	<- fitdist(subset(Y.rawbrl.linked, b4T!='both' & b4T!='none')[,brl], 'gamma')
	mle.ln.one	<- fitdist(subset(Y.rawbrl.linked, b4T!='both' & b4T!='none')[,brl], 'lnorm')
	mle.w.one	<- fitdist(subset(Y.rawbrl.linked, b4T!='both' & b4T!='none')[,brl], 'weibull')
	#	plot
	if(!is.na(plot.file.both))
	{
		pdf(file=plot.file.both, w=12,h=5)
		par(mfrow=c(2, 4))
		qqcomp(mle.exp, addlegend=FALSE, main='exp')
		denscomp(mle.exp, addlegend=FALSE, xlab='exp')		
		qqcomp(mle.ga, addlegend=FALSE, main='gamma')
		denscomp(mle.ga, addlegend=FALSE, xlab='gamma')
		qqcomp(mle.ln, addlegend=FALSE, main='ln')
		denscomp(mle.ln, addlegend=FALSE, xlab='ln')
		qqcomp(mle.w, addlegend=FALSE, main='weibull')
		denscomp(mle.w, addlegend=FALSE, xlab='weibull')
		dev.off()
		###		zero adjusted for ZAGA
		brl				<- subset(Y.rawbrl.linked, b4T=='both')[,brlz]
		brl.cdf			<- data.table( brl=sort(brl), cdf=seq_along(brl)/length(brl) )[, list(cdf=tail(cdf, 1)), by='brl']
		brl.cdf[, cdf.f:= pZAGA( brl.cdf[,brl], mu=exp(ml.zaga[['mu.coefficients']]), sigma=exp(ml.zaga[['sigma.coefficients']]), nu=1/(1+exp(-ml.zaga[['nu.coefficients']])) ) ]				
		ggplot( data=melt(brl.cdf, id='brl'), aes(x=brl, y=value, colour=variable)) + geom_line()
		ggsave(file=paste(substr(plot.file.both,1,nchar(plot.file.both)-4),'_zagacdf','.pdf',sep=''), w=5,h=5)
		###		
		qv 		<- quantile(brl, c(.25, .75))
		qt 		<- qZAGA(c(.25, .75), mu=exp(ml.zaga[['mu.coefficients']]), sigma=exp(ml.zaga[['sigma.coefficients']]), nu=1/(1+exp(-ml.zaga[['nu.coefficients']])) )
		slope 	<- diff(qv)/diff(qt)
		int 	<- qv[1] - slope * qt[1]	
		qplot(sample=brl, geom = "point", stat = "qq", distribution = qZAGA, dparams = list(mu=exp(ml.zaga[['mu.coefficients']]), sigma=exp(ml.zaga[['sigma.coefficients']]), nu=1/(1+exp(-ml.zaga[['nu.coefficients']])))) + geom_abline(slope = slope, intercept = int)
		ggsave(file=paste(substr(plot.file.both,1,nchar(plot.file.both)-4),'_zagaqq','.pdf',sep=''), w=5,h=5)
	}
	if(0)
	{
		brl				<- subset(Y.rawbrl.linked, b4T=='both')[,brl]
		pdf(file=paste(substr(plot.file.both,1,nchar(plot.file.both)-4),'_gam','.pdf',sep=''), w=12,h=5)	
		hist(brl, breaks=100, freq=FALSE)
		dummy		<- sapply(seq_along(ml.gam), function(ml.gam.i)
				{
					print(ml.gam.i)
					dbrl		<- dMX(y=seq(min(brl),max(brl),len=1024), mu=as.list(exp(sapply(ml.gam[[ml.gam.i]]$models,'[[','mu.coefficients'))), sigma=as.list(exp(sapply(ml.gam[[ml.gam.i]]$models,'[[','sigma.coefficients'))), pi=as.list(ml.gam[[ml.gam.i]]$prob), family=ml.gam[[ml.gam.i]]$family )
					lines(seq(min(brl),max(brl),len=1024), dbrl, col=rainbow(length(ml.gam))[ml.gam.i], lty=ml.gam.i)
				})
		legend('topright', fill=rainbow(length(ml.gam)), bty='n', border=NA, legend=paste('GA mixture cmpnts',seq_along(ml.gam)))
		dev.off()		
	}
	if(!is.na(plot.file.one))
	{
		pdf(file=plot.file.one, w=12,h=5)
		par(mfrow=c(2, 4))
		qqcomp(mle.exp.one, addlegend=FALSE, main='exp')
		denscomp(mle.exp.one, addlegend=FALSE, xlab='exp')		
		qqcomp(mle.ga.one, addlegend=FALSE, main='gamma')
		denscomp(mle.ga.one, addlegend=FALSE, xlab='gamma')
		qqcomp(mle.ln.one, addlegend=FALSE, main='ln')
		denscomp(mle.ln.one, addlegend=FALSE, xlab='ln')
		qqcomp(mle.w.one, addlegend=FALSE, main='weibull')
		denscomp(mle.w.one, addlegend=FALSE, xlab='weibull')
		dev.off()
		###		zero adjusted for ZAGA
		brl				<- subset(Y.rawbrl.linked, b4T!='both' & b4T!='none' )[,brlz]
		brl.cdf			<- data.table( brl=sort(brl), cdf=seq_along(brl)/length(brl) )[, list(cdf=tail(cdf, 1)), by='brl']
		brl.cdf[, cdf.f:= pZAGA( brl.cdf[,brl], mu=exp(ml.zaga.one[['mu.coefficients']]), sigma=exp(ml.zaga.one[['sigma.coefficients']]), nu=1/(1+exp(-ml.zaga.one[['nu.coefficients']])) ) ]						
		ggplot( data=melt(brl.cdf, id='brl'), aes(x=brl, y=value, colour=variable)) + geom_line()
		ggsave(file=paste(substr(plot.file.one,1,nchar(plot.file.one)-4),'_zagacdf','.pdf',sep=''), w=5,h=5)		
		###		
		qv 		<- quantile(brl, c(.25, .75))
		qt 		<- qZAGA(c(.25, .75), mu=exp(ml.zaga.one[['mu.coefficients']]), sigma=exp(ml.zaga.one[['sigma.coefficients']]), nu=1/(1+exp(-ml.zaga.one[['nu.coefficients']])) )
		slope 	<- diff(qv)/diff(qt)
		int 	<- qv[1] - slope * qt[1]	
		qplot(sample=brl, geom = "point", stat = "qq", distribution = qZAGA, dparams = list(mu=exp(ml.zaga.one[['mu.coefficients']]), sigma=exp(ml.zaga.one[['sigma.coefficients']]), nu=1/(1+exp(-ml.zaga.one[['nu.coefficients']])))) + geom_abline(slope = slope, intercept = int)
		ggsave(file=paste(substr(plot.file.one,1,nchar(plot.file.one)-4),'_zagaqq','.pdf',sep=''), w=5,h=5)
	}
	if(0)
	{
		pdf(file=paste(substr(plot.file.one,1,nchar(plot.file.one)-4),'_gam','.pdf',sep=''), w=12,h=5)
		brl				<- subset(Y.rawbrl.linked, b4T!='both' & b4T!='none')[,brl]
		hist(brl, breaks=100, freq=FALSE)
		dummy		<- sapply(seq_along(ml.gam.one), function(ml.gam.i)
				{
					dbrl		<- dMX(y=seq(min(brl),max(brl),len=1024), mu=as.list(exp(sapply(ml.gam.one[[ml.gam.i]]$models,'[[','mu.coefficients'))), sigma=as.list(exp(sapply(ml.gam.one[[ml.gam.i]]$models,'[[','sigma.coefficients'))), pi=as.list(ml.gam.one[[ml.gam.i]]$prob), family=ml.gam.one[[ml.gam.i]]$family )
					lines(seq(min(brl),max(brl),len=1024), dbrl, col=rainbow(length(ml.gam.one))[ml.gam.i])
				})
		legend('topright', fill=rainbow(length(ml.gam.one)), bty='n', border=NA, legend=paste('GA mixture cmpnts',seq_along(ml.gam.one)))
		dev.off()		
	}
	test		<- list(	exp= ad.test(subset(Y.rawbrl.linked, b4T=='both')[,brl], pexp, rate=mle.exp$estimate['rate']),								
							gamma= ad.test(subset(Y.rawbrl.linked, b4T=='both')[,brl], pgamma, shape=mle.ga$estimate['shape'], rate=mle.ga$estimate['rate'])		)
	cat(paste('\nAD test for MLE fit in distribution for both before ART'))
	print( sapply(test, '[[', 'p.value') )	
	test		<- list(	exp= ad.test(subset(Y.rawbrl.linked, b4T!='both' & b4T!='none')[,brl], pexp, rate=mle.exp.one$estimate['rate']),								
							gamma= ad.test(subset(Y.rawbrl.linked, b4T!='both' & b4T!='none')[,brl], pgamma, shape=mle.ga.one$estimate['shape'], rate=mle.ga.one$estimate['rate'])		)
	cat(paste('\nAD test for MLE fit in distribution for one before ART'))
	print( sapply(test, '[[', 'p.value') )	
	#
	#	fit Gamma for b4T==both and b4T one in ART
	#
	Y.brl			<- copy(Y.rawbrl)	
	Y.brl[, score.brl.TN:= p.rawbrl.unlinked(brl)]
	
	if(substr(method,1,2)=='3c')	#ignore zeros use Gamma
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
		
		Y.brl[, score.brl.TPp:= pgamma(Y.brl[,brl], shape=2*mle.ga$estimate['shape'], rate=mle.ga$estimate['rate'], lower.tail=FALSE)]
		tmp<- Y.brl[, which(b4T=='one')]
		set(Y.brl, tmp, 'score.brl.TPp', pgamma(Y.brl[tmp,brl], shape=2*mle.ga.one$estimate['shape'], rate=mle.ga.one$estimate['rate'], lower.tail=FALSE) )
		Y.brl[, score.brl.TPd:= score.brl.TPp]		
		Y.brl[, dt:= abs(PosSeqT-t.PosSeqT)]
	}	
	if(substr(method,1,2)=='3d')	#dont ignore zeros - use zero adjusted Gamma
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
		tmp				<- 2*ml.zaga.pa['nu']*(1-ml.zaga.pa['nu'])*pGA( Y.brl[,brlz],  mu=ml.zaga.pa['mu'], sigma=ml.zaga.pa['sigma'] ) + (1-ml.zaga.pa['nu'])*(1-ml.zaga.pa['nu'])*pgamma( Y.brl[,brlz],  shape=2*ml.zaga.pa['shape'], scale=ml.zaga.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.pa['nu']*ml.zaga.pa['nu'])				
		Y.brl[, score.brl.TPp:= 1-tmp]
		tmp				<- Y.brl[, which(b4T=='one')]
		tmp				<- 2*ml.zaga.one.pa['nu']*(1-ml.zaga.one.pa['nu'])*pGA( Y.brl[tmp,brlz],  mu=ml.zaga.one.pa['mu'], sigma=ml.zaga.one.pa['sigma'] ) + (1-ml.zaga.one.pa['nu'])*(1-ml.zaga.one.pa['nu'])*pgamma( Y.brl[tmp,brlz],  shape=2*ml.zaga.one.pa['shape'], scale=ml.zaga.one.pa['scale'] ) 
		tmp				<- tmp / (1-ml.zaga.one.pa['nu']*ml.zaga.one.pa['nu'])				
		set(Y.brl, Y.brl[, which(b4T=='one')], 'score.brl.TPp', 1-tmp )
		
		Y.brl[, score.brl.TPd:= score.brl.TPp]		
		Y.brl[, dt:= abs(PosSeqT-t.PosSeqT)]
	}	
	if(!is.na(plot.file.score))
	{
		pdf(plot.file.score, w=5, h=5)
		cat(paste('\nplot to file',plot.file.score))
		require(RColorBrewer)
		par(mar=c(4,4,0.5,0.5))
		cols		<- brewer.pal(4, 'Set1')[c(1,3,2)]
		legend.txt	<- c('same host', 'potential transmission pairs', 'branch length weight')	
		xlim		<- c(0, max(Y.brl[, max(brl)],  Y.rawbrl.linked[, max(brl)])*1.1 )
		tmp			<- seq(from=xlim[1], to=xlim[2], by=diff(xlim)/200)
		hist( Y.rawbrl.linked[, brl], breaks=tmp , col=my.fade.col(cols[1],0.5), add=0, freq=0, xlab='branch length', main='', ylim=c(0,50) )
		tmp2		<- hist( Y.brl[, brl], breaks=tmp, col=my.fade.col(cols[2],0.5), freq=0, add=1 )
		setkey(Y.brl, brl)
		lines(subset(Y.brl, b4T=='both')[, brl], subset(Y.brl, b4T=='both')[, score.brl.TPp]*max(tmp2$density), col=cols[3], lwd=2)
		lines(subset(Y.brl, b4T=='one')[, brl], subset(Y.brl, b4T=='one')[, score.brl.TPp]*max(tmp2$density), col=cols[3], lwd=2)
		legend('topright', bty='n', fill= cols, border=NA, legend=legend.txt)		
		dev.off()
		Y.brl[, dummy:= score.brl.TPp*Y.brl[, length(which(brl<=0.002))] ]
		ggplot(data=Y.brl, aes(x = brl, fill=b4T))+stat_bin(binwidth=0.002)+geom_line(data=Y.brl, aes(x=brl, y=dummy, colour=b4T))+ scale_color_manual(values=c('black', 'blue'))
		ggsave(file=paste(substr(plot.file.score,1,nchar(plot.file.score)-4),'_b4T','.pdf',sep=''), w=5, h=5)
		Y.brl[, dummy:=NULL]
	}		
	Y.brl		<- subset( Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, score.brl.TPd, score.brl.TPp, score.brl.TN, brl) )
	cat(paste('\nReturn brl score for #t.Patient seqs=',Y.brl[, length(unique(t.FASTASampleCode))]))
	Y.brl		
}
######################################################################################
project.athena.Fisheretal.v3a.Y.brlweight<- function(Y.rawbrl, Y.rawbrl.linked, Y.rawbrl.unlinked, df.all, brl.linked.min.brl=-4, brl.linked.min.dt=1, plot.file=NA, method='3aa')
{
	#brl.linked.min.brl=-4; brl.linked.min.dt=1
	#	check if Exp model would be reasonable
	require(MASS)
	setkey(Y.rawbrl.unlinked, brl)
	#	compute cdf for pairs given they are unlinked
	p.rawbrl.unlinked	<- Y.rawbrl.unlinked[, approxfun(brl , seq_along(brl)/length(brl), yleft=0., yright=1., rule=2)]		 
	#	compute sequence sampling times	to see if time between seq sampling times could be useful
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )	
	tmp					<- merge( data.table(FASTASampleCode=Y.rawbrl.linked[, unique(t.FASTASampleCode)]), unique(subset(df.all, select=c(FASTASampleCode, PosSeqT, AnyT_T1))), by='FASTASampleCode' )
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, tmp, by='t.FASTASampleCode')
	set(Y.rawbrl.linked, NULL, 'PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,PosSeqT]))
	set(Y.rawbrl.linked, NULL, 't.PosSeqT', hivc.db.Date2numeric(Y.rawbrl.linked[,t.PosSeqT]))
	set(Y.rawbrl.linked, NULL, 'AnyT_T1', hivc.db.Date2numeric(Y.rawbrl.linked[,AnyT_T1]))
	set(Y.rawbrl.linked, NULL, 't.AnyT_T1', hivc.db.Date2numeric(Y.rawbrl.linked[,t.AnyT_T1]))
	Y.rawbrl.linked[, b4T:= 'both']
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT<=AnyT_T1)], 'b4T', 'only.RI')
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT<=t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'only.T')
	set(Y.rawbrl.linked, Y.rawbrl.linked[, which(t.PosSeqT>t.AnyT_T1 & PosSeqT>AnyT_T1)], 'b4T', 'none')
	Y.rawbrl.linked		<- merge( Y.rawbrl.linked, data.table(b4T= Y.rawbrl.linked[, unique(b4T)], col=c('green','red','orange','yellow')), by='b4T' )	
	set(Y.rawbrl.linked, NULL, 'b4T', Y.rawbrl.linked[, factor(b4T)])
	Y.rawbrl.linked[, dt:= abs(PosSeqT-t.PosSeqT)]	
	#	subset(Y.rawbrl.linked, dt==0 & brl>1e-4)
	#	M36101 M10149 M15536 M15701 M29207 M28446 M10038 M31096
	Y.rawbrl.linked		<- subset( Y.rawbrl.linked, dt!=0 )
	Y.rawbrl.linked[, brlr:= brl/dt]
	if(0)
	{
		plot(Y.rawbrl.linked[,dt], Y.rawbrl.linked[,log(brl)], pch=18)		#no exponential clock within hosts
		tmp				<- loess(log(brl)~dt, Y.rawbrl.linked)
		tmp				<- predict(tmp, data.frame(dt = seq(0, 20, 0.5)), se=TRUE)
		lines(seq(0, 20, 0.5), tmp$fit, col='red')
		lines(seq(0, 20, 0.5), tmp$fit+2*tmp$se.fit, col='red')
		lines(seq(0, 20, 0.5), tmp$fit-2*tmp$se.fit, col='red')		
	}
	if(0)
	{
		brl.max	<- 0.2
		brl.max	<- 0.4
		setkey(Y.rawbrl.linked, brl, b4T)
		par(mfrow=c(1,5))
		#	plot distribution of BRL, excluding close SeqT
		dummy	<- sapply(c(0,1,2,3,4), function(dt.min)
				{					
					tmp					<- subset(Y.rawbrl.linked, dt>=dt.min &  brl<brl.max)				
					plot( tmp[, log10(brl)], pch=18, main=paste('exclude dt<',dt.min))					 				
				})
		#	plot distribution of BRL, excluding close SeqT & no evolution
		dummy	<- sapply(c(0,1,2,3,4), function(dt.min)
				{					
					tmp					<- subset(Y.rawbrl.linked, dt>=dt.min & log10(brl)>brl.linked.min.brl &  brl<brl.max)				
					plot( tmp[, log10(brl)], pch=18, main=paste('exclude dt<',dt.min))					 				
				})
		#	plot raw data for each b4T group, excluding close SeqT & no evolution
		tmp					<- subset(Y.rawbrl.linked, dt>=1 & log10(brl)>brl.linked.min.brl & brl<brl.max)
		setkey(tmp, b4T, brl)		
		par(mfrow=c(1,1))
		plot(1,1,type='n', ylim=c(-3, -0.5), xlim=c(1, tmp[, max(table(b4T))]), ylab='log10 brl', xlab='index')
		tmp[, points( seq_along(brl), log10(brl), col=col, pch=18 ), by='b4T']		
		legend('bottomright',legend=unique(subset(tmp, select=c(b4T, col)))[, b4T], fill=unique(subset(tmp, select=c(b4T, col)))[, col], bty='n', border=NA) 
		#	exact two sample permutation test for equal location
		library(coin)
		test	<- list( 	both.only.T= 	{	tmp2	<- subset( tmp, b4T%in%c('both','only.T') ); set(tmp2, NULL, 'b4T', tmp2[,factor(as.character(b4T))]); oneway_test( brl~b4T, data=tmp2, distribution = "exact", conf.int = TRUE ) 	},
				both.only.RI= 	{	tmp2	<- subset( tmp, b4T%in%c('both','only.RI') ); set(tmp2, NULL, 'b4T', tmp2[,factor(as.character(b4T))]); oneway_test( brl~b4T, data=tmp2, distribution = "exact", conf.int = TRUE ) 	},
				both.none= 		{	tmp2	<- subset( tmp, b4T%in%c('both','none') ); set(tmp2, NULL, 'b4T', tmp2[,factor(as.character(b4T))]); oneway_test( brl~b4T, data=tmp2, distribution = "exact", conf.int = TRUE ) 	},
				none.only.T= 	{	tmp2	<- subset( tmp, b4T%in%c('only.T','none') ); set(tmp2, NULL, 'b4T', tmp2[,factor(as.character(b4T))]); oneway_test( brl~b4T, data=tmp2, distribution = "exact", conf.int = TRUE ) 	},
				none.only.RI= 	{	tmp2	<- subset( tmp, b4T%in%c('only.RI','none') ); set(tmp2, NULL, 'b4T', tmp2[,factor(as.character(b4T))]); oneway_test( brl~b4T, data=tmp2, distribution = "exact", conf.int = TRUE ) 	}
		)
		sapply(test, pvalue)
		#within host among RI
		#both.only.T both.only.RI    both.none  none.only.T none.only.RI 
		#0.75287863   0.41671419   0.02926871   0.08318482   0.12930924
		#within host among all
		#
		#	KS test
		#
		test	<- list( 	both.only.T= 	ks.test( subset(tmp, b4T=='both')[, brl], subset(tmp, b4T=='only.T')[, brl] ),
				both.only.RI= 	ks.test( subset(tmp, b4T=='both')[, brl], subset(tmp, b4T=='only.RI')[, brl] ),
				both.none= 		ks.test( subset(tmp, b4T=='both')[, brl], subset(tmp, b4T=='none')[, brl] ),
				none.only.T= 	ks.test( subset(tmp, b4T=='none')[, brl], subset(tmp, b4T=='only.T')[, brl] ),
				none.only.RI= 	ks.test( subset(tmp, b4T=='none')[, brl], subset(tmp, b4T=='only.RI')[, brl] )	)
		sapply(test, '[[', 'p.value')
		#within host among RI
		#both.only.T both.only.RI    both.none  none.only.T none.only.RI 
		#0.7271106    0.6860559    0.1520723    0.3264294    0.5416300
		#within host among all
		# both.only.T 	both.only.RI    both.none  		none.only.T 	none.only.RI 
		#6.320585e-06 	1.341786e-04 	3.330669e-15 	5.086774e-04 	2.301351e-05 
		#
		#	plot points by b4T
		#
		tmp				<- subset(Y.rawbrl.linked, dt>=1 & log10(brl)>brl.linked.min.brl & brl<brl.max)
		plot(1,1,type='n', ylim=tmp[,range(log10(brl))], xlim=tmp[,range(dt)], ylab='log10 brl', xlab='delta SeqT')
		tmp[, points( dt, log10(brl), col=col, pch=18 ), by='b4T']
		legend('bottomright',legend=unique(subset(tmp, select=c(b4T, col)))[, b4T], fill=unique(subset(tmp, select=c(b4T, col)))[, col], bty='n', border=NA)
		#	expect rate 2e-3/site/year, so at most 4e-3/site/year pairwise divergence from between host evolution
		#	exclude 20-fold: n=136	out of 4090
		#	exclude 10-fold: n=323	out of 4090
		Y.rawbrl.linked	<- subset(Y.rawbrl.linked,  brlr<=0.04)		
		#	exclude dt>10, n=100
		Y.rawbrl.linked	<- subset(Y.rawbrl.linked,  dt<=10)		
		hist(Y.rawbrl.linked[, brl],breaks=100)
		#	exclude short dt to avoid bias 
		tmp				<- subset(Y.rawbrl.linked, dt>=1)
		plot(1,1,type='n', ylim=tmp[,range(log10(brl))], xlim=tmp[,range(dt)], ylab='log10 brl', xlab='delta SeqT')
		tmp[, points( dt, log10(brl), col=col, pch=18 ), by='b4T']
		legend('bottomright',legend=unique(subset(tmp, select=c(b4T, col)))[, b4T], fill=unique(subset(tmp, select=c(b4T, col)))[, col], bty='n', border=NA)
		setkey(tmp, b4T, brl)		
		plot(1,1,type='n', ylim=c(-3, -0.5), xlim=c(1, tmp[, max(table(b4T))]), ylab='log10 brl', xlab='index')
		tmp[, points( seq_along(brl), log10(brl), col=col, pch=18 ), by='b4T']		
		legend('bottomright',legend=unique(subset(tmp, select=c(b4T, col)))[, b4T], fill=unique(subset(tmp, select=c(b4T, col)))[, col], bty='n', border=NA)
		test	<- list( 	both.only.T= 	ks.test( subset(tmp, b4T=='both')[, brl], subset(tmp, b4T=='only.T')[, brl] ),
				both.only.RI= 	ks.test( subset(tmp, b4T=='both')[, brl], subset(tmp, b4T=='only.RI')[, brl] ),
				both.none= 		ks.test( subset(tmp, b4T=='both')[, brl], subset(tmp, b4T=='none')[, brl] ),
				none.only.T= 	ks.test( subset(tmp, b4T=='none')[, brl], subset(tmp, b4T=='only.T')[, brl] ),
				none.only.RI= 	ks.test( subset(tmp, b4T=='none')[, brl], subset(tmp, b4T=='only.RI')[, brl] )	)
		sapply(test, '[[', 'p.value')
		#	same distribution for 'both b4T' and 'one b4T' rejected		--> expect similar results with coin but run out of mem
		#	both.only.T 	both.only.RI    both.none  		none.only.T 	none.only.RI 
		#	1.211440e-04 	3.670217e-04	 0.000000e+00 	7.075251e-11 	9.402701e-11
		test	<- list( 	both.only.T= 	wilcox.test( subset(tmp, b4T=='both')[, log10(brl)], subset(tmp, b4T=='only.T')[, log10(brl)], alternative='greater' ),
				both.only.RI= 	wilcox.test( subset(tmp, b4T=='both')[, log10(brl)], subset(tmp, b4T=='only.RI')[, log10(brl)] ),
				both.none= 		wilcox.test( subset(tmp, b4T=='both')[, log10(brl)], subset(tmp, b4T=='none')[, log10(brl)] ),
				only.RI.T= 		wilcox.test( subset(tmp, b4T=='only.RI')[, log10(brl)], subset(tmp, b4T=='only.T')[, log10(brl)] ),
				none.only.T= 	wilcox.test( subset(tmp, b4T=='none')[, log10(brl)], subset(tmp, b4T=='only.T')[, log10(brl)] ),
				none.only.RI= 	wilcox.test( subset(tmp, b4T=='none')[, log10(brl)], subset(tmp, b4T=='only.RI')[, log10(brl)] )	)
		sapply(test, '[[', 'p.value')
		# 	both.only.T 	both.only.RI    both.none    	only.RI.T  		none.only.T 	none.only.RI 
		#	9.999306e-01 	6.530813e-04 	2.214453e-26 	7.411552e-01 	1.921079e-12 	2.952566e-12 
		tmp[, list(med.brl=median(log10(brl)), mean.brl=mean(log10(brl))), by='b4T']
		#       b4T   med.brl  mean.brl
		#1:    both -2.303041 -3.198350
		#2:    none -1.665873 -2.105176
		#3: only.RI -2.054468 -2.772468
		#4:  only.T -1.992733 -2.739838
		
		#
		#	fit Gamma and Exp and ZIExp	to data for b4T=='both'
		#
		require(fitdistrplus)
		library(ADGofTest)
		tmp			<- subset(Y.rawbrl.linked, dt>=1 & b4T=='both')		
		mle.exp		<- fitdist(tmp[,brl], 'exp')
		mle.ga		<- fitdist(tmp[,brl], 'gamma')
		mle.ln		<- fitdist(tmp[,brl], 'lnorm')
		mle.w		<- fitdist(tmp[,brl], 'weibull')
		par(mfrow=c(2, 4))
		qqcomp(mle.exp, addlegend=FALSE, main='exp dt>=1 b4T==both')
		denscomp(mle.exp, addlegend=FALSE, xlab='exp dt>=1 b4T==both')		
		qqcomp(mle.ga, addlegend=FALSE, main='gamma dt>=1 b4T==both')
		denscomp(mle.ga, addlegend=FALSE, xlab='gamma dt>=1 b4T==both')
		qqcomp(mle.ln, addlegend=FALSE, main='ln dt>=1 b4T==both')
		denscomp(mle.ln, addlegend=FALSE, xlab='ln dt>=1 b4T==both')
		qqcomp(mle.w, addlegend=FALSE, main='weibull dt>=1 b4T==both')
		denscomp(mle.w, addlegend=FALSE, xlab='weibull dt>=1 b4T==both')
		mle3			<- list()
		tmp2			<- copy(tmp)
		set(tmp2, tmp2[, which(log10(brl)<=brl.linked.min.brl)], 'brl', 0.)
		mle3$estimate	<- c(z=tmp2[, mean(brl==0)], rate=subset(tmp2, brl>0)[, 1/mean(brl)])
		test		<- list(	exp= ad.test(tmp[,brl], pexp, rate=mle.exp$estimate['rate']),								
				gamma= ad.test(tmp[,brl], pgamma, shape=mle.ga$estimate['shape'], rate=mle.ga$estimate['rate']),
				ziexp= ad.test(tmp2[,brl], pziexp, z=mle3$estimate['z'], rate=mle3$estimate['rate'])		)			
		sapply(test, '[[', 'p.value')
		#      exp.AD     gamma.AD     ziexp.AD 
		#	4.195804e-06 1.130774e-05 4.195804e-06
		tmp			<- subset(Y.rawbrl.linked, dt>=1 & b4T=='both' & log10(brl)>brl.linked.min.brl)		
		mle.exp		<- fitdist(tmp[,brl], 'exp')
		mle.ga		<- fitdist(tmp[,brl], 'gamma')
		mle.ln		<- fitdist(tmp[,brl], 'lnorm')
		mle.w		<- fitdist(tmp[,brl], 'weibull')
		par(mfrow=c(2, 4))
		qqcomp(mle.exp, addlegend=FALSE, main='exp dt>=1 b4T==both')
		denscomp(mle.exp, addlegend=FALSE, xlab='exp dt>=1 b4T==both')		
		qqcomp(mle.ga, addlegend=FALSE, main='gamma dt>=1 b4T==both')
		denscomp(mle.ga, addlegend=FALSE, xlab='gamma dt>=1 b4T==both')
		qqcomp(mle.ln, addlegend=FALSE, main='ln dt>=1 b4T==both')
		denscomp(mle.ln, addlegend=FALSE, xlab='ln dt>=1 b4T==both')
		qqcomp(mle.w, addlegend=FALSE, main='weibull dt>=1 b4T==both')
		denscomp(mle.w, addlegend=FALSE, xlab='weibull dt>=1 b4T==both')
		mle3			<- list()
		tmp2			<- copy(tmp)
		set(tmp2, tmp2[, which(log10(brl)<=brl.linked.min.brl)], 'brl', 0.)
		mle3$estimate	<- c(z=tmp2[, mean(brl==0)], rate=subset(tmp2, brl>0)[, 1/mean(brl)])
		test		<- list(	exp= ad.test(tmp[,brl], pexp, rate=mle.exp$estimate['rate']),								
				gamma= ad.test(tmp[,brl], pgamma, shape=mle.ga$estimate['shape'], rate=mle.ga$estimate['rate']),
				ziexp= ad.test(tmp2[,brl], pziexp, z=mle3$estimate['z'], rate=mle3$estimate['rate'])		)			
		sapply(test, '[[', 'p.value')
		#exp.AD  gamma.AD  ziexp.AD 
		#0.0421112 0.3320742 0.0421112
		
		tmp			<- subset(Y.rawbrl.linked, dt>=1 & b4T!='both' & b4T!='none' & log10(brl)>brl.linked.min.brl)		
		mle.exp		<- fitdist(tmp[,brl], 'exp')
		mle.ga		<- fitdist(tmp[,brl], 'gamma')
		mle.ln		<- fitdist(tmp[,brl], 'lnorm')
		mle.w		<- fitdist(tmp[,brl], 'weibull')
		par(mfrow=c(2, 4))
		qqcomp(mle.exp, addlegend=FALSE, main='exp dt>=1 b4T==one')
		denscomp(mle.exp, addlegend=FALSE, xlab='exp dt>=1 b4T==one')		
		qqcomp(mle.ga, addlegend=FALSE, main='gamma dt>=1 b4T==one')
		denscomp(mle.ga, addlegend=FALSE, xlab='gamma dt>=1 b4T==one')
		qqcomp(mle.ln, addlegend=FALSE, main='ln dt>=1 b4T==one')
		denscomp(mle.ln, addlegend=FALSE, xlab='ln dt>=1 b4T==one')
		qqcomp(mle.w, addlegend=FALSE, main='weibull dt>=1 b4T==one')
		denscomp(mle.w, addlegend=FALSE, xlab='weibull dt>=1 b4T==one')
		test		<- list(	exp= ad.test(tmp[,brl], pexp, rate=mle.exp$estimate['rate']),								
				gamma= ad.test(tmp[,brl], pgamma, shape=mle.ga$estimate['shape'], rate=mle.ga$estimate['rate'])		)			
		sapply(test, '[[', 'p.value')
		#exp.AD  	gamma.AD 
		#0.0028749 	0.2870875 
		
		#
		#	fit Gamma and Exp and ZIExp	to data excluding b4T=='none' AND small brl
		#
		tmp			<- subset(Y.rawbrl.linked, dt>=1 & brl<0.2 & b4T!='none' & log10(brl)>brl.linked.min.brl)
		par(mfrow=c(1, 4))
		mle			<- fitdist(tmp[,brl], 'exp')
		qqcomp(mle, addlegend=FALSE, main='exp dt>=1 b4T!=none brl>1e-4')
		denscomp(mle, addlegend=FALSE, xlab='exp dt>=1 b4T!=none brl>1e-4')
		mle.ga		<- fitdist(tmp[,brl], 'gamma')
		qqcomp(mle.ga, addlegend=FALSE, main='gamma dt>=1 b4T!=none brl>1e-4')
		denscomp(mle.ga, addlegend=FALSE, xlab='gamma dt>=1 b4T!=none brl>1e-4')		
		test		<- list(	exp= ad.test(tmp[,brl], pexp), 
				gamma= ad.test(tmp[,brl], pgamma, shape=mle.ga$estimate['shape'], rate=mle.ga$estimate['rate'])		)		
		sapply(test, '[[', 'p.value')
		#exp.AD     	gamma.AD 
		#0.0000122449 	0.3602207894 
		#
		#	estimate mean under Exp model for each b4T group
		#
		rawbrl.exp			<- lapply(c(0,1,2,3,4), function(dt.min)
				{
					tmp					<- subset(Y.rawbrl.linked, dt>=dt.min & log10(brl)>brl.linked.min.brl & brl<0.2)
					tmp					<- tmp[, list(exp.m= mean(brl)), by='b4T']
					#tmp					<- tmp[, list(exp.m= 1/fitdistr(brl, "exponential")$estimate), by='b4T']
					tmp[, dt:=dt.min]
					#rawbrl.l	<- subset( tmp, brl>=0 )[, sort(brl)]		#allow for one mutation
					#plot( tmp[, log10(brl)], col=tmp[,col], pch=18, main=paste('exclude dt<',dt.min))
					tmp 				
				})
		rawbrl.exp			<- do.call('rbind', rawbrl.exp )
		
		plot( tmp[, log10(brl)], col=tmp[,col], pch=18 )
		#rate estimates very similar irrespective of whether the first 1 2 3 years excluded:
		#1/rawbrl.exp	:	0.009735224 0.011242407 0.009690996 0.010396757 0.012383081
	}	
	if(substr(method,1,2)=='3a')
	{
		Y.rawbrl.linked	<- subset(Y.rawbrl.linked, b4T!='none' & log10(brl)>-12 & brl<0.2 & dt>=brl.linked.min.dt)
		#estimated mean= 0.0089
	}
	if(0)
	{
		Y.rawbrl.linked	<- subset(Y.rawbrl.linked, b4T!='none' & log10(brl)>brl.linked.min.brl & brl<0.2 & dt>=brl.linked.min.dt)
		#estimated mean= 0.0138
	}
	if(substr(method,1,2)=='3b')
	{
		Y.rawbrl.linked	<- subset(Y.rawbrl.linked, log10(brl)>brl.linked.min.brl & brl<0.2 & dt>=brl.linked.min.dt)
		#estimated mean= 0.0175
	}
	#plot(log10(Y.rawbrl.linked[, sort(brl)]	))
	rawbrl.exp		<- fitdist( Y.rawbrl.linked[, brl]	, "exp")$estimate
	cat(paste('\nestimated mean=',round(1/rawbrl.exp,d=4)))
	#rawbrl.h	<- hist(Y.rawbrl.linked[,brl], breaks=1e3, freq=0, xlim=c(0,0.05))
	#tmp			<- seq(min(Y.rawbrl.linked[,brl]), max(Y.rawbrl.linked[,brl]), length.out=1024) 
	#lines( tmp, dexp(tmp, rate = rawbrl.exp$estimate), col='red')
	Y.brl			<- Y.rawbrl
	Y.brl[, score.brl.TPp:= pgamma(Y.brl[,brl], shape=2, rate = rawbrl.exp, lower.tail=FALSE)]
	Y.brl[, score.brl.TPd:= pgamma(Y.brl[,brl], shape=2, rate = rawbrl.exp, lower.tail=FALSE)]
	Y.brl[, score.brl.TN:= p.rawbrl.unlinked(brl)]	
	#tmp			<- seq(min(Y.rawbrl.linked[,brl]), max(Y.rawbrl.linked[,brl]), length.out=1024)
	#plot(tmp, pgamma(tmp, shape=2, rate = rawbrl.exp$estimate, lower.tail=FALSE), type='l')
	#lines(tmp, pgamma(tmp, shape=1, rate = rawbrl.exp$estimate, lower.tail=FALSE), col='red')
	if(!is.na(plot.file))
	{
		pdf(plot.file, w=5, h=5)
		cat(paste('\nplot to file',plot.file))
		require(RColorBrewer)
		par(mar=c(4,4,0.5,0.5))
		cols		<- brewer.pal(4, 'Set1')[c(1,3,2)]
		legend.txt	<- c('same host', 'potential transmission pairs', 'branch length weight')	
		xlim		<- c(0, max(Y.rawbrl[, max(brl)],  Y.rawbrl.linked[, max(brl)])*1.1 )
		tmp			<- seq(from=xlim[1], to=xlim[2], by=diff(xlim)/200)
		hist( Y.rawbrl.linked[, brl], breaks=tmp , col=my.fade.col(cols[1],0.5), add=0, freq=0, xlab='branch length', main='', ylim=c(0,100) )
		hist( Y.rawbrl[, brl], breaks=tmp, col=my.fade.col(cols[2],0.5), freq=0, add=1 )
		tmp2		<- pgamma(tmp, shape=2, rate = rawbrl.exp, lower.tail=FALSE) 				
		lines(tmp, tmp2/(sum(tmp2)*diff(tmp)[1]), col=cols[3], lwd=2)		
		legend('topright', bty='n', fill= cols, border=NA, legend=legend.txt)		
		dev.off()
	}	
	Y.brl		<- subset( Y.brl, select=c(FASTASampleCode, t.FASTASampleCode, score.brl.TPd, score.brl.TPp, score.brl.TN, brl) )
	cat(paste('\nReturn brl score for #t.Patient seqs=',Y.brl[, length(unique(t.FASTASampleCode))]))
	Y.brl		
}
######################################################################################
project.athena.Fisheretal.v3.Y.rm.missedtransmitter<- function(df.tpairs, Y.brl, Y.U, Y.coal=NULL, cut.date=0.5, cut.brl=1e-3)
{
	missed					<- merge(df.tpairs, Y.brl, by=c('FASTASampleCode','t.FASTASampleCode'))
	if(!is.null(Y.coal))
		missed				<- merge(missed, Y.coal, by=c('FASTASampleCode','t.FASTASampleCode'))
	missed					<- merge(missed, Y.U, by=c('Patient','t.Patient'), allow.cartesian=TRUE)
	missed					<- subset(missed, select=c(cluster, Patient, t.Patient, FASTASampleCode, t.FASTASampleCode, t, score.brl.TN, coal.after.i.AnyPos_T1, coal.after.t.NegT, score.Inf, score.t.inf, brl))
	if(any(is.na(missed)))	stop('unexpected NA entries in missed')
	#	exclude coalescence within infected
	if(!is.null(Y.coal))
	{
		set(missed, NULL, 'coal.after.t.NegT', missed[,coal.after.t.NegT-coal.after.i.AnyPos_T1])
		set(missed, missed[,which(coal.after.t.NegT<0)], 'coal.after.t.NegT', 0.)		
	}			
	#	simple rule for screening out missed sources, using genetics and dates
	if(1)
	{
		missed[, score.t.inf:=NULL]
		setkey(missed, FASTASampleCode, t.FASTASampleCode)
		missed					<- unique(missed)
		cat(paste('\nprop of tpairs that meet the BRL.TN cutoff:', nrow(subset( missed, score.brl.TN<cut.brl ))/nrow(missed) ))
		cat(paste('\nprop of tpairs that meet the BRL.TN & COAL cutoff:', nrow(subset( missed, score.brl.TN<cut.brl &  coal.after.t.NegT>cut.date))/nrow(missed)  ))
		nmissed					<- subset( missed, score.brl.TN<cut.brl  &  coal.after.t.NegT>cut.date )		
	}
	#	simple rule for screening out missed sources and missed intermediates, using genetics and dates
	#	took this out because of case 1326, M30122->M35631: there is time to infect after coalescence in transmitter and the average is irrelevant
	if(0)
	{
		missed					<- merge(missed, missed[, list(p.nmissed= coal.after.t.NegT*score.Inf*score.t.inf), by=c('t','FASTASampleCode','t.FASTASampleCode')],by=c('t','FASTASampleCode','t.FASTASampleCode'))
		missed					<- missed[, list(p.nmissed= mean(p.nmissed), brl=brl[1], score.brl.TN=score.brl.TN[1]), by=c('FASTASampleCode','t.FASTASampleCode')]
		cat(paste('\nprop of tpairs that meet the BRL.TN cutoff:', nrow(subset( missed, score.brl.TN<cut.brl ))/nrow(missed) ))
		cat(paste('\nprop of tpairs that meet the BRL.TN & COAL cutoff:', nrow(subset( missed, score.brl.TN<cut.brl &  p.nmissed>cut.date))/nrow(missed)  ))
		nmissed					<- subset( missed, score.brl.TN<cut.brl  &  p.nmissed>cut.date )
		nmissed					<- merge( df.tpairs, nmissed, by=c('FASTASampleCode','t.FASTASampleCode') )
	}	
	#	simple rule for screening out missed sources and missed intermediates, using genetics and dates
	if(0)
	{
		missed					<- merge(missed, missed[, list(p.nmissed= coal.after.t.NegT*score.Inf*score.t.inf), by=c('t','FASTASampleCode','t.FASTASampleCode')],by=c('t','FASTASampleCode','t.FASTASampleCode'))
		missed					<- missed[, list(p.nmissed= max(p.nmissed), brl=brl[1], score.brl.TN=score.brl.TN[1]), by=c('FASTASampleCode','t.FASTASampleCode')]
		cat(paste('\nprop of tpairs that meet the BRL.TN cutoff:', nrow(subset( missed, score.brl.TN<cut.brl ))/nrow(missed) ))
		cat(paste('\nprop of tpairs that meet the BRL.TN & COAL cutoff:', nrow(subset( missed, score.brl.TN<cut.brl &  p.nmissed>cut.date))/nrow(missed)  ))
		nmissed					<- subset( missed, score.brl.TN<cut.brl  &  p.nmissed>cut.date )
		nmissed					<- merge( df.tpairs, nmissed, by=c('FASTASampleCode','t.FASTASampleCode') )
	}
	#	effect on distribution of branch lengths
	hist( subset(Y.brl, brl<0.15)[, brl], breaks=seq(0,0.15,0.002), col='blue')		
	hist( nmissed[, brl], breaks=seq(0,0.15,0.002), border=NA, col='red', add=TRUE)
	#
	nmissed	<- subset( nmissed, select=c(cluster, Patient, FASTASampleCode, t.Patient, t.FASTASampleCode))
	cat(paste('\nreturning #infecteds after screening for missed invidivuals, n=',  nmissed[, length(unique(Patient))]))
	cat(paste('\nreturning #potential transmitters after screening for missed invidivuals, n=',  nmissed[, length(unique(t.Patient))]))
	nmissed
}
######################################################################################
project.athena.Fisheretal.Y.rm.missedtransmitter<- function(YX.tpairs, df.all, Y.brl, Y.U, Y.coal=NULL, cut.date=0.5, cut.brl=1e-3, any.pos.grace.yr= 3.5, rm.zero.score=FALSE, plot.file=NA)
{
	missed					<- merge(YX.tpairs, Y.brl, by=c('FASTASampleCode','t.FASTASampleCode'), all.x=TRUE)
	if(!is.null(Y.coal))
		missed				<- merge(missed, Y.coal, by=c('FASTASampleCode','t.FASTASampleCode'), all.x=TRUE)
	missed					<- merge(missed, Y.U, by=c('Patient','t.Patient'), allow.cartesian=TRUE, all.x=TRUE)
	if(!is.null(Y.coal))
		missed				<- subset(missed, select=c(cluster, Patient, t.Patient, FASTASampleCode, t.FASTASampleCode, t, score.brl.TN, coal.after.i.AnyPos_T1, coal.after.t.NegT, score.Inf, score.t.inf, brl, class))
	if(is.null(Y.coal))
		missed				<- subset(missed, select=c(cluster, Patient, t.Patient, FASTASampleCode, t.FASTASampleCode, t, score.brl.TN, score.Inf, score.t.inf, brl, class))
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
		set(missed, missed[, which(coal.after.t.NegT<0)], 'coal.after.t.NegT', 0.)		
	}			
	#	exclude transmitters if date of diagnosis exceeds date of diagnosis of RI by more than 'any.pos.grace.yr'
	if(!is.na(any.pos.grace.yr))
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
		setkey(missed, FASTASampleCode, t.FASTASampleCode)
		missed					<- unique(missed)
		cat(paste('\nprop of tpairs that meet the BRL.TN cutoff:', nrow(subset( missed, score.brl.TN<cut.brl ))/nrow(missed) ))		
		set(missed, missed[, which(score.brl.TN>=cut.brl)], 'score.Y', 0.)
		if(!is.null(Y.coal))
		{
			cat(paste('\nprop of tpairs that meet the BRL.TN & COAL cutoff:', nrow(subset( missed, score.brl.TN<cut.brl &  coal.after.t.NegT>cut.date))/nrow(missed)  ))
			set(missed, missed[, which(!is.na(coal.after.t.NegT) & coal.after.t.NegT<=cut.date)], 'score.Y', 0.)
			cat(paste('\nnumber of pn pairs with score.Y>0:',missed[, length(which(is.na(coal.after.t.NegT) & class=='pn' & score.Y>0))]))			
		}						
	}
	#	effect on distribution of branch lengths
	if(!is.na(plot.file))
	{
		pdf(file=plot.file, 5, 5)
		hist( subset(Y.brl, brl<0.15)[, brl], breaks=seq(0,0.15,0.002), col='blue')		
		hist( subset(missed, score.Y>0)[, brl], breaks=seq(0,0.15,0.002), border=NA, col='red', add=TRUE)
		dev.off()		
		pdf(file=paste( substr(plot.file,1,nchar(plot.file)-4),'ptvspn.pdf',sep='' ), 5, 5)
		hist( subset(missed, class=='pt')[, brl], breaks=seq(0,0.5,0.002), col='green')		
		hist( subset(missed, class=='pn')[, brl], breaks=seq(0,0.5,0.002), border=NA, col='red', add=TRUE)
		legend('topright',bty='n',border=NA,fill=c('green','red'),legend=c('pt','pn'))
		dev.off()
	}
	#
	missed		<- subset( missed, select=c(cluster, Patient, FASTASampleCode, t.Patient, t.FASTASampleCode, score.Y, class))
	cat(paste('\n#infecteds with score.Y>0 after screening for missed invidivuals, n=',  subset(missed, score.Y>0)[, length(unique(Patient))]))
	cat(paste('\n#potential with score.Y>0 transmitters after screening for missed invidivuals, n=',  subset(missed, score.Y>0)[, length(unique(t.Patient))]))
	if(rm.zero.score)
		missed	<- subset(missed, score.Y>0)
	missed
}
######################################################################################
project.athena.Fisheretal.v2.Y.weight<- function(YX)
{
	#	add tpair weight: every transmitter can infect only in one time period
	YX		<- merge(YX, YX[, list(w= 1/length(t)),by=c('t.FASTASampleCode','FASTASampleCode')], by=c('t.FASTASampleCode','FASTASampleCode'))
	#	add infected weight: every infected can only be infected by one transmitter
	tmp		<- subset(YX, select=c(Patient, t.Patient))
	setkey(tmp, Patient, t.Patient)
	tmp		<- unique(tmp)	
	YX		<- merge(YX, tmp[, list(w.i= length(unique(t.Patient))),by=c('Patient')], by='Patient') 
	set(YX, NULL, 'w', YX[, w/w.i])
	YX[, w.i:=NULL]
	YX
}
######################################################################################
project.athena.Fisheretal.YX.weight<- function(YX)
{	
	#	add tpair weight: every transmitter can infect only in one time period
	YX		<- merge(YX, YX[, list(w= 1/length(t)),by=c('t.FASTASampleCode','FASTASampleCode')], by=c('t.FASTASampleCode','FASTASampleCode'))
	#	check weight
	if( abs(nrow(unique( subset(YX, select=c(Patient, t.Patient)) ))-YX[, sum(w)])>5*EPS )	stop('unexpected weight')	
	print( YX[, table(w)] )
	#	add infected weight: every infected can only be infected by one transmitter
	tmp		<- subset(YX, select=c(Patient, t.Patient, score.brl.TPd))
	setkey(tmp, Patient, t.Patient)
	tmp		<- unique(tmp)	
	tmp		<- tmp[,	list(w.i=score.brl.TPd/sum(score.brl.TPd), t.Patient=t.Patient), by='Patient']	
	if( tmp[,sum(w.i)]!=YX[, length(unique(Patient))] )	stop('unexpected weight')	
	YX		<- merge(YX, tmp, by=c('Patient','t.Patient')) 
	set(YX, NULL, 'w', YX[, w*w.i])
	if( abs(YX[,sum(w)]-YX[, length(unique(Patient))])>5*EPS )	stop('unexpected weight')
	#	weights are now normalised per 'Patient' but there are fewer infection EVENTS because Patient may also be a transmitter to t.Patient
	triplet.weight				<- subset(YX, select=c(Patient, t.Patient))	
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
	YX		<- merge(YX, subset(triplet.weight, select=c(Patient, t.Patient, w.t)), by=c('Patient','t.Patient'))
	set(YX, NULL, 'w', YX[, w*w.t])
	YX[,score.brl.TPd:=NULL]
	#	
	YX	
}
######################################################################################
project.athena.Fisheretal.v3.Y.transmitterinfected<- function(Y.infwindow, X.pt)
{
	tmp			<- subset(X.pt, select=c(Patient, t.Patient, t, U.score))
	setkey(tmp, Patient, t.Patient, t)
	tmp			<- unique(tmp)		
	Y.U 		<- merge( subset(Y.infwindow, select=c(Patient, t, score.Inf)), tmp, by=c('Patient','t'))
	set(Y.U, Y.U[, which(is.na(U.score))], 'U.score', 1)
	setnames(Y.U, 'U.score', 'score.t.inf')
	Y.U
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
		
		tmp			<- msm.unlinked.bytime[, quantile(brl, p=c(1e-5, 1e-4, 5e-4, 0.001, 0.01))]
		ltys		<- seq_along(tmp)+1
		abline(v=tmp, col=cols[2], lty=ltys)
		legend('bottomright', bty='n', lty= c(ltys,NA), col=cols[2], border=NA, legend=c('1e-5', '1e-4', '5e-4', '1e-3', '1e-2',''))
		dev.off()
	}
	
	list(tpairs=df.tpairs.brl, linked=msm.linked, unlinked=msm.unlinked.bytime)	
}
######################################################################################
project.athena.Fisheretal.select.denominator<- function(indir, infile, insignat, indircov, infilecov, infile.viro, infile.immu, infile.treatment, infiletree=NULL, adjust.AcuteByNegT=NA, adjust.NegT4Acute=NA, adjust.NegTByDetectability=NA, adjust.minSCwindow=NA, adjust.AcuteSelect=c('Yes'), t.recent.endctime=2013.)
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
	#
	#
	msm.recent		<- subset( df.all, Sex=='M' & !Trm%in%c('OTH','IDU','HET','BLOOD','PREG','HETfa','NEEACC','SXCH') )
	msm.recent		<- subset(msm.recent, AnyPos_T1<t.recent.endctime)
	cat(paste('\nmsm before endctime',t.recent.endctime,': #seq=',nrow(msm.recent),'#patient=',length(msm.recent[,unique(Patient)])))
	setkey(msm.recent, isAcute)
	msm.recent		<- msm.recent[adjust.AcuteSelect,]
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
		clumsm.recent	<- clumsm.info[adjust.AcuteSelect,] 
		clumsm.recent	<- subset( clumsm.recent, Trm%in%c('MSM','BI') )
		set(clumsm.recent, NULL, 'Trm', clumsm.recent[,factor(as.character(Trm))])
		clumsm.recent		<- subset(clumsm.recent, AnyPos_T1<t.recent.endctime)
		cat(paste('\nmsm clustering before endctime',t.recent.endctime,': #seq=',nrow(clumsm.recent),'#patient=',length(clumsm.recent[,unique(Patient)])))		
		# 	recent seroconcerverters in cluster
		clumsm.recentsc	<- subset( clumsm.recent, !is.na(NegT))		
		cat(paste('\nmsm clustering: #seq=',nrow(clumsm.info),'#patient=',length(clumsm.info[,unique(Patient)]),'#cluster=',length(clumsm.info[,unique(cluster)])))		
		cat(paste('\nmsm clustering isAcute: #seq=',nrow(clumsm.recent),'#patient=',length(clumsm.recent[,unique(Patient)]),'#cluster=',length(clumsm.recent[,unique(cluster)])))
		cat(paste('\nmsm clustering isAcute & !is.na(NegT): #seq=',nrow(clumsm.recentsc),'#patient=',length(clumsm.recentsc[,unique(Patient)]),'#cluster=',length(clumsm.recentsc[,unique(cluster)])))
		
		ans<- list(df.all=df.all, df.viro=df.viro, df.immu=df.immu, df.treatment=df.treatment, clumsm.info=clumsm.info, df.select=clumsm.recent, clumsm.subtrees=msm$cluphy.subtrees, clumsm.ph=msm$cluphy)
	}
	else
	{
		ans<- list(df.all=df.all, df.viro=df.viro, df.immu=df.immu, df.treatment=df.treatment, df.select=msm.recent )
	}
	ans	
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
project.athena.Fisheretal.t2inf<- function(indircov, infilecov, adjust.AcuteByNegT=0.75, adjust.dt.CD4=1, adjust.AnyPos_y=2003, adjust.NegT=2)
{
	require(MASS)
	
	load(paste(indircov,'/',infilecov,'.R',sep=''))		
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
					mp.NegT=round(as.numeric(difftime(AnyPos_T1, NegT, units='days'))/adjust.NegT),
					dt.CD4=as.numeric(difftime(PosCD4_T1, AnyPos_T1, units='days'))/365,
					AnyPos_a=as.numeric(difftime(AnyPos_T1, DateBorn, units='days'))/365,
					AnyPos_y=hivc.db.Date2numeric(AnyPos_T1)), by='Patient']
	
	df.sc.negT		<- merge(df.sc.negT, tmp, by='Patient')
	df.scchr.negT	<- subset( df.sc.negT, isAcute=='No' & AnyPos_y>adjust.AnyPos_y)
	df.scchr.cd4	<- subset( df.scchr.negT, dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)
	
	
	m4a	<- glm.nb(mp.NegT ~  1, data=subset(df.scchr.cd4, CD4_T1>850 & dt.NegT<4), link=identity, trace= 1, maxit= 50)	
	m4d	<- glm.nb(mp.NegT ~  AnyPos_a, data=subset(df.scchr.negT, AnyPos_a<=45), link=identity, trace= 1, maxit= 50)		
	m4b	<- glm.nb(mp.NegT ~  AnyPos_a+0+offset(rep(m4d$coefficients["(Intercept)"],length(mp.NegT))), data=subset(df.scchr.cd4, CD4_T1<=850 & CD4_T1>250 & AnyPos_a<=45), link=identity, trace= 1, maxit= 50)
	m4c	<- glm.nb(mp.NegT ~  AnyPos_a+0+offset(rep(m4d$coefficients["(Intercept)"],length(mp.NegT))), data=subset(df.scchr.cd4, CD4_T1<250 & AnyPos_a<=45), link=identity, trace= 1, maxit= 50)
	
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
	
	#		
	predict.t2inf	<- function(q, df, t2inf.args, t2inf.method='glm.nb1')
	{		
		if(t2inf.method!='glm.nb1')	stop('unknown method to predict t2inf')
		df[, {
					if(!is.na(isAcute) & isAcute=='Yes')
					{
						ans	<- pexp(q, 2/365, lower.tail=FALSE)
					}
					else if(!is.na(isAcute) & isAcute=='Maybe')
					{
						ans	<- pexp(q, 1/320, lower.tail=FALSE)
					}
					else if(is.na(PosCD4_T1) | (!is.na(AnyT_T1) & (AnyT_T1<PosCD4_T1 | PosCD4_T1-AnyPos_T1<1)))		#chronic or missing isAcute, first CD4 after ART start or first CD4 too far from PosDiag
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
					list(q=q, score=ans)
				}, by='Patient']			
	}
	
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
project.athena.Fisheretal.plot.selected.transmitters<- function(df.all, df.immu, df.viro, df.treatment, df.tpairs, cluphy, cluphy.info, cluphy.subtrees, cluphy.map.nodectime, file, pdf.height=400)
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
	cat(paste('\nplot to file=',file))
	dummy				<- hivc.beast2out.plot.cluster.trees(df.all, df.immu, df.viro, df.treatment, cluphy, cluphy.root.ctime, cluphy.tip.ctime, ph.prob=NA, df.node.ctime=cluphy.map.nodectime, df.rates=NULL, df.tips=df.tpairs.plot, end.ctime=end.ctime,  cex.nodelabel=0.5,  cex.tiplabel=0.5,  file=file,  pdf.width=7, pdf.height=pdf.height, pdf.xlim=pdf.xlim)	
}
######################################################################################
project.athena.Fisheretal.X.b4care<- function(df.tpairs, clumsm.info, predict.t2inf, t2inf.args, t.period=0.25, ts.min=1980)
{
	b4care	<- merge(data.table(Patient=df.tpairs[, unique(t.Patient)]), unique(subset(clumsm.info, select=c(Patient, DateBorn, AnyPos_T1, isAcute, NegT, PosCD4_T1, CD4_T1, AnyT_T1))), by='Patient')
	tmp		<- b4care[, list(AnyPos_a=AnyPos_T1-DateBorn, ts=ifelse(is.na(NegT), AnyPos_T1-10, NegT), te=AnyPos_T1), by='Patient']
	b4care	<- merge(b4care, tmp, by='Patient')
	#check if ts before 1980 and if so clip
	set(b4care, b4care[,which(ts<ts.min)], 'ts', ts.min )	
	set(b4care, NULL, 'ts', b4care[, floor(ts) + round( (ts%%1)*100 %/% (t.period*100) ) * t.period ] )
	set(b4care, NULL, 'te', b4care[, floor(te) + round( (te%%1)*100 %/% (t.period*100) ) * t.period - t.period] )	#last time undiagnosed is just before first time period diagnosed
	#apply predict.t2inf	
	t2inf	<- predict.t2inf(seq(0,10,t.period)*365, b4care, t2inf.args)
	t2inf	<- merge( subset( t2inf, score>0.01 ), subset(b4care, select=c(Patient,ts,te)), by='Patient' )
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
project.athena.Fisheretal.X.incare<- function(df.tpairs, clumsm.info, df.viro, df.immu, df.treatment, t.period=0.25, t.endctime= 2013.)
{
	#	prepare incare timeline for potential transmitters
	incare		<- merge( data.table(Patient=df.tpairs[, unique(t.Patient)]), unique( subset(clumsm.info, select=c(Patient, AnyPos_T1, DateDied)) ), by='Patient' )
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
	cat(paste('\nincare.t entries, n=',nrow(incare.t)))
	#	prepare treatment variables for potential transmitters
	treat		<- subset(df.treatment, select=c(Patient, AnyT_T1, StartTime, StopTime, TrI, TrCh.failure, TrCh.adherence, TrCh.patrel, NoDrug, NoNRT, NoNNRT, NoPI)) 
	treat		<- merge( data.table(Patient=df.tpairs[, unique(t.Patient)]), treat, by='Patient' )
	set(treat, NULL, 'AnyT_T1', hivc.db.Date2numeric(treat[,AnyT_T1]))
	set(treat, NULL, 'StopTime', hivc.db.Date2numeric(treat[,StopTime]))
	set(treat, NULL, 'StartTime', hivc.db.Date2numeric(treat[,StartTime]))
	set(treat, NULL, 'StartTime', treat[, floor(StartTime) + round( (StartTime%%1)*100 %/% (t.period*100) ) * t.period] )
	set(treat, NULL, 'StopTime', treat[, floor(StopTime) + round( (StopTime%%1)*100 %/% (t.period*100) ) * t.period] )
	treat		<- subset(treat, StartTime!=StopTime)	#delete too small Start to Stop time intervals	-- we can have interrupted at start
	cat(paste('\nTreatment info available for pot transmitters, n=',treat[, length(unique(Patient))]))
	#	remove first ART periods that begin with an interruption
	treat		<- merge( treat, subset( treat[, list(select=!all(TrI=='Yes')), by='Patient'], select, Patient), by='Patient' )
	cat(paste('\nTreatment info for pot transmitters for which not all periods are ART interrupted, n=',treat[, length(unique(Patient))]))
	treat		<- treat[, {
								tmp<- seq.int(which(TrI!='Yes')[1], length(TrI))
								lapply(.SD,'[',tmp)
							},by='Patient']	
	#	get treatment timeline
	treat.t		<- merge(subset(incare.t, select=c(Patient, t)), treat, by='Patient', allow.cartesian=TRUE )	
	treat.t		<- subset( treat.t, StartTime<=t+t.period/2 & t+t.period/2<StopTime )
	cat(paste('\ntreat.t entries, n=',nrow(treat.t)))
	setnames(treat.t, c('TrI','TrCh.failure','TrCh.adherence','TrCh.patrel','NoDrug','NoNRT','NoNNRT','NoPI'), c('ART.I','ART.F','ART.A','ART.P','ART.nDrug','ART.nNRT','ART.nNNRT','ART.nPI'))
	#	merge incare and treatment timelines
	incare.t	<- merge(incare.t, treat.t, all.x=1, by=c('Patient','t'))	
	#	set ever on ART per period t		(avoid AnyT_T1 as it may give periods that would start with treatment interruption)
	set(incare.t, incare.t[, which(!is.na(StartTime))], 'stage', 'ART.started')
	cat(paste("\nNumber entries with ART.I=='Yes' & stage=='Diag', n=",nrow(subset(incare.t, ART.I=='Yes' & stage=='Diag'))))
	#	compute X: viral load of potential transmitter  	
	X.viro		<- project.athena.Fisheretal.X.viro(df.tpairs, df.viro, t.period=t.period, lRNA.cutoff=NA)
	cat(paste('\nVL info available for pot transmitters, n=',X.viro[, length(unique(t.Patient))]))
	#	compute X: CD4 of potential transmitter
	X.cd4		<- project.athena.Fisheretal.X.cd4(df.tpairs, df.immu, t.period=t.period)
	cat(paste('\nCD4 info available for pot transmitters, n=',X.cd4[, length(unique(t.Patient))]))
	#	add viro and cd4 time periods to incare.t
	setnames(incare.t, 'Patient','t.Patient')		
	incare.t	<- merge(incare.t, X.viro, by=c('t.Patient','t'), all.x=1)
	incare.t	<- merge(incare.t, X.cd4, by=c('t.Patient','t'), all.x=1)
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
project.athena.Fisheretal.YX.model4.v2<- function(YX, clumsm.info, vl.suppressed=log10(1e3), acute.select=c('Yes','Maybe'), cd4.cut= c(-1, 350, 550, 5000) )
{
	require(betareg)
	#	Diagnosis risk groups
	#	
	#	PQ: Are CD4 counts at time period or at diagnosis associated with higher transmission probability?		 diagnosis risk groups: age at diagnosis, t2care, t2vlsupp, median follow up?
	#	PQ: Is Acute associated with higher transmission probability? Is high viral load associated with higher transmission probability? 
	#	SQ:	Are the younger ages associated with high transmission prob?
	#		--> this is a more general question, not just after diagnosis
	#	SQ:	Is time to care or time to suppression associated with high transmission prob?
	#	SQ: Are late presenters associated with high transmission prob? (CD4<=220, CDCC at diag, Age at diag > 45)
	#		--> this could also include time before diagnosis
	
	YX		<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(YX.m4, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA, plot.file.or=NA )
	YX.m4	<- copy(YX)
	set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(CD4t))])
	# collect PoslRNA_TL and lRNA_TL
	tmp		<- merge( unique(subset(clumsm.info, select=Patient)), subset(df.viro, select=c(Patient, PosRNA, lRNA)), by='Patient' )
	set( tmp, NULL, 'PosRNA', hivc.db.Date2numeric(tmp[,PosRNA]) )
	tmp		<- subset( tmp[, 	{
						tmp<- which(!is.na(lRNA))
						list(t.Patient=Patient[1], lRNA_TL=lRNA[tail(tmp,1)], PoslRNA_TL=PosRNA[tail(tmp,1)])
					}, by='Patient'], select=c(t.Patient, lRNA_TL, PoslRNA_TL) )
	YX.m4	<- merge(YX.m4, tmp, by= 't.Patient', all.x=TRUE)
	#define lRNA.mxw and lRNA.gm
	YX.m4	<- merge(YX.m4, YX.m4[, {
										tmp<- which(!is.na(lRNA))					
										if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1])
											lRNA.gm		<- lRNA.mxw	<- lRNA_TL[1]						
										else if(!length(tmp))
											lRNA.gm		<- lRNA.mxw	<- NA_real_
										else
										{
											lRNA.gm		<- mean(lRNA[tmp])
											lRNA.mxw	<- max(lRNA[tmp])
										}
										list( lRNA.gm= lRNA.gm, lRNA.mxw= lRNA.mxw )
									}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))		
	#	focus on diagnosis stage
	YX.m4s	<- subset(YX.m4, !stage%in%c('ART.suA.N','ART.suA.Y','ART.vlNA','U','UAy','UAm'))
	set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
	YX.m4s[, table(stage)]				
	#
	YX.m4.fit11.t 	<- betareg(score.Y ~ bs(lRNA, knots=c(3.00001,4.5), degree=1)+stage-1, link='logit', weights=w, data = YX.m4s)
	YX.m4.fit11.gm 	<- betareg(score.Y ~ bs(lRNA.gm, knots=c(3.00001,4.5), degree=1)+stage-1, link='logit', weights=w, data = YX.m4s)
	YX.m4.fit11.mxw <- betareg(score.Y ~ bs(lRNA.mxw, knots=c(3.00001,4.5), degree=1)+stage-1, link='logit', weights=w, data = YX.m4s)
	
	sapply( list(YX.m4.fit11.t, YX.m4.fit11.mxw, YX.m4.fit11.gm), '[[', 'pseudo.r.squared' ) 
	#	0.06925479 0.07393748 0.08227272
	AIC(YX.m4.fit11.t, YX.m4.fit11.mxw, YX.m4.fit11.gm)
	#                df       AIC
	#YX.m4.fit11.t   10 -210.5019
	#YX.m4.fit11.mxw 10 -242.1103
	#YX.m4.fit11.gm  10 -242.8205
	
	
}
######################################################################################
project.athena.Fisheretal.YX.model4.v1<- function(YX, clumsm.info, vl.suppressed=log10(1e3), acute.select=c('Yes','Maybe'), cd4.cut= c(-1, 350, 550, 5000) )
{
	require(betareg)
	#	Diagnosis risk groups
	#	
	#	PQ: Are CD4 counts at time period or at diagnosis associated with higher transmission probability?		 diagnosis risk groups: age at diagnosis, t2care, t2vlsupp, median follow up?
	#	PQ: Is Acute associated with higher transmission probability? Is high viral load associated with higher transmission probability? 
	#	SQ:	Are the younger ages associated with high transmission prob?
	#		--> this is a more general question, not just after diagnosis
	#	SQ:	Is time to care or time to suppression associated with high transmission prob?
	#	SQ: Are late presenters associated with high transmission prob? (CD4<=220, CDCC at diag, Age at diag > 45)
	#		--> this could also include time before diagnosis
	
	YX.m4	<- copy(YX)
	YX.m4[, U.score:=NULL]
	#	if transmitter acute, set undiagnosed to undiagnosed and acute
	set(YX.m4, YX.m4[, which(stage=='U' & t.isAcute=='Yes')], 'stage', 'UAy')
	set(YX.m4, YX.m4[, which(stage=='U' & t.isAcute=='Maybe')], 'stage', 'UAm')
	#
	#	stratify ART into suppressed during infection window and not suppressed during infection window
	#	minimize NAs by using extrapolation for 6 months after last known VL
	#
	VL.cur	<- 3
	tmp		<- unique(subset( clumsm.info, select=c(Patient, PoslRNA_T1) ))
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	YX.m4	<- merge(YX.m4, tmp, by= 't.Patient')	
	# collect PoslRNA_TL and lRNA_TL
	tmp		<- merge( unique(subset(clumsm.info, select=Patient)), subset(df.viro, select=c(Patient, PosRNA, lRNA)), by='Patient' )
	set( tmp, NULL, 'PosRNA', hivc.db.Date2numeric(tmp[,PosRNA]) )
	tmp		<- subset( tmp[, 	{
						tmp<- which(!is.na(lRNA))
						list(t.Patient=Patient[1], lRNA_TL=lRNA[tail(tmp,1)], PoslRNA_TL=PosRNA[tail(tmp,1)])
					}, by='Patient'], select=c(t.Patient, lRNA_TL, PoslRNA_TL) )
	YX.m4	<- merge(YX.m4, tmp, by= 't.Patient', all.x=TRUE)
	set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
	tmp		<- YX.m4[, {
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
	YX.m4	<- merge(YX.m4, tmp, by=c('Patient','t.Patient'))	
	set(YX.m4, YX.m4[,which(stage=='ART.started' & lRNA.c=='SuA.Y')],'stage', 'ART.suA.Y' )
	set(YX.m4, YX.m4[,which(stage=='ART.started' & lRNA.c=='SuA.N')],'stage', 'ART.suA.N' )
	set(YX.m4, YX.m4[,which(stage=='ART.started')],'stage', 'ART.vlNA' )	
	set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])		
					
	YX.m4[, lRNA.c:= NULL]			
	YX.m4[, stage.orig:= stage]
	YX.m4[, table(stage)]
	#		base case: stratify Diag by first CD4 after diagnosis
	if(0)
	{
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		cd4.label	<- c('D1<=350','D1<=550','D1>550')
		cd4.cut		<- c(-1, 350, 550, 5000)
		tmp			<- subset(clumsm.info, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
		setkey(tmp, Patient)
		tmp			<- unique(tmp)	
		tmp[, CD4.c:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
		cat(paste('\ntransmitters with no CD4 after diagnosis\n'))
		print( subset( tmp, is.na(CD4.c) ) )
		setnames(tmp, 'Patient', 't.Patient')
		set(tmp, tmp[, which(is.na(CD4.c))], 'CD4.c', 'D1.NA')
		YX.m4		<- merge(YX.m4, subset(tmp, select=c(t.Patient, CD4.c)), by='t.Patient')
		tmp			<- YX.m4[, which(stage=='Diag')]
		set(YX.m4, tmp, 'stage', YX.m4[tmp, CD4.c])		
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
		YX.m4[, table(stage)]
		
		YX.m4.fit1 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
		YX.m4.fit1.or	<- data.table( 	D1l350= my.or.from.logit(YX.m4.fit1, 'stageD1<=350', 'stageD1>550', subset(YX.m4, stage=='D1<=350')[, sum(w)], subset(YX.m4, stage=='D1>550')[, sum(w)], 1.962),
										D1l550= my.or.from.logit(YX.m4.fit1, 'stageD1<=550', 'stageD1>550', subset(YX.m4, stage=='D1<=550')[, sum(w)], subset(YX.m4, stage=='D1>550')[, sum(w)], 1.962)	)
		# YX.m4.fit1.or	(not adjusted for acute)
      	#	D1l350    D1l550    D1g550
		#1: 1.3039782 0.9266223 0.6534647
		#2: 0.9731111 0.7166847 0.5337555
		#3: 1.7473434 1.1980566 0.8000220
		#	no significant difference between diagnosis groups						
		YX.m4.p			<- YX.m4[, list(n=round(sum(w),d=1),prop= round(sum(w)/YX.m4[,sum(w)],d=3) ),by='stage']
		#
		#	adjusted for acute:
		#
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAm' )				
		tmp			<- YX.m4[, which(stage=='Diag')]
		set(YX.m4, tmp, 'stage', YX.m4[tmp, CD4.c])		
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
		YX.m4[, table(stage)]
		#	ART.suA.N ART.suA.Y  ART.vlNA   D1<=350   D1<=550    D1>550       DAm       DAy         U       UAm       UAy 
     	#	3172      4526        70       	807       2590       3322         702       851      7606      1275      1291 
		YX.m4.fit1 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
		YX.m4.fit1.or	<- data.table( 	D1l350= my.or.from.logit(YX.m4.fit1, 'stageD1<=350', 'stageD1>550', subset(YX.m4, stage=='D1<=350')[, sum(w)], subset(YX.m4, stage=='D1>550')[, sum(w)], 1.962),
				D1l550= my.or.from.logit(YX.m4.fit1, 'stageD1<=550', 'stageD1>550', subset(YX.m4, stage=='D1<=550')[, sum(w)], subset(YX.m4, stage=='D1>550')[, sum(w)], 1.962),
				D1g550= my.or.from.logit(YX.m4.fit1, 'stageD1>550', 'stageU', subset(YX.m4, stage=='D1>550')[, sum(w)], subset(YX.m4, stage=='U')[, sum(w)], 1.962) 		)
		#D1l350    D1l550
		#1: 1.0851541 0.8530707
		#2: 0.7779292 0.6410546
		#3: 1.5137102 1.1352070
		#	plot
		tmp			<- data.table(stage= YX.m4[, levels(stage)], col=sapply( rainbow(YX.m4[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m4		<- merge(YX.m4, tmp, by='stage')
		plot( seq_len(nrow(YX.m4)), YX.m4[, score.Y], pch=18, col= YX.m4[, col], cex=YX.m4[, w^0.4])
		tmp				<- subset(YX.m4, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m4.fit1, tmp, type='response'))	
		start			<- c( YX.m4.fit1$coef$mean + head( 2*sqrt(diag(vcov(YX.m4.fit1))), length(YX.m4.fit1$coef$mean)), YX.m4.fit1$coef$precision)
		YX.m4.fit1.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m4, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m4.fit1$coef$mean - head( 2*sqrt(diag(vcov(YX.m4.fit1))), length(YX.m4.fit1$coef$mean)), YX.m4.fit1$coef$precision)
		YX.m4.fit1.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m4, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
				c( predict(YX.m4.fit1.sup, tmp, type='response'), rev(predict(YX.m4.fit1.slw, tmp, type='response'))), 
				border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m4[, levels(stage)], fill=YX.m4[, unique(col)])
		YX.m4[, col:=NULL]
		#

		YX.m4[, CD4.c:=NULL]
	}
	#		try 	Diag<350, Diag<550, Diag>550 a better explanation?
	if(0)
	{
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAm' )						
		cd4.label	<- c('Dt<=350','Dt<=550','Dt>550')
		set(YX.m4, YX.m4[, which(stage=='Diag' & is.na(CD4))], 'stage', 'Diag.NA')
		for(i in seq_along(cd4.cut)[-1])
		{
			tmp	<- YX.m4[, which(stage=='Diag' & !is.na(CD4) & CD4<=cd4.cut[i] & CD4>cd4.cut[i-1])]
			cat(paste('\nnumber entries with Diag and CD4 ',cd4.cut[i-1], cd4.cut[i],' , n=', length(tmp)))
			set(YX.m4, tmp, 'stage', cd4.label[i-1])
		}
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
		YX.m4[, table(stage)]
		#ART.suA.N ART.suA.Y  ART.vlNA       DAm       	DAy   Diag.NA   Dt<=350   Dt<=550    Dt>550         U       UAm       UAy 
     	#3172      4526        70       	702       	851       766       964      3010      1979      7606      1275      1291 
		YX.m4.fit2 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
		YX.m4.fit2.or	<- data.table( 	D1l350= my.or.from.logit(YX.m4.fit2, 'stageDt<=350', 'stageDt>550', subset(YX.m4, stage=='Dt<=350')[, sum(w)], subset(YX.m4, stage=='Dt>550')[, sum(w)], 1.962),
										D1l550= my.or.from.logit(YX.m4.fit2, 'stageDt<=550', 'stageDt>550', subset(YX.m4, stage=='Dt<=550')[, sum(w)], subset(YX.m4, stage=='Dt>550')[, sum(w)], 1.962),
										D1g550= my.or.from.logit(YX.m4.fit2, 'stageDt>550', 'stageU', subset(YX.m4, stage=='Dt>550')[, sum(w)], subset(YX.m4, stage=='U')[, sum(w)], 1.962) 		)
		#> YX.m4.fit2.or	(not adjusted for acute)
      	#		D1l350    D1l550    
		#1: 	1.0229417 1.0216966 
		#2: 	0.7195683 0.7703738 
		#3: 	1.4542187 1.3550095 
		#D1l350    D1l550    (adjusted for acute)
		#1: 1.0200557 0.9414970 
		#2: 0.6847272 0.6844514 
		#3: 1.5196032 1.2950762 
		# no significant difference between diagnosis groups
	}
	#		late presenters: stratify Diag by first CD4 after diagnosis
	if(1)
	{
		cd4.cut		<- c(-1, 220, 350, 1e4)
		cd4.label	<- c('D1<=220','D1<=350','D1>350')
		YX.m4[, CD4.c:=NULL]
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAm' )		
		#
		tmp			<- subset(clumsm.info, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
		setkey(tmp, Patient)
		tmp			<- unique(tmp)	
		tmp[, CD4.c:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]		
		setnames(tmp, 'Patient', 't.Patient')
		set(tmp, tmp[, which(is.na(CD4.c))], 'CD4.c', 'D1.NA')
		YX.m4		<- merge(YX.m4, subset(tmp, select=c(t.Patient, CD4.c)), by='t.Patient')
		tmp			<- YX.m4[, which(stage=='Diag')]
		set(YX.m4, tmp, 'stage', YX.m4[tmp, CD4.c])
		YX.m4[, CD4.c:=NULL]
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
		YX.m4[, table(stage)]
		#ART.suA.N ART.suA.Y  ART.vlNA   	D1<=220   D1<=350    	D1>350       	DAm       DAy         U       UAm       UAy 
     	#3172      4526        70       	122       685      		5912       		702       851      7606      1275      1291
		YX.m4.fit12 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
		YX.m4.fit12.or		<- data.table( 	D1l220= my.or.from.logit(YX.m4.fit12, 'stageD1<=220', 'stageD1>350', subset(YX.m4, stage=='D1<=220')[, sum(w)], subset(YX.m4, stage=='D1>350')[, sum(w)], 1.962),
											D1l350= my.or.from.logit(YX.m4.fit12, 'stageD1<=350', 'stageD1>350', subset(YX.m4, stage=='D1<=350')[, sum(w)], subset(YX.m4, stage=='D1>350')[, sum(w)], 1.962))		
		#	D1l220   	D1l350
		#1: 0.6272600 	1.417792
		#2: 0.4792174 	1.086661
		#3: 0.8210368 	1.849825
	}
	#		try 	Age at diagnosis
	if(1)
	{
		YX.m4[, t.AnyPos_A:=NULL]
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAm' )				
		tmp		<- subset( clumsm.info, select=c(Patient, DateBorn, AnyPos_T1) )
		setkey(tmp, Patient)
		tmp		<- unique(tmp)
		tmp[, AnyPos_A:= AnyPos_T1-DateBorn]
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		tmp		<- merge( unique( subset(YX.m4, select=t.Patient) ), tmp, by='t.Patient' )		
		YX.m4	<- merge(YX.m4, subset(tmp, select=c(t.Patient, t.AnyPos_A)), by='t.Patient')
		
		age.cut		<- c(-1, 20, 45, 100)
		age.label	<- c('D1<=20','D1<=45','D1>45')
		for(i in seq_along(cd4.cut)[-1])
			set(YX.m4, YX.m4[, which(stage=='Diag' & t.AnyPos_A<=age.cut[i] & t.AnyPos_A>age.cut[i-1])], 'stage', age.label[i-1])
		YX.m4[, t.AnyPos_A:=NULL]
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
		YX.m4[, table(stage)]
		#	ART.suA.N ART.suA.Y  	ART.vlNA    D1<=25   D1<=45    D1>45        U       	UAm       UAy 
     	#	3172      4526        	70       	234      6652      1386      	7606      	1275      1291	
		#	ART.suA.N ART.suA.Y  ART.vlNA    	D1<=20   D1<=45    D1>45       DAm       DAy         U       UAm       UAy 
     	#	3172      4526        70       		155      5436      1128        702       851      7606      1275      1291 
		YX.m4.fit3 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
		YX.m4.fit3.or	<- data.table( 	D1l20= my.or.from.logit(YX.m4.fit3, 'stageD1<=20', 'stageD1<=45', subset(YX.m4, stage=='D1<=20')[, sum(w)], subset(YX.m4, stage=='D1<=45')[, sum(w)], 1.962),
										D1g45= my.or.from.logit(YX.m4.fit3, 'stageD1>45', 'stageD1<=45', subset(YX.m4, stage=='D1>45')[, sum(w)], subset(YX.m4, stage=='D1<=45')[, sum(w)], 1.962)
										)		
		#      D1l20    D1g45
		#1: 2.145616 1.341343
		#2: 1.707319 1.061021
		#3: 2.696432 1.695725
		#      D1l20     D1g45	(adjusting for acute)
		#1: 2.721543 1.1584342
		#2: 2.136278 0.8892229
		#3: 3.467149 1.5091489
		#	young and old age at diagnosis significant; D1l25 was not significant
		#	young age at diagnosis remains significant
	}
	#		try 	mean age in infection window
	if(1)
	{		
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])		
		#YX.m4[, list(t.Age.m= mean(t.Age)),by=c('t.Patient','Patient')]
		age.cut		<- c(-1, 23, 45, 100)
		age.label	<- c('Dt<=20','Dt<=45','Dt>45')
		for(i in seq_along(cd4.cut)[-1])
			set(YX.m4, YX.m4[, which(stage=='Diag' & t.Age<=age.cut[i] & t.Age>age.cut[i-1])], 'stage', age.label[i-1])
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
		YX.m4[, table(stage)]
		#	ART.suA.N ART.suA.Y  ART.vlNA    Dt<=20    Dt<=45     Dt>45         U       	UAm       UAy 
     	#	3172      4526        70       	136      	6351      1785      	7606      	1275      1291 
		YX.m4.fit4 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
		YX.m4.fit4.or	<- data.table( 	Dtl20= my.or.from.logit(YX.m4.fit4, 'stageDt<=20', 'stageDt<=45', subset(YX.m4, stage=='Dt<=20')[, sum(w)], subset(YX.m4, stage=='Dt<=45')[, sum(w)], 1.962),
										Dtg45= my.or.from.logit(YX.m4.fit4, 'stageDt>45', 'stageDt<=45', subset(YX.m4, stage=='Dt>45')[, sum(w)], subset(YX.m4, stage=='Dt<=45')[, sum(w)], 1.962)	)		
		#       Dtl20     	Dtg45
		#1: 	2.669321 	1.0475037
		#2: 	2.149098 	0.8275752
		#3: 	3.315471 	1.3258784
		#	young age at t significant, but generally older is not significant; Dtl25 is just signif 1.33 [1.03, 1.72]; Dtl23 is signif 1.85 [1.45, 1.37]
		#	there seems to be a trend with decreasing age <= 25	
	}
	#		try t2care	
	if(1)
	{
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAm' )				
		set(YX.m4, YX.m4[, which(stage=='Diag' & t2.care.t1=='>=0.25y')], 'stage', 'Dc>=0.25y')
		set(YX.m4, YX.m4[, which(stage=='Diag' & t2.care.t1=='>=0.5y')], 'stage', 'Dc>=0.5y')
		set(YX.m4, YX.m4[, which(stage=='Diag' & t2.care.t1=='other')], 'stage', 'Dc<0.25y')
		
		YX.m4s	<- subset(YX.m4, stage!='Diag')
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
		YX.m4s[, table(stage)]
		#	ART.suA.N ART.suA.Y  ART.vlNA       DAm       DAy  		Dc<0.25y 	Dc>=0.25y  	Dc>=0.5y         U       UAm       UAy 
     	#	3172      4526        70       		702       851      	5514       	340       	858      		7606      1275      1291 
		YX.m4s.fit5 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4s)
		YX.m4s.fit5.or	<- data.table( 	Dcg25= my.or.from.logit(YX.m4s.fit5, 'stageDc>=0.25y', 'stageDc<0.25y', subset(YX.m4s, stage=='Dc>=0.25y')[, sum(w)], subset(YX.m4s, stage=='Dc<0.25y')[, sum(w)], 1.962),
										Dcg50= my.or.from.logit(YX.m4s.fit5, 'stageDc>=0.5y', 'stageDc<0.25y', subset(YX.m4s, stage=='Dc>=0.5y')[, sum(w)], subset(YX.m4s, stage=='Dc<0.25y')[, sum(w)], 1.962)	)
		#      	Dcg25     Dcg50
		#1: 	1.361898 0.5660182
		#2: 	1.062440 0.4443110
		#3: 	1.745762 0.7210638	
		#	this is unexpected ..
		#      Dcg25     Dcg50
		#1: 1.852120 0.6782845
		#2: 1.399815 0.5179190
		#3: 2.450573 0.8883048
		#	remains after adjusting for Acute
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		set(YX.m4, YX.m4[, which(stage=='Diag' & t2.care.t1=='>=0.25y')], 'stage', 'Dc>=0.25y')
		set(YX.m4, YX.m4[, which(stage=='Diag' & t2.care.t1=='>=0.5y')], 'stage', 'Dc>=0.25y')
		set(YX.m4, YX.m4[, which(stage=='Diag' & t2.care.t1=='other')], 'stage', 'Dc<0.25y')
		
		YX.m4s	<- subset(YX.m4, stage!='Diag')
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
		YX.m4s[, table(stage)]
		#	ART.suA.N ART.suA.Y  ART.vlNA  Dc<0.25y  Dc>=0.25y          U       	UAm       UAy 
     	#	3172      4526        70       6952      1313      			7606      	1275      1291
		YX.m4s.fit5 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4s)
		YX.m4s.fit5.or	<- data.table( 	Dcg25= my.or.from.logit(YX.m4s.fit5, 'stageDc>=0.25y', 'stageDc<0.25y', subset(YX.m4s, stage=='Dc>=0.25y')[, sum(w)], subset(YX.m4s, stage=='Dc<0.25y')[, sum(w)], 1.962)	)
		#       Dcg25
		#1: 0.7277992
		#2: 0.5748090
		#3: 0.9215088
		#	this is unexpected ...
	}
	#		try t2supp	
	if(1)
	{
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		set(YX.m4, YX.m4[, which(t2.vl.supp<1)], 'stage', 's<=1' )
		set(YX.m4, YX.m4[, which(t2.vl.supp>1)], 'stage', 's>1' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag')], 'stage', 'UAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag')], 'stage', 'UAm' )
		set(YX.m4, YX.m4[, which(t.Patient%in%subset(YX.m4, ART.pulse=='Yes')[, unique(t.Patient)])], 'stage', 'Ppulse' )		
		
		ggplot(YX.m4, aes(stage, score.Y)) + geom_boxplot() + geom_jitter(size=0.4)	
	}
	#		try t2supp	
	if(1)
	{
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		
		YX.m4[, CD4.c:=NULL]
		cd4.label	<- c('D1<=350','D1<=550','D1>550')
		tmp			<- subset(clumsm.info, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
		setkey(tmp, Patient)
		tmp			<- unique(tmp)	
		tmp[, CD4.c:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
		cat(paste('\ntransmitters with no CD4 after diagnosis\n'))
		print( subset( tmp, is.na(CD4.c) ) )
		setnames(tmp, 'Patient', 't.Patient')
		set(tmp, tmp[, which(is.na(CD4.c))], 'CD4.c', 'D1.NA')
		YX.m4		<- merge(YX.m4, tmp, by='t.Patient')
		#	adjust for Acute and pulse
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAm' )
		set(YX.m4, YX.m4[, which(t.Patient%in%subset(YX.m4, ART.pulse=='Yes')[, unique(t.Patient)])], 'stage', 'Ppulse' )		
				
		
		fw.cut		<- c(-1, 1.5, 2.5, 30)
		fw.label	<- c('Ds<=1','Ds<=2','Dother')
		for(i in seq_along(fw.cut)[-1])
			set(YX.m4, YX.m4[, which(stage=='Diag' & CD4.c=='D1>550' & t2.vl.supp<=fw.cut[i] & t2.vl.supp>fw.cut[i-1])], 'stage', fw.label[i-1])
		set(YX.m4, YX.m4[, which(stage=='Diag' & CD4.c=='D1>550')], 'stage', 'Dother')		
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
		#
		YX.m4s			<- subset(YX.m4, !stage%in%c('ART.suA.N','ART.suA.Y','ART.vlNA','U','UAy','UAm'))
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
		YX.m4s[, table(stage)]
		ggplot(YX.m4s, aes(stage, score.Y)) + geom_boxplot() + geom_jitter(size=0.4)
		
		YX.m4.fit14 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
		YX.m4.fit14.or	<- data.table( 	Dsl1= my.or.from.logit(YX.m4.fit14, 'stageDs<=1', 'stageDother', subset(YX.m4, stage=='Ds<=1')[, sum(w)], subset(YX.m4, stage=='Dother')[, sum(w)], 1.962),
										Dsl2= my.or.from.logit(YX.m4.fit14, 'stageDs<=2', 'stageDother', subset(YX.m4, stage=='Ds<=2')[, sum(w)], subset(YX.m4, stage=='Dother')[, sum(w)], 1.962)	)		
		#	adjusting for Acute. adjusting for pulse: no such individual
		#	not adjusted for high CD4 at diagnosis
		#      Dsl1    Dsl2
		#1: 1.503301 1.449234
		#2: 1.154431 1.112343
		#3: 1.957600 1.888157
		#	adjusted for high CD4 at diagnosis
		#Dsl1      Dsl2
		#1: 2.177977 1.2301539
		#2: 1.669790 0.9291165
		#3: 2.840826 1.6287285
			
		cd4.label	<- c('D1<=350','D1<=550','D1>550')
		tmp			<- subset(clumsm.info, select=c(Patient, AnyPos_T1, PosCD4_T1, CD4_T1))
		setkey(tmp, Patient)
		tmp			<- unique(tmp)	
		tmp[, CD4.c:=cut(tmp[,CD4_T1], breaks=cd4.cut, labels=cd4.label, right=1)]
		cat(paste('\ntransmitters with no CD4 after diagnosis\n'))
		print( subset( tmp, is.na(CD4.c) ) )
		setnames(tmp, 'Patient', 't.Patient')
		set(tmp, tmp[, which(is.na(CD4.c))], 'CD4.c', 'D1.NA')
		YX.m4		<- merge(YX.m4, tmp, by='t.Patient')
		#	adjust for Acute and pulse
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])		
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAm' )
		set(YX.m4, YX.m4[, which(t.Patient%in%subset(YX.m4, ART.pulse=='Yes')[, unique(t.Patient)])], 'stage', 'Ppulse' )		
		#	adjust for young age
		age.cut		<- c(-1, 20, 22, 100)
		age.label	<- c('DC.T<=20','DC.T<=22', 'Diag')
		for(i in seq_along(age.cut)[-1])
			set(YX.m4, YX.m4[, which(stage=='Diag' & t.Age<=age.cut[i] & t.Age>age.cut[i-1])], 'stage', age.label[i-1])
		#	t2supp
		
		fw.cut		<- c(-1, 1, 2, 30)
		fw.label	<- c('Ds<=1','Ds<=2','Dother')
		for(i in seq_along(fw.cut)[-1])
			set(YX.m4, YX.m4[, which(stage=='Diag' & CD4>350 & t2.vl.supp<=fw.cut[i] & t2.vl.supp>fw.cut[i-1])], 'stage', fw.label[i-1])
		set(YX.m4, YX.m4[, which(stage=='Diag')], 'stage', 'Dother')		
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
		YX.m4[, table(stage)]
		#	ignoring CD4(t)
		#ART.suA.N ART.suA.Y  ART.vlNA       DAm       DAy  DC.T<=20  DC.T<=22    Dother     Ds<=1     Ds<=2    Ppulse         U       UAm       UAy 
		#3024      4497        70       	 702       844        72       148      5693       203       600       222      7593      1275
		#
		#ART.suA.N ART.suA.Y  ART.vlNA       DAm       DAy  DC.T<=20  DC.T<=22    Dother     Ds<=1     Ds<=2    Ppulse         U       UAm       UAy 
	    #3024      4497        70      		 702       844        72       148      6084        75       337       222      7593      1275      1269 

		YX.m4.fit14 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
		YX.m4.fit14.or		<- data.table( 	Dsl1= my.or.from.logit(YX.m4.fit14, 'stageDs<=1', 'stageDother', subset(YX.m4, stage=='Ds<=1')[, sum(w)], subset(YX.m4, stage=='Dother')[, sum(w)], 1.962),
											Dsl2= my.or.from.logit(YX.m4.fit14, 'stageDs<=2', 'stageDother', subset(YX.m4, stage=='Ds<=2')[, sum(w)], subset(YX.m4, stage=='Dother')[, sum(w)], 1.962)	)		
		#       Dsl1     Dsl2		ignoring CD4(t)
		#1: 1.505878 1.463572
		#2: 1.154378 1.121267
		#3: 1.964406 1.910376

		#       Dsl1     Dsl2		only CD4(t)>350
		#1: 2.433970 1.516979
		#2: 1.916189 1.173828
		#3: 3.091664 1.960445

	}
	#	how about contact==No ?
	if(1)
	{
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		set(YX.m4, YX.m4[, which(stage=='Diag' & contact=='No')], 'stage', 'D.cN')
		set(YX.m4, YX.m4[, which(stage=='Diag' & contact=='Yes')], 'stage', 'D.cY')
		YX.m4[, table(stage)]
		YX.m4s			<- subset(YX.m4, !stage%in%c('ART.suA.N','ART.suA.Y','ART.vlNA','U','UAy','UAm'))
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
		YX.m4s[, table(stage)]
		
		YX.m4s.fit13 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4s)
		summary(YX.m4s.fit13)
		ggplot(YX.m4s, aes(stage, score.Y)) + geom_boxplot() + geom_jitter(size=0.4)
	}
	#		try fw.up.med	
	if(1)
	{
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		fw.cut		<- c(-1, 2/12, 0.5, 2, 30)
		fw.label	<- c('Df<=0.16','Df<=0.5','Df<=2','Df>2')
		for(i in seq_along(fw.cut)[-1])
			set(YX.m4, YX.m4[, which(stage=='Diag' & fw.up.med<=fw.cut[i] & fw.up.med>fw.cut[i-1])], 'stage', fw.label[i-1])
		YX.m4[, table(stage)]
		
		YX.m4s	<- subset(YX.m4, stage!='Diag')
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
		YX.m4s[, table(stage)]
		#	ART.suA.N ART.suA.Y  ART.vlNA   Df<=0.2     Df<=2     	Df<=5         U       	UAm       UAy 
     	#	3172      4526        70      	2972      	5114       	66      		7606      	1275      1291
		YX.m4s.fit6 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4s)
		YX.m4s.fit6.or	<- data.table( 	Dfl02= my.or.from.logit(YX.m4s.fit6, 'stageDf<=0.16', 'stageDf<=2', subset(YX.m4s, stage=='Df<=0.16')[, sum(w)], subset(YX.m4s, stage=='Df<=2')[, sum(w)], 1.962),
										Dfl05= my.or.from.logit(YX.m4s.fit6, 'stageDf<=0.5', 'stageDf<=2', subset(YX.m4s, stage=='Df<=0.5')[, sum(w)], subset(YX.m4s, stage=='Df<=2')[, sum(w)], 1.962),
										Dfg2= my.or.from.logit(YX.m4s.fit6, 'stageDf>2', 'stageDf<=2', subset(YX.m4s, stage=='Df>2')[, sum(w)], subset(YX.m4s, stage=='Df<=2')[, sum(w)], 1.962)		)
		#       Dfl02    Dfl05      Dfg2
		#1: 2.633420 1.140633 0.6190514
		#2: 1.821406 0.856658 0.3310043
		#3: 3.807445 1.518743 1.1577631
		#	only frequent follow up associated with higher transmissibility. this is similar story to t2.care
		#	is this still true when we adjust for CDCC event ? --> YES
		YX.m4s	<- subset(YX.m4s, !stage%in%c('ART.suA.N','ART.suA.Y','ART.vlNA','U','UAm','UAy'))
		set(YX.m4s, NULL, 'CDCC', YX.m4s[, factor(CDCC)])
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
		YX.m4s.fit6b 	<- betareg(score.Y ~ stage:CDCC-1, link='logit', weights=w, data = YX.m4s)
		#
		#	does this still hold with acute and viral load ?
		#
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		#	set acute
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )		
		#	define lRNA.gm
		YX.m4[, lRNA.gm:=NULL]
		YX.m4	<- merge(YX.m4, YX.m4[, {
							tmp<- which(!is.na(lRNA))					
							if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1])
								lRNA.gm		<- lRNA_TL[1]						
							else if(!length(tmp))
								lRNA.gm		<- NA_real_
							else
								lRNA.gm		<- mean(lRNA[tmp])
							list( lRNA.gm= lRNA.gm )
						}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
		#	
		fw.label	<- c('Df<=0.16','Df<=0.5','Df<=2','Df>2')
		for(i in seq_along(fw.cut)[-1])
			set(YX.m4, YX.m4[, which(stage=='Diag' & fw.up.med<=fw.cut[i] & fw.up.med>fw.cut[i-1])], 'stage', fw.label[i-1])
		YX.m4[, table(stage)]
		#	focus on diagnosis stage
		YX.m4s	<- subset(YX.m4, !stage%in%c('ART.suA.N','ART.suA.Y','ART.vlNA','U','UAy','UAm'))
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
		YX.m4s[, table(stage)]			
		#
		YX.m4s.fit6.fw 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4s)
		YX.m4s.fit6.vl.fw 	<- betareg(score.Y ~ lRNA.gm+stage-1, link='logit', weights=w, data = YX.m4s)
		#	there are some transmitters with long follow up that obtain a high score, but most tend to be low
		ggplot(YX.m4s, aes(stage, score.Y)) + geom_boxplot() + geom_jitter(size=0.4)
		#		
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])		
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )
		YX.m4s	<- subset(YX.m4, !stage%in%c('ART.suA.N','ART.suA.Y','ART.vlNA','U','UAy','UAm'))
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])		
		YX.m4s.fit6.vl 		<- betareg(score.Y ~ lRNA.gm+stage-1, link='logit', weights=w, data = YX.m4s)
		AIC( YX.m4s.fit6.fw, YX.m4s.fit6.vl, YX.m4s.fit6.vl.fw )
		#                  df       AIC
		#YX.m4s.fit6.fw     6 -414.8447		acute + fw
		#YX.m4s.fit6.vl     5 -391.6777		acute + vl
		#YX.m4s.fit6.vl.fw  7 -389.2229		acute + fw + vl
	}
	#		see if single VL thresh stratifies diagnosed	- use max VL per infection window
	if(0)
	{		
		YX.m4[, lRNA.c:=NULL]
		YX.m4[, lRNA_TL:=NULL]
		YX.m4[, PoslRNA_TL:=NULL]		
		tmp		<- unique(subset( clumsm.info, select=c(Patient, PoslRNA_T1) ))
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		YX.m4	<- merge(YX.m4, tmp, by= 't.Patient')
		YX.m4	<- merge(YX.m4, YX.m4[, {	
											tmp<- which(!is.na(lRNA))
											list(lRNA_TL= ifelse(length(tmp)>0, lRNA[tail(tmp,1)], NA_real_), 	PoslRNA_TL=ifelse(length(tmp)>0, t[tail(tmp,1)], NA_real_))
										}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
		
		VL.win		<- c( seq(1e4, 1e5, by=1e4), seq(2e5, 5e5, by=1e5) ) 
		YX.m4.VL	<- sapply(VL.win, function(VL.cur)
				{
					cat(paste('\nprocess VL=', VL.cur))									
					VL.cur	<- log10(VL.cur)
					set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
					tmp		<- YX.m4[, {
										tmp<- which(!is.na(lRNA))
										#if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE)) stop('found')
										if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1]<=VL.cur)
											lRNA.c	<- 'VLw.L'						
										else if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1]>VL.cur)
											lRNA.c	<- 'VLw.H'						
										else if(!length(tmp))
											lRNA.c	<- 'VLw.NA'
										else if(max(lRNA[tmp])>VL.cur)
											lRNA.c	<- 'VLw.H'
										else
											lRNA.c	<- 'VLw.L'
										list( lRNA.c= lRNA.c )
							}, by=c('Patient','t.Patient')]
					YX.m4	<- merge(YX.m4, tmp, by=c('Patient','t.Patient'))
					#	actually, the number of entries with missing VL does not improve based on the contact stuff because we there is no such individual
					set(YX.m4, YX.m4[,which(stage=='Diag' & lRNA.c=='VLw.L')],'stage', 'DVLw.L' )
					set(YX.m4, YX.m4[,which(stage=='Diag' & lRNA.c=='VLw.H')],'stage', 'DVLw.H' )
					set(YX.m4, YX.m4[,which(stage=='Diag')],'stage', 'DVLw.NA' )
					set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
					YX.m4[, lRNA.c:=NULL]
					YX.m4.fit7 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
					
					
					ans			<- c(		YX.m4[, length(which(stage=='DVLw.H'))], YX.m4[, length(which(stage=='DVLw.L'))],
											subset(YX.m4, stage=='DVLw.H')[, sum(w)], subset(YX.m4, stage=='DVLw.L')[, sum(w)],
											my.or.from.logit(YX.m4.fit7, 'stageDVLw.H', 'stageDVLw.L', subset(YX.m4, stage=='DVLw.H')[, sum(w)], subset(YX.m4, stage=='DVLw.L')[, sum(w)], 1.962),						
											my.or.from.logit(YX.m4.fit7, 'stageDVLw.NA', 'stageDVLw.L', subset(YX.m4, stage=='DVLw.NA')[, sum(w)], subset(YX.m4, stage=='DVLw.L')[, sum(w)], 1.962),
											logLik(YX.m4.fit7), YX.m4.fit7$pseudo.r.squared, 10^VL.cur	)
					names(ans)	<- c('n.DVLw.H','n.DVLw.L','w.DVLw.H','w.DVLw.L','or.HL','or.HL.l95','or.HL.u95','or.NAL','or.NAL.l95','or.NAL.u95','lkl','r2','VL.thresh')
					ans
				})
		YX.m4.VL	<- as.data.table(t(YX.m4.VL))
		plot(YX.m4.VL[, VL.thresh], YX.m4.VL[, or.HL], type='p', pch=18, ylim= c( min(YX.m4.VL[, or.HL.l95], na.rm=TRUE),max(YX.m4.VL[, or.HL.u95], na.rm=TRUE))  )
		dummy	<- YX.m4.VL[,	lines( rep(VL.thresh,2), c(or.HL.l95,or.HL.u95))	, by='VL.thresh']
		abline(h=1, lty=2)
		#	suggests to consider Acute separately!
		YX.m4[, PoslRNA_TL:=NULL]
		YX.m4[, lRNA_TL:=NULL]
		YX.m4[, t.PoslRNA_T1:=NULL]
	}	
	#	investigate Acute=='Yes' and Acute=='Maybe' at diagnosis and infection window periods close to diagnosis
	if(0)
	{
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 1)], 'stage', 'DAm' )
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
		YX.m4[, table(stage)]
		#	ART.suA.N ART.suA.Y  ART.vlNA       DAm       DAy      Diag         U       UAm       UAy 
     	#	3172      4526        70       		702       851      6719      7606      1275      1291 
		YX.m4.fit8 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
		YX.m4.fit8.or	<- data.table( 	DAy= my.or.from.logit(YX.m4.fit8, 'stageDAy', 'stageDiag', subset(YX.m4, stage=='DAy')[, sum(w)], subset(YX.m4, stage=='Diag')[, sum(w)], 1.962),
										DAm= my.or.from.logit(YX.m4.fit8, 'stageDAm', 'stageDiag', subset(YX.m4, stage=='DAy')[, sum(w)], subset(YX.m4, stage=='Diag')[, sum(w)], 1.962)	)
		#        DAy      DAm
		#1: 2.453331 3.061342
		#2: 1.930366 2.390187
		#3: 3.117974 3.920955								
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.5)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.5)], 'stage', 'DAm' )
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
		YX.m4[, table(stage)]
		#	ART.suA.N ART.suA.Y  ART.vlNA       DAm       DAy      Diag         U       UAm       UAy 
     	#	3172      4526        70       		436       523      7313      7606      1275      1291  
		YX.m4.fit9 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
		YX.m4.fit9.or	<- data.table( 	DAy= my.or.from.logit(YX.m4.fit9, 'stageDAy', 'stageDiag', subset(YX.m4, stage=='DAy')[, sum(w)], subset(YX.m4, stage=='Diag')[, sum(w)], 1.962),
										DAm= my.or.from.logit(YX.m4.fit9, 'stageDAm', 'stageDiag', subset(YX.m4, stage=='DAy')[, sum(w)], subset(YX.m4, stage=='Diag')[, sum(w)], 1.962)	)
		#               DAy      DAm
		#1: 2.668173 3.266706
		#2: 2.115017 2.585143
		#3: 3.365999 4.127960
		#
		#	YES -- they are strongly associated with transmission at a fairly small time period after diagnosis	
	}
	#	interplay Acute and VL high -- max VL
	if(1)
	{		
		YX.m4[, lRNA.c:=NULL]		
		VL.win		<- c( seq(1e4, 1e5, by=1e4), seq(1.5e5, 5e5, by=5e4) ) 
		YX.m4.VL	<- sapply(VL.win, function(VL.cur)
				{
					cat(paste('\nprocess VL=', VL.cur))									
					VL.cur	<- log10(VL.cur)
					set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
					tmp		<- YX.m4[, {
								tmp<- which(!is.na(lRNA))
								#if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE)) stop('found')
								if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1]<=VL.cur)
									lRNA.c	<- 'VLw.L'						
								else if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1]>VL.cur)
									lRNA.c	<- 'VLw.H'						
								else if(!length(tmp))
									lRNA.c	<- 'VLw.NA'
								else if(max(lRNA[tmp])>VL.cur)
									lRNA.c	<- 'VLw.H'
								else
									lRNA.c	<- 'VLw.L'
								list( lRNA.c= lRNA.c )
							}, by=c('Patient','t.Patient')]
					YX.m4	<- merge(YX.m4, tmp, by=c('Patient','t.Patient'))
					#	actually, the number of entries with missing VL does not improve based on the contact stuff because we there is no such individual
					set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
					set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )					
					set(YX.m4, YX.m4[,which(stage=='Diag' & lRNA.c=='VLw.L')],'stage', 'DVLw.L' )
					set(YX.m4, YX.m4[,which(stage=='Diag' & lRNA.c=='VLw.H')],'stage', 'DVLw.H' )
					set(YX.m4, YX.m4[,which(stage=='Diag')],'stage', 'DVLw.NA' )
					set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
					YX.m4[, lRNA.c:=NULL]
					#YX.m4[, table(stage)]
					YX.m4.fit10 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
					ans			<- c(		YX.m4[, length(which(stage=='DVLw.H'))], YX.m4[, length(which(stage=='DVLw.L'))],
											subset(YX.m4, stage=='DVLw.H')[, sum(w)], subset(YX.m4, stage=='DVLw.L')[, sum(w)],
											my.or.from.logit(YX.m4.fit10, 'stageDVLw.H', 'stageDVLw.L', subset(YX.m4, stage=='DVLw.H')[, sum(w)], subset(YX.m4, stage=='DVLw.L')[, sum(w)], 1.962),						
											my.or.from.logit(YX.m4.fit10, 'stageDVLw.NA', 'stageDVLw.L', subset(YX.m4, stage=='DVLw.NA')[, sum(w)], subset(YX.m4, stage=='DVLw.L')[, sum(w)], 1.962),
											logLik(YX.m4.fit10), YX.m4.fit10$pseudo.r.squared, 10^VL.cur	)
					names(ans)	<- c('n.DVLw.H','n.DVLw.L','w.DVLw.H','w.DVLw.L','or.HL','or.HL.l95','or.HL.u95','or.NAL','or.NAL.l95','or.NAL.u95','lkl','r2','VL.thresh')
					ans
				})
		YX.m4.VL	<- as.data.table(t(YX.m4.VL))
		plot(YX.m4.VL[, VL.thresh], YX.m4.VL[, or.HL], type='p', pch=18, ylim= c( min(YX.m4.VL[, or.HL.l95], na.rm=TRUE),max(YX.m4.VL[, or.HL.u95], na.rm=TRUE))  )
		dummy	<- YX.m4.VL[,	lines( rep(VL.thresh,2), c(or.HL.l95,or.HL.u95))	, by='VL.thresh']
		abline(h=1, lty=2)
		#	Once Acute is considered separately, there is no indication that VL.high is associated with higher transmissibility
		#	typical: 0.9985499 [ 0.7548396  1.320946] for VL 1e5
		YX.m4[, PoslRNA_TL:=NULL]
		YX.m4[, lRNA_TL:=NULL]
		YX.m4[, t.PoslRNA_T1:=NULL]
	}
	#	compare model using lRNA as linear predictor over three possibilities: VL.t VL.mxw VL.gm
	if(1)
	{
		#define lRNA.mxw and lRNA.gm
		YX.m4	<- merge(YX.m4, YX.m4[, {
											tmp<- which(!is.na(lRNA))					
											if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1])
												lRNA.gm		<- lRNA.mxw	<- lRNA_TL[1]						
											else if(!length(tmp))
												lRNA.gm		<- lRNA.mxw	<- NA_real_
											else
											{
												lRNA.gm		<- mean(lRNA[tmp])
												lRNA.mxw	<- max(lRNA[tmp])
											}
											list( lRNA.gm= lRNA.gm, lRNA.mxw= lRNA.mxw )
										}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))	
		#	actually, the number of entries with missing VL does not improve based on the contact stuff because we there is no such individual
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )							
		#	focus on diagnosis stage
		YX.m4s	<- subset(YX.m4, !stage%in%c('ART.suA.N','ART.suA.Y','ART.vlNA','U','UAy','UAm'))
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
		YX.m4s[, table(stage)]				
		#
		YX.m4.fit11.t 	<- betareg(score.Y ~ lRNA+stage-1, link='logit', weights=w, data = YX.m4s)
		YX.m4.fit11.gm 	<- betareg(score.Y ~ lRNA.gm+stage-1, link='logit', weights=w, data = YX.m4s)
		YX.m4.fit11.mxw <- betareg(score.Y ~ lRNA.mxw+stage-1, link='logit', weights=w, data = YX.m4s)
		
		sapply( list(YX.m4.fit11.t, YX.m4.fit11.mxw, YX.m4.fit11.gm), '[[', 'pseudo.r.squared' ) 
		#	0.05427173 0.05646519 0.06212560
		AIC(YX.m4.fit11.t, YX.m4.fit11.mxw, YX.m4.fit11.gm)
		#                df       AIC
		#YX.m4.fit11.t    5 -339.9402
		#YX.m4.fit11.mxw  5 -389.6059
		#YX.m4.fit11.gm   5 -391.6777
		summary( YX.m4.fit11.gm )
		data.table( 	DAy= my.or.from.logit(YX.m4.fit11.gm, 'stageDAy', 'stageDiag', subset(YX.m4, stage=='DAy')[, sum(w)], subset(YX.m4, stage=='Diag')[, sum(w)], 1.962),
						DAm= my.or.from.logit(YX.m4.fit11.gm, 'stageDAm', 'stageDiag', subset(YX.m4, stage=='DAm')[, sum(w)], subset(YX.m4, stage=='Diag')[, sum(w)], 1.962)		)
		#         DAy       DAm
		#1: 2.2460174 2.7237358
		#2: 0.7100023 0.8615171
		#3: 7.1050396 8.6112475
		YX.m4.fit11.noa <- betareg(score.Y ~ lRNA.mxw-1, link='logit', weights=w, data = YX.m4s)		
		sapply( list(YX.m4.fit11.noa, YX.m4.fit11.gm), '[[', 'pseudo.r.squared' )
		#	0.01663442 0.06212560
		AIC(YX.m4.fit11.noa, YX.m4.fit11.gm)
		#					df       AIC
		#YX.m4.fit11.noa  2 -384.0294
		#YX.m4.fit11.gm   5 -391.6777
	}
	#	interplay Acute and VL high using Geometric mean of RNA, or the log Geometric mean ( mean of lRNA )
	if(1)
	{		
		YX.m4[, lRNA.c:=NULL]
		VL.win		<- c( seq(1e4, 1e5, by=1e4), seq(1.5e5, 5e5, by=5e4) ) 
		YX.m4.VL	<- sapply(VL.win, function(VL.cur)
				{
					cat(paste('\nprocess VL=', VL.cur))									
					VL.cur	<- log10(VL.cur)
					set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
					tmp		<- YX.m4[, {
								tmp<- which(!is.na(lRNA))
								#if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE)) stop('found')
								if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1]<=VL.cur)
									lRNA.c	<- 'VLw.L'						
								else if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1]>VL.cur)
									lRNA.c	<- 'VLw.H'						
								else if(!length(tmp))
									lRNA.c	<- 'VLw.NA'
								else if(mean(lRNA[tmp])>VL.cur)
									lRNA.c	<- 'VLw.H'
								else
									lRNA.c	<- 'VLw.L'
								list( lRNA.c= lRNA.c )
							}, by=c('Patient','t.Patient')]
					YX.m4	<- merge(YX.m4, tmp, by=c('Patient','t.Patient'))
					#	actually, the number of entries with missing VL does not improve based on the contact stuff because we there is no such individual
					set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
					set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )					
					set(YX.m4, YX.m4[,which(stage=='Diag' & lRNA.c=='VLw.L')],'stage', 'DVLw.L' )
					set(YX.m4, YX.m4[,which(stage=='Diag' & lRNA.c=='VLw.H')],'stage', 'DVLw.H' )
					set(YX.m4, YX.m4[,which(stage=='Diag')],'stage', 'DVLw.NA' )
					set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
					YX.m4[, lRNA.c:=NULL]
					#YX.m4[, table(stage)]
					YX.m4.fit10 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
					ans			<- c(		YX.m4[, length(which(stage=='DVLw.H'))], YX.m4[, length(which(stage=='DVLw.L'))],
							subset(YX.m4, stage=='DVLw.H')[, sum(w)], subset(YX.m4, stage=='DVLw.L')[, sum(w)],
							my.or.from.logit(YX.m4.fit10, 'stageDVLw.H', 'stageDVLw.L', subset(YX.m4, stage=='DVLw.H')[, sum(w)], subset(YX.m4, stage=='DVLw.L')[, sum(w)], 1.962),						
							my.or.from.logit(YX.m4.fit10, 'stageDVLw.NA', 'stageDVLw.L', subset(YX.m4, stage=='DVLw.NA')[, sum(w)], subset(YX.m4, stage=='DVLw.L')[, sum(w)], 1.962),
							logLik(YX.m4.fit10), YX.m4.fit10$pseudo.r.squared, 10^VL.cur	)
					names(ans)	<- c('n.DVLw.H','n.DVLw.L','w.DVLw.H','w.DVLw.L','or.HL','or.HL.l95','or.HL.u95','or.NAL','or.NAL.l95','or.NAL.u95','lkl','r2','VL.thresh')
					ans
				})
		YX.m4.VL	<- as.data.table(t(YX.m4.VL))
		plot(YX.m4.VL[, VL.thresh], YX.m4.VL[, or.HL], type='p', pch=18, ylim= c( min(YX.m4.VL[, or.HL.l95], na.rm=TRUE),max(YX.m4.VL[, or.HL.u95], na.rm=TRUE)) )
		dummy	<- YX.m4.VL[,	lines( rep(VL.thresh,2), c(or.HL.l95,or.HL.u95))	, by='VL.thresh']
		abline(h=1, lty=2)
		#	with dt=0.25, this is now significant
		#
		#	compare threshold to linear on lRNA model
		VL.cur	<- 5
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		tmp		<- YX.m4[, {
					tmp<- which(!is.na(lRNA))
					#if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE)) stop('found')
					if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1]<=VL.cur)
						lRNA.c	<- 'VLw.L'						
					else if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1]>VL.cur)
						lRNA.c	<- 'VLw.H'						
					else if(!length(tmp))
						lRNA.c	<- 'VLw.NA'
					else if(mean(lRNA[tmp])>VL.cur)
						lRNA.c	<- 'VLw.H'
					else
						lRNA.c	<- 'VLw.L'
					list( lRNA.c= lRNA.c )
				}, by=c('Patient','t.Patient')]
		YX.m4	<- merge(YX.m4, tmp, by=c('Patient','t.Patient'))
		#	actually, the number of entries with missing VL does not improve based on the contact stuff because we there is no such individual
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )					
		set(YX.m4, YX.m4[,which(stage=='Diag' & lRNA.c=='VLw.L')],'stage', 'DVLw.L' )
		set(YX.m4, YX.m4[,which(stage=='Diag' & lRNA.c=='VLw.H')],'stage', 'DVLw.H' )
		set(YX.m4, YX.m4[,which(stage=='Diag')],'stage', 'DVLw.NA' )
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
		YX.m4[, lRNA.c:=NULL]
		YX.m4s	<- subset(YX.m4, !stage%in%c('ART.suA.N','ART.suA.Y','ART.vlNA','U','UAy','UAm','DAy','DAm'))
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
		YX.m4s[, table(stage)]				
		YX.m4.fit10a 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4s)
		summary(YX.m4.fit10a)
		#Log-likelihood: 217.4 on 6 Df
		#Pseudo R-squared: 0.05876		
		#Log-likelihood: 152.3 on 4 Df		(excluding acute Diagnosed)
		#Pseudo R-squared: 0.008702
		
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )					
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])		
		tmp		<- YX.m4[, {
					tmp<- which(!is.na(lRNA))
					#if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE)) stop('found')
					if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) )
						lRNA.c	<- lRNA_TL						
					else if(!length(tmp))
						lRNA.c	<- NA_real_
					else
						lRNA.c	<- mean(lRNA[tmp])
					list( lRNA.c= lRNA.c, lRNA.g=exp(lRNA.c) )
				}, by=c('Patient','t.Patient')]
		YX.m4	<- merge(YX.m4, tmp, by=c('Patient','t.Patient'))
		YX.m4s	<- subset(YX.m4, !stage%in%c('ART.suA.N','ART.suA.Y','ART.vlNA','U','UAy','UAm','DAy','DAm'))
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
		YX.m4s[, table(stage)]
		YX.m4.fit10b 	<- betareg(score.Y ~ lRNA.c, link='logit', weights=w, data = YX.m4s)
		summary(YX.m4.fit10b)
		#Log-likelihood: 197.3 on 5 Df
		#Pseudo R-squared: 0.06579
		#Log-likelihood: 141.3 on 3 Df		(excluding acute Diagnosed)
		#Pseudo R-squared: 0.0236		
		YX.m4.fit10c 	<- betareg(score.Y ~ lRNA.g+stage-1, link='logit', weights=w, data = YX.m4s)
		summary(YX.m4.fit10c)
		#Log-likelihood: 197.5 on 5 Df
		#Pseudo R-squared: 0.06377
		#Log-likelihood: 141.5 on 4 Df		(excluding acute Diagnosed)
		#Pseudo R-squared: 0.02102

		tmp<- subset( YX.m4s, stage=='Diag', select=c(lRNA.c, score.Y))
		setkey(tmp, lRNA.c)
		tmp	<- subset(tmp, !is.na(lRNA.c))
		plot(tmp[, lRNA.c], tmp[, log( score.Y/(1-score.Y) )], pch=18)
		lines( lowess(tmp[, lRNA.c], tmp[, log( score.Y/(1-score.Y) )]), col='blue')	#	linear fit prob not too bad
		
		YX.m4[, lRNA.c:=NULL]
		
		YX.m4[, PoslRNA_TL:=NULL]
		YX.m4[, lRNA_TL:=NULL]
		YX.m4[, t.PoslRNA_T1:=NULL]
	}
	#	interplay Acute and VL high using every time period individually
	if(1)
	{		
		YX.m4[, lRNA.c:=NULL]
		VL.win		<- c( seq(1e4, 1e5, by=1e4), seq(1.5e5, 5e5, by=5e4) ) 
		YX.m4.VL	<- sapply(VL.win, function(VL.cur)
				{
					cat(paste('\nprocess VL=', VL.cur))									
					VL.cur	<- log10(VL.cur)
					set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
					set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
					#set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )										
					set(YX.m4, YX.m4[, which(stage=='Diag' & is.na(lRNA) & !is.na(lRNA_TL) & t>PoslRNA_TL & (t-PoslRNA_TL)<0.5 & lRNA_TL[1]<=VL.cur)],'stage', 'DVLw.L' )
					set(YX.m4, YX.m4[, which(stage=='Diag' & !is.na(lRNA) & lRNA<=VL.cur)],'stage', 'DVLw.L' )
					set(YX.m4, YX.m4[, which(stage=='Diag' & is.na(lRNA) & !is.na(lRNA_TL) & t>PoslRNA_TL & (t-PoslRNA_TL)<0.5 & lRNA_TL[1]>VL.cur)],'stage', 'DVLw.H' )					
					set(YX.m4, YX.m4[, which(stage=='Diag' & !is.na(lRNA) & lRNA>VL.cur)],'stage', 'DVLw.H' )
					set(YX.m4, YX.m4[, which(stage=='Diag')],'stage', 'DVLw.NA' )
					set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])					
					#YX.m4[, table(stage)]
					YX.m4.fit10 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
					ans			<- c(		YX.m4[, length(which(stage=='DVLw.H'))], YX.m4[, length(which(stage=='DVLw.L'))],
							subset(YX.m4, stage=='DVLw.H')[, sum(w)], subset(YX.m4, stage=='DVLw.L')[, sum(w)],
							my.or.from.logit(YX.m4.fit10, 'stageDVLw.H', 'stageDVLw.L', subset(YX.m4, stage=='DVLw.H')[, sum(w)], subset(YX.m4, stage=='DVLw.L')[, sum(w)], 1.962),						
							my.or.from.logit(YX.m4.fit10, 'stageDVLw.NA', 'stageDVLw.L', subset(YX.m4, stage=='DVLw.NA')[, sum(w)], subset(YX.m4, stage=='DVLw.L')[, sum(w)], 1.962),
							logLik(YX.m4.fit10), YX.m4.fit10$pseudo.r.squared, 10^VL.cur	)
					names(ans)	<- c('n.DVLw.H','n.DVLw.L','w.DVLw.H','w.DVLw.L','or.HL','or.HL.l95','or.HL.u95','or.NAL','or.NAL.l95','or.NAL.u95','lkl','r2','VL.thresh')
					ans
				})
		YX.m4.VL	<- as.data.table(t(YX.m4.VL))
		plot(YX.m4.VL[, VL.thresh], YX.m4.VL[, or.HL], type='p', pch=18, ylim= c( min(YX.m4.VL[, or.HL.l95], na.rm=TRUE),max(YX.m4.VL[, or.HL.u95], na.rm=TRUE))  )
		dummy	<- YX.m4.VL[,	lines( rep(VL.thresh,2), c(or.HL.l95,or.HL.u95))	, by='VL.thresh']
		abline(h=1, lty=2)
		#	Once Acute is considered separately, there is no indication that VL.high is associated with higher transmissibility
		#	typical: 0.9985499 [ 0.7548396  1.320946] for VL 1e5
		YX.m4[, PoslRNA_TL:=NULL]
		YX.m4[, lRNA_TL:=NULL]
		YX.m4[, t.PoslRNA_T1:=NULL]
	}
	# combine Acute, lRNA.gm, young age
	if(1)
	{
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		#	add age at diagnosis
		YX.m4[, t.AnyPos_A:=NULL]
		tmp		<- subset( clumsm.info, select=c(Patient, DateBorn, AnyPos_T1) )
		setkey(tmp, Patient)
		tmp		<- unique(tmp)
		tmp[, AnyPos_A:= AnyPos_T1-DateBorn]
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		tmp		<- merge( unique( subset(YX.m4, select=t.Patient) ), tmp, by='t.Patient' )		
		YX.m4	<- merge(YX.m4, subset(tmp, select=c(t.Patient, t.AnyPos_A)), by='t.Patient')	
		#	define lRNA.gm
		YX.m4[, lRNA.gm:=NULL]
		YX.m4	<- merge(YX.m4, YX.m4[, {
							tmp<- which(!is.na(lRNA))					
							if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1])
								lRNA.gm		<- lRNA_TL[1]						
							else if(!length(tmp))
								lRNA.gm		<- NA_real_
							else
								lRNA.gm		<- mean(lRNA[tmp])
							list( lRNA.gm= lRNA.gm )
						}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))	
		#	set acute
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )
		#	set young age groups
		age.cut		<- c(-1, 20, 22, 24, 100)
		age.label	<- c('D.T<=20','D.T<=22', 'D.T<=24', 'D.T>25')
		age.cut		<- c(-1, 23, 100)
		age.label	<- c('D.T<=24', 'D.T>25')
		
		YX.m4[, t.Age.c:= cut(YX.m4[, t.Age], breaks=age.cut, labels=age.label)]
		#	focus on diagnosis stage
		YX.m4s	<- subset(YX.m4, !stage%in%c('ART.suA.N','ART.suA.Y','ART.vlNA','U','UAy','UAm'))
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
		YX.m4s[, table(stage)]			
		#
		YX.m4.fit12.age.stage 	<- betareg(score.Y ~ lRNA.gm+t.Age.c+stage-1, link='logit', weights=w, data = YX.m4s)
		YX.m4.fit12.age		 	<- betareg(score.Y ~ lRNA.gm+t.Age.c-1, link='logit', weights=w, data = YX.m4s)
		YX.m4.fit12.stage 	<- betareg(score.Y ~ lRNA.gm+stage-1, link='logit', weights=w, data = YX.m4s)
		AIC(YX.m4.fit12.age, YX.m4.fit12.stage, YX.m4.fit12.age.stage)
		#	using age<=23 gave largest r2 for YX.m4.fit12.age.stage -- but not lower AIC than stage+lRNA.gw	
		#                      df       AIC
		#YX.m4.fit12.age        4 -384.3063
		#YX.m4.fit12.stage      5 -391.6777
		#YX.m4.fit12.age.stage  6 -390.5392
		sapply( list(YX.m4.fit12.age, YX.m4.fit12.stage, YX.m4.fit12.age.stage), '[[', 'pseudo.r.squared' )
		#	0.04103416 0.06212560 0.07381119
	}
	# combine Acute, lRNA.gm, young age among chronic
	if(1)
	{
		set(YX.m4, NULL, 'stage', YX.m4[, stage.orig])
		#	add age at diagnosis
		YX.m4[, t.AnyPos_A:=NULL]
		tmp		<- subset( clumsm.info, select=c(Patient, DateBorn, AnyPos_T1) )
		setkey(tmp, Patient)
		tmp		<- unique(tmp)
		tmp[, AnyPos_A:= AnyPos_T1-DateBorn]
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		tmp		<- merge( unique( subset(YX.m4, select=t.Patient) ), tmp, by='t.Patient' )		
		YX.m4	<- merge(YX.m4, subset(tmp, select=c(t.Patient, t.AnyPos_A)), by='t.Patient')	
		#	define lRNA.gm
		YX.m4[, lRNA.gm:=NULL]
		YX.m4	<- merge(YX.m4, YX.m4[, {
							tmp<- which(!is.na(lRNA))					
							if(!length(tmp) & all( t>PoslRNA_TL, na.rm=TRUE ) & any( (t-PoslRNA_TL)<0.5, na.rm=TRUE) & !is.na(lRNA_TL[1]) & lRNA_TL[1])
								lRNA.gm		<- lRNA_TL[1]						
							else if(!length(tmp))
								lRNA.gm		<- NA_real_
							else
								lRNA.gm		<- mean(lRNA[tmp])
							list( lRNA.gm= lRNA.gm )
						}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))	
		#	set acute
		set(YX.m4, YX.m4[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAy' )
		set(YX.m4, YX.m4[, which(t.isAcute=='Maybe' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stage', 'DAm' )
		#	set young age groups
		age.cut		<- c(-1, 22, 100)
		age.label	<- c('DC.T<=22', 'Diag')
		for(i in seq_along(age.cut)[-1])
			set(YX.m4, YX.m4[, which(stage=='Diag' & t.Age<=age.cut[i] & t.Age>age.cut[i-1])], 'stage', age.label[i-1])		
		#	focus on diagnosis stage
		YX.m4s	<- subset(YX.m4, !stage%in%c('ART.suA.N','ART.suA.Y','ART.vlNA','U','UAy','UAm'))
		set(YX.m4s, NULL, 'stage', YX.m4s[, factor(as.character(stage))])
		YX.m4s[, table(stage)]			
		#
		YX.m4.fit13.age.stage 	<- betareg(score.Y ~ lRNA.gm+stage-1, link='logit', weights=w, data = YX.m4s)
		AIC(YX.m4.fit12.age, YX.m4.fit12.stage, YX.m4.fit12.age.stage, YX.m4.fit13.age.stage)
		#                      df       AIC		
		#YX.m4.fit12.age        4 -383.4353
		#YX.m4.fit12.stage      5 -391.6777
		#YX.m4.fit12.age.stage  6 -389.7269
		#YX.m4.fit13.age.stage  6 -390.8501		slightly higher age effect among chronic, but still not significant
	}
	# combine Acute, lRNA.gm, time2viralsuppression
	if(1)
	{
		


#	adjust for young age



		summary(YX.m4.fit12)
		
		sapply( list(YX.m4.fit11.t, YX.m4.fit11.mxw, YX.m4.fit11.gm), '[[', 'pseudo.r.squared' ) 
		#	0.05427173 0.05646519 0.06212560
		AIC(YX.m4.fit11.t, YX.m4.fit11.mxw, YX.m4.fit11.gm)
		#                df       AIC
		#YX.m4.fit11.t    5 -339.9402
		#YX.m4.fit11.mxw  5 -389.6059
		#YX.m4.fit11.gm   5 -391.6777
		summary( YX.m4.fit11.gm )
		data.table( 	DAy= my.or.from.logit(YX.m4.fit11.gm, 'stageDAy', 'stageDiag', subset(YX.m4, stage=='DAy')[, sum(w)], subset(YX.m4, stage=='Diag')[, sum(w)], 1.962),
				DAm= my.or.from.logit(YX.m4.fit11.gm, 'stageDAm', 'stageDiag', subset(YX.m4, stage=='DAm')[, sum(w)], subset(YX.m4, stage=='Diag')[, sum(w)], 1.962)		)
		#         DAy       DAm
		#1: 2.2460174 2.7237358
		#2: 0.7100023 0.8615171
		#3: 7.1050396 8.6112475
		YX.m4.fit11.noa <- betareg(score.Y ~ lRNA.mxw-1, link='logit', weights=w, data = YX.m4s)		
		sapply( list(YX.m4.fit11.noa, YX.m4.fit11.gm), '[[', 'pseudo.r.squared' )
		#	0.01663442 0.06212560
		AIC(YX.m4.fit11.noa, YX.m4.fit11.gm)
		
		
		
		#	set age 45 at diagnosis
		set(YX.m4, YX.m4[, which(t.AnyPos_A>45 & stage=='Diag' & t-t.AnyPos_T1 < 0.5)], 'stage', 'DC.T1>45' )
		#			
		set(YX.m4, NULL, 'stage', YX.m4[, factor(as.character(stage))])
		table( YX.m4[, stage])
		#	ART.suA.N ART.suA.Y  ART.vlNA       DAm       DAy  		DC.T<=20  	DC.T<=25   DC.T>25  	DC.T1>45         U       UAm       UAy 
     	#	3172      4526        70       		436       523       75       	427      	6533       278      		7606      1275      1291
		YX.m4.fit11 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m4)
		YX.m4.fit11.or	<- data.table( 	DAy= my.or.from.logit(YX.m4.fit11, 'stageDAy', 'stageDC.T>25', subset(YX.m4, stage=='DAy')[, sum(w)], subset(YX.m4, stage=='DC.T>25')[, sum(w)], 1.962),
										DAm= my.or.from.logit(YX.m4.fit11, 'stageDAm', 'stageDC.T>25', subset(YX.m4, stage=='DAy')[, sum(w)], subset(YX.m4, stage=='DC.T>25')[, sum(w)], 1.962),
										D1g45= my.or.from.logit(YX.m4.fit11, 'stageDC.T1>45', 'stageDC.T>25', subset(YX.m4, stage=='DC.T1>45')[, sum(w)], subset(YX.m4, stage=='DC.T>25')[, sum(w)], 1.962),
										Dl20= my.or.from.logit(YX.m4.fit11, 'stageDC.T<=20', 'stageDC.T>25', subset(YX.m4, stage=='DC.T<=20')[, sum(w)], subset(YX.m4, stage=='DC.T>25')[, sum(w)], 1.962),
										Dl22= my.or.from.logit(YX.m4.fit11, 'stageDC.T<=22', 'stageDC.T>25', subset(YX.m4, stage=='DC.T<=20')[, sum(w)], subset(YX.m4, stage=='DC.T>25')[, sum(w)], 1.962),
										Dl24= my.or.from.logit(YX.m4.fit11, 'stageDC.T<=24', 'stageDC.T>25', subset(YX.m4, stage=='DC.T<=25')[, sum(w)], subset(YX.m4, stage=='DC.T>25')[, sum(w)], 1.962))
		#        DAy      DAm    D1g45
		#1: 2.697307 3.302536 1.443833
		#2: 2.131744 2.605657 1.135020
		#3: 3.412917 4.185795 1.836668
		#	D1>45 remains significant
		#        DAy      DAm    D1g45     Dl20     Dl22      Dl24
		#1: 2.749188 3.366780 1.470059 4.696014 1.379099 1.0650515
		#2: 2.165804 2.647792 1.151648 3.887683 1.144275 0.9277367
		#3: 3.489713 4.281004 1.876507 5.672414 1.662112 1.2226903
		#	young ages also remain significant
	}
}
######################################################################################
#	wrapper to betareg to come up with better starting values due to scores close to 0 and 1
project.athena.Fisheretal.betareg<- function(YX.m3, formula, include.colnames, gamlss.BE.limit.u=c( seq(0.95, 0.99, 0.01),seq(0.991, 1, 0.001)), verbose=0)
{
	require(gamlss)
	gamlss.BE.limit.i	<- gamlss.BE.limit.u[1]
	YX.m3b				<- copy(YX.m3)	
	set(YX.m3b, YX.m3b[, which(score.Y>gamlss.BE.limit.i)], 'score.Y', gamlss.BE.limit.i)
	betafit.or			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0)	
	betafit.rr			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0)	
	tryCatch({
				for(i in seq_along(gamlss.BE.limit.u)[-1])
				{
					#print(gamlss.BE.limit.u[i])			
					YX.m3b	<- copy(YX.m3)	
					set(YX.m3b, YX.m3b[, which(score.Y>gamlss.BE.limit.u[i])], 'score.Y', gamlss.BE.limit.u[i])
					tmp.or				<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0, start.from=betafit.or)
					tmp.rr				<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3b[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0, start.from=betafit.rr)
					gamlss.BE.limit		<- gamlss.BE.limit.u[i]	#not run if error in gamlss
					betafit.rr			<- tmp.rr				#not run if error in gamlss
					betafit.or			<- tmp.or				#not run if error in gamlss
				}
				gamlss.init	<- list(start.from.rr=betafit.rr, start.from.or=betafit.or, BE.limit=gamlss.BE.limit)
			}, error=function(e){ print(e$message); gamlss.init<<- list(start.from.rr=betafit.rr, start.from.or=betafit.or, BE.limit=gamlss.BE.limit)})
	
	if(verbose) cat(paste('\nsuccess: gamlss beta regression for score<=',gamlss.init$BE.limit))
	set(YX.m3, YX.m3[, which(score.Y>gamlss.init$BE.limit)], 'score.Y', gamlss.init$BE.limit)
	betafit.or			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), family=BE(mu.link='logit'), trace=0, start.from=gamlss.init$start.from.or)
	betafit.rr			<- gamlss(as.formula(formula), weights=w, data=na.omit(as.data.frame(YX.m3[, include.colnames, with=FALSE])), family=BE(mu.link='log'), trace=0, start.from=gamlss.init$start.from.rr)
	
	list(betafit.or=betafit.or, betafit.rr=betafit.rr, gamlss.BE.limit=gamlss.init$BE.limit)
}
######################################################################################
project.athena.Fisheretal.estimate.risk.core<- function(YX.m3, X.seq, formula, predict.df, risk.df, include.colnames, bs.n=1e3, gamlss.BE.required.limit=0.99 )
{
	require(gamlss)
	options(warn=0)
	#
	cat(paste('\nbeta regression for formula=',formula))	 
	tmp				<- project.athena.Fisheretal.betareg(YX.m3, formula, include.colnames, gamlss.BE.limit.u=c( c(0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1)), verbose=1 )
	betafit.or		<- tmp$betafit.or
	betafit.rr		<- tmp$betafit.rr
	stopifnot(tmp$gamlss.BE.limit>gamlss.BE.required.limit)
		
	#AIC(YX.m3.fit1, YX.m3.fit2, YX.m3.fit3, YX.m3.fit2b, YX.m3.fit4, YX.m3.fit5, YX.m3.fit6, YX.m3.fit7 )
	#AIC(YX.m3.fit1, YX.m3.fit2, YX.m3.fit3, YX.m3.fit2b, YX.m3.fit4, YX.m3.fit5, YX.m3.fit6, YX.m3.fit7, k=YX.m3[, log(round(sum(w)))])
	#sapply( list(YX.m3.fit1, YX.m3.fit2, YX.m3.fit3, YX.m3.fit2b, YX.m3.fit4, YX.m3.fit5, YX.m3.fit6, YX.m3.fit7), '[[', 'pseudo.r.squared')
	#str(X.seq)
	#str(YX.m3)	
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
				tmp2			<- na.omit(merge(tmp2[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3[which( YX.m3[, risk.ref, with=FALSE][[1]]==factor.ref ), c(risk.ref,corisk,'w'), with=FALSE], by=risk.ref))
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
				tmp2			<- na.omit(merge(tmp2[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3[which( YX.m3[, risk.ref, with=FALSE][[1]]==factor.ref ), c(risk.ref,corisk,'w'), with=FALSE], by=risk.ref))
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
	#	person years in infection window
	cat(paste('\nperson years across infection windows'))
	setkey(risk.df, coef)
	if(is.null(X.seq))
	{
		tmp			<- unique(risk.df)
		tmp[, stat:='PY']
		tmp			<- subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref, factor.ref, PY))
		set(tmp, NULL, c('coef.ref','risk.ref', 'factor.ref'), 'None')
		setnames(tmp,'PY','v')
		risk.ans	<- rbind(risk.ans, tmp)		
	}
	if(!is.null(X.seq))
	{		
		tmp			<- unique(risk.df)[, list(risk=risk, factor=factor, risk.ref="None", factor.ref="None", coef.ref="None", stat="PY", v=length(which(as.character(X.seq[, risk, with=FALSE][[1]])==factor))), by='coef']
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
	#	relative infectiousness
	cat(paste('\nrelative infectiousness'))
	tmp			<- subset(risk.ans, stat=='PY' )[, sum(v)]
	tmp			<- subset(risk.ans, risk.ref=='None')[, list( stat='RI', risk=risk[1], factor=factor[1], risk.ref='None', factor.ref='None', coef.ref='None', v= v[stat=='P'] / ( v[stat=='PY'] / tmp ) ), by='coef']
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
						tmp2			<- project.athena.Fisheretal.betareg(YX.m3.bs, formula, include.colnames, gamlss.BE.limit.u=c( c(0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 1)), verbose=1 )
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
																				tmp2			<- na.omit(merge(tmp2[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3.bs[which( YX.m3.bs[, risk.ref, with=FALSE][[1]]==factor.ref ), c(risk.ref,corisk,'w'), with=FALSE], by=risk.ref))
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
																				tmp2			<- na.omit(merge(tmp2[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3.bs[which( YX.m3.bs[, risk.ref, with=FALSE][[1]]==factor.ref ), c(risk.ref,corisk,'w'), with=FALSE], by=risk.ref))							
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
				#	for proportions		
				setkey(risk.df, coef)
				tmp			<- unique(risk.df)[, 		{
							tmp				<- copy(predict.df)
							set(tmp, NULL, risk, factor(factor, levels=levels( predict.df[[risk]] )))
							corisk			<- intersect(colnames(risk.df), colnames(predict.df))
							if(!length(corisk))	
								corisk<- NULL
							tmp				<- na.omit(merge(tmp[,setdiff( colnames(predict.df), c(colnames(risk.df),'w') ), with=FALSE], YX.m3.bs[which( YX.m3.bs[, risk, with=FALSE][[1]]==factor ), c(risk,corisk,'w'), with=FALSE], by=risk))											
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
	tmp			<- risk.ans.bs[, list(stat='N', risk=risk[1], factor=factor[1], risk.ref='None', factor.ref='None', coef.ref='None', v= v[stat=='prob']*v[stat=='n'] ), by=c('bs','coef')]
	set(tmp, tmp[, which(v==0)],'v', NA_real_)
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, v, bs)))
	tmp			<- tmp[,	list(stat='P', coef=coef, risk=risk, factor=factor, risk.ref='None', factor.ref='None', coef.ref='None', v= v/sum(v, na.rm=TRUE)),	by='bs']		
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, v, bs)))	
	tmp			<- merge( subset(tmp, stat=='P'), subset(risk.ans, stat=='PY', select=c(coef, v)), by='coef' )
	tmp[, v:= v.x/( v.y/subset(risk.ans, stat=='PY')[, sum(v)] )]
	set(tmp, NULL, 'stat', 'RI')
	set(tmp, tmp[,which(is.nan(v) | v==0)],'v',NA_real_)
	risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(coef, coef.ref, stat, risk, factor, risk.ref,  factor.ref, v, bs)))
	tmp			<- risk.ans.bs[,	list(l95.bs=quantile(v, prob=0.025, na.rm=TRUE), u95.bs=quantile(v, prob=0.975, na.rm=TRUE)), by=c('coef','coef.ref','stat')]
	risk.ans	<- merge(risk.ans, tmp, by=c('coef','coef.ref','stat'), all.x=TRUE)
	#	
	setkey(risk.ans, stat, coef.ref, coef)
	list(risk=risk.ans, fit.or=betafit.or, fit.rr=betafit.rr, risk.bs=risk.ans.bs)
}
######################################################################################
project.athena.Fisheretal.estimate.risk.wrap<- function(YX, X.tables, plot.file.or=NA, bs.n=1e3, resume=TRUE, save.file=NA, method=NA)
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
		cat(paste('\nregression on data set by method', method))
		if(method%in%c(	'm21st.cas','m21st.cas.clu','m21st.cas.adj','m21st.cas.clu.adj','m21st.cas.cens','m21st.cas.clu.cens',
						'm2B1st.cas','m2B1st.cas.clu','m2B1st.cas.adj','m2B1st.cas.clu.adj','m2B1st.cas.cens','m2B1st.cas.clu.cens','m2B1st.cas.censp','m2B1st.cas.clu.censp'))
		{  
			#	censoring adjustment
			if(grepl('cens', method))
			{
				tmp			<- ifelse(grepl('m2wmx',method),"CD41st.tperiod","CD4t.tperiod")	
				set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
				tmp			<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')	
				if(grepl('censp', method))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			tmp				<- ifelse(grepl('m2wmx',method),"CD41st","CD4t")
			set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
			#	sequence adjustment (not censoring)
			if(grepl('adj', method))
			{
				tmp			<- ifelse(grepl(method, 'clu'), 'adj.clu', 'adj.seq')
				tmp			<- subset( X.tables[[tmp]], factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')			
			predict.df		<- data.table(stage=factor('ART1.su.Y', levels=YX[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)], risk.ref='stage', factor.ref='ART1.su.Y')
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]			
			tmp				<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')
			risk.df			<- merge(risk.df, subset(X.tables$risk.table, stat==tmp, c(risk, factor, n)), by=c('risk','factor'), all.x=1)
			setnames(risk.df, 'n', 'PY')
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
		}	
		if(method%in%c(	'm2t.cas','m2t.cas.clu','m2t.cas.adj','m2t.cas.clu.adj','m2t.cas.cens','m2t.cas.clu.cens',
						'm2Bt.cas','m2Bt.cas.clu','m2Bt.cas.adj','m2Bt.cas.clu.adj','m2Bt.cas.cens','m2Bt.cas.clu.cens','m2Bt.cas.censp','m2Bt.cas.clu.censp'))
		{  
			#	censoring adjustment
			if(grepl('cens', method))
			{
				tmp			<- ifelse(grepl('m2wmx',method),"CD41st.tperiod","CD4t.tperiod")	
				set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
				tmp			<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')	
				if(grepl('censp', method))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			tmp				<- ifelse(grepl('m2wmx',method),"CD41st","CD4t")
			set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
			#	sequence adjustment (not censoring)
			if(grepl('adj', method))
			{
				tmp			<- ifelse(grepl(method, 'clu'), 'adj.clu', 'adj.seq')
				tmp			<- subset( X.tables[[tmp]], factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')			
			predict.df		<- data.table(stage=factor('ART.su.Y', levels=YX[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)], risk.ref='stage', factor.ref='ART.su.Y')
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]			
			tmp				<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')
			risk.df			<- merge(risk.df, subset(X.tables$risk.table, stat==tmp, c(risk, factor, n)), by=c('risk','factor'), all.x=1)
			setnames(risk.df, 'n', 'PY')
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
		}					
		if(method%in%c(	'm2wmx.cas','m2wmx.cas.clu','m2wmx.cas.adj','m2wmx.cas.clu.adj','m2wmx.cas.cens','m2wmx.cas.clu.cens',
						'm2Bwmx.cas','m2Bwmx.cas.clu','m2Bwmx.cas.adj','m2Bwmx.cas.clu.adj','m2Bwmx.cas.cens','m2Bwmx.cas.clu.cens','m2Bwmx.cas.censp','m2Bwmx.cas.clu.censp'))
		{  
			#	censoring adjustment
			if(grepl('cens', method))
			{
				tmp			<- ifelse(grepl('m2wmx',method),"CD41st.tperiod","CD4t.tperiod")	
				set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
				tmp			<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')	
				if(grepl('censp', method))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			tmp				<- ifelse(grepl('m2wmx',method),"CD41st","CD4t")
			set(YX, NULL, 'stage', factor(as.character(YX[[tmp]])))
			#	sequence adjustment (not censoring)
			if(grepl('adj', method))
			{
				tmp			<- ifelse(grepl(method, 'clu'), 'adj.clu', 'adj.seq')
				tmp			<- subset( X.tables[[tmp]], factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')						
			predict.df		<- data.table(stage=factor('ART.suA.Y', levels=YX[, levels(stage)]), w=1.)
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)], risk.ref='stage', factor.ref='ART.suA.Y')
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			tmp				<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')
			risk.df			<- merge(risk.df, subset(X.tables$risk.table, stat==tmp, c(risk, factor, n)), by=c('risk','factor'), all.x=1)
			setnames(risk.df, 'n', 'PY')
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
		}	
		if(grepl('m2wmx.tp', method) | grepl('m2Bwmx.tp', method))
		{
			tp				<- regmatches(method, regexpr('tp[0-9]', method))
			cat(paste('\nprocess time period',tp))
			tp				<- substr(tp, 3, 3)
			YX				<- subset(YX, t.period==tp)						
			tmp				<- ifelse(grepl('m2wmx',method),"CD41st.tperiod","CD4t.tperiod")
			set(YX, NULL, 'stage', YX[[tmp]])						
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')			
			
			set(YX, NULL, 'stage', YX[,factor(as.character(stage))])
			predict.df		<- data.table(stage=factor(paste('ART.suA.Y',tp,sep='.'), levels=YX[, levels(stage)]), w=1.)									
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.suA.Y',tp,sep='.'))
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref=paste('U',tp,sep='.')) )
			tmp				<- YX[, levels(stage)][ substr(YX[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= YX[, levels(stage)][ YX[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			tmp				<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')
			risk.df			<- merge(risk.df, subset(X.tables$risk.table, stat==tmp, c(risk, factor, n)), by=c('risk','factor'), all.x=1)
			setnames(risk.df, 'n', 'PY')
			#	sequence adjustment
			if(grepl('adj', method))
			{
				tmp			<- ifelse(grepl(method, 'clu'), 'adj.clu', 'adj.seq')
				tmp			<- subset( X.tables[[tmp]], factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			#	sequence and censoring adjustment
			if(grepl('cens', method))
			{
				tmp			<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')	
				if(grepl('censp', method))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}	
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )					
		}
		if(method%in%c('m3.i','m3.i.clu'))
		{
			#	TODO
			#	ART indicators only						
			formula			<- 'score.Y ~ ART.pulse+ART.I+ART.F+ART.P+ART.A-1'
			include.colnames<- c('score.Y','w','ART.pulse','ART.I','ART.F','ART.P','ART.A')			
			set(YX, NULL, 'ART.pulse', YX[,factor(as.character(ART.pulse), levels=X.den[, levels(ART.pulse)])])
			set(YX, NULL, 'ART.I', YX[,factor(as.character(ART.I), levels=X.den[, levels(ART.I)])])
			set(YX, NULL, 'ART.F', YX[,factor(as.character(ART.F), levels=X.den[, levels(ART.F)])])
			set(YX, NULL, 'ART.P', YX[,factor(as.character(ART.P), levels=X.den[, levels(ART.P)])])
			set(YX, NULL, 'ART.A', YX[,factor(as.character(ART.A), levels=X.den[, levels(ART.A)])])
			predict.df		<- data.table( 	ART.pulse=factor('No',levels=c('No','Yes')), ART.I=factor('No',levels=c('No','Yes')), ART.F=factor('No',levels=c('No','Yes')), ART.P=factor('No',levels=c('No','Yes')), ART.A=factor('No',levels=c('No','Yes')),
											w=1.)
			risk.df			<- data.table( risk= c('ART.pulse','ART.I','ART.F','ART.P','ART.A'), factor=rep('Yes',5) )
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			tmp				<- data.table(risk.ref= rep('ART.pulse',nrow(risk.df)), factor.ref= rep('No',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
			ans$risk.table	<- do.call('rbind',list(
										risk.df[,	{
												z	<- table( YX[, risk, with=FALSE], useNA='ifany')
												list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
											},by='risk'],
							 			risk.df[,	{
												z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
												list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
											},by='risk']))			
		}
		if(method%in%c('m3.ni','m3.ni.clu'))
		{
			#	TODO
			#	ART indicators and number of drugs
			#	ART.nDrug.c independent of indicators			
			set(YX, NULL, 'stage', YX[, ART.nDrug.c])
			set(X.den, NULL, 'stage', X.den[, ART.nDrug.c])
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			set(YX, NULL, 'ART.pulse', YX[,factor(as.character(ART.pulse), levels=X.den[, levels(ART.pulse)])])
			set(YX, NULL, 'ART.I', YX[,factor(as.character(ART.I), levels=X.den[, levels(ART.I)])])
			set(YX, NULL, 'ART.F', YX[,factor(as.character(ART.F), levels=X.den[, levels(ART.F)])])
			set(YX, NULL, 'ART.P', YX[,factor(as.character(ART.P), levels=X.den[, levels(ART.P)])])
			set(YX, NULL, 'ART.A', YX[,factor(as.character(ART.A), levels=X.den[, levels(ART.A)])])			
			formula			<- 'score.Y ~ stage+ART.pulse+ART.I+ART.F+ART.P+ART.A-1'
			include.colnames<- c('score.Y','w','stage','ART.pulse','ART.I','ART.F','ART.P','ART.A')
			predict.df		<- data.table(	stage=factor('ART.3', levels=X.den[, levels(stage)]), 
											ART.pulse=factor('No',levels=c('No','Yes')), ART.I=factor('No',levels=c('No','Yes')), ART.F=factor('No',levels=c('No','Yes')), ART.P=factor('No',levels=c('No','Yes')), ART.A=factor('No',levels=c('No','Yes')),
											w=1.)
			risk.df			<- data.table( risk= c(rep('stage',3), 'ART.pulse','ART.I','ART.F','ART.P','ART.A'), factor=c('ART.3','ART.l3','ART.g3', rep('Yes',5)))
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))				
		}
		if(method%in%c('m3.nic','m3.nic.clu','m3.nic.adj','m3.nic.clu.adj','m3.nic.cens','m3.nic.clu.cens','m3.nic.censp','m3.nic.clu.censp'))
		{
			#	number of drugs conditional on no indicators and ART indicators (not overlapping)  
			#	censoring adjustment
			if(grepl('cens', method))
			{				
				YX[, ART.nstage.c.tperiod:= YX[, paste(ART.nstage.c, t.period, sep='.')]]
				set(YX, YX[,which(is.na(ART.nstage.c))], 'ART.nstage.c.tperiod', NA_character_)
				set(YX, NULL, 'stage', YX[, factor(as.character(ART.nstage.c.tperiod))])
				tmp			<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')
				if(grepl('censp', method))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			set(YX, NULL, 'stage', YX[, factor(as.character(ART.nstage.c))])
			#	sequence adjustment (not censoring)
			if(grepl('adj', method))
			{
				tmp			<- ifelse(grepl(method, 'clu'), 'adj.clu', 'adj.seq')
				tmp			<- subset( X.tables[[tmp]], factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			predict.df		<- data.table(stage=factor('ART.3', levels=YX[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)])
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])
			tmp				<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')
			risk.df			<- merge(risk.df, subset(X.tables$risk.table, stat==tmp, c(risk, factor, n)), by=c('risk','factor'), all.x=1)
			setnames(risk.df, 'n', 'PY')
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
		}
		if(method%in%c('m3.nicv','m3.nicv.clu','m3.nicv.adj','m3.nicv.clu.adj','m3.nicv.cens','m3.nicv.clu.cens','m3.nicv.censp','m3.nicv.clu.censp'))
		{
			#	number of drugs conditional on no indicators and ART indicators (not overlapping)  
			#	censoring adjustment
			if(grepl('cens', method))
			{				
				YX[, ART.nstage.c.tperiod:= YX[, paste(ART.nstage.c, t.period, sep='.')]]
				set(YX, YX[,which(is.na(ART.nstage.c))], 'ART.nstage.c.tperiod', NA_character_)
				set(YX, NULL, 'stage', YX[, factor(as.character(ART.nstage.c.tperiod))])
				tmp			<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')				
				if(grepl('censp', method))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			set(YX, NULL, 'stage', YX[, factor(as.character(ART.nstage.c))])
			#	sequence adjustment (not censoring)
			if(grepl('adj', method))
			{
				tmp			<- ifelse(grepl(method, 'clu'), 'adj.clu', 'adj.seq')
				tmp			<- subset( X.tables[[tmp]], factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
			include.colnames<- c('score.Y','w','lRNA.mx','stage')
			predict.df		<- data.table(	stage=factor('ART.3', levels=YX[, levels(stage)]), 
											lRNA.mx=subset(YX, stage=='ART.3')[, mean(lRNA.mx, na.rm=TRUE)],
											w=1.)			
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)])
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			risk.df			<- merge(risk.df, risk.df[, {
															tmp				<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
															list(lRNA.mx=ifelse(is.nan(tmp), YX[, mean( lRNA.mx, na.rm=TRUE )], tmp))
														}, by='coef'], by='coef')			
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])
			tmp				<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')
			risk.df			<- merge(risk.df, subset(X.tables$risk.table, stat==tmp, c(risk, factor, n)), by=c('risk','factor'), all.x=1)
			setnames(risk.df, 'n', 'PY')
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
		}		
		if(method%in%c('m3.tni','m3.tni.clu'))
		{
			#	TODO
			#	ART indicators and number/type of drugs
			#	ART.tnDrug.c independent of indicators
			set(YX, NULL, 'stage', YX[, ART.tnDrug.c])
			set(X.den, NULL, 'stage', X.den[, ART.tnDrug.c])
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			set(YX, NULL, 'ART.pulse', YX[,factor(as.character(ART.pulse), levels=X.den[, levels(ART.pulse)])])
			set(YX, NULL, 'ART.I', YX[,factor(as.character(ART.I), levels=X.den[, levels(ART.I)])])
			set(YX, NULL, 'ART.F', YX[,factor(as.character(ART.F), levels=X.den[, levels(ART.F)])])
			set(YX, NULL, 'ART.P', YX[,factor(as.character(ART.P), levels=X.den[, levels(ART.P)])])
			set(YX, NULL, 'ART.A', YX[,factor(as.character(ART.A), levels=X.den[, levels(ART.A)])])						
			formula			<- 'score.Y ~ stage+ART.pulse+ART.I+ART.F+ART.P+ART.A-1'
			include.colnames<- c('score.Y','w','stage','ART.pulse','ART.I','ART.F','ART.P','ART.A')
			predict.df		<- data.table( 	stage=factor('ART.3.NRT.PI', levels=X.den[, levels(stage)]), 
											ART.pulse=factor('No',levels=c('No','Yes')), ART.I=factor('No',levels=c('No','Yes')), ART.F=factor('No',levels=c('No','Yes')), ART.P=factor('No',levels=c('No','Yes')), ART.A=factor('No',levels=c('No','Yes')),
											w=1.)
			risk.df			<- data.table( risk= c(rep('stage',9), 'ART.pulse','ART.I','ART.F','ART.P','ART.A'), factor=c('ART.3.NRT.PI','ART.l3','ART.g3','ART.3.NRT','ART.3.PI','ART.3.NNRT','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NNRT.PI', rep('Yes',5)))
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))
		}
		if(method%in%c('m3.tniv','m3.tniv.clu'))
		{
			#	TODO
			#	ART indicators and number/type of drugs
			#	ART.tnDrug.c independent of indicators
			set(YX, NULL, 'stage', YX[, ART.tnDrug.c])
			set(X.den, NULL, 'stage', X.den[, ART.tnDrug.c])	
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			set(YX, NULL, 'ART.pulse', YX[,factor(as.character(ART.pulse), levels=X.den[, levels(ART.pulse)])])
			set(YX, NULL, 'ART.I', YX[,factor(as.character(ART.I), levels=X.den[, levels(ART.I)])])
			set(YX, NULL, 'ART.F', YX[,factor(as.character(ART.F), levels=X.den[, levels(ART.F)])])
			set(YX, NULL, 'ART.P', YX[,factor(as.character(ART.P), levels=X.den[, levels(ART.P)])])
			set(YX, NULL, 'ART.A', YX[,factor(as.character(ART.A), levels=X.den[, levels(ART.A)])])						
			formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage+ART.pulse+ART.I+ART.F+ART.P+ART.A-1'
			include.colnames<- c('score.Y','w','stage','lRNA.mx','ART.pulse','ART.I','ART.F','ART.P','ART.A')
			predict.df		<- data.table( 	stage=factor('ART.3.NRT.PI', levels=X.den[, levels(stage)]), 
											lRNA.mx=subset(YX, stage=='ART.3.NRT.PI')[, mean(lRNA.mx, na.rm=TRUE)],
											ART.pulse=factor('No',levels=c('No','Yes')), ART.I=factor('No',levels=c('No','Yes')), ART.F=factor('No',levels=c('No','Yes')), ART.P=factor('No',levels=c('No','Yes')), ART.A=factor('No',levels=c('No','Yes')),
											w=1.)
			risk.df			<- data.table( risk= c(rep('stage',9), 'ART.pulse','ART.I','ART.F','ART.P','ART.A'), factor=c('ART.3.NRT.PI','ART.l3','ART.g3','ART.3.NRT','ART.3.PI','ART.3.NNRT','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NNRT.PI', rep('Yes',5)) )
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]			
			risk.df			<- merge(risk.df, risk.df[, {
															tmp				<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
															list(lRNA.mx=ifelse(is.nan(tmp), YX[, mean( lRNA.mx, na.rm=TRUE )], tmp))
														}, by='coef'], by='coef')			
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))
		}
		if(method%in%c('m3.tnic','m3.tnic.clu','m3.tnic.adj','m3.tnic.clu.adj','m3.tnic.cens','m3.tnic.clu.cens','m3.tnic.censp','m3.tnic.clu.censp'))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			#	censoring adjustment
			if(grepl('cens', method))
			{				
				YX[, ART.ntstage.c.tperiod:= YX[, paste(ART.ntstage.c, t.period, sep='.')]]
				set(YX, YX[,which(is.na(ART.ntstage.c))], 'ART.ntstage.c.tperiod', NA_character_)
				set(YX, NULL, 'stage', YX[, factor(as.character(ART.ntstage.c.tperiod))])
				tmp			<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')				
				if(grepl('censp', method))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			set(YX, NULL, 'stage', YX[, factor(as.character(ART.ntstage.c))])
			#	sequence adjustment (not censoring)
			if(grepl('adj', method))
			{
				tmp			<- ifelse(grepl(method, 'clu'), 'adj.clu', 'adj.seq')
				tmp			<- subset( X.tables[[tmp]], factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			predict.df		<- data.table(stage=factor('ART.3.NRT.PI', levels=YX[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)])
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])
			tmp				<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')
			risk.df			<- merge(risk.df, subset(X.tables$risk.table, stat==tmp, c(risk, factor, n)), by=c('risk','factor'), all.x=1)
			setnames(risk.df, 'n', 'PY')
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )						
		}
		if(method%in%c('m3.tnicv','m3.tnicv.clu','m3.tnicv.adj','m3.tnicv.clu.adj','m3.tnicv.cens','m3.tnicv.clu.cens','m3.tnicv.censp','m3.tnicv.clu.censp'))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			#	censoring adjustment
			if(grepl('cens', method))
			{				
				YX[, ART.ntstage.c.tperiod:= YX[, paste(ART.ntstage.c, t.period, sep='.')]]
				set(YX, YX[,which(is.na(ART.ntstage.c))], 'ART.ntstage.c.tperiod', NA_character_)
				set(YX, NULL, 'stage', YX[, factor(as.character(ART.ntstage.c.tperiod))])
				tmp			<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')				
				if(grepl('censp', method))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			set(YX, NULL, 'stage', YX[, factor(as.character(ART.ntstage.c))])
			#	sequence adjustment (not censoring)
			if(grepl('adj', method))
			{
				tmp			<- ifelse(grepl(method, 'clu'), 'adj.clu', 'adj.seq')
				tmp			<- subset( X.tables[[tmp]], factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
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
			tmp				<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')
			risk.df			<- merge(risk.df, subset(X.tables$risk.table, stat==tmp, c(risk, factor, n)), by=c('risk','factor'), all.x=1)
			setnames(risk.df, 'n', 'PY')
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
		}		
		if(method%in%c('m3.tnicvNo','m3.tnicvNo.clu','m3.tnicvNo.adj','m3.tnicvNo.clu.adj','m3.tnicvNo.cens','m3.tnicvNo.clu.cens','m3.tnicvNo.censp','m3.tnicvNo.clu.censp'))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			x
			#	censoring adjustment
			if(grepl('cens', method))
			{				
				YX[, ART.ntstage.no.c.tperiod:= YX[, paste(ART.ntstage.no.c, t.period, sep='.')]]
				set(YX, YX[,which(is.na(ART.ntstage.no.c))], 'ART.ntstage.no.c.tperiod', NA_character_)
				set(YX, NULL, 'stage', YX[, factor(as.character(ART.ntstage.no.c.tperiod))])
				YX			<- subset(YX, !is.na(stage) & !is.na(lRNA.mx))
				tmp			<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')				
				if(grepl('censp', method))
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyPU[stat=='X.msm']/p.adjbyPU[stat==tmp]),by=c('risk','factor')]
				else
					tmp		<- X.tables$cens.table[, list(w.b= p.adjbyNU[stat=='X.msm']/p.adjbyNU[stat==tmp]),by=c('risk','factor')]
				tmp			<- subset( tmp, factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}
			set(YX, NULL, 'stage', YX[, factor(as.character(ART.ntstage.no.c))])
			YX				<- subset(YX, !is.na(stage) & !is.na(lRNA.mx))
			#	sequence adjustment (not censoring)
			if(grepl('adj', method))
			{
				tmp			<- ifelse(grepl(method, 'clu'), 'adj.clu', 'adj.seq')
				tmp			<- subset( X.tables[[tmp]], factor%in%YX[, levels(stage)] )				
				YX			<- merge( YX, data.table( stage=factor( tmp[, factor], levels=YX[, levels(stage)] ), w.b=tmp[, w.b] ), by='stage' )
				set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			}												
			#
			formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
			include.colnames<- c('score.Y','w','stage','lRNA.mx')
			predict.df		<- data.table( 	stage=factor('ART.3.NRT.PI', levels=YX[, levels(stage)]), 
											lRNA.mx=subset(YX, stage=='ART.3.NRT.PI')[, mean(lRNA.mx, na.rm=TRUE)], w=1.	)
			risk.df			<- data.table( risk= rep('stage',14), factor=c('ART.3.NRT.PI','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A', 'ART.F','ART.l3','ART.g3','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NRT','ART.3.NNRT.PI','ART.3.PI','ART.3.NNRT'))
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]			
			risk.df			<- merge(risk.df, risk.df[, {
								tmp				<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
								list(lRNA.mx=ifelse(is.nan(tmp), YX[, mean( lRNA.mx, na.rm=TRUE )], tmp))
							}, by='coef'], by='coef')			
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])			
			tmp				<- ifelse(grepl(method, 'clu'), 'X.clu', 'X.seq')
			risk.df			<- merge(risk.df, subset(X.tables$risk.table, stat==tmp, c(risk, factor, n)), by=c('risk','factor'), all.x=1)
			setnames(risk.df, 'n', 'PY')
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, NULL, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
		}			
		if(!is.na(save.file))
		{
			ans$YX			<- YX
			ans$risk.table	<- X.tables$risk.table
			ans$cens.table	<- X.tables$cens.table
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
project.athena.Fisheretal.v1.estimate.risk.wrap<- function(YX, X.den, df.all, df.viro, X.msm=NULL, X.clu=NULL, df.all.allmsm=NULL, df.viro.allmsm=NULL, plot.file.or=NA, bs.n=1e3, resume=TRUE, save.file=NA, method=NA)
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
		stopifnot(method%in%c(	'm21st.cas','m2t.cas','m2wmx.cas','m2wmx.tp1','m2wmx.tp2', 'm2wmx.tp3','m2wmx.tp4',
						'm2B1st.cas','m2Bt.cas','m2Bwmx.cas', 'm2Bwmx.tp1', 'm2Bwmx.tp2', 'm2Bwmx.tp3', 'm2Bwmx.tp4',
						'm21st.cas.clu','m2t.cas.clu','m2wmx.cas.clu','m2wmx.tp1.clu','m2wmx.tp2.clu','m2wmx.tp3.clu','m2wmx.tp4.clu',
						'm2B1st.cas.clu','m2Bt.cas.clu','m2Bwmx.cas.clu','m2Bwmx.tp1.clu', 'm2Bwmx.tp2.clu', 'm2Bwmx.tp3.clu', 'm2Bwmx.tp4.clu',
						'm21st.cas.adj','m2t.cas.adj','m2wmx.cas.adj','m2wmx.tp1.adj','m2wmx.tp2.adj','m2wmx.tp3.adj','m2wmx.tp4.adj',
						'm2B1st.cas.adj','m2Bt.cas.adj','m2Bwmx.cas.adj','m2Bwmx.tp1.adj','m2Bwmx.tp2.adj','m2Bwmx.tp3.adj','m2Bwmx.tp4.adj',
						'm21st.cas.clu.adj','m2t.cas.clu.adj','m2wmx.cas.clu.adj','m2wmx.tp1.clu.adj','m2wmx.tp2.clu.adj','m2wmx.tp3.clu.adj','m2wmx.tp4.clu.adj',
						'm2B1st.cas.clu.adj','m2Bt.cas.clu.adj','m2Bwmx.cas.clu.adj','m2Bwmx.tp1.clu.adj', 'm2Bwmx.tp2.clu.adj', 'm2Bwmx.tp3.clu.adj', 'm2Bwmx.tp4.clu.adj',								
						'm2Bwmx.tp1.cens','m2Bwmx.tp2.cens','m2Bwmx.tp3.cens','m2Bwmx.tp4.cens',
						'm2Bwmx.tp1.clu.cens','m2Bwmx.tp2.clu.cens','m2Bwmx.tp3.clu.cens','m2Bwmx.tp4.clu.cens',
						'm3.i','m3.ni','m3.nic','m3.tni','m3.tniv','m3.tnic','m3.nicv','m3.tnicv','m3.tnicvNo',								
						'm3.i.clu','m3.ni.clu','m3.nic.clu','m3.tni.clu','m3.tniv.clu','m3.tnic.clu','m3.nicv.clu','m3.tnicv.clu','m3.tnicvNo.clu',
						'm3.nic.clu.adj','m3.nicv.clu.adj','m3.tnic.clu.adj','m3.tnicv.clu.adj','m3.tnicvNo.clu.adj',
						'm3.nic.adj','m3.nicv.adj','m3.tnic.adj','m3.tnicv.adj','m3.tnicvNo.adj'))
		stopifnot( grepl('clu',method)==!is.null(X.clu) )		
		cat(paste('\nregression on data set by method', method))
		if(method%in%c('m21st.cas','m21st.cas.clu'))
		{  
			set(YX, NULL, 'stage', YX[, CD41st])
			set(X.den, NULL, 'stage', X.den[, CD41st])
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			set(YX, NULL, 'stage', YX[, factor(as.character(stage), levels=X.den[, levels(stage)])])
			predict.df		<- data.table(stage=factor('ART1.su.Y', levels=X.den[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref='ART1.su.Y')
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))	
		}
		if(method%in%c('m21st.cas.adj','m21st.cas.clu.adj'))
		{  
			set(YX, NULL, 'stage', YX[, CD41st])
			set(X.den, NULL, 'stage', X.den[, CD41st])
			set(X.msm, NULL, 'stage', X.msm[, CD41st])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])						
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#	adjust weight by selection bias weights	
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, CD41st])
				set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])												
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.clu[, table(stage)] / nrow(X.clu) )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.den[, table(stage)] / nrow(X.den) )
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#						
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')			
			predict.df		<- data.table(stage=factor('ART1.su.Y', levels=X.den[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref='ART1.su.Y')
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))		
		}
		if(method%in%c('m2B1st.cas','m2B1st.cas.clu'))
		{  
			set(YX, NULL, 'stage', YX[, CD4t])
			set(X.den, NULL, 'stage', X.den[, CD4t])
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			set(YX, NULL, 'stage', YX[, factor(as.character(stage), levels=X.den[, levels(stage)])])
			predict.df		<- data.table(stage=factor('ART1.su.Y', levels=X.den[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref='ART1.su.Y')
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))	
		}
		if(method%in%c('m2B1st.cas.adj','m2B1st.cas.clu.adj'))
		{  
			set(YX, NULL, 'stage', YX[, CD4t])
			set(X.den, NULL, 'stage', X.den[, CD4t])
			set(X.msm, NULL, 'stage', X.msm[, CD4t])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])						
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#	adjust weight by selection bias weights	
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, CD4t])
				set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])												
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.clu[, table(stage)] / nrow(X.clu) )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.den[, table(stage)] / nrow(X.den) )							
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#						
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')			
			predict.df		<- data.table(stage=factor('ART1.su.Y', levels=X.den[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref='ART1.su.Y')
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))		
		}		
		if(method%in%c('m2t.cas','m2t.cas.clu'))
		{  
			set(YX, NULL, 'stage', YX[, CD41st])
			set(X.den, NULL, 'stage', X.den[, CD41st])
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			set(YX, NULL, 'stage', YX[, factor(as.character(stage), levels=X.den[, levels(stage)])])
			predict.df		<- data.table(stage=factor('ART.su.Y', levels=X.den[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref='ART.su.Y')
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))	
		}
		if(method%in%c('m2t.cas.adj','m2t.cas.clu.adj'))
		{  
			set(YX, NULL, 'stage', YX[, CD41st])
			set(X.den, NULL, 'stage', X.den[, CD41st])
			set(X.msm, NULL, 'stage', X.msm[, CD41st])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])						
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#	adjust weight by selection bias weights	
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, CD41st])
				set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])												
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.clu[, table(stage)] / nrow(X.clu) )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.den[, table(stage)] / nrow(X.den) )			
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#			
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')			
			predict.df		<- data.table(stage=factor('ART.su.Y', levels=YX[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref='ART.su.Y')
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))			
		}
		if(method%in%c('m2Bt.cas','m2Bt.cas.clu'))
		{  
			set(YX, NULL, 'stage', YX[, CD4t])
			set(X.den, NULL, 'stage', X.den[, CD4t])
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			set(YX, NULL, 'stage', YX[, factor(as.character(stage), levels=X.den[, levels(stage)])])
			predict.df		<- data.table(stage=factor('ART.su.Y', levels=X.den[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref='ART.su.Y')
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))	
		}
		if(method%in%c('m2Bt.cas.adj','m2Bt.cas.clu.adj'))
		{  
			set(YX, NULL, 'stage', YX[, CD4t])
			set(X.den, NULL, 'stage', X.den[, CD4t])
			set(X.msm, NULL, 'stage', X.msm[, CD4t])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])						
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#	adjust weight by selection bias weights
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, CD4t])
				set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])												
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.clu[, table(stage)] / nrow(X.clu) )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.den[, table(stage)] / nrow(X.den) )			
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#			
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')			
			predict.df		<- data.table(stage=factor('ART.su.Y', levels=YX[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref='ART.su.Y')
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))			
		}		
		if(method%in%c('m2wmx.cas','m2wmx.cas.clu'))
		{  
			set(YX, NULL, 'stage', YX[, CD41st])
			set(X.den, NULL, 'stage', X.den[, CD41st])
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			set(YX, NULL, 'stage', YX[, factor(as.character(stage), levels=X.den[, levels(stage)])])			
			predict.df		<- data.table(stage=factor('ART.suA.Y', levels=X.den[, levels(stage)]), w=1.)
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref='ART.suA.Y')
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))						
		}
		if(method%in%c('m2wmx.cas.adj','m2wmx.cas.clu.adj'))
		{  
			set(YX, NULL, 'stage', YX[, CD41st])
			set(X.den, NULL, 'stage', X.den[, CD41st])
			set(X.msm, NULL, 'stage', X.msm[, CD41st])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])						
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#	adjust weight by selection bias weights			
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, CD41st])
				set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])												
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.clu[, table(stage)] / nrow(X.clu) )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.den[, table(stage)] / nrow(X.den) )
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')						
			predict.df		<- data.table(stage=factor('ART.suA.Y', levels=YX[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref='ART.suA.Y')
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))
		}
		if(method%in%c('m2Bwmx.cas','m2Bwmx.cas.clu'))
		{  
			set(YX, NULL, 'stage', YX[, CD4t])
			set(X.den, NULL, 'stage', X.den[, CD4t])
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			set(YX, NULL, 'stage', YX[, factor(as.character(stage), levels=X.den[, levels(stage)])])			
			predict.df		<- data.table(stage=factor('ART.suA.Y', levels=X.den[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref='ART.suA.Y')
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))						
		}
		if(method%in%c('m2Bwmx.cas.adj','m2Bwmx.cas.clu.adj'))
		{  
			set(YX, NULL, 'stage', YX[, CD4t])
			set(X.den, NULL, 'stage', X.den[, CD4t])
			set(X.msm, NULL, 'stage', X.msm[, CD4t])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])						
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#	adjust weight by selection bias weights			
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, CD4t])
				set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])								
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.clu[, table(stage)] / nrow(X.clu) )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.den[, table(stage)] / nrow(X.den) )
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')						
			predict.df		<- data.table(stage=factor('ART.suA.Y', levels=YX[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref='ART.suA.Y')
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref='U') )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))
		}		
		if(grepl('m2wmx.tp', method) & !grepl('adj', method) & !grepl('cens', method))
		{
			tp				<- regmatches(method, regexpr('tp[0-9]', method))
			cat(paste('\nprocess time period',tp))
			tp				<- substr(tp, 3, 3)
			YX				<- subset(YX, t.period==tp)						
			X.den			<- subset(X.den, t.period==tp)
			
			set(YX, NULL, 'stage', YX[, CD41st.tperiod])
			set(X.den, NULL, 'stage', X.den[, CD41st.tperiod])			
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')			
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage))])
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			predict.df		<- data.table(stage=factor(paste('ART.suA.Y',tp,sep='.'), levels=X.den[, levels(stage)]), w=1.)									
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.suA.Y',tp,sep='.'))
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref=paste('U',tp,sep='.')) )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))		
		}
		if(grepl('m2wmx.tp', method) & grepl('adj', method))
		{
			tp				<- regmatches(method, regexpr('tp[0-9]', method))
			cat(paste('\nprocess time period',tp))
			tp				<- substr(tp, 3, 3)
			YX				<- subset(YX, t.period==tp)						
			X.den			<- subset(X.den, t.period==tp)
			X.msm			<- subset(X.msm, t.period==tp)
			gc()
			set(YX, NULL, 'stage', YX[, CD41st.tperiod])
			set(X.den, NULL, 'stage', X.den[, CD41st.tperiod])	
			set(X.msm, NULL, 'stage', X.msm[, CD41st.tperiod])	
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])
			set(X.den, NULL, 'stage', X.den[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(YX, NULL, 'stage', YX[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#	adjust weight by selection bias weights			
			if(grepl('clu',method))
			{
				X.clu		<- subset(X.clu, t.period==tp)
				set(X.clu, NULL, 'stage', X.clu[, CD41st.tperiod])
				set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])								
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.clu[, table(stage)] / nrow(X.clu) )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.den[, table(stage)] / nrow(X.den) )
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#
			predict.df		<- data.table(stage=factor(paste('ART.suA.Y',tp,sep='.'), levels=X.msm[, levels(stage)]), w=1.)									
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.suA.Y',tp,sep='.'))
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref=paste('U',tp,sep='.')) )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))				
		}
		if(grepl('m2wmx.tp', method) & grepl('cens', method))
		{
			tp				<- regmatches(method, regexpr('tp[0-9]', method))
			cat(paste('\nprocess time period',tp))
			tp				<- substr(tp, 3, 3)
			set(YX, NULL, 'stage', YX[, CD41st.tperiod])
			set(X.den, NULL, 'stage', X.den[, CD41st.tperiod])	
			set(X.msm, NULL, 'stage', X.msm[, CD41st.tperiod])	
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])
			set(X.den, NULL, 'stage', X.den[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(YX, NULL, 'stage', YX[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, CD41st.tperiod])
				set(X.clu, NULL, 'stage', X.clu[, factor(as.character(stage), levels=X.msm[, levels(stage)])])												
			}
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.suA.Y',tp,sep='.'))
			cens.table		<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(stage=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										if(grepl('clu',method))
											z	<- table( X.clu[, risk, with=FALSE], useNA='ifany')
										else
											z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(stage=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(stage=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))		
			cens.table[, t.period:=cens.table[, substr(stage, nchar(stage), nchar(stage))]]
			cens.table[, stage2:=cens.table[, substr(stage, 1, nchar(stage)-2)]]
			for(f in c('U','UAy','UAm'))
				for(z in cens.table[, unique(t.period)])	
					set(cens.table, cens.table[, which(stage2==f & stat=='X.msm' & t.period==z)], 'n',  cens.table[which(stage2==f & stat=='X.msm' & t.period%in%c('2',z)), max(n)])
			cens.table		<- merge(cens.table, cens.table[, list(stage=stage, p= n/sum(n, na.rm=TRUE)), by=c('stat','t.period')], by=c('stat','t.period','stage'))
			YX				<- subset(YX, t.period==tp)						
			X.den			<- subset(X.den, t.period==tp)
			X.msm			<- subset(X.msm, t.period==tp)
			cens.table		<- subset(cens.table, t.period==tp)
			gc()
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			bias.adj		<- cens.table[, list(w.b= p[stat=='X.msm']/p[stat=='X']),by='stage']
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])
			set(X.den, NULL, 'stage', X.den[, factor(as.character(stage), levels=X.msm[, levels(stage)])])						
			set(YX, NULL, 'stage', YX[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#
			predict.df		<- data.table(stage=factor(paste('ART.suA.Y',tp,sep='.'), levels=X.msm[, levels(stage)]), w=1.)									
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.suA.Y',tp,sep='.'))
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref=paste('U',tp,sep='.')) )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$cens.table	<- cens.table
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))				
		}
		if(grepl('m2Bwmx.tp', method) & !grepl('adj', method) & !grepl('cens', method))
		{
			tp				<- regmatches(method, regexpr('tp[0-9]', method))
			cat(paste('\nprocess time period',tp))
			tp				<- substr(tp, 3, 3)
			YX				<- subset(YX, t.period==tp)						
			X.den			<- subset(X.den, t.period==tp)
			
			set(YX, NULL, 'stage', YX[, CD4t.tperiod])
			set(X.den, NULL, 'stage', X.den[, CD4t.tperiod])			
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')			
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage))])
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			predict.df		<- data.table(stage=factor(paste('ART.suA.Y',tp,sep='.'), levels=X.den[, levels(stage)]), w=1.)									
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.suA.Y',tp,sep='.'))
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref=paste('U',tp,sep='.')) )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))		
		}
		if(grepl('m2Bwmx.tp', method) & grepl('adj', method))
		{
			tp				<- regmatches(method, regexpr('tp[0-9]', method))
			cat(paste('\nprocess time period',tp))
			tp				<- substr(tp, 3, 3)
			YX				<- subset(YX, t.period==tp)						
			X.den			<- subset(X.den, t.period==tp)
			X.msm			<- subset(X.msm, t.period==tp)
			set(YX, NULL, 'stage', YX[, CD4t.tperiod])
			set(X.den, NULL, 'stage', X.den[, CD4t.tperiod])
			set(X.msm, NULL, 'stage', X.msm[, CD4t.tperiod])
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])
			set(X.den, NULL, 'stage', X.den[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(YX, NULL, 'stage', YX[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#	adjust weight by selection bias weights			
			if(grepl('clu',method))
			{
				X.clu		<- subset(X.clu, t.period==tp)
				set(X.clu, NULL, 'stage', X.clu[, CD4t.tperiod])
				set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])				
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.clu[, table(stage)] / nrow(X.clu) )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.den[, table(stage)] / nrow(X.den) )
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#
			predict.df		<- data.table(stage=factor(paste('ART.suA.Y',tp,sep='.'), levels=X.msm[, levels(stage)]), w=1.)									
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.suA.Y',tp,sep='.'))
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref=paste('U',tp,sep='.')) )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))				
		}
		if(grepl('m2Bwmx.tp', method) & grepl('cens', method))
		{
			tp				<- regmatches(method, regexpr('tp[0-9]', method))
			cat(paste('\nprocess time period',tp))
			tp				<- substr(tp, 3, 3)
			set(YX, NULL, 'stage', YX[, CD4t.tperiod])
			set(X.den, NULL, 'stage', X.den[, CD4t.tperiod])	
			set(X.msm, NULL, 'stage', X.msm[, CD4t.tperiod])	
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])
			set(X.den, NULL, 'stage', X.den[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(YX, NULL, 'stage', YX[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, CD4t.tperiod])
				set(X.clu, NULL, 'stage', X.clu[, factor(as.character(stage), levels=X.msm[, levels(stage)])])												
			}
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.suA.Y',tp,sep='.'))
			cens.table		<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(stage=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										if(grepl('clu',method))
											z	<- table( X.clu[, risk, with=FALSE], useNA='ifany')
										else
											z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(stage=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(stage=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))		
			cens.table[, t.period:=cens.table[, substr(stage, nchar(stage), nchar(stage))]]
			cens.table[, stage2:=cens.table[, substr(stage, 1, nchar(stage)-2)]]
			for(f in c('U','UAy','UAm'))
				for(z in cens.table[, unique(t.period)])	
					set(cens.table, cens.table[, which(stage2==f & stat=='X.msm' & t.period==z)], 'n',  cens.table[which(stage2==f & stat=='X.msm' & t.period%in%c('2',z)), max(n)])
			cens.table		<- merge(cens.table, cens.table[, list(stage=stage, p= n/sum(n, na.rm=TRUE)), by=c('stat','t.period')], by=c('stat','t.period','stage'))
			YX				<- subset(YX, t.period==tp)						
			X.den			<- subset(X.den, t.period==tp)
			X.msm			<- subset(X.msm, t.period==tp)
			cens.table		<- subset(cens.table, t.period==tp)
			gc()
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			#	adjust weight by selection bias weights			
			bias.adj		<- cens.table[, list(w.b= p[stat=='X.msm']/p[stat=='X']),by='stage']
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])
			set(X.den, NULL, 'stage', X.den[, factor(as.character(stage), levels=X.msm[, levels(stage)])])						
			set(YX, NULL, 'stage', YX[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#
			predict.df		<- data.table(stage=factor(paste('ART.suA.Y',tp,sep='.'), levels=X.msm[, levels(stage)]), w=1.)									
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.suA.Y',tp,sep='.'))
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='D' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref='stage', factor.ref=paste('U',tp,sep='.')) )
			tmp				<- X.den[, levels(stage)][ substr(X.den[, levels(stage)],1,1)=='A' ]
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('l500', levels(stage))] ] )	)
			risk.df			<- rbind(risk.df, data.table(risk='stage',factor=tmp, risk.ref= 'stage', factor.ref= X.den[, levels(stage)][ X.den[, substr(levels(stage),1,1)=='D' & grepl('g500', levels(stage))] ] )	)
			risk.df			<- risk.df[, list(coef=paste(risk, factor,sep=''), coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk','factor','risk.ref','factor.ref')]
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$cens.table	<- cens.table
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))				
		}
		if(method%in%c('m3.i','m3.i.clu'))
		{
			#	ART indicators only						
			formula			<- 'score.Y ~ ART.pulse+ART.I+ART.F+ART.P+ART.A-1'
			include.colnames<- c('score.Y','w','ART.pulse','ART.I','ART.F','ART.P','ART.A')			
			set(YX, NULL, 'ART.pulse', YX[,factor(as.character(ART.pulse), levels=X.den[, levels(ART.pulse)])])
			set(YX, NULL, 'ART.I', YX[,factor(as.character(ART.I), levels=X.den[, levels(ART.I)])])
			set(YX, NULL, 'ART.F', YX[,factor(as.character(ART.F), levels=X.den[, levels(ART.F)])])
			set(YX, NULL, 'ART.P', YX[,factor(as.character(ART.P), levels=X.den[, levels(ART.P)])])
			set(YX, NULL, 'ART.A', YX[,factor(as.character(ART.A), levels=X.den[, levels(ART.A)])])
			predict.df		<- data.table( 	ART.pulse=factor('No',levels=c('No','Yes')), ART.I=factor('No',levels=c('No','Yes')), ART.F=factor('No',levels=c('No','Yes')), ART.P=factor('No',levels=c('No','Yes')), ART.A=factor('No',levels=c('No','Yes')),
					w=1.)
			risk.df			<- data.table( risk= c('ART.pulse','ART.I','ART.F','ART.P','ART.A'), factor=rep('Yes',5) )
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			tmp				<- data.table(risk.ref= rep('ART.pulse',nrow(risk.df)), factor.ref= rep('No',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))			
		}
		if(method%in%c('m3.ni','m3.ni.clu'))
		{
			#	ART indicators and number of drugs
			#	ART.nDrug.c independent of indicators			
			set(YX, NULL, 'stage', YX[, ART.nDrug.c])
			set(X.den, NULL, 'stage', X.den[, ART.nDrug.c])
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			set(YX, NULL, 'ART.pulse', YX[,factor(as.character(ART.pulse), levels=X.den[, levels(ART.pulse)])])
			set(YX, NULL, 'ART.I', YX[,factor(as.character(ART.I), levels=X.den[, levels(ART.I)])])
			set(YX, NULL, 'ART.F', YX[,factor(as.character(ART.F), levels=X.den[, levels(ART.F)])])
			set(YX, NULL, 'ART.P', YX[,factor(as.character(ART.P), levels=X.den[, levels(ART.P)])])
			set(YX, NULL, 'ART.A', YX[,factor(as.character(ART.A), levels=X.den[, levels(ART.A)])])			
			formula			<- 'score.Y ~ stage+ART.pulse+ART.I+ART.F+ART.P+ART.A-1'
			include.colnames<- c('score.Y','w','stage','ART.pulse','ART.I','ART.F','ART.P','ART.A')
			predict.df		<- data.table(	stage=factor('ART.3', levels=X.den[, levels(stage)]), 
					ART.pulse=factor('No',levels=c('No','Yes')), ART.I=factor('No',levels=c('No','Yes')), ART.F=factor('No',levels=c('No','Yes')), ART.P=factor('No',levels=c('No','Yes')), ART.A=factor('No',levels=c('No','Yes')),
					w=1.)
			risk.df			<- data.table( risk= c(rep('stage',3), 'ART.pulse','ART.I','ART.F','ART.P','ART.A'), factor=c('ART.3','ART.l3','ART.g3', rep('Yes',5)))
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )			
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))				
		}
		if(method%in%c('m3.nic','m3.nic.clu'))
		{
			#	number of drugs conditional on no indicators and ART indicators (not overlapping)  
			set(YX, NULL, 'stage', YX[, ART.nstage.c])
			set(X.den, NULL, 'stage', X.den[, ART.nstage.c])
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			predict.df		<- data.table(stage=factor('ART.3', levels=X.den[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)])
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))				
		}
		if(method%in%c('m3.nic.adj','m3.nic.clu.adj'))
		{
			#	number of drugs conditional on no indicators and ART indicators (not overlapping)  
			set(YX, NULL, 'stage', YX[, ART.nstage.c])
			set(X.den, NULL, 'stage', X.den[, ART.nstage.c])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(ART.nstage.c))])						
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#	adjust weight by selection bias weights			
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, ART.nstage.c])
				set(X.clu, NULL, 'stage', X.clu[, factor(as.character(stage), levels=X.msm[, levels(stage)])])				
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.clu[, table(stage)] / nrow(X.clu) )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.den[, table(stage)] / nrow(X.den) )
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			predict.df		<- data.table(stage=factor('ART.3', levels=YX[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)])
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))	
		}
		if(method%in%c('m3.nicv','m3.nicv.clu'))
		{
			#	number of drugs conditional on no indicators and ART indicators (not overlapping)  
			set(YX, NULL, 'stage', YX[, ART.nstage.c])
			set(X.den, NULL, 'stage', X.den[, ART.nstage.c])
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
			include.colnames<- c('score.Y','w','lRNA.mx','stage')
			predict.df		<- data.table(	stage=factor('ART.3', levels=X.den[, levels(stage)]), 
					lRNA.mx=subset(YX, stage=='ART.3')[, mean(lRNA.mx, na.rm=TRUE)],
					w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)])
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			risk.df			<- merge(risk.df, risk.df[, {
								tmp				<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
								list(lRNA.mx=ifelse(is.nan(tmp), YX[, mean( lRNA.mx, na.rm=TRUE )], tmp))
							}, by='coef'], by='coef')			
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))
		}
		if(method%in%c('m3.nicv.adj','m3.nicv.clu.adj'))
		{
			#	number of drugs conditional on no indicators and ART indicators (not overlapping)  
			set(YX, NULL, 'stage', YX[, ART.nstage.c])
			set(X.den, NULL, 'stage', X.den[, ART.nstage.c])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(ART.nstage.c))])
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#	adjust weight by selection bias weights			
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, ART.nstage.c])
				set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])				
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.clu[, table(stage)] / nrow(X.clu) )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.den[, table(stage)] / nrow(X.den) )
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#
			formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
			include.colnames<- c('score.Y','w','lRNA.mx','stage')
			predict.df		<- data.table(	stage=factor('ART.3', levels=X.den[, levels(stage)]), 
					lRNA.mx=subset(YX, stage=='ART.3')[, mean(lRNA.mx, na.rm=TRUE)],
					w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)])
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			risk.df			<- merge(risk.df, risk.df[, {
								tmp				<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
								list(lRNA.mx=ifelse(is.nan(tmp), YX[, mean( lRNA.mx, na.rm=TRUE )], tmp))
							}, by='coef'], by='coef')			
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))	
		}
		if(method%in%c('m3.tni','m3.tni.clu'))
		{
			#	ART indicators and number/type of drugs
			#	ART.tnDrug.c independent of indicators
			set(YX, NULL, 'stage', YX[, ART.tnDrug.c])
			set(X.den, NULL, 'stage', X.den[, ART.tnDrug.c])
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			set(YX, NULL, 'ART.pulse', YX[,factor(as.character(ART.pulse), levels=X.den[, levels(ART.pulse)])])
			set(YX, NULL, 'ART.I', YX[,factor(as.character(ART.I), levels=X.den[, levels(ART.I)])])
			set(YX, NULL, 'ART.F', YX[,factor(as.character(ART.F), levels=X.den[, levels(ART.F)])])
			set(YX, NULL, 'ART.P', YX[,factor(as.character(ART.P), levels=X.den[, levels(ART.P)])])
			set(YX, NULL, 'ART.A', YX[,factor(as.character(ART.A), levels=X.den[, levels(ART.A)])])						
			formula			<- 'score.Y ~ stage+ART.pulse+ART.I+ART.F+ART.P+ART.A-1'
			include.colnames<- c('score.Y','w','stage','ART.pulse','ART.I','ART.F','ART.P','ART.A')
			predict.df		<- data.table( 	stage=factor('ART.3.NRT.PI', levels=X.den[, levels(stage)]), 
					ART.pulse=factor('No',levels=c('No','Yes')), ART.I=factor('No',levels=c('No','Yes')), ART.F=factor('No',levels=c('No','Yes')), ART.P=factor('No',levels=c('No','Yes')), ART.A=factor('No',levels=c('No','Yes')),
					w=1.)
			risk.df			<- data.table( risk= c(rep('stage',9), 'ART.pulse','ART.I','ART.F','ART.P','ART.A'), factor=c('ART.3.NRT.PI','ART.l3','ART.g3','ART.3.NRT','ART.3.PI','ART.3.NNRT','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NNRT.PI', rep('Yes',5)))
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))
		}
		if(method%in%c('m3.tniv','m3.tniv.clu'))
		{
			#	ART indicators and number/type of drugs
			#	ART.tnDrug.c independent of indicators
			set(YX, NULL, 'stage', YX[, ART.tnDrug.c])
			set(X.den, NULL, 'stage', X.den[, ART.tnDrug.c])	
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			set(YX, NULL, 'ART.pulse', YX[,factor(as.character(ART.pulse), levels=X.den[, levels(ART.pulse)])])
			set(YX, NULL, 'ART.I', YX[,factor(as.character(ART.I), levels=X.den[, levels(ART.I)])])
			set(YX, NULL, 'ART.F', YX[,factor(as.character(ART.F), levels=X.den[, levels(ART.F)])])
			set(YX, NULL, 'ART.P', YX[,factor(as.character(ART.P), levels=X.den[, levels(ART.P)])])
			set(YX, NULL, 'ART.A', YX[,factor(as.character(ART.A), levels=X.den[, levels(ART.A)])])						
			formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage+ART.pulse+ART.I+ART.F+ART.P+ART.A-1'
			include.colnames<- c('score.Y','w','stage','lRNA.mx','ART.pulse','ART.I','ART.F','ART.P','ART.A')
			predict.df		<- data.table( 	stage=factor('ART.3.NRT.PI', levels=X.den[, levels(stage)]), 
					lRNA.mx=subset(YX, stage=='ART.3.NRT.PI')[, mean(lRNA.mx, na.rm=TRUE)],
					ART.pulse=factor('No',levels=c('No','Yes')), ART.I=factor('No',levels=c('No','Yes')), ART.F=factor('No',levels=c('No','Yes')), ART.P=factor('No',levels=c('No','Yes')), ART.A=factor('No',levels=c('No','Yes')),
					w=1.)
			risk.df			<- data.table( risk= c(rep('stage',9), 'ART.pulse','ART.I','ART.F','ART.P','ART.A'), factor=c('ART.3.NRT.PI','ART.l3','ART.g3','ART.3.NRT','ART.3.PI','ART.3.NNRT','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NNRT.PI', rep('Yes',5)) )
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]			
			risk.df			<- merge(risk.df, risk.df[, {
								tmp				<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
								list(lRNA.mx=ifelse(is.nan(tmp), YX[, mean( lRNA.mx, na.rm=TRUE )], tmp))
							}, by='coef'], by='coef')			
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))
		}
		if(method%in%c('m3.tnic','m3.tnic.clu'))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)  
			set(YX, NULL, 'stage', YX[, ART.ntstage.c])
			set(X.den, NULL, 'stage', X.den[, ART.ntstage.c])
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			predict.df		<- data.table(stage=factor('ART.3.NRT.PI', levels=X.den[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)])
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))
		}
		if(method%in%c('m3.tnic.adj','m3.tnic.clu.adj'))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			set(YX, NULL, 'stage', YX[, ART.ntstage.c])
			set(X.den, NULL, 'stage', X.den[, ART.ntstage.c])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(ART.ntstage.c))])						
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#	adjust weight by selection bias weights			
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, ART.ntstage.c])
				set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])				
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.clu[, table(stage)] / nrow(X.clu) )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / nrow(X.msm) ) / ( X.den[, table(stage)] / nrow(X.den) )
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#
			formula			<- 'score.Y ~ stage-1'
			include.colnames<- c('score.Y','w','stage')
			predict.df		<- data.table(stage=factor('ART.3.NRT.PI', levels=YX[, levels(stage)]), w=1.)			
			risk.df			<- data.table(risk='stage',factor=YX[, levels(stage)])
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))				
		}
		if(method%in%c('m3.tnicv','m3.tnicv.clu'))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			set(YX, NULL, 'stage', YX[, ART.ntstage.c])
			set(X.den, NULL, 'stage', X.den[, ART.ntstage.c])			
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
			include.colnames<- c('score.Y','w','stage','lRNA.mx')
			predict.df		<- data.table( 	stage=factor('ART.3.NRT.PI', levels=X.den[, levels(stage)]), 
					lRNA.mx=subset(YX, stage=='ART.3.NRT.PI')[, mean(lRNA.mx, na.rm=TRUE)], w=1.	)
			risk.df			<- data.table( risk= rep('stage',14), factor=c('ART.3.NRT.PI','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A', 'ART.F','ART.l3','ART.g3','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NRT','ART.3.NNRT.PI','ART.3.PI','ART.3.NNRT'))
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]			
			risk.df			<- merge(risk.df, risk.df[, {
								tmp				<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
								list(lRNA.mx=ifelse(is.nan(tmp), YX[, mean( lRNA.mx, na.rm=TRUE )], tmp))
							}, by='coef'], by='coef')			
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))
		}
		if(method%in%c('m3.tnicv.adj','m3.tnicv.clu.adj'))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			set(YX, NULL, 'stage', YX[, ART.ntstage.c])
			set(X.den, NULL, 'stage', X.den[, ART.ntstage.c])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(ART.ntstage.c))])						
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			#	adjust weight by selection bias weights			
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, ART.ntstage.c])
				set(X.clu, NULL, 'stage', X.clu[, factor(as.character(stage), levels=X.msm[, levels(stage)])])				
				bias.adj	<- ( X.msm[, table(stage)] / X.msm[, length(which(!is.na(stage)))] ) / ( X.clu[, table(stage)] / X.clu[, length(which(!is.na(stage)))] )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / X.msm[, length(which(!is.na(stage)))] ) / ( X.den[, table(stage)] / X.den[, length(which(!is.na(stage)))] )			
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)			
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#	
			formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
			include.colnames<- c('score.Y','w','stage','lRNA.mx')
			predict.df		<- data.table( 	stage=factor('ART.3.NRT.PI', levels=X.den[, levels(stage)]), 
					lRNA.mx=subset(YX, stage=='ART.3.NRT.PI')[, mean(lRNA.mx, na.rm=TRUE)], w=1.	)
			risk.df			<- data.table( risk= rep('stage',14), factor=c('ART.3.NRT.PI','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A', 'ART.F','ART.l3','ART.g3','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NRT','ART.3.NNRT.PI','ART.3.PI','ART.3.NNRT'))
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]			
			risk.df			<- merge(risk.df, risk.df[, {
								tmp				<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
								list(lRNA.mx=ifelse(is.nan(tmp), YX[, mean( lRNA.mx, na.rm=TRUE )], tmp))
							}, by='coef'], by='coef')			
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))		
		}
		if(method%in%c('m3.tnicvNo','m3.tnicvNo.clu'))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			set(YX, NULL, 'stage', YX[, ART.ntstage.no.c])
			set(X.den, NULL, 'stage', X.den[, ART.ntstage.no.c])
			YX				<- subset(YX, !is.na(stage) & !is.na(lRNA.mx))
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.den[, levels(stage)])])
			formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
			include.colnames<- c('score.Y','w','stage','lRNA.mx')
			predict.df		<- data.table( 	stage=factor('ART.3.NRT.PI', levels=X.den[, levels(stage)]), 
					lRNA.mx=subset(YX, stage=='ART.3.NRT.PI')[, mean(lRNA.mx, na.rm=TRUE)], w=1.	)
			risk.df			<- data.table( risk= rep('stage',14), factor=c('ART.3.NRT.PI','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A', 'ART.F','ART.l3','ART.g3','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NRT','ART.3.NNRT.PI','ART.3.PI','ART.3.NNRT'))
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]			
			risk.df			<- merge(risk.df, risk.df[, {
								tmp				<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
								list(lRNA.mx=ifelse(is.nan(tmp), YX[, mean( lRNA.mx, na.rm=TRUE )], tmp))
							}, by='coef'], by='coef')			
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk']))
		}
		if(method%in%c('m3.tnicvNo.adj','m3.tnicvNo.clu.adj'))
		{
			#	number/type of drugs conditional on no indicators and ART indicators (not overlapping)
			set(YX, NULL, 'stage', YX[, ART.ntstage.no.c])
			set(X.den, NULL, 'stage', X.den[, ART.ntstage.no.c])
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(ART.ntstage.no.c))])						
			set(YX, NULL, 'stage', YX[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.den, NULL, 'stage', X.den[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
			YX				<- subset(YX, !is.na(stage) & !is.na(lRNA.mx))
			#	adjust weight by selection bias weights			
			if(grepl('clu',method))
			{
				set(X.clu, NULL, 'stage', X.clu[, ART.ntstage.no.c])
				set(X.clu, NULL, 'stage', X.clu[,factor(as.character(stage), levels=X.msm[, levels(stage)])])
				bias.adj	<- ( X.msm[, table(stage)] / X.msm[, length(which(!is.na(stage)))] ) / ( X.clu[, table(stage)] / X.clu[, length(which(!is.na(stage)))] )
			}
			else
				bias.adj	<- ( X.msm[, table(stage)] / X.msm[, length(which(!is.na(stage)))] ) / ( X.den[, table(stage)] / X.den[, length(which(!is.na(stage)))] )
			bias.adj		<- cbind( data.table(stage=factor(rownames(bias.adj))), data.table(w.b=unclass(bias.adj)) )
			set(bias.adj, bias.adj[, which(is.infinite(w.b))], 'w.b', 1.)			
			YX				<- merge( YX, bias.adj, by='stage' )
			set(bias.adj, NULL, 'w.b', bias.adj[,w.b] * YX[, sum(w)] / YX[, sum(w*w.b)] )
			set(YX, NULL, 'w', YX[, w*w.b*sum(w)/sum(w*w.b) ] )
			#	
			formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
			include.colnames<- c('score.Y','w','stage','lRNA.mx')
			predict.df		<- data.table( 	stage=factor('ART.3.NRT.PI', levels=X.den[, levels(stage)]), 
					lRNA.mx=subset(YX, stage=='ART.3.NRT.PI')[, mean(lRNA.mx, na.rm=TRUE)], w=1.	)
			risk.df			<- data.table( risk= rep('stage',14), factor=c('ART.3.NRT.PI','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A', 'ART.F','ART.l3','ART.g3','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NRT','ART.3.NNRT.PI','ART.3.PI','ART.3.NNRT'))
			risk.df[, coef:=paste(risk.df[,risk],risk.df[,factor],sep='')]			
			risk.df			<- merge(risk.df, risk.df[, {
								tmp				<- YX[ which(unclass(YX[, risk, with=FALSE])[[1]]==factor), mean( lRNA.mx, na.rm=TRUE )]
								list(lRNA.mx=ifelse(is.nan(tmp), YX[, mean( lRNA.mx, na.rm=TRUE )], tmp))
							}, by='coef'], by='coef')			
			tmp				<- data.table(risk.ref= rep('stage',nrow(risk.df)), factor.ref= rep('ART.3.NRT.PI',nrow(risk.df)))
			risk.df			<- cbind(risk.df, tmp[, list(coef.ref=paste(risk.ref,factor.ref,sep='') ), by=c('risk.ref','factor.ref')])			
			ans				<- project.athena.Fisheretal.estimate.risk.core(YX, X.den, formula, predict.df, risk.df, include.colnames, bs.n=bs.n )
			ans$bias.adj	<- bias.adj
			ans$risk.table	<- do.call('rbind',list(
							risk.df[,	{
										z	<- table( YX[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.den[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X')												
									},by='risk'],
							risk.df[,	{
										z	<- table( X.msm[, risk, with=FALSE], useNA='ifany')
										list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')												
									},by='risk']))
		}		
		if(!is.na(save.file))
		{
			ans$YX			<- YX
			cat(paste('\nsave to file', save.file))
			save(ans, file=save.file)
		}
	}
	ans
}
######################################################################################
project.athena.Fisheretal.estimate.risk.table<- function(YX=NULL, X.den=NULL, X.msm=NULL, X.clu=NULL, resume=TRUE, save.file=NA, method=NA)
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
		if(substr(method,1,2)=='m2' || (substr(method,1,2)=='m3' && grepl('nic',method)))
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
				risktp.col		<- 'CD4t.tperiod'
				risk.col		<- 'CD4t'
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
				risktp.col		<- 'CD4t.tperiod'
				risk.col		<- 'CD4t'
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
				risktp.col		<- 'CD4t.tperiod'
				risk.col		<- 'CD4t'
			}			
			if(grepl('m3',method) & grepl('nic',method) & !grepl('tnic',method))
			{
				factor.ref.v	<- paste('ART.3',tp,sep='')
				risktp.col		<- 'ART.nstage.c.tperiod'
				risk.col		<- 'ART.nstage.c'
				YX[, ART.nstage.c.tperiod:= YX[, paste(ART.nstage.c, t.period, sep='.')]]
				set(YX, YX[,which(is.na(ART.nstage.c))], risktp.col, NA_character_)
				X.den[, ART.nstage.c.tperiod:= X.den[, paste(ART.nstage.c, t.period, sep='.')]]
				set(X.den, X.den[,which(is.na(ART.nstage.c))], risktp.col, NA_character_)
				X.clu[, ART.nstage.c.tperiod:= X.clu[, paste(ART.nstage.c, t.period, sep='.')]]
				set(X.clu, X.clu[,which(is.na(ART.nstage.c))], risktp.col, NA_character_)
				X.msm[, ART.nstage.c.tperiod:= X.msm[, paste(ART.nstage.c, t.period, sep='.')]]
				set(X.msm, X.msm[,which(is.na(ART.nstage.c))], risktp.col, NA_character_)
			}
			if(grepl('m3',method) & grepl('tnic',method) & !grepl('No',method))
			{
				factor.ref.v	<- paste('ART.3.NRT.PI',tp,sep='')
				risktp.col		<- 'ART.ntstage.c.tperiod'
				risk.col		<- 'ART.ntstage.c'
				YX[, ART.ntstage.c.tperiod:= YX[, paste(ART.ntstage.c, t.period, sep='.')]]
				set(YX, YX[,which(is.na(ART.ntstage.c))], risktp.col, NA_character_)
				X.den[, ART.ntstage.c.tperiod:= X.den[, paste(ART.ntstage.c, t.period, sep='.')]]
				set(X.den, X.den[,which(is.na(ART.ntstage.c))], risktp.col, NA_character_)
				X.clu[, ART.ntstage.c.tperiod:= X.clu[, paste(ART.ntstage.c, t.period, sep='.')]]
				set(X.clu, X.clu[,which(is.na(ART.ntstage.c))], risktp.col, NA_character_)
				X.msm[, ART.ntstage.c.tperiod:= X.msm[, paste(ART.ntstage.c, t.period, sep='.')]]
				set(X.msm, X.msm[,which(is.na(ART.ntstage.c))], risktp.col, NA_character_)
			}
			if(grepl('m3',method) & grepl('tnic',method) & grepl('No',method))
			{
				factor.ref.v	<- paste('ART.3.NRT.PI',tp,sep='')
				risktp.col		<- 'ART.ntstage.no.c.tperiod'
				risk.col		<- 'ART.ntstage.no.c'
				YX[, stage:=NA_character_]
				YX[, ART.ntstage.no.c.tperiod:= YX[, paste(ART.ntstage.no.c, t.period, sep='.')]]
				YX				<- subset(YX, !is.na(ART.ntstage.no.c))
				X.den[, ART.ntstage.no.c.tperiod:= X.den[, paste(ART.ntstage.no.c, t.period, sep='.')]]
				X.den			<- subset(X.den, !is.na(ART.ntstage.no.c))
				X.clu[, ART.ntstage.no.c.tperiod:= X.clu[, paste(ART.ntstage.no.c, t.period, sep='.')]]
				X.clu			<- subset(X.clu, !is.na(ART.ntstage.no.c))
				X.msm[, ART.ntstage.no.c.tperiod:= X.msm[, paste(ART.ntstage.no.c, t.period, sep='.')]]
				X.msm			<- subset(X.msm, !is.na(ART.ntstage.no.c))
			}			
			gc()
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
			#	adjust for censoring, keeping number of undiagnosed as in tperiod=='2'
			cens.table[, n.adjbyNU:=n]
			for(f in c('U','UAy','UAm'))
				for(z in cens.table[, unique(t.period)])	
					set(cens.table, cens.table[, which(factor2==f & stat=='X.msm' & t.period==z)], 'n.adjbyNU',  cens.table[which(factor2==f & stat=='X.msm' & t.period%in%c('2',z)), max(n.adjbyNU)])
			#	adjust for censoring, keeping proportion of undiagnosed as in tperiod=='2'
			cens.table[, n.adjbyPU:=n]
			for(z in cens.table[, unique(t.period)])
			{
				tmp		<- rbind(	unadjusted	= sapply(c('U','UAy','UAm'), function(f2)	subset(cens.table, factor2==f2 & stat=='X.msm' & t.period==z)[, p]),
									adjusted	= sapply(c('U','UAy','UAm'), function(f2)	subset(cens.table, factor2==f2 & stat=='X.msm' & t.period%in%c('2',z))[, max(p, na.rm=TRUE) ])		)		
				tmp		<- cens.table[which(stat=='X.msm' & t.period==z), sum[1] * tmp['adjusted',] * (1-sum(tmp['unadjusted',]))/(1-sum(tmp['adjusted',])) ]
				for(f2 in c('U','UAy','UAm'))
					set(cens.table, cens.table[, which(factor2==f2 & stat=='X.msm' & t.period==z)], 'n.adjbyPU',  tmp[f2])				
			}
			cens.table	<- merge(cens.table, cens.table[, list(factor=factor, p.adjbyNU= n.adjbyNU/sum(n.adjbyNU, na.rm=TRUE), p.adjbyPU= n.adjbyPU/sum(n.adjbyPU, na.rm=TRUE)), by=c('stat','t.period')], by=c('stat','t.period','factor'))
			#	if 'tp' is not '', reduce data to t.period
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
			ans$cens.table	<- cens.table
			#	subsetting or resetting: need to readjust levels
			set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])									
			set(X.den, NULL, 'stage', X.den[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(X.clu, NULL, 'stage', X.clu[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
			set(YX,    NULL, 'stage', YX[,    factor(as.character(stage), levels=X.msm[, levels(stage)])])
			risk.df			<- data.table(risk='stage',factor=X.den[, levels(stage)], risk.ref='stage', factor.ref=factor.ref.v)
			#	risk tables
			risk.table		<- do.call('rbind',list(
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
			risk.table		<- merge(risk.table, risk.table[, list(factor=factor, p= n/sum(n, na.rm=TRUE)), by='stat'], by=c('stat','factor'))
			ans$risk.table	<- risk.table
			#	bias + clustering adjustment
			ans$adj.clu			<- risk.table[, list(w.b= p[stat=='X.msm']/p[stat=='X.clu']),by=c('risk','factor')]
			#	bias adjustment
			ans$adj.seq			<- risk.table[, list(w.b= p[stat=='X.msm']/p[stat=='X.seq']),by=c('risk','factor')]			
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
	if('score.Y'%in%colnames(YX.m3))
		set(YX.m3, NULL, 'score.Y', YX.m3[,(score.Y*(nrow(YX.m3)-1)+0.5)/nrow(YX.m3)] )					
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
	}
	YX.m3[, stage.orig:= stage]
	if(return.only.ART)
		YX.m3	<- subset( YX.m3,stage=='ART.started' )
	#
	#	number and type of drugs, comparing those with at least one indicator to those with confirmed No indicator
	#	do this before we collapse the missing ART indicators into 'No'
	#
	YX.m3[, ART.ntstage.no.c:=NULL]
	YX.m3[, ART.ntstage.no.c:=NA_integer_]
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='Yes')], 'ART.ntstage.no.c', 2L )			#ART.I excludes pulsed
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'ART.ntstage.no.c', 3L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'ART.ntstage.no.c', 4L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'ART.ntstage.no.c', 5L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'ART.ntstage.no.c', 6L)	
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' &  ART.I=='No' & ART.F=='No' & ART.P=='No' & ART.A=='No' & is.na(ART.ntstage.no.c) & ART.nDrug<3 & ART.nDrug>0)], 'ART.ntstage.no.c', 7L)		
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' &  ART.I=='No' & ART.F=='No' & ART.P=='No' & ART.A=='No' & is.na(ART.ntstage.no.c) & ART.nDrug>3)], 'ART.ntstage.no.c', 8L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' &  ART.I=='No' & ART.F=='No' & ART.P=='No' & ART.A=='No' & is.na(ART.ntstage.no.c) & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT==0)], 'ART.ntstage.no.c', 1L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' &  ART.I=='No' & ART.F=='No' & ART.P=='No' & ART.A=='No' & is.na(ART.ntstage.no.c) & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT>0)], 'ART.ntstage.no.c', 9L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' &  ART.I=='No' & ART.F=='No' & ART.P=='No' & ART.A=='No' & is.na(ART.ntstage.no.c) & ART.nDrug==3 & ART.nNRT>0 & ART.nPI==0 & ART.nNNRT>0 )], 'ART.ntstage.no.c', 10L)	
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' &  ART.I=='No' & ART.F=='No' & ART.P=='No' & ART.A=='No' & is.na(ART.ntstage.no.c) & ART.nNRT>0 & ART.nPI==0 & ART.nNNRT==0)], 'ART.ntstage.no.c', 11L)				
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' &  ART.I=='No' & ART.F=='No' & ART.P=='No' & ART.A=='No' & ART.nDrug==3 & ART.nNRT==0 & ART.nNNRT>0 & ART.nPI>0)], 'ART.ntstage.no.c', 12L)	#-- these are all on pulsed, need to take out	
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' &  ART.I=='No' & ART.F=='No' & ART.P=='No' & ART.A=='No' & is.na(ART.ntstage.no.c) & ART.nNRT==0 & ART.nPI>0 & ART.nNNRT==0)], 'ART.ntstage.no.c', 13L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' &  ART.I=='No' & ART.F=='No' & ART.P=='No' & ART.A=='No' & is.na(ART.ntstage.no.c) & ART.nNRT==0 & ART.nPI==0 & ART.nNNRT>0)], 'ART.ntstage.no.c', 14L)
	if(!return.only.ART)
	{
		set(YX.m3, YX.m3[, which(stage=='U')], 'ART.ntstage.no.c', 15L)
		set(YX.m3, YX.m3[, which(stage=='UAy')], 'ART.ntstage.no.c', 16L)
		set(YX.m3, YX.m3[, which(stage=='UAm')], 'ART.ntstage.no.c', 17L)
		set(YX.m3, YX.m3[, which(stage=='Diag')], 'ART.ntstage.no.c', 18L)
		set(YX.m3, NULL, 'ART.ntstage.no.c', YX.m3[, factor(ART.ntstage.no.c, levels=1:18, labels=c('ART.3.NRT.PI','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A', 'ART.F','ART.l3','ART.g3','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NRT','ART.3.NNRT.PI','ART.3.PI','ART.3.NNRT','U','UAy','UAm','Diag'))])
	}
	else
		set(YX.m3, NULL, 'ART.ntstage.no.c', YX.m3[, factor(ART.ntstage.no.c, levels=1:14, labels=c('ART.3.NRT.PI','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A', 'ART.F','ART.l3','ART.g3','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NRT','ART.3.NNRT.PI','ART.3.PI','ART.3.NNRT'))])
	#YX.m3[, table(ART.ntstage.no.c)]
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
	set(YX.m3, YX.m3[,which(is.na(ART.A) | stage!='ART.started' | ART.pulse=='Yes')], 'ART.A', 'No')
	set(YX.m3, NULL, 'ART.A', YX.m3[, factor(as.character(ART.A))])
	set(YX.m3, YX.m3[,which(is.na(ART.P) | stage!='ART.started' | ART.pulse=='Yes')], 'ART.P', 'No')
	set(YX.m3, NULL, 'ART.P', YX.m3[, factor(as.character(ART.P))])		
	#	number of drugs -- independent. I can set nDrug==0 to NA but then I get 'contrasts can be applied only to factors with 2' because this is perf corr with ART.I
	#	set reference group for dummy coding to ART.3  
	YX.m3[, ART.nDrug.c:=NULL]
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
	else
		set(YX.m3, NULL, 'ART.nDrug.c', YX.m3[, factor(ART.nDrug.c, levels=1:3, labels=c('ART.3','ART.l3','ART.g3'))])
	if( YX.m3[, any(is.na(ART.nDrug.c))] ) stop('unexpected NA in ART.nDrug.c')
	#	type of drugs for ART.3 
	#	set reference group for dummy coding to ART.3.NRT.PI
	#	include pulse as separate
	#YX.m3[, ART.tDrug.c:=NULL]
	#YX.m3[, ART.tDrug.c:=NA_integer_]
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==0)], 'ART.tDrug.c', 2L)
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT>0)], 'ART.tDrug.c', 3L)
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT==0)], 'ART.tDrug.c', 1L)
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nNRT>0 & ART.nPI==0 & ART.nNNRT>0 )], 'ART.tDrug.c', 4L)		
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nNRT>0 & ART.nPI==0 & ART.nNNRT==0)], 'ART.tDrug.c', 5L)
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nNRT==0 & ART.nNNRT>0 & ART.nPI>0)], 'ART.tDrug.c', 6L)
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nNRT==0 & ART.nPI==0 & ART.nNNRT>0)], 'ART.tDrug.c', 7L)
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nNRT==0 & ART.nPI>0 & ART.nNNRT==0)], 'ART.tDrug.c', 8L)
	#set(YX.m3, YX.m3[, which(ART.P=='Yes' | ART.A=='Yes' | ART.F=='Yes')], 'ART.tnDrug.c', NA_integer_)		#exclude ART with indication		
	#set(YX.m3, NULL, 'ART.tDrug.c', YX.m3[, factor(ART.tDrug.c, levels=1:8, labels=c('ART.NRT.PI','ART.None','ART.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.NRT','ART.PI.NNRT','ART.NNRT','ART.PI'))])										
	#if( YX.m3[, any(is.na(ART.tDrug.c))] ) stop('unexpected NA in ART.tDrug.c')
	#	number and type of drugs for ART.3 
	#	set reference group for dummy coding to ART.3.NRT.PI
	#	do not include ART.I as separate		
	YX.m3[, ART.tnDrug.c:=NULL]
	YX.m3[, ART.tnDrug.c:=NA_integer_]
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug<3)], 'ART.tnDrug.c', 3L)		
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug>3)], 'ART.tnDrug.c', 4L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT>0)], 'ART.tnDrug.c', 5L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT==0)], 'ART.tnDrug.c', 1L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI==0 & ART.nNNRT>0 )], 'ART.tnDrug.c', 6L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT==0 & ART.nNNRT>0 & ART.nPI>0)], 'ART.tnDrug.c', 7L)	
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT==0 & ART.nNNRT==0 & ART.nPI>0)], 'ART.tnDrug.c', 8L)	
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT==0 & ART.nNNRT>0 & ART.nPI==0)], 'ART.tnDrug.c', 9L)	
	#-- these are all on pulsed, need to take out
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nNRT>0 & ART.nPI==0 & ART.nNNRT==0)], 'ART.tnDrug.c', 2L)		
	#set(YX.m3, YX.m3[, which(ART.P=='Yes' | ART.A=='Yes' | ART.F=='Yes')], 'ART.tnDrug.c', NA_integer_)		#exclude ART with indication		
	if(!return.only.ART)
	{
		set(YX.m3, YX.m3[, which(stage=='U')], 'ART.tnDrug.c', 10L)
		set(YX.m3, YX.m3[, which(stage=='UAy')], 'ART.tnDrug.c', 11L)
		set(YX.m3, YX.m3[, which(stage=='UAm')], 'ART.tnDrug.c', 12L)
		set(YX.m3, YX.m3[, which(stage=='Diag')], 'ART.tnDrug.c', 13L)
		set(YX.m3, NULL, 'ART.tnDrug.c', YX.m3[, factor(ART.tnDrug.c, levels=1:13, labels=c('ART.3.NRT.PI','ART.3.NRT','ART.l3','ART.g3','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NNRT.PI','ART.3.PI','ART.3.NNRT','U','UAy','UAm','Diag'))])
	}
	else
		set(YX.m3, NULL, 'ART.tnDrug.c', YX.m3[, factor(ART.tnDrug.c, levels=1:9, labels=c('ART.3.NRT.PI','ART.3.NRT','ART.l3','ART.g3','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NNRT.PI','ART.3.PI','ART.3.NNRT'))])
	if( YX.m3[, any(is.na(ART.tnDrug.c))] ) stop('unexpected NA in ART.tnDrug.c')		
	#	number and type of drugs for ART.3 
	#	set reference group for dummy coding to ART.3.NRT.PI
	#	do not include pulse as separate
	#YX.m3[, ART.tpDrug.c:=NULL]
	#YX.m3[, ART.tpDrug.c:=NA_integer_]
	#set(YX.m3, YX.m3[,which(ART.I=='Yes')], 'ART.tpDrug.c', 2L)
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' & ART.nDrug<3 & ART.nDrug>0)], 'ART.tpDrug.c', 3L)		
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' & ART.nDrug>3)], 'ART.tpDrug.c', 4L)
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT>0)], 'ART.tpDrug.c', 5L)
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT==0)], 'ART.tpDrug.c', 1L)
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI==0 & ART.nNNRT>0 )], 'ART.tpDrug.c', 6L)
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' & ART.nNRT>0 & ART.nPI==0 & ART.nNNRT==0)], 'ART.tpDrug.c', 7L)	
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='No' & ART.nDrug==3 & ART.nNRT==0 & ART.nNNRT>0 & ART.nPI>0)], 'ART.tpDrug.c', 8L)	
	#-- these are all on pulsed, need to take out
	#set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='Yes')], 'ART.tpDrug.c', 9L)	
	#set(YX.m3, YX.m3[, which(ART.P=='Yes' | ART.A=='Yes' | ART.F=='Yes')], 'ART.tpDrug.c', NA_integer_)		#exclude ART with indication		
	#set(YX.m3, NULL, 'ART.tpDrug.c', YX.m3[, factor(ART.tpDrug.c, levels=1:9, labels=c('ART.3.NRT.PI','ART.I','ART.l3','ART.g3','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NRT','ART.3.NNRT.PI','ART.pulse'))])											#still has pulse has NA
	#if( YX.m3[, any(is.na(ART.tpDrug.c))] ) stop('unexpected NA in ART.tpDrug.c')
	#
	#	number of drugs conditional on no indication
	#
	YX.m3[, ART.nstage.c:=NULL]
	YX.m3[, ART.nstage.c:=NA_integer_]
	set(YX.m3, YX.m3[, which(ART.pulse=='Yes')], 'ART.nstage.c', 2L )			#ART.I excludes pulsed
	set(YX.m3, YX.m3[, which(ART.I=='Yes')], 'ART.nstage.c', 3L)
	set(YX.m3, YX.m3[, which(ART.P=='Yes')], 'ART.nstage.c', 4L)
	set(YX.m3, YX.m3[, which(ART.A=='Yes')], 'ART.nstage.c', 5L)
	set(YX.m3, YX.m3[, which(ART.F=='Yes')], 'ART.nstage.c', 6L)				
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.nstage.c) & ART.nDrug<3 & ART.nDrug>0)], 'ART.nstage.c', 7L)		
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.nstage.c) & ART.nDrug>3)], 'ART.nstage.c', 8L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.nstage.c) & ART.nDrug==3)], 'ART.nstage.c', 1L)
	if(!return.only.ART)
	{
		set(YX.m3, YX.m3[, which(stage=='U')], 'ART.nstage.c', 9L)
		set(YX.m3, YX.m3[, which(stage=='UAy')], 'ART.nstage.c', 10L)
		set(YX.m3, YX.m3[, which(stage=='UAm')], 'ART.nstage.c', 11L)
		set(YX.m3, YX.m3[, which(stage=='Diag')], 'ART.nstage.c', 12L)
		set(YX.m3, NULL, 'ART.nstage.c', YX.m3[, factor(ART.nstage.c, levels=1:12, labels=c('ART.3','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A', 'ART.F','ART.l3','ART.g3','U','UAy','UAm','Diag'))])		
	}
	else
		set(YX.m3, NULL, 'ART.nstage.c', YX.m3[, factor(ART.nstage.c, levels=1:8, labels=c('ART.3','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A', 'ART.F','ART.l3','ART.g3'))])
	if( YX.m3[, any(is.na(ART.nstage.c))] ) stop('unexpected ART.nstage.c')
	#
	#	number and type of drugs conditional on no indication
	#
	YX.m3[, ART.ntstage.c:=NULL]
	YX.m3[, ART.ntstage.c:=NA_integer_]
	set(YX.m3, YX.m3[, which(ART.pulse=='Yes')], 'ART.ntstage.c', 2L )			#ART.I excludes pulsed
	set(YX.m3, YX.m3[, which(ART.I=='Yes')], 'ART.ntstage.c', 3L)
	set(YX.m3, YX.m3[, which(ART.P=='Yes')], 'ART.ntstage.c', 4L)
	set(YX.m3, YX.m3[, which(ART.A=='Yes')], 'ART.ntstage.c', 5L)
	set(YX.m3, YX.m3[, which(ART.F=='Yes')], 'ART.ntstage.c', 6L)				
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.ntstage.c) & ART.nDrug<3 & ART.nDrug>0)], 'ART.ntstage.c', 7L)		
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.ntstage.c) & ART.nDrug>3)], 'ART.ntstage.c', 8L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.ntstage.c) & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT==0)], 'ART.ntstage.c', 1L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.ntstage.c) & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT>0)], 'ART.ntstage.c', 9L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.ntstage.c) & ART.nDrug==3 & ART.nNRT>0 & ART.nPI==0 & ART.nNNRT>0 )], 'ART.ntstage.c', 10L)	
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.ntstage.c) & ART.nNRT>0 & ART.nPI==0 & ART.nNNRT==0)], 'ART.ntstage.c', 11L)				
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT==0 & ART.nNNRT>0 & ART.nPI>0)], 'ART.ntstage.c', 12L)	#-- these are all on pulsed, need to take out	
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.ntstage.c) & ART.nNRT==0 & ART.nPI>0 & ART.nNNRT==0)], 'ART.ntstage.c', 13L)
	set(YX.m3, YX.m3[, which(stage=='ART.started' & is.na(ART.ntstage.c) & ART.nNRT==0 & ART.nPI==0 & ART.nNNRT>0)], 'ART.ntstage.c', 14L)	
	if(!return.only.ART)
	{
		set(YX.m3, YX.m3[, which(stage=='U')], 'ART.ntstage.c', 15L)
		set(YX.m3, YX.m3[, which(stage=='UAy')], 'ART.ntstage.c', 16L)
		set(YX.m3, YX.m3[, which(stage=='UAm')], 'ART.ntstage.c', 17L)
		set(YX.m3, YX.m3[, which(stage=='Diag')], 'ART.ntstage.c', 18L)
		set(YX.m3, NULL, 'ART.ntstage.c', YX.m3[, factor(ART.ntstage.c, levels=1:18, labels=c('ART.3.NRT.PI','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A', 'ART.F','ART.l3','ART.g3','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NRT','ART.3.NNRT.PI','ART.3.PI','ART.3.NNRT','U','UAy','UAm','Diag'))])
	}
	else
		set(YX.m3, NULL, 'ART.ntstage.c', YX.m3[, factor(ART.ntstage.c, levels=1:14, labels=c('ART.3.NRT.PI','ART.pulse.Y', 'ART.I', 'ART.P', 'ART.A', 'ART.F','ART.l3','ART.g3','ART.3.NRT.PI.NNRT','ART.3.NRT.NNRT','ART.3.NRT','ART.3.NNRT.PI','ART.3.PI','ART.3.NNRT'))])
	if( YX.m3[, any(is.na(ART.ntstage.c))] ) stop('unexpected ART.ntstage.c')
	#
	#	add mx viral load during infection window
	#
	YX.m3	<- merge(YX.m3, YX.m3[, {
						tmp<- which(!is.na(lRNA))
						list( lRNA.mx= ifelse(length(tmp), max(lRNA[tmp]), NA_real_) )
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))								
	#	focus after ART initiation
	#YX.m3	<- subset( YX.m3, !stage%in%c('UAm','UAy','DAm','DAy','D1<=350','D1<=550','D1>550','D1.NA','U'))
	set(YX.m3, NULL, 'stage', YX.m3[, factor(as.character(stage))])
	#
	YX.m3
}
######################################################################################
project.athena.Fisheretal.YX.model3.stratify.v1.ARTriskgroups<- function(YX, df.all, cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), plot.file.cascade=NA, score.Y.cut=1e-2)
{	
	#	Treatment cascade & ART risk groups
	#	prop and odds by traditional categories 	undiagnosed, 	undiagnosed&evidence acute, 	diagnosed by CD4 at diagnosis
	#	
	#	PQ: are there ART risk groups?
	#	PQ: is VL the causal pathway?		--> currently problematic with uniform infection window approach
	#	SQ: is ART pulse associated with higher transmissibility?
	YX.m3	<- copy(YX)
	YX.m3[, U.score:=NULL]
	if('score.Y'%in%colnames(YX.m3))
		set(YX.m3, NULL, 'score.Y', YX.m3[,(score.Y*(nrow(YX.m3)-1)+0.5)/nrow(YX.m3)] )
	#	transmitter acute
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
	YX.m3[, stage.orig:= stage]
	#
	#	ART indication risk factors
	#
	set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )		
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')		
	#	ART no indication: nDrug
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug<3)], 'stage', 'ART.l3')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug>3)], 'stage', 'ART.g3')
	#	split weight for PI and NNRT
	#if('score.Y'%in%colnames(YX.m3))
	#{
	#	tmp		<- YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT>0)]
	#	tmp2	<- YX.m3[tmp,]
	#	set(tmp2, NULL, 'w', tmp2[, w/2])
	#	set(tmp2, NULL, 'stage', 'ART.3.NRT.PI')
	#	YX.m3	<- rbind(YX.m3, tmp2)
	#	set(YX.m3, tmp, 'w', YX.m3[tmp, w/2])
	#	set(YX.m3, tmp, 'stage', 'ART.3.NRT.NNRT')		
	#}
	#			
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT>0)], 'stage', 'ART.3.NRT.PI.NNRT')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 )], 'stage', 'ART.3.NRT.PI')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nNNRT>0)], 'stage', 'ART.3.NRT.NNRT')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 )], 'stage', 'ART.3.NRT')
	set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT==0 & ART.nNNRT>0 & ART.nPI>0)], 'stage', 'ART.3.NNRT.PI')
	#	add mx viral load during infection window
	YX.m3	<- merge(YX.m3, YX.m3[, {
						tmp<- which(!is.na(lRNA))
						list( lRNA.mx= ifelse(length(tmp), max(lRNA[tmp]), NA_real_) )
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))								
	#	focus after ART initiation
	YX.m3	<- subset( YX.m3, !stage%in%c('UAm','UAy','DAm','DAy','D1<=350','D1<=550','D1>550','D1.NA','U'))
	set(YX.m3, NULL, 'stage', YX.m3[, factor(as.character(stage))])
	#
	if('score.Y'%in%colnames(YX.m3))
		YX.m3	<- subset(YX.m3, select=c(t, t.Patient, Patient, score.Y, stage, lRNA.mx, CDCC, lRNA, contact, fw.up.med, w, stage.orig  ))	
	if(!'score.Y'%in%colnames(YX.m2))
		YX.m3	<- subset(YX.m3, select=c(t, t.Patient, Patient, stage, lRNA.mx, CDCC, lRNA, contact, fw.up.med, stage.orig  ))		
	
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
project.athena.Fisheretal.YX.model3.v1.ARTriskgroups<- function(YX, clumsm.info, cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), plot.file.cascade=NA, score.Y.cut=1e-2)
{
	# cd4.cut= c(-1, 350, 550, 5000); cd4.label=c('D1<=350','D1<=550','D1>550'); plot.file.cascade=NA; score.Y.cut=1e-2
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
	#	deal with acutes
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
	YX.m3[, table(stage)]
	#	
	YX.m3[, stage.orig:= stage]
	#
	#	stratify ART by 	ART.pulse	ART.I	ART.F 	ART.A	ART.P	ART.ni
	#	full model
	#	
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
	tmp				<- data.table(stage= YX.m3[, levels(stage)], col=sapply( rainbow(YX.m3[, nlevels(stage)]), my.fade.col, alpha=0.5) )
	set(tmp, NULL, 'stage', tmp[, factor(stage)])
	YX.m3			<- merge(YX.m3, tmp, by='stage')
	YX.m3.fit4 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m3)				
	YX.m3.fit4.or	<- data.table( 	F= my.or.from.logit(YX.m3.fit4, 'stageART.F', 'stageART.ni', subset(YX.m3, stage=='ART.F')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
									I= my.or.from.logit(YX.m3.fit4, 'stageART.I', 'stageART.ni', subset(YX.m3, stage=='ART.I')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
									P= my.or.from.logit(YX.m3.fit4, 'stageART.P', 'stageART.ni', subset(YX.m3, stage=='ART.P')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
									A= my.or.from.logit(YX.m3.fit4, 'stageART.A', 'stageART.ni', subset(YX.m3, stage=='ART.A')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
									pulse= my.or.from.logit(YX.m3.fit4, 'stageART.pulse.Y', 'stageART.ni', subset(YX.m3, stage=='ART.pulse.Y')[, sum(w)], subset(YX.m3, stage=='ART.ni')[, sum(w)], 1.962),
									ni= my.or.from.logit(YX.m3.fit4, 'stageART.ni', 'stageU', subset(YX.m3, stage=='ART.ni')[, sum(w)], subset(YX.m3, stage=='U')[, sum(w)], 1.962) )	
	YX.m3.fit4.p	<- YX.m3[, list(n=round(sum(w),d=1),prop= round(sum(w)/YX.m3[,sum(w)],d=3) ),by='stage']	
	#
	#	stratify ART by 	ART.pulse	ART.I	ART.F 	ART.A	ART.P	ART.ni
	#	sub model after treatment
	#	
	YX.m3s			<- subset(YX.m3, !stage%in%c('D1<=350','D1<=550','D1>550','U','UAm','UAy'))
	set(YX.m3s, NULL, 'stage', YX.m3s[, factor(as.character(stage))])
	YX.m3s.fit4 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m3s)
	YX.m3s.fit4.or	<- data.table( 	F= my.or.from.logit(YX.m3s.fit4, 'stageART.F', 'stageART.ni', subset(YX.m3s, stage=='ART.F')[, sum(w)], subset(YX.m3s, stage=='ART.ni')[, sum(w)], 1.962),
									I= my.or.from.logit(YX.m3s.fit4, 'stageART.I', 'stageART.ni', subset(YX.m3s, stage=='ART.I')[, sum(w)], subset(YX.m3s, stage=='ART.ni')[, sum(w)], 1.962),
									P= my.or.from.logit(YX.m3s.fit4, 'stageART.P', 'stageART.ni', subset(YX.m3s, stage=='ART.P')[, sum(w)], subset(YX.m3s, stage=='ART.ni')[, sum(w)], 1.962),
									A= my.or.from.logit(YX.m3s.fit4, 'stageART.A', 'stageART.ni', subset(YX.m3s, stage=='ART.A')[, sum(w)], subset(YX.m3s, stage=='ART.ni')[, sum(w)], 1.962),
									pulse= my.or.from.logit(YX.m3s.fit4, 'stageART.pulse.Y', 'stageART.ni', subset(YX.m3s, stage=='ART.pulse.Y')[, sum(w)], subset(YX.m3s, stage=='ART.ni')[, sum(w)], 1.962))
	if(!is.na(plot.file.cascade))
	{
		pdf(file=plot.file.cascade, w=5, h=5)
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
		dev.off()
		#		
	}
	list(m3=YX.m3, m3.fit=YX.m3.fit4, m3.fit.art=YX.m3s.fit4, m3.or=YX.m3.fit4.or, m3.or.art=YX.m3s.fit4.or)
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
project.athena.Fisheretal.YX.model3.v2<- function(YX, clumsm.info, vl.suppressed=log10(1e3), acute.select=c('Yes','Maybe'), cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), score.Y.cut=1e-2)
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
	#	if transmitter acute, set undiagnosed to undiagnosed and acute
	set(YX.m3, YX.m3[, which(t.isAcute=='Yes')], 'stage', 'Ay')
	set(YX.m3, YX.m3[, which(t.isAcute=='Maybe')], 'stage', 'Am')
	YX.m3[, stage.orig:= stage]
	#	number of drugs
	if(1)
	{
		YX.m3s	<- subset( YX.m3, stage=='ART.started') 		
		YX.m3s[, c(any(is.na( ART.nNRT )), any(is.na( ART.nNNRT )), any(is.na( ART.nPI )))]		#no missing entries - great
		YX.m3s[, table(ART.nDrug)]
		
		#	ART indication risk factors
		set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])
		set(YX.m3, YX.m3[, which(ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )		
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')
		
		#	ART no indication: break up by nDrug<3 or nDrug>3
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug<3)], 'stage', 'ART.3l')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug>3)], 'stage', 'ART.3g')
		set(YX.m3, YX.m3[, which(stage=='ART.started')], 'stage', 'ART.ni')
		
		ggplot(YX.m3, aes(stage, score.Y)) + geom_boxplot() + geom_jitter(size=0.4)
		#	some indication that 3l has higher transmissibility
		
		YX.m3s	<- subset( YX.m3, !stage%in%c('Am','Ay','D1<=350','D1<=550','D1>550','U'))
		set(YX.m3s, NULL, 'stage', YX.m3s[, factor(as.character(stage))])
		YX.m3s[, table(stage)]
		#   stage
    	# ART.3g      ART.3l       ART.F       ART.I       ART.P ART.pulse.Y ART.started 
        # 171          51         703         209         191         164        2328  		
		YX.m3s.fit8 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m3s)				
		summary(YX.m3s.fit8)
		#	plot
		tmp			<- data.table(stage= YX.m3s[, levels(stage)], col=sapply( rainbow(YX.m3s[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m3s		<- merge(YX.m3s, tmp, by='stage')		
		plot( seq_len(nrow(YX.m3s)), YX.m3s[, score.Y], pch=18, col= YX.m3s[, col], cex=YX.m3s[, w^0.4])
		tmp				<- subset(YX.m3s, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m3s.fit8, tmp, type='response'))	
		start			<- c( YX.m3s.fit8$coef$mean + head( 2*sqrt(diag(vcov(YX.m3s.fit8))), length(YX.m3s.fit8$coef$mean)), YX.m3s.fit8$coef$precision)
		YX.m3s.fit8.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3s, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m3s.fit8$coef$mean - head( 2*sqrt(diag(vcov(YX.m3s.fit8))), length(YX.m3s.fit8$coef$mean)), YX.m3s.fit8$coef$precision)
		YX.m3s.fit8.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3s, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
				c( predict(YX.m3s.fit8.sup, tmp, type='response'), rev(predict(YX.m3s.fit8.slw, tmp, type='response'))), 
				border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m3s[, levels(stage)], fill=YX.m3s[, unique(col)])
		YX.m3s[, col:=NULL]
		#
		YX.m3s.fit8.or	<- data.table( 	F= my.or.from.logit(YX.m3s.fit8, 'stageART.F', 'stageART.ni', subset(YX.m3s, stage=='ART.F')[, sum(w)], subset(YX.m3s, stage=='ART.ni')[, sum(w)], 1.962),
										I= my.or.from.logit(YX.m3s.fit8, 'stageART.I', 'stageART.ni', subset(YX.m3s, stage=='ART.I')[, sum(w)], subset(YX.m3s, stage=='ART.ni')[, sum(w)], 1.962),
										P= my.or.from.logit(YX.m3s.fit8, 'stageART.P', 'stageART.ni', subset(YX.m3s, stage=='ART.P')[, sum(w)], subset(YX.m3s, stage=='ART.ni')[, sum(w)], 1.962),
										A= my.or.from.logit(YX.m3s.fit8, 'stageART.A', 'stageART.ni', subset(YX.m3s, stage=='ART.A')[, sum(w)], subset(YX.m3s, stage=='ART.ni')[, sum(w)], 1.962),
										g3= my.or.from.logit(YX.m3s.fit8, 'stageART.3g', 'stageART.ni', subset(YX.m3s, stage=='ART.3g')[, sum(w)], subset(YX.m3s, stage=='ART.ni')[, sum(w)], 1.962),
										l3= my.or.from.logit(YX.m3s.fit8, 'stageART.3l', 'stageART.ni', subset(YX.m3s, stage=='ART.3l')[, sum(w)], subset(YX.m3s, stage=='ART.ni')[, sum(w)], 1.962)	)
	}
	#	number drugs != 3, no drugs==3 and which type
	if(1)
	{							
		#	ART indication risk factors
		set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])
		set(YX.m3, YX.m3[, which(ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )		
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')		
		#	ART no indication: nDrug
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug!=3)], 'stage', 'ART.3n')
		#	ART triple therapy - type
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 & ART.nNNRT>0)], 'stage', 'ART.3.NRT.PI.NNRT')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nPI>0 )], 'stage', 'ART.3.NRT.PI')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 & ART.nNNRT>0)], 'stage', 'ART.3.NRT.NNRT')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.nDrug==3 & ART.nNRT>0 )], 'stage', 'ART.3.NRT')
		
		YX.m3s	<- subset( YX.m3, !stage%in%c('Am','Ay','D1<=350','D1<=550','D1>550','U'))
		set(YX.m3s, NULL, 'stage', YX.m3s[, factor(as.character(stage))])
		YX.m3s[, table(stage)]
		#   stage
		# ART.3.NRT    ART.3.NRT.NNRT      ART.3.NRT.PI 	ART.3.NRT.PI.NNRT     ART.3n             ART.F             ART.I		  ART.P       ART.pulse.Y 
        # 72              1423               798                35               	222               703               209 			191               164
		YX.m3s.fit9 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m3s)				
		summary(YX.m3s.fit9)
		#	plot
		tmp			<- data.table(stage= YX.m3s[, levels(stage)], col=sapply( rainbow(YX.m3s[, nlevels(stage)]), my.fade.col, alpha=0.5) )
		set(tmp, NULL, 'stage', tmp[, factor(stage)])
		YX.m3s		<- merge(YX.m3s, tmp, by='stage')		
		plot( seq_len(nrow(YX.m3s)), YX.m3s[, score.Y], pch=18, col= YX.m3s[, col], cex=YX.m3s[, w^0.4])
		tmp				<- subset(YX.m3s, select=stage)
		lines(seq_len(nrow(tmp)), predict(YX.m3s.fit9, tmp, type='response'))	
		start			<- c( YX.m3s.fit9$coef$mean + head( 2*sqrt(diag(vcov(YX.m3s.fit9))), length(YX.m3s.fit9$coef$mean)), YX.m3s.fit9$coef$precision)
		YX.m3s.fit9.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3s, start=start, hessian=TRUE, maxit=0)		
		start			<- c( YX.m3s.fit9$coef$mean - head( 2*sqrt(diag(vcov(YX.m3s.fit9))), length(YX.m3s.fit9$coef$mean)), YX.m3s.fit9$coef$precision)
		YX.m3s.fit9.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX.m3s, start=start, hessian=TRUE, maxit=0)		
		polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
				c( predict(YX.m3s.fit9.sup, tmp, type='response'), rev(predict(YX.m3s.fit9.slw, tmp, type='response'))), 
				border=NA, col=my.fade.col('black',0.5) )
		legend('bottomright', bty='n', border=NA, legend= YX.m3s[, levels(stage)], fill=YX.m3s[, unique(col)])
		YX.m3s[, col:=NULL]
		#
		YX.m3s.fit9.or	<- data.table( 	F= my.or.from.logit(YX.m3s.fit9, 'stageART.F', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.F')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										I= my.or.from.logit(YX.m3s.fit9, 'stageART.I', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.I')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										P= my.or.from.logit(YX.m3s.fit9, 'stageART.P', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.P')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										Pulse= my.or.from.logit(YX.m3s.fit9, 'stageART.pulse.Y', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.pulse.Y')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										A= my.or.from.logit(YX.m3s.fit9, 'stageART.A', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.A')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										n3= my.or.from.logit(YX.m3s.fit9, 'stageART.3n', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.3n')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										NRT= my.or.from.logit(YX.m3s.fit9, 'stageART.3.NRT', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.3.NRT')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										NNRT= my.or.from.logit(YX.m3s.fit9, 'stageART.3.NRT.NNRT', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.3.NRT.NNRT')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										NNRTPI= my.or.from.logit(YX.m3s.fit9, 'stageART.3.NRT.PI.NNRT', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.3.NRT.PI.NNRT')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962))
		#           F        I         P    Pulse  A       n3      NRT      NNRT    NNRTPI
		#1: 1.6683804 2.515669  4.626841 3.291011 NA 3.060213 3.595480 1.4408900 1.3585172
		#2: 0.8347612 1.208841  2.095001 1.427700 NA 1.425840 1.553389 0.8554175 0.6189499
		#3: 3.3344786 5.235258 10.218448 7.586156 NA 6.567993 8.322113 2.4270770 2.9817746
	}
	#	number drugs!=3, NRT alone, split weight for PI and NNRT
	if(1)
	{							
		#	ART indication risk factors
		set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])
		set(YX.m3, YX.m3[, which(ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )		
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
		#	add viral load suppressed during infection window
		VL.cur	<- 3
		YX.m3[,lRNA.c:=NULL]
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
		YX.m3s	<- subset( YX.m3, !stage%in%c('Am','Ay','D1<=350','D1<=550','D1>550','U'))
		set(YX.m3s, NULL, 'stage', YX.m3s[, factor(as.character(stage))])
		#set(YX.m3s, NULL, 'lRNA.c', YX.m3s[, factor(as.character(lRNA.c))])
		YX.m3s[, table(stage)]
		#     ART.3.NRT ART.3.NRT.NNRT   ART.3.NRT.PI         ART.3n          ART.F          ART.I          ART.P    ART.pulse.Y 
        #	    72           1458            833            222            703            209            191            164			
		YX.m3s[, table(lRNA.c)]
		# SuA.N SuA.NA  SuA.Y 
  		#1548     72   2232
		
		ggplot(YX.m3s, aes(x = lRNA, y = score.Y)) + geom_point(size=0.4) + stat_smooth(method = "loess") + facet_grid(. ~ stage, margins=TRUE)	#suggest polynomial of order 3 for lRNA
		YX.m3s.fit10 		<- betareg(score.Y ~ lRNA.c+stage-1, link='logit', weights=w, data = YX.m3s)	#dummy variables act as baseline for lRNA, so that s fine
		#YX.m3s.fit10 		<- betareg(score.Y ~ poly(lRNA.c,3)+stage-1, link='logit', weights=w, data = YX.m3s)	#poly does not like NAs
		
		require(splines)
		YX.m3s.fit10 		<- betareg(score.Y ~ ns(lRNA.c,knot=median(lRNA.c))+stage-1, link='logit', weights=w, data = YX.m3s)	#exactly same r2 as linear model in x
		YX.m3s.fit10 		<- betareg(score.Y ~ ns(lRNA.c,knot=3)+stage-1, link='logit', weights=w, data = YX.m3s)					#unclear what ns does
		
		summary(YX.m3s.fit10)
		YX.m3s.fit10.or	<- data.table( 	F= my.or.from.logit(YX.m3s.fit10, 'stageART.F', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.F')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										I= my.or.from.logit(YX.m3s.fit10, 'stageART.I', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.I')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										P= my.or.from.logit(YX.m3s.fit10, 'stageART.P', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.P')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										Pulse= my.or.from.logit(YX.m3s.fit10, 'stageART.pulse.Y', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.pulse.Y')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										A= my.or.from.logit(YX.m3s.fit10, 'stageART.A', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.A')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										n3= my.or.from.logit(YX.m3s.fit10, 'stageART.3n', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.3n')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										NRT= my.or.from.logit(YX.m3s.fit10, 'stageART.3.NRT', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.3.NRT')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962),
										NNRT= my.or.from.logit(YX.m3s.fit10, 'stageART.3.NRT.NNRT', 'stageART.3.NRT.PI', subset(YX.m3s, stage=='ART.3.NRT.NNRT')[, sum(w)], subset(YX.m3s, stage=='ART.3.NRT.PI')[, sum(w)], 1.962) )
		#	all ORs not significant when I include lRNA
		#          F          I          P      Pulse  A         n3        NRT      NNRT
		#1: 1.361985  1.8438760  3.2542092  2.1139994 NA  2.4075422  4.2927092 1.2845132
		#2: 0.238022  0.2980641  0.5405002  0.3345607 NA  0.4160757  0.7406314 0.2573805
		#3: 7.793409 11.4065339 19.5927344 13.3577941 NA 13.9307793 24.8805967 6.4106403								
	}	
	#
	#	global ART.pulse		stratify after treatment initiation		 by 		ART.I	ART.F 	ART.A	ART.P	ART.ni
	#	
	if(1)
	{
		#	ART indication risk factors
		set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])
		set(YX.m3, YX.m3[, which(ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )		
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & contact=='No')], 'stage', 'ART.cN')
		set(YX.m3, YX.m3[, which(stage=='ART.started')], 'stage', 'ART.ni')				
		YX.m3[, table(stage)]
				
		
		ggplot(YX.m3, aes(stage, score.Y)) + geom_boxplot() + geom_jitter(size=0.4)
		#	contact No has no neg effect
		#	combining pulse has neg effect
	}
	if(1)
	{
		VL.cur	<- 3
		#
		set(YX.m3, NULL, 'stage', YX.m3[,stage.orig])
		set(YX.m3, YX.m3[, which(ART.pulse=='Yes')], 'stage', 'ART.pulse.Y' )		
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.I=='Yes')], 'stage', 'ART.I')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.P=='Yes')], 'stage', 'ART.P')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.A=='Yes')], 'stage', 'ART.A')
		set(YX.m3, YX.m3[, which(stage=='ART.started' & ART.F=='Yes')], 'stage', 'ART.F')
		# collect PoslRNA_TL and lRNA_TL etc
		YX.m3[, t.PoslRNA_T1:=NULL]
		YX.m3[, lRNA_T1:=NULL]
		YX.m3[,PoslRNA_TL:=NULL]
		YX.m3[,lRNA_TL:=NULL]
		tmp		<- unique(subset( clumsm.info, select=c(Patient, PoslRNA_T1) ))
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		YX.m3	<- merge(YX.m3, tmp, by= 't.Patient')						
		tmp		<- merge( unique(subset(clumsm.info, select=Patient)), subset(df.viro, select=c(Patient, PosRNA, lRNA)), by='Patient' )
		set( tmp, NULL, 'PosRNA', hivc.db.Date2numeric(tmp[,PosRNA]) )
		tmp		<- subset( tmp[, 	{
							tmp<- which(!is.na(lRNA))
							list(t.Patient=Patient[1], lRNA_TL=lRNA[tail(tmp,1)], PoslRNA_TL=PosRNA[tail(tmp,1)])
						}, by='Patient'], select=c(t.Patient, lRNA_TL, PoslRNA_TL) )
		YX.m3	<- merge(YX.m3, tmp, by= 't.Patient', all.x=TRUE)
		#	set VL.mx suppressed in infection window
		YX.m3[,lRNA.c:=NULL]
		YX.m3	<- merge(YX.m3, YX.m3[, {
								tmp<- which(!is.na(lRNA))
								list( lRNA.c= ifelse(length(tmp), ifelse( max(lRNA[tmp])>VL.cur, 'SuA.N', 'SuA.Y'), 'SuA.NA') )
							}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
		#	
		set(YX.m3, YX.m3[,which(stage=='ART.started' & lRNA.c=='SuA.Y')],'stage', 'ART.suA.Y' )
		set(YX.m3, YX.m3[,which(stage=='ART.started' & lRNA.c=='SuA.N')],'stage', 'ART.suA.N' )		
		set(YX.m3, YX.m3[,which(stage=='ART.started' & lRNA.c=='SuA.NA')], 'stage', 'ART.suA.NA' )	
		#	additional benefit in fast suppression?
		set(YX.m3, YX.m3[, which(stage%in%c('ART.suA.Y') & t2.vl.supp<1 )], 'stage', 'ARTs<=1')
		#set(YX.m3, YX.m3[, which(stage%in%c('ART.started','D1<=350','D1<=550','D1>550') & t2.vl.supp<1 & ( is.na(lRNA) | max(lRNA)<3 ) )], 'stage', 'ARTs<=1')
						
		YX.m3[, table(stage)]
		
		ggplot(YX.m3, aes(stage, score.Y)) + geom_boxplot() + geom_jitter(size=0.4)
	}	
	if(1)
	{
		YX.m3		<- project.athena.Fisheretal.YX.model3.stratify.ARTriskgroups(YX, df.all, cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'))	
		#X.seq		<- project.athena.Fisheretal.YX.model3.stratify.ARTriskgroups(X.seq, df.all, cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'))				
		#}
		#
		
		
		set(YX.m3, NULL, 'stage', YX.m3[, factor(as.character(stage))])
		
		#	ART indicators without number/type
		YX.m3.fit1 	<- betareg(score.Y ~ ART.pulse+ART.I+ART.F+ART.P+ART.A, link='log', weights=w, data = YX.m3)
		YX.m3.fit1 	<- betareg(score.Y ~ ART.I+ART.P+ART.A+ART.pulse+ART.F, link='log', weights=w, data = YX.m3)
		#	ART indicators and number
		cf			<- coef(betareg(score.Y ~ ART.nDrug.c+ART.pulse+ART.I+ART.F+ART.P+ART.A-1, link='log', weights=w, data = subset(YX.m3, score.Y>0.1)))		
		YX.m3.fit2 	<- betareg(score.Y ~ ART.nDrug.c+ART.pulse+ART.I+ART.F+ART.P+ART.A-1, link='log', weights=w, data = YX.m3, start=list(cf))
		#	variation in variance not significant
		cf			<- coef(betareg(score.Y ~ ART.nDrug.c+ART.pulse+ART.I+ART.F+ART.P+ART.A-1 | ART.nDrug.c+ART.pulse+ART.I+ART.F+ART.P+ART.A-1, link='log', weights=w, data = subset(YX.m3, score.Y>0.1)))
		YX.m3.fit2b <- betareg(score.Y ~ ART.nDrug.c+ART.pulse+ART.I+ART.F+ART.P+ART.A-1 | ART.nDrug.c+ART.pulse+ART.I+ART.F+ART.P+ART.A-1, link='log', weights=w, data = YX.m3, start=list(cf))		
		cf			<- coef(betareg(score.Y ~ ART.nDrug.c+ART.pulse+ART.I+ART.F+ART.P+ART.A-1 | ART.nDrug.c+ART.pulse+ART.A-1, link='log', weights=w, data = subset(YX.m3, score.Y>0.1)))
		YX.m3.fit2c <- betareg(score.Y ~ ART.nDrug.c+ART.pulse+ART.I+ART.F+ART.P+ART.A-1 | ART.nDrug.c+ART.pulse+ART.A-1, link='log', weights=w, data = YX.m3, start=list(cf))		
		cf			<- coef(betareg(score.Y ~ ART.nDrug.c+ART.pulse+ART.I+ART.F+ART.P+ART.A-1 | ART.pulse+ART.A-1, link='log', weights=w, data = subset(YX.m3, score.Y>0.1)))
		YX.m3.fit2d <- betareg(score.Y ~ ART.nDrug.c+ART.pulse+ART.I+ART.F+ART.P+ART.A-1 | ART.pulse+ART.A-1, link='log', weights=w, data = YX.m3, start=list(cf))
		require(lmtest)
		lrtest(YX.m3.fit2, YX.m3.fit2b)
		lrtest(YX.m3.fit2, YX.m3.fit2c)
		lrtest(YX.m3.fit2, YX.m3.fit2d)
		
		#	ART number of drugs and type of drugs independent
		#cf			<- coef(betareg(score.Y ~ ART.tDrug.c+ART.nDrug.c+ART.I+ART.A+ART.F+ART.P+ART.pulse-1, link='log', weights=w, data = subset(YX.m3, score.Y>0.1)))		
		#YX.m3.fit2 	<- betareg(score.Y ~ ART.tDrug.c+ART.nDrug.c+ART.I+ART.A+ART.F+ART.P+ART.pulse-1, link='log', weights=w, data = YX.m3, start=list(cf))
		
		#	ART.tpDrug.c (with pulse) independent of indicators		
		cf			<- coef(betareg(score.Y ~ ART.tpDrug.c+ART.A+ART.F+ART.P-1, link='log', weights=w, data = subset(YX.m3, score.Y>0.1)))		
		YX.m3.fit3 	<- betareg(score.Y ~ ART.tpDrug.c+ART.A+ART.F+ART.P-1, link='log', weights=w, data = YX.m3, start=list(cf))
		
		#	ART.tiDrug.c independent of indicators and of pulse	-- this is better than YX.m3.fit2
		cf			<- coef(betareg(score.Y ~ ART.tiDrug.c+ART.A+ART.F+ART.P+ART.pulse-1, link='log', weights=w, data = subset(YX.m3, score.Y>0.1)))		
		YX.m3.fit4 	<- betareg(score.Y ~ ART.tiDrug.c+ART.A+ART.F+ART.P+ART.pulse-1, link='log', weights=w, data = YX.m3, start=list(cf))
		
		#	ART.tnDrug.c independent of indicators and of pulse and of ART.I
		cf			<- coef(betareg(score.Y ~ ART.tnDrug.c+ART.I+ART.A+ART.F+ART.P+ART.pulse-1, link='log', weights=w, data = subset(YX.m3, score.Y>0.1)))		
		YX.m3.fit5 	<- betareg(score.Y ~ ART.tnDrug.c+ART.I+ART.A+ART.F+ART.P+ART.pulse-1, link='log', weights=w, data = YX.m3, start=list(cf))
		
		#	number of drugs conditional on no indication				
		YX.m3.fit6 	<- betareg(score.Y ~ ART.nstage.c-1, link='log', weights=w, data = YX.m3 )
		
		#	number and type of drugs conditional on no indication				
		YX.m3.fit7 	<- betareg(score.Y ~ ART.ntstage.c-1, link='log', weights=w, data = YX.m3 )
		
	}
	if(1)
	{
		YX.m3		<- project.athena.Fisheretal.YX.model3.stratify.ARTriskgroups(YX, df.all, cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'))	
		
		#non parametric fit
		ggplot(YX.m3, aes(x = lRNA.mx, y = score.Y)) + geom_point(size=0.4) + stat_smooth(method = "loess") + facet_grid(. ~ ART.nDrug.c, margins=TRUE)
		
		#ns spline have no order, only knots or dfs		
		ggplot(YX.m3, aes(x = lRNA.mx, y = score.Y)) + geom_point(size=0.4) + stat_smooth(method = "lm", formula= y~ns(x,2), size=1) + facet_grid(. ~ ART.nDrug.c, margins=TRUE)
		
		#linear on log scale; two change points at 3 and 5
		ggplot(YX.m3, aes(x = lRNA.mx, y = score.Y)) + geom_point(size=0.4) + stat_smooth(method = "lm", formula= y~bs(x,knots=c(3,5), degree=1)) + facet_grid(. ~ ART.nDrug.c, margins=TRUE)
		
		#linear on log scale; two change points at 3.5 and 5 to capture the bump in loess
		ggplot(YX.m3, aes(x = lRNA.mx, y = score.Y)) + geom_point(size=0.4) + stat_smooth(method = "lm", formula= y~bs(x,knots=c(3.5,5), degree=1)) + facet_grid(. ~ ART.nDrug.c, margins=TRUE)
		
		#use betaregression to determine reasonable knots
		require(gamlss)
		YX.m3b	<- na.omit(subset(YX.m3, select=c(score.Y, ART.nstage.c, lRNA.mx, w)))
		#set(YX.m3b, YX.m3b[, which(score.Y<1e-4)], 'score.Y', 1e-4)
		tmp		<- gamlss(formula= score.Y ~ bs(lRNA.mx, knots=c(3.25,4.5), degree=1), family=BE(mu.link='logit'), weights=w, data=YX.m3b)
		YX.m3b[, score.Y.gamlss:= predict(tmp, type='response')]
		tmp		<- loess(formula= score.Y ~ lRNA.mx, weights=w, span=0.8, degree=1, data=YX.m3b)
		YX.m3b[, score.Y.loess:= predict(tmp)]
		ggplot(YX.m3b, aes(x = lRNA.mx, y = score.Y)) + geom_point(aes(size=w)) + geom_smooth(aes(y= score.Y.gamlss), stat='identity') + geom_smooth(aes(y= score.Y.loess), stat='identity')
		
		#	plot relationship of pDT ~ lRNA by ART risk group
		set(YX.m3, NULL, 'stage', YX.m3[, ART.tnDrug.c])			
		formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage+ART.pulse+ART.I+ART.F+ART.P+ART.A-1'
		set(YX.m3, NULL, 'stage', YX.m3[, ART.ntstage.c])			
		formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
		#set(YX.m3, NULL, 'stage', YX.m3[, ART.ntstage.no.c])
		#formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
		#set(YX.m3, NULL, 'stage', YX.m3[, ART.ntstage.c])			
		#formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=2)+stage-1'		
		#
		YX.m3b			<- na.omit(subset(YX.m3, select=c(score.Y, lRNA.mx, stage, ART.pulse, ART.I, ART.F, ART.P, ART.A, CDCC, CD4, fw.up.med, t.period, t.Age, t.AnyPos_T1, t.RegionHospital, w)))
		YX.m3b[, table(stage)]
		
		tmp				<- project.athena.Fisheretal.betareg(YX.m3b, formula, c('score.Y','w','lRNA.mx','stage'), gamlss.BE.limit.u=c( c(0.9, 0.95, 0.975, 0.99, 0.993, 0.996, 0.998, 0.999, 1)), verbose=1 )
		betafit.rr		<- tmp$betafit.rr
		
		tmp				<- predict(betafit.rr, type='response', se.fit=TRUE)		
		YX.m3b[, score.Y.b:= tmp$fit]
		YX.m3b[, score.Y.su:= tmp$fit+tmp$se.fit]
		YX.m3b[, score.Y.sl:= tmp$fit-tmp$se.fit]
		set(YX.m3b, YX.m3b[, which(score.Y.su>1)], 'score.Y.su', 1.)
		set(YX.m3b, YX.m3b[, which(score.Y.sl<0)], 'score.Y.sl', 0.)
		
		tmp				<- YX.m3b[, list(mean.score.Y.b= mean(score.Y.b, na.rm=TRUE)), by='stage']
		setkey(tmp, mean.score.Y.b)
		set(YX.m3b, NULL, 'stage', YX.m3b[, factor( stage, levels=tmp[, stage])])
		ggplot(YX.m3b, aes(x = lRNA.mx, y = score.Y)) + geom_point(aes(size=w)) + geom_smooth(aes(y= score.Y.b), stat='identity') + geom_smooth(aes(y= score.Y.s), stat='identity') + facet_grid(. ~ stage, margins=FALSE)
		#YX.m3b[, c:= t.Age<25]				
		#YX.m3b[, c:= t.AnyPos_T1<2008]
		#YX.m3b[, c:= t.period%in%c('1','2')]
		YX.m3b[, c:= t.period]
		#YX.m3b[, c:= t.Age>45]	
		#YX.m3b[, c:= t.RegionHospital]
		#YX.m3b[, c:= fw.up.med>0.5]
		#YX.m3b[, c:= CDCC=='Yes']
		#YX.m3b[, c:= CD4<300]
		#YX.m3b[, c:= CD4>500]		
		ggplot(YX.m3b, aes(x = lRNA.mx, y = score.Y)) + geom_point(aes(size=w)) + geom_smooth(aes(y= score.Y.b), stat='identity') + geom_smooth(aes(y= score.Y.su), colour="#999999", stat='identity') + geom_smooth(aes(y= score.Y.sl), colour="#999999", stat='identity') + facet_grid(c ~ stage, margins=FALSE)
		setkey(YX.m3b, c, lRNA.mx)
		#ggplot(YX.m3b, aes(x = lRNA.mx, y = score.Y)) + geom_polygon(mapping=aes(x=c(lRNA.mx,rev(lRNA.mx)), y=c(score.Y.su,rev(score.Y.sl))), stat='identity', colour="#999999") + geom_point(aes(size=w)) + geom_smooth(aes(y= score.Y.b), stat='identity') + facet_grid(c ~ stage, margins=FALSE)
		
		
		#	plot relationship of pDT ~ lRNA by ART risk group	-- treat indicators in separate row
		set(YX.m3, NULL, 'stage', YX.m3[, ART.tnDrug.c])
		YX.m3[, ART.ind:=YX.m3[, ART.pulse=='Yes' | ART.I=='Yes' | ART.F=='Yes' | ART.P=='Yes' | ART.A=='Yes']]
		formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage+ART.ind-1'
		YX.m3b			<- na.omit(subset(YX.m3, select=c(score.Y, lRNA.mx, stage, ART.ind, ART.pulse, ART.I, ART.F, ART.P, ART.A, w)))
		tmp				<- coef(betareg(formula=formula, link='log', weights=w, data = subset(YX.m3b, score.Y>0.1)))
		betafit.rr 		<- betareg(formula=formula, link='log', weights=w, data = YX.m3b, start=list(tmp))
		YX.m3b[, score.Y.b:= predict(betafit.rr, type='response')]
		tmp				<- YX.m3b[, list(mean.score.Y.b= mean(score.Y.b, na.rm=TRUE)), by='stage']
		setkey(tmp, mean.score.Y.b)
		set(YX.m3b, NULL, 'stage', YX.m3b[, factor( stage, levels=tmp[, stage])])
		ggplot(YX.m3b, aes(x = lRNA.mx, y = score.Y)) + geom_point(aes(size=w)) + geom_smooth(aes(y= score.Y.b), stat='identity') + facet_grid(. ~ stage+ART.ind, margins=FALSE)
		
	}
	if(1)
	{
		YX.m3		<- project.athena.Fisheretal.YX.model3.stratify.ARTriskgroups(YX, df.all, return.only.ART=FALSE)
		set(YX.m3, NULL, 'stage', YX.m3[, ART.ntstage.c])
		tmp			<- YX.m3[, which(is.na(stage))]
		set(YX.m3, tmp, 'stage', YX.m3[tmp, stage.orig])
		set(YX.m3, NULL, 'stage', YX.m3[, factor(as.character(stage))])
		formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=1)+stage-1'
		#formula			<- 'score.Y ~ bs(lRNA.mx, knots=c(3.00001,4.5), degree=2)+stage-1'
		YX.m3b			<- na.omit(subset(YX.m3, select=c(score.Y, lRNA.mx, stage, ART.pulse, ART.I, ART.F, ART.P, ART.A, CDCC, CD4, fw.up.med, t.period, t.Age, t.AnyPos_T1, t.RegionHospital, w)))
		tmp				<- coef(betareg(formula=formula, link='log', weights=w, data = subset(YX.m3b, score.Y>0.1)))
		betafit.rr 		<- betareg(formula=formula, link='log', weights=w, data = YX.m3b, start=list(tmp))
		YX.m3b[, score.Y.b:= predict(betafit.rr, type='response')]
		tmp				<- YX.m3b[, list(mean.score.Y.b= mean(score.Y.b, na.rm=TRUE)), by='stage']
		setkey(tmp, mean.score.Y.b)		 
		set(YX.m3b, NULL, 'stage', YX.m3b[, factor( stage, levels=c( subset(tmp, !stage%in%c('Diag','DAy','DAm'))[,as.character(stage)],c('Diag','DAy','DAm') ))])
		ggplot(YX.m3b, aes(x = lRNA.mx, y = score.Y)) + geom_point(aes(size=w)) + geom_smooth(aes(y= score.Y.b), stat='identity') + facet_grid(. ~ stage, margins=FALSE)
		
	}
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
project.athena.Fisheretal.YX.model2.stratify.VL1stsu<- function(YX.m2, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA)
{
	#YX.m2	<- copy(YX)
	gc()
	if(is.null(YX.m2))	return(NULL)
	YX.m2[, U.score:=NULL]
	#score.Y.cut<- 1e-8
	#set(YX.m2, YX.m2[,which(score.Y<score.Y.cut)], 'score.Y', score.Y.cut)
	#score.Y.cut<- 1e-5
	#set(YX.m2, YX.m2[,which(score.Y>1-score.Y.cut)], 'score.Y', 1-score.Y.cut)
	if('score.Y'%in%colnames(YX.m2))
		set(YX.m2, NULL, 'score.Y', YX.m2[,(score.Y*(nrow(YX.m2)-1)+0.5)/nrow(YX.m2)] )
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
	#	merge PoslRNA_T1
	cat(paste('\nsetting PoslRNA_T1\n'))
	tmp			<- unique(subset( df.all, select=c(Patient, PoslRNA_T1) ))
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	YX.m2		<- merge(YX.m2, tmp, by= 't.Patient', all.x=TRUE)	
	cat(paste('\nsubset\n'))
	if('score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, score.Y, stage, CDCC, lRNA, t.isAcute, t.PoslRNA_T1, t.AnyT_T1, contact, fw.up.med, w, CD41st, CD4t  ))	
	if(!'score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, stage, CDCC, lRNA, t.isAcute, t.PoslRNA_T1, t.AnyT_T1, contact, fw.up.med, CD41st, CD4t  ))	
	gc()
	#
	#	add suppressed to 'stage'
	#
	YX.m2[, stage.orig:=stage]
	if(!is.na(plot.file.varyvl) && 'score.Y'%in%colnames(YX.m2))
	{			
		require(betareg)
		set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])		
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
		pdf(file=plot.file.varyvl, w=5, h=5)
		plot(YX.m2.VL[, VL.thresh], YX.m2.VL[, or.YN], type='p', pch=18, ylim= c( min(YX.m2.VL[, or.YN.l95], na.rm=TRUE),max(YX.m2.VL[, or.YN.u95], na.rm=TRUE))  )
		dummy	<- YX.m2.VL[,	lines( rep(VL.thresh,2), c(or.YN.l95,or.YN.u95))	, by='VL.thresh']
		dev.off()								
	}
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
				z<- which(t.lRNA<vl.suppressed)
				list(t.ARTSu_T1= ifelse(length(z), t.PosRNA[z[1]], NA_real_) )	
			}, by='t.Patient']
	YX.m2		<- merge(YX.m2, tmp, by= 't.Patient', all.x=TRUE)
	set(YX.m2, YX.m2[, which(stage=='ART.started' & !is.na(t.ARTSu_T1) & t>=t.ARTSu_T1)], 'stage', 'ART1.su.Y')
	set(YX.m2, YX.m2[, which(stage=='ART.started')], 'stage', 'ART1.su.N')
	#	set CD41st and CD4t
	tmp			<- YX.m2[, which(stage.orig=='ART.started')] 
	set(YX.m2, tmp, 'CD4t', YX.m2[tmp, stage])
	set(YX.m2, tmp, 'CD41st', YX.m2[tmp, stage])
	set(YX.m2, NULL, 'CD4t', YX.m2[, factor(as.character(CD4t))])
	set(YX.m2, NULL, 'CD41st', YX.m2[, factor(as.character(CD41st))])
	YX.m2[, t.ARTSu_T1:=NULL]
	YX.m2[, stage.orig:=NULL]
	gc()
	#
	#	add tperiod
	#
	if('t.period'%in%colnames(YX.m2))
	{
		cat(paste('\nadding CD4t.tperiod\n'))	
		YX.m2[, CD41st.tperiod:= paste(CD41st, t.period,sep='.')]
		YX.m2[, CD4t.tperiod:= paste(CD4t, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4t.tperiod', YX.m2[, factor(as.character(CD4t.tperiod))])
		set(YX.m2, NULL, 'CD41st.tperiod', YX.m2[, factor(as.character(CD41st.tperiod))])		
	}
	gc()	
	YX.m2
}
######################################################################################
project.athena.Fisheretal.v1.estimate.risk.wrap<- function(YX, X.seq, df.all, df.viro, plot.file.varyvl=NA, plot.file.or=NA, bs.n=1e3, resume=TRUE, save.file=NA, method=NA)
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
	if(resume & !is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(!resume || is.na(save.file) || inherits(readAttempt, "try-error"))
	{
		cat(paste('\nstratify by', method))
		if(method=='VL1stsu')
		{			
			YX			<- project.athena.Fisheretal.YX.model2.stratify.VL1stsu(YX, df.all, df.viro, vl.suppressed=log10(1e3), cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), plot.file.varyvl=plot.file.varyvl)	
			X.seq		<- project.athena.Fisheretal.YX.model2.stratify.VL1stsu(X.seq, df.all, df.viro, vl.suppressed=log10(1e3), cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), plot.file.varyvl=NA)				
		}
		if(method=='VLmxsu')
		{			
			#YX			<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(YX, df.all, df.viro, vl.suppressed=log10(1e3), cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), plot.file.varyvl=plot.file.varyvl)
			YX			<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(YX, df.all, df.viro, vl.suppressed=log10(1e3), cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), plot.file.varyvl=NA)
			X.seq		<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.seq, df.all, df.viro, vl.suppressed=log10(1e3), cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), plot.file.varyvl=NA)
		}
		if(method=='VLtsu')
		{			
			YX			<- project.athena.Fisheretal.YX.model2.stratify.VLt(YX, df.all, df.viro, vl.suppressed=log10(1e3), cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), plot.file.varyvl=plot.file.varyvl)	
			X.seq		<- project.athena.Fisheretal.YX.model2.stratify.VLt(X.seq, df.all, df.viro, vl.suppressed=log10(1e3), cd4.cut= c(-1, 350, 550, 5000), cd4.label=c('D1<=350','D1<=550','D1>550'), plot.file.varyvl=NA)
		}
		#
		gc()
		cat(paste('\nregression on data set'))
		YX.fit1 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX)	
		YX.fit1b 	<- betareg(score.Y ~ stage-1, link='log', weights=w, data = YX)
		#require(VGAM)
		#YX.fit1c <- vglm(score.Y ~ stage-1, link='log', tobit(Lower=0, Upper=1), data = YX)
		#tmp				<- data.table(idx=seq_len(nrow(YX)),r.logit=YX.fit1$residuals,r.log=YX.fit1b$residuals)	
		#ggplot( data=melt(tmp, id='idx'), aes(x=idx, y=value, colour=variable)) + geom_points()
		
		#	odds ratio and risk ratio
		risk.df		<- subset(data.table(risk=X.seq[,levels(stage)], risk.ref='U'), risk!='U')
		risk.df		<- rbind(risk.df, data.table(risk=X.seq[,levels(stage)][1], risk.ref=X.seq[,levels(stage)][2]))
		setkey(risk.df, risk)
		risk.prefix	<- 'stage'
		risk.ans	<- risk.df[, 	{
										tmp	<- my.or.from.logit(YX.fit1, paste(risk.prefix,risk,sep=''), paste(risk.prefix,risk.ref,sep=''), subset(YX, stage==risk)[, sum(w)], subset(YX, stage==risk.ref)[, sum(w)], 1.962)
										list(stat= 'OR', v=tmp[1], l95.asym=tmp[2], u95.asym=tmp[3])
									},by=c('risk','risk.ref')]
		tmp			<- risk.df[, 	{
										tmp	<- my.rr.from.log(YX.fit1b, paste(risk.prefix,risk,sep=''), paste(risk.prefix,risk.ref,sep=''), subset(YX, stage==risk)[, sum(w)], subset(YX, stage==risk.ref)[, sum(w)], 1.962)
										list(stat= 'RR', v=tmp[1], l95.asym=tmp[2], u95.asym=tmp[3])
									},by=c('risk','risk.ref')]
		risk.ans	<- rbind(risk.ans, tmp)	
		#	person years in infection window	
		tmp			<- X.seq[ , table( stage, useNA='ifany') ]	
		tmp			<- cbind( data.table(risk=rownames(tmp)), data.table(risk.ref='None', stat='PY', v=as.numeric(unclass(tmp)), l95.asym=NA, u95.asym=NA) )
		risk.ans	<- rbind(risk.ans, subset(tmp, select=c(risk, risk.ref, stat, v, l95.asym, u95.asym)))
		#	number of transmissions and proportion of transmissions
		tmp			<- rbind(unique(risk.df), data.table(risk='U', risk.ref='U'))
		tmp			<- tmp[, 		{
												tmp		<- my.prop.from.log(YX.fit1b, paste(risk.prefix,risk,sep=''), 1.962)
												tmp[4]	<- subset(YX, stage==risk)[, sum(w)]
												list(risk.ref='None', stat= 'N', expbeta=ifelse(tmp[4]<2*EPS, 0., tmp[1]), n=tmp[4], l95.asym=NA, u95.asym=NA)
											},by= 'risk']
		tmp[, v:= tmp[,expbeta*n]]	
		risk.ans	<- rbind(risk.ans, subset(tmp, select=c(risk, risk.ref, stat, v, l95.asym, u95.asym)))
		set(tmp, NULL, 'stat', 'P')
		set(tmp, NULL, 'v', tmp[, v/sum(expbeta*n)])
		risk.ans	<- rbind(risk.ans, subset(tmp, select=c(risk, risk.ref, stat, v, l95.asym, u95.asym)))
		#	relative infectiousness
		tmp			<- subset(risk.ans, risk.ref=='None')[, list( stat='RI', risk.ref='None', v= v[stat=='P'] / ( v[stat=='PY'] / nrow(X.seq) ), l95.asym=NA, u95.asym=NA ), by='risk']	
		risk.ans	<- rbind(risk.ans, subset(tmp, select=c(risk, risk.ref, stat, v, l95.asym, u95.asym)))	
		#	Bootstraps	on ratio and on prob
		cat(paste('\nregression on bootstrap data sets bs.n=',bs.n))
		tmp			<- lapply(seq_len(bs.n), function(bs.i)
				{
					if(bs.i%%100==0)	cat(paste('\nregression on bootstrap data sets bs.i=',bs.i))
					#bootstrap over recently infected Patient
					gc()
					tmp				<- unique(subset(YX, select=Patient))
					tmp				<- tmp[ sample( seq_len(nrow(tmp)), nrow(tmp), replace=TRUE ), ]
					YX.bs			<- merge( YX, tmp, by='Patient', allow.cartesian=TRUE )				
					YX.fit1.bs 		<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.bs)	
					YX.fit1b.bs 	<- betareg(score.Y ~ stage-1, link='log', weights=w, data = YX.bs)
					#	odds ratio and risk ratio
					risk.ans.bs		<- risk.df[, 	list(stat='OR', v=my.or.from.logit(YX.fit1.bs, paste(risk.prefix,risk,sep=''), paste(risk.prefix,risk.ref,sep=''), subset(YX.bs, stage==risk)[, sum(w)], subset(YX.bs, stage==risk.ref)[, sum(w)], 1.962)[1]),	by=c('risk','risk.ref')]
					tmp				<- risk.df[, 	list(stat='RR', v=my.rr.from.log(YX.fit1b.bs, paste(risk.prefix,risk,sep=''), paste(risk.prefix,risk.ref,sep=''), subset(YX.bs, stage==risk)[, sum(w)], subset(YX.bs, stage==risk.ref)[, sum(w)], 1.962)[1]),	by=c('risk','risk.ref')]
					risk.ans.bs		<- rbind(risk.ans.bs, tmp)
					#	for proportions
					tmp				<- rbind(unique(risk.df), data.table(risk='U', risk.ref='U'))				 
					risk.ans.bs		<- rbind(risk.ans.bs, tmp[, 	list(risk.ref='None', stat= 'prob', v=my.prop.from.log(YX.fit1b.bs, paste(risk.prefix,risk,sep=''), 1.962)[1]), by= 'risk'])
					risk.ans.bs		<- rbind(risk.ans.bs, tmp[, 	list(risk.ref='None', stat= 'n', v=subset(YX.bs, stage==risk)[, sum(w)]),by= 'risk'])				
					risk.ans.bs[, bs:=bs.i]
				})
		risk.ans.bs	<- do.call('rbind',tmp)
		tmp			<- risk.ans.bs[, list(stat='N', risk.ref='None', v= v[stat=='prob']*v[stat=='n'] ), by=c('bs','risk')]
		risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(risk, risk.ref, stat, v, bs)))
		tmp			<- tmp[,	list(stat='P', risk=risk, risk.ref=risk.ref, v= v/sum(v, na.rm=TRUE)),	by='bs']
		risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(risk, risk.ref, stat, v, bs)))	
		tmp			<- merge( subset(tmp, stat=='P'), subset(risk.ans, stat=='PY', select=c(risk, v)), by='risk' )
		tmp[, v:= v.x/(v.y/nrow(X.seq))]
		set(tmp, NULL, 'stat', 'RI')
		risk.ans.bs	<- rbind(risk.ans.bs, subset(tmp, select=c(risk, risk.ref, stat, v, bs)))
		tmp			<- risk.ans.bs[,	list(l95.bs=quantile(v, prob=0.025, na.rm=TRUE), u95.bs=quantile(v, prob=0.975, na.rm=TRUE)), by=c('risk','risk.ref','stat')]
		risk.ans	<- merge(risk.ans, tmp, by=c('risk','risk.ref','stat'), all.x=TRUE)
		#	
		setkey(risk.ans, stat, risk, risk.ref)
		
		if(0)
		{
			tmp			<- cooks.distance(YX.fit1)
			plot(seq_along(tmp), tmp, type='h', col=YX[, col], xlab='index', ylab='Cooks D')
			legend('topright', bty='n', border=NA, legend= YX[, levels(stage)], fill=YX[, unique(col)])
			tmp			<- residuals(YX.fit1)
			plot(seq_along(tmp), tmp, type='p', pch=18, col=YX[, col], xlab='index', ylab='std residuals')
			legend('topright', bty='n', border=NA, legend= YX[, levels(stage)], fill=YX[, unique(col)])		
		}
		if(!is.na(plot.file.or))
		{
			cat(paste('\nplot to',plot.file.or))
			pdf(file=plot.file.or, w=5, h=5)
			tmp			<- data.table(stage= YX[, levels(stage)], col=sapply( rainbow(YX[, nlevels(stage)]), my.fade.col, alpha=0.5) )
			set(tmp, NULL, 'stage', tmp[, factor(stage)])
			YX		<- merge(YX, tmp, by='stage')
			plot( seq_len(nrow(YX)), YX[, score.Y], pch=18, col= YX[, col], cex=YX[, w^0.4])
			tmp				<- subset(YX, select=stage)
			lines(seq_len(nrow(tmp)), predict(YX.fit1, tmp, type='response'))	
			start			<- c( YX.fit1$coef$mean + head( 2*sqrt(diag(vcov(YX.fit1))), length(YX.fit1$coef$mean)), YX.fit1$coef$precision)
			YX.fit1.sup	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX, start=start, hessian=TRUE, maxit=0)		
			start			<- c( YX.fit1$coef$mean - head( 2*sqrt(diag(vcov(YX.fit1))), length(YX.fit1$coef$mean)), YX.fit1$coef$precision)
			YX.fit1.slw	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data=YX, start=start, hessian=TRUE, maxit=0)		
			polygon(	c( seq_len(nrow(tmp)),rev(seq_len(nrow(tmp))) ), 
					c( predict(YX.fit1.sup, tmp, type='response'), rev(predict(YX.fit1.slw, tmp, type='response'))), 
					border=NA, col=my.fade.col('black',0.5) )
			legend('bottomright', bty='n', border=NA, legend= YX[, levels(stage)], fill=YX[, unique(col)])
			YX[, col:=NULL]
			dev.off()
		}
		if(!is.na(save.file))
		{
			cat(paste('\nsave to file', save.file))
			save(risk.ans, risk.ans.bs, YX.fit1, YX.fit1b, file=save.file)
		}
	}
	list(risk=risk.ans, or.fit=YX.fit1, rr.fit=YX.fit1b, risk.bs=risk.ans.bs)
}
######################################################################################
project.athena.Fisheretal.YX.model2.stratify.VLt<- function(YX.m2, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA, plot.file.or=NA )
{
	#YX.m2	<- copy(YX)
	gc()
	if(is.null(YX.m2))	return(NULL)
	YX.m2[, U.score:=NULL]
	#score.Y.cut<- 1e-8
	#set(YX.m2, YX.m2[,which(score.Y<score.Y.cut)], 'score.Y', score.Y.cut)
	#score.Y.cut<- 1e-5
	#set(YX.m2, YX.m2[,which(score.Y>1-score.Y.cut)], 'score.Y', 1-score.Y.cut)
	if('score.Y'%in%colnames(YX.m2))
		set(YX.m2, NULL, 'score.Y', YX.m2[,(score.Y*(nrow(YX.m2)-1)+0.5)/nrow(YX.m2)] )
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
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, score.Y, stage, CDCC, lRNA, t.isAcute, t.PoslRNA_T1, t.AnyT_T1, contact, fw.up.med, w, lRNA_TL, PoslRNA_TL, CD41st, CD4t  ))	
	if(!'score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, stage, CDCC, lRNA, t.isAcute, t.PoslRNA_T1, t.AnyT_T1, contact, fw.up.med, lRNA_TL, PoslRNA_TL, CD41st, CD4t  ))	
	#
	#	add suppressed to 'stage'
	#
	YX.m2[, stage.orig:=stage]
	gc()
	if(!is.na(plot.file.varyvl))
	{			
		require(betareg)
		print(YX.m2)
		print(YX.m2[,table(stage.orig)])
		set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
		VL.win		<- c( seq(400, 1000, by=100), seq(1250, 1e4, by=250) ) 
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
					YX.m2.fit4 	<- betareg(score.Y ~ stage-1, link='logit', weights=w, data = YX.m2)
					YX.m2[, lRNA.c:=NULL]
					
					ans			<- c(		YX.m2[, length(which(stage=='ART.su.Y'))], YX.m2[, length(which(stage=='ART.su.N'))],
							subset(YX.m2, stage=='ART.su.Y')[, sum(w)], subset(YX.m2, stage=='ART.su.N')[, sum(w)],
							coef(YX.m2.fit4)['stageART.su.Y'], coef(YX.m2.fit4)['stageART.su.N'], coef(YX.m2.fit4)['stageU'],	
							sqrt(diag(vcov(YX.m2.fit4)))[c('stageART.su.Y','stageART.su.N')],
							my.or.from.logit(YX.m2.fit4, 'stageART.su.Y', 'stageART.su.N', subset(YX.m2, stage=='ART.su.Y')[, sum(w)], subset(YX.m2, stage=='ART.su.N')[, sum(w)], 1.962),						
							my.or.from.logit(YX.m2.fit4, 'stageART.su.Y', 'stageU', subset(YX.m2, stage=='ART.su.Y')[, sum(w)], subset(YX.m2, stage=='U')[, sum(w)], 1.962),
							logLik(YX.m2.fit4), YX.m2.fit4$pseudo.r.squared, 10^VL.cur	)
					names(ans)	<- c('n.su.Y','n.su.N','w.su.Y','w.su.N','coef.su.Y','coef.su.N','coef.U','coef.su.Y.sd','coef.su.N.sd','or.YN','or.YN.l95','or.YN.u95','or.YU','or.YU.l95','or.YU.u95','lkl','r2','VL.thresh')
					ans
				})
		YX.m2.VL	<- as.data.table(t(YX.m2.VL))
		ylim		<- c( min(YX.m2.VL[, or.YN.l95], na.rm=TRUE),max(YX.m2.VL[, or.YN.u95], na.rm=TRUE))
		print(YX.m2.VL)
		print(ylim)
		pdf(file=plot.file.varyvl, w=5, h=5)
		plot(YX.m2.VL[, VL.thresh], YX.m2.VL[, or.YN], type='p', pch=18, ylim= ylim  )
		dummy	<- YX.m2.VL[,	lines( rep(VL.thresh,2), c(or.YN.l95,or.YN.u95))	, by='VL.thresh']
		dev.off()
	}
	#
	#	using given VL threshold 
	#		
	set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
	set(YX.m2, YX.m2[,which(stage=='ART.started' & !is.na(lRNA) & lRNA<=vl.suppressed)],'stage', 'ART.su.Y' )
	set(YX.m2, YX.m2[,which(stage=='ART.started' & !is.na(lRNA) & lRNA>vl.suppressed)],'stage', 'ART.su.N' )		
	set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & !is.na(lRNA_TL) & lRNA_TL<=vl.suppressed)],'stage', 'ART.su.Y' )	
	set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & !is.na(lRNA_TL) & lRNA_TL>vl.suppressed)],'stage', 'ART.su.N' )
	set(YX.m2, YX.m2[,which(stage=='ART.started' & is.na(lRNA) & t>t.PoslRNA_T1 & contact=='Yes' & lRNA>vl.suppressed)],'stage', 'ART.su.N' )		
	set(YX.m2, YX.m2[,which(stage=='ART.started')],'stage', 'ART.vlNA' )
	#	set CD41st and CD4t
	tmp			<- YX.m2[, which(stage.orig=='ART.started')] 
	set(YX.m2, tmp, 'CD4t', YX.m2[tmp, stage])
	set(YX.m2, tmp, 'CD41st', YX.m2[tmp, stage])
	set(YX.m2, NULL, 'CD4t', YX.m2[, factor(as.character(CD4t))])
	set(YX.m2, NULL, 'CD41st', YX.m2[, factor(as.character(CD41st))])
	YX.m2[, stage.orig:=NULL]	
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
		YX.m2[, CD4t.tperiod:= paste(CD4t, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4t.tperiod', YX.m2[, factor(as.character(CD4t.tperiod))])
		set(YX.m2, NULL, 'CD41st.tperiod', YX.m2[, factor(as.character(CD41st.tperiod))])		
	}
	gc()	
	YX.m2
}
######################################################################################
project.athena.Fisheretal.YX.model2.stratify.VLgm<- function(YX.m2, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA, plot.file.or=NA )
{
	#YX.m2	<- copy(YX)
	gc()
	if(is.null(YX.m2))	return(NULL)
	YX.m2[, U.score:=NULL]
	#score.Y.cut<- 1e-8
	#set(YX.m2, YX.m2[,which(score.Y<score.Y.cut)], 'score.Y', score.Y.cut)
	#score.Y.cut<- 1e-5
	#set(YX.m2, YX.m2[,which(score.Y>1-score.Y.cut)], 'score.Y', 1-score.Y.cut)
	if('score.Y'%in%colnames(YX.m2))
		set(YX.m2, NULL, 'score.Y', YX.m2[,(score.Y*(nrow(YX.m2)-1)+0.5)/nrow(YX.m2)] )
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
	cat(paste('\nsubset\n'))
	if('score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, score.Y, stage, CDCC, lRNA, t.isAcute, t.AnyT_T1, contact, fw.up.med, t.period, w, CD41st, CD4t  ))	
	if(!'score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, stage, CDCC, lRNA, t.isAcute, t.AnyT_T1, contact, fw.up.med, t.period, CD41st, CD4t  ))	
	gc()
	#
	#	add suppressed to 'stage'
	#
	YX.m2[, stage.orig:=stage]	
	if(!is.na(plot.file.varyvl))
	{			
		require(betareg)
		set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
		VL.win		<- c( seq(400, 1000, by=100), seq(1250, 1e4, by=250) ) 
		YX.m2.VL	<- sapply(VL.win, function(VL.cur)
				{
					cat(paste('\nprocess VL=', VL.cur))									
					VL.cur	<- log10(VL.cur)
					set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
					YX.m2[,lRNA.c:=NULL]
					YX.m2	<- merge(YX.m2, YX.m2[, {
										tmp<- which(!is.na(lRNA))
										list( lRNA.c= ifelse(length(tmp), ifelse( max(lRNA[tmp])>VL.cur, 'SuA.N', 'SuA.Y'), 'SuA.NA') )
									}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))					
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
		tmp			<- c( min(YX.m2.VL[, or.YN.l95], na.rm=TRUE),max(YX.m2.VL[, or.YN.u95], na.rm=TRUE))
		pdf(file=plot.file.varyvl, w=5, h=5)
		plot(YX.m2.VL[, VL.thresh], YX.m2.VL[, or.YN], type='p', pch=18, ylim=tmp  )
		dummy	<- YX.m2.VL[,	lines( rep(VL.thresh,2), c(or.YN.l95,or.YN.u95))	, by='VL.thresh']
		dev.off()
	}
	#
	#	using given VL threshold 
	#	
	gc()
	cat(paste('\nadding lRNA.mx\n'))
	set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
	YX.m2[,lRNA.c:=NULL]
	YX.m2	<- merge(YX.m2, YX.m2[, {
						tmp<- which(!is.na(lRNA))
						list( lRNA.c= ifelse(length(tmp), ifelse( max(lRNA[tmp])>vl.suppressed, 'SuA.N', 'SuA.Y'), 'SuA.NA') )
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)
	gc()
	set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.c=='SuA.Y')],'stage', 'ART.suA.Y' )
	set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.c=='SuA.N')],'stage', 'ART.suA.N' )
	set(YX.m2, YX.m2[,which(stage=='ART.started')],'stage', 'ART.vlNA' )
	#	set CD41st and CD4t
	tmp			<- YX.m2[, which(stage.orig=='ART.started')] 
	set(YX.m2, tmp, 'CD4t', YX.m2[tmp, stage])
	set(YX.m2, tmp, 'CD41st', YX.m2[tmp, stage])
	set(YX.m2, NULL, 'CD4t', YX.m2[, factor(as.character(CD4t))])
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
		YX.m2[, CD4t.tperiod:= paste(CD4t, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4t.tperiod', YX.m2[, factor(as.character(CD4t.tperiod))])
		set(YX.m2, NULL, 'CD41st.tperiod', YX.m2[, factor(as.character(CD41st.tperiod))])		
	}
	gc()
	YX.m2
}
######################################################################################
project.athena.Fisheretal.YX.model2.stratify.VLmxwindow<- function(YX.m2, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA, plot.file.or=NA )
{
	#YX.m2	<- copy(YX)
	gc()
	if(is.null(YX.m2))	return(NULL)
	YX.m2[, U.score:=NULL]
	#score.Y.cut<- 1e-8
	#set(YX.m2, YX.m2[,which(score.Y<score.Y.cut)], 'score.Y', score.Y.cut)
	#score.Y.cut<- 1e-5
	#set(YX.m2, YX.m2[,which(score.Y>1-score.Y.cut)], 'score.Y', 1-score.Y.cut)
	if('score.Y'%in%colnames(YX.m2))
		set(YX.m2, NULL, 'score.Y', YX.m2[,(score.Y*(nrow(YX.m2)-1)+0.5)/nrow(YX.m2)] )
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
	cat(paste('\nsubset\n'))
	if('score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, score.Y, stage, CDCC, lRNA, t.isAcute, t.AnyT_T1, contact, fw.up.med, t.period, w, CD41st, CD4t  ))	
	if(!'score.Y'%in%colnames(YX.m2))
		YX.m2	<- subset(YX.m2, select=c(t, t.Patient, Patient, stage, CDCC, lRNA, t.isAcute, t.AnyT_T1, contact, fw.up.med, t.period, CD41st, CD4t  ))	
	gc()
	#
	#	add suppressed to 'stage'
	#
	YX.m2[, stage.orig:=stage]	
	if(!is.na(plot.file.varyvl))
	{			
		require(betareg)
		set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
		VL.win		<- c( seq(400, 1000, by=100), seq(1250, 1e4, by=250) ) 
		YX.m2.VL	<- sapply(VL.win, function(VL.cur)
				{
					cat(paste('\nprocess VL=', VL.cur))									
					VL.cur	<- log10(VL.cur)
					set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
					YX.m2[,lRNA.c:=NULL]
					YX.m2	<- merge(YX.m2, YX.m2[, {
										tmp<- which(!is.na(lRNA))
										list( lRNA.c= ifelse(length(tmp), ifelse( max(lRNA[tmp])>VL.cur, 'SuA.N', 'SuA.Y'), 'SuA.NA') )
									}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))					
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
		tmp			<- c( min(YX.m2.VL[, or.YN.l95], na.rm=TRUE),max(YX.m2.VL[, or.YN.u95], na.rm=TRUE))
		pdf(file=plot.file.varyvl, w=5, h=5)
		plot(YX.m2.VL[, VL.thresh], YX.m2.VL[, or.YN], type='p', pch=18, ylim=tmp  )
		dummy	<- YX.m2.VL[,	lines( rep(VL.thresh,2), c(or.YN.l95,or.YN.u95))	, by='VL.thresh']
		dev.off()
	}
	#
	#	using given VL threshold 
	#	
	gc()
	cat(paste('\nadding lRNA.mx\n'))
	set(YX.m2, NULL, 'stage', YX.m2[, stage.orig])
	YX.m2[,lRNA.c:=NULL]
	YX.m2	<- merge(YX.m2, YX.m2[, {
						tmp<- which(!is.na(lRNA))
						list( lRNA.c= ifelse(length(tmp), ifelse( max(lRNA[tmp])>vl.suppressed, 'SuA.N', 'SuA.Y'), 'SuA.NA') )
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)
	gc()
	set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.c=='SuA.Y')],'stage', 'ART.suA.Y' )
	set(YX.m2, YX.m2[,which(stage=='ART.started' & lRNA.c=='SuA.N')],'stage', 'ART.suA.N' )
	set(YX.m2, YX.m2[,which(stage=='ART.started')],'stage', 'ART.vlNA' )
	#	set CD41st and CD4t
	tmp			<- YX.m2[, which(stage.orig=='ART.started')] 
	set(YX.m2, tmp, 'CD4t', YX.m2[tmp, stage])
	set(YX.m2, tmp, 'CD41st', YX.m2[tmp, stage])
	set(YX.m2, NULL, 'CD4t', YX.m2[, factor(as.character(CD4t))])
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
		YX.m2[, CD4t.tperiod:= paste(CD4t, t.period,sep='.')]
		set(YX.m2, NULL, 'CD4t.tperiod', YX.m2[, factor(as.character(CD4t.tperiod))])
		set(YX.m2, NULL, 'CD41st.tperiod', YX.m2[, factor(as.character(CD41st.tperiod))])		
	}
	gc()
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
	stopifnot(class(df.immu$Patient)=='character',class(df.tpairs$t.Patient)=='character')
	immu	<- subset( df.immu, select=c(Patient, PosCD4, CD4) )
	setnames(immu, 'Patient','t.Patient')
	immu	<- merge(immu, unique(subset(df.tpairs, select=t.Patient)), by='t.Patient')	
	set(immu, NULL, 'PosCD4', hivc.db.Date2numeric(immu[,PosCD4]))
	tmp		<- subset(immu, select=c(t.Patient, PosCD4))[, list(ts=min(PosCD4), te=max(PosCD4)), by='t.Patient']
	set(tmp, NULL, 'ts', tmp[, floor(ts) + round( (ts%%1)*100 %/% (t.period*100) ) * t.period] )
	set(tmp, NULL, 'te', tmp[, floor(te) + round( (te%%1)*100 %/% (t.period*100) ) * t.period] )
	tmp		<- tmp[, list(t= seq(ts, te, by=t.period)),by='t.Patient']
	setnames(tmp, 't.Patient', 'tt.Patient')
	setkey(tmp, tt.Patient)
	immu	<- immu[, {
				z		<- subset(tmp, tt.Patient==t.Patient[1])	
				if(length(PosCD4)<2)
				{
					y	<- rep(NA, nrow(z))
					y[z[,which( PosCD4<=t+t.period )[1]]]<- CD4
				}
				else
					y	<- approx(PosCD4 , CD4, xout=z[,t]+t.period/2, yleft=NA_real_, yright=NA_real_, rule=2)$y
				list(t=z[,t], CD4=y)
			},by='t.Patient']
	set(immu, NULL, 'CD4', immu[, as.integer(round(CD4))])
	immu		
}
######################################################################################
project.athena.Fisheretal.X.viro<- function(df.tpairs, df.viro, t.period=0.25, lRNA.cutoff= log10(1e3))
{
	stopifnot(class(df.viro$Patient)=='character',class(df.tpairs$t.Patient)=='character') 
	viro	<- subset( df.viro, select=c(Patient, PosRNA, lRNA) )
	setnames(viro, 'Patient','t.Patient')
	viro	<- merge(viro, unique(subset(df.tpairs, select=t.Patient)), by='t.Patient')	
	set(viro, NULL, 'PosRNA', hivc.db.Date2numeric(viro[,PosRNA]))
	tmp		<- subset(viro, select=c(t.Patient, PosRNA))[, list(ts=min(PosRNA), te=max(PosRNA)), by='t.Patient']
	set(tmp, NULL, 'ts', tmp[, floor(ts) + round( (ts%%1)*100 %/% (t.period*100) ) * t.period] )
	set(tmp, NULL, 'te', tmp[, floor(te) + round( (te%%1)*100 %/% (t.period*100) ) * t.period] )
	tmp		<- tmp[, list(t= seq(ts, te, by=t.period)),by='t.Patient']
	
	setnames(tmp, 't.Patient', 'tt.Patient')
	setkey(tmp, tt.Patient)
	#merge(data.table(t.Patient=tmp[, unique(tt.Patient)]), viro, by='t.Patient')	
	viro	<- viro[, {
						z		<- subset(tmp, tt.Patient==t.Patient[1])						#fails if t.Patient is factor because t.Patient[1] gets converted to integer
						if(length(PosRNA)<2)
						{
							y	<- rep(NA, nrow(z))
							y[z[,which( PosRNA<=t+t.period )[1]]]<- lRNA
							y2	<- y
						}
						else
						{
							y	<- approx(PosRNA , lRNA, xout=z[,t]+t.period/2, yleft=NA_real_, yright=NA_real_, rule=2, method='linear')$y
							y2	<- approx(PosRNA , lRNA, xout=z[,t]+t.period/2, yleft=NA_real_, yright=NA_real_, rule=2, method='constant')$y
						}
						list(t=z[,t], lRNA=y, lRNAc=y2)
					},by='t.Patient']	
	if(!is.na(lRNA.cutoff))
	{
		set(viro, viro[, which(lRNA<lRNA.cutoff)], 'lRNA', 0.)
		set(viro, viro[, which(lRNA>=lRNA.cutoff)], 'lRNA', 1.)
		set(viro, NULL, 'lRNA', as.integer(viro[,lRNA]))
		set(viro, viro[, which(lRNAc<lRNA.cutoff)], 'lRNAc', 0.)
		set(viro, viro[, which(lRNAc>=lRNA.cutoff)], 'lRNAc', 1.)
		set(viro, NULL, 'lRNAc', as.integer(viro[,lRNAc]))
	}
	viro
}
######################################################################################
project.athena.Fisheretal.Y.infectiontime<- function(YX.tpairs, df.all, predict.t2inf, t2inf.args, t.period=0.25, ts.min=1980, score.min=0.2, score.set.value=1, method.infectiontime='for.infected')
{
	#ts.min=1980; score.min=0.01; score.set.value=1
	stopifnot(method.infectiontime%in%c('for.infected','for.transmitter'))
	if(method.infectiontime=='for.infected')
		b4care	<- merge(unique(subset( YX.tpairs, select=Patient )), unique(subset(df.all, select=c(Patient, DateBorn, AnyPos_T1, isAcute, NegT, PosCD4_T1, CD4_T1, AnyT_T1))), by='Patient')
	if(method.infectiontime=='for.transmitter')
		b4care	<- merge(data.table(Patient=YX.tpairs[, unique(t.Patient)]), unique(subset(df.all, select=c(Patient, DateBorn, AnyPos_T1, isAcute, NegT, PosCD4_T1, CD4_T1, AnyT_T1))), by='Patient')
	#
	tmp		<- b4care[, list(AnyPos_a=AnyPos_T1-DateBorn, ts=ifelse(is.na(NegT), AnyPos_T1-10, NegT), te=AnyPos_T1), by='Patient']
	b4care	<- merge(b4care, tmp, by='Patient')
	#check if ts before 1980 and if so clip
	set(b4care, b4care[,which(ts<ts.min)], 'ts', ts.min )	
	set(b4care, NULL, 'ts', b4care[, floor(ts) + round( (ts%%1)*100 %/% (t.period*100) ) * t.period ] )
	set(b4care, NULL, 'te', b4care[, floor(te) + round( (te%%1)*100 %/% (t.period*100) ) * t.period - t.period] )	#last time undiagnosed is just before first time period diagnosed
	#apply predict.t2inf	
	t2inf	<- predict.t2inf(seq(0,10,t.period)*365, b4care, t2inf.args)
	t2inf	<- merge( subset( t2inf, score>score.min ), subset(b4care, select=c(Patient,ts,te)), by='Patient' )
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
project.athena.Fisheretal.v3.Y.coal<- function(df.tpairs, clumsm.info, X.pt, cluphy, cluphy.info, cluphy.map.nodectime, coal.t.Uscore.min=0.01, t.period= 0.25, coal.within.inf.grace= 0.25, save.file=NA, resume=0 )
{
	#df.tpairs	<- subset(X.pt, select=c(t.Patient, Patient, cluster, FASTASampleCode, t.FASTASampleCode))
	#setkey(df.tpairs, FASTASampleCode, t.FASTASampleCode)
	#df.tpairs	<- unique(df.tpairs)
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
		setkey(df.tpairs, cluster)	
		#	select tpairs for which dated phylogenies are available
		tmp					<- setdiff( df.tpairs[,unique(cluster)], cluphy.info[, unique(cluster)] )
		cat(paste('\nnumber of clusters missing for analysis of pot transmitters, n=',length(tmp)))
		df.tpairs.mrca		<- df.tpairs[J(cluphy.info[, unique(cluster)]),]	
		cat(paste('\nnumber of pot transmitters for which dated phylogenies are available, n=',df.tpairs.mrca[,length(unique(t.Patient))]))		
		#	compute mrcas of tpairs 
		tmp					<- df.tpairs.mrca[,	list(node=getMRCA(cluphy, c(FASTASampleCode, t.FASTASampleCode))), by=c('FASTASampleCode','t.FASTASampleCode')]
		if(nrow(tmp)!=nrow(df.tpairs.mrca))	stop('unexpected length of tmp')
		df.tpairs.mrca		<- merge(df.tpairs.mrca, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))
		#	compute posterior prob that coalescence is after last NegT of transmitter or after 80% cutoff from U.score 
		#			posterior prob that coalescence is before AnyPos_T1 of individual in denominator population
		tmp					<- unique( subset( clumsm.info, select=c(Patient, NegT)) )
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		df.tpairs.mrca		<- merge(df.tpairs.mrca, tmp, by='t.Patient')		
		
		tmp	<- subset(X.pt, !is.na( U.score ) & U.score>coal.t.Uscore.min)
		setkey(tmp, t.Patient, t)
		tmp					<- tmp[, list(t.UT= min(t)), by='t.Patient']
		df.tpairs.mrca		<- merge(df.tpairs.mrca, tmp, by='t.Patient')
		tmp					<- df.tpairs.mrca[, list(t.queryT= max(t.NegT, t.UT, na.rm=TRUE)), by=c('FASTASampleCode','t.FASTASampleCode')]
		df.tpairs.mrca		<- merge(df.tpairs.mrca, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))
		tmp					<- subset( df.tpairs.mrca, select=c(cluster, node, FASTASampleCode, t.FASTASampleCode, t.queryT) )
		tmp					<- merge(tmp, subset(clumsm.info, select=c(FASTASampleCode, AnyPos_T1)), by='FASTASampleCode')
		set(tmp, NULL, 'AnyPos_T1', tmp[, AnyPos_T1+coal.within.inf.grace])
		tmp					<- merge( cluphy.map.nodectime, tmp, by=c('cluster','node'), allow.cartesian=TRUE)	
		coal				<- tmp[,  {
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
project.athena.Fisheretal.Y.coal<- function(YX.tpairs, df.all, Y.U, cluphy, cluphy.info, cluphy.map.nodectime, coal.t.Uscore.min=0.01, coal.within.inf.grace= 0.25, t.period= 0.25, save.file=NA, resume=0 )
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
		tmp					<- df.tpairs.mrca[, list(t.queryT= max(t.NegT, t.UT, na.rm=TRUE)), by=c('FASTASampleCode','t.FASTASampleCode')]
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
project.athena.Fisheretal.v1.Y.coal.and.inf<- function(df.tpairs, clumsm.info, cluphy, cluphy.info, cluphy.map.nodectime, t.period= 0.25, save.file=NA )
{
	#df.tpairs	<- subset(X.pt, select=c(t.Patient, Patient, cluster, FASTASampleCode, t.FASTASampleCode))
	#setkey(df.tpairs, FASTASampleCode, t.FASTASampleCode)
	#df.tpairs	<- unique(df.tpairs)
	#
	if(!is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(is.na(save.file) || inherits(readAttempt, "try-error"))
	{
		setkey(df.tpairs, cluster)	
		#	select tpairs for which dated phylogenies are available
		tmp					<- setdiff( df.tpairs[,unique(cluster)], cluphy.info[, unique(cluster)] )
		cat(paste('\nnumber of clusters missing for analysis of pot transmitters, n=',length(tmp)))
		df.tpairs			<- df.tpairs[J(cluphy.info[, unique(cluster)]),]	
		cat(paste('\nnumber of pot transmitters for which dated phylogenies are available, n=',df.tpairs[,length(unique(t.Patient))]))		
		#	compute mrcas of tpairs 
		tmp					<- df.tpairs[,	list(node=getMRCA(cluphy, c(FASTASampleCode, t.FASTASampleCode))), by=c('FASTASampleCode','t.FASTASampleCode')]
		if(nrow(tmp)!=nrow(df.tpairs))	stop('unexpected length of tmp')
		df.tpairs			<- merge(df.tpairs, tmp, by=c('FASTASampleCode','t.FASTASampleCode'))
		#	compute posterior prob that coalescence is after last NegT and before TipT of individual in denominator population
		tmp					<- unique( subset( clumsm.info, select=c(Patient, NegT)) )
		setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
		df.tpairs			<- merge(df.tpairs, tmp, by='t.Patient')		
		tmp					<- subset( df.tpairs, select=c(cluster, FASTASampleCode, t.FASTASampleCode, node,  t.NegT) )
		tmp					<- merge(tmp, subset(clumsm.info, select=c(FASTASampleCode, AnyPos_T1)), by='FASTASampleCode')
		tmp					<- merge( cluphy.map.nodectime, tmp, by=c('cluster','node'), allow.cartesian=TRUE)	
		coal				<- tmp[,  {
											z		<- 1-approx(q ,cdf, xout=c(t.NegT[1], AnyPos_T1[1]), yleft=0., yright=1., rule=2)$y
											#z[1] is prob survive in transmitter, z[2] is prob survive in infected
											list(coal.after.t.NegT= z[1], coal.after.i.AnyPos_T1=z[2])
										} , by=c('cluster','FASTASampleCode','t.FASTASampleCode')]
		coal[, score.coal:= coal.after.t.NegT-coal.after.i.AnyPos_T1]
		set(coal, coal[, which(score.coal<0)], 'score.coal', 0.)
		set(coal, NULL, 'score.coal', coal[, round(score.coal, d=3)])						
		coal				<- subset(coal, select=c(FASTASampleCode, t.FASTASampleCode, score.coal))
		#	earliest time of infection: either lowest coal or NegT(infected)
		tmp					<- merge( subset( df.tpairs, select=c(cluster,  node, FASTASampleCode, t.FASTASampleCode) ), subset( clumsm.info, select=c(FASTASampleCode, NegT, PosSeqT, AnyPos_T1)), by='FASTASampleCode' )
		inf					<- merge(tmp, cluphy.map.nodectime[, list(min.coalT= min(q)), by=c('cluster','node')], by=c('cluster','node'))
		tmp					<- inf[,which(is.na(NegT))]
		set(inf, tmp, 'NegT', inf[tmp, min.coalT])
		setnames(inf, 'NegT', 'InfT')	
		#	get list of time periods for every infected
		set(inf, NULL, 'InfT', inf[, floor(InfT) + round( (InfT%%1)*100 %/% (t.period*100) ) * t.period] )								#min infection time period
		set(inf, NULL, 'PosSeqT', inf[, floor(PosSeqT) + round( (PosSeqT%%1)*100 %/% (t.period*100) ) * t.period] )	
		set(inf, NULL, 'AnyPos_T1', inf[, floor(AnyPos_T1) + round( (AnyPos_T1%%1)*100 %/% (t.period*100) ) * t.period] )				#time period in which the recipient is 'diagnosed'
		
		tmp	<- merge(df.tpairs, subset(inf, InfT>=AnyPos_T1), by=c('FASTASampleCode','t.FASTASampleCode'))
		setkey(tmp, cluster.x)
		
		cat(paste('\nnumer of sequences with zero infection probability because the viral lineages coalesce after time of diagnosis of the putative infected, n=',nrow(subset(inf, InfT>=AnyPos_T1)))) 
		inf					<- subset(inf, InfT<AnyPos_T1)
		inf					<- inf[, list(InfT= seq(InfT, AnyPos_T1-t.period, by=t.period), cluster=cluster, node=node), by=c('FASTASampleCode','t.FASTASampleCode')]
		#	evaluated coal at InfT
		setkey(cluphy.map.nodectime, cluster, node)
		setnames(inf, c('cluster','node'), c('i.cluster','i.node'))
		inf					<- inf[, {
										tmp  <- subset(cluphy.map.nodectime, cluster==i.cluster[1] & node==i.node[1])
										InfV <- approx(tmp[,q] ,tmp[,cdf], xout=InfT+t.period/2, yleft=0., yright=1., rule=2)$y		#never includes TMRCA mass after time of diagnosis of infected
										list(t=InfT, score.inf=InfV)
									},by=c('FASTASampleCode','t.FASTASampleCode')]
		set(inf, NULL, 'score.inf', inf[, round(score.inf,d=3)])
		if(!is.na(save.file))
		{
			cat(paste('\nsave to file',save.file))
			save(file= save.file, inf, coal)
		}		
	}
	#
	list(inf=inf, coal=coal)		
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
project.athena.Fisheretal.sensitivity<- function()
{
	require(data.table)
	require(ape)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_140423",sep='/')		
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"
	outdir					<- paste(DATA,"fisheretal_140423",sep='/')
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	method.risk				<- c(	'm2B1st.cas','m2Bwmx.cas','m2Bt.cas','m2Bwmx.tp1','m2Bwmx.tp2','m2Bwmx.tp3','m2Bwmx.tp4',
									'm2B1st.cas.adj','m2Bwmx.cas.adj','m2Bt.cas.adj','m2Bwmx.tp1.adj','m2Bwmx.tp2.adj','m2Bwmx.tp3.adj','m2Bwmx.tp4.adj',
									'm2Bwmx.tp1.cens','m2Bwmx.tp2.cens','m2Bwmx.tp3.cens','m2Bwmx.tp4.cens',
									'm2Bwmx.tp1.clu.cens','m2Bwmx.tp2.clu.cens','m2Bwmx.tp3.clu.cens','m2Bwmx.tp4.clu.cens',
									'm3.i','m3.ni','m3.nic','m3.nicv','m3.tni','m3.tnic','m3.tnicv','m3.tniv','m3.tnicvNo',
									'm3.i.clu','m3.ni.clu','m3.nic.clu','m3.nicv.clu','m3.tni.clu','m3.tnic.clu','m3.tnicv.clu','m3.tniv.clu','m3.tnicvNo.clu',
									'm3.nic.clu.adj','m3.nicv.clu.adj','m3.tnic.clu.adj','m3.tnicv.clu.adj','m3.tnicvNo.clu.adj',
									'm3.nic.adj','m3.nicv.adj','m3.tnic.adj','m3.tnicv.adj','m3.tnicvNo.adj')
							
	
							
	if(resume)
	{
		options(show.error.messages = FALSE)		
		file				<- paste(indir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.R", sep='')		
		readAttempt			<- try(suppressWarnings(load(file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))		
		file				<- paste(indir, '/', infile, '_', gsub('/',':',insignat), '_', "method.fit.R", sep='')
		readAttempt			<- try(suppressWarnings(load(file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))		
		options(show.error.messages = TRUE)		
	}
	if(!resume)
	{
		runs.opt				<- list()
		if(1)
		{
			run.opt							<- list()
			run.opt$method					<- '3c'
			run.opt$method.nodectime		<- 'any'
			run.opt$method.dating			<- 'gmrf'		
			run.opt$infiletree				<- paste(infile,"examlbs500",sep="_")							
			run.opt$clu.infilexml.opt		<- "mph4clutx4tip"
			run.opt$clu.infilexml.template	<- "um192rhU2080"	
			run.opt$outfile					<- paste(infile,'Ac=MY_D=35_gmrf',sep='_')
			runs.opt[[length(runs.opt)+1]]	<- run.opt
		}
		if(1)
		{
			run.opt							<- list()
			run.opt$method					<- '3d'
			run.opt$method.nodectime		<- 'any'
			run.opt$method.dating			<- 'gmrf'
			run.opt$infiletree				<- paste(infile,"examlbs500",sep="_")							
			run.opt$clu.infilexml.opt		<- "mph4clutx4tip"
			run.opt$clu.infilexml.template	<- "um192rhU2080"	
			run.opt$outfile					<- paste(infile,'Ac=MY_D=35_gmrf',sep='_')
			runs.opt[[length(runs.opt)+1]]	<- run.opt
		}
		if(1)
		{
			run.opt							<- list()
			run.opt$method					<- '3c'
			run.opt$method.nodectime		<- 'any'
			run.opt$method.dating			<- 'sasky'
			run.opt$infiletree				<- paste(infile,"examlbs500",sep="_")									
			run.opt$clu.infilexml.opt		<- "clrh80"
			run.opt$clu.infilexml.template	<- "sasky_sdr06fr"	
			run.opt$outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
			runs.opt[[length(runs.opt)+1]]	<- run.opt
		}
		if(1)
		{
			run.opt							<- list()
			run.opt$method					<- '3d'
			run.opt$method.nodectime		<- 'any'
			run.opt$method.dating			<- 'sasky'
			run.opt$infiletree				<- paste(infile,"examlbs500",sep="_")									
			run.opt$clu.infilexml.opt		<- "clrh80"
			run.opt$clu.infilexml.template	<- "sasky_sdr06fr"	
			run.opt$outfile					<- paste(infile,'Ac=MY_D=35_sasky',sep='_')
			runs.opt[[length(runs.opt)+1]]	<- run.opt
		}
		runs.opt	<- lapply( runs.opt, function(run.opt)
				{
					method			<- paste(run.opt$method, ifelse(run.opt$method.nodectime=='any','a','m'), sep='')
					files			<- list.files(indir)		
					files			<- files[ sapply(files, function(x) 	grepl(run.opt$outfile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) &
												grepl(paste('_Yscore',method,'_',sep=''), x, fixed=1) &
												grepl('R$',x) ) ]
					tmp				<- sapply(method.risk, function(x){
								tmp	<- which(grepl(paste('_',x,'.R',sep=''), files, fixed=1))
								ifelse(length(tmp)==0, NA_character_, files[tmp])
							})				
					ans				<- data.table(file=tmp, method.brl=method, method.nodectime=run.opt$method.nodectime, method.dating=run.opt$method.dating, method.risk=method.risk)				
				})
		runs.opt	<- do.call('rbind', runs.opt)
		setkey(runs.opt, method.dating, method.brl)	
		
		#	missing runs
		print( subset(runs.opt, is.na(file)) )
		runs.opt	<- subset(runs.opt, !is.na(file))
		#	load risk estimates
		runs.risk	<- runs.opt[,	{
					tmp	<- paste(indir, file, sep='/')
					cat(paste('\nprocess file=',file))
					tmp	<- load(tmp)
					ans	<- ans$risk
					if(!any(colnames(ans)=='t.period'))
						ans[, t.period:= NA_character_]										
					set(ans, NULL, 'factor', ans[, as.character(factor)])											
					ans[, method.risk:=method.risk]
					ans[, method.dating:=method.dating]
					ans[, method.nodectime:=method.nodectime]
					ans[, method.brl:=method.brl ]
					ans
				},by='file']
		runs.risk[, file:=NULL]	
		file			<- paste(indir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.R", sep='')
		save(runs.risk, file=file)
		#	load regression models
		sapply(c('m2wmx.tp','m2Bwmx.tp'), function(x) grepl, x=method.risk)
		runs.fit	<- subset(runs.opt, !grepl('m2wmx.tp', method.risk) & !grepl('m2Bwmx.tp', method.risk))		
		runs.fit	<- runs.fit[,
					{
						cat(paste('\nprocess file=',file))
						tmp	<- load( paste(indir, file, sep='/') )
						list( AIC= AIC(ans$fit.rr), logLik=logLik(ans$fit.rr), method.risk=method.risk, method.dating=method.dating, method.nodectime=method.nodectime, method.brl=method.brl )
					}, by=file]
		runs.fit[, file:=NULL]
		file			<- paste(indir, '/', infile, '_', gsub('/',':',insignat), '_', "method.fit.R", sep='')
		save(runs.fit, file=file)	
						
		#	load risk tables
		runs.table	<- runs.opt[,	{
					tmp	<- paste(indir, file, sep='/')
					cat(paste('\nprocess file=',file))
					tmp	<- load(tmp)					
					ans	<- ans$risk.table
					if(!any(colnames(ans)=='p'))
						ans	<- merge(ans, ans[, list(risk=risk, factor=factor, p= n/sum(n)  ), by='stat'], by=c('risk','factor','stat'))
					ans[, method.risk:=method.risk]
					ans[, method.dating:=method.dating]
					ans[, method.nodectime:=method.nodectime]
					ans[, method.brl:=method.brl ]
					ans
				},by='file']
		runs.table[, file:=NULL]
		set(runs.table, NULL, 'n', runs.table[, n*t.period])		
		file			<- paste(indir, '/', infile, '_', gsub('/',':',insignat), '_', "method.table.R", sep='')
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
	tmp	<- subset(runs.fit, 	grepl('m2B',method.risk,fixed=1) & grepl('cas',method.risk,fixed=1))
	tmp[, dAIC:= AIC-min(AIC)]
	setkey(tmp, method.brl, method.dating, dAIC)
	tmp
	#	check P cascade m2Bwmx.cas vs m2Bwmx.cas.adj for 3ca
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3ca' & method.dating=='sasky' & method.risk%in%c('m2Bwmx.cas','m2Bwmx.cas.adj') & stat=='P')
	#	check P cascade m2Bwmx.cas across 3ca 3da gmrf sasky
	#	3ca vs 3da 	more sensitive than gmrf vs sasky
	tmp	<- subset(runs.risk, method.nodectime=='any'  & method.risk=='m2Bwmx.cas' & stat=='P')
	#
	#	Table 2
	#
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2Bwmx.cas' & stat=='RR')
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2Bwmx.cas' & stat=='P')
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2Bwmx.cas' & stat=='RI')	
	tmp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2Bwmx.cas' & stat=='PY')	
	tmp	<- subset(runs.table, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m2Bwmx.cas' & stat=='YX')
	#
	#	MODEL 2B time trends
	#
	run.tp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m2Bwmx.tp',method.risk) & grepl('cens',method.risk) & !grepl('clu',method.risk) & stat=='P', c(risk, factor, stat, v, l95.bs, u95.bs, method.risk))
	run.tp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m2Bwmx.tp',method.risk) & grepl('adj',method.risk) & !grepl('clu',method.risk) & stat=='P', c(risk, factor, stat, v, l95.bs, u95.bs, method.risk))
	run.tp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m2Bwmx.tp',method.risk) & grepl('cens',method.risk) & grepl('clu',method.risk) & stat=='P', c(risk, factor, stat, v, l95.bs, u95.bs, method.risk))
	#run.tp	<- subset(runs.risk, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m2Bwmx.tp',method.risk) & stat=='N', c(risk, factor, stat, v, l95.bs, u95.bs, method.risk))	
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])
	run.tp[, cascade.stage:=run.tp[, substr(factor, 1, 1)]]
	set(run.tp, run.tp[,which(cascade.stage=='A')], 'cascade.stage', 'cART initiated')
	set(run.tp, run.tp[,which(cascade.stage=='U')], 'cascade.stage', 'Undiagnosed')
	set(run.tp, run.tp[,which(cascade.stage=='D')], 'cascade.stage', 'Diagnosed')
	set(run.tp, NULL, 'cascade.stage', run.tp[, factor(cascade.stage, levels=c('Undiagnosed','Diagnosed','cART initiated'))])
	ggplot(subset(run.tp, !is.na(v)), aes(x=t.period, y=v, colour=factor, group=factor)) + labs(x="calendar time periods", y="proportion of transmissions") + #scale_colour_brewer(palette="Set3") + 
				geom_point() + geom_line(aes(linetype=factor)) + geom_ribbon(aes(ymin=l95.bs, ymax=u95.bs, linetype=NA, fill=factor), alpha=0.3) + facet_grid(. ~ cascade.stage, margins=FALSE)		
	#	PYIW trends in %	
	run.tp	<- subset(runs.table, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & grepl('m2Bwmx.tp',method.risk), c(risk, factor, stat, n, prop))	
	setkey(run.tp, factor)	
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
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
	tmp		<- c(	'Undiagnosed,\nStrong evidence for acute infection',
					'Undiagnosed,\nSome evidence for acute infection',
					'Undiagnosed,\nNo evidence for acute infection',
					'Diagnosed,\nStrong evidence for acute infection\nuntil 3 months after diagnosis',
					'Diagnosed,\nSome evidence for acute infection\nuntil 3 months after diagnosis',
					'Diagnosed,\nNo evidence for acute infection,\nmin CD4 > 500',
					'Diagnosed,\nNo evidence for acute infection,\nmin CD4 > 350',
					'Diagnosed,\nNo evidence for acute infection,\nmin CD4 < 350',
					'Diagnosed,\nNo evidence for acute infection,\nmissing 1st CD4',
					'ART initiated,\nViral load missing',
					'ART initiated,\nViral load permanently\nsuppressed', 
					'ART initiated,\nViral load not permanently\nsuppressed')
	tmp		<- data.table( factor.legend= factor(tmp), factor=c("UAy","UAm","U","Dtl500","Dtl350","Dtg500","Dt.NA","DAy","DAm","ART.vlNA","ART.suA.Y","ART.suA.N"))
	run.tp	<- merge(run.tp, tmp, by='factor')	
	run.tp	<- merge( run.tp, run.tp[, list(sum=sum(n)),by=c('stat','t.period')],by=c('stat','t.period') )
	subset(run.tp, stat=='cohort' & cascade.stage=='Undiagnosed')
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
	run.tp	<- merge( run.tp, run.tp[, list(factor=factor, propa= na/sum(na, na.rm=TRUE)), by=c('t.period','stat')], by=c('t.period','stat','factor'))
	
	ggplot(subset(run.tp, !is.na(n)), aes(x=t.period, y=prop, colour=factor, group=factor)) + labs(x="calendar time periods", y="proportion PYIW") + #scale_colour_brewer(palette="Set3") + 
			geom_point() + geom_line(aes(linetype=factor)) + facet_grid(stat ~ cascade.stage, margins=FALSE)
	#	PYIW trends in #	
	ggplot(subset(run.tp, !is.na(n)), aes(x=t.period, y=n, colour=factor, group=factor)) + labs(x="calendar time periods", y="PYIW") + #scale_colour_brewer(palette="Set3") + 
			geom_point() + geom_line(aes(linetype=factor)) + facet_grid(stat ~ cascade.stage, scales='free', margins=FALSE)
	file	<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "m2wmx.tp_sasky_3da_PYIWtotal.pdf", sep='')
	ggsave(file=file, w=8,h=8)
	#set(run.tp, run.tp[, which(factor=='ART.suA.Y' & stat=='cohort' & t.period=='4')], 'na', run.tp[which(factor=='ART.suA.Y' & stat=='cohort' & t.period=='3' ), na])
	#	PYIW trends in # with adjustment for censoring
	require(RColorBrewer)
	require(grid)
	ggplot(subset(run.tp, !is.na(n)), aes(x=t.period, y=n, colour=factor.legend, group=factor.legend)) + labs(x="calendar time periods", y="PYIW") + 
			scale_colour_manual(name='Transmission times', values = colorRampPalette(brewer.pal(9, "Set1"),interpolate='spline',space="Lab")(run.tp[,length(unique(factor.legend))])) +			
			theme(legend.key.size=unit(13,'mm')) + guides(colour = guide_legend(override.aes = list(size=5))) + 
			geom_point(data=subset(run.tp, !is.na(n) & n!=na), aes(y=na), shape=8, show_guide= FALSE) + geom_point() +  geom_line(aes(linetype=factor.legend), show_guide= FALSE) + facet_grid(stat ~ cascade.stage, scales='free', margins=FALSE) 
	file	<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "m2wmx.tp_sasky_3da_PYIWtotal_adjcens.pdf", sep='')
	ggsave(file=file, w=8,h=8)
	#	PYIW trends in % with adjustment for censoring
	ggplot(subset(run.tp, !is.na(n)), aes(x=t.period, y=prop, colour=factor.legend, group=factor.legend)) + labs(x="calendar time periods", y="% PYIW") + 
			scale_colour_manual(name='Transmission times', values = colorRampPalette(brewer.pal(9, "Set1"),interpolate='spline',space="Lab")(run.tp[,length(unique(factor.legend))])) +			
			theme(legend.key.size=unit(13,'mm')) + guides(colour = guide_legend(override.aes = list(size=5))) + 			
			geom_point(data=subset(run.tp, !is.na(n) & prop!=propa), aes(y=propa), shape=8, show_guide= FALSE) + geom_point() +  geom_line(data=subset(run.tp, !is.na(n) & prop!=propa), aes(y=propa), colour='black', show_guide= FALSE) + geom_line(aes(linetype=factor.legend), show_guide= FALSE) + facet_grid(stat ~ cascade.stage, margins=FALSE)
	file	<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "m2wmx.tp_sasky_3da_PYIWprop_adjcens.pdf", sep='')
	ggsave(file=file, w=8,h=8)
	#	PYIW Margin in # before adjustment	
	ggplot(subset(run.tp, !is.na(n)), aes(x=t.period, y=n, fill=factor.legend)) + labs(x="calendar time periods", y="PYIW") +
			scale_fill_manual(name='Transmission times', values = colorRampPalette(brewer.pal(9, "Set1"),interpolate='spline',space="Lab")(run.tp[,length(unique(factor.legend))])) +			
			theme(legend.key.size=unit(13,'mm')) + guides(colour = guide_legend(override.aes = list(size=5))) + 						
			geom_bar(stat="identity") + facet_grid(stat ~ ., scales='free', margins=FALSE)
	file	<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "m2wmx.tp_sasky_3da_PYIWtotalmargin.pdf", sep='')
	ggsave(file=file, w=8,h=8)
	#	PYIW Margin in % before adjustment	
	ggplot(subset(run.tp, !is.na(n)), aes(x=t.period, y=prop, fill=factor.legend)) + labs(x="calendar time periods", y="% PYIW") +
			scale_fill_manual(name='Transmission times', values = colorRampPalette(brewer.pal(9, "Set1"),interpolate='spline',space="Lab")(run.tp[,length(unique(factor.legend))])) +			
			theme(legend.key.size=unit(13,'mm')) + guides(colour = guide_legend(override.aes = list(size=5))) + 						
			geom_bar(stat="identity") + facet_grid(stat ~ ., scales='free', margins=FALSE)
	file	<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "m2wmx.tp_sasky_3da_PYIWpropmargin.pdf", sep='')
	ggsave(file=file, w=8,h=8)
	
	
	
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
	

	#	Table 3
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnicvNo') & stat=='RR', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )	
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnicv.adj') & stat=='RR', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnic') & stat=='RR', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnic.adj') & stat=='RR', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnicvNo.adj') & stat=='RR', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	
	
	file	<- subset(runs.opt, method.nodectime=='any' & method.brl=='3da' & method.dating=='sasky' & method.risk=='m3.tnic.adj')[,file]
	tmp		<- load( paste(indir, file, sep='/') )
	tmp		<- ans$YX
	tmp[, table(stage)]
	
	
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnic.adj') & stat=='PY', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	tmp	<- subset(runs.risk,  method.risk%in%c('m3.tnic.adj') & stat=='P', c(factor, v, l95.bs, u95.bs,method.risk, method.dating, method.nodectime, method.brl ) )
	
	tmp	<- subset(runs.risk, method.risk%in%c("m3.i","m3.ni","m3.nic","m3.tnic","m3.tni") & stat=='RR' & (factor=='ART.pulse' | risk=='ART.pulse') & l95.bs<1)	
	setkey(tmp, method.brl, method.dating, method.risk)
	tmp
	
	
	subset(runs.risk,  method.risk=='m3.i' & stat=='RR' & factor=='ART.I')
	
	
	subset(runs.risk, method.risk%in%c("m3.tni") & stat=='RR'	)
	
	subset(runs.risk, method.risk%in%c("m3.tnic","m3.tniv","m3.ni","m3.tni","m3.tnicvNo") & stat=='RR'	)
	subset(runs.risk, method.risk%in%c("m3.i","m3.nic","m3.tnic","m3.tniv","m3.ni","m3.tni","m3.tnicvNo") )
	

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
project.athena.Fisheretal.YX.part1<- function(df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=NULL, df.tpairs=NULL, tperiod.info=NULL, t.period=0.25, t.endctime=2013., sample.n=NA, sample.exclude=NULL, save.file=NA, resume=FALSE)
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
		X.incare				<- project.athena.Fisheretal.X.incare(df.tpairs, df.all, df.viro, df.immu, df.treatment, t.period=t.period, t.endctime=t.endctime)
		X.incare				<- project.athena.Fisheretal.X.nocontact(X.incare, df.viro, df.immu, df.tpairs, df.all, contact.grace=0.5, t.period=t.period, t.endctime= t.endctime)		
		X.incare				<- project.athena.Fisheretal.X.CDCC(X.incare, df.tpairs, df.all, t.period=t.period, t.endctime=t.endctime)
		X.incare				<- project.athena.Fisheretal.X.followup(X.incare, df.all, df.immu, t.period=t.period, t.endctime=t.endctime)
		X.b4care				<- project.athena.Fisheretal.X.b4care(df.tpairs, df.all, predict.t2inf, t2inf.args, t.period=t.period)
		tmp						<- merge(X.incare, X.b4care, by=c('t.Patient','t'))
		cat(paste('\nnumber entries (Patient,t) that overlap in and before care [should be zero], n=',nrow(tmp)))
		X.pt					<- merge(X.incare, X.b4care, by=c('t.Patient','t'), all.x=1, all.y=1)
		set(X.pt, X.pt[,which(is.na(stage))], 'stage', 'U')					
		X.pt[, AnyPos_T1:=NULL]
		X.pt[, AnyT_T1:=NULL]
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
		if(!is.null(tperiod.info))
		{
			X.pt[, t.period:= cut(X.pt[,t], breaks=c(-Inf,tperiod.info[,t.period.min],tperiod.info[nrow(tperiod.info),t.period.max],Inf), labels=seq.int(0,nrow(tperiod.info)+1), right=FALSE)]
			X.pt				<- merge(X.pt, tperiod.info, by='t.period') 			
		}
		#	compute infection window of recipient for direct potential transmitters ri	
		Y.infwindow				<- project.athena.Fisheretal.Y.infectiontime( ri, df.all, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, score.min=0.1, score.set.value=1, method='for.infected')
		Y.infwindow[, AnyPos_a:=NULL]
		if('FASTASampleCode'%in%colnames(df.tpairs))		#mode 1
		{			
			tmp						<- subset(df.tpairs, select=c(t.Patient, t.FASTASampleCode, FASTASampleCode))
			X.pt					<- merge(X.pt, tmp, by='t.Patient', allow.cartesian=TRUE)									#merge all seq pairs
			tmp						<- subset(df.tpairs, select=c(Patient, t.FASTASampleCode, FASTASampleCode))
			Y.infwindow				<- merge(Y.infwindow, tmp, by='Patient', allow.cartesian=TRUE)							#merge all seq pairs
			YX.part1				<- merge(Y.infwindow, X.pt, by=c('FASTASampleCode','t.FASTASampleCode','t'), allow.cartesian=TRUE)
			set(YX.part1,NULL,'FASTASampleCode', YX.part1[, as.character(FASTASampleCode)])
			set(YX.part1,NULL,'t.FASTASampleCode', YX.part1[, as.character(t.FASTASampleCode)])	
			setkey(YX.part1, FASTASampleCode, t.FASTASampleCode, t)
		}			
		else if(is.na(sample.n))							#mode 2
		{			
			cat(paste('\nstart big merge of patient pairs'))
			#	the BIG merge. this is easier if we first reduce to the relevant time periods
			X.b4care				<- X.incare	<- X.ARTpulsed<- X.t2.vlsupp<- X.t2.care<- tmp<- NULL
			gc()
			X.pt					<- merge(X.pt, unique(subset( Y.infwindow, select=t )), by='t')		
			YX.part1				<- merge(Y.infwindow, X.pt, by='t', allow.cartesian=TRUE)
			cat(paste('\ncompleted big merge of patient pairs, nrows=',nrow(YX.part1)))
			setkey(YX.part1, Patient, t.Patient, t)
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
			cat(paste('\ncompleted sampled merge of sequence pairs, nrows=',nrow(YX.part1)))
			set(YX.part1,NULL,'FASTASampleCode', YX.part1[, as.character(FASTASampleCode)])
			set(YX.part1,NULL,'t.FASTASampleCode', YX.part1[, as.character(t.FASTASampleCode)])		
			setkey(YX.part1, FASTASampleCode, t.FASTASampleCode, t)
		}	
		if(!is.na(save.file))
		{			
			YX.part1	<- unique(YX.part1)
			cat(paste('\nsave YX.part1 to file=',save.file))
			save(YX.part1, file=save.file)
		}
		Y.infwindow<- X.b4care<- X.incare<- X.pt<- X.ARTpulsed<- X.t2.vlsupp<- X.t2.care<- tmp<- NULL
		gc()	
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
hivc.prog.betareg.estimaterisks<- function()
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
	t.recent.endctime		<- hivc.db.Date2numeric(as.Date("2013-03-01"))	
	#t.recent.endctime		<- hivc.db.Date2numeric(as.Date("2011-01-01"))
	t.recent.endctime		<- floor(t.recent.endctime) + floor( (t.recent.endctime%%1)*100 %/% (t.period*100) ) * t.period
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
		outfile					<- paste(infile,'_Ac=MY_D=35_gmrf',sep='')
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
		outfile					<- paste(infile,'_Ac=MY_D=35_gmrf',sep='')
	}	
	if(1)
	{		
		method					<- '3d'
		method.nodectime		<- 'any'
		method.risk				<- 'm2B1st.cas'
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
		infiletree				<- paste(infile,"examlbs500",sep="_")
		insignat				<- "Wed_Dec_18_11:37:00_2013"							
		clu.infilexml.opt		<- "clrh80"
		clu.infilexml.template	<- "sasky_sdr06fr"	
		outfile					<- paste(infile,'_Ac=MY_D=35_sasky',sep='')
	}
	if(0)
	{		
		method					<- '3c'
		method.nodectime		<- 'any'
		method.risk				<- 'm21st.cas'
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
	}	
	clu.infile		<- infile
	clu.indir		<- indir
	clu.insignat	<- insignat	
	outfile			<- paste( outfile, ifelse(t.recent.endctime==t.endctime,'',paste('_',t.recent.endctime,sep='')), sep='')
	
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
		files		<- files[ sapply(files, function(x) grepl(outfile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(paste('Yscore',method,sep=''), x, fixed=1) & grepl(paste(method.risk,'.R',sep=''),x, fixed=1)  ) ]		
		stopifnot(length(files)==0)		
	}
	#
	#	check if we have precomputed tables
	#
	X.tables			<- NULL
	if(0)
	{
		save.file		<- NA
		if(grepl('m21st',method.risk))		save.file	<- 'm21st'
		if(grepl('m2B1st',method.risk))		save.file	<- 'm2B1st'
		if(grepl('m2t',method.risk))		save.file	<- 'm2t'
		if(grepl('m2Bt',method.risk))		save.file	<- 'm2Bt'
		if(grepl('m2wmx',method.risk))		save.file	<- 'm2wmx'
		if(grepl('m2Bwmx',method.risk))		save.file	<- 'm2Bwmx'
		if(grepl('m3.i',method.risk) & !grepl('m3.ni',method.risk))								save.file	<- 'm3.i'	
		if(grepl('m3',method.risk) & grepl('ni',method.risk) & !grepl('nic',method.risk))		save.file	<- 'm3.ni'
		if(grepl('m3',method.risk) & grepl('tni',method.risk) & !grepl('tnic',method.risk))		save.file	<- 'm3.tni'
		if(grepl('m3',method.risk) & grepl('nic',method.risk) & !grepl('tnic',method.risk))		save.file	<- 'm3.nic'
		if(grepl('m3',method.risk) & grepl('tnic',method.risk) & !grepl('No',method.risk))		save.file	<- 'm3.tnic'
		if(grepl('m3',method.risk) & grepl('tnic',method.risk) & grepl('No',method.risk))		save.file	<- 'm3.tnicNo'
		if(is.na(save.file))	stop('unknown method.risk')				
		tmp				<- regmatches(method.risk, regexpr('tp[0-9]', method.risk))
		save.file		<- paste(save.file, ifelse(length(tmp), paste('.',tmp,sep=''), ''), sep='')		
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_tables_',save.file,'.R',sep='')
		X.tables		<- project.athena.Fisheretal.estimate.risk.table(YX=NULL, X.den=NULL, X.msm=NULL, X.clu=NULL, resume=TRUE, save.file=save.file, method=method.risk)
		if(!is.null(X.tables))	cat('\nloaded X.tables')
	}
	#
	#	get rough idea about (backward) time to infection from time to diagnosis, taking midpoint of SC interval as 'training data'
	#
	tmp				<- project.athena.Fisheretal.t2inf(indircov, infile.cov.study, adjust.AcuteByNegT=0.75, adjust.dt.CD4=1, adjust.AnyPos_y=2003, adjust.NegT=2)
	predict.t2inf	<- tmp$predict.t2inf
	t2inf.args		<- tmp$t2inf.args
	#
	#	get data relating to full population (MSM including those without seq)
	#
	if(1)
	{
		tmp					<- project.athena.Fisheretal.select.denominator(indir, infile, insignat, indircov, infile.cov.all, infile.viro.all, infile.immu.all, infile.treatment.all, infiletree=NULL, adjust.AcuteByNegT=adjust.AcuteByNegT, adjust.NegT4Acute=NA, adjust.NegTByDetectability=0.25, adjust.minSCwindow=0.25, adjust.AcuteSelect=c('Yes','Maybe'), t.recent.endctime=t.recent.endctime)	
		df.all.allmsm		<- tmp$df.all	
		df.viro.allmsm		<- tmp$df.viro
		df.immu.allmsm		<- tmp$df.immu
		df.treatment.allmsm	<- tmp$df.treatment
		tmp					<- tmp$df.select
		setkey(tmp, Patient)
		ri.allmsm			<- unique(tmp)	
	}
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
	df.tpairs		<- project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos(clumsm.info, df.denom, any.pos.grace.yr= 3.5, select.if.transmitter.seq.unique=FALSE)	
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
		tmp				<- merge( subset(df.tpairs, select=Patient), subset(df.denom, select=c(Patient, AnyPos_T1)), by='Patient' )
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
	#
	#	get timelines for the candidate transmitters in ATHENA.clu to the recently infected RI.PT; remove zero scores
	#
	resume			<- 1
	rm.zero.score	<- TRUE
	save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RICT_',method,'_tATHENAclu','.R',sep='')	
	YX.part1		<- project.athena.Fisheretal.YX.part1(df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=NULL, df.tpairs=df.tpairs, tperiod.info=NULL, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)
	YX.part1		<- merge( YX.part1, subset( df.tpairs, select=c(FASTASampleCode, t.FASTASampleCode, cluster) ), by=c('FASTASampleCode','t.FASTASampleCode'), all.x=1)
	YX.part1[, class:='pt']
	gc()
	save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'YX',method,'.R',sep='')	
	YX				<- project.athena.Fisheretal.YX.part2(YX.part1, df.all, predict.t2inf, t2inf.args, indir, insignat, indircov, infile.cov.study, infiletree, outdir, outfile, cluphy=cluphy, cluphy.info=cluphy.info, cluphy.map.nodectime=cluphy.map.nodectime, df.tpairs.4.rawbrl=df.tpairs, rm.zero.score=rm.zero.score, t.period=t.period, save.file=save.file, resume=resume, method=method)
	gc()
	tperiod.info	<- merge(df.all, unique( subset(YX, select=c(Patient, t.period)) ), by='Patient')
	tperiod.info	<- tperiod.info[, list(t.period.min=min(AnyPos_T1)), by='t.period']
	tperiod.info[, t.period.max:=c(tperiod.info[-1, t.period.min], t.endctime)]
	ri.PT			<- subset(YX, select=c(Patient, t))[, list(n.t.infw= length(unique(t))), by='Patient']
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
		YXS.part1		<- project.athena.Fisheretal.YX.part1(df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=tmp, df.tpairs=NULL, tperiod.info=NULL, sample.n=sample.n, sample.exclude=sample.exclude, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)
		if( !all(YX.part1[,FASTASampleCode]%in%df.all[,FASTASampleCode]) | !all(YX.part1[,t.FASTASampleCode]%in%df.all[,FASTASampleCode]) | !all(YXS.part1[,FASTASampleCode]%in%df.all[,FASTASampleCode]) | !all(YXS.part1[,t.FASTASampleCode]%in%df.all[,FASTASampleCode]) ) stop('unexpected FASTASampleCode')		
		YX.part1[, class:='pt']
		YXS.part1[, class:='pn']
		YXS.part1		<- do.call('rbind',list(YX.part1, subset(YXS.part1, select=colnames(YX.part1))))		#put potential transmitters and potential non-transmitters together		
		YXS.part1		<- merge( YXS.part1, subset( df.tpairs, select=c(FASTASampleCode, t.FASTASampleCode, cluster) ), by=c('FASTASampleCode','t.FASTASampleCode'), all.x=1)	
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'YX',method,'_RICT_tATHENAseq','.R',sep='')
		YXS				<- project.athena.Fisheretal.YX.part2(YXS.part1, df.all, predict.t2inf, t2inf.args, indir, insignat, indircov, infile.cov.study, infiletree, outdir, paste(outfile, 'RICT_tATHENAseq', sep='_'), cluphy=cluphy, cluphy.info=cluphy.info, cluphy.map.nodectime=cluphy.map.nodectime, df.tpairs.4.rawbrl=df.tpairs, rm.zero.score=rm.zero.score, t.period=t.period, save.file=save.file, resume=resume, method=method)		
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
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RIPDT_',method,'_tATHENAclu','.R',sep='')
		X.clu			<- project.athena.Fisheretal.YX.part1(tmp, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=ri.PT, df.tpairs=NULL, tperiod.info=tperiod.info, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)
		gc()		
	}
	#
	#	get timelines for all candidate transmitters in df.all (anyone with subtype B sequ) to the recently infected RI.PT
	#	
	X.seq			<- NULL
	if(is.null(X.tables))
	{
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RIPDT_',method,'_tATHENAseq','.R',sep='')
		X.seq			<- project.athena.Fisheretal.YX.part1(df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, ri=ri.PT, df.tpairs=NULL, tperiod.info=tperiod.info, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)	
		gc()				
	}
	#
	#	get timelines for all candidate transmitters in ATHENA.MSM (anyone with MSM exposure group irrespective of sequ available or not) to the recently infected RI.PT
	#	
	X.msm				<- NULL	
	if(is.null(X.tables) & ( grepl('adj',method.risk) | grepl('cens',method.risk)))
	{
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RIPDT_',method,'_tATHENAmsm','.R',sep='')
		X.msm			<- project.athena.Fisheretal.YX.part1(df.all.allmsm, df.immu.allmsm, df.viro.allmsm, df.treatment.allmsm, predict.t2inf, t2inf.args, ri=ri.PT, df.tpairs=NULL, tperiod.info=tperiod.info, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)
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
		#unique(subset(ri.allmsm, select=Patient)) 
		ri.allmsm.su			<- project.athena.Fisheretal.RI.summarize(unique(subset(ri.allmsm, select=Patient)), df.all.allmsm, tperiod.info, info=NULL, with.seq=FALSE, with.cluster=FALSE, cols=c('lRNA_T1','CD4_T1','AnyPos_A','SC_dt','PosCD4_dt','PoslRNA_dt'))		
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
		plot.file.varyvl		<- NA	#paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_model2_',method.risk,'_VL_adjAym_dt025','.pdf',sep='')
		if(any(sapply(c('m2wmx','m2Bwmx'), grepl, x=method.risk)))
			YX					<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(YX, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=plot.file.varyvl, plot.file.or=NA )			
		if(any(sapply(c('m21st','m2B1st'), grepl, x=method.risk)))
			YX					<- project.athena.Fisheretal.YX.model2.stratify.VL1stsu(YX, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=plot.file.varyvl)
		if(any(sapply(c('m2t','m2Bt'), grepl, x=method.risk)))
			YX					<- project.athena.Fisheretal.YX.model2.stratify.VLt(YX, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=plot.file.varyvl)
		if(substr(method.risk,1,2)=='m3')			
			YX					<- project.athena.Fisheretal.YX.model3.stratify.ARTriskgroups(YX, df.all, return.only.ART=TRUE)					
	}
	if(is.null(X.tables))
	{		
		cat(paste('\nstratify by', method))
		if(substr(method.risk,1,2)=='m2')
		{
			plot.file.varyvl	<- NA	#paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_model2_',method.risk,'_VL_adjAym_dt025','.pdf',sep='')
			if(any(sapply(c('m2wmx','m2Bwmx'), grepl, x=method.risk)))
			{
				YX					<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(YX, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=plot.file.varyvl, plot.file.or=NA )
				X.clu				<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.clu, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA, plot.file.or=NA )
				X.seq				<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.seq, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA, plot.file.or=NA )
				X.msm				<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(X.msm, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, vl.suppressed=log10(1e3), plot.file.varyvl=NA, plot.file.or=NA )
			}
			if(any(sapply(c('m21st','m2B1st'), grepl, x=method.risk)))
			{
				YX					<- project.athena.Fisheretal.YX.model2.stratify.VL1stsu(YX, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=plot.file.varyvl)
				X.clu				<- project.athena.Fisheretal.YX.model2.stratify.VL1stsu(X.clu, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA)
				X.seq				<- project.athena.Fisheretal.YX.model2.stratify.VL1stsu(X.seq, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA)
				X.msm				<- project.athena.Fisheretal.YX.model2.stratify.VL1stsu(X.msm, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, vl.suppressed=log10(1e3), plot.file.varyvl=NA)
			}
			if(any(sapply(c('m2t','m2Bt'), grepl, x=method.risk)))
			{
				YX					<- project.athena.Fisheretal.YX.model2.stratify.VLt(YX, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=plot.file.varyvl)
				X.clu				<- project.athena.Fisheretal.YX.model2.stratify.VLt(X.clu, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA)
				X.seq				<- project.athena.Fisheretal.YX.model2.stratify.VLt(X.seq, df.all, df.viro, df.immu, vl.suppressed=log10(1e3), plot.file.varyvl=NA)
				X.msm				<- project.athena.Fisheretal.YX.model2.stratify.VLt(X.msm, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, vl.suppressed=log10(1e3), plot.file.varyvl=NA)
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
		#	compute tables
		if(grepl('adj',method.risk) & grepl('clu',method.risk))
		{
			save.file		<- NA
			if(grepl('m21st',method.risk))		save.file	<- 'm21st'
			if(grepl('m2B1st',method.risk))		save.file	<- 'm2B1st'
			if(grepl('m2t',method.risk))		save.file	<- 'm2t'
			if(grepl('m2Bt',method.risk))		save.file	<- 'm2Bt'
			if(grepl('m2wmx',method.risk))		save.file	<- 'm2wmx'
			if(grepl('m2Bwmx',method.risk))		save.file	<- 'm2Bwmx'
			if(grepl('m3.i',method.risk) & !grepl('m3.ni',method.risk))								save.file	<- 'm3.i'	
			if(grepl('m3',method.risk) & grepl('ni',method.risk) & !grepl('nic',method.risk))		save.file	<- 'm3.ni'
			if(grepl('m3',method.risk) & grepl('tni',method.risk) & !grepl('tnic',method.risk))		save.file	<- 'm3.tni'
			if(grepl('m3',method.risk) & grepl('nic',method.risk) & !grepl('tnic',method.risk))		save.file	<- 'm3.nic'
			if(grepl('m3',method.risk) & grepl('tnic',method.risk) & !grepl('No',method.risk))		save.file	<- 'm3.tnic'
			if(grepl('m3',method.risk) & grepl('tnic',method.risk) & grepl('No',method.risk))		save.file	<- 'm3.tnicNo'
			if(is.na(save.file))	stop('unknown method.risk')				
			tmp				<- regmatches(method.risk, regexpr('tp[0-9]', method.risk))
			save.file		<- paste(save.file, ifelse(length(tmp), paste('.',tmp,sep=''), ''), sep='')		
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_tables_',save.file,'.R',sep='')
			X.tables		<- project.athena.Fisheretal.estimate.risk.table(YX, X.seq, X.msm, X.clu, resume=TRUE, save.file=save.file, method=method.risk)
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
	if(grepl('m2',method.risk))
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_model2_',method.risk,'.R',sep='')
	if(grepl('m3',method.risk))
		save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_model2_',method.risk,'.R',sep='')	
	tmp					<- project.athena.Fisheretal.estimate.risk.wrap(YX, X.tables, plot.file.or=NA, bs.n=1e3, resume=resume, save.file=save.file, method=method.risk)					
	
	stop()
		
		
		
	#
	#	get data for selection
	#	
	tmp						<- project.athena.Fisheretal.get.dated.phylo.for.selection(df.tpairs, clu.indir, clu.infile, clu.insignat, clu.infilexml.opt, clu.infilexml.template, method.nodectime=method.nodectime)
	#df.all2					<- subset(clumsm.info, select=c(Patient, cluster))
	#setkey(df.all2, Patient)
	#df.all2					<- merge( tmp$df.all, unique(df.all2), by='Patient', all.x=TRUE )
	df.tpairs				<- tmp$df.tpairs
	cluphy					<- tmp$clu$cluphy
	cluphy.subtrees			<- tmp$clu$cluphy.subtrees
	cluphy.info				<- tmp$clu$cluphy.info
	cluphy.map.nodectime	<- tmp$clu$cluphy.map.nodectime	
	#			
	#	plot	
	#	
	outfile					<- paste(indir,'/',infile, '_', clu.infilexml.template, '_', clu.infilexml.opt, '_', gsub('/',':',insignat), '_', 'pt_anypos_3.5_anynodectime', '.pdf',sep='')	
	project.athena.Fisheretal.plot.selected.transmitters(clumsm.info, df.immu, df.viro, df.treatment, df.tpairs, cluphy, cluphy.info, cluphy.subtrees, cluphy.map.nodectime, outfile, pdf.height=600)

	
}
