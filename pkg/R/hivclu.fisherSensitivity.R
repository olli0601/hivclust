######################################################################################
project.athena.Fisheretal.sensitivity<- function()
{
	require(data.table)
	require(ape)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_141108",sep='/')
	outdir					<- paste(DATA,"fisheretal_141108",sep='/')		
	indir					<- paste(DATA,"fisheretal_141221",sep='/')
	outdir					<- paste(DATA,"fisheretal_141221",sep='/')		
	indir					<- paste(DATA,"fisheretal_150105",sep='/')
	outdir					<- paste(DATA,"fisheretal_150105",sep='/')		
	indir					<- paste(DATA,"fisheretal_150216",sep='/')
	outdir					<- paste(DATA,"fisheretal_150216",sep='/')		
	indir					<- paste(DATA,"fisheretal_150303",sep='/')
	outdir					<- paste(DATA,"fisheretal_150303",sep='/')		
	indir					<- paste(DATA,"fisheretal_150308",sep='/')
	outdir					<- paste(DATA,"fisheretal_150308",sep='/')		
	
	
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
		runs.opt	<- subset(runs.opt, !grepl('beforepool',file) & !grepl('Hypo',file))
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
		#	load risk.bs estimates
		tmp			<- lapply(seq_len(nrow(runs.opt)), function(i)
				{
					tmp	<- paste(indir, runs.opt[i,file], sep='/')
					cat(paste('\nprocess file=',runs.opt[i,file]))
					tmp	<- load(tmp)
					ans	<- ans$risk.bs
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
		runs.riskbs	<- do.call('rbind', tmp)
		#	get pooled proportions across all tp
		tmp			<- lapply( runs.risk[, unique(method.brl)], function(METHOD.BRL)
				{
					df		<- subset( runs.risk, !grepl('ARTstarted|GroupsUDA|Hypo', method.risk) & method.brl==METHOD.BRL)
					df.bs	<- subset( runs.riskbs, !grepl('ARTstarted|GroupsUDA|Hypo', method.risk) & method.brl==METHOD.BRL)
					project.athena.Fisheretal.sensitivity.pool.TPALL(df, df.bs)			
				})
		tmp			<- do.call('rbind', tmp)
		tmp[, bs:=NULL]
		runs.risk	<- rbind( runs.risk, tmp, use.names=TRUE )
		file			<- paste(indir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.Rdata", sep='')
		save(runs.risk, file=file)
		#	reduce runs.opt to files for which we have a table
		runs.opt	<- subset(runs.opt, !grepl('ARTstarted', method.risk) & !grepl('GroupsUDA', method.risk) )
		#	load risk tables
		tmp			<- lapply(seq_len(nrow(runs.opt)), function(i)
				{
					tmp	<- paste(indir, runs.opt[i,file], sep='/')
					cat(paste('\nprocess file=',runs.opt[i,file]))
					tmp	<- load(tmp)					
					ans	<- ans$X.tables$risk.table
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
project.athena.Fisheretal.sensitivity.getfigures<- function()
{	
	require(data.table)
	require(ape)
	require(grid)
	require(reshape2)
	require(ggplot2)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_141221",sep='/')
	outdir					<- paste(DATA,"fisheretal_141221",sep='/')		
	indir					<- paste(DATA,"fisheretal_150105",sep='/')
	outdir					<- paste(DATA,"fisheretal_150105",sep='/')		
	indir					<- paste(DATA,"fisheretal_150216",sep='/')
	outdir					<- paste(DATA,"fisheretal_150216",sep='/')		
	indir					<- paste(DATA,"fisheretal_150303",sep='/')
	outdir					<- paste(DATA,"fisheretal_150303",sep='/')		
	indir					<- paste(DATA,"fisheretal_150308",sep='/')
	outdir					<- paste(DATA,"fisheretal_150308",sep='/')		
	
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
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.5, 2006.5, 2008, 2009.5), t.period.max = c(2006.45, 2007.99, 2009.45, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )	
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-7','2008-1','2009-7'), labels=c('96/07','\n\n06/07','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-6','2007-12','2009-6','2010-12'), labels=c('06/06','07/12','09/06','10/12'))])
	
	#	updated stages
	#	set up factor legends
	factor.color	<- c(	"#990000","#EF6548","#FDBB84",
			"#0C2C84","#0570B0","#74A9CF","#41B6C4","#35978F",  
			"#FCC5C0","#F768A1","#7A0177",
			"#1A9850","#A6D96A","grey70"
			)
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
			'ART initiated,\n After first viral suppression\nViral suppression, >1 observations','Not in contact'
			)
	levels			<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.NotYetFirstSu","ART.vlNA","ART.suA.N","ART.suA.Y1","ART.suA.Y2","Lost")
	#levels			<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.NotYetFirstSu","ART.vlNA","ART.suA.N","ART.suA.Y1","ART.suA.Y2")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(levels, levels=levels), factor.color=factor.color, method.risk='m2Cwmx')	
	factors			<- subset(factors, grepl('m2Cwmx',method.risk), select=c(factor, factor.legend, factor.color))
	#
	#	WTN 3ma 	m2Cwmx
	#
	#	WTN 3pa1 1.35 T7 BRL INF
	method.DENOM	<- 'SEQ'	
	
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'	
	stat.select		<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile			<- infile
	method.BRLs		<- c('3pa1H1.48C2V100bInfT7', '3pa1H1.94C2V100bInfT7', '3pa1H1.09C2V100bInfT7','3pa1H1.48C1V100bInfT7','3pa1H1.48C3V100bInfT7')
	method.BRLs		<- c('3pa1H1.48C2V100b0.02T7','3pa1H1.48C2V100b0.04T7')
	dummy			<- sapply(method.BRLs, function(method.BRL)
			{				
				method.RISK		<- 'm2Cwmx.wtn.tp'
				project.athena.Fisheretal.sensitivity.getfigures.m2(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  factors, stat.select, outfile, tperiod.info=tperiod.info)		
				method.RISK		<- "m2Cwmx.wtn"
				project.athena.Fisheretal.sensitivity.getfigures.RR(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING,  factors, stat.select, outfile, tperiod.info=tperiod.info)
				project.athena.Fisheretal.sensitivity.tables.m2.prop(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING, factors, stat.select, outfile, tperiod.info=tperiod.info)				
			})
	
	
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
project.athena.Fisheretal.sensitivity.getfigures.comparetransmissionintervals.notime<- function(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, group)
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
	tmp		<- tmp[, list(n.coh=sum(n)), by='factor']
	tmp[, p.coh:= tmp[, n.coh/sum(n.coh)]]
	set(tmp, NULL, 'n.coh', tmp[, round(n.coh, d=0)])			
	ans		<- subset(tmp, select=c(factor,  p.coh))
	#	prop likely transmission intervals 
	tmp		<- subset(df, stat%in%c('YX','nRecLkl'))
	tmp2	<- subset(df, stat=='YX', select=c(factor, t.period, p))
	tmp		<- dcast.data.table(tmp, t.period+factor+method.risk~stat, value.var='n')
	tmp		<- merge(tmp, tmp2, by=c('factor','t.period'))
	tmp[, n:= YX/nRecLkl*100]	
	tmp		<- tmp[, list(n.lkl=sum(n)), by='factor']
	tmp[, p.lkl:= tmp[, n.lkl/sum(n.lkl)]]	
	set(tmp, NULL, 'n.lkl', tmp[, round(n.lkl, d=0)])		
	ans		<- merge(ans, subset(tmp, select=c(factor, p.lkl)), by='factor')
	#	prop likely transmission intervals + expected missing intervals per 100 'recipient with lkl transmitter' in time period
	tmp		<- subset(df, stat%in%c('YX','Sx.e0cp','nRecLkl'))
	tmp2	<- subset(df, stat=='Sx.e0cp', select=c(factor, t.period, l95.bs, u95.bs))
	tmp		<- dcast.data.table(tmp, t.period+factor+method.risk~stat, value.var='n')
	tmp		<- merge(tmp, tmp2, by=c('factor','t.period'))
	tmp[, n:= (YX+Sx.e0cp/8)/nRecLkl*100]
	tmp[, l95.bs:= (YX+l95.bs/8)/nRecLkl*100]
	tmp[, u95.bs:= (YX+u95.bs/8)/nRecLkl*100]	
	tmp		<- tmp[, list(n.lklm=sum(n), l95.lklm=sum(l95.bs), u95.lklm=sum(u95.bs)), by='factor']
	tmp[, p.lklm:= tmp[, n.lklm/sum(n.lklm)]]
	ans		<- merge(ans, subset(tmp, select=c(factor, p.lklm)), by='factor')
	#	ratio of proportions
	ans[, r.lkl:= p.lkl/p.coh]
	ans[, r.lklm:= p.lklm/p.coh]
	ans		<- melt(ans, measure.vars=c('p.coh','p.lklm'))
	# group
	if(group=='cascade')
	{
		scale.name	<- 'in infection/care stage'
		ans[, group:=ans[, substr(factor, 1, 1)]]
		set(ans, ans[,which(group=='A')], 'group', 'ART started')
		set(ans, ans[,which(group=='U')], 'group', 'Undiagnosed')
		set(ans, ans[,which(group=='D')], 'group', 'Diagnosed')
		set(ans, ans[,which(group=='L')], 'group', 'Not in contact')
		set(ans, NULL, 'group', ans[, factor(group, levels=c('Undiagnosed','Diagnosed','ART started','Not in contact'))])		
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
	#factor.long	
	ans	<- merge(ans, factors, by='factor')	
	ggplot(ans, aes(x=interaction(variable, factor), y=value*100, colour=variable, fill=factor.legend)) + geom_bar(stat='identity', position='dodge') +
			scale_fill_manual(values=ans[, unique(factor.color)], guide=FALSE) +
			scale_colour_manual(values=c('transparent','black'), guide=FALSE) +
			scale_y_continuous(breaks=seq(0,50,10), minor_breaks=seq(0,50,2)) +
			theme_bw() + 
			geom_text(aes(x=seq_along(interaction(variable, factor)), y=-2, label='*'), size=7) +
			geom_text(aes(x=14, y=26, label='* phylogenetically probable transmitters'), size=3) +
			labs(x='', y='Transmission intervals\n(%)') +
			theme(axis.text=element_text(size=14), axis.title=element_text(size=14), legend.key.size=unit(11,'mm'), legend.position = "bottom", legend.box = "vertical", axis.ticks.x=element_blank(), axis.text.x=element_blank(), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major.x=element_blank(), panel.grid.major=element_line(colour="grey70", size=0.7)) 
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','enrich',".pdf", sep='')	
	ggsave(file=file, w=5, h=3)
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
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-5','2008-1','2009-7'), labels=c('96/07','\n\n06/05','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-4','2007-12','2009-6','2010-12'), labels=c('06/04','07/12','09/06','10/12'))])	
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
	indir			<- paste(DATA,"fisheretal",sep='/')
	outdir			<- paste(DATA,"fisheretal_141108",sep='/')	
	infile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXm2_SEQ_3pa1H1C3V100.R'
	outfile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXm2Cwmx.wtn_SEQ_3pa1H1C3V100_lkl.pdf'	
	file			<- paste(indir, '/', infile, sep='')
	load(file)
	YX				<- copy(YX.m2)
	set(YX, NULL, 'score.Y', YX[, score.Y*w.tn])
	YX				<- merge(YX, YX[, list(tprob=score.Y/sum(score.Y), t.Patient=t.Patient, t=t), by=c('Patient')], by=c('Patient','t.Patient','t'))
	
	YX				<- subset(YX, select=c(t.period, t, t.Patient, Patient, score.Y, tprob, telapsed, CD4b, CD4c ))
	setnames( YX, 'CD4c', 'factor')	
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
	#
	#
	ans[, summary(tprob)]
	#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
	#0.000000 0.004108 0.023080 0.048170 0.070660 1.000000 
	ans[, summary(telapsed)]
	#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
	#0.4337  0.8674  1.5510  1.9490  2.6480  9.0890 
	#
	#	boxplot of likelihood
	#
	ggplot(ans, aes(x=factor.legend, y=score.Y, colour=factor.legend, fill=factor.legend), log='y')  + 
			labs(x='', y='phylogenetic, relative probability\nof direct HIV transmission') +				
			scale_colour_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
			scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)]) +
			scale_y_continuous(breaks=c(0.1, 0.5, 1, 2, 5, 10, 20, 100), minor_breaks=c(0.1)) +
			scale_x_discrete(breaks=NULL) +
			coord_trans(ytrans="log10", limy=c(0.01, 1000)) +		
			geom_boxplot(outlier.shape = NA) +
			stat_summary(geom= "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +		
			theme_bw() + theme(legend.key.size=unit(11,'mm'), legend.position = "bottom", legend.box = "vertical", legend.title.align = 0, axis.ticks.margin = unit(c(0,0,0,0), "lines"), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major.x=element_blank(), panel.grid.major=element_line(colour="grey70", size=0.7)) +
			guides(fill=guide_legend(nrow=4))
	file			<- paste(outdir, '/', outfile, sep='')
	ggsave(file=file, w=8, h=6)		
	#
	#	boxplot of transmission prob (observed only)
	#
	ggplot(ans, aes(x=factor.legend, y=tprob, colour=factor.legend, fill=factor.legend))  + 
			labs(x='', y='phylogenetic probability\nof direct HIV transmission') +
			geom_boxplot(outlier.shape = NA) +
			scale_colour_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
			scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
			scale_y_continuous(breaks=seq(0, 0.5, 0.1), minor_breaks=seq(0, 0.5, 0.02)) +
			scale_x_discrete(breaks=NULL) +	
			#coord_trans(ytrans="log10", limy=c(0.001, 0.3)) +
			coord_cartesian(ylim=c(0, 0.31)) +
			stat_summary(geom= "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +		
			theme_bw() + theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.key.size=unit(11,'mm'), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major.x=element_blank(), panel.grid.major=element_line(colour="grey70", size=0.7))
	outfile			<- paste(substr(outfile, 1, nchar(outfile)-7),'tprob.pdf',sep='')
	file			<- paste(outdir, '/', outfile, sep='')
	ggsave(file=file, w=6, h=4)		
	
	
	#
	#	CD4c updated
	#
	indir			<- paste(DATA,"fisheretal",sep='/')
	outdir			<- paste(DATA,"fisheretal_150216",sep='/')
	infile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXm2_SEQ_3pa1H1.35C3V100bInfT7.R'
	outfile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXm2_SEQ_3pa1H1.35C3V100bInfT7_lkl.pdf'
	file			<- paste(indir, '/', infile, sep='')
	load(file)
	YX				<- copy(YX.m2)
	set(YX, NULL, 'score.Y', YX[, score.Y*w.tn])
	YX				<- merge(YX, YX[, list(tprob=score.Y/sum(score.Y), t.Patient=t.Patient, t=t), by=c('Patient')], by=c('Patient','t.Patient','t'))	
	YX				<- subset(YX, select=c(t.period, t, t.Patient, Patient, score.Y, tprob, telapsed, CD4b, CD4c ))
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
	ggplot(ans, aes(x=factor.legend, y=tprob, colour=factor.legend, fill=factor.legend))  + 
			labs(x='', y='phylogenetic probability\nof direct HIV transmission') +
			geom_boxplot(outlier.shape = NA) +
			scale_colour_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
			scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
			scale_y_continuous(breaks=seq(0, 0.5, 0.1), minor_breaks=seq(0, 0.5, 0.02)) +
			scale_x_discrete(breaks=NULL) +	
			#coord_trans(ytrans="log10", limy=c(0.001, 0.3)) +
			coord_cartesian(ylim=c(0, 0.31)) +
			stat_summary(geom= "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +		
			theme_bw() + theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.key.size=unit(11,'mm'), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major.x=element_blank(), panel.grid.major=element_line(colour="grey70", size=0.7))
	outfile			<- paste(substr(outfile, 1, nchar(outfile)-7),'tprob.pdf',sep='')
	file			<- paste(outdir, '/', outfile, sep='')
	ggsave(file=file, w=7, h=3)		
	
	
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
	set(ans, ans[,which(group=='A')], 'group', 'After diagnosis')
	set(ans, ans[,which(group=='U')], 'group', 'Before diagnosis')
	set(ans, ans[,which(group=='D')], 'group', 'After diagnosis')
	set(ans, ans[,which(group=='L')], 'group', 'After diagnosis')
	
	tmp		<- subset(ans, stat=='expected missing', select=c(t.period, factor, stat, group, l95.bs, u95.bs))
	tmp		<- tmp[, list(stat=stat[1], factor=factor[1], l95.bs=sum(l95.bs), u95.bs=sum(u95.bs)), by=c('t.period','group')]	
	ans		<- merge( 	subset(ans, select=c(t.period, factor, group, stat, n)), 
			subset(tmp, select=c(t.period, factor, stat, l95.bs, u95.bs)), by=c('t.period','factor','stat'), all.x=1) 
	
	set(ans, NULL, 'group', ans[, paste(group,'\n(',stat,')',sep='')]) 
	set(ans, NULL, 'stat', NULL)		
	set(ans, NULL, 'group',	ans[, factor(group, levels=c('Before diagnosis\n(observed)','After diagnosis\n(observed)','Before diagnosis\n(expected missing)','After diagnosis\n(expected missing)'))])				
	
	#tperiod.long
	ans	<- merge(ans, tperiod.info, by='t.period')
	ans[, t.period.long:= paste(t.period.min, '-', t.period.max,sep='')]
	set(ans, NULL, 't.period.long', ans[,factor(t.period.long, levels= tperiod.info[, paste(t.period.min, '-', t.period.max,sep='')])])	
	#factor.long	
	ans	<- merge(ans, factors, by='factor')	
	
	ggplot(ans, aes(x=t.period.long, y=n, ymin=l95.bs, ymax=u95.bs, group=group, fill=factor.legend))  + 
			labs(x='', y='probable transmission intervals\nper 100 recipient MSM\n(person-years)') +				
			scale_fill_manual(name=scale.name, values=ans[, unique(factor.color)], guide=FALSE) +
			scale_y_continuous(breaks=seq(0,1000,100), minor_breaks=seq(0,1000,20)) +
			geom_bar(stat='identity') + geom_errorbar(stat='identity',width=.2) +  				
			theme_bw() + theme(strip.text= element_text(size=14), axis.text.x=element_text(size=14),  axis.text.y=element_text(size=14), axis.title.y=element_text(size=14), legend.key.size=unit(12,'mm'), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major.x=element_blank(), panel.grid.major=element_line(colour="grey70", size=0.7)) +
			facet_grid(. ~ group, scales='free_y', margins=FALSE)	
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','lkladj',".pdf", sep='')	
	ggsave(file=file, w=12, h=5)	
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
	set(ans, ans[,which(group=='L')], 'group', 'Not in contact')
	
	
	tmp		<- subset(ans, group=='Undiagnosed\n(observed)', select=c(t.period, factor, X.msm.e0cp, l95.bs, u95.bs))
	setnames(tmp, 'X.msm.e0cp', 'n')
	tmp[, group:='Undiagnosed\n(adjusted for\nright censoring)']
	tmp2	<- tmp[, list(factor='U', l95.bs=sum(l95.bs), u95.bs=sum(u95.bs)), by=c('t.period','group')]
	tmp		<- merge(subset(tmp, select=c(t.period, factor, n, group)), subset(tmp2, select=c(t.period, factor, l95.bs, u95.bs)), by=c('t.period','factor'), all.x=1) 
	
	ans	<- subset(ans, select=c(t.period, factor, n, l95.bs, u95.bs, group))
	set(ans, NULL, c('l95.bs', 'u95.bs'), NA_real_)
	ans	<- rbind(ans, tmp, use.names=TRUE)
	set(ans, NULL, 'group', ans[, factor(group, levels=c('Undiagnosed\n(observed)','Undiagnosed\n(adjusted for\nright censoring)','Diagnosed','ART started','Not in contact'))])
	
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
				scale_y_continuous(breaks=seq(0,20000,1000), minor_breaks=seq(0,20000,250)) +
				geom_bar(stat='identity') + geom_errorbar(stat='identity',width=.2) +  				
				theme_bw() + theme(legend.key.size=unit(12,'mm'), axis.title.x = element_text(vjust=-5), panel.grid.major.x = element_blank(), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.7)) +				
				facet_grid(. ~ group, scales='free_y', margins=FALSE)	
		file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','ptadj',".pdf", sep='')	
		ggsave(file=file, w=12, h=5)
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
				theme_bw() + theme(legend.key.height=unit(12,'mm'), legend.key.width=unit(0,'mm'), panel.grid.major.x = element_blank(), legend.key=element_rect(colour='transparent', fill='transparent'), axis.text.x=element_text(angle = -60, vjust = 0.5, hjust=0.5)) +
				facet_grid(. ~ group, scales='free_y', margins=FALSE)	
		file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_','2011','_',method.DENOM, '_',method.BRL,'_',df[1, method.risk],'_','ptadj',".pdf", sep='')	
		ggsave(file=file, w=9, h=5)	
	}
}
######################################################################################
project.athena.Fisheretal.sensitivity.getfigures.pseq<- function(runs.table, file, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, tperiod.info)
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
	#tperiod.long
	ans	<- merge(ans, tperiod.info, by='t.period')
	ans[, t.period.long:= gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))]
	tmp		<- c("96/07-06/06", "06/07-07/12", "08/01-09/06", "09/07-10/12")
	set(ans, NULL, 't.period.long', ans[, factor(t.period.long, levels=tmp, labels=tmp)])	
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
		ggsave(file=file, w=25, h=8)		
	}
	if(1)
	{
		ggplot(ans, aes(x=p.seq, y=factor.legend, pch=t.period.long, colour=factor.legend)) + geom_point(size=3) +
				scale_x_continuous(breaks=seq(0,100,10), limit=c(0,100), expand=c(0,0)) +
				scale_colour_manual(values=ans[, unique(factor.color)], guide=FALSE) +
				theme_bw() +				
				theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=7), axis.title=element_text(size=10), legend.position='bottom', panel.grid.major.x=element_line(colour="grey70", size=0.4), panel.grid.minor.x=element_line(colour="grey70", size=0.4), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) + 				
				labs(y='', x='potential transmission intervals\nof a potential transmitter with a sequence\n(%)',pch='time of diagnosis\nof recipient MSM') +
				guides(pch=guide_legend(ncol=2))
		cat(paste('plot to',file))
		ggsave(file=file, w=7, h=7)
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
	if(0)
	{
		run.tp			<- subset(runs.risk, (is.na(t.period) | t.period!=0) & method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating==method.DATING & grepl(method.RISK,method.risk)  & (grepl('P.',stat,fixed=1) | stat=='P') )
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
		run.tp[, t.period.long:= gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))]
		set(run.tp, NULL, 't.period.long', run.tp[,factor(t.period.long, levels= tperiod.info[, gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))])])	
		dummy	<- lapply(seq_along(stat.select), function(i)
				{
					cat(paste('\nprocess', stat.select[i]))
					ggplot(subset(run.tp, !is.na(v) & stat==stat.select[i]), aes(x=factor, y=v, fill=factor.legend, colour=factor.legend)) + labs(x="", y=ylab) + 
							scale_y_continuous(breaks=seq(0,0.9,0.1), minor_breaks=seq(0,0.9,0.02), labels=paste(seq(0,0.9,0.1)*100,'%',sep='')) + 
							scale_x_discrete(breaks=NULL, limits=run.tp[, intersect(levels(factor), as.character(factor))]) +
							scale_fill_manual(name='from cascade stage', values=run.tp[, unique(factor.color)], guide=with.guide) + 
							scale_colour_manual(name='from cascade stage', values = rep('black',run.tp[, length(unique(factor))]), guide=with.guide) +
							#scale_fill_brewer(palette='PRGn',name='from cascade stage') + scale_colour_manual(name='from cascade stage', values = rep('black',11)) +					
							#guides(colour=FALSE, fill = guide_legend(override.aes = list(size=5))) +
							theme_bw() + theme(legend.key.size=unit(10.5,'mm'), plot.margin=unit(c(0,0,-5,0),"mm"), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.7)) + #coord_flip() +
							geom_bar(stat='identity',binwidth=1, position='dodge') + geom_errorbar(aes(ymin=l95.bs, ymax=u95.bs), width=0.3, position=position_dodge(width=0.9))	+ 
							facet_grid(. ~ t.period.long, margins=FALSE)
					file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_',run.tp[1, method.recentctime],'_',run.tp[1, method.denom], '_',run.tp[1, method.dating], '_',run.tp[1, method.brl],'_',run.tp[1, method.risk],'_',stat.select[i],".pdf", sep='')
					cat(paste('\nsave to file',file))							
					ggsave(file=file, w=12,h=4)							
				})	
	}
	#
	#
	#
	run.tp			<- subset(runs.risk, (is.na(t.period) | t.period!=0) & method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating==method.DATING & grepl(method.RISK,method.risk)  & (grepl('P.',stat,fixed=1) | stat=='P' | grepl('N.',stat,fixed=1) | stat=='N') )
	run.tp			<- subset(run.tp, !grepl('ARTstarted',method.risk) & !grepl('GroupsUDA',method.risk))
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
	ylab		<- "Proportion of transmissions\n(Number of transmissions)"	
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
	run.tp	<- merge(run.tp, subset(factors, select=factors[, which(colnames(factors)!='method.risk')]), by='factor')	
	run.tp	<- merge(run.tp, tperiod.info, by='t.period')
	run.tp[, t.period.long:= gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))]
	set(run.tp, NULL, 't.period.long', run.tp[,factor(t.period.long, levels= tperiod.info[, gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))])])
	
	yscale	<- subset(run.tp, factor=='UA')
	set(yscale, NULL, 'stat', yscale[, as.character(stat)])
	yscale[, stat1:= yscale[, substr(stat,1,1)]]
	yscale[, stat:= yscale[, substr(stat,3,nchar(stat))]]
	yscale	<- yscale[, list(scale= v[which(stat1=='N')]/v[which(stat1=='P')]), by=c('t.period.long','stat')]
	set(yscale, NULL, 'stat', yscale[, paste('P.',stat,sep='')])
	run.tp	<- merge(run.tp, yscale, by=c('t.period.long','stat'))
	dummy	<- lapply(seq_along(stat.select), function(i)
			{
				cat(paste('\nprocess', stat.select[i]))
				tmp	<- subset(run.tp, !is.na(v) & stat==stat.select[i])
				ggplot(tmp, aes(x=factor, y=v, fill=factor.legend, colour=factor.legend)) + labs(x="", y=ylab) + 
						scale_y_continuous(breaks=seq(0,0.9,0.1), minor_breaks=seq(0,0.9,0.02), labels=paste(seq(0,0.9,0.1)*100,'%\n(n=',round(seq(0,0.9,0.1)*tmp$scale[1], d=1),')',sep='')) + 
						scale_x_discrete(breaks=NULL, limits=run.tp[, intersect(levels(factor), as.character(factor))]) +
						scale_fill_manual(name='from cascade stage', values=run.tp[, unique(factor.color)], guide=with.guide) + 
						scale_colour_manual(name='from cascade stage', values = rep('black',run.tp[, length(unique(factor))]), guide=with.guide) +
						#scale_fill_brewer(palette='PRGn',name='from cascade stage') + scale_colour_manual(name='from cascade stage', values = rep('black',11)) +					
						#guides(colour=FALSE, fill = guide_legend(override.aes = list(size=5))) +
						theme_bw() + theme(legend.key.size=unit(10.5,'mm'), plot.margin=unit(c(0,0,-5,0),"mm"), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.7)) + #coord_flip() +
						geom_bar(stat='identity',binwidth=1, position='dodge') + geom_errorbar(aes(ymin=l95.bs, ymax=u95.bs), width=0.3, position=position_dodge(width=0.9))	+ 
						facet_grid(. ~ t.period.long, margins=FALSE)
				file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_',run.tp[1, method.recentctime],'_',run.tp[1, method.denom], '_',run.tp[1, method.dating], '_',run.tp[1, method.brl],'_',run.tp[1, method.risk],'_',stat.select[i],".pdf", sep='')
				cat(paste('\nsave to file',file))							
				ggsave(file=file, w=12,h=4)							
			})
	
	#
	ylab		<- "Relative transmissibility"
	run.tp		<- subset(runs.risk, (is.na(t.period) | t.period!=0) &  method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating==method.DATING & grepl(method.RISK,method.risk)  & (grepl('RI.',stat,fixed=1) | stat=='RI') )
	run.tp		<- subset(run.tp, !grepl('ARTstarted',method.risk) & !grepl('GroupsUDA',method.risk))
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
	run.tp[, t.period.long:= paste(t.period.min, '-', t.period.max,sep='')]
	set(run.tp, NULL, 't.period.long', run.tp[,factor(t.period.long, levels= tperiod.info[, paste(t.period.min, '-', t.period.max,sep='')])])		
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
							theme_bw() + theme(legend.key.size=unit(10.5,'mm'), panel.grid.minor=element_line(colour="grey70", size=0.2), panel.grid.major=element_line(colour="grey70")) + #coord_flip() +
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
							theme_bw() + theme(legend.key.size=unit(11,'mm'), panel.grid.minor=element_line(colour="grey70", size=0.2), panel.grid.major=element_line(colour="grey70")) +
							facet_wrap( ~ group, scales='free_y', ncol=4)		
					ggsave(file=file, w=12,h=4)
				}									
			})
	#
	ylab		<- "Risk ratio relative to\nDiagnosed, CD4 progression to > 500"
	run.tp		<- subset(runs.risk, (is.na(t.period) | t.period!=0) &  method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating==method.DATING & grepl(method.RISK,method.risk)  & (grepl('RR.',stat,fixed=1) | stat=='RR') & grepl('Dtg500',factor.ref))
	run.tp		<- subset(run.tp, !grepl('ARTstarted',method.risk) & !grepl('GroupsUDA',method.risk))
	stat.select	<- gsub('RI','RR', stat.select)
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
	run.tp[, t.period.long:= paste(t.period.min, '-', t.period.max,sep='')]
	set(run.tp, NULL, 't.period.long', run.tp[,factor(t.period.long, levels= tperiod.info[, paste(t.period.min, '-', t.period.max,sep='')])])		
	color	<- run.tp[, unique(factor.color)]
	dummy	<- lapply(seq_along(stat.select), function(i)
			{
				cat(paste('\nprocess', stat.select[i]))
				file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_',run.tp[1, method.recentctime],'_',run.tp[1, method.denom], '_',run.tp[1, method.dating], '_',run.tp[1, method.brl],'_',run.tp[1, method.risk],'_',stat.select[i],".pdf", sep='')
				cat(paste('\nsave to file',file))				
				tmp		<- subset(run.tp, !is.na(v) & stat==stat.select[i])
				ggplot(tmp, aes(x=t.period.long, y=v, ymin=l95.bs, ymax=u95.bs, fill=factor.legend, colour=factor.legend, group=factor.legend)) + 
						labs(x="", y=ylab) + #coord_trans(ytrans=my_trans, limy=c(0.01, 15)) +	
						#scale_y_continuous(breaks=c(seq(0,1,0.2), seq(1,4,1), seq(4,20,2))) +			
						scale_fill_manual(name='from cascade stage', values = color, guide=with.guide) + 
						scale_colour_manual(name='from cascade stage', values = color, guide=with.guide) +
						geom_hline(yintercept=1) + geom_errorbar(position=position_dodge(.5), width=0.4) + geom_point(position=position_dodge(.5) , shape=18, size=3) +  
						theme_bw() + theme(legend.key.size=unit(11,'mm'), panel.grid.minor=element_line(colour="grey70", size=0.2), panel.grid.major=element_line(colour="grey70")) +
						facet_wrap( ~ group, scales='free_y', ncol=4)
				ggsave(file=file, w=12,h=4)
			})
	
	
	#
	if(0)
	{
		ylab		<- "Number of transmissions"
		run.tp		<- subset(runs.risk, (is.na(t.period) | t.period!=0) &  method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating==method.DATING & grepl(method.RISK,method.risk)  & (grepl('N.',stat,fixed=1) | stat=='N' ) )
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
		run.tp[, t.period.long:= gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))]
		set(run.tp, NULL, 't.period.long', run.tp[,factor(t.period.long, levels= tperiod.info[, gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))])])		
		dummy	<- lapply(seq_along(stat.select), function(i)
				{
					cat(paste('\nprocess', stat.select[i]))
					ggplot(subset(run.tp, !is.na(v) & stat%in%c(stat.select[i])), aes(x=factor, y=v, fill=factor.legend, colour=factor.legend)) + labs(x="", y=ylab) + 
							scale_y_continuous(breaks=seq(0,300,10), minor_breaks=seq(0,300,2)) + scale_x_discrete(breaks=NULL, limits=run.tp[, levels(factor)]) +
							scale_fill_manual(name='from cascade stage', values=run.tp[, unique(factor.color)], guide=with.guide) + 
							scale_colour_manual(name='from cascade stage', values = rep('black',run.tp[, length(unique(factor))]), guide=with.guide) +
							#scale_fill_brewer(palette='PRGn',name='from cascade stage') + scale_colour_manual(name='from cascade stage', values = rep('black',11)) +					
							#guides(colour=FALSE, fill = guide_legend(override.aes = list(size=5))) +
							theme_bw() + theme(legend.key.size=unit(10.5,'mm'), plot.margin=unit(c(0,0,-5,0),"mm"), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.7)) + #coord_flip() +
							geom_bar(stat='identity',binwidth=1, position='dodge')	+ geom_errorbar(aes(ymin=l95.bs, ymax=u95.bs), width=0.3, position=position_dodge(width=0.9))	+ 
							facet_grid(. ~ t.period.long, margins=FALSE)
					file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_',run.tp[1, method.recentctime],'_',run.tp[1, method.denom], '_',run.tp[1, method.dating], '_',run.tp[1, method.brl],'_',run.tp[1, method.risk],'_',stat.select[i],".pdf", sep='')
					cat(paste('\nsave to file',file))							
					ggsave(file=file, w=12,h=4)							
				})		
	}	
}
######################################################################################
project.athena.Fisheretal.sensitivity.getfigures.RR<- function(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING, factors, stat.select, outfile, tperiod.info=NULL, with.guide=FALSE)
{
	#
	ylab		<- "Relative transmissibility"
	run.tp		<- subset(runs.risk, 	(is.na(t.period) | t.period!=0) &  method.denom==method.DENOM & method.nodectime=='any' & method.BRL==method.brl & method.dating==method.DATING & grepl(method.RISK,method.risk)  & (grepl('RI.',stat,fixed=1) | stat=='RI') 
					& grepl('ART',factor)	)
	run.tp		<- subset( run.tp, !grepl('GroupsUDA',method.risk) )											
	stat.select	<- gsub('P','RI', stat.select)
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])
	run.tp[, group:=run.tp[, substr(factor, 1, 1)]]	
	set(run.tp, run.tp[,which(group=='A' & grepl('suA\\.Y',factor))], 'group', 'ART initiated\nViral load < 100 cps/ml')
	set(run.tp, run.tp[,which(group=='A' & !grepl('suA\\.Y',factor))], 'group', 'ART initiated\nOther')
	set(run.tp, run.tp[,which(grepl('ART.started',factor))], 'group', 'ART initiated\n( Overall )')	
	set(run.tp, NULL, 'group', run.tp[, factor(group, levels=c('ART initiated\n( Overall )','ART initiated\nViral load < 100 cps/ml','ART initiated\nOther'))])	
	set(run.tp, NULL, 'factor', run.tp[, factor(factor, levels=factors[, levels(factor)])])	
	set(run.tp, run.tp[, which(grepl('ARTstarted',method.risk))], 'factor', 'ART.started')	
	tmp		<- rbind(factors, data.table(factor='ART.started', factor.legend='ART initiated', factor.color='black'))
	run.tp	<- merge(run.tp, tmp, by='factor')		
	run.tp	<- merge(run.tp, tperiod.info, by='t.period')
	run.tp[, t.period.long:= paste(t.period.min, '-', t.period.max,sep='')]
	set(run.tp, NULL, 't.period.long', run.tp[,factor(t.period.long, levels= tperiod.info[, paste(t.period.min, '-', t.period.max,sep='')])])		
	color	<- run.tp[, unique(factor.color)]
	dummy	<- lapply(seq_along(stat.select), function(i)
			{
				cat(paste('\nprocess', stat.select[i]))
				file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_',run.tp[1, method.recentctime],'_',run.tp[1, method.denom], '_',run.tp[1, method.dating], '_',run.tp[1, method.brl],'_',run.tp[1, method.risk],'_',stat.select[i],"ARTstarted.pdf", sep='')
				cat(paste('\nsave to file',file))				
				ggplot(subset(run.tp, !is.na(v) & stat==stat.select[i]), aes(x=t.period.long, y=v, ymin=l95.bs, ymax=u95.bs, fill=factor.legend, colour=factor.legend, group=factor.legend)) + 
						geom_hline(yintercept=1) + geom_point(position=position_dodge(.5) , shape=18, size=3) + geom_errorbar(position=position_dodge(.5), width=0.4) +
						scale_fill_manual(name='from cascade stage', values = color, guide=FALSE) + 
						scale_colour_manual(name='from cascade stage', values = color, guide=FALSE) +
						theme_bw() + theme(legend.key.size=unit(11,'mm'), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_line(colour="grey70", size=0.2), panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey70", size=0.4)) +
						labs(x="", y=ylab) + facet_grid( ~ group)
				ggsave(file=file, w=8,h=4)													
			})
	#
	ylab		<- "Transmission risk ratio\nof stages after ART start\nversus Diagnosed, CD4>500"
	ylab		<- "Transmission risk ratio"
	run.tp		<- subset(runs.risk, (is.na(t.period) | t.period!=0) &  method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating==method.DATING & grepl(method.RISK,method.risk)  & (grepl('RR.',stat,fixed=1) | stat=='RR') & grepl('Dtg500',factor.ref)
					& grepl('ART',factor)	)
	run.tp		<- subset( run.tp, !grepl('GroupsUDA',method.risk) )						
	stat.select	<- gsub('RI','RR', stat.select)
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])
	run.tp[, group:=run.tp[, substr(factor, 1, 1)]]	
	set(run.tp, run.tp[,which(group=='A' & grepl('suA\\.Y',factor))], 'group', 'ART initiated\nViral load < 100 cps/ml')
	set(run.tp, run.tp[,which(group=='A' & !grepl('suA\\.Y',factor))], 'group', 'ART initiated\nOther')
	set(run.tp, run.tp[,which(grepl('ART.started',factor))], 'group', 'ART initiated\n( Overall )')	
	set(run.tp, NULL, 'group', run.tp[, factor(group, levels=c('ART initiated\n( Overall )','ART initiated\nViral load < 100 cps/ml','ART initiated\nOther'))])	
	set(run.tp, NULL, 'factor', run.tp[, factor(factor, levels=factors[, levels(factor)])])	
	set(run.tp, run.tp[, which(grepl('ARTstarted',method.risk))], 'factor', 'ART.started')	
	tmp		<- rbind(factors, data.table(factor='ART.started', factor.legend='ART initiated', factor.color='black'))
	run.tp	<- merge(run.tp, tmp, by='factor')		
	run.tp	<- merge(run.tp, tperiod.info, by='t.period')
	run.tp[, t.period.long:= paste(t.period.min, '-', t.period.max,sep='')]
	set(run.tp, NULL, 't.period.long', run.tp[,factor(t.period.long, levels= tperiod.info[, paste(t.period.min, '-', t.period.max,sep='')])])		
	color	<- run.tp[, unique(factor.color)]
	if(1)
	{
		dummy	<- lapply(seq_along(stat.select), function(i)
				{
					cat(paste('\nprocess', stat.select[i]))
					file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_',run.tp[1, method.recentctime],'_',run.tp[1, method.denom], '_',run.tp[1, method.dating], '_',run.tp[1, method.brl],'_',run.tp[1, method.risk],'_',stat.select[i],"ARTstarted.pdf", sep='')
					cat(paste('\nsave to file',file))				
					tmp		<- subset(run.tp, !is.na(v) & stat==stat.select[i] & factor.legend!='ART initiated')
					ggplot(tmp, aes(x=t.period.long, y=v, ymin=l95.bs, ymax=u95.bs, fill=factor.legend, colour=factor.legend, group=factor.legend)) + 
							geom_hline(yintercept=1) + geom_point(position=position_dodge(.5) , shape=18, size=4.5) + geom_errorbar(position=position_dodge(.5), width=0.4, size=1) +
							scale_x_discrete(expand=c(0.03,0.03)) +
							scale_y_continuous(breaks=seq(0,2,0.5), minor_breaks=seq(0,2,0.1)) + coord_trans(limy=c(0, 1.1))  +
							scale_fill_manual(name='from cascade stage', values = color, guide=FALSE) + 
							scale_colour_manual(name='from cascade stage', values = color, guide=FALSE) +
							theme_bw() + theme(axis.text=element_text(size=15), axis.title=element_text(size=18), legend.key.size=unit(11,'mm'), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_line(colour="grey70", size=0.4), panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey70", size=0.4)) +
							labs(x="", y=ylab) 
					ggsave(file=file, w=5,h=4)
				})	
	}
	if(0)
	{
		dummy	<- lapply(seq_along(stat.select), function(i)
				{
					cat(paste('\nprocess', stat.select[i]))
					file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_',run.tp[1, method.recentctime],'_',run.tp[1, method.denom], '_',run.tp[1, method.dating], '_',run.tp[1, method.brl],'_',run.tp[1, method.risk],'_',stat.select[i],"ARTstarted.pdf", sep='')
					cat(paste('\nsave to file',file))				
					tmp		<- subset(run.tp, !is.na(v) & stat==stat.select[i])
					ggplot(subset(run.tp, !is.na(v) & stat==stat.select[i]), aes(x=t.period.long, y=v, ymin=l95.bs, ymax=u95.bs, fill=factor.legend, colour=factor.legend, group=factor.legend)) + 
							geom_hline(yintercept=1) + geom_point(position=position_dodge(.5) , shape=18, size=3) + geom_errorbar(position=position_dodge(.5), width=0.4) +
							scale_y_continuous(breaks=seq(0,2,0.5), minor_breaks=seq(0,2,0.1)) + coord_trans(limy=c(0, 1.1))  +
							scale_fill_manual(name='from cascade stage', values = color, guide=FALSE) + 
							scale_colour_manual(name='from cascade stage', values = color, guide=FALSE) +
							theme_bw() + theme(legend.key.size=unit(11,'mm'), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_line(colour="grey70", size=0.4), panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey70", size=0.4)) +
							labs(x="", y=ylab) + facet_grid( ~ group)
					ggsave(file=file, w=7,h=4)
				})	
	}
	
}
######################################################################################
project.athena.Fisheretal.sensitivity.gettables<- function()
{
	require(data.table)
	require(ape)
	require(RColorBrewer)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_141108",sep='/')	
	outdir					<- paste(DATA,"fisheretal_141108",sep='/')	
	indir					<- paste(DATA,"fisheretal_150216",sep='/')	
	outdir					<- paste(DATA,"fisheretal_150216",sep='/')	
	indir					<- paste(DATA,"fisheretal_150308",sep='/')	
	outdir					<- paste(DATA,"fisheretal_150308",sep='/')	
	
	
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.table.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))	
	
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.5, 2006.5, 2008, 2009.5), t.period.max = c(2006.45, 2007.99, 2009.45, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )	
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-7','2008-1','2009-7'), labels=c('96/07','\n\n06/07','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-6','2007-12','2009-6','2010-12'), labels=c('06/06','07/12','09/06','10/12'))])
	#	
	#	updated stages
	#	set up factor legends
	factor.color	<- c(	"#990000","#EF6548","#FDBB84",
			"#0C2C84","#0570B0","#74A9CF","#41B6C4","#35978F",  
			"#FCC5C0","#F768A1","#7A0177",
			"#1A9850","#A6D96A","grey70")
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
			'ART initiated,\n After first viral suppression\n Viral suppression, >1 observations','Not in contact')
	tmp				<- c("UA","U","UAna","DA","Dtg500","Dtl500","Dtl350","Dt.NA","ART.NotYetFirstSu","ART.vlNA","ART.suA.N","ART.suA.Y1","ART.suA.Y2","Lost")
	factors			<- data.table( factor.legend= factor(factor.long, levels=factor.long), factor=factor(tmp, levels=tmp), factor.color=factor.color, method.risk='m2Cwmx')	
	factors			<- subset(factors, grepl('m2Cwmx',method.risk), select=c(factor, factor.legend, factor.color))
	runs.table		<- subset(runs.table, !grepl('ARTstarted|GroupsUDA',method.risk))
	#
	outfile			<- infile	
	method.DENOM	<- 'SEQ'	
	method.RISK		<- 'm2Cwmx.wtn.tp'
	method.WEIGHT	<- ''			
	method.BRLs		<- c('3pa1H1.48C2V100bInfT7','3pa1H1.09C2V100bInfT7','3pa1H1.94C2V100bInfT7')
	dummy			<- lapply(method.BRLs, function(method.BRL)
			{
				project.athena.Fisheretal.sensitivity.getfigures.npotentialtransmissionintervals.adjusted(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, 'cascade')
				project.athena.Fisheretal.sensitivity.getfigures.nlikelytransmissionintervals(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, 'cascade')
				project.athena.Fisheretal.sensitivity.getfigures.comparetransmissionintervals.notime(runs.table, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, outfile, tperiod.info, 'cascade')				
			})
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
project.athena.Fisheretal.sensitivity.pool.TPALL<- function(df, df.bs)
{
	df		<- subset(df, grepl('^N\\.',stat))		
	df.bs	<- subset(df.bs, grepl('^N\\.',stat))
	df[, bs:=0]
	set(df, NULL, c("l95.bs","u95.bs","m50.bs"), NULL)		
	df		<- rbind(df, df.bs, use.names=TRUE)
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	set(df, NULL, 'method.risk', df[, substr(method.risk, 1, nchar(method.risk)-4)])
	set(df, NULL, 'coef', df[, paste(risk,factor,sep='')])
	method.risk			<- df[1, method.risk]
	method.dating		<- df[1, method.dating]
	method.nodectime	<- df[1, method.nodectime]
	method.brl			<- df[1, method.brl]
	method.denom		<- df[1, method.denom]
	method.recentctime	<- df[1, method.recentctime]
	#	pool N's
	ans					<- dcast.data.table(df, coef+coef.ref+risk+risk.ref+factor+factor.ref+bs~stat, value.var='v', fun.aggregate=sum)
	#	get P's
	ans					<- melt(ans, measure.vars=c('N.raw','N.raw.e0','N.raw.e0cp'), value.name='v', variable.name='stat')
	tmp					<- ans[, list(factor=factor, v=v/sum(v)), by=c('stat','bs','risk')]
	ans					<- merge( unique(subset(ans, select=setdiff(names(ans),c('stat','bs','v')))), tmp, by=c('risk','factor') )
	set(ans, NULL, 'stat', ans[,gsub('N\\.','P\\.',stat)])
	#	UA/(U+UA)
	tmp					<- subset(df, grepl('^U$|^UA$', factor))
	tmp					<- dcast.data.table(tmp, risk+risk.ref+stat+bs ~ factor, value.var='v', fun.aggregate=sum)	
	tmp[, v:=UA/(U+UA)]
	tmp[, factor:= paste('UA', sep='')]
	tmp[, coef:= paste(risk,factor,sep='')]
	tmp[, factor.ref:= "None"]
	tmp[, coef.ref:= "None"]
	set(tmp, NULL, 'stat', tmp[, gsub('N.','CUA.',stat)])
	ans					<- rbind(ans, subset(tmp, select=intersect(names(ans),names(tmp))))
	#	groups of factors
	set(df, NULL, 'factor', df[, paste(substr(factor, 1, 1), '_total', sep='')])
	set(df, NULL, 'coef', df[, paste(risk,factor,sep='')])
	df					<- dcast.data.table(df, coef+coef.ref+risk+risk.ref+factor+factor.ref+bs~stat, value.var='v', fun.aggregate=sum)
	df					<- melt(df, measure.vars=c('N.raw','N.raw.e0','N.raw.e0cp'), value.name='v', variable.name='stat')
	tmp					<- df[, list(factor=factor, v=v/sum(v)), by=c('stat','bs','risk')]
	df					<- merge( unique(subset(df, select=setdiff(names(df),c('stat','bs','v')))), tmp, by=c('risk','factor') )
	set(df, NULL, 'stat', df[,gsub('N\\.','P\\.',stat)])
	ans					<- rbind(ans, subset(df, select=intersect(names(ans),names(df))))
	#	return only 95% bootstrap quantiles + central estimate
	tmp					<- subset(ans, bs>0)[,list(l95.bs=quantile(v,p=0.025), u95.bs=quantile(v,p=0.975), m50.bs=quantile(v,p=0.5)), by=c('risk','factor','stat')]		
	ans					<- merge( subset(ans, bs==0), tmp, by=c('risk','factor','stat') )
	ans[, method.risk:=method.risk]
	ans[, method.dating:=method.dating]
	ans[, method.nodectime:=method.nodectime]
	ans[, method.brl:=method.brl ]
	ans[, method.denom:=method.denom]
	ans[, method.recentctime:=method.recentctime ]
	ans[, t.period:=0]
	set(ans, NULL, 'factor', ans[, paste(factor, '.tp0', sep='')])
	set(ans, NULL, 'coef', ans[, paste(risk, factor, sep='')])
	set(ans, NULL, 'method.risk', ans[, paste(method.risk, '.tp0', sep='')])
	ans
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
project.athena.Fisheretal.sensitivity.tables.m2.prop<- function(runs.risk, method.DENOM, method.BRL, method.RISK, method.WEIGHT, method.DATING, factors, stat.select, outfile, tperiod.info=NULL, with.guide=FALSE)
{	
	df		<- subset(runs.risk, method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating==method.DATING & grepl(method.RISK,method.risk))  
	if(method.WEIGHT=='')
		df	<- subset(df, !grepl('wstar',method.risk) & !grepl('now',method.risk))
	if(method.WEIGHT!='')
		df	<- subset(df, grepl(method.WEIGHT,method.risk) )
	dfg		<- subset(df, grepl('GroupsUDA',method.risk))
	dfa		<- subset(df, t.period==0)
	df		<- subset(df, !grepl('ARTstarted',method.risk) & !grepl('GroupsUDA',method.risk) & (is.na(t.period) | t.period!=0))
	#	time period	
	set(df, NULL, 'factor', df[, as.character(factor)])
	setkey(df, factor)
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])	
	set(dfg, NULL, 'factor', dfg[, as.character(factor)])
	setkey(dfg, factor)
	dfg[, t.period:=dfg[, substr(factor, nchar(factor), nchar(factor))]]
	set(dfg, NULL, 'factor', dfg[, substr(factor, 1, nchar(factor)-2)])	
	#	cascade group
	set(df, NULL, 'group', df[, paste(substr(factor, 1, 1),'total',sep='_')])
	set(dfg, NULL, 'group', dfg[, paste(substr(factor, 1, 1),'total',sep='_')])
	#	reset factor order
	set(df, NULL, 'factor', df[, factor(factor, levels=levels, labels=levels)])
	setkey(df, stat, t.period, factor)
	#set(dfg, NULL, 'factor', dfg[, factor(factor, levels=levels, labels=levels)])
	setkey(dfg, stat, t.period, factor)	
	#df		<- unique(df)
	#	prop infections by period
	tmp		<- subset(df, stat%in%c('P.raw.e0cp'))	
	tmp2	<- dcast.data.table(tmp, group+factor~t.period, value.var='v')
	setnames(tmp2, as.character(1:4), paste('v.',1:4,sep=''))
	tmp2	<- merge(tmp2, dcast.data.table(tmp, group+factor~t.period, value.var='l95.bs'), by=c('group','factor'))
	setnames(tmp2, as.character(1:4), paste('l95.bs.',1:4,sep=''))
	tmp2	<- merge(tmp2, dcast.data.table(tmp, group+factor~t.period, value.var='u95.bs'), by=c('group','factor'))
	setnames(tmp2, as.character(1:4), paste('u95.bs.',1:4,sep=''))	
	#	add group prop
	tmp		<- subset(dfg, stat%in%c('P.raw.e0cp'))	
	tmp3	<- dcast.data.table(tmp, group+factor~t.period, value.var='v')
	setnames(tmp3, as.character(1:4), paste('v.',1:4,sep=''))
	tmp3	<- merge(tmp3, dcast.data.table(tmp, group+factor~t.period, value.var='l95.bs'), by=c('group','factor'))
	setnames(tmp3, as.character(1:4), paste('l95.bs.',1:4,sep=''))
	tmp3	<- merge(tmp3, dcast.data.table(tmp, group+factor~t.period, value.var='u95.bs'), by=c('group','factor'))
	setnames(tmp3, as.character(1:4), paste('u95.bs.',1:4,sep=''))	
	tmp3[, factor:= tmp3[, paste(factor,'_total',sep='')]]
	tmp2	<- rbind(tmp3, tmp2, use.names=TRUE)
	tmp2[, stat:='P.raw.e0cp']
	#	add CUA
	tmp		<- unique(subset(df, stat=='CUA.raw.e0cp'))
	tmp3	<- dcast.data.table(tmp, group+factor+stat~t.period, value.var='v')
	setnames(tmp3, as.character(1:4), paste('v.',1:4,sep=''))
	tmp3	<- merge(tmp3, dcast.data.table(tmp, group+factor~t.period, value.var='l95.bs'), by=c('group','factor'))
	setnames(tmp3, as.character(1:4), paste('l95.bs.',1:4,sep=''))
	tmp3	<- merge(tmp3, dcast.data.table(tmp, group+factor~t.period, value.var='u95.bs'), by=c('group','factor'))
	setnames(tmp3, as.character(1:4), paste('u95.bs.',1:4,sep=''))	
	tmp3[, factor:='UA']
	tmp2	<- rbind(tmp3, tmp2, use.names=TRUE)
	#	prop infections overall
	tmp		<- subset(dfa, stat%in%c('P.raw.e0cp','CUA.raw.e0cp'))
	set( tmp, NULL, 'factor', tmp[, substr(factor,1,nchar(factor)-4)] )
	set( tmp, NULL, 'group', tmp[, paste(substr(factor,1,1),'_total',sep='')])
	tmp2	<- merge( subset(tmp, select=c(group, factor, stat, v, l95.bs, u95.bs)), tmp2, by=c('group','factor','stat'))
	#	prepare number of recipient with at least one prob transmitter
	tmp3	<- unique(subset(df, stat=='nRecLkl', select=c(t.period, v)))
	tmp3[, All:=sum(v)]
	set(tmp3, NULL, 'v', tmp3[, paste('=\"(n=', v, ')\"', sep='')])
	set(tmp3, NULL, 'All', tmp3[, paste('=\"(n=', All, ')\"', sep='')])
	tmp3	<- dcast.data.table(tmp3, All~t.period, value.var='v')
	set(tmp3, NULL, c('group','factor','stat'), '')
	#	prepare Excel
	set(tmp2, NULL, 'All', tmp2[, paste('=\"',round(v*100,d=1),' (',round(l95.bs*100,d=1),'-',round(u95.bs*100,d=1),')\"',sep='')])
	set(tmp2, NULL, '1', tmp2[, paste('=\"',round(v.1*100,d=1),' (',round(l95.bs.1*100,d=1),'-',round(u95.bs.1*100,d=1),')\"',sep='')])
	set(tmp2, NULL, '2', tmp2[, paste('=\"',round(v.2*100,d=1),' (',round(l95.bs.2*100,d=1),'-',round(u95.bs.2*100,d=1),')\"',sep='')])
	set(tmp2, NULL, '3', tmp2[, paste('=\"',round(v.3*100,d=1),' (',round(l95.bs.3*100,d=1),'-',round(u95.bs.3*100,d=1),')\"',sep='')])
	set(tmp2, NULL, '4', tmp2[, paste('=\"',round(v.4*100,d=1),' (',round(l95.bs.4*100,d=1),'-',round(u95.bs.4*100,d=1),')\"',sep='')])
	set(tmp2, NULL, paste('v.',1:4,sep=''), NULL)
	set(tmp2, NULL, paste('l95.bs.',1:4,sep=''), NULL)
	set(tmp2, NULL, paste('u95.bs.',1:4,sep=''), NULL)
	set(tmp2, NULL, c('v','l95.bs','u95.bs'), NULL)
	ans		<- copy(tmp2)	
	#	add number of recipient with at least one prob transmitter
	ans		<- rbind(ans, tmp3, use.names=TRUE) 
	#	get dummy to sort
	set(ans, NULL, 'totalsort',  ans[, as.numeric(!grepl('total',factor))])
	#	get factor names
	ans		<- merge(ans, subset(factors, select=c(factor, factor.legend)), by='factor', all.x=1)		
	set(ans, ans[, which(is.na(factor.legend))], 'factor.legend', '')
	set(ans, ans[, which(grepl('total',factor))], 'factor', '')
	#	get column names
	tmp		<- copy(tperiod.info)
	tmp[, window:= tmp[,paste(gsub('\n','',t.period.min),'-',gsub('\n','',t.period.max),sep='')]]
	setnames(ans, tmp[, as.character(t.period)], tmp[, window])
	setnames(ans, 'All', paste(tmp[1,t.period.min],'-',tmp[4,t.period.max],sep=''))
	#	reorder	
	set(ans, NULL, 'group', ans[, factor(group, levels=c('','U_total','D_total','A_total'), labels=c('','U_total','D_total','A_total'))])
	set(ans, NULL, 'stat', ans[, as.character(stat)])
	setkey(ans, group, totalsort, factor, stat)
	
	#	write to file
	file			<- paste(outdir, '/', outfile, '_', gsub('/',':',insignat),'_',df[1, method.recentctime],'_',df[1, method.denom], '_',df[1, method.dating], '_',df[1, method.brl],'_',df[1, method.risk],'_table_PROP',".csv", sep='')
	cat(paste('\nsave to file',file))
	write.csv(subset(ans, select=c(10, 4, 5, 6, 7, 8)), file=file, eol="\r\n", row.names=FALSE)
	#in Excel replace COMMA with ,
	#right click on table or Apple-1 and set wrap text
}

