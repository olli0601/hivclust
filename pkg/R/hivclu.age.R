######################################################################################
hivc.prog.age_props_univariate<- function()
{	
	require(reshape2)
	require(data.table)
	
	require(ape)
	#stop()
	indir					<- paste(DATA,"fisheretal_data",sep='/')		
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	outdir					<- paste(DATA,"tpairs_age",sep='/')
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
		method					<- '3p'
		method.recentctime		<- '2011-01-01'
		method.nodectime		<- 'any'
		method.risk				<- 'm5A.tp4'
		method.Acute			<- 'higher'	#'central'#'empirical'
		method.minQLowerU		<- 0.148
		method.use.AcuteSpec	<- 1
		method.brl.bwhost		<- 2
		method.lRNA.supp		<- 100
		method.thresh.pcoal		<- 0.2
		method.minLowerUWithNegT<- 1
		method.cut.brl			<- Inf		#does not make a difference because compatibility test kills these anyway
		method.tpcut			<- 7
		method.PDT				<- 'SEQ'	# 'PDT'	
		method.thresh.bs		<- 0.8
		method.realloc			<- NA
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
						{	switch(substr(arg,2,15),
									method.cut.brl= return(as.numeric(substr(arg,17,nchar(arg)))),NA)	}))
		if(length(tmp)>0) method.cut.brl<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,17),
									method.thresh.bs= return(as.numeric(substr(arg,19,nchar(arg)))),NA)	}))
		if(length(tmp)>0) method.thresh.bs<- tmp[1]			
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,25),
									method.minLowerUWithNegT= return(as.numeric(substr(arg,27,nchar(arg)))),NA)	}))
		if(length(tmp)>0) method.minLowerUWithNegT<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,15),
									method.realloc= return(substr(arg,17,nchar(arg))),NA)	}))
		if(length(tmp)>0) method.realloc<- tmp[1]
	}	
	clu.infile			<- infile
	clu.indir			<- indir
	clu.insignat		<- insignat	
	t.recent.endctime	<- hivc.db.Date2numeric(as.Date(method.recentctime))	
	t.recent.endctime	<- floor(t.recent.endctime) + floor( (t.recent.endctime%%1)*100 %/% (t.period*100) ) * t.period	
	outfile				<- paste( outfile, ifelse(t.recent.endctime==t.endctime,'',paste('_',t.recent.endctime,sep='')), sep='')	
	#if(method.thresh.bs!=0.8)
	#{		
	clu.infilexml.opt	<- paste(clu.infilexml.opt,'_bs',method.thresh.bs,'_brl1000',sep='')
	cat(paste('\nusing file', clu.infilexml.opt))
	#}	
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
		print(method.cut.brl)
		print(method.realloc)
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
	if(method.cut.brl!=0.08)
		method				<- paste(method,'b',method.cut.brl,sep='')	
	if(method.thresh.bs!=0.8)
		method				<- paste(method,'s',method.thresh.bs,sep='')
	if(method.tpcut==4)
		tp.cut				<- NULL
	if(method.tpcut==7)
	{
		tp.cut				<- c(-Inf, 2006.5, 2008, 2009.5, 2010, 2010.5, 2011)
		method				<- paste(method,'T',method.tpcut,sep='')
	}			
	adjust.AcuteByNegT		<- 1
	any.pos.grace.yr		<- Inf	
	method.lRNA.supp		<- log10(method.lRNA.supp)	
	
	
	#	check if we are done	
	if(0 && resume)
	{		
		files		<- list.files(outdir)			
		if(length(files))
			files	<- files[ sapply(files, function(x) grepl(outfile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(paste('Yscore',method,sep=''), x, fixed=1) & grepl(paste('denom',method.PDT,sep=''), x) & !grepl('tables', x, fixed=1) & grepl(paste(method.risk,'.R',sep=''),x, fixed=1)  ) ]		
		stopifnot(length(files)==0)		
	}		
	#	check if we have precomputed tables
	X.tables				<- age.get.Xtables(method, method.PDT, method.risk, outdir, outfile, insignat)
	#	if no tables, precompute tables and stop (because mem intensive)
	#	otherwise return stratified YX data.table and continue
	tmp	<- age.precompute(	indir, indircov, infile.cov.study, infile.viro.study, infile.immu.study, infile.treatment.study, infile.cov.all, infile.viro.all, infile.immu.all, infile.treatment.all, infile.trm.model,
			clu.indir, clu.insignat, clu.infile,
			infile, infiletree, insignat, clu.infilexml.opt, clu.infilexml.template,
			method, method.recentctime, method.nodectime, method.risk, method.Acute, method.minQLowerU, method.use.AcuteSpec, method.brl.bwhost, method.lRNA.supp, method.thresh.pcoal, method.minLowerUWithNegT, method.tpcut, method.PDT, method.cut.brl, tp.cut, adjust.AcuteByNegT, any.pos.grace.yr, dur.Acute, method.thresh.bs, 
			outdir, outfile,
			t.period, t.recent.startctime, t.endctime, t.recent.endctime,
			X.tables, resume, verbose
	)
	predict.t2inf	<- tmp$predict.t2inf
	t2inf.args		<- tmp$t2inf.args
	df.all			<- copy(tmp$df.all)
	YX				<- copy(tmp$YX)
	Y.brl.bs		<- copy(tmp$Y.brl.bs)
	gc()
#STOP2
#stop()	
	stopifnot(is.null(X.tables)==FALSE)
	#
	#	for each time period, estimate N transmitted etc
	#
	resume				<- 1
	bs.n				<- 1e3
	if(method.PDT=='')
		method.PDT<- 'PDT'
	if(grepl('m[0-9]',method.risk))
	{
		tmp			<- substr(regmatches(method.risk,regexpr('m[0-9]', method.risk)),2,2)
		save.file	<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_denom',method.PDT,'_model',tmp,'_',method.risk,'.R',sep='')
		#	estimate as stratified
		YXe			<- project.athena.Fisheretal.estimate.risk.wrap(YX, Y.brl.bs, X.tables, tperiod.info, plot.file.or=NA, bs.n=bs.n, resume=resume, save.file=save.file, method.risk=method.risk)
		if( as.numeric(substring(regmatches(method.risk,regexpr('tp[0-9]', method.risk)),3))<4 )
		{
			#	pool ART stages, re-estimate
			tmp			<- substr(regmatches(method.risk,regexpr('m[0-9]', method.risk)),2,2)
			save.file	<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_denom',method.PDT,'_model',tmp,'_',sep='')
			save.file	<- paste(save.file, substr(method.risk, 1, regexpr('tp[0-9]', method.risk)-1), 'ARTstarted.', regmatches(method.risk,regexpr('tp[0-9]', method.risk)), '.R', sep='')	
			tmp			<- project.athena.Fisheretal.poolARTstarted( YXe, method.risk, save.file=save.file)
			#	pool into groups U D A, re-estimate
			tmp			<- substr(regmatches(method.risk,regexpr('m[0-9]', method.risk)),2,2)
			save.file	<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_denom',method.PDT,'_model',tmp,'_',sep='')
			save.file	<- paste(save.file, substr(method.risk, 1, regexpr('tp[0-9]', method.risk)-1), 'GroupsUDA.', regmatches(method.risk,regexpr('tp[0-9]', method.risk)), '.R', sep='')	
			tmp			<- project.athena.Fisheretal.poolIntoGroups( YXe, save.file=save.file)	
		}		
	}		
	#	see if we can pool results for tperiod 4
	if(method.tpcut%in%c(7) & as.numeric(substring(regmatches(method.risk,regexpr('tp[0-9]', method.risk)),3))>=4)
	{
		YXe			<- project.athena.Fisheretal.pool.TP4(outdir, outfile, insignat, method, method.PDT, method.risk, resume=1)
		#	if success, 
		if(!is.null(YXe))
		{
			method.risk	<- gsub('tp[0-9]','tp4',method.risk)
			#	pool ART stages, re-estimate
			tmp			<- substr(regmatches(method.risk,regexpr('m[0-9]', method.risk)),2,2)
			save.file	<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_denom',method.PDT,'_model',tmp,'_',sep='')
			save.file	<- paste(save.file, substr(method.risk, 1, regexpr('tp[0-9]', method.risk)-1), 'ARTstarted.', regmatches(method.risk,regexpr('tp[0-9]', method.risk)), '.R', sep='')	
			tmp			<- project.athena.Fisheretal.poolARTstarted( YXe, method.risk, save.file=save.file)
			#	pool into groups U D A, re-estimate
			tmp			<- substr(regmatches(method.risk,regexpr('m[0-9]', method.risk)),2,2)
			save.file	<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_denom',method.PDT,'_model',tmp,'_',sep='')
			save.file	<- paste(save.file, substr(method.risk, 1, regexpr('tp[0-9]', method.risk)-1), 'GroupsUDA.', regmatches(method.risk,regexpr('tp[0-9]', method.risk)), '.R', sep='')	
			tmp			<- project.athena.Fisheretal.poolIntoGroups( YXe, save.file=save.file)
#STOP3
#stop()
			#
			#	TP4 hypothetical scenarios
			#
			if(!is.na(method.realloc))
			{
				tmp			<- substr(regmatches(method.risk,regexpr('m[0-9]', method.risk)),2,2)
				save.file	<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore',method,'_denom',method.PDT,'_model',tmp,'_',sep='')
				save.file	<- paste(save.file, substr(method.risk, 1, regexpr('tp[0-9]', method.risk)-1), 'Hypo', gsub('+','',method.realloc), '.', regmatches(method.risk,regexpr('tp[0-9]', method.risk)), '.R', sep='')		
				tmp			<- project.athena.Fisheretal.Hypo.run(YXe, method.risk, predict.t2inf=predict.t2inf, t2inf.args=t2inf.args, df.all=df.all, method.realloc=method.realloc, t.period=t.period,  use.YXf= 1, bs.n=1e3, save.file=save.file, resume=resume)
			}														
		}		
	}
	
	
#SSS	
	
	quit('no')
	
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
######################################################################################
age.get.Xtables<- function(method, method.PDT, method.risk, outdir, outfile, insignat)
{
	X.tables			<- NULL
	if(1)
	{
		save.file		<- NA
		if(grepl('m5A',method.risk))	save.file	<- 'm5A'
		if(grepl('m5B',method.risk))	save.file	<- 'm5B'
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
		}
	}
	X.tables
}
######################################################################################
censoring.frac.Age_253045.Stage_UA_UC_D_T_F<- function()
{		
	require(Hmisc)
	indir		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3pa1H1.48C2V100bInfT7_CtablesSEQ_m5A.R'	
	plot.file	<- paste(indir,'/',gsub('\\.R','\\.pdf',gsub('CtablesSEQ','CENSORINGFRACTION',infile)),sep='')
	
	load(paste(indir,'/',infile,sep=''))
	ct			<- copy(ans$cens.table)
	setkey(ct, stat, t.period, risk, factor)
	ct			<- unique(ct)
	ctb			<- copy(ans$cens.table.bs)
	setkey(ctb, stat, t.period, risk, factor)
	ctb			<- unique(ctb)
	 	
	#	ct.bs contains counts at times [t1bs, t2bs]
	#	divide these by the counts at times [t1,t2] --> bootstrap fraction
	set(ctb, NULL, 'c', ctb[, nc/n])
	set(ctb, ctb[, which(n==0)], 'c', 1.)
	set(ctb, NULL, 't2c', ctb[, cens.t-t.period.max.bs])
	if(!is.na(plot.file))
	{
		ct.plot		<- subset(ctb, select=c(risk, factor, factor2, c, t2c))	
		#ct.plot		<- merge( ct.plot, factors, by='factor2')
		ct.plot		<- merge(ct.plot, ct.plot[, list(mc=median(c)), by='factor'], by='factor')
		#if(method.group=='cascade')
		ct.plot		<- subset(ct.plot, grepl('^U',factor2))
		set(ct.plot, NULL, 't2cf', ct.plot[, factor(round(t2c,d=2))])
		ct.plot[, factor.legend:= factor2]
		ct.plot[, age:= ct.plot[, factor(sapply(strsplit(factor2,'_'),'[[',2))]]
		ct.plot[, stage:= ct.plot[, factor(sapply(strsplit(factor2,'_'),'[[',1))]]
		setkey(ct.plot, factor.legend)
		ggplot(ct.plot, aes(x=2011-t2c, y=100*c, colour=factor.legend)) + geom_boxplot(aes(position=factor(t2c)), outlier.shape = NA) +
				labs(x=expression('end time of observation period ('*t[2]^k*' )'), y='fraction of non-censored\npotential transmission intervals\n( % )') +
				#scale_colour_manual(values=ct.plot[, unique(factor.color)], guide = FALSE) +
				scale_colour_discrete(guide=FALSE) +
				scale_y_continuous(breaks=seq(0,100,20), lim=c(0,100), expand=c(0,1) ) + theme_bw() +
				theme(panel.grid.major=element_line(colour="grey70", size=0.4), legend.position = "bottom") + 
				facet_grid(stage ~ age, margins=FALSE)
		ggsave(file=plot.file, w=10, h=6)		
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
	ctb		<- subset(ctb, select=c(t.period, risk, factor, factor2, p.cens, p.cens.95l, p.cens.95u))
	ctn		<- subset(ctn, select=c(t.period, risk, factor, factor2, n.adj))
	list( ctb=ctb, ctn=ctn )
}
######################################################################################
sampling.frac.Age_253045.Stage_UA_UC_D_T_F<- function()
{		
	require(Hmisc)
	indir	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile	<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3pa1H1.48C2V100bInfT7_StablesSEQ_m5A.R'	
	
	load(paste(indir,'/',infile,sep=''))
	df		<- ans$risk.table	
	set(df, NULL, 'factor', df[, as.character(factor)])	
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])
	set(df, NULL, 'factor', df[, factor(factor)])
	dfp		<- df[, 	{
				#BINCONF not appropriate because potential transmission intervals are not independent
				p.seq.cis	<- as.double(binconf( sum(n[stat=='X.seq']), sum(n[stat=='X.msm']), alpha=0.05, method= "wilson", include.x=FALSE, include.n=FALSE, return.df=FALSE))
				list( 	p.seq= sum(n[stat=='X.seq']) / sum(n[stat=='X.msm']),
						p.seq.l= p.seq.cis[2],
						p.seq.u= p.seq.cis[3] )
			}, by=c('factor','t.period')]
	dfp		<- merge(dfp, df[, list( factor=factor[stat=='X.msm'], pp.coh=n[stat=='X.msm']/sum(n[stat=='X.msm']), pp.pos=n[stat=='YX']/sum(n[stat=='YX']) ), by='t.period'], by=c('factor','t.period'))	
	set(dfp, NULL, 'p.seq', dfp[, round(p.seq, d=3)]*100)		
	set(dfp, NULL, 't.period.long', dfp[, factor(t.period, levels=c('1','2','3','4'), labels=c("96/07-06/06", "06/07-07/12", "08/01-09/06", "09/07-10/12"))])
	#	separate stage and age
	set(dfp, NULL, 'age', dfp[, sapply(strsplit(as.character(factor),'_'),'[[',2)])
	set(dfp, NULL, 'stage', dfp[, sapply(strsplit(as.character(factor),'_'),'[[',1)])
	set(dfp, NULL, 'stage', dfp[, factor(stage, levels=c('UA','UC','D','T','L'))])
	set(dfp, NULL, 'age', dfp[, factor(age, levels=c('(-1,25]','(25,30]','(30,45]','(45,100]'))])
	#	plot p.seq		
	ggplot(dfp, aes(x=p.seq, y=t.period.long, shape=age, colour=age)) + geom_point(size=3) +
				scale_x_continuous(breaks=seq(0,100,10), minor_breaks=NULL, limit=c(0,100), expand=c(0,0)) +
				scale_shape_discrete(guide=FALSE) +
				#scale_colour_manual(values=dfp[, unique(factor.color)], guide=FALSE) +
				theme_bw() +				
				theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=10), legend.position='bottom', panel.grid.major.x=element_line(colour="grey70", size=0.4), panel.grid.minor.x=element_line(colour="grey70", size=0.4), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), panel.margin = unit(2, "lines")) + 				
				labs(y='', x='overlap intervals\nof a potential transmitter with a sequence\n(%)', colour='age of potential transmitter\nat transmission interval') +
				facet_grid(~stage)  				
	ggsave(file=paste(indir, '/', gsub('\\.R','_STRATSAMPLING_cmp_Age.pdf',infile), sep=''), w=15, h=5)
	
	ggplot(dfp, aes(x=p.seq, y=t.period.long, shape=stage, colour=stage)) + geom_point(size=3) +
			scale_x_continuous(breaks=seq(0,100,10), minor_breaks=NULL, limit=c(0,100), expand=c(0,0)) +
			scale_shape_discrete(guide=FALSE) +
			#scale_colour_manual(values=dfp[, unique(factor.color)], guide=FALSE) +
			theme_bw() +				
			theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=10), legend.position='bottom', panel.grid.major.x=element_line(colour="grey70", size=0.4), panel.grid.minor.x=element_line(colour="grey70", size=0.4), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), panel.margin = unit(2, "lines")) + 				
			labs(y='', x='overlap intervals\nof a potential transmitter with a sequence\n(%)', colour='stage of potential transmitter\nat transmission interval') +
			facet_grid(~age)  				
	ggsave(file=paste(indir, '/', gsub('\\.R','_STRATSAMPLING_cmp_Stage.pdf',infile), sep=''), w=15, h=5)
	
	#
	#	sequence sampling does not vary much by age of potential transmitter
	#
}
######################################################################################
sampling.get.dataset<- function(pt.df, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, df.treatment.allmsm, t.period, t.endctime, method.lRNA.supp=2, method.lRNA.nsupp=4, contact.grace=1.5)
{	
	#pt.df	<- data.table(Patient=X.msm[, unique(t.Patient)])
	ptd.df	<- merge(pt.df, subset( df.all.allmsm, select=c(Patient, PosSeqT, RegionHospital, DateBorn, DateDied, isAcute, AnyPos_T1, AnyT_T1)), by='Patient')
	setkey(ptd.df, Patient, PosSeqT)
	setkey(ptd.df, Patient)
	ptd.df	<- unique(ptd.df)
	#	With sequence	
	ptd.df[, SQD:= as.integer(!is.na(PosSeqT))]
	#	Acute at Diag
	set(ptd.df, ptd.df[, which(isAcute!='Yes' | is.na(isAcute))], 'isAcute', 'NI')
	set(ptd.df, NULL, 'isAcute', ptd.df[, factor(isAcute)])
	#	Acute at Diag in Amsterdam or North
	set(ptd.df, NULL, 'isAcuteAN', ptd.df[, as.character(isAcute)])
	set(ptd.df, ptd.df[, which(isAcute!='Yes' | is.na(isAcute) | RegionHospital%in%c('S','W','E'))], 'isAcuteAN', 'NI')
	set(ptd.df, NULL, 'isAcuteAN', ptd.df[, factor(isAcuteAN)])
	
	#	time since diagnosis
	set(ptd.df, ptd.df[, which(is.na(DateDied))], 'DateDied', t.endctime)
	ptd.df[, DT_Diag:= DateDied-AnyPos_T1]
	stopifnot( ptd.df[, any(DT_Diag>=0)] )
	#	time since ART start
	ptd.df[, DT_ART:= DateDied-AnyT_T1]
	#	ever on ART
	ptd.df[, Ever_ART:= factor(as.numeric(!is.na(AnyT_T1)), levels=c(0,1), labels=c('N','Y'))]	
	#	ever not suppressed for more than 3 mnths after first suppression 
	tmp		<- subset( df.viro.allmsm, select=c(Patient, PosRNA, lRNA) )
	set(tmp, NULL, 'PosRNA', hivc.db.Date2numeric(tmp[,PosRNA]))
	tmp		<- merge(subset(ptd.df, Ever_ART=='Y',c(Patient, AnyT_T1)),tmp,by='Patient')
	#	add VS_T1: time of first viral suppression
	tmp		<- merge(tmp, tmp[, list(VS_T1=PosRNA[ head(which(lRNA<method.lRNA.supp),1)]), by='Patient'], by='Patient', all.x=1)
	#	check for which patient we have two consecutive msrmnts above 10,000 cps
	z		<- subset(tmp, !is.na(VS_T1) & PosRNA>=VS_T1 )[, 
			{				
				list( VNS_ANY= grepl('11',paste( as.numeric(lRNA>method.lRNA.nsupp), collapse='' )), VS_ALL= all(lRNA<=method.lRNA.supp) )				
			}, by='Patient']
	tmp		<- merge(tmp, z, by='Patient',all.x=T)
	setkey(tmp, Patient)
	tmp		<- unique(tmp)
	ptd.df	<- merge(ptd.df, subset(tmp, select=c(Patient, VNS_ANY, VS_ALL)), by='Patient', all.x=T)
	#	ever not in contact	
	incare		<- merge( subset(ptd.df, select=Patient), unique( subset(df.all.allmsm, select=c(Patient, AnyPos_T1, DateDied)) ), by='Patient' )	
	set(incare, NULL, 'AnyPos_T1', incare[, floor(AnyPos_T1) + round( (AnyPos_T1%%1)*100 %/% (t.period*100) ) * t.period] )
	set(incare, incare[,which(is.na(DateDied))], 'DateDied', t.endctime)
	set(incare, NULL, 'DateDied', incare[, floor(DateDied) + round( (DateDied%%1)*100 %/% (t.period*100) ) * t.period] )
	incare		<- subset( incare, AnyPos_T1<DateDied )
	incare.t	<- incare[, list(t= seq(AnyPos_T1, DateDied-t.period, by=t.period)),by='Patient']
	tmp			<- subset(ptd.df, select=Patient)
	setnames(tmp, 'Patient', 't.Patient')
	setnames(incare.t, 'Patient', 't.Patient')
	tmp			<- project.athena.Fisheretal.X.nocontact(incare.t, df.viro.allmsm, df.immu.allmsm, df.treatment.allmsm, tmp, df.all.allmsm, contact.grace=contact.grace, t.period=t.period, t.endctime=t.endctime)
	setnames(tmp, 't.Patient', 'Patient')
	ptd.df		<- merge(ptd.df, tmp[, list(NC_ANY=any(contact=='No')), by='Patient'], by='Patient',all.x=T)	
	set(ptd.df, NULL, 'VNS_ANY', ptd.df[, factor(as.numeric(VNS_ANY), levels=c(0,1), labels=c('N','Y'))])
	set(ptd.df, NULL, 'VS_ALL', ptd.df[, factor(as.numeric(VS_ALL), levels=c(0,1), labels=c('N','Y'))])
	set(ptd.df, NULL, 'NC_ANY', ptd.df[, factor(as.numeric(NC_ANY), levels=c(0,1), labels=c('N','Y'))])
	#	handle NA's
	set(ptd.df, ptd.df[, which(is.na(DT_ART))],'DT_ART',0)
	set(ptd.df, ptd.df[, which(is.na(VNS_ANY))],'VNS_ANY','N')
	set(ptd.df, ptd.df[, which(is.na(VS_ALL))],'VS_ALL','N')
	set(ptd.df, ptd.df[, which(is.na(NC_ANY))],'NC_ANY','N')
	#	Acute at Diag in Amsterdam or North, only when in contact
	set(ptd.df, NULL, 'isAcuteANC', ptd.df[, as.character(isAcute)])
	set(ptd.df, ptd.df[, which(isAcute!='Yes' | is.na(isAcute) | RegionHospital%in%c('S','W','E') | NC_ANY=='Y')], 'isAcuteANC', 'NI')
	set(ptd.df, NULL, 'isAcuteANC', ptd.df[, factor(isAcuteANC)])	
	#	save
	if(!is.na(save.file))		
		save(ptd.df, file=save.file)
	ptd.df
}
######################################################################################
sampling.model.150722<- function(ptd.df)
{	
	#	use model m12	
	set(ptr.df, NULL, 'p.seq', NA_real_)
	z		<- subset(ptr.df, NC_ANY=='Y', select=which(colnames(ptr.df)!='p.seq'))
	ms12a	<- gamlss(formula= SQD ~ RegionHospital*ns(AnyPos_T1, df=1) + VS_ALL  + isAcuteANC - 1, family=BI(), data=z)
	z[, p.seq:=predict(ms12a, type='response')]		
	set(ptr.df, ptr.df[, which(NC_ANY=='Y')], 'p.seq', z[, p.seq])
	z		<- subset(ptr.df, isAcuteANC=='Yes', select=which(colnames(ptr.df)!='p.seq'))
	ms12b	<- gamlss(formula= SQD ~ ns(AnyPos_T1, df=1) + VS_ALL  - 1, family=BI(), data=z)
	z[, p.seq:=predict(ms12b, type='response')]	
	set(ptr.df, ptr.df[, which(isAcuteANC=='Yes')], 'p.seq', z[, p.seq])
	z		<- subset(ptr.df, isAcuteANC=='NI' & NC_ANY=='N', select=which(colnames(ptr.df)!='p.seq'))
	ms12c	<- gamlss(formula= SQD ~ RegionHospital*ns(AnyPos_T1, df=6) + VS_ALL - 1, family=BI(), data=z)	
	z[, p.seq:=predict(ms12c, type='response')]
	set(ptr.df, ptr.df[, which(isAcuteANC=='NI' & NC_ANY=='N')], 'p.seq', z[, p.seq])
	
	ptr.df
}
######################################################################################
sampling.model.exploreoptions<- function(ptd.df)
{	
	#
	#	explore first regression models
	#
	require(zoo)
	outdir	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	oufile	<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011'
	
	#	is age at diagnosis an explanatory variable?
	#	use 2009.2 to kill censoring effects
	pta.df	<- subset(ptd.df, AnyPos_T1>2006.5 & AnyPos_T1<=2009.2 & VNS_ANY=='N' & NC_ANY=='N' & RegionHospital%in%c('Amst','W'))
	pta.df[, Age_T1:= AnyPos_T1-DateBorn]
	setkey(pta.df, Age_T1)
	pta.df[, SQD.rm:=pta.df[, rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE)]]	
	ggplot(pta.df, aes(x=Age_T1, y=SQD.rm)) + geom_line()
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ageatdiag.pdf', sep=''), w=5, h=5)
	#	once censoring is 'excluded' as an effect ie focuse on <2009, there seems to be a 'dip' close to age 30,
	#	but overall, there is no additional age effect

	#	is time since diagnosis an explanatory variable?
	setkey(pta.df, DT_Diag)
	pta.df[, SQD.rm:=pta.df[, rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE)]]	
	ggplot(pta.df, aes(x=DT_Diag, y=SQD.rm)) + geom_line()
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_tsincediag.pdf', sep=''), w=5, h=5)
	#	once censoring is 'excluded' as an effect ie focuse on <2009,
	#	time since diag has no strong impact on rolling mean
	
	ptr.df	<- copy(ptd.df)
	ptr.df[, PosSeqT:=NULL]
	ptr.df[, AnyT_T1:=NULL]			
	
	setkey(ptr.df, AnyPos_T1)	
	ptr.df[, SQD.rm:=ptr.df[, rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE)]]	
	ms1		<- gamlss(formula= SQD ~ ns(AnyPos_T1, df=5), family=BI(), data=ptr.df)
	ptr.df[, MS1:=predict(ms1, type='response')]	
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm','MS1'), id.vars=c('AnyPos_T1')), aes(x=AnyPos_T1, y=value, colour=variable)) + geom_line()
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms1.pdf', sep=''), w=5, h=5)
	#	150721:	ns(AnyPos_T1, df=5) worked better than bs(AnyPos_T1, degree=5)

	#	take time dependent regression as baseline model; 
	#	try to improve based on ACUTE
	setkey(ptr.df, isAcute, AnyPos_T1)
	tmp		<- 	ptr.df[, list( Patient=Patient, SQD.rm2=rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE) ), by='isAcute']
	ptr.df	<- merge(ptr.df, tmp, by=c('Patient','isAcute'))
	ms2		<- gamlss(formula= SQD ~ ns(AnyPos_T1, df=5) + isAcute-1, family=BI(), data=ptr.df)
	ptr.df[, MS2:=predict(ms2, type='response')]	
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm2','MS2'), id.vars=c('AnyPos_T1','isAcute')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(~isAcute)
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms2.pdf', sep=''), w=8, h=5)
	
	#	take time dependent regression as baseline model; 
	#	try to improve based on ACUTE-AN
	setkey(ptr.df, isAcuteAN, AnyPos_T1)
	tmp		<- 	ptr.df[, list( Patient=Patient, SQD.rm2b=rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE) ), by='isAcuteAN']
	ptr.df	<- merge(ptr.df, tmp, by=c('Patient','isAcuteAN'))	
	ms2b	<- gamlss(formula= SQD ~ ns(AnyPos_T1, df=5) + isAcuteAN - 1, family=BI(), data=ptr.df)
	ptr.df[, MS2b:=predict(ms2b, type='response')]	
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm2b','MS2b'), id.vars=c('AnyPos_T1','isAcuteAN')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(~isAcuteAN)
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms2b.pdf', sep=''), w=8, h=5)
	
	
	#	try to improve based on Not Suppressed
	setkey(ptr.df, VNS_ANY, AnyPos_T1)
	tmp		<- 	ptr.df[, list( Patient=Patient, SQD.rm3=rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE) ), by='VNS_ANY']
	ptr.df	<- merge(ptr.df, tmp, by=c('Patient','VNS_ANY'))
	ms3		<- gamlss(formula= SQD ~ ns(AnyPos_T1, df=5) + VNS_ANY-1, family=BI(), data=ptr.df)
	ptr.df[, MS3:=predict(ms3, type='response')]	
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm3','MS3'), id.vars=c('AnyPos_T1','VNS_ANY')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(~VNS_ANY)
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms3.pdf', sep=''), w=8, h=5)
	
	#	try to improve based on Suppressed
	setkey(ptr.df, VS_ALL, AnyPos_T1)
	tmp		<- 	ptr.df[, list( Patient=Patient, SQD.rm4=rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE) ), by='VS_ALL']
	ptr.df	<- merge(ptr.df, tmp, by=c('Patient','VS_ALL'))
	ms4		<- gamlss(formula= SQD ~ ns(AnyPos_T1, df=5) + VS_ALL-1, family=BI(), data=ptr.df)
	ptr.df[, MS4:=predict(ms4, type='response')]	
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm4','MS4'), id.vars=c('AnyPos_T1','VS_ALL')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(~VS_ALL)
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms4.pdf', sep=''), w=8, h=5)
	#	suppressed has lower AIC
	
	#	try to improve based on region hospital
	setkey(ptr.df, RegionHospital, AnyPos_T1)
	tmp		<- 	ptr.df[, list( Patient=Patient, SQD.rm5=rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE) ), by='RegionHospital']
	ptr.df	<- merge(ptr.df, tmp, by=c('Patient','RegionHospital'))
	ms5		<- gamlss(formula= SQD ~ ns(AnyPos_T1, df=5) + RegionHospital-1, family=BI(), data=ptr.df)
	ptr.df[, MS5:=predict(ms5, type='response')]	
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm5','MS5'), id.vars=c('AnyPos_T1','RegionHospital')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(~RegionHospital)
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms5.pdf', sep=''), w=8, h=5)
	#	baseline model seems OK except for South -- try separate model there	
	setkey(ptr.df, RegionHospital, AnyPos_T1)
	#ms5b	<- gamlss(formula= SQD ~ I(RegionHospital=='S')*lo(~AnyPos_T1, span=0.3) + RegionHospital-1, family=BI(), data=ptr.df)
	#ms5b	<- gamlss(formula= SQD ~ I(RegionHospital=='S')*ns(AnyPos_T1, df=7) + RegionHospital-1, family=BI(), data=ptr.df)
	ms5b	<- gamlss(formula= SQD ~ RegionHospital*ns(AnyPos_T1, df=6), family=BI(), data=ptr.df)
	ptr.df[, MS5b:=predict(ms5b, type='response')]	
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm5','MS5b'), id.vars=c('AnyPos_T1','RegionHospital')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(~RegionHospital)
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms5b.pdf', sep=''), w=8, h=5)
	
	#z		<- subset(ptr.df, RegionHospital=='S')
	#ms5b	<- gamlss(formula= SQD ~ lo(~AnyPos_T1, span=0.3), family=BI(), data=z)
	#z[, MS5b:=predict(ms5b, type='response')]
	#ggplot(melt(z, measure.vars=c('SQD.rm5','MS5b'), id.vars=c('AnyPos_T1','RegionHospital')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(~RegionHospital)
	
	#	try to improve based on no contact
	setkey(ptr.df, NC_ANY, AnyPos_T1)
	tmp		<- 	ptr.df[, list( Patient=Patient, SQD.rm6=rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE) ), by='NC_ANY']
	ptr.df	<- merge(ptr.df, tmp, by=c('Patient','NC_ANY'))
	ms6		<- gamlss(formula= SQD ~ ns(AnyPos_T1, df=5) + NC_ANY-1, family=BI(), data=ptr.df)
	ptr.df[, MS6:=predict(ms6, type='response')]	
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm6','MS6'), id.vars=c('AnyPos_T1','NC_ANY')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(~NC_ANY)
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms6.pdf', sep=''), w=8, h=5)
	
	#	combine VS_ALL > region >  no contact
	setkey(ptr.df, VS_ALL, RegionHospital, NC_ANY, AnyPos_T1)
	tmp		<- 	ptr.df[, list( Patient=Patient, SQD.rm7=rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE) ), by=c('VS_ALL','RegionHospital','NC_ANY')]
	ptr.df	<- merge(ptr.df, tmp, by=c('Patient','VS_ALL','RegionHospital','NC_ANY'))
	ms7		<- gamlss(formula= SQD ~ ns(AnyPos_T1, df=5) + VS_ALL + RegionHospital + NC_ANY-1, family=BI(), data=ptr.df)
	ptr.df[, MS7:=predict(ms7, type='response')]	
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm7','MS7'), id.vars=c('AnyPos_T1','VS_ALL','RegionHospital','NC_ANY')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(RegionHospital~VS_ALL+NC_ANY)
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms7.pdf', sep=''), w=10, h=10)

	#	combine VS_ALL > region >  no contact;	use region specific sampling model
	setkey(ptr.df, VS_ALL, RegionHospital, NC_ANY, AnyPos_T1)
	tmp		<- 	ptr.df[, list( Patient=Patient, SQD.rm8=rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE) ), by=c('VS_ALL','RegionHospital','NC_ANY')]
	ptr.df	<- merge(ptr.df, tmp, by=c('Patient','VS_ALL','RegionHospital','NC_ANY'))
	ms8		<- gamlss(formula= SQD ~ RegionHospital*ns(AnyPos_T1, df=6) + VS_ALL + NC_ANY - 1, family=BI(), data=ptr.df)
	ptr.df[, MS8:=predict(ms8, type='response')]	
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm8','MS8'), id.vars=c('AnyPos_T1','VS_ALL','RegionHospital','NC_ANY')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(RegionHospital~VS_ALL+NC_ANY)
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms8.pdf', sep=''), w=10, h=10)
	
	#	combine VS_ALL > region > isAcuteAN >  no contact;	use region specific sampling model
	setkey(ptr.df, VS_ALL, RegionHospital, NC_ANY, isAcuteAN, AnyPos_T1)
	tmp		<- 	ptr.df[, list( Patient=Patient, SQD.rm10=rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE) ), by=c('VS_ALL','RegionHospital','NC_ANY','isAcuteAN')]
	ptr.df	<- merge(ptr.df, tmp, by=c('Patient','VS_ALL','RegionHospital','NC_ANY','isAcuteAN'))
	ms10	<- gamlss(formula= SQD ~ RegionHospital*ns(AnyPos_T1, df=6) + VS_ALL + NC_ANY + isAcuteAN - 1, family=BI(), data=ptr.df)
	ptr.df[, MS10:=predict(ms10, type='response')]	
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm10','MS10'), id.vars=c('AnyPos_T1','VS_ALL','RegionHospital','NC_ANY','isAcuteAN')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(RegionHospital~VS_ALL+NC_ANY+isAcuteAN)
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms10.pdf', sep=''), w=12, h=10)
	
	#	combine VS_ALL > region > isAcuteANC >  no contact;	use region specific sampling model
	setkey(ptr.df, VS_ALL, RegionHospital, NC_ANY, isAcuteANC, AnyPos_T1)
	tmp		<- 	ptr.df[, list( Patient=Patient, SQD.rm11=rollapply(SQD, width=100, FUN=mean, align="center", partial=TRUE) ), by=c('VS_ALL','RegionHospital','NC_ANY','isAcuteANC')]
	ptr.df	<- merge(ptr.df, tmp, by=c('Patient','VS_ALL','RegionHospital','NC_ANY','isAcuteANC'))
	ms11	<- gamlss(formula= SQD ~ RegionHospital*ns(AnyPos_T1, df=6) + VS_ALL + NC_ANY + isAcuteANC - 1, family=BI(), data=ptr.df)
	ptr.df[, MS11:=predict(ms11, type='response')]	
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm11','MS11'), id.vars=c('AnyPos_T1','VS_ALL','RegionHospital','NC_ANY','isAcuteANC')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(RegionHospital~VS_ALL+NC_ANY+isAcuteANC)
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms11.pdf', sep=''), w=12, h=10)
	
	set(ptr.df, NULL, 'MS12', NA_real_)
	z		<- subset(ptr.df, NC_ANY=='Y', select=which(colnames(ptr.df)!='MS12'))
	ms12a	<- gamlss(formula= SQD ~ RegionHospital*ns(AnyPos_T1, df=1) + VS_ALL  + isAcuteANC - 1, family=BI(), data=z)
	z[, MS12:=predict(ms12a, type='response')]		
	set(ptr.df, ptr.df[, which(NC_ANY=='Y')], 'MS12', z[, MS12])
	z		<- subset(ptr.df, isAcuteANC=='Yes', select=which(colnames(ptr.df)!='MS12'))
	ms12b	<- gamlss(formula= SQD ~ ns(AnyPos_T1, df=1) + VS_ALL  - 1, family=BI(), data=z)
	z[, MS12:=predict(ms12b, type='response')]	
	set(ptr.df, ptr.df[, which(isAcuteANC=='Yes')], 'MS12', z[, MS12])
	z		<- subset(ptr.df, isAcuteANC=='NI' & NC_ANY=='N', select=which(colnames(ptr.df)!='MS12'))
	ms12c	<- gamlss(formula= SQD ~ RegionHospital*ns(AnyPos_T1, df=6) + VS_ALL - 1, family=BI(), data=z)	
	z[, MS12:=predict(ms12c, type='response')]
	set(ptr.df, ptr.df[, which(isAcuteANC=='NI' & NC_ANY=='N')], 'MS12', z[, MS12])
	ggplot(melt(ptr.df, measure.vars=c('SQD.rm11','MS12'), id.vars=c('AnyPos_T1','VS_ALL','RegionHospital','NC_ANY','isAcuteANC')), aes(x=AnyPos_T1, y=value, group=variable, colour=variable)) + 
			geom_line() +
			coord_cartesian(xlim=c(1990,2012)) +
			facet_grid(RegionHospital~VS_ALL+NC_ANY+isAcuteANC)
	ggsave(file=paste(outdir, '/', outfile, '_SAMPLINGMODEL_ms12.pdf', sep=''), w=12, h=10)
	
	#	sapply(list(ms1, ms2, ms2b, ms3, ms4, ms5, ms5b, ms6, ms7, ms8, ms9, ms10), Rsq)
	#	0.05869068 0.06241876 0.07190648 0.07340679 0.07538781 0.09400546 0.11839415 0.06397835 0.11540097 0.13952236 0.14305379 0.14303762
	#	sapply(list(ms1, ms2,  ms2b, ms3, ms4, ms5, ms5b, ms6, ms7, ms8, ms9, ms10), AIC)
	#	16169.58 16123.19 15999.18 15979.45 15953.35 15711.33 15428.61 16102.89 15423.94 15136.84 15088.69 15088.92

	#	not much impr of ms9 over ms8
	#	but spline for each region is much better than spline+region 

	#
	#	overall univariate effect: 
	#	region > VS_ALL > VNS_ANY > no contact > isAcute
	#	region > VS_ALL > VNS_ANY > isAcuteAN > no contact 
	#	use ms9
	#
}
######################################################################################
sampling.get.all.tables.args<- function(method)
{	
	tp				<- regmatches(method, regexpr('tp[0-9]', method))
	tp				<- ifelse(length(tp), paste('.',substr(tp, 3, 3),sep=''), '')			
	if(grepl('m5A',method))
	{
			factor.ref.v	<- paste('T_(30,45]',tp,sep='')
			risktp.col		<- 't.stAgeC.prd'
			risk.col		<- 't.stAgeC'			
	}
	if(grepl('m5B',method))
	{
		factor.ref.v	<- paste('T_(30,45]',tp,sep='')
		risktp.col		<- 't.stAgeC.prd'
		risk.col		<- 't.stAgeC'			
	}
	c('factor.ref.v'=factor.ref.v, 'risktp.col'=risktp.col, 'risk.col'=risk.col)
}
######################################################################################
censoring.get.all.tables.by.tperiod<- function(YX=NULL, X.seq=NULL, X.msm=NULL, X.clu=NULL, tperiod.info=NULL, resume=TRUE, save.file=NA, method=NA, risk.col=NA, risktp.col=NA, factor.ref.v=NA, bs.n=100, bs.cdelta.min=2, bs.cdelta.max=3)
{
	stopifnot(!is.na(risk.col), !is.na(risktp.col), !is.na(factor.ref.v))
	
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
		if(is.null(YX) || is.null(X.seq) || is.null(X.msm) || is.null(X.clu))
		{
			cat(paste('\nreturn NULL table', method))
			return(ans)	
		}			
		cat(paste('\ntables by method', method,'\nrisk.col', risk.col,'\nrisk.col.tp', risktp.col,'\nfactor.ref.v',factor.ref.v))
		#YX<- copy(YX.s); X.clu<- copy(X.clu.s); X.seq<- copy(X.seq.s); X.msm<- copy(X.msm.s)
		#	set stratifications. always use 'risktp.col' for censoring
		set(YX, NULL, 'stage', YX[[risktp.col]])				
		set(X.clu, NULL, 'stage', X.clu[[risktp.col]])
		set(X.seq, NULL, 'stage', X.seq[[risktp.col]])
		set(X.msm, NULL, 'stage', X.msm[[risktp.col]])
		set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])									
		set(X.seq, NULL, 'stage', X.seq[, factor(as.character(stage), levels=X.msm[, levels(stage)])])		
		set(X.clu, NULL, 'stage', X.clu[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
		set(YX,    NULL, 'stage', YX[,    factor(as.character(stage), levels=X.msm[, levels(stage)])])
		risk.df			<- data.table(risk='stage',factor=X.seq[, levels(stage)])		
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
									z	<- table( X.seq[, risk, with=FALSE], useNA='ifany')
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
		#
		# 	compute cens.table.bs	
		# 	for every bootstrap cdelta.bs
		#		compute tperiod.bs by	tperiod.start/end - cdelta.bs
		#		get table
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
		gc()
		#
		# 	end boostrapping
		#
		# 	revert factor to tp for X.msm
		set(X.msm, NULL, 'stage', X.msm[[risktp.col]])
		set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])									
		#	save
		ans$cens.table		<- cens.table
		ans$cens.table.bs	<- cens.table.bs
		if(!is.na(save.file))
		{			
			cat(paste('\nsave to file', save.file))
			save(ans, file=save.file)
		}
	}
	ans
}
######################################################################################
censoring.get.all.tables.by.trinterval<- function(YX=NULL, X.seq=NULL, X.msm=NULL, X.clu=NULL, tperiod.info=NULL, resume=TRUE, save.file=NA, method=NA, risk.col=NA, risktp.col=NA, factor.ref.v=NA, c.period=0.125, bs.n=100, bs.cdelta.min=2, bs.cdelta.max=3)
{
	stopifnot(!is.na(risk.col), !is.na(risktp.col), !is.na(factor.ref.v))
	t.period	
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
		if(is.null(YX) || is.null(X.seq) || is.null(X.msm) || is.null(X.clu))
		{
			cat(paste('\nreturn NULL table', method))
			return(ans)	
		}			
		cat(paste('\ntables by method', method,'\nrisk.col', risk.col,'\nrisk.col.tp', risktp.col,'\nfactor.ref.v',factor.ref.v))
		#YX<- copy(YX.s); X.clu<- copy(X.clu.s); X.seq<- copy(X.seq.s); X.msm<- copy(X.msm.s)
		#	set stratifications. always use 'risktp.col' for censoring
		set(YX, NULL, 'stage', YX[[risktp.col]])				
		set(X.clu, NULL, 'stage', X.clu[[risktp.col]])
		set(X.seq, NULL, 'stage', X.seq[[risktp.col]])
		set(X.msm, NULL, 'stage', X.msm[[risktp.col]])
		set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])									
		set(X.seq, NULL, 'stage', X.seq[, factor(as.character(stage), levels=X.msm[, levels(stage)])])		
		set(X.clu, NULL, 'stage', X.clu[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
		set(YX,    NULL, 'stage', YX[,    factor(as.character(stage), levels=X.msm[, levels(stage)])])
		risk.df			<- data.table(risk='stage',factor=X.seq[, levels(stage)])		
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
									z	<- table( X.seq[, risk, with=FALSE], useNA='ifany')
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
		tmp				<- rev(seq( tperiod.info[, max(t.period.max)], tperiod.info[, min(t.period.min)]-c.period, by=-c.period ))
		tperiod.bs		<- data.table(t.period=factor(seq_len(length(tmp)-1)), t.period.min=tmp[-length(tmp)], t.period.max=tmp[-1])
		tmp				<- tperiod.bs[, max(t.period.max)]
		tperiod.bs		<- lapply(seq_along(bs.cdelta), function(b)
				{
					tmp	<- tperiod.bs[, list(t.period.min.bs= t.period.min-bs.cdelta[b], t.period.max.bs= t.period.max-bs.cdelta[b], cens.t=tmp-bs.cdelta[b], cens.delta=bs.cdelta[b]), by='t.period']
					tmp[, BS:=b]
					tmp
				})
		tperiod.bs		<- do.call('rbind', tperiod.bs)		
		#	expand tperiod.bs to include risk and factor
		tmp				<- as.data.table(expand.grid(risk='stage', factor=X.msm[, levels(stage)], t.period= tperiod.bs[, unique(t.period)], stringsAsFactors=F))
		tperiod.bs		<- merge(tperiod.bs, tmp, by='t.period', allow.cartesian=TRUE)
		#
		# 	compute cens.table.bs	
		# 	for every bootstrap cdelta.bs
		#		compute tperiod.bs by	tperiod.start/end - cdelta.bs
		#		get table
		cat(paste('\ncompute cens.table.bs'))
		cens.table.bs	<- tperiod.bs[, {
					bs.breaks	<- c( -Inf, t.period.min.bs, max(t.period.max.bs), Inf)
					bs.labels	<- paste(factor,seq.int(0,length(bs.breaks)-2),sep='.')
					tmp			<- which( X.msm[[risk]]==factor )							
					n.all		<- cut( X.msm[tmp, AnyPos_T1], breaks=bs.breaks, labels=bs.labels, right=FALSE)
					n.all		<- table(n.all, useNA='ifany'  )					
					n.pseudocens<- cut( subset(X.msm[tmp,], t.AnyPos_T1<cens.t[1])[, AnyPos_T1], breaks=bs.breaks, labels=bs.labels, right=FALSE)
					n.pseudocens<- table(n.pseudocens, useNA='ifany'  )
					tmp			<- merge( data.table( factor.tp=rownames(n.all), n=as.numeric(unclass(n.all)) ), data.table( factor.tp=rownames(n.pseudocens), nc=as.numeric(unclass(n.pseudocens)) ), by='factor.tp', all.x=TRUE, all.y=TRUE )
					set(tmp, NULL, 't.period', tmp[, as.numeric(regmatches(factor.tp, regexpr('[0-9]+$',factor.tp)))])
					setkey(tmp, t.period)
					tmp[, t.period.min.bs:=c(NA_real_,t.period.min.bs,NA_real_)]
					tmp[, t.period.max.bs:=c(NA_real_,t.period.max.bs,NA_real_)]
					list(factor.tp=tmp[, factor.tp], n=tmp[, n], nc=tmp[, nc], t.period.min.bs=tmp[, t.period.min.bs], t.period.max.bs=tmp[,t.period.max.bs], cens.delta=cens.delta[1], cens.t=cens.t[1])
				}, by=c('BS','risk','factor')]	
		cens.table.bs	<- subset(cens.table.bs, !is.na(t.period.min.bs))
		setnames(cens.table.bs, c('factor.tp','factor'), c('factor', 'factor2'))			
		set(cens.table.bs, NULL, 'stat',  cens.table.bs[, paste('X.msm.cbs',BS,sep='')])
		cens.table.bs[, t.period:= cens.table.bs[, regmatches(factor,regexpr('[0-9]+$',factor))]]
		gc()
		#
		# 	end boostrapping
		#
		# 	revert factor to tp for X.msm
		set(X.msm, NULL, 'stage', X.msm[[risktp.col]])
		set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])									
		#	save
		ans$cens.table		<- cens.table
		ans$cens.table.bs	<- cens.table.bs
		if(!is.na(save.file))
		{			
			cat(paste('\nsave to file', save.file))
			save(ans, file=save.file)
		}
	}
	ans
}
######################################################################################
sampling.get.all.tables<- function(YX=NULL, X.seq=NULL, X.msm=NULL, X.clu=NULL, tperiod.info=NULL, resume=TRUE, save.file=NA, method=NA, risk.col=NA, risktp.col=NA, factor.ref.v=NA, bs.n=100, bs.cdelta.min=2, bs.cdelta.max=3)
{	
	stopifnot(!is.na(risk.col), !is.na(risktp.col), !is.na(factor.ref.v))
	
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
		if(is.null(YX) || is.null(X.seq) || is.null(X.msm) || is.null(X.clu))
		{
			cat(paste('\nreturn NULL table', method))
			return(ans)	
		}			
		cat(paste('\ntables by method', method,'\nrisk.col', risk.col,'\nrisk.col.tp', risktp.col,'\nfactor.ref.v',factor.ref.v))
		#	YX<- copy(YX.s); X.clu<- copy(X.clu.s); X.seq<- copy(X.seq.s); X.msm<- copy(X.msm.s)
		#	get number of recipients per time period			
		ans$cens.Patient.n	<- do.call('rbind',list( 	
						YX[, list(stat='YX', Patient.n=length(unique(Patient))), by='t.period'],
						X.clu[, list(stat='X.clu', Patient.n=length(unique(Patient))), by='t.period'],
						X.seq[, list(stat='X.seq', Patient.n=length(unique(Patient))), by='t.period'],
						X.msm[, list(stat='X.msm', Patient.n=length(unique(Patient))), by='t.period']	))
		#	get diagnosis times of transmitters for each stage
		set(X.msm, NULL, 'stage', X.msm[[risktp.col]])
		set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])			
		risk.df			<- data.table(risk='stage',factor=X.msm[, levels(stage)])
		ans$cens.AnyPos_T1	<- risk.df[,	{
					tmp	<- subset( X.msm[ which(X.msm[[risk]]==factor), ], select=c(t.Patient, t.AnyPos_T1))
					setkey(tmp, t.Patient)			
					tmp	<- unique(tmp)
					tmp[[risk]]	<- factor[1]
					tmp
				}, by=c('risk','factor')]		
		#	if 'tp' is '', reset 'stage' so we use the global risk.col instead of the one stratified by t.period
		tp				<- regmatches(method, regexpr('tp[0-9]', method))
		tp				<- ifelse(length(tp), paste('.',substr(tp, 3, 3),sep=''), '')		
		tmp				<- ifelse(tp=='', risk.col, risktp.col)
		#	set stratifications to stage
		set(YX, NULL, 'stage', YX[[tmp]])				
		set(X.clu, NULL, 'stage', X.clu[[tmp]])
		set(X.seq, NULL, 'stage', X.seq[[tmp]])
		set(X.msm, NULL, 'stage', X.msm[[tmp]])
		set(X.msm, NULL, 'stage', X.msm[, factor(as.character(stage))])									
		set(X.seq, NULL, 'stage', X.seq[, factor(as.character(stage), levels=X.msm[, levels(stage)])])								
		set(X.clu, NULL, 'stage', X.clu[, factor(as.character(stage), levels=X.msm[, levels(stage)])])
		set(YX,    NULL, 'stage', YX[,    factor(as.character(stage), levels=X.msm[, levels(stage)])])			
		risk.df			<- data.table(risk='stage',factor=X.seq[, levels(stage)], risk.ref='stage', factor.ref=factor.ref.v)
		gc()
		cat(paste('\ncompute nt.table'))
		#	potential transmission intervals by stage and transmitter
		nt.table.pt		<- do.call('rbind', list( 	YX[, list(nt= length(t), stat='YX'), by=c('stage','t.Patient')],
						X.clu[, list(nt= length(t), stat='X.clu'), by=c('stage','t.Patient')],
						X.seq[, list(nt= length(t), stat='X.seq'), by=c('stage','t.Patient')],
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
						X.seq[, list(nt= length(t), stat='X.seq'), by=c('stage','Patient')],
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
									z	<- table( X.seq[, risk, with=FALSE])
									list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.seq')												
								},by='risk'],
						risk.df[,	{
									z	<- table( X.msm[, risk, with=FALSE])
									list(factor=rownames(z), n=as.numeric(unclass(z)), stat='X.msm')
								},by='risk']))
		risk.table		<- merge(risk.table, risk.table[, list(factor=factor, p= n/sum(n, na.rm=TRUE)), by=c('stat','risk')], by=c('stat','risk','factor'))
		ans$risk.table	<- risk.table	
		#	
		if(!is.na(save.file))
		{			
			cat(paste('\nsave to file', save.file))
			save(ans, file=save.file)
		}
	}
	ans
}
######################################################################################
age.precompute<- function(	indir, indircov, infile.cov.study, infile.viro.study, infile.immu.study, infile.treatment.study, infile.cov.all, infile.viro.all, infile.immu.all, infile.treatment.all, infile.trm.model,
		clu.indir, clu.insignat, clu.infile,
		infile, infiletree, insignat, clu.infilexml.opt, clu.infilexml.template,
		method, method.recentctime, method.nodectime, method.risk, method.Acute, method.minQLowerU, method.use.AcuteSpec, method.brl.bwhost, method.lRNA.supp, method.thresh.pcoal, method.minLowerUWithNegT, method.tpcut, method.PDT, method.cut.brl, tp.cut, adjust.AcuteByNegT, any.pos.grace.yr, dur.Acute, method.thresh.bs, 
		outdir, outfile,
		t.period, t.recent.startctime, t.endctime, t.recent.endctime,
		X.tables, resume, verbose
)
{
	#
	#	get data relating to study population (subtype B sequ)
	#
	tmp				<- project.athena.Fisheretal.select.denominator(	indir, infile, insignat, indircov, infile.cov.study, infile.viro.study, infile.immu.study, infile.treatment.study, infiletree=infiletree, 
			adjust.AcuteByNegT=adjust.AcuteByNegT, adjust.NegT4Acute=NA, adjust.NegTByDetectability=1/12, adjust.minSCwindow=0.25, adjust.AcuteSelect=c('Yes','Maybe'), use.AcuteSpec=method.use.AcuteSpec, 
			thresh.bs=method.thresh.bs, t.recent.endctime=t.recent.endctime, t.recent.startctime=t.recent.startctime)	
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
			infiletree=NULL, adjust.AcuteByNegT=adjust.AcuteByNegT, adjust.NegT4Acute=NA, adjust.NegTByDetectability=1/12, adjust.minSCwindow=0.25, adjust.AcuteSelect=c('Yes','Maybe'), use.AcuteSpec=method.use.AcuteSpec,
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
	#	get rough idea about (backward) time to infection from time to diagnosis, taking midpoint of SC interval as 'training data'
	#
	plot.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 't2inf',method.PDT,method,sep='')
	tmp				<- project.athena.Fisheretal.t2inf(	df.all.allmsm,
			method.Acute=method.Acute, method.minQLowerU=method.minQLowerU,
			adjust.AcuteByNegT=0.75, adjust.dt.CD4=1, adjust.AnyPos_y=2003, adjust.NegT=2, dur.AcuteYes=dur.Acute['Yes'], dur.AcuteMaybe=dur.Acute['Maybe'], use.AcuteSpec=method.use.AcuteSpec, t.recent.endctime=t.recent.endctime, 
			plot.file=plot.file)
	predict.t2inf	<- tmp$predict.t2inf
	t2inf.args		<- tmp$t2inf.args	
	#	determine best quantile parameter
	if(is.na(method.minQLowerU))
	{
		plot.file			<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 't2infq_',method.PDT,method,sep='')
		method.minQLowerU	<- project.athena.Fisheretal.t2inf.estimate.quantileparameter(df.all.allmsm, predict.t2inf, t2inf.args, t.recent.endctime, plot.file=plot.file)		
	}
	#
	#	select potential transmitters within same cluster
	#
	df.tpairs		<- project.athena.Fisheretal.select.transmitters.by.B4WindowAnyPos(clumsm.info, df.denom.CLU, any.pos.grace.yr= any.pos.grace.yr, select.if.transmitter.seq.unique=FALSE)		
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
	save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RICT',method.PDT,'_',method,'_tATHENAclu','.R',sep='')	
	YX.part1		<- project.athena.Fisheretal.YX.part1(df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, indircov=indircov, ri=NULL, df.tpairs=df.tpairs, tperiod.info=NULL, lRNA.supp=method.lRNA.supp, method.minLowerUWithNegT=method.minLowerUWithNegT, t.period=t.period, t.endctime=t.endctime, save.file=save.file, resume=resume)
	YX.part1		<- merge( YX.part1, subset( df.tpairs, select=c(FASTASampleCode, t.FASTASampleCode, cluster) ), by=c('FASTASampleCode','t.FASTASampleCode'), all.x=1)
	YX.part1[, class:='pt']
	gc()	
	rm.zero.score	<- TRUE	
	save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'YX',method.PDT,method,'.R',sep='')
	save.all		<- FALSE
	#save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'YX',method.PDT,method,'_all.R',sep='')
	#save.all		<- TRUE
	# df.tpairs.4.rawbrl=df.tpairs; thresh.pcoal=method.thresh.pcoal; brl.bwhost.multiplier=method.brl.bwhost; method.minLowerUWithNegT=method.minLowerUWithNegT; lRNA.supp=method.lRNA.supp; infilecov= infile.cov.study
	tmp				<- project.athena.Fisheretal.YX.part2(	YX.part1, df.all, df.treatment, df.viro, predict.t2inf, t2inf.args, indir, insignat, indircov, infile.cov.study, infiletree, infile.trm.model, outdir, outfile, cluphy=cluphy, cluphy.info=cluphy.info, cluphy.map.nodectime=cluphy.map.nodectime, df.tpairs.4.rawbrl=df.tpairs, dur.Acute=dur.Acute,
			rm.zero.score=rm.zero.score, any.pos.grace.yr=any.pos.grace.yr, thresh.pcoal=method.thresh.pcoal, cut.brl=method.cut.brl, brl.bwhost.multiplier=method.brl.bwhost, method.minLowerUWithNegT=method.minLowerUWithNegT, lRNA.supp=method.lRNA.supp, tp.cut=tp.cut,
			t.period=t.period, save.file=save.file, resume=resume, method=method, save.all=save.all)
	YX				<- copy(tmp$YX)
	Y.brl.bs		<- copy(tmp$Y.brl.bs)
	tmp				<- NULL
	gc()
	tperiod.info	<- merge(df.all, unique( subset(YX, select=c(Patient, t.period)) ), by='Patient')
	tperiod.info	<- tperiod.info[, list(t.period.min=min(AnyPos_T1)), by='t.period']
	set(tperiod.info, 1L,'t.period.min',t.recent.startctime-0.1)	
	tperiod.info[, t.period.max:=c(tperiod.info[-1, t.period.min], t.recent.endctime)]
	ri.PT			<- subset(YX, score.Y>0, select=c(Patient, t))[, list(n.t.infw= length(unique(t))), by='Patient']		
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
		X.clu			<- project.athena.Fisheretal.YX.part1(	tmp, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, indircov=indircov, ri=ri.CLU, df.tpairs=NULL, 
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
			X.seq			<- project.athena.Fisheretal.YX.part1(	df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, indircov=indircov, ri=ri.CLU, df.tpairs=NULL, 
					tperiod.info=tperiod.info, lRNA.supp=method.lRNA.supp, t.period=t.period, t.endctime=t.endctime, method.minLowerUWithNegT=method.minLowerUWithNegT, save.file=save.file, resume=resume)
		}
		if(method.PDT=='SEQ')
		{						
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RISEQ_',method,'_tATHENAseq','.R',sep='')
			X.seq			<- project.athena.Fisheretal.YX.part1(	df.all, df.immu, df.viro, df.treatment, predict.t2inf, t2inf.args, indircov=indircov, ri=ri.SEQ, df.tpairs=NULL, 
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
			X.msm			<- project.athena.Fisheretal.YX.part1(	df.all.allmsm, df.immu.allmsm, df.viro.allmsm, df.treatment.allmsm, predict.t2inf, t2inf.args, indircov=indircov, ri=ri.CLU, df.tpairs=NULL, 
					tperiod.info=tperiod.info, lRNA.supp=method.lRNA.supp, t.period=t.period, t.endctime=t.endctime, method.minLowerUWithNegT=method.minLowerUWithNegT, save.file=save.file, resume=resume)
		}
		if(method.PDT=='SEQ')
		{												
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'RISEQ_',method,'_tATHENAmsm','.R',sep='')
			X.msm			<- project.athena.Fisheretal.YX.part1(	df.all.allmsm, df.immu.allmsm, df.viro.allmsm, df.treatment.allmsm, predict.t2inf, t2inf.args, indircov=indircov, ri=ri.SEQ, df.tpairs=NULL, 
					tperiod.info=tperiod.info, lRNA.supp=method.lRNA.supp, t.period=t.period, t.endctime=t.endctime, method.minLowerUWithNegT=method.minLowerUWithNegT, save.file=save.file, resume=resume)
		}
		gc()		
	}
	#
	#X.clu	<- X.clu[sample(seq_len(nrow(X.clu)), 1e6),]
	#X.seq	<- X.seq[sample(seq_len(nrow(X.seq)), 2e6),]
	#X.msm	<- X.msm[sample(seq_len(nrow(X.msm)), 3e6),]
	#	stratify YX
	if(grepl('m5A',method.risk))	
	{
		YX					<- stratificationmodel.Age_253045.Stage_UA_UC_D_T_F(YX)
		if(is.null(X.tables))
		{
			X.clu			<- stratificationmodel.Age_253045.Stage_UA_UC_D_T_F(X.clu)
			gc()
			X.seq			<- stratificationmodel.Age_253045.Stage_UA_UC_D_T_F(X.seq)
			gc()
			X.msm			<- stratificationmodel.Age_253045.Stage_UA_UC_D_T_F(X.msm)
			gc()
		}
	}
	if(grepl('m5B',method.risk))	
	{
		YX					<- stratificationmodel.Age_253045.Stage_UAE_UAC_UC_DAC_D_T_F(YX)
		if(is.null(X.tables))
		{
			X.clu			<- stratificationmodel.Age_253045.Stage_UAE_UAC_UC_DAC_D_T_F(X.clu)
			gc()
			X.seq			<- stratificationmodel.Age_253045.Stage_UAE_UAC_UC_DAC_D_T_F(X.seq)
			gc()
			X.msm			<- stratificationmodel.Age_253045.Stage_UAE_UAC_UC_DAC_D_T_F(X.msm)
			gc()
		}
	}
	save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', method, 'STRAT_',gsub('\\.clu\\.adj','',gsub('\\.tp[0-9]','',method.risk)),'.R',sep='')
	save(YX, X.clu, X.seq, X.msm, file=save.file)	
#STOP1		
#stop()
	#
	#	compute sampling and censoring tables that are needed for adjustments
	#
	if(grepl('adj',method.risk) & grepl('clu',method.risk))
	{
			#	get args
			save.file		<- tmp	<- NA
			resume			<- 0
			args			<- sampling.get.all.tables.args(method.risk)			
			if(grepl('m5A',method.risk))		
				tmp			<- 'm5A'
			if(grepl('m5B',method.risk))		
				tmp			<- 'm5B'
			if(is.na(tmp))	
				stop('unknown method.risk')				
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore', method,'_Stables',method.PDT,'_',tmp,'.R',sep='')
			#	sampling tables
			Stab			<- sampling.get.all.tables(YX, X.seq, X.msm, X.clu, tperiod.info=tperiod.info, resume=resume, save.file=save.file, method=method.risk, risk.col=args['risk.col'], risktp.col=args['risktp.col'], factor.ref.v=args['factor.ref.v'])
			#	censoring tables
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore', method,'_Ctables',method.PDT,'_',tmp,'.R',sep='')			
			Ctab			<- censoring.get.all.tables.by.trinterval(YX, X.seq, X.msm, X.clu, tperiod.info=tperiod.info, resume=resume, save.file=save.file, method=method.risk, risk.col=args['risk.col'], risktp.col=args['risktp.col'], factor.ref.v=args['factor.ref.v'], c.period=t.period, bs.n=1000, bs.cdelta.min=2, bs.cdelta.max=3)
			#Ctab			<- censoring.get.all.tables(YX, X.seq, X.msm, X.clu, tperiod.info=tperiod.info, resume=resume, save.file=save.file, method=method.risk, risk.col=args['risk.col'], risktp.col=args['risktp.col'], factor.ref.v=args['factor.ref.v'], bs.n=100, bs.cdelta.min=2, bs.cdelta.max=3)
stop()
	}
	X.clu<- X.seq<- X.msm<- NULL
	gc()
stop()	
	ans		<- list(predict.t2inf=predict.t2inf, t2inf.args=t2inf.args, df.all=df.all, YX=YX, Y.brl.bs=Y.brl.bs)
	ans
}

######################################################################################
stratificationmodel.Age_253045.Stage_UA_UC_D_T_F<- function(YX.m5)
{
	#	get age stages
	#YX.m5	<- copy(YX)
	if('U.score'%in%colnames(YX.m5))
		YX.m5[, U.score:=NULL]
	cat(paste('\nsubset to save mem\n'))
	if('score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, FASTASampleCode, t.FASTASampleCode, score.Y, telapsed, brl, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, w, w.i, w.in, w.t, w.tn, t.Age, Age ))	
	if(!'score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, t.Age, Age ))
	gc()	
	if('score.Y'%in%colnames(YX.m5) && YX.m5[, !any(score.Y>1.1)])
	{
		tmp	<- YX.m5[, score.Y>0]
		set(YX.m5, which(tmp), 'score.Y', YX.m5[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )
	}	
	#
	#	PREPARE AGE (Age and t.Age are ages at midpoint of infection interval)
	#
	#	set age group of transmitter				
	age.breaks		<- c(-1, 25, 30, 45, 100)
	age.labels		<- c('<25','<30','<45','<100')	
	YX.m5[, t.AgeC:= YX.m5[, cut(t.Age, breaks=age.breaks)]]	
	set(YX.m5, NULL, 't.AgeC', YX.m5[, factor(as.character(t.AgeC))])
	#	set age group of recipient
	YX.m5[, AgeC:= YX.m5[, cut(Age, breaks=age.breaks)]]	
	set(YX.m5, NULL, 'AgeC', YX.m5[, factor(as.character(AgeC))])
	#	set age group of pair
	YX.m5[, p.AgeC:= YX.m5[, paste(as.character(t.AgeC),'->',as.character(AgeC),sep='')]]
	set(YX.m5, NULL, 'p.AgeC', YX.m5[, factor(p.AgeC)])
	#
	#	PREPARE STAGE
	#	
	YX.m5[, stageC:=NA_character_]
	set(YX.m5, YX.m5[, which(stage=='Diag')], 'stageC', 'D')
	set(YX.m5, YX.m5[, which(stage=='ART.started')], 'stageC', 'T' )
	set(YX.m5, YX.m5[, which(stage=='U')], 'stageC', 'UC')
	#	confirmed acute --> UA
	set(YX.m5, YX.m5[, which(stage=='U' & t.isAcute=='Yes')], 'stageC', 'UA')
	#	confirmed acute and diagnosed in first 3 mnths --> UA
	set(YX.m5, YX.m5[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stageC', 'UA' )
	#	estimated recent
	set(YX.m5, YX.m5[, which(stage=='U' & (t-t.InfT)<1 & (t.isAcute!='Yes' | is.na(t.isAcute)))], 'stageC', 'UA')
	#	set Lost
	stopifnot( nrow(subset(YX.m5, stage!='U' & is.na(contact)))==0 )	
	set(YX.m5, YX.m5[, which(stage!='U' & contact!='Yes')], 'stageC', 'L')
	set(YX.m5, NULL, 'stageC', YX.m5[, factor(stageC)])
	#
	#	PREPARE STAGE:AGE
	#	
	YX.m5[, t.stAgeC:= YX.m5[, factor(paste(as.character(stageC),'_',as.character(t.AgeC),sep=''))]]
	#
	#	PREPARE cross with time period
	#	
	if(1 & 't.period'%in%colnames(YX.m5) & YX.m5[, any(t.period=='6')])
	{
		set(YX.m5, YX.m5[, which(t.period%in%c('5','6'))], 't.period', '4')
		set(YX.m5, NULL, 't.period', YX.m5[, factor(t.period)])
	}
	if('t.period'%in%colnames(YX.m5))
	{
		
		cat(paste('\nadding tperiod columns\n'))
		YX.m5[, AgeC.prd:= factor(paste(as.character(AgeC), t.period,sep='.'))]
		YX.m5[, t.AgeC.prd:= factor(paste(as.character(t.AgeC), t.period,sep='.'))]
		YX.m5[, p.AgeC.prd:= factor(paste(as.character(t.AgeC.prd),'->',as.character(AgeC.prd),sep=''))]
		YX.m5[, t.stAgeC.prd:= factor(paste(as.character(t.stAgeC), t.period,sep='.'))]		
	}	
	gc()	
	cat(paste('\nsubset to save further mem\n'))
	set(YX.m5, NULL, c('t.InfT','t.isAcute','contact','stage','t.AnyT_T1'), NULL)
	gc()
	YX.m5
}
######################################################################################
stratificationmodel.Age_253045.Stage_UAE_UAC_UC_DAC_D_T_F<- function(YX.m5)
{
	#	get age stages
	#	YX.m5	<- copy(YX)
	if('U.score'%in%colnames(YX.m5))
		YX.m5[, U.score:=NULL]
	cat(paste('\nsubset to save mem\n'))
	if('score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, FASTASampleCode, t.FASTASampleCode, score.Y, telapsed, brl, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, w, w.i, w.in, w.t, w.tn, t.Age, Age ))	
	if(!'score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, t.Age, Age ))
	gc()	
	if('score.Y'%in%colnames(YX.m5) && YX.m5[, !any(score.Y>1.1)])
	{
		tmp	<- YX.m5[, score.Y>0]
		set(YX.m5, which(tmp), 'score.Y', YX.m5[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )
	}	
	#
	#	PREPARE AGE (Age and t.Age are ages at midpoint of infection interval)
	#
	#	set age group of transmitter				
	age.breaks		<- c(-1, 25, 30, 45, 100)
	age.labels		<- c('<25','<30','<45','<100')	
	YX.m5[, t.AgeC:= YX.m5[, cut(t.Age, breaks=age.breaks)]]	
	set(YX.m5, NULL, 't.AgeC', YX.m5[, factor(as.character(t.AgeC))])
	#	set age group of recipient
	YX.m5[, AgeC:= YX.m5[, cut(Age, breaks=age.breaks)]]	
	set(YX.m5, NULL, 'AgeC', YX.m5[, factor(as.character(AgeC))])
	#	set age group of pair
	YX.m5[, p.AgeC:= YX.m5[, paste(as.character(t.AgeC),'->',as.character(AgeC),sep='')]]
	set(YX.m5, NULL, 'p.AgeC', YX.m5[, factor(p.AgeC)])
	#
	#	PREPARE STAGE
	#	
	YX.m5[, stageC:=NA_character_]
	set(YX.m5, YX.m5[, which(stage=='Diag')], 'stageC', 'D')
	set(YX.m5, YX.m5[, which(stage=='ART.started')], 'stageC', 'T' )
	set(YX.m5, YX.m5[, which(stage=='U')], 'stageC', 'UC')
	#	confirmed acute --> UA
	set(YX.m5, YX.m5[, which(stage=='U' & t.isAcute=='Yes')], 'stageC', 'UAC')
	#	confirmed acute and diagnosed in first 3 mnths --> UA
	set(YX.m5, YX.m5[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stageC', 'DAC' )
	#	estimated recent
	set(YX.m5, YX.m5[, which(stage=='U' & (t-t.InfT)<1 & (t.isAcute!='Yes' | is.na(t.isAcute)))], 'stageC', 'UAE')
	#	set Lost
	stopifnot( nrow(subset(YX.m5, stage!='U' & is.na(contact)))==0 )	
	set(YX.m5, YX.m5[, which(stage!='U' & contact!='Yes')], 'stageC', 'L')
	set(YX.m5, NULL, 'stageC', YX.m5[, factor(stageC)])
	#
	#	PREPARE STAGE:AGE
	#	
	YX.m5[, t.stAgeC:= YX.m5[, factor(paste(as.character(stageC),'_',as.character(t.AgeC),sep=''))]]
	#
	#	PREPARE cross with time period
	#	
	if(1 & 't.period'%in%colnames(YX.m5) & YX.m5[, any(t.period=='6')])
	{
		set(YX.m5, YX.m5[, which(t.period%in%c('5','6'))], 't.period', '4')
		set(YX.m5, NULL, 't.period', YX.m5[, factor(t.period)])
	}
	if('t.period'%in%colnames(YX.m5))
	{
		
		cat(paste('\nadding tperiod columns\n'))
		YX.m5[, AgeC.prd:= factor(paste(as.character(AgeC), t.period,sep='.'))]
		YX.m5[, t.AgeC.prd:= factor(paste(as.character(t.AgeC), t.period,sep='.'))]
		YX.m5[, p.AgeC.prd:= factor(paste(as.character(t.AgeC.prd),'->',as.character(AgeC.prd),sep=''))]
		YX.m5[, t.stAgeC.prd:= factor(paste(as.character(t.stAgeC), t.period,sep='.'))]		
	}	
	gc()	
	cat(paste('\nsubset to save further mem\n'))
	set(YX.m5, NULL, c('t.InfT','t.isAcute','contact','stage','t.AnyT_T1'), NULL)
	gc()
	YX.m5
}