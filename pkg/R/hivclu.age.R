######################################################################################
age.props_univariate<- function()
{	
	require(reshape2)
	require(data.table)	
	require(survival)
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
		method.risk				<- 'm5D.tp4'
		method.Acute			<- 'higher'	#'central'#'empirical'
		method.minQLowerU		<- 0.148
		method.use.AcuteSpec	<- 1
		method.brl.bwhost		<- 2
		method.lRNA.supp		<- 100
		#method.thresh.pcoal	<- 0.2
		method.thresh.pcoal		<- 0.1		# try lower to see if bias reduced
		method.minLowerUWithNegT<- 1
		method.cut.brl			<- Inf		#does not make a difference because compatibility test kills these anyway
		method.tpcut			<- 7
		method.PDT				<- 'SEQ'	# 'PDT'	
		#method.thresh.bs		<- 0.8
		method.thresh.bs		<- 0.7		# try lower to see if bias reduced
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
	#X.tables				<- age.get.Xtables(method, method.PDT, method.risk, outdir, outfile, insignat)
	X.tables				<- age.get.sampling.censoring.models(method, method.PDT, method.risk, outdir, outfile, insignat, load.bs.id=1)
	#	if no tables, precompute tables and stop (because mem intensive)
	#	otherwise return stratified YX data.table and continue
	tmp	<- age.precompute(	indir, indircov, infile.cov.study, infile.viro.study, infile.immu.study, infile.treatment.study, infile.cov.all, infile.viro.all, infile.immu.all, infile.treatment.all, infile.trm.model,
			clu.indir, clu.insignat, clu.infile,
			infile, infiletree, insignat, clu.infilexml.opt, clu.infilexml.template,
			method, method.recentctime, method.nodectime, method.risk, method.Acute, method.minQLowerU, method.use.AcuteSpec, method.brl.bwhost, method.lRNA.supp, method.thresh.pcoal, method.minLowerUWithNegT, method.tpcut, method.PDT, method.cut.brl, tp.cut, adjust.AcuteByNegT, any.pos.grace.yr, dur.Acute, method.thresh.bs, 
			outdir, outfile,
			t.period, t.recent.startctime, t.endctime, t.recent.endctime,
			is.null(X.tables$st) | is.null(X.tables$sm) | is.null(X.tables$cm), resume, verbose
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
age.explore.Recipient.bias<- function()
{
	YX
	cluphy.info
	ri.ALLMSM
	ri.CLU
	indir		<- '~/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C1V100bInfs0.7T7STRAT_m5D.R'
	outfile.pdf	<- paste(indir, '/', gsub('\\.R',paste('_Rec_by_AgeC_','160217','.pdf',sep=''),infile), sep='')
	outfile.pdf	<- paste(indir, '/', gsub('\\.R',paste('_Rec_by_AgeC_','160223_C1s70','.pdf',sep=''),infile), sep='')
	
	tmp3	<- YX[, list(tna=length(unique(t))), by=c('Patient','AgeC')]
	tmp3	<- merge(tmp3, YX[, list(AnyPos_T1=AnyPos_T1[1], tn=length(unique(t))), by=c('Patient')], by=c('Patient'))
	tmp3[, tmidC2:= cut(AnyPos_T1, breaks=c(-Inf, 2008, Inf), labels=c('2004-2007','2008-2010'))]
	tmp3	<- merge(tmp3,tmp3[, list(nRec=length(unique(Patient))), by='tmidC2'],by='tmidC2')
	set(tmp3, NULL, 'AgeC', tmp3[, factor(as.character(AgeC), levels=c('(-1,28]','(28,33]','(33,38]','(38,43]','(43,48]','(48,100]'), labels=c('16-27','28-32','33-37','38-42','43-47','48-80'))])
	YXar	<- copy(tmp3)
	setkey(YXar, tmidC2, Patient, AgeC)
	YXar	<- unique(subset(YXar, select=c(tmidC2, Patient, AgeC, tn, tna, nRec)))
	set(YXar, NULL, 'tna', YXar[, tna/tn])
	YXar	<- YXar[, list(TYPE='Recipient', N= sum(tna), Ntp=nRec[1]), by=c('tmidC2','AgeC')]
	#	add distribution among new diagnoses
	tmp		<- subset(df.all.allmsm, AnyPos_T1>2004 & AnyPos_T1<2011 & Trm%in%c('MSM','BI'), c(Patient, AnyPos_T1, DateBorn, PosSeqT))
	setkey(tmp, Patient)
	tmp		<- unique(tmp)
	tmp[, Age:=AnyPos_T1-DateBorn]
	tmp[, tmidC2:= cut(AnyPos_T1, breaks=c(-Inf, 2008, Inf), labels=c('2004-2007','2008-2010'))]
	tmp[, AgeC:= cut(Age, breaks=c(-Inf, 28, 33, 38, 43, 48, Inf), labels=c('16-27','28-32','33-37','38-42','43-47','48-80'))]
	tmp2	<- tmp[, list(TYPE='Diagnosed', N= length(Patient)), by=c('tmidC2','AgeC')]
	tmp2	<- merge(tmp2, tmp2[, list(Ntp=sum(N)),by='tmidC2'], by='tmidC2')
	YXar	<- rbind(YXar, tmp2)
	#	add distribution among recently infected
	tmp3	<- copy(ri.ALLMSM)	
	tmp3[, Recent:='Y']	
	tmp		<- merge(tmp, unique(subset(tmp3, select=c(Patient,Recent))), by='Patient', all.x=1)
	set(tmp, tmp[, which(is.na(Recent))],'Recent', 'N')
	tmp2	<- subset(tmp, Recent=='Y')[, list(TYPE='Recently infected', N= length(Patient)), by=c('tmidC2','AgeC')]
	tmp2	<- merge(tmp2, tmp2[, list(Ntp=sum(N)),by='tmidC2'], by='tmidC2')
	YXar	<- rbind(YXar, tmp2)
	#	add distribution among sequenced recently infected
	tmp2	<- subset(tmp, Recent=='Y' & !is.na(PosSeqT))[, list(TYPE='Recently Infected Sequenced', N= length(Patient)), by=c('tmidC2','AgeC')]	
	tmp2	<- merge(tmp2, tmp2[, list(Ntp=sum(N)),by='tmidC2'], by='tmidC2')
	YXar	<- rbind(YXar, tmp2)
	#	add distribution among clustering recently infected
	tmp3	<- copy(ri.CLU)	
	tmp3[, Clustering:='Y']	
	tmp		<- merge(tmp, tmp3, by='Patient', all.x=1)
	set(tmp, tmp[, which(is.na(Clustering))],'Clustering', 'N')
	tmp2	<- subset(tmp, Recent=='Y' & Clustering=='Y')[, list(TYPE='Recently Infected Clustering', N= length(Patient)), by=c('tmidC2','AgeC')]
	tmp2	<- merge(tmp2, tmp2[, list(Ntp=sum(N)),by='tmidC2'], by='tmidC2')
	YXar	<- rbind(YXar, tmp2)	
	YXar[, P:=N/Ntp]
	#	plot
	YXplot	<- copy(YXar)
	#set(YXplot, NULL, 'TYPE', YXplot[, factor(TYPE, levels=c('Recipient','Recently Infected Clustering','Recently Infected Sequenced','Recently infected','Diagnosed'), labels=c('Recipient MSM\n(of whom sources were characterised)','Recently infected\nIn phylo cluster\nBS=80% GD=Inf Coal=80%','Recently infected\nwith a sequence','Recently infected\n(population)','Newly diagnosed MSM\n(population)') )])
	set(YXplot, NULL, 'TYPE', YXplot[, factor(TYPE, levels=c('Recipient','Recently Infected Clustering','Recently Infected Sequenced','Recently infected','Diagnosed'), labels=c('Recipient MSM\n(of whom sources were characterised)','Recently infected\nIn phylo cluster\nBS=70% GD=Inf Coal=90%','Recently infected\nwith a sequence','Recently infected\n(population)','Newly diagnosed MSM\n(population)') )])
	set(YXplot, NULL, 'tmidC2', YXplot[, factor(as.character(tmidC2, levels=c('2004-2007','2008-2020'),labels=c('2004-2007\n(time of diagnosis)','2008-2010\n(time of diagnosis)')))])
	ggplot(YXplot, aes(y=100*P, x=AgeC, fill=TYPE)) +
			geom_bar(stat='identity',position='dodge', colour='black', width=0.8) + 
			facet_grid(~tmidC2) + 
			scale_y_continuous(expand=c(0,0), limit=c(0,25)) +
			#scale_fill_manual(values=c('Newly diagnosed MSM\n(population)'="#FEB24C", 'Recipient MSM\n(of whom sources could be characterised)'="#E31A1C")) +
			scale_fill_brewer(palette='Pastel1') +
			coord_flip() +
			theme_bw() + theme(legend.position='bottom') +
			labs(x='Age group\n(years)\n', y='%', fill='', title='Study population\n')			
	ggsave(file=outfile.pdf, w=8.5, h=6)
	#
	#	superspreaders
	#
	tmp2	<- subset(YX, stageC=='UAC')[, list(AnyPos_T1=mean(AnyPos_T1), iR=length(unique(Patient))), by=c('t.Patient','t.AnyPos_T1','t.AgeC')]
	set(tmp2, NULL, 't.AgeC', tmp2[, factor(as.character(t.AgeC), levels=c('(-1,28]','(28,33]','(33,38]','(38,43]','(43,48]','(48,100]'), labels=c('16-27','28-32','33-37','38-42','43-47','48-80'))])
	tmp2[, tmidC2:= cut(AnyPos_T1, breaks=c(-Inf, 2008, Inf), labels=c('2004-2007','2008-2010'))]
	tmp2	<- tmp2[, list(AverageNoRecipientsFromAcuteTransmitters=mean(iR)), by=c('tmidC2','t.AgeC')]
	setkey(tmp2, tmidC2, t.AgeC)
}
######################################################################################
age.get.sampling.censoring.models<- function(method, method.PDT, method.risk, outdir, outfile, insignat, load.bs.id=100)
{
	ans				<- vector('list')
	tmp				<- NA
	if(grepl('m5A',method.risk))	tmp	<- 'm5A'
	if(grepl('m5B',method.risk))	tmp	<- 'm5B'
	if(grepl('m5C',method.risk))	tmp	<- 'm5C'
	if(grepl('m5D',method.risk))	tmp	<- 'm5D'
	if(grepl('m5E',method.risk))	tmp	<- 'm5E'
	if(grepl('m5F',method.risk))	tmp	<- 'm5F'
	if(is.na(tmp))	stop('unknown method.risk')	
	
	save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore', method,'_Stables',method.PDT,'_',tmp,'.R',sep='')			
	ans$st			<- sampling.get.all.tables(NULL, NULL, NULL, NULL, tperiod.info=NULL, resume=TRUE, save.file=save.file, method=NA, risk.col=NA, risktp.col=NA, factor.ref.v=NA)
	save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore', method,'_Smodel',method.PDT,'_',tmp,'.R',sep='')
	ans$sm			<- sampling.model.calculate(NULL, NULL, NULL, NULL, NULL, resume=TRUE, save.file=save.file)		
	save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore', method,'_Cmodel',method.PDT,'_',tmp,'.R',sep='')
	ans$cm			<- censoring.model.calculate.bs(NULL, NULL, resume=TRUE, save.file=save.file, t.recent.endctime=NA, risk.col=NA, c.period=NA, c.smpl.n=NA, bs.n=load.bs.id, bs.cdelta.min=NA, bs.cdelta.max=NA)		
	ans
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
		if(grepl('m5C',method.risk))	save.file	<- 'm5C'
		if(grepl('m5D',method.risk))	save.file	<- 'm5D'
		if(grepl('m5E',method.risk))	save.file	<- 'm5E'
		if(grepl('m5F',method.risk))	save.file	<- 'm5F'
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
censoring.get.dataset<- function(tp.df, df.all.allmsm)
{	
	#	set up data set	
	if(!any(colnames(tp.df)=='r.Patient'))
		setnames(tp.df, c('t.Patient','Patient'), c('Patient','r.Patient'))
	setkey(tp.df, r.Patient, Patient)
	setkey(df.all.allmsm, Patient)
	tmp			<- unique(df.all.allmsm)
	tmp			<- subset( tmp, select=c(Patient, RegionHospital, DateBorn, isAcute, NegT, AnyPos_T1, PosCD4_T1, CD4_T1))
	#	add age at diagnosis
	tmp[, Age_T1:= AnyPos_T1-DateBorn]
	tmp[, DateBorn:=NULL]
	#	add isAcute
	set(tmp, tmp[, which(isAcute!='Yes' | is.na(isAcute))], 'isAcute', 'NI')
	set(tmp, NULL, 'isAcute', tmp[, factor(isAcute)])
	#	add CD4C
	tmp[, CD4C:= tmp[, cut(CD4_T1, breaks=c(-Inf,250, 350, 500, 750, Inf))]]
	set(tmp, tmp[, which(is.na(CD4C))],'CD4C','CD4.NA')
	set(tmp, tmp[, which(PosCD4_T1-AnyPos_T1>1)],'CD4C','CD4.late')	
	tmp[, AgeC_T1:= tmp[, cut(Age_T1, breaks=c(-Inf,25, 30, 45, Inf))]]	
	#	add ACD4C isAcute=='NI' and CD4C
	tmp[, ACD4C:= tmp[, cut(CD4_T1, breaks=c(-Inf,250, 350, 500, 750, Inf))]]
	set(tmp, tmp[, which(is.na(ACD4C))],'ACD4C','CD4.NA')
	set(tmp, tmp[, which(PosCD4_T1-AnyPos_T1>1)],'ACD4C','CD4.late')	
	set(tmp, tmp[, which(isAcute=='Yes')],'ACD4C','CD4.Acute')
	#	add age at diagnosis
	tmp[, AgeC_T1:= tmp[, cut(Age_T1, breaks=c(-Inf,25, 30, 45, Inf))]]
	#	interaction age and ACD4C
	tmp[, AACD4C:= tmp[, factor(paste(as.character(ACD4C),'-',as.character(AgeC_T1),sep=''))]]
	#
	tpd.df		<- merge(tp.df, tmp, by='Patient')	
	tpd.df
}
######################################################################################
clustering.getrisks<- function(risk.table, indir=NA, infile=NA, rskf.name='t.stAgeC', rskf.baseline='D_(30,45]')
{
	#	p.clu	(clustering prob)
	#	model derived from t.stAgeC.tperiod stratification
	tclu		<- subset(risk.table, stat%in%c('R.clu','R.seq')) 
	tclu[, t.period:= substr(factor, nchar(factor),nchar(factor))]
	tclu[, p:=NULL]	
	tclu[, factor2:= substr(factor, 1,nchar(factor)-2)]	
	tmp			<- dcast.data.table(tclu, risk+factor2~stat, value.var='n', fun.aggregate=sum)
	tmp[, t.period:='Overall']
	tclu		<- dcast.data.table(tclu, risk+t.period+factor2~stat, value.var='n')
	tclu		<- rbind(tclu, tmp, fill=T, use.names=T)
	tclu[, p.clu:= R.clu/R.seq]	
	tmp			<- subset(tclu, factor2==rskf.baseline)
	setnames(tmp, c('factor2','p.clu'), c('factor2r','p.clur'))
	tclu		<- merge(tclu, subset(tmp, select=c('t.period','risk','factor2r','p.clur')), by=c('t.period','risk'))
	tclu[, RR:= p.clu/p.clur]
	#	p.pt	(prob of being a prob transmitter if in same cluster)
	#	model derived from t.stAgeC.tperiod stratification
	tpt			<- subset(risk.table, stat%in%c('R.clu','YX')) 
	tpt[, t.period:= substr(factor, nchar(factor),nchar(factor))]
	tpt[, p:=NULL]	
	tpt[, factor2:= substr(factor, 1,nchar(factor)-2)]
	tmp			<- dcast.data.table(tpt, risk+factor2~stat, value.var='n', fun.aggregate=sum)
	tmp[, t.period:='Overall']	
	tpt			<- dcast.data.table(tpt, risk+t.period+factor2~stat, value.var='n')
	tpt			<- rbind(tpt, tmp, fill=T, use.names=T)
	tpt[, p.pt:= YX/R.clu]
	tmp			<- subset(tpt, factor2==rskf.baseline)
	setnames(tmp, c('factor2','p.pt'), c('factor2r','p.ptr'))
	tpt		<- merge(tpt, subset(tmp, select=c('t.period','risk','factor2r','p.ptr')), by=c('t.period','risk'))
	tpt[, RR:= p.pt/p.ptr]
	
	if(!is.na(indir) & !is.na(infile))
	{
		ggplot(subset(tclu,t.period!='Overall'), aes(y=factor2, x=100*p.clu, pch=factor(t.period), size=R.clu)) + geom_point() + 
				scale_x_continuous(breaks=seq(0,10,0.2)) + theme_bw()
		file	<- paste(indir,'/',gsub('\\.R','_ClusteringByTperiod.pdf',infile),sep='')	
		ggsave(file=file, w=6,h=7)
		ggplot(subset(tclu,t.period!='Overall'), aes(y=factor2, x=RR, pch=factor(t.period), size=R.clu)) + geom_point() + 
				geom_vline(xintercept=1, colour='grey50',lwd=1) +
				scale_x_continuous(breaks=seq(0,10,0.5)) + theme_bw()
		file	<- paste(indir,'/',gsub('\\.R','_ClusteringByTperiodRR.pdf',infile),sep='')	
		ggsave(file=file, w=6,h=7)
		ggplot(subset(tclu,t.period=='Overall'), aes(y=factor2, x=100*p.clu, size=R.clu)) + geom_point() + 
				scale_x_continuous(breaks=seq(0,10,0.2)) + theme_bw()
		file	<- paste(indir,'/',gsub('\\.R','_ClusteringOverall.pdf',infile),sep='')	
		ggsave(file=file, w=6,h=7)
		ggplot(subset(tclu,t.period=='Overall'), aes(y=factor2, x=RR, size=R.clu)) + geom_point() + 
				geom_vline(xintercept=1, colour='grey50',lwd=1) +
				scale_x_continuous(breaks=seq(0,10,0.5)) + theme_bw()
		file	<- paste(indir,'/',gsub('\\.R','_ClusteringOverallRR.pdf',infile),sep='')	
		ggsave(file=file, w=6,h=7)
		#		
		ggplot(subset(tpt,t.period!='Overall'), aes(y=factor2, x=100*p.pt, pch=factor(t.period), size=YX)) + geom_point() + 
				scale_x_continuous(breaks=seq(0,100,10)) + theme_bw()
		file	<- paste(indir,'/',gsub('\\.R','_ProbTrByTperiod.pdf',infile),sep='')	
		ggsave(file=file, w=6,h=7)	
		ggplot(subset(tpt,t.period!='Overall'), aes(y=factor2, x=RR, pch=factor(t.period), size=YX)) + geom_point() + 
				geom_vline(xintercept=1, colour='grey50',lwd=1) + scale_x_continuous(breaks=seq(0,10,0.5)) + theme_bw()
		file	<- paste(indir,'/',gsub('\\.R','_ProbTrByTperiodRR.pdf',infile),sep='')	
		ggsave(file=file, w=6,h=7)	
		ggplot(subset(tpt,t.period=='Overall'), aes(y=factor2, x=100*p.pt, size=YX)) + geom_point() + 
				scale_x_continuous(breaks=seq(0,100,10)) + theme_bw()
		file	<- paste(indir,'/',gsub('\\.R','_ProbTrOverall.pdf',infile),sep='')	
		ggsave(file=file, w=6,h=7)	
		ggplot(subset(tpt,t.period=='Overall'), aes(y=factor2, x=RR, size=YX)) + geom_point() + 
				geom_vline(xintercept=1, colour='grey50',lwd=1) + scale_x_continuous(breaks=seq(0,10,0.5)) + theme_bw()
		file	<- paste(indir,'/',gsub('\\.R','_ProbTrOverallRR.pdf',infile),sep='')	
		ggsave(file=file, w=6,h=7)	
	}
	list(tpt=tpt, tclu=tclu)
}
######################################################################################
prop.consistency<- function(YXc, YXm, indir, infile, suffix)	
{
	YXm[, list(miss.pr=sum(miss.pr)), by=c('stageC','t.period')]
	tmp		<- merge( subset(YXm, AnyPos_T1<2010)[, list(miss.pr=sum(miss.pr)), by=c('stageC','t.period')], subset(YXc, AnyPos_T1<2010)[, list(Patient.n= length(unique(Patient))), by='t.period'], by='t.period')
	tmp[, miss.pr:= miss.pr/Patient.n]
	
	#	consistency OK, replicating exactly previous method. this only has different stages
	tmp2	<- merge(subset(YXc, AnyPos_T1<2010)[, list(me=sum(me), mm= length(t)*(1-median(p.seqnc))/median(p.seqnc), n=length(t), pme=median(p.seqnc), scoreYm=median(score.Y)), by=c('stageC','t.period')], subset(YXc, AnyPos_T1<2010)[, list(Patient.n= length(unique(Patient))), by='t.period'], by='t.period')
	#	subset(tmp2, t.period=='4')[, sum(me)]
	tmp2	<- merge(subset(YXc, AnyPos_T1<2010)[, list(me=sum(me), mm= length(t)*(1-mean(p.seqnc))/mean(p.seqnc), n=length(t), pme=median(p.seqnc), scoreYm=median(score.Y)), by=c('stageC','t.period')], subset(YXc, AnyPos_T1<2010)[, list(Patient.n= length(unique(Patient))), by='t.period'], by='t.period')
	#tmp2	<- merge(subset(YXc, AnyPos_T1<2010)[, list(me=sum(me), mm= length(t)*(1-mean(p.seqnc))/mean(p.seqnc), n=length(t), pme=median(p.seqnc), scoreYm=mean(score.Y)), by=c('stageC','t.period')], subset(YXc, AnyPos_T1<2010)[, list(Patient.n= length(unique(Patient))), by='t.period'], by='t.period')	
	tmp2[, mm:= mm/Patient.n]		
	tmp2	<- merge(as.data.table(expand.grid(t.period='4', stageC= subset(YXc, t.period=='4' & AnyPos_T1<2010)[, unique(stageC)], Patient=subset(YXc, t.period=='4' & AnyPos_T1<2010)[, unique(Patient)], stringsAsFactors=0)), subset(tmp2, select=c(t.period, stageC, mm, scoreYm)), by=c('stageC','t.period'))
	tmp3	<- subset(YXc, t.period==4 & AnyPos_T1<2010)[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','stageC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','stageC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	tmp2[, scoreYm.per.recst:= mm*scoreYm]
	setnames(tmp2, c('mm'), c('miss.per.recst'))
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(stageC=stageC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','stageC'))	
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R','_ptp_asbeforeByStageC.pdf',infile), sep=''), w=10, h=7)
	
	#	work now with new p.seqnc instead of median p.seqnc
	tmp2	<- merge(subset(YXc, AnyPos_T1<2010)[, list(me=sum(me), mm= length(t)*(1-mean(p.seqnc))/mean(p.seqnc), n=length(t), pme=median(p.seqnc), scoreYm=median(score.Y)), by=c('stageC','t.period')], subset(YXc, AnyPos_T1<2010)[, list(Patient.n= length(unique(Patient))), by='t.period'], by='t.period')
	tmp2[, me:= me/Patient.n]
	tmp2	<- merge(as.data.table(expand.grid(t.period='4', stageC= subset(YXc, t.period=='4' & AnyPos_T1<2010)[, unique(stageC)], Patient=subset(YXc, t.period=='4' & AnyPos_T1<2010)[, unique(Patient)], stringsAsFactors=0)), subset(tmp2, select=c(t.period, stageC, me, scoreYm)), by=c('stageC','t.period'))
	tmp3	<- subset(YXc, t.period==4 & AnyPos_T1<2010)[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','stageC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','stageC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	tmp2[, scoreYm.per.recst:= me*scoreYm]
	setnames(tmp2, c('me'), c('miss.per.recst'))
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(stageC=stageC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','stageC'))	
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R','_ptp_pnseqncByStageC.pdf',infile), sep=''), w=10, h=7)
	
	#	work now with score.Y instead of median score.Y
	#	not OK
	if(0)
	{
		tmp2	<- merge(subset(YXc, AnyPos_T1<2010)[, list(me=sum(me*score.Y)), by=c('stageC','t.period')], subset(YXc, AnyPos_T1<2010)[, list(Patient.n= length(unique(Patient))), by='t.period'], by='t.period')
		tmp2[, me:= me/Patient.n]
		tmp2	<- merge(as.data.table(expand.grid(t.period='4', stageC= subset(YXc, t.period=='4' & AnyPos_T1<2010)[, unique(stageC)], Patient=subset(YXc, t.period=='4' & AnyPos_T1<2010)[, unique(Patient)], stringsAsFactors=0)), subset(tmp2, select=c(t.period, stageC, me)), by=c('stageC','t.period'))
		tmp3	<- subset(YXc, t.period==4 & AnyPos_T1<2010)[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','stageC')]
		tmp2	<- merge(tmp2, tmp3, by=c('Patient','stageC','t.period'), all.x=1)
		set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
		setnames(tmp2, c('me'), c('scoreYm.per.recst'))
		YXpr	<- tmp2
		YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
		tmp		<- YXpr[, list(stageC=stageC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
		YXpr	<- merge(YXpr, tmp, by=c('Patient','stageC'))	
		tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('stageC','t.period')]
		setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
		tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
		#tmp[,sum(p.per.tc),by='tmidC']	
		ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
				facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')	
		ggsave(file=paste(indir, '/', gsub('\\.R','_ptp_pnseqncscoreYByStageC.pdf',infile), sep=''), w=10, h=7)		
	}
		
	#	work now with smoothed miss.pr instead of raw p.seqnc, keep median score.Y
	tmp2	<- merge(subset(YXc, AnyPos_T1<2010)[, list(me=sum(me), mm= length(t)*(1-mean(p.seqnc))/mean(p.seqnc), n=length(t), pme=median(p.seqnc), scoreYm=median(score.Y)), by=c('stageC','t.period')], subset(YXc, AnyPos_T1<2010)[, list(Patient.n= length(unique(Patient))), by='t.period'], by='t.period')
	tmp2	<- subset(tmp2, select=c(t.period, stageC, scoreYm))
	tmp2	<- merge(subset(YXm, AnyPos_T1<2010), tmp2, by=c('t.period','stageC'))	
	#subset(tmp2, t.period=='4' & AnyPos_T1<2010)[, list(miss.prsm=sum(miss.prsm), miss.pr=sum(miss.pr))]	
	tmp2	<- tmp2[, list(scoreYm.per.recst=sum(miss.prsm*scoreYm)), by=c('t.period','stageC','Patient')]	
	tmp2	<- subset(tmp2, t.period=='4')
	tmp3	<- subset(YXc, t.period==4 & AnyPos_T1<2010)[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','stageC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','stageC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(stageC=stageC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','stageC'))	
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R','_ptp_missprsmByStageC.pdf',infile), sep=''), w=10, h=7)
	
	#	work now with smoothed miss.pr instead of raw p.seqnc, keep median score.Y and force UAE/UAC same score
	tmp2	<- merge(subset(YXc, AnyPos_T1<2010)[, list(me=sum(me), mm= length(t)*(1-mean(p.seqnc))/mean(p.seqnc), n=length(t), pme=median(p.seqnc)), by=c('stageC','t.period')], subset(YXc, AnyPos_T1<2010)[, list(Patient.n= length(unique(Patient))), by='t.period'], by='t.period')
	tmp		<- subset(YXc, AnyPos_T1<2010)
	set(tmp, NULL, 'stageC', tmp[, gsub('UAC|UAE','UA',stageC)])
	tmp		<- tmp[, list(scoreYm=median(score.Y)), by=c('stageC','t.period')]	
	set(tmp, NULL, 'stageC', tmp[, gsub('UA','UAC',stageC)])
	tmp2	<- merge(tmp2, tmp,by=c('stageC','t.period'), all.x=TRUE)
	set(tmp, NULL, 'stageC', tmp[, gsub('UAC','UAE',stageC)])
	tmp2	<- merge(tmp2, tmp,by=c('stageC','t.period'), all.x=TRUE)
	tmp		<- tmp2[, which(is.na(scoreYm.x))]
	set(tmp2, tmp, 'scoreYm.x', tmp2[tmp,scoreYm.y])
	set(tmp2, NULL, 'scoreYm.y', NULL)
	setnames(tmp2, 'scoreYm.x', 'scoreYm')
	tmp2	<- merge(subset(YXm, AnyPos_T1<2010), tmp2, by=c('t.period','stageC'))		
	#subset(tmp2, t.period=='4' & AnyPos_T1<2010)[, list(miss.prsm=sum(miss.prsm), miss.pr=sum(miss.pr))]	
	tmp2	<- tmp2[, list(scoreYm.per.recst=sum(miss.prsm*scoreYm)), by=c('t.period','stageC','Patient')]	
	tmp2	<- subset(tmp2, t.period=='4')
	tmp3	<- subset(YXc, t.period==4 & AnyPos_T1<2010)[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','stageC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','stageC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(stageC=stageC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','stageC'))	
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R','_ptp_missprsmUAEUACSameByStageC.pdf',infile), sep=''), w=10, h=7)
	
	if(0)
	{
		#	work now with smoothed miss.pr instead of raw p.seqnc, and use mean score.Y
		tmp2	<- merge(subset(YXc, AnyPos_T1<2010)[, list(me=sum(me), mm= length(t)*(1-mean(p.seqnc))/mean(p.seqnc), n=length(t), pme=median(p.seqnc), scoreYm=mean(score.Y)), by=c('stageC','t.period')], subset(YXc, AnyPos_T1<2010)[, list(Patient.n= length(unique(Patient))), by='t.period'], by='t.period')
		tmp2	<- subset(tmp2, select=c(t.period, stageC, scoreYm))
		tmp2	<- merge(subset(YXm, AnyPos_T1<2010), tmp2, by=c('t.period','stageC'))	
		#subset(tmp2, t.period=='4' & AnyPos_T1<2010)[, list(miss.prsm=sum(miss.prsm), miss.pr=sum(miss.pr))]	
		tmp2	<- tmp2[, list(scoreYm.per.recst=sum(miss.prsm*scoreYm)), by=c('t.period','stageC','Patient')]	
		tmp2	<- subset(tmp2, t.period=='4')
		tmp3	<- subset(YXc, t.period==4 & AnyPos_T1<2010)[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','stageC')]
		tmp2	<- merge(tmp2, tmp3, by=c('Patient','stageC','t.period'), all.x=1)
		set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
		YXpr	<- tmp2
		YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
		tmp		<- YXpr[, list(stageC=stageC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
		YXpr	<- merge(YXpr, tmp, by=c('Patient','stageC'))	
		tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('stageC','t.period')]
		setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
		tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
		#tmp[,sum(p.per.tc),by='tmidC']	
		ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
				facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
		ggsave(file=paste(indir, '/', gsub('\\.R','_ptp_missprsmscoreYmeanByStageC.pdf',infile), sep=''), w=10, h=7)		
	}
	
	
	#	work now with smoothed miss.pr instead of raw p.seqnc, use predicted median score.Y
	subset(YXc, t.period==4 & AnyPos_T1<2010)[, list(score.Y=median(score.Y)), by='stageC']
	subset(YXm, t.period==4 & AnyPos_T1<2010)[, list(scoreY.pr=median(scoreY.pr)), by='stageC']
	
	tmp2	<- subset(YXm, AnyPos_T1<2010)[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','stageC','Patient')]
	tmp2	<- subset(tmp2, t.period=='4')
	tmp3	<- subset(YXc, t.period==4 & AnyPos_T1<2010)[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','stageC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','stageC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(stageC=stageC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','stageC'))	
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R','_ptp_missprsmscoreYprByStageC.pdf',infile), sep=''), w=10, h=7)
	
	tmp2	<- YXm[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','stageC','Patient')]
	tmp2	<- subset(tmp2, t.period=='4')
	tmp3	<- subset(YXc, t.period==4)[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','stageC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','stageC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(stageC=stageC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','stageC'))	
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	
}
######################################################################################
censoring.check.censdelta<- function()
{
	indir		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C2V100bInfT7STRAT_m5B.R'
	load(paste(indir,'/',infile,sep=''))
	stage.labels<- c("UAC", "UAE", "UC", "DAC", "D", "L", "T")
	c.period	<- 0.125 
	bs.cdelta.min	<- 2
	bs.cdelta.max	<- 3
	
	require(zoo)
	tp.df		<- subset(X.msm, grepl('^U',stageC))[, list(tm=mean(t)), by=c('t.Patient','Patient')]
	tpd.df		<- censoring.get.dataset(tp.df, df.all.allmsm)	
	
	
	cens.bs		<- c(bs.cdelta.min, (bs.cdelta.min+bs.cdelta.max)/2, bs.cdelta.max)
	cens.bs		<- data.table(cens.t=t.recent.endctime-cens.bs, cens.delta=cens.bs, BS=seq_along(cens.bs))
	tpd.df[, NCNS:= tpd.df[, as.numeric(AnyPos_T1<cens.bs[2, cens.t])]]	#transmitter not censored
	set(tpd.df, NULL, 'tm', tpd.df[, tm-cens.bs[2, cens.t]])
	tpd.df[, NCNST:='2008.5']
	tmp			<- copy(tpd.df)
	tmp[, NCNS:= tmp[, as.numeric(AnyPos_T1<cens.bs[3, cens.t])]]	
	tmp[, NCNST:='2008']
	set(tmp, NULL, 'tm', tmp[, cens.bs[2, cens.t]+tm-cens.bs[3, cens.t]])
	tpd.df		<- rbind(tpd.df, tmp)
	tmp[, NCNS:= tmp[, as.numeric(AnyPos_T1<cens.bs[1, cens.t])]]	
	set(tmp, NULL, 'tm', tmp[, cens.bs[3, cens.t]+tm-cens.bs[1, cens.t]])
	tmp[, NCNST:='2009']
	tpd.df		<- rbind(tpd.df, tmp)
	
	tmp			<- rev(seq(t.recent.endctime, 2000, by=-c.period)) - (t.recent.endctime-(bs.cdelta.min+bs.cdelta.max)/2)
	tmp			<- c( tpd.df[, min(tm)-c.period], tmp)
	tpd.df[, tmC:= tpd.df[, cut(tm, breaks=tmp)]]
	tpd.df[, DUMMY:=seq_len(nrow(tpd.df))]
	tmp			<- tpd.df[, {
				list(DUMMY=sample(DUMMY, min( length(DUMMY),c.smpl.n ))) 
			}, by=c('tmC','ACD4C','AgeC_T1','NCNST')]
	tpds.df		<- merge( tpd.df, subset(tmp, select=DUMMY), by='DUMMY' )
	tpds.df[, DUMMY:=NULL]
	
	#	stratify by isAcute + NCNST
	setkey(tpds.df, NCNST, isAcute, tm)
	tmp		<- merge(tpds.df, data.table(isAcute=c('Yes','NI'), DUMMY=c(1e2, 2e3)), by='isAcute')
	setkey(tmp, NCNST, isAcute, tm)		
	z		<- tmp[, list( Patient=Patient, r.Patient=r.Patient, tm=tm, NCNS.rm5=rollapply(NCNS, width=DUMMY[1], FUN=mean, align="center", partial=TRUE) ), by=c('NCNST','isAcute')]
	tpds.df	<- merge(tpds.df, z, by=c('Patient','r.Patient','tm','NCNST','isAcute'))
	
	
	ggplot(tpds.df, aes(x=tm, y=NCNS.rm5, colour=NCNST)) + geom_line() + facet_grid(~isAcute)
	ggsave(file=paste(outdir, '/', outfile, '_CENSMODEL_censdelta_hasnoeffect.pdf', sep=''), w=10, h=5)
}
######################################################################################
censoring.subsample.tmC_ACD4C_AgeC_T1<- function(tpd.df, t.recent.endctime, bs.cdelta.min, bs.cdelta.max, c.period, c.smpl.n)
{
	tpds.df		<- copy(tpd.df)
	tmp			<- rev(seq(t.recent.endctime, 2000, by=-c.period)) - (t.recent.endctime-(bs.cdelta.min+bs.cdelta.max)/2)
	tmp			<- c( tpds.df[, min(tm)-c.period], tmp)
	tpds.df[, tmC:= tpds.df[, cut(tm, breaks=tmp)]]
	tpds.df[, DUMMY:=seq_len(nrow(tpds.df))]
	tmp			<- tpds.df[, {
				list(DUMMY=sample(DUMMY, min( length(DUMMY),c.smpl.n ))) 
			}, by=c('tmC','ACD4C','AgeC_T1')]
	tpds.df		<- merge( tpds.df, subset(tmp, select=DUMMY), by='DUMMY' )
	tpds.df[, DUMMY:=NULL]
	tpds.df
}
######################################################################################
censoring.model.exploreoptions<- function(tpds.df)
{
	#	explore different censoring models
	#	baseline model: depends on midpoint of transmission interval
	setkey(tpds.df, tm)
	tpds.df[, NCNS.rm:=tpds.df[, rollapply(NCNS, width=1e3, FUN=mean, align="center", partial=TRUE)]]	
	ggplot(tpds.df, aes(x=tm, y=NCNS.rm)) + geom_line()
	#	stratify by isAcute
	setkey(tpds.df, isAcute, tm)
	tmp		<- tpds.df[, list( Patient=Patient, r.Patient=r.Patient, tm=tm, NCNS.rm1=rollapply(NCNS, width=1e3, FUN=mean, align="center", partial=TRUE) ), by='isAcute']
	tpds.df	<- merge(tpds.df, tmp, by=c('Patient','r.Patient','tm','isAcute'))
	ggplot(tpds.df, aes(x=tm, y=NCNS.rm1, colour=isAcute)) + geom_line() + scale_x_continuous(breaks=seq(-20,10,5), minor_breaks=seq(-20,10,1))
	ggsave(file=paste(outdir, '/', outfile, '_CENSRMEAN_rm1.pdf', sep=''), w=5, h=5)
	#	substantial differences in censoring
	
	#	stratify by RegionHospital + isAcute
	setkey(tpds.df, RegionHospital, isAcute, tm)
	tmp		<- merge(tpds.df, data.table(isAcute=c('Yes','NI'), DUMMY=c(1e2, 1e3)), by='isAcute')
	setkey(tmp, RegionHospital, isAcute, tm)	
	z		<- tmp[, list( Patient=Patient, r.Patient=r.Patient, tm=tm, NCNS.rm2=rollapply(NCNS, width=DUMMY[1], FUN=mean, align="center", partial=TRUE) ), by=c('RegionHospital','isAcute')]
	tpds.df	<- merge(tpds.df, z, by=c('tm','Patient','r.Patient','RegionHospital','isAcute'))
	setkey(tpds.df, RegionHospital, isAcute, tm)
	ggplot(tpds.df, aes(x=tm, y=NCNS.rm2, colour=RegionHospital)) + geom_line() + facet_grid(~isAcute) + scale_x_continuous(breaks=seq(-20,10,5), minor_breaks=seq(-20,10,1))
	ggsave(file=paste(outdir, '/', outfile, '_CENSRMEAN_rm2.pdf', sep=''), w=10, h=5)
	#	no particular geographical differences
	
	#	stratify by CD4 at diagnosis
	setkey(tpds.df, ACD4C, tm)
	tmp		<- merge(tpds.df, data.table(ACD4C=c("(-Inf,250]","(250,350]","(350,500]","(500,750]","(750, Inf]","CD4.NA","CD4.late","CD4.Acute"), DUMMY=c(rep(1e3,7),1e2)), by='ACD4C')
	setkey(tmp, ACD4C, tm)		
	z		<- tmp[, list( Patient=Patient, r.Patient=r.Patient, tm=tm, NCNS.rm3=rollapply(NCNS, width=DUMMY[1], FUN=mean, align="center", partial=TRUE) ), by=c('ACD4C')]
	tpds.df	<- merge(tpds.df, z, by=c('tm','Patient','r.Patient','ACD4C'))
	setkey(tpds.df, ACD4C, tm)
	ggplot(tpds.df, aes(x=tm, y=NCNS.rm3, colour=ACD4C)) + geom_line() + scale_x_continuous(breaks=seq(-20,10,5), minor_breaks=seq(-20,10,1))
	ggsave(file=paste(outdir, '/', outfile, '_CENSRMEAN_rm3.pdf', sep=''), w=5, h=5)
	#	for isAcute=Yes, no stratification needed
	#	for isAcute=NI, <250, 250-350, 350-500, 500-750, >750, CD4.NA, CD4.late seem good
	
	#	stratify by age at diagnosis + CD4 at diagnosis
	setkey(tpds.df, AACD4C, tm)
	tmp		<- data.table(AACD4C=tpds.df[, levels(AACD4C)], DUMMY=2e2)
	set(tmp, tmp[, which(grepl('CD4.Acute',AACD4C,fixed=T))], 'DUMMY', 1e2)
	tmp		<- merge(tpds.df, tmp, by='AACD4C')
	setkey(tmp, AACD4C, tm)		
	z		<- tmp[, list( Patient=Patient, r.Patient=r.Patient, tm=tm, NCNS.rm4=rollapply(NCNS, width=DUMMY[1], FUN=mean, align="center", partial=TRUE) ), by=c('AACD4C')]
	tpds.df	<- merge(tpds.df, z, by=c('tm','Patient','r.Patient','AACD4C'))
	setkey(tpds.df, AACD4C, tm)
	ggplot(tpds.df, aes(x=tm, y=NCNS.rm4, colour=AgeC_T1)) + geom_line() + facet_grid(~ACD4C) + scale_x_continuous(breaks=seq(-20,10,5), minor_breaks=seq(-20,10,1))
	ggsave(file=paste(outdir, '/', outfile, '_CENSRMEAN_rm4.pdf', sep=''), w=20, h=5)
	#	for isAcute=NI, Age is an independent effect in addition to CD4C
	
	#
	#	build models
	#
	tprs.df	<- subset(tpds.df, select=c(tm, Patient, r.Patient, AgeC_T1, CD4C, ACD4C, AACD4C, isAcute, RegionHospital, NCNS, NCNS.rm, NCNS.rm1, NCNS.rm2, NCNS.rm3, NCNS.rm4))	
	#	isAcute
	cs1		<- gamlss(formula= NCNS ~ isAcute*ns(tm, df=6), family=BI(), data=tprs.df)
	tprs.df[, CS1:=predict(cs1, type='response')]	
	ggplot(melt(tprs.df, measure.vars=c('NCNS.rm1','CS1'), id.vars=c('tm','isAcute')), aes(x=tm, y=value, group=variable, colour=variable)) +
			scale_x_continuous(breaks=seq(1980,2020,5), minor_breaks=seq(1980,2020,1)) +
			geom_line() + facet_grid(~isAcute)
	ggsave(file=paste(outdir, '/', outfile, '_CENSMODEL_cs1.pdf', sep=''), w=10, h=5)
	#
	#	isAcute and CD4
	#cs3		<- gamlss(formula= NCNS ~ ACD4C*ns(tm, df=6), family=BI(), data=tprs.df)
	#	df=6 results in bumps; df3 looks much better
	cs3		<- gamlss(formula= NCNS ~ ACD4C*ns(tm, df=3), family=BI(), data=tprs.df)	
	tprs.df[, CS3:=predict(cs3, type='response')]
	tmp		<- melt(tprs.df, measure.vars=c('NCNS.rm3','CS3'), id.vars=c('tm','ACD4C'))
	ggplot(tmp, aes(x=tm, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(~ACD4C)
	ggsave(file=paste(outdir, '/', outfile, '_CENSMODEL_cs3.pdf', sep=''), w=20, h=5)
	#
	#	isAcute, CD4, Age at diagnosis
	cs4		<- gamlss(formula= NCNS ~ AACD4C*cs(tm, df=5), family=BI(), data=tprs.df)
	tprs.df[, CS4:=predict(cs4, type='response')]
	tmp		<- melt(tprs.df, measure.vars=c('NCNS.rm4','CS4'), id.vars=c('tm','ACD4C','AgeC_T1'))
	ggplot(tmp, aes(x=tm, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(AgeC_T1~ACD4C)
	ggsave(file=paste(outdir, '/', outfile, '_CENSMODEL_cs4.pdf', sep=''), w=20, h=15)
	#z		<- subset(tprs.df, AgeC_T1=='(-Inf,25]')
	#cs4b	<- gamlss(formula= NCNS ~ ACD4C*cs(tm, df=3), family=BI(), data=z)
	#z[, CS4:=predict(cs4b, type='response')]
	#ggplot(melt(z, measure.vars=c('NCNS.rm4','CS4'), id.vars=c('tm','ACD4C','AgeC_T1')), aes(x=tm, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(AgeC_T1~ACD4C)	
}
######################################################################################
censoring.model.150728<- function(tpds.df)
{
	tprs.df	<- subset(tpds.df, select=c(tm, Patient, r.Patient, AACD4C, NCNS))
	cs4		<- gamlss(formula= NCNS ~ AACD4C*cs(tm, df=5), family=BI(), data=tprs.df)
	tprs.df[, p.nc:=predict(cs4, type='response')]
	list(predict=tprs.df, model=cs4)
	#tmp		<- melt(tprs.df, measure.vars=c('NCNS.rm4','CS4'), id.vars=c('tm','ACD4C','AgeC_T1'))
	#ggplot(tmp, aes(x=tm, y=value, group=variable, colour=variable)) + geom_line() + facet_grid(AgeC_T1~ACD4C)
	#ggsave(file=paste(outdir, '/', outfile, '_CENSMODEL_cs4.pdf', sep=''), w=20, h=15)	
}
######################################################################################
censoring.model.calculate.bs.args<- function(method)
{	
	if(grepl('m5A',method))
	{
		risk.col		<- 'stageC'			
	}
	if(grepl('m5B',method))
	{
		risk.col		<- 'stageC'			
	}
	if(grepl('m5C',method))
	{
		risk.col		<- 'stageC'			
	}
	if(grepl('m5D',method))
	{
		risk.col		<- 'stageC'			
	}
	if(grepl('m5E',method))
	{
		risk.col		<- 'stageC'			
	}
	if(grepl('m5F',method))
	{
		risk.col		<- 'stageC'			
	}
	risk.col
}
######################################################################################
censoring.model.calculate.bs<- function(X.msm, df.all.allmsm, resume=TRUE, save.file=NA, t.recent.endctime=NA, risk.col=NA, c.period=0.125, c.smpl.n=50, bs.n=100, bs.cdelta.min=2, bs.cdelta.max=3)
{
	require(zoo)	
	if(resume & !is.na(save.file))
	{
		options(show.error.messages = FALSE)
		tmp				<- gsub('\\.R',paste('_bs',bs.n,'\\.R',sep=''),save.file)
		readAttempt		<- try(suppressWarnings(load(tmp)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",tmp))					
		options(show.error.messages = TRUE)		
	}
	if(!resume || is.na(save.file) || inherits(readAttempt, "try-error"))
	{		
		if(is.null(X.msm))
			return(NULL)
		stopifnot(!is.na(risk.col),  !is.na(t.recent.endctime))
		cens.bs		<- runif(bs.n, bs.cdelta.min, bs.cdelta.max)
		cens.bs		<- data.table(cens.t=t.recent.endctime-cens.bs, cens.delta=cens.bs, BS=seq_along(cens.bs))
		cens.m		<- vector('list', length=bs.n)
		cens.p		<- vector('list', length=bs.n)
		tp.df		<- subset(X.msm, grepl('^U', X.msm[[risk.col]]))[, list(tm=mean(t)), by=c('t.Patient','Patient')]
		tpd.df		<- censoring.get.dataset(tp.df, df.all.allmsm)
		tpd.df[, tmo:= tpd.df[, tm]]
		#	add pseudo censoring time
		for(bs in seq_len(bs.n))
		{
			tpd.df[, NCNS:= tpd.df[, as.numeric(AnyPos_T1<cens.bs[bs, cens.t])]]	#transmitter not censored
			set(tpd.df, NULL, 'tm', tpd.df[, tmo-cens.bs[bs, cens.t]])
			#	balanced sub sampling so regression is comp feasible
			tpds.df		<- censoring.subsample.tmC_ACD4C_AgeC_T1(tpd.df, t.recent.endctime, bs.cdelta.min, bs.cdelta.max, c.period, c.smpl.n)
			tmp			<- censoring.model.150728(tpds.df)
			cens.m		<- tmp$model
			cens.p		<- copy(tmp$predict)
			cens.p[, BS:=bs]
			#	save
			if(!is.na(save.file))
			{				
				save(tp.df, tpd.df, cens.p, cens.m, file=gsub('\\.R',paste('_bs',bs,'\\.R',sep=''),save.file))
			}							
		}		
	}
	list(cens.p=cens.p, cens.m=cens.m)
}
######################################################################################
adjust.dev.code.for.ntPatient.adjustment<- function()
{
	indir		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C2V100bInfT7STRAT_m5A.R'
	
	YX
	sm		<- X.tables$sm
	YXs		<- merge(YX, sm, by='t.Patient')
	
	#	get w per interval
	cat(paste('\nsetting likelihood to likelihood of pair / number transmission intervals'))
	set(YXs, NULL, 'score.Y.raw', YXs[, score.Y])
	set(YXs, NULL, 'score.Y', YXs[, score.Y*w.tn])				
	YXs[, ntPatient:= factor(YXs[, 1/w.in])]
	
	#	are the w's confounded by cluster size?
	tmp		<- YXs[, list(score.Y.m=mean(score.Y), stageC=stageC[1], t.AgeC=t.AgeC[1]), by=c('ntPatient','t.stAgeC')]
	ggplot( tmp, aes(x=ntPatient, y=score.Y.m)) + geom_bar(stat='identity') + facet_grid(stageC~t.AgeC)
	file	<- paste(indir,'/',gsub('\\.R','_scoreYmeansByStageNtPatient_raw.pdf',infile),sep='')	
	ggsave(file=file, w=20,h=10)
	tmp		<- YXs[, list(score.Y.mx=max(score.Y)), by=c('ntPatient','Patient')]
	tmp		<- tmp[, list(score.Y.mx.me=mean(score.Y.mx)), by='ntPatient']
	ggplot( tmp, aes(x=ntPatient, y=score.Y.mx.me)) + geom_bar(stat='identity')
	tmp		<- YXs[, list(score.Y.me=mean(score.Y),score.Y.md=median(score.Y),n=length(score.Y)), by=c('ntPatient','Patient')][, list(n=as.double(sum(n)),score.Y.me=mean(score.Y.me), score.Y.md=mean(score.Y.md)), by='ntPatient']
	ggplot( melt(tmp, id.vars='ntPatient'), aes(x=ntPatient, y=value)) + geom_bar(stat='identity') + facet_wrap(~variable, scales='free_y')
	file	<- paste(indir,'/',gsub('\\.R','_scoreYmeansByNtPatient_raw.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=7)
	
	ggplot(YXs, aes(x=score.Y.raw, y=score.Y, colour=as.numeric(as.character(ntPatient)))) + geom_point(alpha=0.5) + geom_abline(slope=1, intercept=0) + theme_bw()
	#
	#	no, of course not!! the p's will be confounded by cluster size!!
	#
	
	
	#	need mean or median w for each stage: problem is it s confounded by cluster size..	
	YXr		<- subset(YXs, select=c('t.Patient','t','Patient', 'score.Y','ntPatient','stageC','t.AgeC', rfactor)) 	
	#	Exp has default link log
	mExp	<- gamlss( score.Y~ntPatient, data=YXr, family=EXP() )
	#	Gamma has default links log log
	mGA		<- gamlss( score.Y~ntPatient, data=YXr, family=GA() )
	mGAS	<- gamlss( score.Y~ntPatient, sigma.formula=~ntPatient, data=YXr, family=GA() )
	
	setkey(YXr, ntPatient)
	mpars	<- data.table(ntPatient= YXr[, levels(ntPatient)])
	mpars	<- merge(unique(YXr), mpars, by='ntPatient')
	mpars[, score.mu:= predict(mExp, data=YXr, newdata=mpars, type='response', what='mu')]
	mpars[, TYPE:= 'Exp']
	tmp		<- data.table(ntPatient= YXr[, levels(ntPatient)])
	tmp		<- merge(unique(YXr), tmp, by='ntPatient')
	mu		<- predict(mGA, data=YXr, newdata=tmp, type='response', what='mu')
	sigma	<- predict(mGA, data=YXr, newdata=tmp, type='response', what='sigma')
	tmp[, score.mu:= mu]
	tmp[, score.sigma:= sigma]
	tmp[, TYPE:= 'GA']
	mpars	<- rbind(mpars, tmp, use.names=TRUE, fill=TRUE)
	tmp		<- data.table(ntPatient= YXr[, levels(ntPatient)])
	tmp		<- merge(unique(YXr), tmp, by='ntPatient')
	mu		<- predict(mGAS, data=YXr, newdata=tmp, type='response', what='mu')
	sigma	<- predict(mGAS, data=YXr, newdata=tmp, type='response', what='sigma')
	tmp[, score.mu:= mu]
	tmp[, score.sigma:= sigma]
	tmp[, TYPE:= 'GAS']
	mpars	<- rbind(mpars, tmp, use.names=TRUE, fill=TRUE)
}
######################################################################################
altvtp.explore.time.rollingmean<- function(YXc, indir, infile)
{
	require(zoo)	
	YXr		<- subset(YXc, select=c('t.Patient','t','Patient', 'score.p','score.p.nadj','score.Y','ntPatient','stageC','t.AgeC','t.stAgeC'))
	
	#
	#	get smooth in terms of t.AgeC	6m rolling mean
	#
	YXp		<- copy(YXr)
	setnames(YXp, 't.AgeC','RSKF')
	tmp		<- YXp[, list(nt=length(score.Y)), by=c('RSKF','t')]
	setkey(tmp, RSKF, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=0.5) ])), cnt=cumsum(nt)	), by=c('RSKF')]
	YXp		<- merge(YXp, tmp, by=c('RSKF','t'))		
	YXp		<- melt( YXp, measure.vars=c('score.p','score.p.nadj','missexp.pr','score.Y'), variable.name='STAT', value.name='V' )
	#tmp		<- YXp[, list(V=mean(V), nt=nt[1], cnt=cnt[1], wnt=wnt[1]), by=c('t','RSKF','STAT')]
	#tmp		<- tmp[, list(t=t, V.rm=rollapply(V, width=ceiling(max(cnt)/100), FUN=mean, align="center", partial=TRUE)), by=c('STAT','RSKF')]
	setkey(YXp, STAT, RSKF, t)
	tmp		<- YXp[, {
				z	<- unique(t)
				list(t=z, V.rm= sapply(seq_along(z), function(i)	mean(V[ which(abs(z[i]-t)<=0.5) ]))	)	
			}, by=c('STAT','RSKF')]
	YXp		<- merge(YXp, tmp, by=c('STAT','t','RSKF'))		
	tmp		<- melt(subset(YXp, STAT=='score.p'), id.vars=c('t.Patient','Patient','t','RSKF','stageC','t.stAgeC','ntPatient'), measure.vars=c('wnt'), variable.name='STAT', value.name='V')
	set(YXp, NULL, c('nt','wnt','cnt'), NULL)
	tmp[, V.rm:=V]	
	YXp		<- rbind(YXp, tmp, use.names=TRUE)
	tmp[, ALPHA:=0.2]
	set(tmp, tmp[, which(V>100)],'ALPHA',1)
	YXp		<- merge(YXp, subset(tmp, select=c('t.Patient','Patient','t','ALPHA')), by=c('t.Patient','Patient','t'))
	ggplot(subset(YXp, !(STAT=='score.p.nadj' & V.rm>0.3)), aes(x=t, y=V.rm, colour=RSKF, group=RSKF, alpha=ALPHA)) + geom_line() + geom_point(size=1.2) +
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(STAT~., scales='free') + labs(y='absolute phylogenetic transmission probability\n6m rolling mean\n(%)') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R','_scorePmeantAgeCByTime_6mrollingmean.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=10)	
	#
	#	get smooth in terms of t.AgeC	rolling mean by intervals of length cnt/100 
	#
	YXp		<- copy(YXr)
	setnames(YXp, 't.AgeC','RSKF')
	tmp		<- YXp[, list(nt=length(score.Y)), by=c('RSKF','t')]
	setkey(tmp, RSKF, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=0.5) ])), cnt=cumsum(nt)	), by=c('RSKF')]
	YXp		<- merge(YXp, tmp, by=c('RSKF','t'))		
	YXp		<- melt( YXp, measure.vars=c('score.p','score.p.nadj','missexp.pr','score.Y'), variable.name='STAT', value.name='V' )
	tmp		<- YXp[, list(V=mean(V), nt=nt[1], cnt=cnt[1], wnt=wnt[1]), by=c('t','RSKF','STAT')]
	setkey(YXp, STAT, RSKF, t)
	tmp		<- tmp[, list(t=t, V.rm=rollapply(V, width=ceiling(max(cnt)/100), FUN=mean, align="center", partial=TRUE)), by=c('STAT','RSKF')]	
	YXp		<- merge(YXp, tmp, by=c('STAT','t','RSKF'))		
	tmp		<- melt(subset(YXp, STAT=='score.p'), id.vars=c('t.Patient','Patient','t','RSKF','stageC','t.stAgeC','ntPatient'), measure.vars=c('wnt'), variable.name='STAT', value.name='V')
	set(YXp, NULL, c('nt','wnt','cnt'), NULL)
	tmp[, V.rm:=V]	
	YXp		<- rbind(YXp, tmp, use.names=TRUE)
	tmp[, ALPHA:=0.2]
	set(tmp, tmp[, which(V>100)],'ALPHA',1)
	YXp		<- merge(YXp, subset(tmp, select=c('t.Patient','Patient','t','ALPHA')), by=c('t.Patient','Patient','t'))
	ggplot(subset(YXp, !(STAT=='score.p.nadj' & V.rm>0.3)), aes(x=t, y=V.rm, colour=RSKF, group=RSKF, alpha=ALPHA)) + geom_line() + geom_point(size=1.2) +
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(STAT~., scales='free') + labs(y='absolute phylogenetic transmission probability\ncnt/100 rolling mean\n(%)') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R','_scorePmeantAgeCByTime_cnt100rollingmean.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=10)	
	#
	#	get smooth in terms of t.stAgeC	1y rolling mean
	#
	YXp		<- copy(YXr)
	setnames(YXp, 't.stAgeC','RSKF')
	tmp		<- YXp[, list(nt=length(score.Y)), by=c('RSKF','t')]
	setkey(tmp, RSKF, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ])), cnt=cumsum(nt)	), by=c('RSKF')]
	YXp		<- merge(YXp, tmp, by=c('RSKF','t'))		
	YXp		<- melt( YXp, measure.vars=c('score.p','score.p.nadj','score.Y'), variable.name='STAT', value.name='V' )
	#tmp		<- YXp[, list(V=mean(V), nt=nt[1], cnt=cnt[1], wnt=wnt[1]), by=c('t','RSKF','STAT')]
	#tmp		<- tmp[, list(t=t, V.rm=rollapply(V, width=ceiling(max(cnt)/100), FUN=mean, align="center", partial=TRUE)), by=c('STAT','RSKF')]
	setkey(YXp, STAT, RSKF, t)
	tmp		<- YXp[, {
				z	<- unique(t)
				list(t=z, V.rm= sapply(seq_along(z), function(i)	mean(V[ which(abs(z[i]-t)<=1) ]))	)	
			}, by=c('STAT','RSKF')]
	YXp		<- merge(YXp, tmp, by=c('STAT','t','RSKF'))		
	tmp		<- melt(subset(YXp, STAT=='score.p'), id.vars=c('t.Patient','Patient','t','RSKF','stageC','t.AgeC','ntPatient'), measure.vars=c('wnt'), variable.name='STAT', value.name='V')
	set(YXp, NULL, c('nt','wnt','cnt'), NULL)
	tmp[, V.rm:=V]	
	YXp		<- rbind(YXp, tmp, use.names=TRUE)
	tmp[, ALPHA:=0.2]
	set(tmp, tmp[, which(V>100)],'ALPHA',1)
	YXp		<- merge(YXp, subset(tmp, select=c('t.Patient','Patient','t','ALPHA')), by=c('t.Patient','Patient','t'))
	ggplot(subset(YXp, !(STAT=='score.p.nadj' & V.rm>0.3)), aes(x=t, y=V.rm, colour=t.AgeC, group=RSKF, alpha=ALPHA)) + geom_line() + geom_point(size=1.2) +
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(STAT~stageC, scales='free') + labs(y='absolute phylogenetic transmission probability\n1yr rolling mean\n(%)') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R','_scorePmeanstAgeCByTime_1yrollingmean.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=10)	
	#
	#	get smooth in terms of t.stAgeC	by intervals of length cnt/100
	#
	YXp		<- copy(YXr)
	setnames(YXp, 't.stAgeC','RSKF')
	tmp		<- YXp[, list(nt=length(score.Y)), by=c('RSKF','t')]
	setkey(tmp, RSKF, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=0.5) ])), cnt=cumsum(nt)	), by=c('RSKF')]
	YXp		<- merge(YXp, tmp, by=c('RSKF','t'))		
	YXp		<- melt( YXp, measure.vars=c('score.p','score.p.nadj','missexp.pr','score.Y'), variable.name='STAT', value.name='V' )
	setkey(YXp, STAT, RSKF, t)
	tmp		<- YXp[, list(V=mean(V), nt=nt[1], cnt=cnt[1], wnt=wnt[1]), by=c('t','RSKF','STAT')]
	setkey(YXp, STAT, RSKF, t)
	tmp		<- tmp[, list(t=t, V.rm=rollapply(V, width=ceiling(max(cnt)/100), FUN=mean, align="center", partial=TRUE)), by=c('STAT','RSKF')]
	YXp		<- merge(YXp, tmp, by=c('STAT','t','RSKF'))		
	tmp		<- melt(subset(YXp, STAT=='score.p'), id.vars=c('t.Patient','Patient','t','RSKF','stageC','t.AgeC','ntPatient'), measure.vars=c('wnt'), variable.name='STAT', value.name='V')
	set(YXp, NULL, c('nt','wnt','cnt'), NULL)
	tmp[, V.rm:=V]	
	YXp		<- rbind(YXp, tmp, use.names=TRUE)
	tmp[, ALPHA:=0.2]
	set(tmp, tmp[, which(V>100)],'ALPHA',1)
	YXp		<- merge(YXp, subset(tmp, select=c('t.Patient','Patient','t','ALPHA')), by=c('t.Patient','Patient','t'))
	ggplot(subset(YXp, !(STAT=='score.p.nadj' & V.rm>0.3)), aes(x=t, y=V.rm, colour=t.AgeC, group=RSKF, alpha=ALPHA)) + geom_line() + geom_point(size=1.2) +
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(STAT~stageC, scales='free') + labs(y='absolute phylogenetic transmission probability\ncnt/100 rolling mean\n(%)') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R','_scorePmeanstAgeCByTime_cnt100rollingmean.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=10)	
}
######################################################################################
altvtp.explore.time.models.150904<- function(YXc, indir, infile, suffix='')
{
	YXr		<- subset(YXc, select=c('t', 'score.p','ntPatient','stageC','t.AgeC','t.stAgeC'))
	YXr		<- subset(YXr, !stageC%in%c('L','T'))
	tmp		<- YXr[, list(nt=length(score.p)), by=c('t.stAgeC','t')]
	setkey(tmp, t.stAgeC, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ])), cnt=cumsum(nt)	), by=c('t.stAgeC')]
	YXr		<- merge(YXr, tmp, by=c('t.stAgeC','t'))
	setkey(YXr, t.stAgeC, t)
	tmp		<- YXr[, {
				z	<- unique(t)
				list(t=z, score.p.rm= sapply(seq_along(z), function(i)	mean(score.p[ which(abs(z[i]-t)<=1) ]))	)	
			}, by='t.stAgeC']
	YXr		<- merge(YXr, tmp, by=c('t','t.stAgeC'))
	#	
	am1		<- gamlss( score.p~t.stAgeC+t:t.stAgeC-1, sigma.formula=~t.stAgeC, data=YXr, family=GA() )
	am2		<- gamlss( score.p~t.stAgeC+t:t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient, data=YXr, family=GA() )
	#am3		<- gamlss( score.p~t.stAgeC+t:t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient, data=subset(YXr, score.p>0.001), family=GA(mu.link='identity'), i.control = glim.control(bf.trace=1, glm.trace=1), method=CG() )	
	#am3		<- gamlss( log(score.p)~t.stAgeC+t:t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient, data=YXr, family=NO() )	
	#	get GA(mu.link='identity') model for score.p.rm
	#	need to deal with convergence issues
	YXh				<- subset(YXr, !stageC%in%c('L','T') & !(t.stAgeC=='D_(45,100]' & t>2010))
	set(YXh, NULL, 'score.p', NULL)
	setnames(YXh,'score.p.rm','score.p')	
	gamlss.limit	<- c(0, 0.001, 0.005, 0.006, 0.0061, 0.00623, 0.00623, 0.0063, 0.00638, 0.00639, 0.0064, 0.0065, 0.0066, 0.0067, 0.0068, 0.0069, 0.007, 0.008, 0.009, 0.01 )
	mu.start		<- NULL
	sigma.start		<- NULL
	for(i in rev(seq_along(gamlss.limit)[-1]))
	{
		cat('\nat iteration', gamlss.limit[i],'\n')
		am3h		<- gamlss( 	score.p~t.stAgeC+t:t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient, 
								data=subset(YXh, score.p>gamlss.limit[i]), family=GA(mu.link='identity'), i.control = glim.control(bf.trace=1, glm.trace=1), mu.start=mu.start, sigma.start=sigma.start)		
		tmp			<- subset(YXh, score.p>gamlss.limit[i-1])
		mu.start	<- predict(am3h, data=subset(YXh, score.p>gamlss.limit[i]), newdata=tmp, type='response', what='mu')
		sigma.start	<- predict(am3h, data=subset(YXh, score.p>gamlss.limit[i]), newdata=tmp, type='response', what='sigma')		
	}	
	YXp				<- copy(YXr)
	set(YXp, NULL, 'score.p', NULL)
	setnames(YXp,'score.p.rm','score.p')		
	mu.start		<- predict(am3h, data=YXh, newdata=YXp, type='response', what='mu')
	sigma.start		<- predict(am3h, data=YXh, newdata=YXp, type='response', what='sigma')	
	am3rm			<- gamlss( 	score.p~t.stAgeC+t:t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient, 
								data=YXp, family=GA(mu.link='identity'), i.control = glim.control(bf.trace=1, glm.trace=1), mu.start=mu.start, sigma.start=sigma.start)	
	#
	YXh				<- subset(YXr, score.p>0.03 & !stageC%in%c('L','T') & !(t.stAgeC=='D_(45,100]' & t>2010), select=c(t.stAgeC, t, ntPatient, stageC, t.AgeC, nt, wnt, cnt, score.p))					
	mu.start		<- predict(am3rm, data=YXp, newdata=YXh, type='response', what='mu')
	sigma.start		<- predict(am3rm, data=YXp, newdata=YXh, type='response', what='sigma')	
	gamlss( 	score.p~t.stAgeC+t:t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient, 
								data=subset(YXh, score.p>0.03), family=GA(mu.link='identity'), i.control = glim.control(bf.trace=1, glm.trace=1), mu.start=mu.start, sigma.start=sigma.start)					
	YXpp	<- as.data.table(expand.grid(t=YXr[, sort(unique(t))], t.stAgeC=YXr[, unique(as.character(t.stAgeC))], ntPatient=YXr[, median(ntPatient)]))	
	yam1	<- predict(am1, data=YXr, newdata=YXpp, type='response', what='mu')
	yam2	<- predict(am2, data=YXr, newdata=YXpp, type='response', what='mu')
	yam3	<- predict(am3rm, data=YXp, newdata=YXpp, type='response', what='mu')	
	YXpp	<- merge(YXpp, unique(subset(YXr, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	YXpp[, score.p.am1:=yam1]
	YXpp[, score.p.am2:=yam2]
	YXpp[, score.p.am3:=yam3]
	YXpp	<- melt(YXpp, measure.vars=c('score.p.am1', 'score.p.am2', 'score.p.am3'))
	ggplot(YXpp, aes(x=t, colour=as.character(t.AgeC))) + 
			geom_point(data=YXr, aes(y=score.p.rm), size=1.2) +
			geom_line(aes(y=value, group=as.character(t.stAgeC))) + 
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(variable~stageC, scales='free') + labs(y='absolute phylogenetic transmission probability\n1yr rolling mean\n(%)') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))	
	file	<- paste(indir,'/',gsub('\\.R',paste('_scorePmeanstAgeCByTimeModels_1yrollingmean',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=10,h=10)		
}
######################################################################################
altvtp.explore.time.models.150909<- function(YXc, indir, infile, suffix='')
{
	YXr		<- subset(YXc, select=c('t', 'score.p','ntPatient','stageC','t.AgeC','t.stAgeC'))
	YXr		<- subset(YXr, !stageC%in%c('L','T'))
	tmp		<- YXr[, list(nt=length(score.p)), by=c('t.stAgeC','t')]
	setkey(tmp, t.stAgeC, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ])), cnt=cumsum(nt)	), by=c('t.stAgeC')]
	YXr		<- merge(YXr, tmp, by=c('t.stAgeC','t'))
	setkey(YXr, t.stAgeC, t)
	tmp		<- YXr[, {
				z	<- unique(t)
				list(t=z, score.p.rm= sapply(seq_along(z), function(i)	mean(score.p[ which(abs(z[i]-t)<=1) ]))	)	
			}, by='t.stAgeC']
	YXr		<- merge(YXr, tmp, by=c('t','t.stAgeC'))
	#	
	am1		<- gamlss( score.p~t.stAgeC+t:t.stAgeC-1, sigma.formula=~t.stAgeC, data=YXr, family=GA() )
	am2		<- gamlss( score.p~t.stAgeC+t:t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient, data=YXr, family=GA() )
	am3		<- gamlss( score.p~t.stAgeC+ns(t,df=2):t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient, data=YXr, family=GA() )
	am4		<- gamlss( score.p~t.stAgeC+ns(t,df=3):t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient, data=YXr, family=GA() )
	am5		<- gamlss( score.p~t.stAgeC+ns(t,df=4):t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient, data=YXr, family=GA() )
	
	YXpp	<- as.data.table(expand.grid(t=YXr[, sort(unique(t))], t.stAgeC=YXr[, unique(as.character(t.stAgeC))], ntPatient=YXr[, median(ntPatient)]))	
	yam1	<- predict(am1, data=YXr, newdata=YXpp, type='response', what='mu')
	yam2	<- predict(am2, data=YXr, newdata=YXpp, type='response', what='mu')
	yam3	<- predict(am3, data=YXp, newdata=YXpp, type='response', what='mu')
	yam4	<- predict(am4, data=YXp, newdata=YXpp, type='response', what='mu')
	yam5	<- predict(am5, data=YXp, newdata=YXpp, type='response', what='mu')
	YXpp	<- merge(YXpp, unique(subset(YXr, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	YXpp[, score.p.am1:=yam1]
	YXpp[, score.p.am2:=yam2]
	YXpp[, score.p.am3:=yam3]
	YXpp[, score.p.am4:=yam4]
	YXpp[, score.p.am5:=yam5]
	YXpp	<- melt(YXpp, measure.vars=c('score.p.am1', 'score.p.am2', 'score.p.am3', 'score.p.am4', 'score.p.am5'))
	ggplot(YXpp, aes(x=t, colour=as.character(t.AgeC))) + 
			geom_point(data=YXr, aes(y=score.p.rm), size=1.2) +
			geom_line(aes(y=value, group=as.character(t.stAgeC))) + 
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(variable~stageC, scales='free') + labs(y='absolute phylogenetic transmission probability\n1yr rolling mean\n(%)') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))	
	file	<- paste(indir,'/',gsub('\\.R',paste('_scorePmeanstAgeCByTimeModels_1yrollingmean',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=10,h=10)	
	
	YXr		<- subset(YXc, select=c('t', 'score.p','ntPatient','stageC','t.AgeC','t.stAgeC'))
	ggplot(YXr, aes(x=t, y=100*score.p, colour=as.character(t.AgeC))) + 
			geom_point(size=1.2, alpha=0.5) +
			geom_smooth(colour='black', method='loess', span=0.5, degree=1) + 
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			scale_y_continuous(expand=c(0,0)) +
			facet_grid(t.AgeC~stageC, scales='free') + labs(y='absolute phylogenetic transmission probability\n(%)') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4), legend.position='bottom')
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreP_tAgeCByTime_loess',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=14,h=14)	
}
######################################################################################
altvp.gamlss.start<- function(YXr,mu.formula,sigma.formula,gamlss.limit,gamlss.family,gamlss.control,gamlss.i.control,verbose=FALSE)
{
	mu.start		<- NULL
	sigma.start		<- NULL
	for(i in rev(seq_along(gamlss.limit)[-1]))
	{
		if(verbose)
			cat('\nat iteration', gamlss.limit[i],'\n')
		tmp.mo		<- gamlss( 	as.formula(mu.formula), sigma.formula=as.formula(sigma.formula), 
				data=subset(YXr, score.p>gamlss.limit[i]), family=gamlss.family, control=gamlss.control, i.control=gamlss.i.control, mu.start=mu.start, sigma.start=sigma.start)		
		tmp			<- subset(YXr, score.p>gamlss.limit[i-1])
		mu.start	<- predict(tmp.mo, data=subset(YXr, score.p>gamlss.limit[i]), newdata=tmp, type='response', what='mu')
		sigma.start	<- predict(tmp.mo, data=subset(YXr, score.p>gamlss.limit[i]), newdata=tmp, type='response', what='sigma')		
	}	
	list(mu=mu.start, sigma=sigma.start)		
}
######################################################################################
altvtp.readcoef<- function(mo, tmp2, tmp3, tmp4)
{
	tmp		<- coef(mo)
	#tmp		<- tmp[ grepl(tmp2,names(tmp)) ]
	cf		<- data.table(MO_ID=tmp3, MO=tmp4, ST=tmp2, RSKF=names(tmp), RSKFbaseline=base.df[tmp2,][,RSKFbaseline], MU=tmp)
	tmp		<- confint(mo)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MUL= tmp[,1], MUU= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub(tmp2,'',RSKF)])
	cf
}
######################################################################################
altvtp.model.150911<- function(YXc)
{	
	YXr		<- subset(YXc, score.p>1e-10, select=c('t.Patient','t','Patient', 'score.p','ntPatient','stageC','t.AgeC','p.AgeC','AgeC','t.stAgeC'))
	YXr2	<- subset(YXc, !stageC%in%c('T','L'),select=c('t.Patient','t','Patient', 'score.p','ntPatient','stageC','t.AgeC','p.AgeC','AgeC','t.stAgeC'))
	
	wbm3	<- gamlss( score.p~t.stAgeC-1, sigma.formula=~t.stAgeC-1, data=YXr, family=GA() )	
	wbm2	<- gamlss( score.p~t.stAgeC-1, sigma.formula=~t.stAgeC-1, nu.formula=~t.stAgeC-1, data=YXr, family=GG(), control=gamlss.control(n.cyc=60))	
	wbm4	<- gamlss( score.p~t.stAgeC+ns(ntPatient, df=2)-1, sigma.formula=~t.stAgeC+ns(ntPatient, df=2)-1, data=YXr, family=GA() )
	wbm5	<- gamlss( score.p~t.stAgeC+ns(ntPatient, df=2)-1, sigma.formula=~t.stAgeC+ns(ntPatient, df=2)-1, nu.formula=~t.stAgeC-1, data=YXr, family=GG(), control=gamlss.control(n.cyc=60) )
	wbm6	<- gamlss( score.p~t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient-1, data=YXr, family=GA() )
	wbm7	<- gamlss( score.p~t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient-1, nu.formula=~t.stAgeC-1, data=YXr, family=GG(), control=gamlss.control(n.cyc=200) )
	wbm8	<- gamlss( score.p~t.stAgeC+ntPatient:stageC-1, sigma.formula=~t.stAgeC+ntPatient:stageC-1, data=YXr, family=GA() )	
	wbm9	<- gamlss( score.p~t.stAgeC+ntPatient:t.stAgeC-1, sigma.formula=~t.stAgeC+ntPatient:t.stAgeC-1, data=YXr, family=GA() )
	wbm10	<- gamlss( score.p~t.stAgeC+ns(ntPatient, df=2):t.stAgeC-1, sigma.formula=~t.stAgeC+ns(ntPatient, df=2):t.stAgeC-1, data=YXr, family=GA() )
	wbm11	<- gamlss( score.p~t.stAgeC+ns(ntPatient, df=3):t.stAgeC-1, sigma.formula=~t.stAgeC+ns(ntPatient, df=3):t.stAgeC-1, data=YXr, family=GA() )
	wbm12	<- gamlss( score.p~t.AgeC+stageC+ns(ntPatient, df=2):t.stAgeC-1, sigma.formula=~t.AgeC+stageC+ns(ntPatient, df=2):t.stAgeC-1, data=YXr2, family=GA() )
	#
	#	given there are 'outliers', try understand probabilistic model fit
	#
	file	<- paste(indir,'/',gsub('\\.R','_scorePResidualsModels150911.pdf',infile),sep='')
	pdf(file=file, w=6, h=6)
	qqnorm(resid(wbm3), main='score.p~t.stAgeC-1 GA', xlab = "Theoretical Quantiles", 
			ylab = "Sample Quantiles", plot.it = TRUE, frame.plot = TRUE, 
			col = "darkgreen")
	lines(resid(wbm3), resid(wbm3), col = "red", lwd = 0.4, cex = 0.4)	
	qqnorm(resid(wbm2), main='score.p~t.stAgeC-1 GG', xlab = "Theoretical Quantiles", 
			ylab = "Sample Quantiles", plot.it = TRUE, frame.plot = TRUE, 
			col = "darkgreen")
	lines(resid(wbm2), resid(wbm2), col = "red", lwd = 0.4, cex = 0.4)
	qqnorm(resid(wbm4), main='score.p~t.stAgeC+ns(ntPatient, df=2)-1 GA', xlab = "Theoretical Quantiles", 
			ylab = "Sample Quantiles", plot.it = TRUE, frame.plot = TRUE, 
			col = "darkgreen")
	lines(resid(wbm4), resid(wbm4), col = "red", lwd = 0.4, cex = 0.4)
	qqnorm(resid(wbm5), main='score.p~t.stAgeC+ns(ntPatient, df=2)-1 GG', xlab = "Theoretical Quantiles", 
			ylab = "Sample Quantiles", plot.it = TRUE, frame.plot = TRUE, 
			col = "darkgreen")
	lines(resid(wbm5), resid(wbm5), col = "red", lwd = 0.4, cex = 0.4)	
	qqnorm(resid(wbm6), main='score.p~t.stAgeC+ntPatient-1 GA', xlab = "Theoretical Quantiles", 
			ylab = "Sample Quantiles", plot.it = TRUE, frame.plot = TRUE, 
			col = "darkgreen")
	lines(resid(wbm6), resid(wbm6), col = "red", lwd = 0.4, cex = 0.4)
	qqnorm(resid(wbm7), main='score.p~t.stAgeC+ntPatient-1 GG', xlab = "Theoretical Quantiles", 
			ylab = "Sample Quantiles", plot.it = TRUE, frame.plot = TRUE, 
			col = "darkgreen")
	lines(resid(wbm7), resid(wbm7), col = "red", lwd = 0.4, cex = 0.4)	
	qqnorm(resid(wbm8), main='score.p~t.stAgeC+ntPatient:stageC-1 GA', xlab = "Theoretical Quantiles", 
			ylab = "Sample Quantiles", plot.it = TRUE, frame.plot = TRUE, 
			col = "darkgreen")
	lines(resid(wbm8), resid(wbm8), col = "red", lwd = 0.4, cex = 0.4)		
	dev.off()
	#
	#	plot boxplots on residuals -- did not help much to see what s going on
	#
	YXe		<- copy(YXr)
	YXe[, rm3:= resid(wbm3)]
	YXe[, rm4:= resid(wbm4)]
	YXe[, rm6:= resid(wbm6)]
	ggplot(melt(YXe, measure.vars=c('rm3', 'rm4', 'rm6')), aes(x=variable, y=value)) + 
			geom_hline(yintercept=0, colour='grey80', lwd=1) + 
			geom_boxplot() + facet_grid(t.stAgeC~.) + theme_bw() + coord_flip()
	file	<- paste(indir,'/',gsub('\\.R','_scorePResidualsBoxPlotsModels150911.pdf',infile),sep='')
	ggsave(file=file, w=10, h=15)
	#
	#	plot predictions as function of ntPatient -- useful
	#
	YXe		<- copy(YXr)
	YXe[, p4:= predict(wbm4, what='mu', type='response')]
	YXe[, p8:= predict(wbm8, what='mu', type='response')]
	YXe[, p9:= predict(wbm9, what='mu', type='response')]
	YXe[, p10:= predict(wbm10, what='mu', type='response')]
	YXe[, p11:= predict(wbm11, what='mu', type='response')]
	setkey(YXe, t.stAgeC, ntPatient)
	YXe		<- melt(unique(subset(YXe, select=c(ntPatient, stageC, t.AgeC, t.stAgeC, p4, p8, p9, p10, p11))), measure.vars=c('p4', 'p8','p9', 'p10', 'p11'))
	ggplot(YXr, aes(x=ntPatient)) + #geom_point(aes(y=score.p), alpha=0.2) + 
			geom_smooth(aes(y=score.p), method='loess', degree=1, span=0.75, colour='black', fill='grey80') +
			geom_line(data=YXe, aes(y=value, colour=variable, group=variable, linetype=variable)) +
			facet_grid(t.AgeC~stageC, scales='free') + theme_bw() 
	file	<- paste(indir,'/',gsub('\\.R','_scorePvsNtPatient150911.pdf',infile),sep='')
	ggsave(file=file, w=8, h=8)
	#
	#	read out coefficients
	#
	cf		<- altvtp.readcoef(wbm3, 't.stAgeC', 'wbm3', 'score.p~t.stAgeC-1 GA')
	wb		<- copy(cf)
	cf		<- altvtp.readcoef(wbm4, 't.stAgeC', 'wbm4', 'score.p~t.stAgeC+\nns(ntPatient, df=2)-1 GA')
	wb		<- rbind(wb, cf)	
	cf		<- altvtp.readcoef(wbm8, 't.stAgeC', 'wbm8', 'score.p~t.stAgeC+\nntPatient:stageC-1 GA')
	wb		<- rbind(wb, cf)
	cf		<- altvtp.readcoef(wbm9, 't.stAgeC', 'wbm9', 'score.p~t.stAgeC+\nntPatient:t.stAgeC-1 GA')
	wb		<- rbind(wb, cf)
	cf		<- altvtp.readcoef(wbm10, 't.stAgeC', 'wbm10', 'score.p~t.stAgeC+\nns(ntPatient, df=2):t.stAgeC-1 GA')
	wb		<- rbind(wb, cf)
	tmp		<- melt(wb, id.vars=c('RSKF','MO_ID','MO','ST','RSKFbaseline'))
	set(tmp, NULL, 'value', tmp[, exp(value)])
	wb		<- dcast.data.table(tmp, MO_ID+MO+ST+RSKFbaseline+RSKF~variable, value.var='value')	
	#	for additive models, predict holding ntPatient fixed
	tmp		<- copy(YXr2)	
	tmp2	<- predict(wbm12, what='mu', type='response', se.fit=1)
	tmp[, MU:=tmp2$fit]
	tmp[, MUL:=tmp2$fit-2*tmp2$se.fit]
	tmp[, MUU:=tmp2$fit+2*tmp2$se.fit]
	setkey(tmp, t.stAgeC)
	tmp		<- unique(subset(tmp, ntPatient==4, c(t.stAgeC, MU, MUL, MUU)))
	tmp[, MO_ID:='wbm12']
	tmp[, MO:='score.p~t.AgeC+stageC+\nns(ntPatient, df=2):t.stAgeC-1']
	tmp[, ST:='t.stAgeC']
	tmp[, RSKFbaseline:= base.df['t.stAgeC',][,RSKFbaseline]]
	setnames(tmp, 't.stAgeC', 'RSKF')
	wb		<- rbind(wb, tmp, use.names=TRUE)
	
	tmp		<- subset(wb, !grepl('ntPatient',RSKF))	
	lvls	<- c( 	as.vector(sapply(		c('UA','UC','D'), function(x) paste(x,c('(-1,25]','(25,30]','(30,45]','(45,100]'),sep='_')		)),
			'T_(-1,100]','L_(-1,100]' )	
	set(tmp, NULL, 'RSKF', tmp[, factor(RSKF, levels=rev(lvls))])
	set(tmp, NULL, 'stageC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',1)])
	set(tmp, NULL, 't.AgeC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',2)])
	setnames(tmp, 'RSKF','t.stAgeC')
	tmp		<- merge(tmp, YXr[, list(mean=mean(score.p), ql= quantile(score.p,p=0.025), qu= quantile(score.p,p=0.975)), by=c('stageC','t.AgeC','t.stAgeC')], by=c('stageC','t.AgeC','t.stAgeC'))
	ggplot( tmp, aes(x=t.AgeC) ) +
			scale_y_continuous(breaks=seq(0,50,5), minor_breaks=seq(0,50,1)) +
			coord_cartesian(ylim=c(0,25)) +
			geom_boxplot(data=YXr, aes(y=100*score.p, fill=stageC), colour='grey80', alpha=0.9, outlier.shape=NA) +
			stat_summary(data=YXr, aes(y=100*score.p), fun.y=mean, colour="grey80", geom="point", shape=18, size=3, show_guide = FALSE) +
			geom_point(aes(y=100*MU)) + geom_errorbar(aes(ymin=100*MUL, ymax=100*MUU), width=0.5) + 
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(MO~stageC) 	
	file	<- paste(indir,'/',gsub('\\.R','_scorePMUByModel150911.pdf',infile),sep='')	
	ggsave(file=file, w=14,h=14)	

	#
	#	boxplots
	#	
	ggplot(YXr2, aes(x=t.AgeC, y=100*score.p, fill=t.AgeC)) + geom_boxplot(, colour='grey80', alpha=0.9,outlier.shape=NA) +
			scale_y_continuous(breaks=seq(0,50,5), minor_breaks=seq(0,50,1)) + coord_cartesian(ylim=c(0,15)) +
			stat_summary(data=YXr2, aes(y=100*score.p), fun.y=mean, colour="grey80", geom="point", shape=18, size=3, show_guide = FALSE) +
			facet_grid(~AgeC, scales='free') + theme_bw()
	file	<- paste(indir,'/',gsub('\\.R','_scorePBypAgeC150911.pdf',infile),sep='')	
	ggsave(file=file, w=14,h=8)	
	ggplot(YXr2, aes(x=t.AgeC, y=100*score.p, fill=t.AgeC)) + geom_boxplot(, colour='grey80', alpha=0.9,outlier.shape=NA) +
			scale_y_continuous(breaks=seq(0,50,5), minor_breaks=seq(0,50,1)) + coord_cartesian(ylim=c(0,25)) +
			facet_grid(stageC~AgeC, scales='free') + theme_bw()
	file	<- paste(indir,'/',gsub('\\.R','_scorePBypAgeCStage150911.pdf',infile),sep='')	
	ggsave(file=file, w=14,h=14)	
	
	
	pm	<- gamlss( score.p~p.AgeC+stageC+ns(ntPatient, df=2):t.stAgeC-1, sigma.formula=~p.AgeC+stageC+ns(ntPatient, df=2):t.stAgeC-1, data=YXr2, family=GA() )
			
}
######################################################################################
altvtp.model.150910<- function(YXc)
{
	
	YXr		<- subset(YXc, select=c('t.Patient','t','Patient', 'score.p','ntPatient','stageC','t.AgeC','t.stAgeC'))
	YXr2	<- subset(YXc, !stageC%in%c('T','L'),select=c('t.Patient','t','Patient', 'score.p','ntPatient','stageC','t.AgeC','t.stAgeC'))
	#	relevel so regression coefficients are contrasts of interest
	base.df	<- data.table(	ST=c('stageC','t.AgeC','t.stAgeC'), 
							RSKFbaseline=c('D','(30,45]','D_(30,45]'))
	setkey(base.df, ST)
	set(YXr, NULL, 'stageC', YXr[, relevel(stageC, ref=base.df['stageC',][,RSKFbaseline])])
	set(YXr, NULL, 't.AgeC', YXr[, relevel(t.AgeC, ref=base.df['t.AgeC',][,RSKFbaseline])])
	set(YXr, NULL, 't.stAgeC', YXr[, relevel(t.stAgeC, ref=base.df['t.stAgeC',][,RSKFbaseline])])
	#	fit contrasts 
	#	(we only do this to get the risk ratio and confidence intervals)		
	wcm3	<- gamlss( score.p~t.stAgeC, sigma.formula=~t.stAgeC, data=YXr, family=GA() )
	#tmp		<- altvp.gamlss.start(YXr, mu.formula='score.p~t.stAgeC+ntPatient', sigma.formula='~t.stAgeC+ntPatient', gamlss.limit=c(0, 0.001, 0.005, 0.006, 0.01 ), gamlss.family=GA(mu.link='inverse'),gamlss.control=gamlss.control(trace=0),gamlss.i.control=glim.control(bf.trace=0, glm.trace=0), verbose=FALSE)
	#wcm4	<- gamlss( score.p~t.stAgeC+ntPatient, sigma.formula=~t.stAgeC+ntPatient, data=YXr, family=GA(mu.link='inverse'), mu.start=tmp$mu, sigma.start=tmp$sigma)
	wcm4	<- gamlss( score.p~t.stAgeC+ns(ntPatient, df=2), sigma.formula=~t.stAgeC+ns(ntPatient, df=2), data=YXr, family=GA() )
	wcm5	<- gamlss( score.p~t.AgeC+stageC+ns(ntPatient, df=2), sigma.formula=~t.stAgeC+stageC+ns(ntPatient, df=2), data=YXr2, family=GA() )
	#	fit coefficients that correspond to stages 
	#	(we only do this to get the risk ratio and confidence intervals)	
	wbm3	<- gamlss( score.p~t.stAgeC-1, sigma.formula=~t.stAgeC-1, data=YXr, family=GA() )
	#	I don t get the inverse link to work!? The mu's are negative??
	#tmp		<- altvp.gamlss.start(YXr, mu.formula='score.p~t.stAgeC+ntPatient-1', sigma.formula='~t.stAgeC+ntPatient-1', gamlss.limit=c(0, 0.001, 0.005, 0.006, 0.01 ), gamlss.family=GA(mu.link='inverse'),gamlss.control=gamlss.control(trace=0),gamlss.i.control=glim.control(bf.trace=0, glm.trace=0), verbose=FALSE)
	#wbm4	<- gamlss( score.p~t.stAgeC+ntPatient-1, sigma.formula=~t.stAgeC+ntPatient-1, data=YXr, family=GA(mu.link='inverse'), mu.start=tmp$mu, sigma.start=tmp$sigma)
	wbm4	<- gamlss( score.p~t.stAgeC+ns(ntPatient, df=2)-1, sigma.formula=~t.stAgeC+ns(ntPatient, df=2)-1, data=YXr, family=GA() )
	wbm5	<- gamlss( score.p~t.AgeC+stageC+ns(ntPatient, df=2)-1, sigma.formula=~t.AgeC+stageC+ns(ntPatient, df=2)-1, data=YXr2, family=GA() )
	#	read out coefficients
	cf		<- altvtp.readcoef(wbm3, 't.stAgeC', 'wbm3', 'score.p~t.stAgeC-1')
	wb		<- copy(cf)
	cf		<- altvtp.readcoef(wbm4, 't.stAgeC', 'wbm4', 'score.p~t.stAgeC+ns(ntPatient, df=2)-1')
	wb		<- rbind(wb, cf)
	tmp		<- melt(wb, id.vars=c('RSKF','MO_ID','MO','ST','RSKFbaseline'))
	set(tmp, NULL, 'value', tmp[, exp(value)])
	wb		<- dcast.data.table(tmp, MO_ID+MO+ST+RSKFbaseline+RSKF~variable, value.var='value')
	#	read out contrasts
	cf		<- altvtp.readcoef(wcm3, 't.stAgeC', 'wcm3', 'score.p~t.stAgeC')
	wc		<- copy(cf)
	cf		<- altvtp.readcoef(wcm4, 't.stAgeC', 'wcm4', 'score.p~t.stAgeC+ns(ntPatient, df=2)')
	wc		<- rbind(wc, cf)
	tmp		<- melt(wc, id.vars=c('RSKF','MO_ID','MO','ST','RSKFbaseline'))
	set(tmp, NULL, 'value', tmp[, exp(value)])
	wc		<- dcast.data.table(tmp, MO_ID+MO+ST+RSKFbaseline+RSKF~variable, value.var='value')
	tmp		<- names(wc)[ grepl('MU',names(wc)) ]
	setnames(wc, tmp, gsub('MU','RR',tmp))
	set(wc, wc[, which(RSKF=='(Intercept)')],'RR',1.)
	set(wc, wc[, which(RSKF=='(Intercept)')],c('RRL','RRU'),NA_real_)
	tmp		<- wc[, which(RSKF=='(Intercept)')]	
	set(wc, tmp, 'RSKF', wc[tmp,RSKFbaseline])
	#	for additive models, predict holding ntPatient fixed
	tmp		<- copy(YXr2)	
	tmp2	<- predict(wbm5, what='mu', type='response', se.fit=1)
	tmp[, MU:=tmp2$fit]
	tmp[, MUL:=tmp2$fit-2*tmp2$se.fit]
	tmp[, MUU:=tmp2$fit+2*tmp2$se.fit]
	setkey(tmp, t.stAgeC)
	tmp		<- unique(subset(tmp, ntPatient==4, c(t.stAgeC, MU, MUL, MUU)))
	tmp[, MO_ID:='wbm5']
	tmp[, MO:='score.p~t.AgeC+stageC+ns(ntPatient, df=2)-1']
	tmp[, ST:='t.stAgeC']
	tmp[, RSKFbaseline:= base.df['t.stAgeC',][,RSKFbaseline]]
	setnames(tmp, 't.stAgeC', 'RSKF')
	wb		<- rbind(wb, tmp, use.names=TRUE)
	#	predict quantiles
	tmp		<- copy(YXr)	
	tmp2	<- predict(wbm3, what='mu', type='response', se.fit=1)
	tmp[, MU:=tmp2$fit]
	tmp[, MUL:=tmp2$fit-2*tmp2$se.fit]
	tmp[, MUU:=tmp2$fit+2*tmp2$se.fit]
	tmp2	<- predict(wbm3, what='sigma', type='response', se.fit=1)
	tmp[, SI:=tmp2$fit]
	tmp[, SIL:=tmp2$fit-2*tmp2$se.fit]
	tmp[, SIU:=tmp2$fit+2*tmp2$se.fit]
	tmp		<- unique(subset(tmp, ntPatient==4, c(t.stAgeC, MU, MUL, MUU, SI, SIL, SIU)))
	tmp[, Q50:=qGA(0.5, mu=MU, sigma=SI)]
	tmp[, MO_ID:='wbm3q']
	tmp[, MO:='Q score.p~t.stAgeC-1']
	tmp[, ST:='t.stAgeC']
	tmp[, RSKFbaseline:= base.df['t.stAgeC',][,RSKFbaseline]]
	set(tmp, NULL, c('MU','MUL','MUU','SI','SIL','SIU'), NULL)
	setnames(tmp, c('t.stAgeC','Q50'),c('RSKF','MU'))
	wb		<- rbind(wb, tmp, fill=TRUE, use.names=TRUE)
	#
	
	
	#	plot mus
	tmp		<- subset(wb, grepl('wbm3|wbm4|wbm5',MO_ID) & !grepl('ntPatient',RSKF))	
	lvls	<- c( 	as.vector(sapply(		c('UA','UC','D'), function(x) paste(x,c('(-1,25]','(25,30]','(30,45]','(45,100]'),sep='_')		)),
					'T_(-1,100]','L_(-1,100]' )	
	set(tmp, NULL, 'RSKF', tmp[, factor(RSKF, levels=rev(lvls))])
	set(tmp, NULL, 'stageC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',1)])
	set(tmp, NULL, 't.AgeC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',2)])
	setnames(tmp, 'RSKF','t.stAgeC')
	tmp		<- merge(tmp, YXr[, list(mean=mean(score.p), ql= quantile(score.p,p=0.025), qu= quantile(score.p,p=0.975)), by=c('stageC','t.AgeC','t.stAgeC')], by=c('stageC','t.AgeC','t.stAgeC'))
	ggplot( tmp, aes(x=factor(as.character(t.AgeC), levels=c('(-1,25]','(25,30]','(30,45]','(45,100]','(-1,100]'), labels=c('(-1,25]','(25,30]','(30,45]','(45,100]','(-1,100]')))) +
			scale_y_continuous(breaks=seq(0,50,5), minor_breaks=seq(0,50,1)) +
			coord_cartesian(ylim=c(0,25)) +
			geom_boxplot(data=YXr, aes(y=100*score.p, fill=stageC), colour='grey80', alpha=0.9, outlier.shape=NA) +
			stat_summary(data=YXr, aes(y=100*score.p), fun.y=mean, colour="black", geom="point", shape=18, size=3, show_guide = FALSE) +
			geom_point(aes(y=100*MU)) + geom_errorbar(aes(ymin=100*MUL, ymax=100*MUU), width=0.5) + 
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(MO~stageC) 	
	file	<- paste(indir,'/',gsub('\\.R','_scorePMUByModel150910.pdf',infile),sep='')	
	ggsave(file=file, w=14,h=12)	
	#	plot mu vs empirical mean
	ggplot(subset(tmp, MO_ID=='wbm3'), aes(x=mean, xmin=ql, xmax=qu, y=MU, ymin=MUL, ymax=MUU, colour=t.AgeC)) + 
			geom_abline(intercept=0, slope=1) + 
			geom_point() + geom_errorbar() + geom_errorbarh() +			
			theme_bw() + theme(legend.position='bottom') +
			facet_wrap(~stageC, ncol=5, scales='free')
	file	<- paste(indir,'/',gsub('\\.R','_scorePMUvsMEAN150910.pdf',infile),sep='')	
	ggsave(file=file, w=15,h=5)		
	#	plot contrasts
	tmp		<- subset(wc, grepl('wcm3|wcm4',MO_ID) & !grepl('ntPatient',RSKF))	
	lvls	<- c( 	as.vector(sapply(		c('UA','UC','D'), function(x) paste(x,c('(-1,25]','(25,30]','(30,45]','(45,100]'),sep='_')		)),
					'T_(-1,100]','L_(-1,100]' )	
	set(tmp, NULL, 'RSKF', tmp[, factor(RSKF, levels=rev(lvls))])
	set(tmp, NULL, 'stageC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',1)])
	set(tmp, NULL, 't.AgeC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',2)])
	ggplot( tmp, aes(x=t.AgeC, y=RR, ymin=RRL, ymax=RRU, colour=stageC)) +
			scale_y_continuous(breaks=seq(0,5,0.5), minor_breaks=seq(0,5,0.1)) + geom_hline(y=1, colour='grey50', lwd=1) +
			geom_point() + geom_errorbar(width=0.5) + 
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(MO~stageC) 
	file	<- paste(indir,'/',gsub('\\.R','_scorePRRByModel150910.pdf',infile),sep='')	
	ggsave(file=file, w=12,h=12)		
	ggplot( tmp, aes(y=MO, x=RR, xmin=RRL, xmax=RRU, colour=RSKF) ) +
			scale_x_continuous(breaks=seq(0,5,0.2)) + geom_vline(x=1, colour='grey50', lwd=1) +
			geom_point() + geom_errorbarh(height=0.5) + theme_bw() + theme(legend.position='bottom') +
			facet_grid(RSKF~.) +
			guides(colour=guide_legend(ncol=4))
	file	<- paste(indir,'/',gsub('\\.R','_scorePRRByModel150910_long.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=20)	
	
	list(data=YXr, wb=wb, wc=wc, wbm1=wbm1, wbm2=wbm2, wbm3=wbm3, wbm4=wbm4, wbm5=wbm5, wcm1=wcm1, wcm2=wcm2, wcm3=wcm3, wcm4=wcm4, wcm5=wcm5)
}
######################################################################################
altvtp.model.150902<- function(YXc)
{
	YXr		<- subset(YXc, !stageC%in%c('T','L'),select=c('t.Patient','t','Patient', 'score.p','ntPatient','stageC','t.AgeC','t.stAgeC'))
	YXr		<- subset(YXc, select=c('t.Patient','t','Patient', 'score.p','ntPatient','stageC','t.AgeC','t.stAgeC'))
	#	relevel so regression coefficients are contrasts of interest
	base.df	<- data.table(	ST=c('stageC','t.AgeC','t.stAgeC'), 
							RSKFbaseline=c('D','(30,45]','D_(30,45]'))
	setkey(base.df, ST)
	set(YXr, NULL, 'stageC', YXr[, relevel(stageC, ref=base.df['stageC',][,RSKFbaseline])])
	set(YXr, NULL, 't.AgeC', YXr[, relevel(t.AgeC, ref=base.df['t.AgeC',][,RSKFbaseline])])
	set(YXr, NULL, 't.stAgeC', YXr[, relevel(t.stAgeC, ref=base.df['t.stAgeC',][,RSKFbaseline])])
	#	fit contrasts 
	#	(we only do this to get the risk ratio and confidence intervals)	
	wcm1	<- gamlss( score.p~stageC+t.AgeC, sigma.formula=~stageC+t.AgeC, data=YXr, family=GA() )
	wcm2	<- altvp.gamlss.guide(YXr, mu.formula='score.p~stageC+t.AgeC+ntPatient', sigma.formula='~stageC+t.AgeC+ntPatient', gamlss.limit=c(0, 0.001, 0.005, 0.006, 0.01 ), gamlss.family=GA(mu.link='inverse'),gamlss.control=gamlss.control(trace=0),gamlss.i.control=glim.control(bf.trace=0, glm.trace=0), verbose=FALSE) 
	wcm3	<- gamlss( score.p~t.stAgeC, sigma.formula=~t.stAgeC, data=YXr, family=GA() )
	wcm4	<- altvp.gamlss.guide(YXr, mu.formula='score.p~t.stAgeC+ntPatient', sigma.formula='~t.stAgeC+ntPatient', gamlss.limit=c(0, 0.001, 0.005, 0.006, 0.01 ), gamlss.family=GA(mu.link='inverse'),gamlss.control=gamlss.control(trace=0),gamlss.i.control=glim.control(bf.trace=0, glm.trace=0), verbose=FALSE)	
	#	fit coefficients that correspond to stages 
	#	(we only do this to get the risk ratio and confidence intervals)	
	wbm1	<- gamlss( score.p~t.AgeC+stageC-1, sigma.formula=~t.AgeC+stageC-1, data=YXr, family=GA() )
	wbm2	<- altvp.gamlss.guide(YXr, mu.formula='score.p~t.AgeC+stageC+ntPatient-1', sigma.formula='~t.AgeC+stageC+ntPatient-1', gamlss.limit=c(0, 0.001, 0.005, 0.006, 0.01 ), gamlss.family=GA(mu.link='inverse'),gamlss.control=gamlss.control(trace=0),gamlss.i.control=glim.control(bf.trace=0, glm.trace=0), verbose=FALSE) 
	wbm3	<- gamlss( score.p~t.stAgeC-1, sigma.formula=~t.stAgeC-1, data=YXr, family=GA() )
	wbm4	<- altvp.gamlss.guide(YXr, mu.formula='score.p~t.stAgeC+ntPatient-1', sigma.formula='~t.stAgeC+ntPatient-1', gamlss.limit=c(0, 0.001, 0.005, 0.006, 0.01 ), gamlss.family=GA(mu.link='inverse'),gamlss.control=gamlss.control(trace=0),gamlss.i.control=glim.control(bf.trace=0, glm.trace=0), verbose=FALSE)
	
	
	#	read out coefficients
	tmp		<- coef(wbm1)
	tmp2	<- 't.AgeC'
	tmp		<- tmp[ grepl(tmp2,names(tmp)) ]
	cf		<- data.table(MO_ID='wbm1', MO='score.p~t.AgeC+stageC-1', ST=tmp2, RSKF=names(tmp), RSKFbaseline=base.df[tmp2,][,RSKFbaseline], MU=tmp)
	tmp		<- confint(wbm1)
	tmp		<- tmp[ grepl(tmp2,names(tmp)) ]
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MUL= tmp[,1], MUU= tmp[,2]), by='RSKF')
	tmp		<- confint(wbm1, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MULr= tmp[,1], MUUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.AgeC','',RSKF)])
	wb		<- copy(cf)
	#
	tmp		<- coef(wbm2)
	cf		<- data.table(MO='wbm2', ST='stageC', RSKF=names(tmp), RSKFbaseline=base.df['stageC',][,RSKFbaseline], MU=tmp)
	tmp		<- confint(wbm2)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MUL= tmp[,1], MUU= tmp[,2]), by='RSKF')
	tmp		<- confint(wbm2, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MULr= tmp[,1], MUUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('stageC','',RSKF)])
	wb		<- rbind(wb, cf)
	#
	tmp		<- coef(wbm3)
	cf		<- data.table(MO='wbm3', ST='t.stAgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.stAgeC',][,RSKFbaseline], MU=tmp)
	tmp		<- confint(wbm3)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MUL= tmp[,1], MUU= tmp[,2]), by='RSKF')
	tmp		<- confint(wbm3, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MULr= tmp[,1], MUUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.stAgeC','',RSKF)])
	wb		<- rbind(wb, cf)	
	#
	tmp		<- coef(wbm4)
	cf		<- data.table(MO='wbm4', ST='t.stAgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.stAgeC',][,RSKFbaseline], MU=tmp)
	tmp		<- confint(wbm4)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MUL= tmp[,1], MUU= tmp[,2]), by='RSKF')
	tmp		<- confint(wbm4, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MULr= tmp[,1], MUUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.stAgeC','',RSKF)])
	wb		<- rbind(wb, cf)
	#
	tmp		<- coef(wbm5)
	cf		<- data.table(MO='wbm5', ST='t.stAgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.stAgeC',][,RSKFbaseline], MU=tmp)
	tmp		<- confint(wbm5)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MUL= tmp[,1], MUU= tmp[,2]), by='RSKF')
	tmp		<- confint(wbm5, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MULr= tmp[,1], MUUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.stAgeC','',RSKF)])
	wb		<- rbind(wb, cf)	
	tmp		<- melt(wb, id.vars=c('RSKF','MO','ST','RSKFbaseline'))
	set(tmp, NULL, 'value', tmp[, exp(value)])
	wb		<- dcast.data.table(tmp, MO+ST+RSKFbaseline+RSKF~variable, value.var='value')
	#
	#	read out contrasts
	#
	tmp		<- coef(wcm1)
	cf		<- data.table(MO='wcm1', ST='t.AgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.AgeC',][,RSKFbaseline], CNTR=tmp)
	tmp		<- confint(wcm1)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRL= tmp[,1], CNTRU= tmp[,2]), by='RSKF')
	tmp		<- confint(wcm1, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRLr= tmp[,1], CNTRUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.AgeC','',RSKF)])
	wc		<- copy(cf)
	#
	tmp		<- coef(wcm2)
	cf		<- data.table(MO='wcm2', ST='stageC', RSKF=names(tmp), RSKFbaseline=base.df['stageC',][,RSKFbaseline], CNTR=tmp)
	tmp		<- confint(wcm2)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRL= tmp[,1], CNTRU= tmp[,2]), by='RSKF')
	tmp		<- confint(wcm2, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRLr= tmp[,1], CNTRUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('stageC','',RSKF)])
	wc		<- rbind(wc, cf)
	#
	tmp		<- coef(wcm3)
	cf		<- data.table(MO='wcm3', ST='t.stAgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.stAgeC',][,RSKFbaseline], CNTR=tmp)
	tmp		<- confint(wcm3)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRL= tmp[,1], CNTRU= tmp[,2]), by='RSKF')
	tmp		<- confint(wcm3, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRLr= tmp[,1], CNTRUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.stAgeC','',RSKF)])
	wc		<- rbind(wc, cf)	
	#
	tmp		<- coef(wcm4)
	cf		<- data.table(MO='wcm4', ST='t.stAgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.stAgeC',][,RSKFbaseline], CNTR=tmp)
	tmp		<- confint(wcm4)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRL= tmp[,1], CNTRU= tmp[,2]), by='RSKF')
	tmp		<- confint(wcm4, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRLr= tmp[,1], CNTRUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.stAgeC','',RSKF)])
	wc		<- rbind(wc, cf)
	#
	tmp		<- coef(wcm5)
	cf		<- data.table(MO='wcm5', ST='t.stAgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.stAgeC',][,RSKFbaseline], CNTR=tmp)
	tmp		<- confint(wcm5)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRL= tmp[,1], CNTRU= tmp[,2]), by='RSKF')
	tmp		<- confint(wcm5, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRLr= tmp[,1], CNTRUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.stAgeC','',RSKF)])
	wc		<- rbind(wc, cf)
	#	transform to risk ratios
	tmp		<- melt(wc, id.vars=c('RSKF','MO','ST','RSKFbaseline'))
	set(tmp, NULL, 'value', tmp[, exp(value)])
	wc		<- dcast.data.table(tmp, MO+ST+RSKFbaseline+RSKF~variable, value.var='value')
	tmp		<- names(wc)[ grepl('CNTR',names(wc)) ]
	setnames(wc, tmp, gsub('CNTR','RR',tmp))
	set(wc, wc[, which(RSKF=='(Intercept)')],'RR',1.)
	set(wc, wc[, which(RSKF=='(Intercept)')],c('RRL','RRU','RRLr','RRUr'),NA_real_)
	tmp		<- wc[, which(RSKF=='(Intercept)')]	
	set(wc, tmp, 'RSKF', wc[tmp,RSKFbaseline])
	
	#	plot mus
	tmp		<- subset(wb, grepl('wbm3|wbm4|wbm5',MO) & RSKF!='ntPatient')	
	lvls	<- c( 	as.vector(sapply(		c('UA','UC','D'), function(x) paste(x,c('(-1,25]','(25,30]','(30,45]','(45,100]'),sep='_')		)),
					'T_(-1,100]','L_(-1,100]' )	
	set(tmp, NULL, 'RSKF', tmp[, factor(RSKF, levels=rev(lvls))])
	set(tmp, NULL, 'stageC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',1)])
	set(tmp, NULL, 't.AgeC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',2)])
	setnames(tmp, 'RSKF','t.stAgeC')
	tmp		<- merge(tmp, YXr[, list(mean=mean(score.p), ql= quantile(score.p,p=0.025), qu= quantile(score.p,p=0.975)), by=c('stageC','t.AgeC','t.stAgeC')], by=c('stageC','t.AgeC','t.stAgeC'))
	ggplot( tmp, aes(x=t.AgeC, y=100*MU, ymin=100*MUL, ymax=100*MUU, colour=stageC)) +
			scale_y_continuous(breaks=seq(0,50,5), minor_breaks=seq(0,50,1)) + 
			geom_point() + geom_errorbar(width=0.5) + 
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(MO~stageC) 
	file	<- paste(indir,'/',gsub('\\.R','_scorePMUByModel150902.pdf',infile),sep='')	
	ggsave(file=file, w=12,h=12)	
	#	plot mu vs empirical mean
	ggplot(subset(tmp, MO=='wbm3'), aes(x=mean, xmin=ql, xmax=qu, y=MU, ymin=MUL, ymax=MUU, colour=t.AgeC)) + geom_abline(intercept=0, slope=1) + 
			geom_point() + geom_errorbar() + geom_errorbarh() +			
			theme_bw() + 
			facet_wrap(~stageC, ncol=5, scales='free')
	file	<- paste(indir,'/',gsub('\\.R','_scorePMUvsMEAN150902.pdf',infile),sep='')	
	ggsave(file=file, w=15,h=5)		
	#	plot contrasts
	tmp		<- subset(wc, grepl('wcm3|wcm4|wcm5',MO) & RSKF!='ntPatient')	
	lvls	<- c( 	as.vector(sapply(		c('UA','UC','D'), function(x) paste(x,c('(-1,25]','(25,30]','(30,45]','(45,100]'),sep='_')		)),
					'T_(-1,100]','L_(-1,100]' )	
	set(tmp, NULL, 'RSKF', tmp[, factor(RSKF, levels=rev(lvls))])
	set(tmp, NULL, 'stageC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',1)])
	set(tmp, NULL, 't.AgeC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',2)])
	ggplot( tmp, aes(x=t.AgeC, y=RR, ymin=RRL, ymax=RRU, colour=stageC)) +
			scale_y_continuous(breaks=seq(0,5,0.5), minor_breaks=seq(0,5,0.1)) + geom_hline(y=1, colour='grey50', lwd=1) +
			geom_point() + geom_errorbar(width=0.5) + 
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(MO~stageC) 
	file	<- paste(indir,'/',gsub('\\.R','_scorePRRByModel150902.pdf',infile),sep='')	
	ggsave(file=file, w=12,h=12)	
	
	ggplot( tmp, aes(y=MO, x=RR, xmin=RRL, xmax=RRU, colour=RSKF) ) +
			scale_x_continuous(breaks=seq(0,5,0.2)) + geom_vline(x=1, colour='grey50', lwd=1) +
			geom_point() + geom_errorbarh(height=0.5) + theme_bw() + theme(legend.position='bottom') +
			facet_grid(RSKF~.) +
			guides(colour=guide_legend(ncol=4))
	file	<- paste(indir,'/',gsub('\\.R','_scorePRRByModel150902.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=20)	
	
	list(data=YXr, wb=wb, wc=wc, wbm1=wbm1, wbm2=wbm2, wbm3=wbm3, wbm4=wbm4, wbm5=wbm5, wcm1=wcm1, wcm2=wcm2, wcm3=wcm3, wcm4=wcm4, wcm5=wcm5)
}
######################################################################################
altvtp.model.150730<- function(YXc)
{
	YXr		<- subset(YXc, select=c('t.Patient','t','Patient', 'score.p','ntPatient','stageC','t.AgeC','t.stAgeC'))
	YXr[, ntPatientn:= as.numeric(as.character(ntPatient))]
	#	relevel so regression coefficients are contrasts of interest
	base.df	<- data.table(	ST=c('stageC','t.AgeC','t.stAgeC'), 
							RSKFbaseline=c('D','(30,45]','D_(30,45]'))
	setkey(base.df, ST)
	set(YXr, NULL, 'stageC', YXr[, relevel(stageC, ref=base.df['stageC',][,RSKFbaseline])])
	set(YXr, NULL, 't.AgeC', YXr[, relevel(t.AgeC, ref=base.df['t.AgeC',][,RSKFbaseline])])
	set(YXr, NULL, 't.stAgeC', YXr[, relevel(t.stAgeC, ref=base.df['t.stAgeC',][,RSKFbaseline])])
	#	fit contrasts 
	#	(we only do this to get the risk ratio and confidence intervals)
	wcm1	<- gamlss( score.p~t.AgeC+ntPatientn, sigma.formula=~t.AgeC+ntPatientn, data=YXr, family=GA() )
	wcm2	<- gamlss( score.p~stageC+ntPatientn, sigma.formula=~stageC+ntPatientn, data=YXr, family=GA() )	
	wcm3	<- gamlss( score.p~t.stAgeC, sigma.formula=~t.stAgeC, data=YXr, family=GA() )
	wcm4	<- gamlss( score.p~t.stAgeC+ntPatientn, sigma.formula=~ntPatientn, data=YXr, family=GA() )
	wcm5	<- gamlss( score.p~t.stAgeC+ntPatientn, sigma.formula=~t.stAgeC+ntPatientn, data=YXr, family=GA() )
	#	fit coefficients that correspond to stages 
	#	(we only do this to get the risk ratio and confidence intervals)
	wbm1	<- gamlss( score.p~t.AgeC+ntPatientn-1, sigma.formula=~t.AgeC+ntPatientn-1, data=YXr, family=GA() )
	wbm2	<- gamlss( score.p~stageC+ntPatientn-1, sigma.formula=~stageC+ntPatientn-1, data=YXr, family=GA() )	
	#wbm3	<- gamlss( score.p~t.stAgeC-1, sigma.formula=~t.stAgeC-1, data=YXr, family=GA() )
	wbm3	<- gamlss( score.p~t.stAgeC-1, data=YXr, family=GA() )
	wbm4	<- gamlss( score.p~t.stAgeC+ntPatientn-1, sigma.formula=~ntPatientn-1, data=YXr, family=GA() )
	wbm5	<- gamlss( score.p~t.stAgeC+ntPatientn-1, sigma.formula=~t.stAgeC+ntPatientn-1, data=YXr, family=GA() )
	#	read out coefficients
	tmp		<- coef(wbm1)
	cf		<- data.table(MO='wbm1', ST='t.AgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.AgeC',][,RSKFbaseline], MU=tmp)
	tmp		<- confint(wbm1)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MUL= tmp[,1], MUU= tmp[,2]), by='RSKF')
	tmp		<- confint(wbm1, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MULr= tmp[,1], MUUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.AgeC','',RSKF)])
	wb		<- copy(cf)
	#
	tmp		<- coef(wbm2)
	cf		<- data.table(MO='wbm2', ST='stageC', RSKF=names(tmp), RSKFbaseline=base.df['stageC',][,RSKFbaseline], MU=tmp)
	tmp		<- confint(wbm2)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MUL= tmp[,1], MUU= tmp[,2]), by='RSKF')
	tmp		<- confint(wbm2, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MULr= tmp[,1], MUUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('stageC','',RSKF)])
	wb		<- rbind(wb, cf)
	#
	tmp		<- coef(wbm3)
	cf		<- data.table(MO='wbm3', ST='t.stAgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.stAgeC',][,RSKFbaseline], MU=tmp)
	tmp		<- confint(wbm3)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MUL= tmp[,1], MUU= tmp[,2]), by='RSKF')
	tmp		<- confint(wbm3, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MULr= tmp[,1], MUUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.stAgeC','',RSKF)])
	wb		<- rbind(wb, cf)	
	#
	tmp		<- coef(wbm4)
	cf		<- data.table(MO='wbm4', ST='t.stAgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.stAgeC',][,RSKFbaseline], MU=tmp)
	tmp		<- confint(wbm4)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MUL= tmp[,1], MUU= tmp[,2]), by='RSKF')
	tmp		<- confint(wbm4, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MULr= tmp[,1], MUUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.stAgeC','',RSKF)])
	wb		<- rbind(wb, cf)
	#
	tmp		<- coef(wbm5)
	cf		<- data.table(MO='wbm5', ST='t.stAgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.stAgeC',][,RSKFbaseline], MU=tmp)
	tmp		<- confint(wbm5)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MUL= tmp[,1], MUU= tmp[,2]), by='RSKF')
	tmp		<- confint(wbm5, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), MULr= tmp[,1], MUUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.stAgeC','',RSKF)])
	wb		<- rbind(wb, cf)	
	tmp		<- melt(wb, id.vars=c('RSKF','MO','ST','RSKFbaseline'))
	set(tmp, NULL, 'value', tmp[, exp(value)])
	wb		<- dcast.data.table(tmp, MO+ST+RSKFbaseline+RSKF~variable, value.var='value')
	#
	#	read out contrasts
	#
	tmp		<- coef(wcm1)
	cf		<- data.table(MO='wcm1', ST='t.AgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.AgeC',][,RSKFbaseline], CNTR=tmp)
	tmp		<- confint(wcm1)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRL= tmp[,1], CNTRU= tmp[,2]), by='RSKF')
	tmp		<- confint(wcm1, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRLr= tmp[,1], CNTRUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.AgeC','',RSKF)])
	wc		<- copy(cf)
	#
	tmp		<- coef(wcm2)
	cf		<- data.table(MO='wcm2', ST='stageC', RSKF=names(tmp), RSKFbaseline=base.df['stageC',][,RSKFbaseline], CNTR=tmp)
	tmp		<- confint(wcm2)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRL= tmp[,1], CNTRU= tmp[,2]), by='RSKF')
	tmp		<- confint(wcm2, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRLr= tmp[,1], CNTRUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('stageC','',RSKF)])
	wc		<- rbind(wc, cf)
	#
	tmp		<- coef(wcm3)
	cf		<- data.table(MO='wcm3', ST='t.stAgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.stAgeC',][,RSKFbaseline], CNTR=tmp)
	tmp		<- confint(wcm3)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRL= tmp[,1], CNTRU= tmp[,2]), by='RSKF')
	tmp		<- confint(wcm3, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRLr= tmp[,1], CNTRUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.stAgeC','',RSKF)])
	wc		<- rbind(wc, cf)	
	#
	tmp		<- coef(wcm4)
	cf		<- data.table(MO='wcm4', ST='t.stAgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.stAgeC',][,RSKFbaseline], CNTR=tmp)
	tmp		<- confint(wcm4)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRL= tmp[,1], CNTRU= tmp[,2]), by='RSKF')
	tmp		<- confint(wcm4, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRLr= tmp[,1], CNTRUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.stAgeC','',RSKF)])
	wc		<- rbind(wc, cf)
	#
	tmp		<- coef(wcm5)
	cf		<- data.table(MO='wcm5', ST='t.stAgeC', RSKF=names(tmp), RSKFbaseline=base.df['t.stAgeC',][,RSKFbaseline], CNTR=tmp)
	tmp		<- confint(wcm5)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRL= tmp[,1], CNTRU= tmp[,2]), by='RSKF')
	tmp		<- confint(wcm5, robust=TRUE)
	cf		<- merge(cf, data.table(RSKF= rownames(tmp), CNTRLr= tmp[,1], CNTRUr= tmp[,2]), by='RSKF')
	set(cf, NULL, 'RSKF', cf[,gsub('t.stAgeC','',RSKF)])
	wc		<- rbind(wc, cf)
	#	transform to risk ratios
	tmp		<- melt(wc, id.vars=c('RSKF','MO','ST','RSKFbaseline'))
	set(tmp, NULL, 'value', tmp[, exp(value)])
	wc		<- dcast.data.table(tmp, MO+ST+RSKFbaseline+RSKF~variable, value.var='value')
	tmp		<- names(wc)[ grepl('CNTR',names(wc)) ]
	setnames(wc, tmp, gsub('CNTR','RR',tmp))
	set(wc, wc[, which(RSKF=='(Intercept)')],'RR',1.)
	set(wc, wc[, which(RSKF=='(Intercept)')],c('RRL','RRU','RRLr','RRUr'),NA_real_)
	tmp		<- wc[, which(RSKF=='(Intercept)')]	
	set(wc, tmp, 'RSKF', wc[tmp,RSKFbaseline])
	
	#	plot mus
	tmp		<- subset(wb, grepl('wbm3|wbm4|wbm5',MO) & RSKF!='ntPatientn')	
	lvls	<- as.vector(sapply(		c('UA','UC','D','T','L'), function(x) paste(x,c('(-1,25]','(25,30]','(30,45]','(45,100]'),sep='_')		))
	set(tmp, NULL, 'RSKF', tmp[, factor(RSKF, levels=rev(lvls))])
	set(tmp, NULL, 'stageC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',1)])
	set(tmp, NULL, 't.AgeC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',2)])
	setnames(tmp, 'RSKF','t.stAgeC')
	tmp		<- merge(tmp, YXr[, list(mean=mean(score.p), ql= quantile(score.p,p=0.025), qu= quantile(score.p,p=0.975)), by=c('stageC','t.AgeC','t.stAgeC')], by=c('stageC','t.AgeC','t.stAgeC'))
	ggplot( tmp, aes(x=t.AgeC, y=100*MU, ymin=100*MUL, ymax=100*MUU, colour=stageC)) +
			scale_y_continuous(breaks=seq(0,50,5), minor_breaks=seq(0,50,1)) + 
			geom_point() + geom_errorbar(width=0.5) + 
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(MO~stageC) 
	file	<- paste(indir,'/',gsub('\\.R','_scorePMUByModel2.pdf',infile),sep='')	
	ggsave(file=file, w=12,h=12)	
	#	plot mu vs empirical mean
	ggplot(subset(tmp, MO=='wbm3'), aes(x=mean, xmin=ql, xmax=qu, y=MU, ymin=MUL, ymax=MUU, colour=t.AgeC)) + geom_abline(intercept=0, slope=1) + 
			geom_point() + geom_errorbar() + geom_errorbarh() +			
			theme_bw() + 
			facet_wrap(~stageC, ncol=5, scales='free')
	file	<- paste(indir,'/',gsub('\\.R','_scorePMUvsMEAN.pdf',infile),sep='')	
	ggsave(file=file, w=15,h=5)		
	#	plot contrasts
	tmp		<- subset(wc, grepl('wcm3|wcm4|wcm5',MO) & RSKF!='ntPatientn')	
	lvls	<- as.vector(sapply(		c('UA','UC','D','T','L'), function(x) paste(x,c('(-1,25]','(25,30]','(30,45]','(45,100]'),sep='_')		))
	set(tmp, NULL, 'RSKF', tmp[, factor(RSKF, levels=rev(lvls))])
	set(tmp, NULL, 'stageC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',1)])
	set(tmp, NULL, 't.AgeC', tmp[, sapply(strsplit(as.character(RSKF), '_', fixed=TRUE),'[[',2)])
	ggplot( tmp, aes(x=t.AgeC, y=RR, ymin=RRL, ymax=RRU, colour=stageC)) +
			scale_y_continuous(breaks=seq(0,5,0.5), minor_breaks=seq(0,5,0.1)) + geom_hline(y=1, colour='grey50', lwd=1) +
			geom_point() + geom_errorbar(width=0.5) + 
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(MO~stageC) 
	file	<- paste(indir,'/',gsub('\\.R','_scorePRRByModel2.pdf',infile),sep='')	
	ggsave(file=file, w=12,h=12)	
	
	ggplot( tmp, aes(y=MO, x=RR, xmin=RRL, xmax=RRU, colour=RSKF) ) +
			scale_x_continuous(breaks=seq(0,5,0.2)) + geom_vline(x=1, colour='grey50', lwd=1) +
			geom_point() + geom_errorbarh(height=0.5) + theme_bw() + theme(legend.position='bottom') +
			facet_grid(RSKF~.) +
			guides(colour=guide_legend(ncol=4))
	file	<- paste(indir,'/',gsub('\\.R','_scorePRRByModel.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=20)	
	
	list(data=YXr, wb=wb, wc=wc, wbm1=wbm1, wbm2=wbm2, wbm3=wbm3, wbm4=wbm4, wbm5=wbm5, wcm1=wcm1, wcm2=wcm2, wcm3=wcm3, wcm4=wcm4, wcm5=wcm5)
}
######################################################################################
altvtp.explore.distribution<- function(YXc, indir, infile)
{
	YXr		<- subset(YXc, select=c('score.p','ntPatient'))	
	YXr[, ntPatientf:= factor(ntPatient)]
	setkey(YXr, ntPatient)
	#	Exp has default link log	
	mExp	<- gamlss( score.p~1, data=YXr, family=EXP() )
	mExpNt	<- gamlss( score.p~ntPatientf-1, data=YXr, family=EXP() )
	mExpNtL	<- gamlss( score.p~ntPatient, data=YXr, family=EXP() )
	mExpNtQ	<- gamlss( score.p~ns(ntPatient, df=2), data=YXr, family=EXP() )
	#	Gamma has default links log log
	mGA		<- gamlss( score.p~1, data=YXr, family=GA() )
	mGANt	<- gamlss( score.p~ntPatientf-1, data=YXr, family=GA() )
	mGANtS	<- gamlss( score.p~ntPatientf-1, sigma.formula=~ntPatientf-1, data=YXr, family=GA() )
	mGANtSL	<- gamlss( score.p~ntPatient, sigma.formula=~ntPatient, data=YXr, family=GA() )
	mGANtSLI<- gamlss( score.p~ntPatient, sigma.formula=~ntPatient, data=YXr, family=GA(mu.link='inverse') )
	mGANtSQ	<- gamlss( score.p~ns(ntPatient, df=2), sigma.formula=~ns(ntPatient, df=2), data=YXr, family=GA() )
	mGANtSQI<- gamlss( score.p~ns(ntPatient, df=2), sigma.formula=~ns(ntPatient, df=2), data=YXr, family=GA(mu.link='identity') )
	#	Beta with links log log
	mBEL	<- gamlss( score.p~1, data=YXr, family=BE(mu.link='log'), n.cyc = 40 )
	mBELNt	<- gamlss( score.p~ntPatientf-1, data=YXr, family=BE(mu.link='log'), n.cyc = 40 )
	mBELNtS	<- gamlss( score.p~ntPatientf-1, sigma.formula=~ntPatientf-1, data=YXr, family=BE(mu.link='log', sigma.link='log'), n.cyc = 80 )
	mBELNtSL<- gamlss( score.p~ntPatient, sigma.formula=~ntPatient, data=YXr, family=BE(mu.link='log', sigma.link='log'), n.cyc = 80 )
	mBELNtSQ<- gamlss( score.p~ns(ntPatient, df=2), sigma.formula=~ns(ntPatient, df=2), data=YXr, family=BE(mu.link='log', sigma.link='log'), n.cyc = 80 )
	#	
	z		<- list(mExp, mExpNt, mExpNtL, mExpNtQ, mGA, mGANt, mGANtS, mGANtSL, mGANtSLI, mGANtSQ, mGANtSQI, mBEL, mBELNt, mBELNtS, mBELNtSL, mBELNtSQ)
	names(z)<- c('mExp','mExpNt','mExpNtL','mExpNtQ', 'mGA', 'mGANt', 'mGANtS', 'mGANtSL', 'mGANtSLI', 'mGANtSQ', 'mGANtSQI', 'mBEL', 'mBELNt', 'mBELNtS', 'mBELNtSL', 'mBELNtSQ')	
	#z		<- list(mExpNt)
	#names(z)<- c('mExpNt')	
	mpars	<- do.call('rbind',lapply( seq_along(z), function(i){
				tmp		<- copy(YXr)
				mu		<- predict(z[[i]], type='response', what='mu')
				if(grepl('Exp',names(z)[i]))
					sigma	<- NA_real_
				if(!grepl('Exp',names(z)[i]))
					sigma	<- predict(z[[i]], type='response', what='sigma')				
				tmp[, score.mu:= mu]
				tmp[, score.sigma:= sigma]
				tmp[, TYPE:= names(z)[i]]
				setkey(tmp, ntPatient)
				unique(tmp)
			}))
	#	
	mpars[, score.sd:=NA_real_]
	tmp		<- mpars[, which(grepl('Exp',TYPE))]
	set(mpars, tmp, 'score.sd', mpars[tmp, score.mu])
	tmp		<- mpars[, which(grepl('GA',TYPE))]
	set(mpars, tmp, 'score.sd', mpars[tmp, score.mu*score.sigma])
	tmp		<- mpars[, which(grepl('BE',TYPE))]
	set(mpars, tmp, 'score.sd', mpars[tmp, sqrt(score.mu*(1-score.mu))*score.sigma])
		
	#	check fit in terms of means and sds
	tmp		<- YXr[, list(score.mu=mean(score.p), score.sd=sd(score.p)), by='ntPatient']
	ggplot(mpars, aes(x=as.numeric(ntPatient), y=score.mu)) + 
			geom_point(aes(colour=TYPE)) + geom_line(aes(colour=TYPE, group=TYPE)) +
			geom_point(data=tmp, colour='black') +  
			facet_wrap(~TYPE, ncol=4)
	file	<- paste(indir,'/',gsub('\\.R','_scorePmeansByModel.pdf',infile),sep='')	
	ggsave(file=file, w=15,h=10)
	ggplot(mpars, aes(x=as.numeric(ntPatient), y=score.sd)) + 
			geom_point(aes(colour=TYPE)) + geom_line(aes(colour=TYPE, group=TYPE)) +
			geom_point(data=tmp, colour='black') +  
			facet_wrap(~TYPE, ncol=4)
	file	<- paste(indir,'/',gsub('\\.R','_scorePsdsByModel.pdf',infile),sep='')	
	ggsave(file=file, w=15,h=10)
	#	--> GALNtSLI looks good
	
	#	get samples from fitted distribution
	YXp		<- subset(mpars, TYPE!='OBS')[, {
				if(grepl('Exp',TYPE))
					z<- rEXP(1e4, mu=score.mu)
				if(grepl('GA',TYPE))
					z<- rGA(1e4, mu=score.mu, sigma=score.sigma)
				if(grepl('BE',TYPE))
					z<- rBE(1e4, mu=score.mu, sigma=score.sigma)
				list(score.p=z)
			}, by=c('ntPatientf','TYPE')]
	tmp		<- subset(YXr, select=c(ntPatientf, score.p))
	tmp[, TYPE:= 'OBS']
	YXp		<- rbind(tmp, YXp,use.names=TRUE)
	
	#	check fit in terms of densities
	ggplot(subset(YXp, grepl('BELNtSL|GANtSL|OBS',TYPE)), aes(x=score.p, group=TYPE, colour=TYPE)) +
			scale_x_continuous(limits=c(0,0.1)) +
			scale_colour_brewer(palette='Set1') +
			geom_density() + 
			facet_wrap(~ntPatientf, ncol=5, scales='free_y') 
	file	<- paste(indir,'/',gsub('\\.R','_scorePdensitiesByNtPatient.pdf',infile),sep='')	
	ggsave(file=file, w=15,h=15)	
	#
	file	<- paste(indir,'/',gsub('\\.R','_scorePdensities_mBELNtSL.pdf',infile),sep='')
	pdf(file=file, w=6, h=6)
	plot(mBELNtSL)
	dev.off()
	file	<- paste(indir,'/',gsub('\\.R','_scorePdensities_mGANtSL.pdf',infile),sep='')
	pdf(file=file, w=6, h=6)
	plot(mGANtSL)
	dev.off()
}
######################################################################################
altvtp.plot.150901	<- function(YXc, indir, infile)
{
	ggplot(YXc, aes(x=missexp.pr/DUMMY, fill=stageC)) + geom_histogram(binwidth=0.25) +
			scale_x_continuous(breaks=seq(0,100,5), minor_breaks=seq(0,100,1)) +
			facet_grid(t.AgeC~stageC, scales='free') + 
			theme_bw() 
	file	<- paste(indir,'/',gsub('\\.R','_missexpBytstAgeC150901.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=10)	
	
	ggplot(YXc, aes(x=t.AgeC, y=missexp.pr/DUMMY, fill=stageC)) + geom_boxplot(outlier.shape = NA) +
			scale_y_continuous(breaks=seq(0,100,5), minor_breaks=seq(0,100,1), expand=c(0,0)) +
			coord_cartesian(ylim=c(0,10)) +
			facet_grid(~stageC, space='free', scales='free')  +
			theme_bw() + theme(legend.position='bottom') +
			labs(y='expected missing\n(transmission intervals)') 
	file	<- paste(indir,'/',gsub('\\.R','_missexpBytstAgeC2150901.pdf',infile),sep='')	
	ggsave(file=file, w=13,h=7)	
	
	ggplot(YXc, aes(x=t.AgeC, y=missexp.pr/DUMMY, fill=stageC)) + geom_boxplot(outlier.shape = NA) +
			scale_y_continuous(breaks=seq(0,100,5), minor_breaks=seq(0,100,1), expand=c(0,0)) +
			coord_cartesian(ylim=c(0,10)) +
			facet_grid(stageC~t.period, space='free', scales='free')  +
			theme_bw() + theme(legend.position='bottom') +
			labs(y='expected missing\n(transmission intervals)') 
	file	<- paste(indir,'/',gsub('\\.R','_missexpBytstAgeCPeriod150901.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=10)	
	
	ggplot(YXc, aes(x=t.period, y=missexp.pr/DUMMY, fill=stageC)) + geom_boxplot(outlier.shape = NA) +
			scale_y_continuous(breaks=seq(0,100,5), minor_breaks=seq(0,100,1), expand=c(0,0)) +
			coord_cartesian(ylim=c(0,10)) +
			facet_grid(.~stageC, space='free', scales='free')  +
			theme_bw() + theme(legend.position='bottom') +
			labs(y='expected missing\n(transmission intervals)') 
	file	<- paste(indir,'/',gsub('\\.R','_missexpBystageCPeriod150901.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=6)	
	
	
	ggplot(YXc, aes(x=t.AgeC, y=score.Y, fill=stageC)) + geom_boxplot(outlier.shape = NA) +
			scale_y_continuous(breaks=seq(0,100,5), minor_breaks=seq(0,100,1), expand=c(0,0)) +
			coord_cartesian(ylim=c(0,70)) +
			facet_grid(~stageC, space='free', scales='free')  +
			theme_bw() + theme(legend.position='bottom') +
			labs(y='relative phylogenetic transmission probability\nfor observed transmission intervals\n(scoreY)') 
	file	<- paste(indir,'/',gsub('\\.R','_scoreYBytstAgeC150901.pdf',infile),sep='')	
	ggsave(file=file, w=13,h=7)	
	
	ggplot(YXc, aes(x=t.AgeC, y=100*score.p.nadj, fill=stageC)) + geom_boxplot(outlier.shape = NA) +
			scale_y_continuous(breaks=seq(0,100,5), minor_breaks=seq(0,100,1), expand=c(0,0)) +
			coord_cartesian(ylim=c(0,30)) +
			facet_grid(~stageC, space='free', scales='free')  +
			theme_bw() + theme(legend.position='bottom') +
			labs(y='unadjusted absolute phylogenetic transmission probability\nfor observed transmission intervals\n(%)') 
	file	<- paste(indir,'/',gsub('\\.R','_scorePnadjBytstAgeC150901.pdf',infile),sep='')	
	ggsave(file=file, w=13,h=7)	
	
	ggplot(YXc, aes(x=t.AgeC, y=100*score.p, fill=stageC)) + geom_boxplot(outlier.shape = NA) +
			scale_y_continuous(breaks=seq(0,100,5), minor_breaks=seq(0,100,1), expand=c(0,0)) +
			coord_cartesian(ylim=c(0,10)) +
			facet_grid(~stageC, space='free', scales='free')  +
			theme_bw() + theme(legend.position='bottom') +
			labs(y='absolute phylogenetic transmission probability\nfor observed transmission intervals\n(%)') 
	file	<- paste(indir,'/',gsub('\\.R','_scorePBytstAgeC150901.pdf',infile),sep='')	
	ggsave(file=file, w=13,h=7)		
}
######################################################################################
altvtp.plot.150828	<- function(YXc, indir, infile)
{
	ggplot(YXc, aes(x=missexp, fill=stageC)) + geom_histogram(binwidth=0.25) +
			scale_x_continuous(breaks=seq(0,100,5), minor_breaks=seq(0,100,1)) +
			facet_grid(t.AgeC~stageC, scales='free') + 
			theme_bw() 
	file	<- paste(indir,'/',gsub('\\.R','_missexpBytstAgeC.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=10)	
	
	ggplot(YXc, aes(x=t.AgeC, y=missexp, fill=stageC)) + geom_boxplot(outlier.shape = NA) +
			scale_y_continuous(breaks=seq(0,100,5), minor_breaks=seq(0,100,1), expand=c(0,0)) +
			coord_cartesian(ylim=c(0,10)) +
			facet_grid(~stageC, space='free', scales='free')  +
			theme_bw() + theme(legend.position='bottom') +
			labs(y='expected missing\n(transmission intervals)') 
	file	<- paste(indir,'/',gsub('\\.R','_missexpBytstAgeC2.pdf',infile),sep='')	
	ggsave(file=file, w=13,h=7)	
	
	
	ggplot(YXc, aes(x=t.AgeC, y=score.Y, fill=stageC)) + geom_boxplot(outlier.shape = NA) +
			scale_y_continuous(breaks=seq(0,100,5), minor_breaks=seq(0,100,1), expand=c(0,0)) +
			coord_cartesian(ylim=c(0,70)) +
			facet_grid(~stageC, space='free', scales='free')  +
			theme_bw() + theme(legend.position='bottom') +
			labs(y='relative phylogenetic transmission probability\nfor observed transmission intervals\n(scoreY)') 
	file	<- paste(indir,'/',gsub('\\.R','_scoreYBytstAgeC.pdf',infile),sep='')	
	ggsave(file=file, w=13,h=7)	
	
	ggplot(YXc, aes(x=t.AgeC, y=100*score.p.nadj, fill=stageC)) + geom_boxplot(outlier.shape = NA) +
			scale_y_continuous(breaks=seq(0,100,5), minor_breaks=seq(0,100,1), expand=c(0,0)) +
			coord_cartesian(ylim=c(0,30)) +
			facet_grid(~stageC, space='free', scales='free')  +
			theme_bw() + theme(legend.position='bottom') +
			labs(y='unadjusted absolute phylogenetic transmission probability\nfor observed transmission intervals\n(%)') 
	file	<- paste(indir,'/',gsub('\\.R','_scorePnadjBytstAgeC.pdf',infile),sep='')	
	ggsave(file=file, w=13,h=7)	
	
	ggplot(YXc, aes(x=t.AgeC, y=100*score.p, fill=stageC)) + geom_boxplot(outlier.shape = NA) +
			scale_y_continuous(breaks=seq(0,100,5), minor_breaks=seq(0,100,1), expand=c(0,0)) +
			coord_cartesian(ylim=c(0,10)) +
			facet_grid(~stageC, space='free', scales='free')  +
			theme_bw() + theme(legend.position='bottom') +
			labs(y='absolute phylogenetic transmission probability\nfor observed transmission intervals\n(%)') 
	file	<- paste(indir,'/',gsub('\\.R','_scorePBytstAgeC.pdf',infile),sep='')	
	ggsave(file=file, w=13,h=7)		
}
######################################################################################
altvtp.confounded.ntransmitters<- function(YXc, indir, infile, suffix='')
{
	tmp		<- YXc[, list(score.p.me=mean(score.p)), by=c('ntPatient')]
	ggplot( tmp, aes(x=ntPatient, y=score.p.me)) + geom_bar(stat='identity')
	file	<- paste(indir,'/',gsub('\\.R',paste('_scorePmeansByNtPatient_raw',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=5,h=5)		
	tmp		<- YXc[, list(score.p.m=mean(score.p), stageC=stageC[1], t.AgeC=t.AgeC[1]), by=c('ntPatient','stageC')]
	ggplot( tmp, aes(x=ntPatient, y=score.p.m)) + geom_bar(stat='identity') + facet_grid(stageC~.)
	file	<- paste(indir,'/',gsub('\\.R',paste('_scorePmeansBystageCNtPatient_raw',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=5,h=7)
	tmp		<- YXc[, list(score.p.m=mean(score.p), stageC=stageC[1], t.AgeC=t.AgeC[1]), by=c('ntPatient','t.AgeC')]
	ggplot( tmp, aes(x=ntPatient, y=score.p.m)) + geom_bar(stat='identity') + facet_grid(t.AgeC~.)
	file	<- paste(indir,'/',gsub('\\.R',paste('_scorePmeansBytAgeCNtPatient_raw',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=5,h=7)	
	tmp		<- YXc[, list(score.p.m=mean(score.p), stageC=stageC[1], t.AgeC=t.AgeC[1]), by=c('ntPatient','t.stAgeC')]
	ggplot( tmp, aes(x=ntPatient, y=score.p.m)) + geom_bar(stat='identity') + facet_grid(t.AgeC~stageC)
	file	<- paste(indir,'/',gsub('\\.R',paste('_scorePmeansBytstAgeCNtPatient_raw',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=20,h=10)
}
######################################################################################
adjust.dev.Age_253045.Stage_UA_UC_D_TS_TO_F<- function()
{
	require(survival)
	YX
	X.tables
	t.recent.endctime
	
	indir		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C2V100bInfT7STRAT_m5C.R'
	
	YXo			<- copy(YX)	#before strat
	#YXo2		<-	stratificationmodel.Age_253045.Stage_UA_UC_D_TS_TO_F(YXo, lRNA.supp=2)
	YXo2		<-	stratificationmodel.Age_253045.Stage_UAE_UAC_UC_D_TS_TO_F(YXo, lRNA.supp=2)
	risk.col	<- 't.stAgeC'	
	#
	#	focus on >2003 and before ART start
	#
	#YX			<- subset(YX, t>2004.5 & !stageC%in%c('L','T'))
	YX			<- copy(YXo2)
	YX			<- subset(YX, AnyPos_T1>2005)
	tmp			<- YX[, which(stageC=='L')]
	set(YX, tmp, 't.AgeC', '(-1,100]')
	set(YX, tmp, 't.stAgeC', 'L_(-1,100]')
	tmp			<- YX[, which(stageC=='TS')]
	set(YX, tmp, 't.AgeC', '(-1,100]')
	set(YX, tmp, 't.stAgeC', 'TS_(-1,100]')
	tmp			<- YX[, which(stageC=='TO')]
	set(YX, tmp, 't.AgeC', '(-1,100]')
	set(YX, tmp, 't.stAgeC', 'TO_(-1,100]')		
	#tmp		<- YXc[, which(stageC=='L' | stageC=='T')]
	#set(YXc, tmp, 't.AgeC', '(-1,100]')
	#set(YXc, tmp, 'stageC', 'LT')
	#set(YXc, tmp, 't.stAgeC', 'LT_(-1,100]')	
	set(YX, NULL, 'stageC', YX[, factor(as.character(stageC))])
	set(YX, NULL, 't.stAgeC', YX[, factor(as.character(t.stAgeC))])
	set(YX, NULL, 't.stAgeC.prd', YX[, factor(as.character(t.stAgeC.prd))])
	
	#	get clustering, sampling and non-censoring probabilities 
	#	for each prob transmitter at the transmission interval 
	sm			<- X.tables$sm
	YXs			<- merge(YX, sm, by='t.Patient')
	#	p.nc	(non-censoring prob)
	cp			<- copy(X.tables$cm$cens.p)		#model built after subsampling, so cannot apply straight-away
	set(cp, NULL, c('BS','p.nc'), NULL)
	cm			<- X.tables$cm$cens.m	
	tp.df		<- subset(YXs, grepl('^U', YXs[[risk.col]]))[, list(tm=mean(t)), by=c('t.Patient','Patient')]
	tp.df		<- censoring.get.dataset(tp.df, df.all.allmsm)
	setnames(YXs, c('t.Patient','Patient'), c('Patient','r.Patient'))
	#	YXc so far only contains transmitters with at least one undiagnosed interval
	YXc			<- merge(YXs, subset(tp.df, select=c(Patient, r.Patient, tm, AACD4C)), by=c('Patient','r.Patient'))
	set(YXc, NULL, 'tm', YXc[, tm-t.recent.endctime])	
	YXc[, p.nc:=predict(cm, data=cp, newdata=subset(YXc, select=c(tm, AACD4C)), type='response')]
	#	p.nc is allocated to transmitters with at least one undiagnosed interval
	#	set p.nc to 1 for all intervals after diagnosis
	set(YXc, YXc[, which(!is.na(p.nc) & !grepl('^U', stageC))], 'p.nc', 1)
	#ggplot(YXc, aes(x=tm, y=p.nc, group=AACD4C)) + geom_line() + facet_wrap(~AACD4C, ncol=4)
	#ggplot(YXc, aes(x=tm)) + geom_histogram()	
	setnames(YXs, c('Patient','r.Patient'), c('t.Patient','Patient'))
	setnames(YXc, c('Patient','r.Patient'), c('t.Patient','Patient'))
	YXc			<- merge(YXs, subset(YXc, select=c(t.Patient, Patient, t, p.nc)), by=c('t.Patient','Patient','t'), all.x=TRUE)
	set(YXc, YXc[, which(is.na(p.nc))], 'p.nc', 1)
	#	get p.clu and p.pt
	risk.table	<- X.tables$st$risk.table
	tmp			<- clustering.getrisks(risk.table, indir=indir, infile=infile, rskf.name='t.stAgeC', rskf.baseline='D_(30,45]')
	tpt			<- tmp$tpt
	tclu		<- tmp$tclu	
	tmp			<- subset(tclu, t.period=='Overall', select=c(factor2, p.clu))
	setnames(tmp, 'factor2', 't.stAgeC')	
	set(tmp, NULL, 't.stAgeC', tmp[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)], labels=YXc[, levels(t.stAgeC)])])
	YXc			<- merge(YXc, tmp, by=c('t.stAgeC'), all.x=TRUE)
	tmp			<- subset(tpt, t.period=='Overall', select=c(factor2, p.pt))
	setnames(tmp, 'factor2', 't.stAgeC')		
	set(tmp, NULL, 't.stAgeC', tmp[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)], labels=YXc[, levels(t.stAgeC)])])
	YXc			<- merge(YXc, tmp, by=c('t.stAgeC'), all.x=TRUE)	
	#	get w per interval
	cat(paste('\nsetting likelihood to likelihood of pair / number transmission intervals'))
	set(YXc, NULL, 'score.Y.raw', YXc[, score.Y])
	set(YXc, NULL, 'score.Y', YXc[, score.Y*w.tn])				
	YXc		<- merge(YXc, YXc[, list(t=t, t.Patient=t.Patient, ntPatient= length(unique(t.Patient))), by=c('Patient')], by=c('Patient','t.Patient','t'))	
	YXc[, p.seqnc:= p.seq*p.nc]
	YXc[, me:= 1/(p.seq*p.nc)-1]	
	#	calculate absolute transmission probabilities
	#	need first: predicted score by stage and time t
	#	need second: predicted missing intervals per recipient, by stage and time t
	#	then calculate: exp missing score.Y per recipient by time t= sum_stage missexp.pr(stage,t) * scoreY.pr(stage,t)
	#tmp		<- rltvtp.model.time.150901(YXc)
	tmp		<- rltvtp.model.time.150916(YXc)
	rltvp	<- tmp$predict
	rltva	<- tmp$predict.args
	YXc		<- subset(YXc, AnyPos_T1<2010.8)
	tmp		<- miss.model.150904b(YXc)	
	md		<- tmp$miss.d
	mm		<- tmp$miss.m
	#YXm		<- miss.smoothing.150909(YXc)
	#YXm		<- miss.smoothing.150914(YXc, mm, md, rltvm, rltvd, indir=NA, infile=NA, suffix='150914')
	YXm		<- miss.smoothing.150916(YXc, mm, md, rltvp, rltva, indir=NA, infile=NA, suffix='150916')
	YXm		<- merge(YXm, unique(subset(YXc, select=c(Patient, AnyPos_T1, t.period))), by='Patient')
	if(0)
	{
		#	get expected missing term for denominator 'scoreYm.per.rec'
		tmp		<- YXm[, list(miss.per.rec= sum(miss.prsm), scoreYm.per.rec= sum(miss.prsm*scoreY.pr)), by='Patient']
		YXc		<- merge( YXc, tmp, by='Patient' )
		#	then calculate p = w / ( sum(w) + scoreYm.per.rec)
		YXc		<- merge(YXc, YXc[, list(score.p= score.Y/(sum(score.Y)+scoreYm.per.rec), score.p.nadj= score.Y/sum(score.Y), t.Patient=t.Patient, t=t), by='Patient'], by=c('Patient','t.Patient','t'))
		#ggplot(YXc, aes(x=AnyPos_T1, y=miss.per.rec)) + geom_point()		
	}
	
	#prop.consistency(YXc, YXm)
	
	#	
	YXm		<- YXm[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr), miss.per.recst= sum(miss.prsm)), by=c('t.period','t.stAgeC','Patient')]
	#	merge with observed scores by stage for each transmitter	
	#	DEPENDS ON STRATIFICATION - explore first on t.stAgeC
	tmp		<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.stAgeC')]
	YXpr	<- merge(tmp, YXm, by=c('Patient','t.stAgeC'), all=1)
	set(YXpr, YXpr[, which(is.na(obs.per.recst))], c('obs.per.recst','scoreYo.per.recst'), 0)
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.stAgeC=t.stAgeC, pnadj.per.rec=scoreYo.per.recst/sum(scoreYo.per.recst), p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.stAgeC'))
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))
	YXpr	<- merge(YXpr, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	
	tmp		<- YXpr[, list(p.per.rec=sum(p.per.rec), po.per.rec=sum(p.per.rec*(1-pm.per.rec)), pm.per.rec=sum(p.per.rec*pm.per.rec)), by=c('Patient','stageC','t.period')]
	tmp		<- tmp[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	
	
	
	#	get rolling mean and plot
	setkey(YXpr, t.stAgeC, tmid)	
	tmp		<- YXpr[, {
				z	<- unique(tmid)
				list(tmid=z, prm.per.rec= sapply(seq_along(z), function(i)	mean(p.per.rec[ which(abs(z[i]-tmid)<=1) ]))	)	
			}, by='t.stAgeC']
	YXpr	<- merge(YXpr, tmp, by=c('tmid','t.stAgeC'))
	ggplot(YXpr, aes(x=tmid)) + geom_point(aes(y=p.per.rec, colour=t.AgeC)) + geom_line(aes(y=prm.per.rec), colour='black') + 
			facet_grid(t.AgeC~stageC) + theme_bw() 
	
	
	#	plot proportion of obs+missing across time
	YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009,Inf), labels=c('<2007','2008','2009','>2009'))]
	YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2009,Inf), labels=c('<2007','2008/09','>2009'))]
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.stAgeC','tmidC')]
	tmp		<- merge(tmp, unique(subset(YXpr, select=c(t.stAgeC, t.AgeC, stageC))), by='t.stAgeC')
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=tmidC, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(t.AgeC~stageC) + theme_bw() + theme(legend.position='bottom')	
}
######################################################################################
adjust.dev.Age_253045.Stage_UAC_UAE_UC_D_TS_TO_F<- function()
{
	require(survival)
	YX
	X.tables
	t.recent.endctime
	
	indir		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C2V100bInfT7STRAT_m5B.R'

	YXo			<- copy(YX)	#before strat
	#YXo2		<-	stratificationmodel.Age_253045.Stage_UA_UC_D_TS_TO_F(YXo, lRNA.supp=2)
	YXo2		<-	stratificationmodel.Age_253045.Stage_UAE_UAC_UC_D_TS_TO_F(YXo, lRNA.supp=2)
	risk.col	<- 't.stAgeC'	
	#
	#	focus on >2003 and before ART start
	#
	#YX			<- subset(YX, t>2004.5 & !stageC%in%c('L','T'))
	YX			<- copy(YXo2)
	YX			<- subset(YX, AnyPos_T1>2005)
	
	#	get clustering, sampling and non-censoring probabilities 
	#	for each prob transmitter at the transmission interval 
	sm			<- X.tables$sm
	YXs			<- merge(YX, sm, by='t.Patient')
	#	p.nc	(non-censoring prob)
	cp			<- copy(X.tables$cm$cens.p)		#model built after subsampling, so cannot apply straight-away
	set(cp, NULL, c('BS','p.nc'), NULL)
	cm			<- X.tables$cm$cens.m	
	tp.df		<- subset(YXs, grepl('^U', YXs[[risk.col]]))[, list(tm=mean(t)), by=c('t.Patient','Patient')]
	tp.df		<- censoring.get.dataset(tp.df, df.all.allmsm)
	setnames(YXs, c('t.Patient','Patient'), c('Patient','r.Patient'))
	#	YXc so far only contains transmitters with at least one undiagnosed interval
	YXc			<- merge(YXs, subset(tp.df, select=c(Patient, r.Patient, tm, AACD4C)), by=c('Patient','r.Patient'))
	set(YXc, NULL, 'tm', YXc[, tm-t.recent.endctime])	
	YXc[, p.nc:=predict(cm, data=cp, newdata=subset(YXc, select=c(tm, AACD4C)), type='response')]
	#	p.nc is allocated to transmitters with at least one undiagnosed interval
	#	set p.nc to 1 for all intervals after diagnosis
	set(YXc, YXc[, which(!is.na(p.nc) & !grepl('^U', stageC))], 'p.nc', 1)
	#ggplot(YXc, aes(x=tm, y=p.nc, group=AACD4C)) + geom_line() + facet_wrap(~AACD4C, ncol=4)
	#ggplot(YXc, aes(x=tm)) + geom_histogram()	
	setnames(YXs, c('Patient','r.Patient'), c('t.Patient','Patient'))
	setnames(YXc, c('Patient','r.Patient'), c('t.Patient','Patient'))
	YXc			<- merge(YXs, subset(YXc, select=c(t.Patient, Patient, t, p.nc)), by=c('t.Patient','Patient','t'), all.x=TRUE)
	set(YXc, YXc[, which(is.na(p.nc))], 'p.nc', 1)
	#	get p.clu and p.pt
	risk.table	<- X.tables$st$risk.table
	tmp			<- clustering.getrisks(risk.table, indir=indir, infile=infile, rskf.name='t.stAgeC', rskf.baseline='D_(30,45]')
	tpt			<- tmp$tpt
	tclu		<- tmp$tclu	
	tmp			<- subset(tclu, t.period=='Overall', select=c(factor2, p.clu))
	setnames(tmp, 'factor2', 't.stAgeC')	
	set(tmp, NULL, 't.stAgeC', tmp[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)], labels=YXc[, levels(t.stAgeC)])])
	YXc			<- merge(YXc, tmp, by=c('t.stAgeC'), all.x=TRUE)
	tmp			<- subset(tpt, t.period=='Overall', select=c(factor2, p.pt))
	setnames(tmp, 'factor2', 't.stAgeC')		
	set(tmp, NULL, 't.stAgeC', tmp[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)], labels=YXc[, levels(t.stAgeC)])])
	YXc			<- merge(YXc, tmp, by=c('t.stAgeC'), all.x=TRUE)	
	#	get w per interval
	cat(paste('\nsetting likelihood to likelihood of pair / number transmission intervals'))
	set(YXc, NULL, 'score.Y.raw', YXc[, score.Y])
	set(YXc, NULL, 'score.Y', YXc[, score.Y*w.tn])				
	YXc		<- merge(YXc, YXc[, list(t=t, t.Patient=t.Patient, ntPatient= length(unique(t.Patient))), by=c('Patient')], by=c('Patient','t.Patient','t'))	
	YXc[, p.seqnc:= p.seq*p.nc]
	YXc[, me:= 1/(p.seq*p.nc)-1]	
	#	calculate absolute transmission probabilities
	#	need first: predicted score by stage and time t
	#	need second: predicted missing intervals per recipient, by stage and time t
	#	then calculate: exp missing score.Y per recipient by time t= sum_stage missexp.pr(stage,t) * scoreY.pr(stage,t)
	#tmp		<- rltvtp.model.time.150901(YXc)
	#tmp		<- rltvtp.model.time.150917(YXc, indir=indir, infile=infile, suffix='150918')
	tmp		<- rltvtp.model.time.150918(YXc)
	rltvp	<- tmp$prd
	rltvd	<- tmp$fit
	rltvm	<- tmp$mo
	YXc		<- subset(YXc, AnyPos_T1<2010.8)
	tmp		<- miss.model.150917(YXc)		
	md		<- tmp$miss.d
	mm		<- tmp$miss.m
	#YXm		<- miss.smoothing.150909(YXc)
	#YXm		<- miss.smoothing.150914(YXc, mm, md, rltvm, rltvd, indir=NA, infile=NA, suffix='150914')
	#YXm		<- miss.smoothing.150916(YXc, mm, md, rltvp, rltva, indir=indir, infile=infile, suffix='150916')
	#YXm		<- miss.smoothing.150917(YXc, mm, md, rltvp, rltva, indir=indir, infile=infile, suffix='150917', use.obs.pseqnc=TRUE)
	YXm		<- miss.smoothing.150918(YXc, mm, md, rltvp, rltvd, rltvm, indir=indir, infile=infile, suffix='150918', use.obs.pseqnc=TRUE)
	YXm		<- merge(YXm, unique(subset(YXc, select=c(Patient, AnyPos_T1, t.period))), by='Patient')
	if(0)
	{
		#	get expected missing term for denominator 'scoreYm.per.rec'
		tmp		<- YXm[, list(miss.per.rec= sum(miss.prsm), scoreYm.per.rec= sum(miss.prsm*scoreY.pr)), by='Patient']
		YXc		<- merge( YXc, tmp, by='Patient' )
		#	then calculate p = w / ( sum(w) + scoreYm.per.rec)
		YXc		<- merge(YXc, YXc[, list(score.p= score.Y/(sum(score.Y)+scoreYm.per.rec), score.p.nadj= score.Y/sum(score.Y), t.Patient=t.Patient, t=t), by='Patient'], by=c('Patient','t.Patient','t'))
		#ggplot(YXc, aes(x=AnyPos_T1, y=miss.per.rec)) + geom_point()		
	}
	
	#prop.consistency(YXc, YXm)

	tmp2	<- YXm[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.stAgeC','Patient')]	
	tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','t.stAgeC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','t.stAgeC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.stAgeC=t.stAgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.stAgeC'))	
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))
	YXpr	<- merge(YXpr, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	
	#	consistency	--> aggregating first by t.stAgeC leads to slightly diff results
	tmp		<- YXpr[, list(p.per.rec=sum(p.per.rec), po.per.rec=sum(p.per.rec*(1-pm.per.rec)), pm.per.rec=sum(p.per.rec*pm.per.rec)), by=c('Patient','stageC','t.period')]	
	tmp		<- tmp[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(po.per.rec), pm.per.tc=mean(pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodStageC',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=7)
	#	consistency: combine UAE and UAC		
	tmp		<- copy(YXpr)
	set(tmp, NULL, 'stageC', tmp[, gsub('UAE|UAC','UA', stageC)])
	tmp		<- tmp[, list(p.per.rec=sum(p.per.rec), po.per.rec=sum(p.per.rec*(1-pm.per.rec)), pm.per.rec=sum(p.per.rec*pm.per.rec)), by=c('Patient','stageC','t.period')]
	tmp		<- tmp[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(po.per.rec), pm.per.tc=mean(pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodStageC_UA',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=7)
	
	#	by t.stAgeC	
	tmp		<- YXpr[, list(stageC=stageC[1], t.AgeC=t.AgeC[1], p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.stAgeC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(t.AgeC~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodtstAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=10)

	#	by t.AgeC
	tmp2	<- YXm[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.AgeC','Patient')]	
	tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','t.AgeC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','t.AgeC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.AgeC=t.AgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.AgeC'))	
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~t.AgeC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodtAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=7)
	YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009,Inf), labels=c('<2007','2007','2008','>=2009'))]
	#YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009, 2010, Inf), labels=c('<2007','2007','2008','2009','2010'))]
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','tmidC')]	
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggp		<- ggplot(tmp, aes(x=t.AgeC, y=100*value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') +
			scale_y_continuous(expand=c(0,0), limits=c(0,55)) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			scale_fill_brewer(palette='Set1', guide=FALSE) +
			facet_grid(~tmidC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) +
			labs(x='\nage of probable transmitter', y='%', alpha='')
	z <- ggplot_gtable(ggplot_build(ggp))
	z <- gtable_add_rows(z, z$heights[[3]], 2)
	z <- gtable_add_grob(z, 
			list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
					textGrob("Transmissions by calendar year", gp = gpar(col = 'black'))), 3, 4, 3, 10, name = paste(runif(2)))
	z <- gtable_add_rows(z, unit(1/8, "line"), 3)
	grid.newpage()
	pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByTmid',suffix,'.pdf',sep=''),infile), sep=''), w=12, h=5)
	grid.draw(z)
	dev.off()
	
	#
	#	by p.AgeC
	#
	YXmp	<- merge(YXm, unique(subset(YXc, select=c(Patient, AgeC))), by='Patient')
	tmp2	<- YXmp[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.AgeC','AgeC','Patient')]
	tmp2[, p.AgeC:= paste(as.character(t.AgeC),'->',as.character(AgeC),sep='')]
	tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','p.AgeC','t.AgeC', 'AgeC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','AgeC','t.AgeC','p.AgeC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.AgeC=t.AgeC, AgeC=AgeC, p.AgeC=p.AgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.AgeC','AgeC','p.AgeC'))	
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))	
	#	p.AgeC proportions
	require(gtable)
	tmp		<- YXpr[, list(p.per.tc=100*mean(p.per.rec), po.per.tc=100*mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=100*mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))	
	ggp		<- ggplot(tmp, aes(x=t.AgeC, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~AgeC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) + 
			scale_y_continuous(expand=c(0,0), limits=c(0,55)) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			scale_fill_brewer(palette='Set1', guide=FALSE) +
			labs(x='\nage of probable transmitter', y='%', alpha='')	
	z <- ggplot_gtable(ggplot_build(ggp))
	z <- gtable_add_rows(z, z$heights[[3]], 2)
	z <- gtable_add_grob(z, 
			list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
					textGrob("Transmissions by age of recipient", gp = gpar(col = 'black'))), 3, 4, 3, 10, name = paste(runif(2)))
	z <- gtable_add_rows(z, unit(1/8, "line"), 3)
	grid.newpage()
	pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=12, h=5)
	grid.draw(z)
	dev.off()
	
	#	p.AgeC total
	tmp		<- YXpr[, list(n.per.tc=sum(p.per.rec), no.per.tc=sum(p.per.rec*(1-pm.per.rec)), nm.per.tc=sum(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC')]
	tmp[, p.per.tc:= 100*n.per.tc/sum(n.per.tc)]
	tmp[, text:= as.character(round(100*n.per.tc/sum(n.per.tc),d=1))]
	set(tmp, tmp[, which(p.per.tc<3)], 'text', '')
	tmp[, po.per.tc:= 100*no.per.tc/sum(n.per.tc)]
	tmp[, pm.per.tc:= 100*nm.per.tc/sum(n.per.tc)]	
	setnames(tmp, c('po.per.tc','p.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('expected missing','observed'), labels=c('expected missing','observed'))])
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('observed','expected missing'), labels=c('observed','expected missing'))])
	ggplot(tmp, aes(x=AgeC, y=t.AgeC)) + 			
			geom_point(aes(size=value, alpha=variable, fill=variable), pch=21) +
			geom_text(aes(label=text), colour='black', size=4, hjust=0.5, vjust=0.5) +
			scale_size_continuous(range=c(4,50), guide=FALSE) + 
			scale_fill_manual(values=c('expected missing'='#E41A1C', 'observed'='#E41A1C')) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			theme_bw() + theme(legend.position='bottom') + labs(x='age of recipient', y='age of probable transmitter', fill='transmissions\n(%)', alpha='transmissions\n(%)')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_AgeCtAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=5, h=5)
	
	
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(AgeC~t.AgeC) + theme_bw() + theme(legend.position='bottom')

	
	tmp		<- merge(tmp, unique(subset(YXpr, select=c(t.stAgeC, t.AgeC, stageC))), by='t.stAgeC')
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=tmidC, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(t.AgeC~stageC) + theme_bw() + theme(legend.position='bottom')
	
	
	#	get rolling mean and plot
	setkey(YXpr, t.stAgeC, tmid)	
	tmp		<- YXpr[, {
				z	<- unique(tmid)
				list(tmid=z, prm.per.rec= sapply(seq_along(z), function(i)	mean(p.per.rec[ which(abs(z[i]-tmid)<=1) ]))	)	
			}, by='t.stAgeC']
	YXpr	<- merge(YXpr, tmp, by=c('tmid','t.stAgeC'))
	ggplot(YXpr, aes(x=tmid)) + geom_point(aes(y=p.per.rec, colour=t.AgeC)) + geom_line(aes(y=prm.per.rec), colour='black') + 
			facet_grid(t.AgeC~stageC) + theme_bw() 
	
	
	#	plot proportion of obs+missing across time
	YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009,Inf), labels=c('<2007','2008','2009','>2009'))]
	YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2009,Inf), labels=c('<2007','2008/09','>2009'))]
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.stAgeC','tmidC')]
	tmp		<- merge(tmp, unique(subset(YXpr, select=c(t.stAgeC, t.AgeC, stageC))), by='t.stAgeC')
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=tmidC, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(t.AgeC~stageC) + theme_bw() + theme(legend.position='bottom')	
}
######################################################################################
adjust.dev.Gen_55657585.Stage_UAC_UAE_UC_D_TS_TO_F<- function()
{
	require(survival)
	require(gtable)
	YX
	X.tables
	t.recent.endctime
	suffix		<- '150922'
	indir		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C2V100bInfT7STRAT_m5F.R'
	
	YXo			<- copy(YX)	#before strat
	#YXo2		<-	stratificationmodel.Age_253045.Stage_UA_UC_D_TS_TO_F(YXo, lRNA.supp=2)
	YXo2		<-	stratificationmodel.Gen_55657585.Stage_UAE_UAC_UC_D_TS_TO_F(YXo, lRNA.supp=2)
	risk.col	<- 't.stAgeC'	
	#
	#	focus on >2003 and before ART start
	#
	#YX			<- subset(YX, t>2004.5 & !stageC%in%c('L','T'))
	YX			<- copy(YXo2)
	YX			<- subset(YX, AnyPos_T1>2005)
	
	#	get clustering, sampling and non-censoring probabilities 
	#	for each prob transmitter at the transmission interval 
	sm			<- X.tables$sm
	YXs			<- merge(YX, sm, by='t.Patient')
	#	p.nc	(non-censoring prob)
	cp			<- copy(X.tables$cm$cens.p)		#model built after subsampling, so cannot apply straight-away
	set(cp, NULL, c('BS','p.nc'), NULL)
	cm			<- X.tables$cm$cens.m	
	tp.df		<- subset(YXs, grepl('^U', YXs[[risk.col]]))[, list(tm=mean(t)), by=c('t.Patient','Patient')]
	tp.df		<- censoring.get.dataset(tp.df, df.all.allmsm)
	setnames(YXs, c('t.Patient','Patient'), c('Patient','r.Patient'))
	#	YXc so far only contains transmitters with at least one undiagnosed interval
	YXc			<- merge(YXs, subset(tp.df, select=c(Patient, r.Patient, tm, AACD4C)), by=c('Patient','r.Patient'))
	set(YXc, NULL, 'tm', YXc[, tm-t.recent.endctime])	
	YXc[, p.nc:=predict(cm, data=cp, newdata=subset(YXc, select=c(tm, AACD4C)), type='response')]
	#	p.nc is allocated to transmitters with at least one undiagnosed interval
	#	set p.nc to 1 for all intervals after diagnosis
	set(YXc, YXc[, which(!is.na(p.nc) & !grepl('^U', stageC))], 'p.nc', 1)
	#ggplot(YXc, aes(x=tm, y=p.nc, group=AACD4C)) + geom_line() + facet_wrap(~AACD4C, ncol=4)
	#ggplot(YXc, aes(x=tm)) + geom_histogram()	
	setnames(YXs, c('Patient','r.Patient'), c('t.Patient','Patient'))
	setnames(YXc, c('Patient','r.Patient'), c('t.Patient','Patient'))
	YXc			<- merge(YXs, subset(YXc, select=c(t.Patient, Patient, t, p.nc)), by=c('t.Patient','Patient','t'), all.x=TRUE)
	set(YXc, YXc[, which(is.na(p.nc))], 'p.nc', 1)
	#	get p.clu and p.pt
	risk.table	<- X.tables$st$risk.table
	tmp			<- clustering.getrisks(risk.table, indir=indir, infile=infile, rskf.name='t.stAgeC', rskf.baseline='D_(30,45]')
	tpt			<- tmp$tpt
	tclu		<- tmp$tclu	
	tmp			<- subset(tclu, t.period=='Overall', select=c(factor2, p.clu))
	setnames(tmp, 'factor2', 't.stAgeC')	
	set(tmp, NULL, 't.stAgeC', tmp[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)], labels=YXc[, levels(t.stAgeC)])])
	YXc			<- merge(YXc, tmp, by=c('t.stAgeC'), all.x=TRUE)
	tmp			<- subset(tpt, t.period=='Overall', select=c(factor2, p.pt))
	setnames(tmp, 'factor2', 't.stAgeC')		
	set(tmp, NULL, 't.stAgeC', tmp[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)], labels=YXc[, levels(t.stAgeC)])])
	YXc			<- merge(YXc, tmp, by=c('t.stAgeC'), all.x=TRUE)	
	#	get w per interval
	cat(paste('\nsetting likelihood to likelihood of pair / number transmission intervals'))
	set(YXc, NULL, 'score.Y.raw', YXc[, score.Y])
	set(YXc, NULL, 'score.Y', YXc[, score.Y*w.tn])				
	YXc		<- merge(YXc, YXc[, list(t=t, t.Patient=t.Patient, ntPatient= length(unique(t.Patient))), by=c('Patient')], by=c('Patient','t.Patient','t'))	
	YXc[, p.seqnc:= p.seq*p.nc]
	YXc[, me:= 1/(p.seq*p.nc)-1]	
	#	calculate absolute transmission probabilities
	#	need first: predicted score by stage and time t
	#	need second: predicted missing intervals per recipient, by stage and time t
	#	then calculate: exp missing score.Y per recipient by time t= sum_stage missexp.pr(stage,t) * scoreY.pr(stage,t)
	#tmp		<- rltvtp.model.time.150901(YXc)
	#tmp		<- rltvtp.model.time.150917(YXc, indir=indir, infile=infile, suffix='150918')
	#tmp		<- rltvtp.model.time.150918(YXc)
	tmp		<- rltvtp.model.time.150923Gen(YXc)
	rltvp	<- tmp$prd
	rltvd	<- tmp$fit
	rltvm	<- tmp$mo
	YXc		<- subset(YXc, AnyPos_T1<2010.8)	
	tmp		<- miss.model.150923Gen(YXc)
	mp		<- tmp$prd
	mma		<- tmp$moa
	mmb		<- tmp$mob
	md		<- tmp$fit
	YXm		<- miss.smoothing.150922(YXc, mp, mma, mmb, md, rltvp, rltvd, rltvm, indir=indir, infile=infile, suffix='150923', use.obs.pseqnc=TRUE)
	YXm		<- merge(YXm, unique(subset(YXc, select=c(Patient, AnyPos_T1, t.period))), by='Patient')
	if(0)
	{
		#	get expected missing term for denominator 'scoreYm.per.rec'
		tmp		<- YXm[, list(miss.per.rec= sum(miss.prsm), scoreYm.per.rec= sum(miss.prsm*scoreY.pr)), by='Patient']
		YXc		<- merge( YXc, tmp, by='Patient' )
		#	then calculate p = w / ( sum(w) + scoreYm.per.rec)
		YXc		<- merge(YXc, YXc[, list(score.p= score.Y/(sum(score.Y)+scoreYm.per.rec), score.p.nadj= score.Y/sum(score.Y), t.Patient=t.Patient, t=t), by='Patient'], by=c('Patient','t.Patient','t'))
		#ggplot(YXc, aes(x=AnyPos_T1, y=miss.per.rec)) + geom_point()		
	}
	#
	#	consistency checks
	#
	#prop.consistency(YXc, YXm)	
	tmp2	<- YXm[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.stAgeC','Patient')]	
	tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','t.stAgeC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','t.stAgeC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.stAgeC=t.stAgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.stAgeC'))	
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))
	YXpr	<- merge(YXpr, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')	
	#	consistency	--> aggregating first by t.stAgeC leads to slightly diff results
	tmp		<- YXpr[, list(p.per.rec=sum(p.per.rec), po.per.rec=sum(p.per.rec*(1-pm.per.rec)), pm.per.rec=sum(p.per.rec*pm.per.rec)), by=c('Patient','stageC','t.period')]	
	tmp		<- tmp[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(po.per.rec), pm.per.tc=mean(pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodStageC',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=7)
	#	consistency: combine UAE and UAC		
	tmp		<- copy(YXpr)
	set(tmp, NULL, 'stageC', tmp[, gsub('UAE|UAC','UA', stageC)])
	tmp		<- tmp[, list(p.per.rec=sum(p.per.rec), po.per.rec=sum(p.per.rec*(1-pm.per.rec)), pm.per.rec=sum(p.per.rec*pm.per.rec)), by=c('Patient','stageC','t.period')]
	tmp		<- tmp[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(po.per.rec), pm.per.tc=mean(pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodStageC_UA',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=7)
	
	#	by t.stAgeC	
	tmp		<- YXpr[, list(stageC=stageC[1], t.AgeC=t.AgeC[1], p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.stAgeC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(t.AgeC~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodtstAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=10)
	
	#	transmissions by calendar year / missingness; and DOB of transmitter
	#	
	#	by t.AgeC
	tmp2	<- YXm[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.AgeC','Patient')]	
	tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','t.AgeC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','t.AgeC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.AgeC=t.AgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.AgeC'))	
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') +
			scale_fill_hue(h=c(200, 360), l=50, c=100, guide=FALSE) +			
			facet_grid(~t.AgeC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodtAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=7)
	YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009,Inf), labels=c('<2007','2007','2008','>=2009'))]
	#YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009, 2010, Inf), labels=c('<2007','2007','2008','2009','2010'))]
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','tmidC')]	
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggp		<- ggplot(tmp, aes(x=t.AgeC, y=100*value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') +
			scale_y_continuous(expand=c(0,0), limits=c(0,55)) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			scale_fill_hue(h=c(200, 360), l=50, c=100, guide=FALSE) +
			facet_grid(~tmidC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) +
			labs(x='\nDOB of probable transmitter', y='%', alpha='')
	z <- ggplot_gtable(ggplot_build(ggp))
	z <- gtable_add_rows(z, z$heights[[3]], 2)
	z <- gtable_add_grob(z, 
			list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
					textGrob("Transmissions by calendar year", gp = gpar(col = 'black'))), 3, 4, 3, 10, name = paste(runif(2)))
	z <- gtable_add_rows(z, unit(1/8, "line"), 3)
	grid.newpage()
	pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByTmid',suffix,'.pdf',sep=''),infile), sep=''), w=16, h=5)
	grid.draw(z)
	dev.off()
	
	#	transmissions by calendar year / diagnosis status; and DOB of transmitter
	#	by t.AgeC2
	YXm[, t.AgeC2:= YXm[, gsub('TO|TS|L','D',gsub('U[^_]+','U',t.stAgeC))]]
	tmp2	<- YXm[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.AgeC2','t.AgeC','Patient')]
	YXc[, t.AgeC2:= YXc[, gsub('TO|TS|L','D',gsub('U[^_]+','U',t.stAgeC))]]
	tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','t.AgeC2','t.AgeC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','t.AgeC2','t.AgeC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.AgeC2=t.AgeC2, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.AgeC2'))	
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))
	YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009,Inf), labels=c('<2007','2007','2008','>=2009'))]
	YXpr[, stageC:= factor(substr(t.AgeC2,1,1),levels=c('U','D'),labels=c('undiagnosed','diagnosed'))]	
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.AgeC2','t.AgeC','stageC','tmidC')]	
	ggp		<- ggplot(tmp, aes(x=t.AgeC, y=100*p.per.tc, fill=t.AgeC, alpha=stageC)) + geom_bar(stat='identity',colour='black') +
			scale_y_continuous(expand=c(0,0), limits=c(0,55)) +
			scale_alpha_manual(values=c('diagnosed'=1, 'undiagnosed'=0.25)) +
			scale_fill_hue(h=c(200, 360), l=50, c=100, guide=FALSE) +
			facet_grid(~tmidC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) +
			labs(x='\nDOB of probable transmitter', y='%', alpha='')
	z <- ggplot_gtable(ggplot_build(ggp))
	z <- gtable_add_rows(z, z$heights[[3]], 2)
	z <- gtable_add_grob(z, 
			list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
					textGrob("Transmissions by generation", gp = gpar(col = 'black'))), 3, 4, 3, 10, name = paste(runif(2)))
	z <- gtable_add_rows(z, unit(1/8, "line"), 3)
	grid.newpage()
	pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByTmidUD',suffix,'.pdf',sep=''),infile), sep=''), w=13, h=5)
	grid.draw(z)
	dev.off()
	
	#
	#	by p.AgeC
	#
	YXmp	<- merge(YXm, unique(subset(YXc, select=c(Patient, AgeC))), by='Patient')
	tmp2	<- YXmp[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.AgeC','AgeC','Patient')]
	tmp2[, p.AgeC:= paste(as.character(t.AgeC),'->',as.character(AgeC),sep='')]
	tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','p.AgeC','t.AgeC', 'AgeC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','AgeC','t.AgeC','p.AgeC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.AgeC=t.AgeC, AgeC=AgeC, p.AgeC=p.AgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.AgeC','AgeC','p.AgeC'))	
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))	
	#	p.AgeC proportions	
	tmp		<- YXpr[, list(p.per.tc=100*mean(p.per.rec), po.per.tc=100*mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=100*mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))	
	ggp		<- ggplot(tmp, aes(x=t.AgeC, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~AgeC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) + 
			scale_y_continuous(expand=c(0,0), limits=c(0,55)) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			scale_fill_hue(h=c(200, 360), l=50, c=100, guide=FALSE) +
			labs(x='\nDOB of probable transmitter', y='%', alpha='')	
	z <- ggplot_gtable(ggplot_build(ggp))
	z <- gtable_add_rows(z, z$heights[[3]], 2)
	z <- gtable_add_grob(z, 
			list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
					textGrob("Transmissions by DOB of recipient", gp = gpar(col = 'black'))), 3, 4, 3, 12, name = paste(runif(2)))
	z <- gtable_add_rows(z, unit(1/8, "line"), 3)
	grid.newpage()
	pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=18, h=5)
	grid.draw(z)
	dev.off()
	#		
	tmp		<- YXpr[, list(p.per.tc=sum(p.per.rec), po.per.tc=sum(p.per.rec*(1-pm.per.rec)), pm.per.tc=sum(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))	
	ggp		<- ggplot(tmp, aes(x=t.AgeC, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~AgeC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) + 
			scale_y_continuous(expand=c(0,0), limits=c(0,80)) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			scale_fill_hue(h=c(200, 360), l=50, c=100, guide=FALSE) +
			labs(x='\nDOB of probable transmitter', y='#', alpha='')	
	z <- ggplot_gtable(ggplot_build(ggp))
	z <- gtable_add_rows(z, z$heights[[3]], 2)
	z <- gtable_add_grob(z, 
			list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
					textGrob("Transmissions by DOB of recipient", gp = gpar(col = 'black'))), 3, 4, 3, 12, name = paste(runif(2)))
	z <- gtable_add_rows(z, unit(1/8, "line"), 3)
	grid.newpage()
	pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByAgeCTotals',suffix,'.pdf',sep=''),infile), sep=''), w=18, h=5)
	grid.draw(z)
	dev.off()

	#	p.AgeC total
	tmp		<- YXpr[, list(n.per.tc=sum(p.per.rec), no.per.tc=sum(p.per.rec*(1-pm.per.rec)), nm.per.tc=sum(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC')]
	tmp[, p.per.tc:= 100*n.per.tc/sum(n.per.tc)]
	tmp[, text:= as.character(round(100*n.per.tc/sum(n.per.tc),d=1))]
	set(tmp, tmp[, which(p.per.tc<3)], 'text', '')
	tmp[, po.per.tc:= 100*no.per.tc/sum(n.per.tc)]
	tmp[, pm.per.tc:= 100*nm.per.tc/sum(n.per.tc)]	
	setnames(tmp, c('po.per.tc','p.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('expected missing','observed'), labels=c('expected missing','observed'))])
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('observed','expected missing'), labels=c('observed','expected missing'))])
	ggplot(tmp, aes(x=AgeC, y=t.AgeC)) + 			
			geom_point(aes(size=value, alpha=variable, fill=variable), pch=21) +
			geom_text(aes(label=text), colour='black', size=4, hjust=0.5, vjust=0.5) +
			scale_size_continuous(range=c(4,50), guide=FALSE) + 
			scale_fill_manual(values=c('expected missing'='#E41A1C', 'observed'='#E41A1C')) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			theme_bw() + theme(legend.position='bottom') + labs(x='age of recipient', y='age of probable transmitter', fill='transmissions\n(%)', alpha='transmissions\n(%)')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_AgeCtAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=5, h=5)
}
######################################################################################
adjust.dev.Age_283338etc.Stage_UAC_UAE_UC_D_TS_TO_F<- function()
{
	require(survival)
	require(gtable)
	YX
	X.tables
	t.recent.endctime
	
	indir		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C2V100bInfT7STRAT_m5D.R'
	
	YXo			<- copy(YX)	#before strat
	#YXo2		<-	stratificationmodel.Age_253045.Stage_UA_UC_D_TS_TO_F(YXo, lRNA.supp=2)
	YXo2		<-	stratificationmodel.Age_2833384348.Stage_UAE_UAC_UC_D_TS_TO_F(YXo, lRNA.supp=2)
	risk.col	<- 't.stAgeC'	
	#
	#	focus on >2003 and before ART start
	#
	#YX			<- subset(YX, t>2004.5 & !stageC%in%c('L','T'))
	YX			<- copy(YXo2)
	YX			<- subset(YX, AnyPos_T1>2005)
	
	#	get clustering, sampling and non-censoring probabilities 
	#	for each prob transmitter at the transmission interval 
	sm			<- X.tables$sm
	YXs			<- merge(YX, sm, by='t.Patient')
	#	p.nc	(non-censoring prob)
	cp			<- copy(X.tables$cm$cens.p)		#model built after subsampling, so cannot apply straight-away
	set(cp, NULL, c('BS','p.nc'), NULL)
	cm			<- X.tables$cm$cens.m	
	tp.df		<- subset(YXs, grepl('^U', YXs[[risk.col]]))[, list(tm=mean(t)), by=c('t.Patient','Patient')]
	tp.df		<- censoring.get.dataset(tp.df, df.all.allmsm)
	setnames(YXs, c('t.Patient','Patient'), c('Patient','r.Patient'))
	#	YXc so far only contains transmitters with at least one undiagnosed interval
	YXc			<- merge(YXs, subset(tp.df, select=c(Patient, r.Patient, tm, AACD4C)), by=c('Patient','r.Patient'))
	set(YXc, NULL, 'tm', YXc[, tm-t.recent.endctime])	
	YXc[, p.nc:=predict(cm, data=cp, newdata=subset(YXc, select=c(tm, AACD4C)), type='response')]
	#	p.nc is allocated to transmitters with at least one undiagnosed interval
	#	set p.nc to 1 for all intervals after diagnosis
	set(YXc, YXc[, which(!is.na(p.nc) & !grepl('^U', stageC))], 'p.nc', 1)
	#ggplot(YXc, aes(x=tm, y=p.nc, group=AACD4C)) + geom_line() + facet_wrap(~AACD4C, ncol=4)
	#ggplot(YXc, aes(x=tm)) + geom_histogram()	
	setnames(YXs, c('Patient','r.Patient'), c('t.Patient','Patient'))
	setnames(YXc, c('Patient','r.Patient'), c('t.Patient','Patient'))
	YXc			<- merge(YXs, subset(YXc, select=c(t.Patient, Patient, t, p.nc)), by=c('t.Patient','Patient','t'), all.x=TRUE)
	set(YXc, YXc[, which(is.na(p.nc))], 'p.nc', 1)
	#	get p.clu and p.pt
	risk.table	<- X.tables$st$risk.table
	tmp			<- clustering.getrisks(risk.table, indir=indir, infile=infile, rskf.name='t.stAgeC', rskf.baseline='D_(30,45]')
	tpt			<- tmp$tpt
	tclu		<- tmp$tclu	
	tmp			<- subset(tclu, t.period=='Overall', select=c(factor2, p.clu))
	setnames(tmp, 'factor2', 't.stAgeC')	
	set(tmp, NULL, 't.stAgeC', tmp[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)], labels=YXc[, levels(t.stAgeC)])])
	YXc			<- merge(YXc, tmp, by=c('t.stAgeC'), all.x=TRUE)
	tmp			<- subset(tpt, t.period=='Overall', select=c(factor2, p.pt))
	setnames(tmp, 'factor2', 't.stAgeC')		
	set(tmp, NULL, 't.stAgeC', tmp[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)], labels=YXc[, levels(t.stAgeC)])])
	YXc			<- merge(YXc, tmp, by=c('t.stAgeC'), all.x=TRUE)	
	#	get w per interval
	cat(paste('\nsetting likelihood to likelihood of pair / number transmission intervals'))
	set(YXc, NULL, 'score.Y.raw', YXc[, score.Y])
	set(YXc, NULL, 'score.Y', YXc[, score.Y*w.tn])				
	YXc		<- merge(YXc, YXc[, list(t=t, t.Patient=t.Patient, ntPatient= length(unique(t.Patient))), by=c('Patient')], by=c('Patient','t.Patient','t'))	
	YXc[, p.seqnc:= p.seq*p.nc]
	YXc[, me:= 1/(p.seq*p.nc)-1]	
	#	calculate absolute transmission probabilities
	#	need first: predicted score by stage and time t
	#	need second: predicted missing intervals per recipient, by stage and time t
	#	then calculate: exp missing score.Y per recipient by time t= sum_stage missexp.pr(stage,t) * scoreY.pr(stage,t)
	tmp		<- rltvtp.model.time.150923(YXc)
	rltvp	<- tmp$prd
	rltvd	<- tmp$fit
	rltvm	<- tmp$mo
	YXc		<- subset(YXc, AnyPos_T1<2010.8)
	tmp		<- miss.model.150923(YXc)
	mp		<- tmp$prd
	mma		<- tmp$moa
	mmb		<- tmp$mob
	md		<- tmp$fit
	#YXm		<- miss.smoothing.150909(YXc)
	#YXm		<- miss.smoothing.150914(YXc, mm, md, rltvm, rltvd, indir=NA, infile=NA, suffix='150914')
	#YXm		<- miss.smoothing.150916(YXc, mm, md, rltvp, rltva, indir=indir, infile=infile, suffix='150916')
	#YXm		<- miss.smoothing.150917(YXc, mm, md, rltvp, rltva, indir=indir, infile=infile, suffix='150917', use.obs.pseqnc=TRUE)
	#YXm		<- miss.smoothing.150918(YXc, mm, md, rltvp, rltvd, rltvm, indir=indir, infile=infile, suffix='150918', use.obs.pseqnc=TRUE)
	YXm		<- miss.smoothing.150923(YXc, mp, mma, mmb, md, rltvp, rltvd, rltvm, indir=indir, infile=infile, suffix='150922', use.obs.pseqnc=TRUE)
	YXm		<- merge(YXm, unique(subset(YXc, select=c(Patient, AnyPos_T1, t.period))), by='Patient')
	if(0)
	{
		#	get expected missing term for denominator 'scoreYm.per.rec'
		tmp		<- YXm[, list(miss.per.rec= sum(miss.prsm), scoreYm.per.rec= sum(miss.prsm*scoreY.pr)), by='Patient']
		YXc		<- merge( YXc, tmp, by='Patient' )
		#	then calculate p = w / ( sum(w) + scoreYm.per.rec)
		YXc		<- merge(YXc, YXc[, list(score.p= score.Y/(sum(score.Y)+scoreYm.per.rec), score.p.nadj= score.Y/sum(score.Y), t.Patient=t.Patient, t=t), by='Patient'], by=c('Patient','t.Patient','t'))
		#ggplot(YXc, aes(x=AnyPos_T1, y=miss.per.rec)) + geom_point()		
	}
	
	#prop.consistency(YXc, YXm)
	
	tmp2	<- YXm[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.stAgeC','Patient')]	
	tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','t.stAgeC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','t.stAgeC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.stAgeC=t.stAgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.stAgeC'))	
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))
	YXpr	<- merge(YXpr, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	if(0)
	{
		#	consistency	--> aggregating first by t.stAgeC leads to slightly diff results
		tmp		<- YXpr[, list(p.per.rec=sum(p.per.rec), po.per.rec=sum(p.per.rec*(1-pm.per.rec)), pm.per.rec=sum(p.per.rec*pm.per.rec)), by=c('Patient','stageC','t.period')]	
		tmp		<- tmp[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(po.per.rec), pm.per.tc=mean(pm.per.rec)), by=c('stageC','t.period')]
		setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
		tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
		#tmp[,sum(p.per.tc),by='tmidC']	
		ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
				facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
		ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodStageC',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=7)
		
	}
	if(0)
	{
		#	consistency: combine UAE and UAC		
		tmp		<- copy(YXpr)
		set(tmp, NULL, 'stageC', tmp[, gsub('UAE|UAC','UA', stageC)])
		tmp		<- tmp[, list(p.per.rec=sum(p.per.rec), po.per.rec=sum(p.per.rec*(1-pm.per.rec)), pm.per.rec=sum(p.per.rec*pm.per.rec)), by=c('Patient','stageC','t.period')]
		tmp		<- tmp[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(po.per.rec), pm.per.tc=mean(pm.per.rec)), by=c('stageC','t.period')]
		setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
		tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
		#tmp[,sum(p.per.tc),by='tmidC']	
		ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
				facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
		ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodStageC_UA',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=7)		
	}	
	if(0)
	{
		#	check for odd behaviour in underlying model
		#	by care stage and age; obs/missing
		#	by t.stAgeC	
		tmp		<- YXpr[, list(stageC=stageC[1], t.AgeC=t.AgeC[1], p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.stAgeC','t.period')]
		setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
		tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
		#tmp[,sum(p.per.tc),by='tmidC']	
		ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
				facet_grid(t.AgeC~stageC) + theme_bw() + theme(legend.position='bottom')
		ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodtstAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=10)		
	}
	if(0)
	{
		#	transmissions by age group and the four time periods; obs/exp missing
		#	by t.AgeC
		tmp2	<- YXm[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.AgeC','Patient')]	
		tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','t.AgeC')]
		tmp2	<- merge(tmp2, tmp3, by=c('Patient','t.AgeC','t.period'), all.x=1)
		set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
		YXpr	<- tmp2
		YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
		tmp		<- YXpr[, list(t.AgeC=t.AgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
		YXpr	<- merge(YXpr, tmp, by=c('Patient','t.AgeC'))	
		YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))
		tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','t.period')]
		setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
		tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
		#tmp[,sum(p.per.tc),by='tmidC']	
		ggplot(tmp, aes(x=t.period, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
				facet_grid(~t.AgeC) + theme_bw() + theme(legend.position='bottom')
		ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodtAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=7)		
	}
	if(0)
	{
		#
		#	transmissions by calendar year (<2007,..,2009-2010)
		#	obs/exp missing
		#
		YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009,Inf), labels=c('<2007','2007','2008','>=2009'))]
		#YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009, 2010, Inf), labels=c('<2007','2007','2008','2009','2010'))]
		tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','tmidC')]	
		setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
		tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
		#tmp[,sum(p.per.tc),by='tmidC']	
		ggp		<- ggplot(tmp, aes(x=t.AgeC, y=100*value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') +
				scale_y_continuous(expand=c(0,0), limits=c(0,55)) +
				scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
				scale_fill_brewer(palette='Set1', guide=FALSE) +
				facet_grid(~tmidC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) +
				labs(x='\nage of probable transmitter', y='%', alpha='')
		z <- ggplot_gtable(ggplot_build(ggp))
		z <- gtable_add_rows(z, z$heights[[3]], 2)
		z <- gtable_add_grob(z, 
				list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
						textGrob("Transmissions by calendar year", gp = gpar(col = 'black'))), 3, 4, 3, 10, name = paste(runif(2)))
		z <- gtable_add_rows(z, unit(1/8, "line"), 3)
		grid.newpage()
		pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByTmid',suffix,'.pdf',sep=''),infile), sep=''), w=12, h=5)
		grid.draw(z)
		dev.off()
	}
	#	THIS IS THE LATEST CODE
	#
	#	transmissions by age of recipient and age of transmitter with undiagnosed/diagnosed (prepare)
	#	by p.AgeC
	#	
	YXmp	<- copy(YXm)
	#	add t.age x diagnosis status; 
	#	determine which calendar year each recipient belongs to by midpoint of its transmission window
	#	add age of recipient
	YXmp[, t.AgeC2:= YXmp[, gsub('TO|TS|L','D',gsub('U[^_]+','U',t.stAgeC))]]
	tmp		<- YXc[, list(tmid= mean(t)),  by=c('Patient')]
	YXmp	<- merge(YXmp, tmp, by='Patient')
	YXmp[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009,Inf), labels=c('<2007','2007','2008','>=2009'))]
	YXmp	<- merge(YXmp, unique(subset(YXc, select=c(Patient, AgeC))), by='Patient')
	#	evaluate missing trm scores per patient, age stages, tmidC 	
	tmp2	<- YXmp[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('tmidC','t.AgeC2','t.AgeC','AgeC','Patient')]
	tmp2[, p.AgeC:= paste(as.character(t.AgeC),'->',as.character(AgeC),sep='')]
	set(tmp2, NULL, 't.AgeC2', tmp2[, substr(t.AgeC2,1,1)])
	#	evaluate observed trm scores per patient, age stages, tmidC
	tmp3	<- copy(YXc)	
	tmp3[, t.AgeC2:= tmp3[, gsub('TO|TS|L','D',gsub('U[^_]+','U',t.stAgeC))]]
	tmp		<- YXc[, list(tmid= mean(t)),  by=c('Patient')]
	tmp3	<- merge(tmp3, tmp, by='Patient')
	tmp3[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009,Inf), labels=c('<2007','2007','2008','>=2009'))]
	tmp3	<- tmp3[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('tmidC','t.AgeC2','t.AgeC', 'AgeC','p.AgeC','Patient')]
	set(tmp3, NULL, 't.AgeC2', tmp3[, substr(t.AgeC2,1,1)])
	#	merge observed and missing trm scores
	tmp2	<- merge(tmp2, tmp3, by=c('tmidC','t.AgeC2','t.AgeC','AgeC','p.AgeC','Patient'), all=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	#subset(tmp2, is.na(scoreYm.per.recst))
	set(tmp2, tmp2[, which(is.na(scoreYm.per.recst))], 'scoreYm.per.recst', 0) 	#TODO this should not happen! the recipients transition age classes, so maybe I drop some of these in the construction of YXm??	 
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	setnames(YXpr, 't.AgeC2','stageC')
	#	calculate trm probs
	YXpr	<- dcast.data.table(YXpr, tmidC+t.AgeC+AgeC+p.AgeC+Patient~stageC, value.var='scoreY.per.recst')	
	set(YXpr, YXpr[, which(is.na(D))],'D',0)
	set(YXpr, YXpr[, which(is.na(U))],'U',0)
	YXpr	<- YXpr[, list(	t.AgeC=t.AgeC, AgeC=AgeC, p.AgeC=p.AgeC, 
					p.per.rec= (D+U)/sum(D+U), pD.per.rec= D/sum(D+U), pU.per.rec=U/sum(D+U)), by=c('tmidC','Patient')]	
	#	get the distribution of recipients by age group
	#	tn and tna is number of transmission intervals per recipient and age group	
	tmp3	<- YXc[, list(tna=length(unique(t))), by=c('Patient','AgeC')]
	tmp3	<- merge(tmp3, YXc[, list(tn=length(unique(t))), by=c('Patient')], by=c('Patient'))
	YXpr	<- merge(YXpr, tmp3, by=c('Patient','AgeC'))	
	#	make pretty
	YXpr[, tmidC2:= '2004-2007']
	set(YXpr, YXpr[, which(tmidC%in%c('2008','>=2009'))],'tmidC2','2008-2010')	
	tmp		<- YXpr[, levels(AgeC)]
	set(YXpr, NULL, 'AgeC', YXpr[, factor(as.character(AgeC), levels=tmp, labels=c('16-27','28-32','33-37','38-42','43-47','48-80'))])
	set(YXpr, NULL, 't.AgeC', YXpr[, factor(as.character(t.AgeC), levels=tmp, labels=c('16-27','28-32','33-37','38-42','43-47','48-80'))])	
	YXpr	<- merge(YXpr, YXpr[, list(nRec=length(unique(Patient))), by='tmidC2'], by='tmidC2')
	#
	#	calculate proportion of recipients by age in 2004-2007 and 2008-2010 and
	#	compare to distribution among new diagnoses
	#
	YXar	<- copy(YXpr)	
	setkey(YXar, tmidC2, Patient, AgeC)
	YXar	<- unique(subset(YXar, select=c(tmidC2, Patient, AgeC, tn, tna, nRec)))
	set(YXar, NULL, 'tna', YXar[, tna/tn])
	YXar	<- YXar[, list(TYPE='Recipient', N= sum(tna), Ntp=nRec[1]), by=c('tmidC2','AgeC')]
	#	add distribution among new diagnoses
	tmp		<- subset(df.all.allmsm, AnyPos_T1>2004 & AnyPos_T1<2010 & Trm%in%c('MSM','BI'), c(Patient, AnyPos_T1, DateBorn, PosSeqT))
	setkey(tmp, Patient)
	tmp		<- unique(tmp)
	tmp[, Age:=AnyPos_T1-DateBorn]
	tmp[, tmidC2:= cut(AnyPos_T1, breaks=c(-Inf, 2008, Inf), labels=c('2004-2007','2008-2010'))]
	tmp[, AgeC:= cut(Age, breaks=c(-Inf, 28, 33, 38, 43, 48, Inf), labels=c('16-27','28-32','33-37','38-42','43-47','48-80'))]
	tmp2	<- tmp[, list(TYPE='Diagnosed', N= length(Patient)), by=c('tmidC2','AgeC')]
	tmp2	<- merge(tmp2, tmp2[, list(Ntp=sum(N)),by='tmidC2'], by='tmidC2')
	YXar	<- rbind(YXar, tmp2)
	#	add distribution among sequenced
	tmp2	<- subset(tmp, !is.na(PosSeqT))[, list(TYPE='Sampled', N= length(Patient)), by=c('tmidC2','AgeC')]	
	tmp2	<- merge(tmp2, tmp2[, list(Ntp=sum(N)),by='tmidC2'], by='tmidC2')
	YXar	<- rbind(YXar, tmp2)	
	YXar[, P:=N/Ntp]
	#	plot
	YXplot	<- copy(YXar)
	set(YXplot, NULL, 'TYPE', YXplot[, factor(TYPE, levels=c('Recipient','Sampled','Diagnosed'), labels=c('Recipient MSM\n(of whom sources could be characterised)','Newly diagnosed MSM\nwith a sequence','Newly diagnosed MSM\n(population)') )])
	set(YXplot, NULL, 'tmidC2', YXplot[, factor(as.character(tmidC2, levels=c('2004-2007','2008-2020'),labels=c('2004-2007\n(time of diagnosis)','2008-2010\n(time of diagnosis)')))])
	ggplot(YXplot, aes(y=100*P, x=AgeC, fill=TYPE)) +
			geom_bar(stat='identity',position='dodge', colour='black', width=0.8) + 
			facet_grid(~tmidC2) + 
			scale_y_continuous(expand=c(0,0), limit=c(0,24)) +
			#scale_fill_manual(values=c('Newly diagnosed MSM\n(population)'="#FEB24C", 'Recipient MSM\n(of whom sources could be characterised)'="#E31A1C")) +
			scale_fill_brewer(palette='Pastel1') +
			coord_flip() +
			theme_bw() + theme(legend.position='bottom') +
			labs(x='Age group\n(years)\n', y='%', fill='', title='Study population\n')			
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_Rec_by_AgeC_',suffix,'.pdf',sep=''),infile), sep=''), w=8.5, h=4)
	#	calculate adjustments to transform the age distribution of the recipients to that of the newly diagnosed
	YXar	<- dcast.data.table(YXar, tmidC2+AgeC~TYPE, value.var='P')
	YXar[, wRec:= Diagnosed/Recipient]
	YXpr	<- merge(YXpr, YXar, by=c('tmidC2','AgeC'))	
	#
	#	UNADJUSTED aggregate by p.AgeC (age of transmitters and recipients) and diagnosis status
	#	for years <2007 and 2008-2010
	#	
	YXplot	<- YXpr[, list(p.per.tc=sum(p.per.rec)/nRec[1], pD.per.tc=sum(pD.per.rec)/nRec[1]), by=c('tmidC2','t.AgeC','AgeC','p.AgeC')]
	set(YXplot, NULL, 'tmidC2', YXplot[, factor(tmidC2, levels=c('2004-2007','2008-2010'), labels=c('2004-2007\n(midpoint of seroconversion window of recipient)','2008-2010\n(midpoint of seroconversion window of recipient)'))])
	YXplot[, text:= as.character(round(100*p.per.tc,d=1))]
	YXplot[, text.size:=4]
	set(YXplot, YXplot[, which(p.per.tc<0.03)], 'text.size', 2)
	YXplot	<- melt(YXplot, measure.vars=c('p.per.tc','pD.per.tc'))
	set(YXplot, NULL, 'variable', YXplot[,factor(variable, levels=c('p.per.tc','pD.per.tc'), labels=c('undiagnosed','diagnosed'))])
	ggplot(YXplot, aes(x=AgeC, y=t.AgeC)) + 			
			geom_point(aes(size=100*value, alpha=variable, fill=100*as.numeric(text)), pch=21) +
			geom_text(data=subset(YXplot, text.size==4), aes(label=text), size=4, colour='black', hjust=0.5, vjust=0.5) +
			geom_text(data=subset(YXplot, text.size==2), aes(label=text), size=2, colour='black', hjust=0.5, vjust=0.5) +
			scale_size_continuous(range=c(2,30), guide=FALSE) +
			scale_fill_gradient(low='LightBlue', high='red', guide=FALSE) +			
			scale_alpha_manual(values=c('undiagnosed'=0.4, 'diagnosed'=0.8)) +
			facet_grid(~tmidC2) +
			theme_bw() + theme(panel.margin = unit(4/8, "lines"), legend.position='bottom') + 
			labs(	x='\nage of recipient\n(years at probable transmission interval)', 
					y='age of probable transmitter\n(years at probable transmission interval)\n', 
					fill='transmissions\n(%)', alpha='transmissions\n(%)',
					title='Transmissions between age groups\n')	
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_AgeCtAgeCByTime_DU_',suffix,'.pdf',sep=''),infile), sep=''), w=8.5, h=5.5)
	#
	#	ADJUSTED aggregate by p.AgeC (age of transmitters and recipients) and diagnosis status
	#	for years <2007 and 2008-2010
	#	
	YXplot	<- YXpr[, list(p.per.tc=sum(wRec*p.per.rec)/nRec[1], pD.per.tc=sum(wRec*pD.per.rec)/nRec[1]), by=c('tmidC2','t.AgeC','AgeC','p.AgeC')]
	#set(YXplot, NULL, 'tmidC2', YXplot[, factor(tmidC2, levels=c('2004-2007','2008-2010'), labels=c('2004-2007\n(midpoint of seroconversion window of recipient)','2008-2010\n(midpoint of seroconversion window of recipient)'))])
	YXplot[, text:= as.character(round(100*p.per.tc,d=1))]
	YXplot[, text.size:=4]
	set(YXplot, YXplot[, which(p.per.tc<0.03)], 'text.size', 2)
	YXplot	<- melt(YXplot, measure.vars=c('p.per.tc','pD.per.tc'))
	set(YXplot, NULL, 'variable', YXplot[,factor(variable, levels=c('p.per.tc','pD.per.tc'), labels=c('undiagnosed','diagnosed'))])
	ggplot(YXplot, aes(x=AgeC, y=t.AgeC)) + 			
			geom_point(aes(size=100*value, alpha=variable, fill=100*as.numeric(text)), pch=21) +
			geom_text(data=subset(YXplot, text.size==4), aes(label=text), size=4, colour='black', hjust=0.5, vjust=0.5) +
			geom_text(data=subset(YXplot, text.size==2), aes(label=text), size=2, colour='black', hjust=0.5, vjust=0.5) +
			scale_size_continuous(range=c(2,30), guide=FALSE) +
			scale_fill_gradient(low='LightBlue', high='red', guide=FALSE) +			
			scale_alpha_manual(values=c('undiagnosed'=0.4, 'diagnosed'=0.8)) +
			facet_grid(~tmidC2) +
			theme_bw() + theme(panel.margin = unit(4/8, "lines"), legend.position='bottom') + 
			labs(	x='\nage of infected MSM\n', 
					y='age of transmitter\n', 
					fill='transmissions\n(%)', alpha='transmissions\n(%)',
					title='Transmissions between age groups\n')	
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_AgeCtAgeCByTime_DU_ADJ_',suffix,'.pdf',sep=''),infile), sep=''), w=8.5, h=5.5)
	#
	#ggplot(YXplot, aes(x=tmidC2, weight=tna, fill=AgeC)) +
	#		geom_bar(position='fill', colour='black') +
	#		scale_fill_brewer(palette='Reds') +
	#		scale_x_discrete(expand=c(0,0)) +
	#		scale_y_continuous(expand=c(0,0), labels=function(x) x*100) +
	#		theme_bw() + theme(legend.position='bottom') +
	#		labs(	x='', y='%\n', fill='age', title='Recipients by age\n')
	#ggsave(file=paste(indir, '/', gsub('\\.R',paste('_Rec_by_AgeC_',suffix,'.pdf',sep=''),infile), sep=''), w=3.5, h=5.5)
	#
	#	UNADJUSTED transmissions by age of transmitter in 2004-2007 and 2008-2010 and diagnosis status
	#
	YXplot	<- YXpr[, list(p.per.tc=sum(p.per.rec)/nRec[1], pD.per.tc=sum(pD.per.rec)/nRec[1], pU.per.tc=sum(wRec*pU.per.rec)/nRec[1] ), by=c('tmidC2','t.AgeC')]
	YXplot	<- melt(YXplot, measure.vars=c('pD.per.tc','pU.per.tc'))
	set(YXplot, NULL, 'variable', YXplot[, factor(variable, levels=c('pU.per.tc','pD.per.tc'), labels=c('before diagnosis','after diagnosis'))])
	ggplot(YXplot, aes(x=t.AgeC, y=100*value, fill=variable)) +
			geom_bar(stat='identity', position='stack', colour='black') +
			#scale_x_discrete(expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,30,5), minor_breaks=seq(0,30,1)) +
			scale_fill_manual(values=c('before diagnosis'="#CAB2D6", 'after diagnosis'="#6A3D9A")) +
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(~tmidC2) +
			labs(	x='\nage of probable transmitter\n(years)', y='%\n', fill='diagnosis status of transmitter\n(at probable transmission interval)   ', title='Transmissions from age groups\n')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByTime_DU_',suffix,'.pdf',sep=''),infile), sep=''), w=8.5, h=5.5)
	#
	#	ADJUSTED transmissions by age of transmitter in 2004-2007 and 2008-2010 and diagnosis status
	#
	YXplot	<- YXpr[, list(p.per.tc=sum(wRec*p.per.rec)/nRec[1], pD.per.tc=sum(wRec*pD.per.rec)/nRec[1], pU.per.tc=sum(wRec*pU.per.rec)/nRec[1] ), by=c('tmidC2','t.AgeC')]
	YXplot	<- melt(YXplot, measure.vars=c('pD.per.tc','pU.per.tc'))
	set(YXplot, NULL, 'variable', YXplot[, factor(variable, levels=c('pU.per.tc','pD.per.tc'), labels=c('before diagnosis','after diagnosis'))])
	ggplot(YXplot, aes(x=t.AgeC, y=100*value, fill=variable)) +
			geom_bar(stat='identity', position='stack', colour='black') +
			#scale_x_discrete(expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,30,5), minor_breaks=seq(0,30,1)) +
			scale_fill_manual(values=c('before diagnosis'="#CAB2D6", 'after diagnosis'="#6A3D9A")) +
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(~tmidC2) +
			labs(	x='\nage of transmitter\n(years)', y='%\n', fill='diagnosis status\nof transmitter  ', title='Transmissions from age groups\n')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByTime_DU_ADJ_',suffix,'.pdf',sep=''),infile), sep=''), w=8.5, h=5.5)
	#
	#	transmissions by age of recipient and age of transmitter with observed/exp missing (prepare)
	#	by p.AgeC
	#
	YXmp	<- merge(YXm, unique(subset(YXc, select=c(Patient, AgeC))), by='Patient')
	tmp2	<- YXmp[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.AgeC','AgeC','Patient')]
	tmp2[, p.AgeC:= paste(as.character(t.AgeC),'->',as.character(AgeC),sep='')]
	tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','p.AgeC','t.AgeC', 'AgeC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','AgeC','t.AgeC','p.AgeC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.AgeC=t.AgeC, AgeC=AgeC, p.AgeC=p.AgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.AgeC','AgeC','p.AgeC'))	
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))	
	if(0)
	{
		#	transmissions by age of recipient and age of transmitter (totals)
		#	p.AgeC totals
		tmp		<- YXpr[, list(p.per.tc=sum(p.per.rec), po.per.tc=sum(p.per.rec*(1-pm.per.rec)), pm.per.tc=sum(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC')]
		setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
		tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))	
		ggp		<- ggplot(tmp, aes(x=t.AgeC, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
				facet_grid(~AgeC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) + 
				scale_y_continuous(expand=c(0,0), limits=c(0,40), breaks=seq(0,200,20), minor_breaks=seq(0,200,5)) +
				scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
				scale_fill_brewer(palette='Set1', guide=FALSE) +
				labs(x='\nage of probable transmitter', y='#', alpha='')	
		z <- ggplot_gtable(ggplot_build(ggp))
		z <- gtable_add_rows(z, z$heights[[3]], 2)
		z <- gtable_add_grob(z, 
				list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
						textGrob("Transmissions by age of recipients in cohort", gp = gpar(col = 'black'))), 3, 4, 3, 14, name = paste(runif(2)))
		z <- gtable_add_rows(z, unit(1/8, "line"), 3)
		grid.newpage()
		pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByAgeCTotals',suffix,'.pdf',sep=''),infile), sep=''), w=19, h=4)
		grid.draw(z)
		dev.off()
		#	transmissions by age of recipient and age of transmitter (proportions)
		#	p.AgeC proportions	
		tmp		<- YXpr[, list(p.per.tc=100*mean(p.per.rec), po.per.tc=100*mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=100*mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC')]
		setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
		tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))	
		ggp		<- ggplot(tmp, aes(x=t.AgeC, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
				facet_grid(~AgeC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) + 
				scale_y_continuous(expand=c(0,0), limits=c(0,35)) +
				scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
				scale_fill_brewer(palette='Set1', guide=FALSE) +
				labs(x='\nage of probable transmitter', y='%', alpha='')	
		z <- ggplot_gtable(ggplot_build(ggp))
		z <- gtable_add_rows(z, z$heights[[3]], 2)
		z <- gtable_add_grob(z, 
				list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
						textGrob("Transmissions by age of recipient", gp = gpar(col = 'black'))), 3, 4, 3, 14, name = paste(runif(2)))
		z <- gtable_add_rows(z, unit(1/8, "line"), 3)
		grid.newpage()
		pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=19, h=4)
		grid.draw(z)
		dev.off()
	}
	if(0)
	{
		#	the dot plot by observed/exp missing
		#	p.AgeC total by <=2007, >2007
		tmp		<- copy(YXpr)
		tmp[, tmidC:= cut(tmid, breaks=c(-Inf, 2008, Inf), labels=c('2004-2007','2008-2010'))]
		tmp		<- tmp[, list(n.per.tc=sum(p.per.rec), no.per.tc=sum(p.per.rec*(1-pm.per.rec)), nm.per.tc=sum(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC','tmidC')]
		tmp		<- tmp[, list(t.AgeC=t.AgeC, AgeC=AgeC, p.AgeC=p.AgeC, p.per.tc= 100*n.per.tc/sum(n.per.tc), po.per.tc= 100*no.per.tc/sum(n.per.tc), pm.per.tc= 100*nm.per.tc/sum(n.per.tc)), by='tmidC']
		tmp[, text:= as.character(round(p.per.tc,d=1))]
		tmp[, text.size:=4]
		set(tmp, tmp[, which(p.per.tc<3)], 'text.size', 2)
		tmp[, total:=p.per.tc]
		setnames(tmp, c('po.per.tc','p.per.tc'), c('observed','expected missing'))
		tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
		#set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('expected missing','observed'), labels=c('expected missing','observed'))])
		set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('observed','expected missing'), labels=c('observed','expected missing'))])
		ggp		<- ggplot(tmp, aes(x=AgeC, y=t.AgeC)) + 			
				geom_point(aes(size=value, alpha=variable, fill=total), pch=21) +
				geom_text(data=subset(tmp, text.size==4), aes(label=text), size=4, colour='black', hjust=0.5, vjust=0.5) +
				geom_text(data=subset(tmp, text.size==2), aes(label=text), size=2, colour='black', hjust=0.5, vjust=0.5) +
				scale_size_continuous(range=c(4,50), guide=FALSE) +
				scale_fill_gradient(low='blue', high='red', guide=FALSE) +
				#scale_fill_manual(values=c('expected missing'='#E41A1C', 'observed'='#E41A1C')) +
				scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=0.6)) +
				facet_grid(~tmidC) +
				theme_bw() + theme(panel.margin = unit(4/8, "lines"), legend.position='bottom') + labs(x='\nage of recipient', y='age of probable transmitter', fill='transmissions\n(%)', alpha='transmissions\n(%)')
		z <- ggplot_gtable(ggplot_build(ggp))
		z <- gtable_add_rows(z, z$heights[[3]], 2)
		z <- gtable_add_grob(z, 
				list(	rectGrob(height = unit(2, "npc"), gp = gpar(col='black', fill=gray(0.8))),
						textGrob("Transmissions by age", gp = gpar(col = 'black'))), 3, 4, 3, 6, name = paste(runif(2)))
		z <- gtable_add_rows(z, unit(2/8, "line"), 3)
		grid.newpage()
		pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_AgeCtAgeCByTime',suffix,'.pdf',sep=''),infile), sep=''), w=8.5, h=5.5)
		grid.draw(z)
		dev.off()
	}
	if(0)
	{	
		#	the dot plots, totals	
		#	p.AgeC total
		tmp		<- YXpr[, list(n.per.tc=sum(p.per.rec), no.per.tc=sum(p.per.rec*(1-pm.per.rec)), nm.per.tc=sum(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC')]
		tmp[, p.per.tc:= 100*n.per.tc/sum(n.per.tc)]
		tmp[, text:= as.character(round(100*n.per.tc/sum(n.per.tc),d=1))]
		tmp[, text.size:=4]
		set(tmp, tmp[, which(p.per.tc<3)], 'text.size', 2)
		tmp[, po.per.tc:= 100*no.per.tc/sum(n.per.tc)]
		tmp[, pm.per.tc:= 100*nm.per.tc/sum(n.per.tc)]	
		setnames(tmp, c('po.per.tc','p.per.tc'), c('observed','expected missing'))
		tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
		#set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('expected missing','observed'), labels=c('expected missing','observed'))])
		set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('observed','expected missing'), labels=c('observed','expected missing'))])
		ggplot(tmp, aes(x=AgeC, y=t.AgeC)) + 			
				geom_point(aes(size=value, alpha=variable, fill=variable), pch=21) +
				geom_text(data=subset(tmp, text.size==4), aes(label=text), size=4, colour='black', hjust=0.5, vjust=0.5) +
				geom_text(data=subset(tmp, text.size==2), aes(label=text), size=2, colour='black', hjust=0.5, vjust=0.5) +
				scale_size_continuous(range=c(4,50), guide=FALSE) + 
				scale_fill_manual(values=c('expected missing'='#E41A1C', 'observed'='#E41A1C')) +
				scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
				theme_bw() + theme(legend.position='bottom') + labs(x='age of recipient', y='age of probable transmitter', fill='transmissions\n(%)', alpha='transmissions\n(%)')
		ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_AgeCtAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=5, h=5)	
	}
	
}
######################################################################################
adjust.dev.Age_25354555.Stage_UAC_UAE_UC_D_TS_TO_F<- function()
{
	require(survival)
	YX
	X.tables
	t.recent.endctime
	
	indir		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C2V100bInfT7STRAT_m5E.R'
	
	YXo			<- copy(YX)	#before strat
	#YXo2		<-	stratificationmodel.Age_253045.Stage_UA_UC_D_TS_TO_F(YXo, lRNA.supp=2)
	YXo2		<-	stratificationmodel.Age_25354555.Stage_UAE_UAC_UC_D_TS_TO_F(YXo, lRNA.supp=2)
	risk.col	<- 't.stAgeC'	
	#
	#	focus on >2003 and before ART start
	#
	#YX			<- subset(YX, t>2004.5 & !stageC%in%c('L','T'))
	YX			<- copy(YXo2)
	YX			<- subset(YX, AnyPos_T1>2005)
	
	#	p.nc	(non-censoring prob)
	cp			<- copy(X.tables$cm$cens.p)		#model built after subsampling, so cannot apply straight-away
	set(cp, NULL, c('BS','p.nc'), NULL)
	cm			<- X.tables$cm$cens.m	
	
	load('~/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C2V100bInfT7STRAT_m5E.R')
	tp.df		<- subset(X.msm, grepl('^U', X.msm[[risk.col]]))[, list(tm=mean(t)), by=c('t.Patient','Patient')]
	tp.df		<- censoring.get.dataset(tp.df, df.all.allmsm)
	setnames(X.msm, c('t.Patient','Patient'), c('Patient','r.Patient'))
	#	YXc so far only contains transmitters with at least one undiagnosed interval
	tmp			<- merge(X.msm, subset(tp.df, select=c(Patient, r.Patient, tm, AACD4C)), by=c('Patient','r.Patient'))
	set(tmp, NULL, 'tm', tmp[, tm-t.recent.endctime])	
	tmp[, p.nc:=predict(cm, data=cp, newdata=subset(tmp, select=c(tm, AACD4C)), type='response')]	
	X.msm		<- merge(X.msm, subset(tmp, select=c(Patient, r.Patient, t, p.nc)), by=c('Patient','r.Patient','t'), all.x=1)
	tmp			<- NULL
	gc()
	set(X.msm, X.msm[, which(is.na(p.nc) | !grepl('^U', stageC))], 'p.nc', 1)
	X.msm[, me:= 1/(p.nc)-1]	
	save(X.msm, file=paste(indir,'/',gsub('\\.R','_CENSMSM.R',infile),sep=''))
	setnames(X.msm, c('Patient','r.Patient'), c('t.Patient','Patient'))
	tmp			<- X.msm[, list(tm=mean(t)), by=c('t.Patient','Patient')]
	tmp[, tmidC:= cut(tm, breaks=c(-Inf, 2007, 2008, 2009,Inf), labels=c('<2007','2007','2008','>=2009'))]
	X.msm		<- merge(X.msm, tmp, by=c('t.Patient','Patient'))
	Xmn			<- X.msm[, list(to=length(t), tm=sum(me)), by=c('t.AgeC','stageC','t.AgeC','tmidC')]
	
	
	#	get clustering, sampling and non-censoring probabilities 
	#	for each prob transmitter at the transmission interval 
	sm			<- X.tables$sm
	YXs			<- merge(YX, sm, by='t.Patient')
	tp.df		<- subset(YXs, grepl('^U', YXs[[risk.col]]))[, list(tm=mean(t)), by=c('t.Patient','Patient')]
	tp.df		<- censoring.get.dataset(tp.df, df.all.allmsm)
	setnames(YXs, c('t.Patient','Patient'), c('Patient','r.Patient'))
	#	YXc so far only contains transmitters with at least one undiagnosed interval
	YXc			<- merge(YXs, subset(tp.df, select=c(Patient, r.Patient, tm, AACD4C)), by=c('Patient','r.Patient'))
	set(YXc, NULL, 'tm', YXc[, tm-t.recent.endctime])	
	YXc[, p.nc:=predict(cm, data=cp, newdata=subset(YXc, select=c(tm, AACD4C)), type='response')]
	#	p.nc is allocated to transmitters with at least one undiagnosed interval
	#	set p.nc to 1 for all intervals after diagnosis
	set(YXc, YXc[, which(!is.na(p.nc) & !grepl('^U', stageC))], 'p.nc', 1)
	#ggplot(YXc, aes(x=tm, y=p.nc, group=AACD4C)) + geom_line() + facet_wrap(~AACD4C, ncol=4)
	#ggplot(YXc, aes(x=tm)) + geom_histogram()	
	setnames(YXs, c('Patient','r.Patient'), c('t.Patient','Patient'))
	setnames(YXc, c('Patient','r.Patient'), c('t.Patient','Patient'))
	YXc			<- merge(YXs, subset(YXc, select=c(t.Patient, Patient, t, p.nc)), by=c('t.Patient','Patient','t'), all.x=TRUE)
	set(YXc, YXc[, which(is.na(p.nc))], 'p.nc', 1)
	#	get p.clu and p.pt
	risk.table	<- X.tables$st$risk.table
	tmp			<- clustering.getrisks(risk.table, indir=indir, infile=infile, rskf.name='t.stAgeC', rskf.baseline='D_(30,45]')
	tpt			<- tmp$tpt
	tclu		<- tmp$tclu	
	tmp			<- subset(tclu, t.period=='Overall', select=c(factor2, p.clu))
	setnames(tmp, 'factor2', 't.stAgeC')	
	set(tmp, NULL, 't.stAgeC', tmp[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)], labels=YXc[, levels(t.stAgeC)])])
	YXc			<- merge(YXc, tmp, by=c('t.stAgeC'), all.x=TRUE)
	tmp			<- subset(tpt, t.period=='Overall', select=c(factor2, p.pt))
	setnames(tmp, 'factor2', 't.stAgeC')		
	set(tmp, NULL, 't.stAgeC', tmp[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)], labels=YXc[, levels(t.stAgeC)])])
	YXc			<- merge(YXc, tmp, by=c('t.stAgeC'), all.x=TRUE)	
	#	get w per interval
	cat(paste('\nsetting likelihood to likelihood of pair / number transmission intervals'))
	set(YXc, NULL, 'score.Y.raw', YXc[, score.Y])
	set(YXc, NULL, 'score.Y', YXc[, score.Y*w.tn])				
	YXc		<- merge(YXc, YXc[, list(t=t, t.Patient=t.Patient, ntPatient= length(unique(t.Patient))), by=c('Patient')], by=c('Patient','t.Patient','t'))	
	YXc[, p.seqnc:= p.seq*p.nc]
	YXc[, me:= 1/(p.seq*p.nc)-1]	
	#	calculate absolute transmission probabilities
	#	need first: predicted score by stage and time t
	#	need second: predicted missing intervals per recipient, by stage and time t
	#	then calculate: exp missing score.Y per recipient by time t= sum_stage missexp.pr(stage,t) * scoreY.pr(stage,t)
	tmp		<- rltvtp.model.time.150923(YXc)
	rltvp	<- tmp$prd
	rltvd	<- tmp$fit
	rltvm	<- tmp$mo
	YXc		<- subset(YXc, AnyPos_T1<2010.8)
	tmp		<- miss.model.150923(YXc)
	mp		<- tmp$prd
	mma		<- tmp$moa
	mmb		<- tmp$mob
	md		<- tmp$fit
	#YXm		<- miss.smoothing.150909(YXc)
	#YXm		<- miss.smoothing.150914(YXc, mm, md, rltvm, rltvd, indir=NA, infile=NA, suffix='150914')
	#YXm		<- miss.smoothing.150916(YXc, mm, md, rltvp, rltva, indir=indir, infile=infile, suffix='150916')
	#YXm		<- miss.smoothing.150917(YXc, mm, md, rltvp, rltva, indir=indir, infile=infile, suffix='150917', use.obs.pseqnc=TRUE)
	#YXm		<- miss.smoothing.150918(YXc, mm, md, rltvp, rltvd, rltvm, indir=indir, infile=infile, suffix='150918', use.obs.pseqnc=TRUE)
	YXm		<- miss.smoothing.150923(YXc, mp, mma, mmb, md, rltvp, rltvd, rltvm, indir=indir, infile=infile, suffix='150922', use.obs.pseqnc=TRUE)
	YXm		<- merge(YXm, unique(subset(YXc, select=c(Patient, AnyPos_T1, t.period))), by='Patient')
	if(0)
	{
		#	get expected missing term for denominator 'scoreYm.per.rec'
		tmp		<- YXm[, list(miss.per.rec= sum(miss.prsm), scoreYm.per.rec= sum(miss.prsm*scoreY.pr)), by='Patient']
		YXc		<- merge( YXc, tmp, by='Patient' )
		#	then calculate p = w / ( sum(w) + scoreYm.per.rec)
		YXc		<- merge(YXc, YXc[, list(score.p= score.Y/(sum(score.Y)+scoreYm.per.rec), score.p.nadj= score.Y/sum(score.Y), t.Patient=t.Patient, t=t), by='Patient'], by=c('Patient','t.Patient','t'))
		#ggplot(YXc, aes(x=AnyPos_T1, y=miss.per.rec)) + geom_point()		
	}
	
	#prop.consistency(YXc, YXm)	
	tmp2	<- YXm[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.stAgeC','Patient')]	
	tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','t.stAgeC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','t.stAgeC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.stAgeC=t.stAgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.stAgeC'))	
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))
	YXpr	<- merge(YXpr, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	
	#	consistency	--> aggregating first by t.stAgeC leads to slightly diff results
	tmp		<- YXpr[, list(p.per.rec=sum(p.per.rec), po.per.rec=sum(p.per.rec*(1-pm.per.rec)), pm.per.rec=sum(p.per.rec*pm.per.rec)), by=c('Patient','stageC','t.period')]	
	tmp		<- tmp[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(po.per.rec), pm.per.tc=mean(pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodStageC',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=7)
	#	consistency: combine UAE and UAC		
	tmp		<- copy(YXpr)
	set(tmp, NULL, 'stageC', tmp[, gsub('UAE|UAC','UA', stageC)])
	tmp		<- tmp[, list(p.per.rec=sum(p.per.rec), po.per.rec=sum(p.per.rec*(1-pm.per.rec)), pm.per.rec=sum(p.per.rec*pm.per.rec)), by=c('Patient','stageC','t.period')]
	tmp		<- tmp[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(po.per.rec), pm.per.tc=mean(pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodStageC_UA',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=7)
	
	#	by t.stAgeC	
	tmp		<- YXpr[, list(stageC=stageC[1], t.AgeC=t.AgeC[1], p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.stAgeC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(t.AgeC~stageC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodtstAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=10)
	
	#	by t.AgeC tperiod
	tmp2	<- YXm[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.AgeC','Patient')]	
	tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','t.AgeC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','t.AgeC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.AgeC=t.AgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.AgeC'))	
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~t.AgeC) + theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_ByTperiodtAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=10, h=7)
	#	by t.AgeC tmid
	YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009,Inf), labels=c('<2007','2007','2008','>=2009'))]
	#YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009, 2010, Inf), labels=c('<2007','2007','2008','2009','2010'))]
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','tmidC')]	
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	tmp[, group:='probable transmitter']
	tmp2	<- Xmn[, list(to=sum(to), tm=sum(tm)), by=c('tmidC','t.AgeC')]
	tmp2	<- tmp2[, list(t.AgeC=t.AgeC, value= (to+tm)/sum(to+tm)), by='tmidC']
	tmp2[, variable:='observed']
	tmp2[, group:='population']
	tmp		<- rbind(tmp, tmp2, fill=TRUE, use.names=TRUE)
	set(tmp, NULL, 'group', tmp[, factor(group, labels=c('probable transmitter','population'), levels=c('probable transmitter','population'))])
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggp		<- ggplot(subset(tmp, group=='probable transmitter'), aes(x=t.AgeC)) + geom_bar(aes(y=100*value, fill=t.AgeC, alpha=variable), stat='identity',colour='black') +			
			scale_y_continuous(expand=c(0,0), limits=c(0,45)) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			scale_fill_brewer(palette='Set1', guide=FALSE) +
			facet_grid(~tmidC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) +
			labs(x='\nage', y='%', alpha='')
	z <- ggplot_gtable(ggplot_build(ggp))
	z <- gtable_add_rows(z, z$heights[[3]], 2)
	z <- gtable_add_grob(z, 
			list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
					textGrob("Transmissions by calendar year", gp = gpar(col = 'black'))), 3, 4, 3, 10, name = paste(runif(2)))
	z <- gtable_add_rows(z, unit(1/8, "line"), 3)
	grid.newpage()
	pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByTmid',suffix,'.pdf',sep=''),infile), sep=''), w=12, h=3.5)
	grid.draw(z)
	dev.off()
	ggp		<- ggplot(subset(tmp, group=='population'), aes(x=t.AgeC)) + geom_bar(aes(y=100*value, fill=t.AgeC, alpha=variable), stat='identity',colour='black') +			
			scale_y_continuous(expand=c(0,0), limits=c(0,45)) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			scale_fill_brewer(palette='Set1', guide=FALSE) +
			facet_grid(~tmidC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) +
			labs(x='\nage', y='%', alpha='')
	z <- ggplot_gtable(ggplot_build(ggp))
	z <- gtable_add_rows(z, z$heights[[3]], 2)
	z <- gtable_add_grob(z, 
			list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
					textGrob("Population by calendar year", gp = gpar(col = 'black'))), 3, 4, 3, 10, name = paste(runif(2)))
	z <- gtable_add_rows(z, unit(1/8, "line"), 3)
	grid.newpage()
	pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByTmidPop',suffix,'.pdf',sep=''),infile), sep=''), w=12, h=3.5)
	grid.draw(z)
	dev.off()
	
	
	#YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009,Inf), labels=c('<2007','2007','2008','>=2009'))]
	YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009, 2010, Inf), labels=c('<2007','2007','2008','2009','2010'))]
	tmp		<- subset(YXpr, tmidC=='2009')
	tmp		<- tmp[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','tmidC')]	
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	tmp2	<- data.table(t.AgeC=c('(-1,25]','(25,35]','(35,45]','(45,55]','(55,100]'), mid=c(0, 27.3, 31.8, 29.9, 11.0), lower=c(0, 25.5, 30.1, 28.4, 10.2), upper=c(0, 29, 33.5, 31.4, 11.8))
	ggplot(tmp, aes(x=t.AgeC)) + geom_bar(aes(y=100*value, fill=t.AgeC, alpha=variable), stat='identity',colour='black') +
			geom_errorbar(data=tmp2, aes(x=t.AgeC, ymax=upper, ymin=lower), size=1, width=0.3) +
			geom_point(data=tmp2, aes(x=t.AgeC, y=mid), size=3) +
			scale_y_continuous(expand=c(0,0), limits=c(0,40)) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			scale_fill_brewer(palette='Set1', guide=FALSE) +
			facet_grid(~tmidC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) +
			labs(x='\nage of probable transmitter', y='%', alpha='')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_Skarbinski',suffix,'.pdf',sep=''),infile), sep=''), w=4, h=6)
	
	#
	#	by p.AgeC
	#
	YXmp	<- merge(YXm, unique(subset(YXc, select=c(Patient, AgeC))), by='Patient')
	tmp2	<- YXmp[, list(scoreYm.per.recst=sum(miss.prsm*scoreY.pr)), by=c('t.period','t.AgeC','AgeC','Patient')]
	tmp2[, p.AgeC:= paste(as.character(t.AgeC),'->',as.character(AgeC),sep='')]
	tmp3	<- YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.period','p.AgeC','t.AgeC', 'AgeC')]
	tmp2	<- merge(tmp2, tmp3, by=c('Patient','AgeC','t.AgeC','p.AgeC','t.period'), all.x=1)
	set(tmp2, tmp2[, which(is.na(obs.per.recst))], c('scoreYo.per.recst','obs.per.recst'), 0)
	YXpr	<- tmp2
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.AgeC=t.AgeC, AgeC=AgeC, p.AgeC=p.AgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.AgeC','AgeC','p.AgeC'))	
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by=c('Patient','t.period'))
	#	p.AgeC totals
	require(gtable)
	tmp		<- YXpr[, list(p.per.tc=sum(p.per.rec), po.per.tc=sum(p.per.rec*(1-pm.per.rec)), pm.per.tc=sum(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))	
	ggp		<- ggplot(tmp, aes(x=t.AgeC, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~AgeC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) + 
			scale_y_continuous(expand=c(0,0), limits=c(0,70), breaks=seq(0,200,20), minor_breaks=seq(0,200,5)) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			scale_fill_brewer(palette='Set1', guide=FALSE) +
			labs(x='\nage of probable transmitter', y='#', alpha='')	
	z <- ggplot_gtable(ggplot_build(ggp))
	z <- gtable_add_rows(z, z$heights[[3]], 2)
	z <- gtable_add_grob(z, 
			list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
					textGrob("Transmissions by age of recipients in cohort", gp = gpar(col = 'black'))), 3, 4, 3, 12, name = paste(runif(2)))
	z <- gtable_add_rows(z, unit(1/8, "line"), 3)
	grid.newpage()
	pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByAgeCTotals',suffix,'.pdf',sep=''),infile), sep=''), w=15, h=5)
	grid.draw(z)
	dev.off()
	
	#	p.AgeC proportions
	require(gtable)
	tmp		<- YXpr[, list(p.per.tc=100*mean(p.per.rec), po.per.tc=100*mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=100*mean(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))	
	ggp		<- ggplot(tmp, aes(x=t.AgeC, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~AgeC) + theme_bw() + theme(legend.position='bottom', panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey75", size=0.4), panel.grid.minor.y=element_line(colour="grey75", size=0.4)) + 
			scale_y_continuous(expand=c(0,0), limits=c(0,55)) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			scale_fill_brewer(palette='Set1', guide=FALSE) +
			labs(x='\nage of probable transmitter', y='%', alpha='')	
	z <- ggplot_gtable(ggplot_build(ggp))
	z <- gtable_add_rows(z, z$heights[[3]], 2)
	z <- gtable_add_grob(z, 
			list(	rectGrob(gp = gpar(col='black', fill=gray(0.8))),
					textGrob("Transmissions by age of recipient", gp = gpar(col = 'black'))), 3, 4, 3, 12, name = paste(runif(2)))
	z <- gtable_add_rows(z, unit(1/8, "line"), 3)
	grid.newpage()
	pdf(file=paste(indir, '/', gsub('\\.R',paste('_ptp_tAgeCByAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=12, h=5)
	grid.draw(z)
	dev.off()
	
	#	p.AgeC total
	tmp		<- YXpr[, list(n.per.tc=sum(p.per.rec), no.per.tc=sum(p.per.rec*(1-pm.per.rec)), nm.per.tc=sum(p.per.rec*pm.per.rec)), by=c('t.AgeC','AgeC','p.AgeC')]
	tmp[, p.per.tc:= 100*n.per.tc/sum(n.per.tc)]
	tmp[, text:= as.character(round(100*n.per.tc/sum(n.per.tc),d=1))]
	tmp[, text.size:=4]
	set(tmp, tmp[, which(p.per.tc<3)], 'text.size', 2)
	tmp[, po.per.tc:= 100*no.per.tc/sum(n.per.tc)]
	tmp[, pm.per.tc:= 100*nm.per.tc/sum(n.per.tc)]	
	setnames(tmp, c('po.per.tc','p.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('expected missing','observed'), labels=c('expected missing','observed'))])
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('observed','expected missing'), labels=c('observed','expected missing'))])
	ggplot(tmp, aes(x=AgeC, y=t.AgeC)) + 			
			geom_point(aes(size=value, alpha=variable, fill=variable), pch=21) +
			geom_text(data=subset(tmp, text.size==4), aes(label=text), size=4, colour='black', hjust=0.5, vjust=0.5) +
			geom_text(data=subset(tmp, text.size==2), aes(label=text), size=2, colour='black', hjust=0.5, vjust=0.5) +
			scale_size_continuous(range=c(4,50), guide=FALSE) + 
			scale_fill_manual(values=c('expected missing'='#E41A1C', 'observed'='#E41A1C')) +
			scale_alpha_manual(values=c('expected missing'=0.3, 'observed'=1)) +
			theme_bw() + theme(legend.position='bottom') + labs(x='age of recipient', y='age of probable transmitter', fill='transmissions\n(%)', alpha='transmissions\n(%)')
	ggsave(file=paste(indir, '/', gsub('\\.R',paste('_ptp_AgeCtAgeC',suffix,'.pdf',sep=''),infile), sep=''), w=5, h=5)
}
######################################################################################
adjust.dev.Age_253045.Stage_UA_UC_D_T_F<- function()
{
	require(survival)
	YX
	X.tables
	t.recent.endctime
		
	indir		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C2V100bInfT7STRAT_m5A.R'
	
	YXo			<- copy(YX)	#before strat
	YXo2		<- stratificationmodel.Age_253045.Stage_UA_UC_D_T_F(YX)	
	risk.col	<- 't.stAgeC'	
	#
	#	focus on >2003 and before ART start
	#
	#YX			<- subset(YX, t>2004.5 & !stageC%in%c('L','T'))
	YX			<- copy(YXo2)
	YX			<- subset(YX, AnyPos_T1>2005)
	tmp			<- YX[, which(stageC=='L')]
	set(YX, tmp, 't.AgeC', '(-1,100]')
	set(YX, tmp, 't.stAgeC', 'L_(-1,100]')
	tmp			<- YX[, which(stageC=='T')]
	set(YX, tmp, 't.AgeC', '(-1,100]')
	set(YX, tmp, 't.stAgeC', 'T_(-1,100]')	
	#tmp		<- YXc[, which(stageC=='L' | stageC=='T')]
	#set(YXc, tmp, 't.AgeC', '(-1,100]')
	#set(YXc, tmp, 'stageC', 'LT')
	#set(YXc, tmp, 't.stAgeC', 'LT_(-1,100]')	
	set(YX, NULL, 'stageC', YX[, factor(as.character(stageC))])
	set(YX, NULL, 't.stAgeC', YX[, factor(as.character(t.stAgeC))])
	set(YX, NULL, 't.stAgeC.prd', YX[, factor(as.character(t.stAgeC.prd))])
	
	#	get clustering, sampling and non-censoring probabilities 
	#	for each prob transmitter at the transmission interval 
	sm			<- X.tables$sm
	YXs			<- merge(YX, sm, by='t.Patient')
	#	p.nc	(non-censoring prob)
	cp			<- copy(X.tables$cm$cens.p)		#model built after subsampling, so cannot apply straight-away
	set(cp, NULL, c('BS','p.nc'), NULL)
	cm			<- X.tables$cm$cens.m	
	tp.df		<- subset(YXs, grepl('^U', YXs[[risk.col]]))[, list(tm=mean(t)), by=c('t.Patient','Patient')]
	tp.df		<- censoring.get.dataset(tp.df, df.all.allmsm)
	setnames(YXs, c('t.Patient','Patient'), c('Patient','r.Patient'))
	#	YXc so far only contains transmitters with at least one undiagnosed interval
	YXc			<- merge(YXs, subset(tp.df, select=c(Patient, r.Patient, tm, AACD4C)), by=c('Patient','r.Patient'))
	set(YXc, NULL, 'tm', YXc[, tm-t.recent.endctime])	
	YXc[, p.nc:=predict(cm, data=cp, newdata=subset(YXc, select=c(tm, AACD4C)), type='response')]
	#	p.nc is allocated to transmitters with at least one undiagnosed interval
	#	set p.nc to 1 for all intervals after diagnosis
	set(YXc, YXc[, which(!is.na(p.nc) & !grepl('^U', stageC))], 'p.nc', 1)
	#ggplot(YXc, aes(x=tm, y=p.nc, group=AACD4C)) + geom_line() + facet_wrap(~AACD4C, ncol=4)
	#ggplot(YXc, aes(x=tm)) + geom_histogram()	
	setnames(YXs, c('Patient','r.Patient'), c('t.Patient','Patient'))
	setnames(YXc, c('Patient','r.Patient'), c('t.Patient','Patient'))
	YXc			<- merge(YXs, subset(YXc, select=c(t.Patient, Patient, t, p.nc)), by=c('t.Patient','Patient','t'), all.x=TRUE)
	set(YXc, YXc[, which(is.na(p.nc))], 'p.nc', 1)
	#	get p.clu and p.pt
	risk.table	<- X.tables$st$risk.table
	tmp			<- clustering.getrisks(risk.table, indir=indir, infile=infile, rskf.name='t.stAgeC', rskf.baseline='D_(30,45]')
	tpt			<- tmp$tpt
	tclu		<- tmp$tclu	
	tmp			<- subset(tclu, t.period=='Overall', select=c(factor2, p.clu))
	setnames(tmp, 'factor2', 't.stAgeC')	
	set(tmp, NULL, 't.stAgeC', tmp[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)], labels=YXc[, levels(t.stAgeC)])])
	YXc			<- merge(YXc, tmp, by=c('t.stAgeC'), all.x=TRUE)
	tmp			<- subset(tpt, t.period=='Overall', select=c(factor2, p.pt))
	setnames(tmp, 'factor2', 't.stAgeC')		
	set(tmp, NULL, 't.stAgeC', tmp[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)], labels=YXc[, levels(t.stAgeC)])])
	YXc			<- merge(YXc, tmp, by=c('t.stAgeC'), all.x=TRUE)	
	#	get w per interval
	cat(paste('\nsetting likelihood to likelihood of pair / number transmission intervals'))
	set(YXc, NULL, 'score.Y.raw', YXc[, score.Y])
	set(YXc, NULL, 'score.Y', YXc[, score.Y*w.tn])				
	YXc		<- merge(YXc, YXc[, list(t=t, t.Patient=t.Patient, ntPatient= length(unique(t.Patient))), by=c('Patient')], by=c('Patient','t.Patient','t'))
	if(0)
	{
		#	are the w's confounded by cluster size?
		rltvtp.confounded.ntransmitters(YXc, indir, infile)
		#	no, of course not!! the p's will be confounded by cluster size!!		
		#	mean or median w for each stage:	
		rltvtp.exploredistribution(YXc, indir, infile)		
	}
	YXc[, p.seqnc:= p.seq*p.nc]
	YXc[, me:= 1/(p.seq*p.nc)-1]	
	#	calculate absolute transmission probabilities
	#	need first: predicted score by stage and time t
	#	need second: predicted missing intervals per recipient, by stage and time t
	#	then calculate: exp missing score.Y per recipient by time t= sum_stage missexp.pr(stage,t) * scoreY.pr(stage,t)
	tmp		<- rltvtp.model.time.150901(YXc)
	rltvm	<- tmp$rltvm
	rltvd	<- tmp$rltvd
	YXc		<- subset(YXc, AnyPos_T1<2010.8)
	tmp		<- miss.model.150904b(YXc)	
	md		<- tmp$miss.d
	mm		<- tmp$miss.m
	YXm		<- miss.smoothing.150909(YXc)
	#	get expected missing term for denominator 'scoreYm.per.rec'
	tmp		<- YXm[, list(miss.per.rec= sum(miss.prsm), scoreYm.per.rec= sum(miss.prsm*scoreY.pr)), by='Patient']
	YXc		<- merge( YXc, tmp, by='Patient' )
	#	then calculate p = w / ( sum(w) + scoreYm.per.rec)
	YXc		<- merge(YXc, YXc[, list(score.p= score.Y/(sum(score.Y)+scoreYm.per.rec), score.p.nadj= score.Y/sum(score.Y), t.Patient=t.Patient, t=t), by='Patient'], by=c('Patient','t.Patient','t'))
	ggplot(YXc, aes(x=AnyPos_T1, y=miss.per.rec)) + geom_point()
	
	#	get expected missing by stage for each transmitter
	tmp		<- YXm[, list(miss.per.recst= sum(miss.prsm), scoreYm.per.recst= sum(miss.prsm*scoreY.pr)), by= c('Patient','t.stAgeC','stageC','t.AgeC')]
	#	merge with observed scores by stage for each transmitter	
	#	DEPENDS ON STRATIFICATION - explore first on t.stAgeC
	YXpr	<- merge(YXc[, list(obs.per.recst=length(t), scoreYo.per.recst=sum(score.Y)), by=c('Patient','t.stAgeC')], tmp, by=c('Patient','t.stAgeC'), all=1)
	set(YXpr, YXpr[, which(is.na(obs.per.recst))], c('obs.per.recst','scoreYo.per.recst'), 0)
	YXpr[, scoreY.per.recst:= scoreYo.per.recst+scoreYm.per.recst]
	tmp		<- YXpr[, list(t.stAgeC=t.stAgeC, p.per.rec= scoreY.per.recst/sum(scoreY.per.recst), pm.per.rec=scoreYm.per.recst/scoreY.per.recst), by='Patient']
	YXpr	<- merge(YXpr, tmp, by=c('Patient','t.stAgeC'))
	YXpr	<- merge(YXpr, YXc[, list(tmid= mean(t), AnyPos_T1=AnyPos_T1[1]),  by=c('Patient','t.period')], by='Patient')
	
	#	get rolling mean and plot
	setkey(YXpr, t.stAgeC, tmid)	
	tmp		<- YXpr[, {
				z	<- unique(tmid)
				list(tmid=z, prm.per.rec= sapply(seq_along(z), function(i)	mean(p.per.rec[ which(abs(z[i]-tmid)<=1) ]))	)	
			}, by='t.stAgeC']
	YXpr	<- merge(YXpr, tmp, by=c('tmid','t.stAgeC'))
	ggplot(YXpr, aes(x=tmid)) + geom_point(aes(y=p.per.rec, colour=t.AgeC)) + geom_line(aes(y=prm.per.rec), colour='black') + 
			facet_grid(t.AgeC~stageC) + theme_bw() 
	
	
	#	plot proportion of obs+missing across time
	YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2008, 2009,Inf), labels=c('<2007','2008','2009','>2009'))]
	YXpr[, tmidC:= cut(tmid, breaks=c(-Inf, 2007, 2009,Inf), labels=c('<2007','2008/09','>2009'))]
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('t.stAgeC','tmidC')]
	tmp		<- merge(tmp, unique(subset(YXpr, select=c(t.stAgeC, t.AgeC, stageC))), by='t.stAgeC')
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=tmidC, y=value, fill=t.AgeC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(t.AgeC~stageC) + theme_bw() + theme(legend.position='bottom')
	
	#	consistency with prev analysis
	tmp		<- YXpr[, list(p.per.tc=mean(p.per.rec), po.per.tc=mean(p.per.rec*(1-pm.per.rec)), pm.per.tc=mean(p.per.rec*pm.per.rec)), by=c('stageC','t.period')]
	setnames(tmp, c('po.per.tc','pm.per.tc'), c('observed','expected missing'))
	tmp		<- melt(tmp, measure.vars=c('observed','expected missing'))
	#tmp[,sum(p.per.tc),by='tmidC']	
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC, alpha=variable)) + geom_bar(stat='identity',colour='black') + 
			facet_grid(~stageC) + theme_bw() + theme(legend.position='bottom')
	
	
	
	#	are the p's confounded by cluster size? --> Yes, quite strongly!
	altvtp.confounded.ntransmitters(YXc, indir, infile, suffix='150909')	
	#	explore p across time and stages --> I am not sure, it s not super convincing ..
	#	TODO: get CIs on Gamma mu
	altvtp.explore.time.models.150909(YXc, indir, infile, suffix='150909')
	#	for time-independent analysis: find best depdendence on ntPatient
	altvtp.explore.distribution(YXc, indir, infile)
	#	GANtSLI is best, followed by GANtSL

	#	get RR irrespective of time 


}
######################################################################################
censoring.dev.sampling.model.calculate<- function()
{
	#	take Xmsm potential transmitters to recipients and collect variables
	#	t.isAcute, t.diag time, t.age at diag, t.RegionHospital
	#	
	#	for pseudo censoring time tc, and censor every transmitter with t.AnyPos_T1>tc -> y=0
	#	p is then the prob of not being censored
	#	p depends on as a spline on trintervaltime t. Too much data. Take midpoint of infection window.
	#	spline differs for isAcute. 
	#	Age of transmitter at diagnosis? If old, likely to have progressed further and more likely to be observed. 
	#	CD4 of transmitter at diagnosis. If low CD4, more likely to be oserved.
	#	RegionHospital	
	indir		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C2V100bInfT7STRAT_m5B.R'
	load(paste(indir,'/',infile,sep=''))
	stage.labels<- c("UAC", "UAE", "UC", "DAC", "D", "L", "T")
	
	#need df.all.allmsm
	c.period	<- 0.125
	c.smpl.n	<- 50
	bs.n=10
	bs.cdelta.min=2
	bs.cdelta.max=3
	t.recent.endctime
		
	cens.bs		<- runif(bs.n, bs.cdelta.min, bs.cdelta.max)
	cens.bs		<- data.table(cens.t=t.recent.endctime-cens.bs, cens.delta=cens.bs, BS=seq_along(cens.bs))
	bs			<- 1
	
	cens.m		<- vector('list', length=bs.n)
	cens.p		<- vector('list', length=bs.n)
	require(zoo)
	tp.df		<- subset(X.msm, grepl('^U',stageC))[, list(tm=mean(t)), by=c('t.Patient','Patient')]
	tpd.df		<- censoring.get.dataset(tp.df, df.all.allmsm)
	tpd.df[, tmo:= tpd.df[, tm]]
	#	add pseudo censoring time
	for(bs in seq_len(bs.n))
	{
		tpd.df[, NCNS:= tpd.df[, as.numeric(AnyPos_T1<cens.bs[bs, cens.t])]]	#transmitter not censored
		set(tpd.df, NULL, 'tm', tpd.df[, tmo-cens.bs[bs, cens.t]])
		#	balanced sub sampling so regression is comp feasible
		tpds.df		<- censoring.subsample.tmC_ACD4C_AgeC_T1(tpd.df, t.recent.endctime, bs.cdelta.min, bs.cdelta.max, c.period)
		tmp			<- censoring.model.150728(tpds.df)
		cens.m[[bs]]<- tmp$model
		cens.p[[bs]]<- copy(tmp$predict)
		cens.p[[bs]][, BS:=bs]
	}
	cens.p		<- do.call('rbind', cens.p)
	cens.pm		<- cens.p[, list(p.nc.med= median(p.nc)), by=c('Patient','r.Patient')]
	cens.p		<- merge(cens.p, subset(tpd.df, select=c(Patient, r.Patient, AgeC_T1, ACD4C, AACD4C)), by=c('Patient','r.Patient'))
	ggplot(cens.p, aes(x=tm, y=p.nc, group=factor(BS), colour=factor(BS))) + geom_line() + facet_grid(AgeC_T1~ACD4C)	
	ggsave(file=paste(outdir, '/', outfile, '_CENSMODEL_cs4_bsn10.pdf', sep=''), w=20, h=15)	
	#	save
	save(tp.df, tpd.df, cens.p, cens.m, df.all.allmsm, file=paste(outdir, '/', outfile, '_CENSMODEL_cs4_bsn10.R', sep=''))
	
	#	predict
	z		<- subset(cens.p, BS==1 & AACD4C=='(-Inf,250]-(30,45]')[1,]
	z		<- merge(subset(z, select=which(colnames(z)!='tm')), data.table(Patient=z$Patient[1], tm=seq(-5,0,0.25)), by='Patient')
	predict(cens.m[[4]], data=subset(cens.p, BS==4), newdata=z , type='response')
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
	indir			<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	infile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3pa1H1.48C2V100bInfT7_StablesSEQ_m5A.R'	
	stage.labels	<- c('UA','UC','D','T','L')
	pdf.w			<- 15
	infile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3pa1H1.48C2V100bInfT7_StablesSEQ_m5B.R'
	stage.labels	<- c("UAC", "UAE", "UC", "DAC", "D", "L", "T")
	pdf.w			<- 20
	
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
	set(dfp, NULL, 'stage', dfp[, factor(stage, levels=stage.labels)])
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
	ggsave(file=paste(indir, '/', gsub('\\.R','_STRATSAMPLING_cmp_Age.pdf',infile), sep=''), w=pdf.w, h=5)
	
	ggplot(dfp, aes(x=p.seq, y=t.period.long, colour=stage)) + geom_point(size=3) +
			scale_x_continuous(breaks=seq(0,100,10), minor_breaks=NULL, limit=c(0,100), expand=c(0,0)) +			
			#scale_colour_manual(values=dfp[, unique(factor.color)], guide=FALSE) +
			theme_bw() +				
			theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=10), legend.position='bottom', panel.grid.major.x=element_line(colour="grey70", size=0.4), panel.grid.minor.x=element_line(colour="grey70", size=0.4), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), panel.margin = unit(2, "lines")) + 				
			labs(y='', x='overlap intervals\nof a potential transmitter with a sequence\n(%)', colour='stage of potential transmitter\nat transmission interval') +
			facet_grid(~age)  				
	ggsave(file=paste(indir, '/', gsub('\\.R','_STRATSAMPLING_cmp_Stage.pdf',infile), sep=''), w=15, h=5)	
	#
	#	compare with new sampling model
	#
	tmp		<- subset(dfp, select=c(t.period.long, age, stage, p.seq))
	setnames(tmp, c('p.seq','stage','age'),c('p.seq.prev','stageC','t.AgeC'))
	z		<- YXs[, list(p.seq.me=mean(p.seq), p.seq.md=median(p.seq)), by=c('t.period.long','t.AgeC','stageC')]
	tmp		<- merge(tmp, z, by=c('t.period.long','t.AgeC','stageC'))
	ggplot(tmp, aes(x=p.seq.prev, y=100*p.seq.me, colour=stageC)) + geom_abline(intercept=0, slope=1, colour='black') + 
			geom_point() +
			facet_grid(t.period.long~t.AgeC) +
			labs(colour='stage of potential transmitter\nat transmission interval') +
			theme_bw() + theme(legend.position='bottom')
	ggsave(file=paste(indir, '/', gsub('\\.R','_STRATSAMPLINGMODEL_cmp_Stage.pdf',infile), sep=''), w=7, h=7)
	#
	#	when there s few cases eg (-1,25] then the sampling probs may disagree 
	#	overall, the new model gives higher sampling probabilities
}
######################################################################################
sampling.model.calculate<- function(X.msm, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, df.treatment.allmsm, t.period=NA, t.endctime=NA, method.lRNA.supp=NA, method.lRNA.nsupp=4, contact.grace=NA, resume=TRUE, save.file=NA)	
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
		if(is.null(X.msm))
			return(NULL)
		stopifnot(!is.na(t.period), !is.na(t.endctime), !is.na(method.lRNA.supp), !is.na(contact.grace))
		pt.df		<- data.table(Patient=X.msm[, unique(t.Patient)])
		ptd.df		<- sampling.get.dataset(pt.df, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, df.treatment.allmsm, t.period, t.endctime, method.lRNA.supp=method.lRNA.supp, method.lRNA.nsupp=4, contact.grace=contact.grace)
		ptr.df		<- sampling.model.150722(ptd.df)
		setnames(ptr.df, 'Patient', 't.Patient')
		ptr.df		<- subset( ptr.df, select=c(t.Patient,p.seq) )
		if(!is.na(save.file))
			save(ptr.df, file=save.file)
	}
	ptr.df
}
######################################################################################
miss.model.150904<- function(YXc, indir, infile)
{	
	YXm		<- subset(YXc, AnyPos_T1<2010.8, select=c(t, t.stAgeC, me))
	set(YXm, NULL, 't.stAgeC', YXm[, factor(as.character(t.stAgeC))])
	mm4		<- gamlss(me~t.stAgeC+ns(t, df=3):t.stAgeC, sigma.formula=~t.stAgeC, data=YXm, family=GA())
	list(miss.m=mm4, miss.d=YXm)	
}
######################################################################################
miss.model.150904b<- function(YXc, indir, infile)
{	
	YXm		<- subset(YXc, AnyPos_T1<2010.8, select=c(t, t.stAgeC, p.seqnc))
	set(YXm, NULL, 't.stAgeC', YXm[, factor(as.character(t.stAgeC))])
	mm4		<- gamlss(p.seqnc~t.stAgeC+ns(t, df=3):t.stAgeC, sigma.formula=~t.stAgeC, data=YXm, family=GA())
	list(miss.m=mm4, miss.d=YXm)	
}
######################################################################################
miss.model.150917<- function(YXc, indir, infile)
{	
	YXm		<- subset(YXc, AnyPos_T1<2010.8, select=c(t, t.stAgeC, p.seqnc))
	set(YXm, NULL, 't.stAgeC', YXm[, factor(as.character(t.stAgeC))])
	mm4		<- gamlss(p.seqnc~t.stAgeC+ns(t, df=4):t.stAgeC, sigma.formula=~t.stAgeC, data=YXm, family=GA())
	list(miss.m=mm4, miss.d=YXm)	
}
######################################################################################
miss.model.150923Gen<- function(YXc)
{		
	YXm		<- subset(YXc, AnyPos_T1<2010.8, select=c(t, stageC, t.AgeC, t.stAgeC, p.seqnc))
	YXm[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXm[, unique(gsub('TO_.*','TO_1920-95',gsub('TS_.*','TS_1920-95',gsub('L_.*','L_1920-95',gsub('1956-65','1920-65',gsub('1920-55','1920-65',levels(t.stAgeC),fixed=1),fixed=1)))))]
	set(YXm, NULL, 't.stAgeC2', YXm[, factor(gsub('TO_.*','TO_1920-95',gsub('TS_.*','TS_1920-95',gsub('L_.*','L_1920-95',gsub('1956-65','1920-65',gsub('1920-55','1920-65',t.stAgeC2,fixed=1),fixed=1)))), levels=tmp, labels=tmp)])
	mm5a	<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=3):t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')), family=GA())
	mm5b	<- gamlss(p.seqnc~t.stAgeC2+t:t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, stageC%in%c('UAC','TO','TS','L')), family=GA())
	
	prd<- function(dat, mm5a, mm5b, YXm)
	{		
		dat[, t.stAgeC2:= t.stAgeC]
		tmp		<- dat[, unique(gsub('TO_.*','TO_1920-95',gsub('TS_.*','TS_1920-95',gsub('L_.*','L_1920-95',gsub('1956-65','1920-65',gsub('1920-55','1920-65',levels(t.stAgeC),fixed=1),fixed=1)))))]
		set(dat, NULL, 't.stAgeC2', dat[, factor(gsub('TO_.*','TO_1920-95',gsub('TS_.*','TS_1920-95',gsub('L_.*','L_1920-95',gsub('1956-65','1920-65',gsub('1920-55','1920-65',t.stAgeC2,fixed=1),fixed=1)))), levels=tmp, labels=tmp)])
		me5a	<- predict(mm5a,data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')),newdata=subset(dat,select=c(t, stageC, t.AgeC, t.stAgeC, t.stAgeC2)),type='response',what='mu')
		me5b	<- predict(mm5b,data=subset(YXm, stageC%in%c('UAC','TO','TS','L')),newdata=subset(dat,select=c(t, stageC, t.AgeC, t.stAgeC, t.stAgeC2)),type='response',what='mu')
		dat[, pseqnc.pr:= me5a]
		tmp		<- dat[, which(stageC%in%c('UAC','TO','TS','L'))]
		set(dat, tmp, 'pseqnc.pr', me5b[tmp])
		dat				
	}
	list(prd=prd, moa=mm5a, mob=mm5b, fit=YXm)		
}
######################################################################################
miss.model.150922Gen<- function(YXc)
{		
	YXm		<- subset(YXc, AnyPos_T1<2010.8, select=c(t, stageC, t.AgeC, t.stAgeC, p.seqnc))
	YXm[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXm[, unique(gsub('1956-65','1920-65',gsub('1920-55','1920-65',levels(t.stAgeC),fixed=1),fixed=1))]	
	set(YXm, NULL, 't.stAgeC2', YXm[, factor(gsub('1956-65','1920-65',gsub('1920-55','1920-65',t.stAgeC2,fixed=1),fixed=1), levels=tmp, labels=tmp)])
	mm5a	<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=3):t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')), family=GA())
	mm5b	<- gamlss(p.seqnc~t.stAgeC2+t:t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, stageC%in%c('UAC','TO','TS','L')), family=GA())
	
	prd<- function(dat, mm5a, mm5b, YXm)
	{		
		dat[, t.stAgeC2:= t.stAgeC]
		tmp		<- dat[, unique(gsub('1956-65','1920-65',gsub('1920-55','1920-65',levels(t.stAgeC),fixed=1),fixed=1))]	
		set(dat, NULL, 't.stAgeC2', dat[, factor(gsub('1956-65','1920-65',gsub('1920-55','1920-65',t.stAgeC2,fixed=1),fixed=1), levels=tmp, labels=tmp)])		
		me5a	<- predict(mm5a,data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')),newdata=subset(dat,select=c(t, stageC, t.AgeC, t.stAgeC, t.stAgeC2)),type='response',what='mu')
		me5b	<- predict(mm5b,data=subset(YXm, stageC%in%c('UAC','TO','TS','L')),newdata=subset(dat,select=c(t, stageC, t.AgeC, t.stAgeC, t.stAgeC2)),type='response',what='mu')
		dat[, pseqnc.pr:= me5a]
		tmp		<- dat[, which(stageC%in%c('UAC','TO','TS','L'))]
		set(dat, tmp, 'pseqnc.pr', me5b[tmp])
		dat				
	}
	list(prd=prd, moa=mm5a, mob=mm5b, fit=YXm)		
}
######################################################################################
miss.model.150923<- function(YXc)
{		
	YXm		<- subset(YXc, AnyPos_T1<2010.8, select=c(t, stageC, t.AgeC, t.stAgeC, p.seqnc))
	YXm[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXm[, unique(gsub('TO_.*','TO_(-1,00]',gsub('TS_.*','TS_(-1,00]',gsub('L_.*','L_(-1,00]',gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',levels(t.stAgeC),fixed=1),fixed=1)))))]	
	set(YXm, NULL, 't.stAgeC2', YXm[, factor(gsub('TO_.*','TO_(-1,00]',gsub('TS_.*','TS_(-1,00]',gsub('L_.*','L_(-1,00]',gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',t.stAgeC2,fixed=1),fixed=1)))), levels=tmp, labels=tmp)])
	mm5a	<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=4):t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')), family=GA())
	mm5b	<- gamlss(p.seqnc~t.stAgeC2+t:t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, stageC%in%c('UAC','TO','TS','L')), family=GA())
	
	prd<- function(dat, mm5a, mm5b, YXm)
	{		
		dat[, t.stAgeC2:= t.stAgeC]
		tmp		<- dat[, unique(gsub('TO_.*','TO_(-1,00]',gsub('TS_.*','TS_(-1,00]',gsub('L_.*','L_(-1,00]',gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',levels(t.stAgeC),fixed=1),fixed=1)))))]	
		set(dat, NULL, 't.stAgeC2', dat[, factor(gsub('TO_.*','TO_(-1,00]',gsub('TS_.*','TS_(-1,00]',gsub('L_.*','L_(-1,00]',gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',t.stAgeC2,fixed=1),fixed=1)))), levels=tmp, labels=tmp)])
		me5a	<- predict(mm5a,data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')),newdata=subset(dat,select=c(t, stageC, t.AgeC, t.stAgeC, t.stAgeC2)),type='response',what='mu')
		me5b	<- predict(mm5b,data=subset(YXm, stageC%in%c('UAC','TO','TS','L')),newdata=subset(dat,select=c(t, stageC, t.AgeC, t.stAgeC, t.stAgeC2)),type='response',what='mu')
		dat[, pseqnc.pr:= me5a]
		tmp		<- dat[, which(stageC%in%c('UAC','TO','TS','L'))]
		set(dat, tmp, 'pseqnc.pr', me5b[tmp])
		dat				
	}
	list(prd=prd, moa=mm5a, mob=mm5b, fit=YXm)		
}
######################################################################################
miss.model.150922<- function(YXc)
{		
	YXm		<- subset(YXc, AnyPos_T1<2010.8, select=c(t, stageC, t.AgeC, t.stAgeC, p.seqnc))
	YXm[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXm[, unique(gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',levels(t.stAgeC),fixed=1),fixed=1))]
	set(YXm, NULL, 't.stAgeC2', YXm[, factor(gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',t.stAgeC2,fixed=1),fixed=1), levels=tmp, labels=tmp)])	
	mm5a	<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=4):t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')), family=GA())
	mm5b	<- gamlss(p.seqnc~t.stAgeC2+t:t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, stageC%in%c('UAC','TO','TS','L')), family=GA())
	
	prd<- function(dat, mm5a, mm5b, YXm)
	{		
		dat[, t.stAgeC2:= t.stAgeC]
		tmp		<- dat[, unique(gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',levels(t.stAgeC),fixed=1),fixed=1))]
		set(dat, NULL, 't.stAgeC2', dat[, factor(gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',t.stAgeC2,fixed=1),fixed=1), levels=tmp, labels=tmp)])

		me5a	<- predict(mm5a,data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')),newdata=subset(dat,select=c(t, stageC, t.AgeC, t.stAgeC, t.stAgeC2)),type='response',what='mu')
		me5b	<- predict(mm5b,data=subset(YXm, stageC%in%c('UAC','TO','TS','L')),newdata=subset(dat,select=c(t, stageC, t.AgeC, t.stAgeC, t.stAgeC2)),type='response',what='mu')
		dat[, pseqnc.pr:= me5a]
		tmp		<- dat[, which(stageC%in%c('UAC','TO','TS','L'))]
		set(dat, tmp, 'pseqnc.pr', me5b[tmp])
		dat				
	}
	list(prd=prd, moa=mm5a, mob=mm5b, fit=YXm)		
}
######################################################################################
miss.model.150831<- function(YXc)
{
	YXm		<- subset(YXc, AnyPos_T1<2010.8)
	tmp		<- YXm[, list(ttlme= sum(1/(p.seq*p.nc)-1)/length(unique(Patient)), nR=length(unique(Patient)) ), by='t']
	YXm		<- merge(YXm, tmp,by='t')
	tmp		<- YXm[, list(me= sum(1/(p.seq*p.nc)-1)/nR[1] ), by=c('t','t.stAgeC')]
	YXm		<- merge(YXm, tmp, by=c('t','t.stAgeC'))
	YXm[, nR:=NULL]
	YXm		<- subset(YXm, select=c(t, t.stAgeC, me))
	mm3		<- gamlss(me~t.stAgeC+ns(t, df=4):t.stAgeC, sigma.formula=~t.stAgeC, data=YXm, family=GA())
	list(miss.m=mm3, miss.d=YXm)
}
######################################################################################
miss.smoothing.150909.onns<- function(YXc, indir=NA, infile=NA)
{
	#	get scaffold: for all recipients at time t, all stages
	YXm		<- YXc[, list(nr=length(Patient)), by='t']
	YXm		<- merge(YXm, YXc[, list(tn=length(t.Patient)), by=c('t','Patient')], by='t')
	YXm		<- merge( as.data.table(expand.grid(t.stAgeC=YXc[, levels(t.stAgeC)], t=YXc[, sort(unique(t))])), YXm, by='t', all.x=TRUE, allow.cartesian=TRUE )
	#	add number observed intervals by t and stage
	YXm		<- merge(YXm, YXc[, list(ns=length(t.Patient)), by=c('t','Patient','t.stAgeC')], by=c('t','Patient','t.stAgeC'), all.x=TRUE)
	set(YXm, YXm[, which(is.na(ns))], c('ns'), 0)	
	setkey(YXm, t, Patient)	
	#	get smoothed w's and p's
	YXm[, pseqnc.pr:=predict(mm, data=md, newdata=subset(YXm, select=c(t, t.stAgeC)), type='response', what='mu')]
	YXm[, scoreY.pr:=predict(rltvm, data=rltvd, newdata=subset(YXm, select=c(t, t.stAgeC)), type='response', what='mu')]
	set(YXm, YXm[, which(scoreY.pr<0)], 'scoreY.pr', 0)
	set(YXm, YXm[, which(pseqnc.pr<0)], 'pseqnc.pr', 0)
	set(YXm, YXm[, which(scoreY.pr>1)], 'scoreY.pr', 1)
	set(YXm, YXm[, which(pseqnc.pr>1)], 'pseqnc.pr', 1)		
	#	define adaptive chunks with at least 5 miss.pr>0
	YXm		<- merge(YXm, YXm[, list(ns.per.t=sum(ns)), by=c('t','t.stAgeC')], by=c('t','t.stAgeC'))
	#subset(YXm, t.stAgeC=='UC_(-1,25]' & t>2010.3)[,
	tmp		<- YXm[, 
			{
				z	<- sapply(seq_along(t), function(i)	sort(abs(t[i]-t[ns>0]))		)				
				z	<- z[ ifelse(nrow(z)<20, nrow(z), 20), ]
				#print(z)
				z	<- sapply(seq_along(t), function(i) c(	'tn.per.ch'=sum(tn[ abs(t[i]-t)<=z[i] ]), 
									'ns.per.ch'=sum(ns[ abs(t[i]-t)<=z[i] ])				))
				#print(z)				
				list(t=t, Patient=Patient, tn.per.ch=z['tn.per.ch',], ns.per.ch=z['ns.per.ch',])				
			}, by='t.stAgeC']
	#	smooth ns
	YXm		<- merge(YXm, tmp, by=c('t','t.stAgeC','Patient'))
	YXm[, ns.sm:= ns.per.ch/tn.per.ch]
	#	predict missing based on the smooth n's 
	YXm[, miss.prsm:= ns.sm/pseqnc.pr-ns.sm]
	YXm		<- merge(YXm, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	#	plot smooth
	tmp		<- subset(YXm, select=c(t, stageC, t.AgeC, t.stAgeC, miss.prsm, nr))
	setkey(tmp, t, t.stAgeC)
	tmp		<- unique(tmp)
	if(!is.na(indir) & !is.na(infile))
	{
		ggplot( tmp, aes(x=t, y=miss.prsm, colour=t.AgeC, group=t.stAgeC)) + geom_step() + geom_point() +
				theme_bw() + theme(legend.position='bottom') +
				facet_wrap(~stageC, scales='free', ncol=5)
		file	<- paste(indir,'/',gsub('\\.R','_missexpsmoothBytstAgeC.pdf',infile),sep='')	
		ggsave(file=file, w=15,h=7)
		
		ggplot( YXm, aes(x=t, y=miss.prsm/nr, fill=t.AgeC)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient') +
				facet_grid(t.AgeC~stageC, scales='free')
		file	<- paste(indir,'/',gsub('\\.R','_missexpsmoothtotalBytstAgeC.pdf',infile),sep='')	
		ggsave(file=file, w=12,h=12)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient')
		file	<- paste(indir,'/',gsub('\\.R','_missexpsmoothtotal.pdf',infile),sep='')	
		ggsave(file=file, w=5,h=4)
	}
	YXm
}
######################################################################################
miss.smoothing.150914<- function(YXc, mm, md, rltvm, rltvd, indir=NA, infile=NA, suffix='')
{
	#	get scaffold: for all recipients at time t, all stages
	YXm		<- YXc[, list(nr=length(Patient)), by='t']
	YXm		<- merge(YXm, YXc[, list(tn=length(t.Patient)), by=c('t','Patient')], by='t')
	YXm		<- merge( as.data.table(expand.grid(t.stAgeC=YXc[, levels(t.stAgeC)], t=YXc[, sort(unique(t))])), YXm, by='t', all.x=TRUE, allow.cartesian=TRUE )
	#	add number observed intervals by t and stage
	YXm		<- merge(YXm, YXc[, list(ns=length(t.Patient)), by=c('t','Patient','t.stAgeC')], by=c('t','Patient','t.stAgeC'), all.x=TRUE)
	set(YXm, YXm[, which(is.na(ns))], c('ns'), 0)	
	setkey(YXm, t, Patient)	
	#	get smoothed w's and p's
	YXm[, pseqnc.pr:=predict(mm, data=md, newdata=subset(YXm, select=c(t, t.stAgeC)), type='response', what='mu')]
	YXm[, scoreY.pr:=predict(rltvm, data=rltvd, newdata=subset(YXm, select=c(t, t.stAgeC)), type='response', what='mu')]
	set(YXm, YXm[, which(scoreY.pr<0)], 'scoreY.pr', 0)
	set(YXm, YXm[, which(pseqnc.pr<0)], 'pseqnc.pr', 0)	
	#	if we have an observed sampling prob, use that
	tmp		<- YXc[, list(p.seqnc=median(p.seqnc)), by=c('Patient','t.stAgeC','t')]
	YXm		<- merge(YXm, tmp, by=c('Patient','t.stAgeC','t'), all.x=1)
	tmp		<- YXm[, which(is.na(p.seqnc))]
	set(YXm, tmp, 'p.seqnc', YXm[tmp, pseqnc.pr])
	#	predict missing based on the n's -- allow for large tails in NB
	YXm[, miss.pr:= ns/p.seqnc-ns]
	#	define adaptive chunks with at least 5 miss.pr>0
	YXm		<- merge(YXm, YXm[, list(ns.per.t=sum(ns)), by=c('t','t.stAgeC')], by=c('t','t.stAgeC'))
	#subset(YXm, t.stAgeC=='UC_(-1,25]' & t>2010.3)[,
	tmp		<- YXm[, 
			{
				z	<- sapply(seq_along(t), function(i)	sort(abs(t[i]-t[miss.pr>0]))		)				
				z	<- z[ ifelse(nrow(z)<20, nrow(z), 20), ]
				#print(z)
				z	<- sapply(seq_along(t), function(i) c(	'miss.pr.per.ch'= sum(miss.pr[ abs(t[i]-t)<=z[i] ]),
									't.per.ch'=length(which(abs(t[i]-t)<=z[i])),
									'tn.per.ch'=sum(tn[ abs(t[i]-t)<=z[i] ]), 
									'ns.per.ch'=sum(ns[ abs(t[i]-t)<=z[i] ])				))
				#print(z)				
				list(t=t, Patient=Patient, miss.pr.per.ch=z['miss.pr.per.ch',], t.per.ch=z['t.per.ch',], tn.per.ch=z['tn.per.ch',], ns.per.ch=z['ns.per.ch',])				
			}, by='t.stAgeC']
	YXm		<- merge(YXm, tmp, by=c('t','t.stAgeC','Patient'))
	#YXm[, miss.prsm:= miss.pr.per.ch/tn.per.ch]
	YXm[, miss.prsm:= miss.pr.per.ch/t.per.ch]
	YXm		<- merge(YXm, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	
	#	plot smooth
	if(!is.na(indir) & !is.na(infile))
	{
		ggplot(YXm, aes(x=t, colour=t.AgeC)) + geom_line(aes(y=pseqnc.pr)) + facet_grid(t.AgeC~stageC) + theme_bw() +
				geom_smooth(data=YXc, aes(y=p.seqnc), colour='black', method='loess', span=0.7) +
				theme_bw()
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpovershootBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)			
		tmp		<- subset(YXm, select=c(t, stageC, t.AgeC, t.stAgeC, miss.prsm, nr))
		setkey(tmp, t, t.stAgeC)
		tmp		<- unique(tmp)	
		ggplot( tmp, aes(x=t, y=miss.prsm, colour=t.AgeC, group=t.stAgeC)) + geom_step() + geom_point() +
				theme_bw() + theme(legend.position='bottom') +
				facet_wrap(~stageC, scales='free', ncol=YXm[, length(unique(stageC))])
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=15,h=7)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr, fill=t.AgeC)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient') +
				facet_grid(t.AgeC~stageC, scales='free')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothtotalBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothtotal',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=5,h=4)	
	}
	YXm
}
######################################################################################
miss.smoothing.150916<- function(YXc, mm, md, rltvp, rltva, indir=NA, infile=NA, suffix='')
{
	#	get scaffold: for all recipients at time t, all stages
	YXm		<- YXc[, list(nr=length(Patient)), by='t']
	YXm		<- merge(YXm, YXc[, list(tn=length(t.Patient)), by=c('t','Patient')], by='t')
	YXm		<- merge( as.data.table(expand.grid(t.stAgeC=YXc[, levels(t.stAgeC)], t=YXc[, sort(unique(t))])), YXm, by='t', all.x=TRUE, allow.cartesian=TRUE )
	#	add number observed intervals by t and stage
	YXm		<- merge(YXm, YXc[, list(ns=length(t.Patient)), by=c('t','Patient','t.stAgeC')], by=c('t','Patient','t.stAgeC'), all.x=TRUE)
	set(YXm, YXm[, which(is.na(ns))], c('ns'), 0)	
	setkey(YXm, t, Patient)	
	#	get smoothed w's and p's
	YXm[, pseqnc.pr:=predict(mm, data=md, newdata=subset(YXm, select=c(t, t.stAgeC)), type='response', what='mu')]
	tmp		<- rltvp(subset(YXm, select=c(Patient,t,t.stAgeC)), rltva)
	YXm		<- merge(YXm, tmp, by=c('Patient','t','t.stAgeC'))	
	set(YXm, YXm[, which(scoreY.pr<0)], 'scoreY.pr', 0)
	set(YXm, YXm[, which(pseqnc.pr<0)], 'pseqnc.pr', 0)	
	#	if we have an observed sampling prob, use that
	tmp		<- YXc[, list(p.seqnc=median(p.seqnc)), by=c('Patient','t.stAgeC','t')]
	YXm		<- merge(YXm, tmp, by=c('Patient','t.stAgeC','t'), all.x=1)
	tmp		<- YXm[, which(is.na(p.seqnc))]
	set(YXm, tmp, 'p.seqnc', YXm[tmp, pseqnc.pr])
	#	predict missing based on the n's -- allow for large tails in NB
	YXm[, miss.pr:= ns/p.seqnc-ns]
	#	define adaptive chunks with at least 5 miss.pr>0
	YXm		<- merge(YXm, YXm[, list(ns.per.t=sum(ns)), by=c('t','t.stAgeC')], by=c('t','t.stAgeC'))
	#subset(YXm, t.stAgeC=='UC_(-1,25]' & t>2010.3)[,
	tmp		<- YXm[, 
			{
				z	<- sapply(seq_along(t), function(i)	sort(abs(t[i]-t[miss.pr>0]))		)				
				z	<- z[ ifelse(nrow(z)<20, nrow(z), 20), ]
				#print(z)
				z	<- sapply(seq_along(t), function(i) c(	'miss.pr.per.ch'= sum(miss.pr[ abs(t[i]-t)<=z[i] ]),
									't.per.ch'=length(which(abs(t[i]-t)<=z[i])),
									'tn.per.ch'=sum(tn[ abs(t[i]-t)<=z[i] ]), 
									'ns.per.ch'=sum(ns[ abs(t[i]-t)<=z[i] ])				))
				#print(z)				
				list(t=t, Patient=Patient, miss.pr.per.ch=z['miss.pr.per.ch',], t.per.ch=z['t.per.ch',], tn.per.ch=z['tn.per.ch',], ns.per.ch=z['ns.per.ch',])				
			}, by='t.stAgeC']
	YXm		<- merge(YXm, tmp, by=c('t','t.stAgeC','Patient'))
	#YXm[, miss.prsm:= miss.pr.per.ch/tn.per.ch]
	YXm[, miss.prsm:= miss.pr.per.ch/t.per.ch]
	YXm		<- merge(YXm, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	
	#	plot smooth
	if(!is.na(indir) & !is.na(infile))
	{
		ggplot(YXm, aes(x=t, colour=t.AgeC)) + geom_line(aes(y=pseqnc.pr)) + facet_grid(t.AgeC~stageC) + theme_bw() +
				geom_smooth(data=YXc, aes(y=p.seqnc), colour='black', method='loess', span=0.7) +
				theme_bw()
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpovershootBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)			
		tmp		<- subset(YXm, select=c(t, stageC, t.AgeC, t.stAgeC, miss.prsm, nr))
		setkey(tmp, t, t.stAgeC)
		tmp		<- unique(tmp)	
		ggplot( tmp, aes(x=t, y=miss.prsm, colour=t.AgeC, group=t.stAgeC)) + geom_step() + geom_point() +
				theme_bw() + theme(legend.position='bottom') +
				facet_wrap(~stageC, scales='free', ncol=YXm[, length(unique(stageC))])
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=15,h=7)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr, fill=t.AgeC)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient') +
				facet_grid(t.AgeC~stageC, scales='free')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothtotalBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothtotal',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=5,h=4)	
	}
	YXm
}
######################################################################################
miss.smoothing.150918<- function(YXc, mm, md, rltvp, rltvd, rltvm, indir=NA, infile=NA, suffix='',tperiod=0.125, use.obs.pseqnc=TRUE)
{
	#	get scaffold: for all recipients at time t, all stages
	YXm		<- YXc[, list(nr=length(Patient)), by='t']
	YXm		<- merge(YXm, YXc[, list(tn=length(t.Patient)), by=c('t','Patient')], by='t')
	YXm		<- merge( as.data.table(expand.grid(t.stAgeC=YXc[, levels(t.stAgeC)], t=YXc[, sort(unique(t))])), YXm, by='t', all.x=TRUE, allow.cartesian=TRUE )
	#	add number observed intervals by t and stage
	YXm		<- merge(YXm, YXc[, list(ns=length(t.Patient)), by=c('t','Patient','t.stAgeC')], by=c('t','Patient','t.stAgeC'), all.x=TRUE)
	set(YXm, YXm[, which(is.na(ns))], c('ns'), 0)	
	setkey(YXm, t, Patient)	
	#	get smoothed w's and p's	
	YXm		<- merge(YXm, YXc[, list(tmin=min(t), tmax=max(t)), by='t.stAgeC'], by='t.stAgeC')
	YXm.pr	<- subset(YXm, t>tmin & t<tmax)	
	YXm.pr[, pseqnc.pr:=predict(mm, data=md, newdata=subset(YXm.pr, select=c(t, t.stAgeC)), type='response', what='mu')]
	setnames(YXm.pr, 't', 'tmid')
	YXm.pr	<- rltvp(YXm.pr, rltvm, rltvd)
	setnames(YXm.pr, 'tmid', 't')
	set(YXm.pr, NULL, 't.stAgeC2', NULL)	
	YXm		<- merge(YXm, subset(YXm.pr, select=c(t.stAgeC, t, Patient, pseqnc.pr, scoreY.pr)), by=c('t.stAgeC','t','Patient'), all.x=TRUE)
	#	extrapolate prediction
	tmp		<- unique(subset(YXm, t==tmin+tperiod, select=c(t.stAgeC, pseqnc.pr, scoreY.pr)))
	setnames(tmp, c('pseqnc.pr','scoreY.pr'),c('pseqnc.pr.min','scoreY.pr.min'))
	YXm		<- merge(YXm, tmp, by=c('t.stAgeC'), all.x=TRUE)
	tmp		<- unique(subset(YXm, t==tmax-tperiod, select=c(t.stAgeC, pseqnc.pr, scoreY.pr)))
	setnames(tmp, c('pseqnc.pr','scoreY.pr'),c('pseqnc.pr.max','scoreY.pr.max'))
	YXm		<- merge(YXm, tmp, by=c('t.stAgeC'), all.x=TRUE)
	stopifnot( YXm[, length(which(is.na(pseqnc.pr) & is.na(pseqnc.pr.min) & is.na(pseqnc.pr.max)))==0] )
	tmp		<- YXm[, which(t<=tmin & is.na(pseqnc.pr))]
	set(YXm, tmp, 'pseqnc.pr', YXm[tmp, pseqnc.pr.min])
	set(YXm, tmp, 'scoreY.pr', YXm[tmp, scoreY.pr.min])
	tmp		<- YXm[, which(t>=tmax & is.na(pseqnc.pr))]
	set(YXm, tmp, 'pseqnc.pr', YXm[tmp, pseqnc.pr.max])
	set(YXm, tmp, 'scoreY.pr', YXm[tmp, scoreY.pr.max])
	YXm		<- subset(YXm, select=c(t.stAgeC, t, Patient, nr, tn, ns, pseqnc.pr, scoreY.pr))
	#	check
	cat('\nFound scoreY.pr<0, n=', length(YXm[, which(scoreY.pr<0)]))
	set(YXm, YXm[, which(scoreY.pr<0)], 'scoreY.pr', 0)
	cat('\nFound pseqnc.pr<0, n=', length(YXm[, which(pseqnc.pr<0)]))
	set(YXm, YXm[, which(pseqnc.pr<0)], 'pseqnc.pr', 0)	
	cat('\nFound pseqnc.pr>1, n=', length(YXm[, which(pseqnc.pr>1)]))
	set(YXm, YXm[, which(pseqnc.pr>1)], 'pseqnc.pr', 1)	
	#
	tmp		<- dcast.data.table(YXm, Patient+t~t.stAgeC, value.var='pseqnc.pr')
	setkey(tmp, t)
	tmp		<- unique(tmp)
	#	if we have an observed sampling prob, use that
	if(use.obs.pseqnc)
	{
		tmp		<- YXc[, list(p.seqnc=median(p.seqnc)), by=c('Patient','t.stAgeC','t')]
		YXm		<- merge(YXm, tmp, by=c('Patient','t.stAgeC','t'), all.x=1)
		tmp		<- YXm[, which(is.na(p.seqnc))]
		set(YXm, tmp, 'p.seqnc', YXm[tmp, pseqnc.pr])		
	}
	if(!use.obs.pseqnc)
		YXm[, p.seqnc:= pseqnc.pr]	
	#	predict missing based on the n's -- allow for large tails in NB
	YXm[, miss.pr:= ns/p.seqnc-ns]
	#	define adaptive chunks with at least 5 miss.pr>0
	YXm		<- merge(YXm, YXm[, list(ns.per.t=sum(ns)), by=c('t','t.stAgeC')], by=c('t','t.stAgeC'))
	#subset(YXm, t.stAgeC=='UC_(-1,25]' & t>2010.3)[,
	tmp		<- YXm[, 
			{
				z	<- sapply(seq_along(t), function(i)	sort(abs(t[i]-t[miss.pr>0]))		)				
				z	<- z[ ifelse(nrow(z)<20, nrow(z), 20), ]
				#print(z)
				z	<- sapply(seq_along(t), function(i) c(	'miss.pr.per.ch'= sum(miss.pr[ abs(t[i]-t)<=z[i] ]),
									't.per.ch'=length(which(abs(t[i]-t)<=z[i])),
									'tn.per.ch'=sum(tn[ abs(t[i]-t)<=z[i] ]), 
									'ns.per.ch'=sum(ns[ abs(t[i]-t)<=z[i] ])				))
				#print(z)				
				list(t=t, Patient=Patient, miss.pr.per.ch=z['miss.pr.per.ch',], t.per.ch=z['t.per.ch',], tn.per.ch=z['tn.per.ch',], ns.per.ch=z['ns.per.ch',])				
			}, by='t.stAgeC']
	YXm		<- merge(YXm, tmp, by=c('t','t.stAgeC','Patient'))
	#YXm[, miss.prsm:= miss.pr.per.ch/tn.per.ch]
	YXm[, miss.prsm:= miss.pr.per.ch/t.per.ch]
	YXm		<- merge(YXm, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')	
	#	plot smooth
	if(!is.na(indir) & !is.na(infile))
	{
		ggplot(YXm, aes(x=t, colour=t.AgeC)) + geom_line(aes(y=pseqnc.pr)) + facet_grid(t.AgeC~stageC) + theme_bw() +
				geom_smooth(data=YXc, aes(y=p.seqnc), colour='black', method='loess', span=0.7) +
				theme_bw()
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpovershootBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)			
		tmp		<- subset(YXm, select=c(t, stageC, t.AgeC, t.stAgeC, miss.prsm, nr))
		setkey(tmp, t, t.stAgeC)
		tmp		<- unique(tmp)	
		ggplot( tmp, aes(x=t, y=miss.prsm, colour=t.AgeC, group=t.stAgeC)) + geom_step() + geom_point() +
				theme_bw() + theme(legend.position='bottom') +
				facet_wrap(~stageC, scales='free', ncol=YXm[, length(unique(stageC))])
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=15,h=7)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr, fill=t.AgeC)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient') +
				facet_grid(t.AgeC~stageC, scales='free')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothtotalBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothtotal',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=5,h=4)	
	}
	YXm
}
######################################################################################
miss.smoothing.150923<- function(YXc, mp, mma, mmb, md, rltvp, rltvd, rltvm, indir=NA, infile=NA, suffix='',tperiod=0.125, use.obs.pseqnc=TRUE)
{
	#	get scaffold: for all recipients at time t, all stages
	YXm		<- YXc[, list(nr=length(Patient)), by='t']
	YXm		<- merge(YXm, YXc[, list(tn=length(t.Patient)), by=c('t','Patient')], by='t')
	YXm		<- merge( as.data.table(expand.grid(t.stAgeC=YXc[, levels(t.stAgeC)], t=YXc[, sort(unique(t))])), YXm, by='t', all.x=TRUE, allow.cartesian=TRUE )
	YXm		<- merge(YXm, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')	
	#	add number observed intervals by t and stage
	YXm		<- merge(YXm, YXc[, list(ns=length(t.Patient)), by=c('t','Patient','t.stAgeC')], by=c('t','Patient','t.stAgeC'), all.x=TRUE)
	set(YXm, YXm[, which(is.na(ns))], c('ns'), 0)	
	setkey(YXm, t, Patient)	
	#	get smoothed w's and p's	
	YXm		<- merge(YXm, YXc[, list(tmin=min(t), tmax=max(t)), by='t.stAgeC'], by='t.stAgeC')
	YXm.pr	<- subset(YXm, t>tmin & t<tmax)	
	YXm.pr	<- mp(YXm.pr, mma, mmb, md)
	stopifnot(YXm.pr[, !any(is.na(pseqnc.pr))])
	setnames(YXm.pr, 't', 'tmid')
	YXm.pr	<- rltvp(YXm.pr, rltvm, rltvd)
	stopifnot(YXm.pr[, !any(is.na(scoreY.pr))])
	setnames(YXm.pr, 'tmid', 't')
	set(YXm.pr, NULL, 't.stAgeC2', NULL)	
	YXm		<- merge(YXm, subset(YXm.pr, select=c(t.stAgeC, t, Patient, pseqnc.pr, scoreY.pr)), by=c('t.stAgeC','t','Patient'), all.x=TRUE)
	#	extrapolate prediction
	tmp		<- unique(subset(YXm, t==tmin+tperiod, select=c(t.stAgeC, pseqnc.pr, scoreY.pr)))
	setnames(tmp, c('pseqnc.pr','scoreY.pr'),c('pseqnc.pr.min','scoreY.pr.min'))
	YXm		<- merge(YXm, tmp, by=c('t.stAgeC'), all.x=TRUE)
	tmp		<- unique(subset(YXm, t==tmax-tperiod, select=c(t.stAgeC, pseqnc.pr, scoreY.pr)))
	setnames(tmp, c('pseqnc.pr','scoreY.pr'),c('pseqnc.pr.max','scoreY.pr.max'))
	YXm		<- merge(YXm, tmp, by=c('t.stAgeC'), all.x=TRUE)
	#stopifnot( YXm[, length(which(is.na(pseqnc.pr) & is.na(pseqnc.pr.min) & is.na(pseqnc.pr.max)))==0] )
	tmp		<- YXm[, which(t<=tmin & is.na(pseqnc.pr))]
	set(YXm, tmp, 'pseqnc.pr', YXm[tmp, pseqnc.pr.min])
	set(YXm, tmp, 'scoreY.pr', YXm[tmp, scoreY.pr.min])
	tmp		<- YXm[, which(t>=tmax & is.na(pseqnc.pr))]
	set(YXm, tmp, 'pseqnc.pr', YXm[tmp, pseqnc.pr.max])
	set(YXm, tmp, 'scoreY.pr', YXm[tmp, scoreY.pr.max])
	YXm		<- subset(YXm, select=c(t.stAgeC, t, Patient, nr, tn, ns, pseqnc.pr, scoreY.pr))
	#	fixup 	
	if(YXm[, any(t.stAgeC=='L_(55,100]')])
	{
		tmp	<- dcast.data.table(YXm, t+Patient~t.stAgeC, value.var='pseqnc.pr')
		set(tmp, NULL, 'L_(55,100]', tmp[, 'L_(45,55]', with=0])
		tmp	<- melt(tmp, id.vars=c('t','Patient'), variable.name='t.stAgeC', value.name='pseqnc.pr')
		YXm	<- merge(subset(YXm, select=c(t.stAgeC, t, Patient, nr, tn, ns, scoreY.pr)), tmp, by=c('t','Patient','t.stAgeC'))
		tmp	<- dcast.data.table(YXm, t+Patient~t.stAgeC, value.var='scoreY.pr')
		set(tmp, NULL, 'L_(55,100]', tmp[, 'L_(45,55]', with=0])
		tmp	<- melt(tmp, id.vars=c('t','Patient'), variable.name='t.stAgeC', value.name='scoreY.pr')
		YXm	<- merge(subset(YXm, select=c(t.stAgeC, t, Patient, nr, tn, ns, pseqnc.pr)), tmp, by=c('t','Patient','t.stAgeC'))
	}
	#	check
	cat('\nFound scoreY.pr<0, n=', length(YXm[, which(scoreY.pr<0)]))
	set(YXm, YXm[, which(scoreY.pr<0)], 'scoreY.pr', 0)
	cat('\nFound pseqnc.pr<0, n=', length(YXm[, which(pseqnc.pr<0)]))
	set(YXm, YXm[, which(pseqnc.pr<0)], 'pseqnc.pr', 0)	
	cat('\nFound pseqnc.pr>1, n=', length(YXm[, which(pseqnc.pr>1)]))
	set(YXm, YXm[, which(pseqnc.pr>1)], 'pseqnc.pr', 1)	
	#
	tmp		<- dcast.data.table(YXm, Patient+t~t.stAgeC, value.var='pseqnc.pr')
	setkey(tmp, t)
	tmp		<- unique(tmp)
	#	if we have an observed sampling prob, use that
	if(use.obs.pseqnc)
	{
		tmp		<- YXc[, list(p.seqnc=median(p.seqnc)), by=c('Patient','t.stAgeC','t')]
		YXm		<- merge(YXm, tmp, by=c('Patient','t.stAgeC','t'), all.x=1)
		tmp		<- YXm[, which(is.na(p.seqnc))]
		set(YXm, tmp, 'p.seqnc', YXm[tmp, pseqnc.pr])		
	}
	if(!use.obs.pseqnc)
		YXm[, p.seqnc:= pseqnc.pr]	
	#	predict missing based on the n's -- allow for large tails in NB
	YXm[, miss.pr:= ns/p.seqnc-ns]
	#	define adaptive chunks with at least 5 miss.pr>0
	YXm		<- merge(YXm, YXm[, list(ns.per.t=sum(ns)), by=c('t','t.stAgeC')], by=c('t','t.stAgeC'))
	#subset(YXm, t.stAgeC=='UC_(-1,25]' & t>2010.3)[,
	tmp		<- YXm[, 
			{
				z	<- sapply(seq_along(t), function(i)	sort(abs(t[i]-t[miss.pr>0]))		)				
				z	<- z[ ifelse(nrow(z)<20, nrow(z), 20), ]
				#print(z)
				z	<- sapply(seq_along(t), function(i) c(	'miss.pr.per.ch'= sum(miss.pr[ abs(t[i]-t)<=z[i] ]),
									't.per.ch'=length(which(abs(t[i]-t)<=z[i])),
									'tn.per.ch'=sum(tn[ abs(t[i]-t)<=z[i] ]), 
									'ns.per.ch'=sum(ns[ abs(t[i]-t)<=z[i] ])				))
				#print(z)				
				list(t=t, Patient=Patient, miss.pr.per.ch=z['miss.pr.per.ch',], t.per.ch=z['t.per.ch',], tn.per.ch=z['tn.per.ch',], ns.per.ch=z['ns.per.ch',])				
			}, by='t.stAgeC']
	YXm		<- merge(YXm, tmp, by=c('t','t.stAgeC','Patient'))
	#YXm[, miss.prsm:= miss.pr.per.ch/tn.per.ch]
	YXm[, miss.prsm:= miss.pr.per.ch/t.per.ch]
	YXm		<- merge(YXm, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')	
	#	plot smooth
	if(!is.na(indir) & !is.na(infile))
	{
		ggplot(YXm, aes(x=t, colour=t.AgeC)) + geom_line(aes(y=pseqnc.pr)) + facet_grid(t.AgeC~stageC) + theme_bw() +
				geom_smooth(data=subset(YXc, !t.stAgeC%in%c('L_(55,100]','UAE_(55,100]','UAC_(55,100]','TS_(-1,25]')), aes(y=p.seqnc), colour='black', method='loess', span=0.7) +
				theme_bw()
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpovershootBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)			
		tmp		<- subset(YXm, select=c(t, stageC, t.AgeC, t.stAgeC, miss.prsm, nr))
		setkey(tmp, t, t.stAgeC)
		tmp		<- unique(tmp)	
		ggplot( tmp, aes(x=t, y=miss.prsm, colour=t.AgeC, group=t.stAgeC)) + geom_step() + geom_point() +
				theme_bw() + theme(legend.position='bottom') +
				facet_wrap(~stageC, scales='free', ncol=YXm[, length(unique(stageC))])
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=15,h=7)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr, fill=t.AgeC)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient') +
				facet_grid(t.AgeC~stageC, scales='free')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothtotalBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothtotal',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=5,h=4)	
	}
	YXm
}
######################################################################################
miss.smoothing.150922<- function(YXc, mp, mma, mmb, md, rltvp, rltvd, rltvm, indir=NA, infile=NA, suffix='',tperiod=0.125, use.obs.pseqnc=TRUE)
{
	#	get scaffold: for all recipients at time t, all stages
	YXm		<- YXc[, list(nr=length(Patient)), by='t']
	YXm		<- merge(YXm, YXc[, list(tn=length(t.Patient)), by=c('t','Patient')], by='t')
	YXm		<- merge( as.data.table(expand.grid(t.stAgeC=YXc[, levels(t.stAgeC)], t=YXc[, sort(unique(t))])), YXm, by='t', all.x=TRUE, allow.cartesian=TRUE )
	YXm		<- merge(YXm, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')	
	#	add number observed intervals by t and stage
	YXm		<- merge(YXm, YXc[, list(ns=length(t.Patient)), by=c('t','Patient','t.stAgeC')], by=c('t','Patient','t.stAgeC'), all.x=TRUE)
	set(YXm, YXm[, which(is.na(ns))], c('ns'), 0)	
	setkey(YXm, t, Patient)	
	#	get smoothed w's and p's	
	YXm		<- merge(YXm, YXc[, list(tmin=min(t), tmax=max(t)), by='t.stAgeC'], by='t.stAgeC')
	YXm.pr	<- subset(YXm, t>tmin & t<tmax)	
	YXm.pr	<- mp(YXm.pr, mma, mmb, md)
	stopifnot(YXm.pr[, !any(is.na(pseqnc.pr))])
	setnames(YXm.pr, 't', 'tmid')
	YXm.pr	<- rltvp(YXm.pr, rltvm, rltvd)
	stopifnot(YXm.pr[, !any(is.na(scoreY.pr))])
	setnames(YXm.pr, 'tmid', 't')
	set(YXm.pr, NULL, 't.stAgeC2', NULL)	
	YXm		<- merge(YXm, subset(YXm.pr, select=c(t.stAgeC, t, Patient, pseqnc.pr, scoreY.pr)), by=c('t.stAgeC','t','Patient'), all.x=TRUE)
	#	extrapolate prediction
	tmp		<- unique(subset(YXm, t==tmin+tperiod, select=c(t.stAgeC, pseqnc.pr, scoreY.pr)))
	setnames(tmp, c('pseqnc.pr','scoreY.pr'),c('pseqnc.pr.min','scoreY.pr.min'))
	YXm		<- merge(YXm, tmp, by=c('t.stAgeC'), all.x=TRUE)
	tmp		<- unique(subset(YXm, t==tmax-tperiod, select=c(t.stAgeC, pseqnc.pr, scoreY.pr)))
	setnames(tmp, c('pseqnc.pr','scoreY.pr'),c('pseqnc.pr.max','scoreY.pr.max'))
	YXm		<- merge(YXm, tmp, by=c('t.stAgeC'), all.x=TRUE)
	stopifnot( YXm[, length(which(is.na(pseqnc.pr) & is.na(pseqnc.pr.min) & is.na(pseqnc.pr.max)))==0] )
	tmp		<- YXm[, which(t<=tmin & is.na(pseqnc.pr))]
	set(YXm, tmp, 'pseqnc.pr', YXm[tmp, pseqnc.pr.min])
	set(YXm, tmp, 'scoreY.pr', YXm[tmp, scoreY.pr.min])
	tmp		<- YXm[, which(t>=tmax & is.na(pseqnc.pr))]
	set(YXm, tmp, 'pseqnc.pr', YXm[tmp, pseqnc.pr.max])
	set(YXm, tmp, 'scoreY.pr', YXm[tmp, scoreY.pr.max])
	YXm		<- subset(YXm, select=c(t.stAgeC, t, Patient, nr, tn, ns, pseqnc.pr, scoreY.pr))	
	#	check
	cat('\nFound scoreY.pr<0, n=', length(YXm[, which(scoreY.pr<0)]))
	set(YXm, YXm[, which(scoreY.pr<0)], 'scoreY.pr', 0)
	cat('\nFound pseqnc.pr<0, n=', length(YXm[, which(pseqnc.pr<0)]))
	set(YXm, YXm[, which(pseqnc.pr<0)], 'pseqnc.pr', 0)	
	cat('\nFound pseqnc.pr>1, n=', length(YXm[, which(pseqnc.pr>1)]))
	set(YXm, YXm[, which(pseqnc.pr>1)], 'pseqnc.pr', 1)	
	#
	tmp		<- dcast.data.table(YXm, Patient+t~t.stAgeC, value.var='pseqnc.pr')
	setkey(tmp, t)
	tmp		<- unique(tmp)
	#	if we have an observed sampling prob, use that
	if(use.obs.pseqnc)
	{
		tmp		<- YXc[, list(p.seqnc=median(p.seqnc)), by=c('Patient','t.stAgeC','t')]
		YXm		<- merge(YXm, tmp, by=c('Patient','t.stAgeC','t'), all.x=1)
		tmp		<- YXm[, which(is.na(p.seqnc))]
		set(YXm, tmp, 'p.seqnc', YXm[tmp, pseqnc.pr])		
	}
	if(!use.obs.pseqnc)
		YXm[, p.seqnc:= pseqnc.pr]	
	#	predict missing based on the n's -- allow for large tails in NB
	YXm[, miss.pr:= ns/p.seqnc-ns]
	#	define adaptive chunks with at least 5 miss.pr>0
	YXm		<- merge(YXm, YXm[, list(ns.per.t=sum(ns)), by=c('t','t.stAgeC')], by=c('t','t.stAgeC'))
	#subset(YXm, t.stAgeC=='UC_(-1,25]' & t>2010.3)[,
	tmp		<- YXm[, 
			{
				z	<- sapply(seq_along(t), function(i)	sort(abs(t[i]-t[miss.pr>0]))		)				
				z	<- z[ ifelse(nrow(z)<20, nrow(z), 20), ]
				#print(z)
				z	<- sapply(seq_along(t), function(i) c(	'miss.pr.per.ch'= sum(miss.pr[ abs(t[i]-t)<=z[i] ]),
									't.per.ch'=length(which(abs(t[i]-t)<=z[i])),
									'tn.per.ch'=sum(tn[ abs(t[i]-t)<=z[i] ]), 
									'ns.per.ch'=sum(ns[ abs(t[i]-t)<=z[i] ])				))
				#print(z)				
				list(t=t, Patient=Patient, miss.pr.per.ch=z['miss.pr.per.ch',], t.per.ch=z['t.per.ch',], tn.per.ch=z['tn.per.ch',], ns.per.ch=z['ns.per.ch',])				
			}, by='t.stAgeC']
	YXm		<- merge(YXm, tmp, by=c('t','t.stAgeC','Patient'))
	#YXm[, miss.prsm:= miss.pr.per.ch/tn.per.ch]
	YXm[, miss.prsm:= miss.pr.per.ch/t.per.ch]
	YXm		<- merge(YXm, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')	
	#	plot smooth
	if(!is.na(indir) & !is.na(infile))
	{
		ggplot(YXm, aes(x=t, colour=t.AgeC)) + geom_line(aes(y=pseqnc.pr)) + facet_grid(t.AgeC~stageC) + theme_bw() +
				geom_smooth(data=YXc, aes(y=p.seqnc), colour='black', method='loess', span=0.7) +
				theme_bw()
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpovershootBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)			
		tmp		<- subset(YXm, select=c(t, stageC, t.AgeC, t.stAgeC, miss.prsm, nr))
		setkey(tmp, t, t.stAgeC)
		tmp		<- unique(tmp)	
		ggplot( tmp, aes(x=t, y=miss.prsm, colour=t.AgeC, group=t.stAgeC)) + geom_step() + geom_point() +
				theme_bw() + theme(legend.position='bottom') +
				facet_wrap(~stageC, scales='free', ncol=YXm[, length(unique(stageC))])
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=15,h=7)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr, fill=t.AgeC)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient') +
				facet_grid(t.AgeC~stageC, scales='free')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothtotalBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothtotal',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=5,h=4)	
	}
	YXm
}
######################################################################################
miss.smoothing.150917<- function(YXc, mm, md, rltvp, rltva, indir=NA, infile=NA, suffix='',tperiod=0.125, use.obs.pseqnc=TRUE)
{
	#	get scaffold: for all recipients at time t, all stages
	YXm		<- YXc[, list(nr=length(Patient)), by='t']
	YXm		<- merge(YXm, YXc[, list(tn=length(t.Patient)), by=c('t','Patient')], by='t')
	YXm		<- merge( as.data.table(expand.grid(t.stAgeC=YXc[, levels(t.stAgeC)], t=YXc[, sort(unique(t))])), YXm, by='t', all.x=TRUE, allow.cartesian=TRUE )
	#	add number observed intervals by t and stage
	YXm		<- merge(YXm, YXc[, list(ns=length(t.Patient)), by=c('t','Patient','t.stAgeC')], by=c('t','Patient','t.stAgeC'), all.x=TRUE)
	set(YXm, YXm[, which(is.na(ns))], c('ns'), 0)	
	setkey(YXm, t, Patient)	
	#	get smoothed w's and p's	
	YXm		<- merge(YXm, YXc[, list(tmin=min(t), tmax=max(t)), by='t.stAgeC'], by='t.stAgeC')
	YXm.pr	<- subset(YXm, t>tmin & t<tmax)	
	YXm.pr[, pseqnc.pr:=predict(mm, data=md, newdata=subset(YXm.pr, select=c(t, t.stAgeC)), type='response', what='mu')]
	YXm		<- merge(YXm, subset(YXm.pr, select=c(t.stAgeC, t, Patient, pseqnc.pr)), by=c('t.stAgeC','t','Patient'), all.x=TRUE)
		
	tmp		<- unique(subset(YXm, t==tmin+tperiod, select=c(t.stAgeC, pseqnc.pr)))
	setnames(tmp, 'pseqnc.pr','pseqnc.pr.min')
	YXm		<- merge(YXm, tmp, by=c('t.stAgeC'), all.x=TRUE)
	tmp		<- unique(subset(YXm, t==tmax-tperiod, select=c(t.stAgeC, pseqnc.pr)))
	setnames(tmp, 'pseqnc.pr','pseqnc.pr.max')
	YXm		<- merge(YXm, tmp, by=c('t.stAgeC'), all.x=TRUE)
	stopifnot( YXm[, length(which(is.na(pseqnc.pr) & is.na(pseqnc.pr.min) & is.na(pseqnc.pr.max)))==0] )
	tmp		<- YXm[, which(t<=tmin & is.na(pseqnc.pr))]
	set(YXm, tmp, 'pseqnc.pr', YXm[tmp, pseqnc.pr.min])
	tmp		<- YXm[, which(t>=tmax & is.na(pseqnc.pr))]
	set(YXm, tmp, 'pseqnc.pr', YXm[tmp, pseqnc.pr.max])
	YXm		<- subset(YXm, select=c(t.stAgeC, t, Patient, nr, tn, ns, pseqnc.pr))
	
	tmp		<- rltvp(subset(YXm, select=c(Patient,t,t.stAgeC)), rltva)
	YXm		<- merge(YXm, tmp, by=c('Patient','t','t.stAgeC'))	
	set(YXm, YXm[, which(scoreY.pr<0)], 'scoreY.pr', 0)
	set(YXm, YXm[, which(pseqnc.pr<0)], 'pseqnc.pr', 0)	
	
	tmp		<- dcast.data.table(YXm, Patient+t~t.stAgeC, value.var='pseqnc.pr')
	setkey(tmp, t)
	tmp		<- unique(tmp)
	#	if we have an observed sampling prob, use that
	if(use.obs.pseqnc)
	{
		tmp		<- YXc[, list(p.seqnc=median(p.seqnc)), by=c('Patient','t.stAgeC','t')]
		YXm		<- merge(YXm, tmp, by=c('Patient','t.stAgeC','t'), all.x=1)
		tmp		<- YXm[, which(is.na(p.seqnc))]
		set(YXm, tmp, 'p.seqnc', YXm[tmp, pseqnc.pr])		
	}
	if(!use.obs.pseqnc)
		YXm[, p.seqnc:= pseqnc.pr]	
	#	predict missing based on the n's -- allow for large tails in NB
	YXm[, miss.pr:= ns/p.seqnc-ns]
	#	define adaptive chunks with at least 5 miss.pr>0
	YXm		<- merge(YXm, YXm[, list(ns.per.t=sum(ns)), by=c('t','t.stAgeC')], by=c('t','t.stAgeC'))
	#subset(YXm, t.stAgeC=='UC_(-1,25]' & t>2010.3)[,
	tmp		<- YXm[, 
			{
				z	<- sapply(seq_along(t), function(i)	sort(abs(t[i]-t[miss.pr>0]))		)				
				z	<- z[ ifelse(nrow(z)<20, nrow(z), 20), ]
				#print(z)
				z	<- sapply(seq_along(t), function(i) c(	'miss.pr.per.ch'= sum(miss.pr[ abs(t[i]-t)<=z[i] ]),
									't.per.ch'=length(which(abs(t[i]-t)<=z[i])),
									'tn.per.ch'=sum(tn[ abs(t[i]-t)<=z[i] ]), 
									'ns.per.ch'=sum(ns[ abs(t[i]-t)<=z[i] ])				))
				#print(z)				
				list(t=t, Patient=Patient, miss.pr.per.ch=z['miss.pr.per.ch',], t.per.ch=z['t.per.ch',], tn.per.ch=z['tn.per.ch',], ns.per.ch=z['ns.per.ch',])				
			}, by='t.stAgeC']
	YXm		<- merge(YXm, tmp, by=c('t','t.stAgeC','Patient'))
	#YXm[, miss.prsm:= miss.pr.per.ch/tn.per.ch]
	YXm[, miss.prsm:= miss.pr.per.ch/t.per.ch]
	YXm		<- merge(YXm, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	
	#	plot smooth
	if(!is.na(indir) & !is.na(infile))
	{
		ggplot(YXm, aes(x=t, colour=t.AgeC)) + geom_line(aes(y=pseqnc.pr)) + facet_grid(t.AgeC~stageC) + theme_bw() +
				geom_smooth(data=YXc, aes(y=p.seqnc), colour='black', method='loess', span=0.7) +
				theme_bw()
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpovershootBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)			
		tmp		<- subset(YXm, select=c(t, stageC, t.AgeC, t.stAgeC, miss.prsm, nr))
		setkey(tmp, t, t.stAgeC)
		tmp		<- unique(tmp)	
		ggplot( tmp, aes(x=t, y=miss.prsm, colour=t.AgeC, group=t.stAgeC)) + geom_step() + geom_point() +
				theme_bw() + theme(legend.position='bottom') +
				facet_wrap(~stageC, scales='free', ncol=YXm[, length(unique(stageC))])
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=15,h=7)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr, fill=t.AgeC)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient') +
				facet_grid(t.AgeC~stageC, scales='free')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothtotalBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpmsmoothtotal',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=5,h=4)	
	}
	YXm
}
######################################################################################
miss.smoothing.150909<- function(YXc, mm, md, rltvm, rltvd, indir=NA, infile=NA, suffix='')
{
	#	get scaffold: for all recipients at time t, all stages
	YXm		<- YXc[, list(nr=length(Patient)), by='t']
	YXm		<- merge(YXm, YXc[, list(tn=length(t.Patient)), by=c('t','Patient')], by='t')
	YXm		<- merge( as.data.table(expand.grid(t.stAgeC=YXc[, levels(t.stAgeC)], t=YXc[, sort(unique(t))])), YXm, by='t', all.x=TRUE, allow.cartesian=TRUE )
	#	add number observed intervals by t and stage
	YXm		<- merge(YXm, YXc[, list(ns=length(t.Patient)), by=c('t','Patient','t.stAgeC')], by=c('t','Patient','t.stAgeC'), all.x=TRUE)
	set(YXm, YXm[, which(is.na(ns))], c('ns'), 0)	
	setkey(YXm, t, Patient)	
	#	get smoothed w's and p's
	YXm[, pseqnc.pr:=predict(mm, data=md, newdata=subset(YXm, select=c(t, t.stAgeC)), type='response', what='mu')]
	YXm[, scoreY.pr:=predict(rltvm, data=rltvd, newdata=subset(YXm, select=c(t, t.stAgeC)), type='response', what='mu')]
	set(YXm, YXm[, which(scoreY.pr<0)], 'scoreY.pr', 0)
	set(YXm, YXm[, which(pseqnc.pr<0)], 'pseqnc.pr', 0)	
	#	predict missing based on the n's -- allow for large tails in NB
	YXm[, miss.pr:= ns/pseqnc.pr-ns]
	#	define adaptive chunks with at least 5 miss.pr>0
	YXm		<- merge(YXm, YXm[, list(ns.per.t=sum(ns)), by=c('t','t.stAgeC')], by=c('t','t.stAgeC'))
	#subset(YXm, t.stAgeC=='UC_(-1,25]' & t>2010.3)[,
	tmp		<- YXm[, 
			{
				z	<- sapply(seq_along(t), function(i)	sort(abs(t[i]-t[miss.pr>0]))		)				
				z	<- z[ ifelse(nrow(z)<20, nrow(z), 20), ]
				#print(z)
				z	<- sapply(seq_along(t), function(i) c(	'miss.pr.per.ch'= sum(miss.pr[ abs(t[i]-t)<=z[i] ]), 
									'tn.per.ch'=sum(tn[ abs(t[i]-t)<=z[i] ]), 
									'ns.per.ch'=sum(ns[ abs(t[i]-t)<=z[i] ])				))
				#print(z)				
				list(t=t, Patient=Patient, miss.pr.per.ch=z['miss.pr.per.ch',], tn.per.ch=z['tn.per.ch',], ns.per.ch=z['ns.per.ch',])				
			}, by='t.stAgeC']
	YXm		<- merge(YXm, tmp, by=c('t','t.stAgeC','Patient'))
	YXm[, miss.prsm:= miss.pr.per.ch/tn.per.ch]
	YXm		<- merge(YXm, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	#	plot smooth
	if(!is.na(indir) & !is.na(infile))
	{
		tmp		<- subset(YXm, select=c(t, stageC, t.AgeC, t.stAgeC, miss.prsm, nr))
		setkey(tmp, t, t.stAgeC)
		tmp		<- unique(tmp)	
		ggplot( tmp, aes(x=t, y=miss.prsm, colour=t.AgeC, group=t.stAgeC)) + geom_step() + geom_point() +
				theme_bw() + theme(legend.position='bottom') +
				facet_wrap(~stageC, scales='free', ncol=YXm[, length(unique(stageC))])
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpsmoothBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=15,h=7)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr, fill=t.AgeC)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient') +
				facet_grid(t.AgeC~stageC, scales='free')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpsmoothtotalBytstAgeC',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=12,h=12)	
		ggplot( YXm, aes(x=t, y=miss.prsm/nr)) + geom_bar(stat="identity") + 
				theme_bw() + theme(legend.position='bottom') +
				labs(y='total miss.prsm per t, standardized per recipient')
		file	<- paste(indir,'/',gsub('\\.R',paste('_missexpsmoothtotal',suffix,'.pdf',sep=''),infile),sep='')	
		ggsave(file=file, w=5,h=4)	
	}
	YXm
}
######################################################################################
miss.explore.150907<- function(YXc)
{
	#	get scaffold: for all recipients at time t, all stages
	YXm		<- YXc[, list(tn=length(t.Patient)), by=c('t','Patient')]
	YXm		<- merge( as.data.table(expand.grid(t.stAgeC=YXc[, levels(t.stAgeC)], t=YXc[, sort(unique(t))])), YXm, by='t', all.x=TRUE, allow.cartesian=TRUE )
	#	add number observed intervals by t and stage
	YXm		<- merge(YXm, YXc[, list(ns=length(t.Patient)), by=c('t','Patient','t.stAgeC')], by=c('t','Patient','t.stAgeC'), all.x=TRUE)
	set(YXm, YXm[, which(is.na(ns))], c('ns'), 0)	
	setkey(YXm, t, Patient)	
	#	get smoothed w's and p's
	YXm[, pseqnc.pr:=predict(mm, data=md, newdata=subset(YXm, select=c(t, t.stAgeC)), type='response', what='mu')]
	YXm[, scoreY.pr:=predict(rltvm, data=rltvd, newdata=subset(YXm, select=c(t, t.stAgeC)), type='response', what='mu')]
	set(YXm, YXm[, which(scoreY.pr<0)], 'scoreY.pr', 0)
	set(YXm, YXm[, which(pseqnc.pr<0)], 'pseqnc.pr', 0)	
	#	predict missing based on the n's -- allow for large tails in NB
	YXm[, miss.pr:= ns/pseqnc.pr-ns]
	#	define adaptive chunks with at least 5 miss.pr>0
	YXm		<- merge(YXm, YXm[, list(ns.per.t=sum(ns)), by=c('t','t.stAgeC')], by=c('t','t.stAgeC'))
	#	define adaptive chunks with at least one ns.per.tsm
	tsm.period	<- 0.25
	YXm[, tsm:= YXm[, floor(t) + floor( (t%%1)*100 %/% (tsm.period*100) ) * tsm.period]]	
	YXm		<- merge(YXm, YXm[, list(ns.per.t=sum(ns)), by=c('tsm','t.stAgeC')], by=c('tsm','t.stAgeC'))
	tmp		<- unique(subset(YXm, select=c(tsm, t.stAgeC, ns.per.t)))
	tmp		<- tmp[, list(tsm=tsm, adaptivechunk_id= cumsum(ns.per.t>0) ), by='t.stAgeC']
	YXm		<- merge(YXm, tmp, by=c('tsm','t.stAgeC'))
	#	smooth across time
	YXm		<- merge(YXm, YXm[, list(tn.per.ch= sum(tn), ns.per.ch=sum(ns), miss.pr.per.ch=sum(miss.pr)), by=c('t.stAgeC','adaptivechunk_id')], by=c('t.stAgeC','adaptivechunk_id'))
	YXm[, miss.prsm:= miss.pr.per.ch/tn.per.ch]
	YXm		<- merge(YXm, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	#	plot smooth
	tmp		<- subset(YXm, select=c(t, stageC, t.AgeC, t.stAgeC, miss.prsm))
	setkey(tmp, t, t.stAgeC)
	tmp		<- unique(tmp)	
	ggplot( tmp, aes(x=t, y=miss.prsm, colour=t.AgeC, group=t.stAgeC)) + geom_step() + geom_point() +
			theme_bw() + theme(legend.position='bottom') +
			facet_wrap(~stageC, scales='free', ncol=5)
	file	<- paste(indir,'/',gsub('\\.R','_missexpsmoothBytstAgeC.pdf',infile),sep='')	
	ggsave(file=file, w=15,h=7)	
	#
	#	ALTERNATIVE
	#
	#	define adaptive chunks with at least one ns.per.t
	YXm		<- merge(YXm, YXm[, list(ns.per.t=sum(ns)), by=c('t','t.stAgeC')], by=c('t','t.stAgeC'))
	tmp		<- unique(subset(YXm, select=c(t, t.stAgeC, ns.per.t)))
	tmp		<- tmp[, list(t=t, adaptivechunk_id= cumsum(ns.per.t>0) ), by='t.stAgeC']
	YXm		<- merge(YXm, tmp, by=c('t','t.stAgeC'))
	#	smooth across time
	YXm		<- merge(YXm, YXm[, list(tn.per.ch= sum(tn), ns.per.ch=sum(ns), miss.pr.per.ch=sum(miss.pr)), by=c('t.stAgeC','adaptivechunk_id')], by=c('t.stAgeC','adaptivechunk_id'))
	YXm[, miss.prsm:= miss.pr.per.ch/tn.per.ch]
	YXm		<- merge(YXm, unique(subset(YXc, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
	ggplot(subset(YXm, t.stAgeC=='UC_(-1,25]'), aes(x=t, y=miss.prsm)) + geom_bar(stat="identity")
	file	<- paste(indir,'/',gsub('\\.R','_missexpsmoothtotalBytstAgeC.pdf',infile),sep='')	
	ggsave(file=file, w=10,h=10)	
}
######################################################################################
miss.explore.150922Gen<- function(YXc, indir, infile, suffix='')
{
	tmp		<- melt(YXc, measure.vars=c('p.seq','p.nc','p.seqnc','me'))
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC)) + geom_boxplot() +
			facet_grid(variable~stageC, scales='free') + theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSforcensoringsampling_empirical',suffix,'\\.pdf',sep=''),infile),sep=''), w=10, h=14)	
	ggplot(YXc, aes(x=t, y=me, colour=t.AgeC)) + geom_point(alpha=0.5) + geom_smooth(colour='black', fill='grey50', alpha=1) +
			facet_grid(t.AgeC~stageC, scales='free') + theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSforcensoringsampling_smooth',suffix,'\\.pdf',sep=''),infile),sep=''), w=10, h=14)
	
	
	YXc		<- subset(YXc, AnyPos_T1<2010.8)
	YXm		<- subset(YXc, select=c(t, stageC, t.AgeC, t.stAgeC, p.seqnc))
	YXm[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXm[, unique(gsub('1956-65','1920-65',gsub('1920-55','1920-65',levels(t.stAgeC),fixed=1),fixed=1))]	
	set(YXm, NULL, 't.stAgeC2', YXm[, factor(gsub('1956-65','1920-65',gsub('1920-55','1920-65',t.stAgeC2,fixed=1),fixed=1), levels=tmp, labels=tmp)])
	#YXm[, t.stAgeC3:= t.stAgeC2]
	#YXm[, t.stAgeC4:= t.stAgeC2]
	#set(YXm, YXm[, which(stageC=='UAC')], 't.stAgeC3', NA_character_)
	#set(YXm, YXm[, which(stageC!='UAC')], 't.stAgeC4', NA_character_)
	
	mm1		<- gamlss(p.seqnc~t.stAgeC2+t.stAgeC2:t, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm2		<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=2):t.stAgeC2, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm3		<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=4):t.stAgeC2, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm4		<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=3):t.stAgeC2, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm5a	<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=4):t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')), family=GA())
	mm5b	<- gamlss(p.seqnc~t.stAgeC2+t:t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, stageC%in%c('UAC','TO','TS','L')), family=GA())
	
	YXm2	<- as.data.table(expand.grid(t=YXm[, sort(unique(t))], t.stAgeC2=YXm[, levels(t.stAgeC2)], stringsAsFactors=0))
	YXm2	<- merge(YXm2, unique(subset(YXm, select=c(stageC, t.AgeC, t.stAgeC, t.stAgeC2))), by='t.stAgeC2', allow.cartesian=1)	
	me1		<- predict(mm1,data=YXm,newdata=YXm2,type='response',what='mu')
	me2		<- predict(mm2,data=YXm,newdata=YXm2,type='response',what='mu')
	me3		<- predict(mm3,data=YXm,newdata=YXm2,type='response',what='mu')
	me4		<- predict(mm4,data=YXm,newdata=YXm2,type='response',what='mu')
	me5a	<- predict(mm5a,data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')),newdata=YXm2,type='response',what='mu')
	me5b	<- predict(mm5b,data=subset(YXm, stageC%in%c('UAC','TO','TS','L')),newdata=YXm2,type='response',what='mu')
	YXm2[, p.seqnc1:=me1]
	YXm2[, p.seqnc2:=me2]
	YXm2[, p.seqnc3:=me3]
	YXm2[, p.seqnc4:=me4]
	YXm2[, p.seqnc5:=me5a]
	tmp		<- YXm2[, which(stageC%in%c('UAC','TO','TS','L'))]
	set(YXm2, tmp, 'p.seqnc5', me5b[tmp])
	
	YXm2	<- merge(YXm2, YXm, by=c('t','stageC','t.AgeC','t.stAgeC','t.stAgeC2'),all=1,allow.cartesian=1)
	#YXm		<- melt(YXm, measure.vars=c('me1','me2','me3','me4'))
	ggplot(melt(YXm2, measure.vars=c('p.seqnc1','p.seqnc2','p.seqnc4','p.seqnc5')), aes(x=t, colour=t.AgeC, group=t.stAgeC)) + 
			geom_point(aes(y=p.seqnc), alpha=0.4) + geom_smooth(aes(y=p.seqnc), colour='black', fill='grey50', alpha=1) +
			geom_line(aes(y=value), colour='red') +
			facet_grid(t.AgeC~variable+stageC, scales='free') +
			theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSpmodels',suffix,'\\.pdf',sep=''),infile),sep=''), w=45, h=15)			
}
######################################################################################
miss.explore.150923Gen<- function(YXc, indir, infile, suffix='')
{
	tmp		<- melt(YXc, measure.vars=c('p.seq','p.nc','p.seqnc','me'))
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC)) + geom_boxplot() +
			facet_grid(variable~stageC, scales='free') + theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSforcensoringsampling_empirical',suffix,'\\.pdf',sep=''),infile),sep=''), w=10, h=14)	
	ggplot(YXc, aes(x=t, y=me, colour=t.AgeC)) + geom_point(alpha=0.5) + geom_smooth(colour='black', fill='grey50', alpha=1) +
			facet_grid(t.AgeC~stageC, scales='free') + theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSforcensoringsampling_smooth',suffix,'\\.pdf',sep=''),infile),sep=''), w=10, h=14)
	
	
	YXc		<- subset(YXc, AnyPos_T1<2010.8)
	YXm		<- subset(YXc, select=c(t, stageC, t.AgeC, t.stAgeC, p.seqnc))
	YXm[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXm[, unique(gsub('TO_.*','TO_1920-95',gsub('TS_.*','TS_1920-95',gsub('L_.*','L_1920-95',gsub('1956-65','1920-65',gsub('1920-55','1920-65',levels(t.stAgeC),fixed=1),fixed=1)))))]
	set(YXm, NULL, 't.stAgeC2', YXm[, factor(gsub('TO_.*','TO_1920-95',gsub('TS_.*','TS_1920-95',gsub('L_.*','L_1920-95',gsub('1956-65','1920-65',gsub('1920-55','1920-65',t.stAgeC2,fixed=1),fixed=1)))), levels=tmp, labels=tmp)])
		
	#YXm[, t.stAgeC3:= t.stAgeC2]
	#YXm[, t.stAgeC4:= t.stAgeC2]
	#set(YXm, YXm[, which(stageC=='UAC')], 't.stAgeC3', NA_character_)
	#set(YXm, YXm[, which(stageC!='UAC')], 't.stAgeC4', NA_character_)
	
	mm1		<- gamlss(p.seqnc~t.stAgeC2+t.stAgeC2:t, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm2		<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=2):t.stAgeC2, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm3		<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=4):t.stAgeC2, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm4		<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=3):t.stAgeC2, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm5a	<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=3):t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')), family=GA())
	mm5b	<- gamlss(p.seqnc~t.stAgeC2+t:t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, stageC%in%c('UAC','TO','TS','L')), family=GA())
	
	YXm2	<- as.data.table(expand.grid(t=YXm[, sort(unique(t))], t.stAgeC2=YXm[, levels(t.stAgeC2)], stringsAsFactors=0))
	YXm2	<- merge(YXm2, unique(subset(YXm, select=c(stageC, t.AgeC, t.stAgeC, t.stAgeC2))), by='t.stAgeC2', allow.cartesian=1)	
	me1		<- predict(mm1,data=YXm,newdata=YXm2,type='response',what='mu')
	me2		<- predict(mm2,data=YXm,newdata=YXm2,type='response',what='mu')
	me3		<- predict(mm3,data=YXm,newdata=YXm2,type='response',what='mu')
	me4		<- predict(mm4,data=YXm,newdata=YXm2,type='response',what='mu')
	me5a	<- predict(mm5a,data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')),newdata=YXm2,type='response',what='mu')
	me5b	<- predict(mm5b,data=subset(YXm, stageC%in%c('UAC','TO','TS','L')),newdata=YXm2,type='response',what='mu')
	YXm2[, p.seqnc1:=me1]
	YXm2[, p.seqnc2:=me2]
	YXm2[, p.seqnc3:=me3]
	YXm2[, p.seqnc4:=me4]
	YXm2[, p.seqnc5:=me5a]
	tmp		<- YXm2[, which(stageC%in%c('UAC','TO','TS','L'))]
	set(YXm2, tmp, 'p.seqnc5', me5b[tmp])
	
	YXm2	<- merge(YXm2, YXm, by=c('t','stageC','t.AgeC','t.stAgeC','t.stAgeC2'),all=1,allow.cartesian=1)
	#YXm		<- melt(YXm, measure.vars=c('me1','me2','me3','me4'))
	ggplot(melt(YXm2, measure.vars=c('p.seqnc3','p.seqnc4','p.seqnc5')), aes(x=t, colour=t.AgeC, group=t.stAgeC)) + 
			geom_point(aes(y=p.seqnc), alpha=0.4) + geom_smooth(aes(y=p.seqnc), colour='black', fill='grey50', alpha=1) +
			geom_line(aes(y=value), colour='red') +
			facet_grid(t.AgeC~variable+stageC, scales='free') +
			theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSpmodels',suffix,'\\.pdf',sep=''),infile),sep=''), w=35, h=15)			
}
######################################################################################
miss.explore.150923<- function(YXc, indir, infile, suffix='')
{
	tmp		<- melt(YXc, measure.vars=c('p.seq','p.nc','p.seqnc','me'))
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC)) + geom_boxplot() +
			facet_grid(variable~stageC, scales='free') + theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSforcensoringsampling_empirical',suffix,'\\.pdf',sep=''),infile),sep=''), w=10, h=14)	
	ggplot(YXc, aes(x=t, y=me, colour=t.AgeC)) + geom_point(alpha=0.5) + geom_smooth(colour='black', fill='grey50', alpha=1) +
			facet_grid(t.AgeC~stageC, scales='free') + theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSforcensoringsampling_smooth',suffix,'\\.pdf',sep=''),infile),sep=''), w=10, h=14)
	
	
	YXc		<- subset(YXc, AnyPos_T1<2010.8)
	YXm		<- subset(YXc, select=c(t, stageC, t.AgeC, t.stAgeC, p.seqnc))
	YXm[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXm[, unique(gsub('TO_.*','TO_(-1,00]',gsub('TS_.*','TS_(-1,00]',gsub('L_.*','L_(-1,00]',gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',levels(t.stAgeC),fixed=1),fixed=1)))))]	
	set(YXm, NULL, 't.stAgeC2', YXm[, factor(gsub('TO_.*','TO_(-1,00]',gsub('TS_.*','TS_(-1,00]',gsub('L_.*','L_(-1,00]',gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',t.stAgeC2,fixed=1),fixed=1)))), levels=tmp, labels=tmp)])
	#YXm[, t.stAgeC3:= t.stAgeC2]
	#YXm[, t.stAgeC4:= t.stAgeC2]
	#set(YXm, YXm[, which(stageC=='UAC')], 't.stAgeC3', NA_character_)
	#set(YXm, YXm[, which(stageC!='UAC')], 't.stAgeC4', NA_character_)	
	mm1		<- gamlss(p.seqnc~t.stAgeC2+t.stAgeC2:t, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm2		<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=2):t.stAgeC2, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm3		<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=4):t.stAgeC2, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm4		<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=3):t.stAgeC2, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm5a	<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=4):t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')), family=GA())
	mm5b	<- gamlss(p.seqnc~t.stAgeC2+t:t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, stageC%in%c('UAC','TO','TS','L')), family=GA())
	
	YXm2	<- as.data.table(expand.grid(t=YXm[, sort(unique(t))], t.stAgeC2=YXm[, levels(t.stAgeC2)], stringsAsFactors=0))
	YXm2	<- merge(YXm2, unique(subset(YXm, select=c(stageC, t.AgeC, t.stAgeC, t.stAgeC2))), by='t.stAgeC2', allow.cartesian=1)	
	me1		<- predict(mm1,data=YXm,newdata=YXm2,type='response',what='mu')
	me2		<- predict(mm2,data=YXm,newdata=YXm2,type='response',what='mu')
	me3		<- predict(mm3,data=YXm,newdata=YXm2,type='response',what='mu')
	me4		<- predict(mm4,data=YXm,newdata=YXm2,type='response',what='mu')
	me5a	<- predict(mm5a,data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')),newdata=YXm2,type='response',what='mu')
	me5b	<- predict(mm5b,data=subset(YXm, stageC%in%c('UAC','TO','TS','L')),newdata=YXm2,type='response',what='mu')
	YXm2[, p.seqnc1:=me1]
	YXm2[, p.seqnc2:=me2]
	YXm2[, p.seqnc3:=me3]
	YXm2[, p.seqnc4:=me4]
	YXm2[, p.seqnc5:=me5a]
	tmp		<- YXm2[, which(stageC%in%c('UAC','TO','TS','L'))]
	set(YXm2, tmp, 'p.seqnc5', me5b[tmp])
	
	YXm2	<- merge(YXm2, YXm, by=c('t','stageC','t.AgeC','t.stAgeC','t.stAgeC2'),all=1,allow.cartesian=1)
	#YXm		<- melt(YXm, measure.vars=c('me1','me2','me3','me4'))
	ggplot(melt(YXm2, measure.vars=c('p.seqnc1','p.seqnc2','p.seqnc4','p.seqnc5')), aes(x=t, colour=t.AgeC, group=t.stAgeC)) + 
			geom_point(aes(y=p.seqnc), alpha=0.4) + geom_smooth(aes(y=p.seqnc), colour='black', fill='grey50', alpha=1) +
			geom_line(aes(y=value), colour='red') +
			coord_cartesian(ylim=c(0,1)) +
			facet_grid(t.AgeC~variable+stageC, scales='free') +
			theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSpmodels',suffix,'\\.pdf',sep=''),infile),sep=''), w=45, h=15)			
}
######################################################################################
miss.explore.150922<- function(YXc, indir, infile, suffix='')
{
	tmp		<- melt(YXc, measure.vars=c('p.seq','p.nc','p.seqnc','me'))
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC)) + geom_boxplot() +
			facet_grid(variable~stageC, scales='free') + theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSforcensoringsampling_empirical',suffix,'\\.pdf',sep=''),infile),sep=''), w=10, h=14)	
	ggplot(YXc, aes(x=t, y=me, colour=t.AgeC)) + geom_point(alpha=0.5) + geom_smooth(colour='black', fill='grey50', alpha=1) +
			facet_grid(t.AgeC~stageC, scales='free') + theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSforcensoringsampling_smooth',suffix,'\\.pdf',sep=''),infile),sep=''), w=10, h=14)
	
	
	YXc		<- subset(YXc, AnyPos_T1<2010.8)
	YXm		<- subset(YXc, select=c(t, stageC, t.AgeC, t.stAgeC, p.seqnc))
	YXm[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXm[, unique(gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',levels(t.stAgeC),fixed=1),fixed=1))]
	set(YXm, NULL, 't.stAgeC2', YXm[, factor(gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',t.stAgeC2,fixed=1),fixed=1), levels=tmp, labels=tmp)])
	#YXm[, t.stAgeC3:= t.stAgeC2]
	#YXm[, t.stAgeC4:= t.stAgeC2]
	#set(YXm, YXm[, which(stageC=='UAC')], 't.stAgeC3', NA_character_)
	#set(YXm, YXm[, which(stageC!='UAC')], 't.stAgeC4', NA_character_)
	
	mm1		<- gamlss(p.seqnc~t.stAgeC2+t.stAgeC2:t, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm2		<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=2):t.stAgeC2, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm3		<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=4):t.stAgeC2, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm4		<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=3):t.stAgeC2, sigma.formula=~t.stAgeC2, data=YXm, family=GA())
	mm5a	<- gamlss(p.seqnc~t.stAgeC2+ns(t, df=4):t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')), family=GA())
	mm5b	<- gamlss(p.seqnc~t.stAgeC2+t:t.stAgeC2+t:t.stAgeC2, sigma.formula=~t.stAgeC2, data=subset(YXm, stageC%in%c('UAC','TO','TS','L')), family=GA())
	
	YXm2	<- as.data.table(expand.grid(t=YXm[, sort(unique(t))], t.stAgeC2=YXm[, levels(t.stAgeC2)], stringsAsFactors=0))
	YXm2	<- merge(YXm2, unique(subset(YXm, select=c(stageC, t.AgeC, t.stAgeC, t.stAgeC2))), by='t.stAgeC2', allow.cartesian=1)	
	me1		<- predict(mm1,data=YXm,newdata=YXm2,type='response',what='mu')
	me2		<- predict(mm2,data=YXm,newdata=YXm2,type='response',what='mu')
	me3		<- predict(mm3,data=YXm,newdata=YXm2,type='response',what='mu')
	me4		<- predict(mm4,data=YXm,newdata=YXm2,type='response',what='mu')
	me5a	<- predict(mm5a,data=subset(YXm, !stageC%in%c('UAC','TO','TS','L')),newdata=YXm2,type='response',what='mu')
	me5b	<- predict(mm5b,data=subset(YXm, stageC%in%c('UAC','TO','TS','L')),newdata=YXm2,type='response',what='mu')
	YXm2[, p.seqnc1:=me1]
	YXm2[, p.seqnc2:=me2]
	YXm2[, p.seqnc3:=me3]
	YXm2[, p.seqnc4:=me4]
	YXm2[, p.seqnc5:=me5a]
	tmp		<- YXm2[, which(stageC%in%c('UAC','TO','TS','L'))]
	set(YXm2, tmp, 'p.seqnc5', me5b[tmp])
	
	YXm2	<- merge(YXm2, YXm, by=c('t','stageC','t.AgeC','t.stAgeC','t.stAgeC2'),all=1,allow.cartesian=1)
	#YXm		<- melt(YXm, measure.vars=c('me1','me2','me3','me4'))
	ggplot(melt(YXm2, measure.vars=c('p.seqnc1','p.seqnc2','p.seqnc4','p.seqnc5')), aes(x=t, colour=t.AgeC, group=t.stAgeC)) + 
			geom_point(aes(y=p.seqnc), alpha=0.4) + geom_smooth(aes(y=p.seqnc), colour='black', fill='grey50', alpha=1) +
			geom_line(aes(y=value), colour='red') +
			coord_cartesian(ylim=c(0,1)) +
			facet_grid(t.AgeC~variable+stageC, scales='free') +
			theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSpmodels',suffix,'\\.pdf',sep=''),infile),sep=''), w=45, h=15)			
}
######################################################################################
miss.explore.150904<- function(YXc, indir, infile, suffix='')
{
	tmp		<- melt(YXc, measure.vars=c('p.seq','p.nc','p.seqnc','me'))
	ggplot(tmp, aes(x=t.period, y=value, fill=stageC)) + geom_boxplot() +
			facet_grid(variable~stageC, scales='free') + theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSforcensoringsampling_empirical',suffix,'\\.pdf',sep=''),infile),sep=''), w=10, h=14)	
	ggplot(YXc, aes(x=t, y=me, colour=t.AgeC)) + geom_point(alpha=0.5) + geom_smooth(colour='black', fill='grey50', alpha=1) +
			facet_grid(t.AgeC~stageC, scales='free') + theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSforcensoringsampling_smooth',suffix,'\\.pdf',sep=''),infile),sep=''), w=10, h=14)
	
	
	YXc		<- subset(YXc, AnyPos_T1<2010.8)
	YXm		<- subset(YXc, select=c(t, stageC, t.AgeC, t.stAgeC, p.seqnc))
	set(YXm, NULL, 't.stAgeC', YXm[, factor(as.character(t.stAgeC))])	
	mm1		<- gamlss(p.seqnc~t.stAgeC+t.stAgeC:t, sigma.formula=~t.stAgeC, data=YXm, family=GA())
	mm2		<- gamlss(p.seqnc~t.stAgeC+ns(t, df=2):t.stAgeC, sigma.formula=~t.stAgeC, data=YXm, family=GA())
	mm3		<- gamlss(p.seqnc~t.stAgeC+ns(t, df=4):t.stAgeC, sigma.formula=~t.stAgeC, data=YXm, family=GA())
	mm4		<- gamlss(p.seqnc~t.stAgeC+ns(t, df=3):t.stAgeC, sigma.formula=~t.stAgeC, data=YXm, family=GA())	
	me1		<- predict(mm1,type='response',what='mu')
	me2		<- predict(mm2,type='response',what='mu')
	me3		<- predict(mm3,type='response',what='mu')
	me4		<- predict(mm4,type='response',what='mu')
	YXm[, p.seqnc1:=me1]
	YXm[, p.seqnc2:=me2]
	YXm[, p.seqnc3:=me3]
	YXm[, p.seqnc4:=me4]
	#YXm		<- melt(YXm, measure.vars=c('me1','me2','me3','me4'))
	ggplot(melt(YXm, measure.vars=c('p.seqnc1','p.seqnc2','p.seqnc3','p.seqnc4')), aes(x=t, colour=t.AgeC, group=t.stAgeC)) + 
			geom_point(aes(y=p.seqnc), alpha=0.4) + geom_smooth(aes(y=p.seqnc), colour='black', fill='grey50', alpha=1) +
			geom_line(aes(y=value), colour='red') +
			facet_grid(t.AgeC~variable+stageC, scales='free') +
			theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSpmodels',suffix,'\\.pdf',sep=''),infile),sep=''), w=45, h=15)
	
	
	YXc		<- subset(YXc, AnyPos_T1<2010.8)		
	YXm		<- subset(YXc, select=c(t, stageC, t.AgeC, t.stAgeC, me))
	set(YXm, NULL, 't.stAgeC', YXm[, factor(as.character(t.stAgeC))])
	mm1		<- gamlss(me~t.stAgeC+t.stAgeC:t, data=YXm, family=GA())
	mm2		<- gamlss(me~t.stAgeC+t.stAgeC:t, sigma.formula=~t.stAgeC, data=YXm, family=GA())
	mm3		<- gamlss(me~t.stAgeC+ns(t, df=4):t.stAgeC, sigma.formula=~t.stAgeC, data=YXm, family=GA())
	mm4		<- gamlss(me~t.stAgeC+ns(t, df=3):t.stAgeC, sigma.formula=~t.stAgeC, data=YXm, family=GA())	
	me1		<- predict(mm1,type='response',what='mu')
	me2		<- predict(mm2,type='response',what='mu')
	me3		<- predict(mm3,type='response',what='mu')
	me4		<- predict(mm4,type='response',what='mu')
	YXm[, me1:=me1]
	YXm[, me2:=me2]
	YXm[, me3:=me3]
	YXm[, me4:=me4]
	#YXm		<- melt(YXm, measure.vars=c('me1','me2','me3','me4'))
	ggplot(melt(YXm, measure.vars=c('me1','me2','me3','me4')), aes(x=t, colour=t.AgeC, group=t.stAgeC)) + 
			geom_point(aes(y=me), alpha=0.4) + geom_smooth(aes(y=me), colour='black', fill='grey50', alpha=1) +
			geom_line(aes(y=value), colour='red') +
			facet_grid(t.AgeC~variable+stageC, scales='free') +
			theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R',paste('_EXPMISSmodels',suffix,'\\.pdf',sep=''),infile),sep=''), w=35, h=15)		
}
######################################################################################
miss.explore.150831<- function(YXc, indir, infile)
{
	YXc		<- subset(YXc, AnyPos_T1<2010.8)
	tmp		<- YXc[, list(tm= sum(1/(p.nc)-1)/length(unique(Patient)), rn=length(unique(Patient)) ), by='t']
	ggplot(tmp, aes(x=t, y=tm)) + geom_line() + labs(y='exp missing intervals per recipient')
	ggsave(file=paste(indir,'/',gsub('\\.R','_EXPMISSforcensoring\\.pdf',infile),sep=''), w=5, h=4)
	tmp		<- YXc[, list(tm= sum(1/(p.seq)-1)/length(unique(Patient)), rn=length(unique(Patient)) ), by='t']
	ggplot(tmp, aes(x=t, y=tm)) + geom_line() + labs(y='exp missing intervals per recipient')
	ggsave(file=paste(indir,'/',gsub('\\.R','_EXPMISSforsampling\\.pdf',infile),sep=''), w=5, h=4)	
	tmp		<- YXc[, list(ttlme= sum(1/(p.seq*p.nc)-1)/length(unique(Patient)), nR=length(unique(Patient)) ), by='t']
	setkey(tmp, t)
	ggplot(tmp, aes(x=t, y=ttlme)) + geom_line() + labs(y='exp missing intervals per recipient')
	ggsave(file=paste(indir,'/',gsub('\\.R','_EXPMISSforcensoringsampling\\.pdf',infile),sep=''), w=5, h=4)	
	YXc		<- merge(YXc, tmp,by='t')
	tmp		<- YXc[, list(me= sum(1/(p.seq*p.nc)-1)/nR[1] ), by=c('t','t.stAgeC')]
	YXc		<- merge(YXc, tmp, by=c('t','t.stAgeC'))
	YXm		<- subset(YXc, t<2010.75, select=c(t, stageC, t.AgeC, t.stAgeC, me, ttlme))
	set(YXm, NULL, 't.stAgeC', YXm[, factor(as.character(t.stAgeC))])
	mm1		<- gamlss(me~t.stAgeC+t.stAgeC:t, data=YXm, family=GA())
	mm2		<- gamlss(me~t.stAgeC+t.stAgeC:t, sigma.formula=~t.stAgeC, data=YXm, family=GA())
	mm3		<- gamlss(me~t.stAgeC+ns(t, df=4):t.stAgeC, sigma.formula=~t.stAgeC, data=YXm, family=GA())
	mm4		<- gamlss(me~t.stAgeC+ns(t, knots=c(2005,2009,2010,2010.5)):t.stAgeC, sigma.formula=~t.stAgeC, data=YXm, family=GA())
	me1		<- predict(mm1,type='response',what='mu')
	me2		<- predict(mm2,type='response',what='mu')
	me3		<- predict(mm3,type='response',what='mu')
	me4		<- predict(mm4,type='response',what='mu')
	YXm[, me1:=me1]
	YXm[, me2:=me2]
	YXm[, me3:=me3]
	YXm[, me4:=me4]
	#YXm		<- melt(YXm, measure.vars=c('me1','me2','me3','me4'))
	ggplot(melt(YXm, measure.vars=c('me1','me2','me3','me4')), aes(x=t, colour=t.AgeC, group=t.stAgeC)) + geom_point(aes(y=me)) +
			geom_line(aes(y=value)) +
			facet_wrap(variable~stageC, scales='free', ncol=5) +
			theme_bw()
	ggsave(file=paste(indir,'/',gsub('\\.R','_EXPMISSmodels\\.pdf',infile),sep=''), w=15, h=10)	
	
	tmp2	<- as.data.table(expand.grid(t.stAgeC=YXm[, levels(t.stAgeC)], t=YXm[, sort(unique(t))]))
	tmp2[, mem:=predict(mm3, data=YXm, newdata=tmp2, type='response', what='mu')]
	ggplot(tmp2[, list(ttlmem=sum(mem)), by='t'], aes(x=t)) + 
			geom_line(aes(y=ttlmem), colour='red') +
			geom_line(data=YXm, aes(y=ttlme)) + labs(y='black line: empirical\nred line: predicted')
	ggsave(file=paste(indir,'/',gsub('\\.R','_EXPMISSmodelvsempirical\\.pdf',infile),sep=''), w=10, h=10)
}
######################################################################################
sampling.dev<- function()
{
	indir		<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tpairs_age'
	#
	#	for m5A
	#
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C2V100bInfT7STRAT_m5A.R'
	load(paste(indir,'/',infile,sep=''))
	stage.labels<- c('UA','UC','D','T','L')
	pt.df		<- data.table(Patient=X.msm[, unique(t.Patient)])
	ptd.df		<- sampling.get.dataset(pt.df, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, df.treatment.allmsm, t.period, t.endctime, method.lRNA.supp=2, method.lRNA.nsupp=4, contact.grace=1.5)
	ptr.df		<- sampling.model.150722(ptd.df)
	setnames(ptr.df, 'Patient', 't.Patient')
	riskcl		<- 't.stAgeC'	
	YXs			<- merge( YX, subset(ptr.df, select=c('t.Patient','p.seq')), all.x=T, by='t.Patient' )
	set(YXs, NULL, 'stageC', YXs[, factor(as.character(stageC), levels=stage.labels)])
	ggplot(YXs, aes(y=100*p.seq, x=stageC, fill=stageC)) +
			geom_boxplot() +
			#geom_point(size=3) +
			scale_y_continuous(breaks=seq(0,100,10), minor_breaks=NULL, limit=c(0,100), expand=c(0,0)) +
			scale_shape_discrete(guide=FALSE) +
			#scale_colour_manual(values=dfp[, unique(factor.color)], guide=FALSE) +
			theme_bw() +				
			theme(axis.text=element_text(size=10), axis.title=element_text(size=10), legend.position='bottom', panel.grid.major.x=element_line(colour="grey70", size=0.4), panel.grid.minor.x=element_line(colour="grey70", size=0.4), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), panel.margin = unit(2, "lines")) +
			labs(x='', y='sampling probability\namong prob transmitters\n(%)', colour='stage of potential transmitter\nat transmission interval') +
			coord_flip() +
			facet_grid(~t.AgeC) 
	plot.file	<- paste(indir,'/',gsub('\\.R','_SAMPLINGMODEL_ms12.pdf',infile),sep='')	   				
	ggsave(file=plot.file, w=10, h=5)
		
	set(YXs, NULL, 't.period.long', YXs[, factor(t.period, levels=c('1','2','3','4'), labels=c("96/07-06/06", "06/07-07/12", "08/01-09/06", "09/07-10/12"))])
	ggplot(YXs, aes(y=100*p.seq, x=t.period.long, fill=stageC)) +
			geom_boxplot(size=0.5, weight=0.5, outlier.size=0.5) +
			#geom_point(size=3) +
			scale_y_continuous(breaks=seq(0,100,10), minor_breaks=NULL, limit=c(0,100), expand=c(0,0)) +
			#scale_fill_brewer(palette='Set1') +
			scale_shape_discrete(guide=FALSE) +
			#scale_colour_manual(values=dfp[, unique(factor.color)], guide=FALSE) +
			theme_bw() +				
			theme(axis.text=element_text(size=10), axis.title=element_text(size=10), legend.position='bottom', panel.grid.major.x=element_line(colour="grey70", size=0.4), panel.grid.minor.x=element_line(colour="grey70", size=0.4), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), panel.margin = unit(2, "lines")) +
			labs(x='', y='sampling probability\namong prob transmitters\n(%)', fill='stage of potential transmitter\nat transmission interval') +
			coord_flip() +
			facet_grid(~t.AgeC) 
	plot.file	<- paste(indir,'/',gsub('\\.R','_SAMPLINGMODEL_ms12_by_tperiod.pdf',infile),sep='')	   				
	ggsave(file=plot.file, w=10, h=8)
	#
	#	for m5B
	#
	infile		<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_3pa1H1.48C2V100bInfT7STRAT_m5B.R'
	load(paste(indir,'/',infile,sep=''))
	stage.labels<- c("UAC", "UAE", "UC", "DAC", "D", "L", "T")	   
	pt.df		<- data.table(Patient=X.msm[, unique(t.Patient)])
	ptd.df		<- sampling.get.dataset(pt.df, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, df.treatment.allmsm, t.period, t.endctime, method.lRNA.supp=2, method.lRNA.nsupp=4, contact.grace=1.5)
	ptr.df		<- sampling.model.150722(ptd.df)
	setnames(ptr.df, 'Patient', 't.Patient')
	riskcl		<- 't.stAgeC'	
	YXs			<- merge( YX, subset(ptr.df, select=c('t.Patient','p.seq')), all.x=T, by='t.Patient' )
	set(YXs, NULL, 'stageC', YXs[, factor(as.character(stageC), levels=stage.labels)])
	set(YXs, NULL, 't.period.long', YXs[, factor(t.period, levels=c('1','2','3','4'), labels=c("96/07-06/06", "06/07-07/12", "08/01-09/06", "09/07-10/12"))])
	ggplot(YXs, aes(y=100*p.seq, x=t.period.long, fill=stageC)) +
			geom_boxplot(size=0.5, weight=0.5, outlier.size=0.5) +
			#geom_point(size=3) +
			scale_y_continuous(breaks=seq(0,100,10), minor_breaks=NULL, limit=c(0,100), expand=c(0,0)) +
			#scale_fill_brewer(palette='Set1') +
			scale_shape_discrete(guide=FALSE) +
			#scale_colour_manual(values=dfp[, unique(factor.color)], guide=FALSE) +
			theme_bw() +				
			theme(axis.text=element_text(size=10), axis.title=element_text(size=10), legend.position='bottom', panel.grid.major.x=element_line(colour="grey70", size=0.4), panel.grid.minor.x=element_line(colour="grey70", size=0.4), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), panel.margin = unit(2, "lines")) +
			labs(x='', y='sampling probability\namong prob transmitters\n(%)', fill='stage of potential transmitter\nat transmission interval') +
			coord_flip() +
			facet_grid(~t.AgeC) 
	plot.file	<- paste(indir,'/',gsub('\\.R','_SAMPLINGMODEL_ms12_by_tperiod.pdf',infile),sep='')	   				
	ggsave(file=plot.file, w=10, h=9)
	 
}
######################################################################################
sampling.get.dataset<- function(pt.df, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, df.treatment.allmsm, t.period, t.endctime, method.lRNA.supp=2, method.lRNA.nsupp=4, contact.grace=1.5, save.file=NA)
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
	ptr.df	<- copy(ptd.df)
	if(any(colnames(ptr.df)=='PosSeqT'))
		ptr.df[, PosSeqT:=NULL]
	if(any(colnames(ptr.df)=='AnyT_T1'))
		ptr.df[, AnyT_T1:=NULL]					
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
		factor.ref.v	<- paste('TS_(30,45]',tp,sep='')
		risktp.col		<- 't.stAgeC.prd'
		risk.col		<- 't.stAgeC'			
	}
	if(grepl('m5C',method))
	{
		factor.ref.v	<- paste('TS_(30,45]',tp,sep='')
		risktp.col		<- 't.stAgeC.prd'
		risk.col		<- 't.stAgeC'			
	}
	if(grepl('m5D',method))
	{
		factor.ref.v	<- paste('TS_(38,43]',tp,sep='')
		risktp.col		<- 't.stAgeC.prd'
		risk.col		<- 't.stAgeC'			
	}
	if(grepl('m5E',method))
	{
		factor.ref.v	<- paste('TS_(30,45]',tp,sep='')
		risktp.col		<- 't.stAgeC.prd'
		risk.col		<- 't.stAgeC'			
	}
	if(grepl('m5F',method))
	{
		factor.ref.v	<- paste('TS_1920-95',tp,sep='')
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
sampling.get.all.tables<- function(YX=NULL, X.seq=NULL, X.msm=NULL, X.clu=NULL, clumsm.info=NULL, tperiod.info=NULL, resume=TRUE, save.file=NA, method=NA, risk.col=NA, risktp.col=NA, factor.ref.v=NA, bs.n=100, bs.cdelta.min=2, bs.cdelta.max=3)
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
		if(is.null(YX) || is.null(X.seq) || is.null(X.msm) || is.null(X.clu))
		{
			cat(paste('\nreturn NULL table', method))
			return(ans)	
		}			
		stopifnot(!is.na(risk.col), !is.na(risktp.col), !is.na(factor.ref.v))
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
		tmp				<- subset(clumsm.info, select=c(Patient, cluster))
		setkey(tmp, Patient)
		tmp				<- merge( unique(tmp), unique( subset(YX, select=Patient) ), by='Patient' )
		R.clu			<- merge(X.clu, tmp, by='Patient')		
		R.seq			<- merge(X.seq, tmp, by='Patient')		
		tmp				<- subset(clumsm.info, select=c(Patient, cluster))
		setnames(tmp, c('Patient','cluster'), c('t.Patient','t.cluster'))
		setkey(tmp, t.Patient)
		R.clu			<- merge(R.clu, unique(tmp), by='t.Patient')	
		R.clu			<- subset(R.clu, cluster==t.cluster)
		cat(paste('\ncompute risk.table'))
		risk.table		<- do.call('rbind',list(
						risk.df[,	{
									z	<- table( YX[, risk, with=FALSE])
									list(factor=rownames(z), n=as.numeric(unclass(z)), stat='YX')												
								},by='risk'],
						risk.df[,	{
									z	<- table( R.clu[, risk, with=FALSE])
									list(factor=rownames(z), n=as.numeric(unclass(z)), stat='R.clu')												
								},by='risk'],
						risk.df[,	{
									z	<- table( R.seq[, risk, with=FALSE])
									list(factor=rownames(z), n=as.numeric(unclass(z)), stat='R.seq')												
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
		with.Xmsmetc, resume, verbose
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
	if(with.Xmsmetc && grepl('.clu',method.risk,fixed=1))
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
	if(with.Xmsmetc)
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
	if(with.Xmsmetc & ( grepl('adj',method.risk) | grepl('cens',method.risk)))
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
		if(with.Xmsmetc)
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
		
		YX					<- stratificationmodel.Age_253045.Stage_UAE_UAC_UC_D_TS_TO_F(YX)
		if(with.Xmsmetc)
		{
			X.clu			<- stratificationmodel.Age_253045.Stage_UAE_UAC_UC_D_TS_TO_F(X.clu)
			gc()
			X.seq			<- stratificationmodel.Age_253045.Stage_UAE_UAC_UC_D_TS_TO_F(X.seq)
			gc()
			X.msm			<- stratificationmodel.Age_253045.Stage_UAE_UAC_UC_D_TS_TO_F(X.msm)
			gc()
		}
	}
	if(grepl('m5C',method.risk))	
	{
		YX					<- stratificationmodel.Age_253045.Stage_UA_UC_D_TS_TO_F(YX)
		if(with.Xmsmetc)
		{
			X.clu			<- stratificationmodel.Age_253045.Stage_UA_UC_D_TS_TO_F(X.clu)
			gc()
			X.seq			<- stratificationmodel.Age_253045.Stage_UA_UC_D_TS_TO_F(X.seq)
			gc()
			X.msm			<- stratificationmodel.Age_253045.Stage_UA_UC_D_TS_TO_F(X.msm)
			gc()
		}
	}
	if(grepl('m5D',method.risk))	
	{		
		YX					<- stratificationmodel.Age_2833384348.Stage_UAE_UAC_UC_D_TS_TO_F(YX)
		if(with.Xmsmetc)
		{
			X.clu			<- stratificationmodel.Age_2833384348.Stage_UAE_UAC_UC_D_TS_TO_F(X.clu)
			gc()
			X.seq			<- stratificationmodel.Age_2833384348.Stage_UAE_UAC_UC_D_TS_TO_F(X.seq)
			gc()
			X.msm			<- stratificationmodel.Age_2833384348.Stage_UAE_UAC_UC_D_TS_TO_F(X.msm)
			gc()
		}
	}
	if(grepl('m5E',method.risk))	
	{		
		YX					<- stratificationmodel.Age_25354555.Stage_UAE_UAC_UC_D_TS_TO_F(YX)
		if(with.Xmsmetc)
		{
			X.clu			<- stratificationmodel.Age_25354555.Stage_UAE_UAC_UC_D_TS_TO_F(X.clu)
			gc()
			X.seq			<- stratificationmodel.Age_25354555.Stage_UAE_UAC_UC_D_TS_TO_F(X.seq)
			gc()
			X.msm			<- stratificationmodel.Age_25354555.Stage_UAE_UAC_UC_D_TS_TO_F(X.msm)
			gc()
		}
	}
	if(grepl('m5F',method.risk))	
	{		
		YX					<- stratificationmodel.Gen_55657585.Stage_UAE_UAC_UC_D_TS_TO_F(YX)
		if(with.Xmsmetc)
		{
			X.clu			<- stratificationmodel.Gen_55657585.Stage_UAE_UAC_UC_D_TS_TO_F(X.clu)
			gc()
			X.seq			<- stratificationmodel.Gen_55657585.Stage_UAE_UAC_UC_D_TS_TO_F(X.seq)
			gc()
			X.msm			<- stratificationmodel.Gen_55657585.Stage_UAE_UAC_UC_D_TS_TO_F(X.msm)
			gc()
		}
	}
		
	save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', method, 'STRAT_',gsub('\\.clu\\.adj','',gsub('\\.tp[0-9]','',method.risk)),'.R',sep='')
	save(YX, X.clu, X.seq, X.msm, file=save.file)
#STOP1		
#stop()
	#
	#	compute sampling and censoring tables/models that are needed for adjustments
	#
	if(grepl('adj',method.risk) & grepl('clu',method.risk))
	{
			#	get args
			save.file		<- tmp	<- NA
			resume			<- 0						
			if(grepl('m5A',method.risk))		
				tmp			<- 'm5A'
			if(grepl('m5B',method.risk))		
				tmp			<- 'm5B'
			if(grepl('m5C',method.risk))		
				tmp			<- 'm5C'
			if(grepl('m5D',method.risk))		
				tmp			<- 'm5D'
			if(grepl('m5E',method.risk))		
				tmp			<- 'm5E'						
			if(grepl('m5F',method.risk))		
				tmp			<- 'm5F'									
			if(is.na(tmp))	
				stop('unknown method.risk')
			#	sampling tables
			args			<- sampling.get.all.tables.args(method.risk)
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore', method,'_Stables',method.PDT,'_',tmp,'.R',sep='')			
			Stab			<- sampling.get.all.tables(YX, X.seq, X.msm, X.clu, clumsm.info, tperiod.info=tperiod.info, resume=resume, save.file=save.file, method=method.risk, risk.col=args['risk.col'], risktp.col=args['risktp.col'], factor.ref.v=args['factor.ref.v'])
			#	sampling model
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore', method,'_Smodel',method.PDT,'_',tmp,'.R',sep='')
			sm				<- sampling.model.calculate(X.msm, df.all.allmsm, df.viro.allmsm, df.immu.allmsm, df.treatment.allmsm, t.period=t.period, t.endctime=t.endctime, method.lRNA.supp=method.lRNA.supp, method.lRNA.nsupp=4, contact.grace=1.5, resume=resume, save.file=save.file)			
			#	censoring tables
			#save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore', method,'_Ctables',method.PDT,'_',tmp,'.R',sep='')			
			#Ctab			<- censoring.get.all.tables.by.trinterval(YX, X.seq, X.msm, X.clu, tperiod.info=tperiod.info, resume=resume, save.file=save.file, method=method.risk, risk.col=args['risk.col'], risktp.col=args['risktp.col'], factor.ref.v=args['factor.ref.v'], c.period=t.period, bs.n=1000, bs.cdelta.min=2, bs.cdelta.max=3)
			#	censoring model
			save.file		<- paste(outdir,'/',outfile, '_', gsub('/',':',insignat), '_', 'Yscore', method,'_Cmodel',method.PDT,'_',tmp,'.R',sep='')
			risk.col		<- censoring.model.calculate.bs.args(method.risk)
			cm				<- censoring.model.calculate.bs(X.msm, df.all.allmsm, resume=resume, save.file=save.file, t.recent.endctime=t.recent.endctime, risk.col=risk.col, c.period=t.period, c.smpl.n=50, bs.n=100, bs.cdelta.min=2, bs.cdelta.max=3)			
			stop()
	}
	X.clu<- X.seq<- X.msm<- NULL
	gc()
stop()	
	ans		<- list(predict.t2inf=predict.t2inf, t2inf.args=t2inf.args, df.all=df.all, YX=YX, Y.brl.bs=Y.brl.bs)
	ans
}
######################################################################################
stratificationmodel.Age_253045.Stage_UAE_UAC_UC_D_TS_TO_F<- function(YX.m5, lRNA.supp=2)
{
	#	get age stages
	#YX.m5	<- copy(YX)
	if('U.score'%in%colnames(YX.m5))
		YX.m5[, U.score:=NULL]
	cat(paste('\nsubset to save mem\n'))
	if('score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, FASTASampleCode, t.FASTASampleCode, score.Y, telapsed, brl, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, w, w.i, w.in, w.t, w.tn, t.Age, Age, lRNA_T1.supp, lRNA, nlRNA.supp, nlRNA.nsupp ))	
	if(!'score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, t.Age, Age, lRNA_T1.supp, lRNA, nlRNA.supp, nlRNA.nsupp ))
	gc()	
	if('score.Y'%in%colnames(YX.m5) && YX.m5[, !any(score.Y>1.1)])
	{
		tmp	<- YX.m5[, score.Y>0]
		set(YX.m5, which(tmp), 'score.Y', YX.m5[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )
	}	
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
	set(YX.m5, YX.m5[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stageC', 'UAC' )
	#	estimated recent
	set(YX.m5, YX.m5[, which(stage=='U' & (t-t.InfT)<1 & (t.isAcute!='Yes' | is.na(t.isAcute)))], 'stageC', 'UAE')
	#	set Lost
	stopifnot( nrow(subset(YX.m5, stage!='U' & is.na(contact)))==0 )	
	set(YX.m5, YX.m5[, which(stage!='U' & contact!='Yes')], 'stageC', 'L')
	#	break up treated
	YX.m5	<- merge(YX.m5, YX.m5[, {
					lRNA.c3			<- 'ART.vlNA'
					tmp				<- which(!is.na(lRNA))
					if( any(t<lRNA_T1.supp) | any(is.na(lRNA_T1.supp)))
						lRNA.c3		<- 'ART.NotYetFirstSu'
					if( all(t>=lRNA_T1.supp) && length(tmp) && max(lRNA[tmp])>lRNA.supp)
					{
						tmp2		<- sum(nlRNA.nsupp, na.rm=TRUE)
						if(tmp2<=0)
							lRNA.c3	<- 'ART.vlNA'
						if(tmp2>0)
							lRNA.c3	<- 'ART.suA.N'	
					}
					if( all(t>=lRNA_T1.supp) && length(tmp) &&  max(lRNA[tmp])<=lRNA.supp)
					{
						tmp2		<- sum(nlRNA.supp, na.rm=TRUE)
						if(tmp2<=0)
							lRNA.c3	<- 'ART.vlNA'
						if(tmp2>=1)
							lRNA.c3	<- 'ART.suA.Y'																					
					} 						
					list( lRNA.c3=lRNA.c3 )
				}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)
	set(YX.m5, YX.m5[, which(stageC=='T' & lRNA.c3=='ART.suA.Y')], 'stageC', 'TS')
	set(YX.m5, YX.m5[, which(stageC=='T' & lRNA.c3!='ART.suA.Y')], 'stageC', 'TO')
	set(YX.m5, NULL, 'stageC', YX.m5[, factor(stageC, levels=c('UAC','UAE','UC','D','TO','TS','L'))])
	#
	#	PREPARE AGE (Age and t.Age are ages at midpoint of infection interval)
	#
	#	set age group of transmitter				
	age.breaks		<- c(-1, 25, 30, 45, 100)
	tmp				<- c("(-1,25]","(25,35]","(35,45]","(45,55]","(55,100]")
	YX.m5[, t.AgeC:= YX.m5[, cut(t.Age, breaks=age.breaks, labels=tmp)]]	
	#	set age group of recipient
	YX.m5[, AgeC:= YX.m5[, cut(Age, breaks=age.breaks, labels=tmp)]]	
	#	set age group of pair
	YX.m5[, p.AgeC:= YX.m5[, paste(as.character(t.AgeC),'->',as.character(AgeC),sep='')]]
	#	factors
	set(YX.m5, NULL, 't.AgeC', YX.m5[, factor(as.character(t.AgeC), levels=tmp, labels=tmp)])
	set(YX.m5, NULL, 'AgeC', YX.m5[, factor(as.character(AgeC), levels=tmp, labels=tmp)])
	set(YX.m5, NULL, 'p.AgeC', YX.m5[, factor(p.AgeC)])
	#
	#	PREPARE STAGE:AGE
	#		
	tmp				<- as.vector(sapply( c("UAC","UAE","UC","D","TO","TS","L"), function(x)	paste(x,'_',tmp,sep='')))	
	YX.m5[, t.stAgeC:= YX.m5[, factor(paste(as.character(stageC),'_',as.character(t.AgeC),sep=''),labels=tmp, levels=tmp)]]
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
stratificationmodel.Age_2833384348.Stage_UAE_UAC_UC_D_TS_TO_F<- function(YX.m5, lRNA.supp=2)
{
	#	get age stages
	#YX.m5	<- copy(YX)
	if('U.score'%in%colnames(YX.m5))
		YX.m5[, U.score:=NULL]
	cat(paste('\nsubset to save mem\n'))
	if('score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, FASTASampleCode, t.FASTASampleCode, score.Y, telapsed, brl, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, w, w.i, w.in, w.t, w.tn, t.Age, Age, lRNA_T1.supp, lRNA, nlRNA.supp, nlRNA.nsupp ))	
	if(!'score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, t.Age, Age, lRNA_T1.supp, lRNA, nlRNA.supp, nlRNA.nsupp ))
	gc()	
	if('score.Y'%in%colnames(YX.m5) && YX.m5[, !any(score.Y>1.1)])
	{
		tmp	<- YX.m5[, score.Y>0]
		set(YX.m5, which(tmp), 'score.Y', YX.m5[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )
	}	
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
	set(YX.m5, YX.m5[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stageC', 'UAC' )
	#	estimated recent
	set(YX.m5, YX.m5[, which(stage=='U' & (t-t.InfT)<1 & (t.isAcute!='Yes' | is.na(t.isAcute)))], 'stageC', 'UAE')
	#	set Lost
	stopifnot( nrow(subset(YX.m5, stage!='U' & is.na(contact)))==0 )	
	set(YX.m5, YX.m5[, which(stage!='U' & contact!='Yes')], 'stageC', 'L')
	#	break up treated
	YX.m5	<- merge(YX.m5, YX.m5[, {
						lRNA.c3			<- 'ART.vlNA'
						tmp				<- which(!is.na(lRNA))
						if( any(t<lRNA_T1.supp) | any(is.na(lRNA_T1.supp)))
							lRNA.c3		<- 'ART.NotYetFirstSu'
						if( all(t>=lRNA_T1.supp) && length(tmp) && max(lRNA[tmp])>lRNA.supp)
						{
							tmp2		<- sum(nlRNA.nsupp, na.rm=TRUE)
							if(tmp2<=0)
								lRNA.c3	<- 'ART.vlNA'
							if(tmp2>0)
								lRNA.c3	<- 'ART.suA.N'	
						}
						if( all(t>=lRNA_T1.supp) && length(tmp) &&  max(lRNA[tmp])<=lRNA.supp)
						{
							tmp2		<- sum(nlRNA.supp, na.rm=TRUE)
							if(tmp2<=0)
								lRNA.c3	<- 'ART.vlNA'
							if(tmp2>=1)
								lRNA.c3	<- 'ART.suA.Y'																					
						} 						
						list( lRNA.c3=lRNA.c3 )
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)
	set(YX.m5, YX.m5[, which(stageC=='T' & lRNA.c3=='ART.suA.Y')], 'stageC', 'TS')
	set(YX.m5, YX.m5[, which(stageC=='T' & lRNA.c3!='ART.suA.Y')], 'stageC', 'TO')
	set(YX.m5, NULL, 'stageC', YX.m5[, factor(stageC, levels=c('UAC','UAE','UC','D','TO','TS','L'))])
	#
	#	PREPARE AGE (Age and t.Age are ages at midpoint of infection interval)
	#
	#	set age group of transmitter				
	age.breaks		<- c(-1, 28, 33, 38, 43, 48, 100)
	age.labels		<- c('<28','<33','<38','<43','<48','<100')	
	YX.m5[, t.AgeC:= YX.m5[, cut(t.Age, breaks=age.breaks)]]	
	#	set age group of recipient
	YX.m5[, AgeC:= YX.m5[, cut(Age, breaks=age.breaks)]]	
	#	set age group of pair
	YX.m5[, p.AgeC:= YX.m5[, paste(as.character(t.AgeC),'->',as.character(AgeC),sep='')]]
	#	factors
	tmp				<- c("(-1,28]","(28,33]","(33,38]","(38,43]","(43,48]","(48,100]")
	set(YX.m5, NULL, 't.AgeC', YX.m5[, factor(as.character(t.AgeC), levels=tmp, labels=tmp)])
	set(YX.m5, NULL, 'AgeC', YX.m5[, factor(as.character(AgeC), levels=tmp, labels=tmp)])
	set(YX.m5, NULL, 'p.AgeC', YX.m5[, factor(p.AgeC)])
	#
	#	PREPARE STAGE:AGE
	#
	tmp				<- as.vector(sapply( c("UAC","UAE","UC","D","TO","TS","L"), function(x)	paste(x,'_',tmp,sep='')))
	YX.m5[, t.stAgeC:= YX.m5[, factor(paste(as.character(stageC),'_',as.character(t.AgeC),sep=''),labels=tmp, levels=tmp)]]
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
stratificationmodel.Age_25354555.Stage_UAE_UAC_UC_D_TS_TO_F<- function(YX.m5, lRNA.supp=2)
{
	#	get age stages
	#YX.m5	<- copy(YX)
	if('U.score'%in%colnames(YX.m5))
		YX.m5[, U.score:=NULL]
	cat(paste('\nsubset to save mem\n'))
	if('score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, FASTASampleCode, t.FASTASampleCode, score.Y, telapsed, brl, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, w, w.i, w.in, w.t, w.tn, t.Age, Age, lRNA_T1.supp, lRNA, nlRNA.supp, nlRNA.nsupp ))	
	if(!'score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, t.Age, Age, lRNA_T1.supp, lRNA, nlRNA.supp, nlRNA.nsupp ))
	gc()	
	if('score.Y'%in%colnames(YX.m5) && YX.m5[, !any(score.Y>1.1)])
	{
		tmp	<- YX.m5[, score.Y>0]
		set(YX.m5, which(tmp), 'score.Y', YX.m5[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )
	}	
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
	set(YX.m5, YX.m5[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stageC', 'UAC' )
	#	estimated recent
	set(YX.m5, YX.m5[, which(stage=='U' & (t-t.InfT)<1 & (t.isAcute!='Yes' | is.na(t.isAcute)))], 'stageC', 'UAE')
	#	set Lost
	stopifnot( nrow(subset(YX.m5, stage!='U' & is.na(contact)))==0 )	
	set(YX.m5, YX.m5[, which(stage!='U' & contact!='Yes')], 'stageC', 'L')
	#	break up treated
	YX.m5	<- merge(YX.m5, YX.m5[, {
						lRNA.c3			<- 'ART.vlNA'
						tmp				<- which(!is.na(lRNA))
						if( any(t<lRNA_T1.supp) | any(is.na(lRNA_T1.supp)))
							lRNA.c3		<- 'ART.NotYetFirstSu'
						if( all(t>=lRNA_T1.supp) && length(tmp) && max(lRNA[tmp])>lRNA.supp)
						{
							tmp2		<- sum(nlRNA.nsupp, na.rm=TRUE)
							if(tmp2<=0)
								lRNA.c3	<- 'ART.vlNA'
							if(tmp2>0)
								lRNA.c3	<- 'ART.suA.N'	
						}
						if( all(t>=lRNA_T1.supp) && length(tmp) &&  max(lRNA[tmp])<=lRNA.supp)
						{
							tmp2		<- sum(nlRNA.supp, na.rm=TRUE)
							if(tmp2<=0)
								lRNA.c3	<- 'ART.vlNA'
							if(tmp2>=1)
								lRNA.c3	<- 'ART.suA.Y'																					
						} 						
						list( lRNA.c3=lRNA.c3 )
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)
	set(YX.m5, YX.m5[, which(stageC=='T' & lRNA.c3=='ART.suA.Y')], 'stageC', 'TS')
	set(YX.m5, YX.m5[, which(stageC=='T' & lRNA.c3!='ART.suA.Y')], 'stageC', 'TO')
	set(YX.m5, NULL, 'stageC', YX.m5[, factor(stageC, levels=c('UAC','UAE','UC','D','TO','TS','L'))])
	#
	#	PREPARE AGE (Age and t.Age are ages at midpoint of infection interval)
	#
	#	set age group of transmitter				
	age.breaks		<- c(-1, 25, 35, 45, 55, 100)		
	YX.m5[, t.AgeC:= YX.m5[, cut(t.Age, breaks=age.breaks)]]	
	#	set age group of recipient
	YX.m5[, AgeC:= YX.m5[, cut(Age, breaks=age.breaks)]]	
	#	set age group of pair
	YX.m5[, p.AgeC:= YX.m5[, paste(as.character(t.AgeC),'->',as.character(AgeC),sep='')]]
	#	factors
	tmp				<- c("(-1,25]","(25,35]","(35,45]","(45,55]","(55,100]")
	set(YX.m5, NULL, 't.AgeC', YX.m5[, factor(as.character(t.AgeC), levels=tmp, labels=tmp)])
	set(YX.m5, NULL, 'AgeC', YX.m5[, factor(as.character(AgeC), levels=tmp, labels=tmp)])
	set(YX.m5, NULL, 'p.AgeC', YX.m5[, factor(p.AgeC)])
	#
	#	PREPARE STAGE:AGE
	#		
	tmp				<- as.vector(sapply( c("UAC","UAE","UC","D","TO","TS","L"), function(x)	paste(x,'_',tmp,sep='')))
	YX.m5[, t.stAgeC:= YX.m5[, factor(paste(as.character(stageC),'_',as.character(t.AgeC),sep=''),labels=tmp, levels=tmp)]]
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
stratificationmodel.Gen_55657585.Stage_UAE_UAC_UC_D_TS_TO_F<- function(YX.m5, lRNA.supp=2)
{
	#	get age stages
	#YX.m5	<- copy(YX)
	if('U.score'%in%colnames(YX.m5))
		YX.m5[, U.score:=NULL]
	cat(paste('\nsubset to save mem\n'))
	if('score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, FASTASampleCode, t.FASTASampleCode, score.Y, telapsed, brl, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, w, w.i, w.in, w.t, w.tn, t.Age, Age, lRNA_T1.supp, lRNA, nlRNA.supp, nlRNA.nsupp ))	
	if(!'score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, t.Age, Age, lRNA_T1.supp, lRNA, nlRNA.supp, nlRNA.nsupp ))
	gc()	
	if('score.Y'%in%colnames(YX.m5) && YX.m5[, !any(score.Y>1.1)])
	{
		tmp	<- YX.m5[, score.Y>0]
		set(YX.m5, which(tmp), 'score.Y', YX.m5[tmp,(score.Y*(length(tmp)-1)+0.5)/length(tmp)] )
	}	
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
	set(YX.m5, YX.m5[, which(t.isAcute=='Yes' & stage=='Diag' & t-t.AnyPos_T1 < 0.25)], 'stageC', 'UAC' )
	#	estimated recent
	set(YX.m5, YX.m5[, which(stage=='U' & (t-t.InfT)<1 & (t.isAcute!='Yes' | is.na(t.isAcute)))], 'stageC', 'UAE')
	#	set Lost
	stopifnot( nrow(subset(YX.m5, stage!='U' & is.na(contact)))==0 )	
	set(YX.m5, YX.m5[, which(stage!='U' & contact!='Yes')], 'stageC', 'L')
	#	break up treated
	YX.m5	<- merge(YX.m5, YX.m5[, {
						lRNA.c3			<- 'ART.vlNA'
						tmp				<- which(!is.na(lRNA))
						if( any(t<lRNA_T1.supp) | any(is.na(lRNA_T1.supp)))
							lRNA.c3		<- 'ART.NotYetFirstSu'
						if( all(t>=lRNA_T1.supp) && length(tmp) && max(lRNA[tmp])>lRNA.supp)
						{
							tmp2		<- sum(nlRNA.nsupp, na.rm=TRUE)
							if(tmp2<=0)
								lRNA.c3	<- 'ART.vlNA'
							if(tmp2>0)
								lRNA.c3	<- 'ART.suA.N'	
						}
						if( all(t>=lRNA_T1.supp) && length(tmp) &&  max(lRNA[tmp])<=lRNA.supp)
						{
							tmp2		<- sum(nlRNA.supp, na.rm=TRUE)
							if(tmp2<=0)
								lRNA.c3	<- 'ART.vlNA'
							if(tmp2>=1)
								lRNA.c3	<- 'ART.suA.Y'																					
						} 						
						list( lRNA.c3=lRNA.c3 )
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)
	set(YX.m5, YX.m5[, which(stageC=='T' & lRNA.c3=='ART.suA.Y')], 'stageC', 'TS')
	set(YX.m5, YX.m5[, which(stageC=='T' & lRNA.c3!='ART.suA.Y')], 'stageC', 'TO')
	set(YX.m5, NULL, 'stageC', YX.m5[, factor(stageC, levels=c('UAC','UAE','UC','D','TO','TS','L'))])
	#
	#	PREPARE AGE (Age and t.Age are ages at midpoint of infection interval)
	#
	#	set age group of transmitter				
	age.breaks		<- c(-Inf, 1955,1965,1975,1985,Inf)
	age.labels		<- c('1920-55','1956-65','1966-75','1976-85','1985-95')
	YX.m5[, t.AgeC:= YX.m5[, cut(t-t.Age, breaks=age.breaks, labels=age.labels)]]	
	#	set age group of recipient
	YX.m5[, AgeC:= YX.m5[, cut(t-Age, breaks=age.breaks, labels=age.labels)]]	
	#	set age group of pair
	YX.m5[, p.AgeC:= YX.m5[, paste(as.character(t.AgeC),'->',as.character(AgeC),sep='')]]
	#	factors	
	set(YX.m5, NULL, 't.AgeC', YX.m5[, factor(as.character(t.AgeC), levels=age.labels, labels=age.labels)])
	set(YX.m5, NULL, 'AgeC', YX.m5[, factor(as.character(AgeC), levels=age.labels, labels=age.labels)])
	set(YX.m5, NULL, 'p.AgeC', YX.m5[, factor(p.AgeC)])
	#
	#	PREPARE STAGE:AGE
	#		
	tmp				<- as.vector(sapply( c("UAC","UAE","UC","D","TO","TS","L"), function(x)	paste(x,'_',age.labels,sep='')))
	YX.m5[, t.stAgeC:= YX.m5[, factor(paste(as.character(stageC),'_',as.character(t.AgeC),sep=''),labels=tmp, levels=tmp)]]
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
stratificationmodel.Age_253045.Stage_UA_UC_D_TS_TO_F<- function(YX.m5, lRNA.supp=2)
{
	#	get age stages
	#YX.m5	<- copy(YX)
	if('U.score'%in%colnames(YX.m5))
		YX.m5[, U.score:=NULL]
	cat(paste('\nsubset to save mem\n'))
	if('score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, FASTASampleCode, t.FASTASampleCode, score.Y, telapsed, brl, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, w, w.i, w.in, w.t, w.tn, t.Age, Age, lRNA_T1.supp, lRNA, nlRNA.supp, nlRNA.nsupp ))	
	if(!'score.Y'%in%colnames(YX.m5))
		YX.m5	<- subset(YX.m5, select=c(t, t.Patient, Patient, stage, t.period, t.isAcute, contact, t.AnyT_T1, AnyPos_T1, t.AnyPos_T1, t.InfT, t.Age, Age, lRNA_T1.supp, lRNA, nlRNA.supp, nlRNA.nsupp ))
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
	#	break up treated
	YX.m5	<- merge(YX.m5, YX.m5[, {
						lRNA.c3			<- 'ART.vlNA'
						tmp				<- which(!is.na(lRNA))
						if( any(t<lRNA_T1.supp) | any(is.na(lRNA_T1.supp)))
							lRNA.c3		<- 'ART.NotYetFirstSu'
						if( all(t>=lRNA_T1.supp) && length(tmp) && max(lRNA[tmp])>lRNA.supp)
						{
							tmp2		<- sum(nlRNA.nsupp, na.rm=TRUE)
							if(tmp2<=0)
								lRNA.c3	<- 'ART.vlNA'
							if(tmp2>0)
								lRNA.c3	<- 'ART.suA.N'	
						}
						if( all(t>=lRNA_T1.supp) && length(tmp) &&  max(lRNA[tmp])<=lRNA.supp)
						{
							tmp2		<- sum(nlRNA.supp, na.rm=TRUE)
							if(tmp2<=0)
								lRNA.c3	<- 'ART.vlNA'
							if(tmp2>=1)
								lRNA.c3	<- 'ART.suA.Y'																					
						} 						
						list( lRNA.c3=lRNA.c3 )
					}, by=c('Patient','t.Patient')], by=c('Patient','t.Patient'), all.x=TRUE)
	set(YX.m5, YX.m5[, which(stageC=='T' & lRNA.c3=='ART.suA.Y')], 'stageC', 'TS')
	set(YX.m5, YX.m5[, which(stageC=='T' & lRNA.c3!='ART.suA.Y')], 'stageC', 'TO')
	set(YX.m5, NULL, 'stageC', YX.m5[, factor(stageC, levels=c('UA','UC','D','TO','TS','L'))])
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
######################################################################################
rltvtp.model.time.150901<- function(YXc)
{
	YXr		<- subset(YXc, select=c('t', 'score.Y', 't.stAgeC'))
	#	get model
	wm8		<- gamlss( score.Y~t.stAgeC+t:t.stAgeC-1, sigma.formula=~t.stAgeC, data=YXr, family=GA(mu.link='identity') )
	list(rltvm= wm8, rltvd= YXr)
}
######################################################################################
rltvtp.model.time.150828<- function(YXc)
{
	YXr		<- subset(YXc, t>=2004.5, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))
	#	get rolling mean
	tmp		<- YXr[, list(nt=length(score.Y)), by=c('t.stAgeC','t')]
	setkey(tmp, t.stAgeC, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ]))	), by=c('t.stAgeC')]
	YXr		<- merge(YXr, tmp, by=c('t.stAgeC','t'))			
	setkey(YXr, t.stAgeC, t)
	tmp		<- YXr[, {
				z	<- unique(t)
				list(t=z, score.Y.rm= sapply(seq_along(z), function(i)	mean(score.Y[ which(abs(z[i]-t)<=1) ]))	)	
			}, by=c('t.stAgeC')]
	YXr		<- merge(YXr, tmp, by=c('t','t.stAgeC'))
	#	get model
	wm8		<- gamlss( score.Y~t.stAgeC+t:t.stAgeC-1, sigma.formula=~t.stAgeC, data=YXr, family=GA(mu.link='identity') )
	list(rltvm= wm8, rltvd= YXr)
}
######################################################################################
rltvtp.model.time.150923Gen<- function(YXc)
{		
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))
	set(YXi, YXi[, which(score.Y<1e-15)], 'score.Y', 1e-15)
	YXi		<- merge(YXi,YXi[, list(tmid=mean(t)), by=c('t.Patient','Patient')],by=c('t.Patient','Patient'))
	setkey(YXi, t.Patient, Patient, t.stAgeC, tmid)
	YXi		<- unique(YXi)
	YXi[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXi[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]	
	tmp		<- unique(gsub('TO_.*','TO_1920-95',gsub('TS_.*','TS_1920-95',gsub('L_.*','L_1920-95',gsub('1956-65','1920-65',gsub('1920-55','1920-65',tmp,fixed=1),fixed=1)))))	
	set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('TO_.*','TO_1920-95',gsub('TS_.*','TS_1920-95',gsub('L_.*','L_1920-95',gsub('1956-65','1920-65',gsub('1920-55','1920-65',gsub('UAE|UAC','UA',t.stAgeC2),fixed=1),fixed=1)))), levels=tmp, labels=tmp)])	
	YXi		<- subset(YXi, select=c(tmid, t.stAgeC2, score.Y))
	#	get median models		
	wm9		<- gamlss( score.Y~t.stAgeC2+tmid:t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXi, family=WEI() )	
	
	prd<- function(dat, wm9, YXi)
	{		
		dat[, t.stAgeC2:= t.stAgeC]
		tmp		<- dat[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
		tmp		<- unique(gsub('TO_.*','TO_1920-95',gsub('TS_.*','TS_1920-95',gsub('L_.*','L_1920-95',gsub('1956-65','1920-65',gsub('1920-55','1920-65',tmp,fixed=1),fixed=1)))))	
		set(dat, NULL, 't.stAgeC2', dat[, factor(gsub('TO_.*','TO_1920-95',gsub('TS_.*','TS_1920-95',gsub('L_.*','L_1920-95',gsub('1956-65','1920-65',gsub('1920-55','1920-65',gsub('UAE|UAC','UA',t.stAgeC2),fixed=1),fixed=1)))), levels=tmp, labels=tmp)])	
		ywm9	<- predict(wm9, data=YXi, newdata=subset(dat, select=c(tmid, t.stAgeC2)), type='response', what='mu')*log(2)^(1/predict(wm9, data=YXi, newdata=subset(dat, select=c(tmid, t.stAgeC2)), type='response', what='sigma'))		
		dat[, scoreY.pr:=ywm9]
		dat				
	}
	list(prd=prd, mo=wm9, fit=YXi)
}
######################################################################################
rltvtp.model.time.150922Gen<- function(YXc)
{		
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))
	set(YXi, YXi[, which(score.Y<1e-15)], 'score.Y', 1e-15)
	YXi		<- merge(YXi,YXi[, list(tmid=mean(t)), by=c('t.Patient','Patient')],by=c('t.Patient','Patient'))
	setkey(YXi, t.Patient, Patient, t.stAgeC, tmid)
	YXi		<- unique(YXi)
	YXi[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXi[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]	
	tmp		<- unique(gsub('1956-65','1920-65',gsub('1920-55','1920-65',tmp,fixed=1),fixed=1))
	set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('1956-65','1920-65',gsub('1920-55','1920-65',gsub('UAE|UAC','UA',t.stAgeC2),fixed=1),fixed=1), levels=tmp, labels=tmp)])
	YXi		<- subset(YXi, select=c(tmid, t.stAgeC2, score.Y))
	#	get median models		
	wm9		<- gamlss( score.Y~t.stAgeC2+tmid:t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXi, family=WEI() )	
	
	prd<- function(dat, wm9, YXi)
	{		
		dat[, t.stAgeC2:= t.stAgeC]
		tmp		<- dat[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
		tmp		<- unique(gsub('1956-65','1920-65',gsub('1920-55','1920-65',tmp,fixed=1),fixed=1))
		set(dat, NULL, 't.stAgeC2', dat[, factor(gsub('1956-65','1920-65',gsub('1920-55','1920-65',gsub('UAE|UAC','UA',t.stAgeC2),fixed=1),fixed=1), levels=tmp, labels=tmp)])
		ywm9	<- predict(wm9, data=YXi, newdata=subset(dat, select=c(tmid, t.stAgeC2)), type='response', what='mu')*log(2)^(1/predict(wm9, data=YXi, newdata=subset(dat, select=c(tmid, t.stAgeC2)), type='response', what='sigma'))		
		dat[, scoreY.pr:=ywm9]
		dat				
	}
	list(prd=prd, mo=wm9, fit=YXi)
}
######################################################################################
rltvtp.model.time.150923<- function(YXc)
{		
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))
	set(YXi, YXi[, which(score.Y<1e-15)], 'score.Y', 1e-15)
	YXi		<- merge(YXi,YXi[, list(tmid=mean(t)), by=c('t.Patient','Patient')],by=c('t.Patient','Patient'))
	setkey(YXi, t.Patient, Patient, t.stAgeC, tmid)
	YXi		<- unique(YXi)
	YXi[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXi[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
	tmp		<- unique(gsub('TO_.*','TO_(-1,00]',gsub('TS_.*','TS_(-1,00]',gsub('L_.*','L_(-1,00]',gsub('(33,38]','(33,43]',gsub('(38,43]','(33,43]',gsub('(48,100]','(43,100]',gsub('(43,48]','(43,100]',gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',tmp,fixed=1),fixed=1),fixed=1),fixed=1),fixed=1),fixed=1)))))	
	set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('TO_.*','TO_(-1,00]',gsub('TS_.*','TS_(-1,00]',gsub('L_.*','L_(-1,00]',gsub('(33,38]','(33,43]',gsub('(38,43]','(33,43]',gsub('(48,100]','(43,100]',gsub('(43,48]','(43,100]',gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',gsub('UAE|UAC','UA',t.stAgeC2),fixed=1),fixed=1),fixed=1),fixed=1),fixed=1),fixed=1)))), levels=tmp, labels=tmp)])	
	YXi		<- subset(YXi, select=c(tmid, t.stAgeC2, score.Y))
	#	get median models		
	wm9		<- gamlss( score.Y~t.stAgeC2+tmid:t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXi, family=WEI() )	
	
	prd<- function(dat, wm9, YXi)
	{		
		dat[, t.stAgeC2:= t.stAgeC]
		tmp		<- dat[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
		tmp		<- unique(gsub('TO_.*','TO_(-1,00]',gsub('TS_.*','TS_(-1,00]',gsub('L_.*','L_(-1,00]',gsub('(33,38]','(33,43]',gsub('(38,43]','(33,43]',gsub('(48,100]','(43,100]',gsub('(43,48]','(43,100]',gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',tmp,fixed=1),fixed=1),fixed=1),fixed=1),fixed=1),fixed=1)))))	
		set(dat, NULL, 't.stAgeC2', dat[, factor(gsub('TO_.*','TO_(-1,00]',gsub('TS_.*','TS_(-1,00]',gsub('L_.*','L_(-1,00]',gsub('(33,38]','(33,43]',gsub('(38,43]','(33,43]',gsub('(48,100]','(43,100]',gsub('(43,48]','(43,100]',gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',gsub('UAE|UAC','UA',t.stAgeC2),fixed=1),fixed=1),fixed=1),fixed=1),fixed=1),fixed=1)))), levels=tmp, labels=tmp)])	
		ywm9	<- predict(wm9, data=YXi, newdata=subset(dat, select=c(tmid, t.stAgeC2)), type='response', what='mu')*log(2)^(1/predict(wm9, data=YXi, newdata=subset(dat, select=c(tmid, t.stAgeC2)), type='response', what='sigma'))		
		dat[, scoreY.pr:=ywm9]
		dat				
	}
	list(prd=prd, mo=wm9, fit=YXi)
}
######################################################################################
rltvtp.model.time.150922<- function(YXc)
{		
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))
	set(YXi, YXi[, which(score.Y<1e-15)], 'score.Y', 1e-15)
	YXi		<- merge(YXi,YXi[, list(tmid=mean(t)), by=c('t.Patient','Patient')],by=c('t.Patient','Patient'))
	setkey(YXi, t.Patient, Patient, t.stAgeC, tmid)
	YXi		<- unique(YXi)
	YXi[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXi[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]	
	tmp		<- unique(gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',tmp,fixed=1),fixed=1))
	set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',gsub('UAE|UAC','UA',t.stAgeC2),fixed=1),fixed=1), levels=tmp, labels=tmp)])
	YXi		<- subset(YXi, select=c(tmid, t.stAgeC2, score.Y))
	#	get median models		
	wm9		<- gamlss( score.Y~t.stAgeC2+tmid:t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXi, family=WEI() )	
	
	prd<- function(dat, wm9, YXi)
	{		
		dat[, t.stAgeC2:= t.stAgeC]
		tmp		<- dat[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
		tmp		<- unique(gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',tmp,fixed=1),fixed=1))
		set(dat, NULL, 't.stAgeC2', dat[, factor(gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',gsub('UAE|UAC','UA',t.stAgeC2),fixed=1),fixed=1), levels=tmp, labels=tmp)])
		ywm9	<- predict(wm9, data=YXi, newdata=subset(dat, select=c(tmid, t.stAgeC2)), type='response', what='mu')*log(2)^(1/predict(wm9, data=YXi, newdata=subset(dat, select=c(tmid, t.stAgeC2)), type='response', what='sigma'))		
		dat[, scoreY.pr:=ywm9]
		dat				
	}
	list(prd=prd, mo=wm9, fit=YXi)
}
######################################################################################
rltvtp.model.time.150918<- function(YXc)
{		
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))
	set(YXi, YXi[, which(score.Y<1e-15)], 'score.Y', 1e-15)
	YXi		<- merge(YXi,YXi[, list(tmid=mean(t)), by=c('t.Patient','Patient')],by=c('t.Patient','Patient'))
	setkey(YXi, t.Patient, Patient, t.stAgeC, tmid)
	YXi		<- unique(YXi)
	YXi[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXi[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
	set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('UAE|UAC','UA',t.stAgeC2), levels=tmp, labels=tmp)])
	YXi		<- subset(YXi, select=c(tmid, t.stAgeC2, score.Y))
	#	get median models		
	wm10	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=2):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXi, family=WEI() )
	
	prd<- function(dat, wm10, YXi)
	{		
		dat[, t.stAgeC2:= t.stAgeC]
		tmp		<- dat[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
		set(dat, NULL, 't.stAgeC2', dat[, factor(gsub('UAE|UAC','UA',t.stAgeC2), levels=tmp, labels=tmp)])
		ywm10	<- predict(wm10, data=YXi, newdata=subset(dat, select=c(tmid, t.stAgeC2)), type='response', what='mu')*log(2)^(1/predict(wm10, data=YXi, newdata=subset(dat, select=c(tmid, t.stAgeC2)), type='response', what='sigma'))		
		dat[, scoreY.pr:=ywm10]
		dat				
	}
	list(prd=prd, mo=wm10, fit=YXi)
}
######################################################################################
rltvtp.explore.time.pAgeC.150920<- function(YXc, indir=NA, infile=NA, suffix='')
{
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC','AgeC','p.AgeC'))
	YXi		<- merge(YXi,YXi[, list(tmid=mean(t)), by=c('t.Patient','Patient')],by=c('t.Patient','Patient'))
	setkey(YXi, t.Patient, Patient, t.stAgeC, AgeC, t.AgeC, p.AgeC, tmid)
	YXi		<- unique(YXi)
	
	tmp		<- YXi[, list(score.Y=median(score.Y)), by=c('t.stAgeC','AgeC','stageC','t.AgeC')]
	
	ggplot(tmp, aes(x=AgeC, y=t.AgeC, size=score.Y)) + geom_point() + 
			facet_grid(~stageC) + theme_bw()
	ggplot(YXi, aes(x=AgeC, y=score.Y)) + geom_boxplot() +
			facet_grid(stageC~t.AgeC, scales='free') + theme_bw()
	ggplot(YXi, aes(x=AgeC, y=score.Y)) + geom_boxplot() +
			facet_grid(~t.AgeC, scales='free') + theme_bw()
	
	
	tmp		<- YXi[, list(nt=length(score.Y)), by=c('t.stAgeC','t')]
	setkey(tmp, t.stAgeC, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= as.double(sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ]))), cnt=cumsum(nt)	), by=c('t.stAgeC')]
	YXi		<- merge(YXi, tmp, by=c('t.stAgeC','t'))			
	setkey(YXi, t.stAgeC, t)
	YXi[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXi[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
	set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('UAE|UAC','UA',t.stAgeC2), levels=tmp, labels=tmp)])
	set(YXi, YXi[, which(score.Y<1e-15)], 'score.Y', 1e-15)
	#	get rolling median
	tmp		<- YXi[, {
				#z	<- tmid
				#z2	<- rollapply(score.Y, width=200, FUN=median, align="center", partial=TRUE)
				list(	tmid=tmid, 
						#V.rmd5e2= z2[!duplicated(t)],
						V.rmd12m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=1) ])),						
						V.rmd18m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=1.5) ])),
						V.rmd24m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=2) ]))
				)	
			}, by=c('t.stAgeC2')]	
	setkey(tmp, tmid, t.stAgeC2)
	YXi		<- merge(YXi, unique(tmp), by=c('tmid','t.stAgeC2'))
	#	get median models		
	wm9		<- gamlss( score.Y~t.stAgeC2+tmid:t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXi, family=WEI() )
	ywm9	<- predict(wm9, data=YXi, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXi, type='response', what='sigma'))
	wm10	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=2):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXi, family=WEI() )
	ywm10	<- predict(wm10, data=YXi, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXi, type='response', what='sigma'))
	wm11	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=3):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXi, family=WEI() )
	ywm11	<- predict(wm10, data=YXi, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXi, type='response', what='sigma'))
	
	YXi[, score.Y.wm9:=ywm9]
	YXi[, score.Y.wm10:=ywm10]
	YXi[, score.Y.wm11:=ywm11]
	#
	#	plot
	#	
	tmp2	<- melt(YXi, measure.vars=c('V.rmd12m','V.rmd18m','V.rmd24m'))
	ggplot(tmp2, aes(x=tmid, colour=t.AgeC, group=t.stAgeC)) + 
			geom_point(aes(y=value), size=1.2) +
			geom_line(aes(y=value)) +		
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(variable~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmedianstAgeCByTime_1yrollingmedian',suffix,sep=''),infile),'.pdf',sep='')	
	ggsave(file=file, w=10,h=7)
	#	
	tmp		<- melt(YXi, measure.vars=c('score.Y.wm9','score.Y.wm10','score.Y.wm11'))
	ggplot(tmp, aes(x=tmid, group=t.stAgeC)) + geom_line(aes(y=value, colour=t.AgeC)) +
			geom_point(aes(y=score.Y, colour=t.AgeC), alpha=0.2) +
			geom_point(aes(y=V.rmd12m), pch=2, size=1.2, colour='black') +			
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			scale_y_continuous(limits=c(0,40)) +
			facet_grid(variable+t.AgeC~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmedianstAgeCByTimeModels_1yrollingmedian',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=10,h=20)	
}
######################################################################################
rltvtp.explore.time.150918<- function(YXc, indir=NA, infile=NA, suffix='')
{
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))
	YXi		<- merge(YXi,YXi[, list(tmid=mean(t)), by=c('t.Patient','Patient')],by=c('t.Patient','Patient'))
	setkey(YXi, t.Patient, Patient, t.stAgeC, tmid)
	YXi		<- unique(YXi)
	tmp		<- YXi[, list(nt=length(score.Y)), by=c('t.stAgeC','t')]
	setkey(tmp, t.stAgeC, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= as.double(sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ]))), cnt=cumsum(nt)	), by=c('t.stAgeC')]
	YXi		<- merge(YXi, tmp, by=c('t.stAgeC','t'))			
	setkey(YXi, t.stAgeC, t)
	YXi[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXi[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
	set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('UAE|UAC','UA',t.stAgeC2), levels=tmp, labels=tmp)])
	set(YXi, YXi[, which(score.Y<1e-15)], 'score.Y', 1e-15)
	#	get rolling median
	tmp		<- YXi[, {
				#z	<- tmid
				#z2	<- rollapply(score.Y, width=200, FUN=median, align="center", partial=TRUE)
				list(	tmid=tmid, 
						#V.rmd5e2= z2[!duplicated(t)],
						V.rmd12m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=1) ])),						
						V.rmd18m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=1.5) ])),
						V.rmd24m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=2) ]))
				)	
			}, by=c('t.stAgeC2')]	
	setkey(tmp, tmid, t.stAgeC2)
	YXi		<- merge(YXi, unique(tmp), by=c('tmid','t.stAgeC2'))
	#	get median models		
	wm9		<- gamlss( score.Y~t.stAgeC2+tmid:t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXi, family=WEI() )
	ywm9	<- predict(wm9, data=YXi, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXi, type='response', what='sigma'))
	wm10	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=2):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXi, family=WEI() )
	ywm10	<- predict(wm10, data=YXi, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXi, type='response', what='sigma'))
	wm11	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=3):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXi, family=WEI() )
	ywm11	<- predict(wm10, data=YXi, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXi, type='response', what='sigma'))
	
	YXi[, score.Y.wm9:=ywm9]
	YXi[, score.Y.wm10:=ywm10]
	YXi[, score.Y.wm11:=ywm11]
	#
	#	plot
	#	
	tmp2	<- melt(YXi, measure.vars=c('V.rmd12m','V.rmd18m','V.rmd24m'))
	ggplot(tmp2, aes(x=tmid, colour=t.AgeC, group=t.stAgeC)) + 
			geom_point(aes(y=value), size=1.2) +
			geom_line(aes(y=value)) +		
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(variable~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmedianstAgeCByTime_1yrollingmedian',suffix,sep=''),infile),'.pdf',sep='')	
	ggsave(file=file, w=10,h=7)
	#	
	tmp		<- melt(YXi, measure.vars=c('score.Y.wm9','score.Y.wm10','score.Y.wm11'))
	ggplot(tmp, aes(x=tmid, group=t.stAgeC)) + geom_line(aes(y=value, colour=t.AgeC)) +
			geom_point(aes(y=score.Y, colour=t.AgeC), alpha=0.2) +
			geom_point(aes(y=V.rmd12m), pch=2, size=1.2, colour='black') +			
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			scale_y_continuous(limits=c(0,40)) +
			facet_grid(variable+t.AgeC~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmedianstAgeCByTimeModels_1yrollingmedian',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=10,h=20)	
}
######################################################################################
rltvtp.explore.time.150922Gen<- function(YXc, indir=NA, infile=NA, suffix='')
{
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))
	YXi		<- merge(YXi,YXi[, list(tmid=mean(t)), by=c('t.Patient','Patient')],by=c('t.Patient','Patient'))
	setkey(YXi, t.Patient, Patient, t.stAgeC, tmid)
	YXi		<- unique(YXi)
	tmp		<- YXi[, list(nt=length(score.Y)), by=c('t.stAgeC','t')]
	setkey(tmp, t.stAgeC, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= as.double(sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ]))), cnt=cumsum(nt)	), by=c('t.stAgeC')]
	YXi		<- merge(YXi, tmp, by=c('t.stAgeC','t'))			
	setkey(YXi, t.stAgeC, t)
	YXi[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXi[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
	tmp		<- unique(gsub('1956-65','1920-65',gsub('1920-55','1920-65',tmp,fixed=1),fixed=1))
	set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('1956-65','1920-65',gsub('1920-55','1920-65',gsub('UAE|UAC','UA',t.stAgeC2),fixed=1),fixed=1), levels=tmp, labels=tmp)])
	#set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('UAE|UAC','UA',t.stAgeC2), levels=tmp, labels=tmp)])
	set(YXi, YXi[, which(score.Y<1e-15)], 'score.Y', 1e-15)
	#	get rolling median
	tmp		<- YXi[, {
				#z	<- tmid
				#z2	<- rollapply(score.Y, width=200, FUN=median, align="center", partial=TRUE)
				list(	tmid=tmid, 
						#V.rmd5e2= z2[!duplicated(t)],
						V.rmd12m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=1) ])),						
						V.rmd18m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=1.5) ])),
						V.rmd24m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=2) ]))
				)	
			}, by=c('t.stAgeC2')]	
	setkey(tmp, tmid, t.stAgeC2)
	YXi		<- merge(YXi, unique(tmp), by=c('tmid','t.stAgeC2'))
	YXj		<- subset(YXi, select=c(tmid,  t.AgeC, stageC, t.stAgeC, t.stAgeC2, score.Y))
	YXp		<- as.data.table(expand.grid(tmid=YXi[, sort(unique(tmid))], t.stAgeC2=YXi[, levels(t.stAgeC2)], stringsAsFactors=0))
	YXp		<- merge(YXp, unique(subset(YXi, select=c(t.AgeC, stageC, t.stAgeC, t.stAgeC2))), by='t.stAgeC2', allow.cartesian=1)
	#	get median models		
	wm9		<- gamlss( score.Y~t.stAgeC2+tmid:t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXj, family=WEI() )
	ywm9	<- predict(wm9, data=YXj, newdata=YXp, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXj, newdata=YXp, type='response', what='sigma'))
	wm10	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=2):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXj, family=WEI() )
	ywm10	<- predict(wm10, data=YXj, newdata=YXp, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXj, newdata=YXp, type='response', what='sigma'))
	wm11	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=3):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXj, family=WEI() )
	ywm11	<- predict(wm10, data=YXj, newdata=YXp, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXj, newdata=YXp, type='response', what='sigma'))
	
	YXp[, score.Y.wm9:=ywm9]
	YXp[, score.Y.wm10:=ywm10]
	YXp[, score.Y.wm11:=ywm11]
	#
	#	plot
	#	
	tmp2	<- melt(YXi, measure.vars=c('V.rmd12m','V.rmd18m','V.rmd24m'))
	ggplot(tmp2, aes(x=tmid, colour=t.AgeC, group=t.stAgeC)) + 
			geom_point(aes(y=value), size=1.2) +
			geom_line(aes(y=value)) +		
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(variable~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmedianstAgeCByTime_1yrollingmedian',suffix,sep=''),infile),'.pdf',sep='')	
	ggsave(file=file, w=10,h=7)
	#	
	tmp		<- merge(YXp, YXi, by=c('tmid','t.AgeC','stageC','t.stAgeC','t.stAgeC2'),all=1,allow.cartesian=1)
	tmp		<- melt(tmp, measure.vars=c('score.Y.wm9','score.Y.wm10','score.Y.wm11'))
	ggplot(tmp, aes(x=tmid, group=t.stAgeC)) + geom_line(aes(y=value, colour=t.AgeC)) +
			geom_point(aes(y=score.Y, colour=t.AgeC), alpha=0.2) +
			geom_point(aes(y=V.rmd12m), pch=2, size=1.2, colour='black') +			
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			scale_y_continuous(limits=c(0,40)) +
			facet_grid(variable+t.AgeC~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmedianstAgeCByTimeModels_1yrollingmedian',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=10,h=20)	
}
######################################################################################
rltvtp.explore.time.150923Gen<- function(YXc, indir=NA, infile=NA, suffix='')
{
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))
	YXi		<- merge(YXi,YXi[, list(tmid=mean(t)), by=c('t.Patient','Patient')],by=c('t.Patient','Patient'))
	setkey(YXi, t.Patient, Patient, t.stAgeC, tmid)
	YXi		<- unique(YXi)
	tmp		<- YXi[, list(nt=length(score.Y)), by=c('t.stAgeC','t')]
	setkey(tmp, t.stAgeC, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= as.double(sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ]))), cnt=cumsum(nt)	), by=c('t.stAgeC')]
	YXi		<- merge(YXi, tmp, by=c('t.stAgeC','t'))			
	setkey(YXi, t.stAgeC, t)
	YXi[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXi[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]	
	tmp		<- unique(gsub('TO_.*','TO_1920-95',gsub('TS_.*','TS_1920-95',gsub('L_.*','L_1920-95',gsub('1956-65','1920-65',gsub('1920-55','1920-65',tmp,fixed=1),fixed=1)))))	
	set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('TO_.*','TO_1920-95',gsub('TS_.*','TS_1920-95',gsub('L_.*','L_1920-95',gsub('1956-65','1920-65',gsub('1920-55','1920-65',gsub('UAE|UAC','UA',t.stAgeC2),fixed=1),fixed=1)))), levels=tmp, labels=tmp)])	
	#set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('UAE|UAC','UA',t.stAgeC2), levels=tmp, labels=tmp)])
	set(YXi, YXi[, which(score.Y<1e-15)], 'score.Y', 1e-15)
	#	get rolling median
	tmp		<- YXi[, {
				#z	<- tmid
				#z2	<- rollapply(score.Y, width=200, FUN=median, align="center", partial=TRUE)
				list(	tmid=tmid, 
						#V.rmd5e2= z2[!duplicated(t)],
						V.rmd12m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=1) ])),						
						V.rmd18m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=1.5) ])),
						V.rmd24m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=2) ]))
				)	
			}, by=c('t.stAgeC2')]	
	setkey(tmp, tmid, t.stAgeC2)
	YXi		<- merge(YXi, unique(tmp), by=c('tmid','t.stAgeC2'))
	YXj		<- subset(YXi, select=c(tmid,  t.AgeC, stageC, t.stAgeC, t.stAgeC2, score.Y))
	YXp		<- as.data.table(expand.grid(tmid=YXi[, sort(unique(tmid))], t.stAgeC2=YXi[, levels(t.stAgeC2)], stringsAsFactors=0))
	YXp		<- merge(YXp, unique(subset(YXi, select=c(t.AgeC, stageC, t.stAgeC, t.stAgeC2))), by='t.stAgeC2', allow.cartesian=1)
	#	get median models		
	wm9		<- gamlss( score.Y~t.stAgeC2+tmid:t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXj, family=WEI() )
	ywm9	<- predict(wm9, data=YXj, newdata=YXp, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXj, newdata=YXp, type='response', what='sigma'))
	wm10	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=2):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXj, family=WEI() )
	ywm10	<- predict(wm10, data=YXj, newdata=YXp, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXj, newdata=YXp, type='response', what='sigma'))
	wm11	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=3):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXj, family=WEI() )
	ywm11	<- predict(wm10, data=YXj, newdata=YXp, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXj, newdata=YXp, type='response', what='sigma'))
	
	YXp[, score.Y.wm9:=ywm9]
	YXp[, score.Y.wm10:=ywm10]
	YXp[, score.Y.wm11:=ywm11]
	#
	#	plot
	#	
	tmp2	<- melt(YXi, measure.vars=c('V.rmd12m','V.rmd18m','V.rmd24m'))
	ggplot(tmp2, aes(x=tmid, colour=t.AgeC, group=t.stAgeC)) + 
			geom_point(aes(y=value), size=1.2) +
			geom_line(aes(y=value)) +		
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(variable~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmedianstAgeCByTime_1yrollingmedian',suffix,sep=''),infile),'.pdf',sep='')	
	ggsave(file=file, w=10,h=7)
	#	
	tmp		<- merge(YXp, YXi, by=c('tmid','t.AgeC','stageC','t.stAgeC','t.stAgeC2'),all=1,allow.cartesian=1)
	tmp		<- melt(tmp, measure.vars=c('score.Y.wm9','score.Y.wm10','score.Y.wm11'))
	ggplot(tmp, aes(x=tmid, group=t.stAgeC)) + geom_line(aes(y=value, colour=t.AgeC)) +
			geom_point(aes(y=score.Y, colour=t.AgeC), alpha=0.2) +
			geom_point(aes(y=V.rmd12m), pch=2, size=1.2, colour='black') +			
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			scale_y_continuous(limits=c(0,40)) +
			facet_grid(variable+t.AgeC~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmedianstAgeCByTimeModels_1yrollingmedian',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=10,h=20)	
}
######################################################################################
rltvtp.explore.time.150922<- function(YXc, indir=NA, infile=NA, suffix='')
{
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))
	YXi		<- merge(YXi,YXi[, list(tmid=mean(t)), by=c('t.Patient','Patient')],by=c('t.Patient','Patient'))
	setkey(YXi, t.Patient, Patient, t.stAgeC, tmid)
	YXi		<- unique(YXi)
	tmp		<- YXi[, list(nt=length(score.Y)), by=c('t.stAgeC','t')]
	setkey(tmp, t.stAgeC, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= as.double(sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ]))), cnt=cumsum(nt)	), by=c('t.stAgeC')]
	YXi		<- merge(YXi, tmp, by=c('t.stAgeC','t'))			
	setkey(YXi, t.stAgeC, t)
	YXi[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXi[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
	tmp		<- unique(gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',tmp,fixed=1),fixed=1))
	set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',gsub('UAE|UAC','UA',t.stAgeC2),fixed=1),fixed=1), levels=tmp, labels=tmp)])
	set(YXi, YXi[, which(score.Y<1e-15)], 'score.Y', 1e-15)
	#	get rolling median
	tmp		<- YXi[, {
				#z	<- tmid
				#z2	<- rollapply(score.Y, width=200, FUN=median, align="center", partial=TRUE)
				list(	tmid=tmid, 
						#V.rmd5e2= z2[!duplicated(t)],
						V.rmd12m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=1) ])),						
						V.rmd18m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=1.5) ])),
						V.rmd24m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=2) ]))
				)	
			}, by=c('t.stAgeC2')]	
	setkey(tmp, tmid, t.stAgeC2)
	YXi		<- merge(YXi, unique(tmp), by=c('tmid','t.stAgeC2'))
	YXj		<- subset(YXi, select=c(tmid,  t.AgeC, stageC, t.stAgeC, t.stAgeC2, score.Y))
	YXp		<- as.data.table(expand.grid(tmid=YXi[, sort(unique(tmid))], t.stAgeC2=YXi[, levels(t.stAgeC2)], stringsAsFactors=0))
	YXp		<- merge(YXp, unique(subset(YXi, select=c(t.AgeC, stageC, t.stAgeC, t.stAgeC2))), by='t.stAgeC2', allow.cartesian=1)
	#	get median models		
	wm9		<- gamlss( score.Y~t.stAgeC2+tmid:t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXj, family=WEI() )
	ywm9	<- predict(wm9, data=YXj, newdata=YXp, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXj, newdata=YXp, type='response', what='sigma'))
	wm10	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=2):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXj, family=WEI() )
	ywm10	<- predict(wm10, data=YXj, newdata=YXp, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXj, newdata=YXp, type='response', what='sigma'))
	wm11	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=3):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXj, family=WEI() )
	ywm11	<- predict(wm10, data=YXj, newdata=YXp, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXj, newdata=YXp, type='response', what='sigma'))
	
	YXp[, score.Y.wm9:=ywm9]
	YXp[, score.Y.wm10:=ywm10]
	YXp[, score.Y.wm11:=ywm11]
	#
	#	plot
	#	
	tmp2	<- melt(YXi, measure.vars=c('V.rmd12m','V.rmd18m','V.rmd24m'))
	ggplot(tmp2, aes(x=tmid, colour=t.AgeC, group=t.stAgeC)) + 
			geom_point(aes(y=value), size=1.2) +
			geom_line(aes(y=value)) +		
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(variable~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmedianstAgeCByTime_1yrollingmedian',suffix,sep=''),infile),'.pdf',sep='')	
	ggsave(file=file, w=10,h=7)
	#	
	tmp		<- merge(YXp, YXi, by=c('tmid','t.AgeC','stageC','t.stAgeC','t.stAgeC2'),all=1,allow.cartesian=1)
	tmp		<- melt(tmp, measure.vars=c('score.Y.wm9','score.Y.wm10','score.Y.wm11'))
	ggplot(tmp, aes(x=tmid, group=t.stAgeC)) + geom_line(aes(y=value, colour=t.AgeC)) +
			geom_point(aes(y=score.Y, colour=t.AgeC), alpha=0.2) +
			geom_point(aes(y=V.rmd12m), pch=2, size=1.2, colour='black') +			
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			scale_y_continuous(limits=c(0,40)) +
			facet_grid(variable+t.AgeC~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmedianstAgeCByTimeModels_1yrollingmedian',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=10,h=20)	
}
######################################################################################
rltvtp.explore.time.150923<- function(YXc, indir=NA, infile=NA, suffix='')
{
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))
	YXi		<- merge(YXi,YXi[, list(tmid=mean(t)), by=c('t.Patient','Patient')],by=c('t.Patient','Patient'))
	setkey(YXi, t.Patient, Patient, t.stAgeC, tmid)
	YXi		<- unique(YXi)
	tmp		<- YXi[, list(nt=length(score.Y)), by=c('t.stAgeC','t')]
	setkey(tmp, t.stAgeC, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= as.double(sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ]))), cnt=cumsum(nt)	), by=c('t.stAgeC')]
	YXi		<- merge(YXi, tmp, by=c('t.stAgeC','t'))			
	setkey(YXi, t.stAgeC, t)
	YXi[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXi[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
	tmp		<- unique(gsub('TO_.*','TO_(-1,00]',gsub('TS_.*','TS_(-1,00]',gsub('L_.*','L_(-1,00]',gsub('(33,38]','(33,43]',gsub('(38,43]','(33,43]',gsub('(48,100]','(43,100]',gsub('(43,48]','(43,100]',gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',tmp,fixed=1),fixed=1),fixed=1),fixed=1),fixed=1),fixed=1)))))	
	set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('TO_.*','TO_(-1,00]',gsub('TS_.*','TS_(-1,00]',gsub('L_.*','L_(-1,00]',gsub('(33,38]','(33,43]',gsub('(38,43]','(33,43]',gsub('(48,100]','(43,100]',gsub('(43,48]','(43,100]',gsub('(45,55]','(45,100]',gsub('(55,100]','(45,100]',gsub('UAE|UAC','UA',t.stAgeC2),fixed=1),fixed=1),fixed=1),fixed=1),fixed=1),fixed=1)))), levels=tmp, labels=tmp)])
	set(YXi, YXi[, which(score.Y<1e-15)], 'score.Y', 1e-15)
	#	get rolling median
	tmp		<- YXi[, {
				#z	<- tmid
				#z2	<- rollapply(score.Y, width=200, FUN=median, align="center", partial=TRUE)
				list(	tmid=tmid, 
						#V.rmd5e2= z2[!duplicated(t)],
						V.rmd12m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=1) ])),						
						V.rmd18m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=1.5) ])),
						V.rmd24m= sapply(seq_along(tmid), function(i)	median(score.Y[ which(abs(tmid[i]-tmid)<=2) ]))
				)	
			}, by=c('t.stAgeC2')]	
	setkey(tmp, tmid, t.stAgeC2)
	YXi		<- merge(YXi, unique(tmp), by=c('tmid','t.stAgeC2'))
	YXj		<- subset(YXi, select=c(tmid,  t.AgeC, stageC, t.stAgeC, t.stAgeC2, score.Y))
	YXp		<- as.data.table(expand.grid(tmid=YXi[, sort(unique(tmid))], t.stAgeC2=YXi[, levels(t.stAgeC2)], stringsAsFactors=0))
	YXp		<- merge(YXp, unique(subset(YXi, select=c(t.AgeC, stageC, t.stAgeC, t.stAgeC2))), by='t.stAgeC2', allow.cartesian=1)
	#	get median models		
	wm9		<- gamlss( score.Y~t.stAgeC2+tmid:t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXj, family=WEI() )
	ywm9	<- predict(wm9, data=YXj, newdata=YXp, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXj, newdata=YXp, type='response', what='sigma'))
	wm10	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=2):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXj, family=WEI() )
	ywm10	<- predict(wm10, data=YXj, newdata=YXp, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXj, newdata=YXp, type='response', what='sigma'))
	wm11	<- gamlss( score.Y~t.stAgeC2+ns(tmid, df=3):t.stAgeC2-1, sigma.formula=~t.stAgeC2, data=YXj, family=WEI() )
	ywm11	<- predict(wm10, data=YXj, newdata=YXp, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXj, newdata=YXp, type='response', what='sigma'))
	
	YXp[, score.Y.wm9:=ywm9]
	YXp[, score.Y.wm10:=ywm10]
	YXp[, score.Y.wm11:=ywm11]
	#
	#	plot
	#	
	tmp2	<- melt(YXi, measure.vars=c('V.rmd12m','V.rmd18m','V.rmd24m'))
	ggplot(tmp2, aes(x=tmid, colour=t.AgeC, group=t.stAgeC)) + 
			geom_point(aes(y=value), size=1.2) +
			geom_line(aes(y=value)) +		
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(variable~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmedianstAgeCByTime_1yrollingmedian',suffix,sep=''),infile),'.pdf',sep='')	
	ggsave(file=file, w=10,h=7)
	#	
	tmp		<- merge(YXp, YXi, by=c('tmid','t.AgeC','stageC','t.stAgeC','t.stAgeC2'),all=1,allow.cartesian=1)
	tmp		<- melt(tmp, measure.vars=c('score.Y.wm9','score.Y.wm10','score.Y.wm11'))
	ggplot(tmp, aes(x=tmid, group=t.stAgeC)) + geom_line(aes(y=value, colour=t.AgeC)) +
			geom_point(aes(y=score.Y, colour=t.AgeC), alpha=0.2) +
			geom_point(aes(y=V.rmd12m), pch=2, size=1.2, colour='black') +			
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			scale_y_continuous(limits=c(0,40)) +
			facet_grid(variable+t.AgeC~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmedianstAgeCByTimeModels_1yrollingmedian',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=10,h=20)	
}
######################################################################################
rltvtp.model.time.150916<- function(YXc, indir=NA, infile=NA, suffix='')
{
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))			
	tmp		<- YXi[, list(nt=length(score.Y)), by=c('t.stAgeC','t')]
	setkey(tmp, t.stAgeC, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= as.double(sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ]))), cnt=cumsum(nt)	), by=c('t.stAgeC')]
	YXi		<- merge(YXi, tmp, by=c('t.stAgeC','t'))			
	setkey(YXi, t.stAgeC, t)
	YXi		<- YXi[, {
				z	<- unique(t)
				list(	t=z, 
						V.rmd12m= sapply(seq_along(z), function(i)	median(score.Y[ which(abs(z[i]-t)<=1) ])),
						V.rmd5e2= 				
						#V.rmd18m= sapply(seq_along(z), function(i)	median(score.Y[ which(abs(z[i]-t)<=1.5) ])),
						#V.rmd24m= sapply(seq_along(z), function(i)	median(score.Y[ which(abs(z[i]-t)<=2) ]))
				)	
			}, by=c('t.stAgeC','t.AgeC', 'stageC')]	
	
	setkey(tpds.df, NCNST, isAcute, tm)
	tmp		<- merge(tpds.df, data.table(isAcute=c('Yes','NI'), DUMMY=c(1e2, 2e3)), by='isAcute')
	setkey(tmp, NCNST, isAcute, tm)		
	z		<- tmp[, list( Patient=Patient, r.Patient=r.Patient, tm=tm, NCNS.rm5=rollapply(NCNS, width=DUMMY[1], FUN=mean, align="center", partial=TRUE) ), by=c('NCNST','isAcute')]
	tpds.df	<- merge(tpds.df, z, by=c('Patient','r.Patient','tm','NCNST','isAcute'))
	
	
	prd		 <- function(dat, YXi)
	{
		dat[, {
					z	<- which(YXi[['t.stAgeC']]==t.stAgeC)
					list(t=t, scoreY.pr= approx(x= YXi[['t']][z], y=YXi[['V.rmd12m']][z], xout=t, rule=2, method='linear')$y)			
				}, by=setdiff(colnames(dat),'t')]
	}
	if(!is.na(indir) & !is.na(infile))
	{
		tmp		<- as.data.table(expand.grid(t.stAgeC=YXc[, levels(t.stAgeC)], t=YXc[, sort(unique(t))]))
		tmp		<- prd(tmp, YXi)
		tmp		<- merge(tmp, unique(subset(YXi, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
		ggplot(YXi, aes(x=t, colour=t.AgeC, group=t.stAgeC)) + 
				geom_point(aes(y=V.rmd12m), size=1.2) + 
				geom_line(data=tmp, aes(y=scoreY.pr)) +
				scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
				facet_grid(~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
				theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
		file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmeanstAgeCByTime_1yrollingmedian',suffix,sep=''),infile),'.pdf',sep='')	
		ggsave(file=file, w=10,h=7)	
	}
	
	list(predict=prd, predict.args=YXi)
}
######################################################################################
rltvtp.model.time.150917<- function(YXc, indir=NA, infile=NA, suffix='')
{
	YXi		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','stageC','t.AgeC', 't.stAgeC'))			
	tmp		<- YXi[, list(nt=length(score.Y)), by=c('t.stAgeC','t')]
	setkey(tmp, t.stAgeC, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= as.double(sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ]))), cnt=cumsum(nt)	), by=c('t.stAgeC')]
	YXi		<- merge(YXi, tmp, by=c('t.stAgeC','t'))			
	setkey(YXi, t.stAgeC, t)
	YXi[, t.stAgeC2:= t.stAgeC]
	tmp		<- YXi[, unique(gsub('UAE|UAC','UA',levels(t.stAgeC)))]
	set(YXi, NULL, 't.stAgeC2', YXi[, factor(gsub('UAE|UAC','UA',t.stAgeC2), levels=tmp, labels=tmp)])	
	tmp		<- YXi[, {
				z	<- unique(t)
				list(	t=z, 
						V.rmd12m= sapply(seq_along(z), function(i)	median(score.Y[ which(abs(z[i]-t)<=1) ]))
				#V.rmd18m= sapply(seq_along(z), function(i)	median(score.Y[ which(abs(z[i]-t)<=1.5) ])),
				#V.rmd24m= sapply(seq_along(z), function(i)	median(score.Y[ which(abs(z[i]-t)<=2) ]))
				)	
			}, by=c('t.stAgeC2')]	
	YXi		<- subset(tmp, !grepl('UA',t.stAgeC2))
	tmp		<- subset(tmp, grepl('UA',t.stAgeC2))
	set(tmp, NULL, 't.stAgeC2', tmp[, gsub('UA','UAC', t.stAgeC2)])
	YXi		<- rbind(YXi, tmp)
	set(tmp, NULL, 't.stAgeC2', tmp[, gsub('UAC','UAE', t.stAgeC2)])
	YXi		<- rbind(YXi, tmp)
	setnames(YXi, 't.stAgeC2', 't.stAgeC')
	set(YXi, NULL, 't.stAgeC', YXi[, as.character(t.stAgeC)]) 
	set(YXi, NULL, 't.AgeC', YXi[, sapply(strsplit(t.stAgeC,'_'),'[[',2)])
	set(YXi, NULL, 'stageC', YXi[, sapply(strsplit(t.stAgeC,'_'),'[[',1)])
	set(YXi, NULL, 'stageC', YXi[, factor(stageC, levels=YXc[, levels(stageC)],labels=YXc[, levels(stageC)])])
	set(YXi, NULL, 't.AgeC', YXi[, factor(t.AgeC, levels=YXc[, levels(t.AgeC)],labels=YXc[, levels(t.AgeC)])])
	set(YXi, NULL, 't.stAgeC', YXi[, factor(t.stAgeC, levels=YXc[, levels(t.stAgeC)],labels=YXc[, levels(t.stAgeC)])])
	
	prd		 <- function(dat, YXi)
	{
		dat[, {
					z	<- which(YXi[['t.stAgeC']]==t.stAgeC)
					list(t=t, scoreY.pr= approx(x= YXi[['t']][z], y=YXi[['V.rmd12m']][z], xout=t, rule=2, method='linear')$y)			
				}, by=setdiff(colnames(dat),'t')]
	}
	if(!is.na(indir) & !is.na(infile))
	{
		tmp		<- as.data.table(expand.grid(t.stAgeC=YXc[, levels(t.stAgeC)], t=YXc[, sort(unique(t))]))
		tmp		<- prd(tmp, YXi)
		tmp		<- merge(tmp, unique(subset(YXi, select=c(stageC, t.AgeC, t.stAgeC))), by='t.stAgeC')
		ggplot(YXi, aes(x=t, colour=t.AgeC, group=t.stAgeC)) + 
				geom_point(aes(y=V.rmd12m), size=1.2) + 
				geom_line(data=tmp, aes(y=scoreY.pr)) +
				scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
				facet_grid(~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
				theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
		file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmeanstAgeCByTime_1yrollingmedian',suffix,sep=''),infile),'.pdf',sep='')	
		ggsave(file=file, w=10,h=7)	
	}	
	list(predict=prd, predict.args=YXi)
}
######################################################################################
rltvtp.explore.time.150916<- function(YXc, indir=NA, infile=NA, suffix='')
{
	#	it s very difficult to look at anything before 2003
	#	before 2004.5, the GA ID link model leads to negative values.. 		
	YXr		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','ntPatient','stageC','t.AgeC', 't.stAgeC'))
	YXr[, ntPatientn:= as.numeric(as.character(ntPatient))]
	YXp		<- copy(YXr)
	setnames(YXp, 't.stAgeC','RSKF')
	tmp		<- YXp[, list(nt=length(score.Y)), by=c('RSKF','t')]
	setkey(tmp, RSKF, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= as.double(sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ]))), cnt=cumsum(nt)	), by=c('RSKF')]
	YXp		<- merge(YXp, tmp, by=c('RSKF','t'))		
	YXp		<- melt( YXp, measure.vars=c('score.Y'), variable.name='STAT', value.name='V' )
	setkey(YXp, STAT, RSKF, t)
	#tmp		<- YXp[, list(V=mean(V), nt=nt[1], cnt=cnt[1], wnt=wnt[1]), by=c('t','RSKF','STAT')]
	#tmp		<- tmp[, list(t=t, V.rm=rollapply(V, width=ceiling(max(cnt)/100), FUN=mean, align="center", partial=TRUE)), by=c('STAT','RSKF')]	
	tmp		<- YXp[, {
				z	<- unique(t)
				list(	t=z, 
						V.rm12m= sapply(seq_along(z), function(i)	mean(V[ which(abs(z[i]-t)<=1) ])),
						V.rmd12m= sapply(seq_along(z), function(i)	median(V[ which(abs(z[i]-t)<=1) ])),
						V.rmd18m= sapply(seq_along(z), function(i)	median(V[ which(abs(z[i]-t)<=1.5) ])),
						V.rmd24m= sapply(seq_along(z), function(i)	median(V[ which(abs(z[i]-t)<=2) ]))
						)	
			}, by=c('STAT','RSKF')]
	YXp		<- merge(YXp, tmp, by=c('STAT','t','RSKF'))		
	YXp		<- melt(YXp, id.vars=c('t.Patient','Patient','t','RSKF','stageC','t.AgeC','ntPatient','ntPatientn','V'), measure.vars=c('wnt','V.rm12m','V.rmd12m','V.rmd18m','V.rmd24m'))	
	ggplot(YXp, aes(x=t, y=value, colour=t.AgeC, group=RSKF)) + geom_line() + geom_point(size=1.2) +
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(variable~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmeanstAgeCByTime_1yrollingmean',suffix,sep=''),infile),'.pdf',sep='')	
	ggsave(file=file, w=10,h=10)	
	#
	
	YXr		<- merge(tmp, dcast.data.table(YXp, t.Patient+Patient+t+ntPatientn+RSKF+stageC+t.AgeC+ALPHA~STAT, value.var='V'), by=c('t.Patient','Patient','t','ntPatientn','RSKF','stageC','t.AgeC','ALPHA'))
	
	wm4		<- gamlss( score.Y~t.stAgeC-1, sigma.formula=~t.stAgeC-1, data=YXr, family=GA() )
	wm5		<- gamlss( score.Y~ns(t, df=1)+t.AgeC+stageC-1, sigma.formula=~t.AgeC+stageC-1, data=YXr, family=GA(mu.link='identity') )
	wm6		<- gamlss( score.Y~t.AgeC+t:t.AgeC+stageC-1, sigma.formula=~t.AgeC+stageC-1, data=YXr, family=GA() )
	wm7		<- gamlss( score.Y~t:t.stAgeC, sigma.formula=~t.stAgeC, data=YXr, family=GA(mu.link='identity') )	
	wm8		<- gamlss( score.Y~t.stAgeC+t:t.stAgeC-1, sigma.formula=~t.stAgeC, data=YXr, family=GA(mu.link='identity') )
	wm9st	<- rltvtp.gamlss.start(YXr, mu.formula='score.Y~t.stAgeC+t:t.stAgeC-1', sigma.formula='~t.stAgeC', gamlss.limit=c(0,1e-15,1e-10,1e-8), verbose=0)
	set(YXr, YXr[, which(score.Y<1e-15)], 'score.Y', 1e-15)
	wm9		<- gamlss( score.Y~t.stAgeC+t:t.stAgeC-1, sigma.formula=~t.stAgeC, data=YXr, family=WEI(), mu.start=wm9st$mu.start, sigma.start=wm9st$sigma.start )
	
	ywm4	<- predict(wm4, data=YXr, type='response', what='mu')
	ywm5	<- predict(wm5, data=YXr, type='response', what='mu')
	ywm6	<- predict(wm6, data=YXr, type='response', what='mu')
	ywm7	<- predict(wm7, data=YXr, type='response', what='mu')
	ywm8	<- predict(wm8, data=YXr, type='response', what='mu')
	ywm9	<- predict(wm9, data=YXr, type='response', what='mu')*log(2)^(1/predict(wm9, data=YXr, type='response', what='sigma'))	
	YXr[, score.Y.wm4:=ywm4]
	YXr[, score.Y.wm5:=ywm5]
	YXr[, score.Y.wm6:=ywm6]
	YXr[, score.Y.wm7:=ywm7]
	YXr[, score.Y.wm8:=ywm8]
	YXr[, score.Y.wm9:=ywm9]
	
	tmp		<- dcast.data.table(YXp, t.Patient+Patient+t~variable, value.var='value')
	tmp		<- merge(YXr, tmp, by=c('Patient','t.Patient','t'))	
	tmp		<- melt(tmp, measure.vars=c('score.Y.wm4', 'score.Y.wm5', 'score.Y.wm6', 'score.Y.wm7', 'score.Y.wm8','score.Y.wm9'))
	
	ggplot(subset(tmp, variable!='score.Y.wm9'), aes(x=t, colour=t.AgeC, group=t.stAgeC)) + geom_line(aes(y=value)) + 
			geom_point(aes(y=V.rm12m), size=1.2) +			
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			scale_y_continuous(limits=c(0,40)) +
			facet_grid(variable~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmeanstAgeCByTimeModels_1yrollingmean',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=10,h=10)	
	
	ggplot(subset(tmp, variable=='score.Y.wm9'), aes(x=t, colour=t.AgeC, group=t.stAgeC)) + geom_line(aes(y=value)) + 
			#geom_point(aes(y=V.rm12m), size=1.2, shape=1) +
			geom_point(aes(y=V.rmd18m), size=1.4) +
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			scale_y_continuous(limits=c(0,40)) +
			facet_grid(variable~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmeanstAgeCByTimeModels_1yrollingmedian',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=10,h=7)		
}
######################################################################################
rltvtp.explore.time.150828<- function(YXc, indir=NA, infile=NA, suffix='')
{
	#	it s very difficult to look at anything before 2003
	#	before 2004.5, the GA ID link model leads to negative values.. 		
	YXr		<- subset(YXc, AnyPos_T1>=2005, select=c('t.Patient','t','Patient', 'score.Y','ntPatient','stageC','t.AgeC', 't.stAgeC'))
	YXr[, ntPatientn:= as.numeric(as.character(ntPatient))]
	YXp		<- copy(YXr)
	setnames(YXp, 't.stAgeC','RSKF')
	tmp		<- YXp[, list(nt=length(score.Y)), by=c('RSKF','t')]
	setkey(tmp, RSKF, t)
	tmp		<- tmp[, list(t=t, nt=nt, wnt= sapply(seq_along(t), function(i)	sum(nt[ which(abs(t[i]-t)<=1) ])), cnt=cumsum(nt)	), by=c('RSKF')]
	YXp		<- merge(YXp, tmp, by=c('RSKF','t'))		
	YXp		<- melt( YXp, measure.vars=c('score.Y'), variable.name='STAT', value.name='V' )
	setkey(YXp, STAT, RSKF, t)
	#tmp		<- YXp[, list(V=mean(V), nt=nt[1], cnt=cnt[1], wnt=wnt[1]), by=c('t','RSKF','STAT')]
	#tmp		<- tmp[, list(t=t, V.rm=rollapply(V, width=ceiling(max(cnt)/100), FUN=mean, align="center", partial=TRUE)), by=c('STAT','RSKF')]	
	tmp		<- YXp[, {
				z	<- unique(t)
				list(t=z, V.rm= sapply(seq_along(z), function(i)	mean(V[ which(abs(z[i]-t)<=1) ]))	)	
			}, by=c('STAT','RSKF')]
	YXp		<- merge(YXp, tmp, by=c('STAT','t','RSKF'))		
	tmp		<- melt(subset(YXp, STAT=='score.Y'), id.vars=c('t.Patient','Patient','t','RSKF','stageC','t.AgeC','ntPatient','ntPatientn'), measure.vars=c('wnt'), variable.name='STAT', value.name='V')
	set(YXp, NULL, c('nt','wnt','cnt'), NULL)
	tmp[, V.rm:=V]	
	YXp		<- rbind(YXp, tmp, use.names=TRUE)
	tmp[, ALPHA:=0.2]
	set(tmp, tmp[, which(V>100)],'ALPHA',1)
	YXp		<- merge(YXp, subset(tmp, select=c('t.Patient','Patient','t','ALPHA')), by=c('t.Patient','Patient','t'))
	ggplot(YXp, aes(x=t, y=V.rm, colour=t.AgeC, group=RSKF, alpha=ALPHA)) + geom_line() + geom_point(size=1.2) +
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			facet_grid(STAT~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmeanstAgeCByTime_1yrollingmean',suffix,sep=''),infile),'.pdf',sep='')	
	ggsave(file=file, w=10,h=10)	
	#
	tmp		<- dcast.data.table(YXp, t.Patient+Patient+t+ntPatientn+RSKF+stageC+t.AgeC+ALPHA~STAT, value.var='V.rm')
	set(tmp, NULL, 'wnt', NULL)
	setnames(tmp, 'score.Y', 'score.Y.rm')	
	YXr		<- merge(tmp, dcast.data.table(YXp, t.Patient+Patient+t+ntPatientn+RSKF+stageC+t.AgeC+ALPHA~STAT, value.var='V'), by=c('t.Patient','Patient','t','ntPatientn','RSKF','stageC','t.AgeC','ALPHA'))
	wm4		<- gamlss( score.Y~RSKF-1, sigma.formula=~RSKF-1, data=YXr, family=GA() )
	wm5		<- gamlss( score.Y~ns(t, df=1)+t.AgeC+stageC-1, sigma.formula=~t.AgeC+stageC-1, data=YXr, family=GA(mu.link='identity') )
	wm6		<- gamlss( score.Y~t.AgeC+t:t.AgeC+stageC-1, sigma.formula=~t.AgeC+stageC-1, data=YXr, family=GA() )
	wm7		<- gamlss( score.Y~t:RSKF, sigma.formula=~RSKF, data=YXr, family=GA(mu.link='identity') )	
	wm8		<- gamlss( score.Y~RSKF+t:RSKF-1, sigma.formula=~RSKF, data=YXr, family=GA(mu.link='identity') )	
	ywm4	<- predict(wm4, data=YXr, type='response', what='mu')
	ywm5	<- predict(wm5, data=YXr, type='response', what='mu')
	ywm6	<- predict(wm6, data=YXr, type='response', what='mu')
	ywm7	<- predict(wm7, data=YXr, type='response', what='mu')
	ywm8	<- predict(wm8, data=YXr, type='response', what='mu')
	YXr[, score.Y.wm4:=ywm4]
	YXr[, score.Y.wm5:=ywm5]
	YXr[, score.Y.wm6:=ywm6]
	YXr[, score.Y.wm7:=ywm7]
	YXr[, score.Y.wm8:=ywm8]
	YXp		<- melt(YXr, measure.vars=c('score.Y.wm4', 'score.Y.wm5', 'score.Y.wm6', 'score.Y.wm7', 'score.Y.wm8'))
	ggplot(YXp, aes(x=t, colour=t.AgeC, group=RSKF)) + geom_line(aes(y=value)) + geom_point(aes(y=score.Y.rm), size=1.2) +
			scale_x_continuous(breaks=seq(1995,2020, 5), minor_breaks=seq(1995,2020,0.5)) +
			scale_y_continuous(limits=c(0,40)) +
			facet_grid(variable~stageC, scales='free') + labs(y='relative phylogenetic transmission probability\n1yr rolling mean') +
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey90", size=0.4), panel.grid.major=element_line(colour="grey90", size=0.4))
	file	<- paste(indir,'/',gsub('\\.R',paste('_scoreYmeanstAgeCByTimeModels_1yrollingmean',suffix,'.pdf',sep=''),infile),sep='')	
	ggsave(file=file, w=10,h=10)	
	
}
######################################################################################
rltvtp.model.150729<- function(YXc, indir=NA, infile=NA)
{
	YXr		<- subset(YXc, select=c('t.Patient','t','Patient', 'score.Y','ntPatient','stageC','t.AgeC', 't.stAgeC'))
	wm1		<- gamlss( score.Y~t.AgeC-1, sigma.formula=~t.AgeC-1, data=YXr, family=GA() )
	wm2		<- gamlss( score.Y~stageC-1, sigma.formula=~stageC-1, data=YXr, family=GA() )
	wm3		<- gamlss( score.Y~t.stAgeC-1, sigma.formula=~1, data=YXr, family=GA() )
	wm4		<- gamlss( score.Y~t.stAgeC-1, sigma.formula=~t.stAgeC-1, data=YXr, family=GA() )
	
	#	get model parameters so that we can add densities to histograms and compare means
	tmp		<- as.data.table(expand.grid(t.AgeC=YXr[, levels(t.AgeC)], stageC=YXr[, levels(stageC)]))
	tmp[, t.stAgeC:= tmp[, factor(paste(as.character(stageC),'_',as.character(t.AgeC),sep=''), levels=YXr[, levels(t.stAgeC)])]]	
	mu		<- predict(wm4, data=YXr, newdata=tmp, type='response', what='mu')
	sigma	<- predict(wm4, data=YXr, newdata=tmp, type='response', what='sigma')
	tmp[, mu:=mu]
	tmp[, sigma:=sigma]
	tmp[, TYPE:='wm4']
	wpars	<- copy(tmp)
	tmp		<- as.data.table(expand.grid(t.AgeC=YXr[, levels(t.AgeC)], stageC=YXr[, levels(stageC)]))
	tmp[, t.stAgeC:= tmp[, factor(paste(as.character(stageC),'_',as.character(t.AgeC),sep=''), levels=YXr[, levels(t.stAgeC)])]]	
	mu		<- predict(wm3, data=YXr, newdata=tmp, type='response', what='mu')
	sigma	<- predict(wm3, data=YXr, newdata=tmp, type='response', what='sigma')
	tmp[, mu:=mu]
	tmp[, sigma:=sigma]
	tmp[, TYPE:='wm3']
	wpars	<- rbind(wpars,tmp)
	tmp		<- as.data.table(expand.grid(t.AgeC=YXr[, levels(t.AgeC)], stageC=YXr[, levels(stageC)]))
	tmp[, t.stAgeC:= tmp[, factor(paste(as.character(stageC),'_',as.character(t.AgeC),sep=''), levels=YXr[, levels(t.stAgeC)])]]	
	mu		<- predict(wm2, data=YXr, newdata=tmp, type='response', what='mu')
	sigma	<- predict(wm2, data=YXr, newdata=tmp, type='response', what='sigma')
	tmp[, mu:=mu]
	tmp[, sigma:=sigma]
	tmp[, TYPE:='wm2']
	wpars	<- rbind(wpars,tmp)
	tmp		<- as.data.table(expand.grid(t.AgeC=YXr[, levels(t.AgeC)], stageC=YXr[, levels(stageC)]))
	tmp[, t.stAgeC:= tmp[, factor(paste(as.character(stageC),'_',as.character(t.AgeC),sep=''), levels=YXr[, levels(t.stAgeC)])]]	
	mu		<- predict(wm1, data=YXr, newdata=tmp, type='response', what='mu')
	sigma	<- predict(wm1, data=YXr, newdata=tmp, type='response', what='sigma')
	tmp[, mu:=mu]
	tmp[, sigma:=sigma]
	tmp[, TYPE:='wm1']
	wpars	<- rbind(wpars,tmp)
	
	if(!is.na(indir))
	{
		#	show if means differ from mus
		tmp		<- YXr[, list(score.Y.me= mean(score.Y), TYPE='wm1'), by='t.AgeC']
		wparp	<- merge(wpars, tmp, by=c('TYPE','t.AgeC'))
		tmp		<- YXr[, list(score.Y.me= mean(score.Y), TYPE='wm2'), by='stageC']
		tmp		<- merge(wpars, tmp, by=c('TYPE','stageC'))
		wparp	<- rbind(wparp, tmp, use.names=TRUE)
		tmp		<- YXr[, list(score.Y.me= mean(score.Y), TYPE='wm3'), by='t.stAgeC']
		tmp		<- merge(wpars, tmp, by=c('TYPE','t.stAgeC'))
		wparp	<- rbind(wparp, tmp, use.names=TRUE)
		tmp		<- YXr[, list(score.Y.me= mean(score.Y), TYPE='wm4'), by='t.stAgeC']
		tmp		<- merge(wpars, tmp, by=c('TYPE','t.stAgeC'))
		wparp	<- rbind(wparp, tmp, use.names=TRUE)
		ggplot(wparp, aes(x=score.Y.me, y=mu)) + geom_abline(intercept=0, slope=1,colour='grey50') + 
				geom_point() + facet_wrap(~TYPE,ncol=2) + theme_bw()
		file	<- paste(indir,'/',gsub('\\.R','_scoreYmeansVsMus.pdf',infile),sep='')	
		ggsave(file=file, w=6,h=6)
		#	--> spot on
	}
	list(wm1=wm1, wm2=wm2, wm3=wm3, wm4=wm4, par=wpars, data=YXr)
}
######################################################################################
rltvtp.gamlss.start<- function(YXr, mu.formula, sigma.formula, gamlss.limit=c(0,1e-15,1e-10,1e-8), verbose=0)
{
	mu.start		<- NULL
	sigma.start		<- NULL		
	for(i in rev(seq_along(gamlss.limit)[-1]))
	{
		if(verbose)
			cat('\nat iteration', gamlss.limit[i],'\n')
		tmp.mo		<- gamlss( 	as.formula(mu.formula), sigma.formula=as.formula(sigma.formula), 
				data=subset(YXr, score.Y>gamlss.limit[i]), family=WEI(), mu.start=mu.start, sigma.start=sigma.start)		
		tmp			<- subset(YXr, score.Y>gamlss.limit[i-1])
		mu.start	<- predict(tmp.mo, data=subset(YXr, score.Y>gamlss.limit[i]), newdata=tmp, type='response', what='mu')
		sigma.start	<- predict(tmp.mo, data=subset(YXr, score.Y>gamlss.limit[i]), newdata=tmp, type='response', what='sigma')		
	}		
	list(mu=mu.start, sigma=sigma.start)	
}
######################################################################################
rltvtp.confounded.ntransmitters<- function(YXc, indir, infile)
{
	tmp		<- YXc[, list(score.Y.m=mean(score.Y), stageC=stageC[1], t.AgeC=t.AgeC[1]), by=c('ntPatient','t.stAgeC')]
	ggplot( tmp, aes(x=ntPatient, y=score.Y.m)) + geom_bar(stat='identity') + facet_grid(stageC~t.AgeC)
	file	<- paste(indir,'/',gsub('\\.R','_scoreYmeansByStageNtPatient_raw.pdf',infile),sep='')	
	ggsave(file=file, w=20,h=10)
	tmp		<- YXc[, list(score.Y.me=mean(score.Y)), by=c('ntPatient')]
	ggplot( tmp, aes(x=ntPatient, y=score.Y.me)) + geom_bar(stat='identity')
	file	<- paste(indir,'/',gsub('\\.R','_scoreYmeansByNtPatient_raw.pdf',infile),sep='')	
	ggsave(file=file, w=5,h=5)	
	ggplot(YXc, aes(x=score.Y.raw, y=score.Y, colour=as.numeric(as.character(ntPatient)))) + geom_point(alpha=0.5) + geom_abline(slope=1, intercept=0) + theme_bw()
	file	<- paste(indir,'/',gsub('\\.R','_scoreYvsScoreYRaw.pdf',infile),sep='')	
	ggsave(file=file, w=5,h=5)		
}
######################################################################################
rltvtp.exploredistribution<- function(YXs, indir, infile)
{
	YXr		<- subset(YXs, select=c('t.Patient','t','Patient', 'score.Y','ntPatient','stageC','t.AgeC', rfactor)) 	
	#	Exp has default link log	
	mExp	<- gamlss( score.Y~1, data=YXr, family=EXP() )
	mExpNt	<- gamlss( score.Y~ntPatient, data=YXr, family=EXP() )
	#	Gamma has default links log log
	mGA		<- gamlss( score.Y~1, data=YXr, family=GA() )
	mGANt	<- gamlss( score.Y~ntPatient, data=YXr, family=GA() )
	mGANtS	<- gamlss( score.Y~ntPatient, sigma.formula=~ntPatient, data=YXr, family=GA() )
	#
	setkey(YXr, ntPatient)
	mpars	<- data.table(ntPatient= YXr[, levels(ntPatient)])
	mpars	<- merge(unique(YXr), mpars, by='ntPatient')
	mpars[, score.mu:= predict(mExp, data=YXr, newdata=mpars, type='response', what='mu')]
	mpars[, TYPE:= 'Exp']
	#
	tmp		<- data.table(ntPatient= YXr[, levels(ntPatient)])
	tmp		<- merge(unique(YXr), tmp, by='ntPatient')
	mu		<- predict(mExpNt, data=YXr, newdata=tmp, type='response', what='mu')	
	tmp[, score.mu:= mu]
	tmp[, TYPE:= 'ExpNt']
	mpars	<- rbind(mpars, tmp, use.names=TRUE, fill=TRUE)
	#	
	tmp		<- data.table(ntPatient= YXr[, levels(ntPatient)])
	tmp		<- merge(unique(YXr), tmp, by='ntPatient')
	mu		<- predict(mGA, data=YXr, newdata=tmp, type='response', what='mu')
	sigma	<- predict(mGA, data=YXr, newdata=tmp, type='response', what='sigma')
	tmp[, score.mu:= mu]
	tmp[, score.sigma:= sigma]
	tmp[, TYPE:= 'GA']
	mpars	<- rbind(mpars, tmp, use.names=TRUE, fill=TRUE)
	#
	tmp		<- data.table(ntPatient= YXr[, levels(ntPatient)])
	tmp		<- merge(unique(YXr), tmp, by='ntPatient')
	mu		<- predict(mGANt, data=YXr, newdata=tmp, type='response', what='mu')
	sigma	<- predict(mGANt, data=YXr, newdata=tmp, type='response', what='sigma')
	tmp[, score.mu:= mu]
	tmp[, score.sigma:= sigma]
	tmp[, TYPE:= 'GANt']
	mpars	<- rbind(mpars, tmp, use.names=TRUE, fill=TRUE)
	#
	tmp		<- data.table(ntPatient= YXr[, levels(ntPatient)])
	tmp		<- merge(unique(YXr), tmp, by='ntPatient')
	mu		<- predict(mGANtS, data=YXr, newdata=tmp, type='response', what='mu')
	sigma	<- predict(mGANtS, data=YXr, newdata=tmp, type='response', what='sigma')
	tmp[, score.mu:= mu]
	tmp[, score.sigma:= sigma]
	tmp[, TYPE:= 'GANtS']
	mpars	<- rbind(mpars, tmp, use.names=TRUE, fill=TRUE)
	#
	YXp		<- mpars[, {
				if(grepl('Exp',TYPE))
					z<- rEXP(1e4, mu=score.mu)
				if(grepl('GA',TYPE))
					z<- rGA(1e4, mu=score.mu, sigma=score.sigma)
				list(score.Y=z)
			}, by=c('ntPatient','TYPE')]
	tmp		<- subset(YXr, select=c(ntPatient, score.Y))
	tmp[, TYPE:= 'OBS']
	YXp		<- rbind(tmp, YXp,use.names=TRUE)
	#ggplot(YXp, aes(x=score.Y, group=TYPE, colour=TYPE)) + stat_ecdf() + facet_wrap(~ntPatient, ncol=5)
	ggplot(YXp, aes(x=score.Y, group=TYPE, colour=TYPE)) +
			scale_x_continuous(limits=c(0,50)) +
			scale_colour_brewer(palette='Set1') +
			geom_density() + 
			facet_wrap(~ntPatient, ncol=5, scales='free_y') 
	file	<- paste(indir,'/',gsub('\\.R','_scoreYdensitiesByNtPatient.pdf',infile),sep='')	
	ggsave(file=file, w=15,h=15)
	ggplot(subset(YXp, !grepl('Nt',TYPE)), aes(x=score.Y, group=TYPE, colour=TYPE)) +
			scale_x_continuous(limits=c(0,50)) +
			scale_colour_brewer(palette='Set1') +
			geom_density() 
	file	<- paste(indir,'/',gsub('\\.R','_scoreYdensities.pdf',infile),sep='')	
	ggsave(file=file, w=7,h=5)
#	--> GA model is better than Exp
}