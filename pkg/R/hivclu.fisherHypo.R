######################################################################################
project.athena.Fisheretal.Hypo.evaluate<- function()
{
	require(data.table)
	require(ape)
	#stop()
	resume					<- 1 	
	indir					<- paste(DATA,"fisheretal_150312",sep='/')
	outdir					<- paste(DATA,"fisheretal_150312",sep='/')		
	indir					<- paste(DATA,"fisheretal_150319",sep='/')
	outdir					<- paste(DATA,"fisheretal_150319",sep='/')		
	
	
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
		file			<- paste(indir, '/', infile, '_', gsub('/',':',insignat), '_', "method.Hypo.Rdata", sep='')		
		readAttempt			<- try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))		
	}
	if(!resume)
	{
		files					<- list.files(indir)
		files					<- files[ grepl('.R$',files) & grepl('Hypo',files) & grepl('tp4',files)]
		
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
		print(runs.opt)
		
		#subset(runs.opt, grepl('C3V100',method.brl))
		#	load risk estimates
		tmp			<- lapply(seq_len(nrow(runs.opt)), function(i)
				{
					tmp		<- paste(indir, runs.opt[i,file], sep='/')
					cat(paste('\nprocess file=',runs.opt[i,file]))
					tmp		<- load(tmp)					
					if(!any(colnames(averted)=='t.period'))
						averted[, t.period:= 4]
					if(!any(colnames(averted)=='STAT'))
					{
						averted[, STAT:= 'Pjx.e0cp']
						setnames(averted, c('Pjx.e0cp.sum', 'Pjx.e0cp.sum.h'), c('H0', 'H1'))
					}
					averted	<- averted[, list(AV=mean(1-H1/H0)), by=c('BS','STAT')]
					averted[, method.risk:=runs.opt[i,method.risk]]
					averted[, method.dating:=runs.opt[i,method.dating]]
					averted[, method.nodectime:=runs.opt[i,method.nodectime]]
					averted[, method.brl:=runs.opt[i,method.brl ]]
					averted[, method.denom:=runs.opt[i,method.denom]]
					averted[, method.recentctime:=runs.opt[i,method.recentctime ]]
					averted
				})
		runs.av		<- do.call('rbind', tmp)
		file			<- paste(indir, '/', infile, '_', gsub('/',':',insignat), '_', "method.Hypo.Rdata", sep='')
		save(runs.av, file=file)
	}
	
	runs.av.info	<- copy(runs.av)
	#	get central estimate
	tmp				<- subset(runs.av.info, BS<1)[, list(EST='central', AV=median(AV)), by=c('method.risk','method.dating','method.nodectime','method.brl','method.denom','method.recentctime','STAT')]
	#	get quantile bootstrap estimates
	runs.av.info	<- subset(runs.av.info, BS>=1)[, list(EST=paste('Q',100*c(0.025,0.25,0.5,0.75,0.975),sep=''), AV=quantile(AV, probs=c(0.025,0.25,0.5,0.75,0.975))), by=c('method.risk','method.dating','method.nodectime','method.brl','method.denom','method.recentctime','STAT')]
	runs.av.info	<- rbind(tmp, runs.av.info)	
	set(runs.av.info, runs.av.info[, which(AV<1e-3)], 'AV', 0)
	set(runs.av.info, NULL, 'HYPO', runs.av.info[, gsub('m2Awmx.wtn.','',method.risk)])
	set(runs.av.info, NULL, 'HYPO', runs.av.info[, regmatches(method.risk, regexpr('Hypo.*',method.risk))])
	set(runs.av.info, NULL, 'HYPO', runs.av.info[, gsub('.tp4','',HYPO)])
	set(runs.av.info, NULL, 'method.risk', runs.av.info[, regmatches(method.risk, regexpr('m2[^H]*',method.risk))])
	set(runs.av.info, NULL, 'method.risk', runs.av.info[, substr(method.risk, 1, nchar(method.risk)-1)])
	#runs.av.info[,list(n=length(EST), method.risk=unique(method.risk)),by='HYPO']	
	runs.av.info	<- dcast.data.table(runs.av.info, method.brl+method.risk+STAT+HYPO~EST, value.var='AV')
	runs.av.info[, TEST_TYPE:= 'NONE']
	set(runs.av.info, runs.av.info[, which(grepl('TestC|PROUDC|IPrEXC', HYPO))], 'TEST_TYPE', 'ELISA' )
	set(runs.av.info, runs.av.info[, which(grepl('TestA', HYPO))], 'TEST_TYPE', 'RNA' )
	runs.av.info[, TEST_FREQ:= 'NONE']
	tmp				<- runs.av.info[, which(grepl('A[0-9]+', HYPO))]
	set(runs.av.info, tmp, 'TEST_FREQ', substring(runs.av.info[tmp, regmatches(HYPO, regexpr('A[0-9]+', HYPO))],2) )
	tmp				<- runs.av.info[, which(grepl('C[0-9]+', HYPO))]
	set(runs.av.info, tmp, 'TEST_FREQ', substring(runs.av.info[tmp, regmatches(HYPO, regexpr('C[0-9]+', HYPO))],2) )
	runs.av.info[, TEST_COV:= 'NONE']
	tmp				<- runs.av.info[, which(grepl('m[0-9]+pc', HYPO))]
	set(runs.av.info, tmp, 'TEST_COV', substring(runs.av.info[tmp, regmatches(HYPO, regexpr('m[0-9]+pc', HYPO))],2) )
	set(runs.av.info, tmp, 'TEST_COV', runs.av.info[tmp, substr(TEST_COV, 1, nchar(TEST_COV)-2)] )
	runs.av.info[, PREV:= 'NONE']
	set(runs.av.info, runs.av.info[, which(!grepl('Test', HYPO) & !grepl('Prest', HYPO) & !grepl('Trest', HYPO) & grepl('ART', HYPO))], 'PREV', 'TasP')	
	set(runs.av.info, runs.av.info[, which(grepl('Test', HYPO) & TEST_TYPE!='RNA')], 'PREV', 'test')
	set(runs.av.info, runs.av.info[, which(grepl('Test', HYPO) & TEST_TYPE=='RNA')], 'PREV', 'test (RNA)')
	tmp				<- runs.av.info[, which(grepl('Test', HYPO) & grepl('ART', HYPO))]
	set(runs.av.info, tmp, 'PREV', runs.av.info[tmp, paste(PREV, '-treat',sep='')]) 	
	set(runs.av.info, runs.av.info[, which(grepl('Prest', HYPO) )], 'PREV', 'test-PrEP (all)')
	set(runs.av.info, runs.av.info[, which(grepl('Trest', HYPO) )], 'PREV', 'test-PrEP (<30 yrs)')
	runs.av.info[, PREP_EFF:= 'NONE']
	set(runs.av.info, runs.av.info[, which(grepl('PROUD', HYPO))], 'PREP_EFF', 'PrEP reduction incidence\n86%')
	set(runs.av.info, runs.av.info[, which(grepl('IPrEX', HYPO))], 'PREP_EFF', 'PrEP reduction incidence\n44%')
	tmp				<- runs.av.info[, which(grepl('Prest', HYPO) & grepl('ART', HYPO))]	
	set(runs.av.info, tmp, 'PREV', runs.av.info[tmp, paste(PREV, '-treat',sep='')])	
	tmp				<- runs.av.info[, which(grepl('Trest', HYPO) & grepl('ART', HYPO))]	
	set(runs.av.info, tmp, 'PREV', runs.av.info[tmp, paste(PREV, '-treat',sep='')])		
	tmp				<- runs.av.info[, which(grepl('ARTat500', HYPO))]
	set(runs.av.info, tmp, 'PREV', runs.av.info[tmp, paste(PREV, '(CD4<500)')])
	tmp				<- runs.av.info[, which(grepl('ImmediateART', HYPO))]
	set(runs.av.info, tmp, 'PREV', runs.av.info[tmp, paste(PREV, '(Immediate)')])	
	tmp				<- c(	"test", "test (RNA)",                                               
							"test-treat (CD4<500)", "test-treat (Immediate)", 
							"test (RNA)-treat (CD4<500)", "test (RNA)-treat (Immediate)", 
							"test-PrEP (<30 yrs)", "test-PrEP (<30 yrs)-treat (CD4<500)", "test-PrEP (<30 yrs)-treat (Immediate)",
							"test-PrEP (all)", "test-PrEP (all)-treat (CD4<500)", "test-PrEP (all)-treat (Immediate)"	)
	set( runs.av.info, NULL, 'PREV', runs.av.info[, factor(PREV, levels=rev(tmp), labels=rev(tmp))] )	
	runs.av.info	<- subset(runs.av.info, !is.na(PREV) & STAT=='Pjx.e0cp')
	#
	#	plot for Sciene TM revision
	#
	tmp				<- subset( runs.av.info, grepl('TasP', PREV), select=which(!grepl('TEST_',names(runs.av.info))) )
	tmp				<- merge(tmp, as.data.table(expand.grid(PREV=unique(tmp$PREV), TEST_COV=c(18,seq(30,70,20)), TEST_FREQ=12, TEST_TYPE='ELISA')), by='PREV', allow.cartesian=TRUE)	
	runs.av.plot	<- subset( runs.av.info, !TEST_FREQ%in%c('NONE','6') & !TEST_COV%in%c('NONE', '40','60')  )
	runs.av.plot	<- rbind(runs.av.plot, tmp, use.names=TRUE)	
	set(runs.av.plot, NULL, 'LEGEND', runs.av.plot[, paste('annual testing coverage\nof probable transmitters\n',TEST_COV,'%',sep='')])
	set(runs.av.plot, runs.av.plot[, which(grepl('18',LEGEND))], 'LEGEND', 'no improvements to\nannual testing coverage')
	tmp				<- c("no improvements to\nannual testing coverage","annual testing coverage\nof probable transmitters\n30%","annual testing coverage\nof probable transmitters\n50%","annual testing coverage\nof probable transmitters\n70%")
	set(runs.av.plot, NULL, 'LEGEND', runs.av.plot[, factor(LEGEND, levels=tmp, labels=tmp)])
	#
	#	plot for main text ie 3pa1H1.48C2V100bInfT7
	#
	setkey(runs.av.plot, PREV)
	ggplot( subset(runs.av.plot, grepl('test',PREV) & !grepl('IPrEX',HYPO) & !grepl('test (RNA)-treat',PREV, fixed=1) & STAT=='Pjx.e0cp' & method.risk=='m2Awmx.wtn' & method.brl=='3pa1H1.48C2V100bInfT7'), aes(x=PREV, fill=factor( grepl('test-treat', PREV) + 2*grepl('PrEP', PREV) + grepl('all', PREV)) ) ) +
			geom_hline(yintercept=50, colour="grey70", size=2) +
			scale_fill_manual(values=c("#DFC27D","#F6E8C3","#C7EAE5","#80CDC1"), name='hypothetical interventions\n in time period 09/07-10/12', guide=FALSE) +
			geom_boxplot(aes(ymin=Q2.5*100, ymax=Q97.5*100, lower=Q25*100, middle=Q50*100, upper=Q75*100), stat="identity", fatten=0) +
			geom_errorbar(aes(ymin=Q50*100,ymax=Q50*100), color='black', width=0.9, size=0.7) +
			scale_y_continuous(expand=c(0,0), limits=c(0, 100), breaks=seq(0,100,10), minor_breaks=seq(0,100,5)) +
			labs(x='', y='\nHIV infections amongst MSM in the transmission cohort\nthat could have been averted in 08/07 - 10/12\n(%)') + 
			coord_flip() +			
			theme_bw() + theme(panel.margin=unit(1.25,"lines"), legend.position='bottom', axis.text.y=element_text(size=12), panel.grid.major.x=element_line(colour="grey70", size=0.6), panel.grid.minor.x=element_line(colour="grey70", size=0.6), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
			facet_grid(~LEGEND)
	file			<- '~/Dropbox (Infectious Disease)/2014_MSMtransmission_ATHENA1303/151006_Prevention.pdf'
	ggsave(file=file, w=16, h=6)	
	#
	#	plot comparison 44% vs 86%
	#
	ggplot( subset(runs.av.plot, grepl('PROUD|IPrEX',HYPO) & STAT=='Pjx.e0cp' & method.risk=='m2Awmx.wtn' & method.brl=='3pa1H1.48C2V100bInfT7'), aes(x=PREV, fill=factor( grepl('test-treat', PREV) + 2*grepl('PrEP', PREV) + grepl('all', PREV)) ) ) +
			geom_hline(yintercept=50, colour="grey70", size=2) +
			scale_fill_manual(values=c("#C7EAE5","#80CDC1"), name='hypothetical interventions\n in time period 09/07-10/12', guide=FALSE) +
			geom_boxplot(aes(ymin=Q2.5*100, ymax=Q97.5*100, lower=Q25*100, middle=Q50*100, upper=Q75*100), stat="identity", fatten=0) +
			geom_errorbar(aes(ymin=Q50*100,ymax=Q50*100), color='black', width=0.9, size=0.7) +
			scale_y_continuous(expand=c(0,0), limits=c(0, 100), breaks=seq(0,100,10), minor_breaks=seq(0,100,5)) +
			labs(x='', y='\nHIV infections amongst MSM in the transmission cohort\nthat could have been averted in 08/07 - 10/12\n(%)') + 
			coord_flip() +			
			theme_bw() + theme(panel.margin=unit(1.25,"lines"), legend.position='bottom', axis.text.y=element_text(size=12), panel.grid.major.x=element_line(colour="grey70", size=0.6), panel.grid.minor.x=element_line(colour="grey70", size=0.6), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
			facet_grid(PREP_EFF~LEGEND)
	file			<- '~/Dropbox (Infectious Disease)/2014_MSMtransmission_ATHENA1303/151006_Prevention4486.pdf'
	ggsave(file=file, w=16, h=8)	
	
	
	
	#
	#	plot before Sciene TM revision
	#
	tmp				<- subset( runs.av.info, grepl('TasP', PREV), select=which(!grepl('TEST_',names(runs.av.info))) )
	tmp				<- merge(tmp, as.data.table(expand.grid(PREV=unique(tmp$PREV), TEST_COV=seq(30,70,20), TEST_FREQ=12, TEST_TYPE='ELISA')), by='PREV', allow.cartesian=TRUE)	
	runs.av.plot	<- subset( runs.av.info, !TEST_FREQ%in%c('NONE','6') & !TEST_COV%in%c('NONE', '40','60')  )
	runs.av.plot	<- rbind(runs.av.plot, tmp, use.names=TRUE)	
	set(runs.av.plot, NULL, 'LEGEND', runs.av.plot[, paste('annual testing coverage\nof probable transmitters\n',TEST_COV,'%',sep='')])
	#
	#	plot for main text ie 3pa1H1.48C2V100bInfT7
	#
	setkey(runs.av.plot, PREV)
	ggplot( subset(runs.av.plot, STAT=='Pjx.e0cp' & method.risk=='m2Awmx.wtn' & method.brl=='3pa1H1.48C2V100bInfT7'), aes(x=PREV, fill=factor( grepl('PrEP', PREV) + grepl('86%', PREV)) ) ) +
			geom_hline(yintercept=50, colour="grey70", size=2) +
			scale_fill_manual(values=c("#FED9A6", "#FFFFCC", "#E5D8BD"), name='hypothetical interventions\n in time period 09/07-10/12', guide=FALSE) +
			geom_boxplot(aes(ymin=Q2.5*100, ymax=Q97.5*100, lower=Q25*100, middle=Q50*100, upper=Q75*100), stat="identity", fatten=0) +
			geom_errorbar(aes(ymin=Q50*100,ymax=Q50*100), color='black', width=0.9, size=0.7) +
			scale_y_continuous(expand=c(0,0), limits=c(0, 100), breaks=seq(0,100,10), minor_breaks=seq(0,100,5)) +
			labs(x='', y='\nHIV infections amongst MSM\nthat could have been averted in 09/07 - 10/12\n(%)') + 
			coord_flip() +			
			theme_bw() + theme(panel.margin=unit(1.25,"lines"), legend.position='bottom', axis.text.y=element_text(size=12), panel.grid.major.x=element_line(colour="grey70", size=0.6), panel.grid.minor.x=element_line(colour="grey70", size=0.6), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
			facet_grid(~LEGEND)
	file			<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2014/MSMtransmission_ATHENA1303/150327_Prevention.pdf'
	ggsave(file=file, w=16, h=6)	
	#
	#	plot for SOM: across censoring adjustments etc
	#
	setkey(runs.av.plot, PREV)
	tmp		<- data.table(STAT=c('Pjx.e0cp','Pjx.e0','Pjx'), stat.legend=c('with adjustment for\ncensoring and sampling biases','with adjustment for\nsampling biases','no adjustment for\nsampling and censoring biases'))
	set(tmp, NULL, 'stat.legend', tmp[, factor(stat.legend, levels=stat.legend, labels=stat.legend)])
	tmp		<- merge(runs.av.plot, tmp, by='STAT')	
	ggplot( subset(tmp, method.risk=='m2Awmx.wtn' & method.brl=='3pa1H1.48C2V100bInfT7'), aes(x=PREV, fill=factor( grepl('PrEP', PREV) + grepl('86%', PREV)) ) ) +
			geom_hline(yintercept=50, colour="grey70", size=2) +
			scale_fill_manual(values=c("#FED9A6", "#FFFFCC", "#E5D8BD"), name='hypothetical interventions\n in time period 09/07-10/12', guide=FALSE) +
			geom_boxplot(aes(ymin=Q2.5*100, ymax=Q97.5*100, lower=Q25*100, middle=Q50*100, upper=Q75*100), stat="identity", fatten=0) +
			geom_errorbar(aes(ymin=Q50*100,ymax=Q50*100), color='black', width=0.9, size=0.7) +
			scale_y_continuous(expand=c(0,0), limits=c(0, 100), breaks=seq(0,100,10), minor_breaks=seq(0,100,5)) +
			labs(x='', y='\nHIV infections amongst MSM\nthat could have been averted in 09/07 - 10/12\n(%)') + 
			coord_flip() +			
			theme_bw() + theme(panel.margin=unit(1.25,"lines"), legend.position='bottom', axis.text.y=element_text(size=12), panel.grid.major.x=element_line(colour="grey70", size=0.6), panel.grid.minor.x=element_line(colour="grey70", size=0.6), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
			facet_grid(stat.legend~LEGEND)
	file			<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2014/MSMtransmission_ATHENA1303/150327_PreventionBySamplingAdj.pdf'
	ggsave(file=file, w=16, h=14)
	#
	#	plot for SOM: across likelihood methods
	#
	setkey(runs.av.plot, PREV)
	tmp		<- data.table(	method.risk=c("m2Awmx.wtn","m2Awmx","m2Awmx.noscore","m2Awmx.nophyloscore"), 
							method.legend=c(	'with phylogenetic transmission\nprobability per interval',
												'with phylogenetic transmission\nprobability per probable transmitter',
												'every interval\nequally likely',
												'every probable transmitter\n equally likely'))
	set(tmp, NULL, 'method.legend', tmp[, factor(method.legend, levels=method.legend, labels=method.legend)])
	tmp		<- merge(runs.av.plot, tmp, by='method.risk')	
	ggplot( subset(tmp, STAT=='Pjx.e0cp' &  method.brl=='3pa1H1.48C2V100bInfT7'), aes(x=PREV, fill=factor( grepl('PrEP', PREV) + grepl('86%', PREV)) ) ) +
			geom_hline(yintercept=50, colour="grey70", size=2) +
			scale_fill_manual(values=c("#FED9A6", "#FFFFCC", "#E5D8BD"), name='hypothetical interventions\n in time period 09/07-10/12', guide=FALSE) +
			geom_boxplot(aes(ymin=Q2.5*100, ymax=Q97.5*100, lower=Q25*100, middle=Q50*100, upper=Q75*100), stat="identity", fatten=0) +
			geom_errorbar(aes(ymin=Q50*100,ymax=Q50*100), color='black', width=0.9, size=0.7) +
			scale_y_continuous(expand=c(0,0), limits=c(0, 100), breaks=seq(0,100,10), minor_breaks=seq(0,100,5)) +
			labs(x='', y='\nHIV infections amongst MSM\nthat could have been averted in 09/07 - 10/12\n(%)') + 
			coord_flip() +			
			theme_bw() + theme(panel.margin=unit(1.25,"lines"), legend.position='bottom', axis.text.y=element_text(size=12), panel.grid.major.x=element_line(colour="grey70", size=0.6), panel.grid.minor.x=element_line(colour="grey70", size=0.6), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
			facet_grid(method.legend~LEGEND)
	file			<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2014/MSMtransmission_ATHENA1303/150327_PreventionByPhyloLkl.pdf'
	ggsave(file=file, w=14, h=14)	
	#
	#	plot for SOM: across exclusion criteria
	#
	tmp				<- data.table(	method.brl=c(	"3pa1H1.48C2V100bInfT7", "3pa1H1.94C2V100bInfT7", "3pa1H1.09C2V100bInfT7",					
											"3pa1H1.48C3V100bInfT7", "3pa1H1.48C1V100bInfT7", 
											"3pa1H1.48C3V100bInfs0.85T7", "3pa1H1.48C2V100bInfs0.85T7", "3pa1H1.48C1V100bInfs0.85T7",
											"3pa1H1.48C3V100bInfs0.7T7", "3pa1H1.48C2V100bInfs0.7T7", "3pa1H1.48C1V100bInfs0.7T7"), 
						 			method.legend=c( 'central estimate of HIV infection times\ncentral phylogenetic exclusion criteria',
											'lower estimate of HIV infection times',
											'upper estimate of HIV infection times',
											'phylogenetic exclusion criteria\ncoalescent compatibility < 30%\nclade frequency < 80%',
											'phylogenetic exclusion criteria\ncoalescent compatibility < 10%\nclade frequency < 80%',
											'phylogenetic exclusion criteria\ncoalescent compatibility < 30%\nclade frequency < 85%',
											'phylogenetic exclusion criteria\ncoalescent compatibility < 20%\nclade frequency < 85%',
											'phylogenetic exclusion criteria\ncoalescent compatibility < 10%\nclade frequency < 85%',
											'phylogenetic exclusion criteria\ncoalescent compatibility < 30%\nclade frequency < 70%',
											'phylogenetic exclusion criteria\ncoalescent compatibility < 20%\nclade frequency < 70%',
											'phylogenetic exclusion criteria\ncoalescent compatibility < 10%\nclade frequency < 70%'))										 
	set(tmp, NULL, 'method.legend', tmp[, factor(method.legend, levels=method.legend, labels=method.legend)])
	runs.av.plot	<- merge(runs.av.plot, tmp, by='method.brl')
	setkey(runs.av.plot, PREV)
	ggplot( runs.av.plot, aes(x=PREV, fill=factor( grepl('PrEP', PREV) + grepl('86%', PREV)) ) ) +
			geom_hline(yintercept=50, colour="grey70", size=2) +
			scale_fill_manual(values=c("#FED9A6", "#FFFFCC", "#E5D8BD"), name='hypothetical interventions\n in time period 09/07-10/12', guide=FALSE) +
			geom_boxplot(aes(ymin=Q2.5*100, ymax=Q97.5*100, lower=Q25*100, middle=Q50*100, upper=Q75*100), stat="identity", fatten=0) +
			geom_errorbar(aes(ymin=Q50*100,ymax=Q50*100), color='black', width=0.9, size=0.7) +
			scale_y_continuous(expand=c(0,0), limits=c(0, 100), breaks=seq(0,100,10), minor_breaks=seq(0,100,5)) +
			labs(x='', y='\nHIV infections amongst MSM\nthat could have been averted in 09/07 - 10/12\n(%)') + 
			coord_flip() +			
			theme_bw() + theme(panel.margin=unit(1.25,"lines"), legend.position='bottom', axis.text.y=element_text(size=12), panel.grid.major.x=element_line(colour="grey70", size=0.6), panel.grid.minor.x=element_line(colour="grey70", size=0.6), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
			facet_grid(method.legend~LEGEND)
	file			<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2014/MSMtransmission_ATHENA1303/150327_PreventionByExclusionCriteria.pdf'
	ggsave(file=file, w=14, h=45)
	#
	#	plot for SOM: across exclusion criteria gen distance
	#
	tmp				<- data.table(	method.brl=c(	"3pa1H1.48C2V100b0.02T7", "3pa1H1.48C2V100b0.04T7"	), 
									method.legend=c( 'central estimate of HIV infection times\ncentral phylogenetic exclusion criteria\nand genetic distance <2%',
													 'central estimate of HIV infection times\ncentral phylogenetic exclusion criteria\nand genetic distance <4%'	))										 
	set(tmp, NULL, 'method.legend', tmp[, factor(method.legend, levels=method.legend, labels=method.legend)])
	runs.av.plot	<- merge(runs.av.plot, tmp, by='method.brl')
	tmp				<- subset(runs.av.plot, method.brl=='3pa1H1.48C2V100b0.02T7' & HYPO=='HypoPrestPROUDC12m70pc70pc+ARTat500')
	set(tmp, 1L, 'method.legend','central estimate of HIV infection times\ncentral phylogenetic exclusion criteria\nand genetic distance <4%')
	runs.av.plot	<- rbind(runs.av.plot,tmp)
	setkey(runs.av.plot, PREV)
	
	ggplot( runs.av.plot, aes(x=PREV, fill=factor( grepl('PrEP', PREV) + grepl('86%', PREV)) ) ) +
			geom_hline(yintercept=50, colour="grey70", size=2) +
			scale_fill_manual(values=c("#FED9A6", "#FFFFCC", "#E5D8BD"), name='hypothetical interventions\n in time period 09/07-10/12', guide=FALSE) +
			geom_boxplot(aes(ymin=Q2.5*100, ymax=Q97.5*100, lower=Q25*100, middle=Q50*100, upper=Q75*100), stat="identity", fatten=0) +
			geom_errorbar(aes(ymin=Q50*100,ymax=Q50*100), color='black', width=0.9, size=0.7) +
			scale_y_continuous(expand=c(0,0), limits=c(0, 100), breaks=seq(0,100,10), minor_breaks=seq(0,100,5)) +
			labs(x='', y='\nHIV infections amongst MSM\nthat could have been averted in 09/07 - 10/12\n(%)') + 
			coord_flip() +			
			theme_bw() + theme(panel.margin=unit(1.25,"lines"), legend.position='bottom', axis.text.y=element_text(size=12), panel.grid.major.x=element_line(colour="grey70", size=0.6), panel.grid.minor.x=element_line(colour="grey70", size=0.6), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
			facet_grid(method.legend~LEGEND)
	file			<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2014/MSMtransmission_ATHENA1303/150618_PreventionByGenDist.pdf'
	ggsave(file=file, w=14, h=10)
	
	#
	#
	subset(runs.av.plot, TEST_COV=='50')
	
	#	quick test plot
	ggplot(runs.av.info, aes(x=HYPO))	+ geom_boxplot(aes(ymin=Q2.5, ymax=Q97.5, lower=Q25, middle=Q50, upper=Q75), stat="identity") + geom_point(aes(y=central), colour='red') +
			coord_flip()
	#	keep for figure
	select			<- c(	'HypoARTat500', 'HypoImmediateART', 'HypoTestC06m100pc', 'HypoTestA06m100pc', 'HypoTestC18m100pcARTat500', 'HypoTestC18m100pcImmediateART', 'HypoPrestC18m50pc', 'HypoPrestC18m50pcImmediateART', 'HypoPrestC18m60pcImmediateART', 'HypoPrestC18m60pcARTat500', 'HypoPrestC18m70pcARTat500')
	ggplot(subset(runs.av.info, HYPO%in%select), aes(x=HYPO))	+ geom_boxplot(aes(ymin=Q2.5, ymax=Q97.5, lower=Q25, middle=Q50, upper=Q75), stat="identity") + geom_point(aes(y=central), colour='red') +
			coord_flip()
	
	
	if(0)
	{
		tmp				<- data.table(	HYPO	= rev(select), 
				legend	= rev(c(	'ART at CD4<500', 'immediate ART', 'testing for HIV every 6 mo by 100%', 'testing for acute HIV every 6 mo by 100%',													
								'testing for HIV every 18 mo by 50% + ART at 500', 'testing for HIV every 18 mo by 50% + immediate ART',
								'testing for HIV every 18 mo + oral PrEP by 50%', 
								'testing for HIV every 18 mo + oral PrEP by 60% + ART at 500', 'testing for HIV every 18 mo + oral PrEP by 70% + ART at 500',
								'testing for HIV every 18 mo + oral PrEP by 50% + immediate ART', 'testing for HIV every 18 mo + oral PrEP by 60% + immediate ART'
						)),
				levels	= rev(factor(c( 0, 0, 0, 0,
										1, 2, 
										0, 
										3, 3, 
										4, 4))))
		set(tmp, NULL, 'legend', tmp[, factor(legend, levels=tmp$legend, labels=tmp$legend)])						
		runs.av.plot	<- merge(runs.av.info, tmp, by='HYPO')
		setkey(runs.av.plot, legend)
		ggplot( runs.av.plot, aes(x=legend, fill=levels) ) +
				geom_hline(yintercept=50, colour="grey70", size=2) +
				scale_fill_manual(values=c("white","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C"), name='hypothetical interventions\n in time period 09/07-10/12', guide=FALSE) +
				geom_boxplot(aes(ymin=Q2.5*100, ymax=Q97.5*100, lower=Q25*100, middle=Q50*100, upper=Q75*100), stat="identity", fatten=4) +
				scale_y_continuous(expand=c(0,0), limits=c(0, 100), breaks=seq(0,100,10), minor_breaks=seq(0,100,5)) +
				labs(x='', y='HIV infections among MSM\nthat could have been averted 09/07 - 10/12\n(%)') + 
				coord_flip() +			
				theme_bw() + theme(legend.position='bottom', legend.title=element_text(size=12), legend.text=element_text(size=12), axis.text.y=element_text(size=12), panel.grid.major.x=element_line(colour="grey70", size=0.6), panel.grid.minor.x=element_line(colour="grey70", size=0.6), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) 		
	}
	if(0)	#CROI poster: no acute testing for simplicity
	{
		select			<- c(	'HypoARTat500', 'HypoImmediateART', 'HypoTestC12m59pc', 'HypoTestC06m65pc', 'HypoTestC12m59pcARTat500', 'HypoTestC12m59pcImmediateART', 
								'HypoPrestC12m32pc50pc', 'HypoPrestC12m32pc50pcARTat500', 'HypoPrestC12m32pc50pcImmediateART', 
								'HypoPrestC12m46pc60pcARTat500', 'HypoPrestC12m46pc60pcImmediateART')
		tmp				<- data.table(	HYPO	= rev(select), 
				legend	= rev(c(	'ART at CD4<500', 'immediate ART', 'testing for HIV every 12 mo by 70%', 'testing for HIV every 6 mo by 70%',													
								'testing for HIV every 12 mo by 70% + ART at 500', 'testing for HIV every 12 mo by 70% + immediate ART',
								'testing for HIV every 12 mo + oral PrEP by 50%',
								'testing for HIV every 12 mo + oral PrEP by 50% + ART at 500', 'testing for HIV every 12 mo + oral PrEP by 50% + immediate ART',
								'testing for HIV every 12 mo + oral PrEP by 60% + ART at 500', 'testing for HIV every 12 mo + oral PrEP by 60% + immediate ART'
						)),
				levels	= rev(factor(c( 0, 0, 0, 0,
										0, 0, 
										0, 
										0, 0, 
										0, 0))))
		set(tmp, NULL, 'legend', tmp[, factor(legend, levels=tmp$legend, labels=tmp$legend)])						
		runs.av.plot	<- merge(runs.av.info, tmp, by='HYPO')
		setkey(runs.av.plot, legend)
		ggplot( runs.av.plot, aes(x=legend, fill=levels) ) +
				geom_hline(yintercept=50, colour="grey70", size=2) +
				scale_fill_manual(values=c("black"), name='hypothetical interventions\n in time period 09/07-10/12', guide=FALSE) +
				geom_boxplot(aes(ymin=Q2.5*100, ymax=Q97.5*100, lower=Q25*100, middle=Q50*100, upper=Q75*100), stat="identity", fatten=0) +
				geom_errorbar(aes(ymin=Q50*100,ymax=Q50*100), color="#FF7F00", width=0.9, size=0.7) +
				scale_y_continuous(expand=c(0,0), limits=c(0, 100), breaks=seq(0,100,10), minor_breaks=seq(0,100,5)) +
				labs(x='', y='\nHIV infections among MSM\nthat could have been averted 09/07 - 10/12\n(%)') + 
				coord_flip() +			
				theme_bw() + theme(legend.position='bottom', axis.text.x=element_text(size=14), axis.text.y=element_text(size=4), panel.grid.major.x=element_line(colour="grey70", size=0.6), panel.grid.minor.x=element_line(colour="grey70", size=0.6), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) 
		file			<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.HypoAverted.pdf", sep='')
		ggsave(file=file, w=6, h=6)
	}
	if(1)	#CROI poster: no acute testing for simplicity
	{
		select			<- c(	'HypoARTat500', 'HypoImmediateART', 'HypoTestC12m59pc', 'HypoTestA06m65pc', 'HypoTestC12m59pcARTat500', 'HypoTestC12m59pcImmediateART', 
				'HypoPrestC12m32pc50pc', 'HypoPrestC12m32pc50pcARTat500', 'HypoPrestC12m32pc50pcImmediateART', 
				'HypoPrestC12m46pc60pcARTat500', 'HypoPrestC12m46pc60pcImmediateART')
		tmp				<- data.table(	HYPO	= rev(select), 
				legend	= rev(c(	'ART at CD4<500', 'immediate ART', 'testing for HIV every 12 mo by 70%', 'testing for acute HIV every 6 mo by 70%',													
								'testing for HIV every 12 mo by 70% + ART at 500', 'testing for HIV every 12 mo by 70% + immediate ART',
								'testing for HIV every 12 mo + oral PrEP by 50%',
								'testing for HIV every 12 mo + oral PrEP by 50% + ART at 500', 'testing for HIV every 12 mo + oral PrEP by 50% + immediate ART',
								'testing for HIV every 12 mo + oral PrEP by 60% + ART at 500', 'testing for HIV every 12 mo + oral PrEP by 60% + immediate ART'
						)),
				levels	= rev(factor(c( 0, 0, 0, 0,
										0, 0, 
										0, 
										0, 0, 
										0, 0))))
		set(tmp, NULL, 'legend', tmp[, factor(legend, levels=tmp$legend, labels=tmp$legend)])						
		runs.av.plot	<- merge(runs.av.info, tmp, by='HYPO')
		setkey(runs.av.plot, legend)
		ggplot( runs.av.plot, aes(x=legend, fill=levels) ) +
				geom_hline(yintercept=50, colour="grey70", size=2) +
				scale_fill_manual(values=c("black"), name='hypothetical interventions\n in time period 09/07-10/12', guide=FALSE) +
				geom_boxplot(aes(ymin=Q2.5*100, ymax=Q97.5*100, lower=Q25*100, middle=Q50*100, upper=Q75*100), stat="identity", fatten=0) +
				geom_errorbar(aes(ymin=Q50*100,ymax=Q50*100), color="#FF7F00", width=0.9, size=0.7) +
				scale_y_continuous(expand=c(0,0), limits=c(0, 100), breaks=seq(0,100,10), minor_breaks=seq(0,100,5)) +
				labs(x='', y='\nHIV infections among MSM\nthat could have been averted 09/07 - 10/12\n(%)') + 
				coord_flip() +			
				theme_bw() + theme(legend.position='bottom', axis.text.x=element_text(size=14), axis.text.y=element_text(size=4), panel.grid.major.x=element_line(colour="grey70", size=0.6), panel.grid.minor.x=element_line(colour="grey70", size=0.6), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) 
		file			<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.HypoAverted.pdf", sep='')
		ggsave(file=file, w=6, h=6)
	}
	
}
######################################################################################
project.athena.Fisheretal.Hypo.run.median<- function(YXe, method.risk, predict.t2inf=NULL, t2inf.args=NULL, df.all=NULL, method.realloc='ImmediateART',  t.firstsuppressed=0.3, use.YXf= 1, bs.n=1e3, t.period=0.125, save.file=NA, resume=FALSE)
{
	#	get prop in ART stages after first suppression	
	#	method.realloc<- 'RPrEP100+TestA1y'
	stopifnot(grepl('Test|ImmediateART|ARTat500|RPrEP|RPrEP+ImmediateART',method.realloc))	
	options(warn=0)	
	if(resume & !is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(!resume || is.na(save.file) || inherits(readAttempt, "try-error"))
	{	
		#	censoring adjustment for nt.table
		YXf						<- copy(YXe$YXf)
		tmp						<- copy(YXe$X.tables$cens.table)
		setkey(tmp, stat, t.period, risk, factor)
		ct						<- unique(tmp)
		tmp2					<- copy(YXe$X.tables$cens.table.bs)
		setkey(tmp2, stat, t.period, risk, factor)
		ctb						<- unique(tmp2)
		tmp						<- project.athena.Fisheretal.censoring.model(ct, ctb, plot.file=NA )
		ct						<- copy(tmp$ctn)	
		setkey(ct, t.period, risk, factor)
		tmp						<- ct[, seq_len( length(n.adj)/length(unique(factor)) )]
		ct[, bs:=rep(tmp, nrow(ct)/length(tmp))]		#subset this for BS run					
		
		YX						<- copy(YXe$YX)	
		YX.h					<- copy(YXe$YX)	
		if(grepl('m2Awmx',method.risk))
		{
			set( YX, NULL, 'stage', YX[, CD4a.tperiod] )
			set( YX.h, NULL, 'stage', YX.h[, CD4a.tperiod] )			
		}
		if(grepl('m2Bwmx',method.risk))
		{
			set( YX, NULL, 'stage', YX[, CD4b.tperiod] )
			set( YX.h, NULL, 'stage', YX.h[, CD4b.tperiod] )			
		}
		if(grepl('m2Cwmx',method.risk))
		{
			set( YX, NULL, 'stage', YX[, CD4c.tperiod] )
			set( YX.h, NULL, 'stage', YX.h[, CD4c.tperiod] )			
		}
		X.tables				<- copy(YXe$X.tables)	
		nt.table				<- X.tables$nt.table
		nt.table.h				<- copy(X.tables$nt.table)
		cat(paste('\nusing method',method.realloc))
		df.trinfo				<<- NULL 
		df.uinfo				<<- NULL
		if(grepl('Test',method.realloc))
		{			
			tmp					<- project.athena.Fisheretal.Hypo.ReallocUToDiag.getYXetc( YX.h, nt.table.h, method.risk, predict.t2inf, t2inf.args, df.all, YXf=YXf, th.starttime=2008.5, th.endtime=2011, t.period=t.period, method.realloc=method.realloc, method.sample= 'stage=prop, y=median, t=start')			
			YX.h				<- copy(tmp$YX.h)
			nt.table.h			<- copy(tmp$nt.table.h)
			tmp					<- tmp$df.uinfo		#not NULL -- need this for reallocate.handler.cens
			df.uinfo			<<- copy(tmp)
		}
		if(grepl('ART',method.realloc))
		{
			tmp					<- project.athena.Fisheretal.Hypo.ReallocDiagToART.getYXetc(YX.h, nt.table.h, method.risk, t.firstsuppressed=0, t.delta=t.period, method.realloc=method.realloc, method.sample='pair, stage=prop, y=median')
			YX.h				<- copy(tmp$YX.h)
			nt.table.h			<- copy(tmp$nt.table.h)				
		}
		if(grepl('PrEP',method.realloc))
		{
			p.reachable			<- as.numeric(substring( regmatches(method.realloc,regexpr('RPrEP[0-9]+', method.realloc)), 6))/100
			cat(paste('\nsetting p.reachable=', p.reachable))
			tmp					<- project.athena.Fisheretal.Hypo.ReallocUToNone.getYXetc(YX.h, nt.table.h, method.risk, predict.t2inf, t2inf.args, df.all, p.reachable=p.reachable, th.starttime=2008.5, t.period=t.period, method.sample= 'stage=prop, y=median')
			YX.h				<- copy(tmp$YX.h)
			nt.table.h			<- copy(tmp$nt.table.h)
			tmp					<- tmp$df.trinfo		#not NULL -- need this for reallocate.handler.cens
			df.trinfo			<<- copy(tmp)
		}
		#
		#	prepare nt.table for YX and YX.hypothetical
		#	
		risk.df					<- data.table(risk='stage',factor=YXe$YX[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.suA.Y',4,sep='.'))
		set(risk.df, NULL, c('risk.ref','factor.ref','coef.ref'), 'None')
		nt.table				<- project.athena.Fisheretal.Wallinga.prep.nttable(nt.table, YX=YX)
		nt.table.h				<- project.athena.Fisheretal.Wallinga.prep.nttable(nt.table.h, YX=YX.h)
		stopifnot(!nrow(subset(nt.table.h, X.msm.e0<X.seq)))
		nt.table				<- project.athena.Fisheretal.Wallinga.censoring(ct, nt.table)
		#	use same censoring adjustment for nt.table.h ( we only have cens.table for YX, so that s all we can do )
		reallocate.handler.cens	<- project.athena.Fisheretal.Hypo.ReallocHandler.cens(method.realloc)
		nt.table.h				<- project.athena.Fisheretal.Wallinga.censoring(ct, nt.table.h, reallocate.handler.cens=reallocate.handler.cens)
		#	get censoring and sampling adjustments YX.hypothetical 
		adj.h		<- nt.table.h[, list(YX=sum(YX), X.clu=sum(X.clu), X.seq=sum(X.seq), X.msm.e0=sum(X.msm.e0), X.msm.e0cp=sum(X.msm.e0cp), nRec=length(unique(Patient))), by=c('risk','factor')]
		tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')
		adj.h[, PYs:= adj.h[[tmp]]]
		adj.h[, PTx:= adj.h[, YX] / adj.h[[tmp]]]
		adj.h[, Sx.e0:= PYs/X.msm.e0]					
		adj.h[, Sx.e0cp:= PYs/X.msm.e0cp]
		set(adj.h, adj.h[, which(!is.finite(Sx.e0))], 'Sx.e0', 1)
		set(adj.h, adj.h[, which(!is.finite(Sx.e0cp))], 'Sx.e0cp', 1)
		set(adj.h, adj.h[, which(!is.finite(PTx))], 'PTx', 1)
		risk.df.h	<- merge(risk.df, subset(adj.h, select=c(risk, factor, PYs, PTx, X.msm.e0, X.msm.e0cp)), by=c('risk','factor'))
		nt.table.h	<- merge(nt.table.h, subset(adj.h, select=c(risk, factor, Sx.e0, Sx.e0cp)), by=c('risk','factor'))
		#	get censoring and sampling adjustments YX 
		adj			<- nt.table[, list(YX=sum(YX), X.clu=sum(X.clu), X.seq=sum(X.seq), X.msm.e0=sum(X.msm.e0), X.msm.e0cp=sum(X.msm.e0cp), nRec=length(unique(Patient))), by=c('risk','factor')]
		tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')
		adj[, PYs:= adj[[tmp]]]
		adj[, PTx:= adj[, YX] / adj[[tmp]]]
		adj[, Sx.e0:= PYs/X.msm.e0]					
		adj[, Sx.e0cp:= PYs/X.msm.e0cp]
		set(adj, adj[, which(!is.finite(Sx.e0))], 'Sx.e0', 1)
		set(adj, adj[, which(!is.finite(Sx.e0cp))], 'Sx.e0cp', 1)
		set(adj, adj[, which(!is.finite(PTx))], 'PTx', 1)
		risk.df		<- merge(risk.df, subset(adj, select=c(risk, factor, PYs, PTx, X.msm.e0, X.msm.e0cp)), by=c('risk','factor'))
		nt.table	<- merge(nt.table, subset(adj, select=c(risk, factor, Sx.e0, Sx.e0cp)), by=c('risk','factor'))	
		#
		set(risk.df, NULL, 'risk', risk.df[, as.character(risk)])
		set(risk.df, NULL, 'factor', risk.df[, as.character(factor)])	
		set(nt.table, NULL, 'risk', nt.table[, as.character(risk)])
		set(nt.table, NULL, 'factor', nt.table[, as.character(factor)])
		set(nt.table, NULL, 'Patient', nt.table[, as.character(Patient)])
		set(nt.table, NULL, 'YX', NULL)
		set(risk.df.h, NULL, 'risk', risk.df.h[, as.character(risk)])
		set(risk.df.h, NULL, 'factor', risk.df.h[, as.character(factor)])	
		set(nt.table.h, NULL, 'risk', nt.table.h[, as.character(risk)])
		set(nt.table.h, NULL, 'factor', nt.table.h[, as.character(factor)])
		set(nt.table.h, NULL, 'Patient', nt.table.h[, as.character(Patient)])
		set(nt.table.h, NULL, 'YX', NULL)	
		#				
		missing		<- project.athena.Fisheretal.Wallinga.prep.expmissing(nt.table, risk.df, YX, YXf, use.YXf=use.YXf, method.missingy='y=median')
		missing.h	<- project.athena.Fisheretal.Wallinga.prep.expmissing(nt.table.h, risk.df.h, YX.h, YXf, use.YXf=use.YXf, method.missingy='y=median')
		#if(!is.na(method.reallocate))
		#	stopifnot( !nrow(subset(missing.h, grepl(method.reallocate,factor) & (YX.n>0 | YXm.sum.e0>0 | YXm.sum.e0cp>0))) )	
		#	calculate proportion of recipients averted
		averted		<- missing[, 	list(	Pjx.e0cp.sum= sum((yYX.sum+YXm.sum.e0cp)*YX.w)), by=c('risk','Patient')]
		tmp			<- missing.h[, 	list(	Pjx.e0cp.sum.h= sum((yYX.sum+YXm.sum.e0cp)*YX.w)), by=c('risk','Patient')]
		averted		<- merge(averted, tmp, by=c('risk','Patient'), all.x=1)
		set(averted, averted[, which(is.na(Pjx.e0cp.sum.h))], 'Pjx.e0cp.sum.h', 0)
		tmp			<- averted[, which(Pjx.e0cp.sum<Pjx.e0cp.sum.h)]
		tmp			<- averted[, which(round(Pjx.e0cp.sum,d=2)<round(Pjx.e0cp.sum.h,d=2))]
		cat(paste('\nFound Patients with Pjx.e0cp.sum<Pjx.e0cp.sum.h, n=',length(tmp)))
		set(averted, tmp, 'Pjx.e0cp.sum.h', averted[tmp, Pjx.e0cp.sum])
		averted[, BS:=0]
		#	averted[, BS:=bs.i/bs.n/10]
		#	averted[, median(1-Pjx.e0cp.sum.h/Pjx.e0cp.sum)]
		#	averted[, mean(1-Pjx.e0cp.sum.h/Pjx.e0cp.sum)]
		
		
		#	averted[, list(averted=mean(1-Pjx.e0cp.sum.h/Pjx.e0cp.sum)), by='BS'][, c(summary(averted), quantile(averted, p=c(0.025,0.5,0.975)))]
		#	
		#	need nt.table of YX for censoring		
		cat(paste('\nbootstrap Hypo.run, bs.n=',bs.n))
		averted.bs	<- lapply(seq_len(bs.n), function(bs.i)
				{
#print(bs.i)							
					if(bs.i%%100==0)	cat(paste('\nprocess bootstrap data sets bs.i=',bs.i))
					#	bootstrap over recently infected Patients
					X.tables.bs				<- copy(YXe$X.tables)
					YXf.bs.id				<- unique(subset(X.tables.bs$nt.table, select=Patient))
					YXf.bs.id				<- YXf.bs.id[	sample( seq_len(nrow(YXf.bs.id)), nrow(YXf.bs.id), replace=TRUE ), ]
					set(YXf.bs.id, NULL, 'Patient.bs',  YXf.bs.id[, paste(Patient,'_bs',seq_len(nrow(YXf.bs.id)),sep='')] )
					#	construct X.tables.bs and YX.bs								
					nt.table.bs				<- merge( X.tables.bs$nt.table, YXf.bs.id, by='Patient' )									
					YX.bs					<- as.data.table( merge( as.data.frame(YXe$YX), as.data.frame(YXf.bs.id), by='Patient' ) )
					setnames(nt.table.bs, c('Patient','Patient.bs'), c('Patient.b4bs','Patient'))
					setnames(YX.bs, c('Patient','Patient.bs'), c('Patient.b4bs','Patient'))							
					stopifnot( length(setdiff( YX.bs[, sort(unique(Patient))], subset( nt.table.bs, stat=='YX')[, sort(unique(Patient))] ))==0 )
					YX.h.bs					<- copy(YX.bs)
					nt.table.h.bs			<- copy(nt.table.bs)							
					#	bootstrap over reallocations	
					if(grepl('Test',method.realloc))
					{
						tmp					<- project.athena.Fisheretal.Hypo.ReallocUToDiag.getYXetc( YX.h.bs, nt.table.h.bs, method.risk, predict.t2inf, t2inf.args, df.all, YXf=YXf, th.starttime=2008.5, th.endtime=2011, t.period=t.period, method.realloc=method.realloc, method.sample= 'stage=sample, y=median, t=sample', verbose=FALSE)			
						YX.h.bs				<- copy(tmp$YX.h)
						nt.table.h.bs		<- copy(tmp$nt.table.h)
						tmp					<- copy(tmp$df.uinfo)
#print(tmp)								
						set(df.uinfo, NULL, 'UPT', tmp[, UPT])		#not NULL -- need this for reallocate.handler.cens	
					}					
					if(grepl('ART',method.realloc))
					{
						tmp					<- project.athena.Fisheretal.Hypo.ReallocDiagToART.getYXetc(YX.h.bs, nt.table.h.bs, method.risk, t.firstsuppressed=0, t.delta=t.period, method.realloc=method.realloc, method.sample='pair, stage=sample, y=sample', verbose=FALSE)
						YX.h.bs				<- copy(tmp$YX.h)
						nt.table.h.bs		<- copy(tmp$nt.table.h)									
					}
					if(grepl('PrEP',method.realloc))
					{
						p.reachable			<- as.numeric(substring( regmatches(method.realloc,regexpr('RPrEP[0-9]+', method.realloc)), 6))/100
						tmp					<- project.athena.Fisheretal.Hypo.ReallocUToNone.getYXetc(YX.h.bs, nt.table.h.bs, method.risk, predict.t2inf, t2inf.args, df.all, p.reachable=p.reachable, th.starttime=2008.5, t.period=t.period, method.sample= 'stage=sample, y=median', verbose=FALSE)
						YX.h.bs				<- copy(tmp$YX.h)
						nt.table.h.bs		<- copy(tmp$nt.table.h)
						tmp					<- copy(tmp$df.trinfo)	
						set(df.trinfo, NULL, 'Patient.nztr.h', 	tmp[,Patient.nztr.h])		# the number of bs recipients with a prob transmitter may vary
					}
					stopifnot( length(setdiff( YX.h.bs[, unique(Patient)], YX.bs[, unique(Patient)] ))==0  )
					stopifnot( length(setdiff( nt.table.h.bs[, unique(Patient)], nt.table.bs[, unique(Patient)] ))==0 )
					#	prepare nt.table for bs YX and YX.hypothetical
					nt.table.bs				<- project.athena.Fisheretal.Wallinga.prep.nttable(nt.table.bs, YX=YX.bs)
					nt.table.h.bs			<- project.athena.Fisheretal.Wallinga.prep.nttable(nt.table.h.bs, YX=YX.h.bs)							
					#	bs censoring probability				
					ct.bs					<- subset( ct, bs==sample(ct[, unique(bs)], 1) )
					#	bs censoring adjustment for nt.table
#print('OK1')			
					nt.table.bs				<- project.athena.Fisheretal.Wallinga.censoring(ct.bs, nt.table.bs, nt.table.4c=nt.table)
					#	use same censoring adjustment for nt.table.h ( we only have cens.table for YX, so that s all we can do )
					reallocate.handler.cens	<- project.athena.Fisheretal.Hypo.ReallocHandler.cens(method.realloc)
#print('OK2')							
					nt.table.h.bs			<- project.athena.Fisheretal.Wallinga.censoring(ct.bs, nt.table.h.bs, reallocate.handler.cens=reallocate.handler.cens)
#print('OK3')								
					stopifnot(length(setdiff( subset(nt.table.h.bs, YX>0)[, sort(unique(Patient))], subset(nt.table.bs, YX>0)[, unique(Patient)]))==0)
					#	bootstrap sample missing probabilities
					missing.bs				<- project.athena.Fisheretal.Hypo.prepmissingbs(YX.bs, nt.table.bs, method.risk, YXf=YXf, method.missingy='y=median')
					missing.h.bs			<- project.athena.Fisheretal.Hypo.prepmissingbs(YX.h.bs, nt.table.h.bs, method.risk, YXf=YXf, method.missingy='y=median')							
					stopifnot(length(setdiff( missing.h.bs[, sort(unique(Patient))], missing.bs[, sort(unique(Patient))]))==0)							
					#	calculate proportion of recipients averted	
					averted.bs			<- missing.bs[, 	list(	Pjx.e0cp.sum= sum((yYX.sum+yYXm.sum.e0cp)*YX.w)), by=c('risk','Patient')]
					tmp					<- missing.h.bs[, 	list(	Pjx.e0cp.sum.h= sum((yYX.sum+yYXm.sum.e0cp)*YX.w)), by=c('risk','Patient')]
					averted.bs			<- merge(averted.bs, tmp, by=c('risk','Patient'), all.x=TRUE)
					set(averted.bs, averted.bs[, which(is.na(Pjx.e0cp.sum.h))], 'Pjx.e0cp.sum.h', 0)
					tmp					<- averted.bs[, which(Pjx.e0cp.sum<Pjx.e0cp.sum.h)]
					#cat(paste('\nFound Patients with Pjx.e0cp.sum<Pjx.e0cp.sum.h, n=',length(tmp)))
					set(averted.bs, tmp, 'Pjx.e0cp.sum.h', averted.bs[tmp, Pjx.e0cp.sum])
					averted.bs[, BS:=bs.i]
					#averted.bs[, mean(1-Pjx.e0cp.sum.h/Pjx.e0cp.sum)]
					averted.bs
				})
		averted.bs	<- do.call('rbind',averted.bs)
		averted		<- rbind(averted, averted.bs)
		#tmp			<- averted[, list(averted=median(1-Pjx.e0cp.sum.h/Pjx.e0cp.sum)), by='BS']
		tmp			<- averted[, list(averted=mean(1-Pjx.e0cp.sum.h/Pjx.e0cp.sum)), by='BS']		
		print( subset(tmp, BS<1)[, quantile(averted, prob=c(0.025, 0.5, 0.975))] )
		print( subset( tmp, BS>1 )[, c(summary(averted), quantile(averted, prob=c(0.025, 0.5, 0.975)))] )
		if(!is.na(save.file))
		{
			cat(paste('\nsave to file', save.file))
			save(file=save.file, averted)
		}
	}
	averted
}	
######################################################################################
project.athena.Fisheretal.Hypo.run<- function(YXe, method.risk, predict.t2inf=NULL, t2inf.args=NULL, df.all=NULL, method.realloc='ImmediateART',  use.YXf= 1, bs.n=1e3, t.period=0.125, save.file=NA, resume=FALSE)
{
	#	get prop in ART stages after first suppression	
	#	method.realloc<- 'PrestIPrEXC12m50pc50pc+ImmediateART'
	#	method.realloc<- 'PrestPROUDC12m50pc50pc'
	#	method.realloc<- 'TrestPROUDC12m50pc50pc28y'
	#	method.realloc<- 'TestC12m18pc+ImmediateART'
	#	method.realloc<- 'ImmediateART'
	stopifnot(grepl('Test|ImmediateART|ARTat500|RPrEP|RPrEP|Prest|Trest',method.realloc))	
	options(warn=0)	
	if(resume & !is.na(save.file))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(save.file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",save.file))					
		options(show.error.messages = TRUE)		
	}
	if(!resume || is.na(save.file) || inherits(readAttempt, "try-error"))
	{	
		#	censoring adjustment for nt.table
		YXf						<- copy(YXe$YXf)
		tmp						<- copy(YXe$X.tables$cens.table)
		setkey(tmp, stat, t.period, risk, factor)
		ct						<- unique(tmp)
		tmp2					<- copy(YXe$X.tables$cens.table.bs)
		setkey(tmp2, stat, t.period, risk, factor)
		ctb						<- unique(tmp2)
		tmp						<- project.athena.Fisheretal.censoring.model(ct, ctb, plot.file=NA )
		ct						<- copy(tmp$ctn)	
		setkey(ct, t.period, risk, factor)
		tmp						<- ct[, seq_len( length(n.adj)/length(unique(factor)) )]
		ct[, bs:=rep(tmp, nrow(ct)/length(tmp))]		#subset this for BS run					
		averted	<- lapply(seq_len(bs.n), function(bs.i)
				{					
					YX						<- copy(YXe$YX)	
					YX.h					<- copy(YXe$YX)	
					if(grepl('m2Awmx',method.risk))
					{
						set( YX, NULL, 'stage', YX[, CD4a.tperiod] )
						set( YX.h, NULL, 'stage', YX.h[, CD4a.tperiod] )						
					}
					if(grepl('m2Bwmx',method.risk))
					{
						set( YX, NULL, 'stage', YX[, CD4b.tperiod] )
						set( YX.h, NULL, 'stage', YX.h[, CD4b.tperiod] )						
					}
					if(grepl('m2Cwmx',method.risk))
					{
						set( YX, NULL, 'stage', YX[, CD4c.tperiod] )
						set( YX.h, NULL, 'stage', YX.h[, CD4c.tperiod] )						
					}
					X.tables				<- copy(YXe$X.tables)	
					nt.table				<- X.tables$nt.table
					nt.table.h				<- copy(X.tables$nt.table)
					cat(paste('\nusing method',method.realloc))
					#
					#	we reallocate raw scores in this step. this means we also need to adjust the score.Y field that stores how the transmission probs were derived from the relative probs
					#	this is not necessary for YXf because there is no reallocation. so can just use score.Y
					#
					if(grepl('Test',method.realloc))
					{			
						tmp					<- project.athena.Fisheretal.Hypo.ReallocTest.getYXetc( YX.h, nt.table.h, method.risk, predict.t2inf, t2inf.args, df.all, YXf=YXf, th.starttime=2008.5, th.endtime=2011, t.period=t.period, method.realloc=method.realloc, method.sample= 'stage=sample, y=sample, t=sample')			
						YX.h				<- copy(tmp$YX.h)
						nt.table.h			<- copy(tmp$nt.table.h)						
					}
					if(grepl('Prest',method.realloc))
					{			
						tmp					<- project.athena.Fisheretal.Hypo.ReallocPrEPandTest.getYXetc( YX.h, nt.table.h, method.risk, predict.t2inf, t2inf.args, df.all, YXf=YXf, th.starttime=2008.5, th.endtime=2011, t.period=t.period, method.realloc=method.realloc, method.sample= 'stage=sample, y=sample, t=sample, eff=median')			
						YX.h				<- copy(tmp$YX.h)
						nt.table.h			<- copy(tmp$nt.table.h)						
					}
					if(grepl('Trest',method.realloc))
					{			
						tmp					<- project.athena.Fisheretal.Hypo.ReallocPrEPandTestForTarget.getYXetc( YX.h, nt.table.h, method.risk, predict.t2inf, t2inf.args, df.all, YXf=YXf, th.starttime=2008.5, th.endtime=2011, t.period=t.period, method.realloc=method.realloc, method.sample= 'stage=sample, y=sample, t=sample, eff=median')			
						YX.h				<- copy(tmp$YX.h)
						nt.table.h			<- copy(tmp$nt.table.h)						
					}					
					if(grepl('ART',method.realloc))
					{
						tmp					<- project.athena.Fisheretal.Hypo.ReallocDiagToART.getYXetc(YX.h, nt.table.h, method.risk, t.firstsuppressed=0, t.delta=t.period, method.realloc=method.realloc, method.sample='pair, stage=sample, y=sample')
						YX.h				<- copy(tmp$YX.h)
						nt.table.h			<- copy(tmp$nt.table.h)				
					}
					if(grepl('PrEP',method.realloc))
					{
						tmp					<- project.athena.Fisheretal.Hypo.ReallocPrEP.getYXetc(YX.h, nt.table.h, method.risk, predict.t2inf, t2inf.args, df.all, th.starttime=2008.5, t.period=t.period, method.realloc=method.realloc, method.sample= 'stage=sample, eff=median')
						YX.h				<- copy(tmp$YX.h)
						nt.table.h			<- copy(tmp$nt.table.h)						
					}
					#
					#	prepare nt.table for YX and YX.hypothetical
					#	
					risk.df					<- data.table(risk='stage',factor=YXe$YX[, levels(stage)], risk.ref='stage', factor.ref=paste('ART.suA.Y',4,sep='.'))
					set(risk.df, NULL, c('risk.ref','factor.ref','coef.ref'), 'None')
					nt.table				<- project.athena.Fisheretal.Wallinga.prep.nttable(nt.table, YX=YX)
					nt.table.h				<- project.athena.Fisheretal.Wallinga.prep.nttable(nt.table.h, YX=YX.h)
					stopifnot(!nrow(subset(nt.table.h, X.msm.e0<X.seq)))
					nt.table				<- project.athena.Fisheretal.Wallinga.censoring(ct, nt.table)
					nt.table.h				<- project.athena.Fisheretal.Wallinga.censoring(ct, nt.table.h, nt.table.4c=nt.table)
					#	get censoring and sampling adjustments YX.hypothetical 
					adj.h		<- nt.table.h[, list(YX=sum(YX), X.clu=sum(X.clu), X.seq=sum(X.seq), X.msm.e0=sum(X.msm.e0), X.msm.e0cp=sum(X.msm.e0cp), nRec=length(unique(Patient))), by=c('risk','factor')]
					tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')
					adj.h[, PYs:= adj.h[[tmp]]]
					adj.h[, PTx:= adj.h[, YX] / adj.h[[tmp]]]
					adj.h[, Sx.e0:= PYs/X.msm.e0]					
					adj.h[, Sx.e0cp:= PYs/X.msm.e0cp]
					set(adj.h, adj.h[, which(!is.finite(Sx.e0))], 'Sx.e0', 1)
					set(adj.h, adj.h[, which(!is.finite(Sx.e0cp))], 'Sx.e0cp', 1)
					set(adj.h, adj.h[, which(!is.finite(PTx))], 'PTx', 1)
					risk.df.h	<- merge(risk.df, subset(adj.h, select=c(risk, factor, PYs, PTx, X.msm.e0, X.msm.e0cp)), by=c('risk','factor'))
					nt.table.h	<- merge(nt.table.h, subset(adj.h, select=c(risk, factor, Sx.e0, Sx.e0cp)), by=c('risk','factor'))
					#	get censoring and sampling adjustments YX 
					adj			<- nt.table[, list(YX=sum(YX), X.clu=sum(X.clu), X.seq=sum(X.seq), X.msm.e0=sum(X.msm.e0), X.msm.e0cp=sum(X.msm.e0cp), nRec=length(unique(Patient))), by=c('risk','factor')]
					tmp			<- ifelse(grepl('clu',method.risk), 'X.clu', 'X.seq')
					adj[, PYs:= adj[[tmp]]]
					adj[, PTx:= adj[, YX] / adj[[tmp]]]
					adj[, Sx.e0:= PYs/X.msm.e0]					
					adj[, Sx.e0cp:= PYs/X.msm.e0cp]
					set(adj, adj[, which(!is.finite(Sx.e0))], 'Sx.e0', 1)
					set(adj, adj[, which(!is.finite(Sx.e0cp))], 'Sx.e0cp', 1)
					set(adj, adj[, which(!is.finite(PTx))], 'PTx', 1)
					risk.df		<- merge(risk.df, subset(adj, select=c(risk, factor, PYs, PTx, X.msm.e0, X.msm.e0cp)), by=c('risk','factor'))
					nt.table	<- merge(nt.table, subset(adj, select=c(risk, factor, Sx.e0, Sx.e0cp)), by=c('risk','factor'))	
					#
					set(risk.df, NULL, 'risk', risk.df[, as.character(risk)])
					set(risk.df, NULL, 'factor', risk.df[, as.character(factor)])	
					set(nt.table, NULL, 'risk', nt.table[, as.character(risk)])
					set(nt.table, NULL, 'factor', nt.table[, as.character(factor)])
					set(nt.table, NULL, 'Patient', nt.table[, as.character(Patient)])
					set(risk.df.h, NULL, 'risk', risk.df.h[, as.character(risk)])
					set(risk.df.h, NULL, 'factor', risk.df.h[, as.character(factor)])	
					set(nt.table.h, NULL, 'risk', nt.table.h[, as.character(risk)])
					set(nt.table.h, NULL, 'factor', nt.table.h[, as.character(factor)])
					set(nt.table.h, NULL, 'Patient', nt.table.h[, as.character(Patient)])
					#				
					missing		<- project.athena.Fisheretal.Wallinga.prep.expmissing(nt.table, risk.df, YX, YXf, use.YXf=use.YXf, method.missingy='y=median', verbose=0)
					missing.h	<- project.athena.Fisheretal.Wallinga.prep.expmissing(nt.table.h, risk.df.h, YX.h, YXf, use.YXf=use.YXf, method.missingy='y=median', verbose=0)
					#if(!is.na(method.reallocate))
					#	stopifnot( !nrow(subset(missing.h, grepl(method.reallocate,factor) & (YX.n>0 | YXm.sum.e0>0 | YXm.sum.e0cp>0))) )	
					#	calculate proportion of recipients averted
					averted		<- missing[, 	list(	Pjx= sum(yYX.sum*YX.w),
														Pjx.e0= sum((yYX.sum+YXm.sum.e0)*YX.w),
														Pjx.e0cp= sum((yYX.sum+YXm.sum.e0cp)*YX.w),
														TYPE='H0'), by=c('risk','Patient')]
					tmp			<- missing.h[, 	list(	Pjx= sum(yYX.sum*YX.w),
														Pjx.e0= sum((yYX.sum+YXm.sum.e0)*YX.w),									
														Pjx.e0cp= sum((yYX.sum+YXm.sum.e0cp)*YX.w),
														TYPE='H1'), by=c('risk','Patient')]
					averted		<- melt(averted, id.vars=c('risk','Patient','TYPE'), value.name='SUM', variable.name='STAT')
					tmp			<- melt(tmp, id.vars=c('risk','Patient','TYPE'), value.name='SUM', variable.name='STAT')
					averted		<- merge( dcast.data.table(averted, risk+Patient+STAT~TYPE, value.var='SUM'), dcast.data.table(tmp, risk+Patient+STAT~TYPE, value.var='SUM'), by=c('risk','Patient','STAT'), all.x=1)					
					set(averted, averted[, which(is.na(H1))], 'H1', 0)					
					tmp			<- averted[, which(round(H0,d=2)<round(H1,d=2))]
					cat(paste('\nFound Patients with H0.sum<H1.sum, n=',length(tmp)))
					set(averted, tmp, 'H0', averted[tmp, H1])
					averted[, BS:=bs.i/bs.n/10]
					#	averted[, median(1-Pjx.e0cp.sum.h/Pjx.e0cp.sum)]
					#	averted[, mean(1-Pjx.e0cp.sum.h/Pjx.e0cp.sum)]
				})
		averted		<- do.call('rbind',averted)
		#	averted[, list(averted=mean(1-Pjx.e0cp.sum.h/Pjx.e0cp.sum)), by='BS'][, c(summary(averted), quantile(averted, p=c(0.025,0.5,0.975)))]
		#	
		#	need nt.table of YX for censoring		
		nt.table	<- copy(YXe$X.tables$nt.table)	
		nt.table	<- project.athena.Fisheretal.Wallinga.prep.nttable(nt.table)		
		cat(paste('\nbootstrap Hypo.run, bs.n=',bs.n))
		averted.bs	<- lapply(seq_len(bs.n), function(bs.i)
				{
#print(bs.i)							
					if(bs.i%%100==0)	cat(paste('\nprocess bootstrap data sets bs.i=',bs.i))
					#	bootstrap over recently infected Patients
					X.tables.bs				<- copy(YXe$X.tables)
					YXf.bs.id				<- unique(subset(X.tables.bs$nt.table, select=Patient))
					YXf.bs.id				<- YXf.bs.id[	sample( seq_len(nrow(YXf.bs.id)), nrow(YXf.bs.id), replace=TRUE ), ]
					set(YXf.bs.id, NULL, 'Patient.bs',  YXf.bs.id[, paste(Patient,'_bs',seq_len(nrow(YXf.bs.id)),sep='')] )
					#	construct X.tables.bs and YX.bs								
					nt.table.bs				<- merge( X.tables.bs$nt.table, YXf.bs.id, by='Patient' )									
					YX.bs					<- as.data.table( merge( as.data.frame(YXe$YX), as.data.frame(YXf.bs.id), by='Patient' ) )
					setnames(nt.table.bs, c('Patient','Patient.bs'), c('Patient.b4bs','Patient'))
					setnames(YX.bs, c('Patient','Patient.bs'), c('Patient.b4bs','Patient'))							
					stopifnot( length(setdiff( YX.bs[, sort(unique(Patient))], subset( nt.table.bs, stat=='YX')[, sort(unique(Patient))] ))==0 )
					YX.h.bs					<- copy(YX.bs)
					nt.table.h.bs			<- copy(nt.table.bs)							
					#	bootstrap over reallocations	
					if(grepl('Test',method.realloc))
					{
						tmp					<- project.athena.Fisheretal.Hypo.ReallocTest.getYXetc( YX.h.bs, nt.table.h.bs, method.risk, predict.t2inf, t2inf.args, df.all, YXf=YXf, th.starttime=2008.5, th.endtime=2011, t.period=t.period, method.realloc=method.realloc, method.sample= 'stage=sample, y=sample, t=sample', verbose=TRUE)			
						YX.h.bs				<- copy(tmp$YX.h)
						nt.table.h.bs		<- copy(tmp$nt.table.h)
					}					
					if(grepl('Prest',method.realloc))
					{			
						#YX.h<- YX.h.bs; nt.table<- nt.table.h.bs
						tmp					<- project.athena.Fisheretal.Hypo.ReallocPrEPandTest.getYXetc( YX.h.bs, nt.table.h.bs, method.risk, predict.t2inf, t2inf.args, df.all, YXf=YXf, th.starttime=2008.5, th.endtime=2011, t.period=t.period, method.realloc=method.realloc, method.sample= 'stage=sample, y=sample, t=sample, eff=sample', verbose=FALSE)			
						YX.h.bs				<- copy(tmp$YX.h)
						nt.table.h.bs		<- copy(tmp$nt.table.h)						
					}
					if(grepl('Trest',method.realloc))
					{			
						tmp					<- project.athena.Fisheretal.Hypo.ReallocPrEPandTestForTarget.getYXetc( YX.h.bs, nt.table.h.bs, method.risk, predict.t2inf, t2inf.args, df.all, YXf=YXf, th.starttime=2008.5, th.endtime=2011, t.period=t.period, method.realloc=method.realloc, method.sample= 'stage=sample, y=sample, t=sample, eff=median', verbose=FALSE)			
						YX.h.bs				<- copy(tmp$YX.h)
						nt.table.h.bs		<- copy(tmp$nt.table.h)						
					}										
					if(grepl('ART',method.realloc))
					{
						tmp					<- project.athena.Fisheretal.Hypo.ReallocDiagToART.getYXetc(YX.h.bs, nt.table.h.bs, method.risk, t.firstsuppressed=0, t.delta=t.period,  method.realloc=method.realloc, method.sample='pair, stage=sample, y=sample', verbose=FALSE)
						YX.h.bs				<- copy(tmp$YX.h)
						nt.table.h.bs		<- copy(tmp$nt.table.h)									
					}
					if(grepl('PrEP',method.realloc))
					{						
						tmp					<- project.athena.Fisheretal.Hypo.ReallocPrEP.getYXetc(YX.h.bs, nt.table.h.bs, method.risk, predict.t2inf, t2inf.args, df.all, th.starttime=2008.5, t.period=t.period, method.realloc=method.realloc, method.sample= 'stage=sample, eff=sample', verbose=FALSE)
						YX.h.bs				<- copy(tmp$YX.h)
						nt.table.h.bs		<- copy(tmp$nt.table.h)
					}
					#z<- merge(nt.table.bs, nt.table.h.bs, by=c('Patient','risk','factor','stat'))
					stopifnot( length(setdiff( YX.h.bs[, unique(Patient)], YX.bs[, unique(Patient)] ))==0  )
					stopifnot( length(setdiff( nt.table.h.bs[, unique(Patient)], nt.table.bs[, unique(Patient)] ))==0 )
					#	prepare nt.table for bs YX and YX.hypothetical
					nt.table.bs				<- project.athena.Fisheretal.Wallinga.prep.nttable(nt.table.bs, YX=YX.bs)
					nt.table.h.bs			<- project.athena.Fisheretal.Wallinga.prep.nttable(nt.table.h.bs, YX=YX.h.bs)
					#z<- merge(nt.table.bs, nt.table.h.bs, by=c('Patient','risk','factor'))
					#	bs censoring probability				
					ct.bs					<- subset( ct, bs==sample(ct[, unique(bs)], 1) )
					#	bs censoring adjustment for nt.table
#print('OK1')			
					nt.table.bs				<- project.athena.Fisheretal.Wallinga.censoring(ct.bs, nt.table.bs, nt.table.4c=nt.table)
					nt.table.h.bs			<- project.athena.Fisheretal.Wallinga.censoring(ct.bs, nt.table.h.bs, nt.table.4c=nt.table)
#print('OK3')				
#z<- merge(nt.table.bs, nt.table.h.bs, by=c('Patient','risk','factor'))
					stopifnot(length(setdiff( subset(nt.table.h.bs, YX>0)[, sort(unique(Patient))], subset(nt.table.bs, YX>0)[, unique(Patient)]))==0)
					#	bootstrap sample missing probabilities
					set.seed(bs.i)
					missing.bs				<- project.athena.Fisheretal.Hypo.prepmissingbs(YX.bs, nt.table.bs, method.risk, YXf=YXf, method.missingy='y=median')
					set.seed(bs.i)
					missing.h.bs			<- project.athena.Fisheretal.Hypo.prepmissingbs(YX.h.bs, nt.table.h.bs, method.risk, YXf=YXf, method.missingy='y=median')							
					stopifnot(length(setdiff( missing.h.bs[, sort(unique(Patient))], missing.bs[, sort(unique(Patient))]))==0)							
					#	calculate proportion of recipients averted	
					averted.bs				<- missing.bs[, 	list(	Pjx= sum(yYX.sum*YX.w),
																		Pjx.e0= sum((yYX.sum+yYXm.sum.e0)*YX.w),
																		Pjx.e0cp= sum((yYX.sum+yYXm.sum.e0cp)*YX.w),
																		TYPE='H0'), by=c('risk','Patient')]
					tmp						<- missing.h.bs[, 	list(	Pjx= sum(yYX.sum*YX.w),
																	Pjx.e0= sum((yYX.sum+yYXm.sum.e0)*YX.w),									
																	Pjx.e0cp= sum((yYX.sum+yYXm.sum.e0cp)*YX.w),
																	TYPE='H1'), by=c('risk','Patient')]
					averted.bs				<- melt(averted.bs, id.vars=c('risk','Patient','TYPE'), value.name='SUM', variable.name='STAT')
					tmp						<- melt(tmp, id.vars=c('risk','Patient','TYPE'), value.name='SUM', variable.name='STAT')
					averted.bs				<- merge( dcast.data.table(averted.bs, risk+Patient+STAT~TYPE, value.var='SUM'), dcast.data.table(tmp, risk+Patient+STAT~TYPE, value.var='SUM'), by=c('risk','Patient','STAT'), all.x=1)					
					set(averted.bs, averted.bs[, which(is.na(H1))], 'H1', 0)					
					tmp						<- averted.bs[, which(round(H0,d=2)<round(H1,d=2))]
					#cat(paste('\nFound Patients with H0.sum<H1.sum, n=',length(tmp)))
					set(averted.bs, tmp, 'H0', averted.bs[tmp, H1])
					averted.bs[, BS:=bs.i]
					#averted.bs[, list(AV=mean(1-H1/H0)), by='STAT']
					averted.bs
				})
		averted.bs	<- do.call('rbind',averted.bs)
		averted		<- rbind(averted, averted.bs)
		tmp			<- averted[, list(AV=mean(1-H1/H0)), by=c('STAT','BS')]			
		print( subset(tmp, BS<1)[, list(P=c(0.025, 0.5, 0.975), AV=quantile(AV, prob=c(0.025, 0.5, 0.975))), by='STAT'] )
		print( subset( tmp, BS>1 )[, list(P=c(0.025, 0.5, 0.975), AV=quantile(AV, prob=c(0.025, 0.5, 0.975))), by='STAT'] )
		if(!is.na(save.file))
		{
			cat(paste('\nsave to file', save.file))
			save(file=save.file, averted)
		}
	}
	averted
}	
######################################################################################
project.athena.Fisheretal.Hypo.prepmissingbs<- function(YX.bs, nt.table.bs, method.risk, YX=NULL, YXf=NULL, method.missingy='y=median')
{	
	# 	boostrap sample Xseq and PTx.bs
	missing.bs			<- nt.table.bs[, list(YX=sum(YX), X.clu=sum(X.clu), X.seq=sum(X.seq), X.msm.e0=sum(X.msm.e0), X.msm.e0cp=sum(X.msm.e0cp), nRec=length(unique(Patient))), by=c('risk','factor')]
	tmp					<- ifelse(grepl('clu',method.risk),'X.clu','X.seq')
	set(missing.bs, NULL, 'PYs', rbinom( nrow(missing.bs), sum(missing.bs[[tmp]]), missing.bs[[tmp]]/sum(missing.bs[[tmp]]) ))
	tmp					<- unique(subset(nt.table.bs, select=c(risk,factor)))
	tmp					<- tmp[,	{
				z	<- table( YX.bs[, risk, with=FALSE])
				list(factor=rownames(z), YX.bs=as.numeric(unclass(z)))												
			}, by='risk']
	missing.bs			<- merge(missing.bs, tmp, by=c('risk','factor'), all.x=TRUE)	
	missing.bs[, PTx.bs:= YX.bs/PYs]
	# 	boostrap sample X.msm.e0 and X.msm.e0cp
	set(missing.bs, NULL, 'X.msm.e0.bs', missing.bs[, rbinom(length(X.msm.e0), sum(X.msm.e0), X.msm.e0/sum(X.msm.e0))])
	set(missing.bs, NULL, 'X.msm.e0cp.bs', missing.bs[, rbinom(length(X.msm.e0cp), sum(X.msm.e0cp), X.msm.e0cp/sum(X.msm.e0cp))])	
	missing.bs[, Sx.e0.bs:= PYs/X.msm.e0.bs]					
	missing.bs[, Sx.e0cp.bs:= PYs/X.msm.e0cp.bs]
	#	can be >1 due to bootstrapping
	set(missing.bs, missing.bs[, which(Sx.e0.bs>1 | !is.finite(Sx.e0.bs))], 'Sx.e0.bs', 1)
	set(missing.bs, missing.bs[, which(Sx.e0cp.bs>1 | !is.finite(Sx.e0cp.bs))], 'Sx.e0cp.bs', 1)
	#	zero denominator means missing YX.bs, PTx.bs and NaN Sx.e0.bs, Sx.e0cp.bs
	tmp					<- missing.bs[, which(X.msm.e0==0)]
	set(missing.bs, tmp, c('YX.bs'), 0)
	set(missing.bs, tmp, c('PTx.bs','Sx.e0.bs','Sx.e0cp.bs'), 1)	
	#	prepare risk.df and nt.table.bs as needed
	missing.bs		<- merge(nt.table.bs, subset(missing.bs, select=c(risk, factor, PTx.bs, Sx.e0.bs, Sx.e0cp.bs)), by=c('risk','factor'))
	#	reduce to bootstrap sampled recipients
	missing.bs		<- merge(unique(subset(YX.bs, select=Patient)), missing.bs, by='Patient')				
	#	compute the sum of observed Y's by risk factor for each recipient
	tmp				<- YX.bs[, list(yYX.sum= sum(score.Y), YX.bs=length(score.Y), YX.w=w.t[1]), by=c('stage','Patient')]
	setnames(tmp, 'stage','factor')
	set(tmp, NULL, 'factor', tmp[, as.character(factor)])	
	tmp[, risk:='stage']	
	missing.bs		<- merge(missing.bs, tmp, by=c('risk','factor','Patient'), all.x=TRUE)
	set(missing.bs, missing.bs[, which(is.na(YX.bs))], 'YX.w', 1.) 
	set(missing.bs, missing.bs[, which(is.na(YX.bs))], c('yYX.sum','YX.bs'), 0.)
	#	bootstrap sample number missing: to avoid rounding issues, sample total and then re-allocate to individual recipient MSM
	tmp			<- missing.bs[, {
				#	total missed through sampling
				YXm.r.e0	<- rznbinom( 1, sum(YX.bs), Sx.e0.bs[1] )
				YXm.r.e0cp	<- rznbinom( 1, sum(YX.bs), Sx.e0cp.bs[1] )
				#	distribute among all recipients as fraction
				YXm.r.e0	<- as.vector(rmultinom(1, YXm.r.e0, rep( 1/length(YX.bs), length(YX.bs))))
				YXm.r.e0cp	<- as.vector(rmultinom(1, YXm.r.e0cp, rep( 1/length(YX.bs), length(YX.bs))))
				list(Patient=Patient, YXm.r.e0=YXm.r.e0, YXm.r.e0cp=YXm.r.e0cp)
			}, by=c('risk','factor')]				
	missing.bs	<- merge(missing.bs, tmp, by=c('risk','factor','Patient'))	
	#	draw missing scores from all yijt in that stage	for number missing YXm.r.e0
	if(is.null(YXf))
		YX.used	<- YX
	if(!is.null(YXf))
		YX.used	<- YXf
	missing.bs[, factor2:= substr(factor, 1, nchar(factor)-1)]
	if(grepl('y=mean',method.missingy))
		tmp			<- missing.bs[, {		
					#z	<- YX.used[ which( grepl(factor2[1], YX.used[[risk]], fixed=TRUE) ), ]
					z	<- subset(YX.used, grepl(factor2[1], YX.used[[risk]], fixed=TRUE), c('Patient','t.Patient','score.Y'))
					setkey(z, Patient, t.Patient)
					z	<- unique(z)																																																			
					z	<- mean( z[['score.Y']] )																						
					list(Patient=Patient, yYXm.sum.e0=YXm.r.e0*z, yYXm.sum.e0cp=YXm.r.e0cp*z )
				}, by=c('risk','factor')]	
	if(grepl('y=median',method.missingy))
		tmp			<- missing.bs[, {		
					#z	<- YX.used[ which( grepl(factor2[1], YX.used[[risk]], fixed=TRUE) ), ]
					z	<- subset(YX.used, grepl(factor2[1], YX.used[[risk]], fixed=TRUE), c('Patient','t.Patient','score.Y'))
					setkey(z, Patient, t.Patient)
					z	<- unique(z)
					z	<- median( z[['score.Y']] )																						
					list(Patient=Patient, yYXm.sum.e0=YXm.r.e0*z, yYXm.sum.e0cp=YXm.r.e0cp*z )
				}, by=c('risk','factor')]
	if(grepl('y=sample.i',method.missingy))
		tmp			<- missing.bs[, {		
					yYXm.sum.e0		<- sum(sample(YX.used[ which( grepl(factor2[1], YX.used[[risk]], fixed=TRUE)), ][['score.Y']], YXm.r.e0, replace=TRUE))
					yYXm.sum.e0cp	<- sum(sample(YX.used[ which( grepl(factor2[1], YX.used[[risk]], fixed=TRUE)), ][['score.Y']], YXm.r.e0cp, replace=TRUE))
					list( yYXm.sum.e0=yYXm.sum.e0, yYXm.sum.e0cp=yYXm.sum.e0cp )
				}, by=c('risk','factor','Patient')]
	if(grepl('y=sample.p',method.missingy))
		tmp			<- missing.bs[, {		
					z	<- subset(YX.used, grepl(factor2[1], YX.used[[risk]], fixed=TRUE), c('Patient','t.Patient','score.Y'))
					setkey(z, Patient, t.Patient)
					z	<- my.sample(unique(z)[['score.Y']], length(Patient), replace=TRUE)																															
					list(Patient=Patient, yYXm.sum.e0=YXm.r.e0*z, yYXm.sum.e0cp=YXm.r.e0cp*z )
				}, by=c('risk','factor')]						
	missing.bs[, factor2:=NULL]				
	missing.bs	<- merge(missing.bs, tmp, by=c('Patient','risk','factor'), all.x=TRUE)	
	missing.bs
}
######################################################################################
project.athena.Fisheretal.Hypo.ReallocHandler.cens<- function(method.realloc)
{
	stopifnot(grepl('ImmediateART|ARTat500|RPrEP|Test|Prest',method.realloc))
	if('ImmediateART'==method.realloc)
	{		
		reallocate.handler.cens	<- function(ct.table)
		{
			set(ct.table, ct.table[, which(grepl('Dt|DA',factor))], 'n.adj.med', 0)
			set(ct.table, ct.table[, which(grepl('ART',factor))], 'n.adj.med', ct.table[which(grepl('ART',factor)),X.msm.e0])
			if( ct.table[ grepl('UA.', factor, fixed=TRUE), X.msm.e0>n.adj.med ] )
			{
				cat(paste('\nWARNING: Found UA X.msm.e0>n.adj.med (may happen occasionally)'))
				tmp		<- ct.table[, which(grepl('UA.', factor, fixed=TRUE))]
				set(ct.table, tmp, 'n.adj.med', ct.table[tmp, X.msm.e0])
			}	
			ct.table
		}			
	}
	if('ARTat500'==method.realloc)
	{		
		reallocate.handler.cens	<- function(ct.table)
		{
			set(ct.table, ct.table[, which(grepl('Dtl500|Dtl350',factor))], 'n.adj.med', 0)
			tmp		<- ct.table[, which(grepl('Dtg500|DA|Dt.NA',factor) & X.msm.e0>n.adj.med)]
			set(ct.table, tmp, 'n.adj.med', ct.table[tmp, X.msm.e0])		#occurs half the time because of bs sampling of X.msm.e0
			set(ct.table, ct.table[, which(grepl('ART',factor))], 'n.adj.med', ct.table[which(grepl('ART',factor)),X.msm.e0])
			if( ct.table[ grepl('UA.', factor, fixed=TRUE), X.msm.e0>n.adj.med ] )
			{
				cat(paste('\nWARNING: Found UA X.msm.e0>n.adj.med (may happen occasionally)'))
				tmp		<- ct.table[, which(grepl('UA.', factor, fixed=TRUE))]
				set(ct.table, tmp, 'n.adj.med', ct.table[tmp, X.msm.e0])
			}						
			ct.table
		}	
	}
	if(grepl('Test',method.realloc) & !grepl('RPrEP',method.realloc))	#same with or without ART
	{		
		reallocate.handler.cens	<- function(ct.table)
		{
			ct.table	<- merge(ct.table, df.uinfo, by=c('risk','factor'))
			set( ct.table, NULL, 'n.adj.med', ct.table[, round(n.adj.med*(1-UPT))] )
			tmp			<- ct.table[, which(grepl('ART|Dt|DA', factor))]
			set( ct.table, tmp, 'n.adj.med', ct.table[tmp,X.msm.e0] )			
			if( ct.table[ grepl('UA.', factor, fixed=TRUE), X.msm.e0>n.adj.med ] )
			{
				cat(paste('\nWARNING: Found UA X.msm.e0>n.adj.med (may happen occasionally)'))
				tmp		<- ct.table[, which(grepl('UA.', factor, fixed=TRUE))]
				set(ct.table, tmp, 'n.adj.med', ct.table[tmp, X.msm.e0])
			}						
			ct.table
		}	
	}
	if(grepl('RPrEP',method.realloc) & !grepl('Test',method.realloc))	#same with or without ART
	{
		reallocate.handler.cens	<- function(ct.table)
		{	
			ct.table	<- merge(ct.table, df.trinfo, by=c('risk','factor'))			
			#set(ct.table, NULL, 'n.adj.med', ct.table[, round(n.adj.med*pt*Patient.nztr.h/Patient.nztr)])
			set(ct.table, NULL, 'n.adj.med', ct.table[, round(n.adj.med*pt*Patient.n.h/Patient.n)])
			tmp			<- ct.table[, which(grepl('^ART|^D', factor))]
			set(ct.table, tmp, 'n.adj.med', ct.table[tmp, X.msm.e0] )			
			if( ct.table[ grepl('UA.', factor, fixed=TRUE), X.msm.e0>n.adj.med ] )
			{
				cat(paste('\nWARNING: Found UA X.msm.e0>n.adj.med (may happen occasionally)'))
				tmp		<- ct.table[, which(grepl('UA.', factor, fixed=TRUE))]
				set(ct.table, tmp, 'n.adj.med', ct.table[tmp, X.msm.e0])
			}			
			ct.table							
		}		
	}		
	if(grepl('RPrEP',method.realloc) & grepl('Test',method.realloc))	#same with or without ART
	{
		reallocate.handler.cens	<- function(ct.table)
		{	
			ct.table	<- merge(ct.table, df.uinfo, by=c('risk','factor'))
			set( ct.table, NULL, 'n.adj.med', ct.table[, round(n.adj.med*(1-UPT))] )			
			ct.table	<- merge(ct.table, df.trinfo, by=c('risk','factor'))			
			#set(ct.table, NULL, 'n.adj.med', ct.table[, round(n.adj.med*pt*Patient.nztr.h/Patient.nztr)])
			set(ct.table, NULL, 'n.adj.med', ct.table[, round(n.adj.med*pt*Patient.n.h/Patient.n)])
			tmp			<- ct.table[, which(grepl('^ART|^D', factor))]
			set(ct.table, tmp, 'n.adj.med', ct.table[tmp, X.msm.e0] )			
			if( ct.table[ grepl('UA.', factor, fixed=TRUE), X.msm.e0>n.adj.med ] )
			{
				cat(paste('\nWARNING: Found UA X.msm.e0>n.adj.med (may happen occasionally)'))
				tmp		<- ct.table[, which(grepl('UA.', factor, fixed=TRUE))]
				set(ct.table, tmp, 'n.adj.med', ct.table[tmp, X.msm.e0])
			}			
			ct.table							
		}		
	}	
	if(grepl('Prest',method.realloc))
	{
		reallocate.handler.cens	<- function(ct.table)
		{
			ct.table	<- merge(ct.table, df.uinfo, by=c('risk','factor'))
			#	adjust calculations for fewer recipients due to PrEP
			set( ct.table, NULL, 'n.adj.med', ct.table[, round(n.adj.med*Patient.n.h/Patient.n, d=0)] )
			set( ct.table, NULL, 'X.msm.e0', ct.table[, round(X.msm.e0*Patient.n.h/Patient.n, d=0)] )			
		}
	}
	reallocate.handler.cens
}
######################################################################################
project.athena.Fisheretal.Hypo.ReallocDiagToART.getYXetc<- function(YX, nt.table, method.risk, t.firstsuppressed=0.3, t.delta=0.125, method.realloc='ImmediateART', method.sample= 'stage=prop, y=mean', verbose=TRUE)
{
	stopifnot(grepl('tri|pair',method.sample))
	stopifnot(grepl('stage=sample|stage=prop',method.sample))
	stopifnot(grepl('y=sample|y=median|y=mean',method.sample))
	stopifnot(grepl('ImmediateART|ARTat500',method.realloc))
	
	r.artinfo	<- function(df.artinfo, n, method.sample)
	{
		setkey(df.artinfo, stage)
		tmp		<- unique(df.artinfo)
		setkey(tmp,pt)
		if(grepl('stage=prop',method.sample))
		{
			tmp		<- merge(tmp, data.table(stage=rep(tmp$stage, ceiling(tmp$pt*n))), by='stage')
			tmp		<- tmp[seq_len(n),]			
		}
		if(grepl('stage=sample',method.sample))
			tmp		<- merge(tmp, data.table(stage=cut(runif(n), breaks=c(-Inf,cumsum(tmp$pt)), labels=tmp$stage)), by='stage')	
		set(tmp, NULL, 'stage', tmp[, as.character(stage)])
		if(grepl('y=sample',method.sample))
			ans	<- tmp[ ,	list(score.Y= my.sample(df.artinfo[['score.Y.raw']][ df.artinfo[['stage']]==stage ], length(nt), replace=TRUE)), by='stage']
		if(grepl('y=median',method.sample))
			ans	<- tmp[ ,	list(score.Y= rep(median(df.artinfo[['score.Y.raw']][ df.artinfo[['stage']]==stage ]), length(nt))), by='stage']
		if(grepl('y=mean',method.sample))
			ans	<- tmp[ ,	list(score.Y= rep(mean(df.artinfo[['score.Y.raw']][ df.artinfo[['stage']]==stage ]), length(nt))), by='stage']		
		ans
	}	
	r.artinfo.s	<- function(YXf, tp, n, method.sample)
	{
		tmp	<- subset( YXf, grepl('ART.NotYetFirstSu',stage))$score.Y.raw
		ans	<- data.table(stage=rep(paste('ART.NotYetFirstSu.',tp,sep=''),n))
		if(grepl('y=sample',method.sample))
			ans[, score.Y:= my.sample( tmp, n, replace=TRUE)]
		if(grepl('y=median',method.sample))
			ans[, score.Y:= rep(median(tmp), n)]
		set(YX.h, ARTs.i, 'score.Y.raw', rep(median(subset( YXf, grepl('ART.NotYetFirstSu',stage))$score.Y.raw), length(ARTs.i)))	
		if(grepl('y=mean',method.sample))
			ans[, score.Y:= rep(mean(tmp), n)]		
		ans		
	}
	if(grepl('ImmediateART',method.realloc))
		method.reallocate		<- 'Dt|DA'		
	if(grepl('ARTat500',method.realloc))
		method.reallocate		<- 'Dtl500|Dtl350'		
	tp						<- as.numeric(substr(regmatches(method.risk,regexpr('tp[0-9]', method.risk)),3,3))	
	#
	#	prepare YX.hypothetical
	#		
	#ggplot(subset(YX, grepl('Dt|DA',stage))[, list(n=length(unique(as.character(stage)))), by=c('Patient','t.Patient')], aes(x=n)) + geom_histogram()
	tmp			<- subset(YX, grepl('ART',stage) & !grepl('ART.NotYetFirstSu',stage))[, table(stage)]
	df.artinfo	<- subset( data.table(stage=names(tmp), nt=tmp, pt=tmp/sum(tmp)), nt>0 )
#print(df.artinfo)	
	#	change to YX raw
	tmp			<- subset(YX, select=c(stage, Patient, t.Patient, score.Y.raw))	
	setkey(tmp, Patient, t.Patient, stage)
	tmp			<- unique(tmp)
	set( tmp, NULL, 'stage', tmp[, as.character(stage)] )		#use all ART intervals, not just those from tp4
	set( tmp, NULL, 'stage', tmp[, paste(substr(stage,1,nchar(stage)-1),tp,sep='')] )
	df.artinfo	<- merge(df.artinfo, tmp, by='stage')
	set(df.artinfo, NULL, 'stage', df.artinfo[, as.character(stage)])	
	#	reallocate ART after first viral suppression
	YX.h		<- copy(YX)	
	setkey(YX.h, Patient, t.Patient)
	ARTns.i		<- YX.h[, which( (t-t.AnyPos_T1+t.delta)>=t.firstsuppressed  &  grepl(method.reallocate,stage)  )]
	if(verbose)
		cat(paste('\nreallocate entries, n=', length(ARTns.i)))	
	if(grepl('pair',method.sample))
	{
		ARTns.pair	<- unique(subset(YX.h, t-t.AnyPos_T1+t.delta>=t.firstsuppressed  &  grepl(method.reallocate,stage), select=c(Patient,t.Patient)))	
		ARTns.pair	<- cbind(ARTns.pair, r.artinfo(df.artinfo, nrow(ARTns.pair), method.sample=method.sample))
		setnames(ARTns.pair, c('stage','score.Y'), c('stage2','score.Y2'))
		YX.h 		<- merge(YX.h, ARTns.pair, by=c('Patient', 't.Patient'), all.x=1)				
	}
	if(grepl('tri',method.sample))
	{
		tmp2		<- r.artinfo(df.artinfo, length(ARTns.i), method.sample=method.sample)
		set(YX.h, ARTns.i, 'stage', tmp2$stage )
		set(YX.h, ARTns.i, 'score.Y.raw', tmp2$score.Y )		
	}
	#	reallocate ART before first viral suppression
	ARTs.i			<- YX.h[, which(  (t-t.AnyPos_T1+t.delta)<t.firstsuppressed  &  grepl(method.reallocate,stage)  )]
	if(verbose)
		cat(paste('\nreallocate entries, n=', length(ARTs.i)))
	if(grepl('pair',method.sample))
	{
		ARTs.pair	<- unique(subset(YX.h, t-t.AnyPos_T1+t.delta<t.firstsuppressed  &  grepl(method.reallocate,stage), select=c(Patient,t.Patient)))
		ARTs.pair	<- cbind(ARTs.pair, r.artinfo.s(YX, tp, nrow(ARTs.pair), method.sample))
		setnames(ARTs.pair, c('stage','score.Y'), c('stage3','score.Y3'))
		YX.h 		<- merge(YX.h, ARTs.pair, by=c('Patient', 't.Patient'), all.x=1)
		# some of the patient pairs in ARTs.pair and ARTns.pair also have other intervals than diagnosed -- keep these
		tmp			<- YX.h[, which(   (!is.na(stage2) | !is.na(stage3))  &  !grepl(method.reallocate,stage) )]
		set(YX.h, tmp, c('stage2','stage3'), NA_character_)
		set(YX.h, tmp, c('score.Y2', 'score.Y3'), NA_real_)
		if(verbose)
			cat(paste('\nsetting new entries, n=', nrow(subset(YX.h, !is.na(stage2) | !is.na(stage3))) ))				
		tmp			<- YX.h[, which(!is.na(stage2))]
		set(YX.h, tmp, 'stage', YX.h[tmp, stage2])
		set(YX.h, tmp, 'score.Y.raw', YX.h[tmp, score.Y2])
		tmp			<- YX.h[, which(!is.na(stage3))]
		set(YX.h, tmp, 'stage', YX.h[tmp, stage3])
		set(YX.h, tmp, 'score.Y.raw', YX.h[tmp, score.Y3])		
		#subset(YX.h, !is.na(stage2) | !is.na(stage3))[, sort(table(as.character(stage))/sum(!is.na(stage)))]; unique(subset(df.artinfo, select=c(stage, pt)))
		#set(YX.h, NULL, c('stage2','stage3','score.Y2','score.Y3'), NULL)				
	}
	if(grepl('tri',method.sample))
	{
		tmp2		<- r.artinfo.s(YX, tp, length(ARTs.i), method.sample)
		set(YX.h, ARTs.i, 'stage', tmp2$stage )
		set(YX.h, ARTs.i, 'score.Y.raw', tmp2$score.Y )
	}	
	set(YX.h, NULL, 'score.Y', YX.h[, score.Y.raw])
	set(YX.h, NULL, 'score.Y.raw', NULL)	
	set(YX.h, NULL, 'stage', YX.h[, factor(as.character(stage))])
	#	
	if(grepl('wtn',method.risk))
	{
		if(verbose)
			cat(paste('\nsetting likelihood to likelihood of pair / number transmission intervals'))
		set(YX.h, NULL, 'score.Y', YX.h[, score.Y*w.tn])				
	}								
	if(grepl('noscore',method.risk))
	{
		if(verbose)
			cat(paste('\nsetting likelihood to 1'))
		set(YX.h, NULL, 'score.Y', YX.h[, 1])						
	}
	if(grepl('nophyloscore',method.risk))
	{
		if(verbose)
			cat(paste('\nsetting likelihood to nophyloscore'))
		set(YX.h, NULL, 'score.Y', YX.h[, w.tn])						
	}
	#	don t have t.AnyPos for pot tr. allocate by proportion 
	df.artinfo	<- unique(df.artinfo)		 	
	set(df.artinfo, NULL, 'pt', df.artinfo[, pt] * length(ARTns.i) / (length(ARTns.i)+length(ARTs.i)))
	df.artinfo	<- rbind(df.artinfo, data.table(stage=paste('ART.NotYetFirstSu.',tp,sep=''), pt= ifelse(length(ARTs.i), length(ARTs.i) / (length(ARTns.i)+length(ARTs.i)), 0) ), fill=TRUE)
	if(grepl('pair',method.sample))
	{		
		setnames(df.artinfo, 'pt', 'pt.ideal')
		tmp			<- subset(YX.h, !is.na(stage2) | !is.na(stage3))[, sort(table(as.character(stage))/sum(!is.na(stage)))]
		tmp			<- data.table(stage=names(tmp), pt=as.numeric(tmp))
		if( tmp[, !length(which(grepl('ART.NotYetFirstSu', stage)))] )
			tmp		<- rbind(tmp, data.table(stage=paste('ART.NotYetFirstSu.',tp,sep=''), pt=0))
		df.artinfo	<- merge(df.artinfo, tmp, by='stage')		
		setkey(df.artinfo, pt)
	}
	#stopifnot( df.artinfo[, tail(cumsum(pt),1)==1] )
	#	for Xclu Xseq Xmsm, we don t have NegT etc unless we load the bigmem table X.msm 
	#	avoid this by removing a fraction
	nt.table.h	<- subset(nt.table, stat!='YX')[,{
				realloc.i		<- which(grepl(method.reallocate,factor))
				realloc.sum		<- ifelse(length(realloc.i), sum(nt[realloc.i]), 0)
				tmp				<- nt
				names(tmp)		<- factor
				tmp[realloc.i]	<- 0
				if(grepl('stage=prop',method.sample))
				{
					realloc.n		<- ceiling(realloc.sum*df.artinfo$pt)
					realloc.d		<- sum(realloc.n)-realloc.sum					
					stopifnot(realloc.d<=length(realloc.n))
					realloc.n		<- realloc.n + rev(c( rep(-1, realloc.d), rep(0, length(realloc.n)-realloc.d) ))
					if(!all(realloc.n>=0)) cat(paste(Patient, stat, realloc.n))
					stopifnot(all(realloc.n>=0))
					stopifnot(sum(realloc.n)==realloc.sum)
				}
				if(grepl('stage=sample',method.sample))
					realloc.n		<- as.numeric(rmultinom(1, realloc.sum, df.artinfo$pt))										
				tmp[df.artinfo$stage]	<- tmp[df.artinfo$stage] + realloc.n 
				list(factor=factor, nt=tmp)
			}, by=c('Patient','risk','stat')]
	#	for YX we know easily how many transmission intervals are removed
	tmp				<- merge( unique(subset(nt.table.h, select=c(risk, Patient))), unique(subset(nt.table.h, select=c(risk, factor))), by='risk', allow.cartesian=TRUE)
	stopifnot( nrow(merge(nt.table.h, tmp, all.y=1, by=c('Patient','risk','factor')))==nrow(tmp)*3 )	
	#	for YX we know easily how many transmission intervals are removed
	nt.table.YX		<- YX.h[, 	{
				z				<- table(as.character(stage))
				list(risk='stage', stat='YX', factor=names(z), nt=as.numeric(z))
			}, by='Patient']
	nt.table.YX		<- merge( nt.table.YX, tmp, by=c('Patient','risk','factor'), all.y=TRUE) 		
	set(nt.table.YX, nt.table.YX[, which(is.na(nt))], 'stat', 'YX')
	set(nt.table.YX, nt.table.YX[, which(is.na(nt))], 'nt', 0)
	#
	nt.table.h	<- rbind(nt.table.YX, nt.table.h, use.names=TRUE)
	list(YX.h=YX.h, nt.table.h=nt.table.h)	
}
######################################################################################
project.athena.Fisheretal.Hypo.ReallocTest.getYXetc<- function( YX, nt.table, method.risk, predict.t2inf, t2inf.args, df.all, YXf=NULL, th.starttime=2008.5, th.endtime=2011, t.period=0.125, test.pc.prtr.current= 0.173, method.minLowerUWithNegT=1, method.resolveUAna=0, method.realloc='TestC18m', method.sample= 'stage=prop, y=median, t=mean', verbose=TRUE)
{
	#method.realloc<- 'TestC12m18pc'
	stopifnot(grepl('stage=sample|stage=prop',method.sample))
	stopifnot(grepl('y=sample|y=median|y=mean',method.sample))
	stopifnot(grepl('t=mean|t=sample|t=start',method.sample))
	stopifnot(grepl('Test',method.realloc))
	method.reallocate	<- '^U'
	tmp					<- regmatches(method.realloc,regexpr('Test[[:alnum:]]+', method.realloc))	
	test.repeat			<- as.numeric(regmatches(substring(tmp,6),regexpr('[^m]+', substring(tmp,6)))) / 12	
	test.pc.target		<- as.numeric(substring(regmatches(substring(tmp,6),regexpr('m[0-9]+', substring(tmp,6))),2)) / 100
	test.delay			<- substr(tmp,5,5)
	stopifnot(test.delay%in%c('A','C'), is.finite(test.repeat), is.finite(test.pc.target))
	test.delay			<- ifelse(test.delay=='A',0,1/12)
	test.pc				<- round( (test.pc.target-test.pc.prtr.current)/(1-test.pc.prtr.current), d=3 )
	if(verbose)
	{
		cat(paste('\ntest.repeat=',test.repeat))
		cat(paste('\ntest.pc.target=',test.pc.target))
		cat(paste('\ntest.pc.prtr.current=',test.pc.prtr.current))
		cat(paste('\ntest.pc=',test.pc))		
		cat(paste('\ntest.delay=',test.delay))
	}
	#
	rTest.Realloc<- function(df.tr, method.sample, df.dinfo)
	{		
		setkey(df.dinfo, stage)
		tmp		<- unique(subset(df.dinfo, select=c(stage, pt)))
		setkey(tmp,pt)			
		#	allocate stage
		if(grepl('stage=prop',method.sample))
		{
			tmp		<- merge(tmp, data.table(stage=rep(tmp$stage, ceiling(tmp$pt*nrow(df.tr)))), by='stage')
			tmp		<- tmp[seq_len(nrow(df.tr)),]
		}			
		if(grepl('stage=sample',method.sample))
			tmp		<- merge(tmp, data.table(stage=cut(runif(nrow(df.tr)), breaks=c(-Inf,cumsum(tmp$pt)), labels=tmp$stage)), by='stage')
		tmp			<- cbind(df.tr, tmp)
		#	allocate score by stage
		set(tmp, NULL, 'stage', tmp[, as.character(stage)])
		if(grepl('y=sample',method.sample))
			tmp	<- tmp[ ,	list(Patient=Patient, t.Patient=t.Patient, REALLOC_SCORE_Y_RAW= my.sample(df.dinfo[['score.Y.raw']][ df.dinfo[['stage']]==stage ], length(pt), replace=TRUE)), by='stage']
		if(grepl('y=median',method.sample))
			tmp	<- tmp[ ,	list(Patient=Patient, t.Patient=t.Patient,  REALLOC_SCORE_Y_RAW= rep(median(df.dinfo[['score.Y.raw']][ df.dinfo[['stage']]==stage ]), length(pt))), by='stage']
		if(grepl('y=mean',method.sample))
			tmp	<- tmp[ ,	list(Patient=Patient, t.Patient=t.Patient,  REALLOC_SCORE_Y_RAW= rep(mean(df.dinfo[['score.Y.raw']][ df.dinfo[['stage']]==stage ]), length(pt))), by='stage']		
		setnames(tmp, 'stage', 'REALLOC_STAGE')
		tmp
	}	
	rTest<- function(df.rec, method.sample, test.pc)
	{		
		if(grepl('stage=prop',method.sample))	 
		{
			tmp		<- df.rec[, which(is.na(TESTED))]
			z		<- c( rep(0, ceiling(length(tmp)*(1-test.pc))), rep(1, ceiling(length(tmp)*test.pc)) )
			if(length(z)%%2)
				z	<- c(z, NA)
			z		<- na.omit(as.vector(t(matrix(z, ncol=2))))
			set( df.rec, tmp, 'TESTED', z[seq_along(tmp)] )			
		}
		if(grepl('stage=sample',method.sample))
		{
			tmp		<- df.rec[, which(is.na(TESTED))]
			set( df.rec, tmp, 'TESTED', sample(c(0,1), length(tmp), replace=TRUE, prob=c(1-test.pc, test.pc)) )							
		}		
		df.rec
	}	
	rTest.Time<- function(df.tr, method.sample, test.pc=0.5, test.repeat=1, th.starttime=2008.5, th.endtime=2013.5)
	{		
		ans		<- unique(subset( df.tr, select=c(t.Patient, t.INFECTION_T) ))
		#	first test for every t.Patient after th.starttime
		if(grepl('t=mean',method.sample))
			ans[, REALLOC_T1:=rep(th.starttime+test.repeat/2, nrow(ans))]
		if(grepl('t=start',method.sample))
			ans[, REALLOC_T1:=rep(th.starttime, nrow(ans))]		
		if(grepl('t=sample',method.sample))
			ans[, REALLOC_T1:=runif(nrow(ans), th.starttime-test.repeat/2, th.starttime+test.repeat/2)]
		#	subsequent tests 
		ans		<- ans[, 	{
					tmp	<- c( seq(REALLOC_T1, th.endtime-0.01, test.repeat), Inf )
					list(REALLOC_T1=REALLOC_T1, REALLOC_T=tmp[ which(tmp>=t.INFECTION_T) ][1])					
				}, by='t.Patient']
		ans		<- merge(ans, df.tr, by='t.Patient')
		#	only fraction tests
		ans[, REALLOC_Tc:= floor(REALLOC_T)]
		setkey(ans, t.Patient)
		if(grepl('stage=sample',method.sample))
			tmp	<- unique(ans)[, list(t.Patient=t.Patient, TESTED= sample(c(0,1), length(t.Patient), replace=TRUE, prob=c(1-test.pc, test.pc))), by='REALLOC_Tc']
		if(grepl('stage=prop',method.sample))
			tmp	<- unique(ans)[, {
						z		<- c(rep(0,ceiling(length(t.Patient)*(1-test.pc))), rep(1,ceiling(length(t.Patient)*test.pc)))
						if(length(z)%%2)
							z	<- c(z, NA)
						z		<- na.omit(as.vector(t(matrix(z, ncol=2))))
						list(t.Patient=t.Patient, TESTED= z[seq_along(t.Patient)])		
					}, by='REALLOC_Tc']
		ans		<- merge(ans, subset(tmp, select=c(t.Patient, TESTED)), by='t.Patient')
		ans
	}
	resolve.UAna<- function(df.all, df.select, method.sample, p.UAcond, verbose=TRUE)
	{
		setkey(df.all, Patient)
		df.all	<- unique(df.all)
		tmp3	<- df.all[, which(Patient%in%df.select[, t.Patient] & is.na(isAcute))]
		if(verbose)
			cat(paste('\nFound NA isAcute among transmitters, allocate to Yes or No, n=',length(tmp3)))
		if(grepl('stage=sample',method.sample))
			set(df.all, tmp3, 'isAcute', factor(rbinom(length(tmp3), 1, p.UAcond), levels=c(0,1), labels=c('No','Yes')) )
		if(grepl('stage=prop',method.sample))
		{
			tmp		<- rep(c('No','Yes'), ceiling( c(1-p.UAcond, p.UAcond)*length(tmp3)))
			tmp		<- tmp[seq_along(tmp3)]
			set(df.all, tmp3, 'isAcute', tmp )
		}		
		df.all
	}
	
	#	get proportion of stages to reallocate to
	tp			<- as.numeric(substr(regmatches(method.risk,regexpr('tp[0-9]', method.risk)),3,3))
	tmp			<- subset(YX, grepl('Dt|DA',stage))[, table(stage)]
	df.dinfo	<- subset( data.table(stage=names(tmp), nt=tmp, pt=tmp/sum(tmp)), nt>0 )	
	#	use YXf get set of new scores to sample from
	if(is.null(YXf))
		tmp		<- subset(YX, select=c(stage, Patient, t.Patient, score.Y.raw))
	if(!is.null(YXf))
		tmp		<- subset(YXf, select=c(stage, Patient, t.Patient, score.Y.raw))
	setkey(tmp, Patient, t.Patient, stage)
	tmp			<- unique(tmp)
	set( tmp, NULL, 'stage', tmp[, as.character(stage)] )		#use all Diag intervals, not just those from tp4
	set( tmp, NULL, 'stage', tmp[, paste(substr(stage,1,nchar(stage)-1),tp,sep='')] )
	df.dinfo	<- merge(df.dinfo, tmp, by='stage')
	set(df.dinfo, NULL, 'stage', df.dinfo[, as.character(stage)])	
	#	estimated time of infection of transmitters
	if(verbose)
		cat(paste('\nusing method.resolveUAna=',method.resolveUAna))
	#	get time of infection of transmitters
	#	method.resolveUAna<- 0	
	if(method.resolveUAna)
	{
		tmp		<- subset(YX, grepl(method.reallocate, stage), select=c(t.Patient, t, stage))
		p.UAcond<- subset( tmp, grepl('UA.',stage, fixed=TRUE) | grepl('U.',stage, fixed=TRUE) )[, sum(grepl('UA.',stage, fixed=TRUE))/length(stage)]
		#	of those with isAcute NA, set a proportion p.UAcond to Yes and the rest to No 
		tmp		<- unique(subset(tmp, select=t.Patient))
		tmp2	<- resolve.UAna(df.all, tmp, method.sample, p.UAcond, verbose=verbose)
		df.tr	<- project.athena.Fisheretal.Y.infectiontime(tmp, tmp2, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, method='for.transmitter', method.minLowerUWithNegT=method.minLowerUWithNegT, verbose=FALSE)
	}
	if(!method.resolveUAna)
	{
		tmp			<- unique(subset(YX, grepl(method.reallocate, stage), select=t.Patient))
		df.tr		<- project.athena.Fisheretal.Y.infectiontime(tmp, df.all, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, method='for.transmitter', method.minLowerUWithNegT=method.minLowerUWithNegT, verbose=FALSE)		
	}	
	setkey(df.tr, t.Patient, t)
	df.tr		<- df.tr[, list(t=t[1]), by='t.Patient']
	setnames(df.tr, 't','t.INFECTION_T')
	#	sample test times of transmitters	
	df.tr		<- rTest.Time(df.tr, method.sample, test.pc=test.pc, test.repeat=test.repeat, th.starttime=th.starttime, th.endtime=2013.5)
	if(verbose)
	{
		cat(paste('\nFound transmitters that are tested at baseline and at repeat points, %=', df.tr[, round(mean(TESTED==1),d=3) ]))
	}
	#	
	#	determine which intervals could have been avoided trough testing
	#
	YX.h			<- copy(YX)
	YX.h			<- merge(YX.h, subset(df.tr, select=c(t.Patient, t.INFECTION_T, REALLOC_T1, REALLOC_T, TESTED)), by='t.Patient', all.x=TRUE)	
	if(!method.resolveUAna)
		stopifnot(nrow(subset(YX.h, t.INFECTION_T>t))==0)	
	if(verbose)
		cat(paste('\ntransmission intervals, n=', nrow(YX.h)))
	#	realloc intervals from transmitters that tested positive
	YX.h		<- merge(YX.h, YX.h[, list(t_min=min(t)), by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
	YX.h[, t.grace:= YX.h[,t.INFECTION_T-t_min]]
	set(YX.h, YX.h[, which(!grepl('UAna',stage) | t.grace<0 | is.na(t.grace)) ], 't.grace', 0)
	if(verbose)
		cat(paste('\ntransmission intervals with t.grace>0, n=',YX.h[, length(which(t.grace>0))] ))
	#	determine which intervals to re-alloc
	YX.h[, TR_TEST:= 0]
	#	t.INFECTION_T+test.delay<=REALLOC_T		--> infection time + infection window is before test
	tmp			<- YX.h[, which( 	!is.na(REALLOC_T) & TESTED & t.INFECTION_T+test.delay<=REALLOC_T & 
									t+t.grace>REALLOC_T & grepl(method.reallocate, stage))]
	set(YX.h, tmp, 'TR_TEST', 1)
	if(verbose)
	{		
		cat(paste('\nU transmission intervals, n=', YX.h[, length(which( grepl(method.reallocate, stage)))]	))
		cat(paste('\nU transmission intervals of testing transmitters, n=', YX.h[, length(which( grepl(method.reallocate, stage) & TESTED))]	))
		cat(paste('\nU transmission interval before first test date, n=', YX.h[, length(which( t+t.grace<=REALLOC_T1 & TESTED &  grepl(method.reallocate, stage)))] ))
		cat(paste('\nU transmission interval before test date, n=', YX.h[, length(which( t+t.grace<=REALLOC_T & TESTED &  grepl(method.reallocate, stage)))] ))
		cat(paste('\nU transmission interval after test date, n=', YX.h[, length(which( t+t.grace>REALLOC_T & TESTED &  grepl(method.reallocate, stage)))] ))		
		cat(paste('\nU transmission interval after test date and test is pos, n=', YX.h[, length(which( t.INFECTION_T+test.delay<=REALLOC_T & TESTED & t+t.grace>REALLOC_T &  grepl(method.reallocate, stage)))] ))
		cat(paste('\nU transmission interval after test date and test is neg, n=', YX.h[, length(which( t.INFECTION_T+test.delay>REALLOC_T & TESTED & t+t.grace>REALLOC_T &  grepl(method.reallocate, stage)))] ))
		cat(paste('\nreallocate entries, n=', length(tmp)))
		cat(paste('\nsum of score.Y.raw before reallocate', YX.h[tmp, sum(score.Y.raw)]))		
	}
	#	determine stages and scores to be re-alloced to
	df.tr		<- unique(subset(YX.h, TR_TEST==1, select=c(Patient, t.Patient)))
	df.tr		<- rTest.Realloc(df.tr, method.sample, df.dinfo)
	YX.h		<- merge( YX.h, df.tr, by=c('Patient','t.Patient'), all.x=1 )
	set(YX.h, YX.h[, which(TR_TEST==0)], 'REALLOC_STAGE', NA_character_)
	set(YX.h, YX.h[, which(TR_TEST==0)], 'REALLOC_SCORE_Y_RAW', NA_real_)
	#
	#	add proportion of in and out reallocated stages to df.uinfo
	#
	df.uinfo	<- subset(YX.h, !is.na(REALLOC_STAGE))[, list(REALLOC_NT_IN= length(t)), by='REALLOC_STAGE']
	if(nrow(df.uinfo))
	{
		setnames(df.uinfo, 'REALLOC_STAGE', 'stage')	 
		df.uinfo	<- merge( df.uinfo, YX.h[, list(NT= length(t)), by='stage'], by='stage', all=1 )
		df.uinfo	<- merge( df.uinfo, subset(YX.h, !is.na(REALLOC_STAGE))[, list(REALLOC_NT_OUT= length(t)), by='stage'], by='stage', all=1 )
		set(df.uinfo, df.uinfo[, which(is.na(REALLOC_NT_IN))], 'REALLOC_NT_IN', 0)
		set(df.uinfo, df.uinfo[, which(is.na(REALLOC_NT_OUT))], 'REALLOC_NT_OUT', 0)
		df.uinfo[, U_OUT_AFTER_RM:= REALLOC_NT_OUT/NT]
		set(df.uinfo, df.uinfo[, which(is.na(NT))], c('NT','U_OUT_AFTER_RM'), 0)
		df.uinfo[, D_IN_AFTER_RM:= REALLOC_NT_IN/sum(REALLOC_NT_IN)]		
		setnames(df.uinfo, 'stage', 'factor')
		set(df.uinfo, NULL, 'factor', df.uinfo[, as.character(factor)])
		df.uinfo[, risk:='stage']
		set(df.uinfo, NULL, c('REALLOC_NT_IN','REALLOC_NT_OUT','NT'), NULL)			
	}
	#
	#	reallocate
	#
	tmp			<- YX.h[, which(TR_TEST==1)]
	set(YX.h, tmp, 'stage', YX.h[tmp, REALLOC_STAGE])
	set(YX.h, tmp, c('score.Y','score.Y.raw'), YX.h[tmp, REALLOC_SCORE_Y_RAW])
	set(YX.h, NULL, c('t_min', 't.grace', 't.INFECTION_T', 'REALLOC_T1','REALLOC_T', 'REALLOC_STAGE', 'REALLOC_SCORE_Y_RAW', 'TESTED', 'TR_TEST'), NULL)
	#	update scores
	if(grepl('wtn',method.risk))
	{
		if(verbose)
			cat(paste('\nsetting likelihood to likelihood of pair / number transmission intervals'))
		set(YX.h, NULL, 'score.Y', YX.h[, score.Y.raw*w.tn])				
	}	
	if(grepl('noscore',method.risk))
	{
		if(verbose)
			cat(paste('\nsetting likelihood to 1'))
		set(YX.h, NULL, 'score.Y', YX.h[, 1])						
	}
	if(grepl('nophyloscore',method.risk))
	{
		if(verbose)
			cat(paste('\nsetting likelihood to nophyloscore'))
		set(YX.h, NULL, 'score.Y', YX.h[, w.tn])						
	}
	#
	#	prepare nt.table for hypothetical scenario
	#			
	#	for Xclu Xseq Xmsm, we don t have NegT etc unless we load the bigmem table X.msm 
	#	avoid this by removing a fraction	
	setkey(nt.table, Patient, stat, risk, factor)
	if(!nrow(df.uinfo))
		nt.table.h		<- subset(nt.table, stat!='YX', select=c(Patient, risk, factor, stat, nt))
	if(nrow(df.uinfo))
	{
		setkey(df.uinfo, risk, factor)
		#	move stages according to fractions UPT and DPT
		#	subset(nt.table, Patient=='M41498' & stat=='X.msm')
		if(grepl('stage=sample',method.sample))
			nt.table.h	<- subset(nt.table, stat!='YX')[, {
						flow.out	<- rbinom(length(nt), nt, df.uinfo$U_OUT_AFTER_RM) 
						flow.in		<- as.numeric(rmultinom(1, sum(flow.out), df.uinfo$D_IN_AFTER_RM))
						stopifnot(sum(flow.out)==sum(flow.in))
						#print(nt); print(flow.out); print(flow.in)
						list(factor=factor, nt= nt-flow.out+flow.in)
					}, by=c('Patient','risk','stat')]
		if(grepl('stage=prop',method.sample))
			nt.table.h	<- subset(nt.table, stat!='YX')[, {
						flow.out	<- round( nt * df.uinfo$U_OUT_AFTER_RM )
						flow.in		<- ceiling( sum(flow.out) * df.uinfo$D_IN_AFTER_RM )
						flow.d		<- sum(flow.in)-sum(flow.out)
						stopifnot( flow.d>=0, flow.d<=length(flow.in>0) )
						flow.in[ flow.in>0 ]	<- flow.in[ flow.in>0 ] + c( rep(-1, flow.d), rep(0, length(which(flow.in>0))-flow.d) )
						stopifnot(sum(flow.out)==sum(flow.in))
						list(factor=factor, nt= nt-flow.out+flow.in)					
					}, by=c('Patient','risk','stat')]
	}						
	#	for YX we know easily how many transmission intervals are removed
	nt.table.YX	<- YX.h[, 	{
				z				<- table(as.character(stage))
				list(risk='stage', stat='YX', factor=names(z), nt=as.numeric(z))
			}, by='Patient']
	tmp			<- merge( unique(subset(nt.table.h, select=c(risk, Patient))), unique(subset(nt.table.h, select=c(risk, factor))), by='risk', allow.cartesian=TRUE)
	nt.table.YX	<- merge( nt.table.YX, tmp, by=c('Patient','risk','factor'), all.y=TRUE) 		
	set(nt.table.YX, nt.table.YX[, which(is.na(nt))], 'stat', 'YX')
	set(nt.table.YX, nt.table.YX[, which(is.na(nt))], 'nt', 0)
	#
	nt.table.h	<- rbind(nt.table.YX, nt.table.h, use.names=TRUE)
	stopifnot( nrow(merge(nt.table.h, tmp, all.y=1, by=c('Patient','risk','factor')))==nrow(tmp)*4 )
	nt.table.h	<- dcast.data.table(nt.table.h, Patient+risk+factor~stat, value.var='nt')
	tmp			<- nt.table.h[, which(YX>X.clu)]
	if(length(tmp))
	{
		cat(paste('\nWARNING UtoD Found YX>X.clu, n=',length(tmp)))
		set(nt.table.h, tmp, 'X.clu', nt.table.h[tmp, YX])
	}
	tmp			<- nt.table.h[, which(X.clu>X.seq)]	
	if(length(tmp))
	{
		cat(paste('\nWARNING UtoD Found X.clu>X.seq, n=',length(tmp)))	
		set(nt.table.h, tmp, 'X.seq', nt.table.h[tmp, X.clu])
	}
	tmp			<- nt.table.h[, which(X.seq>X.msm)]	
	if(length(tmp))
	{
		cat(paste('\nWARNING UtoD Found X.seq>X.msm, n=',length(tmp)))
		set(nt.table.h, tmp, 'X.msm', nt.table.h[tmp, X.seq])
	}
	nt.table.h	<- melt(nt.table.h, measure.vars=c('X.clu','X.msm','X.seq','YX'), value.name='nt', variable.name='stat')	
	list(YX.h=YX.h, nt.table.h=nt.table.h)	
}
######################################################################################
#	method.sample<- 't=start, stage=prop, y=median, eff=median'; th.starttime=2008.5; th.endtime=2011;  t.period=0.125; method.minLowerUWithNegT=1; method.resolveUAna=1
project.athena.Fisheretal.Hypo.ReallocPrEPandTestForTarget.getYXetc<- function( YX, nt.table, method.risk, predict.t2inf, t2inf.args, df.all, YXf=NULL, th.starttime=2008.5, th.endtime=2011, t.period=0.125, test.pc.prtr.current= 0.173, method.minLowerUWithNegT=1, method.resolveUAna=0, method.realloc='PrestC18m50pc', method.sample= "t=start, stage=prop, y=median, eff=median", verbose=TRUE)
{
	stopifnot(grepl('stage=sample|stage=prop',method.sample))
	stopifnot(grepl('y=sample|y=median|y=mean',method.sample))
	stopifnot(grepl('t=mean|t=sample|t=start',method.sample))
	stopifnot(grepl('eff=sample|eff=median',method.sample))
	stopifnot(grepl('Trest',method.realloc))
	method.reallocate	<- '^U'
	tmp					<- regmatches(method.realloc,regexpr('Trest[[:alnum:]]+', method.realloc))
	prep.eff			<- substr(tmp,6,10)
	prest.delay			<- substr(tmp,11,11)
	prest.repeat		<- as.numeric(substr(tmp,12,13)) / 12
	test.pc.target		<- regmatches(method.realloc,regexpr('m[0-9]+pc', tmp))
	test.pc.target		<- as.numeric(substr(test.pc.target, 2, nchar(test.pc.target)-2)) / 100
	prep.pc				<- regmatches(method.realloc,regexpr('c[0-9]+pc', tmp))
	prep.pc				<- as.numeric(substr(prep.pc, 2, nchar(prep.pc)-2)) / 100
	prep.tgt			<- regmatches(method.realloc,regexpr('[0-9]+y$', tmp))
	prep.tgt			<- as.numeric(substr(prep.tgt, 1, nchar(prep.tgt)-1))
	stopifnot(prep.eff%in%c('PROUD','IPrEX'), prest.delay%in%c('A','C'), is.finite(prep.tgt), is.finite(prest.repeat), is.finite(test.pc.target), is.finite(prep.pc))
	prest.delay			<- ifelse(prest.delay=='A',0,1/12)	
	#	set test.pc for transmitters
	test.pc				<- round( (test.pc.target-test.pc.prtr.current)/(1-test.pc.prtr.current), d=3 )
	if(verbose)
	{
		cat(paste('\nprest.repeat=',prest.repeat))
		cat(paste('\ntest.pc.target=',test.pc.target))
		cat(paste('\ntest.pc.prtr.current=',test.pc.prtr.current))
		cat(paste('\ntest.pc=',test.pc))
		cat(paste('\nprep.pc=',prep.pc))
		cat(paste('\nprest.delay=',prest.delay))
		cat(paste('\nprep.eff=',prep.eff))
		cat(paste('\nprep.tgt=',prep.tgt))
	}
	#
	rTest.Realloc<- function(df.tr, method.sample, df.dinfo)
	{		
		setkey(df.dinfo, stage)
		tmp		<- unique(subset(df.dinfo, select=c(stage, pt)))
		setkey(tmp,pt)			
		#	allocate stage
		if(grepl('stage=prop',method.sample))
		{
			tmp		<- merge(tmp, data.table(stage=rep(tmp$stage, ceiling(tmp$pt*nrow(df.tr)))), by='stage')
			tmp		<- tmp[seq_len(nrow(df.tr)),]
		}			
		if(grepl('stage=sample',method.sample))
			tmp		<- merge(tmp, data.table(stage=cut(runif(nrow(df.tr)), breaks=c(-Inf,cumsum(tmp$pt)), labels=tmp$stage)), by='stage')
		tmp			<- cbind(df.tr, tmp)
		#	allocate score by stage
		set(tmp, NULL, 'stage', tmp[, as.character(stage)])
		if(grepl('y=sample',method.sample))
			tmp	<- tmp[ ,	list(Patient=Patient, t.Patient=t.Patient, REALLOC_SCORE_Y_RAW= my.sample(df.dinfo[['score.Y.raw']][ df.dinfo[['stage']]==stage ], length(pt), replace=TRUE)), by='stage']
		if(grepl('y=median',method.sample))
			tmp	<- tmp[ ,	list(Patient=Patient, t.Patient=t.Patient,  REALLOC_SCORE_Y_RAW= rep(median(df.dinfo[['score.Y.raw']][ df.dinfo[['stage']]==stage ]), length(pt))), by='stage']
		if(grepl('y=mean',method.sample))
			tmp	<- tmp[ ,	list(Patient=Patient, t.Patient=t.Patient,  REALLOC_SCORE_Y_RAW= rep(mean(df.dinfo[['score.Y.raw']][ df.dinfo[['stage']]==stage ]), length(pt))), by='stage']		
		setnames(tmp, 'stage', 'REALLOC_STAGE')
		tmp
	}	
	rPrEP<- function(df.rec, method.sample, prep.pc, prep.eff='IPrEX')
	{		
		if(grepl('stage=prop',method.sample))	 
		{
			tmp		<- df.rec[, which(is.na(PRESTED))]
			z		<- c( rep(0, ceiling(length(tmp)*(1-prep.pc))), rep(1, ceiling(length(tmp)*prep.pc)) )
			if(length(z)%%2)
				z	<- c(z, NA)
			z		<- na.omit(as.vector(t(matrix(z, ncol=2))))
			set( df.rec, tmp, 'PRESTED', z[seq_along(tmp)] )			
		}
		if(grepl('stage=sample',method.sample))
		{
			tmp		<- df.rec[, which(is.na(PRESTED))]
			set( df.rec, tmp, 'PRESTED', sample(c(0,1), length(tmp), replace=TRUE, prob=c(1-prep.pc, prep.pc)) )							
		}		
		if(!any('PREP_EFF'==names(df.rec)))
			df.rec[, PREP_EFF:=NA_real_]
		if(grepl('eff=median', method.sample))
		{
			if(!any('PREP_EFF_PAR'==names(df.rec)) & prep.eff=='IPrEX' )
				df.rec[, PREP_EFF_PAR:=0.44]
			if(!any('PREP_EFF_PAR'==names(df.rec)) & prep.eff=='PROUD' )
				df.rec[, PREP_EFF_PAR:=0.86]			
			tmp		<- df.rec[, which(is.na(PREP_EFF))]				
			z		<- c( rep(0, ceiling(length(tmp)*(1-df.rec[1, PREP_EFF_PAR]))), rep(1, ceiling(length(tmp)*df.rec[1, PREP_EFF_PAR])) )
			set(df.rec, tmp, 'PREP_EFF', z[seq_along(tmp)])			
		}
		if(grepl('eff=sample', method.sample))
		{			
			if(!any('PREP_EFF_PAR'==names(df.rec)) & prep.eff=='IPrEX' )
				df.rec[, PREP_EFF_PAR:=rbeta(1, 6, (1-0.44)/0.44*6 )]
			if(!any('PREP_EFF_PAR'==names(df.rec)) & prep.eff=='PROUD' )
				df.rec[, PREP_EFF_PAR:=rbeta(1, 2.9, (1-0.86)/0.86*2.9 )]			
			tmp		<- df.rec[, which(is.na(PREP_EFF))]
			set(df.rec, tmp, 'PREP_EFF', sample(c(0,1), length(tmp), prob=c(1-df.rec[1, PREP_EFF_PAR], df.rec[1, PREP_EFF_PAR]), replace=TRUE) )
		}
		df.rec
	}	
	rTest.Time<- function(df.tr, method.sample, test.pc=0.5, prest.repeat=1, th.starttime=2008.5, th.endtime=2013.5)
	{		
		ans		<- unique(subset( df.tr, select=c(t.Patient, t.INFECTION_T) ))
		#	first test for every t.Patient after th.starttime
		if(grepl('t=mean',method.sample))
			ans[, REALLOC_T1:=rep(th.starttime+prest.repeat/2, nrow(ans))]
		if(grepl('t=start',method.sample))
			ans[, REALLOC_T1:=rep(th.starttime, nrow(ans))]		
		if(grepl('t=sample',method.sample))
			ans[, REALLOC_T1:=runif(nrow(ans), th.starttime-prest.repeat/2, th.starttime+prest.repeat/2)]
		#	subsequent tests 
		ans		<- ans[, 	{
					tmp	<- c( seq(REALLOC_T1, th.endtime-0.01, prest.repeat), Inf )
					list(REALLOC_T1=REALLOC_T1, REALLOC_T=tmp[ which(tmp>=t.INFECTION_T) ][1])					
				}, by='t.Patient']
		ans		<- merge(ans, df.tr, by='t.Patient')
		#	only fraction tests
		ans[, REALLOC_Tc:= floor(REALLOC_T)]
		setkey(ans, t.Patient)
		if(grepl('stage=sample',method.sample))
			tmp	<- unique(ans)[, list(t.Patient=t.Patient, PRESTED= sample(c(0,1), length(t.Patient), replace=TRUE, prob=c(1-test.pc, test.pc))), by='REALLOC_Tc']
		if(grepl('stage=prop',method.sample))
			tmp	<- unique(ans)[, {
						z		<- c(rep(0,ceiling(length(t.Patient)*(1-test.pc))), rep(1,ceiling(length(t.Patient)*test.pc)))
						if(length(z)%%2)
							z	<- c(z, NA)
						z		<- na.omit(as.vector(t(matrix(z, ncol=2))))
						list(t.Patient=t.Patient, PRESTED= z[seq_along(t.Patient)])		
					}, by='REALLOC_Tc']
		ans		<- merge(ans, subset(tmp, select=c(t.Patient, PRESTED)), by='t.Patient')
		ans
	}
	resolve.UAna<- function(df.all, df.select, method.sample, p.UAcond, verbose=TRUE)
	{
		setkey(df.all, Patient)
		df.all	<- unique(df.all)
		tmp3	<- df.all[, which(Patient%in%df.select[, t.Patient] & is.na(isAcute))]
		if(verbose)
			cat(paste('\nFound NA isAcute among transmitters, allocate to Yes or No, n=',length(tmp3)))
		if(grepl('stage=sample',method.sample))
			set(df.all, tmp3, 'isAcute', factor(rbinom(length(tmp3), 1, p.UAcond), levels=c(0,1), labels=c('No','Yes')) )
		if(grepl('stage=prop',method.sample))
		{
			tmp		<- rep(c('No','Yes'), ceiling( c(1-p.UAcond, p.UAcond)*length(tmp3)))
			tmp		<- tmp[seq_along(tmp3)]
			set(df.all, tmp3, 'isAcute', tmp )
		}		
		df.all
	}
	
	#	get proportion of stages to reallocate to
	tp			<- as.numeric(substr(regmatches(method.risk,regexpr('tp[0-9]', method.risk)),3,3))
	tmp			<- subset(YX, grepl('Dt|DA',stage))[, table(stage)]
	df.dinfo	<- subset( data.table(stage=names(tmp), nt=tmp, pt=tmp/sum(tmp)), nt>0 )	
	#	use YXf get set of new scores to sample from
	if(is.null(YXf))
		tmp		<- subset(YX, select=c(stage, Patient, t.Patient, score.Y.raw))
	if(!is.null(YXf))
		tmp		<- subset(YXf, select=c(stage, Patient, t.Patient, score.Y.raw))
	setkey(tmp, Patient, t.Patient, stage)
	tmp			<- unique(tmp)
	set( tmp, NULL, 'stage', tmp[, as.character(stage)] )		#use all Diag intervals, not just those from tp4
	set( tmp, NULL, 'stage', tmp[, paste(substr(stage,1,nchar(stage)-1),tp,sep='')] )
	df.dinfo	<- merge(df.dinfo, tmp, by='stage')
	set(df.dinfo, NULL, 'stage', df.dinfo[, as.character(stage)])	
	#	estimated time of infection of transmitters
	if(verbose)
		cat(paste('\nusing method.resolveUAna=',method.resolveUAna))
	#	get time of infection of transmitters
	#	method.resolveUAna<- 0	
	if(method.resolveUAna)
	{
		tmp		<- subset(YX, grepl(method.reallocate, stage), select=c(t.Patient, t, stage))
		p.UAcond<- subset( tmp, grepl('UA.',stage, fixed=TRUE) | grepl('U.',stage, fixed=TRUE) )[, sum(grepl('UA.',stage, fixed=TRUE))/length(stage)]
		#	of those with isAcute NA, set a proportion p.UAcond to Yes and the rest to No 
		tmp		<- unique(subset(tmp, select=t.Patient))
		tmp2	<- resolve.UAna(df.all, tmp, method.sample, p.UAcond, verbose=verbose)
		df.tr	<- project.athena.Fisheretal.Y.infectiontime(tmp, tmp2, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, method='for.transmitter', method.minLowerUWithNegT=method.minLowerUWithNegT, verbose=FALSE)
	}
	if(!method.resolveUAna)
	{
		tmp			<- unique(subset(YX, grepl(method.reallocate, stage), select=t.Patient))
		df.tr		<- project.athena.Fisheretal.Y.infectiontime(tmp, df.all, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, method='for.transmitter', method.minLowerUWithNegT=method.minLowerUWithNegT, verbose=FALSE)		
	}	
	setkey(df.tr, t.Patient, t)
	df.tr		<- df.tr[, list(t=t[1]), by='t.Patient']
	setnames(df.tr, 't','t.INFECTION_T')
	#	sample test times of transmitters	
	df.tr		<- rTest.Time(df.tr, method.sample, test.pc=test.pc, prest.repeat=prest.repeat, th.starttime=th.starttime, th.endtime=2013.5)
	#	those tested are offered PrEP, which may be efficacious from the first test time onwards
	#	so PrEP improves on testing
	#
	tmp			<- subset(df.tr, select=c(t.Patient, PRESTED))
	setnames(tmp, 't.Patient', 'Patient')
	tmp			<- rPrEP(tmp, method.sample, prep.pc, prep.eff=prep.eff)		#prep.pc does not do anything here because all transmitters have !NA PRESTED
	setnames(tmp, 'Patient', 't.Patient')
	df.tr		<- merge(df.tr, unique(subset(tmp, select=c(t.Patient, PREP_EFF, PREP_EFF_PAR))), by='t.Patient')
	#	use REALLOC_T1 to determine if first test before infection time. REALLOC_T is always first test AFTER infection time 
	df.tr[, TR_PREP_ON:= as.numeric(PRESTED & PREP_EFF & t.INFECTION_T+prest.delay>REALLOC_T1)]
	#	set PREP_EFF to 0 if t.Patient older than prep.tgt at first presting time 
	tmp			<- unique(subset(df.all, select=c(Patient, DateBorn)))
	setnames(tmp, 'Patient', 't.Patient')
	df.tr		<- merge(df.tr, tmp, by='t.Patient')
	set(df.tr, df.tr[, which(REALLOC_T1-DateBorn>=prep.tgt)], 'TR_PREP_ON', 0)	
	if(verbose)
	{
		cat(paste('\nFound transmitters that are tested at baseline and at repeat points, n=', df.tr[, round(mean(PRESTED==1),d=3) ]))
		cat(paste('\nFound transmitters with hyp effective PrEP, n=', df.tr[, round(mean(PREP_EFF==1),d=3) ]))
		cat(paste('\nFound transmitters that could be on PrEP if inf time is after neg test, n=', df.tr[, round(mean(PRESTED & PREP_EFF),d=3) ]))
		cat(paste('\nFound transmitters that are on PrEP, n=', df.tr[, round(mean(TR_PREP_ON==1),d=3) ]))
	}
	#
	#	some of the transmitters are also recipients: offer PrEP consistently	
	#
#tmp	<- sapply(1:1000, function(j){
	df.rec		<- unique(subset(nt.table, select=Patient))
	tmp			<- unique(subset(df.tr, select=c(t.Patient, PRESTED, PREP_EFF, PREP_EFF_PAR)))
	setnames(tmp, 't.Patient','Patient')
	df.rec		<- merge(df.rec, tmp, by='Patient', all.x=1)
	if(df.rec[,!all(is.na(PREP_EFF_PAR))])
		set(df.rec, df.rec[, which(is.na(PREP_EFF_PAR))], 'PREP_EFF_PAR', df.rec[which(!is.na(PREP_EFF_PAR))[1],PREP_EFF_PAR])
	if(df.rec[,all(is.na(PREP_EFF_PAR))])
		set(df.rec, NULL, 'PREP_EFF_PAR', NULL)	
	#	for those that are no transmitter, allocate if Prested and if PREP efficacious
	df.rec		<- rPrEP(df.rec, method.sample, prep.pc, prep.eff=prep.eff)
	df.rec[, REC_PREP_ON:= as.numeric(PRESTED & PREP_EFF)]
	#	set PREP_EFF to 0 if t.Patient older than prep.tgt at first presting time
	if(df.rec[1,grepl('_bs[0-9]+',Patient)])
	{
		df.rec[, Patient.b4bs:= sapply(strsplit(Patient,'_'), '[[', 1)]
		tmp			<- unique(subset(df.all, select=c(Patient, DateBorn, AnyPos_T1)))
		setnames(tmp, 'Patient','Patient.b4bs')
		df.rec		<- merge(df.rec, tmp, by='Patient.b4bs')
		set(df.rec, NULL, 'Patient.b4bs', NULL)
	}
	if(!df.rec[1,grepl('_bs[0-9]+',Patient)])
	{
		tmp			<- unique(subset(df.all, select=c(Patient, DateBorn, AnyPos_T1)))	
		df.rec		<- merge(df.rec, tmp, by='Patient')		
	}
	set(df.rec, df.rec[, which(AnyPos_T1-1-DateBorn>=prep.tgt)], 'REC_PREP_ON', 0)		
	if(verbose)
	{
		cat(paste('\nFound recipients that are offered PrEP at baseline, n=', df.rec[, round(length(which(PRESTED==1))/length(PRESTED),d=3) ]))
		cat(paste('\nFound recipients with hyp effective PrEP, n=', df.rec[, round(length(which(PREP_EFF==1))/length(PRESTED),d=3) ]))
		cat(paste('\nFound recipients that are on PrEP, n=', df.rec[, round(length(which(REC_PREP_ON==1))/length(PRESTED),d=3) ]))
	}
	#df.rec[, length(which(REC_PREP_ON=='1'))/length(REC_PREP_ON)]
#})#	summary(tmp)
	#
	#	combine test of transmitters, PrEP of transmitters and PrEP of recipients  
	#
	
	#	
	#	determine which intervals could have been avoided trough testing
	YX.h			<- copy(YX)
	YX.h			<- merge(YX.h, subset(df.tr, select=c(t.Patient, t.INFECTION_T, REALLOC_T1, REALLOC_T, PRESTED, TR_PREP_ON)), by='t.Patient', all.x=TRUE)
	YX.h			<- merge(YX.h, subset(df.rec, select=c(Patient, REC_PREP_ON)), by='Patient', all.x=TRUE)
	if(!method.resolveUAna)
		stopifnot(nrow(subset(YX.h, t.INFECTION_T>t))==0)	
	if(verbose)
		cat(paste('\ntransmission intervals, n=', nrow(YX.h)))	
	#	remove intervals to recipients that are on effective PrEP
	YX.h			<- subset(YX.h, REC_PREP_ON==0)
	if(verbose)
		cat(paste('\ntransmission intervals to recipients not on eff PrEP, n=', nrow(YX.h)))
	#	proportion of transmission intervals for which PrEP works, after we rm recipient on PrEP
	#	need this to reduce nt.table
	df.uinfo		<- merge( subset(YX.h, TR_PREP_ON==1)[, list(TR_RM= length(t)), by='stage'], YX.h[, list(N= length(t)), by='stage'], by='stage', all.y=1)
	set(df.uinfo, df.uinfo[, which(is.na(TR_RM))], 'TR_RM', 0)
	df.uinfo[, TR_NRM:= 1-TR_RM/N]
	set(df.uinfo, NULL, c('TR_RM','N'), NULL)
	#	remove intervals from transmitters that are on effective PrEP
	YX.h			<- subset(YX.h, is.na(TR_PREP_ON) | TR_PREP_ON==0)
	if(verbose)
		cat(paste('\ntransmission intervals from transmitters not on eff PrEP, n=', nrow(YX.h)))
	#	realloc intervals from transmitters that tested positive
	YX.h		<- merge(YX.h, YX.h[, list(t_min=min(t)), by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
	set(YX.h, NULL, 't.grace', YX.h[,t.INFECTION_T-t_min])
	set(YX.h, YX.h[, which(!grepl('UAna',stage) | t.grace<0 | is.na(t.grace)) ], 't.grace', 0)
	if(verbose)
		cat(paste('\ntransmission intervals with t.grace>0, n=', YX.h[, length(which(t.grace>0))]))	
	#	re-alloc intervals of tested transmitters to diag
	#	t.INFECTION_T+prest.delay<=REALLOC_T	--> infected before pos test
	#	t+t.grace>REALLOC_T						--> intervals after pos test
	YX.h[, TR_TEST:= 0]
	tmp			<- YX.h[, which( 	!is.na(REALLOC_T) & PRESTED & t.INFECTION_T+prest.delay<=REALLOC_T & 
									t+t.grace>REALLOC_T & grepl(method.reallocate, stage))]
	set(YX.h, tmp, 'TR_TEST', 1)		
	if(verbose)
		cat(paste('\nrealloc transmission intervals from transmitters that got tested pos, n=', length(tmp)))
	tmp			<- unique(subset(YX.h, TR_TEST==1, select=c(Patient, t.Patient)))
	tmp			<- rTest.Realloc(tmp, method.sample, df.dinfo)
	YX.h		<- merge( YX.h, tmp, by=c('Patient','t.Patient'), all.x=1 )
	set(YX.h, YX.h[, which(TR_TEST==0)], 'REALLOC_STAGE', NA_character_)
	set(YX.h, YX.h[, which(TR_TEST==0)], 'REALLOC_SCORE_Y_RAW', NA_real_)
	#
	#	add proportion of in and out reallocated stages to df.uinfo
	#
	tmp			<- subset(YX.h, !is.na(REALLOC_STAGE))[, list(REALLOC_NT_IN= length(t)), by='REALLOC_STAGE']	
	setnames(tmp, 'REALLOC_STAGE', 'stage')
	if(!nrow(tmp))
		set(df.uinfo,NULL, c('U_OUT_AFTER_RM','D_IN_AFTER_RM'),0)
	if(nrow(tmp))
	{
		df.uinfo	<- merge(df.uinfo, tmp, by='stage', all.x=1)
		df.uinfo	<- merge( df.uinfo, YX.h[, list(NT= length(t)), by='stage'], by='stage', all=1 )
		df.uinfo	<- merge( df.uinfo, subset(YX.h, !is.na(REALLOC_STAGE))[, list(REALLOC_NT_OUT= length(t)), by='stage'], by='stage', all=1 )
		set(df.uinfo, df.uinfo[, which(is.na(REALLOC_NT_IN))], 'REALLOC_NT_IN', 0)
		set(df.uinfo, df.uinfo[, which(is.na(REALLOC_NT_OUT))], 'REALLOC_NT_OUT', 0)
		df.uinfo[, U_OUT_AFTER_RM:= REALLOC_NT_OUT/NT]
		set(df.uinfo, df.uinfo[, which(is.na(NT))], c('NT','U_OUT_AFTER_RM'), 0)
		df.uinfo[, D_IN_AFTER_RM:= REALLOC_NT_IN/sum(REALLOC_NT_IN)]		
		set(df.uinfo, NULL, c('REALLOC_NT_IN','REALLOC_NT_OUT','NT'), NULL)
	}
	setnames(df.uinfo, 'stage', 'factor')
	set(df.uinfo, NULL, 'factor', df.uinfo[, as.character(factor)])
	df.uinfo[, risk:='stage']		
	#
	#	reallocate
	#
	tmp			<- YX.h[, which(TR_TEST==1)]
	set(YX.h, tmp, 'stage', YX.h[tmp, REALLOC_STAGE])
	set(YX.h, tmp, c('score.Y','score.Y.raw'), YX.h[tmp, REALLOC_SCORE_Y_RAW])
	set(YX.h, NULL, c('t_min', 't.grace', 't.INFECTION_T', 'REALLOC_T1','REALLOC_T', 'REALLOC_STAGE', 'REALLOC_SCORE_Y_RAW', 'PRESTED', 'TR_PREP_ON', 'REC_PREP_ON', 'TR_TEST'), NULL)
	#	update scores
	if(grepl('wtn',method.risk))
	{
		if(verbose)
			cat(paste('\nsetting likelihood to likelihood of pair / number transmission intervals'))
		set(YX.h, NULL, 'score.Y', YX.h[, score.Y.raw*w.tn])				
	}									
	if(grepl('noscore',method.risk))
	{
		if(verbose)
			cat(paste('\nsetting likelihood to 1'))
		set(YX.h, NULL, 'score.Y', YX.h[, 1])						
	}
	if(grepl('nophyloscore',method.risk))
	{
		if(verbose)
			cat(paste('\nsetting likelihood to nophyloscore'))
		set(YX.h, NULL, 'score.Y', YX.h[, w.tn])						
	}
	#
	#	prepare nt.table for hypothetical scenario
	#			
	#	for Xclu Xseq Xmsm, we don t have NegT etc unless we load the bigmem table X.msm 
	#	avoid this by removing a fraction
	setkey(df.uinfo, risk, factor)
	setkey(nt.table, Patient, stat, risk, factor)
	#	remove recipients that are on PrEP	
	#	Xmsm includes participants with no transmitter, so need to use df.rec
	nt.table.h		<- merge( nt.table, subset(df.rec, !REC_PREP_ON, Patient), by='Patient' )	
	#	
	#	remove transmitters that are on PrEP and thereafter, move stages according to fractions UPT and DPT
	#
	#	subset(nt.table, Patient=='M41498' & stat=='X.msm')
	if(df.uinfo[, all(U_OUT_AFTER_RM==0)] & grepl('stage=sample',method.sample))
		nt.table.h	<- subset(nt.table.h, stat!='YX')[, {
					nt.nrm		<- rbinom(length(nt), nt, df.uinfo$TR_NRM)
					list(factor=factor, nt= nt.nrm)
				}, by=c('Patient','risk','stat')]
	if(df.uinfo[, all(U_OUT_AFTER_RM==0)] & grepl('stage=prop',method.sample))
		nt.table.h	<- subset(nt.table.h, stat!='YX')[, {
					nt.nrm		<- round( nt * df.uinfo$TR_NRM)
					list(factor=factor, nt= nt.nrm)					
				}, by=c('Patient','risk','stat')]	
	if(df.uinfo[, !all(U_OUT_AFTER_RM==0)] & grepl('stage=sample',method.sample))
		nt.table.h	<- subset(nt.table.h, stat!='YX')[, {
					nt.nrm		<- rbinom(length(nt), nt, df.uinfo$TR_NRM)
					flow.out	<- rbinom(length(nt.nrm), nt.nrm, df.uinfo$U_OUT_AFTER_RM) 
					flow.in		<- as.numeric(rmultinom(1, sum(flow.out), df.uinfo$D_IN_AFTER_RM))
					stopifnot(sum(flow.out)==sum(flow.in))
					#print(nt); print(flow.out); print(flow.in)
					list(factor=factor, nt= nt.nrm-flow.out+flow.in)
				}, by=c('Patient','risk','stat')]
	if(df.uinfo[, !all(U_OUT_AFTER_RM==0)] & grepl('stage=prop',method.sample))
		nt.table.h	<- subset(nt.table.h, stat!='YX')[, {
					nt.nrm		<- round( nt * df.uinfo$TR_NRM)
					flow.out	<- round( nt.nrm * df.uinfo$U_OUT_AFTER_RM )
					flow.in		<- ceiling( sum(flow.out) * df.uinfo$D_IN_AFTER_RM )
					flow.d		<- sum(flow.in)-sum(flow.out)
					stopifnot( flow.d>=0, flow.d<=length(flow.in>0) )
					flow.in[ flow.in>0 ]	<- flow.in[ flow.in>0 ] + c( rep(-1, flow.d), rep(0, length(which(flow.in>0))-flow.d) )
					stopifnot(sum(flow.out)==sum(flow.in))
					list(factor=factor, nt= nt.nrm-flow.out+flow.in)					
				}, by=c('Patient','risk','stat')]		
	#	
	#	for YX, we know easily how many transmission intervals are removed / reallocated
	nt.table.YX	<- YX.h[, 	{
				z				<- table(as.character(stage))
				list(risk='stage', stat='YX', factor=names(z), nt=as.numeric(z))
			}, by='Patient']
	tmp			<- merge( unique(subset(nt.table.h, select=c(risk, Patient))), unique(subset(nt.table.h, select=c(risk, factor))), by='risk', allow.cartesian=TRUE)
	nt.table.YX	<- merge( nt.table.YX, tmp, by=c('Patient','risk','factor'), all.y=TRUE) 		
	set(nt.table.YX, nt.table.YX[, which(is.na(nt))], 'stat', 'YX')
	set(nt.table.YX, nt.table.YX[, which(is.na(nt))], 'nt', 0)
	#
	nt.table.h	<- rbind(nt.table.YX, nt.table.h, use.names=TRUE)
	stopifnot( nrow(merge(nt.table.h, tmp, all.y=1, by=c('Patient','risk','factor')))==nrow(tmp)*4 )
	nt.table.h	<- dcast.data.table(nt.table.h, Patient+risk+factor~stat, value.var='nt')
	tmp			<- nt.table.h[, which(YX>X.clu)]
	if(length(tmp))
	{
		cat(paste('\nWARNING UtoD Found YX>X.clu, n=',length(tmp)))
		set(nt.table.h, tmp, 'X.clu', nt.table.h[tmp, YX])
	}
	tmp			<- nt.table.h[, which(X.clu>X.seq)]	
	if(length(tmp))
	{
		cat(paste('\nWARNING UtoD Found X.clu>X.seq, n=',length(tmp)))	
		set(nt.table.h, tmp, 'X.seq', nt.table.h[tmp, X.clu])
	}
	tmp			<- nt.table.h[, which(X.seq>X.msm)]	
	if(length(tmp))
	{
		cat(paste('\nWARNING UtoD Found X.seq>X.msm, n=',length(tmp)))
		set(nt.table.h, tmp, 'X.msm', nt.table.h[tmp, X.seq])
	}	
	nt.table.h	<- melt(nt.table.h, measure.vars=c('X.clu','X.msm','X.seq','YX'), value.name='nt', variable.name='stat')
	#	
	list(YX.h=YX.h, nt.table.h=nt.table.h)	 
}
######################################################################################
#	method.sample<- 't=start, stage=prop, y=median, eff=median'; th.starttime=2008.5; th.endtime=2011;  t.period=0.125; method.minLowerUWithNegT=1; method.resolveUAna=1
project.athena.Fisheretal.Hypo.ReallocPrEPandTest.getYXetc<- function( YX, nt.table, method.risk, predict.t2inf, t2inf.args, df.all, YXf=NULL, th.starttime=2008.5, th.endtime=2011, t.period=0.125, test.pc.prtr.current= 0.173, method.minLowerUWithNegT=1, method.resolveUAna=0, method.realloc='PrestC18m50pc', method.sample= "t=start, stage=prop, y=median, eff=median", verbose=TRUE)
{
	stopifnot(grepl('stage=sample|stage=prop',method.sample))
	stopifnot(grepl('y=sample|y=median|y=mean',method.sample))
	stopifnot(grepl('t=mean|t=sample|t=start',method.sample))
	stopifnot(grepl('eff=sample|eff=median',method.sample))
	stopifnot(grepl('Prest',method.realloc))
	method.reallocate	<- '^U'
	tmp					<- regmatches(method.realloc,regexpr('Prest[[:alnum:]]+', method.realloc))
	prep.eff			<- substr(tmp,6,10)
	prest.delay			<- substr(tmp,11,11)
	prest.repeat		<- as.numeric(substr(tmp,12,13)) / 12
	test.pc.target		<- regmatches(method.realloc,regexpr('m[0-9]+pc', tmp))
	test.pc.target		<- as.numeric(substr(test.pc.target, 2, nchar(test.pc.target)-2)) / 100
	prep.pc				<- regmatches(method.realloc,regexpr('c[0-9]+pc', tmp))
	prep.pc				<- as.numeric(substr(prep.pc, 2, nchar(prep.pc)-2)) / 100
	stopifnot(prep.eff%in%c('PROUD','IPrEX'), prest.delay%in%c('A','C'), is.finite(prest.repeat), is.finite(test.pc.target), is.finite(prep.pc))
	prest.delay			<- ifelse(prest.delay=='A',0,1/12)	
	#	set test.pc for transmitters
	test.pc				<- round( (test.pc.target-test.pc.prtr.current)/(1-test.pc.prtr.current), d=3 )
	if(verbose)
	{
		cat(paste('\nprest.repeat=',prest.repeat))
		cat(paste('\ntest.pc.target=',test.pc.target))
		cat(paste('\ntest.pc.prtr.current=',test.pc.prtr.current))
		cat(paste('\ntest.pc=',test.pc))
		cat(paste('\nprep.pc=',prep.pc))
		cat(paste('\nprest.delay=',prest.delay))
		cat(paste('\nprep.eff=',prep.eff))
	}
	#
	rTest.Realloc<- function(df.tr, method.sample, df.dinfo)
	{		
		setkey(df.dinfo, stage)
		tmp		<- unique(subset(df.dinfo, select=c(stage, pt)))
		setkey(tmp,pt)			
		#	allocate stage
		if(grepl('stage=prop',method.sample))
		{
			tmp		<- merge(tmp, data.table(stage=rep(tmp$stage, ceiling(tmp$pt*nrow(df.tr)))), by='stage')
			tmp		<- tmp[seq_len(nrow(df.tr)),]
		}			
		if(grepl('stage=sample',method.sample))
			tmp		<- merge(tmp, data.table(stage=cut(runif(nrow(df.tr)), breaks=c(-Inf,cumsum(tmp$pt)), labels=tmp$stage)), by='stage')
		tmp			<- cbind(df.tr, tmp)
		#	allocate score by stage
		set(tmp, NULL, 'stage', tmp[, as.character(stage)])
		if(grepl('y=sample',method.sample))
			tmp	<- tmp[ ,	list(Patient=Patient, t.Patient=t.Patient, REALLOC_SCORE_Y_RAW= my.sample(df.dinfo[['score.Y.raw']][ df.dinfo[['stage']]==stage ], length(pt), replace=TRUE)), by='stage']
		if(grepl('y=median',method.sample))
			tmp	<- tmp[ ,	list(Patient=Patient, t.Patient=t.Patient,  REALLOC_SCORE_Y_RAW= rep(median(df.dinfo[['score.Y.raw']][ df.dinfo[['stage']]==stage ]), length(pt))), by='stage']
		if(grepl('y=mean',method.sample))
			tmp	<- tmp[ ,	list(Patient=Patient, t.Patient=t.Patient,  REALLOC_SCORE_Y_RAW= rep(mean(df.dinfo[['score.Y.raw']][ df.dinfo[['stage']]==stage ]), length(pt))), by='stage']		
		setnames(tmp, 'stage', 'REALLOC_STAGE')
		tmp
	}	
	rPrEP<- function(df.rec, method.sample, prep.pc, prep.eff='IPrEX')
	{		
		if(grepl('stage=prop',method.sample))	 
		{
			tmp		<- df.rec[, which(is.na(PRESTED))]
			z		<- c( rep(0, ceiling(length(tmp)*(1-prep.pc))), rep(1, ceiling(length(tmp)*prep.pc)) )
			if(length(z)%%2)
				z	<- c(z, NA)
			z		<- na.omit(as.vector(t(matrix(z, ncol=2))))
			set( df.rec, tmp, 'PRESTED', z[seq_along(tmp)] )			
		}
		if(grepl('stage=sample',method.sample))
		{
			tmp		<- df.rec[, which(is.na(PRESTED))]
			set( df.rec, tmp, 'PRESTED', sample(c(0,1), length(tmp), replace=TRUE, prob=c(1-prep.pc, prep.pc)) )							
		}		
		if(!any('PREP_EFF'==names(df.rec)))
			df.rec[, PREP_EFF:=NA_real_]
		if(grepl('eff=median', method.sample))
		{
			if(!any('PREP_EFF_PAR'==names(df.rec)) & prep.eff=='IPrEX' )
				df.rec[, PREP_EFF_PAR:=0.44]
			if(!any('PREP_EFF_PAR'==names(df.rec)) & prep.eff=='PROUD' )
				df.rec[, PREP_EFF_PAR:=0.86]			
			tmp		<- df.rec[, which(is.na(PREP_EFF))]				
			z		<- c( rep(0, ceiling(length(tmp)*(1-df.rec[1, PREP_EFF_PAR]))), rep(1, ceiling(length(tmp)*df.rec[1, PREP_EFF_PAR])) )
			set(df.rec, tmp, 'PREP_EFF', z[seq_along(tmp)])			
		}
		if(grepl('eff=sample', method.sample))
		{			
			if(!any('PREP_EFF_PAR'==names(df.rec)) & prep.eff=='IPrEX' )
				df.rec[, PREP_EFF_PAR:=rbeta(1, 6, (1-0.44)/0.44*6 )]
			if(!any('PREP_EFF_PAR'==names(df.rec)) & prep.eff=='PROUD' )
				df.rec[, PREP_EFF_PAR:=rbeta(1, 2.9, (1-0.86)/0.86*2.9 )]			
			tmp		<- df.rec[, which(is.na(PREP_EFF))]
			set(df.rec, tmp, 'PREP_EFF', sample(c(0,1), length(tmp), prob=c(1-df.rec[1, PREP_EFF_PAR], df.rec[1, PREP_EFF_PAR]), replace=TRUE) )
		}
		df.rec
	}	
	rTest.Time<- function(df.tr, method.sample, test.pc=0.5, prest.repeat=1, th.starttime=2008.5, th.endtime=2013.5)
	{		
		ans		<- unique(subset( df.tr, select=c(t.Patient, t.INFECTION_T) ))
		#	first test for every t.Patient after th.starttime
		if(grepl('t=mean',method.sample))
			ans[, REALLOC_T1:=rep(th.starttime+prest.repeat/2, nrow(ans))]
		if(grepl('t=start',method.sample))
			ans[, REALLOC_T1:=rep(th.starttime, nrow(ans))]		
		if(grepl('t=sample',method.sample))
			ans[, REALLOC_T1:=runif(nrow(ans), th.starttime-prest.repeat/2, th.starttime+prest.repeat/2)]
		#	subsequent tests 
		ans		<- ans[, 	{
					tmp	<- c( seq(REALLOC_T1, th.endtime-0.01, prest.repeat), Inf )
					list(REALLOC_T1=REALLOC_T1, REALLOC_T=tmp[ which(tmp>=t.INFECTION_T) ][1])					
				}, by='t.Patient']
		ans		<- merge(ans, df.tr, by='t.Patient')
		#	only fraction tests
		ans[, REALLOC_Tc:= floor(REALLOC_T)]
		setkey(ans, t.Patient)
		if(grepl('stage=sample',method.sample))
			tmp	<- unique(ans)[, list(t.Patient=t.Patient, PRESTED= sample(c(0,1), length(t.Patient), replace=TRUE, prob=c(1-test.pc, test.pc))), by='REALLOC_Tc']
		if(grepl('stage=prop',method.sample))
			tmp	<- unique(ans)[, {
						z		<- c(rep(0,ceiling(length(t.Patient)*(1-test.pc))), rep(1,ceiling(length(t.Patient)*test.pc)))
						if(length(z)%%2)
							z	<- c(z, NA)
						z		<- na.omit(as.vector(t(matrix(z, ncol=2))))
						list(t.Patient=t.Patient, PRESTED= z[seq_along(t.Patient)])		
					}, by='REALLOC_Tc']
		ans		<- merge(ans, subset(tmp, select=c(t.Patient, PRESTED)), by='t.Patient')
		ans
	}
	resolve.UAna<- function(df.all, df.select, method.sample, p.UAcond, verbose=TRUE)
	{
		setkey(df.all, Patient)
		df.all	<- unique(df.all)
		tmp3	<- df.all[, which(Patient%in%df.select[, t.Patient] & is.na(isAcute))]
		if(verbose)
			cat(paste('\nFound NA isAcute among transmitters, allocate to Yes or No, n=',length(tmp3)))
		if(grepl('stage=sample',method.sample))
			set(df.all, tmp3, 'isAcute', factor(rbinom(length(tmp3), 1, p.UAcond), levels=c(0,1), labels=c('No','Yes')) )
		if(grepl('stage=prop',method.sample))
		{
			tmp		<- rep(c('No','Yes'), ceiling( c(1-p.UAcond, p.UAcond)*length(tmp3)))
			tmp		<- tmp[seq_along(tmp3)]
			set(df.all, tmp3, 'isAcute', tmp )
		}		
		df.all
	}
	
	#	get proportion of stages to reallocate to
	tp			<- as.numeric(substr(regmatches(method.risk,regexpr('tp[0-9]', method.risk)),3,3))
	tmp			<- subset(YX, grepl('Dt|DA',stage))[, table(stage)]
	df.dinfo	<- subset( data.table(stage=names(tmp), nt=tmp, pt=tmp/sum(tmp)), nt>0 )	
	#	use YXf get set of new scores to sample from
	if(is.null(YXf))
		tmp		<- subset(YX, select=c(stage, Patient, t.Patient, score.Y.raw))
	if(!is.null(YXf))
		tmp		<- subset(YXf, select=c(stage, Patient, t.Patient, score.Y.raw))
	setkey(tmp, Patient, t.Patient, stage)
	tmp			<- unique(tmp)
	set( tmp, NULL, 'stage', tmp[, as.character(stage)] )		#use all Diag intervals, not just those from tp4
	set( tmp, NULL, 'stage', tmp[, paste(substr(stage,1,nchar(stage)-1),tp,sep='')] )
	df.dinfo	<- merge(df.dinfo, tmp, by='stage')
	set(df.dinfo, NULL, 'stage', df.dinfo[, as.character(stage)])	
	#	estimated time of infection of transmitters
	if(verbose)
		cat(paste('\nusing method.resolveUAna=',method.resolveUAna))
	#	get time of infection of transmitters
	#	method.resolveUAna<- 0	
	if(method.resolveUAna)
	{
		tmp		<- subset(YX, grepl(method.reallocate, stage), select=c(t.Patient, t, stage))
		p.UAcond<- subset( tmp, grepl('UA.',stage, fixed=TRUE) | grepl('U.',stage, fixed=TRUE) )[, sum(grepl('UA.',stage, fixed=TRUE))/length(stage)]
		#	of those with isAcute NA, set a proportion p.UAcond to Yes and the rest to No 
		tmp		<- unique(subset(tmp, select=t.Patient))
		tmp2	<- resolve.UAna(df.all, tmp, method.sample, p.UAcond, verbose=verbose)
		df.tr	<- project.athena.Fisheretal.Y.infectiontime(tmp, tmp2, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, method='for.transmitter', method.minLowerUWithNegT=method.minLowerUWithNegT, verbose=FALSE)
	}
	if(!method.resolveUAna)
	{
		tmp			<- unique(subset(YX, grepl(method.reallocate, stage), select=t.Patient))
		df.tr		<- project.athena.Fisheretal.Y.infectiontime(tmp, df.all, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, method='for.transmitter', method.minLowerUWithNegT=method.minLowerUWithNegT, verbose=FALSE)		
	}	
	setkey(df.tr, t.Patient, t)
	df.tr		<- df.tr[, list(t=t[1]), by='t.Patient']
	setnames(df.tr, 't','t.INFECTION_T')
	#	sample test times of transmitters	
	df.tr		<- rTest.Time(df.tr, method.sample, test.pc=test.pc, prest.repeat=prest.repeat, th.starttime=th.starttime, th.endtime=2013.5)
	#
	#	those tested are offered PrEP, which may be efficacious from the first test time onwards
	#	so PrEP improves on testing
	#
	tmp			<- subset(df.tr, select=c(t.Patient, PRESTED))
	setnames(tmp, 't.Patient', 'Patient')
	tmp			<- rPrEP(tmp, method.sample, prep.pc, prep.eff=prep.eff)		#prep.pc does not do anything here because all transmitters have !NA PRESTED
	setnames(tmp, 'Patient', 't.Patient')
	df.tr		<- merge(df.tr, unique(subset(tmp, select=c(t.Patient, PREP_EFF, PREP_EFF_PAR))), by='t.Patient')
	#	use REALLOC_T1 to determine if first test before infection time. REALLOC_T is always first test AFTER infection time 
	df.tr[, TR_PREP_ON:= as.numeric(PRESTED & PREP_EFF & t.INFECTION_T+prest.delay>REALLOC_T1)]	
	if(verbose)
	{
		cat(paste('\nFound transmitters that are tested at baseline and at repeat points, n=', df.tr[, round(mean(PRESTED==1),d=3) ]))
		cat(paste('\nFound transmitters with hyp effective PrEP, n=', df.tr[, round(mean(PREP_EFF==1),d=3) ]))
		cat(paste('\nFound transmitters that could be on PrEP if inf time is after neg test, n=', df.tr[, round(mean(PRESTED & PREP_EFF),d=3) ]))
		cat(paste('\nFound transmitters that are on PrEP, n=', df.tr[, round(mean(TR_PREP_ON==1),d=3) ]))
	}
	#
	#	some of the transmitters are also recipients: offer PrEP consistently	
	#
#tmp	<- sapply(1:1000, function(j){
	df.rec		<- unique(subset(nt.table, select=Patient))
	tmp			<- unique(subset(df.tr, select=c(t.Patient, PRESTED, PREP_EFF, PREP_EFF_PAR)))
	setnames(tmp, 't.Patient','Patient')
	df.rec		<- merge(df.rec, tmp, by='Patient', all.x=1)
	if(df.rec[,!all(is.na(PREP_EFF_PAR))])
		set(df.rec, df.rec[, which(is.na(PREP_EFF_PAR))], 'PREP_EFF_PAR', df.rec[which(!is.na(PREP_EFF_PAR))[1],PREP_EFF_PAR])
	if(df.rec[,all(is.na(PREP_EFF_PAR))])
		set(df.rec, NULL, 'PREP_EFF_PAR', NULL)	
	#	for those that are no transmitter, allocate if Prested and if PREP efficacious
	df.rec		<- rPrEP(df.rec, method.sample, prep.pc, prep.eff=prep.eff)
	df.rec[, REC_PREP_ON:= as.numeric(PRESTED & PREP_EFF)]
	if(verbose)
	{
		cat(paste('\nFound recipients that are offered PrEP at baseline, n=', df.rec[, round(length(which(PRESTED==1))/length(PRESTED),d=3) ]))
		cat(paste('\nFound recipients with hyp effective PrEP, n=', df.rec[, round(length(which(PREP_EFF==1))/length(PRESTED),d=3) ]))
		cat(paste('\nFound recipients that are on PrEP, n=', df.rec[, round(length(which(REC_PREP_ON==1))/length(PRESTED),d=3) ]))
	}
	#df.rec[, length(which(REC_PREP_ON=='1'))/length(REC_PREP_ON)]
#})#	summary(tmp)
	#
	#	combine test of transmitters, PrEP of transmitters and PrEP of recipients  
	#
	
	#	
	#	determine which intervals could have been avoided trough testing
	YX.h			<- copy(YX)
	YX.h			<- merge(YX.h, subset(df.tr, select=c(t.Patient, t.INFECTION_T, REALLOC_T1, REALLOC_T, PRESTED, TR_PREP_ON)), by='t.Patient', all.x=TRUE)
	YX.h			<- merge(YX.h, subset(df.rec, select=c(Patient, REC_PREP_ON)), by='Patient', all.x=TRUE)
	if(!method.resolveUAna)
		stopifnot(nrow(subset(YX.h, t.INFECTION_T>t))==0)	
	if(verbose)
		cat(paste('\ntransmission intervals, n=', nrow(YX.h)))	
	#	remove intervals to recipients that are on effective PrEP
	YX.h			<- subset(YX.h, REC_PREP_ON==0)
	if(verbose)
		cat(paste('\ntransmission intervals to recipients not on eff PrEP, n=', nrow(YX.h)))
	#	proportion of transmission intervals for which PrEP works, after we rm recipient on PrEP
	#	need this to reduce nt.table
	df.uinfo		<- merge( subset(YX.h, TR_PREP_ON==1)[, list(TR_RM= length(t)), by='stage'], YX.h[, list(N= length(t)), by='stage'], by='stage', all.y=1)
	set(df.uinfo, df.uinfo[, which(is.na(TR_RM))], 'TR_RM', 0)
	df.uinfo[, TR_NRM:= 1-TR_RM/N]
	set(df.uinfo, NULL, c('TR_RM','N'), NULL)
	#	remove intervals from transmitters that are on effective PrEP
	YX.h			<- subset(YX.h, is.na(TR_PREP_ON) | TR_PREP_ON==0)
	if(verbose)
		cat(paste('\ntransmission intervals from transmitters not on eff PrEP, n=', nrow(YX.h)))
	#	realloc intervals from transmitters that tested positive
	YX.h		<- merge(YX.h, YX.h[, list(t_min=min(t)), by=c('Patient','t.Patient')], by=c('Patient','t.Patient'))
	set(YX.h, NULL, 't.grace', YX.h[,t.INFECTION_T-t_min])
	set(YX.h, YX.h[, which(!grepl('UAna',stage) | t.grace<0 | is.na(t.grace)) ], 't.grace', 0)
	if(verbose)
		cat(paste('\ntransmission intervals with t.grace>0, n=', YX.h[, length(which(t.grace>0))]))	
	#	re-alloc intervals of tested transmitters to diag
	#	t.INFECTION_T+prest.delay<=REALLOC_T	--> infected before pos test
	#	t+t.grace>REALLOC_T						--> intervals after pos test
	YX.h[, TR_TEST:= 0]
	tmp			<- YX.h[, which( 	!is.na(REALLOC_T) & PRESTED & t.INFECTION_T+prest.delay<=REALLOC_T & 
							t+t.grace>REALLOC_T & grepl(method.reallocate, stage))]
	set(YX.h, tmp, 'TR_TEST', 1)		
	if(verbose)
		cat(paste('\nrealloc transmission intervals from transmitters that got tested pos, n=', length(tmp)))
	tmp			<- unique(subset(YX.h, TR_TEST==1, select=c(Patient, t.Patient)))
	tmp			<- rTest.Realloc(tmp, method.sample, df.dinfo)
	YX.h		<- merge( YX.h, tmp, by=c('Patient','t.Patient'), all.x=1 )
	set(YX.h, YX.h[, which(TR_TEST==0)], 'REALLOC_STAGE', NA_character_)
	set(YX.h, YX.h[, which(TR_TEST==0)], 'REALLOC_SCORE_Y_RAW', NA_real_)
	#
	#	add proportion of in and out reallocated stages to df.uinfo
	#
	tmp			<- subset(YX.h, !is.na(REALLOC_STAGE))[, list(REALLOC_NT_IN= length(t)), by='REALLOC_STAGE']	
	setnames(tmp, 'REALLOC_STAGE', 'stage')
	if(!nrow(tmp))
		set(df.uinfo,NULL, c('U_OUT_AFTER_RM','D_IN_AFTER_RM'),0)
	if(nrow(tmp))
	{
		df.uinfo	<- merge(df.uinfo, tmp, by='stage', all.x=1)
		df.uinfo	<- merge( df.uinfo, YX.h[, list(NT= length(t)), by='stage'], by='stage', all=1 )
		df.uinfo	<- merge( df.uinfo, subset(YX.h, !is.na(REALLOC_STAGE))[, list(REALLOC_NT_OUT= length(t)), by='stage'], by='stage', all=1 )
		set(df.uinfo, df.uinfo[, which(is.na(REALLOC_NT_IN))], 'REALLOC_NT_IN', 0)
		set(df.uinfo, df.uinfo[, which(is.na(REALLOC_NT_OUT))], 'REALLOC_NT_OUT', 0)
		df.uinfo[, U_OUT_AFTER_RM:= REALLOC_NT_OUT/NT]
		set(df.uinfo, df.uinfo[, which(is.na(NT))], c('NT','U_OUT_AFTER_RM'), 0)
		df.uinfo[, D_IN_AFTER_RM:= REALLOC_NT_IN/sum(REALLOC_NT_IN)]		
		set(df.uinfo, NULL, c('REALLOC_NT_IN','REALLOC_NT_OUT','NT'), NULL)
	}
	setnames(df.uinfo, 'stage', 'factor')
	set(df.uinfo, NULL, 'factor', df.uinfo[, as.character(factor)])
	df.uinfo[, risk:='stage']		
	#
	#	reallocate
	#
	tmp			<- YX.h[, which(TR_TEST==1)]
	set(YX.h, tmp, 'stage', YX.h[tmp, REALLOC_STAGE])
	set(YX.h, tmp, c('score.Y','score.Y.raw'), YX.h[tmp, REALLOC_SCORE_Y_RAW])
	set(YX.h, NULL, c('t_min', 't.grace', 't.INFECTION_T', 'REALLOC_T1','REALLOC_T', 'REALLOC_STAGE', 'REALLOC_SCORE_Y_RAW', 'PRESTED', 'TR_PREP_ON', 'REC_PREP_ON', 'TR_TEST'), NULL)
	#	update scores
	if(grepl('wtn',method.risk))
	{
		if(verbose)
			cat(paste('\nsetting likelihood to likelihood of pair / number transmission intervals'))
		set(YX.h, NULL, 'score.Y', YX.h[, score.Y.raw*w.tn])				
	}									
	if(grepl('noscore',method.risk))
	{
		if(verbose)
			cat(paste('\nsetting likelihood to 1'))
		set(YX.h, NULL, 'score.Y', YX.h[, 1])						
	}
	if(grepl('nophyloscore',method.risk))
	{
		if(verbose)
			cat(paste('\nsetting likelihood to nophyloscore'))
		set(YX.h, NULL, 'score.Y', YX.h[, w.tn])						
	}
	#
	#	prepare nt.table for hypothetical scenario
	#			
	#	for Xclu Xseq Xmsm, we don t have NegT etc unless we load the bigmem table X.msm 
	#	avoid this by removing a fraction
	setkey(df.uinfo, risk, factor)
	setkey(nt.table, Patient, stat, risk, factor)
	#	remove recipients that are on PrEP	
	#	Xmsm includes participants with no transmitter, so need to use df.rec
	nt.table.h		<- merge( nt.table, subset(df.rec, !REC_PREP_ON, Patient), by='Patient' )	
	#	
	#	remove transmitters that are on PrEP and thereafter, move stages according to fractions UPT and DPT
	#
	#	subset(nt.table, Patient=='M41498' & stat=='X.msm')
	if(df.uinfo[, all(U_OUT_AFTER_RM==0)] & grepl('stage=sample',method.sample))
		nt.table.h	<- subset(nt.table.h, stat!='YX')[, {
					nt.nrm		<- rbinom(length(nt), nt, df.uinfo$TR_NRM)
					list(factor=factor, nt= nt.nrm)
				}, by=c('Patient','risk','stat')]
	if(df.uinfo[, all(U_OUT_AFTER_RM==0)] & grepl('stage=prop',method.sample))
		nt.table.h	<- subset(nt.table.h, stat!='YX')[, {
					nt.nrm		<- round( nt * df.uinfo$TR_NRM)
					list(factor=factor, nt= nt.nrm)					
				}, by=c('Patient','risk','stat')]	
	if(df.uinfo[, !all(U_OUT_AFTER_RM==0)] & grepl('stage=sample',method.sample))
		nt.table.h	<- subset(nt.table.h, stat!='YX')[, {
					nt.nrm		<- rbinom(length(nt), nt, df.uinfo$TR_NRM)
					flow.out	<- rbinom(length(nt.nrm), nt.nrm, df.uinfo$U_OUT_AFTER_RM) 
					flow.in		<- as.numeric(rmultinom(1, sum(flow.out), df.uinfo$D_IN_AFTER_RM))
					stopifnot(sum(flow.out)==sum(flow.in))
					#print(nt); print(flow.out); print(flow.in)
					list(factor=factor, nt= nt.nrm-flow.out+flow.in)
				}, by=c('Patient','risk','stat')]
	if(df.uinfo[, !all(U_OUT_AFTER_RM==0)] & grepl('stage=prop',method.sample))
		nt.table.h	<- subset(nt.table.h, stat!='YX')[, {
					nt.nrm		<- round( nt * df.uinfo$TR_NRM)
					flow.out	<- round( nt.nrm * df.uinfo$U_OUT_AFTER_RM )
					flow.in		<- ceiling( sum(flow.out) * df.uinfo$D_IN_AFTER_RM )
					flow.d		<- sum(flow.in)-sum(flow.out)
					stopifnot( flow.d>=0, flow.d<=length(flow.in>0) )
					flow.in[ flow.in>0 ]	<- flow.in[ flow.in>0 ] + c( rep(-1, flow.d), rep(0, length(which(flow.in>0))-flow.d) )
					stopifnot(sum(flow.out)==sum(flow.in))
					list(factor=factor, nt= nt.nrm-flow.out+flow.in)					
				}, by=c('Patient','risk','stat')]		
	#	
	#	for YX, we know easily how many transmission intervals are removed / reallocated
	nt.table.YX	<- YX.h[, 	{
				z				<- table(as.character(stage))
				list(risk='stage', stat='YX', factor=names(z), nt=as.numeric(z))
			}, by='Patient']
	tmp			<- merge( unique(subset(nt.table.h, select=c(risk, Patient))), unique(subset(nt.table.h, select=c(risk, factor))), by='risk', allow.cartesian=TRUE)
	nt.table.YX	<- merge( nt.table.YX, tmp, by=c('Patient','risk','factor'), all.y=TRUE) 		
	set(nt.table.YX, nt.table.YX[, which(is.na(nt))], 'stat', 'YX')
	set(nt.table.YX, nt.table.YX[, which(is.na(nt))], 'nt', 0)
	#
	nt.table.h	<- rbind(nt.table.YX, nt.table.h, use.names=TRUE)
	stopifnot( nrow(merge(nt.table.h, tmp, all.y=1, by=c('Patient','risk','factor')))==nrow(tmp)*4 )
	nt.table.h	<- dcast.data.table(nt.table.h, Patient+risk+factor~stat, value.var='nt')
	tmp			<- nt.table.h[, which(YX>X.clu)]
	if(length(tmp))
	{
		cat(paste('\nWARNING UtoD Found YX>X.clu, n=',length(tmp)))
		set(nt.table.h, tmp, 'X.clu', nt.table.h[tmp, YX])
	}
	tmp			<- nt.table.h[, which(X.clu>X.seq)]	
	if(length(tmp))
	{
		cat(paste('\nWARNING UtoD Found X.clu>X.seq, n=',length(tmp)))	
		set(nt.table.h, tmp, 'X.seq', nt.table.h[tmp, X.clu])
	}
	tmp			<- nt.table.h[, which(X.seq>X.msm)]	
	if(length(tmp))
	{
		cat(paste('\nWARNING UtoD Found X.seq>X.msm, n=',length(tmp)))
		set(nt.table.h, tmp, 'X.msm', nt.table.h[tmp, X.seq])
	}	
	nt.table.h	<- melt(nt.table.h, measure.vars=c('X.clu','X.msm','X.seq','YX'), value.name='nt', variable.name='stat')
	#	
	list(YX.h=YX.h, nt.table.h=nt.table.h)	 
}
######################################################################################
project.athena.Fisheretal.Hypo.ReallocPrEP.getYXetc<- function( YX, nt.table, method.risk, predict.t2inf, t2inf.args, df.all, th.starttime=2008.5, t.period=0.125, method.minLowerUWithNegT=1, method.realloc='PrEP33', method.sample= 'stage=prop, y=mean', verbose=FALSE)
{
	stopifnot(grepl('stage=sample|stage=prop',method.sample))
	stopifnot(grepl('eff=sample|eff=median',method.sample))
	p.reachable			<- as.numeric(substring( regmatches(method.realloc,regexpr('RPrEP[0-9]+', method.realloc)), 6))/100
	if(verbose)
		cat(paste('\nsetting p.reachable=', p.reachable))
	
	rEfficacy<- function(method.sample)
	{
		#tmp		<- data.table(x=seq(0,1,0.001))
		#tmp[, a:= 6]
		#tmp[, b:= (1-0.44)/0.44*a]
		#tmp[,y:=dbeta(x,a,b)]
		#ggplot(tmp, aes(x=x, y=y)) + geom_line()
		#quantile( rbeta(1e4, tmp[1,a],tmp[1,b] ), probs=c(0.025,0.975) )		
		if(grepl('eff=median',method.sample))
			ans	<- rep(0.44, 1)
		if(grepl('eff=sample',method.sample))
			ans	<- rbeta(1, 6, (1-0.44)/0.44*6 )
		ans
	}				
	#tp				<- as.numeric(substr(regmatches(method.risk,regexpr('tp[0-9]', method.risk)),3,3))
	#subset(unique(subset(merge(data.table(Patient=YX[, unique(t.Patient)]), df.all, by='Patient'), select=c(Patient, AnyPos_T1, NegT))), AnyPos_T1-NegT<3)
	#subset( df.all, AnyPos_T1>2009.5, select=c(Patient, AnyPos_T1, NegT) )[, table(!is.na(NegT) & AnyPos_T1-NegT<1)/length(NegT)]
	tmp				<- unique(subset(YX, t.AnyPos_T1>=th.starttime, select=t.Patient))
	df.tr			<- project.athena.Fisheretal.Y.infectiontime(tmp, df.all, predict.t2inf, t2inf.args, t.period=t.period, ts.min=1980, score.set.value=NA, method='for.transmitter', method.minLowerUWithNegT=method.minLowerUWithNegT, verbose=FALSE)
	df.tr			<- df.tr[, list(t=tail(t,1)), by='t.Patient']
	#
	#	prepare removable transmitters
	#			
	#	select transmitters for which RPrEP in TP4 might have worked
	df.tr			<- subset(df.tr, t>=th.starttime, select=t.Patient)	
	tmp				<- merge( subset(YX, select=c(t.Patient, stage)), subset(df.tr, select=t.Patient), by='t.Patient' )
	df.uinfo		<- data.table(stage= tmp[, names(table(as.character(stage)))], nt=tmp[, table(as.character(stage))])
	df.uinfo		<- subset(df.uinfo, nt>5)
	df.uinfo[, pt:= nt/sum(nt)]	
	#	select those that could have been reached before infected: take actual NegT within past year
	tmp				<- unique(subset(df.all, select=c(Patient, NegT)))
	setnames(tmp, colnames(tmp), paste('t.',colnames(tmp),sep=''))
	df.tr			<- merge( df.tr, tmp, by='t.Patient' )
	#	offer PrEP to those with a neg test since  th.starttime
	df.tr			<- subset(df.tr, !is.na(t.NegT) & th.starttime<t.NegT)	
	#	select proportion for whom PrEP works	
	df.tr[, pEfficacy:=rEfficacy(method.sample)]
	if(grepl('stage=prop',method.sample))	#select first 
	{
		set(df.tr, NULL, 'rEfficacy', 1)
		set(df.tr, seq_len(round(nrow(df.tr)*df.tr[1, pEfficacy])), 'rEfficacy', 0)		
	}
	if(grepl('stage=sample',method.sample))
		df.tr[, rEfficacy:=runif(nrow(df.tr))]
	df.tr			<- subset(df.tr, rEfficacy<=pEfficacy)
	#
	#	prepare removable recipients
	#			
	df.rec		<- unique(subset(nt.table, select=Patient))
	df.rec[, pEfficacy:=rEfficacy(method.sample)]	
	df.rec[, pReachable:=p.reachable]	
	if(grepl('stage=prop',method.sample))	#select first 
	{
		set(df.rec, NULL, c('rReachable','rEfficacy'), 0)		
		set(df.rec, seq_len(round(nrow(df.rec)*df.rec[1, pEfficacy*pReachable])), c('rReachable','rEfficacy'), 1)
		df.rec[, REC_PREP_ON:= as.numeric(rReachable & rEfficacy)]
	}
	if(grepl('stage=sample',method.sample))
	{
		df.rec[, rEfficacy:=runif(nrow(df.rec))]
		df.rec[, rReachable:=runif(nrow(df.rec))]
		df.rec[, REC_PREP_ON:= as.numeric(rReachable<=pReachable & rEfficacy<=pEfficacy)]
	}		
	#	
	#	get YX.hypothetical
	#
	if(verbose)
	{
		tmp			<- merge(YX, df.tr, by='t.Patient')[, paste(unique(as.character(stage)), collapse=' ')]
		cat(paste('\nrm stages because transmitters lost=',tmp))		
		tmp			<- merge(YX, df.rec, by='Patient')[, paste(unique(as.character(stage)), collapse=' ')]
		cat(paste('\nrm stages because recipients lost=',tmp))
	}	
	YX.h			<- copy(YX)
	if(verbose)
		cat(paste('\ntransmission intervals, n=', nrow(YX.h)))	
	YX.h			<- merge( YX.h, subset( df.rec,  REC_PREP_ON==0, Patient ), by='Patient')
	if(verbose)
		cat(paste('\ntransmission intervals to recipients not on eff PrEP, n=', nrow(YX.h)))	
	tmp				<- unique( subset( YX.h, select=t.Patient ) )
	df.tr[, TR_PREP_ON:= 1]
	tmp				<- merge(tmp, df.tr, by='t.Patient', all.x=TRUE)
	set(tmp, tmp[,which(is.na(TR_PREP_ON))],'TR_PREP_ON',0)	
	YX.h			<- merge(YX.h, subset(tmp, select=c(t.Patient, TR_PREP_ON)), by='t.Patient')
	#
	df.uinfo		<- merge( subset(YX.h, TR_PREP_ON==1)[, list(TR_RM= length(t)), by='stage'], YX.h[, list(N= length(t)), by='stage'], by='stage', all.y=1)
	set(df.uinfo, df.uinfo[, which(is.na(TR_RM))], 'TR_RM', 0)
	df.uinfo[, TR_NRM:= 1-TR_RM/N]
	set(df.uinfo, NULL, c('TR_RM','N'), NULL)
	setnames(df.uinfo, 'stage', 'factor')
	set(df.uinfo, NULL, 'factor', df.uinfo[, as.character(factor)])
	df.uinfo[, risk:='stage']	
	#	remove intervals from transmitters that are on effective PrEP
	YX.h			<- subset(YX.h, TR_PREP_ON==0)
	if(verbose)
		cat(paste('\ntransmission intervals from transmitters not on eff PrEP, n=', nrow(YX.h)))	
	#	prepare scores				 
	set(YX.h, NULL, 'score.Y', YX.h[, score.Y.raw])
	set(YX.h, NULL, 'score.Y.raw', NULL)
	set(YX.h, NULL, 'stage', YX.h[, factor(as.character(stage))])	
	if(grepl('wtn',method.risk))
	{				
		set(YX.h, NULL, 'score.Y', YX.h[, score.Y*w.tn])				
	}									
	if(grepl('noscore',method.risk))
	{
		if(verbose) cat(paste('\nsetting likelihood to 1'))
		set(YX.h, NULL, 'score.Y', YX.h[, 1])						
	}
	if(grepl('nophyloscore',method.risk))
	{
		if(verbose) cat(paste('\nsetting likelihood to nophyloscore'))
		set(YX.h, NULL, 'score.Y', YX.h[, w.tn])						
	}	
	if(verbose)
	{
		cat(paste('\nBefore PrEP, nRec=',YX[, length(unique(Patient))],' nTrans=', YX[, length(unique(t.Patient))], ' nInt=', nrow(YX)))
		cat(paste('\nAfter PrEP, nRec=',YX.h[, length(unique(Patient))],' nTrans=',YX.h[, length(unique(t.Patient))], ' nInt=', nrow(YX.h)))
	}
	#
	#	prepare nt.table for hypothetical scenario
	#			
	#	for Xclu Xseq Xmsm, we don t have NegT etc unless we load the bigmem table X.msm 
	#	avoid this by removing a fraction
	setkey(df.uinfo, risk, factor)
	setkey(nt.table, Patient, stat, risk, factor)
	#	remove recipients that are on PrEP	
	#	Xmsm includes participants with no transmitter, so need to use df.rec
	nt.table.h		<- merge( nt.table, subset(df.rec, REC_PREP_ON==0, Patient), by='Patient' )	
	#	
	#	remove transmitters that are on PrEP and thereafter, move stages according to fractions UPT and DPT
	#
	#	subset(nt.table, Patient=='M41498' & stat=='X.msm')
	if(grepl('stage=sample',method.sample))
		nt.table.h	<- subset(nt.table.h, stat!='YX')[, {
					nt.nrm		<- rbinom(length(nt), nt, df.uinfo$TR_NRM)
					list(factor=factor, nt= nt.nrm)
				}, by=c('Patient','risk','stat')]
	if(grepl('stage=prop',method.sample))
		nt.table.h	<- subset(nt.table.h, stat!='YX')[, {
					nt.nrm		<- round( nt * df.uinfo$TR_NRM)
					list(factor=factor, nt= nt.nrm)					
				}, by=c('Patient','risk','stat')]	
	#	
	#	for YX, we know easily how many transmission intervals are removed / reallocated
	nt.table.YX	<- YX.h[, 	{
				z				<- table(as.character(stage))
				list(risk='stage', stat='YX', factor=names(z), nt=as.numeric(z))
			}, by='Patient']
	tmp			<- merge( unique(subset(nt.table.h, select=c(risk, Patient))), unique(subset(nt.table.h, select=c(risk, factor))), by='risk', allow.cartesian=TRUE)
	nt.table.YX	<- merge( nt.table.YX, tmp, by=c('Patient','risk','factor'), all.y=TRUE) 		
	set(nt.table.YX, nt.table.YX[, which(is.na(nt))], 'stat', 'YX')
	set(nt.table.YX, nt.table.YX[, which(is.na(nt))], 'nt', 0)
	#
	nt.table.h	<- rbind(nt.table.YX, nt.table.h, use.names=TRUE)
	stopifnot( nrow(merge(nt.table.h, tmp, all.y=1, by=c('Patient','risk','factor')))==nrow(tmp)*4 )
	nt.table.h	<- dcast.data.table(nt.table.h, Patient+risk+factor~stat, value.var='nt')
	tmp			<- nt.table.h[, which(YX>X.clu)]
	if(length(tmp))
	{
		cat(paste('\nWARNING UtoD Found YX>X.clu, n=',length(tmp)))
		set(nt.table.h, tmp, 'X.clu', nt.table.h[tmp, YX])
	}
	tmp			<- nt.table.h[, which(X.clu>X.seq)]	
	if(length(tmp))
	{
		cat(paste('\nWARNING UtoD Found X.clu>X.seq, n=',length(tmp)))	
		set(nt.table.h, tmp, 'X.seq', nt.table.h[tmp, X.clu])
	}
	tmp			<- nt.table.h[, which(X.seq>X.msm)]	
	if(length(tmp))
	{
		cat(paste('\nWARNING UtoD Found X.seq>X.msm, n=',length(tmp)))
		set(nt.table.h, tmp, 'X.msm', nt.table.h[tmp, X.seq])
	}	
	nt.table.h	<- melt(nt.table.h, measure.vars=c('X.clu','X.msm','X.seq','YX'), value.name='nt', variable.name='stat')
	#	
	list(YX.h=YX.h, nt.table.h=nt.table.h)
}