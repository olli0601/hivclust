######################################################################################
project.athena.Fisheretal.revision.subtypeB<- function()
{
	tmp		<- subset(df.all, Trm%in%c('MSM','BI') & PosSeqT<'2013-04-01')
	setkey(tmp, Patient)
	tmp		<- unique(tmp)
	tmp[, table(Subtype)]		#94% are B
	
	tmp[, tC:= cut(hivc.db.Date2numeric(PosSeqT), breaks=c(-Inf,2003, 2006, 2009, Inf))]
	
	tmp[, table(Subtype, tC)]
	tmp[, mean(Subtype=='B'), by='tC']
}
######################################################################################
project.athena.Fisheretal.revision.coclustering.recentlate<- function()
{
	crl		<- clumsm.info[, {
				
				LATn	<- ADVn	<- 0L
				z		<- which(CD4_T1<200 & Trm%in%c('BI','MSM'))
				if(length(z))
					ADVn= length(setdiff(Patient[z], ri.CLU$Patient))
				z		<- which(CD4_T1<350 & Trm%in%c('BI','MSM'))
				if(length(z))
					LATn= length(setdiff(Patient[z], ri.CLU$Patient))				
				list( LATn=LATn, ADVn=ADVn, RECn=length(intersect(Patient, ri.CLU$Patient)), clu.npat=clu.npat[1], clu.nFrgnInfection=clu.nFrgnInfection[1], clu.AnyPos_T1=clu.AnyPos_T1[1] )
			}, by='cluster']
	crl[, LATp:= LATn/clu.npat]
	crl[, RECp:= RECn/clu.npat]
	crl[, RECc:= RECn>0]
	crl[, tC:= cut(clu.AnyPos_T1, breaks=c(1985, 2000, 2007, 2013))]
	
	crlm1	<- gamlss(formula=LATn~RECn, family=NO(), data=subset(crl, select=c(LATn, RECn, clu.npat, clu.nFrgnInfection, clu.AnyPos_T1)))
	crl[, m1p:=predict(crlm1)]	
	ggplot(crl, aes(x=RECn, colour=clu.AnyPos_T1)) +
			geom_abline(slope=1, intercept=0) +
			geom_line(aes(y=m1p)) +
			geom_jitter(aes(y=LATn), position = position_jitter(width=1/3, height=1/3), alpha=0.8) + theme_bw()
	
	crlm5	<- gamlss(formula=LATn~RECn:tC, family=NO(), data=subset(crl, select=c(LATn, RECn, tC, clu.npat, clu.nFrgnInfection, clu.AnyPos_T1)))
	crl[, m5p:=predict(crlm5)]
	ggplot(crl, aes(x=RECn, colour=clu.AnyPos_T1)) +
			geom_abline(slope=1, intercept=0) +
			geom_line(aes(y=m5p)) +			
			geom_jitter(aes(y=LATn), position = position_jitter(width=1/3, height=1/3), alpha=0.8) + theme_bw() +
			facet_grid(~tC)
	
	crlm6	<- gamlss(formula=LATn~RECc:tC, family=NO(), data=subset(crl, select=c(LATn, RECc, tC, clu.npat, clu.nFrgnInfection, clu.AnyPos_T1)))
	crl[, m6p:=predict(crlm6)]
	crl[, tC:= cut(clu.AnyPos_T1, breaks=c(1985, 2000, 2007, 2011, 2013))]
	crl[, clu.npatC:= cut(clu.npat, breaks=c(-1,3,7,60), labels=c('<4','4-7','>7'))]
	box.stat <- function(x)
	{
		v <- c(max(min(x), mean(x)-diff(quantile(x, p=c(0.25, 0.75)))), mean(x) - 2*sd(x)/sqrt(length(x)), mean(x), mean(x) + 2*sd(x)/sqrt(length(x)), min(max(x),mean(x)+diff(quantile(x, p=c(0.25, 0.75)))))
		names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
		v
	}
	tmp		<- melt(crl, measure.vars=c('LATn','ADVn'))
	set(tmp, tmp[, which(variable=='LATn')],'variable','late presenters\nin phylogenetic cluster')
	set(tmp, tmp[, which(variable=='ADVn')],'variable','men with advanced disease\nin phylogenetic cluster')
	ggp		<- ggplot(tmp, aes(x=RECc, y=value)) +
			stat_summary(fun.data=box.stat, geom="boxplot", outlier.shape =NA) +
			geom_jitter(aes(y=value), colour='grey70', position = position_jitter(width=1/3, height=1/10), alpha=0.7) + 							 
			theme_bw() +
			facet_grid(variable~tC, scales='free') + labs(x='\npresence of recipient MSM\nin cluster', y='#') +
			theme(legend.position='bottom')
	z <- ggplot_gtable(ggplot_build(ggp))
	z <- gtable_add_rows(z, z$heights[[3]], 2)
	z <- gtable_add_grob(z, 
			list(	rectGrob(height = unit(2, "npc"), gp = gpar(col='black', fill=gray(0.8))),
					textGrob("phylogenetic clusters by first date of diagnosis", gp = gpar(col = 'black'))), 3, 4, 3, 8, name = paste(runif(2)))
	z <- gtable_add_rows(z, unit(2/8, "line"), 3)
	grid.newpage()
	pdf(file='~/Dropbox (SPH Imperial College)/2014_MSMtransmission_ATHENA1303/150626_revisionlate.pdf', w=10, h=7)
	grid.draw(z)
	dev.off()
	
	ggp		<- ggplot(tmp, aes(x=RECc, y=value)) +
			stat_summary(fun.data=box.stat, geom="boxplot", outlier.shape =NA) +
			geom_jitter(aes(y=value), colour='grey70', position = position_jitter(width=1/3, height=1/10), alpha=0.7) + 							 
			theme_bw() +
			facet_grid(variable~clu.npatC, scales='free') + labs(x='\npresence of recipient MSM\nin cluster', y='#') +
			theme(legend.position='bottom')
	z <- ggplot_gtable(ggplot_build(ggp))
	z <- gtable_add_rows(z, z$heights[[3]], 2)
	z <- gtable_add_grob(z, 
			list(	rectGrob(height = unit(2, "npc"), gp = gpar(col='black', fill=gray(0.8))),
					textGrob("number of men in phylogenetic cluster", gp = gpar(col = 'black'))), 3, 4, 3, 8, name = paste(runif(2)))
	z <- gtable_add_rows(z, unit(2/8, "line"), 3)
	grid.newpage()
	pdf(file='~/Dropbox (SPH Imperial College)/2014_MSMtransmission_ATHENA1303/150626_revisionlate_clusize.pdf', w=10, h=7)
	grid.draw(z)
	dev.off()
	
	ggplot(crl, aes(x=clu.npat, y=LATn)) + geom_jitter(aes(colour=RECc), position = position_jitter(width=1/5, height=1/10), alpha=0.7) +			 
			geom_smooth(method='lm', colour='black', fill='black', alpha=0.4) +
			scale_colour_brewer(palette='Set2') +
			coord_trans(xtrans="log10")  +			
			geom_smooth(colour='blue', fill='blue', alpha=0.2) +
			theme_bw() +
			labs(colour='recipient MSM\nin cluster', x='men in cluster',y='late presenters in cluster') +
			theme(legend.position='bottom')
	ggsave(file='~/Dropbox (SPH Imperial College)/2014_MSMtransmission_ATHENA1303/150626_revisionlate_LATnvsSize.pdf', h=10, w=6)
	
	
	crlm2	<- gamlss(formula=LATn~RECc+clu.npat, family=NO(), data=subset(crl, select=c(LATn, RECc, RECn, clu.npat, clu.nFrgnInfection, clu.AnyPos_T1)))
	#	RECcTRUE      -0.4739    0.081008   -5.850   7.742e-09
	crlm2b	<- gamlss(formula=LATn~RECc+clu.npat, sigma.formula=~clu.npat, family=NO(), data=subset(crl, select=c(LATn, RECc, RECn, clu.npat, clu.nFrgnInfection, clu.AnyPos_T1)))
	#	RECcTRUE      -0.5747     0.06584   -8.729  2.119e-17
	crlm3	<- gamlss(formula=LATn~RECc+clu.npat, family=PO(), data=subset(crl, select=c(LATn, RECc, RECn, clu.npat, clu.nFrgnInfection, clu.AnyPos_T1)))
	#	RECcTRUE     -0.14073    0.078411   -1.795   7.269e-02
	crlm4	<- gamlss(formula=LATn~RECc+clu.npat, family=PO(mu.link='identity'), data=subset(crl, select=c(LATn, RECc, RECn, clu.npat, clu.nFrgnInfection, clu.AnyPos_T1)))
	#	RECcTRUE      -0.5525     0.07182   -7.693  1.433e-14
	crlp4	<- data.table(RECc=TRUE, clu.npat=crl[, sort(unique(clu.npat))])
	crlp4[, LATnp:= predict(crlm4, data=subset(crl, select=c(LATn, RECc, RECn, clu.npat, clu.nFrgnInfection, clu.AnyPos_T1)), newdata=crlp4)]
	crlp4	<- merge(crl, subset(crlp4, select=c(clu.npat, LATnp)), by='clu.npat')
	ggplot(crlp4, aes(x=RECc, y=LATn-LATnp)) +
			stat_summary(fun.data=box.stat, geom="boxplot", outlier.shape =NA) +
			geom_jitter(colour='grey70', position = position_jitter(width=1/3, height=1/10), alpha=0.7) +
			geom_hline(yintercept=0, colour='blue') +
			theme_bw() +
			#facet_grid(~tC, scales='free') + 
			labs(x='\npresence of recipient MSM\nin cluster', y='late presenters in cluster minus\nexpected number of late presenters in clusters with same size that have a recipient') +
			theme(legend.position='bottom')
	ggsave(file='~/Dropbox (SPH Imperial College)/2014_MSMtransmission_ATHENA1303/150626_revisionlate_contrastLATn.pdf', h=10, w=6)
	
	
	library(coin)
	wilcox_test(LATn~factor(RECc), data=subset(crl, tC=='(2007,2013]'), distribution='exact', alternative = "greater")
	#	factors: 1=FALSE, 2=TRUE. contrast is mu=Y_1-Y_2. so Null is mu<=0 (absence does not results in more late representers)
	oneway_test(LATn~factor(RECc), data=subset(crl, tC=='(2007,2013]'), distribution='exact', alternative = "greater")
	#	p-value = 0.03548 n=175
	oneway_test(ADVn~factor(RECc), data=subset(crl, tC=='(2007,2013]'), distribution='exact', alternative = "greater")
	#	p-value = 0.2035
	oneway_test(LATn~factor(RECc), data=crl, distribution='exact', alternative = "greater")
	#	p-value = 0.9996 n=659
	oneway_test(ADVn~factor(RECc), data=crl, distribution='exact', alternative = "greater")
	#	p-value = 0.8958
	crl[, tC:= cut(clu.AnyPos_T1, breaks=c(1985, 2006.6, 2008, 2009.5, 2011, 2013))]
	ggplot(crl, aes(x=RECc, y=LATn)) +
		geom_boxplot(outlier.shape =NA) +
		geom_jitter(aes(y=LATn), colour='grey70', position = position_jitter(width=1/3, height=1/3), alpha=0.7) + 				
		stat_summary(fun.y=mean, colour="red", geom="point", shape=18, size=3,show_guide = FALSE) + 
		theme_bw() +
		facet_grid(~tC) + labs(x='\npresence of recipient MSM\nin cluster', y='late presenters in phylogenetic cluster\n(#)') +
		theme(legend.position='bottom')
	
		
	
	crlm2	<- gamlss(formula=LATn~RECn+clu.npat-1, family=NO(), data=subset(crl, select=c(LATn, RECn, clu.npat, clu.nFrgnInfection, clu.AnyPos_T1)))
	crl[, m1p:=predict(crlm1)]
	
	crlm3	<- gamlss(formula=LATp~RECp-1, family=NO(), data=subset(crl, select=c(LATp, RECp, clu.npat, clu.nFrgnInfection, clu.AnyPos_T1)))
	crl[, m3p:=predict(crlm3)]
	ggplot(crl, aes(x=RECp, colour=clu.AnyPos_T1)) +
			geom_abline(slope=1, intercept=0) +
			geom_line(aes(y=m3p)) +
			geom_jitter(aes(y=LATp), position = position_jitter(width=1/100, height=1/100), alpha=0.8) + theme_bw()
	
	crlm4	<- gamlss(formula=LATp~RECp:tC, family=NO(), data=subset(crl, select=c(LATp, RECp, tC, clu.npat, clu.nFrgnInfection, clu.AnyPos_T1)))
	crl[, m4p:=predict(crlm4)]
	ggplot(crl, aes(x=RECp, colour=clu.AnyPos_T1)) +
			geom_line(aes(y=m4p)) +			
			geom_jitter(aes(y=LATp), position = position_jitter(width=1/100, height=1/100), alpha=0.8) + theme_bw() +
			facet_grid(~tC)
	
	
}
######################################################################################
project.athena.Fisheretal.composition.censoringfraction<- function()
{
	require(data.table)
	require(ape)	
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal",sep='/')			
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	outfile					<- infile
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period	
	#	read files
	files					<- list.files(indir, pattern='*tables.*R$')
	select					<- 'Yscore3pa1H1.48C2V100bInfT7_'
	files					<- files[  grepl(select, files) & grepl('m2A', files) ]	
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
	#
	cens.AnyPos_T1<- runs.opt[,	{
				tmp	<- paste(indir,'/',file,sep='')
				cat(paste('\nprocess file=',tmp))
				tmp	<- load(tmp)
				ans	<- ans$cens.AnyPos_T1
				ans[, method.risk:=method.risk]
				ans[, method.brl:=method.brl ]								
				ans
			},by='file']
	#
	#	set up factor legends
	#
	factors				<- project.athena.Fisheretal.sensitivity.factor.legend('m2Awmx')
	#
	#	set up t.period
	#
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	#
	#	CASCADE
	#		
	ctb			<- subset(cens.tables.bs, method.brl=='3pa1H1.48C2V100bInfT7')
	setkey(ctb, stat, t.period, risk, factor)
	ctb			<- unique(ctb)
	ct			<- subset(cens.tables, method.brl=='3pa1H1.48C2V100bInfT7')
	setkey(ct, stat, t.period, risk, factor)
	ct			<- unique(ct)	
	ctn			<- subset(cens.Patient.n, grepl('X.msm',stat) & method.brl=='3pa1H1.48C2V100bInfT7')
	setkey(ctn, stat, t.period)
	ctn			<- unique(ctn)
	plot.file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150304_CensoringFraction.pdf'
	ct1			<- project.athena.Fisheretal.censoring.model(ct, ctb, ctn, plot.file=plot.file, factors=factors)	
	#
	# illustrate extent of censoring
	#	
	tmp				<- subset(cens.AnyPos_T1, stage%in%c('U.2','UA.2'))
	set(tmp, tmp[, which(stage=='UA.2')], 'stage', 'Undiagnosed,\n Confirmed recent infection\n at diagnosis')
	set(tmp, tmp[, which(stage=='U.2')], 'stage', 'Undiagnosed,\n Unconfirmed chronic infection')
	ggplot( tmp, aes(x=t.AnyPos_T1, fill=stage)) + geom_histogram(binwidth=0.25) +
			labs(y='potential transmitters\nto recipient MSM diagnosed in 06/06-07/12\n(#)', x='time of diagnosis of potential transmitter', fill='', colour='') +
			scale_fill_manual(values=c("#990000","#FDBB84"), guide = FALSE) + geom_vline(xintercept = c(2006.5, 2008.)) +
			scale_y_continuous(breaks=seq(0,2000,100), minor_breaks=NULL, expand=c(0,5)) + 
			scale_x_continuous(breaks=seq(2004, 2015, 2), minor_breaks=seq(2004, 2015, 0.5), expand=c(0,0)) +
			facet_grid(. ~ stage, margins=FALSE) + theme_bw() + theme(panel.grid.major= element_line(colour="grey70", size=0.4), panel.grid.minor= element_line(colour="grey70", size=0.4))	
	plot.file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150305_censoringUandUA.pdf'
	ggsave(file=plot.file, w=10, h=4)
	#
	dfc			<- data.table(expand.grid(tC=seq(2008,2013,0.05), factor=unique(tmp$stage), stringsAsFactors=FALSE))
	dfc			<- dfc[,  {
				z	<- tmp$t.AnyPos_T1[ tmp$stage==factor]
				list( pnc=mean(z<=tC) )	
			}, by=c('factor','tC')]
	ggplot(dfc, aes(x=tC, ymax=100*(1-pnc), ymin=0, group=factor, fill=factor)) + 
			geom_ribbon() + 
			scale_fill_manual(values=c("#990000","#FDBB84"), guide = FALSE) +
			scale_y_continuous(breaks=seq(0,100,10), minor_breaks=NULL) + 
			scale_x_continuous(breaks=seq(2004, 2015, 2), minor_breaks=seq(2004, 2015, 0.5)) +			
			facet_grid(. ~ factor, margins=FALSE) + theme_bw() +
			labs(y='censored\npotential transmission intervals\n(%)', x='hypothetical time of database closure') +
			theme(panel.grid.major= element_line(colour="grey70", size=0.4), panel.grid.minor= element_line(colour="grey70", size=0.4))
	plot.file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150305_censoringUandUApc.pdf'
	ggsave(file=plot.file, w=10, h=3)
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
project.athena.Fisheretal.composition.potentialintervals.table<- function()
{
	require(data.table)
	require(ape)
	require(RColorBrewer)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_150319",sep='/')	
	outdir					<- paste(DATA,"fisheretal_150319",sep='/')	
		
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.table.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))	
	dft					<- subset(runs.table, !grepl('ARTstarted|GroupsUDA',method.risk) & stat%in%c('X.msm','nRec'))	
	dft					<- subset(dft, grepl('bInf', method.brl))
	dft[, t.period:= dft[, as.character(factor)]]	
	set(dft, NULL, 't.period', dft[, substring(t.period, nchar(t.period))] )
	set(dft, NULL, 'factor', dft[, as.character(factor)] )
	set(dft, NULL, 'factor', dft[, substr(factor, 1, nchar(factor)-2)] )
	set(dft, NULL, 'p', dft[, round(100*p,d=1)])
	#	proportions
	dftp				<- dcast.data.table(subset(dft, stat=='X.msm'), factor+method.brl~t.period, value.var='p')	
	#	total per recipient MSM
	tmp					<- subset(dft, stat=='X.msm')[, list(n=sum(n)), by=c('method.brl','t.period')]
	tmp2				<- unique(subset(dft, stat=='nRec', select=c(method.brl, t.period, n)))
	tmp					<- merge(tmp, tmp2, by=c('method.brl','t.period'))
	set(tmp, NULL, 'n', tmp[, n.x/n.y])
	tmp					<- dcast.data.table(tmp, method.brl~t.period, value.var='n')
	tmp[, factor:='']
	#	put everything together
	dftp				<- rbind(dftp, tmp, use.names=TRUE)
	tmp					<- c("","UA","UAna","U","DA","Dt.NA","Dtg500","Dtl500","Dtl350","ART.NotYetFirstSu","ART.vlNA","ART.suA.N","ART.suA.Y1","ART.suA.Y2","Lost")
	set( dftp, NULL, 'factor', dftp[, factor(factor, levels=tmp, labels=tmp)] )
	setkey(dftp, method.brl, factor)
	ans	<- rbind(	dcast.data.table( subset(dftp, select=c('method.brl', 'factor', '1')), factor~method.brl, value.var='1'),
					dcast.data.table( subset(dftp, select=c('method.brl', 'factor', '2')), factor~method.brl, value.var='2'),	
					dcast.data.table( subset(dftp, select=c('method.brl', 'factor', '3')), factor~method.brl, value.var='3'),
					dcast.data.table( subset(dftp, select=c('method.brl', 'factor', '4')), factor~method.brl, value.var='4')	)
	#	inf estimates	
	tmp			<- subset(ans, select=c('factor','3pa1H1.48C2V100bInfT7','3pa1H1.94C2V100bInfT7','3pa1H1.09C2V100bInfT7'))
	file		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150303_PoIntervals_InfTime.csv'
	write.csv(tmp, file=file, eol='\r\n')
}
######################################################################################
project.athena.Fisheretal.composition.problintervals.enrichment<- function()
{
	require(data.table)
	require(ape)
	require(RColorBrewer)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_150319",sep='/')	
	outdir					<- paste(DATA,"fisheretal_150319",sep='/')		
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"
	method.RISK				<- 'm2Awmx'
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	file					<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.table.Rdata", sep='')		
	readAttempt				<- try(suppressWarnings(load(file)))	
	#
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.5, 2006.5, 2008, 2009.5), t.period.max = c(2006.45, 2007.99, 2009.45, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )	
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-7','2008-1','2009-7'), labels=c('96/07','\n\n06/07','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-6','2007-12','2009-6','2010-12'), labels=c('06/06','07/12','09/06','10/12'))])
	#		
	runs.table			<- subset(runs.table, !grepl('ARTstarted|GroupsUDA',method.risk) & grepl(method.RISK, method.risk))
	#
	tmp					<- data.table(	method.brl=c(	"3pa1H1.48C2V100bInfT7", "3pa1H1.94C2V100bInfT7", "3pa1H1.09C2V100bInfT7",					
														"3pa1H1.48C3V100bInfT7", "3pa1H1.48C1V100bInfT7",
														"3pa1H1.48C2V100b0.02T7", "3pa1H1.48C2V100b0.04T7", 
														"3pa1H1.48C3V100bInfs0.85T7", "3pa1H1.48C2V100bInfs0.85T7", "3pa1H1.48C1V100bInfs0.85T7",
														"3pa1H1.48C3V100bInfs0.7T7", "3pa1H1.48C2V100bInfs0.7T7", "3pa1H1.48C1V100bInfs0.7T7"), 
										method.legend=c( 'central estimate of HIV infection times\ncentral phylogenetic exclusion criteria',
														'lower estimate of HIV infection times',
														'upper estimate of HIV infection times',
														'central phylogenetic exclusion criteria\nand genetic distance < 2%',
														'central phylogenetic exclusion criteria\nand genetic distance < 4%',
														'phylogenetic exclusion criteria\ncoalescent compatibility < 30%\nclade frequency < 80%',
														'phylogenetic exclusion criteria\ncoalescent compatibility < 10%\nclade frequency < 80%',
														'phylogenetic exclusion criteria\ncoalescent compatibility < 30%\nclade frequency < 85%',
														'phylogenetic exclusion criteria\ncoalescent compatibility < 20%\nclade frequency < 85%',
														'phylogenetic exclusion criteria\ncoalescent compatibility < 10%\nclade frequency < 85%',
														'phylogenetic exclusion criteria\ncoalescent compatibility < 30%\nclade frequency < 70%',
														'phylogenetic exclusion criteria\ncoalescent compatibility < 20%\nclade frequency < 70%',
														'phylogenetic exclusion criteria\ncoalescent compatibility < 10%\nclade frequency < 70%'))										 
	tmp				<- tmp[c(1:3,6:13),]
	set(tmp, NULL, 'method.legend', tmp[, factor(method.legend, levels=method.legend, labels=method.legend)])
	df				<- merge(runs.table, tmp, by='method.brl')	
	df[, t.period:= df[, as.character(factor)]]		
	set(df, NULL, 'factor', df[, as.character(factor)])	
	df[, t.period:=df[, substr(factor, nchar(factor), nchar(factor))]]
	set(df, NULL, 'factor', df[, substr(factor, 1, nchar(factor)-2)])		
	set(df, NULL, 'factor', df[, factor(factor, levels=factors[, levels(factor)])])
	
	#	prop potential transmission intervals 
	tmp		<- subset(df, stat%in%c('X.msm.e0cp','nRec'))	
	tmp		<- dcast.data.table(tmp, method.legend+method.brl+t.period+factor+method.risk~stat, value.var='n')
	tmp[, n:= X.msm.e0cp/8/nRec*100]	
	tmp		<- tmp[, list(n.coh=sum(n)), by=c('factor','method.legend','method.brl')]
	tmp		<- merge(tmp, tmp[, list(factor=factor, p.coh=n.coh/sum(n.coh)), by=c('method.legend','method.brl')], by=c('factor','method.legend','method.brl'))
	ans		<- subset(tmp, select=c(method.brl, method.legend, factor,  p.coh))
	#	prop likely transmission intervals 
	tmp		<- subset(df, stat%in%c('YX','nRecLkl'))	
	tmp		<- dcast.data.table(tmp, method.legend+method.brl+t.period+factor+method.risk~stat, value.var='n')
	tmp[, n:= YX/nRecLkl*100]
	tmp		<- tmp[, list(n.lkl=sum(n)), by=c('factor','method.legend','method.brl')]
	#	add prop likely transmission intervals + expected missing intervals per 100 'recipient with lkl transmitter' in time period
	tmp2	<- subset(df, stat%in%c('YX','Sx.e0cp','nRecLkl'))	
	tmp2	<- dcast.data.table(tmp2, method.legend+method.brl+t.period+factor+method.risk~stat, value.var='n')	
	tmp2[, nm:= (YX+Sx.e0cp/8)/nRecLkl*100]
	tmp2	<- tmp2[, list(n.lklm=sum(nm)), by=c('factor','method.legend','method.brl')]	
	tmp		<- merge(subset(tmp, select=c(method.legend, method.brl, factor, n.lkl)), subset(tmp2, select=c(method.legend, method.brl, factor, n.lklm)), by=c('factor','method.legend','method.brl'))	
	tmp		<- tmp[, list(factor=factor, p.lkl=(n.lkl+n.lklm)/sum(n.lkl+n.lklm)), by=c('method.legend','method.brl')]
	#	
	ans		<- merge(ans, subset(tmp, select=c(method.legend, method.brl, factor, p.lkl)), by=c('factor','method.legend','method.brl'))
	#	plot
	ans		<- melt(ans, measure.vars=c('p.coh','p.lkl'))
	factors	<- project.athena.Fisheretal.sensitivity.factor.legend(method.RISK)
	ans		<- merge(ans, factors, by='factor')	
	setkey(ans, factor.legend)
	ggplot(ans, aes(x=interaction(variable, factor), y=value*100, colour=variable, fill=factor.legend)) + 
		geom_bar(stat='identity', size=1, width=0.7, position = position_dodge(width = 0.8)) +
		scale_fill_manual(values=ans[, unique(factor.color)], guide=FALSE) +
		scale_colour_manual(values=c('transparent','black'), guide=FALSE) +
		scale_y_continuous(breaks=seq(0,50,10), minor_breaks=seq(0,50,2)) +
		theme_bw() + 
		geom_text( aes(y=-2, label='*'), size=7) +		
		labs(x='stage in HIV infection and care continuum', y='Transmission intervals\n(%)') +
		facet_wrap(~method.legend, ncol=3) +
		theme(axis.text=element_text(size=14), axis.title=element_text(size=14), legend.key.size=unit(11,'mm'), legend.position = "bottom", legend.box = "vertical", axis.ticks.x=element_blank(), axis.text.x=element_blank(), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major.x=element_blank(), panel.grid.major=element_line(colour="grey70", size=0.7)) 
	file		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150309_PrMIntervalsEnrichment.pdf'
	ggsave(file=file, width=14, height=14)


}
######################################################################################
project.athena.Fisheretal.composition.problintervals.withmissing.table.generate<- function(runs.table, file)
{
	dfm					<- subset(runs.table, stat%in%c('Sx.e0cp','nRecLkl'))
	tmp					<- subset(dfm, stat=='Sx.e0cp', c(method.brl, t.period, factor, l95.bs,  u95.bs))
	dfm					<- dcast.data.table(dfm, t.period+factor+method.risk+method.brl~stat, value.var='n')
	dfm					<- merge(dfm, tmp, by=c('method.brl','t.period','factor'))
	dfm[, Sx.e0cp:= Sx.e0cp/8/nRecLkl]
	dfm[, l95.bs:= l95.bs/8/nRecLkl]
	dfm[, u95.bs:= u95.bs/8/nRecLkl]
	#	observed
	dft					<- subset(runs.table, stat%in%c('YX','nRecLkl'))
	dft					<- dcast.data.table(dft, t.period+factor+method.risk+method.brl~stat, value.var='n')
	set(dft, NULL, 'YX', dft[, YX/nRecLkl])
	#	both
	dfb					<- merge( 	subset(dft, select=c(t.period, factor, method.risk, method.brl, YX)), 
			subset(dfm, select=c(t.period, factor, method.risk, method.brl, Sx.e0cp, l95.bs, u95.bs)), by=c('t.period','factor','method.risk','method.brl') )
	set(dfb, NULL, 'n', dfb[, YX+Sx.e0cp])
	set(dfb, NULL, 'l95.bs', dfb[, YX+l95.bs])
	set(dfb, NULL, 'u95.bs', dfb[, YX+u95.bs])
	#	proportions per time period			
	tmp					<- dfb[, list(factor=factor, stat='p', central= round(100*n/sum(n),d=1), l95.bs= round(100*l95.bs/sum(l95.bs),d=1), u95.bs= round(100*u95.bs/sum(u95.bs),d=1)), by=c('t.period','method.risk','method.brl')]
	#	total per recipient MSM
	dfb					<- dfb[, list(factor='', stat='n', central=round(sum(n),d=3), l95.bs=round(sum(l95.bs),d=3), u95.bs=round(sum(u95.bs),d=3)), by=c('t.period','method.risk','method.brl')]	
	dfb					<- rbind(dfb, tmp, use.names=TRUE)
	dfb[, char:= paste(central,' (',l95.bs,'-',u95.bs,')',sep='')]
	tmp					<- dfb[, which(stat=='p')]
	set(dfb, tmp, 'char', dfb[tmp, as.character(central)])
	#	put everything together
	tmp					<- c("","UA","UAna","U","DA","Dt.NA","Dtg500","Dtl500","Dtl350","ART.NotYetFirstSu","ART.vlNA","ART.suA.N","ART.suA.Y1","ART.suA.Y2","Lost")
	set( dfb, NULL, 'factor', dfb[, factor(factor, levels=tmp, labels=tmp)] )
	setkey(dfb, method.brl, factor)
	
	ans	<- rbind(	dcast.data.table( subset(dfb, t.period==1, select=c('method.brl', 'factor', 'char')), factor~method.brl, value.var='char'),
			dcast.data.table( subset(dfb, t.period==2, select=c('method.brl', 'factor', 'char')), factor~method.brl, value.var='char'),	
			dcast.data.table( subset(dfb, t.period==3, select=c('method.brl', 'factor', 'char')), factor~method.brl, value.var='char'),
			dcast.data.table( subset(dfb, t.period==4, select=c('method.brl', 'factor', 'char')), factor~method.brl, value.var='char')	)
	ans	<- subset(ans, select=c(1,3,4,2))	
	write.csv(ans, file=file, eol='\r\n')
}
######################################################################################
project.athena.Fisheretal.composition.potentialintervals.withmissing.table.generate<- function(runs.table, file)
{
	#	central
	dft					<- subset(runs.table, stat%in%c('X.msm.e0cp','nRec'))
	dft					<- dcast.data.table(dft, t.period+factor+method.risk+method.brl~stat, value.var='n')
	set(dft, NULL, 'X.msm.e0cp', dft[, X.msm.e0cp/8/nRec])
	#	lower
	tmp					<- subset(runs.table, stat=='X.msm.e0cp')
	tmp					<- dcast.data.table(tmp, t.period+factor+method.risk+method.brl~stat, value.var='l95.bs')
	setnames(tmp, 'X.msm.e0cp', 'l95.bs')
	dft					<- merge(dft, tmp, by=c('t.period','factor','method.risk','method.brl'))
	set(dft, NULL, 'l95.bs', dft[, l95.bs/8/nRec])
	#	upper
	tmp					<- subset(runs.table, stat=='X.msm.e0cp')
	tmp					<- dcast.data.table(tmp, t.period+factor+method.risk+method.brl~stat, value.var='u95.bs')
	setnames(tmp, 'X.msm.e0cp', 'u95.bs')
	dft					<- merge(dft, tmp, by=c('t.period','factor','method.risk','method.brl'))
	set(dft, NULL, 'u95.bs', dft[, u95.bs/8/nRec])
	#	proportions per time period			
	tmp					<- dft[, list(factor=factor, stat='p', central= round(100*X.msm.e0cp/sum(X.msm.e0cp),d=1)), by=c('t.period','method.risk','method.brl')]
	#	total per recipient MSM
	dft					<- dft[, list(factor='', stat='n', central=round(sum(X.msm.e0cp),d=0), l95.bs=round(sum(l95.bs),d=0), u95.bs=round(sum(u95.bs),d=0)), by=c('t.period','method.risk','method.brl')]
	dft					<- rbind(dft, tmp, use.names=TRUE, fill=TRUE)
	dft[, char:= paste(central,' (',l95.bs,'-',u95.bs,')',sep='')]
	tmp					<- dft[, which(stat=='p')]
	set(dft, tmp, 'char', dft[tmp, as.character(central)])
	#	put everything together
	tmp					<- c("","UA","UAna","U","DA","Dt.NA","Dtg500","Dtl500","Dtl350","ART.NotYetFirstSu","ART.vlNA","ART.suA.N","ART.suA.Y1","ART.suA.Y2","Lost")
	set( dft, NULL, 'factor', dft[, factor(factor, levels=tmp, labels=tmp)] )
	setkey(dft, method.brl, factor)
	
	ans	<- rbind(	dcast.data.table( subset(dft, t.period==1, select=c('method.brl', 'factor', 'char')), factor~method.brl, value.var='char'),
			dcast.data.table( subset(dft, t.period==2, select=c('method.brl', 'factor', 'char')), factor~method.brl, value.var='char'),	
			dcast.data.table( subset(dft, t.period==3, select=c('method.brl', 'factor', 'char')), factor~method.brl, value.var='char'),
			dcast.data.table( subset(dft, t.period==4, select=c('method.brl', 'factor', 'char')), factor~method.brl, value.var='char')	)
	ans	<- subset(ans, select=c(1,3,4,2))	
	write.csv(ans, file=file, eol='\r\n')
}
######################################################################################
project.athena.Fisheretal.composition.potentialintervals.withmissing.table<- function()
{
	require(data.table)
	require(ape)
	require(RColorBrewer)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_150319",sep='/')	
	outdir					<- paste(DATA,"fisheretal_150319",sep='/')	
	
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.table.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))		
	runs.table[, t.period:= runs.table[, as.character(factor)]]	
	set(runs.table, NULL, 't.period', runs.table[, substring(t.period, nchar(t.period))] )
	set(runs.table, NULL, 'factor', runs.table[, as.character(factor)] )
	set(runs.table, NULL, 'factor', runs.table[, substr(factor, 1, nchar(factor)-2)] )
	set(runs.table, NULL, 'p', runs.table[, round(100*p,d=1)])
	#
	#	first table
	#
	runs.table			<- subset(runs.table, !grepl('ARTstarted|GroupsUDA',method.risk))
	tmp					<- subset(runs.table, grepl('bInfT7', method.brl) & grepl('C2', method.brl))	
	file				<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150327_PotAdjIntervalsByCoalComp.csv'
	project.athena.Fisheretal.composition.potentialintervals.withmissing.table.generate(tmp, file)
	
}
######################################################################################
project.athena.Fisheretal.composition.problintervals.withmissing.table<- function()
{
	require(data.table)
	require(ape)
	require(RColorBrewer)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_150319",sep='/')	
	outdir					<- paste(DATA,"fisheretal_150319",sep='/')	
	
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.table.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))		
	runs.table[, t.period:= runs.table[, as.character(factor)]]	
	set(runs.table, NULL, 't.period', runs.table[, substring(t.period, nchar(t.period))] )
	set(runs.table, NULL, 'factor', runs.table[, as.character(factor)] )
	set(runs.table, NULL, 'factor', runs.table[, substr(factor, 1, nchar(factor)-2)] )
	set(runs.table, NULL, 'p', runs.table[, round(100*p,d=1)])
	#
	#	first table
	#
	runs.table			<- subset(runs.table, !grepl('ARTstarted|GroupsUDA',method.risk))
	tmp					<- subset(runs.table, grepl('bInfT7', method.brl) & grepl('H1.48', method.brl))	
	file				<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150308_PrMIntervals.csv'
	project.athena.Fisheretal.composition.problintervals.withmissing.table.generate(tmp, file)
	#
	#	second table
	#	
	tmp					<- subset(runs.table, grepl('C2', method.brl) & grepl('bInf', method.brl) & grepl('H1.48', method.brl))	
	file				<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150308_PrMIntervalsByCladeFreq.csv'
	project.athena.Fisheretal.composition.problintervals.withmissing.table.generate(tmp, file)	
}
######################################################################################
project.athena.Fisheretal.composition.transprob.uncertainty<- function()	
{
	method.RISK		<- 'm2Awmx.wtn'
	#	t.period
	tperiod.info	<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-5','2008-1','2009-7'), labels=c('96/07','\n\n06/05','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-4','2007-12','2009-6','2010-12'), labels=c('06/04','07/12','09/06','10/12'))])			
	#	get trm.p for data
	indir			<- paste(DATA,"fisheretal_150319",sep='/')
	file			<- paste(indir,'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Wed_Dec_18_11:37:00_2013_method.trmp.BS.Rdata',sep='/')
	load(file)		
	trm.p			<- subset(trm.p, grepl(method.RISK,method.risk) & grepl('3pa1H1.48C2V100bInfT7',method.brl))
	trm.p[, t.period:= trm.p[, substring(factor, nchar(factor))]]
	set(trm.p, NULL, 'factor', trm.p[, substr(factor, 1, nchar(factor)-2)])
	#	normalize tp per patient
	#subset(trm.p, BS==0 & Patient=='M12884')[, list(factor=factor, tp= ys/sum(yn+ym.e0cp)/yn  ), by='Patient'  ]
	trm.pu			<- trm.p[, list(factor=factor, tp= ys/sum(ys+ys.e0cp)/yn  ), by=c('BS','Patient')  ]
	trm.pu			<- subset(trm.pu, !is.na(tp))
	#set(trm.pu, trm.pu[, which(is.na(tp))], 'tp', 0)		don t count zeros -- we look at pos transmission intervals throughout
	#	match recipients in bootstrap to recipient in BS==0
	tmp				<- subset(trm.pu, BS>0)
	setnames(tmp, 'Patient', 'Patient.bs')	
	tmp[, Patient:= tmp[, regmatches(Patient.bs, regexpr('.*_bs',Patient.bs))]]
	set(tmp, NULL, 'Patient', tmp[, substr(Patient, 1, nchar(Patient)-3)])	
	tmp				<- tmp[, list(Patient.bs=Patient.bs[1], tp=tp[1]), by=c('BS','factor','Patient')]	
	trm.pu			<- rbind(trm.pu, subset(tmp, select=c(BS, Patient, factor, tp)))
	#	take selection of recipients
	tmp				<- trm.pu[, list(nBS=length(BS)), by=c('factor','Patient')]
	setkey(tmp, factor, nBS)
	tmp				<- tmp[,  list(Patient=rev(Patient)[1:min(40,length(Patient))]), by='factor']
	trm.pu			<- merge(tmp, trm.pu, by=c('factor','Patient'))
	#	do boxplots
	tmp				<- subset(trm.pu, factor%in%c('UA',"UAna","ART.suA.N","ART.suA.Y2"))
	factors			<- project.athena.Fisheretal.sensitivity.factor.legend(method.RISK)
	tmp				<- merge(tmp, factors, by='factor')
	#	something going wrong, rescale for now -- need to fix
	tmp[, list(tpmed= mean(tp)), by='factor']	
	setkey(tmp, factor.legend)
	ggplot(tmp, aes(x=factor(Patient), y=100*tp, colour=factor.legend, fill=factor.legend)) + geom_boxplot(data=subset(tmp,BS>0), outlier.size = 0.4) + 
			#geom_point(data=subset(tmp,BS==0), size=1.2, colour="#E41A1C") +			
			scale_y_continuous(expand=c(0,0)) +
			scale_colour_manual(name='', values=tmp[, unique(factor.color)], guide=FALSE) +
			scale_fill_manual(name='', values=tmp[, unique(factor.color)], guide=FALSE) +
			stat_summary(geom= "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
			theme_bw() +
			theme(strip.text=element_blank(), strip.background=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), panel.grid.major.x=element_blank(), panel.grid.major=element_line(colour="grey60", size=0.4), panel.grid.minor=element_line(colour="grey60", size=0.4) ) +			
			facet_wrap(~factor.legend, ncol=2, scales='free') +
			labs(	x='transmission intervals', y='transmission probability\n(%)')
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150311_TrProbUncertainty.pdf'
	ggsave(file=file, w=10, h=8)
}
######################################################################################
project.athena.Fisheretal.composition.transprob.byexclusioncriteria<- function()	
{
	method.RISK		<- 'm2Awmx'
	if(grepl('m2A', method.RISK))
		stage.id	<- 'CD4a'
	if(grepl('m2C', method.RISK))
		stage.id	<- 'CD4c'	
	#	t.period
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-5','2008-1','2009-7'), labels=c('96/07','\n\n06/05','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-4','2007-12','2009-6','2010-12'), labels=c('06/04','07/12','09/06','10/12'))])		
	#	load YX
	indir			<- paste(DATA,"fisheretal_150319",sep='/')
	infiles			<- list.files(indir, pattern='wtn.tp1.R$')
	infiles			<- infiles[ !grepl('b0.02|b0.04',infiles) ]
	infiles			<- data.table(FILE=infiles)
	infiles[, method.brl:= infiles[, regmatches(FILE, regexpr('Yscore[^_]+',FILE))]]
	set(infiles, NULL, 'method.brl', infiles[, substring(method.brl, 7)])
	tmp				<- lapply(seq_len(nrow(infiles)), function(i)
			{
				file			<- paste(indir, '/', infiles[i,FILE], sep='')
				load(file)
				YX				<- subset(ans$YXf, select=c(t.period, t, t.Patient, Patient, score.Y, CD4a, CD4b, CD4c))
				setnames( YX, stage.id, 'factor')	
				set(YX, YX[, which(as.numeric(t.period)>4)],'t.period','4')			
				YX[, method.brl:=infiles[i,method.brl]]
				YX
			})
	YX				<- do.call('rbind', tmp)
	#	get trm.p for data
	file			<- paste(indir,'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Wed_Dec_18_11:37:00_2013_method.trmp.BS0.Rdata',sep='/')
	load(file)	
	trm.p			<- subset(trm.p, BS==0 & grepl(method.RISK,method.risk))
	trm.p[, t.period:= trm.p[, substring(factor, nchar(factor))]]
	set(trm.p, NULL, 'factor', trm.p[, substr(factor, 1, nchar(factor)-2)])
	#	get denominator for every patient
	tmp		<- trm.p[, list( ys.e0cp=sum(ys.e0cp)+sum(ys) ) , by=c('method.brl','Patient','t.period')]
	#	merge with YXf and plot probs
	#	YX[, list(Patient=setdiff(  Patient, tmp[['Patient']][ tmp[['method.brl']]==method.brl ] )), by='method.brl']
	# 	something wrong with "3pa1H1.48C1V100bInfs0.7T7" -- exclude
	YX		<- merge(YX, tmp, by=c('method.brl','t.period','Patient'))
	YX[, tp:= score.Y/ys.e0cp]
	#	check
	#merge(YX, YX[, list(yscheck=sum(score.Y)), by=c('t.period','Patient')], by=c('COAL','BRL','t.period','Patient'))
	#	merge tperiod.long
	YX		<- merge(YX, tperiod.info, by='t.period')
	YX[, t.period.long:= paste(t.period.min, '-', t.period.max,sep='')]
	#	merge factor legends
	factors	<- project.athena.Fisheretal.sensitivity.factor.legend(method.RISK)
	YX		<- merge(YX, factors, by='factor')	
	#	merge run legends
	tmp		<- data.table(	method.brl=c(	"3pa1H1.48C2V100bInfT7", "3pa1H1.94C2V100bInfT7", "3pa1H1.09C2V100bInfT7",					
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
	YX		<- merge(YX, tmp, by='method.brl')
	setkey(YX, method.brl, factor.legend, t.period.long)
	#
	#	boxplot of transmission prob (observed only)
	#
	ggplot(YX, aes(x=factor.legend, y=tp*100, colour=factor.legend, fill=factor.legend))  + 
			labs(x='stage in HIV infection and care continuum', y='phylogenetic transmission probability\nper observed, probable interval\n(%)') +
			geom_boxplot(outlier.shape = NA) +
			scale_colour_manual(name='', values=YX[, unique(factor.color)], guide=FALSE) +
			scale_fill_manual(name='', values=YX[, unique(factor.color)], guide=FALSE) +
			scale_y_continuous(breaks=seq(0, 50, 5), minor_breaks=seq(0, 50, 1)) +
			scale_x_discrete(breaks=NULL) +	
			#coord_trans(ytrans="log10", limy=c(0.001, 0.3)) +
			coord_cartesian(ylim=c(0, 20)) +
			stat_summary(geom= "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +		
			theme_bw() + theme(panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major.x=element_blank(), panel.grid.major=element_line(colour="grey70", size=0.7)) +
			facet_wrap(~method.legend, ncol=3)
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150312_TrProbByCriteria.pdf'
	ggsave(file=file, w=10, h=10)
			
}
######################################################################################
project.athena.Fisheretal.composition.riskratio.exclusioncriteria<- function()	
{
	#	t.period
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-5','2008-1','2009-7'), labels=c('96/07','\n06/07','08/01','\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-4','2007-12','2009-6','2010-12'), labels=c('06/06','07/12','09/06','10/12'))])
	tperiod.info[, t.period.long:= tperiod.info[,paste(t.period.min, '-', t.period.max,sep='')]]
	set(tperiod.info, NULL, 't.period.long', tperiod.info[, factor(t.period.long, levels=t.period.long, labels=t.period.long)])		
	#	load risk tables
	indir				<- paste(DATA,"fisheretal_150319",sep='/')
	outdir				<- paste(DATA,"fisheretal_150319",sep='/')			
	infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	insignat			<- "Wed_Dec_18_11:37:00_2013"	
	
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))
	runs.risk			<- subset(runs.risk, !grepl('ARTstarted',method.risk) & !grepl('GroupsUDA',method.risk) )
	method.DENOM		<- 'SEQ'		
	method.WEIGHT		<- ''
	method.BRL			<- c(  "3pa1H1.09C2V100bInfT7", "3pa1H1.48C1V100bInfT7", "3pa1H1.48C1V100bInfs0.7T7", "3pa1H1.48C1V100bInfs0.85T7", 
			"3pa1H1.48C2V100b0.02T7", "3pa1H1.48C2V100b0.04T7", "3pa1H1.48C2V100bInfT7", "3pa1H1.48C2V100bInfs0.7T7",  
			"3pa1H1.48C2V100bInfs0.85T7", "3pa1H1.48C3V100bInfT7", "3pa1H1.48C3V100bInfs0.7T7", "3pa1H1.48C3V100bInfs0.85T7",
			"3pa1H1.94C2V100bInfT7")
	method.DATING		<- 'sasky'	
	stat.select			<- 'RR.raw.e0cp'
	outfile				<- infile
	method.RISK			<- "m2Awmx.wtn.tp"	
	run.tp				<- subset(runs.risk, 	(is.na(t.period) | t.period!=0) & method.denom==method.DENOM & method.nodectime=='any' & 
									grepl(method.RISK,method.risk) & method.dating==method.DATING & 
									stat==stat.select & grepl('Dtg500', factor.ref) & grepl('ART',factor)	)
	#	merge methods legend
	tmp					<- data.table(	method.brl=c(	"3pa1H1.48C2V100bInfT7", "3pa1H1.94C2V100bInfT7", "3pa1H1.09C2V100bInfT7",					
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
	run.tp				<- merge(run.tp, tmp, by='method.brl')
	#	prepare
	run.tp[, m50.bs:=NULL]
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])	
	#	merge factors
	factors				<- project.athena.Fisheretal.sensitivity.factor.legend(method.RISK)
	run.tp				<- merge(run.tp, subset(factors, select=factors[, which(colnames(factors)!='method.risk')]), by='factor')			
	set(run.tp, NULL, 'factor', run.tp[, factor(factor, levels=factors[, levels(factor)])])	
	#	merge t.period
	run.tp				<- merge(run.tp, tperiod.info, by='t.period')	
	#
	#	plot
	#	
	setkey(run.tp, method.legend, factor.legend, t.period.long)
	ggplot(run.tp, aes(x=t.period.long, y=100*v, ymin=100*l95.bs, ymax=100*u95.bs, fill=factor.legend, colour=factor.legend, group=factor.legend)) + 
			geom_hline(yintercept=100) + 
			geom_point(position=position_dodge(.5) , shape=18) + 
			geom_errorbar(position=position_dodge(.5), width=0.5) +
			scale_x_discrete(expand=c(0.03,0.03)) +
			scale_y_continuous(breaks=seq(0,200,50), minor_breaks=seq(0,200,10), limit=c(0, 101)) + 
			scale_fill_manual(name='from cascade stage', values= run.tp[, unique(factor.color)], guide=FALSE) + 
			scale_colour_manual(name='from cascade stage', values= run.tp[, unique(factor.color)], guide=FALSE) +
			theme_bw() + theme(panel.margin.y=unit(1.25, "lines"), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_line(colour="grey70", size=0.4), panel.grid.major.x=element_blank(), panel.grid.major.y=element_line(colour="grey70", size=0.4)) +
			labs(x="", y="Transmission risk ratio\ncompared to diagnosed, untreated men with CD4>500\n(%)") +
			facet_wrap(~method.legend, ncol=3, scales='free')
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150316_RRByExclusionCriteria.pdf'
	ggsave(file=file, w=10, h=10)
	#
	

}
project.athena.Fisheretal.composition.prop.overall<- function()
{
	require(data.table)
	require(ape)
	require(grid)
	require(reshape2)
	require(ggplot2)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_150319",sep='/')
	outdir					<- paste(DATA,"fisheretal_150319",sep='/')			
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period	
	file					<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.Rdata", sep='')		
	readAttempt				<- try(suppressWarnings(load(file)))	
	#
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.5, 2006.5, 2008, 2009.5), t.period.max = c(2006.45, 2007.99, 2009.45, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )	
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-7','2008-1','2009-7'), labels=c('96/07','\n\n06/07','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-6','2007-12','2009-6','2010-12'), labels=c('06/06','07/12','09/06','10/12'))])	
	#
	method.DENOM	<- 'SEQ'		
	method.WEIGHT	<- ''
	method.DATING	<- 'sasky'	
	stat.select		<- 'P.raw.e0cp'
	outfile			<- infile
	method.BRL		<- '3pa1H1.94C2V100bInfT7'
	method.BRL		<- '3pa1H1.48C2V100bInfT7'
	factors			<- project.athena.Fisheretal.sensitivity.factor.legend('m2Awmx.wtn.tp')			
	method.RISK		<- 'm2Awmx.wtn.tp'
	#
	run.tp			<- subset(runs.risk,  t.period==0 & method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating==method.DATING & stat.select==stat  )	
	run.tp			<- subset(run.tp, !grepl('total', factor))
	run.tp[, m50.bs:=NULL]
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-4)])	
	set(run.tp, NULL, 'factor', run.tp[, factor(factor, levels=factors[, levels(factor)])])
	run.tp	<- merge(run.tp, subset(factors, select=factors[, which(colnames(factors)!='method.risk')]), by='factor')	
	run.tp[, FCT:='from phylogenetically probable\ntransmitters']
	
	ggplot(run.tp, aes(x=factor, y=100*v, fill=factor.legend, colour=factor.legend)) + labs(x="", y="Transmissions\n(%)") + 
			scale_y_continuous(breaks=seq(0,90,10), minor_breaks=seq(0,90,2), limits=c(0,35), expand=c(0,0)) + 
			scale_x_discrete(breaks=NULL, limits=run.tp[, intersect(levels(factor), as.character(factor))]) +
			scale_fill_manual(name='from cascade stage', values=run.tp[, unique(factor.color)], guide=FALSE) + 
			scale_colour_manual(name='from cascade stage', values = rep('black',run.tp[, length(unique(factor))]), guide=FALSE) +
			#scale_fill_brewer(palette='PRGn',name='from cascade stage') + scale_colour_manual(name='from cascade stage', values = rep('black',11)) +					
			#guides(colour=FALSE, fill = guide_legend(override.aes = list(size=5))) +
			theme_bw() + theme(strip.text=element_text(size=14), axis.text=element_text(size=14), axis.title=element_text(size=14), legend.key.size=unit(10.5,'mm'), plot.margin=unit(c(0,0,-5,0),"mm"), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.7)) + #coord_flip() +
			geom_bar(stat='identity', position='dodge') + geom_errorbar(aes(ymin=100*l95.bs, ymax=100*u95.bs), width=0.3, position=position_dodge(width=0.9)) 
	
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2014_MSMtransmission_ATHENA1303/160603_PropOverall.pdf'		
	ggsave(file=file, w=5, h=3.5)		
}
######################################################################################
project.athena.Fisheretal.composition.prop.exclusiongendistance<- function()	
{
	#	t.period
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-5','2008-1','2009-7'), labels=c('96/07','\n\n06/05','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-4','2007-12','2009-6','2010-12'), labels=c('06/04','07/12','09/06','10/12'))])		
	#	load risk tables
	indir				<- paste(DATA,"fisheretal_150319",sep='/')
	outdir				<- paste(DATA,"fisheretal_150319",sep='/')			
	infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	insignat			<- "Wed_Dec_18_11:37:00_2013"	
	
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))
	runs.risk			<- subset(runs.risk, !grepl('ARTstarted',method.risk) & !grepl('GroupsUDA',method.risk))
	method.DENOM		<- 'SEQ'		
	method.WEIGHT		<- ''
	method.BRL			<- c(  "3pa1H1.48C2V100b0.02T7", "3pa1H1.48C2V100b0.04T7")
	method.DATING		<- 'sasky'	
	stat.select			<- 'P.raw.e0cp'
	stat.select			<- 'P.raw'
	outfile				<- infile
	method.RISK			<- "m2Awmx.wtn.tp"	
	#	set up factor legends
	factors				<- project.athena.Fisheretal.sensitivity.factor.legend(method.RISK)
	
	run.tp				<- subset(runs.risk, 	(is.na(t.period) | t.period!=0) & method.denom==method.DENOM & method.nodectime=='any' & 
					grepl(method.RISK,method.risk) & method.dating==method.DATING & 
					(grepl('P.',stat,fixed=1) | stat=='P' | grepl('N.',stat,fixed=1) | stat=='N') )
	run.tp				<- subset(run.tp, method.brl%in%method.BRL)
	run.tp				<- subset(run.tp, !is.na(v) & stat%in%stat.select )
	#	prepare
	run.tp[, m50.bs:=NULL]
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])
	set(run.tp, NULL, 'factor', run.tp[, factor(factor, levels=factors[, levels(factor)])])	
	#	merge t.period
	run.tp				<- merge(run.tp, tperiod.info, by='t.period')
	run.tp[, t.period.long:= gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))]
	set(run.tp, NULL, 't.period.long', run.tp[,factor(t.period.long, levels= tperiod.info[, gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))])])	
	#	merge factors
	run.tp				<- merge(run.tp, subset(factors, select=factors[, which(colnames(factors)!='method.risk')]), by='factor')	
	#	merge methods legend
	tmp					<- data.table(	method.brl=c(	"3pa1H1.48C2V100b0.02T7", "3pa1H1.48C2V100b0.04T7"), 
										method.legend=c( 'central phylogenetic exclusion criteria\nand genetic distance < 2%',
														 'central phylogenetic exclusion criteria\nand genetic distance < 4%'))										 
	tmp					<- data.table(	method.brl=c(	"3pa1H1.48C2V100b0.02T7", "3pa1H1.48C2V100b0.04T7"), 
										method.legend=c( 'central phylogenetic exclusion criteria\nno censoring/sampling bias adjustments\nand genetic distance < 2%',
														 'central phylogenetic exclusion criteria\nno censoring/sampling bias adjustments\nand genetic distance < 4%'))										 										 
	set(tmp, NULL, 'method.legend', tmp[, factor(method.legend, levels=method.legend, labels=method.legend)])
	run.tp	<- merge(run.tp, tmp, by='method.brl')
	setkey(run.tp, factor.legend, t.period.long)
	
	tmp2	<- subset(tmp[1:2, ], select=method.brl)
	tmp2	<- merge(run.tp, tmp2, by='method.brl')
	ggplot(tmp2, aes(x=factor, y=v*100, fill=factor.legend, colour=factor.legend)) + 
			labs(x="stage in HIV infection and care continuum", y="Proportion of transmissions\n(%)") + 
			scale_y_continuous(breaks=seq(0,90,10), minor_breaks=seq(0,90,2), limit=c(0,110*run.tp[, round(max(u95.bs), d=2)]), expand=c(0,0)) + 
			scale_x_discrete(breaks=NULL, limits=run.tp[, intersect(levels(factor), as.character(factor))]) +
			scale_fill_manual(name='from cascade stage', values=run.tp[, unique(factor.color)], guide=FALSE) + 
			scale_colour_manual(name='from cascade stage', values = rep('black',run.tp[, length(unique(factor))]), guide=FALSE) +
			theme_bw() + theme(panel.margin = unit(0.8, "lines"), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.7)) + #coord_flip() +
			geom_bar(stat='identity',binwidth=1, position='dodge') + 
			geom_errorbar(aes(ymin=l95.bs*100, ymax=u95.bs*100), width=0.3, position=position_dodge(width=0.9))	+ 
			facet_grid(method.legend ~ t.period.long, margins=FALSE)
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150618_PropByInfectionTime.pdf'
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150618_PropByInfectionTimePraw.pdf'
	ggsave(file=file, w=10, h=7)		
}
######################################################################################
project.athena.Fisheretal.composition.prop.exclusioncriteria<- function()	
{
	#	t.period
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-5','2008-1','2009-7'), labels=c('96/07','\n\n06/05','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-4','2007-12','2009-6','2010-12'), labels=c('06/04','07/12','09/06','10/12'))])		
	#	load risk tables
	indir				<- paste(DATA,"fisheretal_150319",sep='/')
	outdir				<- paste(DATA,"fisheretal_150319",sep='/')			
	infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	insignat			<- "Wed_Dec_18_11:37:00_2013"	
	
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))
	runs.risk			<- subset(runs.risk, !grepl('ARTstarted',method.risk) & !grepl('GroupsUDA',method.risk))
	method.DENOM		<- 'SEQ'		
	method.WEIGHT		<- ''
	method.BRL			<- c(  "3pa1H1.09C2V100bInfT7", "3pa1H1.48C1V100bInfT7", "3pa1H1.48C1V100bInfs0.7T7", "3pa1H1.48C1V100bInfs0.85T7", 
 								"3pa1H1.48C2V100b0.02T7", "3pa1H1.48C2V100b0.04T7", "3pa1H1.48C2V100bInfT7", "3pa1H1.48C2V100bInfs0.7T7",  
 								"3pa1H1.48C2V100bInfs0.85T7", "3pa1H1.48C3V100bInfT7", "3pa1H1.48C3V100bInfs0.7T7", "3pa1H1.48C3V100bInfs0.85T7",
								"3pa1H1.94C2V100bInfT7")
	method.DATING		<- 'sasky'	
	stat.select			<- 'P.raw.e0cp'
	outfile				<- infile
	method.RISK			<- "m2Awmx.wtn.tp"	
	run.tp				<- subset(runs.risk, 	(is.na(t.period) | t.period!=0) & method.denom==method.DENOM & method.nodectime=='any' & 
												grepl(method.RISK,method.risk) & method.dating==method.DATING & 
												(grepl('P.',stat,fixed=1) | stat=='P' | grepl('N.',stat,fixed=1) | stat=='N') )
	run.tp				<- subset(run.tp, method.brl%in%method.BRL)
	run.tp				<- subset(run.tp, !is.na(v) & stat%in%stat.select )
	#	prepare
	run.tp[, m50.bs:=NULL]
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])
	set(run.tp, NULL, 'factor', run.tp[, factor(factor, levels=factors[, levels(factor)])])	
	#	merge t.period
	run.tp				<- merge(run.tp, tperiod.info, by='t.period')
	run.tp[, t.period.long:= gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))]
	set(run.tp, NULL, 't.period.long', run.tp[,factor(t.period.long, levels= tperiod.info[, gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))])])	
	#	merge factors
	#	set up factor legends
	factors				<- project.athena.Fisheretal.sensitivity.factor.legend(method.RISK)
	run.tp				<- merge(run.tp, subset(factors, select=factors[, which(colnames(factors)!='method.risk')]), by='factor')	
	#	merge methods legend
	tmp					<- data.table(	method.brl=c(	"3pa1H1.48C2V100bInfT7", "3pa1H1.94C2V100bInfT7", "3pa1H1.09C2V100bInfT7",					
														"3pa1H1.48C3V100bInfT7", "3pa1H1.48C1V100bInfT7",
														"3pa1H1.48C2V100b0.02T7", "3pa1H1.48C2V100b0.04T7", 
														"3pa1H1.48C3V100bInfs0.85T7", "3pa1H1.48C2V100bInfs0.85T7", "3pa1H1.48C1V100bInfs0.85T7",
														"3pa1H1.48C3V100bInfs0.7T7", "3pa1H1.48C2V100bInfs0.7T7", "3pa1H1.48C1V100bInfs0.7T7"), 
										method.legend=c( 'central estimate of HIV infection times\ncentral phylogenetic exclusion criteria',
														 'lower estimate of HIV infection times',
														 'upper estimate of HIV infection times',
														 'central phylogenetic exclusion criteria\nand genetic distance < 2%',
														 'central phylogenetic exclusion criteria\nand genetic distance < 4%',
														 'phylogenetic exclusion criteria\ncoalescent compatibility < 30%\nclade frequency < 80%',
														 'phylogenetic exclusion criteria\ncoalescent compatibility < 10%\nclade frequency < 80%',
														 'phylogenetic exclusion criteria\ncoalescent compatibility < 30%\nclade frequency < 85%',
														 'phylogenetic exclusion criteria\ncoalescent compatibility < 20%\nclade frequency < 85%',
														 'phylogenetic exclusion criteria\ncoalescent compatibility < 10%\nclade frequency < 85%',
														 'phylogenetic exclusion criteria\ncoalescent compatibility < 30%\nclade frequency < 70%',
														 'phylogenetic exclusion criteria\ncoalescent compatibility < 20%\nclade frequency < 70%',
														 'phylogenetic exclusion criteria\ncoalescent compatibility < 10%\nclade frequency < 70%'))										 
	set(tmp, NULL, 'method.legend', tmp[, factor(method.legend, levels=method.legend, labels=method.legend)])
	run.tp	<- merge(run.tp, tmp, by='method.brl')
	setkey(run.tp, factor.legend, t.period.long)
	
	tmp2	<- subset(tmp[1:3, ], select=method.brl)
	tmp2	<- merge(run.tp, tmp2, by='method.brl')
	ggplot(tmp2, aes(x=factor, y=v*100, fill=factor.legend, colour=factor.legend)) + 
			labs(x="stage in HIV infection and care continuum", y="Proportion of transmissions\n(%)") + 
			scale_y_continuous(breaks=seq(0,90,10), minor_breaks=seq(0,90,2), limit=c(0,100*run.tp[, round(max(u95.bs), d=2)]), expand=c(0,0)) + 
			scale_x_discrete(breaks=NULL, limits=run.tp[, intersect(levels(factor), as.character(factor))]) +
			scale_fill_manual(name='from cascade stage', values=run.tp[, unique(factor.color)], guide=FALSE) + 
			scale_colour_manual(name='from cascade stage', values = rep('black',run.tp[, length(unique(factor))]), guide=FALSE) +
			theme_bw() + theme(panel.margin = unit(0.8, "lines"), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.7)) + #coord_flip() +
			geom_bar(stat='identity',binwidth=1, position='dodge') + 
			geom_errorbar(aes(ymin=l95.bs*100, ymax=u95.bs*100), width=0.3, position=position_dodge(width=0.9))	+ 
			facet_grid(method.legend ~ t.period.long, margins=FALSE)
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150316_PropByInfectionTime.pdf'
	ggsave(file=file, w=10, h=10)
	
	tmp2	<- subset(tmp[c(1,6:13), ], select=method.brl)
	tmp2	<- merge(run.tp, tmp2, by='method.brl')
	ggplot(tmp2, aes(x=factor, y=v*100, fill=factor.legend, colour=factor.legend)) + 
			labs(x="stage in HIV infection and care continuum", y="Proportion of transmissions\n(%)") + 
			scale_y_continuous(breaks=seq(0,90,10), minor_breaks=seq(0,90,2), limit=c(0,run.tp[, 100*round(max(u95.bs), d=2)]), expand=c(0,0)) + 
			scale_x_discrete(breaks=NULL, limits=run.tp[, intersect(levels(factor), as.character(factor))]) +
			scale_fill_manual(name='from cascade stage', values=run.tp[, unique(factor.color)], guide=FALSE) + 
			scale_colour_manual(name='from cascade stage', values = rep('black',run.tp[, length(unique(factor))]), guide=FALSE) +
			theme_bw() + theme(panel.margin = unit(0.8, "lines"), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.7)) + #coord_flip() +
			geom_bar(stat='identity',binwidth=1, position='dodge') + 
			geom_errorbar(aes(ymin=l95.bs*100, ymax=u95.bs*100), width=0.3, position=position_dodge(width=0.9))	+ 
			facet_grid(method.legend ~ t.period.long, margins=FALSE)
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150316_PropByExclusionCriteria.pdf'
	ggsave(file=file, w=10, h=25)	
}
######################################################################################
project.athena.Fisheretal.composition.prop.phylolikelihoods<- function()	
{
	#	t.period
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-5','2008-1','2009-7'), labels=c('96/07','\n\n06/05','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-4','2007-12','2009-6','2010-12'), labels=c('06/04','07/12','09/06','10/12'))])		
	#	load risk tables
	indir				<- paste(DATA,"fisheretal_150312",sep='/')
	outdir				<- paste(DATA,"fisheretal_150312",sep='/')			
	infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	insignat			<- "Wed_Dec_18_11:37:00_2013"	
	
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))
	runs.risk			<- subset(runs.risk, !grepl('ARTstarted',method.risk) & !grepl('GroupsUDA',method.risk))
	method.DENOM		<- 'SEQ'		
	method.WEIGHT		<- ''
	method.RISK			<- c("m2Awmx.nophyloscore.tp","m2Awmx.noscore.tp","m2Awmx.wtn.tp","m2Awmx.tp")
	method.DATING		<- 'sasky'	
	stat.select			<- 'P.raw.e0cp'
	outfile				<- infile
	method.BRL			<- '3pa1H1.48C2V100bInfT7'
	
	run.tp				<- subset(runs.risk, 	(is.na(t.period) | t.period!=0) & method.denom==method.DENOM & method.nodectime=='any' & 
												method.brl==method.BRL & method.dating==method.DATING & 
												(grepl('P.',stat,fixed=1) | stat=='P' | grepl('N.',stat,fixed=1) | stat=='N') )
	set(run.tp, NULL, 'method.risk', run.tp[, substr(method.risk, 1, nchar(method.risk)-1)])
	run.tp				<- subset(run.tp, method.risk%in%method.RISK)
	run.tp				<- subset(run.tp, !is.na(v) & stat%in%stat.select )
	#	set up factor legends
	factors				<- project.athena.Fisheretal.sensitivity.factor.legend("m2Awmx")	
	#	prepare
	run.tp[, m50.bs:=NULL]
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])
	set(run.tp, NULL, 'factor', run.tp[, factor(factor, levels=factors[, levels(factor)])])	
	#	merge t.period
	run.tp				<- merge(run.tp, tperiod.info, by='t.period')
	run.tp[, t.period.long:= gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))]
	set(run.tp, NULL, 't.period.long', run.tp[,factor(t.period.long, levels= tperiod.info[, gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))])])	
	#	merge factors
	run.tp				<- merge(run.tp, subset(factors, select=factors[, which(colnames(factors)!='method.risk')]), by='factor')	
	#	merge methods legend
	tmp		<- data.table(	method.risk=c("m2Awmx.wtn.tp","m2Awmx.tp","m2Awmx.noscore.tp","m2Awmx.nophyloscore.tp"), 
							method.legend=c(	'with phylogenetic transmission\nprobability per interval',
											'with phylogenetic transmission\nprobability per probable transmitter',
											'every interval\nequally likely',
											'every probable transmitter\n equally likely'))
	set(tmp, NULL, 'method.legend', tmp[, factor(method.legend, levels=method.legend, labels=method.legend)])
	run.tp	<- merge(run.tp, tmp, by='method.risk')

	setkey(run.tp, factor.legend, t.period.long)
	ggplot(run.tp, aes(x=factor, y=v*100, fill=factor.legend, colour=factor.legend)) + 
			labs(x="stage in HIV infection and care continuum", y="Proportion of transmissions\n(%)") + 
			scale_y_continuous(breaks=seq(0,90,10), minor_breaks=seq(0,90,2), limit=c(0,100*run.tp[, round(max(u95.bs), d=2)]), expand=c(0,0)) + 
			scale_x_discrete(breaks=NULL, limits=run.tp[, intersect(levels(factor), as.character(factor))]) +
			scale_fill_manual(name='from cascade stage', values=run.tp[, unique(factor.color)], guide=FALSE) + 
			scale_colour_manual(name='from cascade stage', values = rep('black',run.tp[, length(unique(factor))]), guide=FALSE) +
			theme_bw() + theme(panel.margin = unit(0.8, "lines"), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.7)) + #coord_flip() +
			geom_bar(stat='identity',binwidth=1, position='dodge') + 
			geom_errorbar(aes(ymin=l95.bs*100, ymax=u95.bs*100), width=0.3, position=position_dodge(width=0.9))	+ 
			facet_grid(method.legend ~ t.period.long, margins=FALSE)
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150312_PropByLklMethod.pdf'
	ggsave(file=file, w=10, h=13)	
}
######################################################################################
project.athena.Fisheretal.composition.prop.samplingcensoring<- function()	
{
	#	t.period
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-5','2008-1','2009-7'), labels=c('96/07','\n\n06/05','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-4','2007-12','2009-6','2010-12'), labels=c('06/04','07/12','09/06','10/12'))])		
	#	load risk tables
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))	
	method.DENOM		<- 'SEQ'		
	method.WEIGHT		<- ''
	method.DATING		<- 'sasky'	
	stat.select			<- c(	'P.raw','P.raw.e0','P.raw.e0cp'	)
	outfile				<- infile
	method.BRL			<- '3pa1H1.48C2V100bInfT7'
	method.RISK			<- "m2Awmx.wtn.tp"
	run.tp				<- subset(runs.risk, (is.na(t.period) | t.period!=0) & method.denom==method.DENOM & method.nodectime=='any' & method.brl==method.BRL & method.dating==method.DATING & grepl(method.RISK,method.risk)  & (grepl('P.',stat,fixed=1) | stat=='P' | grepl('N.',stat,fixed=1) | stat=='N') )
	run.tp				<- subset(run.tp, !grepl('ARTstarted',method.risk) & !grepl('GroupsUDA',method.risk))
	if(method.WEIGHT=='')
		run.tp			<- subset(run.tp, !grepl('wstar',method.risk) & !grepl('now',method.risk))
	if(method.WEIGHT!='')
		run.tp			<- subset(run.tp, grepl(method.WEIGHT,method.risk) )
	if(is.null(tperiod.info))
	{
		tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
		set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
		set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )		
	}		
	run.tp[, m50.bs:=NULL]
	setkey(run.tp, factor)
	run.tp[, t.period:=run.tp[, substr(factor, nchar(factor), nchar(factor))]]
	set(run.tp, NULL, 'factor', run.tp[, substr(factor, 1, nchar(factor)-2)])
	set(run.tp, NULL, 'factor', run.tp[, factor(factor, levels=factors[, levels(factor)])])
	factors				<- project.athena.Fisheretal.sensitivity.factor.legend(method.RISK)
	run.tp	<- merge(run.tp, subset(factors, select=factors[, which(colnames(factors)!='method.risk')]), by='factor')	
	run.tp	<- merge(run.tp, tperiod.info, by='t.period')
	run.tp[, t.period.long:= gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))]
	set(run.tp, NULL, 't.period.long', run.tp[,factor(t.period.long, levels= tperiod.info[, gsub('\n','',paste(t.period.min, '-', t.period.max,sep=''))])])
	run.tp	<- subset(run.tp, !is.na(v) & stat%in%c(	'P.raw','P.raw.e0','P.raw.e0cp'	) )
	setkey(run.tp, factor.legend, t.period.long)
	tmp		<- data.table(stat=c('P.raw.e0cp','P.raw.e0','P.raw'), stat.legend=c('with adjustment for\ncensoring and sampling biases','with adjustment for\nsampling biases','no adjustment for\nsampling and censoring biases'))
	set(tmp, NULL, 'stat.legend', tmp[, factor(stat.legend, levels=stat.legend, labels=stat.legend)])
	run.tp	<- merge(run.tp, tmp, by='stat')
	ggplot(run.tp, aes(x=factor, y=v*100, fill=factor.legend, colour=factor.legend)) + 
			labs(x="stage in HIV infection and care continuum", y="Proportion of transmissions\n(%)") + 
			scale_y_continuous(breaks=seq(0,90,10), minor_breaks=seq(0,90,2), limit=c(0,1+100*run.tp[, round(max(u95.bs), d=2)]), expand=c(0,0)) + 
			scale_x_discrete(breaks=NULL, limits=run.tp[, intersect(levels(factor), as.character(factor))]) +
			scale_fill_manual(name='from cascade stage', values=run.tp[, unique(factor.color)], guide=FALSE) + 
			scale_colour_manual(name='from cascade stage', values = rep('black',run.tp[, length(unique(factor))]), guide=FALSE) +
			theme_bw() + theme(panel.margin = unit(0.8, "lines"), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey70", size=0.7)) + #coord_flip() +
			geom_bar(stat='identity',binwidth=1, position='dodge') + 
			geom_errorbar(aes(ymin=l95.bs*100, ymax=u95.bs*100), width=0.3, position=position_dodge(width=0.9))	+ 
			facet_grid(stat.legend ~ t.period.long, margins=FALSE)
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150312_PropBySamplingCens.pdf'
	ggsave(file=file, w=10, h=10)	
}
######################################################################################
project.athena.Fisheretal.composition.transprob<- function()	
{
	#	t.period
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.503, 2006.408, 2008.057, 2009.512), t.period.max = c(2006.308, 2007.957, 2009.49, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-5','2008-1','2009-7'), labels=c('96/07','\n\n06/05','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-4','2007-12','2009-6','2010-12'), labels=c('06/04','07/12','09/06','10/12'))])		
	#	load YX
	method.RISK		<- 'm2Awmx.wtn'
	indir			<- paste(DATA,"fisheretal_150319",sep='/')
	infile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Yscore3pa1H1.48C2V100bInfT7_denomSEQ_model2_m2Awmx.wtn.tp1.R'
	#outfile			<- 'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXm2_SEQ_3la2H1C2V100_lkl.pdf'	
	file			<- paste(indir, '/', infile, sep='')
	load(file)
	YX				<- copy(ans$YXf)
	#set(YX, NULL, 'score.Y', YX[, score.Y*w.tn])
	YX				<- subset(YX, select=c(t.period, t, t.Patient, Patient, score.Y, CD4b, CD4c ))	
	setnames( YX, 'CD4c', 'factor')	
	set(YX, YX[, which(as.numeric(t.period)>4)],'t.period','4')
	#	get trm.p for data
	file			<- paste(indir,'ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Wed_Dec_18_11:37:00_2013_method.trmp.BS0.Rdata',sep='/')
	load(file)	
	trm.p			<- subset(trm.p, BS==0 & grepl(method.RISK,method.risk) & grepl('3pa1H1.48C2V100bInfT7',method.brl))
	trm.p[, t.period:= trm.p[, substring(factor, nchar(factor))]]
	set(trm.p, NULL, 'factor', trm.p[, substr(factor, 1, nchar(factor)-2)])
	#	get denominator for every patient
	tmp		<- trm.p[, list( ys.e0cp=sum(ys.e0cp)+sum(ys), ys=sum(ys), ysr= sum(ys)/(sum(ys.e0cp)+sum(ys)) ) , by=c('Patient','t.period')]
	#	for some patients have quite a lot of data ..
	#ggplot(tmp, aes(x=ysr)) + geom_histogram()	
	#	merge with YXf and plot probs
	YX		<- merge(YX, tmp, by=c('t.period','Patient'))
	YX[, tp:= score.Y/ys.e0cp]
	#	check
	#merge(YX, YX[, list(yscheck=sum(score.Y)), by=c('t.period','Patient')], by=c('t.period','Patient'))
	#	merge tperiod.long
	YX		<- merge(YX, tperiod.info, by='t.period')
	YX[, t.period.long:= paste(t.period.min, '-', t.period.max,sep='')]
	#	merge legends
	factors	<- project.athena.Fisheretal.sensitivity.factor.legend(method.RISK)
	YX		<- merge(YX, factors, by='factor')
	setkey(YX, factor)
	#	overall
	ggplot(YX, aes(x=tp*100)) + geom_histogram(binwidth=0.5) + theme_bw() +
			scale_x_continuous(breaks=seq(0,100,10), minor_breaks=seq(0,100,5), expand=c(0,0)) +
			labs(x='phylogenetically derived transmission probability\nper observed transmission interval\n(%)')
			#theme(panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major=element_line(colour="grey60", size=0.4))
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150311_TrProbHist.pdf'
	ggsave(file=file, w=10, h=3)
	#	overall stats
	YX[, c(mean(tp), quantile(tp, prob=c(0.025, 0.25, 0.5, 0.75, 0.975, 1)))]
	#	1.919032e-02 2.167467e-09 9.322279e-04 7.260539e-03 2.254388e-02 1.080519e-01 6.290796e-01 
	YX[, c(mean(tp), quantile(tp, prob=c(0.025, 0.25, 0.5, 0.75, 0.975, 1))), by='t.period']
	#
	#	boxplot of transmission prob (observed only)
	#
	setkey(YX, factor.legend)
	ggplot(YX, aes(x=factor.legend, y=tp*100, colour=factor.legend, fill=factor.legend))  + 
			labs(x='', y='phylogenetic\ntransmission probability\nduring observed intervals\n(%)') +
			geom_boxplot(outlier.shape = NA) +
			scale_colour_manual(name='', values=YX[, unique(factor.color)], guide=FALSE) +
			scale_fill_manual(name='', values=YX[, unique(factor.color)], guide=FALSE) +
			scale_y_continuous(breaks=seq(0, 50, 4), minor_breaks=seq(0, 50, 1)) +
			scale_x_discrete(breaks=NULL) +	
			#coord_trans(ytrans="log10", limy=c(0.001, 0.3)) +
			coord_cartesian(ylim=c(0, 14)) +
			stat_summary(geom= "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +		
			theme_bw() + theme(axis.title=element_text(size=12), axis.text=element_text(size=12), legend.key.size=unit(11,'mm'), panel.grid.minor=element_line(colour="grey60", size=0.2), panel.grid.major.x=element_blank(), panel.grid.major=element_line(colour="grey60", size=0.2))	
	file			<- paste(indir, '/', gsub('tp1.R','tprob.pdf',infile), sep='')
	ggsave(file=file, w=5.5, h=3)	
	
	YX[, mean(tp), by='factor']
}

######################################################################################
project.athena.Fisheretal.composition.problintervals.table.generate<- function(dft, file)
{
	dft[, t.period:= dft[, as.character(factor)]]	
	set(dft, NULL, 't.period', dft[, substring(t.period, nchar(t.period))] )
	set(dft, NULL, 'factor', dft[, as.character(factor)] )
	set(dft, NULL, 'factor', dft[, substr(factor, 1, nchar(factor)-2)] )
	set(dft, NULL, 'p', dft[, round(100*p,d=1)])
	#	proportions
	dftp				<- dcast.data.table(subset(dft, stat=='YX'), factor+method.brl~t.period, value.var='p')	
	#	total per recipient MSM
	tmp					<- subset(dft, stat=='YX')[, list(n=sum(n)), by=c('method.brl','t.period')]
	tmp2				<- unique(subset(dft, stat=='nRecLkl', select=c(method.brl, t.period, n)))
	tmp					<- merge(tmp, tmp2, by=c('method.brl','t.period'))
	set(tmp, NULL, 'n', tmp[, round(n.x/n.y, d=3)])
	tmp					<- dcast.data.table(tmp, method.brl~t.period, value.var='n')
	tmp[, factor:='']
	#	put everything together
	dftp				<- rbind(dftp, tmp, use.names=TRUE)
	tmp					<- c("","UA","UAna","U","DA","Dt.NA","Dtg500","Dtl500","Dtl350","ART.NotYetFirstSu","ART.vlNA","ART.suA.N","ART.suA.Y1","ART.suA.Y2","Lost")
	set( dftp, NULL, 'factor', dftp[, factor(factor, levels=tmp, labels=tmp)] )
	setkey(dftp, method.brl, factor)
	ans	<- rbind(	dcast.data.table( subset(dftp, select=c('method.brl', 'factor', '1')), factor~method.brl, value.var='1'),
			dcast.data.table( subset(dftp, select=c('method.brl', 'factor', '2')), factor~method.brl, value.var='2'),	
			dcast.data.table( subset(dftp, select=c('method.brl', 'factor', '3')), factor~method.brl, value.var='3'),
			dcast.data.table( subset(dftp, select=c('method.brl', 'factor', '4')), factor~method.brl, value.var='4')	)
	ans	<- subset(ans, select=c(1,3,4,2))	
	write.csv(ans, file=file, eol='\r\n')
}
######################################################################################
project.athena.Fisheretal.composition.problintervals.table<- function()
{
	require(data.table)
	require(ape)
	require(RColorBrewer)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_150319",sep='/')	
	outdir					<- paste(DATA,"fisheretal_150319",sep='/')		
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.table.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))	
	dft.all				<- subset(runs.table, !grepl('ARTstarted|GroupsUDA',method.risk) & stat%in%c('YX','nRecLkl'))
	#
	#	first table
	#
	dft					<- subset(dft.all, grepl('bInfT7', method.brl) & grepl('H1.48', method.brl))
	file				<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150304_PrIntervalsByCoalComp.csv'
	project.athena.Fisheretal.composition.problintervals.table.generate(dft, file)
	#
	#	second table
	#
	dft					<- subset(dft.all, grepl('C2', method.brl) & grepl('bInf', method.brl) & grepl('H1.48', method.brl))
	file				<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150304_PrIntervalsByCladeFreq.csv'
	project.athena.Fisheretal.composition.problintervals.table.generate(dft, file)
}
######################################################################################
project.athena.Fisheretal.composition.acute<- function()
{
	tmp<- df.all
	tmp<- df.all.allmsm
	setkey(tmp, Patient)
	tmp<- unique(tmp)
	tmp<- subset(tmp, AnyPos_T1>=1996.6 & AnyPos_T1<2011 & Trm%in%c('MSM','BI'), select=c(Patient, isAcute, AnyPos_T1))
	set(tmp, tmp[, which(is.na(isAcute))],'isAcute', 'Missing')
	
	ggplot(tmp, aes(x=AnyPos_T1, fill=isAcute, colour=isAcute)) + geom_bar(aes(y=..count../sum(..count..)), binwidth=1, position='dodge')
	ggplot(tmp, aes(x=AnyPos_T1, fill=isAcute, colour=isAcute)) + geom_bar(binwidth=0.5, position='dodge')
	ggplot(tmp, aes(x=AnyPos_T1, fill=isAcute, colour=isAcute)) + geom_bar(binwidth=1, position='fill')
	
	#tmp<- subset(YX, select=c(Patient, t.Patient, w.i, score.Y))
	#setkey(tmp, Patient, t.Patient)
	#tmp<- unique(tmp)
	
	#
	tmp	<- df.all.allmsm
	setkey(tmp, Patient)
	tmp	<- unique(tmp)
	tmp	<- subset(tmp, AnyPos_T1>=1996.5 & AnyPos_T1<2011 & Trm%in%c('MSM','BI'), select=c(Patient, isAcute, Acute_Spec, AnyPos_T1, Trm))
	#
	subset(tmp, isAcute=='Yes')[, table(Acute_Spec, useNA='if')/length(Acute_Spec)]
	#
	set(tmp, tmp[, which(is.na(Acute_Spec))],'Acute_Spec', 'Missing')
	set(tmp, tmp[, which(isAcute=='No')],'Acute_Spec', 'Recent HIV infection not indicated')
	set(tmp, tmp[, which(isAcute=='Yes')],'Acute_Spec', 'Confirmed recent HIV infection')
	set(tmp, tmp[, which(Acute_Spec%in%c('CLIN','LAB','NEGTEST'))],'Acute_Spec','Confirmed recent HIV infection')
	set(tmp, tmp[, which(Acute_Spec=='SYM')],'Acute_Spec','Recent HIV infection not indicated')
	tmp	<- merge(tmp, data.table(Acute_Spec=c('Confirmed recent HIV infection', 'Recent HIV infection not indicated', 'Missing' ), Acute.col=c("#FF7F00","#FDBF6F",'grey70')), by='Acute_Spec')
	set(tmp, NULL, 'Acute_Spec', tmp[,factor(Acute_Spec, levels=c('Confirmed recent HIV infection', 'Recent HIV infection not indicated', 'Missing' ), labels=c('Confirmed recent HIV infection', 'Recent HIV infection not indicated', 'Missing' ))])
	tmp[, t.period:= tmp[, cut(AnyPos_T1, breaks=c(-Inf, 2006.5, 2008, 2009.5, 2011))]]		
	#
	setkey(tmp, Acute_Spec)
	ggplot(tmp, aes(x=AnyPos_T1, fill=Acute_Spec)) + 
			geom_bar(binwidth=0.25, position='fill',alpha=0.8) +
			scale_y_continuous(breaks=seq(0,1,0.2),labels=seq(0,1,0.2)*100) +
			scale_x_continuous(breaks=seq(1997,2014,2), minor_breaks=NULL, expand=c(0,0)) +
			scale_fill_manual(name='', values=tmp[, unique(Acute.col)]) +
			labs(x='', y='Infection status at diagnosis\n( % )') + 
			facet_grid(.~t.period, scales='free_x', space='free_x') +					
			theme_bw() +
			theme(legend.position='bottom', legend.text=element_text(size=14), legend.key.size=unit(7,'mm'), axis.title=element_text(size=14), axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), strip.background = element_blank(), strip.text = element_blank(), panel.grid.major.y = element_line(colour="black", size=0.4), panel.grid.minor.y = element_line(colour="black", size=0.4), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), plot.margin=unit(c(0,2,0,0),"cm")) +	
			guides(fill=guide_legend(ncol=2))
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150223_AcuteSpec_Time.pdf'
	ggsave(file=file, w=10, h=4)
	#
	tmp[, table(t.period, Acute_Spec)]
	
}

######################################################################################
project.athena.Fisheretal.composition.totalmissing<- function()	
{
	require(data.table)
	require(ape)
	require(RColorBrewer)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_150319",sep='/')	
	outdir					<- paste(DATA,"fisheretal_150319",sep='/')	
	#
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.table.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))	
	
	
	method.DENOM	<- 'SEQ'
	method.BRL		<- '3pa1H1.48C2V100bInfT7'
	method.RISK		<- 'm2Awmx.wtn.tp'
	method.WEIGHT	<- ''	
	
	df		<- subset(runs.table, method.denom==method.DENOM & method.BRL==method.brl & grepl(method.RISK,method.risk), c(stat, method.risk, factor, n, p, l95.bs, u95.bs))
	tmp		<- subset(df, stat%in%c('Sx.e0cp'))	
	tmp[, list(n=sum(n), py=sum(n/8), l95.bs=sum(l95.bs), u95.bs=sum(u95.bs), pyl95.bs=sum(l95.bs/8), pyu95.bs=sum(u95.bs/8))]
}
######################################################################################
project.athena.Fisheretal.composition.contact<- function()
{
	t.endctime		<- 2013.3
	contact.grace	<- 1.5
	t.period		<- 0.125
	
	contact		<- subset(df.all.allmsm, Trm%in%c('MSM','BI'), select=c(Patient, AnyPos_T1, DateLastContact, DateDied, ReasonStopRegistration))
	setkey(contact, Patient)
	contact		<- unique(contact)
	#	set last contact for NA last contact
	set(contact, contact[, which(is.na(DateDied) | DateDied>t.endctime)], 'DateDied', t.endctime)
	tmp			<- contact[, which(DateLastContact>=DateDied)]
	set(contact, tmp, 'DateLastContact', contact[tmp, DateDied])
	#	set DateDied for those that moved
	tmp			<- contact[, which(!is.na(DateLastContact) & ReasonStopRegistration=='Moved')]
	set(contact, tmp, 'DateDied', contact[tmp, DateLastContact])
	#
	set(df.treatment.allmsm, NULL, 'StopTime', hivc.db.Date2numeric(df.treatment.allmsm[,StopTime]))
	set(df.treatment.allmsm, NULL, 'StartTime', hivc.db.Date2numeric(df.treatment.allmsm[,StartTime]))
	df.treat	<- melt( df.treatment.allmsm, id.vars=c('Patient'), measure.vars=c('StartTime','StopTime'), value.name='AnyT_T' )
	setkey(df.treat, Patient, AnyT_T)
	df.treat	<- unique(df.treat)
	df.treat	<- merge( df.treat, subset(contact, select=c(Patient, DateDied)),  by='Patient', all.x=1 )
	#	keep only last AnyT_T before death. If not, adjusting DateLastContact precisely removes all not in contact periods
	df.treat	<- subset(df.treat, AnyT_T<DateDied, c(Patient, AnyT_T))
	tmp			<- merge( unique(subset(contact, select=Patient)), subset(df.treat, select=c(Patient, AnyT_T)), by='Patient')
	tmp			<- tmp[, list(AnyT_TL=max(AnyT_T)), by='Patient']	
	contact		<- merge(contact, tmp, by='Patient', all.x=1)	
	tmp			<- merge(unique(subset(contact, select=Patient)), subset(df.viro.allmsm, select=c(Patient, PosRNA)), by='Patient')
	set(tmp, NULL, 'PosRNA', hivc.db.Date2numeric(tmp[,PosRNA]))
	tmp			<- tmp[, list(PosRNA_TL=max(PosRNA)), by='Patient']	
	contact		<- merge(contact, tmp, by='Patient', all.x=1)
	tmp			<- merge(unique(subset(contact, select=Patient)), subset(df.immu.allmsm, select=c(Patient, PosCD4)), by='Patient')
	set(tmp, NULL, 'PosCD4', hivc.db.Date2numeric(tmp[,PosCD4]))
	tmp			<- tmp[, list(PosCD4_TL=max(PosCD4)), by='Patient']	
	contact		<- merge(contact, tmp, by='Patient', all.x=1)		
	tmp			<- contact[, which(DateLastContact<PosRNA_TL & PosRNA_TL<t.endctime)]
	if(length(tmp))
	{
		cat(paste('\nsetting relaxed DateLastContact because PosRNA_TL for n=',length(tmp)))
		set(contact, tmp, 'DateLastContact', contact[tmp,PosRNA_TL])
	}
	tmp			<- contact[, which(DateLastContact<PosCD4_TL & PosCD4_TL<t.endctime)]
	if(length(tmp))
	{
		cat(paste('\nsetting relaxed DateLastContact because PosCD4_TL for n=',length(tmp)))
		set(contact, tmp, 'DateLastContact', contact[tmp,PosCD4_TL])
	}
	tmp			<- contact[, which(DateLastContact<AnyT_TL & AnyT_TL<t.endctime)]	
	if(length(tmp))
	{
		cat(paste('\nsetting relaxed DateLastContact because AnyT_TL for n=',length(tmp)))
		set(contact, tmp, 'DateLastContact', contact[tmp,AnyT_TL])
	}	
	#	allow for grace at end
	#tmp			<- contact[, which(DateLastContact<t.endctime & DateLastContact+contact.grace>=t.endctime)]
	#if(length(tmp))
	#{
	#	cat(paste('\nsetting relaxed DateLastContact because DateLastContact+contact.grace<t.endctime for n=',length(tmp)))
	#	set(contact, tmp, 'DateLastContact', contact[tmp, t.endctime])
	#}
	#	discretize
	set(contact, NULL, 'DateLastContact', contact[, floor(DateLastContact) + ceiling( (DateLastContact%%1)*100 %/% (t.period*100) ) * t.period] )
	set(contact, NULL, 'DateDied', contact[, ceiling(DateDied) + round( (DateDied%%1)*100 %/% (t.period*100) ) * t.period] )
	set(contact, NULL, 'AnyPos_T1', contact[, floor(AnyPos_T1) + round( (AnyPos_T1%%1)*100 %/% (t.period*100) ) * t.period] )
	#	set up patient timelines
	X.incare	<- contact[, list(t=seq(AnyPos_T1, DateDied-t.period, by=t.period)), by='Patient']
	tmp			<- subset(contact, DateLastContact<DateDied)	
	tmp			<- tmp[, list(t= seq(DateLastContact, DateDied-t.period, by=t.period), contact='No'),by='Patient']	
	setnames(tmp, 'Patient','t.Patient')
	setnames(X.incare, 'Patient','t.Patient')	
	X.incare	<- merge(X.incare, tmp, by=c('t.Patient','t'), all.x=1)
	set(X.incare, X.incare[,which(is.na(contact))],'contact','Yes')
	#set(X.incare, X.incare[, which(contact=='No' & t+contact.grace>t.endctime)], 'contact', 'Yes')
	#
	#	for each time t, check if there is at least 1 VL or CD4 measurement or Treatment visit in the last year/2 or next year/2
	#
	set(df.immu.allmsm, NULL, 'PosCD4', hivc.db.Date2numeric(df.immu.allmsm[,PosCD4]))
	set(df.viro.allmsm, NULL, 'PosRNA', hivc.db.Date2numeric(df.viro.allmsm[,PosRNA]))
	tmp			<- X.incare[, {
				cntct					<- rep(0, length(t))
				endgrace				<- t+contact.grace-t.endctime
				endgrace[endgrace<0]	<- 0	
				endgrace[endgrace>0]	<- 0	
				z						<- df.immu.allmsm$PosCD4[ which(df.immu.allmsm$Patient==t.Patient) ]
				if(length(z))
				{
					tmp			<- sapply(seq_along(t), 	function(i) min(abs(t[i]-z))<=(contact.grace/2+endgrace[i])		)
					cntct[tmp]	<- 1
				}									
				z		<- df.viro.allmsm$PosRNA[ which(df.viro.allmsm$Patient==t.Patient) ]
				if(!all(cntct==1) & length(z))
				{
					tmp			<- sapply(seq_along(t), 	function(i) min(abs(t[i]-z))<=(contact.grace/2+endgrace[i])		)
					cntct[tmp]	<- 1
				}
				z						<- df.treat$AnyT_T[ which(df.treat$Patient==t.Patient) ]
				if(length(z))
				{
					tmp			<- sapply(seq_along(t), 	function(i) min(abs(t[i]-z))<=(contact.grace/2+endgrace[i])		)
					cntct[tmp]	<- 1
				}																		
				list(t=t, incontact=cntct)
			}, by=c('t.Patient')]	
	X.incare	<- merge(X.incare, tmp, by=c('t.Patient','t'))	
	
	tmp			<- subset(X.incare, t>=1996 & t<2013.125)[, list(NOCON= mean(contact=='No' | incontact==0)), by='t']
	#tmp			<- subset(X.incare, t>=1996 & t<2013.125)[, list(NOCON= mean(contact=='No')), by='t']
	ggplot(tmp, aes(x=t, ymax=100*NOCON, ymin=0)) + geom_ribbon() + 
			scale_x_continuous(breaks=seq(1980, 2020, 5), minor_breaks=seq(1980, 2020, 1)) + 
			scale_y_continuous(expand=c(0,0), limits=c(0,25), breaks=seq(0,100,5)) +
			coord_trans(limx=c(1996.5, 2013.125)) +
			theme(panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.4), axis.text.x=element_text(angle=0, vjust=0, hjust=0)) +
			labs(x='', y='MSM with no contact to care\nfor at least 18 months\n(%)') +
			theme_bw()
	subset(tmp, t==2012 | t==2003)
	
}
######################################################################################
project.athena.Fisheretal.composition.age<- function()
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
project.athena.Fisheretal.composition.clock<- function()
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
project.athena.Fisheretal.composition.coalvsbrl<- function()
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
project.athena.Fisheretal.composition.coal<- function()
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
project.athena.Fisheretal.composition.seqdata<- function()
{
	load('~/duke/2013_HIV_NL/ATHENA_2013/data/derived/ATHENA_2013_03_Sequences.R')
	#	total 8754 with date, 8380 longer than 21
	verbose			<- 1
	resume			<- 1 
	indircov		<- paste(DATA,"derived",sep='/')
	infilecov		<- "ATHENA_2013_03_AllSeqPatientCovariates"
	indir			<- paste(DATA,"tmp",sep='/')	
	infile			<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500"			
	insignat		<- "Wed_Dec_18_11/37/00_2013"
	#	total used in analysis
	argv			<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	nsh.clu.pre		<- hivc.prog.get.clustering.precompute()	
	nshs			<- subset(nsh.clu.pre$df.seqinfo, select=c(FASTASampleCode, Patient, PosSeqT, AnyPos_T1, AnyT_T1, Trm, Sex))
	#	total in data set
	load('~/duke/2013_HIV_NL/ATHENA_2013/data/derived/ATHENA_2013_03_-DR-RC+LANL_Sequences_Fri_Nov_01_16:07:23_2013.R')	
	dfs		<- data.table( FASTASampleCode=rownames(seq.PROT.RT), LEN=seq.length(seq.PROT.RT) )
	dfs[, FRG:= grepl('TN|PROT',FASTASampleCode)]
	#	working out differences..
	all( nshs[, FASTASampleCode]%in%dfs[, FASTASampleCode] )
	nshs[, DATASET:='nsh']
	dfs		<- merge(dfs, nshs, by='FASTASampleCode', all.x=1)	
	dfsNL	<- subset(dfs, DATASET=='nsh' & !FRG)
	
	set(dfsNL, NULL, 'PosSeqT', hivc.db.Date2numeric(dfsNL[,PosSeqT]))
	set(dfsNL, NULL, 'AnyPos_T1', hivc.db.Date2numeric(dfsNL[,AnyPos_T1]))
	set(dfsNL, NULL, 'AnyT_T1', hivc.db.Date2numeric(dfsNL[,AnyT_T1]))	
	
	dfsNL[, dts:= PosSeqT-AnyPos_T1]
	dfsNL[, dtt:= 'After ART start']
	set(dfsNL, dfsNL[, which(is.na(AnyT_T1) | PosSeqT<AnyT_T1)], 'dtt','Before ART start') 
	
	dfpNL	<- subset(dfsNL, !is.na(PosSeqT))[, list(SEQN=length(FASTASampleCode), PosSeqT1=min(PosSeqT), AnyPos_T1=min(AnyPos_T1)), by='Patient']
	subset(dfsNL, !is.na(dts))[, quantile(LEN, p=c(0.025, 0.05, 0.95, 0.975))]
	dfpNL[, summary(PosSeqT1-AnyPos_T1)]
	dfpNL[, quantile(PosSeqT1-AnyPos_T1, p=c(0.025, 0.05, 0.95, 0.975))]
	
	setkey(dfsNL, Patient)
	unique(subset(dfsNL, !is.na(PosSeqT)))[, table(Trm, Sex, useNA='if')]
	#	plot numbers by year
	dfsNL[, PosSeqTc:=floor(PosSeqT)]
	ggplot(dfsNL, aes(x=PosSeqTc, fill=dtt)) + geom_bar(binwidth=1) + theme_bw() +
			labs(x='Sampling year', y='ATHENA sequences') +
			scale_fill_brewer(name='sampling time', palette='Set1') +
			scale_x_continuous(breaks=seq(1980,2020,5), minor_breaks=seq(1980,2020,1)) +
			scale_y_continuous(breaks=seq(0,2000,200)) +
			geom_vline(xintercept=2011, size=1) +
			geom_vline(xintercept=2013.3, size=1, linetype=2) +
			theme(legend.justification=c(0,1), legend.position=c(0,1), panel.grid.major = element_line(colour='grey70'), panel.grid.minor = element_line(colour='grey80'))
	file	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2014/MSMtransmission_ATHENA1303/150226_SeqNo.pdf'
	ggsave(file=file, w=8, h=3)
	
}
######################################################################################
project.athena.Fisheretal.composition.seqcoverage<- function()
{
	require(Hmisc)
	pts			<- unique(subset(X.seq, select=t.Patient))
	pt			<- unique(subset(X.msm, select=t.Patient))
	prob.pt.seq	<- nrow(pts)/nrow(pt)
	
	#	calculate sequence coverage over time
	setnames(pts, 't.Patient','Patient')
	setnames(pt, 't.Patient','Patient')
	tmp	<- subset(df.all.allmsm, select=c(Patient, AnyPos_T1, DateDied, PosSeqT))
	setkey(tmp, Patient)
	tmp	<- unique(tmp)
	pts	<- merge(pts, tmp, by='Patient')
	pt	<- merge(pt, tmp, by='Patient')
	
	dfd	<- data.table(t=seq(1996.5, 2011, 0.125))
	dfd	<- dfd[, list( 	NPT=subset(pt, AnyPos_T1>=t & (is.na(DateDied) | DateDied<t))[, length(unique(Patient))],
					NPTS=subset(pts, AnyPos_T1>=t & (is.na(DateDied) | DateDied<t))[, length(unique(Patient))]		)
			, by='t']
	dfd[, SIFD:=NPTS/NPT]	
	dfr	<- YX[, list(NR=length(unique(Patient))), by='t.Patient']	
	ggplot(dfd, aes(x=t, y=SIFD)) + geom_line()
	#	because of recent decline in coverage, the average prob is quite generous
	
	#	load a few YX to see how many recipients retained
	files	<- list( 	data.table(b=Inf, c=0.3, FILE='ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXSEQ3pa1H1.35C3V100bInfT7.R'),
			data.table(b=0.08, c=0.3, FILE='ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXSEQ3pa1H1.35C3V100T7.R'),	
			data.table(b=Inf, c=0.4, FILE='ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXSEQ3pa1H1.35C4V100bInfT7.R'),
			data.table(b=Inf, c=0.5, FILE='ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXSEQ3pa1H1.35V100bInfT7.R')	
	)
	files		<- do.call('rbind',files)
	
	#X.msm.prs	<- unique(subset(X.msm, select=c(Patient, t.Patient)))	
	#load( paste(outdir,'/',files[1,FILE],sep='') )	
	#tmp	<- merge( X.msm.prs, unique(subset(YX, select=Patient)), by='Patient' )
	#tmp[, length(unique(t.Patient))]
	
	alpha	<- 0.01
	YX.info	<- files[, {
				load( paste(outdir,'/',FILE,sep='') )
				list(	nRexp= nrow(ri.SEQ)*prob.pt.seq, nRexp.ql=qbinom(alpha/2, nrow(ri.SEQ), prob.pt.seq), nRexp.qu=qbinom(1-alpha/2, nrow(ri.SEQ), prob.pt.seq), 
						nR= YX[, length(unique(Patient))], nPT=YX[, length(unique(t.Patient))], nI=nrow(YX))
			}, by=c('b','c')]
}
######################################################################################
project.athena.Fisheretal.composition.vrancken.brl<- function()
{
	file				<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_YXSEQ3pa1H1.35C3V100T7.R"
	load(file)					#YX
	infile.trm.model	<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_data/TchainBelgium_set7_pol_GA11emodel_nA_INFO.R"
	load(infile.trm.model)		#expect "trm.pol.GA" "trm.pol.nA"
	
	df.vr	<- subset(YX, select=c(Patient, t.Patient, t, brl, score.Y, telapsed, t.period))
	setkey( df.vr, brl )
	
	
	setnames(df.vr,'telapsed', 'd_TSeqT')
	df.vr[, mu:= predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(df.vr, select=d_TSeqT)), what='mu', type='link')]
	df.vr[, sigma:= predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(df.vr, select=d_TSeqT)), what='sigma', type='link')]		
	df.vr[, p:= df.vr[, pGA(brl, mu=mu, sigma=sigma, lower.tail=FALSE)]]
	
	setkey( df.vr, Patient, t.Patient )
	df.vr	<- unique(df.vr)
	
	ggplot( df.vr, aes(x=brl)) + geom_histogram()
	ggplot(df.vr, aes(x=p)) + geom_histogram()
	
	df.vr[, table(brl<0.07)/length(brl)]
	#     FALSE       TRUE 
	#	0.06855576 0.93144424 
	df.vr[, table(brl<0.045)/length(brl)]
	#	0.1553931 0.8446069
	df.vr[, table(brl<0.015)/length(brl)]
	#	0.558958 0.441042
	
	df.vr[, table(brl>0.015 & brl<0.07)/length(brl)]
	#	FALSE      TRUE 
	#	0.5095978 0.4904022
	
	df.vr[, c(quantile(d_TSeqT, p=c(0.8, 0.95, 0.99)), max(d_TSeqT))]
	#	2.921351  4.859400  6.925940 10.174000
	
	df.vr[, table(p>0.05)/length(p)]	
	#	0.4570384 0.5429616
}
######################################################################################
project.athena.Fisheretal.composition.intervals.data<- function()
{	
	indir	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal'
	files	<- data.table(	file.MSM=list.files(indir, pattern='*tATHENAmsm.R$'),
			file.SEQ=list.files(indir, pattern='*tATHENAseq.R$'),
			file.YX=list.files(indir, pattern='*YXSEQ*'))
	files	<- files[9,]				
	files[, method:= files[,regmatches(file.MSM,regexpr('3pa[^_]+', file.MSM))]]
	
	file	<- paste(indir, '/', files[i,file.MSM], sep='')
	load(file)
	set(YX.part1, YX.part1[, which(as.numeric(t.period)>4)], 't.period', '4')
	setkey(YX.part1, t.Patient)	
	dfpot	<- subset(YX.part1, select=c(t.Patient,t.period)) 
	setnames(dfpot, 't.Patient', 'Patient')
	setkey(dfpot, Patient, t.period)
	dfpot	<- unique(dfpot)
	setkey(df.all.allmsm, Patient)
	dfpot	<- merge(unique(df.all.allmsm), dfpot, by='Patient')
	set(dfpot, dfpot[, which(is.na(DateDied))], 'DateDied', 2013.3)
	#
	#	Proportions
	#
	ans		<- dfpot[, list(AcuteNA= mean(is.na(isAcute))), by='t.period']
	tmp		<- subset(dfpot, !is.na(AnyT_T1))[, list(PoslRNA_T1= mean(is.na(PoslRNA_T1))), by='t.period']
	ans		<- merge(ans, tmp, by='t.period')
	tmp		<- dfpot[, list(PosCD4_T1= mean(is.na(PosCD4_T1))), by='t.period']
	ans		<- merge(ans, tmp, by='t.period')
	
	dfpot[, list(PosCD4_T1= length(which(is.na(PosCD4_T1))), PosCD4_2ART= length(which(is.na(PosCD4_T1) & AnyT_T1-AnyPos_T1<1))   ), by='t.period']
	#
	#	Frequency
	#
	tmp		<- merge(df.viro.allmsm, subset(dfpot, select=c(Patient, t.period, AnyT_T1, DateDied)), by='Patient', allow.cartesian=TRUE)
	set(tmp, NULL, 'PosRNA', hivc.db.Date2numeric(tmp[,PosRNA]))
	tmp		<- subset(tmp, PosRNA>=AnyT_T1)
	tmp		<- tmp[, list(PosRNAnpy= length(PosRNA)/(DateDied[1]-AnyT_T1[1])), by=c('Patient','t.period')]
	ansf	<- tmp[, list(DATA='RNA', STAT=c('mean',c(0,0.025,0.25,0.5,0.75,0.975,1)), npy= c(mean(PosRNAnpy), quantile(PosRNAnpy, prob=c(0,0.025,0.25,0.5,0.75,0.975,1))) ), by='t.period']	
	#
	tmp		<- merge(df.immu.allmsm, subset(dfpot, select=c(Patient, t.period, AnyT_T1, AnyPos_T1, DateDied)), by='Patient', allow.cartesian=TRUE)
	set(tmp, NULL, 'PosCD4', hivc.db.Date2numeric(tmp[,PosCD4]))
	tmp		<- subset(tmp, PosCD4<=AnyT_T1)
	tmp2	<- tmp[, which(!is.na(AnyT_T1))]
	set(tmp, tmp2, 'DateDied', tmp[tmp2, AnyT_T1])
	tmp		<- subset(tmp, ( DateDied-AnyPos_T1 )>1 )
	tmp		<- tmp[, list(PosCD4npy= length(PosCD4)/(DateDied[1]-AnyPos_T1[1])), by=c('Patient','t.period')]
	tmp		<- tmp[, list(DATA='CD4', STAT=c('mean',c(0,0.025,0.25,0.5,0.75,0.975,1)), npy= c(mean(PosCD4npy), quantile(PosCD4npy, prob=c(0,0.025,0.25,0.5,0.75,0.975,1))) ), by='t.period']
	ansf	<- rbind(ansf, tmp)
	#
	#	table format
	#	
	ansf	<- dcast.data.table(subset(ansf, STAT%in%c('0.025','0.25','0.5','0.75','0.975')), DATA+t.period~STAT, value.var='npy')
	#	check -- those with no PosCD4 are dropped?
	tmp		<- subset(YX.part1, t.period==1 & t<=t.AnyT_T1 & t>=t.AnyPos_T1)
	tmp2	<- merge(tmp, unique(subset(df.all.allmsm, select=c(Patient, PosCD4_T1, DateDied))), by='Patient', all=TRUE)
	subset(tmp2, is.na(PosCD4_T1))
	
	setkey(df.all.allmsm, Patient)
	dfch	<- unique(df.all.allmsm)
	subset(dfch, is.na(PosCD4_T1))
	merge(dfch, data.table(Patient=df.tpairs[, unique(t.Patient)]), by='Patient' )[, length(unique(Patient))]
	
	
	dfall.ch	<- subset(merge(dfch, data.table(Patient=df.tpairs[, unique(t.Patient)]), by='Patient' ),is.na(PosCD4_T1))[, Patient] 
	X.incare.ch	<- subset( subset(X.incare, stage=='Diag' | stage=='ART.started')[, list(s=all(is.na(CD4))), by='t.Patient'], s)[, t.Patient]	
	check		<- data.table(Patient=setdiff(X.incare.ch, dfall.ch))
	merge(df.immu.allmsm, check, by='Patient')[, length(CD4), by='Patient']
	
	length( subset( subset(X.pt, stage=='Diag' | stage=='ART.started')[, list(s=all(is.na(CD4))), by='t.Patient'], s)[, t.Patient] )
	length( subset( subset(YX.part1, stage=='Diag' | stage=='ART.started')[, list(s=all(is.na(CD4))), by='t.Patient'], s)[, t.Patient] )
	YXch		<- unique(subset( subset(YX.part1, stage=='Diag')[, list(s=all(is.na(CD4))), by='t.Patient'], s, t.Patient))
	Ymch		<- subset(YX.m2[, list(s=all(CD4t=='Dt.NA')), by=c('t.Patient','Patient')], s)[, unique(t.Patient)]
	check		<- setdiff(YXch[,t.Patient], Ymch)
	subset(YX.m2, t.Patient%in%check[1])
	
	length( subset( subset(YX.part1, stage=='Diag')[, list(s=all(is.na(CD4))), by='t.Patient'], s)[, t.Patient] )
	#merge(X.incare, data.table(t.Patient=c('M15922','M26930','M27268')), by='t.Patient')
	merge(X.incare, data.table(t.Patient=c('M15922')), by='t.Patient')
	subset(df.immu.allmsm, Patient=='M15922')
	subset(X.cd4, Patient=='M15922')
	
	merge(data.table(Patient=subset(X.viro, is.na(lRNA_T1.supp), t.Patient)[, unique(t.Patient)]), subset(df.all, select=c(Patient, AnyT_T1)), by='Patient')
	
	length( subset( subset(YX.m2, stage=='Diag' | stage=='ART.started')[, list(s=all(is.na(CD4))), by='t.Patient'], s)[, t.Patient] )
	
	YX.part1x	<- project.athena.Fisheretal.YX.model2.stratify.VLmxwindow(YX.part1, df.all, df.viro, df.immu, indircov, lRNA.supp=method.lRNA.supp )
	length( subset( subset(YX.part1x, stage=='Diag' | stage=='ART.started')[, list(s=all(is.na(CD4))), by='t.Patient'], s)[, t.Patient] )
	
	length(subset(YX.m2[, list(s=all(CD4t=='Dt.NA')), by=c('t.Patient','Patient')], s)[, unique(t.Patient)] )
	
	length( df.tpairs[, unique(t.Patient)] )
	length( X.incare[, unique(t.Patient)] )
}
######################################################################################
project.athena.Fisheretal.composition.potentialtranspairs.table<- function()
{	
	indir	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal'
	files	<- data.table(	file.MSM=list.files(indir, pattern='*tATHENAmsm.R$'),
			file.SEQ=list.files(indir, pattern='*tATHENAseq.R$'),
			file.YX=list.files(indir, pattern='*YXSEQ*'))
	files	<- files[c(1,6,9),]				
	files[, method:= files[,regmatches(file.MSM,regexpr('3pa[^_]+', file.MSM))]]
	
	ans.n	<- lapply(seq_len(nrow(files)), function(i)
			{				
				file	<- paste(indir, '/', files[i,file.YX], sep='')
				load(file)
				set(YX, YX[, which(as.numeric(t.period)>4)], 't.period', '4')
				setkey(YX, Patient, t.Patient)	
				tmp		<- subset(unique(YX), select=c(Patient, t.Patient, FASTASampleCode, t.FASTASampleCode, t.Trm, t.period))
				ans.n	<- tmp[, list(nRwPr= length(unique(Patient)), nPrT= length(unique(t.Patient)), nPrP= length(Patient)), by='t.period']				
				ans.n	<- rbind(ans.n, tmp[, list(t.period='All', nRwPr= length(unique(Patient)), nPrT= length(unique(t.Patient)), nPrP= length(Patient))])
				#	
				file	<- paste(indir, '/', files[i,file.SEQ], sep='')
				load(file)
				set(YX.part1, YX.part1[, which(as.numeric(t.period)>4)], 't.period', '4')
				setkey(YX.part1, Patient, t.Patient)	
				tmp		<- subset(unique(YX.part1), select=c(Patient, t.Patient, t.period))
				tmp2	<- tmp[, list(nRs= length(unique(Patient)), nPoTs= length(unique(t.Patient)), nPoPs= length(Patient)), by='t.period']				
				tmp2	<- rbind(tmp2, tmp[, list(t.period='All', nRs= length(unique(Patient)), nPoTs= length(unique(t.Patient)), nPoPs= length(Patient))])				
				ans.n	<- merge(ans.n, tmp2, by='t.period')
				#	
				file	<- paste(indir, '/', files[i,file.MSM], sep='')
				load(file)
				set(YX.part1, YX.part1[, which(as.numeric(t.period)>4)], 't.period', '4')
				setkey(YX.part1, Patient, t.Patient)	
				tmp		<- subset(unique(YX.part1), select=c(Patient, t.Patient, t.period))
				tmp2	<- tmp[, list(nR= length(unique(Patient)), nPoT= length(unique(t.Patient)), nPoP= length(Patient)), by='t.period']
				tmp2	<- rbind(tmp2, tmp[, list(t.period='All', nR= length(unique(Patient)), nPoT= length(unique(t.Patient)), nPoP= length(Patient))])
				ans.n	<- merge(ans.n, tmp2, by='t.period')
				ans.n[, method:=files[i,method]]
				ans.n
			})
	ans.n	<- do.call('rbind', ans.n)
	
	tmp		<- merge(tmp, subset(df.all.allmsm, select=c(FASTASampleCode, AnyPos_T1, PosSeqT)), by='FASTASampleCode')
	setnames(tmp, c('FASTASampleCode','t.FASTASampleCode','PosSeqT'),  c('r.FASTASampleCode','FASTASampleCode','r.PosSeqT'))
	tmp		<- merge(tmp, subset(df.all.allmsm, select=c(FASTASampleCode, PosSeqT)), by='FASTASampleCode')
	setnames(tmp, c('FASTASampleCode','r.FASTASampleCode','r.PosSeqT','PosSeqT'),  c('t.FASTASampleCode','FASTASampleCode','PosSeqT','t.PosSeqT'))
	tmp[, t.elapsed.RPrT:= abs(PosSeqT-AnyPos_T1+0.5)+abs(t.PosSeqT-AnyPos_T1+0.5) ]
	ans.s	<- subset(tmp, select=c(t.elapsed.RPrT))
}
######################################################################################
project.athena.Fisheretal.composition.recipients.repr.info<- function(df.ar, quantiles)
{
	#	Age
	tmp			<- df.ar[, list( STAT= c('mean', quantiles), AgeT= round(c(mean(AgeT), quantile(AgeT, prob=quantiles)), d=1)), by='t.period']
	tmp2		<- df.ar[, list( t.period='Overall', STAT= c('mean', quantiles), AgeT= round(c(mean(AgeT), quantile(AgeT, prob=quantiles)), d=1))]
	df.ari		<- rbind(tmp, tmp2)	
	#	First CD4
	tmp			<- subset(df.ar, !is.na(CD4_T1) & (PosCD4_T1-AnyPos_T1)<1)[, list( STAT= c('mean', quantiles), CD4_T1= round(c(mean(CD4_T1), quantile(CD4_T1, prob=quantiles)))), by='t.period']
	tmp2		<- subset(df.ar, !is.na(CD4_T1) & (PosCD4_T1-AnyPos_T1)<1)[, list( t.period='Overall', STAT= c('mean', quantiles), CD4_T1= round(c(mean(CD4_T1), quantile(CD4_T1, prob=quantiles))))]
	tmp			<- rbind(tmp, tmp2)
	df.ari		<- merge(df.ari, tmp, by=c('t.period','STAT'))
	#	First Viral load
	tmp			<- subset(df.ar, !is.na(lRNA_T1) & (PoslRNA_T1-AnyPos_T1)<1)[, list( STAT= c('mean', quantiles), lRNA_T1= round(c(mean(lRNA_T1), quantile(lRNA_T1, prob=quantiles)), d=1)), by='t.period']
	tmp2		<- subset(df.ar, !is.na(lRNA_T1) & (PoslRNA_T1-AnyPos_T1)<1)[, list( t.period='Overall', STAT= c('mean', quantiles), lRNA_T1= round(c(mean(lRNA_T1), quantile(lRNA_T1, prob=quantiles)), d=1))]
	tmp			<- rbind(tmp, tmp2)
	df.ari		<- merge(df.ari, tmp, by=c('t.period','STAT'))
	#	% with last neg test in 12mo
	tmp			<- df.ar[, list( STAT= c('mean', quantiles), NegT12P= round(100*c(mean(!is.na(NegT) & (AnyPos_T1-NegT)<1), quantile(!is.na(NegT) & (AnyPos_T1-NegT)<1, prob=quantiles)), d=1) ), by='t.period']
	tmp2		<- df.ar[, list( t.period='Overall', STAT= c('mean', quantiles), NegT12P= round(100*c(mean(!is.na(NegT) & (AnyPos_T1-NegT)<1), quantile(!is.na(NegT) & (AnyPos_T1-NegT)<1, prob=quantiles)), d=1) )]
	tmp			<- rbind(tmp, tmp2)
	df.ari		<- merge(df.ari, tmp, by=c('t.period','STAT'))
	#	% from Amsterdam
	tmp			<- df.ar[, list( STAT= c('mean', quantiles), HAmst= round(100*c(mean(!is.na(RegionHospital) & RegionHospital=='Amst'), quantile(!is.na(RegionHospital) & RegionHospital=='Amst', prob=quantiles)), d=1) ), by='t.period']
	tmp2		<- df.ar[, list( t.period='Overall', STAT= c('mean', quantiles), HAmst= round(100*c(mean(!is.na(RegionHospital) & RegionHospital=='Amst'), quantile(!is.na(RegionHospital) & RegionHospital=='Amst', prob=quantiles)), d=1) )]
	tmp			<- rbind(tmp, tmp2)
	df.ari		<- merge(df.ari, tmp, by=c('t.period','STAT'))		
	#	% country infection NL
	tmp			<- df.ar[, {
								z	<- CountryInfection[!is.na(CountryInfection)]
								list( STAT= c('mean', quantiles), INFNL= round(100*c(mean(z=='NL'), quantile(z=='NL', prob=quantiles)), d=1) )	
							}, by='t.period']
	tmp2		<- df.ar[, {
								z	<- CountryInfection[!is.na(CountryInfection)]
								list( t.period='Overall', STAT= c('mean', quantiles), INFNL= round(100*c(mean(z=='NL'), quantile(z=='NL', prob=quantiles)), d=1) )		
							}]
	tmp			<- rbind(tmp, tmp2)
	df.ari		<- merge(df.ari, tmp, by=c('t.period','STAT'))
	df.ari
}
######################################################################################
project.athena.Fisheretal.composition.recipients.repr<- function()
{	
	quantiles	<- c(0,0.025,0.25,0.5,0.75,0.975,1)
	#
	#	 I re-set inaccurate testing dates -- get df.all.allmsm with original testing dates first
	#
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
	method.recentctime		<- '2011-01-01'
	method.Acute			<- 'higher'	#'central'#'empirical'
	method.minQLowerU		<- 0.135
	method.use.AcuteSpec	<- 1
	#
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	infiletree				<- paste(infile,"examlbs500",sep="_")
	insignat				<- "Wed_Dec_18_11:37:00_2013"							
	t.recent.endctime		<- hivc.db.Date2numeric(as.Date(method.recentctime))	
	t.recent.endctime		<- floor(t.recent.endctime) + floor( (t.recent.endctime%%1)*100 %/% (t.period*100) ) * t.period	
	
	adjust.AcuteByNegT			<- 1
	any.pos.grace.yr			<- Inf	
	adjust.minSCwindow			<- NA		#differs from main study
	adjust.NegTByDetectability	<- NA		#differs from main study
	
	tmp				<- project.athena.Fisheretal.select.denominator(indir, infile, insignat, indircov, infile.cov.study, infile.viro.study, infile.immu.study, infile.treatment.study, 
			infiletree=infiletree, adjust.AcuteByNegT=adjust.AcuteByNegT, adjust.NegT4Acute=NA, adjust.NegTByDetectability=adjust.NegTByDetectability, adjust.minSCwindow=adjust.minSCwindow, adjust.AcuteSelect=c('Yes','Maybe'), use.AcuteSpec=method.use.AcuteSpec, 
			t.recent.endctime=t.recent.endctime, t.recent.startctime=t.recent.startctime)	
	df.all			<- tmp$df.all	
	df.denom.CLU	<- tmp$df.select
	df.denom.SEQ	<- tmp$df.select.SEQ	
	ri.CLU			<- unique(subset(df.denom.CLU, select=Patient))
	ri.SEQ			<- unique(subset(df.denom.SEQ, select=Patient))
	df.viro			<- tmp$df.viro
	df.immu			<- tmp$df.immu
	df.treatment	<- tmp$df.treatment	
	#
	tmp					<- project.athena.Fisheretal.select.denominator(	indir, infile, insignat, indircov, infile.cov.all, infile.viro.all, infile.immu.all, infile.treatment.all, 
			infiletree=NULL, adjust.AcuteByNegT=adjust.AcuteByNegT, adjust.NegT4Acute=NA, adjust.NegTByDetectability=adjust.NegTByDetectability, adjust.minSCwindow=adjust.minSCwindow, adjust.AcuteSelect=c('Yes','Maybe'), use.AcuteSpec=method.use.AcuteSpec,
			t.recent.endctime=t.recent.endctime, t.recent.startctime=t.recent.startctime,
			df.viro.part=df.viro, df.immu.part=df.immu, df.treatment.part=df.treatment, df.all.part=df.all)	
	df.all.allmsm		<- tmp$df.all
	setkey(df.all.allmsm, Patient)
	#
	#	all diag
	#
	df.ad	<- subset(df.all.allmsm, AnyPos_T1<2011 & AnyPos_T1>=1996.6 & Trm%in%c('MSM','BI')) 
	setkey(df.ad, Patient)
	df.ad	<- unique(df.ad)	
	df.ad[, AgeT:= AnyPos_T1-DateBorn]
	df.ad[, t.period:= cut(AnyPos_T1, breaks=c(1996, 2006.5, 2008, 2009.5, 2011))]	
	#	all rec
	df.ar	<- subset(df.all.allmsm, isAcute=='Yes' & AnyPos_T1<2011 & AnyPos_T1>=1996.6 & Trm%in%c('MSM','BI')) 
	setkey(df.ar, Patient)
	df.ar	<- unique(df.ar)	
	df.ar[, AgeT:= AnyPos_T1-DateBorn]
	df.ar[, t.period:= cut(AnyPos_T1, breaks=c(1996, 2006.5, 2008, 2009.5, 2011))]
	#	rec w pr
	indir	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal'
	files	<- data.table(	file.MSM=list.files(indir, pattern='*tATHENAmsm.R$'),
			file.SEQ=list.files(indir, pattern='*tATHENAseq.R$'),
			file.YX=list.files(indir, pattern='*YXSEQ*'))
	file	<- paste(indir, '/', files[9,file.YX], sep='')
	load(file)
	df.rr	<- merge(df.ar, unique(subset(YX, select=Patient)), by='Patient')
	df.rr[, round(table(Acute_Spec, useNA='if')/nrow(df.rr), d=2)]
	
	df.adi	<- project.athena.Fisheretal.composition.recipients.repr.info(df.ad, quantiles)
	df.ari	<- project.athena.Fisheretal.composition.recipients.repr.info(df.ar, quantiles)
	df.rri	<- project.athena.Fisheretal.composition.recipients.repr.info(df.rr, quantiles)
	#	re-format all recent
	ans			<- dcast.data.table(df.ari, t.period~STAT, value.var='AgeT')
	ans[, STAT:='AgeT']
	ans[, DATA:='All']
	tmp			<- dcast.data.table(df.ari, t.period~STAT, value.var='CD4_T1')
	tmp[, STAT:='CD4_T1']
	tmp[, DATA:='All']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.ari, t.period~STAT, value.var='lRNA_T1')
	tmp[, STAT:='lRNA_T1']
	tmp[, DATA:='All']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.ari, t.period~STAT, value.var='NegT12P')
	tmp[, STAT:='NegT12P']
	tmp[, DATA:='All']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.ari, t.period~STAT, value.var='INFNL')
	tmp[, STAT:='INFNL']
	tmp[, DATA:='All']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.ari, t.period~STAT, value.var='HAmst')
	tmp[, STAT:='HAMST']
	tmp[, DATA:='All']
	ans			<- rbind(ans, tmp)		
	#	re-format those with prob transm
	tmp			<- dcast.data.table(df.rri, t.period~STAT, value.var='AgeT')
	tmp[, STAT:='AgeT']
	tmp[, DATA:='wProb']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.rri, t.period~STAT, value.var='CD4_T1')
	tmp[, STAT:='CD4_T1']
	tmp[, DATA:='wProb']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.rri, t.period~STAT, value.var='lRNA_T1')
	tmp[, STAT:='lRNA_T1']
	tmp[, DATA:='wProb']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.rri, t.period~STAT, value.var='NegT12P')
	tmp[, STAT:='NegT12P']
	tmp[, DATA:='wProb']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.rri, t.period~STAT, value.var='INFNL')
	tmp[, STAT:='INFNL']
	tmp[, DATA:='wProb']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.rri, t.period~STAT, value.var='HAmst')
	tmp[, STAT:='HAMST']
	tmp[, DATA:='wProb']
	ans			<- rbind(ans, tmp)		
	#	re-format diagnosed
	tmp			<- dcast.data.table(df.adi, t.period~STAT, value.var='AgeT')
	tmp[, STAT:='AgeT']
	tmp[, DATA:='Diag']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.adi, t.period~STAT, value.var='CD4_T1')
	tmp[, STAT:='CD4_T1']
	tmp[, DATA:='Diag']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.adi, t.period~STAT, value.var='lRNA_T1')
	tmp[, STAT:='lRNA_T1']
	tmp[, DATA:='Diag']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.adi, t.period~STAT, value.var='NegT12P')
	tmp[, STAT:='NegT12P']
	tmp[, DATA:='Diag']
	ans			<- rbind(ans, tmp)
	tmp			<- dcast.data.table(df.adi, t.period~STAT, value.var='INFNL')
	tmp[, STAT:='INFNL']
	tmp[, DATA:='Diag']
	tmp			<- dcast.data.table(df.adi, t.period~STAT, value.var='HAmst')
	tmp[, STAT:='HAMST']
	tmp[, DATA:='Diag']
	ans			<- rbind(ans, tmp)	
	#levels(ans$t.period)
	tmp			<- c("Overall","(1996,2006]","(2006,2008]","(2008,2010]","(2010,2011]")
	set(ans, NULL, 't.period', ans[, factor(as.character(t.period), levels=tmp, labels=tmp)])
	#ans[, unique(DATA)]
	tmp			<- c("wProb","All","Diag")
	set(ans, NULL, 'DATA', ans[, factor(DATA, levels=tmp, labels=tmp)])
	setkey(ans, STAT, t.period, DATA)
	
	ans			<- subset(ans, select=c("STAT","t.period","DATA","0.025","0.25","mean","0.75","0.975"))
	file		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150303_RecipMSMRepresentative.csv'
	write.csv(ans, file=file, eol="\r\n", row.names=FALSE)
}
######################################################################################
project.athena.Fisheretal.composition.probtranspairs.Ard<- function()
{	
	#load('~/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Ac=MY_D=35_sasky_2011_Wed_Dec_18_11:37:00_2013_Xmsmm2_SEQ_3pa1H1.48C2V100bInfT7.R')
	#tmp		<- subset(YX.m2, CD4c=='Lost')
	#df.nc	<- tmp[, list(WindowStart=min(t)), by='t.Patient']
	#setnames(df.nc, 't.Patient','Patient')
	tmp		<- subset(df.all.allmsm, select=c(Patient, DateDied, DateLastContact))
	setkey(tmp, Patient)	
	df.nc	<- data.table( Patient= subset(X.incare, contact=='No' & t>2010 & t<2010.5)[, unique(t.Patient)] )	
	df.nc	<- merge(df.nc, unique(tmp), by='Patient')
	
	#	
	set(df.immu.allmsm, NULL, 'PosCD4', hivc.db.Date2numeric(df.immu.allmsm[,PosCD4]))
	set(df.viro.allmsm, NULL, 'PosRNA', hivc.db.Date2numeric(df.viro.allmsm[,PosRNA]))
	#
	df.nc	<- merge(df.nc, df.immu.allmsm, by='Patient', all.x=1)
	df.nc[, T2CD4:= abs(WindowStart-PosCD4)]
	set(df.nc, df.nc[, which(is.na(T2CD4))], 'T2CD4', 1000)
	setkey(df.nc, Patient, WindowStart, T2CD4)
	df.nc[, list(WindowStart=WindowStart[1], DateDied=DateDied[1], DateLastContact=DateLastContact[1], T2CD4=T2CD4[1]), by=c('Patient','WindowStart')]
	
	
	df.immu.allmsm
	file		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150319_NoContactPotentialTransmitters.csv'
	write.csv(df.nc, file=file, eol='\r\n')	
}
######################################################################################
project.athena.Fisheretal.composition.conftranspairs.divergence<- function()
{	
	indir	<- '/Users/Oliver/duke/2014_HIV_TChainBelgium'
	infile	<- '140921_set7_INFO_TRM.R'
	file	<- paste(indir, '/',infile, sep='')
	load(file)
	#	exclude B->A, same as A->B
	trm.pol.nA		<- subset(trm.pol, withA==FALSE, select=c(d_SeqT, d_TSeqT, BRL))			
	trm.pol.nA		<- rbind(trm.pol.nA, data.table(BRL=0.0025, d_SeqT=seq(2,4,len=200), d_TSeqT=seq(2,8,len=200)), fill=TRUE)
	#trm.pol.nA		<- rbind(trm.pol.nA, data.table(BRL=0.001, d_SeqT=rep(0.1,40), d_TSeqT=rep(0.1,40)), fill=TRUE)			
	#trm.pol.GA32	<- gamlss(as.formula('BRL ~ bs(d_TSeqT, degree=4)'), sigma.formula=as.formula('~ bs(d_TSeqT, degree=2)'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity',sigma.link='identity'), n.cyc = 40)			
	trm.pol.GA11e	<- gamlss(as.formula('BRL ~ d_TSeqT-1'), sigma.formula=as.formula('~ d_TSeqT'), data=as.data.frame(trm.pol.nA), family=GA(mu.link='identity',sigma.link='identity'), n.cyc = 40)			
	trm.pol.GA		<- trm.pol.GA11e
	trm.pol.p2		<- data.table(d_TSeqT=seq(0,18,0.1))
	trm.pol.p2		<- data.table(d_TSeqT=seq(0.1,18,0.1))	
	trm.pol.p2[, BRL_p:=predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(trm.pol.p2), type='response', se.fit=FALSE)]
	trm.pol.p2[, mu:=predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(trm.pol.p2, select=d_TSeqT)), type='link', what='mu',se.fit=FALSE)]
	trm.pol.p2[, sig:=predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(trm.pol.p2, select=d_TSeqT)), type='link', what='sigma',se.fit=FALSE)]	
	trm.pol.p2		<- merge(trm.pol.p2, as.data.table(expand.grid(d_TSeqT=trm.pol.p2[, d_TSeqT], p=c(1, 2.5, 10, 25, 75, 90, 97.5, 99)/100 )), by='d_TSeqT')
	trm.pol.p2[, q:=trm.pol.p2[, qGA(p, mu=mu, sigma=sig)]]	
	trm.pol.p3		<- copy(trm.pol.p2)
	set(trm.pol.p3, NULL, 'p', trm.pol.p3[, paste('q',100*p,sep='')])
	trm.pol.p3		<- dcast.data.table(trm.pol.p3, d_TSeqT~p, value.var='q' )
		
	ggplot(trm.pol, aes(x=d_TSeqT)) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=q10, ymax=q90), fill="#E41A1C", alpha=0.5) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=q2.5, ymax=q10), fill="#E41A1C", alpha=0.3) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=q90, ymax=q97.5), fill="#E41A1C", alpha=0.3) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=q1, ymax=q2.5), fill="#E41A1C", alpha=0.1) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=q97.5, ymax=q99), fill="#E41A1C", alpha=0.1) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=q99, ymax=0.15), fill="#377EB8", alpha=0.3) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=0, ymax=q1), fill="#377EB8", alpha=0.3) +
			geom_jitter(aes(y=BRL), size=1.2, alpha=0.5, position = position_jitter(width = .1), data=subset(trm.pol, withA==FALSE & BRL>0.003)) +									
			geom_line(data=trm.pol.p2, aes(y=BRL_p)) +
			#geom_line(data=trm.pol.p3, aes(x=x, y=q2.5), lty='dotted', color=colors[1], size=line.size) + geom_line(data=trm.pol.p3, aes(x=x, y=q97.5), lty='dotted', color=colors[1], size=line.size) +
			#geom_line(data=trm.pol.p3, aes(x=x, y=q1), lty='twodash') + geom_line(data=trm.pol.p3, aes(x=x, y=q99), lty='twodash') +
			#geom_line(data=trm.pol.p3, aes(x=x, y=q25), lty='dashed') + geom_line(data=trm.pol.p3, aes(x=x, y=q75), lty='dashed') +
			#geom_line(data=trm.pol.p3, aes(x=x, y=q10), lty='twodash', color=colors[1], size=line.size) + geom_line(data=trm.pol.p3, aes(x=x, y=q90), lty='twodash', color=colors[1], size=line.size) +
			#geom_vline(xintercept=10.17, color=colors[2] ) +					
			scale_x_continuous(breaks=seq(0,20,2), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,0.2,0.02), expand=c(0,0)) +
			coord_trans(limy=c(0,0.1), limx=c(0,13) ) +
			theme_bw() + 
			#labs(x='time elapsed\n(years)', y='genetic distance\n(subst/site)', title='sequence pairs\nbetween transmitters and recipients\nin epidemiologically confirmed pairs\n') +
			labs(x='time elapsed\n(years)', y='genetic distance\n(subst/site)', title='epidemiologically confirmed\ntransmission pairs\n') +
			theme(panel.grid.major=element_line(colour="grey70", size=0.4), panel.grid.minor=element_line(colour="grey70", size=0.2))
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150320_ConfPairsBrl2.pdf'
	ggsave(file=file, w=5, h=5)	
}
######################################################################################
project.athena.Fisheretal.composition.probtranspairs.divergence<- function()
{	
	indir	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal'
	files	<- data.table(	file.MSM=list.files(indir, pattern='*tATHENAmsm.R$'),
			file.SEQ=list.files(indir, pattern='*tATHENAseq.R$'),
			file.YX=list.files(indir, pattern='*YXSEQ*'))
	file	<- paste(indir, '/', files[3,file.YX], sep='')
	load(file)
	tmp		<- unique(subset(YX, select=c(Patient, t.Patient, cluster, score.Y,  telapsed, brl, t.period)))
	
	file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_data/TchainBelgium_set7_pol_GA11emodel_nA_INFO.R'
	load(file)
	trm.pol.p2		<- data.table(d_TSeqT=seq(0.1,18,0.1))	
	trm.pol.p2[, BRL_p:=predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(trm.pol.p2), type='response', se.fit=FALSE)]
	trm.pol.p2[, mu:=predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(trm.pol.p2, select=d_TSeqT)), type='link', what='mu',se.fit=FALSE)]
	trm.pol.p2[, sig:=predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(trm.pol.p2, select=d_TSeqT)), type='link', what='sigma',se.fit=FALSE)]	
	trm.pol.p2		<- merge(trm.pol.p2, as.data.table(expand.grid(d_TSeqT=trm.pol.p2[, d_TSeqT], p=c(1, 2.5, 10, 25, 75, 90, 97.5, 99)/100 )), by='d_TSeqT')
	trm.pol.p2[, q:=trm.pol.p2[, qGA(p, mu=mu, sigma=sig)]]
	
	trm.pol.p3		<- copy(trm.pol.p2)
	set(trm.pol.p3, NULL, 'p', trm.pol.p3[, paste('q',100*p,sep='')])
	trm.pol.p3		<- dcast.data.table(trm.pol.p3, d_TSeqT~p, value.var='q' )
	
	ggplot(tmp) + 
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=q10, ymax=q90), fill="#E41A1C", alpha=0.5) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=q2.5, ymax=q10), fill="#E41A1C", alpha=0.3) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=q90, ymax=q97.5), fill="#E41A1C", alpha=0.3) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=q1, ymax=q2.5), fill="#E41A1C", alpha=0.1) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=q97.5, ymax=q99), fill="#E41A1C", alpha=0.1) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=q99, ymax=0.15), fill="#377EB8", alpha=0.3) +
			geom_ribbon(data=trm.pol.p3, aes(x=d_TSeqT, ymin=0, ymax=q1), fill="#377EB8", alpha=0.3) +
			geom_point(size=1.2, alpha=0.75, aes(x=telapsed, y=brl)) +
			#geom_hline(yintercept=c(0.02, 0.04), colour="#E41A1C") +
			scale_x_continuous(breaks=seq(0,20,2), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,0.2,0.02), expand=c(0,0)) +
			coord_trans(limy=c(0,0.1), limx=c(0,13) ) +
			theme_bw() + 
			#labs(x='time elapsed\n(years)', y='genetic distance\n(subst/site)', title='sequence pairs\nbetween transmitters and recipients\nin phylogenetically probable pairs\n') +
			labs(x='time elapsed\n(years)', y='genetic distance\n(subst/site)', title='phylogenetically probable\ntransmission pairs\n') +
			theme(panel.grid.major=element_line(colour="grey70", size=0.4), panel.grid.minor=element_line(colour="grey70", size=0.2))
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150320_ProbPairsBrl2.pdf'
	ggsave(file=file, w=5, h=5)
	
	
	ggplot(tmp) + 
			geom_ribbon(data=trm.pol.p3, aes(x=x, ymin=q10, ymax=q90), fill="#377EB8") +
			geom_ribbon(data=trm.pol.p3, aes(x=x, ymin=q2.5, ymax=q10), fill="#377EB8", alpha=0.6) +
			geom_ribbon(data=trm.pol.p3, aes(x=x, ymin=q90, ymax=q97.5), fill="#377EB8", alpha=0.6) +
			geom_ribbon(data=trm.pol.p3, aes(x=x, ymin=q1, ymax=q2.5), fill="#377EB8", alpha=0.3) +
			geom_ribbon(data=trm.pol.p3, aes(x=x, ymin=q97.5, ymax=q99), fill="#377EB8", alpha=0.3) +			
			geom_point(size=1.2, alpha=0.75, aes(x=telapsed, y=brl)) +
			geom_hline(yintercept=c(0.02, 0.04), colour="#E41A1C") +
			scale_x_continuous(breaks=seq(0,20,2), limit=c(0,11), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,0.2,0.02), limit=c(0,0.15), expand=c(0,0)) +
			theme_bw() + labs(x='time elapsed\n(years)', y='evolutionary divergence\n(subst/site)') +
			theme(panel.grid.major=element_line(colour="grey70", size=0.4), panel.grid.minor=element_line(colour="grey70", size=0.4))
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150302_ProbPairsBrl.pdf'
	ggsave(file=file, w=8, h=8)
	
	require(scales)
	tmp2	<- sqsq_trans<- function(){		trans_new("sqsq", transform = function(x){x^(1/4)}, inverse = function(x){x^4})	}
	ggplot(tmp, aes(x=telapsed, y=brl, colour=score.Y, size=score.Y)) + 
			geom_point(alpha=0.75) +
			scale_x_continuous(breaks=seq(0,20,2), limit=c(0,11), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,0.2,0.02), limit=c(0,0.15), expand=c(0,0)) +
			scale_colour_gradientn(name='relative\ntransmission\nprobability', colours=c("#377EB8","#FF7F00"), trans = tmp2(), breaks=c(1, 10, 100, 400)) +
			scale_size_continuous(range = c(1,4), trans = sqrt_trans(), guide=FALSE) +
			theme_bw() + labs(x='time elapsed\n(years)', y='evolutionary divergence\n(subst/site)') +
			theme( legend.position=c(1,1), legend.justification=c(1,1), panel.grid.major=element_line(colour="grey70", size=0.4), panel.grid.minor=element_line(colour="grey70", size=0.4))
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150302_ProbPairsPhyloLkl.pdf'
	ggsave(file=file, w=5, h=5)
	
	
	
	Y.brl	<- data.table(d_TSeqT=c(1, 2, 3, 5, 7, 10))
	Y.brl[, mu:= predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='mu', type='link')]
	Y.brl[, sigma:= predict(trm.pol.GA, data=trm.pol.nA, newdata=as.data.frame(subset(Y.brl, select=d_TSeqT)), what='sigma', type='link') ]
	Y.brl	<- merge(Y.brl, as.data.table(expand.grid(d_TSeqT=c(1, 2, 3, 5, 7, 10), d=c(seq(0.001, 0.01, 0.0001),seq(0.011, 0.1, 0.001)))), by='d_TSeqT')	
	Y.brl[, lkl:= Y.brl[, dGA(d, mu=mu, sigma=sigma)]]
	
	ggplot(Y.brl, aes(x=d, y=lkl, colour= factor(d_TSeqT), group=factor(d_TSeqT))) + geom_line(size=0.8) +
			scale_x_continuous(breaks=seq(0, 0.2, 0.02), expand=c(0,0)) +
			scale_y_continuous(expand=c(0,0)) +
			scale_colour_brewer(name='time elapsed\n(years)', palette='Dark2') +
			labs(x='evolutionary divergence\n(subst/site)', y='relative transmission probability') +
			theme_bw() + theme(legend.position=c(1,1), legend.justification=c(1,1),panel.grid.major=element_line(colour="grey70", size=0.4)) 
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150302_PhyloLkl.pdf'
	ggsave(file=file, w=5, h=5)
	
	
	
	
	#tmp[, mean(brl<=0.07)]
	#[1] 0.9765733
	#> tmp[, mean(brl<=0.04)]
	#[1] 0.8700046
	#> tmp[, mean(brl<=0.02)]
	#[1] 0.597152
	#> tmp[,c(mean(telapsed), quantile(telapsed, prob=c(0.25, 0.5, 0.75, 1)))]
	#[1] 2.0122714  0.8706756  1.6013511  2.6621756 10.8690000 
	
	#geom_line(data=trm.pol.p2, aes(y=BRL_p)) +
			
}
######################################################################################
project.athena.Fisheretal.composition.probtranspairs.table<- function()
{	
	indir	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal'
	files	<- data.table(	file.MSM=list.files(indir, pattern='*tATHENAmsm.R$'),
							file.SEQ=list.files(indir, pattern='*tATHENAseq.R$'),
							file.YX=list.files(indir, pattern='*YXSEQ.*T7.R$'))
					
	#files	<- files[grepl('YXSEQ3pa1H1.48C1V100bInfT7', file.YX) | grepl('H1.48C3.*bInfT7', file.YX) | grepl('H1.48C1.*bInfT7', file.YX) | grepl('H1.48C2.*bInfs0.7T7', file.YX) | grepl('H1.48C2.*bInfs0.85T7', file.YX), ]
	files	<- files[grepl('YXSEQ3pa1H1.48C[1-3]V100bInfT7', file.YX) ]
	files[, method:= files[,regmatches(file.MSM,regexpr('3pa[^_]+', file.MSM))]]
	
	ans.n	<- lapply(seq_len(nrow(files)), function(i)
			{				
				file	<- paste(indir, '/', files[i,file.YX], sep='')
				load(file)
				set(YX, YX[, which(as.numeric(t.period)>4)], 't.period', '4')
				setkey(YX, Patient, t.Patient)	
				tmp		<- subset(unique(YX), select=c(Patient, t.Patient, FASTASampleCode, t.FASTASampleCode, t.Trm, t.period))
				ans.n	<- tmp[, list(nRwPr= length(unique(Patient)), nPrT= length(unique(t.Patient)), nPrP= length(Patient)), by='t.period']				
				ans.n	<- rbind(ans.n, tmp[, list(t.period='All', nRwPr= length(unique(Patient)), nPrT= length(unique(t.Patient)), nPrP= length(Patient))])
				#	
				file	<- paste(indir, '/', files[i,file.SEQ], sep='')
				load(file)
				set(YX.part1, YX.part1[, which(as.numeric(t.period)>4)], 't.period', '4')
				setkey(YX.part1, Patient, t.Patient)	
				tmp		<- subset(unique(YX.part1), select=c(Patient, t.Patient, t.period))
				tmp2	<- tmp[, list(nRs= length(unique(Patient)), nPoTs= length(unique(t.Patient)), nPoPs= length(Patient)), by='t.period']				
				tmp2	<- rbind(tmp2, tmp[, list(t.period='All', nRs= length(unique(Patient)), nPoTs= length(unique(t.Patient)), nPoPs= length(Patient))])				
				ans.n	<- merge(ans.n, tmp2, by='t.period')
				#	
				file	<- paste(indir, '/', files[i,file.MSM], sep='')
				load(file)
				set(YX.part1, YX.part1[, which(as.numeric(t.period)>4)], 't.period', '4')
				setkey(YX.part1, Patient, t.Patient)	
				tmp		<- subset(unique(YX.part1), select=c(Patient, t.Patient, t.period))
				tmp2	<- tmp[, list(nR= length(unique(Patient)), nPoT= length(unique(t.Patient)), nPoP= length(Patient)), by='t.period']
				tmp2	<- rbind(tmp2, tmp[, list(t.period='All', nR= length(unique(Patient)), nPoT= length(unique(t.Patient)), nPoP= length(Patient))])
				ans.n	<- merge(ans.n, tmp2, by='t.period')
				ans.n[, method:=files[i,method]]
				ans.n
			})
	ans.n	<- do.call('rbind', ans.n)
	ans.p	<- copy(ans.n)
	
		
	ans.p[, pRwPr:= round(nRwPr/nRs*100,d=2)]
	ans.p[, pPrT:= round(nPrT/nPoTs*100,d=2)]
	ans.p[, pPrP:= round(nPrP/nPoPs*100,d=2)]
	
	tmp		<- dcast.data.table(subset(ans.p, select=c(t.period, method, nRwPr)), t.period~method, value.var='nRwPr')
	setnames(tmp, names(tmp)[-1], paste('n',names(tmp)[-1],sep=''))
	tmp2	<- dcast.data.table(subset(ans.p, select=c(t.period, method, pRwPr)), t.period~method, value.var='pRwPr')
	setnames(tmp2, names(tmp2)[-1], paste('p',names(tmp2)[-1],sep=''))	
	ans		<- merge(tmp, tmp2, by=c('t.period'))
	ans[, stat:='RwPr']
	tmp		<- dcast.data.table(subset(ans.p, select=c(t.period, method, nPrT)), t.period~method, value.var='nPrT')
	setnames(tmp, names(tmp)[-1], paste('n',names(tmp)[-1],sep=''))
	tmp2	<- dcast.data.table(subset(ans.p, select=c(t.period, method, pPrT)), t.period~method, value.var='pPrT')
	setnames(tmp2, names(tmp2)[-1], paste('p',names(tmp2)[-1],sep=''))	
	tmp		<- merge(tmp, tmp2, by=c('t.period'))
	tmp[, stat:='PrT']
	ans		<- rbind(ans, tmp)
	tmp		<- dcast.data.table(subset(ans.p, select=c(t.period, method, nPrP)), t.period~method, value.var='nPrP')
	setnames(tmp, names(tmp)[-1], paste('n',names(tmp)[-1],sep=''))
	tmp2	<- dcast.data.table(subset(ans.p, select=c(t.period, method, pPrP)), t.period~method, value.var='pPrP')
	setnames(tmp2, names(tmp2)[-1], paste('p',names(tmp2)[-1],sep=''))	
	tmp		<- merge(tmp, tmp2, by=c('t.period'))
	tmp[, stat:='PrP']
	ans		<- rbind(ans, tmp)
	set(ans, NULL, 't.period', ans[, factor(as.character(t.period), levels=c('All','1','2','3','4'), labels=c('All','1','2','3','4'))])
	set(ans, NULL, 'stat', ans[, factor(stat, levels=c('RwPr','PrT','PrP'), labels=c('RwPr','PrT','PrP'))])
	setkey(ans, stat, t.period)
	ans		<- subset(ans, select=c('stat','t.period','n3pa1H1.48C2V100bInfT7','p3pa1H1.48C2V100bInfT7','n3pa1H1.48C3V100bInfT7','p3pa1H1.48C3V100bInfT7','n3pa1H1.48C1V100bInfT7','p3pa1H1.48C1V100bInfT7',
									'n3pa1H1.48C2V100bInfs0.7T7','p3pa1H1.48C2V100bInfs0.7T7','n3pa1H1.48C2V100bInfs0.85T7','p3pa1H1.48C2V100bInfs0.85T7'))
	file		<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150310_PrPairs.csv'
	write.csv(ans, file=file, eol='\r\n')	
	#
	#	multiple transmitters
	#
	file	<- paste(indir, '/', files[2,file.YX], sep='')
	load(file)
	tmp		<- unique(subset(YX, select=c(Patient, t.Patient, t.period)))	
	set(tmp, tmp[, which(as.numeric(t.period)>4)], 't.period', '4')
	
	tmp		<- tmp[, list(NPR= length(unique(t.Patient))), by=c('t.period','Patient')]
	tmp2	<- data.table(t.period=c('1','2','3','4'), legend=c('96/07-06/06','06/07-07/12','08/01-09/06','09/07-10/12'))
	set(tmp2, NULL, 'legend', tmp2[,factor(legend, levels=legend, labels=legend)])
	tmp		<- merge(tmp, tmp2, by='t.period')
	
	ggplot(tmp, aes(x=factor(NPR), fill=legend)) + geom_bar(width=0.7, position = position_dodge(width = 0.8)) +
			theme_bw() + labs(x='Phylogenetically probable transmitters per recipient MSM\n(#)', y='') +
			scale_y_continuous(breaks=seq(0,100,10)) +
			scale_fill_brewer(name='time of diagnosis of\nrecipient MSM', palette='Set2') +
			theme(legend.position=c(1,1), legend.justification=c(1,1), panel.grid.major.y=element_line(colour="grey70", size=0.4), panel.grid.minor.y=element_line(colour="grey70", size=0.4), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150302_MultiPairs.pdf'
	ggsave(file=file, w=10, h=3)
	
	
	tmp		<- merge(tmp, subset(df.all.allmsm, select=c(FASTASampleCode, AnyPos_T1, PosSeqT)), by='FASTASampleCode')
	setnames(tmp, c('FASTASampleCode','t.FASTASampleCode','PosSeqT'),  c('r.FASTASampleCode','FASTASampleCode','r.PosSeqT'))
	tmp		<- merge(tmp, subset(df.all.allmsm, select=c(FASTASampleCode, PosSeqT)), by='FASTASampleCode')
	setnames(tmp, c('FASTASampleCode','r.FASTASampleCode','r.PosSeqT','PosSeqT'),  c('t.FASTASampleCode','FASTASampleCode','PosSeqT','t.PosSeqT'))
	tmp[, t.elapsed.RPrT:= abs(PosSeqT-AnyPos_T1+0.5)+abs(t.PosSeqT-AnyPos_T1+0.5) ]
	ans.s	<- subset(tmp, select=c(t.elapsed.RPrT))
}
######################################################################################
project.athena.Fisheretal.composition.cluster.typeI<- function()
{	
	verbose		<- 1
	resume		<- 1
	patient.n	<- 15700; 	thresh.brl<- 0.096; 	thresh.bs<- 0.8;	opt.brl<- "dist.brl.casc" 
	indircov	<- paste(DATA,"derived",sep='/')
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"
	indir		<- paste(DATA,"tmp",sep='/')
	#
	# get clusters for No Recombination + No Drug resistance mutations + No short sequences, single linkage criterion		
	#						
	infile			<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500"			
	insignat		<- "Wed_Dec_18_11/37/00_2013"
	#
	argv			<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	nsh.clu.pre		<- hivc.prog.get.clustering.precompute()
	#
	resume			<- 0
	argv			<<- hivc.cmd.clustering.tptn(indir, infile, insignat, indircov, infilecov, opt.brl="dist.brl.casc", patient.n=patient.n, resume=resume)
	argv			<<- unlist(strsplit(argv,' '))
	nsh.clu.tptn	<- hivc.prog.get.clustering.TPTN(clu.pre=nsh.clu.pre, with.plot=0)
	#			
	tnn.by.sum		<- nsh.clu.tptn$tnn.by.sum
	fnn.by.sum		<- nsh.clu.tptn$fpn.by.sum
	tn.by.sum		<- nsh.clu.tptn$tn.by.sum
	fp.by.all		<- nsh.clu.tptn$fp.by.all
	fp.by.sum		<- nsh.clu.tptn$fp.by.sum
	colnames(fp.by.all)<- colnames(fp.by.sum)
	fp.by.sum2		<- nsh.clu.tptn$fp.by.sum2
	
	
	fp.by.all		<- data.table(fp.by.all= as.numeric(fp.by.all), bs=factor(rownames(fp.by.all), ordered=1), brl=rep(factor(colnames(fp.by.all), ordered=1), each=nrow(fp.by.all)))
	fp.by.sum		<- data.table(fp.by.sum= as.numeric(fp.by.sum), bs=factor(rownames(fp.by.sum), ordered=1), brl=rep(factor(colnames(fp.by.sum), ordered=1), each=nrow(fp.by.sum)))
	fp.by.sum2		<- data.table(fp.by.sum2= as.numeric(fp.by.sum2), bs=factor(rownames(fp.by.sum2), ordered=1), brl=rep(factor(colnames(fp.by.sum2), ordered=1), each=nrow(fp.by.sum2)))
	tnn.by.sum		<- data.table(tnn= as.numeric(tnn.by.sum), bs=factor(rownames(tnn.by.sum), ordered=1), brl=rep(factor(colnames(tnn.by.sum), ordered=1), each=nrow(tnn.by.sum	)))
	fnn.by.sum		<- data.table(fnn= as.numeric(fnn.by.sum), bs=factor(rownames(fnn.by.sum), ordered=1), brl=rep(factor(colnames(fnn.by.sum), ordered=1), each=nrow(fnn.by.sum	)))
	tn.by.sum		<- data.table(tn= as.numeric(tn.by.sum), bs=factor(rownames(tn.by.sum), ordered=1), brl=rep(factor(colnames(tn.by.sum), ordered=1), each=nrow(tn.by.sum	)))
	df			<- merge(tnn.by.sum, fp.by.all, by=c('bs','brl'))
	df			<- merge(df, tn.by.sum, by=c('bs','brl'))
	df			<- merge(df, fp.by.sum, by=c('bs','brl'))
	df			<- merge(df, fp.by.sum2, by=c('bs','brl'))
	df			<- merge(df, fnn.by.sum, by=c('bs','brl'))	
	df			<- subset(df, brl==0.24 | brl==0.02 | brl==0.04)
	set(df, df[, which(brl==0.24)], 'brl', 'None')
	set(df, df[, which(brl==0.02)], 'brl', '2% subst/site')
	set(df, df[, which(brl==0.04)], 'brl', '4% subst/site')
	ggplot(df, aes(x=bs, y=100*fp.by.sum, group=brl, colour=brl)) + geom_line(size=0.8) +
			geom_hline(yintercept=5) +
			scale_y_continuous(limits=60*c(0,1), expand=c(0,0), breaks=seq(0,100,10), minor_breaks=seq(0,100,5)) +
			scale_colour_brewer(palette='Set2', name='Distance threshold') + 
			theme_bw() + labs(x='clade frequency threshold',y='probability that sequences from\nthe same individual do not co-cluster\n(%)') +			
			theme(legend.position=c(0,1), legend.justification=c(0,1), panel.grid.major=element_line(colour="grey70", size=0.4))
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150228_ClusterTypeI.pdf'
	ggsave(file=file, w=5, h=5)
	
}
######################################################################################
project.athena.Fisheretal.composition.cluster.thresholds<- function()
{
	verbose			<- 1
	resume			<- 1
	patient.n		<- 15700
	indircov		<- paste(DATA,"derived",sep='/')
	infilecov		<- "ATHENA_2013_03_AllSeqPatientCovariates"
	indir			<- paste(DATA,"tmp",sep='/')
	infile			<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500"			
	insignat		<- "Wed_Dec_18_11/37/00_2013"
	
	argv			<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	nsh.clu.pre		<- hivc.prog.get.clustering.precompute()	
	ph				<- nsh.clu.pre$ph
	ph.node.bs		<- nsh.clu.pre$ph.node.bs
	dist.brl		<- nsh.clu.pre$dist.brl.casc
	df.seqinfo		<- copy(nsh.clu.pre$df.seqinfo)
	setkey(df.seqinfo, Node)
	
	
	thresh.bs		<- 0.7
	clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, nodesupport=ph.node.bs, retval="all")
	df.seqinfo[, cluster.70I:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]	
	thresh.bs		<- 0.8
	clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, nodesupport=ph.node.bs, retval="all")
	df.seqinfo[, cluster.80I:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]	
	thresh.bs		<- 0.85
	clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, nodesupport=ph.node.bs, retval="all")
	df.seqinfo[, cluster.85I:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]		
	#thresh.bs		<- 0.9
	#clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, nodesupport=ph.node.bs, retval="all")
	#df.seqinfo[, cluster.90I:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]	
	thresh.bs		<- 0.95
	clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, nodesupport=ph.node.bs, retval="all")
	df.seqinfo[, cluster.95I:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]		
	thresh.bs		<- 0.8
	clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=0.02, dist.brl=dist.brl, nodesupport=ph.node.bs, retval="all")
	df.seqinfo[, cluster.802:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]		
	#clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=0.04, dist.brl=dist.brl, nodesupport=ph.node.bs, retval="all")
	#df.seqinfo[, cluster.804:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]	
	
		
	#
	#	plot individuals in cluster
	#
	df.plot			<- melt( df.seqinfo, id.vars=c('FASTASampleCode','Patient'), measure.vars=names(df.seqinfo)[ grepl('cluster',names(df.seqinfo)) ] )
	setkey(df.plot, variable, Patient)
	df.plot			<- unique(df.plot)
	df.plot			<- subset(df.plot, !is.na(value))[,  	{
				tmp<- table(value)
				list(CLU_ID= names(tmp), CLU_N=as.numeric(tmp))
			} , by='variable']
	
	df.plot			<- df.plot[, list(CLU_NP= length(CLU_ID)*CLU_N), by=c('variable','CLU_N')]	
	setkey(df.plot, variable, CLU_N)
	df.plot			<- merge( df.plot, subset(df.plot, CLU_N>1)[, list(CLU_N=CLU_N, CLU_NPC=cumsum(CLU_NP)), by='variable'], by=c('variable','CLU_N'), all.x=1)
	
	df.plot			<- merge(df.plot, data.table(variable=c('cluster.70I','cluster.80I','cluster.85I','cluster.95I','cluster.802'), legend=c('70%','80%','85%','95%','80% and 2%\ndistance threshold')), by='variable')
	ggplot(subset(df.plot, CLU_N>1), aes(x=CLU_N, y=CLU_NPC, group=legend, colour=legend)) + geom_line(size=0.6)  +
			scale_x_continuous(breaks=c(2,5,seq(10,100,10)), minor_breaks= c(seq(2,10,1),seq(15,100,5))) +
			scale_y_continuous(breaks=c(100, 500, seq(1000,5000,1000)), minor_breaks=seq(100, 8000,100)) +
			scale_colour_brewer(name='clustering\nthreshold', palette='Dark2') +
			theme_bw() + labs(x='phylogenetic cluster size', y='individuals in cluster\n(cumulative number)') +
			theme(legend.position=c(1,0), legend.justification=c(1,0), panel.grid.major=element_line(colour="grey70", size=0.5))
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150228_ClusterInd.pdf'
	ggsave(file=file, w=5, h=6)
	#
	#	plot prob transmitters in cluster; need YX.tpairs for correct clustering criteria .. TODO
	#
	df.plot			<- melt( df.seqinfo, id.vars=c('FASTASampleCode','Patient'), measure.vars=names(df.seqinfo)[ grepl('cluster',names(df.seqinfo)) ], value.name='CLU_ID' )
	setkey(df.plot, variable, Patient)
	df.plot			<- subset(unique(df.plot), !is.na(CLU_ID))	
	tmp				<- df.plot[,  	list(CLU_N=length(Patient)), by=c('variable','CLU_ID')]
	df.plot			<- merge(df.plot, tmp, by=c('variable','CLU_ID'))
	setkey(df.tpairs, cluster, Patient, t.Patient)
	df.plot			<- merge(df.plot, subset(unique(df.tpairs), select=c(Patient, t.Patient)), by='Patient', allow.cartesian=T)
	tmp				<- df.plot[,  	list(CLU_P=length(t.Patient)), by=c('variable','CLU_ID')]
	df.plot			<- merge(df.plot, tmp, by=c('variable','CLU_ID'))
	df.plot			<- unique(subset(df.plot, select=c(variable, CLU_ID, CLU_N, CLU_P)))
	
	# copy from above
	
}
######################################################################################
project.athena.Fisheretal.composition.cluster.brlimpact<- function()
{
	verbose			<- 1
	resume			<- 1
	patient.n		<- 15700
	indircov		<- paste(DATA,"derived",sep='/')
	infilecov		<- "ATHENA_2013_03_AllSeqPatientCovariates"
	indir			<- paste(DATA,"tmp",sep='/')
	infile			<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_examlbs500"			
	insignat		<- "Wed_Dec_18_11/37/00_2013"
	
	argv			<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
	argv			<<- unlist(strsplit(argv,' '))
	nsh.clu.pre		<- hivc.prog.get.clustering.precompute()	
	ph				<- nsh.clu.pre$ph
	ph.node.bs		<- nsh.clu.pre$ph.node.bs
	dist.brl		<- nsh.clu.pre$dist.brl.casc
	df.seqinfo		<- copy(nsh.clu.pre$df.seqinfo)
	setkey(df.seqinfo, Node)
	
	thresh.bs		<- 0.8
	clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=0.015, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	df.seqinfo[, cluster.15:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]	
	
	clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=0.045, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	df.seqinfo[, cluster.45:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]	
	
	clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=0.096, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	df.seqinfo[, cluster.80:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]
	
	clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=0.096, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	df.seqinfo[, cluster.96:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]
	
	clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=0.2, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	df.seqinfo[, cluster.200:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]	
	
	clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	df.seqinfo[, cluster.nobrl:= clustering[["clu.mem"]][seq_len(Ntip(ph))]]
	
	df.plot			<- melt( df.seqinfo, id.vars=c('FASTASampleCode','Patient'), measure.vars=names(df.seqinfo)[ grepl('cluster',names(df.seqinfo)) ] )	
	df.plot			<- subset(df.plot, !is.na(value))[,  	{
				tmp<- table(value)
				list(value= names(tmp), count=as.numeric(tmp))
			} , by='variable']
	setkey(df.plot, variable, count)
	df.plot			<- merge( df.plot, df.plot[, list(id=seq_along(value), value=value), by='variable'], by=c('variable', 'value') )
	setkey(df.plot, variable, count)
	ggplot(df.plot, aes(x=id, y=count, colour=variable, group=variable)) + geom_line()
	#	essentially no difference betw
	df.plot[, list(mx=max(count), n=length(count)), by='variable']
	#df.seqinfo[, table(cluster.96, useNA='if')]
	#ggplot( df.seqinfo, aes(x=cluster.96)) + geom_histogram(binwidth = 1)
}
######################################################################################
project.athena.Fisheretal.composition.prop.recent<- function()
{
	require(data.table)
	require(ape)
	require(RColorBrewer)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_150319",sep='/')	
	outdir					<- paste(DATA,"fisheretal_150319",sep='/')	
	
	
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	indircov				<- paste(DATA,"fisheretal_data",sep='/')
	insignat				<- "Wed_Dec_18_11:37:00_2013"	
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	t.period				<- 1/8
	t.endctime				<- hivc.db.Date2numeric(as.Date("2013-03-01"))
	t.endctime				<- floor(t.endctime) + floor( (t.endctime%%1)*100 %/% (t.period*100) ) * t.period
	
	file					<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.recent.Rdata", sep='')		
	readAttempt				<- try(suppressWarnings(load(file)))	
	#	merge tperiod legends
	tperiod.info			<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.5, 2006.5, 2008, 2009.5), t.period.max = c(2006.45, 2007.99, 2009.45, 2010.999)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))
	set(tperiod.info, NULL, 't.period.min', tperiod.info[,  paste(floor(t.period.min), floor( 1+(t.period.min%%1)*12 ), sep='-')] )
	set(tperiod.info, NULL, 't.period.max', tperiod.info[,  paste(floor(t.period.max), floor( 1+(t.period.max%%1)*12 ), sep='-')] )	
	set(tperiod.info, NULL, 't.period.min', tperiod.info[, factor(t.period.min, levels=c('1996-7','2006-7','2008-1','2009-7'), labels=c('96/07','\n\n06/07','08/01','\n\n09/07'))])
	set(tperiod.info, NULL, 't.period.max', tperiod.info[, factor(t.period.max, levels=c('2006-6','2007-12','2009-6','2010-12'), labels=c('06/06','07/12','09/06','10/12'))])
	tperiod.info[, t.period.long:= paste(gsub('\n','',t.period.min), '-', gsub('\n','',t.period.max),sep='')]
	set(tperiod.info,NULL,c('t.period.min','t.period.max'),NULL)
	tperiod.info			<- rbind(data.table(t.period='0', t.period.long='Overall'),tperiod.info)
	set(tperiod.info, NULL, 't.period.long', tperiod.info[, factor(t.period.long, levels=t.period.long, labels=t.period.long)])
	runs.UA					<- merge(runs.UA, tperiod.info, by='t.period')	
	#	merge run legends
	tmp			<- data.table(	method.brl=c(	"3pa1H1.48C2V100bInfT7", "3pa1H1.94C2V100bInfT7", "3pa1H1.09C2V100bInfT7",					
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
	set(tmp, NULL, 'method.legend', tmp[, factor(method.legend, levels=rev(method.legend), labels=rev(method.legend))])
	runs.UA		<- merge(runs.UA, tmp, by='method.brl')
	#	select
	runs.UA		<- subset(runs.UA, stat=='P.raw.e0cp')
	#	plot
	ggplot(runs.UA, aes(y=method.legend, x=100*v, xmin=100*l95.bs, xmax=100*u95.bs)) + 
			geom_point(size=2.2) +
			geom_errorbarh(height=0.5, size=0.5) +
			labs(x='\nTransmissions from men in their first year of HIV infection\n(%)', y='') +
			scale_x_continuous(breaks=seq(0,80,10), minor_breaks=seq(0,80,1)) +
			theme_bw() +
			theme(panel.margin=unit(1, "lines"), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.x=element_line(colour='grey70', size=0.5), panel.grid.minor.x=element_line(colour='grey70', size=0.2)) +
			facet_grid(~t.period.long)
	file			<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150325_PropRecent.pdf'
	ggsave(file, w=12, h=10)
}
######################################################################################
project.athena.Fisheretal.composition.seqfraction<- function()
{
	require(data.table)
	require(ape)
	require(RColorBrewer)
	#stop()
	resume					<- 1 
	indir					<- paste(DATA,"fisheretal_150319",sep='/')	
	outdir					<- paste(DATA,"fisheretal_150319",sep='/')	
	
	
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
	runs.table		<- subset(runs.table, !grepl('ARTstarted|GroupsUDA',method.risk))
	#
	outfile			<- infile	
	method.DENOM	<- 'SEQ'	
	method.RISK		<- 'm2Awmx.wtn.tp'
	method.WEIGHT	<- ''			
	method.BRL		<- '3pa1H1.48C2V100bInfT7'		
	file			<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150303_SeqFraction.pdf'
	factors			<- project.athena.Fisheretal.sensitivity.factor.legend(method.RISK)
	project.athena.Fisheretal.sensitivity.getfigures.pseq(runs.table, file, method.DENOM, method.BRL, method.RISK, method.WEIGHT, factors, tperiod.info)
}
######################################################################################
project.athena.Fisheretal.composition.testing.of.ProbTr<- function()
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
	method.recentctime		<- '2011-01-01'
	method.Acute			<- 'higher'	#'central'#'empirical'
	method.minQLowerU		<- 0.135
	method.use.AcuteSpec	<- 1
	#
	#	 I re-set inaccurate testing dates -- get df.all.allmsm with original testing dates first
	#
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	infiletree				<- paste(infile,"examlbs500",sep="_")
	insignat				<- "Wed_Dec_18_11:37:00_2013"							
	t.recent.endctime	<- hivc.db.Date2numeric(as.Date(method.recentctime))	
	t.recent.endctime	<- floor(t.recent.endctime) + floor( (t.recent.endctime%%1)*100 %/% (t.period*100) ) * t.period	
	
	adjust.AcuteByNegT			<- 1
	any.pos.grace.yr			<- Inf	
	adjust.minSCwindow			<- NA		#differs from main study
	adjust.NegTByDetectability	<- NA		#differs from main study
	
	tmp				<- project.athena.Fisheretal.select.denominator(indir, infile, insignat, indircov, infile.cov.study, infile.viro.study, infile.immu.study, infile.treatment.study, 
			infiletree=infiletree, adjust.AcuteByNegT=adjust.AcuteByNegT, adjust.NegT4Acute=NA, adjust.NegTByDetectability=adjust.NegTByDetectability, adjust.minSCwindow=adjust.minSCwindow, adjust.AcuteSelect=c('Yes','Maybe'), use.AcuteSpec=method.use.AcuteSpec, 
			t.recent.endctime=t.recent.endctime, t.recent.startctime=t.recent.startctime)	
	df.all			<- tmp$df.all	
	df.denom.CLU	<- tmp$df.select
	df.denom.SEQ	<- tmp$df.select.SEQ	
	ri.CLU			<- unique(subset(df.denom.CLU, select=Patient))
	ri.SEQ			<- unique(subset(df.denom.SEQ, select=Patient))
	df.viro			<- tmp$df.viro
	df.immu			<- tmp$df.immu
	df.treatment	<- tmp$df.treatment	
	#
	tmp					<- project.athena.Fisheretal.select.denominator(	indir, infile, insignat, indircov, infile.cov.all, infile.viro.all, infile.immu.all, infile.treatment.all, 
			infiletree=NULL, adjust.AcuteByNegT=adjust.AcuteByNegT, adjust.NegT4Acute=NA, adjust.NegTByDetectability=adjust.NegTByDetectability, adjust.minSCwindow=adjust.minSCwindow, adjust.AcuteSelect=c('Yes','Maybe'), use.AcuteSpec=method.use.AcuteSpec,
			t.recent.endctime=t.recent.endctime, t.recent.startctime=t.recent.startctime,
			df.viro.part=df.viro, df.immu.part=df.immu, df.treatment.part=df.treatment, df.all.part=df.all)	
	df.all.allmsm		<- tmp$df.all
	setkey(df.all.allmsm, Patient)
	#
	#	load all probable transmitters in TP4 and determine testing behaviour of those
	#
	indir			<- paste(DATA,"fisheretal_150319",sep='/')
	infiles			<- list.files(indir, pattern='wtn.tp1.R$')
	infiles			<- infiles[ !grepl('b0.02|b0.04',infiles) ]
	infiles			<- data.table(FILE=infiles)
	infiles[, method.brl:= infiles[, regmatches(FILE, regexpr('Yscore[^_]+',FILE))]]
	set(infiles, NULL, 'method.brl', infiles[, substring(method.brl, 7)])
	tmp				<- lapply(seq_len(nrow(infiles)), function(i)
			{
				file			<- paste(indir, '/', infiles[i,FILE], sep='')
				load(file)
				YX				<- subset(ans$YXf, select=c(t.period, t, t.Patient, Patient, score.Y, CD4b, CD4c))
				setnames( YX, 'CD4c', 'factor')	
				set(YX, YX[, which(as.numeric(t.period)>4)],'t.period','4')			
				YX[, method.brl:=infiles[i,method.brl]]
				YX
			})
	YX				<- do.call('rbind', tmp)
	#	testing of probable transmitters
	df.prt	<- subset(YX, t.period==4)[, list(t.period=t.period[1], Patient=unique(t.Patient)), by='method.brl']	
	df.prt	<- merge(df.prt, unique(df.all.allmsm), by='Patient')
	df.prt[, TestT:= AnyPos_T1-NegT]
	df.test	<- data.table(TESTL= seq(0.5,3,0.05))
	df.test	<- df.prt[, list(TESTL=df.test$TESTL, TESTLP=sapply(df.test$TESTL, function(x)  mean(!is.na(TestT) & TestT<=x))), by='method.brl']
	df.test[, DATA:='Probable transmitters\nto recipients diagnosed in 09/07-10/12']
	#	testing of recipients
	if(0)
	{
		df.rec	<- subset(YX, t.period==4)[, list(t.period=t.period[1], Patient=unique(Patient)), by='method.brl']	
		df.rec	<- merge(df.rec, unique(df.all.allmsm), by='Patient')
		df.rec[, TestT:= AnyPos_T1-NegT]
		tmp		<- data.table(TESTL= seq(0.5,3,0.05))
		tmp		<- df.rec[, list(TESTL=tmp$TESTL, TESTLP=sapply(tmp$TESTL, function(x)  mean(!is.na(TestT) & TestT<=x))), by='method.brl']
		tmp[, DATA:='Rec']
		df.test	<- rbind(df.test, tmp)		
	}
	#	testing of newly diagnosed MSM
	df.dm	<- subset(unique(df.all.allmsm), AnyPos_T1<2011 & AnyPos_T1>=2009.5 & Trm%in%c('MSM','BI')) 
	df.dm[, TestT:= AnyPos_T1-NegT]
	tmp		<- data.table(TESTL= seq(0.5,3,0.05))
	tmp		<- df.dm[, list(TESTL=tmp$TESTL, TESTLP=sapply(tmp$TESTL, function(x)  mean(!is.na(TestT) & TestT<=x)))]
	tmp[, DATA:='Diagnosed in 09/07-10/12']	
	tmp		<- merge(data.table(expand.grid(TESTL= seq(0.5,3,0.05), method.brl=df.test[, unique(method.brl)])), tmp, by='TESTL')
	df.test	<- rbind(df.test, tmp, use.names=TRUE)
	#	method legends
	tmp					<- data.table(	method.brl=c(	"3pa1H1.48C2V100bInfT7", "3pa1H1.94C2V100bInfT7", "3pa1H1.09C2V100bInfT7",					
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
	df.test				<- merge(df.test, tmp, by='method.brl')	
	
	#	plot
	ggplot(df.test, aes(x=TESTL, y=100*TESTLP, colour=DATA, group=DATA)) + geom_line() +
			facet_wrap(~method.legend, ncol=3, scales='free') +
			scale_x_continuous(breaks=seq(1,3,1), minor_breaks=seq(0.5,3,0.5), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,60,20), limits=c(0,60), expand=c(0,0)) +
			theme_bw() +
			theme(legend.position=c(1,0), legend.justification=c(1,0), panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.4)) +
			labs(x= 'time between last negative HIV test and diagnosis\n(years)', y='Proportion with a last negative test', colour='')
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150318_Testing.pdf'
	ggsave(file=file, w=10, h=10)
	#
	subset(df.test, TESTL==1 & method.brl=='3pa1H1.48C2V100bInfT7')
	
}
######################################################################################
project.athena.Fisheretal.composition.nocontact<- function()
{	
	contact.grace	<- 1.5 
	t.period		<- 0.125
	t.endctime		<- 2013.
	
	contact			<- subset( df.all.allmsm, Trm%in%c('MSM','BI'), select=c(Patient, DateLastContact, DateDied, ReasonStopRegistration))
	setkey(contact, Patient)
	contact			<- unique(contact)
	#	set last contact for NA last contact
	set(contact, contact[, which(is.na(DateDied) | DateDied>t.endctime)], 'DateDied', t.endctime)
	tmp			<- contact[, which(DateLastContact>=DateDied)]
	set(contact, tmp, 'DateLastContact', contact[tmp, DateDied])
	#	set DateDied for those that moved
	tmp			<- contact[, which(!is.na(DateLastContact) & ReasonStopRegistration=='Moved')]
	set(contact, tmp, 'DateDied', contact[tmp, DateLastContact])
	#
	set(df.treatment.allmsm, NULL, 'StopTime', hivc.db.Date2numeric(df.treatment.allmsm[,StopTime]))
	set(df.treatment.allmsm, NULL, 'StartTime', hivc.db.Date2numeric(df.treatment.allmsm[,StartTime]))
	df.treat	<- melt( df.treatment.allmsm, id.vars=c('Patient'), measure.vars=c('StartTime','StopTime'), value.name='AnyT_T' )
	setkey(df.treat, Patient, AnyT_T)
	df.treat	<- unique(df.treat)
	df.treat	<- merge( df.treat, subset(contact, select=c(Patient, DateDied)),  by='Patient', all.x=1 )
	#	keep only last AnyT_T before death. If not, adjusting DateLastContact precisely removes all not in contact periods
	df.treat	<- subset(df.treat, AnyT_T<DateDied, c(Patient, AnyT_T))
	tmp			<- merge(data.table(Patient=df.tpairs[, unique(t.Patient)]), subset(df.treat, select=c(Patient, AnyT_T)), by='Patient')
	tmp			<- tmp[, list(AnyT_TL=max(AnyT_T)), by='Patient']	
	contact		<- merge(contact, tmp, by='Patient', all.x=1)	
	tmp			<- merge(data.table(Patient=df.tpairs[, unique(t.Patient)]), subset(df.viro, select=c(Patient, PosRNA)), by='Patient')
	set(tmp, NULL, 'PosRNA', hivc.db.Date2numeric(tmp[,PosRNA]))
	tmp			<- tmp[, list(PosRNA_TL=max(PosRNA)), by='Patient']	
	contact		<- merge(contact, tmp, by='Patient', all.x=1)
	tmp			<- merge(data.table(Patient=df.tpairs[, unique(t.Patient)]), subset(df.immu, select=c(Patient, PosCD4)), by='Patient')
	set(tmp, NULL, 'PosCD4', hivc.db.Date2numeric(tmp[,PosCD4]))
	tmp			<- tmp[, list(PosCD4_TL=max(PosCD4)), by='Patient']	
	contact		<- merge(contact, tmp, by='Patient', all.x=1)		
	tmp			<- contact[, which(DateLastContact<PosRNA_TL & PosRNA_TL<t.endctime)]
	if(length(tmp))
	{
		cat(paste('\nsetting relaxed DateLastContact because PosRNA_TL for n=',length(tmp)))
		set(contact, tmp, 'DateLastContact', contact[tmp,PosRNA_TL])
	}
	tmp			<- contact[, which(DateLastContact<PosCD4_TL & PosCD4_TL<t.endctime)]
	if(length(tmp))
	{
		cat(paste('\nsetting relaxed DateLastContact because PosCD4_TL for n=',length(tmp)))
		set(contact, tmp, 'DateLastContact', contact[tmp,PosCD4_TL])
	}
	tmp			<- contact[, which(DateLastContact<AnyT_TL & AnyT_TL<t.endctime)]	
	if(length(tmp))
	{
		cat(paste('\nsetting relaxed DateLastContact because AnyT_TL for n=',length(tmp)))
		set(contact, tmp, 'DateLastContact', contact[tmp,AnyT_TL])
	}	
	#	allow for grace at end
	tmp			<- contact[, which(DateLastContact<t.endctime & DateLastContact+contact.grace>=t.endctime)]
	if(length(tmp))
	{
		cat(paste('\nsetting relaxed DateLastContact because DateLastContact+contact.grace<t.endctime for n=',length(tmp)))
		set(contact, tmp, 'DateLastContact', contact[tmp, t.endctime])
	}
	#	set No contact at end
	set(contact, NULL, 'DateLastContact', contact[, floor(DateLastContact) + ceiling( (DateLastContact%%1)*100 %/% (t.period*100) ) * t.period] )
	set(contact, NULL, 'DateDied', contact[, floor(DateDied) + round( (DateDied%%1)*100 %/% (t.period*100) ) * t.period] )
	tmp			<- subset(contact, DateLastContact<DateDied)	
	tmp			<- tmp[, list(t= seq(DateLastContact, DateDied-t.period, by=t.period), contact='No'),by='Patient']	
	setnames(tmp, 'Patient','t.Patient')
	X.incare	<- merge(X.incare, tmp, by=c('t.Patient','t'), all.x=1)
	set(X.incare, X.incare[,which(is.na(contact))],'contact','Yes')
	set(X.incare, X.incare[, which(contact=='No' & t+contact.grace>t.endctime)], 'contact', 'Yes')
	#
	#	for each time t, check if there is at least 1 VL or CD4 measurement or Treatment visit in the last year/2 or next year/2
	#
	set(df.immu, NULL, 'PosCD4', hivc.db.Date2numeric(df.immu[,PosCD4]))
	set(df.viro, NULL, 'PosRNA', hivc.db.Date2numeric(df.viro[,PosRNA]))
	tmp			<- X.incare[, {
				cntct					<- rep(0, length(t))
				endgrace				<- t+contact.grace-t.endctime
				endgrace[endgrace<0]	<- 0									
				z						<- df.immu$PosCD4[ which(df.immu$Patient==t.Patient) ]
				if(length(z))
				{
					tmp			<- sapply(seq_along(t), 	function(i) min(abs(t[i]-z))<=(contact.grace/2+endgrace[i])		)
					cntct[tmp]	<- 1
				}									
				z		<- df.viro$PosRNA[ which(df.viro$Patient==t.Patient) ]
				if(!all(cntct==1) & length(z))
				{
					tmp			<- sapply(seq_along(t), 	function(i) min(abs(t[i]-z))<=(contact.grace/2+endgrace[i])		)
					cntct[tmp]	<- 1
				}
				z						<- df.treat$AnyT_T[ which(df.treat$Patient==t.Patient) ]
				if(length(z))
				{
					tmp			<- sapply(seq_along(t), 	function(i) min(abs(t[i]-z))<=(contact.grace/2+endgrace[i])		)
					cntct[tmp]	<- 1
				}																		
				list(t=t, incontact=cntct)
			}, by=c('t.Patient')]	
	X.incare	<- merge(X.incare, tmp, by=c('t.Patient','t'))	
	if(!is.null(plot.file))
	{
		tmp			<- subset(X.incare, t>1996.5 & t<2011 & (contact=='No' | incontact==0))[, list(FIR=t[1], DUR=length(t)*t.period), by=c('t.Patient')]
		setnames(tmp, 't.Patient', 'Patient')
		tmp2		<- subset(df.all, select=c(Patient, PosSeqT))
		setkey(tmp2, Patient)
		tmp			<- merge(tmp, unique(tmp2), by='Patient', all.x=1)
		ggplot(tmp, aes(x=FIR, y=DUR)) + geom_jitter(size=1) +
				scale_x_continuous(breaks=seq(1980, 2020, 5), minor_breaks=seq(1980, 2020, 1)) +
				scale_y_continuous(breaks=seq(0, 20, 4), minor_breaks=seq(0, 20, 1), expand=c(0.01,0.01)) +
				labs(x='first time not in contact', y='duration of loss to contact\n(years)') +
				theme(panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.4), axis.text.x=element_text(angle=0, vjust=0, hjust=0)) +
				theme_bw()
		file		<- paste(plot.file, '_contactDur.pdf',sep='')
		ggsave(file=file, h=6, w=6)
		#
		tmp			<- subset(tmp, !is.na(PosSeqT))
		set(tmp, NULL, 'PosSeqT', tmp[, floor(PosSeqT) + round( (PosSeqT%%1)*100 %/% (t.period*100) ) * t.period] )
		tmp[, COL:=tmp[, factor(as.numeric(FIR<=PosSeqT), levels=c(0,1), labels=c('before loss of contact', 'after loss of contact'))]]		
		ggplot(tmp, aes(x=FIR, y=PosSeqT, colour=COL)) + geom_point() + geom_abline() +
				labs(colour='sequenced', x='first time not in contact', y='sequence sampling time') +
				theme_bw() +
				theme(legend.position='bottom')
		file		<- paste(plot.file, '_contactBySeq.pdf',sep='')
		ggsave(file=file, h=6, w=6)		
		#ggplot(subset(X.incare, contact=='No' | incontact==0), aes(x=t, fill=incontact==0)) + geom_histogram(binwidth=0.125) + scale_x_continuous(breaks=seq(1980, 2020, 5), minor_breaks=seq(1980, 2020, 1)) + coord_trans(limx=c(1996.5, 2011))
		tmp			<- subset(X.incare, contact=='No' | incontact==0)
		setnames(tmp, 't.Patient', 'Patient')
		tmp			<- merge(tmp, unique(tmp2), by='Patient', all.x=1)
		ggplot(tmp, aes(x=t, fill=factor(as.numeric(!is.na(PosSeqT)), levels=c(0,1),labels=c("No","Yes")))) + geom_histogram(binwidth=0.125) + 
				scale_x_continuous(breaks=seq(1980, 2020, 5), minor_breaks=seq(1980, 2020, 1)) + 
				scale_y_continuous(expand=c(0,0)) +
				coord_trans(limx=c(1996.5, 2011)) +
				labs(x='', y='potential transmitters with no contact\n(#)', fill='with a sequence') +
				theme_bw() +
				theme(panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.4), axis.text.x=element_text(angle=0, vjust=0, hjust=0), legend.position='bottom')
		file		<- paste(plot.file, '_contactHist.pdf',sep='')
		ggsave(file=file, h=6, w=6)		
		tmp			<- subset(X.incare, t>1996 & t<2011.5)[, list(NOCON= mean(contact=='No' | incontact==0)), by='t']
		ggplot(tmp, aes(x=t, ymax=100*NOCON, ymin=0)) + geom_ribbon() + 
				scale_x_continuous(breaks=seq(1980, 2020, 5), minor_breaks=seq(1980, 2020, 1)) + 
				scale_y_continuous(expand=c(0,0), limits=c(0,100)) +
				coord_trans(limx=c(1996.5, 2011)) +
				theme(panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.4), axis.text.x=element_text(angle=0, vjust=0, hjust=0)) +
				labs(x='', y='potential transmitters with no contact\n(%)') +
				theme_bw()
		file		<- paste(plot.file, '_contactProp.pdf',sep='')
		ggsave(file=file, h=6, w=6)
		
	}
	#
	set(X.incare, X.incare[, which(!incontact & contact=='Yes')], 'contact', 'No')	
	set(X.incare, NULL, 'incontact', NULL )
	if(X.incare[, length(which(is.na(stage)))])	stop('unexpected NA stage')
	X.incare
}
######################################################################################
project.athena.Fisheretal.composition.viralload.of.treated<- function()
{
	dfv		<- subset( df.viro.allmsm, select=c(Patient, PosRNA, lRNA) )
	set(dfv, NULL, 'PosRNA', hivc.db.Date2numeric(dfv[,PosRNA]))
	
	df.cov	<- subset( df.all.allmsm, Trm%in%c('MSM','BI'))
	setkey(df.cov, Patient)
	df.cov	<- unique(df.cov)
	tmp		<- subset( df.cov, AnyT_T1<2013 & DateLastContact>2013 & (is.na(DateDied) | DateDied>2013))		
	tmp		<- merge( subset(dfv, PosRNA>=2012 & PosRNA<2013), subset(tmp, select=c(Patient, AnyPos_T1, AnyT_T1)), by='Patient' )
	tmp		<- subset(tmp, PosRNA>AnyT_T1)
	setkey(tmp, Patient, PosRNA)
	tmp		<- tmp[, {
				z	<- which.max(PosRNA)
				list(lRNA= lRNA[z])
			}, by='Patient']
	tmp[, mean(lRNA<log10(100))]
	#	0.9459634
	
	tmp		<- subset( df.cov, AnyT_T1<2004 & DateLastContact>2004 & (is.na(DateDied) | DateDied>2004))		
	tmp		<- merge( subset(dfv, PosRNA>=2003 & PosRNA<2004), subset(tmp, select=c(Patient, AnyPos_T1, AnyT_T1)), by='Patient' )
	tmp		<- subset(tmp, PosRNA>AnyT_T1)
	setkey(tmp, Patient, PosRNA)
	tmp		<- tmp[, {
				z	<- which.max(PosRNA)
				list(lRNA= lRNA[z])
			}, by='Patient']
	tmp[, mean(lRNA<log10(100))]
	#	0.796652
}
######################################################################################
project.athena.Fisheretal.composition.viralload<- function()
{
	dfv		<- subset( df.viro.allmsm, select=c(Patient, PosRNA, lRNA) )
	set(dfv, NULL, 'PosRNA', hivc.db.Date2numeric(dfv[,PosRNA]))
	
	dfs		<- data.table(Patient=c('M10724','M30724','M30740','M30783','M30799'))
	dfv		<- merge(dfs, dfv, by='Patient')
	dfv[,lRNAi:= lRNA]
	
	#dfid	<- merge(dfs, immu, by='Patient')
	dfvp	<- merge(dfs, df.all.allmsm, by='Patient')
	set(dfvp, dfvp[, which(is.na(DateDied))],'DateDied', 2013.3)
	#	I think this is not too bad at all
	#	only difficulty is for patients with large gaps in surveillance
	#	replace times with gap>1 yr with NA 
	ggplot(dfis, aes(group=Patient)) + 
			geom_rect(data=dfvp, aes(xmin=AnyPos_T1, xmax=AnyT_T1, ymin=0, ymax=10), fill='orange', alpha=0.1) +
			geom_rect(data=dfvp, aes(xmin=AnyT_T1, xmax=DateDied, ymin=0, ymax=10), fill='green', alpha=0.1) +
			geom_hline(yintercept=log10(100), colour='blue', alpha=0.3) +
			geom_point(data=dfv, aes(x=PosRNA, y=lRNA), colour='grey50') +			
			geom_line(data=dfv, aes(x=PosRNA, y=lRNA, group=Patient), colour='black') + 				
			scale_y_continuous(breaks=seq(0, 6, 2), expand=c(0,0)) +				
			scale_x_continuous(breaks=seq(1995, 2015, 5), minor_breaks=seq(1995, 2015, 1), expand=c(0,0)) +
			#coord_trans(limx=c(1996, 2013), limy=c(0,1500)) + 
			theme_bw() +
			facet_wrap(~Patient, scales='free') + labs(x='', y='log10 viral load\n(cps/ml)') +
			theme(strip.text=element_blank(), strip.background=element_blank())
	
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150303_VLtraj.pdf'
	ggsave(file=file, w=10, h=7)
	
}
######################################################################################
project.athena.Fisheretal.composition.CD4atDiagnosis<- function()
{
	t.period				<- 0.125
	t.recent.startctime		<- hivc.db.Date2numeric(as.Date("1996-07-15"))
	t.recent.startctime		<- floor(t.recent.startctime) + floor( (t.recent.startctime%%1)*100 %/% (t.period*100) ) * t.period
	t.recent.endctime		<- 2011
	#tp.cut					<- c(-Inf, 2006.5, 2008, 2009.5, 2011)
	tperiod.info			<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.5, 2006.5, 2008, 2009.5)+0.001, t.period.max = c(2006.5, 2008, 2009.5, 2011)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))		
	file					<- paste(DATA,'/derived/','ATHENA_2014_06_Patient_AllMSM_CD4.R',sep='')
	load(file)		
		
	df.cov	<- subset( df.all.allmsm, Trm%in%c('MSM','BI'))
	setkey(df.cov, Patient)
	df.cov	<- unique(df.cov)
	tmp		<- subset( df.cov, AnyPos_T1<2004 & DateLastContact>2004 & (is.na(DateDied) | DateDied>2004))	
	tmp[, mean(!is.na(AnyT_T1) & AnyT_T1<2004)]
	#	0.7852525
	tmp		<- subset(df.cov, !is.na(AnyT_T1) & AnyT_T1<2004 & AnyT_T1>=2003, select=c(Patient, AnyT_T1))	
	tmp		<- merge(immu.sm, tmp, by='Patient')
	tmp		<- tmp[, {
			z	<- which.min(abs(t-AnyT_T1))
			list(t=t[z], CD4=CD4[z], AnyT_T1=AnyT_T1[z])
			}, by='Patient']
	subset(tmp, !is.na(CD4))[, median(CD4)]
	#	238.64
	tmp		<- subset(df.cov, !is.na(AnyT_T1) & AnyT_T1<2004 & AnyT_T1>=2003, select=c(Patient, AnyT_T1))
	tmp		<- merge(df.immu.allmsm, tmp, by='Patient')
	set(tmp, NULL, 'PosCD4', tmp[, hivc.db.Date2numeric(PosCD4)])
	tmp		<- tmp[, {
				z	<- tail(which(AnyT_T1>=PosCD4),1)
				list(PosCD4=PosCD4[z], CD4=CD4[z], AnyT_T1=AnyT_T1[z])
			}, by='Patient']
	subset(tmp, !is.na(CD4))[, median(CD4)]
	#	200
	
	
	
	tmp		<- subset( df.cov, AnyPos_T1<2013 & (is.na(DateLastContact) | DateLastContact>2013) & (is.na(DateDied) | DateDied>2013))
	nrow(tmp)
	#	9799
	tmp[, mean(!is.na(AnyT_T1) & AnyT_T1<2013)]
	#	0.8543729	??
	tmp		<- subset(df.cov, !is.na(AnyT_T1) & AnyT_T1<2013 & AnyT_T1>=2012, select=c(Patient, AnyT_T1))	
	tmp		<- merge(immu.sm, tmp, by='Patient')
	tmp		<- tmp[, {
				z	<- tail(which(AnyT_T1>=t),1)
				list(t=t[z], CD4=CD4[z], AnyT_T1=AnyT_T1[z])
			}, by='Patient']
	subset(tmp, !is.na(CD4))[, median(CD4)]
	#	390
	tmp		<- subset(df.cov, !is.na(AnyT_T1) & AnyT_T1<2013 & AnyT_T1>=2012, select=c(Patient, AnyT_T1))
	nrow(tmp)
	#	499
	tmp		<- merge(df.immu.allmsm, tmp, by='Patient')
	set(tmp, NULL, 'PosCD4', tmp[, hivc.db.Date2numeric(PosCD4)])
	setkey(tmp, Patient, PosCD4)
	tmp		<- tmp[, {
						z	<- tail(which(AnyT_T1>=PosCD4),1)
						list(PosCD4=PosCD4[z], CD4=CD4[z], AnyT_T1=AnyT_T1[z])
					}, by='Patient']
	subset(tmp, !is.na(CD4))[, median(CD4)]
	#	354.5
	tmp		<- subset(df.cov, !is.na(AnyT_T1) & AnyT_T1<2013 & AnyT_T1>=2012, select=c(Patient, AnyT_T1))
	tmp		<- merge(df.immu.allmsm, tmp, by='Patient')
	set(tmp, NULL, 'PosCD4', tmp[, hivc.db.Date2numeric(PosCD4)])
	setkey(tmp, Patient, PosCD4)
	tmp		<- tmp[, {
				z	<- which(AnyT_T1>=PosCD4 & AnyT_T1-PosCD4<0.25)
				list(CD4=mean(CD4[z]), AnyT_T1=AnyT_T1[z[1]])
			}, by='Patient']
	subset(tmp, !is.nan(CD4))[, median(CD4)]
	#	347	
}
######################################################################################
project.athena.Fisheretal.composition.DiagCD4350<- function()
{
	t.period				<- 0.125
	t.recent.startctime		<- hivc.db.Date2numeric(as.Date("1996-07-15"))
	t.recent.startctime		<- floor(t.recent.startctime) + floor( (t.recent.startctime%%1)*100 %/% (t.period*100) ) * t.period
	t.recent.endctime		<- 2011
	#tp.cut				<- c(-Inf, 2006.5, 2008, 2009.5, 2011)
	tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.5, 2006.5, 2008, 2009.5)+0.001, t.period.max = c(2006.5, 2008, 2009.5, 2011)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))		
	
	df.cov	<- copy(df.all.allmsm)
	setkey(df.cov, Patient)
	df.cov	<- unique(df.cov)
	#
	# get first time CD4<350
	#
	file	<- paste(DATA,'/derived/','ATHENA_2014_06_Patient_AllMSM_CD4.R',sep='')
	load(file)		
	tmp		<- immu.sm[, {
				CD4_350_T1	<- which(CD4<=350)
				CD4_500_T1	<- which(CD4<=500)
				list(CD4sm='Y', CD4_350_T1=ifelse(length(CD4_350_T1), t, NA_real_), CD4_500_T1=ifelse(length(CD4_500_T1), t, NA_real_))
			}, by='Patient']
	df.cov	<- merge(df.cov, tmp, all.x=TRUE, by='Patient')
	df.cov	<- subset(df.cov, Trm=='MSM' | Trm=='BI')
	cat(paste('\nNumber of patients for whom we don t have a smooth, n=', subset(df.cov, is.na(CD4sm))[, length(unique(Patient))] ))
	#	for those for whom we have one CD4 count, set according to that
	tmp		<- df.cov[, which(is.na(CD4sm) & !is.na(CD4_T1) & CD4_T1<=500)]	
	set( df.cov, tmp, 'CD4_500_T1', df.cov[tmp, PosCD4_T1] )
	tmp		<- df.cov[, which(is.na(CD4sm) & !is.na(CD4_T1) & CD4_T1<=350)]	
	set( df.cov, tmp, 'CD4_350_T1', df.cov[tmp, PosCD4_T1] )
	#	for all remaining ones without a smooth, we don t have a single CD4 count
	cat(paste('\nNumber of patients for whom we don t have a single CD4 count, n=', nrow(unique(subset(df.cov, is.na(CD4_T1), Patient))) ))
	#	13 patients missing that have 1 CD4 count and not interpolated
	#	1044 patients without CD4		
	df.cd4cov	<- data.table(t= seq(t.recent.startctime, t.recent.endctime, 1/12))
	df.cd4cov	<- df.cd4cov[, list(	N_CD4_350=		nrow(subset(df.cov, t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & t>=CD4_350_T1 )),
					N_CD4_l500=		nrow(subset(df.cov, t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & t>=CD4_500_T1 & (is.na(CD4_350_T1) |  t<CD4_350_T1))),
					N_CD4_g500=		nrow(subset(df.cov, t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & t>=PosCD4_T1 & (t<CD4_500_T1 | is.na(CD4_500_T1)))),
					N_CD4_NA= 		nrow(subset(df.cov, t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & (t<PosCD4_T1 | is.na(PosCD4_T1)) )),
					N_CD4_350_ART=	nrow(subset(df.cov, t>=AnyT_T1 & t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & t>=CD4_350_T1)),	
					N_CD4_l500_ART=	nrow(subset(df.cov, t>=AnyT_T1 & t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & t>=CD4_500_T1 & (is.na(CD4_350_T1) |  t<CD4_350_T1))),
					N_CD4_g500_ART=	nrow(subset(df.cov, t>=AnyT_T1 & t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & t>=PosCD4_T1 &  (t<CD4_500_T1 | is.na(CD4_500_T1)))),
					N_CD4_NA_ART= 	nrow(subset(df.cov, t>=AnyT_T1 & t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & (t<PosCD4_T1 | is.na(PosCD4_T1)) ))										
			), by='t']
	df.cd4cov[, N_CD4_350_NoART:= N_CD4_350-N_CD4_350_ART]
	df.cd4cov[, N_CD4_l500_NoART:= N_CD4_l500-N_CD4_l500_ART]
	df.cd4cov[, N_CD4_g500_NoART:= N_CD4_g500-N_CD4_g500_ART]
	df.cd4cov[, N_CD4_NA_NoART:= N_CD4_NA-N_CD4_NA_ART]
	df.cd4cov[, N_NoART:=N_CD4_350_NoART+N_CD4_l500_NoART+N_CD4_g500_NoART+N_CD4_NA_NoART]
	set(df.cd4cov, NULL, 't.period', df.cd4cov[, cut(t, breaks=c(-Inf, 2006.5, 2008, 2009.5, 2011), labels=1:4)])
	#
	#
	#
	method.RISK			<- 'm2Awmx.wtn.tp'
	indir				<- paste(DATA,"fisheretal_150319",sep='/')
	infile				<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	insignat			<- "Wed_Dec_18_11:37:00_2013"		
	file				<- paste(outdir, '/', infile, '_', gsub('/',':',insignat), '_', "method.risks.Rdata", sep='')		
	readAttempt			<- try(suppressWarnings(load(file)))	
	df					<- subset(runs.risk, is.na(t.period) & method.nodectime=='any' & method.denom=='SEQ' & method.recentctime==2011 & grepl(method.RISK, method.risk) & method.dating=='sasky' & stat=='P.raw.e0cp' & grepl('stageDtl350', coef))
	tmp					<- data.table(	method.brl=c(	"3pa1H1.48C2V100bInfT7", "3pa1H1.94C2V100bInfT7", "3pa1H1.09C2V100bInfT7",					
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
	df					<- merge(df, tmp, by='method.brl')	
	for(method.BRL in df[,unique(method.brl)])
	{
		df.Praw.e0cp		<- subset(df, method.brl==method.BRL)
		set(df.Praw.e0cp, NULL, 't.period', df.Praw.e0cp[, as.character(factor)])
		set(df.Praw.e0cp, NULL, 't.period', df.Praw.e0cp[, factor(substr(t.period, nchar(t.period), nchar(t.period)))])
		#df.cd4cov	<- subset(df.cd4cov, select=c(t, N_CD4_350_NoART))
		df.cd4me	<- tperiod.info[, list(M_CD4_350_NoART=subset(df.cd4cov, t>=t.period.min & t<t.period.max)[, mean(N_CD4_350_NoART)]), by=t.period]
		df.cd4me	<- merge(tperiod.info, df.cd4me, by='t.period')
		df.cd4me	<- merge(df.cd4me, df.Praw.e0cp, by='t.period')
		tmp			<- lm(M_CD4_350_NoART~v-1, data=df.cd4me)
		ans			<- lm(M_CD4_350_NoART~v, data=df.cd4me)
		set(df.cd4me, NULL, 'tv', predict(tmp))
		set(df.cd4me, NULL, 'tl95.bs', predict(tmp, data.frame(v= df.cd4me$l95.bs)))
		set(df.cd4me, NULL, 'tu95.bs', predict(tmp, data.frame(v= df.cd4me$u95.bs)))
		set(df.cd4cov, NULL, 't.period', df.cd4cov[, cut(t, breaks=c(-Inf, 2006.5, 2008, 2009.5, 2011), labels=1:4)])
		df.cd4cov[, method.legend:= df.Praw.e0cp[1, method.legend] ]
		#
		library(gtable)
		library(grid)
		p1	<- ggplot(df.cd4cov) + geom_step(aes(x=t, y=N_CD4_350_NoART), stat='identity', colour="#41B6C4", width=1/12, alpha=0.4) + theme_bw() +
				scale_x_continuous(breaks=seq(1997,2014,2), minor_breaks=NULL, expand=c(0,0)) + 
				scale_y_continuous(breaks=seq(0,3000,250)) +
				geom_segment(data=df.cd4me, aes(x=t.period.min, xend=t.period.max, y=M_CD4_350_NoART, yend=M_CD4_350_NoART), col="#41B6C4", size=1.4) +
				geom_point(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, y=tv), col='black', shape=18, size=4) +				
				geom_errorbar(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, ymin=tl95.bs, ymax=tu95.bs), col='black', width=0.5, size=1) +
				labs(x='', y='Untreated HIV infected MSM\nwith CD4<350\n( # )', title=df.Praw.e0cp[1, method.legend]) + 
				facet_grid(.~t.period, scales='free_x', space='free_x') +					
				theme(strip.background=element_blank(), strip.text = element_blank(), axis.text.y=element_text(colour="#41B6C4"), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x=element_blank())
		p2	<- ggplot(df.cd4cov) + geom_step(aes(x=t, y=N_CD4_350_NoART), stat='identity', colour="#41B6C4", width=1/12, alpha=0.4) + theme_bw() +
				scale_x_continuous(breaks=seq(1997,2014,2), minor_breaks=NULL, expand=c(0,0)) + 
				scale_y_continuous(breaks=predict(tmp, data.frame(v= seq(0, 14, 2)/100)), label=seq(0, 14, 2)) +					
				geom_segment(data=df.cd4me, aes(x=t.period.min, xend=t.period.max, y=M_CD4_350_NoART, yend=M_CD4_350_NoART), col="#41B6C4", size=1.4) +
				geom_point(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, y=tv), col='black', shape=18, size=4) +				
				geom_errorbar(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, ymin=tl95.bs, ymax=tu95.bs), col='black', width=0.5, size=1) +
				labs(x='', y='Transmissions from\nuntreated men with CD4<350\n( % )', title=df.Praw.e0cp[1, method.legend]) +
				facet_grid(.~t.period, scales='free_x', space='free_x') +					
				theme(strip.background=element_blank(), strip.text = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x=element_blank())				
		#	prepare y axis label 	
		p2	<- p2 + theme(axis.text.y = element_text(hjust = 0), axis.title.y = element_text(angle = 270)) 
		#extract gtable
		g1	<- ggplot_gtable(ggplot_build(p1))
		g2	<- ggplot_gtable(ggplot_build(p2))		
		#overlap the panel of the 2nd plot on that of the 1st plot		
		pp	<- c(subset(g1$layout, name=="panel", se=t:r)[1, ])
		g	<- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")[1]]], pp$t, pp$l, pp$b, pp$l)
		#extract left axis
		ia <- which(g2$layout$name == "axis-l")
		ga <- g2$grobs[[ia]]
		ax <- ga$children[[2]]
		ax$widths <- rev(ax$widths)
		ax$grobs <- rev(ax$grobs)
		ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
		g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
		g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
		#extract left axis label		
		g 	<- gtable_add_grob(g, g2$grobs[[which(g2$layout$name == "ylab")]], pp$t, length(g$widths), pp$b)
		#draw the whole thing	
		file	<- paste('/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150316_DiagCD4350_',df.Praw.e0cp[1,method.brl],'.pdf', sep='')		
		pdf(file=file, w=6, h=4) # plot saved by default to Rplots.pdf
		grid.newpage()
		grid.draw(g)
		dev.off() 	
	}
}
######################################################################################
project.athena.Fisheretal.composition.CD4model<- function()
{
	if(0)
	{
		#	model CD4 counts
		t.period<- 0.125
		df.cov	<- copy(df.all.allmsm)
		setkey(df.cov, Patient)
		df.cov	<- unique(df.cov)		
		immu	<- subset( df.immu.allmsm, select=c(Patient, PosCD4, CD4) )
		set(immu, NULL, 'PosCD4', hivc.db.Date2numeric(immu[,PosCD4]))
		tmp		<- subset(df.cov, select=c(Patient, AnyT_T1))
		set(tmp, tmp[, which(is.na(AnyT_T1))], 'AnyT_T1', 2030.)
		tmp2	<- immu[, list(CD4.gap=ifelse(length(PosCD4)==1, 0, max(diff(PosCD4)))), by='Patient']
		tmp		<- merge(tmp, tmp2, by='Patient')		
		immu	<- merge(immu, tmp, by='Patient')	
		require(gamlss)
		stopifnot(class(immu$Patient)=='character',class(df.cov$Patient)=='character')
		#	add time to next PosCD4
		tmp		<- immu[, {
					z	<- c(diff(PosCD4),0)
					if(length(PosCD4)==1)
						z<- 0
					list(PosCD4=PosCD4, PosCD4d=z)	
				}, by='Patient']
		immu	<- merge(immu, tmp, by=c('Patient','PosCD4'))
		#			
		tmp		<- subset(immu, select=c(Patient, PosCD4))[, list(ts=min(PosCD4), te=max(PosCD4)), by='Patient']
		set(tmp, NULL, 'ts', tmp[, floor(ts) + round( (ts%%1)*100 %/% (t.period*100) ) * t.period] )
		set(tmp, NULL, 'te', tmp[, floor(te) + round( (te%%1)*100 %/% (t.period*100) ) * t.period] )
		tmp		<- tmp[, list(PosCD4= seq(ts, te, by=t.period)),by='Patient']	
		setkey(tmp, Patient)
		immu.pa	<- unique(subset(immu, select=Patient))
		immu.pa[, BATCH:=floor(seq_len(nrow(immu.pa))/100)]
		#	do 100 at a time
		immu.sm<- lapply( immu.pa[, unique(BATCH)], function(b)
				{
					cat(paste('\nprocess Batch',b))
					immu.sm	<- lapply( subset( immu.pa, BATCH==b)[, unique(Patient)], function(x)
							{								
								ART.st	<- subset(immu, Patient==x)[1, AnyT_T1]
								ART.st	<- ifelse(is.na(ART.st), 2030., ART.st)	
								#	before ART start
								z		<- subset(tmp, Patient==x & PosCD4<=ART.st, PosCD4)		#	times to predict at
								z2		<- subset(immu, Patient==x & PosCD4<=ART.st)			#	data to build model
								if(!nrow(z2) || !nrow(z))
									ans	<- data.table(Patient=x, t=NA_real_, CD4=NA_real_)
								else
								{
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
											tryCatch({
														cd4.ml	<- gamlss(CD4 ~ PosCD4, data=z2, family='NO', trace = FALSE)
														cd4.m	<- gamlss(CD4 ~ bs(PosCD4, degree=cd4.d), data=z2, family='NO', trace = FALSE)										
														#	gamlss fit may go wild occasionally, in this case fall back to linear interpolation		
														if(deviance(cd4.ml)-deviance(cd4.m)>10)
															cd4.s	<- predict(cd4.m, type='response', newdata=z, data=z2) 
														if(deviance(cd4.ml)-deviance(cd4.m)<=10)
															cd4.s	<- predict(cd4.ml, type='response', newdata=z, data=z2) 	
													}, 
													warning=function(w)
													{ 
														cat(paste('\nWarning: fall back to approx for patient',x))
														cd4.s	<<- approx(z2[,PosCD4], z2[,CD4], xout=z[,PosCD4]+t.period/2, rule=2)$y	
													}, 
													error=function(e)
													{ 
														cat(paste('\nError: fall back to approx for patient',x))
														cd4.s	<<- approx(z2[,PosCD4], z2[,CD4], xout=z[,PosCD4]+t.period/2, rule=2)$y	
													})						
										}
									}
									ans	<- data.table(Patient=x, t=z[,PosCD4], CD4=cd4.s)	#answer
									z2	<- subset(z2, PosCD4d>2)				
									for(i in seq_len(nrow(z2)))
									{
										set(ans, ans[,which( t>z2[i, PosCD4] & t<z2[i, PosCD4+PosCD4d])], 'CD4', NA_real_)
									}
								}							
								#	after ART start
								z		<- subset(tmp, Patient==x & PosCD4>ART.st, PosCD4)		#	times to predict at
								z2		<- subset(immu, Patient==x & PosCD4>ART.st)				#	data to build model
								if(nrow(z) && nrow(z2))								
								{
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
											tryCatch({
														cd4.ml	<- gamlss(CD4 ~ PosCD4, data=z2, family='NO', trace = FALSE)
														cd4.m	<- gamlss(CD4 ~ bs(PosCD4, degree=cd4.d), data=z2, family='NO', trace = FALSE)													
														#	gamlss fit may go wild occasionally, in this case fall back to linear interpolation		
														if(deviance(cd4.ml)-deviance(cd4.m)>10)
															cd4.s	<- predict(cd4.m, type='response', newdata=z, data=z2) 
														if(deviance(cd4.ml)-deviance(cd4.m)<=10)
															cd4.s	<- predict(cd4.ml, type='response', newdata=z, data=z2) 	
													}, 
													warning=function(w)
													{ 
														cat(paste('\nWarning: fall back to approx for patient',x))
														cd4.s	<<- approx(z2[,PosCD4], z2[,CD4], xout=z[,PosCD4]+t.period/2, rule=2)$y	
													},
													error=function(e)
													{ 
														cat(paste('\nError: fall back to approx for patient',x))
														cd4.s	<<- approx(z2[,PosCD4], z2[,CD4], xout=z[,PosCD4]+t.period/2, rule=2)$y	
													})						
										}
									}
									z3	<- data.table(Patient=x, t=z[,PosCD4], CD4=cd4.s)	#answer
									z2	<- subset(z2, PosCD4d>2)				
									for(i in seq_len(nrow(z2)))
									{
										set(z3, z3[,which( t>z2[i, PosCD4] & t<z2[i, PosCD4+PosCD4d])], 'CD4', NA_real_)
									}
									ans	<- subset(rbind(ans, z3), !is.na(t))
								}
								ans				
							})
					immu.sm	<- do.call('rbind',immu.sm)
					save(file=paste(outdir,'/','ATHENA_composition_CD4_batch',b,'.R',sep=''),immu.sm)
					immu.sm				
				})
		immu.sm	<- do.call('rbind',immu.sm)
		immu.sm	<- subset(immu.sm, !is.na(t))
		
		file	<- paste(DATA,'/derived/','ATHENA_2014_06_Patient_AllMSM_CD4.R',sep='')
		save(file=paste(outdir,'/','ATHENA_2014_06_Patient_AllMSM_CD4.R',sep=''),immu.sm)
	}
	if(0)
	{
		file			<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_data/ATHENA_2014_06_Patient_AllMSM_CD4.R"
		load(file)
		#	select a few patients
		immu.sm			<- subset(immu.sm, !is.na(t))
		immu	<- subset( df.immu.allmsm, select=c(Patient, PosCD4, CD4) )
		set(immu, NULL, 'PosCD4', hivc.db.Date2numeric(immu[,PosCD4]))
		
		dfs		<- immu.sm[, list(CD4N=length(t), CD41=min(t)), by='Patient']
		#subset(dfs, CD41>2005)
		#dfs		<- data.table(Patient=c('M10559','M10694','M10710','M10724','M10758','M10790','M17184','M30538','M36429'))
		#dfs		<- data.table(Patient=c('M10724','M10758','M17203','M17066','M17081','M17176','M17184','M17189','M17209','M17217','M17220','M17221','M17223','M10790','M17184','M30538','M36429'))
		dfs		<- data.table(Patient=c('M10724','M30724','M30740','M30783','M30799'))
		dfis	<- merge(dfs, immu.sm, by='Patient')
		dfid	<- merge(dfs, immu, by='Patient')
		dfip	<- merge(dfs, df.all.allmsm, by='Patient')
		set(dfip, dfip[, which(is.na(DateDied))],'DateDied', 2013.3) 
		#	I think this is not too bad at all
		#	only difficulty is for patients with large gaps in surveillance
		#	replace times with gap>1 yr with NA 
		ggplot(dfis, aes(group=Patient)) + 
				geom_rect(data=dfip, aes(xmin=AnyPos_T1, xmax=AnyT_T1, ymin=0, ymax=1400), fill='orange', alpha=0.1) +
				geom_rect(data=dfip, aes(xmin=AnyT_T1, xmax=DateDied, ymin=0, ymax=1400), fill='green', alpha=0.1) +
				geom_hline(yintercept=c(350,500), colour='blue', alpha=0.3) +
				geom_point(data=dfid, aes(x=PosCD4, y=CD4), colour='grey50') +
				geom_line(colour='black', aes(x=t, y=CD4)) + 				
				scale_y_continuous(breaks=seq(200, 2000, 200), expand=c(0,0)) +				
				scale_x_continuous(breaks=seq(1995, 2015, 5), minor_breaks=seq(1995, 2015, 1), expand=c(0,0)) +
				#coord_trans(limx=c(1996, 2013), limy=c(0,1500)) + 
				theme_bw() +
				facet_wrap(~Patient, scales='free') + labs(x='', y='CD4 count\n(cells/mm3)') +
				theme(strip.text=element_blank(), strip.background=element_blank())
		file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150303_CD4.pdf'
		ggsave(file=file, w=10, h=7)
		
		#file			<- paste(outdir, '/ATHENA0312_CD4endpoint_check.pdf',sep='')	
		#ggsave(file=file, w=10, h=3*tmp[, length(unique(Patient))], limitsize=FALSE )		
	}
	if(1)
	{
		t.period				<- 0.125
		t.recent.startctime		<- hivc.db.Date2numeric(as.Date("1996-07-15"))
		t.recent.startctime		<- floor(t.recent.startctime) + floor( (t.recent.startctime%%1)*100 %/% (t.period*100) ) * t.period
		t.recent.endctime		<- 2011
		#tp.cut				<- c(-Inf, 2006.5, 2008, 2009.5, 2011)
		tperiod.info<- as.data.table(structure(list(t.period = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), t.period.min = c(1996.5, 2006.5, 2008, 2009.5)+0.001, t.period.max = c(2006.5, 2008, 2009.5, 2011)), row.names = c(NA, -4L), class = "data.frame", .Names = c("t.period", "t.period.min", "t.period.max")))		
		
		df.cov	<- copy(df.all.allmsm)
		setkey(df.cov, Patient)
		df.cov	<- unique(df.cov)
		#
		# get first time CD4<350
		#
		file	<- paste(DATA,'/derived/','ATHENA_2014_06_Patient_AllMSM_CD4.R',sep='')
		load(file)		
		tmp		<- immu.sm[, {
					CD4_350_T1	<- which(CD4<=350)
					CD4_500_T1	<- which(CD4<=500)
					list(CD4sm='Y', CD4_350_T1=ifelse(length(CD4_350_T1), t, NA_real_), CD4_500_T1=ifelse(length(CD4_500_T1), t, NA_real_))
				}, by='Patient']
		df.cov	<- merge(df.cov, tmp, all.x=TRUE, by='Patient')
		df.cov	<- subset(df.cov, Trm=='MSM' | Trm=='BI')
		cat(paste('\nNumber of patients for whom we don t have a smooth, n=', subset(df.cov, is.na(CD4sm))[, length(unique(Patient))] ))
		#	for those for whom we have one CD4 count, set according to that
		tmp		<- df.cov[, which(is.na(CD4sm) & !is.na(CD4_T1) & CD4_T1<=500)]	
		set( df.cov, tmp, 'CD4_500_T1', df.cov[tmp, PosCD4_T1] )
		tmp		<- df.cov[, which(is.na(CD4sm) & !is.na(CD4_T1) & CD4_T1<=350)]	
		set( df.cov, tmp, 'CD4_350_T1', df.cov[tmp, PosCD4_T1] )
		#	for all remaining ones without a smooth, we don t have a single CD4 count
		cat(paste('\nNumber of patients for whom we don t have a single CD4 count, n=', nrow(unique(subset(df.cov, is.na(CD4_T1), Patient))) ))
		#	13 patients missing that have 1 CD4 count and not interpolated
		#	1044 patients without CD4		
		df.cd4cov	<- data.table(t= seq(t.recent.startctime, t.recent.endctime, 1/12))
		df.cd4cov	<- df.cd4cov[, list(	N_CD4_350=		nrow(subset(df.cov, t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & t>=CD4_350_T1 )),
						N_CD4_l500=		nrow(subset(df.cov, t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & t>=CD4_500_T1 & (is.na(CD4_350_T1) |  t<CD4_350_T1))),
						N_CD4_g500=		nrow(subset(df.cov, t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & t>=PosCD4_T1 & (t<CD4_500_T1 | is.na(CD4_500_T1)))),
						N_CD4_NA= 		nrow(subset(df.cov, t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & (t<PosCD4_T1 | is.na(PosCD4_T1)) )),
						N_CD4_350_ART=	nrow(subset(df.cov, t>=AnyT_T1 & t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & t>=CD4_350_T1)),	
						N_CD4_l500_ART=	nrow(subset(df.cov, t>=AnyT_T1 & t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & t>=CD4_500_T1 & (is.na(CD4_350_T1) |  t<CD4_350_T1))),
						N_CD4_g500_ART=	nrow(subset(df.cov, t>=AnyT_T1 & t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & t>=PosCD4_T1 &  (t<CD4_500_T1 | is.na(CD4_500_T1)))),
						N_CD4_NA_ART= 	nrow(subset(df.cov, t>=AnyT_T1 & t>=AnyPos_T1 & (is.na(DateDied) | DateDied<t) & (t<PosCD4_T1 | is.na(PosCD4_T1)) ))										
				), by='t']
		df.cd4cov[, N_CD4_350_NoART:= N_CD4_350-N_CD4_350_ART]
		df.cd4cov[, N_CD4_l500_NoART:= N_CD4_l500-N_CD4_l500_ART]
		df.cd4cov[, N_CD4_g500_NoART:= N_CD4_g500-N_CD4_g500_ART]
		df.cd4cov[, N_CD4_NA_NoART:= N_CD4_NA-N_CD4_NA_ART]
		df.cd4cov[, N_NoART:=N_CD4_350_NoART+N_CD4_l500_NoART+N_CD4_g500_NoART+N_CD4_NA_NoART]
		set(df.cd4cov, NULL, 't.period', df.cd4cov[, cut(t, breaks=c(-Inf, 2006.5, 2008, 2009.5, 2011), labels=1:4)])
		#
		#	some plotting
		#
		if(0)
		{
			tmp	<- melt(df.cd4cov, id.vars=c('t'), measure.vars=c('N_CD4_350_NoART','N_CD4_l500_NoART','N_CD4_g500_NoART','N_CD4_NA_NoART','N_NoART'))			
			tmp2	<- data.table(	variable= c('N_CD4_350_NoART','N_CD4_l500_NoART','N_CD4_g500_NoART','N_CD4_NA_NoART','N_NoART'), 
									colour=c("#41B6C4", "#74A9CF","#0570B0","#35978F", "grey70"))
			tmp		<- merge(tmp, tmp2, by='variable')			
			set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('N_CD4_g500_NoART','N_CD4_l500_NoART','N_CD4_350_NoART','N_CD4_NA_NoART','N_NoART'), labels=c('>500','350-500','<350', 'no CD4 measured','(All)')) ])
			setkey(tmp, variable)
			ggplot(tmp, aes(x=t, ymax=value, ymin=0, group=variable, fill=variable)) + 
					geom_ribbon(size=1)  + 
					facet_wrap(~variable, scales='free', nrow=3) +
					labs(x='', y='diagnosed but untreated MSM\n(#)') +
					scale_x_continuous(breaks=seq(1996,2020,2)) +
					scale_y_continuous(limit=c(0,NA)) +
					scale_fill_manual(name='CD4 progression to', values=tmp[, unique(colour)], guide=T) +
					theme_bw() + theme(legend.position='bottom',  strip.background = element_blank(), strip.text = element_blank()) +
					guides(fill=guide_legend(ncol=2))
			file	<- file				<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150123_ATHENA_2014_06_Patient_AllMSM_ARTno_by_CD4.pdf'
			ggsave(w=10,h=6,file=file)
			
			ggplot(df.cd4cov, aes(x=t, y=N_CD4_350_ART/N_CD4_350)) + geom_bar(stat='identity') + theme_bw() +
					scale_y_continuous(breaks=seq(0,1,0.1)) + scale_x_continuous(breaks=seq(1996,2014,2)) +
					labs(x='', y='ART coverage\namong diagnosed with CD4 < 350')
			file		<- paste(DATA, '/tmp/', 'ATHENA_2014_06_Patient_AllMSM_ARTcoverage_CD4350.pdf',sep='')
			ggsave(file=file, w=6, h=4)
			
			tmp			<- melt(df.cd4cov, measure.vars=c('N_CD4_350_ART','N_CD4_350_NoART'), id.vars='t', variable.name='CD4_350', value.name='N')
			set(tmp, NULL, 'CD4_350', tmp[, sapply(strsplit(as.character(CD4_350),'_',fixed=1),'[[',4)])
			set(tmp, NULL, 'CD4_350', tmp[, factor(CD4_350, labels=c('NoART','ART'), levels=c('NoART','ART'))])
			ggplot(tmp, aes(x=t, y=N, fill=CD4_350)) + geom_bar(stat='identity') + scale_fill_brewer(palette='Set1') + theme_bw() +
					scale_y_continuous(breaks=seq(0,1e4,500)) + scale_x_continuous(breaks=seq(1996,2014,2)) +
					labs(x='', y='Patients with CD4<350')
			file		<- paste(DATA, '/tmp/', 'ATHENA_2014_06_Patient_AllMSM_ARTcoverage_CD4350n.pdf',sep='')
			ggsave(file=file, w=6, h=4)			
		}
		#	load estimated prop and Nt
		file	<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/fisheretal_150216/ATHENA_2013_03_-DR-RC-SH+LANL_Sequences_Wed_Dec_18_11:37:00_2013_method.risks.Rdata'
		load(file)
		#
		#	< 350
		#
		if(1)
		{
			df.Praw.e0cp	<- subset(runs.risk, is.na(t.period) & method.nodectime=='any' & method.brl=='3pa1H1.35C3V100bInfT7' & method.denom=='SEQ' & method.recentctime==2011 & method.dating=='sasky' & stat=='P.raw.e0cp' & grepl('stageDtl350', coef))
			set(df.Praw.e0cp, NULL, 't.period', df.Praw.e0cp[, as.character(factor)])
			set(df.Praw.e0cp, NULL, 't.period', df.Praw.e0cp[, factor(substr(t.period, nchar(t.period), nchar(t.period)))])
			#df.cd4cov	<- subset(df.cd4cov, select=c(t, N_CD4_350_NoART))
			df.cd4me	<- tperiod.info[, list(M_CD4_350_NoART=subset(df.cd4cov, t>=t.period.min & t<t.period.max)[, mean(N_CD4_350_NoART)]), by=t.period]
			df.cd4me	<- merge(tperiod.info, df.cd4me, by='t.period')
			df.cd4me	<- merge(df.cd4me, df.Praw.e0cp, by='t.period')
			tmp			<- lm(M_CD4_350_NoART~v-1, data=df.cd4me)
			ans			<- lm(M_CD4_350_NoART~v, data=df.cd4me)
			set(df.cd4me, NULL, 'tv', predict(tmp))
			set(df.cd4me, NULL, 'tl95.bs', predict(tmp, data.frame(v= df.cd4me$l95.bs)))
			set(df.cd4me, NULL, 'tu95.bs', predict(tmp, data.frame(v= df.cd4me$u95.bs)))
			set(df.cd4cov, NULL, 't.period', df.cd4cov[, cut(t, breaks=c(-Inf, 2006.5, 2008, 2009.5, 2011), labels=1:4)])
			
			summary(ans)
			#lm(formula = M_CD4_350_NoART ~ v, data = df.cd4me)
			#
			#Coefficients:
			#           Estimate Std. Error t value Pr(>|t|)  
			#(Intercept)    64.55     105.74   0.610    0.604  
			#	16389.67    3052.63   5.369    0.033 *
			summary(tmp)
			#lm(formula = M_CD4_350_NoART ~ v - 1, data = df.cd4me)
			#Coefficients:
			#		Estimate Std. Error t value Pr(>|t|)    
			#v  17040.5      640.1   26.62 0.000116 ***
			
			#
			#	final plot 350
			#
			library(gtable)
			library(grid)
			p1	<- ggplot(df.cd4cov) + geom_step(aes(x=t, y=N_CD4_350_NoART), stat='identity', colour="#41B6C4", width=1/12, alpha=0.4) + theme_bw() +
					scale_x_continuous(breaks=seq(1997,2014,2), minor_breaks=NULL, expand=c(0,0)) + 
					scale_y_continuous(breaks=seq(0,3000,250)) +
					geom_segment(data=df.cd4me, aes(x=t.period.min, xend=t.period.max, y=M_CD4_350_NoART, yend=M_CD4_350_NoART), col="#41B6C4", size=1.4) +
					geom_point(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, y=tv), col='black', shape=18, size=4) +				
					geom_errorbar(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, ymin=tl95.bs, ymax=tu95.bs), col='black', width=0.5, size=1) +
					labs(x='', y='Untreated HIV infected MSM\nwith CD4<350\n( # )') + 
					facet_grid(.~t.period, scales='free_x', space='free_x') +					
					theme(axis.title=element_text(size=18), axis.text.x=element_text(size=14), axis.text.y=element_text(colour="#41B6C4", size=14), strip.background = element_blank(), strip.text = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x=element_blank(), plot.margin=unit(c(0,2,0,0),"cm"))
			p2	<- ggplot(df.cd4cov) + geom_step(aes(x=t, y=N_CD4_350_NoART), stat='identity', colour="#41B6C4", width=1/12, alpha=0.4) + theme_bw() +
					scale_x_continuous(breaks=seq(1997,2014,2), minor_breaks=NULL, expand=c(0,0)) + 
					scale_y_continuous(breaks=predict(tmp, data.frame(v= seq(0, 14, 2)/100)), label=seq(0, 14, 2)) +					
					geom_segment(data=df.cd4me, aes(x=t.period.min, xend=t.period.max, y=M_CD4_350_NoART, yend=M_CD4_350_NoART), col="#41B6C4", size=1.4) +
					geom_point(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, y=tv), col='black', shape=18, size=4) +				
					geom_errorbar(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, ymin=tl95.bs, ymax=tu95.bs), col='black', width=0.5, size=1) +
					labs(x='', y='Transmissions from\nuntreated men with CD4<350\n( % )') + 
					facet_grid(.~t.period, scales='free_x', space='free_x') +										
					theme(axis.title=element_text(size=18), axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), strip.background = element_blank(), strip.text = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x=element_blank(), plot.margin=unit(c(0,2,0,0),"cm"))
			#	prepare y axis label 	
			p2	<- p2 + theme(axis.text.y = element_text(hjust = 0), axis.title.y = element_text(angle = 270)) 
			#extract gtable
			g1	<- ggplot_gtable(ggplot_build(p1))
			g2	<- ggplot_gtable(ggplot_build(p2))		
			#overlap the panel of the 2nd plot on that of the 1st plot		
			pp	<- c(subset(g1$layout, name=="panel", se=t:r)[1, ])
			g	<- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")[1]]], pp$t, pp$l, pp$b, pp$l)
			#extract left axis
			ia <- which(g2$layout$name == "axis-l")
			ga <- g2$grobs[[ia]]
			ax <- ga$children[[2]]
			ax$widths <- rev(ax$widths)
			ax$grobs <- rev(ax$grobs)
			ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
			g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
			g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
			#extract left axis label		
			g 	<- gtable_add_grob(g, g2$grobs[[which(g2$layout$name == "ylab")]], pp$t, length(g$widths), pp$b)
			#draw the whole thing		
			file<- paste(DATA, '/fisheretal_150216/', 'ATHENA_2014_06_Patient_AllMSM_ARTno_CD4350.pdf',sep='')
			pdf(file=file, w=6, h=4) # plot saved by default to Rplots.pdf
			grid.newpage()
			grid.draw(g)
			dev.off() 	
		}
		#
		#	350-500
		#
		if(1)
		{
			df.Praw.e0cp	<- subset(runs.risk, is.na(t.period) & method.nodectime=='any' & method.brl=='3pa1H1.35C3V100bInfT7' & method.denom=='SEQ' & method.recentctime==2011 & method.dating=='sasky' & stat=='P.raw.e0cp' & grepl('stageDtl500', coef))
			set(df.Praw.e0cp, NULL, 't.period', df.Praw.e0cp[, as.character(factor)])
			set(df.Praw.e0cp, NULL, 't.period', df.Praw.e0cp[, factor(substr(t.period, nchar(t.period), nchar(t.period)))])
			df.cd4me	<- tperiod.info[, list(M_CD4_l500_NoART=subset(df.cd4cov, t>=t.period.min & t<t.period.max)[, mean(N_CD4_l500_NoART)]), by=t.period]
			df.cd4me	<- merge(tperiod.info, df.cd4me, by='t.period')
			df.cd4me	<- merge(df.cd4me, df.Praw.e0cp, by='t.period')
			tmp			<- lm(M_CD4_l500_NoART~v-1, data=df.cd4me)
			ans			<- lm(M_CD4_l500_NoART~v, data=df.cd4me)
			set(df.cd4me, NULL, 'tv', predict(tmp))
			set(df.cd4me, NULL, 'tl95.bs', predict(tmp, data.frame(v= df.cd4me$l95.bs)))
			set(df.cd4me, NULL, 'tu95.bs', predict(tmp, data.frame(v= df.cd4me$u95.bs)))
			
			
			summary(ans)
			#lm(formula = M_CD4_350_NoART ~ v, data = df.cd4me)
			#
			#Coefficients:
			#             Estimate Std. Error t value Pr(>|t|)
			#(Intercept)   -16.28     369.52  -0.044    0.969
			#v            7871.68    5866.39   1.342    0.312
			summary(tmp)
			#	Estimate Std. Error t value Pr(>|t|)  
			#	v     7624       1346   5.663   0.0109 *
			#
			#	final plot 350-500
			#
			library(gtable)
			library(grid)
			p1	<- ggplot(df.cd4cov) + geom_step(aes(x=t, y=N_CD4_l500_NoART), stat='identity', colour="#74A9CF", width=1/12, alpha=0.4) + theme_bw() +
					scale_x_continuous(breaks=seq(1997,2014,2), minor_breaks=NULL, expand=c(0,0)) + 
					scale_y_continuous(breaks=seq(0,3000,250)) +
					geom_segment(data=df.cd4me, aes(x=t.period.min, xend=t.period.max, y=M_CD4_l500_NoART, yend=M_CD4_l500_NoART), col="#74A9CF", size=1.4) +
					geom_point(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, y=tv), col='black', shape=18, size=4) +				
					geom_errorbar(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, ymin=tl95.bs, ymax=tu95.bs), col='black', width=0.5, size=1) +
					labs(x='', y='Untreated HIV infected MSM\nwith CD4 350-500\n( # )') + 
					facet_grid(.~t.period, scales='free_x', space='free_x') +					
					theme(axis.title=element_text(size=18), axis.text.x=element_text(size=14), axis.text.y=element_text(colour="#74A9CF", size=14), strip.background = element_blank(), strip.text = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x=element_blank(), plot.margin=unit(c(0,2,0,0),"cm"))
			p2	<- ggplot(df.cd4cov) + geom_step(aes(x=t, y=N_CD4_l500_NoART), stat='identity', colour="#74A9CF", width=1/12, alpha=0.4) + theme_bw() +
					scale_x_continuous(breaks=seq(1997,2014,2), minor_breaks=NULL, expand=c(0,0)) + 
					scale_y_continuous(breaks=predict(tmp, data.frame(v= seq(0, 14, 2)/100)), label=seq(0, 14, 2)) +					
					geom_segment(data=df.cd4me, aes(x=t.period.min, xend=t.period.max, y=M_CD4_l500_NoART, yend=M_CD4_l500_NoART), col="#74A9CF", size=1.4) +
					geom_point(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, y=tv), col='black', shape=18, size=4) +				
					geom_errorbar(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, ymin=tl95.bs, ymax=tu95.bs), col='black', width=0.5, size=1) +
					labs(x='', y='Transmissions from\nuntreated men with CD4 350-500\n( % )') + 
					facet_grid(.~t.period, scales='free_x', space='free_x') +										
					theme(axis.title=element_text(size=18), axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), strip.background = element_blank(), strip.text = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x=element_blank(), plot.margin=unit(c(0,2,0,0),"cm"))
			#	prepare y axis label 	
			p2	<- p2 + theme(axis.text.y = element_text(hjust = 0), axis.title.y = element_text(angle = 270)) 
			#extract gtable
			g1	<- ggplot_gtable(ggplot_build(p1))
			g2	<- ggplot_gtable(ggplot_build(p2))		
			#overlap the panel of the 2nd plot on that of the 1st plot		
			pp	<- c(subset(g1$layout, name=="panel", se=t:r)[1, ])
			g	<- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")[1]]], pp$t, pp$l, pp$b, pp$l)
			#extract left axis
			ia <- which(g2$layout$name == "axis-l")
			ga <- g2$grobs[[ia]]
			ax <- ga$children[[2]]
			ax$widths <- rev(ax$widths)
			ax$grobs <- rev(ax$grobs)
			ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
			g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
			g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
			#extract left axis label		
			g 	<- gtable_add_grob(g, g2$grobs[[which(g2$layout$name == "ylab")]], pp$t, length(g$widths), pp$b)
			#draw the whole thing		
			file<- paste(DATA, '/fisheretal_150105/', 'ATHENA_2014_06_Patient_AllMSM_ARTno_CD4350500.pdf',sep='')
			pdf(file=file, w=6, h=4) # plot saved by default to Rplots.pdf
			grid.newpage()
			grid.draw(g)
			dev.off() 	
		}
		#
		#	> 500
		#
		if(1)
		{
			df.Praw.e0cp	<- subset(runs.risk, is.na(t.period) & method.nodectime=='any' & method.brl=='3pa1H1.35C3V100bInfT7' & method.denom=='SEQ' & method.recentctime==2011 & method.dating=='sasky' & stat=='P.raw.e0cp' & grepl('stageDtg500', coef))
			set(df.Praw.e0cp, NULL, 't.period', df.Praw.e0cp[, as.character(factor)])
			set(df.Praw.e0cp, NULL, 't.period', df.Praw.e0cp[, factor(substr(t.period, nchar(t.period), nchar(t.period)))])
			df.cd4me	<- tperiod.info[, list(M_CD4_g500_NoART=subset(df.cd4cov, t>=t.period.min & t<t.period.max)[, mean(N_CD4_g500_NoART)]), by=t.period]
			df.cd4me	<- merge(tperiod.info, df.cd4me, by='t.period')
			df.cd4me	<- merge(df.cd4me, df.Praw.e0cp, by='t.period')
			tmp			<- lm(M_CD4_g500_NoART~v-1, data=df.cd4me)
			ans			<- lm(M_CD4_g500_NoART~v, data=df.cd4me)
			set(df.cd4me, NULL, 'tv', predict(tmp))
			set(df.cd4me, NULL, 'tl95.bs', predict(tmp, data.frame(v= df.cd4me$l95.bs)))
			set(df.cd4me, NULL, 'tu95.bs', predict(tmp, data.frame(v= df.cd4me$u95.bs)))
			
			summary(ans)
			#lm(formula = M_CD4_350_NoART ~ v, data = df.cd4me)
			#
			#Coefficients:
			#             Estimate Std. Error t value Pr(>|t|)
			#(Intercept)    658.3      907.5   0.725    0.544
			#v            -3784.7    10203.6  -0.371    0.746
			summary(tmp)		
			#Estimate Std. Error t value Pr(>|t|)  
			#v     3559       1171    3.04   0.0559 .
			
			#
			#	final plot > 500
			#
			library(gtable)
			library(grid)
			p1	<- ggplot(df.cd4cov) + geom_step(aes(x=t, y=N_CD4_g500_NoART), stat='identity', colour="#0570B0", width=1/12, alpha=0.4) + theme_bw() +
					scale_x_continuous(breaks=seq(1997,2014,2), minor_breaks=NULL, expand=c(0,0)) + 
					scale_y_continuous(breaks=seq(0,3000,250)) +
					geom_segment(data=df.cd4me, aes(x=t.period.min, xend=t.period.max, y=M_CD4_g500_NoART, yend=M_CD4_g500_NoART), col="#0570B0", size=1.4) +
					geom_point(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, y=tv), col='black', shape=18, size=4) +				
					geom_errorbar(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, ymin=tl95.bs, ymax=tu95.bs), col='black', width=0.5, size=1) +
					labs(x='', y='Untreated HIV infected MSM\nwith CD4>500\n( # )') + 
					facet_grid(.~t.period, scales='free_x', space='free_x') +					
					theme(axis.title=element_text(size=18), axis.text.x=element_text(size=14), axis.text.y=element_text(colour="#0570B0", size=14), strip.background = element_blank(), strip.text = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x=element_blank(), plot.margin=unit(c(0,2,0,0),"cm"))
			p2	<- ggplot(df.cd4cov) + geom_step(aes(x=t, y=N_CD4_g500_NoART), stat='identity', colour="#0570B0", width=1/12, alpha=0.4) + theme_bw() +
					scale_x_continuous(breaks=seq(1997,2014,2), minor_breaks=NULL, expand=c(0,0)) + 
					scale_y_continuous(breaks=predict(tmp, data.frame(v= seq(0, 14, 2)/100)), label=seq(0, 14, 2)) +					
					geom_segment(data=df.cd4me, aes(x=t.period.min, xend=t.period.max, y=M_CD4_g500_NoART, yend=M_CD4_g500_NoART), col="#0570B0", size=1.4) +
					geom_point(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, y=tv), col='black', shape=18, size=4) +				
					geom_errorbar(data=df.cd4me, aes(x=(t.period.min+t.period.max)/2, ymin=tl95.bs, ymax=tu95.bs), col='black', width=0.5, size=1) +
					labs(x='', y='Transmissions from\nuntreated men with CD4>500\n( % )') + 
					facet_grid(.~t.period, scales='free_x', space='free_x') +										
					theme(axis.title=element_text(size=18), axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), strip.background = element_blank(), strip.text = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x=element_blank(), plot.margin=unit(c(0,2,0,0),"cm"))
			#	prepare y axis label 	
			p2	<- p2 + theme(axis.text.y = element_text(hjust = 0), axis.title.y = element_text(angle = 270)) 
			#extract gtable
			g1	<- ggplot_gtable(ggplot_build(p1))
			g2	<- ggplot_gtable(ggplot_build(p2))		
			#overlap the panel of the 2nd plot on that of the 1st plot		
			pp	<- c(subset(g1$layout, name=="panel", se=t:r)[1, ])
			g	<- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")[1]]], pp$t, pp$l, pp$b, pp$l)
			#extract left axis
			ia <- which(g2$layout$name == "axis-l")
			ga <- g2$grobs[[ia]]
			ax <- ga$children[[2]]
			ax$widths <- rev(ax$widths)
			ax$grobs <- rev(ax$grobs)
			ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
			g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
			g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
			#extract left axis label		
			g 	<- gtable_add_grob(g, g2$grobs[[which(g2$layout$name == "ylab")]], pp$t, length(g$widths), pp$b)
			#draw the whole thing		
			file<- paste(DATA, '/fisheretal_150105/', 'ATHENA_2014_06_Patient_AllMSM_ARTno_CD4g500.pdf',sep='')
			pdf(file=file, w=6, h=4) # plot saved by default to Rplots.pdf
			grid.newpage()
			grid.draw(g)
			dev.off() 	
		}
	}
}
######################################################################################
project.athena.Fisheretal.composition.putativeinfectionwindow<- function()
{
	tmp		<- subset(ri.ALLMSM, AnyPos_T1>=1996.5 & AnyPos_T1<2011 & isAcute=='Yes' & Trm%in%c('MSM','BI'), select=c(Patient, isAcute, Acute_Spec, AnyPos_T1, NegT, Trm))
	
	tmp[, IPWd:= 1]
	tmp2	<- tmp[,which(!is.na(NegT) & AnyPos_T1-NegT<11/12)]
	set(tmp, tmp2, 'IPWd', tmp[tmp2, AnyPos_T1-NegT+1/12])
	tmp[, t.period:= tmp[, cut(AnyPos_T1, breaks=c(-Inf, 2006.5, 2008, 2009.5, 2011), labels=c('96/07-06/06','06/07-07/12','08/01-09/06','09/07-10/12'))]]
	
	ggplot(tmp, aes(x=IPWd*12)) + geom_histogram(binwidth=0.5) + facet_grid(.~t.period, scales='free', space='free') +
			scale_x_continuous(breaks=seq(4.25,12.25,2), labels=seq(4,12,2), lim=c(4,12.5)) +
			scale_y_log10(breaks=c(10, seq(20,100,20), 200))+
			labs(x='duration of putative infection window\n(months)', y='') +
			theme_bw()
	file	<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2014/MSMtransmission_ATHENA1303/150223_PutInfWindow_Duration.pdf'
	ggsave(file=file, w=10, h=4)
	#
	nrow(subset(tmp, IPWd<1)) / nrow(tmp)
	
}
######################################################################################
project.athena.Fisheretal.composition.t2inf.cd4.350etc<- function(df.sc.negT, adjust.dt.CD4, plot.file)
{
	df.sc.negT[, CD4_T1c:=  cut(CD4_T1, breaks=c(-Inf,250,350,500,850,Inf))]
	df.scay.cd4		<- subset( df.sc.negT, isAcute=='Yes' & dt.CD4<=adjust.dt.CD4 & PosCD4_T1<=AnyT_T1)	
	df.scaym.cd4	<- subset( df.sc.negT, is.na(isAcute) & Acute_Spec=='SYM' & dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)
	df.scana.cd4	<- subset( df.sc.negT, is.na(isAcute) & is.na(Acute_Spec) & dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)
	df.scchr.cd4	<- subset( df.sc.negT, isAcute=='No' & dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)
	#	loess fit
	tmp<- subset(df.sc.negT, dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)
	set(tmp, tmp[, which(is.na(isAcute))], 'isAcute','Missing infection status')
	set(tmp, tmp[, which(isAcute=='No')],'isAcute', 'Recent HIV infection\nnot indicated')
	set(tmp, tmp[, which(isAcute=='Yes')],'isAcute', 'Confirmed\nrecent HIV infection')
	ggplot(tmp, aes(x=AnyPos_a, y=mpy.NegT)) + 
			geom_point(size=0.85) + scale_y_continuous(limits=c(0,13), breaks=seq(0,20,2)) + xlim(18, 61) +
			geom_smooth(fill='blue', alpha=0.2) +
			labs(x="age at diagnosis", y='time between midpoint of\nseroconversion interval and diagnosis\n(years)') +
			facet_grid(isAcute ~ CD4_T1c) + theme_bw()
	ggsave(paste(plot.file,'sc_exploreloess.pdf'), w=8, h=8)
	#	Gamma fit
	#	confirmed recent
	tmp	<- subset(df.scay.cd4, CD4_T1c=='(-Inf,250]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- copy(tmp)
	tmp	<- subset(df.scay.cd4, CD4_T1c=='(250,350]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	tmp	<- subset(df.scay.cd4, CD4_T1c=='(350,500]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)		
	tmp	<- subset(df.scay.cd4, CD4_T1c=='(500,850]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	tmp	<- subset(df.scay.cd4, CD4_T1c=='(850, Inf]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=3)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	# 	unknown status / symptomatic
	tmp	<- subset(rbind(df.scana.cd4,df.scaym.cd4), CD4_T1c=='(-Inf,250]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	tmp	<- subset(df.scana.cd4, CD4_T1c=='(250,350]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	tmp	<- subset(df.scana.cd4, CD4_T1c=='(350,500]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)		
	tmp	<- subset(df.scana.cd4, CD4_T1c=='(500,850]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	tmp	<- subset(rbind(df.scana.cd4,df.scaym.cd4), CD4_T1c=='(850, Inf]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=3)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	# 	chronic
	tmp	<- subset(df.scchr.cd4, CD4_T1c=='(-Inf,250]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	tmp	<- subset(df.scchr.cd4, CD4_T1c=='(250,350]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	tmp	<- subset(df.scchr.cd4, CD4_T1c=='(350,500]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)		
	tmp	<- subset(df.scchr.cd4, CD4_T1c=='(500,850]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	tmp	<- subset(df.scchr.cd4, CD4_T1c=='(850, Inf]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=3)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)		
	ans		<- rbind(ans,tmp)
	#
	set(ans, ans[, which(isAcute=='No')],'isAcute', 'Chronic HIV infection')
	set(ans, ans[, which(isAcute=='Yes')],'isAcute', 'Confirmed recent HIV infection')
	set(ans, ans[, which(is.na(isAcute))],'isAcute', 'Unconfirmed infection status')
	#		
	ggplot(ans, aes(x=AnyPos_a, y=mpy.NegT)) + 
			geom_point(size=0.85) + scale_y_continuous(limits=c(0,13), breaks=seq(0,20,2)) + xlim(18, 61) +			
			geom_line(aes(y=y.b), colour='blue') +
			geom_ribbon(aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2, fill='blue', alpha=0.5) +
			labs(x="age at diagnosis", y='time between midpoint of\nseroconversion interval and diagnosis\n(years)') +
			facet_grid(isAcute ~ CD4_T1c) + theme_bw()
	ggsave(paste(plot.file,'sc_explore.pdf'), w=8, h=8)
}
######################################################################################
project.athena.Fisheretal.composition.t2inf.cd4.850etc<- function(df.sc.negT, adjust.dt.CD4, plot.file)
{
	df.sc.negT[, CD4_T1c:=  cut(CD4_T1, breaks=c(-Inf,250,850,Inf))]
	df.scay.cd4		<- subset( df.sc.negT, isAcute=='Yes' & dt.CD4<=adjust.dt.CD4 & PosCD4_T1<=AnyT_T1)	
	df.scaym.cd4	<- subset( df.sc.negT, is.na(isAcute) & Acute_Spec=='SYM' & dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)
	df.scana.cd4	<- subset( df.sc.negT, is.na(isAcute) & is.na(Acute_Spec) & dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)
	df.scchr.cd4	<- subset( df.sc.negT, isAcute=='No' & dt.CD4<adjust.dt.CD4 & PosCD4_T1<AnyT_T1)		
	#	confirmed recent
	tmp	<- subset(df.scay.cd4, CD4_T1c=='(-Inf,250]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- copy(tmp)
	tmp	<- subset(df.scay.cd4, CD4_T1c=='(250,850]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	tmp	<- subset(df.scay.cd4, CD4_T1c=='(850, Inf]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=3)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	# 	unknown status / symptomatic
	tmp	<- subset(rbind(df.scana.cd4,df.scaym.cd4), CD4_T1c=='(-Inf,250]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	tmp	<- subset(rbind(df.scana.cd4,df.scaym.cd4), CD4_T1c=='(250,850]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	tmp	<- subset(rbind(df.scana.cd4,df.scaym.cd4), CD4_T1c=='(850, Inf]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=3)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	# 	chronic
	tmp	<- subset(df.scchr.cd4, CD4_T1c=='(-Inf,250]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]	
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)
	ans		<- rbind(ans,tmp)
	tmp	<- subset(df.scchr.cd4, CD4_T1c=='(250,850]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=5)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)		
	ans		<- rbind(ans,tmp)
	tmp	<- subset(df.scchr.cd4, CD4_T1c=='(850, Inf]', select=c(Patient, mpy.NegT, AnyPos_a, CD4_T1c, isAcute, Acute_Spec))		
	setkey(tmp, AnyPos_a)
	tmp2	<- gamlss(as.formula('mpy.NegT ~ bs(AnyPos_a, degree=3)'), sigma.formula=as.formula('~ bs(AnyPos_a, degree=3)'), data=as.data.frame(subset(tmp, select=c(mpy.NegT, AnyPos_a,CD4_T1c))), family=GA)				 
	tmp[, y.b:= predict(tmp2, type='response', se.fit=FALSE)]
	tmp3	<- gamlss.centiles.get(tmp2, tmp$AnyPos_a, cent = c(2.5, 97.5), with.ordering=FALSE )
	tmp		<- cbind(tmp, tmp3)		
	ans		<- rbind(ans,tmp)
	#
	set(ans, ans[, which(isAcute=='No')],'isAcute', 'Chronic HIV infection')
	set(ans, ans[, which(isAcute=='Yes')],'isAcute', 'Confirmed recent HIV infection')
	set(ans, ans[, which(is.na(isAcute))],'isAcute', 'Unconfirmed infection status')
	
	
	ggplot(ans, aes(x=AnyPos_a, y=mpy.NegT)) + 
			geom_point(size=0.85) + scale_y_continuous(limits=c(0,13), breaks=seq(0,20,2)) + xlim(18, 61) +			
			geom_line(aes(y=y.b), colour='blue') +
			geom_ribbon(aes(x=x, ymin=q2.5, ymax=q97.5), alpha=0.2, fill='blue', alpha=0.5) +
			labs(x="age at diagnosis", y='time between midpoint of\nseroconversion interval and diagnosis\n(years)') +
			facet_grid(isAcute ~ CD4_T1c) + theme_bw()
	ggsave(paste(plot.file,'sc_explore.pdf'), w=8, h=8)
}
