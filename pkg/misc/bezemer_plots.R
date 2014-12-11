db.clusterR.plot<- function()
{
	require(data.table)
	indir	<- '/Users/Oliver/Dropbox (Infectious Disease)/material for paper/data_SA_Rt'
	infile	<- 'R-Estimates_SA80p10.csv'
	indir	<- '/Users/Oliver/Dropbox (Infectious Disease)/material for paper'
	infile	<- 'R-Estimates_MSM_only.csv'
	file	<- paste(indir, '/', infile, sep='')
	read.csv(file)	
	df		<- as.data.table(read.csv(file))
	setnames(df, c('ClusterName','ClusterSizeWithDiagDate','TimeLastREstimate'), c('CID','CN','T'))
	
	indir	<- '/Users/Oliver/Dropbox (Infectious Disease)/material for paper/data_SA_Rt'
	infile	<- 'WT-smoothed_SA80p10.RData'	
	file	<- paste(indir, '/', infile, sep='')
	z		<- load(file)
	#Transmission networks are shown in blue if their first sequence appeared before 1991, 
	# in green if their first sequence appeared between 1991 and 2000, and 
	# in yellow if it appeared after 2000. The red horizontal dotted line represents the threshold value
	
	df		<- lapply(seq_along(results_width_5), function(i)
			{
				tmp<- as.data.table( results_width_5[[i]][['R']] )
				tmp[, CID:= names(results_width_5)[i]]
				tmp
			})	
	df		<- do.call('rbind', df)
	setnames(df, c('T.Start','T.End','Mean(R)','Std(R)','Quantile.0.025(R)','Quantile.0.975(R)'),
					c('TS','TE','RM','RS','RQL','RQU'))
	df		<- subset(df, !is.nan(RM))
	set(df, NULL, 'TS', df[,TS+1980])
	set(df, NULL, 'TE', df[,TE+1980])
	df		<- merge( df, df[, list(TMID= mean(c(TS, TE))), by=c('CID', 'TS') ], by=c('CID', 'TS') )
	
	tmp		<- df[, list(TMIN=min(TS)), by='CID']
	tmp[, CID.ERA:= tmp[, cut(TMIN, breaks=c(-Inf, 1990, 2000, Inf), labels=c("before 1991","1991 - 2000","after 2000"))]]	
	df		<- merge(df, tmp, by='CID')
	
	dfq		<- df[, list(RQLM= mean(RQL), RQUM= mean(RQU), RMQL=quantile(RM, p=0.025), RMM=quantile(RM, p=0.5), RMQU=quantile(RM, p=0.975)), by='TMID']
	dfq		<- subset(dfq, RMQL<RMQU)
	setkey(dfq, TMID)
	
	dfq2		<- df[, list(RQLM= mean(RQL), RQUM= mean(RQU), RMQL=quantile(RM, p=0.025), RMM=quantile(RM, p=0.5), RMQU=quantile(RM, p=0.975)), by=c('CID.ERA','TMID')]
	dfq2		<- subset(dfq2, RMQL<RMQU)
	setkey(dfq2, TMID)
	#
	#	new code: taking Anne's comment on board
	#	
	ggplot(df, aes(x=TMID)) + theme_bw() +	labs(x='', y='Case reproduction number within cluster\n(clusters of size > 10)\n') +		
			geom_ribbon(alpha=0.01,  fill='black', aes(y=RM, ymax=RQU, ymin=RQL, group=CID)) +
			geom_line(alpha=0.5, aes(y=RM, group=CID, colour=CID.ERA)) +
			geom_line(data=dfq, linetype='dashed', aes(y=RQLM)) +
			geom_line(data=dfq, linetype='dashed', aes(y=RQUM)) +
			geom_line(alpha=0.5, aes(y=RM, group=CID, colour=CID.ERA)) +
			geom_line(data=dfq, aes(y=RMM),size=2, colour='black') +
			geom_hline(yintercept=1, linetype='dotdash') +			
			scale_colour_brewer(palette='Set1', name='Era') +
			coord_cartesian(ylim=c(0,2.5), xlim=c(1985, 2008)) +
			scale_y_continuous(breaks=seq(0,3,0.2)) + scale_x_continuous(breaks=seq(1980,2020,5), minor_breaks=seq(1980,2020,1)) +
			theme(legend.position=c(0,1), legend.justification=c(0,1), legend.direction='horizontal', legend.key.size=unit(10,'mm'))
	file	<- paste(indir, '/', substr(infile,1,nchar(infile)-6),'_RAllEras.pdf', sep='')
	ggsave(file=file, w=8, h=6)	
	
	ggplot(df, aes(x=TMID)) + theme_bw() +	labs(x='', y='Case reproduction number within cluster\n(clusters of size > 10)\n') +			
			geom_ribbon(alpha=0.01,  fill='black', aes(y=RM, ymax=RQU, ymin=RQL, group=CID)) +
			geom_line(alpha=0.5, aes(y=RM, group=CID, colour=CID.ERA)) +
			geom_line(data=dfq2, linetype='dashed', aes(y=RQLM)) +
			geom_line(data=dfq2, linetype='dashed', aes(y=RQUM)) +
			geom_line(alpha=0.5, aes(y=RM, group=CID, colour=CID.ERA)) +
			geom_line(data=dfq2, aes(y=RMM),size=2, colour='black') +
			geom_hline(yintercept=1, linetype='dotdash') +			
			scale_colour_brewer(palette='Set1', name='Era', guide=FALSE) +
			coord_cartesian(ylim=c(0,2.5), xlim=c(1985, 2008)) +
			scale_y_continuous(breaks=seq(0,3,0.2)) + scale_x_continuous(breaks=seq(1980,2020,5), minor_breaks=seq(1980,2020,1)) +
			theme(legend.position=c(0,1), legend.justification=c(0,1), legend.direction='horizontal', legend.key.size=unit(10,'mm')) +
			facet_grid(.~CID.ERA)
	file	<- paste(indir, '/', substr(infile,1,nchar(infile)-6),'_RByEra.pdf', sep='')
	ggsave(file=file, w=12, h=5)
	
	#
	#	Rmean for all eras
	#
	ggplot(df, aes(x=TMID)) + theme_bw() + labs(x='', y='mean R per cluster') +
			geom_line(data=dfq, aes(y=RMM),size=2, colour='black') +
			geom_ribbon(data=dfq, aes(ymin=RMQL, ymax=RMQU), fill='black', alpha=0.2) +
			geom_line(alpha=0.5, aes(y=RM, group=CID, colour=CID.ERA)) +
			geom_hline(yintercept=1, linetype='dashed') +
			#geom_errorbar(aes(ymin=RQL, ymax=RQU), colour="black", alpha=0.2, width=.1, position=position_dodge(.5)) +
			scale_colour_brewer(palette='Set1', name='era') +
			scale_y_continuous(breaks=seq(0,3,0.2)) + scale_x_continuous(breaks=seq(1980,2020,5), minor_breaks=seq(1980,2020,1)) +
			theme(legend.position='bottom')
	file	<- paste(indir, '/', substr(infile,1,nchar(infile)-6),'_RmeanAllEras.pdf', sep='')
	ggsave(file=file, w=8, h=6)	
	#
	#	Rmean by eras
	#	
	ggplot(df, aes(x=TMID)) + theme_bw() + labs(x='', y='mean R per cluster') +
			geom_line(data=dfq2, aes(y=RMM),size=2, colour='black') +
			geom_ribbon(data=dfq2, aes(ymin=RMQL, ymax=RMQU), fill='black', alpha=0.2) +
			geom_line(alpha=0.5, aes(y=RM, group=CID, colour=CID.ERA)) +
			geom_hline(yintercept=1, linetype='dashed') +
			#geom_errorbar(aes(ymin=RQL, ymax=RQU), colour="black", alpha=0.2, width=.1, position=position_dodge(.5)) +
			scale_colour_brewer(palette='Set1', name='era') +
			scale_y_continuous(breaks=seq(0,3,0.2)) + scale_x_continuous(breaks=seq(1980,2020,5), minor_breaks=seq(1980,2020,1)) +
			theme(legend.position='bottom') +
			facet_grid(.~CID.ERA)
	file	<- paste(indir, '/', substr(infile,1,nchar(infile)-6),'_RmeanByEra.pdf', sep='')
	ggsave(file=file, w=12, h=6)
	
	
}