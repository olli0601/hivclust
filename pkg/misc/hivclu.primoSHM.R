primoSHM.plotViralLoad<- function()
{
	dfv		<- subset( df.viro.allmsm, select=c(Patient, PosRNA, lRNA) )
	set(dfv, NULL, 'PosRNA', hivc.db.Date2numeric(dfv[,PosRNA]))
	dfpr	<- subset(df.all.allmsm, WithPrimoSHMID=='Yes')
	set(dfpr, dfpr[, which(is.na(DateDied))], 'DateDied', 2013.3)
	set(dfpr, dfpr[, which(is.na(AnyT_T1))], 'AnyT_T1', 2013.3)
	setkey(dfpr, Patient)	
	
	setdiff(dfpr[, Patient], df.viro.allmsm[, Patient])
	#[1] "M11453" "M12803" "M15053" "M28707" "M29294" "M30202" "M31104" "M31555" "M31972" "M32892"
	#[11] "M33098" "M33508" "M34539" "M35680" "M39384"

	dfvp	<- merge(subset(dfpr, select=c(Patient, NegT, AnyPos_T1, AnyT_T1, PosSeqT, DateDied, DateBorn)), dfv, by='Patient')
	setkey(dfvp, Patient)
	dfvpu	<- unique(dfvp)
	setkey(dfvp, AnyPos_T1, Patient)
	
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_MSMtransmission_Primo/150417_Primo.R'
	save(dfvpu, dfvp, dfpr, file=file)
	
	ggplot(dfvpu, aes(x=AnyT_T1-AnyPos_T1)) + geom_histogram(binwidth=1/12) + scale_x_continuous(breaks=seq(0,15,1))
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_MSMtransmission_Primo/150417_TimeToARTStart.pdf'
	ggsave(file=file, w=8, h=6)
	
	subset(dfvpu, AnyPos_T1>2003.33 & AnyPos_T1<=2013.16)[, table(AnyT_T1-AnyPos_T1<.5, useNA='if')]
	
	#	ART within 6mo
	tmp		<- subset(dfvp, AnyT_T1-AnyPos_T1<.5)
	ggplot(tmp, aes(group=Patient)) +
			scale_y_continuous(breaks=seq(0, 7, 2), expand=c(0,0)) +				
			scale_x_continuous(breaks=seq(1980, 2015, 5), minor_breaks=seq(1980, 2015, 1), expand=c(0,0)) +
			coord_trans(limy=c(0,7)) +
			labs(title='AnyT_T1-AnyPos_T1<.5',x='', y='log10 viral load\n(cps/ml)') +
			geom_rect(data=merge(dfvpu, unique(subset(tmp, select=Patient)), by='Patient'), aes(xmin=AnyPos_T1, xmax=AnyT_T1, ymin=0, ymax=10), fill='orange', alpha=0.1) +
			geom_rect(data=merge(dfvpu, unique(subset(tmp, select=Patient)), by='Patient'), aes(xmin=AnyT_T1, xmax=DateDied, ymin=0, ymax=10), fill='green', alpha=0.1) +			
			geom_line(aes(x=PosRNA, y=lRNA, group=Patient), colour='black') +
			geom_point(aes(x=PosRNA, y=lRNA), colour='grey50') +
			theme_bw() +
			theme(panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.3)) +
			facet_grid(Patient~., scales='free')			
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_MSMtransmission_Primo/150417_VLtraj_ARTwithin6mo.pdf'
	ggsave(file=file, w=8, h=1*tmp[,length(unique(Patient))],limitsize=FALSE)
	
	#	ART after 6mo
	tmp		<- subset(dfvp, AnyT_T1-AnyPos_T1>=.5)
	ggplot(tmp, aes(group=Patient)) +
			scale_y_continuous(breaks=seq(0, 7, 2), expand=c(0,0)) +				
			scale_x_continuous(breaks=seq(1980, 2015, 5), minor_breaks=seq(1980, 2015, 1), expand=c(0,0)) +
			coord_trans(limy=c(0,7)) +
			labs(title='AnyT_T1-AnyPos_T1>=.5',x='', y='log10 viral load\n(cps/ml)') +
			geom_rect(data=merge(dfvpu, unique(subset(tmp, select=Patient)), by='Patient'), aes(xmin=AnyPos_T1, xmax=AnyT_T1, ymin=0, ymax=10), fill='orange', alpha=0.1) +
			geom_rect(data=merge(dfvpu, unique(subset(tmp, select=Patient)), by='Patient'), aes(xmin=AnyT_T1, xmax=DateDied, ymin=0, ymax=10), fill='green', alpha=0.1) +			
			geom_line(aes(x=PosRNA, y=lRNA, group=Patient), colour='black') +
			geom_point(aes(x=PosRNA, y=lRNA), colour='grey50') +
			theme_bw() +
			theme(panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.3)) +
			facet_grid(Patient~., scales='free')			
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_MSMtransmission_Primo/150417_VLtraj_ARTafter6mo.pdf'
	ggsave(file=file, w=8, h=1*tmp[,length(unique(Patient))],limitsize=FALSE)
}

primoSHM.probableTransmitters<- function()
{
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_MSMtransmission_Primo/150417_Primo.R'
	load(file)
		
	setnames(dfvpu, 'Patient', 't.Patient')
	YXpr	<- merge( YX, subset(dfvpu, select=c(t.Patient)), by='t.Patient' )
	setkey(YXpr, t.Patient)
	
	unique(YXpr)[, table(t.AnyT_T1-t.AnyPos_T1<.5, useNA='if')]
	#FALSE  TRUE 
	#13    12
	subset(dfpr, !is.na(PosSeqT))[, table(AnyT_T1-AnyPos_T1<.5, useNA='if')]
	#FALSE  TRUE 
	#76   136 
}