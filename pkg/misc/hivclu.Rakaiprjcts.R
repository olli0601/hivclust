######################################################################################
project.Rakai.processRegion1PANGEAalignment.151129<- function()
{
	susa.f	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151129.fasta'
	susa.s	<- read.dna(susa.f, format='fasta')
	susa.d	<- data.table(TAXA=rownames(susa.s), DATA= factor(grepl('^PG', rownames(susa.s)), levels=c(TRUE,FALSE), label=c('PNG','LNL')))
	
	png.f	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/PANGEA_orig/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	load(png.f)	#sq, sqi, si
	
	tmp		<- subset(sqi, is.na(SITE))[, PANGEA_ID]
	susa.d[, tmp%in%TAXA]	#Susanna did not keep the references?
	tmp		<- sapply(strsplit(tmp, '.', fixed=1), function(x) tail(x,1))
	sapply(tmp, function(x) susa.d[, any(grepl(x, TAXA))]	)	#no...
	
	#
	# find PANGEA sequence with greatest coverage in the PANGEA and Su alignment overlap
	#
	sqiu	<- subset(sqi, SITE=='UG')	
	tmp		<- subset(sqiu, COV==max(COV))[1, TAXA]
	susa.s[tmp,]	# empty!	
	#	3628 UG seqs in PANGEA alignment
	subset(susa.d, DATA=='PNG')
	#	2923 UG PANGEA in SU alignment	
	tmp		<- merge(susa.d, sqiu, by='TAXA')
	#	2903 of 3628 in SU alignment
	setdiff( subset(susa.d, DATA=='PNG')[, TAXA], sqiu[, TAXA])	# Susanna gave a different name to the duplicates TODO: Chris to change taxa names with ( and space
	setdiff( sqiu[, TAXA], subset(susa.d, DATA=='PNG')[, TAXA])	# quite a few more missing!
	#	find start of SU alignment in PANGEA alignment
	tmp2	<- susa.s[tmp[, TAXA],]
	tmp2	<- as.character(tmp2[,1:50])
	tmp2	<- data.table(TAXA=rownames(tmp2), COV50=apply(tmp2, 1, function(x) length(which(!x%in%c('-','?')))	))	
	mxcv	<- tmp2[which.max(COV50),TAXA]						# mxcov for key
	tmp2	<- gsub('?','',gsub('-','',paste(as.vector(as.character(susa.s[mxcv,1:50])), collapse='')),fixed=1)		#build SU alignment key
	tmp2	<- paste(strsplit(tmp2, '')[[1]],'-*',sep='',collapse='')
	startat	<- regexpr(tmp2, gsub('?','-',paste(as.vector(as.character(sq[mxcv,])), collapse=''),fixed=1))	
	#	TODO: double check that start of region1 is at PANGEA 1
	#	TODO: startat ignores - at start, so pos will not be right if startat != 1
	#	code below assumes that startat is 1
	
	#
	# calculate offset
	#	
	lanl.s	<- susa.s[subset(susa.d, DATA!='PNG')[, TAXA],]
	tmp2	<- as.character(lanl.s)
	tmp2	<- data.table(	TAXA= rownames(lanl.s), 
			FIRST= apply( tmp2, 1, function(x) which(x!='-')[1] ),
			LAST= ncol(lanl.s)-apply( tmp2, 1, function(x) which(rev(x)!='-')[1] )+1L		)
	lanl.d	<- merge(susa.d, tmp2, by='TAXA')	
	lanl.l	<- lanl.d[, max(LAST)]			#this is the final pos with partial seqs in the new alignment
	sqm		<- sq[tmp[, TAXA],seq.int(startat,lanl.l)]
	tmp2	<- as.character(sqm)
	tmp2	<- data.table(TAXA=rownames(tmp2), COV=apply(tmp2, 1, function(x) length(which(!x%in%c('-','?')))	))
	mxcv	<- tmp2[which.max(COV),TAXA]						# mxcov for offset calculations
	# calculate offset	
	y		<- susa.s[mxcv,]		#y must be the longer sequence with extra gaps
	x		<- sq[mxcv,]	
	x		<- as.character(x)	
	y		<- as.character(y)	
	stopat	<- 9e3
	#substring(paste(as.vector(x), collapse=''),280,400)
	#k<- 285
	#rbind( z[1, ((k-5):(k+20))+offset.in.x[(k-5):(k+20)]], z[2, ((k-5):(k+20))+offset.in.y[(k-5):(k+20)]] )
	#z[1, ((k-5):(k+20))]
	if(0)	#DEV
	{
		x<- matrix(strsplit('tcaa---------aaattt','')[[1]], nrow=1)
		y<- matrix(strsplit('-tcaaaa---a-ttt----','')[[1]], nrow=1)		
	}	
	stopifnot( gsub('?','',gsub('-','',paste(as.vector(x), collapse='')),fixed=1)==gsub('?','',gsub('-','',paste(as.vector(y), collapse='')),fixed=1) )
	z							<- rbind.fill.matrix(x,y)
	z[1,which(is.na(z[1,]))]	<- '-'
	z[2,which(is.na(z[2,]))]	<- '-'
	zz							<- unname(z)	
	offset.to.x	<- rep(0, ncol(z))
	offset.to.y	<- rep(0, ncol(z))			
	k			<- seq_len(ncol(z))+offset.to.y
	k			<- k[ k>0 & k<=ncol(z)]
	k			<- which( z[1, seq_along(k)]!=z[2, k] )[1]
	while(!is.na(k) & k<stopat)
	{
		#print(k)		
		done	<- 0
		if( !zz[1,k]%in%c('-','?') )
		{
			kn		<- k-offset.to.x[min(k,length(offset.to.x))]
			offset.to.x[ seq.int(kn,length(offset.to.x)) ] <- offset.to.x[ seq.int(kn,length(offset.to.x)) ]+1
			if(k==1)
				zz	<- rbind( c('-',zz[1,]), c(zz[2,],'-') )
			if(k>1)
				zz	<- rbind( c(zz[1, 1:(k-1)],'-',zz[1,k:ncol(zz)]), c(zz[2,],'-')	)
			done	<- 1
		}
		if( !done & zz[1,k]%in%c('-','?') )
		{
			kn		<- k-offset.to.y[min(k,length(offset.to.x))]
			offset.to.y[ seq.int(kn,length(offset.to.x)) ] <- offset.to.y[ seq.int(kn,length(offset.to.y)) ]+1
			if(k==1)
				zz	<- rbind( c(zz[1,],'-'), c('-',zz[2,]) )
			if(k>1)
				zz	<- rbind( c(zz[1,],'-'), c(zz[2, 1:(k-1)],'-',zz[2,k:ncol(zz)])	)							
		}
		k	<- unname(which(zz[1,]!=zz[2,])[1])		
	}
	if(0)
	{
		#	DEV
		# to x add offset.to.x
		k			<- offset.to.x+seq_along(offset.to.x)
		xx			<- matrix('-',ncol=max(k),nrow=nrow(x))
		xx[,k]		<- x[, seq_along(k)]
		# to y add offset.to.y
		k			<- offset.to.y+seq_along(offset.to.y)
		yy			<- matrix('-',ncol=max(k),nrow=nrow(y))
		yy[,k]		<- y[, seq_along(k)]
		rbind(xx, yy)
	}
	# add gaps to SU alignment
	tmp			<- as.character(lanl.s[, seq_len(lanl.l)])
	lanl.s2		<- rbind(tmp,as.character(susa.s[mxcv, seq_len(lanl.l)]))	#ONLY TO CHECK		
	k			<- offset.to.y[seq_len(lanl.l)]+seq_len(lanl.l)
	tmp			<- matrix('-',ncol=max(k),nrow=nrow(lanl.s2), dimnames=dimnames(lanl.s2))
	tmp[,k]		<- lanl.s2[, seq_along(k)]
	lanl.s2		<- as.DNAbin(tmp)
	# add gaps to PNGEA alignment 
	tmp			<- offset.to.x+seq_along(offset.to.x)
	k			<- tmp[tmp<=ncol(lanl.s2)]
	tmp			<- as.character(sq[,seq_along(k)])			#take PANGEA alignment and stretch
	sqn			<- matrix('-',nrow=nrow(tmp), ncol=max(k), dimnames=dimnames(tmp))
	sqn[,k]		<- tmp[,seq_along(k)]							
	sqn			<- as.DNAbin(sqn)
	# now add non-stretched part of PANGEA alignment
	sqn			<- cbind(sqn,sq[,seq.int(length(k)+1,ncol(sq))])
	# now bring LANL to same size as new alignment
	tmp			<- as.DNAbin( matrix('-',nrow=nrow(lanl.s2), ncol=ncol(sqn)-ncol(lanl.s2), dimnames=dimnames(lanl.s2)) )
	lanl.s2		<- cbind(lanl.s2, tmp)
	# now merge!
	sqn			<- rbind(sqn,lanl.s2)
	if(0)	
	{	
		#CHECK
		write.dna( sqn[which(mxcv==rownames(sqn)),], format='fasta', file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/check.fasta' )
		# regmatches(paste(as.vector(x), collapse=''),regexpr('ttttttagggagactt.*',paste(as.vector(x), collapse='')))		
	}
	#	rm Check
	sqn			<- sqn[-which(mxcv==rownames(sqn))[2],]
	#	get info
	sqni		<- data.table(TAXA=rownames(sqn), IDX=seq_len(nrow(sqn)))		
	sqni[, SRC:= factor(grepl('^PG',TAXA),levels=c(TRUE,FALSE),labels=c('PNG','LANL'))]
	tmp			<- sqni[, which(SRC=='LANL' & IDX<300)]
	set(sqni, tmp, 'SRC', 'COMPENDIUM')	
	tmp			<- as.character(sqn)
	tmp			<- data.table(	TAXA=rownames(tmp), 
								FIRST= apply( tmp, 1, function(x) which(!x%in%c('-','?'))[1] ),
								LAST= ncol(tmp)-apply( tmp, 1, function(x) which(!rev(x)%in%c('-','?'))[1] ) + 1L,						
								COV=apply(tmp, 1, function(x) length(which(!x%in%c('-','?')))	))
	sqni		<- merge(sqni,tmp,by='TAXA')
	tmp			<- subset(sqni, SRC=='LANL')[, max(LAST)]
	tmp			<- as.character(sqn[,seq_len(tmp)])
	tmp			<- data.table(	TAXA=rownames(tmp), 
								COV_REGION=apply(tmp, 1, function(x) length(which(!x%in%c('-','?')))	))
	sqni		<- merge(sqni,tmp,by='TAXA')
	sqni[, SUBT:=NA_character_]
	tmp			<- sqni[, which(SRC=='COMPENDIUM')]
	set(sqni, tmp, 'SUBT', sqni[tmp, sapply(strsplit(TAXA,'.',fixed=1),'[[',1)])
	tmp			<- sqni[, which(SRC=='LANL')]
	set(sqni, tmp, 'SUBT', sqni[tmp, sapply(strsplit(TAXA,'*',fixed=1),'[[',3)])
	tmp			<- sqni[, which(SUBT%in%c('-','U'))]
	set(sqni, tmp, 'SUBT', NA_character_)
	sqni[, ACCN:=NA_character_]
	tmp			<- sqni[, which(SRC=='COMPENDIUM')]
	set(sqni, tmp, 'ACCN', sqni[tmp, sapply(strsplit(TAXA,'.',fixed=1),function(x) rev(x)[1])])
	tmp			<- sqni[, which(SRC=='LANL')]
	set(sqni, tmp, 'ACCN', sqni[tmp, sapply(strsplit(TAXA,'*',fixed=1),'[[',1)])
	write.dna( sqn, format='fasta', file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201.fasta', colsep='', nbcol=-1)
	save(sqni, sqn, file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201.rda')
	
	sqn				<- read.dna(file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201_rm99gps.fasta',format='fasta')
	rownames(sqn)	<- gsub(' (stripped)','',rownames(sqn),fixed=1)
	#	get info
	sqni		<- data.table(TAXA=rownames(sqn), IDX=seq_len(nrow(sqn)))		
	sqni[, SRC:= factor(grepl('^PG',TAXA),levels=c(TRUE,FALSE),labels=c('PNG','LANL'))]
	tmp			<- sqni[, which(SRC=='LANL' & IDX<300)]
	set(sqni, tmp, 'SRC', 'COMPENDIUM')	
	tmp			<- as.character(sqn)
	tmp			<- data.table(	TAXA=rownames(tmp), 
			FIRST= apply( tmp, 1, function(x) which(!x%in%c('-','?'))[1] ),
			LAST= ncol(tmp)-apply( tmp, 1, function(x) which(!rev(x)%in%c('-','?'))[1] ) + 1L,						
			COV=apply(tmp, 1, function(x) length(which(!x%in%c('-','?')))	))
	sqni		<- merge(sqni,tmp,by='TAXA')
	tmp			<- subset(sqni, SRC=='LANL')[, max(LAST)]
	tmp			<- as.character(sqn[,seq_len(tmp)])
	tmp			<- data.table(	TAXA=rownames(tmp), 
			COV_REGION=apply(tmp, 1, function(x) length(which(!x%in%c('-','?')))	))
	sqni		<- merge(sqni,tmp,by='TAXA')
	sqni[, SUBT:=NA_character_]
	tmp			<- sqni[, which(SRC=='COMPENDIUM')]
	set(sqni, tmp, 'SUBT', sqni[tmp, sapply(strsplit(TAXA,'.',fixed=1),'[[',1)])
	tmp			<- sqni[, which(SRC=='LANL')]
	set(sqni, tmp, 'SUBT', sqni[tmp, sapply(strsplit(TAXA,'*',fixed=1),'[[',3)])
	tmp			<- sqni[, which(SUBT%in%c('-','U'))]
	set(sqni, tmp, 'SUBT', NA_character_)
	sqni[, ACCN:=NA_character_]
	tmp			<- sqni[, which(SRC=='COMPENDIUM')]
	set(sqni, tmp, 'ACCN', sqni[tmp, sapply(strsplit(TAXA,'.',fixed=1),function(x) rev(x)[1])])
	tmp			<- sqni[, which(SRC=='LANL')]
	set(sqni, tmp, 'ACCN', sqni[tmp, sapply(strsplit(TAXA,'*',fixed=1),'[[',1)])
	save(sqni, sqn, file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201_rm99gps.rda')	
	write.dna( sqn, format='fasta', file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201_rm99gps.fasta', colsep='', nbcol=-1)
}
######################################################################################
project.Rakai.checkforARVs.150910<- function()
{
	require(data.table)
	#	Kate, I generally use data.table rather than data.frame since it s faster and more versatile. 
	#	warning: it s a bit of a learning curve though to use it
	
	require(big.phylo)
	#	library(devtools)
	#	install_github("olli0601/big.phylo")
	
	f.arv	<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_150909.rda'
	f.rccsid<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/150625_PangeaBoxes.csv'
	f.sid	<- '~/Dropbox (Infectious Disease)/PANGEA_data/SangerUpdates/2015-07-20_PANGEA_3.csv'
	f.seq	<- '~/Dropbox (Infectious Disease)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908.fasta'
	#
	#	get IDs into OK format
	#
	load(f.arv)	#	loads arvdat
	arvdat	<- as.data.table(arvdat)
	d.rccsid<- as.data.table(read.csv(f.rccsid, stringsAsFactors=0))
	d.rccsid<- subset(d.rccsid, select=c(Pangea.ID, Rakai.Study.ID))
	setnames(d.rccsid, c('Rakai.Study.ID','Pangea.ID'), c('RCCS_ID','PNG_ID'))
	set(d.rccsid, NULL, 'RCCS_ID', d.rccsid[, substr(RCCS_ID, 1, nchar(RCCS_ID)-3)])
	set(d.rccsid, NULL, 'PNG_ID', d.rccsid[, gsub('-Sy*','',PNG_ID)])
	setnames(arvdat, 'RCCS_studyid', 'RCCS_ID')
	d.sid	<- as.data.table(read.csv(f.sid, stringsAsFactors=0, header=0))
	setnames(d.sid, c('V1','V2','V3'), c('PNG_ID_FULL', 'SNG_ID', 'x'))
	d.sid[, PNG_ID:= gsub('-S[0-9]+','',PNG_ID_FULL)]
	set(d.sid, NULL, 'SNG_ID', d.sid[, gsub('#','_',SNG_ID)])
	d.sid[, x:=NULL]
	#
	#	merge
	#	
	subset(arvdat, any.pangea==1)	
	#	126		
	arvdat	<- merge( arvdat, d.rccsid, by='RCCS_ID' )
	#	85
	arvdat	<- merge( arvdat, d.sid, by='PNG_ID')
	#	76
	
	#
	#	check for DRMs
	#
	seq				<- read.dna(file=f.seq, format='fasta')	
	seq				<- seq[c("B.FR.83.HXB2_LAI_IIIB_BRU.K03455",arvdat[, PNG_ID_FULL]),]
	rownames(seq)[1]	<- 'HXB2'
	outfile			<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_PotentialDRMs_OR_150910.R'
	tmp				<- big.phylo:::seq.rm.drugresistance(seq, outfile=outfile)
	nodr.info		<- tmp$nodr.info
	nodr.seq		<- tmp$nodr.seq
}
######################################################################################
project.Rakai.checkforARVs.150911<- function()
{
	require(data.table)
	#	Kate, I generally use data.table rather than data.frame since it s faster and more versatile. 
	#	warning: it s a bit of a learning curve though to use it
	
	require(big.phylo)
	#	library(devtools)
	#	install_github("olli0601/big.phylo")
	
	f.arv			<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_150911.rda'
	f.rccsid		<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/150625_PangeaBoxes.csv'
	f.sid			<- '~/Dropbox (Infectious Disease)/PANGEA_data/SangerUpdates/2015-07-20_PANGEA_3.csv'
	f.seq			<- '~/Dropbox (Infectious Disease)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908.fasta'
	#
	load(f.arv)	#	loads rakaidat
	rakaidat		<- as.data.table(rakaidat)
	setnames(rakaidat, 'ProjectID', 'PNG_ID')
	#	482
	seq				<- read.dna(file=f.seq, format='fasta')	
	rownames(seq)	<- gsub('-S[0-9]+','',rownames(seq))
	d.seq			<- merge( data.table(PNG_ID=rownames(seq)), rakaidat, by='PNG_ID' )
	#	459 with PANGEA sequence
	
	
	#
	#	check for DRMs
	#
	seq				<- seq[c("B.FR.83.HXB2_LAI_IIIB_BRU.K03455",d.seq[, PNG_ID]),]		
	rownames(seq)[1]	<- 'HXB2'
	outfile			<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_PotentialDRMs_OR_150911.R'
	tmp				<- big.phylo:::seq.rm.drugresistance(seq, outfile=outfile)
	nodr.info		<- tmp$nodr.info
	nodr.seq		<- tmp$nodr.seq
}