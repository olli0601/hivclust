seattle.wrapper<- function()
{
	seattle.170621.fastree()
}

seattle.170601.rm.drug.resistance.mutations<- function()
{
	require(big.phylo)
	#	handle DRMs
	infile.fasta		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.fasta'
	outfile				<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm.fasta'
	sq					<- read.dna(infile.fasta, format='fa')	
	tmp					<- which(grepl("HXB2",rownames(sq)))
	rownames(sq)[tmp]	<- 'HXB2'
	tmp					<- big.phylo:::seq.rm.drugresistance(sq)
	sq					<- tmp$nodr.seq
	write.dna(sq, file= outfile, format='fasta')
}

seattle.170607.rm.drug.resistance.mutations<- function()
{
	require(big.phylo)
	#	handle DRMs
	infile.fasta		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.plus.LANL.USA.refse.fasta'
	outfile				<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.plus.LANL.USA.refse.ndrm.fasta'	
	infile.fasta		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.PHSKCLANLe.refse.fasta'
	outfile				<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.PHSKCLANLe.refse.ndrm.fasta'
	sq					<- read.dna(infile.fasta, format='fa')	
	tmp					<- which(grepl("HXB2",rownames(sq)))
	rownames(sq)[tmp]	<- 'HXB2'
	tmp					<- big.phylo:::seq.rm.drugresistance(sq)
	sq					<- tmp$nodr.seq
	write.dna(sq, file= outfile, format='fasta')
}

seattle.170601.fastree<- function()
{	
	require(big.phylo)
	
	#	first tree
	infile.fasta	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.alignment.fasta'
	outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
	tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
	#	second tree, DRMs removed + including references for rooting
	infile.fasta	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm.fasta'
	outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
	tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
	#	run on command line
	cat(tmp)		
	
	#
	# 	bootstrap trees
	
	#	second tree, DRMs removed + including references for rooting
	infile.fasta	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm.fasta'
	bs.id			<- 2
	outfile.ft		<- gsub('\\.fasta',paste0('_ft.',sprintf("%03d",bs.id),'.newick'),infile.fasta)
	tmp				<- cmd.fasttree.one.bootstrap(infile.fasta, bs.id, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
	#	run on command line
	cat(tmp)		
	
	infile.fasta	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm.fasta'
	bs.n			<- 3
	bs.dir			<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm_bootstrap_trees'
	outfile			<- paste0('/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm_ft_bs',bs.n,'.newick')
	tmp				<- cmd.fasttree.many.bootstraps(infile.fasta, bs.dir, bs.n, outfile)
	cat(paste0('#!/bin/sh\n',tmp), file=paste0(outfile,'.sh'))
	Sys.chmod(paste0(outfile,'.sh'), mode = "777")	
}

seattle.170607.fastree<- function()
{	
	require(big.phylo)
	
	infile.fasta	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.plus.LANL.USA.refse.ndrm.fasta'
	outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
	tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma', check.binary=TRUE)
	#	run on command line
	cat(tmp)	
}

seattle.170621.fastree<- function()
{	
	require(big.phylo)
		
	#indir			<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle'
	indir			<- '/work/or105/SEATTLE'
	infile.fasta	<- file.path(indir,'Seattle.June2017.PHSKCLANLe.refse.ndrm.fasta')
	
	bs.n			<- 100	
	bs.dir			<- file.path(indir,'Seattle.June2017.PHSKCLANLe.refse.ndrm_bootstrap_trees')	
	outfile.ft		<- gsub('\\.fasta',paste0('_ft_bs',bs.n,'.newick'),infile.fasta)
	tmp				<- cmd.fasttree.many.bootstraps(infile.fasta, bs.dir, bs.n, outfile.ft, pr.args='-nt -gtr -gamma', opt.bootstrap.by='nucleotide')
	
	#	run on HPC
	cmd				<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=998, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=1, hpc.load='module load R/3.2.0')
	cmd				<- paste(cmd,tmp,sep='\n')
	cat(cmd)					
	outfile.cmd		<- paste("sea",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
	cmd.hpccaller(indir, outfile.cmd, cmd)	
}

seattle.170814.metadata.regionorigin<- function()
{
	infile	<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.2017_08_04.metadata/sequences_meta.RData'
	outfile	<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.2017_08_04.metadata/sequences_first_OR.rda'
	load(infile)
	ds	<- as.data.table(sequences_meta)
	ds	<- subset(ds, select=c(seqID, newnum, seqy, seqm))
	setnames(ds, colnames(ds), toupper(colnames(ds)))
	set(ds, NULL, 'NEWNUM', ds[, as.integer(NEWNUM)])
	set(ds, NULL, 'SEQID', ds[, as.character(SEQID)])
	set(ds, NULL, 'SEQY', ds[, as.integer(SEQY)])
	set(ds, NULL, 'SEQM', ds[, as.integer(SEQM)])
	set(ds, NULL, 'SEQ_DATE', ds[, as.numeric(SEQY)+(as.numeric(SEQM)-1+.5)/12])
	
	dsf	<- ds[, {
					z<- which(SEQ_DATE==min(SEQ_DATE))
					list(SEQID=SEQID[z], SEQ_DATE1=min(SEQ_DATE))
				}, by='NEWNUM']
	dsf[, GENE:= dsf[,gsub('(.*)-[0-9]+$','\\1',SEQID)]]	
	save(ds, dsf, file=outfile)
	
	infile	<- '~/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.2017_08_04.metadata/person.RData'
	outfile	<- gsub('\\.RData','_OR.rda',infile)
	load(infile)
	dp	<- as.data.table(person)
	dp	<- subset(dp, select=c(newnum, race, b_yr, gender, transm, birthCountry, rsh_county_name, rsa_county_name, cur_county_name, hivyr, hivmo))
	setnames(dp, colnames(dp), toupper(colnames(dp)))
	
	set(dp, NULL, 'NEWNUM', dp[, as.integer(NEWNUM)])
	set(dp, NULL, 'RACE', dp[, as.character(RACE)])
	set(dp, NULL, 'B_YR', dp[, as.integer(B_YR)])
	set(dp, NULL, 'GENDER', dp[, as.character(GENDER)])
	set(dp, NULL, 'TRANSM', dp[, as.character(TRANSM)])
	set(dp, NULL, 'BIRTHCOUNTRY', dp[, as.character(BIRTHCOUNTRY)])
	set(dp, NULL, 'RSH_COUNTY_NAME', dp[, as.character(RSH_COUNTY_NAME)])
	set(dp, NULL, 'RSA_COUNTY_NAME', dp[, as.character(RSA_COUNTY_NAME)])
	set(dp, NULL, 'CUR_COUNTY_NAME', dp[, as.character(CUR_COUNTY_NAME)])
	set(dp, NULL, 'HIVYR', dp[, as.character(HIVYR)])
	set(dp, NULL, 'HIVMO', dp[, as.character(HIVMO)])
	
	set(dp, NULL, 'RACE', dp[, gsub('\\([0-9]+\\)(.*)','\\1',RACE)])
	set(dp, NULL, 'RACE', dp[, gsub('-','_',gsub('/','_',gsub(' ','_',gsub('Not Hispanic, |, All races','',RACE))))])
	set(dp, NULL, 'TRANSM', dp[, gsub('[0-9]+\\.(.*)','\\1',TRANSM)])	
	#set(dp, NULL, 'HIV_DIAG', dp[, as.numeric(HIVYR)+(as.numeric(HIVMO)-1+.5)/12])
	
	# use line below to create the following vector
	# dp[, cat(paste(unique(BIRTHCOUNTRY), collapse='\", \"'))]
	tmp	<- c( "(USA)United States"="USA", 
			"(ZMB)Zambia"="AFR", 
			"(SLV)El Salvador"="CAR", 
			"(X99)Not Specified/Unknown"="UNKNOWN", 
			"(MEX)Mexico"="MEX", 
			"(VNM)Viet Nam"="SASIA", 
			"(TZA)Tanzania, United Republic of"="AFR", 
			"(MWI)Malawi"="AFR", 
			"(PHL)Philippines"="SASIA", 
			"(IRL)Ireland"="EU", 
			"(PER)Peru"="SAM", 
			"(ETH)Ethiopia"="AFR", 
			"(KEN)Kenya"="AFR", 
			"(PRI)Puerto Rico"="CAR", 
			"(GRC)Greece"="EU", 
			"(GTM)Guatemala"="CAR", 
			"(CUB)Cuba"="CAR", 
			"(X98)Other"="UNKNOWN", 
			"(NIC)Nicaragua"="CAR", 
			"(CIV)Cote d'Ivoire (Ivory Coast)"="AFR", 
			"(ARG)Argentina"="SAM", 
			"(ESP)Spain"="EU", 
			"(ERI)Eritrea"="AFR", 
			"(HKG)Hong Kong"="SASIA", 
			"(VEN)Venezuela"="SAM", 
			"(LBR)Liberia"="AFR", 
			"(SOM)Somalia"="AFR", 
			"(DOM)Dominican Republic"="CAR", 
			"(ITA)Italy"="EU", 
			"(THA)Thailand"="SASIA", 
			"(UGA)Uganda"="AFR", 
			"(GUF)French Guiana"="SAM", 
			"(JAM)Jamaica"="CAR", 
			"(CAN)Canada"="CAN", 
			"(DMA)Dominica"="CAR", 
			"(MRT)Mauritania"="AFR", 
			"(EGY)Egypt"="AFR", 
			"(CHN)China"="NASIA", 
			"(JPN)Japan"="NASIA", 
			"(DEU)Germany"="EU", 
			"(ZWE)Zimbabwe"="AFR", 
			"(PAN)Panama"="CAR", 
			"(TWN)Taiwan, Province of China"="NASIA", 
			"(IND)India"="NASIA", 
			"(GUM)Guam"="PAC", 
			"(KOR)Korea, Republic of (South)"="NASIA", 
			"(LAO)Lao People's Democratic Republic"="SASIA", 
			"(KHM)Cambodia"="SASIA", 
			"(GBR)United Kingdom"="EU", 
			"(ASM)American Samoa"="PAC", 
			"(FJI)Fiji"="PAC", 
			"(RUS)Russian Federation"="NASIA", 
			"(SAU)Saudi Arabia"="MIDEAST", 
			"(TGO)Togo"="AFR", 
			"(COD)Congo, Democratic Republic of Zaire"="AFR", 
			"(SLE)Sierra Leone"="AFR", 
			"(MMR)Myanmar (Burma)"="SASIA", 
			"(NGA)Nigeria"="AFR", 
			"(COL)Colombia"="SAM", 
			"(COG)Congo"="AFR", 
			"(NOR)Norway"="EU", 
			"(SDN)Sudan"="AFR", 
			"(SVK)Slovakia"="EU", 
			"(TON)Tonga"="PAC", 
			"(BRA)Brazil"="SAM", 
			"(POL)Poland"="EU", 
			"(LUX)Luxembourg"="EU", 
			"(NZL)New Zealand"="PAC", 
			"(AFG)Afghanistan"="NASIA", 
			"(MLI)Mali"="AFR", 
			"(RWA)Rwanda"="AFR", 
			"(CHL)Chile"="SAM", 
			"(BDI)Burundi"="AFR", 
			"(ECU)Ecuador"="CAR", 
			"(GHA)Ghana"="AFR", 
			"(FRA)France"="EU", 
			"(AUS)Australia"="PAC", 
			"(HND)Honduras"="CAR", 
			"(SYR)Syrian Arab Republic"="MIDEAST", 
			"(HTI)Haiti"="CAR", 
			"(MNP)Northern Mariana Islands"="PAC", 
			"(IDN)Indonesia"="SASIA", 
			"(ZAF)South Africa"="AFR", 
			"(CMR)Cameroon"="AFR", 
			"(LKA)Sri Lanka"="NASIA", 
			"(KAZ)Kazakhstan"="NASIA", 
			"(IRN)Iran, Islamic Republic of"="NASIA", 
			"(VIR)Virgin Islands, U.S."="CAR", 
			"(MHL)Marshall Islands"="PAC", 
			"(GIN)Guinea"="AFR", 
			"(MAR)Morocco"="AFR", 
			"(ROU)Romania"="EU", 
			"(YUG)Yugoslavia"="EU", 
			"(CRI)Costa Rica"="SAM", 
			"(MYS)Malaysia"="SASIA", 
			"(GMB)Gambia"="AFR", 
			"(YEM)Yemen"="MIDEAST", 
			"(SWZ)Swaziland"="AFR", 
			"(UKR)Ukraine"="EU", 
			"(AUT)Austria"="EU", 
			"(TCD)Chad"="AFR", 
			"(NLD)Netherlands"="EU", 
			"(BFA)Burkina Faso"="AFR", 
			"(URY)Uruguay"="SAM", 
			"(NAM)Namibia"="AFR", 
			"(NER)Niger"="AFR", 
			"(SGP)Singapore"="SASIA", 
			"(LBN)Lebanon"="MIDEAST", 
			"(PRT)Portugal"="EU", 
			"(BIH)Bosnia and Herzegovina"="EU",
			"(BEN)Benin"="AFR", 
			"(ISR)Israel"="MIDEAST", 
			"(HUN)Hungary"="EU", 
			"(MDA)Moldova, Republic of"="EU", 
			"(BGR)Bulgaria"="EU", 
			"(SEN)Senegal"="AFR", 
			"(DNK)Denmark"="EU", 
			"(LSO)Lesotho"="AFR", 
			"(BHR)Bahrain"="MIDEAST", 
			"(ALB)Albania"="EU", 
			"(UZB)Uzbekistan"="NASIA", 
			"(TTO)Trinidad and Tobago"="CAR", 
			"(BWA)Botswana"="AFR", 
			"(PAK)Pakistan"="NASIA", 
			"(IRQ)Iraq"="NASIA", 
			"(DJI)Djibouti"="AFR", 
			"(BOL)Bolivia"="SAM", 
			"(ATG)Antigua and Barbuda"="CAR", 
			"(MOZ)Mozambique"="AFR", 
			"(GAB)Gabon"="AFR", 
			"(BEL)Belgium"="EU", 
			"(BMU)Bermuda"="CAR", 
			"(LVA)Latvia"="EU", 
			"(LBY)Libyan Arab Jamahiriya"="AFR", 
			"(FSM)Micronesia, Federated States of"="PAC", 
			"(MNG)Mongolia"="NASIA", 
			"(SWE)Sweden"="EU", 
			"(TUR)Turkey"="MIDEAST", 
			"(GEO)Georgia"="EU", 
			"(WSM)Samoa"="PAC", 
			"(AGO)Angola"="AFR", 
			"(PNG)Papua New Guinea"="SASIA", 
			"(GUY)Guyana"="SAM", 
			"(AZE)Azerbaijan"="NASIA", 
			"(NPL)Nepal"="NASIA", 
			"(CZE)Czech Republic"="EU", 
			"(DZA)Algeria"="AFR", 
			"(PRY)Paraguay"="SAM", 
			"(QAT)Qatar"="MIDEAST", 
			"(BLZ)Belize"="SAM", 
			"(MTQ)Martinique"="CAR", 
			"(FIN)Finland"="EU", 
			"(BRB)Barbados"="CAR", 
			"(ISL)Iceland"="EU", 
			"(PRK)Korea, Dem People's Rep of (North)"="NASIA", 
			"(CHE)Switzerland"="EU", 
			"(UMI)U.S. Minor Outlying Areas"="PAC")
	dp	<- merge(dp, data.table(BIRTHCOUNTRY=names(tmp), BIRTHLOC=tmp), all.x=1, by='BIRTHCOUNTRY')
	set(dp, dp[, which(is.na(BIRTHLOC))], 'BIRTHLOC', 'UNKNOWN')
	save(dp, file=outfile)
}

seattle.170814.LSD<- function()
{	
	require(data.table)
	require(big.phylo)
	#
	#	write the overall dates file
	#
	if(0)
	{
		infile.seq	<- '~/Box Sync/OR_Work/2017/2017_Seattle/Seattle.2017_08_04.metadata/sequences_first_OR.rda'
		infile.pers	<- "~/Box Sync/OR_Work/2017/2017_Seattle/Seattle.2017_08_04.metadata/person_OR.rda"
		infile.tree	<- '~/Box Sync/OR_Work/2017/2017_Seattle/Seattle.plus.LANL.USA/Seattle.plus.LANL.USA.refse.ndrm_ft_reroot.newick'
		load(infile.seq)
		load(infile.pers)
		
		ph			<- read.tree(infile.tree)
		dph			<- data.table(TAXA=ph$tip.label)
		dph[, SEA:= as.numeric(grepl('^SEA',TAXA))]
		tmp			<- dph[, which(SEA==1)]
		set(dph, tmp, 'NEWNUM', as.integer(dph[tmp, gsub('SEATTLE\\.([0-9]+)\\|.*','\\1',TAXA)]))
		#	there is one sequence per individual
		dph			<- merge(dph, unique(subset(dsf, select=c(NEWNUM, SEQ_DATE1))), by='NEWNUM', all.x=TRUE)
		
		write.csv(subset(dph, !is.na(NEWNUM) & is.na(SEQ_DATE1)), file='~/Box Sync/OR_Work/2017/2017_Seattle/Seattle.plus.LANL.USA/Seattle.plus.LANL.USA.refse.ndrm_ft_reroot_missingsamplingdates.csv')
		
		
		infile.dates	<- '~/Dropbox (Infectious Disease)/2017_NL_Introductions/seq_info/Geneflow/NONB_flowinfo_OR.csv'
		df				<- as.data.table(read.csv(infile.dates))
		df				<- subset(df, select=c(ID, SUBTYPE, SAMPLING_DATE))
		setnames(df, c('ID','SAMPLING_DATE'), c('TAXA','DATE'))
		tmp				<- copy(df)
		set(tmp, NULL, 'TAXA', tmp[,gsub('[0-9]$','',paste0(TAXA,'_subtype',SUBTYPE))])
		df				<- rbind(df, tmp)	
		infile.dates	<- gsub('\\.csv','_lsddates.csv',infile.dates)
		write.csv(df, row.names=FALSE, infile.dates)
		#
		#	write subtype specific dates files
		#
		indir.ft		<- '~/Dropbox (Infectious Disease)/2017_NL_Introductions/trees_ft'
		infiles			<- data.table(F=list.files(indir.ft, pattern='*newick$', full.names=TRUE))
		infiles[, DATES_FILE:= gsub('_ft_bs100\\.newick','_lsddates.csv',F)]
		invisible(infiles[, {
							cmd		<- cmd.lsd.dates(infile.dates, F, DATES_FILE, run.lsd=FALSE)
							system(cmd)
							NULL
						}, by='F'])
	}	
	#
	#	run LSD 
	#
	#indir.ft		<- '~/Dropbox (Infectious Disease)/2017_NL_Introductions/trees_ft'
	#indir.dates		<- '~/Dropbox (Infectious Disease)/2017_NL_Introductions/trees_ft'
	#outdir			<- '~/Dropbox (Infectious Disease)/2017_NL_Introductions/trees_ft'
	
	indir.ft		<- '/work/or105/ATHENA_2016/vlintros/trees_ft'
	indir.dates		<- '/work/or105/ATHENA_2016/vlintros/trees_lsd'
	outdir			<- '/work/or105/ATHENA_2016/vlintros/trees_lsd'
	
	ali.len			<- 1000				
	lsd.args		<- '-v 1 -c -b 10'			# no rooting, no re-estimation of rates				
	#	get files	
	infiles	<- data.table(F=list.files(indir.ft, pattern='*newick$', full.names=TRUE, recursive=TRUE))
	infiles	<- subset(infiles, !grepl('RAxML',F))
	infiles[, SUBTYPE:= gsub('_.*','',basename(F))]
	tmp		<- data.table(FD=list.files(indir.dates, pattern='*_lsddates.csv$', full.names=TRUE, recursive=TRUE))
	tmp[, SUBTYPE:= gsub('_.*','',basename(FD))]
	infiles	<- merge(infiles, tmp, by='SUBTYPE')
	infiles[, FL:= file.path(outdir,gsub('_ft\\.|_ft_','_lsd_',gsub('\\.newick','',basename(F))))]	
	#	build LSD commands
	dlsd	<- infiles[, {
				cmd		<- cmd.lsd(F, FD, ali.len, outfile=FL, pr.args=lsd.args)
				list(CMD=cmd)
			}, by='F']
	tmp		<- dlsd[, paste(CMD, collapse='\n')]
	#	qsub
	cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=50, hpc.q="pqeelab", hpc.mem="5800mb",  hpc.nproc=1, hpc.load='module load R/3.3.3')							
	cmd		<- paste(cmd,tmp,sep='\n')
	#cat(cmd)					
	outfile	<- paste("scRAr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
	cmd.hpccaller(indir.ft, outfile, cmd)	
}


seattle.170601.clustering<- function()
{	
	require(big.phylo)
	#	first tree
	infile.ft		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm_ft.newick'
	ph				<- read.tree(infile.ft)
	#	reroot at AY371169
	tmp				<- which(grepl('AY371169',ph$tip.label))
	ph				<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
	ph 				<- ladderize( ph )
	#	process bootstrap support values	
	tmp				<- as.numeric( ph$node.label )		
	tmp[is.na(tmp)]	<- 0
	ph$node.label	<- tmp
	#	
	write.tree(ph, gsub('\\.newick','_reroot.newick',infile.ft))
	pdf(gsub('newick','pdf',infile.ft), w=15, h=350)
	plot(ph, cex=0.4, show.node.label=TRUE)
	dev.off()
	#
	#read patristic distances
	require(hivclust)
	stat.fun						<- hivc.clu.min.transmission.cascade	
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)	
	thresh.brl						<- 0.045
	thresh.bs						<- 0.8
	#produce clustering 
	clustering						<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph$node.label, retval="all")	
	outfile							<- gsub('\\.newick',paste0('_reroot.clu',thresh.brl*1e2,'-',thresh.bs*1e2,'.pdf'),infile.ft)	
	hivc.clu.plot(ph, 
					clustering[["clu.mem"]], 
					file=outfile, 
					pdf.off=0, pdf.width=15, pdf.height=1000, 
					show.node.label=TRUE, show.edge.label=TRUE,  
					cex.edge.incluster=2)
	dev.off()
	
	dfc		<- data.table(TAXA=ph$tip.label, CLU=clustering$clu.mem[1:Ntip(ph)])
	outfile	<- gsub('\\.newick',paste0('_reroot.clu',thresh.brl*1e2,'-',thresh.bs*1e2,'.rda'),infile.ft)
	save(ph, clustering, dfc, file, file=outfile)	
}

seattle.170607.clustering<- function()
{	
	require(big.phylo)
	#	first tree
	infile.ft		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.plus.LANL.USA.refse.ndrm_ft.newick'
	ph				<- read.tree(infile.ft)
	#	reroot at AY371169
	tmp				<- which(grepl('AY371169',ph$tip.label))
	ph				<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
	ph 				<- ladderize( ph )
	#	process bootstrap support values	
	tmp				<- as.numeric( ph$node.label )		
	tmp[is.na(tmp)]	<- 0
	ph$node.label	<- tmp
	#	
	write.tree(ph, gsub('\\.newick','_reroot.newick',infile.ft))
	pdf(gsub('\\.newick','_reroot.pdf',infile.ft), w=15, h=2e3)
	plot(ph, cex=0.4, show.node.label=TRUE)
	dev.off()
	#
	#read patristic distances
	require(hivclust)
	stat.fun						<- hivc.clu.min.transmission.cascade	
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)	
	thresh.brl						<- 0.045
	thresh.bs						<- 0.8
	#produce clustering 
	clustering						<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph$node.label, retval="all")	
	outfile							<- gsub('\\.newick',paste0('_reroot.clu',thresh.brl*1e2,'-',thresh.bs*1e2,'.pdf'),infile.ft)	
	hivc.clu.plot(ph, 
			clustering[["clu.mem"]], 
			file=outfile, 
			pdf.off=0, pdf.width=15, pdf.height=1000, 
			show.node.label=TRUE, show.edge.label=TRUE,  
			cex.edge.incluster=2)
	dev.off()
	
	dfc		<- data.table(TAXA=ph$tip.label, CLU=clustering$clu.mem[1:Ntip(ph)])
	outfile	<- gsub('\\.newick',paste0('_reroot.clu',thresh.brl*1e2,'-',thresh.bs*1e2,'.rda'),infile.ft)
	save(ph, clustering, dfc, file, file=outfile)	
}

seattle.170601.clustering.plot<- function()
{
	
	require(data.table)
	infile	<- "/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm_ft_reroot.clu4.5-80.rda"
	load(infile)
	#	exclude reference sequences non-Seattle
	dfc		<- subset(dfc, !grepl('AY371169|HXB2',TAXA))
	#	read meta-data 
	dfc[, ID:= sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',1)]
	dfc[, SXO:= sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',2)]
	dfc[, RACE:= sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',3)]
	dfc[, SEQ_YR:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',4))]
	dfc[, SEQ_Q:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',5))]
	dfc[, HIV_YR:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',6))]
	dfc[, BIRTH_YR:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',7))]
	dfc[, AGE_AT_SEQ:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',8))]
	dfc[, AGE_AT_HIV:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',9))]
	#	fixup
	set(dfc, dfc[, which(SXO=='UNKNOWN')], 'SXO', NA_character_)
	set(dfc, dfc[, which(SXO=='')], 'SXO', NA_character_)
	set(dfc, dfc[, which(SXO=='BLOOD')], 'SXO', 'Other')
	set(dfc, dfc[, which(SXO=='PERINAT')], 'SXO', 'Other')
	set(dfc, dfc[, which(RACE=='FALSE')], 'RACE', NA_character_)
	#	give each singleton a new cluster ID (cluster of size 1)
	tmp	<- dfc[, which(is.na(CLU))]
	set(dfc, tmp, 'CLU', max(dfc$CLU, na.rm=TRUE)+seq_along(tmp))
	
	#	compute cluster statistics 
	tmp	<- dfc[, list(	CLU_TS= min(HIV_YR, na.rm=TRUE),
				CLU_TE= max(HIV_YR, na.rm=TRUE),
				CLU_TM=mean(HIV_YR, na.rm=TRUE),
				CLU_N=length(ID),
				CLU_N_MSM=length(which(grepl('MSM',SXO))),
				CLU_N_IDU=length(which(grepl('IDU',SXO))),
				CLU_N_HET=length(which(grepl('HET',SXO)))
				), by='CLU']
	dfc	<- merge(dfc, tmp, by='CLU') 
	#	determine max lkl SXO of cluster
	tmp	<- melt(dfc, id.vars='CLU',measure.vars=c('CLU_N_MSM','CLU_N_IDU','CLU_N_HET'))
	tmp	<- tmp[, list(CLU_MX_SXO=variable[which.max(value)]), by='CLU']
	set(tmp, NULL, 'CLU_MX_SXO', tmp[, gsub('CLU_N_','',CLU_MX_SXO)])
	dfc	<- merge(dfc, tmp, by='CLU')
	#	set factors 
	set(dfc, NULL, 'CLU_MX_SXO', dfc[, factor(CLU_MX_SXO, levels=c('MSM','HET','IDU'), labels=c('MSM','heterosexual','injecting\ndrug\nuser'))])
	set(dfc, dfc[, which(is.na(SXO))], 'SXO', 'Unknown')
	set(dfc, dfc[, which(SXO=='HETEROp')], 'SXO', 'HETERO')
	set(dfc, NULL, 'SXO', dfc[, factor(SXO,  levels=c('MSM','MSM_IDU','IDU','HETERO','Other','Unknown'), 
											 labels=c('MSM','MSM/injecting drug user','injecting drug user','heterosexual','other','unknown'))])
	#
	#	plot clusters only
	dfp	<- subset(dfc, CLU_N>1 & !is.na(CLU_TS))
	#	order by CLU_MX_SXO and start date of cluster
	tmp	<- unique(subset(dfp, select=c(CLU, CLU_MX_SXO, CLU_TS)))
	tmp	<- tmp[order(CLU_MX_SXO, CLU_TS), ]	
	set(tmp, NULL, 'CLU_ORDERED', tmp[, factor(CLU, levels=tmp$CLU)])
	dfp	<- merge(dfp, subset(tmp, select=c(CLU, CLU_ORDERED)), by='CLU')
	tmp	<- melt(unique(subset(dfp, select=c(CLU_ORDERED, CLU_TS, CLU_TE, CLU_MX_SXO))), id.vars=c('CLU_ORDERED','CLU_MX_SXO'))
	ggplot(dfp, aes(y=CLU_ORDERED)) +
			scale_colour_brewer(palette='Set1') +
			scale_x_continuous(breaks=seq(1950,2030,5), expand=c(0,0)) +
			geom_line(data=tmp, aes(x=value, y=CLU_ORDERED), colour='black') +
			geom_point(aes(x=HIV_YR, colour=SXO)) +
			facet_grid(CLU_MX_SXO~., scales='free_y', space='free_y') +
			labs(x='\nyear of diagnosis', y='phylogenetic clusters\n', colour='sexual orientation') +
			theme_bw() +
			theme(	strip.text.y=element_text(angle=0),
					panel.grid.minor.x=element_line(colour="grey70", size=0.25), 
					panel.grid.major.x=element_line(colour="grey70", size=0.5),
					panel.grid.minor.y=element_blank(), 
					panel.grid.major.y=element_blank(),					
					axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					legend.position='bottom')
	ggsave(gsub('\\.rda','_cluevolution_by_sxo.pdf',infile), w=12,h=18)	
}