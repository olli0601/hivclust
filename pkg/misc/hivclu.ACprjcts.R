######################################################################################
project.ACpolext.rmDRM<- function()
{
	indir	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_SA/ACpolext150831'
	infile	<- 'SATURN150831.csv'
	dfstr	<- data.table(read.csv(file=paste(indir, '/', infile, sep=''), stringsAsFactors=FALSE))
	infile	<- 'Eduan_DRT_170815.csv'
	dfdrt	<- data.table(read.csv(file=paste(indir, '/', infile, sep=''), stringsAsFactors=FALSE))
	infile	<- 'Comet150831.csv'
	dfcmt	<- data.table(read.csv(file=paste(indir, '/', infile, sep=''), stringsAsFactors=FALSE))	
	#	what subtypes do we have? 
	dfcmt[, table(subtype)]
	dfdrt[, table(assignment)]
	#	OK separate all C and the non-C sequences.  
	dfc		<- subset(dfdrt, assignment=='HIV-1 Subtype C', name) 
	dfnc	<- subset(dfdrt, assignment!='HIV-1 Subtype C', name)
	#	clarify: not all C? do we have subtype assignment for all others that are not in SATURN?
	
	
	indir	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_SA/ACpolext150830_original'
	infile	<- 'ZA.SubC.12432.aln.fasta'
	acxs	<- read.dna(file=paste(indir, '/', infile,sep=''), format='fasta')
	#	determine alignment start relative to HXB2
	#	pos 60 is pos 2300 in HXB2 ie pos 1 is pos 2241 in HXB2
	paste(as.character(acxs[ which(grepl("B.FR.K03455.1983",rownames(acxs))), ]), collapse='')
	acxs.st	<- 2241
	#	rm primary DRMs
	load( system.file(package="hivclust", "data",'IAS_primarydrugresistance_201303.rda') )		
	dr				<- as.data.table(IAS_primarydrugresistance_201303)				
	set(dr, NULL, "Alignment.nuc.pos", dr[,Alignment.nuc.pos]-acxs.st+1)		
	seq				<- seq.rm.drugresistance(as.character(acxs), dr, verbose=1, rtn.DNAbin=1 )
	#	removed DR mutations, n= 19340 n per patient= 1.55566280566281
	#	rm any additional DMRs that are indicated in the Excel sheet
	dfdrm			<- subset(dfstr, select=c(Seq_Code, PI_MAJOR, PI_MINOR, NRTI, NNRTI))
	dfdrm			<- subset(melt(dfdrm, id.vars=c('Seq_Code'), variable.name='TYPE', value.name='DRM'), !DRM%in%c('None','No PR seq','No RT seq'))
	dfdrm			<- dfdrm[, list(DRM=unlist(strsplit(DRM, ','))), by=c('Seq_Code','TYPE')]
	dfdrm[, Gene.codon.number:=dfdrm[, regmatches(DRM,regexpr('[0-9]+', DRM))]]
	dfdrm[, Wild.type:=dfdrm[, regmatches(DRM,regexpr('^[^0-9]*', DRM))]]
	dfdrm[, Mutant:=dfdrm[, regmatches(DRM,regexpr('[^0-9]*$', DRM))]]
	dfdrm			<- subset(dfdrm, !grepl('insertion|deletion',Mutant))
	set(dfdrm, NULL, 'Mutant', dfdrm[, gsub('*','',Mutant,fixed=1)])
	dfdrm			<- dfdrm[, list(Mutant= unlist(strsplit(Mutant,''))), by=c('Seq_Code','TYPE','DRM','Gene.codon.number','Wild.type')]	
	nt2aa			<- as.data.table( read.csv( system.file(package="hivclust", "data",'standard_nt_code.csv.gz'), stringsAsFactors=F ) )
	setnames(nt2aa,c("AA","NTs"),c("Mutant","Mutant.NTs"))
	nt2aa			<- subset(nt2aa, select=c(Mutant,Mutant.NTs))
	dfdrm			<- merge(dfdrm, nt2aa, all.x=1, by="Mutant", allow.cartesian=TRUE)
	dfdrm			<- subset(dfdrm, Mutant!=Wild.type)
	set(dfdrm, NULL, "Mutant.NTs", tolower(dfdrm[,Mutant.NTs]))
	#	get pos of wild type in alignment: align HXB2s
	load( system.file(package="hivclust", "data",'refseq_hiv1_hxb2.rda') )
	tmp				<- hxb2[acxs.st:(acxs.st+ncol(acxs)+100-1L), as.list(as.DNAbin(as.character(HXB2.K03455)))]
	names(tmp)	<- 'HXB2'
	tmp				<- c(as.list(acxs[ which(grepl("B.FR.K03455.1983",rownames(acxs))), ]), tmp)
	write.dna(tmp, file=paste(indir, 'ZA.SubC.HXB2s.fasta', sep='/'), format='fasta')
	#	OK HXB2 is aligned in one stretch :-)
	#	HXB2 pos for PI mutations
	set(dfdrm, NULL,'HXB2.pos', dfdrm[,(as.numeric(Gene.codon.number)-15)*3+2295])
	#	HXB2 pos for RT mutations
	tmp			<- dfdrm[, which(grepl('RT', TYPE))]
	set(dfdrm, tmp, 'HXB2.pos', dfdrm[tmp,(as.numeric(Gene.codon.number)-1)*3+2550])
	#	save all DRM locations for future use
	dr			<- unique(subset(dfdrm, select=c(DRM, Gene.codon.number, Wild.type, Mutant.NTs, HXB2.pos)))
	setnames(dr, c('DRM'), c('DR.name'))
	save(dr, file=paste( CODE.HOME,"/data/AC_drugresistance_201508.rda",sep=''))
	#,'Alignment.nuc.pos'
	set(dr, NULL,'Alignment.nuc.pos', dr[, HXB2.pos-acxs.st+1L])
	tmp			<- subset(dr, select=c(DR.name, Mutant.NTs, Alignment.nuc.pos))
	seq			<- seq.rm.drugresistance(as.character(acxs), tmp, verbose=1, rtn.DNAbin=1 )
	outdir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_SA/ACextpol150831'
	outfile		<- 'ZA_SubC_12432_nDRM.fasta'
	write.dna(seq, paste(outdir, '/', outfile, sep=''), format='fasta', colsep='', nbcol=-1)
	save(seq, file=paste(outdir, '/', gsub('\\.fasta','\\.R',outfile), sep=''))	
}
######################################################################################
project.ACpolext.trees.inspect<- function()
{	  	
	indir		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_SA/ACpolext150831/ExaML'
	infiles		<- data.table(FILE=list.files(indir, pattern='^ExaML_result', full.names=1))
	require(phytools)
	i			<- 1
	ph			<- read.tree(infiles[i,FILE])	
	tmp			<- which(ph$tip.label=="B.FR.K03455.1983")
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
	
	ph			<- ladderize(ph)
	phd			<- data.table(TAXA= ph$tip.label)
	phd[, AC:= grepl('AC_',TAXA)]
	phd[, TIP.CLR:= 'black']
	set(phd, phd[, which(AC)], 'TIP.CLR','red')
	
	pdf(file= infiles[i,paste(FILE,'.pdf',sep='')], width=40, height=250)
	plot(ph, tip.color= phd[, TIP.CLR], cex=0.5, adj=1)
	dev.off()
	
	infile		<- "ZA_SubC_12432_nDRM"
	signat.in	<- '150831' 
	signat.out	<- '150831'
	#	ExaML bootstrap args
	bs.from		<- 1
	bs.to		<- 1
	bs.n		<- 500
	outdir		<- indir
	
	#args.parser	<- paste("-m DNA -q",PARTITION)
	args.parser	<- paste("-m DNA")
	args.examl	<- "-f d -m GAMMA"	
	cmd			<- hivc.cmd.examl.bootstrap(indir, infile, signat.in, signat.out, bs.from=bs.from, bs.to=bs.to, opt.bootstrap.by="codon", args.parser=args.parser, args.examl=args.examl, tmpdir.prefix="examl")					
	dummy		<- lapply(cmd, function(x)
			{				
				#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
				x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=34, hpc.q="pqeph", hpc.mem="450mb", hpc.nproc=1)
				signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
				outfile	<- paste("exa",signat,sep='.')
				#cat(x)
				hivc.cmd.hpccaller(outdir, outfile, x)
				Sys.sleep(1)
			})	
}
######################################################################################
project.ACpolext.examl<- function()
{	  	
	indir		<- DATA
	infile		<- "ZA_SubC_12432_nDRM"
	signat.in	<- '150831' 
	signat.out	<- '150831'
	signat.out	<- '150906'
	#	ExaML bootstrap args
	bs.from		<- 0
	bs.to		<- 0
	bs.n		<- 500
	outdir		<- indir
	
	#args.parser	<- paste("-m DNA -q",PARTITION)
	args.parser	<- paste("-m DNA")
	args.examl	<- "-f d -D -m GAMMA"
	args.examl	<- "-f o -D -m GAMMA"
	cmd			<- hivc.cmd.examl.bootstrap(indir, infile, signat.in, signat.out, bs.from=bs.from, bs.to=bs.to, opt.bootstrap.by="codon", args.parser=args.parser, args.examl=args.examl, tmpdir.prefix="examl")					
	dummy		<- lapply(cmd, function(x)
			{				
				#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
				x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=999, hpc.q="pqeph", hpc.mem="800mb", hpc.nproc=1)
				signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
				outfile	<- paste("exa",signat,sep='.')
				cat(x)
				#hivc.cmd.hpccaller(outdir, outfile, x)
				Sys.sleep(1)
			})	
}