#' @export
HIVC.COUNTRY.TABLE<- data.table(matrix(c(	"AG","ANTIGUA AND BARBUDA",		"AL","ALBANIA",		"AN", "ANTILLES", "AO","ANGOLA",	"AR","ARGENTINA",	"AT","AUSTRIA",		"AU","AUSTRALIA",				
						"BB","BARBADOS",	"BE","BELGIUM",	"BG","BULGARIA",	"BR","BRAZIL",	"CA","CANADA",	"CH","SWITZERLAND",		"CI","COTE D'IVOIRE",	"CL","CHILE",	
						"CN","CHINA",	"CO","COLOMBIA",	"CU","CUBA",	"CW", "CURACAO", "CV","CAPE VERDE",	"CY","CYPRUS",	"CZ","CZECH REPUBLIC",	"DE","GERMANY",	"DK","DENMARK",	
						"DO","DOMINICAN REPUBLIC",	"DZ","ALGERIA",		"EC","ECUADOR",	"ER", "ERITREA", "ES","SPAIN",	"FI","FINLAND",	"FR","FRANCE",	"GA","GABON",	"GB","UNITED KINGDOM",	
						"GD","GRENADA",	"GE","GIORGIA",	"GH","GHANA",	"GR","GREECE",	"GW","GUINEA-BISSAU",	"GY","GUYANA",	"HN","HONDURAS", "HT","HAITI",	
						"HU","HUNGARY",	"ID", "INDONESIA", "IL","ISRAEL",	"IN","INDIA",	"IT","ITALY",	"JM","JAMAICA",	"JP","JAPAN",	"KE","KENYA",	"KR","SOUTH KOREA",	"LB","LEBANON",				
						"LU","LUXEMBOURG",	"LK", "SRILANKA", "LV","LATVIA",	"MA","MOROCCO",	"ME","MONTENEGRO",	"MG","MADAGASCAR",	"ML","MALI",	"MN","MONGOLIA",	"MX","MEXICO",	"MY","MALAYSIA",
						"NG","NIGERIA",	"NL","NETHERLANDS",	"NO","NORWAY",	"PA","PANAMA",	"PE","PERU",	"PH","PHILIPPINES",	"PL","POLAND",	"PR","PUERTO RICO",	"PT","PORTUGAL",			
						"PY","PARAGUAY",	"RO","ROMANIA",	"RS","SERBIA",	"RU","RUSSIAN FEDERATION",	"SC","SEYCHELLES",	"SD","SUDAN",	"SE","SWEDEN",	"SG","SINGAPORE",	
						"SI","SLOVENIA",	"SK","SLOVAKIA",	"SN","SENEGAL",	"SR","SURINAME",	"SV","EL SALVADOR",	"TH","THAILAND",	"TT","TRINIDAD AND TOBAGO",	"TW","TAIWAN",			
						"UA","UKRAINE",	"UG","UGANDA",	"US","UNITED STATES",	"UY","URUGUAY",	"VE","VENEZUELA",	"VN","VIETNAM",	"YE","YEMEN","ZA","SOUTH AFRICA"),ncol=2,byrow=1,dimnames=list(c(),c("key","country"))))

######################################################################################
#create PROT+RT data set of first sequences from all patients
hivc.prog.get.bootstrapseq<- function(check.any.bs.identical=0)
{	
	library(ape)
	library(data.table)
	library(hivclust)
	
	indir				<- outdir		<- paste(DATA,"tmp",sep='/')
	infile				<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
	signat.out			<- signat.in	<- "Sat_May_11_14/23/46_2013"
	verbose				<- resume		<- 1
	opt.bootstrap.by	<- "codon"
	bs					<- 0
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.in<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.out<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									bootstrap= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) bs<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,3),
									by= return(substr(arg,5,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.bootstrap.by<- tmp[1]
	}
	if(1)
	{
		print( indir ) 
		print(outdir)
		print(infile)
		print(signat.in)
		print(signat.out)
		print(verbose)
		print(resume)
		print(bs)
		print(opt.bootstrap.by)
		print(signat.in)
		print(signat.out)
	}
	if(!opt.bootstrap.by%in%c("nucleotide","codon"))	stop("Unexpected opt.bootstrap.by")		
	pattern 	<- paste(infile,"_",gsub('/',':',signat.out),".phylip.",sprintf("%03d",bs),sep='')
	file		<- list.files(path=outdir, pattern=pattern, full.names=1)
	if(!resume || !length(file))	
	{					
		file		<- paste(outdir,"/",infile,"_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nread",file))
		tmp			<- load(file)
		if(length(tmp)!=1)	stop("Unexpected lenght of loaded objects")
		eval(parse(text=paste("seq<- ",tmp,sep='')))		
		print(seq)
		#print(bs)
		if(bs)		#keep bs0 intact
		{
			dummy			<- 0
			any.eq			<- 1
			j				<- 0
			while(any.eq)
			{
				j			<- j+1
				if(opt.bootstrap.by=="codon")
				{
					bs.blocks.n	<- floor( ncol(seq )/3)
					bs.blocks.s	<- sample(seq_len(bs.blocks.n),bs.blocks.n,replace=T)-1
					bs.seq.s	<- as.numeric( sapply(bs.blocks.s,function(x)		3*x+c(1,2,3)		) )
				}
				else if(opt.bootstrap.by=="nucleotide")
				{
					bs.seq.s	<- sample(seq_len(ncol(seq )),ncol(seq ),replace=T)
				}
				seq.BS		<- seq[,bs.seq.s]				
				if(check.any.bs.identical)
				{
					if(verbose) cat(paste("\ncheck for identity proposed boostrap seq alignment no",j))
					#check no seqs identical								
					for(i1 in seq_len(nrow(seq.BS)-1))
					{
						seq1		<- seq.BS[i1,]
						tmp			<- 1-sapply(seq.int(i1+1,nrow(seq.BS)),function(i2)
													{		
														.C("hivc_dist_ambiguous_dna", seq1, seq.BS[i2,], ncol(seq1), dummy )[[4]]			
													})
						#print(tmp)
						if(any(tmp==0))
						{
							print(tmp)
							break
						}									
						if(i1==nrow(seq.BS)-1)
							any.eq	<- 0
					}
					if(verbose) cat(paste("\nchecked for identity proposed boostrap seq alignment no",j,"is any identical",any.eq))
				}
				else
					any.eq	<- 0
			}					
		}
		else
		{
			cat(paste("\nkeep boostrap seq alignment no",bs,"as original"))
			seq.BS	<- seq
		}
		file		<- paste(outdir,"/",infile,"_",gsub('/',':',signat.out),".phylip.",sprintf("%03d",bs),sep='')
		cat(paste("\nsave boostrap seq alignment to",file))
		seq.write.dna.phylip(seq.BS, file=file)
	}
	else
		cat("\nfound boostrap sequence alignment")
}
######################################################################################
#create PROT+RT data set of first sequences from all patients
hivc.prog.get.allseq<- function()
{	
	library(ape)
	library(data.table)
	library(RFLPtools)
	library(hivclust)
	
	indir		<- paste(DATA,"derived",sep='/')
	outdir		<- paste(DATA,"tmp",sep='/')
	infile		<- "ATHENA_2013_03_All1stPatientCovariates.R"
	signat.out	<- "Sat_Jun_16_17/23/46_2013"
	signat.in	<- "Wed_May__1_17/08/15_2013"
	verbose		<- 1 
	resume		<- 0
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.in<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.out<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
	}
	if(1)
	{
		print(indir)
		print(outdir)
		print(infile)
		print(signat.in)
		print(signat.out)
		print(verbose)
		print(resume)
	}
	pattern 	<- gsub('/',':',paste("ATHENA_2013_03_All+LANL_Sequences_",signat.out,".R$",sep=''))
	file		<- list.files(path=outdir, pattern=pattern, full.names=1)
	if(!resume || !length(file))	
	{						
		#create file of sequences that are in preliminary cluster. use these as seed to enrich data set to get more reliable boostrap
		dir.name		<- DATA
		infile.tree		<- "ATHENA_2013_03_FirstCurSequences_PROTRT_examlbs100"
		infile.seq		<- "ATHENA_2013_03_FirstCurSequences_PROTRT"
		signat.in		<- "Sat_May_11_14/23/46_2013"
		
		
		infile.enrich	<- "LosAlamos_HIV1B_Prot_P51_n4243_gapsyes"
		outfile.enrich	<- "LosAlamos_HIV1B_Prot_P51_n4243_gapsno"
		infile.enrich	<- "LosAlamos_HIV1B_Prot_n85928_gapsyes"
		outfile.enrich	<- "LosAlamos_HIV1B_Prot_n85928_gapsno"

		if(0) 
		{
			insignat				<- "Sat_Jun_15_18/23/46_2013"
			#build BLAST database
			#remove all gaps from FASTA file
			file<- paste(paste(dir.name,"original",sep='/'),'/',paste(infile.enrich,'_',gsub('/',':',insignat),".fasta", sep=''), sep='')						
			seq.enrich				<- read.dna( file, format="fa" )
			seq.enrich				<- seq.rmgaps(seq.enrich, rm.only.col.gaps=0)
			names(seq.enrich)		<- gsub('-','NA',names(seq.enrich))					
			file<- paste(paste(dir.name,"derived",sep='/'),'/',paste(outfile.enrich,'_',gsub('/',':',insignat),".fasta", sep=''), sep='')			
			write.dna(seq.enrich, file=file, format="fasta", append=0, colsep='', colw=80, blocksep=0)
			#build BLAST database
			cmd						<- hivc.cmd.blast.makedb(paste(dir.name,"derived",sep='/'), outfile.enrich, signat=gsub('/',':',insignat), with.mask=0, verbose=0)
			cat ( cmd )
			#cat( system(cmd, intern=TRUE) )
		}
		if(0)
		{	
			plot							<- 0
			#identify query sequences for BLAST search from large trial clusters
			signat.in						<- "Fri_May_24_12/59/06_2013"						
			#read tree to construct large trial clusters
			file							<- paste(dir.name,"tmp",paste(infile.tree,'_',gsub('/',':',signat.in),".newick",sep=''),sep='/')
			cat(paste("read file",file))
			ph								<- ladderize( read.tree(file) )
			#read bootstrap support values		
			ph.node.bs						<- as.numeric( ph$node.label )
			ph.node.bs[is.na(ph.node.bs)]	<- 0
			ph.node.bs						<- ph.node.bs/100
			ph$node.label					<- ph.node.bs
			dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf")			
			#produce trial clustering
			thresh.bs						<- 0.9
			thresh.brl						<- 0.105		#subst rate 3.4 10^-3  so for 15 years either way have 10.5 * 10^-2
			clustering						<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
			#get tip names in trial cluster
			ph.tips.in.cluster				<- which( !is.na(clustering[["clu.mem"]][seq_len(Ntip(ph))]) )
			ph.tips.in.cluster.names		<- ph$tip.label[ph.tips.in.cluster]
			#get sequences in trial cluster			
			file							<- paste(dir.name,"tmp",paste(infile.seq,'_',gsub('/',':',signat.in),".R",sep=''),sep='/')
			load(file)
			if(verbose) cat(paste("\nload trial sequences from",file,sep=''))
			seq.Trial.PROT.RT				<- seq.PROT.RT[ph.tips.in.cluster.names,]
			seq.Trial.PROT.RT				<- seq.replace(seq.Trial.PROT.RT, code.from='?', code.to='n')			
			if(verbose) cat(paste("\nfound trial sequences, n=",nrow(seq.Trial.PROT.RT),sep=''))			
			outsignat						<- "Sat_Jun_15_18/23/46_2013"
			file							<- paste(paste(dir.name,"tmp",sep='/'),'/',paste(infile.tree,"_intrialcluster_",gsub('/',':',outsignat),".fasta", sep=''), sep='')
			if(verbose) cat(paste("\nwrite trial sequences to",file,sep=''))
			write.dna(seq.Trial.PROT.RT, file=file, format="fasta", append=0, colsep='', colw=ncol(seq.Trial.PROT.RT), blocksep=0)
			
			#BLAST search against PROT+P51 database
			indir							<- paste(dir.name,"tmp",sep='/')
			infile							<- paste(infile.tree,"_intrialcluster",sep='')
			insignat						<- "Sat_Jun_15_18/23/46_2013"
			dbdir							<- paste(dir.name,"derived",sep='/')
			dbfile							<- "LosAlamos_HIV1B_Prot_P51_n4243_gapsno"			
			outsignat						<- "Sat_Jun_15_18/23/46_2013"
			cmd<- hivc.cmd.blast(indir, infile, gsub('/',':',insignat), dbdir, dbfile, gsub('/',':',insignat), outdir=indir, outsignat=gsub('/',':',outsignat), blast.max_target_seqs=20)			
			cat(cmd)
			system(cmd, intern=TRUE)	#wait until finished	
			#read BLAST file against PROT+P51 database and take unique best 10 hits
			file							<- paste(paste(dir.name,"tmp",sep='/'),'/',paste(infile.tree,"_intrialcluster_",gsub('/',':',signat.out),".blast", sep=''), sep='')
			seq.enrich						<- seq.blast.read(file=file)
			if(plot)
			{
				tmp							<- 1-seq.enrich[,identity]/100
				summary(tmp)
				hist(tmp)
			}
			seq.enrich.unique				<- data.table(subject.id=unique(seq.enrich[,subject.id]), key="subject.id")			
			tmp								<- seq.enrich.unique[, substr(subject.id,3,4)]
			tmp[tmp=="NA"]					<- NA
			seq.enrich.unique[,subject.location:=tmp]
			tmp								<- sapply(seq.enrich.unique[, strsplit(subject.id,'.',fixed=1) ], function(x) x[3])
			tmp[tmp=="NA"]					<- NA
			seq.enrich.unique[,subject.PosSeq:=as.numeric(tmp)]
			if(plot)
			{
				hist(seq.enrich.unique[,subject.PosSeq])
				seq.enrich.unique.loc		<- sort(table(seq.enrich.unique[,subject.location]),decreasing=TRUE)
				barplot(seq.enrich.unique.loc, cex.names=.75)
			}			
			#best PROT+p51 hits to enrich NL data set
			indir				<- paste(dir.name,"derived",sep='/')
			infile				<- "LosAlamos_HIV1B_Prot_P51_n4243_gapsno"
			insignat			<- "Sat_Jun_15_18/23/46_2013"
			file				<- paste(indir,'/',paste(infile,'_',gsub('/',':',insignat),".fasta", sep=''), sep='')
			seq.enrich.PROT.P51	<- read.dna( file, format="fa" )
			tmp					<- which( names(seq.enrich.PROT.P51)%in%seq.enrich.unique[,subject.id] )
			tmp2				<- lapply( tmp, function(i)	seq.enrich.PROT.P51[[i]] )
			names(tmp2)			<- paste("PROT+P51",names(seq.enrich.PROT.P51)[tmp],sep="_")
			class(tmp2)			<- "DNAbin"
			seq.enrich.PROT.P51	<- tmp2
			
			
			#BLAST search against PROT database
			indir							<- paste(dir.name,"tmp",sep='/')
			infile							<- paste(infile.tree,"_intrialcluster",sep='')			
			dbdir							<- paste(dir.name,"derived",sep='/')			
			dbfile							<- "LosAlamos_HIV1B_Prot_n85928_gapsno"
			outsignat						<- "Sat_Jun_16_17/23/46_2013"
			cmd<- hivc.cmd.blast(indir, infile, gsub('/',':',insignat), dbdir, dbfile, gsub('/',':',insignat), outdir=indir, outsignat=gsub('/',':',outsignat), blast.max_target_seqs=20)			
			cat(cmd)
			#cat( system(cmd, intern=TRUE) )
			
			#read BLAST file against PROT database and take unique best 10 hits that are geographically distant as TRUE NEGATIVES
			file							<- paste(paste(dir.name,"tmp",sep='/'),'/',paste(infile.tree,"_intrialcluster_",gsub('/',':',outsignat),".blast", sep=''), sep='')
			seq.enrich						<- seq.blast.read(file=file)
			if(plot)
			{
				tmp							<- 1-seq.enrich[,identity]/100
				summary(tmp)
				hist(tmp)
			}
			seq.enrich.unique				<- data.table(subject.id=unique(seq.enrich[,subject.id]), key="subject.id")			
			tmp								<- seq.enrich.unique[, substr(subject.id,3,4)]
			tmp[tmp=="NA"]					<- NA
			seq.enrich.unique[,subject.location:=tmp]
			tmp								<- sapply(seq.enrich.unique[, strsplit(subject.id,'.',fixed=1) ], function(x) x[3])
			tmp[tmp=="NA"]					<- NA
			seq.enrich.unique[,subject.PosSeq:=as.numeric(tmp)]
			if(plot)
			{
				hist(seq.enrich.unique[,subject.PosSeq])
				seq.enrich.unique.loc			<- sort(table(seq.enrich.unique[,subject.location]),decreasing=TRUE)
				barplot(seq.enrich.unique.loc, cex.names=.75)
			}
			
			#to build HIVC.COUNTRY.TABLE, use this to search LANL
			#tmp								<- seq.enrich.unique[,min(subject.id),by=subject.location]
			#tmp								<- sapply(tmp[, strsplit(V1,'.',fixed=1) ], function(x) x[5])			
			#extract seq.enrich.unique that are very unlikely to transmit to NL		seq.enrich.unique.loc.excl taken from Bezemer2013
			seq.enrich.unique.loc.level		<- "0.01"			
			seq.enrich.unique.loc.excl		<- structure(list(	`0.01` = c("GREENLAND", "CHILE", "ALBANIA", "Suriname", "SINGAPORE", "PANAMA", "SENEGAL", "FINLAND", "TAIWAN", "JAPAN", 
																			"SLOVAKIA", "LUXEMBOURG", "AUSTRALIA", "THAILAND", "SERBIA", "ROMANIA", "NORWAY", "VENEZUELA", "SLOVENIA", "RUSSIAN FEDERATION", 
																			"AUSTRIA", "SWEDEN", "FRANCE", "SOUTH KOREA", "CYPRUS", "MONTENEGRO", "CUBA", "HONDURAS", "DENMARK", "POLAND", "PORTUGAL", "CHINA", 
																			"SWITZERLAND", "BRAZIL", "ARGENTINA", "GERMANY", "BELGIUM", "CANADA","CZECH REPUBLIC", "UNITED STATES", "SPAIN", "ITALY", "UNITED KINGDOM","NETHERLANDS"), 
																`0.05` = c("SERBIA", "ROMANIA", "NORWAY", "VENEZUELA", "SLOVENIA","RUSSIAN FEDERATION", "AUSTRIA", "SWEDEN", "FRANCE", "SOUTH KOREA", 
																			"CYPRUS", "MONTENEGRO", "CUBA", "HONDURAS", "DENMARK", "POLAND","PORTUGAL", "CHINA", "SWITZERLAND", "BRAZIL", "ARGENTINA", "GERMANY",
																			"BELGIUM", "CANADA", "CZECH REPUBLIC", "UNITED STATES", "SPAIN","ITALY", "UNITED KINGDOM","NETHERLANDS")), .Names = c("0.01", "0.05"))																										
			seq.enrich.unique.loc.incl		<- subset(HIVC.COUNTRY.TABLE, !country%in%seq.enrich.unique.loc.excl[[seq.enrich.unique.loc.level]] )
			seq.enrich.unique.loc.incl		<- subset(seq.enrich.unique.loc.incl, !country%in%c("BULGARIA","GRENADA","GIORGIA","GREECE","HUNGARY","UKRAINE","LATVIA") )						
			seq.enrich.unlikelytransmission	<- subset(seq.enrich.unique,subject.location%in%as.character(seq.enrich.unique.loc.incl[,key]) )
			if(plot)
			{				
				hist(seq.enrich.unlikelytransmission[,subject.PosSeq])
				seq.enrich.unique.loc			<- sort(table(seq.enrich.unlikelytransmission[,subject.location]),decreasing=TRUE)
				barplot(seq.enrich.unique.loc, cex.names=.75)
			}
			#best geographically distant PROT hits to enrich NL data set
			indir				<- paste(dir.name,"derived",sep='/')
			infile				<- "LosAlamos_HIV1B_Prot_n85928_gapsno"
			insignat			<- "Sat_Jun_15_18/23/46_2013"
			file				<- paste(indir,'/',paste(infile,'_',gsub('/',':',insignat),".fasta", sep=''), sep='')
			seq.enrich.TN		<- read.dna( file, format="fa" )
			tmp					<- which( names(seq.enrich.TN)%in%seq.enrich.unlikelytransmission[,subject.id] )
			tmp2				<- lapply( tmp, function(i)	seq.enrich.TN[[i]] )
			names(tmp2)			<- paste("TN",names(seq.enrich.TN)[tmp],sep="_")
			class(tmp2)			<- "DNAbin"
			seq.enrich.TN		<- tmp2

			#combine all sequence data sets and save
			seq.enrich			<- c(seq.enrich.PROT.P51,seq.enrich.TN)
			outdir				<- paste(dir.name,"derived",sep='/')
			outsignat			<- "Sat_Jun_16_17/23/46_2013"			
			file				<- paste(outdir,"/LosAlamos_EnrichSequences_For_ATHENA201303_",gsub('/',':',outsignat),".R",sep='')
			if(verbose) cat(paste("\nwrite to",file))
			save(seq.enrich, file=file)						
		}
		if(0)
		{
			#get all ATHENA sequences into PROT+RT format, and add foreign sequences
			indir		<- paste(dir.name,"tmp",sep='/')
			outdir		<- paste(dir.name,"tmp",sep='/')
			insignat	<- "Wed_May__1_17/08/15_2013"
			outsignat	<- "Sat_Jun_16_17/23/46_2013"
			outfile		<- "ATHENA_2013_03_All+LANL_Sequences"
			
			if(verbose)	cat(paste("\ncreate LosAlamos_EnrichSequences_For_ATHENA201303 file"))						
			#create matrix of all PROT+RT sequences				
			pattern 	<- gsub('/',':',paste(insignat,".clustalo$",sep=''))
			files		<- list.files(path=indir, pattern=pattern, full.names=1)
			#read all sequences and add a missing one with name "NA" if not in union of all possible sequences take at a sampling date		
			seq.PROT	<- read.dna( files[ grep("PROT",files) ], format="fa", as.matrix=1 )
			if(verbose) cat(paste("\nfound PROT sequences, n=",nrow(seq.PROT),"\n"))
			seq.RT		<- read.dna( files[ grep("RT",files) ], format="fa", as.matrix=1 )
			if(verbose) cat(paste("\nfound RT sequences, n=",nrow(seq.RT),"\n"))
			seq.nam.all	<- union( rownames(seq.PROT), rownames(seq.RT) )
			if(verbose) cat(paste("\nfound unique sampling IDs, n=",length(seq.nam.all),"\n"))
			seq.nam.PROT<- seq.nam.all
			seq.nam.RT	<- seq.nam.all
			seq.nam.PROT[ !seq.nam.all%in%rownames(seq.PROT) ]	<- "NA"		#those sampling dates missing among PROT get name NA
			seq.nam.RT[ !seq.nam.all%in%rownames(seq.RT) ]		<- "NA"
			#prepare PROT and RT sequences: add a missing one with name "NA" 		
			seq.PROT	<- read.dna( files[ grep("PROT",files) ], format="fa", as.matrix=1 )
			tmp			<- as.DNAbin( matrix(rep('-',ncol(seq.PROT)),1,ncol(seq.PROT), dimnames=list(c("NA"),c())) )
			seq.PROT	<- rbind(seq.PROT,tmp)
			seq.RT		<- read.dna( files[ grep("RT",files) ], format="fa", as.matrix=1 )						 				
			tmp			<- as.DNAbin( matrix(rep('-',ncol(seq.RT)),1,ncol(seq.RT), dimnames=list(c("NA"),c())) )
			seq.RT		<- rbind(seq.RT,tmp)
			#add missing sequence where needed to combine PROT and RT 
			seq.PROT	<- seq.PROT[seq.nam.PROT,]
			seq.RT		<- seq.RT[seq.nam.RT,]
			rownames(seq.PROT)	<- seq.nam.all
			rownames(seq.RT)	<- seq.nam.all
			#combine PROT and RT
			seq.PROT.RT	<- cbind(seq.PROT,seq.RT)
			if(verbose) cat(paste("\ncombined PROT and RT sequences, n=",nrow(seq.PROT.RT),"\n"))
			print(seq.PROT.RT)
			
			#load EnrichSequences
			indir				<- paste(dir.name,"derived",sep='/')
			insignat			<- "Sat_Jun_16_17/23/46_2013"
			file				<- paste(indir,"/LosAlamos_EnrichSequences_For_ATHENA201303_",gsub('/',':',insignat),".R",sep='')
			if(verbose) cat(paste("\nloading file",file))
			load(file)	
			print(seq.enrich)
			#load HXB2 reference sequence
			data( refseq_hiv1_hxb2 )
			hxb2				<- as.character( data.table( hxb2 )[, HXB2.K03455 ] )
			hxb2				<- hxb2[seq_len(length(hxb2)-2)]
			
			#write all sequences to fasta file
			file		<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".fasta",sep='')			
			if(verbose) cat(paste("\nwrite reference sequence to",file,"\n"))
			cat( paste(">HXB2\n",paste( hxb2, collapse='',sep='' ),"\n",sep=''), file=file, append=0)			
			if(verbose) cat(paste("\nappend ATHENA combined PROT and RT sequences to",file,"\n"))
			write.dna(seq.PROT.RT, file=file, format="fasta", append=1, colsep='', colw=length(hxb2), blocksep=0)
			if(verbose) cat(paste("\nappend LosAlamos_EnrichSequences to",file,"\n"))
			write.dna(seq.enrich, file=file, format="fasta", append=1, colsep='', colw=length(hxb2), blocksep=0)				 
		}		
		if(1)	#curate alignment with reference 
		{
			indir								<- paste(dir.name,"tmp",sep='/')
			infile								<- "ATHENA_2013_03_All+LANL_Sequences"
			outdir								<- paste(dir.name,"tmp",sep='/')
			outfile								<- "ATHENA_2013_03_CurAll+LANL_Sequences"
			insignat							<- "Sat_Jun_16_17/23/46_2013"
			outsignat							<- "Sat_Jun_16_17/23/46_2013"
			file								<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo",sep='')			
			if(verbose) cat(paste("\nread ",file))
			seq.PROT.RT							<- read.dna(file, format="fasta", as.character=1)
			query.yes							<- seq.find(seq.PROT.RT, 2273, c("g","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2273:2275]<- matrix( c("-","-","g"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 2348, c("?","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2348:2349]<- matrix( c("g","?"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- seq.find(seq.PROT.RT, 2350, c("-","-","-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2350:2353]<- matrix( c("t","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			#always delete drug resistance from pos 2356 because this is really cryptic; leave g at 2355 always OK
			seq.PROT.RT[,2356:2363]				<- matrix( c("-","-","-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=8, byrow=1 )			
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo1",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)			
			#manual edits seq 3151 3243 9049
			query.yes							<- seq.find(seq.PROT.RT, 2624, c("t","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2624:2625]<- matrix( c("-","t"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- seq.find(seq.PROT.RT, 2634, c("-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2634:2635]<- matrix( c("g","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			seq.PROT.RT.sort.by					<- apply(seq.PROT.RT,1,function(x)		which(rev(x)!="-" )[1]  )
			seq.PROT.RT							<- seq.PROT.RT[sort(seq.PROT.RT.sort.by, index.return=1)$ix,]
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo2",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			
			query.yes							<- seq.find(seq.PROT.RT, 2759, c("-","g","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("g","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 2759, c("-","a","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("a","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 2760, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- seq.find(seq.PROT.RT, 2760, c("-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("g","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- seq.find(seq.PROT.RT, 2760, c("c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("-","c"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- seq.find(seq.PROT.RT, 2761, c("-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2761:2762]<- matrix( c("c","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			#always delete drug resistance from pos 2752 because this is really cryptic; leave a at 2751 always OK
			seq.PROT.RT[,2752:2759]				<- matrix( c("-","-","-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=8, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo3",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			#double fixup needed
			query.yes							<- seq.find(seq.PROT.RT, 2759, c("-","g","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("g","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 2759, c("-","g","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("g","g","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 2759, c("-","a","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("a","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 2759, c("-","r","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2759:2761]<- matrix( c("r","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			seq.PROT.RT[,2752:2760]				<- matrix( c("-","-","-","-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=9, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 2760, c("-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("g","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- seq.find(seq.PROT.RT, 2760, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,2760:2761]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			seq.PROT.RT[,2752:2760]				<- matrix( c("-","-","-","-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=9, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo4",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- seq.find(seq.PROT.RT, 3149, c("-"))
			if(length(query.yes))
			{	
				seq.PROT.RT[query.yes,3149:3156]<- seq.PROT.RT[query.yes,3150:3157]
				seq.PROT.RT[query.yes,3157]		<- "-"
			}	
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo5",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)				
			query.yes							<- seq.find(seq.PROT.RT, 3156, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3156:3157]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- seq.find(seq.PROT.RT, 3157, c("a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3157:3158]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
			#rm col 3213 - not in HB2 and we keep the standard numbering
			seq.PROT.RT							<- seq.PROT.RT[,-3213]
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo6",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- seq.find(seq.PROT.RT, 3243, c("-"))
			if(length(query.yes))
			{	
				seq.PROT.RT[query.yes,3241:3243]<- seq.PROT.RT[query.yes,3240:3242]
				seq.PROT.RT[query.yes,3240]		<- "-"
			}	
			#fixup: rm col 2671 to 2676 - not in HB2 and we return to the standard numbering
			seq.PROT.RT							<- seq.PROT.RT[,c(1:2670,2677:ncol(seq.PROT.RT))]
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo7",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)						
			query.yes							<- seq.find(seq.PROT.RT, 3433, c("-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3433:3434]<- matrix( c("t","-"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- seq.find(seq.PROT.RT, 3433, c("-","-","-","-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3433:3437]<- matrix( c("t","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3434, c("-","-","-","-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3434:3438]<- matrix( c("a","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3434, c("g","-","-","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3434:3438]<- matrix( c("-","-","-","-","g"), nrow=length(query.yes), ncol=5, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo8",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- seq.find(seq.PROT.RT, 3434, c("-","-","-","-","g","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3434:3439]<- matrix( c("g","-","-","-","-","g"), nrow=length(query.yes), ncol=6, byrow=1 )						
			query.yes							<- seq.find(seq.PROT.RT, 3435, c("-","-","-","-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3435:3439]<- matrix( c("a","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3435, c("-","-","-","-","r"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3435:3439]<- matrix( c("r","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )			
			query.yes							<- seq.find(seq.PROT.RT, 3436, c("-","-","-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3439]<- matrix( c("c","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3436, c("-","-","-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3439]<- matrix( c("t","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3436, c("-","-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3438]<- matrix( c("c","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3436, c("-","-","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3438]<- matrix( c("t","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3436, c("-","-","y"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3436:3438]<- matrix( c("y","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )			
			query.yes							<- seq.find(seq.PROT.RT, 3437, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3437:3438]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3438, c("-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3438:3439]<- matrix( c("g","-"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3438, c("a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3437:3438]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
			
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo9",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- seq.find(seq.PROT.RT, 3491, c("c","a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3491:3493]<- matrix( c("-","c","a"), nrow=length(query.yes), ncol=3, byrow=1 )															
			query.yes							<- seq.find(seq.PROT.RT, 3501, c("-","t","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3501:3503]<- matrix( c("t","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3501, c("a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3501:3502]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )			
			query.yes							<- seq.find(seq.PROT.RT, 3503, c("-"))
			for(i in query.yes)
				seq.PROT.RT[i,3501:3503]		<- c("-",seq.PROT.RT[i,3501:3502])			
			query.yes							<- seq.find(seq.PROT.RT, 3504, c("-","-","-","t"))	
			for(i in query.yes)
				seq.PROT.RT[i,3504:3524]		<- c(seq.PROT.RT[i,3507:3524],"-","-","-")
			query.yes							<- seq.find(seq.PROT.RT, 3506, c("-","-","-","t","g"))
			for(i in query.yes)
				seq.PROT.RT[i,3506:3526]		<- c(seq.PROT.RT[i,3509:3526],"-","-","-")
			query.yes							<- seq.find(seq.PROT.RT, 3507, c("-","-","-","g","a"))
			for(i in query.yes)
				seq.PROT.RT[i,3507:3537]		<-	c(seq.PROT.RT[i,3510:3537],"-","-","-")			
			query.yes							<- seq.find(seq.PROT.RT, 3543, c("-","c"))
			if(length(query.yes))
			{	
				seq.PROT.RT[query.yes,3543:3550]<- seq.PROT.RT[query.yes,3544:3551]
				seq.PROT.RT[query.yes,3551]		<- "-"
			}
			query.yes							<- seq.find(seq.PROT.RT, 3549, c("-","c","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3549, c("-","g","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("g","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3549, c("-","g","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("g","t","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3551, c("g","g","c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3551:3554]<- matrix( c("-","g","g","c"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3495, c("g","-","-","-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3495:3499]<- matrix( c("g","g","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo10",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			
			query.yes							<- seq.find(seq.PROT.RT, 3497, c("-","-","-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3497:3500]<- matrix( c("a","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )						
			query.yes							<- seq.find(seq.PROT.RT, 3497, c("-","-","-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3497:3500]<- matrix( c("g","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )									
			query.yes							<- seq.find(seq.PROT.RT, 3499, c("t","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3499:3501]<- matrix( c("-","-","t"), nrow=length(query.yes), ncol=3, byrow=1 )			
			query.yes							<- seq.find(seq.PROT.RT, 3578, c("g","a","g","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3578:3581]<- matrix( c("-","g","a","g"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3578, c("g","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3578:3579]<- matrix( c("-","g"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3551, c("g","g"))
			for(i in query.yes)
				seq.PROT.RT[i,3551:3600]		<-	c("-",seq.PROT.RT[i,3551:3599])
			query.yes							<- seq.find(seq.PROT.RT, 3552, c("g","c"))
			for(i in query.yes)
				seq.PROT.RT[i,3552:3600]		<-	c("-",seq.PROT.RT[i,3552:3599])
			query.yes							<- seq.find(seq.PROT.RT, 3498, c("c","g","t"))
			for(i in query.yes)
				seq.PROT.RT[i,3498:3520]		<-	c("-",seq.PROT.RT[i,3498:3519])
			query.yes							<- seq.find(seq.PROT.RT, 3549, c("a","a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3551, c("g","a","c"))
			for(i in query.yes)
				seq.PROT.RT[i,3551:3560]		<-	c("-",seq.PROT.RT[i,3551:3559])
			query.yes							<- seq.find(seq.PROT.RT, 3551, c("g","a","y","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3551:3554]<- matrix( c("-","g","a","y"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3549, c("-","c","c","g","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","-","c"), nrow=length(query.yes), ncol=3, byrow=1 )
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo11",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			
			#manual curation on the end
			query.yes							<- seq.find(seq.PROT.RT, 3422, c("t","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3422:3424]<- matrix( c("-","-","t"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3494, c("g","g","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3494:3496]<- matrix( c("-","g","g"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3498, c("-","-","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3498:3500]<- matrix( c("g","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3578, c("g","a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3578:3580]<- matrix( c("-","g","a"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3618, c("g","-","-","-","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3618:3622]<- matrix( c("-","-","-","-","g"), nrow=length(query.yes), ncol=5, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3499, c("-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3499:3500]<- matrix( c("c","-"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3504, c("a","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3504:3505]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3437, c("-","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3437:3438]<- matrix( c("a","-"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3438, c("c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3438:3439]<- matrix( c("-","c"), nrow=length(query.yes), ncol=2, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 3444, c("c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,3444:3445]<- matrix( c("-","c"), nrow=length(query.yes), ncol=2, byrow=1 )
			#cut before PROT pol start and at max len in database
			seq.PROT.RT							<- seq.PROT.RT[,2253:ncol(seq.PROT.RT)] 
			seq.PROT.RT							<- seq.PROT.RT[,1:1624]
			seq.PROT.RT.sort.by					<- apply(seq.PROT.RT,1,function(x)		which(rev(x)!="-" )[1]  )
			seq.PROT.RT							<- seq.PROT.RT[sort(seq.PROT.RT.sort.by, index.return=1)$ix,]
			#always delete drug resistance at G333D/E because this is at the end of the alignment and alignment unreliable - could have picked up other ends
			seq.PROT.RT[,1294:1299]				<- matrix( c("-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=6, byrow=1 )			
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo12",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
			query.yes							<- seq.find(seq.PROT.RT, 1300, c("g","g","g"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1302]<- matrix( c("g","g","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 1300, c("g","g","t"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1302]<- matrix( c("g","g","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 1300, c("g","g","a"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1302]<- matrix( c("g","g","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 1300, c("g","a","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1302]<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 1300, c("g","g","-","c"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("g","g","g","c"), nrow=length(query.yes), ncol=4, byrow=1 )			
			write.dna(seq.PROT.RT, file=paste(indir,'/',infile,'_',gsub('/',':',insignat),".fasta.clustalo14",sep=''), format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)						
			query.yes							<- seq.find(seq.PROT.RT, 1300, c("a","a","c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("-","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 1300, c("g","a","t","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("g","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 1300, c("t","g","t","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("t","g","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 1300, c("a","t","c","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("-","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 1300, c("-","g","t","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("-","g","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			query.yes							<- seq.find(seq.PROT.RT, 1300, c("g","a","y","-"))
			if(length(query.yes))
				seq.PROT.RT[query.yes,1300:1303]<- matrix( c("g","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
			
			seq.PROT.RT							<- as.DNAbin( seq.PROT.RT )
			seq.PROT.RT							<- seq.replace(seq.PROT.RT, code.from='?', code.to='n')
			seq.PROT.RT							<- seq.PROT.RT[-1,]												#remove HXB2
			rownames(seq.PROT.RT)				<- gsub(' ','',rownames(seq.PROT.RT))
			
			file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".phylip",sep='')
			if(verbose)	cat(paste("\nwrite phylip file to",file))
			seq.write.dna.phylip(seq.PROT.RT, file=file)			
			file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
			if(verbose)	cat(paste("\nwrite R file to",file))
			save(seq.PROT.RT, file=file)
			file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".fasta",sep='')			
			if(verbose)	cat(paste("\nwrite fasta file to",file))
			write.dna(seq.PROT.RT, file=file, format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)						
		}
		quit("no")				
	}
	else
	{
		if(verbose)	cat(paste("\nload",file))
		load(file)
	}	
}
######################################################################################
#create PROT+RT data set of first sequences from all patients
hivc.prog.get.firstseq<- function()
{	
	library(ape)
	library(data.table)
	library(hivclust)
	
	indir		<- outdir		<- paste(DATA,"tmp",sep='/')
	infile		<- "ATHENA_2013_03_SeqMaster.R"
	signat.out	<- signat.in	<- "Wed_May__1_17/08/15_2013"
	verbose		<- resume		<- 1
	
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									insignat= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.in<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat.out<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
	}
	if(1)
	{
		print(indir)
		print(outdir)
		print(infile)
		print(signat.in)
		print(signat.out)
		print(verbose)
		print(resume)
	}
	pattern 	<- gsub('/',':',paste("FirstAliSequences_PROTRT_",signat.in,".R$",sep=''))
	file		<- list.files(path=outdir, pattern=pattern, full.names=1)
	if(!resume || !length(file))	
	{		
		if(verbose)	cat(paste("\ncreate FirstAliSequences_PROTRT file"))			
		file		<- paste(indir,infile,sep='/')
		if(verbose)	cat(paste("\nload",file,"\n"))
		load(file)
		#str(df.all)
		
		#get correct order of sequence SampleCodes corresponding to first seq of Patient
		seq.PROT.nam<- as.character( df.all[,SeqPROT] )
		seq.PROT.nam[ which(is.na(seq.PROT.nam)) ]	<- "NA"
		seq.RT.nam	<- as.character( df.all[,SeqRT] )
		seq.RT.nam[ which(is.na(seq.RT.nam)) ]		<- "NA"
					
		pattern 	<- gsub('/',':',paste(signat.in,".clustalo$",sep=''))
		files		<- list.files(path=indir, pattern=pattern, full.names=1)
		#read all sequences and add a missing one with name "NA" 		
		seq.PROT	<- read.dna( files[ grep("PROT",files) ], format="fa", as.matrix=1 )
		tmp			<- as.DNAbin( matrix(rep('-',ncol(seq.PROT)),1,ncol(seq.PROT), dimnames=list(c("NA"),c())) )
		seq.PROT	<- rbind(seq.PROT,tmp)
		seq.RT		<- read.dna( files[ grep("RT",files) ], format="fa", as.matrix=1 )						 				
		tmp			<- as.DNAbin( matrix(rep('-',ncol(seq.RT)),1,ncol(seq.RT), dimnames=list(c("NA"),c())) )
		seq.RT		<- rbind(seq.RT,tmp)
		
		seq.PROT	<- seq.PROT[seq.PROT.nam,]
		seq.RT		<- seq.RT[seq.RT.nam,]
		rownames(seq.PROT)<- as.character( df.all[,Patient] )
		rownames(seq.RT)<- as.character( df.all[,Patient] )
		
		seq.PROT.RT	<- cbind(seq.PROT,seq.RT)
		if(verbose) print(seq.PROT.RT)
		file		<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_PROTRT_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nwrite to",file))
		save(seq.PROT.RT, file=file)	
	}
	else
	{
		if(verbose)	cat(paste("\nload",file))
		load(file)
	}		
	if(0){	print("DEBUG FirstAliSequences"); seq.PROT.RT<- seq.PROT.RT[1:10,]	}
	if(0)	#create full fasta file with reference sequence and align once more
	{
		file		<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.out),".fasta",sep='')
		if(verbose) cat(paste("\nwrite to ",file))		
		data( refseq_hiv1_hxb2 )
		hxb2		<- data.table( hxb2 )
		hxb2		<- as.character( hxb2[, HXB2.K03455 ] )
		cat( paste(">HXB2\n",paste( hxb2[ seq.int(1,length(hxb2)-2) ], collapse='',sep='' ),"\n",sep=''), file=file, append=1)
		write.dna(seq.PROT.RT, file=file, format="fasta", append=1, colsep='', colw=length(hxb2)-2, blocksep=0)
				
		file<- paste("ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.out),".fasta",sep='')
		cmd<- hivc.cmd.clustalo(outdir, file, signat='')
		system(cmd)		
	}
	if(0)	#curate alignment with reference 
	{
		file								<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.out),".fasta.clustalo",sep='')
		if(verbose) cat(paste("\nread ",file))
		seq.PROT.RT							<- read.dna(file, format="fasta", as.character=1)
		query.yes							<- seq.find(seq.PROT.RT, 2550, c("-","-","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,2550:2552]<- matrix( c("c","c","y"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 2550, c("c","c","y","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,2550:2553]<- matrix( c("-","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 2751, c("-","-","-","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,2751:2754]<- matrix( c("a","-","-","-"), nrow=length(query.yes), ncol=4, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3143, c("-"))
		for(i in query.yes)	
			seq.PROT.RT[i,3143:3151]		<-	c(seq.PROT.RT[i,3144:3151],"-") 	#align pos 3143 and move gap into 3rd codon
		query.yes							<- seq.find(seq.PROT.RT, 3151, c("a","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3151:3152]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
		#always delete drug resistance insertion at pos 69 because this is really cryptic
		seq.PROT.RT[,2752:2754]				<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3237, c("-"))
		for(i in query.yes)	
			seq.PROT.RT[i,3233:3237]<- c("-",seq.PROT.RT[i,3233:3236]) 		#align pos 3237 and move gap into 3rd codon
		query.yes							<- seq.find(seq.PROT.RT, 3433, c("-","t"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3433:3434]<- matrix( c("t","-"), nrow=length(query.yes), ncol=2, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3433, c("-","-","-","-","t"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3433:3437]<- matrix( c("t","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3434, c("-","-","-","-","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3434:3438]<- matrix( c("a","-","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3503, c("-"))
		for(i in query.yes)
			seq.PROT.RT[i,3501:3503]		<- c("-",seq.PROT.RT[i,3501:3502])
		query.yes							<- seq.find(seq.PROT.RT, 3504, c("-","-","-","t"))	
		for(i in query.yes)
			seq.PROT.RT[i,3504:3524]		<- c(seq.PROT.RT[i,3507:3524],"-","-","-")
		query.yes							<- seq.find(seq.PROT.RT, 3506, c("-","-","-","t","g"))
		for(i in query.yes)
			seq.PROT.RT[i,3506:3526]		<- c(seq.PROT.RT[i,3509:3526],"-","-","-")
		query.yes							<- seq.find(seq.PROT.RT, 3507, c("-","-","-","g","a"))
		for(i in query.yes)
			seq.PROT.RT[i,3507:3537]		<-	c(seq.PROT.RT[i,3510:3537],"-","-","-")
		query.yes							<- seq.find(seq.PROT.RT, 3504, c("g","a","c"))
		for(i in query.yes)
			seq.PROT.RT[i,3504:3530]		<-	c("-","-","-",seq.PROT.RT[i,3504:3527])
		query.yes							<- seq.find(seq.PROT.RT, 3495, c("g","-","-","-","g"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3495:3499]<- matrix( c("g","g","-","-","-"), nrow=length(query.yes), ncol=5, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3501, c("-","t","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3501:3503]<- matrix( c("t","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3501, c("a","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3501:3502]<- matrix( c("-","a"), nrow=length(query.yes), ncol=2, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3551, c("g","g"))
		for(i in query.yes)
			seq.PROT.RT[i,3551:3600]		<-	c("-",seq.PROT.RT[i,3551:3599])
		query.yes							<- seq.find(seq.PROT.RT, 3552, c("g","c"))
		for(i in query.yes)
			seq.PROT.RT[i,3552:3600]		<-	c("-",seq.PROT.RT[i,3552:3599])
		query.yes							<- seq.find(seq.PROT.RT, 3501, c("-","t","t"))
		for(i in query.yes)
			seq.PROT.RT[i,3501:3600]		<-	c("-",seq.PROT.RT[i,3501:3599])
		query.yes							<- seq.find(seq.PROT.RT, 3498, c("c","g","t"))
		for(i in query.yes)
			seq.PROT.RT[i,3498:3520]		<-	c("-",seq.PROT.RT[i,3498:3519])
		query.yes							<- seq.find(seq.PROT.RT, 3549, c("-","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","a","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3543, c("-","c","a","g","g","g","g","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","g","g","g","g","c","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3543, c("-","c","a","g","g","a","g","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","g","g","a","g","c","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3543, c("-","c","a","k","g","g","a","c","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","k","g","g","a","c","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3543, c("-","c","a","g","g","g","g","g","a"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3543:3551]<- matrix( c("c","a","g","g","g","g","g","a","-"), nrow=length(query.yes), ncol=9, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3549, c("a","a","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3549:3551]<- matrix( c("-","-","-"), nrow=length(query.yes), ncol=3, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3551, c("g","a","c"))
		for(i in query.yes)
			seq.PROT.RT[i,3551:3560]		<-	c("-",seq.PROT.RT[i,3551:3559])
		query.yes							<- seq.find(seq.PROT.RT, 3551, c("g","a","y","-"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3551:3554]<- matrix( c("-","g","a","y"), nrow=length(query.yes), ncol=4, byrow=1 )
		query.yes							<- seq.find(seq.PROT.RT, 3549, c("-","c","c","g","g"))
		if(length(query.yes))
			seq.PROT.RT[query.yes,3549:3551]<- matrix( c("c","-","c"), nrow=length(query.yes), ncol=3, byrow=1 )
		#some more manual editing at the end
		#cut before PROT pol start and at max len in database
		seq.PROT.RT							<- seq.PROT.RT[,2253:ncol(seq.PROT.RT)] 
		seq.PROT.RT							<- seq.PROT.RT[,1:1624]
		seq.PROT.RT.sort.by					<- apply(seq.PROT.RT,1,function(x)		which(rev(x)!="-" )[1]  )
		seq.PROT.RT							<- seq.PROT.RT[sort(seq.PROT.RT.sort.by, index.return=1)$ix,]
		#some more manual editing at all ends after sorting
		#always delete drug resistance at G333D/E because this is at the end of the alignment and alignment unreliable - could have picked up other ends
		seq.PROT.RT[,1294:1299]				<- matrix( c("-","-","-","-","-","-"), nrow=nrow(seq.PROT.RT), ncol=6, byrow=1 )

		seq.PROT.RT							<- as.DNAbin( seq.PROT.RT )
		file								<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.out),".fasta.clustalo2",sep='')		
		write.dna(seq.PROT.RT, file=file, format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
	}
	if(0)	#create final HXB2PROTRT R and phylip files; need phylip for ExaML
	{
		signat.in	<- "Wed_May__1_17/08/15_2013"
		signat.out	<- "Sat_May_11_14:23:46_2013"
		file<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_HXB2PROTRT_",gsub('/',':',signat.in),".fasta.clustalo4",sep='')
		seq.PROT.RT	<- read.dna(file, format="fasta")
		print(seq.PROT.RT)
						
		seq.PROT.RT	<- seq.PROT.RT[ rownames(seq.PROT.RT)!="HXB2", ]
		
		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRT_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nwrite to",file))
		save(seq.PROT.RT, file=file)
		
		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRT_",gsub('/',':',signat.out),".phylip",sep='')
		if(verbose) cat(paste("\nwrite to ",file))
		seq.write.dna.phylip(seq.PROT.RT, file=file)				
	}
	if(1)	#retain only third codon positions
	{
		signat.in	<- "Sat_May_11_14:23:46_2013"
		signat.out	<- "Sat_May_11_14:23:46_2013"

		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRT_",gsub('/',':',signat.in),".R",sep='')
		if(verbose) cat(paste("\nload ",file))			
		load(file)
		
		codon3.idx	<- seq.int(3,ncol(seq.PROT.RT),3)
		seq.PROT.RT3<- seq.PROT.RT[, codon3.idx]
		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRTCD3_",gsub('/',':',signat.out),".R",sep='')
		if(verbose) cat(paste("\nwrite to ",file))
		save(seq.PROT.RT3, file=file)
		
		file		<- paste(outdir,"/ATHENA_2013_03_FirstCurSequences_PROTRTCD3_",gsub('/',':',signat.out),".phylip",sep='')
		if(verbose) cat(paste("\nwrite to ",file))
		seq.write.dna.phylip(seq.PROT.RT3, file=file)
	}
}
######################################################################################
hivc.prog.get.clustering.MSM<- function(clu.pre= NULL)
{
	verbose		<- 1
	resume		<- 1
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	indircov	<- paste(DATA,"derived",sep='/')
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"
	opt.brl		<- "dist.brl.casc" 
	thresh.brl	<- 0.096
	thresh.bs	<- 0.8
	
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
		if(length(tmp)>0) infilecov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									opt.brl= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.brl<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									thresh.bs= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.bs<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									thresh.brl= return(as.numeric(substr(arg,13,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.brl<- tmp[1]		
	}	
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infilecov)
		print(resume)
		print(opt.brl)
		print(thresh.brl)
		print(thresh.bs)		
	}	
	#stop()
	if(resume)												#//load if there is R Master data.table
	{
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr",sep='')
		outsignat		<- insignat	
		file			<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')		
		options(show.error.messages = FALSE)		
		readAttempt		<-try(suppressWarnings(load(file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))			
		options(show.error.messages = TRUE)		
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{			
		#
		# precompute clustering stuff		
		#
		if(is.null(clu.pre))
		{
			argv		<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
			argv		<<- unlist(strsplit(argv,' '))
			clu.pre		<- hivc.prog.get.clustering.precompute()
		}
		if(0)
		{
			clu.pre$df.seqinfo<- merge(df.all, subset(clu.pre$df.seqinfo,select=c(FASTASampleCode, Node)), all.y=1, by="FASTASampleCode")
			ph<- clu.pre$ph; dist.brl.max<- clu.pre$dist.brl.max; dist.brl.med<- clu.pre$dist.brl.med; dist.brl.casc<- clu.pre$dist.brl.casc; ph.node.bs<- clu.pre$ph.node.bs; ph.linked<- clu.pre$ph.linked; ph.unlinked.info<- clu.pre$ph.unlinked.info; ph.unlinked<- clu.pre$ph.unlinked; df.seqinfo<- clu.pre$df.seqinfo; unlinked.byspace<- clu.pre$unlinked.byspace; unlinked.bytime<- clu.pre$unlinked.bytime; linked.bypatient<- clu.pre$linked.bypatient	
			save(ph, dist.brl.max, dist.brl.med, dist.brl.casc, ph.node.bs, ph.linked, ph.unlinked.info, ph.unlinked, df.seqinfo, unlinked.byspace, unlinked.bytime, linked.bypatient, file="/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tmp/ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100_preclust_Thu_Aug_01_17:05:23_2013.R")
		}
		#
		#precompute clustering stuff for particular thresholds etc	
		argv		<<- hivc.cmd.clustering(indir, infile, insignat, opt.brl, thresh.brl, thresh.bs, resume=resume)				 
		argv		<<- unlist(strsplit(argv,' '))
		clu			<- hivc.prog.get.clustering()
		if(0)
		{
			clu$df.seqinfo<- merge(clu.pre$df.seqinfo,subset(clu$df.seqinfo,select=c(FASTASampleCode, cluster)), all.y=1, by="FASTASampleCode")
			set(clu$df.seqinfo, which( clu$df.seqinfo[,substr(FASTASampleCode,1,2)=="TN"] ), "CountryInfection", "FRGNTN")			
			set(clu$df.seqinfo, which( clu$df.seqinfo[,substr(FASTASampleCode,1,8)=="PROT+P51"] ), "CountryInfection", "FRGN")
			df.seqinfo<- clu$df.seqinfo; clustering<- clu$clustering; clusters.tp<- clu$clusters.tp; clusters.tn<- clu$clusters.tn
			save(df.seqinfo, clustering, clusters.tp, clusters.tn,file="/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tmp/ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100_clust_dist.brl.casc_bs80_brl9.6_Thu_Aug_01_17:05:23_2013.R")		
		}
		#
		# remove singletons
		#
		if(verbose) cat(paste("\nnumber of seq in tree is n=", nrow(clu$df.cluinfo)))
		df.cluinfo	<- subset(clu$df.seqinfo, !is.na(cluster) )
		if(verbose) cat(paste("\nnumber of seq in clusters is n=", nrow(df.cluinfo)))
		if(verbose) cat(paste("\nnumber of clusters is n=", length(unique(df.cluinfo[,cluster]))))
		#
		# remove within patient clusters
		#
		tmp			<- subset(df.cluinfo[,list(clu.is.bwpat=length(unique(Patient))>1),by="cluster"], clu.is.bwpat, cluster )
		df.cluinfo	<- merge(tmp, df.cluinfo, by="cluster", all.x=1)
		if(verbose) cat(paste("\nnumber of seq in clusters between patients is n=", nrow(df.cluinfo)))
		if(verbose) cat(paste("\nnumber of clusters between patients is n=", length(unique(df.cluinfo[,cluster]))))				
		if(0)	#plot clusters that have multiple foreign infections
		{
			outdir			<- indir
			outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"multifrgninfection",sep='')
			outsignat		<- insignat									
			hivc.clu.getplot.multifrgninfection(clu.pre$ph, clu$clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''))
		}					
		if(0)	#plot clusters with mixed exposure group
		{			
			outdir		<- indir
			outfile		<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"mixedexpgr",sep='')
			outsignat	<- insignat		
			hivc.clu.getplot.mixedexposuregroup( clu.pre$ph, clu$clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep='') )
		}	
		#
		# remove clusters with all seq coming from frgn infection
		#
		tmp			<- hivc.clu.getplot.excludeallmultifrgninfection(clu.pre$ph, clu$clustering, df.cluinfo )
		ph			<- tmp$cluphy
		df.cluinfo	<- tmp$cluphy.df
		clustering	<- tmp$cluphy.clustering		
		#
		#get in-country clusters. this splits clusters with a foreign sequence
		#		
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"incountry",sep='')
		outsignat		<- insignat									
		#ph			<- clu.pre$ph; clustering	<- clu$clustering; plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''); char.frgn  	='CountryInfection=="FRGN"'; char.frgntn	='CountryInfection=="FRGNTN"'; 
		incountry		<- hivc.clu.getplot.incountry(ph, clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''))
		ph				<- incountry$cluphy
		clustering		<- incountry$clustering
		df.cluinfo		<- incountry$df.cluinfo		
		#
		# get msm exposure group clusters. this splits clusters with HET-F
		#
		set(df.cluinfo, which( df.cluinfo[,Trm%in%c("BLOOD","BREAST","PREG","NEEACC")] ), "Trm", "OTH" )
		set(df.cluinfo, which( df.cluinfo[,Trm=="HETfa"] ), "Trm", "HET" )		
		set(df.cluinfo, NULL, "Trm", factor(df.cluinfo[,Trm]) )		
		#ph<- incountry$cluphy; 		plot.file	<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''); levels.msm=c("BI","MSM","IDU","NA"); levels.het=c("BI","HET","IDU","NA"); levels.mixed=c("BI","MSM","HET","IDU","NA"); levels.oth="OTH"
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr",sep='')
		outsignat		<- insignat							
		msm				<- hivc.clu.getplot.msmexposuregroup(ph, clustering, df.cluinfo, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''))		
		#
		# collapse within patient subclades in each subtree
		#
		msm				<- hivc.clu.collapse.monophyletic.withinpatientseq(msm$cluphy.subtrees, msm$cluphy.df )
		ph				<- msm$cluphy
		df.cluinfo		<- msm$cluphy.df
		cluphy.subtrees	<- msm$cluphy.subtrees		
		clustering		<- msm$cluphy.clustering
		#
		# get statistics of clusters that describe acute/early infection, branch lengths (#mutations)
		#				
		msm.brl.bwpat				<- hivc.clu.brl.bwpat(cluphy.subtrees, df.cluinfo)
		tmp							<- sapply(msm.brl.bwpat, function(x) c(median(na.omit(x)), sd(na.omit(x))))
		if(any(is.na(tmp[1,])))	stop("unexpected NA in median branch length - all within patient clusters removed ?")
		msm.clusu.df				<- data.table(cluster= as.numeric( names(cluphy.subtrees) ), clu.bwpat.medbrl= tmp[1,])		
		# add number of patients in cluster to 'msm$cluphy.df'
		tmp							<- df.cluinfo[, list(	cluster			= unique(cluster),
															nseq			= length(FASTASampleCode),
															FrgnInfection	= (!is.na(CountryInfection) & CountryInfection!="NL")[1],
															PossAcute		= (isAcute%in%c("Yes","Maybe"))[1],
															AnyPos_T1		= AnyPos_T1[1]
													), by="Patient"]
		tmp							<- tmp[, list(		clu.npat			= length(Patient), 
														clu.ntip			= sum(nseq),
														clu.nFrgnInfection	= length(which(FrgnInfection)),
														clu.fPossAcute		= length(which(PossAcute)) / length(Patient),
														clu.AnyPos_T1		= min(AnyPos_T1)						
												), by="cluster"]
		msm.clusu.df				<- merge(tmp, msm.clusu.df, by="cluster")
		#
		# re-order cluphy by branch length and plot
		#
		df.cluinfo				<- merge( df.cluinfo, msm.clusu.df, by="cluster" )
		setkey(df.cluinfo,clu.bwpat.medbrl)											#sort clusters by median branch length
		cluphy.subtrees			<- lapply( as.character(unique(df.cluinfo[,cluster])), function(name) cluphy.subtrees[[name]] )
		names(cluphy.subtrees)	<- as.character(unique(df.cluinfo[,cluster]))
		#
		# plot
		#		
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr_bybwpatmedbrl",sep='')
		outsignat		<- insignat									
		tmp				<- hivc.clu.polyphyletic.clusters(df.cluinfo, cluphy.subtrees=cluphy.subtrees, plot.file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".pdf",sep=''), pdf.scaley=35, adj.tiplabel= c(-0.05,0.5), cex.tiplabel=0.3, pdf.xlim=0.36)
		cluphy			<- tmp$cluphy
		#
		# save
		#
		outdir			<- indir
		outfile			<- paste(infile,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,'_',"msmexpgr",sep='')
		outsignat		<- insignat	
		file			<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
		if(verbose)	cat(paste("\nsave msm output to",file))
		save(df.cluinfo, cluphy.subtrees, clustering, cluphy, file=file)
	}
	list(df.cluinfo=df.cluinfo, cluphy=cluphy, cluphy.subtrees=cluphy.subtrees, clustering=clustering)
}
######################################################################################
hivc.prog.get.clustering.TPTN<- function(clu.pre= NULL, with.plot=TRUE)
{
	require(RColorBrewer)
	require(colorspace)
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_examlbs100"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	indircov	<- paste(DATA,"derived",sep='/')
	infilecov	<- "ATHENA_2013_03_AllSeqPatientCovariates"							
		
	patient.n	<- 15700
	opt.brl		<- "dist.brl.casc"
	thresh.brl	<- c(seq(0.02,0.05,0.01),seq(0.06,0.12,0.02),seq(0.16,0.24,0.04))
	#thresh.bs	<- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
	thresh.bs	<- c(0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
	resume		<- 1
	verbose		<- 1
	
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
		if(length(tmp)>0) infilecov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									opt.brl= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.brl<- tmp[1]
		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									patient.n= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) patient.n<- tmp[1]
	}	
	outdir			<- indir
	outsignat		<- insignat
	outfile			<- paste(infile,"_clusttptn_",opt.brl,sep='')	
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infilecov)
		print(resume)
		print(opt.brl)
		print(patient.n)
	}
	#precompute clustering stuff	
	if(is.null(clu.pre))
	{
		argv	<<- hivc.cmd.preclustering(indir, infile, insignat, indircov, infilecov, resume=resume)				 
		argv	<<- unlist(strsplit(argv,' '))
		clu.pre	<- hivc.prog.get.clustering.precompute()
	}
	else
		cat(paste("\nuse input clu.pre"))
	dist.brl	<- switch(	opt.brl, 
							"dist.brl.max"		= clu.pre$dist.brl.max,
							"dist.brl.med"		= clu.pre$dist.brl.med,
							"dist.brl.casc"		= clu.pre$dist.brl.casc,
							NA)
	if(any(is.na(dist.brl)))	stop("unexpected NA in dist.brl")
	#
	# evaluate TP and TN for several clustering thresholds	 
	#
	file		<- paste(outdir, '/', outfile, '_', gsub('/',':',outsignat),".R", sep='' )
	if(resume)												#//load if there is R Master data.table
	{
		options(show.error.messages = FALSE)		
		readAttempt<-try(suppressWarnings(load(file)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))			
		options(show.error.messages = TRUE)		
	}
	if(!resume || inherits(readAttempt, "try-error"))		#else generate R Master data.table
	{
		cat(paste("\ncreate file",file))
		clusters		<- lapply(thresh.bs, function(bs)
							{
								clusters	<- lapply(thresh.brl, function(brl)
										{
											hivc.clu.clusterbythresh(clu.pre$ph,thresh.brl=brl,dist.brl=dist.brl,thresh.nodesupport=bs,nodesupport=clu.pre$ph.node.bs,retval="all")
										})
								clusters.tp	<- lapply(clusters,function(x)
										{
											hivc.clu.truepos(x, clu.pre$ph.linked, Ntip(clu.pre$ph), verbose=0)
										})		
								clusters.fp	<- lapply(clusters,function(x)
										{
											hivc.clu.trueneg(x, clu.pre$ph.unlinked.info, clu.pre$ph.unlinked, Ntip(clu.pre$ph), verbose=0)
										})
								list(clu= clusters, tp= clusters.tp, fp= clusters.fp)							
							})	
		names(clusters)	<- thresh.bs	
		#
		# get number of patients in clustering	 
		#
		df.cluinfo						<- copy(clu.pre$df.seqinfo)
		df.cluinfo[,"cluster":= NA]
		setkey(df.cluinfo, Node)
		clusters.nbwpat					<- t( sapply(seq_along(clusters),function(i)
											{				
												sapply(seq_along(clusters[[i]][["clu"]]), function(j)
														{
															clustering			<- clusters[[i]][["clu"]][[j]]
															set(df.cluinfo, NULL, "cluster", clustering[["clu.mem"]][seq_len(Ntip(clu.pre$ph))] )
															df.bwpatclu			<- subset(df.cluinfo, !is.na(cluster) )[ , list(n.pat= length(unique(na.omit(Patient))) ), by="cluster"]
															df.bwpatclu			<- subset(df.bwpatclu, n.pat>1)
															sum(df.bwpatclu[,n.pat]) 							
														})
											}) )
		rownames(clusters.nbwpat)		<- thresh.bs
		colnames(clusters.nbwpat)		<- thresh.brl	
		#
		# get number of sequences in clustering	 
		#	
		clusters.nseq					<- t( sapply(seq_along(clusters),function(i)
											{				
												sapply(seq_along(clusters[[i]][["clu"]]), function(j)
														{
															clustering			<- clusters[[i]][["clu"]][[j]]
															set(df.cluinfo, NULL, "cluster", clustering[["clu.mem"]][seq_len(Ntip(clu.pre$ph))] )
															nrow( subset( df.cluinfo, 	!is.na(cluster) & substr(FASTASampleCode, 1, 2)!="TN" & substr(FASTASampleCode, 1, 8)!="PROT+P51" ) ) 							
														})
											}) )
		rownames(clusters.nseq)		<- thresh.bs
		colnames(clusters.nseq)		<- thresh.brl	
		#
		# get distribution of patients in clustering	 
		#
		clusters.dbwpat				<- lapply(seq_along(clusters),function(i)
										{				
											tmp			<- lapply(seq_along(clusters[[i]][["clu"]]), function(j)
															{
																clustering			<- clusters[[i]][["clu"]][[j]]
																set(df.cluinfo, NULL, "cluster", clustering[["clu.mem"]][seq_len(Ntip(clu.pre$ph))] )
																df.bwpatclu			<- subset(df.cluinfo, !is.na(cluster) )[ , list(n.pat= length(unique(na.omit(Patient))), size=length(FASTASampleCode) ), by="cluster"]
																
																df.bwpatclu			<- df.cluinfo[ , list(n.pat= length(unique(na.omit(Patient))), size=length(FASTASampleCode) ), by="cluster"]													
																t.bwpatclu			<- table( subset(df.bwpatclu, !is.na(cluster) & n.pat>0)[, n.pat] )		#n.pat>0 excludes clusters of foreign seq
																#add non-clustering ATHENA seqs to '1'; these have cluster==NA and Patient!=NA 
																t.bwpatclu[1]		<- t.bwpatclu[1] + subset(df.bwpatclu, is.na(cluster))[,n.pat]
																t.bwpatclu
															})
											names(tmp)	<- thresh.brl
											tmp
										})
		names(clusters.dbwpat)		<- thresh.bs
	
		file			<- paste(outdir, '/', outfile, '_', gsub('/',':',outsignat),".R", sep='')
		save(clusters,clusters.nbwpat,clusters.dbwpat,clusters.nseq,file=file)
	}
	else
		df.cluinfo						<- copy(clu.pre$df.seqinfo)
	
	clusters.cov.epidemic	<- clusters.nbwpat / patient.n
	clusters.cov.patindb	<- clusters.nbwpat / length(unique(subset(df.cluinfo,!is.na(Patient))[,Patient]))	
	clusters.cov.seqindb	<- clusters.nseq / nrow( subset( df.cluinfo, 	substr(FASTASampleCode, 1, 2)!="TN" & substr(FASTASampleCode, 1, 8)!="PROT+P51" ) )
	
	if(with.plot)
	{	
		#
		# plot cluster size distributions per bootstrap
		#
		with.singleton			<- 0
		thresh.brl.select		<- which(thresh.brl %in% c(0.02, 0.04, 0.06, 0.1, 0.16) )
		thresh.bs.select		<- which(thresh.bs %in% c(0.7, 0.8, 0.95) )
		pch						<- 15 + seq_along(thresh.bs.select)
		lapply(thresh.bs.select,function(i)
								{				
									file				<- paste(outdir, paste(infile,"_clustpdf_",opt.brl,"_bs",thresh.bs[i],'_',gsub('/',':',outsignat),".pdf",sep=''),sep='/')				
									pdf(file=file, width=5,height=5)				
									if( with.singleton )
									{
										x					<- lapply(thresh.brl.select, function(j){	as.numeric(names(clusters.dbwpat[[i]][[j]]))	})
										y					<- lapply(thresh.brl.select, function(j){	clusters.dbwpat[[i]][[j]]	})					
									}
									else
									{
										x					<- lapply(thresh.brl.select, function(j){	as.numeric(names(clusters.dbwpat[[i]][[j]]))[-1]	})
										y					<- lapply(thresh.brl.select, function(j){	clusters.dbwpat[[i]][[j]][-1]	})
									}
									xlim				<- range( c(70, unlist(x) ) )
									cols				<- brewer.pal( length(x), "Accent")
									plot(1,1,bty='n',type='n',xlim=xlim, ylim=range(unlist(y)), xlab="patients", ylab="frequency", log='y')
									dummy<- lapply(seq_along(x),function(j)
											{
												points(x[[j]],y[[j]],type='l',col=cols[j],pch=pch[i])
												points(x[[j]],y[[j]],col=cols[j],pch=pch[i], cex=0.5)
											})				
									legend("topright", 		fill=cols, 	legend= paste("BRL", thresh.brl[ thresh.brl.select ]), bty='n', border=NA)
									legend("bottomright", 	pch=pch[i], legend= paste("BS", thresh.bs[i] ), bty='n', border=NA)
									dev.off()				
								})
		#
		# plot cluster size distributions per brl
		#
		clusters.dbwpat			<- lapply(seq_along(clusters.dbwpat[[1]]), function(i)
										{
											tmp			<- lapply(seq_along(clusters.dbwpat), function(j) 		clusters.dbwpat[[j]][[i]]		)
											names(tmp)	<- thresh.bs
											tmp
										})
		names(clusters.dbwpat)	<- thresh.brl
		thresh.brl.select	<- which(thresh.brl %in% c(0.06, 0.08, 0.1) )
		thresh.bs.select	<- which(thresh.bs %in% c(0.7, 0.75, 0.8, 0.85, 0.9, 0.95) )
		pch					<- 15 + seq_along(thresh.brl.select)
		lapply(thresh.brl.select,function(i)
								{												
									file				<- paste(outdir, paste(infile,"_clustpdf_",opt.brl,"_brl",thresh.brl[i],'_',gsub('/',':',outsignat),".pdf",sep=''),sep='/')				
									pdf(file=file, width=5,height=5)				
									if( with.singleton )
									{
										x					<- lapply(thresh.bs.select, function(j){	as.numeric(names(clusters.dbwpat[[i]][[j]]))	})
										y					<- lapply(thresh.bs.select, function(j){	clusters.dbwpat[[i]][[j]]	})					
									}
									else
									{
										x					<- lapply(thresh.bs.select, function(j){	as.numeric(names(clusters.dbwpat[[i]][[j]]))[-1]	})
										y					<- lapply(thresh.bs.select, function(j){	clusters.dbwpat[[i]][[j]][-1]	})
									}
									xlim				<- range( c(70, unlist(x) ) )
									cols				<- brewer.pal( length(x), "Accent")
									plot(1,1,bty='n',type='n',xlim=xlim, ylim=range(unlist(y)), xlab="patients", ylab="frequency", log='y')
									dummy<- lapply(seq_along(x),function(j)
											{
												points(x[[j]],y[[j]],type='l',col=cols[j],pch=pch[i])
												points(x[[j]],y[[j]],col=cols[j], pch=pch[i], cex=0.5)
											})				
									legend("topright", 		pch=pch[i],	legend= paste("BRL", thresh.brl[ i ]), bty='n', border=NA)
									legend("bottomright", 	fill=cols,	legend= paste("BS", thresh.bs[ thresh.bs.select ] ), bty='n', border=NA)
									dev.off()				
								})	
	}
	#
	# get FP, TP and FN, TN for several clustering thresholds
	#
	#	false positives
	fpclu.by.all			<- t( sapply(seq_along(clusters),function(i){				sapply(clusters[[i]][["tp"]], function(x) x["fpclu.by.all"] )		}))
	rownames(fpclu.by.all)	<- thresh.bs
	colnames(fpclu.by.all)	<- thresh.brl	
	fp.by.sum				<- t( sapply(seq_along(clusters),function(i){				sapply(clusters[[i]][["tp"]], function(x) x["fp.by.sum"] )		}))
	rownames(fp.by.sum)		<- thresh.bs
	colnames(fp.by.sum)		<- thresh.brl		
	fp.by.sum2				<- t( sapply(seq_along(clusters),function(i){				sapply(clusters[[i]][["tp"]], function(x) x["fp.by.sum2"] )		}))
	rownames(fp.by.sum2)	<- thresh.bs
	colnames(fp.by.sum2)	<- thresh.brl		
	#	true positives
	tp.by.all				<- t( sapply(seq_along(clusters),function(i){				sapply(clusters[[i]][["tp"]], function(x) x["tp.by.all"] )		}))
	rownames(tp.by.all)		<- thresh.bs
	colnames(tp.by.all)		<- thresh.brl	
	tpclu.by.all			<- t( sapply(seq_along(clusters),function(i){				sapply(clusters[[i]][["tp"]], function(x) x["tpclu.by.all"] )		}))
	rownames(tpclu.by.all)	<- thresh.bs
	colnames(tpclu.by.all)	<- thresh.brl
	#	true negatives
	tnn.by.sum				<- t( sapply(seq_along(clusters),function(i){				sapply(clusters[[i]][["fp"]], function(x) x["tnn.by.sum"] )	}) )
	rownames(tnn.by.sum)	<- thresh.bs
	colnames(tnn.by.sum)	<- thresh.brl
	tn.by.sum				<- t( sapply(seq_along(clusters),function(i){				sapply(clusters[[i]][["fp"]], function(x) x["tn.by.sum"] )	}) )
	rownames(tn.by.sum)		<- thresh.bs
	colnames(tn.by.sum)		<- thresh.brl	
	#	false negatives
	fpn.by.sum				<- t( sapply(seq_along(clusters),function(i){		sapply(clusters[[i]][["fp"]], function(x) x["fpn.by.sum"] )	}) )
	rownames(fpn.by.sum)	<- thresh.bs
	colnames(fpn.by.sum)	<- thresh.brl
	fp.by.all				<- t( sapply(seq_along(clusters),function(i){		sapply(clusters[[i]][["fp"]], function(x) x["fp.by.all"] )	}) )
	rownames(fp.by.all)		<- thresh.bs
	if(with.plot)
	{			
		#
		# plot TP and TN for several clustering thresholds
		#	
		colnames(fp.by.all)		<- thresh.brl	
		cols					<- diverge_hcl(nrow(fpn.by.sum), h = c(246, 40), c = 96, l = c(65, 90))
		names(cols)				<- thresh.bs
		#	 
		hivc.clu.plot.tptn(	fpn.by.sum, tp.by.all, paste(outdir, paste(outfile,"_tpbyall_fpnbysum_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
							cols, xlab= "#FP (among all)", ylab= "%TP (among all)", xlim= range(c(0,16,fpn.by.sum)))
		#
		hivc.clu.plot.tptn(	fp.by.all, tp.by.all, paste(outdir, paste(outfile,"_tpbyall_fpbyall_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
							cols, xlab= "%FP (among all)", ylab= "%TP (among all)", xlim= range(c(0,0.01,fp.by.all)))
		#
		hivc.clu.plot.tptn(	fpn.by.sum, tpclu.by.all, paste(outdir, paste(outfile,"_tpclubyall_fpnbysum_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
							cols, xlab= "#FP (among all)", ylab= "%TP (among clu)", xlim= range(c(0,16,fpn.by.sum)))
		#
		hivc.clu.plot.tptn(	fp.by.all, tpclu.by.all, paste(outdir, paste(outfile,"_tpclubyall_fpbyall_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
							cols, xlab= "%FP (among all)", ylab= "%TP (among clu)", xlim= range(c(0,0.01,fp.by.all)))
		#
		hivc.clu.plot.tptn(	fpn.by.sum, clusters.cov.epidemic, paste(outdir, paste(outfile,"_covepi_fpnbysum_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
							cols, xlab= "#FP (among all)", ylab= "%coverage (of epi)", xlim= range(c(0,16,fpn.by.sum)))
		#
		hivc.clu.plot.tptn(	fp.by.all, clusters.cov.epidemic, paste(outdir, paste(outfile,"_covepi_fpbyall_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
							cols, xlab= "%FP (among all)", ylab= "%coverage (of epi)", xlim= range(c(0,0.01,fp.by.all)))
		#
		hivc.clu.plot.tptn(	fpn.by.sum, clusters.cov.patindb, paste(outdir, paste(outfile,"_covpat_fpnbysum_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
							cols, xlab= "#FP (among all)", ylab= "%coverage (of db)", xlim= range(c(0,16,fpn.by.sum)))
		#
		hivc.clu.plot.tptn(	fp.by.all, clusters.cov.patindb, paste(outdir, paste(outfile,"_covpat_fpbyall_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
							cols, xlab= "%FP (among all)", ylab= "%coverage (of db)", xlim= range(c(0,0.01,fp.by.all)))
		#
		hivc.clu.plot.tptn(	fpn.by.sum, clusters.cov.seqindb, paste(outdir, paste(outfile,"_covseq_fpnbysum_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
							cols, xlab= "#FP (among all)", ylab= "%seq coverage (of db)", xlim= range(c(0,16,fpn.by.sum)))
		#
		hivc.clu.plot.tptn(	fp.by.all, clusters.cov.seqindb, paste(outdir, paste(outfile,"_covseq_fpbyall_",gsub('/',':',outsignat),".pdf",sep=''),sep='/'), 
							cols, xlab= "%FP (among all)", ylab= "%seq coverage (of db)", xlim= range(c(0,0.01,fp.by.all)))
	}
	
	list(	clusters				= clusters,
			clusters.cov.epidemic	= clusters.cov.epidemic, 
			clusters.cov.patindb	= clusters.cov.patindb, 
			clusters.cov.seqindb	= clusters.cov.seqindb,
			clusters.dbwpat			= clusters.dbwpat,
			fpn.by.sum				= fpn.by.sum,
			fp.by.all				= fp.by.all,
			tp.by.all				= tp.by.all, 
			tpclu.by.all			= tpclu.by.all,
			tnn.by.sum				= tnn.by.sum,
			tn.by.sum				= tn.by.sum,
			fp.by.sum				= fp.by.sum,
			fp.by.sum2				= fp.by.sum2,
			fpclu.by.all			= fpclu.by.all
			)
}
######################################################################################
hivc.prog.recombination.process.3SEQ.output<- function()
{	
	verbose		<- 1
	resume		<- 1
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	
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
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}	
	if(verbose)
	{
		print(verbose)
		print(resume)
		print(indir)		
		print(infile)			
		print(insignat)			
	}
	file		<- paste(indir,'/',infile,'_', gsub('/',':',insignat),".R",sep='')
	if(verbose)	cat(paste("\nload file ",file))			
	load(file)
	#loaded seq.PROT.RT
	if(resume)
	{
		file		<- paste(indir,'/',infile,"_3seq_", gsub('/',':',insignat),".R",sep='')		
		options(show.error.messages = FALSE)		
		if(verbose)	cat(paste("\ntry to resume file ",file))
		readAttempt<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)	cat(paste("\nresumed file ",file))			
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{
		if(verbose)	cat(paste("\ngenerate file ",file))			
		tmp			<- list.files(indir, pattern="3seq$")
		files		<- tmp[grepl(infile, tmp, fixed=1) & grepl(gsub('/',':',insignat), tmp)]
		tmp			<- sapply(strsplit(files,'_'), function(x) rev(x)[6] )		#see if filename has '_xx-xx_' string at pos 6 from end  
		files		<- files[grepl('-', tmp)]		
		if(verbose)	cat(paste("\nFound 3seq output files matching infile and insignat, n=",length(files)))
		#	figure out if any consecutive files are missing
		tmp			<- tmp[grepl('-', tmp)]
		df.seqchunks<- as.data.table( t( sapply( strsplit(tmp, '-'), as.numeric) ) )
		setnames(df.seqchunks, c("V1","V2"), c("start","end"))
		setkey(df.seqchunks, start)
		tmp			<- subset(df.seqchunks, df.seqchunks[-1,start]-1 != df.seqchunks[-nrow(df.seqchunks),end] )
		if(nrow(tmp)){		print(tmp); stop("Found missing sequence chunks")		}
		#	determine non-empty files and number of columns		
		tmp			<- sapply(files, function(x)
				{
					tmp	<- count.fields(paste(indir, '/', x, sep=''), skip=1, sep='\t')
					ifelse(length(tmp), max(tmp), 0)			
				})
		col.max		<- max(tmp)		
		files		<- files[tmp>0]
		if(verbose)	cat(paste("\nFound non-empty 3seq files, n=",length(files)))
		#	read non-empty files
		col.names		<- c("parent1","parent2","child","m","n","k","p","[p_max]","HS?","log(p)","DS(p)","DS(p)","min_rec_length")
		if(col.max>length(col.names))
			col.names	<- c(col.names, paste("bp",seq_len( col.max-length(col.names) ),sep=''))
		df.recomb		<- lapply(files, function(x)
				{	
					if(verbose)	cat(paste("\nprocess file", x))
					tmp					<- read.delim(paste(indir, '/', x, sep=''), skip=1, header=0, fill=1, col.names=col.names, stringsAsFactors=0, strip.white=T)					
					as.data.table(tmp)			
				})
		df.recomb		<- rbindlist(df.recomb)
		#	extract FASTASampleCodes from job output
		tmp				<- df.recomb[, list(m1= regexpr("-[",parent1,fixed=1), m2=regexpr("-[",parent2,fixed=1), mc=regexpr("-[",child,fixed=1)) ]
		if(any(tmp<1))	stop("Unexpected parent1 or parent2. Parsing error?")
		df.recomb		<- cbind( df.recomb, tmp )
		df.recomb[,dummy:=seq_len(nrow(df.recomb))]
		set(df.recomb, NULL, "parent1", df.recomb[,substr(parent1, 1, m1-1), by="dummy"][,V1])
		set(df.recomb, NULL, "parent2", df.recomb[,substr(parent2, 1, m2-1), by="dummy"][,V1])
		set(df.recomb, NULL, "child", df.recomb[,substr(child, 1, mc-1), by="dummy"][,V1])
		#	extract unique recombinants from job output
		setkey(df.recomb,child)
		df.recomb		<- unique(df.recomb)
		if(verbose)	cat(paste("\nFound recombinant sequences, n=",nrow(df.recomb)))
		#	order parents
		tmp				<- df.recomb[, {
					z<- sort(c(parent1, parent2))
					list(parent1=z[1], parent2=z[2])
				} ,by="dummy"]
		df.recomb		<- df.recomb[, setdiff( colnames(df.recomb), c("parent1","parent2","p","X.p_max.","HS.","DS.p.","m1","m2","mc") ), with=F]
		df.recomb		<- merge(tmp, df.recomb, by="dummy")
		#	check if triplets are unique
		tmp				<- df.recomb[, {
					z<- sort(c(parent1, parent2, child))
					list(parent1=z[1], parent2=z[2], child=z[3])
				} ,by="dummy"]
		setkeyv(tmp, c("parent1","parent2","child"))
		tmp				<- unique(tmp)
		if(nrow(tmp)<nrow(df.recomb))	warning("triplets in df.recomb not unique")
		#
		tmp				<- subset(df.recomb, select=c(parent1,parent2,child))
		setkey(tmp, parent1, parent2)
		tmp				<- unique(tmp)
		tmp[, parentpair:=seq_len(nrow(tmp))]				
		if(verbose)	cat(paste("\nNumber of unique parent pairs, n=",nrow(tmp)))
		#
		tmp				<- unique(c(df.recomb[, parent1], df.recomb[, parent2]))
		if(verbose)	cat(paste("\nNumber of unique parent sequences, n=",length(tmp)))
		#
		if(length(unique(df.recomb[,child]))!=nrow(df.recomb))	stop("Unexpected duplicate children - select parents with smallest p?")
		#
		tmp				<- subset(df.recomb, select=c(parent1,parent2))
		setkey(tmp, parent1, parent2)
		tmp				<- unique(tmp)
		tmp[, parentpair:=seq_len(nrow(tmp))]				
		if(verbose)	cat(paste("\nNumber of unique parent pairs, n=",nrow(tmp)))
		df.recomb		<- merge(df.recomb, tmp, by=c("parent1","parent2"))
		setnames(df.recomb,c("log.p.","DS.p..1"),c("logp","adjp"))
		#	candidate breakpoints bp1 etc are all overlapping, consider only bp1 for simplicity	
		#	evaluate midpoint of breapoint regions		
		tmp<- strsplit(df.recomb[,bp1], ',')
		if(any(sapply(tmp, length )!=2)) 	stop("\nunexpected bp1: missing ','")		
		df.recomb[, bp1.1:= sapply(tmp, "[", 1)]
		df.recomb[, bp1.2:= sapply(tmp, "[", 2)]
		set(df.recomb, NULL, "bp1.1", round( sapply(strsplit(df.recomb[,bp1.1], '-'), function(x)	mean(as.numeric(x))	) ) )
		set(df.recomb, NULL, "bp1.2", round( sapply(strsplit(df.recomb[,bp1.2], '-'), function(x)	min(as.numeric(x))	) ) )
		df.recomb	<- merge(df.recomb, df.recomb[, list(child.len=max(which(as.character(seq.PROT.RT[child,])!='-'))), by="dummy"], by="dummy" )
		if(any(df.recomb[, child.len-bp1.2]<0))	stop("\nunexpected breakpoint past end of child sequence")
		df.recomb[, child.start:= 1]
		tmp			<- which( df.recomb[, bp1.1<10] )
		set(df.recomb, tmp, "child.start", df.recomb[tmp, bp1.1])
		tmp			<- which( df.recomb[, child.len-bp1.2<25] )
		set(df.recomb, tmp, "child.len", as.integer(df.recomb[tmp, bp1.2]))
		#	save potential recombinants
		file		<- paste(indir,'/',infile,"_3seq_", gsub('/',':',insignat),".R",sep='')
		if(verbose)	cat(paste("\nSave candidate triplets to",file))
		save(df.recomb, file=file)		
	}
	if(0)
	{
		setnames(df.recomb, "dummy", "triplet.id")
		#
		#	get candidate recombinant sequences
		#
		df.recombseq	<- data.table( FASTASampleCode= unique( c(df.recomb[, parent1], df.recomb[, parent2], df.recomb[, child]) ) )
		df.recombseq	<- df.recombseq[ , 	{
												tmp<- subset(df.recomb, parent1==FASTASampleCode | parent2==FASTASampleCode | child==FASTASampleCode )
												list( n.triplets=nrow(tmp), min.p=min(tmp[,adjp]), med.p=median(tmp[,adjp]), triplet.id=tmp[,triplet.id] )
											} , by="FASTASampleCode"]
		setkey(df.recombseq, n.triplets, FASTASampleCode)
		#plot( df.recombseq[,n.triplets], df.recombseq[,min.p], pch=18 )
		
		#
		#determine how many other triplet sequences there are for an m candidate
		#
		df.mrecombseq	<- subset(df.recombseq, n.triplets>2)[, {
					tmp			<- triplet.id
					tmp			<- subset(df.recomb, triplet.id%in%tmp)
					triplet.seq	<- setdiff( unique( c(tmp[, parent1], tmp[, parent2], tmp[, child]) ), FASTASampleCode)
					list( triplet.seq=triplet.seq, triplet.seq.n=length(triplet.seq)  )
				}, by="FASTASampleCode"]
		df.mrecombbp	<- subset(df.recombseq, n.triplets>2)[, {
					tmp			<- triplet.id
					tmp			<- subset(df.recomb, triplet.id%in%tmp)
					list( bp1= tmp[,bp1], bp1.1= tmp[,bp1.1], bp1.2= tmp[,bp1.2] )																	
				}, by="FASTASampleCode"]														
		setkeyv(df.mrecombbp, c("FASTASampleCode", "bp1.1", "bp1.2"))												
		#	breakpoints among mrecombinants are not necessarily the same
		print(df.mrecombbp, n=250)										
	}
	df.recomb
}
######################################################################################
hivc.prog.recombination.plot.incongruence<- function()
{
	require(RColorBrewer)
	require(ape)
	
	verbose		<- 1
	resume		<- 1
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	
	id			<- NA		
	bs.n		<- 500
	select		<- ''
	#select		<- 'ng2'
	
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
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									tripletid= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) id<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									select= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) select<- tmp[1]
	}	
	if(verbose)
	{
		print(verbose)		
		print(indir)		
		print(infile)			
		print(insignat)
		print(id)		
		print(bs.n)		
		print(select)
	}
	cols				<- brewer.pal(6,"Paired")
	pattern				<- c("^in_parent1","^out_parent1","^in_parent2","^out_parent2","^in_child","^out_child")
	cex					<- 0.6
	thresh.nodesupport	<- 0.6
	edge.length.outliers<- 0.5
	
	if(is.na(id))
	{
		#	read candidate triplets
		argv				<<-	hivc.cmd.recombination.process.3SEQ.output(indir, infile, insignat, resume=1, verbose=1) 
		argv				<<- unlist(strsplit(argv,' '))
		df.recomb			<- hivc.prog.recombination.process.3SEQ.output()
		setnames(df.recomb, "dummy", "triplet.id")
		setkey(df.recomb, triplet.id)
		#	get candidate recombinant sequences
		df.recombseq	<- data.table( FASTASampleCode= unique( c(df.recomb[, parent1], df.recomb[, parent2], df.recomb[, child]) ) )
		df.recombseq	<- df.recombseq[ , 	{
					tmp<- subset(df.recomb, parent1==FASTASampleCode | parent2==FASTASampleCode | child==FASTASampleCode )
					list( n.triplets=nrow(tmp), min.p=min(tmp[,adjp]), med.p=median(tmp[,adjp]), triplet.id=tmp[,triplet.id] )
				} , by="FASTASampleCode"]
		setkey(df.recombseq, FASTASampleCode)
		#	select unclear triplets	
		if(select=="ng2")
		{
			tmp				<- unique(subset(df.recombseq, n.triplets>2)[, triplet.id])			
			df.recomb		<- df.recomb[J(setdiff(df.recomb[,triplet.id], tmp )),]
			if(verbose)	cat(paste("\nSelected unclear triplets, n=", nrow(df.recomb)))
		}
		#	select triplets that have a given FASTASampleCode
		if(select %in% unique(subset(df.recombseq, n.triplets>2)[, FASTASampleCode]))
		{
			tmp				<- subset(df.recombseq, n.triplets>2)[J(select),]
			df.recomb		<- df.recomb[J(tmp[,triplet.id]),]
		}
		#	analyze FASTASampleCodes that occur in >2 triplets
		if(select=="g2")
		{
			tmp				<- unique(subset(df.recombseq, n.triplets>2)[, FASTASampleCode])
			dummy			<- lapply(tmp, function(x)
						{
							argv			<<- hivc.cmd.recombination.plot.incongruence(indir, infile, insignat, prog= PR.RECOMB.PLOTINCONGRUENCE, opt.select=x,verbose=1)
							argv			<<- unlist(strsplit(argv,' '))
							hivc.prog.recombination.plot.incongruence()			
						})
			stop()
		}
		
		#	read available checks for triplets
		files				<- list.files(indir, pattern=".newick$")
		files				<- files[ grepl(paste('_3seqcheck_',sep=''), files) & grepl(infile, files, fixed=1) & grepl(gsub('/',':',insignat), files) ]
		tmp					<- regexpr("3seqcheck_",files)
		tmp					<- sapply(seq_along(files), function(i){		substr(files[i],tmp[i],nchar(files[i]))		})
		files.df			<- data.table(	file=files, 
											triplet.id= sapply(strsplit(tmp, '_'), function(x)		substr(x[[2]],3,nchar(x[[2]]))	),
											region= sapply(strsplit(tmp, '_'), function(x)		substr(x[[3]],2,nchar(x[[3]]))	)		)
		set(files.df, NULL, "triplet.id", as.numeric(files.df[,triplet.id]))							
		setkey(files.df, triplet.id, region)
		if(verbose)	cat(paste("\nFound files, n=", nrow(files.df)))
		files.df			<- merge(files.df, files.df[,list(region.n= length(region)), by="triplet.id"], by="triplet.id")
		files.df			<- subset(files.df, region.n>1)
		if(verbose)	cat(paste("\nFound checked triplets, n=", nrow(files.df)/2))
		if(verbose) cat(paste("\nTriplets still to check=",paste(setdiff( df.recomb[,triplet.id],unique(files.df[,triplet.id]) ), collapse=', ')))
		#	select candidate triplets for which checks available
		files.df			<- merge(files.df, df.recomb,by="triplet.id")
		if(verbose)	cat(paste("\nSelected triplets for plotting, n=", nrow(files.df)/2))
		setkey(files.df, parentpair, adjp)
		#setkey(files.df, adjp)
		
		file				<- paste(indir,'/',infile,"_3seqcheck_examlbs",bs.n,'_',select,'_',gsub('/',':',insignat),".pdf",sep='')
		if(verbose)	cat(paste("\nPlot both phylogenies to", file))		
		pdf(file, width=12, height=6)
		def.par 			<- par(no.readonly = TRUE)
		par(mar=c(2,0.5,1,0))
		layout(matrix(c(1,1,2,2), 2, 2))					
		dummy				<- lapply( unique( files.df[, triplet.id] ), function(z)
				{
					#x<- subset(files.df, triplet.id==58)
					x<- subset(files.df, triplet.id==z)
					#	read tree corresponding to 'in' recombinant region
					tmp					<- paste(indir,'/',x[1, file],sep='')
					if(verbose)	cat(paste("\nRead file", tmp))
					ph.in				<- ladderize( read.tree(tmp) )		
					tmp					<- as.numeric( ph.in$node.label )
					tmp[is.na(tmp)]		<- 0 
					ph.in$node.label	<- tmp/100
					#	read tree corresponding to 'out' recombinant region
					tmp					<- paste(indir,'/',x[2, file],sep='')		
					if(verbose)	cat(paste("\nRead file", tmp))
					ph.out				<- ladderize( read.tree(tmp) )		
					tmp					<- as.numeric( ph.out$node.label )
					tmp[is.na(tmp)]		<- 0 
					ph.out$node.label	<- tmp/100
					#	remove outlier filler sequences
					outliers			<- c()
					tmp					<- which( ph.in$edge.length>edge.length.outliers )
					if(length(tmp))
					{
						tmp				<- ph.in$edge[tmp,2]
						outliers		<- ph.in$tip.label[ tmp[tmp<Ntip(ph.in)] ]
					}
					if(length(tmp))
					{						
						tmp				<- which( ph.out$edge.length>1 )
						tmp				<- ph.out$edge[tmp,2]
						outliers		<- c(outliers, ph.out$tip.label[ tmp[tmp<Ntip(ph.out)] ])
					}
					if(length(outliers))
					{
						ph.in			<- drop.tip(ph.in, outliers)
						ph.out			<- drop.tip(ph.out, outliers)
					}
					#	show half way stable subtrees -- ideally would contain triplet sequences
					clustering.in		<- hivc.clu.clusterbythresh(ph.in, thresh.nodesupport=thresh.nodesupport, nodesupport=ph.in$node.label, retval="all")
					clustering.out		<- hivc.clu.clusterbythresh(ph.out, thresh.nodesupport=thresh.nodesupport, nodesupport=ph.out$node.label, retval="all")
					#	show filler sequences in different color
					tip.color.in		<- rep("black", Ntip(ph.in))	
					for(i in seq_along(pattern))
					{
						tmp						<- which( grepl(pattern[i],ph.in$tip) )
						if(length(tmp))	
							tip.color.in[ tmp ]	<- cols[i]
					}
					tip.color.out		<- rep("black", Ntip(ph.out))	
					for(i in seq_along(pattern))
					{
						tmp						<- which( grepl(pattern[i],ph.out$tip) )
						if(length(tmp))	
							tip.color.out[ tmp ]<- cols[i]
					}	
					#	plot both phylogenies side by side	
					dummy<- hivc.clu.plot(ph.in, clustering.in[["clu.mem"]], show.tip.label=T, cex.nodelabel=cex, cex.edge.incluster=2*cex, no.margin=F, tip.color=tip.color.in)
					mtext(paste("parent.pair=",x[,parentpair]," triplet.id=",x[1,triplet.id]," log10p=",round(log10(x[1,adjp]),d=2)," region 'in' length=",x[1,bp1.2-bp1.1]), side = 3, cex=cex)
					axisPhylo(cex=cex)
					dummy<- hivc.clu.plot(ph.out, clustering.out[["clu.mem"]], show.tip.label=T, cex.nodelabel=cex, cex.edge.incluster=2*cex, no.margin=F, tip.color=tip.color.out)		
					mtext(paste("parent.pair=",x[,parentpair]," triplet.id=",x[1,triplet.id]," log10p=",round(log10(x[1,adjp]),d=2)," region 'out' length=",x[1,child.len-child.start-bp1.2+bp1.1]), side = 3, cex=cex)
					axisPhylo(cex=cex)
				})
		par(def.par)
		dev.off()
	
		#x<- "R12-15108"
		#subset(df.recomb, parent1==x | parent2==x | child==x )
	}
	if(!is.na(id))
	{
		#	read tree corresponding to 'in' recombinant region
		file				<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rIn_examlbs",bs.n,'_',gsub('/',':',insignat),".newick",sep='')
		if(verbose)	cat(paste("\nRead file", file))
		ph.in				<- ladderize( read.tree(file) )		
		tmp					<- as.numeric( ph.in$node.label )
		tmp[is.na(tmp)]		<- 0 
		ph.in$node.label	<- tmp/100
		#	read tree corresponding to 'out' recombinant region
		file				<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rOut_examlbs",bs.n,'_',gsub('/',':',insignat),".newick",sep='')
		if(verbose)	cat(paste("\nRead file", file))
		ph.out				<- ladderize( read.tree(file) )		
		tmp					<- as.numeric( ph.out$node.label )
		tmp[is.na(tmp)]		<- 0 
		ph.out$node.label	<- tmp/100
		#	remove outlier filler sequences
		outliers			<- c()
		tmp					<- which( ph.in$edge.length>edge.length.outliers )
		if(length(tmp))
		{
			tmp				<- ph.in$edge[tmp,2]
			outliers		<- ph.in$tip.label[ tmp[tmp<Ntip(ph.in)] ]
		}
		if(length(tmp))
		{						
			tmp				<- which( ph.out$edge.length>1 )
			tmp				<- ph.out$edge[tmp,2]
			outliers		<- c(outliers, ph.out$tip.label[ tmp[tmp<Ntip(ph.out)] ])
		}
		if(length(outliers))
		{
			ph.in			<- drop.tip(ph.in, outliers)
			ph.out			<- drop.tip(ph.out, outliers)
		}
		#	show half way stable subtrees -- ideally would contain triplet sequences
		clustering.in		<- hivc.clu.clusterbythresh(ph.in, thresh.nodesupport=thresh.nodesupport, nodesupport=ph.in$node.label, retval="all")
		clustering.out		<- hivc.clu.clusterbythresh(ph.out, thresh.nodesupport=thresh.nodesupport, nodesupport=ph.out$node.label, retval="all")
		#	show filler sequences in different color
		tip.color.in		<- rep("black", Ntip(ph.in))	
		for(i in seq_along(pattern))
		{
			tmp						<- which( grepl(pattern[i],ph.in$tip) )
			if(length(tmp))	
				tip.color.in[ tmp ]	<- cols[i]
		}
		tip.color.out		<- rep("black", Ntip(ph.out))	
		for(i in seq_along(pattern))
		{
			tmp						<- which( grepl(pattern[i],ph.out$tip) )
			if(length(tmp))	
				tip.color.out[ tmp ]<- cols[i]
		}	
		#	plot both phylogenies side by side
		file				<- paste(indir,'/',infile,"_3seqcheck_id",id,"_examlbs",bs.n,'_',gsub('/',':',insignat),".pdf",sep='')
		if(verbose)	cat(paste("\nPlot both phylogenies to", file))
		pdf(file, width=8, height=6)
		def.par 			<- par(no.readonly = TRUE)
		layout(matrix(c(1,1,2,2), 2, 2))		
		dummy<- hivc.clu.plot(ph.in, clustering.in[["clu.mem"]], show.tip.label=T, cex.nodelabel=cex, cex.edge.incluster=2*cex, no.margin=F, tip.color=tip.color.in)
		mtext(paste("parent.pair=",x[,parentpair]," triplet.id=",x[1,triplet.id]," log10p=",round(log10(x[1,adjp]),d=2)," region 'in' length=",x[1,bp1.2-bp1.1]), side = 3, cex=cex)
		axisPhylo(cex=cex)
		dummy<- hivc.clu.plot(ph.out, clustering.out[["clu.mem"]], show.tip.label=T, cex.nodelabel=cex, cex.edge.incluster=2*cex, no.margin=F, tip.color=tip.color.out)		
		mtext(paste("parent.pair=",x[,parentpair]," triplet.id=",x[1,triplet.id]," log10p=",round(log10(x[1,adjp]),d=2)," region 'out' length=",x[1,child.len-child.start-bp1.2+bp1.1]), side = 3, cex=cex)
		axisPhylo(cex=cex)
		par(def.par)
		dev.off()
	}
}
######################################################################################
hivc.prog.recombination.check.candidates<- function()
{	
	require(ape)
	verbose		<- 1
	resume		<- 0
	indir		<- paste(DATA,"tmp",sep='/')		
	infile		<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"			
	insignat	<- "Thu_Aug_01_17/05/23_2013"
	
	id			<- 51
	seq.select.n<- 10
	bs.from		<- 0
	bs.to		<- 499
	bs.n		<- 500
	
	hpc.walltime<- 36
	hpc.mem		<- "600mb"
	hpc.nproc	<- 1		
	hpc.q		<- "pqeph"
	
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
						{	switch(substr(arg,2,10),
									tripletid= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) id<- tmp[1]
	}	
	if(verbose)
	{
		print(verbose)
		print(resume)
		print(indir)		
		print(infile)			
		print(insignat)
		print(id)
		print(seq.select.n)
		print(bs.from)
		print(bs.to)
		print(bs.n)		
	}
	if(resume)
	{
		tmp			<- 1
		file		<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rIn_",gsub('/',':',insignat),".R",sep='')
		options(show.error.messages = FALSE)		
		if(verbose)	cat(paste("\ntry to load file ",file))
		readAttempt	<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)		cat(paste("\nloaded file=",file))
		if(inherits(readAttempt, "try-error"))					tmp		<- 0
		
		file		<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rOut_",gsub('/',':',insignat),".R",sep='')
		options(show.error.messages = FALSE)		
		if(verbose)	cat(paste("\ntry to load file ",file))
		readAttempt	<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)		cat(paste("\nloaded file=",file))
		if(inherits(readAttempt, "try-error"))					tmp		<- 0		
	}
	if(!resume || tmp==0)
	{
		file		<- paste(indir,'/',infile,'_', gsub('/',':',insignat),".R",sep='')
		if(verbose)	cat(paste("\nload file ",file))			
		load(file)
		#	loaded seq.PROT.RT
		file		<- paste(indir,'/',infile,"_3seq_", gsub('/',':',insignat),".R",sep='')		
		options(show.error.messages = FALSE)		
		if(verbose)	cat(paste("\ntry to load file ",file))
		readAttempt	<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(inherits(readAttempt, "try-error"))	stop(paste("\nCannot find 3SEQ file, run hivc.prog.recombination.process.3SEQ.output?, file=",file))			
		#	loaded df.recomb
		
		#
		#	process triplet for dummy id	
		#	
		df.recomb	<- subset(df.recomb, dummy==id)
		if(verbose)	cat(paste("\nprocess triplet number",id))
		if(verbose)	print(df.recomb)
		#	create sequence matrices corresponding to the two breakpoint regions 
		seq.in		<- seq.PROT.RT[,seq.int(df.recomb[,bp1.1],df.recomb[,bp1.2])]
		seq.out		<- if(df.recomb[,child.start]<df.recomb[,bp1.1]-1) seq.int(df.recomb[,child.start],df.recomb[,bp1.1]-1) else numeric(0) 
		seq.out		<- if(df.recomb[,bp1.2]+1<df.recomb[,child.len]) c(seq.out,seq.int(df.recomb[,bp1.2]+1, df.recomb[,child.len]))	else 	seq.out
		seq.out		<- seq.PROT.RT[,seq.out]
		seq.select.f<- ifelse(min(ncol(seq.out),ncol(seq.in))<150, 10, 10)
		if(verbose)	cat(paste("\nsetting inflation factor to",seq.select.f))
		seq.select.n<- seq.select.n * seq.select.f
		#	select background sequences for child based on sequence similarity
		if(verbose)	cat(paste("\ncompute genetic distances for parent1 parent2 child"))
		tmp				<- which( rownames(seq.PROT.RT)==df.recomb[,child] )		
		dummy			<- 0				
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.in))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.in[tmp,], seq.in[i,], ncol(seq.in), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- data.table( FASTASampleCode=rownames(seq.in)[-tmp] , dist=seq.dist, group="child", region="in" ) 		
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.out))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.out[tmp,], seq.out[i,], ncol(seq.out), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df, data.table( FASTASampleCode=rownames(seq.out)[-tmp] , dist=seq.dist, group="child", region="out" )) 
		#	select background sequences for parent1 based on sequence similarity
		tmp				<- which( rownames(seq.PROT.RT)==df.recomb[,parent1] )				
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.in))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.in[tmp,], seq.in[i,], ncol(seq.in), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df,data.table( FASTASampleCode=rownames(seq.in)[-tmp] , dist=seq.dist, group="parent1", region="in" )) 		
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.out))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.out[tmp,], seq.out[i,], ncol(seq.out), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df, data.table( FASTASampleCode=rownames(seq.out)[-tmp] , dist=seq.dist, group="parent1", region="out" )) 
		#	select background sequences for parent2 based on sequence similarity
		tmp				<- which( rownames(seq.PROT.RT)==df.recomb[,parent2] )				
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.in))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.in[tmp,], seq.in[i,], ncol(seq.in), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df,data.table( FASTASampleCode=rownames(seq.in)[-tmp] , dist=seq.dist, group="parent2", region="in" )) 		
		seq.dist					<- 1 - sapply(seq_len(nrow(seq.out))[-tmp],function(i){		.C("hivc_dist_ambiguous_dna", seq.out[tmp,], seq.out[i,], ncol(seq.out), dummy )[[4]]			})	
		seq.dist[is.nan(seq.dist)]	<- Inf
		seq.df			<- rbind(seq.df, data.table( FASTASampleCode=rownames(seq.out)[-tmp] , dist=seq.dist, group="parent2", region="out" ))
		#	first pass:
		#	select closest FASTASampleCode by group and region to verify recombination breakpoint by phylogenetic incongruence
		setkeyv(seq.df, c("group","region","dist"))
		seq.df			<- subset(seq.df, dist>0)
		if(verbose)	cat(paste("\nFound related sequences with dist>0, n=",nrow(seq.df)))
		#	for each group and region, select the n closest sequence names (non-unique)
		seq.df			<- seq.df[	,	list(FASTASampleCode=FASTASampleCode[seq_len(seq.select.n)], dist=dist[seq_len(seq.select.n)]), by=c("group","region")]
		#	keep each sequence name once, for the group it is closest to		
		seq.df			<- seq.df[, {
										tmp<- which.min(dist)
										list(dist= dist[tmp], group=group[tmp], region=region[tmp])
									}, by=c("FASTASampleCode")]
		setkeyv(seq.df, c("group","region","dist"))					
		if(verbose)	cat(paste("\ndetermined candidates for balancing filler sequences, n=",nrow(seq.df)))
		#	select unique sequences
		tmp				<- c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child], seq.df[, FASTASampleCode] )
		tmp				<- seq.unique(seq.in[tmp,])
		seq.df			<- merge( data.table(FASTASampleCode=rownames(tmp)), seq.df, by="FASTASampleCode" )
		tmp				<- c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child], seq.df[, FASTASampleCode] )
		tmp				<- seq.unique(seq.out[tmp,])
		seq.df			<- merge( data.table(FASTASampleCode=rownames(tmp)), seq.df, by="FASTASampleCode" )
		seq.df			<- subset(seq.df, !FASTASampleCode%in%c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ) )		
		if(verbose)	cat(paste("\nfound candidates for filler sequences that are unique on both recombinant regions, n=",nrow(seq.df)))
		if(verbose)	print( seq.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )
		#		
		seq.select.n	<- seq.select.n/seq.select.f 
		#	select 'seq.select.n' unique filler sequences for 'in' region, balancing by group as much as possible
		if(verbose)	cat(paste("\nSelect filler sequences for recombinant region 'in'"))
		tmp						<- subset( seq.df, region=='in' )[,FASTASampleCode]
		tmp						<- seq.unique(seq.in[tmp ,])
		seq.in.df				<- merge( data.table(FASTASampleCode=rownames(tmp)), subset(seq.df,region=="in"), by="FASTASampleCode" )
		setkey(seq.in.df, dist)
		if(nrow(seq.in.df)<seq.select.n)	cat(paste("\ncan only select less than the requested number of sequences, n=",nrow(seq.in.df)))
		tmp						<- rbind( seq.in.df, data.table(FASTASampleCode=NA, dist=NA, group=c("child","parent1","parent2"), region=NA) )
		seq.in.order			<- tmp[	,	list(n=length(na.omit(FASTASampleCode))) ,by=c("group")]		
		seq.in.order			<- seq.in.order[order(n),]
		overflow				<- 0
		ans						<- data.table(FASTASampleCode=NA, dist=NA, group=NA, region=NA)
		for(x in seq.in.order[,group])
		{			
			#print(x)
			tmp					<- subset(seq.in.df, group==x)
			#print(tmp)
			ans					<- rbind(tmp[seq_len( min(seq.select.n+overflow, nrow(tmp)) ),], ans )
			overflow			<- ifelse(seq.select.n+overflow<nrow(tmp), 0, seq.select.n+overflow-nrow(tmp))
		}
		if(overflow>0)	cat(paste("\nNot as many filler sequences as requested for recombinant region 'in', n=",nrow(seq.out.df)))
		seq.in.df				<- ans[-nrow(ans),]
		if(verbose)	cat(paste("\nSelected balancing sequences for recombinant region 'in', n=",nrow(seq.in.df)))
		if(verbose)	print( seq.in.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )
		#
		#	select 'seq.select.n' unique filler sequences for 'out' region, balancing by group as much as possible
		#
		if(verbose)	cat(paste("\nSelect sequences for recombinant region 'out'"))
		tmp						<- subset( seq.df, region=='out' )[,FASTASampleCode] 
		tmp						<- seq.unique(seq.out[tmp,])		
		seq.out.df				<- merge( data.table(FASTASampleCode=rownames(tmp)), subset(seq.df,region=="out"), by="FASTASampleCode" )
		setkey(seq.out.df, dist)
		if(nrow(seq.out.df)<seq.select.n)	cat(paste("\ncan only select less than the requested number of sequences, n=",nrow(seq.out.df)))
		tmp						<- rbind( seq.out.df, data.table(FASTASampleCode=NA, dist=NA, group=c("child","parent1","parent2"), region=NA) )
		seq.out.order			<- tmp[	,	list(n=length(FASTASampleCode)) ,by=c("group")]		
		seq.out.order			<- seq.out.order[order(n),]
		overflow				<- 0
		ans						<- data.table(FASTASampleCode=NA, dist=NA, group=NA, region=NA)		
		for(x in seq.out.order[,group])
		{			
			#print(x)
			tmp					<- subset(seq.out.df, group==x)
			#print(tmp)
			ans					<- rbind(tmp[seq_len( min(seq.select.n+overflow, nrow(tmp)) ),], ans )
			overflow			<- ifelse(seq.select.n+overflow<nrow(tmp), 0, seq.select.n+overflow-nrow(tmp))
		}
		if(overflow>0)	cat(paste("\nNot as many filler sequences as requested for recombinant region 'out', n=",nrow(seq.out.df)))
		seq.out.df				<- ans[-nrow(ans),]
		if(verbose)	cat(paste("\nSelected balancing sequences for recombinant region 'out', n=",nrow(seq.out.df)))
		if(verbose)	print( seq.out.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )
		#	combine unique filler sequences
		seq.df					<- rbind( seq.in.df, seq.out.df )
		if(verbose)	cat(paste("\nSelected balancing set of closest filler sequences"))
		if(verbose)	print( seq.df[	,	list(n=length(FASTASampleCode)) ,by=c("group","region")] )		
		if( length(unique( seq.df[, FASTASampleCode] ))!=nrow(seq.df) )		stop("Unexpected non-unique sequence names")
		if( nrow(seq.unique( seq.in[ seq.df[, FASTASampleCode], ] ))!=nrow(seq.df) ) stop("Unexpected non-unique 'in' sequences")
		if( nrow(seq.unique( seq.out[ seq.df[, FASTASampleCode], ] ))!=nrow(seq.df) ) stop("Unexpected non-unique 'out' sequences")
		#	could be that the triplet sequences are not unique among each other 
		#	if so, make change to one of the triplet sequences
		seq.select				<- c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child], seq.df[, FASTASampleCode] )
		tmp						<- seq.unique( seq.in[ seq.select, ] )
		if( !length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )
			seq.in				<- tmp
		if( length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )
		{
			tmp					<- setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp) )		#name of sequence in triplet that is identical with one other sequence in triplet
			if(verbose)	cat(paste("\nFound identical triplet sequence for region 'in'", tmp))
			seq.in				<- as.character(seq.in)
			seq.in[tmp,1]		<- ifelse(seq.in[tmp,1]=='t','c',ifelse(seq.in[tmp,1]=='c','t',ifelse(seq.in[tmp,1]=='a','g','a')))
			seq.in				<- as.DNAbin(seq.in)
			tmp					<- seq.unique( seq.in[ seq.select, ] )
			if( length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )	stop("Unexpected duplicate for 'in'")
			seq.in				<- tmp
		}			
		seq.out					<- seq.out[seq.select,]
		tmp						<- seq.unique(seq.out)
		if( !length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )
			seq.out				<- tmp
		if( length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )
		{
			tmp					<- setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp) )		#name of sequence in triplet that is identical with one other sequence in triplet
			if(verbose)	cat(paste("\nFound identical triplet sequence for region 'out'", tmp))
			seq.out				<- as.character(seq.out)
			seq.out[tmp,1]		<- ifelse(seq.out[tmp,1]=='t','c',ifelse(seq.out[tmp,1]=='c','t',ifelse(seq.out[tmp,1]=='a','g','a')))
			seq.out				<- as.DNAbin(seq.out)
			tmp					<- seq.unique(seq.out[seq.select,])
			if( length(setdiff(c( df.recomb[,parent1],df.recomb[,parent2],df.recomb[,child] ),rownames(tmp))) )	stop("Unexpected duplicate for 'out'")
			seq.out				<- tmp
		}
		if(any(rownames(seq.in)!=rownames(seq.out)))	stop("Unexpected unequal sequences selected")
		#	reset rownames
		tmp						<- c( paste(c("tparent1","tparent2","tchild"),seq.select[1:3],sep='_'), seq.df[,list(label= paste(region,group,FASTASampleCode,sep='_')), by="FASTASampleCode"][,label] )
		rownames(seq.in)		<- tmp
		rownames(seq.out)		<- tmp
		#	save
		file		<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rIn_",gsub('/',':',insignat),".R",sep='')
		if(verbose) cat(paste("\nsave to ",file))
		save(seq.in, file=file)		
		file		<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rOut_",gsub('/',':',insignat),".R",sep='')
		if(verbose) cat(paste("\nsave to ",file))
		save(seq.out, file=file)
	}
	if(1)
	{
		#
		#	run bootstrap ExaML for region 'in', all boostraps on one processor
		#			
		cmd				<- NULL
		file			<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rIn_examlbs",bs.n,'_',gsub('/',':',insignat),".newick",sep='')		
		if(!resume || !file.exists(file))
		{
			infile.exa	<- paste(infile,"_3seqcheck_id",id,"_rIn",sep='')		
			cmd			<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.exa, gsub('/',':',insignat),gsub('/',':',insignat), bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, outdir=indir, opt.bootstrap.by="nucleotide", resume=1, verbose=1)
		}
		#
		#	run bootstrap ExaML for region 'out', all boostraps on one processor
		#						
		file			<- paste(indir,'/',infile,"_3seqcheck_id",id,"_rOut_examlbs",bs.n,'_',gsub('/',':',insignat),".newick",sep='')
		if(!resume || !file.exists(file))
		{
			infile.exa	<- paste(infile,"_3seqcheck_id",id,"_rOut",sep='')		
			cmd			<- c(cmd, hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.exa, gsub('/',':',insignat),gsub('/',':',insignat), bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, outdir=indir, opt.bootstrap.by="nucleotide", resume=1, verbose=1))
		}
		#
		if(verbose) cat(paste("\ncreated ExaML bootstrap runs, n=",length(cmd)))
		if(!is.null(cmd))
		{
			cmd			<- paste(cmd,collapse='\n')
			#cmd		<- paste(cmd,hivc.cmd.recombination.plot.incongruence(indir, infile, gsub('/',':',insignat), triplet.id=id, verbose=1),sep='')				
			#cat(cmd)
			if(verbose) cat(paste("\nqsub ExaML bootstrap runs, hpc.walltime=",hpc.walltime," hpc.mem=",hpc.mem," hpc.nproc=",hpc.nproc," hpc.q=",hpc.q))
			cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem, hpc.nproc=hpc.nproc)
			signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
			outdir		<- paste(DATA,"tmp",sep='/')
			outfile		<- paste("3sc",signat,sep='.')
			#cat(cmd)			
			hivc.cmd.hpccaller(outdir, outfile, cmd)
			Sys.sleep(1)
		}
	}
}		
######################################################################################
hivc.prog.eval.clustering.bias<- function()
{
	indir				<- paste(DATA,"tmp",sep='/')		
	infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
	insignat			<- "Thu_Aug_01_17/05/23_2013"
	indircov			<- paste(DATA,"derived",sep='/')
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infiletree			<- paste(infile,"examlbs100",sep="_")
	
	opt.brl				<- "dist.brl.casc" 
	thresh.brl			<- 0.096
	thresh.bs			<- 0.8
	resume				<- 1
	verbose				<- 1		
	#
	#	load msm clusters
	#
	argv			<<- hivc.cmd.clustering.msm(indir, infiletree, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=resume)
	argv			<<- unlist(strsplit(argv,' '))		
	msm				<- hivc.prog.get.clustering.MSM()	
	#
	df.all			<- msm$df.cluinfo
	file.out.name	<- paste( infiletree,"_clust_",opt.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,"_msmexpgr", sep='')
	#
	#	new diagnoses by CD4
	#
	df.newdiag			<- copy(subset(df.all, select=c(Patient,AnyPos_T1, CD4_T1)))
	setkey(df.newdiag,Patient)
	df.newdiag			<- unique(df.newdiag)
	msmclu.newdiagCD4 	<- hivc.db.getplot.newdiagnosesbyCD4(	df.newdiag, 
			plot.file= paste(dir.name,"/tmp/",file.out.name,"_NewDiagByCD4_",gsub('/',':',insignat),".pdf",sep=''),
			plot.file.p= paste(dir.name,"/tmp/",file.out.name,"_NewDiagByCD4_prop_",gsub('/',':',insignat),".pdf",sep=''),
			plot.ylab= "New diagnoses HIV-1 subtype B,\nin MSM cluster")
	#
	#	seem in care by risk group
	#
	df.living			<- copy(subset(df.all, select=c(Patient, AnyPos_T1, Trm, DateDied, DateLastContact)))
	setkey(df.living,Patient)
	df.living			<- unique(df.living)		
	msmclu.livExpGr 	<- hivc.db.getplot.livingbyexposure(	df.living, 
			plot.file=paste(dir.name,"/tmp/",file.out.name,"_Seen4CareByExpGroup_",gsub('/',':',insignat),".pdf",sep=''),
			plot.file.p=paste(dir.name,"/tmp/",file.out.name,"_Seen4CareByExpGroup_prop_",gsub('/',':',insignat),".pdf",sep=''),
			plot.ylab="Seen for care with HIV-1 subtype B,\nin MSM cluster", db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)
	#
	#	seen in care by recent, untreated/CD4, treated
	#
	file.immu			<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
	load(file.immu)
	df.immu				<- df
	df.immu				<- subset(df.immu, select=c(Patient,PosCD4, CD4) )
	set(df.immu, NULL, "PosCD4", hivc.db.Date2numeric(df.immu[,PosCD4]))
	
	df.living			<- copy(subset(df.all, (is.na(df.all[,Trm]) | Trm%in%c("MSM","BI","HET")) & Sex!='F', select=c(Patient, AnyPos_T1, CD4_T1, isAcute, AnyT_T1, DateDied, DateLastContact)))
	setkey(df.living,Patient)
	df.living			<- unique(df.living)
	msmclu.livCD4		<- hivc.db.getplot.livingbyCD4(	df.living, df.immu, 
			plot.file=paste(dir.name,"/tmp/",file.out.name,"_MSMSeen4CareByCD4_",gsub('/',':',insignat),".pdf",sep=''),
			plot.file.p=paste(dir.name,"/tmp/",file.out.name,"_MSMSeen4CareByCD4_prop_",gsub('/',':',insignat),".pdf",sep=''), 
			plot.ylab="MSM seen for care with HIV-1 subtype B,\nin MSM cluster", db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)
	#
	#
	#	load ATHENA_03 data, subset to MSM and BI
	#
	#
	file			<- paste(dir.name,"/derived/",infilecov,".R",sep='')
	load(file)		
	#
	#	seem in care by risk group
	#		
	df.newdiag		<- copy(subset(df.all, (is.na(df.all[,Trm]) | Trm%in%c("MSM","BI")) & Sex!='F', select=c(Patient,AnyPos_T1, CD4_T1)))
	setkey(df.newdiag,Patient)
	df.newdiag		<- unique(df.newdiag)
	msm.newdiagCD4 	<- hivc.db.getplot.newdiagnosesbyCD4(	df.newdiag, 
															plot.file= paste(dir.name,"/derived/",infilecov,"_MSMNewDiagByCD4.pdf",sep=''),
															plot.file.p= paste(dir.name,"/derived/",infilecov,"_MSMNewDiagByCD4_prop.pdf",sep=''),
															plot.ylab= "MSM new diagnoses HIV-1 subtype B,\n seq available")
	#
	df.living		<- copy(subset(df.all, (is.na(df.all[,Trm]) | Trm%in%c("MSM","BI")) & Sex!='F', select=c(Patient, AnyPos_T1, Trm, DateDied, DateLastContact)))
	setkey(df.living,Patient)
	df.living		<- unique(df.living)
	msm.livExpGr	<- hivc.db.getplot.livingbyexposure(	df.living, 
															plot.file=paste(dir.name,"/derived/",infilecov,"_MSMSeen4CareByExpGroup.pdf",sep=''),
															plot.file.p=paste(dir.name,"/derived/",infilecov,"_MSMSeen4CareByExpGroup_prop.pdf",sep=''),
															plot.ylab="MSM seen for care with HIV-1 subtype B,\n seq available", db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)
	#
	#	seen in care by recent, untreated/CD4, treated
	#
	file.immu		<- paste(dir.name,"derived/ATHENA_2013_03_Immu.R",sep='/')
	load(file.immu)
	df.immu			<- df
	df.immu			<- subset(df.immu, select=c(Patient,PosCD4, CD4) )
	set(df.immu, NULL, "PosCD4", hivc.db.Date2numeric(df.immu[,PosCD4]))
	
	df.living		<- copy(subset(df.all, (is.na(df.all[,Trm]) | Trm%in%c("MSM","BI")) & Sex!='F', select=c(Patient, AnyPos_T1, CD4_T1, isAcute, AnyT_T1, DateDied, DateLastContact)))
	setkey(df.living,Patient)
	df.living		<- unique(df.living)	
	msm.livCD4		<- hivc.db.getplot.livingbyCD4(			df.living, df.immu, 
															plot.file=paste(dir.name,"/derived/",infilecov,"_MSMSeen4CareByCD4.pdf",sep=''),
															plot.file.p=paste(dir.name,"/derived/",infilecov,"_MSMSeen4CareByCD4_prop.pdf",sep=''),
															plot.ylab="MSM seen for care with HIV-1 subtype B,\n seq available", db.endtime=2013.3, db.diff.lastcontact2died=0.5, db.diff.lastcontact2now= 2.3, verbose=1)						
	#	new diagnoses by year	
	clu		<- msmclu.newdiagCD4$t.newdiag
	ds		<- msm.newdiagCD4$t.newdiag
	tmp		<- apply(clu, 2, sum)[as.character(1985:2012)] / apply(ds, 2, sum)[as.character(1985:2012)]
	tmp		<- round(tmp, d=2)
	#
	clu		<- msmclu.newdiagCD4$p.newdiag
	ds		<- msm.newdiagCD4$p.newdiag	
	tmp		<- round( clu[, as.character(1999:2012)] - ds[, as.character(1999:2012)], d=2)
	apply(tmp,1,mean)
	#	seen in care
	clu		<- msmclu.livCD4$p.living
	ds		<- msm.livCD4$p.living
	tmp		<- round( clu[, as.character(1996:2012)] - ds[, as.character(1996:2012)], d=2)
	
	#
}			
######################################################################################
hivc.prog.get.clustering<- function()
{
	library(ape)
	library(data.table)	
			
	opt.dist.brl	<- "dist.brl.casc"
	thresh.bs		<- 0.8
	thresh.brl		<- 0.096	#with mut rate 6e-3/yr expect 2*0.048 brl for 8 yrs apart	
	verbose			<- 1
	resume			<- 0
	indir			<- paste(DATA,"tmp",sep='/')
	infile			<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
	insignat		<- "Sat_Jun_16_17/23/46_2013"

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
						{	switch(substr(arg,2,10),
									thresh.bs= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.bs<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									thresh.brl= return(as.numeric(substr(arg,13,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.brl<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									opt.brl= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.dist.brl<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}	
	
	if(verbose)
	{
		cat(paste("\nopt.dist.brl",opt.dist.brl))
		cat(paste("\nthresh.bs",thresh.bs))
		cat(paste("\nthresh.brl",thresh.brl))
		cat(paste("\nindir",indir))
		cat(paste("\ninfile",infile))
		cat(paste("\ninsignat",insignat))				
	}
	outdir			<- indir
	outfile			<- paste(infile,"_clust_",opt.dist.brl,"_bs",thresh.bs*100,"_brl",thresh.brl*100,sep='')
	outsignat		<- insignat	
	
	if(resume)
	{
		file		<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')		
		options(show.error.messages = FALSE)		
		if(verbose)
			cat(paste("\ntry to resume file ",file))
		readAttempt<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)
			cat(paste("\nresumed file ",file))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{		
		#
		#	load preclustlink files
		#
		file		<- paste(indir,paste(infile,"_preclust_",gsub('/',':',insignat),".R",sep=''),sep='/')
		if(verbose) cat(paste("load file",file))
		options(show.error.messages = FALSE)
		readAttempt<-	try(suppressWarnings(load(file)))
		if(inherits(readAttempt, "try-error"))	stop(paste("cannot load required file", file))
		options(show.error.messages = TRUE)
		#
		#	generate clustering
		#
		dist.brl		<- switch(	opt.dist.brl, 
									"dist.brl.max"		= dist.brl.max,
									"dist.brl.med"		= dist.brl.med,
									"dist.brl.casc"		= dist.brl.casc,
									NA)
		if(any(is.na(dist.brl)))	stop("unexpected NA in dist.brl")					
		clustering		<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
		#
		#	add cluster membership to df.seqinfo
		#
		setkey(df.seqinfo, Node)
		df.seqinfo[,"cluster":= clustering[["clu.mem"]][seq_len(Ntip(ph))]]
		#	set CountryInfection for non-ATHENA seqs
		tmp				<- which( df.seqinfo[,substr(FASTASampleCode,1,2)=="TN"] )
		set(df.seqinfo, tmp, "CountryInfection", "FRGNTN")
		tmp				<- which( df.seqinfo[,substr(FASTASampleCode,1,8)=="PROT+P51"] )
		set(df.seqinfo, tmp, "CountryInfection", "FRGN")
		#
		#	compute TP clusters
		#
		clusters.tp			<- hivc.clu.truepos(clustering, ph.linked, Ntip(ph), verbose= 0)								
		missed.df			<- copy( df.seqinfo )
		plotfile			<- paste(outdir,paste(outfile,"_missedtp_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
		hivc.clu.plot.withinpatientseq.not.samecluster(ph, clustering, clusters.tp, missed.df, plotfile)
		#
		#	compute TN clusters
		#		
		clusters.tn		<- hivc.clu.trueneg(clustering, ph.unlinked.info, ph.unlinked, Ntip(ph), verbose=0)
		#
		#	save
		#		
		file			<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
		cat(paste("\nwrite cluster info to file",file))
		save(df.seqinfo, clustering, clusters.tp, clusters.tn, file=file)
	}
	
	ans	<- list(thresh.bs=thresh.bs, thresh.brl=thresh.brl, opt.dist.brl=opt.dist.brl, df.seqinfo=df.seqinfo, clustering=clustering, clusters.tp=clusters.tp, clusters.tn=clusters.tn)
	ans
}
######################################################################################
hivc.prog.get.clustering.precompute<- function()
{
	library(ape)
	#library(adephylo)
	library(data.table)	
	
	verbose				<- 1
	resume				<- 0
	use.seroneg.as.is	<- 1	#use with updated "ATHENA_2013_03_AllSeqPatientCovariates"
	indir				<- paste(DATA,"tmp",sep='/')	
	infile				<- "ATHENA_2013_03_CurAll+LANL_Sequences_examlbs100"
	insignat			<- "Sat_Jun_16_17/23/46_2013"
	indircov			<- paste(DATA,"derived",sep='/')	
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"	
	
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
		if(length(tmp)>0) infilecov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}	
	outdir			<- indir
	outsignat		<- insignat
	outfile			<- paste(infile,"preclust",sep='_')
	outfile.dtips	<- paste(infile,"predtips",sep='_')
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infilecov)
		print(resume)
	}
	if(resume)
	{
		file		<- paste(outdir,paste(outfile,'_',gsub('/',':',outsignat),".R",sep=''),sep='/')
		options(show.error.messages = FALSE)		
		if(verbose)
			cat(paste("\ntry to resume file ",file))
		readAttempt<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)
			cat(paste("\nresumed file ",file))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{	
		file							<- paste(indir,paste(infile,'_',gsub('/',':',insignat),".newick",sep=''),sep='/')
		if(verbose) cat(paste("\nread phylo from file",file))		
		ph								<- ladderize( read.tree(file) )
		gc()
		#
		#	easy: extract bootstrap support
		#
		ph$node.label[2]				<- 0								#little hack so that clustering works
		ph.node.bs						<- as.numeric( ph$node.label )
		ph.node.bs[is.na(ph.node.bs)]	<- 0
		ph.node.bs						<- ph.node.bs/100
		ph$node.label					<- ph.node.bs
		if(1)
		{
			#
			#	memory consuming: extract branch length statistic of subtree
			#
			if(verbose) cat("\ncompute dist.brl.med")
			dist.brl.med					<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=median)
			gc()
			if(verbose) cat("\ncompute dist.brl.max")
			dist.brl.max					<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=max)		#read patristic distances -- this is the expensive step but still does not take very long
			gc()
			if(verbose) cat("\ncompute dist.brl.casc")
			dist.brl.casc					<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
			gc()		
		}
		#
		#	easy: extract tree specific TP and FN data sets
		#		
		file							<- paste(indircov,"/",infilecov,".R",sep='')
		load(file)
		if(verbose) cat(paste("\nfound covariates for patients, n=",nrow(df.all)))
		df.seqinfo						<- subset(df.all, !is.na(PosSeqT) )
		if(verbose) cat(paste("\nfound covariates for patients with non-NA PosSeqT, n=",nrow(df.seqinfo)))
		if(verbose) cat(paste("\nstart: compute TP and TN data tables for phylogeny"))
		tmp								<- hivc.phy.get.TP.and.TN(ph, df.seqinfo, verbose=verbose, use.seroneg.as.is= use.seroneg.as.is)		
		if(verbose) cat(paste("\nend: compute TP and TN data tables for phylogeny"))
		unlinked.byspace				<- tmp[["unlinked.byspace"]]
		unlinked.bytime					<- tmp[["unlinked.bytime"]]
		linked.bypatient				<- tmp[["linked.bypatient"]]	
		ph.linked						<- tmp[["ph.linked"]]
		ph.unlinked.info				<- tmp[["ph.unlinked.info"]]
		ph.unlinked						<- tmp[["ph.unlinked"]]
		#
		#	add node number to df.seqinfo
		#
		df.seqinfo						<- merge( df.all, data.table(Node=seq_len(Ntip(ph)), FASTASampleCode=ph$tip.label), all.y=1, by="FASTASampleCode")
		#
		#	get and plot distribution of bootstrap values for TP and TN pairs
		#
		ph.mrca				<- mrca(ph)
		#	prepare data.tables with mrca: dead/seroneg TN pairs		
		bs.unlinkedpairs	<- lapply(unlinked.bytime, function(x)
								{						
									set(x, NULL, "query.FASTASampleCode", as.character(x[,query.FASTASampleCode]))
									set(x, NULL, "FASTASampleCode", as.character(x[,FASTASampleCode]))
									ans					<- merge(x, x[, list(mrca= ph.mrca[query.FASTASampleCode,FASTASampleCode]) ,by=FASTASampleCode],by="FASTASampleCode")
									setnames(ans, c("FASTASampleCode","query.FASTASampleCode"), c("tip2","tip1"))
									subset(ans, select=c(tip1, tip2, mrca))
								})
		bs.unlinkedpairs	<- rbindlist(bs.unlinkedpairs)
		#	prepare data.tables with mrca: geographically distant seqs and any indb seq
		unlinked.byspace[,dummy:=seq_len(nrow(unlinked.byspace))]
		set(unlinked.byspace, NULL, "FASTASampleCode", as.character(unlinked.byspace[,FASTASampleCode]))	
		seq.indb			<- colnames(ph.mrca)[ which( substr(colnames(ph.mrca),1,2)!="TN" & substr(colnames(ph.mrca),1,8)!="PROT+P51" ) ]
		if(nrow(unlinked.byspace))
			bs.unlinked.byspace	<- unlinked.byspace[,	list(tip1=FASTASampleCode,tip2=seq.indb, mrca= ph.mrca[seq.indb,FASTASampleCode]), by="dummy"]
		else
			bs.unlinked.byspace	<- NULL
		
		#	prepare data.tables with mrca: within patient TP pairs
		setkey(linked.bypatient,Patient)
		bs.linked.bypatient	<- linked.bypatient[, {
														tmp					<- match(FASTASampleCode, ph$tip.label)
														ans					<- t(combn(tmp, 2 ))
														ans					<- cbind(ans, apply(ans, 1, function(z)  ph.mrca[z[1],z[2]]))
														data.table(tip1=ans[,1], tip2=ans[,2], mrca=ans[,3])
													}, by="Patient"]		
		tmp								<- hivc.phy.get.TP.and.TN.bootstrapvalues(ph, bs.linked.bypatient, ph.mrca=ph.mrca, df.seqinfo=df.seqinfo, bs.unlinkedpairs=bs.unlinkedpairs, bs.unlinked.byspace=bs.unlinked.byspace, dist.brl=dist.brl.casc, thresh.brl=0.096, plot.file=paste(outdir,'/',outfile,"_distbs_",gsub('/',':',outsignat),".pdf",sep=''), verbose=verbose)		
		bs.linked.bypatient				<- tmp[["bs.linked.bypatient"]] 
		bs.unlinkedpairs				<- tmp[["bs.unlinkedpairs"]]
		bs.unlinked.byspace				<- tmp[["bs.unlinked.byspace"]]
		#
		#	memory consuming: compute branch length matrix between tips
		#
		#if(verbose) cat("\ncompute dist.root")
		#dist.root						<-  distRoot(ph, method= "patristic")
		#gc()
		#if(verbose) cat("\ncompute dist.tips.mat")
		#dist.tips.mat					<-  distTips(ph, method= "patristic")
		#gc()		
		#
		#	save output
		#
		file							<- paste(outdir,paste(outfile,'_',gsub('/',':',outsignat),".R",sep=''),sep='/')	
		if(verbose)	cat(paste("write to file",file))
		save(	ph, dist.brl.max, dist.brl.med, dist.brl.casc, ph.node.bs, ph.linked, ph.unlinked.info, ph.unlinked, 
				df.seqinfo, unlinked.byspace, unlinked.bytime, linked.bypatient,  
				bs.linked.bypatient, bs.unlinkedpairs, bs.unlinked.byspace, file=file )
		#save(ph, dist.brl.max, dist.brl.med, dist.brl.casc, ph.node.bs, ph.linked, ph.unlinked.info, ph.unlinked, df.seqinfo, unlinked.byspace, unlinked.bytime, linked.bypatient,  file=file )		
	}
	
	ans	<- list(	ph=ph, #dist.tips.mat=dist.tips.mat, #dist.root=dist.root,  
					dist.brl.max=dist.brl.max, dist.brl.med=dist.brl.med, dist.brl.casc=dist.brl.casc, 
					ph.node.bs=ph.node.bs, ph.linked=ph.linked, ph.unlinked.info=ph.unlinked.info, ph.unlinked=ph.unlinked, 
					df.seqinfo=df.seqinfo, unlinked.byspace=unlinked.byspace, unlinked.bytime=unlinked.bytime, linked.bypatient=linked.bypatient,
					bs.linked.bypatient=bs.linked.bypatient, bs.unlinkedpairs=bs.unlinkedpairs, bs.unlinked.byspace
					)
	ans				
}



hivc.prog.remove.resistancemut<- function()
{
	library(ape)
	library(data.table)
	
	#load drug resistance mutations and select unique mutants by codon
	load( paste( CODE.HOME,"/data/IAS_primarydrugresistance_201303.rda",sep='' ) )
	dr				<- as.data.table(IAS_primarydrugresistance_201303)
	
	indir			<- paste(DATA,"tmp",sep='/')
	infile			<- "ATHENA_2013_03_CurAll+LANL_Sequences"
	insignat		<- "Sat_Jun_16_17/23/46_2013"
	outdir			<- paste(DATA,"tmp",sep='/')
	outfile			<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"
	outsignat		<- "Thu_Aug_01_17/05/23_2013"
	alignment.start	<- 2253	
	verbose			<- 1
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
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									outfile= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) outsignat<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,16),
									alignment.start= return(substr(arg,18,nchar(arg))),NA)	}))
		if(length(tmp)>0) alignment.start<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
	}
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(outdir)
		print(outfile)		
		print(outsignat)
		print(alignment.start)
		print(verbose)
	}	
	#load alignment
	file				<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
	load(file)	
	#modify dr table for particular alignment	
	set(dr, NULL, "Alignment.nuc.pos", dr[,Alignment.nuc.pos]-alignment.start+1)
	
	#remove	likely.nonB.outliers
	cat("\nchange infile: remove	likely.nonB.outliers")
	likely.nonB.outliers	<- c("R03-07193","2006G206","PROT+P51_B.AU.1995.C92.AF538307","2008G084")
	likely.nonB.outliers	<- which(rownames(seq.PROT.RT) %in% likely.nonB.outliers)
	seq.PROT.RT				<- seq.PROT.RT[-likely.nonB.outliers,]
	
	#if alignment equals any of the drug resistance mutants, replace with NNN	
	seq.PROT.RT			<- seq.rm.drugresistance(as.character(seq.PROT.RT), dr, verbose=verbose, rtn.DNAbin=1 )	
	
	#save No Drug resistance alignment to file
	file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".R",sep='')
	if(verbose)	cat(paste("\nwrite R file to",file))
	save(seq.PROT.RT, file=file)	
	file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".phylip",sep='')
	if(verbose)	cat(paste("\nwrite phylip file to",file))
	seq.write.dna.phylip(seq.PROT.RT, file=file)					
	file								<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".fasta",sep='')			
	if(verbose)	cat(paste("\nwrite fasta file to",file))
	write.dna(seq.PROT.RT, file=file, format="fasta", colsep='', colw=ncol(seq.PROT.RT), blocksep=0)
	
	seq.PROT.RT
}
######################################################################################
hivc.prog.BEAST.evalpoolrun<- function()
{	
	require(phangorn)
	indircov			<- paste(DATA,"derived",sep='/')	
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"			
	file.cov			<- paste(indircov,"/",infilecov,".R",sep='')
	file.viro			<- paste(indircov,"/ATHENA_2013_03_Viro.R",sep='/')
	file.immu			<- paste(indircov,"/ATHENA_2013_03_Immu.R",sep='/')
	file.treatment		<- paste(indircov,"/ATHENA_2013_03_Regimens.R",sep='/')
		
	indir				<- paste(DATA,"beast/beast_131011",sep='/')
	#indir				<- paste(DATA,"tmp",sep='/')	
	infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_beast_seroneg"
	insignat			<- "Tue_Aug_26_09/13/47_2013"
	infilexml.opt		<- "txs4clu"
	infilexml.opt		<- "mph4clu"
	infilexml.opt		<- "mph4clutx4tip"
	infilexml.template	<- "um232rhU2045"
	infilexml.template	<- "um182rhU2045"
	#infilexml.opt		<- "mph4cluLdTd"
	#infilexml.template	<- "um22rhG202018"
	#infilexml.template	<- "um182rhU2045ay"	
	
	plot						<- 1
	verbose						<- 1
	resume						<- 0
	pool.n						<- 3	
	beastlabel.idx.clu			<- 1
	beastlabel.idx.hivs			<- ifelse(any(sapply(c("LdTd","LsTd"), function(x) grepl(x,infilexml.opt))),	5,	4)
	beastlabel.idx.samplecode	<- 6	
	beastlabel.idx.rate			<- 7
	
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
		#		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,14),
									infilexml.opt= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.opt<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									infilexml.template= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.template<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									pool.n= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) pool.n<- tmp[1]
	}	
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(infilexml.opt)
		print(infilexml.template)
		print(pool.n)
	}
	if(resume)
	{
		file		<- paste(indir,'/',infile,'_',infilexml.template,'_',infilexml.opt,"_mcc_",gsub('/',':',insignat),".R",sep='')		
		options(show.error.messages = FALSE)		
		if(verbose)
			cat(paste("\ntry to resume file ",file))
		readAttempt<-	try(suppressWarnings(load(file)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)
			cat(paste("\nresumed file ",file))
		return(list(cluphy=cluphy, cluphy.trmca=cluphy.trmca, cluphy.df=cluphy.df))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{
		#	
		#	read annotated mcc trees
		#
		tmp			<- list.files(indir, pattern=paste(".nex$",sep=''))
		file.nex	<- tmp[ grepl(paste('_',infilexml.template,'_',sep=''), tmp) & grepl(paste('_',infilexml.opt,'_',sep=''), tmp) & grepl(gsub('/',':',insignat), tmp) ]
		if(verbose)	cat(paste("\nFound .nex files matching input args, n=", length(file.nex)))
		tmp			<- list.files(indir, pattern=paste(".log$",sep=''))
		file.log	<- tmp[ grepl(paste('_',infilexml.template,'_',sep=''), tmp) & grepl(paste('_',infilexml.opt,'_',sep=''), tmp) & grepl(gsub('/',':',insignat), tmp) ]
		if(verbose)	cat(paste("\nFound .log files matching input args, n=", length(file.log)))		
		#		
		if(length(file.nex)==pool.n)
		{
			#	read treeannotator .nex file
			ph.beast		<- lapply(file.nex, function(x)		hivc.treeannotator.read(paste(indir,x,sep='/'), add.to.tiplabel=c("rate_median"), rate.multiplier=1e3, round.digit=2, verbose=verbose)		)
			if(verbose)	cat(paste("\nRead trees matching input args, n=", length(ph.beast)))
			#	read length of tip stems
			file.log		<- sapply(file.log, function(x)					paste(indir,x,sep='/') )
			file.xml		<- sapply(file.log, function(x)					paste(substr(x,1,nchar(x)-3), "xml", sep='') )			
			df.tstem		<- lapply(seq_along(file.log), function(i)		hivc.beast.read.log2tstem(file.log[i], file.xml[i], beastlabel.idx.samplecode=6, burn.in= 5e6, breaks.n= 30, verbose=1)		)
			df.tstem		<- rbindlist( df.tstem )
			#	load all patient covariates
			load(file.cov)								
			#	build single tree from pooled runs
			tmp				<- hivc.treeannotator.get.phy(ph.beast, beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate, debug=0)
			cluphy			<- tmp$cluphy 
			ph.tip.ctime	<- tmp$ph.tip.ctime 				
			ph.root.ctime	<- tmp$ph.root.ctime
			#	extract rates		
			rates.df		<- hivc.treeannotator.get.rates(cluphy, tmp$ph.tip.df, nodelabel.idx.edgewidth=5)
			#	convert tstem time into calendar time
			tmp				<- data.table( FASTASampleCode=cluphy$tip.label, tip=seq_along(cluphy$tip.label), mrca= Ancestors(cluphy, seq_along(cluphy$tip.label), type="parent")-Ntip(cluphy) )			 
			df.tstem		<- merge( df.tstem, tmp, by="FASTASampleCode" )			
			set(df.tstem, NULL, "tstem", df.tstem[, max(ph.tip.ctime)-tstem])
			df.tstem		<- subset(df.tstem, select=c(tip, mrca, tstem, density))
			#	get length of 95% TMRCAs of tip stems
			cluphy.tstem	<- df.tstem[,	{
												x<- data.table(tstem, density)
												x[,dummy:=-density]
												setkey(x,dummy)			
												x[,cdensity:= cumsum(x[,density]/sum(x[,density]))]							
												list(height_95_diff= diff(range(subset(x, cdensity<0.95, tstem))) )							
											},	by="tip"]			
			#
			cluphy.df		<- hivc.treeannotator.get.clusterprob(ph.beast, beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.samplecode=beastlabel.idx.samplecode)
			cluphy.df		<- merge(cluphy.df, df.all, by="FASTASampleCode")
			cluphy.tmrca	<- hivc.treeannotator.get.tmrcas(ph.beast, beastlabel.idx.hivs=beastlabel.idx.hivs) 						
			#
			#	save mcc tree in R format
			#
			file		<- paste(indir,'/',infile,'_',infilexml.template,'_',infilexml.opt,"_mcc_",gsub('/',':',insignat),".R",sep='')
			if(verbose)	cat(paste("\nSave mcc tree objects to file",file))
			save(cluphy, cluphy.tmrca, cluphy.df, cluphy.tstem, file=file)			
			#
			#	plot cluster TMRCAs
			#
			if(plot)
			{
				tmp			<- sort(cluphy.tmrca[, height_95_diff])	
				file		<- paste(indir,'/',infile,'_',infilexml.template,'_',infilexml.opt,"_mcc_clutmrca_",gsub('/',':',insignat),".pdf",sep='')
				if(verbose)	cat(paste("\nPlot cluster TMRCAs to file",file))
				pdf(file,width=7,height=7)
				plot( tmp, seq_len(nrow(cluphy.tmrca)), type="l", xlab='width of 95% interval of TMRCA of clusters', ylab='#MRCA of clusters', bty='n' )
				tmp			<- quantile(tmp, prob=c(0.25,0.5,0.8,0.9,0.95))
				legend("bottomright",bty='n',border=NA,legend=paste( apply(rbind(names(tmp), round(tmp,d=2)),2,function(x) paste(x,collapse='=',sep='')), collapse=', ',sep=''))
				legend("topleft",bty='n',border=NA,legend=paste(infilexml.template,infilexml.opt,sep='_'))
				dev.off()
			}
			#
			#	plot tip stem TMRCAs
			#
			if(plot)
			{
				tmp			<- sort(cluphy.tstem[, height_95_diff])	
				file		<- paste(indir,'/',infile,'_',infilexml.template,'_',infilexml.opt,"_mcc_tiptmrca_",gsub('/',':',insignat),".pdf",sep='')
				if(verbose)	cat(paste("\nPlot tip TMRCAs to file",file))
				pdf(file,width=7,height=7)
				plot( tmp, seq_len(nrow(cluphy.tstem)), type="l", xlab='width of 95% interval of TMRCA of tips', ylab='#MRCA of tips', bty='n' )
				tmp			<- quantile(tmp, prob=c(0.25,0.5,0.8,0.9,0.95))
				legend("bottomright",bty='n',border=NA,legend=paste( apply(rbind(names(tmp), round(tmp,d=2)),2,function(x) paste(x,collapse='=',sep='')), collapse=', ',sep=''))
				legend("topleft",bty='n',border=NA,legend=paste(infilexml.template,infilexml.opt,sep='_'))
				dev.off()
			}
			#
			#	plot clusters
			#
			if(plot)
			{
				#load patient RNA
				load(file.viro)
				df.viro				<- df				
				#load patient CD4				
				load(file.immu)
				df.immu				<- df
				#load patient regimen
				load(file.treatment)
				df.treatment		<- df
				
				youngest.tip.ctime	<- max(ph.tip.ctime)
				#youngest.tip.ctime	<- 2010.46
				file				<- paste(indir,'/',infile,'_',infilexml.template,'_',infilexml.opt,"_mcc_",gsub('/',':',insignat),".pdf",sep='')
				if(verbose)	cat(paste("\nplotting dated clusters to file", file ))
				dummy				<- hivc.treeannotator.plot(cluphy, ph.root.ctime, youngest.tip.ctime, df.all, df.viro, df.immu, df.treatment=df.treatment, df.tstem=df.tstem, df.rates=rates.df, end.ctime=2013.3, cex.nodelabel=0.5, cex.tiplabel=0.5, file=file, pdf.width=7, pdf.height=150)
			}
		}
		else if(verbose)	
			cat(paste("\nNothing evaluated - waiting for BEAST runs to complete"))
	}
}
######################################################################################
hivc.prog.BEAST2.plot.cluster.trees<- function()
{
	require(ape)
	require(data.table)
	require(hivclust)
	#	program options
	opt.topo.posterior.cprob	<- 0.75
	beastlabel.idx.clu			<- 1
	beastlabel.idx.hivn			<- 2
	beastlabel.idx.hivd			<- 3
	beastlabel.idx.hivs			<- 4
	beastlabel.idx.samplecode	<- 6
	beastlabel.idx.rate			<- NA
	#	input files
	indir					<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/beast/beast2_140201'
	infile					<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_seroneg-130"
	insignat				<- "Tue_Aug_26_09:13:47_2013"
	infilexml.opt			<- "rsu815"
	infilexml.template		<- "sasky_sdr06"
	indircov				<- paste(DATA,"derived",sep='/')	
	infilecov				<- "ATHENA_2013_03_AllSeqPatientCovariates"			
	#
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
						{	switch(substr(arg,2,14),
									infilexml.opt= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.opt<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									infilexml.template= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.template<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indircov= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indircov<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infilecov= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilecov<- tmp[1]		
		#
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) outsignat<- tmp[1]				
		#		
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}
	#
	file.cov				<- paste(indircov,"/",infilecov,".R",sep='')
	file.viro				<- paste(indircov,"/ATHENA_2013_03_Viro.R",sep='/')
	file.immu				<- paste(indircov,"/ATHENA_2013_03_Immu.R",sep='/')
	file.treatment			<- paste(indircov,"/ATHENA_2013_03_Regimens.R",sep='/')
	outdir					<- indir
	outfile					<- paste(infile,infilexml.template,infilexml.opt, sep='_')
	outsignat				<- insignat	
	#
	if(verbose)
	{
		print(resume)
		print(indir)		
		print(infile)
		print(insignat)
		print(infilexml.opt)
		print(infilexml.template)
		print(indircov)
		print(infilecov)
		print(outdir)
		print(outfile)
		print(outsignat)
	}
	#	load patient demographic data	(loads df.all)
	load(file.cov)		
	#	load patient RNA
	load(file.viro)
	df.viro				<- df				
	#	load patient CD4				
	load(file.immu)
	df.immu				<- df
	#	load patient regimen
	load(file.treatment)
	df.treatment		<- df		
	#	collect files containing dated cluster phylogenies
	files		<- list.files(indir)
	files		<- files[ sapply(files, function(x) grepl(infile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(infilexml.opt, x, fixed=1) & grepl(infilexml.template,x, fixed=1) & grepl('_cluposterior_[0-9]+',x) & grepl('R$',x) ) ]		
	if(!length(files))	stop('no input files matching criteria')
	tmp			<- regmatches( files, regexpr('_cluposterior_[0-9]+',files)) 
	cluster		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
	file.info	<- data.table(file=files, cluster=cluster)
	setkey(file.info, cluster)
	#
	#	for each cluster
	#	plot dated cluster phylogenies that together have 'opt.topo.posterior.cprob' posterior prob with clinical data points 
	#
	dummy<- lapply(seq_len(nrow(file.info)), function(i)
			{
				#	load dated cluster phylogenies
				file			<- paste(indir, file.info[i,file], sep='/')
				cat(paste('\nload file',file))
				tmp				<- load(file)
				#	start: remove in updated saves		
				#set(mph.node.ctime,NULL,'q',round(mph.node.ctime[,q],d=3))
				#mph.node.ctime	<- mph.node.ctime[,	{
				#			tmp<- setdiff(seq_along(q),which(duplicated(q))-1)
				#			if(length(tmp)==1)	
				#				ans	<- list(cluster= rep(cluster[tmp],2), q=q[tmp]+c(-0.001,0), cdf=c(0.,1.))
				#			else
				#				ans	<- list(cluster= cluster[tmp], q=q[tmp], cdf=cdf[tmp])
				#			ans
				#		},by=c('equal.to','node')]
				#	end:remove in updated saves	
				#	get pdf from cdf of node calendar times
				#	mph.node.ctime	<- mph.node.ctime[, list(cluster=cluster, q=q, pdf=c(0,diff(cdf)/diff(q)), cdf=cdf),by=c('equal.to','node')]
				#
				topo.info		<- lapply( strsplit(names(ph.consensus),'_'), function(x)		data.table(cluster=x[1], mph.i=x[2], prob=x[3])	)
				topo.info		<- do.call('rbind',topo.info)
				topo.info[,ph.i:=seq_along(ph.consensus)]
				set(topo.info, NULL, 'cluster', as.numeric( regmatches(topo.info[, cluster], regexpr('[0-9]+',topo.info[, cluster])) ))
				set(topo.info, NULL, 'mph.i', as.numeric( regmatches(topo.info[, mph.i], regexpr('[0-9]+',topo.info[, mph.i])) ))
				set(topo.info, NULL, 'prob', as.numeric( regmatches(topo.info[, prob], regexpr('[0-1].?[0-9]*',topo.info[, prob])) ))
				topo.info[, sort:= -prob]
				setkey(topo.info, sort)		
				topo.info[, cprob:= cumsum( topo.info[, prob] )]						
				#	select a subset of topologies for plotting
				if(nrow(topo.info)>5)
					topo.info	<- rbind( subset( topo.info, cprob<opt.topo.posterior.cprob ), topo.info[which(cprob>=opt.topo.posterior.cprob)[1], ])
				#	plot all topologies into a single pdf
				file				<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),'_cluposterior+clinical_',topo.info[1,cluster],".pdf",sep='')
				cat(paste('\nplot dated cluster phylogenies with clinical data to file',file))
				pdf(file=file, w=5, h=5)				
				dummy		<- sapply( seq_len(nrow(topo.info)), function(topo.i)
						{
							#
							#	for each dated cluster phylogeny
							#
							cat(paste('\nplot topology',topo.i))
							cluphy				<- ph.consensus[[ topo.info[topo.i, ph.i] ]]
							cluphy.prob			<- topo.info[topo.i, prob]
							cluphy.info			<- hivc.treeannotator.tiplabel2df(cluphy, beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
							cluphy.info[, tip:=match( cluphy.info[,BEASTlabel],cluphy$tip.label )]
							cluphy$tip.label	<- cluphy.info[, FASTASampleCode]							
							cluphy.node.ctime	<- subset(mph.node.ctime, equal.to==topo.info[topo.i, mph.i] & node!=0)
							cluphy.root.ctime	<- cluphy.info[1,TipT]-node.depth.edgelength(cluphy)[1] 							
							cluphy.tip.ctime	<- cluphy.info[, TipT]											
							dummy				<- hivc.beast2out.plot.cluster.trees(df.all, df.immu, df.viro, df.treatment, cluphy, cluphy.root.ctime, cluphy.tip.ctime, ph.prob=cluphy.prob, df.node.ctime=copy(cluphy.node.ctime), df.rates=NULL, end.ctime=2013.3,  cex.nodelabel=0.5,  cex.tiplabel=0.5,  file=NULL,  pdf.width=7, pdf.height=20)					
						})
				dev.off()
			})
}
######################################################################################
hivc.prog.BEAST2.process.cluster.trees<- function()
{
	require(ape)
	require(data.table)
	require(hivclust)
	#	program options
	resume						<- TRUE
	verbose						<- TRUE
	clu							<- NA
	opt.max.runs.per.cluster	<- 4
	beastlabel.idx.clu			<- 1
	beastlabel.idx.hivn			<- 2
	beastlabel.idx.hivd			<- 3
	beastlabel.idx.hivs			<- 4
	beastlabel.idx.samplecode	<- 6
	beastlabel.idx.rate			<- NA
	cdf.by						<- 0.025
	#	input files		
	indir					<- paste(DATA,"tmp",sep='/')
	indir					<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/beast/beast2_140201'
	infile					<- "ATHENA_2013_03_NoDRAll+LANL_Sequences_seroneg-130"
	insignat				<- "Tue_Aug_26_09:13:47_2013"
	infilexml.opt			<- "rsu815"
	infilexml.template		<- "sasky_sdr06"
	#
	indir					<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/sasky_sdr06_-DR-RC-SH+LANL_alrh160_pool'
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	insignat				<- "Wed_Dec_18_11:37:00_2013"
	infilexml.opt			<- "alrh160"
	infilexml.template		<- "sasky_sdr06fr"
	#
	if(0)
	{
		indir					<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/beast1pool"		
		infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"		
		insignat				<- "Wed_Dec_18_11:37:00_2013"				
		infilexml.template		<- "um192rhU2080"
		infilexml.opt			<- "mph4clutx4tip"			
	}
	#
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
						{	switch(substr(arg,2,14),
									infilexml.opt= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.opt<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									infilexml.template= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.template<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									cluster= return(as.numeric(substr(arg,10,nchar(arg)))),NA)	}))
		if(length(tmp)>0) clu<- tmp[1]
		#
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) outsignat<- tmp[1]				
		#		
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}	
	#
	outdir					<- indir
	outfile					<- paste(infile,infilexml.template,infilexml.opt, sep='_')
	outsignat				<- insignat
	#
	if(verbose)
	{
		print(resume)
		print(indir)		
		print(infile)
		print(insignat)
		print(infilexml.opt)
		print(infilexml.template)
		print(outdir)
		print(outfile)
		print(outsignat)
		print(clu)
	}
	#	pool overlapping cluster trees  
	files		<- list.files(indir)
	files		<- files[ sapply(files, function(x) grepl(infile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(paste('_',infilexml.opt,'_',sep=''), x, fixed=1) & grepl(paste('_',infilexml.template,'_',sep=''), x, fixed=1) & grepl('_pool_[0-9]+',x) & grepl('_clu_[0-9]+',x) & grepl('R$',x) ) ]
	#
	tmp			<- regexpr('_clu_[0-9]+',files)
	if(any(tmp<0))	stop('unexpected _clu_ files')
	tmp			<- regmatches( files, tmp) 
	cluster		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
	tmp			<- regmatches( files, regexpr('_pool_[0-9]+',files)) 
	run.id		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
	file.info	<- data.table(file=files, cluster=cluster, run.id=run.id )
	setkey(file.info, cluster)
	if(!is.na(clu))
		file.info	<- subset( file.info, cluster==clu )
	cat(paste('\nnumber of files to process, n=',nrow(file.info)))
	#
	dummy	<- lapply( file.info[, unique(cluster)], function(clu)
			{
				file			<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),'_','clutrees_',clu,'.R',sep='')
				if(resume)
				{
					options(show.error.messages = FALSE)		
					readAttempt		<- try(suppressWarnings(load(file)))
					if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))					
					options(show.error.messages = TRUE)		
				}
				if(!resume || inherits(readAttempt, "try-error"))
				{
					tmp				<- subset(file.info, cluster==clu)[, file]					
					if(length(tmp)>opt.max.runs.per.cluster)	
					{
						cat(paste('\nFound more than',opt.max.runs.per.cluster,'runs for cluster',clu,'Keep only first',opt.max.runs.per.cluster,'of these.'))
						tmp			<- tmp[seq_len(opt.max.runs.per.cluster)]
					}
					mph.clu			<- hivc.beast2out.pool.cluster.trees( paste(indir, tmp, sep='/') )					
					cat(paste("\nsave mph.clu to file",file))
					save(mph.clu, file=file)
				}
			})
	
	#
	files		<- list.files(indir)
	files		<- files[ sapply(files, function(x) grepl(infile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(paste('_',infilexml.opt,'_',sep=''), x, fixed=1) & grepl(paste('_',infilexml.template,'_',sep=''),x, fixed=1) & grepl('_clutrees_[0-9]+',x) & grepl('R$',x) ) ]		
	if(!length(files))	stop('no input files matching criteria')
	tmp			<- regmatches( files, regexpr('_clutrees_[0-9]+',files)) 
	cluster		<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
	file.info	<- data.table(file=files, cluster=cluster)
	setkey(file.info, cluster)
	if(!is.na(clu))
		file.info	<- subset( file.info, cluster==clu )
	cat(paste('\nnumber of files to process, n=',nrow(file.info)))
	#
	dummy		<- lapply(seq_len(nrow(file.info)), function(file.i)
			{			
				clu				<- file.info[file.i,cluster]
				if(resume)
				{
					file	<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),'_cluposterior_',clu,'.R',sep='')
					options(show.error.messages = FALSE)		
					readAttempt		<- try(suppressWarnings(load(file)))
					if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))					
					options(show.error.messages = TRUE)							
				}
				if(!resume || inherits(readAttempt, "try-error"))
				{					
					file			<- paste(indir,'/',file.info[file.i,file],sep='')
					cat(paste("\nload file",file))
					options(show.error.messages = FALSE)		
					readAttempt		<- try(suppressWarnings(load(file)))
					if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))
					if(inherits(readAttempt, "try-error"))	stop(paste("\ncould not read file",file))
					options(show.error.messages = TRUE)		
					#
					#	for each cluster, extract indices of sequentially distinct topologies and their frequency
					#
					cat(paste('\nprocess cluster=',clu,'with #tips=',Ntip(mph.clu[[1]])))
					tmp				<- hivc.beast2.extract.distinct.topologies(mph.clu)
					mph.clu.dtopo	<- tmp$dtopo
					mph.clu.itopo	<- tmp$itopo					
					suppressWarnings( mph.clu.dtopo[, cluster:= clu ] )
					suppressWarnings( mph.clu.itopo[, cluster:= clu ] )
					setkey(mph.clu.dtopo, freq)
					#
					#for each topology, check if tips are in same order and if edges are in same order					
					#
					cat(paste('\ncheck tip, node and edge order'))
					mph.topo.info	<- mph.clu.itopo[,	{
								if(length(mph.i)==1)
								{
									tips.in.order	<- edges.in.order<- nodes.in.order<- TRUE
									tmp				<- mph.i[1]
								}
								else
								{
									tmp				<- mph.clu[[ mph.i[1] ]]$tip.label
									tips.in.order	<- sapply( mph.i[-1], function(i)	identical( tmp, mph.clu[[ i ]]$tip.label) )
									tmp				<- mph.clu[[ mph.i[1] ]]$edge[,2]
									edges.in.order	<- sapply( mph.i[-1], function(i)	identical( tmp, mph.clu[[ i ]]$edge[,2]) )
									tmp				<- mph.clu[[ mph.i[1] ]]$edge[,1]
									nodes.in.order	<- sapply( mph.i[-1], function(i)	identical( tmp, mph.clu[[ i ]]$edge[,1]) )																
									tmp				<- mph.i[-1]
								}
								list(tips.in.order= tips.in.order, edges.in.order=edges.in.order, nodes.in.order=nodes.in.order, mph.i=tmp)
							}, by='equal.to']
					mph.topo.info	<- subset( mph.topo.info, !tips.in.order | !edges.in.order | !nodes.in.order )
					cat(paste('\ntip, node or edges not in order for n=',nrow(mph.topo.info)))
					for(i in seq_len(nrow(mph.topo.info)))
						mph.clu[[ mph.topo.info[i,mph.i] ]]<- hivc.phy.reorder( mph.clu[[ mph.topo.info[i,equal.to] ]], mph.clu[[ mph.topo.info[i,mph.i] ]] )
					#
					#	tips and edges in same order, so very easy to summarize node heights of each topology
					#
					cat(paste('\ncalculate node calendar times'))
					mph.info		<- hivc.treeannotator.tiplabel2df(mph.clu[[1]], beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
					mph.info[, tip:=match( mph.info[, BEASTlabel], mph.clu[[1]]$tip.label)]
					if(1)
					{
						tmp	<- node.depth.edgelength( mph.clu[[ mph.clu.dtopo[1, mph.i] ]] )
						tmp	<- mph.info[1,TipT]-tmp[1] + tmp
						tmp	<- tmp[seq_len(nrow(mph.info))] - mph.info[,TipT]
						if( any( abs(tmp)>EPS ) ) warning('branch lengths do not add up to tip times')
					}
					mph.node.ctime	<- mph.clu.itopo[,	{					
															mph.info		<- hivc.treeannotator.tiplabel2df(mph.clu[[ mph.i[1] ]], beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
															mph.info[, tip:=match( mph.info[, BEASTlabel], mph.clu[[ mph.i[1] ]]$tip.label)]															
															node			<- c(seq.int(from=Ntip(mph.clu[[1]])+1, to=Nnode(mph.clu[[1]], internal=FALSE)),0)
															tmp				<- sapply( mph.i, function(i)
																	{
																		depth	<- c( node.depth.edgelength( mph.clu[[ i ]] ), -mph.clu[[ i ]]$root.edge )
																		tmp		<- which.max(depth)
																		depth	<- depth-depth[tmp]+subset(mph.info, tip==tmp)[, TipT]
																		depth[ seq.int(from=Ntip(mph.clu[[ i ]])+1, to=length(depth)) ]																		
																	})
															tmp				<- apply(tmp, 1, quantile, probs=seq(0,1,by=0.1))		
															list( q= round(as.numeric(tmp),d=3), cdf=seq(0,1,by=0.1), node=rep(node, each=11) )					
														}, by='equal.to']
					mph.node.ctime	<- mph.node.ctime[,	{
															tmp<- setdiff(seq_along(q),which(duplicated(q))-1)
															if(tmp[1]!=1)
																tmp[1]<- 1
															if(length(tmp)==1)	
																ans	<- list(cluster= rep(cluster[tmp],2), q=q[tmp]+c(-0.001,0), cdf=c(0.,1.))
															else
																ans	<- list(cluster= cluster[tmp], q=q[tmp], cdf=cdf[tmp])
															ans
														},by=c('equal.to','node')]
					mph.node.ctime		<- mph.node.ctime[, {
																tmp<- c(0,diff(cdf)/diff(q))
																list(cluster=cluster, q=q, pdf=tmp/sum(tmp), cdf=cdf)
															},by=c('equal.to','node')]					
					mph.node.ctime[, cluster:=clu]
					#
					#	tips and edges in same order, so not difficult to summarize posterior node heights across all topologies
					#
					tmp			<- mph.clu.dtopo[, mph.i[which.max(freq)]]
					tmp			<- as.matrix( mrca( mph.clu[[tmp]], full = FALSE ) )					
					ref.mrca	<- do.call('rbind', lapply( seq_len(ncol(tmp)), function(j)		data.table( tip1=rownames(tmp)[-j], tip2=colnames(tmp)[j], node=tmp[-j,j]) ) )
					
					mph.mapnode.pctime	<- mph.clu.itopo[,	{					
								
										tmp			<- as.matrix( mrca( mph.clu[[ mph.i[1] ]], full = FALSE ) )					
										mph.mrca	<- do.call('rbind', lapply( seq_len(ncol(tmp)), function(j)		data.table( tip1=rownames(tmp)[-j], tip2=colnames(tmp)[j], mph.node=tmp[-j,j]) ) )
										mph.mrca	<- merge(mph.mrca, ref.mrca, by=c('tip1','tip2'))
										mph.mrca	<- rbind(mph.mrca, data.table(tip1='root',tip2='root', mph.node=0, node=0))
										mph.info	<- hivc.treeannotator.tiplabel2df(mph.clu[[ mph.i[1] ]], beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
										mph.info[, tip:=match( mph.info[, BEASTlabel], mph.clu[[ mph.i[1] ]]$tip.label)]															
										node		<- c(seq.int(from=Ntip(mph.clu[[1]])+1, to=Nnode(mph.clu[[1]], internal=FALSE)),0)
										#	extract ctimes for nodes in mph topology
										tmp			<- sapply( mph.i, function(i)
												{
													depth	<- c( node.depth.edgelength( mph.clu[[ i ]] ), -mph.clu[[ i ]]$root.edge )
													tmp		<- which.max(depth)
													depth	<- depth-depth[tmp]+subset(mph.info, tip==tmp)[, TipT]
													depth[ seq.int(from=Ntip(mph.clu[[ i ]])+1, to=length(depth)) ]																		
												})
										tmp	<- data.table(mph.node=rep(node, each=ncol(tmp)), ctime= as.numeric(t(tmp)))
										#	every node in map topology can correspond to multipe mrca s in the mph topology
										#	to weigh the ctimes correctly, count nodes multiple times according to how often they are the mrca s in the map topology
										#	this may be mem intensive
										tmp	<- merge(tmp, subset(mph.mrca, select=c(node, mph.node)), by='mph.node', allow.cartesian=TRUE)
										list( node=tmp[,node], ctime=tmp[, ctime] )															
									}, by='equal.to']
					mph.mapnode.pctime	<- mph.mapnode.pctime[, list( q= round(quantile(ctime, probs=seq(0,1,by=cdf.by)),d=3), cdf=seq(0,1,by=cdf.by) ), by='node']		
					mph.mapnode.pctime	<- mph.mapnode.pctime[,	{
																	tmp<- setdiff(seq_along(q),which(duplicated(q))-1)
																	if(tmp[1]!=1)
																		tmp[1]	<- 1
																	if(length(tmp)==1)	
																		ans		<- list(q=q[tmp]+c(-0.001,0), cdf=c(0.,1.))
																	else
																		ans		<- list(q=q[tmp], cdf=cdf[tmp])
																	ans
																}, by='node']
					mph.mapnode.pctime	<- mph.mapnode.pctime[, {
																	tmp	<- c(0,diff(cdf)/diff(q))
																	list(q=q, pdf=tmp/sum(tmp), cdf=cdf)	
																}, by='node']							
					mph.mapnode.pctime[, cluster:= clu]
					mph.mapnode.pctime[, mph.i:= mph.clu.dtopo[, mph.i[which.max(freq)]]]					
					#
					#	tips and edges in same order, so very easy to summarize brls of each topology
					#
					cat(paste('\ncalculate brl.cdf'))
					mph.brl.cdf		<- mph.clu.itopo[,	{					
															edge.of	<- c(mph.clu[[ mph.i[1] ]]$edge[,2], 0)
															tmp		<- sapply( mph.i, function(i)	c( mph.clu[[ i ]]$edge.length, mph.clu[[ i ]]$root.edge)  )
															tmp		<- apply(tmp, 1, quantile, probs=seq(0,1,by=cdf.by))																
															list( q= round(as.numeric(tmp),d=4), cdf=seq(0,1,by=cdf.by), edge.of=rep(edge.of, each=nrow(tmp)) )					
														}, by='equal.to']	
					mph.brl.cdf		<- mph.brl.cdf[,	{
															tmp<- setdiff(seq_along(q),which(duplicated(q))-1)
															if(tmp[1]!=1)
																tmp[1]	<- 1
															if(length(tmp)==1)	
																ans		<- list(cluster= rep(cluster[tmp],2), q=q[tmp]+c(-0.001,0), cdf=c(0.,1.))
															else
																ans		<- list(cluster= cluster[tmp], q=q[tmp], cdf=cdf[tmp])
															ans
														},by=c('equal.to','edge.of')]
					mph.brl.cdf		<- mph.brl.cdf[, {
														tmp	<- c(0,diff(cdf)/diff(q))
														list(cluster=cluster, q=q, pdf=tmp/sum(tmp), cdf=cdf)	
													},by=c('equal.to','edge.of')]							
					mph.brl.cdf[, cluster:= clu]
					#
					#	tips and edges in same order, so very easy to compute the prob that a tip is an internal node
					#
					cat(paste('\ncalculate SA.cnt'))
					mph.SA.cnt		<- mph.clu.itopo[,	{
															mph.topo.ntip	<- Ntip(mph.clu[[ mph.i[1] ]])
															mph.topo.tipn	<- mph.clu[[ mph.i[1] ]]$tip.label
															mph.topo.itip	<- which( mph.clu[[ mph.i[1] ]]$edge[,2] <= mph.topo.ntip )
															tmp				<- sapply( mph.i, function(i)	mph.clu[[ i ]]$edge.length[mph.topo.itip]  )
															list(tip=mph.topo.tipn, SA.freq= apply(tmp, 1, function(x) length(which(x==0.))), n= length(mph.i) )																						
														}, by='equal.to']
					mph.SA.cnt[, cluster:= clu]
					#
					#	generate consensus tree: take the posterior tree such that all its branch lengths are closest to the marginal median branch lengths
					#	if we took the median branch lengths, then the branch lengths might not any longer be consistent with the sampling dates
					#
					mph.brl			<- mph.clu.itopo[,	{					
															edge.of			<- c(mph.clu[[ mph.i[1] ]]$edge[,2], 0)
															tmp				<- sapply( mph.i, function(i)	c( mph.clu[[ i ]]$edge.length, mph.clu[[ i ]]$root.edge)  )
															tmp				<- lapply( seq_len(nrow(tmp)), function(i)		tmp[i,])	
															names(tmp)		<- paste('edge.of.',edge.of,sep='')
															tmp																				
														}, by='equal.to']	
					tmp				<- mph.brl[,  lapply(.SD, rank, ties.method='random'), by='equal.to']
					tmp				<- tmp[, 	{
													z		<- do.call('cbind',lapply( .SD, function(x) abs(x-length(edge.of.0)/2)))
													list(dist.from.midrank= apply(z, 1, sum) )								
												},by='equal.to']					
					mph.brl.mid		<- tmp[, list(mph.ii=which.min(dist.from.midrank)), by='equal.to']
												
					ph.consensus	<- lapply(seq_len(nrow(mph.clu.dtopo)), function(j)
							{								
								mph.ii					<- subset(mph.brl.mid, equal.to==mph.clu.dtopo[j, mph.i])[, mph.ii]
								mph.i					<- subset(mph.clu.itopo, equal.to==mph.clu.dtopo[j, mph.i])[, mph.i]
								mph.clu[[ mph.i[mph.ii] ]]								
							})					 
					names(ph.consensus)	<- paste('cluster=',clu,'_mph.i=', mph.clu.dtopo[,mph.i],'_prob=',mph.clu.dtopo[,dens], sep='')
					if(0)
					{
						mph.info		<- hivc.treeannotator.tiplabel2df(ph.consensus[[1]], beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
						mph.info[, tip:=match( mph.info[, BEASTlabel], ph.consensus[[1]]$tip.label)]
						
						tmp	<- node.depth.edgelength( ph.consensus[[1]] )
						tmp	<- mph.info[1,TipT]-tmp[1] + tmp
						tmp	<- tmp[seq_len(nrow(mph.info))] - mph.info[,TipT]
						if( any( abs(tmp)>EPS ) ) warning('branch lengths do not add up to tip times')
					}
					#
					#	save
					#
					file	<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),'_cluposterior_',clu,'.R',sep='')					
					cat(paste("\nsave ph.consensus, mph.SA.cnt, mph.brl.cdf, mph.node.ctime, mph.clu.dtopo, mph.clu.itopo, mph.mapnode.pctime to file",file))
					save(ph.consensus, mph.SA.cnt, mph.brl.cdf, mph.clu.dtopo, mph.clu.itopo, mph.node.ctime, mph.mapnode.pctime, file=file)
				}								
			})
}
######################################################################################
hivc.prog.BEAST2.get.cluster.trees<- function()
{
	require(ape)
	require(data.table)
	require(hivclust)
	#	program options
	opt.pool					<- NA
	opt.burnin					<- 5e6		
	opt.save.cluster.nexfile	<- FALSE
	resume						<- TRUE
	verbose						<- TRUE
	beastlabel.idx.clu			<- 1
	beastlabel.idx.hivn			<- 2
	beastlabel.idx.hivd			<- 3
	beastlabel.idx.hivs			<- 4
	beastlabel.idx.samplecode	<- 6
	beastlabel.idx.rate			<- NA
	#	input files		
	indir					<- paste(DATA,"tmp",sep='/')
	indir					<- '/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/zip'
	infile					<- "ATHENA_2013_03_-DR-RC-SH+LANL_Sequences"
	insignat				<- "Wed_Dec_18_11:37:00_2013"
	infilexml.opt			<- "clrh80"
	infilexml.template		<- "sasky_sdr06fr"
	#
	#indir				<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/beast/beast_131011"		
	#infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
	#insignat			<- "Tue_Aug_26_09:13:47_2013"				
	#infilexml.template	<- "um182rhU2045"
	#infilexml.opt		<- "mph4clutx4tip"
	#opt.burnin			<- 2e7
	
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
						{	switch(substr(arg,2,14),
									infilexml.opt= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.opt<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									infilexml.template= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.template<- tmp[1]
		#
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) outsignat<- tmp[1]				
		#		
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									burnin= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) opt.burnin<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,5),
									pool= return(as.numeric(substr(arg,7,nchar(arg)))),NA)	}))
		if(length(tmp)>0) opt.pool<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]		
	}
	#
	outdir					<- indir
	outfile					<- paste(infile,infilexml.template,infilexml.opt, sep='_')
	outsignat				<- insignat
	#
	if(verbose)
	{
		print(opt.pool)
		print(opt.burnin)		
		print(opt.save.cluster.nexfile)
		print(resume)
		print(beastlabel.idx.clu)
		print(beastlabel.idx.hivn)
		print(beastlabel.idx.hivd)
		print(beastlabel.idx.hivs)
		print(beastlabel.idx.samplecode)
		print(beastlabel.idx.rate)
		print(indir)		
		print(infile)
		print(insignat)
		print(infilexml.opt)
		print(infilexml.template)
		print(outdir)
		print(outfile)
		print(outsignat)
	}
	# read log files and return mean TMRCA and file specific TMRCA
	#files		<- list.files(indir)
	#files		<- files[ sapply(files, function(x) grepl(infile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(paste('_',infilexml.opt,'_',sep=''), x, fixed=1) & grepl(paste('_',infilexml.template,'_',sep=''), x, fixed=1) & grepl('log$',x) ) ]
	#files		<- paste(indir, files, sep='/')
	#if(!length(files))	stop('cannot find files matching criteria')
	#df.th		<- hivc.beast2out.get.log.evolrate(files, opt.burnin=opt.burnin)
	
	#	read tree files	- this takes a while
	#	do this one by one as reading all in one go is too mem intensive
	files		<- list.files(indir)
	files		<- files[ sapply(files, function(x) grepl(infile, x, fixed=1) & grepl(gsub('/',':',insignat), x, fixed=1) & grepl(paste('_',infilexml.opt,'_',sep=''), x, fixed=1) & grepl(paste('_',infilexml.template,'_',sep=''), x, fixed=1) & grepl('trees$',x) ) ]	
	if(!length(files))	stop('cannot find files matching criteria')
	files		<- paste(indir, files, sep='/')
	tmp			<- grepl('timetrees$',files)
	if(any(tmp))
		files	<- files[tmp]
	if(!is.na(opt.pool))
		files	<- files[ grepl(paste('_pool_',opt.pool,'_',sep=''), files) ]
	cat(paste("\nfound files, n=",length(files)))
	dummy		<- lapply(seq_along(files), function(i)
			{
				file					<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),'_pool_',i,'.R',sep='')
				if(resume)
				{
					options(show.error.messages = FALSE)		
					readAttempt		<- try(suppressWarnings(load(file)))
					if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))					
					options(show.error.messages = TRUE)							
				}
				if(!resume || inherits(readAttempt, "try-error"))
				{
					#opt.rescale.edge.length	<- df.th[i, ucldMean.mean / ucldMean.mean.i]
					cat(paste("\nread file ",files[i]))
					mph						<- hivc.beast2out.read.trees(files[i], opt.rescale.edge.length=1, opt.burnin=opt.burnin)					
					cat(paste('\nsave mph to file=',file))
					tmp						<- hivc.beast2out.tip.date.check(mph[[length(mph)]], fun=hivc.treeannotator.tiplabel2df, beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
					print(tmp)
					if(tmp>2*EPS )	stop('hivc.beast2out.read.trees does not pass random tip date check')						
					save(mph, file=file)					
				}
				mph						<- NULL
				gc()				
			})
	
	#	for each cluster, extract monophyletic subtrees corresponding to clusters
	dummy		<- lapply(seq_along(files), function(i)
			{
				if(resume)
				{
					tmp			<- list.files(outdir)
					tmp			<- tmp[ sapply(tmp, function(x) grepl(outfile, x, fixed=1) & grepl(gsub('/',':',outsignat), x, fixed=1) & grepl(infilexml.opt, x, fixed=1) & grepl(infilexml.template,x, fixed=1) & grepl(paste('_pool_',i,'_clu_[0-9]+',sep=''),x) & grepl('.R$',x) ) ]
					tmp			<- length(tmp)
				}
				if(!resume || !tmp )
				{
					file			<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),'_pool_',i,'.R',sep='')
					cat(paste("\nload file",file))
					options(show.error.messages = FALSE)		
					readAttempt		<- try(suppressWarnings(load(file)))
					if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",file))
					if(inherits(readAttempt, "try-error"))	stop(paste("\ncould not read file",file))
					options(show.error.messages = TRUE)		
					mph.info		<- hivc.treeannotator.tiplabel2df(mph[[1]], beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
					mph.by.clu		<- hivc.beast2out.extract.cluster.trees(mph, mph.info)
					mph				<- NULL					
					#for each cluster, produce R and possibly nexus tree file
					dummy			<- lapply(seq_along(mph.by.clu), function(clu)
							{
								file	<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),'_pool_',i,'_clu_',names(mph.by.clu)[clu],'.R',sep='')
								mph.clu	<- mph.by.clu[[clu]]
								cat(paste("\nsave mph.clu to file",file))
								save(mph.clu, file=file)
								if(opt.save.cluster.nexfile)
								{
									file	<- paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),'_pool_',i,'_clu_',names(mph.by.clu)[clu],'.nex',sep='')
									cat(paste("\nsave mph.clu to file",file))
									write.nexus(mph.clu, file=file)					
								}
							})
					dummy			<- NULL
					gc()	
				}					
			})		
}
######################################################################################
hivc.prog.BEAST2.generate.xml<- function()
{	
	require(XML)
	
	indir				<- paste(DATA,"tmp",sep='/')		
	infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
	insignat			<- "Thu_Aug_01_17/05/23_2013"
	indircov			<- paste(DATA,"derived",sep='/')
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infiletree			<- paste(infile,"examlbs100",sep="_")
	infilexml			<- paste(infile,'_',"seroneg",sep='')	
	outdir				<- indir
	outsignat			<- "Tue_Aug_26_09/13/47_2013"
	
	opt.burnin				<- 5e6
	opt.brl					<- "dist.brl.casc" 
	thresh.brl				<- 0.096
	thresh.bs				<- 0.8
	pool.ntip				<- NA
	beast.mcmc.length		<- 25e7
	infilexml.opt			<- "rbe420"
	infilexml.template		<- "sasky_sdr06"	
	
	resume				<- 1
	verbose				<- 1
	hpc.walltime		<- 800
	hpc.ncpu			<- 1
	hpc.mem				<- "1200mb"
	
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
		if(length(tmp)>0) infilecov<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infiletree= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infiletree<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infilexml= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) outsignat<- tmp[1]				
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
						{	switch(substr(arg,2,8),
									opt.brl= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.brl<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									thresh.bs= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.bs<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									thresh.brl= return(as.numeric(substr(arg,13,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.brl<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									pool.ntip= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) pool.ntip<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,14),
									infilexml.opt= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.opt<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									infilexml.template= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.template<- tmp[1]
	}	
	if(grepl("sdr06",infilexml.template))
		alignment.filter<- c("1::3,2::3", "3::3")
	else
		alignment.filter<- NA
	#	modify beast2.spec depending on infilexml.opt
	if(grepl("S4p",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=4, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.)
	}
	else if(grepl("S5p",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 7.596, 5.596, 1.596, 0.)
	}
	else if(grepl("S8p",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=8, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 8.596, 7.596, 6.596, 5.596, 1.596, 0.596, 0.)
	}
	else if(grepl("s0106",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.2, 0.7, 0.6, 0.6)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.1/0","Uniform/0.1/1.0","Uniform/0.6/1.0","Uniform/0.5/1.0","Uniform/0.5/1.0")	
	}
	else if(grepl("s00106",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.01, 0.2, 0.7, 0.6, 0.6)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Uniform/0.1/1.0","Uniform/0.6/1.0","Uniform/0.5/1.0","Uniform/0.5/1.0")	
	}	
	else if(grepl("s0108",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.8, 0.8)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Uniform/0.4/1.0","Uniform/0.8/1.0","Uniform/0.7/1.0","Uniform/0.7/1.0")	
	}
	else if(grepl("s124",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.8, 0.8)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.8/1.0","Uniform/0.2/1.0","Beta/2.5/4.0/0")	
	}	
	else if(grepl("s424",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.8, 0.8)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.4/0","Uniform/0.8/1.0","Uniform/0.2/1.0","Beta/2.5/4.0/0")	
	}
	else if(grepl("s024",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.8, 0.8)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.01/0","Uniform/0.8/1.0","Uniform/0.2/1.0","Beta/2.5/4.0/0")	
	}
	else if(grepl("s184",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")	
	}
	else if(grepl("sartest",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$sasky.r.value					<- c(0.1, 0.2, 0.5, 0.7, 0.7)
		beast2.spec$sasky.r.prior					<- c("Uniform/0.0/0.5","Uniform/0.0/0.5","Uniform/0.0/1.0","Uniform/0.5/0.8","Uniform/0.5/1.0")
	}
	else if(grepl("d999",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(9, 9, 9, 9, 9)
		beast2.spec$bdsky.notInf.prior				<- c("Exponential/0.11/0","Exponential/0.11/0","Exponential/0.11/0","Exponential/0.11/0","Exponential/0.11/0")
	}
	else if(grepl("d774",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(7, 7, 7, 4, 4)
		beast2.spec$bdsky.notInf.prior				<- c("Exponential/0.14/0","Exponential/0.14/0","Exponential/0.14/0","Exponential/0.25/0","Exponential/0.25/0")
	}
	else if(grepl("d543",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)
		beast2.spec$bdsky.notInf.prior				<- c("Exponential/0.2/0","Exponential/0.25/0","Exponential/0.25/0","Exponential/0.33/0","Exponential/0.33/0")
	}
	else if(grepl("dg543",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)
		beast2.spec$bdsky.notInf.prior				<- c("Gamma/5/0.03/0.1","Gamma/5/0.03/0.1","Exponential/0.25/0","Gamma/5/0.05/0.1","Exponential/0.33/0")
		beast2.spec$sasky.r.value					<- c(0.9, 0.5, 0.5, 0.5, 0.5)
		beast2.spec$sasky.r.prior					<- c("Beta/2/1/0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0")
	}
	else if(grepl("r543",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)
		#beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.2/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Beta/4.0/3.0/0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)
		#beast2.spec$bdsky.notInf.prior				<- c("Gamma/5/0.03/0.1","Gamma/5/0.03/0.1","Exponential/0.25/0","Gamma/5/0.05/0.1","Exponential/0.33/0")
		beast2.spec$bdsky.notInf.prior				<- c("LogNormal/0.2/0.6/0.1/true","LogNormal/0.25/0.8/0.1/true","LogNormal/0.3/0.8/0.1/true","LogNormal/0.5/1.2/0.1/true","LogNormal/0.5/1.2/0.1/true")
		beast2.spec$sasky.r.prior					<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt))) / 100
		beast2.spec$sasky.r.value					<- rep( beast2.spec$sasky.r.prior / 2, 5 )
		beast2.spec$sasky.r.prior					<- rep( paste("Uniform/0.0/",beast2.spec$sasky.r.prior,sep=''), 5 )
	}
	else if(grepl("ori",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.1, 0.6, 0.6, 0.3)
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Beta/4.0/3.0/0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)
		beast2.spec$bdsky.notInf.prior				<- c("LogNormal/0.2/0.6/0.1/true","LogNormal/0.25/0.8/0.1/true","LogNormal/0.3/0.8/0.1/true","LogNormal/0.5/1.2/0.1/true","LogNormal/0.5/1.2/0.1/true")
		beast2.spec$sasky.r.value					<- rep(0.05, 5)
		beast2.spec$sasky.r.prior					<- rep("Uniform/0.0/0.1",5)		
		beast2.spec$bdsky.origin.prior				<- as.numeric(substr(infilexml.opt,4,nchar(infilexml.opt)))		
		beast2.spec$bdsky.origin.prior				<- paste("Uniform/",min(20,beast2.spec$bdsky.origin.prior-20),"/",beast2.spec$bdsky.origin.prior,sep='')		
	}
	else if(grepl("rbe9",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)		
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Beta/9.0/5.0/0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)		
		beast2.spec$bdsky.notInf.prior				<- c("LogNormal/0.2/0.6/0.1/true","LogNormal/0.25/0.8/0.1/true","LogNormal/0.3/0.8/0.1/true","LogNormal/0.5/1.2/0.1/true","LogNormal/0.5/1.2/0.1/true")
		beast2.spec$sasky.r.prior					<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt))) / 100
		beast2.spec$sasky.r.value					<- rep( beast2.spec$sasky.r.prior / 2, 5 )
		beast2.spec$sasky.r.prior					<- rep( paste("Uniform/0.0/",beast2.spec$sasky.r.prior,sep=''), 5 )
	}
	else if(grepl("run2",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.9, 0.6, 0.3)		
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.25/1.0","Beta/4.0/3.0/0","Beta/2.5/4.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)
		#beast2.spec$bdsky.notInf.prior				<- c("Gamma/5/0.03/0.1","Gamma/5/0.03/0.1","Exponential/0.25/0","Gamma/5/0.05/0.1","Exponential/0.33/0")
		beast2.spec$bdsky.notInf.prior				<- c("LogNormal/0.2/0.6/0.1/true","LogNormal/0.25/0.8/0.1/true","LogNormal/0.3/0.8/0.1/true","LogNormal/0.5/1.2/0.1/true","LogNormal/0.5/1.2/0.1/true")
		beast2.spec$sasky.r.prior					<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt))) / 100
		beast2.spec$sasky.r.value					<- rep( beast2.spec$sasky.r.prior / 2, 5 )
		beast2.spec$sasky.r.prior					<- rep( paste("Uniform/0.0/",beast2.spec$sasky.r.prior,sep=''), 5 )
	}
	else if(grepl("rbe4",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.2, 0.2, 0.2)		
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Beta/9.0/5.0/0","Beta/4.0/16.0/0","Beta/4.0/16.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)
		beast2.spec$bdsky.notInf.prior				<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt))) / 100
		beast2.spec$bdsky.notInf.prior				<- c("LogNormal/0.2/0.6/0.1/true", rep( paste("LogNormal/",beast2.spec$bdsky.notInf.prior,"/0.8/0.1/true",sep=''), 4) )				
		beast2.spec$sasky.r.value					<- rep( 0.5, 5 )
		beast2.spec$sasky.r.prior					<- rep( "Uniform/0.0/1.0", 5 )
	}
	else if(grepl("pmct",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.2, 0.2, 0.2)		
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Beta/9.0/5.0/0","Beta/4.0/16.0/0","Beta/4.0/16.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)		
		beast2.spec$bdsky.notInf.prior				<- rep( "LogNormal/0.2/0.6/0.1/true",5)				
		beast2.spec$sasky.r.value					<- rep( 0.5, 5 )
		beast2.spec$sasky.r.prior					<- rep( "Uniform/0.0/1.0", 5 )
		beast2.spec$pool.cnts.requested				<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt))) 
		beast2.spec$pool.cnts.requested				<- c(NA, beast2.spec$pool.cnts.requested, 70, 50, NA)				
	}
	else if(grepl("pfn6",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$pool.fNegT						<- 0.6
		beast2.spec$pool.ntip						<- 200
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.2, 0.2, 0.2)		
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Beta/9.0/5.0/0","Beta/4.0/16.0/0","Beta/4.0/16.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)		
		beast2.spec$bdsky.notInf.prior				<- rep( "LogNormal/0.2/0.6/0.1/true",5)				
		beast2.spec$sasky.r.value					<- rep( 0.5, 5 )
		beast2.spec$sasky.r.prior					<- rep( "Uniform/0.0/1.0", 5 )
		beast2.spec$pool.cnts.requested				<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt))) 
		beast2.spec$pool.cnts.requested				<- c(NA, beast2.spec$pool.cnts.requested, max(beast2.spec$pool.cnts.requested,70), max(50,beast2.spec$pool.cnts.requested), NA)				
	}
	else if(grepl("piv4",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=4, alignment.filter=alignment.filter)
		beast2.spec$pool.fNegT						<- 0.6
		beast2.spec$pool.ntip						<- 200
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.2, 0.2)		
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Beta/9.0/5.0/0","Beta/4.0/16.0/0")
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3)		
		beast2.spec$bdsky.notInf.prior				<- rep( "LogNormal/0.2/0.6/0.1/true",4)				
		beast2.spec$sasky.r.value					<- rep( 0.5, 4 )
		beast2.spec$sasky.r.prior					<- rep( "Uniform/0.0/1.0", 4 )
		beast2.spec$pool.cnts.requested				<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt))) 
		beast2.spec$pool.cnts.requested				<- c(max(50,beast2.spec$pool.cnts.requested), max(beast2.spec$pool.cnts.requested,70), beast2.spec$pool.cnts.requested,  NA )				
	}
	else if(grepl("rfn8",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.2, 0.2, 0.2)		
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Beta/9.0/5.0/0","Beta/4.0/16.0/0","Beta/4.0/16.0/0")
		beast2.spec$bdsky.R0.prior					<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt))) / 100
		beast2.spec$bdsky.R0.prior					<- c("Gamma/1.5/1.5/0", rep(paste("Gamma/2.3/",beast2.spec$bdsky.R0.prior,"/0.7",sep=''),4))
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)		
		beast2.spec$bdsky.notInf.prior				<- rep( "LogNormal/0.2/0.6/0.1/true",5)				
		beast2.spec$sasky.r.value					<- rep( 0.5, 5 )
		beast2.spec$sasky.r.prior					<- rep( "Uniform/0.0/1.0", 5 )		 
		beast2.spec$pool.cnts.requested				<- c(NA, 35, 70, 50, NA)				
	}
	else if(grepl("rsu8",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter, cluster.monophyletic=TRUE, tip.log.stem=TRUE)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.2, 0.2, 0.2)		
		beast2.spec$bdsky.sprop.prior				<- c("Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0")
		beast2.spec$bdsky.R0.prior					<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt))) / 100
		beast2.spec$bdsky.R0.prior					<- c("Gamma/1.5/1.5/0", rep(paste("Gamma/2.3/",beast2.spec$bdsky.R0.prior,"/0.7",sep=''),4))
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)		
		beast2.spec$bdsky.notInf.prior				<- rep( "LogNormal/0.2/0.6/0.1/true",5)				
		beast2.spec$sasky.r.value					<- rep( 0.5, 5 )
		beast2.spec$sasky.r.prior					<- rep( "Uniform/0.0/1.0", 5 )		 
		beast2.spec$pool.cnts.requested				<- c(NA, 35, 70, 50, NA)				
	}
	else if(grepl("rse8",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter, cluster.monophyletic=TRUE, tip.log.stem=TRUE)
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.2, 0.2, 0.2)		
		beast2.spec$bdsky.sprop.prior				<- c("Exponential/0.01/0","Exponential/0.1/0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0")
		beast2.spec$bdsky.R0.prior					<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt))) / 100
		beast2.spec$bdsky.R0.prior					<- c("Gamma/1.5/1.5/0", rep(paste("Gamma/2.3/",beast2.spec$bdsky.R0.prior,"/0.7",sep=''),4))
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)		
		beast2.spec$bdsky.notInf.prior				<- rep( "LogNormal/0.2/0.6/0.1/true",5)				
		beast2.spec$sasky.r.value					<- rep( 0.5, 5 )
		beast2.spec$sasky.r.prior					<- rep( "Uniform/0.0/1.0", 5 )		 
		beast2.spec$pool.cnts.requested				<- c(NA, 35, 70, 50, NA)				
	}
	else if(grepl("alsu",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter, cluster.monophyletic=TRUE, cluster.log=FALSE, tip.log.stem=FALSE)
		beast2.spec$pool.fNegT						<- 0
		beast2.spec$pool.ntip						<- 150		
		beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.596, 0.)
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.2, 0.2, 0.2)		
		beast2.spec$bdsky.sprop.prior				<- c("Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0")
		beast2.spec$bdsky.R0.prior					<- c("Gamma/1.5/1.5/0", rep(paste("Gamma/2.3/0.15/0.7",sep=''),4))
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)		
		beast2.spec$bdsky.notInf.prior				<- rep( "LogNormal/0.2/0.6/0.1/true",5)				
		beast2.spec$sasky.r.value					<- rep( 0.5, 5 )
		beast2.spec$sasky.r.prior					<- rep( "Uniform/0.0/1.0", 5 )
		beast2.spec$pool.cnts.requested				<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt)))
		beast2.spec$pool.cnts.requested				<- c(NA, max(50,beast2.spec$pool.cnts.requested), 70, NA, NA)				
	}	
	else if(grepl("alrh",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter, cluster.monophyletic=TRUE, cluster.log=FALSE, tip.log.stem=FALSE)
		beast2.spec$pool.fNegT						<- 0
		beast2.spec$pool.ntip						<- 150		
		beast2.spec$sasky.r.changepoint.value		<- beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(0.596, 1.596, 5.596, 9.596, 0.)
		beast2.spec$bdsky.reverseTimeArrays.value	<- c('true', 'true', 'true', 'false', 'true')
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.2, 0.2, 0.2)		
		beast2.spec$bdsky.sprop.prior				<- c("Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0")
		beast2.spec$bdsky.R0.prior					<- c("Gamma/1.5/1.5/0", rep(paste("Gamma/2.3/0.15/0.7",sep=''),4))
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)		
		beast2.spec$bdsky.notInf.prior				<- rep( "LogNormal/0.2/0.6/0.1/true",5)				
		beast2.spec$sasky.r.value					<- rep( 0.5, 5 )
		beast2.spec$sasky.r.prior					<- rep( "Uniform/0.0/1.0", 5 )		
		beast2.spec$pool.cnts.requested				<- c(NA, 50, 70, NA, NA)	
		beast2.spec$bdsky.origin.prior				<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt)))		
		beast2.spec$bdsky.origin.prior				<- paste("Uniform/",min(20,beast2.spec$bdsky.origin.prior-20),"/",beast2.spec$bdsky.origin.prior,sep='')
	}	
	else if(grepl("clrh",infilexml.opt))
	{
		beast2.spec		<- hivc.beast2.get.specifications(mcmc.length=beast.mcmc.length, bdsky.intervalNumber=5, alignment.filter=alignment.filter, cluster.monophyletic=TRUE, cluster.log=FALSE, tip.log.stem=FALSE)
		beast2.spec$pool.fNegT						<- 0
		beast2.spec$pool.ntip						<- 150		
		beast2.spec$sasky.r.changepoint.value		<- beast2.spec$bdsky.sprop.changepoint.value	<- beast2.spec$bdsky.R0.changepoint.value		<- beast2.spec$bdsky.notInf.changepoint.value	<- c(1.596, 3.596, 5.596, 9.596, 0.)
		beast2.spec$bdsky.reverseTimeArrays.value	<- c('true', 'true', 'true', 'false', 'true')
		beast2.spec$bdsky.sprop.value				<- c(0.1, 0.5, 0.2, 0.2, 0.2)		
		beast2.spec$bdsky.sprop.prior				<- c("Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0","Uniform/0.0/1.0")
		beast2.spec$bdsky.R0.prior					<- c("Gamma/1.5/1.5/0", rep(paste("Gamma/2.3/0.15/0.7",sep=''),4))
		beast2.spec$bdsky.notInf.value				<- 1/c(5, 4, 4, 3, 3)		
		beast2.spec$bdsky.notInf.prior				<- rep( "LogNormal/0.2/0.6/0.1/true",5)				
		beast2.spec$sasky.r.value					<- rep( 0.5, 5 )
		beast2.spec$sasky.r.prior					<- rep( "Uniform/0.0/1.0", 5 )		
		beast2.spec$pool.cnts.requested				<- c(50, 70, 70, NA, NA)	
		beast2.spec$bdsky.origin.prior				<- as.numeric(substr(infilexml.opt,5,nchar(infilexml.opt)))		
		beast2.spec$bdsky.origin.prior				<- paste("Uniform/",min(20,beast2.spec$bdsky.origin.prior-20),"/",beast2.spec$bdsky.origin.prior,sep='')
	}	
	else stop("unknown infilexml.opt")
	#		 
	#
	#
	#
	if(!is.na(pool.ntip))
		beast2.spec$pool.ntip			<- pool.ntip
	if(grepl("sasky",infilexml.template))
	{
		beast2.spec$template			<- infilexml.template
		beast2.spec$treemodel			<- "SampledAncestorSkylineModel"
		prog.beast						<- PR.BEAST2SA
	}
	if(grepl("bdsky",infilexml.template))
	{
		beast2.spec$template			<- infilexml.template
		beast2.spec$treemodel			<- "BirthDeathSkylineModel"
		prog.beast						<- PR.BEAST2
	}
	#	
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infilecov)
		print(infiletree)
		print(infilexml)
		print(outdir)		
		print(outsignat)
		print(resume)
		print(opt.brl)
		print(thresh.brl)
		print(thresh.bs)
		print(pool.ntip)		 
		print(infilexml.opt)
		print(infilexml.template)
	}	
	#
	#	load complete tree to generate starting tree
	#
	file					<- paste(indir,'/',infiletree,'_',gsub('/',':',insignat),".R",sep='')
	if(verbose)	cat(paste("\nload complete tree to generate starting tree from file",file))
	load(file)	#load object 'ph'	
	#
	#	load sequences
	#
	file		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
	if(verbose)	cat(paste("\nload sequences from file",file))
	load( file )
	#print( seq.PROT.RT )	
	#
	#	load msm clusters
	#
	argv		<<- hivc.cmd.clustering.msm(indir, infiletree, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=resume)
	argv		<<- unlist(strsplit(argv,' '))		
	msm			<- hivc.prog.get.clustering.MSM()	
	#
	#	select seroconverters
	#
	df.cluinfo				<- msm$df.cluinfo
	tmp						<- df.cluinfo[,	list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], fNegT=length(which(!is.na(NegT))) / clu.ntip[1]),by="cluster"]										
	tmp						<- subset(tmp, fNegT>=quantile(tmp[,fNegT], probs=beast2.spec$pool.fNegT) )
	cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
	if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))
	cluphy.df				<- hivc.beast.addBEASTLabel( cluphy.df )
	#
	#	create sets of cluster pools for BEAST
	#
	if(any(!is.na(beast2.spec$pool.cnts.requested)))
		df.clupool			<- hivc.beast2.poolclusters.mincnts(cluphy.df, beast2.spec, verbose=1)
	else
		df.clupool			<- hivc.beast.poolclusters(cluphy.df, pool.ntip= beast2.spec$pool.ntip, verbose=1)
	#print( sapply(seq_along(df.clupool$pool.df), function(i)	nrow(df.clupool$pool.df[[i]])) )
	#
	#	load xml template file
	#	
	file			<- paste(CODE.HOME,"/data/BEAST2_template_",beast2.spec$template,".xml",sep='')
	bxml.template	<- xmlTreeParse(file, useInternalNodes=TRUE, addFinalizer = TRUE)			
	#
	#	create BEAST2 XML file	
	#
	bfile			<- lapply(seq_len(length(df.clupool$pool.df)), function(pool.id)
						{
							df							<- df.clupool$pool.df[[pool.id]]
							setkey(df, cluster)							
							beast2.spec$xml.dir			<- indir
							beast2.spec$xml.filename	<- paste(infilexml,'-',beast2.spec$pool.ntip,'_pool_',pool.id,'_',infilexml.template,'_',infilexml.opt,'_',gsub('/',':',outsignat),sep='')
							beast2.spec$starttree.newick<- hivc.beast2.get.startingtree(ph, df, beast2.spec, verbose=verbose)
							beast2.xml					<- hivc.beast2.get.xml( bxml.template, seq.PROT.RT, df, beast2.spec, verbose=1)			
							file						<- paste(beast2.spec$xml.dir,'/',beast2.spec$xml.filename,".xml", sep='')
							if(verbose)	cat(paste("\nwrite xml file to",file))
							saveXML(beast2.xml, file=file)
							paste(infilexml,'-',beast2.spec$pool.ntip,'_pool_',pool.id,'_',infilexml.template,'_',infilexml.opt,sep='')
						})
	#
	#	generate BEAST commands and run
	#
	sapply(seq_along(bfile), function(pool.id)
			{
				cmd			<- hivc.cmd.beast2.runxml(indir, bfile[[pool.id]], outsignat, hpc.ncpu=hpc.ncpu, prog.beast=prog.beast, prog.opt.Xmx="1200m", hpc.tmpdir.prefix="beast2")
				#cmd		<- paste(cmd, hivc.cmd.beast2.getclustertrees.pipe(indir, infilexml, outsignat, infilexml.opt, infilexml.template, opt.burnin, outdir=indir, outsignat=outsignat, opt.pool=pool.id, verbose=verbose, resume=resume), sep='' )
				#cmd		<- paste(cmd,hivc.cmd.beast.evalrun(outdir, infilexml, outsignat, infilexml.opt, infilexml.template, length(bfile), verbose=1),sep='')				
				cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=hpc.walltime, hpc.q="pqeph", hpc.mem=hpc.mem,  hpc.nproc=hpc.ncpu)					
				cat(cmd)
				outfile		<- paste("b2",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
				hivc.cmd.hpccaller(outdir, outfile, cmd)
				stop()
			})	
}
######################################################################################
hivc.prog.BEAST.read.nexus.and.stats<- function()
{	
	require(ape)
	indir				<- paste(DATA,"tmp",sep='/')		
	infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"			
	outdir				<- indir	
	
	resume				<- 1
	verbose				<- 1
	tree.id				<- NA 
	method.node.stat	<- 'any.node'
	
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
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									tree.id= return(as.numeric(substr(arg,10,nchar(arg)))),NA)	}))
		if(length(tmp)>0) tree.id<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,17),
									method.node.stat= return(substr(arg,19,nchar(arg))),NA)	}))
		if(length(tmp)>0) method.node.stat<- tmp[1]
	}	
	outfile		<- unlist(strsplit(infile, '\\.'))
	outfile		<- paste( c(outfile[-length(outfile)],'R'),sep='',collapse='.')
	if(verbose)
	{
		print(indir)
		print(infile)
		print(outdir)
		print(outfile)
		print(resume)
		print(method.node.stat)
		print(tree.id)
	}
	
	file		<- paste(indir,'/',infile, sep='')
	tmp			<- hivc.beast2out.read.nexus.and.stats(file, tree.id=tree.id, method.node.stat=method.node.stat)
	node.stat	<- tmp$node.stat
	tree		<- tmp$tree
	file		<- paste(outdir, '/', outfile, sep='')
	cat(paste('\nsave to file=', file))
	save(file=file, node.stat, tree)	
}
######################################################################################
hivc.prog.BEAST.generate.xml<- function()
{	
	require(XML)
	
	indir				<- paste(DATA,"tmp",sep='/')		
	infile				<- "ATHENA_2013_03_NoDRAll+LANL_Sequences"		
	insignat			<- "Thu_Aug_01_17/05/23_2013"
	indircov			<- paste(DATA,"derived",sep='/')
	infilecov			<- "ATHENA_2013_03_AllSeqPatientCovariates"
	infiletree			<- paste(infile,"examlbs100",sep="_")
	infilexml			<- paste(infile,'_',"beast",'_',"seroneg",sep='')
	
	outdir				<- indir
	outsignat			<- "Tue_Aug_26_09/13/47_2013"
		
	opt.fNegT				<- 0 #0.8
	opt.burnin				<- 2e7
	opt.brl					<- "dist.brl.casc" 
	thresh.brl				<- 0.096
	thresh.bs				<- 0.8
	pool.ntip				<- 130
	beast.mcmc.chainLength	<- 100000000
	infilexml.opt			<- "mph4cluLdTd"
	infilexml.template		<- "standard"	
	
	resume				<- 1
	verbose				<- 1
	hpc.walltime		<- 171
	hpc.ncpu			<- 4
	hpc.mem				<- "1800mb"

	
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
		if(length(tmp)>0) infilecov<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infiletree= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infiletree<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infilexml= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									outsignat= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) outsignat<- tmp[1]				
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
						{	switch(substr(arg,2,8),
									opt.brl= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) opt.brl<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									thresh.bs= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.bs<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									thresh.brl= return(as.numeric(substr(arg,13,nchar(arg)))),NA)	}))
		if(length(tmp)>0) thresh.brl<- tmp[1]
		#
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									pool.ntip= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) pool.ntip<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,14),
									infilexml.opt= return(substr(arg,16,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.opt<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,19),
									infilexml.template= return(substr(arg,21,nchar(arg))),NA)	}))
		if(length(tmp)>0) infilexml.template<- tmp[1]
	}	
	
	#	load complete tree to generate starting tree
	file					<- paste(indir,'/',infiletree,'_',gsub('/',':',insignat),".R",sep='')
	if(verbose)	cat(paste("\nload complete tree to generate starting tree from file",file))
	load(file)	#load object 'ph'
	
	xml.monophyly4clusters			<- ifelse(grepl("mph4clu",infilexml.opt),1,0)		
	pool.includealwaysbeforeyear	<- ifelse(grepl("fx03",infilexml.opt),2003, NA)
	xml.prior4tipstem				<- ifelse(grepl("u4tip",infilexml.opt),"uniform",NA)
	xml.taxon4tipstem				<- ifelse(	any(sapply(c("u4tip","tx4tip"), function(x) grepl(x,infilexml.opt))),	1, 0)
	df.resetTipDate					<- NA
	if(grepl("LdTd",infilexml.opt))				df.resetTipDate<- "LdTd"
	else if(grepl("LsTd",infilexml.opt))		df.resetTipDate<- "LsTd"
	else if(grepl("UmTd",infilexml.opt))		df.resetTipDate<- "UmTd"
	#case NoTd is dealt with below
	#
	#
	if(verbose)
	{
		print(indir)
		print(infile)
		print(insignat)
		print(indircov)
		print(infilecov)
		print(infiletree)
		print(infilexml)
		print(outdir)		
		print(outsignat)
		print(resume)
		print(opt.brl)
		print(thresh.brl)
		print(thresh.bs)
		print(pool.ntip)
		print(pool.includealwaysbeforeyear)
		print(infilexml.opt)
		print(infilexml.template)
		print(xml.monophyly4clusters)
		print(xml.taxon4tipstem)
		print(xml.prior4tipstem)
		print(df.resetTipDate)
	}	
	#
	#	load sequences
	#
	file		<- paste(indir,'/',infile,'_',gsub('/',':',insignat),".R",sep='')
	if(verbose)	cat(paste("\nload sequences from file",file))
	load( file )
	#print( seq.PROT.RT )	
	#
	#	load msm clusters
	#
	argv		<<- hivc.cmd.clustering.msm(indir, infiletree, insignat, indircov, infilecov, opt.brl, thresh.brl, thresh.bs, resume=resume)
	argv		<<- unlist(strsplit(argv,' '))		
	msm			<- hivc.prog.get.clustering.MSM()	
	#
	#	select seroconverters
	#
	df.cluinfo				<- msm$df.cluinfo
	tmp						<- df.cluinfo[,	list(clu.bwpat.medbrl=clu.bwpat.medbrl[1],clu.npat=clu.npat[1], clu.fPossAcute=clu.fPossAcute[1], fNegT=length(which(!is.na(NegT))) / clu.ntip[1]),by="cluster"]										
	tmp						<- subset(tmp, fNegT>=quantile(tmp[,fNegT], probs=opt.fNegT) )
	cluphy.df				<- merge( subset(tmp,select=cluster), df.cluinfo, all.x=1, by="cluster" )
	if(verbose) cat(paste("\nnumber of selected sequences is n=",nrow(cluphy.df)))
	cluphy.df				<- hivc.beast.addBEASTLabel( cluphy.df, df.resetTipDate=df.resetTipDate )
	#
	#	create sets of cluster pools for BEAST
	#
	df.clupool				<- hivc.beast.poolclusters(cluphy.df, pool.ntip= pool.ntip, pool.includealwaysbeforeyear=pool.includealwaysbeforeyear, verbose=1)			
	#
	if(0)		#used to generate standard xml file
	{
		outfile		<- paste(infile,"beast","seroneg",sep='_')
		outsignat	<- "Thu_Aug_01_17/05/23_2013"
		hivc.beast.writeNexus4Beauti(seq.PROT.RT, cluphy.df, file=paste(outdir,'/',outfile,'_',gsub('/',':',outsignat),".nex",sep=''))
	}
	#
	#	load xml template file
	#	
	file		<- paste(CODE.HOME,"/data/BEAST_template_",infilexml.template,".xml",sep='')
	btemplate	<- xmlTreeParse(file, useInternalNodes=TRUE, addFinalizer = TRUE)			
	#
	#	produce xml files for each cluster pool from template
	#	
	bfile		<- lapply(seq_len(length(df.clupool$pool.df)), function(pool.id)
			{
				df					<- df.clupool$pool.df[[pool.id]]
				setkey(df, cluster)
				#print( unique(df[,cluster]) )
				#	get xml file 
				outfile				<- paste(infilexml,'-',pool.ntip,'-',pool.id,'_',infilexml.template,'_',infilexml.opt,'_',gsub('/',':',outsignat),sep='')
				beast.label.datepos	<- ifelse(!is.na(df.resetTipDate), 5, 4)
				beast.usingDates	<- ifelse(grepl("NoTd",infilexml.opt), "false", "true")
				bxml				<- hivc.beast.get.xml(btemplate, seq.PROT.RT, df, outfile, ph=ph, xml.monophyly4clusters=xml.monophyly4clusters, xml.taxon4tipstem=xml.taxon4tipstem, xml.prior4tipstem=xml.prior4tipstem, beast.label.datepos=beast.label.datepos, beast.label.sep= '_', beast.date.direction= "forwards", beast.date.units= "years", beast.usingDates=beast.usingDates, beast.mcmc.chainLength=beast.mcmc.chainLength, verbose=1)
				getNodeSet(bxml, "//*[@id='tmrca(c1)']")		
				#	write xml file
				file				<- paste(outdir,'/',outfile,".xml",sep='')
				if(verbose)	cat(paste("\nwrite xml file to",file))
				saveXML(bxml, file=file)		
				#	freed through R garbage collection (I hope!) since addFinalizer=TRUE throughout
				paste(infilexml,'-',pool.ntip,'-',pool.id,'_',infilexml.template,'_',infilexml.opt,sep='')
			})
	#
	#	generate BEAST commands and run
	#
	sapply(seq_along(bfile), function(pool.id)
			{
				cmd			<- hivc.cmd.beast.runxml(outdir, bfile[[pool.id]], outsignat, hpc.tmpdir.prefix="beast", hpc.ncpu=hpc.ncpu)
				#cmd			<- paste(cmd, hivc.cmd.beast2.getclustertrees.pipe(outdir, infilexml, outsignat, infilexml.opt, infilexml.template, opt.burnin, outdir=outdir, outsignat=outsignat, opt.pool=pool.id, verbose=verbose, resume=resume), sep='' )
				cmd			<- paste(cmd,hivc.cmd.beast.evalrun(outdir, infilexml, outsignat, infilexml.opt, infilexml.template, length(bfile), verbose=1),sep='')				
				cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=hpc.walltime, hpc.q="pqeph", hpc.mem=hpc.mem,  hpc.nproc=hpc.ncpu)					
				cat(cmd)
				outfile		<- paste("bea",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
				hivc.cmd.hpccaller(outdir, outfile, cmd)
			})		
}
######################################################################################
hivc.prog.get.geneticdist<- function()
{
	library(bigmemory)
	library(ape)
	
	indir		<- outdir<- paste(DATA,"tmp",sep='/')
	infile		<- "ATHENA_2013_03_FirstAliSequences_PROTRT"
	resume		<- verbose <- 1	
	signat 		<- "Wed_May__1_17/08/15_2013"
	gd.max		<- NA
	out.phylip	<- 0
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
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))		
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									signat= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) signat<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									maxgd= return(as.numeric(substr(arg,8,nchar(arg)))),NA)	}))
		if(length(tmp)>0) gd.max<- tmp[1]
	}
	if(1)
	{
		print(indir)
		print(infile)
		print(outdir)
		print(signat)
		print(resume)
		print(gd.max)
	}	
	
	tmp			<- list.files(outdir, pattern=paste(".gdm$",sep=''))		
	file		<- tmp[ grepl(paste(infile,'_',sep=''), tmp) & grepl(gsub('/',':',signat), tmp) ]
	if(verbose)	cat(paste("\nFound .gdm files matching input args, n=", length(file)))	
	if(resume)
	{
		if(verbose)	cat(paste("\nloading",file))
		gd.bigmat			<- read.big.matrix(file, has.row.names=1, sep=',', type="char")		
	}
	if(!resume || !length(file))
	{		
		if(verbose)	cat(paste("\ncreate gdm file"))
		file				<- paste(indir,"/",infile,"_",gsub('/',':',signat),".R",sep='')
		if(verbose)	cat(paste("\nload",file))
		load(file)
		str(seq.PROT.RT)		
			#tmp				<- tmp[1:10,]
		gd.bigmat			<- seq.dist(  seq.PROT.RT )
		file				<- paste(outdir,"/",infile,"_",gsub('/',':',signat),".gdm",sep='')
		if(verbose) cat(paste("\nwrite to",file))
		write.big.matrix(gd.bigmat, file, row.names= 1, col.names=0, sep=',')		
		#gd.bigmat.d 		<- describe( gd.bigmat )
		#file				<- paste(outdir,"/ATHENA_2013_03_FirstAliSequences_PROTRT_",gsub('/',':',signat),".gdm",sep='')		
		#dput(gd.bigmat.d, file=file)
	}
	if(!is.na(gd.max))
	{
		if(verbose)	cat(paste("\nselect sequences with at least one other sequence within gd.max",gd.max))
		#	now have pairwise genetic distances in gd.bigmat
		gd.bigmat.min	<- sapply(seq.int(1,nrow(gd.bigmat)-1),function(i)
								{
									tmp<- which.min(gd.bigmat[i,])
									c(i,tmp, gd.bigmat[i,tmp])
								})
		rownames(gd.bigmat.min)<- c("seq.idx1","seq.idx2","gd")
		gd.bigmat.seq.i	<- which( gd.bigmat.min["gd",] < (gd.max * 1e3) )
		gd.seqs			<- gd.bigmat.min[c("seq.idx1","seq.idx2"),gd.bigmat.seq.i]
		gd.seqs			<- sort(unique(as.vector(gd.seqs)))
		if(verbose) cat(paste("\ncomputed sequences sth at least one pair with less than ",gd.max*100,"% genetic distance, n=",length(gd.seqs)))		
		file			<- paste(indir,"/",infile,"_",gsub('/',':',signat),".R",sep='')
		if(verbose)	cat(paste("\nload",file))
		load(file)
		seq.PROT.RT.gd	<- seq.PROT.RT[ rownames(gd.bigmat)[gd.seqs], ]
		#	save selected sequences in R
		file 			<- paste(outdir,"/",infile,"_Gd",gd.max*1000,"_",gsub('/',':',signat),".R",sep='')		
		save(seq.PROT.RT.gd, file=file)		
		#	save selected sequences in phylip
		if(out.phylip)
		{
			file 			<- paste(outdir,"/",infile,"_Gd",gd.max*1000,"_",gsub('/',':',signat),".phylip",sep='')
			if(verbose) cat(paste("\nwrite to ",file))
			seq.write.dna.phylip(seq.PROT.RT.gd, file=file)
		}
	}
	else
	{
		seq.PROT.RT.gd<- NULL		
	}
	list(gd.bigmat=gd.bigmat, seq.PROT.RT.gd=seq.PROT.RT.gd)						
}
