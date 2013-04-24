project.hivc.clustalo<- function(dir.name= DATA, min.seq.len=21)
{		
	if(0)
	{
		#read csv data file and preprocess		
		verbose		<- 1
		analyze		<- 1
		
		file		<- paste(dir.name,"derived/ATHENA_2013_03_Sequences.csv",sep='/')
		df			<- read.csv(file, stringsAsFactors=FALSE)		
		df.names	<- paste("GeneCode",unique(df[,"GeneCode"]),sep='')
		df			<- lapply(unique(df[,"GeneCode"]),function(gene)
						{
							cat(paste("\nprocess GeneCode", gene))
							tmp<- df[ df[,"GeneCode"]==gene, c("SampleCode","Sequence"), drop=0 ]
							cat(paste("\nfound n=", nrow(tmp)))
							seq.len		<- nchar(tmp[,"Sequence"])
							seq.nok.idx	<- which(seq.len<min.seq.len)
							seq.ok.idx	<- which(seq.len>=min.seq.len)
							if(verbose)	cat(paste("\ndiscarding sequences with insufficient length, n=",length(seq.nok.idx),'\n'))
							if(verbose)	cat(paste("\ndiscarding sequence with insufficient length, SampleCode",paste(tmp[seq.nok.idx,"SampleCode"],collapse=', ')))
							tmp			<- tmp[seq.ok.idx,]
							if(verbose) cat(paste("\nfinal n=", nrow(tmp)))
							tmp
						})		
		names(df)	<- df.names
		file		<- paste(dir.name,"tmp/ATHENA_2013_03_Sequences.R",sep='/')		
		save(df, file=file)
	}
	if(0)
	{
		#generate fasta files per gene
		verbose<- 1
				
		file<- paste(dir.name,"tmp/ATHENA_2013_03_Sequences.R",sep='/')
		load(file)		
		lapply(seq_along(df),function(i)
				{
					df.gene<- df[[i]]
					df.gene[,"SampleCode"]	<- paste('>',df.gene[,"SampleCode"],sep='')
					seq.len		<- nchar(df.gene[,"Sequence"])	
					if(analyze)
					{
						cat(paste("\nnumber of sequences:",nrow(df.gene),"\nnumber of unique SampleCodes:", length(unique(df.gene[,"SampleCode"])),"\n") )
						hist(seq.len, breaks=100)	
						if(nrow(df.gene)!=length(unique(df.gene[,"SampleCode"])))
							o<- sapply(unique(df.gene[,"SampleCode"]),function(x)
									{
										tmp<- which(df.gene[,"SampleCode"]==x)
										if(length(tmp)>1)
										{ 
											cat("\n found duplicates\n")
											print(df.gene[tmp,])	
										}
									})
						seq.nok.idx	<- which(seq.len<min.seq.len)
						if(verbose)	cat(paste("discarding sequences with insufficient length, n=",length(seq.nok.idx),'\n'))			
					}
					tmp			<- as.vector(t(df.gene))
					tmp			<- paste(tmp,collapse='\n')					
					file		<- paste("derived/ATHENA_2013_03_Sequences_",names(df)[i],".fa",sep='')
					cat(paste("\nwrite file to",file))
					cat(tmp,file=file)	
				})			
	}
	if(1)
	{
		#generate clustalo command		 
		indir	<- paste(dir.name,"derived",sep='/')
		outdir	<- paste(dir.name,"tmp",sep='/')
		infiles	<- list.files(path=indir, pattern=".fa$", full.names=0)		
		signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
		cmd		<- hiv.cmd.clustalo(indir, infiles, signat=signat, outdir=outdir)
		cmd		<- hiv.cmd.hpcwrapper(cmd)
		
		lapply(cmd, cat)
		stop()
		lapply(cmd,function(x)
				{
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("clustalo",signat,"qsub",sep='.')					
					hiv.cmd.hpccaller(outdir, outfile, x)			
				})
		
		
	}
}