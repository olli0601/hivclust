PR.CLUSTALO	<<- "clustalo"
HPC.NPROC	<<- 1
HPC.SYS		<<- "CX1"
HPC.LOAD	<<- "module load intel-suite/10.0 R/2.13.0"

#generate clustalo command
hiv.cmd.clustalo<- function(indir, infiles, signat=paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''), outdir=indir, prog= PR.CLUSTALO, nproc=HPC.NPROC)
{
	lapply(infiles,function(x)
			{
				cat(paste("\nprocess",x,"\n"))				
				#verbose stuff
				cmd<- "#######################################################
# run clustalo
#######################################################"
				cmd<- paste(cmd,paste("\necho \'run ",PR.CLUSTALO,"\'\n",sep=''))
				
				#default commands
				cmd<- paste(cmd,PR.CLUSTALO,sep=" ")								
				cmd<- paste(cmd, " --infmt=fa --outfmt=fa --force --threads=",nproc," ", sep='' )
				#file in/out
				tmp<- paste(indir,x,sep='/')
				cmd<- paste(cmd, paste("--in",tmp,sep=' ') )				
				tmp<- paste(outdir,paste(x,signat,"clustalo",sep='.'),sep='/')
				cmd<- paste(cmd, paste("--out",tmp,sep=' ') )
				
				#verbose stuff
				cmd<- paste(cmd,paste("\necho \'end ",PR.CLUSTALO,"\'\n\n",sep=''))
				cmd
			})
}

#add additional high performance computing information 
hiv.cmd.hpcwrapper<- function(cmd, hpc.sys=HPC.SYS, hpc.walltime=24, hpc.select= "1:ncpus=1:mem=400mb", hpc.q=NA)
{
	wrap<- "#!/bin/sh"
	if(hpc.sys=="CX1")
	{				
		tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
		wrap<- paste(wrap, tmp, sep='\n')
		tmp	<- paste("#PBS -l select=",hpc.select,sep='')
		wrap<- paste(wrap, tmp, sep='\n')
		wrap<- paste(wrap, "#PBS -j oe", sep='\n')
		if(!is.na(hpc.q))
			wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
		wrap<- paste(wrap, HPC.LOAD, sep='\n\n')

	}
	else
		stop(paste("unknown hpc.sys",hpc.sys))
	
	cmd<- lapply(seq_along(cmd),function(i){	paste(wrap,cmd[[i]],sep='\n')	}) 
	cmd	
}

#create high performance computing qsub file and submit
hiv.cmd.hpccaller<- function(outdir, outfile, cmd)
{
	file<- paste(outdir,outfile,sep='/')
	cat(paste("\nwrite cmd to",file))
	cat(cmd,file=file)
	cmd<- paste("qsub",file)
	cat( cmd )
	cat( system(cmd, intern=TRUE) )	
}