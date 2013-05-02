if(!exists("HIVC.CODE.HOME"))	
{	
	HIVC.CODE.HOME	<- getwd()
	INST			<- paste(HIVC.CODE.HOME,"inst",sep='/')
}

#' @export
PR.CLUSTALO		<- "clustalo"

#' @export
PR.CLUSTALO.HMM	<- paste(INST,"align_HIV-1_pol_DNA.hmm",sep='/')

#' @export
PR.GENDISTMAT	<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeGENDISTMAT",sep='/')

#' @export
HPC.NPROC		<- 1

#' @export
HPC.SYS			<- "CX1"

#' @export
HPC.LOAD		<- "module load intel-suite/10.0 R/2.13.0"

#generate clustalo command
#' @export
hivc.cmd.clustalo<- function(indir, infiles, signat=paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''), outdir=indir, prog= PR.CLUSTALO, hmm=PR.CLUSTALO.HMM, nproc=HPC.NPROC)
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
				cmd<- paste(cmd, " --infmt=fa --outfmt=fa --force", sep='')
				cmd<- paste(cmd, " --threads=",nproc," ", sep='' )
				#cmd<- paste(cmd, " --hmm-in=",hmm," ", sep='' )		not supported in clustalo v1.1.0
				
				#file in/out
				tmp<- paste(indir,x,sep='/')
				cmd<- paste(cmd, paste("--in",tmp,sep=' ') )				
				tmp<- paste(outdir,paste(x,ifelse(nchar(signat),'.',''),signat,".clustalo",sep=''),sep='/')
				cmd<- paste(cmd, paste("--out",tmp,sep=' ') )
				
				#verbose stuff
				cmd<- paste(cmd,paste("\necho \'end ",PR.CLUSTALO,"\'\n\n",sep=''))
				cmd
			})
}

hivc.cmd.get.geneticdist<- function(indir, signat, outdir=indir, prog= PR.GENDISTMAT)
{
	cmd<- "#######################################################
# run geneticdist
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",PR.GENDISTMAT,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,PR.GENDISTMAT,"-v=1 -resume=1",sep=" ")
	cmd<- paste(cmd," -indir=",indir," -outdir=",outdir," -signat=",signat,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",PR.GENDISTMAT,"\'\n\n",sep=''))
	cmd
}
	
#add additional high performance computing information 
#' @export
hivc.cmd.hpcwrapper<- function(cmd, hpc.sys=HPC.SYS, hpc.walltime=24, hpc.select= "1:ncpus=1:mem=400mb", hpc.q=NA)
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
#' @export
hivc.cmd.hpccaller<- function(outdir, outfile, cmd)
{
	file<- paste(outdir,outfile,sep='/')
	cat(paste("\nwrite cmd to",file,"\n"))
	cat(cmd,file=file)
	cmd<- paste("qsub",file)
	cat( cmd )
	cat( system(cmd, intern=TRUE) )	
}