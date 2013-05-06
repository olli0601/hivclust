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
PR.FIRSTSEQ		<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeFIRSTSEQ",sep='/')

#' @export
PR.GENDISTMAT	<- paste(HIVC.CODE.HOME,"pkg/misc/hivclu.startme.R -exeGENDISTMAT",sep='/')

#' @export
PR.EXAML.PARSER	<- "ExaML-parser"

#' @export
PR.EXAML.STARTTREE	<- "ExaML-parsimonator"

#' @export
PR.EXAML.EXAML	<- "ExaML-examl"

#' @export
HPC.NPROC		<- {tmp<- c(1,6); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}

#' @export
HPC.MPIRUN		<- {tmp<- c("mpirun","mpiexec"); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}

#' @export
HPC.CX1.IMPERIAL<- "cx1.hpc.ic.ac.uk"		#this is set to system('domainname',intern=T) for the hpc cluster of choice

#' @export
HPC.MEM			<- "3600mb"

#' @export
HPC.LOAD		<- "module load intel-suite/10.0 mpi R/2.13.0"

#generate clustalo command
#' @export
hivc.cmd.clustalo<- function(indir, infiles, signat=paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''), outdir=indir, prog= PR.CLUSTALO, hmm=PR.CLUSTALO.HMM, nproc=1)
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

#' @export
hivc.cmd.get.geneticdist<- function(indir, signat, outdir=indir, prog= PR.GENDISTMAT, resume=1, verbose=1)
{
	cmd<- "#######################################################
# compute geneticdist
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -outdir=",outdir," -signat=",signat,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n\n",sep=''))
	cmd
}

#' @export
hivc.cmd.get.firstseq<- function(indir, infile, signat.in, signat.out, outdir=indir, prog= PR.FIRSTSEQ, resume=1, verbose=1)
{
	cmd<- "#######################################################
# extract firstseq
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog,"\'\n",sep=''))
	#default commands
	cmd<- paste(cmd,prog," -v=",verbose," -resume=",resume,sep='')
	cmd<- paste(cmd," -indir=",indir," -infile=",infile," -outdir=",outdir," -insignat=",signat.in," -outsignat=",signat.out,sep='')
	#verbose stuff
	cmd<- paste(cmd,paste("\necho \'end ",prog,"\'\n\n",sep=''))
	cmd
}

#' @export
hivc.cmd.examl<- function(indir, infile, signat.in, signat.out, outdir=indir, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, args.starttree="-p 12345", prog.examl= PR.EXAML.EXAML, args.examl="-m GAMMA -B 20 -D", resume=1, verbose=1)
{
	cmd<- "#######################################################
# compute ExaML tree
#######################################################"
	cmd<- paste(cmd,paste("\necho \'run ",prog.parser,"\'\n",sep=''))
	#default commands for parser	
				
	curr.dir	<- getwd()
	cmd			<- paste(cmd,"cd ",outdir,'\n',sep='')
	tmp			<- paste(indir,paste(infile,'_',signat.in,".phylip",sep=''),sep='/')
	cmd			<- paste(cmd,prog.parser," -m DNA -s ",tmp,sep='')
	tmp			<- paste(infile,'_',signat.out,".phylip.examl",sep='')
	cmd			<- paste(cmd," -n ",tmp,sep='')
	#verbose stuff for parser	
	cmd			<- paste(cmd,paste("\necho \'end ",prog.parser,"\'",sep=''))
		
	cmd			<- paste(cmd,paste("\necho \'run ",prog.starttree,"\'\n",sep=''))
	#default commands for parser
	tmp			<- paste(indir,paste(infile,'_',signat.in,".phylip",sep=''),sep='/')
	cmd			<- paste(cmd,prog.starttree,' ',args.starttree," -s ",tmp,sep='')	
	tmp			<- paste(infile,'_',signat.out,".examlstarttree",sep='')
	cmd			<- paste(cmd," -n ",tmp,sep='')
	#verbose stuff
	cmd			<- paste(cmd,paste("\necho \'end ",prog.starttree,"\'",sep=''))
	
	cmd			<- paste(cmd,paste("\necho \'run ",prog.examl,"\'\n",sep=''))
	#default commands for parser
	tmp			<- hivc.get.hpcsys()
	if(tmp=="debug")
		cmd		<- paste(cmd,HPC.MPIRUN[tmp]," -np ",HPC.NPROC[tmp],' ',prog.examl,' ',args.examl,sep='')
	else if(tmp==HPC.CX1.IMPERIAL)
		cmd		<- paste(cmd,HPC.MPIRUN[tmp],' ',prog.examl,' ',args.examl,sep='')
	else	
		stop("unknown hpc sys")
	tmp			<- paste(infile,'_',signat.out,".phylip.examl.binary",sep='')
	cmd			<- paste(cmd," -s ",tmp,sep='')
	tmp			<- paste("RAxML_parsimonyTree.",infile,'_',signat.out,".examlstarttree.0",sep='')
	cmd			<- paste(cmd," -t ",tmp,sep='')
	tmp			<- paste(infile,'_',signat.out,".examltree",sep='')
	cmd			<- paste(cmd," -n ",tmp,sep='')		
	
	cmd			<- paste(cmd,paste("\necho \'end ",prog.examl,"\'\n",sep=''))	
	cmd			<- paste(cmd,"cd ",curr.dir,'\n',sep='')
	
	cmd
}

#' @export
hivc.cmd.examl.cleanup<- function(outdir, prog= PR.EXAML.EXAML)
{
	cmd<- "#######################################################
# clean up after ExaML tree
#######################################################"
	cmd<- paste(cmd,paste("\necho \'clean after ",prog,"\'\n",sep=''))	
	
	tmp<- list.files(outdir, full.names=1)
	#tmp<- tmp[c(grep("examlstarttree",tmp), grep("examlbin",tmp),grep("examl.binary",tmp),grep("phylip",tmp))]
	#cmd<- paste(cmd, "\nrm ", paste(tmp,collapse=' ',sep=''), '\n', sep='')	
	cmd<- paste(cmd,"rm ",paste(paste(outdir,c("RAxML_info*","RAxML_parsimonyTree*","ExaML_binaryCheckpoint*"),sep='/'),collapse=' ',sep=''),sep=' ')
	
	cmd<- paste(cmd,paste("\necho \'cleaned up after ",prog,"\'\n",sep=''))
	cmd
}

#add additional high performance computing information 
#' @export
hivc.cmd.hpcwrapper<- function(cmd, domain.name= hivc.get.hpcsys(), hpc.walltime=24, hpc.mem=HPC.MEM, hpc.nproc=HPC.NPROC[domain.name], hpc.q=NA)
{
	wrap<- "#!/bin/sh"	
	if(1 || domain.name==HPC.CX1.IMPERIAL)
	{				
		tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
		wrap<- paste(wrap, tmp, sep='\n')		
		tmp	<- paste("#PBS -l select=1:ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
		wrap<- paste(wrap, tmp, sep='\n')
		wrap<- paste(wrap, "#PBS -j oe", sep='\n')
		if(!is.na(hpc.q))
			wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
		wrap<- paste(wrap, HPC.LOAD, sep='\n')
		#tmp	<- paste("export PATH=$PATH:",HPC.BIN,sep='')
		#wrap<- paste(wrap, tmp, sep='\n')

	}
	else if(tmp=='')
		cat(paste("\ndetected no HPC system and no hpcwrapper generated, domain name is",tmp))
	else
		stop(paste("unknown hpc system with domain name",tmp))
	
	cmd<- lapply(seq_along(cmd),function(i){	paste(wrap,cmd[[i]],sep='\n')	})
	if(length(cmd)==1)
		cmd<- unlist(cmd)
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