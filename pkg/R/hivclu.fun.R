#' this file contains all R functions of the hivclust package
#' @useDynLib hivc

#' @export
hivc.seq.read.GenBank<- function (access.nb, seq.names = access.nb, species.names = TRUE, gene.names = FALSE, as.character = FALSE, attributes= c("origin")) 
{
	require(ape)
	N <- length(access.nb)
	nrequest <- N%/%400 + as.logical(N%%400)
	X <- character(0)
	for (i in 1:nrequest) {
			a <- (i - 1) * 400 + 1
			b <- 400 * i
			if (i == nrequest) 
				b <- N
			URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", 
					paste(access.nb[a:b], collapse = ","), "&rettype=gb&retmode=text", 
					sep = "")
			X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
		}
		
	FI <- grep("^ {0,}ORIGIN", X) + 1
	LA <- which(X == "//") - 1
	obj <- vector("list", N)
	for (i in 1:N) {
		tmp <- gsub("[[:digit:] ]", "", X[FI[i]:LA[i]])
		obj[[i]] <- unlist(strsplit(tmp, NULL))
	}
	names(obj) <- seq.names
	if (!as.character) 
		obj <- as.DNAbin(obj)
	if(length(attributes)) 
	{
		attr.lines	<- lapply(attributes, function(attr) grep(paste("^ {0,}/",attr,sep=''), X) 	)
		attr.fields	<- lapply(seq_along(attr.lines),function(i)		gsub("\"$","",gsub(paste("^ {0,}/",attributes[i],"=\"",sep=''),"",X[attr.lines[[i]]]))		)						
		for(i in seq_along(attributes))
			attr(obj, attributes[[i]])<- attr.fields[[i]]
	}
	if (gene.names) {
		tmp <- character(N)
		sp <- grep(" +gene +<", X)
		for (i in 1:N) tmp[i] <- unlist(strsplit(X[sp[i + 1L]], 
							" +/gene=\""))[2]
		attr(obj, "gene") <- gsub("\"$", "", tmp)
	}
	obj
}

#' @export
or.dist.dna<- function (x, model = "K80", variance = FALSE, gamma = FALSE, 
		pairwise.deletion = FALSE, base.freq = NULL, as.matrix = FALSE) 
{
	MODELS <- c("RAW", "JC69", "K80", "F81", "K81", "F84", "T92", 
			"TN93", "GG95", "LOGDET", "BH87", "PARALIN", "N", "TS", 
			"TV", "INDEL", "INDELBLOCK")
	imod <- pmatch(toupper(model), MODELS)
	if (is.na(imod)) 
		stop(paste("'model' must be one of:", paste("\"", MODELS, 
								"\"", sep = "", collapse = " ")))
	if (imod == 11 && variance) {
		warning("computing variance not available for model BH87")
		variance <- FALSE
	}
	if (gamma && imod %in% c(1, 5:7, 9:17)) {
		warning(paste("gamma-correction not available for model", 
						model))
		gamma <- FALSE
	}
	if (is.list(x)) 
		x <- as.matrix(x)
	nms <- dimnames(x)[[1]]
	n <- dim(x)
	s <- n[2]
	n <- n[1]
	if (imod %in% c(4, 6:8)) {
		BF <- if (is.null(base.freq)) 
					base.freq(x)
				else base.freq
	}
	else BF <- 0
	if (imod %in% 16:17) 
		pairwise.deletion <- TRUE
	if (!pairwise.deletion) {
		keep <- .C("GlobalDeletionDNA", x, n, s, rep(1L, s), 
				PACKAGE = "ape")[[4]]
		x <- x[, as.logical(keep)]
		s <- dim(x)[2]
	}
	Ndist <- if (imod == 11) 
				n * n
			else n * (n - 1)/2
	var <- if (variance) 
				double(Ndist)
			else 0
	if (!gamma) 
		gamma <- alpha <- 0
	else alpha <- gamma <- 1
print(x)	
	d <- .C("dist_dna", x, n, s, imod, double(Ndist), BF, as.integer(pairwise.deletion), 
			as.integer(variance), var, as.integer(gamma), alpha, 
			DUP = FALSE, NAOK = TRUE, PACKAGE = "ape")
	if (variance) 
		var <- d[[9]]
	d <- d[[5]]
	if (imod == 11) {
		dim(d) <- c(n, n)
		dimnames(d) <- list(nms, nms)
	}
	else {
		attr(d, "Size") <- n
		attr(d, "Labels") <- nms
		attr(d, "Diag") <- attr(d, "Upper") <- FALSE
		attr(d, "call") <- match.call()
		attr(d, "method") <- model
		class(d) <- "dist"
		if (as.matrix) 
			d <- as.matrix(d)
	}
	if (variance) 
		attr(d, "variance") <- var
	d
}

#' @export
hivc.seq.create.referencepairs<- function(dir.name= DATA)
{
	if(0)	#generate ATHENA_2013_hptn052.rda
	{
		gban				<- c( paste("JN",seq.int(247047,247075),sep=''),paste("JN",seq.int(634296,634492),sep='') )	
		tmp					<- hivc.read.GenBank(gban, as.character=0, attributes= c("isolate","country","collection_date"))		
		file				<- "ATHENA_2013_hptn052.fa"
		write.dna(tmp, paste(DATA,file,sep='/'), format = "fasta" )
		cmd					<- hiv.cmd.clustalo(paste(dir.name,"tmp",sep='/'), file, signat='', outdir=paste(dir.name,"tmp",sep='/'))
		cat(cmd)
		if(0) system(cmd)		
		file				<- paste(DATA,"tmp/ATHENA_2013_hptn052.fa.clustalo",sep='/')
		hptn052				<- read.dna(file, format = "fasta" )		
		rownames(hptn052)	<- attr(tmp,"isolate")		
		file				<- paste(DATA,"tmp/ATHENA_2013_hptn052.rda",sep='/')
		save(hptn052, file= file)
	}
	file<- paste(DATA,"tmp/ATHENA_2013_hptn052.rda",sep='/')		
	load(file)
	#select control sequences and compute pairwise distances between them
	hptn052.control		<- grep("control", rownames(hptn052))	
	hptn052.control.d	<- hivc.pwdist( hptn052[hptn052.control,] )
	hptn052.control.d	<- hptn052.control.d[ upper.tri(hptn052.control.d) ]
	#select positive sequences and compute pairwise distance between them	
	hptn052.sdc			<- c(grep("A", rownames(hptn052)), grep("B", rownames(hptn052)))
	hptn052.sdc			<- hptn052[hptn052.sdc,]
	hptn052.sdc.ind		<- sapply(strsplit(rownames(hptn052.sdc),'-'), function(x)
			{ 
				as.numeric(x[2]) + ifelse(x[3]=='I',0,0.5)
			})
	hptn052.sdc.ind		<- sapply( unique( hptn052.sdc.ind ), function(x){		hptn052.sdc[hptn052.sdc.ind==x,]	})
	hptn052.sdc.ind		<- lapply( which( sapply(hptn052.sdc.ind,nrow)==2 ), function(i)  hptn052.sdc.ind[[i]] )
	hptn052.sdc.ind.d	<- sapply( hptn052.sdc.ind, function(x) hivc.pwdist(x)[1,2] )
	
	list( gd.islink= hptn052.sdc.ind.d, gd.unlinked= hptn052.control.d )
}

#' @export
hivc.seq.write.dna.phylip<- function(seq.DNAbin.mat, file)
{		
	tmp<- cbind( rownames(seq.DNAbin.mat), apply( as.character( seq.DNAbin.mat ), 1, function(x) paste(x,sep='',collapse='')  ) )
	tmp<- paste(t(tmp),collapse='\n',sep='')	
	tmp<- paste( paste(c(nrow(seq.DNAbin.mat),ncol(seq.DNAbin.mat)),sep='',collapse=' '),'\n',tmp,'\n',collapse='',sep='' )
	cat(tmp, file=file)
}


#' @export
hivc.seq.dist<- function(seq.DNAbin.matrix, verbose=1)
{
	if(0)
	{
		require(ape)
		ans<- dist.dna(seq.DNAbin.matrix, model="raw", as.matrix=1)
	}
	if(1)
	{
		library(bigmemory)
		options(bigmemory.typecast.warning=FALSE)
		big.matrix.charmax<- 127
		dummy	<- 0
		ans<- big.matrix(nrow(seq.DNAbin.matrix),nrow(seq.DNAbin.matrix), dimnames=list(rownames(seq.DNAbin.matrix),c()), type="char", init=NA)		
		ans[1,1:6]<- 2^seq.int(6,11)-1
		if(is.na(ans[1,2]))		stop("unexpected behaviour of bigmemory")
		ans[1,1:6]<- NA
								
		for(i1 in seq.int(1,nrow(seq.DNAbin.matrix)-1))
		{			
			seq1<- seq.DNAbin.matrix[i1,]
			time<- system.time	(
							tmp	<- 1 - sapply(seq.int(i1+1,nrow(seq.DNAbin.matrix)),function(i2){		.C("hivc_dist_ambiguous_dna", seq1, seq.DNAbin.matrix[i2,], ncol(seq1), dummy )[[4]]			})
						)[3]
			if(verbose)	cat(paste("\ncompute distance of row",i1,"entries",nrow(seq.DNAbin.matrix)-i1,"took",time))			
			tmp												<- round(tmp*1e3,d=0)			
			tmp[tmp>big.matrix.charmax]						<- big.matrix.charmax
			ans[i1, seq.int(i1+1,nrow(seq.DNAbin.matrix))]	<- tmp
		}		
	}
	ans
}

#' @export
hivc.seq.rmgaps<- function(seq.DNAbin.matrix, verbose=0)
{
	seq.DNAbin.matrix	<- as.character(seq.DNAbin.matrix)
	nogap				<- which( !apply(seq.DNAbin.matrix,2,function(x) all(x=="-" || x=="?")) )
	if(verbose)	cat(paste("\nremove gaps, n=",ncol(seq.DNAbin.matrix)-length(nogap)))
	as.DNAbin( seq.DNAbin.matrix[,nogap] )
}

#' @export
hivc.clu.geneticdist.cutoff<- function(dir.name= DATA, plot=1, verbose=1, level.retain.unlinked=0.05)
{	
	refpairs			<- hivc.seq.create.referencepairs(dir.name)
	refpairs			<- lapply(refpairs,function(x) x*100)
	
	xlim				<- range(c(refpairs[[1]], refpairs[[2]]))
	width				<- c(1.5,0.75)
	breaks				<- seq(xlim[1],xlim[2],len=40)	
	refpairs.h			<- lapply(seq_along(refpairs),function(i)
			{ 				
				tmp<- hist(refpairs[[i]],breaks=breaks,plot=0)
				tmp$kde<- density(refpairs[[i]], kernel="biweight",from=breaks[1],to=breaks[length(breaks)],width = max(EPS,width[i]*diff(summary(refpairs[[i]])[c(2,5)])))
				tmp
			})
	names(refpairs.h)	<- names(refpairs)
	ylim				<- range(c(refpairs.h[[1]]$intensities,refpairs.h[[2]]$intensities))
	cols				<- c(my.fade.col("black",0.2),my.fade.col("black",0.6))
	plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab="% genetic distance",ylab="density")		
	lapply(seq_along(refpairs.h),function(i)
			{
				plot(refpairs.h[[i]],add=1,freq=0,col=cols[i],border=NA)				
			})
	lapply(seq_along(refpairs.h),function(i)
			{
				lines(refpairs.h[[i]][["kde"]]$x,refpairs.h[[i]][["kde"]]$y,col=cols[i])								
			})
	
	ref.gd.unlinked.cdf	<- cumsum( refpairs.h[["gd.unlinked"]][["kde"]]$y / sum(refpairs.h[["gd.unlinked"]][["kde"]]$y) )		
	gd.cutoff			<- refpairs.h[["gd.unlinked"]][["kde"]]$x[ which.max( ref.gd.unlinked.cdf[ref.gd.unlinked.cdf<=level.retain.unlinked] ) ]
	abline(v=gd.cutoff,lty=2)
	gd.specificity<- refpairs.h[["gd.islink"]][["kde"]]$y / sum(refpairs.h[["gd.islink"]][["kde"]]$y)
	gd.specificity<- sum( gd.specificity[ refpairs.h[["gd.islink"]][["kde"]]$x<=gd.cutoff ] )
	
	if(verbose)	cat(paste("\ncutoff for high sensitivity of",1-level.retain.unlinked,"of detecting unlinked is",gd.cutoff))
	if(verbose)	cat(paste("\nspecificity for detecting linked at cutoff",gd.cutoff," is",gd.specificity))
	
	list(gd.cutoff= gd.cutoff, gd.specificity=gd.specificity)	
}
