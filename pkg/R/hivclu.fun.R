#' this file contains all R functions of the hivclust package
#' @useDynLib hivc
#' @import ape
#' @import phytools
#' @import igraph
#' @import geiger
#' @import data.table

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

hivc.seq.find<- function(char.matrix, pos0= NA, from= c(), verbose=1)
{
	if(is.na(pos0)) 	stop("start position of token to be replaced is missing")
	if(!length(from))	stop("token to be replaced is missing")
	query.colidx	<- seq.int(pos0,pos0+length(from)-1)
	query.yes		<- which( apply(char.matrix, 1, function(x)	all(x[query.colidx]==from) ) )
	query.yes	
}

hivc.seq.gc.content<- function(seq.DNAbin.mat)
{	
	rna.gc.fraction.n		<- c('a','c','g','t',	'r','m','w','s',	'k','y','v','h',		'd','b','n','-','?')
	rna.gc.fraction			<- c( 0, 1, 1, 0,		0.5, 0.5, 0, 1, 	1/2, 1/2, 2/3, 1/3,		1/3,2/3, 1/4, 0, 0)		#this fraction assumes that U is synonymous with T and that U does not occur in the code
	rna.gc.sum				<- c( T, T, T, T,       T, T, T, T,         T, T, T, T,				T, T, T, F, F )	
	counts					<- apply(seq.DNAbin.mat,1,function(x) base.freq(x, freq=1, all=1))
	apply(counts*rna.gc.fraction,2,sum) / apply(counts[rna.gc.sum,],2,sum)
}

#slight modification of blastSequences() in pkg annotate
hivc.seq.blast<- function (x, database = "nr", hitListSize = "10", filter = "L", expect = "10", program = "blastn", organism= "HIV-1") 
{
	baseUrl <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
	query <- paste("QUERY=", as.character(x), "&DATABASE=", database, "&ORGANISM=", organism,
			"&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=", 
			expect, "&PROGRAM=", program, sep = "")
	url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
	results <- tempfile()
	Sys.sleep(5)
	require(XML)
	post <- htmlTreeParse(url0, useInternalNodes = TRUE)
	x <- post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
	rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
	rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", 
					x))
	url1 <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, 
			rid)
	Sys.sleep(rtoe)
	result <- annotate:::.tryParseResult(url1)
	qseq <- xpathApply(result, "//Hsp_qseq", xmlValue)
	hseq <- xpathApply(result, "//Hsp_hseq", xmlValue)
	require(Biostrings)
	res <- list()
	for (i in seq_len(length(qseq))) {
		res[i] <- DNAMultipleAlignment(c(hseq[[i]], qseq[[i]]), 
				rowmask = as(IRanges(), "NormalIRanges"), colmask = as(IRanges(), 
						"NormalIRanges"))
	}
	res
}

#slight modification of read.blast() in pkg RFLPtools; expects blast was run with -outfmt 6
hivc.seq.blast.read<- function (file, sep = "\t") 
{
	require(data.table)
	x <- read.table(file = file, header = FALSE, sep = sep, quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
	if (ncol(x) != 12) 
		stop("Data in given file", basename(file), "has wrong dimension!")
	names(x) <- c("query.id", "subject.id", "identity", "alignment.length","mismatches", "gap.opens", "q.start", "q.end", "s.start","s.end", "evalue", "bit.score")
	data.table(x, key="query.id")
}

#' @export
hivc.seq.rm.drugresistance<- function(char.matrix, dr, verbose=1, rtn.DNAbin=1)
{
	if(verbose)	cat(paste("\nchecking for drug resistance mutations, n=",nrow(dr)))
	tmp	<- rep(0, nrow(dr))
	for(i in seq_len(nrow(dr)))
	{		
		query.yes	<- hivc.seq.find(char.matrix, dr[i,Alignment.nuc.pos], unlist(strsplit(unlist(dr[i,Mutant.NTs]),'')))
		if(length(query.yes))
		{
			if(verbose)	
				cat(paste("\nfound DR mutation",dr[i,DR.name],"NT code",dr[i,Mutant.NTs],"at pos",dr[i,Alignment.nuc.pos],"n=",length(query.yes), ". Replacing NT seq with nnn."))			
			char.matrix[query.yes,	seq.int(dr[i,Alignment.nuc.pos], length.out=3) ]<- matrix("n", nrow=length(query.yes), ncol=3)			
			#print( char.matrix[query.yes, seq.int(dr[i,Alignment.nuc.pos]-3, length.out=9)] ); stop()
		}	
		tmp[i]		<- length(query.yes)
	}
	if(verbose)
		cat(paste("\nremoved DR mutations, n=",sum(tmp), "n per patient=",sum(tmp)/nrow(char.matrix)))
	if(rtn.DNAbin)
		return( as.DNAbin(char.matrix) )
	else
		return( char.matrix )
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

hivc.seq.replace<- function(seq.DNAbin.matrix, code.from='?', code.to='n', verbose=0)
{
	seq.DNAbin.matrix	<- as.character(seq.DNAbin.matrix)		
	seq.DNAbin.matrix	<- apply(seq.DNAbin.matrix, 2, function(col) 		gsub(code.from,code.to,col,fixed=1)			)	
	as.DNAbin( seq.DNAbin.matrix )
}

#' @export
hivc.seq.rmgaps<- function(seq.DNAbin.matrix, rm.only.col.gaps=1, verbose=0)
{
	seq.DNAbin.matrix		<- as.character(seq.DNAbin.matrix)		
	if(!rm.only.col.gaps)
	{	
		if(is.matrix(seq.DNAbin.matrix))
		{
			tmp					<- lapply(seq_len(nrow(seq.DNAbin.matrix)), function(i){	seq.DNAbin.matrix[i, seq.DNAbin.matrix[i,]!="-" & seq.DNAbin.matrix[i,]!="?"]	})
			names(tmp)			<- rownames(seq.DNAbin.matrix)
		}
		else
		{
			tmp					<- lapply(seq_along(seq.DNAbin.matrix), function(i){	seq.DNAbin.matrix[[i]][ seq.DNAbin.matrix[[i]]!="-" & seq.DNAbin.matrix[[i]]!="?"]	})
			names(tmp)			<- names(seq.DNAbin.matrix)
		}		
		seq.DNAbin.matrix	<- tmp
	}
	else
	{		
		nogap				<- which( !apply(seq.DNAbin.matrix,2,function(x) all(x=="-" || x=="?")) )
		if(verbose)	cat(paste("\nremove gaps, n=",ncol(seq.DNAbin.matrix)-length(nogap)))
		seq.DNAbin.matrix	<- seq.DNAbin.matrix[,nogap]	
	}
	as.DNAbin( seq.DNAbin.matrix )
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

#' Compute a set of phylogenetic clusters according to several criteria; taken from http://dx.doi.org/10.6084/m9.figshare.97225
#' Perform DFS. Label a subtree a cluster if both
#'   -its median pairwise patristic ditance (MPPD) is below a percentile of the whole-tree p-distance and
#'   -it is not a leaf. 
#' @export
#' @param ph					phylogenetic tree with branch lengths
#' @param thresh.brl 			MPPD threshold
#' @param dist.brl 				precomputed MPPD values for each node
#' @param thresh.nodesupport	boostrap reliability threshold
#' @param nodesupport			precomputed bootstrap values
#' @param retval				evaluate MPPD by tip nodes only or for all nodes
#' @return vector indicating, for each internal node, which cluster it belongs to
hivc.clu.clusterbythresh<- function(ph,thresh.brl=NULL,dist.brl=NULL,thresh.nodesupport=NULL,nodesupport=NULL,retval="all")
{
	require(ape)
	require(igraph)
	require(geiger)	
	if(all( is.null(c(thresh.brl, thresh.nodesupport))) ) 		stop("all threshold criteria NULL - provide at least one")
	if(!is.null(thresh.nodesupport) && is.null(nodesupport))	stop("node support threshold set but no node support values provided")
	if(!is.null(thresh.brl) && is.null(dist.brl))				stop("branch length threshold set but no branch length distances provided")
	if(!is.null(nodesupport) && any(nodesupport>1+EPS))			warning("Found nodesupport values above 1")
	if(is.null(thresh.brl))
	{		
		dist.brl			<- rep(1,Nnode(ph))
		thresh.brl			<- 1
	}	
	if(is.null(thresh.nodesupport))
	{
		nodesupport			<- rep(1,Nnode(ph))
		thresh.nodesupport	<- 1		
	}
#print(nodesupport[1:10]); print(thresh.nodesupport)
	## set up clustering
	ntips		<- Ntip(ph)	
	clu.i 		<- 0 ## cluster number
	clu.mem 	<- rep(NA,ntips+ph$Nnode) ## cluster member assignment
	clu.idx		<- rep(NA,ntips+ph$Nnode) ## cluster index assignment
	igraph.ph	<- graph.edgelist(ph$edge) ## ph in igraph form
	dfs 		<- graph.dfs(igraph.ph,root=ntips+1,neimode='out',order=TRUE,dist=TRUE)
#print(c(length(dfs$order),Ntip(ph),Nnode(ph)))
	## travese the ph in depth first order
	for(i in 1:length(dfs$order))
	{
		node <- dfs$order[i]
#print( c(node,node-ntips, is.na(clu.mem[node]), thresh.brl, dist.brl[node-ntips], dist.brl[node-ntips]<=thresh.brl, thresh.nodesupport, nodesupport[node-ntips], nodesupport[node-ntips]<thresh.nodesupport) )
		if(	node > ntips	&&											## skip leaves
			is.na(clu.mem[node]) &&										## only consider unassigned nodes
			nodesupport[node-ntips]>=thresh.nodesupport-EPS &&
			dist.brl[node-ntips]<=thresh.brl+EPS
			)	 
		{				
#print(nodesupport[node-ntips]);print(c(node,i))			
			clu.i 			<- clu.i+1
			subtree 		<- graph.dfs(igraph.ph,node,neimode='out',unreachable=FALSE)$order
			subtree 		<- subtree[! is.nan(subtree)]
			clu.mem[subtree]<- clu.i
			clu.idx[node]	<- clu.i
		}
	}
	ans <- list( clu.mem=clu.mem, clu.idx=which(!is.na(clu.idx)), size.all=table(clu.mem), size.tips=table(clu.mem[1:ntips]), ntips=ntips, thresh.brl=thresh.brl, thresh.nodesupport=thresh.nodesupport)
	if(retval=="tips")
		ans$clu.mem 		<- ans$clu.mem[1:ntips]
	return(ans)
}
######################################################################################
hivc.clu.exp.typeIerror.randomclu<- function(ph, dist.brl, nodesupport, ph.unlinked.seroneg, ph.unlinked.dead, thresh.brl, thresh.bs)
{
	ph.tips.n				<- Ntip(ph)
	#seroneg.n				<- nrow(ph.unlinked.seroneg)
	#ph.unlinked.seroneg.v	<- ph.unlinked.seroneg[,PhNode]
	clustering				<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=nodesupport, retval="all")
	clusters				<- unique(clustering[["clu.mem"]][!is.na(clustering[["clu.mem"]])])
	clu.fp.exp	<- sapply(clusters,function(clu)		
			{		
				tmp							<- which(clustering[["clu.mem"]]==clu)
				clu.leaves					<- tmp[ tmp<=ph.tips.n ]
				clu.leaves.seroneg			<- sapply(clu.leaves, function(x)
						{
							as.numeric( subset(ph.unlinked.seroneg, x==PhNode, c(PhNodeUnlinked,NegT) ) )	 	#O(logn) search; returns NA if subset dt is empty																 
						})				
#print(clu); print(tmp); print(clu.leaves); print(clu.leaves.seroneg)
				#if all leaves not seroneg, return 0								
				if(all(is.na(clu.leaves.seroneg[1,])))	return(0)
				#if at least one seroneg, compute prob that 1-P(X=0) where X is number of seroconverters dead at SC of latest seroconverter				
				clu.leaves.seroneg			<- rbind(clu.leaves,clu.leaves.seroneg)[,!is.na(clu.leaves.seroneg[2,]),drop=0]								
				clu.leaves.latestseroneg.ix	<- which.max( clu.leaves.seroneg[3,] )
				tmp							<- sapply(clu.leaves.seroneg[2,],function(x) length(ph.unlinked.dead[[x]]) )
				latest.seroneg.unlinked		<- ph.unlinked.dead[[ clu.leaves.seroneg[2,clu.leaves.latestseroneg.ix] ]]				
				if( max(tmp)!=length(latest.seroneg.unlinked) )	
					stop("expect latest seroneg to have most unlinked HIV+ dead associated")	#which.max does not work because .ix is not necessarily the first																
				#m is "number false pos" is length(latest.seroneg.unlinked)
				#n is "number false pos not known" is ph.tips.n-m
				#k is "number of draws" is length(clu.leaves)-1 
				ans<- 1-dhyper(0, length(latest.seroneg.unlinked), ph.tips.n-length(latest.seroneg.unlinked), length(clu.leaves)-1)
				#c(ans, length(clu.leaves)-1, length(latest.seroneg.unlinked), ph.tips.n)
				ans
			})
	list(fp.data=clu.fp.exp, fp.n=sum(clu.fp.exp), fp.rate= sum(clu.fp.exp)/length(clusters), clustering=clustering)
}
######################################################################################
hivc.clu.truepos<- function(clustering, ph.linked, ph.tips.n, verbose=0)
{
	clusters	<- unique(clustering[["clu.mem"]][!is.na(clustering[["clu.mem"]])])
	tmp			<- which(!is.na(clustering[["clu.mem"]][seq_len(ph.tips.n)]))
	clu.linked	<- merge( data.table(Node=tmp), ph.linked, by="Node" )	
	clu.linked	<- clu.linked[,list(nTP.clu= length(FASTASampleCode)), by="Patient"]
	clu.linked	<- merge(ph.linked, clu.linked, all.x=1, by="Patient")
	set(clu.linked, which(is.na(clu.linked[,nTP.clu])) ,"nTP.clu", 1L )
	setkey(clu.linked, Node)
	
	if(!length(clusters))
	{
		clu.tp		<- cbind( subset(clu.linked,select=c(Patient,nTP,nTP.clu)), TP.clu=0)
		setnames(clu.tp, c("nTP","nTP.clu"),c("TP.n","TP.n.clu"))				
	}
	else
	{				
		clu.tp	<- lapply(c(-1,clusters), function(clu)		
			{		
				#clu<- 48
				if(clu==-1)					#bug in data.table; if first data.table empty, rbind fails		http://lists.r-forge.r-project.org/pipermail/datatable-help/2012-November/001372.html
					return( as.data.table(data.frame(Patient="XXX",TP.clu=0,TP.n=0,TP.n.clu=0,clu=0,clu.tips=0,stringsAsFactors=F)) )
				tmp							<- which(clustering[["clu.mem"]]==clu)
				clu.tips					<- tmp[ tmp<=ph.tips.n ]				
				clu.tips.linked				<- na.omit(merge( data.table(Node=clu.tips, key="Node"), clu.linked, all.x=1))	#extract clu.tips that should be linked
				if(!nrow(clu.tips.linked))	
					ans						<- as.data.table(na.omit(data.frame(Patient=NA,TP.clu=NA,TP.n=NA, TP.n.clu=NA, clu=NA, clu.tips=NA)))
				else	
				{
					ans						<- clu.tips.linked[, list(TP.clu=length(FASTASampleCode), TP.n=nTP[1], TP.n.clu=nTP.clu[1], clu=clu, clu.tips=length(clu.tips) ), by=Patient]
					set(ans, NULL, "Patient", as.character(ans[,Patient]))
				}
				ans
			})
		#rbind and remove bugfix first row
		clu.tp		<- rbindlist(clu.tp)[-1,]	
		
	}	
	# account for different seq numbers per patient with a combinatorial argument 
	# -- only important when the number of correctly clustered within-patient sequences is evaluated
	clu.tp[, pairs.TP.n:= choose(TP.n,2)]
	clu.tp[, pairs.TP.n.clu:= choose(TP.n.clu,2)]
	clu.tp[, pairs.TP.clu:= choose(TP.clu,2)]
	# patients in wrong cluster
	clu.diffclu								<- subset( clu.tp[,list(nclu=length(na.omit(unique(clu)))),by="Patient"], nclu>1, Patient )
	clu.diffclu								<- merge( clu.diffclu, clu.tp, all.x=1, by="Patient" )

	# evaluate among all sequences in tree
	patients.allseqinsameclu				<- nrow( subset(clu.tp, TP.clu==TP.n) )							#if TP.clu==TP.n then Patient can occur only once in subset
	patients.total							<- length(unique(clu.linked[,Patient]))
	tmp										<- clu.tp[, list(TP.clu= max(TP.clu), pairs.TP.clu= max(pairs.TP.clu), TP.n=TP.n[1], pairs.TP.n=pairs.TP.n[1]),by=Patient]	#if TP.clu!=TP.n then Patient occurs multiple times
	patients.freqinsameclu					<- sum(tmp[, TP.clu]) / sum(tmp[, TP.n])
	# preferred way to quantify  the number of correctly clustered within-patient sequences
	pairs.freqinsameclu						<- sum(tmp[, pairs.TP.clu]) / sum(tmp[, pairs.TP.n])
	
	# evaluate among all sequences that cluster
	patients.clustering.allseqindiffclu		<- length(unique(clu.wc[,Patient]))
	patients.clustering.allseqinsameclu		<- nrow( subset(clu.tp, TP.clu==TP.n.clu & TP.n.clu>1) )
	patients.clustering.total				<- length(unique(subset(clu.linked,nTP.clu>1)[,Patient]))
	tmp										<- subset(clu.tp, TP.n.clu>1)[, list(TP.clu= max(TP.clu), pairs.TP.clu=max(pairs.TP.clu), TP.n.clu=TP.n.clu[1], pairs.TP.n.clu=pairs.TP.n.clu[1]),by=Patient]
	patients.clustering.freqinsameclu		<- sum(tmp[, TP.clu]) / sum(tmp[,TP.n.clu])
	# preferred way to quantify  the number of correctly clustered within-patient sequences	
	pairs.clustering.freqinsameclu			<- sum(tmp[, pairs.TP.clu]) / sum(tmp[,pairs.TP.n.clu])
	
	list(	clu.n				= length(clusters),	
			clu.onlytp			= subset(clu.tp, TP.clu==TP.n),
			clu.missedtp		= subset(clu.tp, TP.clu!=TP.n),
			clu.diffclu 		= clu.diffclu,
			tp.by.all			= patients.allseqinsameclu / patients.total,
			tp.by.sum			= pairs.freqinsameclu,
			tpclu.by.all		= patients.clustering.allseqinsameclu / patients.clustering.total,			
			tpclu.by.diffclu	= (patients.clustering.total - patients.clustering.allseqindiffclu) / patients.clustering.total,
			tpclu.by.sum		= pairs.clustering.freqinsameclu			
			)	
}
######################################################################################
hivc.clu.trueneg<- function(clustering, ph.unlinked.info, ph.unlinked, ph.tips.n, verbose=0)
{
	clusters			<- unique(clustering[["clu.mem"]][!is.na(clustering[["clu.mem"]])])
	clu.fp				<- sapply(clusters, function(clu)		
							{		
								tmp							<- which(clustering[["clu.mem"]]==clu)
								clu.tips					<- tmp[ tmp<=ph.tips.n ]
								clu.tips					<- merge( data.table(Node=clu.tips, key="Node"), subset(ph.unlinked.info,select=c(Node,NodeIdx,NegT)) )						
								if(all( is.na(clu.tips[,NodeIdx]) ))	return(rep(NA,2))
								if(!nrow(clu.tips))	return(rep(NA,2))				
								tmp							<- clu.tips[,NegT]
								tmp							<- ifelse(all(is.na(tmp)),1,which.max(tmp))	
								clu.tips.unlinked			<- ph.unlinked[[ clu.tips[tmp,NodeIdx] ]]			
								#print(clu.tips)
								c(nrow( merge(clu.tips, clu.tips.unlinked, by="Node") ), nrow(clu.tips.unlinked))																					
							})
	colnames(clu.fp)	<- clusters		
	clu.fp				<-	clu.fp[, apply(!is.na(clu.fp),2,all),drop=0]
	list( 	clu.n		= length(clusters), 
			clu.fp		= clu.fp, 
			fpn.by.all	= length(which(clu.fp[1,]>0)), 
			fpn.by.sum	= sum(clu.fp[1,]), 
			fp.by.all	= length(which(clu.fp[1,]>0))/length(clusters)
			)
}
######################################################################################
hivc.clu.clusterbytruepos<- function(ph, dist.brl, nodesupport, ph.linked, thresh.brl=NULL, thresh.bs=NULL, level= 0.95, tol= 0.005, mxit= 10, thresh.bs.lower= min(nodesupport), thresh.brl.upper=max(dist.brl), method.tp="tp.rate.tpclu",verbose=0)
{		
	if(is.null(thresh.brl) && is.null(thresh.bs))	stop("either thresh.bs or thresh.brl required")
	if(!is.null(thresh.brl) && !is.null(thresh.bs))	stop("nothing to do")
	
	error								<- 1
	ph.tips.n							<- Ntip(ph)
	nit									<- 0	
	if(is.null(thresh.bs))
	{
		srch.find.thresh.bs				<- 1
		srch.upper						<- 1
		srch.lower						<- thresh.bs	<- thresh.bs.lower		
	}
	else
	{
		srch.find.thresh.bs				<- 0
		srch.lower						<- thresh.brl 	<- 0
		srch.upper						<- thresh.brl.upper		
	}
	if(verbose)
		cat(paste("\ncalibrate bs?",srch.find.thresh.bs,"\nsrch.lower",srch.lower,"srch.upper",srch.upper))
	#binary search algorithm
	while(	abs(error)>tol && srch.lower<srch.upper && nit<mxit)
	{	
		if(nit)
		{
			if(srch.find.thresh.bs)	
					thresh.bs			<- ( srch.upper + srch.lower ) / 2
			else
					thresh.brl			<- ( srch.upper + srch.lower ) / 2
		}
		clustering					<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=nodesupport, retval="all")
		tmp							<- hivc.clu.truepos(clustering, ph.linked, ph.tips.n)
print(tmp)		
		if(verbose)	
			cat(paste("\nit=",nit,", #clu=",tmp[["clu.n"]],", BS=",thresh.bs,", BRL=",thresh.brl,", avg of %TP in cluster=",round(tmp[[method.tp]],digits=4),sep=''))			
		error						<- tmp[[method.tp]]	- level	
#print(c(srch.upper,srch.lower, thresh.bs, thresh.brl))
		if(srch.find.thresh.bs)	
		{
			if(error<0)
				srch.upper			<- thresh.bs
			else
				srch.lower			<- thresh.bs							
		}
		else
		{
			if(error<0)
				srch.lower			<- thresh.brl				
			else
				srch.upper			<- thresh.brl				
		}		
		nit							<- nit+1
	}
	ans<- list(thresh.bs=thresh.bs, thresh.brl=thresh.brl, clu=clustering, clu.fp.n=error+level, srch.nit= nit, srch.error=error, srch.tol=tol)
	ans		
}
######################################################################################
hivc.clu.clusterbytrueneg<- function(ph, dist.brl, nodesupport, ph.unlinked.info, ph.unlinked, thresh.brl=NULL, thresh.bs=NULL, level= 0.01, tol= 0.005, mxit= 20, thresh.bs.lower= min(nodesupport), thresh.brl.upper=max(dist.brl), verbose=0)
{		
	if(is.null(thresh.brl) && is.null(thresh.bs))	stop("either thresh.bs or thresh.brl required")
	if(!is.null(thresh.brl) && !is.null(thresh.bs))	stop("nothing to do")
	
	error								<- 1
	ph.tips.n							<- Ntip(ph)
	nit									<- 0	
	if(is.null(thresh.bs))
	{
		srch.find.thresh.bs				<- 1
		srch.upper						<- 1
		srch.lower						<- thresh.bs	<- thresh.bs.lower		
	}
	else
	{
		srch.find.thresh.bs				<- 0
		srch.lower						<- 0
		srch.upper						<- thresh.brl	<- thresh.brl.upper		
	}
	if(verbose)
		cat(paste("\ncalibrate bs?",srch.find.thresh.bs,"\nsrch.lower",srch.lower,"srch.upper",srch.upper,"thresh.bs",thresh.bs,"thresh.brl",thresh.brl))
	#binary search algorithm
	while(abs(error)>tol && srch.lower<srch.upper && nit<mxit)
	{	
		if(nit)
		{
			if(srch.find.thresh.bs)	
				thresh.bs			<- ( srch.upper + srch.lower ) / 2
			else
				thresh.brl			<- ( srch.upper + srch.lower ) / 2				
		}
			
		clustering					<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=nodesupport, retval="all")
		tmp							<- hivc.clu.trueneg(clustering, ph.unlinked.info, ph.unlinked, ph.tips.n, verbose=0)
		clustering.fp				<- tmp$clu.fp
		if(verbose)	
			cat(paste("\nit=",nit,", #clu=",tmp[["clu.n"]],", BS=",thresh.bs,", BRL=",thresh.brl,", #unlinked in any cluster=",tmp[["fp.n.sum"]],", #FP=",tmp[["fp.n"]],", %FP=",round(tmp[["fp.rate"]],d=5),sep=''))			
		error						<- tmp[["fp.rate"]]	- level	
#print(c(srch.upper,srch.lower, thresh.bs, thresh.brl))
		if(srch.find.thresh.bs)	
		{
			if(error<0)
				srch.upper			<- thresh.bs
			else
				srch.lower			<- thresh.bs							
		}
		else
		{
			if(error<0)
				srch.lower			<- thresh.brl				
			else
				srch.upper			<- thresh.brl				
		}
		nit							<- nit+1
	}
	ans<- list(thresh.bs=thresh.bs, thresh.brl=thresh.brl, clu=clustering, clu.fp=clustering.fp, srch.nit= nit, srch.error=error, srch.tol=tol)
	ans		
}
######################################################################################
hivc.clu.polyphyletic.clusters<- function(ph, clustering, cluphy.df, verbose=1, plot.file=NA, pdf.scaley=25, cex.nodelabel=0.2, cex.tiplabel=0.2, adj.tiplabel= c(-0.15,0.5))
{
	#get node corresponding to index of selected clusters		
	cluphy.cluidx				<- clustering[["clu.idx"]][ unique( cluphy.df[,cluster] ) ]		
	if(any(is.na(cluphy.cluidx)))	stop("unexcpected NA in cluphy.cluidx")
	#get selected clusters into single phylogeny
	cluphy.subtrees				<- lapply(cluphy.cluidx, function(x)		extract.clade(ph, x, root.edge= 1, interactive = FALSE) 		)
	names(cluphy.subtrees)		<- unique( cluphy.df[,cluster] )
	cluphy						<- eval(parse(text=paste('cluphy.subtrees[[',seq_along(cluphy.subtrees),']]', sep='',collapse='+')))	
	#plot selected clusters
	if(!is.na(plot.file))
	{		
		cluphy.tiplabels		<- hivc.clu.get.tiplabels( cluphy, cluphy.df )
		if(verbose) cat(paste("\nwrite tree to file",plot.file))
		hivc.clu.plot(cluphy, cluphy.df[,cluster], file=plot.file, pdf.scaley=pdf.scaley, pdf.off=0, cex.nodelabel=cex.nodelabel )											
		hivc.clu.plot.tiplabels( seq_len(Ntip(cluphy)), cluphy.tiplabels$text, cluphy.tiplabels$col, cex=cex.tiplabel, adj=adj.tiplabel, add.xinch=0, add.yinch=0 )
		dev.off()
	}	
	list(cluphy=cluphy, cluphy.subtrees=cluphy.subtrees)
}	
######################################################################################
hivc.clu.brl.bwpat<- function(cluphy.subtrees, cluphy.df, rtn.val.for.na.patient=FALSE)
{
	library(adephylo)
	#get branch length matrices for each subtree 
	cluphy.disttips				<- lapply(cluphy.subtrees, distTips)
	#select branch lengths between seqs of different patients 
	setkey(cluphy.df, FASTASampleCode)
	cluphy.brl.bwpat			<- lapply(cluphy.disttips, function(x)
			{
				x.patient		<- subset( cluphy.df[attr(x, "Labels"),], select=Patient )
				tmp				<- CJ(p1= seq_len(nrow(x.patient)), p2= seq_len(nrow(x.patient)))
				tmp[,lower.tri:=tmp[,p1<p2]]
				x.bwpat.matidx	<- tmp[, lower.tri & x.patient[p1,Patient]!=x.patient[p2,Patient]]		#if Patient ID not know, matidx is NA
				if(rtn.val.for.na.patient)					
					x.bwpat.matidx[is.na(x.bwpat.matidx)]<- TRUE										
				as.matrix(x)[x.bwpat.matidx]															#if Patient ID not know, returns NA
			}) 
	names(cluphy.brl.bwpat)		<- names(cluphy.subtrees)		
	cluphy.brl.bwpat
}
######################################################################################
hivc.clu.getplot.incountry<- function(ph, clustering, df.cluinfo, verbose=1, plot.file= NA, char.frgn='CountryInfection=="FRGN"', char.frgntn='CountryInfection=="FRGNTN"')
{
	clut						<- parse(text=paste("!is.na(CountryInfection) & (",char.frgn," | ",char.frgntn,") ", sep=''))		
	clut						<- table( df.cluinfo[,cluster, eval(clut) ] )
	#select clusters containing only in-country seq
	clut.nofrgn					<- apply(clut, 2, function(x)	x[1]>0 & x[2]==0		)
	clut.nofrgn					<- as.numeric( colnames(clut[,clut.nofrgn]) )
	if(verbose) cat(paste("\nnumber of clusters with only in-country seq is n=", length(clut.nofrgn)))
	cluphy.nofrgn				<- merge(data.table(cluster=clut.nofrgn), df.cluinfo, by="cluster")
	if(verbose) cat(paste("\nnumber of seq in clusters with only in-country seq is n=", nrow(cluphy.nofrgn)))
	#select clusters containing at least one in-country seq
	clut.mixed					<- apply(clut, 2, function(x)	all(x>0)		)		
	clut.mixed					<- as.numeric( colnames(clut[,clut.mixed]) )
	if(verbose) cat(paste("\nnumber of clusters with at least one in-country seq is n=", length(clut.mixed)))
	cluphy.mixed				<- merge(data.table(cluster=clut.mixed), df.cluinfo, by="cluster")
	if(verbose) cat(paste("\nnumber of seq in clusters with at least one in-country seq is n=", nrow(cluphy.mixed)))
	#
	# split mixed clusters
	#
	#exclude pairs -- split would result in singleton
	tmp							<- cluphy.mixed[, list(size=length(unique(FASTASampleCode))), by="cluster"]
	cluphy.mixed.rm				<- subset(tmp,size<=2, cluster)[,cluster]		
	cluphy.mixed				<- merge( subset(tmp,size>2, cluster), df.cluinfo, all.x=1, by="cluster" )
	if(verbose) cat(paste("\nnumber of splitting clusters is n=", length(unique(cluphy.mixed[,cluster])) ))
	if(length(intersect( unique(cluphy.mixed[,cluster]),unique(cluphy.nofrgn[,cluster]) )))
		stop("unexpected overlap between 'cluphy.mixed' and 'cluphy.nofrgn'")		
	#extract selected clusters		
	cluphy.cluidx				<- clustering[["clu.idx"]][ unique( cluphy.mixed[,cluster] ) ]					
	cluphy.subtrees				<- lapply(cluphy.cluidx, function(x)		extract.clade(ph, x, root.edge= 1, interactive = FALSE) 		)
	names(cluphy.subtrees)		<- unique( cluphy.mixed[,cluster] )				
	#split selected clusters				
	#subset(cluphy.mixed, select=c(cluster,FASTASampleCode,Patient,CountryInfection))
	#splitexpr<- parse(text=paste( char.frgn, char.frgntn, sep=" | " ))
	#cluphy.df<- df.cluinfo
	tmp							<- hivc.clu.splitcluster(ph, clustering, cluphy.subtrees, df.cluinfo, parse(text=paste( char.frgn, char.frgntn, sep=" | " )) )
	#z<- subset(tmp$cluphy.df[cluphy.mixed[,FASTASampleCode],], !is.na(cluster),select=c(cluster,FASTASampleCode,Patient,CountryInfection))
	#length(unique(z[,cluster]))	
	cluphy.subtrees				<- tmp$cluphy.subtrees
	if(verbose) cat(paste("\nnumber of splitted clusters is n=", length(cluphy.subtrees)))
	clustering					<- tmp$clustering 
	df.cluinfo					<- tmp$cluphy.df
	#
	# remove previously excluded pairs
	#
	set(df.cluinfo,which(df.cluinfo[,cluster]==cluphy.mixed.rm),"cluster",NA)
	df.cluinfo																<- subset(df.cluinfo, !is.na(cluster) )
	clustering[["clu.mem"]][ clustering[["clu.mem"]]%in%cluphy.mixed.rm ]	<- NA
	clustering[["clu.idx"]][ cluphy.mixed.rm ]								<- NA
	clustering[["size.tips"]]												<- table( clustering[["clu.mem"]][1:Ntip(ph)] )
	# update clustering after splitting
	cluphy.mixed				<- subset(df.cluinfo[cluphy.mixed[,FASTASampleCode],], !is.na(cluster))				
	# combine 'cluphy.mixed', 'cluphy.nofrgn'
	setcolorder(cluphy.mixed, colnames(cluphy.nofrgn))
	cluphy.df					<- rbind(cluphy.mixed, cluphy.nofrgn)	
	if(verbose) cat(paste("\nnumber of in-country clusters is n=", length(unique(cluphy.df[,cluster]))))
	if(verbose) cat(paste("\nnumber of seq in in-country clusters is n=", nrow(cluphy.df)))
	if(verbose) cat(paste("\nCHECK: number of in-country clusters in df.cluinfo is n=", length(unique(df.cluinfo[,cluster]))))
	#
	# build polyphyletic tree from clusters
	#		
	tmp							<- hivc.clu.polyphyletic.clusters(ph, clustering, cluphy.df, plot.file=plot.file, verbose=verbose )
	cluphy						<- tmp$cluphy
	cluphy.subtrees				<- tmp$cluphy.subtrees
	
	list(clustering=clustering, df.cluinfo=df.cluinfo, cluphy=cluphy, cluphy.subtrees=cluphy.subtrees, cluphy.df=cluphy.df)
}
######################################################################################
hivc.clu.getplot.msmexposuregroup<- function(ph, clustering, df.cluinfo, verbose=1, plot.file= NA, levels.msm=c("BI","MSM","IDU","NA"), levels.het=c("BI","HET","IDU","NA"), levels.mixed=c("BI","MSM","HET","IDU","NA"),levels.oth="OTH", split.clusters=0)
{	
	clut						<- table( df.cluinfo[,cluster,Trm], useNA= "always" )
	rownames(clut)[nrow(clut)]	<- "NA"
	clut						<- clut[,-ncol(clut)]
	clut						<- rbind(clut, apply(clut, 2, sum) )
	rownames(clut)[nrow(clut)]	<- "sum"
	#
	# clusters with MSM only
	#
	clut.onlymsm				<- apply(clut, 2, function(x)		sum( x[levels.msm] )==x["sum"]		)
	clut.onlymsm				<- as.numeric(colnames(clut[,clut.onlymsm]))
	cluphy.onlymsm				<- merge(data.table(cluster=clut.onlymsm), df.cluinfo, by="cluster")	
	if(verbose)	cat(paste("\nnumber of clusters with MSM only seq is n=",length(unique(cluphy.onlymsm[,cluster]))))		
	#
	# clusters with MSM/HET
	#
	clut.het					<- apply(clut, 2, function(x)		sum( x[levels.msm] )!=x["sum"] & sum( x[levels.het] )!=x["sum"] & sum( x[levels.mixed] )==x["sum"]		)
	clut.het					<- as.numeric(colnames(clut[,clut.het]))
	cluphy.het					<- merge(data.table(cluster=clut.het), df.cluinfo, by="cluster")
	if(verbose)	cat(paste("\nnumber of clusters with MSM/HET seq is n=",length(unique(cluphy.het[,cluster]))))
	#	
	# clusters with MSM/HET-M
	#	
	cluphy.hetM 				<- cluphy.het[,list( HetM.only= all(Sex=="M") ),by=cluster]
	if(verbose)	cat(paste("\nnumber of clusters with MSM/HET-M seq is n=",nrow(subset(cluphy.hetM,HetM.only))))
	cluphy.hetM					<- merge(subset(cluphy.hetM,HetM.only,select=cluster), df.cluinfo, all.x=1, by="cluster")				
	#
	# clusters with MSM/HET-F
	#
	cluphy.hetF 				<- cluphy.het[,list( with.HETF= any(Sex=="F") ),by=cluster]
	if(verbose)	cat(paste("\nnumber of clusters with MSM/HET-F seq is n=",nrow(subset(cluphy.hetF,with.HETF))))
	cluphy.hetF					<- merge(subset(cluphy.hetF,with.HETF,select=cluster), df.cluinfo, all.x=1, by="cluster")
	if(length(intersect( unique(cluphy.hetF[,cluster]),unique(cluphy.hetM[,cluster]) )))
		stop("unexpected overlap between cluphy.hetF and cluphy.hetM")
	if(length(intersect( unique(cluphy.onlymsm[,cluster]),unique(cluphy.hetM[,cluster]) )))
		stop("unexpected overlap between cluphy.onlymsm and cluphy.hetM")		
	#
	# clusters with OTH
	#	
	clut.oth					<- apply(clut, 2, function(x)		x["MSM"]>0 & sum( x[levels.msm] )!=x["sum"] & x[levels.oth]>0 & sum( x[levels.oth] )!=x["sum"] & sum( x[levels.het] )!=x["sum"] & sum( x[unique(c(levels.msm,levels.het,levels.oth))] )==x["sum"]		)
	clut.oth					<- as.numeric(colnames(clut[,clut.oth]))
	cluphy.oth					<- merge(data.table(cluster=clut.oth), df.cluinfo, by="cluster")
	if(verbose)	cat(paste("\nnumber of clusters with OTH seq is n=",length(unique(cluphy.oth[,cluster]))))
	# combine cluphy.oth and cluphy.hetF for splitting
	if(length(intersect( unique(cluphy.hetF[,cluster]),unique(cluphy.oth[,cluster]) )))
		stop("unexpected overlap between cluphy.hetF and cluphy.oth")
	cluphy.split				<- rbind(cluphy.hetF, cluphy.oth)
	# exclude pairs -- split would result in singleton
	tmp							<- cluphy.split[, list(size=length(unique(FASTASampleCode))), by="cluster"]
	cluphy.split.rm				<- subset(tmp,size<=2, cluster)[,cluster]		
	cluphy.split				<- merge( subset(tmp,size>2, cluster), df.cluinfo, all.x=1, by="cluster" )	
	# extract subtrees for MSM/HET-F clusters		
	cluphy.cluidx				<- clustering[["clu.idx"]][ unique( cluphy.split[,cluster] ) ]					
	cluphy.subtrees				<- lapply(cluphy.cluidx, function(x)		extract.clade(ph, x, root.edge= 1, interactive = FALSE) 		)
	names(cluphy.subtrees)		<- unique( cluphy.split[,cluster] )				
	if(split.clusters)	#split cluster at HETF or OTH sequence
	{
		# split HET-F, OTH clusters					
		#subset(cluphy.hetF, select=c(cluster,FASTASampleCode,Patient,Sex))
		#splitexpr<- parse(text='Sex=="F"')
		#cluphy.df<- df.cluinfo
		if(verbose)	cat(paste("\nnumber of splitting clusters is n=",length(unique(cluphy.split[,cluster]))))
		tmp							<- hivc.clu.splitcluster(ph, clustering, cluphy.subtrees, df.cluinfo, parse(text='Sex=="F" | Trm=="OTH"'))
		#z<- subset(tmp$cluphy.df[cluphy.hetF[,FASTASampleCode],], !is.na(cluster),select=c(cluster,FASTASampleCode,Patient,Sex,Trm))	
		cluphy.subtrees				<- tmp$cluphy.subtrees
		if(verbose)	cat(paste("\nnumber of splitted clusters is n=",length(cluphy.subtrees)))		
	}
	else				#drop HETF or OTH sequence from cluster
	{		
		if(verbose)	cat(paste("\nnumber of droptip clusters is n=",length(unique(cluphy.split[,cluster]))))
		tmp							<- hivc.clu.droptipincluster(ph, clustering, cluphy.subtrees, df.cluinfo, parse(text='Sex=="F" | Trm=="OTH"') )
		cluphy.subtrees				<- tmp$cluphy.subtrees
		if(verbose)	cat(paste("\nnumber of droptip clusters is n=",length(cluphy.subtrees)))						
	}
	clustering					<- tmp$clustering 
	df.cluinfo					<- tmp$cluphy.df
	# update clustering after splitting
	cluphy.split				<- subset(df.cluinfo[cluphy.split[,FASTASampleCode],], !is.na(cluster))
	#
	# remove previously excluded pairs
	#
	set(df.cluinfo,which(df.cluinfo[,cluster]%in%cluphy.split.rm),"cluster",NA)
	df.cluinfo																<- subset(df.cluinfo, !is.na(cluster) )
	clustering[["clu.mem"]][ clustering[["clu.mem"]]%in%cluphy.split.rm ]	<- NA
	clustering[["clu.idx"]][ cluphy.split.rm ]								<- NA
	clustering[["size.tips"]]												<- table( clustering[["clu.mem"]][1:Ntip(ph)] )
	#
	# combine 'cluphy.split', 'cluphy.hetM', 'cluphy.onlymsm'
	#
	setcolorder(cluphy.split, colnames(cluphy.hetM))
	cluphy.df					<- rbind( rbind(cluphy.split, cluphy.hetM), cluphy.onlymsm )
	# may have lost all MSM in cluster due to splitting 
	cluphy.df					<- subset( cluphy.df[,list(with.MSM= any(Trm=="MSM")),by="cluster"], with.MSM, cluster )
	if(verbose) cat(paste("\nnumber of MSM clusters is n=", nrow(cluphy.df)))
	cluphy.df					<- merge(cluphy.df, df.cluinfo, all.x=1, by="cluster")
	if(verbose) cat(paste("\nnumber of seq in MSM clusters is n=", nrow(cluphy.df)))	
	#
	# build polyphyletic tree from clusters
	#
	tmp							<- hivc.clu.polyphyletic.clusters(ph, clustering, cluphy.df, plot.file=plot.file )
	cluphy						<- tmp$cluphy
	cluphy.subtrees				<- tmp$cluphy.subtrees
	
	list(clustering=clustering, df.cluinfo=df.cluinfo, cluphy=cluphy, cluphy.subtrees=cluphy.subtrees, cluphy.df=cluphy.df)
}	
######################################################################################
hivc.clu.getplot.mixedexposuregroup<- function(ph, clustering, df.cluinfo, verbose=1, plot.file= NA, levels.msm=c("BI","MSM","IDU","NA"), levels.het=c("BI","HET","IDU","NA"), levels.oth=c("OTH","NA"))
{	
	clut						<- table( df.cluinfo[,cluster,Trm], useNA= "always" )
	rownames(clut)[nrow(clut)]	<- "NA"
	clut						<- clut[,-ncol(clut)]
	clut						<- rbind(clut, apply(clut, 2, sum) )
	rownames(clut)[nrow(clut)]	<- "sum"
	clut.sametype				<- apply(clut, 2, function(x)		any( c( sum( x[levels.msm] ),sum( x[levels.het] ),sum( x[levels.oth] )	)==x["sum"])		)
	clut.onlymsm				<- apply(clut, 2, function(x)		sum( x[levels.msm] )==x["sum"]		)	
	if(verbose) print( table(clut.sametype) )
	if(verbose) print( table(clut.onlymsm) )
	#get selected clusters into single phylogeny
	cluphy.df					<- merge(data.table(cluster=which(!clut.sametype)), df.cluinfo, by="cluster")
	tmp							<- hivc.clu.polyphyletic.clusters(ph, clustering, cluphy.df, plot.file=plot.file )
	
	ans<- list(cluphy=tmp$cluphy, cluphy.df=cluphy.df )
	ans
}	
######################################################################################
hivc.clu.getplot.female2female<- function( ph, clustering, df.cluinfo, plot.file=NA )
{	
	#select clusters for HET transmission
	clut						<- table( df.cluinfo[,cluster,Trm], useNA= "always")
	rownames(clut)[nrow(clut)]	<- "NA"
	clut						<- clut[,-ncol(clut)]
	clut						<- rbind(clut, apply(clut, 2, sum) )
	rownames(clut)[nrow(clut)]	<- "sum"		
	onlyhet.clut				<- apply(clut, 2, function(x)		sum( x["HET"] )==x["sum"]		)		
	onlyhet.df					<- subset(df.cluinfo, cluster%in%which(onlyhet.clut) )
	#select clusters for gender: all F
	tmp							<- onlyhet.df[, { 
				tmp<- length(unique(Patient)); 
				list(clu.bwpat=tmp>1, clu.genderF=all(Sex=="F") ) 
			}, by="cluster"]
	onlyhet.select				<- subset( tmp, clu.bwpat & clu.genderF, select=cluster )		
	cluphy.df					<- merge( onlyhet.select, onlyhet.df, by="cluster" )		
	#select clusters for country infection: at most one out of NL origin		
	tmp							<- cluphy.df[, list(cluster=cluster[1], Patient=Patient[1], CountryInfection=CountryInfection[1]), by="Patient"]
	tmp							<- table(tmp[,cluster,CountryInfection], useNA= "always")
	tmp							<- tmp[, -ncol(tmp)]
	rownames(tmp)[nrow(tmp)]	<- "NA"
	tmp							<- rbind(tmp, apply(tmp, 2, sum))
	rownames(tmp)[nrow(tmp)]	<- "sum"
	tmp							<- as.numeric( colnames(tmp)[ tmp["NL",] + tmp["NA",] + 1 >= tmp["sum",] ] )
	cluphy.df					<- subset( cluphy.df, cluster%in%tmp )
	tmp							<- hivc.clu.polyphyletic.clusters(ph, clustering, cluphy.df, plot.file=plot.file, pdf.scaley=4, cex.nodelabel=0.4, cex.tiplabel=0.4)
		
	ans<- list(cluphy=tmp$cluphy, cluphy.df=cluphy.df, cluphy.subtrees=tmp$cluphy.subtrees )
	ans
}	
######################################################################################
hivc.clu.getplot.mixedfrgngroup<- function(ph, clustering, df.cluinfo, verbose=1, plot.file= NA, char.frgn='CountryInfection=="FRGN"', char.frgntn='CountryInfection=="FRGNTN"' )							
{	
	clut						<- parse(text=paste("!is.na(CountryInfection) & (",char.frgn," | ",char.frgntn,") ", sep=''))		
	clut						<- table( df.cluinfo[,cluster, eval(clut) ] )
	clut.mixed					<- apply(clut, 2, function(x)	all(x>0)		)
	clut.mixed					<- as.numeric( colnames(clut[,clut.mixed]) )
	if(verbose) cat(paste("\nnumber of mixed clusters (in-country / frgn) is n=", length(clut.mixed) ))
	cluphy.df					<- merge(data.table(cluster=clut.mixed), df.cluinfo, by="cluster")		
	cluphy			<- NULL
	
	if(!is.na(plot.file))
	{
		tmp			<- hivc.clu.polyphyletic.clusters(ph, clustering, cluphy.df, plot.file=plot.file, pdf.scaley=4, cex.nodelabel=0.4, cex.tiplabel=0.4)
		cluphy		<- tmp$cluphy
	}
	list(cluphy=cluphy, cluphy.df=cluphy.df )			
}	
######################################################################################
#prepare standard format of tip labels -- requires df.info to be sorted along the tips as they appear in a phylogeny
hivc.clu.get.tiplabels<- function(ph, df.info)
{
	#set colors CountryInfection
	tmp					<- rep("transparent",nrow(df.info))
	df.info[,CountryInfection.col:=tmp]
	set(df.info, which(df.info[,CountryInfection=="NL"]), "CountryInfection.col", "#EF9708")
	#set colors Patient
	tmp					<- unique( df.info[,Patient] )
	tmp2				<- diverge_hcl(length(tmp), h = c(246, 40), c = 96, l = c(85, 90))		
	tmp					<- data.table(Patient=tmp, Patient.col=tmp2, key="Patient")
	df.info				<- merge(  df.info, tmp, all.x=1, by="Patient" )
	#set colors Sex
	tmp					<- rep("transparent",nrow(df.info))
	df.info[,Sex.col:=tmp]						
	set(df.info, which(df.info[,Sex=="M"]), "Sex.col", "#EF9708")
	#set colors MSM or BI
	tmp					<- rep("transparent",nrow(df.info))
	df.info[,Trm.col:=tmp]						
	set(df.info, which(df.info[,Trm%in%c("MSM","BI")]), "Trm.col", "#EF9708")
	#set colors RegionHospital
	tmp					<- levels( df.info[,RegionHospital] )		
	tmp2				<- brewer.pal(length(tmp), "Dark2")
	tmp					<- data.table(RegionHospital=tmp, RegionHospital.col=tmp2, key="RegionHospital")
	df.info				<- merge(  df.info, tmp, all.x=1, by="RegionHospital" )		
	#set colors time		
	tmp					<- range( range(df.info[, NegT],na.rm=1), range(df.info[, AnyPos_T1],na.rm=1) )
	tmp					<- as.POSIXlt( seq.Date(tmp[1],tmp[2]+365,by="years") )$year
	tmp2				<- heat_hcl(length(tmp), h = c(0, -100), l = c(75, 40), c = c(40, 80), power = 1)
	yearcols			<- data.table(Year=tmp, Year.col=tmp2)
	#set colors PosSeqT
	tmp					<- data.table(PosSeqT= unique( df.info[,PosSeqT] ), key="PosSeqT" )
	tmp2				<- tmp[, as.POSIXlt(PosSeqT)$year]
	tmp[,"Year":=tmp2]
	tmp					<- subset( merge(tmp, yearcols, all.x=1, by="Year"), select=c(PosSeqT,Year.col))
	setnames(tmp,"Year.col","PosSeqT.col")
	df.info			<- merge(  df.info, tmp, all.x=1, by="PosSeqT" )
	#set colors NegT
	tmp					<- data.table(NegT= unique( df.info[,NegT] ), key="NegT" )
	tmp2				<- tmp[, as.POSIXlt(NegT)$year]
	tmp[,"Year":=tmp2]
	tmp					<- subset( merge(tmp, yearcols, all.x=1, by="Year"), select=c(NegT,Year.col))
	setnames(tmp,"Year.col","NegT.col")
	df.info			<- merge(  df.info, tmp, all.x=1, by="NegT" )
	#set colors AnyPos_T1
	tmp					<- data.table(AnyPos_T1= unique( df.info[,AnyPos_T1] ), key="AnyPos_T1" )
	tmp2				<- tmp[, as.POSIXlt(AnyPos_T1)$year]
	tmp[,"Year":=tmp2]
	tmp					<- subset( merge(tmp, yearcols, all.x=1, by="Year"), select=c(AnyPos_T1,Year.col))
	setnames(tmp,"Year.col","AnyPos_T1.col")
	df.info			<- merge(  df.info, tmp, all.x=1, by="AnyPos_T1" )									
	
	#convert time to string
	set(df.info,NULL,"AnyPos_T1",substr(as.character( df.info[,AnyPos_T1] ),1,7))
	set(df.info,NULL,"NegT",substr(as.character( df.info[,NegT] ),1,7))
	set(df.info,NULL,"PosSeqT",substr(as.character( df.info[,PosSeqT] ),1,7))
	
	#handle missing entries
	setkey(df.info, Patient)
	tmp					<- which(is.na(df.info[,Patient]))
	set(df.info, tmp, "Patient", '')
	set(df.info, tmp, "Patient.col", "transparent")
	tmp					<- which(is.na(df.info[,RegionHospital]))
	set(df.info, tmp, "RegionHospital", '-')
	set(df.info, tmp, "RegionHospital.col", "transparent")		
	tmp					<- which(is.na(df.info[,CountryInfection]))
	set(df.info, tmp, "CountryInfection", "--")
	set(df.info, tmp, "CountryInfection.col", "transparent")				
	tmp					<- which(is.na(df.info[,Trm]))
	set(df.info, tmp, "Trm", '')
	set(df.info, tmp, "Trm.col", "transparent")
	tmp					<- which(is.na(df.info[,Sex]))
	set(df.info, tmp, "Sex", '')
	set(df.info, tmp, "Sex.col", "transparent")
	tmp					<- which(is.na(df.info[,AnyPos_T1]))
	set(df.info, tmp, "AnyPos_T1", "-------")
	set(df.info, NULL, "AnyPos_T1", paste("HIV+:",df.info[,AnyPos_T1],sep=''))
	set(df.info, tmp, "AnyPos_T1.col", "transparent")
	tmp					<- which(is.na(df.info[,NegT]))
	set(df.info, tmp, "NegT", "-------")
	set(df.info, NULL, "NegT", paste("HIV-:",df.info[,NegT],sep=''))
	set(df.info, tmp, "NegT.col", "transparent")
	tmp					<- which(is.na(df.info[,PosSeqT]))
	set(df.info, tmp, "PosSeqT", "-------")
	set(df.info, NULL, "PosSeqT", paste("HIVS:",df.info[,PosSeqT],sep=''))
	set(df.info, tmp, "PosSeqT.col", "transparent")		
	#get df.info into order of tips
	setkey(df.info, FASTASampleCode)
	df.info				<- df.info[ph$tip.label,]
	#select text and col matrix 
	text				<- t( as.matrix( subset(df.info,select=c(CountryInfection, Trm, Sex, NegT, AnyPos_T1, PosSeqT, Patient, RegionHospital)) ) )
	colnames(text)		<- ph$tip.label
	col					<- t( as.matrix( subset(df.info,select=c(CountryInfection.col, Trm.col, Sex.col, NegT.col, AnyPos_T1.col, PosSeqT.col, Patient.col, RegionHospital.col)) ) )
	colnames(col)		<- ph$tip.label
	ans<- list(text=text, col=col)
	ans
}	
######################################################################################
hivc.clu.plot.tiplabels<- function (tip, text, col, adj = c(-0.05, 0.5), cex=1,add.xinch= 0.03, add.yinch= 0.02) 
{		
	lastPP 			<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	xx 				<- lastPP$xx[tip]
	yy 				<- lastPP$yy[tip]
	if(length(tip)==1)
	{
		#assume vector of text and col of same length
		if(!is.vector(text) || !is.vector(col))	stop("expect vector text and vector col")
		if(length(text)!=length(col))			stop("expect same length of text and col")
		wh				<- sapply(text, function(x){	 c(xinch(strwidth(x, units = "inches", cex = cex)), yinch(strheight(x, units = "inches", cex = cex))) 		})
		wh.total		<- c(sum(wh[1,]),max(wh[2,]))
		coord			<- matrix(NA,6,length(text),dimnames=list(c("xl","xr","yb","yt","xx","yy"),c()))
		coord["xl",]	<- xx - wh.total[1] * adj[1] - xinch(add.xinch)
		tmp				<- cumsum(wh[1,]) + xinch(add.xinch) / ncol(wh)
		coord["xl",]	<- coord["xl",] + c(0, tmp[-ncol(coord)])
		coord["xr",]	<- coord["xl",] + wh[1,] + xinch(add.xinch)/ ncol(wh)
		coord["yb",]	<- yy - wh.total[2] * adj[2] - yinch(add.yinch)
		coord["yt",]	<- coord["yb",] + wh.total[2] + yinch(add.yinch)
		coord["xx",]	<- coord["xl",] + c( diff( coord[1,] ), diff(coord[1:2,ncol(wh)])  ) / 2
		coord["yy",]	<- rep(yy, ncol(coord))
#print(wh); print(wh.total); print(coord)		
	}
	else
	{
		#assume matrix of text and col of same ncol
		if(!is.matrix(text) || !is.matrix(col))				stop("expect matrix text and vector col")
		if(ncol(text)!=ncol(col) || nrow(text)!=nrow(col))	stop("expect same dimensions of text and col")
		if(ncol(text)!=length(tip))							stop("expect cols of text and col to correspond to tips")
		coord	<- lapply( seq_along(xx),function(i)
				{
					wh				<- sapply(text[,i], function(x){	 c(xinch(strwidth(x, units = "inches", cex = cex)), yinch(strheight(x, units = "inches", cex = cex))) 		})
					wh.total		<- c(sum(wh[1,]),max(wh[2,]))
					coord			<- matrix(NA,6,nrow(text),dimnames=list(c("xl","xr","yb","yt","xx","yy"),c()))
					coord["xl",]	<- xx[i] - wh.total[1] * adj[1] - xinch(add.xinch)
					tmp				<- cumsum(wh[1,]) + xinch(add.xinch) / ncol(wh)
					coord["xl",]	<- coord["xl",] + c(0, tmp[-ncol(coord)])
					coord["xr",]	<- coord["xl",] + wh[1,] + xinch(add.xinch)/ ncol(wh)
					coord["yb",]	<- yy[i] - wh.total[2] * adj[2] - yinch(add.yinch)
					coord["yt",]	<- coord["yb",] + wh.total[2] + yinch(add.yinch)
					coord["xx",]	<- coord["xl",] + c( diff( coord[1,] ), diff(coord[1:2,ncol(wh)])  ) / 2
					coord["yy",]	<- rep(yy[i], ncol(coord))
					coord
				})
		coord	<- do.call("cbind",coord)
		text	<- as.vector(text)
		col		<- as.vector(col)
#print(coords); print(text); print(col)
	}
	rect(coord["xl",], coord["yb",], coord["xr",], coord["yt",], col = col, border=NA)
	text(coord["xx",], coord["yy",], text, cex=cex)
}
######################################################################################
hivc.clu.plot<- function(	ph, clu, edge.col.basic="black", show.node.label= T, file=NULL,  
							highlight.edge.of.tiplabel=NULL, highlight.edge.of.tiplabel.col="red", 
							highlight.cluster=NULL, highlight.cluster.col="red",							
							pdf.scaley=10, pdf.width= 7, pdf.height=pdf.scaley*7, pdf.off=1, 
							cex.nodelabel=0.5, cex.edge.incluster=1, cex.edge.outcluster= cex.edge.incluster/3)
{
	require(colorspace)	
	clu.edge							<- clu[ ph$edge[,1] ]
	clu.edge							<- clu.edge+1				#set col index for non-clustering edges to 1
	clu.edge[is.na(clu.edge)]			<- 1
	#cols.n								<- length(unique(clu.edge))-1	
	#cols								<- c("black",rainbow_hcl(cols.n, start = 30, end = 300))	
	clu.edge.col						<- rep(edge.col.basic, nrow(ph$edge))	#cols[clu.edge]
	clu.edge.col[clu.edge==1]			<- "grey50" 
	if(!is.null(highlight.cluster))
	{
		if(length(highlight.cluster.col)==1)
			highlight.cluster.col<- rep(highlight.cluster.col,length(highlight.cluster))
		for(i in seq_along(highlight.cluster))
			clu.edge.col[clu.edge%in%(highlight.cluster[[i]]+1)]<- 	highlight.cluster.col[i]	
	}	
	if(!is.null(highlight.edge.of.tiplabel))
	{
		highlight.edge<- lapply(highlight.edge.of.tiplabel, function(x)		which( ph$edge[,2]%in%which( substr(ph$tip.label, 1, nchar(x))==x ) )		)					
		if(length(highlight.edge.of.tiplabel.col)==1)
			highlight.edge.of.tiplabel.col	<- rep(highlight.edge, length(highlight.edge.of.tiplabel.col) )		
		for(i in seq_along(highlight.edge))
			clu.edge.col[highlight.edge[[i]]]		<- highlight.edge.of.tiplabel.col[i]		
	}	
	clu.edge.width						<- rep(cex.edge.outcluster, length(clu.edge))
	clu.edge.width[clu.edge!=1]			<- cex.edge.incluster
	clu.edge.lty						<- rep(1, length(clu.edge))
	clu.edge.lty[clu.edge!=1]			<- 1
	if(class(file)=="character")
		pdf(file,width=pdf.width,height=pdf.height)
	plot(ph, show.tip.label=F, show.node.label=show.node.label, cex=cex.nodelabel, edge.color=clu.edge.col, edge.width=clu.edge.width, edge.lty=clu.edge.lty)
	if(class(file)=="character" && pdf.off)
		dev.off()
}
######################################################################################
hivc.clu.plot.tptn<- function(stat.cmp.x, stat.cmp.y, plotfile, cols, xlab= "#FP (among all)", ylab= "%TP (among all)", xlim= range(stat.cmp.x), labels= as.numeric(colnames(stat.cmp.x)), verbose=1)
{		
	if(verbose) cat(paste("\nplot file",plotfile))
	pdf(file=plotfile,width=5,height=5)
	plot(1,1,type='n',bty='n',xlim=xlim,ylim=range(stat.cmp.y),xlab=xlab,ylab=ylab)
	dummy	<- sapply(seq_len(nrow(stat.cmp.x)),function(i)
			{
				points(stat.cmp.x[i,],stat.cmp.y[i,],col=cols[i],type='b')
				text(stat.cmp.x[i,],stat.cmp.y[i,],col=cols[i],labels=labels,adj=c(-0.8,0.5),cex=0.5)
			})
	legend("bottomright",border=NA,bty='n',fill=cols,legend=names(cols))
	dev.off()
}	
######################################################################################
hivc.clu.plot.withinpatientseq.not.samecluster<- function(missed.ph, clustering, clusters.tp, missed.df, plotfile, verbose=1)
{
	clusters.tp.missed									<- merge( data.table(Patient=unique( clusters.tp$clu.missedtp[,Patient] ) ), subset(missed.df,select=c(Patient, FASTASampleCode, cluster, Node)), by="Patient" )				
	# highlight sequences of patients that are not clustering in blue
	node.missed											<- subset(clusters.tp.missed,is.na(cluster),Node)[,Node]
	missed.ph$tip.label[ node.missed ]					<- paste("TPm",missed.ph$tip.label[ node.missed ],sep='-') 
	tmp													<- which( missed.df[,Node]%in%node.missed )
	set(missed.df, tmp, "FASTASampleCode", paste("TPm",missed.df[tmp,FASTASampleCode],sep='-') )		
	# highlight sequences of patients that are not in the same cluster in red
	node.wrong											<- data.table(Patient=unique(clusters.tp$clu.diffclu[,Patient]))
	node.wrong											<- merge(node.wrong , subset(missed.df,select=c(Patient, FASTASampleCode, cluster, Node)), by="Patient" )
	node.wrong											<- subset(node.wrong, !is.na(cluster), Node)[,Node]
	missed.ph$tip.label[ node.wrong ]					<- paste("TPw",missed.ph$tip.label[ node.wrong ],sep='-') 
	tmp													<- which( missed.df[,Node]%in%node.wrong )
	set(missed.df, tmp, "FASTASampleCode", paste("TPw",missed.df[tmp,FASTASampleCode],sep='-') )
	# highlight clusters with missed within patient seq 
	missed.clumem										<- clustering[["clu.mem"]]
	tmp													<- na.omit( unique(clusters.tp.missed[, cluster]) )
	missed.clumem[ !missed.clumem%in%tmp ]	 			<- NA
	#plot tree to file	
	tiplabels.missed									<- hivc.clu.get.tiplabels( missed.ph, missed.df )					
	if(verbose) cat(paste("write tree to file",plotfile))
	hivc.clu.plot(missed.ph, missed.clumem, file=plotfile, pdf.scaley=25, pdf.off=0, highlight.edge.of.tiplabel=c("TPw","TPm"), highlight.edge.of.tiplabel.col= c("red","blue"), cex.nodelabel=0.1 )
	hivc.clu.plot.tiplabels( seq_len(Ntip(missed.ph)), tiplabels.missed$text, tiplabels.missed$col, cex=0.12, adj=c(-0.15,0.5), add.xinch=0, add.yinch=0 )
	dev.off()
}	
######################################################################################
hivc.clu.min.transmission.cascade<- function(brlmat)
{
	if(is.matrix(brlmat))
		ans	<- .Call("hivc_clu_mintransmissioncascade", brlmat[upper.tri(brlmat)])  						
	else
		ans	<- .Call("hivc_clu_mintransmissioncascade", brlmat)
	ans		
}
######################################################################################
#' Compute a statistic of patristic distances for the subtree starting at \code{node}
#' @param node 				internal node at which the subtree starts
#' @param tree 				phylogenetic tree
#' @param distmat 			matrix of patristic distances, either precomputed or NULL
#' @param eval.dist.btw 	string, either "leaf" or "all". Specifies the nodes of the subtree on which patristic distances are considered.
#' @param stat.fun 			function, statistic to be computed from a set of patristic distances.
#' @return statistic of patristic distances		
hivc.clu.brdist.stats.subtree<- function(node, tree, distmat, eval.dist.btw="leaf", stat.fun=max)
{
	require(ape)
	require(igraph)
	require(geiger)	
	if(eval.dist.btw=="leaf")
	{
		nlist	<- tips(tree,node)
		foo 	<- distmat[nlist,nlist] 		
	}
	else if(eval.dist.btw=="all")
	{
		nlist	<- tips(tree,node)
		elist 	<- tree$edge[which.edge(tree,nlist),2]
		foo 	<- distmat[elist,elist] 	
	}
	else	
		stop("unrecognized eval.dist.btw")
	return( stat.fun(foo[upper.tri(foo,diag=FALSE)]) )
}
######################################################################################
#' Compute a statistic of patristic distances for each subtree in the tree
hivc.clu.brdist.stats<- function(tree, distmat=NULL, eval.dist.btw="leaf", stat.fun=max)
{
	require(phytools)
	if(is.null(distmat))
	{
		if(eval.dist.btw=='leaf')
			distmat	<- cophenetic.phylo(tree)
		else if(eval.dist.btw=='all')
			distmat <- dist.nodes(tree)
		else	stop("unrecognized eval.dist.btw")
	}
	ntips			<- Ntip(tree)
	nint 			<- tree$Nnode 		## number of internal nodes
	return(sapply(seq.int(ntips+1,ntips+nint), function(x) 		hivc.clu.brdist.stats.subtree(x,tree,distmat,eval.dist.btw=eval.dist.btw, stat.fun=stat.fun)		))	
}
######################################################################################
hivc.clu.droptipincluster<- function(ph, clustering, cluphy.subtrees, df.cluinfo, splitexpr)
{	
	setkey(df.cluinfo, FASTASampleCode)
	for(i in seq_along(cluphy.subtrees))
	{
		x																	<- cluphy.subtrees[[i]]
		x.cluster															<- as.numeric(names(cluphy.subtrees)[i])
		#drop selected tips from cluster					
		x.droplabel															<- subset( df.cluinfo[x$tip.label,], eval(splitexpr), FASTASampleCode )[,FASTASampleCode]
		if(length(x.droplabel)+1==length(x$tip.label))
		{			
			x.droplabel														<- x$tip.label
			x.tip.mrca														<- NA
			x																<- numeric()					#cannot do NULL which would automatically shorten 'cluphy.subtrees' and give an outofbounds error			
		}
		else	
		{
			x																<- drop.tip(x, x.droplabel, root.edge=1)			
			x.tip															<- sapply(x$tip.label, function(z) match(z, ph$tip.label) )				
			x.tip.anc														<- lapply(x.tip, function(z) Ancestors(ph, z) )
			x.tip.jnt														<- my.intersect.n( x.tip.anc )			
			x.tip.mrca														<- x.tip.anc[[1]][min(which(x.tip.anc[[1]] %in% x.tip.jnt))]							
		}
		#drop selected tips from data.table
		set(df.cluinfo, which(df.cluinfo[,FASTASampleCode%in%x.droplabel ]), "cluster", NA)												
		#drop selected tips from clustering
		clustering[["clu.mem"]][ clustering[["clu.mem"]]==x.cluster ]		<- NA				
		if(!is.na(x.tip.mrca))		
			clustering[["clu.mem"]][ c(x.tip.mrca, x.tip) ]					<- x.cluster					#WARNING: might miss inner nodes, so size.all would not be accurate
		clustering[["clu.idx"]][ x.cluster ]								<- x.tip.mrca		
		cluphy.subtrees[[i]]												<- x							
	}	
	cluphy.subtrees	<- lapply( which(sapply(cluphy.subtrees, length)>0), function(i) 	cluphy.subtrees[[i]] 	)
	
	list(cluphy.subtrees=cluphy.subtrees, clustering=clustering, cluphy.df=df.cluinfo)
}
######################################################################################
#modify clustering by splitting existing clusters based on 'splitexpr'. This is useful when additional meta-information is available.
hivc.clu.splitcluster<- function(ph, clustering, cluphy.subtrees, cluphy.df, splitexpr)
{
	require(phangorn)
	setkey(cluphy.df, FASTASampleCode)
	tmp	<- list()
	for(j in seq_along(cluphy.subtrees))		#need for loop because we change 'cluphy.df' and 'clustering'
	{
		#print(j)
		#j<- 17
		x					<- cluphy.subtrees[[j]]			
		x.tiplabel			<- x$tip.label
		x.cluster			<- cluphy.df[x.tiplabel[1],cluster][,cluster]
		x.splitlabel		<- subset( cluphy.df[x.tiplabel,], eval(splitexpr), FASTASampleCode )[,FASTASampleCode]
		#print(x); print(x.splitlabel); print(cluphy.df[x$tip.label,], )
		splitlabel	<- x.splitlabel
		x.splittrees		<- hivc.phy.splittip(x, x.splitlabel)
		#from x.splittrees, select between patient clusters
		#print(x.splittrees)
		x.clusters			<- c()
		if(length(x.splittrees))
			x.clusters		<- which( sapply( seq_along(x.splittrees), function(i) length(unique(cluphy.df[x.splittrees[[i]][["tip.label"]],Patient][,Patient]))>1 	) ) 	
		x.splittrees		<- lapply( x.clusters, function(i)	x.splittrees[[i]]	)		
		#update 'clustering' and 'cluphy.df'
		if(!length(x.splittrees))			#Splitting breaks x into singletons. In 'df', set 'cluster' to NA. In 'clustering', update 'clu.idx' and 'clu.mem'.
		{
			x.clusters		<- c()
			clustering[["clu.idx"]][ x.cluster ]							<- NA
			clustering[["clu.mem"]][ clustering[["clu.mem"]]==x.cluster ]	<- NA			
			set(cluphy.df, which(cluphy.df[,FASTASampleCode%in%x.tiplabel]), "cluster", NA)
		}
		else	
		{				
			x.clusters		<- c( x.cluster, seq.int(length(clustering[["clu.idx"]])+1, len=length(x.splittrees)-1) )
			#print("HERE"); print(x.clusters); print(x.splittrees)
			for(i in seq_along(x.splittrees))
			{
				#print(i)
				#i<- 1
				#single smaller cluster. In 'df', set 'cluster' for lost x.tiplabels to NA. In 'clustering', update 'clu.idx' and 'clu.mem'.
				x.tip																<- sapply(x.splittrees[[i]][["tip.label"]], function(z) match(z, ph$tip.label) )				
				x.tip.anc															<- lapply(x.tip, function(z) Ancestors(ph, z) )
				x.tip.jnt															<- my.intersect.n( x.tip.anc )			
				x.tip.mrca															<- x.tip.anc[[1]][min(which(x.tip.anc[[1]] %in% x.tip.jnt))]				
				clustering[["clu.mem"]][ clustering[["clu.mem"]]==x.clusters[i] ]	<- NA				
				clustering[["clu.mem"]][ c(x.tip.mrca, x.tip) ]						<- x.clusters[i]												#WARNING: might miss inner nodes, so size.all would not be accurate
				clustering[["clu.idx"]][ x.clusters[i] ]							<- x.tip.mrca													#ok to assign beyond current length
				x.tiplabel															<- setdiff(x.tiplabel, x.splittrees[[i]][["tip.label"]])
				set(cluphy.df, which(cluphy.df[,FASTASampleCode%in%x.tiplabel ]), "cluster", NA)													#set lost tip labels to NA
				set(cluphy.df, which(cluphy.df[,FASTASampleCode%in%x.splittrees[[i]][["tip.label"]]]), "cluster", x.clusters[i])					#reset tip labels to x.clusters[i], as this might be a new cluster					
			}
		}
		tmp[[j]]<- list(clu= x.splittrees, clu.names=x.clusters)					
	}
	cluphy.cluidx				<- which(sapply(seq_along(tmp), function(i)		length(tmp[[i]][["clu"]])>0	))	
	cluphy.subtrees				<- eval(parse(text=paste("c(",paste('tmp[[',seq_along(tmp),']][["clu"]]', sep='',collapse=','),")",sep='')))	
	names(cluphy.subtrees)		<- eval(parse(text=paste("c(",paste('tmp[[',seq_along(tmp),']][["clu.names"]]', sep='',collapse=','),")",sep='')))	
	clustering[["size.all"]]	<- NULL
	clustering[["size.tips"]]	<- table( clustering[["clu.mem"]][1:Ntip(ph)] )
	
	list(cluphy.subtrees=cluphy.subtrees, clustering=clustering, cluphy.df=cluphy.df)
}	
######################################################################################
hivc.phy.splittip<- function(x, splitlabel)
{
	#str(x)
	x.splittip			<- which(x$tip.label%in%splitlabel)								#take first matching label as this 'x.splittip' and put rest on stack
	if(!length(x.splittip))		
		return(list(x))
	x.splittip			<- x.splittip[1]
	splitlabel.stack	<- setdiff(splitlabel, x$tip.label[ x.splittip ])	
	#print(x$tip.label); print(x.splittip); print(splitlabel.stack)	
	if(Ntip(x)<3)																		#tip is among tiplabels, so splitting will break pair
		return(list())
	x.splitnode			<- x$edge[x$edge[,2]==x.splittip, 1]							#ancestor of HETF tip		
	#x$tip.label[x.splittip]
	x.splittree2		<- extract.clade(x, x.splitnode, root.edge = 1)
	
	#replace 'x.splitnode' with 'x.splittip' and drop dummy tip	
	x.rmnode							<- x.splitnode
	tmp									<- x$edge[x$edge[,1]==x.splitnode,2] 			
	while(any(tmp>Ntip(x)))
	{
		x.rmnode	<- c(x.rmnode, tmp[ tmp>Ntip(x) ])
		tmp			<- x$edge[x$edge[,1]%in%tmp,2]
	}
	x$node.label						<- x$node.label[ -(x.rmnode-Ntip(x)) ]
	tmp									<- x$edge[,1]%in%x.rmnode
	x$edge.length						<- x$edge.length[ !tmp ]	
	x.rmtip								<- x$edge[tmp,2]
	x.rmtip								<- x.rmtip[x.rmtip!=x.splittip & x.rmtip<=Ntip(x)]	 
	x$tip.label							<- x$tip.label[-x.rmtip]
	x$Nnode								<- x$Nnode-length(x.rmnode)	
	x$edge								<- x$edge[!x$edge[,1]%in%x.rmnode,]				#drop clade
	x.splitedge							<- x$edge[,2]==x.splitnode						#find edge leading to clade	
	x$edge[x.splitedge,2]				<- x.splittip									#set edge to dummy tip
	if(nrow(x$edge))
	{
		tmp								<- as.vector(x$edge)							#need to get node numbers right again
		tmp.rm							<- which(!seq_len(max(tmp))%in%tmp)
		for(z in rev(tmp.rm))
			tmp[tmp>=z]					<- tmp[tmp>=z]-1			
		x$edge							<- matrix(tmp, nrow=nrow(x$edge), ncol=ncol(x$edge) )
		x.splittip						<- x$edge[x.splitedge,2]						#renumbering might have changed dummy tip
	}
	else
	{
		x.splittip						<- 1
	}
	ans							<- list()
	#print("HERE"); str(x); print(x.splittip); str(x.splittree2); print(x$tip.label[x.splittip]); print("HERE2")
	if(Ntip(x)>2)
		ans[[1]]				<- drop.tip(x, x.splittip, root.edge=1)
	if(Ntip(x.splittree2)>2)
		ans[[length(ans)+1]]	<- drop.tip(x.splittree2, x$tip.label[x.splittip], root.edge=1)
	if(length(splitlabel.stack) && length(ans))
	{
		tmp	<- paste("c(",paste('hivc.phy.splittip(ans[[',seq_along(ans),']], splitlabel.stack)', sep='',collapse=','),")",sep='')		
		tmp	<- eval(parse(text=tmp))																	#annoyingly, sapply etc is not the same as c( etc )		
		#tmp	<- sapply(seq_along(ans), function(i)	hivc.phy.splittip(ans[[i]], splitlabel.stack) )		#c(list()) is [[1]] list(), need a hack
		#ans	<- list()
		#ans	<- lapply( which(sapply(tmp, length)>0), function(i) tmp[[i]] )
		ans <- tmp
	}
	ans
}		
######################################################################################
hivc.phy.get.TP.and.TN<- function(ph, df.all, verbose= 1)
{
	#
	# exctract geographically distant seqs that are assumed to be truly unlinked to NL seqs
	#
	tmp					<-  ph$tip.label[ substr(ph$tip.label,1,2)=="TN" ]
	unlinked.byspace	<- data.table(FASTASampleCode=tmp, key="FASTASampleCode")		
	#
	# extract list of truly linked sample codes
	#
	tmp					<- subset(df.all[, length(FASTASampleCode)>1, by=Patient], V1==T, select=Patient)
	linked.bypatient	<- merge(tmp, df.all, all.x=1, by="Patient")
	setkey(linked.bypatient, "FASTASampleCode")
	linked.bypatient	<- subset(linked.bypatient, select=c(Patient, FASTASampleCode, PosSeqT))					
	#
	# extract list of truly unlinked sample codes by temporal separation
	#
	df.dead						<- subset(df.all, !is.na(DateDied))	
	setkey(df.dead, DateDied)	
	#extract seroconverters
	df.serocon.acc				<- subset(df.all, NegT_Acc=="Yes" & NegT>=df.dead[1,DateDied],)
	#add seroconverters with inaccurate info
	df.serocon.nacc				<- subset(df.all, NegT_Acc=="No" & !is.na(NegT) & NegT>=df.dead[1,DateDied], )
	df.serocon.nacc.dy			<- subset(df.serocon.nacc, as.POSIXlt(NegT)$mday==15, )						#for inaccurate days, we (conservatively) assume the patient was only seronegative at the start of the month
	tmp							<- as.POSIXlt(df.serocon.nacc.dy[,NegT] )
	tmp$mday					<- 1
	df.serocon.nacc.dy[,NegT:=as.Date(tmp)]
	df.serocon.nacc.mody		<- subset(df.serocon.nacc, as.POSIXlt(NegT)$mon==6 & as.POSIXlt(NegT)$mday==1, )		#for inaccurate months and days, we (conservatively) assume the patient was only seronegative at the start of the year
	tmp							<- as.POSIXlt(df.serocon.nacc.mody[,NegT] )
	tmp$mon						<- 0
	df.serocon.nacc.mody[,NegT:=as.Date(tmp)]
	df.serocon					<- rbind(df.serocon.acc, df.serocon.nacc.dy, df.serocon.nacc.mody)		
	if(verbose) cat(paste("\nnumber of seroconverters with at least 1 preceeding dead HIV+",nrow(df.serocon)))
	#for each seq of seroconverter, extract HIV+ seqs that are dead before seroconversion
	if( any(as.logical(df.serocon[,is.na(NegT)])) )	warning("Found accurate seroconverters with missing NegT")		
	setkey(df.serocon,NegT)
	unlinked.bytime				<- lapply(seq_len(nrow(df.serocon)), function(i)
				{
					tmp				<- subset(df.dead, DateDied<=df.serocon[i,NegT],select=c(Patient,FASTASampleCode,DateDied))
					tmp2			<- rep(df.serocon[i,FASTASampleCode],nrow(tmp))
					tmp[,query.FASTASampleCode:= tmp2]
					setkey(tmp,"FASTASampleCode")
					tmp					
				})
	names(unlinked.bytime)	<- df.serocon[,FASTASampleCode]											
	#
	#plot number of unlinked HIV+ by date (this date is when someone else is still seroneg)
	#
	if(0)
	{
			#
			y		<- sapply(unlinked.bytime, function(x) nrow(x) )
			y2		<- nrow(unlinked.byspace)
			x		<- df.serocon[,NegT]
			xlim	<- range(x)
			ylim	<- c(0, y2+max(y))
			tmp		<- as.POSIXlt(xlim[1])
			tmp$mday<- 1
			tmp$mon	<- 1
			xlim[1]	<- as.Date(tmp)
			plot(1,1,type='n',bty='n',xlab="time of seroconversion",ylab="number unlinked",xaxt='n',xlim=xlim,ylim=ylim)
			axis.Date(1,Year,at=seq(xlim[1], xlim[2],by="12 months"))			
			polygon(c(x[1],x[length(x)],x[length(x)],x[1]), c(y2,y2,0,0), border=NA, col="blue" )
			polygon(c(x,c(x[length(x)],x[1])), c(y+y2,y2,y2), border=NA, col="grey60" )
			legend("topleft",bty='n',fill=c("blue","grey60"),legend=c("by space","by time"),border=NA)						
	}	
	#
	# get TN for seroneg: convert unlinked.bytime from SampleCode to ph node index and add truly unlinked.byspace
	#
	df.tips							<- data.table(Node=seq_len(Ntip(ph)), FASTASampleCode=ph$tip.label, key="FASTASampleCode" )
	unlinked.byspace				<- merge(unlinked.byspace, df.tips, all.x=1, by="FASTASampleCode")
	setkey(unlinked.byspace,"Node")
	ph.unlinked.bytime				<- lapply(unlinked.bytime, function(x)
				{
					tmp	<- merge(x, df.tips, all.x=1)
					tmp2<- tmp[1,query.FASTASampleCode]
					tmp	<- rbind( subset(tmp, select=c(FASTASampleCode,Node)), unlinked.byspace )
					tmp2<- rep(tmp2, nrow(tmp))
					tmp[,query.FASTASampleCode:= tmp2]
					setkey(tmp, Node)
					tmp
				})	
	names(ph.unlinked.bytime)		<- names(unlinked.bytime)
	#
	# get TN for seropos: use only unlinked.byspace
	#
	ph.seroneg						<- merge( data.table(FASTASampleCode=names(ph.unlinked.bytime)), df.tips, by="FASTASampleCode", all.x=1)		
	ph.seropos						<- subset( df.tips[,c(2,1), with=F], !FASTASampleCode%in%names(ph.unlinked.bytime) 
														& 	substr(FASTASampleCode,1,2)!="TN"
														&	substr(FASTASampleCode,1,8)!="PROT+P51"
														)			
	ph.unlinked.byspace				<- lapply(seq_len(nrow(ph.seropos)), function(i)
			{
				tmp	<- unlinked.byspace
				tmp2<- rep(ph.seropos[i,FASTASampleCode], nrow(tmp))
				tmp[,query.FASTASampleCode:= tmp2]
				setkey(tmp, Node)
				tmp
			})
	names(ph.unlinked.byspace)		<- ph.seropos[,FASTASampleCode]
	#
	# put all unlinked ATHENA seqs together and sort by Node
	#
	ph.unlinked						<- c(ph.unlinked.byspace, ph.unlinked.bytime)
	names(ph.unlinked)				<- c(names(ph.unlinked.byspace),names(ph.unlinked.bytime))
	ph.unlinked.info				<- rbind(ph.seropos,ph.seroneg)		
	setkey(ph.unlinked.info, "Node")
	ph.unlinked						<- lapply(ph.unlinked.info[,FASTASampleCode], function(i){		ph.unlinked[[i]]	})
	names(ph.unlinked)				<- ph.unlinked.info[,FASTASampleCode]
	tmp								<- seq_len(nrow(ph.unlinked.info))
	ph.unlinked.info[,NodeIdx:= tmp]
	ph.unlinked.info				<- merge(ph.unlinked.info, df.all, all.x=1, by="FASTASampleCode")
	setkey(ph.unlinked.info, "Node")		
	#
	#	get data table of all linked ATHENA seqs
	#
	ph.linked						<- merge(linked.bypatient, df.tips, all.x=1)
	ph.linked						<- subset(ph.linked, select=c(FASTASampleCode, Patient, Node))
	ph.linked						<- merge(ph.linked, ph.linked[, length(FASTASampleCode), by=Patient], all.x=1, by="Patient")
	setnames(ph.linked, "V1","nTP")
	setkey(ph.linked, "Node")
	ans								<- list(	unlinked.byspace=unlinked.byspace, unlinked.bytime=unlinked.bytime, linked.bypatient=linked.bypatient,
												ph.linked=ph.linked, ph.unlinked.info=ph.unlinked.info, ph.unlinked=ph.unlinked)
	ans									
}
######################################################################################
hivc.get.hpcsys<- function()
{
	tmp<- system('domainname',intern=T)
	if(!nchar(tmp))	tmp<- "debug"
	tmp
}

