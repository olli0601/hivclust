
######################################################################################
#' @export
seq.write.dna.nexus<- function(seq.DNAbin.mat, ph=NULL, file=NULL, nexus.format="DNA",nexus.gap='-', nexus.missing='?', nexus.interleave="NO")
{		
	tmp		<- cbind( rownames(seq.DNAbin.mat), apply( as.character( seq.DNAbin.mat ), 1, function(x) paste(x,sep='',collapse='')  ) )
	tmp		<- apply(tmp, 1, function(x) paste(x, collapse='\t', sep=''))
	tmp		<- paste(tmp, collapse='\n',sep='')
	header	<- paste( "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=",nrow(seq.DNAbin.mat)," NCHAR=",ncol(seq.DNAbin.mat),";\nFORMAT DATATYPE=",nexus.format," MISSING=",nexus.missing," GAP=",nexus.gap," INTERLEAVE=",nexus.interleave,";\nMATRIX\n", collapse='',sep='')
	tmp		<- paste(header, tmp, "\n;\nEND;\n", sep='')
	if(!is.null(ph))
	{
		tmp	<- paste(tmp,"#BEGIN TAXA;\nTAXLABELS ", paste(ph$tip.label, sep='',collapse=' '), ';\nEND;\n\n',sep='')
		tmp	<- paste(tmp, "BEGIN TREES;\nTREE tree1 = ", write.tree(ph), "\nEND;\n", sep='')
	}
	if(!is.null(file))
		cat(tmp, file=file)
	tmp
}		
######################################################################################
#' @export
seq.write.dna.phylip<- function(seq.DNAbin.mat, file)
{		
	tmp<- cbind( rownames(seq.DNAbin.mat), apply( as.character( seq.DNAbin.mat ), 1, function(x) paste(x,sep='',collapse='')  ) )
	tmp<- paste(t(tmp),collapse='\n',sep='')	
	tmp<- paste( paste(c(nrow(seq.DNAbin.mat),ncol(seq.DNAbin.mat)),sep='',collapse=' '),'\n',tmp,'\n',collapse='',sep='' )
	cat(tmp, file=file)
}
######################################################################################
seq.singleton2bifurcatingtree<- function(ph.s, dummy.label=NA)
{	
	if(!Nnode(ph.s))
	{
		stopifnot(!nrow(ph.s$edge))
		if(is.na(dummy.label))
			dummy.label		<- paste('DUMMY',ph.s$tip.label, sep='_')
		ph.s$edge			<- matrix(c(3,1,3,2), nrow=2, ncol=2, byrow=TRUE) 	
		ph.s$edge.length	<- c(0,0)
		ph.s$tip.label		<- c(ph.s$tip.label, dummy.label)
		ph.s$Nnode			<- 1
	}
	ph.s
}
######################################################################################
seq.find<- function(seq, pos0= NA, from= c(), verbose=1)
{
	if(is.na(pos0)) 	stop("start position of token to be replaced is missing")
	if(!length(from))	stop("token to be replaced is missing")
	query.colidx	<- seq.int(pos0,pos0+length(from)-1)
	if(class(seq)=='matrix')	
		query.yes		<- which( apply(seq, 1, function(x)	all(x[query.colidx]==from) ) )
	else if(class(seq)=='list')
		query.yes		<- which( sapply(seq, function(x)	all(x[query.colidx]==from) ) )
	query.yes	
}
######################################################################################
seq.length<- function(seq.DNAbin.mat, exclude=c('-','?'))
{
	counts	<- sapply(seq_len(nrow(seq.DNAbin.mat)), function(i) base.freq(seq.DNAbin.mat[i,], freq=1, all=1))
	apply(counts[ !rownames(counts)%in%exclude, ],2,sum)
}
######################################################################################
seq.proportion.ambiguous<- function(seq.DNAbin.mat, exclude=c('-','?'))
{
	counts	<- apply(seq.DNAbin.mat,1,function(x) base.freq(x, freq=1, all=1))
	len		<- apply(counts[ !rownames(counts)%in%exclude, ],2,sum)
	pa		<- apply(counts[c("r", "m", "w", "s", "k", "y", "v", "h", "d", "b"),],2,sum)
	pa/len
}
######################################################################################
seq.gc.content<- function(seq.DNAbin.mat)
{	
	rna.gc.fraction.n		<- c('a','c','g','t',	'r','m','w','s',	'k','y','v','h',		'd','b','n','-','?')
	rna.gc.fraction			<- c( 0, 1, 1, 0,		0.5, 0.5, 0, 1, 	1/2, 1/2, 2/3, 1/3,		1/3,2/3, 1/4, 0, 0)		#this fraction assumes that U is synonymous with T and that U does not occur in the code
	rna.gc.sum				<- c( T, T, T, T,       T, T, T, T,         T, T, T, T,				T, T, T, F, F )	
	counts					<- apply(seq.DNAbin.mat,1,function(x) base.freq(x, freq=1, all=1))
	apply(counts*rna.gc.fraction,2,sum) / apply(counts[rna.gc.sum,],2,sum)
}
######################################################################################
#slight modification of blastSequences() in pkg annotate
seq.blast<- function (x, database = "nr", hitListSize = "10", filter = "L", expect = "10", program = "blastn", organism= "HIV-1") 
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
######################################################################################
#slight modification of read.blast() in pkg RFLPtools; expects blast was run with -outfmt 6
seq.blast.read<- function (file, sep = "\t") 
{
	require(data.table)
	x <- read.table(file = file, header = FALSE, sep = sep, quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
	if (ncol(x) != 12) 
		stop("Data in given file", basename(file), "has wrong dimension!")
	names(x) <- c("query.id", "subject.id", "identity", "alignment.length","mismatches", "gap.opens", "q.start", "q.end", "s.start","s.end", "evalue", "bit.score")
	data.table(x, key="query.id")
}
######################################################################################
#' @export
seq.rm.drugresistance<- function(char.matrix, dr, verbose=1, rtn.DNAbin=1)
{
	if(verbose)	cat(paste("\nchecking for drug resistance mutations, n=",nrow(dr)))
	tmp	<- rep(0, nrow(dr))
	for(i in seq_len(nrow(dr)))
	{		
		query.yes	<- seq.find(char.matrix, dr[i,Alignment.nuc.pos], unlist(strsplit(unlist(dr[i,Mutant.NTs]),'')))
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
######################################################################################
#' @export
seq.unique<- function(seq.DNAbin.matrix)
{
	x<- as.character(seq.DNAbin.matrix)
	x<- apply(x, 1, function(z) paste(z,collapse=''))
	seq.DNAbin.matrix[!duplicated(x),]			
}
######################################################################################
#' @export
seq.dist.pairwise<- function(x, y)
{
	stopifnot(class(x)=='DNAbin', class(y)=='DNAbin', length(dim(x))==2, length(dim(y))==2, nrow(x)==1, nrow(y)==1)
	dummy		<- 0
	1-.C("hivc_dist_ambiguous_dna", x, y, ncol(x), dummy )[[4]]
}
######################################################################################
#' @export
seq.dist<- function(seq.DNAbin.matrix, verbose=1)
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
		if(nrow(seq.DNAbin.matrix)>5)
		{
			ans[1,1:6]<- 2^seq.int(6,11)-1
			if(is.na(ans[1,2]))		stop("unexpected behaviour of bigmemory")
			ans[1,1:6]<- NA
		}						
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
######################################################################################
#' Collapse singleton nodes (bugfix to collapse.singles in ape)
#' @export
seq.collapse.singles<- function (tree) 
{
	elen 		<- tree$edge.length
	xmat 		<- tree$edge
	node.lab 	<- tree$node.label
	nnode 		<- tree$Nnode
	ntip 		<- length(tree$tip.label)
	root		<- 0
	singles 	<- NA
	while (length(singles) > 0) 
	{
		tx <- tabulate(xmat[, 1])
		singles <- which(tx == 1)
		if (length(singles) > 0) 
		{
			i 					<- singles[1]
			prev.node 			<- which(xmat[, 2] == i)
			next.node 			<- which(xmat[, 1] == i)
			xmat[prev.node, 2] 	<- xmat[next.node, 2]
			xmat 				<- xmat[xmat[, 1] != i, , drop=0]
			xmat[xmat > i] 		<- xmat[xmat > i] - 1L
			if(!length(prev.node))
				root			<- root + elen[next.node]
			if(length(prev.node))
				elen[prev.node] <- elen[prev.node] + elen[next.node]
			if (!is.null(node.lab)) 
				node.lab <- node.lab[-c(i - ntip)]
			nnode <- nnode - 1L
			elen <- elen[-next.node]
		}
	}
	tree$edge 			<- xmat
	tree$edge.length 	<- elen
	tree$node.label 	<- node.lab
	tree$Nnode 			<- nnode
	tree$root.edge		<- root
	tree
}
######################################################################################
#' Reads newick file with singletons and node labels (extends read.newick in phytools)
#' @export
seq.read.newick<- function (file = "", text) 
{
	if (file != "") 
		text <- scan(file, sep = "\n", what = "character")	
	Nnode		<- length(gregexpr("\\(", text)[[1]])
	Ntip		<- 1+length(gregexpr(",", text)[[1]])	
	tree 		<- unlist(strsplit(text, NULL))
	tip.label 	<- vector(mode = "character")
	node.label	<- vector(mode = 'character')
	edge 		<- matrix(data = 0, Nnode + Ntip - 1, 2)
	edge.length <- rep(0, Nnode + Ntip - 1)
	ntip 		<- vector(mode = "numeric")
	currnode 	<- Ntip + 1
	nodecount 	<- currnode
	i 	<- 1
	j 	<- 1
	k 	<- 1
	while(tree[i] != ";") 
	{
		if(tree[i] == "(") 
		{
			edge[j, 1] <- currnode
			i <- i + 1
			if(is.na(match(tree[i], c("(", ")", ",", ":", ";")))) 
			{
				l				<- gregexpr(",|:|\\)", substr(text, i, nchar(text)))[[1]][1]
				stopifnot(l>0)
				tip.label[k] 	<- substr(text, i, i+l-2)	
				i				<- i+l-1
				edge[j, 2] 		<- k
				k 				<- k + 1
				ntip[j] 		<- 1
				if (tree[i] == ":") 
				{
					i				<- i + 1
					l				<- gregexpr(",|\\)", substr(text, i, nchar(text)))[[1]][1]
					stopifnot(l>0)
					edge.length[j]	<- as.numeric(substr(text, i, i+l-2))
					i				<- i+l-1					
				}
			}
			else if(tree[i] == "(") 
			{
				nodecount 	<- nodecount + 1
				currnode 	<- nodecount
				edge[j, 2] 	<- currnode
			}
			j <- j + 1
		}
		else if(tree[i] == ")") 
		{
			i <- i + 1			
			if(is.na(match(tree[i], c(":", ")")))) 
			{
				l	<- gregexpr(":|;|\\)", substr(text, i, nchar(text)))[[1]][1]
				stopifnot(l>0)
				node.label[currnode-Ntip]	<- substr(text, i, i+l-2)
				i	<- i+l-1				
			}
			if(tree[i] == ":") 
			{
				i 	<- i + 1
				l	<- gregexpr(",|\\)", substr(text, i, nchar(text)))[[1]][1]
				stopifnot(l>0)
				edge.length[match(currnode, edge[, 2])] <- as.numeric(substr(text, i, i+l-2))
				i	<- i+l-1	
			}
			ntip[match(currnode, edge[, 2])] 	<- sum(ntip[which(edge[, 1] == currnode)])
			currnode 							<- edge[match(currnode, edge[, 2]), 1]
		}
		else if(tree[i]==",") 
		{
			edge[j, 1] 	<- currnode
			i 			<- i + 1
			if(is.na(match(tree[i], c("(", ")", ",", ":", ";")))) 
			{
				l				<- gregexpr(",|:|\\)", substr(text, i, nchar(text)))[[1]][1]
				stopifnot(l>0)
				tip.label[k] 	<- substr(text, i, i+l-2)	
				i				<- i+l-1
				edge[j, 2] 		<- k
				k 				<- k + 1
				ntip[j] 		<- 1
				if (tree[i] == ":") 
				{
					i				<- i + 1
					l				<- gregexpr(",|\\)", substr(text, i, nchar(text)))[[1]][1]
					stopifnot(l>0)
					edge.length[j]	<- as.numeric(substr(text, i, i+l-2))
					i				<- i+l-1
				}
			}
			else if (tree[i] == "(") 
			{
				nodecount 	<- nodecount + 1
				currnode 	<- nodecount
				edge[j, 2] 	<- currnode
			}
			j <- j + 1
		}
	}
	tmp	<- which(edge[,1]==0)
	if(length(tmp))
	{
		edge		<- edge[-tmp,]
		edge.length	<- edge.length[-tmp]		
		tmp			<- sort( unique( as.numeric( edge ) ) )
		tmp			<- rbind(tmp, seq_along(tmp))		
		tmp			<- sapply( as.numeric( edge ), function(j)	tmp[2, match(j, tmp[1,])] )
		edge		<- matrix(tmp, ncol=2)
	}
	phy <- list(edge = edge, Nnode = as.integer(Nnode), tip.label = tip.label, Ndesc = ntip)
	if(sum(edge.length) > 1e-08) 
		phy$edge.length	<- edge.length
	if(length(node.label))
		phy$node.label	<- node.label
	class(phy) <- "phylo"
	return(phy)
}
######################################################################################
#' @export
seq.read.GenBank<- function (access.nb, seq.names = access.nb, species.names = TRUE, gene.names = FALSE, as.character = FALSE, attributes= c("origin")) 
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
######################################################################################
#' @export
seq.dist.dna<- function (x, model = "K80", variance = FALSE, gamma = FALSE, 
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
######################################################################################
#' @export
seq.create.referencepairs<- function(dir.name= DATA)
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
######################################################################################
seq.replace<- function(seq.DNAbin.matrix, code.from='?', code.to='n', verbose=0)
{
	seq.DNAbin.matrix	<- as.character(seq.DNAbin.matrix)		
	seq.DNAbin.matrix	<- apply(seq.DNAbin.matrix, 2, function(col) 		gsub(code.from,code.to,col,fixed=1)			)	
	as.DNAbin( seq.DNAbin.matrix )
}
######################################################################################
#' @export
seq.rmgaps<- function(seq.DNAbin.matrix, rm.only.col.gaps=1, rm.char='-', verbose=0)
{
	if(class(seq.DNAbin.matrix)=='DNAbin')
		seq.DNAbin.matrix		<- as.character(seq.DNAbin.matrix)		
	if(!rm.only.col.gaps)
	{	
		if(is.matrix(seq.DNAbin.matrix))
		{
			tmp					<- lapply(seq_len(nrow(seq.DNAbin.matrix)), function(i){	seq.DNAbin.matrix[i, !seq.DNAbin.matrix[i,]%in%rm.char]	})
			names(tmp)			<- rownames(seq.DNAbin.matrix)
		}
		else
		{
			tmp					<- lapply(seq_along(seq.DNAbin.matrix), function(i){	seq.DNAbin.matrix[[i]][ !seq.DNAbin.matrix[[i]]%in%rm.char]	})
			names(tmp)			<- names(seq.DNAbin.matrix)
		}		
		seq.DNAbin.matrix	<- tmp
	}
	else
	{		
		gap					<- apply(seq.DNAbin.matrix,2,function(x) all(x%in%rm.char)) 
		if(verbose)		cat(paste("\nremove gaps, n=",length(which(gap))))
		if(verbose>1)	cat(paste("\nremove gaps, at pos=",which(gap)))
		seq.DNAbin.matrix	<- seq.DNAbin.matrix[,!gap]	
	}
	as.DNAbin( seq.DNAbin.matrix )
}
######################################################################################
#' @export
seq.rmallchar<- function(seq, rm.char='-', verbose=0)
{
	if(class(seq)=='DNAbin')
		seq		<- as.character(seq)	
	tmp			<- apply(seq, 2, function(x) all(x==rm.char)) 	
	if(verbose)	cat(paste("\nremove gaps, n=",length(which(tmp))))
	as.DNAbin( seq[, !tmp] )	
}