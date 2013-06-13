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

hivx.seq.find<- function(char.matrix, pos0= NA, from= c(), verbose=1)
{
	if(is.na(pos0)) 	stop("start position of token to be replaced is missing")
	if(!length(from))	stop("token to be replaced is missing")
	query.colidx	<- seq.int(pos0,pos0+length(from)-1)
	query.yes		<- which( apply(char.matrix, 1, function(x)	all(x[query.colidx]==from) ) )
	query.yes	
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
	## set up clustering
	ntips		<- Ntip(ph)	
	clu.i 		<- 0 ## cluster number
	clu.mem 	<- rep(NA,ntips+ph$Nnode) ## cluster member assignment
	clu.idx		<- rep(NA,ntips+ph$Nnode) ## cluster index assignment
	igraph.ph	<- graph.edgelist(ph$edge) ## ph in igraph form
	dfs 		<- graph.dfs(igraph.ph,root=ntips+1,neimode='out',order=TRUE,dist=TRUE)
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
hivc.clu.typeIerror<- function(clustering, ph.unlinked.seroneg, ph.unlinked.dead, verbose=0)
{
	clusters	<- unique(clustering[["clu.mem"]][!is.na(clustering[["clu.mem"]])])
	clu.fp.n	<- sapply(clusters, function(clu)		
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
				length(which( latest.seroneg.unlinked%in%clu.leaves ))																	
			})		
	list(clu.n=length(clusters), fp.n= length(which(clu.fp.n>0)), fp.unlinked.n= sum(clu.fp.n), fp.rate= length(which(clu.fp.n>0))/length(clusters))
}
######################################################################################
hivc.clu.clusterbytypeIerror<- function(ph, dist.brl, nodesupport, ph.unlinked.seroneg, ph.unlinked.dead, thresh.brl=NULL, thresh.bs=NULL, level= 0.01, tol= level, mxit= 100, thresh.bs.lower= min(nodesupport), thresh.brl.upper=max(dist.brl), verbose=0)
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
	while(abs(error)>tol && nit<mxit)
	{	
		if(nit)
		{
			if(srch.find.thresh.bs)	
				thresh.bs			<- ( srch.upper + srch.lower ) / 2
			else
				thresh.brl			<- ( srch.upper + srch.lower ) / 2				
		}
			
		clustering					<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=nodesupport, retval="all")
		tmp							<- hivc.clu.typeIerror(clustering, ph.unlinked.seroneg, ph.unlinked.dead, verbose=0)		
		if(verbose)	
			cat(paste("\nit=",nit,", #clu=",tmp[["clu.n"]],", BS=",thresh.bs,", BRL=",thresh.brl,", #unlinked in any cluster=",tmp[["fp.unlinked.n"]],", #FP=",tmp[["fp.n"]],", %FP=",round(tmp[["fp.rate"]],d=5),sep=''))			
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
	ans<- list(thresh.bs=thresh.bs, thresh.brl=thresh.brl, clu=clustering, clu.fp.n=clu.fp.n, srch.nit= nit, srch.error=error, srch.tol=tol)
	ans		
}
######################################################################################
hivc.clu.plot<- function(ph, clu, show.node.label= T, file=NULL, pdf.scaley=10, pdf.width= 7, pdf.height=pdf.scaley*7)
{
	require(colorspace)	
	clu.edge							<- clu[ ph$edge[,1] ]
	clu.edge							<- clu.edge+1				#set col index for non-clustering edges to 1
	clu.edge[is.na(clu.edge)]			<- 1
	cols.n								<- length(unique(clu.edge))-1	
	#cols.palette 						<- colorRampPalette(c("red", "orange", "blue"),space = "Lab")
	#cols								<- c("black",cols.palette(cols.n))
	cols								<- c("black",rainbow_hcl(cols.n, start = 30, end = 300))	
	clu.edge.col						<- cols[clu.edge]
	clu.edge.width						<- rep(1, length(clu.edge))
	clu.edge.width[clu.edge!=1]			<- 2
	if(class(file)=="character")
		pdf(file,width=pdf.width,height=pdf.height)
	plot(ph, show.tip.label=F, show.node.label=show.node.label, cex=0.5, edge.color=clu.edge.col, edge.width=clu.edge.width)
	if(class(file)=="character")
		dev.off()
	#tiplabels(rep("  ",Ntip(ph)), bg=clu)	
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
hivc.clu.brdist.stats<- function(tree, distmat=NULL, eval.dist.btw="leaf")
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
	return(sapply(seq.int(ntips+1,ntips+nint), function(x) 		hivc.clu.brdist.stats.subtree(x,tree,distmat,eval.dist.btw=eval.dist.btw)		))	
}
######################################################################################
hivc.get.hpcsys<- function()
{
	tmp<- system('domainname',intern=T)
	if(!nchar(tmp))	tmp<- "debug"
	tmp
}

