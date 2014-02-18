
######################################################################################
#	get BEAST taxon labels:		cluster	FASTASampleCode	NegT	AnyPosT	SeqT  -> turn date into numerical format
hivc.beast.addBEASTLabel<- function( df, df.resetTipDate=NA )
{
	if(is.na(df.resetTipDate))
		df[,dummy:=NA]
	else if(df.resetTipDate=="LdTd")
		df	<- merge(df, df[,list(dummy=max(AnyPos_T1)),by="cluster"])
	else if(df.resetTipDate=="LsTd")
		df	<- merge(df, df[,list(dummy=max(PosSeqT)),by="cluster"])
	else if(df.resetTipDate=="UmTd")
		df[,dummy:=max(df[,PosSeqT])]
	else
		stop("Unexpected df.resetTipDate")
	df	<- merge(df, df[,list(PosSeqTadj= max(PosSeqT, dummy, na.rm=T)),by="FASTASampleCode"], by="FASTASampleCode")
	
	tmp	<- df[,		{
				z	<- as.POSIXlt(c(NegT, AnyPos_T1, PosSeqT, PosSeqTadj))
				tmp	<- z$year + 1900
				z	<- tmp + round( z$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )												
				list(BEASTlabel= paste(c(cluster, z, FASTASampleCode), collapse='_', sep=''))
			}, by="FASTASampleCode"]
	df	<- merge(df, tmp, by="FASTASampleCode")
	subset(df,select=-which(colnames(df)=="dummy"))	
}
######################################################################################
#	pool clusters into sets containing roughly 'pool.ntip' sequences
hivc.beast.poolclusters<- function(cluphy.df, pool.ntip= 130, pool.includealwaysbeforeyear=NA, verbose=1)
{	
	df			<- cluphy.df[, list(clu.ntip=clu.ntip[1], clu.AnyPos_T1=clu.AnyPos_T1[1]), by="cluster"]	
	if(!is.na(pool.includealwaysbeforeyear))
	{
		if(verbose) cat(paste("\nalways include clusters starting before ",pool.includealwaysbeforeyear,"and then pool evenly across clu.AnyPos_T1"))
		df.always	<- subset(df,as.POSIXlt(clu.AnyPos_T1)$year<(pool.includealwaysbeforeyear-1900))
		df			<- subset(df,as.POSIXlt(clu.AnyPos_T1)$year>=(pool.includealwaysbeforeyear-1900))
		pool.ntip	<- pool.ntip - sum(df.always[,clu.ntip])		
		setkey(df, clu.AnyPos_T1)
		pool.n		<- ceiling( sum( df[,clu.ntip] ) / pool.ntip )
		tmp			<- lapply( seq_len(pool.n), function(x)	seq.int(x,nrow(df),by=pool.n) )		
		pool.df		<- lapply(seq_along(tmp), function(i) merge(subset(rbind(df.always, df[tmp[[i]],]), select=cluster), cluphy.df, by="cluster") )		
	}
	else
	{
		if(verbose) cat(paste("\npool evenly across clu.AnyPos_T1"))
		setkey(df, clu.AnyPos_T1)
		pool.n	<- ceiling( sum( df[,clu.ntip] ) / pool.ntip )
		tmp		<- lapply( seq_len(pool.n), function(x)	seq.int(x,nrow(df),by=pool.n) )
		pool.df	<- lapply(seq_along(tmp), function(i) merge(subset(df[tmp[[i]],], select=cluster), cluphy.df, by="cluster") )		
	}
	if(verbose) cat(paste("\nnumber of pools is n=",pool.n))		
	if(verbose) cat(paste("\nnumber of seq in pools is n=",paste( sapply(pool.df, nrow), sep='', collapse=', ' )))
	list(pool.df=pool.df, pool.ntip=pool.ntip)
}	
######################################################################################
#	add taxa and alignment in bxml from BEASTlabels in df and alignment in seq.PROT.RT
#	beast.label.datepos= 4; beast.label.sep= '_'; beast.date.direction= "forwards"; beast.date.units= "years"; beast.alignment.dataType= "nucleotide"; xml.resetTipDate2LastDiag=1
hivc.beast.add.seq<- function(bxml, df, seq.PROT.RT, beast.label.datepos= 4, beast.label.sep= '_', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.dataType= "nucleotide", verbose=1)
{			
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	
	#	add sequence taxa
	tmp.label	<- df[,BEASTlabel]	
	if(verbose)	cat(paste("\nsetting tip date to time at pos x in label, x=", beast.label.datepos))
	tmp.date	<- sapply( strsplit(tmp.label, beast.label.sep, fixed=1), function(x) x[beast.label.datepos] )
	dummy		<- newXMLCommentNode(text="The list of taxa to be analysed (can also include dates/ages).", parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(df),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	seqtaxa		<- newXMLNode("taxa", attrs= list(id="taxa"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_along(tmp.label), function(i)
			{
				taxon	<- newXMLNode("taxon", attrs= list(id=tmp.label[i]), parent=seqtaxa, doc=bxml, addFinalizer=T )
				dummy	<- newXMLNode("date", attrs= list(value=tmp.date[i], direction=beast.date.direction, units=beast.date.units), parent=taxon, doc=bxml, addFinalizer=T )
				taxon
			})	
	if(verbose)	cat(paste("\nadded new seq taxa, n=", xmlSize(seqtaxa)))
	#	add alignment	
	dummy		<- newXMLCommentNode(text="The sequence alignment (each sequence refers to a taxon above).", parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(df)," nchar=",ncol(seq.PROT.RT),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	seqalign	<- newXMLNode("alignment", attrs= list(id="alignment", dataType=beast.alignment.dataType), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- lapply( seq_len(nrow(df)), function(i)
			{
				seq		<- newXMLNode("sequence", parent=seqalign, doc=bxml, addFinalizer=T)
				dummy	<- newXMLNode("taxon", attrs= list(idref= df[i, BEASTlabel]), parent=seq, doc=bxml, addFinalizer=T)
				tmp		<- which( rownames(seq.PROT.RT)==df[i, FASTASampleCode] )
				if(length(tmp)!=1)	stop("unexpected row in seq.PROT.RT selected")
				tmp		<- paste(as.character(seq.PROT.RT[tmp,])[1,],collapse='',sep='')						
				dummy	<- newXMLTextNode(text=tmp, parent=seq, doc=bxml,addFinalizer=T)
				seq
			})	
	if(verbose)	cat(paste("\nadded new alignments, n=", xmlSize(seqalign)))	
	bxml
}	
######################################################################################
#	Read a BEAST treeannotator file. The node.label is set to a data.table that contains the SIMMAP annotation for the interior nodes in the newick tree.
hivc.treeannotator.read<- function(file, add.to.tiplabel=NA, rate.multiplier=NA, round.digit=NA, verbose=1) 
{
	require(data.table)
	require(ape)
	ph 			<- read.nexus(file)
	if(verbose)	cat(paste("\nReading BEAST file",file,sep=''))
	X 			<- scan(file = file, what = "", sep = "\n", quiet = TRUE)	
	#	read annotations of node in order as they appear in 'file'
	tab			<- X[grep("tree TREE1[[:space:]]+=", X)]
	tab 		<- gsub("tree TREE1[[:space:]]+= \\[&R\\] ", "", tab)
	tab 		<- unlist(strsplit(tab, "\\["))[-1]
	tab 		<- gsub("&|;|\\]", "", tab)
	tab 		<- gsub(":.+$", "", tab)
	tab 		<- lapply(tab, function(x) unlist(strsplit(x, ","))	)	
	tab			<- lapply(tab, function(x)
			{
				x			<- gsub('%','',x,fixed=1)
				ind 		<- grep("[{]", x)
				names 		<- gsub("=.+$", "", x[ind])
				x[ind] 		<- gsub("[{]", "", x[ind])
				x[ind] 		<- gsub("=", "_MIN=", x[ind])
				x[ind + 1] 	<- gsub("[}]", "", x[ind + 1])
				x[ind + 1] 	<- paste(paste(names, "MAX=", sep = "_"), x[ind + 1],sep='')
				x
			})
	colnames	<- unique(gsub("=.+$", "", unlist(tab)))
	if(verbose)	cat(paste("\nFound BEAST variables ",paste(colnames,collapse=', '),sep=''))
	tab			<- c(list(paste(colnames,-1,sep='=')), tab)		#rbindlist bug fix
	tab			<- lapply(tab, function(x)
			{
				ans									<- rep(NA, length(colnames))
				names(ans)							<- colnames				
				x									<- strsplit(x,'=',fixed=1)
				ans[ sapply(x, function(z) z[1]) ]	<- sapply(x, function(z) z[2])
				ans									<- paste(apply( rbind( names(ans), ans ), 2, function(z) paste(z,collapse='=',sep='')),collapse=',')
				eval(parse(text=paste("data.table(",ans,")",sep='')))				
			})
	df.beast	<- rbindlist(tab)[-1,]
	if(!is.na(rate.multiplier))
	{
		tmp	<- colnames(df.beast)[grepl("rate",colnames(df.beast))]
		sapply(tmp, function(x)		set(df.beast, NULL, x, as.numeric(unlist(df.beast[,x,with=F]))*rate.multiplier)			)		
	}
	if(!any(is.na(round.digit)))
	{
		if(length(round.digit)!=ncol(df.beast))
			round.digit<- rep(round.digit[1], ncol(df.beast))
		tmp	<- colnames(df.beast)
		sapply(seq_along(tmp), function(i)  
				{
					if(class(df.beast[[i]])=="numeric")			
						set(df.beast, NULL, tmp[i], round(as.numeric(unlist(df.beast[,tmp[i],with=F])), d=round.digit[i]))				
				})
	}
	
	tmp			<- length(which(df.beast[,!is.na(posterior)]))
	if(verbose)	cat(paste("\nFound annotated nodes, n=", tmp))
	if(verbose)	cat(paste("\nFound annotated tips, n=", nrow(df.beast)-tmp))	
	#	determine node index for 'df.beast':
	#
	#	- delete SIMMAP information from 'X'	
	LEFT 		<- grep("\\[", X)
	RIGHT 		<- grep("\\]", X)
	if (length(LEFT)) 
	{
		w 	<- LEFT == RIGHT
		if (any(w)) 
		{
			s 		<- LEFT[w]
			X[s] 	<- gsub("\\[[^]]*\\]", "", X[s])
		}
		w <- !w
		if(any(w)) 
		{
			s 		<- LEFT[w]
			X[s] 	<- gsub("\\[.*", "", X[s])
			sb	 	<- RIGHT[w]
			X[sb] 	<- gsub(".*\\]", "", X[sb])
			if(any(s < sb - 1)) 
				X 	<- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
		}
	}	
	#	- read tree block
	endblock 			<- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	semico 				<- grep(";", X)
	i1 					<- grep("BEGIN TREES;", X, ignore.case = TRUE)
	i2 					<- grep("TRANSLATE", X, ignore.case = TRUE)
	tree 				<- X[(semico[semico > i2][1] + 1):(endblock[endblock > i1][1] - 1)]
	tree 				<- gsub("^.*= *", "", tree)
	tree				<- substr(tree, 1, nchar(tree)-1)	 
	# 	- For each node, add a dummy node label that is the index in 'df.beast'	
	tmp					<- unlist(strsplit(tree, ":"))
	interiorm1			<- which( df.beast[,!is.na(posterior)] )
	tmp					<- sapply(seq_along(tmp), function(i) 			ifelse(i %in% interiorm1, 	paste(tmp[i],i,sep=''),	tmp[i])				)
	tmp					<- paste(paste(tmp, collapse=':'),';',sep='')
	#	- read this newick string and determine the node index in 'df.beast'
	tmp					<- read.tree(text=tmp)
	ph[["node.label"]]	<- cbind(data.table(node=Ntip(ph) + seq_len(Nnode(ph))), df.beast[as.numeric( tmp$node.label ),])
	setkey(ph[["node.label"]], node)
	
	if(!any(is.na(add.to.tiplabel)))
	{
		if(length(intersect(add.to.tiplabel,colnames(df.beast)))!=length(add.to.tiplabel))
			stop("Cannot find add.to.tiplabel")
		tmp				<- as.matrix(subset( df.beast[as.numeric( tmp$tip.label ),], select=add.to.tiplabel, with=F))
		tmp				<- cbind(ph$tip.label, tmp)
		ph$tip.label	<- apply(tmp, 1, function(x) paste(x,collapse='_'))		
	}	
	
	ph
}
######################################################################################
#	write nexus file for all sequences specified in df. assumes df has BEASTlabel. assumes seq.DNAbin.matrix and ph contain FASTASampleCode in df.
hivc.beast.writeNexus4Beauti<- function( seq.DNAbin.matrix, df, ph=NULL, file=NULL )
{
	# 	select sequences and relabel					
	seq.DNAbin.matrix			<- seq.DNAbin.matrix[ df[,FASTASampleCode], ]	
	rownames(seq.DNAbin.matrix)	<- df[,BEASTlabel]
	#	generate starting tree and relabel
	if(!is.null(ph))
	{			
		tmp					<- setdiff( ph$tip.label, df[,FASTASampleCode] )
		tmp					<- match( tmp, ph$tip.label)
		ph.start			<- drop.tip(ph, tmp)
		ph.start$node.label	<- NULL
		setkey(df, FASTASampleCode)
		ph.start$tip.label	<- df[ph.start$tip.label,][,BEASTlabel]		
	}
	#	produce nexus text			
	ans<- seq.write.dna.nexus(seq.DNAbin.matrix, ph=ph.start, file=file)
	ans
}
######################################################################################
#	read tip stem samples after burn in from log file and return upper left points of a histogram of the Monte Carlo sample
hivc.beast.read.log2tstem<- function(file.log, file.xml, beastlabel.idx.samplecode=6, burn.in= 5e6, breaks.n= 30, verbose=0)
{
	require(data.table)
	require(XML)	
	if(verbose)	cat(paste("\nReading file ",file.log))
	df.log	<- read.delim2(file.log, header = TRUE, sep = "\t", quote="\"", dec=".", fill = TRUE, comment.char="#")
	df.log	<- as.data.table(df.log)		
	if(verbose)	cat(paste("\nRead log file with ncol=", ncol(df.log)))
	tmp		<- c("state", colnames(df.log)[ grepl("tstem", colnames(df.log)) ] )
	if(length(tmp)==1)
	{
		if(verbose)	cat(paste("\nNo tstem found, return NA"))
		return( data.table( FASTASampleCode=NA, tstem= NA, density= NA ) )
	}
	df.log	<- subset( df.log, state>burn.in, select=tmp )
	if(verbose)	cat(paste("\nFound tstem data for n=", ncol(df.log)-1))
	#	translate tips back into FASTASampleCode
	if(verbose)	cat(paste("\nReading file ",file.xml))
	bxml				<- xmlTreeParse(file.xml, useInternalNodes=TRUE, addFinalizer = TRUE)
	tmp					<- regexpr("tip[0-9]+",colnames(df.log))
	if(!length(which(tmp>0)))
	{
		if(verbose)	cat(paste("\nNo tip tstem found, return NA"))
		return( data.table( FASTASampleCode=NA, tstem= NA, density= NA ) )
	}	
	log.tips			<- sapply(seq_along(tmp)[tmp>0], function(i)	substr(colnames(df.log)[i],tmp[i],tmp[i]+attr(tmp,"match.length")[i]-1)		)
	log.FASTASampleCode	<- sapply(log.tips, function(x) unlist( xpathApply(bxml, paste("//taxa[@id='",x,"']/taxon",sep=''), xmlGetAttr, "idref") )	)
	log.FASTASampleCode	<- sapply( strsplit(log.FASTASampleCode, '_'), function(x) x[beastlabel.idx.samplecode])		
	setnames(df.log, colnames(df.log)[tmp>0], log.FASTASampleCode)
	#	compute histograms for each tip stem
	ans		<- lapply(log.FASTASampleCode, function(x)
			{
				#x<- "R11-11357"
				tmp	<- hist( as.numeric(unlist( df.log[, x, with=0] )), breaks=breaks.n, plot=0 )
				data.table( FASTASampleCode=x, tstem= tmp$breaks, density= c(tmp$density,0) )										
			})
	if(verbose)	cat(paste("\nProcessed tstem data for tips n=", length(log.FASTASampleCode)))
	ans	<- rbindlist(ans)
}
######################################################################################
hivc.beast2.add.data<- function(bxml, seq.PROT.RT, df, beast2.spec, verbose=1)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	dummy		<- newXMLCommentNode(text="The sequences.", parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(df)," nchar=",ncol(seq.PROT.RT),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	seqalign	<- newXMLNode("data", attrs= list(id=beast2.spec$data.id, dataType=beast2.spec$data.dataType, missing=beast2.spec$data.missing), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- lapply( seq_len(nrow(df)), function(i)
			{
				tmp		<- which( rownames(seq.PROT.RT)==df[i, FASTASampleCode] )
				if(length(tmp)!=1)	stop("unexpected row in seq.PROT.RT selected")
				tmp		<- paste(as.character(seq.PROT.RT[tmp,])[1,],collapse='',sep='')										
				seq		<- newXMLNode("sequence", attrs=list(id= paste('seq',df[i, BEASTlabel],sep=''), taxon=df[i, BEASTlabel], value=tmp), parent=seqalign, doc=bxml, addFinalizer=T)
				seq
			})	
	if(verbose)	cat(paste("\nadded new sequences, n=", xmlSize(seqalign)))
	bxml
}
######################################################################################
#	df<- cluphy.df
hivc.beast2.poolclusters.mincnts<- function(df, beast2.spec, verbose=1)
{
	if(verbose) cat(paste("\npool evenly across clu.AnyPos_T1 and fill sequences so that min cnts requested per time period is met"))
	pool.ntip		<- beast2.spec$pool.ntip
	cnts.requested	<- beast2.spec$pool.cnts.requested
	breaks			<- c(Inf, beast2.spec$bdsky.sprop.changepoint.value)
	#	add tip heights and tip periods to cluphy.df
	tmp				<- df[, max(PosSeqT)]
	tmp				<- as.numeric( df[, difftime(tmp, PosSeqT, units='days') / 365] )	
	df				<- merge( df, data.table(TipHeight=tmp, TipPeriod=as.character(cut(tmp, breaks, right=FALSE)), FASTASampleCode=df[,FASTASampleCode]), by='FASTASampleCode')	
	#	first attempt of pooling
	clu.df			<- df[, list(clu.ntip=clu.ntip[1], clu.AnyPos_T1=clu.AnyPos_T1[1]), by="cluster"]
	df.mxclu		<- clu.df[, { tmp<- c(which.min(clu.AnyPos_T1),which.max(clu.AnyPos_T1)); list(cluster= cluster[tmp])}]		
	setkey(df, clu.AnyPos_T1)
	pool.n			<- ceiling( sum( clu.df[,clu.ntip] ) / pool.ntip )
	tmp				<- lapply( seq_len(pool.n), function(x)	seq.int(x,nrow(clu.df),by=pool.n) )
	#	select the clusters evenly across clu.AnyPos_T1 and always add the earliest and latest cluster
	pool.df			<- lapply(seq_along(tmp), function(i) unique(rbind(subset(clu.df[tmp[[i]],], select=cluster), df.mxclu)) )		
	pool.df			<- lapply(seq_along(tmp), function(i) merge(pool.df[[i]], df, by="cluster") )
	#	compute the intial numbers of selected sequences per time period
	cnts			<- sapply(seq_along(pool.df), function(i)	table( pool.df[[i]][,TipPeriod])		)
	if(verbose)		print(cnts)
	if( any( na.omit(rowSums(cnts)<cnts.requested) ) )	
		stop('total number of sequences is smaller than cnts.requested')
	#	compute the number of sequences to be added
	setkey(df, FASTASampleCode)	
	cnts							<- cnts.requested - cnts
	cnts[is.na(cnts)|cnts<=0]		<- 0
	#	add sequences to pool	
	for(i in seq_along(pool.df))
		for(j in seq_len(nrow(cnts)))
			if(cnts[j,i]>0)
			{
				cat(paste('\nadding sequences to pool',i,'for period',rownames(cnts)[j]))
				#	only if cnts>0
				#	determine clusters that can be added
				pool.df.notin			<- setdiff( subset(df,TipPeriod==rownames(cnts)[j])[,FASTASampleCode], pool.df[[i]][,FASTASampleCode] )
				clu.df.notin			<- unique( subset( df[pool.df.notin,], select=cluster) )
				clu.df.notin			<- merge(df, clu.df.notin, by='cluster')
				clu.df.notin			<- clu.df.notin[, list(clu.ntipperiod=length(which(TipPeriod==rownames(cnts)[j]))), by=cluster]	
				#	determine clusters that will be added
				tmp						<- clu.df.notin[,tail(which(cumsum(clu.ntipperiod)<cnts[j,i]),1),]
				tmp						<- ifelse(length(tmp), tmp[1]+1, 1)
				clu.df.notin			<- clu.df.notin[ seq_len( min( nrow(clu.df.notin), tmp ) ), ]
				#	add clusters
				pool.df[[i]]			<- rbind( pool.df[[i]], merge( df, subset(clu.df.notin, select=cluster), by='cluster') )
				tmp						<- cnts.requested - as.numeric( table( pool.df[[i]][,TipPeriod]) )
				tmp[is.na(tmp)|tmp<0]	<- 0 
				cnts[,i]				<- tmp
				cat(paste('\nnew number of sequences in pool',i,'is',nrow(pool.df[[i]])))
			}
	#	compute the final numbers of selected sequences per time period
	cnts	<- sapply(seq_along(pool.df), function(i)	table( pool.df[[i]][,TipPeriod])		)
	if(verbose)		print(cnts)
	list(pool.df=pool.df, pool.ntip=pool.ntip)
}
######################################################################################
hivc.beast2.extract.distinct.topologies<- function(mph.clu)					
{
	#compute if topology between retained mph.clu's is identical (branch lengths may differ)
	tmp				<- sapply( seq_along(mph.clu)[-1], function(mph.i) all.equal(mph.clu[[mph.i-1]], mph.clu[[mph.i]], use.edge.length=FALSE, use.tip.label=TRUE, index.return=FALSE ) )
	tmp				<- data.table( mph.i= seq_along(mph.clu)[-1], equal.to.previous= tmp )				
	#	dtopo -> distinct topologies
	mph.clu.dtopo	<- data.table( mph.i= c(1, subset(tmp, !equal.to.previous)[, mph.i]), freq= as.vector(table(cumsum(!c(TRUE,tmp[,equal.to.previous])))), collapsed=FALSE	)					
	#	itopo -> identical topologies 
	mph.clu.itopo	<- cbind(tmp, equal.to= 1+cumsum(!tmp[,equal.to.previous]))
	mph.clu.itopo	<- rbind( data.table(mph.i=1, equal.to=1), subset(mph.clu.itopo, select=c(mph.i, equal.to)) )
	#	select by index
	set(mph.clu.itopo, NULL, 'equal.to', mph.clu.dtopo[mph.clu.itopo[,equal.to],][,mph.i])
	#	now select by key				
	setkey(mph.clu.dtopo, mph.i)
	#for each sequentially different topology, check if subsequent ones are identical and if yes collapse					
	while(mph.clu.dtopo[,any(!collapsed)])
	{
		cat(paste('\nprogress: number of seq distinct topologies =',nrow(mph.clu.dtopo)))
		mph.indexrow	<- mph.clu.dtopo[, which(!collapsed)][1]
		if(mph.indexrow<nrow(mph.clu.dtopo))
		{
			mph.index		<- mph.clu.dtopo[mph.indexrow, mph.i]
			mph.is			<- mph.clu.dtopo[seq.int(mph.indexrow+1,nrow(mph.clu.dtopo)), ][,mph.i]						
			tmp				<- sapply( mph.is, function(mph.i) all.equal(mph.clu[[mph.index]], mph.clu[[mph.i]], use.edge.length=FALSE, use.tip.label=TRUE, index.return=FALSE ) )
			topo.duplicates	<- data.table( mph.index= mph.index, mph.i=mph.is, equal.to.index= tmp )
			topo.duplicates	<- merge( mph.clu.dtopo, subset(topo.duplicates, equal.to.index), by='mph.i' )
			if(nrow(topo.duplicates))
			{
				tmp				<- topo.duplicates[, list(index=seq_len(freq)), by='mph.i']
				set(mph.clu.itopo, as.integer( tmp[, mph.i+index-1] ), 'equal.to', mph.index)
				set(mph.clu.dtopo, mph.indexrow, 'freq', mph.clu.dtopo[mph.indexrow, freq] + topo.duplicates[, sum(freq)])							
				collapsed.i		<- setdiff(mph.clu.dtopo[,mph.i], topo.duplicates[,mph.i])
				mph.clu.dtopo	<- mph.clu.dtopo[J(collapsed.i),]	
			}							
		}
		set(mph.clu.dtopo, mph.indexrow, 'collapsed', TRUE)
	}
	if( mph.clu.dtopo[, sum(freq)]!=length(mph.clu) )	stop('unexpected freq in mph.clu.dtopo')
	if( !setequal(mph.clu.dtopo[, mph.i], mph.clu.itopo[, unique(equal.to)]) )  stop('unexpected mph.i difference between mph.clu.dtopo and mph.clu.itopo')
	cat(paste('\nnumber of distinct topologies in mph.clu.dtopo=',nrow(mph.clu.dtopo)))
	cat(paste('\nnumber of distinct topologies in mph.clu.itopo=',length(unique(mph.clu.itopo[, equal.to]))))
	
	mph.clu.dtopo	<- subset(mph.clu.dtopo, select=c(mph.i, freq))	
	mph.clu.dtopo[, dens:= round(mph.clu.dtopo[, freq]/length(mph.clu),d=3) ]	
	mph.clu.dtopo[, topo.n:= nrow(mph.clu.dtopo) ]
	mph.clu.dtopo[, tip.n:= Ntip(mph.clu[[1]]) ]			
	#	creates warnings because list(DT) creates a copy of the data.table
	list(dtopo= mph.clu.dtopo, itopo=mph.clu.itopo)
}
######################################################################################
hivc.beast2out.combine.clu.trees<- function(indir, file.info, beastlabel.idx.clu=1, beastlabel.idx.hivn=2, beastlabel.idx.hivd=3, beastlabel.idx.hivs=4, beastlabel.idx.samplecode= 6, beastlabel.idx.rate= NA)
{
	#	collect consensus tree and further info for plotting
	tmp			<- lapply(seq_len(nrow(file.info)), function(i)
			{				
				#	load dated cluster phylogenies
				file				<- paste(indir, file.info[i,file], sep='/')
				cat(paste('\nload file=',file,'i=',i))
				tmp					<- load(file)
				topo.map			<- mph.clu.dtopo[which.max(freq),]					
				tmp					<- which( grepl( paste('mph.i=',topo.map[,mph.i],'_',sep=''), names(ph.consensus) ) )
				if(length(tmp)!=1)	stop('unexpected ph.consensus index')
				topo.map.ph			<- ph.consensus[[ tmp ]]
				topo.map.SA			<- subset(mph.SA.cnt, equal.to==topo.map[,mph.i])
				set(topo.map.SA, NULL, 'SA.freq', topo.map.SA[,SA.freq/n])													
				topo.map.nodectime	<- subset(mph.node.ctime, equal.to==topo.map[,mph.i])
				topo.map			<- merge(subset(topo.map, select=c(cluster, dens)), subset(topo.map.SA, select=c(cluster, tip, SA.freq)), by='cluster')
				topo.map.nodectime	<- subset(topo.map.nodectime, select=c(cluster, node, q, cdf, pdf))					
				list(map= topo.map, nodectime=topo.map.nodectime, ph=topo.map.ph)			
			})
	cluphy.map					<- do.call('rbind',lapply(tmp, function(x) x$map))
	cluphy.map.nodectime		<- do.call('rbind',lapply(tmp, function(x) x$nodectime))
	cluphy.subtrees				<- lapply(tmp, function(x) x$ph)
	names(cluphy.subtrees)		<- file.info[,cluster]
	#	determine subtree root times and set the root time of the combined phylogeny to the mean
	cluphy.subtrees.root.ctime	<- sapply(seq_along(cluphy.subtrees), function(i)
			{
				tmp		<- hivc.treeannotator.tiplabel2df(cluphy.subtrees[[i]], beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
				tmp[1,TipT]-node.depth.edgelength(cluphy.subtrees[[i]])[1]-cluphy.subtrees[[i]]$root.edge
			})
	cluphy.root.ctime			<- mean(cluphy.subtrees.root.ctime)
	for(i in seq_along(cluphy.subtrees))
		cluphy.subtrees[[i]]$root.edge	<- cluphy.subtrees[[i]]$root.edge+cluphy.subtrees.root.ctime[i]-cluphy.root.ctime	
	brl.error	<- max( sapply(seq_along(cluphy.subtrees), function(i)
			{
				tmp		<- hivc.treeannotator.tiplabel2df(cluphy.subtrees[[i]], beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
				tmp		<- node.depth.edgelength(cluphy.subtrees[[i]])[seq_len(Ntip(cluphy.subtrees[[i]]))]+cluphy.subtrees[[i]]$root.edge+cluphy.root.ctime-tmp[, TipT]
				max(abs(tmp))
			}) )
	if(brl.error>2*EPS)	warning(paste('found brl error, max error=',brl.error))
	#	combine subtree phylogenies
	cluphy						<- hivc.clu.polyphyletic.clusters(cluphy.subtrees=cluphy.subtrees)$cluphy
	cluphy.info					<- hivc.treeannotator.tiplabel2df(cluphy, beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)		
	cluphy$tip.label			<- cluphy.info[, FASTASampleCode]
	cluphy.tip.ctime			<- cluphy.info[, TipT]
	if(cluphy.root.ctime!=cluphy.info[1,TipT]-node.depth.edgelength(cluphy)[1])	stop('unexpected root time of the combined phylogeny')
	#	need to map nodes from subtrees to nodes in combined phylogeny 
	cluphy.info					<- merge(cluphy.info, cluphy.info[,	list(mrca= hivc.clu.mrca(cluphy, FASTASampleCode)$mrca),	by='cluster'], by='cluster')		
	cluphy.vertexmap			<- cluphy.info[,	{
				tmp		<- cluphy.subtrees[[as.character(cluster)]]
				list( 	ph.vertex=c(mrca[1], Descendants(cluphy, mrca[1], type="all")), clu.vertex=c(Ntip(tmp)+1, Descendants(tmp, Ntip(tmp)+1, type="all")) )
			},by='cluster']
	#	cluphy.vertexmap does not contain time of root to mrca (node 0) and we don t need it.
	#	cluphy.vertexmap contains tip indices and we don t need those.
	setnames(cluphy.vertexmap,'clu.vertex','node')
	cluphy.map.nodectime		<- merge(cluphy.map.nodectime, cluphy.vertexmap, by=c('cluster','node'))
	set(cluphy.map.nodectime, NULL, 'node', cluphy.map.nodectime[,ph.vertex]) 
	cluphy.map.nodectime[, ph.vertex:=NULL]				
	#	set node.label to MAP prob
	tmp											<- merge( unique(subset(cluphy.info, select=c(cluster, mrca))), unique(subset(cluphy.map, select=c(cluster, dens))), by='cluster' )
	cluphy$node.label							<- rep('',Nnode(cluphy))
	cluphy$node.label[tmp[,mrca-Ntip(cluphy)]]	<- as.character(tmp[,dens])		
	list(cluphy=cluphy, cluphy.subtrees=cluphy.subtrees, cluphy.info=cluphy.info, cluphy.map.nodectime=cluphy.map.nodectime)
}
######################################################################################
hivc.beast2out.plot.cluster.trees<- function(df.all, df.immu, df.viro, df.treatment, ph, ph.root.ctime, ph.tip.ctime, ph.prob=NA, df.node.ctime=NULL, df.rates=NULL, df.tips=NULL, end.ctime=2013.3,  cex.nodelabel=0.5,  cex.tiplabel=0.5,  file=NULL,  pdf.width=7, pdf.height=20, pdf.xlim=NULL)
{
	#df.all, df.immu, df.viro, df.treatment, 
	#df.tips	<- df.tpairs.plot; ph<- cluphy; ph.root.ctime=cluphy.root.ctime; ph.tip.ctime=cluphy.tip.ctime; df.node.ctime=cluphy.map.nodectime;  cex.nodelabel=0.5;  cex.tiplabel=0.5;  file=NULL;  pdf.width=7; pdf.height=120; pdf.xlim=pdf.xlim
	require(RColorBrewer)
	if(class(file)=="character")
		pdf(file, width=pdf.width, height=pdf.height)
	par(mar=rep(0,4))		
	youngest.tip.ctime	<- max(ph.tip.ctime)
	cols				<- brewer.pal(12,"Paired")
	cols[1]				<- cols[8]			
	#	get tip labels	
	ph.tiplabel			<- hivc.clu.get.tiplabels(ph, 	df.all, col.notmsm="#4EB3D3", col.Early="#EF9708", col.highVL="#FEE391", col.AfterTreat="#D4B9DA", col.green="#D9F0A3", col.latePres="#FA9FB5", select=c("CountryInfection","Trm","Sex","isAcute","lRNA.early","cluster","Patient","RegionHospital") )
	if(is.null(pdf.xlim))
	{
		tmp				<- max( apply(ph.tiplabel$text, 2, function(x)  sum(xinch(strwidth(x, units="inches", cex=cex.tiplabel)))  ) )
		pdf.xlim		<- c(-3, ceiling(end.ctime)-ph.root.ctime+max(2.5,tmp))		
	}
	pdf.ylim			<- c(1,Ntip(ph)) + c(-1,1)							
	plot(ph, x.lim=pdf.xlim, y.lim= pdf.ylim, show.tip.label=0, edge.color = 1, tip.color = 0)
	# add calendar timeline
	hivc.treeannotator.plot.ctimeline(ph, youngest.tip.ctime, end.ctime, add.yinch= 0.5)
	# add NegT and AnyPos_T1
	ph.seronodeheight	<- hivc.treeannotator.sero.getnodeheight.range(ph, df.all, youngest.tip.ctime)
	hivc.treeannotator.plot.seronodeheightrange(ph, ph.seronodeheight, add.yinch= -0.03, width.yinch= 0.03, width.yinch.past.AnyPos_T1= 0, col=cols[2])
	# add lRNA timeline
	ph.viro.timeline	<- hivc.treeannotator.get.viro.timeline(ph, df.all, df.viro, youngest.tip.ctime, df.treatment=df.treatment)
	hivc.treeannotator.plot.viro.timeline(ph, ph.viro.timeline, viro.min= log10(300), width.yinch= 0.15, add.yinch= 0.005, col.bg= cols[c(5,10,12)], col.legend= cols[6], cex.txt= 0.2, lines.lwd=0.1)
	# add BEAST posterior density of nodes where available
	if(!is.null(df.node.ctime))
		hivc.treeannotator.plot.node.ctime(df.node.ctime, ph.root.ctime, width.yinch=0.1, add.yinch=0.005, col.bg=cols[1] )
	# add CD4 timeline
	ph.immu.timeline	<- hivc.treeannotator.get.immu.timeline(ph, df.all, df.immu, youngest.tip.ctime, end.ctime=2013.3)
	hivc.treeannotator.plot.immu.timeline(ph, ph.immu.timeline, immu.min= 150, immu.max= 800, width.yinch= 0.15, add.yinch= -0.005, col.bg= cols[3], col.legend= cols[4], cex.txt= 0.2, lines.lwd=0.1)
	# re-plot phylogeny
	if(!is.null(ph$node.label))
		ph$node.label	<- as.numeric(sapply( strsplit( ph$node.label, '_' ), function(x)	x[1] ))
	edge.width			<- hivc.treeannotator.get.edgewidth(ph, df.rates, scale.edgewidth= 8)
	hivc.phy.plotupon(ph, show.tip.label=0, show.node.label=ifelse(is.null(ph$node.label),0,1), cex=cex.nodelabel, edge.width=edge.width[,width])
	# add root edge
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	lines(c(-ph$root.edge,0),rep(lastPP$yy[Ntip(ph)+1],2), lwd=edge.width[1,width])
	# plot tick marks
	if(!is.null(df.tips))
		hivc.treeannotator.plot.tipmarks(ph, df.tips, add.xinch=0, add.yinch=0)			
	# add rate labels			
	if(!is.null(df.rates))	
		hivc.treeannotator.plot.rates(ph, edge.width, add.xinch=-0.1, cex.rate=0.3)
	# add tip labels								
	tmp					<- rep( max(node.depth.edgelength(ph)) - (youngest.tip.ctime-ceiling(end.ctime)), Ntip(ph))
	hivc.clu.plot.tiplabels(seq_len(Ntip(ph)), ph.tiplabel$text, ph.tiplabel$col, xx=tmp, adj = c(-0.05, 0.5), cex=cex.tiplabel, add.xinch= 0.03, add.yinch= 0.02)
	# add legend	
	legend("topright", fill= cols[c(1,2,3,5,10,12)], legend=c("SA-BEAST2 TMRCA pdf", "interval [last HIV-, diagnosis]", "CD4 timeline", "VL timeline", "VL timeline under treatment", "VL timeline under treatment"), bty='n', border=NA, cex=cex.tiplabel)
	if(!is.na(ph.prob))
		legend("bottomright", legend=paste('prob=',ph.prob),bty='n', border=NA, cex=cex.tiplabel*2)
	if(class(file)=="character")
		dev.off()	
}
######################################################################################
hivc.beast2out.pool.cluster.trees<- function(files)
{
	tmp				<- regmatches( files, regexpr('_pool_[0-9]+',files)) 
	run				<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
	mph.clu.pool	<- NULL
	for(i in seq_along(run))
	{
		options(show.error.messages = FALSE)		
		readAttempt		<- try(suppressWarnings(load(files[i])))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",files[i]))
		if(inherits(readAttempt, "try-error"))	stop(paste("\ncould not read file",files[i]))		
		names(mph.clu)	<- paste('RUN_',run[i],'_',names(mph.clu),sep='')
		mph.clu.pool	<- c(mph.clu.pool, mph.clu)			
	}
	mph.clu.pool			
}
######################################################################################
hivc.beast2out.extract.cluster.trees<- function(mph, mph.info)
{
	mph.clusters		<- unique(mph.info[,cluster])		
	mph.by.clu			<- lapply(mph.clusters, function(clu)
			{
				cat(paste('\nprocess cluster',clu))
				mph.clu			<- lapply( seq_along(mph), function(mph.i)
						{
							clu.seq			<- subset(mph.info, cluster==clu)[,BEASTlabel]					
							clu.mrca		<- hivc.clu.mrca(mph[[mph.i]], clu.seq)
							extract.clade( mph[[mph.i]], clu.mrca$mrca, root.edge=clu.mrca$mrca.height, interactive=FALSE )					
						})	
				names(mph.clu)	<- names(mph)
				mph.clu
			})
	names(mph.by.clu)	<- mph.clusters
	mph.by.clu
}
######################################################################################
hivc.beast2out.tip.date.check<- function(ph, fun, ...)
{	
	ph.info	<- fun(ph, ...) 
	ph.info[, tip:= match(ph.info[,BEASTlabel], ph$tip.label)]
	setkey(ph.info, tip)
	tmp		<- node.depth.edgelength(ph)[seq_len(Ntip(ph))] 
	max( abs( ph.info[1,TipT]-tmp[1]+tmp  -  ph.info[,TipT] ) ) 		
}
######################################################################################
hivc.beast2out.read.trees<- function(file, opt.rescale.edge.length= 1., opt.burnin=0 )
{
	mph			<- read.nexus(file)
	#	remove burn in 
	tmp			<- regexpr('[0-9]+',names(mph))
	if(any(tmp<0))	stop('unexpected nexus file without STATE iteration numbers')
	mph.it		<- as.numeric( regmatches( names(mph), tmp) )
	mph			<- lapply( which( mph.it>opt.burnin), function(j)	mph[[j]]	)
	mph.it		<- mph.it[ mph.it > opt.burnin ] 
	names(mph)	<- paste('STATE_',mph.it,sep='')
	#	rescale edge lengths
	if(opt.rescale.edge.length!=1.)	
		for(j in seq_along(mph)) 
			mph[[j]]$edge.length	<- mph[[j]]$edge.length * opt.rescale.edge.length
	mph
}
######################################################################################
hivc.beast2out.get.log.evolrate<- function(files, opt.burnin=opt.burnin)
{		
	if(!opt.burnin)	warning('burn in equals zero')
	#	collect posterior TreeHeight samples
	df.log		<- lapply(seq_along(files), function(i)
			{
				file		<- files[i]
				cat(paste("\nReading file ",file))
				df.log		<- read.delim2(file, header = TRUE, sep = "\t", quote="\"", dec=".", fill = TRUE, comment.char="#")
				df.log		<- as.data.table(df.log)
				is.BEAST2	<- ifelse('Sample'%in%colnames(df.log), 1, 0)
				cat(paste("\ndetected is.BEAST2=",is.BEAST2))
				if(is.BEAST2)
				{					
					if(!'Sample'%in%colnames(df.log)) stop('expect column Sample in log file')
					if(!'ucldMean'%in%colnames(df.log)) stop('expect column ucldMean in log file')					
					df.log	<- subset(df.log, Sample>opt.burnin, select=c(Sample, ucldMean))					
				}
				else
				{
					if(!'state'%in%colnames(df.log)) stop('expect column state in log file')
					if(!'ucld.mean'%in%colnames(df.log)) stop('expect column ucld.mean in log file')					
					df.log	<- subset(df.log, state>opt.burnin, select=c(state, ucld.mean))
					setnames(df.log, 'ucld.mean','ucldMean')
				}
				df.log[, file.i:= i]
				df.log
			})
	df.log		<- do.call('rbind', df.log)
	#	compute overall mean and file specific mean
	ans			<- df.log[, list(ucldMean.mean.i=mean(ucldMean)), by='file.i']
	ans[,ucldMean.mean:= mean(ucldMean.mean.i)]
	ans
}		
######################################################################################
hivc.beast2.add.alignment<- function(bxml, beast2.spec, verbose=1)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	dummy		<- newXMLCommentNode(text="The sequence alignments.", parent=bxml.beast, doc=bxml, addFinalizer=T)
	if(length(beast2.spec$alignment.filter)==1 && is.na(beast2.spec$alignment.filter))
	{		
		dummy	<- newXMLCommentNode(text="The <alignment> is the <data> block if there is no alignment filter.", parent=bxml.beast, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("taxa", attrs= list(id=paste("TaxonSet.t",beast2.spec$alignment.id,sep=':'), spec=beast2.spec$alignment.taxa.spec, alignment=paste('@',beast2.spec$alignment.id,sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)
	}
	else	
		dummy	<- lapply(seq_along(beast2.spec$alignment.filter), function(i)
				{
					dummy	<- newXMLNode("alignment", attrs= list(id=beast2.spec$alignment.id[i], filter=beast2.spec$alignment.filter[i], data=paste('@',beast2.spec$data.id,sep=''), spec=beast2.spec$alignment.spec), parent=bxml.beast, doc=bxml, addFinalizer=T)
					dummy	<- newXMLNode("taxa", attrs= list(id=paste("TaxonSet.t",beast2.spec$alignment.id[i],sep=':'), spec=beast2.spec$alignment.taxa.spec, alignment=paste('@',beast2.spec$alignment.id[i],sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)
				})		
	if(verbose)	cat(paste("\nadded new alignments and taxonsets, n=", length(beast2.spec$alignment.id)))	
	bxml
}
######################################################################################
hivc.beast2.add.datetrait<- function(bxml, df, beast2.spec, verbose=1)	
{			
	tmp			<- df[,BEASTlabel]	
	if(verbose)	cat(paste("\nsetting tip date to time at pos x in label, x=", beast2.spec$beast.label.datepos))
	tmp			<- sapply( strsplit(tmp, beast2.spec$beast.label.sep, fixed=1), function(x) paste(paste(x,collapse=beast2.spec$beast.label.sep,sep=''),'=',x[beast2.spec$beast.label.datepos],sep='') )
	
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	dummy		<- newXMLCommentNode(text="The tip dates.", parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(df),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	dummy		<- newXMLNode("trait", attrs= list(id=beast2.spec$datetrait.id, 
					spec= beast2.spec$datetrait.spec,
					units= beast2.spec$datetrait.units,
					taxa= paste('@',beast2.spec$datetrait.taxa,sep=''),
					traitname= beast2.spec$datetrait.traitname, 
					value=paste(tmp, collapse=", ",sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadded new date trait for sequences, n=", length(tmp)))
	#define taxonset in here
	bxml
}
######################################################################################
hivc.beast2.add.satree<- function(bxml, beast2.spec, verbose=verbose)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	if(1 || is.na(beast2.spec$starttree.newick))
	{
		if(verbose)	cat(paste('\nadd SA tree, initialized with UPGMA clustering'))
		dummy		<- newXMLNode("tree", attrs= list(	id=beast2.spec$tree.id, 													 
						trait=paste('@',beast2.spec$datetrait.id,sep=''),
						nodetype=beast2.spec$sasky.tree.nodetype,
						clusterType=beast2.spec$sasky.tree.cluster.type,
						spec=beast2.spec$sasky.tree.cluster.spec,
						taxa=paste('@',beast2.spec$tree.taxonset,sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)			
	}
	else
	{
		if(verbose)	cat(paste('\nadd SA tree, initialized with newick tree'))
		bxml.tree	<- newXMLNode("tree", attrs= list(	id=beast2.spec$tree.id, 													 
						trait=paste('@',beast2.spec$datetrait.id,sep=''),
						nodetype=beast2.spec$sasky.tree.nodetype,						
						spec=beast2.spec$sasky.tree.parser.spec,
						taxa=paste('@',beast2.spec$tree.taxonset,sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)
		dummy		<- newXMLCommentNode(text=paste("start: Newick starting tree"), parent=bxml.tree, doc=bxml, addFinalizer=T)				
		tmp			<- newXMLNode("input", attrs= list(name='newick'), parent=bxml.tree, doc=bxml, addFinalizer=T)
		dummy		<- newXMLTextNode(text=beast2.spec$starttree.newick, parent=tmp, doc=bxml, addFinalizer=T)
		dummy		<- newXMLCommentNode(text=paste("end: Newick starting tree"), parent=bxml.tree, doc=bxml, addFinalizer=T)		
	}
	if(verbose)	cat(paste("\nadded SA trees for taxonsets, n=", length(beast2.spec$tree.taxonset)))
	bxml
}
######################################################################################
hivc.beast2.add.tree<- function(bxml, beast2.spec, verbose=verbose)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	dummy		<- newXMLNode("tree", attrs= list(	id=beast2.spec$tree.id, 													 
					trait=paste('@',beast2.spec$datetrait.id,sep=''),
					taxonset=paste('@TaxonSet.t:',beast2.spec$tree.taxonset,sep='')), 
			parent=bxml.beast, doc=bxml, addFinalizer=T)	
	if(verbose)	cat(paste("\nadded trees for taxonsets, n=", length(beast2.spec$tree.taxonset)))
	bxml
}
######################################################################################
hivc.beast2.add.treemodel.bdsky<- function(bxml, beast2.spec, verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	beast.treemodel	<- newXMLNode("BirthDeathSkylineModel", attrs= list(	id=beast2.spec$treemodel.id, 
					name=beast2.spec$treemodel.id, 
					tree=paste('@',beast2.spec$tree.id,sep=''),
					spec=beast2.spec$bdsky.spec,
					intervalNumber=as.character(beast2.spec$bdsky.intervalNumber)), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.sprop.id, 
					name="samplingProportion",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.sprop.lower),
					upper=as.character(beast2.spec$bdsky.sprop.upper),
					value=paste(beast2.spec$bdsky.sprop.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.R0.id, 
					name="R0",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.R0.lower),
					upper=as.character(beast2.spec$bdsky.R0.upper),
					value=paste(beast2.spec$bdsky.R0.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)											
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.notInf.id, 
					name="becomeUninfectiousRate",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.notInf.lower),
					upper=as.character(beast2.spec$bdsky.notInf.upper),
					value=paste(beast2.spec$bdsky.notInf.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.origin.id, 
					name="origin",
					lower=as.character(beast2.spec$bdsky.origin.lower),
					upper=as.character(beast2.spec$bdsky.origin.upper),
					value=paste(beast2.spec$bdsky.origin.value)), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	tmp				<- rep("false",4)
	tmp[which(!sapply(c(beast2.spec$bdsky.R0.changepoint.id[1],beast2.spec$bdsky.notInf.changepoint.id[1],beast2.spec$bdsky.sprop.changepoint.id[1]),is.null))]		<- "true"
	dummy			<- newXMLNode("reverseTimeArrays", attrs= list(	id=beast2.spec$bdsky.reverseTimeArrays.id, 
					spec=beast2.spec$bdsky.reverseTimeArrays.spec,
					value=paste(tmp, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.R0.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.R0.changepoint.id, 
						name="birthRateChangeTimes",
						value=paste(beast2.spec$bdsky.R0.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	
	if(!is.null(beast2.spec$bdsky.notInf.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.notInf.changepoint.id, 
						name="deathRateChangeTimes",
						value=paste(beast2.spec$bdsky.notInf.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.sprop.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.sprop.changepoint.id, 
						name="samplingRateChangeTimes",
						value=paste(beast2.spec$bdsky.sprop.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadded BDSKY tree models for taxonsets, n=", length(beast2.spec$tree.taxonset)))
	bxml	
}
######################################################################################
hivc.beast2.add.treemodel.sasky<- function(bxml, beast2.spec, verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	beast.treemodel	<- newXMLNode("BirthDeathSkylineModel", attrs= list(	id=beast2.spec$treemodel.id, 
					name=beast2.spec$treemodel.id, 
					tree=paste('@',beast2.spec$tree.id,sep=''),
					spec=beast2.spec$sasky.spec,
					intervalNumber=as.character(beast2.spec$bdsky.intervalNumber)), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.sprop.id, 
					name="samplingProportion",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.sprop.lower),
					upper=as.character(beast2.spec$bdsky.sprop.upper),
					value=paste(beast2.spec$bdsky.sprop.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.R0.id, 
					name="R0",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.R0.lower),
					upper=as.character(beast2.spec$bdsky.R0.upper),
					value=paste(beast2.spec$bdsky.R0.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)											
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.notInf.id, 
					name="becomeUninfectiousRate",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.notInf.lower),
					upper=as.character(beast2.spec$bdsky.notInf.upper),
					value=paste(beast2.spec$bdsky.notInf.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$sasky.r.id, 
					name="becomeNoninfectiousAfterSamplingProbability",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$sasky.r.lower),
					upper=as.character(beast2.spec$sasky.r.upper),
					value=paste(beast2.spec$sasky.r.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)	
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.origin.id, 
					name="origin",
					lower=as.character(beast2.spec$bdsky.origin.lower),
					upper=as.character(beast2.spec$bdsky.origin.upper),
					value=paste(beast2.spec$bdsky.origin.value)), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	tmp				<- rep("false",4)
	tmp[which(!sapply(c(beast2.spec$bdsky.R0.changepoint.id[1],beast2.spec$bdsky.notInf.changepoint.id[1],beast2.spec$bdsky.sprop.changepoint.id[1]),is.null))]		<- "true"
	dummy			<- newXMLNode("reverseTimeArrays", attrs= list(	id=beast2.spec$bdsky.reverseTimeArrays.id, 
					spec=beast2.spec$bdsky.reverseTimeArrays.spec,
					value=paste(tmp, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.R0.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.R0.changepoint.id, 
						name="birthRateChangeTimes",
						value=paste(beast2.spec$bdsky.R0.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)	
	if(!is.null(beast2.spec$bdsky.notInf.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.notInf.changepoint.id, 
						name="deathRateChangeTimes",
						value=paste(beast2.spec$bdsky.notInf.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.sprop.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.sprop.changepoint.id, 
						name="samplingRateChangeTimes",
						value=paste(beast2.spec$bdsky.sprop.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadded SASKY tree models for taxonsets, n=", length(beast2.spec$tree.taxonset)))
	bxml	
}
######################################################################################
hivc.beast2.init.xml<- function( beast2.spec=NULL, verbose=1)
{
	bxml		<- newXMLDoc(addFinalizer=T)
	bxml.beast	<- newXMLNode("beast", attrs=list(beautitemplate='HIVCLUST', beautistatus='', version="2.0", namespace=paste(beast2.spec$namespace,sep='',collapse=':')), doc=bxml, addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Beta"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Beta, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="ExcludablePrior"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.ExcludablePrior, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Exponential"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Exponential, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="InverseGamma"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.InverseGamma, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="LogNormal"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.LogNormal, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Gamma"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Gamma, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Uniform"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Uniform, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="prior"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.prior, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="LaplaceDistribution"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Laplace, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="OneOnX"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.OneOnX, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Normal"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Normal, parent=dummy, doc=bxml,addFinalizer=T)
	bxml
}	
######################################################################################
hivc.beast2.add.tiplogstem<- function(bxml, beast2.spec, verbose=0)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	taxons		<- xpathApply(bxml.beast, "//sequence", xmlGetAttr, "taxon")
	if(verbose)	cat(paste("\nadd taxonsets for tips, n=", length(taxons)))
	#	add taxonsets	
	tmp			<- newXMLCommentNode(text="Tip Taxonsets, used to log the tip stem height", doc=bxml, addFinalizer=T)
	tmp			<- c(tmp, lapply(seq_along(taxons), function(i)
					{				
						taxonset	<- newXMLNode("taxonset", attrs= list(	id=paste(beast2.spec$tip.taxonset.id.prefix,i,sep=''), spec=beast2.spec$taxonset.spec ), doc=bxml, addFinalizer=T)
						if(length(getNodeSet(bxml, paste("//taxon[@id='",taxons[[i]],"']",sep=''))))
							dummy	<- newXMLNode("taxon", attrs= list(	idref=taxons[[i]], spec=beast2.spec$taxon.spec ), parent=taxonset, doc=bxml, addFinalizer=T)
						else
							dummy	<- newXMLNode("taxon", attrs= list(	id=taxons[[i]], spec=beast2.spec$taxon.spec ), parent=taxonset, doc=bxml, addFinalizer=T)
						taxonset
					}))	
	dummy		<- addChildren(bxml.beast, tmp, at= tail(which(xmlSApply(bxml.beast, xmlName)=="data"),1) )
	#	add MRCA priors	
	if(verbose)	cat(paste("\nadd MRCA priors for tips, n=", length(taxons)))
	bxml.prior	<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior	<- bxml.prior[[1]]	
	dummy		<- newXMLCommentNode(text="start: Tip Taxonset priors, used to log the tip stem height", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_along(taxons), function(i)
			{
				newXMLNode("distribution", attrs= list(	id=paste(beast2.spec$tip.prior.id.prefix,i,sep=''), taxonset=paste('@',beast2.spec$tip.taxonset.id.prefix,i,sep=''), useOriginate='true', monophyletic='false', tree=paste('@',beast2.spec$tree.id,sep=''), spec=beast2.spec$mrca.prior.spec ), parent=bxml.prior, doc=bxml, addFinalizer=T)
			})
	dummy		<- newXMLCommentNode(text="end: Tip Taxonset priors, used to log the tip stem height", parent=bxml.prior, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadd tracelog for stem of tips, n=", length(taxons)))
	#	add log entry to tracelog
	bxml.log	<- getNodeSet(bxml, "//logger[@id='tracelog']")
	if(length(bxml.log)!=1)	stop("unexpected length of bxml.log")
	bxml.log	<- bxml.log[[1]]
	dummy		<- newXMLCommentNode(text="start: Tip Taxonset log entries, used to log the tip stem height", parent=bxml.log, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_along(taxons), function(i)
			{
				newXMLNode("log", attrs= list(	idref=paste(beast2.spec$tip.prior.id.prefix,i,sep='')	), parent=bxml.log, doc=bxml, addFinalizer=T)
			})
	dummy		<- newXMLCommentNode(text="end: Tip Taxonset log entries, used to log the tip stem height", parent=bxml.log, doc=bxml, addFinalizer=T)
	bxml.beast
}
######################################################################################
hivc.beast.add.taxonsets4clusters<- function(bxml, df, xml.monophyly4clusters=1, verbose=1)
{
	bxml.treeModel.id			<- unlist(xpathApply(bxml, "//treeModel[@id]", xmlGetAttr, "id"))
	
	#get taxon sets for each cluster
	btaxonsets.clusters			<- hivc.beast.get.taxonsets4clusters(bxml, df)	
	#get tmrcaStatistic for each cluster
	btmrcaStatistics.clusters	<- hivc.beast.get.tmrcaStatistic(bxml, btaxonsets.clusters, bxml.treeModel.id, includeStem="false") 		
	
	#	modify from template
	bxml.beast					<- getNodeSet(bxml, "//beast")[[1]]
	#	add 'btaxonsets.clusters' after last taxa
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)=="taxa")		
	addChildren(bxml.beast, btaxonsets.clusters, at=bxml.idx[length(bxml.idx)] )
	if(verbose) cat(paste("\nadded taxon sets comprising each cluster, n=",length(btaxonsets.clusters)))
	#	add 'btmrcaStatistics.clusters' after last treeModel
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)=="treeModel")
	addChildren(bxml.beast, btmrcaStatistics.clusters, at=bxml.idx[length(bxml.idx)] )
	if(verbose) cat(paste("\nadded tmrcaStatistics for mrca of each cluster, n=",length(btmrcaStatistics.clusters)))
	#add tmrcaStatistics to log 
	bxml.fileLog				<- getNodeSet(bxml, "//log[@id='fileLog']")
	tmrcaStatistics.id			<- sapply(btmrcaStatistics.clusters, function(x)	xmlGetAttr(x,"id")	)
	tmp							<- lapply(tmrcaStatistics.id, function(x)	newXMLNode("tmrcaStatistic", attrs= list(idref=x), parent=bxml.fileLog, doc=bxml, addFinalizer=T )	)
	if(verbose) cat(paste("\nadded tmrcaStatistics to fileLog"))
	#
	#	get monophylyStatistic and add after last 'tmrcaStatistic'
	#	populate booleanLikelihood with monophylyStatistic constraints
	#
	if(xml.monophyly4clusters)
	{
		bmStatistics.clusters	<- hivc.beast.get.monophylyStatistic(bxml, btaxonsets.clusters, bxml.treeModel.id)
		bxml.idx				<- which(xmlSApply(bxml.beast, xmlName)=="tmrcaStatistic")
		addChildren(bxml.beast, bmStatistics.clusters, at=bxml.idx[length(bxml.idx)] )
		if(verbose) cat(paste("\nadded monophylyStatistics for mrca of each cluster, n=",length(bmStatistics.clusters)))
		hivc.beast.add.monophylylkl(bxml, bmStatistics.clusters)		
		if(verbose) cat(paste("\nadded monophylyStatistics to booleanLikelihood, n=",length(bmStatistics.clusters)))
	}
	
	bxml
}

######################################################################################
#	get a list of monopylyStatistics. assumes each taxonset in 'btaxonsets' is monophyletic
hivc.beast.get.monophylyStatistic<- function(bxml, btaxonsets, treeModel.id) 
{		
	btaxonset.id	<- sapply(btaxonsets, function(x)	xmlGetAttr(x,"id")	)
	ans				<- lapply(btaxonset.id, function(x)
			{
				monophylyStatistic	<- newXMLNode("monophylyStatistic", attrs= list(id=paste("monophyly(",x,")",sep='')), doc=bxml, addFinalizer=T )
				mrca			<- newXMLNode("mrca", parent=monophylyStatistic, doc=bxml, addFinalizer=T)
				newXMLNode("taxa", attrs= list(idref=x), parent=mrca, doc=bxml, addFinalizer=T )
				newXMLNode("treeModel", attrs= list(idref=treeModel.id), parent=monophylyStatistic, doc=bxml, addFinalizer=T )
				monophylyStatistic					
			})
	ans
}	
######################################################################################
#	add a list of monopylyStatistics to the BEAST prior. assumes all monophylyStatistics are referenced in the list 'monophylyStatistics'
hivc.beast.add.monophylylkl<- function(bxml, monophylyStatistics)
{
	bxml.prior				<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior				<- bxml.prior[[1]]
	#see if there is a booleanLikelihood, and if yes use it, otherwise add a new XML node
	bxml.bool.lkl			<- getNodeSet(bxml.prior,"//booleanLikelihood")
	if(length(bxml.bool.lkl)>1)	stop("unexpected length of bxml.bool.lkl")
	if(length(bxml.bool.lkl)<1)
		bxml.bool.lkl		<- newXMLNode("booleanLikelihood", parent=bxml.prior, doc=bxml, addFinalizer=T )
	else
		bxml.bool.lkl		<- bxml.bool.lkl[[1]]
	
	monophylyStatistics.id	<- sapply(monophylyStatistics, function(x)	xmlGetAttr(x,"id")	)
	dummy					<- lapply(monophylyStatistics.id, function(x)
			{
				dummy		<- newXMLNode("monophylyStatistic", attrs= list(idref=x), parent=bxml.bool.lkl, doc=bxml, addFinalizer=T )
			})
	bxml
}
######################################################################################
#	For each taxonset, get a list of tmrcaStatistics. Assumes all tmrcaStatistics share a treeModel.id and the same includeStem attribute
hivc.beast.get.tmrcaStatistic<- function(bxml, btaxonsets, treeModel.id, includeStem="false") 
{		
	btaxonset.id	<- sapply(btaxonsets, function(x)	xmlGetAttr(x,"id")	)
	prefix.id		<- ifelse(includeStem=="false","tmrca","tstem")
	ans				<- lapply(btaxonset.id, function(x)
			{
				tmrcaStatistic	<- newXMLNode("tmrcaStatistic", attrs= list(id=paste(prefix.id,"(",x,")",sep=''), includeStem=includeStem), doc=bxml, addFinalizer=T )
				mrca			<- newXMLNode("mrca", parent=tmrcaStatistic, doc=bxml, addFinalizer=T)
				newXMLNode("taxa", attrs= list(idref=x), parent=mrca, doc=bxml, addFinalizer=T )
				newXMLNode("treeModel", attrs= list(idref=treeModel.id), parent=tmrcaStatistic, doc=bxml, addFinalizer=T )
				tmrcaStatistic					
			})
	ans
}	
######################################################################################
#	For each tip: construct a prior for the corresponding tmrcaStatistics
hivc.beast.get.tipPrior<- function(bxml, df, btmrcaStatistics.tips, xml.prior4stem="uniform", beast.label.negpos=2, beast.label.diagpos=3, beast.label.datepos=4, verbose=1)
{
	if(xml.prior4stem!="uniform")	stop("unexpected xml.tipprior")
	#
	if(verbose)	cat(paste("\nuse tip date found at pos x in label, x=", beast.label.datepos))
	df.height	<- t( sapply( strsplit(df[,BEASTlabel],'_',fixed=1), function(x)  as.numeric( x[c(beast.label.negpos, beast.label.diagpos, beast.label.datepos)]) ) )
	df.height	<- data.table(BEASTlabel=df[,BEASTlabel], NegT=df.height[,1], AnyPos_T1=df.height[,2], PosSeqT=df.height[,3])
	tmp			<- max( df.height[, PosSeqT])		#TODO should this be height or length ?
	df.height	<- df.height[, list(BEASTlabel=BEASTlabel, NegT=tmp-NegT, AnyPos_T1=tmp-AnyPos_T1, PosSeqT=tmp-PosSeqT)]
	#add uniform prior according to last NegT and first PosT
	ans			<- lapply(btmrcaStatistics.tips, function(x)
			{						
				tmrcaStatistics.id	<- xmlGetAttr(x,"id")
				tmp					<- xpathApply(x, "mrca/taxa", xmlGetAttr, "idref" )
				if(length(tmp)!=1)	stop("unexpected length of idref for mrca/taxa") 
				tip					<- as.numeric( substr(tmp[[1]],4,nchar(tmp[[1]])) )												
				bxml.tipprior		<- newXMLNode("uniformPrior", attrs= list(lower=df.height[tip,AnyPos_T1], upper=df.height[tip,NegT]), doc=bxml, addFinalizer=T )
				dummy				<- newXMLNode("statistic", attrs= list(idref=tmrcaStatistics.id), parent=bxml.tipprior, doc=bxml, addFinalizer=T )
				bxml.tipprior
			})
	ans			
}
######################################################################################
#	For each tip: add a taxonset, tmrcaStatistic, and potentially a reference to fileLog to 'bxml'
#	The aim of this as standalone is to log the length of tip stems. 
hivc.beast.add.taxonsets4tips<- function(bxml, df, log=1, verbose=1)
{
	bxml.treeModel.id			<- unlist(xpathApply(bxml, "//treeModel[@id]", xmlGetAttr, "id"))
	
	#get taxon sets for each tip
	btaxonsets.tips				<- hivc.beast.get.taxonsets4tips(bxml, df)	
	#get tmrcaStatistic for each tip
	btmrcaStatistics.tips		<- hivc.beast.get.tmrcaStatistic(bxml, btaxonsets.tips, bxml.treeModel.id, includeStem="true") 			
	#	modify from template
	bxml.beast					<- getNodeSet(bxml, "//beast")[[1]]
	#	add 'btaxonsets.tips' after last taxa
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)=="taxa")		
	addChildren(bxml.beast, btaxonsets.tips, at=bxml.idx[length(bxml.idx)] )
	if(verbose) cat(paste("\nadded taxon sets for each tip, n=",length(btaxonsets.tips)))
	#	add 'btmrcaStatistics.tips' after last treeModel or tmrcaStatistic
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)%in%c("treeModel","tmrcaStatistic"))
	addChildren(bxml.beast, btmrcaStatistics.tips, at=tail(bxml.idx,1) )
	if(verbose) cat(paste("\nadded tmrcaStatistics for stem of each tip, n=",length(btmrcaStatistics.tips)))
	#add tmrcaStatistics to log 
	if(log)
	{
		bxml.fileLog				<- getNodeSet(bxml, "//log[@id='fileLog']")
		tmrcaStatistics.id			<- sapply(btmrcaStatistics.tips, function(x)	xmlGetAttr(x,"id")	)
		tmp							<- lapply(tmrcaStatistics.id, function(x)	newXMLNode("tmrcaStatistic", attrs= list(idref=x), parent=bxml.fileLog, doc=bxml, addFinalizer=T )	)
		if(verbose) cat(paste("\nadded tmrcaStatistics for tip stems to fileLog, n=",length(tmrcaStatistics.id)))
	}
	
	bxml
}
######################################################################################
#	For each tip: add a taxonset, tmrcaStatistic, prior for the tmrcaStatistics and potentially a reference to fileLog to 'bxml'
hivc.beast.add.prior4tips<- function(bxml, df, xml.prior4stem="uniform", beast.label.datepos=4, verbose=1)
{
	#find list of tmrcaStatistics with id containing 'tip'
	btmrcaStatistics.tips		<- getNodeSet(bxml, "//tmrcaStatistic[starts-with(@id,'tstem(tip')]")
	print(btmrcaStatistics.tips)
	#get prior for each tip stem
	bprior.tips					<- hivc.beast.get.tipPrior(bxml, df, btmrcaStatistics.tips, xml.prior4stem=xml.prior4stem, beast.label.datepos=beast.label.datepos, verbose=verbose)		
	#	add 'bprior.tips' to prior
	bxml.prior					<- getNodeSet(bxml, "//prior[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected number of //prior[@id='prior']")
	bxml.prior					<- bxml.prior[[1]]
	addChildren(bxml.prior, bprior.tips)
	if(verbose) cat(paste("\nadded priors for stem of each tip, n=",length(bprior.tips)))	
	bxml
}
######################################################################################
hivc.beast2.add.startingtree.random<- function(bxml, beast2.spec, verbose=0)
{
	bxml.run	<- getNodeSet(bxml, "//run")[[1]]
	if(length(bxml.run)!=1)	stop("unexpected length of bxml.run")
	if(verbose)	cat(paste('\nadd random starting tree'))
	dummy		<- newXMLCommentNode(text=paste("Random starting tree"), parent=bxml.run, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("init", attrs= list(	id=beast2.spec$starttree.id, estimate='false', initial=beast2.spec$tree.id, 
					taxa=paste('@',beast2.spec$rndstarttree.taxonset,sep=''), taxonset=paste('@TaxonSet.t:',beast2.spec$rndstarttree.taxonset,sep='')), parent=bxml.run, doc=bxml, addFinalizer=T)		
	tmp			<- newXMLNode("populationModel", attrs= list(	id=paste('ConstantPopulation',beast2.spec$tree.id,sep=''), spec='ConstantPopulation'	), parent=tmp, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("parameter", attrs= list(	id=paste('PopSize',beast2.spec$tree.id,sep=''), name='popSize', value='1.0'), parent=tmp, doc=bxml, addFinalizer=T)
	bxml
}
######################################################################################
hivc.beast2.add.startingtree.newick<- function(bxml, beast2.spec, verbose=0)
{
	bxml.run	<- getNodeSet(bxml, "//run")[[1]]
	if(length(bxml.run)!=1)	stop("unexpected length of bxml.run")
	if(verbose)	cat(paste('\nadd newick starting tree'))
	dummy		<- newXMLCommentNode(text=paste("start: Newick starting tree"), parent=bxml.run, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("init", attrs= list(	id=beast2.spec$starttree.id, initial=beast2.spec$tree.id, IsLabelledNewick=beast2.spec$starttree.islabelledtree, spec=beast2.spec$starttree.spec), parent=bxml.run, doc=bxml, addFinalizer=T)		
	tmp			<- newXMLNode("input", attrs= list(name='newick'), parent=tmp, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(text=beast2.spec$starttree.newick, parent=tmp, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("end: Newick starting tree"), parent=bxml.run, doc=bxml, addFinalizer=T)
	bxml
}
######################################################################################
hivc.beast2.add.sasky.serialpriors<- function(bxml, beast2.spec, verbose=1)
{
	bxml.prior	<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior	<- bxml.prior[[1]]	
	dummy		<- newXMLCommentNode(text="start: serial SASKY priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text="The SASKY model prior.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("distribution", attrs= list(	id=paste("prior",beast2.spec$treemodel.id,sep='-'),spec=beast2.spec$compoundprior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)	
	dummy		<- newXMLNode("distribution", attrs= list(	idref=beast2.spec$treemodel.id), parent=tmp, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text="The origin prior.", parent=bxml.prior, doc=bxml, addFinalizer=T)	
	tmp			<- newXMLNode("prior", attrs= list(	id=paste("sprior",beast2.spec$bdsky.origin.id,sep='-'),
					name="distribution",
					x= paste('@',beast2.spec$bdsky.origin.id,sep='')), parent=bxml.prior, doc=bxml, addFinalizer=T)	
	dummy		<- hivc.beast2.get.prior(beast2.spec$bdsky.origin.prior, tmp, xmlAttrs(tmp)["id"], bxml)
	dummy		<- newXMLCommentNode(text="The serial samplingProb priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.sprop.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.sprop.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$bdsky.sprop.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})
	if(verbose)	cat(paste("\nadded serial samplingProb priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial R0 priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.R0.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.R0.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$bdsky.R0.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})												
	if(verbose)	cat(paste("\nadded serial R0 priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial becomingUninfectious priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.notInf.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.notInf.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$bdsky.notInf.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})	
	if(verbose)	cat(paste("\nadded serial becomingUninfectious priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial becomeNoninfectiousAfterSamplingProbability priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$sasky.r.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$sasky.r.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$sasky.r.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})	
	if(verbose)	cat(paste("\nadded serial becomeNoninfectiousAfterSamplingProbability priors, n=", length(dummy)))	
	dummy		<- newXMLCommentNode(text="end: serial SASKY priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	bxml.prior
}
######################################################################################
hivc.beast2.get.startingtree<- function(ph, df, beast2.spec, verbose=1)
{
	require(adephylo)
	if(verbose) cat(paste("\ncreate startingTree with root height=",beast2.spec$starttree.rootHeight))
	tmp					<- match( setdiff( ph$tip.label, df[,FASTASampleCode] ), ph$tip.label)
	ph.start			<- drop.tip(ph, tmp)		
	ph.start$node.label	<- NULL
	setkey(df, FASTASampleCode)
	ph.start$tip.label	<- df[ph.start$tip.label,][,BEASTlabel]
	if(verbose) cat(paste("\nselected tips for startingTree, n=",Ntip(ph.start)))
	#	adjust rootHeight to 'beast.rootHeight'	
	ph.start$edge.length<- ph.start$edge.length	*	beast2.spec$starttree.rootHeight / max(distRoot(ph.start))
	if(verbose) cat(paste("\nadjusted root height=",max(distRoot(ph.start))))	
	if(beast2.spec$starttree.rootHeight>=beast2.spec$bdsky.origin.value)	stop('Detected invalid rootHeight. Must be smaller than origin.value.')
	write.tree( ph.start )		
}
######################################################################################
hivc.beast2.get.specifications	<- function(xml.dir=NA, xml.filename=NA, mcmc.length=20e6, bdsky.intervalNumber=4, alignment.filter=NA, tip.log.stem=FALSE, cluster.log=FALSE, cluster.monophyletic=FALSE)
{
	beast2.spec<- list()	
	beast2.spec$namespace						<- c("beast.core","beast.evolution.alignment","beast.evolution.tree.coalescent","beast.core.util","beast.evolution.nuc","beast.evolution.operators","beast.evolution.sitemodel","beast.evolution.substitutionmodel","beast.evolution.likelihood","beast.evolution.speciation","beast.core.parameter")
	beast2.spec$xml.dir							<- xml.dir		
	beast2.spec$xml.filename					<- xml.filename
	beast2.spec$pool.cnts.requested				<- rep(NA, bdsky.intervalNumber)
	beast2.spec$pool.ntip						<- 130
	beast2.spec$pool.fNegT						<- 0.8
	beast2.spec$map.Beta						<- "beast.math.distributions.Beta"
	beast2.spec$map.Exponential					<- "beast.math.distributions.Exponential"
	beast2.spec$map.ExcludablePrior				<- "beast.math.distributions.ExcludablePrior"
	beast2.spec$map.InverseGamma				<- "beast.math.distributions.InverseGamma"
	beast2.spec$map.LogNormal					<- "beast.math.distributions.LogNormalDistributionModel"
	beast2.spec$map.Gamma						<- "beast.math.distributions.Gamma"
	beast2.spec$map.Uniform						<- "beast.math.distributions.Uniform"
	beast2.spec$map.prior						<- "beast.math.distributions.Prior"
	beast2.spec$map.Laplace						<- "beast.math.distributions.LaplaceDistribution"
	beast2.spec$map.OneOnX						<- "beast.math.distributions.OneOnX"
	beast2.spec$map.Normal						<- "beast.math.distributions.Normal"
	beast2.spec$mcmc.length						<- mcmc.length
	beast2.spec$beast.label.datepos				<- 4
	beast2.spec$beast.label.sep					<- '_'		
	beast2.spec$data.missing					<- "-?"
	beast2.spec$data.dataType					<- 'nucleotide'
	beast2.spec$alignment.spec					<- "FilteredAlignment"
	beast2.spec$alignment.taxa.spec				<- "TaxonSet"
	beast2.spec$alignment.filter				<- alignment.filter
	beast2.spec$taxon.spec						<- "Taxon"
	beast2.spec$taxonset.spec					<- "TaxonSet"	
	beast2.spec$mrca.prior.spec					<- "beast.math.distributions.MRCAPrior"
	if(length(alignment.filter)==1 && is.na(alignment.filter))
	{
		beast2.spec$data.id						<- 'ds'
		beast2.spec$alignment.id				<- 'ds'
	}
	else
	{
		beast2.spec$data.id						<- 'data'
		beast2.spec$alignment.id				<- paste('ds', seq_along(beast2.spec$alignment.filter), sep='_')
	}
	beast2.spec$tree.taxonset					<- beast2.spec$alignment.id[1]
	beast2.spec$tree.id							<- paste('Tree',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$sequence.totalcount				<- 4
	beast2.spec$datetrait.spec					<- "beast.evolution.tree.TraitSet"
	beast2.spec$datetrait.taxa.spec				<- "TaxonSet"
	beast2.spec$datetrait.taxa					<- paste("TaxonSet.t",beast2.spec$alignment.id[1],sep=':')
	beast2.spec$datetrait.id					<- paste("dateTrait.t",beast2.spec$alignment.id[1],sep=':')
	beast2.spec$datetrait.units					<- "year"
	beast2.spec$datetrait.traitname				<- "date"
	beast2.spec$treemodel						<- "BirthDeathSkylineModel"
	beast2.spec$treemodel.id					<- paste("birthDeath",beast2.spec$tree.taxonset,sep='.t:')	
	beast2.spec$bdsky.spec						<- "beast.evolution.speciation.BirthDeathSkylineModel"
	beast2.spec$bdsky.prior.spec				<- "beast.math.distributions.ExcludablePrior"
	beast2.spec$bdsky.intervalNumber			<- bdsky.intervalNumber
	beast2.spec$bdsky.origin.id					<- paste('originS',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.origin.value				<- rep(38, length(beast2.spec$bdsky.origin.id))
	beast2.spec$bdsky.origin.lower				<- 0.0
	beast2.spec$bdsky.origin.upper				<- 1000.0
	beast2.spec$bdsky.origin.prior				<- "Uniform/20.0/40.0"
	beast2.spec$bdsky.sprop.id					<- paste('samplingProportionS',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.sprop.value				<- rep(0.4, beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.sprop.lower				<- 0.0
	beast2.spec$bdsky.sprop.upper				<- 1.0
	beast2.spec$bdsky.sprop.prior				<- rep("Uniform/0.2/1.0",beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.sprop.changepoint.id		<- paste('samplingRateChangeTimes',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.sprop.changepoint.value	<- c(9.596, 5.596, 1.596, 0.)	
	beast2.spec$bdsky.R0.id						<- paste('R0S',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.R0.value					<- rep(1.2, beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.R0.lower					<- 0.0
	beast2.spec$bdsky.R0.upper					<- 10.0
	beast2.spec$bdsky.R0.prior					<- rep("Gamma/1.5/1.5/0",beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.R0.changepoint.id			<- paste('birthRateChangeTimes',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.R0.changepoint.value		<- c(9.596, 5.596, 1.596, 0.)	
	beast2.spec$bdsky.notInf.id					<- paste('becomeUninfectiousRateS',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.notInf.value				<- rep(0.1, beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.notInf.lower				<- 0.0
	beast2.spec$bdsky.notInf.upper				<- 10.0
	beast2.spec$bdsky.notInf.prior				<- rep("OneOnX/0",beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.notInf.changepoint.id		<- paste('deathRateChangeTimes',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.notInf.changepoint.value	<- c(9.596, 5.596, 1.596, 0.)
	beast2.spec$bdsky.reverseTimeArrays.id		<- paste('reverseTimeArrays',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.reverseTimeArrays.spec	<- "parameter.BooleanParameter"	
	beast2.spec$sasky.spec						<- "beast.evolution.speciation.SABDSkylineModel"		
	beast2.spec$sasky.tree.nodetype				<- "beast.evolution.tree.ZeroBranchSANode"
	beast2.spec$sasky.tree.parser.spec			<- "beast.util.ZeroBranchSATreeParser"
	beast2.spec$sasky.tree.cluster.spec			<- "beast.util.ClusterZBSATree"
	beast2.spec$sasky.tree.cluster.type			<- "upgma"
	beast2.spec$sasky.r.id						<- paste('r',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$sasky.r.value					<- rep(0.5, beast2.spec$bdsky.intervalNumber)
	beast2.spec$sasky.r.lower					<- 0.0
	beast2.spec$sasky.r.upper					<- 1.0
	beast2.spec$sasky.r.prior					<- rep("Uniform/0.0/1.0",beast2.spec$bdsky.intervalNumber)
	beast2.spec$sasky.r.changepoint.id			<- paste('rChangeTimes',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$sasky.r.changepoint.value		<- c(9.596, 5.596, 1.596, 0.)		
	beast2.spec$compoundprior.spec				<- "util.CompoundDistribution"
	beast2.spec$tip.log.stem					<- tip.log.stem
	beast2.spec$tip.taxonset.id.prefix			<- 'tip'
	beast2.spec$tip.prior.id.prefix				<- 'tipprior'
	beast2.spec$cluster.taxonset.id.prefix		<- 'c'
	beast2.spec$cluster.prior.id.prefix			<- 'cprior'
	beast2.spec$cluster.monophyletic			<- ifelse(cluster.monophyletic, 'true', 'false')
	beast2.spec$cluster.log						<- cluster.log	
	beast2.spec$starttree.rootHeight			<- 35
	beast2.spec$starttree.usingDates			<- 'true'
	beast2.spec$starttree.brlunits				<- 'years'
	beast2.spec$starttree.islabelledtree		<- 'true'
	beast2.spec$starttree.id					<- paste('Starting',beast2.spec$tree.id,sep='')	
	beast2.spec$starttree.newick				<- NA
	beast2.spec$starttree.spec					<- 'beast.util.TreeParser'	
	beast2.spec$rndstarttree.taxonset			<- beast2.spec$tree.taxonset
	beast2.spec$rndstarttree.spec				<- 'beast.evolution.tree.RandomTree'
	beast2.spec
} 
######################################################################################
hivc.beast2.get.xml<- function(	bxml.template, seq.PROT.RT, df, beast2.spec, ph=NULL, verbose=1)
{	
	require(XML)
	#	init XML
	bxml		<- hivc.beast2.init.xml( beast2.spec, verbose=verbose)
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]	
	#	add data	
	dummy		<- hivc.beast2.add.data(bxml, seq.PROT.RT, df, beast2.spec, verbose=verbose)
	#	add alignment, alignment filters, alignment taxonsets
	dummy		<- hivc.beast2.add.alignment(bxml, beast2.spec, verbose=verbose)
	#	add tip dates
	dummy		<- hivc.beast2.add.datetrait(bxml, df, beast2.spec, verbose=verbose)
	#	add tree and starting tree for alignment
	if(beast2.spec$treemodel=="BirthDeathSkylineModel")
		dummy	<- hivc.beast2.add.tree(bxml, beast2.spec, verbose=verbose)
	if(beast2.spec$treemodel=="SampledAncestorSkylineModel")
		dummy	<- hivc.beast2.add.satree(bxml, beast2.spec, verbose=verbose)
	#	add tree model -- TODO move starting tree into here
	if(beast2.spec$treemodel=="BirthDeathSkylineModel")
		dummy	<- hivc.beast2.add.treemodel.bdsky(bxml, beast2.spec, verbose=verbose)
	#	add tree model and starting tree
	if(beast2.spec$treemodel=="SampledAncestorSkylineModel")
		dummy	<- hivc.beast2.add.treemodel.sasky(bxml, beast2.spec, verbose=verbose)	
	#	copy branchRateModel from template
	tmp			<- getNodeSet(bxml.template, "//branchRateModel")
	if(length(tmp)!=1) stop("unexpected number of //branchRateModel")
	dummy<- addChildren( bxml.beast, xmlClone( tmp[[1]], addFinalizer=T, doc=bxml ) )
	if(verbose)	cat(paste("\nadded branchRateModel from template, size=", xmlSize(tmp[[1]])))	
	#	copy siteModel from template
	tmp			<- getNodeSet(bxml.template, "//siteModel")
	if(length(tmp)<1) stop("not found //siteModel")
	dummy		<- lapply(seq_along(tmp),function(i)	addChildren( bxml.beast, xmlClone( tmp[[i]], addFinalizer=T, doc=bxml ) )		) 		 
	if(verbose)	cat(paste("\nadded siteModel from template, numer=", length(tmp)))
	#	copy run from template
	tmp			<- getNodeSet(bxml.template, "//run")
	if(length(tmp)!=1) stop("unexpected number of //run")
	dummy		<- addChildren( bxml.beast, xmlClone( tmp[[1]], addFinalizer=T, doc=bxml ) )
	if(verbose)	cat(paste("\nadded run from template, size=", xmlSize(tmp[[1]])))
	#	add initial tree now if bdsky
	if(beast2.spec$treemodel=="BirthDeathSkylineModel" && is.na(beast2.spec$starttree.newick))
		dummy	<- hivc.beast2.add.startingtree.random(bxml, beast2.spec, verbose=verbose)
	if(beast2.spec$treemodel=="BirthDeathSkylineModel" && !is.na(beast2.spec$starttree.newick))
		dummy	<- hivc.beast2.add.startingtree.newick(bxml, beast2.spec, verbose=verbose)		
	#	add tree model prior
	if(beast2.spec$treemodel=="BirthDeathSkylineModel")
		dummy	<- hivc.beast2.add.bdsky.serialpriors(bxml, beast2.spec, verbose=verbose)
	if(beast2.spec$treemodel=="SampledAncestorSkylineModel")
		dummy	<- hivc.beast2.add.sasky.serialpriors(bxml, beast2.spec, verbose=verbose)
	
	#	reset output fileNames
	bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
	tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))	
	tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(beast2.spec$xml.filename, '.', rev(x)[1], sep=''))		
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
	if(verbose)	cat(paste("\nchanged trunk filename to", beast2.spec$xml.filename))
	#	reset chain length and logEvery
	bxml.onodes	<- getNodeSet(bxml, "//*[@chainLength]")
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["chainLength"]<- sprintf("%d",beast2.spec$mcmc.length)		})
	if(verbose)	cat(paste("\nchanged chain length to", beast2.spec$mcmc.length))
	bxml.onodes	<- getNodeSet(bxml, "//*[@logEvery]")
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["logEvery"]<- sprintf("%d",round(beast2.spec$mcmc.length/1e4))		})
	if(verbose)	cat(paste("\nchanged logEvery to", round(beast2.spec$mcmc.length/1e4)))	
	#	add logs for height of tip stem
	if(verbose)	cat(paste("\nadd log for height of tip stems=", beast2.spec$tip.log.stem))
	if(beast2.spec$tip.log.stem)
		dummy	<- hivc.beast2.add.tiplogstem(bxml, beast2.spec, verbose=verbose)
	#	add cluster taxonsets
	if(verbose)	cat(paste("\nadd taxonsets for clusters=", beast2.spec$cluster.monophyletic))
	if(beast2.spec$cluster.monophyletic)
		dummy	<- hivc.beast2.add.cluster.taxonsets(bxml, df, beast2.spec, verbose=verbose)
	
	bxml	
}
######################################################################################
#	create xml file from btemplate and seq.PROT.RT, using seq in df 
# 	beast.label.datepos= 4; beast.label.sep= '_'; beast.date.direction= "forwards"; beast.date.units= "years"; verbose=1; xml.prior4tipstem="uniform"; xml.resetTipDate2LastDiag=1
hivc.beast.get.xml<- function(	btemplate, seq.PROT.RT, df, file, ph=NULL, xml.monophyly4clusters=0, xml.taxon4tipstem=0, xml.prior4tipstem=NA, 
		beast.label.datepos= 4, beast.label.sep= '_', beast.date.direction= "forwards", beast.usingDates="true", beast.date.units= "years", beast.mcmc.chainLength=50000000, verbose=1)
{
	
	bxml		<- newXMLDoc(addFinalizer=T)
	bxml.beast	<- newXMLNode("beast", doc=bxml, addFinalizer=T)
	newXMLCommentNode(text=paste("Generated by HIVCLUST from template",file), parent=bxml.beast, doc=bxml, addFinalizer=T)
	#	add new set of sequences
	dummy		<- hivc.beast.add.seq(bxml, df, seq.PROT.RT, beast.label.datepos=beast.label.datepos, beast.label.sep=beast.label.sep, beast.date.direction=beast.date.direction, beast.date.units=beast.date.units, verbose=verbose)
	#	copy everything after alignment up to but not including constantSize from template
	bt.beast	<- getNodeSet(btemplate, "//beast")[[1]]
	dummy		<- sapply(seq.int( which( xmlSApply(bt.beast, xmlName)=="alignment" )+1, which( xmlSApply(bt.beast, xmlName)=="patterns" ) ), function(i)
			{
				if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
					dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
				else
					dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
			})
	#	add startingTree
	if(!is.null(ph))		#	add startingTree if provided
		dummy	<- hivc.beast.add.startingtree(bxml, ph, df, beast.usingDates=beast.usingDates)		
	else					#	otherwise copy upgmaTree
	{
		dummy		<- sapply(seq.int( which( xmlSApply(bt.beast, xmlName)=="constantSize" )-2, which( xmlSApply(bt.beast, xmlName)=="upgmaTree" ) ), function(i)
				{
					if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
						dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
					else
						dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
				})
	}
	#	copy everything from 'treeModel'-1 from template
	dummy		<- sapply(seq.int( which( xmlSApply(bt.beast, xmlName)=="treeModel" )-1, xmlSize(bt.beast) ), function(i)
			{
				if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
					dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
				else
					dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
			})
	#	if user-provided startingTree, reset <upgmaTree idref="startingTree"/> to <newick idref="startingTree"/>
	if(!is.null(ph))
	{
		tmp					<- getNodeSet(bxml, "//*[@idref='startingTree']")
		if(length(tmp)!=1) stop("unexpected number of //*[@idref='startingTree']")
		xmlName(tmp[[1]])	<- "newick"
	}
	# 	reset dimension of GMRF likelihood
	tmp			<- getNodeSet(bxml, "//*[@id='skyride.logPopSize']")
	if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.logPopSize'")
	tmp			<- tmp[[1]]
	xmlAttrs(tmp)["dimension"]	<-	nrow(df)-1  
	tmp			<- getNodeSet(bxml, "//*[@id='skyride.groupSize']")
	if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.groupSize'")
	tmp			<- tmp[[1]]
	xmlAttrs(tmp)["dimension"]	<-	nrow(df)-1
	#	if uniform prior for rootheight, set minimum to earliest sample time in data set
	dummy		<- hivc.beast.adjust.rootheightprior(bxml, df, verbose=verbose)
	dummy		<- hivc.beast.adjust.mcmc(bxml, beast.mcmc.chainLength=beast.mcmc.chainLength, verbose=verbose)
	#	for tips, add taxon sets and tmrcaStatistics
	if(xml.taxon4tipstem)
		dummy	<- hivc.beast.add.taxonsets4tips(bxml, df, log=1, verbose=1)
	if(!is.na(xml.prior4tipstem))
		dummy	<- hivc.beast.add.prior4tips(bxml, df, xml.prior4stem=xml.prior4tipstem, beast.label.datepos=beast.label.datepos, verbose=verbose)
	
	#	for clusters, add taxon sets and tmrcaStatistics 
	dummy		<- hivc.beast.add.taxonsets4clusters(bxml, df, xml.monophyly4clusters=xml.monophyly4clusters, verbose=verbose)		
	#	reset output fileNames
	bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
	tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
	tmp			<- gsub("(time).","time",tmp,fixed=1)
	tmp			<- gsub("(subst).","subst",tmp,fixed=1)
	tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(file, '.', x[2], sep=''))		
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
	#
	bxml
}	
######################################################################################
#	adjust mcmc BEAST XML element
hivc.beast.adjust.mcmc<- function(bxml, beast.mcmc.chainLength=50000000, verbose=1)
{
	tmp			<- getNodeSet(bxml, "//*[@id='mcmc']")
	if(length(tmp)!=1)	stop("unexpected number of *[@id='mcmc']")
	tmp			<- tmp[[1]]
	if(verbose)	cat(paste("\nSet MCMC chainLength to",beast.mcmc.chainLength))
	xmlAttrs(tmp)["chainLength"]	<-	sprintf("%d",beast.mcmc.chainLength)
	bxml
}
######################################################################################
#	if rootheight prior uniform, sets lower bound to earliest sampling time in data set
hivc.beast.adjust.rootheightprior<- function(bxml, df, verbose=1)
{
	bxml.prior	<- getNodeSet(bxml, "//uniformPrior[descendant::parameter[@idref='treeModel.rootHeight']]")
	if(length(bxml.prior)>1)	stop("unexpected length of //uniformPrior[descendant::parameter[@idref='treeModel.rootHeight']]")
	if(length(bxml.prior)==1)
	{
		tmp								<- df[,range(AnyPos_T1)]
		if(verbose)	cat(paste("\nfound uniformPrior for treeModel.rootHeight. Range of tip dates is",tmp[1],tmp[2]))
		tmp								<- difftime(tmp[2],tmp[1], units="days")	
		bxml.prior						<- bxml.prior[[1]]
		if(verbose)	cat(paste("\nset lower bound of uniformPrior for treeModel.rootHeight to", floor( tmp/365 )))
		xmlAttrs(bxml.prior)["lower"]	<- floor( tmp/365 )	
	}
	bxml
}
######################################################################################
#	For each cluster, create a taxonset. Assumes df has BEASTlabel and cluster
hivc.beast.get.taxonsets4clusters	<- function(bxml, df)
{	
	ans	<- lapply( unique( df[,cluster] ), function(clu)
			{							
				taxonset	<- newXMLNode("taxa", attrs= list(id=paste("c",clu,sep='')), doc=bxml, addFinalizer=T )
				tmp<- df[cluster==clu,][,BEASTlabel]
				tmp<- lapply(tmp, function(x)		newXMLNode("taxon", attrs= list(idref=x), parent=taxonset, doc=bxml, addFinalizer=T )	)
				taxonset
			})
	ans		
}
######################################################################################
#	For each tip, create a taxonset. Assumes df has BEASTlabel
hivc.beast.get.taxonsets4tips	<- function(bxml, df)
{	
	ans	<- lapply( seq_len(nrow(df)), function(i)
			{							
				taxonset	<- newXMLNode("taxa", attrs= list(id=paste("tip",i,sep='')), doc=bxml, addFinalizer=T )
				newXMLNode("taxon", attrs= list(idref=df[i,BEASTlabel]), parent=taxonset, doc=bxml, addFinalizer=T )
				taxonset
			})
	ans		
}
######################################################################################
# 	extract starting tree from 'ph' by tips in 'df'. Only keeps tree topology and resets branch lengths so that the maximum root distance is 'beast.rootHeight'
#	beast.rootHeight= 35; beast.usingDates= "false"; beast.newickid= "startingTree"
hivc.beast.add.startingtree<- function(bxml, ph, df, beast.rootHeight= 35, beast.usingDates="true", beast.newickid= "startingTree", beast.brlunits="years", verbose=1)
{
	require(adephylo)
	if(verbose) cat(paste("\ncreate startingTree with root height=",beast.rootHeight))
	tmp					<- setdiff( ph$tip.label, df[,FASTASampleCode] )
	tmp					<- match( tmp, ph$tip.label)
	ph.start			<- drop.tip(ph, tmp)		
	ph.start$node.label	<- NULL
	setkey(df, FASTASampleCode)
	ph.start$tip.label	<- df[ph.start$tip.label,][,BEASTlabel]
	if(verbose) cat(paste("\nselected tips for startingTree, n=",Ntip(ph.start)))
	#	adjust rootHeight to 'beast.rootHeight'
	tmp					<- beast.rootHeight / max(distRoot(ph.start))
	ph.start$edge.length<- ph.start$edge.length*tmp
	#	compute adjusted branch lengths for each tip: midpoint within NegT and AnyPos_T1
	if(0)
	{
		df.length			<- suppressWarnings( t( sapply( strsplit(df[,BEASTlabel],'_',fixed=1), function(x)  as.numeric( x[2:4]) ) ) )
		df.length			<- data.table(BEASTlabel=df[,BEASTlabel], NegT=df.length[,1], AnyPos_T1=df.length[,2], PosSeqT=df.length[,3])
		tmp					<- max( df.length[, PosSeqT])		#TODO should this be height or length ?
		df.length			<- df.length[, list(BEASTlabel=BEASTlabel, brl=(AnyPos_T1-NegT)/2+PosSeqT-AnyPos_T1)]
		setkey(df.length, BEASTlabel)
		#	adjust stem of each tip to be within NegT and AnyPos_T1
		tmp							<- sapply(seq_len(Ntip(ph.start)), function(x) which( ph.start$edge[,2]==x ) )
		ph.start$edge.length[ tmp ]	<- df.length[ph.start$tip.label,][,brl]
		if(verbose) cat(paste("\nadjusted branch lengths of tips to be within NegT and AnyPos_T1. New root height is",max(distRoot(ph.start))))
	}
	#	write ph.start as newick tree to bxml
	tmp					<- write.tree( ph.start )	
	bxml.beast			<- getNodeSet(bxml, "//beast")[[1]]
	dummy				<- newXMLCommentNode(text="The user-specified starting tree in a newick tree format", parent=bxml.beast, doc=bxml, addFinalizer=T)
	bxml.startingTree	<- newXMLNode("newick", attrs= list(id=beast.newickid, usingDates=beast.usingDates, units=beast.brlunits), parent= bxml.beast, doc=bxml, addFinalizer=T)
	dummy				<- newXMLTextNode(text=tmp, parent=bxml.startingTree, doc=bxml, addFinalizer=T) 
	bxml
}
######################################################################################
hivc.beast2.add.cluster.taxonsets<- function(bxml, df, beast2.spec, verbose=0)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	clusters	<- df[,unique(cluster)]	
	if(verbose)	cat(paste("\nadd taxonsets for clusters, n=", length( clusters )))
	#	add taxonsets	
	tmp			<- newXMLCommentNode(text="start: The cluster taxonsets, used to log their TMRCA or to enforce monophyly", doc=bxml, addFinalizer=T)
	tmp			<- c(tmp, lapply(clusters, function(clu)
					{				
						taxonset	<- newXMLNode("taxonset", attrs= list(	id=paste(beast2.spec$cluster.taxonset.id.prefix,clu,sep=''), spec=beast2.spec$taxonset.spec ), doc=bxml, addFinalizer=T)
						dummy		<- sapply( subset(df, cluster==clu)[, BEASTlabel], function(x)
								{
									if(length(getNodeSet(bxml, paste("//taxon[@id='",x,"']",sep=''))))									
										newXMLNode("taxon", attrs= list(	idref=x, spec=beast2.spec$taxon.spec ), parent=taxonset, doc=bxml, addFinalizer=T)
									else
										newXMLNode("taxon", attrs= list(	id=x, spec=beast2.spec$taxon.spec ), parent=taxonset, doc=bxml, addFinalizer=T)
								})										
						taxonset
					}))	
	tmp			<- c(tmp, newXMLCommentNode(text="end: The cluster taxonsets, used to log their TMRCA or to enforce monophyly", doc=bxml, addFinalizer=T))
	dummy		<- addChildren(bxml.beast, tmp, at= tail(which(xmlSApply(bxml.beast, xmlName)=="data"),1) )
	#	add MRCA priors	
	if(verbose)	cat(paste("\nadd MRCA priors for clusters, n=", length(clusters)))
	bxml.prior	<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior	<- bxml.prior[[1]]	
	dummy		<- newXMLCommentNode(text="start: Cluster taxonset priors, used to log their TMRCA or to enforce monophyly", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(clusters, function(clu)
			{
				newXMLNode("distribution", attrs= list(	id=paste(beast2.spec$cluster.prior.id.prefix,clu,sep=''), taxonset=paste('@',beast2.spec$cluster.taxonset.id.prefix,clu,sep=''), useOriginate='false', monophyletic=beast2.spec$cluster.monophyletic, tree=paste('@',beast2.spec$tree.id,sep=''), spec=beast2.spec$mrca.prior.spec ), parent=bxml.prior, doc=bxml, addFinalizer=T)
			})
	dummy		<- newXMLCommentNode(text="end: Cluster taxonset priors, used to log their TMRCA or to enforce monophyly", parent=bxml.prior, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadd tracelog for cluster tmrca=", beast2.spec$cluster.log))
	if(beast2.spec$cluster.log)
	{
		if(verbose)	cat(paste("\nadd tracelog for cluster tmrca, n=", length(clusters)))				
		#	add log entry to tracelog
		bxml.log	<- getNodeSet(bxml, "//logger[@id='tracelog']")
		if(length(bxml.log)!=1)	stop("unexpected length of bxml.log")
		bxml.log	<- bxml.log[[1]]
		dummy		<- newXMLCommentNode(text="start: Cluster taxonset log entries, used to log their TMRCA", parent=bxml.log, doc=bxml, addFinalizer=T)
		dummy		<- lapply(clusters, function(clu)
				{
					newXMLNode("log", attrs= list(	idref=paste(beast2.spec$cluster.prior.id.prefix,clu,sep='')	), parent=bxml.log, doc=bxml, addFinalizer=T)
				})
		dummy		<- newXMLCommentNode(text="end: Cluster taxonset log entries, used to log their TMRCA", parent=bxml.log, doc=bxml, addFinalizer=T)		
	}
	bxml.beast	
}
######################################################################################
hivc.beast2.get.prior<- function(args, parent, parentid, bxml)
{
	args	<- strsplit(args,'/')[[1]]
	if(args[1]=="Uniform")
		prior	<- newXMLNode("Uniform", attrs= list( id= paste('U',parentid,sep='-'), name="distr", lower=args[2], upper=args[3]), parent=parent, doc=bxml, addFinalizer=T)											
	else if(args[1]=="OneOnX")
		prior	<- newXMLNode("OneOnX", attrs= list( id= paste('OneOnX.',parentid,sep=''), name="distr"), parent=parent, doc=bxml, addFinalizer=T)
	else if(args[1]=="Exponential")
	{
		prior	<- newXMLNode("Exponential", attrs= list( id= paste('Exponential',parentid,sep='-'), name="distr", offset=args[3]), parent=parent, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pExponential',parentid,sep='-'), name="mean", value=args[2], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
	}
	else if(args[1]=="Beta")
	{
		prior	<- newXMLNode("Beta", attrs= list( id= paste('Beta',parentid,sep='-'), name="distr", offset=args[4]), parent=parent, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pBeta1',parentid,sep='-'), lower="0.0", upper="10.0", name="alpha", value=args[2], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pBeta2',parentid,sep='-'), lower="0.0", upper="10.0", name="beta", value=args[3], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
	}
	else if(args[1]=="Gamma")
	{
		prior	<- newXMLNode("Gamma", attrs= list( id= paste('Gamma',parentid,sep='-'), name="distr", offset=args[4]), parent=parent, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pGamma1',parentid,sep='-'), name="alpha",  value=args[2], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pGamma2',parentid,sep='-'), name="beta",  value=args[3], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
	}
	else if(args[1]=="LogNormal")
	{
		prior	<- newXMLNode("LogNormal", attrs= list( id= paste('LogNormal',parentid,sep='-'), name="distr", offset=args[4], meanInRealSpace=args[5]), parent=parent, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pLogNormal1',parentid,sep='-'), name="M",  value=args[2], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pLogNormal2',parentid,sep='-'), name="S",  value=args[3], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
	}
	else	stop("prior not implemented")	
	prior
}
######################################################################################
hivc.beast2.add.bdsky.serialpriors<- function(bxml, beast2.spec, verbose=1)
{
	bxml.prior	<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior	<- bxml.prior[[1]]	
	dummy		<- newXMLCommentNode(text="start: serial BDSKY priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text="The BDSKY model prior.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("distribution", attrs= list(	id=paste("prior",beast2.spec$treemodel.id,sep='-'),spec=beast2.spec$compoundprior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)	
	dummy		<- newXMLNode("distribution", attrs= list(	idref=beast2.spec$treemodel.id), parent=tmp, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text="The origin prior.", parent=bxml.prior, doc=bxml, addFinalizer=T)	
	tmp			<- newXMLNode("prior", attrs= list(	id=paste("sprior",beast2.spec$bdsky.origin.id,sep='-'),
					name="distribution",
					x= paste('@',beast2.spec$bdsky.origin.id,sep='')), parent=bxml.prior, doc=bxml, addFinalizer=T)	
	dummy		<- hivc.beast2.get.prior(beast2.spec$bdsky.origin.prior, tmp, xmlAttrs(tmp)["id"], bxml)
	dummy		<- newXMLCommentNode(text="The serial samplingProb priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.sprop.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.sprop.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$bdsky.sprop.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})
	if(verbose)	cat(paste("\nadded serial samplingProb priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial R0 priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.R0.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.R0.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$bdsky.R0.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})												
	if(verbose)	cat(paste("\nadded serial R0 priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial becomingUninfectious priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.notInf.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.notInf.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- hivc.beast2.get.prior(beast2.spec$bdsky.notInf.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})	
	if(verbose)	cat(paste("\nadded serial becomingUninfectious priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="end: serial BDSKY priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	bxml.prior
}

######################################################################################
hivc.treeannotator.tiplabel2df<- function(x, beastlabel.idx.clu=1, beastlabel.idx.hivn=2, beastlabel.idx.hivd=3, beastlabel.idx.hivs=4, beastlabel.idx.samplecode=5, beastlabel.idx.rate=6)
{
	tmp		<- t( sapply(strsplit(x$tip.label,'_'), function(z)	z[c(beastlabel.idx.clu,beastlabel.idx.hivn, beastlabel.idx.hivd, beastlabel.idx.hivs, beastlabel.idx.samplecode, beastlabel.idx.rate)] ) )
	data.table(cluster=as.numeric(tmp[,1]), NegT= suppressWarnings(as.numeric(tmp[,2])), AnyPos_T1= as.numeric(tmp[,3]), TipT= as.numeric(tmp[,4]), FASTASampleCode=tmp[,5], rate=as.numeric(tmp[,6]), BEASTlabel=x$tip.label )				
}	
######################################################################################
hivc.treeannotator.get.clusterprob<- function(ph.beast, beastlabel.idx.clu=1, beastlabel.idx.samplecode=5, verbose=1)
{
	#	for each of the clusters, compute the posterior probability of a common MRCA
	clu.df	<- lapply(ph.beast, function(x)
			{							
				tmp		<- data.table(	tip=seq_len(Ntip(x)), 
						cluster=as.numeric( sapply( strsplit(x$tip.label,'_'),function(z)  z[beastlabel.idx.clu] ) ),
						FASTASampleCode=sapply( strsplit(x$tip.label,'_'),function(z)  z[beastlabel.idx.samplecode] )										
				)
				ans		<- merge( tmp[, list(node=hivc.clu.mrca(x, x.tip=tip)$mrca, FASTASampleCode=FASTASampleCode), by=cluster], subset(x$node.label, select=c(node, posterior)), by="node" )
				subset(ans, select=c(cluster, FASTASampleCode, posterior))
			})
	clu.df	<- rbindlist(clu.df)		
	if(verbose)	cat(paste("\nRange of posterior probabilities that each of the putative clusters each has a common MRCA, min=",min(clu.df[,posterior])," max=",max(clu.df[,posterior]) ))
	clu.df
}
######################################################################################
hivc.treeannotator.get.edgewidth<- function(ph, rates.df, scale.edgewidth= 12)
{
	if(is.null(rates.df))
		rates.df<- data.table(node=Nnode(ph,internal.only=0), rate=NA)
	
	edge.width	<- merge(rates.df, data.table(node=ph$edge[,2], edge=seq_len(nrow(ph$edge))), all.y=1, by="node")
	edge.width[, width:=edge.width[,rate]/mean(edge.width[,rate], na.rm=1)]
	set(edge.width, NULL, "width", (edge.width[,width]-1)*scale.edgewidth + 1)
	set(edge.width, which(edge.width[,is.na(width)]), "width", 1)
	set(edge.width, which(edge.width[,width<0.1]), "width", 0.1)
	setkey(edge.width, edge)
	edge.width
}
######################################################################################
#	ph<- cluphy; end.ctime=2013.3; cex.nodelabel=0.5; cex.tiplabel=0.5; file=NULL; pdf.width=7; pdf.height=20
hivc.treeannotator.plot<- function(ph, ph.root.ctime, youngest.tip.ctime, df.all, df.viro, df.immu, df.treatment=NULL, df.tstem=NULL, df.rates=NULL, end.ctime=2013.3, cex.nodelabel=0.5, cex.tiplabel=0.5, file=NULL, pdf.width=7, pdf.height=20)
{		
	require(RColorBrewer)
	if(class(file)=="character")
		pdf(file, width=pdf.width, height=pdf.height)
	par(mar=c(0,0,0,0))
	
	cols			<- brewer.pal(12,"Paired")
	ph.xlim			<- end.ctime-ph.root.ctime+ c(-22,6)
	ph.ylim			<- c(1,Ntip(ph)) + c(-1,1)			
	
	plot(ph, x.lim=ph.xlim, y.lim= ph.ylim, show.tip.label=0, edge.color = 0, tip.color = 0)
	# add calendar timeline			
	hivc.treeannotator.plot.ctimeline(ph, youngest.tip.ctime, end.ctime, add.yinch= 0.5)
	# add BEAST TMRCA 95% credibility interval
	ph.nodeheighthpds	<- hivc.treeannotator.nodelabels.getnodeheightHPD(ph, youngest.tip.ctime)
	#do not plot the BEAST TMRCA 95% credibility interval for those nodes for which more data than the 95% interval is available
	if(!is.null(df.tstem))	
		ph.nodeheighthpds	<- subset( ph.nodeheighthpds, !node%in%unique(df.tstem[,mrca]) )
	hivc.treeannotator.plot.hpdbars(ph, ph.nodeheighthpds, col=cols[1], lwd=4)
	# add NegT and AnyPos_T1
	ph.seronodeheight	<- hivc.treeannotator.sero.getnodeheight.range(ph, df.all, youngest.tip.ctime)
	hivc.treeannotator.plot.seronodeheightrange(ph, ph.seronodeheight, add.yinch= -0.03, width.yinch= 0.03, width.yinch.past.AnyPos_T1= 0, col=cols[2])
	# add lRNA timeline
	ph.viro.timeline	<- hivc.treeannotator.get.viro.timeline(ph, df.all, df.viro, youngest.tip.ctime, df.treatment=df.treatment)
	hivc.treeannotator.plot.viro.timeline(ph, ph.viro.timeline, viro.min= log10(300), width.yinch= 0.15, add.yinch= 0.005, col.bg= cols[c(5,10,12)], col.legend= cols[6], cex.txt= 0.2)
	# add BEAST posterior density of TMRCAs where available
	if(!is.null(df.tstem))
		hivc.treeannotator.plot.tipstem.timeline(ph, youngest.tip.ctime, df.tstem, width.yinch=0.1, add.yinch=0.005, col.bg=cols[1] )		
	# add CD4 timeline
	ph.immu.timeline	<- hivc.treeannotator.get.immu.timeline(ph, df.all, df.immu, youngest.tip.ctime, end.ctime=2013.3)
	hivc.treeannotator.plot.immu.timeline(ph, ph.immu.timeline, immu.min= 150, immu.max= 800, width.yinch= 0.15, add.yinch= -0.005, col.bg= cols[3], col.legend= cols[4], cex.txt= 0.2)
	# re-plot phylogeny
	ph$node.label		<- as.numeric(sapply( strsplit( ph$node.label, '_' ), function(x)	x[1] ))
	edge.width			<- hivc.treeannotator.get.edgewidth(ph, df.rates, scale.edgewidth= 8)
	hivc.phy.plotupon(ph, show.tip.label=0, show.node.label=1, cex=cex.nodelabel, edge.width=edge.width[,width],)
	# add rate labels			
	if(!is.null(df.rates))	
		hivc.treeannotator.plot.rates(ph, edge.width, add.xinch=-0.1, cex.rate=0.3)
	# add tip labels			
	ph.tiplabel			<- hivc.clu.get.tiplabels(ph, 	df.all, col.notmsm="#4EB3D3", col.Early="#EF9708", col.highVL="#FEE391", col.AfterTreat="#D4B9DA", col.green="#D9F0A3", col.latePres="#FA9FB5", select=c("CountryInfection","Trm","Sex","isAcute","lRNA.early","Patient","RegionHospital") )				
	tmp					<- rep( max(node.depth.edgelength(ph)) - (youngest.tip.ctime-ceiling(end.ctime)), Ntip(ph))
	hivc.clu.plot.tiplabels(seq_len(Ntip(ph)), ph.tiplabel$text, ph.tiplabel$col, xx=tmp, adj = c(-0.05, 0.5), cex=cex.tiplabel, add.xinch= 0.03, add.yinch= 0.02)
	# add legend	
	legend("topright", fill= cols[c(1,2,3,5,10,12)], legend=c("BEAST 95% TMRCA", "interval [last HIV-, diagnosis]", "CD4 timeline", "VL timeline", "VL timeline under treatment", "VL timeline under treatment"), bty='n', border=NA, cex=cex.tiplabel)
	
	if(class(file)=="character")
		dev.off()				
}
######################################################################################
hivc.treeannotator.plot.tipmarks<- function(ph, df.tips, add.xinch=0, add.yinch=0)
{
	tmp	<- match(df.tips[, FASTASampleCode], ph$tip.label)
	if(any(is.na(tmp)))		stop('unexpected missing tip names in ph for df.tips')				
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)	
	tmp		<- data.table(xx= lastPP$xx[ tmp ]+xinch(add.xinch), yy= lastPP$yy[ tmp ], tips=tmp)
	tmp		<- cbind(df.tips, tmp)
	points(tmp[,xx], tmp[,yy], col=tmp[,col], pch=tmp[,pch], cex=tmp[,cex] )	
}
######################################################################################
hivc.treeannotator.plot.rates<- function(ph, edge.width, add.xinch=-0.1, cex.rate=0.5)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)	
	tmp	<- edge.width[, list(xx= lastPP$xx[ ph$edge[edge,1] ]+xinch(add.xinch), yy= mean( lastPP$yy[ ph$edge[edge,] ] ), rate=rate), by="edge"]
	text(tmp[,xx], tmp[,yy], tmp[,rate], cex=cex.rate)
}
######################################################################################
hivc.treeannotator.plot.node.ctime<- function(df.node.ctime, ph.root.ctime, width.yinch=0.1, add.yinch=0.001, col.bg="black", density=30, lwd=0.5)		
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)	
	df.node.ctime[, yyi:=pdf]
	df.node.ctime[, xx:=q-ph.root.ctime]	
	setkey(df.node.ctime, node)
	df.node.ctime	<- merge(df.node.ctime, df.node.ctime[,	list(scale=yinch(width.yinch)/max(yyi))	,by="node"], by="node")
	set(df.node.ctime, NULL, "yyi",  df.node.ctime[,yyi*scale])
	df.node.ctime[, yy:= yinch(add.yinch)+lastPP$yy[df.node.ctime[,node]]]
	df.node.ctime[, col:=col.bg]
	setkey(df.node.ctime, node)		
	dummy<- sapply( unique(df.node.ctime[,node]), function(x)
			{
				z		<- df.node.ctime[J(x)]
				setkey(z, xx)
				polygon( c( z[,xx],z[nrow(z),xx] ), z[1,yy]+c( z[,yyi],0 ), border=z[1,col], col=z[1,col], density=density, lwd=lwd )  				
			})	
}
######################################################################################
hivc.treeannotator.plot.tipstem.timeline<- function(ph, youngest.tip.ctime, df.tstem, width.yinch=0.15, add.yinch=0, col.bg="grey75", density=30, lwd=0.5)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	set(df.tstem, NULL, "tstem", youngest.tip.ctime - df.tstem[,tstem])
	set(df.tstem, NULL, "tstem", max(node.depth.edgelength(ph)) - df.tstem[,tstem])
	df.tstem[, yyi:= density]
	setkey(df.tstem, tip)
	df.tstem<- merge(df.tstem, df.tstem[,	list(scale=yinch(width.yinch)/max(yyi))	,by="tip"], by="tip")	
	set(df.tstem, NULL, "yyi",  df.tstem[,yyi*scale])
	df.tstem[, yy:= yinch(add.yinch)+lastPP$yy[df.tstem[,tip]]]
	df.tstem[, col:=col.bg]
	setkey(df.tstem, tip)
	dummy<- sapply( unique(df.tstem[,tip]), function(x)
			{
				z		<- df.tstem[J(x)]
				setkey(z, tstem)
				polygon( c( z[,tstem],z[nrow(z),tstem] ), z[1,yy]+c( z[,yyi],0 ), border=z[1,col], col=z[1,col], density=density, lwd=lwd )  				
			})
}
######################################################################################
hivc.treeannotator.plot.ctimeline<- function(ph, youngest.tip.ctime, end.ctime, col.bg= c(my.fade.col("black",0.15),my.fade.col("black",0.05)), col.txt= c(my.fade.col("black",1),"transparent"), cex.txt= 0.5, add.yinch= 0.5)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)			
	tmp 	<- seq( max(lastPP$xx) - ( youngest.tip.ctime-floor(youngest.tip.ctime) )+1, lastPP$x.lim[1], -1 )
	df.time	<- data.table(ctime= seq( floor(youngest.tip.ctime), by=-1, len= length(tmp) ), xx.u=tmp)
	#	if tips don t reach to end.ctime, add more bars
	if(1+floor(youngest.tip.ctime)<floor(end.ctime))
	{
		tmp		<- seq( 1+floor(youngest.tip.ctime),floor(end.ctime),by=1 )
		tmp		<- data.table( ctime=rev(tmp), xx.u=rev(seq(df.time[1,xx.u]+1, len=length(tmp), by=1)))
		df.time	<- rbind(tmp, df.time)
	}
	df.time	<- cbind( 	df.time[-nrow(df.time), ], 
			data.table(	xx.l	= df.time[-1,xx.u], 
					col.bg	= rep(col.bg,ceiling(nrow(df.time)/length(col.bg)))[seq_len(nrow(df.time)-1)],
					col.txt = rep(col.txt,ceiling(nrow(df.time)/length(col.txt)))[seq_len(nrow(df.time)-1)]		) )
	
	rect(df.time[,xx.l], lastPP$y.lim[1]-yinch(add.yinch), df.time[,xx.u], lastPP$y.lim[2]+yinch(add.yinch), col=df.time[,col.bg], border=NA)				
	text(df.time[,xx.l+(xx.u-xx.l)/2.1], lastPP$y.lim[1]-yinch(add.yinch)/4, df.time[,ctime], cex=cex.txt, col=df.time[,col.txt], offset=0)
	text(df.time[,xx.l+(xx.u-xx.l)/2.1], lastPP$y.lim[2]+yinch(add.yinch)/4, df.time[,ctime], cex=cex.txt, col=df.time[,col.txt], offset=0)	
}
######################################################################################
hivc.treeannotator.plot.immu.timeline<- function(ph, ph.immu.timeline, immu.min= 150, immu.max= 800, immu.legend= c(200, 350, 500, immu.max), width.yinch= 0.15, add.yinch= -0.005, col.bg= cols[3], col.legend= cols[4], cex.txt= 0.2, lines.lwd=0.2)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	set(ph.immu.timeline, NULL, "PosCD4", max(node.depth.edgelength(ph)) - ph.immu.timeline[,PosCD4])
	ph.immu.timeline[, yyCD4:= CD4]
	set(ph.immu.timeline, which(ph.immu.timeline[,yyCD4]>immu.max), "yyCD4", immu.max)
	set(ph.immu.timeline, NULL, "yyCD4", ph.immu.timeline[,yyCD4]-immu.min)
	set(ph.immu.timeline, which(ph.immu.timeline[,yyCD4]<0), "yyCD4", 0.)			
	scale	<- yinch(width.yinch) / max( ph.immu.timeline[,yyCD4])
	set(ph.immu.timeline, NULL, "yyCD4", ph.immu.timeline[,yyCD4] * scale )
	ph.immu.timeline[, yy:= yinch(add.yinch)+lastPP$yy[ph.immu.timeline[,tip]]]
	
	dummy<- sapply( unique(ph.immu.timeline[,tip]), function(x)
			{
				z<- ph.immu.timeline[J(x)]
				polygon( c( z[,PosCD4], z[nrow(z),PosCD4], z[1,PosCD4] ), c( z[,yy-yyCD4], z[nrow(z),yy], z[1,yy] ), border=NA, col=col.bg	)
				sapply(z[1,yy]-(immu.legend-immu.min)*scale,function(i)		lines(z[c(1,nrow(z)),PosCD4], rep(i,2), col=col.legend, lty=3, lwd=lines.lwd)		)		
				text(rep(z[nrow(z),PosCD4],3),z[1,yy]-(immu.legend-immu.min)*scale,immu.legend,cex=cex.txt, col=col.legend)
			})
}
######################################################################################
hivc.treeannotator.get.immu.timeline<- function(ph, df, df.immu, youngest.tip.ctime, end.ctime=2013.3)
{				
	setkey(df, FASTASampleCode)
	tmp				<- df[J(ph$tip.label)][, list(tip=seq_along(ph$tip.label), FASTASampleCode=FASTASampleCode, Patient=Patient, DateDied=DateDied)]		
	ans				<- merge(subset(df.immu, select=c(Patient, PosCD4, CD4)), tmp, all.y=1, by="Patient")
	set(ans,NULL,"PosCD4",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,PosCD4]))
	set(ans,NULL,"DateDied",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,DateDied]))				
	set(ans,which(is.na(ans[,DateDied])), "DateDied", youngest.tip.ctime - end.ctime)				
	ans				<- ans[,list(PosCD4=c(PosCD4,DateDied[1]), CD4=c(CD4,tail(CD4,1))),by="tip"]
	setkey(ans,tip)
	ans
}
######################################################################################
#	viro.min= log10(300); width.yinch= 0.15; add.yinch= 0.005; col.bg= cols[c(5,9,10)]; col.legend= cols[6]; cex.txt= 0.2
hivc.treeannotator.plot.viro.timeline<- function(ph, ph.viro.timeline, viro.min= log10(300), viro.legend=c(3,4,5,6), width.yinch= 0.2, add.yinch= 0.005, col.bg= "red", col.legend="red", cex.txt= 0.2, lines.lwd=0.2)
{	
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	set(ph.viro.timeline, NULL, "PosRNA", max(node.depth.edgelength(ph)) - ph.viro.timeline[,PosRNA])
	ph.viro.timeline[, yylRNA:= lRNA]
	set(ph.viro.timeline, NULL, "yylRNA", ph.viro.timeline[,yylRNA]-viro.min)
	scale	<- yinch(width.yinch) / max( ph.viro.timeline[,yylRNA])
	set(ph.viro.timeline, NULL, "yylRNA", ph.viro.timeline[,yylRNA] * scale )
	ph.viro.timeline[, yy:= yinch(add.yinch)+lastPP$yy[ph.viro.timeline[,tip]]]
	#define treatment periods if there
	if(!any(ph.viro.timeline[, NoDrug!=0]))
		ph.viro.timeline[, col:=col.bg[1]]
	else
	{
		if(length(col.bg)==1)	stop("for TPeriod != NA, expect more than one 'col.bg'")
		col.bg.nNA	<- rep(col.bg[-1], ceiling(max(ph.viro.timeline[, TPeriod])/(length(col.bg)-1)))
		ph.viro.timeline[, col:=""]
		set(ph.viro.timeline, which(ph.viro.timeline[,NoDrug==0]), "col", col.bg[1])
		set(ph.viro.timeline, which(ph.viro.timeline[,NoDrug!=0]), "col", col.bg.nNA[ subset(ph.viro.timeline, NoDrug!=0)[,TPeriod] ])		
	}	
	setkey(ph.viro.timeline, tip)
	dummy<- sapply( unique(ph.viro.timeline[,tip]), function(x)
			{
				z		<- ph.viro.timeline[J(x)]	
				#reset TPeriod because there can be multiple off treatment periods (ie TPeriod==0) and we cannot lump them together in the next line
				#NOT NEEDED ANY LONGER	set(z, NULL, "TPeriod",cumsum(c(0,as.numeric(abs(diff(z[,TPeriod]))>0))))
				dummy	<- z[,	{								
							polygon( c( PosRNA, PosRNA[length(PosRNA)], PosRNA[1] ), c( yylRNA+yy, yy[length(yy)], yy[1] ), border=NA, col=col[1]	)	
						}, by="TPeriod"]								
				sapply(z[1,yy]+(viro.legend-viro.min)*scale,function(i)		lines(z[c(1,nrow(z)),PosRNA], rep(i,2), col=col.legend, lty=3, lwd=lines.lwd)		)		
				text(rep(z[nrow(z),PosRNA],3),z[1,yy]+(viro.legend-viro.min)*scale,paste("1e",viro.legend,sep=''),cex=cex.txt, col=col.legend)
				#stop()
			})
}
######################################################################################
hivc.treeannotator.get.viro.timeline<- function(ph, df, df.viro, youngest.tip.ctime, df.treatment=NULL, end.ctime=2013.3)
{
	setkey(df, FASTASampleCode)
	
	if(is.null(df.treatment))		#prepare a single viral load timeline without treatment periods
	{
		tmp				<- df[J(ph$tip.label)][, list(tip=seq_along(ph$tip.label), FASTASampleCode=FASTASampleCode, Patient=Patient, DateDied=DateDied)]		
		ans				<- merge(subset(df.viro, select=c(Patient, PosRNA, lRNA)), tmp, all.y=1, by="Patient")
		set(ans,NULL,"PosRNA",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,PosRNA]))
		set(ans,NULL,"DateDied",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,DateDied]))				
		set(ans,which(is.na(ans[,DateDied])), "DateDied", youngest.tip.ctime - end.ctime)				
		ans				<- ans[,list(PosRNA=c(PosRNA,DateDied[1]), lRNA=c(lRNA,tail(lRNA,1), TPeriod=0, NoDrug=0)),by="tip"]
		setkey(ans,tip)
	}
	else							#prepare a single viral load timeline with treatment periods
	{
		#as above except TPeriod=NA
		tmp				<- df[J(ph$tip.label)][, list(tip=seq_along(ph$tip.label), FASTASampleCode=FASTASampleCode, Patient=Patient, DateDied=DateDied)]		
		ans				<- merge(subset(df.viro, select=c(Patient, PosRNA, lRNA)), tmp, all.y=1, by="Patient")
		set(ans,NULL,"PosRNA",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,PosRNA]))
		set(ans,NULL,"DateDied",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,DateDied]))				
		set(ans,which(is.na(ans[,DateDied])), "DateDied", youngest.tip.ctime - end.ctime)				
		ans				<- ans[,list(Patient= rep(Patient[1], length(Patient)+1), PosRNA=c(PosRNA,DateDied[1]), lRNA=c(lRNA,tail(lRNA,1))),by="tip"]
		#prepare subset of df.treatment
		tmp				<- merge(subset(df.treatment, select=c(Patient, StartTime, StopTime, NoDrug)), unique(subset(ans,select=Patient)), all.y=1, by="Patient")	
		set(tmp,NULL,"StartTime",	youngest.tip.ctime - hivc.db.Date2numeric(tmp[,StartTime]))
		set(tmp,NULL,"StopTime",	youngest.tip.ctime - hivc.db.Date2numeric(tmp[,StopTime]))
		set(tmp, which( tmp[,StopTime==min(StopTime, na.rm=1)] ), "StopTime", youngest.tip.ctime - 2013.3)
		ans				<- merge(tmp, ans, by="Patient", allow.cartesian=1)
		tmp				<- which(ans[,is.na(StartTime)])
		set(ans,tmp,"StartTime",youngest.tip.ctime - 2013.3)
		set(ans,tmp,"StopTime",youngest.tip.ctime - 2013.3)
		set(ans,tmp,"NoDrug",0L)
		#now have all treatments and viral loads measures togethers
		ans					<- ans[,	{
					x		<- data.table(StartTime, StopTime, NoDrug, tip, PosRNA,  lRNA)											
					#select outside any treatment period
					tmp		<- subset(x, PosRNA>max(StartTime))								#handle no drug before first treatment
					tmp		<- rbind(tmp, subset(x, PosRNA<min(StopTime)) )					#handle no drug after stop treatment
					if(nrow(tmp))
					{
						setkey(tmp, PosRNA)
						tmp		<- unique(tmp)
						set(tmp,NULL,"NoDrug",0L)
						#select only those viral loads within a particular treatment period
						tmp		<- rbind(tmp, subset(x, StartTime>=PosRNA & PosRNA>=StopTime))												
					}
					else
						tmp		<- subset(x, StartTime>=PosRNA & PosRNA>=StopTime)
					#assign integer value to different treatment periods
					setkey(tmp, NoDrug, StartTime)
					tmp2	<- subset(unique(tmp), select=c(StartTime, NoDrug))[order(-StartTime, NoDrug)]
					tmp2	<- tmp2[, list(StartTime=StartTime, NoDrug=NoDrug, TPeriod=seq_along(NoDrug))]
					tmp		<- merge(tmp, tmp2, all.x=1, by=c("StartTime", "NoDrug"))											
					tmp		<- subset(tmp[order(-PosRNA)], select=c(PosRNA,  lRNA, NoDrug, TPeriod))
					#add endpoints for each treatment period
					tmp2	<- subset(tmp, TPeriod>min(TPeriod))[,  list(PosRNA=PosRNA[1], lRNA=lRNA[1] ),by=TPeriod]
					set(tmp2, NULL, "TPeriod", tmp2[,TPeriod]-1)
					tmp2	<- merge(tmp2, unique(subset(tmp, select=c(TPeriod, NoDrug))), by="TPeriod")
					tmp		<- rbind(tmp, subset(tmp2, select=c(PosRNA,  lRNA, NoDrug, TPeriod)))
					tmp
				},by="tip"]
		ans					<- ans[order(tip, -PosRNA, TPeriod)]		
	}
	ans
}
######################################################################################
hivc.treeannotator.sero.getnodeheight.range<- function(ph, df, youngest.tip.ctime)
{
	setkey(df, FASTASampleCode)
	ans					<- cbind( data.table(tip=seq_along(ph$tip.label)), subset(df[J(ph$tip.label)], select=c(NegT, AnyPos_T1, DateDied)) )
	set(ans,NULL,"NegT",		youngest.tip.ctime - hivc.db.Date2numeric(ans[,NegT]))
	set(ans,NULL,"AnyPos_T1",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,AnyPos_T1]))				
	set(ans,NULL,"DateDied",	youngest.tip.ctime - hivc.db.Date2numeric(ans[,DateDied]))
	set(ans,which(is.na(ans[,DateDied])), "DateDied", 0.)
	ans
}
######################################################################################
hivc.treeannotator.plot.seronodeheightrange<- function(ph, ph.seronodeheight, add.yinch= -0.05, width.yinch= 0.1, width.yinch.past.AnyPos_T1= 0.02, col="red")
{						
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	if (!lastPP$use.edge.length) 
		stop("function needs edge length information")
	if (lastPP$type != "phylogram") 
		stop("currently only 'type == phylogram' supported")
	if (lastPP$dir!="rightwards")
		stop("currently only rightwards supported")			
	
	
	tmp				<- max(node.depth.edgelength(ph)) - subset(ph.seronodeheight, select=c(NegT, AnyPos_T1, DateDied))		#from node heights to root heights, assuming root is plotted at 0
	yy.l 			<- lastPP$yy[ph.seronodeheight[,tip]] + yinch(add.yinch)			
	rect(tmp[,NegT], yy.l, tmp[,AnyPos_T1], yy.l+yinch(width.yinch), col = col, border=NA)			
	rect(tmp[,AnyPos_T1], yy.l, tmp[,DateDied], yy.l+yinch(width.yinch.past.AnyPos_T1), col=col, border=NA)	
}
######################################################################################
hivc.treeannotator.nodelabels.getnodeheightHPD<- function(ph, youngest.tip.ctime, node.label.hpd.l= 3, node.label.hpd.u= 4)
{
	nodes			<- which(!is.na(ph$node.label))
	node.hpd		<- t( sapply(strsplit(ph$node.label[nodes],'_'), function(x)	as.numeric(x[c(node.label.hpd.l, node.label.hpd.u)])) )
	node.hpd		<- youngest.tip.ctime - node.hpd				#from calendar time to node heights
	data.table(node=nodes, hpd.l=node.hpd[,1], hpd.u=node.hpd[,2])
}
######################################################################################
hivc.treeannotator.plot.hpdbars<- function(ph, ph.nodeheighthpds, col="grey75", lwd=4 ) 
{
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	if (!lastPP$use.edge.length) 
		stop("function needs edge length information")
	if (lastPP$type != "phylogram") 
		stop("currently only 'type == phylogram' supported")
	if (lastPP$dir!="rightwards")
		stop("currently only rightwards supported")				
	tmp				<- max(node.depth.edgelength(ph)) - subset(ph.nodeheighthpds, select=c(hpd.l,hpd.u))		#from node heights to root heights, assuming root is plotted at 0
	yy 				<- lastPP$yy[Ntip(ph)+ph.nodeheighthpds[,node]]
	segments(tmp[,hpd.l], yy, tmp[,hpd.u], yy, col = col, lwd = lwd)					
}	
######################################################################################
hivc.treeannotator.get.rates<- function(ph, tip.df, nodelabel.idx.edgewidth=5)
{
	tmp				<- as.numeric( sapply( strsplit(ph$node.label,'_'),function(x)	x[nodelabel.idx.edgewidth] ) )	
	rates.df		<- data.table(node=seq_len(Nnode(ph))+Ntip(ph), rate=tmp)
	ans				<- rbind(rates.df, subset(tip.df, select=c(node, rate)))
	#set(ans, which(is.na(ans[,rate])),"rate",mean(ans[,rate],na.rm=1))
	setkey(ans,node)
	ans
}
######################################################################################
hivc.treeannotator.get.phy<- function(ph.beast, beastlabel.idx.clu=1, beastlabel.idx.hivs=4, beastlabel.idx.samplecode=5, beastlabel.idx.rate=6, verbose=1, debug=0)
{
	#	get root height for final tree in calendar time
	ph.tip.ctime	<- sapply(ph.beast, function(x) max( as.numeric( sapply(strsplit(x$tip.label,'_'), function(x)	x[beastlabel.idx.hivs] ) ) ))			
	ph.root.ctime	<- min( sapply(seq_along(ph.beast), function(i)	ph.tip.ctime[i]-max(node.depth.edgelength(ph.beast[[i]]))	) )
	#	extract tip label information
	tip.df			<- lapply(ph.beast, function(x) hivc.treeannotator.tiplabel2df(x, beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate) )
	#	prepare cluster subtrees
	clu.subtrees	<- lapply( seq_along(ph.beast), function(i)
			{				
				x<- ph.beast[[i]]
				#	convert heights into calendar time and collapse node.label
				tmp					<- ph.tip.ctime[i]								
				#subset(x$node.label, node==147, select=c(node, height_median, height_95_HPD_MIN, height_95_HPD_MAX, posterior))
				tmp					<- x$node.label[, list(	node=node, 
								height_median=ph.tip.ctime[i]-height_median, 
								height_95_HPD_MIN=ph.tip.ctime[i]-height_95_HPD_MAX, 
								height_95_HPD_MAX=ph.tip.ctime[i]-height_95_HPD_MIN, 
								rate_median=rate_median,
								posterior=posterior)]			
				tmp					<- tmp[, list(node.label= paste(posterior,height_median, height_95_HPD_MIN, height_95_HPD_MAX, rate_median, sep='_')),by="node"]
				x$node.label		<- tmp[, node.label]
				x$node.label.format	<- "posterior height_median height_95_HPD_MIN height_95_HPD_MAX rate_median"				
				#
				#	extract rooted ExaML clusters
				#								
				x$tip.label	<- tip.df[[i]][, FASTASampleCode]
				tmp			<- tip.df[[i]][, list(node=getMRCA(x,FASTASampleCode)),by=cluster]
				clu.subtrees<- lapply(tmp[,node], function(z)
						{	
							ans						<- extract.clade(x, z, root.edge= 1, interactive = FALSE)
							ans$root.edge			<- as.numeric(strsplit(ans$node.label[1],'_')[[1]][2])-ph.root.ctime		#reset root edge against root of all runs combined
							ans$node.label.format	<- x$node.label.format
							ans
						})	
				names(clu.subtrees)<- tmp[,cluster]	
				clu.subtrees
			})
	clu.subtrees	<- eval(parse(text= paste("c(",paste('clu.subtrees[[',seq_along(clu.subtrees),']]', sep='',collapse=','),")",sep='') ))
	if(verbose)	cat(paste("\nFound ExaML clusters in treeannotator files, number of clusters is n=", length(clu.subtrees) ))
	if(debug)
		clu.subtrees	<- lapply(1:3, function(i) clu.subtrees[[i]] )
	#	join all clusters 
	cluphy				<- eval(parse(text=paste('clu.subtrees[[',seq_along(clu.subtrees),']]', sep='',collapse='+')))	
	if(verbose)	cat(paste("\nFound ExaML clusters in treeannotator files, number of sequences is n=", Ntip(cluphy) ))
	#	retain tip info for those tip labels in cluphy
	tip.df			<- rbindlist(tip.df)	
	tip.df			<- merge(data.table(FASTASampleCode=cluphy$tip.label), tip.df, by="FASTASampleCode")
	tip.df			<- cbind(tip.df, node=seq_len(Ntip(cluphy)))
	
	list(cluphy=cluphy, ph.tip.df=tip.df, ph.tip.ctime=ph.tip.ctime, ph.root.ctime=ph.root.ctime)
}
######################################################################################
hivc.treeannotator.get.tmrcas<- function(ph.beast, beastlabel.idx.hivs=4)
{
	#	for each tree, return 	height_median height_95_HPD_MIN height_95_HPD_MAX for the parent of each tip 		
	ph.trmca	<- lapply(ph.beast, function(x)
			{
				#select heights and convert heights into calendar time
				tip.latest			<- max( as.numeric( sapply(strsplit(x$tip.label,'_'), function(x)	x[beastlabel.idx.hivs] ) ) )
				tip.parents			<- unique( x$edge[x$edge[,2]<=Ntip(x),1] )		
				ans					<- x$node.label[J(tip.parents)][, list(node=node, height_median=tip.latest-height_median, height_95_HPD_MIN=tip.latest-height_95_HPD_MAX, height_95_HPD_MAX=tip.latest-height_95_HPD_MIN)]
				tmp					<- sapply(ans[,node], function(z)	x$edge[ x$edge[,1]==z, 2 ] )
				tmp[tmp>Ntip(x)]	<- NA
				tmp					<- t( apply(tmp,2,function(z) x$tip.label[ sort(z,na.last=1) ]) )
				ans[,height_95_diff:= height_95_HPD_MAX-height_95_HPD_MIN]
				ans[,tip1:= tmp[,1]]		
				ans[,tip2:= tmp[,2]]					
				ans
			})
	ph.trmca	<- rbindlist( ph.trmca )
}