seattle.170601.rm.drug.resistance.mutations<- function()
{
	require(big.phylo)
	#	handle DRMs
	infile.fasta		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.fasta'
	outfile				<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm.fasta'
	sq					<- read.dna(infile.fasta, format='fa')	
	tmp					<- which(grepl("HXB2",rownames(sq)))
	rownames(sq)[tmp]	<- 'HXB2'
	tmp					<- big.phylo:::seq.rm.drugresistance(sq)
	sq					<- tmp$nodr.seq
	write.dna(sq, file= outfile, format='fasta')
}

seattle.170601.fastree<- function()
{	
	require(big.phylo)
	#	first tree
	infile.fasta	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.alignment.fasta'
	outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
	tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
	#	second tree, DRMs removed + including references for rooting
	infile.fasta	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm.fasta'
	outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
	tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
	#	run on command line
	cat(tmp)		
}

seattle.170601.clustering<- function()
{	
	require(big.phylo)
	#	first tree
	infile.ft		<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm_ft.newick'
	ph				<- read.tree(infile.ft)
	#	reroot at AY371169
	tmp				<- which(grepl('AY371169',ph$tip.label))
	ph				<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
	ph 				<- ladderize( ph )
	#	process bootstrap support values	
	tmp				<- as.numeric( ph$node.label )		
	tmp[is.na(tmp)]	<- 0
	ph$node.label	<- tmp
	#	
	write.tree(ph, gsub('\\.newick','_reroot.newick',infile.ft))
	pdf(gsub('newick','pdf',infile.ft), w=15, h=350)
	plot(ph, cex=0.4, show.node.label=TRUE)
	dev.off()
	#
	#read patristic distances
	require(hivclust)
	stat.fun						<- hivc.clu.min.transmission.cascade	
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)	
	thresh.brl						<- 0.045
	thresh.bs						<- 0.8
	#produce clustering 
	clustering						<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph$node.label, retval="all")	
	outfile							<- gsub('\\.newick',paste0('_reroot.clu',thresh.brl*1e2,'-',thresh.bs*1e2,'.pdf'),infile.ft)	
	hivc.clu.plot(ph, 
					clustering[["clu.mem"]], 
					file=outfile, 
					pdf.off=0, pdf.width=15, pdf.height=1000, 
					show.node.label=TRUE, show.edge.label=TRUE,  
					cex.edge.incluster=2)
	dev.off()
	
	dfc		<- data.table(TAXA=ph$tip.label, CLU=clustering$clu.mem[1:Ntip(ph)])
	outfile	<- gsub('\\.newick',paste0('_reroot.clu',thresh.brl*1e2,'-',thresh.bs*1e2,'.rda'),infile.ft)
	save(ph, clustering, dfc, file, file=outfile)	
}

seattle.170601.clustering.plot<- function()
{
	
	require(data.table)
	infile	<- "/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2017/2017_Seattle/Seattle.June2017.refse.ndrm_ft_reroot.clu4.5-80.rda"
	load(infile)
	#	exclude reference sequences non-Seattle
	dfc		<- subset(dfc, !grepl('AY371169|HXB2',TAXA))
	#	read meta-data 
	dfc[, ID:= sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',1)]
	dfc[, SXO:= sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',2)]
	dfc[, RACE:= sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',3)]
	dfc[, SEQ_YR:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',4))]
	dfc[, SEQ_Q:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',5))]
	dfc[, HIV_YR:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',6))]
	dfc[, BIRTH_YR:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',7))]
	dfc[, AGE_AT_SEQ:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',8))]
	dfc[, AGE_AT_HIV:= as.numeric(sapply(strsplit(dfc[, TAXA],'|',fixed=TRUE),'[[',9))]
	#	fixup
	set(dfc, dfc[, which(SXO=='UNKNOWN')], 'SXO', NA_character_)
	set(dfc, dfc[, which(SXO=='')], 'SXO', NA_character_)
	set(dfc, dfc[, which(SXO=='BLOOD')], 'SXO', 'Other')
	set(dfc, dfc[, which(SXO=='PERINAT')], 'SXO', 'Other')
	set(dfc, dfc[, which(RACE=='FALSE')], 'RACE', NA_character_)
	#	give each singleton a new cluster ID (cluster of size 1)
	tmp	<- dfc[, which(is.na(CLU))]
	set(dfc, tmp, 'CLU', max(dfc$CLU, na.rm=TRUE)+seq_along(tmp))
	
	#	compute cluster statistics 
	tmp	<- dfc[, list(	CLU_TS= min(HIV_YR, na.rm=TRUE),
				CLU_TE= max(HIV_YR, na.rm=TRUE),
				CLU_TM=mean(HIV_YR, na.rm=TRUE),
				CLU_N=length(ID),
				CLU_N_MSM=length(which(grepl('MSM',SXO))),
				CLU_N_IDU=length(which(grepl('IDU',SXO))),
				CLU_N_HET=length(which(grepl('HET',SXO)))
				), by='CLU']
	dfc	<- merge(dfc, tmp, by='CLU') 
	#	determine max lkl SXO of cluster
	tmp	<- melt(dfc, id.vars='CLU',measure.vars=c('CLU_N_MSM','CLU_N_IDU','CLU_N_HET'))
	tmp	<- tmp[, list(CLU_MX_SXO=variable[which.max(value)]), by='CLU']
	set(tmp, NULL, 'CLU_MX_SXO', tmp[, gsub('CLU_N_','',CLU_MX_SXO)])
	dfc	<- merge(dfc, tmp, by='CLU')
	#	set factors 
	set(dfc, NULL, 'CLU_MX_SXO', dfc[, factor(CLU_MX_SXO, levels=c('MSM','HET','IDU'), labels=c('MSM','heterosexual','injecting\ndrug\nuser'))])
	set(dfc, dfc[, which(is.na(SXO))], 'SXO', 'Unknown')
	set(dfc, dfc[, which(SXO=='HETEROp')], 'SXO', 'HETERO')
	set(dfc, NULL, 'SXO', dfc[, factor(SXO,  levels=c('MSM','MSM_IDU','IDU','HETERO','Other','Unknown'), 
											 labels=c('MSM','MSM/injecting drug user','injecting drug user','heterosexual','other','unknown'))])
	#
	#	plot clusters only
	dfp	<- subset(dfc, CLU_N>1 & !is.na(CLU_TS))
	#	order by CLU_MX_SXO and start date of cluster
	tmp	<- unique(subset(dfp, select=c(CLU, CLU_MX_SXO, CLU_TS)))
	tmp	<- tmp[order(CLU_MX_SXO, CLU_TS), ]	
	set(tmp, NULL, 'CLU_ORDERED', tmp[, factor(CLU, levels=tmp$CLU)])
	dfp	<- merge(dfp, subset(tmp, select=c(CLU, CLU_ORDERED)), by='CLU')
	tmp	<- melt(unique(subset(dfp, select=c(CLU_ORDERED, CLU_TS, CLU_TE, CLU_MX_SXO))), id.vars=c('CLU_ORDERED','CLU_MX_SXO'))
	ggplot(dfp, aes(y=CLU_ORDERED)) +
			scale_colour_brewer(palette='Set1') +
			scale_x_continuous(breaks=seq(1950,2030,5), expand=c(0,0)) +
			geom_line(data=tmp, aes(x=value, y=CLU_ORDERED), colour='black') +
			geom_point(aes(x=HIV_YR, colour=SXO)) +
			facet_grid(CLU_MX_SXO~., scales='free_y', space='free_y') +
			labs(x='\nyear of diagnosis', y='phylogenetic clusters\n', colour='sexual orientation') +
			theme_bw() +
			theme(	strip.text.y=element_text(angle=0),
					panel.grid.minor.x=element_line(colour="grey70", size=0.25), 
					panel.grid.major.x=element_line(colour="grey70", size=0.5),
					panel.grid.minor.y=element_blank(), 
					panel.grid.major.y=element_blank(),					
					axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					legend.position='bottom')
	ggsave(gsub('\\.rda','_cluevolution_by_sxo.pdf',infile), w=12,h=18)	
}