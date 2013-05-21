project.bezemer2013a.figs<- function()
{
	require(data.table)
	require(RColorBrewer)
	
	dir.name		<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/derived"
	f.name			<- "Bezemer2014_clusters.csv"
	df				<- read.csv(paste(dir.name,f.name,sep='/'), stringsAsFactors=FALSE)
	
	#deal with missing data and get into reasonable shape
	df[df[,"fpos1"]=="","fpos1"]			<- NA
	df[,"fpos1"]							<- as.Date(df[,"fpos1"], format="%d/%m/%Y")	
	df[,"patient"]							<- factor(df[,"patient"])	
	df[which(df[,"Country"]==""),"Country"]	<- NA
	df[,"Country"]							<- factor(df[,"Country"])
	df[which(df[,"hregion"]==""),"hregion"]	<- "unknown"
	df[,"hregion"]							<- factor(df[,"hregion"], levels=c("CU","Ea","NH","No","So","UT","ZH","unknown"))
	df[which(df[,"sexy"]==""),"sexy"]		<- NA
	df[,"sexy"]								<- factor(df[,"sexy"])
	df[which(df[,"rout"]==""),"rout"]		<- "unknown"
	df[which(df[,"rout"]=="U"),"rout"]		<- "unknown"
	df[,"rout"]								<- factor(df[,"rout"])
	#print(levels(df[,"Country"]))
	#print(levels(df[,"hregion"]))
	#print(levels(df[,"sexy"]))
	#print(levels(df[,"rout"]))
	#print(summary(df[,"fpos1"]))
	df	<- data.table(df, key="clustername")
	str(df)
	#I wondered if you'd like to make a figure by date of fpos1 and then A. coloured by rout, and B. coloured by hregion.				
	#The los alamos sequences do not have this information so you could colour by country and date of sequence (yearLA) instead, and in a second option leave them out.	
	
	dir.name			<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tmp"
	CLUSTER.N.THRESHOLD	<- 10
	ylab				<- "transmission network ID"
	yheight				<- 1.5
	
	if(1)	#los alamos left out, figure x: fpos1 y: clusters, coloured by hregion
	{
		#input
		f.name						<- "BezemerFig1_hregion_noSingletons.pdf"
		df[,covariate:=hregion]		#set which covariate to use
		tmp							<- df[,covariate]
		tmp[ which(tmp=="unknown") ]<- NA
		df[,covariateNA:=tmp]
		
		#create plot
		xlim			<- range( df[,fpos1], na.rm=1 )
		cluster.n		<- df[,length(fpos1),by=clustername]		
		#select clusters with more than CLUSTER.N.THRESHOLD members
		df				<- df[ subset(cluster.n,V1>=CLUSTER.N.THRESHOLD,clustername) ]
		#sort clusters by median covariate and then by median fpos	
		clusters.sortby1<- as.vector( df[,median(as.numeric(covariateNA), na.rm=1), by=clustername]$V1 )
		clusters.sortby2<- as.vector( df[,median(fpos1, na.rm=1), by=clustername]$V1 )
		clusters.sortby	<- order(clusters.sortby1,clusters.sortby2)
		clusters		<- unique( df[,clustername] )[clusters.sortby]				
		clusters.levels	<- levels(df[,covariate])			
		clusters		<- lapply(clusters,function(x)
				{
					subset(df,clustername==x,c(fpos1,covariate))
				})
		ylim			<- c(0,length(clusters))
		cols			<- c(brewer.pal(length(clusters.levels)-1, "Set3"), "grey30")		#add grey for "unknown"
		cols			<- sapply(cols, function(x) my.fade.col(x,0.7))
		pch				<- c(rep(19,length(clusters.levels)),3)		
		xlim[1]			<- xlim[1]-700
		file			<- paste(dir.name,f.name,sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=yheight*7)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xaxt='n',ylab=ylab, xlab="")
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		sapply(seq_along(clusters),function(i)
				{							
					points( as.numeric(clusters[[i]][,fpos1]), rep(i,nrow(clusters[[i]])), col=cols[ as.numeric(clusters[[i]][,covariate]) ], pch=pch[1], cex=0.75 )					
					points( median(as.numeric(clusters[[i]][,fpos1]), na.rm=1), i, col="black", pch=pch[length(pch)], cex=0.5 )
				})
		pch[pch==19]	<- 21
		col				<- c(rep("transparent",length(clusters.levels)),"black")
		legend("topleft",bty='n',pt.bg=c(cols,"black"),pch=pch,legend=c(levels(df[,covariate]),"median 1st HIV+"), col=col)
		dev.off()		
	}
	if(1)	#los alamos left out, figure x: fpos1 y: clusters, coloured by rout
	{
		#input
		f.name				<- "BezemerFig1_rout_noSingletons.pdf"
		df[,covariate:=rout]		#set which covariate to use
		
		#create plot
		xlim			<- range( df[,fpos1], na.rm=1 )
		cluster.n		<- df[,length(fpos1),by=clustername]		
		#select clusters with more than CLUSTER.N.THRESHOLD members
		df				<- df[ subset(cluster.n,V1>=CLUSTER.N.THRESHOLD,clustername) ]
		#sort clusters by median fpos	
		clusters.sortby	<- sort( as.vector( df[,median(fpos1, na.rm=1), by=clustername]$V1 ), index.return=T)$ix		
		clusters		<- unique( df[,clustername] )[clusters.sortby]				
		clusters.levels	<- levels(df[,covariate])			
		clusters		<- lapply(clusters,function(x)
				{
					subset(df,clustername==x,c(fpos1,covariate))
				})
		ylim			<- c(0,length(clusters))
		cols			<- c(brewer.pal(length(clusters.levels)-1, "Set3"), "grey30")		#add grey for "unknown"
		cols			<- sapply(cols, function(x) my.fade.col(x,0.7))
		pch				<- c(rep(19,length(clusters.levels)),3)		
		xlim[1]			<- xlim[1]-700
		file			<- paste(dir.name,f.name,sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=yheight*7)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xaxt='n',ylab=ylab, xlab="")
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		sapply(seq_along(clusters)[-1],function(i)
				{							
					points( as.numeric(clusters[[i]][,fpos1]), rep(i,nrow(clusters[[i]])), col=cols[ as.numeric(clusters[[i]][,covariate]) ], pch=pch[1], cex=0.75 )					
					points( median(as.numeric(clusters[[i]][,fpos1]), na.rm=1), i, col="black", pch=pch[length(pch)], cex=0.5 )
				})
		pch[pch==19]	<- 21
		col				<- c(rep("transparent",length(clusters.levels)),"black")
		legend("topleft",bty='n',pt.bg=c(cols,"black"),pch=pch,legend=c(levels(df[,covariate]),"median 1st HIV+"), col=col)
		dev.off()
		#save
		
	}
	if(0)
	{
		#los alamos left out, figure x: fpos1 y: clusters, coloured by rout
	}
	
}