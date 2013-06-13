project.bezemer2013a.figs.v130528<- function()
{
	require(data.table)
	require(RColorBrewer)	 
	
	dir.name		<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"	
	signat			<- "130528"
	f.name			<- paste("Bezemer2014_clusters_",signat,".csv",sep='')
	df				<- read.csv(paste(dir.name,f.name,sep='/'), stringsAsFactors=FALSE)
	
	#deal with missing data and get into reasonable shape
	#@TODO Daniela, see if I interpret all the fields correctly, especially what is NA
	df[df[,"fpos1"]=="","fpos1"]				<- NA
	df[,"fpos1"]								<- as.Date(df[,"fpos1"], format="%d/%m/%Y")	
	df[,"patient"]								<- factor(df[,"patient"])
	tmp											<- !is.na(df[,"RegionOrigin"]) & (df[,"RegionOrigin"]=="" | df[,"RegionOrigin"]=="U") 		
	df[tmp,"RegionOrigin"]						<- NA
	df[,"RegionOrigin"]							<- factor(df[,"RegionOrigin"])
	df[which(df[,"Country"]==""),"Country"]		<- NA
	df[,"Country"]								<- factor(df[,"Country"])
	df[which(df[,"hregion"]==""),"hregion"]		<- "unknown"
	df[,"hregion"]								<- factor(df[,"hregion"], levels=c("CU","Ea","NH","No","So","UT","ZH","unknown"))
	df[which(df[,"sexy"]==""),"sexy"]			<- NA
	df[,"sexy"]									<- factor(df[,"sexy"])
	df[which(df[,"rout"]==""),"rout"]			<- "unknown"
	df[which(df[,"rout"]=="U"),"rout"]			<- "unknown"
	df[,"rout"]									<- factor(df[,"rout"])
	#print(levels(df[,"Country"]))
	#print(levels(df[,"hregion"]))
	#print(levels(df[,"sexy"]))
	#print(levels(df[,"rout"]))
	#print(summary(df[,"fpos1"]))
	df	<- data.table(df, key="clustername")
	str(df)
	
	dir.name			<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"
	CLUSTER.N.THRESHOLD	<- 10
	ylab				<- "transmission network ID"
	yheight				<- 1.5
	cex.points			<- 0.75
	#networks with >= 10 sequences by date of diagnosis, ordered by duration of the network (added)
	#bottom line: singletons and smaller clusters as two summary networks on the bottom line.
	 
	#C, by country origin, blue NL, orange NLseAnt, green LAtin, black Europe. 
	#
	if(0)
	{
		#input
		f.name					<- paste("BezemerFig1_RegionOrigin_wSingletons_",signat,".pdf",sep='')
		#reduce risk group as desired and set as covariate
		tmp													<- unclass(df[,RegionOrigin]	)
		code												<- seq_along(attr(tmp,"levels"))
		names(code)											<- attr(tmp,"levels")
		print(code)
		tmp[ tmp==code[c("NL")] ]							<- 20
		tmp[ tmp==code[c("NLseAnt")] ]						<- 21
		tmp[ tmp==code[c("Lat")] ]							<- 22
		tmp[ tmp%in%code[c("EUC","EUO","EUW")] ]			<- 23
		tmp[ tmp%in%code[c("AUS","NAM","OAP","SSA","ZAz")] ]<- 24
		tmp													<- factor(tmp, levels=as.character(20:24), labels=c("NL","NLseAnt","Latin","Europe","other"))
		df[,covariate:=tmp ]
		df[,covariateNA:=tmp]
		
		
		print( attr(unclass(df[,RegionOrigin]),"levels") )
		print( attr(unclass(df[,covariate]),"levels") )
		
		xlim			<- range( df[,fpos1], na.rm=1 )
		zlim			<- range( df[,duration], na.rm=1 )
		#select singletons
		singletons		<- subset(df,clustername==0,c(fpos1,covariate,duration))
		#rm singletons from data table
		df				<- subset(df,clustername!=0,)
		#select clusters with less than CLUSTER.N.THRESHOLD non-NA fpos1 members
		cluster.n		<- df[,length(which(!is.na(fpos1))),by=clustername]		
		clusters.small	<- df[ subset(cluster.n,V1<CLUSTER.N.THRESHOLD,clustername) ]		
		#select clusters with more than CLUSTER.N.THRESHOLD non-NA fpos1 members
		df				<- df[ subset(cluster.n,V1>=CLUSTER.N.THRESHOLD,clustername) ]
		#sort clusters by duration -- this is the same as median duration
		
		clusters.sortby1<- as.vector( df[,as.numeric(any(hregion=="CU", na.rm=1)), by=clustername]$V1 )
		clusters.sortby2<- as.vector( df[,median(duration, na.rm=1), by=clustername]$V1 )
		clusters.sortby	<- order(clusters.sortby1,clusters.sortby2, decreasing=T)		
		#extract clusters one by one
		clusters		<- unique( df[,clustername] )[clusters.sortby]		
		clusters.levels	<- levels(df[,covariate])			
		clusters		<- lapply(clusters,function(x)
				{
					subset(df,clustername==x,c(fpos1,covariate,duration,hregion))
				})
		ylim			<- c(-15,length(clusters))
		cols			<- c(brewer.pal(length(clusters.levels)-1, "Set1"), "grey30")[c(2,1,3,4,5)]		#add grey for "unknown"
		#cols			<- c("green","red","blue","grey30")
		cols			<- sapply(cols, function(x) my.fade.col(x,0.7))
		pch				<- c(rep(19,length(clusters.levels)))
		xlim[1]			<- xlim[1]-700
		file			<- paste(dir.name,f.name,sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=yheight*7)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',xlab="")
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		#axis(3, )
		sapply(seq_along(clusters),function(i)
				{							
					cluster.cex	<- (1 - nrow(subset(clusters[[i]],hregion=="unknown")) / nrow(clusters[[i]])) * cex.points
					cluster.ix	<- order(as.numeric(clusters[[i]][,fpos1]))
					cluster.x	<- as.numeric(clusters[[i]][,fpos1])[cluster.ix]										
					cluster.y	<- rep(i,nrow(clusters[[i]]))
					cluster.z	<- as.numeric(clusters[[i]][,covariate])[cluster.ix]
					#print(clusters[[i]][,covariate])
					#if(cluster.z[1]!=4)
					lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
					points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=pch[1], cex=cluster.cex )										
				})		
		#plot small clusters
		abline(h=-5, lty=3,col="grey50")
		points( as.numeric(clusters.small[,fpos1]), rep(-5,nrow(clusters.small)), col=cols[ as.numeric(clusters.small[,covariate]) ], pch=pch[1], cex=0.75 )
		#plot singletons
		abline(h=-10, lty=4,col="grey50")
		points( as.numeric(singletons[,fpos1]), rep(-10,nrow(singletons)), col=cols[ as.numeric(singletons[,covariate]) ], pch=pch[1], cex=0.75 )
		
		pch[pch==19]	<- 21
		col				<- c(rep("transparent",length(clusters.levels)),"black")
		legend("topleft",bty='n',pt.bg=cols,pch=pch,legend=c(levels(df[,covariate])), col=col)
		legend("bottomleft", bty='n', lty=c(4), col="grey50", legend=c("singletons"))
		legend("bottomright", bty='n', lty=c(3), col="grey50", legend=c("networks<10"))
		dev.off()			
	}
	#B, in blue hregion 'NH', in orange in 'CU', in green in other provinces.
	if(0)	#los alamos left out, figure x: fpos1 y: clusters, coloured by hregion
	{
		#input
		f.name						<- paste("BezemerFig1_hregion_wSingletons_",signat,".pdf",sep='')
		#reduce risk group as desired and set as covariate
		tmp												<- unclass(df[,hregion]	)
		code											<- seq_along(attr(tmp,"levels"))
		names(code)										<- attr(tmp,"levels")						
		tmp[ tmp%in%code[c("Ea","No","So","UT","ZH")] ]	<- 12
		tmp[ tmp==code[c("NH")] ]						<- 10
		tmp[ tmp==code[c("CU")] ]						<- 11
		tmp[ tmp==code[c("unknown")] ]					<- 13	
		tmp												<- factor(tmp, levels=as.character(10:13), labels=c("NH","CU","other","unknown"))
		df[,covariate:=tmp ]
		tmp[ which(tmp=="unknown") ]					<- NA
		df[,covariateNA:=tmp]
		#print(df)
		
		xlim			<- range( df[,fpos1], na.rm=1 )
		zlim			<- range( df[,duration], na.rm=1 )
		#select singletons
		singletons		<- subset(df,clustername==0,c(fpos1,covariate,duration))
		#rm singletons from data table
		df				<- subset(df,clustername!=0,)
		#select clusters with less than CLUSTER.N.THRESHOLD non-NA fpos1 members
		cluster.n		<- df[,length(which(!is.na(fpos1))),by=clustername]		
		clusters.small	<- df[ subset(cluster.n,V1<CLUSTER.N.THRESHOLD,clustername) ]		
		#select clusters with more than CLUSTER.N.THRESHOLD non-NA fpos1 members
		df				<- df[ subset(cluster.n,V1>=CLUSTER.N.THRESHOLD,clustername) ]
		#sort clusters by duration -- this is the same as median duration
		
		clusters.sortby1<- as.vector( df[,as.numeric(any(covariateNA=="CU", na.rm=1)), by=clustername]$V1 )
		clusters.sortby2<- as.vector( df[,median(duration, na.rm=1), by=clustername]$V1 )
		clusters.sortby	<- order(clusters.sortby1,clusters.sortby2, decreasing=T)		
		#extract clusters one by one
		clusters		<- unique( df[,clustername] )[clusters.sortby]		
		clusters.levels	<- levels(df[,covariate])			
		clusters		<- lapply(clusters,function(x)
				{
					subset(df,clustername==x,c(fpos1,covariate,duration,hregion))
				})
		ylim			<- c(-15,length(clusters))
		cols			<- c(brewer.pal(length(clusters.levels)-1, "Set1"), "grey30")[c(2,1,3,4)]		#add grey for "unknown"
		#cols			<- c("green","red","blue","grey30")
		cols			<- sapply(cols, function(x) my.fade.col(x,0.7))
		pch				<- c(rep(19,length(clusters.levels)))
		xlim[1]			<- xlim[1]-700
		file			<- paste(dir.name,f.name,sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=yheight*7)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',xlab="")
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		#axis(3, )
		sapply(seq_along(clusters),function(i)
				{							
					cluster.cex	<- (1 - nrow(subset(clusters[[i]],hregion=="unknown")) / nrow(clusters[[i]])) * cex.points
					cluster.ix	<- order(as.numeric(clusters[[i]][,fpos1]))
					cluster.x	<- as.numeric(clusters[[i]][,fpos1])[cluster.ix]										
					cluster.y	<- rep(i,nrow(clusters[[i]]))
					cluster.z	<- as.numeric(clusters[[i]][,covariate])[cluster.ix]
					#if(cluster.z[1]!=4)
					lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
					points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=pch[1], cex=cluster.cex )										
				})		
		#plot small clusters
		abline(h=-5, lty=3,col="grey50")
		points( as.numeric(clusters.small[,fpos1]), rep(-5,nrow(clusters.small)), col=cols[ as.numeric(clusters.small[,covariate]) ], pch=pch[1], cex=0.75 )
		#plot singletons
		abline(h=-10, lty=4,col="grey50")
		points( as.numeric(singletons[,fpos1]), rep(-10,nrow(singletons)), col=cols[ as.numeric(singletons[,covariate]) ], pch=pch[1], cex=0.75 )
		
		pch[pch==19]	<- 21
		col				<- c(rep("transparent",length(clusters.levels)),"black")
		legend("topleft",bty='n',pt.bg=cols,pch=pch,legend=c(levels(df[,covariate])), col=col)
		legend("bottomleft", bty='n', lty=c(4), col="grey50", legend=c("singletons"))
		legend("bottomright", bty='n', lty=c(3), col="grey50", legend=c("networks<10"))
		dev.off()					
	}
	if(1)	#A, By riskgroup, blue MSM, red HT, green DU, black other. Annotate somehow those with high percentage / number of online retrieved seqs.
	{
		#output file
		f.name						<- paste("BezemerFig1_rout_wSingletons_",signat,".pdf",sep='')				
		#reduce risk group as desired and set as covariate
		tmp											<- unclass(df[,rout]	)
		code										<- seq_along(attr(tmp,"levels"))
		names(code)									<- attr(tmp,"levels")
		tmp[ tmp%in%code[c("BL","CH","unknown")] ]	<- code["unknown"]
		tmp											<- factor(tmp, levels=as.character(code[c("MSM","HT","DU","unknown")]), labels=c("MSM","HT","DU","other/unknown"))
		df[,covariate:=tmp ]
		tmp											<- df[,covariate]
		tmp[ which(tmp=="other/unknown") ]			<- NA
		df[,covariateNA:=tmp]
		
		xlim			<- range( df[,fpos1], na.rm=1 )
		zlim			<- range( df[,duration], na.rm=1 )
		#select singletons
		singletons		<- subset(df,clustername==0,c(fpos1,covariate,duration))
		#rm singletons from data table
		df				<- subset(df,clustername!=0,)
		#select clusters with less than CLUSTER.N.THRESHOLD non-NA fpos1 members
		cluster.n		<- df[,length(which(!is.na(fpos1))),by=clustername]		
		clusters.small	<- df[ subset(cluster.n,V1<CLUSTER.N.THRESHOLD,clustername) ]		
		#select clusters with more than CLUSTER.N.THRESHOLD non-NA fpos1 members
		df				<- df[ subset(cluster.n,V1>=CLUSTER.N.THRESHOLD,clustername) ]
		#sort clusters by duration -- this is the same as median duration
		clusters.sortby1<- as.vector( df[,as.numeric(any(hregion=="CU", na.rm=1)), by=clustername]$V1 )
		#clusters.sortby1<- as.vector( df[,median(as.numeric(covariateNA), na.rm=1), by=clustername]$V1 )
		clusters.sortby2<- as.vector( df[,median(duration, na.rm=1), by=clustername]$V1 )
		clusters.sortby	<- order(clusters.sortby1,clusters.sortby2, decreasing=T)
		
		#clusters.sortby	<- sort( as.vector( df[,median(duration, na.rm=1), by=clustername]$V1 ), index.return=T, decreasing=T)$ix
		clusters		<- unique( df[,clustername] )[clusters.sortby]
		#extract clusters one by one
		clusters.levels	<- levels(df[,covariate])			
		clusters		<- lapply(clusters,function(x)
				{
					subset(df,clustername==x,c(fpos1,covariate,duration,hregion))
				})
		ylim			<- c(-15,length(clusters))
		cols			<- c(brewer.pal(length(clusters.levels)-1, "Set1"), "grey30")[c(2,1,3,4)]		#add grey for "unknown"
		#cols			<- c("green","red","blue","grey30")
		cols			<- sapply(cols, function(x) my.fade.col(x,0.7))
		pch				<- c(rep(19,length(clusters.levels)))
		xlim[1]			<- xlim[1]-700
		file			<- paste(dir.name,f.name,sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=yheight*7)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',xlab="")
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		#axis(3, )
		sapply(seq_along(clusters),function(i)
				{							
					cluster.cex<- (1 - nrow(subset(clusters[[i]],hregion=="unknown")) / nrow(clusters[[i]])) * cex.points
					points( as.numeric(clusters[[i]][,fpos1]), rep(i,nrow(clusters[[i]])), col=cols[ as.numeric(clusters[[i]][,covariate]) ], pch=pch[1], cex=cluster.cex )										
				})		
		
		#plot small clusters
		abline(h=-5, lty=3,col="grey50")
		points( as.numeric(clusters.small[,fpos1]), rep(-5,nrow(clusters.small)), col=cols[ as.numeric(clusters.small[,covariate]) ], pch=pch[1], cex=0.75 )
		#plot singletons
		abline(h=-10, lty=4,col="grey50")
		points( as.numeric(singletons[,fpos1]), rep(-10,nrow(singletons)), col=cols[ as.numeric(singletons[,covariate]) ], pch=pch[1], cex=0.75 )
		
		pch[pch==19]	<- 21
		col				<- c(rep("transparent",length(clusters.levels)),"black")
		legend("topleft",bty='n',pt.bg=cols,pch=pch,legend=c(levels(df[,covariate])), col=col)
		legend("bottomleft", bty='n', lty=c(4), col="grey50", legend=c("singletons"))
		legend("bottomright", bty='n', lty=c(3), col="grey50", legend=c("networks<10"))
		dev.off()
		#save		
	}			
}

project.bezemer2013a.figs.v130517<- function()
{
	require(data.table)
	require(RColorBrewer)
	
	dir.name		<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"	
	signat			<- "130517"
	f.name			<- paste("Bezemer2014_clusters_",signat,".csv",sep='')
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
	
	dir.name			<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"
	CLUSTER.N.THRESHOLD	<- 10
	ylab				<- "transmission network ID"
	yheight				<- 1.5
	
	if(1)	#los alamos left out, figure x: fpos1 y: clusters, coloured by hregion
	{
		#input
		f.name						<- paste("BezemerFig1_hregion_noSingletons_",signat,".pdf",sep='')
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
		f.name						<- paste("BezemerFig1_rout_noSingletons_",signat,".pdf",sep='')
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
}

project.bezemer2013a.figs<- function()
{
	#project.bezemer2013a.figs.v130517()
	project.bezemer2013a.figs.v130528()
}
